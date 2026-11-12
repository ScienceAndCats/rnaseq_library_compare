"""Real SAM/GFF fixtures exercising the category analyses end to end."""
import csv
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest
import yaml

from rnaseq_compare.cli import main
from rnaseq_compare.io import InputError, load_sources, Catalog, comparison_pairs
from rnaseq_compare.de import prepare_plans
from rnaseq_compare.rna_annotation import RNAIndex, SUBTYPES, subtype
from rnaseq_compare.rna_categories import (CATEGORIES, Stage, classify_stage, member,
                                          prepare_categories, discover_stages, CategoryData,
                                          analyze_categories, MissingCategoryInput)


def sam_file(path, reads, length=3000, sort="unsorted"):
    path.parent.mkdir(parents=True, exist_ok=True)
    with pysam.AlignmentFile(str(path), "w", header={"HD": {"VN": "1.6", "SO": sort}, "SQ": [{"SN": "chrom", "LN": length}]}) as handle:
        for name, start, cigar, flag in reads:
            read = pysam.AlignedSegment(handle.header)
            read.query_name, read.flag = name, flag
            read.reference_id, read.reference_start = 0, start
            read.mapping_quality = 42
            read.cigarstring = cigar
            read.query_sequence = "A" * read.infer_query_length()
            if flag & 1:
                read.next_reference_id, read.next_reference_start = 0, start
            handle.write(read)
    return path


@pytest.fixture
def annotation(tmp_path):
    path = tmp_path / "annotation.gff3"
    path.write_text("##gff-version 3\n"
        "chrom\ttest\tCDS\t101\t200\t.\t+\t0\tlocus_tag=geneA\n"
        "chrom\ttest\trRNA\t151\t180\t.\t+\t.\tlocus_tag=rrna\n"
        "chrom\ttest\tncRNA\t301\t380\t.\t+\t.\tlocus_tag=small\n"
        "chrom\ttest\ttRNA\t401\t480\t.\t-\t.\tlocus_tag=transfer\n"
        "chrom\ttest\tncRNA\t501\t580\t.\t+\t.\tlocus_tag=antisense;ncRNA_class=antisense_RNA\n"
        "chrom\ttest\tmisc_RNA\t601\t680\t.\t+\t.\tlocus_tag=misc\n"
        "chrom\ttest\tfive_prime_UTR\t701\t780\t.\t+\t.\tlocus_tag=utr\n"
        "chrom\ttest\tmRNA\t1001\t1500\t.\t+\t.\tID=tx\n"
        "chrom\ttest\texon\t1001\t1100\t.\t+\t.\tParent=tx\n"
        "chrom\ttest\texon\t1401\t1500\t.\t+\t.\tParent=tx\n")
    return path


def test_blocks_priority_strand_and_subtypes(annotation):
    idx = RNAIndex(annotation, {"chrom": 3000})
    assert idx.classify("chrom", [(155, 165)], "+")[0] == "rRNA"
    assert idx.classify("chrom", [(100, 110)], "+")[0] == "CDS"
    assert idx.classify("chrom", [(300, 310)], "+")[0] == "sRNA_ncRNA"
    assert idx.classify("chrom", [(400, 410)], "-")[0] == "tRNA_tmRNA"
    assert idx.classify("chrom", [(400, 410)], "+")[0] == "intergenic_unannotated"
    assert idx.classify("chrom", [(500, 510)], "+")[0] == "antisense_RNA"
    assert idx.classify("chrom", [(600, 610)], "+")[0] == "other_noncoding_RNA"
    assert idx.classify("chrom", [(700, 710)], "+")[0] == "UTR_transcript_regions"
    assert idx.classify("chrom", [(1050, 1060)], "+")[0] == "UTR_transcript_regions"
    assert idx.classify("chrom", [(1200, 1210)], "+")[0] == "intergenic_unannotated"
    # A spliced read's unaligned gap overlaps CDS/rRNA; the blocks do not.
    assert idx.classify("chrom", [(80, 90), (210, 220)], "+")[0] == "intergenic_unannotated"
    assert idx.classify("chrom", [(300, 310)], None)[0] == "sRNA_ncRNA"
    assert idx.relative_position("tRNA_tmRNA::transfer", "chrom", [(400, 410)]) > .9
    assert idx.relative_position("sRNA_ncRNA::small", "chrom", [(300, 310)]) < .1


def test_annotation_validation_and_rna_coordinate_lengths(annotation, tmp_path):
    with pytest.raises(InputError, match="alignment header"):
        RNAIndex(annotation, {"different": 3000})
    with pytest.raises(InputError, match="alignment header"):
        RNAIndex(annotation, {"chrom": 100})
    idx = RNAIndex(annotation, {"chrom": 3000}, kind="rna")
    assert all(not spans for spans in idx.intervals.values())
    assert idx.classify("chrom", [(1200, 1210)], "+")[0] == "UTR_transcript_regions"
    empty = tmp_path / "empty.gff3"
    empty.write_text("##gff-version 3\n")
    with pytest.raises(InputError, match="No eligible"):
        RNAIndex(empty, {"chrom": 3000})


def test_reverse_proper_pairs_and_unsequenced_gap(annotation, tmp_path):
    # Coordinate order separates the proper mates; sort by name before counting.
    path = sam_file(tmp_path / "paired.sam", [
        ("cds_gap", 80, "10M", 99), ("small", 300, "10M", 163),
        ("cds_gap", 210, "10M", 147), ("small", 340, "10M", 83)])
    index = RNAIndex(annotation, {"chrom": 3000}, strand="reverse")
    total, counts, ambiguous = classify_stage(Stage("S", "dna", path, annotation, "reverse"), True, index)
    assert total["sRNA_ncRNA"] == 1
    assert total["intergenic_unannotated"] == 1
    assert total["CDS"] == total["rRNA"] == 0
    assert sum(total.values()) == 2


def test_ambiguous_counts_retained_but_not_invented_features(tmp_path):
    gff = tmp_path / "ambiguous.gff3"
    gff.write_text("chrom\tt\tncRNA\t101\t200\t.\t+\t.\tID=a\n"
                   "chrom\tt\tncRNA\t151\t250\t.\t+\t.\tID=b\n")
    path = sam_file(tmp_path / "amb.sam", [("multi", 160, "10M", 0), ("secondary", 100, "10M", 256), ("supp", 100, "10M", 2048)])
    total, counts, ambiguous = classify_stage(Stage("S", "dna", path, gff, "forward"), False, RNAIndex(gff, {"chrom": 3000}))
    assert total == {"sRNA_ncRNA": 1}
    assert not counts
    assert ambiguous == {"sRNA_ncRNA": 1}
    with pytest.raises(InputError, match="merges distinct RNA features"):
        RNAIndex(gff, {"chrom": 3000}, gene_map={"a": "b"})


def test_upstream_totals_must_match(annotation, tmp_path):
    path = sam_file(tmp_path / "s.sam", [("small", 300, "10M", 0)])
    expected = tmp_path / "s.categories.json"
    expected.write_text(json.dumps({"aligned": 1, "categories": {"rRNA": 0, "CDS": 1, "other_RNA": 0, "intergenic_unannotated": 0}}))
    with pytest.raises(InputError, match="Reconstructed categories disagree"):
        classify_stage(Stage("S", "dna", path, annotation, "forward", expected), False, RNAIndex(annotation, {"chrom": 3000}))


@pytest.fixture
def category_config(tmp_path, fc_writer, annotation):
    root = tmp_path / "source"
    fc_writer(root, ["A", "B"], [[12, 25], [30, 35], [10, 11]])
    gff = root / "annotation.gff3"
    gff.write_text(annotation.read_text())
    stages = []
    for sample, extra in (("A", 0), ("B", 1)):
        specs = []
        for i, (start, flag) in enumerate(((100, 0), (155, 0), (300, 0), (400, 16), (500, 0), (600, 0), (700, 0), (1900, 0))):
            specs.extend((f"{i}_{n}", start, "10M", flag) for n in range(12+i+extra))
        alignment = sam_file(root / "bowtie_alignments/bacteria" / f"BACTERIA_{sample}_22bp.sam", specs)
        category = alignment.with_suffix(".categories.json")
        category.write_text(json.dumps({"aligned": len(specs), "categories": {"CDS": 12+extra, "rRNA": 13+extra,
             "other_RNA": sum(12+i+extra for i in range(2, 7)), "intergenic_unannotated": 19+extra}}))
        stages.append({"sample": sample, "stages": [{"organism": "bacteria", "kind": "dna", "categories": str(category)}]})
    (root / "run_manifest.json").write_text(json.dumps({"read_layout": "single", "samples": stages}))
    (root / "reference_manifest.json").write_text(json.dumps({"bacteria": {"dna": {"gff": str(gff)}, "rna": None}}))
    return {"version": 1, "_base": str(tmp_path), "output_dir": "results", "reference": {"genome_fasta": None},
            "inputs": [{"id": "short", "path": "source", "pipeline": "shortread", "library_prep": "method"}],
            "comparison_groups": {"one": ["A", "B"]}, "analysis": {"min_count": 1},
            "rna_categories": {"write_read_assignments": True}}


def test_all_category_analyses_and_csvs(category_config, tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(yaml.safe_dump({k: v for k, v in category_config.items() if not k.startswith("_")}))
    assert main([str(path), "--check"]) == 0
    assert not (tmp_path / "results").exists()
    assert main([str(path), "--debug"]) == 0
    output = tmp_path / "results"
    status = pd.read_csv(output / "rna_category_analysis_status.csv")
    assert set(status.category) == set(CATEGORIES)
    assert (status.status == "completed").all()
    totals = pd.read_csv(output / "rna_category_counts.csv")
    for sample in totals.sample_id.unique():
        part = totals[totals.sample_id == sample].set_index("category")["count"]
        assert part.bacteria_aligned == part.bacteria_rRNA + part.bacteria_CDS + part.bacteria_other_RNA + part.bacteria_intergenic_unannotated
        assert part.bacteria_non_rRNA == part.bacteria_aligned - part.bacteria_rRNA
        assert part.other_RNA_total == sum(part[c] for c in SUBTYPES)
        assert part.other_RNA_total == part.bacteria_other_RNA
    for cat in CATEGORIES:
        directory = output / "rna_categories" / cat
        assert (output / "rna_categories" / f"{cat}.csv").is_file()
        assert len(pd.read_csv(directory / "pairwise_metrics.csv")) == 1
        assert len(pd.read_csv(directory / "comparison_index.csv")) == 1
        assert (directory / "differential_expression_status.csv").is_file()
        assert (directory / "technical_bias_associations.csv").is_file()
        cpm = pd.read_csv(directory / "counts_cpm.csv", index_col=0)
        np.testing.assert_allclose(cpm.sum(), 1e6)
    assignments = pd.read_csv(output / "rna_categories/read_assignments.csv")
    assert len(assignments) == totals[totals.category == "bacteria_aligned"]["count"].sum()
    profiles = pd.read_csv(output / "rna_categories/tRNA_tmRNA/category_metagene_profiles.csv")
    assert profiles.loc[profiles.mean_fraction > 0, "bin_start_percent"].min() > 80
    manifest = json.loads((output / "run_manifest.json").read_text())
    assert manifest["status"] == "completed"
    assert len(manifest["rna_category_analyses"]) == 12


def test_missing_inputs_totals_only_and_strict_mode(category_config, tmp_path):
    root = tmp_path / "source"
    for alignment in root.glob("bowtie_alignments/bacteria/*.sam"):
        alignment.unlink()
    with (root / "output_read_breakdown.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Sample", "Bacteria aligned (reads)", "Bacteria other RNA (reads)"])
        writer.writeheader()
        writer.writerows({"Sample": sample, "Bacteria aligned (reads)": "50", "Bacteria other RNA (reads)": "20"} for sample in ("A", "B"))
    sources = load_sources(category_config, [])
    data = prepare_categories(sources, ["short::A", "short::B"], category_config, [])
    assert len(data.missing) == 1
    totals = pd.DataFrame(data.totals).set_index(["sample_id", "category"])
    assert totals.loc[("short::A", "bacteria_aligned"), "count"] == 50
    assert totals.loc[("short::A", "other_RNA_total"), "status"] == "totals_only"
    assert totals.loc[("short::A", "sRNA_ncRNA"), "status"] == "unavailable"
    category_config["rna_categories"]["require_alignments"] = True
    with pytest.raises(InputError, match="retained alignment is unavailable"):
        prepare_categories(sources, ["short::A", "short::B"], category_config, [])


def test_zero_category_is_not_fabricated_normalization(category_config, tmp_path):
    path = tmp_path / "source/bowtie_alignments/bacteria/BACTERIA_B_22bp.sam"
    with pysam.AlignmentFile(str(path)) as handle:
        reads = [(r.query_name, r.reference_start, r.cigarstring, r.flag) for r in handle if not r.query_name.startswith("2_")]
    sam_file(path, reads)
    path.with_suffix(".categories.json").unlink()
    category_config["rna_categories"]["categories"] = ["sRNA_ncRNA"]
    config = tmp_path / "zero.yaml"
    config.write_text(yaml.safe_dump({k: v for k, v in category_config.items() if not k.startswith("_")}))
    assert main([str(config), "--debug"]) == 0
    directory = tmp_path / "results/rna_categories/sRNA_ncRNA"
    cpm = pd.read_csv(directory / "counts_cpm.csv", index_col=0)
    assert cpm["short::B"].isna().all()
    assert pd.read_csv(directory / "comparison_status.csv").iloc[0].status == "zero_category_counts"


def test_explicit_stages_and_no_molecule_read_substitution(category_config, tmp_path):
    source = load_sources(category_config, [])[0]
    override = category_config["inputs"][0]
    override["category_stages"] = [{"sample": "A", "kind": "dna", "gff": "source/annotation.gff3",
                                   "alignment": "source/bowtie_alignments/bacteria/BACTERIA_A_22bp.sam"}]
    assert discover_stages(source, override, category_config)[0].sample == "A"
    source.layer = "molecules"
    del override["category_stages"]
    with pytest.raises(MissingCategoryInput, match="retained alignment is unavailable"):
        discover_stages(source, override, category_config)


def test_rna_stage_without_cds_table_and_selected_samples_only(category_config, tmp_path):
    root = tmp_path / "source"
    gff = root / "rna.gff3"
    gff.write_text("chrom\tt\ttRNA\t401\t480\t.\t-\t.\tlocus_tag=transfer\n")
    refs = json.loads((root / "reference_manifest.json").read_text())
    refs["bacteria"]["rna"] = {"gff": str(gff)}
    (root / "reference_manifest.json").write_text(json.dumps(refs))
    manifest = json.loads((root / "run_manifest.json").read_text())
    for entry in manifest["samples"]:
        path = sam_file(root / "bowtie_alignments/bacteria/rna" / f"BACTERIA_{entry['sample']}_22bp.sam", [("rna_trna", 400, "10M", 16)])
        category = path.with_suffix(".categories.json")
        category.write_text(json.dumps({"aligned": 1, "categories": {"rRNA": 0, "CDS": 0, "other_RNA": 1, "intergenic_unannotated": 0}}))
        entry["stages"].append({"organism": "bacteria", "kind": "rna", "categories": str(category)})
    (root / "run_manifest.json").write_text(json.dumps(manifest))
    # Unselected B has no retained DNA SAM; this must not block A.
    (root / "bowtie_alignments/bacteria/BACTERIA_B_22bp.sam").unlink()
    sources = load_sources(category_config, [])
    data = prepare_categories(sources, ["short::A"], category_config, [])
    assert not data.missing
    assert {p["stage"] for p in data.provenance} == {"rna", "dna"}
    total = next(t for t in data.totals if t["category"] == "tRNA_tmRNA")
    assert total["count"] == 15 + 1


def test_same_category_alignment_cannot_supply_two_libraries(category_config):
    spec = category_config["inputs"][0]
    spec["category_stages"] = [{"sample": sample, "kind": "dna", "gff": "source/annotation.gff3",
                               "alignment": "source/bowtie_alignments/bacteria/BACTERIA_A_22bp.sam"} for sample in ("A", "B")]
    with pytest.raises(InputError, match="same category alignment"):
        prepare_categories(load_sources(category_config, []), ["short::A", "short::B"], category_config, [])


def test_category_specific_actual_de_and_untestable_panel(simple_config, tmp_path):
    simple_config["rna_categories"] = {"categories": ["sRNA_ncRNA"]}
    simple_config["analysis"] = {"min_count": 10, "import_positional_profiles": False}
    simple_config["differential_expression"] = {"enabled": True, "analyses": [{"name": "paired_small_RNA",
        "reference": ["A1", "A2", "A3"], "query": ["B1", "B2", "B3"], "paired_by": "biological_sample"}]}
    sources = load_sources(simple_config, [])
    catalog = Catalog(sources, simple_config)
    plans = prepare_plans(simple_config, catalog, [])
    selected = plans[0].samples
    rng = np.random.default_rng(123)
    genes = [f"sRNA_ncRNA::g{i:03d}" for i in range(250)]
    means = rng.uniform(80, 500, size=(250, 1)) * np.ones((1, 6))
    means[:20, 3:] *= 8
    counts = pd.DataFrame(rng.negative_binomial(25, 25/(25+means)), index=pd.Index(genes, name="gene_id"), columns=selected)
    data = CategoryData()
    for source in sources:
        keys = [f"{source.id}::{s}" for s in source.counts]
        frame = counts[keys].rename(columns={k: k.split("::", 1)[1] for k in keys})
        data.sources.append(replace(source, counts=frame, intervals={}))
        data.classes[source.id] = dict.fromkeys(genes, "sRNA_ncRNA")
    data.totals = [{"sample_id": key, "category": "sRNA_ncRNA", "status": "reconstructed", "count": int(counts[key].sum()),
                    "bacteria_aligned_count": int(counts[key].sum())} for key in selected]
    output = tmp_path / "de_results"
    output.mkdir()
    analyze_categories(data, selected, simple_config, [], [], plans, catalog.metadata(selected), None, 1, output, [])
    category = output / "rna_categories/sRNA_ncRNA"
    status = pd.read_csv(category / "differential_expression_status.csv")
    assert status.iloc[0].status == "completed"
    results = pd.read_csv(next((category / "differential_expression").glob("*/results.csv")), index_col=0)
    assert results.index.str.startswith("sRNA_ncRNA::").all()
    assert results.iloc[:20].log2FoldChange.median() > 1.5
    assert results.pvalue.notna().any()
    assert results.padj.dropna().between(0, 1).all()
    metadata = pd.read_csv(next((category / "differential_expression").glob("*/sample_metadata.csv")))
    assert set(metadata.counting_policy) == {"exclusive_category_unique_feature"}
    # The same design cannot fabricate a model on a one-feature RNA panel.
    for source in data.sources:
        source.counts = source.counts.iloc[:1]
    output = tmp_path / "insufficient_panel"
    output.mkdir()
    analyze_categories(data, selected, simple_config, [], [], plans, catalog.metadata(selected), None, 1, output, [])
    status = pd.read_csv(output / "rna_categories/sRNA_ncRNA/differential_expression_status.csv")
    assert status.iloc[0].status == "not_testable"
    assert "fewer than two genes" in status.iloc[0].reason


@pytest.mark.parametrize("typ,attrs,expected", [
    ("sRNA", {}, "sRNA_ncRNA"), ("tmRNA", {}, "tRNA_tmRNA"),
    ("ncRNA", {"ncRNA_class": "antisense_RNA"}, "antisense_RNA"),
    ("misc_RNA", {}, "other_noncoding_RNA"), ("five_prime_UTR", {}, "UTR_transcript_regions")])
def test_subtype_mapping(typ, attrs, expected):
    assert subtype(typ.lower(), attrs) == expected
