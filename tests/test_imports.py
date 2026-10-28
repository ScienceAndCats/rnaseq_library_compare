import csv
import json
from pathlib import Path

import numpy as np
import pytest

from rnaseq_compare.io import (Catalog, InputError, build_matrix, comparison_pairs, count_sample_name,
    gc_from_fasta, load_config, load_sources, read_count_table, select_root)


def test_stage_sum_not_duplicate(tmp_path, fc_writer):
    fc_writer(tmp_path, ["S"], [[2], [3], [4]])
    fc_writer(tmp_path, ["S"], [[5], [7], [11]], stage="rna")
    config = {"_base": str(tmp_path), "inputs": [{"id": "a", "pipeline": "shortread", "path": str(tmp_path)}]}
    source = load_sources(config, [])[0]
    assert source.counts.S.tolist() == [7, 10, 15]
    assert sum(e-s+1 for _,s,e,_ in source.intervals["geneA"]) == 100


def test_longread_combined_once_and_length_not_doubled(tmp_path, fc_writer):
    path = fc_writer(tmp_path, ["S_R1", "S_R2"], longread=True)
    table = read_count_table(path, "combined", "bacteria")
    assert list(table.counts) == ["S_R1", "S_R2"]
    assert table.intervals["geneA"] == (("chrom", 1, 100, "+"),)
    config = {"_base": str(tmp_path), "inputs": [{"id": "long", "pipeline": "longread", "path": str(tmp_path)}]}
    sources = load_sources(config, [])
    catalog = Catalog(sources, config)
    assert catalog.resolve("S_R1.fastq.gz") == "long::S_R1"
    with pytest.raises(InputError):
        catalog.resolve("S")


def test_paired_r1_alias(tmp_path, fc_writer):
    fc_writer(tmp_path, ["MW147-J-EM0330_S2_L001"])
    config = {"_base": str(tmp_path), "inputs": [{"id": "a", "pipeline": "agilent_pe", "path": str(tmp_path)}]}
    catalog = Catalog(load_sources(config, []), config)
    assert catalog.resolve("MW147-J-EM0330_S2_L001_R1_001.fastq.gz") == "a::MW147-J-EM0330_S2_L001"


def test_ambiguous_name_requires_source(tmp_path, fc_writer):
    for sub in ("a", "b"):
        fc_writer(tmp_path/sub, ["same"])
    config = {"_base": str(tmp_path), "inputs": [{"id": x, "pipeline": "shortread", "path": str(tmp_path/x)} for x in ("a", "b")]}
    catalog = Catalog(load_sources(config, []), config)
    with pytest.raises(InputError, match="Ambiguous"):
        catalog.resolve("same")
    assert catalog.resolve("b::same") == "b::same"


def test_molecule_tables_not_summed_with_reads(tmp_path, fc_writer):
    fc_writer(tmp_path, ["S"], [[100], [200], [300]])
    fc_writer(tmp_path, ["S"], [[1], [2], [3]], layer="molecules")
    config = {"_base": str(tmp_path), "inputs": [{"id": "a", "pipeline": "agilent_pe", "path": str(tmp_path), "count_layer": "molecules"}]}
    source = load_sources(config, [])[0]
    assert source.counts.S.tolist() == [1,2,3]


def test_missing_stage_sample_fails(tmp_path, fc_writer):
    fc_writer(tmp_path, ["S", "T"])
    fc_writer(tmp_path, ["S"], stage="rna")
    with pytest.raises(InputError, match="sample columns differ"):
        load_sources({"_base": str(tmp_path), "inputs": [{"id":"a","pipeline":"shortread","path":str(tmp_path)}]}, [])


def test_missing_rna_cds_table_detected(tmp_path, fc_writer):
    fc_writer(tmp_path, ["S"])
    (tmp_path/"run_manifest.json").write_text(json.dumps({"read_layout":"single", "samples":[{"sample":"S"}]}))
    gff = tmp_path / "rna.gff3"
    gff.write_text("tx\tx\tCDS\t1\t100\t.\t+\t0\tID=g\n")
    (tmp_path/"reference_manifest.json").write_text(json.dumps({"bacteria":{"rna":{"gff":str(gff)}}}))
    with pytest.raises(InputError, match="CDS count table is missing"):
        load_sources({"_base":str(tmp_path),"inputs":[{"id":"a","pipeline":"shortread","path":str(tmp_path)}]}, [])


def test_multiple_runs_refuse_latest_guess(tmp_path, fc_writer):
    for name in ("old", "new"):
        fc_writer(tmp_path/"slurm_runs"/name, ["S"], longread=True)
    with pytest.raises(InputError, match="Multiple run directories"):
        select_root(tmp_path, "bacteria", "longread", "reads")


def test_unfinished_longread_aggregation_fails(tmp_path, fc_writer):
    fc_writer(tmp_path, ["S"], longread=True)
    (tmp_path/"run.json").write_text("{}")
    with pytest.raises(InputError, match="not marked complete"):
        select_root(tmp_path,"bacteria","longread","reads")
    (tmp_path/"aggregation.completed.json").write_text("{}")
    assert select_root(tmp_path,"bacteria","longread","reads") == tmp_path


def test_fractional_counts_fail(tmp_path, fc_writer):
    path = fc_writer(tmp_path, ["S"])
    path.write_text(path.read_text().replace("\t10\n", "\t10.5\n"))
    with pytest.raises(InputError, match="Noninteger"):
        read_count_table(path, "dna", "bacteria")


def test_duplicate_genes_fail(tmp_path, fc_writer):
    path = fc_writer(tmp_path, ["S"], genes=["same","same"])
    with pytest.raises(InputError, match="Duplicate"):
        read_count_table(path, "dna", "bacteria")


def test_duplicate_input_matrix_rejected(tmp_path, fc_writer):
    fc_writer(tmp_path, ["S"])
    config = {"_base":str(tmp_path), "inputs":[{"id":x,"pipeline":"shortread","path":str(tmp_path)} for x in ("a","b")]}
    with pytest.raises(InputError, match="same count matrix"):
        load_sources(config, [])


def test_gene_intersection_not_zero_fill(tmp_path, fc_writer):
    fc_writer(tmp_path/"a", ["S"], genes=["geneA","geneB","only_A"])
    fc_writer(tmp_path/"b", ["T"], genes=["geneA","geneB","only_B"])
    config = {"_base":str(tmp_path),"inputs":[{"id":x,"pipeline":"shortread","path":str(tmp_path/x)} for x in ("a","b")],
              "analysis":{"min_shared_gene_fraction":.6}}
    sources = load_sources(config, [])
    matrix, *_ = build_matrix(sources,["a::S","b::T"],config,[])
    assert list(matrix.index) == ["geneA","geneB"]
    config["analysis"]["min_shared_gene_fraction"] = .8
    with pytest.raises(InputError, match="share exact IDs"):
        build_matrix(sources,["a::S","b::T"],config,[])


def test_coordinate_mismatch_fails(simple_config):
    sources = load_sources(simple_config, [])
    sources[1].intervals["geneA"] = (("different_assembly",1,100,"+"),)
    with pytest.raises(InputError, match="coordinates disagree"):
        build_matrix(sources,["a::A1","b::B1"],simple_config,[])


def test_long_yaml_lists_and_named_groups(simple_config):
    cat = Catalog(load_sources(simple_config, []),simple_config)
    simple_config["comparison_groups"] = [["A1","B1","B2"],["A1","B3"]]
    pairs = comparison_pairs(simple_config,cat)
    assert [p["query"] for p in pairs] == ["b::B1","b::B2","b::B3"]
    assert {p["reference"] for p in pairs} == {"a::A1"}


def test_duplicate_yaml_key_fails(tmp_path):
    path = tmp_path/"config.yaml"
    path.write_text("version: 1\nversion: 2\n")
    with pytest.raises(InputError, match="Duplicate YAML"):
        load_config(path)


def test_gc_streaming_including_ambiguous_bases(tmp_path):
    path = tmp_path/"ref.fa"
    path.write_text(">chrom\nACGTNN\nGGCCAA\n")
    gc = gc_from_fasta(path, {"a":(("chrom",1,4,"+"),),"b":(("chrom",5,12,"-"),)})
    assert gc["a"] == 50
    assert gc["b"] == pytest.approx(400/6)
    with pytest.raises(InputError, match="shorter"):
        gc_from_fasta(path, {"a":(("chrom",1,100,"+"),)})
    with pytest.raises(InputError, match="absent"):
        gc_from_fasta(path, {"a":(("missing",1,2,"+"),)})


def test_large_featurecounts_fields(tmp_path):
    path=tmp_path/"long.txt"
    n=23000
    with path.open("w",newline="") as handle:
        writer=csv.writer(handle,delimiter="\t")
        writer.writerow(["Geneid","Chr","Start","End","Strand","Length","BACTERIA_S_22bp.sam"])
        writer.writerow(["g",";".join(["chromosome"]*n),";".join(["1"]*n),";".join(["100"]*n),";".join(["+"]*n),100,10])
    table=read_count_table(path,"dna","bacteria")
    assert table.counts.loc["g","S"] == 10
    assert table.intervals["g"] == (("chromosome",1,100,"+"),)


def test_zero_library_not_normalized(simple_config):
    sources=load_sources(simple_config,[])
    sources[0].counts["A1"]=0
    with pytest.raises(InputError,match="zero counts"):
        build_matrix(sources,["a::A1","b::B1"],simple_config,[])


def test_auto_fasta_from_short_manifest(tmp_path, fc_writer):
    from rnaseq_compare.io import detect_genome_fasta
    fc_writer(tmp_path,["S"])
    ref=tmp_path/".reference_cache/hash"
    ref.mkdir(parents=True)
    (ref/"dna.fa").write_text(">chrom\nACGT\n")
    (tmp_path/"reference_manifest.json").write_text(json.dumps({"bacteria":{"dna":{"gff":str(ref/"dna.gff3"),"index":str(ref/"dna")}}}))
    sources=load_sources({"_base":str(tmp_path),"inputs":[{"id":"a","pipeline":"shortread","path":str(tmp_path)}]},[])
    path,sid=detect_genome_fasta(sources,"bacteria")
    assert path == ref/"dna.fa"
    assert sid == "a"


def test_one_to_one_gene_mapping(tmp_path, fc_writer):
    path=fc_writer(tmp_path,["S"])
    table=read_count_table(path,"dna","bacteria",{"geneA":"newA"})
    assert "newA" in table.counts.index
    with pytest.raises(InputError,match="Duplicate"):
        read_count_table(path,"dna","bacteria",{"geneA":"geneB"})


def test_alias_metadata_cannot_conflict(simple_config):
    simple_config["samples"]["second_alias"]={"source":"a","sample":"A1","biological_sample":"DIFFERENT"}
    with pytest.raises(InputError,match="Conflicting metadata"):
        Catalog(load_sources(simple_config,[]),simple_config)


def test_host_table_import(tmp_path, fc_writer):
    fc_writer(tmp_path,["S"],organism="host")
    config={"_base":str(tmp_path),"organism":"host","inputs":[{"id":"a","pipeline":"shortread","path":str(tmp_path)}]}
    source=load_sources(config,[])[0]
    assert list(source.counts.columns) == ["S"]
