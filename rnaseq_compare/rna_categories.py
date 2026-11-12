"""Reconstruct RNA categories from retained alignments and analyze each panel.

Disposition totals alone cannot reconstruct per-feature expression. They are
exported separately when the alignments or normalized GFFs were not retained.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field, replace
from itertools import groupby
import csv
from pathlib import Path
import re
import tempfile

import numpy as np
import pandas as pd

from .analysis import (depth_analysis, global_qc, pairwise_analysis,
                       reproducibility_analysis, write_csv, save)
from .de import count_filter, run_de
from .io import (InputError, build_matrix, config_path, file_sha256, gc_from_fasta,
                 merge_intervals, read_disposition, read_gene_maps, read_json)
from .rna_annotation import RNAIndex, SUBTYPES


CATEGORIES = ("bacteria_aligned", "bacteria_non_rRNA", "bacteria_rRNA", "bacteria_CDS",
              "bacteria_other_RNA", "bacteria_intergenic_unannotated", "other_RNA_total",
              "sRNA_ncRNA", "tRNA_tmRNA", "antisense_RNA", "other_noncoding_RNA", "UTR_transcript_regions")
LABELS = dict(zip(CATEGORIES[:6], ("Bacteria aligned", "Bacteria non-rRNA", "Bacteria rRNA", "Bacteria CDS",
                                 "Bacteria other RNA", "Bacteria intergenic/unannotated")))


def member(cls, category):
    if category == "bacteria_aligned":
        return True
    if category == "bacteria_non_rRNA":
        return cls != "rRNA"
    if category in {"bacteria_other_RNA", "other_RNA_total"}:
        return cls in SUBTYPES
    return cls == {"bacteria_rRNA": "rRNA", "bacteria_CDS": "CDS",
                   "bacteria_intergenic_unannotated": "intergenic_unannotated"}.get(category, category)


def settings(config):
    block = config.get("rna_categories", {})
    if block is None:
        block = {}
    if not isinstance(block, dict):
        raise InputError("rna_categories must be a mapping.")
    unknown = set(block) - {"enabled", "require_alignments", "categories", "intergenic_bin_bp", "write_read_assignments"}
    if unknown:
        raise InputError(f"Unknown rna_categories settings: {sorted(unknown)}")
    for key in ("enabled", "require_alignments", "write_read_assignments"):
        if key in block and not isinstance(block[key], bool):
            raise InputError(f"rna_categories.{key} must be true or false.")
    chosen = block.get("categories", list(CATEGORIES))
    if not isinstance(chosen, list) or not chosen or any(c not in CATEGORIES for c in chosen) or len(chosen) != len(set(chosen)):
        raise InputError(f"rna_categories.categories must list distinct names from {list(CATEGORIES)}")
    bin_bp = block.get("intergenic_bin_bp", 1000)
    if isinstance(bin_bp, bool) or not isinstance(bin_bp, int) or bin_bp < 1:
        raise InputError("rna_categories.intergenic_bin_bp must be a positive integer.")
    return dict(block, categories=chosen, intergenic_bin_bp=bin_bp)


def resolve_retained(value, source, base=None):
    """Resolve recorded or recognizably relocated files, refusing conflicting copies."""
    if not value:
        return None
    original = Path(value)
    candidates = [original if original.is_absolute() else (base or source.root) / original]
    for anchor in ("samples", "bowtie_alignments", "bbsplit_alignments", ".reference_cache", "shared"):
        if anchor in original.parts:
            offset = original.parts.index(anchor)
            candidates.append(source.root.joinpath(*original.parts[offset:]))
    if base:
        candidates.append(base / original.name)
    found = list(dict.fromkeys(p.resolve() for p in candidates if p.is_file()))
    if len(found) > 1 and len({file_sha256(p) for p in found}) > 1:
        raise InputError(f"Conflicting retained files for {source.id}: {found}")
    return found[0] if found else None


@dataclass
class Stage:
    sample: str
    kind: str
    alignment: Path
    gff: Path
    strand: str
    expected: Path | None = None


class MissingCategoryInput(Exception):
    """Optional category reconstruction is unavailable, rather than invalid."""


def discover_stages(source, spec, config, samples=None):
    samples = set(source.counts.columns if samples is None else samples)
    recorded_strands = {match for t in source.tables for line in t.comments
                        for match in re.findall(r"(?:^|\s)-s\s+([012])(?:\s|$)", line)}
    if len(recorded_strands) > 1 and "category_strand" not in spec:
        raise InputError(f"{source.id}: inconsistent featureCounts strands; set category_strand explicitly.")
    inferred = {"0": "unstranded", "1": "forward", "2": "reverse"}.get(next(iter(recorded_strands), ""))
    strand = spec.get("category_strand", inferred or ("reverse" if source.pipeline == "agilent_pe" else "forward"))
    if strand not in {"forward", "reverse", "unstranded"}:
        raise InputError(f"{source.id}: category_strand must be forward, reverse, or unstranded.")
    if source.layout == "unknown":
        raise MissingCategoryInput("Set read_layout: single or paired to classify reads safely.")
    overrides = spec.get("category_stages")
    if overrides is not None:
        if not isinstance(overrides, list) or not overrides:
            raise InputError(f"{source.id}: category_stages must be a nonempty list.")
        result = []
        for item in overrides:
            if not isinstance(item, dict) or not {"sample", "kind", "alignment", "gff"} <= item.keys():
                raise InputError("Each category_stages entry needs sample, kind, alignment, and gff.")
            if item["sample"] not in source.counts.columns or item["kind"] not in {"dna", "rna"}:
                raise InputError(f"{source.id}: invalid category stage sample/kind: {item}")
            if item["sample"] not in samples:
                continue
            paths = [config_path(item[k], Path(config["_base"])) for k in ("alignment", "gff")]
            if not all(p.is_file() for p in paths):
                raise MissingCategoryInput(f"Explicit category alignment/GFF is missing: {paths}")
            result.append(Stage(item["sample"], item["kind"], *paths, strand))
        return result
    ref_candidates = [source.root / "reference_manifest.json",
                      source.root / "bbsplit_alignments/references/manifest.json",
                      source.root / "shared/bbsplit_alignments/references/manifest.json"]
    ref_path = next((p for p in ref_candidates if p.is_file()), None)
    if ref_path is None:
        raise MissingCategoryInput("No reference manifest; provide category_stages with retained alignments and normalized GFFs.")
    refs = read_json(ref_path).get("bacteria") or {}
    gffs = {kind: resolve_retained((refs.get(kind) or {}).get("gff"), source, ref_path.parent) for kind in ("rna", "dna")}
    for kind in ("rna", "dna"):
        if refs.get(kind) and gffs[kind] is None:
            raise MissingCategoryInput(f"Normalized {kind} GFF is unavailable; restore it or provide category_stages.")
    result = []
    entries = source.metadata.get("samples") or []
    if entries and source.pipeline != "longread":
        for sample in entries:
            if sample["sample"] not in samples:
                continue
            for stage in sample.get("stages", []):
                if stage.get("organism") != "bacteria":
                    continue
                kind = stage.get("kind")
                if kind not in gffs or gffs[kind] is None:
                    raise MissingCategoryInput(f"{sample['sample']}: normalized GFF is unavailable for stage {kind}.")
                category_path = stage.get("molecule_categories" if source.layer == "molecules" else "categories")
                original = Path(category_path) if category_path else None
                recorded = stage.get("molecules") if source.layer == "molecules" else (str(original).removesuffix(".categories.json") + ".sam" if original else None)
                alignment = resolve_retained(recorded, source)
                if alignment is None and recorded and source.layer == "reads":
                    alignment = resolve_retained(str(Path(recorded).with_suffix(".sorted.bam")), source)
                if alignment is None:
                    raise MissingCategoryInput(f"{sample['sample']} {kind}: retained alignment is unavailable ({recorded}).")
                expected = resolve_retained(category_path, source)
                result.append(Stage(sample["sample"], kind, alignment, gffs[kind], strand, expected))
        if result:
            return result
    # Legacy stage matrices or long-read combined matrices: use their recorded
    # sample alignment paths, never recursively glob old sample-worker runs.
    for sample in sorted(samples):
        tables = {t.stage: t for t in source.tables if sample in t.columns}
        for kind in ("rna", "dna"):
            if not refs.get(kind):
                if kind in tables:
                    raise MissingCategoryInput(f"{sample}: {kind} count table exists but its reference annotation is not recorded.")
                continue
            table = tables.get(kind) or next(iter(tables.values()))
            path = Path(table.columns[sample])
            if table.stage == "combined":
                path = path.parent / kind / path.name
            elif table.stage != kind:
                parent = path.parent.parent if table.stage == "rna" else path.parent
                path = (parent / "rna" if kind == "rna" else parent) / path.name
            alignment = resolve_retained(str(path), source)
            if alignment is None and source.layer == "reads":
                alignment = resolve_retained(str(path.with_suffix(".sorted.bam")), source)
            if alignment is None:
                raise MissingCategoryInput(f"{sample} {kind}: retained stage alignment is unavailable ({path}).")
            expected = resolve_retained(str(path.with_suffix(".categories.json")), source)
            result.append(Stage(sample, kind, alignment, gffs[kind], strand, expected))
    if not result:
        raise MissingCategoryInput("No explicit bacterial alignment stages could be resolved; provide category_stages.")
    return result


@dataclass
class CategoryData:
    sources: list = field(default_factory=list)
    classes: dict = field(default_factory=dict)
    totals: list = field(default_factory=list)
    missing: dict = field(default_factory=dict)
    provenance: list = field(default_factory=list)
    positions: dict = field(default_factory=dict)
    matrices: dict = field(default_factory=dict)
    assignments: object = None


def classify_stage(stage, paired, index, assignment_writer=None, sample_id="", positions=None):
    import pysam

    totals, counts, ambiguous = Counter(), Counter(), Counter()
    with tempfile.TemporaryDirectory(prefix="rnaseq-categories-") as temp:
        path = stage.alignment
        # Proper mates must be grouped even if only coordinate-sorted BAMs
        # remain. Name sorting is local; upstream alignments are never changed.
        if paired:
            sorted_path = Path(temp) / "name_sorted.bam"
            pysam.sort("-n", "-m", "256M", "-o", str(sorted_path), str(path))
            path = sorted_path
        with pysam.AlignmentFile(str(path), "rb" if path.suffix == ".bam" else "r") as handle:
            lengths = dict(zip(handle.references, handle.lengths))
            for name, group in groupby(handle, key=lambda read: read.query_name):
                primary = [r for r in group if not (r.is_unmapped or r.is_secondary or r.is_supplementary)]
                if paired:
                    if any(not r.is_paired for r in primary):
                        raise InputError(f"Unpaired alignment in paired input: {stage.alignment}, {name}")
                    primary = [r for r in primary if r.is_proper_pair and not r.mate_is_unmapped]
                    if not primary:
                        continue
                    if len(primary) != 2 or sum(r.is_read1 for r in primary) != 1 or sum(r.is_read2 for r in primary) != 1:
                        raise InputError(f"Incomplete/duplicate proper mates: {stage.alignment}, {name}")
                    observations = [sorted(primary, key=lambda r: not r.is_read1)]
                else:
                    if any(r.is_paired for r in primary):
                        raise InputError(f"Paired alignments require read_layout: paired: {stage.alignment}")
                    observations = [[r] for r in primary]
                for reads in observations:
                    first = reads[0]
                    if any(r.reference_start < 0 or r.reference_end is None or r.reference_end > lengths[r.reference_name] for r in reads):
                        raise InputError(f"Alignment exceeds reference bounds: {stage.alignment}, {name}")
                    if any(r.reference_name != first.reference_name for r in reads):
                        raise InputError(f"Mates on different references: {stage.alignment}, {name}")
                    selected_strand = None if stage.strand == "unstranded" else ("-" if first.is_reverse else "+")
                    if stage.strand == "reverse":
                        selected_strand = "+" if selected_strand == "-" else "-"
                    blocks = [b for r in reads for b in r.get_blocks() if b[0] < b[1]]
                    if not blocks:
                        raise InputError(f"Mapped alignment has no aligned blocks: {stage.alignment}, {name}")
                    cls, features = index.classify(first.reference_name, blocks, selected_strand)
                    totals[cls] += 1
                    if len(features) == 1:
                        fid = next(iter(features))
                        counts[fid] += 1
                        if positions is not None:
                            position = index.relative_position(fid, first.reference_name, blocks)
                            if position is not None:
                                positions.setdefault(fid, np.zeros(100, dtype=np.int64))[min(99, int(position*100))] += 1
                    else:
                        ambiguous[cls] += 1
                    if assignment_writer:
                        assignment_writer.writerow([sample_id, stage.kind, name, cls, ";".join(sorted(features)),
                                                    "assigned" if len(features) == 1 else "ambiguous_feature"])
    if stage.expected:
        expected = read_json(stage.expected)
        actual = {"rRNA": totals["rRNA"], "CDS": totals["CDS"],
                  "other_RNA": sum(totals[c] for c in SUBTYPES), "intergenic_unannotated": totals["intergenic_unannotated"]}
        if expected.get("aligned") != sum(totals.values()) or expected.get("categories") != actual:
            raise InputError(f"Reconstructed categories disagree with {stage.expected}; check strand, GFF version, and stage alignments.")
    return totals, counts, ambiguous


def prepare_categories(sources, selected, config, warnings):
    options = settings(config)
    if config.get("organism", "bacteria") != "bacteria" or not options.get("enabled", True):
        return None
    try:
        import pysam
    except ImportError as exc:
        raise InputError("RNA category analysis requires pysam. Install requirements.txt or update the supplied conda environment.") from exc

    data = CategoryData()
    if options.get("write_read_assignments", False):
        data.assignments = tempfile.TemporaryFile(mode="w+", newline="", encoding="utf-8")
    writer = csv.writer(data.assignments) if data.assignments else None
    if writer:
        writer.writerow(["sample_id", "stage", "read_name", "rna_class", "feature_ids", "assignment_status"])
    maps = read_gene_maps(config)
    cache = {}
    used_alignments = {}
    specs = {str(s["id"]): s for s in config["inputs"]}
    for source in sources:
        samples = [k.split("::", 1)[1] for k in selected if k.startswith(source.id + "::")]
        if not samples:
            continue
        try:
            stages = discover_stages(source, specs[source.id], config, samples)
            seen = set()
            kinds_by_sample = defaultdict(set)
            for stage in stages:
                key = (stage.sample, stage.kind)
                if key in seen:
                    raise InputError(f"{source.id}: duplicated category sample/stage: {key}")
                seen.add(key)
                kinds_by_sample[stage.sample].add(stage.kind)
                alignment = stage.alignment.resolve()
                if alignment in used_alignments:
                    raise InputError(f"The same category alignment is imported more than once: {alignment} ({used_alignments[alignment]} and {source.id}::{stage.sample}, {stage.kind})")
                used_alignments[alignment] = f"{source.id}::{stage.sample}, {stage.kind}"
            if set(kinds_by_sample) != set(samples) or len({tuple(sorted(v)) for v in kinds_by_sample.values()}) != 1:
                raise MissingCategoryInput("Category stages must cover every selected sample with the same RNA/DNA stages.")
            for sample in samples:
                recorded = {st["kind"] for entry in source.metadata.get("samples", []) if entry["sample"] == sample
                            for st in entry.get("stages", []) if st.get("organism") == "bacteria"}
                recorded.update(t.stage for t in source.tables if t.stage in {"rna", "dna"})
                if recorded - kinds_by_sample[sample]:
                    raise MissingCategoryInput(f"{sample}: category_stages omit recorded bacterial stages: {sorted(recorded - kinds_by_sample[sample])}")
            values = {sample: Counter() for sample in samples}
            total_counts = {sample: Counter() for sample in samples}
            ambiguity = {sample: Counter() for sample in samples}
            classes, intervals = {}, defaultdict(list)
            signatures = {}
            for stage in stages:
                with pysam.AlignmentFile(str(stage.alignment), "rb" if stage.alignment.suffix == ".bam" else "r") as handle:
                    references = dict(zip(handle.references, handle.lengths))
                cache_key = (stage.gff, tuple(sorted(references.items())), stage.kind, stage.strand, options["intergenic_bin_bp"], source.id)
                if cache_key not in cache:
                    cache[cache_key] = RNAIndex(stage.gff, references, stage.kind, stage.strand,
                                              options["intergenic_bin_bp"], maps.get(source.id))
                index = cache[cache_key]
                signature = (stage.gff, index.reference_signature, stage.strand)
                if stage.kind in signatures and signatures[stage.kind] != signature:
                    raise InputError(f"{source.id}: samples have different GFFs, references or strand settings for {stage.kind}.")
                signatures[stage.kind] = signature
                key = f"{source.id}::{stage.sample}"
                positions = {} if config.get("analysis", {}).get("import_positional_profiles", True) else None
                totals, counts, ambiguous = classify_stage(stage, source.layout == "paired", index, writer, key, positions)
                if positions is not None:
                    data.positions[(key, stage.kind)] = positions
                values[stage.sample].update(counts)
                total_counts[stage.sample].update(totals)
                ambiguity[stage.sample].update(ambiguous)
                classes.update(index.features)
                for fid, spans in index.intervals.items():
                    intervals[fid].extend(spans)
                data.provenance.append({"sample_id": key, "stage": stage.kind, "strand": stage.strand,
                    "alignment": str(stage.alignment), "alignment_sha256": file_sha256(stage.alignment),
                    "gff": str(stage.gff), "gff_sha256": file_sha256(stage.gff),
                    "aligned_observations": sum(totals.values()), "ambiguous_feature_observations": sum(ambiguous.values())})
            frame = pd.DataFrame({s: [values[s][f] for f in sorted(classes)] for s in samples},
                                 index=pd.Index(sorted(classes), name="gene_id"), dtype=np.int64)
            category_source = replace(source, counts=frame, intervals={f: merge_intervals(v) for f, v in intervals.items()})
            data.sources.append(category_source)
            data.classes[source.id] = classes
            for sample in samples:
                for cat in options["categories"]:
                    data.totals.append({"sample_id": f"{source.id}::{sample}", "category": cat, "status": "reconstructed",
                        "count": sum(n for c, n in total_counts[sample].items() if member(c, cat)),
                        "ambiguous_feature_observations": sum(n for c, n in ambiguity[sample].items() if member(c, cat)),
                        "bacteria_aligned_count": sum(total_counts[sample].values()),
                        "observation_unit": source.observation_unit, "provenance": "retained_alignments"})
        except MissingCategoryInput as exc:
            data.missing[source.id] = str(exc)
            warnings.append(f"{source.id}: RNA category feature analyses unavailable: {exc}")
            if options.get("require_alignments", False):
                raise InputError(warnings[-1]) from exc
    # A missing category is NA, never a biological zero. Modern disposition
    # tables can still supply main-category counts when alignments were removed.
    disposition = read_disposition([s for s in sources if s.id in data.missing], selected, [])
    for source in sources:
        if source.id not in data.missing:
            continue
        for key in selected:
            if not key.startswith(source.id + "::"):
                continue
            rows = disposition.loc[disposition.sample_id == key] if not disposition.empty else pd.DataFrame()
            if len(rows) > 1:
                raise InputError(f"Duplicate read-disposition row for {key}")
            for cat in options["categories"]:
                # Disposition is read-based, even if expression selects molecules.
                label = LABELS.get(cat, "Bacteria other RNA" if cat == "other_RNA_total" else None)
                value = rows.iloc[0].get(f"{label} (reads)") if label and not rows.empty and source.layer == "reads" else None
                if value is not None and str(value).strip() not in {"", "NA", "nan"}:
                    try:
                        count = int(value)
                    except (ValueError, TypeError) as exc:
                        raise InputError(f"Invalid category count for {key}, {cat}: {value}") from exc
                    if count < 0:
                        raise InputError(f"Negative category count for {key}, {cat}")
                else:
                    count = np.nan
                data.totals.append({"sample_id": key, "category": cat,
                    "status": "totals_only" if pd.notna(count) else "unavailable", "count": count,
                    "bacteria_aligned_count": rows.iloc[0].get("Bacteria aligned (reads)", np.nan) if not rows.empty and source.layer == "reads" else np.nan,
                    "ambiguous_feature_observations": np.nan, "observation_unit": source.observation_unit,
                    "provenance": "upstream_read_disposition" if pd.notna(count) else "unavailable"})
    if not data.missing:
        for cat in options["categories"]:
            data.matrices[cat] = category_matrix(data, cat, selected, config)
    return data


def category_matrix(data, cat, selected, config):
    filtered = [replace(s, counts=s.counts.loc[[f for f in s.counts.index if member(data.classes[s.id][f], cat)]]) for s in data.sources]
    if any(s.counts.empty for s in filtered):
        return None
    category_warnings = []
    category_config = dict(config, analysis={k: v for k, v in config.get("analysis", {}).items() if k != "gene_panel"})
    matrix = build_matrix(filtered, selected, category_config, category_warnings, allow_zero_counts=True)
    matrix[1]["rna_class"] = [data.classes[filtered[0].id][f] for f in matrix[1].index]
    return (*matrix, category_warnings)


def category_positions(data, genes, pairs, directory, minimum):
    """Recompute class-specific profiles; never reuse a CDS-only average."""
    from scipy.spatial.distance import jensenshannon
    import matplotlib.pyplot as plt

    rows, summaries, curves = [], [], {}
    for (sample, stage), features in data.positions.items():
        included = {f: v for f, v in features.items() if f in genes and v.sum() >= minimum}
        if not included:
            continue
        fractions = np.stack([v / v.sum() for v in included.values()]).mean(axis=0)
        curves[(sample, stage)] = fractions
        centers = (np.arange(100) + .5)
        summaries.append({"sample_id": sample, "stage": stage, "features_contributing": len(included),
            "mean_position_percent": float(np.sum(centers * fractions)),
            "first_20pct": float(fractions[:20].sum()), "last_20pct": float(fractions[80:].sum()),
            "midpoint_definition": "median_overlapping_aligned_base_on_spliced_feature"})
        for i, value in enumerate(fractions):
            rows.append({"sample_id": sample, "stage": stage, "bin_start_percent": i,
                         "bin_end_percent": i+1, "mean_fraction": value, "features_contributing": len(included)})
    write_csv(pd.DataFrame(rows), directory / "category_metagene_profiles.csv")
    write_csv(pd.DataFrame(summaries), directory / "positional_bias_summary.csv")
    comparisons = []
    for pair in pairs:
        stages = sorted({stage for sample, stage in curves if sample == pair["reference"]} &
                        {stage for sample, stage in curves if sample == pair["query"]})
        if not stages:
            comparisons.append(dict(pair, status="unavailable", reason="No matching stage with sufficient unique feature assignments."))
        for stage in stages:
            a, b = curves[(pair["reference"], stage)], curves[(pair["query"], stage)]
            comparisons.append(dict(pair, stage=stage, status="descriptive_only",
                midpoint_js_distance_base2=float(jensenshannon(a, b, base=2)),
                query_minus_reference_mean_position=float(np.sum((b-a)*(np.arange(100)+.5))),
                same_expressed_feature_panel_verified=False))
    write_csv(pd.DataFrame(comparisons), directory / "positional_pairwise.csv")
    for stage in sorted({stage for sample, stage in curves}):
        fig, ax = plt.subplots(figsize=(9, 6))
        for (sample, kind), fractions in curves.items():
            if kind == stage:
                ax.plot(np.arange(100)+.5, fractions, label=sample)
        ax.set(xlabel="Spliced feature position (5′ → 3′, %)", ylabel="Mean per-feature fraction",
               title=f"Category aligned-base positions: {stage} stage")
        if len(curves) <= 15:
            ax.legend(fontsize=6)
        save(fig, directory / f"metagene_{stage}.png")


def analyze_categories(data, selected, config, pairs, groups, plans, metadata, fasta, threads, output, warnings):
    if data is None:
        return []
    destination = output / "rna_categories"
    destination.mkdir()
    options = settings(config)
    totals = pd.DataFrame(data.totals)
    denominators = pd.to_numeric(totals.bacteria_aligned_count, errors="raise")
    totals["fraction_of_bacteria_aligned"] = totals["count"] / denominators.where(denominators > 0)
    write_csv(totals, output / "rna_category_counts.csv")
    write_csv(pd.DataFrame(data.provenance), destination / "alignment_provenance.csv")
    if data.assignments:
        import shutil
        data.assignments.seek(0)
        with (destination / "read_assignments.csv").open("w", encoding="utf-8", newline="") as handle:
            shutil.copyfileobj(data.assignments, handle)
        data.assignments.close()
    summaries = []
    for cat in options["categories"]:
        print(f"[RNA category] {cat}", flush=True)
        directory = destination / cat
        directory.mkdir()
        part = totals.loc[totals.category == cat].copy()
        write_csv(part, destination / f"{cat}.csv")
        # Per-category total comparisons are meaningful even without gene-level
        # data. Aligned/non-rRNA/other-total are nested, not disjoint components.
        by_sample = part.set_index("sample_id")
        total_pairs = []
        for pair in pairs:
            a, b = by_sample.loc[pair["reference"]], by_sample.loc[pair["query"]]
            valid = pd.notna(a["count"]) and pd.notna(b["count"])
            total_pairs.append(dict(pair, category=cat, reference_count=a["count"], query_count=b["count"],
                reference_fraction_of_bacteria_aligned=a.fraction_of_bacteria_aligned,
                query_fraction_of_bacteria_aligned=b.fraction_of_bacteria_aligned,
                status="available" if valid else "unavailable", normalization="raw_category_observations",
                log2_ratio_count_plus_one=float(np.log2((b["count"]+1)/(a["count"]+1))) if valid else np.nan))
        write_csv(pd.DataFrame(total_pairs), directory / "category_total_comparisons.csv")
        summary = {"category": cat, "results_directory": f"rna_categories/{cat}",
                   "status": "unavailable", "reason": "", "features": 0}
        if data.missing:
            summary["reason"] = "; ".join(f"{sid}: {why}" for sid, why in data.missing.items())
            summaries.append(summary)
            continue
        matrix = data.matrices[cat] if cat in data.matrices else category_matrix(data, cat, selected, config)
        if matrix is None:
            summary.update(status="no_features", reason="This category has no annotated features in at least one selected source.")
            summaries.append(summary)
            continue
        # CDS panels/control genes cannot implicitly restrict a non-CDS panel.
        category_config = dict(config, analysis={k: v for k, v in config.get("analysis", {}).items() if k != "gene_panel"})
        counts, annotation, genomic, audit, retained, category_warnings = matrix
        analysis_metadata = metadata.assign(native_counting_policy=metadata.counting_policy,
            counting_policy="exclusive_category_unique_feature", rna_category=cat)
        if fasta:
            annotation["gc_percent"] = pd.Series(gc_from_fasta(fasta, genomic))
        write_csv(counts, directory / "counts_raw.csv", index=True)
        write_csv(annotation, directory / "gene_annotation.csv", index=True)
        write_csv(audit, directory / "gene_set_audit.csv")
        write_csv(retained, directory / "retained_assignments.csv")
        write_csv(analysis_metadata, directory / "sample_metadata.csv", index=True)
        positive = counts.columns[counts.sum(axis=0) > 0].tolist()
        excluded = [s for s in selected if s not in positive]
        pair_status = [dict(p, status="completed" if p["reference"] in positive and p["query"] in positive else "zero_category_counts") for p in pairs]
        write_csv(pd.DataFrame(pair_status), directory / "comparison_status.csv")
        group_status = [dict(group=g["name"], status="completed" if set(g["samples"]) <= set(positive) else "zero_category_counts") for g in groups]
        write_csv(pd.DataFrame(group_status), directory / "reproducibility_status.csv")
        valid_pairs = [p for p in pairs if p["reference"] in positive and p["query"] in positive]
        valid_groups = [g for g in groups if set(g["samples"]) <= set(positive)]
        if positive:
            active = counts[positive]
            analysis_settings = category_config.get("analysis", {})
            cpm, log_cpm, _ = pairwise_analysis(active, annotation, valid_pairs, analysis_settings, directory)
            # Keep zero-count libraries in raw tables and mark their normalized
            # values NA; a zero category has no defined within-category CPM.
            write_csv(counts, directory / "counts_raw.csv", index=True)
            write_csv(cpm.reindex(columns=selected), directory / "counts_cpm.csv", index=True)
            write_csv(log_cpm.reindex(columns=selected), directory / "counts_log2cpm.csv", index=True)
            global_qc(active, cpm, log_cpm, analysis_metadata.loc[positive], analysis_settings, directory, category_warnings)
            depth_analysis(active, directory)
            reproducibility_analysis(active, cpm, log_cpm, valid_groups, analysis_settings, directory)
        else:
            empty_normalized = counts.astype(float) * np.nan
            write_csv(empty_normalized, directory / "counts_cpm.csv", index=True)
            write_csv(empty_normalized, directory / "counts_log2cpm.csv", index=True)
        valid_plans, skipped_de = [], []
        for plan in plans:
            try:
                count_filter(counts, plan)
            except InputError as exc:
                skipped_de.append({"analysis": plan.name, "status": "not_testable", "reason": str(exc)})
                continue
            if plan.settings.get("control_genes"):
                skipped_de.append({"analysis": plan.name, "status": "not_testable",
                    "reason": "Root CDS control_genes cannot normalize every RNA category. Use control_genes: null for category DE."})
                continue
            plan_metadata = plan.metadata.assign(native_counting_policy=plan.metadata.counting_policy,
                counting_policy="exclusive_category_unique_feature", rna_category=cat)
            valid_plans.append(replace(plan, metadata=plan_metadata))
        run_de(counts, annotation, valid_plans, category_config, threads, directory, category_warnings)
        if skipped_de:
            status_path = directory / "differential_expression_status.csv"
            existing = pd.read_csv(status_path) if valid_plans else pd.DataFrame()
            write_csv(pd.concat([existing, pd.DataFrame(skipped_de)], ignore_index=True), status_path)
        category_warnings.append("Feature counts use unique assignment within the winning RNA class; ambiguous overlaps are excluded from feature matrices but retained in category totals. Intergenic features are reference-anchored annotation-free bins. CPM is normalized separately within this category.")
        if config.get("analysis", {}).get("import_positional_profiles", True):
            category_positions(data, set(counts.index), pairs, directory, int(config.get("analysis", {}).get("min_count", 10)))
            category_warnings.append("Category positional profiles use the median overlapping aligned base on the spliced annotated feature, not fragment spans, per-base coverage, or upstream CDS-only averages. RNA/DNA stages remain separate. Expressed features contributing to each average may differ.")
        if config.get("analysis", {}).get("gene_panel"):
            category_warnings.append("The root CDS gene_panel is not applied to RNA category panels.")
        (directory / "NOTES.txt").write_text("\n".join(dict.fromkeys(category_warnings)) + "\n", encoding="utf-8")
        summary.update(status="completed" if not excluded else "partial_zero_counts", features=len(counts),
                       reason="Zero category assignments in: " + ", ".join(excluded) if excluded else "")
        summaries.append(summary)
    result = pd.DataFrame(summaries)
    write_csv(result, output / "rna_category_analysis_status.csv")
    warnings.append("RNA categories are analyzed separately under rna_categories/. Aligned and non-rRNA are nested totals; other_RNA_total repeats bacteria_other_RNA. Category feature assignments use aligned blocks and strand-aware annotation; the original root CDS results retain native featureCounts semantics.")
    if any(r["status"] != "completed" for r in summaries):
        warnings.append("Some RNA category feature analyses are unavailable, empty, or include zero-count libraries; consult rna_category_analysis_status.csv and each category's comparison_status.csv.")
    return summaries
