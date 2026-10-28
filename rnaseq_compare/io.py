"""Validated adapters for the three ScienceAndCats pipelines.

Only explicit published count matrices are imported. We never recursively glob
all historical worker outputs, sum read and molecule tables, or infer zeros for
genes missing from another run's annotation.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

FC_HEADER = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
PIPELINES = {"shortread", "agilent_pe", "longread"}
CSV_FIELD_LIMIT = 64 * 1024 * 1024
csv.field_size_limit(CSV_FIELD_LIMIT)


class InputError(ValueError):
    """An input or study design cannot be interpreted safely."""


class UniqueLoader(yaml.SafeLoader):
    """Reject duplicate YAML keys instead of silently replacing settings."""


def _unique_mapping(loader, node, deep=False):
    result = {}
    for key_node, value_node in node.value:
        key = loader.construct_object(key_node, deep=deep)
        if key in result:
            raise InputError(f"Duplicate YAML key: {key!r}")
        result[key] = loader.construct_object(value_node, deep=deep)
    return result


UniqueLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _unique_mapping)


def load_config(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        config = yaml.load(handle, Loader=UniqueLoader)
    if not isinstance(config, dict) or config.get("version") != 1:
        raise InputError("The YAML must be a mapping with version: 1.")
    if not isinstance(config.get("inputs"), list) or not config["inputs"]:
        raise InputError("inputs must be a nonempty list of output-directory definitions.")
    if config.get("organism", "bacteria") not in {"bacteria", "host"}:
        raise InputError("organism must be bacteria or host; run them separately.")
    allowed = {"version", "output_dir", "organism", "inputs", "samples", "comparison_groups",
               "reproducibility_groups", "analysis", "reference", "differential_expression"}
    unknown = set(config) - allowed
    if unknown:
        raise InputError(f"Unknown top-level YAML keys: {sorted(unknown)}")
    config["_config_path"] = str(path.resolve())
    config["_base"] = str(path.resolve().parent)
    return config


def config_path(value: str, base: Path) -> Path:
    path = Path(str(value)).expanduser()
    return (path if path.is_absolute() else base / path).resolve()


def read_json(path: Path) -> dict:
    try:
        with path.open(encoding="utf-8") as handle:
            result = json.load(handle)
    except (ValueError, OSError) as exc:
        raise InputError(f"Cannot read JSON {path}: {exc}") from exc
    if not isinstance(result, dict):
        raise InputError(f"Expected a JSON object in {path}")
    return result


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fastq_stem(value: str) -> str:
    return re.sub(r"\.(?:fastq|fq)(?:\.gz)?$", "", Path(str(value)).name, flags=re.I)


def count_sample_name(value: str, organism: str) -> str:
    name = Path(value.replace("\\", "/")).name
    name = re.sub(r"(?:\.sorted|\.molecules)?\.(?:sam|bam|cram)$", "", name, flags=re.I)
    prefix = organism.upper() + "_"
    if name.startswith(prefix):
        name = name[len(prefix):]
    return re.sub(r"_\d+bp$", "", name)


def merge_intervals(intervals):
    """Return nonoverlapping, 1-based inclusive intervals, keyed by contig/strand."""
    groups = {}
    for chrom, start, end, strand in intervals:
        groups.setdefault((chrom, strand), []).append((int(start), int(end)))
    merged = []
    for (chrom, strand), spans in sorted(groups.items()):
        current = []
        for start, end in sorted(set(spans)):
            if current and start <= current[-1][1] + 1:
                current[-1][1] = max(current[-1][1], end)
            else:
                current.append([start, end])
        merged.extend((chrom, start, end, strand) for start, end in current)
    return tuple(merged)


def genomic_intervals(row: list[str], stage: str):
    # RNA-stage CDS coordinates are not genomic gene lengths. Long-read combined
    # tables may list the same locus on both RNA__ targets and genomic contigs.
    if stage == "rna":
        return ()
    chroms, starts, ends, strands = [row[i].split(";") for i in (1, 2, 3, 4)]
    if not (len(chroms) == len(starts) == len(ends) == len(strands)):
        raise InputError(f"Mismatched featureCounts coordinate lists for {row[0]!r}")
    spans = []
    for chrom, start, end, strand in zip(chroms, starts, ends, strands):
        if chrom.startswith("RNA__"):
            continue
        try:
            lo, hi = int(start), int(end)
        except ValueError as exc:
            raise InputError(f"Invalid coordinates for {row[0]!r}") from exc
        if lo < 1 or hi < lo:
            raise InputError(f"Invalid coordinates for {row[0]!r}: {lo}-{hi}")
        spans.append((chrom, lo, hi, strand))
    return merge_intervals(spans)


@dataclass
class CountTable:
    path: Path
    stage: str
    counts: pd.DataFrame
    intervals: dict
    columns: dict  # canonical sample -> original featureCounts alignment column
    comments: list[str]


def read_count_table(path: Path, stage: str, organism: str, gene_map: dict | None = None) -> CountTable:
    comments = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        first = handle.readline()
        while first and (first.startswith("#") or not first.strip()):
            comments.append(first.rstrip())
            first = handle.readline()
        if not first:
            raise InputError(f"Empty count table: {path}")
        delimiter = "\t" if "\t" in first else ","
        import itertools
        rows = csv.reader(itertools.chain([first], handle), delimiter=delimiter)
        header = next(rows)
        if header[:6] != FC_HEADER or len(header) < 7:
            raise InputError(f"Not a gene-level featureCounts matrix: {path}. Use the .txt table, not its .summary report.")
        sample_names = [count_sample_name(x, organism) for x in header[6:]]
        if len(sample_names) != len(set(sample_names)) or any(not x for x in sample_names):
            raise InputError(f"Duplicate/empty resolved sample columns in {path}")
        values, genes, intervals, seen = [], [], {}, set()
        mapping = gene_map or {}
        for number, row in enumerate(rows, 2):
            if not row:
                continue
            if len(row) != len(header):
                raise InputError(f"Wrong column count in {path}, record {number}; prefer the original tab-delimited .txt.")
            gene = mapping.get(row[0], row[0])
            if not gene or gene in seen:
                raise InputError(f"Duplicate/empty gene ID after mapping in {path}: {gene!r}. Many-to-one ID mappings are not supported.")
            seen.add(gene)
            if any(re.fullmatch(r"[0-9]+", x.strip()) is None for x in row[6:]):
                raise InputError(f"Noninteger, negative, or missing counts in {path}, gene {gene!r}; normalized or fractional counts are not accepted.")
            counts = [int(x) for x in row[6:]]
            if any(x > np.iinfo(np.int64).max for x in counts):
                raise InputError(f"Count exceeds int64 in {path}: {gene}")
            genes.append(gene)
            values.append(counts)
            intervals[gene] = genomic_intervals(row, stage)
    frame = pd.DataFrame(values, index=pd.Index(genes, name="gene_id"), columns=sample_names, dtype=np.int64)
    return CountTable(path, stage, frame, intervals, dict(zip(sample_names, header[6:])), comments)


def table_path(directory: Path, organism: str) -> Path | None:
    stem = f"featurecounts_{organism.upper()}_summary"
    for suffix in (".txt", ".csv"):
        candidate = directory / (stem + suffix)
        if candidate.is_file():
            return candidate
    return None


def select_root(path: Path, organism: str, pipeline: str, layer: str) -> Path:
    def has_counts(root):
        root = root / "molecule_counts" if layer == "molecules" else root
        if organism == "bacteria":
            return table_path(root, organism) is not None
        return any(table_path(root / p, organism) for p in ("bowtie_alignments/host", "bbsplit_alignments/host", "host"))
    if not path.is_dir():
        raise InputError(f"Input output directory does not exist: {path}")
    if has_counts(path):
        if pipeline == "longread" and (path / "run.json").exists() and not (path / "aggregation.completed.json").exists():
            raise InputError(f"Long-read Slurm aggregation is not marked complete: {path}")
        return path
    runs = path / "slurm_runs"
    candidates = [p for p in sorted(runs.iterdir()) if p.is_dir() and has_counts(p)] if runs.is_dir() else []
    if len(candidates) == 1:
        return select_root(candidates[0], organism, pipeline, layer)
    if len(candidates) > 1:
        raise InputError("Multiple run directories contain count tables; set the exact path instead of choosing the newest:\n" + "\n".join(map(str, candidates)))
    raise InputError(f"No published {organism} CDS count table at {path}. Finalize the upstream run or point to its exact run directory. count_layer={layer}.")


@dataclass
class Source:
    id: str
    root: Path
    pipeline: str
    layer: str
    layout: str
    library_prep: str
    tables: list[CountTable]
    counts: pd.DataFrame
    intervals: dict
    metadata: dict = field(default_factory=dict)

    @property
    def counting_policy(self):
        return "multiple_overlapping_genes" if self.pipeline == "longread" else "unique_gene_assignment"

    @property
    def observation_unit(self):
        if self.layer == "molecules":
            return "barcode_deduplicated_pairs"
        return "read_pairs" if self.layout == "paired" else "reads"


def read_gene_maps(config) -> dict:
    value = config.get("reference", {}).get("gene_id_map")
    if not value:
        return {}
    path = config_path(value, Path(config["_base"]))
    maps = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        required = {"source", "gene_id", "canonical_gene_id"}
        if not required <= set(reader.fieldnames or []):
            raise InputError(f"Gene mapping CSV needs {sorted(required)}: {path}")
        for row in reader:
            source, old, new = (row[k].strip() for k in ("source", "gene_id", "canonical_gene_id"))
            if not source or not old or not new:
                raise InputError(f"Empty ID in {path}")
            mapping = maps.setdefault(source, {})
            if old in mapping or new in mapping.values():
                raise InputError(f"Gene mapping must be one-to-one within a source: {source}, {old}, {new}")
            mapping[old] = new
    return maps


def load_sources(config: dict, warnings: list[str]) -> list[Source]:
    sources, ids, used_tables = [], set(), set()
    used_alignments = {}
    organism = config.get("organism", "bacteria")
    maps = read_gene_maps(config)
    for spec in config["inputs"]:
        if not isinstance(spec, dict) or not {"id", "path", "pipeline"} <= spec.keys():
            raise InputError("Each input needs id, path, and pipeline.")
        sid = str(spec["id"])
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", sid) or sid in ids:
            raise InputError(f"Invalid or duplicate input id: {sid}")
        ids.add(sid)
        pipeline, layer = spec["pipeline"], spec.get("count_layer", "reads")
        if pipeline not in PIPELINES or layer not in {"reads", "molecules"}:
            raise InputError(f"Invalid pipeline/count_layer for {sid}")
        if layer == "molecules" and pipeline != "agilent_pe":
            raise InputError("Only agilent_pe provides the supported molecule_counts tree.")
        root = select_root(config_path(spec["path"], Path(config["_base"])), organism, pipeline, layer)
        manifest = read_json(root / "run_manifest.json") if (root / "run_manifest.json").is_file() else {}
        layout = spec.get("read_layout", manifest.get("read_layout", "paired" if pipeline == "agilent_pe" else "single" if pipeline == "longread" else "unknown"))
        if layout not in {"single", "paired", "unknown"}:
            raise InputError(f"Invalid read_layout for {sid}: {layout}")
        if pipeline == "longread" and layout != "single":
            raise InputError("The long-read pipeline treats each FASTQ record independently, not as a paired read.")
        if layout == "unknown":
            warnings.append(f"{sid}: read layout is not recorded; exact sample names work, but paired FASTQ-name aliases require read_layout: paired.")
        count_root = root / "molecule_counts" if layer == "molecules" else root
        if pipeline == "longread":
            dirs = [(count_root if organism == "bacteria" else count_root / "bbsplit_alignments/host", "combined")]
            if organism == "host" and not table_path(dirs[0][0], organism):
                dirs = [(count_root / "host", "combined")]
        else:
            dirs = [(count_root if organism == "bacteria" else count_root / "bowtie_alignments/host", "dna"),
                    (count_root / f"bowtie_alignments/{organism}/rna", "rna")]
        tables = []
        for directory, stage in dirs:
            path = table_path(directory, organism)
            if path is None:
                continue
            resolved = path.resolve()
            if resolved in used_tables:
                raise InputError(f"The same count matrix is imported more than once: {path}")
            used_tables.add(resolved)
            table = read_count_table(path, stage, organism, maps.get(sid))
            for recorded in table.columns.values():
                alignment = Path(recorded)
                alignment = str((alignment if alignment.is_absolute() else root / alignment).resolve())
                previous = used_alignments.get(alignment)
                if previous is not None and previous != sid:
                    raise InputError(f"The same recorded alignment is imported through two inputs ({previous}, {sid}): {alignment}")
                used_alignments[alignment] = sid
            tables.append(table)
        if not tables or not any(len(t.counts) for t in tables):
            raise InputError(f"No CDS features available for {sid}")
        samples = list(tables[0].counts.columns)
        for table in tables[1:]:
            if set(table.counts.columns) != set(samples):
                raise InputError(f"RNA/DNA stage sample columns differ for {sid}; incomplete/stale outputs cannot be combined.")
        if manifest.get("samples"):
            expected = {entry["sample"] for entry in manifest["samples"]}
            if expected != set(samples):
                raise InputError(f"Count table samples disagree with run_manifest.json for {sid}.")
            # If a normalized RNA annotation declares CDS, its table must exist.
            # Noncoding-only RNA stages correctly have no CDS matrix.
            refs_path = root / "reference_manifest.json"
            if refs_path.is_file() and not any(t.stage == "rna" for t in tables):
                refs = read_json(refs_path)
                rna_gff = ((refs.get(organism) or {}).get("rna") or {}).get("gff")
                if rna_gff:
                    gff = Path(rna_gff)
                    if not gff.is_absolute():
                        gff = root / gff
                    if gff.is_file():
                        opener = gzip.open if gff.name.endswith(".gz") else open
                        with opener(gff, "rt") as handle:
                            if any(not line.startswith("#") and len(line.split("\t")) == 9 and line.split("\t")[2] == "CDS" for line in handle):
                                raise InputError(f"{sid}: RNA reference has CDS annotations, but its CDS count table is missing.")
                    else:
                        warnings.append(f"{sid}: RNA CDS-table completeness could not be checked because the referenced GFF is unavailable: {gff}")
        genes = sorted(set().union(*(set(t.counts.index) for t in tables)))
        counts = pd.DataFrame(0, index=pd.Index(genes, name="gene_id"), columns=samples, dtype=np.int64)
        intervals = {}
        for table in tables:
            counts.loc[table.counts.index, samples] += table.counts.loc[:, samples]
            for gene, spans in table.intervals.items():
                if spans:
                    if gene in intervals and intervals[gene] != spans:
                        raise InputError(f"Inconsistent genomic coordinates across stages for {sid}, {gene}")
                    intervals[gene] = spans
        if pipeline == "longread":
            warnings.append(f"{sid}: native long-read CDS counts use -O (multiple overlapping genes); reads >6000 bp are excluded by the current upstream pipeline. Count units, mapping, and positional definitions differ from short reads.")
        sources.append(Source(sid, root, pipeline, layer, layout, str(spec.get("library_prep", sid)), tables, counts, intervals, manifest))
    unknown_maps = set(maps) - ids
    if unknown_maps:
        raise InputError(f"Unknown source IDs in gene_id_map: {sorted(unknown_maps)}")
    return sources


class Catalog:
    def __init__(self, sources: list[Source], config: dict):
        self.sources = {s.id: s for s in sources}
        self.records, self.aliases = {}, {}
        assigned_metadata = {}
        for source in sources:
            for sample in source.counts.columns:
                key = f"{source.id}::{sample}"
                self.records[key] = {"sample_id": key, "sample": sample, "source": source.id,
                    "pipeline": source.pipeline, "library_prep": source.library_prep,
                    "count_layer": source.layer, "observation_unit": source.observation_unit,
                    "counting_policy": source.counting_policy, "biological_sample": ""}
        for alias, spec in (config.get("samples") or {}).items():
            if not isinstance(spec, dict):
                raise InputError(f"samples.{alias} must be a metadata mapping.")
            if "source" in spec:
                key = self.resolve(f"{spec['source']}::{spec.get('sample', alias)}")
            else:
                key = self.resolve(str(spec.get("sample", alias)))
            if alias in self.aliases:
                raise InputError(f"Duplicate sample alias: {alias}")
            self.aliases[str(alias)] = key
            protected = {"sample_id", "source", "sample", "pipeline", "count_layer", "observation_unit", "counting_policy"}
            for name, value in spec.items():
                if name not in protected:
                    if isinstance(value, (list, dict)):
                        raise InputError(f"Sample metadata values must be scalar: {alias}.{name}")
                    normalized = "" if value is None else str(value)
                    prior = assigned_metadata.setdefault(key, {})
                    if name in prior and prior[name] != normalized:
                        raise InputError(f"Conflicting metadata for {key}: {name}")
                    prior[name] = normalized
                    self.records[key][name] = normalized

    def resolve(self, value: str) -> str:
        value = str(value)
        if value in self.aliases:
            return self.aliases[value]
        if value in self.records:
            return value
        source_id, raw = value.split("::", 1) if "::" in value else (None, value)
        stem = fastq_stem(raw)
        candidates = [key for key, rec in self.records.items() if (source_id is None or rec["source"] == source_id) and rec["sample"] == stem]
        if not candidates:
            # ONLY paired sources can normalize an R1/R2 FASTQ suffix. Long-read
            # _R1 and _R2 are independent sample names and must not be collapsed.
            paired_stem = re.sub(r"_R[12](?:_\d{3})?$", "", stem)
            if paired_stem != stem:
                candidates = [key for key, rec in self.records.items()
                              if (source_id is None or rec["source"] == source_id)
                              and self.sources[rec["source"]].layout == "paired" and rec["sample"] == paired_stem]
        if len(candidates) != 1:
            if candidates:
                raise InputError(f"Ambiguous sample {value!r}; use one of {candidates} or define a samples alias.")
            raise InputError(f"Sample {value!r} was not found. Run --list-samples to inspect the imported names.")
        return candidates[0]

    def metadata(self, selected):
        return pd.DataFrame([self.records[k] for k in selected]).set_index("sample_id")


def comparison_pairs(config, catalog):
    groups = config.get("comparison_groups") or []
    if isinstance(groups, dict):
        entries = list(groups.items())
    elif isinstance(groups, list):
        entries = [(f"group_{i:02d}", group) for i, group in enumerate(groups, 1)]
    else:
        raise InputError("comparison_groups must be a mapping of group names to lists, or a list of lists.")
    pairs = []
    for name, values in entries:
        if not isinstance(values, list) or len(values) < 2:
            raise InputError(f"Comparison group {name} must contain a reference followed by at least one query.")
        ids = [catalog.resolve(v) for v in values]
        if len(ids) != len(set(ids)):
            raise InputError(f"Duplicate sample (possibly a second alias) in comparison group {name}")
        pairs.extend({"group": str(name), "reference": ids[0], "query": query} for query in ids[1:])
    return pairs


def reproducibility_groups(config, catalog):
    groups = []
    for spec in config.get("reproducibility_groups") or []:
        if not isinstance(spec, dict) or not {"name", "samples", "kind"} <= spec.keys():
            raise InputError("Each reproducibility group needs name, samples, and kind (technical or biological).")
        ids = [catalog.resolve(v) for v in spec["samples"]]
        if len(ids) < 2 or len(ids) != len(set(ids)):
            raise InputError(f"Reproducibility group {spec['name']} requires at least two distinct libraries.")
        if spec["kind"] not in {"technical", "biological"}:
            raise InputError("Reproducibility kind must be technical or biological.")
        records = [catalog.records[k] for k in ids]
        if len({r["library_prep"] for r in records}) != 1:
            raise InputError(f"Reproducibility group {spec['name']} mixes library_prep values; use comparison_groups for method agreement.")
        bios = [r.get("biological_sample", "") for r in records]
        if spec["kind"] == "technical" and len({b for b in bios if b}) > 1:
            raise InputError(f"Technical group {spec['name']} contains different biological samples.")
        if spec["kind"] == "biological" and (any(not b for b in bios) or len(set(bios)) != len(bios)):
            raise InputError(f"Biological reproducibility group {spec['name']} needs distinct biological_sample IDs.")
        groups.append(dict(spec, samples=ids))
    names = [g["name"] for g in groups]
    if len(names) != len(set(names)):
        raise InputError("Reproducibility group names must be unique.")
    return groups


def read_gene_list(path: Path) -> list[str]:
    genes = []
    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            text = line.strip()
            if text and not text.startswith("#") and text.lower() not in {"gene_id", "geneid"}:
                if any(ch in text for ch in (",", "\t")):
                    raise InputError(f"Gene lists must contain one canonical ID per line: {path}")
                genes.append(text)
    if not genes or len(genes) != len(set(genes)):
        raise InputError(f"Empty or duplicate IDs in gene list: {path}")
    return genes


def build_matrix(sources, selected, config, warnings):
    used_sources = [s for s in sources if any(k.startswith(s.id + "::") for k in selected)]
    sets = [set(s.counts.index) for s in used_sources]
    common = set.intersection(*sets)
    union = set.union(*sets)
    settings = config.get("analysis", {})
    threshold = float(settings.get("min_shared_gene_fraction", 0.8))
    if not 0 < threshold <= 1:
        raise InputError("min_shared_gene_fraction must be in (0, 1].")
    join = settings.get("gene_join", "intersection")
    if join not in {"intersection", "strict"}:
        raise InputError("gene_join must be intersection or strict; union-zero-filling is intentionally unsupported.")
    if join == "strict" and common != union:
        raise InputError("Gene sets differ between outputs under gene_join: strict.")
    audit = []
    for source in used_sources:
        fraction = len(common) / len(source.counts.index)
        if fraction < threshold:
            raise InputError(f"Only {fraction:.1%} of {source.id} genes share exact IDs across outputs. Check genome versions and provide a one-to-one gene_id_map; do not blindly lower the threshold.")
        audit.append({"source": source.id, "annotated_genes": len(source.counts), "shared_genes": len(common), "shared_fraction": fraction})
    if common != union:
        warnings.append(f"Using {len(common)} shared genes out of {len(union)} distinct IDs; absent annotations are excluded, not treated as biological zero expression.")
    if settings.get("gene_panel"):
        panel = set(read_gene_list(config_path(settings["gene_panel"], Path(config["_base"]))))
        missing = panel - common
        if missing:
            raise InputError(f"Requested panel contains {len(missing)} genes absent from the shared annotation, e.g. {sorted(missing)[:8]}")
        common &= panel
        warnings.append(f"Analysis restricted to {len(common)} configured panel genes. CPM denominators and DE normalization use this panel.")
    if not common:
        raise InputError("No shared genes remain.")
    genes = sorted(common)
    source_by_id = {s.id: s for s in sources}
    data = {}
    coverage_rows = []
    for key in selected:
        sid, sample = key.split("::", 1)
        source = source_by_id[sid]
        values = source.counts[sample]
        selected_values = values.loc[genes]
        data[key] = selected_values
        total, retained = int(values.sum()), int(selected_values.sum())
        coverage_rows.append({"sample_id": key, "all_source_gene_assignments": total,
            "retained_gene_assignments": retained, "retained_assignment_fraction": retained / total if total else np.nan})
        if retained == 0:
            raise InputError(f"{key} has zero counts on the analysis gene set. Remove it from the requested comparisons; normalization is undefined.")
        if total and retained / total < float(settings.get("warn_retained_assignment_fraction", 0.8)):
            warnings.append(f"{key}: only {retained/total:.1%} of source gene assignments remain after gene-set/panel selection.")
    matrix = pd.DataFrame(data, index=pd.Index(genes, name="gene_id"), dtype=np.int64)
    genomic, annotation_rows, mismatch = {}, [], []
    for gene in genes:
        known = [(s.id, s.intervals[gene]) for s in used_sources if s.intervals.get(gene)]
        signatures = {x[1] for x in known}
        if len(signatures) > 1:
            mismatch.append(gene)
        spans = known[0][1] if len(signatures) == 1 else ()
        genomic[gene] = spans
        annotation_rows.append({"gene_id": gene, "length_bp": sum(end-start+1 for _, start, end, _ in spans) if spans else np.nan,
            "genomic_coordinates_verified_across_sources": len(known) == len(used_sources) and len(signatures) == 1,
            "coordinate_sources": ";".join(sid for sid, _ in known), "gc_percent": np.nan})
    if mismatch:
        if settings.get("require_matching_genomic_coordinates", True):
            raise InputError(f"Genomic CDS coordinates disagree for {len(mismatch)} shared IDs, e.g. {mismatch[:8]}. Use the same assembly/annotation. An explicit override can retain counts but disables length/GC for conflicting genes.")
        warnings.append(f"{len(mismatch)} genes have conflicting genomic coordinates: length/GC bias values are unavailable for them. Matching IDs do not prove comparable feature definitions.")
    annotation = pd.DataFrame(annotation_rows).set_index("gene_id")
    if not annotation.genomic_coordinates_verified_across_sources.all():
        warnings.append("Some genomic feature coordinates could not be verified across every source; see gene_annotation.csv.")
    return matrix, annotation, genomic, pd.DataFrame(audit), pd.DataFrame(coverage_rows)


def gc_from_fasta(path: Path, genomic: dict) -> dict[str, float]:
    """Stream FASTA and collect GC/ACGT prefix counts only at CDS boundaries.

    Memory scales with annotation boundaries, not genome size. Strand does not
    change GC. Ambiguous bases are excluded from the denominator.
    """
    needed = {}
    for spans in genomic.values():
        for chrom, start, end, _ in spans:
            needed.setdefault(chrom, set()).update((start - 1, end))
    bounds = {chrom: sorted(points) for chrom, points in needed.items()}
    prefix, seen = {}, set()
    chrom, position, gc, acgt, cursor = None, 0, 0, 0, 0
    def finish():
        if chrom in bounds and cursor != len(bounds[chrom]):
            raise InputError(f"FASTA contig {chrom} is shorter than the gene coordinates in the count table.")
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                chrom = line[1:].split()[0]
                if chrom in seen:
                    raise InputError(f"Duplicate FASTA contig: {chrom}")
                seen.add(chrom)
                position, gc, acgt, cursor = 0, 0, 0, 0
                if chrom in bounds:
                    prefix[chrom] = {}
                continue
            if chrom is None:
                raise InputError(f"Sequence before first FASTA header: {path}")
            if chrom not in bounds:
                continue
            sequence = line.upper()
            stop, offset = position + len(sequence), 0
            while cursor < len(bounds[chrom]) and bounds[chrom][cursor] <= stop:
                boundary = bounds[chrom][cursor]
                cut = boundary - position
                section = sequence[offset:cut]
                gc += section.count("G") + section.count("C")
                acgt += sum(section.count(base) for base in "ACGT")
                prefix[chrom][boundary] = (gc, acgt)
                offset = cut
                cursor += 1
            section = sequence[offset:]
            gc += section.count("G") + section.count("C")
            acgt += sum(section.count(base) for base in "ACGT")
            position = stop
    finish()
    missing = set(bounds) - seen
    if missing:
        raise InputError(f"Genomic contig IDs absent from genome_fasta: {sorted(missing)[:10]}")
    result = {}
    for gene, spans in genomic.items():
        total_gc = total_acgt = 0
        for chrom, start, end, _ in spans:
            hi, lo = prefix[chrom][end], prefix[chrom][start-1]
            total_gc += hi[0] - lo[0]
            total_acgt += hi[1] - lo[1]
        result[gene] = 100 * total_gc / total_acgt if total_acgt else np.nan
    return result


def optional_sidecar_paths(source: Source, table: CountTable, sample: str, suffix: str):
    """Resolve recorded paths and recognizable relocated output suffixes only."""
    original = Path(table.columns[sample])
    stem = re.sub(r"(?:\.sorted|\.molecules)?\.(?:sam|bam|cram)$", "", original.name, flags=re.I)
    original = original.with_name(stem + suffix)
    candidates = [original if original.is_absolute() else source.root / original]
    parts = original.parts
    for anchor in ("slurm_runs", "samples", "bowtie_alignments", "bbsplit_alignments"):
        if anchor in parts:
            index = len(parts) - 1 - tuple(reversed(parts)).index(anchor)
            candidates.append(source.root.joinpath(*parts[index:]))
    candidates.append(table.path.parent / original.name)
    found = list(dict.fromkeys(p.resolve() for p in candidates if p.is_file()))
    if len(found) > 1:
        # Multiple paths can point to different copies; never guess which one.
        hashes = {file_sha256(p) for p in found}
        if len(hashes) > 1:
            raise InputError(f"Conflicting positional sidecars for {source.id}::{sample}: {found}")
    return found[0] if found else None


def read_disposition(sources: list[Source], selected: list[str], warnings: list[str]) -> pd.DataFrame:
    frames = []
    for source in sources:
        candidates = [source.root / "output_read_breakdown.csv", source.root / f"BACTERIA_{source.root.name}_read_breakdown.csv"]
        path = next((p for p in candidates if p.is_file()), None)
        if path is None:
            warnings.append(f"{source.id}: read-disposition CSV is unavailable; expression analysis remains available.")
            continue
        # csv.reader rather than pandas' C parser supports 64 MiB source fields.
        with path.open(newline="", encoding="utf-8-sig") as handle:
            rows = list(csv.DictReader(handle))
        for row in rows:
            sample = row.get("Sample", row.get("sample"))
            if sample is None:
                raise InputError(f"No Sample column in {path}")
            key = f"{source.id}::{sample}"
            if key in selected:
                frames.append({"sample_id": key, "source": source.id, "observation_unit": source.observation_unit,
                               "count_layer_for_expression": source.layer,
                               **{k: v for k, v in row.items() if k not in {"Sample", "sample"}}})
    return pd.DataFrame(frames)


def detect_genome_fasta(sources, organism):
    """Locate a DNA-only FASTA through the known upstream reference manifests.

    Uses exact prepared paths, or recognized relocated directory suffixes. Never
    searches arbitrary FASTAs or substitutes a combined RNA/genome reference.
    """
    for source in sources:
        candidates = [source.root / "reference_manifest.json",
                      source.root / "bbsplit_alignments/references/manifest.json",
                      source.root / "shared/bbsplit_alignments/references/manifest.json"]
        for manifest_path in candidates:
            if not manifest_path.is_file():
                continue
            dna = ((read_json(manifest_path).get(organism) or {}).get("dna") or {})
            recorded = []
            if dna.get("fasta"):
                recorded.append(Path(dna["fasta"]))
            if source.pipeline != "longread" and dna.get("gff"):
                recorded.append(Path(dna["gff"]).with_suffix(".fa"))
            for original in recorded:
                paths = [original if original.is_absolute() else manifest_path.parent / original]
                for anchor in (".reference_cache", "bbsplit_alignments", "shared"):
                    if anchor in original.parts:
                        offset = original.parts.index(anchor)
                        paths.append(source.root.joinpath(*original.parts[offset:]))
                paths.append(manifest_path.parent / original.name)
                for path in paths:
                    if path.is_file():
                        return path.resolve(), source.id
    return None, None
