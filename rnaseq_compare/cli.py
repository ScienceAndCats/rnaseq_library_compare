"""Command-line orchestration and reproducible reporting."""
from __future__ import annotations

import argparse
import html
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
import traceback
from datetime import datetime, timezone

import pandas as pd
import yaml

from . import __version__
from .analysis import (depth_analysis, global_qc, pairwise_analysis, positional_profiles,
                       reproducibility_analysis, write_csv)
from .de import prepare_plans, count_filter, run_de
from .io import (Catalog, InputError, build_matrix, comparison_pairs, config_path, file_sha256,
                 gc_from_fasta, load_config, load_sources, read_disposition, reproducibility_groups, detect_genome_fasta)


LIMITATIONS = [
    "Gene expression here means native CDS gene-assignment counts, not transcript isoforms, all ncRNAs, or absolute RNA molecules.",
    "A sample-to-sample log2 ratio is descriptive. A low correlation is not evidence that a library-prep method is worse when biological samples differ.",
    "P-values require independent biological samples. Repeated lanes or library preps from the same RNA do not supply biological replication.",
    "Cross-platform/pipeline differences can include alignment, counting rules, read length, strandedness, capture targets, PCR duplication, and RNA/DNA-stage routing—not only library prep.",
    "CPM is compositional and uses the selected shared-gene/panel assignment total. DESeq2 normalization assumes a suitable stable reference across genes; targeted capture or broad shifts may violate that assumption.",
    "Assignment-depth rarefaction samples gene-count entries, not original reads, molecules, or library complexity. Long-read -O can assign one read to several genes.",
    "Upstream metagene curves describe read/fragment midpoints on potentially different gene panels and coordinate systems. They are not per-base 5′/3′ coverage or a harmonized coverage assay.",
    "No remapping, recounting, duplicate removal, or automatic technical-replicate collapse is performed. No upstream outputs are modified.",
]


def versions():
    result = {"rnaseq_library_compare": __version__, "python": platform.python_version()}
    for name in ("numpy", "pandas", "scipy", "matplotlib", "PyYAML", "pydeseq2", "anndata"):
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name] = "not installed"
    return result


def report(output, warnings, metrics, metadata, plans):
    warning_text = "\n".join(f"- {x}" for x in dict.fromkeys(warnings)) or "- No additional run-specific warnings."
    body = ["# RNA-seq library comparison", "", f"Samples: {len(metadata)}", f"Requested reference/query comparisons: {len(metrics)}",
            f"Replicate-aware DE analyses requested: {len(plans)}", "", "## Interpretation", "",
            *[f"- {x}" for x in LIMITATIONS], "", "## Run-specific warnings", "", warning_text, "",
            "## Main results", "", "- `pairwise_metrics.csv`: agreement, detection, and descriptive bias metrics.",
            "- `comparison_index.csv`: maps requested comparisons to per-gene CSVs and plots.",
            "- `reproducibility_summary.csv`: explicitly declared replicate groups.",
            "- `differential_expression_status.csv`: what significance testing did or did not run.",
            "- `sample_inventory.csv`, `sample_metadata.csv`, `gene_set_audit.csv`, `retained_assignments.csv`: input/selection audits.",
            "- `run_manifest.json`: input hashes, paths, package versions, and completion status.", ""]
    (output / "REPORT.md").write_text("\n".join(body), encoding="utf-8")
    blocks = ["<!doctype html><html><head><meta charset='utf-8'><title>RNA-seq library comparison</title>",
              "<style>body{font-family:system-ui,sans-serif;max-width:1100px;margin:2em auto;padding:0 1em;line-height:1.5}table{border-collapse:collapse;display:block;overflow-x:auto}td,th{border:1px solid;padding:.35em}img{max-width:100%}code{white-space:normal}</style></head><body>",
              "<h1>RNA-seq library comparison</h1>",
              f"<p>{len(metadata)} samples; {len(metrics)} descriptive comparisons; {len(plans)} replicate-aware DE analyses.</p>",
              "<h2>Interpretation</h2><ul>" + "".join(f"<li>{html.escape(x)}</li>" for x in LIMITATIONS) + "</ul>",
              "<h2>Run-specific warnings</h2><ul>" + "".join(f"<li>{html.escape(x)}</li>" for x in dict.fromkeys(warnings)) + "</ul>",
              "<h2>Pairwise results</h2>", metrics.to_html(index=False, escape=True, float_format=lambda x: f"{x:.4g}"),
              "<h2>Samples</h2>", metadata.to_html(escape=True)]
    for name in ("pca.png", "pearson_heatmap.png", "assignment_rarefaction.png", "metagene_dna.png", "metagene_rna.png", "metagene_combined.png"):
        if (output / name).is_file():
            blocks.append(f"<h2>{html.escape(name)}</h2><img src='{name}' alt='{html.escape(name)}'>")
    blocks.append("<h2>Result files</h2><ul>")
    for path in sorted(output.rglob("*.csv")):
        relative = path.relative_to(output).as_posix()
        blocks.append(f"<li><a href='{html.escape(relative, quote=True)}'>{html.escape(relative)}</a></li>")
    blocks.append("</ul></body></html>")
    (output / "report.html").write_text("\n".join(blocks), encoding="utf-8")


def run(args):
    config = load_config(args.config.resolve())
    warnings = []
    sources = load_sources(config, warnings)
    catalog = Catalog(sources, config)
    if args.list_samples:
        print(catalog.metadata(list(catalog.records)).reset_index().to_csv(index=False), end="")
        return 0
    pairs = comparison_pairs(config, catalog)
    replicate_groups = reproducibility_groups(config, catalog)
    plans = prepare_plans(config, catalog, warnings)
    requested = {pair[key] for pair in pairs for key in ("reference", "query")}
    requested.update(sample for group in replicate_groups for sample in group["samples"])
    requested.update(sample for plan in plans for sample in plan.samples)
    if not requested:
        raise InputError("No samples selected. Define comparison_groups, reproducibility_groups, or enabled differential_expression analyses.")
    selected = [sample for sample in catalog.records if sample in requested]
    metadata = catalog.metadata(selected)
    print(f"Loaded {len(sources)} output directories; selected {len(selected)} libraries.", flush=True)
    counts, annotation, genomic, gene_audit, retained = build_matrix(sources, selected, config, warnings)
    print(f"Analysis gene set: {len(counts)} CDS genes.", flush=True)
    for plan in plans:
        count_filter(counts, plan)
    # Warn when differing biological material makes method comparison confounded.
    for pair in pairs:
        a, b = metadata.loc[pair["reference"]], metadata.loc[pair["query"]]
        if not a.biological_sample or not b.biological_sample:
            warnings.append(f"{pair['group']}: biological identity is unspecified for at least one library; the comparison is agreement, not demonstrated technical reproducibility.")
        elif a.biological_sample != b.biological_sample:
            warnings.append(f"{pair['group']}: {pair['query']} and {pair['reference']} are different biological samples; biology and method effects are not separated by this pairwise comparison.")
        if a.counting_policy != b.counting_policy or a.count_layer != b.count_layer:
            warnings.append(f"{pair['group']}: {pair['query']} versus {pair['reference']} uses different native quantification policies/layers. Descriptive agreement is reported with this limitation; no harmonization has occurred.")
    fasta = config.get("reference", {}).get("genome_fasta", "auto")
    fasta_path = None
    if fasta == "auto":
        used_sources = [s for s in sources if s.id in set(metadata.source)]
        fasta_path, fasta_source = detect_genome_fasta(used_sources, config.get("organism", "bacteria"))
        if fasta_path:
            warnings.append(f"GC sequence reference automatically selected from {fasta_source}: {fasta_path}")
    elif fasta:
        fasta_path = config_path(fasta, Path(config["_base"]))
    if fasta_path and not fasta_path.is_file():
        raise InputError(f"genome_fasta not found: {fasta_path}")
    if plans and versions()["pydeseq2"] == "not installed":
        raise InputError("Differential-expression testing is enabled but PyDESeq2 is not installed. Activate the supplied environment.")
    if args.check:
        print(f"CHECK PASSED: {len(pairs)} descriptive comparisons, {len(replicate_groups)} reproducibility groups, {len(plans)} DE analyses.")
        print("This checks input/design validity, not statistical model convergence or PACE execution.")
        for warning in dict.fromkeys(warnings):
            print("WARNING:", warning)
        return 0
    if fasta_path:
        print("Computing genomic CDS GC fractions from the supplied FASTA...", flush=True)
        annotation["gc_percent"] = pd.Series(gc_from_fasta(fasta_path, genomic))
    else:
        warnings.append("No genome_fasta supplied: GC bias is unavailable (NA), not assumed to be zero. Length bias uses genomic CDS intervals from count tables, excluding RNA__ targets.")
    value = args.output_dir or config.get("output_dir", "results/library_comparison")
    output = config_path(str(value), Path(config["_base"]))
    for source in sources:
        if output == source.root or output in source.root.parents:
            raise InputError("Output directory cannot equal an upstream output directory or be its ancestor.")
    if output.exists() and any(output.iterdir()):
        raise InputError(f"Result directory is nonempty: {output}. Choose a new output_dir to preserve the previous analysis.")
    output.mkdir(parents=True, exist_ok=True)
    # Atomic claim prevents two simultaneously submitted jobs from publishing
    # over one another even when both initially observed an empty directory.
    with (output / ".analysis.lock").open("x") as lock:
        lock.write(f"pid={os.getpid()} slurm_job={os.environ.get('SLURM_JOB_ID', '')}\n")
    manifest = {"program": "rnaseq_library_compare", "status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
                "config": str(args.config.resolve()), "config_sha256": file_sha256(args.config.resolve()),
                "versions": versions(), "threads": args.threads, "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                "selected_samples": selected, "gene_count": len(counts), "inputs": [],
                "gc_genome_fasta": {"path": str(fasta_path), "sha256": file_sha256(fasta_path)} if fasta_path else None}
    for source in sources:
        manifest["inputs"].append({"id": source.id, "path": str(source.root), "pipeline": source.pipeline,
            "count_layer": source.layer, "native_counting_policy": source.counting_policy,
            "tables": [{"path": str(t.path), "stage": t.stage, "sha256": file_sha256(t.path),
                        "genes": len(t.counts), "comments": t.comments} for t in source.tables]})
    manifest_path = output / "run_manifest.json"
    def persist():
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    persist()
    try:
        snapshot = {key: value for key, value in config.items() if not key.startswith("_")}
        (output / "config.used.yaml").write_text(yaml.safe_dump(snapshot, sort_keys=False), encoding="utf-8")
        write_csv(catalog.metadata(list(catalog.records)), output / "sample_inventory.csv", index=True)
        write_csv(metadata, output / "sample_metadata.csv", index=True)
        write_csv(gene_audit, output / "gene_set_audit.csv")
        used_sources = [s for s in sources if s.id in set(metadata.source)]
        all_genes = set().union(*(set(s.counts.index) for s in used_sources))
        exclusions = []
        for gene in sorted(all_genes - set(counts.index)):
            present = [s.id for s in used_sources if gene in s.counts.index]
            missing = [s.id for s in used_sources if gene not in s.counts.index]
            exclusions.append({"gene_id": gene, "present_in_sources": ";".join(present),
                               "missing_from_sources": ";".join(missing),
                               "reason": "not_shared" if missing else "outside_selected_panel"})
        write_csv(pd.DataFrame(exclusions, columns=["gene_id", "present_in_sources", "missing_from_sources", "reason"]), output / "excluded_genes.csv")
        write_csv(retained, output / "retained_assignments.csv")
        write_csv(annotation, output / "gene_annotation.csv", index=True)
        write_csv(read_disposition(sources, selected, warnings), output / "upstream_read_disposition.csv")
        settings = config.get("analysis", {})
        cpm, log_cpm, metrics = pairwise_analysis(counts, annotation, pairs, settings, output)
        global_qc(counts, cpm, log_cpm, metadata, settings, output, warnings)
        depth_analysis(counts, output)
        reproducibility_analysis(counts, cpm, log_cpm, replicate_groups, settings, output)
        if settings.get("import_positional_profiles", True):
            positional_profiles(sources, selected, pairs, output, warnings)
        else:
            warnings.append("Positional-profile import disabled by configuration.")
        run_de(counts, annotation, plans, config, args.threads, output, warnings)
        report(output, warnings, metrics, metadata, plans)
        manifest.update(status="completed", completed_utc=datetime.now(timezone.utc).isoformat(), warnings=list(dict.fromkeys(warnings)))
        persist()
        print(f"COMPLETED: {output}\nOpen report.html; primary tables are pairwise_metrics.csv and differential_expression_status.csv.", flush=True)
        return 0
    except Exception as exc:
        manifest.update(status="failed", error=str(exc), warnings=list(dict.fromkeys(warnings)))
        persist()
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare expression, technical bias, and reproducibility across ScienceAndCats RNA-seq pipeline outputs.")
    parser.add_argument("config", type=Path, help="YAML configuration; relative paths resolve against this file's directory")
    parser.add_argument("--threads", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    parser.add_argument("--output-dir", type=Path, help="Override YAML output_dir; must be empty/new")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true", help="Validate counts, names, shared genes, and DE design without writing results")
    mode.add_argument("--list-samples", action="store_true", help="Print imported sample inventory as CSV; no analysis")
    parser.add_argument("--debug", action="store_true", help="Include Python traceback on failure")
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads must be positive")
    try:
        return run(args)
    except (InputError, OSError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        if args.debug:
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
