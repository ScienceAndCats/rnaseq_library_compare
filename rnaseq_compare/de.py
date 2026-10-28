"""Replicate-aware PyDESeq2 analysis with explicit biological sample identities."""
from __future__ import annotations

import json
import re
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .analysis import save, slug, write_csv
from .io import InputError, config_path, read_gene_list
from pathlib import Path


@dataclass
class DEPlan:
    name: str
    purpose: str
    reference: list[str]
    query: list[str]
    metadata: pd.DataFrame
    design: pd.DataFrame
    settings: dict

    @property
    def samples(self):
        return self.reference + self.query


def prepare_plans(config, catalog, warnings):
    block = config.get("differential_expression") or {}
    if not block.get("enabled", False):
        if block.get("analyses"):
            warnings.append("Differential-expression analyses are configured but enabled is false; no significance tests will run.")
        return []
    analyses = block.get("analyses") or []
    if not analyses:
        raise InputError("differential_expression.enabled is true, but no analyses are defined.")
    plans = []
    for spec in analyses:
        if not isinstance(spec, dict) or not {"name", "reference", "query"} <= spec.keys():
            raise InputError("Each differential analysis needs name, reference: [...], and query: [...].")
        settings = {k: v for k, v in block.items() if k not in {"enabled", "analyses"}}
        settings.update(spec)
        name = str(spec["name"])
        if not isinstance(spec["reference"], list) or not isinstance(spec["query"], list):
            raise InputError(f"{name}: reference and query must be lists of sample names.")
        ref = [catalog.resolve(x) for x in spec["reference"]]
        query = [catalog.resolve(x) for x in spec["query"]]
        ids = ref + query
        if min(len(ref), len(query)) < 2:
            raise InputError(f"{name}: at least TWO independent biological samples per arm are required; a one-versus-one list is descriptive only. Three or more per arm is preferable.")
        if len(ids) != len(set(ids)):
            raise InputError(f"{name}: duplicated libraries or overlap between reference/query lists.")
        metadata = catalog.metadata(ids).copy()
        bios = metadata.biological_sample
        if bios.isna().any() or (bios.astype(str).str.strip() == "").any():
            raise InputError(f"{name}: biological_sample metadata is required for every library. FASTQ files are not automatically biological replicates.")
        for arm in (ref, query):
            if bios.loc[arm].nunique() != len(arm):
                raise InputError(f"{name}: an arm contains repeated libraries/lanes from one biological_sample. They cannot inflate biological n; explicitly combine compatible technical chunks upstream or select one library per biological sample.")
        purpose = settings.get("purpose", "technical")
        if purpose not in {"technical", "biological"}:
            raise InputError(f"{name}: purpose must be technical or biological.")
        if metadata.count_layer.nunique() > 1:
            raise InputError(f"{name}: DE cannot mix read counts and barcode-deduplicated molecule counts. Use the same count_layer for every input in this analysis.")
        if metadata.counting_policy.nunique() > 1:
            if purpose != "technical" or not settings.get("allow_mixed_native_counting", False):
                raise InputError(f"{name}: native long/short-read gene-counting policies differ (-O versus unique assignment). Keep this comparison descriptive, harmonize upstream counting, or explicitly set purpose: technical and allow_mixed_native_counting: true to test the COMBINED prep/platform/pipeline effect, not isolated biological expression.")
            warnings.append(f"{name}: mixed native counting explicitly accepted. The tested effect includes library prep, platform, alignment, length selection, and assignment rules; it does not isolate library-prep efficiency.")
        if metadata.pipeline.nunique() > 1:
            warnings.append(f"{name}: multiple upstream pipelines are included. Pipeline effects are not automatically removed by size-factor normalization.")
        if metadata.observation_unit.nunique() > 1:
            warnings.append(f"{name}: input observations mix single reads and read pairs. Results describe relative gene representation, not absolute molecule numbers.")
        if min(len(ref), len(query)) < 3:
            warnings.append(f"{name}: only two biological replicates in at least one arm; uncertainty and power can be poor.")
        metadata["comparison_arm"] = ["reference"] * len(ref) + ["query"] * len(query)
        paired_by = settings.get("paired_by")
        covariates = settings.get("covariates") or []
        if not isinstance(covariates, list) or len(covariates) != len(set(covariates)):
            raise InputError(f"{name}: covariates must be a list of distinct categorical metadata columns.")
        factors = ([paired_by] if paired_by else []) + covariates
        if len(factors) != len(set(factors)):
            raise InputError(f"{name}: do not repeat paired_by in covariates.")
        if not paired_by and set(bios.loc[ref]) & set(bios.loc[query]):
            raise InputError(f"{name}: biological samples occur in both arms. Specify paired_by: biological_sample rather than treating matched libraries as independent.")
        if paired_by:
            if paired_by not in metadata:
                raise InputError(f"{name}: missing pairing column {paired_by}")
            a, b = metadata.loc[ref, paired_by], metadata.loc[query, paired_by]
            if a.nunique() != len(ref) or b.nunique() != len(query) or set(a) != set(b):
                raise InputError(f"{name}: paired analysis requires the same complete set of blocks, exactly one library per block per arm.")
        design = pd.DataFrame({"Intercept": np.ones(len(ids))}, index=metadata.index)
        for factor in factors:
            if not isinstance(factor, str) or not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", factor):
                raise InputError(f"{name}: invalid metadata factor {factor!r}")
            if factor not in metadata or metadata[factor].isna().any() or (metadata[factor].astype(str).str.strip() == "").any():
                raise InputError(f"{name}: incomplete covariate/pairing metadata for {factor}")
            values = metadata[factor].astype(str)
            levels = sorted(values.unique())
            if len(levels) == 1:
                warnings.append(f"{name}: constant covariate {factor} omitted from the model.")
                continue
            for index, level in enumerate(levels[1:], 1):
                design[f"{factor}__level_{index}"] = (values == level).astype(float)
        design["query_vs_reference"] = (metadata.comparison_arm == "query").astype(float)
        if np.linalg.matrix_rank(design.values) != design.shape[1]:
            raise InputError(f"{name}: the design is not full rank. Group, pairing, and/or batch are confounded; a library-prep effect cannot be separated by this model.")
        if design.shape[0] <= design.shape[1]:
            raise InputError(f"{name}: no residual degrees of freedom. Add independent biological replicates or remove unsupported covariates.")
        if not 0 < float(settings.get("alpha", .05)) < 1:
            raise InputError(f"{name}: alpha must be in (0,1).")
        if settings.get("size_factors_fit_type", "poscounts") not in {"ratio", "poscounts", "iterative"}:
            raise InputError(f"{name}: unsupported size_factors_fit_type.")
        plans.append(DEPlan(name, purpose, ref, query, metadata, design, settings))
    if len({p.name for p in plans}) != len(plans):
        raise InputError("Differential-analysis names must be unique.")
    return plans


def count_filter(counts, plan):
    minimum = int(plan.settings.get("min_count", 10))
    samples = int(plan.settings.get("min_samples", 2))
    if minimum < 1 or not 1 <= samples <= len(plan.samples):
        raise InputError(f"{plan.name}: invalid min_count/min_samples for differential testing.")
    subset = counts[plan.samples]
    keep = (subset >= minimum).sum(axis=1) >= samples
    if keep.sum() < 2:
        raise InputError(f"{plan.name}: fewer than two genes pass the DE low-count filter.")
    subset = subset.loc[keep]
    if (subset.sum(axis=0) == 0).any():
        raise InputError(f"{plan.name}: a library has zero counts after the DE filter.")
    return subset, keep


def run_de(counts, annotation, plans, config, threads, output, warnings):
    if not plans:
        write_csv(pd.DataFrame([{"status": "not_requested", "reason": "No replicate-aware significance tests requested. Pairwise ratios are descriptive, not DE p-values."}]), output / "differential_expression_status.csv")
        return
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:
        raise InputError("PyDESeq2 is required for enabled differential analyses. Create/activate the supplied environment.yml; do not run them in the upstream pipeline environment.") from exc
    summaries = []
    destination = output / "differential_expression"
    destination.mkdir()
    for plan in plans:
        print(f"[DE] {plan.name}: {len(plan.reference)} reference and {len(plan.query)} query biological samples", flush=True)
        directory = destination / slug(plan.name)
        directory.mkdir()
        subset, keep = count_filter(counts, plan)
        controls = None
        if plan.settings.get("control_genes"):
            controls = read_gene_list(config_path(plan.settings["control_genes"], Path(config["_base"])))
            missing = set(controls) - set(subset.index)
            if missing:
                raise InputError(f"{plan.name}: normalization control genes absent or filtered: {sorted(missing)[:10]}")
        # An explicit numeric design avoids ambiguous category/reference coding.
        dds = DeseqDataSet(counts=subset.T, metadata=plan.metadata.copy(), design=plan.design.copy(),
                          size_factors_fit_type=plan.settings.get("size_factors_fit_type", "poscounts"),
                          control_genes=controls, refit_cooks=True, n_cpus=threads)
        dds.deseq2()
        contrast = np.zeros(plan.design.shape[1])
        contrast[-1] = 1  # QUERY minus REFERENCE, including paired/batch adjustment.
        alpha = float(plan.settings.get("alpha", .05))
        ds = DeseqStats(dds, contrast=contrast, alpha=alpha, cooks_filter=True,
                       independent_filter=True, n_cpus=threads)
        ds.summary()
        results = annotation.join(ds.results_df, how="left")
        results["passed_count_filter"] = keep
        results["significant_fdr"] = results.padj.le(alpha).fillna(False)
        results["absolute_log2fc_at_least_1"] = results.log2FoldChange.abs().ge(1).fillna(False)
        results["status"] = np.where(~keep, "excluded_low_counts",
                                    np.where(results.pvalue.isna(), "not_testable_or_cooks_filtered",
                                             np.where(results.padj.isna(), "independent_filter_excluded", "tested")))
        results["interpretation"] = "differential_representation" if plan.purpose == "technical" else "differential_expression"
        results.index.name = "gene_id"
        write_csv(results, directory / "results.csv", index=True)
        write_csv(results.loc[results.significant_fdr].sort_values("padj"), directory / "significant_genes.csv", index=True)
        write_csv(subset, directory / "counts_used.csv", index=True)
        write_csv(plan.metadata, directory / "sample_metadata.csv", index=True)
        write_csv(plan.design, directory / "design_matrix.csv", index=True)
        normalized = pd.DataFrame(np.asarray(dds.layers["normed_counts"]).T, index=subset.index, columns=subset.columns)
        write_csv(normalized, directory / "deseq2_normalized_counts.csv", index=True)
        factors = dds.obs["size_factors"] if "size_factors" in dds.obs else np.asarray(dds.obsm["size_factors"]).ravel()
        write_csv(pd.DataFrame({"sample_id": plan.samples, "size_factor": np.asarray(factors).ravel()}), directory / "size_factors.csv")
        write_csv(dds.var, directory / "fit_diagnostics.csv", index=True)
        (directory / "model.json").write_text(json.dumps({"name": plan.name, "purpose": plan.purpose,
            "direction": "query / reference", "reference": plan.reference, "query": plan.query,
            "design_columns": list(plan.design.columns), "contrast": contrast.tolist(),
            "settings": plan.settings, "lfc_shrinkage": "none; reported log2FoldChange is the unshrunk estimate",
            "pvalue_null": "log2 fold change = 0", "multiple_testing": "BH per analysis, with independent filtering"}, indent=2) + "\n")
        good = results.pvalue.notna() & np.isfinite(results.log2FoldChange)
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.scatter(results.loc[good, "log2FoldChange"], -np.log10(results.loc[good, "pvalue"].clip(lower=1e-300)), s=10, alpha=.4)
        ax.set(xlabel="Adjusted query/reference log2 fold change (unshrunk)", ylabel="−log10 raw p-value", title=plan.name)
        save(fig, directory / "volcano.png")
        summaries.append({"analysis": plan.name, "purpose": plan.purpose, "status": "completed",
            "reference_biological_n": len(plan.reference), "query_biological_n": len(plan.query),
            "genes_passing_count_filter": int(keep.sum()), "genes_with_adjusted_pvalue": int(results.padj.notna().sum()),
            "significant_fdr": int(results.significant_fdr.sum()), "alpha": alpha,
            "results_csv": str((directory / "results.csv").relative_to(output))})
    write_csv(pd.DataFrame(summaries), output / "differential_expression_status.csv")
    warnings.append("DE adjusted p-values control BH FDR separately for each analysis, not across every contrast in this run. Log2 fold changes are unshrunk; absolute-log2FC flags are descriptive effect-size filters, not threshold-null tests.")
