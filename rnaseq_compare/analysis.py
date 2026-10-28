"""Descriptive expression agreement, reproducibility, and technical bias.

These outputs are explicitly not differential-expression significance tests.
"""
from __future__ import annotations

import hashlib
import itertools
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import jensenshannon

from .io import InputError, optional_sidecar_paths


def slug(text: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("._")[:70] or "item"
    return f"{safe}_{hashlib.sha256(str(text).encode()).hexdigest()[:8]}"


def write_csv(frame: pd.DataFrame, path: Path, index: bool = False):
    path.parent.mkdir(parents=True, exist_ok=True)
    if len(frame.columns) == 0:
        frame = pd.DataFrame([{"status": "no_records", "detail": "No records for this optional output; see REPORT.md and configuration."}])
        index = False
    frame.to_csv(path, index=index, na_rep="NA", float_format="%.10g")


def save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def correlation(x, y, kind="pearson") -> float:
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    if kind == "spearman":
        return float(stats.spearmanr(x, y).statistic)
    return float(np.corrcoef(x, y)[0, 1])


def concordance(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 2:
        return np.nan
    covariance = np.mean((x - x.mean()) * (y - y.mean()))
    denominator = x.var() + y.var() + (x.mean() - y.mean()) ** 2
    return float(2 * covariance / denominator) if denominator else np.nan


def pair_statistics(counts, cpm, log_cpm, reference, query, min_count):
    a, b = counts[reference], counts[query]
    keep = (a >= min_count) | (b >= min_count)
    x, y = log_cpm.loc[keep, reference].values, log_cpm.loc[keep, query].values
    difference = y - x
    da, db = a >= min_count, b >= min_count
    both, union = int((da & db).sum()), int((da | db).sum())
    selected_a, selected_b = cpm[reference].values, cpm[query].values
    slope = intercept = np.nan
    if len(x) >= 3 and np.ptp(x) > 0:
        slope, intercept = np.polyfit(x, y, 1)
    return {
        "reference": reference, "query": query, "genes_total": len(counts),
        "genes_in_agreement_metrics": int(keep.sum()),
        "pearson_log2cpm": correlation(x, y),
        "spearman_cpm": correlation(cpm.loc[keep, reference], cpm.loc[keep, query], "spearman"),
        "lins_concordance_log2cpm": concordance(x, y),
        "regression_slope_log2cpm": slope, "regression_intercept_log2cpm": intercept,
        "rmse_log2cpm": float(np.sqrt(np.mean(difference ** 2))) if len(x) else np.nan,
        "median_log2_ratio": float(np.median(difference)) if len(x) else np.nan,
        "median_absolute_log2_ratio": float(np.median(np.abs(difference))) if len(x) else np.nan,
        "fraction_within_twofold": float(np.mean(np.abs(difference) <= 1)) if len(x) else np.nan,
        "reference_detected_genes": int(da.sum()), "query_detected_genes": int(db.sum()),
        "detected_in_both": both, "reference_only": int((da & ~db).sum()),
        "query_only": int((db & ~da).sum()), "detection_jaccard": both / union if union else np.nan,
        "reference_detection_retained_fraction": both / int(da.sum()) if da.sum() else np.nan,
        "jensen_shannon_distance_base2": float(jensenshannon(selected_a, selected_b, base=2)),
    }


def bias_associations(frame, group, reference, query):
    rows, bins = [], []
    selected = frame[frame["included_in_agreement_metrics"]]
    covariates = {"log2_length_bp": np.log2(selected.length_bp.where(selected.length_bp > 0)),
                  "gc_percent": selected.gc_percent,
                  "mean_log2cpm": selected.mean_log2cpm}
    for name, series in covariates.items():
        good = np.isfinite(series) & np.isfinite(selected.log2_ratio_cpm)
        x, y = series[good].values, selected.loc[good, "log2_ratio_cpm"].values
        slope = np.polyfit(x, y, 1)[0] if len(x) >= 3 and np.ptp(x) > 0 else np.nan
        rows.append({"group": group, "reference": reference, "query": query, "covariate": name,
                     "genes": len(x), "spearman_with_log2_ratio": correlation(x, y, "spearman"),
                     "linear_slope_log2_ratio_per_unit": slope})
        if len(x) >= 5 and len(np.unique(x)) >= 2:
            groups = pd.qcut(pd.Series(x), min(10, len(np.unique(x))), labels=False, duplicates="drop")
            for number in sorted(groups.dropna().unique()):
                chosen = np.asarray(groups == number)
                bins.append({"group": group, "reference": reference, "query": query, "covariate": name,
                             "quantile_bin": int(number) + 1, "genes": int(chosen.sum()),
                             "covariate_median": float(np.median(x[chosen])),
                             "log2_ratio_median": float(np.median(y[chosen])),
                             "log2_ratio_q25": float(np.quantile(y[chosen], .25)),
                             "log2_ratio_q75": float(np.quantile(y[chosen], .75))})
    return rows, bins


def make_pair_plots(frame, reference, query, destination, prefix):
    selected = frame[frame.included_in_agreement_metrics]
    title = f"{query}\nversus {reference}"
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(selected.reference_log2cpm, selected.query_log2cpm, s=9, alpha=.35, rasterized=True)
    if len(selected):
        lo = min(selected.reference_log2cpm.min(), selected.query_log2cpm.min())
        hi = max(selected.reference_log2cpm.max(), selected.query_log2cpm.max())
        ax.plot([lo, hi], [lo, hi], linestyle="--")
    ax.set(xlabel="Reference log2(CPM + pseudocount)", ylabel="Query log2(CPM + pseudocount)", title=title)
    save(fig, destination / f"{prefix}.scatter.png")
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(selected.mean_log2cpm, selected.log2_ratio_cpm, s=9, alpha=.35, rasterized=True)
    ax.axhline(0, linestyle="--")
    ax.set(xlabel="Mean log2(CPM + pseudocount)", ylabel="Query / reference log2 ratio", title="Descriptive MA plot\n" + title)
    save(fig, destination / f"{prefix}.ma.png")
    for column, xlabel in (("length_bp", "Genomic CDS length (bp; log scale)"), ("gc_percent", "Genomic CDS GC (%)")):
        good = np.isfinite(selected[column])
        if column == "length_bp":
            good &= selected[column] > 0
        if good.sum() < 3:
            continue
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.scatter(selected.loc[good, column], selected.loc[good, "log2_ratio_cpm"], s=9, alpha=.35, rasterized=True)
        if column == "length_bp":
            ax.set_xscale("log")
        ax.axhline(0, linestyle="--")
        ax.set(xlabel=xlabel, ylabel="Query / reference log2 ratio", title="Technical bias diagnostic\n" + title)
        save(fig, destination / f"{prefix}.{column}_bias.png")


def pairwise_analysis(counts, annotation, pairs, settings, output):
    pseudocount = float(settings.get("cpm_pseudocount", 1.0))
    minimum = int(settings.get("min_count", 10))
    if pseudocount <= 0 or minimum < 1:
        raise InputError("cpm_pseudocount must be positive and min_count at least 1.")
    cpm = counts.div(counts.sum(axis=0), axis=1) * 1_000_000
    log_cpm = np.log2(cpm + pseudocount)
    write_csv(counts, output / "counts_raw.csv", index=True)
    write_csv(cpm, output / "counts_cpm.csv", index=True)
    write_csv(log_cpm, output / "counts_log2cpm.csv", index=True)
    destination = output / "pairwise"
    destination.mkdir()
    metrics, bias, binned, index = [], [], [], []
    for pair in pairs:
        reference, query, group = pair["reference"], pair["query"], pair["group"]
        prefix = slug(group + "__" + query + "__vs__" + reference)
        frame = annotation.copy()
        frame["reference_count"] = counts[reference]
        frame["query_count"] = counts[query]
        frame["reference_cpm"], frame["query_cpm"] = cpm[reference], cpm[query]
        frame["reference_log2cpm"], frame["query_log2cpm"] = log_cpm[reference], log_cpm[query]
        frame["log2_ratio_cpm"] = log_cpm[query] - log_cpm[reference]
        frame["mean_log2cpm"] = (log_cpm[query] + log_cpm[reference]) / 2
        frame["reference_detected"] = counts[reference] >= minimum
        frame["query_detected"] = counts[query] >= minimum
        frame["included_in_agreement_metrics"] = frame.reference_detected | frame.query_detected
        write_csv(frame, destination / f"{prefix}.genes.csv", index=True)
        metrics.append(dict(group=group, **pair_statistics(counts, cpm, log_cpm, reference, query, minimum)))
        associations, bins = bias_associations(frame, group, reference, query)
        bias.extend(associations)
        binned.extend(bins)
        make_pair_plots(frame, reference, query, destination, prefix)
        index.append(dict(pair, gene_csv=f"pairwise/{prefix}.genes.csv", plot_prefix=f"pairwise/{prefix}"))
    metric_frame = pd.DataFrame(metrics)
    write_csv(metric_frame, output / "pairwise_metrics.csv")
    write_csv(pd.DataFrame(bias), output / "technical_bias_associations.csv")
    write_csv(pd.DataFrame(binned), output / "technical_bias_bins.csv")
    write_csv(pd.DataFrame(index), output / "comparison_index.csv")
    return cpm, log_cpm, metric_frame


def global_qc(counts, cpm, log_cpm, metadata, settings, output, warnings):
    minimum = int(settings.get("min_count", 10))
    eligible = (counts >= minimum).any(axis=1)
    matrix = log_cpm.loc[eligible]
    qc = metadata.copy()
    qc["gene_assignments_in_analysis"] = counts.sum(axis=0)
    qc["genes_detected_at_min_count"] = (counts >= minimum).sum(axis=0)
    qc["genes_with_any_count"] = (counts > 0).sum(axis=0)
    qc["top_10_genes_fraction_of_assignments"] = [counts[s].nlargest(10).sum() / counts[s].sum() for s in counts]
    write_csv(qc, output / "sample_qc.csv", index=True)
    for kind in ("pearson", "spearman"):
        corr = matrix.corr(method=kind, min_periods=3)
        write_csv(corr, output / f"all_sample_{kind}_log2cpm.csv", index=True)
        n = len(corr)
        if n >= 2:
            fig, ax = plt.subplots(figsize=(min(18, max(6, n * .6)), min(18, max(5, n * .6))))
            im = ax.imshow(corr.values, vmin=-1, vmax=1)
            ax.set_xticks(range(n), corr.columns, rotation=90, fontsize=7)
            ax.set_yticks(range(n), corr.index, fontsize=7)
            ax.set_title(f"{kind.capitalize()} agreement: log2(CPM + pseudocount)")
            fig.colorbar(im, ax=ax)
            save(fig, output / f"{kind}_heatmap.png")
    if matrix.shape[1] < 2 or matrix.shape[0] < 2:
        warnings.append("PCA unavailable: at least two samples and two expressed genes are required.")
        return
    top_n = int(settings.get("pca_top_genes", 2000))
    if top_n < 2:
        raise InputError("pca_top_genes must be at least 2.")
    top = matrix.var(axis=1).nlargest(min(top_n, len(matrix))).index
    x = matrix.loc[top].T.values
    x = x - x.mean(axis=0, keepdims=True)
    u, singular, _ = np.linalg.svd(x, full_matrices=False)
    if np.sum(singular ** 2) == 0:
        warnings.append("PCA unavailable: all normalized sample profiles are identical.")
        return
    coordinates = u * singular
    proportion = singular ** 2 / np.sum(singular ** 2)
    pcs = pd.DataFrame(coordinates[:, :min(5, len(singular))], index=counts.columns,
                       columns=[f"PC{i+1}" for i in range(min(5, len(singular)))])
    pcs.index.name = "sample_id"
    pcs["plot_label"] = np.arange(1, len(pcs) + 1)
    write_csv(pcs.join(metadata), output / "pca_scores.csv", index=True)
    write_csv(pd.DataFrame({"component": [f"PC{i+1}" for i in range(len(proportion))], "variance_fraction": proportion}), output / "pca_variance.csv")
    write_csv(pd.DataFrame({"gene_id": top}), output / "pca_genes.csv")
    fig, ax = plt.subplots(figsize=(9, 7))
    for method in metadata.library_prep.unique():
        ids = metadata.index[metadata.library_prep == method]
        ax.scatter(pcs.loc[ids, "PC1"], pcs.loc[ids, "PC2"], label=method)
    for sample in pcs.index:
        ax.annotate(str(int(pcs.loc[sample, "plot_label"])), (pcs.loc[sample, "PC1"], pcs.loc[sample, "PC2"]), fontsize=8, xytext=(4, 5), textcoords="offset points")
    ax.set(xlabel=f"PC1 ({proportion[0]:.1%})", ylabel=f"PC2 ({proportion[1]:.1%})",
           title="PCA of top-variable log2-CPM genes\nSample numbers: pca_scores.csv; not batch-corrected")
    ax.legend(fontsize=8)
    save(fig, output / "pca.png")


def expected_detection(counts, depth):
    """Exact expected nonzero categories in a hypergeometric assignment subsample.

    This is sampling GENE ASSIGNMENTS, not unique reads/molecules. In particular,
    native long-read -O counts are not independent physical reads across genes.
    """
    values = np.asarray(counts, dtype=np.int64)
    total, depth = int(values.sum()), int(depth)
    if depth < 0 or depth > total:
        raise InputError("Invalid assignment rarefaction depth.")
    if depth == total:
        return float((values > 0).sum())
    log_p0 = stats.hypergeom.logpmf(0, total, values, depth)
    return float(np.sum(-np.expm1(np.minimum(log_p0, 0))))


def depth_analysis(counts, output):
    common_depth = int(counts.sum(axis=0).min())
    rows, curve = [], []
    for sample in counts:
        total = int(counts[sample].sum())
        rows.append({"sample_id": sample, "gene_assignments": total, "matched_assignment_depth": common_depth,
                     "expected_genes_with_at_least_one_assignment": expected_detection(counts[sample], common_depth)})
        depths = sorted(set([max(1, int(total * f)) for f in np.linspace(.05, 1, 15)] + [common_depth]))
        for depth in depths:
            curve.append({"sample_id": sample, "assignment_depth": depth, "expected_detected_genes": expected_detection(counts[sample], depth)})
    write_csv(pd.DataFrame(rows), output / "depth_matched_detection.csv")
    data = pd.DataFrame(curve)
    write_csv(data, output / "assignment_rarefaction.csv")
    fig, ax = plt.subplots(figsize=(9, 6))
    for sample, part in data.groupby("sample_id", sort=False):
        ax.plot(part.assignment_depth, part.expected_detected_genes, label=sample)
    ax.set(xlabel="Sampled gene-assignment counts (NOT sequencing reads)", ylabel="Expected genes with ≥1 assignment",
           title="Assignment-depth sensitivity; not molecular complexity")
    if len(counts.columns) <= 15:
        ax.legend(fontsize=6)
    save(fig, output / "assignment_rarefaction.png")


def reproducibility_analysis(counts, cpm, log_cpm, groups, settings, output):
    destination = output / "reproducibility"
    destination.mkdir()
    minimum = int(settings.get("min_count", 10))
    summaries, all_pairs = [], []
    for group in groups:
        ids = group["samples"]
        subset = cpm[ids]
        # Per-gene CV is descriptive; libraries are not counted as biological
        # replicates for inference merely because they occur in this group.
        enough = (counts[ids] >= minimum).sum(axis=1) >= max(2, int(np.ceil(len(ids)/2)))
        means = subset.mean(axis=1)
        sds = subset.std(axis=1, ddof=1)
        cv = sds / means.replace(0, np.nan)
        frame = pd.DataFrame({"mean_cpm": means, "sd_cpm": sds, "cv_cpm": cv,
                              "included_in_group_summary": enough})
        write_csv(frame, destination / f"{slug(group['name'])}.gene_cv.csv", index=True)
        pair_rows = [dict(group=group["name"], kind=group["kind"], **pair_statistics(counts, cpm, log_cpm, a, b, minimum))
                     for a, b in itertools.combinations(ids, 2)]
        all_pairs.extend(pair_rows)
        rhos = [p["spearman_cpm"] for p in pair_rows if np.isfinite(p["spearman_cpm"])]
        summaries.append({"group": group["name"], "kind": group["kind"], "libraries": len(ids),
                          "genes_for_cv_summary": int(enough.sum()), "median_gene_cv_cpm": cv[enough].median(),
                          "median_pairwise_spearman": np.median(rhos) if rhos else np.nan})
    write_csv(pd.DataFrame(summaries), output / "reproducibility_summary.csv")
    write_csv(pd.DataFrame(all_pairs), output / "reproducibility_pairwise.csv")


def positional_profiles(sources, selected, pairs, output, warnings):
    """Import upstream midpoint metagene curves, preserving stage and provenance.

    These are not per-base coverage. Different upstream feature panels, strand
    handling, and midpoint definitions prevent treating curve differences as an
    isolated library-prep effect. Do not pool RNA/DNA means or compare combined
    long-read curves with DNA-only short-read curves.
    """
    rows, metrics = [], []
    for source in sources:
        for table in source.tables:
            for sample in table.counts.columns:
                key = f"{source.id}::{sample}"
                if key not in selected:
                    continue
                if source.layer == "molecules":
                    warnings.append(f"{key}: upstream positional profiles are read-based, so no molecule-level profile is substituted.")
                    continue
                path = optional_sidecar_paths(source, table, sample, ".metagene_CDS_profile.csv")
                if path is None:
                    continue
                frame = pd.read_csv(path)
                required = {"bin_start_percent", "bin_end_percent", "mean_fraction", "features_contributing"}
                if not required <= set(frame.columns):
                    warnings.append(f"Unsupported positional CSV schema, skipped: {path}")
                    continue
                if frame.empty:
                    continue
                for column in required:
                    frame[column] = pd.to_numeric(frame[column], errors="raise")
                if frame[list(required)].isna().any().any() or (frame.mean_fraction < 0).any():
                    raise InputError(f"Invalid positional fractions in {path}")
                frame = frame.sort_values("bin_start_percent")
                if frame.bin_start_percent.duplicated().any():
                    raise InputError(f"Duplicate bins in {path}")
                weights = frame.mean_fraction.values
                total = weights.sum()
                if total <= 0:
                    continue
                if not np.isclose(total, 1, atol=.02):
                    raise InputError(f"Metagene mean fractions do not sum to approximately 1: {path}")
                weights = weights / total
                centers = (frame.bin_start_percent.values + frame.bin_end_percent.values) / 2
                fractions = {}
                # Fractional bin overlap avoids boundary rounding for arbitrary bins.
                starts, ends = frame.bin_start_percent.values, frame.bin_end_percent.values
                if np.any(ends <= starts) or starts.min() < 0 or ends.max() > 100:
                    raise InputError(f"Invalid positional bin bounds: {path}")
                for label, lo, hi in (("first_20pct", 0, 20), ("last_20pct", 80, 100)):
                    overlap = np.maximum(0, np.minimum(ends, hi) - np.maximum(starts, lo)) / (ends-starts)
                    fractions[label] = float(np.sum(weights * overlap))
                metrics.append({"sample_id": key, "stage": table.stage,
                    "midpoint_definition": "fragment_midpoint" if source.layout == "paired" else "alignment_midpoint",
                    "features_contributing": int(frame.features_contributing.max()),
                    "mean_position_percent": float(np.sum(weights * centers)), **fractions,
                    "last20_over_first20": fractions["last_20pct"] / fractions["first_20pct"] if fractions["first_20pct"] else np.nan,
                    "source_file": str(path)})
                for row, fraction in zip(frame.to_dict("records"), weights):
                    rows.append({"sample_id": key, "stage": table.stage, "pipeline": source.pipeline,
                                 "source_file": str(path), **row, "mean_fraction": fraction})
    write_csv(pd.DataFrame(metrics), output / "positional_bias_summary.csv")
    data = pd.DataFrame(rows)
    write_csv(data, output / "upstream_metagene_profiles.csv")
    comparisons = []
    for pair in pairs:
        reference, query = pair["reference"], pair["query"]
        ref_stages = set(data.loc[data.sample_id == reference, "stage"]) if not data.empty else set()
        query_stages = set(data.loc[data.sample_id == query, "stage"]) if not data.empty else set()
        common = ref_stages & query_stages
        if not common:
            comparisons.append(dict(pair, status="unavailable", reason="No matching available coordinate stage; RNA/DNA/combined curves are not interchanged."))
        for stage in sorted(common):
            a = data[(data.sample_id == reference) & (data.stage == stage)].sort_values("bin_start_percent")
            b = data[(data.sample_id == query) & (data.stage == stage)].sort_values("bin_start_percent")
            grid = ["bin_start_percent", "bin_end_percent"]
            if a.shape[0] != b.shape[0] or not np.allclose(a[grid].values, b[grid].values):
                comparisons.append(dict(pair, stage=stage, status="unavailable", reason="Different bin grids; no automatic rebinning."))
                continue
            av, bv = a.mean_fraction.values, b.mean_fraction.values
            mid = (a.bin_start_percent.values + a.bin_end_percent.values) / 2
            comparisons.append(dict(pair, stage=stage, status="descriptive_only", same_gene_panel_verified=False,
                midpoint_js_distance_base2=float(jensenshannon(av, bv, base=2)),
                midpoint_l1_distance=float(np.abs(av-bv).sum()),
                query_minus_reference_mean_position=float(np.sum((bv-av)*mid)),
                reason="Imported upstream curves; contributing gene panels and midpoint definitions may differ."))
    write_csv(pd.DataFrame(comparisons), output / "positional_pairwise.csv")
    if data.empty:
        warnings.append("No available CDS metagene sidecars: positional bias was not evaluated. Gene-length/GC and abundance-bias diagnostics are separate.")
        return
    warnings.append("Positional summaries reuse upstream MIDPOINT profiles, not per-base coverage. RNA, DNA, and combined stages stay separate; contributing gene panels may differ. Agilent profiles are not UMI-deduplicated.")
    for stage, part in data.groupby("stage", sort=False):
        fig, ax = plt.subplots(figsize=(9, 6))
        for sample, curve in part.groupby("sample_id", sort=False):
            ax.plot((curve.bin_start_percent + curve.bin_end_percent)/2, curve.mean_fraction, label=sample)
        ax.set(xlabel="Normalized feature position (5′ → 3′)", ylabel="Mean per-feature midpoint fraction",
               title=f"Upstream CDS midpoint profiles: {stage} stage (descriptive)")
        if part.sample_id.nunique() <= 15:
            ax.legend(fontsize=6)
        save(fig, output / f"metagene_{stage}.png")
