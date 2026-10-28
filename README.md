# RNA-seq library-preparation comparison

Compare gene representation, technical bias, and reproducibility across multiple
outputs of these ScienceAndCats pipelines:

- `bacteria_RNAseq_pipeline` (`pipeline: shortread`)
- `agilentPE_bacteria_RNAseq_pipeline` (`pipeline: agilent_pe`)
- `longread_bacteria_RNAseq_pipeline` (`pipeline: longread`)

The program reads **native CDS featureCounts matrices**, not FASTQs. It does not
rerun alignment, change upstream results, or require the three repositories to
be installed alongside it. A PACE Slurm launcher is included.

**Two different questions remain separate:** single-library comparisons report
expression ratios and agreement; replicate-aware analyses use PyDESeq2 to test
systematic differences. A list with one reference library is not an experiment
with biological replication.

## Quick start on PACE

From this repository directory:

```bash
conda env create -f environment.yml
conda activate rnaseq-library-compare

# Edit paths, preparation labels, and comparison lists before continuing.
cp config.yaml config.local.yaml
nano config.local.yaml

# Optional: see the exact imported names, then validate the requested analysis.
python compare_rnaseq.py config.local.yaml --list-samples
python compare_rnaseq.py config.local.yaml --check

sbatch run_comparison.sbatch config.local.yaml
```

The job defaults to account `gts-mwhiteley3-biocluster`, partition `cpu-small`,
QOS `inferno`, four CPUs, 32 GB RAM, and 12 hours. These are editable resource
requests, not predicted resource consumption. Activate the **new** environment
before submission. No R installation is required.

From another directory, use the helper, which captures absolute paths:

```bash
bash /path/to/rnaseq_library_compare/submit.sh /path/to/config.local.yaml
```

The batch script intentionally uses the real repository path or
`SLURM_SUBMIT_DIR`, **not the directory of Slurm's temporary script copy**.
Direct `sbatch` submission should be made from the repository directory.

For a local/interactive run:

```bash
python compare_rnaseq.py config.local.yaml --threads 4
```

All YAML paths, including `output_dir`, are relative to the YAML's own directory
unless absolute. Results require an empty/new directory: reruns never overwrite
an earlier analysis. Use a different `output_dir` or `--output-dir` for a new run.

## Configuration

### Multiple input directories

```yaml
inputs:
  - id: neb_run_1
    path: /path/to/first_output
    pipeline: shortread
    library_prep: NEB
    count_layer: reads
  - id: agilent_run_1
    path: /path/to/second_output
    pipeline: agilent_pe
    library_prep: Agilent
    count_layer: reads
  - id: ont_run_1
    path: /path/to/longread_output/slurm_runs/EXACT_RUN_ID
    pipeline: longread
    library_prep: ONT_method
    count_layer: reads
```

Add any number of inputs, including several runs of the same pipeline. Input IDs
must be unique. Use the completed, published output directory. For long-read
Slurm outputs, this is normally `OUTPUT_DIR/slurm_runs/<run-id>/`. If a parent
contains exactly one usable run, it can be resolved automatically; with multiple
runs the program refuses to guess which one is intended.

`count_layer: molecules` selects Agilent's `molecule_counts/` tree from a run
made with `MBC_MODE=count`. The default is ordinary read-pair counts. The trees
are never added together. DE comparisons cannot mix read and deduplicated
molecule layers. There is no automatic UMI conversion for other pipelines.

### First sample is always the reference

```yaml
comparison_groups:
  first_set:
    - MW147-J-EM0330_S2_L001_R1_001
    - MW155_3_R1
    - MW155_5_R1
  second_set:
    - MW147-J-EM0330_S2_L001_R1_001
    - MW155_6_R1
    - MW155_7_R1
    - MW155_8_R1
```

This produces five reference/query comparisons, not every possible pair. The
same reference can occur in several lists. Literal YAML lists of lists also work:

```yaml
comparison_groups:
  - [reference_library, test_library_1, test_library_2]
  - [reference_library, test_library_3]
```

`log2_ratio_cpm = log2((query CPM + pseudocount)/(reference CPM + pseudocount))`.
Positive values mean greater relative representation in the query. The default
pseudocount is **one CPM**, not one raw read. It is not used by the DE count model.

### Names and metadata

A unique FASTQ stem normally works, with or without `.fastq[.gz]`/`.fq[.gz]`.
For paired outputs, an `_R1` or `_R1_001` FASTQ stem can resolve to the pipeline's
pair-level sample name. **Long-read `_R1` and `_R2` stems remain separate
samples.** When a legacy short-read output lacks its run manifest, specify
`read_layout: paired` to enable paired filename aliases.

When a name exists in more than one output, qualify it as `input_id::sample`,
or define a short alias and metadata:

```yaml
samples:
  reference_A:
    source: neb_run_1
    sample: original_fastq_stem
    biological_sample: specimen_A
    library_prep: NEB
    batch: batch_1
```

Use `reference_A` in lists. CSV matrices use unambiguous `input_id::sample`
identifiers; `sample_inventory.csv` records source and original names. Do not
assign different biological IDs to lanes or repeated preps of the same material.

## What is calculated?

### Quantification and agreement

Raw CDS gene counts; shared-gene/panel CPM and log2-CPM; per-gene query/reference
ratios; Pearson correlation on log2-CPM; Spearman correlation; Lin's concordance;
regression slope; RMSE; median absolute log2 ratio; and the fraction within a
2-fold descriptive ratio. Low-count genes are excluded from agreement summaries
using a configurable threshold (default: at least 10 counts in either library),
but remain in the exported gene tables.

Detection overlap, baseline gene retention, query-only and reference-only genes,
and Jensen–Shannon distance are also reported. Global heatmaps and PCA use the
selected samples. PCA is on centered top-variable log2-CPM, not DESeq2 VST, and is
not batch-corrected.

### Technical bias

The workflow relates the pairwise ratio to **genomic CDS length**, **GC content**,
and **average expression**, with per-gene scatterplots, correlation/slope
summaries, and quantile-bin summaries. These are associations, not causal proof
of a biochemical bias. Length comes from unioned genomic CDS intervals; the
long-read table's combined RNA/genome `Length` value is deliberately not used
when it would count the same locus in two coordinate systems.

GC uses `reference.genome_fasta: auto` by default: the first available DNA-only
FASTA in the selected outputs' known reference manifests. The selected path and
hash are recorded. An explicit same-assembly genomic FASTA overrides this;
`.gz` files work. Missing references make GC unavailable (`NA`), not zero.
No external genomes are downloaded. Use `null` to disable GC analysis.

Expected detected genes at matched **gene-assignment depth** and assignment
rarefaction curves help distinguish simple depth effects. These are not physical
read subsampling, UMI complexity, or an estimate of PCR duplicates. With native
long-read `-O` counts, one read can contribute more than one gene assignment.

When present, existing CDS metagene CSVs supply midpoint-position summaries,
first/last-20% fractions, curves, and same-stage reference/query differences.
They remain explicitly **upstream midpoint profiles**, not per-base coverage.
RNA, DNA, and combined stages are never silently interchanged or averaged.
Curves with different bin grids are not automatically rebinned. Contributing gene
panels can differ, so positional comparisons are descriptive and not a matched-
gene coverage assay. Barcode-deduplicated positional profiles are not invented
from Agilent's ordinary read-based curves.

Read-disposition summaries are imported separately when present (rRNA, CDS,
other RNA, host, leftovers, etc.). Their denominators remain the upstream
read/read-pair totals, not CPM gene-assignment totals.

### Reproducibility

Declare actual replicate sets, not just methods to compare:

```yaml
reproducibility_groups:
  - name: NEB_technical_replicates
    kind: technical
    samples: [NEB_A1, NEB_A2, NEB_A3]
```

All libraries in a group must have the same `library_prep`. Technical sets should
represent the same biological material. A `kind: biological` set requires
unique biological IDs. Outputs include every replicate pair's agreement and
per-gene CPM mean, sample SD, and CV. CV summaries require the gene to meet the
count threshold in at least two and at least half of the group's libraries.
These descriptive metrics do not create extra biological degrees of freedom.

## Replicate-aware differential testing

The default config leaves significance testing **disabled** because the sample
lists alone do not establish replication. `examples/paired_methods.yaml` is a
complete template for three distinct biological specimens, each split between
a reference and a test preparation:

```yaml
differential_expression:
  enabled: true
  alpha: 0.05
  min_count: 10
  min_samples: 2
  size_factors_fit_type: poscounts
  analyses:
    - name: test_method_vs_reference
      purpose: technical
      reference: [reference_A, reference_B, reference_C]
      query: [test_A, test_B, test_C]
      paired_by: biological_sample
      covariates: []
```

Provide `biological_sample: A/B/C` in each matching library's metadata. The model
is the equivalent of `~ biological_sample + comparison_arm`. For independent
biological condition groups, use `purpose: biological`, distinct biological IDs,
and omit `paired_by`. Optional `covariates: [batch]` adds categorical adjustment.
This initial version supports additive categorical designs, not continuous
covariates, interactions, nested random effects, or incomplete matched blocks.

The program checks for at least **two independent biological samples per arm**,
complete one-library-per-block pairing, repeated biological samples within an
arm, zero-count libraries, full-rank designs, and residual degrees of freedom.
Three or more biological replicates per arm are preferable. Repeated technical
libraries/lanes must not masquerade as biological replicates; compatible lane
chunks should be combined explicitly upstream. There is no automatic collapse
of independent library preps.

PyDESeq2 0.5.4 fits its negative-binomial DESeq2-style count model to **raw integer
counts**, never CPM/TPM. Results include fitted log2 fold changes, SEs, Wald
statistics, p-values, BH-adjusted p-values, fit diagnostics, size factors, and
normalized counts. Cook's filtering and independent filtering are enabled.
Filtered/untestable values remain `NA`. All genes, including genes excluded by
the low-count filter, remain in the complete results with a status column.

FDR is controlled **per analysis**, not over all contrasts jointly. Log2 fold
changes are **unshrunk**. The `absolute_log2fc_at_least_1` flag is an effect-size
filter, not a formal test against a one-log2-unit threshold. Volcano plots use
raw p-values on their vertical axis; significance calls use adjusted p-values.

For matched library-prep tests, significant genes indicate **differential
representation by method/pipeline**, not biological regulation: the biological
material may be identical. Capture enrichment, broad composition changes, or
many systematically affected genes can violate the normalization assumptions.
A declared `analysis.gene_panel` restricts all comparisons to a common target
panel; optional `control_genes` provides independently justified invariant genes
for DE size-factor estimation. Neither setting is a universal bias correction.

## Long-read compatibility: supported, with explicit limits

The long-read output format is usable: its root featureCounts table already
combines sequential RNA/DNA alignments. It is imported **once**, whereas the
short-read/Agilent DNA and RNA tables are added across their disjoint stages.

However, the inspected long-read pipeline assigns overlapping reads with `-O`,
while the short-read pipelines use unique gene assignment. It also uses a
different mapper/read-length regime and excludes post-trimming reads >6000 bp.
This software does not harmonize those choices. See `UPSTREAM_COMPATIBILITY.md`.

Therefore:

- Long-read descriptive comparisons and within-long-read replicated DE are
  supported without pretending the assays are identical.
- Native long-versus-short inferential comparisons are blocked by default.
- An intentional **technical benchmark of the combined preparation/platform/
  pipeline effect** can opt in per analysis with `purpose: technical` and
  `allow_mixed_native_counting: true`. This is not a repair or a harmonization
  switch, and it is not accepted as an isolated biological-expression test.

To isolate library prep from mapping/counting, harmonize upstream references,
read assignment, filtering, and study design first. Normalization alone cannot
remove gene-specific protocol/platform effects.

## Annotation and input safeguards

Sources should use the same organism, assembly, and annotation. Exact gene IDs
are intersected by default, with a minimum 80% overlap relative to **each**
selected source. Missing annotation is never converted into zero expression.
`gene_join: strict` requires identical gene sets. The 80% threshold is a sanity
check, not proof of assembly identity. Contradictory genomic CDS coordinates
cause an error by default. An explicit coordinate-check override drops
length/GC values for conflicting genes; it does not make their counts equivalent.

An optional `reference.gene_id_map` CSV has columns
`source,gene_id,canonical_gene_id`. Only one-to-one mappings are accepted; no
fuzzy gene-symbol matching, automatic ID-prefix stripping, or many-to-one
aggregation of potentially overlapping features is performed.

The importer prefers `.txt` over `.csv`, supports large featureCounts annotation
fields, rejects fractional/negative counts and duplicate genes/samples, and
checks manifest sample lists when present. It does not recursively import every
historical run or both merged and worker matrices. Incomplete RNA/CDS-table
checks use retained annotation manifests when available; old outputs without
those annotations cannot be exhaustively verified. Original paths embedded in
count headers are retained for provenance and optional profile discovery;
alignment files themselves do not need to exist for gene-count comparisons.

The current adapters are for the inspected pipeline layouts. An older output
with different gene IDs, missing stage tables, or changed counting rules may
need regeneration or explicit mapping; the program does not guess a repair.

## Outputs

Start with `report.html` (offline, with PNGs and CSV links) or `REPORT.md`.

| File/directory | Contents |
|---|---|
| `counts_raw.csv`, `counts_cpm.csv`, `counts_log2cpm.csv` | Selected sample/gene matrices |
| `comparison_index.csv`, `pairwise_metrics.csv`, `pairwise/` | Requested first-sample comparisons, gene ratios, plots |
| `sample_qc.csv`, `pca*.csv`, `*_heatmap.png` | Library summaries and global structure |
| `technical_bias_associations.csv`, `technical_bias_bins.csv` | Length, GC, and abundance associations |
| `depth_matched_detection.csv`, `assignment_rarefaction.csv` | Assignment-depth sensitivity |
| `reproducibility_summary.csv`, `reproducibility_pairwise.csv`, `reproducibility/` | Declared replicate agreement/CVs |
| `positional_bias_summary.csv`, `positional_pairwise.csv`, `upstream_metagene_profiles.csv` | Optional upstream midpoint diagnostics |
| `upstream_read_disposition.csv` | Imported read/pair-category summaries |
| `differential_expression_status.csv`, `differential_expression/` | What significance testing ran, designs, and results |
| `sample_inventory.csv`, `sample_metadata.csv`, `gene_annotation.csv` | Resolved identities and covariates |
| `gene_set_audit.csv`, `excluded_genes.csv`, `retained_assignments.csv` | Shared-gene and count-retention audit |
| `config.used.yaml`, `run_manifest.json` | Configuration snapshot, paths, hashes, software versions, status |

Unavailable optional outputs have an explicit no-records marker or warnings;
missing statistics are not fabricated. A failed run is marked `status: failed`;
a killed/interrupted run may remain `running`. Only `completed` confirms success.

## Tests and a reproducible demonstration

```bash
python -m pytest -q
python examples/make_demo.py /path/to/new/synthetic_demo
python compare_rnaseq.py /path/to/new/synthetic_demo/config.yaml
```

The demonstration generates small **synthetic** outputs in all three layouts.
It is a software check, not evidence about any real library method. Tests cover
import, sample names, stage/count-layer separation, annotations, design guards,
statistics, reporting, and the Slurm spool-path issue. The actual PyDESeq2 engine
test runs when that dependency is installed and is explicitly skipped otherwise.
See `TESTING.md` for checks executed on the delivered version. CI installs the
full environment's Python requirements before running the same suite.

## Methods references

- [PyDESeq2 workflow](https://pydeseq2.readthedocs.io/en/stable/auto_examples/plot_minimal_pydeseq2_pipeline.html)
- [PyDESeq2 dataset API](https://pydeseq2.readthedocs.io/en/stable/api/docstrings/pydeseq2.dds.DeseqDataSet.html)
- [PyDESeq2 statistical tests API](https://pydeseq2.readthedocs.io/en/stable/api/docstrings/pydeseq2.ds.DeseqStats.html)
- [DESeq2 vignette: counts, replication, normalization, and designs](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [Slurm sbatch documentation](https://slurm.schedmd.com/sbatch.html)
