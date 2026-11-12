# Upstream compatibility snapshot — 2026-10-27

## RNA category extension

Category reconstruction was added after inspecting the following complete
upstream revisions:

| Repository | Commit |
|---|---|
| bacteria_RNAseq_pipeline | `b488fe8fc1450baafe1996a7aebf3203fd14890d` |
| agilentPE_bacteria_RNAseq_pipeline | `e7163280b778027ddeb3859f7e8cd233f179d9be` |
| longread_bacteria_RNAseq_pipeline | `e607d2978f768684f0f8aa10e74913c979be8648` |

`alignment_categories.py` in the short-read pipeline defines aligned-block,
strand-aware, exclusive rRNA > CDS > other RNA > intergenic/unannotated totals.
Its JSONs contain totals and a CDS feature set, **not read-level category
membership**. Gene-level CDS/rRNA featureCounts matrices are not sufficient to
reconstruct all other RNA or unannotated alignments. The extension therefore
reads retained stage alignments and their normalized annotations. Generic
genomic transcript spans do not label introns as RNA; exon and UTR handling
follows the short-read category rules. Other-RNA types/biotypes are further
split as documented in README.md.

Short-read/Agilent current-run manifests identify stage category JSONs (and
deduplicated molecule BAMs for Agilent), even when a stage has no CDS table.
Recognizable relocated paths are supported. Legacy stage-count headers provide
alignment paths; custom `category_stages` can supply missing path metadata.
Long-read discovery follows the recorded combined alignment's `rna/` and
`dna/` sibling stages and the prepared reference manifest, without counting
both the combined SAM and its components. RNA-only coordinates never supply
genomic lengths or GC values.

All category **feature** matrices use unique assignment in the winning class,
with ambiguous-feature observations retained in category totals. Category
metadata distinguishes this from each pipeline's native gene-counting policy.
This common classification follows the short-read rules across all adapters;
it does not undo prior alignment/routing/filter choices or make reference
assemblies interchangeable. In particular, modern long-read annotation rules
can differ from short-read handling of UTRs and strandless features; retained
long-read TSV totals are not relabeled as reclassified feature counts.

Native-CDS import and its limitations below are retained. Real upstream runs
were not available for validation; synthetic retained SAM/GFF runs exercise all
three adapters, both stages, proper-pair reversal, all twelve panels, and
category-specific PyDESeq2 testing.

## Original native-CDS adapter snapshot

The adapters were written from direct inspection of the three named repositories,
not inferred from their project names. Upstream `main` branches can change.
These are **file blob SHAs**, not whole-repository commit IDs.

| Repository / inspected file | Blob SHA |
|---|---|
| bacteria_RNAseq_pipeline / merge_featurecounts.py | 11a4e16c302404d1066b753bd6dc7d6c3406f622 |
| bacteria_RNAseq_pipeline / gene_position_profile.py | 8d6cb7d5319cbeb005a7fdd05a2c50caf906e6e6 |
| bacteria_RNAseq_pipeline / read_breakdown_csv.py | bc3b5164fd9a1777fd62f40d45b497264896444d |
| agilentPE_bacteria_RNAseq_pipeline / map_bacteria_with_decoys.sh | 9674d6d4e8a3464ec9f716ca6f8811f3c0fda0d7 |
| agilentPE_bacteria_RNAseq_pipeline / paired_featurecounts.py | c15a276ba796cd652b04208b65acd48a9047639c |
| longread_bacteria_RNAseq_pipeline / count_features.py | 4cc4aa04f7fc0006e55614e4a68efbd79fc9fc24 |
| longread_bacteria_RNAseq_pipeline / map_bacteria_with_decoys.sh | 3fc3761332a26a77d9a12be22ae7966b63f92bb3 |
| longread_bacteria_RNAseq_pipeline / prepare_references.py | 459b93d72c57919309e4992ac8109111ee550b9d |
| longread_bacteria_RNAseq_pipeline / README.md | 017636b8a737f7ecf272e5f8ee729191b3fa0ac7 |

## Short-read and Agilent adapters

Bacterial genomic CDS: `featurecounts_BACTERIA_summary.txt` at the output root.
RNA-stage CDS, when present:
`bowtie_alignments/bacteria/rna/featurecounts_BACTERIA_summary.txt`.
Host equivalents live under `bowtie_alignments/host/`, with RNA under `rna/`.
The program adds disjoint RNA/DNA stages by exact gene ID, never by row position.

Modern finalizers publish merged matrices at the project root, keeping original
alignment columns that identify sample-worker output paths. `run_manifest.json`
is used for pairedness/sample checks. There is no glob of all historical worker
matrices, and therefore no automatic double-counting of merged and worker copies.

Agilent retains regular paired counts and, for `MBC_MODE=count`, a separate
`molecule_counts/` tree mirroring the stage layout. Its paired count wrapper uses
`-p --countReadPairs -B -C -s 2`. The selected count layer is explicit. Agilent's
metagene profiles remain ordinary read-pair profiles rather than molecule counts.

## Long-read adapter

`featurecounts_BACTERIA_summary.txt` is already the combined RNA/DNA table.
Serial host tables are under `bbsplit_alignments/host/`; finalized Slurm host
tables are under `host/`. Finalized Slurm run roots have `run.json` and
`aggregation.completed.json`, rather than the short-read run-manifest layout.

The current long-read program uses BBSplit's long-read mapping settings, keeps
one primary alignment per input record through sequential stages, excludes
post-trimming reads >6000 bp, and calls `count_features.py --allow-overlap` for
CDS. That wrapper uses forward-stranded featureCounts (`-s 1`) and `-O`.
It does **not** apply the same unique-gene assignment policy as the inspected
short-read counting commands.

RNA contigs receive `RNA__` prefixes. Combined featureCounts annotations can
contain a locus in both genomic and RNA coordinates; the comparison program
excludes `RNA__` segments when determining genomic CDS lengths/GC boundaries.
It does not add the combined table to stage-specific tables a second time.

The DNA-only reference FASTA is recorded as
`manifest[organism]["dna"]["fasta"]` under the prepared reference directory.
The short-read normalized reference instead has a `dna.gff3` / `dna.fa` cache
pair. Automatic GC lookup follows those known formats only.

## What native import cannot harmonize

The count matrix does not reveal every filtering parameter or reverse a prior
routing/assignment choice. The program records the declared pipeline type and
known native policy, and preserves available featureCounts comments. It does
not reconstruct the full executed command for a historical run when the final
merge removed it. Use the upstream frozen configs/logs when verifying
strandedness, sampling, annotation versions, capture targets, and count rules.

Native matrices are sufficient for quantified agreement and appropriately
replicated representation models, but not for proving that a difference was
caused solely by library preparation. Cross-pipeline statistical tests remain a
study-design question, not merely a file-format conversion.
