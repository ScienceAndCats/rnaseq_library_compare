"""Strand-aware aligned-block classification and annotated RNA feature indexing.

The main-category rules follow bacteria_RNAseq_pipeline/alignment_categories.py:
rRNA > CDS > other RNA > intergenic/unannotated. Generic genomic gene and
transcript spans do not make introns RNA. Coordinates here are zero-based,
half-open; exported genomic intervals are one-based, inclusive.
"""
from __future__ import annotations

from collections import defaultdict
import gzip
from urllib.parse import unquote

from .io import InputError, merge_intervals


TRANSCRIPTS = {"mrna", "transcript", "primary_transcript"}
UTRS = {"utr", "five_prime_utr", "three_prime_utr", "5utr", "3utr", "5'utr", "3'utr"}
SUBTYPES = ("antisense_RNA", "tRNA_tmRNA", "sRNA_ncRNA", "other_noncoding_RNA", "UTR_transcript_regions")
PRIORITY = ("rRNA", "CDS", *SUBTYPES, "intergenic_unannotated")
BIN_SIZE = 16384


def subtype(typ, attrs):
    """Use explicit types/biotypes, never infer antisense from read orientation."""
    labels = {typ}
    labels.update(attrs.get(k, "").lower().replace("-", "_") for k in
                  ("ncRNA_class", "gene_biotype", "transcript_biotype", "biotype"))
    if any("antisense" in x for x in labels):
        return "antisense_RNA"
    if labels & {"trna", "tmrna"}:
        return "tRNA_tmRNA"
    if labels & {"srna", "ncrna", "nc_rna", "small_rna"}:
        return "sRNA_ncRNA"
    if typ in TRANSCRIPTS or typ in UTRS or typ == "exon":
        return "UTR_transcript_regions"
    return "other_noncoding_RNA"


def category(typ, attrs):
    if typ == "rrna":
        return "rRNA"
    if typ == "cds":
        return "CDS"
    if "rna" in typ or typ in TRANSCRIPTS | UTRS | {"ribozyme"}:
        return subtype(typ, attrs)
    return None


def read_gff(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise InputError(f"Expected nine GFF columns: {path}:{number}")
            attrs = {}
            for item in fields[8].strip(";").split(";"):
                if "=" in item:
                    key, value = item.strip().split("=", 1)
                    attrs[key] = unquote(value)
            try:
                start, end = int(fields[3]) - 1, int(fields[4])
            except ValueError as exc:
                raise InputError(f"Invalid GFF coordinates: {path}:{number}") from exc
            if start < 0 or end <= start or fields[6] not in {"+", "-", ".", "?"}:
                raise InputError(f"Invalid GFF coordinates/strand: {path}:{number}")
            yield fields[0], start, end, fields[6], fields[2].lower(), attrs


def merged(spans):
    result = []
    for start, end in sorted(spans):
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(end, result[-1][1]))
        else:
            result.append((start, end))
    return result


def complement(spans, length):
    start = 0
    for lo, hi in merged(spans):
        if start < lo:
            yield start, lo
        start = max(start, hi)
    if start < length:
        yield start, length


class RNAIndex:
    def __init__(self, gff, references, kind="dna", strand="forward", bin_bp=1000, gene_map=None):
        self.kind = kind
        self.bins = defaultdict(list)
        self.features = {}
        self.intervals = {}
        self.relative_intervals = {}
        self.references = dict(references)
        self.reference_signature = tuple(sorted(references.items()))
        rows = list(read_gff(gff))
        by_id = {r[5]["ID"]: r for r in rows if r[5].get("ID")}
        exon_parents = {parent for r in rows if r[4] == "exon" for parent in r[5].get("Parent", "").split(",")}
        coverage = defaultdict(list)
        spans_by_feature = defaultdict(list)
        canonical_origins = {}

        def ancestors(row):
            seen, pending = set(), row[5].get("Parent", "").split(",")
            while pending:
                ident = pending.pop()
                if ident in seen or ident not in by_id:
                    continue
                seen.add(ident)
                parent = by_id[ident]
                yield parent
                pending.extend(parent[5].get("Parent", "").split(","))

        def add(ref, start, end, selected_strand, cls, ident):
            fid = f"{cls}::{ident}"
            self.features[fid] = cls
            spans_by_feature[fid].append((ref, start, end, selected_strand))
            coverage[(ref, selected_strand)].append((start, end))
            for b in range(start // BIN_SIZE, (end - 1) // BIN_SIZE + 1):
                self.bins[(ref, selected_strand, b)].append((start, end, cls, fid))

        for row in rows:
            ref, start, end, selected_strand, typ, attrs = row
            if ref not in references or end > references[ref]:
                raise InputError(f"GFF reference/coordinates do not match the alignment header: {gff}, {ref}:{start+1}-{end}")
            chain = list(ancestors(row))
            cls = category(typ, attrs)
            if typ == "exon":
                classes = {category(p[4], p[5]) for p in chain}
                cls = next((c for c in PRIORITY if c in classes and c != "CDS"), "UTR_transcript_regions")
            if cls is None:
                continue
            if kind == "dna" and (attrs.get("ID") in exon_parents or typ in TRANSCRIPTS):
                continue
            # Prefer an explicit locus ID; exon children inherit their RNA
            # parent's identity so multiple exons do not become separate RNAs.
            owners = ([p for p in chain if category(p[4], p[5]) == cls] + [row]) if typ == "exon" else [row]
            ident = next((p[5][key] for p in owners for key in ("locus_tag", "gene_id", "locus", "gene", "ID", "Parent") if p[5].get(key)),
                         f"{kind}:{ref}:{start+1}-{end}:{selected_strand}")
            mapped = (gene_map or {}).get(ident, ident)
            if (cls, mapped) in canonical_origins and canonical_origins[(cls, mapped)] != ident:
                raise InputError(f"gene_id_map merges distinct RNA features into {cls}::{mapped}")
            canonical_origins[(cls, mapped)] = ident
            add(ref, start, end, selected_strand, cls, mapped)

        if not self.features:
            raise InputError(f"No eligible CDS, rRNA, exon/UTR or other RNA annotation features in {gff}")

        # Reference-anchored bins on annotation-free intervals provide stable,
        # explicit features for intergenic analysis, rather than fake gene IDs.
        strands = (".",) if strand == "unstranded" else ("+", "-")
        for ref, length in sorted(references.items()):
            for selected_strand in strands:
                eligible = ("+", "-", ".", "?") if strand == "unstranded" else (selected_strand, ".", "?")
                annotated = [s for st in eligible for s in coverage.get((ref, st), ())]
                for lo, hi in complement(annotated, length):
                    start = lo
                    while start < hi:
                        end = min(hi, (start // bin_bp + 1) * bin_bp)
                        ident = f"{kind}:{ref}:{start+1}-{end}:{selected_strand}"
                        add(ref, start, end, selected_strand, "intergenic_unannotated", ident)
                        start = end
        for fid, spans in spans_by_feature.items():
            self.relative_intervals[fid] = merge_intervals((ref, start+1, end, st) for ref, start, end, st in spans)
            self.intervals[fid] = self.relative_intervals[fid] if kind == "dna" else ()

    def relative_position(self, fid, ref, blocks):
        """Median overlapping aligned base on the spliced annotated feature.

        Uses sequenced bases only, including for pairs; never the unsequenced
        insert. Minus-strand positions run from the feature's biological 5' end.
        """
        spans = self.relative_intervals[fid]
        if len({(chrom, st) for chrom, _, _, st in spans}) != 1 or spans[0][3] not in {"+", "-"}:
            return None  # A multi-locus feature has no single 5'-to-3' profile.
        length = sum(end-start+1 for _, start, end, _ in spans)
        selected, offset = [], 0
        for chrom, start, end, _ in spans:
            for lo, hi in merged(blocks):
                left, right = max(start-1, lo), min(end, hi)
                if chrom == ref and right > left:
                    selected.append((offset + left-(start-1), offset + right-(start-1)))
            offset += end-start+1
        bases = sum(hi-lo for lo, hi in selected)
        if not bases:
            return None
        target = (bases-1) / 2
        for lo, hi in selected:
            if target < hi-lo:
                fraction = (lo + target + .5) / length
                return 1-fraction if spans[0][3] == "-" else fraction
            target -= hi-lo
        return None

    def classify(self, ref, blocks, strand):
        strands = ("+", "-", ".", "?") if strand is None else (strand, ".", "?")
        hits = defaultdict(set)
        for lo, hi in blocks:
            for st in strands:
                for b in range(lo // BIN_SIZE, (hi - 1) // BIN_SIZE + 1):
                    for start, end, cls, fid in self.bins.get((ref, st, b), ()):
                        if start < hi and end > lo:
                            hits[cls].add(fid)
        for cls in PRIORITY:
            if hits[cls]:
                return cls, hits[cls]
        raise InputError(f"Alignment has no usable aligned blocks on {ref}; cannot classify it.")
