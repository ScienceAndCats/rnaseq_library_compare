#!/usr/bin/env python3
"""Write SYNTHETIC retained RNA/DNA alignments for all three pipeline layouts.

python examples/make_category_demo.py /path/to/new/demo
python compare_rnaseq.py /path/to/new/demo/config.yaml --threads 2
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pysam
import yaml


def write_alignment(path, reference, length, observations, paired):
    path.parent.mkdir(parents=True, exist_ok=True)
    with pysam.AlignmentFile(str(path), "w", header={"HD": {"VN": "1.6", "SO": "unsorted"},
                              "SQ": [{"SN": reference, "LN": length}]}) as handle:
        for i, start in enumerate(observations):
            for flag, position in ((83, start), (163, start+15)) if paired else ((0, start),):
                read = pysam.AlignedSegment(handle.header)
                read.query_name, read.flag = f"observation_{i}", flag
                read.reference_id, read.reference_start = 0, position
                read.mapping_quality, read.cigarstring = 42, "20M"
                read.query_sequence = "A" * 20
                if paired:
                    read.next_reference_id, read.next_reference_start = 0, start
                handle.write(read)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    root = args.directory.resolve()
    if root.exists() and any(root.iterdir()):
        parser.error("Use a new/empty demo directory.")
    root.mkdir(parents=True, exist_ok=True)
    types = ("CDS", "rRNA", "sRNA", "tRNA", "antisense_RNA", "misc_RNA", "five_prime_UTR")
    features = [(f"{typ}_{i}", typ, 1+(j*3+i)*200, 100+(j*3+i)*200)
                for j, typ in enumerate(types) for i in range(3)]
    length = 5000
    reference = root / "references"
    reference.mkdir()
    (reference / "genome.fa").write_text(">chrom\n" + ("ACGT" * (length//4)) + "\n")
    for kind, contig in (("dna", "chrom"), ("rna", "RNA__demo")):
        (reference / f"{kind}.gff3").write_text("##gff-version 3\n" + "".join(
            f"{contig}\tdemo\t{typ}\t{start}\t{end}\t.\t+\t.\tlocus_tag={fid};gene_id={fid}\n"
            for fid, typ, start, end in features))
    refs = {"bacteria": {"dna": {"gff": str(reference/"dna.gff3"), "fasta": str(reference/"genome.fa")},
                         "rna": {"gff": str(reference/"rna.gff3")}}}
    config = {"version": 1, "output_dir": "results", "organism": "bacteria", "inputs": [], "samples": {},
              "comparison_groups": {}, "analysis": {"min_count": 5}, "reproducibility_groups": [],
              "rna_categories": {"require_alignments": True, "write_read_assignments": True}}
    rng = np.random.default_rng(461)
    for j, (sid, pipeline) in enumerate((("short", "shortread"), ("agilent", "agilent_pe"), ("long", "longread"))):
        source = root / sid
        source.mkdir()
        base = source / ("bbsplit_alignments/bacteria" if pipeline == "longread" else "bowtie_alignments/bacteria")
        samples = [f"{sid}_{i+1}" for i in range(3)]
        native_counts, columns, entries = {}, {}, []
        for kind, contig in (("dna", "chrom"), ("rna", "RNA__demo")):
            native_counts[kind], columns[kind] = [], []
            for i, sample in enumerate(samples):
                directory = base / kind if pipeline == "longread" or kind == "rna" else base
                alignment = directory / f"BACTERIA_{sample}_22bp.sam"
                values = rng.integers(12, 35, size=len(features))
                values[j::4] *= 2
                observations = [start+10 for (_, _, start, _), count in zip(features, values) for _ in range(int(count))]
                observations.extend([4610] * (15+i+j))  # Unannotated reference interval.
                write_alignment(alignment, contig, length, observations, pipeline == "agilent_pe")
                totals = {"CDS": int(values[:3].sum()), "rRNA": int(values[3:6].sum()),
                          "other_RNA": int(values[6:].sum()), "intergenic_unannotated": 15+i+j}
                category = alignment.with_suffix(".categories.json")
                category.write_text(json.dumps({"aligned": len(observations), "categories": totals}))
                if len(entries) <= i:
                    entries.append({"sample": sample, "stages": []})
                entries[i]["stages"].append({"organism": "bacteria", "kind": kind, "categories": str(category)})
                native_counts[kind].append(values[:3])
                columns[kind].append(str(alignment))
        for kind in (("combined",) if pipeline == "longread" else ("dna", "rna")):
            directory = source if kind != "rna" else base / "rna"
            table = directory / "featurecounts_BACTERIA_summary.txt"
            header_columns = [str(base / f"BACTERIA_{s}_22bp.sam") for s in samples] if kind == "combined" else columns[kind]
            counts = np.array(native_counts["dna"]) + np.array(native_counts["rna"]) if kind == "combined" else np.array(native_counts[kind])
            with table.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["Geneid", "Chr", "Start", "End", "Strand", "Length", *header_columns])
                for gene, (fid, typ, start, end) in enumerate(features[:3]):
                    fields = ["chrom;RNA__demo", f"{start};{start}", f"{end};{end}", "+;+", 200] if kind == "combined" else [
                        "chrom" if kind == "dna" else "RNA__demo", start, end, "+", 100]
                    writer.writerow([fid, *fields, *counts[:, gene]])
        if pipeline == "longread":
            ref_manifest = source / "bbsplit_alignments/references/manifest.json"
            ref_manifest.parent.mkdir(parents=True)
        else:
            ref_manifest = source / "reference_manifest.json"
            (source / "run_manifest.json").write_text(json.dumps({"read_layout": "paired" if pipeline == "agilent_pe" else "single", "samples": entries}))
        ref_manifest.write_text(json.dumps(refs))
        config["inputs"].append({"id": sid, "path": sid, "pipeline": pipeline, "library_prep": f"synthetic_method_{j}"})
        for i, sample in enumerate(samples):
            config["samples"][sample] = {"source": sid, "sample": sample, "biological_sample": f"specimen_{i}"}
        config["reproducibility_groups"].append({"name": sid, "kind": "biological", "samples": samples})
    for i in range(1, 4):
        config["comparison_groups"][f"specimen_{i}"] = [f"short_{i}", f"agilent_{i}", f"long_{i}"]
    (root / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    print(root / "config.yaml")


if __name__ == "__main__":
    main()
