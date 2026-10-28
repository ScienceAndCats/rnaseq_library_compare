from pathlib import Path
import csv
import json
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


@pytest.fixture
def fc_writer():
    def write(root, samples, counts=None, genes=None, stage="dna", organism="bacteria", layer="reads", longread=False):
        root = Path(root)
        base = root / "molecule_counts" if layer == "molecules" else root
        if organism == "bacteria":
            directory = base / "bowtie_alignments/bacteria/rna" if stage == "rna" else base
        else:
            directory = base / ("bbsplit_alignments/host" if longread else "bowtie_alignments/host")
            if stage == "rna":
                directory /= "rna"
        directory.mkdir(parents=True, exist_ok=True)
        genes = genes or ["geneA", "geneB", "geneC"]
        data = np.asarray(counts if counts is not None else [[10*(i+1)+j for j in range(len(samples))] for i in range(len(genes))])
        path = directory / f"featurecounts_{organism.upper()}_summary.txt"
        columns = []
        for sample in samples:
            align_dir = root / ("bbsplit_alignments" if longread else "bowtie_alignments") / organism
            if stage == "rna":
                align_dir /= "rna"
            align_dir.mkdir(parents=True, exist_ok=True)
            columns.append(str(align_dir / f"{organism.upper()}_{sample}_22bp{' .molecules.bam'.strip() if layer == 'molecules' else '.sam'}"))
        with path.open("w", newline="") as handle:
            handle.write("# synthetic featureCounts fixture, not real sequencing data\n")
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["Geneid", "Chr", "Start", "End", "Strand", "Length", *columns])
            for i, gene in enumerate(genes):
                start, end = i * 300 + 1, i * 300 + 100 + i
                chrom = f"transcript_{gene}" if stage == "rna" else "chrom"
                if longread:
                    coords = [f"chrom;RNA__tx_{gene}", f"{start};1", f"{end};{end-start+1}", "+;+", str(2*(end-start+1))]
                elif stage == "rna":
                    coords = [chrom, "1", str(end-start+1), "+", str(end-start+1)]
                else:
                    coords = [chrom, str(start), str(end), "+", str(end-start+1)]
                writer.writerow([gene, *coords, *map(int, data[i])])
        return path
    return write


@pytest.fixture
def simple_config(tmp_path, fc_writer):
    a, b = tmp_path / "a", tmp_path / "b"
    fc_writer(a, ["A1", "A2", "A3"])
    fc_writer(b, ["B1", "B2", "B3"])
    return {"version": 1, "_base": str(tmp_path), "output_dir": str(tmp_path / "results"),
            "inputs": [{"id": "a", "path": str(a), "pipeline": "shortread", "read_layout": "single", "library_prep": "reference"},
                       {"id": "b", "path": str(b), "pipeline": "agilent_pe", "library_prep": "test"}],
            "comparison_groups": {"one": ["A1", "B1"]},
            "samples": {**{f"A{i}": {"biological_sample": f"bio{i}"} for i in range(1,4)},
                        **{f"B{i}": {"biological_sample": f"bio{i}"} for i in range(1,4)}}}
