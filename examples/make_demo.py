#!/usr/bin/env python3
"""Generate small SYNTHETIC outputs in all three supported layouts.

Usage: python examples/make_demo.py /path/to/new/demo_directory
Then:  python compare_rnaseq.py /path/to/new/demo_directory/config.yaml
No FASTQs, alignments, or biological conclusions are involved.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import yaml


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    root = args.directory.resolve()
    if root.exists() and any(root.iterdir()):
        parser.error("Use an empty/new demo directory.")
    root.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(54)
    genes = [f"demo_gene_{i:03d}" for i in range(100)]
    sequences, spans, position = [], [], 1
    for i, gene in enumerate(genes):
        length = 150 + 5 * i
        gc = .25 + .5 * i / len(genes)
        seq = "".join(rng.choice(list("ACGT"), size=length, p=[(1-gc)/2,gc/2,gc/2,(1-gc)/2]))
        spans.append((position, position+length-1))
        sequences.extend([seq,"N"*25])
        position += length+25
    reference = root/"references"
    reference.mkdir()
    fasta = reference/"genome.fa"
    sequence = "".join(sequences)
    fasta.write_text(">demo_chromosome\n"+"\n".join(sequence[i:i+80] for i in range(0,len(sequence),80))+"\n")
    gff = reference/"genome.gff3"
    gff.write_text("##gff-version 3\n"+"".join(f"demo_chromosome\tdemo\tCDS\t{a}\t{b}\t.\t+\t0\tID={gene};locus_tag={gene}\n" for gene,(a,b) in zip(genes,spans)))
    base = rng.lognormal(4.5,1.0,size=(len(genes),1))*rng.uniform(.8,1.2,size=(1,3))
    config = {"version":1,"output_dir":"results","organism":"bacteria","inputs":[],"samples":{},
              "comparison_groups":{},"reproducibility_groups":[],"reference":{"genome_fasta":"auto"},
              "analysis":{"min_count":10},"differential_expression":{"enabled":False,"analyses":[]}}
    for j,(sid,pipeline) in enumerate((("short","shortread"),("agilent","agilent_pe"),("long","longread"))):
        source = root/sid
        source.mkdir()
        sample_names = [f"{sid}_{i+1}" for i in range(3)]
        effects = np.ones((len(genes),1))
        if j:
            effects[:,0] = np.exp(np.linspace(-.7,.7,len(genes)) * j)
        values = rng.negative_binomial(30,30/(30+base*effects))
        stage_values = [("combined",values)] if pipeline=="longread" else [("dna",np.floor(values*.7).astype(int)),("rna",values-np.floor(values*.7).astype(int))]
        for stage, data in stage_values:
            directory = source/"bowtie_alignments/bacteria/rna" if stage=="rna" else source
            directory.mkdir(parents=True,exist_ok=True)
            alignment_dir = source/("bbsplit_alignments" if pipeline=="longread" else "bowtie_alignments")/"bacteria"
            if stage=="rna":
                alignment_dir /= "rna"
            alignment_dir.mkdir(parents=True,exist_ok=True)
            headers = [str(alignment_dir/f"BACTERIA_{name}_22bp.sam") for name in sample_names]
            with (directory/"featurecounts_BACTERIA_summary.txt").open("w",newline="") as handle:
                handle.write("# SYNTHETIC DEMONSTRATION COUNTS, NOT REAL DATA\n")
                writer=csv.writer(handle,delimiter="\t",lineterminator="\n")
                writer.writerow(["Geneid","Chr","Start","End","Strand","Length",*headers])
                for gene,(start,end),row in zip(genes,spans,data):
                    length=end-start+1
                    if stage=="combined":
                        coords=[f"demo_chromosome;RNA__{gene}",f"{start};1",f"{end};{length}","+;+",length*2]
                    elif stage=="rna":
                        coords=[f"RNA__{gene}",1,length,"+",length]
                    else:
                        coords=["demo_chromosome",start,end,"+",length]
                    writer.writerow([gene,*coords,*row])
            for header in headers:
                bins=np.arange(10)
                fractions=np.exp((bins-4.5)*(.05*j))
                fractions/=fractions.sum()
                path=Path(header).with_suffix(".metagene_CDS_profile.csv")
                with path.open("w",newline="") as handle:
                    writer=csv.writer(handle)
                    writer.writerow(["sample","feature_class","bin","bin_start_percent","bin_end_percent","mean_fraction","median_fraction","features_contributing","raw_read_count"])
                    writer.writerows([Path(header).stem,"CDS",i+1,10*i,10*(i+1),f,f,len(genes),100] for i,f in enumerate(fractions))
        ref_manifest={"bacteria":{"dna":{"fasta":str(fasta),"gff":str(gff)},"rna":None}}
        if pipeline=="longread":
            target=source/"bbsplit_alignments/references"
            target.mkdir()
            (target/"manifest.json").write_text(json.dumps(ref_manifest))
        else:
            (source/"reference_manifest.json").write_text(json.dumps(ref_manifest))
            (source/"run_manifest.json").write_text(json.dumps({"read_layout":"paired" if pipeline=="agilent_pe" else "single", "samples":[{"sample":x} for x in sample_names]}))
        config["inputs"].append({"id":sid,"path":sid,"pipeline":pipeline,"library_prep":f"synthetic_method_{j+1}"})
        for i,sample in enumerate(sample_names):
            config["samples"][sample]={"source":sid,"sample":sample,"biological_sample":f"synthetic_specimen_{i+1}"}
        config["reproducibility_groups"].append({"name":sid+"_between_specimens","kind":"biological","samples":sample_names})
    for i in range(1,4):
        config["comparison_groups"][f"specimen_{i}"]=[f"short_{i}",f"agilent_{i}",f"long_{i}"]
    (root/"config.yaml").write_text(yaml.safe_dump(config,sort_keys=False))
    print(root/"config.yaml")


if __name__ == "__main__":
    main()
