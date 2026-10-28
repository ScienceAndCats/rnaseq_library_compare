import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
import yaml

from rnaseq_compare.cli import main
from rnaseq_compare.io import Catalog, load_sources
from rnaseq_compare.de import prepare_plans, run_de


def test_full_descriptive_workflow(tmp_path, fc_writer):
    rng=np.random.default_rng(5)
    genes=[f"gene{i:03d}" for i in range(80)]
    counts=rng.integers(5,1000,size=(80,4))
    fc_writer(tmp_path/"source",["A","B","C","D"],counts,genes)
    config={"version":1,"output_dir":"results","inputs":[{"id":"short","path":"source","pipeline":"shortread","read_layout":"single","library_prep":"one"}],
            "samples":{"A":{"biological_sample":"bio1"},"B":{"biological_sample":"bio1"}},
            "comparison_groups":{"one":["A","B","C"]},
            "reproducibility_groups":[{"name":"technical","kind":"technical","samples":["A","B"]}],
            "analysis":{"min_count":10,"import_positional_profiles":True}}
    path=tmp_path/"config.yaml"
    path.write_text(yaml.safe_dump(config))
    assert main([str(path),"--check"]) == 0
    assert not (tmp_path/"results").exists()
    assert main([str(path)]) == 0
    result=tmp_path/"results"
    assert json.loads((result/"run_manifest.json").read_text())["status"] == "completed"
    assert len(pd.read_csv(result/"pairwise_metrics.csv")) == 2
    assert len(pd.read_csv(result/"reproducibility_summary.csv")) == 1
    assert pd.read_csv(result/"differential_expression_status.csv").iloc[0]["status"] == "not_requested"
    assert (result/"report.html").is_file()
    assert (result/"pca.png").is_file()
    # A rerun must preserve the first result instead of silently replacing it.
    assert main([str(path)]) == 1


def test_sbatch_spool_location_safe(tmp_path):
    repo=Path(__file__).resolve().parents[1]
    config=tmp_path/"config.yaml"
    config.write_text("version: 1\n")
    copied=tmp_path/"slurm_spool"/"slurm_script"
    copied.parent.mkdir()
    copied.write_text((repo/"run_comparison.sbatch").read_text())
    # A stub Python records arguments, allowing wrapper behavior to be checked
    # without Slurm or the actual configured data. Its path contains a space.
    stub=tmp_path/"python stub"
    calls=tmp_path/"calls.txt"
    stub.write_text('#!/usr/bin/env bash\nprintf "%s\\n" "$*" >> "$CALLS"\n')
    stub.chmod(0o755)
    import os
    env=dict(os.environ,SLURM_SUBMIT_DIR=str(repo),RNA_COMPARE_PYTHON=str(stub),CALLS=str(calls),SLURM_TMPDIR=str(tmp_path),SLURM_JOB_ID="test",SLURM_CPUS_PER_TASK="4")
    completed=subprocess.run(["bash",str(copied),str(config),"--check"],env=env,capture_output=True,text=True)
    assert completed.returncode == 0, completed.stderr
    assert str(repo/"compare_rnaseq.py") in calls.read_text()
    assert "slurm_spool/compare_rnaseq.py" not in calls.read_text()


def test_actual_pydeseq2_engine_when_installed(tmp_path, simple_config):
    pytest.importorskip("pydeseq2", reason="Full DE engine requires the supplied environment; not replaced by a mock significance test.")
    simple_config["differential_expression"]={"enabled":True,"analyses":[{"name":"synthetic_paired",
        "reference":["A1","A2","A3"],"query":["B1","B2","B3"],"paired_by":"biological_sample"}]}
    catalog=Catalog(load_sources(simple_config,[]),simple_config)
    warnings=[]
    plans=prepare_plans(simple_config,catalog,warnings)
    rng=np.random.default_rng(123)
    genes=[f"g{i:03d}" for i in range(250)]
    means=rng.uniform(80,500,size=(250,1))*np.ones((1,6))
    means[:20,3:]*=8
    shape=25
    counts=pd.DataFrame(rng.negative_binomial(shape,shape/(shape+means)),index=pd.Index(genes,name="gene_id"),columns=plans[0].samples)
    annotation=pd.DataFrame({"length_bp":np.full(250,400),"gc_percent":np.full(250,50)},index=counts.index)
    run_de(counts,annotation,plans,simple_config,1,tmp_path,warnings)
    path=next((tmp_path/"differential_expression").glob("*/results.csv"))
    results=pd.read_csv(path,index_col=0)
    assert results.iloc[:20].log2FoldChange.median() > 1.5
    assert results.pvalue.notna().any()
    assert (results.padj.dropna().between(0,1)).all()
