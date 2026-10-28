import numpy as np
import pandas as pd
import pytest

from rnaseq_compare.analysis import expected_detection, correlation, concordance, pair_statistics
from rnaseq_compare.de import prepare_plans
from rnaseq_compare.io import Catalog, InputError, load_sources


def make_de(config):
    config["differential_expression"] = {"enabled":True,"analyses":[{"name":"paired",
        "reference":["A1","A2","A3"],"query":["B1","B2","B3"],"paired_by":"biological_sample"}]}
    return config["differential_expression"]["analyses"][0]


def test_depth_expected_detection_exact():
    assert expected_detection([2,1],0) == 0
    assert expected_detection([2,1],1) == pytest.approx(1)
    assert expected_detection([2,1],2) == pytest.approx(5/3)
    assert expected_detection([2,1],3) == 2


def test_constant_correlation_not_fake_one():
    assert np.isnan(correlation([1,1,1],[1,2,3]))
    assert np.isnan(concordance([1,1],[1,1]))
    assert concordance([1,2,3],[1,2,3]) == 1


def test_counts_depth_invariance():
    frame=pd.DataFrame({"a":[10,20,40,80],"b":[100,200,400,800]})
    cpm=frame/frame.sum()*1e6
    result=pair_statistics(frame,cpm,np.log2(cpm+1),"a","b",10)
    assert result["median_absolute_log2_ratio"] == 0
    assert result["pearson_log2cpm"] == pytest.approx(1)


def test_paired_de_design(simple_config):
    make_de(simple_config)
    cat=Catalog(load_sources(simple_config,[]),simple_config)
    plan=prepare_plans(simple_config,cat,[])[0]
    assert plan.design.shape == (6,4)
    assert plan.design.query_vs_reference.tolist() == [0,0,0,1,1,1]


def test_no_unreplicated_de(simple_config):
    spec=make_de(simple_config)
    spec["reference"]=["A1"]
    spec["query"]=["B1"]
    cat=Catalog(load_sources(simple_config,[]),simple_config)
    with pytest.raises(InputError,match="TWO independent"):
        prepare_plans(simple_config,cat,[])


def test_missing_biological_ids_fail(simple_config):
    make_de(simple_config)
    simple_config["samples"]={}
    cat=Catalog(load_sources(simple_config,[]),simple_config)
    with pytest.raises(InputError,match="biological_sample metadata"):
        prepare_plans(simple_config,cat,[])


def test_technical_libraries_not_independent(simple_config):
    make_de(simple_config)
    simple_config["samples"]["A2"]["biological_sample"]="bio1"
    cat=Catalog(load_sources(simple_config,[]),simple_config)
    with pytest.raises(InputError,match="repeated libraries"):
        prepare_plans(simple_config,cat,[])


def test_paired_samples_not_unpaired(simple_config):
    spec=make_de(simple_config)
    del spec["paired_by"]
    cat=Catalog(load_sources(simple_config,[]),simple_config)
    with pytest.raises(InputError,match="both arms"):
        prepare_plans(simple_config,cat,[])


def test_confounded_batch_fails(simple_config):
    spec=make_de(simple_config)
    spec["covariates"]=["batch"]
    for name, meta in simple_config["samples"].items():
        meta["batch"]="one" if name.startswith("A") else "two"
    cat=Catalog(load_sources(simple_config,[]),simple_config)
    with pytest.raises(InputError,match="not full rank"):
        prepare_plans(simple_config,cat,[])


def test_mixed_native_longread_default_refused(simple_config):
    make_de(simple_config)
    simple_config["inputs"][1]["pipeline"]="longread"
    cat=Catalog(load_sources(simple_config,[]),simple_config)
    with pytest.raises(InputError,match="policies differ"):
        prepare_plans(simple_config,cat,[])
    simple_config["differential_expression"]["analyses"][0]["allow_mixed_native_counting"]=True
    warnings=[]
    assert prepare_plans(simple_config,cat,warnings)
    assert any("mixed native counting" in x for x in warnings)


def test_mixed_native_not_biological_even_with_override(simple_config):
    spec=make_de(simple_config)
    spec.update(allow_mixed_native_counting=True,purpose="biological")
    simple_config["inputs"][1]["pipeline"]="longread"
    cat=Catalog(load_sources(simple_config,[]),simple_config)
    with pytest.raises(InputError,match="policies differ"):
        prepare_plans(simple_config,cat,[])


def test_positional_coordinate_stages_not_mixed(tmp_path, fc_writer):
    from rnaseq_compare.analysis import positional_profiles
    source_paths=[tmp_path/"short",tmp_path/"long"]
    short=fc_writer(source_paths[0],["S"])
    long=fc_writer(source_paths[1],["L"],longread=True)
    config={"_base":str(tmp_path),"inputs":[{"id":"short","pipeline":"shortread","path":str(source_paths[0])},
        {"id":"long","pipeline":"longread","path":str(source_paths[1])}]}
    sources=load_sources(config,[])
    from pathlib import Path
    for source in sources:
        original=Path(next(iter(source.tables[0].columns.values())))
        path=original.with_suffix(".metagene_CDS_profile.csv")
        pd.DataFrame({"bin_start_percent":[0,50],"bin_end_percent":[50,100],"mean_fraction":[.5,.5],"features_contributing":[3,3]}).to_csv(path,index=False)
    output=tmp_path/"output"
    output.mkdir()
    positional_profiles(sources,["short::S","long::L"],[{"group":"test","reference":"short::S","query":"long::L"}],output,[])
    frame=pd.read_csv(output/"positional_pairwise.csv")
    assert frame.iloc[0].status == "unavailable"
    assert "No matching" in frame.iloc[0].reason
