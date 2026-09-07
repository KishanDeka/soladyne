import numpy as np
import pytest
from solar_dynamo import Model,Parameters

def analytic_model():
    p=Parameters();p.relaxed_initial_state=False;m=Model(p);m.initialize("ignored.dat");return m
def test_exact_grid_and_initialization():
    m=analytic_model();s=m.snapshot();assert s.poloidal.shape==(257,257);assert s.toroidal.shape==(257,257)
    assert s.theta[0]==pytest.approx(np.pi) and s.theta[-1]==0;assert np.isfinite(s.omega).all()
def test_exact_time_and_boundaries():
    m=analytic_model();m.step(2);s=m.snapshot();assert m.time==pytest.approx(.0001);assert m.step_number==2
    assert np.allclose(s.toroidal[[0,-1]],0) and np.allclose(s.toroidal[:,[0,-1]],0)
def test_missing_relaxed_input_is_rejected(tmp_path):
    with pytest.raises(RuntimeError):Model(Parameters()).initialize(tmp_path/"missing.dat")
def test_final_legacy_formats(tmp_path):
    m=analytic_model();m.step(1);m.write_final_outputs(tmp_path)
    assert np.loadtxt(tmp_path/"diffrot.dat").shape==(66049,3)
    assert np.loadtxt(tmp_path/"final.dat").shape==(66049,5)
def test_legacy_output_cadence_and_preincrement_time(tmp_path):
    m=analytic_model();m.step_with_output(200,tmp_path)
    lat=np.loadtxt(tmp_path/"omg_lat.dat");diff=np.loadtxt(tmp_path/"diffrot.dat")
    assert lat.shape==(257,3) and diff.shape==(66049,4)
    assert lat[0,1]==pytest.approx(.00995) and diff[0,3]==pytest.approx(.00995)
