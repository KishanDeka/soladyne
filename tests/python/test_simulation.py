import numpy as np
import pytest
from solar_dynamo import Parameters
from solar_dynamo.simulation import load_result,run,save_result
def test_frames_and_roundtrip(tmp_path):
    p=Parameters();result=run(p,None,frames=3,steps_per_frame=1);assert result.omega.shape==(3,257,257)
    assert np.diff(result.time)==pytest.approx([p.dt,p.dt]);restored=load_result(save_result(result,tmp_path/"run.npz"));np.testing.assert_array_equal(restored.toroidal,result.toroidal)
@pytest.mark.parametrize("frames,steps",[(1,1),(2,0)])
def test_controls(frames,steps):
    with pytest.raises(ValueError):run(Parameters(),None,frames,steps)
