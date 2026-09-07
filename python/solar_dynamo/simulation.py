"""Orchestration for the exact C++ v3 solver."""
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from . import Model, Parameters

@dataclass(frozen=True)
class RunResult:
    time: np.ndarray
    radius: np.ndarray
    theta: np.ndarray
    poloidal: np.ndarray
    toroidal: np.ndarray
    omega: np.ndarray
    energy: np.ndarray

def run(parameters: Parameters, init_file: str | Path | None = None,
        frames: int = 10, steps_per_frame: int = 5) -> RunResult:
    if frames < 2 or steps_per_frame < 1:
        raise ValueError("frames must be >= 2 and steps_per_frame >= 1")
    if init_file is None:
        parameters.relaxed_initial_state = False
        init_file = "init.dat"
    model=Model(parameters);model.initialize(Path(init_file));snapshots=[model.snapshot()]
    for _ in range(frames-1):model.step(steps_per_frame);snapshots.append(model.snapshot())
    poloidal=np.stack([np.asarray(s.poloidal) for s in snapshots]);toroidal=np.stack([np.asarray(s.toroidal) for s in snapshots])
    energy=np.mean(poloidal*poloidal+toroidal*toroidal,axis=(1,2))
    return RunResult(np.array([s.time for s in snapshots]),np.asarray(snapshots[0].radius),np.asarray(snapshots[0].theta),
      poloidal,toroidal,np.stack([np.asarray(s.omega) for s in snapshots]),energy)

def save_result(result:RunResult,path:str|Path)->Path:
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);np.savez_compressed(path,**result.__dict__);return path
def load_result(path:str|Path)->RunResult:
    with np.load(path) as data:return RunResult(**{k:data[k] for k in RunResult.__annotations__})
