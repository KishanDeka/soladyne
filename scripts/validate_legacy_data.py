"""Validate legacy SURYA text files without changing their format."""
import argparse
from pathlib import Path
import numpy as np

def load(path:Path,columns:int,rows:int|None=None):
    values=np.loadtxt(path)
    if values.ndim!=2 or values.shape[1]!=columns:raise ValueError(f"{path}: expected {columns} columns")
    if rows is not None and len(values)!=rows:raise ValueError(f"{path}: expected {rows} rows, got {len(values)}")
    if not np.isfinite(values).all():raise ValueError(f"{path}: contains non-finite values")
    return values
def main():
    p=argparse.ArgumentParser();p.add_argument("--init",type=Path);p.add_argument("--diffrot",type=Path);p.add_argument("--omg-lat",type=Path);a=p.parse_args()
    if a.init:print("init",load(a.init,4,257*257).shape)
    if a.diffrot:print("diffrot",load(a.diffrot,4,257*257).shape)
    if a.omg_lat:
        x=load(a.omg_lat,3);assert len(x)%257==0;times=x.reshape(-1,257,3)[:,0,1];print("omg_lat snapshots",len(times),"resets",(np.flatnonzero(np.diff(times)<0)+1).tolist())
if __name__=="__main__":main()
