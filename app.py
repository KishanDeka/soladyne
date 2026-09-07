from io import BytesIO
from pathlib import Path
import tempfile
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
from solar_dynamo import BACKEND,Parameters
from solar_dynamo.simulation import run
from solar_dynamo.visualization import FIELDS,make_figure

st.set_page_config(page_title="SURYA Dynamo v3",page_icon="☀️",layout="wide")
st.title("SURYA Dynamo v3 · exact C++ port")
st.caption("Original 257×257 ADI scheme · object-oriented C++17 · pybind11")
with st.sidebar:
    st.header("Legacy parameters")
    mode=st.radio("Initial condition",["Analytic (irelax=0)","Upload init.dat"])
    uploaded=st.file_uploader("Four-column init.dat",type=["dat"]) if mode.startswith("Upload") else None
    v0=st.number_input("v0",value=-20.0,step=1.0);et0=st.number_input("et0",value=3.0,step=.1)
    et1=st.number_input("et1",value=.04,format="%.5f");al0=st.number_input("al0",value=30.0,step=1.0)
    frames=st.slider("Frames",2,30,8);steps=st.slider("Steps per frame",1,50,5)
    field=st.selectbox("Displayed field",list(FIELDS));clicked=st.button("Run exact solver",type="primary",use_container_width=True)
    st.caption(f"Backend: {BACKEND}; dt=0.00005")

@st.cache_data(show_spinner=False)
def simulate(v0,et0,et1,al0,frames,steps,init_bytes):
    p=Parameters();p.v0=v0;p.et0=et0;p.et1=et1;p.al0=al0
    if init_bytes is None:
        p.relaxed_initial_state=False;return run(p,None,frames,steps)
    p.relaxed_initial_state=True
    with tempfile.NamedTemporaryFile(suffix=".dat") as f:
        f.write(init_bytes);f.flush();return run(p,Path(f.name),frames,steps)

if clicked or "result" not in st.session_state:
    if mode.startswith("Upload") and uploaded is None:
        st.info("Upload init.dat, or choose the analytic Fortran initialization.");st.stop()
    with st.spinner("Building legacy coefficients and advancing ADI steps…"):
        st.session_state.result=simulate(v0,et0,et1,al0,frames,steps,uploaded.getvalue() if uploaded else None)
result=st.session_state.result
left,right=st.columns([2,1])
with left:
    frame=st.slider("Frame",0,len(result.time)-1,len(result.time)-1)
    fig=make_figure(result,frame,field);st.pyplot(fig,use_container_width=True);plt.close(fig)
with right:
    st.metric("Model time",f"{result.time[frame]:.7f}");st.metric("Magnetic energy",f"{result.energy[frame]:.3e}")
    st.line_chart({"magnetic energy":result.energy[:frame+1]})
    payload=BytesIO();np.savez_compressed(payload,**result.__dict__)
    st.download_button("Download snapshots (.npz)",payload.getvalue(),"surya_v3_snapshots.npz")
st.warning("The uploaded reference omg_lat.dat contains three appended runs with different time-step spacing; it is retained as format evidence, not a single-run parity baseline.")
