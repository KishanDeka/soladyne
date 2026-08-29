#  Real-Time 2D Solar Dynamo Simulator (C++20 & Plotly WebGL)

High-performance 2D kinematic solar dynamo simulator. The compute engine is written in **C++20** and compiled into a native Python extension using **Pybind11**. Results are rendered via a dynamic **Streamlit** dashboard using **Plotly WebGL (`go.Heatmapgl`)** for hardware-accelerated client-side rendering.

---

##  Key Features

* **Accelerated Physics Engine:** Full Alternating Direction Implicit (ADI) semi-implicit solver written in pure C++ with $O(N)$ tridiagonal matrix solutions.
* **Seamless Python Interop:** Zero-copy C++/Python array handoffs via `pybind11::stl`.
* **Hardware-Accelerated WebGL Visualization:** Plotly WebGL heatmaps running on client GPUs to maintain smooth frame rates during real-time parameter tuning.
* **Interactive Control:** Adjust physical diffusivity ($\eta_0$), differential rotation ($\Omega_0$), and alpha effect ($\alpha_0$) dynamically from the UI.

---

##  Architecture Overview

```text
┌─────────────────────────┐        Pybind11 Binding        ┌────────────────────────┐
│   C++ DynamoSolver      │ ─────────────────────────────> │ Streamlit Dashboard    │
│ (ADI Solvers & Grid)    │ <───────────────────────────── │ (Interactive Controls) │
└─────────────────────────┘      Parameters & Controls     └────────────────────────┘
│                                                         │
│ Vector Arrays                                           │ WebGL Canvas
▼                                                         ▼
Poloidal (A) & Toroidal (B)                                GPU Accelerated Plotly
```

---

##  Prerequisites

* **C++ Compiler:** `g++` (10+) or `clang` supporting **C++20**
* **CMake:** Version 3.15+
* **Python:** 3.8 or higher

---

##  Quick Start & Installation

### Clone the Repository

```bash
git clone https://github.com/your-username/solar-dynamo-sim.git
cd solar-dynamo-sim
```

### Set Up Python Environment & Dependencies

**Linux / macOS:**

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

**Windows:**

```powershell
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

### Compile the C++ Engine Module

Build the C++ shared object library using CMake:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

This places the compiled `dynamo_cpp` binary (`.so` on Linux/macOS or `.pyd` on Windows) directly into the repository root.

### Launch the Streamlit Dashboard

```bash
streamlit run app.py
```

---

##  Testing

Run the integration test suite to verify the native C++ engine binding and solver updates:

```bash
pytest tests/
```

---

##  License

Distributed under the MIT License. See `LICENSE` for details.
