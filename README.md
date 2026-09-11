# 🧪 lookatit: ORCA Trajectory & Geometry Visualizer

**lookatit** is a lightweight, interactive Streamlit application designed for Quantum Chemistry practitioners to parse, visualize, analyze, and manipulate ORCA trajectory files (`.trj`) and coordinate files (`.xyz`).

---

## Key Features

- **Interactive 3D Viewport:** Click individual atoms directly inside the WebGL viewport to queue them for measurement or editing.
- **Camera-Persistent Coordinate Nudger:** Adjust atomic positions along the X, Y, and Z axes using step increments without resetting your 3D view angle or zoom level.
- **Geometric Analysis:** Compute real-time distances (2 atoms), angles (3 atoms), and dihedral angles (4 atoms) across trajectory scan steps.
- **Dual-Axis Energy Profile:** Plot relative energies ($\text{kcal/mol}$) alongside active geometric parameters across all scan frames.
- **Dynamic Frame Export:** Download active or modified coordinate frames directly as standalone `.xyz` files.

---

## Installation & Quick Start

### 1. Clone the Repository

```bash
git clone https://github.com/beeblo/lookatit.git
cd lookatit
```

### 2. Install Dependencies

Ensure you have Python 3.8+ installed, then install the required libraries:

```bash
pip install streamlit py3Dmol numpy matplotlib
```

### 3. Run the Application

```bash
streamlit run lookatit.py
```

---

## How to Use

### 1. Upload File

Open the sidebar and upload an ORCA `.trj` or multi-frame/single-frame `.xyz` file.

### 2. Measure Geometry

- Click **2 atoms** for **Distance** ($\text{\AA}$).
- Click **3 atoms** for **Angle** ($^\circ$).
- Click **4 atoms** for **Dihedral Angle** ($^\circ$).
- Click **Clear** in the sidebar to reset selections.

### 3. Nudge Coordinates

- Select a target atom by clicking it in the 3D viewer or entering its index manually.
- Set your desired step size ($\text{\AA}$).
- Use the `-X`, `+X`, `-Y`, `+Y`, `-Z`, and `+Z` buttons to translate the atom.

### 4. Export

Click **Download Frame (.xyz)** in the sidebar to save the current frame.

---

## Requirements

- `python >= 3.8`
- `streamlit`
- `py3Dmol`
- `numpy`
- `matplotlib`
