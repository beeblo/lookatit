# lookatit 🧪

**ORCA Trajectory & Geometry Visualizer**

`lookatit` is a lightweight, interactive web application built with Streamlit, `py3Dmol`, and `matplotlib`. It provides an intuitive GUI for post-processing ORCA relaxed surface scan trajectories (`.trj` / `.xyz`) and static single-frame coordinates, offering real-time synchronization between 3D molecular structures and dual-axis energy/geometry plots.

---

## Key Features

- **Flexible File Parsing:** Reads both multi-frame scan trajectories (`.trj` / `.xyz`) and static single-frame coordinate files.
- **Smart Energy Extraction:** Automatically parses negative total electronic energies from comment headers (converting Hartrees to relative $\text{kcal/mol}$). Safely defaults to `0.0 kcal/mol` if energy headers are missing.
- **Dynamic Geometry Detection:** Automatically measures structural parameters based on 0-indexed atom input:
  - **2 Indices (`0, 1`):** Measures and plots **Bond Distance** ($\text{Å}$).
  - **3 Indices (`0, 1, 2`):** Measures and plots **Bond Angle** ($^\circ$) at the vertex atom.
  - **4 Indices (`0, 1, 2, 3`):** Measures and plots **Dihedral Angle** ($^\circ$).
- **Interactive 3D Overlays:** Highlights selected atoms, renders 3D cylinder bonds for measured parameters, toggles atom index overlays, and displays dynamic measurement text directly on the WebGL canvas.
- **Single-Frame Compatibility:** Handles single-frame `.xyz` files gracefully without slider UI errors.
- **Structure Export:** Export any active scan step or structure directly as a standalone `.xyz` file.

---

## Installation Guide

### Prerequisites

- Python 3.8 or higher installed on your system.

### 1. Set Up Environment

Create a project folder and set up a clean Python virtual environment:

```bash
mkdir lookatit
cd lookatit
python3 -m venv venv
```

Activate the virtual environment.

**Linux / macOS**

```bash
source venv/bin/activate
```

**Windows**

```dos
venv\Scripts\activate
```

### 2. Install Required Dependencies

Install all required libraries via `pip`:

```bash
pip install streamlit py3Dmol matplotlib numpy
```

---

## Usage Instructions

### 1. Launching the App

Save your application code as `lookatit.py` in your project folder, then start Streamlit:

```bash
streamlit run lookatit.py
```

Streamlit will launch a local server and print access URLs in your terminal:

- **Local URL:** `http://localhost:8501`
- **Network URL:** `http://192.168.x.x:8501`

### 2. Step-by-Step Workflow

#### Upload Input File

Open the web browser interface and use the left sidebar to upload an ORCA surface scan file (`.trj`) or coordinate file (`.xyz`).

#### Define Atom Indices

In the sidebar under **Geometrical Measurement**, enter 0-indexed atom numbers separated by commas:

- **Distance between Atom 0 and Atom 1:** `0, 1`
- **Angle centered at Atom 1:** `0, 1, 2`
- **Dihedral angle across 4 atoms:** `0, 1, 2, 3`

#### Toggle View Options

Check or uncheck **Show Atom Indexes** to toggle numeric overlays on individual atoms in the 3D viewport.

#### Scrub Scan Steps

For multi-frame trajectories, drag the **Select Scan Step / Frame** slider to navigate reaction coordinates. The red marker on the relative energy plot and blue marker on the geometric parameter curve will update in sync with the 3D view.

#### Export Frame

Click **Download Frame X (`.xyz`)** under **Export Structure** in the sidebar to download the currently visible 3D geometry as an `.xyz` file.

---

## Troubleshooting & FAQ

### Issue: `StreamlitInvalidMinMaxError` on single-frame files

**Fix:** Ensure you are using the latest version of `lookatit.py`. The slider is automatically disabled when only 1 frame is detected.

### Issue: Energy curve is completely flat (`0.0 kcal/mol`)

**Cause:** Your `.xyz` or coordinate file lacks ORCA comment headers containing the negative total electronic energy in Hartrees (line 2 of each coordinate block).

**Behavior:** `lookatit` automatically assigns a dummy energy of `0.0 kcal/mol` so you can still:

- View 3D structures.
- Inspect index labels.
- Measure geometries.
- Export frames.

### Issue: Cannot connect from a remote workstation or HPC node

**Fix:** Bind the app to all network interfaces when launching on a remote machine:

```bash
streamlit run lookatit.py --server.address 0.0.0.0 --server.port 8501
```

Then access the interface via:

```text
http://YOUR_SERVER_IP:8501
```

### Issue: Atom labels clutter the 3D view

**Fix:** Uncheck **Show Atom Indexes** in the left sidebar to hide atom numbering labels on the WebGL canvas.
