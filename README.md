# lookatit 🧪
**ORCA Trajectory & Geometry Visualizer**

`lookatit` is a lightweight, interactive web application built with Streamlit, `py3Dmol`, and `matplotlib`. It provides an intuitive GUI for post-processing ORCA relaxed surface scan trajectories (`.trj` / `.xyz`) and static single-frame coordinates, offering real-time synchronization between 3D molecular structures and dual-axis energy/geometry plots.

---

## Key Features

* **Flexible File Parsing:** Reads both multi-frame scan trajectories (`.trj` / `.xyz`) and static single-frame coordinate files.
* **Smart Energy Extraction:** Automatically parses negative total electronic energies from comment headers (converting Hartrees to relative $\text{kcal/mol}$). Safely defaults to `0.0 kcal/mol` if energy headers are missing.
* **Dynamic Geometry Detection:** Automatically measures structural parameters based on 0-indexed atom input:
  * **2 Indices (`0, 1`):** Measures and plots **Bond Distance** ($\text{Å}$).
  * **3 Indices (`0, 1, 2`):** Measures and plots **Bond Angle** ($^\circ$) at the vertex atom.
  * **4 Indices (`0, 1, 2, 3`):** Measures and plots **Dihedral Angle** ($^\circ$).
* **Interactive 3D Overlays:** Highlights selected atoms, renders 3D cylinder bonds for measured parameters, toggles atom index overlays, and displays dynamic measurement text right on the WebGL canvas.
* **Single-Frame Compatibility:** Handles single-frame `.xyz` files gracefully without slider UI errors.

---

## Installation Guide

### Prerequisites
Python 3.8 or higher.

### 1. Set Up Environment
```bash
mkdir lookatit
cd lookatit
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
