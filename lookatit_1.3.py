import math
import matplotlib.pyplot as plt
import numpy as np
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components

# Set page layout and browser tab title
st.set_page_config(layout="wide", page_title="lookatit - ORCA Visualizer")

# GUI Header
st.title("🧪 lookatit: ORCA Trajectory & Geometry Visualizer")


# =============================================================================
# 1. VECTOR MATH & PARSING FUNCTIONS
# =============================================================================
def calc_distance(p1, p2):
    return np.linalg.norm(p1 - p2)


def calc_angle(p1, p2, p3):
    v1 = p1 - p2
    v2 = p3 - p2
    cos_theta = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return math.degrees(math.acos(cos_theta))


def calc_dihedral(p1, p2, p3, p4):
    b0 = -1.0 * (p2 - p1)
    b1 = p3 - p2
    b2 = p4 - p3
    b1_norm = np.linalg.norm(b1)
    if b1_norm > 0:
        b1 /= b1_norm
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return math.degrees(math.atan2(y, x))


def parse_trj_content(file_content):
    lines = file_content.splitlines()
    frame_xyz_list = []
    frame_coords = []
    energies = []
    i = 0

    while i < len(lines):
        line = lines[i].strip()
        if line.isdigit():
            natoms = int(line)
            frame_block = "\n".join(lines[i : i + natoms + 2])
            frame_xyz_list.append(frame_block)

            # Parse Energy from comment line (if present)
            comment = lines[i + 1] if (i + 1) < len(lines) else ""
            found_e = False
            for token in comment.split():
                try:
                    val = float(token)
                    if val < 0:  # Isolates Hartree energy
                        energies.append(val)
                        found_e = True
                        break
                except ValueError:
                    continue

            if not found_e:
                energies.append(0.0)  # Default dummy energy if missing

            # Parse Atomic Coordinates
            coords = []
            for j in range(i + 2, i + 2 + natoms):
                parts = lines[j].split()
                coords.append(
                    np.array(
                        [float(parts[1]), float(parts[2]), float(parts[3])],
                        dtype=float,
                    )
                )
            frame_coords.append(coords)

            i += natoms + 2
        else:
            i += 1

    return frame_xyz_list, frame_coords, energies, natoms


# =============================================================================
# 2. SIDEBAR CONTROLS
# =============================================================================
st.sidebar.title("lookatit")
st.sidebar.header("📁 Input & Configuration")

uploaded_file = st.sidebar.file_uploader(
    "Upload Trajectory or Single Frame (.trj / .xyz)", type=["trj", "xyz"]
)

if uploaded_file is not None:
    file_content = uploaded_file.getvalue().decode("utf-8")
    frame_xyz_list, frame_coords, energies, natoms = parse_trj_content(
        file_content
    )

    if not frame_xyz_list:
        st.error("Could not parse any valid XYZ coordinate blocks.")
        st.stop()

    # Relative energy calculation (kcal/mol)
    has_real_energies = any(e != 0.0 for e in energies)
    if has_real_energies:
        e_min = min(energies)
        rel_energies = [(e - e_min) * 627.509 for e in energies]
    else:
        rel_energies = [0.0] * len(energies)

    frames = list(range(1, len(frame_xyz_list) + 1))

    st.sidebar.success(f"Loaded {len(frames)} frame(s) ({natoms} atoms/frame).")

    # Geometry Selection
    st.sidebar.subheader("📐 Geometrical Measurement")
    indices_str = st.sidebar.text_input(
        "Atom Indices (0-indexed, comma separated)", value="0, 1"
    )
    show_labels = st.sidebar.checkbox("Show Atom Indexes", value=True)

    # Parse Indices
    try:
        indices = [
            int(x.strip())
            for x in indices_str.split(",")
            if x.strip() != ""
        ]
        valid_indices = [i for i in indices if 0 <= i < natoms]
    except ValueError:
        valid_indices = []

    num_pts = len(valid_indices)
    geom_mode = "None"
    geom_values = []
    unit = ""

    if num_pts == 2:
        geom_mode = "Distance"
        unit = "Å"
        a1, a2 = valid_indices
        geom_values = [
            calc_distance(frame_coords[f][a1], frame_coords[f][a2])
            for f in range(len(frames))
        ]
    elif num_pts == 3:
        geom_mode = "Angle"
        unit = "°"
        a1, a2, a3 = valid_indices
        geom_values = [
            calc_angle(
                frame_coords[f][a1], frame_coords[f][a2], frame_coords[f][a3]
            )
            for f in range(len(frames))
        ]
    elif num_pts == 4:
        geom_mode = "Dihedral"
        unit = "°"
        a1, a2, a3, a4 = valid_indices
        geom_values = [
            calc_dihedral(
                frame_coords[f][a1],
                frame_coords[f][a2],
                frame_coords[f][a3],
                frame_coords[f][a4],
            )
            for f in range(len(frames))
        ]

    # --- Handle Single-Frame vs Multi-Frame Navigation ---
    if len(frames) > 1:
        frame_num = st.slider(
            "Select Scan Step / Frame",
            min_value=1,
            max_value=len(frames),
            value=1,
            step=1,
        )
    else:
        frame_num = 1
        st.info("Single-frame structure loaded.")

    idx = frame_num - 1

    # --- Export Frame Feature ---
    st.sidebar.subheader("💾 Export Structure")
    current_frame_xyz = frame_xyz_list[idx]
    
    st.sidebar.download_button(
        label=f"Download Frame {frame_num} (.xyz)",
        data=current_frame_xyz,
        file_name=f"lookatit_frame_{frame_num}.xyz",
        mime="chemical/x-xyz",
        help="Export the active frame as a standalone XYZ coordinate file."
    )

    # =========================================================================
    # 3. RENDER DASHBOARD LAYOUT
    # =========================================================================
    col1, col2 = st.columns([1, 1])

    # --- Column 1: 3D Molecule Viewer ---
    with col1:
        st.subheader("3D Structure View")
        view = py3Dmol.view(width=500, height=400)
        view.addModel(frame_xyz_list[idx], "xyz")
        view.setStyle({"stick": {"radius": 0.15}, "sphere": {"scale": 0.25}})

        if show_labels:
            for atom_i, coord in enumerate(frame_coords[idx]):
                view.addLabel(
                    str(atom_i),
                    {
                        "position": {
                            "x": float(coord[0]),
                            "y": float(coord[1]),
                            "z": float(coord[2]),
                        },
                        "backgroundColor": "black",
                        "fontColor": "white",
                        "fontSize": 10,
                        "backgroundOpacity": 0.6,
                    },
                )

        if geom_mode != "None":
            pts = [frame_coords[idx][i] for i in valid_indices]
            for k in range(len(pts) - 1):
                view.addCylinder(
                    {
                        "start": {
                            "x": float(pts[k][0]),
                            "y": float(pts[k][1]),
                            "z": float(pts[k][2]),
                        },
                        "end": {
                            "x": float(pts[k + 1][0]),
                            "y": float(pts[k + 1][1]),
                            "z": float(pts[k + 1][2]),
                        },
                        "radius": 0.08,
                        "color": "yellow",
                        "fromCap": 1,
                        "toCap": 1,
                    }
                )

            label_pos = (
                pts[1] if geom_mode in ["Angle", "Dihedral"] else (pts[0] + pts[1]) / 2
            )
            val_str = f"{geom_values[idx]:.2f} {unit}"
            view.addLabel(
                val_str,
                {
                    "position": {
                        "x": float(label_pos[0]),
                        "y": float(label_pos[1]),
                        "z": float(label_pos[2]),
                    },
                    "backgroundColor": "yellow",
                    "fontColor": "black",
                    "fontSize": 12,
                },
            )

        view.setBackgroundColor("white")
        view.zoomTo()

        html_content = view._make_html()
        components.html(html_content, height=410, width=510)

    # --- Column 2: Dual-Axis Energy Plot ---
    with col2:
        st.subheader("Energy & Geometry Profile")
        fig, ax1 = plt.subplots(figsize=(6, 4))

        color_e = "#d62728"
        ax1.set_xlabel("Scan Step / Frame")
        ax1.set_ylabel("Relative Energy (kcal/mol)", color=color_e)
        ax1.plot(frames, rel_energies, color=color_e, linewidth=1.5, alpha=0.4)
        ax1.plot(
            frames[idx],
            rel_energies[idx],
            marker="o",
            color=color_e,
            markersize=8,
        )
        ax1.tick_params(axis="y", labelcolor=color_e)
        ax1.grid(True, linestyle="--", alpha=0.3)

        if geom_mode != "None":
            ax2 = ax1.twinx()
            color_g = "#1f77b4"
            idx_str = "-".join(map(str, valid_indices))
            ax2.set_ylabel(f"{geom_mode} [{idx_str}] ({unit})", color=color_g)
            ax2.plot(
                frames,
                geom_values,
                color=color_g,
                linestyle="-.",
                linewidth=1.5,
                alpha=0.6,
            )
            ax2.plot(
                frames[idx],
                geom_values[idx],
                marker="s",
                color=color_g,
                markersize=8,
            )
            ax2.tick_params(axis="y", labelcolor=color_g)

            title_str = f"Step {frame_num} | E: {rel_energies[idx]:.2f} kcal/mol | {geom_mode}: {geom_values[idx]:.2f}{unit}"
        else:
            title_str = (
                f"Step {frame_num} | Energy: {rel_energies[idx]:.2f} kcal/mol"
            )

        plt.title(title_str)
        plt.tight_layout()
        st.pyplot(fig)

else:
    st.info(
        "👈 Upload an ORCA `.trj` or single/multi-frame `.xyz` file in the sidebar to launch lookatit."
    )
