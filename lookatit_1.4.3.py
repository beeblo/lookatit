import math
import os
import matplotlib.pyplot as plt
import numpy as np
import streamlit as st
import streamlit.components.v1 as components

# Set page layout and browser tab title
st.set_page_config(layout="wide", page_title="lookatit - ORCA Visualizer")

st.title("🧪 lookatit: ORCA Trajectory & Geometry Visualizer")

# =============================================================================
# 0. DECLARE CUSTOM 3DMOL STREAMLIT COMPONENT (BIDIRECTIONAL CLICK BRIDGE)
# =============================================================================
COMPONENT_DIR = os.path.join(os.path.dirname(__file__), ".mol_viewer_component")
os.makedirs(COMPONENT_DIR, exist_ok=True)

INDEX_HTML = """<!DOCTYPE html>
<html>
<head>
    <script src="https://code.jquery.com/jquery-3.6.0.min.js"></script>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <style>
        body, html { margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden; }
        #gviewer { width: 100%; height: 400px; position: relative; }
    </style>
</head>
<body>
    <div id="gviewer"></div>
    <script>
        function sendToStreamlit(type, data) {
            window.parent.postMessage(Object.assign({isStreamlitMessage: true, type: type}, data), "*");
        }

        var viewer = null;

        window.addEventListener("message", function(event) {
            if (event.data.type === "streamlit:render") {
                var args = event.data.args;
                var xyz = args.xyz;
                var labels = args.labels || [];
                var cylinders = args.cylinders || [];
                var geomLabel = args.geomLabel || null;
                var zoomTo = args.zoomTo;

                if (!viewer) {
                    viewer = $3Dmol.createViewer($("#gviewer"), {backgroundColor: "white"});
                    sendToStreamlit("streamlit:setFrameHeight", {height: 410});
                }

                viewer.clear();
                viewer.addModel(xyz, "xyz");
                viewer.setStyle({}, {stick: {radius: 0.15}, sphere: {scale: 0.25}});

                labels.forEach(function(l) {
                    viewer.addLabel(l.text, l.spec);
                });

                cylinders.forEach(function(c) {
                    viewer.addCylinder(c);
                });

                if (geomLabel) {
                    viewer.addLabel(geomLabel.text, geomLabel.spec);
                }

                viewer.setClickable({}, true, function(atom, viewer, event, container) {
                    if (atom && atom.index !== undefined) {
                        sendToStreamlit("streamlit:setComponentValue", {value: atom.index});
                    }
                });

                if (zoomTo) {
                    viewer.zoomTo();
                }
                
                viewer.render();
            }
        });

        sendToStreamlit("streamlit:componentReady", {apiVersion: 1});
    </script>
</body>
</html>
"""

index_path = os.path.join(COMPONENT_DIR, "index.html")
with open(index_path, "w", encoding="utf-8") as f:
    f.write(INDEX_HTML)

mol_viewer = components.declare_component("mol_viewer", path=COMPONENT_DIR)


# =============================================================================
# 1. VECTOR MATH & PARSING FUNCTIONS
# =============================================================================
def calc_distance(p1, p2):
    return np.linalg.norm(p1 - p2)


def calc_angle(p1, p2, p3):
    v1 = p1 - p2
    v2 = p3 - p2
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
    if norm1 == 0 or norm2 == 0:
        return 0.0
    cos_theta = np.dot(v1, v2) / (norm1 * norm2)
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
    frame_symbols = []
    frame_coords = []
    energies = []
    i = 0

    while i < len(lines):
        line = lines[i].strip()
        if line.isdigit():
            natoms = int(line)
            if i + 1 + natoms >= len(lines):
                break  # Guard against truncated trailing content

            comment = lines[i + 1] if (i + 1) < len(lines) else ""
            found_e = False
            for token in comment.split():
                try:
                    val = float(token)
                    if val < 0:
                        energies.append(val)
                        found_e = True
                        break
                except ValueError:
                    continue

            if not found_e:
                energies.append(0.0)

            coords = []
            symbols = []
            for j in range(i + 2, i + 2 + natoms):
                parts = lines[j].split()
                if len(parts) >= 4:
                    symbols.append(parts[0])
                    coords.append(
                        np.array(
                            [float(parts[1]), float(parts[2]), float(parts[3])],
                            dtype=float,
                        )
                    )

            if len(coords) == natoms:
                frame_coords.append(coords)
                frame_symbols.append(symbols)

            i += natoms + 2
        else:
            i += 1

    natoms_count = natoms if frame_coords else 0
    return frame_symbols, frame_coords, energies, natoms_count


def rebuild_xyz(symbols, coords, comment=""):
    lines = [str(len(symbols)), comment]
    for sym, c in zip(symbols, coords):
        lines.append(f"{sym:2s} {c[0]:12.6f} {c[1]:12.6f} {c[2]:12.6f}")
    return "\n".join(lines)


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
    frame_symbols, frame_coords, energies, natoms = parse_trj_content(
        file_content
    )

    if not frame_coords:
        st.error("Could not parse any valid XYZ coordinate blocks.")
        st.stop()

    # Session State Initialization
    file_id = f"{uploaded_file.name}_{uploaded_file.size}"
    if (
        "current_file_id" not in st.session_state
        or st.session_state.current_file_id != file_id
    ):
        st.session_state.current_file_id = file_id
        st.session_state.editable_coords = [
            [np.array(c, copy=True) for c in frame] for frame in frame_coords
        ]
        st.session_state.editable_symbols = frame_symbols
        st.session_state.selected_indices = [0, 1] if natoms >= 2 else [0]
        st.session_state.last_clicked = None
        st.session_state.should_zoom = True
        st.session_state.target_atom = 0

    # Relative energy calculation (kcal/mol)
    has_real_energies = any(e != 0.0 for e in energies)
    if has_real_energies:
        e_min = min(energies)
        rel_energies = [(e - e_min) * 627.509 for e in energies]
    else:
        rel_energies = [0.0] * len(energies)

    frames = list(range(1, len(frame_coords) + 1))
    st.sidebar.success(f"Loaded {len(frames)} frame(s) ({natoms} atoms/frame).")

    # Handle Frame Navigation
    if len(frames) > 1:
        frame_num = st.sidebar.slider(
            "Select Scan Step / Frame",
            min_value=1,
            max_value=len(frames),
            value=1,
            step=1,
        )
    else:
        frame_num = 1
        st.sidebar.info("Single-frame structure loaded.")

    idx = frame_num - 1

    current_symbols = st.session_state.editable_symbols[idx]
    current_coords = st.session_state.editable_coords[idx]

    # --- Geometrical Measurement Controls ---
    st.sidebar.subheader("📐 Geometrical Measurement")

    col_sel1, col_sel2 = st.sidebar.columns([3, 1])
    with col_sel2:
        if st.button("Clear"):
            st.session_state.selected_indices = []
            st.session_state.last_clicked = None
            st.session_state.should_zoom = False
            st.rerun()

    indices_str = st.sidebar.text_input(
        "Atom Indices (click 3D view or type)",
        value=", ".join(map(str, st.session_state.selected_indices)),
    )

    try:
        indices = [
            int(x.strip())
            for x in indices_str.split(",")
            if x.strip() != ""
        ]
        valid_indices = [i for i in indices if 0 <= i < natoms]
        st.session_state.selected_indices = valid_indices
    except ValueError:
        valid_indices = []

    show_labels = st.sidebar.checkbox("Show Atom Indexes", value=True)

    num_pts = len(valid_indices)
    geom_mode = "None"
    geom_values = []
    unit = ""

    if num_pts == 2:
        geom_mode = "Distance"
        unit = "Å"
        a1, a2 = valid_indices
        geom_values = [
            calc_distance(
                st.session_state.editable_coords[f][a1],
                st.session_state.editable_coords[f][a2],
            )
            for f in range(len(frames))
        ]
    elif num_pts == 3:
        geom_mode = "Angle"
        unit = "°"
        a1, a2, a3 = valid_indices
        geom_values = [
            calc_angle(
                st.session_state.editable_coords[f][a1],
                st.session_state.editable_coords[f][a2],
                st.session_state.editable_coords[f][a3],
            )
            for f in range(len(frames))
        ]
    elif num_pts == 4:
        geom_mode = "Dihedral"
        unit = "°"
        a1, a2, a3, a4 = valid_indices
        geom_values = [
            calc_dihedral(
                st.session_state.editable_coords[f][a1],
                st.session_state.editable_coords[f][a2],
                st.session_state.editable_coords[f][a3],
                st.session_state.editable_coords[f][a4],
            )
            for f in range(len(frames))
        ]

    # --- Atom Nudger Engine ---
    st.sidebar.subheader("🕹️ Atom Coordinate Nudger")

    target_atom = st.sidebar.number_input(
        "Target Atom Index to Move",
        min_value=0,
        max_value=natoms - 1,
        value=int(st.session_state.target_atom),
        step=1,
        help="Select any atom to translate along XYZ axes.",
    )
    st.session_state.target_atom = target_atom

    step_size = st.sidebar.number_input(
        "Step Size (Å)",
        min_value=0.001,
        max_value=2.0,
        value=0.100,
        step=0.025,
        format="%.3f",
    )

    col_x, col_y, col_z = st.sidebar.columns(3)

    with col_x:
        if st.button("-X"):
            st.session_state.editable_coords[idx][target_atom][0] -= step_size
            st.session_state.should_zoom = False
            st.rerun()
        if st.button("+X"):
            st.session_state.editable_coords[idx][target_atom][0] += step_size
            st.session_state.should_zoom = False
            st.rerun()

    with col_y:
        if st.button("-Y"):
            st.session_state.editable_coords[idx][target_atom][1] -= step_size
            st.session_state.should_zoom = False
            st.rerun()
        if st.button("+Y"):
            st.session_state.editable_coords[idx][target_atom][1] += step_size
            st.session_state.should_zoom = False
            st.rerun()

    with col_z:
        if st.button("-Z"):
            st.session_state.editable_coords[idx][target_atom][2] -= step_size
            st.session_state.should_zoom = False
            st.rerun()
        if st.button("+Z"):
            st.session_state.editable_coords[idx][target_atom][2] += step_size
            st.session_state.should_zoom = False
            st.rerun()

    # Rebuild XYZ string for current active frame
    active_xyz = rebuild_xyz(current_symbols, current_coords, f"Frame {frame_num}")

    # Export Structure
    st.sidebar.subheader("💾 Export Structure")
    st.sidebar.download_button(
        label=f"Download Frame {frame_num} (.xyz)",
        data=active_xyz,
        file_name=f"lookatit_frame_{frame_num}.xyz",
        mime="chemical/x-xyz",
    )

    # =========================================================================
    # 3. RENDER DASHBOARD LAYOUT
    # =========================================================================
    col1, col2 = st.columns([1, 1])

    # --- Column 1: 3D Molecule Viewer ---
    with col1:
        st.subheader("3D Structure View")

        atom_labels = []
        if show_labels:
            for atom_i, coord in enumerate(current_coords):
                atom_labels.append(
                    {
                        "text": str(atom_i),
                        "spec": {
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
                    }
                )

        cylinders_data = []
        geom_label_data = None

        if geom_mode != "None":
            pts = [current_coords[i] for i in valid_indices]
            for k in range(len(pts) - 1):
                cylinders_data.append(
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
            geom_label_data = {
                "text": val_str,
                "spec": {
                    "position": {
                        "x": float(label_pos[0]),
                        "y": float(label_pos[1]),
                        "z": float(label_pos[2]),
                    },
                    "backgroundColor": "yellow",
                    "fontColor": "black",
                    "fontSize": 12,
                },
            }

        # Render 3Dmol viewer component & catch click events
        clicked_atom = mol_viewer(
            xyz=active_xyz,
            labels=atom_labels,
            cylinders=cylinders_data,
            geomLabel=geom_label_data,
            zoomTo=st.session_state.should_zoom,
            key="3d_canvas",
        )

        # Handle Click Selection Events
        if (
            clicked_atom is not None
            and clicked_atom != st.session_state.last_clicked
        ):
            st.session_state.last_clicked = clicked_atom
            st.session_state.target_atom = clicked_atom
            st.session_state.should_zoom = False
            current_sel = list(st.session_state.selected_indices)

            if clicked_atom not in current_sel:
                if len(current_sel) >= 4:
                    st.session_state.selected_indices = [clicked_atom]
                else:
                    st.session_state.selected_indices.append(clicked_atom)
                st.rerun()

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
        plt.close(fig)  # Prevent Matplotlib memory leaks across reruns

else:
    st.info(
        "👈 Upload an ORCA `.trj` or single/multi-frame `.xyz` file in the sidebar to launch lookatit."
    )
