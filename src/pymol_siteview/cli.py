#!/usr/bin/env python
"""
Automated PyMOL site-view renderer for PDB/mmCIF/MOL2 and related structures.

Run with PyMOL, for example:
    pymol -cq pymol_siteview.py -- protein_ligand.pdb -o out

Or with a Python interpreter that can import the PyMOL module:
    python pymol_siteview.py protein_ligand.pdb -o out
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple


Vec3 = Tuple[float, float, float]

SUPPORTED_EXTENSIONS = {
    ".pdb",
    ".ent",
    ".cif",
    ".mmcif",
    ".mol2",
    ".sdf",
    ".mol",
}

PROTEIN_RING_ATOMS: Dict[str, List[Tuple[str, ...]]] = {
    "PHE": [("CG", "CD1", "CD2", "CE1", "CE2", "CZ")],
    "TYR": [("CG", "CD1", "CD2", "CE1", "CE2", "CZ")],
    "HIS": [("CG", "ND1", "CD2", "CE1", "NE2")],
    "HID": [("CG", "ND1", "CD2", "CE1", "NE2")],
    "HIE": [("CG", "ND1", "CD2", "CE1", "NE2")],
    "HIP": [("CG", "ND1", "CD2", "CE1", "NE2")],
    "TRP": [
        ("CG", "CD1", "NE1", "CE2", "CD2"),
        ("CD2", "CE2", "CE3", "CZ3", "CH2", "CZ2"),
    ],
}

PROTEIN_CATION_ATOMS: Set[Tuple[str, str]] = {
    ("LYS", "NZ"),
    ("ARG", "CZ"),
    ("HIP", "NE2"),
    ("HIP", "ND1"),
}

SOLVENT_NAMES = {"HOH", "WAT", "DOD", "H2O"}
AROMATIC_LIGAND_ELEMENTS = {"C", "N", "O", "S"}
LIGAND_CATION_ELEMENTS = {"N", "P", "S"}


@dataclass(frozen=True)
class Ring:
    source: str
    label: str
    selection: str
    center: Vec3
    normal: Vec3
    residue_key: Optional[Tuple[str, str, str, str]]
    residue_atom_index: Optional[int]


@dataclass(frozen=True)
class Cation:
    source: str
    label: str
    selection: str
    center: Vec3
    residue_key: Optional[Tuple[str, str, str, str]]
    residue_atom_index: Optional[int]


@dataclass
class RenderResult:
    input_file: Path
    pse_file: Path
    png_file: Optional[Path]
    ligand_atoms: int
    site_residues: int
    hbond_object: Optional[str]
    pi_pi_count: int
    pi_cation_count: int


def strip_pymol_args(argv: Sequence[str]) -> List[str]:
    """Return only script args when called as `pymol script.py -- args`."""
    args = list(argv)
    if "--" in args:
        return args[args.index("--") + 1 :]
    return args


def import_pymol(headless: bool):
    try:
        import pymol  # type: ignore
        from pymol import cmd  # type: ignore
    except ImportError as exc:
        raise SystemExit(
            "PyMOL Python API is required. Run with `pymol -cq pymol_siteview.py -- ...` "
            "or install a Python environment that provides `import pymol`."
        ) from exc

    try:
        cmd.get_version()
    except Exception:
        launch_args = ["pymol", "-cq"] if headless else ["pymol"]
        pymol.finish_launching(launch_args)

    return cmd


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create PyMOL site-view PNG/PSE files with protein-ligand interactions."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        type=Path,
        help="Structure files, e.g. .pdb, .cif, .mol2, .sdf.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=Path("siteview_out"),
        help="Output directory. Default: siteview_out",
    )
    parser.add_argument(
        "--ligand-selection",
        default=None,
        help="Optional PyMOL selection for ligand atoms, e.g. `resn LIG`.",
    )
    parser.add_argument(
        "--protein-selection",
        default="polymer.protein",
        help="PyMOL protein selection inside each object. Default: polymer.protein",
    )
    parser.add_argument(
        "--site-cutoff",
        type=float,
        default=4.5,
        help="Residues within this distance from ligand are shown/labeled. Default: 4.5",
    )
    parser.add_argument(
        "--hbond-cutoff",
        type=float,
        default=3.6,
        help="PyMOL polar-contact cutoff for hydrogen bonds. Default: 3.6",
    )
    parser.add_argument(
        "--pi-cutoff",
        type=float,
        default=5.5,
        help="Ring-center cutoff for pi-pi interactions. Default: 5.5",
    )
    parser.add_argument(
        "--pi-cation-cutoff",
        type=float,
        default=6.0,
        help="Ring-center to cation cutoff for pi-cation interactions. Default: 6.0",
    )
    parser.add_argument(
        "--cartoon-transparency",
        type=float,
        default=0.8,
        help="Protein cartoon transparency when ligand exists. Default: 0.8",
    )
    parser.add_argument(
        "--stick-radius",
        type=float,
        default=0.2,
        help="Stick radius for ligand and site residues. Default: 0.2",
    )
    parser.add_argument(
        "--label-font-id",
        type=int,
        default=7,
        help="PyMOL label_font_id. Default: 7, commonly used for Arial-like labels.",
    )
    parser.add_argument(
        "--label-size",
        type=float,
        default=18.0,
        help="PyMOL label size. Default: 18",
    )
    parser.add_argument(
        "--width",
        type=int,
        default=1800,
        help="PNG width in pixels. Default: 1800",
    )
    parser.add_argument(
        "--height",
        type=int,
        default=1400,
        help="PNG height in pixels. Default: 1400",
    )
    parser.add_argument(
        "--no-png",
        action="store_true",
        help="Only save PSE, do not render PNG.",
    )
    parser.add_argument(
        "--no-ray",
        action="store_true",
        help="Use OpenGL snapshot for PNG instead of ray tracing.",
    )
    parser.add_argument(
        "--keep-solvent",
        action="store_true",
        help="Keep solvent molecules in the session.",
    )
    parser.add_argument(
        "--no-quit",
        action="store_true",
        help="Do not call cmd.quit() after processing.",
    )
    return parser.parse_args(strip_pymol_args(argv))


def sanitize_name(value: str, fallback: str = "obj") -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_]+", "_", value).strip("_")
    if not cleaned:
        cleaned = fallback
    if cleaned[0].isdigit():
        cleaned = f"{fallback}_{cleaned}"
    return cleaned[:60]


def unique_name(base: str, used: Set[str]) -> str:
    name = sanitize_name(base)
    candidate = name
    counter = 2
    while candidate in used:
        candidate = f"{name}_{counter}"
        counter += 1
    used.add(candidate)
    return candidate


def add_vec(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def sub_vec(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def dot_vec(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def norm_vec(a: Vec3) -> float:
    return math.sqrt(dot_vec(a, a))


def scale_vec(a: Vec3, value: float) -> Vec3:
    return (a[0] * value, a[1] * value, a[2] * value)


def dist(a: Vec3, b: Vec3) -> float:
    return norm_vec(sub_vec(a, b))


def unit_vec(a: Vec3) -> Vec3:
    length = norm_vec(a)
    if length < 1.0e-6:
        return (0.0, 0.0, 1.0)
    return (a[0] / length, a[1] / length, a[2] / length)


def center_of(coords: Sequence[Vec3]) -> Vec3:
    total = (0.0, 0.0, 0.0)
    for coord in coords:
        total = add_vec(total, coord)
    count = float(len(coords) or 1)
    return (total[0] / count, total[1] / count, total[2] / count)


def plane_normal(coords: Sequence[Vec3]) -> Vec3:
    if len(coords) < 3:
        return (0.0, 0.0, 1.0)

    nx = ny = nz = 0.0
    for index, coord in enumerate(coords):
        nxt = coords[(index + 1) % len(coords)]
        nx += (coord[1] - nxt[1]) * (coord[2] + nxt[2])
        ny += (coord[2] - nxt[2]) * (coord[0] + nxt[0])
        nz += (coord[0] - nxt[0]) * (coord[1] + nxt[1])
    return unit_vec((nx, ny, nz))


def angle_between_normals(a: Vec3, b: Vec3) -> float:
    cosine = max(-1.0, min(1.0, abs(dot_vec(unit_vec(a), unit_vec(b)))))
    return math.degrees(math.acos(cosine))


def atom_coord(atom) -> Vec3:
    return (float(atom.coord[0]), float(atom.coord[1]), float(atom.coord[2]))


def atom_element(atom) -> str:
    elem = (getattr(atom, "elem", "") or "").strip().upper()
    if elem:
        return elem
    letters = "".join(ch for ch in getattr(atom, "name", "") if ch.isalpha())
    return letters[:1].upper()


def residue_key(atom) -> Tuple[str, str, str, str]:
    return (
        (getattr(atom, "segi", "") or "").strip(),
        (getattr(atom, "chain", "") or "").strip(),
        (getattr(atom, "resi", "") or "").strip(),
        (getattr(atom, "resn", "") or "").strip().upper(),
    )


def residue_label_from_key(key: Tuple[str, str, str, str]) -> str:
    _segi, chain, resi, resn = key
    label = f"{resn}{resi}"
    if chain:
        label = f"{label}:{chain}"
    return label


def atom_index_selection(obj_name: str, atom_index: int) -> str:
    return f"({obj_name} and index {int(atom_index)})"


def atom_indices_selection(obj_name: str, atom_indices: Iterable[int]) -> str:
    indices = sorted({int(index) for index in atom_indices})
    if not indices:
        return f"({obj_name} and none)"
    return f"({obj_name} and index {'+'.join(str(index) for index in indices)})"


def count_residues(cmd, selection: str) -> int:
    model = cmd.get_model(selection)
    return len({residue_key(atom) for atom in model.atom})


def set_if_exists(cmd, setting: str, value, selection: Optional[str] = None) -> None:
    try:
        if selection is None:
            cmd.set(setting, value)
        else:
            cmd.set(setting, value, selection)
    except Exception:
        pass


def choose_ligand_selection(cmd, obj_name: str, args: argparse.Namespace) -> str:
    if args.ligand_selection:
        return f"({obj_name} and ({args.ligand_selection}))"

    candidates = [
        f"({obj_name} and organic and not solvent and not polymer)",
        f"({obj_name} and organic and not solvent)",
        f"({obj_name} and hetatm and not solvent and not polymer)",
    ]
    for selection in candidates:
        try:
            if cmd.count_atoms(selection) > 0:
                return selection
        except Exception:
            continue

    if cmd.count_atoms(f"({obj_name} and not solvent)") > 0 and cmd.count_atoms(
        f"({obj_name} and polymer.protein)"
    ) == 0:
        return f"({obj_name} and not solvent)"

    return f"({obj_name} and none)"


def style_scene(cmd, obj_name: str, protein_sel: str, ligand_sel: str, site_sel: str, args) -> None:
    cmd.bg_color("white")
    cmd.hide("everything", obj_name)
    set_if_exists(cmd, "orthoscopic", "on")
    set_if_exists(cmd, "ray_opaque_background", "off")
    set_if_exists(cmd, "antialias", 2)
    set_if_exists(cmd, "ambient", 0.45)
    set_if_exists(cmd, "spec_reflect", 0.15)
    set_if_exists(cmd, "stick_radius", args.stick_radius)
    set_if_exists(cmd, "dash_gap", 0.28)
    set_if_exists(cmd, "dash_width", 3.0)
    set_if_exists(cmd, "dash_radius", 0.04)
    set_if_exists(cmd, "label_font_id", args.label_font_id)
    set_if_exists(cmd, "label_size", args.label_size)
    set_if_exists(cmd, "label_color", "black")
    set_if_exists(cmd, "label_outline_color", "white")

    try:
        cmd.set_color("site_carbon", [0.20, 0.48, 0.82])
        cmd.set_color("ligand_carbon", [0.95, 0.52, 0.18])
        cmd.set_color("protein_cartoon", [0.78, 0.80, 0.82])
    except Exception:
        pass

    if cmd.count_atoms(protein_sel) > 0:
        cmd.show("cartoon", protein_sel)
        cmd.color("protein_cartoon", protein_sel)
        if cmd.count_atoms(ligand_sel) > 0:
            set_if_exists(cmd, "cartoon_transparency", args.cartoon_transparency, protein_sel)

    if cmd.count_atoms(site_sel) > 0:
        cmd.show("sticks", site_sel)
        cmd.color("site_carbon", f"({site_sel}) and elem C")

    if cmd.count_atoms(ligand_sel) > 0:
        cmd.show("sticks", ligand_sel)
        cmd.color("ligand_carbon", f"({ligand_sel}) and elem C")

    for color, element in [
        ("blue", "N"),
        ("red", "O"),
        ("yelloworange", "S"),
        ("orange", "P"),
    ]:
        cmd.color(color, f"({obj_name}) and elem {element}")

    if not args.keep_solvent:
        cmd.hide("everything", f"({obj_name} and solvent)")


def label_site_residues(cmd, obj_name: str, site_sel: str) -> None:
    model = cmd.get_model(site_sel)
    residues: Dict[Tuple[str, str, str, str], int] = {}
    ca_indices: Dict[Tuple[str, str, str, str], int] = {}
    for atom in model.atom:
        key = residue_key(atom)
        residues.setdefault(key, int(atom.index))
        if getattr(atom, "name", "").strip().upper() == "CA":
            ca_indices[key] = int(atom.index)

    for key, fallback_index in sorted(residues.items(), key=lambda item: residue_label_from_key(item[0])):
        atom_index = ca_indices.get(key, fallback_index)
        label = residue_label_from_key(key)
        label_sel = atom_index_selection(obj_name, atom_index)
        try:
            cmd.label(label_sel, json.dumps(label))
        except Exception:
            continue


def get_protein_rings(cmd, obj_name: str, protein_sel: str) -> List[Ring]:
    model = cmd.get_model(protein_sel)
    residues: Dict[Tuple[str, str, str, str], Dict[str, object]] = {}
    ca_by_residue: Dict[Tuple[str, str, str, str], int] = {}

    for atom in model.atom:
        key = residue_key(atom)
        name = getattr(atom, "name", "").strip().upper()
        residues.setdefault(key, {})[name] = atom
        if name == "CA":
            ca_by_residue[key] = int(atom.index)

    rings: List[Ring] = []
    for key, atoms_by_name in residues.items():
        resn = key[3]
        for ring_atoms in PROTEIN_RING_ATOMS.get(resn, []):
            if not all(name in atoms_by_name for name in ring_atoms):
                continue
            atoms = [atoms_by_name[name] for name in ring_atoms]
            coords = [atom_coord(atom) for atom in atoms]
            atom_indices = [int(atom.index) for atom in atoms]
            rings.append(
                Ring(
                    source="protein",
                    label=residue_label_from_key(key),
                    selection=atom_indices_selection(obj_name, atom_indices),
                    center=center_of(coords),
                    normal=plane_normal(coords),
                    residue_key=key,
                    residue_atom_index=ca_by_residue.get(key, atom_indices[0]),
                )
            )
    return rings


def build_bond_graph(model) -> Dict[int, Set[int]]:
    graph: Dict[int, Set[int]] = {index: set() for index in range(len(model.atom))}
    for bond in getattr(model, "bond", []):
        try:
            left, right = int(bond.index[0]), int(bond.index[1])
        except Exception:
            continue
        if left in graph and right in graph:
            graph[left].add(right)
            graph[right].add(left)
    return graph


def canonical_cycle(path: Sequence[int]) -> Tuple[int, ...]:
    return tuple(sorted(path))


def find_cycles_of_size(graph: Dict[int, Set[int]], sizes: Set[int]) -> Set[Tuple[int, ...]]:
    if not graph:
        return set()

    max_size = max(sizes)
    cycles: Set[Tuple[int, ...]] = set()

    for start in sorted(graph):
        stack: List[Tuple[int, List[int]]] = [(start, [start])]
        while stack:
            current, path = stack.pop()
            if len(path) > max_size:
                continue
            for neighbor in graph[current]:
                if neighbor == start and len(path) in sizes:
                    cycles.add(canonical_cycle(path))
                elif neighbor not in path and neighbor >= start and len(path) < max_size:
                    stack.append((neighbor, path + [neighbor]))
    return cycles


def get_ligand_rings(cmd, obj_name: str, ligand_sel: str) -> List[Ring]:
    model = cmd.get_model(ligand_sel)
    if not getattr(model, "bond", None):
        return []

    graph = build_bond_graph(model)
    cycles = find_cycles_of_size(graph, {5, 6})
    rings: List[Ring] = []

    for cycle in sorted(cycles):
        atoms = [model.atom[index] for index in cycle]
        elements = {atom_element(atom) for atom in atoms}
        if not elements.issubset(AROMATIC_LIGAND_ELEMENTS):
            continue
        coords = [atom_coord(atom) for atom in atoms]
        atom_indices = [int(atom.index) for atom in atoms]
        labels = [residue_label_from_key(residue_key(atom)) for atom in atoms]
        label = labels[0] if labels else "ligand"
        rings.append(
            Ring(
                source="ligand",
                label=label,
                selection=atom_indices_selection(obj_name, atom_indices),
                center=center_of(coords),
                normal=plane_normal(coords),
                residue_key=None,
                residue_atom_index=None,
            )
        )
    return rings


def atom_formal_or_partial_charge(atom) -> float:
    for attr in ("formal_charge", "partial_charge", "charge"):
        value = getattr(atom, attr, None)
        if value in (None, ""):
            continue
        try:
            return float(value)
        except (TypeError, ValueError):
            continue
    return 0.0


def get_protein_cations(cmd, obj_name: str, protein_sel: str) -> List[Cation]:
    model = cmd.get_model(protein_sel)
    cations: List[Cation] = []
    for atom in model.atom:
        resn = (getattr(atom, "resn", "") or "").strip().upper()
        name = (getattr(atom, "name", "") or "").strip().upper()
        if (resn, name) not in PROTEIN_CATION_ATOMS:
            continue
        key = residue_key(atom)
        cations.append(
            Cation(
                source="protein",
                label=f"{residue_label_from_key(key)} {name}",
                selection=atom_index_selection(obj_name, int(atom.index)),
                center=atom_coord(atom),
                residue_key=key,
                residue_atom_index=int(atom.index),
            )
        )
    return cations


def get_ligand_cations(cmd, obj_name: str, ligand_sel: str) -> List[Cation]:
    model = cmd.get_model(ligand_sel)
    cations: List[Cation] = []
    for atom in model.atom:
        elem = atom_element(atom)
        charge = atom_formal_or_partial_charge(atom)
        if charge <= 0.2 and elem not in LIGAND_CATION_ELEMENTS:
            continue
        if elem == "N" and charge < -0.1:
            continue
        label = residue_label_from_key(residue_key(atom))
        cations.append(
            Cation(
                source="ligand",
                label=f"{label} {getattr(atom, 'name', '').strip()}",
                selection=atom_index_selection(obj_name, int(atom.index)),
                center=atom_coord(atom),
                residue_key=None,
                residue_atom_index=None,
            )
        )
    return cations


def create_center_pseudoatom(cmd, name: str, center: Vec3, color: str) -> str:
    cmd.pseudoatom(name, pos=list(center), color=color)
    cmd.hide("everything", name)
    return name


def create_dashed_distance(cmd, name: str, sel1: str, sel2: str, color: str) -> str:
    cmd.distance(name, sel1, sel2)
    cmd.hide("labels", name)
    cmd.color(color, name)
    set_if_exists(cmd, "dash_color", color, name)
    set_if_exists(cmd, "dash_width", 3.0, name)
    set_if_exists(cmd, "dash_gap", 0.25, name)
    set_if_exists(cmd, "dash_radius", 0.045, name)
    return name


def draw_hbonds(cmd, obj_prefix: str, ligand_sel: str, site_sel: str, cutoff: float) -> Optional[str]:
    if cmd.count_atoms(ligand_sel) == 0 or cmd.count_atoms(site_sel) == 0:
        return None
    name = f"{obj_prefix}_hbonds"
    try:
        cmd.distance(name, ligand_sel, site_sel, cutoff=cutoff, mode=2)
        cmd.hide("labels", name)
        cmd.color("yellow", name)
        set_if_exists(cmd, "dash_color", "yellow", name)
        set_if_exists(cmd, "dash_width", 3.2, name)
        set_if_exists(cmd, "dash_gap", 0.22, name)
        set_if_exists(cmd, "dash_radius", 0.045, name)
        return name
    except Exception:
        return None


def draw_pi_pi(
    cmd,
    obj_prefix: str,
    ligand_rings: Sequence[Ring],
    protein_rings: Sequence[Ring],
    cutoff: float,
) -> Tuple[int, Set[Tuple[str, str, str, str]]]:
    count = 0
    residues: Set[Tuple[str, str, str, str]] = set()

    for ligand_ring in ligand_rings:
        for protein_ring in protein_rings:
            center_distance = dist(ligand_ring.center, protein_ring.center)
            if center_distance > cutoff:
                continue
            angle = angle_between_normals(ligand_ring.normal, protein_ring.normal)
            if not (angle <= 35.0 or 55.0 <= angle <= 90.0):
                continue

            count += 1
            left = create_center_pseudoatom(
                cmd, f"{obj_prefix}_pipi_lig_{count}", ligand_ring.center, "blue"
            )
            right = create_center_pseudoatom(
                cmd, f"{obj_prefix}_pipi_pro_{count}", protein_ring.center, "blue"
            )
            create_dashed_distance(cmd, f"{obj_prefix}_pipi_{count}", left, right, "blue")
            if protein_ring.residue_key is not None:
                residues.add(protein_ring.residue_key)

    return count, residues


def cation_projects_near_ring(ring: Ring, cation: Cation, max_plane_distance: float = 2.8) -> bool:
    ring_to_cation = sub_vec(cation.center, ring.center)
    plane_distance = abs(dot_vec(ring_to_cation, ring.normal))
    return plane_distance <= max_plane_distance


def draw_pi_cation(
    cmd,
    obj_prefix: str,
    ligand_rings: Sequence[Ring],
    protein_rings: Sequence[Ring],
    ligand_cations: Sequence[Cation],
    protein_cations: Sequence[Cation],
    cutoff: float,
) -> Tuple[int, Set[Tuple[str, str, str, str]]]:
    count = 0
    residues: Set[Tuple[str, str, str, str]] = set()

    pairs: List[Tuple[Ring, Cation]] = []
    pairs.extend((ring, cation) for ring in ligand_rings for cation in protein_cations)
    pairs.extend((ring, cation) for ring in protein_rings for cation in ligand_cations)

    for ring, cation in pairs:
        if dist(ring.center, cation.center) > cutoff:
            continue
        if not cation_projects_near_ring(ring, cation):
            continue

        count += 1
        ring_atom = create_center_pseudoatom(
            cmd, f"{obj_prefix}_picat_ring_{count}", ring.center, "green"
        )
        create_dashed_distance(cmd, f"{obj_prefix}_picat_{count}", ring_atom, cation.selection, "green")

        if ring.residue_key is not None:
            residues.add(ring.residue_key)
        if cation.residue_key is not None:
            residues.add(cation.residue_key)

    return count, residues


def label_extra_residues(
    cmd,
    obj_name: str,
    residues: Iterable[Tuple[str, str, str, str]],
    representative_indices: Dict[Tuple[str, str, str, str], int],
) -> None:
    for key in residues:
        atom_index = representative_indices.get(key)
        if atom_index is None:
            continue
        try:
            cmd.label(atom_index_selection(obj_name, atom_index), json.dumps(residue_label_from_key(key)))
        except Exception:
            continue


def representative_residue_indices(cmd, protein_sel: str) -> Dict[Tuple[str, str, str, str], int]:
    model = cmd.get_model(protein_sel)
    reps: Dict[Tuple[str, str, str, str], int] = {}
    for atom in model.atom:
        key = residue_key(atom)
        name = (getattr(atom, "name", "") or "").strip().upper()
        if key not in reps or name == "CA":
            reps[key] = int(atom.index)
    return reps


def orient_site(cmd, view_sel: str) -> None:
    try:
        cmd.orient(view_sel)
        cmd.zoom(view_sel, buffer=3.5)
        cmd.center(view_sel)
        cmd.clip("near", -8)
        cmd.clip("far", 20)
    except Exception:
        cmd.zoom("all")


def save_outputs(cmd, pse_file: Path, png_file: Optional[Path], args: argparse.Namespace) -> None:
    cmd.scene("siteview", "store")
    cmd.save(str(pse_file))
    if png_file is not None:
        cmd.png(str(png_file), width=args.width, height=args.height, dpi=300, ray=0 if args.no_ray else 1)


def process_structure(cmd, structure_file: Path, outdir: Path, args, used_names: Set[str]) -> RenderResult:
    if not structure_file.exists():
        raise FileNotFoundError(structure_file)
    if structure_file.suffix.lower() not in SUPPORTED_EXTENSIONS:
        print(f"Warning: {structure_file} has an uncommon extension; PyMOL will still try to load it.")

    cmd.reinitialize()

    obj_name = unique_name(structure_file.stem, used_names)
    obj_prefix = sanitize_name(obj_name, "site")
    cmd.load(str(structure_file.resolve()), obj_name)
    cmd.remove(f"({obj_name} and alt B+C+D+E+F)")

    protein_sel = f"({obj_name} and ({args.protein_selection}))"
    ligand_sel = choose_ligand_selection(cmd, obj_name, args)
    ligand_atoms = cmd.count_atoms(ligand_sel)

    if ligand_atoms > 0 and cmd.count_atoms(protein_sel) > 0:
        site_sel = f"(byres (({protein_sel}) within {args.site_cutoff} of ({ligand_sel})))"
        view_sel = f"(({ligand_sel}) or ({site_sel}))"
    elif ligand_atoms > 0:
        site_sel = f"({obj_name} and none)"
        view_sel = ligand_sel
    else:
        site_sel = protein_sel
        view_sel = obj_name

    style_scene(cmd, obj_name, protein_sel, ligand_sel, site_sel, args)
    site_residue_count = count_residues(cmd, site_sel) if cmd.count_atoms(site_sel) > 0 else 0

    hbond_object = draw_hbonds(cmd, obj_prefix, ligand_sel, site_sel, args.hbond_cutoff)

    representative_indices = representative_residue_indices(cmd, protein_sel)
    label_site_residues(cmd, obj_name, site_sel)

    ligand_rings = get_ligand_rings(cmd, obj_name, ligand_sel) if ligand_atoms > 0 else []
    protein_rings = get_protein_rings(cmd, obj_name, protein_sel)
    ligand_cations = get_ligand_cations(cmd, obj_name, ligand_sel) if ligand_atoms > 0 else []
    protein_cations = get_protein_cations(cmd, obj_name, protein_sel)

    pi_pi_count, pi_pi_residues = draw_pi_pi(
        cmd, obj_prefix, ligand_rings, protein_rings, args.pi_cutoff
    )
    pi_cation_count, pi_cation_residues = draw_pi_cation(
        cmd,
        obj_prefix,
        ligand_rings,
        protein_rings,
        ligand_cations,
        protein_cations,
        args.pi_cation_cutoff,
    )
    label_extra_residues(cmd, obj_name, pi_pi_residues | pi_cation_residues, representative_indices)

    orient_site(cmd, view_sel)

    outdir.mkdir(parents=True, exist_ok=True)
    pse_file = outdir / f"{structure_file.stem}_siteview.pse"
    png_file = None if args.no_png else outdir / f"{structure_file.stem}_siteview.png"
    save_outputs(cmd, pse_file, png_file, args)

    return RenderResult(
        input_file=structure_file,
        pse_file=pse_file,
        png_file=png_file,
        ligand_atoms=ligand_atoms,
        site_residues=site_residue_count,
        hbond_object=hbond_object,
        pi_pi_count=pi_pi_count,
        pi_cation_count=pi_cation_count,
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    cmd = import_pymol(headless=True)

    results: List[RenderResult] = []
    used_names: Set[str] = set()

    for input_file in args.inputs:
        try:
            result = process_structure(cmd, input_file, args.outdir, args, used_names)
        except Exception as exc:
            print(f"[ERROR] {input_file}: {exc}", file=sys.stderr)
            continue
        results.append(result)
        png_text = f", PNG={result.png_file}" if result.png_file else ""
        print(
            f"[OK] {result.input_file} -> PSE={result.pse_file}{png_text}; "
            f"ligand_atoms={result.ligand_atoms}, site_residues={result.site_residues}, "
            f"pi_pi={result.pi_pi_count}, pi_cation={result.pi_cation_count}"
        )

    if not args.no_quit:
        try:
            cmd.quit()
        except Exception:
            pass

    return 0 if results else 1


if __name__ == "__main__":
    raise SystemExit(main())
