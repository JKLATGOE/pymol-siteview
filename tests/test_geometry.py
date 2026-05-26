from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pymol_siteview.cli import angle_between_normals, center_of, dist, sanitize_name


def test_center_and_distance() -> None:
    assert center_of([(0.0, 0.0, 0.0), (2.0, 2.0, 2.0)]) == (1.0, 1.0, 1.0)
    assert dist((0.0, 0.0, 0.0), (0.0, 3.0, 4.0)) == 5.0


def test_angle_between_normals_uses_plane_equivalence() -> None:
    assert angle_between_normals((0.0, 0.0, 1.0), (0.0, 0.0, -1.0)) == 0.0


def test_sanitize_name() -> None:
    assert sanitize_name("3QIP ligand-site.pdb") == "obj_3QIP_ligand_site_pdb"
