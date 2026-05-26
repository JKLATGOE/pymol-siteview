# Examples

Put your own `pdb`, `cif`, `mol2`, `sdf` or `mol` files in a working folder, then run:

```bash
pymol-siteview complex.pdb -o siteview_out
```

If automatic ligand detection chooses the wrong molecule, specify it manually:

```bash
pymol-siteview complex.pdb -o siteview_out --ligand-selection "resn LIG"
```

Batch run:

```bash
pymol-siteview *.pdb *.cif *.mol2 -o siteview_out
```
