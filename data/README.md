# SPDX-FileCopyrightText: Copyright © Timo Koch
# SPDX-License-Identifier: CC-BY-4.0

ukb-atlas data
-----------------

Using https://github.com/ComputationalPhysiology/ukb-atlas

@misc{Finsberg2025,
  doi = {10.5281/ZENODO.14840907},
  url = {https://zenodo.org/doi/10.5281/zenodo.14840907},
  author = {Finsberg,  Henrik and Pankewitz,  Lisa R},
  title = {UK Biobank atlas - mesh generation},
  publisher = {Zenodo},
  year = {2025},
  copyright = {MIT License}
}

Cite as:
--------------

Finsberg, H., & Pankewitz, L. R. (2025). UK Biobank atlas - mesh generation (v0.2.0). Zenodo. https://doi.org/10.5281/zenodo.14840907

Meshes from https://www.cardiacatlas.org/biventricular-modes/

Generate meshes:
-------------------

```
ukb-atlas data --mesh
```


Remeshed with gmsh:
------------------------

```
gmsh -2 -format msh2 heart.geo
```

Origin of coronary arteries:
-------------------------------

The left and right ventricle originate at the base of the aorta.
Using ParaView, we roughly selected the positions coordinates as

```
-14.6428 21.1399 37.5829
-23.9962 9.63808 23.1206
```
