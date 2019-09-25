# lammps-fixes

The fixes used in the paper here: https://pubs.acs.org/doi/full/10.1021/acs.jpcc.9b06635

Contact kdwelle@mit.edu with questions

Fix_imagecharges
Syntax
fix ID group-ID imagecharges px py pz nx ny nz itype keyword value…

*ID, group-ID are documented in fix command
*imagecharges = style name of this fix command
*px, py, pz = coordinates of a point on the image plane
*nx, ny, nz = define a vector normal to the image plane
*itype = atom type to be used as the image charges (note, any existing atoms of this type will be used first then new atoms of this type will be created)
*keyword = region or scale


