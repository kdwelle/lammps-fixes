# lammps-fixes

The fixes used in the paper here: https://pubs.acs.org/doi/full/10.1021/acs.jpcc.9b06635

Contact kdwelle@mit.edu with questions

**Fix_imagecharges**

Syntax

_fix ID group-ID imagecharges px py pz nx ny nz itype keyword value..._

* ID, group-ID are documented in fix command
* imagecharges = style name of this fix command
* px, py, pz = coordinates of a point on the image plane
* nx, ny, nz = define a vector normal to the image plane
* itype = atom type to be used as the image charges (note, any existing atoms of this type will be used first then new atoms of this type will be created)
* keyword = region or scale


**Fix_electrodeboundaries**

Syntax

_fix ID group-ID electrodeboundaries xlo dist v0 dv etype seed keyword value..._

* ID, group-ID are documented in fix command
* electrodeboundaries = style name of this fix command
* xlo = x coordinate of left electrode
* dist = distance between left and right electrodes
* v0
* dv
* etype
* seed
* zero or more keyword/value pairs may be appended
- keyword = ncycles or intercalation or porus
- ncycles value = float
  - float = number of redox attempts per function call, can be less than 1
- couple = neutralType
  - neutralType = atom type to be used as the neutral redox couple species
- porusLeft = xstart
  - xstart = start of the uniform region of redox
- porusRight = xend
	- xend = end of the uniform region of redox
- intercalation = neutralType
  - neutralType = atom type to be used as the neutral redox couple species, -1 for implicit intercalation
- overpotential = float
  - float = energy added to the dE term before a hopping probability is computed
- occupation = spacing
  - Use this keyword to add an occupation matrix which will not let an intercalation occur in the same spot twice without a deintercalation. 
  - spacing = distance between occupation sites, i.e. lattice spacing for the occupation matrix
- v0Increment = float
  - This is the amount by which v0 is changed for each electron tranfer. Higher values will punish deviations from nuetrality more severely. A negative value will lead to an unstable system with runaway charge buildup.
