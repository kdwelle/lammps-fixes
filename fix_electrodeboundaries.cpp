/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Kaitlyn Dwelle (MIT)
------------------------------------------------------------------------- */

#include <math.h>
#include "fix_electrodeboundaries.h"
#include "fix.h"


using namespace std;
using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixElectrodeBoundaries::FixElectrodeBoundaries(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
	pxstr(NULL), pystr(NULL), pzstr(NULL), 
	nxstr(NULL), nystr(NULL), nzstr(NULL), 
	idregion(NULL), diststr(NULL){

  dx = 0.5; //plus/minus search for ion in vicinity
  zCut = 2.0; //distance from electrode to check for electrochem

  if (narg < 13) error->all(FLERR,"Illegal fix electrodeboundaries command -- not enough arguments");

	// first three arguments define a point on the plane
  pxvalue = force->numeric(FLERR,arg[3]); // else just convert input to a number and store
  pyvalue = force->numeric(FLERR,arg[4]); 
  pzvalue = force->numeric(FLERR,arg[5]); 

  // next three arguments define a vector normal to the plane
  nxvalue = force->numeric(FLERR,arg[6]); 
  nyvalue = force->numeric(FLERR,arg[7]); 
  nzvalue = force->numeric(FLERR,arg[8]); 

  // Next argument is a distance between electrodes 
  distvalue = force->numeric(FLERR,arg[9]); 
  // Then voltage@defined plane and voltage difference between electrodes
  vvalue = force->numeric(FLERR,arg[10]); 
  dvvalue = force->numeric(FLERR,arg[11]); 
  active_type = force->inumeric(FLERR,arg[12]); //type of atom that is electrochemically active

	// optional arguments
	iregion = -1;
	idregion = NULL;
	scale = 1.0;

	int iarg = 13;	//start after madatory arguments
	while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) { //keyword = region
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix electrodeboundaries command"); //check to make sure there is a next argument
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix electrodeboundaries does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    
    } else error->all(FLERR,"Illegal fix electrodeboundaries command"); // not a recognized keyword
  }
	atom->add_callback(0);

}

FixElectrodeBoundaries::~FixElectrodeBoundaries(){
  // destructor -- free us pointer arrays

  delete [] pxstr;
  delete [] pystr;
  delete [] pzstr;
  delete [] nxstr;
  delete [] nystr;
  delete [] nzstr;
  delete [] vstr;
  delete [] dvstr;
  delete [] idregion;

  atom->delete_callback(id,0);

}

int FixElectrodeBoundaries::setmask(){
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

void FixElectrodeBoundaries::init(){
  // initial setup -- not sure what should go here
}

void FixElectrodeBoundaries::pre_exchange(){
  // for some number of attempts
  // pick an x(uniform),y(uniform),z(uniform to some cutoff) position
  // Check if ion is within dx --> If so attempt reduction
  // --> if not, attempt oxidation



  for (int i=0; i<)
}

void FixGCMC::attempt_atomic_insertion_full(double *coord){

  ninsertion_attempts += 1.0;
  double energy_before = energy_stored;

  // add atom
  atom->avec->create_atom(ngcmc_type,coord);
  int m = atom->nlocal - 1;

  // add to groups
  // optionally add to type-based groups

  atom->mask[m] = groupbitall;
  for (int igroup = 0; igroup < ngrouptypes; igroup++) {
    if (
      _type == grouptypes[igroup])
atom->mask[m] |= grouptypebits[igroup];
  }

  atom->v[m][0] = random_unequal->gaussian()*sigma;
  atom->v[m][1] = random_unequal->gaussian()*sigma;
  atom->v[m][2] = random_unequal->gaussian()*sigma;
  if (charge_flag) atom->q[m] = charge;
  modify->create_attribute(m);


  atom->natoms++;
  if (atom->tag_enable) {
    atom->tag_extend();
    if (atom->map_style) atom->map_init();
  }
  atom->nghost = 0;

  if (force->kspace) force->kspace->qsum_qsq();
  double energy_after = energy_full();

  if (random_equal->uniform() <
      zz*volume*exp(beta*(energy_before - energy_after))/(ngas+1)) {

    ninsertion_successes += 1.0;
    energy_stored = energy_after;
  } else {
    atom->natoms--;
    if (proc_flag) atom->nlocal--;
    if (force->kspace) force->kspace->qsum_qsq();
    energy_stored = energy_before;
  }
  update_gas_atoms_list();
}


/* ----------------------------------------------------------------------
   compute particle's interaction energy with the rest of the system
------------------------------------------------------------------------- */

double FixElectrodeBoundaries::energy(int i, int itype, tagint imolecule, double *coord)
{
  double delx,dely,delz,rsq;

  double **x = atom->x;
  int *type = atom->type;
  tagint *molecule = atom->molecule;
  int nall = atom->nlocal + atom->nghost;
  pair = force->pair;
  cutsq = force->pair->cutsq;

  double fpair = 0.0;
  double factor_coul = 1.0;
  double factor_lj = 1.0;

  double total_energy = 0.0;

  for (int j = 0; j < nall; j++) {

    if (i == j) continue;
    if (mode == MOLECULE)
      if (imolecule == molecule[j]) continue;

    delx = coord[0] - x[j][0];
    dely = coord[1] - x[j][1];
    delz = coord[2] - x[j][2];
    rsq = delx*delx + dely*dely + delz*delz;
    int jtype = type[j];

    if (rsq < cutsq[itype][jtype])
      total_energy +=
        pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
  }

  return total_energy;
}

