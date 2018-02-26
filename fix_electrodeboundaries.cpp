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
	idregion(NULL),{

  dr = 0.5; //plus/minus search for ion in vicinity
  xcut = 2.0; //distance from electrode to check for electrochem
  ncycles = 100; //number of attempts per timestep

  if (narg < 8) error->all(FLERR,"Illegal fix electrodeboundaries command -- not enough arguments");

  // electrodes have to lie along the x-axis
  // Next argument is a distance between electrodes 
  xlo = force->numeric(FLERR,arg[3]); 
  dist = force->numeric(FLERR,arg[4]); 
  xhi = xlo + dist;
  // Then voltage@defined plane and voltage difference between electrodes
  v0 = force->numeric(FLERR,arg[5]); 
  dv = force->numeric(FLERR,arg[6]); 
  etype = force->inumeric(FLERR,arg[7]); //type of atom that is electrochemically active

	// optional arguments
	iregion = -1;
	idregion = NULL;
	scale = 1.0;

	int iarg = 8;	//start after madatory arguments
	while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) { //keyword = region
        error->all(FLERR,"fix electrodeboundaries does not support regions");
    iarg += 2;
    } else error->all(FLERR,"Illegal fix electrodeboundaries command"); // not a recognized keyword
  }

  // zero out counters
  leftOx=0;
  leftOxAttempts=0;
  leftRed=0;
  leftRedAttempts=0;
  rightOx=0;
  rightOxAttempts=0;
  rightRed=0;
  rightRedAttempts=0;

  // random number generator, same for all procs
  random_equal = new RanPark(lmp,seed);

	atom->add_callback(0);

}

FixElectrodeBoundaries::~FixElectrodeBoundaries(){
  // destructor -- free us pointer arrays

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

  //check that active ion type exists
  int *type = atom->type;

    if (etype <= 0 || etype > atom->ntypes)
      error->all(FLERR,"Invalid atom type in fix electrodeboundaries command");

}

void FixElectrodeBoundaries::pre_exchange(){
  // for some number of attempts
  // pick an x(uniform),y(uniform),z(uniform to some cutoff) position
  // Check if ion is within dx --> If so attempt reduction
  // --> if not, attempt oxidation

  double coords[3];
  ylo = domain->boxlo[1];
  yhi = domain->boxhi[1];
  zlo = domain->boxlo[2];
  zhi = domain->boxhi[2];

  for (int i=0; i<ncycle; ++i){
    coord[0] = random_equal->uniform() * (xcut*2); //only want to sample near electrodes
    coord[1] = ylo + random_equal->uniform() * (yhi-ylo);
    coord[2] = zlo + random_equal->uniform() * (zhi-zlo);

    //translate coord[0] into x position
    if (coord[0] < xcut){ //left side
      coord[0] = coord[0] + xlo;

    }else{ //right side
      coord[0] = xhi - coord[0];
    }

    if (is_particle(coord)){
      //attempt reduction

    }else{
      //attempt oxidation
      
    }


  }
}

void FixElectrodeBoundaries::attempt_atomic_insertion_full(double *coord){

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

