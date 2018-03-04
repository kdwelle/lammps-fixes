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
  charge = 1;

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

  if (etype <= 0 || etype > atom->ntypes){
    error->all(FLERR,"Invalid atom type in fix electrodeboundaries command");
  }
  // set energy
  energy_stored = energy_full();

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
    int side
    if (coord[0] < xcut){ //left side
      coord[0] = coord[0] + xlo;
      side = 0;

    }else{ //right side
      coord[0] = xhi - coord[0];
      side =1;
    }

    int index = is_particle(coord);
    if (index == -1){
      //attempt reduction
      attempt_reduction(index, side);

    }else{
      //attempt oxidation
      attempt_oxidation(coord, side);

    }

  }
}

int is_particle(double *coords){
  // checks to see if there is a particle within dr of coords
  // returns the index of the atom or -1 if not found
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  xmin=coords[0]-dr;
  xmax=coorder[0]+dr;
  ymin=coords[1]-dr;
  ymax=coords[1]+dr;
  zmin=coords[2]-dr;
  zmax=coords[2]+dr;

  //probably going to be very slow -- oh well let's try it
  for (int i=0; i<nlocal; ++i){
    if(mask[i] & groupbit){
      if (tpye[i]==etype){   
        if (x[i][0] > xmin){
          if (x[i][0] < xmax){
            if (x[i][1] > ymin){
              if (x[i][1] < ymax){
                if (x[i][2] > zmin){
                  if (x[i][2] < zmax){
                    return i;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return -1;
}

void FixElectrodeBoundaries::attempt_oxidation(double *coord, side){
  side? rightOxAttempts++ : leftOxAttempts++;

  double energy_before = energy_stored;

  // add atom
  atom->avec->create_atom(etype,coord);
  int m = atom->nlocal - 1; //This is -1 because it's after atom was created

  // add to groups
  // optionally add to type-based groups

  atom->mask[m] = groupbitall;
  for (int igroup = 0; igroup < ngrouptypes; igroup++) {
    if (type == grouptypes[igroup])
    atom->mask[m] |= grouptypebits[igroup];
  }

  atom->v[m][0] = random_unequal->gaussian()*sigma;
  atom->v[m][1] = random_unequal->gaussian()*sigma;
  atom->v[m][2] = random_unequal->gaussian()*sigma;
  if (charge_flag) atom->q[m] = charge;
  modify->create_attribute(m); //what does this do?

  atom->natoms++;
  if (atom->tag_enable) {
    atom->tag_extend();
    if (atom->map_style) atom->map_init();
  }
  atom->nghost = 0;

  if (force->kspace) force->kspace->qsum_qsq();
  double energy_after = energy_full();

  if (random_equal->uniform() > get_transfer_probability(energy_after-energy_before,side) ){ 
  // metropolis condition -- greater than becaude get_transfer probability return p(x) for reduction, oxidation = 1-P(x)
    energy_stored = energy_after;
    (side)? rightOx++ : leftOx++;
  }else{ //not accepted
    atom->natoms--;
    if (proc_flag) atom->nlocal--;
    if (force->kspace) force->kspace->qsum_qsq();
    energy_stored = energy_before;
  }

}

void FixElectrodeBoundaries::attempt_reduction(int i, int side){
  double q_tmp;
  const int q_flag = atom->q_flag;

  side? rightRedAttempts++ : leftRedAttempts++;
  double energy_before = energy_stored;

  int tmpmask;
  if (i >= 0) {
    tmpmask = atom->mask[i];
    atom->mask[i] = exclusion_group_bit;
    if (q_flag) {
      q_tmp = atom->q[i];
      atom->q[i] = 0.0;
    }
  }
  if (force->kspace) force->kspace->qsum_qsq();
  double energy_after = energy_full();

  if (random_equal->uniform() < get_transfer_probability(energy_after-energy_before,side)) {
    atom->avec->copy(atom->nlocal-1,i,1);
    atom->nlocal--;
    atom->natoms--;
    if (atom->map_style) atom->map_init();
    side? rightRed++ : leftRed++;
    energy_stored = energy_after;
  } else { //not accepted
    atom->mask[i] = tmpmask;
    if (q_flag) atom->q[i] = q_tmp;
    if (force->kspace) force->kspace->qsum_qsq();
    energy_stored = energy_before;
  }

}

float FixElectrodeBoundaries::get_transfer_probability(float dE, int side){
  //get P(x) from madeleung potential
  // Note for positive ion, madelueng potential = dE because q=+1
  //this will return oxidation potential, reduction is 1-P(x)

  //TODO: include temperature in this

  float x = dE + v0 + side*(dv); 
  return 1/(1+exp(x));
}




/* ----------------------------------------------------------------------
   compute total system energy incuding fixes
------------------------------------------------------------------------- */

double ElectrodeBoundaries::energy_full()
{
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (modify->n_pre_neighbor) modify->pre_neighbor();
  neighbor->build();
  int eflag = 1;
  int vflag = 0;

  // clear forces so they don't accumulate over multiple
  // calls within fix gcmc timestep, e.g. for fix shake
  
  size_t nbytes = sizeof(double) * (atom->nlocal + atom->nghost);
  if (nbytes) memset(&atom->f[0][0],0,3*nbytes);

  if (modify->n_pre_force) modify->pre_force(vflag);
  if (force->pair) force->pair->compute(eflag,vflag);
  if (force->kspace) force->kspace->compute(eflag,vflag);

  // unlike Verlet, not performing a reverse_comm() or forces here
  // b/c GCMC does not care about forces
  // don't think it will mess up energy due to any post_force() fixes

  if (modify->n_post_force) modify->post_force(vflag);
  if (modify->n_end_of_step) modify->end_of_step();

  // NOTE: all fixes with THERMO_ENERGY mask set and which
  //   operate at pre_force() or post_force() or end_of_step()
  //   and which user has enable via fix_modify thermo yes,
  //   will contribute to total MC energy via pe->compute_scalar()

  update->eflag_global = update->ntimestep;
  double total_energy = c_pe->compute_scalar();

  return total_energy;
}


