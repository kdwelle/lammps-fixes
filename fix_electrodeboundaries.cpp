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
#include <stdlib.h>
#include <string.h>
#include <vector>

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "kspace.h"
#include "modify.h"
#include "neighbor.h"
#include "pair.h"
#include "random_park.h"
#include "update.h"

#include "fix_electrodeboundaries.h"
#include "fix.h"


using namespace std;
using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixElectrodeBoundaries::FixElectrodeBoundaries(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  idregion(NULL){

  dr = 2.0; //plus/minus search for ion in vicinity
  xcut = 0.5; //mean distance from electrode to check for electrochem
  ncycles = 1; //number of attempts per timestep
  pOxidation = 0.50; //probability of oxidation vs. reduction
  charge = 1.0;
  charge_flag = true;
  sigma = sqrt(force->boltz/force->mvv2e);
  intercalation = true; //default is to add/remove ions when reduce/oxidized
  neutralIndex = -1;
  porusLeft = false;
  porusRight = false;
  overpotential = 0;
  occupation = 0;

  if (narg < 8) error->all(FLERR,"Illegal fix electrodeboundaries command -- not enough arguments");

  // electrodes have to lie along the x-axis
  // Next argument is a distance between electrodes
  Vxlo = force->numeric(FLERR,arg[3]);
  dist = force->numeric(FLERR,arg[4]);
  Vxhi = Vxlo + dist;
  // Then voltage@defined plane and voltage difference between electrodes
  v0 = force->numeric(FLERR,arg[5]);
  dv = force->numeric(FLERR,arg[6]);
  etype = force->inumeric(FLERR,arg[7]); //type of atom that is electrochemically active
  seed = force->inumeric(FLERR,arg[8]);

  if (seed <= 0) error->all(FLERR,"Illegal fix electrodeboundaries command -- seed cannot be zero or negative");

  // fprintf(screen,"v0 is %f and dv is %f.\n", v0,dv);

  // optional arguments
  iregion = -1;
  idregion = NULL;

  int iarg = 9; //start after madatory arguments
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) { //keyword = region
      error->all(FLERR,"fix electrodeboundaries does not support regions");
    }
    else if (strcmp(arg[iarg],"ncycles") == 0) { //keyword = ncycles 
      ncycles = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    }else if (strcmp(arg[iarg],"couple") == 0) { //keyword = intercalation true/false neutralIndex
      intercalation = false; // if -1 then use no redox couple
      neutralIndex = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    }else if (strcmp(arg[iarg],"intercalation") == 0) { //keyword = intercalation true/false neutralIndex
      intercalation = true; // if -1 then use no redox couple
      neutralIndex = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    }else if (strcmp(arg[iarg],"porusLeft") == 0){ //keyword = porus ; triggers a porus electrode framework
      porusLeft = true;
      xend = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    }else if (strcmp(arg[iarg],"porusRight") == 0){ //keyword = porus ; triggers a porus electrode framework
      porusRight = true;
      xstart = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    }else if (strcmp(arg[iarg],"overpotential") == 0){ //keyword = overpotential ; adds a term to the de before prob. compute
      overpotential = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    }else if (strcmp(arg[iarg],"occupation") == 0){ //keyword = occupation ; includes a check against an occuptation lattice
      occupation = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
  	}else if (strcmp(arg[iarg],"pOxidation") == 0){ //keyword = occupation ; includes a check against an occuptation lattice
      pOxidation = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    }else if (strcmp(arg[iarg],"xcut") == 0){ //keyword = occupation ; includes a check against an occuptation lattice
      xcut = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    }else error->all(FLERR,"Illegal fix electrodeboundaries command"); // not a recognized keyword
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
  vector_flag = 1;
  size_vector = 8;

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

  // create a new group for interaction exclusions
  // used for attempted atom or molecule deletions
  if (intercalation){ //only needed if doind deletions
  if (!exclusion_group_bit) {
    char **group_arg = new char*[4];

    // create unique group name for atoms to be excluded
    int len = strlen(id) + 30;
    group_arg[0] = new char[len];
    sprintf(group_arg[0],"FixGCMC:gcmc_exclusion_group:%s",id);
    group_arg[1] = (char *) "subtract";
    group_arg[2] = (char *) "all";
    group_arg[3] = (char *) "all";
    group->assign(4,group_arg);
    exclusion_group = group->find(group_arg[0]);
    if (exclusion_group == -1)
      error->all(FLERR,"Could not find fix gcmc exclusion group ID");
    exclusion_group_bit = group->bitmask[exclusion_group];

    // neighbor list exclusion setup
    // turn off interactions between group all and the exclusion group
    int narg = 4;
    char **arg = new char*[narg];;
    arg[0] = (char *) "exclude";
    arg[1] = (char *) "group";
    arg[2] = group_arg[0];
    arg[3] = (char *) "all";
    neighbor->modify_params(narg,arg);
    delete [] group_arg[0];
    delete [] group_arg;
    delete [] arg;
  }
  }

  char *id_pe = (char *) "thermo_pe";
  int ipe = modify->find_compute(id_pe);
  c_pe = modify->compute[ipe];

  // Build occupation array if needed
  if (occupation){
    ylo = domain->boxlo[1];
    yhi = domain->boxhi[1];
    zlo = domain->boxlo[2];
    zhi = domain->boxhi[2];
    xhi = domain->boxhi[0];
    xlo = domain->boxlo[0];
    ny = round((yhi-ylo)*occupation);
    nx = round((xhi-xlo)*occupation);
    nz = round((zhi-zlo)*occupation);
    occ.resize(nx);
    for (int i=0; i<nx; ++i){
      occ[i].resize(ny);
      for (int j=0; j<ny; ++j){
        occ[i][j].resize(nz,false); //all start at false
      }
    }

  }
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

  int numRepeat;
  if (random_equal->uniform() < ncycles){
    numRepeat = ceil(ncycles);
    energy_stored = energy_full();
  }else{
    numRepeat = 0;
  }

  for (int i=0; i<numRepeat; ++i){
    int index;
    int side; // left = 0, right = 1;
    if (porusLeft == porusRight){ //both porus or neither porus
    	side = random_equal->uniform() > 0.5; // 50/50
    }else if (porusRight){
    	double dx = Vxhi - xstart; // length of porus regime
    	side = random_equal->uniform() < (dx/(dx+xcut)); //way more likely to ping right
	}else{ //porusleft is true
		double dx = xend - Vxlo;
		side = random_equal->uniform() > (dx/(dx+xcut)); //more likely for left
	}
	
    // left = 0, right = 1;
    coords[0] = get_x(side);
    coords[1] = ylo + random_equal->uniform() * (yhi-ylo);
    coords[2] = zlo + random_equal->uniform() * (zhi-zlo);

    //random chance for reduction vs oxidation
    if (random_equal->uniform() < pOxidation){
      //oxidation
      fprintf(screen, "oxidation %.2f %.2f %.2f ",coords[0],coords[1],coords[2]); 
      if (porus_side(side)){ //checks if this side is porus and if occupancy is turned on
        if (is_occupied(coords)){
          attempt_oxidation(coords, side);
        }else{ 
          fprintf(screen, "-1 -1 0 \n");
        }
      }else{
        attempt_oxidation(coords, side);
      }
    
    }else{
      //reduction
      fprintf(screen, "reduction ");
      if (porus_side(side)){
        if (!is_occupied(coords)){ //have to reduce to an empty site
          index = is_particle(coords,etype);
        }else{
          index = -1;
        }
      }else{ // no need to worry about occ
        index = is_particle(coords,etype);
      }
      if (index > 0){ //hack because image charges messes up when excluded atom is index 0
        //attempt reduction
        attempt_reduction(index, side);
      }else{
        fprintf(screen, "-1 -1 -1 -1 -1 0 \n");
      }
    }
  }

  // fprintf(screen, "Left oxidations: %d / %d \nRight oxidations: %d / %d \n", leftOx, leftOxAttempts, rightOx, rightOxAttempts);
  // fprintf(screen, "Left reductions: %d / %d \nRight reductions: %d / %d \n", leftRed, leftRedAttempts, rightRed, rightRedAttempts);
  
}

float FixElectrodeBoundaries::get_x(int side){
  // function that determines the x-coordinate of the oxidation/reduction
  if (side){ //right
    if (porusRight){
      return xstart + ((Vxhi - xstart)*random_equal->uniform());
    }else{
      float xtemp = -1*xcut*log(random_equal->uniform()); //exponential distribution centered at xcut
      return Vxhi - xtemp;
    }

  }else{ //left
    if (porusLeft){
      return Vxlo + ((xend-Vxlo)*random_equal->uniform());
    }else{
      float xtemp = -1*xcut*log(random_equal->uniform()); //exponential distribution centered at xcut
      return xtemp + Vxlo;
    }
  }

}

int FixElectrodeBoundaries::is_particle(double *coords, int typeI){
  // checks to see if there is a particle within dr of coords
  // returns the index of the atom or -1 if not found
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  int yperiodic = domain->yperiodic;
  int zperiodic = domain->zperiodic;
  double ylen = yhi-ylo;
  double zlen = zhi-zlo;

  float xmin=max(coords[0]-xcut, Vxlo);
  float xmax=min(coords[0]+xcut, Vxhi);
  float ymin=coords[1]-dr;
  float ymax=coords[1]+dr;
  float zmin=coords[2]-dr;
  float zmax=coords[2]+dr;

  //probably going to be very slow -- oh well let's try it
  for (int i=0; i<nlocal; ++i){
    if(mask[i] & groupbit){
      if (type[i]==typeI){
        if (x[i][0] > xmin){
          if (x[i][0] < xmax){
            if (x[i][1] > ymin || (yperiodic && ymin < ylo && x[i][1] > ymin+ylen)){
              if (x[i][1] < ymax || (yperiodic && ymax > yhi && x[i][1] < ymax-ylen)){
                if (x[i][2] > zmin || (zperiodic && zmin < zlo && x[i][2] > zmin+zlen)){
                  if (x[i][2] < zmax || (zperiodic && zmax > zhi && x[i][2] < zmax-zlen)){
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

bool FixElectrodeBoundaries::is_occupied(double* coords){
  int x = round((coords[0]-xlo)/occupation);
  int y = round((coords[1]-ylo)/occupation);
  int z = round((coords[2]-zlo)/occupation);
  if(x<0){x += nx; }else if(x>nx){x -= nx; }
  if(y<0){y += ny; }else if(y>nx){y -= ny; }
  if(z<0){z += nz; }else if(z>nz){z -= nz; }
  return occ[x][y][z];
}

void FixElectrodeBoundaries::set_occupation(double* coords, bool set){
  int x = round((coords[0]-xlo)/occupation);
  int y = round((coords[1]-ylo)/occupation);
  int z = round((coords[2]-zlo)/occupation);
  if(x<0){x += nx; }else if(x>nx){x -= nx; }
  if(y<0){y += ny; }else if(y>nx){y -= ny; }
  if(z<0){z += nz; }else if(z>nz){z -= nz; }
  occ[x][y][z] = set;
}

bool FixElectrodeBoundaries::porus_side(int side){
  //returns true if side is porus AND occ is turned on
  if (occupation){
    if (!side && porusLeft){
      return true;
    }else if (side && porusRight){
      return true;
    }
  }
  return false;
}

void FixElectrodeBoundaries::attempt_oxidation(double *coord, int side){

  side? rightOxAttempts++ : leftOxAttempts++;

  double energy_before = energy_stored;
  double original_energy = energy_stored;
  int m;
  double energy_after;
  double de;
  bool reject = false;
  bool found = false;
  int oldmask;
  double oldcoords[3];

  // Step 0: check to see if there is an occupied site to "draw" from
  
  // Step 1: Check LJ interactions using Boltzmann factors if intercaltion
  if (intercalation) {
    if (neutralIndex == -1){
      // add atom
      atom->avec->create_atom(etype,coord);
      m = atom->nlocal - 1; //This is -1 because it's after atom was created
      // add to groups
      // optionally add to type-based groups
      // don't add charge yet
      atom->mask[m] = groupbit;
      atom->v[m][0] = random_equal->gaussian()*sigma;
      atom->v[m][1] = random_equal->gaussian()*sigma;
      atom->v[m][2] = random_equal->gaussian()*sigma;
      modify->create_attribute(m); //what does this do? --> adds it to be modified by this fix

      atom->natoms++;
      if (atom->tag_enable) {
        atom->tag_extend();
        if (atom->map_style) atom->map_init();
      }
      atom->nghost = 0;

    } else {
      //find an unused neutral atom
      int *type = atom->type;
      int *mask = atom->mask;

      for (int i=0; i < atom->nlocal; ++i){
        if(mask[i] & groupbit){
          if (type[i]==neutralIndex){
            m=i;
            found = true;
            break;
          }
        }
      }
      if (found){
        atom->type[m] = etype;
        atom->q[m] = 0;
        // oldmask = atom->mask[m];
        // atom->mask[m] = groupbit;
        oldcoords[0] = atom->x[m][0];
        oldcoords[1] = atom->x[m][1];
        oldcoords[2] = atom->x[m][2];
        atom->x[m][0] = coord[0];
        atom->x[m][1] = coord[1];
        atom->x[m][2] = coord[2];

        // atom->v[m][0] = random_equal->gaussian()*sigma;
        // atom->v[m][1] = random_equal->gaussian()*sigma;
        // atom->v[m][2] = random_equal->gaussian()*sigma;
      } else{
        reject = true;
        fprintf(screen, "-2 -2 0 \n"); //out of neutral atoms
      }
    }
    energy_after = energy_full();
    de = energy_after-energy_before;

  } else { //not intercalation a//search for a redox couple
    m = is_particle(coord,neutralIndex);
    if (m == -1){
      reject = true;
      de = 10000;
      fprintf(screen, "-3 -3 0 \n"); // couldn't find a couple
    }else{
      de = -10;
    }
  }

  // TODO: edit so that uses correct kT for non-lj units!
  if(!reject){
    if(de < 0 || exp(-de) > random_equal->uniform()){ //accept, check electrostatics
      energy_before = energy_after; // Now want to compare charge and no charge
      if (charge_flag) atom->q[m] = charge;
      energy_after = energy_full();
      de = energy_after-energy_before;
      double prob = get_transfer_probability(de,side,1);
      fprintf(screen, "%f %f ",de, prob);
      if (random_equal->uniform() < prob ){
        energy_stored = energy_after;
        (side)? rightOx++ : leftOx++;
        fprintf(screen, "1 \n"); // accepted
        if(porus_side(side)){
          set_occupation(coord,0);
        }
      }else{  // charge transfer move rejected
        reject = true;
        fprintf(screen, "0 \n");
      }
    }else{ // insertion move rejected
      reject = true;
      fprintf(screen, "-1 -1 0 \n"); // insertion move rejected
    }
  }

  if (reject){
    if (neutralIndex == -1) {
      int nlocal = atom->nlocal;
      // fprintf(screen, "%s %d %s %d %s", "not accepted, m is ", m, " nlocal is ",nlocal, "\n");
      
      while (m < atom->nlocal-1){
        atom->natoms--;
        atom->nlocal--;
      }
      //delete atom
      atom->natoms--;
      atom->nlocal--;
    
    } else if (found){
      if (charge_flag) atom->q[m] = 0.0;
      // atom->mask[m] = oldmask;
      atom->type[m] = neutralIndex;
      atom->x[m][0] = oldcoords[0];
      atom->x[m][1] = oldcoords[1];
      atom->x[m][2] = oldcoords[2];
    }

    if (modify->n_pre_force) modify->pre_force(0);
    if (force->kspace) force->kspace->qsum_qsq();
    energy_stored = original_energy;
  }
}

void FixElectrodeBoundaries::attempt_reduction(int i, int side){
  double q_tmp=0;
  bool ctAccepted = false;

  double* coord=atom->x[i];
  fprintf(screen, "%f %f %f", coord[0],coord[1],coord[2]);

  side? rightRedAttempts++ : leftRedAttempts++;
  double energy_before = energy_stored;
  double original_energy = energy_before;
  double energy_after;
  double de;
  

  // First, check electrostatics
  if(charge_flag){
    q_tmp = atom->q[i];
    atom->q[i] = 0.0;
    energy_after = energy_full();
    de = energy_after-energy_before;
    double prob = get_transfer_probability(de,side,0);

    fprintf(screen, "%f %f ", de, prob);

    if (random_equal->uniform() < prob) {  // check to see if can remove atom
      ctAccepted = true;
      energy_before = energy_after;  // no charge is new baseline energy
    } else{  // charge transfer move rejected
      ctAccepted = false;
      fprintf(screen, "0 \n");
      atom->q[i] = q_tmp;
      energy_stored = original_energy;
    }
  } else {  // no charge, ct automatically accepted
    ctAccepted = true;
    fprintf(screen, "no charge on this reduction ion");
  }

  if (ctAccepted){  // check to see if noncharged interactions okay
    if (intercalation){ //only need to check for short-range interactions if intercalation
      int tmpmask;
      if (i >= 0) {  // exclude atom
        tmpmask = atom->mask[i];
        atom->mask[i] = exclusion_group_bit;
      }
      energy_after = energy_full();
      de = energy_after-energy_before;

      // TODO: edit so that uses correct kT for non-lj units!
      if(de < 0 || exp(-de) > random_equal->uniform()){ //accept boltzmann and move
        atom->mask[i] = tmpmask;
        remove_atom(i);
        if (atom->map_style) atom->map_init();  //what does this do?
        side? rightRed++ : leftRed++;
        energy_stored = energy_after;
        fprintf(screen, "1 \n");
        if (occupation){
          set_occupation(coord,1);
        }

      } else { //not accepted
        fprintf(screen, "0 \n");
        // reset everything
        atom->mask[i] = tmpmask;
        if (charge_flag) atom->q[i] = q_tmp;
        if (modify->n_pre_force) modify->pre_force(0);
        if (force->kspace) force->kspace->qsum_qsq();
        energy_stored = original_energy;
      }
    }else{
      side? rightRed++ : leftRed++;
      remove_atom(i); //no check needed for redox couple, we already checked charge
      if (occupation){
        set_occupation(coord,1);
      }
    }
  }
}

void FixElectrodeBoundaries::remove_atom(int i){
  if (neutralIndex == -1){ //no neutral atoms or hidden atom group, just delete it
    atom->avec->copy(atom->nlocal-1,i,1);
    atom->nlocal--;
    atom->natoms--;
    if (modify->n_pre_force) modify->pre_force(0);
    if (force->kspace) force->kspace->qsum_qsq();
  }else{
    atom->type[i] = neutralIndex; //else turn it into the neutral species
    if (modify->n_pre_force) modify->pre_force(0);
    if (force->kspace) force->kspace->qsum_qsq(); 
  }
}

float FixElectrodeBoundaries::get_transfer_probability(float dE, int side, int redox){
  // get P(x) from madeleung potential
  // side is 0 for left and 1 for right
  // redox is 0 for reduction and 1 for oxidation

  //TODO: include temperature in this
  float x;
  float fermi = v0 + (side*dv);
  if (redox){
    x = dE - fermi; // oxidation
  }else{
    x = dE + fermi; //fermi - (-dE) = 1-P(ox)
  }
  x = x+overpotential;
  return 1/(1+exp(x));
}
  


/* ----------------------------------------------------------------------
   compute total system energy incuding fixes
------------------------------------------------------------------------- */

double FixElectrodeBoundaries::energy_full()
{
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (modify->n_pre_neighbor) modify->pre_neighbor();
  neighbor->build(1);
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

/* ----------------------------------------------------------------------
   return stats about oxidized and reduced atoms
------------------------------------------------------------------------- */

double FixElectrodeBoundaries::compute_vector(int n){
  //n=  0 --> Left Ox
  //    1 --> Left Red
  //    2 --> Right Ox
  //    3 --> Right Red
  //    4 --> Left Ox Attemps
  //    5 --> Left Red Attemps
  //    6 --> Right Ox Attemps
  //    7 --> Right Red Attemps
  
  switch (n){
    case 0: return leftOx;
    case 1: return leftRed;
    case 2: return rightOx;
    case 3: return rightRed;
    case 4: return leftOxAttempts;
    case 5: return leftRedAttempts;
    case 6: return rightOxAttempts;
    case 7: return rightRedAttempts;
  }
  return -1; //something went wrong here
}

