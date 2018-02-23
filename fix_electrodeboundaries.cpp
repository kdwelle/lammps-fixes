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

  if (narg < 12) error->all(FLERR,"Illegal fix electrodeboundaries command -- not enough arguments");

	// first three arguments define a point on the plane
	if (strstr(arg[3],"v_") == arg[3]) {   // if starts with the string v_
    int n = strlen(&arg[3][2]) + 1;      // length of variable name + 1 (calls strlen on a pointer to the third character in the string)
    pxstr = new char[n];                 // allocate enough space for variable
    strcpy(pxstr,&arg[3][2]);            // copy varible name to class
	} else {
    pxvalue = force->numeric(FLERR,arg[3]); // else just convert input to a number and store
    pxstyle = CONSTANT;
	}

	if (strstr(arg[4],"v_") == arg[4]) { 
    int n = strlen(&arg[4][2]) + 1;    
    pystr = new char[n];                
    strcpy(pystr,&arg[4][2]);           
  } else {
    pyvalue = force->numeric(FLERR,arg[4]); 
    pystyle = CONSTANT;
  }
	
	if (strstr(arg[5],"v_") == arg[5]) { 
    int n = strlen(&arg[5][2]) + 1;    
    pzstr = new char[n];                
    strcpy(pzstr,&arg[5][2]);           
  } else {
    pzvalue = force->numeric(FLERR,arg[5]); 
    pzstyle = CONSTANT;
  }

  // next three arguments define a vector normal to the plane
  if (strstr(arg[6],"v_") == arg[6]) { 
    int n = strlen(&arg[6][2]) + 1;    
    nxstr = new char[n];                
    strcpy(nxstr,&arg[6][2]);           
  } else {
    nxvalue = force->numeric(FLERR,arg[6]); 
    nxstyle = CONSTANT;
  }

  if (strstr(arg[7],"v_") == arg[7]) { 
    int n = strlen(&arg[7][2]) + 1;    
    nystr = new char[n];                
    strcpy(nystr,&arg[7][2]);           
  } else {
    nyvalue = force->numeric(FLERR,arg[7]); 
    nystyle = CONSTANT;
  }

  if (strstr(arg[8],"v_") == arg[8]) { 
    int n = strlen(&arg[8][2]) + 1;
    nzstr = new char[n];
    strcpy(nzstr,&arg[8][2]);
  } else {
    nzvalue = force->numeric(FLERR,arg[8]); 
    nzstyle = CONSTANT;
  }

  // Next argument is a distance between electrodes 
  if (strstr(arg[9],"v_") == arg[9]) { 
    int n = strlen(&arg[9][2]) + 1;
    nzstr = new char[n];
    strcpy(diststr,&arg[9][2]);
  } else {
    distvalue = force->numeric(FLERR,arg[9]); 
    diststyle = CONSTANT;
  }
  // Then left voltage and voltage difference
  if (strstr(arg[10],"v_") == arg[10]) { 
    int n = strlen(&arg[10][2]) + 1;
    vstr = new char[n];
    strcpy(diststr,&arg[10][2]);
  } else {
    vvalue = force->numeric(FLERR,arg[10]); 
    vstyle = CONSTANT;
  }
  if (strstr(arg[11],"v_") == arg[11]) { 
    int n = strlen(&arg[11][2]) + 1;
    dvstr = new char[n];
    strcpy(diststr,&arg[11][2]);
  } else {
    dvvalue = force->numeric(FLERR,arg[11]); 
    dvstyle = CONSTANT;
  }

	// optional arguments
	iregion = -1;
	idregion = NULL;
	scale = 1.0;

	int iarg = 12;	//start after madatory arguments
	while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) { //keyword = region
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix imagecharges command"); //check to make sure there is a next argument
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1)
        error->all(FLERR,"Region ID for fix imagescharges does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    
    } else error->all(FLERR,"Illegal fix imagecharges command"); // not a recognized keyword
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
  // for some number of atoms
  // pick an ion
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

