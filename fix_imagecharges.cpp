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

#include <string.h>
#include <stdlib.h>
#include "fix_imagecharges.h"
#include "fix_electrodeboundaries.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "region.h"
#include "variable.h"



using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixImageCharges::FixImageCharges(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  pxstr(NULL), pystr(NULL), pzstr(NULL),
  nxstr(NULL), nystr(NULL), nzstr(NULL),
  idregion(NULL), scalestr(NULL), imagei(NULL), imageid(NULL){
    if (narg < 10) error->all(FLERR,"Illegal fix imagecharges command -- not enough arguments");

    // initialize the array to keep track of image charge associations
    memory->create(imagei, atom->nmax, "FixImageCharges::imagei");
    memory->create(imageid, atom->nmax, "FixImageCharges::imagei");

    // first three arguments define a point on the plane
    if (strstr(arg[3],"v_") == arg[3]) { // if starts with the string v_
      int n = strlen(&arg[3][2]) + 1;    // length of variable name + 1 (calls strlen on a pointer to the third character in the string)
      pxstr = new char[n];                // allocate enough space for variable
      strcpy(pxstr,&arg[3][2]);           // copy varible name to class
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

  // itype -- index for the image charge types
  // go ahead and evaluate itype only once
  if (strstr(arg[9],"v_") == arg[9]) {
    int itypevar = input->variable->find(&arg[8][2]);
    if (itypevar < 0){
      error->all(FLERR,"Variable name for fix imagecharges (itype) does not exist");
    }
    if (input->variable->equalstyle(itypevar)){
      itype = input->variable->compute_equal(itypevar); //doesn't check before converting to int
    } else {
      error->all(FLERR,"Variable for fix imagecharges (itype) is invalid style");
    }
  } else {
    itype = force->inumeric(FLERR,arg[9]);
  }


  // optional arguments
  iregion = -1;
  idregion = NULL;
  scale = 1.0;

  int iarg = 10;  //start after madatory arguments
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

    } else if (strcmp(arg[iarg],"scale") == 0) { //keyword = scale
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix imagecharges command"); //check to make sure there is a next argument
      if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) { //variable
        int n = strlen(&arg[iarg+1][2]) + 1;
        scalestr = new char[n];
        strcpy(scalestr,&arg[iarg+1][2]);
        } else {
          scale = force->numeric(FLERR,arg[iarg+1]);
          scalestyle = CONSTANT;
        }
        iarg += 2;

    } else error->all(FLERR,"Illegal fix imagecharges command"); // not a recognized keyword
  }

  //This fix produces per-atom information
  peratom_flag = 1;
  peratom_freq = 1;
  vector_atom = imageid;

  atom->add_callback(0);
}


FixImageCharges::~FixImageCharges(){
  // destructor -- free up any arrays of pointers
  delete [] pxstr;
  delete [] pystr;
  delete [] pzstr;
  delete [] nxstr;
  delete [] nystr;
  delete [] nzstr;
  delete [] scalestr;
  delete [] idregion;
  memory->destroy(imagei);
  memory->destroy(imageid);
  atom->delete_callback(id,0);
}

int FixImageCharges::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  return mask;
}

void FixImageCharges::init(){

  // check variables
  if (pxstr) {
    pxvar = input->variable->find(pxstr);
    if (pxvar < 0)
      error->all(FLERR,"Variable name for fix imagecharges (px) does not exist");
    if (input->variable->equalstyle(pxvar)) pxstyle = EQUAL;
    else if (input->variable->atomstyle(pxvar)) pxstyle = ATOM;
    else error->all(FLERR,"Variable for fix imagecharges (px) is invalid style");
  }
  if (pystr) {
    pyvar = input->variable->find(pystr);
    if (pyvar < 0)
      error->all(FLERR,"Variable name for fix imagecharges (py) does not exist");
    if (input->variable->equalstyle(pyvar)) pystyle = EQUAL;
    else if (input->variable->atomstyle(pyvar)) pystyle = ATOM;
    else error->all(FLERR,"Variable for fix imagecharges (py) is invalid style");
  }
  if (pzstr) {
    pzvar = input->variable->find(pzstr);
    if (pzvar < 0)
      error->all(FLERR,"Variable name for fix imagecharges (pz) does not exist");
    if (input->variable->equalstyle(pzvar)) pzstyle = EQUAL;
    else if (input->variable->atomstyle(pzvar)) pzstyle = ATOM;
    else error->all(FLERR,"Variable for fix imagecharges (pz) is invalid style");
  }
  if (nxstr) {
    nxvar = input->variable->find(nxstr);
    if (nxvar < 0)
      error->all(FLERR,"Variable name for fix imagecharges does not exist");
    if (input->variable->equalstyle(nxvar)) nxstyle = EQUAL;
    else if (input->variable->atomstyle(nxvar)) nxstyle = ATOM;
    else error->all(FLERR,"Variable for fix imagecharges (nx) is invalid style");
  }
  if (nystr) {
    nyvar = input->variable->find(nystr);
    if (nyvar < 0)
      error->all(FLERR,"Variable name for fix imagecharges (ny) does not exist");
    if (input->variable->equalstyle(nyvar)) nystyle = EQUAL;
    else if (input->variable->atomstyle(nyvar)) nystyle = ATOM;
    else error->all(FLERR,"Variable for fix imagecharges (ny) is invalid style");
  }
  if (nzstr) {
    nzvar = input->variable->find(nzstr);
    if (nzvar < 0)
      error->all(FLERR,"Variable name for fix imagecharges (nz) does not exist");
    if (input->variable->equalstyle(nzvar)) nzstyle = EQUAL;
    else if (input->variable->atomstyle(nzvar)) nzstyle = ATOM;
    else error->all(FLERR,"Variable for fix imagecharges (nz) is invalid style");
  }


  // set index and check validity of region
  if (iregion >= 0) {
    iregion = domain->find_region(idregion);
    if (iregion == -1)
      error->all(FLERR,"Region ID for fix imagecharges does not exist");
  }

  // find scale variable if set
  if (scalestr) {
    scalevar = input->variable->find(scalestr);
    if (scalevar < 0)
      error->all(FLERR,"Variable name for fix setforce (scale) does not exist");
    if (input->variable->equalstyle(scalevar)) scalestyle = EQUAL;
    else if (input->variable->atomstyle(scalevar)) scalestyle = ATOM;
    else error->all(FLERR,"Variable for fix setforce (scale) is invalid style");
  }

  if (pxstyle == ATOM || pystyle == ATOM || pzstyle == ATOM || nxstyle == ATOM || nystyle == ATOM || nzstyle == ATOM || scalestyle == ATOM)
    varflag = ATOM;
  else if (pxstyle == EQUAL || pystyle == EQUAL || pzstyle == EQUAL || nxstyle == EQUAL || nystyle == EQUAL || nzstyle == EQUAL || scalestyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  // get the index of the "all" groupbit
  // char *tempFind = new char[3];
  // std::sprintf(tempFind,"%s","all");
  // int temp_ind = group->find(tempFind);
  // group_bit_all = group->bitmask[temp_ind];
  // std::fprintf(screen, "%s %d %s", "group_bit_all is: ", group_bit_all, "\n");

}

void FixImageCharges::setup_pre_force(int vflag){
  // add the image charges as new atoms
  // todo: will probably have to link images to atoms at some point

  double **x = atom->x;      //positions
  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;
  int nadded = 0;
  int atomIndex = nlocal;


  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])){
          imagei[i] = -2;
          imageid[i] = -2;
          continue;
        }
        // transform coordinates across plane
        double nnorm = sqrt(nxvalue*nxvalue + nyvalue*nyvalue + nzvalue*nzvalue);
        double prefactor = 2*(nxvalue/nnorm*x[i][0] + nyvalue/nnorm*x[i][1] + nzvalue/nnorm*x[i][2]);
        double delta = 2*(nxvalue/nnorm*pxvalue + nyvalue/nnorm*pyvalue + nzvalue/nnorm*pzvalue);
        double r[3];
        r[0] = x[i][0] - (prefactor-delta)*nxvalue;
        r[1] = x[i][1] - (prefactor-delta)*nyvalue;
        r[2] = x[i][2] - (prefactor-delta)*nzvalue;
        //add a new atom
        nadded++;
        atom->avec->create_atom(itype,r);
        atom->q[atomIndex] = -1*scale*q[i];
        atom->mask[atomIndex] = groupbit;
        imagei[i] = atomIndex;
        imageid[i] = atomIndex;
        imagei[atomIndex] = -1;
        imageid[atomIndex] = -1;
        atomIndex++;
      }
      vector_atom = imageid;

  }
  //MPI adding here TODO
  if(nadded){
      atom->natoms += nadded;
  if (atom->natoms < 0)
    error->all(FLERR,"Too many total atoms");
  if (atom->tag_enable) atom->tag_extend();
  if (atom->map_style) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
    }
  }

  // //same for exclusion
  // int len = strlen(id) + 30;
  // char *group_arg = new char[len];
  // std::sprintf(group_arg,"FixGCMC:gcmc_exclusion_group:%s",id);
  // exclusion_group = group->find(group_arg);
  // std::fprintf(screen, "%s %d %s", "exclusion_bit is: ", exclusion_group, "\n");

}

void FixImageCharges::pre_force(int vflag){
  // Move all image charges to proper locations before calculating forces
  // Check to make sure all atoms in region/group have images and no extra images
  // imagei = -1 if is an image charge, -2 if not currently in region, undefined for not in group
  double **x = atom->x;      //positions
  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;
  int nadded = 0;
  int atomIndex = nlocal;
  int seenCount = 0;
  int reqCount = 0;
  int tmpmask;

  bool toDelete = false;
  int dlist[nlocal]; //list to use as a mask for atoms that need to be deleted
  int imList[nlocal];
  int imReqList[nlocal];

  for (int i=0; i<nlocal; i++){ //initialize
    dlist[i] = 0;
    imList[i] = 0;
    imReqList[i] = 0;
  }

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }

  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++){
      if (mask[i] & groupbit) {
        // if (region && !region->match(x[i][0],x[i][1],x[i][2])){
        //   // check to see if there's an existing image charge to be deleted
        //   if (imagei[i] >= 0) {
        //       toDelete = true;
        //       dlist[imagei[i]] = 1;
        //       nadded--; //got rid of an atom so keep accounting constant
        //       imagei[i] = -2;
        //       imageid[i] = -2;
        //   }
        //   continue;
        // }
        // check to see if an image charge already exists
        int j = imagei[i];

        if(j == -1){ // this is an image charge, will get taken care of later or deleted
          dlist[i] = !dlist[i];
          seenCount++;
        }else{
          // get new position -- transform coordinates across plane
          double nnorm = sqrt(nxvalue*nxvalue + nyvalue*nyvalue + nzvalue*nzvalue);
          double prefactor = 2*(nxvalue/nnorm*x[i][0] + nyvalue/nnorm*x[i][1] + nzvalue/nnorm*x[i][2]);
          double delta = 2*(nxvalue/nnorm*pxvalue + nyvalue/nnorm*pyvalue + nzvalue/nnorm*pzvalue);
          double r[3];
          r[0] = x[i][0] - (prefactor-delta)*nxvalue;
          r[1] = x[i][1] - (prefactor-delta)*nyvalue;
          r[2] = x[i][2] - (prefactor-delta)*nzvalue;

        
          if(j == -2 || j==0 || j >= nlocal){ //used to not be in region or is new atom
            // probably won't fail even if j was supposed to be zero
            j=atomIndex;

            fprintf(screen,"%s %d %s %d %s", "New atom ", i, " gets image ", j, "\n");
            atomIndex++;
            nadded++;
            atom->avec->create_atom(itype,r); //add a new atom
            atom->mask[atomIndex] = groupbit;
          
          }else{
            // mark that we updated/saw this image
            dlist[j] = !dlist[j];
            reqCount++;
            
            // update image coordinates
            for (int k=0; k<3; ++k){
              x[j][k] = r[k];
            }
          }
          atom->q[j] = -1*scale*q[i]; //update charge
          imagei[i] = j;
          imageid[i] = j;
          imagei[j] = -1;
          imageid[j] = -1;
        }
      }
    }
    // } else { TODO add the atom and equal style interpretations
  }

  // deal with the deleteList
  if (seenCount > reqCount){
    fprintf(screen,"%s", "in imList: ");
    for (int i=0; i<nlocal; i++){
      if (dlist[i]) fprintf(screen,"%d %s", i, ", ");
    }
    toDelete = true;
    fprintf(screen,"%s", "\n");
  }

  if (toDelete){
    int i = 0;
    while (i<nlocal){
      if (dlist[i]){
        atom -> avec -> copy(nlocal-1, i, 1);
        dlist[i] = dlist[nlocal-1];
        atom->nlocal--;
        nadded--;
      } else i++;
    }
    
  }

  if(nadded){
    atom->natoms += nadded;
    if (atom->natoms < 0)
      error->all(FLERR,"Too many total atoms");
    if (atom->tag_enable) atom->tag_extend();
    if (atom->map_style) {
      atom->nghost = 0;
      atom->map_init();
      atom->map_set();
    }
  }
}

void FixImageCharges::post_force(int vflag){
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // update region if necessary

  Region *region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
  }


  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++){
      // if (mask[i] & groupbit) {
        // if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        // check whether this is an image charge
        if (imagei[i] == -1){
          f[i][0] = 0; //if so zero forces
          f[i][1] = 0;
          f[i][2] = 0;
        }
      // }
    }
  }

}

double FixImageCharges::memory_usage(){
  int nmax = atom->nmax;
  double bytes = 0.0;
  bytes += nmax * sizeof(int);
  return bytes;
}

void FixImageCharges::grow_arrays(int nmax){
  memory->grow(this->imagei,nmax,"FixImageCharges::imagei");
  memory->grow(this->imageid,nmax,"FixImageCharges::imageid");
  vector_atom = imageid;
}

void FixImageCharges::copy_arrays(int i, int j, int delflag){
  int i1 = imagei[i];
  imagei[j] = imagei[i];
  imageid[j] = imageid[i];

  if (i1 == -1){ //this is an image charge
    for (int x=0; x < atom->nlocal+1; ++x){ //need to loop over empty space at end used for sorting
      if (imagei[x] == i){
        imagei[x] = j; //now points to new location
        imageid[x] = j;
      }
    }
  }

  // imageid[j] = imageid[i];
  // memcpy(&imagei[j], &imagei[i], sizeof(int));
  // memcpy(&imageid[j], &imageid[i], sizeof(double));
}

void FixImageCharges::set_arrays(int i){
  memset(&imagei[i], 0, sizeof(int));
  memset(&imageid[i], 0, sizeof(int));
}
