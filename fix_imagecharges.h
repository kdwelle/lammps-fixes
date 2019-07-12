/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(imagecharges,FixImageCharges)

#else

#ifndef LMP_FIX_IMAGE_CHARGES_H
#define LMP_FIX_IMAGE_CHARGES_H

#include "fix.h"

namespace LAMMPS_NS {

class FixImageCharges : public Fix {
 public:
  FixImageCharges(class LAMMPS *, int, char **);
  virtual ~FixImageCharges();
  int setmask();
  virtual void init();
  void min_setup_pre_force(int);
  void setup_pre_force(int);
  void min_pre_force(int);
  void pre_force(int);
  void min_post_force(int); 
  void post_force(int);
  void post_run();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int,int);
  void set_arrays(int);


 protected: 
  double pxvalue,pyvalue,pzvalue,nxvalue,nyvalue,nzvalue,scale;
  double energy_stored;
  int varflag,iregion,itype;
  int exclusionAtom; 
  char *pxstr,*pystr,*pzstr,*nxstr,*nystr,*nzstr;
  char *idregion, *scalestr;
  int pxvar,pyvar,pzvar,nxvar,nyvar,nzvar,scalevar; 
  int pxstyle,pystyle,pzstyle,nxstyle,nystyle,nzstyle,scalestyle;
  int *imagei;
  double *imageid;


  double foriginal[3],foriginal_all[3];
  int force_flag;
  int nlevels_respa,ilevel_respa; 

  int maxatom;
  double **sforce;
};

}

#endif
#endif

/* ERROR/WARNING messages:



*/
