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

FixStyle(electrodeboundaries,FixElectrodeBoundaries)

#else

#ifndef LMP_FIX_ELECTRODEBOUNDARIES_H
#define LMP_FIX_ELECTRODEBOUNDARIES_H

#include "fix.h"

namespace LAMMPS_NS {

class FixElectrodeBoundaries : public Fix {
 public:
  FixElectrodeBoundaries(class LAMMPS *, int, char **);
  ~FixElectrodeBoundaries();
  int setmask();
  void init();
  void pre_exchange();

  //internal subroutines
  void attempt_charge_transfer();
  int is_particle(double*);
  void attempt_atomic_insertion_full(double*)

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  void set_arrays(int);


 protected: 
  double xlo, xhi, dist,v0,dv;
  double ylo,yhi,zlo,zhi;
  int varflag,iregion,etype,ncycles;
  int leftOx, leftOxAttempts
  int leftRed, leftRedAttempts;
  int rightOx, rightOxAttempts;
  int rightRed, rightRedAttempts;
  char *idregion;

  double dr,xcut;

  int force_flag;
  int nlevels_respa,ilevel_respa;

  int maxatom;
  double **sforce;

  class RanPark *random_equal;


};
}

#endif
#endif