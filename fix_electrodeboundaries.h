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
/* ----------------------------------------------------------------------
   Contributing author: Kaitlyn Dwelle (MIT)
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
  double compute_vector(int);

  // double memory_usage();
  // void grow_arrays(int);
  // void copy_arrays(int, int);
  // void set_arrays(int);


 protected: 
  double xlo,xhi,dist,v0,dv,sigma,charge;
  double ylo,yhi,zlo,zhi;
  bool charge_flag, intercalation;
  double energy_stored,ncycles,pOxidation;
  int exclusion_group, neutralIndex;
  int varflag,iregion,etype,seed;
  int leftOx, leftOxAttempts;
  int leftRed, leftRedAttempts;
  int rightOx, rightOxAttempts;
  int rightRed, rightRedAttempts;
  char *idregion;

  double dr,xcut;

  int force_flag;
  int nlevels_respa,ilevel_respa;
  int exclusion_group_bit;


  int maxatom;
  double **sforce;

  class RanPark *random_equal;
  class Compute *c_pe;

  //internal subroutines
  int is_particle(double*,int);
  void attempt_oxidation(double*, int);
  void attempt_reduction(double*, int);
  void remove_atom(int);

  float get_transfer_probability(float, int, int);
  double energy_full();



};
}

#endif
#endif