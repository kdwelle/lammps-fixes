
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
#include <vector>

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
  double Vxlo,Vxhi; //left electrode coordinate, right electrode coordinate
  double dist;      //distace between electrodes
  double v0,dv;     //left electrode voltage, volage difference between electrodes
  double sigma,charge; //width for random distributions for atom velocity, chage of active ion
  double ylo,yhi,zlo,zhi,xlo,xhi; //total box boundaries
  bool charge_flag, intercalation, porusLeft,porusRight;
  double energy_stored,ncycles,pOxidation, xstart, xend;
  double occupation,v0Increment;
  int exclusion_group, neutralIndex;
  int varflag,iregion,etype,seed,overpotential;
  int leftOx, leftOxAttempts;
  int leftRed, leftRedAttempts;
  int rightOx, rightOxAttempts;
  int rightRed, rightRedAttempts;
  char *idregion;

  double dr,xcut;

  int force_flag;
  int nlevels_respa,ilevel_respa;
  int exclusion_group_bit;

  std::vector< std::vector< std::vector<bool> > >  occ; //vector that holds info about site occupation
  int nx,ny,nz; //lengths of the occ vector
  int maxatom;
  double **sforce;


  class RanPark *random_equal;
  class Compute *c_pe;

  //internal subroutines
  int is_particle(double*,int);
  void attempt_oxidation(double*, int);
  void attempt_reduction(int, int);
  void remove_atom(int);
  float get_x(int);
  bool is_occupied(double*);
  void set_occupation(double*, bool);
  bool porus_side(int);

  float get_transfer_probability(float, int, int);
  double energy_full();



};
}

#endif
#endif