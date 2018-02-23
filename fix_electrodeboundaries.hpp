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

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  void set_arrays(int);


 protected: 
  double pxvalue,pyvalue,pzvalue,nxvalue,nyvalue,nzvalue,distvalue,vvalue,dvvalue;
  int varflag,iregion,itype;
  char *pxstr,*pystr,*pzstr,*nxstr,*nystr,*nzstr,*diststr,*vstr,*dvstr;
  char *idregion, *scalestr;
  int pxvar,pyvar,pzvar,nxvar,nyvar,nzvar,distvar; 
  int pxstyle,pystyle,pzstyle,nxstyle,nystyle,nzstyle,diststyle,vstyle,dvstyle;


  double foriginal[3],foriginal_all[3];
  int force_flag;
  int nlevels_respa,ilevel_respa;

  int maxatom;
  double **sforce;

};
}

#endif
#endif