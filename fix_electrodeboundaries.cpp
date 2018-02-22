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

	idregion(NULL), scalestr(NULL){

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


	// // optional arguments
	// iregion = -1;
 //  	idregion = NULL;
 //  	scale = 1.0;

 //  	int iarg = 10;	//start after madatory arguments
 //  	while (iarg < narg) {

 //    if (strcmp(arg[iarg],"region") == 0) { //keyword = region
 //      if (iarg+2 > narg) error->all(FLERR,"Illegal fix imagecharges command"); //check to make sure there is a next argument
 //      iregion = domain->find_region(arg[iarg+1]);
 //      if (iregion == -1)
 //        error->all(FLERR,"Region ID for fix imagescharges does not exist");
 //      int n = strlen(arg[iarg+1]) + 1;
 //      idregion = new char[n];
 //      strcpy(idregion,arg[iarg+1]);
 //      iarg += 2;

 //    } else if (strcmp(arg[iarg],"scale") == 0) { //keyword = scale
 //    	if (iarg+2 > narg) error->all(FLERR,"Illegal fix imagecharges command"); //check to make sure there is a next argument
 //    	if (strstr(arg[iarg+1],"v_") == arg[iarg+1]) { //variable
 //    		int n = strlen(&arg[iarg+1][2]) + 1;
	//     	scalestr = new char[n];
	//     	strcpy(scalestr,&arg[iarg+1][2]);
	// 		} else {
	// 			scale = force->numeric(FLERR,arg[iarg+1]);
	// 			scalestyle = CONSTANT;
	// 		}
	// 		iarg += 2;
    
 //    } else error->all(FLERR,"Illegal fix imagecharges command"); // not a recognized keyword
 //  }
	// atom->add_callback(0);

}