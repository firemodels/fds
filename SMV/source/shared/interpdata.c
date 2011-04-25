// $Date: 2011-03-30 12:04:57 -0400 (Wed, 30 Mar 2011) $ 
// $Revision: 8021 $
// $Author: gforney $

#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "interpdata.h"

// svn revision character string
char interpdata_revision[]="$Revision: 8021 $";

/* ----------------------- integrate3d ----------------------------- */

float getalpha(xyzdata *xyzinfo, float ksoot, float xyz0[3], float xyz1[3]){
  float lpath, dlpath, xyz[3];
  float dx, dy, dz;
  int i, nl;
  float xnl;
  float sum=0.0;
  float val;
  float alpha;

  dx = xyz1[0]-xyz0[0];
  dy = xyz1[1]-xyz0[1];
  dz = xyz1[2]-xyz0[2];

  lpath = (float)sqrt(dx*dx+dy*dy+dz*dz);
  dlpath = lpath/xyzinfo->dxyzmin;
  nl = lpath/dlpath+1;
  xnl = (float)nl;
  dlpath=lpath/(float)nl;
  for(i=0;i<nl;i++){
    float ii;

    ii = (float)i + (float)0.5;
    xyz[0] = (xyz0[0]*(xnl-ii) + ii*xyz1[0])/xnl;
    xyz[1] = (xyz0[1]*(xnl-ii) + ii*xyz1[1])/xnl;
    xyz[2] = (xyz0[2]*(xnl-ii) + ii*xyz1[2])/xnl;
    val = interp3d(xyzinfo,xyz);
    sum+=val;
  }
  sum*=ksoot/xnl;
  alpha = (float)exp(-sum);
  return val;
}

/* ----------------------- interp3d ----------------------------- */

float interp3d(xyzdata *xyzinfo, float xyz[3]){
  int i, j, k;
  int ijk;

  i = (xyz[0]-xyzinfo->xbar0)/xyzinfo->dx;
  j = (xyz[1]-xyzinfo->ybar0)/xyzinfo->dy;
  k = (xyz[2]-xyzinfo->zbar0)/xyzinfo->dz;

  ijk = i + j*xyzinfo->ibar + k*xyzinfo->ibar*xyzinfo->jbar;

  return xyzinfo->data[ijk];
}
