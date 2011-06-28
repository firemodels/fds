// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "flowfiles.h"
#include "interp.h"
#include "smokeviewdefs.h"

char interp_revision[]="$Revision$";

/* ------------------ get_z_interp_factors ------------------------ */

void get_z_interp_factors(float *zplt, int nz, float z, int *k1, int *k2, float *f1, float *f2){
  float dz;
  int ileft, iright;

  dz = zplt[1] - zplt[0];

  ileft = (z-zplt[0])/dz;
  if(ileft<0)ileft=0;
  if(ileft>nz-1)ileft=nz-1;
  iright = ileft + 1;

  *k1 = ileft;
  *k2 = iright;
  *f1 = (z-zplt[ileft])/dz;
  *f2 = (zplt[iright]-z)/dz;
  return;
}

/* ------------------ interp3dsliceindex ------------------------ */

int interp3dsliceindex(unsigned char *data, float *zplt, int nz, int n0, float z){
  int k1, k2;
  int n1, n2;
  float dz;
  float val1, val2;
  float dz1, dz2;
  float z1, z2;
  int ival;

  dz = zplt[1] - zplt[0];

  k1 = (z-zplt[0])/dz;
  if(k1<0)k1=0;
  if(k1>nz-1)k1=nz-1;
  k2 = k1 + 1;

  val1 = data[n0+k1];
  val2 = data[n0+k2];
  z1 = zplt[k1];
  z2 = zplt[k2];
  ival = ((z-z1)*val2 + (z2-z)*val1)/dz;
  if(ival<0)ival=0;
  if(ival>255)ival=255;
  return ival;
}

