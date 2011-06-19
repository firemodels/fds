// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "interp.h"
#include "smokeviewdefs.h"

char interp_revision[]="$Revision$";


/* ----------------------- interp3d ----------------------------- */

#define INTERP1D(f0,f1,dx) (float)((f0) + ((f1)-(f0))*(dx))
float interp3d(float xyz[3], float *vals, float *xplt, float *yplt, float *zplt, int ibar, int jbar, int kbar){
  int i, j, k;
  int ijk;
  float val000,val100,val010,val110;
  float val001,val101,val011,val111;
  float val00,val01,val10,val11;
  float val0, val1, val;
  int nx, ny, nz, nxy;
  float dx, dy, dz;
  float dxbar, dybar, dzbar;
  float *vv;

  dxbar = xplt[1]-xplt[0];
  dybar = yplt[1]-yplt[0];
  dzbar = zplt[1]-zplt[0];

  nx = ibar + 1;
  ny = jbar + 1;
  nz = kbar + 1;
  nxy = nx*ny;

  GETINDEX(i,xyz[0],xplt[0],dxbar,ibar);
  GETINDEX(j,xyz[1],yplt[0],dybar,jbar);
  GETINDEX(k,xyz[2],zplt[0],dzbar,kbar);

  //ijk = i + j*nx + k*nxy;
  ijk = k + (j + i*ny)*nz; 

  return vals[ijk];
  /*
  dx = (xyz[0] - i*dxbar)/dxbar;
  dy = (xyz[1] - j*dybar)/dybar;
  dz = (xyz[2] - k*dzbar)/dzbar;
  val000 = vv[0];
  val100 = vv[1];
  val010 = vv[nx];
  val110 = vv[1+nx];
  vv += nxy;
  val001 = vv[0];
  val101 = vv[1];
  val011 = vv[nx];
  val111 = vv[1+nx];
  val00 = INTERP1D(val000,val100,dx);
  val10 = INTERP1D(val010,val110,dx);
  val01 = INTERP1D(val001,val101,dx);
  val11 = INTERP1D(val011,val111,dx);
   val0 = INTERP1D( val00, val10,dy);
   val1 = INTERP1D( val01, val11,dy);
    val = INTERP1D(  val0,  val1,dz);

  return val;
  */
}

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

