// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

char interp_revision[]="$Revision$";


/* ----------------------- interp3d ----------------------------- */

#define INTERP1D(f0,f1,dx) (float)((f0) + ((f1)-(f0))*(dx))
float interp3d(float *xplt, float *yplt, float *zplt, int ibar, int jbar, int kbar, float *vals, float xyz[3]){
  int i, j, k;
  int ijk;
  float val000,val100,val010,val110;
  float val001,val101,val011,val111;
  float val00,val01,val10,val11;
  float val0, val1, val;
  int nx, ny, nxy;
  float dx, dy, dz;
  float dxbar, dybar, dzbar;
  float *vv;

  dxbar = xplt[1]-xplt[0];
  dybar = yplt[1]-yplt[0];
  dzbar = zplt[1]-zplt[0];

  i = (xyz[0]-xplt[0])/dxbar;
  j = (xyz[1]-yplt[0])/dybar;
  k = (xyz[2]-zplt[0])/dzbar;

  dx = (xyz[0] - i*dxbar)/dxbar;
  dy = (xyz[1] - j*dybar)/dybar;
  dz = (xyz[2] - k*dzbar)/dzbar;

  nx = ibar;
  ny = jbar;
  nxy = nx*ny;

  ijk = i + j*nx + k*nxy;

  vv = vals + ijk;
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
}
