// $Date: 2011-03-24 16:22:30 -0400 (Thu, 24 Mar 2011) $ 
// $Revision: 7970 $
// $Author: gforney $

/* --------------------------  xyzdata ------------------------------------ */

typedef struct {
  int ibar, jbar, kbar;
  float *data, xbar0, xbar, ybar0, ybar, zbar0, zbar;
  float dx, dy, dz, dxyzmin;
} xyzdata;

float interp3d(xyzdata *xyzinfo, float xyz[3]);
