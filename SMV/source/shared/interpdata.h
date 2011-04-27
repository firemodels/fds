// $Date: 2011-03-24 16:22:30 -0400 (Thu, 24 Mar 2011) $ 
// $Revision: 7970 $
// $Author: gforney $

/* --------------------------  xyzdata ------------------------------------ */

typedef struct {
  int ibar, jbar, kbar;
  unsigned char *lightmap, *opacity;
  float *flightmap;
  float xbar0, xbar, ybar0, ybar, zbar0, zbar;
  float dx, dy, dz, dxyzmin;
} lightdata;

void create_lightmap(lightdata *lightinfo);
unsigned char interp3d(lightdata *lightinfo, float xyz[3]);
