// $Date: 2011-03-24 16:22:30 -0400 (Thu, 24 Mar 2011) $ 
// $Revision: 7970 $
// $Author: gforney $

/* --------------------------  xyzdata ------------------------------------ */

typedef struct {
  int *ijkbar;
  unsigned char *radiance, *opacity;
  float *xyzbar0, *xyzbar;
  float *dxyz;
} radiancedata;

void setup_radiancemap(radiancedata *radianceinfo, int ijkbar[3], float xyzbar0[3], float xyzbar[3], float dxyz[3], 
                       unsigned char *radiance, unsigned char *opacity);
void build_radiancemap(radiancedata *radianceinfo);
float interp3d(float *xplt, float *yplt, float *zplt, int ibar, int jbar, int kbar, float *vals, float xyz[3]);
