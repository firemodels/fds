// $Date: 2007-10-07 22:08:47 -0400 (Sun, 07 Oct 2007) $ 
// $Revision: 800 $
// $Author: gforney $

#ifndef SETCSPHERE
#define SETCSPHERE

typedef struct {
  int n;
  float rad;
  float dphi;
  unsigned int npoints;
  float *dtheta;
  int *nlong;
  unsigned int *vallist;
  float *normals;
  short *snormals;
  float pi;
  float maxerr_deg;
} spherepoints;

#ifdef pp_DRAWISO
void drawspherepoints(spherepoints *spherei);
#endif
void initspherepoints(spherepoints *sphereinfo, int n);
void freespherepoints(spherepoints *sphereinfo);
unsigned int getnormalindex(spherepoints *sphereinfo, float *normal);
void getnormalvector(spherepoints *sphereinfo, unsigned int index, float *normal);
float *getnormalvectorptr(spherepoints *sphereinfo, unsigned int index);

#endif
