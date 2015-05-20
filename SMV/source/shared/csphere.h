// $Date: 2012-08-20 19:39:01 -0400 (Mon, 20 Aug 2012) $ 
// $Revision: 12156 $
// $Author: koverholt $

#ifndef CSPHERE_H_DEFINED
#define CSPHERE_H_DEFINED

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
