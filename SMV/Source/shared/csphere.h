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
