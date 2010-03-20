// $Date$ 
// $Revision$
// $Author$

#ifndef DEF_CONTOURTEST
#define DEF_CONTOURTEST

#if defined(WIN32)
#include <windows.h>
#endif
#ifdef pp_OSX
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif


typedef struct {
  float *levels;
  int nlevels, *nnodes, *npolys, *nlines;
  int **polysize;
  float **xnode, **ynode, **xlines, **ylines, xyzval;
  float **rgbptr;
  int idir;
} contour;

void getcontournodes(int n, int levels, const double x[4], const double y[4], const double val[4],
                     double contlow, int modelow, double conthigh, int modehigh,
                     int *nnode, float *xnode, float *ynode,
                     int *nnode2,float *xline, float *yline,
                     int *casen,int blankit);
void getlinecontournodes(double linelevel, const double x[4], const double y[4], const double val[4],
                     int *nline_nodes,float *xline, float *yline,
                     int blankit);

void DrawContours(const contour *ci);
void setcontourslice(contour *ci,int idir,float xyz);
void getcontours(const float *xgrid, const float *ygrid, int nx, int ny,  
                 const float *vals, const char *iblank, const float *levels,  
                 const contour *ci);
void initcontour(contour *ci, float **rgbptr, int nlevels);
void freecontour(contour *ci);

#endif





