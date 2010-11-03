// $Date$ 
// $Revision$
// $Author$

#ifndef DEF_ISOTEST2
#define DEF_ISOTEST2
#define INCPOINTS 100000

#if defined(WIN32)
#include <windows.h>
#endif

#ifdef pp_DRAWISO
#ifdef pp_OSX
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif
#endif

#ifdef IN_ISOBOX
#define SV_EXTERN
#else
#define SV_EXTERN extern
#endif

/* iso-surface definitions */

typedef struct {
  float level;
  float *color;
  int dataflag;
  int  nvertices, ntriangles, nnorm;      /* actual number */
  int maxvertices, maxtriangles, maxnorm; /* space reserved */
  int colorindex;
  int normtype, defined, plottype;
  int *triangles;
  int *closestnodes;
  unsigned short *triangles2;
  unsigned char *triangles1;
  unsigned short *vertices, *tvertices;
  unsigned char *color8;
  int *sortedlist,*rank;
  float xmin, ymin, zmin, xyzmaxdiff;
  float *xvert, *yvert, *zvert, *tvert;
  float tmin, tmax;
  float *xnorm, *ynorm, *znorm;
  short *norm, *vertexnorm;
  unsigned char *comp_bufferframe, *full_bufferframe;
  int ncomp_bufferframe, nfull_bufferframe;
  unsigned char *s_norm;
  int cullfaces;
  int compression_type;
} isosurface;

typedef struct {
  int index;
  unsigned short vertex[3];
} sortdata;

typedef struct {
  int index;
  int sortedlist;
} rankdata;

typedef struct {
  int index;
  int closest;
} orderdata;



#ifndef pp_DRAWISO

#ifndef pp_noappend
#define CCisosurface2file iso2file_
#define CCisosurfacet2file isot2file_
#define CCisoheader isoheader_
#define CCtisoheader tisoheader_
#else
#define CCisosurface2file iso2file
#define CCisosurfacet2file isot2file
#define CCtisoheader tisoheader
#define CCisoheader isoheader
#endif
SV_EXTERN void CCisoheader(char *isofile, 
                 char *isolonglabel, char *isoshortlabel, char *isounits,
                 float *levels, int *nlevels, int *error);
SV_EXTERN void CCtisoheader(char *isofile, 
                 char *isolonglabel, char *isoshortlabel, char *isounits,
                 float *levels, int *nlevels, int *error);
SV_EXTERN void isoout(FILE *isostream,float t, int timeindex, isosurface *surface, int *error);
SV_EXTERN void CCisosurface2file(char *isofile, float *t, float *data, int *iblank, float *level, int *nlevels,
     float *xplt, int *nx, float *yplt, int *ny, float *zplt, int *nz,int *reduce_triangles, int *error);
SV_EXTERN void CCisosurfacet2file(char *isofile, float *t, float *data, int *data2flag, float *data2, int *iblank, 
						float *level, int *nlevels,
                   float *xplt, int *nx, 
                   float *yplt, int *ny, 
                   float *zplt, int *nz,
                   int *reduce_triangles, int *error);
#endif
SV_EXTERN int CompressIsosurface(isosurface *surface, int reduce_triangles, 
                        float xmin, float xmax,
                        float ymin, float ymax,
                        float zmin, float zmax);
SV_EXTERN int UpdateIsosurface(isosurface *surface, 
                      const float *xvert, 
                      const float *yvert, 
                      const float *zvert, 
                      const float *tvert, 
                      const int *closestnodes, 
                      int nvert, 
                      const int *triangles, 
                      int ntriangles);

SV_EXTERN void DrawIsosurface(const isosurface *isodata);
SV_EXTERN void freesurface(isosurface *surfacedata);
SV_EXTERN void InitIsosurface(isosurface *surfacedata, float level, float *color, int colorindex);
SV_EXTERN int ResizeSurface(isosurface *surfacedata, int incvert, int inctrilist, int incnorm);
SV_EXTERN void GetIsobox(const float *x, 
               const float *y, 
               const float *z, 
               const float *vals, 
               const float *tvals, 
               const int *nodeindexes, 
               float level,
               float *xvert, 
               float *yvert, 
               float *zvert, 
               float *tvert, int *closestnodes, int *nvert,
               int *triangles, int *ntriangles);
int GetIsosurface(isosurface *surface, 
                  const float *data, 
                  const float *tdata, 
                  const char *iblank_cell, 
                  float level,
                  const float *xplt, int nx, 
                  const float *yplt, int ny, 
                  const float *zplt, int nz
                   );

SV_EXTERN void ReduceToUnit(float v[3]);
SV_EXTERN void calcNormal(const float *v1, const float *v2, const float *v3, float *out);
SV_EXTERN void calcNormal2(const unsigned short *v1, 
                           const unsigned short *v2, 
                           const unsigned short *v3, 
                           float *out, float *area);
SV_EXTERN int GetNormalSurface(isosurface *surfacedata);
#endif





