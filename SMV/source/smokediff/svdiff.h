#ifndef SVDIFF_H_DEFINED
#define SVDIFF_H_DEFINED
#include "histogram.h"

//************************** pre-processing directives ****************************************

#ifdef INMAIN
#define EXTERN
#else
#define EXTERN extern
#endif

#include "string_util.h"


#ifdef WIN32
#define STDCALLF extern void _stdcall
#else
#define STDCALLF extern void
#endif

#ifdef X64
#ifndef STRUCTSTAT
#define STRUCTSTAT struct __stat64
#endif
#ifndef STAT
#define STAT _stat64
#endif
#else
#ifndef STRUCTSTAT
#define STRUCTSTAT struct stat
#endif
#ifndef STAT
#define STAT stat
#endif
#endif

#ifndef FILE_SIZE
#define FILE_SIZE unsigned long long
#endif

//************************** data structures ****************************************

typedef struct {
  int ibar, jbar, kbar;
  float xbar0, xbar, ybar0, ybar, zbar0, zbar;
  float dx, dy, dz;
  float *xplt, *yplt, *zplt;
} meshdata;

typedef struct _boundary {
  char *file;
  int version;
  struct _boundary *boundary2;
  FILE_SIZE filesize;
  int npatches;
  int *pi1, *pi2, *pj1, *pj2, *pk1, *pk2, *patchdir;
  int *patch2index, *patchsize, *qoffset;
  char keyword[255];
  int boundarytype;
  histogramdata *histogram;
  meshdata *boundarymesh;
  flowlabels label;
} boundary;

typedef struct _slice {
  char *file;
  int is1, is2, js1, js2, ks1, ks2;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  FILE_SIZE filesize;
  int factor[3];
  int version;
  int volslice;
  struct _slice *slice2;
  char keyword[255];
  int slicetype;
  meshdata *slicemesh;
  histogramdata *histogram;
  flowlabels label;
} slice;

typedef struct _plot3d {
  char keyword[255];
  char *file;
  float time;
  struct _plot3d *plot3d2;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  histogramdata *histogram[5];
  meshdata *plot3dmesh;
  flowlabels labels[5];
} plot3d;

typedef struct {
  slice *sliceinfo;
  meshdata *meshinfo;
  plot3d *plot3dinfo;
  boundary *boundaryinfo;
  char *dir;
  int endian;
  int nmeshes;
  int nsliceinfo, nplot3dinfo, nboundary_files;
} casedata;

//************************** headers ****************************************

int getendian(void);
void Usage(void);
int mesh_match(meshdata *mesh1, meshdata *mesh2);
int ReadSMV(FILE *streamsmv, FILE *stream_out, casedata *smvcase);
void setup_boundary(FILE *stream_out);
void setup_slice(FILE *stream_out);
void setup_plot3d(FILE *stream_out);
plot3d *getplot3d(plot3d *plot3din, casedata *case2);
slice *getslice(slice *slicein, casedata *case2);
boundary *getboundary(boundary *boundaryin, casedata *case2);
void diff_boundaryes(FILE *stream_out);
void diff_slices(FILE *stream_out);
void diff_plot3ds(FILE *stream_out);
int similar_grid(meshdata *mesh1, meshdata *mesh2, int *factor);
int exact_grid(meshdata *mesh1, meshdata *mesh2, int *factor);
int getpatchindex(int in1, boundary *boundaryin, boundary *boundaryout);

#define FORTgetsliceparms _F(getsliceparms)
#define FORTclosefortranfile _F(closefortranfile)
#define FORTopenslice _F(openslice)
#define FORTgetsliceframe _F(getsliceframe)
#define FORToutsliceframe _F(outsliceframe)
#define FORToutsliceheader _F(outsliceheader)
#define FORTgetplot3dq _F(getplot3dq)
#define FORTplot3dout _F(plot3dout)
#define FORTgetboundaryheader1 _F(getboundaryheader1)
#define FORTgetboundaryheader2 _F(getboundaryheader2)
#define FORTopenboundary _F(openboundary)
#define FORTgetpatchdata _F(getpatchdata)
#define FORToutboundaryheader _F(outboundaryheader)
#define FORToutpatchframe _F(outpatchframe)
#define FORTendianout _F(endianout)
#define FORTget_file_unit _F(get_file_unit)

STDCALLF FORTget_file_unit(int *file_unit, int *file_unit_start);
STDCALLF FORToutpatchframe(int *lunit, int *npatch,
                          int *pi1, int *pi2, int *pj1, int *pj2, int *pk1, int *pk2,
                          float *patchtime, float *pqq, int *error);
STDCALLF FORToutboundaryheader(char *outfile, int *unit3, int *npatches,
                              int *pi1, int *pi2, int *pj1, int *pj2, int *pk1, int *pk2,
                              int *patchdir, int *error1, FILE_SIZE len);
STDCALLF FORTgetpatchdata(int *lunit, int *npatch,int *pi1,int *pi2,int *pj1,int *pj2,int *pk1,int *pk2,
                         float *patch_times,float *pqq, int *npqq, int *error);
STDCALLF FORTopenboundary(char *boundaryfilename, int *boundaryunitnumber,
                         int *version, int *error, FILE_SIZE len);
STDCALLF FORTgetboundaryheader1(char *boundaryfilename, int *boundaryunitnumber,
                               int *npatch, int *error, FILE_SIZE lenfile);
STDCALLF FORTgetboundaryheader2(int *boundaryunitnumber, int *version, int *npatches,
                               int *pi1, int *pi2, int *pj1, int *pj2, int *pk1, int *pk2, int *patchdir);
STDCALLF FORTgetsliceframe(int *lu11,
                          int *is1,int *is2,int *js1,int *js2,int *ks1,int *ks2,
                          float *time,float *qframe,int *slicetest, int *error);
STDCALLF FORTgetsliceparms(char *file,
                          int *is1,int *is2,int *js1,int *js2,int *ks1, int *ks2,
                          int *ni, int *nj, int *nk,
                          int *slice3d, int *error,FILE_SIZE lenfile);
STDCALLF FORTopenslice(char *slicefilename, int *unit,
                      int *is1, int *is2, int *js1, int *js2, int *ks1, int *ks2,
                      int *error, FILE_SIZE lenfile);
STDCALLF FORTclosefortranfile(int *unit);
STDCALLF FORToutsliceframe(int *unit3,
                          int *is1a,int *is1b,int *js1a,int *js1b,int *ks1a,int *ks1b,
                          float *time1,float *qframeout, int *error);
STDCALLF FORToutsliceheader(char *outfile,int *unit3,
                             int *is1a,int *is2a,int *js1a,int *js2a,int *ks1a,int *ks2a,
                             int *error1,int len);
STDCALLF FORTgetplot3dq(char *qfilename, int *nx, int *ny, int *nz, float *qq, int *error, int *isotest, int len);
STDCALLF FORTplot3dout(char *outfile,int *nx,int *ny,int *nz,float *qout,int *error3,int lenout);
STDCALLF FORTendianout(char *endian_filename,int lenout);

//************************** global variables ****************************************

EXTERN char pp[2];
EXTERN casedata *caseinfo;
EXTERN char *sourcedir1, *sourcedir2, *destdir;
EXTERN int test_mode, display_warnings;
EXTERN char type_label[1024];
EXTERN FILE *LOG_FILENAME;

#endif
