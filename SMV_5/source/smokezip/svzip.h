// $Date$ 
// $Revision$
// $Author$

//***********************
//************* #definess
//***********************
#ifdef INMAIN
#define EXTERN
#else
#define EXTERN extern
#endif
#include "csphere.h"

#define UNLINK unlink

#include "histogram.h"
#include "threader.h"

#ifdef pp_PART
#define rgb_white 12
#define rgb_yellow 13
#define rgb_blue 14
#define rgb_red 15
#define rgb_green 16
#define rgb_magenta 17
#define rgb_cyan 18
#define rgb_black 19
#endif

#ifdef X64
#define STRUCTSTAT struct __stat64
#define STAT _stat64
#else
#define STRUCTSTAT struct stat
#define STAT stat
#endif

//***********************
//************* structures
//***********************

/* --------------------------  mesh ------------------------------------ */


typedef struct {
  int ibar, jbar, kbar;
  float *xplt, *yplt, *zplt;
  float *xpltcell, *ypltcell, *zpltcell;
  float xbar0, xbar, ybar0, ybar, zbar0, zbar;
  float dx, dy, dz;
  float dxx, dyy, dzz;
} mesh;

/* --------------------------  flowlabels ------------------------------------ */

typedef struct {
  char *longlabel, *shortlabel, *unit;
} flowlabels;

/* --------------------------  patch ------------------------------------ */

typedef struct {
  char *file,*filebase;
  int unit_start;
  char summary[1024];
  int compressed;
  int filesize;
  int inuse,inuse_getbounds;
  int seq_id, autozip;
  int doit, done;
  int *pi1, *pi2, *pj1, *pj2, *pk1, *pk2, *patchdir, *patchsize;
  int npatches;
  int setvalmin, setvalmax;
  float valmin, valmax;
  int version;
  histogramdata *histogram;
  flowlabels label;
  int dup;
} patch;

/* --------------------------  iso ------------------------------------ */

typedef struct {
  int blocknumber;
  int unit_start;
  char summary[1024];
  int compressed;
  int inuse;
  char *file, *filebase;
  flowlabels label;
  int filesize;
  int seq_id, autozip;
  int dup;
  int version;
  int dataflag;
  float tmin, tmax;
  float *isolevels;
  int nisolevels,nisosteps;
} iso;

/* --------------------------  slice ------------------------------------ */

typedef struct {
  char *file,*filebase;
  int unit_start;
  char summary[1024];
  int compressed;
  int inuse,inuse_getbounds;
  int filesize;
  int seq_id,autozip;
  int doit, done, count;
  int setvalmin, setvalmax;
  float valmin, valmax;
  int setchopvalmin, setchopvalmax;
  float chopvalmin, chopvalmax;
  int version;
  histogramdata *histogram;
  flowlabels label;
  int dup;
  int rle;
} slice;

/* --------------------------  bound ------------------------------------ */

typedef struct {
  int setvalmin, setvalmax;
  float valmin, valmax;
} bound;

/* --------------------------  plot3d ------------------------------------ */

typedef struct {
  char *file,*filebase;
  int unit_start;
  char summary[1024];
  int compressed;
  int inuse;
  float time;
  int blocknumber;
  mesh *plot3d_mesh;
  int filesize;
  int seq_id,autozip;
  int doit, done, count;
  bound bounds[5];
  int version;
  flowlabels labels[5];
  int dup;
} plot3d;

/* --------------------------  vert ------------------------------------ */

typedef struct {
  float normal[3];
} vert;

/* --------------------------  smoke3d ------------------------------------ */

typedef struct {
  char *file,*filebase;
  int unit_start;
  char summary[1024];
  int compressed;
  int inuse;
  int seq_id, autozip;
  int nx, ny, nz, filesize;
  unsigned char *compressed_lightingbuffer;
  uLongf ncompressed_lighting_zlib;
} smoke3d;

#ifdef pp_PART

/* --------------------------  part5prop ------------------------------------ */

typedef struct {
  int used;
  char isofilename[1024];
  float *partvals;
  flowlabels label;
  float valmin, valmax;
  histogramdata *histogram;
  int setvalmin, setvalmax;
} part5prop;

/* --------------------------  partclass ------------------------------------ */

typedef struct {
  char *name;
  int ntypes;
  flowlabels *labels;
} part5class;

/* --------------------------  part5data ------------------------------------ */

typedef struct {
  int npoints,n_rtypes, n_itypes;
  int *tags,*sort_tags;
  float *rvals;
  unsigned char *irvals;
  unsigned char **cvals;
} part5data;

/* --------------------------  part ------------------------------------ */

typedef struct {
  char *file,*filebase;
  char summary[1024], summary2[1024];
  int unit_start;
  char **summaries;
  int nsummaries;
  int compressed,compressed2;
  int inuse,inuse_part2iso;
  int filesize;
  int seq_id, autozip;
  int setvalmin, setvalmax;
  float valmin, valmax;
  flowlabels *label;
  mesh *partmesh;

  int nclasses;
  part5class **classptr;
  part5data *data5;
} part;
#endif

#define PDFMAX 100000
typedef struct {
  int ncount;
  int buckets[PDFMAX];
  float pdfmin,pdfmax;
} pdfdata;

#define BOUND(x,xmin,xmax) (x<xmin)?xmin:(x>xmax)?xmax:x
#define GET_INTERVAL(xyz,xyz0,dxyz) ((xyz)-(xyz0))/(dxyz)


//***********************
//************* headers
//***********************

void print_summary(void);
void *compress_all(void *arg);
void mt_compress_all(void);
void rand_absdir(float xyz[3], int dir);
void rand_cone_dir(float xyz[3], float dir[3], float mincosangle);
void rand_sphere_dir(float xyz[3]);
float rand_1d(float xmin, float xmax);
void rand_2d(float xy[2], float xmin, float xmax, float ymin, float ymax);
void rand_3d(float xyz[3], float xmin, float xmax, float ymin, float ymax, float zmin, float zmax);
void get_startup_slice(int seq_id);
void get_startup_iso(int seq_id);
void get_startup_smoke(int seq_id);
void get_startup_patch(int seq_id);
int getrevision(char *svn);
void getSMZversion(char *SMZversion);
int getmaxrevision(void);
void version(void);
unsigned int irle(unsigned char *buffer_in, int nchars_in, unsigned char *buffer_out);
int readsmv(char *file);
int getendian(void);
int convert_slice(slice *slicei, int *thread_index);
slice *getslice(char *string);
void *compress_slices(void *arg);
void *compress_isos(void *arg);
int plot3ddup(plot3d *plot3dj, int iplot3d);
int slicedup(slice *slicej, int islice);
void *compress_plot3ds(void *arg);
void getfilesizelabel(int size, char *sizelabel);
void initpdf(pdfdata *pdf);
int getfileinfo(char *filename, char *sourcedir, int *filesize);
void filecopy(char *destdir, char *file, char *filebase);
void copyfile(char *destfile, char *sourcefile);
void makesvd(char *destdir, char *smvfile);
void trimzeros(char *line);
void usage(char *prog);
void getpdf(float *vals, int nvals, pdfdata *pdf);
void mergepdf(pdfdata *pdf1, pdfdata *pdf2, pdfdata *pdfmerge);
void smoothlabel(float *a, float *b, int n);
#ifdef pp_PART
void compress_parts(void *arg);
void *convert_parts2iso(void *arg);
part *getpart(char *string);
part5prop *getpartprop(char *string);
int getpartprop_index(char *string);
void convert_part(part *parti, int *thread_index);
int convertable_part(part *parti);
#endif
void *compress_patches(void *arg);
patch *getpatch(char *string);
char *trim_front(char *line);
int patchdup(patch *patchj, int ipatch);
int readlabels(flowlabels *flowlabel, FILE *stream);
void readini(char *file);
void readini2(char *file2);
void Get_Boundary_Bounds(void);
void Get_Slice_Bounds(void);
#ifdef pp_PART
void Get_Part_Bounds(void);
#endif
void convert_3dsmoke(smoke3d *smoke3di, int *thread_index);
void *compress_smoke3ds(void *arg);
int match(const char *buffer, const char *key, unsigned int lenkey);
void trim(char *line);
void Normal(unsigned short *v1, unsigned short *v2, unsigned short *v3, float *normal, float *area);
float atan3(float y, float x);

#ifdef pp_noappend
#define FORTgetpartheader1 getpartheader1
#define FORTgetpartheader2 getpartheader2
#define FORTgetpartdataframe getpartdataframe
#define FORTclosefortranfile closefortranfile
#define FORTgetboundaryheader1 getboundaryheader1
#define FORTgetboundaryheader2 getboundaryheader2
#define FORTopenboundary openboundary
#define FORTgetpatchdata getpatchdata
#define FORTopenslice openslice
#define FORTopenpart openpart
#define FORTgetsliceframe getsliceframe
#define FORTget_file_unit get_file_unit
#else
#define FORTgetpartheader1 getpartheader1_
#define FORTgetpartheader2 getpartheader2_
#define FORTgetpartdataframe getpartdataframe_
#define FORTclosefortranfile closefortranfile_
#define FORTgetboundaryheader1 getboundaryheader1_
#define FORTgetboundaryheader2 getboundaryheader2_
#define FORTopenboundary openboundary_
#define FORTgetpatchdata getpatchdata_
#define FORTopenslice openslice_
#define FORTopenpart openpart_
#define FORTgetsliceframe getsliceframe_
#define FORTget_file_unit get_file_unit_
#endif
#ifdef WIN32
#define STDCALL extern void _stdcall
#else
#define STDCALL extern void
#endif

STDCALL FORTget_file_unit(int *file_unit,int *file_unit_start);
STDCALL FORTopenpart(char *partfilename, int *unit, int *endian, int *error, FILE_SIZE lenfile);
STDCALL FORTgetpartheader1(int *unit, int *nclasses, int *fdsversion, int *size);
STDCALL FORTgetpartheader2(int *unit, int *nclasses, int *nquantities, int *size);
STDCALL FORTgetpartdataframe(int *unit, int *nclasses, int *nquantities, int *npoints, float *time, int *tagdata, float *pdata, int *size, int *error);

STDCALL FORTclosefortranfile(int *lunit);

STDCALL FORTgetpatchdata(int *lunit, int *npatch,int *pi1,int *pi2,int *pj1,int *pj2,int *pk1,int *pk2,
                         float *patchtimes,float *pqq, int *ndummy, int *error);
STDCALL FORTopenboundary(char *boundaryfilename, int *boundaryunitnumber, 
                         int *endian, int *version, int *error, FILE_SIZE len);
STDCALL FORTgetboundaryheader1(char *boundaryfilename, int *boundaryunitnumber, 
                               int *endian, int *npatch, int *error, FILE_SIZE lenfile);
STDCALL FORTgetboundaryheader2(int *boundaryunitnumber, int *version, int *npatches,
                               int *pi1, int *pi2, int *pj1, int *pj2, int *pk1, int *pk2, int *patchdir);

STDCALL FORTgetsliceframe(int *lu11,
                          int *is1,int *is2,int *js1,int *js2,int *ks1,int *ks2,
                          float *time,float *qframe,int *slicetest, int *error);
STDCALL FORTopenslice(char *slicefilename, int *unit, int *endian, 
                      int *is1, int *is2, int *js1, int *js2, int *ks1, int *ks2,
                      int *error, FILE_SIZE lenfile);

//***********************
//************* variables
//***********************

EXTERN int doit_smoke3d, doit_boundary, doit_slice, doit_plot3d, doit_iso;
#ifdef pp_PART2
EXTERN int doit_particle;
#endif

EXTERN int first_initsphere,first_slice,first_patch,first_plot3d,first_part2iso,first_part2iso_smvopen;
EXTERN int frameskip;
EXTERN int no_chop;
EXTERN patch *patchinfo;
EXTERN mesh *meshinfo;
EXTERN smoke3d *smoke3dinfo;
EXTERN int npatch_files;

EXTERN int nslice_files, niso_files, nplot3d_files;

EXTERN slice *sliceinfo;
EXTERN iso *isoinfo;
EXTERN plot3d *plot3dinfo;

EXTERN int nmeshes;
EXTERN int overwrite_slice;
EXTERN int overwrite_plot3d;
#ifdef pp_PART
EXTERN int overwrite_part;
EXTERN part *partinfo;
EXTERN int npart_files;
EXTERN int npartclassinfo;
EXTERN part5class *partclassinfo;
EXTERN part5prop *part5propinfo;
EXTERN int maxpart5propinfo, npart5propinfo;
#endif
EXTERN int nsmoke3d_files;
EXTERN int endianswitch,overwrite_b,overwrite_s;
EXTERN int overwrite_iso;
EXTERN int cleanfiles;
EXTERN char *destdir,*sourcedir;
EXTERN int lensourcedir,lendestdir;
EXTERN char dirseparator[3];
EXTERN char pp[2];
EXTERN int smoke3dzipstep, boundzipstep, slicezipstep;
EXTERN int isozipstep,doiso;
EXTERN int filesremoved;
EXTERN int endf, syst;
EXTERN char endianfilebase[1024];
EXTERN char *endianfile;
EXTERN spherepoints sphereinfo;
EXTERN int autozip, make_demo;
EXTERN int get_bounds, get_slice_bounds, get_plot3d_bounds, get_boundary_bounds;
#ifdef pp_PART
EXTERN int get_part_bounds;
EXTERN int partfile2iso;
#endif
EXTERN char smvisofile[1024];
EXTERN threaddata *threadinfo;

