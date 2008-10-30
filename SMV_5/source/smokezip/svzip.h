// $Date$ 
// $Revision$
// $Author$
#ifdef INMAIN
#define EXTERN
#else
#define EXTERN extern
#endif
#include "csphere.h"
#define PERCENTILE_MIN 0
#define SET_MIN 1
#define GLOBAL_MIN 2

#define PERCENTILE_MAX 0
#define SET_MAX 1
#define GLOBAL_MAX 2

#define NBUCKETS 100000

#define UNLINK unlink

#ifdef pp_LIGHT
#define NPHOTONS 100000
#define NRAD 10
#define NTHETA 10
#define NPSI 10
#endif

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


/* --------------------------  flowlabels ------------------------------------ */

typedef struct {
  char *longlabel, *shortlabel, *unit;
} flowlabels;

/* --------------------------  patch ------------------------------------ */

typedef struct {
  char *file,*filebase;
  int filesize;
  int seq_id, autozip;
  int doit, done;
  int setvalmin, setvalmax;
  float valmin, valmax;
  int version;
  flowlabels label;
  int dup;
} patch;

typedef struct {
  int blocknumber;
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


typedef struct {
  char *file,*filebase;
  int filesize;
  int seq_id,autozip;
  int doit, done, count;
  int setvalmin, setvalmax;
  float valmin, valmax;
  int version;
  flowlabels label;
  int dup;
  int rle;
} slice;

typedef struct {
  int ibar, jbar, kbar;
  float *xplt, *yplt, *zplt;
  float xbar0, xbar, ybar0, ybar, zbar0, zbar;
  float dx, dy, dz;
  float dxx, dyy, dzz;
#ifdef pp_LIGHT
  float cell_volume, cell_surface_area, cell_cross_sectional_area;
  float *photon_cell;
  float *light_cell_radiance;
  float dxyzmax;
#endif
} mesh;

typedef struct {
  float normal[3];
} vert;

/* --------------------------  smoke3d ------------------------------------ */

typedef struct {
  char *file,*filebase;
  int seq_id, autozip;
  int nx, ny, nz, filesize;
#ifdef pp_LIGHT
  mesh *smoke_mesh;
  float *light_q_rect;
  int type;
  flowlabels label;
#endif
  unsigned char *compressed_lightingbuffer;
  uLongf ncompressed_lighting_zlib;
} smoke3d;

#ifdef pp_PART
/* --------------------------  part ------------------------------------ */

typedef struct {
  char *file,*filebase;
  int filesize;
  int seq_id, autozip;
  int setvalmin, setvalmax;
  float valmin, valmax;
  flowlabels label;
  int dup;
} part;
#endif

#ifdef pp_LIGHT
typedef struct {
  int type,dir;
  int nstep;
  int move;
  float t1, t2;
  float radius, area;
  float xyz1[3], xyz2[3], q, qflux;
} lightdata;
#endif

#define PDFMAX 100000
typedef struct {
  int ncount;
  int buckets[PDFMAX];
  float pdfmin,pdfmax;
} pdfdata;

#define BOUND(x,xmin,xmax) (x<xmin)?xmin:(x>xmax)?xmax:x
#define GET_INTERVAL(xyz,xyz0,dxyz) ((xyz)-(xyz0))/(dxyz)


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
int convert_slice(slice *slicei);
slice *getslice(char *string);
void compress_slices(void);
void compress_isos(void);
int slicedup(slice *slicej, int islice);
void compress_slices(void);
void getfilesizelabel(int size, char *sizelabel);
void initpdf(pdfdata *pdf);
int getfileinfo(char *filename, char *sourcedir, int *filesize);
void filecopy(char *destdir, char *file, char *filebase);
void makesvd(char *destdir, char *smvfile);
void trimzeros(char *line);
void usage(char *prog);
void getpdf(float *vals, int nvals, pdfdata *pdf);
void mergepdf(pdfdata *pdf1, pdfdata *pdf2, pdfdata *pdfmerge);
void smoothlabel(float *a, float *b, int n);
#ifdef pp_PART
unsigned char getpartcolor(float val);
void compress_parts(void);
part *getpart(char *string);
int partdup(part *partj, int ipart);
void convert_part(part *parti);
#endif
void compress_patches(void);
patch *getpatch(char *string);
char *trim_front(char *line);
int patchdup(patch *patchj, int ipatch);
int readlabels(flowlabels *flowlabel, FILE *stream);
void readini(char *file);
void readini2(char *file2);
int convert_boundary(patch *patchi, int pass);
void convert_3dsmoke(smoke3d *smoke3di);
void compress_smoke3ds(void);
int match(const char *buffer, const char *key, unsigned int lenkey);
void trim(char *line);
void Normal(unsigned short *v1, unsigned short *v2, unsigned short *v3, float *normal, float *area);
float atan3(float y, float x);
#ifdef pp_LIGHT
void light_smoke(smoke3d *smoke3di,unsigned char *full_lightingbuffer, float *val_buffer, unsigned char *alpha_buffer);
void set_lightfield(smoke3d *smoke3di,float xyz[3], float hrr);
void update_lightfield(float time, smoke3d *smoke3di, unsigned char *lightingbuffer);
#endif

EXTERN int frameskip;
EXTERN unsigned char *full_alphabuffer;
#ifdef pp_LIGHT
EXTERN float *full_logalphabuffer;
EXTERN int nphotons;
#endif
EXTERN patch *patchinfo;
EXTERN mesh *meshinfo;
EXTERN smoke3d *smoke3dinfo;
EXTERN int npatch_files;
#ifdef pp_LIGHT
EXTERN int make_lighting_file;
EXTERN int nlightinfo;
EXTERN lightdata *lightinfo;
EXTERN float light_delta;
EXTERN float *light_cdf;
#endif
EXTERN slice *sliceinfo;
EXTERN int nslice_files;
EXTERN int niso_files;
EXTERN iso *isoinfo;
EXTERN int nmeshes;
EXTERN int overwrite_slice;
#ifdef pp_PART
EXTERN part *partinfo;
EXTERN int npart_files;
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
#ifdef pp_LIGHT
EXTERN float albedo, light_min, light_max;
#endif
EXTERN int endf, syst;
EXTERN char endianfilebase[1024];
EXTERN char *endianfile;
EXTERN spherepoints sphereinfo;
EXTERN int autozip, make_demo;

