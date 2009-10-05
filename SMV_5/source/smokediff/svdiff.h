// $Date: 2009-07-09 11:08:46 -0400 (Thu, 09 Jul 2009) $ 
// $Revision: 4355 $
// $Author: gforney $
#ifdef INMAIN
#define EXTERN
#else
#define EXTERN extern
#endif


/* --------------------------  flowlabels ------------------------------------ */

typedef struct {
  char *longlabel, *shortlabel, *unit;
} flowlabels;

typedef struct {
  int ibar, jbar, kbar;
  float xbar0, xbar, ybar0, ybar, zbar0, zbar;
  float dx, dy, dz;
  float *xplt, *yplt, *zplt;
} mesh;

typedef struct _slice {
  char *file;
  int is1, is2, js1, js2, ks1, ks2;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  int filesize;
  int version;
  int volslice;
  struct _slice *slice2;
  char keyword[255];
  int slicetype;
  mesh *slicemesh;
  flowlabels label;
} slice;

#ifdef WIN32
#define STDCALL extern void _stdcall
#else
#define STDCALL extern void
#endif

#ifdef X64
#define STRUCTSTAT struct __stat64
#define STAT _stat64
#else
#define STRUCTSTAT struct stat
#define STAT stat
#endif

#ifdef X64
#undef BIT64
#define BIT64
#endif

#ifdef pp_LINUX64
#undef BIT64
#define BIT64
#endif

#ifdef BIT64
#define FILE_SIZE unsigned long long
#else
#define FILE_SIZE unsigned int
#endif

typedef struct {
  slice *sliceinfo;
  mesh *meshinfo;
  char *dir;
  int nmeshes;
  int nslice_files;
} casedata;

void getSMDiffversion(char *SMDiffversion);
int getmaxrevision(void);
int imax(int a, int b);
int getrevision(char *svn);
void version(void);
void usage(char *prog);
int getfileinfo(char *filename, char *source_dir, int *filesize);
int match(const char *buffer, const char *key, unsigned int lenkey);
void trim(char *line);
char *trim_front(char *line);
int readlabels(flowlabels *flowlabel, FILE *stream);
char *setdir(char *argdir);
int readsmv(FILE *streamsmv, FILE *stream_out, casedata *smvcase);
void setup_slice(FILE *stream_out);
slice *getslice(slice *slicein, casedata *case2);
void diff_slices(void);
void fullfile(char *fileout, char *dir, char *file);

#ifdef pp_noappend
#define FORTgetsliceparms getsliceparms
#else
#define FORTgetsliceparms getsliceparms_
#endif
STDCALL FORTgetsliceparms(char *file,int *endian,
                          int *is1,int *is2,int *js1,int *js2,int *ks1, int *ks2,int *slice3d, int *error,FILE_SIZE lenfile);

EXTERN char dirseparator[3];
EXTERN casedata caseinfo[2];
EXTERN char *sourcedir1, *sourcedir2, *destdir;

