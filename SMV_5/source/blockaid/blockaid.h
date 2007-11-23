// $Date$ 
// $Revision$
// $Author$

#define FFALSE 0
#define TTRUE 1
#define MAXRECURSE 100
#define MAXPOS 100000000.0
#define MINPOS -MAXPOS

#ifdef INMAIN
#define EXTERN
#else
#define EXTERN extern
#endif

typedef struct _fdsdata {
  char *line,*linecopy,*line_before, *line_after;
  float xb[6];
  int ibeg, iend;
  int type,is_obst;
  struct _fdsdata *prev, *next;
} fdsdata;

typedef struct _blockaiddata {
  char *id;
  float orig[4],xyzmax[3], dxy[3];
  int in_use;
  float bb_min[3], bb_dxyz[3], bb_max[3];
  int bb_box_defined;
  struct _fdsdata *first_line, *last_line;
  struct _fdsdata f_line, l_line;
  struct _blockaiddata *prev, *next;
  struct _blockaiddata **assemblylist[MAXRECURSE];
} blockaiddata;

int getrevision(char *svn);
void getASMversion(char *SMZversion);
int getmaxrevision(void);
void version(void);
char *trim_front(char *line);
void trim(char *line);
void usage(char *prog);
int getfileinfo(char *filename, char *source_dir, int *filesize);
int readfds(char *fdsfile);
int match(const char *buffer, const char *key, unsigned int lenkey);
int get_fds_line(FILE *stream, char *fdsbuffer, unsigned int len_fdsbuffer);
void expand_assembly(char *buffer, int recurse_level);
blockaiddata *create_assembly(char *buffer);
void init_assemdata(char *id, float *orig, blockaiddata *prev, blockaiddata *next);
void update_assembly(blockaiddata *assembly,char *buffer);
void remove_assembly(blockaiddata *assembly);
char *getkeyid(char *source, const char *key);
int get_irvals(char *line, char *key, int nvals, float *ivals, float *rvals, int *ibeg, int *iend);
void startup(void);
blockaiddata *get_assembly(char *id);
void trimzeros(char *line);
void trimmzeros(char *line);
void rotatexy(float *dx, float *dy, float *orig, float rotate, float *dxy);
void get_boundbox(blockaiddata *assem, int recurse_level);
void init_bb(void);
void reorder(float *xy);

EXTERN blockaiddata *blockaidinfo, *blockaid_first, *blockaid_last, ba_first, ba_last;
EXTERN blockaiddata **assemblylist, **assembly_sorted_list;
EXTERN int nassembly;
EXTERN float *offset_rotate;
EXTERN int nblockaid;