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

typedef struct {
  char **keyword_list, **val_list;
  int nkeywordlist;
} keyvaldata;

typedef struct _fdsdata {
  char *line,*linecopy,*line_before, *line_after;
  float xb[6];
  int ibeg, iend;
  int type,is_obst;
  struct _blockaiddata *blockaid;
  struct _fdsdata *prev, *next;
} fdsdata;

typedef struct _blockaiddata {
  char *id;
  char **keyword_list, **val_list;
  int nkeywordlist;
  float xyz0[4];
  float bb_min[3], bb_dxyz[3], bb_max[3];
  int bb_box_defined;
  struct _fdsdata *first_line, *last_line;
  struct _fdsdata f_line, l_line;
  struct _blockaiddata *prev, *next;
} blockaiddata;

int getrevision(char *svn);
void getASMversion(char *SMZversion);
int getmaxrevision(void);
void version(void);
char *trim_front(char *line);
void trim(char *line);
void usage(void);
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
void rotatexy(float *dx, float *dy, float rotate, float *dxy);
void get_boundbox(blockaiddata *assem, int recurse_level);
void init_bb(void);
void reorder(float *xy);
void get_keywords(blockaiddata *blockaidi, char *buffer);
void subst_keys(char *line_after,int recurse_level);
char *get_val(char *key, int recurse_level);

EXTERN blockaiddata *blockaid_first, *blockaid_last, ba_first, ba_last;
EXTERN blockaiddata **assemblylist;
EXTERN keyvaldata keyvalstack[MAXRECURSE];
EXTERN nkeyvalstack;
EXTERN int nassembly;
EXTERN float *offset_rotate;
EXTERN int nblockaid;
