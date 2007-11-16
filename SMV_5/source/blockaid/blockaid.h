// $Date: 2007-10-09 15:35:36 -0400 (Tue, 09 Oct 2007) $ 
// $Revision: 824 $
// $Author: gforney $

#define FFALSE 0
#define TTRUE 1

#ifdef INMAIN
#define EXTERN
#else
#define EXTERN extern
#endif

typedef struct _fdsdata {
  char *line,*line_before, *line_after;
  float xb[6];
  int ibeg, iend;
  int type;
  struct _fdsdata *prev, *next;
} fdsdata;

typedef struct _blockaiddata {
  char *id;
  float orig[3];
  int in_use;
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
void usage(char *prog);
int getfileinfo(char *filename, char *source_dir, int *filesize);
int readfds(char *fdsfile);
int match(const char *buffer, const char *key, unsigned int lenkey);
int get_fds_line(FILE *stream, char *fdsbuffer, unsigned int len_fdsbuffer);
void expand_assembly(char *buffer, int first_time);
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
void rotatexy(float *dx, float *dy, float *orig, float rotate);

EXTERN blockaiddata *blockaidinfo, *blockaid_first, *blockaid_last, ba_first, ba_last;
EXTERN int nblockaid;