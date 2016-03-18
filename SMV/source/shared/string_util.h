#ifndef STRING_UTIL_H_DEFINED
#define STRING_UTIL_H_DEFINED

#ifdef IN_STRING_UTIL
#define STREXTERN
#define STRDECL(var,val)  var=val
#else
#define STREXTERN extern CCC
#define STRDECL(var,val)  var
#endif

#ifdef __MINGW32__
#include <stdio.h>
#include "options.h"
#endif

#define MATCH 1
#define NOTMATCH 0

/* --------------------------  flowlabels ------------------------------------ */

typedef struct {
  char *longlabel, *shortlabel, *unit;
} flowlabels;

EXTERNCPP void init_rand_ab(int size);
EXTERNCPP float rand_ab(int seed, float minval, float maxval);
EXTERNCPP void to_lower(char *string);
EXTERNCPP char *STRCHRR(char *strbeg, char *searchbeg, int c);
EXTERNCPP unsigned int diffdate(char *token, char *tokenbase);
EXTERNCPP unsigned int time2sec(char *tokenorig);
EXTERNCPP unsigned int date2sec(char *tokenorig);
EXTERNCPP unsigned int date2sec2(char *tokenorig);
EXTERNCPP unsigned int date2day(char *tokenorig);
EXTERNCPP int setlabels(flowlabels *flowlabel, char *longlabel, char *shortlabel, char *unit);
EXTERNCPP int setlabels_iso(flowlabels *flowlabel, char *longlabel, char *shortlabel, char *unit, float *levels, int nlevels);
EXTERNCPP int readlabels_facecenter(flowlabels *flowlabel, FILE *stream);
EXTERNCPP int readlabels_cellcenter(flowlabels *flowlabel, FILE *stream);
EXTERNCPP int readlabels_terrain(flowlabels *flowlabel, FILE *stream);
EXTERNCPP int readlabels(flowlabels *label, FILE *stream);
EXTERNCPP void getPROGversion(char *PROGversion);
EXTERNCPP int match_wild(char *pTameText, char *pWildText);
EXTERNCPP int match(char *buffer, const char *key);
EXTERNCPP int match_upper(char *buffer, const char *key);
EXTERNCPP int randint(int min, int max);
EXTERNCPP void fparsecsv(char *buffer, float *vals, int *valids, int ncols, int *ntokens);
EXTERNCPP void parsecsv(char *buffer, char **tokens, int ncols, int *ntokens);
EXTERNCPP void stripquotes(char *buffer);
EXTERNCPP void stripcommas(char *buffer);
EXTERNCPP int getrowcols(FILE *stream, int *nrows, int *ncols);

EXTERNCPP char *remove_comment(char *buffer);
EXTERNCPP void trim_back(char *line);
EXTERNCPP void trim_commas(char *line);
EXTERNCPP char *trim_front(char *line);
EXTERNCPP void trimzeros(char *line);
EXTERNCPP void trimmzeros(char *line);
EXTERNCPP char *Strstr(char *c, char *key);
EXTERNCPP char *STRSTR(char *c, const char *key);
EXTERNCPP void scalestring(const char *stringfrom, char *stringto, const float *scale, float range);
EXTERNCPP void scalefloat2string(float floatfrom, char *stringto, const float *scale, float range);
EXTERNCPP void num2string(char *string, float tval,float range);
EXTERNCPP char *trim_frontback(char *buffer);
EXTERNCPP int STRCMP(const char *s1, const char *s2);
EXTERNCPP char *get_chid(char *file, char *buffer);
#ifdef pp_GPU
EXTERNCPP int log_base2(float xx);
#endif
EXTERNCPP void array2string(float *vals, int nvals, char *string);
EXTERNCPP float frexp10(float x, int *exp10);
EXTERNCPP void getGitInfo(char *githash, char *gitdate);
EXTERNCPP char *getstring(char *buffer);
EXTERNCPP char *time2timelabel(float time, float dt, char *timelabel);
EXTERNCPP char *randstr(char* str, int length);
EXTERNCPP void getBaseTitle(char *progname, char *title_base);
EXTERNCPP void getTitle(char *progname, char *fulltitle);
EXTERNCPP void version(char *progname);

#ifdef WIN32
STREXTERN char STRDECL(dirseparator[],"\\");
#else
STREXTERN char STRDECL(dirseparator[],"/");
#endif
#endif
