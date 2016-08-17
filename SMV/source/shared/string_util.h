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

EXTERNCPP void InitRandAB(int size);
EXTERNCPP float RandAB(int seed, float minval, float maxval);
EXTERNCPP void ToLower(char *string);
EXTERNCPP char *STRCHRR(char *strbeg, char *searchbeg, int c);
EXTERNCPP unsigned int DiffDate(char *token, char *tokenbase);
EXTERNCPP unsigned int Time2Sec(char *tokenorig);
EXTERNCPP unsigned int Date2Sec(char *tokenorig);
EXTERNCPP unsigned int Date2Sec2(char *tokenorig);
EXTERNCPP unsigned int Date2Day(char *tokenorig);
EXTERNCPP int SetLabels(flowlabels *flowlabel, char *longlabel, char *shortlabel, char *unit);
EXTERNCPP int SetLabelsIso(flowlabels *flowlabel, char *longlabel, char *shortlabel, char *unit, float *levels, int nlevels);
EXTERNCPP int ReadLabelsFaceCenter(flowlabels *flowlabel, FILE *stream);
EXTERNCPP int ReadLabelsCellCenter(flowlabels *flowlabel, FILE *stream);
EXTERNCPP int ReadLabelsTerrain(flowlabels *flowlabel, FILE *stream);
EXTERNCPP int ReadLabels(flowlabels *label, FILE *stream);
EXTERNCPP void GetProgVersion(char *PROGversion);
EXTERNCPP int MatchWild(char *pTameText, char *pWildText);
EXTERNCPP int Match(char *buffer, const char *key);
EXTERNCPP int MatchUpper(char *buffer, const char *key);
EXTERNCPP int RandInt(int min, int max);
EXTERNCPP void FParseCSV(char *buffer, float *vals, int *valids, int ncols, int *ntokens);
EXTERNCPP void ParseCSV(char *buffer, char **tokens, int *ntokens);
EXTERNCPP void StripQuotes(char *buffer);
EXTERNCPP void StripCommas(char *buffer);
EXTERNCPP int GetRowCols(FILE *stream, int *nrows, int *ncols);

EXTERNCPP char *RemoveComment(char *buffer);
EXTERNCPP void TrimBack(char *line);
EXTERNCPP void TrimCommas(char *line);
EXTERNCPP char *TrimFront(char *line);
EXTERNCPP void TrimZeros(char *line);
EXTERNCPP void TrimMZeros(char *line);
EXTERNCPP char *Strstr(char *c, char *key);
EXTERNCPP char *STRSTR(char *c, const char *key);
EXTERNCPP void ScaleString(const char *stringfrom, char *stringto, const float *scale);
EXTERNCPP void ScaleFloat2String(float floatfrom, char *stringto, const float *scale);
EXTERNCPP void Num2String(char *string, float tval);
EXTERNCPP char *TrimFrontBack(char *buffer);
EXTERNCPP int STRCMP(const char *s1, const char *s2);
EXTERNCPP char *GetChid(char *file, char *buffer);
#ifdef pp_GPU
EXTERNCPP int LogBase2(float xx);
#endif
EXTERNCPP void Array2String(float *vals, int nvals, char *string);
EXTERNCPP float FrExp10(float x, int *exp10);
EXTERNCPP void GetGitInfo(char *githash, char *gitdate);
EXTERNCPP char *GetString(char *buffer);
EXTERNCPP char *Time2TimeLabel(float time, float dt, char *timelabel);
EXTERNCPP char *RandStr(char* str, int length);
EXTERNCPP void GetBaseTitle(char *progname, char *title_base);
EXTERNCPP void GetTitle(char *progname, char *fulltitle);
EXTERNCPP void PRINTversion(char *progname);

#ifdef WIN32
STREXTERN char STRDECL(dirseparator[],"\\");
#else
STREXTERN char STRDECL(dirseparator[],"/");
#endif
#endif
