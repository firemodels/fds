#ifndef TRANSLATE_H_DEFINED
#define TRANSLATE_H_DEFINED

#ifdef IN_TRANSLATE
#define TREXTERN
#define TRDECL(var,val)  var=val
#else
#define TREXTERN extern CCC
#define TRDECL(var,val)  var
#endif

#define _(String) Translate((char *)String,1)
#define _d(String) Translate((char *)String,0)
// #define _(String) (String)

/* --------------------------  structs ------------------------------------ */

typedef struct {
  char *key, *value;
} trdata;

//************************** headers ****************************************

TREXTERN int CompareTrdata( const void *arg1, const void *arg2 );
TREXTERN char *Translate(char *string, int option);
TREXTERN void InitTranslate(char *bindir, char *tr_name);
TREXTERN int ParseLang(char *file, trdata **trinfoptr, int *ntrinfoptr);

#define TR_STRING_MAX_LENGTH 1024

TREXTERN char tr_string[TR_STRING_MAX_LENGTH];
TREXTERN char tr_string_before[TR_STRING_MAX_LENGTH];
TREXTERN char tr_string_in[TR_STRING_MAX_LENGTH];
TREXTERN char tr_string_after[TR_STRING_MAX_LENGTH];
TREXTERN int TRDECL(tr_otherlang,0);
TREXTERN char TRDECL(*smokeview_lang,NULL);
TREXTERN trdata TRDECL(*trinfo,NULL);
TREXTERN int TRDECL(ntrinfo,0);
#endif

