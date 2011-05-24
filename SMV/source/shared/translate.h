// $Date$ 
// $Revision$
// $Author$

#ifdef IN_TRANSLATE
#define TREXTERN
#define TRDECL(var,val)  var=val
#else
#define TREXTERN extern CCC
#define TRDECL(var,val)  var
#endif

/* --------------------------  structs ------------------------------------ */

typedef struct {
  char *key, *value;
} trdata;


//************************** headers ****************************************

TREXTERN char *translate(char *string);
void init_translate(char *bindir, char *tr_name);

TREXTERN char tr_string[1024];
TREXTERN int TRDECL(tr_lang,0);
TREXTERN char *smokeview_lang;
TREXTERN trdata TRDECL(*trinfo,NULL);
TREXTERN int TRDECL(ntrinfo,0);
