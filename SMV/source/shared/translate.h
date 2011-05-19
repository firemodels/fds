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


//************************** headers ****************************************

char *translate(char *string);

TREXTERN char tr_string[1024];
TREXTERN int TRDECL(tr_lang,0);
