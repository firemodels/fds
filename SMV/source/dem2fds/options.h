#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNING
#define PROGVERSION "1.0"


 /* ------------------ options ------------------------ */


int reg_path(int setget, int pathtype, char *path);

#define FILE_SIZE unsigned long long

#ifdef CPP
#define CCC "C"
#else
#define CCC
#endif

#ifdef INMAIN
#define SVEXTERN
#define SVDECL(var,val)  var=val
#else
#define SVEXTERN extern CCC
#define SVDECL(var,val)  var
#endif

#ifdef CPP
#define EXTERNCPP extern "C"
#else
#define EXTERNCPP
#endif

#ifdef X64
#define STRUCTSTAT struct __stat64
#define STAT _stat64
#else
#define STRUCTSTAT struct stat
#define STAT stat
#endif

#include "lint.h"
