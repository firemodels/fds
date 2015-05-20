#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNING
#define PROGVERSION "1.1"


 /* ------------------ options ------------------------ */


#pragma warning (disable:1419)		/* external declaration in primary source file */

#define REG_GET 0
#define REG_SET 1
#define REG_USER_PATH 2
#define REG_SYSTEM_PATH 3

int reg_path(int setget, int pathtype, char *path);

//*** turn on BIT64 if compiled on a 64 bit platform

#ifdef X64
#undef BIT64
#define BIT64
#endif

#ifdef BIT64
#define FILE_SIZE unsigned long long
#else
#define FILE_SIZE unsigned int
#endif

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

