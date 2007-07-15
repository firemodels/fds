
 /* ------------------ options ------------------------ */

#ifndef pp_RELEASE
#undef _DEBUG
#define _DEBUG
#endif
#ifdef WIN32
#define pp_noappend
#define pp_cvf
#endif

#define pp_ISO
#define EGZ
#define USE_ZLIB

#ifdef pp_TEST
#define pp_RLETEST
#endif
