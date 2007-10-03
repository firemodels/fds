
 /* ------------------ options ------------------------ */
#ifdef WIN32
#pragma warning (disable:4996)		
#pragma warning (disable:4701)		
#pragma warning (disable:4310)		
#pragma warning (disable:4127)		
#pragma warning (disable:4267)		
#pragma warning (disable:4244)		
#endif



#ifndef pp_RELEASE
#undef _DEBUG
#define _DEBUG
#endif
#ifdef WIN32
#define pp_noappend
#define pp_cvf
#endif

#define EGZ
#define USE_ZLIB

#ifdef pp_TEST
#define pp_RLETEST
#endif
