// $Date$ 
// $Revision$
// $Author$

// build Smokeview as a standard release unless the pp_ALPHA or pp_BETA directives are defined

#define pp_release

// uncomment the following line if using LINT to check code
// #define pp_LINT
#ifdef pp_LINT
#include "lint.h"
#endif

#ifdef pp_BETA
#define SMVVERSION "Test"
#undef pp_release
#endif

#ifdef pp_ALPHA
#define SMVVERSION "Experimental"
#undef pp_release
#endif

#ifdef pp_ALPHA2
#define pp_ALPHA
#define SMVVERSION "Experimental Full"
#undef pp_release
#endif

#ifdef pp_release
#define SMVVERSION "5.4.8"
#endif


#define _CRT_SECURE_NO_DEPRECATE

#ifdef WIN32
#ifndef X64
#define pp_memstatus
#endif
#define pp_COMPRESS
#define pp_noappend
#include "pragmas.h"

#ifndef X64
#define pp_GPU
#define pp_CULL
#endif
// #define pp_GPU_BLANK
#endif

#ifdef X64
#undef BIT64
#define BIT64
#endif

#ifdef pp_LINUX64
#undef BIT64
#define BIT64
#endif

#ifdef BIT64
#define FILE_SIZE unsigned long long
#else
#define FILE_SIZE unsigned int
#endif

#define pp_RENDER

#ifndef X64
#define pp_THREAD
#endif

#define pp_MESSAGE
#define pp_SPHERE
#define pp_ISOOUT
#define pp_DRAWISO
#define EGZ
#define USE_ZLIB
//#define NO_GLUTPOSTDISPLAY

#ifdef NO_GLUTPOSTDISPLAY
#define GLUTPOSTREDISPLAY
#else
#define GLUTPOSTREDISPLAY glutPostRedisplay();
#endif


#ifdef pp_ALPHA
#undef USE_ZLIB
#undef pp_RENDER
#undef pp_THREAD
//#define pp_LIGHT
//#define pp_SHOWLIGHT
//#define pp_COLOR
#endif

#ifdef pp_BETA
#define pp_LIGHT
#define pp_SHOOTER
#define pp_TRANSFORMxxx
#endif

#ifdef pp_CULL
#define pp_GPU
#endif

#ifndef pp_OSX
#define pp_JPEG
#endif

#undef pp_OPEN
#ifdef WIN32
#ifdef pp_BETA
#define pp_OPEN
#endif
#endif
