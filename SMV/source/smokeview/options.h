#ifndef OPTIONS_H_DEFINED
#define OPTIONS_H_DEFINED
// build Smokeview as a standard release unless the pp_BETA directive is defined

#define pp_release

//*** uncomment the following two lines to force all versions to be beta
//#undef pp_BETA
//#define pp_BETA

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVV  define smokeview title VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
#ifdef pp_BETA
#define PROGVERSION "Test"
#undef pp_release
#endif

// comment the following line when building an unofficial release
#define pp_OFFICIAL_RELEASE

#ifdef pp_release
#ifdef pp_OFFICIAL_RELEASE
#define PROGVERSION "6.3.6"
#else
#define PROGVERSION "Unofficial release"
#endif
#endif

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVV  turn on options available on all platforms VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

#define pp_GPU
#define pp_ffmpeg

#ifdef pp_GPU
#define pp_CULL
#define pp_GPUTHROTTLE
#endif

#ifdef VS2015
#define HAVE_SNPRINTF
#define HAVE_STRUCT_TIMESPEC
#endif
#define pp_DRAWISO
#define pp_LANG
#define pp_DEG
#define pp_SLICEDUP
#define pp_SLICECOLORDEFER

#define pp_THREAD
#define pp_THREADIBLANK // test iblank computation in background.
#ifdef pp_THREADIBLANK  // if pp_THREADIBLANK is set then pp_THREAD also has to be set
#define pp_THREAD
#endif

#define _CRT_SECURE_NO_DEPRECATE

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVV  turn on options that are being tested VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

#ifdef pp_BETA
#define pp_PILOT
//#define pp_WINDROSE  2d histogram, variation of pilot data
#define pp_GEOMTEST
#define pp_HAZARD
//#define pp_GPUDEPTH
#define pp_MEMPRINT
#endif

// for debugging, set particle values to 100*parti->seq_id + small random number
//#define pp_PARTTEST

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#ifdef _DEBUG  // comment the following line when debugging REALLY large cases (to avoid memory checks)
#define pp_MEMDEBUG
#endif

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVV  turn on windows only options VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

#ifdef WIN32
#define pp_memstatus
#define pp_COMPRESS
#define pp_noappend
#include "pragmas.h"
#endif
//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVV  used to access fortran routines from C VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

#ifndef _F
#ifdef pp_noappend
#define _F(name) name
#else
#define _F(name) name ## _
#endif
#endif

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


// VVVVVVVVVVVVVVVVVVVVVVVVV  set platform defines VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

#define FILE_SIZE unsigned long long

#ifdef X64
#define STRUCTSTAT struct __stat64
#define STAT _stat64
#else
#define STRUCTSTAT struct stat
#define STAT stat
#endif

#define LINT long int
#ifdef X64
#undef LINT
#ifdef WIN32
#define LINT __int64
#else
#define LINT long long int
#endif
#endif

//*** turn off features on the Mac that don't work there

#ifdef pp_OSX
#undef pp_LANG
#undef pp_DEG
#endif

// VVVVVVVVVVVVVVVVVVVVVVVVV  set defines used by various headers VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
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

#ifdef pp_OSX
#define GLUT_H <GLUT/glut.h>
#else
#define GLUT_H <GL/glut.h>
#endif


#include "lint.h"
#endif

