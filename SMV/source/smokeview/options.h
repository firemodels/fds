// $Date: 2015-04-28 14:58:18 -0400 (Tue, 28 Apr 2015) $ 
// $Revision: 22541 $
// $Author: gforney $

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
#define PROGVERSION "6.2.3"
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

#define pp_DRAWISO
#define pp_THREAD
#define pp_LANG
#define pp_DEG

#define _CRT_SECURE_NO_DEPRECATE

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

//VVVVVVVVVVVVVVVVVVVVVVVVVVVVV  turn on options that are being tested VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

#ifdef pp_BETA
#define pp_SETTIME
#define pp_PILOT
#define pp_GEOMTEST
#define pp_ZONEVENT
#define pp_HAZARD
//#define pp_GPUDEPTH
#define pp_MEMPRINT
#endif

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

#ifdef pp_LINUX64
#define pp_LINUX
#endif

#ifdef pp_OSX64
#define pp_OSX
#endif

//*** turn on BIT64 if compiled on a 64 bit platform

#ifdef X64
#undef BIT64
#define BIT64
#endif

#ifdef pp_LINUX64
#undef BIT64
#define BIT64
#endif

#ifdef pp_OSX64
#undef pp_OSX
#define pp_OSX
#undef BIT64
#define BIT64
#endif

#ifdef BIT64
#define FILE_SIZE unsigned long long
#else
#define FILE_SIZE unsigned int
#endif

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

