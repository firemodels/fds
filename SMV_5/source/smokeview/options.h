// $Date$ 
// $Revision$
// $Author$

// build Smokeview as a standard release unless the pp_ALPHA or pp_BETA directives are defined

#define pp_STANDARD

// uncomment the following line if using LINT to check code
// #define pp_LINT
#ifdef pp_LINT
#include "lint.h"
#endif

#ifdef pp_BETA
#define SMVVERSION "Test"
#undef pp_STANDARD
#endif

#ifdef pp_ALPHA
#define SMVVERSION "Experimental"
#undef pp_STANDARD
#endif

#ifdef pp_ALPHA2
#define pp_ALPHA
#define SMVVERSION "Experimental Full"
#undef pp_STANDARD
#endif

#ifdef pp_STANDARD
#define SMVVERSION "5.0.7"
#endif


#define _CRT_SECURE_NO_DEPRECATE
#ifdef WIN32
#define pp_memstatus
#define pp_COMPRESS
#define pp_noappend
#include "pragmas.h"
#endif

#ifndef _DEBUG
#define pp_RELEASE
#endif

#define pp_DRAWISO
#define EGZ
#define USE_ZLIB
#define pp_THREADS2

#ifdef pp_ALPHA2
#define pp_TOUR
#define pp_LIGHT
#endif
#ifdef pp_ALPHA
#define pp_COMPRESS_AUTOLOADED
#define pp_AVATAR
#define pp_SHOWLIGHT
#define pp_COLOR
#define pp_GPU
#define pp_MEM2
#define pp_CULL
#endif

#undef pp_OPEN
#ifdef WIN32
#ifdef _DEBUG
#define pp_OPEN
#endif
#ifdef pp_MEM2
#define pp_OPEN
#endif
#endif








