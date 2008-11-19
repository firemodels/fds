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
#define SMVVERSION "5.2.7"
#endif


#define _CRT_SECURE_NO_DEPRECATE
#ifdef WIN32
#define pp_memstatus
#define pp_COMPRESS
#define pp_noappend
#include "pragmas.h"
#define pp_GPU
#define pp_CULL
// #define pp_GPU_BLANK
#endif

#define pp_SPHERE
#define pp_ISOOUT
#define pp_DRAWISO
#define EGZ
#define USE_ZLIB
//#define NO_GLUTPOSTDISPLAY

#ifdef pp_ALPHA2
#define pp_TOUR
#endif

#ifdef NO_GLUTPOSTDISPLAY
#define GLUTPOSTREDISPLAY
#else
#define GLUTPOSTREDISPLAY glutPostRedisplay();
#endif


#ifdef pp_ALPHA
#define pp_LIGHT
#define pp_SHOWLIGHT
#define pp_COLOR
#endif

#ifdef pp_CULL
#define pp_GPU
#endif

#undef pp_OPEN
#ifdef WIN32
#ifdef pp_BETA
#define pp_OPEN
#endif
#endif








