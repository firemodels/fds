// $Date$ 
// $Revision$
// $Author$

#define pp_STANDARD

#ifdef pp_BETA
#define BLOCKAIDVERSION "Test"
#undef pp_STANDARD
#endif

#ifdef pp_STANDARD
#define BLOCKAIDVERSION "1.0.0"
#endif

#ifndef _DEBUG
#define pp_RELEASE
#endif

#ifndef pp_RELEASE
#undef _DEBUG
#define _DEBUG
#endif
#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE
