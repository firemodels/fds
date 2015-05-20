// $Date: 2012-08-20 19:39:01 -0400 (Mon, 20 Aug 2012) $ 
// $Revision: 12156 $
// $Author: koverholt $

#define pp_STANDARD

#ifdef pp_BETA
#define BLOCKAIDVERSION "Test"
#undef pp_STANDARD
#endif

#ifdef pp_STANDARD
#define BLOCKAIDVERSION "1.0.0"
#endif

#define _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_DEPRECATE
