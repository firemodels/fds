// $Date: 2012-08-20 19:39:01 -0400 (Mon, 20 Aug 2012) $ 
// $Revision: 12156 $
// $Author: koverholt $

/* ------------------ ASSERT ------------------------ */
#ifndef ASSERT_DEFINED
#define ASSERT_DEFINED

#ifdef _DEBUG

#ifdef CPP
#define ASSERT_EXTERN extern "C"
#else
#define ASSERT_EXTERN
#endif

ASSERT_EXTERN void _Assert(char *file, unsigned linenumber);
ASSERT_EXTERN void _WAssert(char *comment, char *file, unsigned linenumber);
#define WASSERT(f,g) if(f){}else{_WAssert(g,__FILE__,__LINE__);}
#define ASSERT(f) if(f){}else{_Assert(__FILE__,__LINE__);}
#define ASSERTFLAG(f) ASSERT((f)==0||(f)==1)
#else
  #define WASSERT(f,g)
  #define ASSERT(f)
  #define ASSERTFLAG(f)
#endif
#endif

