// $Date: 2013-01-15 14:08:59 -0500 (Tue, 15 Jan 2013) $ 
// $Revision: 14445 $
// $Author: gforney $

/* ------------------ ASSERT ------------------------ */
#ifndef ASSERT_H_DEFINED
#define ASSERT_H_DEFINED

#ifdef _DEBUG

#ifdef CPP
#define ASSERT_EXTERN extern "C"
#else
#define ASSERT_EXTERN
#endif

ASSERT_EXTERN void _Assert(char *file, unsigned linenumber);
ASSERT_EXTERN void _WAssert(char *comment, char *file, unsigned linenumber);
#define ASSERT(f) if(f){}else{_Assert(__FILE__,__LINE__);}
#else
  #define ASSERT(f)
#endif
#endif

