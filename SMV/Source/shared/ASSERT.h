
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

