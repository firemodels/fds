// $Date$ 
// $Revision$
// $Author$

#ifdef INMAIN
#define EXTERN
#else
#define EXTERN extern
#endif

#define VERSION "1.0"
EXTERN int getmaxrevision(void);

extern char assert_revision[];
extern char dmalloc_revision[];
extern char main_revision[];
