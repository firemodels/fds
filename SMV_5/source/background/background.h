// $Date: 2010-04-11 20:48:06 -0400 (Sun, 11 Apr 2010) $ 
// $Revision: 6048 $
// $Author: gforney $
#ifdef INMAIN
#define EXTERN
#else
#define EXTERN extern
#endif

void getPROGversion(char *PROGversion);
int getmaxrevision(void);
void version(void);
int getrevision(char *svn);





