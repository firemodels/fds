// $Date: 2015-05-20 16:50:12 -0400 (Wed, 20 May 2015) $ 
// $Revision: 22691 $
// $Author: gforney $
#ifdef INMAIN
#define EXTERN
#else
#define EXTERN extern
#endif

void version(void);
int getrevision(char *svn);
char *hostlistfile;
char *host;
#ifdef pp_LINUX
char  pid[20];
#endif
#ifdef pp_OSX
char  pid[20];
#endif
int nhostinfo;

typedef struct {
  char *hostname;
  int ncores;
} hostdata;

hostdata *hostinfo;
