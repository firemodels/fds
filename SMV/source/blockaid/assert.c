// $Date: 2012-08-20 19:39:01 -0400 (Mon, 20 Aug 2012) $ 
// $Revision: 12156 $
// $Author: koverholt $

#include "options.h"
#include <stdio.h>
#include <stdlib.h>

// svn revision character string
char assert_revision[]="$Revision: 12156 $";

void _Assert(char *filename, unsigned linenumber){
  int dummy;
#ifdef _DEBUG
  float x=0.0, y;
#endif

#ifdef _DEBUG
   y=(float)1.0/(float)x;
#endif
   fflush(NULL);
   fprintf(stderr, "\nAssertion failed: %s, line %u\n",filename, linenumber);
   fflush(stderr);
   scanf("%i",&dummy);
   abort();
}

void _WAssert(char *comment, char *filename, unsigned linenumber){

  fflush(NULL);
  fprintf(stderr, "\nWarning: %s\nAssertion failed: %s, line %u\n",comment,filename, linenumber);
  fflush(stderr);
}
