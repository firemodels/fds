// $Date: 2009-02-25 10:22:50 -0500 (Wed, 25 Feb 2009) $ 
// $Revision: 3416 $
// $Author: gforney $

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
char assert_revision[]="$Revision: 3416 $";

void _Assert(char *filename, unsigned linenumber){
  int dummy;
#ifdef _DEBUG
  float x=0.0, y;
#endif

#ifdef _DEBUG
   y=1.0/x;
   printf("y=%f\n",y);
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
