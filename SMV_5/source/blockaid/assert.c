// $Date: 2007-10-08 21:39:01 -0400 (Mon, 08 Oct 2007) $ 
// $Revision: 812 $
// $Author: gforney $

#include "options.h"
#include <stdio.h>
#include <stdlib.h>

// svn revision character string
char assert_revision[]="$Revision: 824 $";

void _Assert(char *filename, unsigned linenumber){
  int dummy;
#ifdef _DEBUG
  float x=0.0, y;
#endif

#ifdef _DEBUG
   y=1.0/x;
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
