// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>
#include <stdlib.h>

// svn revision character string
char assert_revision[]="$Revision$";


void _Assert(char *filename, unsigned linenumber){
  int dummy;

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
