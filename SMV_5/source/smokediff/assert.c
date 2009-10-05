// $Date: 2007-10-09 15:35:36 -0400 (Tue, 09 Oct 2007) $ 
// $Revision: 824 $
// $Author: gforney $

#include "options.h"
#include <stdio.h>
#include <stdlib.h>

// svn revision character string
char assert_revision[]="$Revision: 824 $";


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
