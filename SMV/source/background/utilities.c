// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include "background.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "svn_revision.h"
#include "datadefs.h"

// svn revision character string
char utilities_revision[]="$Revision$";

/* ------------------ version ------------------------ */

void version(void){
    char smv_version[100];
    int svn_num;

    getPROGversion(smv_version);  // get Smokeview version (ie 5.x.z)
    svn_num=getmaxrevision();    // get svn revision number
    printf("\n");
    printf("background\n\n");
    printf("Version: %s\n",smv_version);
    printf("SVN Revision Number: %i\n",svn_num);
    printf("Compile Date: %s\n",__DATE__);
}

/* ------------------ getmaxrev ------------------------ */

#define MAXREV(cval) max_revision=MAX(getrevision(cval),max_revision)
int getmaxrevision(void){
  int max_revision=0;

  MAXREV(main_revision);
  MAXREV(utilities_revision);
  return max_revision;
}

/* ------------------ getPROGversion ------------------------ */

void getPROGversion(char *PROGversion){
  strcpy(PROGversion,PROGVERSION);
}
