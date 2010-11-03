// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#define INMAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pcat.h"
#include "MALLOC.h"

//dummy change to bump version number to 1.2.5

// svn revision character string
char main_revision[]="$Revision$";

/* ------------------ main ------------------------ */

int main(int argc, char **argv){

  char *prog=NULL, *arg;
  char *filebase=NULL;
  int i;
  int return_code;

  prog=argv[0];
  if(argc==1){
    version();
    return 1;
  }

  for(i=1;i<argc;i++){
    int lenarg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      switch(arg[1]){
      case 'h':
        usage(prog);
        return 1;
        break;
      default:
        usage(prog);
        return 1;
      }
    }
    else{
      if(filebase==NULL){
        filebase=argv[i];
      }
    }
  }
  if(filebase==NULL){
    usage(prog);
    return 1;
  }

  return 0;
}
       
/* ------------------ usage ------------------------ */

void usage(char *prog){
  char pp[2];
  char spr_version[100];
  int svn_num;

  getSPRversion(spr_version);  // get spreadpose version (ie 1.x.y)
  svn_num=getmaxrevision();    // get svn revision number

  strcpy(pp,"%");
  printf("\n");
  printf("  pcat %s(%i) - %s\n\n",spr_version,svn_num,__DATE__);
  printf("  Transposes an FDS spreadsheet so that sensor data grouped by position\n");
  printf("  (eg. thermocouple trees) may be more easily plotted\n\n");
  printf("  %s",prog);
  printf(" [-c config file] ");
  printf(" [-t time skip] csv_in\n\n");
  printf(" -c config file - specify a file containing:\n");
  printf("      line 1: list of comma delimited field names to be tranposed\n");
  printf("      line 2: list of corresponding values (usually elevations) for each field\n");
  printf(" -t skip - copy data every skip seconds (default %i s).  skip must be\n",DT_SKIP);
  printf("           larger than the time intervals in the file to be transposed.\n");
  printf(" csv_in - FDS spreadsheet file (csv_in.csv) to be transposed.  The transposed\n");
  printf("          file is output in csv_in_tran.csv\n");
}

