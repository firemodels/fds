// $Date: 2008-10-30 14:40:35 -0400 (Thu, 30 Oct 2008) $ 
// $Revision: 2576 $
// $Author: gforney $

#include "options.h"
#define INMAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "transpose.h"
#include "MALLOC.h"

//dummy change to bump version number to 1.2.5

// svn revision character string
char main_revision[]="$Revision: 2576 $";

/* ------------------ main ------------------------ */

int main(int argc, char **argv){

  char *prog=NULL, *arg;
  char *filebase=NULL;
  char *fileoutbase=NULL;
  int i;
  int return_code;

  dt_skip=DT_SKIP;
  key_label=NULL;
  key_unit=NULL;
  strcpy(config_file,"spreadtran.cfg");

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
      case 't':
        if(i+1<argc&&argv[i]!=NULL){
          dt_skip=DT_SKIP;
          sscanf(argv[i+1],"%i",&dt_skip);
          if(dt_skip<=0)dt_skip=1;
          i++;
        }
        break;
      case 'c':
        if(i+1<argc&&argv[i]!=NULL){
          strcpy(config_file,argv[i+1]);
          i++;
        }
        break;
      case 'o':
        if(i+1<argc&&argv[i]!=NULL){
          fileoutbase=argv[i+1];
          i++;
        }
        break;
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

  strcpy(csv_in,filebase);
  strcat(csv_in,".csv");
  strcpy(csv_out,filebase);
  if(fileoutbase!=NULL){
    strcat(csv_out,"_");
    strcat(csv_out,fileoutbase);
  }
  strcat(csv_out,"_tran.csv");

  return_code = getfileinfo(csv_in,NULL,NULL);
  if(return_code!=0){
    printf("*** warning: The file %s could not be opened\n",csv_in);
    return 1;
  }
  {
    FILE *stream;

    stream=fopen(csv_out,"w");
    if(stream==NULL){
      printf("*** warning: The file %s could not be opened for output.\n",csv_out);
      return 1;
    }
    fclose(stream);
  }

  convert_csv(csv_in,csv_out);
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
  printf("  transpose %s(%i) - %s\n\n",spr_version,svn_num,__DATE__);
  printf("  Transposes an FDS spreadsheet so that sensor data grouped by position\n");
  printf("  (eg. thermocouple trees) may be more easily plotted\n\n");
  printf("  %s",prog);
  printf(" [-c config file] ");
  printf(" [-t time skip] csv_in\n\n");
  printf(" -c config file - specify a file containing:\n");
  printf("      line 1: list of comma delimited field names to be tranposed\n");
  printf("      line 2: list of corresponding values (usually elevations) for each field\n");
  printf(" -o file mod - add string used modify output file name\n");
  printf(" -t skip - copy data every skip seconds (default %i s).  skip must be\n",DT_SKIP);
  printf("           larger than the time intervals in the file to be transposed.\n");
  printf(" csv_in - FDS spreadsheet file (csv_in.csv) to be transposed.  The transposed\n");
  printf("          file is output in csv_in_tran.csv\n");
}

