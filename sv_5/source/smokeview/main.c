#include "options.h"
#define INMAIN
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "MALLOC.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"
#include "smokeviewapi.h"

// svn revision character string
char main_revision[]="$Revision$";


/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char **argv_sv;
#ifdef WIN32
  char *smv_file;
#endif

#ifdef pp_MEM2
  initMM();
#endif
  initvars0();
  initcolors();
  initvars1();
  if(argc==1){
    version();
  }
#ifdef WIN32
  copy_args(&argc, argv, &argv_sv);
  if(argc==1){
    exit(0);
  }
  smv_file=argv_sv[1];
#else
  argv_sv=argv;
#endif
  if(argc!=0){
    int startup_flag;
    size_t len;
    char *progname;

    progname=argv_sv[0];
    len=strlen(progname);
    if(len>2){
      NewMemory((void **)&smvprogdir,len+1);
    }
    else{
      smvprogdir=NULL;
    }
    if(smvprogdir!=NULL){
      strcpy(smvprogdir,progname);
      getdir(smvprogdir);
    }
    get_smokezippath(smvprogdir,&smokezippath);

    if(smokezippath!=NULL)printf(" Smokezip file: %s found\n",smokezippath);

    CheckMemory;
    Args(argc, argv_sv);
    sv_startup_c(argc,argv_sv);
    CheckMemory;
    startup_flag=initcase_c(argc,argv_sv);
    if(startup_flag==0){
      sv_update();
    }
    else{
      exit(1);
    }
  }
  return 0;

}	 

