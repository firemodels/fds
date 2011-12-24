// $Date$ 
// $Revision$
// $Author$

#define INMAIN
#include "options.h"
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
#include "translate.h"

// svn revision character string
char main_revision[]="$Revision$";
//  dummy change to update version to  5.6.3
//  dummy change to force revision update

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char **argv_sv;
  int startup_flag;
  size_t len;
  char *progname;

  //listdir(".");
  initMALLOC();
  initcolors();
  initvars();
  if(argc==1){
    version();
  }
#ifdef WIN32
  copy_args(&argc, argv, &argv_sv);
  if(argc==1){
    exit(0);
  }
#else
  argv_sv=argv;
#endif

  if(argc==0)return 0;

  progname=argv_sv[0];
  smvprogdir=getprogdir(progname);
  init_texturedir();
  get_smokezippath(smvprogdir,&smokezippath);

  CheckMemory;
  Args(argc, argv_sv);
  version();
  printf("\n");
  if(smokezippath!=NULL)printf("Smokezip file: %s found\n",smokezippath);
  sv_startup_c(argc,argv_sv);
  init_translate(smvprogdir,tr_name);
  CheckMemory;
  startup_flag=initcase_c(argc,argv_sv);
  if(update_bounds==1){
    Update_All_Patch_Bounds();
#ifdef pp_THREAD
    pthread_join(update_all_patch_bounds_id,NULL);
#endif
    return 0;
  }
  if(startup_flag==0){
    glutMainLoop();
  }
  else{
    exit(1);
  }
  return 0;
}	 

