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
#include "smokeviewapi.h"
#include "translate.h"
#include "string_util.h"

// svn revision character string
char main_revision[]="$Revision$";
//  dummy change to update version to  5.6.3
//  dummy change to  force revision update

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char **argv_sv;

  initMM();
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
#else
  argv_sv=argv;
#endif
  if(argc!=0){
    int startup_flag;
    size_t len;
    char *progname;

    progname=argv_sv[0];
    smvprogdir=getprogdir(progname);
    if(texturedir==NULL){
      char *texture_buffer;
      size_t texture_len;

      texture_buffer=getenv("texturedir");
      if(texture_buffer!=NULL){
        texture_len=strlen(texture_buffer);
        NewMemory((void **)&texturedir,texture_len+1);
        strcpy(texturedir,texture_buffer);
      }
      if(texturedir==NULL&&smvprogdir!=NULL){
        texture_len=strlen(smvprogdir)+strlen("textures");
        NewMemory((void **)&texturedir,texture_len+2);
        strcpy(texturedir,smvprogdir);
        //strcat(texturedir,dirseparator);
        strcat(texturedir,"textures");
      }
    }

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
      sv_update();
    }
    else{
      exit(1);
    }
  }
  return 0;

}	 

