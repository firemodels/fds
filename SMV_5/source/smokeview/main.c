// $Date$ 
// $Revision$
// $Author$

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
//   dummy change to update version to 5.5.3
//   dummy change to force revision update

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
    len=strlen(progname);
    if(len>2){
      smvprogdir=getdir(progname);
    }
    else{
      smvprogdir=NULL;
    }
    if(smvprogdir==NULL||strlen(smvprogdir)==0){
      char *temp_smvbindir=NULL, SMVBINDIR[1024];

      strcpy(SMVBINDIR,"");
      FREEMEMORY(smvprogdir);
      temp_smvbindir=getenv("SMVBINDIR");
      if(temp_smvbindir==NULL)temp_smvbindir=getenv("smvbindir");
      if(temp_smvbindir!=NULL)strcpy(SMVBINDIR,temp_smvbindir);
#ifdef WIN32
      if(strlen(SMVBINDIR)==0){
        strcpy(SMVBINDIR,"c:\\program files\\fds\\fds5\\bin\\");
      }
#endif
      if(strlen(SMVBINDIR)>0){
        NewMemory((void **)&smvprogdir,strlen(SMVBINDIR)+2);
        strcpy(smvprogdir,SMVBINDIR);
      }
    }
    if(smvprogdir!=NULL){
      len=strlen(smvprogdir);
      if(smvprogdir[len-1]!=dirseparator[0]){
        strcat(smvprogdir,dirseparator);
      }
    }
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
//    if(texturedir!=NULL){
//      printf("Texture directory: %s\n",texturedir);
//    }
    version();
    printf("\n");
    if(smokezippath!=NULL)printf("Smokezip file: %s found\n",smokezippath);
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

