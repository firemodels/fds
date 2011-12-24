// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#ifdef WIN32
#include <Winsock2.h>
#endif
#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <time.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include "flowfiles.h"
#include "MALLOC.h"
#include "smokeviewvars.h"

// svn revision character string
char sv_api_revision[]="$Revision$";

/* ------------------ copy_args ------------------------ */

#ifdef WIN32
void copy_args(int *argc, char **aargv, char ***argv_sv){
  char *filename=NULL;
  char **argv=NULL;
  int filelength=1024,openfile;
  int i;

  if(NewMemory((void **)&argv,(*argc+1)*sizeof(char **))!=0){
    *argv_sv=argv;
    for(i=0;i<*argc;i++){
      argv[i]=aargv[i];
    }
    if(*argc==1){
      if(NewMemory((void **)&filename,(unsigned int)(filelength+1))!=0){
        openfile=0;
#ifdef WIN32
        OpenSMVFile(filename,filelength,&openfile);
#endif
        if(openfile==1&&ResizeMemory((void **)&filename,strlen(filename)+1)!=0){
          *argc=2;
          argv[1]=filename;
        }
        else{
          FREEMEMORY(filename);
        }
      }
    }
  }
  else{
    *argc=0;
  }
}
#endif

/* ------------------ get_smokezippath ------------------------ */

void get_smokezippath(char *progdir, char **zippath){
  size_t len;
  STRUCTSTAT statbuffer;

  if(progdir!=NULL){
    len=strlen(progdir);
    NewMemory((void **)zippath,len+20);
    strcpy(*zippath,progdir);
#ifdef WIN32
#ifdef X64
    strcat(*zippath,"smokezip_win_64.exe");
#else
    strcat(*zippath,"smokezip.exe");
#endif
#else
    strcat(*zippath,"smokezip");
#endif
  }
  if(STAT(*zippath,&statbuffer)==0)return;
  FREEMEMORY(*zippath);
#ifdef WIN32
  NewMemory((void **)zippath,28);
  strcpy(*zippath,"c:\\nist\\fds\\smokezip.exe");
  if(STAT(*zippath,&statbuffer)==0)return;
  FREEMEMORY(*zippath);
#else
  NewMemory((void **)zippath,28);
  strcpy(*zippath,"~");
  strcat(*zippath,dirseparator);
  strcat(*zippath,"bin");
  strcat(*zippath,dirseparator);
  strcat(*zippath,"smokezip");
  if(STAT(*zippath,&statbuffer)==0)return;
  FREEMEMORY(*zippath);
  NewMemory((void **)zippath,28);
  strcpy(*zippath,dirseparator);
  strcat(*zippath,"usr");
  strcpy(*zippath,dirseparator);
  strcat(*zippath,"local");
  strcpy(*zippath,dirseparator);
  strcat(*zippath,"bin");
  strcat(*zippath,dirseparator);
  strcat(*zippath,"smokezip");
  if(STAT(*zippath,&statbuffer)==0)return;
  FREEMEMORY(*zippath);
#endif
  return;
}

/* ------------------ abortSV ------------------------ */

void abortSV(char *message){
  int i;
  if(message!=NULL&&strlen(message)>0){
    abort_message(message);
  }
  scanf("%i",&i);
}

/* ------------------ pauseSV ------------------------ */

void pauseSV(void){
  int i;
  printf("program paused - press <CTRL> c to close window\n");
  scanf("%i",&i);
}

/* ------------------ sv_startup ------------------------ */

void sv_startup(char *file, int showpart){

  int argc=2;
  char cprog[]="smokeview";
  char cnopart[]="-nopart";
  char *prog, *nopart;
  char *argv[3];
  prog=cprog;
  nopart=cnopart;
  CheckMemory;
  if(showpart==0)argc++;
  argv[0]=prog;
  argv[1]=file;
  if(showpart==0)argv[2]=nopart;
  CheckMemory;
  if(startup_pass==1)sv_startup_c(argc,argv);
  initcase_c(argc,argv);
}

/* ------------------ sv_init0 ------------------------ */

void sv_init0(void){
/*  int i;*/
  initcolors();
  setsmokeviewvars();
}
/* ------------------ sv_hide ------------------------ */

void sv_unload(void){
    int errorcode,i;
  /* first unload file data */
    for(i=0;i<nterraininfo;i++){
      readterrain("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<nsliceinfo;i++){
      readslice("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<nplot3dinfo;i++){
      readplot3d("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<npatchinfo;i++){
      readpatch(i,UNLOAD,&errorcode);
    }
    for(i=0;i<npartinfo;i++){
      readpart("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<nisoinfo;i++){
      readiso("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<nzone;i++){
      readzone(i,UNLOAD,&errorcode);
    }
    for(i=0;i<nsmoke3dinfo;i++){
      readsmoke3d(i,UNLOAD,&errorcode);
    }
    InitMenus(UNLOAD);
    readsmv(NULL,NULL);
    glutSetWindow(mainwindow_id);
    glutHideWindow();
    FreeAllMemory();
    return;
}


