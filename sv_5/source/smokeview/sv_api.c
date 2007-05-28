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
#include "isodefs.h"
#include "smokeheaders.h"
#include "smokeviewapi.h"
#include "flowfiles.h"
#include "smokeviewdefs.h"
#include "MALLOC.h"
#include "smokeviewvars.h"

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
        OpenSMVFile(filename,filelength,&openfile);
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

/* ------------------ sv_grid ------------------------ */

void svWINAPI sv_grid(int val){
  visGrid = 0;
  if(val==1)visGrid=1;
}

/* ------------------ get_smokezippath ------------------------ */

void get_smokezippath(char *progdir, char **zippath){
  size_t len;
  struct stat statbuffer;

  if(progdir!=NULL){
    len=strlen(progdir);
    NewMemory((void **)zippath,len+13);
    strcpy(*zippath,progdir);
#ifdef WIN32
    strcat(*zippath,"smokezip.exe");
#else
    strcat(*zippath,"smokezip");
#endif
  }
  if(stat(*zippath,&statbuffer)==0)return;
  FREEMEMORY(*zippath);
#ifdef WIN32
  NewMemory((void **)zippath,28);
  strcpy(*zippath,"c:\\nist\\fds\\smokezip.exe");
  if(stat(*zippath,&statbuffer)==0)return;
  FREEMEMORY(*zippath);
#else
  NewMemory((void **)zippath,28);
  strcpy(*zippath,"~");
  strcat(*zippath,dirseparator);
  strcat(*zippath,"bin");
  strcat(*zippath,dirseparator);
  strcat(*zippath,"smokezip");
  if(stat(*zippath,&statbuffer)==0)return;
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
  if(stat(*zippath,&statbuffer)==0)return;
  FREEMEMORY(*zippath);
#endif
  return;
}

/* ------------------ pause ------------------------ */

void pauseSV(void){
  int i;
  printf("program paused - hit any key and \"enter\" to continue\n");
  scanf("%i",&i);
}

/* ------------------ sv_startup ------------------------ */

void svWINAPI sv_startup(char *file, int showpart){

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
  if(first==1)sv_startup_c(argc,argv);
  initcase_c(argc,argv);
}

/* ------------------ ResizeWindow ------------------------ */

void ResizeWindow(int width, int height){
#ifdef pp_RENDER
  if(render_double!=0)return;
#endif
  glutSetWindow(mainwindow_id);
  glutReshapeWindow(width,height);
  glutPostRedisplay();
}
    
/* ------------------ sv_update ------------------------ */

void svWINAPI sv_update(){
  glutMainLoop();
}

/* ------------------ sv_init0 ------------------------ */

void svWINAPI sv_init0(void){
/*  int i;*/
  initcolors();
  setsmokeviewvars();
}
/* ------------------ sv_hide ------------------------ */

void svWINAPI sv_unload(void){
    int errorcode,i;
  /* first unload file data */
#ifdef pp_WUI
    for(i=0;i<nterraininfo;i++){
      readterrain("",i,UNLOAD,&errorcode);
    }
#endif
    for(i=0;i<nslice;i++){
      readslice("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<nplot3d;i++){
      readplot("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<npatch_files;i++){
      readpatch(i,UNLOAD,&errorcode);
    }
    for(i=0;i<npartinfo;i++){
      readpart("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<niso;i++){
      readiso("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<nzone;i++){
      readzone("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<nsmoke3d;i++){
      readsmoke3d(i,UNLOAD,&errorcode);
    }
    glutDetachMenu(GLUT_RIGHT_BUTTON);
    InitMenus(UNLOAD);
    readsmv(NULL);
    glutSetWindow(mainwindow_id);
    glutHideWindow();

    return;
}


