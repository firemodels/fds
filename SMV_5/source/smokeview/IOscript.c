// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>  
#include <stdlib.h>
#include <sys/stat.h>
#include "flowfiles.h"
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "MALLOC.h"
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

// svn revision character string
char IOscript_revision[]="$Revision$";
#ifdef pp_SCRIPT
//
// script commands
//
// RENDERONE
// RENDERALL 
//  skip (int)
// RENDERTYPE
//  type (char - jpeg or png)
// LOADFILE 
//  file (char)
// SETTIME
//  time (float)
// SETFRAME
//  frame (int)
// VIEWPOINT
//  viewpoint (char)

#define RENDERONE 101
#define RENDERALL 102
#define RENDERTYPE 103
#define LOADFILE 104
#define SETTIME 105
#define SETFRAME 106
#define VIEWPOINT 107

/* ------------------ start_script ------------------------ */

void start_script(void){
  printf("starting script\n");
}

/* ------------------ free_script ------------------------ */

void free_script(void){
  scriptdata *scripti;
  int i;

  if(nscriptinfo>0){
    for(i=0;i<nscriptinfo;i++){
      scripti = scriptinfo + i;

      FREEMEMORY(scripti->cval);  
    }
    FREEMEMORY(scriptinfo);
    nscriptinfo=0;
  }

}

/* ------------------ init_scripti ------------------------ */

void init_scripti(scriptdata *scripti, int command){
  scripti->command=command;
  scripti->cval=NULL;
  scripti->fval=0.0;
  scripti->ival=0;
}

/* ------------------ compile_script ------------------------ */

void compile_script(char *scriptfile){
  FILE *stream;
  char buffer[1024];
  scriptdata *scripti;

  if(scriptfile==NULL)return;
  stream=fopen(scriptfile,"r");
  if(stream==NULL)return;
  
  /* 
   ************************************************************************
   ************************ start of pass 1 ********************************* 
   ************************************************************************
 */

  free_script();

  while(!feof(stream)){
    scriptdata *scripti;

    if(fgets(buffer,255,stream)==NULL)break;
    if(strncmp(buffer," ",1)==0)continue;

    if(match(buffer,"RENDERONE",9) == 1){
      nscriptinfo++;
      continue;
    }
    if(match(buffer,"RENDERALL",9) == 1){
      nscriptinfo++;
      continue;
    }
    if(match(buffer,"LOADFILE",8) == 1){
      nscriptinfo++;
      continue;
    }
    if(match(buffer,"SETTIME",7) == 1){
      nscriptinfo++;
      continue;
    }
    if(match(buffer,"SETFRAME",8) == 1){
      nscriptinfo++;
      continue;
    }
    if(match(buffer,"VIEWPOINT",9) == 1){
      nscriptinfo++;
      continue;
    }
  }

  if(nscriptinfo==0)return;

  NewMemory((void **)&scriptinfo,nscriptinfo*sizeof(scriptdata));

  /* 
   ************************************************************************
   ************************ start of pass 2 ********************************* 
   ************************************************************************
 */

  nscriptinfo=0;
  rewind(stream);
  while(!feof(stream)){
    if(fgets(buffer,255,stream)==NULL)break;
    if(strncmp(buffer," ",1)==0)continue;

    if(match(buffer,"RENDERONE",9) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,RENDERONE);

      nscriptinfo++;
      continue;
    }
    if(match(buffer,"RENDERALL",9) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,RENDERALL);
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%i",&scripti->ival);
      if(scripti->ival<1)scripti->ival=1;
      if(scripti->ival>10)scripti->ival=10;

      nscriptinfo++;
      continue;
    }
    if(match(buffer,"LOADFILE",8) == 1){
      int len;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,LOADFILE);
      if(fgets(buffer,255,stream)==NULL)break;
      trim(buffer);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      nscriptinfo++;
      continue;
    }
    if(match(buffer,"SETTIME",7) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SETTIME);
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%f",&scripti->fval);
      if(scripti->fval<0.0)scripti->fval=0.0;


      nscriptinfo++;
      continue;
    }
    if(match(buffer,"SETFRAME",8) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SETFRAME);
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%i",&scripti->ival);
      if(scripti->ival<0)scripti->ival=0.0;

      nscriptinfo++;
      continue;
    }
    if(match(buffer,"VIEWPOINT",9) == 1){
      int len;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,VIEWPOINT);
      if(fgets(buffer,255,stream)==NULL)break;
      trim(buffer);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      nscriptinfo++;
      continue;
    }
  }
}

/* ------------------ run_script ------------------------ */

void run_script(void){
}

#endif
