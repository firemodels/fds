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
// RENDERONCE
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

/* ------------------ start_script ------------------------ */

void start_script(void){
  current_script_command=scriptinfo;
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


/* ------------------ get_file_type ------------------------ */

int get_file_type(char *file){
  int i;

  for(i=0;i<nfileinfo;i++){
    filedata *filei;

    filei = fileinfo + i;
    if(strcmp(filei->file,file)==0)return filei->type;
  }
  return -1;
}

/* ------------------ init_scripti ------------------------ */

void init_scripti(scriptdata *scripti, int command){
  scripti->command=command;
  scripti->cval=NULL;
  scripti->fval=0.0;
  scripti->ival=0;
}

/* ------------------ compile_script ------------------------ */

int compile_script(char *scriptfile){
  FILE *stream;
  char buffer[1024];
  scriptdata *scripti;
  int return_val;

  return_val=1;
  if(scriptfile==NULL){
    printf("*** internal smokeview error, scriptfile name is NULL\n");
    return return_val;
  }
  stream=fopen(scriptfile,"r");
  if(stream==NULL){
    printf("*** scriptfile, %s, could not be opened for input\n",scriptfile);
    return return_val;
  }

  return_val=0;
  
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

    if(match(buffer,"RENDERONCE",10) == 1){
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

  if(nscriptinfo==0){
    printf("*** warning: scriptfile has no usable commands\n");
    return 1;
  }

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

    if(match(buffer,"RENDERONCE",10) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,RENDERONCE);

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
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,LOADFILE);
      if(fgets(buffer,255,stream)==NULL)break;
      trim(buffer);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);
      filetype = get_file_type(buffer);
      scripti->ival = filetype;

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
  return return_val;
}


/* ------------------ run_script ------------------------ */

void run_script(void){
  scriptdata *scripti;

  if(current_script_command>scriptinfo+nscriptinfo-1){
    current_script_command=NULL;
    return;
  }
  scripti = current_script_command;
#ifdef _DEBUG
  printf("executing script command %i\n",scripti->command);
#endif
  switch (scripti->command){
    case RENDERONCE:
      keyboard('r',0,0);
      break;
    case RENDERALL:
      RenderMenu(1);
      break;
    case RENDERTYPE:
      break;
    case LOADFILE:
      break;
    case SETTIME:
      break;
    case SETFRAME:
      break;
    case VIEWPOINT:
      break;
  }
  current_script_command++;
}

#endif
