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
// SETTIMEVAL
//  time (float)
// SETIMEFRAME
//  frame (int)
// SETVIEWPOINT
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
    if(match(buffer,"SETTIMEVAL",10) == 1){
      nscriptinfo++;
      continue;
    }
    if(match(buffer,"SETTIMEFRAME",12) == 1){
      nscriptinfo++;
      continue;
    }
    if(match(buffer,"SETVIEWPOINT",12) == 1){
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
    if(match(buffer,"SETTIMEVAL",10) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SETTIMEVAL);
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%f",&scripti->fval);
      if(scripti->fval<0.0)scripti->fval=0.0;


      nscriptinfo++;
      continue;
    }
    if(match(buffer,"SETTIMEFRAME",12) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SETTIMEFRAME);
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%i",&scripti->ival);
      if(scripti->ival<0)scripti->ival=0.0;

      nscriptinfo++;
      continue;
    }
    if(match(buffer,"SETVIEWPOINT",12) == 1){
      int len;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SETVIEWPOINT);
      if(fgets(buffer,255,stream)==NULL)break;
      trim(buffer);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      nscriptinfo++;
      continue;
    }
  }
  fclose(stream);
  return return_val;
}

/* ------------------ run_renderall ------------------------ */

void script_renderall(scriptdata *scripti){
  int skip;

  skip=scripti->ival;
  if(skip<1)skip=1;
  printf("Script: Rendering every %i frames",skip);
  printf("\n");
  RenderMenu(skip);
}

/* ------------------ script_rendertype ------------------------ */

void script_rendertype(scriptdata *scripti){
  printf("Script: rendertype=%s",scripti->cval);
  printf(" *** not implemented ***");
  printf("\n");
}

/* ------------------ script_loadfile ------------------------ */

void script_loadfile(scriptdata *scripti){
  int i;
  int errorcode;

  printf("Script: loading file %s",scripti->cval);
  printf("\n");
  for(i=0;i<nslice;i++){
    slice *sd;

    sd = sliceinfo + i;
    if(strcmp(sd->reg_file,scripti->cval)==0){
      readslice(sd->reg_file,i,LOAD,&errorcode);
      return;
    }

  }
}

/* ------------------ script_settimeval ------------------------ */

void script_settimeval(scriptdata *scripti){
  float timeval;
  int i;

  timeval = scripti->fval;
  printf("Script: set time to %f",timeval);
  printf("\n");
  if(times!=NULL&&ntimes>0){
    if(timeval<times[0])timeval=times[0];
    if(timeval>times[ntimes-1])timeval=times[ntimes-1];
    for(i=0;i<ntimes-1;i++){
      if(times[i]<=timeval&&timeval<=times[i+1]){
        itime=i;
        stept=1;
        force_redisplay=1;

      }
    }
  }
}

/* ------------------ script_settimeframe ------------------------ */

void script_settimeframe(scriptdata *scripti){
  int timeframe;

  timeframe = scripti->ival;
  printf("Script: set time frame to %i",timeframe);
  printf(" *** not implemented ***");
  printf("\n");
}

/* ------------------ script_setviewpoint ------------------------ */

void script_setviewpoint(scriptdata *scripti){
  char *viewpoint;
  camera *ca;

  viewpoint = scripti->cval;
  printf("Script: set viewpoint to %s",viewpoint);
  printf("\n");
  for(ca=camera_list_first.next;ca->next!=NULL;ca=ca->next){
    if(strcmp(scripti->cval,ca->name)==0){
      ResetMenu(ca->view_id);
      break;
    }
  }
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
      script_renderall(scripti);
      break;
    case RENDERTYPE:
      script_rendertype(scripti);
      break;
    case LOADFILE:
      script_loadfile(scripti);
      break;
    case SETTIMEVAL:
      script_settimeval(scripti);
      break;
    case SETTIMEFRAME:
      script_settimeframe(scripti);
      break;
    case SETVIEWPOINT:
      script_setviewpoint(scripti);
      break;
  }
  current_script_command++;
}

#endif
