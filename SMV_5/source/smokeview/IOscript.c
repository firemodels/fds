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
void remove_comment(char *buffer);
//
// script commands
//
// RENDERONCE

// RENDERDOUBLEONCE

// RENDERALL 
//  skip (int)

// SETJPEG

// SETPNG

// LOADFILE 
//  file (char)

// LOADBOUNDARY
//   type (char)

// LOAD3DSMOKE
//  type (char)

// LOADPARTICLES

// LOADISO
//  type (char)

// LOADSLICE
//  type (char)
//  1/2/3 (int)  val (float)

// LOADVSLICE
//  type (char)
//  1/2/3 (int)  val (float)

// SETIMEFRAME
//  frame (int)

// SETVIEWPOINT
//  viewpoint (char)

// UNLOADALL

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
  char buffer[1024], buffer2[1024];
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
    char *buff_ptr;

    if(fgets(buffer2,255,stream)==NULL)break;
    remove_comment(buffer2);
    buff_ptr = trim_front(buffer2);
    trim(buff_ptr);
    strcpy(buffer,buff_ptr);

    if(strncmp(buffer," ",1)==0)continue;

    if(match_upper(buffer,"UNLOADALL",9) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"RENDERONCE",10) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"RENDERDOUBLEONCE",16) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"RENDERALL",9) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADFILE",8) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADBOUNDARY",12) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOAD3DSMOKE",11) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADSLICE",9) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADVSLICE",10) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADISO",7) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADPARTICLES",13) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"SETTIMEVAL",10) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"SETVIEWPOINT",12) == 1){
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
    char *buff_ptr;

    if(fgets(buffer2,255,stream)==NULL)break;
    remove_comment(buffer2);
    buff_ptr = trim_front(buffer2);
    trim(buff_ptr);
    strcpy(buffer,buff_ptr);
    if(strncmp(buffer," ",1)==0)continue;

    if(match_upper(buffer,"UNLOADALL",9) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_UNLOADALL);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"RENDERONCE",10) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_RENDERONCE);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"RENDERDOUBLEONCE",16) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_RENDERDOUBLEONCE);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"RENDERALL",9) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_RENDERALL);
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%i",&scripti->ival);
      if(scripti->ival<1)scripti->ival=1;
      if(scripti->ival>10)scripti->ival=10;

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADFILE",8) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_LOADFILE);
      if(fgets(buffer,255,stream)==NULL)break;
      remove_comment(buffer);
      trim(buffer);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADBOUNDARY",12) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_LOADBOUNDARY);
      if(fgets(buffer,255,stream)==NULL)break;
      remove_comment(buffer);
      trim(buffer);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOAD3DSMOKE",11) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_LOAD3DSMOKE);
      if(fgets(buffer,255,stream)==NULL)break;
      remove_comment(buffer);
      trim(buffer);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADSLICE",9) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_LOADSLICE);
      if(fgets(buffer,255,stream)==NULL)break;
      remove_comment(buffer);
      trim(buffer);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      if(fgets(buffer,255,stream)==NULL)break;
      remove_comment(buffer);
      trim(buffer);
      sscanf(buffer,"%i %f",&scripti->ival,&scripti->fval);
      if(scripti->ival<1)scripti->ival=1;
      if(scripti->ival>3)scripti->ival=3;

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADVSLICE",10) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_LOADVSLICE);
      if(fgets(buffer,255,stream)==NULL)break;
      remove_comment(buffer);
      trim(buffer);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      if(fgets(buffer,255,stream)==NULL)break;
      remove_comment(buffer);
      trim(buffer);
      sscanf(buffer,"%i %f",&scripti->ival,&scripti->fval);
      if(scripti->ival<1)scripti->ival=1;
      if(scripti->ival>3)scripti->ival=3;

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADISO",7) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_LOADISO);
      if(fgets(buffer,255,stream)==NULL)break;
      remove_comment(buffer);
      trim(buffer);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADPARTICLES",13) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_LOADPARTICLES);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"SETTIMEVAL",10) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_SETTIMEVAL);
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%f",&scripti->fval);
      if(scripti->fval<0.0)scripti->fval=0.0;

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"SETVIEWPOINT",12) == 1){
      int len;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_SETVIEWPOINT);
      if(fgets(buffer,255,stream)==NULL)break;
      remove_comment(buffer);
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

/* ------------------ script_loadparticles ------------------------ */

void script_loadparticles(scriptdata *scripti){
  int i;
  int errorcode;

  printf("Script: loading particles files");
  printf("\n");

  for(i=0;i<npartinfo;i++){
    particle *parti;

    parti = partinfo + i;
    if(parti->evac==1)continue;
    if(parti->version==1){
      readpart(parti->file,i,LOAD,&errorcode);
    }
  }
  force_redisplay=1;
  update_framenumber(0);
  updatemenu=1;
}

/* ------------------ script_loadiso ------------------------ */

void script_loadiso(scriptdata *scripti){
  int i;
  int errorcode;

  printf("Script: loading isosurface files of type: %s",scripti->cval);
  printf("\n");

  for(i=0;i<niso;i++){
    iso *isoi;

    isoi = isoinfo + i;
    if(match_upper(isoi->label.longlabel,scripti->cval,strlen(scripti->cval))==1){
      readiso(isoi->file,i,LOAD,&errorcode);
    }
  }
  force_redisplay=1;
  updatemenu=1;

}

/* ------------------ script_load3dsmoke ------------------------ */

void script_load3dsmoke(scriptdata *scripti){
  int i;
  int errorcode;

  printf("Script: loading smoke3d files of type: %s",scripti->cval);
  printf("\n");

  for(i=0;i<nsmoke3d;i++){
    smoke3d *smoke3di;

    smoke3di = smoke3dinfo + i;
    if(match_upper(smoke3di->label.longlabel,scripti->cval,strlen(scripti->cval))==1){
      readsmoke3d(i,LOAD,&errorcode);
    }
  }
  force_redisplay=1;
  updatemenu=1;

}

/* ------------------ script_loadslice ------------------------ */

void script_loadslice(scriptdata *scripti){
  int i;
  int errorcode;

  printf("Script: loading slice files of type: %s",scripti->cval);
  printf("\n");

  for(i=0;i<nmultislices;i++){
    multislice *mslicei;
    slice *slicei;
    int j;
    float delta;

    mslicei = multisliceinfo + i;
    if(mslicei->nslices<=0)continue;
    slicei = sliceinfo + mslicei->islices[0];
    if(match_upper(slicei->label.longlabel,scripti->cval,strlen(scripti->cval))==0)continue;
    if(slicei->idir!=scripti->ival)continue;
    delta = slicei->position - scripti->fval;
    if(delta<0.0)delta = -delta;
    if(delta>slicei->delta)continue;

    for(j=0;j<mslicei->nslices;j++){
      LoadSliceMenu(mslicei->islices[j]);
    } 
    break;
  }


}


/* ------------------ script_loadvslice ------------------------ */

void script_loadvslice(scriptdata *scripti){
  int i;
  int errorcode;
  float delta;

  printf("Script: loading vector slice files of type: %s",scripti->cval);
  printf("\n");

  for(i=0;i<nmultivslices;i++){
    multivslice *mvslicei;
    int j;
    slice *slicei;

    mvslicei = multivsliceinfo + i;
    if(mvslicei->nvslices<=0)continue;
    slicei = sliceinfo + mvslicei->ivslices[0];
    if(match_upper(slicei->label.longlabel,scripti->cval,strlen(scripti->cval))==0)continue;
    if(slicei->idir!=scripti->ival)continue;
    delta = slicei->position - scripti->fval;
    if(delta<0.0)delta = -delta;
    if(delta>slicei->delta)continue;

    for(j=0;j<mvslicei->nvslices;j++){
      LoadVSliceMenu(mvslicei->ivslices[j]);
    } 
    break;
  }
}

/* ------------------ script_loadboundary ------------------------ */

void script_loadboundary(scriptdata *scripti){
  int i;
  int errorcode;

  printf("Script: loading boundary files of type: %s",scripti->cval);
  printf("\n");

  for(i=0;i<npatch_files;i++){
    patch *patchi;

    patchi = patchinfo + i;
    if(strcmp(patchi->label.longlabel,scripti->cval)==0){
      LOCK_COMPRESS
      readpatch(i,LOAD,&errorcode);
      UNLOCK_COMPRESS
    }
  }
  force_redisplay=1;
  updatemenu=1;
  update_framenumber(0);

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
    case SCRIPT_UNLOADALL:
      LoadUnloadMenu(UNLOADALL);
      break;
    case SCRIPT_RENDERONCE:
      keyboard('r',0,0);
      break;
    case SCRIPT_RENDERDOUBLEONCE:
      keyboard('R',0,0);
      break;
    case SCRIPT_RENDERALL:
      script_renderall(scripti);
      break;
    case SCRIPT_SETJPEG:
      RenderMenu(RenderJPEG);
      break;
    case SCRIPT_SETPNG:
      RenderMenu(RenderPNG);
      break;
    case SCRIPT_LOADFILE:
      script_loadfile(scripti);
      break;
    case SCRIPT_LOADBOUNDARY:
      script_loadboundary(scripti);
      break;
    case SCRIPT_LOADISO:
      script_loadiso(scripti);
      break;
    case SCRIPT_LOADPARTICLES:
      script_loadparticles(scripti);
      break;
    case SCRIPT_LOADSLICE:
      script_loadslice(scripti);
      break;
    case SCRIPT_LOADVSLICE:
      script_loadvslice(scripti);
      break;
    case SCRIPT_SETTIMEVAL:
      script_settimeval(scripti);
      break;
    case SCRIPT_SETVIEWPOINT:
      script_setviewpoint(scripti);
      break;
  }
  glutPostRedisplay();
  current_script_command++;
}

#endif
