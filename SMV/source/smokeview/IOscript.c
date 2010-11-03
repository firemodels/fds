// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include "flowfiles.h"
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "MALLOC.h"
#include "smokeviewvars.h"

// svn revision character string
char IOscript_revision[]="$Revision$";
void remove_comment(char *buffer);
void ParticlePropShowMenu(int var);
void update_menu(void);
void glui_script_disable(void);

//
// script commands
//

// ---------- rendering images -----------

// default image file names are CHID_seq
//   where CHID is the base file name and seq is the frame number.
//   Smokeview will automatically add an .jpg or .png extension
//   depending on what kind of files are rendered.

// RENDERDIR
//  directory name (char) (where rendered files will go)

// RENDERONCE
// file name base (char) (or blank to use smokeview default)

// RENDERDOUBLEONCE
// file name base (char) (or blank to use smokeview default)

// RENDERALL 
//  skip (int)
// file name base (char) (or blank to use smokeview default)

// ---------- loading, unloading files -----------
//
//  Use LOADFILE to load a particular file.  Smokeview will figure
//  out what kind of file it is (3d smoke, slice etc.)  and call the
//  appropriate routine.
//
//  Use other LOAD commands to load files of the specified type for 
//  all meshes.

// LOADINIFILE 
//  file (char)

// LOADFILE 
//  file (char)

// LOADVFILE
//  file (char)

// LOADBOUNDARY
//   type (char)

// LOAD3DSMOKE
//  type (char)

// LOADPARTICLES

// PARTCLASSCOLOR
//   color (char)

// PARTCLASSTYPE
//   type (char)

// LOADISO
//  type (char)

// LOADSLICE
//  type (char)
//  1/2/3 (int)  val (float)

// LOADVSLICE
//  type (char)
//  1/2/3 (int)  val (float)

// LOADPLOT3D
//  mesh number (int) time (float)

// PLOT3DPROPS
//  plot3d type (int) showvector (0/1) (int) vector length index (int) plot3d display type (int)

// SHOWPLOT3DDATA
//  mesh number (int) orientation (int)  value (0/1) (int) position (float)

// UNLOADALL

// ---------- controlling scene -----------
//
// tours and viewpoints are referenced using names defined
// previously in a Smokeview session.  These names are
// stored in the .ini file.

// LOADTOUR
//  type (char)

// UNLOADTOUR

// SETTIMEVAL
//  time (float)

// SETVIEWPOINT
//  viewpoint (char)

// EXIT

/* ------------------ insert_scriptfile ------------------------ */

void get_newscriptfilename(char *newscriptfilename){
  char buffer[1024];
  char filebase[1024];
  int i;
  int nexti;
  scriptfiledata *scriptfile;

  for(i=0;i<1000;i++){
    if(i==0){
      strcpy(buffer,fdsprefix);
      strcat(buffer,".ssf");
    }
    else{
      sprintf(buffer,"%s_%03i.ssf",fdsprefix,i);
    }
    nexti=0;
    for(scriptfile=first_scriptfile.next;scriptfile->next!=NULL;scriptfile=scriptfile->next){
      if(strcmp(scriptfile->file,buffer)==0){
        nexti=1;
        break;
      }
    }
    if(nexti==0){
      strcpy(newscriptfilename,buffer);
      return;
    }
  }
  strcpy(newscriptfilename,"");
}

/* ------------------ get_scriptfilename ------------------------ */

char *get_scriptfilename(int id){
  scriptfiledata *thisptr,*prevptr,*nextptr;
  int len;
  scriptfiledata *scriptfile;

  for(scriptfile=first_scriptfile.next;scriptfile->next!=NULL;scriptfile=scriptfile->next){
    if(scriptfile->id==id)return scriptfile->file;
    if(scriptfile->file==NULL)continue;
  }
  return NULL;
}

/* ------------------ get_scriptfilename ------------------------ */

char *get_inifilename(int id){
  inifiledata *thisptr,*prevptr,*nextptr;
  int len;
  inifiledata *inifile;

  for(inifile=first_inifile.next;inifile->next!=NULL;inifile=inifile->next){
    if(inifile->id==id)return inifile->file;
    if(inifile->file==NULL)continue;
  }
  return NULL;
}

/* ------------------ insert_inifile ------------------------ */

inifiledata *insert_inifile(char *file){
  inifiledata *thisptr,*prevptr,*nextptr;
  int len;
  inifiledata *inifile;
  int idmax=-1;

  for(inifile=first_inifile.next;inifile->next!=NULL;inifile=inifile->next){
    if(inifile->id>idmax)idmax=inifile->id;
    if(inifile->file==NULL)continue;
    if(strcmp(file,inifile->file)==0)return NULL;
  }

  NewMemory((void **)&thisptr,sizeof(inifiledata));
  nextptr = &last_inifile;
  prevptr = nextptr->prev;
  nextptr->prev=thisptr;
  prevptr->next=thisptr;

  thisptr->next=nextptr;
  thisptr->prev=prevptr;
  thisptr->file=NULL;
  thisptr->id=idmax+1;

  if(file!=NULL){
    len = strlen(file);
    if(len>0){
      NewMemory((void **)&thisptr->file,len+1);
      strcpy(thisptr->file,file);
    }
  }
  return thisptr;
}

/* ------------------ insert_scriptfile ------------------------ */

scriptfiledata *insert_scriptfile(char *file){
  scriptfiledata *thisptr,*prevptr,*nextptr;
  int len;
  scriptfiledata *scriptfile;
  int idmax=-1;

  for(scriptfile=first_scriptfile.next;scriptfile->next!=NULL;scriptfile=scriptfile->next){
    if(scriptfile->id>idmax)idmax=scriptfile->id;
    if(scriptfile->file==NULL)continue;
    if(strcmp(file,scriptfile->file)==0)return NULL;
  }

  NewMemory((void **)&thisptr,sizeof(scriptfiledata));
  nextptr = &last_scriptfile;
  prevptr = nextptr->prev;
  nextptr->prev=thisptr;
  prevptr->next=thisptr;

  thisptr->next=nextptr;
  thisptr->prev=prevptr;
  thisptr->file=NULL;
  thisptr->recording=0;
  thisptr->id=idmax+1;

  if(file!=NULL){
    len = strlen(file);
    if(len>0){
      NewMemory((void **)&thisptr->file,len+1);
      strcpy(thisptr->file,file);
    }
  }
  return thisptr;
}

/* ------------------ cleanbuffer ------------------------ */

void cleanbuffer(char *buffer, char *buffer2){
  char *buff_ptr;

  remove_comment(buffer2);
  buff_ptr = trim_front(buffer2);
  trim(buff_ptr);
  strcpy(buffer,buff_ptr);
}

/* ------------------ start_script ------------------------ */

void start_script(void){
  if(scriptinfo==NULL){
    printf("*** warning: Smokeview script does not exist\n");
    return;
  }
  glui_script_disable();
  current_script_command=scriptinfo-1;
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
  scripti->ival2=0;
  scripti->ival3=0;
  scripti->ival4=0;
  scripti->ival5=0;
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
    if(fgets(buffer2,255,stream)==NULL)break;
    cleanbuffer(buffer,buffer2);

    if(strncmp(buffer," ",1)==0)continue;

    if(match_upper(buffer,"UNLOADALL",9) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"EXIT",4) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"RENDERDIR",9) == 1){
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
    if(match_upper(buffer,"LOADINIFILE",11) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADVFILE",9) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADBOUNDARY",12) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"PARTCLASSCOLOR",14) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"PARTCLASSTYPE",13) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"SHOWPLOT3DDATA",14) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"PLOT3DPROPS",11) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADTOUR",8) == 1){
      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"UNLOADTOUR",10) == 1){
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
    if(match_upper(buffer,"LOADPLOT3D",10) == 1){
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
    if(fgets(buffer2,255,stream)==NULL)break;
    cleanbuffer(buffer,buffer2);
    if(strncmp(buffer," ",1)==0)continue;

    if(match_upper(buffer,"UNLOADALL",9) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_UNLOADALL);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"RENDERDIR",9) == 1){
      int len;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_RENDERDIR);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      len = strlen(buffer);
      if(len>0){
        int i;

        NewMemory((void **)&scripti->cval,len+2);
        for(i=0;i<len;i++){
#ifdef WIN32
          if(buffer[i]=='/')buffer[i]='\\';
#else
          if(buffer[i]=='\\')buffer[i]='/';
#endif
        }
#ifdef WIN32
        if(buffer[len-1]!='\\')strcat(buffer,dirseparator);        
#else
        if(buffer[len-1]!='/')strcat(buffer,dirseparator);        
#endif
        strcpy(scripti->cval,buffer);
      }

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"RENDERONCE",10) == 1){
      int len;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_RENDERONCE);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      len = strlen(buffer);
      if(len>0){
        NewMemory((void **)&scripti->cval,len+1);
        strcpy(scripti->cval,buffer);
      }

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"RENDERDOUBLEONCE",16) == 1){
      int len;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_RENDERDOUBLEONCE);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      len = strlen(buffer);
      if(len>0){
        NewMemory((void **)&scripti->cval,len+1);
        strcpy(scripti->cval,buffer);
      }

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"RENDERALL",9) == 1){
      int len;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_RENDERALL);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      sscanf(buffer,"%i",&scripti->ival);
      if(scripti->ival<1)scripti->ival=1;
      if(scripti->ival>20)scripti->ival=20;

      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      len = strlen(buffer);
      if(len>0){
        NewMemory((void **)&scripti->cval,len+1);
        strcpy(scripti->cval,buffer);
      }

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADINIFILE",11) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_LOADINIFILE);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADFILE",8) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_LOADFILE);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADVFILE",9) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_LOADVFILE);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"EXIT",4) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_EXIT);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADBOUNDARY",12) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_LOADBOUNDARY);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"PARTCLASSCOLOR",14) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_PARTCLASSCOLOR);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"PARTCLASSTYPE",13) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_PARTCLASSTYPE);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"PLOT3DPROPS",11) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_PLOT3DPROPS);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      sscanf(buffer2,"%i %i %i %i",&scripti->ival,&scripti->ival2,&scripti->ival3,&scripti->ival4);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"SHOWPLOT3DDATA",14) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_SHOWPLOT3DDATA);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      sscanf(buffer2,"%i %i %i %i %f",&scripti->ival,&scripti->ival2,&scripti->ival3,&scripti->ival4,&scripti->fval);
      if(scripti->ival2==4){
        sscanf(buffer2,"%i %i %i %i %i",&scripti->ival,&scripti->ival2,&scripti->ival3,&scripti->ival4,&scripti->ival5);
      }

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOADTOUR",8) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_LOADTOUR);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"UNLOADTOUR",10) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_UNLOADTOUR);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"LOAD3DSMOKE",11) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_LOAD3DSMOKE);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
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
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
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
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      len=strlen(buffer);
      NewMemory((void **)&scripti->cval,len+1);
      strcpy(scripti->cval,buffer);

      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
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
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
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
    if(match_upper(buffer,"LOADPLOT3D",10) == 1){
      int len;
      int filetype;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_LOADPLOT3D);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      sscanf(buffer," %i %f",&scripti->ival,&scripti->fval);

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"SETTIMEVAL",10) == 1){
      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_SETTIMEVAL);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      sscanf(buffer,"%f",&scripti->fval);
      if(scripti->fval<0.0)scripti->fval=0.0;

      nscriptinfo++;
      continue;
    }
    if(match_upper(buffer,"SETVIEWPOINT",12) == 1){
      int len;

      scripti = scriptinfo + nscriptinfo;
      init_scripti(scripti,SCRIPT_SETVIEWPOINT);
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
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

  for(i=0;i<npart_files;i++){
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

  for(i=0;i<niso_files;i++){
    iso *isoi;
    int len;
    int imatch;

    isoi = isoinfo + i;
    len = strlen(scripti->cval);
    imatch = match_upper(isoi->surface_label.longlabel,scripti->cval,len);
    if(imatch==1||imatch==2){
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

  for(i=0;i<nsmoke3d_files;i++){
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

/* ------------------ script_loadtour ------------------------ */

void script_loadtour(scriptdata *scripti){
  int i;
  int errorcode;

  printf("Script: loading tour %s",scripti->cval);
  printf("\n");

  for(i=0;i<ntours;i++){
    tourdata *touri;

    touri = tourinfo + i;
    if(strcmp(touri->label,scripti->cval)==0){
      TourMenu(i);
      viewtourfrompath=0;
      TourMenu(-5);
      break;
    }
  }

  force_redisplay=1;
  updatemenu=1;
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

/* ------------------ script_partclasscolor ------------------------ */

void script_partclasscolor(scriptdata *scripti){
  int i;

  for(i=0;i<npart5prop;i++){
    part5prop *propi;

    propi = part5propinfo + i;
    if(propi->particle_property==0)continue;
    if(strcmp(propi->label->longlabel,scripti->cval)==0){
      ParticlePropShowMenu(i);
    }
  }
}


/* ------------------ script_plot3dprops ------------------------ */

void script_plot3dprops(scriptdata *scripti){
  int i, p_index;

  p_index = scripti->ival;
  if(p_index<1)p_index=1;
  if(p_index>5)p_index=5;

  visVector = scripti->ival2;
  if(visVector!=1)visVector=0;

  iveclengths = scripti->ival3;

  plotn = p_index;
  if(plotn<1){
    plotn=numplot3dvars;
  }
  if(plotn>numplot3dvars){
    plotn=1;
  }
  updateallplotslices();
  if(visiso==1)updatesurface();
  updateplot3dlistindex();

  vecfactor = get_vecfactor(&iveclengths);
  printf("iveclengths=%i\n",iveclengths);

  p3cont2d = scripti->ival4;
  if(p3cont2d<0)p3cont2d=0;
  if(p3cont2d>2)p3cont2d=2;
  update_plot3d_display();

  if(visVector==1&&ReadPlot3dFile==1){
    mesh *gbsave,*gbi;

    gbsave=current_mesh;
    for(i=0;i<nmeshes;i++){
      gbi = meshinfo + i;
      if(gbi->plot3dfilenum==-1)continue;
      update_current_mesh(gbi);
      updateplotslice(1);
      updateplotslice(2);
      updateplotslice(3);
    }
    update_current_mesh(gbsave);
  }
}

/* ------------------ script_showplot3ddata ------------------------ */

void script_showplot3ddata(scriptdata *scripti){
  int i;
  mesh *meshi;
  int imesh, dir, p_index, showhide;
  float val;
  int isolevel;

  imesh = scripti->ival-1;
  if(imesh<0||imesh>nmeshes-1)return;

  meshi = meshinfo + imesh;
  update_current_mesh(meshi);

  dir = scripti->ival2;
  if(dir<1)dir=1;
  if(dir>4)dir=4;

  plotn=scripti->ival3;

  showhide = scripti->ival4;
  val = scripti->fval;

  switch (dir){
    case 1:
      updateshowstep(showhide,DIRX);
      meshi->plotx=get_plot3d_index(meshi, dir, val);
      break;
    case 2:
      updateshowstep(showhide,DIRY);
      meshi->ploty=get_plot3d_index(meshi, dir, val);
      break;
    case 3:
      updateshowstep(showhide,DIRZ);
      meshi->plotz=get_plot3d_index(meshi, dir, val);
      break;
    case 4:
      isolevel=scripti->ival5;
      plotiso[plotn-1]=isolevel;
      updateshowstep(showhide,ISO);
      updatesurface();
      updatemenu=1;  
      break;
  }
  updateplotslice(dir);

}

/* ------------------ script_partclasstype ------------------------ */

void script_partclasstype(scriptdata *scripti){
  int i;

  for(i=0;i<npart5prop;i++){
    part5prop *propi;
    int j;

    propi = part5propinfo + i;
    if(propi->display==0)continue;
    for(j=0;j<npartclassinfo;j++){
      part5class *partclassj;

      if(propi->class_present[j]==0)continue;
      partclassj = partclassinfo + j;
      if(partclassj->kind==HUMANS)continue;
      if(strcmp(partclassj->name,scripti->cval)==0){
        ParticlePropShowMenu(-10-j);
      }
    }
  }
}

/* ------------------ script_loadinifile ------------------------ */

void script_loadinifile(scriptdata *scripti){
  int i;
  int errorcode;

  printf("Script: loading ini file %s",scripti->cval);
  printf("\n");
  scriptinifilename2=scripti->cval;
  windowresized=0;
  readini(2);
  scriptinifilename2=NULL;

}

/* ------------------ script_loadfile ------------------------ */

void script_loadfile(scriptdata *scripti){
  int i;
  int errorcode;

  printf("Script: loading file %s",scripti->cval);
  printf("\n");
  for(i=0;i<nslice_files;i++){
    slice *sd;

    sd = sliceinfo + i;
    if(strcmp(sd->file,scripti->cval)==0){
      readslice(sd->file,i,LOAD,&errorcode);
      return;
    }

  }
  for(i=0;i<npatch_files;i++){
    patch *patchi;

    patchi = patchinfo + i;
    if(strcmp(patchi->file,scripti->cval)==0){
      readpatch(i,LOAD,&errorcode);
      return;
    }
  }
  for(i=0;i<npart_files;i++){
    particle *parti;

    parti = partinfo + i;
    if(strcmp(parti->file,scripti->cval)==0){
      readpart(parti->file,i,LOAD,&errorcode);
      return;
    }
  }
  for(i=0;i<niso_files;i++){
    iso *isoi;

    isoi = isoinfo + i;
    if(strcmp(isoi->file,scripti->cval)==0){
      readiso(isoi->file,i,LOAD,&errorcode);
      return;
    }
  }
  for(i=0;i<nsmoke3d_files;i++){
    smoke3d *smoke3di;

    smoke3di = smoke3dinfo + i;
    if(strcmp(smoke3di->file,scripti->cval)==0){
      readsmoke3d(i,LOAD,&errorcode);
      return;
    }
  }
  for(i=0;i<nzone;i++){
    zone *zonei;
    char *file;

    zonei = zoneinfo + i;
    file = zonei->file;
    if(strcmp(file,scripti->cval)==0){
      readzone(file,i,LOAD,&errorcode);
      return;
    }
  }
  for(i=0;i<nplot3d_files;i++){
    plot3d *plot3di;

    plot3di = plot3dinfo + i;
    if(strcmp(plot3di->file,scripti->cval)==0){
      ReadPlot3dFile=1;
      readplot3d(plot3di->file,i,LOAD,&errorcode);
      update_menu();
      return;
    }
  }

  printf("file %s was not loaded\n",scripti->cval);

}


/* ------------------ script_loadplot3d ------------------------ */

void script_loadplot3d(scriptdata *scripti){
  int i;
  float time;
  int blocknum;

  time = scripti->fval;
  blocknum = scripti->ival-1;

  for(i=0;i<nplot3d_files;i++){
    plot3d *plot3di;

    plot3di = plot3dinfo + i;
    if(plot3di->blocknumber==blocknum&&fabs(plot3di->time-time)<0.5){
      LoadPlot3dMenu(i);
    }
  }
  updatecolors(-1);
  set_labels_controls();

  //update_menu();
}

/* ------------------ script_loadvfile ------------------------ */

void script_loadvfile(scriptdata *scripti){
  int i;
  int errorcode;

  printf("Script: loading vector slice file %s",scripti->cval);
  printf("\n");
  for(i=0;i<nvslice;i++){
    slice *val;
    vslice *vslicei;

    vslicei = vsliceinfo + i;
    val = sliceinfo + vslicei->ival;
    if(val==NULL)continue;
    if(strcmp(val->reg_file,scripti->cval)==0){
      LoadVSliceMenu(i);
      return;
    }
  }
  printf("vector slice file %s was not loaded\n",scripti->cval);

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
    if(timeval>times[ntimes-1]-0.0001)timeval=times[ntimes-1]-0.0001;
    for(i=0;i<ntimes-1;i++){
      if(times[i]<=timeval&&timeval<times[i+1]){
        itime=i;
        script_itime=i;
        stept=1;
        force_redisplay=1;
        update_framenumber(0);
        UpdateTimeLabels();
        break;
      }
    }
  }
}

/* ------------------ settimeval ------------------------ */

void settimeval(float timeval){
  int i;

  if(times!=NULL&&ntimes>0){
    if(timeval<times[0])timeval=times[0];
    if(timeval>times[ntimes-1]-0.0001)timeval=times[ntimes-1]-0.0001;
    for(i=0;i<ntimes;i++){
      float tlow, thigh;

      if(i==0){
        tlow = times[i];
        thigh = (times[i]+times[i+1])/2.0;
      }
      else if(i==ntimes-1){
        tlow = (times[i-1]+times[i])/2.0;
        thigh = times[i];
      }
      else{
        tlow=(times[i-1]+times[i])/2.0;
        thigh=(times[i]+times[i+1])/2.0;
      }
      if(tlow<=timeval&&timeval<thigh){
        itime=i;
        stept=1;
        force_redisplay=1;
        update_framenumber(0);
        UpdateTimeLabels();
        keyboard('t',0,0);
        break;
      }
    }
  }
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

int run_script(void){

// This procedure should return 1 if the smokeview frame should not be advanced.
// (to ensure images are rendered at the right time step)

  scriptdata *scripti;
  int returnval;

  returnval=0;

  if(current_script_command>scriptinfo+nscriptinfo-1){
    current_script_command=NULL;
    return returnval;
  }
  scripti = current_script_command;
#ifdef _DEBUG
  printf("executing script command %i\n",scripti->command);
#endif
  switch (scripti->command){
    case SCRIPT_UNLOADALL:
      LoadUnloadMenu(UNLOADALL);
      break;
    case SCRIPT_RENDERDIR:
      if(scripti->cval!=NULL&&strlen(scripti->cval)>0){
        script_dir_path=scripti->cval;
      }
      else{
        script_dir_path=NULL;
      }
      break;
    case SCRIPT_RENDERONCE:
      keyboard('r',0,0);
      returnval=1;
      break;
    case SCRIPT_RENDERDOUBLEONCE:
      keyboard('R',0,0);
      returnval=1;
      break;
    case SCRIPT_RENDERALL:
      script_renderall(scripti);
      break;
    case SCRIPT_LOADFILE:
      script_loadfile(scripti);
      break;
    case SCRIPT_LOADINIFILE:
      script_loadinifile(scripti);
      break;
    case SCRIPT_LOADVFILE:
      script_loadvfile(scripti);
      break;
    case SCRIPT_LOADBOUNDARY:
      script_loadboundary(scripti);
      break;
    case SCRIPT_PARTCLASSCOLOR:
      script_partclasscolor(scripti);
      break;
    case SCRIPT_SHOWPLOT3DDATA:
      script_showplot3ddata(scripti);
      break;
    case SCRIPT_PLOT3DPROPS:
      script_plot3dprops(scripti);
      break;
    case SCRIPT_PARTCLASSTYPE:
      script_partclasstype(scripti);
      break;
    case SCRIPT_LOADTOUR:
      script_loadtour(scripti);
      break;
    case SCRIPT_UNLOADTOUR:
      TourMenu(-2);
      break;
    case SCRIPT_EXIT:
#ifndef _DEBUG
      exit(0);
#endif
      break;
    case SCRIPT_LOADISO:
      script_loadiso(scripti);
      break;
    case SCRIPT_LOAD3DSMOKE:
      script_load3dsmoke(scripti);
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
    case SCRIPT_LOADPLOT3D:
      script_loadplot3d(scripti);
      break;
    case SCRIPT_SETTIMEVAL:
      returnval=1;
      script_settimeval(scripti);
      break;
    case SCRIPT_SETVIEWPOINT:
      script_setviewpoint(scripti);
      break;
  }
  GLUTPOSTREDISPLAY
  return returnval;
}
