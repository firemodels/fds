#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include GLUT_H

#include "update.h"
#include "smokeviewvars.h"
#include "IOvolsmoke.h"

#define RENDER_START 3

void update_menu(void);

/* ------------------ get_newscriptfilename ------------------------ */

void get_newscriptfilename(char *newscriptfilename){
  char buffer[1024];
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
  scriptfiledata *scriptfile;

  for(scriptfile=first_scriptfile.next;scriptfile->next!=NULL;scriptfile=scriptfile->next){
    if(scriptfile->id==id)return scriptfile->file;
    if(scriptfile->file==NULL)continue;
  }
  return NULL;
}

/* ------------------ get_inifilename ------------------------ */

char *get_inifilename(int id){
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
    if(strcmp(file,scriptfile->file)==0)return scriptfile;
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

/* ------------------ insert_luascriptfile ------------------------ */

#ifdef pp_LUA
luascriptfiledata *insert_luascriptfile(char *file){
  luascriptfiledata *thisptr,*prevptr,*nextptr;
  int len;
  luascriptfiledata *luascriptfile;
  int idmax=-1;

  for(luascriptfile=first_luascriptfile.next;luascriptfile->next!=NULL;luascriptfile=luascriptfile->next){
    if(luascriptfile->id>idmax)idmax=luascriptfile->id;
    if(luascriptfile->file==NULL)continue;
    if(strcmp(file,luascriptfile->file)==0)return luascriptfile;
  }

  NewMemory((void **)&thisptr,sizeof(luascriptfiledata));
  nextptr = &last_luascriptfile;
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
#endif

/* ------------------ start_script ------------------------ */

void start_script(void){
  if(scriptinfo==NULL){
    fprintf(stderr,"*** Error: Smokeview script does not exist\n");
    return;
  }
  glui_script_disable();
  current_script_command=scriptinfo-1;
}

/* ------------------ get_pointer ------------------------ */

char *get_pointer(char *buffer2){
  char *cval=NULL, *buffptr;
  int len;

  buffptr = remove_comment(buffer2);
  len = strlen(buffptr);
  if(len>0){
    NewMemory((void **)&cval,len+1);
    strcpy(cval,buffptr);
  }
  return cval;
}

/* ------------------ free_script ------------------------ */

void free_script(void){
  scriptdata *scripti;
  int i;

  if(nscriptinfo>0){
    for(i=0;i<nscriptinfo;i++){
      scripti = scriptinfo + i;

      FREEMEMORY(scripti->cval);
      FREEMEMORY(scripti->cval2);
    }
    FREEMEMORY(scriptinfo);
    nscriptinfo=0;
  }

}

/* ------------------ init_scripti ------------------------ */

void init_scripti(scriptdata *scripti, int command,char *label){
  char *label2;

  trim_back(label);
  label2 = trim_front(label);
  strcpy(scripti->command_label,label2);
  scripti->command=command;
  scripti->cval=NULL;
  scripti->cval2=NULL;
  scripti->fval=0.0;
  scripti->ival=0;
  scripti->ival2=0;
  scripti->ival3=0;
  scripti->ival4=0;
  scripti->ival5=0;
}

/* ------------------ get_script_keyword_index ------------------------ */

int get_script_keyword_index(char *keyword){
  if(keyword==NULL||strlen(keyword)==0)return SCRIPT_UNKNOWN;

  if(match_upper(keyword,"CBARFLIP") == MATCH)return SCRIPT_CBARFLIP;
  if(match_upper(keyword,"CBARNORMAL") == MATCH)return SCRIPT_CBARNORMAL;
  if(match_upper(keyword,"EXIT") == MATCH)return SCRIPT_EXIT;
  if(match_upper(keyword,"KEYBOARD") == MATCH)return SCRIPT_KEYBOARD;
  if(match_upper(keyword,"GSLICEORIEN") == MATCH)return SCRIPT_GSLICEORIEN;
  if(match_upper(keyword,"GSLICEPOS") == MATCH)return SCRIPT_GSLICEPOS;
  if(match_upper(keyword,"GSLICEVIEW") == MATCH)return SCRIPT_GSLICEVIEW;
  if(match_upper(keyword,"LOAD3DSMOKE") == MATCH)return SCRIPT_LOAD3DSMOKE;
  if(match_upper(keyword,"LOADBOUNDARY") == MATCH)return SCRIPT_LOADBOUNDARY;
  if(match_upper(keyword,"LOADBOUNDARYM") == MATCH)return SCRIPT_LOADBOUNDARYM;
  if(match_upper(keyword,"LOADFILE") == MATCH)return SCRIPT_LOADFILE;
  if(match_upper(keyword,"LABEL") == MATCH)return SCRIPT_LABEL;
  if(match_upper(keyword,"LOADINIFILE") == MATCH)return SCRIPT_LOADINIFILE;
  if(match_upper(keyword,"LOADISO") == MATCH)return SCRIPT_LOADISO;
  if(match_upper(keyword,"LOADISOM") == MATCH)return SCRIPT_LOADISOM;
  if(match_upper(keyword,"LOADPARTICLES") == MATCH)return SCRIPT_LOADPARTICLES;
  if(match_upper(keyword,"LOADPLOT3D") == MATCH)return SCRIPT_LOADPLOT3D;
  if(match_upper(keyword,"LOADSLICE") == MATCH)return SCRIPT_LOADSLICE;
  if(match_upper(keyword,"LOADSLICEM") == MATCH)return SCRIPT_LOADSLICEM;
  if(match_upper(keyword,"LOADTOUR") == MATCH)return SCRIPT_LOADTOUR;
  if(match_upper(keyword,"LOADVOLSMOKE") == MATCH)return SCRIPT_LOADVOLSMOKE;
  if(match_upper(keyword,"LOADVOLSMOKEFRAME") == MATCH)return SCRIPT_LOADVOLSMOKEFRAME;
  if(match_upper(keyword,"LOADVFILE") == MATCH)return SCRIPT_LOADVFILE;
  if(match_upper(keyword,"LOADVSLICE") == MATCH)return SCRIPT_LOADVSLICE;
  if(match_upper(keyword,"LOADVSLICEM") == MATCH)return SCRIPT_LOADVSLICEM;
  if(match_upper(keyword,"MAKEMOVIE") == MATCH)return SCRIPT_MAKEMOVIE;
  if(match_upper(keyword,"PARTCLASSCOLOR") == MATCH)return SCRIPT_PARTCLASSCOLOR;
  if(match_upper(keyword,"PARTCLASSTYPE") == MATCH)return SCRIPT_PARTCLASSTYPE;
  if(match_upper(keyword,"PLOT3DPROPS") == MATCH)return SCRIPT_PLOT3DPROPS;
  if(match_upper(keyword,"RENDERALL") == MATCH)return SCRIPT_RENDERALL;
  if (match_upper(keyword, "RENDER360ALL") == MATCH)return SCRIPT_RENDER360ALL;
  if(match_upper(keyword,"RENDERCLIP") == MATCH)return SCRIPT_RENDERCLIP;
  if(match_upper(keyword,"RENDERDIR") == MATCH)return SCRIPT_RENDERDIR;
  if(match_upper(keyword,"RENDERTYPE") == MATCH)return SCRIPT_RENDERTYPE;
  if(match_upper(keyword,"MOVIETYPE") == MATCH)return SCRIPT_MOVIETYPE;
  if(match_upper(keyword, "RENDERSIZE") == MATCH)return SCRIPT_RENDERSIZE;
  if(match_upper(keyword,"RENDERDOUBLEONCE") == MATCH)return SCRIPT_RENDERDOUBLEONCE;
  if(match_upper(keyword,"RENDERONCE") == MATCH)return SCRIPT_RENDERONCE;
  if(match_upper(keyword,"RENDERSTART") == MATCH)return SCRIPT_RENDERSTART;
  if(match_upper(keyword,"SCENECLIP") == MATCH)return SCRIPT_SCENECLIP;
  if(match_upper(keyword,"SETTOURKEYFRAME") == MATCH)return SCRIPT_SETTOURKEYFRAME;
  if(match_upper(keyword,"SETTOURVIEW") == MATCH)return SCRIPT_SETTOURVIEW;
  if(match_upper(keyword,"SETTIMEVAL") == MATCH)return SCRIPT_SETTIMEVAL;
  if(match_upper(keyword,"SETVIEWPOINT") == MATCH)return SCRIPT_SETVIEWPOINT;
  if(match_upper(keyword,"SHOWPLOT3DDATA") == MATCH)return SCRIPT_SHOWPLOT3DDATA;
  if(match_upper(keyword,"UNLOADALL") == MATCH)return SCRIPT_UNLOADALL;
  if(match_upper(keyword,"UNLOADTOUR") == MATCH)return SCRIPT_UNLOADTOUR;
  if(match_upper(keyword,"VOLSMOKERENDERALL") == MATCH)return SCRIPT_VOLSMOKERENDERALL;
  if(match_upper(keyword, "ISORENDERALL")==1)return SCRIPT_ISORENDERALL;
  if(match_upper(keyword, "XSCENECLIP")==1)return SCRIPT_XSCENECLIP;
  if(match_upper(keyword,"YSCENECLIP") == MATCH)return SCRIPT_YSCENECLIP;
  if(match_upper(keyword,"ZSCENECLIP") == MATCH)return SCRIPT_ZSCENECLIP;

  return SCRIPT_UNKNOWN;
}

/* ------------------ get_xyz ------------------------ */

void get_xyz(char *buffer,int *ival){
  int i;

  for(i=0;i<strlen(buffer);i++){
    char *c;

    c = buffer+i;
    if(*c!=' '){
      if(*c=='x'||*c=='X'||*c=='1'){
        *ival=1;
        return;
      }
      if(*c=='y'||*c=='Y'||*c=='2'){
        *ival=2;
        return;
      }
      if(*c=='z'||*c=='Z'||*c=='3'){
        *ival=3;
        return;
      }
      if(*c == 'a' || *c == 'A' || *c == '0'){
        *ival = 0;
        return;
      }
      *ival = 1;
      break;
    }
  }
}

/* ------------------ script_error_check ------------------------ */

void script_error_check(char *keyword, char *data){
  if(get_script_keyword_index(data)!=SCRIPT_UNKNOWN){
    fprintf(stderr,"*** Error: While parsing the Smokeview script entry: %s ,\n",keyword);
    fprintf(stderr,"           a keyword was found in \"%s\", data was expected.\n",data);
  }
}

#define SETbuffer \
if(fgets(buffer, 255, stream) == NULL){\
scriptEOF = 1; \
break; \
}\
buffptr = remove_comment(buffer); \
buffptr = trim_front(buffptr); \
script_error_check(keyword, buffptr)

#define SETcval \
SETbuffer;\
scripti->cval=get_pointer(buffptr)

#define SETcval2 \
SETbuffer;\
scripti->cval2 = get_pointer(buffptr)

#define SETfval \
SETbuffer;\
sscanf(buffptr, "%f", &scripti->fval)

#define SETival \
SETbuffer;\
sscanf(buffptr, "%i", &scripti->ival)

#define SETival2 \
SETbuffer;\
sscanf(buffptr, "%i", &scripti->ival2)

#ifndef pp_DEG
/* ------------------ removeDEG ------------------------ */

void removeDEG(char *string){
  int i,ii;

  if(string == NULL)return;
  for(i = 0,ii=0; i < strlen(string);i++){
    unsigned char c;

    c = (unsigned char)string[i];
    if(c == 176)continue;
    string[ii] = string[i];
    ii++;
  }
  string[ii] = 0;
}
#endif

/* ------------------ compile_script ------------------------ */

int compile_script(char *scriptfile){
  FILE *stream;
  int return_val;

  return_val=1;
  if(scriptfile==NULL){
    fprintf(stderr,"*** Error: scriptfile name is NULL\n");
    return return_val;
  }
  stream=fopen(scriptfile,"r");
  if(stream==NULL){
    fprintf(stderr,"*** Error: scriptfile, %s, could not be opened for input\n",scriptfile);
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
    char buffer[1024], buffer2[1024], *buffptr;

    if(fgets(buffer2,255,stream)==NULL)break;
    buffptr = remove_comment(buffer2);
    strcpy(buffer, buffptr);


    if(get_script_keyword_index(buffer)!=SCRIPT_UNKNOWN)nscriptinfo++;
  }

  if(nscriptinfo==0){
    fclose(stream);
    fprintf(stderr,"*** Error: scriptfile has no usable commands\n");
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
    int keyword_index;
    int scriptEOF;
    char keyword[255];
    char buffer[1024], buffer2[1024], *buffptr;
    scriptdata *scripti;

    if(fgets(buffer2,255,stream)==NULL)break;
    buffptr = remove_comment(buffer2);
    strcpy(buffer, buffptr);

    if(strlen(buffer)==0)continue;

    keyword_index = get_script_keyword_index(buffer);
    if(keyword_index==SCRIPT_UNKNOWN)continue;
    strcpy(keyword,buffer);

    scripti = scriptinfo + nscriptinfo;
    init_scripti(scripti,keyword_index,buffer);

    scriptEOF=0;
    switch(keyword_index){

// UNLOADALL
      case SCRIPT_UNLOADALL:

// LOADPARTICLES
      case SCRIPT_LOADPARTICLES:

// CBARFLIP:
      case SCRIPT_CBARFLIP:

// CBARNORMAL:
      case SCRIPT_CBARNORMAL:
        break;

// RENDERSIZE
// width height (int)
      case SCRIPT_RENDERSIZE:
        SETbuffer;
        sscanf(buffer, "%i %i", &scripti->ival, &scripti->ival2);
        break;

// RENDERTYPE
//  jpg or png  (char)
      case SCRIPT_RENDERTYPE:
        SETcval;
        if(STRCMP(scripti->cval, "PNG") == 0){
          scripti->ival = PNG;
        }
        else{
          scripti->ival = JPEG;
        }
        break;

// MOVIETYPE
//  jpg, png or wmv  (char)
      case SCRIPT_MOVIETYPE:
        SETcval;
        if(STRCMP(scripti->cval, "WMV") == 0){
          scripti->ival = WMV;
        }
        if(STRCMP(scripti->cval, "MP4") == 0){
          scripti->ival = MP4;
        }
        else{
          scripti->ival = AVI;
        }
        break;

// RENDERDIR
//  directory name (char) (where rendered files will go)
      case SCRIPT_RENDERDIR:
        {
        int len;
        int i;

        SETbuffer;
        len = strlen(buffer);
        if(len>0){
#ifdef WIN32
          for(i=0;i<len;i++){
            if(buffer[i]=='/')buffer[i]='\\';
          }
          if(buffer[len-1]!='\\')strcat(buffer,dirseparator);
#else
          for(i=0;i<len;i++){
            if(buffer[i]=='\\')buffer[i]='/';
          }
          if(buffer[len-1]!='/')strcat(buffer,dirseparator);
#endif
          scripti->cval=get_pointer(buffer);
        }
        }
        break;

// SCENECLIP
//  clip mode (int)
      case SCRIPT_SCENECLIP:

// LOADVOLSMOKE
//  mesh number (-1 for all meshes) (int)
      case SCRIPT_LOADVOLSMOKE:
        SETival;
        break;

// X/y/ZSCENECLIP
// imin (int) min (float) imax (int) max (float)
      case SCRIPT_XSCENECLIP:
      case SCRIPT_YSCENECLIP:
      case SCRIPT_ZSCENECLIP:
        SETbuffer;
        sscanf(buffer,"%i %f %i %f",&scripti->ival,&scripti->fval,&scripti->ival2,&scripti->fval2);
        break;

// RENDERCLIP
// flag left right bottom top indentations in pixels, clip if flag==1
      case SCRIPT_RENDERCLIP:
        SETbuffer;
        sscanf(buffer,"%i %i %i %i %i",&scripti->ival,&scripti->ival2,&scripti->ival3,&scripti->ival4, &scripti->ival5);
        break;

// RENDERONCE
// file name base (char) (or blank to use smokeview default)
      case SCRIPT_RENDERONCE:

// RENDERDOUBLEONCE
// file name base (char) (or blank to use smokeview default)
      case SCRIPT_RENDERDOUBLEONCE:
        SETcval2;
        break;

// RENDERSTART
//  start_frame (int) skip_frame (int)
      case SCRIPT_RENDERSTART:
        SETbuffer;
        sscanf(buffer,"%i %i",&scripti->ival,&scripti->ival2);
        break;

// RENDERALL
//  skip (int)
// file name base (char) (or blank to use smokeview default)
      case SCRIPT_RENDERALL:
        SETbuffer;
        scripti->ival=1;   // skip
        scripti->ival3=0;  // first frame
        sscanf(buffer,"%i %i",&scripti->ival,&scripti->ival3);
        scripti->ival = MAX(scripti->ival, 1);
        scripti->ival3 = MAX(scripti->ival3, 0);
        first_frame_index=scripti->ival3;

        SETcval2;
        break;

// RENDER360ALL
//  skip (int)
// file name base (char) (or blank to use smokeview default)
      case SCRIPT_RENDER360ALL:
        SETbuffer;
        scripti->ival = 1;   // skip
        scripti->ival3 = 0;  // first frame
        sscanf(buffer, "%i %i", &scripti->ival, &scripti->ival3);
        scripti->ival = MAX(scripti->ival, 1);
        scripti->ival3 = MAX(scripti->ival3, 0);
        first_frame_index = scripti->ival3;

        SETcval2;
        break;
        
        // VOLSMOKERENDERALL
//  skip (int) start_frame (int)
// file name base (char) (or blank to use smokeview default)
      case SCRIPT_VOLSMOKERENDERALL:
        SETbuffer;
        scripti->ival3=0;  // first frame
        scripti->ival=1;
        sscanf(buffer,"%i %i",&scripti->ival,&scripti->ival3);
        scripti->ival=CLAMP(scripti->ival,1,20); // skip
        scripti->exit=0;
        scripti->first=1;
        scripti->remove_frame=-1;
        first_frame_index=scripti->ival3;

        SETcval2;
        break;

// LOADISOM
//  type (char)
//  mesh number (int)
      case SCRIPT_LOADISOM:
        SETcval;
#ifndef pp_DEG
        removeDEG(scripti->cval);
#endif
        scripti->ival = 1;
        SETival;
        break;

// ISORENDERALL
//  skip (int) start_frame (int) iso file index (int) ( index of &ISOF line in .fds input file)
// file name base (char) (or blank to use smokeview default)
      case SCRIPT_ISORENDERALL:
        SETbuffer;
        scripti->ival3 = 0;  // first frame
        scripti->ival = 1;
        sscanf(buffer, "%i %i %i", &scripti->ival, &scripti->ival3, &scripti->ival4);
        scripti->ival=CLAMP(scripti->ival,1,20); // skip
        scripti->exit = 0;
        scripti->first = 1;
        scripti->remove_frame = -1;
        first_frame_index = scripti->ival3;

        SETcval2;  //render file base name
        break;

// MAKEMOVIE
//  movie_name (char)
//  frame_prefix (char)
//  framerate (float)
      case SCRIPT_MAKEMOVIE:
        SETcval;
        SETcval2;
        SETfval;
        break;

// EXIT
      case SCRIPT_EXIT:

// UNLOADTOUR
      case SCRIPT_UNLOADTOUR:
        scripti->cval=NULL;
        break;

      case SCRIPT_KEYBOARD:

// LOADINIFILE
//  file (char)
      case SCRIPT_LOADINIFILE:

 // LOADFILE
//  file (char)
      case SCRIPT_LOADFILE:

// LOADVFILE
//  file (char)
      case SCRIPT_LOADVFILE:

// PARTCLASSCOLOR
//   color (char)
      case SCRIPT_PARTCLASSCOLOR:

// PARTCLASSTYPE
//   type (char)
      case SCRIPT_PARTCLASSTYPE:

// LOADTOUR
//  type (char)
      case SCRIPT_LOADTOUR:

// LOAD3DSMOKE
//  type (char)
      case SCRIPT_LOAD3DSMOKE:

// LOADISO
//  type (char)
      case SCRIPT_LOADISO:

// SETVIEWPOINT
//  viewpoint (char)
      case SCRIPT_SETVIEWPOINT:

// LABEL
//   text
      case SCRIPT_LABEL:

// LOADBOUNDARY
//   type (char)
      case SCRIPT_LOADBOUNDARY:
        SETcval;
        break;

// LOADBOUNDARYM
//  type (char)
//  mesh number (int)
      case SCRIPT_LOADBOUNDARYM:
        SETcval;
        SETival;
        break;

// PLOT3DPROPS
//  plot3d type (int) showvector (0/1) (int) vector length index (int) plot3d display type (int)
      case SCRIPT_PLOT3DPROPS:
        SETbuffer;
        {
          float fv=-1;

          sscanf(buffer,"%i %i %i %i %f",&scripti->ival,&scripti->ival2,&scripti->ival3,&scripti->ival4,&fv);
          if(scripti->ival3<0&&fv>=0.0){
            scripti->fval=fv;
          }
          else{
            scripti->fval=-1.0;
          }
        }
        break;

// SHOWPLOT3DDATA
//  mesh number (int) orientation (int)  value (0/1) (int) position (float)
      case SCRIPT_SHOWPLOT3DDATA:
        SETbuffer;
        sscanf(buffer,"%i %i %i %i %f",&scripti->ival,&scripti->ival2,&scripti->ival3,&scripti->ival4,&scripti->fval);
        if(scripti->ival2==4){
          sscanf(buffer,"%i %i %i %i %i",&scripti->ival,&scripti->ival2,&scripti->ival3,&scripti->ival4,&scripti->ival5);
        }
        break;

// LOADVOLSMOKEFRAME
//  mesh index, frame (int)
      case SCRIPT_LOADVOLSMOKEFRAME:
        SETbuffer;
        sscanf(buffer,"%i %i",&scripti->ival,&scripti->ival2);
        break;

// LOADSLICE
//  type (char)
//  1/2/3 (int)  val (float)
      case SCRIPT_LOADSLICE:

// LOADVSLICE
//  type (char)
//  1/2/3 (int)  val (float)
      case SCRIPT_LOADVSLICE:
        SETcval;

        SETbuffer;
        sscanf(buffer, "%i %f", &scripti->ival, &scripti->fval);
        scripti->ival = CLAMP(scripti->ival, 0, 3);
        break;

// LOADVSLICEM
//  type (char)
//  1/2/3 (int)  val (float)
//  mesh number (int)
      case SCRIPT_LOADVSLICEM:
        SETcval;

        SETbuffer;
        sscanf(buffer, "%i %f", &scripti->ival, &scripti->fval);
        scripti->ival = CLAMP(scripti->ival, 0, 3);
        SETival2;
        break;

// LOADSLICEM
//  type (char)
//  1/2/3 (int)  val (float)
//  mesh number (int)
      case SCRIPT_LOADSLICEM:
        SETcval;

        SETbuffer;
        sscanf(buffer, "%i %f", &scripti->ival, &scripti->fval);
        scripti->ival = CLAMP(scripti->ival, 0, 3);
        SETival2;
        break;

// LOADPLOT3D
//  mesh number (int) time (float)
      case SCRIPT_LOADPLOT3D:
        SETbuffer;
        sscanf(buffer," %i %f",&scripti->ival,&scripti->fval);
        break;

// SETTIMEVAL
//  time (float)
      case SCRIPT_SETTIMEVAL:
        SETfval;
        if(scripti->fval<0.0)scripti->fval=0.0;
        break;

// SETTOURVIEW
//   viewtype  showpath showtour_locus tension
      case SCRIPT_SETTOURVIEW:
        SETbuffer;
        sscanf(buffer,"%i %i %i %f",&scripti->ival,&scripti->ival2,&scripti->ival3,&scripti->fval);
        break;

// SETTOURKEYFRAME
//  time (float)
      case SCRIPT_SETTOURKEYFRAME:
        SETfval;
        break;

// GSLICEVIEW
// show_gslice (int) show_triangles (int)  show_triangulation (int) show_normals (int)
      case SCRIPT_GSLICEVIEW:
        SETbuffer;
        sscanf(buffer,"%i %i %i %i",&scripti->ival,&scripti->ival2,&scripti->ival3,&scripti->ival4);
        break;

// GSLICEPOS
// x (float) y (float) z (float)
      case SCRIPT_GSLICEPOS:
        SETbuffer;
        sscanf(buffer,"%f %f %f",&scripti->fval,&scripti->fval2,&scripti->fval3);
        break;

// GSLICEORIEN
// azimuth (float) elevation (float)
      case SCRIPT_GSLICEORIEN:
        SETbuffer;
        sscanf(buffer,"%f %f",&scripti->fval,&scripti->fval2);
        break;
    }
    if(scriptEOF==1)break;
    if(keyword_index!=SCRIPT_UNKNOWN)nscriptinfo++;
  }
  fclose(stream);
  return return_val;
}

/* ------------------ script_renderstart ------------------------ */

void script_renderstart(scriptdata *scripti){
  script_startframe=scripti->ival;
  script_skipframe=scripti->ival2;
}

/* ------------------ script_renderall ------------------------ */

void script_renderall(scriptdata *scripti){
  int skip_local;


  if(script_startframe>0)scripti->ival3=script_startframe;
  if(startframe0>=0)scripti->ival3=startframe0;
  first_frame_index=scripti->ival3;
  itimes=first_frame_index;

  if(script_skipframe>0)scripti->ival=script_skipframe;
  if(skipframe0>0)scripti->ival=skipframe0;
  skip_local=MAX(1,scripti->ival);

  PRINTF("script: Rendering every %i frame(s) starting at frame %i\n\n",skip_local,scripti->ival3);
  skip_render_frames=1;
  RenderMenu(skip_local);
}

/* ------------------ script_render360all ------------------------ */

void script_render360all(scriptdata *scripti) {
  int skip_local;


  if (script_startframe>0)scripti->ival3 = script_startframe;
  if (startframe0 >= 0)scripti->ival3 = startframe0;
  first_frame_index = scripti->ival3;
  itimes = first_frame_index;

  if (script_skipframe>0)scripti->ival = script_skipframe;
  if (skipframe0>0)scripti->ival = skipframe0;
  skip_local = MAX(1, scripti->ival);

  PRINTF("script: Rendering every %i frame(s) starting at frame %i\n\n", skip_local, scripti->ival3);
  skip_render_frames = 1;
  //RenderMenu(skip_local);
  render_360=1;
  Render_CB(RENDER_START);
}

/* ------------------ script_loadvolsmokeframe ------------------------ */

void script_loadvolsmokeframe(scriptdata *scripti, int flag){
  int framenum, index;
  int first = 1;
  int i;

  index = scripti->ival;
  framenum = scripti->ival2;
  if(index > nmeshes - 1)index = -1;
  for(i = 0; i < nmeshes; i++){
    if(index == i || index < 0){
      meshdata *meshi;
      volrenderdata *vr;

      meshi = meshinfo + i;
      vr = &meshi->volrenderinfo;
      free_volsmoke_frame(vr, framenum);
      read_volsmoke_frame(vr, framenum, &first);
      if(vr->times_defined == 0){
        vr->times_defined = 1;
        get_volsmoke_all_times(vr);
      }
      vr->loaded = 1;
      vr->display = 1;
    }
  }
  plotstate = GetPlotState(DYNAMIC_PLOTS);
  stept = 1;
  UpdateTimes();
  force_redisplay = 1;
  UpdateFrameNumber(framenum);
  i = framenum;
  itimes = i;
  script_itime = i;
  stept = 1;
  force_redisplay = 1;
  UpdateFrameNumber(0);
  UpdateTimeLabels();
  keyboard('r', FROM_SMOKEVIEW);
  if(flag == 1)script_render = 1;// called when only rendering a single frame
}

/* ------------------ script_loadvolsmokeframe2 ------------------------ */

void script_loadvolsmokeframe2(void){
  scriptdata scripti;

  scripti.ival = -1;
  scripti.ival2 = itimes;
  script_loadvolsmokeframe(&scripti, 0);
}

/* ------------------ script_volsmokerenderall ------------------------ */

void script_volsmokerenderall(scriptdata *scripti){
  int skip_local;

  if(nvolrenderinfo==0){
    PRINTF("*** Error: there is no volume rendered smoke data to render\n");
    ScriptMenu(SCRIPT_CANCEL);
    return;
  }
  script_loadvolsmokeframe2();

  if(script_startframe>0)scripti->ival3=script_startframe;
  if(vol_startframe0>0)scripti->ival3=vol_startframe0;
  // check first_frame_index
  first_frame_index=scripti->ival3;
  itimes=first_frame_index;

  if(script_skipframe>0)scripti->ival=script_skipframe;
  if(vol_skipframe0>0)scripti->ival=vol_skipframe0;
  skip_local=MAX(1,scripti->ival);

  PRINTF("script: Rendering every %i frame(s) starting at frame %i\n\n",skip_local,scripti->ival3);
  skip_render_frames=1;
  scripti->ival=skip_local;
  RenderMenu(skip_local);
}

/* ------------------ script_loadisoframe ------------------------ */

void script_loadisoframe(scriptdata *scripti, int flag){
  int framenum, index;
  int i, count = 0;
  int fileindex;

  index = scripti->ival;
  framenum = scripti->ival2;
  fileindex = scripti->ival4;
  if(index > nmeshes - 1)index = -1;

  update_readiso_geom_wrapup = UPDATE_ISO_START_ALL;
  for(i = 0; i < nisoinfo; i++){
    int errorcode;
    isodata *isoi;

    isoi = isoinfo + i;
    if(isoi->isof_index + 1 == fileindex){
      readiso(isoi->file, i, LOAD, &framenum, &errorcode);
      count++;
    }
  }
  if(update_readiso_geom_wrapup == UPDATE_ISO_ALL_NOW)readiso_geom_wrapup();
  update_readiso_geom_wrapup = UPDATE_ISO_OFF;

  plotstate = GetPlotState(DYNAMIC_PLOTS);
  stept = 1;
  UpdateTimes();
  force_redisplay = 1;
  UpdateFrameNumber(framenum);
  i = framenum;
  itimes = i;
  script_itime = i;
  stept = 1;
  force_redisplay = 1;
  UpdateFrameNumber(0);
  UpdateTimeLabels();
  keyboard('r', FROM_SMOKEVIEW);
  if(flag == 1)script_render = 1;// called when only rendering a single frame
}

/* ------------------ script_loadisoframe2 ------------------------ */

void script_loadisoframe2(scriptdata *scripti){
  scripti->ival2 = itimes;
  script_loadisoframe(scripti, 0);
}

/* ------------------ script_isorenderall ------------------------ */

void script_isorenderall(scriptdata *scripti){
  int skip_local;

  //  skip (int) start_frame (int) iso file index (int) ( index of &ISOF line in .fds input file)
  //  skip == scripti->ival
  //  start_frame == scripti->ival3
  //  &ISOF index == scripti->ival4

  if(nisoinfo == 0){
    PRINTF("*** Error: there is no isosurface data to render\n");
    ScriptMenu(SCRIPT_CANCEL);
    return;
  }
  script_loadisoframe2(scripti);

  if(script_startframe>0)scripti->ival3 = script_startframe;
  if(startframe0>0)scripti->ival3 = startframe0;
  // check first_frame_index
  first_frame_index = scripti->ival3;
  itimes = first_frame_index;

  if(script_skipframe>0)scripti->ival = script_skipframe;
  if(skipframe0>0)scripti->ival = skipframe0;
  skip_local = MAX(1, scripti->ival);

  PRINTF("script: Rendering every %i frame(s) starting at frame %i\n\n", skip_local, scripti->ival3);
  skip_render_frames = 1;
  scripti->ival = skip_local;
  RenderMenu(skip_local);
}

/* ------------------ script_makemovie ------------------------ */

void script_makemovie(scriptdata *scripti){
  strcpy(movie_name, scripti->cval);
  strcpy(render_file_base,scripti->cval2);
  movie_framerate=scripti->fval;
  Render_CB(MAKE_MOVIE);
}

/* ------------------ script_loadparticles ------------------------ */

void script_loadparticles(scriptdata *scripti){
  int i;
  int errorcode;
  int count=0;

  FREEMEMORY(loaded_file);

  PRINTF("script: loading particles files\n\n");

  npartframes_max=get_min_partframes();
  for(i=0;i<npartinfo;i++){
    partdata *parti;

    parti = partinfo + i;
    if(parti->evac==1)continue;
    readpart(parti->file,i,UNLOAD,PARTDATA,&errorcode);
    count++;
  }
  for(i=0;i<npartinfo;i++){
    partdata *parti;

    parti = partinfo + i;
    if(parti->evac==1)continue;
    readpart(parti->file,i,LOAD,PARTDATA,&errorcode);
    if(scripti->cval!=NULL&&strlen(scripti->cval)>0){
      FREEMEMORY(loaded_file);
      NewMemory((void **)&loaded_file,strlen(scripti->cval)+1);
      strcpy(loaded_file,scripti->cval);
    }
    count++;
  }
  if(count==0)fprintf(stderr,"*** Error: Particles files failed to load\n");
  force_redisplay=1;
  UpdateFrameNumber(0);
  updatemenu=1;
}

/* ------------------ script_loadiso ------------------------ */

void script_loadiso(scriptdata *scripti, int meshnum){
  int i;
  int count=0;

  FREEMEMORY(loaded_file);
  PRINTF("script: loading isosurface files of type: %s\n\n",scripti->cval);

  update_readiso_geom_wrapup = UPDATE_ISO_START_ALL;
  for(i = 0; i<nisoinfo; i++){
    int errorcode;
    isodata *isoi;
    char label2[100];
    int lencval, lenlabel;

    isoi = isoinfo + i;
    if(meshnum != -1 && isoi->blocknumber+1 != meshnum)continue;
    lencval = strlen(scripti->cval);
    lenlabel = strlen(isoi->surface_label.longlabel);
    if(lencval<=lenlabel){
      strncpy(label2, isoi->surface_label.longlabel, lencval);
      label2[lencval] = 0;
      if(STRCMP(label2, scripti->cval)==0){
        readiso(isoi->file, i, LOAD, NULL, &errorcode);
        if(scripti->cval!=NULL&&strlen(scripti->cval)>0){
          FREEMEMORY(loaded_file);
          NewMemory((void **)&loaded_file, strlen(scripti->cval)+1);
          strcpy(loaded_file, scripti->cval);
        }
        count++;
      }
    }
  }
  if(update_readiso_geom_wrapup == UPDATE_ISO_ALL_NOW)readiso_geom_wrapup();
  update_readiso_geom_wrapup = UPDATE_ISO_OFF;
  if(count == 0)fprintf(stderr, "*** Error: Isosurface files of type %s failed to load\n", scripti->cval);
  force_redisplay=1;
  updatemenu=1;
}

/* ------------------ script_loadvolsmoke ------------------------ */

void script_loadvolsmoke(scriptdata *scripti){
  int imesh;

  imesh = scripti->ival;
  if(imesh==-1){
    read_vol_mesh=VOL_READALL;
    read_volsmoke_allframes_allmeshes2(NULL);
  }
  else if(imesh>=0&&imesh<nmeshes){
    meshdata *meshi;
    volrenderdata *vr;

    meshi = meshinfo + imesh;
    vr = &meshi->volrenderinfo;
    read_volsmoke_allframes(vr);
  }
}

/* ------------------ script_load3dsmoke ------------------------ */

void script_load3dsmoke(scriptdata *scripti){
  int i;
  int errorcode;
  int count=0;

  FREEMEMORY(loaded_file);
  PRINTF("script: loading smoke3d files of type: %s\n\n",scripti->cval);

  for(i=0;i<nsmoke3dinfo;i++){
    smoke3ddata *smoke3di;

    smoke3di = smoke3dinfo + i;
    if(match_upper(smoke3di->label.longlabel,scripti->cval) == MATCH){
      readsmoke3d(i,LOAD,&errorcode);
      if(scripti->cval!=NULL&&strlen(scripti->cval)>0){
        FREEMEMORY(loaded_file);
        NewMemory((void **)&loaded_file,strlen(scripti->cval)+1);
        strcpy(loaded_file,scripti->cval);
      }
      count++;
    }
  }
  if(count==0)fprintf(stderr,"*** Error: Smoke3d files of type %s failed to load\n",scripti->cval);
  force_redisplay=1;
  updatemenu=1;

}

/* ------------------ script_loadslice ------------------------ */

void script_loadslice(scriptdata *scripti){
  int i;
  int count=0;

  PRINTF("script: loading slice files of type: %s\n\n",scripti->cval);

  for(i=0;i<nmultisliceinfo;i++){
    multislicedata *mslicei;
    slicedata *slicei;
    int j;

    mslicei = multisliceinfo + i;
    if(mslicei->nslices<=0)continue;
    slicei = sliceinfo + mslicei->islices[0];
    if(match_upper(slicei->label.longlabel,scripti->cval) == NOTMATCH)continue;
    if(scripti->ival==0){
      if(slicei->volslice==0)continue;
    }
    else{
      if(slicei->idir != scripti->ival)continue;
      if(ABS(slicei->position_orig - scripti->fval) > slicei->delta_orig)continue;
    }
    for(j=0;j<mslicei->nslices;j++){
      LoadSliceMenu(mslicei->islices[j]);
      count++;
    }
    break;
  }
  if(count==0)fprintf(stderr,"*** Error: Slice files of type %s failed to load\n",scripti->cval);
}


/* ------------------ script_loadslicem ------------------------ */

void script_loadslicem(scriptdata *scripti, int meshnum){
  int i;

  PRINTF("script: loading slice files of type: %s in mesh %i\n\n", scripti->cval,meshnum);

  for(i = 0; i < nsliceinfo; i++){
    slicedata *slicei;

    slicei = sliceinfo + i;
    if(slicei->blocknumber + 1 != meshnum)continue;
    if(match_upper(slicei->label.longlabel, scripti->cval) == NOTMATCH)continue;
    if(scripti->ival == 0){
      int *min, *max;
      meshdata *meshi;

      if(slicei->volslice == 0)continue;
      min = slicei->ijk_min;
      max = slicei->ijk_max;
      if(min[0] != 0 || min[1] != 0 || min[2] != 0)continue;
      meshi = meshinfo + slicei->blocknumber;
      if(max[0] != meshi->ibar || max[1] != meshi->jbar || max[2] != meshi->kbar)continue;
    }
    else{
      if(slicei->idir != scripti->ival)continue;
      if(ABS(slicei->position_orig - scripti->fval) > slicei->delta_orig)continue;
    }
    LoadSliceMenu(i);
  }
}

/* ------------------ script_loadvslice ------------------------ */

void script_loadvslice(scriptdata *scripti){
  int i;
  int count=0;

  PRINTF("script: loading vector slice files of type: %s\n\n",scripti->cval);

  for(i=0;i<nmultivsliceinfo;i++){
    multivslicedata *mvslicei;
    vslicedata *vslicei;
    int j;
    slicedata *slicei;

    mvslicei = multivsliceinfo + i;
    if(mvslicei->nvslices<=0)continue;
    vslicei = vsliceinfo + mvslicei->ivslices[0];
    slicei = sliceinfo + vslicei->ival;
    if(match_upper(slicei->label.longlabel,scripti->cval) == NOTMATCH)continue;
    if(scripti->ival == 0){
      if(slicei->volslice == 0)continue;
    }
    else{
      if(slicei->idir != scripti->ival)continue;
      if(ABS(slicei->position_orig - scripti->fval) > slicei->delta_orig)continue;
    }
    for(j=0;j<mvslicei->nvslices;j++){
      LoadVSliceMenu(mvslicei->ivslices[j]);
      count++;
    }
    break;
  }
  if(count==0)fprintf(stderr,"*** Error: Vector slice files of type %s failed to load\n",scripti->cval);
}

/* ------------------ script_loadvslicem ------------------------ */

void script_loadvslicem(scriptdata *scripti, int meshnum){
  int i;
  int count=0;

  PRINTF("script: loading vector slice files of type: %s in mesh %i\n\n", scripti->cval,meshnum);

  for(i=0;i<nmultivsliceinfo;i++){
    multivslicedata *mvslicei;
    vslicedata *vslicei;
    int j;
    slicedata *slicei;

    mvslicei = multivsliceinfo + i;
    if(mvslicei->nvslices<=0)continue;
    vslicei = vsliceinfo + mvslicei->ivslices[0];
    slicei = sliceinfo + vslicei->ival;
    if(slicei->blocknumber + 1 != meshnum)continue;
    if(match_upper(slicei->label.longlabel,scripti->cval) == NOTMATCH)continue;
    if(scripti->ival == 0){
      if(slicei->volslice == 0)continue;
    }
    else{
      if(slicei->idir != scripti->ival)continue;
      if(ABS(slicei->position_orig - scripti->fval) > slicei->delta_orig)continue;
    }
    for(j=0;j<mvslicei->nvslices;j++){
      LoadVSliceMenu(mvslicei->ivslices[j]);
      count++;
    }
    break;
  }
  if(count==0)fprintf(stderr,"*** Error: Vector slice files of type %s in mesh %i failed to load\n",scripti->cval,meshnum);
}

/* ------------------ script_loadtour ------------------------ */

void script_loadtour(scriptdata *scripti){
  int i;
  int count=0;

  PRINTF("script: loading tour %s\n\n",scripti->cval);

  for(i=0;i<ntours;i++){
    tourdata *touri;

    touri = tourinfo + i;
    if(strcmp(touri->label,scripti->cval)==0){
      TourMenu(i);
      viewtourfrompath=0;
      TourMenu(MENU_TOUR_VIEWFROMROUTE);
      count++;
      break;
    }
  }

  if(count==0)fprintf(stderr,"*** Error: The tour %s failed to load\n",scripti->cval);
  force_redisplay=1;
  updatemenu=1;
}

/* ------------------ script_loadboundary ------------------------ */

void script_loadboundary(scriptdata *scripti, int meshnum){
  int i;
  int errorcode;
  int count=0;

  FREEMEMORY(loaded_file);
  PRINTF("Script: loading boundary files of type: %s\n\n",scripti->cval);

  for(i=0;i<npatchinfo;i++){
    patchdata *patchi;

    patchi = patchinfo + i;
    if (meshnum == -1 || patchi->blocknumber + 1 == meshnum) {
      if (strcmp(patchi->label.longlabel, scripti->cval) == 0) {
        LOCK_COMPRESS
          readpatch(i, LOAD, &errorcode);
        if (scripti->cval != NULL&&strlen(scripti->cval) > 0) {
          FREEMEMORY(loaded_file);
          NewMemory((void **)&loaded_file, strlen(scripti->cval) + 1);
          strcpy(loaded_file, scripti->cval);
        }
        count++;
        UNLOCK_COMPRESS
          if(meshnum == -1)break;
      }
    }
  }
  if(count==0)fprintf(stderr,"*** Error: Boundary files of type %s failed to load\n",scripti->cval);
  force_redisplay=1;
  updatemenu=1;
  UpdateFrameNumber(0);

}

/* ------------------ script_partclasscolor ------------------------ */

void script_partclasscolor(scriptdata *scripti){
  int i;
  int count=0;

  for(i=0;i<npart5prop;i++){
    partpropdata *propi;

    propi = part5propinfo + i;
    if(propi->particle_property==0)continue;
    if(strcmp(propi->label->longlabel,scripti->cval)==0){
      ParticlePropShowMenu(i);
      count++;
    }
  }
  if(count==0)fprintf(stderr,"*** Error: particle class color: %s failed to be set\n",scripti->cval);
}


/* ------------------ script_plot3dprops ------------------------ */

void script_plot3dprops(scriptdata *scripti){
  int i, p_index;

  p_index = scripti->ival;
  if(p_index<1)p_index=1;
  if(p_index>5)p_index=5;

  visVector = scripti->ival2;
  if(visVector!=1)visVector=0;

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

  vecfactor=1.0;
  if(scripti->fval>=0.0)vecfactor=scripti->fval;
  update_vector_widgets();

  PRINTF("script: vecfactor=%f\n",vecfactor);

  contour_type=CLAMP(scripti->ival4,0,2);
  update_plot3d_display();

  if(visVector==1&&ReadPlot3dFile==1){
    meshdata *gbsave,*gbi;

    gbsave=current_mesh;
    for(i=0;i<nmeshes;i++){
      gbi = meshinfo + i;
      if(gbi->plot3dfilenum==-1)continue;
      update_current_mesh(gbi);
      updateplotslice(XDIR);
      updateplotslice(YDIR);
      updateplotslice(ZDIR);
    }
    update_current_mesh(gbsave);
  }
}

/* ------------------ script_showplot3ddata ------------------------ */

void script_showplot3ddata(scriptdata *scripti){
  meshdata *meshi;
  int imesh, dir, showhide;
  float val;
  int isolevel;

  imesh = scripti->ival-1;
  if(imesh<0||imesh>nmeshes-1)return;

  meshi = meshinfo + imesh;
  update_current_mesh(meshi);

  dir = CLAMP(scripti->ival2,XDIR,ISO);

  plotn=scripti->ival3;

  showhide = scripti->ival4;
  val = scripti->fval;

  switch(dir){
    case XDIR:
      visx_all=showhide;
      iplotx_all=get_index(val,XDIR,plotx_all,nplotx_all);
      next_xindex(1,0);
      next_xindex(-1,0);
      break;
    case YDIR:
      visy_all=showhide;
      iploty_all=get_index(val,YDIR,ploty_all,nploty_all);
      next_yindex(1,0);
      next_yindex(-1,0);
      break;
    case ZDIR:
      visz_all=showhide;
      iplotz_all=get_index(val,ZDIR,plotz_all,nplotz_all);
      next_zindex(1,0);
      next_zindex(-1,0);
      break;
    case ISO:
      isolevel=scripti->ival5;
      plotiso[plotn-1]=isolevel;
      updateshowstep(showhide,ISO);
      updatesurface();
      updatemenu=1;
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
  updateplotslice(dir);

}

/* ------------------ script_partclasstype ------------------------ */

void script_partclasstype(scriptdata *scripti){
  int i;
  int count=0;

  for(i=0;i<npart5prop;i++){
    partpropdata *propi;
    int j;

    propi = part5propinfo + i;
    if(propi->display==0)continue;
    for(j=0;j<npartclassinfo;j++){
      partclassdata *partclassj;

      if(propi->class_present[j]==0)continue;
      partclassj = partclassinfo + j;
      if(partclassj->kind==HUMANS)continue;
      if(strcmp(partclassj->name,scripti->cval)==0){
        ParticlePropShowMenu(-10-j);
        count++;
      }
    }
  }
  if(count==0)fprintf(stderr,"*** Error: particle class type %s failed to be set\n",scripti->cval);
}

/* ------------------ script_loadinifile ------------------------ */

void script_loadinifile(scriptdata *scripti){
  PRINTF("script: loading ini file %s\n\n",scripti->cval);
  windowresized=0;
  readini(scripti->cval);
}

/* ------------------ script_loadfile ------------------------ */

void script_loadfile(scriptdata *scripti){
  int i;
  int errorcode;

  FREEMEMORY(loaded_file);
  PRINTF("script: loading file %s\n\n",scripti->cval);
  if(scripti->cval!=NULL&&strlen(scripti->cval)>0){
    NewMemory((void **)&loaded_file,strlen(scripti->cval)+1);
    strcpy(loaded_file,scripti->cval);
  }
  for(i=0;i<nsliceinfo;i++){
    slicedata *sd;

    sd = sliceinfo + i;
    if(strcmp(sd->file,scripti->cval)==0){
      if(i<nsliceinfo-nfedinfo){
        readslice(sd->file,i,LOAD,SET_SLICECOLOR,&errorcode);
      }
      else{
        readfed(i,LOAD,FED_SLICE,&errorcode);
      }
      return;
    }
  }
  for(i=0;i<npatchinfo;i++){
    patchdata *patchi;

    patchi = patchinfo + i;
    if(strcmp(patchi->file,scripti->cval)==0){
      readpatch(i,LOAD,&errorcode);
      return;
    }
  }
  npartframes_max=get_min_partframes();
  for(i=0;i<npartinfo;i++){
    partdata *parti;

    parti = partinfo + i;
    if(strcmp(parti->file,scripti->cval)==0){
      readpart(parti->file,i,LOAD,PARTDATA,&errorcode);
      return;
    }
  }
  for(i=0;i<nisoinfo;i++){
    isodata *isoi;

    isoi = isoinfo + i;
    if(strcmp(isoi->file,scripti->cval)==0){
      readiso(isoi->file,i,LOAD,NULL,&errorcode);
      if(update_readiso_geom_wrapup == UPDATE_ISO_ONE_NOW)readiso_geom_wrapup();
      return;
    }
  }
  for(i=0;i<nsmoke3dinfo;i++){
    smoke3ddata *smoke3di;

    smoke3di = smoke3dinfo + i;
    if(strcmp(smoke3di->file,scripti->cval)==0){
      readsmoke3d(i,LOAD,&errorcode);
      return;
    }
  }
  for(i=0;i<nzoneinfo;i++){
    zonedata *zonei;

    zonei = zoneinfo + i;
    if(strcmp(zonei->file,scripti->cval)==0){
      readzone(i,LOAD,&errorcode);
      return;
    }
  }
  for(i=0;i<nplot3dinfo;i++){
    plot3ddata *plot3di;

    plot3di = plot3dinfo + i;
    if(strcmp(plot3di->file,scripti->cval)==0){
      ReadPlot3dFile=1;
      readplot3d(plot3di->file,i,LOAD,&errorcode);
      update_menu();
      return;
    }
  }

  fprintf(stderr,"*** Error: file %s failed to load\n",scripti->cval);
}

/* ------------------ script_label ------------------------ */

void script_label(scriptdata *scripti){

  FREEMEMORY(script_labelstring);
  if(scripti->cval!=NULL&&strlen(scripti->cval)>0){
    NewMemory((void **)&script_labelstring,strlen(scripti->cval)+1);
    strcpy(script_labelstring,scripti->cval);
    PRINTF("*******************************\n");
    PRINTF("*** %s ***\n",script_labelstring);
    PRINTF("*******************************\n");
  }
}

/* ------------------ script_loadplot3d ------------------------ */

void script_loadplot3d(scriptdata *scripti){
  int i;
  float time_local;
  int blocknum;
  int count=0;

  time_local = scripti->fval;
  blocknum = scripti->ival-1;

  for(i=0;i<nplot3dinfo;i++){
    plot3ddata *plot3di;

    plot3di = plot3dinfo + i;
    if(plot3di->blocknumber==blocknum&&ABS(plot3di->time-time_local)<0.5){
      count++;
      LoadPlot3dMenu(i);
    }
  }
  UpdateRGBColors(COLORBAR_INDEX_NONE);
  set_labels_controls();
  if(count==0)fprintf(stderr,"*** Error: Plot3d file failed to load\n");

  //update_menu();
}

/* ------------------ script_loadvfile ------------------------ */

void script_loadvfile(scriptdata *scripti){
  int i;

  FREEMEMORY(loaded_file);
  PRINTF("script: loading vector slice file %s\n\n",scripti->cval);
  for(i=0;i<nvsliceinfo;i++){
    slicedata *val;
    vslicedata *vslicei;

    vslicei = vsliceinfo + i;
    val = sliceinfo + vslicei->ival;
    if(val==NULL)continue;
    if(strcmp(val->reg_file,scripti->cval)==0){
      LoadVSliceMenu(i);
      if(scripti->cval!=NULL&&strlen(scripti->cval)>0){
        NewMemory((void **)&loaded_file,strlen(scripti->cval)+1);
        strcpy(loaded_file,scripti->cval);
      }
      return;
    }
  }
  fprintf(stderr,"*** Error: Vector slice file %s was not loaded\n",scripti->cval);

}

/* ------------------ script_settourkeyframe ------------------------ */

void script_settourkeyframe(scriptdata *scripti){
  float keyframe_time;
  keyframe *keyj,*minkey=NULL;
  tourdata *touri;
  float minkeytime=1000000000.0;

  if(selected_tour==NULL)return;
  touri=selected_tour;
  keyframe_time=scripti->fval;
  for(keyj=(touri->first_frame).next;keyj->next!=NULL;keyj=keyj->next){
    float diff_time;

    if(keyj==(touri->first_frame).next){
      minkey=keyj;
      minkeytime=ABS(keyframe_time-keyj->nodeval.time);
      continue;
    }
    diff_time=ABS(keyframe_time-keyj->nodeval.time);
    if(diff_time<minkeytime){
      minkey=keyj;
      minkeytime=diff_time;
    }
  }
  if(minkey!=NULL){
    new_select(minkey);
    set_glui_keyframe();
    update_tourcontrols();
  }
}

/* ------------------ script_settourview ------------------------ */

void script_settourview(scriptdata *scripti){
  edittour=scripti->ival;
  show_tourlocus=scripti->ival3;
  tour_global_tension_flag=1;
  tour_global_tension=scripti->fval;
  switch(scripti->ival2){
    case 0:
      viewtourfrompath=0;
      keyframe_snap=0;
      break;
    case 1:
      viewtourfrompath=1;
      keyframe_snap=0;
      break;
    case 2:
      viewtourfrompath=0;
      keyframe_snap=1;
      break;
    default:
      viewtourfrompath=0;
      keyframe_snap=0;
      break;
  }
  update_tour_state();
}

/* ------------------ script_settimeval ------------------------ */

void script_settimeval(scriptdata *scripti){
  float timeval;
  int i,imin;
  float valmin;

  timeval = scripti->fval;
  PRINTF("script: setting time to %f\n\n",timeval);
  if(global_times!=NULL&&nglobal_times>0){
    if(timeval<global_times[0])timeval=global_times[0];
    if(timeval>global_times[nglobal_times-1]-0.0001){
      float dt;

      dt = timeval-global_times[nglobal_times-1]-0.0001;
      if(nglobal_times>1&&dt>global_times[1]-global_times[0]){
        fprintf(stderr,"*** Error: data not available at time requested\n");
        fprintf(stderr,"           time: %f s, min time: %f, max time: %f s\n",
          timeval,global_times[0],global_times[nglobal_times-1]);
        if(loaded_file!=NULL)fprintf(stderr,"           loaded file: %s\n",loaded_file);
        if(script_labelstring!=NULL)fprintf(stderr,"                 label: %s\n",script_labelstring);
      }
      timeval=global_times[nglobal_times-1]-0.0001;
    }
    valmin=ABS(global_times[0]-timeval);
    imin=0;
    for(i=1;i<nglobal_times-1;i++){
      float val;

      val = ABS(global_times[i]-timeval);
      if(val<valmin){
        valmin=val;
        imin=i;
      }
    }
    itimes=imin;
    script_itime=imin;
    stept=0;
    force_redisplay=1;
    UpdateFrameNumber(0);
    UpdateTimeLabels();
  }
}

//    sscanf(buffer,"%i %i %i %i",&vis_gslice_data, &show_gslice_triangles, &show_gslice_triangulation, &show_gslice_normal);
//    sscanf(buffer,"%f %f %f",gslice_xyz,gslice_xyz+1,gslice_xyz+2);
//    sscanf(buffer,"%f %f",gslice_normal_azelev,gslice_normal_azelev+1);

/* ------------------ script_gsliceview ------------------------ */

void script_gsliceview(scriptdata *scripti){
  vis_gslice_data = scripti->ival;
  show_gslice_triangles = scripti->ival2;
  show_gslice_triangulation = scripti->ival3;
  show_gslice_normal = scripti->ival4;
  update_gslice=1;
}

/* ------------------ script_gslicepos ------------------------ */

void script_gslicepos(scriptdata *scripti){
  gslice_xyz[0]=scripti->fval;
  gslice_xyz[1]=scripti->fval2;
  gslice_xyz[2]=scripti->fval3;
  update_gslice=1;
}

/* ------------------ script_gsliceorien ------------------------ */

void script_gsliceorien(scriptdata *scripti){
  gslice_normal_azelev[0]=scripti->fval;
  gslice_normal_azelev[1]=scripti->fval2;
  update_gslice=1;
}

/* ------------------ settimeval ------------------------ */

void settimeval(float timeval){
  int i;

  if(global_times!=NULL&&nglobal_times>0){
    if(timeval<global_times[0])timeval=global_times[0];
    if(timeval>global_times[nglobal_times-1]-0.0001)timeval=global_times[nglobal_times-1]-0.0001;
    for(i=0;i<nglobal_times;i++){
      float tlow, thigh;

      if(i==0){
        tlow = global_times[i];
        thigh = (global_times[i]+global_times[i+1])/2.0;
      }
      else if(i==nglobal_times-1){
        tlow = (global_times[i-1]+global_times[i])/2.0;
        thigh = global_times[i];
      }
      else{
        tlow=(global_times[i-1]+global_times[i])/2.0;
        thigh=(global_times[i]+global_times[i+1])/2.0;
      }
      if(tlow<=timeval&&timeval<thigh){
        itimes=i;
        stept=1;
        force_redisplay=1;
        UpdateFrameNumber(0);
        UpdateTimeLabels();
        keyboard('t',FROM_SMOKEVIEW);
        break;
      }
    }
  }
}

/* ------------------ script_setviewpoint ------------------------ */

void script_setviewpoint(scriptdata *scripti){
  char *viewpoint;
  cameradata *ca;
  int count=0;

  viewpoint = scripti->cval;
  PRINTF("script: set viewpoint to %s\n\n",viewpoint);
  for(ca=camera_list_first.next;ca->next!=NULL;ca=ca->next){
    if(strcmp(scripti->cval,ca->name)==0){
      ResetMenu(ca->view_id);
      count++;
      break;
    }
  }
  if(count==0)fprintf(stderr,"*** Error: The viewpoint %s was not found\n",viewpoint);
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
  PRINTF("\n");
  PRINTF("script: %s\n",scripti->command_label);
  if(scripti->cval!=NULL){
    PRINTF("script:  %s\n",scripti->cval);
  }
  if(scripti->cval2!=NULL){
    PRINTF("script:  %s\n",scripti->cval2);
  }
  PRINTF("\n");
  switch(scripti->command){
    case SCRIPT_UNLOADALL:
      LoadUnloadMenu(UNLOADALL);
      break;
    case SCRIPT_RENDERSIZE:
      script_render_width = scripti->ival;
      script_render_height = scripti->ival2;
      RenderMenu(RenderCustom);
      break;
    case SCRIPT_RENDERTYPE:
      if(STRCMP(scripti->cval, "JPG")==0){
        update_render_type(JPEG);
      }
      else{
        update_render_type(PNG);
      }
      break;
    case SCRIPT_MOVIETYPE:
      if(STRCMP(scripti->cval, "WMV") == 0){
        update_movie_type(WMV);
      }
      if(STRCMP(scripti->cval, "MP4") == 0){
        update_movie_type(MP4);
      }
      else{
        update_movie_type(AVI);
      }
      break;
    case SCRIPT_RENDERDIR:
      if(scripti->cval!=NULL&&strlen(scripti->cval)>0){
        script_dir_path=scripti->cval;
        if(can_write_to_dir(script_dir_path)==0){
          fprintf(stderr,"*** Error: Cannot write to the RENDERDIR directory: %s\n",script_dir_path);
        }
        PRINTF("script: setting render path to %s\n",script_dir_path);
      }
      else{
        script_dir_path=NULL;
      }
      break;
    case SCRIPT_KEYBOARD:
      {
        char *key;

        script_keystate=0;
        key = scripti->cval + strlen(scripti->cval) - 1;
        if(strncmp(scripti->cval,"ALT",3)==0)script_keystate=GLUT_ACTIVE_ALT;

        keyboard(*key,FROM_SCRIPT);
        returnval=1;
      }
      break;
    case SCRIPT_SCENECLIP:
      clip_mode=scripti->ival;
      updatefacelists=1;
      break;
    case SCRIPT_XSCENECLIP:
      clipinfo.clip_xmin=scripti->ival;
      clipinfo.xmin = scripti->fval;

      clipinfo.clip_xmax=scripti->ival2;
      clipinfo.xmax = scripti->fval2;
      updatefacelists=1;
      break;
    case SCRIPT_YSCENECLIP:
      clipinfo.clip_ymin=scripti->ival;
      clipinfo.ymin = scripti->fval;

      clipinfo.clip_ymax=scripti->ival2;
      clipinfo.ymax = scripti->fval2;
      updatefacelists=1;
      break;
    case SCRIPT_ZSCENECLIP:
      clipinfo.clip_zmin=scripti->ival;
      clipinfo.zmin = scripti->fval;

      clipinfo.clip_zmax=scripti->ival2;
      clipinfo.zmax = scripti->fval2;
      updatefacelists=1;
      break;
    case SCRIPT_RENDERCLIP:
      clip_rendered_scene=scripti->ival;
      render_clip_left=scripti->ival2;
      render_clip_right=scripti->ival3;
      render_clip_bottom=scripti->ival4;
      render_clip_top=scripti->ival5;
      break;
    case SCRIPT_RENDERONCE:
      keyboard('r',FROM_SMOKEVIEW);
      returnval=1;
      break;
    case SCRIPT_RENDERDOUBLEONCE:
      keyboard('R',FROM_SMOKEVIEW);
      returnval=1;
      break;
    case SCRIPT_RENDERSTART:
      script_renderstart(scripti);
      break;
    case SCRIPT_RENDERALL:
      script_renderall(scripti);
      break;
    case SCRIPT_RENDER360ALL:
      script_render360all(scripti);
      break;
    case SCRIPT_VOLSMOKERENDERALL:
      script_volsmokerenderall(scripti);
      break;
    case SCRIPT_ISORENDERALL:
      script_isorenderall(scripti);
      break;
    case SCRIPT_MAKEMOVIE:
      script_makemovie(scripti);
      break;
    case SCRIPT_LOADFILE:
      script_loadfile(scripti);
      break;
    case SCRIPT_LABEL:
      script_label(scripti);
      break;
    case SCRIPT_LOADINIFILE:
      script_loadinifile(scripti);
      break;
    case SCRIPT_LOADVFILE:
      script_loadvfile(scripti);
      break;
    case SCRIPT_LOADBOUNDARY:
      script_loadboundary(scripti,-1);
      break;
    case SCRIPT_LOADBOUNDARYM:
      script_loadboundary(scripti,scripti->ival);
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
      TourMenu(MENU_TOUR_MANUAL);
      break;
    case SCRIPT_EXIT:
#ifndef _DEBUG
      exit(0);
#endif
      break;
    case SCRIPT_LOADISO:
      script_loadiso(scripti,-1); // load isosurface for all meshes
      break;
    case SCRIPT_LOADISOM:
      script_loadiso(scripti, scripti->ival); // load isosurface for mesh script->ival
      break;
    case SCRIPT_LOAD3DSMOKE:
      script_load3dsmoke(scripti);
      break;
    case SCRIPT_LOADVOLSMOKE:
      script_loadvolsmoke(scripti);
      break;
    case SCRIPT_LOADVOLSMOKEFRAME:
      script_loadvolsmokeframe(scripti,1);
      returnval=1;
      break;
    case SCRIPT_LOADPARTICLES:
      script_loadparticles(scripti);
      break;
    case SCRIPT_LOADSLICE:
      script_loadslice(scripti);
      break;
    case SCRIPT_LOADSLICEM:
      script_loadslicem(scripti, scripti->ival2);
      break;
    case SCRIPT_LOADVSLICE:
      script_loadvslice(scripti);
      break;
    case SCRIPT_LOADVSLICEM:
      script_loadvslicem(scripti,scripti->ival2);
      break;
    case SCRIPT_LOADPLOT3D:
      script_loadplot3d(scripti);
      break;
    case SCRIPT_SETTIMEVAL:
      returnval=1;
      script_settimeval(scripti);
      break;
    case SCRIPT_SETTOURVIEW:
      script_settourview(scripti);
      break;
    case SCRIPT_SETTOURKEYFRAME:
      script_settourkeyframe(scripti);
      break;
    case SCRIPT_SETVIEWPOINT:
      script_setviewpoint(scripti);
      break;
    case SCRIPT_GSLICEVIEW:
      script_gsliceview(scripti);
      break;
    case SCRIPT_GSLICEPOS:
      script_gslicepos(scripti);
      break;
    case SCRIPT_GSLICEORIEN:
      script_gsliceorien(scripti);
      break;
    case SCRIPT_CBARFLIP:
      colorbarflip=0;
      ColorBarMenu(COLORBAR_FLIP);
      break;
    case SCRIPT_CBARNORMAL:
      colorbarflip=1;
      ColorBarMenu(COLORBAR_FLIP);
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
  glutPostRedisplay();
  return returnval;
}
