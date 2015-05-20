// $Date: 2012-08-20 19:39:01 -0400 (Mon, 20 Aug 2012) $ 
// $Revision: 12156 $
// $Author: koverholt $

#include "options.h"
#include <stdio.h>
#include <ctype.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "svn_revision.h"
#include "blockaid.h"
#include "datadefs.h"
#include "ASSERT.h"

// svn revision character string
char utilities_revision[]="$Revision: 12156 $";

/* ------------------ usage ------------------------ */

void usage(void){
  char group_version[100];
  int svn_num;

  getPROGversion(group_version);  // get Blockaid version (ie 5.x.z)
  svn_num=getmaxrevision();    // get svn revision number

    printf("\n");
    printf("Usage:\n");
    printf(" blockaid [-f] [-h] [-l libpath] casename \n");
    printf("          -f - overwrite output file casename.fds\n");
    printf("          -h - display this message\n");
    printf("          -l libdir - optional directory path containing &INCL files\n\n");
    printf("  This program reads in casename.fof and outputs casename.fds where specified\n");
    printf("  groups of blockages, holes and vents are replicated, translated and rotated.\n");
    printf("  These groups are surrounded with &BGRP and &EGRP as in:\n");
    printf("      &BGRP ID='group label' ORIGIN=x,y,z /\n");
    printf("          one or more &OBST, &VENT, &HOLE \n");
    printf("      &EGRP /\n");
    printf("  where ORIGIN is the group's origin.  This group may be replicated using:\n");
    printf("     &GRP ID='group label' XYZ=x,y,z ROTATE=angle /\n");
    printf("  where x,y,z is the translation amount and angle is the rotation amount\n");
    printf("  snapped to the nearest 90 degrees.\n");
    printf("  See the blockaid Wiki at http://fire.nist.gov/fds for more information.\n");
}

/* ------------------ version ------------------------ */

void version(void){
    char blockaid_version[100];
    int svn_num;

    getPROGversion(blockaid_version);  // get Sblockaid verson (ie x,y,z)
    svn_num=getmaxrevision();    // get svn revision number
    printf("\n");
    printf("Blockaid\n\n");
    printf("Version: %s\n",blockaid_version);
    printf("SVN Revision Number: %i\n",svn_num);
    printf("Compile Date: %s\n",__DATE__);
#ifdef WIN32
    printf("Platform: WIN32\n");
#endif
#ifdef pp_OSX
    printf("Platform: OS X\n");
#endif
#ifdef pp_LINUX
    printf("Platform: LINUX\n");
#endif
    printf("\nUse blockaid -h for more information\n");
}

/* ------------------ getmaxrevision ------------------------ */

int getmaxrevision(void){
  int max_revision=0, rev;

  MAXREV(assert_revision);
  MAXREV(dmalloc_revision);
  MAXREV(main_revision);
  MAXREV(readfds_revision);
  MAXREV(utilities_revision);
  return max_revision;
}

/* ------------------ getrevision ------------------------ */

int getrevision(char *svn){
  char svn_string[256];
  char *svn_ptr;
  int return_val;

  svn_ptr=svn_string;
  svn=strchr(svn,':');
  if(svn==NULL||strlen(svn)<=4)return 0;
  
  svn++;
  strcpy(svn_ptr,svn);
  svn_ptr=trim_front(svn_ptr);
  svn_ptr[strlen(svn_ptr)-1]=0;
  trim(svn_ptr);
  sscanf(svn_ptr,"%i",&return_val);
  return return_val;

}

/* ------------------ trim_front ------------------------ */

char *trim_front(char *line){
  size_t i;

  for(i=0;i<strlen(line);i++){
    if(line[i]!=' ')return line+i;
  }
  return line;
}

/* ------------------ trim ------------------------ */

void trim(char *line){
  char *blank=" ";
  const char *c;
  const char *lf="\n";
  size_t len, i;

  len = strlen(line);
  c = line+len-1;
  for(i=0; i<len; i++){
    if(strncmp(c,blank,1)!=0&&strncmp(c,lf,1)!=0){
      c++; 
      line[c-line]='\0';
      return;
    }
    c--;
  }
  *line='\0';
}

  /* ------------------ getfileinfo ------------------------ */

int getfileinfo(char *filename, char *source_dir, int *filesize){
  struct stat statbuffer;
  int statfile;
  char buffer[MAXLINE];

  if(source_dir==NULL){
    strcpy(buffer,filename);
  }
  else{
    strcpy(buffer,source_dir);
    strcat(buffer,filename);
  }
  if(filesize!=NULL)*filesize=0;
  statfile=stat(buffer,&statbuffer);
  if(statfile!=0)return statfile;
  if(filesize!=NULL)*filesize=statbuffer.st_size;
  return statfile;
}

/* ------------------ match ------------------------ */

int match(const char *buffer, const char *key, unsigned int lenkey){
  if(strncmp(buffer,key,lenkey) == 0)return(1);
  return(0);
}

/* ------------------ getkeyparam ------------------------ */

char *get_keyid(char *source, const char *key){
  char *keyptr,*s1,*e1;

// given:  key='xxxyyy'
// returns  location of beginning of xxxyyy

  if(source==NULL||key==NULL||strlen(key)==0)return NULL;
  keyptr=strstr(source,key);
  if(keyptr==NULL)return NULL;
  s1=strstr(keyptr,"'");
  if(s1==NULL)return NULL;
  s1++;
  e1=strstr(s1,"'");
  if(e1==NULL)return NULL;
  e1[0]='\0';
  return s1;
}


/* ------------------ parseobst_xb ------------------------ */

int get_irvals(char *line, char *key, int nvals, int *ivals, float *rvals, int *ibeg, int *iend){
  char *keyptr, *keystart, *xslash, *c;
  char *cvals;
  char line2[10000];
  size_t len;
  int mode,token;
  size_t i;
  int exitloop;
#define INBLANK 0
#define INTOKEN 1

  if(ibeg!=NULL)*ibeg=-1;
  if(iend!=NULL)*iend=-2;
  if(line==NULL||strlen(line)<1)return 0;
  if(key==NULL||strlen(key)<1)return 0;
  if(nvals<1)return 0;

  len=strlen(line);
  if(len+1>10000)return 0;

  for(i=0;i<len;i++){
    line2[i]=toupper(line[i]);
  }
  line2[len]='\0';

  keyptr=NULL;
  keystart=line2;
  while(keyptr==NULL){
    keyptr = strstr(keystart,key);
    if(keyptr==NULL)return 0;
    if(!isspace(keyptr[-1])&&keyptr[-1]!=','){
      keystart=keyptr+1;
      keyptr=NULL;
    }
  }

  xslash = strstr(line2,"/");
  if(xslash==NULL||keyptr>xslash)return 0;

  mode = INBLANK;
  token=0;
  cvals=NULL;
  exitloop=0;
  for(c=keyptr+strlen(key);*c!='\0';c++){
    if(exitloop==1)break;
    switch (mode) {
    case INBLANK:
      if(*c==' '||*c==','||*c=='\n'){
        *c=' ';
        break;
      }
      if(token==0&&*c=='='){
        cvals=c+1;
        if(ibeg!=NULL)*ibeg=(int)(cvals-line2);
        break;
      }
      token++;
      mode = INTOKEN;
      break;
    case INTOKEN:
      if(*c==' '||*c==','||*c=='/'){
        if(token==nvals){
          exitloop=1;
//          *c='\0';
          if(iend!=NULL)*iend=(int)(c-line2);
          break;
        }
        if(*c=='/'){
          return 0;
        }
        *c=' ';
        mode=INBLANK;
      }
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
  }
  if(token!=nvals)return 0;
  if(rvals!=NULL){
    switch (nvals){
      case 1:
        sscanf(cvals,"%f",rvals);
        break;
      case 3:
        sscanf(cvals,"%f %f %f",rvals,rvals+1,rvals+2);
        break;
      case 6:
        sscanf(cvals,"%f %f %f %f %f %f",rvals,rvals+1,rvals+2,rvals+3,rvals+4,rvals+5);
        break;
      default:
        return 0;
        break;
    }
  }
  if(ivals!=NULL){
    switch (nvals){
      case 1:
        sscanf(cvals,"%i",ivals);
        break;
      case 3:
        sscanf(cvals,"%i %i %i",ivals,ivals+1,ivals+2);
        break;
      case 6:
        sscanf(cvals,"%i %i %i %i %i %i",ivals,ivals+1,ivals+2,ivals+3,ivals+4,ivals+5);
        break;
      default:
        return 0;
        break;
    }
  }
  return nvals;
}


/* ------------------ trimzeros ------------------------ */

void trimzeros(char *line){
  size_t i,len;
  char *c;

  len = strlen(line);
  c = line + len-1;
  for(i=len-1;i>0;i--){
    if(*c=='0'){
      c--;
      if(*c=='.'){
        line[i+1]='\0';
        return;
      }
      continue;
    }
    line[i+1]='\0';
    return;
  }
  line[0]='\0';
}

/* ------------------ trimmzeros ------------------------ */

void trimmzeros(char *line){
  size_t i,j,lenline;
  char *linecopy, c;
  char buffer[MAXLINE];
  size_t ibeg;
  size_t start=1;
  size_t bufstart=0;
  size_t lenbuf;

  lenline=strlen(line);
  linecopy=line;
  for(i=0;i<lenline;i++){
    c=*linecopy++;
    if(start==1){
      if(c==' ')continue;
      ibeg=i;
      start=0;
      buffer[0]=c;
      continue;
    }
    buffer[i-ibeg]=c;
    if(c==' '){
      buffer[i-ibeg]='\0';
      trimzeros(buffer);
      lenbuf = strlen(buffer);
      for(j=bufstart;j<bufstart+lenbuf;j++){
        line[j]=buffer[j-bufstart];
      }
      bufstart+=strlen(buffer);
      line[bufstart]=' ';
      bufstart++;
      start=1;
    }
  }
  line[bufstart]='\0';
}

/* ------------------ rotatexy ------------------------ */

  void rotatexy(float *dx, float *dy, float rotate, float *dxy){
    int irotate;
    float dxy2[2];

    irotate = (int)(rotate/90.0+0.5);
    irotate = irotate%4;

    switch (irotate){
      case 0:
        dxy2[0]=*dx;
        dxy2[1]=*dy;
        break;
      case 1:
        dxy2[0]=-(*dy);
        dxy2[1]= *dx;
        if(dxy!=NULL){
          dxy2[0]+=dxy[1];
        }
        break;
      case 2:
        dxy2[0]=-(*dx);
        dxy2[1]=-(*dy);
        if(dxy!=NULL){
          dxy2[0]+=dxy[0];
          dxy2[1]+=dxy[1];
        }
        break;
      case 3:
        dxy2[0]= *dy;
        dxy2[1]=-(*dx);
        if(dxy!=NULL){
          dxy2[1]+=dxy[0];
        }
        break;
    }
    *dx=dxy2[0];
    *dy=dxy2[1];
  }

/* ------------------ float2string ------------------------ */

void float2string(float *xb, int nxb, char *xbstring){
  char charxb[100];
  int i;

  strcpy(xbstring,"");
  for(i=0;i<nxb;i++){
    sprintf(charxb,"%f",xb[i]);
    trimzeros(charxb);
    strcat(xbstring,charxb);
    strcat(xbstring,", ");
  }
}

/* ------------------ int2string ------------------------ */

void int2string(int *ib, int nib, char *ibstring){
  char charib[100];
  int i;

  strcpy(ibstring,"");
  for(i=0;i<nib;i++){

    sprintf(charib,"%i",ib[i]);
    trimzeros(charib);
    strcat(ibstring,charib);
    if(i!=nib-1){
      strcat(ibstring,", ");
    }
  }
}

/* ------------------ subst_string ------------------------ */

void subst_string(char *string, int ibeg, int iend, char *replace){
  char buffer[MAXLINE];
  int ii,i;

  if(string==NULL)return;
  ii=0;
  for(i=0;i<ibeg;i++){
    buffer[ii++]=string[i];
  }
  if(replace!=NULL){
    for(i=0;i<strlen(replace);i++){
      buffer[ii++]=replace[i];
    }
  }
  for(i=iend+1;i<strlen(string);i++){
    buffer[ii++]=string[i];
  }
  buffer[ii]='\0';
  strcpy(string,buffer);
}

