// $Date: 2008-03-05 16:31:28 -0500 (Wed, 05 Mar 2008) $ 
// $Revision: 1409 $
// $Author: gforney $

#include "options.h"
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "transpose.h"
#include "MALLOC.h"
#include "svn_revision.h"

// svn revision character string
char utilities_revision[]="$Revision: 1409 $";

/* ------------------ match ------------------------ */

int match(const char *buffer, const char *key, unsigned int lenkey){
  if(strncmp(buffer,key,lenkey) == 0)return(1);
  return(0);
}

/* ------------------ trim_front ------------------------ */

char *trim_front(char *line){
  char *blank=" ";
  const char *c;
  unsigned int i,len;

  c = line;
  len=strlen(line);
  for(i=0;i<len;i++){
    if(strncmp(c++,blank,1)!=0)return line+i;
  }
  return line;
}

/* ------------------ trim ------------------------ */

void trim(char *line){
  char *blank=" ";
  const char *c;
  const char *lf="\n";
  unsigned int len;
  unsigned int i;
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

/* ------------------ trim_both ------------------------ */

char *trim_both(char *line){
  char *return_val;

  trim(line);
  return_val = trim_front(line);
  return return_val;
}

/* ------------------ readlabels ------------------------ */

int readlabels(flowlabels *flowlabel, FILE *stream){
  char buffer[255];
  unsigned int len;

  flowlabel->longlabel=NULL;
  flowlabel->shortlabel=NULL;
  flowlabel->unit=NULL;

  if(fgets(buffer,255,stream)==NULL)return 2;

  len=strlen(buffer);
  buffer[len-1]='\0';
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void **)&flowlabel->longlabel,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->longlabel,buffer);


  if(fgets(buffer,255,stream)==NULL)return 2;
  len=strlen(buffer);
  buffer[len-1]='\0';
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void **)&flowlabel->shortlabel,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->shortlabel,buffer);

  if(fgets(buffer,255,stream)==NULL)return 2;
  len=strlen(buffer);
  buffer[len-1]='\0';
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void *)&flowlabel->unit,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->unit,buffer);
  return 0;
}

  /* ------------------ trimzeros ------------------------ */

void trimzeros(char *line){
  unsigned int i;
  unsigned int len;
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

  /* ------------------ getfileinfo ------------------------ */

int getfileinfo(char *filename, char *source_dir, int *filesize){
  struct stat statbuffer;
  int statfile;
  char buffer[1024];

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

  /* ------------------ getfilesizelabel ------------------------ */

void getfilesizelabel(int size, char *sizelabel){
  int leftsize,rightsize;

#define sizeGB   1000000000
#define size100MB 100000000
#define size10MB   10000000
#define sizeMB     1000000
#define size100KB    100000
#define size10KB      10000

  if(size>=sizeGB){
    size/=size10MB;
    leftsize=size/100;
    rightsize=size-100*leftsize;
    sprintf(sizelabel,"%i.%02i GB",leftsize,rightsize);
  }
  else if(size>=size100MB&&size<sizeGB){
    size/=sizeMB;
    leftsize=size;
    sprintf(sizelabel,"%i MB",leftsize);
  }
  else if(size>=size10MB&&size<size100MB){
    size/=size100KB;
    leftsize=size/10;
    rightsize=size-10*leftsize;
    sprintf(sizelabel,"%i.%i MB",leftsize,rightsize);
  }
  else if(size>=sizeMB&&size<size10MB){
    size/=size10KB;
    leftsize=size/100;
    rightsize=size-100*leftsize;
    sprintf(sizelabel,"%i.%02i MB",leftsize,rightsize);
  }
  else if(size>=size100KB&&size<sizeMB){
    size/=1000;
    leftsize=size;
    sprintf(sizelabel,"%i KB",leftsize);
  }
  else{
    size/=10;
    leftsize=size/100;
    rightsize=size-100*leftsize;
    sprintf(sizelabel,"%i.%02i KB",leftsize,rightsize);
  }
}
/* ------------------ version ------------------------ */

void version(void){
    char spr_version[100];
    int svn_num;

    getSPRversion(spr_version);  // get Smokeview version (ie 5.x.z)
    svn_num=getmaxrevision();    // get svn revision number
    printf("\n");
    printf("Spreadpose\n\n");
    printf("Version: %s\n",spr_version);
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

}

/* ------------------ getmaxrev ------------------------ */

#define MAXREV(cval) max_revision=imax(getrevision(cval),max_revision)
int getmaxrevision(void){
  int max_revision=0;

  MAXREV(assert_revision);
  MAXREV(dmalloc_revision);
  MAXREV(main_revision);
  MAXREV(utilities_revision);
  MAXREV(csv_convert_revision);

  return max_revision;
}

/* ------------------ getSMVversion ------------------------ */

void getSPRversion(char *SPRversion){
  strcpy(SPRversion,SPRVERSION);
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

/* ------------------ imax ------------------------ */

int imax(int a, int b){
  if(a>b){
    return a;
  }
  else{
    return b;
  }
}
