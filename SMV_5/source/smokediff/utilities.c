// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "svdiff.h"
#include "svn_revision.h"
#include "MALLOC.h"

// svn revision character string
char utilities_revision[]="$Revision$";

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

  /* ------------------ getfileinfo ------------------------ */

int getfileinfo(char *filename, char *source_dir, int *filesize){
  STRUCTSTAT statbuffer;
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
  statfile=STAT(buffer,&statbuffer);
  if(statfile!=0)return statfile;
  if(filesize!=NULL)*filesize=statbuffer.st_size;
  return statfile;
}

/* ------------------ version ------------------------ */

void version(void){
    char smv_version[100];
    int svn_num;

    getSMDiffversion(smv_version);  // get Smokeview version (ie 5.x.z)
    svn_num=getmaxrevision();    // get svn revision number
    printf("\n");
    printf("Smokediff\n\n");
    printf("Version: %s\n",smv_version);
    printf("SVN Revision Number: %i\n",svn_num);
    printf("Compile Date: %s\n",__DATE__);
#ifdef X64
    printf("Platform: WIN64\n");
#endif
#ifdef WIN32
#ifndef X64
    printf("Platform: WIN32\n");
#endif
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

  MAXREV(main_revision);
  MAXREV(utilities_revision);
  return max_revision;
}

/* ------------------ getSMVversion ------------------------ */

void getSMDiffversion(char *SMDiffversion){
  strcpy(SMDiffversion,SMDiffVERSION);
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

/* ------------------ fullfile ------------------------ */

void fullfile(char *fileout, char *dir, char *file){
  strcpy(fileout,"");
  if(dir!=NULL)strcat(fileout,dir);
  strcat(fileout,file);
}

/* ------------------ getendian ------------------------ */

int getendian(void){
  short val;
  char *cval;
  val=1;
  cval = (char *)&val+1;
  return (int)(*cval);
}


