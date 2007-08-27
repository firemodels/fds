#include "options.h"
#ifdef pp_HRR
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

void stripcommas(char *buffer){
  int i;
  char *c;

  for(i=0;i<strlen(buffer);i++){
    c=buffer+i;
    if(*c==',')*c=' ';
  }
}

/* ------------------ printhrr ------------------------ */

void printhrr(void){
  int i;
  float *hrrtime, *hrrval;

  if(hrrinfo==NULL){
    printf("hrr data not available\n");
    return;
  }
  if(hrrinfo->ntimes==0){
    printf("hrr data not loaded\n");
    return;
  }
  hrrtime=hrrinfo->times;
  hrrval=hrrinfo->hrrval;
  for(i=0;i<hrrinfo->ntimes;i++){
    printf(" time=%f hrr=%f\n",*hrrtime,*hrrval);
    hrrtime++;
    hrrval++;
  }
}

/* ------------------ readhrr ------------------------ */

void readhrr(int flag, int *errorcode){
  FILE *HRRFILE;
  int ntimes;
  char buffer[1024];
  float *hrrtime, *hrrval;

  *errorcode=0;
  if(hrrinfo!=NULL){
    FREEMEMORY(hrrinfo->times);
    FREEMEMORY(hrrinfo->hrrval);
    FREEMEMORY(hrrinfo->timeslist);
  }
  FREEMEMORY(hrrinfo);
  if(flag==UNLOAD)return;

  NewMemory((void **)&hrrinfo,sizeof(hrrdata));
  hrrinfo->file=hrrfilename;
  hrrinfo->times=NULL;
  hrrinfo->timeslist=NULL;
  hrrinfo->hrrval=NULL;
  hrrinfo->ntimes=0;
  hrrinfo->loaded=1;
  hrrinfo->display=0;
  hrrinfo->itime=0;

  HRRFILE=fopen(hrrinfo->file,"r");
  if(HRRFILE==NULL){
    readhrr(UNLOAD,errorcode);
    return;
  }

// read first three lines

  if(fgets(buffer,1024,HRRFILE)==NULL)*errorcode=1;
  if(*errorcode==0){
    if(fgets(buffer,1024,HRRFILE)==NULL)*errorcode=1;
  }
  if(*errorcode==0){
    if(fgets(buffer,1024,HRRFILE)==NULL)*errorcode=1;
  }
  if(*errorcode==1){
    readhrr(UNLOAD,errorcode);
    *errorcode=1;
    return;
  }

// size data

  ntimes=0;
  while(!feof(HRRFILE)){
    if(fgets(buffer,1024,HRRFILE)==NULL)break;
    ntimes++;
  }

  rewind(HRRFILE);
  NewMemory((void **)&hrrinfo->times,ntimes*sizeof(float));
  NewMemory((void **)&hrrinfo->hrrval,ntimes*sizeof(float));

// read data
  
  hrrtime=hrrinfo->times;
  hrrval=hrrinfo->hrrval;
  ntimes=0;

  if(*errorcode==0&&fgets(buffer,1024,HRRFILE)==NULL)*errorcode=1;
  if(*errorcode==0&&fgets(buffer,1024,HRRFILE)==NULL)*errorcode=1;
  if(*errorcode==0&&fgets(buffer,1024,HRRFILE)==NULL)*errorcode=1;
  if(*errorcode==1){
    readhrr(UNLOAD,errorcode);
    *errorcode=1;
    return;
  }
  while(!feof(HRRFILE)){
    if(fgets(buffer,10245,HRRFILE)==NULL)break;
    stripcommas(buffer);
    sscanf(buffer,"%f %f",hrrtime,hrrval);
    hrrtime++;
    hrrval++;
    ntimes++;
  }
  hrrinfo->ntimes=ntimes;
}

/* ------------------ drawhrr ------------------------ */

void drawhrr(void){
}
#endif
