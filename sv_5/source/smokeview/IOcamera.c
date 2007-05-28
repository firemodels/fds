#include "options.h"
#ifdef pp_CAMERA
#include <stdio.h>  
#include <stdlib.h>
#include <math.h>
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

/* ------------------ init_camview ------------------------ */

void init_camview(camdata *cv, char *file, char *label){
  int len;

  cv->ncamviews=0;
  cv->camviews=NULL;
  cv->file=NULL;
  cv->time=NULL;
  if(file!=NULL){
    trim(file);
    len=strlen(file);
    if(len>0){
      NewMemory((void **)&cv->file,len+1);
      STRCPY(cv->file,file);
    }
  }
  cv->label=NULL;
  if(label!=NULL){
    trim(label);
    len=strlen(label);
    if(len>0){
      NewMemory((void **)&cv->label,len+1);
      STRCPY(cv->label,label);
    }
  }
}

/* ------------------ free_camview ------------------------ */

void free_camview(camdata *cv){
  FREEMEMORY(cv->camviews);
  FREEMEMORY(cv->file);
  FREEMEMORY(cv->label);
  FREEMEMORY(cv->time);
  cv->ncamviews=0;
}

/* ------------------ read_camview ------------------------ */

void read_camview(int ifile, int flag, int *errorcode){
  FILE *stream;
  char buffer[255];
  int n;
  float *time, *exyz0, *vxyz0, *ap;
  camviewdata *cvd;
  camdata *cv;

  *errorcode=1;
  if(ifile<0||ifile>=ncaminfo)return;

  cv=caminfo+ifile;
  cv->ncamviews=0;
  cv->available=0;
 
  stream=fopen(cv->file,"r");
  if(stream==NULL){
    *errorcode=1;
    return;
  }
  n=0;
  while(!feof(stream)){
    if(fgets(buffer,255,stream)==NULL)break;
    n++;
  }
  if(n==0)return;
 
  cv->ncamviews=n;
  NewMemory((void **)&cv->camviews,n*sizeof(camviewdata));
  NewMemory((void **)&cv->time,n*sizeof(float));
  rewind(stream);

  n=0;
  while(!feof(stream)){
    if(fgets(buffer,255,stream)==NULL)break;
    cvd = cv->camviews+n;
    time=&cvd->time;
    exyz0=cvd->eye0;
    vxyz0=cvd->view0;
    ap=&cvd->aperture;
    sscanf(buffer,"%f %f %f %f %f %f %f %f",time, exyz0,exyz0+1,exyz0+2,vxyz0,vxyz0+1,vxyz0+2,ap);
    cv->time[n]=*time;
    cvd->eye[0]=(cvd->eye0[0]-xbar0)/xyzmaxdiff;
    cvd->eye[1]=(cvd->eye0[1]-ybar0)/xyzmaxdiff;
    cvd->eye[2]=(cvd->eye0[2]-zbar0)/xyzmaxdiff;
    cvd->view[0]=(cvd->view0[0]-xbar0)/xyzmaxdiff;
    cvd->view[1]=(cvd->view0[1]-ybar0)/xyzmaxdiff;
    cvd->view[2]=(cvd->view0[2]-zbar0)/xyzmaxdiff;
    n++;
  }
  cv->available=1;
}


#endif
