#include "options.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "flowfiles.h"
#include "MALLOC.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"


#define EXPMIN -2
#define EXPMAX 3


/* ------------------ getBoundaryColors ------------------------ */

void getBoundaryColors(float *t, int nt, unsigned char *it, 
              int settmin, float *ttmin, int settmax, float *ttmax,
              float *tmin_global, float *tmax_global,
              int ndatalevel, int nlevel,
              char **labels, char *scale, float *tvals256){
  int n;
  float *tcopy, factor, tval, range;
  int expmin, expmax;
  int itt;
  float local_tmin, local_tmax, tmin2, tmax2;
  int local_skip;

  tmin2 = *t;
  tmax2 = *t;
    
  STRCPY(scale,"");
  tcopy = t+1;
  for(n=1;n<nt;n++){
    if(*tcopy<tmin2)tmin2=*tcopy;
    if(*tcopy>tmax2)tmax2=*tcopy;
    tcopy++;
  }
  *tmin_global = tmin2;
  *tmax_global = tmax2;
  local_skip=0;
  adjustdatabounds(t,local_skip,nt,settmin,&tmin2,settmax,&tmax2);
  if(settmin!=SET_MIN){
    *ttmin=tmin2;
  }
  if(settmax!=SET_MAX){
    *ttmax=tmax2;
  }
  local_tmin = *ttmin;
  local_tmax = *ttmax;

  range = local_tmax - local_tmin;
  factor = 0.0f;
  if(range!=0.0f)factor = ndatalevel/range;
  for(n=0;n<nt;n++){
    itt=(int)(factor*(*t++-local_tmin));
    if(itt<0)itt=0;
    if(itt>ndatalevel-1)itt=ndatalevel-1;
    *it=itt;
    it++;
  }
  frexp10(local_tmax, &expmax);
  frexp10(local_tmin, &expmin);
  if(expmin!=0&&expmax!=0&&expmax-expmin<=2&&(expmin<-2||expmin>2)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(scale,"*10^%i",expmin);
  }
  if(expmin==0&&(expmax<EXPMIN||expmax>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmax);
    local_tmax *= pow((double)10.0,(double)-expmax);
    sprintf(scale,"*10^%i",expmax);
  }
  if(expmax==0&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(scale,"*10^%i",expmin);
  }
  range = local_tmax - local_tmin;
  factor = range/(nlevel-2);
  for (n=1;n<nlevel-2;n++){
    tval = local_tmin + (n-1)*factor;
    num2string(&labels[n][0],tval,range);
  }
  tval = local_tmin + (nlevel-3)*factor;
  for(n=0;n<256;n++){
    tvals256[n] = (local_tmin*(255-n) + n*local_tmax)/255.;
  }
  num2string(&labels[nlevel-2][0],tval,range);
  tval = local_tmax;
  num2string(&labels[nlevel-1][0],tval,range);
}

/* ------------------ getBoundaryColors2 ------------------------ */

void getBoundaryColors2(float *t, int nt, unsigned char *it, 
              int settmin, float *ttmin, int settmax, float *ttmax,
              float *tmin_global, float *tmax_global,
              int ndatalevel){
  int n;
  float *tcopy, factor, range;
  int itt;
  float local_tmin, local_tmax, tmin2, tmax2;
  int local_skip;

  tmin2 = *t;
  tmax2 = *t;
    
  tcopy = t+1;
  for(n=1;n<nt;n++){
    if(*tcopy<tmin2)tmin2=*tcopy;
    if(*tcopy>tmax2)tmax2=*tcopy;
    tcopy++;
  }
  *tmin_global = tmin2;
  *tmax_global = tmax2;
  local_skip=0;
  adjustdatabounds(t,local_skip,nt,settmin,&tmin2,settmax,&tmax2);
  if(settmin!=SET_MIN){
    *ttmin=tmin2;
  }
  if(settmax!=SET_MAX){
    *ttmax=tmax2;
  }
  local_tmin = *ttmin;
  local_tmax = *ttmax;

  range = local_tmax - local_tmin;
  factor = 0.0f;
  if(range!=0.0f)factor = ndatalevel/range;
  for(n=0;n<nt;n++){
    itt=(int)(factor*(*t++-local_tmin));
    if(itt<0)itt=0;
    if(itt>ndatalevel-1)itt=ndatalevel-1;
    *it=itt;
    it++;
  }
}

/* ------------------ getBoundaryLabels ------------------------ */

void getBoundaryLabels(
              float local_tmin, float local_tmax,
              char **labels, char *scale, float *tvals256, int nlevel){
  int n;
  float factor, tval, range;
  int expmin, expmax;

  STRCPY(scale,"");

  frexp10(local_tmax, &expmax);
  frexp10(local_tmin, &expmin);
  if(expmin!=0&&expmax!=0&&expmax-expmin<=2&&(expmin<-2||expmin>2)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(scale,"*10^%i",expmin);
  }
  if(expmin==0&&(expmax<EXPMIN||expmax>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmax);
    local_tmax *= pow((double)10.0,(double)-expmax);
    sprintf(scale,"*10^%i",expmax);
  }
  if(expmax==0&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(scale,"*10^%i",expmin);
  }
  range = local_tmax - local_tmin;
  factor = range/(nlevel-2);
  for (n=1;n<nlevel-2;n++){
    tval = local_tmin + (n-1)*factor;
    num2string(&labels[n][0],tval,range);
  }
  tval = local_tmin + (nlevel-3)*factor;
  for(n=0;n<256;n++){
    tvals256[n] = (local_tmin*(255-n) + n*local_tmax)/255.;
  }
  num2string(&labels[nlevel-2][0],tval,range);
  tval = local_tmax;
  num2string(&labels[nlevel-1][0],tval,range);
}


/* ------------------ getPart5Colors ------------------------ */

void getPart5Colors(particle *parti, int nlevel){
  int i,j,k,m;
  part5data *datacopy;

  datacopy = parti->data5;
  for(i=0;i<parti->nframes;i++){
    for(j=0;j<parti->nclasses;j++){
      float valmin, valmax, dval;
      part5class *partclassi;
      float *rvals;
      unsigned char *irvals;

      partclassi = partclassinfo + j;
      rvals = datacopy->rvals;
      irvals = datacopy->irvals;

      for(k=2;k<partclassi->ntypes;k++){
        part5prop *prop_id;

        prop_id = get_part5prop(partclassi->labels[k].longlabel);
        if(prop_id==NULL)continue;

        valmin = prop_id->valmin;
        valmax = prop_id->valmax;
        dval = valmax - valmin;
        if(dval<=0.0)dval=1.0;

        for(m=0;m<datacopy->npoints;m++){
          float val;
          int irval;

          val=*rvals++;
          irval = 255*(val-valmin)/dval;
          if(irval<0)irval=0;
          if(irval>255)irval=255;
          *irvals++=irval;

        }
      }
      FREEMEMORY(datacopy->rvals);
      datacopy++;
    }
  }
  for(i=0;i<npart5prop;i++){
    int n;

    part5prop *propi;
    float local_tmin, local_tmax;
    int expmin, expmax;
    float factor,range,tval;
    char *scale,**labels;
    float *ppartlevels256;

    propi = part5propinfo + i;

    local_tmin = propi->valmin;
    local_tmax = propi->valmax;
    scale = propi->scale;
    labels=propi->partlabels;
    ppartlevels256=propi->ppartlevels256;

    frexp10(local_tmax, &expmax);
    frexp10(local_tmin, &expmin);
    if(expmin!=0&&expmax!=0&&expmax-expmin<=2&&(expmin<EXPMIN||expmin>EXPMAX)){
      local_tmin *= pow((double)10.0,(double)-expmin);
      local_tmax *= pow((double)10.0,(double)-expmin);
      sprintf(scale,"*10^%i",expmin);
    }
    if(expmin==0&&(expmax<EXPMIN||expmax>EXPMAX)){
      local_tmin *= pow((double)10.0,(double)-expmax);
      local_tmax *= pow((double)10.0,(double)-expmax);
      sprintf(scale,"*10^%i",expmax);
    }
    if(expmax==0&&(expmin<EXPMIN||expmin>EXPMAX)){
      local_tmin *= pow((double)10.0,(double)-expmin);
      local_tmax *= pow((double)10.0,(double)-expmin);
      sprintf(scale,"*10^%i",expmin);
    }
    range = local_tmax - local_tmin;

    factor = range/(nlevel-2);
    for (n=1;n<nlevel-2;n++){
      tval = local_tmin + (n-1)*factor;
      num2string(&labels[n][0],tval,range);
    }
    for(n=0;n<256;n++){
      ppartlevels256[n] = (local_tmin*(255-n) + n*local_tmax)/255.;
    }
    tval = local_tmin + (nlevel-3)*factor;
    num2string(&labels[nlevel-2][0],tval,range);
    tval = local_tmax;
    num2string(&labels[nlevel-1][0],tval,range);
    CheckMemory;
  }

}
/* ------------------ getPartColors ------------------------ */

void getPartColors(const float *t, int local_skip, int nt, 
                   unsigned char *it,const unsigned char *isprink, int particle_type, int droplet_type, 
              const float *ttmin, const float *ttmax, int nlevel,
              char **labels, char *scale, float *ppartlevels256){
  int n;
  float factor, tval, range;
  float local_tmin, local_tmax;
  int expmax, expmin;
  int itt;
  const unsigned char *isprinkcopy;
  const float *tcopy;
  unsigned char *itcopy;

  isprinkcopy=isprink;
  itcopy = it;
  tcopy = t;

  STRCPY(scale,"");
  t += local_skip;
  it += local_skip;
  isprink += local_skip;
  if(particle_type==0){
    for (n=local_skip;n<nt;n++){
      if(*isprink==0){  /* its not a sprinkler drop so color it according to below ... */
        if(     *t>-0.1&&*t<0.1 ){
          *it=rgb_white;
        }
        else if(*t>0.9 &&*t<1.1 ){
          *it=rgb_yellow;
        }
        else if(*t>1.9&&*t<2.1){
          *it=rgb_blue;
        }
        else if(*t>2.9&&*t<3.1){
          *it=rgb_red;
        }
        else if(*t>3.9&&*t<4.1){
          *it=rgb_green;
        }
        else if(*t>4.9&&*t<5.1){
          *it=rgb_magenta;
        }
        else if(*t>5.9&&*t<6.1){
          *it=rgb_cyan;
        }
        else if(*t>6.9&&*t<7.1){
          *it=rgb_black;
        }
        else{
           *it=rgb_white;
         }
      }
      else {
        *it=rgb_blue;
      }  /* its a sprinkler drop (*isprink=1) so color it blue */
      t++;
      it++;
      isprink++;
    }
    STRCPY(&labels[1][0],"0.0");
    for (n=2;n<nlevel-1;n++){
      STRCPY(&labels[n][0]," ");
    }
    STRCPY(&labels[nlevel-1][0],"1.0");
    if(particle_type==0&&droplet_type==0)return;
  }

  isprink=isprinkcopy;
  it=itcopy;
  t=tcopy;

  local_tmin=*ttmin;
  local_tmax=*ttmax;
  range = local_tmax - local_tmin;
  factor=0.0f;
  if(range!=0.0f){
    factor = 255./range;
  }
  for(n=local_skip;n<nt;n++){
    if((particle_type==1&&*isprink==0)||(*isprink==1&&droplet_type==1)){
      itt=(int)(factor*(*t-local_tmin));
      if(itt<0)itt=0;
      if(itt>255){
        itt=255;
      }
      *it = itt;
    }
    if(droplet_type==0&&*isprink==1){
      *it=rgb_blue;
    }
    it++;
    t++;
    isprink++;
  }

  frexp10(local_tmax, &expmax);
  frexp10(local_tmin, &expmin);
  if(expmin!=0&&expmax!=0&&expmax-expmin<=2&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(scale,"*10^%i",expmin);
  }
  if(expmin==0&&(expmax<EXPMIN||expmax>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmax);
    local_tmax *= pow((double)10.0,(double)-expmax);
    sprintf(scale,"*10^%i",expmax);
  }
  if(expmax==0&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(scale,"*10^%i",expmin);
  }
  range = local_tmax - local_tmin;

  factor = range/(nlevel-2);
  for (n=1;n<nlevel-2;n++){
    tval = local_tmin + (n-1)*factor;
    num2string(&labels[n][0],tval,range);
  }
  for(n=0;n<256;n++){
    ppartlevels256[n] = (local_tmin*(255-n) + n*local_tmax)/255.;
  }
  tval = local_tmin + (nlevel-3)*factor;
  num2string(&labels[nlevel-2][0],tval,range);
  tval = local_tmax;
  num2string(&labels[nlevel-1][0],tval,range);
}

/* ------------------ getZoneColor ------------------------ */

int getZoneColor(float t, float local_tmin, float local_tmax, int nlevel){
  int level;

  if(t<=local_tmin)return 0;
  if(t>=local_tmax)return nlevel-1;
  if(local_tmin==local_tmax)return 0;
  level=nlevel*(t-local_tmin)/(local_tmax-local_tmin);
  if(level<0)return 0;
  if(level>nlevel-1)return nlevel-1;
  return level;
}


/* ------------------ getZoneColors ------------------------ */

void getZoneColors(const float *t, int nt, unsigned char *it,
               const float *ttmin, const float *ttmax, int nlevel, int nlevel_full,
               char **labels, char *scale, float *tvals256
               ){
  int n;
  float dt, factor;
  int itt;
  int expmin, expmax;
  float local_tmin, local_tmax;
  float range;
  float tval;

  local_tmin = *ttmin;
  local_tmax = *ttmax;

  dt = local_tmax - local_tmin;
  factor=0.0f;
  if(dt!=0.0f)factor = nlevel_full/dt;
  for(n=0;n<nt;n++){
    itt=(int)(factor*(*t-local_tmin));
    if(itt<0)itt=0;
    if(itt>nlevel_full-1)itt=nlevel_full-1;
    *it=itt;
    it++;t++;
	/* iblank++;*/

  }

  STRCPY(scale,"");

  frexp10(local_tmax, &expmax);
  frexp10(local_tmin, &expmin);
  if(expmin!=0&&expmax!=0&&expmax-expmin<=2&&(expmin<-2||expmin>2)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(scale,"*10^%i",expmin);
  }
  if(expmin==0&&(expmax<EXPMIN||expmax>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmax);
    local_tmax *= pow((double)10.0,(double)-expmax);
    sprintf(scale,"*10^%i",expmax);
  }
  if(expmax==0&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(scale,"*10^%i",expmin);
  }
  range = local_tmax - local_tmin;
  factor = range/(nlevel-2);
  for (n=1;n<nlevel-2;n++){
    tval = local_tmin + (n-1)*factor;
    num2string(&labels[n][0],tval,range);
  }
  tval = local_tmin + (nlevel-3)*factor;
  for(n=0;n<256;n++){
    tvals256[n] = (local_tmin*(255-n) + n*local_tmax)/255.;
  }
  num2string(&labels[nlevel-2][0],tval,range);
  tval = local_tmax;
  num2string(&labels[nlevel-1][0],tval,range);

}

/* ------------------ getDataColors ------------------------ */

void getPlot3DColors(int plot3dvar, int settmin, float *ttmin, int settmax, float *ttmax, 
              int ndatalevel, int nlevel,
              char **labels,char **scale, float *tlevels, float *tlevels256){
  int n;
  float dt, factor, tval;
  float local_tmin, local_tmax, tmin2, tmax2;
  float range;
  int expmax,expmin;
  float tminorig, tmaxorig, dtorig;
  int itt;
  float *q;
  int *iq;
  plot3d *p;
  mesh *meshi;
  int *iblank,i;
  int ntotal;

  tmin2=*ttmin;
  tmax2=*ttmax;

  tmin2= 1000000000.;
  tmax2=-1000000000.;

  for(i=0;i<nplot3d;i++){
    p = plot3dinfo+i;
    if(p->loaded==0||p->display==0)continue;
    meshi = selected_case->meshinfo+p->blocknumber;
    ntotal=(meshi->ibar+1)*(meshi->jbar+1)*(meshi->kbar+1);
    iblank=meshi->iblank;
    q=meshi->qdata+plot3dvar*ntotal;
    for(n=0;n<ntotal;n++){
      if(*iblank++==1){
        if(*q<tmin2)tmin2=*q;
        if(*q>tmax2)tmax2=*q;
      }
      q++;
    }
  }
  if(settmin==PERCENTILE_MIN||settmin==GLOBAL_MIN){
    local_tmin=tmin2;
  }
  else{
    local_tmin=*ttmin;
  }
  if(settmax==PERCENTILE_MAX||settmax==GLOBAL_MAX){
    local_tmax=tmax2;
  }
  else{
    local_tmax=*ttmax;
  }
  adjustPlot3Dbounds(plot3dvar,settmin,&local_tmin,settmax,&local_tmax);

  *ttmin=local_tmin;
  *ttmax=local_tmax;
  range = local_tmax-local_tmin;
  tminorig=local_tmin;
  tmaxorig=local_tmax;
  if(range!=0.0f){
    factor = (float)ndatalevel/range;
  }
  else{
    factor = 0.0f;
  }

  for(i=0;i<nplot3d;i++){
    p = plot3dinfo+i;
    if(p->loaded==0||p->display==0)continue;
    meshi = selected_case->meshinfo+p->blocknumber;
    ntotal=(meshi->ibar+1)*(meshi->jbar+1)*(meshi->kbar+1);
    q=meshi->qdata+plot3dvar*ntotal;
    iq=meshi->iqdata+plot3dvar*ntotal;
    for(n=0;n<ntotal;n++){
      itt=(int)(factor*(*q++-local_tmin));
      if(itt<0)itt=0;
      if(itt>ndatalevel-1)itt=ndatalevel-1;
      *iq++=itt;
    }
  }

  STRCPY(*scale,"");
  frexp10(local_tmax, &expmax);
  frexp10(local_tmin, &expmin);
  if(expmin!=0&&expmax!=0&&expmax-expmin<=2&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(*scale,"*10^%i",expmin);
  }
  if(expmin==0&&(expmax<EXPMIN||expmax>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmax);
    local_tmax *= pow((double)10.0,(double)-expmax);
    sprintf(*scale,"*10^%i",expmax);
  }
  if(expmax==0&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(*scale,"*10^%i",expmin);
  }

  range = local_tmax-local_tmin;
  dt = range/(float)(nlevel-2);
  dtorig = (tmaxorig-tminorig)/(float)(nlevel-2);
  for (n=1;n<nlevel-1;n++){
    tval = local_tmin + (n-1)*dt;
    num2string(&labels[n][0],tval,range);
  }
  for(n=0;n<256;n++){
    tlevels256[n] = (local_tmin*(255-n) + local_tmax*n)/255.;
  }
  tval = local_tmax;
  num2string(&labels[nlevel-1][0],tval,range);

  for(n=0;n<nlevel+1;n++){tlevels[n]=tminorig+(float)(n-1)*dtorig;}
}

/* ------------------ getSliceColors ------------------------ */

void getSliceColors(const float *t, int nt, unsigned char *it,
              float local_tmin, float local_tmax, 
              int ndatalevel, int nlevel,
              char labels[12][11],char **scale, float *tlevels256){
  int n;
  float dt, factor, tval;
  float range;
  int expmax,expmin;
  int itt;

  range = local_tmax-local_tmin;
  if(range!=0.0f){
    factor = (float)ndatalevel/range;
  }
   else{
     factor = 0.0f;
   }
  for(n=0;n<nt;n++){
    itt=(int)(factor*(*t++-local_tmin));
    if(itt<0)itt=0;
    if(itt>ndatalevel-1)itt=ndatalevel-1;
    *it=itt;
    it++;
  }

  STRCPY(*scale,"");
  frexp10(local_tmax, &expmax);
  frexp10(local_tmin, &expmin);
  if(expmin!=0&&expmax!=0&&expmax-expmin<=2&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(*scale,"*10^%i",expmin);
  }
  if(expmin==0&&(expmax<EXPMIN||expmax>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmax);
    local_tmax *= pow((double)10.0,(double)-expmax);
    sprintf(*scale,"*10^%i",expmax);
  }
  if(expmax==0&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(*scale,"*10^%i",expmin);
  }

  range = local_tmax-local_tmin;
  dt = range/(float)(nlevel-2);
  for (n=1;n<nlevel-1;n++){
    tval = local_tmin + (n-1)*dt;
    num2string(&labels[n][0],tval,range);
  }
  for(n=0;n<256;n++){
    tlevels256[n] = (local_tmin*(255-n) + local_tmax*n)/255.;
  }
  tval = local_tmax;
  num2string(&labels[nlevel-1][0],tval,range);
}

/* ------------------ getSliceLabelels ------------------------ */

void getSliceLabels(float local_tmin, float local_tmax, int nlevel,
              char labels[12][11],char **scale, float *tlevels256){
  int n;
  float dt, tval;
  float range;
  int expmax,expmin;

  range = local_tmax-local_tmin;

  STRCPY(*scale,"");
  frexp10(local_tmax, &expmax);
  frexp10(local_tmin, &expmin);
  if(expmin!=0&&expmax!=0&&expmax-expmin<=2&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(*scale,"*10^%i",expmin);
  }
  if(expmin==0&&(expmax<EXPMIN||expmax>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmax);
    local_tmax *= pow((double)10.0,(double)-expmax);
    sprintf(*scale,"*10^%i",expmax);
  }
  if(expmax==0&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(*scale,"*10^%i",expmin);
  }

  range = local_tmax-local_tmin;
  dt = range/(float)(nlevel-2);
  for (n=1;n<nlevel-1;n++){
    tval = local_tmin + (n-1)*dt;
    num2string(&labels[n][0],tval,range);
  }
  for(n=0;n<256;n++){
    tlevels256[n] = (local_tmin*(255-n) + local_tmax*n)/255.;
  }
  tval = local_tmax;
  num2string(&labels[nlevel-1][0],tval,range);
}

/* ------------------ frexp10 ------------------------ */

float frexp10(float x, int *exp10){
  float xabs, mantissa;

  xabs = fabs((double)x);
  if(x==0.0f){
    *exp10=0;
    return 0.0f;
  }
  mantissa = log10((double)xabs);
  *exp10 = (int)floor((double)mantissa);
      
  mantissa = pow((double)10.0f,(double)mantissa-(double)*exp10);
  if(x<0)mantissa = -mantissa;
  return mantissa;
}

/* ------------------ drawColorBars ------------------------ */

void drawColorBars(void){
  int i,i2,i3;
  float labeltop;
  float barleft;
  float dbar=.10f, dtext=.3f, dspace=.05f, space;
  float barbot=-0.5f;
  float right[3],left[3],bottom[4];
  int ileft=0;
  int leftzone, leftsmoke, leftslice, leftpatch;
  float yy,yy2;
  int iposition;
  float position;
  float tttval, tttmin, tttmax;
  float val;
  char slicelabel[256], p3dlabel[256], boundarylabel[256], partlabel[256], zonelabel[256];
  char unitlabel[256];
  float *p3lev;
  databounds *sb;
  patch *patchi;

  int sliceunitclass,sliceunittype;
  int sliceflag=0;
  float *slicefactor=NULL;
  float slicerange;
  char slicecolorlabel[256];
  char *slicecolorlabel_ptr=NULL;

  int plot3dunitclass, plot3dunittype;
  int plot3dflag=0;
  float *plot3dfactor=NULL;
  float plot3drange;
  char plot3dcolorlabel[256];
  char *plot3dcolorlabel_ptr=NULL;

  int patchunitclass, patchunittype;
  int zoneunitclass, zoneunittype;
  int patchflag=0;
  int zoneflag=0;
  float *patchfactor=NULL;
  float *zonefactor=NULL;
  float patchrange=0.0;
  float zonerange;
  char patchcolorlabel[256];
  char *patchcolorlabel_ptr=NULL;

  char zonecolorlabel[256];
  char *zonecolorlabel_ptr=NULL;

  int partunitclass, partunittype;
  int partflag=0;
  float *partfactor=NULL;
  float partrange=0.0;
  char partcolorlabel[256];
  char *partcolorlabel_ptr=NULL;

  char partunitlabel2[256], partshortlabel2[256];

  GLfloat *color1, *color2;
  color1=&(foregroundcolor[0]);
  color2=&(redcolor[0]);

  barleft = barright-dbar;
  right[0] = barleft;
  space = dtext + dspace;
  if(fontindex==LARGE_FONT)space *= 1.5;
  left[0] = barleft - space;
  if(left[0]<0.0f)left[0]=0.0;

  right[1]=left[0];
  left[1] = left[0] - space;
  if(left[1]<0.0f)left[1]=0.0;

  right[2]=left[1];
  left[2] = left[1]-space;
  if(left[2]<0.0f)left[2]=0.0;

  labeltop=nrgb+.7;
  switch (fontindex){
  case SMALL_FONT:
    bottom[0]=labeltop-0.0;
    bottom[1]=labeltop-0.4;
    bottom[2]=labeltop-0.8;
    bottom[3]=labeltop-1.2;
    break;
  case LARGE_FONT:
    bottom[0]=labeltop-0.2;
    bottom[1]=labeltop-0.8;
    bottom[2]=labeltop-1.4;
    bottom[3]=labeltop-2.0;
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  
  if((showsmoke==1&&parttype!=0)||showslice==1||showvslice==1||showpatch==1||(showzone==1&&sethazardcolor==0)||showplot3d==1){
    sniffErrors("before colorbar");
    CheckMemory;
    glBegin(GL_QUADS);
    if(showplot3d==1&&p3cont2d!=SHADED_CONTOURS){
      {
        {
          int icolmin, icolmax;

          for (i = 0; i < nrgb; i++){
            i2 = i+1;i3 = i+1;
            icolmin = i*((nrgb_full-1)/(nrgb-1));
            icolmax = i*((nrgb_full-1)/(nrgb-1));
            if(i3>nrgb_full-1)i3=nrgb_full;

            if(showplot3d==1&&(i==0&&minfill==0||i==nrgb-1&&maxfill==0)){}
             else{
              //glColor4fv(rgb[i]); 
              glColor4fv(rgb_full[icolmin]);
              glVertex2f(barleft, (float)i+barbot); 
              glVertex2f(barright,(float)i+barbot);
       
              //glColor4fv(rgb[i3]);
              glColor4fv(rgb_full[icolmax]);
              glVertex2f(barright,(float)i2+barbot);
              glVertex2f(barleft, (float)i2+barbot);
             }
          }
        }
      }
    }
    else{
      for (i = 0; i < nrgb_full; i++){
        yy = (barbot*(255-i)+i*(nrgb+barbot))/255.;
        yy2 = (barbot*(254-i)+(i+1)*(nrgb+barbot))/255.;
        i3=i+1;
        if(i==nrgb_full-1)i3=i;
        glColor4fv(rgb_full[i]); 
        glVertex2f(barleft, yy);
        glVertex2f(barright,yy);

        glColor4fv(rgb_full[i3]);
        glVertex2f(barright,yy2);
        glVertex2f(barleft,yy2);
      }
    }
    /*  code to put in chop triangle
    glColor3f(0.0,1.0,0.0);
    glVertex2f(barleft-0.25,1.0);
    glVertex2f(barleft,1.0);

    glVertex2f(barleft-0.125,1.125);
    glVertex2f(barleft-0.125,1.125);
    */
    glEnd();
  }
  leftsmoke=0;
  leftslice=0;
  leftpatch=0;
  leftzone=0;
  if(showsmoke==1){
    if(parttype!=0){
      leftsmoke=ileft;
      ileft++;
    }
  }
  if(showslice==1||showvslice==1){
    leftslice=ileft;
    ileft++;
  }
  if(showpatch==1){
    leftpatch=ileft;
  }
  strcpy(partshortlabel2,"");
  strcpy(partunitlabel2,"");
  if(showsmoke==1){
    if(parttype!=0)outputBarText(right[leftsmoke],bottom[0],color1,"Part");
    if(parttype==-1){
      strcpy(partshortlabel2,"temp");
      strcpy(partunitlabel2,"C");
    }
    else if(parttype==-2){
      strcpy(partshortlabel2,"HRRPUV");
      strcpy(partunitlabel2,"kW/m^3");
    }
    else{
      if(partshortlabel!=NULL)strcpy(partshortlabel2,partshortlabel);
      if(partunitlabel!=NULL)strcpy(partunitlabel2,partunitlabel);
    }
    if(parttype!=0){
      getunitinfo(partshortlabel2,&partunitclass,&partunittype);
      if(partunitclass>=0&&partunitclass<nunitclasses){
        if(partunittype>=0){
          partflag=1;
          partfactor=unitclasses[partunitclass].units[partunittype].scale;
          strcpy(partunitlabel,unitclasses[partunitclass].units[partunittype].unit);
        }
      }
      outputBarText(right[leftsmoke],bottom[1],color1,partshortlabel);
      outputBarText(right[leftsmoke],bottom[2],color1,partunitlabel);
      outputBarText(right[leftsmoke],bottom[3],color1,partscale);
    }
  }
  if(showslice==1||showvslice==1){
    sb = slicebounds + islicetype;
    strcpy(unitlabel,sb->label->unit);
    getunitinfo(sb->label->shortlabel,&sliceunitclass,&sliceunittype);
    if(sliceunitclass>=0&&sliceunitclass<nunitclasses){
      if(sliceunittype>0){
        sliceflag=1;
        slicefactor=unitclasses[sliceunitclass].units[sliceunittype].scale;
        strcpy(unitlabel,unitclasses[sliceunitclass].units[sliceunittype].unit);
      }
    }
    outputBarText(right[leftslice],bottom[0],color1,"Slice");
    outputBarText(right[leftslice],bottom[1],color1,sb->label->shortlabel);
    outputBarText(right[leftslice],bottom[2],color1,unitlabel);
    outputBarText(right[leftslice],bottom[3],color1,sb->scale);
  }
  if(showpatch==1){
    patchi = patchinfo + patchtypes[ipatchtype];
    strcpy(unitlabel,patchi->label.unit);
    getunitinfo(patchi->label.shortlabel,&patchunitclass,&patchunittype);
    if(patchunitclass>=0&&patchunitclass<nunitclasses){
      if(patchunittype>0){
        patchflag=1;
        patchfactor=unitclasses[patchunitclass].units[patchunittype].scale;
        strcpy(unitlabel,unitclasses[patchunitclass].units[patchunittype].unit);
      }
    }
    outputBarText(right[leftpatch],bottom[0],color1,"Bndry");
    outputBarText(right[leftpatch],bottom[1],color1,patchi->label.shortlabel);
    outputBarText(right[leftpatch],bottom[2],color1,unitlabel);
    outputBarText(right[leftpatch],bottom[3],color1,patchi->scale);
  }
  if(showplot3d==1){
    strcpy(unitlabel,unitp3label[plotn-1]);
    getunitinfo(shortp3label[plotn-1],&plot3dunitclass,&plot3dunittype);
    if(plot3dunitclass>=0&&plot3dunitclass<nunitclasses){
      if(plot3dunittype>0){
        plot3dflag=1;
        plot3dfactor=unitclasses[plot3dunitclass].units[plot3dunittype].scale;
        strcpy(unitlabel,unitclasses[plot3dunitclass].units[plot3dunittype].unit);
      }
    }
    outputBarText(right[0],bottom[0],color1,"Plot3d");
    outputBarText(right[0],bottom[1],color1,shortp3label[plotn-1]);
    outputBarText(right[0],bottom[2],color1,unitlabel);
    outputBarText(right[0],bottom[3],color1,scalep3[plotn-1]);
  }
  if(showzone==1&&sethazardcolor==0){
    strcpy(unitlabel,"C");
    getunitinfo("Temp",&zoneunitclass,&zoneunittype);
    if(zoneunitclass>=0&&zoneunitclass<nunitclasses){
      if(zoneunittype>0){
        zoneflag=1;
        zonefactor=unitclasses[zoneunitclass].units[zoneunittype].scale;
        strcpy(unitlabel,unitclasses[zoneunitclass].units[zoneunittype].unit);
      }
    }
    outputBarText(right[leftzone],bottom[0],color1,"Zone");
    outputBarText(right[leftzone],bottom[1],color1,"Temp");
    outputBarText(right[leftzone],bottom[2],color1,unitlabel);
    outputBarText(right[leftzone],bottom[3],color1,zonescale);

  }
  if(showsmoke==1&&parttype!=0){
    float *partlevels256_ptr;

    partlevels256_ptr=partlevels256;
    if(prop_index>=0&&prop_index<npart5prop){
      partlevels256_ptr=part5propinfo[prop_index].ppartlevels256;
    }

    iposition=-1;
    tttmin = partlevels256_ptr[0];
    tttmax = partlevels256_ptr[255];
    partrange = tttmax - tttmin;
    if(global_changecolorindex!=-1){
      tttval = partlevels256_ptr[valindex];
      num2string(partlabel,tttval,partrange);
      partcolorlabel_ptr=partlabel;
      if(partflag==1){
        scalefloat2string(tttval,partcolorlabel, partfactor, partrange);
        partcolorlabel_ptr=partcolorlabel;
      }
      position = (float)global_changecolorindex/255.*nrgb-0.5;
      iposition=(int)position;
      outputBarText(right[leftsmoke],position,color2,partcolorlabel_ptr);
    }
    for (i=1; i<nrgb; i++){
      if(iposition==i-1)continue;
      if(prop_index>=0&&prop_index<npart5prop){
        partcolorlabel_ptr=&part5propinfo[prop_index].partlabels[i][0];
      }
      else{
        partcolorlabel_ptr=&colorlabelpart[i][0];
      }
      if(partflag==1){
        val = tttmin + (i-1)*partrange/(nrgb-2);
        scalefloat2string(val,partcolorlabel, partfactor, partrange);
        scalestring(partcolorlabel_ptr,partcolorlabel, partfactor, partrange);
        partcolorlabel_ptr=partcolorlabel;
      }
      outputBarText(right[leftsmoke],(float)i+barbot,color1,partcolorlabel_ptr);
    }
  }
  if(showslice==1||showvslice==1){
    sb=slicebounds+islicetype;
    tttmin = sb->levels256[0];
    tttmax = sb->levels256[255];
    slicerange=tttmax-tttmin;
    iposition=-1;
    if(global_changecolorindex!=-1){
      tttval = sb->levels256[valindex];
      num2string(slicelabel,tttval,slicerange);
      slicecolorlabel_ptr=slicelabel;
      if(sliceflag==1){
        scalefloat2string(tttval,slicecolorlabel, slicefactor, slicerange);
        slicecolorlabel_ptr=slicecolorlabel;
      }
      position = (float)global_changecolorindex/255.*nrgb-0.5;
      iposition=(int)position;
      outputBarText(right[leftslice],position,color2,slicecolorlabel_ptr);
    }
    for (i=1; i<nrgb; i++){
      if(iposition==i-1)continue;
      slicecolorlabel_ptr=&(sb->colorlabels[i][0]);
      if(sliceflag==1){
        val = tttmin + (i-1)*slicerange/(nrgb-2);
        scalefloat2string(val,slicecolorlabel, slicefactor, slicerange);
        slicecolorlabel_ptr=slicecolorlabel;
      }
      outputBarText(right[leftslice],(float)i+barbot,color1,slicecolorlabel_ptr);
    }
  }
  if(showpatch==1){
    iposition=-1;
    tttmin = boundarylevels256[0];
    tttmax = boundarylevels256[255];
    patchrange=tttmax-tttmin;
    if(global_changecolorindex!=-1){
      tttval = boundarylevels256[valindex];
      num2string(boundarylabel,tttval,tttmax-tttmin);
      position = (float)global_changecolorindex/255.*nrgb-0.5;
      iposition=(int)position;
      patchcolorlabel_ptr=&(boundarylabel[0]);
      if(patchflag==1){
        scalefloat2string(tttval,patchcolorlabel, patchfactor, patchrange);
        patchcolorlabel_ptr=patchcolorlabel;
      }
      outputBarText(right[leftpatch],position,color2,patchcolorlabel_ptr);
    }
    for (i=1; i<nrgb; i++){
      if(iposition==i-1)continue;
      patchcolorlabel_ptr=&colorlabelpatch[i][0];
      if(patchflag==1){
        val = tttmin + (i-1)*patchrange/(nrgb-2);
        scalefloat2string(val,patchcolorlabel, patchfactor, patchrange);
        patchcolorlabel_ptr=patchcolorlabel;
      }
      outputBarText(right[leftpatch],(float)i+barbot,color1,patchcolorlabel_ptr);
    }
  }
  if(showzone==1&&sethazardcolor==0){
    iposition=-1;
    tttmin = zonelevels256[0];
    tttmax = zonelevels256[255];
    zonerange=tttmax-tttmin;
    if(global_changecolorindex!=-1){
      tttval = zonelevels256[valindex];
      num2string(zonelabel,tttval,tttmax-tttmin);
      position = (float)global_changecolorindex/255.*nrgb-0.5;
      iposition=(int)position;
      zonecolorlabel_ptr=&(zonelabel[0]);
      if(zoneflag==1){
        scalefloat2string(tttval,zonecolorlabel, zonefactor, zonerange);
        zonecolorlabel_ptr=zonecolorlabel;
      }
      outputBarText(right[leftzone],position,color2,zonecolorlabel_ptr);
    }
    for (i=1; i<nrgb; i++){
      if(iposition==i-1)continue;
      zonecolorlabel_ptr=&colorlabelzone[i][0];
      if(zoneflag==1){
        val = tttmin + (i-1)*zonerange/(nrgb-2);
        scalefloat2string(val,zonecolorlabel, zonefactor, zonerange);
        zonecolorlabel_ptr=zonecolorlabel;
      }
      outputBarText(right[leftzone],(float)i+barbot,color1,zonecolorlabel_ptr);
    }
  }
  if(showplot3d==1){
    iposition=-1;
    p3lev = p3levels256[plotn-1];
    tttmin = p3lev[0];
    tttmax = p3lev[255];
    plot3drange = tttmax - tttmin;
    if(global_changecolorindex!=-1){
      tttval = p3lev[valindex];
      num2string(p3dlabel,tttval,tttmax-tttmin);
      plot3dcolorlabel_ptr = p3dlabel;
      if(plot3dflag==1){
        scalefloat2string(tttval,plot3dcolorlabel, plot3dfactor, plot3drange);
        plot3dcolorlabel_ptr=plot3dcolorlabel;
      }
      position = (float)global_changecolorindex/255.*nrgb-0.5;
      iposition=(int)position;
      outputBarText(right[0],position,color2,plot3dcolorlabel_ptr);
    }
    if(visiso==0){
      for (i=1; i<nrgb; i++){
        if(iposition==i-1)continue;
        plot3dcolorlabel_ptr=&colorlabelp3[plotn-1][i][0];
        if(plot3dflag==1){
          val = tttmin + (i-1)*plot3drange/(nrgb-2);
          scalefloat2string(val,plot3dcolorlabel, plot3dfactor, plot3drange);
          plot3dcolorlabel_ptr=plot3dcolorlabel;
        }
        outputBarText(right[0],(float)i+barbot,color1,plot3dcolorlabel_ptr);
      }
    }
    else
    {
      for (i=1; i<nrgb; i++){
        if(iposition==i-1)continue;
        plot3dcolorlabel_ptr=&colorlabelp3[plotn-1][i][0];
        if(plot3dflag==1){
          val = tttmin + (i-1)*plot3drange/(nrgb-2);
          scalefloat2string(val,plot3dcolorlabel, plot3dfactor, plot3drange);
          plot3dcolorlabel_ptr=plot3dcolorlabel;
        }
        if(isolevelindex+1==i||isolevelindex2+1==i){
          outputBarText(right[0],(float)i+barbot,color2,plot3dcolorlabel_ptr);
        }
        else{
          outputBarText(right[0],(float)i+barbot,color1,plot3dcolorlabel_ptr);
        }
      }
    }
  }


}

/* ------------------ initcolors ------------------------ */

void initcolors(){
  mat_specular_orig[0]=0.5f;
  mat_specular_orig[1]=0.5f;
  mat_specular_orig[2]=0.2f;
  mat_specular_orig[3]=1.0f;

  mat_ambient_orig[0] = 0.5f;
  mat_ambient_orig[1] = 0.5f;
  mat_ambient_orig[2] = 0.2f;
  mat_ambient_orig[3] = 1.0f;

  ventcolor_orig[0]=1.0;
  ventcolor_orig[1]=0.0;
  ventcolor_orig[2]=1.0;
  ventcolor_orig[3]=1.0;

  block_ambient_orig[0] = 1.0;
  block_ambient_orig[1] = 0.8;
  block_ambient_orig[2] = 0.4;
  block_ambient_orig[3] = 1.0;

  ventcolor=getcolorptr(ventcolor_orig);
  block_ambient2=getcolorptr(block_ambient_orig);
  mat_ambient2=getcolorptr(mat_ambient_orig);
  mat_specular2=getcolorptr(mat_specular_orig);
}


/* ------------------ initcadcolors ------------------------ */

void initcadcolors(void){
  int n, i1, i2, i;
  float xx, f1, f2, sum;
  switch (setbw){
   case 0:
    for(n=0;n<nrgb_cad;n++){
      xx = (float)n/(float)nrgb_cad * (float)(nrgb-1);
      i1 = (int)xx;
      i2 = (int)(xx+1);
      f2 = xx - (float)i1;
      f1 = 1.0f - f2;
      sum=0;
      for(i=0;i<3;i++){
        rgb_cad[n][i] = f1*rgb[i1][i] + f2*rgb[i2][i];
        sum += rgb_cad[n][i]*rgb_cad[n][i];
      }
      sum=sqrt((double)sum);
      for(i=0;i<3;i++){
        rgb_cad[n][i] /= sum;
      }
      rgb_cad[n][3]=1.0;
    }
    break;
   case 1:
    for(n=0;n<nrgb_cad;n++){
      xx = (float)n/(float)nrgb_cad;
      for(i=0;i<3;i++){
        rgb_cad[n][i] = xx;
      }
      rgb_cad[n][3]=1.0;
    }
    break;
    default:
      ASSERT(FFALSE);
      break;
  }
}

/* ------------------ update_texturebar ------------------------ */

void update_texturebar(void){
  glBindTexture(GL_TEXTURE_1D,texture_colorbar_id);
  glTexImage1D(GL_TEXTURE_1D,0,4,256,0,GL_RGBA,GL_FLOAT,rgb_full);
  glBindTexture(GL_TEXTURE_1D,texture_slice_colorbar_id);
  glTexImage1D(GL_TEXTURE_1D,0,4,256,0,GL_RGBA,GL_FLOAT,rgb_slice);
  glBindTexture(GL_TEXTURE_1D,texture_plot3d_colorbar_id);
  glTexImage1D(GL_TEXTURE_1D,0,4,256,0,GL_RGBA,GL_FLOAT,rgb_plot3d);
}

/* ------------------ initrgb ------------------------ */

void initrgb(void){
  float *rgb_ptr;
  int n,nn;

  if(setbw==0){
    colorconvert(TO_COLOR);
    if(nrgb_ini !=0){
      nrgb = nrgb_ini;
      rgb_ptr=rgb_ini;
      for(n=0;n<nrgb_ini;n++){
        nn=(n+colorbarcycle)%nrgb_ini;
        rgb[n][0] = rgb_ptr[nn*3];
        rgb[n][1] = rgb_ptr[nn*3+1];
        rgb[n][2] = rgb_ptr[nn*3+2];
        rgb[n][3] = transparentlevel;
      }
    }   
    else{
      for(n=0;n<nrgb;n++){
        nn=(n+colorbarcycle)%nrgb;
        rgb[n][0] = rgb_base[nn][0];
        rgb[n][1] = rgb_base[nn][1];
        rgb[n][2] = rgb_base[nn][2];
        rgb[n][3] = transparentlevel;
      }
    }
  }
  else{
    colorconvert(TO_BW);
    for(n=0;n<nrgb;n++){
      nn=(n+colorbarcycle)%nrgb;
      rgb[n][0] = bw_base[nn][0];
      rgb[n][1] = bw_base[nn][1];
      rgb[n][2] = bw_base[nn][2];
      rgb[n][3] = transparentlevel;
    }
  }
}

/* ------------------ updatecolors ------------------------ */

void updatecolors(int changecolorindex){

  int n,nn;
  int i,j;
  float *rgb2ptr;
  int cci;
  mesh *meshi;
  int vent_offset, outline_offset;
  facedata *facej;

  initcadcolors();
  initrgb();
  nrgb_full = MAXRGB;
  for(n=0;n<nrgb_full;n++){
    rgb_trans[4*n+0]=0.0;
    rgb_trans[4*n+1]=0.0;
    rgb_trans[4*n+2]=0.0;
    rgb_trans[4*n+3]=(float)n/(float)(nrgb_full-1);
  }
  if(colorbarinfo!=NULL){
    for(n=0;n<nrgb_full;n++){
      rgb_full[n][0]=colorbarinfo[colorbartype].colorbar[3*n];
      rgb_full[n][1]=colorbarinfo[colorbartype].colorbar[3*n+1];
      rgb_full[n][2]=colorbarinfo[colorbartype].colorbar[3*n+2];
      rgb_full[n][3]=transparentlevel;
    } 
  }
  else{
    for(n=0;n<nrgb_full;n++){
      rgb_full[n][0]=(float)n/(float)(nrgb_full);
      rgb_full[n][1]=(float)n/(float)(nrgb_full);
      rgb_full[n][2]=(float)n/(float)(nrgb_full);
      rgb_full[n][3]=transparentlevel;
    } 
  }
  if(colorbarcycle!=0){
    {
      int icolor,nnn;

      for(n=0;n<nrgb_full;n++){
        rgb_full2[n][0]=rgb_full[n][0];
        rgb_full2[n][1]=rgb_full[n][1];
        rgb_full2[n][2]=rgb_full[n][2];
        rgb_full2[n][3]=rgb_full[n][3];
      }
      icolor=colorbarcycle*nrgb_full/nrgb;
      for(n=0;n<nrgb_full;n++){
        nnn=(n+icolor)%nrgb_full;
        rgb_full[nnn][0]=rgb_full2[n][0];
        rgb_full[nnn][1]=rgb_full2[n][1];
        rgb_full[nnn][2]=rgb_full2[n][2];
        rgb_full[nnn][3]=rgb_full2[n][3];
      }
    }
  }
  if(colorbarflip==1){
    {
      int nnn;

      for(n=0;n<nrgb_full;n++){
        rgb_full2[n][0]=rgb_full[n][0];
        rgb_full2[n][1]=rgb_full[n][1];
        rgb_full2[n][2]=rgb_full[n][2];
        rgb_full2[n][3]=rgb_full[n][3];
      }
      for(n=0;n<nrgb_full;n++){
        nnn=nrgb_full-1-n;
        rgb_full[nnn][0]=rgb_full2[n][0];
        rgb_full[nnn][1]=rgb_full2[n][1];
        rgb_full[nnn][2]=rgb_full2[n][2];
        rgb_full[nnn][3]=rgb_full2[n][3];
      }
    }
  }
  global_changecolorindex=changecolorindex;
  if(changecolorindex>=0){
    valindex = global_changecolorindex;
    valindex = (int)(6.*(valindex-255./12.)/5.+1);
    if(valindex<0)valindex=0;
    if(valindex>255)valindex=255;
    cci = changecolorindex;
    if(setbw==1){
      for(n=-5;n<6;n++){
        if(cci+n>255)break;
        if(cci+n<0)continue;
        rgb_full[cci+n][0]=1.;
        rgb_full[cci+n][1]=0.;
        rgb_full[cci+n][2]=0.;
        rgb_full[cci+n][3]=transparentlevel;
      }
    }
    else{
      for(n=-5;n<6;n++){
        if(cci+n>255)break;
        if(cci+n<0)continue;
        rgb_full[cci+n][0]=0.;
        rgb_full[cci+n][1]=0.;
        rgb_full[cci+n][2]=0.;
        rgb_full[cci+n][3]=transparentlevel;
      }
    }
  }
  if(rgb2_ini!=NULL){
    rgb2ptr=rgb2_ini;
  }
  else{
    rgb2ptr=&(rgb2[0][0]);
  }
  if(changecolorindex!=0){
    for(n=0;n<nrgb;n++){
      nn=n*(nrgb_full-1)/(nrgb-1);
      rgb[n][0] = rgb_full[nn][0];
      rgb[n][1] = rgb_full[nn][1];
      rgb[n][2] = rgb_full[nn][2];
      rgb[n][3] = transparentlevel;
    }
  }
  for(n=nrgb;n<nrgb+nrgb2;n++){
    rgb[n][0]=rgb2ptr[3*(n-nrgb)];
    rgb[n][1]=rgb2ptr[3*(n-nrgb)+1];
    rgb[n][2]=rgb2ptr[3*(n-nrgb)+2];
    rgb[n][3]=transparentlevel;
  }
  rgb_white=nrgb;
  rgb_yellow=nrgb+1;
  rgb_blue=nrgb+2;
  rgb_red=nrgb+3;

  if(flip==0){
    for(i=0;i<3;i++){
      foregroundcolor[i]=foregroundbasecolor[i];
      backgroundcolor[i]=backgroundbasecolor[i];
    }
    rgb[rgb_white][0]=1.0;
    rgb[rgb_white][1]=1.0;
    rgb[rgb_white][2]=1.0;
    rgb[rgb_black][0]=0.0;
    rgb[rgb_black][1]=0.0;
    rgb[rgb_black][2]=0.0;
  }
  else{
    for(i=0;i<3;i++){
      foregroundcolor[i]=backgroundbasecolor[i];
      backgroundcolor[i]=foregroundbasecolor[i];
    }
    rgb[rgb_white][0]=0.0;  //xxx fix or take ouit
    rgb[rgb_white][1]=0.0;
    rgb[rgb_white][2]=0.0;
    rgb[rgb_black][0]=1.0;
    rgb[rgb_black][1]=1.0;
    rgb[rgb_black][2]=1.0;
  }
  if(selected_case!=NULL){
    for(i=0;i<selected_case->nmeshes;i++){
      meshi=selected_case->meshinfo + i;
      vent_offset = 6*meshi->nbptrs;
      outline_offset = vent_offset + meshi->nvents;
      for(j=outline_offset;j<outline_offset+6;j++){
        facej = meshi->faceinfo + j;
        facej->color=foregroundcolor;
      }
    }
    updatechopcolors();
  }
  initcadcolors();
  update_texturebar();
}        

/* ------------------ updatechopcolors ------------------------ */

void updatechopcolors(){
  int i;
  float smin, smax;
  int ichopmin=0,ichopmax=nrgb_full;
#define NCHOP 8  
  int ii;

  last_particle_type=current_particle_type;
  for(i=0;i<nrgb_full;i++){
    rgb_slice[4*i+0]=rgb_full[i][0];
    rgb_slice[4*i+1]=rgb_full[i][1];
    rgb_slice[4*i+2]=rgb_full[i][2];
    rgb_slice[4*i+3]=rgb_full[i][3];
    rgb_part[4*i+0]=rgb_full[i][0];
    rgb_part[4*i+1]=rgb_full[i][1];
    rgb_part[4*i+2]=rgb_full[i][2];
    rgb_part[4*i+3]=rgb_full[i][3];
    rgb_plot3d[4*i+0]=rgb_full[i][0];
    rgb_plot3d[4*i+1]=rgb_full[i][1];
    rgb_plot3d[4*i+2]=rgb_full[i][2];
    rgb_plot3d[4*i+3]=rgb_full[i][3];
  }
  if(slicebounds!=NULL&&islicetype!=-1){
    smin=slicebounds[islicetype].valmin;
    smax=slicebounds[islicetype].valmax;
    if(setslicechopmin==1){
      ichopmin=nrgb_full*(slicechopmin-smin)/(smax-smin);
      if(ichopmin<0)ichopmin=0;
      if(ichopmin>nrgb_full-1)ichopmin=nrgb_full-1;
      for(i=0;i<ichopmin;i++){
        rgb_slice[4*i+3]=0.0;
      }
      for(i=ichopmin-NCHOP;i<ichopmin;i++){
        if(i<=0)continue;
        if(i>nrgb_full-1)continue;
        ii = i - (ichopmin-NCHOP);
        if(ii>NCHOP-1)continue;
        rgb_slice[4*i+3]=transparentlevel*(float)ii/(float)(NCHOP-1);
      }
    }
    if(setslicechopmax==1){
      ichopmax=nrgb_full*(slicechopmax - smin)/(smax-smin);
      if(ichopmax<0)ichopmax=0;
      if(ichopmax>nrgb_full-1)ichopmax=nrgb_full-1;
      for(i=ichopmax;i<nrgb_full;i++){
        rgb_slice[4*i+3]=0.0;
      }
      for(i=ichopmax;i<ichopmax+NCHOP;i++){
        if(i<=0)continue;
        if(i>nrgb_full-1)continue;
        ii = NCHOP-1-(i - ichopmax);
        if(ii>NCHOP-1)continue;
        rgb_slice[4*i+3]=transparentlevel*(float)ii/(float)(NCHOP-1);
      }
    }
  }

  if(partmax>partmin){
    if(setpartchopmin==1){
      ichopmin=nrgb_full*(partchopmin-partmin)/(partmax-partmin);
      if(ichopmin<0)ichopmin=0;
      if(ichopmin>nrgb_full-1)ichopmin=nrgb_full-1;
      for(i=0;i<ichopmin;i++){
        rgb_part[4*i+3]=0.0;
      }
      for(i=ichopmin-NCHOP;i<ichopmin;i++){
        if(i<=0)continue;
        if(i>nrgb_full-1)continue;
        ii = i - (ichopmin-NCHOP);
        if(ii>NCHOP-1)continue;
        rgb_part[4*i+3]=transparentlevel*(float)ii/(float)(NCHOP-1);
      }
    }
    if(setpartchopmax==1){
      ichopmax=nrgb_full*(partchopmax - partmin)/(partmax-partmin);
      if(ichopmax<0)ichopmax=0;
      if(ichopmax>nrgb_full-1)ichopmax=nrgb_full-1;
      for(i=ichopmax;i<nrgb_full;i++){
        rgb_part[4*i+3]=0.0;
      }
      for(i=ichopmax;i<ichopmax+NCHOP;i++){
        if(i<=0)continue;
        if(i>nrgb_full-1)continue;
        ii = NCHOP-1-(i - ichopmax);
        if(ii>NCHOP-1)continue;
        rgb_part[4*i+3]=transparentlevel*(float)ii/(float)(NCHOP-1);
      }
    } 
  }
  if(p3max_temp>p3min_temp){
    if(setp3chopmin_temp==1){
      ichopmin=nrgb_full*(p3chopmin_temp-p3min_temp)/(p3max_temp-p3min_temp);
      if(ichopmin<0)ichopmin=0;
      if(ichopmin>nrgb_full-1)ichopmin=nrgb_full-1;
      for(i=0;i<ichopmin;i++){
        rgb_plot3d[4*i+3]=0.0;
      }
      for(i=ichopmin-NCHOP;i<ichopmin;i++){
        if(i<=0)continue;
        if(i>nrgb_full-1)continue;
        ii = i - (ichopmin-NCHOP);
        if(ii>NCHOP-1)continue;
        rgb_plot3d[4*i+3]=transparentlevel*(float)ii/(float)(NCHOP-1);
      }
    }
    if(setp3chopmax_temp==1){
      ichopmax=nrgb_full*(p3chopmax_temp - p3min_temp)/(p3max_temp-p3min_temp);
      if(ichopmax<0)ichopmax=0;
      if(ichopmax>nrgb_full-1)ichopmax=nrgb_full-1;
      for(i=ichopmax;i<nrgb_full;i++){
        rgb_plot3d[4*i+3]=0.0;
      }
      for(i=ichopmax;i<ichopmax+NCHOP;i++){
        if(i<=0)continue;
        if(i>nrgb_full-1)continue;
        ii = NCHOP-1-(i - ichopmax);
        if(ii>NCHOP-1)continue;
        rgb_plot3d[4*i+3]=transparentlevel*(float)ii/(float)(NCHOP-1);
      }
    } 
  }

  update_texturebar();
}

/* ------------------ getrgb ------------------------ */

void getrgb(unsigned int val, unsigned char *rr, unsigned char *gg, unsigned char *bb){
  unsigned char r, g, b;

  r = val >> (ngreenbits+nbluebits);
  r = r << nredshift;

  g = val >> nbluebits;
  g = g&rgbmask[ngreenbits-1];
  g = g << ngreenshift;

  b = val&rgbmask[nbluebits-1];
  b = b << nblueshift;
  *rr=r; *gg=g; *bb=b;

  return;
}


/* ------------------ getcolorptr ------------------------ */

float *getcolorptr(const float *color){
  colordata *colorptr,*oldlastcolor,*lastcolor;

  int i;

  if(firstcolor==NULL){
    NewMemory((void *)&firstcolor,sizeof(colordata));
    for(i=0;i<4;i++){
      firstcolor->color[i]=color[i];
      firstcolor->full_color[i]=color[i];
    }
    firstcolor->bw_color[0] = 0.299*color[0]+0.587*color[1]+0.114*color[2];
    firstcolor->bw_color[1] = firstcolor->bw_color[0];
    firstcolor->bw_color[2] = firstcolor->bw_color[0];
    firstcolor->bw_color[3] = color[3];
    firstcolor->nextcolor=NULL;
    return firstcolor->color;
  }
  for(colorptr=firstcolor;colorptr!=NULL;colorptr=colorptr->nextcolor){
    oldlastcolor=colorptr;
    if(flabs(colorptr->color[0]-color[0])>0.0001)continue;
    if(flabs(colorptr->color[1]-color[1])>0.0001)continue;
    if(flabs(colorptr->color[2]-color[2])>0.0001)continue;
    if(flabs(colorptr->color[3]-color[3])>0.0001)continue;
    return colorptr->color;
  }
  lastcolor=NULL;
  NewMemory((void *)&lastcolor,sizeof(colordata));
  oldlastcolor->nextcolor=lastcolor;
  for(i=0;i<4;i++){
    lastcolor->color[i]=color[i];
    lastcolor->full_color[i]=color[i];
  }
  lastcolor->bw_color[0] = 0.299*color[0]+0.587*color[1]+0.114*color[2];
  lastcolor->bw_color[1] = lastcolor->bw_color[0];
  lastcolor->bw_color[2] = lastcolor->bw_color[0];
  lastcolor->bw_color[3] = color[3];
  lastcolor->nextcolor=NULL;
  return lastcolor->color;
}

/* ------------------ freecolors ------------------------ */
/*
void freecolors(void){
  colordata *colorptr,*nextcolor;
  extern colordata *firstcolor;

  firstcolor=NULL;
  for(colorptr=firstcolor;;){
    if(colorptr==NULL)return;
    nextcolor=colorptr->nextcolor;
    FREEMEMORY(colorptr);
    colorptr=nextcolor;
  }
}
*/

/* ------------------ colorconvert ------------------------ */

void colorconvert(int flag){
  colordata *colorptr;
  extern colordata *firstcolor;

  switch (flag){
   case TO_BW:
    for(colorptr=firstcolor;colorptr!=NULL;colorptr=colorptr->nextcolor){
      colorptr->color[0]=colorptr->bw_color[0];
      colorptr->color[1]=colorptr->bw_color[1];
      colorptr->color[2]=colorptr->bw_color[2];
    }
    break;
   case TO_COLOR:
    for(colorptr=firstcolor;colorptr!=NULL;colorptr=colorptr->nextcolor){
      colorptr->color[0]=colorptr->full_color[0];
      colorptr->color[1]=colorptr->full_color[1];
      colorptr->color[2]=colorptr->full_color[2];
    }
    break;
   default:
     ASSERT(FFALSE);
     break;
  }
}
