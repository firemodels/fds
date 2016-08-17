#include "options.h"
#include "glew.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "smokeviewvars.h"

#define EXPMIN -1
#define EXPMAX 3

/* ------------------ getBoundaryColors ------------------------ */

void getBoundaryColors(float *t, int nt, unsigned char *it,
              int settmin, float *ttmin, int settmax, float *ttmax,
              float *tmin_arg, float *tmax_arg,
              int ndatalevel, int nlevel,
              char **labels, char *scale, float *tvals256,
              int *extreme_min, int *extreme_max
              ){
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
  *tmin_arg = tmin2;
  *tmax_arg = tmax2;
  *extreme_min=0;
  *extreme_max=0;
  local_skip=0;
  AdjustDataBounds(t,local_skip,nt,settmin,&tmin2,settmax,&tmax2);
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
  if(range!=0.0f)factor = (ndatalevel-2*extreme_data_offset)/range;
  for(n=0;n<nt;n++){
    float val;

    val = *t;

    if(val<local_tmin){
      itt=0;
      *extreme_min=1;
    }
    else if(val>local_tmax){
      itt=ndatalevel-1;
      *extreme_max=1;
    }
    else{
      itt=extreme_data_offset+(int)(factor*(val-local_tmin));
    }
    *it++ = CLAMP(itt, colorbar_offset, ndatalevel - 1 - colorbar_offset);
    t++;
  }
  FrExp10(local_tmax, &expmax);
  FrExp10(local_tmin, &expmin);
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
    Num2String(&labels[n][0],tval);
  }
  tval = local_tmin + (nlevel-3)*factor;
  for(n=0;n<256;n++){
    tvals256[n] = (local_tmin*(255-n) + n*local_tmax)/255.;
  }
  Num2String(&labels[nlevel-2][0],tval);
  tval = local_tmax;
  Num2String(&labels[nlevel-1][0],tval);
}

/* ------------------ getBoundaryColors2 ------------------------ */

void getBoundaryColors2(float *t, int nt, unsigned char *it,
              int settmin, float *ttmin, int settmax, float *ttmax,
              float *tmin_arg, float *tmax_arg,
              int ndatalevel,
              int *extreme_min, int *extreme_max
              ){
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
  *tmin_arg = tmin2;
  *tmax_arg = tmax2;
  local_skip=0;
  AdjustDataBounds(t,local_skip,nt,settmin,&tmin2,settmax,&tmax2);
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
  if(range!=0.0f){
    factor = (ndatalevel-2*extreme_data_offset)/range;
  }
  for(n=0;n<nt;n++){
    float val;

    val=*t;
    if(val<local_tmin){
      *extreme_min=1;
      itt=0;
    }
    else if(val>local_tmax){
      *extreme_max=1;
      itt=ndatalevel-1;
    }
    else{
      itt=extreme_data_offset+(int)(0.5+(factor*(*t-local_tmin)));
    }
    *it++ = CLAMP(itt, colorbar_offset, ndatalevel - 1 - colorbar_offset);
    t++;
  }
}

/* ------------------ WriteBoundINI ------------------------ */

void WriteBoundINI(void){
  FILE *stream = NULL;
  char *fullfilename = NULL;
  int i;

  if(boundini_filename == NULL)return;
  fullfilename = get_filename(smokeviewtempdir, boundini_filename, tempdir_flag);

  if(fullfilename == NULL)return;

  for(i = 0; i < npatchinfo; i++){
    bounddata *boundi;
    patchdata *patchi;
    int skipi;
    int j;

    skipi = 0;
    patchi = patchinfo + i;
    if(patchi->bounds.defined == 0)continue;
    for(j = 0; j < i - 1; j++){
      patchdata *patchj;

      patchj = patchinfo + j;
      if(patchi->type == patchj->type&&patchi->filetype == patchj->filetype){
        skipi = 1;
        break;
      }
    }
    if(skipi == 1)continue;

    boundi = &patchi->bounds;
    if(stream == NULL){
      stream = fopen(fullfilename, "w");
      if(stream == NULL){
        FREEMEMORY(fullfilename);
        return;
      }
    }
    fprintf(stream, "B_BOUNDARY\n");
    fprintf(stream, " %f %f %f %f %i %s\n", boundi->global_min, boundi->percentile_min, boundi->percentile_max, boundi->global_max, patchi->filetype, patchi->label.shortlabel);
  }
  if(stream != NULL)fclose(stream);
  FREEMEMORY(fullfilename);
}

/* ------------------ UpdatePatchBounds ------------------------ */

void UpdatePatchBounds(patchdata *patchi){
  histogramdata full_histogram;
  bounddata *boundi;
  int j;

  boundi = &patchi->bounds;
  if(boundi->defined==1)return;
  init_histogram(&full_histogram,NHIST_BUCKETS);

  for(j=0;j<npatchinfo;j++){
    patchdata *patchj;

    patchj=patchinfo+j;
    if(patchj->type!=patchi->type||patchj->filetype!=patchi->filetype)continue;
    merge_histogram(&full_histogram,patchj->histogram);
  }

  boundi->global_min = get_histogram_value(&full_histogram, 0.0);
  boundi->global_max = get_histogram_value(&full_histogram, 1.0);
  boundi->percentile_min = get_histogram_value(&full_histogram, percentile_level);
  boundi->percentile_max = get_histogram_value(&full_histogram, 1.0-percentile_level);
  boundi->defined=1;

  for(j=0;j<npatchinfo;j++){
    bounddata *boundj;
    patchdata *patchj;

    patchj=patchinfo+j;
    if(patchi==patchj||patchj->type!=patchi->type||patchj->filetype!=patchi->filetype)continue;

    boundj = &patchj->bounds;
    memcpy(boundj,boundi,sizeof(bounddata));
  }
  WriteBoundINI();
  free_histogram(&full_histogram);
}

/* ------------------ getBoundaryColors3 ------------------------ */

void getBoundaryColors3(patchdata *patchi, float *t, int nt, unsigned char *it,
              int settmin, float *ttmin, int settmax, float *ttmax,
              float *tmin_arg, float *tmax_arg,
              int nlevel,
              char **labels, char *scale, float *tvals256,
              int *extreme_min, int *extreme_max
              ){
  int n;
  float factor, tval, range;
  int expmin, expmax;
  int itt;
  float new_tmin, new_tmax, tmin2, tmax2;

  UpdatePatchBounds(patchi);

  CheckMemory;
  tmin2=patchi->bounds.global_min;
  tmax2=patchi->bounds.global_max;

  *tmin_arg = tmin2;
  *tmax_arg = tmax2;
  *extreme_min=0;
  *extreme_max=0;
  if(settmin==PERCENTILE_MIN){
    tmin2=patchi->bounds.percentile_min;
  }
  if(settmax==PERCENTILE_MAX){
    tmax2=patchi->bounds.percentile_max;
  }
  if(axislabels_smooth==1){
    SmoothLabel(&tmin2,&tmax2,nrgb);
  }
  if(settmin!=SET_MIN){
    *ttmin=tmin2;
  }
  if(settmax!=SET_MAX){
    *ttmax=tmax2;
  }
  new_tmin = *ttmin;
  new_tmax = *ttmax;

  patchi->local_valmin=new_tmin;
  patchi->local_valmax=new_tmax;

  CheckMemory;
  range = new_tmax - new_tmin;
  factor = 0.0f;
  if(range!=0.0f)factor = (float)(255-2*extreme_data_offset)/range;
  for(n=0;n<nt;n++){
    float val;

    val = *t;

    if(val<new_tmin){
      itt=0;
      *extreme_min=1;
    }
    else if(val>new_tmax){
      itt=255;
      *extreme_max=1;
    }
    else{
      itt=extreme_data_offset+(int)(factor*(val-new_tmin));
    }
    *it++=CLAMP(itt,colorbar_offset,255-colorbar_offset);
    t++;
  }
  CheckMemory;
  STRCPY(scale,"");
  FrExp10(new_tmax, &expmax);
  FrExp10(new_tmin, &expmin);
  if(expmin!=0&&expmax!=0&&expmax-expmin<=2&&(expmin<-2||expmin>2)){
    new_tmin *= pow((double)10.0,(double)-expmin);
    new_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(scale,"*10^%i",expmin);
  }
  if(expmin==0&&(expmax<EXPMIN||expmax>EXPMAX)){
    new_tmin *= pow((double)10.0,(double)-expmax);
    new_tmax *= pow((double)10.0,(double)-expmax);
    sprintf(scale,"*10^%i",expmax);
  }
  if(expmax==0&&(expmin<EXPMIN||expmin>EXPMAX)){
    new_tmin *= pow((double)10.0,(double)-expmin);
    new_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(scale,"*10^%i",expmin);
  }
  range = new_tmax - new_tmin;
  factor = range/(nlevel-2);
  for (n=1;n<nlevel-2;n++){
    tval = new_tmin + (n-1)*factor;
    Num2String(&labels[n][0],tval);
  }
  tval = new_tmin + (nlevel-3)*factor;
  for(n=0;n<256;n++){
    tvals256[n] = (new_tmin*(255-n) + n*new_tmax)/255.;
  }
  Num2String(&labels[nlevel-2][0],tval);
  tval = new_tmax;
  Num2String(&labels[nlevel-1][0],tval);
}

/* ------------------ UpdateAllPatchColors ------------------------ */

void UpdateAllPatchColors(void){
  int i;

  for(i=0;i<nmeshes;i++){
    meshdata *meshi;
    patchdata *patchi;
    int npatchvals;
    float patchmin_global, patchmax_global;

    meshi = meshinfo + i;
    if(meshi->patchval==NULL||meshi->cpatchval==NULL||meshi->patchfilenum<0)continue;
    patchi = patchinfo + meshi->patchfilenum;
    if(patchi->loaded==0)continue;

    npatchvals = meshi->npatch_times*meshi->npatchsize;

    getBoundaryColors3(patchi,meshi->patchval, npatchvals, meshi->cpatchval,
    setpatchmin,&patchmin, setpatchmax,&patchmax,
    &patchmin_global, &patchmax_global,
    nrgb, colorlabelpatch,patchi->scale,boundarylevels256,
    &patchi->extreme_min,&patchi->extreme_max);
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

  FrExp10(local_tmax, &expmax);
  FrExp10(local_tmin, &expmin);
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
    Num2String(&labels[n][0],tval);
  }
  tval = local_tmin + (nlevel-3)*factor;
  for(n=0;n<256;n++){
    tvals256[n] = (local_tmin*(255-n) + n*local_tmax)/255.;
  }
  Num2String(&labels[nlevel-2][0],tval);
  tval = local_tmax;
  Num2String(&labels[nlevel-1][0],tval);
}

/* ------------------ UpdatePart5Extremes ------------------------ */

void UpdatePart5Extremes(void){
  int ii,i,j,k,m;
  part5data *datacopy;

  for(i=0;i<npart5prop;i++){
    partpropdata *propi;

    propi = part5propinfo + i;
    propi->extreme_max=0;
    propi->extreme_min=0;
  }


  for(ii=0;ii<npartinfo;ii++){
    partdata *parti;

    parti = partinfo + ii;
    if(parti->loaded==0||parti->display==0)continue;
    datacopy = parti->data5;
    for(i=0;i<parti->ntimes;i++){
      for(j=0;j<parti->nclasses;j++){
        partclassdata *partclassi;
        unsigned char *irvals;

        partclassi = parti->partclassptr[j];
        irvals = datacopy->irvals;
        for(k=2;k<partclassi->ntypes;k++){
          partpropdata *prop_id;

          prop_id = get_partprop(partclassi->labels[k].longlabel);
          if(prop_id==NULL)continue;

          if(strcmp(partclassi->labels[k].longlabel,"HUMAN_COLOR")==0){
          }
          else{
            for(m=0;m<datacopy->npoints;m++){
              int irval;

              irval=*irvals++;
              if(irval==0)prop_id->extreme_min=1;
              if(irval==255)prop_id->extreme_max=1;
            }
          }
        }
        datacopy++;
      }
    }
  }
}

/* ------------------ GetPart5Colors ------------------------ */

void GetPart5Colors(partdata *parti, int nlevel, int convert_flag){
  int i;
  part5data *datacopy;
  // float *diameter_data;
  float *length_data, *azimuth_data, *elevation_data;
  float *u_vel_data, *v_vel_data, *w_vel_data;

  datacopy = parti->data5;
  for(i=0;i<parti->ntimes;i++){
    int j;

    for(j=0;j<parti->nclasses;j++){
      float valmin, valmax, dval;
      partclassdata *partclassi;
      float *rvals;
      unsigned char *irvals;
      float *dsx, *dsy, *dsz;
      int flag, k;

      partclassi = parti->partclassptr[j];
      rvals = datacopy->rvals;
      irvals = datacopy->irvals;
      for(k=2;k<partclassi->ntypes;k++){
        partpropdata *prop_id;

        prop_id = get_partprop(partclassi->labels[k].longlabel);
        if(prop_id==NULL)continue;

        if(strcmp(partclassi->labels[k].longlabel,"HUMAN_COLOR")==0){
          if(convert_flag==PARTFILE_MAP){
            int m;

            for(m = 0; m<datacopy->npoints; m++){
              float val;

              val = *rvals++;
              *irvals++ = CLAMP(val+0.5, 0, navatar_colors-1);
            }
          }
        }
        else{
          int prop_id_index;
          float partimin, partimax;

          valmin = prop_id->valmin;
          valmax = prop_id->valmax;
          dval = valmax - valmin;
          if(dval<=0.0)dval=1.0;
          prop_id_index = prop_id-part5propinfo;
          partimin = parti->valmin[prop_id_index];
          partimax = parti->valmax[prop_id_index];

          if(convert_flag==PARTFILE_MAP){
            int m;

            for(m = 0; m<datacopy->npoints; m++){
              float val;
              int irval;

              val = *rvals++;
              if(val<valmin){
                irval = 0;
              }
              else if(val>valmax){
                irval = 255;
              }
              else{
                irval = extreme_data_offset+(float)(255-2*extreme_data_offset)*(val-valmin)/dval;
              }
              *irvals++ = CLAMP(irval, 0, 255);
            }
          }
          else if(convert_flag==PARTFILE_REMAP){
            int m;

            for(m = 0; m<datacopy->npoints; m++){
              float val;
              int irval;

              irval = *irvals;
              val = partimin+(float)(irval-extreme_data_offset)*(partimax-partimin)/(255.0-2.0*extreme_data_offset);
              if(val<valmin){
                irval = 0;
              }
              else if(val>valmax){
                irval = 255;
              }
              else{
                irval = extreme_data_offset+(float)(255-2*extreme_data_offset)*(val-valmin)/dval;
              }
              *irvals++ = CLAMP(irval, 0, 255);
            }
          }
          else{
            ASSERT(FFALSE);
          }
        }
      }
      //** do some data conversion if the right data columns are present
      azimuth_data=NULL;
//      diameter_data=NULL;
      elevation_data=NULL;
      length_data=NULL;
      u_vel_data=NULL;
      v_vel_data=NULL;
      w_vel_data=NULL;

      if(partclassi->col_azimuth>=0){
        azimuth_data=datacopy->rvals+partclassi->col_azimuth*datacopy->npoints;
      }
      if(partclassi->col_diameter>=0){
       // diameter_data=datacopy->rvals+partclassi->col_diameter*datacopy->npoints;
      }
      if(partclassi->col_elevation>=0){
        elevation_data=datacopy->rvals+partclassi->col_elevation*datacopy->npoints;
      }
      if(partclassi->col_length>=0){
        length_data=datacopy->rvals+partclassi->col_length*datacopy->npoints;
      }
      if(partclassi->col_u_vel>=0){
        u_vel_data=datacopy->rvals+partclassi->col_u_vel*datacopy->npoints;
      }
      if(partclassi->col_v_vel>=0){
        v_vel_data=datacopy->rvals+partclassi->col_v_vel*datacopy->npoints;
      }
      if(partclassi->col_w_vel>=0){
        w_vel_data=datacopy->rvals+partclassi->col_w_vel*datacopy->npoints;
      }
      flag=0;
      if(azimuth_data!=NULL&&elevation_data!=NULL&&length_data!=NULL){
        int m;

        flag=1;
        dsx = datacopy->dsx;
        dsy = datacopy->dsy;
        dsz = datacopy->dsz;
        for(m=0;m<datacopy->npoints;m++){
          float az, elev, length;

          az= azimuth_data[m]*DEG2RAD;
          elev = elevation_data[m]*DEG2RAD;
          length=SCALE2SMV(length_data[m]);
          dsx[m] = cos(az)*cos(elev)*length/2.0;
          dsy[m] = sin(az)*cos(elev)*length/2.0;
          dsz[m] =         sin(elev)*length/2.0;
        }
      }
      if(u_vel_data!=NULL&&v_vel_data!=NULL&&w_vel_data!=NULL){
        float denom;
        int m;
        partpropdata *prop_U, *prop_V, *prop_W;

        prop_U = get_partprop(partclassi->labels[partclassi->col_u_vel+2].longlabel);
        prop_V = get_partprop(partclassi->labels[partclassi->col_v_vel+2].longlabel);
        prop_W = get_partprop(partclassi->labels[partclassi->col_w_vel+2].longlabel);
        if(prop_U!=NULL&&prop_V!=NULL&&prop_W!=NULL){
          float umax, vmax, wmax;

          umax = MAX(ABS(prop_U->valmin),ABS(prop_U->valmax));
          vmax = MAX(ABS(prop_V->valmin),ABS(prop_V->valmax));
          wmax = MAX(ABS(prop_W->valmin),ABS(prop_W->valmax));

          denom = sqrt(umax*umax+vmax*vmax+wmax*wmax);
          if(denom==0.0)denom=1.0;
        }
        else{
          denom=1.0;
        }

        flag=1;
        dsx = datacopy->dsx;
        dsy = datacopy->dsy;
        dsz = datacopy->dsz;
        for(m=0;m<datacopy->npoints;m++){
          dsx[m] = 0.05*u_vel_data[m]/denom;
          dsy[m] = 0.05*v_vel_data[m]/denom;
          dsz[m] = 0.05*w_vel_data[m]/denom;
        }
      }
      if(flag==0){
        FREEMEMORY(datacopy->dsx);
        FREEMEMORY(datacopy->dsy);
        FREEMEMORY(datacopy->dsz);
      }
      datacopy++;
    }
  }
// erase data memory in a separate loop (so all "columns" are available when doing any conversions)
  datacopy = parti->data5;
  if(parti->data_type == PARTDATA){
    for(i = 0; i < parti->ntimes; i++){
      int j;

      for(j = 0; j < parti->nclasses; j++){
        FREEMEMORY(datacopy->rvals);
        datacopy++;
      }
    }
  }
  for(i=0;i<npart5prop;i++){
    int n;
    partpropdata *propi;
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

    strcpy(scale,"");

    FrExp10(local_tmax, &expmax);
    FrExp10(local_tmin, &expmin);
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
      Num2String(&labels[n][0],tval);
    }
    for(n=0;n<256;n++){
      ppartlevels256[n] = (local_tmin*(255-n) + n*local_tmax)/255.;
    }
    tval = local_tmin + (nlevel-3)*factor;
    Num2String(&labels[nlevel-2][0],tval);
    tval = local_tmax;
    Num2String(&labels[nlevel-1][0],tval);
    CheckMemory;
  }

}

/* ------------------ GetZoneColor ------------------------ */

int GetZoneColor(float t, float local_tmin, float local_tmax, int nlevel){
  int level;

  if(t<=local_tmin)return 0;
  if(t>=local_tmax)return nlevel-1;
  if(local_tmin==local_tmax)return 0;
  level=nlevel*(t-local_tmin)/(local_tmax-local_tmin);
  if(level<0)return 0;
  if(level>nlevel-1)return nlevel-1;
  return level;
}


/* ------------------ GetZoneColors ------------------------ */

void GetZoneColors(const float *t, int nt, unsigned char *it,
               float ttmin, float ttmax, int nlevel, int nlevel_full,
               char **labels, char *scale, float *tvals256
               ){
  int n;
  float dt, factor;
  int itt;
  int expmin, expmax;
  float local_tmin, local_tmax;
  float range;
  float tval;

  local_tmin = ttmin;
  local_tmax = ttmax;

  dt = local_tmax - local_tmin;
  factor=0.0f;
  if(dt!=0.0f)factor = (nlevel_full-2*extreme_data_offset)/dt;
  for(n=0;n<nt;n++){
    float val;

    val=*t;
    if(val<local_tmin){
      itt=0;
    }
    else if(val>local_tmax){
      itt=nlevel_full-1;
    }
    else{
      itt=extreme_data_offset+(int)(factor*(val-local_tmin));
    }
    *it++=CLAMP(itt,colorbar_offset,nlevel_full-1-colorbar_offset);
    t++;
  }

  STRCPY(scale,"");

  FrExp10(local_tmax, &expmax);
  FrExp10(local_tmin, &expmin);
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
    Num2String(&labels[n][0],tval);
  }
  tval = local_tmin + (nlevel-3)*factor;
  for(n=0;n<256;n++){
    tvals256[n] = (local_tmin*(255-n) + n*local_tmax)/255.;
  }
  Num2String(&labels[nlevel-2][0],tval);
  tval = local_tmax;
  Num2String(&labels[nlevel-1][0],tval);

}

/* ------------------ GetPlot3DColors ------------------------ */

void GetPlot3DColors(int plot3dvar, int settmin, float *ttmin, int settmax, float *ttmax,
              int ndatalevel, int nlevel,
              char **labels,char **labelsiso,char **scale, float *fscale, float *tlevels, float *tlevels256,
              int *extreme_min, int *extreme_max
              ){
  int n;
  float dt, factor, tval;
  float local_tmin, local_tmax, tmin2, tmax2;
  float range;
  int expmax,expmin;
  float tminorig, tmaxorig, dtorig;
  int itt;
  float *q;
  unsigned char *iq;
  plot3ddata *p;
  meshdata *meshi;
  char *iblank;
  int i;
  int ntotal;

  tmin2=*ttmin;
  tmax2=*ttmax;

  tmin2= 1000000000.;
  tmax2=-1000000000.;
  *extreme_min=0;
  *extreme_max=0;
  for(i=0;i<nplot3dinfo;i++){
    p = plot3dinfo+i;
    if(p->loaded==0||p->display==0)continue;
    meshi = meshinfo+p->blocknumber;
    ntotal=(meshi->ibar+1)*(meshi->jbar+1)*(meshi->kbar+1);
    iblank=meshi->c_iblank_node;
    if(cache_qdata==1||meshi->qdata!=NULL){
      q=meshi->qdata+plot3dvar*ntotal;
      for(n=0;n<ntotal;n++){
        if(iblank==NULL||*iblank++==GAS){
          if(*q<tmin2)tmin2=*q;
          if(*q>tmax2)tmax2=*q;
        }
        q++;
      }
    }
    else{
      float qval, *qvals;

      qvals=p3levels256[plot3dvar];
      iq=meshi->iqdata+plot3dvar*ntotal;
      for(n=0;n<ntotal;n++){
        qval=qvals[*iq++];
        if(*iblank++==GAS){
          if(qval<tmin2)tmin2=qval;
          if(qval>tmax2)tmax2=qval;
        }
      }
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
  AdjustPlot3DBounds(plot3dvar,settmin,&local_tmin,settmax,&local_tmax);

  *ttmin=local_tmin;
  *ttmax=local_tmax;
  range = local_tmax-local_tmin;
  tminorig=local_tmin;
  tmaxorig=local_tmax;
  if(range!=0.0f){
    factor = (float)(ndatalevel-2*extreme_data_offset)/range;
  }
  else{
    factor = 0.0f;
  }

  for(i=0;i<nplot3dinfo;i++){
    p = plot3dinfo+i;
    if(p->loaded==0||p->display==0)continue;
    meshi = meshinfo+p->blocknumber;
    ntotal=(meshi->ibar+1)*(meshi->jbar+1)*(meshi->kbar+1);

    if(cache_qdata==1||meshi->qdata!=NULL){
      q=meshi->qdata+plot3dvar*ntotal;
      iq=meshi->iqdata+plot3dvar*ntotal;
      for(n=0;n<ntotal;n++){
        float val;

        val=*q;
        if(val<local_tmin){
          itt=0;
          *extreme_min=1;
        }
        else if(val>local_tmax){
          itt=ndatalevel-1;
          *extreme_max=1;
        }
        else{
          itt=extreme_data_offset+(int)(factor*(val-local_tmin));
        }
        *iq++=CLAMP(itt,colorbar_offset,ndatalevel-1-colorbar_offset);
        q++;
      }
    }
  }

  STRCPY(*scale,"");
  FrExp10(local_tmax, &expmax);
  FrExp10(local_tmin, &expmin);
  if(expmin!=0&&expmax!=0&&expmax-expmin<=2&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(*scale,"*10^%i",expmin);
    *fscale=pow(10.0,(float)expmin);
  }
  if(expmin==0&&(expmax<EXPMIN||expmax>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmax);
    local_tmax *= pow((double)10.0,(double)-expmax);
    sprintf(*scale,"*10^%i",expmax);
    *fscale=pow(10.0,(float)expmax);
  }
  if(expmax==0&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(*scale,"*10^%i",expmin);
    *fscale=pow(10.0,(float)expmin);
  }

  range = local_tmax-local_tmin;
  dt = range/(float)(nlevel-1);
  dtorig = (tmaxorig-tminorig)/(float)(nlevel-1);
  for (n=0;n<nlevel-1;n++){
    tval = local_tmin + n*dt;
    Num2String(&labels[n][0],tval);
    Num2String(&labelsiso[n][0],tval+dt/2.0);
  }
  for(n=0;n<256;n++){
    tlevels256[n] = (local_tmin*(255-n) + local_tmax*n)/255.;
  }
  tval = local_tmax;
  Num2String(&labels[nlevel-1][0],tval);
  Num2String(&labelsiso[nlevel-1][0],tval);

  for(n=0;n<nlevel;n++){
    tlevels[n]=tminorig+(float)n*dtorig;
  }

  for(i=0;i<nplot3dinfo;i++){
    p = plot3dinfo+i;
    if(p->loaded==0||p->display==0)continue;
    meshi = meshinfo+p->blocknumber;
    ntotal=(meshi->ibar+1)*(meshi->jbar+1)*(meshi->kbar+1);

    if(cache_qdata==0&&meshi->qdata==NULL){
      float qval, *qvals;

      qvals=p3levels256[plot3dvar];
      iq=meshi->iqdata+plot3dvar*ntotal;
      for(n=0;n<ntotal;n++){
        qval=qvals[*iq];
        itt=(int)(factor*(qval-local_tmin));
        if(itt<0)itt=0;
        if(itt>ndatalevel-1)itt=ndatalevel-1;
        *iq++=itt;
      }
    }
  }
}

/* ------------------ GetSliceColors ------------------------ */

void GetSliceColors(const float *t, int nt, unsigned char *it,
              float local_tmin, float local_tmax,
              int ndatalevel, int nlevel,
              char labels[12][11],char **scale, float *fscale, float *tlevels256,
              int *extreme_min, int *extreme_max
              ){
  int n;
  float dt, factor, tval;
  float range;
  int expmax,expmin;
  int itt;

  range = local_tmax-local_tmin;
  *extreme_min=0;
  *extreme_max=0;
  if(range!=0.0f){
    factor = (float)(ndatalevel-2*extreme_data_offset)/range;
  }
   else{
     factor = 0.0f;
   }
  for(n=0;n<nt;n++){
    float val;

    val = *t;

    if(val<local_tmin){
      itt=0;
      *extreme_min=1;
    }
    else if(val>local_tmax){
      itt=ndatalevel-1;
      *extreme_max=1;
    }
    else{
      itt=extreme_data_offset+(int)(factor*(val-local_tmin));
    }
    *it++ = CLAMP(itt, colorbar_offset, ndatalevel - 1 - colorbar_offset);
    t++;
  }

  STRCPY(*scale,"");
  FrExp10(local_tmax, &expmax);
  FrExp10(local_tmin, &expmin);
  *fscale=1.0;
  if(expmin!=0&&expmax!=0&&expmax-expmin<=2&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(*scale,"*10^%i",expmin);
    *fscale=pow(10.0,(float)expmin);
  }
  if(expmin==0&&(expmax<EXPMIN||expmax>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmax);
    local_tmax *= pow((double)10.0,(double)-expmax);
    sprintf(*scale,"*10^%i",expmax);
    *fscale=pow(10.0,(float)expmax);
  }
  if(expmax==0&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(*scale,"*10^%i",expmin);
    *fscale=pow(10.0,(float)expmin);
  }

  range = local_tmax-local_tmin;
  dt = range/(float)(nlevel-2);
  for (n=1;n<nlevel-1;n++){
    tval = local_tmin + (n-1)*dt;
    Num2String(&labels[n][0],tval);
  }
  for(n=0;n<256;n++){
    tlevels256[n] = (local_tmin*(255-n) + local_tmax*n)/255.;
  }
  tval = local_tmax;
  Num2String(&labels[nlevel-1][0],tval);
}

/* ------------------ getSliceLabelels ------------------------ */

void GetSliceLabels(float local_tmin, float local_tmax, int nlevel,
              char labels[12][11],char **scale, float *fscale, float *tlevels256){
  int n;
  float dt, tval;
  float range;
  int expmax,expmin;

  range = local_tmax-local_tmin;

  STRCPY(*scale,"");
  *fscale=1.0;
  FrExp10(local_tmax, &expmax);
  FrExp10(local_tmin, &expmin);
  if(expmin!=0&&expmax!=0&&expmax-expmin<=2&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(*scale,"*10^%i",expmin);
    *fscale=pow(10.0,(float)expmin);
  }
  if(expmin==0&&(expmax<EXPMIN||expmax>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmax);
    local_tmax *= pow((double)10.0,(double)-expmax);
    sprintf(*scale,"*10^%i",expmax);
    *fscale=pow(10.0,(float)expmax);
  }
  if(expmax==0&&(expmin<EXPMIN||expmin>EXPMAX)){
    local_tmin *= pow((double)10.0,(double)-expmin);
    local_tmax *= pow((double)10.0,(double)-expmin);
    sprintf(*scale,"*10^%i",expmin);
    *fscale=pow(10.0,(float)expmin);
  }

  range = local_tmax-local_tmin;
  dt = range/(float)(nlevel-2);
  for (n=1;n<nlevel-1;n++){
    tval = local_tmin + (n-1)*dt;
    Num2String(&labels[n][0],tval);
  }
  for(n=0;n<256;n++){
    tlevels256[n] = (local_tmin*(255-n) + local_tmax*n)/255.;
  }
  tval = local_tmax;
  Num2String(&labels[nlevel-1][0],tval);
}


/* ------------------ GetIsoLabels ------------------------ */

void GetIsoLabels(float local_tmin, float local_tmax, int nlevel,
              char labels[12][11],char **scale, float *tlevels256){
  int n;
  float dt, tval;
  float range;
  int expmax,expmin;

  range = local_tmax-local_tmin;

  STRCPY(*scale,"");
  FrExp10(local_tmax, &expmax);
  FrExp10(local_tmin, &expmin);
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
    Num2String(&labels[n][0],tval);
  }
  for(n=0;n<256;n++){
    tlevels256[n] = (local_tmin*(255-n) + local_tmax*n)/255.;
  }
  tval = local_tmax;
  Num2String(&labels[nlevel-1][0],tval);
}

/* ------------------ InitCadColors ------------------------ */

void InitCadColors(void){
  int n, i1, i2, i;
  float xx, f1, f2, sum;
  switch(setbw){
   case 0:
    for(n=0;n<nrgb_cad;n++){
      xx = (float)n/(float)nrgb_cad * (float)(nrgb-1);
      i1 = (int)xx;
      i2 = (int)(xx+1);
      f2 = xx - (float)i1;
      f1 = 1.0f - f2;
      sum=0.0;
      for(i=0;i<3;i++){
        rgb_cad[n][i] = f1*rgb[i1][i] + f2*rgb[i2][i];
        sum += rgb_cad[n][i]*rgb_cad[n][i];
      }
      sum=sqrt((double)sum);
      if(sum>0.0){
        for(i=0;i<3;i++){
          rgb_cad[n][i] /= sum;
        }
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

/* ------------------ UpdateTexturebar ------------------------ */

void UpdateTexturebar(void){
  if(use_graphics==0)return;
  glBindTexture(GL_TEXTURE_1D, terrain_colorbar_id);
  glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, 256, 0, GL_RGBA, GL_FLOAT, rgb_terrain2);
  SNIFF_ERRORS("UpdateTexturebar - glTexImage1D (rgb_terrain2) ");

  glBindTexture(GL_TEXTURE_1D, texture_colorbar_id);
  glTexImage1D(GL_TEXTURE_1D,0,GL_RGBA,256,0,GL_RGBA,GL_FLOAT,rgb_full);
  SNIFF_ERRORS("UpdateTexturebar - glTexImage1D (rgb_full) ");

  glBindTexture(GL_TEXTURE_1D,texture_slice_colorbar_id);
  glTexImage1D(GL_TEXTURE_1D,0,GL_RGBA,256,0,GL_RGBA,GL_FLOAT,rgb_slice);
  SNIFF_ERRORS("UpdateTexturebar - glTexImage1D (rgb_slice) ");

  glBindTexture(GL_TEXTURE_1D,texture_patch_colorbar_id);
  glTexImage1D(GL_TEXTURE_1D,0,GL_RGBA,256,0,GL_RGBA,GL_FLOAT,rgb_patch);
  SNIFF_ERRORS("UpdateTexturebar - glTexImage1D (rgb_patch) ");

  glBindTexture(GL_TEXTURE_1D,texture_plot3d_colorbar_id);
  glTexImage1D(GL_TEXTURE_1D,0,GL_RGBA,256,0,GL_RGBA,GL_FLOAT,rgb_plot3d);
  SNIFF_ERRORS("UpdateTexturebar - glTexImage1D (rgb_plot3d) ");

  glBindTexture(GL_TEXTURE_1D,texture_iso_colorbar_id);
  glTexImage1D(GL_TEXTURE_1D,0,GL_RGBA,256,0,GL_RGBA,GL_FLOAT,rgb_iso);
  SNIFF_ERRORS("UpdateTexturebar - glTexImage1D (rgb_iso) ");

  glBindTexture(GL_TEXTURE_1D,slicesmoke_colormap_id);
  glTexImage1D(GL_TEXTURE_1D,0,GL_RGBA,MAXSMOKERGB,0,GL_RGBA,GL_FLOAT,rgb_slicesmokecolormap);
  SNIFF_ERRORS("UpdateTexturebar - glTexImage1D (rgb_slicesmokecolormap) ");

  glBindTexture(GL_TEXTURE_1D,volsmoke_colormap_id);
  glTexImage1D(GL_TEXTURE_1D,0,GL_RGBA,MAXSMOKERGB,0,GL_RGBA,GL_FLOAT,rgb_volsmokecolormap);
  SNIFF_ERRORS("UpdateTexturebar - glTexImage1D (rgb_volsmokecolormap) ");

#ifdef pp_GPU
  if(gpuactive==1&&nvolrenderinfo>0&&showvolrender==1){
    glActiveTexture(GL_TEXTURE2);
    glTexSubImage1D(GL_TEXTURE_1D,0,0,MAXSMOKERGB,GL_RGBA,GL_FLOAT, rgb_volsmokecolormap);
    SNIFF_ERRORS("UpdateTexturebar - glTexSubImage1D (rgb_volsmokecolormap) ");
    glActiveTexture(GL_TEXTURE0);
  }
  if(gpuactive==1&&SHOW_gslice_data==1){
    glActiveTexture(GL_TEXTURE4);
    glTexSubImage1D(GL_TEXTURE_1D,0,0,256,GL_RGBA,GL_FLOAT, rgb_slice);
    SNIFF_ERRORS("updatecolors after glTexSubImage1D (rgb_slice)");
    glActiveTexture(GL_TEXTURE0);
  }
#endif

}

/* ------------------ InitRGB ------------------------ */

void InitRGB(void){
  int n;
  float transparent_level_local=1.0;

  if(use_transparency_data==1)transparent_level_local=transparent_level;

  if(setbw==0){
    colorconvert(TO_COLOR);
    if(nrgb_ini !=0){
      nrgb = nrgb_ini;
      for(n=0;n<nrgb_ini;n++){
        rgb[n][0] = rgb_ini[n*3];
        rgb[n][1] = rgb_ini[n*3+1];
        rgb[n][2] = rgb_ini[n*3+2];
        rgb[n][3] = transparent_level_local;
      }
    }
    else{
      for(n=0;n<nrgb;n++){
        rgb[n][0] = rgb_base[n][0];
        rgb[n][1] = rgb_base[n][1];
        rgb[n][2] = rgb_base[n][2];
        rgb[n][3] = transparent_level_local;
      }
    }
  }
  else{
    colorconvert(TO_BW);
    for(n=0;n<nrgb;n++){
      rgb[n][0] = bw_base[n][0];
      rgb[n][1] = bw_base[n][1];
      rgb[n][2] = bw_base[n][2];
      rgb[n][3] = transparent_level_local;
    }
  }
}

/* ------------------ UpdateSmokeColormap ------------------------ */

void UpdateSmokeColormap(int option){
  int n;
  float transparent_level_local=1.0;
  unsigned char *alpha;
  float *fire_cb;
  float val, valmin, valmax, valcut;
  int icut;
  float *rgb_colormap;
  int have_fire;

  have_fire = HaveFire();
  if(option==RENDER_SLICE){
    valmin=global_hrrpuv_min;
    valcut=global_hrrpuv_cutoff;
    valmax=global_hrrpuv_max;
    rgb_colormap = rgb_slicesmokecolormap;
  }
  else{
    valmin=temperature_min;
    valcut=temperature_cutoff;
    valmax=temperature_max;
    rgb_colormap = rgb_volsmokecolormap;
  }
  icut = (MAXSMOKERGB-1)*((valcut-valmin)/(valmax-valmin));
  icut = CLAMP(icut,2,(MAXSMOKERGB-3));

  if(use_transparency_data==1)transparent_level_local=transparent_level;

  alpha = colorbarinfo[colorbartype].alpha;
  fire_cb = colorbarinfo[fire_colorbar_index].colorbar;

  switch(firecolormap_type){
    case FIRECOLORMAP_DIRECT:
      for(n=0;n<MAXSMOKERGB;n++){
        if(n<icut||have_fire==0){
          rgb_colormap[4*n+0] = (float)smoke_red / 255.0;
          rgb_colormap[4*n+1] = (float)smoke_green / 255.0;
          rgb_colormap[4*n+2] = (float)smoke_blue / 255.0;
        }
        else{
          rgb_colormap[4*n+0]=(float)fire_red/255.0;
          rgb_colormap[4*n+1]=(float)fire_green/255.0;
          rgb_colormap[4*n+2]=(float)fire_blue/255.0;
        }
        if(alpha[n]==0){
          rgb_colormap[4*n+3]=0.0;
        }
        else{
          rgb_colormap[4*n+3]=transparent_level_local;
        }
      }
      UpdateTexturebar();
      break;
    case FIRECOLORMAP_NOCONSTRAINT:
    case FIRECOLORMAP_CONSTRAINT:
      for(n=0;n<MAXSMOKERGB;n++){
        float n2,factor;
        int nn2;
        float *fire1, *fire2;

        val = valmin + (float)n*(valmax-valmin)/(float)(MAXSMOKERGB-1);
        if(firecolormap_type==FIRECOLORMAP_CONSTRAINT){
          if(val<=valcut){
            n2 = 1+127*(val-valmin)/(valcut-valmin);
          }
          else{
            n2 = 128 + 126*(val-valcut)/(valmax-valcut);
          }
        }
        else{
          n2 = 1.0+253.0*(val-valmin)/(valmax-valmin);
        }
        nn2 = (int)n2;
        nn2 = CLAMP(nn2,1,253);
        factor = n2 - nn2;
        factor = CLAMP(factor,0.0,1.0);
        fire1 = fire_cb + 3*nn2;
        fire2 = fire1 + 3;
        rgb_colormap[4*n]  =(1.0-factor)*fire1[0]+factor*fire2[0];
        rgb_colormap[4*n+1]=(1.0-factor)*fire1[1]+factor*fire2[1];
        rgb_colormap[4*n+2]=(1.0-factor)*fire1[2]+factor*fire2[2];
        if(alpha[n]==0){
          rgb_colormap[4*n+3]=0.0;
        }
        else{
          rgb_colormap[4*n+3]=transparent_level_local;
        }
      }
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
  UpdateTexturebar();
}

/* ------------------ UpdateRGBColors ------------------------ */

void UpdateRGBColors(int colorbar_index){

  int n,nn;
  int i,j;
  float *rgb2ptr;
  int cci;
  meshdata *meshi;
  int vent_offset, outline_offset;
  facedata *facej;
  float transparent_level_local=1.0;

  if(use_transparency_data==1)transparent_level_local=transparent_level;

  InitCadColors();
  InitRGB();
  nrgb_full = MAXRGB;
  for(n=0;n<nrgb_full;n++){
    rgb_trans[4*n]=0.0;
    rgb_trans[4*n+1]=0.0;
    rgb_trans[4*n+2]=0.0;
    rgb_trans[4*n+3]=(float)n/(float)(nrgb_full-1);
  }
  if(colorbarinfo!=NULL){
    unsigned char *alpha;
    colorbardata *cbi;

    cbi = colorbarinfo + colorbartype;

    alpha = colorbarinfo[colorbartype].alpha;
    for(n=0;n<nrgb_full;n++){
      rgb_full[n][0]=cbi->colorbar[3*n];
      rgb_full[n][1]=cbi->colorbar[3*n+1];
      rgb_full[n][2]=cbi->colorbar[3*n+2];
      if(alpha[n]==0){
        rgb_full[n][3]=0.0;
      }
      else{
        rgb_full[n][3]=transparent_level_local;
      }
    }
    UpdateSmokeColormap(RENDER_SLICE);
    UpdateSmokeColormap(RENDER_VOLUME);
  }
  else{
    for(n=0;n<nrgb_full;n++){
      rgb_full[n][0]=(float)n/(float)(nrgb_full);
      rgb_full[n][1]=(float)n/(float)(nrgb_full);
      rgb_full[n][2]=(float)n/(float)(nrgb_full);
      rgb_full[n][3]=transparent_level_local;
    }
  }
  if(contour_type==LINE_CONTOURS){
    for(n=0;n<nrgb_full;n++){
      rgb_full2[n][3]=rgb_full[n][3];
      rgb_full[n][3]=0;
    }
    for(n=0;n<11;n++){
      int nnm1,nnp0,nnp1;

      if(n==0){
        nnp0=1;
      }
      else if(n==10){
        nnp0=254;
      }
      else{
        nnp0=1+n*25.4;
      }
      nnm1=nnp0-1;
      nnp1=nnp0+1;
      rgb_full[nnm1][3]=rgb_full2[nnm1][3];
      rgb_full[nnp0][3]=rgb_full2[nnp0][3];
      rgb_full[nnp1][3]=rgb_full2[nnp1][3];
    }
  }
  if(contour_type==STEPPED_CONTOURS){
    int index[11];

    for(n=0;n<10;n++){
      index[n]=n*25.4;
    }
    index[10]=nrgb_full;
    for(n=0;n<10;n++){
      int mid;

      mid = (index[n]+index[n+1])/2;
      for(i=index[n];i<index[n+1];i++){
        rgb_full[i][0]=rgb_full[mid][0];
        rgb_full[i][1]=rgb_full[mid][1];
        rgb_full[i][2]=rgb_full[mid][2];
        rgb_full[i][3]=rgb_full[mid][3];
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
  global_colorbar_index=colorbar_index;
  if(colorbar_index>=0){
    float highlight_black[3]={0.0,0.0,0.0},highlight_red[3]={1.0,0.0,0.0},*highlight_color;
    int cbmin, cbmax;

    valindex = global_colorbar_index;
    if(valindex<0)valindex=0;
    if(valindex>255)valindex=255;
    cci = colorbar_index;
    if(setbw==1){
      highlight_color=highlight_red;
    }
    else{
      highlight_color=highlight_black;
    }
    cbmin = cci-colorband;
    cbmax = cci+colorband;
    if(cbmin<0){
      cbmax = cbmax - cbmin;
      cbmin = 0;
    }
    if(cbmax>255){
      cbmin = cbmin - (cbmax-255);
      cbmax = 255;
    }
    for(n=cbmin;n<cbmax+1;n++){
      rgb_full[n][0]=highlight_color[0];
      rgb_full[n][1]=highlight_color[1];
      rgb_full[n][2]=highlight_color[2];
      rgb_full[n][3]=transparent_level_local;
    }
  }
  if(show_extreme_mindata==1){
    rgb_full[0][0]=rgb_below_min[0]/255.0;
    rgb_full[0][1]=rgb_below_min[1]/255.0;;
    rgb_full[0][2]=rgb_below_min[2]/255.0;;
  }
  if(show_extreme_maxdata==1){
    rgb_full[255][0]=rgb_above_max[0]/255.0;
    rgb_full[255][1]=rgb_above_max[1]/255.0;
    rgb_full[255][2]=rgb_above_max[2]/255.0;
  }
  if(rgb2_ini!=NULL){
    rgb2ptr=rgb2_ini;
  }
  else{
    rgb2ptr=&(rgb2[0][0]);
  }
  if(colorbar_index!=0){
    for(n=0;n<nrgb;n++){
      nn=n*(nrgb_full-1)/(nrgb-1);
      rgb[n][0] = rgb_full[nn][0];
      rgb[n][1] = rgb_full[nn][1];
      rgb[n][2] = rgb_full[nn][2];
      rgb[n][3] = transparent_level_local;
    }
  }
  for(n=nrgb;n<nrgb+nrgb2;n++){
    rgb[n][0]=rgb2ptr[3*(n-nrgb)];
    rgb[n][1]=rgb2ptr[3*(n-nrgb)+1];
    rgb[n][2]=rgb2ptr[3*(n-nrgb)+2];
    rgb[n][3]=transparent_level_local;
  }
  rgb_white=nrgb;
  rgb_yellow=nrgb+1;
  rgb_blue=nrgb+2;
  rgb_red=nrgb+3;

  if(background_flip==0){
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
    rgb[rgb_white][0]=0.0;  //xxx fix or take out
    rgb[rgb_white][1]=0.0;
    rgb[rgb_white][2]=0.0;
    rgb[rgb_black][0]=1.0;
    rgb[rgb_black][1]=1.0;
    rgb[rgb_black][2]=1.0;
  }
  for(i=0;i<nmeshes;i++){
    meshi=meshinfo + i;
    vent_offset = 6*meshi->nbptrs;
    outline_offset = vent_offset + meshi->nvents;
    for(j=outline_offset;j<outline_offset+6;j++){
      facej = meshi->faceinfo + j;
      facej->color=foregroundcolor;
    }
  }
  UpdateChopColors();
  InitCadColors();
  UpdateTexturebar();
}

/* ------------------ UpdateChopColors ------------------------ */

void UpdateChopColors(void){
  int i;
  int ichopmin=0,ichopmax=nrgb_full;
#define NCHOP 8
  int ii;
  float transparent_level_local=1.0;

  if(use_transparency_data==1)transparent_level_local=transparent_level;

  for(i=0;i<nrgb_full;i++){
    rgb_iso[4*i]=rgb_full[i][0];
    rgb_iso[4*i+1]=rgb_full[i][1];
    rgb_iso[4*i+2]=rgb_full[i][2];
    if(rgb_full[i][3]>0.001){
      rgb_iso[4*i+3]=transparent_level_local;
    }
    else{
      rgb_iso[4*i+3]=0.0;
    }

    rgb_slice[4*i]=rgb_full[i][0];
    rgb_slice[4*i+1]=rgb_full[i][1];
    rgb_slice[4*i+2]=rgb_full[i][2];
    if(rgb_full[i][3]>0.001){
      rgb_slice[4*i+3]=transparent_level_local;
    }
    else{
      rgb_slice[4*i+3]=0.0;
    }

    rgb_terrain2[4 * i] = rgb_full[i][0];
    rgb_terrain2[4 * i + 1] = rgb_full[i][1];
    rgb_terrain2[4 * i + 2] = rgb_full[i][2];
    if(rgb_full[i][3] > 0.001){
      rgb_terrain2[4 * i + 3] = transparent_level_local;
    }
    else{
      rgb_terrain2[4 * i + 3] = 0.0;
    }
    if(show_zlevel == 1){
      int ilevel;
      float dz;

      dz = (terrain_zmax - terrain_zmin)*geom_vert_exag;
      if(ABS(dz)<0.01)dz=1;

      ilevel = 255 * (terrain_zlevel - terrain_zmin) / dz;
      if(ABS(ilevel - i) < 3){
        rgb_terrain2[4 * i] = 0;
        rgb_terrain2[4 * i + 1] = 0;
        rgb_terrain2[4 * i + 2] = 0;
      }
    }

    rgb_part[4 * i] = rgb_full[i][0];
    rgb_part[4*i+1]=rgb_full[i][1];
    rgb_part[4*i+2]=rgb_full[i][2];
    rgb_part[4*i+3]=rgb_full[i][3];
    if(rgb_full[i][3]>0.001){
      rgb_part[4*i+3]=transparent_level_local;
    }
    else{
      rgb_part[4*i+3]=0.0;
    }

    rgb_plot3d[4*i]=rgb_full[i][0];
    rgb_plot3d[4*i+1]=rgb_full[i][1];
    rgb_plot3d[4*i+2]=rgb_full[i][2];
    rgb_plot3d[4*i+3]=rgb_full[i][3];
    if(rgb_full[i][3]>0.001){
      rgb_plot3d[4*i+3]=transparent_level_local;
    }
    else{
      rgb_plot3d[4*i+3]=0.0;
    }

    rgb_patch[4*i]=rgb_full[i][0];
    rgb_patch[4*i+1]=rgb_full[i][1];
    rgb_patch[4*i+2]=rgb_full[i][2];
    if(rgb_full[i][3]>0.001){
      rgb_patch[4*i+3]=transparent_level_local;
    }
    else{
      rgb_patch[4*i+3]=0.0;
    }
  }
  {
    float smin, smax;

    smin = boundarylevels256[0];
    smax = boundarylevels256[255];

    if(setpatchchopmin==1){
      ichopmin=nrgb_full*(patchchopmin-smin)/(smax-smin);
      if(ichopmin<0)ichopmin=0;
      if(ichopmin>nrgb_full-1)ichopmin=nrgb_full-1;
      for(i=0;i<ichopmin;i++){
        rgb_patch[4*i+3]=0.0;
      }
      for(i=ichopmin-NCHOP;i<ichopmin;i++){
        if(i<=0)continue;
        if(i>nrgb_full-1)continue;
        ii = i - (ichopmin-NCHOP);
        if(ii>NCHOP-1)continue;
        rgb_patch[4*i+3]=transparent_level_local*(float)ii/(float)(NCHOP-1);
      }
    }
    if(setpatchchopmax==1){
      ichopmax=nrgb_full*(patchchopmax - smin)/(smax-smin);
      if(ichopmax<0)ichopmax=0;
      if(ichopmax>nrgb_full-1)ichopmax=nrgb_full-1;
      for(i=ichopmax;i<nrgb_full;i++){
        rgb_patch[4*i+3]=0.0;
      }
      for(i=ichopmax;i<ichopmax+NCHOP;i++){
        if(i<=0)continue;
        if(i>nrgb_full-1)continue;
        ii = NCHOP-1-(i - ichopmax);
        if(ii>NCHOP-1)continue;
        rgb_patch[4*i+3]=transparent_level_local*(float)ii/(float)(NCHOP-1);
      }
    }
  }
  if(slicebounds!=NULL&&islicetype!=-1){
    float smin, smax;

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
        rgb_slice[4*i+3]=transparent_level_local*(float)ii/(float)(NCHOP-1);
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
        rgb_slice[4*i+3]=transparent_level_local*(float)ii/(float)(NCHOP-1);
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
        rgb_part[4*i+3]=transparent_level_local*(float)ii/(float)(NCHOP-1);
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
        rgb_part[4*i+3]=transparent_level_local*(float)ii/(float)(NCHOP-1);
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
        rgb_plot3d[4*i+3]=transparent_level_local*(float)ii/(float)(NCHOP-1);
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
        rgb_plot3d[4*i+3]=transparent_level_local*(float)ii/(float)(NCHOP-1);
      }
    }
  }
  for(i=0;i<npartinfo;i++){
    partdata *parti;

    parti = partinfo + i;
    if(parti->loaded==0)continue;
    AdjustPart5Chops(parti);
  }
  UpdateTexturebar();
}

/* ------------------ GetRGB ------------------------ */

void GetRGB(unsigned int val, unsigned char *rr, unsigned char *gg, unsigned char *bb){
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

/* ------------------ color2bw ------------------------ */

float color2bw(const float *color){
  float returnval;

  returnval = 0.299*color[0] + 0.587*color[1] + 0.114*color[2];
  return returnval;
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
    firstcolor->bw_color[0] = color2bw(color);
    firstcolor->bw_color[1] = firstcolor->bw_color[0];
    firstcolor->bw_color[2] = firstcolor->bw_color[0];
    firstcolor->bw_color[3] = color[3];
    firstcolor->nextcolor=NULL;
    return firstcolor->color;
  }
  oldlastcolor = firstcolor;
  for(colorptr = firstcolor; colorptr!=NULL; colorptr = colorptr->nextcolor){
    oldlastcolor=colorptr;
    if(ABS(colorptr->color[0]-color[0])>0.0001)continue;
    if(ABS(colorptr->color[1]-color[1])>0.0001)continue;
    if(ABS(colorptr->color[2]-color[2])>0.0001)continue;
    if(ABS(colorptr->color[3]-color[3])>0.0001)continue;
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

/* ------------------ colorconvert ------------------------ */

void colorconvert(int flag){
  colordata *colorptr;
  extern colordata *firstcolor;

  switch(flag){
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
