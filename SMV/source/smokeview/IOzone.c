// $Date: 2015-04-16 22:07:09 -0400 (Thu, 16 Apr 2015) $ 
// $Revision: 22416 $
// $Author: gforney $

// svn revision character string
char IOzone_revision[]="$Revision: 22416 $";

#include "options.h"
#include <stdio.h>  
#include <stdlib.h> 
#include <string.h>
#include <math.h>
#include "glew.h"

#include "string_util.h"
#include "update.h"
#include "smokeviewvars.h"

void drawtrunccone(float d1, float d2, float height, unsigned char *rgbcolor);

/* ------------------ getzonesizecsv ------------------------ */

void getzonesizecsv(int *nzone_times_local, int *nroom, int *nfires_local, int *nhvents, int *nvvents, int *error){
   devicedata *dev;
   int nr,nf,nv;
   int i;
   char label[100];

   
   *error=0;
   nr=0;
   for(i=0;i<ndeviceinfo;i++){
     sprintf(label,"ULT_%i",i+1);
     dev=getdevice(label,-1);
     if(dev==NULL)break;
     *nzone_times_local=dev->nvals;
     nr++;
   }
   *nroom=nr;

   nv=0;
   for(i=0;i<ndeviceinfo;i++){
     sprintf(label,"HVENT_%i",i+1);
     dev=getdevice(label,-1);
     if(dev==NULL)break;
     nv++;
   }
   *nhvents=nv;

   nv=0;
   for(i=0;i<ndeviceinfo;i++){
     sprintf(label,"VVENT_%i",i+1);
     dev=getdevice(label,-1);
     if(dev==NULL)break;
     nv++;
   }
   *nvvents=nv;

   nf=0;
   for(i=0;i<ndeviceinfo;i++){
     sprintf(label,"HRR_%i",i+1);
     dev=getdevice(label,-1);
     if(dev==NULL)break;
     nf++;
   }
   *nfires_local=nf;
}

#ifdef pp_ZONEVENT
#define FREEZONEMEM \
  FREEMEMORY(zonehvents_devs); \
  FREEMEMORY(zonehvents_devs); \
  FREEMEMORY(zoneslab_n_devs); \
  FREEMEMORY(zoneslab_T_devs); \
  FREEMEMORY(zoneslab_F_devs); \
  FREEMEMORY(zoneslab_YB_devs); \
  FREEMEMORY(zoneslab_YT_devs); \
  FREEMEMORY(zonevvents_devs); \
  FREEMEMORY(zoneqfire_devs); \
  FREEMEMORY(zonepr_devs); \
  FREEMEMORY(zoneylay_devs); \
  FREEMEMORY(zonetl_devs); \
  FREEMEMORY(zonetu_devs); \
  FREEMEMORY(zoneodl_devs); \
  FREEMEMORY(zoneodu_devs); \
  FREEMEMORY(zonefheight_devs); \
  FREEMEMORY(zonefarea_devs); \
  FREEMEMORY(zonefbase_devs)

#define GETZONEDEV(zonedev)\
      zonedev = getdevice(label, -1);\
      if(zonedev == NULL || zonedev->nvals != nzone_times_local){\
        *error=1;\
        FREEZONEMEM;\
        return;\
      }\
      zonedev->in_zone_csv = 1
#else
#define FREEZONEMEM \
  FREEMEMORY(zonehvents_devs); \
  FREEMEMORY(zonevvents_devs); \
  FREEMEMORY(zoneqfire_devs); \
  FREEMEMORY(zonepr_devs); \
  FREEMEMORY(zoneylay_devs); \
  FREEMEMORY(zonetl_devs); \
  FREEMEMORY(zonetu_devs); \
  FREEMEMORY(zoneodl_devs); \
  FREEMEMORY(zoneodu_devs); \
  FREEMEMORY(zonefheight_devs); \
  FREEMEMORY(zonefarea_devs); \
  FREEMEMORY(zonefbase_devs)
#endif

#ifdef pp_ZONEVENT
/* ------------------ getzonedatacsv ------------------------ */

#define GET_VENTSLAB(slabtype,ventbeg,ventend,have_ventslab_flow)\
  have_ventslab_flow = 0;\
  for(i = ventbeg; i < ventend; i++){\
    char label[100];\
    int islab;\
\
    sprintf(label, "%s_%i",#slabtype, i + 1);\
    zoneslab_n_devs[i] = getdevice(label, -1);\
    if(zoneslab_n_devs[i] != NULL){\
      have_ventslab_flow = 1;\
      break;\
    }\
    for(islab = 0; islab < MAXSLABS; islab++){\
      int idev;\
\
      idev = MAXSLABS * i + islab;\
      sprintf(label, "%s_%i_%i", #slabtype,i + 1, islab + 1);\
      zoneslab_T_devs[idev] = getdevice(label, -1);\
      if(zoneslab_T_devs[idev] != NULL)have_ventslab_flow = 1;\
\
      sprintf(label, "%sF_%i_%i",#slabtype, i + 1, islab + 1);\
      zoneslab_F_devs[idev] = getdevice(label, -1);\
      if(zoneslab_T_devs[idev] != NULL)have_ventslab_flow = 1;\
\
      sprintf(label, "%sYB_%i_%i",#slabtype, i + 1, islab + 1);\
      zoneslab_YB_devs[idev] = getdevice(label, -1);\
      if(zoneslab_T_devs[idev] != NULL)have_ventslab_flow = 1;\
\
      sprintf(label, "%sYT_%i_%i",#slabtype, i + 1, islab + 1);\
      zoneslab_YT_devs[idev] = getdevice(label, -1);\
      if(zoneslab_T_devs[idev] != NULL)have_ventslab_flow = 1;\
    }\
    if(have_ventslab_flow == 1)break;\
  }\
  if(have_ventslab_flow == 1){\
    for(i = ventbeg; ventend < nzhvents; i++){\
      char label[100];\
      int islab;\
\
      sprintf(label, "%s_%i", #slabtype,i + 1);\
      GETZONEDEV(zoneslab_n_devs[i]);\
      for(islab = 0; islab < MAXSLABS; islab++){\
        int idev;\
\
        idev = MAXSLABS * i + islab;\
        sprintf(label, "%sT_%i_%i", #slabtype,i + 1, islab + 1);\
        GETZONEDEV(zoneslab_T_devs[idev]);\
\
        sprintf(label, "%sF_%i_%i", #slabtype,i + 1, islab + 1);\
        GETZONEDEV(zoneslab_F_devs[idev]);\
\
        sprintf(label, "%sYB_%i_%i",#slabtype, i + 1, islab + 1);\
        GETZONEDEV(zoneslab_YB_devs[idev]);\
\
        sprintf(label, "%sYT_%i_%i",#slabtype, i + 1, islab + 1);\
        GETZONEDEV(zoneslab_YT_devs[idev]);\
      }\
    }\
  }
#endif  

/* ------------------ getzonedatacsv ------------------------ */

void getzonedatacsv(int nzone_times_local, int nrooms_local, int nfires_local, 
                    float *zone_times_local, float *zoneqfire_local, float *zonefheight_local, float *zonefbase_local, float *zonefdiam_local,
                    float *zonepr_local, float *zoneylay_local,  float *zonetl_local, float *zonetu_local,
                    float **zoneodlptr, float **zoneoduptr, float *zonehvents_local, float *zonevvents_local,
#ifdef pp_ZONEVENT                    
                    int *zoneslab_n_local, float *zoneslab_T_local, float *zoneslab_F_local, float *zoneslab_YB_local, float *zoneslab_YT_local,
#endif
                    int *error){
  int i,ii,iif, use_od=1, iihv, iivv;
  devicedata **zoneqfire_devs=NULL;
  devicedata **zonepr_devs=NULL, **zoneylay_devs=NULL, **zonetl_devs=NULL, **zonetu_devs=NULL, **zoneodl_devs=NULL, **zoneodu_devs=NULL;
  devicedata **zonefheight_devs=NULL, **zonefbase_devs=NULL, **zonefarea_devs=NULL;
  devicedata **zonehvents_devs=NULL, **zonevvents_devs=NULL;
#ifdef pp_ZONEVENT  
  devicedata **zoneslab_n_devs = NULL, **zoneslab_T_devs = NULL, **zoneslab_F_devs = NULL, **zoneslab_YB_devs = NULL, **zoneslab_YT_devs = NULL;
#endif  
  float *zoneodl_local, *zoneodu_local;
  float *times_local;

  *error=0;
  if(nfires_local>0){
    NewMemory((void **)&zoneqfire_devs,nfires_local*sizeof(devicedata *));
    NewMemory((void **)&zonefheight_devs,nfires_local*sizeof(devicedata *));
    NewMemory((void **)&zonefbase_devs,nfires_local*sizeof(devicedata *));
    NewMemory((void **)&zonefarea_devs,nfires_local*sizeof(devicedata *));
  }

  if(nrooms_local>0){
    NewMemory((void **)&zonepr_devs,nrooms_local*sizeof(devicedata *));
    NewMemory((void **)&zoneylay_devs,nrooms_local*sizeof(devicedata *));
    NewMemory((void **)&zonetl_devs,nrooms_local*sizeof(devicedata *));
    NewMemory((void **)&zonetu_devs,nrooms_local*sizeof(devicedata *));
    NewMemory((void **)&zoneodl_devs,nrooms_local*sizeof(devicedata *));
    NewMemory((void **)&zoneodu_devs,nrooms_local*sizeof(devicedata *));
  }

  if(nzhvents>0){
    NewMemory((void **)&zonehvents_devs,nzhvents*sizeof(devicedata *));
#ifdef pp_ZONEVENT    
    NewMemory((void **)&zoneslab_n_devs,  nzhvents*sizeof(devicedata *));
    NewMemory((void **)&zoneslab_T_devs,  MAXSLABS*nzhvents*sizeof(devicedata *));
    NewMemory((void **)&zoneslab_F_devs,  MAXSLABS*nzhvents*sizeof(devicedata *));
    NewMemory((void **)&zoneslab_YB_devs, MAXSLABS*nzhvents*sizeof(devicedata *));
    NewMemory((void **)&zoneslab_YT_devs, MAXSLABS*nzhvents*sizeof(devicedata *));
#endif    
  }

  if(nzvvents>0){
    NewMemory((void **)&zonevvents_devs,nzvvents*sizeof(devicedata *));
  }

  for(i=0;i<nrooms_local;i++){
    char label[100];

    sprintf(label,"PRS_%i",i+1);
    zonepr_devs[i]=getdevice(label,-1);
    if(zonepr_devs[i]==NULL||zonepr_devs[i]->nvals!=nzone_times_local){
      *error = 1;
      FREEZONEMEM;
      return;
    }
    zonepr_devs[i]->in_zone_csv=1;

    sprintf(label,"HGT_%i",i+1);
    zoneylay_devs[i]=getdevice(label,-1);
    if(zoneylay_devs[i]==NULL||zoneylay_devs[i]->nvals!=nzone_times_local){
      zoneylay_devs[i]=NULL;
    }
    else{
      zoneylay_devs[i]->in_zone_csv=1;
    }

    sprintf(label,"LLT_%i",i+1);
    zonetl_devs[i]=getdevice(label,-1);
    if(zonetl_devs[i]==NULL||zonetl_devs[i]->nvals!=nzone_times_local){
      zonetl_devs[i]=NULL;
    }
    else{
      zonetl_devs[i]->in_zone_csv=1;
    }

    sprintf(label,"ULT_%i",i+1);
    zonetu_devs[i]=getdevice(label,-1);
    if(zonetu_devs[i]==NULL||zonetu_devs[i]->nvals!=nzone_times_local){
      *error = 1;
      FREEZONEMEM;
      return;
    }
    zonetu_devs[i]->in_zone_csv=1;

    sprintf(label,"ULOD_%i",i+1);
    zoneodu_devs[i]=getdevice(label,-1);
    if(zoneodu_devs[i]==NULL){
      use_od=0;
    }
    else{
      zoneodu_devs[i]->in_zone_csv=1;
    }

    sprintf(label,"LLOD_%i",i+1);
    zoneodl_devs[i]=getdevice(label,-1);
    if(zoneodl_devs[i]==NULL)zoneodl_devs[i]=zoneodu_devs[i];
    if(zoneodl_devs[i]==NULL){
      use_od=0;
    }
    else{
      zoneodl_devs[i]->in_zone_csv=1;
    }
  }
  if(use_od==1){
    NewMemory((void **)&zoneodl_local     ,nrooms_local*nzone_times_local*sizeof(float));
    NewMemory((void **)&zoneodu_local     ,nrooms_local*nzone_times_local*sizeof(float));
  }
  else{
    zoneodl_local=NULL;
    zoneodu_local=NULL;
  }
  *zoneodlptr=zoneodl_local;
  *zoneoduptr=zoneodu_local;

  for(i=0;i<nfires_local;i++){
    char label[100];

    sprintf(label,"HRR_%i",i+1);
    zoneqfire_devs[i]=getdevice(label,-1);
    if(zoneqfire_devs[i]==NULL||zoneqfire_devs[i]->nvals!=nzone_times_local){
      *error=1;
      FREEZONEMEM;
      return;
    }
    zoneqfire_devs[i]->in_zone_csv=1;

    sprintf(label,"FLHGT_%i",i+1);
    zonefheight_devs[i]=getdevice(label,-1);
    if(zonefheight_devs[i]==NULL||zonefheight_devs[i]->nvals!=nzone_times_local){
      *error=1;
      FREEZONEMEM;
      return;
    }
    zonefheight_devs[i]->in_zone_csv=1;

    sprintf(label,"FBASE_%i",i+1);
    zonefbase_devs[i]=getdevice(label,-1);
    if(zonefbase_devs[i]==NULL||zonefbase_devs[i]->nvals!=nzone_times_local){
      *error=1;
      FREEZONEMEM;
      return;
    }
    zonefbase_devs[i]->in_zone_csv=1;

    sprintf(label,"FAREA_%i",i+1);
    zonefarea_devs[i]=getdevice(label,-1);
    if(zonefarea_devs[i]==NULL||zonefarea_devs[i]->nvals!=nzone_times_local){
      *error=1;
      FREEZONEMEM;
      return;
    }
    zonefarea_devs[i]->in_zone_csv=1;
  }

  for(i = 0; i < nzhvents; i++){
    char label[100];

    sprintf(label, "HVENT_%i", i + 1);
    zonehvents_devs[i] = getdevice(label, -1);
    if(zonehvents_devs[i] == NULL || zonehvents_devs[i]->nvals != nzone_times_local){
      *error=1;
      FREEZONEMEM;
      return;
    }
    zonehvents_devs[i]->in_zone_csv = 1;
  }

#ifdef pp_ZONEVENT
  GET_VENTSLAB("HSLAB",0,nzhvents,have_hventslab_flow);
#endif    

  for(i=0;i<nzvvents;i++){
    char label[100];

    sprintf(label,"VVENT_%i",i+1);
    zonevvents_devs[i]=getdevice(label,-1);
    if(zonevvents_devs[i]==NULL||zonevvents_devs[i]->nvals!=nzone_times_local){
      *error = 1;
      FREEZONEMEM;
      return;
    }
    zonevvents_devs[i]->in_zone_csv=1;
  }

  ii=0;
  iif=0;
  iihv=0;
  iivv=0;
  times_local = zonepr_devs[0]->times;

#ifdef pp_ZONEVENT
  maxslabflow = 0.0;
#endif
  for(i=0;i<nzone_times_local;i++){
    int j, ivent;

    zone_times_local[i]=times_local[i];
    for(j=0;j<nrooms_local;j++){
      zonepr_local[ii]=zonepr_devs[j]->vals[i];
      if(zoneylay_devs[j]==NULL){
        zoneylay_local[ii]=0.0;
      }
      else{
        zoneylay_local[ii]=zoneylay_devs[j]->vals[i];
      }
      zonetu_local[ii]=zonetu_devs[j]->vals[i];
      if(zonetl_devs[j]!=NULL){
        zonetl_local[ii]=zonetl_devs[j]->vals[i];
      }
      else{
        zonetl_local[ii]=zonetu_devs[j]->vals[i];
      }

      if(use_od==1){
        zoneodu_local[ii]=zoneodu_devs[j]->vals[i];
        if(zoneodl_devs[j]!=NULL){
          zoneodl_local[ii]=zoneodl_devs[j]->vals[i];
        }
        else{
          zoneodl_local[ii]=zoneodu_devs[j]->vals[i];
        }
      }
      ii++;
    }
    for(j=0;j<nfires_local;j++){
      float area, diam;

      zoneqfire_local[iif]=1000.0*zoneqfire_devs[j]->vals[i];
      zonefheight_local[iif]=zonefheight_devs[j]->vals[i];
      area=zonefarea_devs[j]->vals[i];
      diam=2.0*sqrt(area/PI);
      if(diam<0.0001){
        diam=SCALE2SMV(0.1);
      }
      zonefdiam_local[iif]=diam;
      zonefbase_local[iif]=zonefbase_devs[j]->vals[i];
      iif++;
    }
    for(ivent=0;ivent<nzhvents;ivent++){
#ifdef pp_ZONEVENT    
      int islab;
#endif      

      zonehvents_local[iihv] = zonehvents_devs[ivent]->vals[i];
#ifdef pp_ZONEVENT      
      if(zoneslab_n_devs[ivent]!=NULL)zoneslab_n_local[iihv] = (int)(zoneslab_n_devs[ivent]->vals[i]+0.1);
      for(islab = 0; islab<MAXSLABS; islab++){
        int idev, ival;

        idev = MAXSLABS * ivent + islab;
        ival = MAXSLABS * iihv + islab;
        if(zoneslab_T_devs[idev]!=NULL)zoneslab_T_local[ival] =  zoneslab_T_devs[idev]->vals[i];
        if(zoneslab_F_devs[idev] != NULL){
          float slabflow;

          slabflow = zoneslab_F_devs[idev]->vals[i];
          maxslabflow = MAX(ABS(slabflow), maxslabflow);
          zoneslab_F_local[ival] = slabflow;
        }
        if(zoneslab_YB_devs[idev]!=NULL)zoneslab_YB_local[ival] = zoneslab_YB_devs[idev]->vals[i];
        if(zoneslab_YT_devs[idev]!=NULL)zoneslab_YT_local[ival] = zoneslab_YT_devs[idev]->vals[i];
      }
#endif      
      iihv++;
    }
    for(j=0;j<nzvvents;j++){
      zonevvents_local[iivv] = zonevvents_devs[j]->vals[i];
      iivv++;
    }
  }
}

/* ------------------ getsmokedir ------------------------ */

void getzonesmokedir(float *mm){
    /*
      ( m0 m4 m8  m12 ) (x)    (0)
      ( m1 m5 m9  m13 ) (y)    (0)
      ( m2 m6 m10 m14 ) (z)  = (0)
      ( m3 m7 m11 m15 ) (1)    (1)

       ( m0 m4  m8 )      (m12)
   Q=  ( m1 m5  m9 )  u = (m13)
       ( m2 m6 m10 )      (m14)
      
      (Q   u) (x)     (0)      
      (v^T 1) (y)   = (1)
       
      m3=m7=m11=0, v^T=0, y=1   Qx+u=0 => x=-Q^Tu
    */
  int i,ii,j;
  float norm[3];
  float eyedir[3];
  float cosdir;
  float angles[7];

  xyzeyeorig[0] = -(mm[0]*mm[12]+mm[1]*mm[13]+ mm[2]*mm[14])/mscale[0];
  xyzeyeorig[1] = -(mm[4]*mm[12]+mm[5]*mm[13]+ mm[6]*mm[14])/mscale[1];
  xyzeyeorig[2] = -(mm[8]*mm[12]+mm[9]*mm[13]+mm[10]*mm[14])/mscale[2];
  
  for(j=0;j<nrooms;j++){
    roomdata *roomj;
    
    roomj = roominfo + j;

    roomj->zoneinside=0;
    if(
      xyzeyeorig[0]>roomj->x0&&xyzeyeorig[0]<roomj->x1&&
      xyzeyeorig[1]>roomj->y0&&xyzeyeorig[1]<roomj->y1&&
      xyzeyeorig[2]>roomj->z0&&xyzeyeorig[2]<roomj->z1
      ){
      for(i=-3;i<=3;i++){
        if(i==0)continue;
        roomj->drawsides[i+3]=1;
      }
      roomj->zoneinside=1;
      continue;
    }

    for(i=-3;i<=3;i++){
      if(i==0)continue;
      ii = i;
      if(i<0)ii=-i;
      norm[0]=0.0;
      norm[1]=0.0;
      norm[2]=0.0;
      switch(ii){
      case 1:
        if(i<0){
          norm[0]=-1.0;
          eyedir[0]=roomj->x0;
        }
        else{
          norm[0]=1.0;
          eyedir[0]=roomj->x0+roomj->dx;
        }
        eyedir[1]=roomj->y0+roomj->dy/2.0;
        eyedir[2]=roomj->z0+roomj->dz/2.0;
        break;
      case 2:
        eyedir[0]=roomj->x0+roomj->dx/2.0;
        if(i<0){
          norm[1]=-1.0;
          eyedir[1]=roomj->y0;
        }
        else{
          norm[1]=1.0;
          eyedir[1]=roomj->y0+roomj->dy;
        }
        eyedir[2]=roomj->z0+roomj->dz/2.0;
        break;
      case 3:
        eyedir[0]=roomj->x0+roomj->dx/2.0;
        eyedir[1]=roomj->y0+roomj->dy/2.0;
        if(i<0){
          norm[2]=-1.0;
          eyedir[2]=roomj->z0;
        }
        else{
          norm[2]=1.0;
          eyedir[2]=roomj->z0+roomj->dz;
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
      }
      eyedir[0]=xyzeyeorig[0]-eyedir[0];
      eyedir[1]=xyzeyeorig[1]-eyedir[1];
      eyedir[2]=xyzeyeorig[2]-eyedir[2];
      normalize(eyedir,3);
      cosdir = (eyedir[0]*norm[0]+eyedir[1]*norm[1]+eyedir[2]*norm[2]);
      if(cosdir>1.0)cosdir=1.0;
      if(cosdir<-1.0)cosdir=-1.0;
      cosdir=acos(cosdir)*RAD2DEG;
      if(cosdir<0.0)cosdir=-cosdir;
      angles[3+i]=cosdir;
    }
    for(i=-3;i<=3;i++){
      if(i==0)continue;
      if(angles[i+3]<90.0){
        roomj->drawsides[i+3]=1;
      }
      else{
        roomj->drawsides[i+3]=0;
      }
    }
  }
}

/* ------------------ readzone ------------------------ */

void readzone(int ifile, int flag, int *errorcode){
  int error,ntotal,i,j,ii;
  int nrooms2,nfires2,nzhvents2,nzvvents2;
  size_t zonefilelen;
  zonedata *zonei;
  char *file;

  *errorcode=0;

  ASSERT(ifile>=0&&ifile<nzoneinfo);
  zonei = zoneinfo + ifile;
  file = zonei->file;
  if(zonei->loaded==0&&flag==UNLOAD)return;
  FREEMEMORY(zonehvents);
  FREEMEMORY(zonevvents);
  FREEMEMORY(zone_times);
  FREEMEMORY(zoneylay);
  FREEMEMORY(zonetl);
  FREEMEMORY(zonetu);
  FREEMEMORY(zonepr);
  FREEMEMORY(hazardcolor);
  FREEMEMORY(zoneqfire);
  FREEMEMORY(zonefheight);
  FREEMEMORY(zonefdiam);
  FREEMEMORY(zonefbase);
  FREEMEMORY(izonetu);
  FREEMEMORY(zoneodl);
  FREEMEMORY(zoneodu);
#ifdef pp_ZONEVENT  
  FREEMEMORY(zoneslab_n);
  FREEMEMORY(zoneslab_T);
  FREEMEMORY(zoneslab_F);
  FREEMEMORY(zoneslab_YB);
  FREEMEMORY(zoneslab_YT);
#endif  

  activezone=NULL;
  if(flag==UNLOAD){
    zonei->loaded=0;
    zonei->display=0;
    nzone_times=0;
    ReadZoneFile=0;
    showzone=0;
    plotstate=getplotstate(DYNAMIC_PLOTS);
    Update_Times();
    updatemenu=1;
    return;
  }
  zonefilelen = strlen(file);
  if(zonei->csv==1){
    read_device_data(zonei->file,CSV_CFAST,UNLOAD);
    read_device_data(zonei->file,CSV_CFAST,LOAD);
    getzonesizecsv(&nzone_times,&nrooms2,&nfires2,&nzhvents2,&nzvvents2,&error);
  }
  else{
    FORTgetzonesize(file,&nzone_times,&nrooms2,&nfires2,&error,zonefilelen);
    nzhvents2=nzhvents;
    nzvvents2=nzvvents;
  }
  CheckMemory;
  if(error!=0||nrooms!=nrooms2||nzone_times==0||nzhvents!=nzhvents2||nzvvents!=nzvvents2){
    showzone=0;
    Update_Times();
    ReadZoneFile=0;
    if(nrooms!=nrooms2){
      fprintf(stderr,"*** Error: number of rooms specified in the smv file (%i)\n",nrooms);
      fprintf(stderr,"    not consistent with the number specified in the zone file (%i)\n",nrooms2);
    }
    if(nzhvents!=nzhvents2){
      fprintf(stderr,"*** Error: number of horizontal flow vents specified in the smv file (%i)\n",nzhvents);
      fprintf(stderr,"    not consistent with the number specified in the data file (%i)\n",nzhvents2);
    }
    if(nzvvents!=nzvvents2){
      fprintf(stderr,"*** Error: number of vertical flow vents specified in the smv file (%i)\n",nzvvents);
      fprintf(stderr,"    not consistent with the number specified in the data file (%i)\n",nzvvents2);
    }
    if(nzone_times<=0)fprintf(stderr,"*** Error: The file, %s, contains no data\n",file);
    return;
  }
  FREEMEMORY(zonelonglabels);
  FREEMEMORY(zoneshortlabels);
  FREEMEMORY(zoneunits);
  FREEMEMORY(zonelevels);
  if(NewMemory((void **)&zonelonglabels  ,LABELLEN)==0||
     NewMemory((void **)&zoneshortlabels ,LABELLEN)==0||
     NewMemory((void **)&zoneunits       ,LABELLEN)==0||
     NewMemory((void **)&zonelevels      ,nrgb*sizeof(float))==0){
    *errorcode=1;
    return;
  }
  CheckMemory;

  PRINTF("Loading zone data: %s\n",file);

  ntotal = nrooms*nzone_times;
  nzonetotal=ntotal;

  if(ntotal>0){
    FREEMEMORY(zone_times); 
    FREEMEMORY(zoneylay); 
    FREEMEMORY(zonetl); 
    FREEMEMORY(zonetu); 
    FREEMEMORY(zonepr);
    FREEMEMORY(zoneodl);
    FREEMEMORY(zoneodu);
    if(NewMemory((void **)&zone_times      ,ntotal*sizeof(float))==0||
       NewMemory((void **)&zoneylay   ,ntotal*sizeof(float))==0||
       NewMemory((void **)&zonetu     ,ntotal*sizeof(float))==0||
       NewMemory((void **)&zonetl     ,ntotal*sizeof(float))==0||
       NewMemory((void **)&zonepr     ,ntotal*sizeof(float))==0||
       NewMemory((void **)&hazardcolor,ntotal*sizeof(unsigned char))==0){
      *errorcode=1;
      return;
    }
    FREEMEMORY(zonehvents); 
#ifdef pp_ZONEVENT    
    FREEMEMORY(zoneslab_n);
    FREEMEMORY(zoneslab_T);
    FREEMEMORY(zoneslab_F);
    FREEMEMORY(zoneslab_YB);
    FREEMEMORY(zoneslab_YT);
#endif    
    if(nzhvents>0){
      NewMemory((void **)&zonehvents,nzhvents*nzone_times*sizeof(float));
#ifdef pp_ZONEVENT      
      NewMemory((void **)&zoneslab_n, nzone_times*nzhvents*sizeof(int));
      NewMemory((void **)&zoneslab_T, nzone_times*nzhvents*MAXSLABS*sizeof(float));
      NewMemory((void **)&zoneslab_YT, nzone_times*nzhvents*MAXSLABS*sizeof(float));
      NewMemory((void **)&zoneslab_F, nzone_times*nzhvents*MAXSLABS*sizeof(float));
      NewMemory((void **)&zoneslab_YB, nzone_times*nzhvents*MAXSLABS*sizeof(float));
#endif      
    }
    FREEMEMORY(zonevvents); 
    if(nzvvents>0){
      NewMemory((void **)&zonevvents,nzvvents*nzone_times*sizeof(float));
    }
    FREEMEMORY(zoneqfire);
    FREEMEMORY(zonefheight);
    FREEMEMORY(zonefdiam);
    FREEMEMORY(zonefbase);
    if(nfires!=0){
      if(
        NewMemory((void **)&zoneqfire,nfires*nzone_times*sizeof(float))==0||
        NewMemory((void **)&zonefheight,nfires*nzone_times*sizeof(float))==0||
        NewMemory((void **)&zonefdiam,nfires*nzone_times*sizeof(float))==0||
        NewMemory((void **)&zonefbase,nfires*nzone_times*sizeof(float))==0
        ){
        *errorcode=1;
        return;
      }
    }
    else{
      zoneqfire=NULL;
    }
    FREEMEMORY(izonetu);
    if(NewMemory((void **)&izonetu,ntotal*sizeof(unsigned char))==0){
      *errorcode=1;
      return;
    }
  }
  else{
    return;
  }
  CheckMemory;
  if(zonei->csv==1){
#ifdef pp_ZONEVENT
    getzonedatacsv(nzone_times,nrooms,  nfires, zone_times,zoneqfire, zonefheight, zonefbase, zonefdiam,
                   zonepr,zoneylay,zonetl,zonetu,&zoneodl,&zoneodu, zonehvents, zonevvents, 
                   zoneslab_n, zoneslab_T, zoneslab_F, zoneslab_YB, zoneslab_YT,
                   &error);
#else
    getzonedatacsv(nzone_times,nrooms,  nfires, zone_times,zoneqfire, zonefheight, zonefbase, zonefdiam,
                   zonepr,zoneylay,zonetl,zonetu,&zoneodl,&zoneodu, zonehvents, zonevvents, &error);
#endif                   
  }
  else{
    FORTgetzonedata(file,&nzone_times,&nrooms, &nfires, zone_times,zoneqfire,zonepr,zoneylay,zonetl,zonetu,&error,zonefilelen);
  }

  if(zonei->csv==0){
    ii=0;
    for(i=0;i<nzone_times;i++){
      for(j=0;j<nrooms;j++){
        zonetu[ii] = K2C(zonetu[ii]);
        zonetl[ii] = K2C(zonetl[ii]);
        ii++;
      }
    }
  }
  CheckMemory;
  ii = 0;
  for(i=0;i<nzone_times;i++){
    for(j=0;j<nrooms;j++){
      if(zonetu[ii]>=500.0){
        hazardcolor[ii]=RED;
      }
      else{
        if(zonetu[ii]>=50.0){
          if(zoneylay[ii]>1.5){
            hazardcolor[ii]=YELLOW;
          }
          else{
            hazardcolor[ii]=PINK;
          }
        }
        else{
          if(zoneylay[ii]>2.0){
            hazardcolor[ii]=BLUE;
          }
          else{
            hazardcolor[ii]=GREEN;
          }
        }
      }
      zoneylay[ii]=SCALE2SMV(zoneylay[ii]);
      ii++;
    }
  }

  PRINTF("computing zone color levels \n");

  getzoneglobalbounds(zonetu,ntotal,&zoneglobalmin,&zoneglobalmax);
  if(setzonemin==GLOBAL_MIN)zonemin = zoneglobalmin;
  if(setzonemax==GLOBAL_MAX)zonemax = zoneglobalmax;
  if(setzonemin==SET_MIN)zonemin = zoneusermin;
  if(setzonemax==SET_MAX)zonemax = zoneusermax;
  update_glui_zonebounds();
  getZoneColors(zonetu, ntotal, izonetu, zonemin, zonemax, nrgb, nrgb_full, 
    colorlabelzone, zonescale, zonelevels256);

  ReadZoneFile=1;
  visZone=1;
  showzone=1;
  zonei->loaded=1;
  zonei->display=1;
  plotstate=getplotstate(DYNAMIC_PLOTS);
  Update_Times();
  updatemenu=1;
  activezone = zoneinfo + ifile;
  if(nzhvents>0||nzvvents>0){
    PRINTF("computing vent bounds\n");
    getzoneventbounds();
  }
  Idle_CB();

}

/* ------------------ fill_zonedata ------------------------ */

void fill_zonedata(int izone_index){
  float *pr0, *tl0, *tu0, *ylay0, *odl0, *odu0, *hvent0;
  int iroom,ivent;
  roomdata *roomi;
#ifdef pp_ZONEVENT
  int *zoneslab_n0;
  float *zoneslab_T0, *zoneslab_F0, *zoneslab_YB0, *zoneslab_YT0;
#endif

  float CP=1004.0;
  float R=0.4*CP/1.4;

  if(ReadZoneFile==0)return;
  pr0 = zonepr + izone_index*nrooms;
  ylay0 = zoneylay + izone_index*nrooms;
  tl0 = zonetl + izone_index*nrooms;
  tu0 = zonetu + izone_index*nrooms;
  hvent0 = zonehvents + izone_index*nzhvents;
#ifdef pp_ZONEVENT
  zoneslab_n0  = zoneslab_n  + izone_index*nzhvents;
  zoneslab_T0  = zoneslab_T  + izone_index*MAXSLABS*nzhvents;
  zoneslab_F0  = zoneslab_F  + izone_index*MAXSLABS*nzhvents;
  zoneslab_YB0 = zoneslab_YB + izone_index*MAXSLABS*nzhvents;
  zoneslab_YT0 = zoneslab_YT + izone_index*MAXSLABS*nzhvents;
#endif
  if(zoneodl!=NULL)odl0 = zoneodl + izone_index*nrooms;
  if(zoneodu!=NULL)odu0 = zoneodu + izone_index*nrooms;
#ifdef pp_ZONEVENT
  for(ivent=0;ivent<nzhvents;ivent++){
    zvent *zventi;
    int islab;

    zventi = zventinfo + ivent;
    zventi->area_fraction=CLAMP(hvent0[ivent]/zventi->area,0.0,1.0);

    zventi->nslab = zoneslab_n0[ivent];
    for(islab = 0; islab<zventi->nslab; islab++){
      int iislab;

      iislab = MAXSLABS*ivent+islab;
      zventi->slab_bot[islab] = zoneslab_YB0[iislab];
      zventi->slab_top[islab] = zoneslab_YT0[iislab];
      zventi->slab_vel[islab] =  zoneslab_F0[iislab];
      zventi->slab_temp[islab] = zoneslab_T0[iislab];
    }
  }
#else
  for(ivent=0;ivent<nzhvents;ivent++){
    zvent *zventi;

    zventi = zventinfo + ivent;
    zventi->area_fraction=CLAMP(hvent0[ivent]/zventi->area,0.0,1.0);
  }
#endif
  for(iroom=0;iroom<nrooms;iroom++){
    roomi = roominfo + iroom;
    roomi->pfloor=pr0[iroom];
    roomi->ylay=ylay0[iroom];
    roomi->tl=C2K(tl0[iroom]);
    roomi->tu=C2K(tu0[iroom]);
    roomi->itl=getZoneColor(tl0[iroom],zonemin,zonemax,nrgb_full);
    roomi->itu=getZoneColor(tu0[iroom],zonemin,zonemax,nrgb_full);
    roomi->rho_L=(pref+pr0[iroom])/R/roomi->tl;
    roomi->rho_U=(pref+pr0[iroom])/R/roomi->tu;
    if(zoneodl!=NULL)roomi->od_L=1.0/MAX(odl0[iroom],0.0001);
    if(zoneodu!=NULL)roomi->od_U=1.0/MAX(odu0[iroom],0.0001);
  }
  roomi=roominfo+nrooms;
  roomi->pfloor=0.0;
  roomi->ylay=99999.0;
  roomi->tl=tamb;
  roomi->tu=tamb;
  roomi->itl=getZoneColor(K2C(tamb),zonemin,zonemax,nrgb_full);
  roomi->itu=getZoneColor(K2C(tamb),zonemin,zonemax,nrgb_full);
  roomi->rho_L=(pref+pamb)/R/roomi->tl;
  roomi->rho_U=(pref+pamb)/R/roomi->tu;
  roomi->z0=0.0;
  roomi->z1=100000.0;
}

/* ------------------ get_p ------------------------ */

float get_p(float y, float pfloor, float ylay, float rho_L, float rho_U){
  float g=9.80;
  float p;

  if(y<ylay){
    p = pfloor - SCALE2FDS(rho_L*g*y);
  }
  else{
    p = pfloor - SCALE2FDS(rho_L*g*ylay);
    p -= SCALE2FDS(rho_U*g*(y-ylay));
  }
  return p;
}

/* ------------------ get_dpT ------------------------ */

void get_dpT(float *yy, int n, roomdata *r1, roomdata *r2, float *delp, float *dpmin, float *dpmax, int *iT){
  float p1, p2;
  int itslab;
  int fsign;
  int i;
  float y;

  for(i=0;i<n;i++){

    y=yy[i];

    if(y<r1->z0||y<r2->z0||y>r1->z1||y>r2->z1){
      delp[i]=0.0;
      iT[i]=r1->itl;
    }

    p1=get_p(y,              r1->pfloor,r1->ylay,r1->rho_L,r1->rho_U);
    p2=get_p(y+r1->z0-r2->z0,r2->pfloor,r2->ylay,r2->rho_L,r2->rho_U);

    if(p1>p2){
      fsign=1.0;
      if(y>r1->ylay){
        itslab=r1->itu;
      }
      else{
        itslab=r1->itl;
      }
    }
    else{
      fsign=-1.0;
      if(y>r2->ylay){
        itslab=r2->itu;
      }
      else{
        itslab=r2->itl;
      }
    }
    delp[i]=fsign*sqrt(ABS(p1-p2));
    iT[i]=itslab;
  }
  *dpmin=delp[0];
  *dpmax=delp[0];
  for(i=1;i<n;i++){
    if(delp[i]<*dpmin)*dpmin=delp[i];
    if(delp[i]>*dpmax)*dpmax=delp[i];
  }
  return;
}

/* ------------------ drawroomgeom ------------------------ */

void drawroomgeom(void){
  float xroom0, yroom0, zroom0, xroom, yroom, zroom;
  float x1,x2,yy1,yy2,z1,z2;
  int i;
  int idir;
  float yy,zz;
  
  fill_zonedata(izone);

/* draw the frame */

  antialias(ON);
  glBegin(GL_LINES);

  for(i=0;i<nrooms;i++){
    roomdata *roomi;

    if(zone_highlight==1&&zone_highlight_room==i){
      glEnd();
      glLineWidth(5.0*linewidth);
      glBegin(GL_LINES);
      glColor3f(1.0,0.0,0.0);
    }
    else{
      glEnd();
      glLineWidth(linewidth);
      glBegin(GL_LINES);
      glColor4fv(foregroundcolor);
    }

    roomi = roominfo + i;
    xroom0 = roomi->x0;
    yroom0 = roomi->y0;
    zroom0 = roomi->z0;
    xroom = roomi->x1;
    yroom = roomi->y1;
    zroom = roomi->z1;

 
    glVertex3f(xroom0,yroom0,zroom);
    glVertex3f(xroom,yroom0,zroom);

    glVertex3f(xroom,yroom0,zroom);
    glVertex3f(xroom,yroom,zroom);

    glVertex3f(xroom,yroom,zroom);
    glVertex3f(xroom0,yroom,zroom);

    glVertex3f(xroom0,yroom,zroom);
    glVertex3f(xroom0,yroom0,zroom);

    glVertex3f(xroom0,yroom0,zroom0);
    glVertex3f(xroom,yroom0,zroom0);

    glVertex3f(xroom,yroom0,zroom0);
    glVertex3f(xroom,yroom,zroom0);

    glVertex3f(xroom,yroom,zroom0);
    glVertex3f(xroom0,yroom,zroom0);

    glVertex3f(xroom0,yroom,zroom0);
    glVertex3f(xroom0,yroom0,zroom0);

    glVertex3f(xroom0,yroom0,zroom0);
    glVertex3f(xroom0,yroom0,zroom);

    glVertex3f(xroom,yroom0,zroom0);
    glVertex3f(xroom,yroom0,zroom);

    glVertex3f(xroom,yroom,zroom0);
    glVertex3f(xroom,yroom,zroom);

    glVertex3f(xroom0,yroom,zroom0);
    glVertex3f(xroom0,yroom,zroom);
  }
  glEnd();
  antialias(OFF);

  if(visVents==1){
    glLineWidth(ventlinewidth);
    for(i=0;i<nzvents;i++){
      zvent *zvi;

      zvi = zventinfo + i;

      glColor4fv(zvi->color);
      idir=zvi->dir;
      x1=zvi->x1;
      x2=zvi->x2;
      z1=zvi->z1;
      z2=zvi->z2;
      yy=zvi->yy;
      glBegin(GL_LINE_LOOP);
      switch(idir){
      case 1:
      case 3:
        glVertex3f(x1,yy,z1);
        glVertex3f(x2,yy,z1);
        glVertex3f(x2,yy,z2);
        glVertex3f(x1,yy,z2);
        glVertex3f(x1,yy,z1);
        break;
      case 2:
      case 4:
        glVertex3f(yy,x1,z1);
        glVertex3f(yy,x2,z1);
        glVertex3f(yy,x2,z2);
        glVertex3f(yy,x1,z2);
        glVertex3f(yy,x1,z1);
        break;
      case 5:
      case 6:
        yy1=zvi->y1;
        yy2=zvi->y2;
        zz=zvi->zz;
        glVertex3f(x1,yy1,zz);
        glVertex3f(x2,yy1,zz);
        glVertex3f(x2,yy2,zz);
        glVertex3f(x1,yy2,zz);
        glVertex3f(x1,yy1,zz);
        break;
      default:
        ASSERT(FFALSE);
        break;
      }
      glEnd();
    }
  }
}

/* ------------------ getzoneventbounds ------------------------ */

void getzoneventbounds(void){
  int i;

  for(i=0;i<nzvents;i++){
    zvent *zvi;

    zvi = zventinfo + i;
    zvi->g_dpmax=-1000000000.0;
    zvi->g_dpmin=1000000000.0;
  }
  for(izone=0;izone<nzone_times;izone++){
    fill_zonedata(izone);
    for(i=0;i<nzvents;i++){
      int j;
      zvent *zvi;
      float yelev[20];

      zvi = zventinfo + i;
      if(zvi->vent_orien==VFLOW_VENT||zvi->vent_orien==HVAC_VENT)continue;
      for(j=0;j<20;j++){
        yelev[j]=(zvi->z1*(19-j)+zvi->z2*j)/19.0;
      }
      get_dpT(yelev, 20, zvi->room1, zvi->room2, zvi->vdata, &zvi->dpmin, &zvi->dpmax, zvi->itempdata);
      if(zvi->dpmin<zvi->g_dpmin)zvi->g_dpmin=zvi->dpmin;
      if(zvi->dpmax>zvi->g_dpmax)zvi->g_dpmax=zvi->dpmax;
    }
  }
  zone_maxventflow=0.0;
  for(i=0;i<nzvents;i++){
    zvent *zvi;

    zvi = zventinfo + i;
    if(zvi->vent_orien==VFLOW_VENT||zvi->vent_orien==HVAC_VENT)continue;
    if(ABS(zvi->g_dpmin)>zone_maxventflow)zone_maxventflow=ABS(zvi->g_dpmin);
    if(ABS(zvi->g_dpmax)>zone_maxventflow)zone_maxventflow=ABS(zvi->g_dpmax);
  }
}

/* ------------------ drawventdataORIG ------------------------ */

#ifdef pp_ZONEVENT
void drawventdataORIG(void){
#else
void drawventdata(void){
#endif
  float factor;
  int i;
  int idir;
  float x1, yy;

  if(visVentFlow==0)return;

  if(cullfaces==1)glDisable(GL_CULL_FACE);

  for(i=0;i<nzvents;i++){
    int j;
    zvent *zvi;
    float yelev[20];

    zvi = zventinfo + i;
    if(zvi->vent_orien==VFLOW_VENT||zvi->vent_orien==HVAC_VENT)continue;
    for(j=0;j<20;j++){
      yelev[j]=(zvi->z1*(19-j)+zvi->z2*j)/19.0;
    }
    get_dpT(yelev, 20, zvi->room1, zvi->room2, zvi->vdata, &zvi->dpmin, &zvi->dpmax, zvi->itempdata);
  }
  factor = 0.1*zone_ventfactor/zone_maxventflow;
  for(i=0;i<nzvents;i++){
    zvent *zvi;
    int j;
    float yelev[20];
    float *vcolor1,*vcolor2;

    zvi = zventinfo + i;

    if(zvi->vent_orien==VFLOW_VENT||zvi->vent_orien==HVAC_VENT)continue;
    for(j=0;j<20;j++){
      yelev[j]=(zvi->z1*(19-j)+zvi->z2*j)/19.0;
    }
    idir=zvi->dir;
    x1=(zvi->x1+zvi->x2)/2.0;
    yy=zvi->yy;
    glBegin(GL_QUADS);
    for(j=0;j<19;j++){
      float dy1,dy2;

      dy1 = factor*zvi->area_fraction*zvi->vdata[j];
      dy2 = factor*zvi->area_fraction*zvi->vdata[j+1];
      if(idir==1||idir==4){
        dy1=-dy1;
        dy2=-dy2;
      }
      vcolor1=rgb_full[zvi->itempdata[j]];
      vcolor2=rgb_full[zvi->itempdata[j+1]];
      vcolor2=vcolor1;
      switch(idir){
      case 4:
      case 2:
        if(dy1*dy2>=0.0){
          glColor3fv(vcolor1);
          glVertex3f(yy,    x1,yelev[j]);
          glVertex3f(yy+dy1,x1,yelev[j]);
     
          glColor3fv(vcolor2);
          glVertex3f(yy+dy2,x1,yelev[j+1]);
          glVertex3f(yy,    x1,yelev[j+1]);
        }
        else{
          float dyy;

          dyy =  yelev[j] - dy1*(yelev[j+1]-yelev[j])/(dy2-dy1);
          glColor3fv(vcolor1);
          glVertex3f(yy,    x1,yelev[j]);
          glVertex3f(yy+dy1,x1,yelev[j]);
          glVertex3f(yy,    x1,dyy);
          glVertex3f(yy,x1,dyy);

          glColor3fv(vcolor2);
          glVertex3f(yy,    x1,dyy);
          glVertex3f(yy,    x1,dyy);
          glVertex3f(yy+dy2,x1,yelev[j+1]);
          glVertex3f(yy,    x1,yelev[j+1]);
        }
        break;
      case 3:
      case 1:
        if(dy1*dy2>=0.0){
          glColor3fv(vcolor1);
          glVertex3f(x1,yy,    yelev[j]);
          glVertex3f(x1,yy+dy1,yelev[j]);

          glColor3fv(vcolor2);
          glVertex3f(x1,yy+dy2,yelev[j+1]);
          glVertex3f(x1,yy,    yelev[j+1]);
        }
        else{
          float dyy;

          dyy =  yelev[j] - dy1*(yelev[j+1]-yelev[j])/(dy2-dy1);
          glColor3fv(vcolor1);
          glVertex3f(x1,yy,    yelev[j]);
          glVertex3f(x1,yy+dy1,yelev[j]);
          glVertex3f(x1,yy,dyy);
          glVertex3f(x1,yy,dyy);

          glColor3fv(vcolor2);
          glVertex3f(x1,yy,dyy);
          glVertex3f(x1,yy,dyy);
          glVertex3f(x1,yy+dy2,yelev[j+1]);
          glVertex3f(x1,yy,    yelev[j+1]);
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
      }
    }
    glEnd();
  }
  if(cullfaces==1)glEnable(GL_CULL_FACE);

}

#ifdef pp_ZONEVENT
/* ------------------ drawventslabdata ------------------------ */

void drawventslabdata(void){
  float factor;
  int i;
  int idir;
  float x1, yy, dyy;

  if(visVentFlow==0)return;

  if(cullfaces==1)glDisable(GL_CULL_FACE);

  for(i = 0; i<nzvents; i++){
    zvent *zvi;
    int j;
    float yelev[20];
    float *vcolor1, *vcolor2;
    float *slab_vel;
    int islab;

    zvi = zventinfo+i;

    if(zvi->vent_orien==HVAC_VENT)continue;
    idir = zvi->dir;
    x1 = (zvi->x1+zvi->x2)/2.0;
    yy = zvi->yy;

    slab_vel = zvi->slab_vel;
    glBegin(GL_QUADS);
    for(islab = 0; islab<zvi->nslab;islab++){
      float slab_bot, slab_top, tslab, *tcolor;
      int itslab;

      slab_bot = NORMALIZE_Z(zvi->slab_bot[islab]);
      slab_top = NORMALIZE_Z(zvi->slab_top[islab]);
      tslab = zvi->slab_temp[islab];
      itslab = getZoneColor(K2C(tslab), zonemin, zonemax, nrgb_full);
      tcolor = rgb_full[itslab];
      glColor3fv(tcolor);

      dyy = -0.1*(slab_vel[islab]/maxslabflow);
      switch(idir){
      case 5:
      case 6:
        break;
      case 4:
      case 2:
        glVertex3f(yy,     x1, slab_bot);
        glVertex3f(yy+dyy, x1, slab_bot);

        glVertex3f(yy+dyy, x1, slab_top);
        glVertex3f(yy,     x1, slab_top);
        break;
      case 3:
      case 1:
        glVertex3f(x1,     yy, slab_bot);
        glVertex3f(x1, yy+dyy, slab_bot);

        glVertex3f(x1, yy+dyy, slab_top);
        glVertex3f(x1,     yy, slab_top);
        break;
      default:
        ASSERT(FFALSE);
        break;
      }
    }
    glEnd();
  }
  if(cullfaces==1)glEnable(GL_CULL_FACE);
}

/* ------------------ drawventdata ------------------------ */

void drawventdata(void){
  if(have_hventslab_flow==1){
    drawventslabdata();
  }
  else{
    drawventdataORIG();
  }
}
#endif

/* ------------------ getzonethick ------------------------ */

float getzonethick(int dir, roomdata *roomi, float xyz[3]){
  float dx, dy, dz, L;
  float alpha, alpha_min, alpha_ylay;
  float x0, x1, yy0, yy1, z0, z1;
  float factor_L, factor_U, factor;
  float ylay;
  float thick;
  float odl, odu;

  x0 = roomi->x0;
  x1 = roomi->x1;
  yy0 = roomi->y0;
  yy1 = roomi->y1;
  z0 = roomi->z0;
  z1 = roomi->z1;
  odl = roomi->od_L;
  odu = roomi->od_U;
  alpha_min = 100000.0;
  ylay = roomi->z0 + roomi->ylay;

  dx = xyz[0] - xyzeyeorig[0];
  dy = xyz[1] - xyzeyeorig[1];
  dz = xyz[2] - xyzeyeorig[2];
  L = sqrt(dx*dx+dy*dy+dz*dz);

  alpha_ylay = (ylay - xyz[2])/dz;
  if(roomi->zoneinside==0){
  if(dir!=-1){
    alpha =  (x0 - xyz[0])/dx;
    if(alpha>0.0&&alpha<alpha_min){
      alpha_min=alpha;
    }
  }
  if(dir!=1){
    alpha =  (x1 - xyz[0])/dx;
    if(alpha>0.0&&alpha<alpha_min){
      alpha_min=alpha;
    }
  }
  if(dir!=-2){
    alpha =  (yy0 - xyz[1])/dy;
    if(alpha>0.0&&alpha<alpha_min){
      alpha_min=alpha;
    }
  }
  if(dir!=2){
    alpha =  (yy1 - xyz[1])/dy;
    if(alpha>0.0&&alpha<alpha_min){
      alpha_min=alpha;
    }
  }
  if(dir!=-3){
    alpha =  (z0 - xyz[2])/dz;
    if(alpha>0.0&&alpha<alpha_min){
      alpha_min=alpha;
    }
  }
  if(dir!=3){
    alpha =  (z1 - xyz[2])/dz;
    if(alpha>0.0&&alpha<alpha_min){
      alpha_min=alpha;
    }
  }
  if(xyzeyeorig[2]>ylay&&xyz[2]>ylay){
    if(alpha_ylay>0.0&&alpha_ylay<alpha_min){
      factor_U=alpha_ylay/odu;
      factor_L=(alpha_min-alpha_ylay)/odl;
    }
    else{
      factor_U=alpha_min/odu;
      factor_L=0.0;
    }
  }
  if(xyzeyeorig[2]>ylay&&xyz[2]<=ylay){
    factor_U=0.0;
    factor_L=alpha_min/odl;
  }
  if(xyzeyeorig[2]<=ylay&&xyz[2]>ylay){
    factor_U=alpha_min/odu;
    factor_L=0.0;
  }
  if(xyzeyeorig[2]<=ylay&&xyz[2]<=ylay){
    if(alpha_ylay>0.0&&alpha_ylay<alpha_min){
      factor_U=(alpha_min-alpha_ylay)/odu;
      factor_L=alpha_ylay/odl;
    }
    else{
      factor_U=0.0;
      factor_L=alpha_min/odl;
    }
  }
  }
  else{
    if(xyzeyeorig[2]>ylay&&xyz[2]>ylay){
      factor_U=1.0/odu;
      factor_L=0.0;
    }
    if(xyzeyeorig[2]>ylay&&xyz[2]<=ylay){
      factor_U=(1.0+alpha_ylay)/odu;
      factor_L=-alpha_ylay/odl;
    }
    if(xyzeyeorig[2]<=ylay&&xyz[2]>ylay){
      factor_U=-alpha_ylay/odu;
      factor_L=(1.0+alpha_ylay)/odl;
    }
    if(xyzeyeorig[2]<=ylay&&xyz[2]<=ylay){
      factor_U=0.0;
      factor_L=1.0/odl;
    }
  }

  factor = SCALE2FDS((factor_U+factor_L)*L);
  thick = 1.0-exp(-factor);
  return thick;
}

#ifdef pp_GPU
/* ------------------ drawzonesmokeGPU ------------------------ */

void drawzonesmokeGPU(roomdata *roomi){
#define NROWS_GPU 2
#define NCOLS_GPU 2
  int iwall;
  float dx, dy, dz;
  
  glUniform3f(GPUzone_eyepos,xyzeyeorig[0],xyzeyeorig[1],xyzeyeorig[2]);
  glUniform1i(GPUzone_zoneinside,roomi->zoneinside);
  glUniform1f(GPUzone_xyzmaxdiff,xyzmaxdiff);
  glUniform3f(GPUzone_boxmin,roomi->x0,roomi->y0,roomi->z0);
  glUniform3f(GPUzone_boxmax,roomi->x1,roomi->y1,roomi->z1);
  glUniform1f(GPUzone_zlay,roomi->z0+roomi->ylay);
  glUniform1f(GPUzone_odl,roomi->od_L);
  glUniform1f(GPUzone_odu,roomi->od_U);

  for(iwall=-3;iwall<=3;iwall++){
    int i,j;
    float x1, x2, yy1, yy2, z1, z2;

    if(iwall==0||roomi->drawsides[iwall+3]==0)continue;

    glUniform1i(GPUzone_zonedir,iwall);
    glBegin(GL_TRIANGLES);

    switch(iwall){
      case 1:
      case -1:
        dy = roomi->dy/(NCOLS_GPU-1);
        dz = roomi->dz/(NROWS_GPU-1);
        if(iwall<0){
          x1 = roomi->x0;
        }
        else{
          x1=roomi->x1;
        }
        for(i=0;i<NCOLS_GPU-1;i++){
          yy1 = roomi->y0 + i*dy;
          yy2 = yy1 + dy;
          for(j=0;j<NROWS_GPU-1;j++){
            z1 = roomi->z0 + j*dz;
            z2 = z1 + dz;

            glVertex3f(x1,yy1,z1);
            glVertex3f(x1,yy2,z1);
            glVertex3f(x1,yy2,z2);

            glVertex3f(x1,yy1,z1);
            glVertex3f(x1,yy2,z2);
            glVertex3f(x1,yy1,z2);
          }
        }
        break;
      case 2:
      case -2:
        dx = roomi->dx/(NCOLS_GPU-1);
        dz = roomi->dz/(NROWS_GPU-1);
        if(iwall<0){
          yy1=roomi->y0;
        }
        else{
          yy1=roomi->y1;
        }
        for(i=0;i<NCOLS_GPU-1;i++){
          x1 = roomi->x0 + i*dx;
          x2 = x1 + dx;
          for(j=0;j<NROWS_GPU-1;j++){
            z1 = roomi->z0 + j*dz;
            z2 = z1 + dz;

            if(roomi->zoneinside==0){
              glVertex3f(x1,yy1,z1);
              glVertex3f(x2,yy1,z1);
              glVertex3f(x2,yy1,z2);

              glVertex3f(x1,yy1,z1);
              glVertex3f(x2,yy1,z2);
              glVertex3f(x1,yy1,z2);
            }
            else{
              glVertex3f(x1,yy1,z1);
              glVertex3f(x2,yy1,z2);
              glVertex3f(x2,yy1,z1);

              glVertex3f(x1,yy1,z1);
              glVertex3f(x1,yy1,z2);
              glVertex3f(x2,yy1,z2);
            }
          }
        }
        break;
      case 3:
      case -3:
        dx = roomi->dx/(NCOLS_GPU-1);
        dy = roomi->dy/(NROWS_GPU-1);
        if(iwall<0){
          z1=roomi->z0;
        }
        else{
          z1=roomi->z1;
        }
        for(i=0;i<NCOLS_GPU-1;i++){
          x1 = roomi->x0 + i*dx;
          x2 = x1 + dx;
          for(j=0;j<NROWS_GPU-1;j++){
            yy1 = roomi->y0 + j*dy;
            yy2 = yy1 + dy;

            if(roomi->zoneinside==0){
              glVertex3f(x1,yy1,z1);
              glVertex3f(x2,yy1,z1);
              glVertex3f(x2,yy2,z1);

              glVertex3f(x1,yy1,z1);
              glVertex3f(x2,yy2,z1);
              glVertex3f(x1,yy2,z1);
            }
            else{
              glVertex3f(x1,yy1,z1);
              glVertex3f(x2,yy2,z1);
              glVertex3f(x2,yy1,z1);

              glVertex3f(x1,yy1,z1);
              glVertex3f(x1,yy2,z1);
              glVertex3f(x2,yy2,z1);
            }
          }
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
    glEnd();
  }

}

#endif
/* ------------------ drawzonesmoke ------------------------ */

void drawzonesmoke(roomdata *roomi){
#define NROWS 100
#define NCOLS 100
  float vxyz[4][NROWS][NCOLS];
  int iwall;
  float xyz[3];
  float dx, dy, dz;
  
  for(iwall=-3;iwall<=3;iwall++){
    int i,j;
  
    if(iwall==0)continue;
    if(roomi->drawsides[iwall+3]==0)continue;

    switch(iwall){
      case 1:
      case -1:
        dy = roomi->dy/(NCOLS-1);
        dz = roomi->dz/(NROWS-1);
        if(iwall<0){
          xyz[0]=roomi->x0;
        }
        else{
          xyz[0]=roomi->x1;
        }
        for(i=0;i<NCOLS;i++){
          xyz[1] = roomi->y0 + i*dy;
          for(j=0;j<NROWS;j++){
            xyz[2] = roomi->z0 + j*dz;
            vxyz[0][i][j]=xyz[0];
            vxyz[1][i][j]=xyz[1];
            vxyz[2][i][j]=xyz[2];
            vxyz[3][i][j]=getzonethick(iwall,roomi,xyz);
          }
        }
        break;
      case 2:
      case -2:
        dx = roomi->dx/(NCOLS-1);
        dz = roomi->dz/(NROWS-1);
        if(iwall<0){
          xyz[1]=roomi->y0;
        }
        else{
          xyz[1]=roomi->y1;
        }
        for(i=0;i<NCOLS;i++){
          xyz[0] = roomi->x0 + i*dx;
          for(j=0;j<NROWS;j++){
            xyz[2] = roomi->z0 + j*dz;
            vxyz[0][i][j]=xyz[0];
            vxyz[1][i][j]=xyz[1];
            vxyz[2][i][j]=xyz[2];
            vxyz[3][i][j]=getzonethick(iwall,roomi,xyz);
          }
        }
        break;
      case 3:
      case -3:
        dx = roomi->dx/(NCOLS-1);
        dy = roomi->dy/(NROWS-1);
        if(iwall<0){
          xyz[2]=roomi->z0;
        }
        else{
          xyz[2]=roomi->z1;
        }
        for(i=0;i<NCOLS;i++){
          xyz[0] = roomi->x0 + i*dx;
          for(j=0;j<NROWS;j++){
            xyz[1] = roomi->y0 + j*dy;
            vxyz[0][i][j]=xyz[0];
            vxyz[1][i][j]=xyz[1];
            vxyz[2][i][j]=xyz[2];
            vxyz[3][i][j]=getzonethick(iwall,roomi,xyz);
          }
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
    }

    glBegin(GL_TRIANGLES);
    for(i=0;i<NCOLS-1;i++){
      for(j=0;j<NROWS-1;j++){
        float x11, x12, x22, x21;
        float y11, y12, y22, y21;
        float z11, z12, z22, z21;
        float a11, a12, a22, a21;
        float grey=0.0;

        x11 = vxyz[0][i][j];
        x12 = vxyz[0][i][j+1];
        x21 = vxyz[0][i+1][j];
        x22 = vxyz[0][i+1][j+1];
        y11 = vxyz[1][i][j];
        y12 = vxyz[1][i][j+1];
        y21 = vxyz[1][i+1][j];
        y22 = vxyz[1][i+1][j+1];
        z11 = vxyz[2][i][j];
        z12 = vxyz[2][i][j+1];
        z21 = vxyz[2][i+1][j];
        z22 = vxyz[2][i+1][j+1];
        a11 = vxyz[3][i][j];
        a12 = vxyz[3][i][j+1];
        a21 = vxyz[3][i+1][j];
        a22 = vxyz[3][i+1][j+1];

        glColor4f(grey,grey,grey,a11);
        glVertex3f(x11,y11,z11);
        glColor4f(grey,grey,grey,a12);
        glVertex3f(x12,y12,z12);
        glColor4f(grey,grey,grey,a22);
        glVertex3f(x22,y22,z22);

        glColor4f(grey,grey,grey,a11);
        glVertex3f(x11,y11,z11);
        glColor4f(grey,grey,grey,a22);
        glVertex3f(x22,y22,z22);
        glColor4f(grey,grey,grey,a21);
        glVertex3f(x21,y21,z21);
      }
    }
    glEnd();
  }
}

/* ------------------ drawfiredata ------------------------ */

void drawfiredata(void){
  int i;
  float *zoneqfirebase, *zonefheightbase, *zonefdiambase, *zonefbasebase;

  if(zone_times[0]>global_times[itimes])return;
  if(cullfaces==1)glDisable(GL_CULL_FACE);

  zoneqfirebase = zoneqfire + izone*nfires;
  zonefheightbase = zonefheight + izone*nfires;
  zonefdiambase = zonefdiam + izone*nfires;
  zonefbasebase = zonefbase + izone*nfires;

  if(viszonefire==1){
    for(i=0;i<nfires;i++){
      float qdot;
      float diameter, flameheight, maxheight;

      qdot = zoneqfirebase[i]/1000.0f;
      if(zonecsv==1){
        if(qdot>0.0f){
          firedata *firei;
          roomdata *roomi;
          float deltaz;
          mesh *meshi;

          // radius/plumeheight = .268 = atan(15 degrees)
          firei = fireinfo + i;
          roomi = roominfo + firei->roomnumber-1;
          meshi = meshinfo + firei->roomnumber-1;
          diameter = SCALE2SMV(zonefdiambase[i]);
          deltaz = SCALE2SMV(zonefbasebase[i]);
          maxheight=roomi->z1-roomi->z0-deltaz;
          flameheight = SCALE2SMV(zonefheightbase[i]);
          setClipPlanes(&(meshi->box_clipinfo),CLIP_ON);
          glPushMatrix();
          glTranslatef(firei->absx,firei->absy,roomi->z0+deltaz);
          DrawFirePlume(diameter,flameheight,maxheight);
          glPopMatrix();
          setClipPlanes(&(meshi->box_clipinfo),CLIP_OFF);
        }
      }
      else{
        if(qdot>0.0f){
          firedata *firei;
          roomdata *roomi;
          mesh *meshi;

          // radius/plumeheight = .268 = atan(15 degrees)
          firei = fireinfo + i;
          roomi = roominfo + firei->roomnumber-1;
          meshi = meshinfo + firei->roomnumber-1;
          maxheight=roomi->z1-firei->absz;
          flameheight = SCALE2SMV((0.23f*pow((double)qdot,(double)0.4)/(1.0f+2.0f*0.268f)));
          diameter = 2.0*flameheight*0.268f;
          setClipPlanes(&(meshi->box_clipinfo),CLIP_ON);
          glPushMatrix();
          glTranslatef(firei->absx,firei->absy,firei->absz);
          DrawFirePlume(diameter,flameheight,maxheight);
          glPopMatrix();
          setClipPlanes(&(meshi->box_clipinfo),CLIP_OFF);
        }
      }
    }
  }
  if(use_transparency_data==1)transparentoff();
  if(cullfaces==1)glEnable(GL_CULL_FACE);
}

/* ------------------ drawroomdata ------------------------ */

void drawroomdata(void){
  float xroom0, yroom0, zroom0, xroom, yroom, zroom;
  float *zoneylaybase,dy;
  unsigned char *hazardcolorbase, *zonecolorbase;
  float ylay;
  float *colorv;
  unsigned char color;
  unsigned char *izonetubase;
  int i;

  if(zone_times[0]>global_times[itimes])return;

  if(cullfaces==1)glDisable(GL_CULL_FACE);
  if(use_transparency_data==1)transparenton();

  izonetubase = izonetu + izone*nrooms;
  hazardcolorbase = hazardcolor + izone*nrooms;
  zoneylaybase = zoneylay + izone*nrooms;

  if(zonecolortype==ZONEHAZARD_COLOR){
    zonecolorbase=hazardcolorbase;
  }
  else{
    zonecolorbase=izonetubase;
  }

#ifdef pp_GPU
   if(usegpu==1){
     LoadZoneSmokeShaders();
   }
#endif
  for(i=0;i<nrooms;i++){
    roomdata *roomi;

    roomi = roominfo + i;

    ylay = *(zoneylaybase+i);
    color = *(zonecolorbase+i);
    if(zonecolortype==ZONEHAZARD_COLOR){
      colorv = rgbhazard[color];
    }
    else{
      colorv = rgb_full[color];
    }
    xroom0 = roomi->x0;
    yroom0 = roomi->y0;
    xroom = roomi->x1;
    yroom = roomi->y1;
    zroom0 = roomi->z0;
    zroom = roomi->z1;
    dy = roomi->dy/2.;

    if(zonecolortype==ZONESMOKE_COLOR&&visSZone==1){
#ifdef pp_GPU
      if(usegpu==1){
        drawzonesmokeGPU(roomi);
      }
      else{
        drawzonesmoke(roomi);
      }
#else
      drawzonesmoke(roomi);
#endif
    }
    else{
      if(visHZone==1){
        glColor4fv(colorv);
        glBegin(GL_QUADS);
        glVertex3f(xroom0,yroom0,ylay+zroom0);
        glVertex3f(xroom,yroom0,ylay+zroom0);
        glVertex3f(xroom,yroom,ylay+zroom0);
        glVertex3f(xroom0,yroom,ylay+zroom0);
        glEnd();
      }
      if(visVZone==1){
        glColor4fv(colorv);
        glBegin(GL_QUADS);
        glVertex3f(xroom0,yroom0+dy,ylay+zroom0);
        glVertex3f(xroom,yroom0+dy,ylay+zroom0);
        glVertex3f(xroom, yroom0+dy,zroom);
        glVertex3f(xroom0,yroom0+dy,zroom);
        glEnd();
      }
    }
  }
#ifdef pp_GPU
  if(usegpu==1){
    UnLoadShaders();
  }
#endif

  if(use_transparency_data==1)transparentoff();
  if(cullfaces==1)glEnable(GL_CULL_FACE);
}

/* ------------------ DrawFirePlume ------------------------ */

void DrawFirePlume(float diameter, float height, float maxheight){
  unsigned char firecolor[3]={255,128,0};

  if(height<=maxheight){
    float dlower, dupper;

    dlower=diameter;
    dupper=0.0;
    drawtrunccone(dlower,dupper,height,firecolor);
  }
  else{
    float dh;
    float dlower1, dupper1;
    float dlower2, dupper2;

    dh = height-maxheight;
    dlower1=diameter;
    dupper1=diameter*dh/height;
    drawtrunccone(dlower1,dupper1,maxheight,firecolor);

    glPushMatrix();
    glTranslatef(0.0,0.0,maxheight-dupper1/2.0);
    dlower2=0.0;
    dupper2=2.0*dh;
    drawtrunccone(dlower2,dupper2,dupper1/2.0,firecolor);
    glPopMatrix();
  }

}

