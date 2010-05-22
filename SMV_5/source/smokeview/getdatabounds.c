// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include "MALLOC.h"
#include "flowfiles.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"

// svn revision character string
char getdatabounds_revision[]="$Revision$";

/* ------------------ adjustdatabounds ------------------------ */

void adjustdatabounds(const float *pdata, int local_skip, int ndata,
                   int setpmin, float *pmin, int setpmax, float *pmax)
{
    int nsmall, nbig, *buckets=NULL, n, level, total, alpha05;
    float dp, pmin2, pmax2;

    if(setpmin==PERCENTILE_MIN||setpmax==PERCENTILE_MAX){
      dp = (*pmax - *pmin)/NBUCKETS;
      nsmall=0;
      nbig=NBUCKETS;
      if(NewMemory((void **)&buckets,NBUCKETS*sizeof(int))==0){
        printf("*** Warning: Unable to allocate memory in getdatabounds\n");
        return;
      }

      for (n=0;n<NBUCKETS;n++){
        buckets[n]=0;
      }
      for (n=local_skip;n<ndata;n++){
        level=0;
        if(dp!=0.0f)level = (int)((pdata[n] - *pmin)/dp);
        if(level<0){
          level=0;
        }
        if(level>NBUCKETS-1){
          level=NBUCKETS-1;
        }
        buckets[level]++;
      }
      alpha05 = (int)(percentile_level*ndata);
      total = 0;
      for (n=0;n<NBUCKETS;n++){
        total += buckets[n];
        if(total>alpha05){
          nsmall=n;
          break;
        }
      }
      total = 0;
      for (n=NBUCKETS;n>0;n--){
        total += buckets[n-1];
        if(total>alpha05){
          nbig=n;
          break;
        }
      }
      pmin2 = *pmin + (nsmall-1)*dp;
      pmax2 = *pmin + (nbig+1)*dp;
      if(setpmin==PERCENTILE_MIN)*pmin = pmin2; 
      if(setpmax==PERCENTILE_MAX)*pmax = pmax2;
      FreeMemory(buckets);
    }
    if(axissmooth==1){
      smoothlabel(pmin,pmax,nrgb);
    }
}

/* ------------------ adjustpart5bounds ------------------------ */

void adjustpart5bounds(particle *parti){
  int i,j,k,m;
  part5data *datacopy;
  int alpha05;

  for(i=0;i<npart5prop;i++){
    part5prop *propi;
    int n;
    int *buckets;

    propi = part5propinfo + i;

    if(propi->set_global_bounds==1){
      propi->global_min =  1000000000.0;
      propi->global_max = -1000000000.0;
    }

    NewMemory((void **)&propi->buckets,NBUCKETS*sizeof(int));
    buckets = propi->buckets;
    for (n=0;n<NBUCKETS;n++){
      buckets[n]=0;
    }
  }

  // compute global min and max

  datacopy = parti->data5;
  for(i=0;i<parti->nframes;i++){
    for(j=0;j<parti->nclasses;j++){
      part5class *partclassi;
      float *rvals;

      partclassi=parti->partclassptr[j];
      rvals = datacopy->rvals;

      for(k=2;k<partclassi->ntypes;k++){
        part5prop *prop_id;
        float *valmin, *valmax;

        prop_id = get_part5prop(partclassi->labels[k].longlabel);
        if(prop_id==NULL)continue;

        valmin = &prop_id->global_min;
        valmax = &prop_id->global_max;

        if(prop_id->set_global_bounds==1&&(valmin!=NULL||valmax!=NULL)){
          for(m=0;m<datacopy->npoints;m++){
            float val;

            val=*rvals++;
            if(valmin!=NULL&&val<*valmin)*valmin=val;
            if(valmax!=NULL&&val>*valmax)*valmax=val;
          }
        }
      }
      datacopy++;
    }
  }

  // generate data histogram (buckets) in order to determine percentile min and max

  datacopy = parti->data5;
  for(i=0;i<parti->nframes;i++){
    for(j=0;j<parti->nclasses;j++){
      part5class *partclassi;
      float *rvals;

      partclassi = parti->partclassptr[j];
      rvals = datacopy->rvals;

      for(k=2;k<partclassi->ntypes;k++){
        part5prop *prop_id;
        float *valmin, *valmax, dg;
        int *buckets;

        prop_id = get_part5prop(partclassi->labels[k].longlabel);
        if(prop_id==NULL)continue;

        valmin = &prop_id->global_min;
        valmax = &prop_id->global_max;
        dg = (*valmax-*valmin)/(float)NBUCKETS;
        if(dg==0.0)dg=1.0;
        buckets=prop_id->buckets;

        for(m=0;m<datacopy->npoints;m++){
          float val;
          int ival;

          val=*rvals++;
          ival = (val-*valmin)/dg;
          if(ival<0)ival=0;
          if(ival>NBUCKETS-1)ival=NBUCKETS-1;
          buckets[ival]++;
        }
      }
      datacopy++;
    }
  }

  // calculate percentile min and max

  for(i=0;i<npart5prop;i++){
    part5prop *propi;
    int total;
    int *buckets;
    int nsmall, nbig;
    int n;
    float gmin, gmax, dg;

    propi = part5propinfo + i;
    buckets = propi->buckets;

    if(propi->set_global_bounds==1){
      total = 0;
      for(n=0;n<NBUCKETS;n++){
        total+=buckets[n];
      }
      alpha05 = (int)(percentile_level*total);

      total = 0;
      nsmall=0;
      nbig = NBUCKETS-1;
      for (n=0;n<NBUCKETS;n++){
        total += buckets[n];
        if(total>alpha05){
          nsmall=n;
          break;
        }
      }
      total = 0;
      for (n=NBUCKETS;n>0;n--){
        total += buckets[n-1];
        if(total>alpha05){
          nbig=n-1;
          break;
        }
      }
      gmin = propi->global_min;
      gmax = propi->global_max;
      dg = (gmax-gmin)/(float)NBUCKETS;
      propi->percentile_min = gmin + nsmall*dg;
      propi->percentile_max = gmin + nbig*dg;
    }
    if(propi->global_min<1000000000.0)propi->set_global_bounds=0;

    FREEMEMORY(propi->buckets);
  }
  for(i=0;i<npart5prop;i++){
    part5prop *propi;

    propi = part5propinfo + i;

    switch (propi->setvalmin){
    case PERCENTILE_MIN:
      propi->valmin=propi->percentile_min;
  //    propi->user_min=propi->valmin;
      break;
    case GLOBAL_MIN:
      propi->valmin=propi->global_min;
      break;
    case SET_MIN:
      propi->valmin=propi->user_min;
      break;
    }
    switch (propi->setvalmax){
    case PERCENTILE_MAX:
      propi->valmax=propi->percentile_max;
//      propi->user_max=propi->valmax;
      break;
    case GLOBAL_MAX:
      propi->valmax=propi->global_max;
      break;
    case SET_MAX:
      propi->valmax=propi->user_max;
      break;
    }
  }
#ifdef _DEBUG
  print_part5prop();
#endif
}

/* ------------------ adjustpartbounds ------------------------ */

void adjustpartbounds(const float *pdata, int particle_type, int droplet_type, const unsigned char *isprink, 
                      int local_skip, int ndataloop, int setpmin, float *pmin, int setpmax, float *pmax)
{
    int nsmall, nbig, *buckets=NULL, n, level, total, alpha05;
    float dp, pmin2, pmax2;
    int ndata;

    if(setpmin==PERCENTILE_MIN||setpmax==PERCENTILE_MAX){
      dp = (*pmax - *pmin)/NBUCKETS;
      nsmall=0;
      nbig=NBUCKETS;
      if(NewMemory((void **)&buckets,NBUCKETS*sizeof(int))==0){
        printf("*** Warning: Unable to allocate memory in getdatabounds\n");
        return;
      }

      for (n=0;n<NBUCKETS;n++){
        buckets[n]=0;
      }
      ndata=0;
      for (n=local_skip;n<ndataloop;n++){
        level=0;
        if(isprink[n]==1){
          if(droplet_type==0)continue;
        }
        else{
          if(particle_type==0)continue;
        }
        if(dp!=0.0f)level = (int)((pdata[n] - *pmin)/dp);
        if(level<0){
          level=0;
        }
        if(level>NBUCKETS-1){
          level=NBUCKETS-1;
        }
        ndata++;
        buckets[level]++;
      }
      alpha05 = (int)(percentile_level*ndata);
      total = 0;
      for (n=0;n<NBUCKETS;n++){
        total += buckets[n];
        if(total>alpha05){
          nsmall=n;
          break;
        }
      }
      total = 0;
      for (n=NBUCKETS;n>0;n--){
        total += buckets[n-1];
        if(total>alpha05){
          nbig=n;
          break;
        }
      }
      pmin2 = *pmin + (nsmall-1)*dp;
      pmax2 = *pmin + (nbig+1)*dp;
      if(setpmin==PERCENTILE_MIN)*pmin = pmin2; 
      if(setpmax==PERCENTILE_MAX)*pmax = pmax2;
      FreeMemory(buckets);
    }
    if(axissmooth==1){
      smoothlabel(pmin,pmax,nrgb);
    }
}

/* ------------------ adjustPlot3Dbounds ------------------------ */

void adjustPlot3Dbounds(int plot3dvar, int setpmin, float *pmin, int setpmax, float *pmax)
{
    int nsmall, nbig, *buckets=NULL, n, level, total, alpha05;
    float dp, pmin2, pmax2;
    plot3d *p;
    mesh *meshi;
    int i;
    char *iblank;
    float *q;
    int ndata=0;
    int ntotal;

    if(setpmin==PERCENTILE_MIN||setpmax==PERCENTILE_MAX){
      dp = (*pmax - *pmin)/NBUCKETS;
      nsmall=0;
      nbig=NBUCKETS;
      if(NewMemory((void **)&buckets,NBUCKETS*sizeof(int))==0){
        printf("*** Warning: Unable to allocate memory in getdatabounds\n");
        return;
      }

      for (n=0;n<NBUCKETS;n++){
        buckets[n]=0;
      }
      for(i=0;i<nplot3d_files;i++){
        p = plot3dinfo+i;
        if(p->loaded==0||p->display==0)continue;
        meshi=meshinfo+p->blocknumber;
        ntotal=(meshi->ibar+1)*(meshi->jbar+1)*(meshi->kbar+1);
        ndata += ntotal;
        iblank=meshi->c_iblank;
        q=meshi->qdata+plot3dvar*ntotal;
        for(n=0;n<ntotal;n++){
          if(iblank==NULL||*iblank++==1){
            level=0;
            if(dp!=0.0f)level = (int)((q[n] - *pmin)/dp);
            if(level<0){
              level=0;
            }
            if(level>NBUCKETS-1){
              level=NBUCKETS-1;
            }
            buckets[level]++;
          }
        }
      }
      alpha05 = (int)(percentile_level*ndata);
      total = 0;
      for (n=0;n<NBUCKETS;n++){
        total += buckets[n];
        if(total>alpha05){
          nsmall=n;break;
        }
      }
      total = 0;
      for (n=NBUCKETS;n>0;n--){
        total += buckets[n-1];
        if(total>alpha05){
          nbig=n;break;
        }
      }
      pmin2 = *pmin + (nsmall-1)*dp;
      pmax2 = *pmin + (nbig+1)*dp;
      if(setpmin==PERCENTILE_MIN)*pmin = pmin2; 
      if(setpmax==PERCENTILE_MAX)*pmax = pmax2;
      FreeMemory(buckets);
    }
    if(axissmooth==1&&setpmin!=SET_MIN&&setpmax!=SET_MAX){
      smoothlabel(pmin,pmax,nrgb);
    }

}

/* ------------------ getzonebounds ------------------------ */

void getzonebounds(const float *pdata, int ndata, 
                   int setpmin, float *pmin, int setpmax, float *pmax)
{
    int n;
    float pmin2, pmax2, val;

    if(setpmin==SET_MIN&&setpmax==SET_MAX)return;
    pmin2 = *pdata;
    pmax2 = pmin2;
    for (n=0;n<ndata;n++){
      val=*pdata++;
      if(val<pmin2)pmin2=val;
      if(val>pmax2)pmax2=val;
    }
    if(setpmin!=SET_MIN)*pmin = pmin2; 
    if(setpmax!=SET_MAX)*pmax = pmax2;
}

