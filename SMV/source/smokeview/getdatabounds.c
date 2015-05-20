#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "smokeviewvars.h"

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
        fprintf(stderr,"*** Error: Unable to allocate memory in getdatabounds\n");
        return;
      }

      for (n=0;n<NBUCKETS;n++){
        buckets[n]=0;
      }
      for (n=local_skip;n<ndata;n++){
        level=0;
        if(dp!=0.0f)level = CLAMP((int)((pdata[n] - *pmin)/dp),0,NBUCKETS-1);
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
    if(axislabels_smooth==1){
      smoothlabel(pmin,pmax,nrgb);
    }
}


/* ------------------ adjustpart5bounds ------------------------ */

void adjustpart5chops(partdata *parti){
  int i;
  
  for(i=0;i<npart5prop;i++){
    part5prop *propi;

    propi = part5propinfo + i;
    propi->imin=0;
    propi->imax=255;
    if(propi->setchopmin==1){
      float dval;

      dval = propi->valmax-propi->valmin;
      if(dval<=0.0)dval=1.0;
      propi->imin=CLAMP(255*(propi->chopmin-propi->valmin)/dval,0,255);
    }
    if(propi->setchopmax==1){
      float dval;

      dval = propi->valmax-propi->valmin;
      if(dval<=0.0)dval=1.0;
      propi->imax=CLAMP(255*(propi->chopmax-propi->valmin)/dval,0,255);
    }
  }
}

/* ------------------ adjustpart5bounds ------------------------ */

void adjustpart5bounds(partdata *parti){
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
  for(i=0;i<parti->ntimes;i++){
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
            if(valmin!=NULL)*valmin=MIN(*valmin,val);
            if(valmax!=NULL)*valmax=MAX(*valmax,val);
          }
        }
      }
      datacopy++;
    }
  }

  // generate data histogram (buckets) in order to determine percentile min and max

  datacopy = parti->data5;
  for(i=0;i<parti->ntimes;i++){
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
          ival = CLAMP((val-*valmin)/dg,0,NBUCKETS-1);
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
    FREEMEMORY(propi->buckets);
  }
  for(i=0;i<npart5prop;i++){
    part5prop *propi;

    propi = part5propinfo + i;

    switch(propi->setvalmin){
    case PERCENTILE_MIN:
      propi->valmin=propi->percentile_min;
      break;
    case GLOBAL_MIN:
      propi->valmin=propi->global_min;
      break;
    case SET_MIN:
      propi->valmin=propi->user_min;
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
    switch(propi->setvalmax){
    case PERCENTILE_MAX:
      propi->valmax=propi->percentile_max;
      break;
    case GLOBAL_MAX:
      propi->valmax=propi->global_max;
      break;
    case SET_MAX:
      propi->valmax=propi->user_max;
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
  }
  adjustpart5chops(parti);
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
        fprintf(stderr,"*** Error: Unable to allocate memory in getdatabounds\n");
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
        if(dp!=0.0f)level = CLAMP((int)((pdata[n] - *pmin)/dp),0,NBUCKETS-1);
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
    if(axislabels_smooth==1){
      smoothlabel(pmin,pmax,nrgb);
    }
}

/* ------------------ adjustPlot3Dbounds ------------------------ */

void adjustPlot3Dbounds(int plot3dvar, int setpmin, float *pmin, int setpmax, float *pmax)
{
    int nsmall, nbig, *buckets=NULL, n, level, total, alpha05;
    float dp, pmin2, pmax2;
    plot3ddata *p;
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
        fprintf(stderr,"*** Error: Unable to allocate memory in getdatabounds\n");
        return;
      }

      for (n=0;n<NBUCKETS;n++){
        buckets[n]=0;
      }
      for(i=0;i<nplot3dinfo;i++){
        p = plot3dinfo+i;
        if(p->loaded==0||p->display==0)continue;
        meshi=meshinfo+p->blocknumber;
        ntotal=(meshi->ibar+1)*(meshi->jbar+1)*(meshi->kbar+1);
        ndata += ntotal;
        iblank=meshi->c_iblank_node;
        q=meshi->qdata+plot3dvar*ntotal;
        for(n=0;n<ntotal;n++){
          if(iblank==NULL||*iblank++==GAS){
            level=0;
            if(dp!=0.0f)level = CLAMP((int)((q[n] - *pmin)/dp),0,NBUCKETS-1);
            buckets[level]++;
          }
        }
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
    if(axislabels_smooth==1&&setpmin!=SET_MIN&&setpmax!=SET_MAX){
      smoothlabel(pmin,pmax,nrgb);
    }

}

/* ------------------ getzonebounds ------------------------ */

void getzoneglobalbounds(const float *pdata, int ndata, float *pglobalmin, float *pglobalmax)
{
    int n;
    float pmin2, pmax2, val;

    pmin2 = pdata[0];
    pmax2 = pmin2;
    for (n=0;n<ndata;n++){
      val=*pdata++;
      pmin2=MIN(val,pmin2);
      pmax2=MAX(val,pmax2);
    }
    *pglobalmin = pmin2; 
    *pglobalmax = pmax2;
}

/* ------------------ smoothlabel ------------------------ */

void smoothlabel(float *a, float *b, int n){
  double delta, factor, logdelta;
  int ndigits;
  double half;

  half=0.5;

  delta = ((double)*b-(double)*a)/(double)(n-2);
  if(delta==0.0)return;
  logdelta = log10((double)delta);
  ndigits=logdelta-1;
  if(logdelta<=1)ndigits--;
  factor = 5*pow(10,ndigits);
  delta = (int)(delta/factor + half)*factor;

  *a = factor*(int)(*a/factor+half);
  *b = *a + (n-2)*delta;

}
