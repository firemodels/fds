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


/* ------------------ adjustpart5chops ------------------------ */

void adjustpart5chops(partdata *parti){
  int i;

  for(i=0;i<npart5prop;i++){
    partpropdata *propi;

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
  int i;

  if(parti->valmin==NULL){
    NewMemory((void **)&parti->valmin, npart5prop*sizeof(float));
  }
  if(parti->valmax==NULL){
    NewMemory((void **)&parti->valmax, npart5prop*sizeof(float));
  }

  for(i=0;i<npart5prop;i++){
    partpropdata *propi;
    histogramdata *histi;

    propi = part5propinfo + i;
    histi = &propi->histogram;

    propi->global_min = histi->valmin;
    propi->global_max = histi->valmax;

    propi->percentile_min = get_histogram_value(histi, percentile_level);
    propi->percentile_max = get_histogram_value(histi, 1.0 - percentile_level);

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
    parti->valmin[i] = propi->valmin;
    parti->valmax[i] = propi->valmax;
  }
  adjustpart5chops(parti);
#ifdef _DEBUG
  print_partprop();
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
    meshdata *meshi;
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

/* ------------------ getzoneglobalbounds ------------------------ */

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
