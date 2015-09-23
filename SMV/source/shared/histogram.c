#include "options.h"
#include "lint.h"

#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "histogram.h"
#include "pragmas.h"
#include "MALLOC.h"
#include "datadefs.h"

/* ------------------ get_histogram_value ------------------------ */

float get_histogram_value(histogramdata *histogram, float cdf){

// get value of histogram for value cdf

  int cutoff, count;
  int i;
  float returnval;

  if(cdf<=0.0){
    return histogram->valmin;
  }
  if(cdf>=1.0){
    return histogram->valmax;
  }
  cutoff = cdf*histogram->ntotal;
  count=0;
  for(i=0;i<NHIST_BUCKETS;i++){
    count+=histogram->buckets[i];
    if(count>cutoff){
      returnval = histogram->valmin + (float)(i+0.5)*(histogram->valmax-histogram->valmin)/(float)NHIST_BUCKETS;
      return returnval;
    }
  }
  return histogram->valmax;
}

/* ------------------ complete_histogram ------------------------ */

void complete_histogram(histogramdata *histogram){

// set variable indicating that histogram is complete

  histogram->complete=1;
}

/* ------------------ init_histogram ------------------------ */

void init_histogram(histogramdata *histogram){

// initialize histogram data structures

  int i,nbuckets;

  nbuckets = NHIST_BUCKETS;
  histogram->ndim = 1;
  histogram->nbuckets = nbuckets;
  for(i = 0; i<nbuckets; i++){
    histogram->buckets[i]=0;
  }
  histogram->defined=0;
  histogram->ntotal=0;
  histogram->valmin=(float)pow(10.0,20.0);
  histogram->valmax=-histogram->valmin;
  histogram->complete=0;
}

/* ------------------ init_histogram2d ------------------------ */

void init_histogram2d(histogramdata *histogram, int nx, int ny){

// initialize histogram data structures

  int i, nbuckets;

  nbuckets = nx*ny;
  histogram->ndim = 2;
  histogram->nbuckets = nbuckets;
  NewMemory((void **)&histogram->buckets2, histogram->nbuckets*sizeof(int));
  for(i = 0; i<nbuckets; i++){
    histogram->buckets2[i]=0;
  }
  histogram->defined=0;
  histogram->ntotal=0;
  histogram->valxmin=(float)pow(10.0,20.0);
  histogram->valxmax=-histogram->valxmin;
  histogram->valymin=(float)pow(10.0,20.0);
  histogram->valymax=-histogram->valymin;
  histogram->complete=0;
}

/* ------------------ copy_data2histogram ------------------------ */

void copy_data2histogram(float *vals, int nvals, histogramdata *histogram){

// copy nvals of the floating point array, vals, into the histogram histogram

  int i;
  float valmin, valmax;
  float dbucket;

  histogram->defined=1;
  for(i=0;i<NHIST_BUCKETS;i++){
    histogram->buckets[i]=0;
  }
  if(nvals==0){
    valmin=(float)pow(10.0,20.0);
    valmax=-valmin;
  }
  else{
    valmin=vals[0];
    valmax=vals[0];
    for(i=1;i<nvals;i++){
      valmin=MIN(vals[i],valmin);
      valmax=MAX(vals[i],valmax);
    }
    dbucket=(valmax-valmin)/NHIST_BUCKETS;
    if(dbucket==0.0){
      histogram->buckets[0]=nvals;
    }
    else{
      for(i=0;i<nvals;i++){
        int ival;

        ival = (vals[i]-valmin)/dbucket;
        ival=MAX(0,ival);
        ival=MIN(NHIST_BUCKETS-1,ival);
        histogram->buckets[ival]++;
      }
    }
  }
  histogram->ntotal=nvals;
  histogram->valmax=valmax;
  histogram->valmin=valmin;
}

/* ------------------ copy_uvdata2histogram ------------------------ */

void copy_uvdata2histogram(float *uvals, float *vvals, int nvals, histogramdata *histogram){
  int i;
  float rmin, rmax;

  NewMemory((void **)&histogram->rvals,nvals*sizeof(float));
  for(i = 0; i < nvals; i++){
    float r, theta;
    float u, v;
    int ix, iy, ixy;

    u = uvals[i];
    v = vvals[i];
    histogram->rvals[i] = sqrt(u*u + v*v);
  }

  copy_data2histogram(histogram->rvals, nvals, histogram);
  rmin = histogram->valmin;
  rmax = histogram->valmax;

  for(i = 0; i < nvals; i++){
    float r, theta;
    float u, v;
    int ix, iy, ixy;

    u = uvals[i];
    v = vvals[i];
    r = histogram->rvals[i];
    ix = 0;
    if(rmax>rmin)ix = CLAMP(histogram->nx*(r - rmin) / (rmax - rmin),0,histogram->nx-1);
    theta = RAD2DEG*atan2(v, u);
    if(theta < 0.0)theta += 360.0;
    iy = CLAMP(histogram->ny*(theta / 360.0),0,histogram->ny-1);
    ixy = ix + iy*histogram->nx;
    histogram->buckets2[ixy]++;
  }
  FREEMEMORY(histogram->rvals);
}

/* ------------------ update_histogram ------------------------ */

void update_histogram(float *vals, int nvals, histogramdata *histogram){

// merge nvals of the floating point array, vals, into the histogram histogram

  histogramdata histogramval;

  if(nvals<=0)return;
  copy_data2histogram(vals,nvals,&histogramval);
  merge_histogram(histogram,&histogramval);
}

/* ------------------ merge_histogram ------------------------ */

void merge_histogram(histogramdata *histogram1, histogramdata *histogram2){
  
  // merge histogram histogram2 into histogram1

  int i;
  float dbucket1, dbucket2, dbucket_new;
  int bucket1copy[NHIST_BUCKETS];
  float valmin_new, valmax_new;

  histogram1->defined=1;
  valmin_new=MIN(histogram1->valmin,histogram2->valmin);
  valmax_new=MAX(histogram1->valmax,histogram2->valmax);

  for(i=0;i<NHIST_BUCKETS;i++){
    bucket1copy[i]=histogram1->buckets[i];
    histogram1->buckets[i]=0;
  }
  dbucket1 = (histogram1->valmax-histogram1->valmin)/NHIST_BUCKETS;
  dbucket2 = (histogram2->valmax-histogram2->valmin)/NHIST_BUCKETS;
  dbucket_new=(valmax_new-valmin_new)/NHIST_BUCKETS;

  if(dbucket_new==0.0){
    histogram1->buckets[0]=histogram1->ntotal+histogram2->ntotal;
    histogram1->ntotal=histogram1->buckets[0];
  }
  else{
    for(i=0;i<NHIST_BUCKETS;i++){
      float val;
      int ival;

      if(bucket1copy[i]!=0){
        val = (float)(histogram1->valmin + (float)(i+0.5)*dbucket1);
        valmin_new=MIN(valmin_new,val);
        valmax_new=MAX(valmax_new,val);
        ival = (val-valmin_new)/dbucket_new;
        if(ival<0)ival=0;
        if(ival>NHIST_BUCKETS-1)ival=NHIST_BUCKETS-1;
        histogram1->buckets[ival]+=bucket1copy[i];
      }
      if(histogram2->buckets[i]!=0){
        val = (float)(histogram2->valmin + (i+0.5)*dbucket2);
        valmin_new=MIN(valmin_new,val);
        valmax_new=MAX(valmax_new,val);
        ival = (val-valmin_new)/dbucket_new;
        if(ival<0)ival=0;
        if(ival>NHIST_BUCKETS-1)ival=NHIST_BUCKETS-1;
        histogram1->buckets[ival]+=histogram2->buckets[i];
      }
    }
  }
  histogram1->valmin=valmin_new;
  histogram1->valmax=valmax_new;
  histogram1->ntotal+=histogram2->ntotal;
}

#ifdef pp_CHECK
/* ------------------ check_histogram ------------------------ */

void check_histogram(void){
#define NVALS 1000000
  float *vals;
  int i;
  histogramdata histogram;
  float v01, v50, v99;

  NewMemory((void **)&vals,NVALS*sizeof(float));

  for(i=0;i<NVALS;i++){
    vals[i]=10.0*rand()/(float)RAND_MAX;
  }
  vals2histogram(vals,NVALS,&histogram);
  v01=get_histogram_value(&histogram, .01);
  v50=get_histogram_value(&histogram, .50);
  v99=get_histogram_value(&histogram, .99);
  FREEMEMORY(vals);

}
#endif
