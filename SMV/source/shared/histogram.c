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

float get_histogram_value(histogramdata *histgram, float cdf){

// get value of histogram for value cdf

  int cutoff, count;
  int i;
  float returnval;

  if(cdf<=0.0){
    return histgram->valmin;
  }
  if(cdf>=1.0){
    return histgram->valmax;
  }
  cutoff = cdf*histgram->ntotal;
  count=0;
  for(i=0;i<NHIST_BUCKETS;i++){
    count+=histgram->buckets[i];
    if(count>cutoff){
      returnval = histgram->valmin + (float)(i+0.5)*(histgram->valmax-histgram->valmin)/(float)NHIST_BUCKETS;
      return returnval;
    }
  }
  return histgram->valmax;
}

/* ------------------ complete_histogram ------------------------ */

void complete_histogram(histogramdata *histgram){

// set variable indicating that histogram is complete

  histgram->complete=1;
}

/* ------------------ init_histogram ------------------------ */

void init_histogram(histogramdata *histgram){

// initialize histogram data structures

  int i,nbuckets;

  nbuckets = NHIST_BUCKETS;
  histgram->ndim = 1;
  histgram->nbuckets = nbuckets;
  for(i = 0; i<nbuckets; i++){
    histgram->buckets[i]=0;
  }
  histgram->defined=0;
  histgram->ntotal=0;
  histgram->valmin=(float)pow(10.0,20.0);
  histgram->valmax=-histgram->valmin;
  histgram->complete=0;
}

/* ------------------ init_histogram2d ------------------------ */

void init_histogram2d(histogramdata *histgram, int nx, int ny){

// initialize histogram data structures

  int i, nbuckets;

  nbuckets = nx*ny;
  histgram->ndim = 2;
  histgram->nbuckets = nbuckets;
  NewMemory((void **)&histgram->buckets2, histgram->nbuckets*sizeof(int));
  for(i = 0; i<nbuckets; i++){
    histgram->buckets2[i]=0;
  }
  histgram->defined=0;
  histgram->ntotal=0;
  histgram->valxmin=(float)pow(10.0,20.0);
  histgram->valxmax=-histgram->valxmin;
  histgram->valymin=(float)pow(10.0,20.0);
  histgram->valymax=-histgram->valymin;
  histgram->complete=0;
}

/* ------------------ copy_data2histogram ------------------------ */

void copy_data2histogram(float *vals, int nvals, histogramdata *histgram){

// copy nvals of the floating point array, vals, into the histogram histgram

  int i;
  float valmin, valmax;
  float dbucket;

  histgram->defined=1;
  for(i=0;i<NHIST_BUCKETS;i++){
    histgram->buckets[i]=0;
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
      histgram->buckets[0]=nvals;
    }
    else{
      for(i=0;i<nvals;i++){
        int ival;

        ival = (vals[i]-valmin)/dbucket;
        ival=MAX(0,ival);
        ival=MIN(NHIST_BUCKETS-1,ival);
        histgram->buckets[ival]++;
      }
    }
  }
  histgram->ntotal=nvals;
  histgram->valmax=valmax;
  histgram->valmin=valmin;
}

/* ------------------ copy_uvdata2histogram ------------------------ */

void copy_uvdata2histogram(float *uvals, float *vvals, int nvals, histogramdata *histgram){
  int i;

  for(i = 0; i<nvals; i++){
    float r, theta;
    float u, v;

    u = uvals[i];
    v = vvals[i];
    r = sqrt(u*u+v*v);
    theta = RAD2DEG*atan2(v, u);
  }
}


/* ------------------ update_histogram ------------------------ */

void update_histogram(float *vals, int nvals, histogramdata *histgram){

// merge nvals of the floating point array, vals, into the histogram histgram

  histogramdata histgramval;

  if(nvals<=0)return;
  copy_data2histogram(vals,nvals,&histgramval);
  merge_histogram(histgram,&histgramval);
}

/* ------------------ merge_histogram ------------------------ */

void merge_histogram(histogramdata *histgram1, histogramdata *histgram2){
  
  // merge histogram histgram2 into histgram1

  int i;
  float dbucket1, dbucket2, dbucket_new;
  int bucket1copy[NHIST_BUCKETS];
  float valmin_new, valmax_new;

  histgram1->defined=1;
  valmin_new=MIN(histgram1->valmin,histgram2->valmin);
  valmax_new=MAX(histgram1->valmax,histgram2->valmax);

  for(i=0;i<NHIST_BUCKETS;i++){
    bucket1copy[i]=histgram1->buckets[i];
    histgram1->buckets[i]=0;
  }
  dbucket1 = (histgram1->valmax-histgram1->valmin)/NHIST_BUCKETS;
  dbucket2 = (histgram2->valmax-histgram2->valmin)/NHIST_BUCKETS;
  dbucket_new=(valmax_new-valmin_new)/NHIST_BUCKETS;

  if(dbucket_new==0.0){
    histgram1->buckets[0]=histgram1->ntotal+histgram2->ntotal;
    histgram1->ntotal=histgram1->buckets[0];
  }
  else{
    for(i=0;i<NHIST_BUCKETS;i++){
      float val;
      int ival;

      if(bucket1copy[i]!=0){
        val = (float)(histgram1->valmin + (float)(i+0.5)*dbucket1);
        valmin_new=MIN(valmin_new,val);
        valmax_new=MAX(valmax_new,val);
        ival = (val-valmin_new)/dbucket_new;
        if(ival<0)ival=0;
        if(ival>NHIST_BUCKETS-1)ival=NHIST_BUCKETS-1;
        histgram1->buckets[ival]+=bucket1copy[i];
      }
      if(histgram2->buckets[i]!=0){
        val = (float)(histgram2->valmin + (i+0.5)*dbucket2);
        valmin_new=MIN(valmin_new,val);
        valmax_new=MAX(valmax_new,val);
        ival = (val-valmin_new)/dbucket_new;
        if(ival<0)ival=0;
        if(ival>NHIST_BUCKETS-1)ival=NHIST_BUCKETS-1;
        histgram1->buckets[ival]+=histgram2->buckets[i];
      }
    }
  }
  histgram1->valmin=valmin_new;
  histgram1->valmax=valmax_new;
  histgram1->ntotal+=histgram2->ntotal;
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
