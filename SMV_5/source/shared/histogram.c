// $Date$ 
// $Revision$
// $Author$

#define IN_BUCKET 1

#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "histogram.h"

// svn revision character string
char histogram_revision[]="$Revision$";

#define HMIN(X, Y)  ((X) < (Y) ? (X) : (Y))
#define HMAX(X, Y)  ((X) > (Y) ? (X) : (Y))

/* ------------------ get_histogram_value ------------------------ */

float get_histogram_value(histogramdata *histgram, float cdf){
  int cutoff, count;
  int i;
  float returnval;

  cutoff = cdf*histgram->ntotal;
  count=0;
  for(i=0;i<NBUCKETS;i++){
    count+=histgram->buckets[i];
    if(count>cutoff){
      returnval = histgram->valmin + (float)(i+0.5)*(histgram->valmax-histgram->valmin)/(float)NBUCKETS;
      return returnval;
    }
  }
  return histgram->valmax;
}

/* ------------------ init_histogram ------------------------ */

void init_histogram(float *vals, int nvals, histogramdata *histgram){
  int i;
  float valmin, valmax;
  float dbucket;

  for(i=0;i<NBUCKETS;i++){
    histgram->buckets[i]=0;
  }
  if(vals==NULL||nvals==0){
    histgram->ntotal=0;
    histgram->valmin=pow(10.0,20.0);
    histgram->valmax=-histgram->valmin;
  }
  else{
    valmin=vals[0];
    valmax=vals[0];
    for(i=1;i<nvals;i++){
      valmin=HMIN(vals[i],valmin);
      valmax=HMAX(vals[i],valmax);
    }
    dbucket=(valmax-valmin)/NBUCKETS;
    if(dbucket==0.0){
      histgram->buckets[0]=nvals;
    }
    else{
      for(i=0;i<nvals;i++){
        int ival;

        ival = (vals[i]-valmin)/dbucket;
        ival=HMAX(0,ival);
        ival=HMIN(NBUCKETS-1,ival);
        histgram->buckets[ival]++;
      }
    }
    histgram->ntotal=nvals;
    histgram->valmax=valmax;
    histgram->valmin=valmin;
  }
}

/* ------------------ update_histogram ------------------------ */

void update_histogram(float *vals, int nvals, histogramdata *histgram){
  histogramdata histgramval;

  init_histogram(vals,nvals,&histgramval);
  merge_histogram(histgram,&histgramval);
}

/* ------------------ merge_histogram ------------------------ */

void merge_histogram(histogramdata *histgram1, histogramdata *histgram2){
  
  // merge histogram histgram2 into histgram1

  int i;
  float dbucket1, dbucket2, dbucket_new;
  int bucket1copy[NBUCKETS];
  float valmin_new, valmax_new;

  valmin_new=HMIN(histgram1->valmin,histgram2->valmin);
  valmax_new=HMAX(histgram1->valmax,histgram2->valmax);

  for(i=0;i<NBUCKETS;i++){
    bucket1copy[i]=histgram1->buckets[i];
    histgram1->buckets[i]=0;
  }
  dbucket1 = (histgram1->valmax-histgram1->valmin)/NBUCKETS;
  dbucket2 = (histgram2->valmax-histgram2->valmin)/NBUCKETS;
  dbucket_new=(valmax_new-valmin_new)/NBUCKETS;

  if(dbucket_new==0.0){
    histgram1->buckets[0]=histgram1->ntotal+histgram2->ntotal;
    histgram1->ntotal=histgram1->buckets[0];
  }
  else{
    for(i=0;i<NBUCKETS;i++){
      float val;
      int ival;

      if(bucket1copy[i]!=0){
        val = histgram1->valmin + (i+0.5)*dbucket1;
        ival = (val-valmin_new)/dbucket_new;
        if(ival<0)ival=0;
        if(ival>NBUCKETS-1)ival=NBUCKETS-1;
        histgram1->buckets[ival]+=bucket1copy[i];
      }
      if(histgram2->buckets[i]!=0){
        val = histgram2->valmin + (i+0.5)*dbucket2;
        ival = (val-valmin_new)/dbucket_new;
        if(ival<0)ival=0;
        if(ival>NBUCKETS-1)ival=NBUCKETS-1;
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

  vals=malloc(NVALS*sizeof(float));

  for(i=0;i<NVALS;i++){
    vals[i]=10.0*rand()/(float)RAND_MAX;
  }
  vals2histogram(vals,NVALS,&histogram);
  v01=get_histogram_value(&histogram, .01);
  v50=get_histogram_value(&histogram, .50);
  v99=get_histogram_value(&histogram, .99);
  free(vals);

}
#endif
