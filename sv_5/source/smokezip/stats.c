#include "options.h"
#include "zlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "egz_stdio.h"
#include "svzip.h"
#include "MALLOC.h"

#ifdef WIN32
#pragma warning (disable:4244)		/* disable bogus conversion warnings */
#endif

// svn revision character string
char stats_revision[]="$Revision: 624 $";

pdfdata pdftemp;

float minval(float f1, float f2);
float maxval(float f1, float f2);


/* ------------------ getpdf ------------------------ */

void getpdf(float *vals, int nvals, pdfdata *pdf){
  int i;

  float pmin, pmax;
  float dval;

  initpdf(pdf);
  if(nvals<=0)return;
  pdf->ncount=nvals;

  pmin=vals[0];
  pmax=pmin;
  for(i=1;i<nvals;i++){
    if(vals[i]<pmin)pmin=vals[i];
    if(vals[i]>pmax)pmax=vals[i];
  }
  pdf->pdfmax=pmax;
  pdf->pdfmin=pmin;
  dval=(pmax-pmin)/PDFMAX;
  if(dval==0.0){
    pdf->buckets[0]+=nvals;
    return;
  }
  for(i=0;i<nvals;i++){
    int ival;

    ival = (vals[i]-pmin)/dval;
    if(ival<0)ival=0;
    if(ival>PDFMAX-1)ival=PDFMAX-1;
    pdf->buckets[ival]++;
  }
  return;
}

/* ------------------ mergepdf ------------------------ */

void mergepdf(pdfdata *pdf1, pdfdata *pdf2, pdfdata *pdfmerge){
  int i;
  float dvalmerge,dval1,dval2;
  float val1, val2;
  int ival1, ival2;
  int count1, count2;

  if(pdf1->ncount==0&&pdf2->ncount==0){
    initpdf(pdfmerge);
    return;
  }
  initpdf(&pdftemp);
  if(pdf1->ncount>0&&pdf2->ncount>0){
    pdftemp.pdfmin=minval(pdf1->pdfmin,pdf2->pdfmin);
    pdftemp.pdfmax=maxval(pdf1->pdfmax,pdf2->pdfmax);
    pdftemp.ncount=pdf1->ncount+pdf2->ncount;
  }
  else{
    if(pdf1->ncount>0){
      memcpy(&pdftemp,pdf1,sizeof(pdfdata));
    }
    else{
      memcpy(&pdftemp,pdf2,sizeof(pdfdata));
    }
    memcpy(pdfmerge,&pdftemp,sizeof(pdfdata));
    return;
  }
  dvalmerge = (pdftemp.pdfmax-pdftemp.pdfmin)/(PDFMAX-1);
  dval1 = (pdf1->pdfmax-pdf1->pdfmin)/(PDFMAX-1);
  dval2 = (pdf2->pdfmax-pdf2->pdfmin)/(PDFMAX-1);

  for(i=0;i<PDFMAX;i++){
    count1=pdf1->buckets[i];
    if(count1!=0){
      val1 = pdf1->pdfmin+i*dval1;
      if(dvalmerge!=0.0){
        ival1 = (val1-pdftemp.pdfmin)/dvalmerge;
      }
      else{
        ival1=0;
      }
      if(ival1<0)ival1=0;
      if(ival1>PDFMAX-1)ival1=PDFMAX-1;
      pdftemp.buckets[ival1]+=count1;
    }
    count2=pdf2->buckets[i];
    if(count2!=0){
      val2 = pdf2->pdfmin+i*dval2;
      if(dvalmerge!=0.0){
        ival2 = (val2-pdftemp.pdfmin)/dvalmerge;
      }
      else{
        ival2=0;
      }
      if(ival2<0)ival2=0;
      if(ival2>PDFMAX-1)ival2=PDFMAX-1;
      pdftemp.buckets[ival2]+=count2;
    }
  }
  memcpy(pdfmerge,&pdftemp,sizeof(pdfdata));
}

/* ------------------ initpdf ------------------------ */

void initpdf(pdfdata *pdf){
  int i;

  pdf->ncount=0;
  for(i=0;i<PDFMAX;i++){
    pdf->buckets[i]=0;
  }
  pdf->pdfmax=0.0;
  pdf->pdfmin=0.0;
}

/* ------------------ maxval ------------------------ */

float maxval(float val1,float val2){
  if(val1>val2)return val1;
  return val2;
}

/* ------------------ minval ------------------------ */

float minval(float val1, float val2){
  if(val1<val2)return val1;
  return val2;
}
