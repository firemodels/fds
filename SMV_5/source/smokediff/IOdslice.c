// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "svdiff.h"
#include "MALLOC.h"

// svn revision character string
char IOdslice_revision[]="$Revision$";

/* ------------------ setup_slice ------------------------ */

void setup_slice(FILE *stream_out){
  casedata *case1, *case2;
  int i;

  case1 = caseinfo;
  case2 = caseinfo + 1;

  for(i=0;i<case1->nslice_files;i++){
    slice *slicei;

    slicei = case1->sliceinfo + i;
    slicei->slice2 = getslice(slicei,case2);
    if(slicei->slice2!=NULL&&stream_out!=NULL){
      char outfile[1024];

      fprintf(stream_out,"%s\n",slicei->keyword);
      make_outfile(outfile,NULL,slicei->file,".sf");
      fprintf(stream_out,"%s\n",outfile);
      fprintf(stream_out,"%s\n",slicei->label.longlabel);
      fprintf(stream_out,"%s\n",slicei->label.shortlabel);
      fprintf(stream_out,"%s\n",slicei->label.unit);
    }
  }
}

/* ------------------ getslice ------------------------ */

slice *getslice(slice *slicein, casedata *case2){
  int i;
  float dx, dy, dz;

  dx = slicein->slicemesh->dx/2.0;
  dy = slicein->slicemesh->dy/2.0;
  dz = slicein->slicemesh->dz/2.0;

  for(i=0;i<case2->nslice_files;i++){
    slice *sliceout;

    sliceout = case2->sliceinfo + i;
    if(slicein->slicetype!=sliceout->slicetype)continue;
    if(strcmp(slicein->label.longlabel,sliceout->label.longlabel)!=0)continue;
    if(fabs(slicein->xmin-sliceout->xmin)>dx)continue;
    if(fabs(slicein->xmax-sliceout->xmax)>dx)continue;
    if(fabs(slicein->ymin-sliceout->ymin)>dy)continue;
    if(fabs(slicein->ymax-sliceout->ymax)>dy)continue;
    if(fabs(slicein->zmin-sliceout->zmin)>dz)continue;
    if(fabs(slicein->zmax-sliceout->zmax)>dz)continue;
    return sliceout;
  }
  return NULL;
}

/* ------------------ diff_slices ------------------------ */

void diff_slices(FILE *stream_out){
  int j;

  for(j=0;j<caseinfo->nslice_files;j++){
    float valmin, valmax;
    char *file1, *file2;
    char fullfile1[1024], fullfile2[1024], outfile[1024],  outfile2[1024];
    slice *slicei, *slice1, *slice2;
    FILE *stream;
    int unit1, unit2, unit3;
    FILE_SIZE len1,len2;
    int is1a, is2a, js1a, js2a, ks1a, ks2a;
    int is1b, is2b, js1b, js2b, ks1b, ks2b;
    int error1=0,error2a=0,error2b=0;
    float time1, *qframe1;
    int nqframe1;
    float time2a, *qframe2a;
    float time2b, *qframe2b;
    int nqframe2;
    float *qframeout;
    int i;
    int len;
    float f1, f2, dt;
    int slicetest1,slicetest2;
    float fraction_complete;
    FILE_SIZE size_sofar;
    int percent_complete;
    int nvals;
    float valmin_percentile, valmax_percentile, cdf01, cdf99;

    slicei = caseinfo->sliceinfo+j;
    slice1 = slicei;
    if(slicei->slice2==NULL)continue;
    slice2 = slicei->slice2;
    file1 = slicei->file;
    file2 = slicei->slice2->file;
    fullfile(fullfile1,sourcedir1,file1);
    fullfile(fullfile2,sourcedir2,file2);

    stream=fopen(fullfile1,"r");
    if(stream==NULL)continue;
    fclose(stream);

    stream=fopen(fullfile2,"r");
    if(stream==NULL)continue;
    fclose(stream);

    make_outfile(outfile,destdir,file1,".sf");
    if(strlen(outfile)==0)continue;
    stream=fopen(outfile,"w");
    if(stream==NULL)continue;
    fclose(stream);

    make_outfile(outfile2,NULL,slice1->file,".sf");

    unit1=11;
    len1=strlen(fullfile1);
    slicetest1=0;
    slicetest2=0;
    if(test_mode==1){
      slicetest1=1;
      slicetest2=2;
    }
    FORTopenslice(fullfile1,&unit1,&caseinfo->endian,&is1a,&is2a,&js1a,&js2a,&ks1a,&ks2a,&error1,len1);
    unit2=12;
    len2=strlen(fullfile2);
    FORTopenslice(fullfile2,&unit2,&caseinfo->endian,&is1b,&is2b,&js1b,&js2b,&ks1b,&ks2b,&error2a,len2);
    if(is1a!=is1b||js1a!=js1b||ks1a!=ks1b||
       is2a!=is2b||js2a!=js2b||ks2a!=ks2b||
       error1!=0||error2a!=0){
      FORTclosefortranfile(&unit1);
      FORTclosefortranfile(&unit2);
      if(error1!=0||error2a!=0){
        if(error1==0)printf("*** problem opening %s\n",fullfile1);
        if(error2a==0)printf("*** problem opening %s\n",fullfile2);
      }
      if(is1a!=is1b||js1a!=js1b||ks1a!=ks1b||
         is2a!=is2b||js2a!=js2b||ks2a!=ks2b){
        printf("*** integer slice bounds do not match for %s\n",fullfile1);
        printf("    %i %i %i %i %i %i\n",is1a, is2a, js1a, js2a, ks1a, ks2a);
        printf("    %i %i %i %i %i %i\n",is1b, is2b, js1b, js2b, ks1b, ks2b);
        printf(" %f %f %f %f %f %f\n",slice1->xmin,slice1->xmax,slice1->ymin,slice1->ymax,slice1->zmin,slice1->zmax);
        printf(" %f %f %f %f %f %f\n",slice2->xmin,slice2->xmax,slice2->ymin,slice2->ymax,slice2->zmin,slice2->zmax);
      }
      continue;
    }

    nqframe1 = (is2a+1-is1a)*(js2a+1-js1a)*(ks2a+1-ks1a);
    NewMemory((void **)&qframe1,nqframe1*sizeof(float));
    NewMemory((void **)&qframeout,nqframe1*sizeof(float));
    nqframe2 = (is2b+1-is1b)*(js2b+1-js1b)*(ks2b+1-ks1b);
    NewMemory((void **)&qframe2a,nqframe2*sizeof(float));
    NewMemory((void **)&qframe2b,nqframe2*sizeof(float));

    len=strlen(outfile);
    unit3=13;
    FORToutsliceheader(outfile,&unit3,&is1a,&is2a,&js1a,&js2a,&ks1a,&ks2a,&error1,len);
    if(error1!=0){
      FORTclosefortranfile(&unit1);
      FORTclosefortranfile(&unit2);
      printf("*** problem writing out header for %s\n",fullfile1);
      continue;
    }
    printf("Subtracting %s from %s\n",fullfile1,fullfile2);
    error1=1;
    error2a=1;
    error2b=1;
    nvals=0;
    update_data_hist(NULL, nvals, slice1->bucket, INIT_HISTOGRAM);
    FORTgetsliceframe(&unit1,&is1a,&is2a,&js1a,&js2a,&ks1a,&ks2a,&time1,qframe1,&slicetest1,&error1);
    if(error1==0)FORTgetsliceframe(&unit2,&is1b,&is2b,&js1b,&js2b,&ks1b,&ks2b,&time2a,qframe2a,&slicetest2,&error2a);
    if(error2a==0)FORTgetsliceframe(&unit2,&is1b,&is2b,&js1b,&js2b,&ks1b,&ks2b,&time2b,qframe2b,&slicetest2,&error2b);
    if(error1!=0||error2a!=0||error2b!=0){
      FORTclosefortranfile(&unit1);
      FORTclosefortranfile(&unit2);
      FORTclosefortranfile(&unit3);
      FREEMEMORY(qframe1);
      FREEMEMORY(qframe2a);
      FREEMEMORY(qframe2b);
      FREEMEMORY(qframeout);
      continue;
    }
    update_data_hist(qframe1, nqframe1, slice1->bucket, UPDATE_HISTOGRAM);
    printf("  Progress: ");
    fflush(stdout);

    percent_complete=0;
    size_sofar=0;
    valmin = 1000000000.0;
    valmax = -valmin;
    for(;;){

      size_sofar+=nqframe1*sizeof(float);
      fraction_complete=(float)size_sofar/(float)slice1->filesize;
      if((int)(fraction_complete*100)>percent_complete+10){
        if(percent_complete<100)percent_complete+=10;
        printf("%i%s ",percent_complete,pp);
        fflush(stdout);
      }
      while(time1>time2b){
        for(i=0;i<nqframe1;i++){
          qframe2a[i]=qframe2b[i];
        }
        time2a=time2b;
        FORTgetsliceframe(&unit2,&is1b,&is2b,&js1b,&js2b,&ks1b,&ks2b,&time2b,qframe2b,&slicetest2,&error2a);
        if(error2a!=0)break;
      }
      if(error2a!=0)break;
      dt = time2b - time2a;
      f1 = 1.0;
      f2 = 0.0;
      if(dt!=0.0){
        f1 = (time2b - time1)/dt;
        f2 = (time1-time2a)/dt;
      }
      for(i=0;i<nqframe1;i++){
        qframeout[i]=f1*qframe2a[i]+f2*qframe2b[i]-qframe1[i];
        if(qframe1[i]<valmin)valmin=qframe1[i];
        if(qframe1[i]>valmax)valmax=qframe1[i];
      }
      FORToutsliceframe(&unit3,&is1a,&is2a,&js1a,&js2a,&ks1a,&ks2a,&time1,qframeout,&error1);
      if(error1!=0)break;
      FORTgetsliceframe(&unit1,&is1a,&is2a,&js1a,&js2a,&ks1a,&ks2a,&time1,qframe1,&slicetest1,&error1);
      if(error1!=0)break;
      update_data_hist(qframe1, nqframe1, slice1->bucket, UPDATE_HISTOGRAM);
    }
    printf("\n");
    fflush(stdout);

    cdf01=0.01;
    cdf99=0.99;
    valmin_percentile = get_hist_val(slice1->bucket, cdf01);
    valmax_percentile = get_hist_val(slice1->bucket, cdf99);
    fprintf(stream_out,"MINMAXSLCF\n");
    fprintf(stream_out,"  %s\n",outfile2);
    fprintf(stream_out,"  %f %f %f %f\n",valmin,valmax,valmin_percentile,valmax_percentile);

    FORTclosefortranfile(&unit1);
    FORTclosefortranfile(&unit2);
    FORTclosefortranfile(&unit3);
    FREEMEMORY(qframe1);
    FREEMEMORY(qframe2a);
    FREEMEMORY(qframe2b);
    FREEMEMORY(qframeout);
  }
}
