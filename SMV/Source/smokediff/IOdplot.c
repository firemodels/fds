#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "svdiff.h"
#include "MALLOC.h"
#include "datadefs.h"
#include "file_util.h"

/* ------------------ setup_plot3d ------------------------ */

void setup_plot3d(FILE *stream_out){
  casedata *case1, *case2;
  int i;

  case1 = caseinfo;
  case2 = caseinfo + 1;

  for(i=0;i<case1->nplot3dinfo;i++){
    plot3d *plot3di;

    plot3di = case1->plot3dinfo + i;
    plot3di->plot3d2 = getplot3d(plot3di,case2);
    if(plot3di->plot3d2!=NULL&&stream_out!=NULL){
      char outfile[1024];
      int j;

      fprintf(stream_out,"%s\n",plot3di->keyword);
      make_outfile(outfile,NULL,plot3di->file,".q");
      fprintf(stream_out," %s\n",outfile);
      for(j=0;j<5;j++){
        flowlabels *label;

        label = plot3di->labels + j;
        fprintf(stream_out," %s\n",label->longlabel);
        fprintf(stream_out," %s\n",label->shortlabel);
        fprintf(stream_out," %s\n",label->unit);
      }
    }
  }
}

/* ------------------ getplot3d ------------------------ */

plot3d *getplot3d(plot3d *plot3din, casedata *case2){
  int i;
  float dx, dy, dz;
  meshdata *meshin;

  meshin = plot3din->plot3dmesh;
  dx = meshin->dx/2.0;
  dy = meshin->dy/2.0;
  dz = meshin->dz/2.0;

  for(i=0;i<case2->nplot3dinfo;i++){
    plot3d *plot3dout_local;
    meshdata *meshout;

    plot3dout_local = case2->plot3dinfo + i;
    meshout = plot3dout_local->plot3dmesh;
    if(ABS(plot3din->time-plot3dout_local->time)>0.05)continue;
    if(strcmp(plot3din->labels[0].longlabel,plot3dout_local->labels[0].longlabel)!=0)continue;
    if(strcmp(plot3din->labels[1].longlabel,plot3dout_local->labels[1].longlabel)!=0)continue;
    if(strcmp(plot3din->labels[2].longlabel,plot3dout_local->labels[2].longlabel)!=0)continue;
    if(strcmp(plot3din->labels[3].longlabel,plot3dout_local->labels[3].longlabel)!=0)continue;
    if(strcmp(plot3din->labels[4].longlabel,plot3dout_local->labels[4].longlabel)!=0)continue;
    if(ABS(meshin->xbar0-meshout->xbar0)>dx)continue;
    if(ABS(meshin->xbar-meshout->xbar)>dx)continue;
    if(ABS(meshin->ybar0-meshout->ybar0)>dy)continue;
    if(ABS(meshin->ybar-meshout->ybar)>dy)continue;
    if(ABS(meshin->zbar0-meshout->zbar0)>dz)continue;
    if(ABS(meshin->zbar-meshout->zbar)>dz)continue;
    return plot3dout_local;
  }
  return NULL;
}

/* ------------------ diff_plot3ds ------------------------ */

void diff_plot3ds(FILE *stream_out){
  int j;
  char *file1, *file2;
  char fullfile1[1024], fullfile2[1024], outfile[1024], outfile2[1024];
  float valmin, valmax;

  for(j=0;j<caseinfo->nplot3dinfo;j++){
    plot3d *plot3di, *plot3d1, *plot3d2;
    float *qframe1, *qframe2, *qout;
    meshdata *plot3dmesh;
    int nx, ny, nz, nq;
    int len1, len2, lenout;
    int isotest=0;
    FILE *stream;
    int error1, error2, error3;
    int endian;
    int i;
    int n;
    int nn;

    endian = getendian();
    plot3di = caseinfo->plot3dinfo+j;
    plot3d1 = plot3di;
    if(plot3di->plot3d2==NULL)continue;
    plot3d2 = plot3di->plot3d2;
    file1 = plot3d1->file;
    file2 = plot3d2->file;
    fullfile(fullfile1,sourcedir1,file1);
    fullfile(fullfile2,sourcedir2,file2);

    stream=fopen(fullfile1,"r");
    if(stream==NULL)continue;
    fclose(stream);

    stream=fopen(fullfile2,"r");
    if(stream==NULL)continue;
    fclose(stream);

    make_outfile(outfile,destdir,file1,".q");
    if(strlen(outfile)==0)continue;
    stream=fopen(outfile,"w");
    if(stream==NULL)continue;
    fclose(stream);

    make_outfile(outfile2,NULL,plot3d1->file,".q");

    plot3dmesh = plot3d1->plot3dmesh;
    nx = plot3dmesh->ibar+1;
    ny = plot3dmesh->jbar+1;
    nz = plot3dmesh->kbar+1;

    nq = 5*nx*ny*nz;
    NewMemory((void **)&qframe1,nq*sizeof(float));
    NewMemory((void **)&qframe2,nq*sizeof(float));
    NewMemory((void **)&qout,nq*sizeof(float));

    len1=strlen(fullfile1);
    len2=strlen(fullfile2);
    lenout=strlen(outfile);
    isotest=0;
    PRINTF("Subtracting %s from %s\n",fullfile2,fullfile1);
    FFLUSH();

    if(test_mode==1)isotest=1;
    PRINTF("  Progress: reading %s,",fullfile1);
    FFLUSH();

    FORTgetplot3dq(fullfile1,&nx,&ny,&nz,qframe1,&error1,&isotest,len1);
    if(test_mode==1)isotest=2;
    PRINTF(" reading %s,",fullfile2);
    FFLUSH();

    FORTgetplot3dq(fullfile2,&nx,&ny,&nz,qframe2,&error2,&isotest,len2);
    PRINTF(" differencing data,");
    FFLUSH();

    valmin=1000000000.0;
    valmax=-valmin;
    for(i=0;i<nq;i++){
      qout[i]=qframe1[i]-qframe2[i];
    }
    nn=0;
    fprintf(stream_out,"MINMAXPL3D\n");
    fprintf(stream_out,"  %s\n",outfile2);
    for(n=0;n<5;n++){
      int nvals;
      float valmin_percentile, valmax_percentile;

      valmin=1000000000.0;
      valmax=-valmin;

      for(i=0;i<nq/5;i++){
        if(qframe1[nn]<valmin)valmin=qframe1[nn];
        if(qframe1[nn]>valmax)valmax=qframe1[nn];
        nn++;
      }

      nvals=0;
      reset_histogram(plot3d1->histogram[n]);
      nvals=nq/5;
      update_histogram(qframe1+n*nvals, nvals, plot3d1->histogram[n]);
      valmin_percentile = get_histogram_value(plot3d1->histogram[n], 0.01);
      valmax_percentile = get_histogram_value(plot3d1->histogram[n], 0.99);
      fprintf(stream_out,"  %f %f %f %f\n",valmin,valmax,valmin_percentile,valmax_percentile);
    }

    FORTplot3dout(outfile,&nx,&ny,&nz,qout,&error3,lenout);
    PRINTF(" completed.\n");
    FFLUSH();

    FREEMEMORY(qframe1);
    FREEMEMORY(qframe2);
    FREEMEMORY(qout);
  }
}
