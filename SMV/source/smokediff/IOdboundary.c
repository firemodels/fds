#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "svdiff.h"
#include "MALLOC.h"
#include "file_util.h"

/* ------------------ setup_slice ------------------------ */

void setup_boundary(FILE *stream_out){
  casedata *case1, *case2;
  int i;

  case1 = caseinfo;
  case2 = caseinfo + 1;

  for(i=0;i<case1->nboundary_files;i++){
    boundary *boundaryi;

    boundaryi = case1->boundaryinfo + i;
    boundaryi->boundary2 = getboundary(boundaryi,case2);
    if(boundaryi->boundary2!=NULL&&stream_out!=NULL){
      char outfile[1024];

      fprintf(stream_out,"%s\n",boundaryi->keyword);
      make_outfile(outfile,NULL,boundaryi->file,".bf");
      fprintf(stream_out," %s\n",outfile);
      fprintf(stream_out," %s\n",boundaryi->label.longlabel);
      fprintf(stream_out," %s\n",boundaryi->label.shortlabel);
      fprintf(stream_out," %s\n",boundaryi->label.unit);
    }
  }
}

/* ------------------ getboundary ------------------------ */

boundary *getboundary(boundary *boundaryin, casedata *case2){
  int i,j;

  if(strlen(type_label)>0&&strcmp(type_label,boundaryin->label.shortlabel)!=0){
    return NULL;
  }
  for(i=0;i<case2->nboundary_files;i++){
    boundary *boundaryout;
    int *patch2index;
    int nmatches=0;

    boundaryout = case2->boundaryinfo + i;
    if(boundaryin->boundarytype!=boundaryout->boundarytype)continue;
    if(strcmp(boundaryin->label.longlabel,boundaryout->label.longlabel)!=0)continue;
    if(mesh_match(boundaryin->boundarymesh,boundaryout->boundarymesh)==0)continue;
    patch2index=boundaryin->patch2index;
    for(j=0;j<boundaryin->npatches;j++){
      patch2index[j]=getpatchindex(j,boundaryin,boundaryout);
      if(patch2index[j]!=-1)nmatches++;
    }
    if(nmatches==0)continue;

    return boundaryout;
  }
  return NULL;
}

/* ------------------ getpatchindex ------------------------ */

int getpatchindex(int in1, boundary *boundaryin, boundary *boundaryout){
  int j;
  int i1, i2, j1, j2, k1, k2;
  int *pi1, *pi2, *pj1, *pj2, *pk1, *pk2;
  int pd1, *patchdir2;

  i1 = boundaryin->pi1[in1];
  i2 = boundaryin->pi2[in1];
  j1 = boundaryin->pj1[in1];
  j2 = boundaryin->pj2[in1];
  k1 = boundaryin->pk1[in1];
  k2 = boundaryin->pk2[in1];
  pd1 = boundaryin->patchdir[in1];

  pi1 = boundaryout->pi1;
  pi2 = boundaryout->pi2;
  pj1 = boundaryout->pj1;
  pj2 = boundaryout->pj2;
  pk1 = boundaryout->pk1;
  pk2 = boundaryout->pk2;
  patchdir2 = boundaryout->patchdir;

  for(j=0;j<boundaryout->npatches;j++){
    if(i1!=pi1[j]||i2!=pi2[j])continue;
    if(j1!=pj1[j]||j2!=pj2[j])continue;
    if(k1!=pk1[j]||k2!=pk2[j])continue;
    if(pd1!=patchdir2[j])continue;
    return j;
  }
  return -1;
}

/* ------------------ diff_slices ------------------------ */

void diff_boundaryes(FILE *stream_out){
  int j;

  for(j=0;j<caseinfo->nboundary_files;j++){
    char *file1, *file2;
    char fullfile1[1024], fullfile2[1024], outfile[1024], outfile2[1024];
    boundary *boundaryi, *boundary1, *boundary2;
    FILE *stream;
    int unit1, unit2, unit3;
    int error1, error2, error3;
    FILE_SIZE len1,len2;
    int *p1i1, *p1i2, *p1j1, *p1j2, *p1k1, *p1k2;
    int *p2i1, *p2i2, *p2j1, *p2j2, *p2k1, *p2k2;
    int *p3i1, *p3i2, *p3j1, *p3j2, *p3k1, *p3k2;
    int *patchdir3;
    float *pqq1, *pqq2a, *pqq2b, *pqq2out, *pqq3;
    float patchtime1, patchtime2a, patchtime2b;
    int nsize1, nsize2, nsize3;
    int npatches3;
    int i;
    int npqq1, npqq2a, npqq2b;

    boundaryi = caseinfo->boundaryinfo+j;
    boundary1 = boundaryi;
    if(boundaryi->boundary2==NULL)continue;
    boundary2 = boundaryi->boundary2;
    file1 = boundary1->file;
    file2 = boundary2->file;
    fullfile(fullfile1,sourcedir1,file1);
    fullfile(fullfile2,sourcedir2,file2);

    p1i1 = boundary1->pi1;
    p1i2 = boundary1->pi2;
    p1j1 = boundary1->pj1;
    p1j2 = boundary1->pj2;
    p1k1 = boundary1->pk1;
    p1k2 = boundary1->pk2;
    nsize1=0;
    nsize3=0;
    npatches3=0;
    for(i=0;i<boundary1->npatches;i++){
      nsize1+=boundary1->patchsize[i];
      if(boundary1->patch2index[i]!=-1){
        npatches3++;
        nsize3+=boundary1->patchsize[i];
      }
    }
    p2i1 = boundary2->pi1;
    p2i2 = boundary2->pi2;
    p2j1 = boundary2->pj1;
    p2j2 = boundary2->pj2;
    p2k1 = boundary2->pk1;
    p2k2 = boundary2->pk2;

    nsize2=0;
    for(i=0;i<boundary2->npatches;i++){
      nsize2+=boundary2->patchsize[i];
    }

    NewMemory((void **)&pqq1,nsize1*sizeof(float));
    NewMemory((void **)&pqq2a,nsize2*sizeof(float));
    NewMemory((void **)&pqq2b,nsize2*sizeof(float));
    NewMemory((void **)&pqq2out,nsize2*sizeof(float));
    NewMemory((void **)&pqq3,nsize3*sizeof(float));

    NewMemory((void **)&p3i1,npatches3*sizeof(int));
    NewMemory((void **)&p3i2,npatches3*sizeof(int));
    NewMemory((void **)&p3j1,npatches3*sizeof(int));
    NewMemory((void **)&p3j2,npatches3*sizeof(int));
    NewMemory((void **)&p3k1,npatches3*sizeof(int));
    NewMemory((void **)&p3k2,npatches3*sizeof(int));
    NewMemory((void **)&patchdir3,npatches3*sizeof(int));

    stream=fopen(fullfile1,"r");
    if(stream==NULL)continue;
    fclose(stream);

    stream=fopen(fullfile2,"r");
    if(stream==NULL)continue;
    fclose(stream);

    make_outfile(outfile,destdir,file1,".bf");
    if(strlen(outfile)==0)continue;
    stream=fopen(outfile,"w");
    if(stream==NULL)continue;
    fclose(stream);
    make_outfile(outfile2,NULL,boundary1->file,".bf");

    PRINTF("Subtracting %s from %s\n",fullfile2,fullfile1);

    unit1=11;
    FORTget_file_unit(&unit1,&unit1);
    len1=strlen(fullfile1);
    FORTopenboundary(fullfile1,&unit1,&boundary1->version,&error1,len1);

    unit2=12;
    FORTget_file_unit(&unit2,&unit2);
    len2=strlen(fullfile2);
    FORTopenboundary(fullfile2,&unit2,&boundary2->version,&error2,len2);

    if(error1==0&&error2==0){
      int ii,len3;
      FILE_SIZE size_sofar;
      float fraction_complete;
      int percent_complete;
      float valmin, valmax;
      float valmin_percentile, valmax_percentile;

      ii=0;
      for(i=0;i<boundary1->npatches;i++){
        int jj;

        jj = boundary1->patch2index[i];
        if(jj==-1)continue;

        p3i1[ii]=p1i1[i];
        p3i2[ii]=p1i2[i];
        p3j1[ii]=p1j1[i];
        p3j2[ii]=p1j2[i];
        p3k1[ii]=p1k1[i];
        p3k2[ii]=p1k2[i];
        patchdir3[ii]=boundary1->patchdir[i];

        ii++;
      }
      unit3=15;
      FORTget_file_unit(&unit3,&unit3);
      len3=strlen(outfile);
      size_sofar=0;
      FORToutboundaryheader(outfile,&unit3,&npatches3,
        p3i1,p3i2,p3j1,p3j2,p3k1,p3k2,patchdir3,&error1,len3);
      PRINTF("  Progress: ");
      FFLUSH();
      percent_complete=0;
      valmin=1000000000.0;
      valmax=-valmin;

      reset_histogram(boundary1->histogram);

      FORTgetpatchdata(&unit1, &boundary1->npatches,
        p1i1, p1i2, p1j1, p1j2, p1k1, p1k2, &patchtime1, pqq1, &npqq1, &error1);
      FORTgetpatchdata(&unit2, &boundary2->npatches,
        p2i1, p2i2, p2j1, p2j2, p2k1, p2k2, &patchtime2a, pqq2a, &npqq2a, &error2);
      if(error2==0)FORTgetpatchdata(&unit2, &boundary2->npatches,
        p2i1, p2i2, p2j1, p2j2, p2k1, p2k2, &patchtime2b, pqq2b, &npqq2b, &error2);
      for(;;){
        int iq;
        float f1, f2, dt;

        if(error1!=0||error2!=0)break;

        while(patchtime1>patchtime2b){
          for(i=0;i<nsize2;i++){
            pqq2a[i]=pqq2b[i];
          }
          patchtime2a=patchtime2b;
          FORTgetpatchdata(&unit2, &boundary2->npatches,
            p2i1, p2i2, p2j1, p2j2, p2k1, p2k2, &patchtime2b, pqq2b, &npqq2b, &error2);
          if(error2!=0)break;
        }
        if(error2!=0)break;
        dt = patchtime2b - patchtime2a;
        f1 = 1.0;
        f2 = 0.0;
        if(dt!=0.0){
          f1 = (patchtime2b - patchtime1)/dt;
          f2 = (patchtime1-patchtime2a)/dt;
        }
        for(i=0;i<nsize2;i++){
          pqq2out[i]=f1*pqq2a[i]+f2*pqq2b[i];
        }

        iq=0;
        for(i=0;i<boundary1->npatches;i++){
          int jj, kk;
          float *pq1, *pq2;

          jj = boundary1->patch2index[i];
          if(jj==-1)continue;

          pq1 = pqq1 + boundary1->qoffset[i];
          pq2 = pqq2out + boundary2->qoffset[jj];
          for(kk=0;kk<boundary1->patchsize[i];kk++){
            pqq3[iq++]=pq1[kk]-pq2[kk];
            if(pq1[kk]<valmin)valmin=pq1[kk];
            if(pq1[kk]>valmax)valmax=pq1[kk];
          }
        }
        update_histogram(pqq1, nsize1, boundary1->histogram);
        FORToutpatchframe(&unit3, &npatches3,
                        p3i1, p3i2, p3j1, p3j2, p3k1, p3k2,
                        &patchtime1, pqq3, &error3);
        size_sofar+=nsize1*sizeof(float);
        fraction_complete=(float)size_sofar/(float)boundary1->filesize;
        if((int)(fraction_complete*100)>percent_complete+10){
          if(percent_complete<100)percent_complete+=10;
          PRINTF("%i%s ",percent_complete,pp);
          FFLUSH();
        }

        FORTgetpatchdata(&unit1, &boundary1->npatches,
          p1i1, p1i2, p1j1, p1j2, p1k1, p1k2, &patchtime1, pqq1, &npqq1, &error1);
      }
      PRINTF("\n");
      FFLUSH();
      valmin_percentile = get_histogram_value(boundary1->histogram, 0.01);
      valmax_percentile = get_histogram_value(boundary1->histogram, 0.99);

      fprintf(stream_out,"MINMAXBNDF\n");
      fprintf(stream_out,"  %s\n",outfile2);
      fprintf(stream_out,"  %f %f %f %f\n",valmin,valmax,valmin_percentile,valmax_percentile);
    }

    FREEMEMORY(pqq1);
    FREEMEMORY(pqq2a);
    FREEMEMORY(pqq2b);
    FREEMEMORY(pqq2out);
    FREEMEMORY(pqq3);
    FREEMEMORY(p3i1);
    FREEMEMORY(p3i2);
    FREEMEMORY(p3j1);
    FREEMEMORY(p3j2);
    FREEMEMORY(p3k1);
    FREEMEMORY(p3k2);
    FREEMEMORY(patchdir3);

    if(error1!=0)FORTclosefortranfile(&unit1);
    if(error2!=0)FORTclosefortranfile(&unit2);
  }
}
