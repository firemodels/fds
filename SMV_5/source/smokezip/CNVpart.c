// $Date$ 
// $Revision$
// $Author$

#include "options.h"
// svn revision character string
char CNVpart_revision[]="$Revision$";
#ifdef pp_PART
#include "zlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "egz_stdio.h"
#include "svzip.h"
#include "MALLOC.h"



void endian_switch(void *val, int nval);
#define FORTREAD(var,size) fseek(PARTFILE,4,SEEK_CUR);\
                           returncode=fread(var,4,size,PARTFILE);\
                           if(returncode!=size)returncode=0;\
                           if(endianswitch==1&&returncode!=0)endian_switch(var,size);\
                           fseek(PARTFILE,4,SEEK_CUR)
#define READ(var,size) returncode=fread(var,4,size,PARTFILE);\
                       if(returncode!=size)returncode=0;\
                       if(endianswitch==1&&returncode!=0)endian_switch(var,size);
#define READBR(var,size) READ(var,size);if(returncode==0)break;
#define SEEK returncode=fseek(PARTFILE,4,SEEK_CUR)
#define SEEKBR SEEK;if(returncode!=0)break


/* ------------------ getpartprop ------------------------ */

part5prop *getpartprop(char *string){
  int i;
  part5prop *partpropi;

  for(i=0;i<npart5propinfo;i++){
    partpropi = part5propinfo + i;
    if(strcmp(partpropi->label.shortlabel,string)==0)return partpropi;
  }
  return NULL;
}

/* ------------------ compress_patches ------------------------ */

void compress_parts(void){
  int i;
  
  if(get_part_bounds==1){
    Get_Part_Bounds();
  }
  for(i=0;i<npart_files;i++){
    part *parti;

    parti = partinfo + i;
    if(autozip==1&&parti->autozip==0)continue;
    convert_part(parti);
  }
}

/* ------------------ compress_part ------------------------ */

void convert_part(part *parti){
  FILE *PARTFILEstream,*partstream,*partsizestream;
  char *partfile, partfile_svz[256], partsizefile_svz[256];
  int lenfile;
  int unit=15;
  int nclasses;
  int *nquantities, *npoints;
  float time;
  int error;
  int endiandata;
  int *tagdata;
  float *pdata;
  int i;

  endiandata=getendian();
  if(endianswitch==1)endiandata=1-endiandata;

  partfile=parti->file;

  // check if part file is accessible

  if(getfileinfo(partfile,NULL,NULL)!=0){
    printf("Particle file %s does not exist\n",partfile);
    return;
  }

  PARTFILEstream=fopen(partfile,"rb");
  if(PARTFILEstream==NULL){
    printf("Particle file %s could not be opened\n",partfile);
    return;
  }

  // set up compressed particle file

  strcpy(partfile_svz,"");
  if(destdir!=NULL){
    strcpy(partfile_svz,destdir);
  }
  strcat(partfile_svz,partfile);
  strcat(partfile_svz,".svz");
  partstream=fopen(partfile_svz,"wb");
  if(partstream==NULL){
    printf("The file %s could not be opened for writing\n",partfile_svz);
    fclose(PARTFILEstream);
    return;
  }

  // set up compressed particle size file

  strcpy(partsizefile_svz,"");
  if(destdir!=NULL){
    strcpy(partsizefile_svz,destdir);
  }
  strcpy(partsizefile_svz,partfile);
  strcat(partsizefile_svz,".sz");
  partsizestream=fopen(partsizefile_svz,"w");

  if(partsizestream==NULL){
    printf("The file %s could not be opened for writing\n",partsizefile_svz);
    fclose(PARTFILEstream);
    fclose(partstream);
    return;
  }

  printf("  Compressing %s\n",parti->file);

  NewMemory((void **)&pdata,1000000*sizeof(float));
  NewMemory((void **)&tagdata,1000000*sizeof(int));

  lenfile=strlen(parti->file);
  FORTopenpart(parti->file,&unit,&endiandata,&error,lenfile);

  FORTgetpartheader1(&unit,&nclasses);
  NewMemory((void **)&nquantities,nclasses*sizeof(int));
  NewMemory((void **)&npoints,nclasses*sizeof(int));

  FORTgetpartheader2(&unit,&nclasses,nquantities);

  error=0;
  for(;;){
    float *x, *y, *z, *vals;
    int j,k;

    FORTgetpartdataframe(&unit,&nclasses,nquantities,npoints,&time,tagdata,pdata,&error);
    if(error!=0)break;

    vals=pdata;
    for(j=0;j<nclasses;j++){
      part5class *classi;

      if(npoints[j]==0)continue;
      classi=parti->classptr[i];
      x = vals;
      y = x + npoints[j];
      z = y + npoints[j];
      vals = z + npoints[j];
      for(k=0;k<nquantities[j];k++){
        part5prop *propi;
          
        propi=getpartprop(classi->labels[k].shortlabel);

        vals += npoints[j];
      }
    }
  }

  FREEMEMORY(nquantities);
  FREEMEMORY(npoints);
  FREEMEMORY(pdata);
  FREEMEMORY(tagdata);

  FORTclosefortranfile(&unit);
}

/* ------------------ Get_Part_Bounds ------------------------ */

void Get_Part_Bounds(void){
  int i;
  float *pdata;
  int *tagdata;

  int endiandata;

  endiandata=getendian();
  if(endianswitch==1)endiandata=1-endiandata;

  printf("Determining particle file bounds\n");

  for(i=0;i<npart5propinfo;i++){
    part5prop *propi;

    propi = part5propinfo + i;

    NewMemory((void **)&propi->histogram,sizeof(histogramdata));
    init_histogram(propi->histogram);

  }

  NewMemory((void **)&pdata,1000000*sizeof(float));
  NewMemory((void **)&tagdata,1000000*sizeof(int));

  for(i=0;i<npart_files;i++){
    part *parti;
    FILE_SIZE lenfile;
    int unit=15;
    int error1;
    int nclasses;
    int *nquantities, *npoints;
    float time;
    int error;


    parti = partinfo + i;
    printf("  Examining %s\n",parti->file);
    lenfile=strlen(parti->file);
    FORTopenpart(parti->file,&unit,&endiandata,&error1,lenfile);

    FORTgetpartheader1(&unit,&nclasses);
    NewMemory((void **)&nquantities,nclasses*sizeof(int));
    NewMemory((void **)&npoints,nclasses*sizeof(int));

    FORTgetpartheader2(&unit,&nclasses,nquantities);

    error=0;
    for(;;){
      float *x, *y, *z, *vals;
      int j,k;

      FORTgetpartdataframe(&unit,&nclasses,nquantities,npoints,&time,tagdata,pdata,&error);
      if(error!=0)break;

      vals=pdata;
      for(j=0;j<nclasses;j++){
        part5class *classi;

        if(npoints[j]==0)continue;
        classi=parti->classptr[i];
        x = vals;
        y = x + npoints[j];
        z = y + npoints[j];
        vals = z + npoints[j];
        for(k=0;k<nquantities[j];k++){
          part5prop *propi;
          
          propi=getpartprop(classi->labels[k].shortlabel);
          update_histogram(vals,npoints[j],propi->histogram);

          vals += npoints[j];
        }
      }
    }

    FREEMEMORY(nquantities);
    FREEMEMORY(npoints);

    FORTclosefortranfile(&unit);
  }

  FREEMEMORY(pdata);
  FREEMEMORY(tagdata);

  for(i=0;i<npart5propinfo;i++){
    part5prop *propi;

    propi = part5propinfo + i;

    propi->valmax=get_histogram_value(propi->histogram,0.99);
    propi->valmin=get_histogram_value(propi->histogram,0.01);
    propi->setvalmax=1;
    propi->setvalmin=1;
    printf(" %s min: %f max: %f\n",propi->label.shortlabel,propi->valmin,propi->valmax);
    FREEMEMORY(propi->histogram);
  }

}
#endif
