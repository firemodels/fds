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
  int sizebefore=0, sizeafter=0, size;
  int one=1, zero=0, fileversion=0, fdsversion;
  unsigned int *int_buffer_uncompressed;
  unsigned char *int_buffer_compressed;
  unsigned char *char_buffer_uncompressed,*char_buffer_compressed;
  mesh *partmesh;
  float xmin, xmax,  ymin, ymax, zmin, zmax;
  uLongf ncompressed_zlib;
  int percent_done;
  int percent_next=10;
  int data_loc;
  int count=0;
  int compression_level;

   //*** PARTICLE FILE FORMATS

  //*** FDS FORMAT (FORTRAN - each FORTRAN record has a 4 byte header and a 4 byte trailer surrounding the data)



  //*** ZLIB format (C - no extra bytes surrounding data) 

  //*** header
  // endian
  // completion (0/1)
  // fileversion (compressed particle file format version)
  // fds version 
  // compression level
  // global min max (used to perform conversion)
  // nclasses
  // quantities_1, ..., quantities_nclasses

  //*** data frame
  // time
  // points_1, ..., points_nclasses
  // ntotal_int, ntotal_char
  // ncompressedpoints
  // xyz_1, ..., xyz_ncompressedpoints
  // ncompresseddata
  // data_1, ..., data_ncompresseddata


  endiandata=getendian();
  if(endianswitch==1)endiandata=1-endiandata;

  partfile=parti->file;
  partmesh=parti->partmesh;
  xmin = partmesh->xbar0;
  xmax = partmesh->xbar;
  ymin = partmesh->ybar0;
  ymax = partmesh->ybar;
  zmin = partmesh->zbar0;
  zmax = partmesh->zbar;

#define COMPRESSION_LEVEL 1024
  compression_level=COMPRESSION_LEVEL;

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
  strcat(partsizefile_svz,".sz2");
  partsizestream=fopen(partsizefile_svz,"w");

  if(partsizestream==NULL){
    printf("The file %s could not be opened for writing\n",partsizefile_svz);
    fclose(PARTFILEstream);
    fclose(partstream);
    return;
  }

  printf("  Compressing %s\n",parti->file);

#define BUFFER_SIZE 1000000
  NewMemory((void **)&pdata,BUFFER_SIZE*sizeof(float));
  NewMemory((void **)&tagdata,BUFFER_SIZE*sizeof(int));
  NewMemory((void **)&int_buffer_uncompressed,BUFFER_SIZE*sizeof(unsigned int));
  NewMemory((void **)&int_buffer_compressed,BUFFER_SIZE*sizeof(unsigned char));
  NewMemory((void **)&char_buffer_uncompressed,BUFFER_SIZE*sizeof(unsigned char));
  NewMemory((void **)&char_buffer_compressed,BUFFER_SIZE*sizeof(unsigned char));

  lenfile=strlen(parti->file);
  FORTopenpart(parti->file,&unit,&endiandata,&error,lenfile);

  FORTgetpartheader1(&unit,&nclasses,&fdsversion,&size);
  NewMemory((void **)&nquantities,nclasses*sizeof(int));
  NewMemory((void **)&npoints,nclasses*sizeof(int));
  sizebefore+=size;

  FORTgetpartheader2(&unit,&nclasses,nquantities,&size);
  sizebefore+=size;

  fwrite(&one,4,1,partstream);           // write out a 1 to determine "endianness" when file is read in later
  fwrite(&zero,4,1,partstream);          // write out a zero now, then a one just before file is closed
  fwrite(&fileversion,4,1,partstream);   // write out compressed fileversion in case file format changes later
  fwrite(&fdsversion,4,1,partstream);    // fds file version
  fwrite(&compression_level,4,1,partstream);
  sizeafter=20;
  fwrite(&nclasses,4,1,partstream);      // compression level
  fwrite(nquantities,4,nclasses,partstream);
  sizeafter+=4*(1+nclasses);

  error=0;
  for(;;){
    float *x, *y, *z, *vals;
    int j;
    unsigned int *xbuff, *ybuff, *zbuff, *tbuff, *tagdataptr, ntotal_int, ntotal_char;
    unsigned char *cbuff;
    int ncompressed_int, ncompressed_char;

    FORTgetpartdataframe(&unit,&nclasses,nquantities,npoints,&time,tagdata,pdata,&size,&error);
    if(error!=0)break;

    fwrite(&time,4,1,partstream);
    fwrite(npoints,4,nclasses,partstream);
    sizeafter+=4*(1+nclasses);
    sizebefore+=size;

    fprintf(partsizestream,"%f\n",time);

    vals=pdata;
    tagdataptr=(unsigned int *)tagdata;
    xbuff = int_buffer_uncompressed;
    cbuff = char_buffer_uncompressed;
    ntotal_int=0;
    ntotal_char=0;
    for(j=0;j<nclasses;j++){
      part5class *classi;
      int k;

      fprintf(partsizestream," %i ",npoints[j]);
      if(j<nclasses-1){
        fprintf(partsizestream," \n");
      }

      if(npoints[j]==0)continue;
      ybuff = xbuff + npoints[j];
      zbuff = ybuff + npoints[j];
      tbuff = zbuff + npoints[j];
      classi=parti->classptr[j];
      x = vals;
      y = x + npoints[j];
      z = y + npoints[j];

      for(k=0;k<npoints[j];k++){
        xbuff[k] = compression_level*(x[k]-xmin)/(xmax-xmin);
        ybuff[k] = compression_level*(y[k]-ymin)/(ymax-ymin);
        zbuff[k] = compression_level*(z[k]-zmin)/(zmax-zmin);
        tbuff[k] = tagdataptr[k];
      }
      ntotal_int+=4*npoints[j];
      xbuff = tbuff + npoints[j];
      vals = z + npoints[j];
      tagdataptr += npoints[j];
      ntotal_char+=nquantities[j]*npoints[j];

      for(k=0;k<nquantities[j];k++){
        part5prop *propi;
        int cval,kk;
        float denom;
          
        propi=getpartprop(classi->labels[k].shortlabel);
        denom=propi->valmax-propi->valmin;
        if(denom<=0.0)denom=1.0;

        for(kk=0;kk<npoints[j];kk++){
          if(vals[kk]<propi->valmin){
            cval=0;
          }
          else if(vals[kk]>propi->valmax){
            cval=255;
          }
          else{
            cval = 1+253*(vals[kk]-propi->valmin)/denom;
            if(cval<1){
              cval=1;
            }
            else if(cval>254){
              cval=254;
            }
          }
          cbuff[kk]=cval;
        }
        cbuff += npoints[j];
        vals += npoints[j];
      }
    }

    // compress data

    fwrite(&ntotal_int,4,1,partstream);
    fwrite(&ntotal_char,4,1,partstream);

    ncompressed_int=0;
    if(ntotal_int>0){
      ncompressed_zlib=BUFFER_SIZE;
      compress(int_buffer_compressed, &ncompressed_zlib, (unsigned char *)int_buffer_uncompressed, 4*ntotal_int);
      ncompressed_int = ncompressed_zlib;
      sizeafter+=(4+ncompressed_int);
      fwrite(&ncompressed_int,4,1,partstream);
      fwrite(&int_buffer_compressed,1,ncompressed_int,partstream);
    }

    ncompressed_char=0;
    if(ntotal_char>0){
      ncompressed_zlib=BUFFER_SIZE;
      compress(char_buffer_compressed, &ncompressed_zlib, char_buffer_uncompressed, ntotal_char);
      ncompressed_char = ncompressed_zlib;
      sizeafter+=(4+ncompressed_char);
      fwrite(&ncompressed_char,4,1,partstream);
      fwrite(&char_buffer_compressed,1,ncompressed_char,partstream);
    }

    fprintf(partsizestream," %i %i %i %i\n",(int)ntotal_int,(int)ntotal_char,ncompressed_int,ncompressed_char);

    data_loc=sizebefore;
    percent_done=100.0*(float)data_loc/(float)parti->filesize;
    if(percent_done>percent_next){
      printf(" %i%s",percent_next,pp);
      fflush(stdout);
      percent_next+=10;
    }
    count++;
  }

  printf(" 100%s completed\n",pp);
  FREEMEMORY(nquantities);
  FREEMEMORY(npoints);
  FREEMEMORY(pdata);
  FREEMEMORY(tagdata);
  FREEMEMORY(int_buffer_uncompressed);
  FREEMEMORY(int_buffer_compressed);
  FREEMEMORY(char_buffer_uncompressed);
  FREEMEMORY(char_buffer_compressed);

  FORTclosefortranfile(&unit);
  fclose(partstream);
  fclose(partsizestream);
  {
    char before_label[256],after_label[256],xxx[2];

    strcpy(xxx,"X");
    getfilesizelabel(sizebefore,before_label);
    getfilesizelabel(sizeafter,after_label);
    printf("    records=%i, ",count);
    printf("Sizes: original=%s, ",before_label);
    printf("compressed=%s (%4.1f%s reduction)\n",after_label,(float)sizebefore/(float)sizeafter,xxx);
  }

}

/* ------------------ Get_Part_Bounds ------------------------ */

void Get_Part_Bounds(void){
  int i;
  float *pdata;
  int *tagdata;
  int fdsversion;

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
    int error, size;
    int nquantities_total;
    int j;

    parti = partinfo + i;
    printf("  Examining %s\n",parti->file);
    lenfile=strlen(parti->file);
    FORTopenpart(parti->file,&unit,&endiandata,&error1,lenfile);

    FORTgetpartheader1(&unit,&nclasses,&fdsversion,&size);
    NewMemory((void **)&nquantities,nclasses*sizeof(int));
    NewMemory((void **)&npoints,nclasses*sizeof(int));

    FORTgetpartheader2(&unit,&nclasses,nquantities,&size);
    nquantities_total=0;
    for(j=0;j<nclasses;j++){
      nquantities_total+=nquantities[j];
    }
    if(nquantities_total==0){
      FREEMEMORY(nquantities);
      FREEMEMORY(npoints);
      FORTclosefortranfile(&unit);
      continue;
    }

    error=0;
    for(;;){
      float *x, *y, *z, *vals;
      int k;

      FORTgetpartdataframe(&unit,&nclasses,nquantities,npoints,&time,tagdata,pdata,&size,&error);
      if(error!=0)break;

      vals=pdata;
      for(j=0;j<nclasses;j++){
        part5class *classj;

        if(npoints[j]==0)continue;
        classj=parti->classptr[j];
        x = vals;
        y = x + npoints[j];
        z = y + npoints[j];
        vals = z + npoints[j];
        for(k=0;k<nquantities[j];k++){
          part5prop *propi;
          
          propi=getpartprop(classj->labels[k].shortlabel);
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
