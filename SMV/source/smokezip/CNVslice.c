#include "options.h"
#include "zlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <zlib.h>
#include "svzip.h"
#include "MALLOC.h"
#include "compress.h"

void mt_update_slice_hist(void);

#define FORTSLICEREAD(var,size) FSEEK(SLICEFILE,4,SEEK_CUR);\
                           returncode=fread(var,4,size,SLICEFILE);\
                           if(endianswitch==1)endian_switch(var,size);\
                           FSEEK(SLICEFILE,4,SEEK_CUR)

/* ------------------ convert_volslice ------------------------ */

int convert_volslice(slice *slicei, int *thread_index){
  char slicefile_svz[1024];
  char *slice_file;
  char filetype[1024];
  char *shortlabel;
  int ijkbar[6];
  uLong framesize;
  float *sliceframe_data=NULL;
  int sizebefore, sizeafter;
  int returncode;
  LINT data_loc;
  int percent_done;
  int percent_next=10;
#ifndef pp_THREAD
  int count=0;
#endif
  FILE *SLICEFILE;
  FILE *slicestream;

#ifdef pp_THREAD
  if(GLOBcleanfiles==0){
    int fileindex;

    fileindex = slicei + 1 - sliceinfo;
    sprintf(threadinfo[*thread_index].label,"vsf %i",fileindex);
  }
#endif

  slice_file=slicei->file;

  // check if slice file is accessible

  strcpy(filetype,"");
  shortlabel=slicei->label.shortlabel;
  if(strlen(shortlabel)>0)strcat(filetype,shortlabel);
  TrimBack(filetype);

  if(getfileinfo(slice_file,NULL,NULL)!=0){
    fprintf(stderr,"*** Warning: The file %s does not exist\n",slice_file);
    return 0;
  }

  SLICEFILE=fopen(slice_file,"rb");
  if(SLICEFILE==NULL){
    fprintf(stderr,"*** Warning: The file %s could not be opened\n",slice_file);
    return 0;
  }

  // set up slice compressed file

  if(GLOBdestdir!=NULL){
    strcpy(slicefile_svz,GLOBdestdir);
    strcat(slicefile_svz,slicei->filebase);
  }
  else{
    strcpy(slicefile_svz,slicei->file);
  }

  if(strlen(slicefile_svz)>4)strcat(slicefile_svz,".svv");

  if(GLOBcleanfiles==1){
    slicestream=fopen(slicefile_svz,"rb");
    if(slicestream!=NULL){
      fclose(slicestream);
      PRINTF("  Removing %s\n",slicefile_svz);
      UNLINK(slicefile_svz);
      LOCK_COMPRESS;
      GLOBfilesremoved++;
      UNLOCK_COMPRESS;
    }
    fclose(SLICEFILE);
    return 0;
  }

  if(GLOBoverwrite_slice==0){
    slicestream=fopen(slicefile_svz,"rb");
    if(slicestream!=NULL){
      fclose(slicestream);
      fprintf(stderr,"*** Warning: The file %s exists.\n",slicefile_svz);
      fprintf(stderr,"     Use the -f option to overwrite smokezip compressed files\n");
      fclose(SLICEFILE);
      return 0;
    }
  }

  slicestream=fopen(slicefile_svz,"wb");
  if(slicestream==NULL){
    fprintf(stderr,"*** Warning: The file %s could not be opened for writing\n",slicefile_svz);
    fclose(SLICEFILE);
    return 0;
  }

  // read and write slice header

#ifndef pp_THREAD
  if(GLOBcleanfiles==0){
    PRINTF("Compressing %s (%s)\n",slice_file,filetype);
  }
#endif


  {
    int skip;

    skip = 3*(4+30+4);  // skip over 3 records each containing a 30 byte FORTRAN character string
    returncode=FSEEK(SLICEFILE,skip,SEEK_CUR);
    sizebefore=skip;
  }

  FORTSLICEREAD(ijkbar,6);
  sizebefore+=4+6*4+4;
  sizeafter=0;

  {
    int one=1, version_local=0, completion=0;

    fwrite(&one,4,1,slicestream);
    fwrite(&version_local,4,1,slicestream);
    fwrite(&completion,4,1,slicestream);
  }


  {
    int ni, nj, nk;

    ni = ijkbar[1]+1-ijkbar[0];
    nj = ijkbar[3]+1-ijkbar[2];
    nk = ijkbar[5]+1-ijkbar[4];
    framesize = ni*nj*nk;
    NewMemory((void **)&sliceframe_data,framesize*sizeof(float));

    for(;;){
      float vmin, vmax;
      float *valmin, *valmax;
      unsigned char *compressed_data_out;
      uLongf ncompressed_data_out;
      float time_local;

      FORTSLICEREAD(&time_local,1);
      if(returncode==0)break;
      CheckMemory;
      sizebefore+=12;

      FORTSLICEREAD(sliceframe_data,framesize);    //---------------
      if(returncode==0)break;
      CheckMemory;
      sizebefore+=(4+framesize*sizeof(float)+4);

      valmin=NULL;
      valmax=NULL;
      if(slicei->voltype==1){
        vmin=0.0;
        valmin=&vmin;
      }
      else if(slicei->voltype==2){
        vmin=20.0;
        valmin=&vmin;
        vmax=1400.0;
        valmax=&vmax;
      }
      else{
        ASSERT(0);
      }
      CheckMemory;
      compress_volsliceframe(sliceframe_data, framesize, time_local, valmin, valmax,
                &compressed_data_out, &ncompressed_data_out);
      CheckMemory;
      sizeafter+=ncompressed_data_out;
      if(ncompressed_data_out>0){
        fwrite(compressed_data_out,1,ncompressed_data_out,slicestream);
      }
      CheckMemory;
      FREEMEMORY(compressed_data_out);

#ifndef pp_THREAD
      count++;
#endif

      data_loc=FTELL(SLICEFILE);
      percent_done=100.0*(float)data_loc/(float)slicei->filesize;
#ifdef pp_THREAD
      threadinfo[*thread_index].stat=percent_done;
      if(percent_done>percent_next){
        LOCK_PRINT;
        print_thread_stats();
        UNLOCK_PRINT;
        percent_next+=10;
      }
#else
      if(percent_done>percent_next){
        PRINTF(" %i%s",percent_next,GLOBpp);
        FFLUSH();
        percent_next+=10;
      }
#endif

    }
    if(returncode!=0){
      fprintf(stderr,"*** Error: compress returncode=%i\n",returncode);
    }
    FREEMEMORY(sliceframe_data);
  }

#ifndef pp_THREAD
    PRINTF(" 100%s completed\n",GLOBpp);
#endif

  {
    int completion=1;

    FSEEK(slicestream,4,SEEK_SET);
    fwrite(&completion,4,1,slicestream);
  }
  fclose(SLICEFILE);
  fclose(slicestream);

  {
    char before_label[256],after_label[256];

    getfilesizelabel(sizebefore,before_label);
    getfilesizelabel(sizeafter,after_label);
#ifdef pp_THREAD
    slicei->vol_compressed=1;
    sprintf(slicei->volsummary,"compressed from %s to %s (%4.1f%s reduction)",before_label,after_label,(float)sizebefore/(float)sizeafter,GLOBx);
    threadinfo[*thread_index].stat=-1;
#else
    PRINTF("  records=%i, ",count);
    PRINTF("Sizes: original=%s, ",before_label);
    PRINTF("compressed=%s (%4.1f%s reduction)\n\n",after_label,(float)sizebefore/(float)sizeafter,GLOBx);
#endif
  }

  return 1;

}

/* ------------------ convert_slice ------------------------ */

// unsigned int uncompress_rle(unsigned char *buffer_in, int nchars_in, unsigned char *buffer_out)

int convert_slice(slice *slicei, int *thread_index){

  char slicefile_svz[1024], slicesizefile_svz[1024];
  int fileversion, one, zero;
  char *slice_file;
  int version_local;
  char filetype[1024];
  char *shortlabel, *unit;
  char units[256];
  int ijkbar[6];
  uLong framesize;
  float *sliceframe_data=NULL;
  unsigned char *sliceframe_compressed=NULL, *sliceframe_uncompressed=NULL;
  unsigned char *sliceframe_compressed_rle=NULL, *sliceframe_uncompressed_rle=NULL;
  char cval[256];
  int sizebefore, sizeafter;
  int returncode;
  float minmax[2];
  float time_local;
  LINT data_loc;
  int percent_done;
  int percent_next=10;
  float valmin, valmax, denom;
  int chop_min, chop_max;
  uLongf ncompressed_zlib;
  int ncompressed_save;
#ifndef pp_THREAD
  int count=0;
#endif
  int ncol, nrow, idir;
  float time_max;
  int itime;
  LINT file_loc;

  FILE *SLICEFILE;
  FILE *slicestream,*slicesizestream;

#ifdef pp_THREAD
  if(GLOBcleanfiles==0){
    int fileindex;

    fileindex = slicei + 1 - sliceinfo;
    sprintf(threadinfo[*thread_index].label,"sf %i",fileindex);
  }
#endif

  slice_file=slicei->file;
  version_local=slicei->version;

  fileversion = 1;
  one = 1;
  zero=0;

  // check if slice file is accessible

  strcpy(filetype,"");
  shortlabel=slicei->label.shortlabel;
  if(strlen(shortlabel)>0)strcat(filetype,shortlabel);
  TrimBack(filetype);

  if(getfileinfo(slice_file,NULL,NULL)!=0){
    fprintf(stderr,"*** Warning: The file %s does not exist\n",slice_file);
    return 0;
  }

  SLICEFILE=fopen(slice_file,"rb");
  if(SLICEFILE==NULL){
    fprintf(stderr,"*** Warning: The file %s could not be opened\n",slice_file);
    return 0;
  }

  // set up slice compressed file

  if(GLOBdestdir!=NULL){
    strcpy(slicefile_svz,GLOBdestdir);
    strcat(slicefile_svz,slicei->filebase);
  }
  else{
    strcpy(slicefile_svz,slicei->file);
  }
  {

    char *ext;
    int lensvz;

    lensvz = strlen(slicefile_svz);

    if(lensvz>4){
      ext = slicefile_svz + lensvz - 4;
      if(strcmp(ext,".rle")==0){
        slicefile_svz[lensvz-4]=0;
      }
      strcat(slicefile_svz,".svz");
    }
  }

  if(GLOBdestdir!=NULL){
    strcpy(slicesizefile_svz,GLOBdestdir);
    strcat(slicesizefile_svz,slicei->filebase);
  }
  else{
    strcpy(slicesizefile_svz,slicei->file);
  }
  {

    char *ext;
    int lensvz;

    lensvz = strlen(slicesizefile_svz);

    if(lensvz>4){
      ext = slicesizefile_svz + lensvz - 4;
      if(strcmp(ext,".rle")==0){
        slicesizefile_svz[lensvz-4]=0;
      }
      strcat(slicesizefile_svz,".sz");
    }
  }

  if(GLOBcleanfiles==1){
    slicestream=fopen(slicefile_svz,"rb");
    if(slicestream!=NULL){
      fclose(slicestream);
      PRINTF("  Removing %s\n",slicefile_svz);
      UNLINK(slicefile_svz);
      LOCK_COMPRESS;
      GLOBfilesremoved++;
      UNLOCK_COMPRESS;
    }
    slicesizestream=fopen(slicesizefile_svz,"rb");
    if(slicesizestream!=NULL){
      fclose(slicesizestream);
      PRINTF("  Removing %s\n",slicesizefile_svz);
      UNLINK(slicesizefile_svz);
      LOCK_COMPRESS;
      GLOBfilesremoved++;
      UNLOCK_COMPRESS;
    }
    fclose(SLICEFILE);
    return 0;
  }

  if(GLOBoverwrite_slice==0){
    slicestream=fopen(slicefile_svz,"rb");
    if(slicestream!=NULL){
      fclose(slicestream);
      fprintf(stderr,"*** Warning:  %s exists.\n",slicefile_svz);
      fprintf(stderr,"     Use the -f option to overwrite smokezip compressed files\n");
      fclose(SLICEFILE);
      return 0;
    }
  }

  slicestream=fopen(slicefile_svz,"wb");
  slicesizestream=fopen(slicesizefile_svz,"w");
  if(slicestream==NULL||slicesizestream==NULL){
    if(slicestream==NULL){
      fprintf(stderr,"*** Warning: The file %s could not be opened for writing\n",slicefile_svz);
    }
    if(slicesizestream==NULL){
      fprintf(stderr,"  %s could not be opened for writing\n",slicesizefile_svz);
    }
    if(slicestream!=NULL)fclose(slicestream);
    if(slicesizestream!=NULL)fclose(slicesizestream);
    fclose(SLICEFILE);
    return 0;
  }

  // read and write slice header

  strcpy(units,"");
  unit=slicei->label.unit;
  if(strlen(unit)>0)strcat(units,unit);
  TrimBack(units);
  sprintf(cval,"%f",slicei->valmin);
  TrimZeros(cval);
#ifndef pp_THREAD
  if(GLOBcleanfiles==0){
    PRINTF("Compressing %s (%s)\n",slice_file,filetype);
    PRINTF("  using min=%s %s",cval,units);
  }
#endif
  sprintf(cval,"%f",slicei->valmax);
  TrimZeros(cval);
#ifndef pp_THREAD
  if(GLOBcleanfiles==0){
    PRINTF(" max=%s %s\n",cval,units);
    PRINTF(" ");
  }
#endif
  valmin=slicei->valmin;
  valmax=slicei->valmax;
  denom = valmax-valmin;
  if(denom==0.0)denom=1.0;

  chop_min=0;
  chop_max=255;
  if(GLOBno_chop==0){
    if(slicei->setchopvalmax==1){
        chop_max = 255*(slicei->chopvalmax-valmin)/denom;
        if(chop_max<0)chop_max=0;
        if(chop_max>255)chop_max=255;
    }
    if(slicei->setchopvalmin==1){
       chop_min = 255*(slicei->chopvalmin-valmin)/denom;
       if(chop_min<0)chop_min=0;
       if(chop_min>255)chop_min=255;
    }
  }


  fwrite(&one,4,1,slicestream);           // write out a 1 to determine "endianness" when file is read in later
  fwrite(&zero,4,1,slicestream);          // write out a zero now, then a one just before file is closed
  fwrite(&fileversion,4,1,slicestream);   // write out compressed fileversion in case file format changes later
  fwrite(&version_local,4,1,slicestream);       // fds slice file version
  sizeafter=16;

  //*** SLICE FILE FORMATS

  //*** FDS FORMAT (FORTRAN - each FORTRAN record has a 4 byte header and a 4 byte trailer surrounding the data)

  // 30 byte long label
  // 30 byte short label
  // 30 byte unit
  // i1,i2,j1,j2,k1,k2

  // for each time step:

  // time, compressed frame size
  // qq(1,nbuffer)              where nbuffer = (i2+1-i1)*(j2+1-j1)*(k2+1-k1)



  //*** ZLIB format (C - no extra bytes surrounding data)

  //*** header
  // endian
  // completion (0/1)
  // fileversion (compressed format)
  // version_local  (slicef version)
  // global min max (used to perform conversion)
  // i1,i2,j1,j2,k1,k2


  //*** frame
  // time, compressed frame size                        for each frame
  // compressed buffer


  //*** RLE format (FORTRAN)

  //*** header
  // endian
  // fileversion, slice version
  // global min max (used to perform conversion)
  // i1,i2,j1,j2,k1,k2


  //*** frame
  // time
  // compressed frame size                        for each frame
  // compressed buffer

  {
    int skip;

    skip = 3*(4+30+4);  // skip over 3 records each containing a 30 byte FORTRAN character string
    returncode=FSEEK(SLICEFILE,skip,SEEK_CUR);
    sizebefore=skip;
  }

  FORTSLICEREAD(ijkbar,6);
  sizebefore+=8+6*4;

  framesize =  (ijkbar[1]+1-ijkbar[0]);
  framesize *= (ijkbar[3]+1-ijkbar[2]);
  framesize *= (ijkbar[5]+1-ijkbar[4]);

  minmax[0]=slicei->valmin;
  minmax[1]=slicei->valmax;
  fwrite(minmax,4,2,slicestream);    // min max vals
  fwrite(ijkbar,4,6,slicestream);
  sizeafter+=(8+24);


  ncompressed_save=1.02*framesize+600;
  if(NewMemory((void **)&sliceframe_data,ncompressed_save*sizeof(float))==0)goto wrapup;
  if(NewMemory((void **)&sliceframe_compressed,ncompressed_save*sizeof(unsigned char))==0)goto wrapup;
  if(NewMemory((void **)&sliceframe_uncompressed,ncompressed_save*sizeof(unsigned char))==0)goto wrapup;

  fprintf(slicesizestream,"%i %i %i %i %i %i\n",ijkbar[0],ijkbar[1],ijkbar[2],ijkbar[3],ijkbar[4],ijkbar[5]);
  fprintf(slicesizestream,"%f %f\n",minmax[0],minmax[1]);

  idir=0;
  if(ijkbar[0]==ijkbar[1]){
    idir=1;
    ncol = ijkbar[3] + 1 - ijkbar[2];
    nrow = ijkbar[5] + 1 - ijkbar[4];
  }
  else if(ijkbar[2]==ijkbar[3]){
    idir=2;
    ncol = ijkbar[1] + 1 - ijkbar[0];
    nrow = ijkbar[5] + 1 - ijkbar[4];
  }
  else if(ijkbar[4]==ijkbar[5]){
    idir=3;
    ncol = ijkbar[1] + 1 - ijkbar[0];
    nrow = ijkbar[3] + 1 - ijkbar[2];
  }
  if(idir==0){
    idir=1;
    ncol = ijkbar[3] + 1 - ijkbar[2];
    nrow = ijkbar[5] + 1 - ijkbar[4];
  }


  {
    int ni, nj, nk;

    ni = ijkbar[1]+1-ijkbar[0];
    nj = ijkbar[3]+1-ijkbar[2];
    nk = ijkbar[5]+1-ijkbar[4];

    time_max=-1000000.0;
    itime=-1;
    for(;;){
      int i;

      FORTSLICEREAD(&time_local,1);
      sizebefore+=12;
      if(returncode==0)break;
      FORTSLICEREAD(sliceframe_data,framesize);    //---------------
      if(returncode==0)break;

      sizebefore+=(8+framesize*4);
      if(time_local<time_max)continue;
      time_max=time_local;

#ifndef pp_THREAD
      count++;
#endif

      data_loc=FTELL(SLICEFILE);
      percent_done=100.0*(float)data_loc/(float)slicei->filesize;
#ifdef pp_THREAD
      threadinfo[*thread_index].stat=percent_done;
      if(percent_done>percent_next){
        LOCK_PRINT;
        print_thread_stats();
        UNLOCK_PRINT;
        percent_next+=10;
      }
#else
      if(percent_done>percent_next){
        PRINTF(" %i%s",percent_next,GLOBpp);
        FFLUSH();
        percent_next+=10;
      }
#endif
      for(i=0;i<framesize;i++){
        int ival;
        int icol, jrow, index2;
        int ii,jj,kk;

        // val_in(i,j,k) = i + j*ni + k*ni*nj

        if(framesize<=ncol*nrow){  // only one slice plane

          // i = jrow*ncol + icol;

          icol = i%ncol;
          jrow = i/ncol;

          index2 = icol*nrow + jrow;
        }
        else{
          ii = i%ni;
          jj = (i/ni)%nj;
          kk = i/(ni*nj);

          index2 = ii*nj*nk + jj*nk + kk;
        }

        {
          float val;

          val = sliceframe_data[i];
          if(val<valmin){
            ival=0;
          }
          else if(val>valmax){
            ival=255;
          }
          else{
            ival = 1 + 253*(val-valmin)/denom;
          }
          if(ival<chop_min)ival=0;
          if(ival>chop_max)ival=255;
          sliceframe_uncompressed[index2] = ival;
        }
      }
      itime++;
      if(itime%GLOBslicezipstep!=0)continue;

      //int compress (Bytef *dest,   uLongf *destLen, const Bytef *source, uLong sourceLen);
      ncompressed_zlib=ncompressed_save;
      returncode=compress_zlib(sliceframe_compressed,&ncompressed_zlib,sliceframe_uncompressed,framesize);

      file_loc=FTELL(slicestream);
      fwrite(&time_local,4,1,slicestream);
      fwrite(&ncompressed_zlib,4,1,slicestream);
      fwrite(sliceframe_compressed,1,ncompressed_zlib,slicestream);
      sizeafter+=(8+ncompressed_zlib);
      fprintf(slicesizestream,"%f %i, %li\n",time_local,(int)ncompressed_zlib,(long)file_loc);
    }
    if(returncode!=0){
      fprintf(stderr,"*** Error: compress returncode=%i\n",returncode);
    }
  }

wrapup:
#ifndef pp_THREAD
    PRINTF(" 100%s completed\n",GLOBpp);
#endif
  FREEMEMORY(sliceframe_data);
  FREEMEMORY(sliceframe_compressed);
  FREEMEMORY(sliceframe_uncompressed);
  FREEMEMORY(sliceframe_compressed_rle);
  FREEMEMORY(sliceframe_uncompressed_rle);

  fclose(SLICEFILE);
  FSEEK(slicestream,4,SEEK_SET);
  fwrite(&one,4,1,slicestream);  // write completion code
  fclose(slicestream);
  fclose(slicesizestream);

  {
    char before_label[256],after_label[256];

    getfilesizelabel(sizebefore,before_label);
    getfilesizelabel(sizeafter,after_label);
#ifdef pp_THREAD
    slicei->compressed=1;
    sprintf(slicei->summary,"compressed from %s to %s (%4.1f%s reduction)",before_label,after_label,(float)sizebefore/(float)sizeafter,GLOBx);
    threadinfo[*thread_index].stat=-1;
#else
    PRINTF("  records=%i, ",count);
    PRINTF("Sizes: original=%s, ",before_label);
    PRINTF("compressed=%s (%4.1f%s reduction)\n\n",after_label,(float)sizebefore/(float)sizeafter,GLOBx);
#endif
  }

  return 1;
}

/* ------------------ getslice ------------------------ */

slice *getslice(char *string){
  int i;
  slice *slicei;

  for(i=0;i<nsliceinfo;i++){
    slicei = sliceinfo + i;
    if(slicei->dup==1)continue;
    if(strcmp(slicei->label.shortlabel,string)==0)return slicei;
  }
  return NULL;
}


/* ------------------ compress_volslices ------------------------ */

void *compress_volslices(void *arg){
  int *thread_index;
  int i;

  thread_index = (int *)arg;
  if(nvolrenderinfo<=0)return NULL;

  if(GLOBcleanfiles==1)return NULL;

  // convert and compress files

  for(i=0;i<nsliceinfo;i++){
    slice *slicei;

    slicei = sliceinfo + i;

    if(slicei->isvolslice==0)continue;

    LOCK_VOLSLICE;
    if(slicei->involuse==1){
      UNLOCK_VOLSLICE;
      continue;
    }
    slicei->involuse=1;
    UNLOCK_VOLSLICE;

    convert_volslice(slicei,thread_index);
  }
  return NULL;
}

/* ------------------ compress_slices ------------------------ */

void *compress_slices(void *arg){
  int i;
  slice *slicei, *sb;
  int *thread_index;

  thread_index = (int *)arg;


  if(nsliceinfo<=0)return NULL;
  LOCK_SLICE;
  if(GLOBfirst_slice==1){
    GLOBfirst_slice=0;
    if(GLOBcleanfiles==1){
      for(i=0;i<nsliceinfo;i++){
        slicei = sliceinfo + i;
        convert_slice(slicei,thread_index);
      }
      UNLOCK_SLICE;
      return NULL;
    }
    for(i=0;i<nsliceinfo;i++){
      slicei = sliceinfo + i;
      if(GLOBautozip==1&&slicei->autozip==0)continue;
      slicei->count=0;
    }
    if(GLOBget_slice_bounds==1){
      Get_Slice_Bounds();
    }
    for(i=0;i<nsliceinfo;i++){
      char *label;

      slicei = sliceinfo + i;
      if(GLOBautozip==1&&slicei->autozip==0)continue;
      slicei->doit=1;

      sb=getslice(slicei->label.shortlabel);
      if(sb==NULL)slicei->doit=0;
      label = slicei->label.longlabel;
      if(GLOBmake_demo==1&&(strcmp(label,"TEMPERATURE")==0||strcmp(label,"oxygen")==0)){
        slicei->setvalmax=1;
        slicei->setvalmin=1;
        if(strcmp(label,"TEMPERATURE")==0){
          slicei->valmax=620.0;
          slicei->valmin=20.0;
        }
        else{
          slicei->valmax=0.23;
          slicei->valmin=0.0;
        }
      }
      else{
        if(sb!=NULL){
          if(sb->setvalmax!=1||sb->setvalmin!=1)slicei->doit=0;
        }
        slicei->setvalmax=sb->setvalmax;
        slicei->setvalmin=sb->setvalmin;
        slicei->valmax=sb->valmax;
        slicei->valmin=sb->valmin;
      }
      sb->count++;
    }
  }
  UNLOCK_SLICE;

  if(GLOBcleanfiles==1)return NULL;

  // convert and compress files

  for(i=0;i<nsliceinfo;i++){
    slicei = sliceinfo + i;
    if(GLOBautozip==1&&slicei->autozip==0)continue;
    LOCK_SLICE;
    if(slicei->inuse==1){
      UNLOCK_SLICE;
      continue;
    }
    slicei->inuse=1;
    UNLOCK_SLICE;

    if(slicei->doit==1){
      convert_slice(slicei,thread_index);
    }
    else{
      PRINTF("%s not compressed\n",slicei->file);
      PRINTF("  Min and Max for %s not set in .ini file\n",slicei->label.shortlabel);
    }
  }
  return NULL;
}

/* ------------------ slicedup ------------------------ */

int slicedup(slice *slicej, int islice){
  int i;
  slice *slicei;

  for(i=0;i<islice;i++){
    slicei = sliceinfo + i;
    if(slicei->dup==1)continue;
    if(strcmp(slicei->label.shortlabel,slicej->label.shortlabel)==0){
      slicej->dup=1;
      return 1;
    }
  }
  return 0;
}

/* ------------------ update_slice_hist ------------------------ */

void update_slice_hist(void){
  int i;

  for(i=0;i<nsliceinfo;i++){
    slice *slicei;
    int unit1;
    FILE_SIZE lenfile;
    int error1;
    float slicetime1, *sliceframe;
    int sliceframesize;
    int is1, is2, js1, js2, ks1, ks2;
    int testslice;

    slicei = sliceinfo + i;

    LOCK_SLICE_BOUND;
    if(slicei->inuse_getbounds==1){
      UNLOCK_SLICE_BOUND;
      continue;
    }
    slicei->inuse_getbounds=1;
    UNLOCK_SLICE_BOUND;
    PRINTF("  Examining %s\n",slicei->file);

    lenfile=strlen(slicei->file);

    LOCK_COMPRESS;
    FORTget_file_unit(&unit1,&slicei->unit_start);
    FORTopenslice(slicei->file,&unit1,&is1,&is2,&js1,&js2,&ks1,&ks2,&error1,lenfile);
    UNLOCK_COMPRESS;

    sliceframesize=(is2+1-is1)*(js2+1-js1)*(ks2+1-ks1);
    NewMemory((void **)&sliceframe,sliceframesize*sizeof(float));
    ResetHistogram(slicei->histogram);
    testslice=0;
    while(error1==0){
      FORTgetsliceframe(&unit1, &is1, &is2, &js1, &js2, &ks1, &ks2, &slicetime1, sliceframe, &testslice,&error1);
      UpdateHistogram(sliceframe,sliceframesize,slicei->histogram);
    }
    FREEMEMORY(sliceframe);

    LOCK_COMPRESS;
    FORTclosefortranfile(&unit1);
    UNLOCK_COMPRESS;
  }
}
#ifdef pp_THREAD
/* ------------------ MT_update_slice_hist ------------------------ */

void *MT_update_slice_hist(void *arg){
  update_slice_hist();
  return NULL;
}

/* ------------------ mt_update_slice_hist ------------------------ */

void mt_update_slice_hist(void){
  pthread_t *thread_ids;
  int i;

  NewMemory((void **)&thread_ids,mt_nthreads*sizeof(pthread_t));

  for(i=0;i<mt_nthreads;i++){
    pthread_create(&thread_ids[i],NULL,MT_update_slice_hist,NULL);
  }

  for(i=0;i<mt_nthreads;i++){
    pthread_join(thread_ids[i],NULL);
  }
  FREEMEMORY(thread_ids);
}
#endif

/* ------------------ Get_Slice_Bounds ------------------------ */

void Get_Slice_Bounds(void){
  int i;

  int endiandata;

  endiandata=getendian();
  if(endianswitch==1)endiandata=1-endiandata;

  PRINTF("Determining slice file bounds\n");
  for(i=0;i<nsliceinfo;i++){
    slice *slicei;

    slicei = sliceinfo + i;
    slicei->inuse_getbounds=0;
  }
#ifdef pp_THREAD
  mt_update_slice_hist();
#else
  update_slice_hist();
#endif
  for(i=0;i<nsliceinfo;i++){
    slice *slicei;
    int j;

    slicei = sliceinfo + i;
    if(slicei->dup==1)continue;
    for(j=i+1;j<nsliceinfo;j++){
      slice *slicej;

      slicej = sliceinfo + j;
      if(strcmp(slicei->label.shortlabel,slicej->label.shortlabel)!=0)continue;
      MergeHistogram(slicei->histogram,slicej->histogram);
    }
    slicei->valmax=GetHistogramVal(slicei->histogram,0.99);
    slicei->valmin=GetHistogramVal(slicei->histogram,0.01);
    slicei->setvalmax=1;
    slicei->setvalmin=1;
    for(j=i+1;j<nsliceinfo;j++){
      slice *slicej;

      slicej = sliceinfo + j;
      if(strcmp(slicei->label.shortlabel,slicej->label.shortlabel)!=0)continue;
      slicej->valmax=slicei->valmax;
      slicej->valmin=slicei->valmin;
      slicej->setvalmax=1;
      slicej->setvalmin=1;
    }
  }
  for(i=0;i<nsliceinfo;i++){
    slice *slicei;

    slicei = sliceinfo + i;
    FREEMEMORY(slicei->histogram);
  }

}

/* ------------------ getsliceparms_c ------------------------ */

void getsliceparms_c(char *file, int *ni, int *nj, int *nk){
    int skip,ijkbar[6];
    FILE *stream;

    *ni=0;
    *nj=0;
    *nk=0;

    stream=fopen(file,"rb");
    if(stream==NULL)return;

    skip = 3*(4+30+4);  // skip over 3 records each containing a 30 byte FORTRAN character string
    FSEEK(stream,skip,SEEK_CUR);

    skip=4;
    FSEEK(stream,skip,SEEK_CUR);
    fread(ijkbar,sizeof(int),6,stream);
    *ni=ijkbar[1]+1-ijkbar[0];
    *nj=ijkbar[3]+1-ijkbar[2];
    *nk=ijkbar[5]+1-ijkbar[4];
    fclose(stream);
}
