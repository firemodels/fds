// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include "zlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "egz_stdio.h"
#include "svzip.h"
#include "MALLOC.h"

// svn revision character string
char CNVslice_revision[]="$Revision$";

void mt_update_slice_hist(void);

#ifdef pp_RLETEST

#ifndef pp_cvf
#define STDCALL extern void
#else
#define STDCALL extern void _stdcall
#endif

#ifndef pp_noappend
#define FORTrleheader rleheader_
#define FORTrleframe rleframe_
#else
#define FORTrleheader rleheader
#define FORTrleframe rleframe
#endif

// SUBROUTINE RLEHEADER(RLEFILE, MINMAX, IJKMINMAX)
STDCALL FORTrleheader(char *rlefile, float *minmax, int *ijkminmax, int lenrlefile);

// SUBROUTINE RLEFRAME(RLEFILE,T,CVAL,NCVAL)
STDCALL FORTrleframe(char *rlefile, float *t, float *q, int *nq, float *minmax, unsigned char *cval, int *ncval, int lenrlefile, int lencval);
#endif



void endian_switch(void *val, int nval);

#define FORTSLICEREAD(var,size) fseek(SLICEFILE,4,SEEK_CUR);\
                           returncode=fread(var,4,size,SLICEFILE);\
                           if(endianswitch==1)endian_switch(var,size);\
                           fseek(SLICEFILE,4,SEEK_CUR)

#define FORTRLESLICEREAD(var,size) fseek(SLICEFILE,4,SEEK_CUR);\
                           returncode=fread(var,4,size,SLICEFILE);\
                           if(endian_rle_switch==1)endian_switch(var,size);\
                           fseek(SLICEFILE,4,SEEK_CUR)

/* ------------------ convert_slice ------------------------ */

// unsigned int irle(unsigned char *buffer_in, int nchars_in, unsigned char *buffer_out)

int convert_slice(slice *slicei, int *thread_index){

  char slicefile_svz[1024], slicesizefile_svz[1024];
#ifdef pp_RLETEST
  char slicefile_rle[1024];
#endif
  int fileversion, one, zero;
  char *slice_file;
  int version;
  char filetype[1024];
  char *shortlabel, *unit;
  char units[256];
  int ijkbar[6];
  uLong framesize;
  float *sliceframe_data=NULL;
#ifdef pp_RLETEST
  unsigned char *sliceframe_rle=NULL;
  float *sliceframe_float_rle=NULL;
#endif
  unsigned char *sliceframe_compressed=NULL, *sliceframe_uncompressed=NULL;
  unsigned char *sliceframe_compressed_rle=NULL, *sliceframe_uncompressed_rle=NULL;
  char pp[2];
  char xxx[2];
  char cval[256];
  int sizebefore, sizeafter;
  int returncode;
  float minmax[2];
  float time;
  long data_loc;
  int percent_done;
  int percent_next=10;
  float valmin, valmax, denom;
  int chop_min, chop_max;
  uLongf ncompressed_zlib;
  int ncompressed_save;
  int count=0;
  int ncol, nrow, idir;
  int endian_rle_switch;
  float time_max;
  int itime;

  FILE *SLICEFILE;
  FILE *slicestream,*slicesizestream;

#ifdef pp_THREAD
  if(cleanfiles==0){
    int fileindex;

    fileindex = slicei + 1 - sliceinfo;
    sprintf(threadinfo[*thread_index].label,"sf %i",fileindex);
  }
#endif

  slice_file=slicei->file;
  version=slicei->version;

  strcpy(pp,"%");
  strcpy(xxx,"X");

  fileversion = 1;
  one = 1;
  zero=0;

  // check if slice file is accessible

  strcpy(filetype,"");
  shortlabel=slicei->label.shortlabel;
  if(strlen(shortlabel)>0)strcat(filetype,shortlabel);
  trim(filetype);

  if(getfileinfo(slice_file,NULL,NULL)!=0){
    printf("  %s does not exist\n",slice_file);
    return 0;
  }

  SLICEFILE=fopen(slice_file,"rb");
  if(SLICEFILE==NULL){
    printf("  %s could not be opened\n",slice_file);
    return 0;
  }

  // set up slice compressed file

  if(destdir!=NULL){
    strcpy(slicefile_svz,destdir);
    strcat(slicefile_svz,slicei->filebase);
  }
  else{
    strcpy(slicefile_svz,slicei->file);
  }
#ifdef pp_RLETEST
    strcpy(slicefile_rle,slicefile_svz);
    strcat(slicefile_rle,"test.rle");
#endif
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

  if(destdir!=NULL){
    strcpy(slicesizefile_svz,destdir);
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

  if(cleanfiles==1){
    slicestream=fopen(slicefile_svz,"rb");
    if(slicestream!=NULL){
      fclose(slicestream);
      printf("  Removing %s.\n",slicefile_svz);
      UNLINK(slicefile_svz);
      LOCK_COMPRESS;
      filesremoved++;
      UNLOCK_COMPRESS;
    }
    slicesizestream=fopen(slicesizefile_svz,"rb");
    if(slicesizestream!=NULL){
      fclose(slicesizestream);
      printf("  Removing %s.\n",slicesizefile_svz);
      UNLINK(slicesizefile_svz);
      LOCK_COMPRESS;
      filesremoved++;
      UNLOCK_COMPRESS;
    }
    return 0;
  }

  if(overwrite_slice==0){
    slicestream=fopen(slicefile_svz,"rb");
    if(slicestream!=NULL){
      fclose(slicestream);
      printf("  %s exists.\n",slicefile_svz);
      printf("     Use the -f option to overwrite smokezip compressed files\n");
      return 0;
    }
  }

  slicestream=fopen(slicefile_svz,"wb");
  slicesizestream=fopen(slicesizefile_svz,"w");
  if(slicestream==NULL||slicesizestream==NULL){
    if(slicestream==NULL){
      printf("  %s could not be opened for writing\n",slicefile_svz);
    }
    if(slicesizestream==NULL){
      printf("  %s could not be opened for writing\n",slicesizefile_svz);
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
  trim(units);
  sprintf(cval,"%f",slicei->valmin);
  trimzeros(cval);
#ifndef pp_THREAD
  if(cleanfiles==0){
    printf("Compressing slice file (%s) %s\n",filetype,slice_file);
    printf("    using min=%s %s",cval,units);
  }
#endif
  sprintf(cval,"%f",slicei->valmax);
  trimzeros(cval);
#ifndef pp_THREAD
  if(cleanfiles==0){
    printf(" max=%s %s\n\n",cval,units);
  }
#endif
  valmin=slicei->valmin;
  valmax=slicei->valmax;
  denom = valmax-valmin;
  if(denom==0.0)denom=1.0;
  
  chop_min=0;
  chop_max=255;
  if(slicei->rle==0&&no_chop==0){
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
  fwrite(&version,4,1,slicestream);       // fds slice file version
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
  // version  (slicef version)
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

  if(slicei->rle==0){
    int skip;

    skip = 3*(4+30+4);  // skip over 3 records each containing a 30 byte FORTRAN character string
    returncode=fseek(SLICEFILE,skip,SEEK_CUR);
    sizebefore=skip;
  }

  if(slicei->rle==0){
    FORTSLICEREAD(ijkbar,6);
    sizebefore+=8+6*4;
  }
  else{
    int one;

    fseek(SLICEFILE,4,SEEK_CUR);fread(&one,4,1,SLICEFILE);fseek(SLICEFILE,4,SEEK_CUR);
    sizebefore = 12;
    
    endian_rle_switch=0;                                
    if(one!=1)endian_rle_switch=1;

    fseek(SLICEFILE,4*4,SEEK_CUR);
    sizebefore += 16;

    FORTRLESLICEREAD(minmax,2);
    sizebefore += 8 + 2*4;

    FORTRLESLICEREAD(ijkbar,6);
    sizebefore += 8 + 6*4;

  }

  framesize =  (ijkbar[1]+1-ijkbar[0]);
  framesize *= (ijkbar[3]+1-ijkbar[2]);
  framesize *= (ijkbar[5]+1-ijkbar[4]);

  if(slicei->rle==0){
    minmax[0]=slicei->valmin;
    minmax[1]=slicei->valmax;
  }
  fwrite(minmax,4,2,slicestream);    // min max vals
  fwrite(ijkbar,4,6,slicestream);
  sizeafter+=(8+24);

#ifdef pp_RLETEST
  {
    int lenslicefile_rle;

    lenslicefile_rle=strlen(slicefile_rle);
    FORTrleheader(slicefile_rle,minmax,ijkbar,lenslicefile_rle);
  }
#endif

  ncompressed_save=1.02*framesize+600;
#ifdef pp_RLETEST
  if(NewMemory((void **)&sliceframe_float_rle,ncompressed_save*sizeof(float))==0)goto wrapup;
#endif
  if(NewMemory((void **)&sliceframe_data,ncompressed_save*sizeof(float))==0)goto wrapup;
  if(NewMemory((void **)&sliceframe_compressed,ncompressed_save*sizeof(unsigned char))==0)goto wrapup;
  if(slicei->rle==1){
    if(NewMemory((void **)&sliceframe_compressed_rle,ncompressed_save*sizeof(unsigned char))==0)goto wrapup;
    if(NewMemory((void **)&sliceframe_uncompressed_rle,ncompressed_save*sizeof(unsigned char))==0)goto wrapup;
  }
  if(NewMemory((void **)&sliceframe_uncompressed,ncompressed_save*sizeof(unsigned char))==0)goto wrapup;
#ifdef pp_RLETEST
  if(NewMemory((void **)&sliceframe_rle,framesize*sizeof(unsigned char))==0)goto wrapup;
#endif

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
      int ncomp_rle;

      if(slicei->rle==0){
        FORTSLICEREAD(&time,1);
      }
      else{
        FORTRLESLICEREAD(&time,1);
      }
      sizebefore+=12;
      if(returncode==0)break;
      if(slicei->rle==0){
        FORTSLICEREAD(sliceframe_data,framesize);    //---------------
        if(returncode==0)break;

        sizebefore+=(8+framesize*4);
      }
      else{
        FORTSLICEREAD(&ncomp_rle,1);
        sizebefore+=12;

        returncode=fseek(SLICEFILE,4,SEEK_CUR);
        if(returncode!=0)break;

        returncode=fread(sliceframe_compressed_rle,1,ncomp_rle,SLICEFILE);
        if(returncode==0)break;

        returncode=fseek(SLICEFILE,4,SEEK_CUR);
        if(returncode!=0)break;
        sizebefore+=(8+ncomp_rle);
        
        irle(sliceframe_compressed_rle,ncomp_rle,sliceframe_uncompressed_rle);


      }
      if(time<time_max)continue;
      time_max=time;

      count++;
   
      data_loc=ftell(SLICEFILE);
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
        printf(" %i%s",percent_next,pp);
        fflush(stdout);
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

        if(slicei->rle==0){
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
        else{
          sliceframe_uncompressed[i] = sliceframe_uncompressed_rle[i];
        }
#ifdef pp_RLETEST
        sliceframe_float_rle[index2] = sliceframe_data[i];
#endif
      }
      itime++;
      if(itime%slicezipstep!=0)continue;

      //int compress (Bytef *dest,   uLongf *destLen, const Bytef *source, uLong sourceLen);
      ncompressed_zlib=ncompressed_save;
      returncode=compress(sliceframe_compressed,&ncompressed_zlib,sliceframe_uncompressed,framesize);

      fwrite(&time,4,1,slicestream);
      fwrite(&ncompressed_zlib,4,1,slicestream);
      fwrite(sliceframe_compressed,1,ncompressed_zlib,slicestream);
      sizeafter+=(8+ncompressed_zlib);
#ifdef pp_RLETEST
      {
        int lenslicefile_rle;
        lenslicefile_rle = strlen(slicefile_rle);
        ncompressed_rle=framesize;
        FORTrleframe(slicefile_rle,&time,sliceframe_float_rle,&framesize, minmax, sliceframe_rle,&ncompressed_rle,lenslicefile_rle,ncompressed_rle);
        fprintf(slicesizestream,"%f %i %i\n",time,ncompressed_zlib,nrle);
      }
#else
      if(slicei->rle==0){
        fprintf(slicesizestream,"%f %i\n",time,(int)ncompressed_zlib);
      }
      else{
        fprintf(slicesizestream,"%f %i %i\n",time,(int)ncompressed_zlib,(int)ncomp_rle);
      }
#endif

    }
    if(returncode!=0){
      printf("*** error: compress returncode=%i\n",returncode);
    }
  }

wrapup:
#ifndef pp_THREAD
    printf(" 100%s completed\n",pp);
#endif
    FREEMEMORY(sliceframe_data);
    FREEMEMORY(sliceframe_compressed);
    FREEMEMORY(sliceframe_uncompressed);
    FREEMEMORY(sliceframe_compressed_rle);
    FREEMEMORY(sliceframe_uncompressed_rle);
#ifdef pp_RLETEST
    FREEMEMORY(sliceframe_float_rle);
    FREEMEMORY(sliceframe_rle);
#endif

  fclose(SLICEFILE);
  fseek(slicestream,4,SEEK_SET);
  fwrite(&one,4,1,slicestream);  // write completion code
  fclose(slicestream);
  fclose(slicesizestream);

  {
    char before_label[256],after_label[256];

    getfilesizelabel(sizebefore,before_label);
    getfilesizelabel(sizeafter,after_label);
#ifdef pp_THREAD
    slicei->compressed=1;
    sprintf(slicei->summary,"compressed from %s to %s (%4.1f%s reduction)",before_label,after_label,(float)sizebefore/(float)sizeafter,xxx);
    threadinfo[*thread_index].stat=-1;
#else
    printf("    records=%i, ",count);
    printf("Sizes: original=%s, ",before_label);
    printf("compressed=%s (%4.1f%s reduction)\n",after_label,(float)sizebefore/(float)sizeafter,xxx);
#endif
  }

  return 1;

}

/* ------------------ getslice ------------------------ */

slice *getslice(char *string){
  int i;
  slice *slicei;

  for(i=0;i<nslice_files;i++){
    slicei = sliceinfo + i;
    if(slicei->dup==1)continue;
    if(strcmp(slicei->label.shortlabel,string)==0)return slicei;
  }
  return NULL;
}

/* ------------------ compress_slices ------------------------ */

void *compress_slices(void *arg){
  int i;
  slice *slicei, *sb;
  int *thread_index;

  thread_index = (int *)arg;


  if(nslice_files<=0)return NULL;
  LOCK_SLICE;
  if(first_slice==1){
    first_slice=0;
    if(cleanfiles==1){
      for(i=0;i<nslice_files;i++){
        slicei = sliceinfo + i;
        convert_slice(slicei,thread_index);
      }
      return NULL;
    }
    for(i=0;i<nslice_files;i++){
      slicei = sliceinfo + i;
      if(autozip==1&&slicei->autozip==0)continue;
      slicei->count=0;
    }
    if(get_slice_bounds==1){
      Get_Slice_Bounds();
    }
    for(i=0;i<nslice_files;i++){
      char *label;

      slicei = sliceinfo + i;
      if(autozip==1&&slicei->autozip==0)continue;
      slicei->doit=1;

      sb=getslice(slicei->label.shortlabel);
      if(sb==NULL)slicei->doit=0;
      label = slicei->label.longlabel;
      if(make_demo==1&&(strcmp(label,"TEMPERATURE")==0||strcmp(label,"oxygen")==0)){
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

  // convert and compress files

  for(i=0;i<nslice_files;i++){
    slicei = sliceinfo + i;
    if(autozip==1&&slicei->autozip==0)continue;
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
      printf("%s not compressed\n",slicei->file);
      printf("  Min and Max for %s not set in .ini file\n",slicei->label.shortlabel);
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

  int endiandata;

  endiandata=getendian();
  for(i=0;i<nslice_files;i++){
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

    lenfile=strlen(slicei->file);
    
    LOCK_COMPRESS;
    FORTget_file_unit(&unit1,&slicei->unit_start);
    FORTopenslice(slicei->file,&unit1,&endiandata,&is1,&is2,&js1,&js2,&ks1,&ks2,&error1,lenfile);
    UNLOCK_COMPRESS;

    sliceframesize=(is2+1-is1)*(js2+1-js1)*(ks2+1-ks1);
    NewMemory((void **)&sliceframe,sliceframesize*sizeof(float));
    init_histogram(slicei->histogram);
    testslice=0;
    while(error1==0){
      FORTgetsliceframe(&unit1, &is1, &is2, &js1, &js2, &ks1, &ks2, &slicetime1, sliceframe, &testslice,&error1);
      update_histogram(sliceframe,sliceframesize,slicei->histogram);
    }
    FREEMEMORY(sliceframe);
    
    LOCK_COMPRESS;
    FORTclosefortranfile(&unit1);
    UNLOCK_COMPRESS;
  }
}

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

/* ------------------ Get_Slice_Bounds ------------------------ */

void Get_Slice_Bounds(void){
  int i;

  int endiandata;

  endiandata=getendian();
  if(endianswitch==1)endiandata=1-endiandata;

#ifndef pp_THREAD
  printf("Determining slice file bounds\n");
#endif
  for(i=0;i<nslice_files;i++){
    slice *slicei;

    slicei = sliceinfo + i;
    slicei->inuse_getbounds=0;
  }
#ifdef pp_THREAD
  if(mt_compress==1){
    mt_update_slice_hist();
  }
  else{
    update_slice_hist();
  }
#else
  update_slice_hist();
#endif
  for(i=0;i<nslice_files;i++){
    slice *slicei;
    int j;

    slicei = sliceinfo + i;
    if(slicei->dup==1)continue;
    for(j=i+1;j<nslice_files;j++){
      slice *slicej;

      slicej = sliceinfo + j;
      if(strcmp(slicei->label.shortlabel,slicej->label.shortlabel)!=0)continue;
      merge_histogram(slicei->histogram,slicej->histogram);
    }
    slicei->valmax=get_histogram_value(slicei->histogram,0.99);
    slicei->valmin=get_histogram_value(slicei->histogram,0.01);
    slicei->setvalmax=1;
    slicei->setvalmin=1;
    for(j=i+1;j<nslice_files;j++){
      slice *slicej;

      slicej = sliceinfo + j;
      if(strcmp(slicei->label.shortlabel,slicej->label.shortlabel)!=0)continue;
      slicej->valmax=slicei->valmax;
      slicej->valmin=slicei->valmin;
      slicej->setvalmax=1;
      slicej->setvalmin=1;
    }
  }
  for(i=0;i<nslice_files;i++){
    slice *slicei;

    slicei = sliceinfo + i;
    FREEMEMORY(slicei->histogram);
  }

}
