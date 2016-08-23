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

pdfdata pdfmerge,pdfframe;

#define FORTREAD(var,size) FSEEK(BOUNDARYFILE,4,SEEK_CUR);\
                           returncode=fread(var,4,size,BOUNDARYFILE);\
                           if(endianswitch==1)endian_switch(var,size);\
                           FSEEK(BOUNDARYFILE,4,SEEK_CUR)

/* ------------------ clean_boundary ------------------------ */

int clean_boundary(patch *patchi){
  FILE *BOUNDARYFILE=NULL;
  char boundaryfile_svz[1024], boundarysizefile_svz[1024];
  FILE *boundarystream=NULL,*boundarysizestream=NULL;
  char *boundary_file;
  char filetype[256];
  char *shortlabel;

  boundary_file=patchi->file;

  // check if boundary file is accessible

  strcpy(filetype,"");
  shortlabel=patchi->label.shortlabel;
  if(strlen(shortlabel)>0)strcat(filetype,shortlabel);
  TrimBack(filetype);

  BOUNDARYFILE=fopen(boundary_file,"rb");
  if(BOUNDARYFILE==NULL){
    return 0;
  }
  fclose(BOUNDARYFILE);

  // set up boundary compressed file

  if(GLOBdestdir!=NULL){
    strcpy(boundaryfile_svz,GLOBdestdir);
    strcat(boundaryfile_svz,patchi->filebase);
  }
  else{
    strcpy(boundaryfile_svz,patchi->file);
  }
  strcat(boundaryfile_svz,".svz");

  if(GLOBdestdir!=NULL){
    strcpy(boundarysizefile_svz,GLOBdestdir);
    strcat(boundarysizefile_svz,patchi->filebase);
  }
  else{
    strcpy(boundarysizefile_svz,patchi->file);
  }
  strcat(boundarysizefile_svz,".szz");

  boundarystream=fopen(boundaryfile_svz,"rb");
  if(boundarystream!=NULL){
    fclose(boundarystream);
    PRINTF("  Removing %s\n",boundaryfile_svz);
    UNLINK(boundaryfile_svz);
    LOCK_COMPRESS;
    GLOBfilesremoved++;
    UNLOCK_COMPRESS;
  }
  boundarysizestream=fopen(boundarysizefile_svz,"rb");
  if(boundarysizestream!=NULL){
    fclose(boundarysizestream);
    PRINTF("  Removing %s\n",boundarysizefile_svz);
    UNLINK(boundarysizefile_svz);
    LOCK_COMPRESS;
    GLOBfilesremoved++;
    UNLOCK_COMPRESS;
  }
  return 0;
}

/* ------------------ convert_boundary ------------------------ */

int convert_boundary(patch *patchi, int *thread_index){
  FILE *BOUNDARYFILE=NULL;
  char boundaryfile_svz[1024], boundarysizefile_svz[1024];
  FILE *boundarystream=NULL,*boundarysizestream=NULL;
  int npatch;
  int ijkbounds[9];
  int i;
  int fileversion,one;
  float time_local;
  int i1, i2, j1, j2, k1, k2;
  float *patchvals=NULL,*patchvalscopy;
  unsigned char *full_boundarybuffer=NULL,*compressed_boundarybuffer=NULL;
  int returncode;
  int ncompressed_zlibSAVE;
  uLongf ncompressed_zlib;
  uLong npatchfull;
  unsigned int sizebefore=0, sizeafter=0;
  int count=-1;
  int version_local;
  char *boundary_file;
  char filetype[256];
  char *shortlabel, *unit;
  char units[256];
  float minmax[2];
  char cval[256];
  int percent_done;
  int percent_next=10;
  LINT data_loc;
  int zero=0;
  float time_max;

  boundary_file=patchi->file;
  version_local=patchi->version;
  patchi->compressed=0;

#ifdef pp_THREAD
  {
    int fileindex;

    fileindex = patchi + 1 - patchinfo;
    sprintf(threadinfo[*thread_index].label,"bf %i",fileindex);
  }
#endif
  fileversion = 1;
  one = 1;
  zero=0;

  // check if boundary file is accessible

  strcpy(filetype,"");
  shortlabel=patchi->label.shortlabel;
  if(strlen(shortlabel)>0)strcat(filetype,shortlabel);
  TrimBack(filetype);

  if(getfileinfo(boundary_file,NULL,NULL)!=0){
    fprintf(stderr,"*** Warning: The file %s does not exist\n",boundary_file);
    return 0;
  }

  BOUNDARYFILE=fopen(boundary_file,"rb");
  if(BOUNDARYFILE==NULL){
    fprintf(stderr,"*** Warning: The file %s could not be opened\n",boundary_file);
    return 0;
  }

  // set up boundary compressed file

  if(GLOBdestdir!=NULL){
    strcpy(boundaryfile_svz,GLOBdestdir);
    strcat(boundaryfile_svz,patchi->filebase);
  }
  else{
    strcpy(boundaryfile_svz,patchi->file);
  }
  strcat(boundaryfile_svz,".svz");

  if(GLOBdestdir!=NULL){
    strcpy(boundarysizefile_svz,GLOBdestdir);
    strcat(boundarysizefile_svz,patchi->filebase);
  }
  else{
    strcpy(boundarysizefile_svz,patchi->file);
  }
  strcat(boundarysizefile_svz,".szz");

  if(GLOBoverwrite_b==0){
    boundarystream=fopen(boundaryfile_svz,"rb");
    boundarysizestream=fopen(boundarysizefile_svz,"r");
    if(boundarystream!=NULL||boundarysizestream!=NULL){
      if(boundarystream!=NULL){
        fclose(boundarystream);
        fprintf(stderr,"*** Warning: The file %s exists.\n",boundaryfile_svz);
        fprintf(stderr,"     Use the -f option to overwrite smokezip compressed files\n");
      }
      fclose(BOUNDARYFILE);
      return 0;
    }
  }

  boundarystream=fopen(boundaryfile_svz,"wb");
  boundarysizestream=fopen(boundarysizefile_svz,"w");
  if(boundarystream==NULL||boundarysizestream==NULL){
    if(boundarystream==NULL){
      fprintf(stderr,"*** Warning: The file %s could not be opened for writing\n",boundaryfile_svz);
    }
    if(boundarysizestream==NULL){
      fprintf(stderr,"*** Warning: The file %s could not be opened for writing\n",boundarysizefile_svz);
    }
    if(boundarystream!=NULL)fclose(boundarystream);
    if(boundarysizestream!=NULL)fclose(boundarysizestream);
    fclose(BOUNDARYFILE);
    return 0;
  }


  // read and write boundary header
#ifndef pp_THREAD
  PRINTF("Compressing %s (%s)\n",boundary_file,filetype);
#endif

  strcpy(units,"");
  unit=patchi->label.unit;
  if(strlen(unit)>0)strcat(units,unit);
  TrimBack(units);
  sprintf(cval,"%f",patchi->valmin);
  TrimZeros(cval);
#ifndef pp_THREAD
  PRINTF("  using min=%s %s",cval,units);
#endif
  sprintf(cval,"%f",patchi->valmax);
  TrimZeros(cval);
#ifndef pp_THREAD
  PRINTF(" max=%s %s\n",cval,units);
#endif


  fwrite(&one,4,1,boundarystream);           // write out a 1 to determine "endianness" when file is read in later
  fwrite(&zero,4,1,boundarystream);          // write out a zero now, then a one just before file is closed
  fwrite(&fileversion,4,1,boundarystream);   // write out compressed fileversion in case file format changes later
  fwrite(&version_local,4,1,boundarystream);       // fds boundary file version
  sizeafter=16;

  // endian
  // completion (0/1)
  // fileversion (compressed format)
  // version  (bndf version)
  // global min max (used to perform conversion)
  // local min max  (min max found for this file)
  // npatch
  // i1,i2,j1,j2,k1,k2,idir,dummy,dummy (npatch times)
  // time_local
  // compressed size of frame
  // compressed buffer


  {
    int skip;

    skip = 3*(4+30+4);  // skip over 3 records each containing a 30 byte FORTRAN character string
    returncode=FSEEK(BOUNDARYFILE,skip,SEEK_CUR);
    sizebefore=skip;
  }

  npatch=0;
  if(returncode==0){
    FORTREAD(&npatch,1);
    if(returncode==0)npatch=0;
    sizebefore+=12;

    minmax[0]=patchi->valmin;
    minmax[1]=patchi->valmax;
    fwrite(minmax,4,2,boundarystream);    // conversion min max vals
    FSEEK(boundarystream,8,SEEK_CUR);       // skip over local min max vals (we're set in pass 1);
    fwrite(&npatch,4,1,boundarystream);   // write out npatch
    sizeafter+=20;
  }

  if(npatch>0){

    int nbounds=6;
    int *ijks=NULL,*ijkscopy;

    if(NewMemory((void **)&ijks,6*npatch*sizeof(int))==0)goto wrapup;
    CheckMemory;
    ijkscopy=ijks;
    if(version_local==1)nbounds=9;

    npatchfull=0;
    for(i=0;i<npatch;i++){
      int j;

      FORTREAD(ijkbounds,nbounds);
      sizebefore+=(nbounds+2)*4;
      if(returncode==0)goto wrapup;
      fwrite(ijkbounds,4,nbounds,boundarystream);         // write out i1,i2,j1,j2,k1,k2,idir,dummy,dummy
      sizeafter+=4*nbounds;                               // note:  data can be read into one block of size 9*nblocks

      i1 = ijkbounds[0];
      i2 = ijkbounds[1];
      j1 = ijkbounds[2];
      j2 = ijkbounds[3];
      k1 = ijkbounds[4];
      k2 = ijkbounds[5];
      for(j=0;j<6;j++){
        *ijkscopy++=ijkbounds[j];
      }
      npatchfull+=(i2+1-i1)*(j2+1-j1)*(k2+1-k1);
      CheckMemory;
    }

    ncompressed_zlibSAVE=1.01*npatchfull+600;
    if(NewMemory((void **)&patchvals,npatchfull*sizeof(float))==0)goto wrapup;
    if(NewMemory((void **)&full_boundarybuffer,npatchfull)==0)goto wrapup;
    if(NewMemory((void **)&compressed_boundarybuffer,ncompressed_zlibSAVE)==0)goto wrapup;
#ifndef pp_THREAD
    PRINTF(" ");
#endif
    time_max=-1000000.0;
    while(feof(BOUNDARYFILE)==0){
      int j ;

      FORTREAD(&time_local,1);
      sizebefore+=12;
      if(returncode==0)break;

      patchvalscopy=patchvals;
      for(j=0;j<npatch;j++){
        int size;

        i1 = ijks[6*j];
        i2 = ijks[6*j+1];
        j1 = ijks[6*j+2];
        j2 = ijks[6*j+3];
        k1 = ijks[6*j+4];
        k2 = ijks[6*j+5];
        size = (i2+1-i1)*(j2+1-j1)*(k2+1-k1);

        FORTREAD(patchvalscopy,size);
        sizebefore+=(size+2)*4;
        if(returncode==0)goto wrapup;
        patchvalscopy+=size;
      }

//      patchi = patchinfo + i;

      if(time_local<time_max)continue;
      count++;

      if(count%GLOBboundzipstep!=0)continue;
      time_max=time_local;

      for(i=0;i<npatchfull;i++){
        unsigned char ival;
        float val;

        val = patchvals[i];

        if(val<patchi->valmin){
          ival=0;
        }
        else if(val>patchi->valmax){
          ival=255;
        }
        else{
          ival=1+253*(val-patchi->valmin)/(patchi->valmax-patchi->valmin);
        }
        full_boundarybuffer[i]=ival;
      }

      //int compress (Bytef *dest,   uLongf *destLen, const Bytef *source, uLong sourceLen);
      ncompressed_zlib=ncompressed_zlibSAVE;
      returncode=compress_zlib(compressed_boundarybuffer, &ncompressed_zlib, full_boundarybuffer, npatchfull);
      if(returncode!=0){
        fprintf(stderr,"*** Error: compress returncode=%i\n",returncode);
      }
//      PRINTF("time=%f before %i after=%i\n",time_local,npatchfull,ncompressed_zlib);

      fprintf(boundarysizestream,"%f %i %i\n",time_local,(int)npatchfull,(int)ncompressed_zlib);
      fwrite(&time_local,4,1,boundarystream);                                       // write out time_local
      fwrite(&ncompressed_zlib,4,1,boundarystream);                           // write out compressed size of frame
      fwrite(compressed_boundarybuffer,1,ncompressed_zlib,boundarystream);    // write out compressed buffer
      sizeafter+=ncompressed_zlib+8;

      data_loc=FTELL(BOUNDARYFILE);
      percent_done=100.0*(float)data_loc/(float)patchi->filesize;
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
wrapup:
#ifndef pp_THREAD
    PRINTF(" 100%s completed\n",GLOBpp);
#endif
    FREEMEMORY(ijks);
    FREEMEMORY(patchvals);
    FREEMEMORY(full_boundarybuffer);
    FREEMEMORY(compressed_boundarybuffer);
  }

  fclose(BOUNDARYFILE);
  FSEEK(boundarystream,8,SEEK_SET);
  fwrite(&one,4,1,boundarystream);  // write completion code
  fclose(boundarystream);
  fclose(boundarysizestream);
  {
    char before_label[256],after_label[256];
    getfilesizelabel(sizebefore,before_label);
    getfilesizelabel(sizeafter,after_label);
#ifdef pp_THREAD
    patchi->compressed=1;
    sprintf(patchi->summary,"compressed from %s to %s (%4.1f%s reduction)",before_label,after_label,(float)sizebefore/(float)sizeafter,GLOBx);
    threadinfo[*thread_index].stat=-1;
#else
    PRINTF("  records=%i, ",count);
    PRINTF("Sizes: original=%s, ",before_label);
    PRINTF("compressed=%s (%4.1f%s reduction)\n\n",after_label,(float)sizebefore/(float)sizeafter,GLOBx);
#endif
  }

  return 1;

}

/* ------------------ patchdup ------------------------ */

int patchdup(patch *patchj, int ipatch){
  int i;
  patch *patchi;

  for(i=0;i<ipatch;i++){
    patchi = patchinfo + i;
    if(patchi->dup==1)continue;
    if(strcmp(patchi->label.shortlabel,patchj->label.shortlabel)==0){
      patchj->dup=1;
      return 1;
    }
  }
  return 0;
}

/* ------------------ getpatch ------------------------ */

patch *getpatch(char *string){
  int i;
  patch *patchi;

  for(i=0;i<npatchinfo;i++){
    patchi = patchinfo + i;
    if(patchi->dup==1)continue;
    if(strcmp(patchi->label.shortlabel,string)==0)return patchi;
  }
  return NULL;
}

/* ------------------ compress_patches ------------------------ */

void *compress_patches(void *arg){
  int i;
  patch *patchi;
  patch *pb;
  int *thread_index;

  thread_index = (int *)arg;

  if(npatchinfo<=0)return NULL;
  LOCK_PATCH;
  if(GLOBfirst_patch==1){
    GLOBfirst_patch=0;

    if(GLOBcleanfiles==1){
      for(i=0;i<npatchinfo;i++){
        patchi = patchinfo + i;
        clean_boundary(patchi);
      }
      UNLOCK_PATCH;
      return NULL;
    }
    for(i=0;i<npatchinfo;i++){
      patchi = patchinfo + i;
      if(GLOBautozip==1&&patchi->autozip==0)continue;

      pb=getpatch(patchi->label.shortlabel);
      if(pb!=NULL){
        patchi->setvalmax=pb->setvalmax;
        patchi->setvalmin=pb->setvalmin;
        patchi->valmax=pb->valmax;
        patchi->valmin=pb->valmin;
      }
      else{
        patchi->setvalmax=1;
        patchi->setvalmin=1;
      }
    }

  // find bounds

    if(GLOBget_boundary_bounds==1){
      Get_Boundary_Bounds();
    }
    for(i=0;i<npatchinfo;i++){
      patchi = patchinfo + i;
      if(patchi->setvalmin==1&&patchi->setvalmax==1){
        patchi->doit=1;
      }
      else{
        patchi->doit=0;
      }
    }
  }
  UNLOCK_PATCH;

  if(GLOBcleanfiles==1)return NULL;

  // convert and compress files

  for(i=0;i<npatchinfo;i++){
    patchi = patchinfo + i;
    if(GLOBautozip==1&&patchi->autozip==0)continue;

    if(patchi->doit==1){
      LOCK_PATCH;
      if(patchi->inuse==1){
        UNLOCK_PATCH;
        continue;
      }
      patchi->inuse=1;
      UNLOCK_PATCH;

      convert_boundary(patchi,thread_index);
    }
    else{
      PRINTF("%s not compressed\n",patchi->file);
      PRINTF("  Min and Max for %s not set in .ini file\n",patchi->label.shortlabel);
    }
  }
  return NULL;
}

/* ------------------ update_patch_hist ------------------------ */

void update_patch_hist(void){
  int i;
  int endiandata;

  endiandata=getendian();
  if(endianswitch==1)endiandata=1-endiandata;

  for(i=0;i<npatchinfo;i++){
    patch *patchi;
    int unit1;
    FILE_SIZE lenfile;
    int error1;
    int *pi1, *pi2, *pj1, *pj2, *pk1, *pk2;
    float patchtime1, *patchframe;
    int patchframesize;
    int j;

    patchi = patchinfo + i;
    LOCK_PATCH_BOUND;
    if(patchi->inuse_getbounds==1){
      UNLOCK_PATCH_BOUND;
      continue;
    }
    patchi->inuse_getbounds=1;
    UNLOCK_PATCH_BOUND;

    PRINTF("  Examining %s\n",patchi->file);
    lenfile=strlen(patchi->file);
    pi1 = patchi->pi1;
    pi2 = patchi->pi2;
    pj1 = patchi->pj1;
    pj2 = patchi->pj2;
    pk1 = patchi->pk1;
    pk2 = patchi->pk2;

    LOCK_COMPRESS;
    FORTget_file_unit(&unit1,&patchi->unit_start);
    FORTopenboundary(patchi->file,&unit1,&patchi->version,&error1,lenfile);
    UNLOCK_COMPRESS;

    patchframesize=0;
    for(j=0;j<patchi->npatches;j++){
      patchframesize+=patchi->patchsize[j];
    }
    NewMemory((void **)&patchframe,patchframesize*sizeof(float));
    ResetHistogram(patchi->histogram);
    while(error1==0){
      int ndummy;

      FORTgetpatchdata(&unit1, &patchi->npatches,
        pi1, pi2, pj1, pj2, pk1, pk2, &patchtime1, patchframe, &ndummy,&error1);
      UpdateHistogram(patchframe,patchframesize,patchi->histogram);
    }
    LOCK_COMPRESS;
    FORTclosefortranfile(&unit1);
    UNLOCK_COMPRESS;
    FREEMEMORY(patchframe);
  }
#ifndef pp_THREAD
  PRINTF("\n");
#endif
}

#ifdef pp_THREAD
/* ------------------ MT_update_slice_hist ------------------------ */

void *MT_update_patch_hist(void *arg){
  update_patch_hist();
  return NULL;
}

/* ------------------ mt_update_slice_hist ------------------------ */

void mt_update_patch_hist(void){
  pthread_t *thread_ids;
  int i;

  NewMemory((void **)&thread_ids,mt_nthreads*sizeof(pthread_t));

  for(i=0;i<mt_nthreads;i++){
    pthread_create(&thread_ids[i],NULL,MT_update_patch_hist,NULL);
  }

  for(i=0;i<mt_nthreads;i++){
    pthread_join(thread_ids[i],NULL);
  }
  FREEMEMORY(thread_ids);
}
#endif

/* ------------------ Get_Boundary_Bounds ------------------------ */

void Get_Boundary_Bounds(void){
  int i;

  int endiandata;

  endiandata=getendian();
  if(endianswitch==1)endiandata=1-endiandata;

  PRINTF("Determining boundary file bounds\n");
  for(i=0;i<npatchinfo;i++){
    patch *patchi;

    patchi = patchinfo + i;
    patchi->inuse_getbounds=0;
  }
#ifdef pp_THREAD
  mt_update_patch_hist();
#else
  update_patch_hist();
#endif
  for(i=0;i<npatchinfo;i++){
    patch *patchi;
    int j;

    patchi = patchinfo + i;
    if(patchi->dup==1)continue;
    for(j=i+1;j<npatchinfo;j++){
      patch *patchj;

      patchj = patchinfo + j;
      if(strcmp(patchi->label.shortlabel,patchj->label.shortlabel)!=0)continue;
      MergeHistogram(patchi->histogram,patchj->histogram);
    }
    patchi->valmax=GetHistogramVal(patchi->histogram,0.99);
    patchi->valmin=GetHistogramVal(patchi->histogram,0.01);
    patchi->setvalmax=1;
    patchi->setvalmin=1;
    for(j=i+1;j<npatchinfo;j++){
      patch *patchj;

      patchj = patchinfo + j;
      if(strcmp(patchi->label.shortlabel,patchj->label.shortlabel)!=0)continue;
      patchj->valmax=patchi->valmax;
      patchj->valmin=patchi->valmin;
      patchj->setvalmax=1;
      patchj->setvalmin=1;
    }
  }
  for(i=0;i<npatchinfo;i++){
    patch *patchi;

    patchi = patchinfo + i;
    FREEMEMORY(patchi->histogram->buckets);
    FREEMEMORY(patchi->histogram->buckets_2d);
    FREEMEMORY(patchi->histogram);
  }

}
