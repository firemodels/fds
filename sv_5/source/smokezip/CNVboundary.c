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
char CNVboundary_revision[]="$Revision: 624 $";

pdfdata pdfmerge,pdfframe;


void endian_switch(void *val, int nval);
#define FORTREAD(var,size) fseek(BOUNDARYFILE,4,SEEK_CUR);\
                           returncode=fread(var,4,size,BOUNDARYFILE);\
                           if(endianswitch==1)endian_switch(var,size);\
                           fseek(BOUNDARYFILE,4,SEEK_CUR)

/* ------------------ convert_boundary ------------------------ */

int convert_boundary(patch *patchi,int pass){
  FILE *BOUNDARYFILE=NULL;
  char boundaryfile_svz[1024], boundarysizefile_svz[1024];
  FILE *boundarystream=NULL,*boundarysizestream=NULL;
  int npatch;
  int ijkbounds[9];
  int i;
  int fileversion,one;
  float time;
  int i1, i2, j1, j2, k1, k2;
  float *patchvals=NULL,*patchvalscopy;
  unsigned char *full_boundarybuffer=NULL,*compressed_boundarybuffer=NULL;
  int returncode;
  int ncompressed_zlibSAVE;
  uLongf ncompressed_zlib;
  uLong npatchfull;
  unsigned int sizebefore=0, sizeafter=0;
  int count=-1;
  char pp[2];
  char xxx[2];
  int version;
  char *boundary_file;
  char filetype[256];
  char *shortlabel, *unit;
  char units[256];
  float minmax[2];
  char cval[256];
  int percent_done;
  int percent_next=10;
  long data_loc;
  int zero=0;

  boundary_file=patchi->file;
  version=patchi->version;

  strcpy(pp,"%");
  strcpy(xxx,"X");

  fileversion = 1;
  one = 1;
  zero=0;

  // check if boundary file is accessible

  strcpy(filetype,"");
  shortlabel=patchi->label.shortlabel;
  if(strlen(shortlabel)>0)strcat(filetype,shortlabel);
  trim(filetype);

  if(getfileinfo(boundary_file,NULL,NULL)!=0){
    printf("  %s does not exist\n",boundary_file);
    return 0;
  }

  BOUNDARYFILE=fopen(boundary_file,"rb");
  if(BOUNDARYFILE==NULL){
    printf("  %s could not be opened\n",boundary_file);
    return 0;
  }
  if(pass==1)fclose(BOUNDARYFILE);

  // set up boundary compressed file

  if(destdir!=NULL){
    strcpy(boundaryfile_svz,destdir);
    strcat(boundaryfile_svz,patchi->filebase);
  }
  else{
    strcpy(boundaryfile_svz,patchi->file);
  }
  strcat(boundaryfile_svz,".svz");

  if(destdir!=NULL){
    strcpy(boundarysizefile_svz,destdir);
    strcat(boundarysizefile_svz,patchi->filebase);
  }
  else{
    strcpy(boundarysizefile_svz,patchi->file);
  }
  strcat(boundarysizefile_svz,".szz");

  if(cleanfiles==1){
    boundarystream=fopen(boundaryfile_svz,"rb");
    if(boundarystream!=NULL){
      fclose(boundarystream);
      printf("  Removing %s.\n",boundaryfile_svz);
      unlink(boundaryfile_svz);
      filesremoved++;
    }
    boundarysizestream=fopen(boundarysizefile_svz,"rb");
    if(boundarysizestream!=NULL){
      fclose(boundarysizestream);
      printf("  Removing %s.\n",boundarysizefile_svz);
      unlink(boundarysizefile_svz);
      filesremoved++;
    }
    return 0;
  }

  if(overwrite_b==0){
    boundarystream=fopen(boundaryfile_svz,"rb");
    boundarysizestream=fopen(boundarysizefile_svz,"r");
    if(boundarystream!=NULL||boundarysizestream!=NULL){
      if(boundarystream!=NULL){
        getdatabounds(patchi,boundarystream); // even though file exists we'll retrieve local bound data
        fclose(boundarystream);
        printf("  %s exists.\n",boundaryfile_svz);
        printf("     Use the -f option to overwrite boundary or 3d smoke files\n");
      }
   //   if(boundarysizestream!=NULL){
   //     fclose(boundarysizestream);
   //     printf("  %s exists.\n  Use the -f option if you wish to overwrite it.\n",boundarysizefile_svz);
   //   }
      return 0;
    }
  }

  if(pass==2){
    boundarystream=fopen(boundaryfile_svz,"wb");
    boundarysizestream=fopen(boundarysizefile_svz,"w");
    if(boundarystream==NULL||boundarysizestream==NULL){
      if(boundarystream==NULL){
        printf("  %s could not be opened for writing\n",boundaryfile_svz);
      }
      if(boundarysizestream==NULL){
        printf("  %s could not be opened for writing\n",boundarysizefile_svz);
      }
      if(boundarystream!=NULL)fclose(boundarystream);
      if(boundarysizestream!=NULL)fclose(boundarysizestream);
      fclose(BOUNDARYFILE);
      return 0;
    }
  }
  else{
    returncode = getdatabounds(patchi,NULL);
    strcpy(units,"");
    unit=patchi->label.unit;
    if(strlen(unit)>0)strcat(units,unit);
    trim(units);
    sprintf(cval,"%f",patchi->valmin);
    trimzeros(cval);
    printf("    local bounds min=%s %s",cval,units);
    sprintf(cval,"%f",patchi->valmax);
    trimzeros(cval);
    printf(" max=%s %s\n",cval,units);
    return returncode;
  }

  // read and write boundary header

  if(cleanfiles==0)printf("Compressing boundary file (%s) %s\n",filetype,boundary_file);

  strcpy(units,"");
  unit=patchi->label.unit;
  if(strlen(unit)>0)strcat(units,unit);
  trim(units);
  sprintf(cval,"%f",patchi->valmin);
  trimzeros(cval);
  printf("    using min=%s %s",cval,units);
  sprintf(cval,"%f",patchi->valmax);
  trimzeros(cval);
  printf(" max=%s %s\n",cval,units);


  fwrite(&one,4,1,boundarystream);           // write out a 1 to determine "endianness" when file is read in later
  fwrite(&zero,4,1,boundarystream);          // write out a zero now, then a one just before file is closed
  fwrite(&fileversion,4,1,boundarystream);   // write out compressed fileversion in case file format changes later
  fwrite(&version,4,1,boundarystream);       // fds boundary file version
  sizeafter=16;

  // endian
  // completion (0/1)
  // fileversion (compressed format)
  // version  (bndf version)
  // global min max (used to perform conversion)
  // local min max  (min max found for this file)
  // npatch
  // i1,i2,j1,j2,k1,k2,idir,dummy,dummy (npatch times)
  // time
  // compressed size of frame
  // compressed buffer
  

  {
    int skip;

    skip = 3*(4+30+4);  // skip over 3 records each containing a 30 byte FORTRAN character string
    returncode=fseek(BOUNDARYFILE,skip,SEEK_CUR);
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
    fseek(boundarystream,8,SEEK_CUR);       // skip over local min max vals (we're set in pass 1);
    fwrite(&npatch,4,1,boundarystream);   // write out npatch
    sizeafter+=20;
  }

  if(npatch>0){

    int nbounds=6;
    int *ijks=NULL,*ijkscopy;

    if(NewMemory((void **)&ijks,6*npatch*sizeof(int))==0)goto wrapup;
    CheckMemory;
    ijkscopy=ijks;
    if(version==1)nbounds=9;

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
    printf("  Compressing: ");
    while(feof(BOUNDARYFILE)==0){
      int j ;

      FORTREAD(&time,1);
      sizebefore+=12;
      if(returncode==0)break;

      patchvalscopy=patchvals;
      for(j=0;j<npatch;j++){
        int size;

        i1 = ijks[6*j+0];
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

      count++;

      if(count%boundzipstep!=0)continue;

      for(i=0;i<npatchfull;i++){
        unsigned char ival;
    
        if(patchvals[i]<patchi->valmin)patchvals[i]=patchi->valmin;
        if(patchvals[i]>patchi->valmax)patchvals[i]=patchi->valmax;
        ival=255*(patchvals[i]-patchi->valmin)/(patchi->valmax-patchi->valmin);
        full_boundarybuffer[i]=ival;
      }

      //int compress (Bytef *dest,   uLongf *destLen, const Bytef *source, uLong sourceLen);
      ncompressed_zlib=ncompressed_zlibSAVE;
      returncode=compress(compressed_boundarybuffer, &ncompressed_zlib, full_boundarybuffer, npatchfull);
      if(returncode!=0){
        printf("*** error: compress returncode=%i\n",returncode);
      }
//      printf("time=%f before %i after=%i\n",time,npatchfull,ncompressed_zlib);

      fprintf(boundarysizestream,"%f %i %i\n",time,npatchfull,ncompressed_zlib);
      fwrite(&time,4,1,boundarystream);                                       // write out time
      fwrite(&ncompressed_zlib,4,1,boundarystream);                           // write out compressed size of frame
      fwrite(compressed_boundarybuffer,1,ncompressed_zlib,boundarystream);    // write out compressed buffer
      sizeafter+=ncompressed_zlib+8;

      data_loc=ftell(BOUNDARYFILE);
      percent_done=100.0*(float)data_loc/(float)patchi->filesize;
      if(percent_done>percent_next){
        printf(" %i%s",percent_next,pp);
        fflush(stdout);
        percent_next+=10;
      }
    }
wrapup:
    printf(" 100%s completed\n",pp);
    FREEMEMORY(ijks);
    FREEMEMORY(patchvals);
    FREEMEMORY(full_boundarybuffer);
    FREEMEMORY(compressed_boundarybuffer);
  }

  fclose(BOUNDARYFILE);
  fseek(boundarystream,8,SEEK_SET);
  fwrite(&one,4,1,boundarystream);  // write completion code
  fclose(boundarystream);
  fclose(boundarysizestream);
  {
    char before_label[256],after_label[256];
    getfilesizelabel(sizebefore,before_label);
    getfilesizelabel(sizeafter,after_label);
    printf("    records=%i, ",count);
    printf("Sizes: original=%s, ",before_label);
    printf("compressed=%s (%4.1f%s reduction)\n",after_label,(float)sizebefore/(float)sizeafter,xxx);
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

  for(i=0;i<npatch_files;i++){
    patchi = patchinfo + i;
    if(patchi->dup==1)continue;
    if(strcmp(patchi->label.shortlabel,string)==0)return patchi;
  }
  return NULL;
}

/* ------------------ compress_patches ------------------------ */

void compress_patches(void){
  int i,j;
  patch *patchi,*patchj;
  patch *pb;
  float valmin, valmax;
  int consolidate=0;

  printf("\n");
  for(i=0;i<npatch_files;i++){
    patchi = patchinfo + i;
    patchi->count=0;
  }
  for(i=0;i<npatch_files;i++){
    patchi = patchinfo + i;

    pb=getpatch(patchi->label.shortlabel);
//    if(pb!=NULL&&pb->setvalmin==1&&pb->setvalmax==1&&pb->valmax>pb->valmin){
    if(pb!=NULL){
      patchi->setvalmax=pb->setvalmax;
      patchi->setvalmin=pb->setvalmin;
      patchi->valmax=pb->valmax;
      patchi->valmin=pb->valmin;
      pb->count++;
      if(pb->count>1)consolidate=1;
    }
    else{
      patchi->setvalmax=1;
      patchi->setvalmin=1;
    }
  }

  // find bounds

  if(npatch_files>0&&cleanfiles==0)printf("Finding bounds.\n");
  for(i=0;i<npatch_files;i++){
    patchi = patchinfo + i;
    patchi->done=0;

    patchi->doit=convert_boundary(patchi,1);
  }

  // consolidate bounds

  if(consolidate==1){
    printf("Consolidating bounds.\n");
    for(i=0;i<npatch_files;i++){
      patchi = patchinfo + i;
      if(patchi->doit==0||patchi->done==1)continue;  // if we can't doit or this patch is done then skip it
      if(patchi->setvalmax!=PERCENTILE_MAX&&
         patchi->setvalmin!=PERCENTILE_MIN)continue; // only need to consolidate percentile bounds

    // find bounds

      valmin=patchi->valmin;
      valmax=patchi->valmax;
      for(j=i+1;j<npatch_files;j++){
        patchj = patchinfo + j;
        if(strcmp(patchi->label.shortlabel,patchj->label.shortlabel)!=0)continue;
        patchi->done=1;

        if(patchj->valmin<valmin)valmin=patchj->valmin;
        if(patchj->valmax>valmax)valmax=patchj->valmax;
      }

    // copy bounds back
      smoothlabel(&valmin,&valmax,12);

      for(j=i;j<npatch_files;j++){
        patchj = patchinfo + j;
        if(strcmp(patchi->label.shortlabel,patchj->label.shortlabel)!=0)continue;

        patchj->valmin=valmin;
        patchj->valmax=valmax;
      }
    }
  }

  // convert and compress files

  for(i=0;i<npatch_files;i++){
    patchi = patchinfo + i;

    if(patchi->doit==0)continue;
    convert_boundary(patchi,2);
  }


}

/* ------------------ getdatabounds ------------------------ */

int getdatabounds(patch *patchi,FILE *stream){
  FILE *BOUNDARYFILE;
  int count;
  int i;
  int npatch;
  int npatchfull;
  int *patchsize;
  int returncode;
  float *patchvals;
  float time;
  float *patchvalscopy;
  float pmin, pmax, gmin, gmax;
  int firsttime=1;
  int *buckets;
  float dbucket;
  int percent_done;
  int percent_next=10;
  long data_loc;
  float minmax[2];
  int one;
  int completion;

  // endian
  // completion (0/1)
  // fileversion (compressed format)
  // version  (bndf version)
  // global min max (used to perform conversion)
  // local min max  (min max found for this file)
  // npatch
  // i1,i2,j1,j2,k1,k2,idir,dummy,dummy (npatch times)
  // time
  // compressed size of frame
  // compressed buffer

  if(stream!=NULL){
    int endswitch=0;

    fread(&one,4,1,stream);
    if(one!=1)endswitch=1;
    fread(&completion,4,1,stream);
    if(completion==0)return 0;  // don't need to endian switch a 0
    fseek(stream,16,SEEK_CUR);  // skip over next 4 values
    fread(minmax,4,2,stream);
    if(endswitch==1)endian_switch(minmax,2);
    patchi->valmin=minmax[0];
    patchi->valmax=minmax[1];
    return 1;
  }
  initpdf(&pdfmerge);
  initpdf(&pdfframe);
  if(patchi->setvalmin==1&&patchi->setvalmin==1)return 1;
  BOUNDARYFILE = fopen(patchi->file,"rb");
  if(BOUNDARYFILE==NULL)return 0;
  printf("  for: %s ",patchi->file);

  {
    int skip;

    skip = 3*(4+30+4);  // skip over 3 records each containing a 30 byte FORTRAN character string
    returncode=fseek(BOUNDARYFILE,skip,SEEK_CUR);
  }

  FORTREAD(&npatch,1);

  if(npatch>0){
    int nbounds=6;

    NewMemory((void **)&patchsize,npatch*sizeof(int));
    npatchfull=0;
    if(patchi->version==1)nbounds=9;
    for(i=0;i<npatch;i++){
      int size;
      int ijkbounds[9];
      int i1, i2, j1, j2, k1, k2;

      FORTREAD(ijkbounds,nbounds);
      i1 = ijkbounds[0];
      i2 = ijkbounds[1];
      j1 = ijkbounds[2];
      j2 = ijkbounds[3];
      k1 = ijkbounds[4];
      k2 = ijkbounds[5];
      size = (i2+1-i1)*(j2+1-j1)*(k2+1-k1);
      npatchfull += size;
      patchsize[i]=size;
    }
    patchvals=NULL;
    NewMemory((void **)&patchvals,npatchfull*sizeof(float));

    for(i=1;;i++){
      int j ;

      FORTREAD(&time,1);
      if(returncode==0)break;

      patchvalscopy=patchvals;

      for(j=0;j<npatch;j++){
        FORTREAD(patchvalscopy,patchsize[j]);
        if(returncode==0)goto percentile_loop;
        patchvalscopy+=patchsize[j];
      }
      if(firsttime==1){
        gmin=patchvals[0];
        gmax=gmin;
        firsttime=0;
      }
      data_loc=ftell(BOUNDARYFILE);
      percent_done=100.0*(float)data_loc/(float)patchi->filesize;
      if(percent_done>percent_next){
        printf(" %i%s",percent_next,pp);
        fflush(stdout);
        percent_next+=10;
      }
      getpdf(patchvals,npatchfull,&pdfframe);
      mergepdf(&pdfframe,&pdfmerge,&pdfmerge);
    }
percentile_loop:
    printf(" completed\n");
    buckets=pdfmerge.buckets;
    count=pdfmerge.ncount;
    gmin=pdfmerge.pdfmin;
    gmax=pdfmerge.pdfmax;

    if(buckets!=NULL){
      int alpha05;
      int nsmall=0, nbig=0;
      int total;
      int n;

      alpha05 = (int)(0.01*count);
      total = 0;
      for (n=0;n<NBUCKETS;n++){
        total += buckets[n];
        if(total>alpha05){
          nsmall=n;
          break;
        }
      }
      total = 0;
      for (n=NBUCKETS;n>0;n--){
        total += buckets[n-1];
        if(total>alpha05){
          nbig=n;
          break;
        }
      }
      dbucket=(gmax-gmin)/(float)NBUCKETS;
      pmin = gmin + (nsmall-1)*dbucket;
      pmax = gmin + (nbig+1)*dbucket;
    }
    else{
      pmin=gmin;
      pmax=gmax;
    }
//    FREEMEMORY(buckets);
    smoothlabel(&pmin,&pmax,12);
    smoothlabel(&gmin,&gmax,12);
    switch(patchi->setvalmin){
    case PERCENTILE_MIN:
      patchi->valmin=pmin;
      break;
    case GLOBAL_MIN:
      patchi->valmin=gmax;
      break;
    }
    switch(patchi->setvalmax){
    case PERCENTILE_MAX:
      patchi->valmax=pmax;
      break;
    case GLOBAL_MAX:
      patchi->valmax=gmax;
      break;
    }
    FREEMEMORY(patchvals);
    FREEMEMORY(patchsize);
  }
  fclose(BOUNDARYFILE);
  return 1;
}
