#include "options.h"
// svn revision character string
char CNVpart_revision[]="$Revision: 747 $";
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

/* ------------------ getpatch ------------------------ */

part *getpart(char *string){
  int i;
  part *parti;

  for(i=0;i<npart_files;i++){
    parti = partinfo + i;
    if(parti->dup==1)continue;
    if(strcmp(parti->label.shortlabel,string)==0)return parti;
  }
  return NULL;
}


/* ------------------ patchdup ------------------------ */

int partdup(part *partj, int ipart){
  int i;
  part *parti;

  for(i=0;i<ipart;i++){
    parti = partinfo + i;
    if(parti->dup==1)continue;
    if(strcmp(parti->label.shortlabel,partj->label.shortlabel)==0){
      partj->dup=1;
      return 1;
    }
  }
  return 0;
}

/* ------------------ compress_patches ------------------------ */

void compress_parts(void){
  int i;
  part *parti;
  part *pb;

  for(i=0;i<npart_files;i++){
    parti = partinfo + i;

    pb=getpart(parti->label.shortlabel);
    if(pb!=NULL){
      parti->setvalmax=pb->setvalmax;
      parti->setvalmin=pb->setvalmin;
      parti->valmax=pb->valmax;
      parti->valmin=pb->valmin;
    }
//      getdatabounds(parti);

    convert_part(parti);

  }

}

/* ------------------ compress_part ------------------------ */

void convert_part(part *parti){
  FILE *PARTFILE,*partsizestream;
  EGZ_FILE *partstream;
  char *part_file;
  float valmin, valmax;
  char pp[2];
  char part_file_svz[256], partsizefile_svz[256];
  char filetype[256];
  char *shortlabel, *unit;
  char units[256];
  int idummy[7];
  float rdummy[3];
  int returncode;
  int mxframepts;
  int ibar, jbar, kbar;
  int nb, nv, nspr;
  int i;
  int *ispr;
  int np;
  float time;
  float *xp, *yp, *zp, *bp;
  float *xsp, *ysp, *zsp, *bsp;
  unsigned short *ixyzbp;
  unsigned short *ixp, *iyp, *izp, *ibp;
  unsigned short *ixsp, *iysp, *izsp, *ibsp;
  int nsp;
  float xmin,xmax,ymin,ymax,zmin,zmax;
  float xyzbounds[6];
  float dx, dy, dz;
  int sizebefore=0, sizeafter=0;
  int ipart,idrop;
  unsigned char *compressed_partbuffer, *full_partbuffer;
  int ncompressed_zlib, ncompressed_zlibSAVE, npartfull;



  part_file=parti->file;
  valmin=parti->valmin;
  valmax=parti->valmax;

  strcpy(pp,"%");

  // check if part file is accessible

  if(getfileinfo(part_file,NULL,NULL)!=0){
    printf("Particle file %s does not exist\n",part_file);
    return;
  }

  PARTFILE=fopen(part_file,"rb");
  if(PARTFILE==NULL){
    printf("Particle file %s could not be opened\n",part_file);
    return;
  }

  // set up part compressed file

  if(destdir!=NULL){
    strcpy(part_file_svz,destdir);
    strcat(part_file_svz,partii->filebase);
  }
  else{
    strcpy(part_file_svz,partii->file);
  }
  strcat(part_file_svz,".svz");


//  partstream=fopen(part_file_svz,"wb");
  partstream=EGZ_FOPEN(part_file_svz,"w",1,2);
  
  strcpy(partsizefile_svz,part_file);
  strcat(partsizefile_svz,".sz");
  partsizestream=fopen(partsizefile_svz,"w");

  if(partstream==NULL){
    printf("The file %s could not be opened for writing\n",part_file_svz);
    fclose(PARTFILE);
    return;
  }

  // read and write part header

  strcpy(filetype,"");
  shortlabel=parti->label.shortlabel;
  if(strlen(shortlabel)>0)strcat(filetype,shortlabel);
  trim(filetype);

  strcpy(units,"");
  unit=parti->label.unit;
  if(strlen(unit)>0)strcat(units,unit);
  trim(units);
  printf("Compressing particle file (%s) %s\n",filetype,part_file);
  printf("  using min=%f %s max=%f %s\n",valmin,units,valmax,units);

 // header

 // rr rr rr ii ii(mxframepts)
 // ii(ibar) ii(jbar) ii(kbar)
 // ii(1->ibar+jbar+kbar+3)
 // ii(nb)
 // ii ii ii ii ii ii ii (repeated nb times)
 // ii(nv)
 // ii ii ii ii ii ii ii (repeated nv times)
 // ii(nspr)
 // rr,rr,rr (repeated nspr times)

 // frame
 // rr(time),ii(np),ii(idummy),(ii,...,ii) nspr times
 // (rr,..,rr)(rr,...,rr)(rr,...,rr)(rr,...,rr) np times
 // ii(nsp)
 // rr,rr,rr,(rr if nsp<0)   (abs(nsp( times)

  SEEK;READ(rdummy,3);READ(idummy,2);SEEK;
  sizeafter+=EGZ_FWRITE(rdummy,4,3,partstream);
  sizeafter+=EGZ_FWRITE(idummy,4,2,partstream);

  sizebefore+=4*7;

  ipart=idummy[0];
  mxframepts=idummy[1];

  FORTREAD(idummy,3);
  sizebefore+=20;

  ibar=idummy[0];
  jbar=idummy[1];
  kbar=idummy[2];
//  fseek(PARTFILE,4+4*(ibar+jbar+kbar+3)+4,SEEK_CUR);

  SEEK;READ(&xmin,1);fseek(PARTFILE,4*(ibar-1),SEEK_CUR);READ(&xmax,1);
       READ(&ymin,1);fseek(PARTFILE,4*(jbar-1),SEEK_CUR);READ(&ymax,1);
       READ(&zmin,1);fseek(PARTFILE,4*(kbar-1),SEEK_CUR);READ(&zmax,1);SEEK;
  sizebefore+=4*(ibar+jbar+kbar+5);
  xyzbounds[0]=xmin;
  xyzbounds[1]=xmax;
  xyzbounds[2]=ymin;
  xyzbounds[3]=ymax;
  xyzbounds[4]=zmin;
  xyzbounds[5]=zmax;

  sizeafter+=EGZ_FWRITE(xyzbounds,4,6,partstream);

  FORTREAD(&nb,1);sizebefore+=12;
  if(nb>0){
    fseek(PARTFILE,nb*(4+7*4+4),SEEK_CUR);
    sizebefore+=nb*36;
  }

  FORTREAD(&nv,1);sizebefore+=12;
  if(nv>0){
    fseek(PARTFILE,nv*(4+7*4+4),SEEK_CUR);
    sizebefore+=nv*36;
  }

  FORTREAD(&nspr,1);sizebefore+=12;
  sizeafter+=EGZ_FWRITE(&nspr,4,1,partstream);
  for(i=0;i<nspr;i++){
    FORTREAD(rdummy,3);
    EGZ_FWRITE(rdummy,4,3,partstream);
  }

  if(nspr>0)ispr=malloc(sizeof(int)*nspr);
  xp=malloc(sizeof(float)*mxframepts);
  yp=malloc(sizeof(float)*mxframepts);
  zp=malloc(sizeof(float)*mxframepts);
  bp=malloc(sizeof(float)*mxframepts);
  xsp=malloc(sizeof(float)*mxframepts);
  ysp=malloc(sizeof(float)*mxframepts);
  zsp=malloc(sizeof(float)*mxframepts);
  bsp=malloc(sizeof(float)*mxframepts);

  ixyzbp=malloc(sizeof(unsigned short)*4*mxframepts);
  /*
  ixp=malloc(sizeof(unsigned short)*mxframepts);
  iyp=malloc(sizeof(unsigned short)*mxframepts);
  izp=malloc(sizeof(unsigned short)*mxframepts);
  ibp=malloc(sizeof(unsigned char)*mxframepts);
  */

  ixsp=malloc(sizeof(unsigned short)*mxframepts);
  iysp=malloc(sizeof(unsigned short)*mxframepts);
  izsp=malloc(sizeof(unsigned short)*mxframepts);
  ibsp=malloc(sizeof(unsigned char)*mxframepts);

  npartfull = 7*mxframepts;
  ncompressed_zlib = 1.01*npartfull+1000;
  ncompressed_zlibSAVE=ncompressed_zlib;
  compressed_partbuffer=malloc(ncompressed_zlib);
  full_partbuffer=(unsigned char *)ixyzbp;


  // frames
#define MAXVAL 65535
  dx = (xmax-xmin)/MAXVAL;
  dy = (ymax-ymin)/MAXVAL;
  dz = (zmax-zmin)/MAXVAL;

  // frame
 // rr(time),ii(np),ii(idummy),(ii,...,ii) nspr times
 // (rr,..,rr)(rr,...,rr)(rr,...,rr)(rr,...,rr) np times
 // ii(nsp)
 // rr,rr,rr,(rr if nsp<0)   (abs(nsp( times)

  for(;;){

    int naspr;
    int abs_nsp;


    SEEKBR;
    READBR(&time,1);
    READBR(idummy,2);
    READBR(ispr,1);
    if(nspr>1){
      READBR(ispr+1,nspr-1)
    };
    SEEKBR;

    np=idummy[0];
    idrop=idummy[1];
    if(np>0){
      SEEKBR;READBR(xp,np);READBR(yp,np);READBR(zp,np);READBR(bp,np);SEEKBR;
    }
    else{
      SEEKBR;
      SEEKBR;
    }
    naspr=0;
    for(i=0;i<nspr;i++){
      if(ispr[i]==1)naspr++;
    }
    nsp=0;
    if(naspr>0){

      abs_nsp=nsp;
      if(nsp<0)abs_nsp=-abs_nsp;

      SEEKBR;READBR(&nsp,1);SEEKBR;

      if(nsp>0){
        SEEKBR;READBR(xsp,nsp);READBR(ysp,nsp);READBR(zsp,nsp);SEEKBR;
      }
      else if(nsp<0){
        SEEKBR;READBR(xsp,-nsp);READBR(ysp,-nsp);READBR(zsp,-nsp);READBR(bsp,-nsp);SEEKBR;
      }
      else{
        SEEKBR;SEEKBR;
      }
      for(i=0;i<abs_nsp;i++){
        unsigned short ix, iy, iz;

        ix=0;
        iy=0;
        iz=0;
        if(dx!=0.0)ix=(unsigned short)((xsp[i]-xmin)/dx);
        if(dy!=0.0)iy=(unsigned short)((ysp[i]-ymin)/dy);
        if(dz!=0.0)iz=(unsigned short)((zsp[i]-zmin)/dz);
        ixsp[i]=ix;
        iysp[i]=iy;
        izsp[i]=iz;
        if(ipart==0){
          ibsp[i]=rgb_blue;
        }
        else{
          ibsp[i]=getpartcolor(bsp[i]);
        }
      }
    }
    EGZ_FWRITE(&time,4,1,partstream);
    printf("time=%f\n",time);
    if(nspr>0)EGZ_FWRITE(ispr,4,nspr,partstream);

    if(np>0){

      npartfull=7*np;

      ixp = ixyzbp;
      iyp = ixp+np;
      izp = iyp+np;
      ibp = izp+np;

      for(i=0;i<np;i++){
        unsigned short ix, iy, iz;

        ix=0;
        iy=0;
        iz=0;
        if(dx!=0.0)ix=(unsigned short)((xp[i]-xmin)/dx);
        if(dy!=0.0)iy=(unsigned short)((yp[i]-ymin)/dy);
        if(dz!=0.0)iz=(unsigned short)((zp[i]-zmin)/dz);
        if(ix<0)ix=0;if(ix>MAXVAL)ix=MAXVAL;
        if(iy<0)iy=0;if(iy>MAXVAL)iy=MAXVAL;
        if(iz<0)iz=0;if(iz>MAXVAL)iz=MAXVAL;
        ixp[i]=ix;
        iyp[i]=iy;
        izp[i]=iz;
        ibp[i]=getpartcolor(bp[i]);
      }

      ncompressed_zlib=ncompressed_zlibSAVE;
      npartfull=7*np;

      returncode=compress(compressed_partbuffer, &ncompressed_zlib, full_partbuffer, npartfull);
      EGZ_FWRITE(compressed_partbuffer,1,ncompressed_zlib,partstream);
    }
  }
  fclose(PARTFILE);
  EGZ_FCLOSE(partstream);

  if(nspr>0)free(ispr);
  free(xp);
  free(yp);
  free(zp);
  free(bp);
  free(xsp);
  free(ysp);
  free(zsp);
  free(bsp);
  free(ixsp);
  free(iysp);
  free(izsp);
  free(ibsp);
  free(compressed_partbuffer);




}

/* ------------------ getpartcolor ------------------------ */

unsigned char getpartcolor(float val){
  unsigned char cval;

  if(     val>-0.1&&val<0.1 ){
    cval=rgb_white;
  }
  else if(val>0.9 &&val<1.1 ){
    cval=rgb_yellow;
  }
  else if(val>1.9&&val<2.1){
    cval=rgb_blue;
  }
  else if(val>2.9&&val<3.1){
    cval=rgb_red;
  }
  else if(val>3.9&&val<4.1){
    cval=rgb_green;
  }
  else if(val>4.9&&val<5.1){
    cval=rgb_magenta;
  }
  else if(val>5.9&&val<6.1){
    cval=rgb_cyan;
  }
  else if(val>6.9&&val<7.1){
    cval=rgb_black;
  }
  else{
     cval=rgb_white;
  }
  return cval;

}
#endif
