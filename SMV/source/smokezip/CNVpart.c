#include "options.h"

#ifdef pp_PART
#include "zlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <zlib.h>
#include "svzip.h"
#include "MALLOC.h"
#include "isodefs.h"



void part2iso(part *parti,int *thread_index);


/* ------------------ getpartprop_index ------------------------ */

int getpartprop_index(char *string){
  int i;
  partpropdata *partpropi;

  for(i=0;i<npart5propinfo;i++){
    partpropi = part5propinfo + i;
    if(strcmp(partpropi->label.shortlabel,string)==0)return i;
  }
  return -1;
}

/* ------------------ getpartprop ------------------------ */

partpropdata *getpartprop(char *string){
  int i;
  partpropdata *partpropi;

  for(i=0;i<npart5propinfo;i++){
    partpropi = part5propinfo + i;
    if(strcmp(partpropi->label.shortlabel,string)==0)return partpropi;
  }
  return NULL;
}

/* ------------------ compress_parts ------------------------ */

void compress_parts(void *arg){
  int i;
  int needbounds=0;
  int convertable=0;
  int *thread_index;

  thread_index=(int *)arg;

  for(i=0;i<npartinfo;i++){
    part *parti;

    parti = partinfo + i;
    if(convertable_part(parti)==1){
      convertable=1;
      break;
    }
  }

  for(i=0;i<npart5propinfo;i++){
    partpropdata *propi;

    propi = part5propinfo + i;
    if(propi->setvalmax!=1||propi->setvalmin!=1){
      needbounds=1;
      break;
    }
  }

  if(convertable==1&&needbounds==1&&GLOBget_part_bounds==1){
    Get_Part_Bounds();
  }

  for(i=0;i<npartinfo;i++){
    part *parti;

    parti = partinfo + i;
    if(GLOBautozip==1&&parti->autozip==0)continue;
    convert_part(parti,thread_index);
  }
}

/* ------------------ compress_part2iso ------------------------ */

void *convert_parts2iso(void *arg){
  int i;
  int *thread_index;

  thread_index=(int *)arg;

  LOCK_PART2ISO;
  if(GLOBfirst_part2iso==1){
    GLOBfirst_part2iso=0;

    if(GLOBcleanfiles==1){
      FILE *stream;
      int j;

      stream=fopen(GLOBsmvisofile,"rb");
      if(stream!=NULL){
        fclose(stream);
        PRINTF("  Removing %s\n",GLOBsmvisofile);
        UNLINK(GLOBsmvisofile);
        LOCK_COMPRESS;
        GLOBfilesremoved++;
        UNLOCK_COMPRESS;
      }
      for(j=0;j<npartinfo;j++){
        part *parti;

        parti = partinfo + j;

        for(i=0;i<npart5propinfo;i++){
          partpropdata *propi;
          flowlabels *labels;
          char isofilename[1024];

          propi = part5propinfo + i;
          labels = &propi->label;
          strcpy(isofilename,parti->file);
          strcat(isofilename,"_");
          strcat(isofilename,labels->shortlabel);
          strcat(isofilename,".tiso");

          stream=fopen(isofilename,"rb");
          if(stream!=NULL){
            fclose(stream);
            PRINTF("  Removing %s\n",isofilename);
            UNLINK(isofilename);
            LOCK_COMPRESS;
            GLOBfilesremoved++;
            UNLOCK_COMPRESS;
          }
        }
      }
      UNLOCK_PART2ISO;
      return NULL;
    }
  }
  UNLOCK_PART2ISO;

  if(GLOBcleanfiles==1)return NULL;

  if(GLOBpartfile2iso==1){
    for(i=0;i<npartinfo;i++){
      part *parti;

      parti = partinfo + i;
      LOCK_PART2ISO;
      if(parti->inuse_part2iso==1){
        UNLOCK_PART2ISO;
        continue;
      }
      parti->inuse_part2iso=1;
      UNLOCK_PART2ISO;

      part2iso(parti,thread_index);
    }
  }
  return NULL;
}

/* ------------------ convertable_part ------------------------ */

int convertable_part(part *parti){
  char partfile_svz[256];
  char *partfile;
  FILE *stream;

  partfile=parti->file;

  if(getfileinfo(partfile,NULL,NULL)!=0)return 0;

  stream=fopen(partfile,"rb");
  if(stream==NULL)return 0;
  fclose(stream);

  // set up compressed particle file

  strcpy(partfile_svz,"");
  if(GLOBdestdir!=NULL){
    strcpy(partfile_svz,GLOBdestdir);
  }
  strcat(partfile_svz,partfile);
  strcat(partfile_svz,".svz");
  stream=fopen(partfile_svz,"wb");
  if(stream==NULL)return 0;
  return 1;
}

/* ------------------ compress_part ------------------------ */

void convert_part(part *parti, int *thread_index){
  FILE *PARTFILEstream,*partstream,*partsizestream;
  char *partfile, partfile_svz[256], partsizefile_svz[256];
  FILE_SIZE lenfile;
  int unit;
  int nclasses;
  int *nquantities, *npoints;
  float time_local;
  int error;
  int endiandata;
  int *tagdata;
  float *pdata;
  int sizebefore=0, sizeafter=0, size;
  int one=1, zero=0, fileversion=0, fdsversion;
  unsigned int *int_buffer_uncompressed;
  unsigned char *int_buffer_compressed;
  unsigned char *char_buffer_uncompressed,*char_buffer_compressed;
  meshdata *partmesh;
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
  // time_local
  // points_1, ..., points_nclasses
  // ntotal_int, ntotal_char
  // ncompressedpoints
  // xyz_1, ..., xyz_ncompressedpoints
  // ncompresseddata
  // data_1, ..., data_ncompresseddata


#ifdef pp_THREAD
  {
    int fileindex;

    fileindex = parti + 1 - partinfo;
    sprintf(threadinfo[*thread_index].label,"prt5 %i",fileindex);
  }
#endif

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
    fprintf(stderr,"*** Warning: The particle file %s does not exist\n",partfile);
    return;
  }

  PARTFILEstream=fopen(partfile,"rb");
  if(PARTFILEstream==NULL){
    fprintf(stderr,"*** Warning: The particle file %s could not be opened\n",partfile);
    return;
  }

  // set up compressed particle file

  strcpy(partfile_svz,"");
  if(GLOBdestdir!=NULL){
    strcpy(partfile_svz,GLOBdestdir);
  }
  strcat(partfile_svz,partfile);
  strcat(partfile_svz,".svz");
  partstream=fopen(partfile_svz,"wb");
  if(partstream==NULL){
    fprintf(stderr,"*** Warning: The file %s could not be opened for writing\n",partfile_svz);
    fclose(PARTFILEstream);
    return;
  }

  // set up compressed particle size file

  strcpy(partsizefile_svz,"");
  if(GLOBdestdir!=NULL){
    strcpy(partsizefile_svz,GLOBdestdir);
  }
  strcpy(partsizefile_svz,partfile);
  strcat(partsizefile_svz,".szz");
  partsizestream=fopen(partsizefile_svz,"w");

  if(partsizestream==NULL){
    fprintf(stderr,"*** Warning: The file %s could not be opened for writing\n",partsizefile_svz);
    fclose(PARTFILEstream);
    fclose(partstream);
    return;
  }

  if(GLOBcleanfiles==1){
    partstream=fopen(partfile_svz,"rb");
    if(partstream!=NULL){
      fclose(partstream);
      PRINTF("  Removing %s\n",partfile_svz);
      UNLINK(partfile_svz);
      LOCK_COMPRESS;
      GLOBfilesremoved++;
      UNLOCK_COMPRESS;
    }
    partsizestream=fopen(partsizefile_svz,"rb");
    if(partsizestream!=NULL){
      fclose(partsizestream);
      PRINTF("  Removing %s\n",partsizefile_svz);
      UNLINK(partsizefile_svz);
      LOCK_COMPRESS;
      GLOBfilesremoved++;
      UNLOCK_COMPRESS;
    }
    return;
  }

  if(GLOBoverwrite_part==0){
    partstream=fopen(partfile_svz,"rb");
    if(partstream!=NULL){
      fclose(partstream);
      fprintf(stderr,"** Warning: The file %s exists.\n",partfile_svz);
      fprintf(stderr,"     Use the -f option to overwrite smokezip compressed files\n");
      return;
    }
  }

  PRINTF("  Compressing %s\n\n",parti->file);

#define BUFFER_SIZE 1000000
  NewMemory((void **)&pdata,BUFFER_SIZE*sizeof(float));
  NewMemory((void **)&tagdata,BUFFER_SIZE*sizeof(int));
  NewMemory((void **)&int_buffer_uncompressed,BUFFER_SIZE*sizeof(unsigned int));
  NewMemory((void **)&int_buffer_compressed,BUFFER_SIZE*sizeof(unsigned char));
  NewMemory((void **)&char_buffer_uncompressed,BUFFER_SIZE*sizeof(unsigned char));
  NewMemory((void **)&char_buffer_compressed,BUFFER_SIZE*sizeof(unsigned char));

  lenfile=strlen(parti->file);
  LOCK_COMPRESS;
  FORTget_file_unit(&unit,&parti->unit_start);
  FORTopenpart(parti->file,&unit,&error,lenfile);
  UNLOCK_COMPRESS;

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

    FORTgetpartdataframe(&unit,&nclasses,nquantities,npoints,&time_local,tagdata,pdata,&size,&error);
    if(error!=0)break;

    fwrite(&time_local,4,1,partstream);
    fwrite(npoints,4,nclasses,partstream);
    sizeafter+=4*(1+nclasses);
    sizebefore+=size;

//    fprintf(partsizestream,"%f\n",time_local);

    vals=pdata;
    tagdataptr=(unsigned int *)tagdata;
    xbuff = int_buffer_uncompressed;
    cbuff = char_buffer_uncompressed;
    ntotal_int=0;
    ntotal_char=0;
    for(j=0;j<nclasses;j++){
      partclassdata *classi;
      int k;

//      fprintf(partsizestream," %i ",npoints[j]);
  //    if(j<nclasses-1){
  //      fprintf(partsizestream," \n");
  //    }

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
        partpropdata *propi;
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
      compress_zlib(int_buffer_compressed, &ncompressed_zlib, (unsigned char *)int_buffer_uncompressed, 4*ntotal_int);
      ncompressed_int = ncompressed_zlib;
      sizeafter+=(4+ncompressed_int);
      fwrite(&ncompressed_int,4,1,partstream);
      fwrite(&int_buffer_compressed,1,ncompressed_int,partstream);
    }

    ncompressed_char=0;
    if(ntotal_char>0){
      ncompressed_zlib=BUFFER_SIZE;
      compress_zlib(char_buffer_compressed, &ncompressed_zlib, char_buffer_uncompressed, ntotal_char);
      ncompressed_char = ncompressed_zlib;
      sizeafter+=(4+ncompressed_char);
      fwrite(&ncompressed_char,4,1,partstream);
      fwrite(&char_buffer_compressed,1,ncompressed_char,partstream);
    }

//    fprintf(partsizestream," %i %i %i %i\n",(int)ntotal_int,(int)ntotal_char,ncompressed_int,ncompressed_char);

    data_loc=sizebefore;
    percent_done=100.0*(float)data_loc/(float)parti->filesize;
    if(percent_done>percent_next){
      PRINTF(" %i%s",percent_next,GLOBpp);
      LOCK_COMPRESS;
      FFLUSH();
      UNLOCK_COMPRESS;
      percent_next+=10;
    }
    count++;
  }

  PRINTF(" 100%s completed\n",GLOBpp);
  FREEMEMORY(nquantities);
  FREEMEMORY(npoints);
  FREEMEMORY(pdata);
  FREEMEMORY(tagdata);
  FREEMEMORY(int_buffer_uncompressed);
  FREEMEMORY(int_buffer_compressed);
  FREEMEMORY(char_buffer_uncompressed);
  FREEMEMORY(char_buffer_compressed);

  LOCK_COMPRESS;
  FORTclosefortranfile(&unit);
  UNLOCK_COMPRESS;
  fclose(partstream);
//  fclose(partsizestream);
  {
    char before_label[256],after_label[256];

    getfilesizelabel(sizebefore,before_label);
    getfilesizelabel(sizeafter,after_label);
    PRINTF("    records=%i, ",count);
    PRINTF("Sizes: original=%s, ",before_label);
    PRINTF("compressed=%s (%4.1f%s reduction)\n",after_label,(float)sizebefore/(float)sizeafter,GLOBx);
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

  PRINTF("Determining particle file bounds\n");

  for(i=0;i<npart5propinfo;i++){
    partpropdata *propi;

    propi = part5propinfo + i;

    NewMemory((void **)&propi->histogram,sizeof(histogramdata));
    init_histogram(propi->histogram,NHIST_BUCKETS);
  }

  NewMemory((void **)&pdata,1000000*sizeof(float));
  NewMemory((void **)&tagdata,1000000*sizeof(int));

  for(i=0;i<npartinfo;i++){
    part *parti;
    FILE_SIZE lenfile;
    int unit;
    int error1;
    int nclasses;
    int *nquantities, *npoints;
    float time_local;
    int error, size;
    int nquantities_total;
    int j;

    parti = partinfo + i;
    PRINTF("  Examining %s\n",parti->file);
    lenfile=strlen(parti->file);
    LOCK_COMPRESS;
    unit=15;
    FORTget_file_unit(&unit,&parti->unit_start);
    FORTopenpart(parti->file,&unit,&error1,lenfile);
    UNLOCK_COMPRESS;

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
      LOCK_COMPRESS;
      FORTclosefortranfile(&unit);
      UNLOCK_COMPRESS;
      continue;
    }

    error=0;
    for(;;){
      float *x, *y, *z, *vals;
      int k;

      FORTgetpartdataframe(&unit,&nclasses,nquantities,npoints,&time_local,tagdata,pdata,&size,&error);
      if(error!=0)break;

      vals=pdata;
      for(j=0;j<nclasses;j++){
        partclassdata *classj;

        if(npoints[j]==0)continue;
        classj=parti->classptr[j];
        x = vals;
        y = x + npoints[j];
        z = y + npoints[j];
        vals = z + npoints[j];
        for(k=0;k<nquantities[j];k++){
          partpropdata *propi;

          propi=getpartprop(classj->labels[k].shortlabel);
          update_histogram(vals,npoints[j],propi->histogram);

          vals += npoints[j];
        }
      }
    }

    FREEMEMORY(nquantities);
    FREEMEMORY(npoints);

    LOCK_COMPRESS;
    FORTclosefortranfile(&unit);
    UNLOCK_COMPRESS;
  }

  FREEMEMORY(pdata);
  FREEMEMORY(tagdata);

  for(i=0;i<npart5propinfo;i++){
    partpropdata *propi;

    propi = part5propinfo + i;

    propi->valmax=get_histogram_value(propi->histogram,0.99);
    propi->valmin=get_histogram_value(propi->histogram,0.01);
    propi->setvalmax=1;
    propi->setvalmin=1;
    PRINTF(" %s min: %f max: %f\n",propi->label.shortlabel,propi->valmin,propi->valmax);
    FREEMEMORY(propi->histogram->buckets);
    FREEMEMORY(propi->histogram->buckets_2d);
    FREEMEMORY(propi->histogram);
  }
}

#define IJKVAL(ix,iy,iz) ((ix) + (iy)*nx2 + (iz)*nx2*ny2)

  /* ------------------ part2iso ------------------------ */

void part2iso(part *parti, int *thread_index){
  float *pdata;
  int *tagdata;
  int fdsversion;

  int endiandata;

  int blocknumber;
  FILE_SIZE len_partfile;
  int unit;
  int error1;
  int nclasses;
  int *nquantities, *npoints, *partindex;
  float time_local;
  int error, size;
  int j;
  int npartcount, i;
  meshdata *partmesh;
  int nx, ny, nz;
  char *isofile, tisofile[1024];
  char isolonglabel[32], isoshortlabel[32], isounits[32];
  int nlevels;
  float levels[1];
  int reduce_triangles=1;
  float *xpltcell, *ypltcell, *zpltcell;
  int data2flag=1;
  float *partcount;
  FILE *SMVISOFILE=NULL;
  int nx2, ny2, nz2;
  float xmin, ymin, zmin;
  partpropdata *part5propinfo_copy;
  int percent_done;
  float file_size;
  int percent_next=10;

  parti->compressed2=0;
#ifdef pp_THREAD
  if(GLOBcleanfiles==0){
    int fileindex;

    fileindex = parti + 1 - partinfo;
    sprintf(threadinfo[*thread_index].label,"prt2iso %i",fileindex);
  }
#else
  PRINTF("Converting %s to\n",parti->file);
#endif

  endiandata=getendian();
  if(endianswitch==1)endiandata=1-endiandata;

  NewMemory((void **)&pdata,1000000*sizeof(float));
  NewMemory((void **)&tagdata,1000000*sizeof(int));
  NewMemory((void **)&partindex,1000000*sizeof(int));

  len_partfile=strlen(parti->file);
  LOCK_COMPRESS;
  FORTget_file_unit(&unit,&parti->unit_start);
  FORTopenpart(parti->file,&unit,&error1,len_partfile);
  UNLOCK_COMPRESS;

  FORTgetpartheader1(&unit,&nclasses,&fdsversion,&size);
  NewMemory((void **)&nquantities,nclasses*sizeof(int));
  NewMemory((void **)&npoints,nclasses*sizeof(int));

  FORTgetpartheader2(&unit,&nclasses,nquantities,&size);

  partmesh = parti->partmesh;

  blocknumber = partmesh-meshinfo + 1;

  nx = partmesh->ibar;
  ny = partmesh->jbar;
  nz = partmesh->kbar;

  nx2 = nx+2;
  ny2 = ny+2;
  nz2 = nz+2;

  npartcount = nx2*ny2*nz2;

  xmin = partmesh->xbar0-partmesh->dx;
  ymin = partmesh->ybar0-partmesh->dy;
  zmin = partmesh->zbar0-partmesh->dz;

  xpltcell = partmesh->xpltcell;
  ypltcell = partmesh->ypltcell;
  zpltcell = partmesh->zpltcell;

  NewMemory((void **)&isofile,strlen(parti->file)+5);
  strcpy(isofile,parti->file);
  strcat(isofile,".iso");

  NewMemory((void **)&tisofile,strlen(parti->file)+6);
  strcpy(tisofile,parti->file);
  strcat(tisofile,".tiso");

  strcpy(isolonglabel,"particle boundary");
  strcpy(isoshortlabel,"pbound");
  strcpy(isounits,"");

  nlevels=1;
  levels[0]=0.5;

  if(npart5propinfo>0)NewMemory((void **)&part5propinfo_copy,npart5propinfo*sizeof(partpropdata));

  for(i=0;i<npart5propinfo;i++){
    partpropdata *propi;

    propi = part5propinfo_copy + i;
    propi->used=0;
  }
  for(j=0;j<nclasses;j++){
    int k;

    for(k=0;k<nquantities[j];k++){
      partclassdata *classj;
      partpropdata *propi;

      classj=parti->classptr[j];
      propi=part5propinfo_copy+getpartprop_index(classj->labels[k].shortlabel);
      propi->used=1;
    }
  }

  NewMemory((void **)&partcount,npartcount*sizeof(float));

  for(i=0;i<npart5propinfo;i++){
    partpropdata *propi;

    propi = part5propinfo_copy + i;
    if(propi->used==0)continue;

    NewMemory((void **)&propi->partvals,npartcount*sizeof(float));
  }

  CCisoheader(isofile,isolonglabel,isoshortlabel,isounits,levels,&nlevels,&error);

  LOCK_PART2ISO;
  if(GLOBfirst_part2iso_smvopen==1){
    GLOBfirst_part2iso_smvopen=0;
    SMVISOFILE=fopen(GLOBsmvisofile,"w");
  }
  else{
    SMVISOFILE=fopen(GLOBsmvisofile,"a");
  }

  fprintf(SMVISOFILE,"ISOF %i\n",blocknumber);
  fprintf(SMVISOFILE," %s\n",isofile);
  fprintf(SMVISOFILE," %s\n",isolonglabel);

  fprintf(SMVISOFILE," %s\n",isoshortlabel);
  fprintf(SMVISOFILE," %s\n",isounits);
  fprintf(SMVISOFILE,"\n");

#ifndef pp_THREAD
  PRINTF("  %s\n",isofile);
#endif

  for(i=0;i<npart5propinfo;i++){
    partpropdata *propi,*propi_ro;
    flowlabels *labels;

    propi_ro = part5propinfo + i;
    propi = part5propinfo_copy + i;
    if(propi->used==0)continue;

    labels = &propi_ro->label;
    strcpy(propi->isofilename,parti->file);
    strcat(propi->isofilename,"_");
    strcat(propi->isofilename,labels->shortlabel);
    strcat(propi->isofilename,".tiso");
    CCtisoheader(propi->isofilename, labels->longlabel, labels->shortlabel, labels->unit, levels, &nlevels, &error);
#ifndef pp_THREAD
    PRINTF("  %s\n",propi->isofilename);
#endif

    fprintf(SMVISOFILE,"TISOF %i\n",blocknumber);
    fprintf(SMVISOFILE," %s\n",propi->isofilename);
    fprintf(SMVISOFILE," %s\n",isolonglabel);
    fprintf(SMVISOFILE," %s\n",isoshortlabel);
    fprintf(SMVISOFILE," %s\n",isounits);
    fprintf(SMVISOFILE," %s\n",labels->longlabel);
    fprintf(SMVISOFILE," %s\n",labels->shortlabel);
    fprintf(SMVISOFILE," %s\n",labels->unit);
    fprintf(SMVISOFILE," \n");
  }
  fclose(SMVISOFILE);
  UNLOCK_PART2ISO;

#ifndef pp_THREAD
  PRINTF(" ");
#endif
  error=0;
  file_size=0.0;
  for(;;){
    float *x, *y, *z, *vals;
    int k;

    FORTgetpartdataframe(&unit,&nclasses,nquantities,npoints,&time_local,tagdata,pdata,&size,&error);

    file_size+=size;

    percent_done=100.0*(float)file_size/(float)parti->filesize;
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
    if(error!=0)break;

    for(j=0;j<npartcount;j++){
      partcount[j]=0.0;
    }
    for(i=0;i<npart5propinfo;i++){
      partpropdata *propi;

      propi = part5propinfo_copy + i;
      if(propi->used==0)continue;

      for(j=0;j<npartcount;j++){
        propi->partvals[j]=0.0;
      }
    }

    vals=pdata;
    for(j=0;j<nclasses;j++){
      partclassdata *classj;

      if(npoints[j]==0)continue;
      classj=parti->classptr[j];
      x = vals;
      y = x + npoints[j];
      z = y + npoints[j];
      vals = z + npoints[j];

// construct 3D particle density array

      for(i=0;i<npoints[j];i++){
        int ix, iy, iz;
        int ijkval;

        GETINDEX(ix,x[i],xmin,partmesh->dx,nx2);
        GETINDEX(iy,y[i],ymin,partmesh->dy,ny2);
        GETINDEX(iz,z[i],zmin,partmesh->dz,nz2);
        ijkval = IJKVAL(ix,iy,iz);
        partindex[i]=ijkval;
        partcount[ijkval]++;
      }
      for(k=0;k<nquantities[j];k++){
        partpropdata *propi;

        propi=part5propinfo_copy+getpartprop_index(classj->labels[k].shortlabel);

        for(i=0;i<npoints[j];i++){
          int ijkval;

          ijkval = partindex[i];
          propi->partvals[ijkval]+=vals[i];
        }
        for(i=0;i<npartcount;i++){
          if(partcount[i]>0){
            propi->partvals[i]/=partcount[i];
          }
        }
        vals += npoints[j];
      }
    }
    CCisosurface2file(isofile, &time_local, partcount, NULL, levels, &nlevels,
        xpltcell, &nx2, ypltcell, &ny2, zpltcell, &nz2,
        &reduce_triangles, &error);

    for(i=0;i<npart5propinfo;i++){
      partpropdata *propi;

      propi = part5propinfo_copy + i;
      if(propi->used==0)continue;

      CCisosurfacet2file(propi->isofilename, &time_local, partcount, &data2flag, propi->partvals, NULL, levels, &nlevels,
            xpltcell, &nx2, ypltcell, &ny2, zpltcell, &nz2,
            &reduce_triangles, &error);
    }
  }

#ifdef pp_THREAD
  {
    int nconv=1;
    int lenfile;
    char **summaries;

    for(i=0;i<npart5propinfo;i++){
      partpropdata *propi;

      propi = part5propinfo_copy + i;
      if(propi->used==1)nconv++;
    }

    NewMemory((void **)&summaries,nconv*sizeof(char *));

    lenfile=strlen(isofile);
    NewMemory((void **)&summaries[0],lenfile+1);
    strcpy(summaries[0],isofile);

    parti->summaries=summaries;
    parti->nsummaries=nconv;

    nconv=1;
    for(i=0;i<npart5propinfo;i++){
      partpropdata *propi;

      propi = part5propinfo_copy + i;
      if(propi->used==0)continue;
      lenfile=strlen(propi->isofilename);
      NewMemory((void **)&summaries[nconv],lenfile+1);
      strcpy(summaries[nconv],propi->isofilename);
      nconv++;
    }
    parti->compressed2=1;
    threadinfo[*thread_index].stat=-1;
  }
#else
  PRINTF(" 100%s completed\n",GLOBpp);
#endif

  FREEMEMORY(nquantities);
  FREEMEMORY(npoints);
  FREEMEMORY(partcount);
  FREEMEMORY(isofile);
  LOCK_COMPRESS;
  FORTclosefortranfile(&unit);
  UNLOCK_COMPRESS;

  FREEMEMORY(pdata);
  FREEMEMORY(tagdata);
  FREEMEMORY(partindex);

  for(i=0;i<npart5propinfo;i++){
    partpropdata *propi;

    propi = part5propinfo_copy + i;
    if(propi->used==0)continue;

    FREEMEMORY(propi->partvals);
    FREEMEMORY(propi->partvals);
  }
  if(npart5propinfo>0){
    FREEMEMORY(part5propinfo_copy);
  }
}

/* ------------------ part2object ------------------------ */

void part2object(part *parti, int *thread_index){
  float *pdata;
  int *tagdata;
  int fdsversion;

  int endiandata;

  int blocknumber;
  FILE_SIZE len_partfile;
  int unit;
  int error1;
  int nclasses;
  int *nquantities, *npoints, *partindex;
  float time_local;
  int error, size;
  int j;
  int npartcount, i;
  meshdata *partmesh;
  int nx, ny, nz;
  char *isofile, tisofile[1024];
  char isolonglabel[32], isoshortlabel[32], isounits[32];
  int nlevels;
  float levels[1];
  int reduce_triangles=1;
  float *xpltcell, *ypltcell, *zpltcell;
  int data2flag=1;
  float *partcount;
  FILE *SMVISOFILE=NULL;
  int nx2, ny2, nz2;
  float xmin, ymin, zmin;
  partpropdata *part5propinfo_copy;
  int percent_done;
  float file_size;
  int percent_next=10;

  parti->compressed2=0;
#ifdef pp_THREAD
  if(GLOBcleanfiles==0){
    int fileindex;

    fileindex = parti + 1 - partinfo;
    sprintf(threadinfo[*thread_index].label,"prt2iso %i",fileindex);
  }
#else
  PRINTF("Converting %s to\n",parti->file);
#endif

  endiandata=getendian();
  if(endianswitch==1)endiandata=1-endiandata;

  NewMemory((void **)&pdata,1000000*sizeof(float));
  NewMemory((void **)&tagdata,1000000*sizeof(int));
  NewMemory((void **)&partindex,1000000*sizeof(int));

  len_partfile=strlen(parti->file);
  LOCK_COMPRESS;
  FORTget_file_unit(&unit,&parti->unit_start);
  FORTopenpart(parti->file,&unit,&error1,len_partfile);
  UNLOCK_COMPRESS;

  FORTgetpartheader1(&unit,&nclasses,&fdsversion,&size);
  NewMemory((void **)&nquantities,nclasses*sizeof(int));
  NewMemory((void **)&npoints,nclasses*sizeof(int));

  FORTgetpartheader2(&unit,&nclasses,nquantities,&size);

  partmesh = parti->partmesh;

  blocknumber = partmesh-meshinfo + 1;

  nx = partmesh->ibar;
  ny = partmesh->jbar;
  nz = partmesh->kbar;

  nx2 = nx+2;
  ny2 = ny+2;
  nz2 = nz+2;

  npartcount = nx2*ny2*nz2;

  xmin = partmesh->xbar0-partmesh->dx;
  ymin = partmesh->ybar0-partmesh->dy;
  zmin = partmesh->zbar0-partmesh->dz;

  xpltcell = partmesh->xpltcell;
  ypltcell = partmesh->ypltcell;
  zpltcell = partmesh->zpltcell;

  NewMemory((void **)&isofile,strlen(parti->file)+5);
  strcpy(isofile,parti->file);
  strcat(isofile,".iso");

  NewMemory((void **)&tisofile,strlen(parti->file)+6);
  strcpy(tisofile,parti->file);
  strcat(tisofile,".tiso");

  strcpy(isolonglabel,"particle boundary");
  strcpy(isoshortlabel,"pbound");
  strcpy(isounits,"");

  nlevels=1;
  levels[0]=0.5;

  if(npart5propinfo>0)NewMemory((void **)&part5propinfo_copy,npart5propinfo*sizeof(partpropdata));

  for(i=0;i<npart5propinfo;i++){
    partpropdata *propi;

    propi = part5propinfo_copy + i;
    propi->used=0;
  }
  for(j=0;j<nclasses;j++){
    int k;

    for(k=0;k<nquantities[j];k++){
      partclassdata *classj;
      partpropdata *propi;

      classj=parti->classptr[j];
      propi=part5propinfo_copy+getpartprop_index(classj->labels[k].shortlabel);
      propi->used=1;
    }
  }

  NewMemory((void **)&partcount,npartcount*sizeof(float));

  for(i=0;i<npart5propinfo;i++){
    partpropdata *propi;

    propi = part5propinfo_copy + i;
    if(propi->used==0)continue;

    NewMemory((void **)&propi->partvals,npartcount*sizeof(float));
  }

  CCisoheader(isofile,isolonglabel,isoshortlabel,isounits,levels,&nlevels,&error);

  LOCK_PART2ISO;
  if(GLOBfirst_part2iso_smvopen==1){
    GLOBfirst_part2iso_smvopen=0;
    SMVISOFILE=fopen(GLOBsmvisofile,"w");
  }
  else{
    SMVISOFILE=fopen(GLOBsmvisofile,"a");
  }

  fprintf(SMVISOFILE,"ISOF %i\n",blocknumber);
  fprintf(SMVISOFILE," %s\n",isofile);
  fprintf(SMVISOFILE," %s\n",isolonglabel);

  fprintf(SMVISOFILE," %s\n",isoshortlabel);
  fprintf(SMVISOFILE," %s\n",isounits);
  fprintf(SMVISOFILE,"\n");

#ifndef pp_THREAD
  PRINTF("  %s\n",isofile);
#endif

  for(i=0;i<npart5propinfo;i++){
    partpropdata *propi,*propi_ro;
    flowlabels *labels;

    propi_ro = part5propinfo + i;
    propi = part5propinfo_copy + i;
    if(propi->used==0)continue;

    labels = &propi_ro->label;
    strcpy(propi->isofilename,parti->file);
    strcat(propi->isofilename,"_");
    strcat(propi->isofilename,labels->shortlabel);
    strcat(propi->isofilename,".tiso");
    CCtisoheader(propi->isofilename, labels->longlabel, labels->shortlabel, labels->unit, levels, &nlevels, &error);
#ifndef pp_THREAD
    PRINTF("  %s\n",propi->isofilename);
#endif

    fprintf(SMVISOFILE,"TISOF %i\n",blocknumber);
    fprintf(SMVISOFILE," %s\n",propi->isofilename);
    fprintf(SMVISOFILE," %s\n",isolonglabel);
    fprintf(SMVISOFILE," %s\n",isoshortlabel);
    fprintf(SMVISOFILE," %s\n",isounits);
    fprintf(SMVISOFILE," %s\n",labels->longlabel);
    fprintf(SMVISOFILE," %s\n",labels->shortlabel);
    fprintf(SMVISOFILE," %s\n",labels->unit);
    fprintf(SMVISOFILE," \n");
  }
  fclose(SMVISOFILE);
  UNLOCK_PART2ISO;

#ifndef pp_THREAD
  PRINTF(" ");
#endif
  error=0;
  file_size=0.0;
  for(;;){
    float *x, *y, *z, *vals;
    int k;

    FORTgetpartdataframe(&unit,&nclasses,nquantities,npoints,&time_local,tagdata,pdata,&size,&error);

    file_size+=size;

    percent_done=100.0*(float)file_size/(float)parti->filesize;
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
    if(error!=0)break;

    for(j=0;j<npartcount;j++){
      partcount[j]=0.0;
    }
    for(i=0;i<npart5propinfo;i++){
      partpropdata *propi;

      propi = part5propinfo_copy + i;
      if(propi->used==0)continue;

      for(j=0;j<npartcount;j++){
        propi->partvals[j]=0.0;
      }
    }

    vals=pdata;
    for(j=0;j<nclasses;j++){
      partclassdata *classj;

      if(npoints[j]==0)continue;
      classj=parti->classptr[j];
      x = vals;
      y = x + npoints[j];
      z = y + npoints[j];
      vals = z + npoints[j];

      // construct 3D particle density array

      for(i=0;i<npoints[j];i++){
        int ix, iy, iz;
        int ijkval;

        GETINDEX(ix,x[i],xmin,partmesh->dx,nx2);
        GETINDEX(iy,y[i],ymin,partmesh->dy,ny2);
        GETINDEX(iz,z[i],zmin,partmesh->dz,nz2);
        ijkval = IJKVAL(ix,iy,iz);
        partindex[i]=ijkval;
        partcount[ijkval]++;
      }
      for(k=0;k<nquantities[j];k++){
        partpropdata *propi;

        propi=part5propinfo_copy+getpartprop_index(classj->labels[k].shortlabel);

        for(i=0;i<npoints[j];i++){
          int ijkval;

          ijkval = partindex[i];
          propi->partvals[ijkval]+=vals[i];
        }
        for(i=0;i<npartcount;i++){
          if(partcount[i]>0){
            propi->partvals[i]/=partcount[i];
          }
        }
        vals += npoints[j];
      }
    }
    CCisosurface2file(isofile, &time_local, partcount, NULL, levels, &nlevels,
      xpltcell, &nx2, ypltcell, &ny2, zpltcell, &nz2,
      &reduce_triangles, &error);

    for(i=0;i<npart5propinfo;i++){
      partpropdata *propi;

      propi = part5propinfo_copy + i;
      if(propi->used==0)continue;

      CCisosurfacet2file(propi->isofilename, &time_local, partcount, &data2flag, propi->partvals, NULL, levels, &nlevels,
        xpltcell, &nx2, ypltcell, &ny2, zpltcell, &nz2,
        &reduce_triangles, &error);
    }
  }

#ifdef pp_THREAD
  {
    int nconv=1;
    int lenfile;
    char **summaries;

    for(i=0;i<npart5propinfo;i++){
      partpropdata *propi;

      propi = part5propinfo_copy + i;
      if(propi->used==1)nconv++;
    }

    NewMemory((void **)&summaries,nconv*sizeof(char *));

    lenfile=strlen(isofile);
    NewMemory((void **)&summaries[0],lenfile+1);
    strcpy(summaries[0],isofile);

    parti->summaries=summaries;
    parti->nsummaries=nconv;

    nconv=1;
    for(i=0;i<npart5propinfo;i++){
      partpropdata *propi;

      propi = part5propinfo_copy + i;
      if(propi->used==0)continue;
      lenfile=strlen(propi->isofilename);
      NewMemory((void **)&summaries[nconv],lenfile+1);
      strcpy(summaries[nconv],propi->isofilename);
      nconv++;
    }
    parti->compressed2=1;
    threadinfo[*thread_index].stat=-1;
  }
#else
  PRINTF(" 100%s completed\n",GLOBpp);
#endif

  FREEMEMORY(nquantities);
  FREEMEMORY(npoints);
  FREEMEMORY(partcount);
  FREEMEMORY(isofile);
  LOCK_COMPRESS;
  FORTclosefortranfile(&unit);
  UNLOCK_COMPRESS;

  FREEMEMORY(pdata);
  FREEMEMORY(tagdata);
  FREEMEMORY(partindex);

  for(i=0;i<npart5propinfo;i++){
    partpropdata *propi;

    propi = part5propinfo_copy + i;
    if(propi->used==0)continue;

    FREEMEMORY(propi->partvals);
    FREEMEMORY(propi->partvals);
  }
  if(npart5propinfo>0){
    FREEMEMORY(part5propinfo_copy);
  }
}

#endif
