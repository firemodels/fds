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

#define FORTgetplot3dq _F(getplot3dq)

STDCALLF FORTgetplot3dq(char *qfilename, int *nx, int *ny, int *nz, float *qq, int *error, int *isotest, FILE_SIZE filelen);

/* ------------------ convert_plot3d ------------------------ */

int convert_plot3d(plot3d *plot3di){

  char plot3dfile_svz[1024];
  int fileversion, one, zero;
  char *plot3d_file;
  int version_local;
  char filetype[1024];
  char *shortlabel;
  FILE *PLOT3DFILE, *plot3dstream;
  int sizebefore, sizeafter;
  float *plot3dframe_data=NULL;
  unsigned char *plot3dframe_compressed=NULL, *plot3dframe_uncompressed=NULL;
  uLongf ncompressed_zlib;
  float minmax[10];
  int ijkbar[3];

  plot3di->compressed=0;
  plot3d_file=plot3di->file;
  version_local=plot3di->version;

  fileversion = 1;
  one = 1;
  zero=0;

  // check if plot3d file is accessible

  strcpy(filetype,"");
  shortlabel=plot3di->labels[0].shortlabel;
  if(strlen(shortlabel)>0)strcat(filetype,shortlabel);
  trim_back(filetype);

  if(getfileinfo(plot3d_file,NULL,NULL)!=0){
    fprintf(stderr,"*** Warning: The file %s does not exist\n",plot3d_file);
    return 0;
  }

  PLOT3DFILE=fopen(plot3d_file,"rb");
  if(PLOT3DFILE==NULL){
    fprintf(stderr,"*** Warning: The file %s could not be opened\n",plot3d_file);
    return 0;
  }

  // set up plot3d compressed file

  if(GLOBdestdir!=NULL){
    strcpy(plot3dfile_svz,GLOBdestdir);
    strcat(plot3dfile_svz,plot3di->filebase);
  }
  else{
    strcpy(plot3dfile_svz,plot3di->file);
  }
  {
    int lensvz;

    lensvz = strlen(plot3dfile_svz);

    if(lensvz>4){
      strcat(plot3dfile_svz,".svz");
    }
  }

  if(GLOBcleanfiles==1){
    plot3dstream=fopen(plot3dfile_svz,"rb");
    if(plot3dstream!=NULL){
      fclose(plot3dstream);
      PRINTF("  Removing %s\n",plot3dfile_svz);
      UNLINK(plot3dfile_svz);
      LOCK_COMPRESS;
      GLOBfilesremoved++;
      UNLOCK_COMPRESS;
    }
    return 0;
  }

  if(GLOBoverwrite_plot3d==0){
    plot3dstream=fopen(plot3dfile_svz,"rb");
    if(plot3dstream!=NULL){
      fclose(plot3dstream);
      fprintf(stderr,"*** Warning: The file %s exists.\n",plot3dfile_svz);
      fprintf(stderr,"     Use the -f option to overwrite smokezip compressed files\n");
      return 0;
    }
  }

  plot3dstream=fopen(plot3dfile_svz,"wb");
  if(plot3dstream==NULL){
    if(plot3dstream==NULL){
      fprintf(stderr,"*** Warning: The file %s could not be opened for writing\n",plot3dfile_svz);
    }
    fclose(PLOT3DFILE);
    return 0;
  }
  {
    int nx, ny, nz;
    FILE_SIZE len;
    int error, isotest;
    int j;
    int framesize;
    int k,kk;

    len=strlen(plot3d_file);
    CheckMemory;
    nx = plot3di->plot3d_mesh->ibar+1;
    ny = plot3di->plot3d_mesh->jbar+1;
    nz = plot3di->plot3d_mesh->kbar+1;
    ijkbar[0]=nx;
    ijkbar[1]=ny;
    ijkbar[2]=nz;
    isotest=0;
    framesize=nx*ny*nz;
    NewMemory((void **)&plot3dframe_data,5*framesize*sizeof(float));
    NewMemory((void **)&plot3dframe_compressed,1.1*5*framesize*sizeof(unsigned char));
    NewMemory((void **)&plot3dframe_uncompressed,5*framesize*sizeof(unsigned char));

    FORTgetplot3dq(plot3d_file, &nx, &ny, &nz, plot3dframe_data, &error, &isotest, len);
    kk=0;
    for(j=0;j<5;j++){
      float valmin, valmax;
      float dv;

      valmin=plot3dinfo[0].bounds[j].valmin;
      valmax=plot3dinfo[0].bounds[j].valmax;
      minmax[2*j]=valmin;
      minmax[2*j+1]=valmax;
      dv = valmax-valmin;
      if(dv==0.0)dv=1.0;
      for(k=0;k<framesize;k++){
        float val;
        int ival;

        val = plot3dframe_data[kk];
        if(val<valmin){
          ival=0;
        }
        else if(val>valmax){
          ival=255;
        }
        else{
          ival = 1 + 253*(val-valmin)/dv;
        }
        plot3dframe_uncompressed[kk++]=ival;
      }
    }
    compress_zlib(plot3dframe_compressed,&ncompressed_zlib,plot3dframe_uncompressed,5*framesize);
    sizeafter=16+ncompressed_zlib;
    sizebefore=5*framesize*sizeof(float);
  }

  fwrite(&one,4,1,plot3dstream);           // write out a 1 to determine "endianness" when file is read in later
  fwrite(&zero,4,1,plot3dstream);          // write out a zero now, then a one just before file is closed
  fwrite(&fileversion,4,1,plot3dstream);   // write out compressed fileversion in case file format changes later
  fwrite(&version_local,4,1,plot3dstream);       // fds plot3d file version
  fwrite(minmax,4,10,plot3dstream);        // write out min/max bounds for each plot3d file component
  fwrite(ijkbar,4,3,plot3dstream);         // write out ibar, jbar, kbar
  fwrite(plot3dframe_compressed,1,ncompressed_zlib,plot3dstream);  // write out compressed plot3d data

  //*** PLOT3D FILE FORMAT

  //*** FDS FORMAT (FORTRAN - each FORTRAN record has a 4 byte header and a 4 byte trailer surrounding the data)

  // nx, ny, nz
  // dum1, dum2, dum3, dum4
  // (((qq(i,j,k),i=1,nx),j=1,ny),k=1,nz)

  //*** ZLIB format (C - no extra bytes surrounding data)

  //*** header
  // endian
  // completion (0/1)
  // fileversion (compressed format)
  // version  (plot3df version)
  // (xmin(i),ymin(i),i=1,5)
  // nx, ny, nz


  //*** frame
  // compressed qq buffer

  FREEMEMORY(plot3dframe_data);
  FREEMEMORY(plot3dframe_compressed);
  FREEMEMORY(plot3dframe_uncompressed);

  fclose(PLOT3DFILE);
  FSEEK(plot3dstream,4,SEEK_SET);
  fwrite(&one,4,1,plot3dstream);  // write completion code
  fclose(plot3dstream);

  {
    char before_label[256],after_label[256];

    getfilesizelabel(sizebefore,before_label);
    getfilesizelabel(sizeafter,after_label);
#ifdef pp_THREAD
    LOCK_PRINT;
    PRINTF("\n%s\n  compressed from %s to %s (%4.1f%s reduction)\n\n",plot3di->file,before_label,after_label,(float)sizebefore/(float)sizeafter,GLOBx);
    UNLOCK_PRINT;
#else
    PRINTF("Sizes: original=%s, ",before_label);
    PRINTF("compressed=%s (%4.1f%s reduction)\n",after_label,(float)sizebefore/(float)sizeafter,GLOBx);
#endif
  }

  return 1;

}

/* ------------------ getplot3d ------------------------ */

plot3d *getplot3d(char *string){
  int i;
  plot3d *plot3di;

  for(i=0;i<nplot3dinfo;i++){
    plot3di = plot3dinfo + i;
    if(plot3di->dup==1)continue;
    if(strcmp(plot3di->labels[0].shortlabel,string)==0)return plot3di;
  }
  return NULL;
}

/* ------------------ compress_plot3ds ------------------------ */

void *compress_plot3ds(void *arg){
  int i, j;
  plot3d *plot3di, *pb;

  LOCK_PLOT3D;
  if(GLOBfirst_plot3d==1){
    GLOBfirst_plot3d=0;
    for(i=0;i<nplot3dinfo;i++){
      plot3di = plot3dinfo + i;
      if(GLOBautozip==1&&plot3di->autozip==0)continue;
      plot3di->count=0;
    }
    for(i=0;i<nplot3dinfo;i++){
      plot3di = plot3dinfo + i;
      if(GLOBautozip==1&&plot3di->autozip==0)continue;
      plot3di->doit=1;

      pb=plot3dinfo;
      for(j=0;j<5;j++){
        plot3di->bounds[j].setvalmin=pb->bounds[j].setvalmin;
        plot3di->bounds[j].setvalmax=pb->bounds[j].setvalmax;
        plot3di->bounds[j].valmin=pb->bounds[j].valmin;
        plot3di->bounds[j].valmax=pb->bounds[j].valmax;
      }
      pb->count++;
    }
  }
  UNLOCK_PLOT3D;

  // convert and compress files

  for(i=0;i<nplot3dinfo;i++){
    plot3di = plot3dinfo + i;
    if(GLOBautozip==1&&plot3di->autozip==0)continue;
    LOCK_PLOT3D;
    if(plot3di->inuse==1){
      UNLOCK_PLOT3D;
      continue;
    }
    plot3di->inuse=1;
    UNLOCK_PLOT3D;

    if(plot3di->doit==1){
      convert_plot3d(plot3di);
    }
    else{
      if(GLOBcleanfiles==0){
        PRINTF("%s not compressed\n",plot3di->file);
        PRINTF("  Min and Max for %s not set in .ini file\n",plot3di->labels[0].shortlabel);
      }
    }
  }
  return NULL;
}

/* ------------------ plot3ddup ------------------------ */

int plot3ddup(plot3d *plot3dj, int iplot3d){
  int i;
  plot3d *plot3di;

  for(i=0;i<iplot3d;i++){
    plot3di = plot3dinfo + i;
    if(plot3di->dup==1)continue;
    if(strcmp(plot3di->labels[0].shortlabel,plot3dj->labels[0].shortlabel)==0){
      plot3dj->dup=1;
      return 1;
    }
  }
  return 0;
}
