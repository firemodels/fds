// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "zlib.h"
#include "egz_stdio.h"
#include "svzip.h"
#include "MALLOC.h"

// svn revision character string
char CNV3dsmoke_revision[]="$Revision$";


#define FORTREAD(read) fseek(BOUNDARYFILE,4,SEEK_CUR);returncode=read;fseek(BOUNDARYFILE,4,SEEK_CUR);

/* ------------------ convert_3dsmoke ------------------------ */

void convert_3dsmoke(smoke3d *smoke3di){
  unsigned char *compressed_alphabuffer;
  FILE *smoke3dstream=NULL,*smoke3dsizestream=NULL;
  EGZ_FILE *SMOKE3DFILE=NULL;
  char smoke3dfile_svz[1024], smoke3dsizefile_svz[1024];
  unsigned char *full_alphabuffer;
  int nxyz[9];
  int nx, ny, nz;
  int version;
  int buffersize;
  unsigned int sizebefore, sizeafter;
  int count;
  float time;
  int nchars[2];
  int nfull,nfull2;
  int ncompressed_rle;
  uLongf ncompressed_zlib;
  int returncode;
  char pp[2];
  char xxx[2];
  int percent_done;
  int percent_next=10;
  long data_loc;
  char *smoke3dfile;
  float time_max;
 
  smoke3dfile=smoke3di->file;


//int compress (Bytef *dest,   uLongf *destLen, const Bytef *source, uLong sourceLen);
//int uncompress (Bytef *dest,   uLongf *destLen, const Bytef *source, uLong sourceLen);

  if(cleanfiles==0)printf("Compressing 3D smokefile %s\n",smoke3dfile);

  strcpy(pp,"%");
  strcpy(xxx,"X");
  full_alphabuffer=NULL;
  compressed_alphabuffer=NULL;
  
  if(getfileinfo(smoke3dfile,NULL,NULL)!=0){
    printf("  %s does not exist\n",smoke3dfile);
    return;
  }
  SMOKE3DFILE=EGZ_FOPEN(smoke3dfile,"rb",0,2);
  if(SMOKE3DFILE==NULL){
    printf("  %s could not be opened\n",smoke3dfile);
    return;
  }

  // name 3d smoke flie

  if(destdir!=NULL){
    strcpy(smoke3dfile_svz,destdir);
    strcat(smoke3dfile_svz,smoke3di->filebase);
  }
  else{
    strcpy(smoke3dfile_svz,smoke3di->file);
  }
  strcat(smoke3dfile_svz,".svz");

  // name size file

  if(destdir!=NULL){
    strcpy(smoke3dsizefile_svz,destdir);
    strcat(smoke3dsizefile_svz,smoke3di->filebase);
  }
  else{
    strcpy(smoke3dsizefile_svz,smoke3di->file);
  }
  strcat(smoke3dsizefile_svz,".szz");

  // remove files if clean option is set

  if(cleanfiles==1){
    smoke3dstream=fopen(smoke3dfile_svz,"rb");
    if(smoke3dstream!=NULL){
      fclose(smoke3dstream);
      printf("  Removing %s.\n",smoke3dfile_svz);
      UNLINK(smoke3dfile_svz);
      filesremoved++;
    }
    smoke3dsizestream=fopen(smoke3dsizefile_svz,"r");
    if(smoke3dsizestream!=NULL){
      fclose(smoke3dsizestream);
      printf("  Removing %s.\n",smoke3dsizefile_svz);
      UNLINK(smoke3dsizefile_svz);
      filesremoved++;
    }
    return;
  }

  if(overwrite_s==0){
    smoke3dstream=fopen(smoke3dfile_svz,"rb");
    if(smoke3dstream!=NULL){
      fclose(smoke3dstream);
      printf("  %s exists.\n",smoke3dfile_svz);
      printf("     Use the -f option to overwrite smokezip compressed files\n");
      return;
    }
  }

  smoke3dsizestream=fopen(smoke3dsizefile_svz,"w");
  smoke3dstream=fopen(smoke3dfile_svz,"wb");
  if(smoke3dstream==NULL||smoke3dsizestream==NULL
    ){
    if(smoke3dstream==NULL){
      printf("  3dsmoke file, %s, could not be opened for output\n",smoke3dfile_svz);
    }
    if(smoke3dsizestream==NULL){
      printf("  3dsmoke size file, %s, could not be opened for output\n",smoke3dsizefile_svz);
    }
    if(smoke3dsizestream!=NULL)fclose(smoke3dsizestream);
    if(smoke3dstream!=NULL)fclose(smoke3dstream);
    EGZ_FCLOSE(SMOKE3DFILE);
    return;
  }

  EGZ_FREAD(nxyz,4,8,SMOKE3DFILE);

  nxyz[0] = 1;
  version = nxyz[1];
  if(version==1){
    printf("  already compressed\n");
    EGZ_FCLOSE(SMOKE3DFILE);
    fclose(smoke3dstream);
    fclose(smoke3dsizestream);
    return;
  }

  version=1;

  nxyz[1]=version;
  fwrite(nxyz,4,8,smoke3dstream);
  fprintf(smoke3dsizestream,"%i\n",version);

  nx = nxyz[3]-nxyz[2]+1;
  ny = nxyz[5]-nxyz[4]+1;
  nz = nxyz[7]-nxyz[6]+1;
  buffersize=1.01*nx*ny*nz+600;
  smoke3di->nx=nx;
  smoke3di->ny=ny;
  smoke3di->nz=nz;
  smoke3di->ncompressed_lighting_zlib=buffersize;

  full_alphabuffer=NULL;
  NewMemory((void **)&full_alphabuffer,buffersize);
  compressed_alphabuffer=NULL;
  NewMemory((void **)&compressed_alphabuffer,buffersize);

  count=-1;
  sizebefore=8;
  sizeafter=8;
  printf("  Compressing: ");
  time_max=-1000000.0;
  for(;;){
    EGZ_FREAD(&time,4,1,SMOKE3DFILE);
    if(EGZ_FEOF(SMOKE3DFILE)!=0)break;

    EGZ_FREAD(nchars,4,2,SMOKE3DFILE);
    nfull=nchars[0];
    ncompressed_rle=nchars[1];

    // read compressed frame

    EGZ_FREAD(compressed_alphabuffer,ncompressed_rle,1,SMOKE3DFILE);

    if(time<time_max)continue;
    count++;

    sizebefore+=12+ncompressed_rle;

    if(count%smoke3dzipstep!=0)continue;
    time_max=time;

    // uncompress frame data (from RLE format)

    nfull2=irle(compressed_alphabuffer, ncompressed_rle, full_alphabuffer);
    CheckMemory;
    if(nfull!=nfull2){
      printf("  ***warning frame size expected=%i frame size found=%i\n",nfull,nfull2);
    }

    // compress frame data (into ZLIB format)

    ncompressed_zlib=buffersize;
    returncode=compress(compressed_alphabuffer, &ncompressed_zlib, full_alphabuffer, nfull2);
    CheckMemory;
    if(returncode!=0){
      printf("  ***warning zlib compressor failed - frame %f\n",time);
    }

    data_loc=EGZ_FTELL(SMOKE3DFILE);
    percent_done=100.0*(float)data_loc/(float)smoke3di->filesize;
    if(percent_done>percent_next){
        printf(" %i%s",percent_next,pp);
        LOCK_SMOKE;
        fflush(stdout);
        UNLOCK_SMOKE;
      percent_next+=10;
    }

    // write out new entries in the size (sz) file

    nchars[0]=nfull2;
    nchars[1]=ncompressed_zlib;
    fwrite(&time,4,1,smoke3dstream);
    fwrite(nchars,4,2,smoke3dstream);
    if(ncompressed_zlib>0)fwrite(compressed_alphabuffer,1,ncompressed_zlib,smoke3dstream);
    sizeafter+=12+ncompressed_zlib;

    fprintf(smoke3dsizestream,"%f %i %i %i\n",time,nfull,ncompressed_rle,(int)ncompressed_zlib);
  }
  printf(" 100%s completed\n",pp);
  printf("  records=%i, ",count);
  {
    char before_label[256],after_label[256];
  
    getfilesizelabel(sizebefore,before_label);
    getfilesizelabel(sizeafter,after_label);

    printf("Sizes: original=%s, ",before_label);

    printf("compressed=%s (%4.1f%s reduction)\n",after_label,(float)sizebefore/(float)sizeafter,xxx);
    LOCK_SMOKE;
    fflush(stdout);
    UNLOCK_SMOKE;
  }

  // close files and free buffers

  EGZ_FCLOSE(SMOKE3DFILE);

  fclose(smoke3dstream);
  fclose(smoke3dsizestream);
  FREEMEMORY(full_alphabuffer);
  FREEMEMORY(compressed_alphabuffer);
}

/* ------------------ convert_smoke3ds ------------------------ */

void compress_smoke3ds(void){
  int i;
  smoke3d *smoke3di;

  printf("\n");
  for(i=0;i<nsmoke3d_files;i++){
    smoke3di = smoke3dinfo + i;
    if(autozip==1&&smoke3di->autozip==0)continue;
    convert_3dsmoke(smoke3di);
    CheckMemory;
  }
}
#ifdef pp_THREAD
/* ------------------ convert_smoke3ds ------------------------ */

void MT_compress_smoke3ds(void){
  int i;
  smoke3d *smoke3di;
  pthread_t *thread_ids;

  if(nsmoke3d_files<=0)return;

  CheckMemory;
  NewMemory((void **)&thread_ids,nsmoke3d_files*sizeof(pthread_t));
  CheckMemory;

  for(i=0;i<nsmoke3d_files;i++){
    smoke3di = smoke3dinfo + i;
    if(autozip==1&&smoke3di->autozip==0)continue;
    pthread_create(&thread_ids[i],NULL,MT_convert_3dsmoke,(void *)(smoke3di));
  }

  printf("\n");
  for(i=0;i<nsmoke3d_files;i++){
    smoke3di = smoke3dinfo + i;
    if(autozip==1&&smoke3di->autozip==0)continue;
    pthread_join(thread_ids[i],NULL);
  }
  FREEMEMORY(thread_ids);
}
#endif
