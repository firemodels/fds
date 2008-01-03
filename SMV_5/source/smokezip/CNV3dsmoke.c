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

#define ijknode(i,j,k) ((i)+(j)*nx+(k)*nxy)

#ifdef pp_LIGHT
void light_smoke(smoke3d *smoke3di,unsigned char *full_lightingbufer, float *val_buffer, unsigned char *alpha_buffer);
#endif

#define FORTREAD(read) fseek(BOUNDARYFILE,4,SEEK_CUR);returncode=read;fseek(BOUNDARYFILE,4,SEEK_CUR);

/* ------------------ convert_3dsmoke ------------------------ */

void convert_3dsmoke(smoke3d *smoke3di){
  unsigned char *full_smokebuffer,*compressed_smokebuffer;
  FILE *smoke3dstream=NULL,*smoke3dsizestream=NULL;
  EGZ_FILE *SMOKE3DFILE=NULL;
  char smoke3dfile_svz[1024], smoke3dsizefile_svz[1024];
#ifdef pp_LIGHT
  unsigned char *full_lightingbuffer,*compressed_lightingbuffer;
  float *val_buffer;
  char smoke3dfile_lvz[1024];
  FILE *light3dstream=NULL;
#endif
  int nxyz[9];
  int nx, ny, nz;
  int version;
  int buffersize;
  unsigned int sizebefore, sizeafter;
#ifdef pp_LIGHT
  unsigned int lightbefore, lightafter;
#endif
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
#ifdef pp_LIGHT
  int light_info[2]={1,0};
  unsigned char mapalpha[256];
  int i;
#endif
 
#ifdef pp_LIGHT
  for(i=0;i<256;i++){
    float xx;

    xx = 1.0 - (float)i/255.0;
    mapalpha[i]=i;
//    mapalpha[i]=255*(1.0-pow(xx,1.0-albedo));
  }
#endif

  smoke3dfile=smoke3di->file;


//int compress (Bytef *dest,   uLongf *destLen, const Bytef *source, uLong sourceLen);
//int uncompress (Bytef *dest,   uLongf *destLen, const Bytef *source, uLong sourceLen);

  if(cleanfiles==0)printf("Compressing 3D smokefile %s\n",smoke3dfile);

  strcpy(pp,"%");
  strcpy(xxx,"X");
  full_smokebuffer=NULL;
  compressed_smokebuffer=NULL;
#ifdef pp_LIGHT
  full_lightingbuffer=NULL;
  val_buffer=NULL;
  compressed_lightingbuffer=NULL;
#endif
  
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

#ifdef pp_LIGHT

  // name lighting file

  if(destdir!=NULL){
    strcpy(smoke3dfile_lvz,destdir);
    strcat(smoke3dfile_lvz,smoke3di->filebase);
  }
  else{
    strcpy(smoke3dfile_lvz,smoke3di->file);
  }
  strcat(smoke3dfile_lvz,".lvz");
#endif

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
      unlink(smoke3dfile_svz);
      filesremoved++;
    }
    smoke3dsizestream=fopen(smoke3dsizefile_svz,"r");
    if(smoke3dsizestream!=NULL){
      fclose(smoke3dsizestream);
      printf("  Removing %s.\n",smoke3dsizefile_svz);
      unlink(smoke3dsizefile_svz);
      filesremoved++;
    }
#ifdef pp_LIGHT
    light3dstream=fopen(smoke3dfile_lvz,"rb");
    if(light3dstream!=NULL){
      fclose(light3dstream);
      printf("  Removing %s.\n",smoke3dfile_lvz);
      unlink(smoke3dfile_lvz);
      filesremoved++;
    }
#endif
    return;
  }

  if(overwrite_s==0){
    smoke3dstream=fopen(smoke3dfile_svz,"rb");
    if(smoke3dstream!=NULL){
      fclose(smoke3dstream);
      printf("  %s exists.\n",smoke3dfile_svz);
      printf("     Use the -f option to overwrite boundary or 3d smoke files\n");
      return;
    }
//  smoke3dsizestream=fopen(smoke3dsizefile_svz,"r");
//    if(smoke3dsizestream!=NULL){
//      fclose(smoke3dsizestream);
//      printf("  %s exists.\n  Use the -f option if you wish to overwrite it.\n",smoke3dsizefile_svz);
//      return;
//    }
  }

  smoke3dsizestream=fopen(smoke3dsizefile_svz,"w");
  smoke3dstream=fopen(smoke3dfile_svz,"wb");
#ifdef pp_LIGHT
  if(make_lighting_file==1)light3dstream=fopen(smoke3dfile_lvz,"wb");
#endif
  if(smoke3dstream==NULL||smoke3dsizestream==NULL
#ifdef pp_LIGHT
    ||(make_lighting_file==1&&light3dstream==NULL)
#endif
    ){
    if(smoke3dstream==NULL){
      printf("  %s could not be opened for output\n",smoke3dfile_svz);
    }
    if(smoke3dsizestream==NULL){
      printf("  %s could not be opened for output\n",smoke3dsizefile_svz);
    }
#ifdef pp_LIGHT
    if(make_lighting_file==1&&light3dstream==NULL){
      printf("  %s could not be opened for output\n",smoke3dfile_lvz);
    }
#endif
    if(smoke3dsizestream!=NULL)fclose(smoke3dsizestream);
    if(smoke3dstream!=NULL)fclose(smoke3dstream);
#ifdef pp_LIGHT
    if(make_lighting_file==1&&light3dstream!=NULL)fclose(light3dstream);
#endif
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

#ifdef pp_LIGHT
  if(make_lighting_file==1){
    light_info[0]=1;
    light_info[1]=1;
    fwrite(light_info,4,2,light3dstream);
  }
#endif

  nx = nxyz[3]-nxyz[2]+1;
  ny = nxyz[5]-nxyz[4]+1;
  nz = nxyz[7]-nxyz[6]+1;
  buffersize=1.01*nx*ny*nz+600;
  smoke3di->nx=nx;
  smoke3di->ny=ny;
  smoke3di->nz=nz;
  smoke3di->ncompressed_lighting_zlib=buffersize;

  full_smokebuffer=NULL;
  NewMemory((void **)&full_smokebuffer,buffersize);
  compressed_smokebuffer=NULL;
  NewMemory((void **)&compressed_smokebuffer,buffersize);

#ifdef pp_LIGHT
  if(make_lighting_file==1){
    full_lightingbuffer=NULL;
    NewMemory((void **)&full_lightingbuffer,buffersize);
    val_buffer=NULL;
    NewMemory((void **)&val_buffer,buffersize*sizeof(float));
    compressed_lightingbuffer=NULL;
    NewMemory((void **)&compressed_lightingbuffer,buffersize);
    smoke3di->compressed_lightingbuffer=compressed_lightingbuffer;
  }
#endif

  count=-1;
  sizebefore=8;
  sizeafter=8;
#ifdef pp_LIGHT
  lightbefore=0;
  lightafter=0;
#endif
  printf("  Compressing: ");
  for(;;){
    EGZ_FREAD(&time,4,1,SMOKE3DFILE);
    if(EGZ_FEOF(SMOKE3DFILE)!=0)break;

    EGZ_FREAD(nchars,4,2,SMOKE3DFILE);
    nfull=nchars[0];
    ncompressed_rle=nchars[1];

    // read compressed frame

    EGZ_FREAD(compressed_smokebuffer,ncompressed_rle,1,SMOKE3DFILE);

    count++;

    sizebefore+=12+ncompressed_rle;
#ifdef pp_LIGHT
    if(make_lighting_file==1)lightbefore+=nx*ny*nz;
#endif

    if(count%smoke3dzipstep!=0)continue;

    // uncompress frame data (from RLE format)

    nfull2=irle(compressed_smokebuffer, ncompressed_rle, full_smokebuffer);
    if(nfull!=nfull2){
      printf("  ***warning frame size expected: %i actual: %i\n",nfull,nfull2);
    }

#ifdef pp_LIGHT
    for(i=0;i<nfull2;i++){
      full_smokebuffer[i]=mapalpha[full_smokebuffer[i]];
    }
#endif

    // compress frame data (into ZLIB format)

    ncompressed_zlib=buffersize;
    returncode=compress(compressed_smokebuffer, &ncompressed_zlib, full_smokebuffer, nfull2);
    if(returncode!=0){
      printf("  ***warning zlib compressor failed - frame %f\n",time);
    }

    data_loc=EGZ_FTELL(SMOKE3DFILE);
    percent_done=100.0*(float)data_loc/(float)smoke3di->filesize;
    if(percent_done>percent_next){
      printf(" %i%s",percent_next,pp);
      fflush(stdout);
      percent_next+=10;
    }

#ifdef pp_LIGHT
    // compress frame of lighting data (into ZLIB format)
    if(make_lighting_file==1){
      int return_code;

      light_smoke(smoke3di,full_lightingbuffer,val_buffer, full_smokebuffer);

      buffersize=1.01*nx*ny*nz+600;
      smoke3di->ncompressed_lighting_zlib=buffersize;
      return_code=compress(smoke3di->compressed_lightingbuffer, &smoke3di->ncompressed_lighting_zlib, full_lightingbuffer, nx*ny*nz);
      if(return_code!=0){
        printf("  ***warning zlib compressor failed - frame %f\n",time);
      }
      lightafter+=smoke3di->ncompressed_lighting_zlib;
    }

#endif

    // write out new entries in the size (sz) file

    nchars[0]=nfull2;
    nchars[1]=ncompressed_zlib;
    fwrite(&time,4,1,smoke3dstream);
    fwrite(nchars,4,2,smoke3dstream);
    if(ncompressed_zlib>0)fwrite(compressed_smokebuffer,1,ncompressed_zlib,smoke3dstream);
    sizeafter+=12+ncompressed_zlib;

#ifdef pp_LIGHT
    if(make_lighting_file==1){
      fwrite(&time,4,1,light3dstream);
      fwrite(&smoke3di->ncompressed_lighting_zlib,4,1,light3dstream);
      if(smoke3di->ncompressed_lighting_zlib>0)fwrite(compressed_lightingbuffer,1,smoke3di->ncompressed_lighting_zlib,light3dstream);
    }
    fprintf(smoke3dsizestream,"%f %i %i %i",time,nfull,ncompressed_rle,ncompressed_zlib);
    if(make_lighting_file==1)fprintf(smoke3dsizestream," %i",smoke3di->ncompressed_lighting_zlib);
    fprintf(smoke3dsizestream,"\n");
#else
    fprintf(smoke3dsizestream,"%f %i %i %i\n",time,nfull,ncompressed_rle,(int)ncompressed_zlib);
#endif
  }
  printf(" 100%s completed\n",pp);
  printf("  records=%i, ",count);
  {
    char before_label[256],after_label[256];
  
    getfilesizelabel(sizebefore,before_label);
    getfilesizelabel(sizeafter,after_label);

    printf("Sizes: original=%s, ",before_label);

    printf("compressed=%s (%4.1f%s reduction)\n",after_label,(float)sizebefore/(float)sizeafter,xxx);
#ifdef pp_LIGHT
    if(make_lighting_file==1){
      getfilesizelabel(lightbefore,before_label);
      getfilesizelabel(lightafter,after_label);

      printf("     Sizes(lighting): original=%s, ",before_label);

      printf("compressed=%s (%4.1f%s reduction)\n",after_label,(float)lightbefore/(float)lightafter,xxx);
    }
#endif
  }


  // close files and free buffers

  EGZ_FCLOSE(SMOKE3DFILE);
  fclose(smoke3dstream);
  fclose(smoke3dsizestream);
#ifdef pp_LIGHT
  if(make_lighting_file==1)fclose(light3dstream);
#endif
  FREEMEMORY(full_smokebuffer);
  FREEMEMORY(compressed_smokebuffer);
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
  }

}

#ifdef pp_LIGHT
void light_smoke(smoke3d *smoke3di, unsigned char *lightingbuffer, float *val_buffer, unsigned char *alpha_buffer){
  int n_lightingbuffer;
  int nx, ny, nz, nxy;
  int i, j, k, ijk, ijkm1;
  float factor;
  float val;

  nx = smoke3di->nx;
  ny = smoke3di->ny;
  nz = smoke3di->nz;
  nxy = nx*ny;

  // set everything to 1.0

  n_lightingbuffer=nx*ny*nz;
  for(i=0;i<n_lightingbuffer;i++){  // initialize data
    lightingbuffer[i]=0;
    val_buffer[i]=1.0;
  }

  // set interior to  0.0 (hence outside boundary is 1.0)

  for(k=1;k<nz-1;k++){
    for(j=1;j<ny-1;j++){
      for(i=1;i<nx-1;i++){
        ijk = ijknode(i,j,k);
        val_buffer[ijk]=0.0;
      }
    }
  }

  // light from left

  for(k=0;k<nz;k++){
    for(j=0;j<ny;j++){
      for(i=1;i<nx;i++){
        ijk = ijknode(i,j,k);
        ijkm1 = ijknode(i-1,j,k);
        factor = (255.0-(float)alpha_buffer[ijkm1])/255.0;
        val=val_buffer[ijkm1]*factor;
        if(val>val_buffer[ijk])val_buffer[ijk]=val;
      }
    }
  }

  // light from right

  for(k=0;k<nz;k++){
    for(j=0;j<ny;j++){
      for(i=nx-1;i>0;i--){
        ijk = ijknode(i,j,k);
        ijkm1 = ijknode(i-1,j,k);
        factor = (255.0-(float)alpha_buffer[ijk])/255.0;
        val=val_buffer[ijk]*factor;
        if(val>val_buffer[ijkm1])val_buffer[ijkm1]=val;
      }
    }
  }

  // light from front

  for(k=0;k<nz;k++){
    for(j=1;j<ny;j++){
      for(i=0;i<nx;i++){
        ijk = ijknode(i,j,k);
        ijkm1 = ijknode(i,j-1,k);
        factor = (255.0-(float)alpha_buffer[ijkm1])/255.0;
        val = factor*val_buffer[ijkm1];
        if(val>val_buffer[ijk])val_buffer[ijk]=val;
      }
    }
  }

  // light from back

  for(k=0;k<nz;k++){
    for(j=ny-1;j>0;j--){
      for(i=0;i<nx;i++){
        ijk = ijknode(i,j,k);
        ijkm1 = ijknode(i,j-1,k);
        factor = (255.0-(float)alpha_buffer[ijk])/255.0;
        val = factor*val_buffer[ijk];
        if(val>val_buffer[ijkm1])val_buffer[ijkm1]=val;
      }
    }
  }


  // light from above

  for(k=nz-1;k>0;k--){
    for(j=0;j<ny;j++){
      for(i=1;i<nx;i++){
        ijk = ijknode(i,j,k);
        ijkm1 = ijknode(i,j,k-1);
        factor = (255.0-(float)alpha_buffer[ijk])/255.0;
        val = factor*val_buffer[ijk];
        if(val>val_buffer[ijkm1])val_buffer[ijkm1]=val;
      }
    }
  }

  for(i=0;i<n_lightingbuffer;i++){
    lightingbuffer[i]=255*albedo*val_buffer[i];
  }

}
#endif
