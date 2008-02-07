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
unsigned char *full_alphabuffer;

#define IJKNODE(i,j,k) ((i)+(j)*nx+(k)*nxy)

#ifdef pp_LIGHT
void light_smoke(smoke3d *smoke3di,unsigned char *full_lightingbuffer, float *val_buffer, unsigned char *alpha_buffer);
void set_lightfield(smoke3d *smoke3di,float xyz[3], float hrr);
#endif

#define FORTREAD(read) fseek(BOUNDARYFILE,4,SEEK_CUR);returncode=read;fseek(BOUNDARYFILE,4,SEEK_CUR);

/* ------------------ convert_3dsmoke ------------------------ */

void convert_3dsmoke(smoke3d *smoke3di){
  unsigned char *compressed_alphabuffer;
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
//    float xx;

  //  xx = 1.0 - (float)i/255.0;
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
  full_alphabuffer=NULL;
  compressed_alphabuffer=NULL;
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
      printf("  3dsmoke file, %s, could not be opened for output\n",smoke3dfile_svz);
    }
    if(smoke3dsizestream==NULL){
      printf("  3dsmoke size file, %s, could not be opened for output\n",smoke3dsizefile_svz);
    }
#ifdef pp_LIGHT
    if(make_lighting_file==1&&light3dstream==NULL){
      printf("  3dsmoke lighting file, %s, could not be opened for output\n",smoke3dfile_lvz);
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

  full_alphabuffer=NULL;
  NewMemory((void **)&full_alphabuffer,buffersize);
  compressed_alphabuffer=NULL;
  NewMemory((void **)&compressed_alphabuffer,buffersize);

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

    EGZ_FREAD(compressed_alphabuffer,ncompressed_rle,1,SMOKE3DFILE);

    count++;

    sizebefore+=12+ncompressed_rle;
#ifdef pp_LIGHT
    if(make_lighting_file==1)lightbefore+=nx*ny*nz;
#endif

    if(count%smoke3dzipstep!=0)continue;

    // uncompress frame data (from RLE format)

    nfull2=irle(compressed_alphabuffer, ncompressed_rle, full_alphabuffer);
    if(nfull!=nfull2){
      printf("  ***warning frame size expected: %i actual: %i\n",nfull,nfull2);
    }

#ifdef pp_LIGHT
    for(i=0;i<nfull2;i++){
      full_alphabuffer[i]=mapalpha[full_alphabuffer[i]];
    }
#endif

    // compress frame data (into ZLIB format)

    ncompressed_zlib=buffersize;
    returncode=compress(compressed_alphabuffer, &ncompressed_zlib, full_alphabuffer, nfull2);
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

      light_smoke(smoke3di,full_lightingbuffer,val_buffer, full_alphabuffer);

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
    if(ncompressed_zlib>0)fwrite(compressed_alphabuffer,1,ncompressed_zlib,smoke3dstream);
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
  }

}

#ifdef pp_LIGHT

/* ------------------ init_lightfield ------------------------ */

void init_lightfield(void){
  NewMemory((void **)&light_q_polar,NRAD*NTHETA*NPSI);
}

/* ------------------ update_lightfield ------------------------ */

void update_lightfield(smoke3d *smoke3di, unsigned char *lightingbuffer){
  int ilight;
  int nlight_q_rect;
  int i, j, k;

  if(smoke3di->smoke_mesh==NULL)return;
  nlight_q_rect = smoke3di->nx*smoke3di->ny*smoke3di->nz;
  if(smoke3di->light_q_rect==NULL){
    if(nlight_q_rect>0){
      NewMemory((void **)&smoke3di->light_q_rect,nlight_q_rect*sizeof(float));
    }
    else{
      return;
    }
  }
  for(i=0;i<nlight_q_rect;i++){
    smoke3di->light_q_rect[i]=0.0;
  }

  // accumulate hrr for each light

  for(ilight=0;ilight<nlightinfo;ilight++){
    lightdata *lighti;

    lighti = lightinfo + ilight;
    switch (lighti->type){
      float *xyz1, *xyz2;
      float dx,dy,dz,length;
      int npoint, nx, ny, nz;
      float dxx, dyy, dzz;
      float xyz[3];

      case 0:      // point
        set_lightfield(smoke3di,lighti->xyz1,lighti->q);
        break;
      case 1:      // line
        xyz1 = lighti->xyz1;
        xyz2 = lighti->xyz2;
        dx = xyz1[0]-xyz2[0];        
        dy = xyz1[1]-xyz2[1];        
        dz = xyz1[2]-xyz2[2];        
        length = sqrt(dx*dx+dy*dy+dz*dz);
        npoint = length/light_delta+1.5;
        if(npoint<2)npoint=2;
        for(i=0;i<npoint;i++){
          xyz[0] = ((float)(npoint-1-i)*xyz1[0] + (float)i*xyz2[0])/(float)(npoint-1);
          xyz[1] = ((float)(npoint-1-i)*xyz1[1] + (float)i*xyz2[1])/(float)(npoint-1);
          xyz[2] = ((float)(npoint-1-i)*xyz1[2] + (float)i*xyz2[2])/(float)(npoint-1);
          set_lightfield(smoke3di,xyz,lighti->q/(float)npoint);
        }
        break;
      case 2:      // region
        xyz1 = lighti->xyz1;
        xyz2 = lighti->xyz2;
        dx = abs(xyz1[0]-xyz2[0]);
        dy = abs(xyz1[1]-xyz2[1]);
        dz = abs(xyz1[2]-xyz2[2]);
        nx = dx/light_delta+1.5;
        ny = dy/light_delta+1.5;
        nz = dz/light_delta+1.5;
        dxx = 0.0;
        dyy = 0.0;
        dzz = 0.0;
        if(nx>1){
          dxx = (xyz2[0]-xyz1[0])/(nx-1);
        }
        if(ny>1){
          dyy = (xyz2[1]-xyz1[1])/(ny-1);
        }
        if(nz>1){
          dzz = (xyz2[2]-xyz1[2])/(nz-1);
        }
        for(k=0;k<nz;k++){
          xyz[2] = xyz1[2] + k*dzz;
          for(j=0;j<ny;j++){
            xyz[1] = xyz1[1] + j*dyy;
            for(i=0;i<nx;i++){
              xyz[0] = xyz1[0] + i*dxx;
              set_lightfield(smoke3di,xyz,lighti->q/(float)(nx*ny*nz));
            }
          }
        }
        break;
    }
  }

  // convert hrr field to colors
}

/* ------------------ getldist ------------------------ */

float getldist(float ldist,float *xyz, float x, float y, float z){
  float dx, dy, dz;
  float dist;

  dx = x-xyz[0];
  dy = y-xyz[1];
  dz = z-xyz[2];
  dist = sqrt(dx*dx+dy*dy+dz*dz);
  if(ldist>dist)dist=ldist;
  return dist;
}

/* ------------------ getalpha ------------------------ */

float getalpha(mesh *smoke_mesh, float *xyz2){
  int i, j, k;
  int nx, ny, nxy;
  int ialpha;

  if(xyz2[0]<smoke_mesh->xbar0||xyz2[0]>smoke_mesh->xbar)return 1.0;
  if(xyz2[1]<smoke_mesh->ybar0||xyz2[1]>smoke_mesh->ybar)return 1.0;
  if(xyz2[2]<smoke_mesh->zbar0||xyz2[2]>smoke_mesh->zbar)return 1.0;

  nx = smoke_mesh->ibar;
  ny = smoke_mesh->jbar;
  nxy = nx*ny;

  i = (xyz2[0]-smoke_mesh->xbar0)/(smoke_mesh->xbar-smoke_mesh->xbar0)*smoke_mesh->ibar;
  if(i<0)i=0;
  if(i>smoke_mesh->ibar)i=smoke_mesh->ibar;

  j = (xyz2[1]-smoke_mesh->ybar0)/(smoke_mesh->ybar-smoke_mesh->ybar0)*smoke_mesh->jbar;
  if(j<0)j=0;
  if(j>smoke_mesh->jbar)j=smoke_mesh->jbar;

  k = (xyz2[2]-smoke_mesh->zbar0)/(smoke_mesh->zbar-smoke_mesh->zbar0)*smoke_mesh->kbar;
  if(k<0)k=0;
  if(k>smoke_mesh->kbar)k=smoke_mesh->kbar;

  ialpha = full_alphabuffer[IJKNODE(i,j,k)];
  return (float)ialpha/255.0;
}

/* ------------------ atan3 ------------------------ */

float atan3(float dy,float dx){
  if(dx!=0.0)return atan(dy/dx);

  // dx is zero so atan(dy/dx) is PI/2 or -PI/2 depending on sign dy

  if(dy>0.0)return 2.0*atan(1.0);
  if(dy<0.0)return -2.0*atan(1.0);
  return 0.0;
}

/* ------------------ set_lightfield ------------------------ */

void set_lightfield(smoke3d *smoke3di,float xyz[3], float light_q_source){
  int i,j,k;
  mesh *smoke_mesh;
  float rads[NRAD+1], area[NRAD+1];
  float cos_theta[NTHETA+1], sin_theta[NTHETA+1];
  float cos_psi[NPSI+1], sin_psi[NPSI+1];
  float ldist;
  float PI, theta, psi;
  int ipsi, irad, itheta;
  int nrad, nradtheta;
  int inode;
  float xbar0, xbar;
  float ybar0, ybar;
  float zbar0, zbar;
  int ibar, jbar, kbar;
  int nx, ny, nxy;

  smoke_mesh = smoke3di->smoke_mesh;
  
  nrad=NRAD;
  nradtheta=NRAD*NTHETA;

  xbar0 = smoke_mesh->xbar0;
  xbar =  smoke_mesh->xbar;
  ybar0 = smoke_mesh->ybar0;
  ybar =  smoke_mesh->ybar;
  zbar0 = smoke_mesh->zbar0;
  zbar =  smoke_mesh->zbar;

  ibar = smoke_mesh->ibar;
  jbar = smoke_mesh->jbar;
  kbar = smoke_mesh->kbar;

  nx = smoke3di->nx;
  ny = smoke3di->ny;
  nxy = nx*ny;

  ldist=0.0;
  ldist=getldist(ldist,xyz,xbar0,ybar0,zbar0);
  ldist=getldist(ldist,xyz, xbar,ybar0,zbar0);
  ldist=getldist(ldist,xyz,xbar0, ybar,zbar0);
  ldist=getldist(ldist,xyz, xbar, ybar,zbar0);
  ldist=getldist(ldist,xyz,xbar0,ybar0, zbar);
  ldist=getldist(ldist,xyz, xbar,ybar0, zbar);
  ldist=getldist(ldist,xyz,xbar0, ybar, zbar);
  ldist=getldist(ldist,xyz, xbar, ybar, zbar);

  PI=4.0*atan(1.0);
  for(i=0;i<NRAD;i++){
    float rad;
    rad=(float)(i+1)*ldist/(float)NRAD;
    rads[i]=rad;
    area[i]=4.0*PI*rad*rad;
  }
  for(i=0;i<NTHETA+1;i++){
    theta = (float)i*2.0*PI/(float)NTHETA;
    cos_theta[i]=cos(theta);
    sin_theta[i]=sin(theta);
  }
  for(i=0;i<NPSI+1;i++){
    psi = (float)i*2.0*PI/(float)NPSI;
    cos_psi[i]=cos(psi);
    sin_psi[i]=sin(psi);
  }

  // set polar hrr field to zero

  for(i=0;i<NRAD*NTHETA*NPSI;i++){
    light_q_polar[i]=0.0;
  }

  // set center of field to hrr/area then
  //   set each successive shell using new_shell_hrrpua = old_shell_hrrpua*(r/(r+dr))^2 (1-alpha)

#define GETPOLARNODE(ipsi,itheta,irad) ((irad)+(itheta)*nrad+(ipsi)*nradtheta)
  inode=0;
  for(ipsi=0;ipsi<NPSI;ipsi++){
    for(itheta=0;itheta<NTHETA;itheta++){
      inode=GETPOLARNODE(ipsi,itheta,0);
      light_q_polar[inode]=light_q_source/area[0];
      for(irad=1;irad<NRAD;irad++){
        float xyz2[3];
        float alpha;

        if(light_q_polar[inode-1]==0.0){
          light_q_polar[inode]=0.0;
        }
        else{
          xyz2[0] = xyz[0] + rads[irad]*cos_theta[itheta]*cos_psi[ipsi];
          xyz2[1] = xyz[1] + rads[irad]*sin_theta[itheta]*cos_psi[ipsi];
          xyz2[2] = xyz[2] + rads[irad]*sin_psi[ipsi];
          alpha = getalpha(smoke_mesh,xyz2);
          light_q_polar[inode]=light_q_polar[inode-1]*area[irad-1]/area[irad]*(1.0-alpha);
        }
        inode++;
      }
    }
  }

  // add polar field to rectangular field
  
  for(k=0;k<kbar;k++){
    float dx, dy, dz;
    float r, theta, psi;
    float xy_length;
    int irad, ipsi, itheta;
    int polarnode, ijknode;

    dz = (zbar0*(float)(kbar-1-k) + (float)k*zbar)/(float)(kbar-1)-xyz[2];
    for(j=0;j<smoke_mesh->jbar;j++){
      dy = (ybar0*(float)(jbar-1-j) + (float)j*ybar)/(float)(jbar-1)-xyz[1];
      for(i=0;i<smoke_mesh->ibar;i++){
        dx = (xbar0*(float)(ibar-1-i) + (float)i*xbar)/(float)(ibar-1)-xyz[0];
        r = sqrt(dx*dx+dy*dy+dz*dz);
        irad=(NRAD-1)*(r/ldist);
        if(irad<0)irad=0;
        if(irad>NRAD-1)irad=NRAD-1;

        xy_length = sqrt(dx*dx+dy*dy);
        psi = atan3(dz,xy_length)+PI/2.0;
        ipsi = (NPSI-1)*psi/PI;
        if(ipsi<0)ipsi=0;
        if(ipsi>NPSI-1)ipsi=NPSI-1;

        theta = 0.0;
        if(dx!=0.0||dy!=0.0){
          theta = atan2(dy,dx);
        }
        theta+=PI;
        itheta = theta*(NTHETA-1)*theta/(2.0*PI);
        if(itheta<0)itheta=0;
        if(itheta>NTHETA-1)itheta=NTHETA-1;

        polarnode = GETPOLARNODE(ipsi,itheta,irad);
        ijknode = IJKNODE(i,j,k);
        smoke3di->light_q_rect[ijknode]+=light_q_polar[polarnode];
      }
    }
  }
}

/* ------------------ light_smoke ------------------------ */

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

  // light from below

  // set planes above k=0 to 0.0

  for(k=1;k<nz;k++){
    for(j=0;j<ny;j++){
      ijk = IJKNODE(0,j,k);
      for(i=0;i<nx;i++){
        val_buffer[ijk]=0.0;
        ijk++;
      }
    }
  }

  // light from below

  for(k=1;k<nz;k++){
    for(j=0;j<ny;j++){
      ijk = IJKNODE(0,j,k);
      for(i=0;i<nx;i++){
        ijkm1 = ijk - nxy;
        factor = (255.0-(float)alpha_buffer[ijkm1])/255.0;
        val = factor*val_buffer[ijkm1];
        if(val>val_buffer[ijk])val_buffer[ijk]=val;
        ijk++;
      }
    }
  }

  for(i=0;i<n_lightingbuffer;i++){
    lightingbuffer[i]=255*albedo*val_buffer[i];
  }

}
#endif
