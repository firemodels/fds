#include "options.h"
#include "glew.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include GLUT_H
#include <float.h>

#include "update.h"
#include "smokeviewvars.h"
#include "IOvolsmoke.h"
#include "compress.h"

#define SKIP FSEEK( SMOKE3DFILE, fortran_skip, SEEK_CUR)

int cull_count=0;

#ifdef pp_CULL
int cullplane_compare( const void *arg1, const void *arg2 );
#endif
void Smoke3d_CB(int var);

//              alphaf_out[n]=adjustalpha(ALPHAIN, xyzeyeorig, xp, ASPECTRATIO, NORM, NORMTYPE);

// -------------------------- ADJUSTALPHA ----------------------------------

#define ADJUSTALPHA(ALPHAIN,ASPECTRATIO,NORM,NORMTYPE) \
            alphaf_out[n]=0;\
            if(ALPHAIN==0)continue;\
            if((adjustalphaflag==2||adjustalphaflag==3)&&iblank_smoke3d!=NULL&&iblank_smoke3d[n]==SOLID)continue;\
            if(adjustalphaflag==2){\
              alphaf_out[n]=ALPHAIN;\
            }\
            else{\
              alphaf_out[n]=adjustalpha(ALPHAIN, ASPECTRATIO);\
            }

// -------------------------- DRAWVERTEX ----------------------------------
#define DRAWVERTEX(XX,YY,ZZ)        \
if(show_smoketest==0){\
  if(show_smoke_lighting==1&&have_light==1){\
    light_value[0]=light_in[n11]*smoke_albedo; \
    light_value[1]=light_in[n12]*smoke_albedo; \
    light_value[2]=light_in[n22]*smoke_albedo; \
    light_value[3]=light_in[n21]*smoke_albedo; \
  }\
  value[0]=alphaf_ptr[n11]; \
  value[1]=alphaf_ptr[n12]; \
  value[2]=alphaf_ptr[n22]; \
  value[3]=alphaf_ptr[n21]; \
  if(value[0]==0&&value[1]==0&&value[2]==0&&value[3]==0)continue;\
  ivalue[0]=n11<<2;  \
  ivalue[1]=n12<<2;  \
  ivalue[2]=n22<<2;  \
  ivalue[3]=n21<<2;  \
  if(ABS(value[0]-value[2])<ABS(value[1]-value[3])){     \
    xyzindex=xyzindex1;                                  \
  }                                                      \
  else{                                                  \
    xyzindex=xyzindex2;                                  \
  }                                                      \
  for(node=0;node<6;node++){                             \
    unsigned char alphabyte;\
    float alphaval;\
    int mm;\
    mm = xyzindex[node];                                 \
    alphabyte = value[mm];                               \
    if(skip_global==2){\
      alphaval=alphabyte/255.0; \
      alphaval=alphaval*(2.0-alphaval);                  \
      alphabyte=alphaval*255.0;\
    }\
    else if(skip_global==3){\
      alphaval=alphabyte/255.0;                              \
      alphaval = alphaval*(3.0-alphaval*(3.0-alphaval));\
      alphabyte = 255*alphaval; \
    }\
    if(show_smoke_lighting==1&&have_light==1){\
      glColor4ub(light_value[mm],light_value[mm],light_value[mm],alphabyte);                                \
    }\
    else{\
      colorptr=mergecolorptr+ivalue[mm];\
      colorptr[3]=alphabyte;                                   \
      glColor4ubv(colorptr);                                \
    }\
    glVertex3f(XX,YY,ZZ);                                \
  }\
}\
else{\
  for(node=0;node<6;node++){                             \
    int mm;\
    mm = xyzindex1[node];                                 \
    glColor4ub(0,0,0,(unsigned char)smoke_alpha);\
    glVertex3f(XX,YY,ZZ);                                \
  }\
}

// -------------------------- DRAWVERTEXTERRAIN ----------------------------------

#define DRAWVERTEXTERRAIN(XX,YY,ZZ)        \
if(show_smoketest==0){\
  value[0]=alphaf_ptr[n11]; \
  value[1]=alphaf_ptr[n12]; \
  value[2]=alphaf_ptr[n22]; \
  value[3]=alphaf_ptr[n21]; \
  if(value[0]==0&&value[1]==0&&value[2]==0&&value[3]==0)continue;\
  z_offset[XXX]=znode_offset[m11];\
  z_offset[YYY]=znode_offset[m12];\
  z_offset[ZZZ]=znode_offset[m22];\
  z_offset[3]=znode_offset[m21];\
  ivalue[0]=n11<<2;  \
  ivalue[1]=n12<<2;  \
  ivalue[2]=n22<<2;  \
  ivalue[3]=n21<<2;  \
  if(ABS(value[0]-value[2])<ABS(value[1]-value[3])){     \
    xyzindex=xyzindex1;                                  \
  }                                                      \
  else{                                                  \
    xyzindex=xyzindex2;                                  \
  }                                                      \
  for(node=0;node<6;node++){                             \
    unsigned char alphabyte;\
    float alphaval;\
    int mm;\
    mm = xyzindex[node];                                 \
    alphabyte = value[mm];                               \
    if(skip_global==2){\
      alphaval=alphabyte/255.0; \
      alphaval=alphaval*(2.0-alphaval);                  \
      alphabyte=alphaval*255.0;\
    }\
    else if(skip_global==3){\
      alphaval=alphabyte/255.0;                              \
      alphaval = alphaval*(3.0-alphaval*(3.0-alphaval));\
      alphabyte = 255*alphaval; \
    }\
    colorptr=mergecolorptr+ivalue[mm];\
    colorptr[3]=alphabyte;                                   \
    glColor4ubv(colorptr);                                \
    glVertex3f(XX,YY,ZZ+z_offset[mm]);                                \
  }\
}\
else{\
  z_offset[XXX]=znode_offset[m11];\
  z_offset[YYY]=znode_offset[m12];\
  z_offset[ZZZ]=znode_offset[m22];\
  z_offset[3]=znode_offset[m21];\
  for(node=0;node<6;node++){                             \
    int mm;\
    mm = xyzindex1[node];                                 \
    glColor4ub(0,0,0,(unsigned char)smoke_alpha);\
    glVertex3f(XX,YY,ZZ+z_offset[mm]);                                \
  }\
}

// -------------------------- DRAWVERTEXGPU ----------------------------------

#define DRAWVERTEXGPU(XX,YY,ZZ) \
  value[0]=alphaf_in[n11];\
  value[1]=alphaf_in[n12];\
  value[2]=alphaf_in[n22];\
  value[3]=alphaf_in[n21];\
  if(firecolor==NULL&&value[0]==0&&value[1]==0&&value[2]==0&&value[3]==0)continue;\
  if((adjustalphaflag==2||adjustalphaflag==3)&&iblank_smoke3d!=NULL){\
    if(iblank_smoke3d[n11]==SOLID)value[0]=0;\
    if(iblank_smoke3d[n12]==SOLID)value[1]=0;\
    if(iblank_smoke3d[n22]==SOLID)value[2]=0;\
    if(iblank_smoke3d[n21]==SOLID)value[3]=0;\
  }\
  if(firecolor==NULL&&value[0]==0&&value[1]==0&&value[2]==0&&value[3]==0)continue;\
  if(ABS(value[0]-value[2])<ABS(value[1]-value[3])){     \
    xyzindex=xyzindex1;                                  \
  }                                                      \
  else{                                                  \
    xyzindex=xyzindex2;                                  \
  }                                                      \
  if(firecolor!=NULL){\
    fvalue[0]=firecolor[n11];\
    fvalue[1]=firecolor[n12];\
    fvalue[2]=firecolor[n22];\
    fvalue[3]=firecolor[n21];\
    if(value[0]==0&&value[1]==0&&value[2]==0&&value[3]==0&&\
       fvalue[0]==0&&fvalue[1]==0&&fvalue[2]==0&&fvalue[3]==0)continue;\
  }\
  else{\
    fvalue[0]=0.0;\
    fvalue[1]=0.0;\
    fvalue[2]=0.0;\
    fvalue[3]=0.0;\
  }\
  for(node=0;node<6;node++){                             \
    mm = xyzindex[node];                                 \
    glVertexAttrib1f(GPU_hrr,fvalue[mm]);\
    glVertexAttrib1f(GPU_smokealpha,value[mm]); \
    glVertex3f(XX,YY,ZZ);                                \
  }

// -------------------------- DRAWVERTEXGPUTERRAIN ----------------------------------

#define DRAWVERTEXGPUTERRAIN(XX,YY,ZZ) \
  z_offset[XXX]=znode_offset[m11];\
  z_offset[YYY]=znode_offset[m12];\
  z_offset[ZZZ]=znode_offset[m22];\
  z_offset[3]=znode_offset[m21];\
  value[0]=alphaf_in[n11];\
  value[1]=alphaf_in[n12];\
  value[2]=alphaf_in[n22];\
  value[3]=alphaf_in[n21];\
  if((adjustalphaflag==2||adjustalphaflag==3)&&iblank_smoke3d!=NULL){\
    if(iblank_smoke3d[n11]==SOLID)value[0]=0;\
    if(iblank_smoke3d[n12]==SOLID)value[1]=0;\
    if(iblank_smoke3d[n22]==SOLID)value[2]=0;\
    if(iblank_smoke3d[n21]==SOLID)value[3]=0;\
  }\
  if(value[0]==0&&value[1]==0&&value[2]==0&&value[3]==0)continue;\
  if(ABS(value[0]-value[2])<ABS(value[1]-value[3])){     \
    xyzindex=xyzindex1;                                  \
  }                                                      \
  else{                                                  \
    xyzindex=xyzindex2;                                  \
  }                                                        \
  if(firecolor!=NULL){\
    fvalue[0]=firecolor[n11];\
    fvalue[1]=firecolor[n12];\
    fvalue[2]=firecolor[n22];\
    fvalue[3]=firecolor[n21];\
  }\
  else{\
    fvalue[0]=0.0;\
    fvalue[1]=0.0;\
    fvalue[2]=0.0;\
    fvalue[3]=0.0;\
  }\
  for(node=0;node<6;node++){                             \
    mm = xyzindex[node];                                 \
    glVertexAttrib1f(GPU_smokealpha,(float)value[mm]);\
    glVertexAttrib1f(GPU_hrr,(float)fvalue[mm]);\
    glVertex3f(XX,YY,ZZ+z_offset[mm]);                                \
  }

/* ------------------ is_fire_or_soot ------------------------ */

int is_fire_or_soot(smoke3ddata *smoke3di){
  smoke3ddata *soot=NULL, *fire=NULL;
  int is_fire, is_soot;
  int return_val;

  is_soot=0;
  if(smoke3di->soot_index!=-1)soot = smoke3dinfo + smoke3di->soot_index;
  if(soot!=NULL&&soot->loaded==1&&soot->display==1){
    is_soot=1;
    if(soot->frame_all_zeros[soot->ismoke3d_time]==1)is_soot=0;
  }

  is_fire=0;
  if(smoke3di->hrrpuv_index!=-1)fire = smoke3dinfo + smoke3di->hrrpuv_index;
  if(fire!=NULL&&fire->loaded==1&&fire->display==1){
    is_fire=1;
    if(fire->frame_all_zeros[fire->ismoke3d_time]==1)is_fire=0;
  }
  return_val=0;
  if(is_soot==1||is_fire==1)return_val=1;
  return return_val;
}

/* ------------------ draw_smokeframe ------------------------ */

void draw_smokeframe(void){
  CheckMemory;
  if(showvolrender==1&&smoke3dVoldebug==1){
    drawsmoke3dVOLdebug();
  }
#ifdef pp_GPU
  if(usegpu==1){
    if(showvolrender==1){
      LoadVolSmokeShaders();
    }
    else{
      LoadSmokeShaders();
    }
  }
#endif
#ifdef pp_CULL
  if(usegpu==1&&cullsmoke==1&&showvolrender==0){
    drawsmoke3dCULL();
  }
  else{
#else
  {
#endif
    int i;

    if(showvolrender==0){
      for(i=0;i<nsmoke3dinfo;i++){
        smoke3ddata *smoke3di;

        smoke3di = smoke3dinfo_sorted[i];
        if(smoke3di->loaded==0||smoke3di->display==0)continue;
        if(smoke3di->d_display==0)continue;
        if(is_fire_or_soot(smoke3di)==0)continue;

#ifdef pp_GPU
        if(usegpu==1){
          drawsmoke3dGPU(smoke3di);
        }
        else{
          drawsmoke3d(smoke3di);
        }
#else
        drawsmoke3d(smoke3di);
#endif
      }
    }
    if(showvolrender==1){
#ifdef pp_GPU
      if(usegpu==1){
      //  glBlendFunc(GL_ONE,GL_ONE_MINUS_SRC_ALPHA);
        SNIFF_ERRORS("before drawsmoke3dGPUVOL");
        drawsmoke3dGPUVOL();
        SNIFF_ERRORS("after drawsmoke3dGPUVOL");
      //  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
      }
      else{
        drawsmoke3dVOL();
        SNIFF_ERRORS("after drawsmoke3dVOL");
      }
#else
      drawsmoke3dVOL();
#endif
    }
  }
#ifdef pp_GPU
  if(usegpu==1){
    UnLoadShaders();
  }
#endif
#ifdef pp_CULL
  if(cullsmoke==1&&stereotype==STEREO_NONE&&nsmoke3dinfo>0&&show3dsmoke==1){
    setPixelCount();
  }
#endif
  SNIFF_ERRORS("after drawsmoke");
}

/* ------------------ readsmoke3d ------------------------ */

void readsmoke3d(int ifile,int flag, int *errorcode){
  smoke3ddata *smoke3di;
  FILE *SMOKE3DFILE;
  int error;
  int ncomp_smoke_total;
  int ncomp_smoke_total_skipped;
  int local_starttime=0, local_stoptime=0;
  FILE_SIZE file_size=0;
  int local_starttime0=0, local_stoptime0=0;
  float delta_time, delta_time0;
  int iii,ii,i,j;
  int nxyz[8];
  int nchars[2];
  int nframes_found=0;
  int have_light;

  float time_local;
  char compstring[128];
  meshdata *meshi;
  int fortran_skip;

  local_starttime0 = glutGet(GLUT_ELAPSED_TIME);
  ASSERT(ifile>=0&&ifile<nsmoke3dinfo);
  smoke3di = smoke3dinfo + ifile;
  if(smoke3di->filetype==FORTRAN_GENERATED){
    fortran_skip=4;
  }
  else{
    fortran_skip=0;
  }


  if(smoke3di->loaded==1){
    freesmoke3d(smoke3di);
    smoke3di->loaded=0;
    smoke3di->display=0;
    smoke3di->d_display=0;
  }

  meshi=meshinfo+smoke3di->blocknumber;
  FREEMEMORY(meshi->merge_alpha);
  FREEMEMORY(meshi->merge_color);

  if(flag==UNLOAD){
    plotstate=GetPlotState(DYNAMIC_PLOTS);
    UpdateTimes();
    Read3DSmoke3DFile=0;
    setsmokecolorflags();

    hrrpuv_loaded=0;
    for(j=0;j<nsmoke3dinfo;j++){
      smoke3ddata *smoke3dj;

      smoke3dj = smoke3dinfo + j;
      if(smoke3dj->loaded==1){
        Read3DSmoke3DFile=1;
        if(strcmp(smoke3dj->label.longlabel,"HRRPUV")==0){
          hrrpuv_loaded=1;
        }
        break;
      }
    }
    if(meshi->iblank_smoke3d != NULL){
      int free_iblank_smoke3d;

      free_iblank_smoke3d = 1;
      for (j = 0; j < nsmoke3dinfo; j++){
        smoke3ddata *smoke3dj;
        meshdata *meshj;

        smoke3dj = smoke3dinfo + j;
        meshj = meshinfo + smoke3dj->blocknumber;
        if(smoke3dj != smoke3di && smoke3dj->loaded == 1 && meshj == meshi){
          free_iblank_smoke3d = 0;
          break;
        }
      }
      if(free_iblank_smoke3d == 1){
        FREEMEMORY(meshi->iblank_smoke3d);
        meshi->iblank_smoke3d_defined = 0;
        update_makeiblank_smoke3d=1;
      }
    }
    update_makeiblank_smoke3d=1;
#ifdef pp_CULL
    if(cullactive==1)update_initcull=1;
#endif
    return;
  }

  if(smoke3di->compression_type==UNKNOWN){
    smoke3di->compression_type=getsmoke3d_version(smoke3di);
    update_smoke3d_menulabels();
  }

  CheckMemory;
  if(getsmoke3d_sizes(skip_global,smoke3di->file,smoke3di->compression_type,&smoke3di->times, &smoke3di->use_smokeframe,
                 &smoke3di->nchars_uncompressed,
                 &smoke3di->nchars_compressed_smoke,
                 &smoke3di->nchars_compressed_smoke_full,
                 &smoke3di->ntimes,
                 &smoke3di->ntimes_full,&have_light)==1){
    readsmoke3d(ifile,UNLOAD,&error);
    *errorcode=1;
    fprintf(stderr,"*** Error: problems sizing 3d smoke data for %s\n",smoke3di->file);
    return;
  }
  smoke3di->have_light=have_light;
  CheckMemory;
  if(
     NewMemory((void **)&smoke3di->smokeframe_comp_list,smoke3di->ntimes_full*sizeof(unsigned char *))==0||
     NewMemory((void **)&smoke3di->frame_all_zeros,smoke3di->ntimes_full*sizeof(unsigned char))==0||
     NewMemory((void **)&smoke3di->smokeframe_in,smoke3di->nchars_uncompressed*sizeof(unsigned char))==0||
     NewMemory((void **)&smoke3di->smokeview_tmp,smoke3di->nchars_uncompressed*sizeof(unsigned char))==0||
     NewMemory((void **)&smoke3di->smokeframe_out,smoke3di->nchars_uncompressed*sizeof(unsigned char))==0||
     NewMemory((void **)&meshi->merge_color,4*smoke3di->nchars_uncompressed*sizeof(unsigned char))==0||
     NewMemory((void **)&meshi->merge_alpha,smoke3di->nchars_uncompressed*sizeof(unsigned char))==0){
     readsmoke3d(ifile,UNLOAD,&error);
     *errorcode=1;
     fprintf(stderr,"*** Error: problems allocating memory for 3d smoke file: %s\n",smoke3di->file);
     return;
  }
  if(have_light==1){
    smoke3di->lightframe_in=smoke3di->smokeframe_in+smoke3di->nchars_uncompressed/2;
  }
  else{
    smoke3di->lightframe_in=NULL;
  }
  for(i=0;i<smoke3di->ntimes_full;i++){
    smoke3di->frame_all_zeros[i]=2;
  }

  ncomp_smoke_total=0;
  ncomp_smoke_total_skipped=0;
  for(i=0;i<smoke3di->ntimes_full;i++){
    ncomp_smoke_total+=smoke3di->nchars_compressed_smoke_full[i];
    if(smoke3di->use_smokeframe[i]==1){
      ncomp_smoke_total_skipped+=smoke3di->nchars_compressed_smoke_full[i];
    }
  }
  if(NewMemory((void **)&smoke3di->smoke_comp_all,ncomp_smoke_total_skipped*sizeof(unsigned char))==0){
    readsmoke3d(ifile,UNLOAD,&error);
    *errorcode=1;
     fprintf(stderr,"*** Error: problems allocating memory for 3d smoke file: %s\n",smoke3di->file);
    return;
  }
  smoke3di->ncomp_smoke_total=ncomp_smoke_total_skipped;

  ncomp_smoke_total=0;
  i=0;
  for(ii=0;ii<smoke3di->ntimes_full;ii++){
    if(smoke3di->use_smokeframe[ii]==1){
      smoke3di->smokeframe_comp_list[i]=smoke3di->smoke_comp_all+ncomp_smoke_total;
      ncomp_smoke_total+=smoke3di->nchars_compressed_smoke[i];
      i++;
    }
  }

// read in data

  file_size=get_filesize(smoke3di->file);
  SMOKE3DFILE=fopen(smoke3di->file,"rb");
  if(SMOKE3DFILE==NULL){
    readsmoke3d(ifile,UNLOAD,&error);
    *errorcode=1;
    return;
  }

  SKIP;fread(nxyz,4,8,SMOKE3DFILE);SKIP;
  smoke3di->is1=nxyz[2];
  smoke3di->is2=nxyz[3];
  smoke3di->js1=nxyz[4];
  smoke3di->js2=nxyz[5];
  smoke3di->ks1=nxyz[6];
  smoke3di->ks2=nxyz[7];
  smoke3di->compression_type=nxyz[1];

  // read smoke data

  iii=0;
  nframes_found=0;
  local_starttime = glutGet(GLUT_ELAPSED_TIME);
  for(i=0;i<smoke3di->ntimes_full;i++){
	  SKIP;fread(&time_local,4,1,SMOKE3DFILE);SKIP;
    if(feof(SMOKE3DFILE)!=0||(use_tload_end==1&&time_local>tload_end)){
      smoke3di->ntimes_full=i;
      smoke3di->ntimes=nframes_found;
      break;
    }
    if(use_tload_begin==1&&time_local<tload_begin)smoke3di->use_smokeframe[i]=0;
    if(smoke3di->use_smokeframe[i]==1){
      PRINTF("3D smoke/fire time=%.2f",time_local);
    }
    SKIP;fread(nchars,4,2,SMOKE3DFILE);SKIP;
    if(feof(SMOKE3DFILE)!=0){
      smoke3di->ntimes_full=i;
      smoke3di->ntimes=nframes_found;
      break;
    }
    if(smoke3di->use_smokeframe[i]==1){
      float complevel;

      nframes_found++;
      SKIP;fread(smoke3di->smokeframe_comp_list[iii],1,smoke3di->nchars_compressed_smoke[iii],SMOKE3DFILE);SKIP;
      iii++;
      CheckMemory;
      if(feof(SMOKE3DFILE)!=0){
        smoke3di->ntimes_full=i;
        smoke3di->ntimes=nframes_found;
        break;
      }

      complevel=40.0*(float)nchars[0]/(float)nchars[1];
      complevel=(int)complevel;
      complevel/=10.0;
      if(complevel<0.0)complevel=-complevel;
      sprintf(compstring," compression ratio: %.1f",complevel);
      trim_back(compstring);
      trimzeros(compstring);
      PRINTF("%s\n",compstring);
    }
    else{
      SKIP;FSEEK(SMOKE3DFILE,smoke3di->nchars_compressed_smoke_full[i],SEEK_CUR);SKIP;
      if(feof(SMOKE3DFILE)!=0){
        smoke3di->ntimes_full=i;
        smoke3di->ntimes=nframes_found;
        break;
      }
    }
  }
  local_stoptime = glutGet(GLUT_ELAPSED_TIME);
  delta_time = (local_stoptime-local_starttime)/1000.0;

  if(SMOKE3DFILE!=NULL){
    fclose(SMOKE3DFILE);
  }

  smoke3di->loaded=1;
  smoke3di->display=1;
  setsmokecolorflags();
  if(strcmp(smoke3di->label.longlabel,"HRRPUV")==0){
    hrrpuv_loaded=1;
  }

  Read3DSmoke3DFile=1;
  update_makeiblank_smoke3d=1;
  plotstate=GetPlotState(DYNAMIC_PLOTS);
  UpdateTimes();
#ifdef pp_CULL
    if(cullactive==1)initcull(cullsmoke);
#endif
  Smoke3d_CB(UPDATE_SMOKEFIRE_COLORS);
  Idle_CB();
  local_stoptime0 = glutGet(GLUT_ELAPSED_TIME);
  delta_time0=(local_stoptime0-local_starttime0)/1000.0;

  if(file_size!=0&&delta_time>0.0){
    float loadrate;

    loadrate = ((float)file_size*8.0/1000000.0)/delta_time;
    PRINTF(" %.1f MB loaded in %.2f s - rate: %.1f Mb/s (overhead: %.2f s)\n",
    (float)file_size/1000000.,delta_time,loadrate,delta_time0-delta_time);
  }
  else{
    PRINTF(" %.1f MB downloaded in %.2f s (overhead: %.2f s)",
    (float)file_size/1000000.,delta_time,delta_time0-delta_time);
  }
#ifdef pp_MEMPRINT
  PRINTF("After 3D Smoke load: \n");
  PrintMemoryInfo;
#endif
  *errorcode = 0;
}

/* ------------------ setsmokecolorflags ------------------------ */

void setsmokecolorflags(void){
  int i;

  for(i=0;i<nsmoke3dinfo;i++){
    smoke3ddata *smoke3di;

    smoke3di = smoke3dinfo + i;
    smoke3di->soot_loaded=0;
    smoke3di->water_loaded=0;
    smoke3di->hrrpuv_loaded=0;
  }

  for(i=0;i<nsmoke3dinfo;i++){
    smoke3ddata *smoke3di;
    int j;

    smoke3di=smoke3dinfo + i;
    smoke3di->water_color=NULL;
    smoke3di->hrrpuv_color=NULL;
    smoke3di->soot_color=NULL;
    smoke3di->water_index=-1;
    smoke3di->hrrpuv_index=-1;
    smoke3di->soot_index=-1;
    if(smoke3di->loaded==0)continue;

    switch(smoke3di->type){
    case SOOT:
      smoke3di->soot_color=smoke3di->smokeframe_in;
      smoke3di->soot_index=i;
      break;
    case FIRE:
      smoke3di->hrrpuv_color=smoke3di->smokeframe_in;
      smoke3di->hrrpuv_index=i;
      break;
    case WATER:
      smoke3di->water_color=smoke3di->smokeframe_in;
      smoke3di->water_index=i;
      break;
    default:
      ASSERT(FFALSE);
      break;
    }

    for(j=0;j<nsmoke3dinfo;j++){
      smoke3ddata *smoke3dj;

      if(i==j)continue;
      smoke3dj = smoke3dinfo + j;
      if(smoke3dj->loaded==0)continue;
      if(smoke3di->blocknumber!=smoke3dj->blocknumber)continue;
      if(smoke3di->is1!=smoke3dj->is1)continue;
      if(smoke3di->is2!=smoke3dj->is2)continue;
      if(smoke3di->js1!=smoke3dj->js1)continue;
      if(smoke3di->js2!=smoke3dj->js2)continue;
      if(smoke3di->ks1!=smoke3dj->ks1)continue;
      if(smoke3di->ks2!=smoke3dj->ks2)continue;

      switch(smoke3dj->type){
      case SOOT:
        smoke3di->soot_color=smoke3dj->smokeframe_in;
        smoke3di->soot_index=j;
        break;
      case FIRE:
        smoke3di->hrrpuv_color=smoke3dj->smokeframe_in;
        smoke3di->hrrpuv_index=j;
        break;
      case WATER:
        smoke3di->water_color=smoke3dj->smokeframe_in;
        smoke3di->water_index=j;
        break;
      default:
        ASSERT(FFALSE);
        break;
      }
      if(smoke3di->soot_color!=NULL)smoke3dj->soot_loaded=1;
      if(smoke3di->water_color!=NULL)smoke3dj->water_loaded=1;
      if(smoke3di->hrrpuv_color!=NULL)smoke3dj->hrrpuv_loaded=1;
    }
  }
  for(i=0;i<nsmoke3dinfo;i++){
    int j;
    smoke3ddata *smoke3di;

    smoke3di = smoke3dinfo + i;
    if(smoke3di->loaded==0||smoke3di->display==0||smoke3di->frame_all_zeros==NULL)continue;
    for(j=0;j<smoke3di->ntimes_full;j++){
      smoke3di->frame_all_zeros[j]=2;
    }
  }
}

/* ------------------ getsmoke3d_sizes ------------------------ */

int getsmoke3d_sizes(int fortran_skip,char *smokefile, int version, float **timelist_found, int **use_smokeframe,
                      int *nchars_smoke_uncompressed,
                      int **nchars_smoke_compressed_found,
                      int **nchars_smoke_compressed_full,
                      int *ntimes_found,int *ntimes_full, int *have_light){
  char smoke_sizefilename[1024],smoke_sizefilename2[1024],buffer[255];
  char *textptr;
  FILE *SMOKE_SIZE=NULL;
  FILE *SMOKE3DFILE;
  int nframes_found;
  float time_local,time_max,*time_found=NULL;
  int *use_smokeframe_full;
  int nch_uncompressed,nch_smoke_compressed;
  int *nch_smoke_compressed_full=NULL;
  int *nch_smoke_compressed_found=NULL;
  int nxyz[8];
  int nchars[2];
  int skip_local;
  int iframe_local;
  int dummy;
  int ntimes_full2;
  size_t lentext;
  int iii;

  *have_light=0;
  // try .sz
  strcpy(smoke_sizefilename,smokefile);
  lentext=strlen(smoke_sizefilename);
  if(lentext>4){
    textptr=smoke_sizefilename + lentext - 4;
    if(strcmp(textptr,".svz")==0){
      smoke_sizefilename[lentext-4]=0;
    }
  }
  strcat(smoke_sizefilename,".szz");
  SMOKE_SIZE=fopen(smoke_sizefilename,"r");

  if(SMOKE_SIZE==NULL){
    // not .szz so try .sz
    strcpy(smoke_sizefilename,smokefile);
    strcat(smoke_sizefilename,".sz");
    SMOKE_SIZE=fopen(smoke_sizefilename,"r");
  }

  if(SMOKE_SIZE==NULL){
    // neither .sz or .szz so try in tempdir
    strcpy(smoke_sizefilename,smokefile);
    strcat(smoke_sizefilename,".sz");
    SMOKE_SIZE=fopen(smoke_sizefilename,"w");
    if(SMOKE_SIZE==NULL&&smokeviewtempdir!=NULL){
      strcpy(smoke_sizefilename2,smokeviewtempdir);
      strcat(smoke_sizefilename2,smoke_sizefilename);
      strcpy(smoke_sizefilename,smoke_sizefilename2);
      SMOKE_SIZE=fopen(smoke_sizefilename,"w");
    }
    if(SMOKE_SIZE==NULL)return 1;  // can't write size file in temp directory so give up
    SMOKE3DFILE=fopen(smokefile,"rb");
    if(SMOKE3DFILE==NULL){
      fclose(SMOKE_SIZE);
      return 1;
    }

    SKIP;fread(nxyz,4,8,SMOKE3DFILE);SKIP;

    if(version!=1)version=0;
    fprintf(SMOKE_SIZE,"%i\n",version);

    time_max=-1000000.0;
    for(;;){
      int nframeboth;

      SKIP;fread(&time_local,4,1,SMOKE3DFILE);SKIP;
      if(feof(SMOKE3DFILE)!=0)break;
      SKIP;fread(nchars,4,2,SMOKE3DFILE);SKIP;
      if(version==0){ // uncompressed data
        fprintf(SMOKE_SIZE,"%f %i %i\n",time_local,nchars[0],nchars[1]);
      }
      else{  // compressed data
        int nlightdata;

        // time, nframeboth, ncompressed_rle, ncompressed_zlib, nlightdata
        // ncompessed_zlib and nlightdata are negative if there is radiance data present

        if(nchars[1]<0){  // light data present
          nframeboth = nchars[0];
          nlightdata = -nchars[0]/2;
          fprintf(SMOKE_SIZE,"%f %i %i %i %i \n",time_local,nframeboth,-1,nchars[1],nlightdata);
        }
        else{
          nframeboth = nchars[0];
          nlightdata = 0;
          fprintf(SMOKE_SIZE,"%f %i %i %i %i \n",time_local,nframeboth,-1,nchars[1],nlightdata);
        }
      }
      skip_local = ABS(nchars[1]);
      SKIP;FSEEK(SMOKE3DFILE,skip_local,SEEK_CUR);SKIP;
    }

    fclose(SMOKE3DFILE);
    fclose(SMOKE_SIZE);
    SMOKE_SIZE=fopen(smoke_sizefilename,"r");
  }
  if(SMOKE_SIZE==NULL)return 1;

  nframes_found=0;
  iframe_local=-1;
  time_max=-1000000.0;
  fgets(buffer,255,SMOKE_SIZE);
  iii=0;
  while(!feof(SMOKE_SIZE)){
    if(fgets(buffer,255,SMOKE_SIZE)==NULL)break;
    sscanf(buffer,"%f",&time_local);
    iframe_local++;
    if(time_local<=time_max)continue;
    if(use_tload_end==1&&time_local>tload_end)break;
    if(iii%smoke3dframestep==0&&(use_tload_begin==0||time_local>=tload_begin)){
      nframes_found++;
      time_max=time_local;
    }
    iii++;
  }
  rewind(SMOKE_SIZE);
  if(nframes_found<=0){
    *ntimes_found=0;
    *ntimes_full=0;
    fclose(SMOKE_SIZE);
    return 1;
  }
  *ntimes_found=nframes_found;
  *ntimes_full=iframe_local+1;

  NewMemory((void **)&use_smokeframe_full,*ntimes_full*sizeof(float));
  NewMemory((void **)&time_found,nframes_found*sizeof(float));
  NewMemory((void **)&nch_smoke_compressed_full,(*ntimes_full)*sizeof(int));
  NewMemory((void **)&nch_smoke_compressed_found,(*ntimes_found)*sizeof(int));

  *use_smokeframe=use_smokeframe_full;
  *timelist_found=time_found;
  *nchars_smoke_compressed_full=nch_smoke_compressed_full;
  *nchars_smoke_compressed_found=nch_smoke_compressed_found;

  fgets(buffer,255,SMOKE_SIZE);
  ntimes_full2=0;
  time_max=-1000000.0;
  iii=0;
  while(!feof(SMOKE_SIZE)){
    int nch_light;

    if(fgets(buffer,255,SMOKE_SIZE)==NULL)break;
    ntimes_full2++;
    if(ntimes_full2>*ntimes_full)break;
    if(version==0){
      sscanf(buffer,"%f %i %i",&time_local,&nch_uncompressed,&nch_smoke_compressed);
    }
    else{
      nch_light=0;
      sscanf(buffer,"%f %i %i %i %i",&time_local,&nch_uncompressed,&dummy,&nch_smoke_compressed,&nch_light);
      if(nch_light<0)*have_light=1;
    }
    *nch_smoke_compressed_full++=nch_smoke_compressed;
    *use_smokeframe_full=0;
    if(use_tload_end==1&&time_local>tload_end)break;
    if(time_local<=time_max){
      use_smokeframe_full++;
      continue;
    }

    if(iii%smoke3dframestep==0&&(use_tload_begin==0||time_local>=tload_begin)){
      *use_smokeframe_full=1;
      *time_found++=time_local;
      time_max=time_local;
      *nch_smoke_compressed_found++=nch_smoke_compressed;
    }
    use_smokeframe_full++;
    iii++;
  }
  *nchars_smoke_uncompressed=nch_uncompressed;
  fclose(SMOKE_SIZE);
  return 0;
}

/* ------------------ freesmoke3d ------------------------ */

void freesmoke3d(smoke3ddata *smoke3di){

  smoke3di->lastiframe=-999;
  FREEMEMORY(smoke3di->smokeframe_in);
  FREEMEMORY(smoke3di->smokeframe_out);
  FREEMEMORY(smoke3di->timeslist);
  FREEMEMORY(smoke3di->times);
  FREEMEMORY(smoke3di->use_smokeframe);
  FREEMEMORY(smoke3di->nchars_compressed_smoke_full);
  FREEMEMORY(smoke3di->nchars_compressed_smoke);
  FREEMEMORY(smoke3di->frame_all_zeros);
  FREEMEMORY(smoke3di->smoke_comp_all);
  FREEMEMORY(smoke3di->smokeframe_comp_list);
  FREEMEMORY(smoke3di->smokeview_tmp);
}

/* ------------------ updatesmoke3d ------------------------ */

void updatesmoke3d(smoke3ddata *smoke3di){
  int iframe_local;
  int countin;
  uLongf countout;

  iframe_local = smoke3di->ismoke3d_time;
  countin = smoke3di->nchars_compressed_smoke[iframe_local];
  countout=smoke3di->nchars_uncompressed;
  switch(smoke3di->compression_type){
  case RLE:
    countout = uncompress_rle(smoke3di->smokeframe_comp_list[iframe_local],countin,smoke3di->smokeframe_in);
    break;
  case ZLIB:
    uncompress_zlib(
      smoke3di->smokeframe_in,&countout,
      smoke3di->smokeframe_comp_list[iframe_local],countin);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  if(smoke3di->frame_all_zeros[iframe_local]==2){
    int i;
    unsigned char *smokeframe_in;

    smokeframe_in = smoke3di->smokeframe_in;
    smoke3di->frame_all_zeros[iframe_local]=1;
    for(i=0;i<smoke3di->nchars_uncompressed;i++){
      if(smokeframe_in[i]!=0){
        smoke3di->frame_all_zeros[iframe_local]=0;
        break;
      }
    }
  }
  ASSERT(countout==smoke3di->nchars_uncompressed);
}

/* ------------------ mergesmokecolors ------------------------ */

void mergesmoke3dcolors(smoke3ddata *smoke3dset){
  int i,j;
  int i_hrrpuv_cutoff;

  i_hrrpuv_cutoff=254*global_hrrpuv_cutoff/hrrpuv_max_smv;

#ifdef pp_CULL
  for(i=0;i<nmeshes;i++){
    meshdata *meshi;

    meshi = meshinfo + i;
    meshi->cull_smoke3d=NULL;
  }
#endif

  for(i=0;i<nsmoke3dinfo;i++){
    smoke3ddata *smoke3di;
    smoke3ddata *smoke_soot;

    smoke3di=smoke3dinfo + i;
    if(smoke3dset!=NULL&&smoke3dset!=smoke3di)continue;
    smoke3di->d_display=0;
    if(smoke3di->loaded==0||smoke3di->display==0)continue;
    smoke_soot=NULL;
    switch(smoke3di->type){
    case SOOT:
      smoke3di->d_display=1;
      break;
    case FIRE:
      if(smoke3di->soot_index!=-1)smoke_soot=smoke3dinfo + smoke3di->soot_index;
      if(smoke3di->soot_loaded==0||(smoke_soot!=NULL&&smoke_soot->display==0)){
        smoke3di->d_display=1;
      }
      break;
    case WATER:
      if(smoke3di->soot_loaded==0&&smoke3di->hrrpuv_loaded==0){
        smoke3di->d_display=1;
      }
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
  }

  for(i=0;i<nsmoke3dinfo;i++){
    smoke3ddata *smoke3di;
    meshdata *meshi;
    float fire_alpha;
    unsigned char *firecolor,*sootcolor;
    unsigned char *mergecolor,*mergealpha;

    smoke3di = smoke3dinfo + i;
    if(smoke3dset!=NULL&&smoke3dset!=smoke3di)continue;
    if(smoke3di->loaded==0||smoke3di->d_display==0)continue;
    meshi=meshinfo+smoke3di->blocknumber;
#ifdef pp_CULL
    meshi->cull_smoke3d=smoke3di;
#endif
    i_hrrpuv_cutoff=254*global_hrrpuv_cutoff/hrrpuv_max_smv;

    if(fire_halfdepth<=0.0){
      smoke3di->fire_alpha=255;
    }
    else{
      smoke3di->fire_alpha=255*(1.0-pow(0.5,meshi->dx/fire_halfdepth));
    }
    fire_alpha=smoke3di->fire_alpha;

    firecolor=NULL;
    sootcolor=NULL;
    if(smoke3di->hrrpuv_color!=NULL){
      firecolor=smoke3di->hrrpuv_color;
      if(smoke3di->hrrpuv_index!=-1){
        smoke3ddata *smoke3dref;

        smoke3dref=smoke3dinfo + smoke3di->hrrpuv_index;
        if(smoke3dref->display==0)firecolor=NULL;
      }
    }
    if(smoke3di->soot_color!=NULL){
      sootcolor=smoke3di->soot_color;
      if(smoke3di->soot_index!=-1){
        smoke3ddata *smoke3dref;

        smoke3dref=smoke3dinfo + smoke3di->soot_index;
        if(smoke3dref->display==0){
          sootcolor=NULL;
        }
      }
    }
#ifdef pp_GPU
    if(usegpu==1)continue;
#endif
    if(meshi->merge_color==NULL){
      NewMemory((void **)&meshi->merge_color,4*smoke3di->nchars_uncompressed*sizeof(unsigned char));
    }
    if(meshi->merge_alpha==NULL){
      NewMemory((void **)&meshi->merge_alpha,smoke3di->nchars_uncompressed*sizeof(unsigned char));
    }

    mergecolor=meshi->merge_color;
    mergealpha=meshi->merge_alpha;
    if(sootcolor!=NULL){
      for(j=0;j<smoke3di->nchars_uncompressed;j++){
        float *firesmoke;

        if(firecolor!=NULL&&firecolor[j]>i_hrrpuv_cutoff){
          firesmoke=rgb_slicesmokecolormap+4*firecolor[j];
          *mergecolor++ = 255*firesmoke[0];
          *mergecolor++ = 255*firesmoke[1];
          *mergecolor++ = 255*firesmoke[2];
          *mergealpha++=fire_alpha;
        }
        else{
          firesmoke=rgb_slicesmokecolormap+4*sootcolor[j];
          *mergecolor++ = 255*smoke_albedo;
          *mergecolor++ = 255*smoke_albedo;
          *mergecolor++ = 255*smoke_albedo;
          *mergealpha++=(sootcolor[j]>>smoke3d_thick);
        }
        mergecolor++;
      }
      continue;
    }
    else{
      if(firecolor!=NULL){
        for(j=0;j<smoke3di->nchars_uncompressed;j++){
          float *fire;

          fire=rgb_slicesmokecolormap+4*firecolor[j];
          *mergecolor++=255*fire[0];
          *mergecolor++=255*fire[1];
          *mergecolor++=255*fire[2];
          mergecolor++;
          if(firecolor[j]>i_hrrpuv_cutoff){
            *mergealpha++=fire_alpha;
          }
          else{
            *mergealpha++=0;
          }
        }
        continue;
      }
      else{
        ASSERT(FFALSE);
        // (fire==NULL&&soot==NULL) this should never happen!!!
      }
    }
  }
}

/* ------------------ drawsmoke3d ------------------------ */

void drawsmoke3d(smoke3ddata *smoke3di){
  int i,j,k,n;
  float constval,x1,x3,z1,z3, yy1, y3;
  int is1, is2, js1, js2, ks1, ks2;
  int ii, jj, kk;
  int ibeg, iend, jbeg, jend, kbeg, kend;
  float *znode_offset, z_offset[4];
  int have_light;

  float *xplt, *yplt, *zplt;
  unsigned char mergealpha,*mergealphaptr,*mergecolorptr;
  int nx,ny,nz;
  unsigned char *alphaf_in,*alphaf_out,*alphaf_ptr,*light_in;
  unsigned char *colorptr;
  int xyzindex1[6],xyzindex2[6],*xyzindex,node;
  float xnode[4],znode[4],ynode[4];
  int skip_local;
    int iterm, jterm, kterm,nxy;
  float x11[3], x12[3], x22[3], x21[3];
  unsigned int n11, n12, n22, n21;
  int ipj,jpk,ipk,jmi,kmi,kmj;
  int iii, jjj, kkk;
  int slice_end,slice_beg;
  float aspectratio;
  int ssmokedir;
  unsigned char *iblank_smoke3d;

  unsigned char value[4],light_value[4];
  int ivalue[4];

  meshdata *meshi;

  have_light = smoke3di->have_light;
  meshi = meshinfo + smoke3di->blocknumber;
  if(meshvisptr[meshi-meshinfo]==0)return;

  if(meshi->merge_alpha==NULL||meshi->update_firehalfdepth==1){
    meshi->update_firehalfdepth=0;
    mergesmoke3dcolors(smoke3di);
  }
  mergealphaptr = meshi->merge_alpha;
  mergecolorptr = meshi->merge_color;
  value[0]=255;
  value[1]=255;
  value[2]=255;
  value[3]=255;

  if(nterraininfo>0){
    znode_offset = meshi->terrain->znode_offset;
  }

  xplt=meshi->xplt;
  yplt=meshi->yplt;
  zplt=meshi->zplt;
  iblank_smoke3d = meshi->iblank_smoke3d;
  alphaf_in=smoke3di->smokeframe_in;
  light_in=smoke3di->lightframe_in;
  alphaf_out=smoke3di->smokeframe_out;

  switch(demo_mode){
  case 0:
    is1 = smoke3di->is1;
    is2 = smoke3di->is2;
    break;
  case 1:
    is1 = (smoke3di->is1+smoke3di->is2)/2;
    is2 = is1+1;
    if(is1<smoke3di->is1)is1=smoke3di->is1;
    if(is2>smoke3di->is2)is2=smoke3di->is2;
    break;
  default:
    is1 = (smoke3di->is1+smoke3di->is2)/2-demo_mode;
    is2 = is1+2*demo_mode;
    if(is1<smoke3di->is1)is1=smoke3di->is1;
    if(is2>smoke3di->is2)is2=smoke3di->is2;
    break;
  }
  js1 = smoke3di->js1;
  js2 = smoke3di->js2;
  ks1 = smoke3di->ks1;
  ks2 = smoke3di->ks2;

  nx = smoke3di->is2 + 1 - smoke3di->is1;
  ny = js2 + 1 - js1;
  nz = ks2 + 1 - ks1;
  nxy = nx*ny;

  ssmokedir=meshi->smokedir;
  skip_local=smokeskipm1+1;

  xyzindex1[0]=0;
  xyzindex1[1]=1;
  xyzindex1[2]=2;
  xyzindex1[3]=0;
  xyzindex1[4]=2;
  xyzindex1[5]=3;

  xyzindex2[0]=0;
  xyzindex2[1]=1;
  xyzindex2[2]=3;
  xyzindex2[3]=1;
  xyzindex2[4]=2;
  xyzindex2[5]=3;

  if(cullfaces==1)glDisable(GL_CULL_FACE);

  TransparentOn();
  switch(ssmokedir){

    // +++++++++++++++++++++++++++++++++++ DIR 1 +++++++++++++++++++++++++++++++++++++++

  case 1:
  case -1:

    // ++++++++++++++++++  adjust transparency +++++++++++++++++

    if(adjustalphaflag!=0){

      aspectratio=meshi->dx;
      for(i=is1;i<=is2;i++){
        iterm=(i-smoke3di->is1);

        if(smokecullflag==1){
          x11[0]=xplt[i];
          x12[0]=xplt[i];
          x22[0]=xplt[i];
          x21[0]=xplt[i];

          x11[1]=yplt[js1];
          x12[1]=yplt[js2];
          x22[1]=yplt[js2];
          x21[1]=yplt[js1];

          x11[2]=zplt[ks1];
          x12[2]=zplt[ks1];
          x22[2]=zplt[ks2];
          x21[2]=zplt[ks2];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(k=ks1;k<=ks2;k++){
          kterm=(k-ks1)*nxy;

          if(smokecullflag==1&&k!=ks2){
            x11[2]=zplt[k];
            x12[2]=zplt[k];
            x22[2]=zplt[k+1];
            x21[2]=zplt[k+1];
            if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
          }

          for(j=js1;j<=js2;j++){
            jterm = (j-js1)*nx;
          //  jterm = (j-js1)*nx;
            n = iterm + jterm + kterm;
            ASSERT(n>=0&&n<smoke3di->nchars_uncompressed);
            mergealpha = mergealphaptr[n];
            ADJUSTALPHA(mergealpha,aspectratio,NULL,1);
          }
        }
      }
      alphaf_ptr=alphaf_out;
    }
    else{
      alphaf_ptr=alphaf_in;
    }

    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);
    if(smokedrawtest==1||show_smoketest==1){
      slice_end=smokedrawtest_nummax;
      if(slice_end>is2)slice_end=is2;
      slice_beg=smokedrawtest_nummin;
      if(slice_beg<is1)slice_beg=is1;
      if(show_smoketest==1){
        slice_end=is2-1;
        slice_beg=slice_end-1;
      }
    }
    else{
      slice_beg=is1;
      slice_end=is2;
    }
    for(iii=slice_beg;iii<slice_end;iii+=skip_local){
      i=iii;
      if(ssmokedir<0)i = is1+is2-iii-1;
      iterm = (i-smoke3di->is1);

      if(smokecullflag==1){
        x11[0]=xplt[i];
        x12[0]=xplt[i];
        x22[0]=xplt[i];
        x21[0]=xplt[i];

        x11[1]=yplt[js1];
        x12[1]=yplt[js2];
        x22[1]=yplt[js2];
        x21[1]=yplt[js1];

        x11[2]=zplt[ks1];
        x12[2]=zplt[ks1];
        x22[2]=zplt[ks2];
        x21[2]=zplt[ks2];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      constval = xplt[i]+0.001;
      for(k=ks1; k<ks2; k++){
        kterm = (k-ks1)*nxy;
        z1 = zplt[k];
        z3 = zplt[k+1];
        znode[0]=z1;
        znode[1]=z1;
        znode[2]=z3;
        znode[3]=z3;

        if(smokecullflag==1&&k!=ks2){
          x11[2]=zplt[k];
          x12[2]=zplt[k];
          x22[2]=zplt[k+1];
          x21[2]=zplt[k+1];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }
        for(j=js1; j<js2; j++){
          jterm = (j-js1)*nx;
          yy1 = yplt[j];
          y3 = yplt[j+1];
          ynode[0]=yy1;
          ynode[1]=y3;
          ynode[2]=y3;
          ynode[3]=yy1;

          n = iterm + jterm + kterm;

          n11 = n;              //n
          n12 = n11 + nx;       //n+nx
          n22 = n12 + nxy;      //n+nx+nxy
          n21 = n22 - nx;       //n+nxy

//        n11 = (i-is1)   + (j-js1)*nx   + (k-ks1)*nx*ny;
//        n12 = (i-is1)   + (j+1-js1)*nx + (k-ks1)*nx*ny;
//        n22 = (i-is1)   + (j+1-js1)*nx + (k+1-ks1)*nx*ny;
//        n21 = (i-is1)   + (j-js1)*nx   + (k+1-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + nx;
            m22 = m12;
            m21 = m22 - nx;
            DRAWVERTEXTERRAIN(constval,ynode[mm],znode[mm])
          }
          else{
            DRAWVERTEX(constval,ynode[mm],znode[mm])
          }

        }
      }
    }
    glEnd();

    break;

  // +++++++++++++++++++++++++++++++++++ DIR 2 +++++++++++++++++++++++++++++++++++++++

  case 2:
  case -2:

    // ++++++++++++++++++  adjust transparency +++++++++++++++++

    if(adjustalphaflag!=0){

      aspectratio=meshi->dy;
      for(j=js1;j<=js2;j++){
        jterm = (j-js1)*nx;
    //    xp[1]=yplt[j];

        if(smokecullflag==1){
          x11[0]=xplt[is1];
          x12[0]=xplt[is2];
          x22[0]=xplt[is2];
          x21[0]=xplt[is1];

          x11[1]=yplt[j];
          x12[1]=yplt[j];
          x22[1]=yplt[j];
          x21[1]=yplt[j];

          x11[2]=zplt[ks1];
          x12[2]=zplt[ks1];
          x22[2]=zplt[ks2];
          x21[2]=zplt[ks2];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(k=ks1;k<=ks2;k++){
          kterm = (k-ks1)*nxy;

          if(smokecullflag==1&&k!=ks2){
            x11[2]=zplt[k];
            x12[2]=zplt[k];
            x22[2]=zplt[k+1];
            x21[2]=zplt[k+1];

            if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
          }

          for(i=is1;i<=is2;i++){
            iterm = (i-is1);
            n = iterm + jterm + kterm;
            ASSERT(n>=0&&n<smoke3di->nchars_uncompressed);
            mergealpha = mergealphaptr[n];
            ADJUSTALPHA(mergealpha,aspectratio,NULL,2);
          }
        }
      }
      alphaf_ptr=alphaf_out;
    }
    else{
      alphaf_ptr=alphaf_in;
    }

    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);
    if(smokedrawtest==1||show_smoketest==1){
      slice_end=smokedrawtest_nummax;
      if(slice_end>js2)slice_end=js2;
      slice_beg=smokedrawtest_nummin;
      if(slice_beg<js1)slice_beg=js1;
      if(show_smoketest==1){
        slice_end=js2-1;
        slice_beg=slice_end-1;
      }
    }
    else{
      slice_beg=js1;
      slice_end=js2;
    }
    for(jjj=slice_beg;jjj<slice_end;jjj+=skip_local){
      j=jjj;
      if(ssmokedir<0)j = js1+js2-jjj-1;
      constval = yplt[j]+0.001;
      jterm = (j-js1)*nx;

      if(smokecullflag==1){
        x11[0]=xplt[is1];
        x12[0]=xplt[is2];
        x22[0]=xplt[is2];
        x21[0]=xplt[is1];

        x11[1]=yplt[j];
        x12[1]=yplt[j];
        x22[1]=yplt[j];
        x21[1]=yplt[j];

        x11[2]=zplt[ks1];
        x12[2]=zplt[ks1];
        x22[2]=zplt[ks2];
        x21[2]=zplt[ks2];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(k=ks1; k<ks2; k++){
        kterm = (k-ks1)*nxy;
        z1 = zplt[k];
        z3 = zplt[k+1];

        znode[0]=z1;
        znode[1]=z1;
        znode[2]=z3;
        znode[3]=z3;

        if(smokecullflag==1){
          x11[2]=zplt[k];
          x12[2]=zplt[k];
          x22[2]=zplt[k+1];
          x21[2]=zplt[k+1];
          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(i=is1; i<is2; i++){
          iterm = (i-is1);
          x1 = xplt[i];
          x3 = xplt[i+1];

          xnode[0]=x1;
          xnode[1]=x3;
          xnode[2]=x3;
          xnode[3]=x1;

          n = iterm + jterm + kterm;
          n11 = n;            //n
          n12 = n11+1;;       //n+1
          n22 = n12+nxy;      //n+1+nxy
          n21 = n22-1;        //n+nxy

//        n11 = (i-is1)   + (j-js1)*nx   + (k-ks1)*nx*ny;
//        n12 = (i+1-is1) + (j-js1)*nx   + (k-ks1)*nx*ny;
//        n22 = (i+1-is1) + (j-js1)*nx   + (k+1-ks1)*nx*ny;
//        n21 = (i-is1)   + (j-js1)*nx   + (k+1-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + 1;
            m22 = m12;
            m21 = m22 - 1;
            DRAWVERTEXTERRAIN(xnode[mm],constval,znode[mm])
          }
          else{
            DRAWVERTEX(xnode[mm],constval,znode[mm])
          }
        }
      }
    }
    glEnd();
    break;

  // +++++++++++++++++++++++++++++++++++ DIR 3 +++++++++++++++++++++++++++++++++++++++

  case 3:
  case -3:
    // ++++++++++++++++++  adjust transparency +++++++++++++++++

    aspectratio=meshi->dz;
    if(adjustalphaflag!=0){
      for(k=ks1;k<=ks2;k++){
        kterm = (k-ks1)*nxy;

        if(smokecullflag==1){
          x11[0]=xplt[is1];
          x12[0]=xplt[is2];
          x22[0]=xplt[is2];
          x21[0]=xplt[is1];

          x11[1]=yplt[js1];
          x12[1]=yplt[js1];
          x22[1]=yplt[js2];
          x21[1]=yplt[js2];

          x11[2]=zplt[k];
          x12[2]=zplt[k];
          x22[2]=zplt[k];
          x21[2]=zplt[k];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(j=js1;j<=js2;j++){
          jterm = (j-js1)*nx;

          if(smokecullflag==1&&j!=js2){
            x11[1]=yplt[j];
            x12[1]=yplt[j];
            x22[1]=yplt[j+1];
            x21[1]=yplt[j+1];
            if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
          }

          for(i=is1;i<=is2;i++){
            iterm = (i-is1);
            n = iterm + jterm + kterm;
            ASSERT(n>=0&&n<smoke3di->nchars_uncompressed);
            mergealpha = mergealphaptr[n];
            ADJUSTALPHA(mergealpha,aspectratio,NULL,3);
          }
        }
      }
      alphaf_ptr=alphaf_out;
    }
    else{
      alphaf_ptr=alphaf_in;
    }

    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);
    if(smokedrawtest==1||show_smoketest==1){
      slice_end=smokedrawtest_nummax;
      if(slice_end>ks2)slice_end=ks2;
      slice_beg=smokedrawtest_nummin;
      if(slice_beg<ks1)slice_beg=ks1;
      if(show_smoketest==1){
        slice_end=ks2-1;
        slice_beg=slice_end-1;
      }
    }
    else{
      slice_beg=ks1;
      slice_end=ks2;
    }
    for(kkk=slice_beg;kkk<slice_end;kkk+=skip_local){
      k=kkk;
      if(ssmokedir<0)k = ks1+ks2-kkk-1;
      constval = zplt[k]+0.001;
      kterm = (k-ks1)*nxy;

      if(smokecullflag==1){
        x11[0]=xplt[is1];
        x12[0]=xplt[is2];
        x22[0]=xplt[is2];
        x21[0]=xplt[is1];

        x11[1]=yplt[js1];
        x12[1]=yplt[js1];
        x22[1]=yplt[js2];
        x21[1]=yplt[js2];

        x11[2]=zplt[k];
        x12[2]=zplt[k];
        x22[2]=zplt[k];
        x21[2]=zplt[k];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(j=js1; j<js2; j++){
        jterm = (j-js1)*nx;

        yy1 = yplt[j];
        y3 = yplt[j+1];

        ynode[0]=yy1;
        ynode[1]=yy1;
        ynode[2]=y3;
        ynode[3]=y3;

        if(smokecullflag==1){
          x11[1]=yplt[j];
          x12[1]=yplt[j];
          x22[1]=yplt[j+1];
          x21[1]=yplt[j+1];
          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(i=is1; i<is2; i++){
          iterm = (i-is1);
          x1 = xplt[i];
          x3 = xplt[i+1];

          xnode[0]=x1;
          xnode[1]=x3;
          xnode[2]=x3;
          xnode[3]=x1;

          n = iterm + jterm + kterm;
          n11 = n;
          n12 = n11+1;;
          n22 = n12+nx;
          n21 = n22-1;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + 1;
            m22 = m12 + nx;
            m21 = m22 - 1;
            DRAWVERTEXTERRAIN(xnode[mm],ynode[mm],constval)
          }
          else{
            DRAWVERTEX(xnode[mm],ynode[mm],constval)
          }
        }
      }
    }
    glEnd();
    break;

  // +++++++++++++++++++++++++++++++++++ DIR 4 +++++++++++++++++++++++++++++++++++++++

  case 4:
  case -4:

    // ++++++++++++++++++  adjust transparency +++++++++++++++++

    aspectratio=meshi->dxy;
    if(adjustalphaflag!=0){

      for(iii=1;iii<nx+ny-2;iii+=skip_local){
        ipj = iii;
        if(ssmokedir<0)ipj = nx+ny-2-iii;
        ibeg=0;
        jbeg=ipj;
        if(jbeg>ny-1){
          jbeg=ny-1;
          ibeg=ipj-jbeg;
        }
        iend=nx-1;
        jend=ipj-iend;
        if(jend<0){
          jend=0;
          iend=ipj-jend;
        }


        if(smokecullflag==1){
          x11[0]=xplt[is1+ibeg];
          x12[0]=xplt[is1+iend];
          x22[0]=xplt[is1+iend];
          x21[0]=xplt[is1+ibeg];

          x11[1]=yplt[js1+jbeg];
          x12[1]=yplt[js1+jend];
          x22[1]=yplt[js1+jend];
          x21[1]=yplt[js1+jbeg];

          x11[2]=zplt[ks1];
          x12[2]=zplt[ks1];
          x22[2]=zplt[ks2];
          x21[2]=zplt[ks2];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(k=ks1;k<=ks2;k++){
          kterm = (k-ks1)*nxy;

          if(smokecullflag==1&&k!=ks2){
            x11[2]=zplt[k];
            x12[2]=zplt[k];
            x22[2]=zplt[k+1];
            x21[2]=zplt[k+1];

            if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
          }

          for(ii=ibeg;ii<=iend;ii++){
            i=is1+ii;
            iterm = (i-is1);

            jj = ipj-ii;
            j=js1+jj;
            jterm = (j-js1)*nx;

            n = iterm + jterm + kterm;
            ASSERT(n>=0&&n<smoke3di->nchars_uncompressed);
            mergealpha = mergealphaptr[n];
            ADJUSTALPHA(mergealpha,aspectratio,norm,4);
          }
        }
      }
      alphaf_ptr=alphaf_out;
    }
    else{
      alphaf_ptr=alphaf_in;
    }

    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);
    if(smokedrawtest==1||show_smoketest==1){
      slice_end=smokedrawtest_nummax;
      if(slice_end>nx+ny-2)slice_end=nx+ny-2;
      slice_beg=smokedrawtest_nummin;
      if(slice_beg<1)slice_beg=1;
      if(show_smoketest==1){
        slice_end=nx+ny-2;
        slice_beg=slice_end-1;
      }
    }
    else{
      slice_beg=1;
      slice_end=nx+ny-2;
    }
    for(iii=slice_beg;iii<slice_end;iii+=skip_local){
      ipj = iii;
      if(ssmokedir<0)ipj = nx+ny-2-iii;
      ibeg=0;
      jbeg=ipj;
      if(jbeg>ny-1){
        jbeg=ny-1;
        ibeg=ipj-jbeg;
      }
      iend=nx-1;
      jend=ipj-iend;
      if(jend<0){
        jend=0;
        iend=ipj-jend;
      }

      if(smokecullflag==1){
        x11[0]=xplt[is1+ibeg];
        x12[0]=xplt[is1+iend];
        x22[0]=xplt[is1+iend];
        x21[0]=xplt[is1+ibeg];

        x11[1]=yplt[js1+jbeg];
        x12[1]=yplt[js1+jend];
        x22[1]=yplt[js1+jend];
        x21[1]=yplt[js1+jbeg];

        x11[2]=zplt[ks1];
        x12[2]=zplt[ks1];
        x22[2]=zplt[ks2];
        x21[2]=zplt[ks2];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(k=ks1; k<ks2; k++){
        kterm = (k-ks1)*nxy;
        z1 = zplt[k];
        z3 = zplt[k+1];
        znode[0]=z1;
        znode[1]=z1;
        znode[2]=z3;
        znode[3]=z3;

        if(smokecullflag==1){
          x11[2]=zplt[k];
          x12[2]=zplt[k];
          x22[2]=zplt[k+1];
          x21[2]=zplt[k+1];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(ii=ibeg;ii<iend;ii++){
          i=is1+ii;
          iterm = (i-is1);
          x1 = xplt[i];
          x3 = xplt[i+1];

          xnode[0]=x1;
          xnode[1]=x3;
          xnode[2]=x3;
          xnode[3]=x1;

          jj = ipj-ii;
          j=js1+jj;
          jterm = (j-js1)*nx;

          yy1=yplt[j];
          y3=yplt[j-1];

          ynode[0]=yy1;
          ynode[1]=y3;
          ynode[2]=y3;
          ynode[3]=yy1;

          n11 = iterm+jterm+kterm;
          n12 = n11 - nx + 1;
          n22 = n12 + nxy;
          n21 = n11 + nxy;

//        n11 = (j-js1)*nx   + (i-is1)   + (k-ks1)*nx*ny;
//        n12 = (j-1-js1)*nx + (i+1-is1) + (k-ks1)*nx*ny;
//        n22 = (j-1-js1)*nx + (i+1-is1) + (k+1-ks1)*nx*ny;
//        n21 = (j-js1)*nx   + (i-is1)   + (k+1-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 - nx + 1;
            m22 = m12;
            m21 = m11;
            DRAWVERTEXTERRAIN(xnode[mm],ynode[mm],znode[mm])
          }
          else{
            DRAWVERTEX(xnode[mm],ynode[mm],znode[mm])
          }
        }
      }
    }
    glEnd();
    break;

  // +++++++++++++++++++++++++++++++++++ DIR 5 +++++++++++++++++++++++++++++++++++++++

    case 5:
    case -5:

    // ++++++++++++++++++  adjust transparency +++++++++++++++++

    aspectratio=meshi->dxy;
    if(adjustalphaflag!=0){

      for(iii=1;iii<nx+ny-2;iii+=skip_local){
        jmi=iii;
        if(ssmokedir<0)jmi = nx+ny-2-iii;

        ibeg=0;
        jbeg=ibeg-nx+1+jmi;
        if(jbeg<0){
          jbeg=0;
          ibeg=jbeg+nx-1-jmi;
        }
        iend=nx-1;
        jend=iend+jmi+1-nx;
        if(jend>ny-1){
          jend=ny-1;
          iend=jend+nx-1-jmi;
        }

        if(smokecullflag==1){
          x11[0]=xplt[is1+ibeg];
          x12[0]=xplt[is1+iend];
          x22[0]=xplt[is1+iend];
          x21[0]=xplt[is1+ibeg];

          x11[1]=yplt[js1+jbeg];
          x12[1]=yplt[js1+jend];
          x22[1]=yplt[js1+jend];
          x21[1]=yplt[js1+jbeg];

          x11[2]=zplt[ks1];
          x12[2]=zplt[ks1];
          x22[2]=zplt[ks2];
          x21[2]=zplt[ks2];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(k=ks1;k<=ks2;k++){
          kterm = (k-ks1)*nxy;

          if(smokecullflag==1&&k!=ks2){
            x11[2]=zplt[k];
            x12[2]=zplt[k];
            x22[2]=zplt[k+1];
            x21[2]=zplt[k+1];
            if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
          }

          for(ii=ibeg;ii<=iend;ii++){
            i = is1 + ii;
            iterm = (i-is1);

            jj = ii + jmi + 1 - nx;
            j = js1 + jj;
            jterm = (j-js1)*nx;

            n = iterm + jterm + kterm;
            ASSERT(n>=0&&n<smoke3di->nchars_uncompressed);
            mergealpha = mergealphaptr[n];
            ADJUSTALPHA(mergealpha,aspectratio,norm,4);
          }
        }
      }
      alphaf_ptr=alphaf_out;
    }
    else{
      alphaf_ptr=alphaf_in;
    }

    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);

    if(smokedrawtest==1||show_smoketest==1){
      slice_end=smokedrawtest_nummax;
      if(slice_end>nx+ny-2)slice_end=nx+ny-2;
      slice_beg=smokedrawtest_nummin;
      if(slice_beg<1)slice_beg=1;
      if(show_smoketest==1){
        slice_end=nx+ny-2;
        slice_beg=slice_end-1;
      }
    }
    else{
      slice_beg=1;
      slice_end=nx+ny-2;
    }
    for(iii=slice_beg;iii<slice_end;iii+=skip_local){
      jmi=iii;
      if(ssmokedir<0)jmi = nx+ny-2-iii;

      ibeg=0;
      jbeg=ibeg-nx+1+jmi;
      if(jbeg<0){
        jbeg=0;
        ibeg=jbeg+nx-1-jmi;
      }
      iend=nx-1;
      jend=iend+jmi+1-nx;
      if(jend>ny-1){
        jend=ny-1;
        iend=jend+nx-1-jmi;
      }

      if(smokecullflag==1){
        x11[0]=xplt[is1+ibeg];
        x12[0]=xplt[is1+iend];
        x22[0]=xplt[is1+iend];
        x21[0]=xplt[is1+ibeg];

        x11[1]=yplt[js1+jbeg];
        x12[1]=yplt[js1+jend];
        x22[1]=yplt[js1+jend];
        x21[1]=yplt[js1+jbeg];

        x11[2]=zplt[ks1];
        x12[2]=zplt[ks1];
        x22[2]=zplt[ks2];
        x21[2]=zplt[ks2];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(k=ks1; k<ks2; k++){
        kterm = (k-ks1)*nxy;
        z1 = zplt[k];
        z3 = zplt[k+1];
        znode[0]=z1;
        znode[1]=z1;
        znode[2]=z3;
        znode[3]=z3;


        if(smokecullflag==1){
          x11[2]=z1;
          x12[2]=z1;
          x22[2]=z3;
          x21[2]=z3;

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(ii=ibeg;ii<iend;ii++){
          i = is1 + ii;
          iterm = (i-is1);

          jj = ii + jmi + 1 - nx;
          j = js1 + jj;
          jterm = (j-js1)*nx;


          yy1=yplt[j];
          y3=yplt[j+1];

          ynode[0]=yy1;
          ynode[1]=y3;
          ynode[2]=y3;
          ynode[3]=yy1;

          x1 = xplt[i];
          x3 = xplt[i+1];
          xnode[0]=x1;
          xnode[1]=x3;
          xnode[2]=x3;
          xnode[3]=x1;


          n11 = jterm + iterm + kterm;
          n12 = n11 + nx + 1;
          n22 = n12 + nxy;
          n21 = n11 + nxy;

        //    n11 = (j-js1)*nx + (i-is1) + (k-ks1)*nx*ny;
        //    n12 = (j+1-js1)*nx + (i+1-is1) + (k-ks1)*nx*ny;
        //    n22 = (j+1-js1)*nx + (i+1-is1) + (k+1-ks1)*nx*ny;
        //    n21 = (j-js1)*nx + (i-is1) + (k+1-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + nx + 1;
            m22 = m12;
            m21 = m11 - nx;
            DRAWVERTEXTERRAIN(xnode[mm], ynode[mm],znode[mm])
          }
          else{
            DRAWVERTEX(xnode[mm], ynode[mm], znode[mm])
          }
        }
      }
    }
    glEnd();
    break;

  // +++++++++++++++++++++++++++++++++++ DIR 6 +++++++++++++++++++++++++++++++++++++++

  case 6:
  case -6:

    // ++++++++++++++++++  adjust transparency +++++++++++++++++

    aspectratio=meshi->dyz;
    if(adjustalphaflag!=0){

      for(iii=1;iii<ny+nz-2;iii+=skip_local){
        jpk = iii;
        if(ssmokedir<0)jpk = ny+nz-2-iii;
        jbeg=0;
        kbeg=jpk;
        if(kbeg>nz-1){
          kbeg=nz-1;
          jbeg=jpk-kbeg;
        }
        jend=ny-1;
        kend=jpk-jend;
        if(kend<0){
          kend=0;
          jend=jpk-kend;
        }


        if(smokecullflag==1){
          x11[0]=xplt[is1];
          x12[0]=xplt[is1];
          x22[0]=xplt[is2];
          x21[0]=xplt[is2];

          x11[1]=yplt[js1+jbeg];
          x12[1]=yplt[js1+jend];
          x22[1]=yplt[js1+jend];
          x21[1]=yplt[js1+jbeg];

          x11[2]=zplt[ks1+kbeg];
          x12[2]=zplt[ks1+kend];
          x22[2]=zplt[ks1+kend];
          x21[2]=zplt[ks1+kbeg];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(i=is1;i<=is2;i++){
          iterm = (i-is1);

          if(smokecullflag==1&&i!=is2){
            x11[0]=xplt[i];
            x12[0]=xplt[i];
            x22[0]=xplt[i+1];
            x21[0]=xplt[i+1];

            if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
          }

          for(jj=jbeg;jj<=jend;jj++){
            j=js1+jj;
            jterm = (j-js1)*nx;

            kk = jpk-jj;
            k=ks1+kk;
            kterm = (k-ks1)*nxy;

            n = iterm + jterm + kterm;
            ASSERT(n>=0&&n<smoke3di->nchars_uncompressed);
            mergealpha = mergealphaptr[n];
            ADJUSTALPHA(mergealpha,aspectratio,norm,4);
          }
        }
      }
      alphaf_ptr=alphaf_out;
    }
    else{
      alphaf_ptr=alphaf_in;
    }

    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);
    if(smokedrawtest==1||show_smoketest==1){
      slice_end=smokedrawtest_nummax;
      if(slice_end>ny+nz-2)slice_end=ny+nz-2;
      slice_beg=smokedrawtest_nummin;
      if(slice_beg<1)slice_beg=1;
      if(show_smoketest==1){
        slice_end=ny+nz-2;
        slice_beg=slice_end-1;
      }
    }
    else{
      slice_beg=1;
      slice_end=ny+nz-2;
    }
    for(iii=slice_beg;iii<slice_end;iii+=skip_local){
      jpk = iii;
      if(ssmokedir<0)jpk = ny+nz-2-iii;
      jbeg=0;
      kbeg=jpk;
      if(kbeg>nz-1){
        kbeg=nz-1;
        jbeg=jpk-kbeg;
      }
      jend=ny-1;
      kend=jpk-jend;
      if(kend<0){
        kend=0;
        jend=jpk-kend;
      }

      if(smokecullflag==1){
        x11[0]=xplt[is1];
        x12[0]=xplt[is1];
        x22[0]=xplt[is2];
        x21[0]=xplt[is2];

        x11[1]=yplt[js1+jbeg];
        x12[1]=yplt[js1+jend];
        x22[1]=yplt[js1+jend];
        x21[1]=yplt[js1+jbeg];

        x11[2]=zplt[ks1+kbeg];
        x12[2]=zplt[ks1+kend];
        x22[2]=zplt[ks1+kend];
        x21[2]=zplt[ks1+kbeg];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(i=is1; i<is2; i++){
        iterm = (i-is1);
        x1 = xplt[i];
        x3 = xplt[i+1];
        xnode[0]=x1;
        xnode[1]=x1;
        xnode[2]=x3;
        xnode[3]=x3;

        if(smokecullflag==1){
          x11[0]=xplt[i];
          x12[0]=xplt[i];
          x22[0]=xplt[i+1];
          x21[0]=xplt[i+1];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(jj=jbeg;jj<jend;jj++){
          j=js1+jj;
          jterm = (j-js1)*nx;
          yy1 = yplt[j];
          y3 = yplt[j+1];

          ynode[0]=yy1;
          ynode[1]=y3;
          ynode[2]=y3;
          ynode[3]=yy1;

          kk = jpk-jj;
          k=ks1+kk;
          kterm = (k-ks1)*nxy;

          z1=zplt[k];
          z3=zplt[k-1];

          znode[0]=z1;
          znode[1]=z3;
          znode[2]=z3;
          znode[3]=z1;

          n11 = iterm+jterm+kterm;
          n12 = n11 + nx - nxy;
          n22 = n12 + 1;
          n21 = n22 - nx +  nxy;

//        n11 = (i-is1)   + (j-js1)*nx   + (k-ks1)*nx*ny;
//        n12 = (i-is1)   + (j+1-js1)*nx + (k-1-ks1)*nx*ny;
//        n22 = (i+1-is1) + (j+1-js1)*nx + (k-1-ks1)*nx*ny;
//        n21 = (i+1-is1) + (j-js1)*nx   + (k-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + nx;
            m22 = m12 + 1;
            m21 = m22 - nx;
            DRAWVERTEXTERRAIN(xnode[mm],ynode[mm],znode[mm])
          }
          else{
            DRAWVERTEX(xnode[mm],ynode[mm],znode[mm])
          }
        }
      }
    }
    glEnd();
    break;

  // +++++++++++++++++++++++++++++++++++ DIR 7 +++++++++++++++++++++++++++++++++++++++

    case 7:
    case -7:
    aspectratio=meshi->dyz;

    if(adjustalphaflag!=0){

      for(iii=1;iii<ny+nz-2;iii+=skip_local){
        kmj=iii;
        if(ssmokedir<0)kmj = ny+nz-2-iii;

        jbeg=0;
        kbeg=jbeg-ny+1+kmj;
        if(kbeg<0){
          kbeg=0;
          jbeg=kbeg+ny-1-kmj;
        }
        jend=ny-1;
        kend=jend+kmj+1-ny;
        if(kend>nz-1){
          kend=nz-1;
          jend=kend+ny-1-kmj;
        }

        if(smokecullflag==1){
          x11[0]=xplt[is1];
          x12[0]=xplt[is1];
          x22[0]=xplt[is2];
          x21[0]=xplt[is2];

          x11[1]=yplt[js1+jbeg];
          x12[1]=yplt[js1+jend];
          x22[1]=yplt[js1+jend];
          x21[1]=yplt[js1+jbeg];

          x11[2]=zplt[ks1+kbeg];
          x12[2]=zplt[ks1+kend];
          x22[2]=zplt[ks1+kend];
          x21[2]=zplt[ks1+kbeg];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(i=is1;i<=is2;i++){
          iterm = (i-is1);

          if(smokecullflag==1&&i!=is2){
            x11[0]=xplt[i];
            x12[0]=xplt[i];
            x22[0]=xplt[i+1];
            x21[0]=xplt[i+1];
            if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
          }

          for(jj=jbeg;jj<=jend;jj++){
            j = js1 + jj;
            jterm = (j-js1)*nx;

            kk = jj + kmj + 1 - ny;
            k = ks1 + kk;
            kterm = (k-ks1)*nxy;

            n = iterm + jterm + kterm;
            ASSERT(n>=0&&n<smoke3di->nchars_uncompressed);
            mergealpha = mergealphaptr[n];
            ADJUSTALPHA(mergealpha,aspectratio,norm,4);
          }
        }
      }
      alphaf_ptr=alphaf_out;
    }
    else{
      alphaf_ptr=alphaf_in;
    }

    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);

    if(smokedrawtest==1||show_smoketest==1){
      slice_end=smokedrawtest_nummax;
      if(slice_end>ny+nz-2)slice_end=ny+nz-2;
      slice_beg=smokedrawtest_nummin;
      if(slice_beg<1)slice_beg=1;
      if(show_smoketest==1){
        slice_end=ny+nz-2;
        slice_beg=slice_end-1;
      }
    }
    else{
      slice_beg=1;
      slice_end=ny+nz-2;
    }
    for(iii=slice_beg;iii<slice_end;iii+=skip_local){
      kmj=iii;
      if(ssmokedir<0)kmi = ny+nz-2-iii;

      jbeg=0;
      kbeg=jbeg-ny+1+kmj;
      if(kbeg<0){
        kbeg=0;
        jbeg=kbeg+ny-1-kmj;
      }
      jend=ny-1;
      kend=jend+kmj+1-ny;
      if(kend>nz-1){
        kend=nz-1;
        jend=kend+ny-1-kmj;
      }

      if(smokecullflag==1){
        x11[0]=xplt[is1];
        x12[0]=xplt[is1];
        x22[0]=xplt[is2];
        x21[0]=xplt[is2];

        x11[1]=yplt[js1+jbeg];
        x12[1]=yplt[js1+jend];
        x22[1]=yplt[js1+jend];
        x21[1]=yplt[js1+jbeg];

        x11[2]=zplt[ks1+kbeg];
        x12[2]=zplt[ks1+kend];
        x22[2]=zplt[ks1+kend];
        x21[2]=zplt[ks1+kbeg];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(i=is1; i<is2; i++){
        iterm = (i-is1);
        x1 = xplt[i];
        x3 = xplt[i+1];
        xnode[0]=x1;
        xnode[1]=x1;
        xnode[2]=x3;
        xnode[3]=x3;


        if(smokecullflag==1){
          x11[0]=x1;
          x12[0]=x1;
          x22[0]=x3;
          x21[0]=x3;

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(jj=jbeg;jj<jend;jj++){
          j = js1 + jj;
          jterm = (j-js1)*nx;

          kk = jj + kmj + 1 - ny;
          k = ks1 + kk;
          kterm = (k-ks1)*nxy;


          z1=zplt[k];
          z3=zplt[k+1];

          znode[0]=z1;
          znode[1]=z3;
          znode[2]=z3;
          znode[3]=z1;

          yy1 = yplt[j];
          y3 = yplt[j+1];
          ynode[0]=yy1;
          ynode[1]=y3;
          ynode[2]=y3;
          ynode[3]=yy1;


          n11 = jterm + iterm + kterm;
          n12 = n11 + nxy + nx;
          n22 = n12 + 1;
          n21 = n22 - nx - nxy;

        //    n11 = (i-is1)   + (j-js1)*nx    + (k-ks1)*nx*ny;
        //    n12 = (i-is1)   + (j+1-js1)*nx  + (k+1-ks1)*nx*ny;
        //    n22 = (i+1-is1) + (j+1-js1)*nx  + (k+1-ks1)*nx*ny;
        //    n21 = (i+1-is1) + (j-js1)*nx    + (k-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + nx;
            m22 = m12 + 1;
            m21 = m22 - nx;
            DRAWVERTEXTERRAIN(xnode[mm],ynode[mm],znode[mm])
          }
          else{
            DRAWVERTEX(xnode[mm], ynode[mm], znode[mm])
          }
        }
      }
    }
    glEnd();
    break;


  // +++++++++++++++++++++++++++++++++++ DIR 8 +++++++++++++++++++++++++++++++++++++++

  case 8:
  case -8:

    // ++++++++++++++++++  adjust transparency +++++++++++++++++

    aspectratio=meshi->dxz;
    if(adjustalphaflag!=0){

      for(iii=1;iii<nx+nz-2;iii+=skip_local){
        ipk = iii;
        if(ssmokedir<0)ipk = nx+nz-2-iii;
        ibeg=0;
        kbeg=ipk;
        if(kbeg>nz-1){
          kbeg=nz-1;
          ibeg=ipk-kbeg;
        }
        iend=nx-1;
        kend=ipk-iend;
        if(kend<0){
          kend=0;
          iend=ipk-kend;
        }


        if(smokecullflag==1){
          x11[0]=xplt[is1+ibeg];
          x12[0]=xplt[is1+iend];
          x22[0]=xplt[is1+iend];
          x21[0]=xplt[is1+ibeg];

          x11[1]=yplt[js1];
          x12[1]=yplt[js1];
          x22[1]=yplt[js2];
          x21[1]=yplt[js2];

          x11[2]=zplt[ks1+kbeg];
          x12[2]=zplt[ks1+kend];
          x22[2]=zplt[ks1+kend];
          x21[2]=zplt[ks1+kbeg];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(j=js1;j<=js2;j++){
          jterm = (j-js1)*nx;

          if(smokecullflag==1&&j!=js2){
            x11[1]=yplt[j];
            x12[1]=yplt[j];
            x22[1]=yplt[j+1];
            x21[1]=yplt[j+1];

            if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
          }

          for(ii=ibeg;ii<=iend;ii++){
            i=is1+ii;
            iterm = (i-is1);

            kk = ipk-ii;
            k=ks1+kk;
            kterm = (k-ks1)*nxy;

            n = iterm + jterm + kterm;
            ASSERT(n>=0&&n<smoke3di->nchars_uncompressed);
            mergealpha = mergealphaptr[n];
            ADJUSTALPHA(mergealpha,aspectratio,norm,4);
          }
        }
      }
      alphaf_ptr=alphaf_out;
    }
    else{
      alphaf_ptr=alphaf_in;
    }

    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);
    if(smokedrawtest==1||show_smoketest==1){
      slice_end=smokedrawtest_nummax;
      if(slice_end>nx+nz-2)slice_end=nx+nz-2;
      slice_beg=smokedrawtest_nummin;
      if(slice_beg<1)slice_beg=1;
      if(show_smoketest==1){
        slice_end=nx+nz-2;
        slice_beg=slice_end-1;
      }
    }
    else{
      slice_beg=1;
      slice_end=nx+nz-2;
    }
    for(iii=slice_beg;iii<slice_end;iii+=skip_local){
      ipk = iii;
      if(ssmokedir<0)ipk = nx+nz-2-iii;
      ibeg=0;
      kbeg=ipk;
      if(kbeg>nz-1){
        kbeg=nz-1;
        ibeg=ipk-kbeg;
      }
      iend=nx-1;
      kend=ipk-iend;
      if(kend<0){
        kend=0;
        iend=ipk-kend;
      }

      if(smokecullflag==1){
        x11[0]=xplt[is1+ibeg];
        x12[0]=xplt[is1+iend];
        x22[0]=xplt[is1+iend];
        x21[0]=xplt[is1+ibeg];

        x11[1]=yplt[js1];
        x12[1]=yplt[js1];
        x22[1]=yplt[js2];
        x21[1]=yplt[js2];

        x11[2]=zplt[ks1+kbeg];
        x12[2]=zplt[ks1+kend];
        x22[2]=zplt[ks1+kend];
        x21[2]=zplt[ks1+kbeg];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(j=js1; j<js2; j++){
        jterm = (j-js1)*nx;
        yy1 = yplt[j];
        y3 = yplt[j+1];
        ynode[0]=yy1;
        ynode[1]=yy1;
        ynode[2]=y3;
        ynode[3]=y3;

        if(smokecullflag==1){
          x11[1]=yplt[j];
          x12[1]=yplt[j];
          x22[1]=yplt[j+1];
          x21[1]=yplt[j+1];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(ii=ibeg;ii<iend;ii++){
          i=is1+ii;
          iterm = (i-is1);
          x1 = xplt[i];
          x3 = xplt[i+1];

          xnode[0]=x1;
          xnode[1]=x3;
          xnode[2]=x3;
          xnode[3]=x1;

          kk = ipk-ii;
          k=ks1+kk;
          kterm = (k-ks1)*nxy;

          z1=zplt[k];
          z3=zplt[k-1];

          znode[0]=z1;
          znode[1]=z3;
          znode[2]=z3;
          znode[3]=z1;

          n11 = iterm+jterm+kterm;
          n12 = n11 + 1 - nxy;
          n22 = n12 + nx;
          n21 = n22 - 1 +  nxy;

//        n11 = (i-is1)   + (j-js1)*nx   + (k-ks1)*nx*ny;
//        n12 = (i+1-is1) + (j-js1)*nx   + (k-1-ks1)*nx*ny;
//        n22 = (i+1-is1) + (j+1-js1)*nx + (k-1-ks1)*nx*ny;
//        n21 = (i-is1)   + (j+1-js1)*nx + (k-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + 1;
            m22 = m12 + nx;
            m21 = m22 - 1;
            DRAWVERTEXTERRAIN(xnode[mm],ynode[mm],znode[mm])
          }
          else{
            DRAWVERTEX(xnode[mm],ynode[mm],znode[mm])
          }
        }
      }
    }
    glEnd();
    break;

  // +++++++++++++++++++++++++++++++++++ DIR 9 +++++++++++++++++++++++++++++++++++++++

      /* interchange y and z, j and z */
    case 9:
    case -9:
    aspectratio=meshi->dxz;
    if(adjustalphaflag!=0){

      for(iii=1;iii<nx+nz-2;iii+=skip_local){
        kmi=iii;
        if(ssmokedir<0)kmi = nx+nz-2-iii;

        ibeg=0;
        kbeg=ibeg-nx+1+kmi;
        if(kbeg<0){
          kbeg=0;
          ibeg=kbeg+nx-1-kmi;
        }
        iend=nx-1;
        kend=iend+kmi+1-nx;
        if(kend>nz-1){
          kend=nz-1;
          iend=kend+nx-1-kmi;
        }

        if(smokecullflag==1){
          x11[0]=xplt[is1+ibeg];
          x12[0]=xplt[is1+iend];
          x22[0]=xplt[is1+iend];
          x21[0]=xplt[is1+ibeg];

          x11[1]=yplt[js1];
          x12[1]=yplt[js1];
          x22[1]=yplt[js2];
          x21[1]=yplt[js2];

          x11[2]=zplt[ks1+kbeg];
          x12[2]=zplt[ks1+kend];
          x22[2]=zplt[ks1+kend];
          x21[2]=zplt[ks1+kbeg];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(j=js1;j<=js2;j++){
          jterm = (j-js1)*nx;

          if(smokecullflag==1&&j!=js2){
            x11[1]=yplt[j];
            x12[1]=yplt[j];
            x22[1]=yplt[j+1];
            x21[1]=yplt[j+1];
            if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
          }

          for(ii=ibeg;ii<=iend;ii++){
            i = is1 + ii;
            iterm = (i-is1);

            kk = ii + kmi + 1 - nx;
            k = ks1 + kk;
            kterm = (k-ks1)*nxy;

            n = iterm + jterm + kterm;
            ASSERT(n>=0&&n<smoke3di->nchars_uncompressed);
            mergealpha = mergealphaptr[n];
            ADJUSTALPHA(mergealpha,aspectratio,norm,4);
          }
        }
      }
      alphaf_ptr=alphaf_out;
    }
    else{
      alphaf_ptr=alphaf_in;
    }

    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);

    if(smokedrawtest==1||show_smoketest==1){
      slice_end=smokedrawtest_nummax;
      if(slice_end>nx+nz-2)slice_end=nx+nz-2;
      slice_beg=smokedrawtest_nummin;
      if(slice_beg<1)slice_beg=1;
      if(show_smoketest==1){
        slice_end=nx+nz-2;
        slice_beg=slice_end-1;
      }
    }
    else{
      slice_beg=1;
      slice_end=nx+nz-2;
    }
    for(iii=slice_beg;iii<slice_end;iii+=skip_local){
      kmi=iii;
      if(ssmokedir<0)kmi = nx+nz-2-iii;

      ibeg=0;
      kbeg=ibeg-nx+1+kmi;
      if(kbeg<0){
        kbeg=0;
        ibeg=kbeg+nx-1-kmi;
      }
      iend=nx-1;
      kend=iend+kmi+1-nx;
      if(kend>nz-1){
        kend=nz-1;
        iend=kend+nx-1-kmi;
      }

      if(smokecullflag==1){
        x11[0]=xplt[is1+ibeg];
        x12[0]=xplt[is1+iend];
        x22[0]=xplt[is1+iend];
        x21[0]=xplt[is1+ibeg];

        x11[1]=yplt[js1];
        x12[1]=yplt[js1];
        x22[1]=yplt[js2];
        x21[1]=yplt[js2];

        x11[2]=zplt[ks1+kbeg];
        x12[2]=zplt[ks1+kend];
        x22[2]=zplt[ks1+kend];
        x21[2]=zplt[ks1+kbeg];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(j=js1; j<js2; j++){
        jterm = (j-js1)*nx;
        yy1 = yplt[j];
        y3 = yplt[j+1];
        ynode[0]=yy1;
        ynode[1]=yy1;
        ynode[2]=y3;
        ynode[3]=y3;


        if(smokecullflag==1){
          x11[1]=yy1;
          x12[1]=yy1;
          x22[1]=y3;
          x21[1]=y3;

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(ii=ibeg;ii<iend;ii++){
          i = is1 + ii;
          iterm = (i-is1);

          kk = ii + kmi + 1 - nx;
          k = ks1 + kk;
          kterm = (k-ks1)*nxy;


          z1=zplt[k];
          z3=zplt[k+1];

          znode[0]=z1;
          znode[1]=z3;
          znode[2]=z3;
          znode[3]=z1;

          x1 = xplt[i];
          x3 = xplt[i+1];
          xnode[0]=x1;
          xnode[1]=x3;
          xnode[2]=x3;
          xnode[3]=x1;


          n11 = jterm + iterm + kterm;
          n12 = n11 + nxy + 1;
          n22 = n12 + nx;
          n21 = n22 - 1 - nxy;

        //    n11 = (i-is1)   + (j-js1)*nx   + (k-ks1)*nx*ny;
        //    n12 = (i+1-is1) + (j-js1)*nx   + (k+1-ks1)*nx*ny;
        //    n22 = (i+1-is1) + (j+1-js1)*nx + (k+1-ks1)*nx*ny;
        //    n21 = (i-is1)   + (j+1-js1)*nx + (k-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + 1;
            m22 = m12 + nx;
            m21 = m22 - 1;
            DRAWVERTEXTERRAIN(xnode[mm],ynode[mm],znode[mm])
          }
          else{
            DRAWVERTEX(xnode[mm], ynode[mm], znode[mm])
          }
        }
      }
    }
    glEnd();
    break;
    default:
      ASSERT(FFALSE);
      break;
  }
  TransparentOff();
  if(cullfaces==1)glEnable(GL_CULL_FACE);
//  PRINTF("majorcull=%i minorcull=%i\n",majorcull,minorcull);
}

#ifdef pp_GPU
/* ------------------ drawsmoke3dGPU ------------------------ */

void drawsmoke3dGPU(smoke3ddata *smoke3di){
  int i,j,k,n;
  float constval,x1,x3,z1,z3, yy1, y3;
  int is1, is2, js1, js2, ks1, ks2;
  int ii, jj, kk;
  int ibeg, iend, jbeg, jend, kbeg, kend;

  float *xplt, *yplt, *zplt;
  float *znode_offset, z_offset[4];

  int nx,ny,nz;
  int xyzindex1[6],xyzindex2[6],*xyzindex,node,mm;
  float xnode[4],znode[4],ynode[4];
  int skip_local;
  int iterm, jterm, kterm,nxy;
  float x11[3], x12[3], x22[3], x21[3];
  int n11, n12, n22, n21;
  int ipj,jpk,ipk,jmi,kmi,kmj;
  int iii, jjj, kkk;
  int slice_end,slice_beg;
  float aspectratio;
  int ssmokedir;
  unsigned char *iblank_smoke3d;
  int have_smoke;

  unsigned char *firecolor, *alphaf_in;
  float value[4], fvalue[4];

  meshdata *meshi;

  meshi = meshinfo + smoke3di->blocknumber;
  if(meshvisptr[meshi-meshinfo]==0)return;

  firecolor=smoke3di->hrrpuv_color;

  {
    smoke3ddata *sooti=NULL;

	if(smoke3di->soot_index>=0){
	  sooti = smoke3dinfo + smoke3di->soot_index;
	}
	if(sooti!=NULL&&sooti->display==1){
      have_smoke=1;
    }
    else{
      have_smoke=0;
    }
  }
  iblank_smoke3d = meshi->iblank_smoke3d;

  // meshi->hrrpuv_cutoff
  // hrrpuv_max_smv;

  meshi = meshinfo + smoke3di->blocknumber;
  if(meshvisptr[meshi-meshinfo]==0)return;
  value[0]=255;
  value[1]=255;
  value[2]=255;
  value[3]=255;

  if(nterraininfo>0){
    znode_offset = meshi->terrain->znode_offset;
  }

  xplt=meshi->xplt;
  yplt=meshi->yplt;
  zplt=meshi->zplt;
  alphaf_in=smoke3di->smokeframe_in;

  is1 = smoke3di->is1;
  is2 = smoke3di->is2;
  js1 = smoke3di->js1;
  js2 = smoke3di->js2;
  ks1 = smoke3di->ks1;
  ks2 = smoke3di->ks2;

  nx = smoke3di->is2 + 1 - smoke3di->is1;
  ny = js2 + 1 - js1;
  nz = ks2 + 1 - ks1;
  nxy = nx*ny;

  ssmokedir=meshi->smokedir;
  skip_local=smokeskipm1+1;

  xyzindex1[0]=0;
  xyzindex1[1]=1;
  xyzindex1[2]=2;
  xyzindex1[3]=0;
  xyzindex1[4]=2;
  xyzindex1[5]=3;

  xyzindex2[0]=0;
  xyzindex2[1]=1;
  xyzindex2[2]=3;
  xyzindex2[3]=1;
  xyzindex2[4]=2;
  xyzindex2[5]=3;

  if(cullfaces==1)glDisable(GL_CULL_FACE);

  glUniform1i(GPU_adjustalphaflag,adjustalphaflag);
  glUniform1i(GPU_have_smoke,have_smoke);
  glUniform1i(GPU_smokecolormap,0);
  glUniform1f(GPU_smoke3d_rthick,smoke3d_rthick);
  glUniform1f(GPU_hrrpuv_max_smv,hrrpuv_max_smv);
  glUniform1f(GPU_hrrpuv_cutoff,global_hrrpuv_cutoff);
  glUniform1f(GPU_fire_alpha, smoke3di->fire_alpha);
  glUniform1f(GPU_smoke_albedo, smoke_albedo);

  TransparentOn();
  switch(ssmokedir){
  // +++++++++++++++++++++++++++++++++++ DIR 1 +++++++++++++++++++++++++++++++++++++++


  case 1:
  case -1:

    aspectratio=meshi->dx;
    glUniform1f(GPU_aspectratio,aspectratio);

    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);
    slice_beg=is1;
    slice_end=is2;
    for(iii=slice_beg;iii<slice_end;iii+=skip_local){
      i=iii;
      if(ssmokedir<0)i = is1+is2-iii-1;
      iterm = (i-smoke3di->is1);

      if(smokecullflag==1){
        x11[0]=xplt[i];
        x12[0]=xplt[i];
        x22[0]=xplt[i];
        x21[0]=xplt[i];

        x11[1]=yplt[js1];
        x12[1]=yplt[js2];
        x22[1]=yplt[js2];
        x21[1]=yplt[js1];

        x11[2]=zplt[ks1];
        x12[2]=zplt[ks1];
        x22[2]=zplt[ks2];
        x21[2]=zplt[ks2];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      constval = xplt[i]+0.001;
      for(k=ks1; k<ks2; k++){
        kterm = (k-ks1)*nxy;
        z1 = zplt[k];
        z3 = zplt[k+1];
        znode[0]=z1;
        znode[1]=z1;
        znode[2]=z3;
        znode[3]=z3;

        if(smokecullflag==1&&k!=ks2){
          x11[2]=zplt[k];
          x12[2]=zplt[k];
          x22[2]=zplt[k+1];
          x21[2]=zplt[k+1];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }
        for(j=js1; j<js2; j++){
          jterm = (j-js1)*nx;
          yy1 = yplt[j];
          y3 = yplt[j+1];
          ynode[0]=yy1;
          ynode[1]=y3;
          ynode[2]=y3;
          ynode[3]=yy1;

          n = iterm + jterm + kterm;

          n11 = n;              //n
          n12 = n11 + nx;       //n+nx
          n22 = n12 + nxy;      //n+nx+nxy
          n21 = n22 - nx;       //n+nxy

//        n11 = (i-is1)   + (j-js1)*nx   + (k-ks1)*nx*ny;
//        n12 = (i-is1)   + (j+1-js1)*nx + (k-ks1)*nx*ny;
//        n22 = (i-is1)   + (j+1-js1)*nx + (k+1-ks1)*nx*ny;
//        n21 = (i-is1)   + (j-js1)*nx   + (k+1-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + nx;
            m22 = m12;
            m21 = m22 - nx;
            DRAWVERTEXGPUTERRAIN(constval,ynode[mm],znode[mm])
          }
          else{
            DRAWVERTEXGPU(constval,ynode[mm],znode[mm])
          }

        }
      }
    }
    glEnd();

    break;

  // +++++++++++++++++++++++++++++++++++ DIR 2 +++++++++++++++++++++++++++++++++++++++

  case 2:
  case -2:

    aspectratio=meshi->dy;
    glUniform1f(GPU_aspectratio,aspectratio);


  // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);
    slice_beg=js1;
    slice_end=js2;
    for(jjj=slice_beg;jjj<slice_end;jjj+=skip_local){
      j=jjj;
      if(ssmokedir<0)j = js1+js2-jjj-1;
      constval = yplt[j]+0.001;
      jterm = (j-js1)*nx;

      if(smokecullflag==1){
        x11[0]=xplt[is1];
        x12[0]=xplt[is2];
        x22[0]=xplt[is2];
        x21[0]=xplt[is1];

        x11[1]=yplt[j];
        x12[1]=yplt[j];
        x22[1]=yplt[j];
        x21[1]=yplt[j];

        x11[2]=zplt[ks1];
        x12[2]=zplt[ks1];
        x22[2]=zplt[ks2];
        x21[2]=zplt[ks2];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(k=ks1; k<ks2; k++){
        kterm = (k-ks1)*nxy;
        z1 = zplt[k];
        z3 = zplt[k+1];

        znode[0]=z1;
        znode[1]=z1;
        znode[2]=z3;
        znode[3]=z3;

        if(smokecullflag==1){
          x11[2]=zplt[k];
          x12[2]=zplt[k];
          x22[2]=zplt[k+1];
          x21[2]=zplt[k+1];
          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(i=is1; i<is2; i++){
          iterm = (i-is1);
          x1 = xplt[i];
          x3 = xplt[i+1];

          xnode[0]=x1;
          xnode[1]=x3;
          xnode[2]=x3;
          xnode[3]=x1;

          n = iterm + jterm + kterm;
          n11 = n;            //n
          n12 = n11+1;;       //n+1
          n22 = n12+nxy;      //n+1+nxy
          n21 = n22-1;        //n+nxy

//        n11 = (i-is1)   + (j-js1)*nx   + (k-ks1)*nx*ny;
//        n12 = (i+1-is1) + (j-js1)*nx   + (k-ks1)*nx*ny;
//        n22 = (i+1-is1) + (j-js1)*nx   + (k+1-ks1)*nx*ny;
//        n21 = (i-is1)   + (j-js1)*nx   + (k+1-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + 1;
            m22 = m12;
            m21 = m22 - 1;
            DRAWVERTEXGPUTERRAIN(xnode[mm],constval,znode[mm])
          }
          else{
            DRAWVERTEXGPU(xnode[mm],constval,znode[mm])
          }

        }
      }
    }
    glEnd();
    break;

  // +++++++++++++++++++++++++++++++++++ DIR 3 +++++++++++++++++++++++++++++++++++++++

  case 3:
  case -3:

    aspectratio=meshi->dz;
    glUniform1f(GPU_aspectratio,aspectratio);


    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);
    slice_beg=ks1;
    slice_end=ks2;
    for(kkk=slice_beg;kkk<slice_end;kkk+=skip_local){
      k=kkk;
      if(ssmokedir<0)k = ks1+ks2-kkk-1;
      constval = zplt[k]+0.001;
      kterm = (k-ks1)*nxy;

      if(smokecullflag==1){
        x11[0]=xplt[is1];
        x12[0]=xplt[is2];
        x22[0]=xplt[is2];
        x21[0]=xplt[is1];

        x11[1]=yplt[js1];
        x12[1]=yplt[js1];
        x22[1]=yplt[js2];
        x21[1]=yplt[js2];

        x11[2]=zplt[k];
        x12[2]=zplt[k];
        x22[2]=zplt[k];
        x21[2]=zplt[k];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(j=js1; j<js2; j++){
        jterm = (j-js1)*nx;

        yy1 = yplt[j];
        y3 = yplt[j+1];

        ynode[0]=yy1;
        ynode[1]=yy1;
        ynode[2]=y3;
        ynode[3]=y3;

        if(smokecullflag==1){
          x11[1]=yplt[j];
          x12[1]=yplt[j];
          x22[1]=yplt[j+1];
          x21[1]=yplt[j+1];
          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(i=is1; i<is2; i++){
          iterm = (i-is1);
          x1 = xplt[i];
          x3 = xplt[i+1];

          xnode[0]=x1;
          xnode[1]=x3;
          xnode[2]=x3;
          xnode[3]=x1;

          n = iterm + jterm + kterm;
          n11 = n;
          n12 = n11+1;;
          n22 = n12+nx;
          n21 = n22-1;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + 1;
            m22 = m12 + nx;
            m21 = m22 - 1;
            DRAWVERTEXGPUTERRAIN(xnode[mm],ynode[mm],constval)
          }
          else{
            DRAWVERTEXGPU(xnode[mm],ynode[mm],constval)
          }
        }
      }
    }
    glEnd();
    break;

  // +++++++++++++++++++++++++++++++++++ DIR 4 +++++++++++++++++++++++++++++++++++++++

  case 4:
  case -4:

    aspectratio=meshi->dxy;
    glUniform1f(GPU_aspectratio,aspectratio);

    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);
    slice_beg=1;
    slice_end=nx+ny-2;
    for(iii=slice_beg;iii<slice_end;iii+=skip_local){
      ipj = iii;
      if(ssmokedir<0)ipj = nx+ny-2-iii;
      ibeg=0;
      jbeg=ipj;
      if(jbeg>ny-1){
        jbeg=ny-1;
        ibeg=ipj-jbeg;
      }
      iend=nx-1;
      jend=ipj-iend;
      if(jend<0){
        jend=0;
        iend=ipj-jend;
      }

      if(smokecullflag==1){
        x11[0]=xplt[is1+ibeg];
        x12[0]=xplt[is1+iend];
        x22[0]=xplt[is1+iend];
        x21[0]=xplt[is1+ibeg];

        x11[1]=yplt[js1+jbeg];
        x12[1]=yplt[js1+jend];
        x22[1]=yplt[js1+jend];
        x21[1]=yplt[js1+jbeg];

        x11[2]=zplt[ks1];
        x12[2]=zplt[ks1];
        x22[2]=zplt[ks2];
        x21[2]=zplt[ks2];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(k=ks1; k<ks2; k++){
        kterm = (k-ks1)*nxy;
        z1 = zplt[k];
        z3 = zplt[k+1];
        znode[0]=z1;
        znode[1]=z1;
        znode[2]=z3;
        znode[3]=z3;

        if(smokecullflag==1){
          x11[2]=zplt[k];
          x12[2]=zplt[k];
          x22[2]=zplt[k+1];
          x21[2]=zplt[k+1];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(ii=ibeg;ii<iend;ii++){
          i=is1+ii;
          iterm = (i-is1);
          x1 = xplt[i];
          x3 = xplt[i+1];

          xnode[0]=x1;
          xnode[1]=x3;
          xnode[2]=x3;
          xnode[3]=x1;

          jj = ipj-ii;
          j=js1+jj;
          jterm = (j-js1)*nx;

          yy1=yplt[j];
          y3=yplt[j-1];

          ynode[0]=yy1;
          ynode[1]=y3;
          ynode[2]=y3;
          ynode[3]=yy1;

          n11 = iterm+jterm+kterm;
          n12 = n11 - nx + 1;
          n22 = n12 + nxy;
          n21 = n11 + nxy;

//        n11 = (j-js1)*nx   + (i-is1)   + (k-ks1)*nx*ny;
//        n12 = (j-1-js1)*nx + (i+1-is1) + (k-ks1)*nx*ny;
//        n22 = (j-1-js1)*nx + (i+1-is1) + (k+1-ks1)*nx*ny;
//        n21 = (j-js1)*nx   + (i-is1)   + (k+1-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 - nx + 1;
            m22 = m12;
            m21 = m11 - nx;
            DRAWVERTEXGPUTERRAIN(xnode[mm],ynode[mm],znode[mm])
          }
          else{
            DRAWVERTEXGPU(xnode[mm],ynode[mm],znode[mm])
          }
        }
      }
    }
    glEnd();
    break;

  // +++++++++++++++++++++++++++++++++++ DIR 5 +++++++++++++++++++++++++++++++++++++++

    case 5:
    case -5:

    aspectratio=meshi->dxy;
    glUniform1f(GPU_aspectratio,aspectratio);


    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);

    slice_beg=1;
    slice_end=nx+ny-2;
    for(iii=slice_beg;iii<slice_end;iii+=skip_local){
      jmi=iii;
      if(ssmokedir<0)jmi = nx+ny-2-iii;

      ibeg=0;
      jbeg=ibeg-nx+1+jmi;
      if(jbeg<0){
        jbeg=0;
        ibeg=jbeg+nx-1-jmi;
      }
      iend=nx-1;
      jend=iend+jmi+1-nx;
      if(jend>ny-1){
        jend=ny-1;
        iend=jend+nx-1-jmi;
      }

      if(smokecullflag==1){
        x11[0]=xplt[is1+ibeg];
        x12[0]=xplt[is1+iend];
        x22[0]=xplt[is1+iend];
        x21[0]=xplt[is1+ibeg];

        x11[1]=yplt[js1+jbeg];
        x12[1]=yplt[js1+jend];
        x22[1]=yplt[js1+jend];
        x21[1]=yplt[js1+jbeg];

        x11[2]=zplt[ks1];
        x12[2]=zplt[ks1];
        x22[2]=zplt[ks2];
        x21[2]=zplt[ks2];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(k=ks1; k<ks2; k++){
        kterm = (k-ks1)*nxy;
        z1 = zplt[k];
        z3 = zplt[k+1];
        znode[0]=z1;
        znode[1]=z1;
        znode[2]=z3;
        znode[3]=z3;


        if(smokecullflag==1){
          x11[2]=z1;
          x12[2]=z1;
          x22[2]=z3;
          x21[2]=z3;

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(ii=ibeg;ii<iend;ii++){
          i = is1 + ii;
          iterm = (i-is1);

          jj = ii + jmi + 1 - nx;
          j = js1 + jj;
          jterm = (j-js1)*nx;


          yy1=yplt[j];
          y3=yplt[j+1];

          ynode[0]=yy1;
          ynode[1]=y3;
          ynode[2]=y3;
          ynode[3]=yy1;

          x1 = xplt[i];
          x3 = xplt[i+1];
          xnode[0]=x1;
          xnode[1]=x3;
          xnode[2]=x3;
          xnode[3]=x1;


          n11 = jterm + iterm + kterm;
          n12 = n11 + nx + 1;
          n22 = n12 + nxy;
          n21 = n11 + nxy;

        //    n11 = (j-js1)*nx + (i-is1) + (k-ks1)*nx*ny;
        //    n12 = (j+1-js1)*nx + (i+1-is1) + (k-ks1)*nx*ny;
        //    n22 = (j+1-js1)*nx + (i+1-is1) + (k+1-ks1)*nx*ny;
        //    n21 = (j-js1)*nx + (i-is1) + (k+1-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + nx + 1;
            m22 = m12;
            m21 = m11;
            DRAWVERTEXGPUTERRAIN(xnode[mm],ynode[mm],znode[mm])
          }
          else{
            DRAWVERTEXGPU(xnode[mm], ynode[mm], znode[mm])
          }
        }
      }
    }
    glEnd();
    break;

  // +++++++++++++++++++++++++++++++++++ DIR 6 +++++++++++++++++++++++++++++++++++++++

  case 6:
  case -6:

    aspectratio=meshi->dyz;
    glUniform1f(GPU_aspectratio,aspectratio);


    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);
    slice_beg=1;
    slice_end=ny+nz-2;
    for(iii=slice_beg;iii<slice_end;iii+=skip_local){
      jpk = iii;
      if(ssmokedir<0)jpk = ny+nz-2-iii;
      jbeg=0;
      kbeg=jpk;
      if(kbeg>nz-1){
        kbeg=nz-1;
        jbeg=jpk-kbeg;
      }
      jend=ny-1;
      kend=jpk-jend;
      if(kend<0){
        kend=0;
        jend=jpk-kend;
      }

      if(smokecullflag==1){
        x11[0]=xplt[is1];
        x12[0]=xplt[is1];
        x22[0]=xplt[is2];
        x21[0]=xplt[is2];

        x11[1]=yplt[js1+jbeg];
        x12[1]=yplt[js1+jend];
        x22[1]=yplt[js1+jend];
        x21[1]=yplt[js1+jbeg];

        x11[2]=zplt[ks1+kbeg];
        x12[2]=zplt[ks1+kend];
        x22[2]=zplt[ks1+kend];
        x21[2]=zplt[ks1+kbeg];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(i=is1; i<is2; i++){
        iterm = (i-is1);
        x1 = xplt[i];
        x3 = xplt[i+1];
        xnode[0]=x1;
        xnode[1]=x1;
        xnode[2]=x3;
        xnode[3]=x3;

        if(smokecullflag==1){
          x11[0]=xplt[i];
          x12[0]=xplt[i];
          x22[0]=xplt[i+1];
          x21[0]=xplt[i+1];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(jj=jbeg;jj<jend;jj++){
          j=js1+jj;
          jterm = (j-js1)*nx;
          yy1 = yplt[j];
          y3 = yplt[j+1];

          ynode[0]=yy1;
          ynode[1]=y3;
          ynode[2]=y3;
          ynode[3]=yy1;

          kk = jpk-jj;
          k=ks1+kk;
          kterm = (k-ks1)*nxy;

          z1=zplt[k];
          z3=zplt[k-1];

          znode[0]=z1;
          znode[1]=z3;
          znode[2]=z3;
          znode[3]=z1;

          n11 = iterm+jterm+kterm;
          n12 = n11 + nx - nxy;
          n22 = n12 + 1;
          n21 = n22 - nx +  nxy;

//        n11 = (i-is1)   + (j-js1)*nx   + (k-ks1)*nx*ny;
//        n12 = (i-is1)   + (j+1-js1)*nx + (k-1-ks1)*nx*ny;
//        n22 = (i+1-is1) + (j+1-js1)*nx + (k-1-ks1)*nx*ny;
//        n21 = (i+1-is1) + (j-js1)*nx   + (k-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + nx;
            m22 = m12 + 1;
            m21 = m22 - nx;
            DRAWVERTEXGPUTERRAIN(xnode[mm],ynode[mm],znode[mm])
          }
          else{
            DRAWVERTEXGPU(xnode[mm],ynode[mm],znode[mm])
          }
        }
      }
    }
    glEnd();
    break;

  // +++++++++++++++++++++++++++++++++++ DIR 7 +++++++++++++++++++++++++++++++++++++++

    case 7:
    case -7:
    aspectratio=meshi->dyz;
    glUniform1f(GPU_aspectratio,aspectratio);


    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);

    slice_beg=1;
    slice_end=ny+nz-2;
    for(iii=slice_beg;iii<slice_end;iii+=skip_local){
      kmj=iii;
      if(ssmokedir<0)kmi = ny+nz-2-iii;

      jbeg=0;
      kbeg=jbeg-ny+1+kmj;
      if(kbeg<0){
        kbeg=0;
        jbeg=kbeg+ny-1-kmj;
      }
      jend=ny-1;
      kend=jend+kmj+1-ny;
      if(kend>nz-1){
        kend=nz-1;
        jend=kend+ny-1-kmj;
      }

      if(smokecullflag==1){
        x11[0]=xplt[is1];
        x12[0]=xplt[is1];
        x22[0]=xplt[is2];
        x21[0]=xplt[is2];

        x11[1]=yplt[js1+jbeg];
        x12[1]=yplt[js1+jend];
        x22[1]=yplt[js1+jend];
        x21[1]=yplt[js1+jbeg];

        x11[2]=zplt[ks1+kbeg];
        x12[2]=zplt[ks1+kend];
        x22[2]=zplt[ks1+kend];
        x21[2]=zplt[ks1+kbeg];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(i=is1; i<is2; i++){
        iterm = (i-is1);
        x1 = xplt[i];
        x3 = xplt[i+1];
        xnode[0]=x1;
        xnode[1]=x1;
        xnode[2]=x3;
        xnode[3]=x3;


        if(smokecullflag==1){
          x11[0]=x1;
          x12[0]=x1;
          x22[0]=x3;
          x21[0]=x3;

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(jj=jbeg;jj<jend;jj++){
          j = js1 + jj;
          jterm = (j-js1)*nx;

          kk = jj + kmj + 1 - ny;
          k = ks1 + kk;
          kterm = (k-ks1)*nxy;


          z1=zplt[k];
          z3=zplt[k+1];

          znode[0]=z1;
          znode[1]=z3;
          znode[2]=z3;
          znode[3]=z1;

          yy1 = yplt[j];
          y3 = yplt[j+1];
          ynode[0]=yy1;
          ynode[1]=y3;
          ynode[2]=y3;
          ynode[3]=yy1;


          n11 = jterm + iterm + kterm;
          n12 = n11 + nxy + nx;
          n22 = n12 + 1;
          n21 = n22 - nx - nxy;

        //    n11 = (i-is1)   + (j-js1)*nx    + (k-ks1)*nx*ny;
        //    n12 = (i-is1)   + (j+1-js1)*nx  + (k+1-ks1)*nx*ny;
        //    n22 = (i+1-is1) + (j+1-js1)*nx  + (k+1-ks1)*nx*ny;
        //    n21 = (i+1-is1) + (j-js1)*nx    + (k-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + nx;
            m22 = m12 + 1;
            m21 = m22 - nx;
            DRAWVERTEXGPUTERRAIN(xnode[mm],ynode[mm],znode[mm])
          }
          else{
            DRAWVERTEXGPU(xnode[mm], ynode[mm], znode[mm])
          }
        }
      }
    }
    glEnd();
    break;


  // +++++++++++++++++++++++++++++++++++ DIR 8 +++++++++++++++++++++++++++++++++++++++

  case 8:
  case -8:
    aspectratio=meshi->dxz;
    glUniform1f(GPU_aspectratio,aspectratio);

    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);
    slice_beg=1;
    slice_end=nx+nz-2;
    for(iii=slice_beg;iii<slice_end;iii+=skip_local){
      ipk = iii;
      if(ssmokedir<0)ipk = nx+nz-2-iii;
      ibeg=0;
      kbeg=ipk;
      if(kbeg>nz-1){
        kbeg=nz-1;
        ibeg=ipk-kbeg;
      }
      iend=nx-1;
      kend=ipk-iend;
      if(kend<0){
        kend=0;
        iend=ipk-kend;
      }

      if(smokecullflag==1){
        x11[0]=xplt[is1+ibeg];
        x12[0]=xplt[is1+iend];
        x22[0]=xplt[is1+iend];
        x21[0]=xplt[is1+ibeg];

        x11[1]=yplt[js1];
        x12[1]=yplt[js1];
        x22[1]=yplt[js2];
        x21[1]=yplt[js2];

        x11[2]=zplt[ks1+kbeg];
        x12[2]=zplt[ks1+kend];
        x22[2]=zplt[ks1+kend];
        x21[2]=zplt[ks1+kbeg];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(j=js1; j<js2; j++){
        jterm = (j-js1)*nx;
        yy1 = yplt[j];
        y3 = yplt[j+1];
        ynode[0]=yy1;
        ynode[1]=yy1;
        ynode[2]=y3;
        ynode[3]=y3;

        if(smokecullflag==1){
          x11[1]=yplt[j];
          x12[1]=yplt[j];
          x22[1]=yplt[j+1];
          x21[1]=yplt[j+1];

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(ii=ibeg;ii<iend;ii++){
          i=is1+ii;
          iterm = (i-is1);
          x1 = xplt[i];
          x3 = xplt[i+1];

          xnode[0]=x1;
          xnode[1]=x3;
          xnode[2]=x3;
          xnode[3]=x1;

          kk = ipk-ii;
          k=ks1+kk;
          kterm = (k-ks1)*nxy;

          z1=zplt[k];
          z3=zplt[k-1];

          znode[0]=z1;
          znode[1]=z3;
          znode[2]=z3;
          znode[3]=z1;

          n11 = iterm+jterm+kterm;
          n12 = n11 + 1 - nxy;
          n22 = n12 + nx;
          n21 = n22 - 1 +  nxy;

//        n11 = (i-is1)   + (j-js1)*nx   + (k-ks1)*nx*ny;
//        n12 = (i+1-is1) + (j-js1)*nx   + (k-1-ks1)*nx*ny;
//        n22 = (i+1-is1) + (j+1-js1)*nx + (k-1-ks1)*nx*ny;
//        n21 = (i-is1)   + (j+1-js1)*nx + (k-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + 1;
            m22 = m12 +nx;
            m21 = m22 - 1;
            DRAWVERTEXGPUTERRAIN(xnode[mm],ynode[mm],znode[mm])
          }
          else{
            DRAWVERTEXGPU(xnode[mm],ynode[mm],znode[mm])
          }
        }
      }
    }
    glEnd();
    break;

  // +++++++++++++++++++++++++++++++++++ DIR 9 +++++++++++++++++++++++++++++++++++++++

      /* interchange y and z, j and z */
    case 9:
    case -9:
    aspectratio=meshi->dxz;
    glUniform1f(GPU_aspectratio,aspectratio);

    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);

    slice_beg=1;
    slice_end=nx+nz-2;
    for(iii=slice_beg;iii<slice_end;iii+=skip_local){
      kmi=iii;
      if(ssmokedir<0)kmi = nx+nz-2-iii;

      ibeg=0;
      kbeg=ibeg-nx+1+kmi;
      if(kbeg<0){
        kbeg=0;
        ibeg=kbeg+nx-1-kmi;
      }
      iend=nx-1;
      kend=iend+kmi+1-nx;
      if(kend>nz-1){
        kend=nz-1;
        iend=kend+nx-1-kmi;
      }

      if(smokecullflag==1){
        x11[0]=xplt[is1+ibeg];
        x12[0]=xplt[is1+iend];
        x22[0]=xplt[is1+iend];
        x21[0]=xplt[is1+ibeg];

        x11[1]=yplt[js1];
        x12[1]=yplt[js1];
        x22[1]=yplt[js2];
        x21[1]=yplt[js2];

        x11[2]=zplt[ks1+kbeg];
        x12[2]=zplt[ks1+kend];
        x22[2]=zplt[ks1+kend];
        x21[2]=zplt[ks1+kbeg];

        if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
      }

      for(j=js1; j<js2; j++){
        jterm = (j-js1)*nx;
        yy1 = yplt[j];
        y3 = yplt[j+1];
        ynode[0]=yy1;
        ynode[1]=yy1;
        ynode[2]=y3;
        ynode[3]=y3;


        if(smokecullflag==1){
          x11[1]=yy1;
          x12[1]=yy1;
          x22[1]=y3;
          x21[1]=y3;

          if(RectangleInFrustum(x11,x12,x22,x21)==0)continue;
        }

        for(ii=ibeg;ii<iend;ii++){
          i = is1 + ii;
          iterm = (i-is1);

          kk = ii + kmi + 1 - nx;
          k = ks1 + kk;
          kterm = (k-ks1)*nxy;


          z1=zplt[k];
          z3=zplt[k+1];

          znode[0]=z1;
          znode[1]=z3;
          znode[2]=z3;
          znode[3]=z1;

          x1 = xplt[i];
          x3 = xplt[i+1];
          xnode[0]=x1;
          xnode[1]=x3;
          xnode[2]=x3;
          xnode[3]=x1;


          n11 = jterm + iterm + kterm;
          n12 = n11 + nxy + 1;
          n22 = n12 + nx;
          n21 = n22 - 1 - nxy;

        //    n11 = (i-is1)   + (j-js1)*nx   + (k-ks1)*nx*ny;
        //    n12 = (i+1-is1) + (j-js1)*nx   + (k+1-ks1)*nx*ny;
        //    n22 = (i+1-is1) + (j+1-js1)*nx + (k+1-ks1)*nx*ny;
        //    n21 = (i-is1)   + (j+1-js1)*nx + (k-ks1)*nx*ny;

          if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
            int m11, m12, m22, m21;

            m11 = iterm + jterm;
            m12 = m11 + 1;
            m22 = m12 + nx;
            m21 = m22 - 1;
            DRAWVERTEXGPUTERRAIN(xnode[mm],ynode[mm],znode[mm])
          }
          else{
            DRAWVERTEXGPU(xnode[mm], ynode[mm], znode[mm])
          }
        }
      }
    }
    glEnd();
    break;
    default:
      ASSERT(FFALSE);
      break;
  }
  TransparentOff();
  if(cullfaces==1)glEnable(GL_CULL_FACE);
}

#endif
#ifdef pp_CULL

/* ------------------ init_cull_exts ------------------------ */

int init_cull_exts(void){
  char version_label[256];
  char version_label2[256];
  int i, major,  minor;
  const GLubyte *version_string;
  int err;

  cullactive=0;
  version_string=glGetString(GL_VERSION);
  if(version_string==NULL){
    fprintf(stderr,"*** Warning: GL_VERSION string is NULL in init_cull_exts()\n");
    err = 1;
    return err;
  }
  strcpy(version_label,(char *)version_string);
  strcpy(version_label2,version_label);
  for(i=0;i<strlen(version_label);i++){
    if(version_label[i]=='.')version_label[i]=' ';
  }
  sscanf(version_label,"%i %i",&major,&minor);
  if(major>1){
    cullactive=1;
    err=0;
  }
  else{
    trim_back(version_label);
    PRINTF("Smokeview is running on a system using OpenGL %s\n",version_label2);
    PRINTF("Smoke culling is not not supported, it requires OpenGL 2.0 or later.\n");
    cullsmoke=0;
    err=1;
  }
  return err;
}

/* ------------------ drawsmoke3dCULL ------------------------ */

void drawsmoke3dCULL(void){
  int i,j,k,n,nn;
  float constval,x1,x3,z1,z3, yy1, y3;
  int is1, is2, js1, js2, ks1;
  float *znode_offset, z_offset[4];

  float *xplt, *yplt, *zplt;
  int nx,ny;
  unsigned char *alphaf_in;
  int xyzindex1[6],xyzindex2[6],*xyzindex,node,mm;
  float xnode[4],znode[4],ynode[4];
  int iterm, jterm, kterm,nxy;
  int n11, n12, n22, n21;
  //int slice_end,slice_beg;
  float aspectratio;
  unsigned char *firecolor;
  unsigned char *iblank_smoke3d;
  int have_smoke;
  smoke3ddata *smoke3di;
  meshdata *meshi;
  cullplanedata *culli;

  unsigned char value[4];
  unsigned char fvalue[4];

  meshdata *mesh_old;

  CheckMemory;
  if(cullfaces==1)glDisable(GL_CULL_FACE);
  TransparentOn();

  SNIFF_ERRORS("before drawsmoke3dcull");
  value[0]=255;
  value[1]=255;
  value[2]=255;
  value[3]=255;

  xyzindex1[0]=0;
  xyzindex1[1]=1;
  xyzindex1[2]=2;
  xyzindex1[3]=0;
  xyzindex1[4]=2;
  xyzindex1[5]=3;

  xyzindex2[0]=0;
  xyzindex2[1]=1;
  xyzindex2[2]=3;
  xyzindex2[3]=1;
  xyzindex2[4]=2;
  xyzindex2[5]=3;

  CheckMemory;
  SNIFF_ERRORS("in drawsmoke3dcull 1");
  glBegin(GL_TRIANGLES);
  mesh_old=NULL;
  for(nn=0;nn<ncullplaneinfo;nn++){
    culli = sort_cullplaneinfo[ncullplaneinfo-1-nn];

    meshi = culli->cull_mesh;
    if(meshvisptr[meshi-meshinfo]==0)continue;

    if(meshi!=mesh_old){
      mesh_old=meshi;

      glEnd();

      if(nterraininfo>0){
        znode_offset = meshi->terrain->znode_offset;
      }
      SNIFF_ERRORS("in drawsmoke3dcull 4");
      smoke3di=meshi->cull_smoke3d;

      firecolor=smoke3di->hrrpuv_color;

      xplt=meshi->xplt;
      yplt=meshi->yplt;
      zplt=meshi->zplt;
      alphaf_in=smoke3di->smokeframe_in;
      iblank_smoke3d = meshi->iblank_smoke3d;

      is1 = smoke3di->is1;
      is2 = smoke3di->is2;
      js1 = smoke3di->js1;
      js2 = smoke3di->js2;
      ks1 = smoke3di->ks1;

      nx = is2 + 1 - is1;
      ny = js2 + 1 - js1;
      nxy = nx*ny;

      switch(meshi->smokedir){
        case 1:
        case -1:
          aspectratio=meshi->dx;
          break;
        case 2:
        case -2:
          aspectratio=meshi->dy;
          break;
        case 3:
        case -3:
          aspectratio=meshi->dz;
          break;
        case 4:
        case -4:
        case 5:
        case -5:
          aspectratio=meshi->dxy;
          break;
        case 6:
        case -6:
        case 7:
        case -7:
          aspectratio = meshi->dyz;
          break;
        case 8:
        case -8:
        case 9:
        case -9:
          aspectratio = meshi->dxz;
          break;
        default:
          ASSERT(FFALSE);
          break;
        }
	    {
	      smoke3ddata *sooti=NULL;

        if(smoke3di->soot_index>=0){
          sooti = smoke3dinfo + smoke3di->soot_index;
        }
        if(sooti!=NULL&&sooti->display==1){
          have_smoke=1;
        }
        else{
          have_smoke=0;
        }
	    }
      glUniform1i(GPU_adjustalphaflag,adjustalphaflag);
      glUniform1i(GPU_have_smoke,have_smoke);
      glUniform1i(GPU_smokecolormap,0);
      glUniform1f(GPU_smoke3d_rthick,smoke3d_rthick);
      glUniform1f(GPU_hrrpuv_max_smv,hrrpuv_max_smv);
      glUniform1f(GPU_hrrpuv_cutoff,global_hrrpuv_cutoff);
      glUniform1f(GPU_aspectratio, aspectratio);
      glUniform1f(GPU_fire_alpha,smoke3di->fire_alpha);
      glUniform1f(GPU_smoke_albedo, smoke_albedo);

      glBegin(GL_TRIANGLES);
    }

    if(is_fire_or_soot(smoke3di)==0)continue;
    switch(meshi->smokedir){

  // +++++++++++++++++++++++++++++++++++ DIR 1 +++++++++++++++++++++++++++++++++++++++

      case 1:
      case -1:
        iterm = (culli->ibeg-is1);
        constval = xplt[culli->ibeg]+0.001;
        for(k=culli->kbeg; k<culli->kend; k++){
          kterm = (k-ks1)*nxy;
          z1 = zplt[k];
          z3 = zplt[k+1];
          znode[0]=z1;
          znode[1]=z1;
          znode[2]=z3;
          znode[3]=z3;

          for(j=culli->jbeg; j<culli->jend; j++){
            jterm = (j-js1)*nx;
            yy1 = yplt[j];
            y3 = yplt[j+1];
            ynode[0]=yy1;
            ynode[1]=y3;
            ynode[2]=y3;
            ynode[3]=yy1;

            n = iterm + jterm + kterm;

            n11 = n;              //n
            n12 = n11 + nx;       //n+nx
            n22 = n12 + nxy;      //n+nx+nxy
            n21 = n22 - nx;       //n+nxy

//        n11 = (i-is1)   + (j-js1)*nx   + (k-ks1)*nx*ny;
//        n12 = (i-is1)   + (j+1-js1)*nx + (k-ks1)*nx*ny;
//        n22 = (i-is1)   + (j+1-js1)*nx + (k+1-ks1)*nx*ny;
//        n21 = (i-is1)   + (j-js1)*nx   + (k+1-ks1)*nx*ny;

            if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
              int m, m11, m12, m22, m21;

              m = iterm+jterm;
              m11 = m;
              m12 = m11 + nx;
              m22 = m12;
              m21 = m22 - nx;
              DRAWVERTEXGPUTERRAIN(constval,ynode[mm],znode[mm])
            }
            else{
              DRAWVERTEXGPU(constval,ynode[mm],znode[mm])
            }

          }
        }
        break;

  // +++++++++++++++++++++++++++++++++++ DIR 2 +++++++++++++++++++++++++++++++++++++++

      case 2:
      case -2:
        constval = yplt[culli->jbeg]+0.001;
        jterm = (culli->jbeg-js1)*nx;

        for(k=culli->kbeg; k<culli->kend; k++){
          kterm = (k-ks1)*nxy;
          z1 = zplt[k];
          z3 = zplt[k+1];

          znode[0]=z1;
          znode[1]=z1;
          znode[2]=z3;
          znode[3]=z3;

          for(i=culli->ibeg; i<culli->iend; i++){
            iterm = (i-is1);
            x1 = xplt[i];
            x3 = xplt[i+1];

            xnode[0]=x1;
            xnode[1]=x3;
            xnode[2]=x3;
            xnode[3]=x1;

            n = iterm + jterm + kterm;
            n11 = n;            //n
            n12 = n11+1;;       //n+1
            n22 = n12+nxy;      //n+1+nxy
            n21 = n22-1;        //n+nxy

//        n11 = (i-is1)   + (j-js1)*nx   + (k-ks1)*nx*ny;
//        n12 = (i+1-is1) + (j-js1)*nx   + (k-ks1)*nx*ny;
//        n22 = (i+1-is1) + (j-js1)*nx   + (k+1-ks1)*nx*ny;
//        n21 = (i-is1)   + (j-js1)*nx   + (k+1-ks1)*nx*ny;

            if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
              int m, m11, m12, m22, m21;

              m = iterm+jterm;
              m11 = m;
              m12 = m11 + 1;
              m22 = m12;
              m21 = m22 - 1;
              DRAWVERTEXGPUTERRAIN(xnode[mm],constval,znode[mm])
            }
            else{
              DRAWVERTEXGPU(xnode[mm],constval,znode[mm])
            }
          }
        }
        break;

  // +++++++++++++++++++++++++++++++++++ DIR 3 +++++++++++++++++++++++++++++++++++++++

      case 3:
      case -3:
        constval = zplt[culli->kbeg]+0.001;
        kterm = (culli->kbeg-ks1)*nxy;

        for(j=culli->jbeg; j<culli->jend; j++){
          jterm = (j-js1)*nx;

          yy1 = yplt[j];
          y3 = yplt[j+1];

          ynode[0]=yy1;
          ynode[1]=yy1;
          ynode[2]=y3;
          ynode[3]=y3;

          for(i=culli->ibeg; i<culli->iend; i++){
            iterm = (i-is1);
            x1 = xplt[i];
            x3 = xplt[i+1];

            xnode[0]=x1;
            xnode[1]=x3;
            xnode[2]=x3;
            xnode[3]=x1;

            n = iterm + jterm + kterm;
            n11 = n;
            n12 = n11+1;
            n22 = n12+nx;
            n21 = n22-1;

//        n11 = (i-is1)   + (j-js1)*nx     + (k-ks1)*nx*ny;
//        n12 = (i+1-is1) + (j-js1)*nx     + (k-ks1)*nx*ny;
//        n22 = (i+1-is1) + (j+1-js1)*nx   + (k-ks1)*nx*ny;
//        n21 = (i-is1)   + (j+1-js1)*nx   + (k-ks1)*nx*ny;

            if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
              int m, m11, m12, m22, m21;

              m = iterm+jterm;
              m11 = m;
              m12 = m11 + 1;
              m22 = m12 + nx;
              m21 = m22 - 1;
              DRAWVERTEXGPUTERRAIN(xnode[mm],ynode[mm],constval)
            }
            else{
              DRAWVERTEXGPU(xnode[mm],ynode[mm],constval)
            }
          }
        }
        break;

  // +++++++++++++++++++++++++++++++++++ DIR 4 +++++++++++++++++++++++++++++++++++++++

      case 4:
      case -4:
        for(k=culli->kbeg; k<culli->kend; k++){
          kterm = (k-ks1)*nxy;
          z1 = zplt[k];
          z3 = zplt[k+1];
          znode[0]=z1;
          znode[1]=z1;
          znode[2]=z3;
          znode[3]=z3;

          for(i=culli->ibeg;i<culli->iend;i++){
            iterm = (i-is1);
            x1 = xplt[i];
            x3 = xplt[i+1];

            xnode[0]=x1;
            xnode[1]=x3;
            xnode[2]=x3;
            xnode[3]=x1;

            j = culli->jend-(i-culli->ibeg);
            jterm = (j-js1)*nx;

            yy1=yplt[j];
            y3=yplt[j-1];

            ynode[0]=yy1;
            ynode[1]=y3;
            ynode[2]=y3;
            ynode[3]=yy1;

            n11 = iterm+jterm+kterm;
            n12 = n11 - nx + 1;
            n22 = n12 + nxy;
            n21 = n11 + nxy;

//        n11 = (j-js1)*nx   + (i-is1)   + (k-ks1)*nx*ny;
//        n12 = (j-1-js1)*nx + (i+1-is1) + (k-ks1)*nx*ny;
//        n22 = (j-1-js1)*nx + (i+1-is1) + (k+1-ks1)*nx*ny;
//        n21 = (j-js1)*nx   + (i-is1)   + (k+1-ks1)*nx*ny;

            if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
              int m, m11, m12, m22, m21;

              m = iterm+jterm;
              m11 = m;
              m12 = m11 - nx + 1;
              m22 = m12;
              m21 = m11;
              DRAWVERTEXGPUTERRAIN(xnode[mm],ynode[mm],znode[mm])
            }
            else{
              DRAWVERTEXGPU(xnode[mm],ynode[mm],znode[mm])
            }
          }
        }
        break;

  // +++++++++++++++++++++++++++++++++++ DIR 5 +++++++++++++++++++++++++++++++++++++++

      case 5:
      case -5:
        for(k=culli->kbeg; k<culli->kend; k++){
          kterm = (k-ks1)*nxy;
          z1 = zplt[k];
          z3 = zplt[k+1];
          znode[0]=z1;
          znode[1]=z1;
          znode[2]=z3;
          znode[3]=z3;

          for(i=culli->ibeg;i<culli->iend;i++){
            iterm = (i-is1);
            x1 = xplt[i];
            x3 = xplt[i+1];

            xnode[0]=x1;
            xnode[1]=x3;
            xnode[2]=x3;
            xnode[3]=x1;

            j = culli->jbeg+(i-culli->ibeg);
            jterm = (j-js1)*nx;

            yy1=yplt[j];
            y3=yplt[j+1];

            ynode[0]=yy1;
            ynode[1]=y3;
            ynode[2]=y3;
            ynode[3]=yy1;

            n11 = jterm + iterm + kterm;
            n12 = n11 + nx + 1;
            n22 = n12 + nxy;
            n21 = n11 + nxy;

        //    n11 = (j-js1)*nx + (i-is1) + (k-ks1)*nx*ny;
        //    n12 = (j+1-js1)*nx + (i+1-is1) + (k-ks1)*nx*ny;
        //    n22 = (j+1-js1)*nx + (i+1-is1) + (k+1-ks1)*nx*ny;
        //    n21 = (j-js1)*nx + (i-is1) + (k+1-ks1)*nx*ny;

            if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
              int m, m11, m12, m22, m21;

              m = iterm+jterm;
              m11 = m;
              m12 = m11 + nx + 1;
              m22 = m12;
              m21 = m11;
              DRAWVERTEXGPUTERRAIN(xnode[mm],ynode[mm],znode[mm])
            }
            else{
              DRAWVERTEXGPU(xnode[mm],ynode[mm],znode[mm])
            }
          }
        }
        break;

  // +++++++++++++++++++++++++++++++++++ DIR 6 +++++++++++++++++++++++++++++++++++++++

      case 6:
      case -6:
        for(i=culli->ibeg; i<culli->iend; i++){
          iterm = (i-is1);
          x1 = xplt[i];
          x3 = xplt[i+1];

          xnode[0]=x1;
          xnode[1]=x3;
          xnode[2]=x3;
          xnode[3]=x1;

          for(j=culli->jend;j>culli->jbeg;j--){

            k = culli->kbeg-(j-culli->jend);
            jterm = (j-js1)*nx;

            yy1=yplt[j];
            y3=yplt[j-1];

            ynode[0]=yy1;
            ynode[1]=yy1;
            ynode[2]=y3;
            ynode[3]=y3;

            kterm = (k-ks1)*nxy;
            z1 = zplt[k];
            z3 = zplt[k+1];
            znode[0]=z1;
            znode[1]=z1;
            znode[2]=z3;
            znode[3]=z3;


            n11 = jterm + iterm + kterm;
            n12 = n11 + 1;
            n22 = n12 -nx + nxy;
            n21 = n11 -nx +  nxy;

        //    n11 = (j-js1)*nx + (i-is1) + (k-ks1)*nx*ny;
        //    n12 = (j-js1)*nx + (i+1-is1) + (k-ks1)*nx*ny;
        //    n22 = (j-1-js1)*nx + (i+1-is1) + (k+1-ks1)*nx*ny;
        //    n21 = (j-1-js1)*nx + (i-is1) + (k+1-ks1)*nx*ny;

            if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
              int m, m11, m12, m22, m21;

              m = iterm+jterm;
              m11 = m;
              m12 = m11 + 1;
              m22 = m12 - nx;
              m21 = m11 - nx;
              DRAWVERTEXGPUTERRAIN(xnode[mm],ynode[mm],znode[mm])
            }
            else{
              DRAWVERTEXGPU(xnode[mm],ynode[mm],znode[mm])
            }
          }
        }
        break;

  // +++++++++++++++++++++++++++++++++++ DIR 7 +++++++++++++++++++++++++++++++++++++++

      case 7:
      case -7:
        for(i=culli->ibeg; i<culli->iend; i++){
          iterm = (i-is1);
          x1 = xplt[i];
          x3 = xplt[i+1];

          xnode[0]=x1;
          xnode[1]=x3;
          xnode[2]=x3;
          xnode[3]=x1;

          for(j=culli->jbeg;j<culli->jend;j++){

            k = culli->kbeg+(j-culli->jbeg);
            jterm = (j-js1)*nx;

            yy1=yplt[j];
            y3=yplt[j+1];

            ynode[0]=yy1;
            ynode[1]=yy1;
            ynode[2]=y3;
            ynode[3]=y3;

            kterm = (k-ks1)*nxy;
            z1 = zplt[k];
            z3 = zplt[k+1];
            znode[0]=z1;
            znode[1]=z1;
            znode[2]=z3;
            znode[3]=z3;


            n11 = jterm + iterm + kterm;
            n12 = n11 + 1;
            n22 = n12 +nx + nxy;
            n21 = n11 +nx +  nxy;

        //    n11 = (j-js1)*nx + (i-is1) + (k-ks1)*nx*ny;
        //    n12 = (j-js1)*nx + (i+1-is1) + (k-ks1)*nx*ny;
        //    n22 = (j+1-js1)*nx + (i+1-is1) + (k+1-ks1)*nx*ny;
        //    n21 = (j+1-js1)*nx + (i-is1) + (k+1-ks1)*nx*ny;

            if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
              int m, m11, m12, m22, m21;

              m = iterm+jterm;
              m11 = m;
              m12 = m11 + 1;
              m22 = m12 + nx;
              m21 = m11 + nx;
              DRAWVERTEXGPUTERRAIN(xnode[mm],ynode[mm],znode[mm])
            }
            else{
              DRAWVERTEXGPU(xnode[mm],ynode[mm],znode[mm])
            }
          }
        }
        break;

  // +++++++++++++++++++++++++++++++++++ DIR 8 +++++++++++++++++++++++++++++++++++++++

      case 8:
      case -8:
        for(j=culli->jbeg; j<culli->jend; j++){
          jterm = (j-js1)*nx;

          yy1=yplt[j];
          y3=yplt[j+1];

          ynode[0]=yy1;
          ynode[1]=yy1;
          ynode[2]=y3;
          ynode[3]=y3;

          for(i=culli->ibeg;i<culli->iend;i++){
            iterm = (i-is1);
            x1 = xplt[i];
            x3 = xplt[i+1];

            xnode[0]=x1;
            xnode[1]=x3;
            xnode[2]=x3;
            xnode[3]=x1;

            k = culli->kend-(i-culli->ibeg);
            kterm = (k-ks1)*nxy;
            z1 = zplt[k];
            z3 = zplt[k-1];
            znode[0]=z1;
            znode[1]=z3;
            znode[2]=z3;
            znode[3]=z1;

            n11 = jterm + iterm + kterm;
            n12 = n11 + 1-nx;
            n22 = n12 + nx;
            n21 = n11 +nx ;

        //    n11 = (j-js1)*nx + (i-is1) + (k-ks1)*nx*ny;
        //    n12 = (j-js1)*nx + (i+1-is1) + (k-1-ks1)*nx*ny;
        //    n22 = (j+1-js1)*nx + (i+1-is1) + (k-1-ks1)*nx*ny;
        //    n21 = (j+1-js1)*nx + (i-is1) + (k-ks1)*nx*ny;

            if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
              int m, m11, m12, m22, m21;

              m = iterm+jterm;
              m11 = m;
              m12 = m11 + 1 - nx;
              m22 = m12 + nx;
              m21 = m11 + nx;
              DRAWVERTEXGPUTERRAIN(xnode[mm],ynode[mm],znode[mm])
            }
            else{
              DRAWVERTEXGPU(xnode[mm],ynode[mm],znode[mm])
            }
          }
        }
        break;

  // +++++++++++++++++++++++++++++++++++ DIR 9 +++++++++++++++++++++++++++++++++++++++

      case 9:
      case -9:
        for(j=culli->jbeg; j<culli->jend; j++){
          jterm = (j-js1)*nx;

          yy1=yplt[j];
          y3=yplt[j+1];

          ynode[0]=yy1;
          ynode[1]=y3;
          ynode[2]=y3;
          ynode[3]=yy1;

          for(i=culli->ibeg;i<culli->iend;i++){
            iterm = (i-is1);
            x1 = xplt[i];
            x3 = xplt[i+1];

            xnode[0]=x1;
            xnode[1]=x1;
            xnode[2]=x3;
            xnode[3]=x3;

            k = culli->kbeg+(i-culli->ibeg);
            kterm = (k-ks1)*nxy;
            z1 = zplt[k];
            z3 = zplt[k+1];
            znode[0]=z1;
            znode[1]=z1;
            znode[2]=z3;
            znode[3]=z3;

            n11 = jterm + iterm + kterm;
            n12 = n11 + nx;
            n22 = n12 + 1 + nxy;
            n21 = n11 + 1 + nxy ;

        //    n11 = (j-js1)*nx + (i-is1) + (k-ks1)*nx*ny;
        //    n12 = (j+1-js1)*nx + (i-is1) + (k-ks1)*nx*ny;
        //    n22 = (j+1-js1)*nx + (i+1-is1) + (k+1-ks1)*nx*ny;
        //    n21 = (j-js1)*nx + (i+1-is1) + (k+1-ks1)*nx*ny;

            if(nterraininfo>0&&ABS(vertical_factor-1.0)>0.01){
              int m, m11, m12, m22, m21;

              m = iterm+jterm;
              m11 = m;
              m12 = m11 + nx;
              m22 = m12 + 1;
              m21 = m11 + 1;
              DRAWVERTEXGPUTERRAIN(xnode[mm],ynode[mm],znode[mm])
            }
            else{
              DRAWVERTEXGPU(xnode[mm],ynode[mm],znode[mm])
            }
          }
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
  }
  glEnd();
  SNIFF_ERRORS("in drawsmoke3dcull 12");
  TransparentOff();
  if(cullfaces==1)glEnable(GL_CULL_FACE);
  SNIFF_ERRORS("after drawsmokecull");
}

#endif

/* ------------------ update_smoke3d_menulabels ------------------------ */

void update_smoke3d_menulabels(void){
  int i;
  smoke3ddata *smoke3di;
  char meshlabel[128];

  PRINTF("  updating smoke3d menu labels");
  for(i=0;i<nsmoke3dinfo;i++){
    smoke3di = smoke3dinfo + i;
    STRCPY(smoke3di->menulabel,smoke3di->label.longlabel);
    if(showfiles==1){
      STRCAT(smoke3di->menulabel,", ");
      STRCAT(smoke3di->menulabel,smoke3di->file);
    }
    switch(smoke3di->compression_type){
    case UNKNOWN:
      // compression type not determined yet
      break;
    case RLE:
      STRCAT(smoke3di->menulabel," (RLE) ");
      break;
    case ZLIB:
      STRCAT(smoke3di->menulabel," (ZLIB) ");
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
//    len=strlen(smoke3di->menulabel);
    if(nmeshes>1){
	  meshdata *smokemesh;

	  smokemesh = meshinfo + smoke3di->blocknumber;
      sprintf(meshlabel,"%s",smokemesh->label);
      STRCAT(smoke3di->menulabel," - ");
      STRCAT(smoke3di->menulabel,meshlabel);
    }
  }
  PRINTF(" - completed\n");
}

/* ------------------ adjustalpha ------------------------ */

unsigned char adjustalpha(unsigned char alpha, float factor){
  double val, rr;
  float falpha;
  float term1, term2, term3, term4;

  rr = factor;
  falpha = alpha/255.0;

  //val = 1.0 - pow(1.0-falpha,rr);
  term1 = falpha*rr;
  term2 = falpha*(rr-1.0)/2.0;
  term3 = falpha*(rr-2.0)/3.0;
  term4 = falpha*(rr-3.0)/4.0;
  val = term1*(1.0 - term2*(1.0 -  term3*(1.0 - term4)));

  val = 255*val+0.5;
  if(val>255)val=255;
  alpha=val;
  return alpha;
}

/* ------------------ makeiblanksmoke3d ------------------------ */

void makeiblank_smoke3d(void){
  int i, ii;
  int ic;

  PRINTF("  initializing 3D smoke/hrrpuv blanking array - ");
  update_makeiblank_smoke3d=0;
  for(i=0;i<nsmoke3dinfo;i++){
    smoke3ddata *smoke3di;
    meshdata *smokemesh;
    int ibar, jbar, kbar;
    int ijksize;

    smoke3di = smoke3dinfo + i;
    smokemesh = meshinfo + smoke3di->blocknumber;

    ibar = smokemesh->ibar;
    jbar = smokemesh->jbar;
    kbar = smokemesh->kbar;
    ijksize=(ibar+1)*(jbar+1)*(kbar+1);

    if(use_iblank==1&&smoke3di->loaded==1&&smokemesh->iblank_smoke3d==NULL){
      unsigned char *iblank_smoke3d;

      NewMemory((void **)&iblank_smoke3d, ijksize*sizeof(unsigned char));
      smokemesh->iblank_smoke3d=iblank_smoke3d;
    }
  }

#define ALLMESHES 0
#define LOWERMESHES 1

  for(ic=nmeshes-1;ic>=0;ic--){
    meshdata *smokemesh;
    unsigned char *iblank_smoke3d;
    float *xplt, *yplt, *zplt;
    float dx, dy, dz;
    int nx, ny, nxy;
    int ibar, jbar, kbar;
    int ijksize;
    int j, k;

    smokemesh = meshinfo + ic;
    iblank_smoke3d = smokemesh->iblank_smoke3d;
   // if(iblank_smoke3d==NULL||smokemesh->iblank_smoke3d_defined==1)continue;
    if(iblank_smoke3d==NULL)continue;
    smokemesh->iblank_smoke3d_defined = 1;

    xplt=smokemesh->xplt;
    yplt=smokemesh->yplt;
    zplt=smokemesh->zplt;
    dx = xplt[1]-xplt[0];
    dy = yplt[1]-yplt[0];
    dz = zplt[1]-zplt[0];

    ibar = smokemesh->ibar;
    jbar = smokemesh->jbar;
    kbar = smokemesh->kbar;
    ijksize=(ibar+1)*(jbar+1)*(kbar+1);
    nx = ibar + 1;
    ny = jbar + 1;
    nxy = nx*ny;

    for(ii=0;ii<ijksize;ii++){
      *iblank_smoke3d++=GAS;
    }

    iblank_smoke3d=smokemesh->iblank_smoke3d;

    for(i=0;i<=smokemesh->ibar;i++){
      for(j=0;j<=smokemesh->jbar;j++){
        for(k=0;k<=smokemesh->kbar;k++){
          float x, y, z;
          int ijk;

          ijk = IJKNODE(i, j, k);
          x = xplt[i];
          y = yplt[j];
          z = zplt[k];
          if(inmesh_smoke(x,y,z,ic-1,LOWERMESHES)>=0)iblank_smoke3d[ijk]=SOLID;
        }
      }
    }

    for(ii=0;ii<smokemesh->nbptrs;ii++){
      blockagedata *bc;

      bc = smokemesh->blockageinfoptrs[ii];
      if(bc->invisible==1||bc->hidden==1||bc->nshowtime!=0)continue;
      for(i=bc->ijk[IMIN];i<=bc->ijk[IMAX];i++){
        for(j=bc->ijk[JMIN];j<=bc->ijk[JMAX];j++){
          for(k=bc->ijk[KMIN];k<=bc->ijk[KMAX];k++){
            int ijk;

            ijk = IJKNODE(i, j, k);
            iblank_smoke3d[ijk]=SOLID;
          }
        }
      }
    }

    for(j=0;j<=jbar;j++){
      for(k=0;k<=kbar;k++){
        float x, y, z;
        int ijk;

        ijk = IJKNODE(0, j, k);
        x = xplt[0];
        y = yplt[j];
        z = zplt[k];
        if(inmesh_smoke(x-dx,y,z,ic,ALLMESHES)<0){
          iblank_smoke3d[ijk]=SOLID;
        }
        else{
          iblank_smoke3d[ijk]=GAS;
        }

        ijk = IJKNODE(ibar,j,k);
        x = xplt[ibar];
        y = yplt[j];
        z = zplt[k];
        if(inmesh_smoke(x+dx,y,z,ic,ALLMESHES)<0){
          iblank_smoke3d[ijk]=SOLID;
        }
        else{
          iblank_smoke3d[ijk]=GAS;
        }
      }
    }

    for(i=0;i<=ibar;i++){
      for(k=0;k<=kbar;k++){
        float x, y, z;
        int ijk;

        ijk = IJKNODE(i, 0, k);
        x = xplt[i];
        y = yplt[0];
        z = zplt[k];
        if(inmesh_smoke(x,y-dy,z,ic,ALLMESHES)<0){
          iblank_smoke3d[ijk]=SOLID;
        }
        else{
          iblank_smoke3d[ijk]=GAS;
        }

        ijk = IJKNODE(i,jbar,k);
        x = xplt[i];
        y = yplt[jbar];
        z = zplt[k];
        if(inmesh_smoke(x,y+dy,z,ic,ALLMESHES)<0){
          iblank_smoke3d[ijk]=SOLID;
        }
        else{
          iblank_smoke3d[ijk]=GAS;
        }
      }
    }

    for(i=0;i<=ibar;i++){
      for(j=0;j<=jbar;j++){
        float x, y, z;
        int ijk;

        ijk = IJKNODE(i, j, 0);
        x = xplt[i];
        y = yplt[j];
        z = zplt[0];
        if(inmesh_smoke(x,y,z-dz,ic,ALLMESHES)<0){
          iblank_smoke3d[ijk]=SOLID;
        }
        else{
          iblank_smoke3d[ijk]=GAS;
        }

        ijk = IJKNODE(i,j,kbar);
        x = xplt[i];
        y = yplt[j];
        z = zplt[kbar];
        if(inmesh_smoke(x,y,z+dz,ic,ALLMESHES)<0){
          iblank_smoke3d[ijk]=SOLID;
        }
        else{
          iblank_smoke3d[ijk]=GAS;
        }
      }
    }
  }
  PRINTF("complete\n");
}

/* ------------------ inmesh_smoke ------------------------ */

 int inmesh_smoke(float x, float y, float z, int nm, int flag){
  int i;
  int n;

  if(flag==LOWERMESHES){
    n = nm;
  }
  else{
    n = nmeshes;
  }
  for(i=0;i<n;i++){
    meshdata *meshi;

    meshi = meshinfo + i;
    if(flag==ALLMESHES&&i==nm)continue;
    if(meshi->iblank_smoke3d == NULL)continue;

    if(x<meshi->xplt[0]||x>meshi->xplt[meshi->ibar])continue;
    if(y<meshi->yplt[0]||y>meshi->yplt[meshi->jbar])continue;
    if(z<meshi->zplt[0]||z>meshi->zplt[meshi->kbar])continue;
    return i;
  }
  return -1;
}

/* ------------------ getsmoke3d_version ------------------------ */

 int getsmoke3d_version(smoke3ddata *smoke3di){

   FILE *SMOKE3DFILE=NULL,*SMOKE3D_REGFILE=NULL, *SMOKE3D_COMPFILE=NULL;
   int nxyz[8];
   char *file;
   int fortran_skip=0;

   if(smoke3di->filetype==FORTRAN_GENERATED)fortran_skip=4;

   smoke3di->have_light=0;
   file=smoke3di->comp_file;
   SMOKE3D_COMPFILE=fopen(file,"rb");
   if(SMOKE3D_COMPFILE==NULL){
     file=smoke3di->reg_file;
     SMOKE3D_REGFILE=fopen(file,"rb");
   }
   if(SMOKE3D_REGFILE==NULL&&SMOKE3D_COMPFILE==NULL)return -1;
   if(SMOKE3D_COMPFILE!=NULL)SMOKE3DFILE=SMOKE3D_COMPFILE;
   if(SMOKE3D_REGFILE!=NULL)SMOKE3DFILE=SMOKE3D_REGFILE;
   smoke3di->file=file;

   SKIP;fread(nxyz,4,8,SMOKE3DFILE);SKIP;
   {
     float time_local;
     int nchars[2];

     SKIP;fread(&time_local,4,1,SMOKE3DFILE);SKIP;
     if(feof(SMOKE3DFILE)==0){
       SKIP;fread(nchars,4,2,SMOKE3DFILE);SKIP;
       if(feof(SMOKE3DFILE)==0){
         if(nchars[1]<0){
           smoke3di->have_light=1;
         }
       }
     }
   }

   fclose(SMOKE3DFILE);


   return nxyz[1];
 }

#ifdef pp_CULL

    /* ------------------ initcullplane ------------------------ */

void initcullplane(int cullflag){
  int iii;
  int iskip, jskip, kskip;
  int nx, ny, nz;
  int nxx, nyy, nzz;
  int ibeg, iend, jbeg, jend, kbeg, kend;
  int i, j, k, ij;
  cullplanedata *cp;
  smoke3ddata *smoke3di;
  float norm[3];
  float dx, dy, dz, factor;

#ifdef _DEBUG
  PRINTF("updating initcullplane %i\n",cull_count++);
#endif
  update_initcullplane=0;
  cp = cullplaneinfo;
  ncullplaneinfo=0;

  for(iii=0;iii<nmeshes;iii++){
    meshdata *meshi;
    culldata *culli;

    meshi=meshinfo+iii;

    meshi->culldefined=0;
    if(meshi->cull_smoke3d==NULL)continue;
    smoke3di = meshi->cull_smoke3d;
    if(smoke3di->loaded==0||smoke3di->display==0)continue;
    meshi->culldefined=1;

    culli=meshi->cullinfo;
    GetCullSkips(meshi,cullflag, cull_portsize, &iskip,&jskip,&kskip);

    nx = meshi->ibar/iskip + 1;
    ny = meshi->jbar/jskip + 1;
    nz = meshi->kbar/kskip + 1;

    for(k=0;k<nz;k++){
      kbeg = k*kskip;
      kend = kbeg + kskip;
      if(kend>meshi->kbar)kend=meshi->kbar;
      for(j=0;j<ny;j++){
        jbeg = j*jskip;
        jend = jbeg + jskip;
        if(jend>meshi->jbar)jend=meshi->jbar;
        for(i=0;i<nx;i++){
          ibeg = i*iskip;
          iend = ibeg + iskip;
          if(iend>meshi->ibar)iend=meshi->ibar;

          if(culli->npixels>0||cullflag==0){
            cp->norm[0]=0.0;
            cp->norm[1]=0.0;
            cp->norm[2]=0.0;
            switch(meshi->smokedir){
              int ii, jj, kk;

              case 1:
              case -1:
                for(ii=ibeg;ii<iend;ii++){
                  cp->dir=meshi->smokedir;
                  cp->cull=culli;
                  cp->cull_mesh=meshi;

                  cp->ibeg=ii;
                  cp->iend=ii;
                  cp->jbeg=jbeg;
                  cp->jend=jend;
                  cp->kbeg=kbeg;
                  cp->kend=kend;
                  if(meshi->smokedir>0){
                    cp->norm[0]=1.0;
                  }
                  else{
                    cp->norm[0]=-1.0;
                  }
                  cp->xmin=meshi->xplt[cp->ibeg];
                  cp->xmax=meshi->xplt[cp->iend];
                  cp->ymin=meshi->yplt[cp->jbeg];
                  cp->ymax=meshi->yplt[cp->jend];
                  cp->zmin=meshi->zplt[cp->kbeg];
                  cp->zmax=meshi->zplt[cp->kend];

                  cp++;
                  ncullplaneinfo++;
                }
                break;
              case 2:
              case -2:
                for(jj=jbeg;jj<jend;jj++){
                  cp->dir=meshi->smokedir;
                  cp->cull=culli;
                  cp->cull_mesh=meshi;

                  cp->ibeg=ibeg;
                  cp->iend=iend;
                  cp->jbeg=jj;
                  cp->jend=jj;
                  cp->kbeg=kbeg;
                  cp->kend=kend;
                  if(meshi->smokedir>0){
                    cp->norm[1]=1.0;
                  }
                  else{
                    cp->norm[1]=-1.0;
                  }

                  cp->xmin=meshi->xplt[cp->ibeg];
                  cp->xmax=meshi->xplt[cp->iend];
                  cp->ymin=meshi->yplt[cp->jbeg];
                  cp->ymax=meshi->yplt[cp->jend];
                  cp->zmin=meshi->zplt[cp->kbeg];
                  cp->zmax=meshi->zplt[cp->kend];

                  ncullplaneinfo++;
                  cp++;
                }
                break;
              case 3:
              case -3:
                for(kk=kbeg;kk<kend;kk++){
                  cp->dir=meshi->smokedir;
                  cp->cull=culli;
                  cp->cull_mesh=meshi;

                  cp->ibeg=ibeg;
                  cp->iend=iend;
                  cp->jbeg=jbeg;
                  cp->jend=jend;
                  cp->kbeg=kk;
                  cp->kend=kk;
                  if(meshi->smokedir>0){
                    cp->norm[2]=1.0;
                  }
                  else{
                    cp->norm[2]=-1.0;
                  }

                  cp->xmin=meshi->xplt[cp->ibeg];
                  cp->xmax=meshi->xplt[cp->iend];
                  cp->ymin=meshi->yplt[cp->jbeg];
                  cp->ymax=meshi->yplt[cp->jend];
                  cp->zmin=meshi->zplt[cp->kbeg];
                  cp->zmax=meshi->zplt[cp->kend];

                  ncullplaneinfo++;
                  cp++;
                }
                break;
              case 4:
              case -4:
                nxx = iend - ibeg;
                nyy = jend - jbeg;
                norm[2]=0.0;
                dx = meshi->xplt_orig[1]-meshi->xplt_orig[0];
                dy = meshi->yplt_orig[1]-meshi->yplt_orig[0];
                factor= dx*dx+dy*dy;
                if(factor==0.0){
                  factor=1.0;
                }
                else{
                  factor=1.0/sqrt(factor);
                }
                norm[0]=dy*factor;
                norm[1]=dx*factor;
                for(ij=1;ij<nxx+nyy+1;ij++){
                  cp->dir=meshi->smokedir;
                  cp->cull=culli;
                  cp->cull_mesh=meshi;

                  if(ij<=nyy){
                    cp->ibeg=ibeg;
                  }
                  else{
                    cp->ibeg=ibeg+(ij-nyy);
                  }
                  cp->iend = ibeg + ij;
                  if(cp->iend>iend)cp->iend=iend;

                  if(ij<=nxx){
                    cp->jbeg=jbeg;
                  }
                  else{
                    cp->jbeg=jbeg+(ij-nxx);
                  }
                  cp->jend = jbeg + ij;
                  if(cp->jend>jend)cp->jend=jend;

                  cp->kbeg=kbeg;
                  cp->kend=kend;

                  if(meshi->smokedir>0){
                    cp->norm[0]=norm[0];
                    cp->norm[1]=norm[1];
                    cp->norm[2]=norm[2];
                  }
                  else{
                    cp->norm[0]=-norm[0];
                    cp->norm[1]=-norm[1];
                    cp->norm[2]=-norm[2];
                  }

                  cp->xmin=meshi->xplt[cp->ibeg];
                  cp->xmax=meshi->xplt[cp->iend];
                  cp->ymin=meshi->yplt[cp->jbeg];
                  cp->ymax=meshi->yplt[cp->jend];
                  cp->zmin=meshi->zplt[cp->kbeg];
                  cp->zmax=meshi->zplt[cp->kend];

                  ncullplaneinfo++;
                  cp++;
                }
                break;
              case 5:
              case -5:
                nxx = iend - ibeg;
                nyy = jend - jbeg;
                norm[2]=0.0;
                dx = meshi->xplt_orig[1]-meshi->xplt_orig[0];
                dy = meshi->yplt_orig[1]-meshi->yplt_orig[0];
                factor= dx*dx+dy*dy;
                if(factor==0.0){
                  factor=1.0;
                }
                else{
                  factor=1.0/sqrt(factor);
                }
                norm[0]=-dy*factor;
                norm[1]=dx*factor;
                for(ij=1;ij<nxx+nyy+1;ij++){
                  cp->dir=meshi->smokedir;
                  cp->cull=culli;
                  cp->cull_mesh=meshi;

                  cp->ibeg=iend-ij;
                  if(ij>nxx){
                    cp->ibeg=ibeg;
                  }
                  cp->iend = iend;
                  if(ij>nyy)cp->iend=iend-(ij-nyy);

                  if(ij<=nxx){
                    cp->jbeg=jbeg;
                  }
                  else{
                    cp->jbeg=jbeg+(ij-nxx);
                  }
                  cp->jend = jbeg + ij;
                  if(cp->jend>jend)cp->jend=jend;

                  cp->kbeg=kbeg;
                  cp->kend=kend;

                  if(meshi->smokedir>0){
                    cp->norm[0]=norm[0];
                    cp->norm[1]=norm[1];
                    cp->norm[2]=norm[2];
                  }
                  else{
                    cp->norm[0]=-norm[0];
                    cp->norm[1]=-norm[1];
                    cp->norm[2]=-norm[2];
                  }

                  cp->xmin=meshi->xplt[cp->ibeg];
                  cp->xmax=meshi->xplt[cp->iend];
                  cp->ymin=meshi->yplt[cp->jbeg];
                  cp->ymax=meshi->yplt[cp->jend];
                  cp->zmin=meshi->zplt[cp->kbeg];
                  cp->zmax=meshi->zplt[cp->kend];

                  ncullplaneinfo++;
                  cp++;
                }
                break;
              case 6:
              case -6:
                nyy = jend - jbeg;
                nzz = kend - kbeg;
                norm[0]=0.0;
                dy = meshi->yplt_orig[1]-meshi->yplt_orig[0];
                dz = meshi->zplt_orig[1]-meshi->zplt_orig[0];
                factor= dz*dz+dy*dy;
                if(factor==0.0){
                  factor=1.0;
                }
                else{
                  factor=1.0/sqrt(factor);
                }
                norm[1]=dz*factor;
                norm[2]=dy*factor;
                for(ij=1;ij<nzz+nyy+1;ij++){
                  cp->dir=meshi->smokedir;
                  cp->cull=culli;
                  cp->cull_mesh=meshi;

                  cp->ibeg=ibeg;
                  cp->iend=iend;

                  cp->jbeg=jbeg;
                  if(ij>nzz){
                    cp->jbeg=jbeg+(ij-nzz);
                  }
                  cp->jend=jbeg+ij;
                  if(ij>nyy){
                    cp->jend=jend;
                  }

                  cp->kbeg = kbeg;
                  if(ij>nyy){
                    cp->kbeg=kbeg+(ij-nyy);
                  }
                  cp->kend = kbeg+ij;
                  if(ij>nzz){
                    cp->kend=kend;
                  }

                  if(meshi->smokedir>0){
                    cp->norm[0]=norm[0];
                    cp->norm[1]=norm[1];
                    cp->norm[2]=norm[2];
                  }
                  else{
                    cp->norm[0]=-norm[0];
                    cp->norm[1]=-norm[1];
                    cp->norm[2]=-norm[2];
                  }

                  cp->xmin=meshi->xplt[cp->ibeg];
                  cp->xmax=meshi->xplt[cp->iend];
                  cp->ymin=meshi->yplt[cp->jbeg];
                  cp->ymax=meshi->yplt[cp->jend];
                  cp->zmin=meshi->zplt[cp->kbeg];
                  cp->zmax=meshi->zplt[cp->kend];

                  ncullplaneinfo++;
                  cp++;
                }
                break;
              case 7:
              case -7:
                nyy = jend - jbeg;
                nzz = kend - kbeg;
                norm[0]=0.0;
                dy = meshi->yplt_orig[1]-meshi->yplt_orig[0];
                dz = meshi->zplt_orig[1]-meshi->zplt_orig[0];
                factor= dz*dz+dy*dy;
                if(factor==0.0){
                  factor=1.0;
                }
                else{
                  factor=1.0/sqrt(factor);
                }
                norm[1]=-dz*factor;
                norm[2]=dy*factor;
                for(ij=1;ij<nzz+nyy+1;ij++){
                  cp->dir=meshi->smokedir;
                  cp->cull=culli;
                  cp->cull_mesh=meshi;

                  cp->ibeg=ibeg;
                  cp->iend=iend;

                  cp->jbeg=jbeg;
                  if(ij>nzz){
                    cp->jbeg=jbeg+(ij-nzz);
                  }
                  cp->jend=jbeg+ij;
                  if(ij>nyy){
                    cp->jend=jend;
                  }

                  cp->kbeg = kend-ij;
                  if(ij>nzz){
                    cp->kbeg=kbeg;
                  }
                  cp->kend = kend;
                  if(ij>nyy){
                    cp->kend=kend-(ij-nyy);
                  }

                  if(meshi->smokedir>0){
                    cp->norm[0]=norm[0];
                    cp->norm[1]=norm[1];
                    cp->norm[2]=norm[2];
                  }
                  else{
                    cp->norm[0]=-norm[0];
                    cp->norm[1]=-norm[1];
                    cp->norm[2]=-norm[2];
                  }

                  cp->xmin=meshi->xplt[cp->ibeg];
                  cp->xmax=meshi->xplt[cp->iend];
                  cp->ymin=meshi->yplt[cp->jbeg];
                  cp->ymax=meshi->yplt[cp->jend];
                  cp->zmin=meshi->zplt[cp->kbeg];
                  cp->zmax=meshi->zplt[cp->kend];

                  ncullplaneinfo++;
                  cp++;
                }
                break;
              case 8:
              case -8:
                nxx = iend - ibeg;
                nzz = kend - kbeg;
                norm[1]=0.0;
                dx = meshi->xplt_orig[1]-meshi->xplt_orig[0];
                dz = meshi->zplt_orig[1]-meshi->zplt_orig[0];
                factor= dz*dz+dx*dx;
                if(factor==0.0){
                  factor=1.0;
                }
                else{
                  factor=1.0/sqrt(factor);
                }
                norm[0]=dz*factor;
                norm[2]=dx*factor;
                for(ij=1;ij<nzz+nxx+1;ij++){
                  cp->dir=meshi->smokedir;
                  cp->cull=culli;
                  cp->cull_mesh=meshi;


                  cp->ibeg=iend-ij;
                  if(ij>nxx){
                    cp->ibeg=ibeg;
                  }
                  ASSERT(cp->ibeg>=ibeg);
                  cp->iend=iend;
                  if(ij>nzz){
                    cp->iend=iend-(ij-nzz);
                  }

                  cp->jbeg=jbeg;
                  cp->jend=jend;

                  cp->kbeg = kbeg-ij;
                  if(ij>nzz){
                    cp->kbeg=kbeg;
                  }
                  cp->kend = kend;
                  if(ij>nxx){
                    cp->kend=kend-(ij-nxx);
                  }

                  if(meshi->smokedir>0){
                    cp->norm[0]=norm[0];
                    cp->norm[1]=norm[1];
                    cp->norm[2]=norm[2];
                  }
                  else{
                    cp->norm[0]=-norm[0];
                    cp->norm[1]=-norm[1];
                    cp->norm[2]=-norm[2];
                  }

                  cp->xmin=meshi->xplt[cp->ibeg];
                  cp->xmax=meshi->xplt[cp->iend];
                  cp->ymin=meshi->yplt[cp->jbeg];
                  cp->ymax=meshi->yplt[cp->jend];
                  cp->zmin=meshi->zplt[cp->kbeg];
                  cp->zmax=meshi->zplt[cp->kend];

                  ncullplaneinfo++;
                  cp++;
                }
                break;
              case 9:
              case -9:
                nxx = iend - ibeg;
                nzz = kend - kbeg;
                norm[1]=0.0;
                dx = meshi->xplt_orig[1]-meshi->xplt_orig[0];
                dz = meshi->zplt_orig[1]-meshi->zplt_orig[0];
                factor= dz*dz+dx*dx;
                if(factor==0.0){
                  factor=1.0;
                }
                else{
                  factor=1.0/sqrt(factor);
                }
                norm[0]=-dz*factor;
                norm[2]=dx*factor;
                for(ij=1;ij<nzz+nxx+1;ij++){
                  cp->dir=meshi->smokedir;
                  cp->cull=culli;
                  cp->cull_mesh=meshi;


                  cp->ibeg=iend-ij;
                  if(ij>nxx){
                    cp->ibeg=ibeg;
                  }
                  ASSERT(cp->ibeg>=ibeg);
                  cp->iend=iend;
                  if(ij>nzz){
                    cp->iend=iend-(ij-nzz);
                  }

                  cp->jbeg=jbeg;
                  cp->jend=jend;

                  cp->kbeg = kbeg;
                  if(ij>nxx){
                    cp->kbeg=kbeg+(ij-nxx);
                  }
                  cp->kend = kbeg+ij;
                  if(ij>nzz){
                    cp->kend=kend;
                  }

                  if(meshi->smokedir>0){
                    cp->norm[0]=norm[0];
                    cp->norm[1]=norm[1];
                    cp->norm[2]=norm[2];
                  }
                  else{
                    cp->norm[0]=-norm[0];
                    cp->norm[1]=-norm[1];
                    cp->norm[2]=-norm[2];
                  }

                  cp->xmin=meshi->xplt[cp->ibeg];
                  cp->xmax=meshi->xplt[cp->iend];
                  cp->ymin=meshi->yplt[cp->jbeg];
                  cp->ymax=meshi->yplt[cp->jend];
                  cp->zmin=meshi->zplt[cp->kbeg];
                  cp->zmax=meshi->zplt[cp->kend];

                  ncullplaneinfo++;
                  cp++;
                }
                break;
              default:
                ASSERT(FFALSE);
                break;
            }
          }
          culli++;
        }
      }
    }
  }

  for(iii=0;iii<ncullplaneinfo;iii++){
    sort_cullplaneinfo[iii]=cullplaneinfo+iii;
  }
  if(ncullplaneinfo>1){
    qsort((cullplanedata *)sort_cullplaneinfo,(size_t)ncullplaneinfo,
       sizeof(cullplanedata *),cullplane_compare);
  }
  glutPostRedisplay();
}

/* ------------------ cullplane_compare ------------------------ */

int cullplane_compare( const void *arg1, const void *arg2 ){
  cullplanedata *cpi, *cpj;

  cpi = *(cullplanedata **)arg1;
  cpj = *(cullplanedata **)arg2;

  if(cpi->dir==cpj->dir){
    float di, dj;

    di = cpi->norm[0]*cpi->xmin + cpi->norm[1]*cpi->ymin + cpi->norm[2]*cpi->zmin;
    dj = cpj->norm[0]*cpj->xmin + cpj->norm[1]*cpj->ymin + cpj->norm[2]*cpj->zmin;

    if(di>dj)return -1;  // sort in descending order
    if(di<dj)return 1;
  }
  else{
    if(cpi->xmax>cpj->xmin)return -1;
    if(cpi->ymax>cpj->ymin)return -1;
    if(cpi->zmax>cpj->zmin)return -1;
    return 1;
  }
  return 0;
}

/* ------------------ initcull ------------------------ */

void initcull(int cullflag){
  culldata *culli;
  int nx, ny, nz;
  int i, j, k, ii;
  int ibeg, iend, jbeg, jend, kbeg, kend;
  float xbeg, xend, ybeg, yend, zbeg, zend;
  int iskip, jskip, kskip;

  update_initcull=0;
  update_initcullplane=1;
  FREEMEMORY(cullplaneinfo);
  FREEMEMORY(sort_cullplaneinfo);
  ncullplaneinfo=0;

  for(ii=0;ii<nmeshes;ii++){
    meshdata *meshi;

    meshi=meshinfo+ii;

    GetCullSkips(meshi,cullflag,cull_portsize,&iskip,&jskip,&kskip);
    nx = meshi->ibar/iskip + 1;
    ny = meshi->jbar/jskip + 1;
    nz = meshi->kbar/kskip + 1;
    meshi->ncullinfo = nx*ny*nz;

    FREEMEMORY(meshi->cullinfo);
    NewMemory( (void **)&meshi->cullinfo,nx*ny*nz*sizeof(culldata));

    FREEMEMORY(meshi->cullQueryId);
    NewMemory( (void **)&meshi->cullQueryId,nx*ny*nz*sizeof(GLuint));
    culli=meshi->cullinfo;

    for(k=0;k<nz;k++){
      kbeg = k*kskip;
      kend = kbeg + kskip;
      if(kend>meshi->kbar)kend=meshi->kbar;
      zbeg = meshi->zplt[kbeg];
      zend = meshi->zplt[kend];
      for(j=0;j<ny;j++){
        jbeg = j*jskip;
        jend = jbeg + jskip;
        if(jend>meshi->jbar)jend=meshi->jbar;
        ybeg = meshi->yplt[jbeg];
        yend = meshi->yplt[jend];
        for(i=0;i<nx;i++){
          ibeg = i*iskip;
          iend = ibeg + iskip;
          if(iend>meshi->ibar)iend=meshi->ibar;
          xbeg = meshi->xplt[ibeg];
          xend = meshi->xplt[iend];

          culli->ibeg=ibeg;
          culli->iend=iend;

          culli->jbeg=jbeg;
          culli->jend=jend;

          culli->kbeg=kbeg;
          culli->kend=kend;

          culli->xbeg=xbeg;
          culli->xend=xend;

          culli->ybeg=ybeg;
          culli->yend=yend;

          culli->zbeg=zbeg;
          culli->zend=zend;

          culli->iskip=iskip;
          culli->jskip=jskip;
          culli->kskip=kskip;

          culli->npixels=0;
          culli->npixels_old=-1;

          ncullplaneinfo+=(iend-ibeg)+(jend-jbeg)+(kend-kbeg);

          culli++;
        }
      }
    }
  }
  NewMemory( (void **)&cullplaneinfo,ncullplaneinfo*sizeof(cullplanedata));
  NewMemory( (void **)&sort_cullplaneinfo,ncullplaneinfo*sizeof(cullplanedata *));
}

/* ------------------ setPixelCount ------------------------ */

void setPixelCount(void){
  int imesh,n;

  have_setpixelcount=1;
  if(show_cullports==0){
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_NORMALIZE);
    glDepthMask(GL_FALSE);
    glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);
  }

  for(imesh=0;imesh<nmeshes;imesh++){
    meshdata *meshi;
    meshi = meshinfo + imesh;

    meshi->culldefined=0;
  }
  for(n=0;n<nsmoke3dinfo;n++){
    meshdata *meshi;
    smoke3ddata *smoke3di;

    smoke3di = smoke3dinfo + n;
    if(show_cullports==0&&(smoke3di->loaded==0||smoke3di->display==0))continue;
    meshi = meshinfo + smoke3di->blocknumber;
    if(meshi->culldefined==1)continue;

    setPixelCountOrthog(meshi);
  }
  if(show_cullports==0){
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);
    glDepthMask(GL_TRUE);
    glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
  }
}

/* ------------------ setPixelCountOrthog ------------------------ */

void setPixelCountOrthog(meshdata *meshi){
  int icull;
  float x0[3], x1[3], x2[3], x3[3];
  float x4[3], x5[3], x6[3], x7[3];



    meshi->culldefined=1;
    if(show_cullports==0)glGenQueries(meshi->ncullinfo,meshi->cullQueryId);

    for(icull=0;icull<meshi->ncullinfo;icull++){
      culldata *culli;

      culli = meshi->cullinfo + icull;

/*
  stuff min and max grid data into a more convenient form
  assuming the following grid numbering scheme

       5-------6
     / |      /|
   /   |     / |
  4 -------7   |
  |    |   |   |
  Z    1---|---2
  |  Y     |  /
  |/       |/
  0--X-----3

  */
      x0[0]=culli->xbeg;
      x1[0]=culli->xbeg;
      x2[0]=culli->xend;
      x3[0]=culli->xend;
      x4[0]=culli->xbeg;
      x5[0]=culli->xbeg;
      x6[0]=culli->xend;
      x7[0]=culli->xend;

      x0[1]=culli->ybeg;
      x1[1]=culli->yend;
      x2[1]=culli->yend;
      x3[1]=culli->ybeg;
      x4[1]=culli->ybeg;
      x5[1]=culli->yend;
      x6[1]=culli->yend;
      x7[1]=culli->ybeg;

      x0[2]=culli->zbeg;
      x1[2]=culli->zbeg;
      x2[2]=culli->zbeg;
      x3[2]=culli->zbeg;
      x4[2]=culli->zend;
      x5[2]=culli->zend;
      x6[2]=culli->zend;
      x7[2]=culli->zend;

      if(show_cullports==0){
        glBeginQuery(GL_SAMPLES_PASSED,meshi->cullQueryId[icull]);
        glBegin(GL_QUADS);

        glVertex3fv(x0);
        glVertex3fv(x3);
        glVertex3fv(x7);
        glVertex3fv(x4);

        glVertex3fv(x3);
        glVertex3fv(x2);
        glVertex3fv(x6);
        glVertex3fv(x7);

        glVertex3fv(x2);
        glVertex3fv(x1);
        glVertex3fv(x5);
        glVertex3fv(x6);

        glVertex3fv(x1);
        glVertex3fv(x0);
        glVertex3fv(x4);
        glVertex3fv(x5);

        glVertex3fv(x1);
        glVertex3fv(x2);
        glVertex3fv(x3);
        glVertex3fv(x0);

        glVertex3fv(x4);
        glVertex3fv(x7);
        glVertex3fv(x6);
        glVertex3fv(x5);
        glEnd();
      }
      else{
        glBegin(GL_LINES);
        glVertex3fv(x0);
        glVertex3fv(x1);
        glVertex3fv(x4);
        glVertex3fv(x5);
        glVertex3fv(x3);
        glVertex3fv(x2);
        glVertex3fv(x7);
        glVertex3fv(x6);

        glVertex3fv(x3);
        glVertex3fv(x0);
        glVertex3fv(x7);
        glVertex3fv(x4);
        glVertex3fv(x6);
        glVertex3fv(x5);
        glVertex3fv(x2);
        glVertex3fv(x1);

        glVertex3fv(x0);
        glVertex3fv(x4);
        glVertex3fv(x3);
        glVertex3fv(x7);
        glVertex3fv(x2);
        glVertex3fv(x6);
        glVertex3fv(x1);
        glVertex3fv(x5);

        glEnd();
      }
      if(show_cullports==0)glEndQuery(GL_SAMPLES_PASSED);
    }
}

/* ------------------ getPixelCount ------------------------ */

void getPixelCount(void){
  meshdata *meshi;
  int i;
  int icull;

  for(i=0;i<nmeshes;i++){
    meshi = meshinfo + i;

    if(meshi->culldefined==0)continue;
    if(show_cullports==0&&have_setpixelcount==1){
      for(icull=0;icull<meshi->ncullinfo;icull++){
        culldata *culli;

        culli = meshi->cullinfo + icull;
        glGetQueryObjectiv(meshi->cullQueryId[icull],GL_QUERY_RESULT,&culli->npixels);
        if(culli->npixels_old!=culli->npixels){
          update_initcullplane=1;
          culli->npixels_old=culli->npixels;
        }
      }
      glDeleteQueries(meshi->ncullinfo,meshi->cullQueryId);
    }
    else{
      for(icull=0;icull<meshi->ncullinfo;icull++){
        culldata *culli;

        culli = meshi->cullinfo + icull;
        culli->npixels=1;
        if(update_initcullplane==0&&culli->npixels!=culli->npixels_old){
          culli->npixels_old=culli->npixels;
          update_initcullplane=1;
        }
      }
    }
  }
#ifdef _DEBUG
  if(update_initcullplane==1)PRINTF("pixel count has changed - 1\n");
#endif
}
#endif

/* ------------------ HaveFire ------------------------ */

int HaveFire(void){
  int i;

  for(i = 0; i<nsmoke3dinfo; i++){
    smoke3ddata *smoke3di;

    smoke3di = smoke3dinfo+i;
    if(smoke3di->loaded==1&&smoke3di->type==FIRE)return 1;
  }
  return 0;
}
