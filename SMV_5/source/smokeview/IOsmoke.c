// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#ifdef pp_GPU
#include <GL/glew.h>
#endif
#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "flowfiles.h"
#include <sys/types.h>
#include <sys/stat.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "egz_stdio.h"
#include "MALLOC.h"
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

int cull_count=0;

#ifdef pp_CULL
#define NCULLS 5
int cullplane_compare( const void *arg1, const void *arg2 );
#endif
// svn revision character string
char IOsmoke_revision[]="$Revision$";

char *textFileRead(char *fn);

//              alphaf_out[n]=adjustalpha(ALPHAIN, xyzeyeorig, xp, ASPECTRATIO, NORM, NORMTYPE);\

#define ADJUSTALPHA(ALPHAIN,ASPECTRATIO,NORM,NORMTYPE) \
            alphaf_out[n]=0;\
            if(ALPHAIN==0)continue;\
            if(adjustalphaflag==2&&iblank_smoke3d!=NULL&&iblank_smoke3d[n]==0)continue;\
            if(adjustalphaflag==3){\
              alphaf_out[n]=ALPHAIN;\
            }\
            else{\
              alphaf_out[n]=adjustalpha(ALPHAIN, ASPECTRATIO);\
            }

#define DRAWVERTEX(XX,YY,ZZ)        \
if(show_smoketest==0){\
  value[0]=alphaf_ptr[n11]; \
  value[1]=alphaf_ptr[n12]; \
  value[2]=alphaf_ptr[n22]; \
  value[3]=alphaf_ptr[n21]; \
  if(value[0]==0&&value[1]==0&&value[2]==0&&value[3]==0)continue;\
  ivalue[0]=n11<<2;  \
  ivalue[1]=n12<<2;  \
  ivalue[2]=n22<<2;  \
  ivalue[3]=n21<<2;  \
  if(abs(value[0]-value[2])<abs(value[1]-value[3])){     \
    xyzindex=xyzindex1;                                  \
  }                                                      \
  else{                                                  \
    xyzindex=xyzindex2;                                  \
  }                                                      \
  for(node=0;node<6;node++){                             \
    mm = xyzindex[node];                                 \
    alphabyte = value[mm];                               \
    switch (skip){                                       \
      case 1:\
       break;\
      case 2:                                              \
        alphaval=alphabyte/255.0; \
        alphaval=alphaval*(2.0-alphaval);                  \
        alphabyte=alphaval*255.0;\
       break;                                             \
      case 3:                                              \
        alphaval=alphabyte/255.0;                              \
        alphaval = 3*alphaval*(1.0-alphaval*(1.0-alphaval/3.0));\
        alphabyte = 255*alphaval; \
        break;                                                   \
      default:\
        ASSERT(FFALSE);\
        break;\
    }                                                          \
    colorptr=mergecolorptr+ivalue[mm];\
    colorptr[3]=alphabyte;                                   \
    glColor4ubv(colorptr);                                \
    glVertex3f(XX,YY,ZZ);                                \
  }\
}\
    else{\
      for(node=0;node<6;node++){                             \
        mm = xyzindex1[node];                                 \
        glColor4ub(0,0,0,(unsigned char)smoke_alpha);\
        glVertex3f(XX,YY,ZZ);                                \
}\
  }

#define DRAWVERTEXTERRAIN(XX,YY,ZZ)        \
if(show_smoketest==0){\
  value[0]=alphaf_ptr[n11]; \
  value[1]=alphaf_ptr[n12]; \
  value[2]=alphaf_ptr[n22]; \
  value[3]=alphaf_ptr[n21]; \
  if(value[0]==0&&value[1]==0&&value[2]==0&&value[3]==0)continue;\
  z_offset[0]=znode_offset[m11];\
  z_offset[1]=znode_offset[m12];\
  z_offset[2]=znode_offset[m22];\
  z_offset[3]=znode_offset[m21];\
  ivalue[0]=n11<<2;  \
  ivalue[1]=n12<<2;  \
  ivalue[2]=n22<<2;  \
  ivalue[3]=n21<<2;  \
  if(abs(value[0]-value[2])<abs(value[1]-value[3])){     \
    xyzindex=xyzindex1;                                  \
  }                                                      \
  else{                                                  \
    xyzindex=xyzindex2;                                  \
  }                                                      \
  for(node=0;node<6;node++){                             \
    mm = xyzindex[node];                                 \
    alphabyte = value[mm];                               \
    switch (skip){                                       \
      case 1:\
       break;\
      case 2:                                              \
        alphaval=alphabyte/255.0; \
        alphaval=alphaval*(2.0-alphaval);                  \
        alphabyte=alphaval*255.0;\
       break;                                             \
      case 3:                                              \
        alphaval=alphabyte/255.0;                              \
        alphaval = 3*alphaval*(1.0-alphaval*(1.0-alphaval/3.0));\
        alphabyte = 255*alphaval; \
        break;                                                   \
      default:\
        ASSERT(FFALSE);\
        break;\
    }                                                          \
    colorptr=mergecolorptr+ivalue[mm];\
    colorptr[3]=alphabyte;                                   \
    glColor4ubv(colorptr);                                \
    glVertex3f(XX,YY,ZZ+z_offset[mm]);                                \
  }\
}\
else{\
  z_offset[0]=znode_offset[m11];\
  z_offset[1]=znode_offset[m12];\
  z_offset[2]=znode_offset[m22];\
  z_offset[3]=znode_offset[m21];\
  for(node=0;node<6;node++){                             \
    mm = xyzindex1[node];                                 \
    glColor4ub(0,0,0,(unsigned char)smoke_alpha);\
    glVertex3f(XX,YY,ZZ+z_offset[mm]);                                \
  }\
}

#ifdef pp_GPU_BLANK
#define GPU_BLANK glVertexAttrib1f(GPU_blank,(float)[mm]);
#define SETBVALS bvalue[0]=(float)iblank_smoke3d[n11];\
  bvalue[1]=(float)iblank_smoke3d[n12];\
  bvalue[2]=(float)iblank_smoke3d[n22];\
  bvalue[3]=(float)iblank_smoke3d[n21];
#else
#define GPU_BLANK
#define SETBVALS
#endif

#define DRAWVERTEXGPU(XX,YY,ZZ) \
  value[0]=alphaf_in[n11];\
  value[1]=alphaf_in[n12];\
  value[2]=alphaf_in[n22];\
  value[3]=alphaf_in[n21];\
  if(value[0]==0&&value[1]==0&&value[2]==0&&value[3]==0)continue;\
  SETBVALS \
  if((adjustalphaflag==2||adjustalphaflag==3)&&iblank_smoke3d!=NULL){\
    if(iblank_smoke3d[n11]==0)value[0]=0;\
    if(iblank_smoke3d[n12]==0)value[1]=0;\
    if(iblank_smoke3d[n22]==0)value[2]=0;\
    if(iblank_smoke3d[n21]==0)value[3]=0;\
  }\
  if(value[0]==0&&value[1]==0&&value[2]==0&&value[3]==0)continue;\
  if(abs(value[0]-value[2])<abs(value[1]-value[3])){     \
    xyzindex=xyzindex1;                                  \
  }                                                      \
  else{                                                  \
    xyzindex=xyzindex2;                                  \
}                                                        \
  if(firecolor==NULL){\
    for(node=0;node<6;node++){                             \
      mm = xyzindex[node];                                 \
      GPU_BLANK \
      glVertexAttrib1f(GPU_smokealpha,(float)value[mm]); \
      glVertex3f(XX,YY,ZZ);                                \
    }\
  }\
  else{\
    fvalue[0]=firecolor[n11];\
    fvalue[1]=firecolor[n12];\
    fvalue[2]=firecolor[n22];\
    fvalue[3]=firecolor[n21];\
    for(node=0;node<6;node++){                             \
      mm = xyzindex[node];                                 \
      GPU_BLANK \
      glVertexAttrib1f(GPU_smokealpha,(float)value[mm]);\
      glVertexAttrib1f(GPU_hrr,(float)fvalue[mm]);\
      glVertex3f(XX,YY,ZZ);                                \
    }\
}
#define DRAWVERTEXGPUTERRAIN(XX,YY,ZZ) \
  z_offset[0]=znode_offset[m11];\
  z_offset[1]=znode_offset[m12];\
  z_offset[2]=znode_offset[m22];\
  z_offset[3]=znode_offset[m21];\
  value[0]=alphaf_in[n11];\
  value[1]=alphaf_in[n12];\
  value[2]=alphaf_in[n22];\
  value[3]=alphaf_in[n21];\
  SETBVALS \
  if((adjustalphaflag==2||adjustalphaflag==3)&&iblank_smoke3d!=NULL){\
    if(iblank_smoke3d[n11]==0)value[0]=0;\
    if(iblank_smoke3d[n12]==0)value[1]=0;\
    if(iblank_smoke3d[n22]==0)value[2]=0;\
    if(iblank_smoke3d[n21]==0)value[3]=0;\
  }\
  if(value[0]==0&&value[1]==0&&value[2]==0&&value[3]==0)continue;\
  if(abs(value[0]-value[2])<abs(value[1]-value[3])){     \
    xyzindex=xyzindex1;                                  \
  }                                                      \
  else{                                                  \
    xyzindex=xyzindex2;                                  \
}                                                        \
  if(firecolor==NULL){\
    for(node=0;node<6;node++){                             \
      mm = xyzindex[node];                                 \
      GPU_BLANK \
      glVertexAttrib1f(GPU_smokealpha,(float)value[mm]); \
      glVertex3f(XX,YY,ZZ+z_offset[mm]);                                \
    }\
  }\
  else{\
    fvalue[0]=firecolor[n11];\
    fvalue[1]=firecolor[n12];\
    fvalue[2]=firecolor[n22];\
    fvalue[3]=firecolor[n21];\
    for(node=0;node<6;node++){                             \
      mm = xyzindex[node];                                 \
      GPU_BLANK \
      glVertexAttrib1f(GPU_smokealpha,(float)value[mm]);\
      glVertexAttrib1f(GPU_hrr,(float)fvalue[mm]);\
      glVertex3f(XX,YY,ZZ+z_offset[mm]);                                \
    }\
}

/* ------------------ readsmoke3d ------------------------ */

void readsmoke3d(int ifile,int flag, int *errorcode){
  smoke3d *smoke3di,*smoke3dj;
  EGZ_FILE *SMOKE3DFILE;
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

  float time;
  char compstring[128];
  mesh *meshi;

  local_starttime0 = glutGet(GLUT_ELAPSED_TIME);
  smoke3di = smoke3dinfo + ifile;

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
    plotstate=getplotstate(DYNAMIC_PLOTS);
    updatetimes();
    Read3DSmoke3DFile=0;
    setsmokecolorflags();

    hrrpuv_loaded=0;
    for(j=0;j<nsmoke3d_files;j++){
      smoke3dj = smoke3dinfo + j;
      if(smoke3dj->loaded==1){
        Read3DSmoke3DFile=1;
        if(strcmp(smoke3dj->label.longlabel,"HRRPUV")==0){
          hrrpuv_loaded=1;
        }
        break;
      }
    }
    update_makeiblank_smoke3d=1;
#ifdef pp_CULL
    if(cullactive==1)update_initcull=1;
#endif
    return;
  }
  CheckMemory;
  if(getsmoke3d_sizes(smoke3di->file,smoke3di->version,&smoke3di->times, &smoke3di->use_smokeframe,
                 &smoke3di->nchars_uncompressed, 
                 &smoke3di->nchars_compressed_smoke,
                 &smoke3di->nchars_compressed_smoke_full,
                 &smoke3di->n_times,
                 &smoke3di->n_times_full)==1){
    readsmoke3d(ifile,UNLOAD,&error);
    *errorcode=1;
    printf("*** error: problems sizing 3d smoke data for %s\n",smoke3di->file);
    return;
  }
  CheckMemory;
  if(
     NewMemory((void **)&smoke3di->smokeframe_comp_list,smoke3di->n_times_full*sizeof(unsigned char *))==0||
     NewMemory((void **)&smoke3di->smoke_state_list,smoke3di->n_times_full*sizeof(unsigned char))==0||
     NewMemory((void **)&smoke3di->smokeframe_in,smoke3di->nchars_uncompressed*sizeof(unsigned char))==0||
     NewMemory((void **)&smoke3di->smokeview_tmp,smoke3di->nchars_uncompressed*sizeof(unsigned char))==0||
     NewMemory((void **)&smoke3di->smokeframe_out,smoke3di->nchars_uncompressed*sizeof(unsigned char))==0||
     NewMemory((void **)&meshi->merge_color,4*smoke3di->nchars_uncompressed*sizeof(unsigned char))==0||
     NewMemory((void **)&meshi->merge_alpha,smoke3di->nchars_uncompressed*sizeof(unsigned char))==0){
     readsmoke3d(ifile,UNLOAD,&error);
     *errorcode=1;
     printf("*** error: problems allocating memory for 3d smoke file: %s\n",smoke3di->file);
     return;
  }
  for(i=0;i<smoke3di->n_times_full;i++){
    smoke3di->smoke_state_list[i]=2;
  }

  ncomp_smoke_total=0;
  ncomp_smoke_total_skipped=0;
  for(i=0;i<smoke3di->n_times_full;i++){
    ncomp_smoke_total+=smoke3di->nchars_compressed_smoke_full[i];
    if(smoke3di->use_smokeframe[i]==1){
      ncomp_smoke_total_skipped+=smoke3di->nchars_compressed_smoke_full[i];
    }
  }
  if(NewMemory((void **)&smoke3di->smoke_comp_all,ncomp_smoke_total_skipped*sizeof(unsigned char))==0){
    readsmoke3d(ifile,UNLOAD,&error);
    *errorcode=1;
     printf("*** error: problems allocating memory for 3d smoke file: %s\n",smoke3di->file);
    return;
  }
  smoke3di->ncomp_smoke_total=ncomp_smoke_total_skipped;

  ncomp_smoke_total=0;
  i=0;
  for(ii=0;ii<smoke3di->n_times_full;ii++){
    if(smoke3di->use_smokeframe[ii]==1){
      smoke3di->smokeframe_comp_list[i]=smoke3di->smoke_comp_all+ncomp_smoke_total;
      ncomp_smoke_total+=smoke3di->nchars_compressed_smoke[i];
      i++;
    }
  }

// read in data

  getfile_size(smoke3di->file,&file_size);
#ifdef EGZ
  SMOKE3DFILE=EGZ_FOPEN(smoke3di->file,"rb",0,2);
#endif
#ifndef EGZ
  SMOKE3DFILE=EGZ_FOPEN(smoke3di->file,"rb");
#endif
  if(SMOKE3DFILE==NULL){
    readsmoke3d(ifile,UNLOAD,&error);
    *errorcode=1;
    return;
  }

  EGZ_FREAD(nxyz,4,8,SMOKE3DFILE);
  smoke3di->is1=nxyz[2];
  smoke3di->is2=nxyz[3];
  smoke3di->js1=nxyz[4];
  smoke3di->js2=nxyz[5];
  smoke3di->ks1=nxyz[6];
  smoke3di->ks2=nxyz[7];
  smoke3di->version=nxyz[1];

  // read smoke data

  iii=0;
  nframes_found=0;
  local_starttime = glutGet(GLUT_ELAPSED_TIME);
  for(i=0;i<smoke3di->n_times_full;i++){
	EGZ_FREAD(&time,4,1,SMOKE3DFILE);
    if(EGZ_FEOF(SMOKE3DFILE)!=0||(use_tload_end==1&&time>tload_end)){
      smoke3di->n_times_full=i;
      smoke3di->n_times=nframes_found;
      break;
    }
    if(use_tload_begin==1&&time<tload_begin)smoke3di->use_smokeframe[i]=0;
    if(smoke3di->use_smokeframe[i]==1){
      printf("3D smoke/fire time=%.2f",time);
    }
    EGZ_FREAD(nchars,4,2,SMOKE3DFILE);
    if(EGZ_FEOF(SMOKE3DFILE)!=0){
      smoke3di->n_times_full=i;
      smoke3di->n_times=nframes_found;
      break;
    }
    if(smoke3di->use_smokeframe[i]==1){
      float complevel;

      nframes_found++;
      EGZ_FREAD(smoke3di->smokeframe_comp_list[iii],1,smoke3di->nchars_compressed_smoke[iii],SMOKE3DFILE);
      iii++;
      CheckMemory;
      if(EGZ_FEOF(SMOKE3DFILE)!=0){
        smoke3di->n_times_full=i;
        smoke3di->n_times=nframes_found;
        break;
      }

      complevel=40.0*(float)nchars[0]/(float)nchars[1];
      complevel=(int)complevel;
      complevel/=10.0;
      sprintf(compstring," compression ratio: %.1f",complevel);
      trim(compstring);
      trimzeros(compstring);
      printf("%s\n",compstring);
    }
    else{
      EGZ_FSEEK(SMOKE3DFILE,smoke3di->nchars_compressed_smoke_full[i],SEEK_CUR);
      if(EGZ_FEOF(SMOKE3DFILE)!=0){
        smoke3di->n_times_full=i;
        smoke3di->n_times=nframes_found;
        break;
      }
    }
  }
  local_stoptime = glutGet(GLUT_ELAPSED_TIME);
  delta_time = (local_stoptime-local_starttime)/1000.0;

  if(SMOKE3DFILE!=NULL){
    EGZ_FCLOSE(SMOKE3DFILE);
  }

  smoke3di->loaded=1;
  smoke3di->display=1;
  setsmokecolorflags();
  if(strcmp(smoke3di->label.longlabel,"HRRPUV")==0){
    hrrpuv_loaded=1;
  }

  Read3DSmoke3DFile=1;
  makeiblank_smoke3d();
  plotstate=getplotstate(DYNAMIC_PLOTS);
  updatetimes();
#ifdef pp_CULL
    if(cullactive==1)initcull(cullsmoke);
#endif
  IDLE();
  local_stoptime0 = glutGet(GLUT_ELAPSED_TIME);
  delta_time0=(local_stoptime0-local_starttime0)/1000.0;

  if(file_size!=0&&delta_time>0.0){
    float loadrate;

    loadrate = ((float)file_size*8.0/1000000.0)/delta_time;
    printf(" %.1f MB loaded in %.2f s - rate: %.1f Mb/s (overhead: %.2f s)\n",
    (float)file_size/1000000.,delta_time,loadrate,delta_time0-delta_time);
  }
  else{
    printf(" %.1f MB downloaded in %.2f s (overhead: %.2f s)",
    (float)file_size/1000000.,delta_time,delta_time0-delta_time);
  }
}

/* ------------------ setsmokecolorflags ------------------------ */

void setsmokecolorflags(void){
  int i,j;
  smoke3d *smoke3dj,*smoke3di;

  for(i=0;i<nsmoke3d_files;i++){
    smoke3di = smoke3dinfo + i;
    smoke3di->soot_loaded=0;
    smoke3di->water_loaded=0;
    smoke3di->hrrpuv_loaded=0;
  }
  
  for(i=0;i<nsmoke3d_files;i++){
    smoke3di=smoke3dinfo + i;
    smoke3di->water_color=NULL;
    smoke3di->hrrpuv_color=NULL;
    smoke3di->soot_color=NULL;
    smoke3di->water_index=-1;
    smoke3di->hrrpuv_index=-1;
    smoke3di->soot_index=-1;
    if(smoke3di->loaded==0)continue;

    switch (smoke3di->type){
    case 1:
      smoke3di->soot_color=smoke3di->smokeframe_in;
      smoke3di->soot_index=i;
      break;
    case 2:
      smoke3di->hrrpuv_color=smoke3di->smokeframe_in;
      smoke3di->hrrpuv_index=i;
      break;
    case 3:
      smoke3di->water_color=smoke3di->smokeframe_in;
      smoke3di->water_index=i;
      break;
    default:
      ASSERT(FFALSE);
      break;
    }

    for(j=0;j<nsmoke3d_files;j++){
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

      switch (smoke3dj->type){
      case 1:
        smoke3di->soot_color=smoke3dj->smokeframe_in;
        smoke3di->soot_index=j;
        break;
      case 2:
        smoke3di->hrrpuv_color=smoke3dj->smokeframe_in;
        smoke3di->hrrpuv_index=j;
        break;
      case 3:
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
}

/* ------------------ getsmoke3d_sizes ------------------------ */
int getsmoke3d_sizes(char *smokefile, int version, float **timelist_found, int **use_smokeframe,
                      int *nchars_smoke_uncompressed, 
                      int **nchars_smoke_compressed_found,
                      int **nchars_smoke_compressed_full,
                      int *n_times_found,int *n_times_full){
  char textfilename[1024],textfilename2[1024],buffer[255];
  char *textptr;
  FILE *TEXTFILE=NULL;
  EGZ_FILE *SMOKE3DFILE;
  int nframes_found;
  float time,time_max,*time_found=NULL;
  int *use_smokeframe_full;
  int nch_uncompressed,nch_smoke_compressed;
  int *nch_smoke_compressed_full=NULL;
  int *nch_smoke_compressed_found=NULL;
  int nxyz[8];
  int nchars[2];
  int skip;
  int iframe;
  int dummy;
  int n_times_full2;
  size_t lentext;
  int iii;

  // try .sz
  strcpy(textfilename,smokefile);
  lentext=strlen(textfilename);
  if(lentext>4){
    textptr=textfilename + lentext - 4;
    if(strcmp(textptr,".svz")==0)textfilename[lentext-4]=0;
  }
  strcat(textfilename,".szz");
  TEXTFILE=fopen(textfilename,"r");

  if(TEXTFILE==NULL){
    // not .szz so try .sz
    strcpy(textfilename,smokefile);
    strcat(textfilename,".sz");
    TEXTFILE=fopen(textfilename,"r");
  }

  if(TEXTFILE==NULL){
    // neither .sz or .szz so try in tempdir
    strcpy(textfilename,smokefile);
    strcat(textfilename,".sz");
    TEXTFILE=fopen(textfilename,"w");
    if(TEXTFILE==NULL&&smokeviewtempdir!=NULL){
      strcpy(textfilename2,smokeviewtempdir);
      strcat(textfilename2,textfilename);
      strcpy(textfilename,textfilename2);
      TEXTFILE=fopen(textfilename,"w");
    }
    if(TEXTFILE==NULL)return 1;  // can't write size file in temp directory so give up
#ifdef EGZ
    SMOKE3DFILE=EGZ_FOPEN(smokefile,"rb",0,2);
#endif
#ifndef EGZ
    SMOKE3DFILE=EGZ_FOPEN(smokefile,"rb");
#endif
    if(SMOKE3DFILE==NULL){
      fclose(TEXTFILE);
      return 1;
    }

    EGZ_FREAD(nxyz,4,8,SMOKE3DFILE);

    if(version!=1)version=0;
    fprintf(TEXTFILE,"%i\n",version);

    time_max=-1000000.0;
    for(;;){
      EGZ_FREAD(&time,4,1,SMOKE3DFILE);
      if(EGZ_FEOF(SMOKE3DFILE)!=0)break;
      fprintf(TEXTFILE,"%f",time);
      EGZ_FREAD(nchars,4,2,SMOKE3DFILE);
      if(version==0){
        fprintf(TEXTFILE," %i %i\n",nchars[0],nchars[1]);
      }
      else{
        fprintf(TEXTFILE," %i -1 %i \n",nchars[0],nchars[1]);
      }
      skip = nchars[1];
      EGZ_FSEEK(SMOKE3DFILE,skip,SEEK_CUR);
    }

    EGZ_FCLOSE(SMOKE3DFILE);
    fclose(TEXTFILE);
    TEXTFILE=fopen(textfilename,"r");
  }
  if(TEXTFILE==NULL)return 1;

  nframes_found=0;
  iframe=-1;
  time_max=-1000000.0;
  fgets(buffer,255,TEXTFILE);
  iii=0;
  while(!feof(TEXTFILE)){
    if(fgets(buffer,255,TEXTFILE)==NULL)break;
    sscanf(buffer,"%f",&time);
    iframe++;
    if(time<=time_max)continue;
    if(use_tload_end==1&&time>tload_end)break;
    if(iii%smoke3dframestep==0&&(use_tload_begin==0||time>=tload_begin)){
      nframes_found++;
      time_max=time;
    }
    iii++;
  }
  rewind(TEXTFILE);
  if(nframes_found<=0){
    *n_times_found=0;
    *n_times_full=0;
    fclose(TEXTFILE);
    return 1;
  }
  *n_times_found=nframes_found;
  *n_times_full=iframe+1;
  if(nframes_found>0){
    NewMemory((void **)&use_smokeframe_full,*n_times_full*sizeof(float));
    NewMemory((void **)&time_found,nframes_found*sizeof(float));
    NewMemory((void **)&nch_smoke_compressed_full,(*n_times_full)*sizeof(int));
    NewMemory((void **)&nch_smoke_compressed_found,(*n_times_found)*sizeof(int));
  }
  *use_smokeframe=use_smokeframe_full;
  *timelist_found=time_found;
  *nchars_smoke_compressed_full=nch_smoke_compressed_full;
  *nchars_smoke_compressed_found=nch_smoke_compressed_found;

  fgets(buffer,255,TEXTFILE);
  n_times_full2=0;
  time_max=-1000000.0;
  iii=0;
  while(!feof(TEXTFILE)){

    if(fgets(buffer,255,TEXTFILE)==NULL)break;
    n_times_full2++;
    if(n_times_full2>*n_times_full)break;
    if(version==0){
      sscanf(buffer,"%f %i %i",&time,&nch_uncompressed,&nch_smoke_compressed);
    }
    else{
      sscanf(buffer,"%f %i %i %i",&time,&nch_uncompressed,&dummy,&nch_smoke_compressed);
    }
    *nch_smoke_compressed_full++=nch_smoke_compressed;
    *use_smokeframe_full=0;
    if(use_tload_end==1&&time>tload_end)break;
    if(time<=time_max){
      use_smokeframe_full++;
      continue;
    }

    if(iii%smoke3dframestep==0&&(use_tload_begin==0||time>=tload_begin)){
      *use_smokeframe_full=1;
      *time_found++=time;
      time_max=time;
      *nch_smoke_compressed_found++=nch_smoke_compressed;
    }
    use_smokeframe_full++;
    iii++;
  }
  *nchars_smoke_uncompressed=nch_uncompressed;
  fclose(TEXTFILE);
  return 0;
}

/* ------------------ freesmoke3d ------------------------ */

void freesmoke3d(smoke3d *smoke3di){

  smoke3di->lastiframe=-999;
  FREEMEMORY(smoke3di->smokeframe_in);
  FREEMEMORY(smoke3di->smokeframe_out);
  FREEMEMORY(smoke3di->timeslist);
  FREEMEMORY(smoke3di->times);
  FREEMEMORY(smoke3di->use_smokeframe);
  FREEMEMORY(smoke3di->nchars_compressed_smoke_full);
  FREEMEMORY(smoke3di->nchars_compressed_smoke);
  FREEMEMORY(smoke3di->smoke_state_list);
  FREEMEMORY(smoke3di->smoke_comp_all);
  FREEMEMORY(smoke3di->smokeframe_comp_list);
  FREEMEMORY(smoke3di->smokeview_tmp);
}

/* ------------------ updatesmoke3d ------------------------ */

void updatesmoke3d(smoke3d *smoke3di){
  int iframe;
  int countin;
#ifdef USE_ZLIB
  uLongf countout;
#else
  int countout;
#endif

  iframe = smoke3di->iframe;
  countin = smoke3di->nchars_compressed_smoke[iframe];
  countout=smoke3di->nchars_uncompressed;
  switch(smoke3di->version){
  case 0:
    countout = irle(smoke3di->smokeframe_comp_list[iframe],countin,smoke3di->smokeframe_in);
    break;
  case 1:
#ifdef USE_ZLIB
    uncompress(
      smoke3di->smokeframe_in,&countout,
      smoke3di->smokeframe_comp_list[iframe],countin);
#endif
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  if(smoke3di->smoke_state_list[iframe]==2){
    int i;
    unsigned char *smokeframe_in;

    smokeframe_in = smoke3di->smokeframe_in;
    smoke3di->smoke_state_list[iframe]=0;
    for(i=0;i<smoke3di->nchars_uncompressed;i++){
      if(smokeframe_in[i]==0)continue;
      smoke3di->smoke_state_list[iframe]=1;
      break;
    }
  }
  ASSERT(countout==smoke3di->nchars_uncompressed);
}

/* ------------------ mergesmokecolors ------------------------ */

void mergesmoke3dcolors(void){
  int i,j;
  smoke3d *smoke3di,*smoke3dref;
  unsigned char smoke[3]; //{0,0,0};
  unsigned char fire[3]; //{255,128,0};
  unsigned char water[3]={205,205,255};

  unsigned char *firecolor,*watercolor,*sootcolor;
  unsigned char *mergecolor,*mergealpha;
  int i_hrrpuv_cutoff;
  float fire_alpha;

  smoke3d *smoke_soot;

  mesh *meshi;

  smoke[0]=smoke_shade;
  smoke[1]=smoke_shade;
  smoke[2]=smoke_shade;
  fire[0]=fire_red;
  fire[1]=fire_green;
  fire[2]=fire_blue;
  i_hrrpuv_cutoff=254*hrrpuv_cutoff/hrrpuv_max_smv;

#ifdef pp_CULL
  for(i=0;i<nmeshes;i++){
    meshi = meshinfo + i;
    meshi->cull_smoke3d=NULL;
  }
#endif

  for(i=0;i<nsmoke3d_files;i++){
    smoke3di=smoke3dinfo + i;
    smoke3di->d_display=0;
    if(smoke3di->loaded==0||smoke3di->display==0)continue;
    smoke_soot=NULL;
    switch (smoke3di->type){
    case 1:
      smoke3di->d_display=1;
      break;
    case 2:
      if(smoke3di->soot_index!=-1)smoke_soot=smoke3dinfo + smoke3di->soot_index;
      if(smoke3di->soot_loaded==0||(smoke_soot!=NULL&&smoke_soot->display==0)){
        smoke3di->d_display=1;
      }
      break;
    case 3:
      if(smoke3di->soot_loaded==0&&smoke3di->hrrpuv_loaded==0){
        smoke3di->d_display=1;
      }
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
  }

  for(i=0;i<nsmoke3d_files;i++){
    smoke3di = smoke3dinfo + i;
    if(smoke3di->loaded==0)continue;
    if(smoke3di->d_display==0)continue;
    meshi=meshinfo+smoke3di->blocknumber;
#ifdef pp_CULL
    meshi->cull_smoke3d=smoke3di;
#endif
    i_hrrpuv_cutoff=254*meshi->hrrpuv_cutoff/hrrpuv_max_smv;

    if(fire_halfdepth<=0.0){
      smoke3di->fire_alpha=255;
    }
    else{
      smoke3di->fire_alpha=255*(1.0-pow(0.5,meshi->dx/fire_halfdepth));
    }
    fire_alpha=smoke3di->fire_alpha;

    firecolor=NULL;
    watercolor=NULL;
    sootcolor=NULL;
    if(smoke3di->hrrpuv_color!=NULL){
      firecolor=smoke3di->hrrpuv_color;
      if(smoke3di->hrrpuv_index!=-1){
        smoke3dref=smoke3dinfo + smoke3di->hrrpuv_index;
        if(smoke3dref->display==0)firecolor=NULL;
      }
    }
    if(smoke3di->water_color!=NULL){
      watercolor=smoke3di->water_color;
      if(smoke3di->water_index!=-1){
        smoke3dref=smoke3dinfo + smoke3di->water_index;
        if(smoke3dref->display==0)watercolor=NULL;
      }
    }
    if(smoke3di->soot_color!=NULL){
      sootcolor=smoke3di->soot_color;
      if(smoke3di->soot_index!=-1){
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
      if(firecolor==NULL&&watercolor==NULL){
        for(j=0;j<smoke3di->nchars_uncompressed;j++){
          *mergecolor++=smoke[0];
          *mergecolor++=smoke[1];
          *mergecolor++=smoke[2];
          mergecolor++;
          *mergealpha=sootcolor[j]>>smoke3d_thick;
          mergealpha++;
        }
        continue;
      }
      if(firecolor!=NULL&&watercolor==NULL){
        for(j=0;j<smoke3di->nchars_uncompressed;j++){
          if(firecolor[j]>i_hrrpuv_cutoff){
            *mergecolor++=fire[0];
            *mergecolor++=fire[1];
            *mergecolor++=fire[2];
            mergecolor++;
            *mergealpha++=fire_alpha;
          }
          else{
            *mergecolor++=smoke[0];
            *mergecolor++=smoke[1];
            *mergecolor++=smoke[2];
            mergecolor++;
            *mergealpha++=sootcolor[j]>>smoke3d_thick;
          }
        }
        continue;
      }
      if(firecolor==NULL&&watercolor!=NULL){
        for(j=0;j<smoke3di->nchars_uncompressed;j++){
          if(watercolor[j]>20){
            *mergecolor++=water[0];
            *mergecolor++=water[1];
            *mergecolor++=water[2];
            mergecolor++;
            *mergealpha++=watercolor[j];
          }
          else{
            *mergecolor++=smoke[0];
            *mergecolor++=smoke[1];
            *mergecolor++=smoke[2];
            mergecolor++;
            *mergealpha++=sootcolor[j]>>smoke3d_thick;
          }
        }
        continue;
      } 
      if(firecolor!=NULL&&watercolor!=NULL){
        for(j=0;j<smoke3di->nchars_uncompressed;j++){
          if(watercolor[j]>20){
            *mergecolor++=water[0];
            *mergecolor++=water[1];
            *mergecolor++=water[2];
            mergecolor++;
            *mergealpha++=watercolor[j];
            continue;
          }
          if(firecolor[j]>i_hrrpuv_cutoff){
            *mergecolor++=fire[0];
            *mergecolor++=fire[1];
            *mergecolor++=fire[2];
            mergecolor++;
            *mergealpha++=fire_alpha;
            continue;
          }
          *mergecolor++=smoke[0];
          *mergecolor++=smoke[1];
          *mergecolor++=smoke[2];
          mergecolor++;
          *mergealpha++=sootcolor[j]>>smoke3d_thick;
          continue;
        }
      }
    }
    if(sootcolor==NULL){
      if(firecolor==NULL&&watercolor==NULL){
        ASSERT(FFALSE);
        // this should never happen!!!
      }
      if(firecolor!=NULL&&watercolor==NULL){
        for(j=0;j<smoke3di->nchars_uncompressed;j++){
          *mergecolor++=fire[0];
          *mergecolor++=fire[1];
          *mergecolor++=fire[2];
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
      if(firecolor==NULL&&watercolor!=NULL){
        for(j=0;j<smoke3di->nchars_uncompressed;j++){
          *mergecolor++=water[0];
          *mergecolor++=water[1];
          *mergecolor++=water[2];
          mergecolor++;
          if(watercolor[j]>20){
            *mergealpha++=watercolor[j];
          }
          else{
            *mergealpha++=0;
          }
        }
        continue;
      } 
      if(firecolor!=NULL&&watercolor!=NULL){
        for(j=0;j<smoke3di->nchars_uncompressed;j++){
          mergecolor[0]=0;
          mergecolor[1]=0;
          mergecolor[2]=0;
          *mergealpha=0;
          if(watercolor[j]>20){
            *mergecolor++=water[0];
            *mergecolor++=water[1];
            *mergecolor++=water[2];
            mergecolor++;
            *mergealpha++=watercolor[j];
            continue;
          }
          if(firecolor[j]>i_hrrpuv_cutoff){
            *mergecolor++=fire[0];
            *mergecolor++=fire[1];
            *mergecolor++=fire[2];
            mergecolor++;
            *mergealpha++=fire_alpha;
            continue;
          }
          mergecolor+=4;
          mergealpha++;

        }
      }
    }
  }
}

/* ------------------ drawsmoke3d ------------------------ */

void drawsmoke3d(smoke3d *smoke3di){
  int i,j,k,n;
  float constval,x1,x3,z1,z3, yy1, y3;
  int is1, is2, js1, js2, ks1, ks2;
  int ii, jj, kk;
  int ibeg, iend, jbeg, jend, kbeg, kend;
  float *znode_offset, z_offset[4];
  
  float *xplt, *yplt, *zplt;
  unsigned char mergealpha,*mergealphaptr,*mergecolorptr;
  int nx,ny,nz;
  unsigned char *alphaf_in,*alphaf_out,*alphaf_ptr;
  float alphaval;
  unsigned char alphabyte;
  unsigned char *colorptr;
  int xyzindex1[6],xyzindex2[6],*xyzindex,node,mm;
  float xnode[4],znode[4],ynode[4];
  int skip;
    int iterm, jterm, kterm,nxy;
  float x11[3], x12[3], x22[3], x21[3];
  int n11, n12, n22, n21;
  int ipj,jpk,ipk,jmi,kmi,kmj;
  int iii, jjj, kkk;
  int slice_end,slice_beg;
  float aspectratio;
  int ssmokedir;
  unsigned char *iblank_smoke3d;

  unsigned char value[4];
  int ivalue[4];

  mesh *meshi;

  meshi = meshinfo + smoke3di->blocknumber;
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
  alphaf_out=smoke3di->smokeframe_out;

  switch (demo_mode){
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
  skip=smokeskipm1+1;

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

  transparenton();
  switch (ssmokedir){
 
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
    for(iii=slice_beg;iii<slice_end;iii+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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
    for(jjj=slice_beg;jjj<slice_end;jjj+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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
    for(kkk=slice_beg;kkk<slice_end;kkk+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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

      for(iii=1;iii<nx+ny-2;iii+=skip){
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
    for(iii=slice_beg;iii<slice_end;iii+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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

      for(iii=1;iii<nx+ny-2;iii+=skip){
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
    for(iii=slice_beg;iii<slice_end;iii+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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

      for(iii=1;iii<ny+nz-2;iii+=skip){
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
    for(iii=slice_beg;iii<slice_end;iii+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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

      for(iii=1;iii<ny+nz-2;iii+=skip){
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
    for(iii=slice_beg;iii<slice_end;iii+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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

      for(iii=1;iii<nx+nz-2;iii+=skip){
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
    for(iii=slice_beg;iii<slice_end;iii+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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

      for(iii=1;iii<nx+nz-2;iii+=skip){
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
    for(iii=slice_beg;iii<slice_end;iii+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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
  transparentoff();
  if(cullfaces==1)glEnable(GL_CULL_FACE);
//  printf("majorcull=%i minorcull=%i\n",majorcull,minorcull);


}
#ifdef pp_GPU
/* ------------------ drawsmoke3dGPU ------------------------ */

void drawsmoke3dGPU(smoke3d *smoke3di){
  int i,j,k,n;
  float constval,x1,x3,z1,z3, yy1, y3;
  int is1, is2, js1, js2, ks1, ks2;
  int ii, jj, kk;
  int ibeg, iend, jbeg, jend, kbeg, kend;

  float *xplt, *yplt, *zplt;
  float *znode_offset, z_offset[4];

  int nx,ny,nz;
  unsigned char *alphaf_in;
  int xyzindex1[6],xyzindex2[6],*xyzindex,node,mm;
  float xnode[4],znode[4],ynode[4];
  int skip;
  int iterm, jterm, kterm,nxy;
  float x11[3], x12[3], x22[3], x21[3];
  int n11, n12, n22, n21;
  int ipj,jpk,ipk,jmi,kmi,kmj;
  int iii, jjj, kkk;
  int slice_end,slice_beg;
  float aspectratio;
  int ssmokedir;
  unsigned char *firecolor;
  unsigned char *iblank_smoke3d;
  int is_smoke;

  unsigned char value[4];
  float bvalue[4];
  unsigned char fvalue[4];

  mesh *meshi;
  float fire_alpha;

  meshi = meshinfo + smoke3di->blocknumber;
  firecolor=smoke3di->hrrpuv_color;

  {
    smoke3d *sooti=NULL;

	if(smoke3di->soot_index>=0){
	  sooti = smoke3dinfo + smoke3di->soot_index;
	}
	if(sooti!=NULL&&sooti->display==1){
      is_smoke=1;
    }
    else{
      is_smoke=0;
    }
  }
  iblank_smoke3d = meshi->iblank_smoke3d;
  if(fire_halfdepth<=0.0){
    fire_alpha=256.0;
  }
  else{
    fire_alpha=256*(1.0-pow(0.5,meshi->dx/fire_halfdepth));
  }

  meshi = meshinfo + smoke3di->blocknumber;
  value[0]=255;
  value[1]=255;
  value[2]=255;
  value[3]=255;
  smoke_shade4[0]=smoke_shade/255.0;
  smoke_shade4[1]=smoke_shade/255.0;
  smoke_shade4[2]=smoke_shade/255.0;
  smoke_shade4[3]=1.0;

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
  skip=smokeskipm1+1;

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
  glUniform1i(GPU_is_smoke,is_smoke);
  {
    float fire_color[4];

    fire_color[0]=fire_red/256;
    fire_color[1]=fire_green/256;
    fire_color[2]=fire_blue/256;
    fire_color[3]=fire_alpha/256.0;
    glUniform4fv(GPU_firecolor,1,fire_color);
  }
  glUniform1f(GPU_smokeshade,(float)smoke_shade);
 // glUniform1i(GPU_skip,skip);
  glUniform1f(GPU_smoke3d_rthick,smoke3d_rthick);
  if(firecolor==NULL){
    i_hrrcutoff=-1;
  }
  else{
    i_hrrcutoff=254*meshi->hrrpuv_cutoff/hrrpuv_max_smv;
    if(i_hrrcutoff<0)i_hrrcutoff=0;
    if(i_hrrcutoff>254)i_hrrcutoff=254;
  }
  glUniform1f(GPU_hrrcutoff,(float)i_hrrcutoff);
  transparenton();
  switch (ssmokedir){
  // +++++++++++++++++++++++++++++++++++ DIR 1 +++++++++++++++++++++++++++++++++++++++


  case 1:
  case -1:

    aspectratio=meshi->dx;
    glUniform1f(GPU_aspectratio,aspectratio);

    // ++++++++++++++++++  draw triangles +++++++++++++++++

    glBegin(GL_TRIANGLES);
    slice_beg=is1;
    slice_end=is2;
    for(iii=slice_beg;iii<slice_end;iii+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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
    for(jjj=slice_beg;jjj<slice_end;jjj+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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
    for(kkk=slice_beg;kkk<slice_end;kkk+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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
    for(iii=slice_beg;iii<slice_end;iii+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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
    for(iii=slice_beg;iii<slice_end;iii+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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
    for(iii=slice_beg;iii<slice_end;iii+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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
    for(iii=slice_beg;iii<slice_end;iii+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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
    for(iii=slice_beg;iii<slice_end;iii+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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
    for(iii=slice_beg;iii<slice_end;iii+=skip){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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
  transparentoff();
  if(cullfaces==1)glEnable(GL_CULL_FACE);
//  printf("majorcull=%i minorcull=%i\n",majorcull,minorcull);


}
#endif
#ifdef pp_CULL

/* ------------------ init_cull_exts ------------------------ */

int init_cull_exts(void){
  char version_label[256];
  char version_label2[256];
  int i, major,  minor;
  float version;
  const GLubyte *version_string;
  int err;

  cullactive=0;
  version_string=glGetString(GL_VERSION);
  if(version_string==NULL){
    printf("*** warning: GL_VERSION string is NULL in init_cull_exts()\n");
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
    trim(version_label);
    printf("Smokeview is running on a system using OpenGL %s\n",version_label2);
    printf("Smoke culling is not not supported, it requires OpenGL 2.0 or later.\n");
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
  int ii, jj, kk;
  int ibeg, iend, jbeg, jend, kbeg, kend;
  float *znode_offset, z_offset[4];

  float *xplt, *yplt, *zplt;
  int nx,ny;
  unsigned char *alphaf_in;
  int xyzindex1[6],xyzindex2[6],*xyzindex,node,mm;
  float xnode[4],znode[4],ynode[4];
  int iterm, jterm, kterm,nxy;
  float x11[3], x12[3], x22[3], x21[3];
  int n11, n12, n22, n21;
  int ipj,jpk,ipk,jmi,kmi,kmj;
  int iii, jjj, kkk;
  //int slice_end,slice_beg;
  float aspectratio;
  unsigned char *firecolor;
  int ntemp;
  unsigned char *iblank_smoke3d;
  int is_smoke;

  unsigned char value[4];
#ifdef pp_GPU_BLANK
  unsigned char bvalue[4];
#endif
  unsigned char fvalue[4];

  mesh *mesh_old;
  float fire_alpha;

  CheckMemory;
  if(cullfaces==1)glDisable(GL_CULL_FACE);
  transparenton();

  sniffErrors("before drawsmoke3dcull");
  value[0]=255;
  value[1]=255;
  value[2]=255;
  value[3]=255;
  smoke_shade4[0]=smoke_shade/255.0;
  smoke_shade4[1]=smoke_shade/255.0;
  smoke_shade4[2]=smoke_shade/255.0;
  smoke_shade4[3]=1.0;

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

  glUniform1i(GPU_adjustalphaflag,adjustalphaflag);
  glUniform1f(GPU_smokeshade,(float)smoke_shade);
  glUniform1f(GPU_smoke3d_rthick,smoke3d_rthick);

  CheckMemory;
  sniffErrors("in drawsmoke3dcull 1");
  glBegin(GL_TRIANGLES);
  mesh_old=NULL;
  for(nn=0;nn<ncullplaneinfo;nn++){
    mesh *meshi;
    smoke3d *smoke3di;

    cullplanedata *culli;

    culli = sort_cullplaneinfo[ncullplaneinfo-1-nn];

    meshi = culli->cull_mesh;


  
    if(meshi!=mesh_old){
      mesh_old=meshi;
  
      glEnd();

      if(nterraininfo>0){
        znode_offset = meshi->terrain->znode_offset;
      }
      sniffErrors("in drawsmoke3dcull 4");
      smoke3di=meshi->cull_smoke3d;

      firecolor=smoke3di->hrrpuv_color;
      if(fire_halfdepth<=0.0){
        fire_alpha=256.0;
      }
      else{
        fire_alpha=256*(1.0-pow(0.5,meshi->dx/fire_halfdepth));
      }
      {
        float fire_color[4];
        float *ffirecolor;

        fire_color[0]=(float)fire_red/256.0;
        fire_color[1]=(float)fire_green/256.0;
        fire_color[2]=(float)fire_blue/256.0;
        fire_color[3]=(float)fire_alpha/256.0;
        ffirecolor=fire_color;
        glUniform4fv(GPU_firecolor,1,ffirecolor);
      }

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

      if(firecolor==NULL){
        i_hrrcutoff=-1;
      }
      else{
        i_hrrcutoff=254*meshi->hrrpuv_cutoff/hrrpuv_max_smv;
        if(i_hrrcutoff<0)i_hrrcutoff=0;
        if(i_hrrcutoff>254)i_hrrcutoff=254;
      }
      glUniform1f(GPU_hrrcutoff,(float)i_hrrcutoff);
      switch (meshi->smokedir){
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
        }
	    {
	      smoke3d *sooti=NULL;

        if(smoke3di->soot_index>=0){
          sooti = smoke3dinfo + smoke3di->soot_index;
        }
        if(sooti!=NULL&&sooti->display==1){
          is_smoke=1;
        }
        else{
          is_smoke=0;
        }
	    }
      glUniform1f(GPU_aspectratio,aspectratio);
      glUniform1i(GPU_is_smoke,is_smoke);
      glBegin(GL_TRIANGLES);
    }

    if(smoke3di->smoke_state_list[smoke3di->iframe]==0)continue;
    switch (meshi->smokedir){

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

            if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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

            if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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
      //  sniffErrors("after drawsmokecull case 2");
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

            if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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

          if(nterraininfo>0&&fabs(vertical_factor-1.0)>0.01){
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
  sniffErrors("in drawsmoke3dcull 12");
  transparentoff();
  if(cullfaces==1)glEnable(GL_CULL_FACE);
  sniffErrors("after drawsmokecull");
}

#endif

/* ------------------ updatesmoke3dmenulabels ------------------------ */

void updatesmoke3dmenulabels(void){
  int i;
  smoke3d *smoke3di;
  char meshlabel[128];

  for(i=0;i<nsmoke3d_files;i++){
    smoke3di = smoke3dinfo + i;
    STRCPY(smoke3di->menulabel,smoke3di->label.longlabel);
    smoke3di->version=getsmoke3dversion(smoke3di);
    if(showfiles==1){
      STRCAT(smoke3di->menulabel,", ");
      STRCAT(smoke3di->menulabel,smoke3di->file);
    }
    switch(smoke3di->version){
    case 0:
      STRCAT(smoke3di->menulabel," (RLE) ");
      break;
    case 1:
      STRCAT(smoke3di->menulabel," (ZLIB) ");
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
//    len=strlen(smoke3di->menulabel);
    if(nmeshes>1){
      sprintf(meshlabel,"Mesh: *%i",smoke3di->blocknumber+1);
      STRCAT(smoke3di->menulabel," - ");
      STRCAT(smoke3di->menulabel,meshlabel);
    }
  } 
}

/* ------------------ adjustalpha ------------------------ */

unsigned char adjustalpha(unsigned char alpha, float factor){
  double val, term, top, bottom, rr;
  int i;
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

/* ------------------ getsmokedir ------------------------ */

void getsmokedir(float *mm){
    /*
      ( m0 m4 m8  m12 ) (x)    (0)
      ( m1 m5 m9  m13 ) (y)    (0)
      ( m2 m6 m10 m14 ) (z)  = (0)
      ( m3 m7 m11 m15 ) (1)    (1)

       ( m0 m4  m8 )      (m12)
   Q=  ( m1 m5  m9 )  u = (m13)
       ( m2 m6 m10 )      (m14)
      
      (Q   u) (x)     (0)      
      (v^T 1) (y)   = (1)
       
      m3=m7=m11=0, v^T=0, y=1   Qx+u=0 => x=-Q^Tu
    */
  int i,ii,j;
  mesh *meshj;
  float norm[3],scalednorm[3];
  float normdir[3];
  float absangle,cosangle,minangle;
  int iminangle;
  float dx, dy, dz;
  float factor;
  float pi;

  pi=4.0*atan(1.0);

  xyzeyeorig[0] = -(mm[0]*mm[12]+mm[1]*mm[13]+ mm[2]*mm[14])/mscale[0];
  xyzeyeorig[1] = -(mm[4]*mm[12]+mm[5]*mm[13]+ mm[6]*mm[14])/mscale[1];
  xyzeyeorig[2] = -(mm[8]*mm[12]+mm[9]*mm[13]+mm[10]*mm[14])/mscale[2];
  
  for(j=0;j<nmeshes;j++){
    meshj = meshinfo + j;

    minangle=1000.0;
    iminangle=-10;
    meshj->dx=meshj->xplt_orig[1]-meshj->xplt_orig[0];
    meshj->dy=meshj->yplt_orig[1]-meshj->yplt_orig[0];
    meshj->dz=meshj->zplt_orig[1]-meshj->zplt_orig[0];
    meshj->dxy=meshj->dx*meshj->dx+meshj->dy*meshj->dy;
    meshj->dxy=sqrt(meshj->dxy)/2.0;
    meshj->dxz=meshj->dx*meshj->dx+meshj->dz*meshj->dz;
    meshj->dxz=sqrt(meshj->dxz)/2.0;
    meshj->dyz=meshj->dy*meshj->dy+meshj->dz*meshj->dz;
    meshj->dyz=sqrt(meshj->dyz)/2.0;



    meshj->dy/=meshj->dx;
    meshj->dz/=meshj->dx;
    meshj->dxy/=meshj->dx;
    meshj->dxz/=meshj->dx;
    meshj->dyz/=meshj->dx;
    meshj->dx=1.0;

    if(smokedrawtest2==1){
      meshj->norm[0]=1.0;
       meshj->norm[1]=0.0;
       meshj->norm[2]=0.0;
       meshj->smokedir=1;
       continue;
    }

    for(i=-9;i<=9;i++){
      if(i==0)continue;
      ii = i;
      if(i<0)ii=-i;
      norm[0]=0.0;
      norm[1]=0.0;
      norm[2]=0.0;
      switch (ii){
      case 1:
        if(i<0)norm[0]=-1.0;
        if(i>0)norm[0]=1.0;
        break;
      case 2:
        if(i<0)norm[1]=-1.0;
        if(i>0)norm[1]=1.0;
        break;
      case 3:
        if(i<0)norm[2]=-1.0;
        if(i>0)norm[2]=1.0;
        break;
      case 4:
        dx = meshj->xplt_orig[1]-meshj->xplt_orig[0];
        dy = meshj->yplt_orig[1]-meshj->yplt_orig[0];
        factor= dx*dx+dy*dy;
        if(factor==0.0){
          factor=1.0;
        }
        else{
          factor=1.0/sqrt(factor);
        }
        if(i<0){
          norm[0]=-dy*factor;
          norm[1]=-dx*factor;
        }
        else{
          norm[0]=dy*factor;
          norm[1]=dx*factor;
        }
        break;
      case 5:
        dx = meshj->xplt_orig[1]-meshj->xplt_orig[0];
        dy = meshj->yplt_orig[1]-meshj->yplt_orig[0];
        factor= dx*dx+dy*dy;
        if(factor==0.0){
          factor=1.0;
        }
        else{
          factor=1.0/sqrt(factor);
        }
        if(i<0){
          norm[0]= dy*factor;
          norm[1]=-dx*factor;
        }
        else{
          norm[0]=-dy*factor;
          norm[1]= dx*factor;
        }
        break;
      case 6:
        dy = meshj->yplt_orig[1]-meshj->yplt_orig[0];
        dz = meshj->zplt_orig[1]-meshj->zplt_orig[0];
        factor= dz*dz+dy*dy;
        if(factor==0.0){
          factor=1.0;
        }
        else{
          factor=1.0/sqrt(factor);
        }
        if(i<0){
          norm[1]=-dz*factor;
          norm[2]=-dy*factor;
        }
        else{
          norm[1]=dz*factor;
          norm[2]=dy*factor;
        }      
        break;
      case 7:
        dy = meshj->yplt_orig[1]-meshj->yplt_orig[0];
        dz = meshj->zplt_orig[1]-meshj->zplt_orig[0];
        factor= dz*dz+dy*dy;
        if(factor==0.0){
          factor=1.0;
        }
        else{
          factor=1.0/sqrt(factor);
        }
        if(i<0){
          norm[1]= dz*factor;
          norm[2]=-dy*factor;
        }
        else{
          norm[1]=-dz*factor;
          norm[2]= dy*factor;
        }
        break;
      case 8:
        dx = meshj->xplt_orig[1]-meshj->xplt_orig[0];
        dz = meshj->zplt_orig[1]-meshj->zplt_orig[0];
        factor= dz*dz+dx*dx;
        if(factor==0.0){
          factor=1.0;
        }
        else{
          factor=1.0/sqrt(factor);
        }
        if(i<0){
          norm[0]=-dz*factor;
          norm[2]=-dx*factor;
        }
        else{
          norm[0]=dz*factor;
          norm[2]=dx*factor;
        }      
        break;
      case 9:
        dx = meshj->xplt_orig[1]-meshj->xplt_orig[0];
        dz = meshj->zplt_orig[1]-meshj->zplt_orig[0];
        factor= dx*dx+dz*dz;
        if(factor==0.0){
          factor=1.0;
        }
        else{
          factor=1.0/sqrt(factor);
        }
        if(i<0){
          norm[0]= dz*factor;
          norm[2]=-dx*factor;
        }
        else{
          norm[0]=-dz*factor;
          norm[2]= dx*factor;
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
      }
      scalednorm[0]=norm[0]*mscale[0];
      scalednorm[1]=norm[1]*mscale[1];
      scalednorm[2]=norm[2]*mscale[2];

      normdir[0] = mm[0]*scalednorm[0] + mm[4]*scalednorm[1] + mm[8]*scalednorm[2];
      normdir[1] = mm[1]*scalednorm[0] + mm[5]*scalednorm[1] + mm[9]*scalednorm[2];
      normdir[2] = mm[2]*scalednorm[0] + mm[6]*scalednorm[1] + mm[10]*scalednorm[2];

      cosangle = normdir[2]/sqrt(normdir[0]*normdir[0]+normdir[1]*normdir[1]+normdir[2]*normdir[2]);
      if(cosangle>1.0)cosangle=1.0;
      if(cosangle<-1.0)cosangle=-1.0;
      absangle=acos(cosangle)*180.0/pi;
      if(absangle<0.0)absangle=-absangle;
      if(absangle<minangle){
        iminangle=i;
        minangle=absangle;
        meshj->norm[0]=norm[0];
        meshj->norm[1]=norm[1];
        meshj->norm[2]=norm[2];
      }
    }
    meshj->smokedir=iminangle;
#ifdef pp_CULL
    if(meshj->smokedir!=meshj->smokedir_old){
      meshj->smokedir_old=meshj->smokedir;
      update_initcullplane=1;
#ifdef _DEBUG
      printf("mesh dir has changed\n");
#endif
    }
#endif
    if(demo_mode!=0){
      meshj->smokedir=1;
    }
  }
}

/* ------------------ filter_smoke ------------------------ */

void filter_smoke(smoke3d *smoke3di){
  unsigned char *smoke_data, *smoke_tmp;
  int ndata;
  int is1, is2, js1, js2, ks1, ks2;
  int i, j, k, m;
  int nx, ny, nxy;
  int ijk;
  unsigned char minval;
  int ncount=0,ncount2=0;

  int skipi[4], skipj[4], skipk[4];
  unsigned int minvals[6];

  smoke_data=smoke3di->smokeframe_in;
  smoke_tmp=smoke3di->smokeview_tmp;
  ndata = smoke3di->nchars_uncompressed;
  is1 = smoke3di->is1;
  is2 = smoke3di->is2;
  js1 = smoke3di->js1;
  js2 = smoke3di->js2;
  ks1 = smoke3di->ks1;
  ks2 = smoke3di->ks2;
  nx = is2 + 1 - is1;
  ny = js2 + 1 - js1;
  nxy = nx*ny;
  skipi[0]=1;
  skipi[1]=2;
  skipi[2]=3;
  skipi[3]=4;
  skipj[0]=nx;
  skipj[1]=2*nx;
  skipj[2]=3*nx;
  skipj[3]=4*nx;
  skipk[0]=nxy;
  skipk[1]=2*nxy;
  skipk[2]=3*nxy;
  skipk[3]=4*nxy;



  memcpy(smoke_tmp,smoke_data,ndata);

  for(k=ks1+4;k<ks2-3;k++){
    for(j=js1+4;j<js2-3;j++){
      for(i=is1+4;i<is2-3;i++){
        ijk = (i-is1) + (j-js1)*nx + (k-ks1)*nx*ny;
        for(m=0;m<6;m++){
          minvals[m]=1;
        }
        for(m=0;m<4;m++){
          minvals[0]*=(255-smoke_tmp[ijk+skipi[m]]);
          minvals[1]*=(255-smoke_tmp[ijk-skipi[m]]);
          minvals[2]*=(255-smoke_tmp[ijk+skipj[m]]);
          minvals[3]*=(255-smoke_tmp[ijk-skipj[m]]);
          minvals[4]*=(255-smoke_tmp[ijk+skipk[m]]);
          minvals[5]*=(255-smoke_tmp[ijk-skipk[m]]);
        }
        minval=255;
        for(m=0;m<6;m++){
          minvals[m]=minvals[m]>>24;
          minvals[m]=255-minvals[m];
          if(minval<minvals[m])minval=minvals[m];
        }
        if(minval>240){
          smoke_data[ijk]=0;
          ncount++;
        }
        if(smoke_data[ijk]!=0)ncount2++;
      }
    }
  }

  printf("ncount=%i ncount2=%i\n",ncount,ncount2);

}

/* ------------------ makeiblanksmoke3d ------------------------ */

#define ijknode(i,j,k) ((i)+(j)*nx+(k)*nxy)

void makeiblank_smoke3d(void){
  smoke3d *smoke3di;
  mesh *smokemesh;
  int ibar, jbar, kbar, ijksize;
  unsigned char *iblank_smoke3d;
  blockagedata *bc;
  int ii;
  int i, j, k;
  int ic;
  int nx, ny, nxy;
  int ijk;
  float *xplt, *yplt, *zplt;
  float x, y, z;
  float dx, dy, dz;

  update_makeiblank_smoke3d=0;
  for(i=0;i<nmeshes;i++){
    smokemesh = meshinfo + i;
    smokemesh->smokeloaded=0;
  }
  for(i=0;i<nsmoke3d_files;i++){
    smoke3di = smoke3dinfo + i;
    smokemesh = meshinfo + smoke3di->blocknumber;

    if(smoke3di->loaded==1){
      smokemesh->smokeloaded=1;
    }

    ibar = smokemesh->ibar;
    jbar = smokemesh->jbar;
    kbar = smokemesh->kbar;
    ijksize=(ibar+1)*(jbar+1)*(kbar+1);
    nx = ibar + 1;
    ny = jbar + 1;
    nxy = nx*ny;

    if(use_iblank==1&&smoke3di->loaded==1&&smokemesh->iblank_smoke3d==NULL){
      NewMemory( (void **)&iblank_smoke3d,ijksize*sizeof(unsigned char));
      smokemesh->iblank_smoke3d=iblank_smoke3d;
    }
  }

#define ALLMESHES 0
#define LOWERMESHES 1

  for(ic=nmeshes-1;ic>=0;ic--){
    smokemesh = meshinfo + ic;
    if(smokemesh->smokeloaded==0)continue;
    iblank_smoke3d = smokemesh->iblank_smoke3d;
    if(iblank_smoke3d==NULL)continue;

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
      *iblank_smoke3d++=1;
    }

    iblank_smoke3d=smokemesh->iblank_smoke3d;

    for(i=0;i<=smokemesh->ibar;i++){
    for(j=0;j<=smokemesh->jbar;j++){
    for(k=0;k<=smokemesh->kbar;k++){
      ijk = ijknode(i,j,k);
      x = xplt[i];
      y = yplt[j];
      z = zplt[k];
      if(inmesh_smoke(x,y,z,ic-1,LOWERMESHES)>=0)iblank_smoke3d[ijk]=0;
    }
    }
    }

    for(ii=0;ii<smokemesh->nbptrs;ii++){
      bc=smokemesh->blockageinfoptrs[ii];
      if(bc->invisible==1||bc->hidden==1||bc->nshowtime!=0)continue;
      for(i=bc->ijk[IMIN];i<bc->ijk[IMAX];i++){
      for(j=bc->ijk[JMIN];j<bc->ijk[JMAX];j++){
      for(k=bc->ijk[KMIN];k<bc->ijk[KMAX];k++){
        ijk = ijknode(i,j,k);
        iblank_smoke3d[ijk]=0;
      }
      }
      }
    }

    for(j=0;j<=jbar;j++){
    for(k=0;k<=kbar;k++){
      ijk = ijknode(0,j,k);
      x = xplt[0];
      y = yplt[j];
      z = zplt[k];
      if(inmesh_smoke(x-dx,y,z,ic,ALLMESHES)<0){
        iblank_smoke3d[ijk]=0;
      }
      else{
        iblank_smoke3d[ijk]=1;
      }

      ijk = ijknode(ibar,j,k);
      x = xplt[ibar];
      y = yplt[j];
      z = zplt[k];
      if(inmesh_smoke(x+dx,y,z,ic,ALLMESHES)<0){
        iblank_smoke3d[ijk]=0;
      }
      else{
        iblank_smoke3d[ijk]=1;
      }

    }
    }

    for(i=0;i<=ibar;i++){
    for(k=0;k<=kbar;k++){
      ijk = ijknode(i,0,k);
      x = xplt[i];
      y = yplt[0];
      z = zplt[k];
      if(inmesh_smoke(x,y-dy,z,ic,ALLMESHES)<0){
        iblank_smoke3d[ijk]=0;
      }
      else{
        iblank_smoke3d[ijk]=1;
      }



      ijk = ijknode(i,jbar,k);
      x = xplt[i];
      y = yplt[jbar];
      z = zplt[k];
      if(inmesh_smoke(x,y+dy,z,ic,ALLMESHES)<0){
        iblank_smoke3d[ijk]=0;
      }
      else{
        iblank_smoke3d[ijk]=1;
      }


    }
    }
    for(i=0;i<=ibar;i++){
    for(j=0;j<=jbar;j++){
      ijk = ijknode(i,j,0);
      x = xplt[i];
      y = yplt[j];
      z = zplt[0];
      if(inmesh_smoke(x,y,z-dz,ic,ALLMESHES)<0){
        iblank_smoke3d[ijk]=0;
      }
      else{
        iblank_smoke3d[ijk]=1;
      }


      ijk = ijknode(i,j,kbar);
      x = xplt[i];
      y = yplt[j];
      z = zplt[kbar];
      if(inmesh_smoke(x,y,z+dz,ic,ALLMESHES)<0){
        iblank_smoke3d[ijk]=0;
      }
      else{
        iblank_smoke3d[ijk]=1;
      }


    }
    }
    
  }
}

/* ------------------ inmesh_smoke ------------------------ */

 int inmesh_smoke(float x, float y, float z, int nm, int flag){
  int i;
  mesh *meshi;
  int n;
  float xmin, ymin, zmin;
  float xmax, ymax, zmax;
  int ibar, jbar, kbar;

  if(flag==LOWERMESHES){
    n = nm;
  }
  else{
    n = nmeshes;
  }
  for(i=0;i<n;i++){
    meshi = meshinfo + i;
    if(flag==ALLMESHES&&i==nm)continue;
    if(meshi->smokeloaded==0)continue;

    ibar=meshi->ibar;
    jbar=meshi->jbar;
    kbar=meshi->kbar;
    xmin=meshi->xplt[0];
    ymin=meshi->yplt[0];
    zmin=meshi->zplt[0];
    xmax=meshi->xplt[ibar];
    ymax=meshi->yplt[jbar];
    zmax=meshi->zplt[kbar];

    if(x<xmin)continue;
    if(y<ymin)continue;
    if(z<zmin)continue;

    if(x>xmax)continue;
    if(y>ymax)continue;
    if(z>zmax)continue;
    return i;
  }
  return -1;
}

/* ------------------ getsmoke3dversion ------------------------ */

 int getsmoke3dversion(smoke3d *smoke3di){

   STRUCTSTAT statbuffer;
   EGZ_FILE *SMOKE3DFILE=NULL;
   int nxyz[8];
   char *file;

   if(STAT(smoke3di->comp_file,&statbuffer)==0){
     smoke3di->file=smoke3di->comp_file;
   }
   else{
     if(STAT(smoke3di->reg_file,&statbuffer)==0){
       smoke3di->file=smoke3di->reg_file;
     }
     else{
       return -1;
     }
   }
   file = smoke3di->file;

   SMOKE3DFILE=EGZ_FOPEN(file,"rb",0,2);
   if(SMOKE3DFILE==NULL)return -1;

   EGZ_FREAD(nxyz,4,8,SMOKE3DFILE);
   EGZ_FCLOSE(SMOKE3DFILE);


   return nxyz[1];
 }

#ifdef pp_CULL

/* ------------------ get_cullskips ------------------------ */

void get_cullskips(mesh *meshi, int cullflag, int *iiskip, int *jjskip, int *kkskip){
  int iskip, jskip, kskip;

  if(cullflag==1){
    iskip = cull_portsize;
    if(iskip<3)iskip=3;
    if(iskip>meshi->ibar+1)iskip=meshi->ibar+1;
    jskip = cull_portsize;
    if(jskip<3)jskip=3;
    if(jskip>meshi->jbar+1)jskip=meshi->jbar+1;

    kskip = cull_portsize;
    if(kskip<3)kskip=3;
    if(kskip>meshi->kbar+1)kskip=meshi->kbar+1;
  }
  else{
    iskip = meshi->ibar+1;
    jskip = meshi->jbar+1;
    kskip = meshi->kbar+1;
  }
  *iiskip=iskip;
  *jjskip=jskip;
  *kkskip=kskip;
}

    /* ------------------ initcullplane ------------------------ */

void initcullplane(int cullflag){
  int ii;
  int iskip, jskip, kskip;
  int nx, ny, nz;
  int nxx, nyy, nzz;
  int ibeg, iend, jbeg, jend, kbeg, kend;
  int i, j, k, ik, ij;
  cullplanedata *cp;
  smoke3d *smoke3di;
  float norm[3];
  float dx, dy, dz, factor;

#ifdef _DEBUG
  printf("updating initcullplane %i\n",cull_count++);
#endif
  update_initcullplane=0; 
  cp = cullplaneinfo;
  ncullplaneinfo=0;

  for(ii=0;ii<nmeshes;ii++){
    mesh *meshi;
    culldata *culli;

    meshi=meshinfo+ii;

    meshi->culldefined=0;
    if(meshi->cull_smoke3d==NULL)continue;
    smoke3di = meshi->cull_smoke3d;
    if(smoke3di->loaded==0||smoke3di->display==0)continue;
    meshi->culldefined=1;

    culli=meshi->cullinfo;
    get_cullskips(meshi,cullflag, &iskip,&jskip,&kskip);

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
          int jj, kk;

          ibeg = i*iskip;
          iend = ibeg + iskip;
          if(iend>meshi->ibar)iend=meshi->ibar;

          if(culli->npixels>0||cullflag==0){
            cp->norm[0]=0.0;
            cp->norm[1]=0.0;
            cp->norm[2]=0.0;
            switch (meshi->smokedir) {
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

  for(ii=0;ii<ncullplaneinfo;ii++){
    sort_cullplaneinfo[ii]=cullplaneinfo+ii;
  }
  if(ncullplaneinfo>1){
    qsort((cullplanedata *)sort_cullplaneinfo,(size_t)ncullplaneinfo,
       sizeof(cullplanedata *),cullplane_compare);
  }
  GLUTPOSTREDISPLAY
}

/* ------------------ cullplane_compare ------------------------ */

int cullplane_compare( const void *arg1, const void *arg2 ){
  cullplanedata *cpi, *cpj;

  cpi = *(cullplanedata **)arg1;
  cpj = *(cullplanedata **)arg2;

  if(cpi->dir==cpj->dir){
    float di, dj;
    float dx, dy, dz;

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
  cullplanedata *cpx, *cpy, *cpz;

  update_initcull=0;
  update_initcullplane=1;
  FREEMEMORY(cullplaneinfo);
  FREEMEMORY(sort_cullplaneinfo);
  ncullplaneinfo=0;

  for(ii=0;ii<nmeshes;ii++){
    mesh *meshi;

    meshi=meshinfo+ii;

    get_cullskips(meshi,cullflag,&iskip,&jskip,&kskip);
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
  int imesh,icull,n;

  have_setpixelcount=1;
  if(show_cullports==0){
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_NORMALIZE);
    glDepthMask(GL_FALSE);
    glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);
  }

  for(imesh=0;imesh<nmeshes;imesh++){
    mesh *meshi;
    meshi = meshinfo + imesh;

    meshi->culldefined=0;
  }
  for(n=0;n<nsmoke3d_files;n++){
    mesh *meshi;
    smoke3d *smoke3di;

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

void setPixelCountOrthog(mesh *meshi){
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
  mesh *meshi;
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
  if(update_initcullplane==1)printf("pixel count has changed - 1\n");
#endif
}
#endif
