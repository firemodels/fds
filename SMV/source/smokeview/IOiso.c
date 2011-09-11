// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include "flowfiles.h"
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <math.h>
#include "egz_stdio.h"
#include "MALLOC.h"
#include "smokeviewvars.h"
#include "update.h"

// svn revision character string
char IOiso_revision[]="$Revision$";

void sync_isobounds(int isottype);
void unloadiso(mesh *gb);

/* ------------------ getisoheader ------------------------ */

void getisoheader(const char *isofile, EGZ_FILE **isostreamptr,
                  const char *isosizefile, int *nvertices, int *ntriangles, int *nbuffer, int *maxfullbuffer, int *maxcompbuffer,
                  float **levelsptr, int *nisolevels, int *nisosteps, int isoframestep, short **normaltable, int *nnormaltable){
  FILE *isosizestream=NULL;
  float time, time_max;
  int nvert, ntri, nbuf;
  char buffer[1024];
  int istep=0;
  int one,nnn;
  int mfullbuffer;
  int mcompbuffer;
  int nfull;

  isosizestream=fopen(isosizefile,"r");
  *nisosteps=0;
  *nisolevels=0;
  *nvertices=0;
  *ntriangles=0;
  *nbuffer=0;
  *maxfullbuffer=0;
  if(isosizestream==NULL)return;
  if(fgets(buffer,255,isosizestream)==NULL)return;
  sscanf(buffer,"%i",nisolevels);

  time_max=-1000000.0;
  for(;;){
    int i;
    int skip_frame;
    
    if(fgets(buffer,255,isosizestream)==NULL)break;
    sscanf(buffer,"%f",&time);
    skip_frame=1;
    if(time>time_max){
      skip_frame=0;
      time_max=time;
    }

    mfullbuffer=0;
    mcompbuffer=0;
    for(i=0;i<*nisolevels;i++){
      if(fgets(buffer,255,isosizestream)==NULL)break;
      if(istep%isoframestep!=0||(settmin_i==1&&time<tmin_i)||(settmax_i==1&&time>tmax_i)||skip_frame==1){
      }
      else{
        sscanf(buffer,"%i %i %i %i ",&nvert,&ntri,&nfull,&nbuf);
        *nvertices+=nvert;
        *ntriangles+=ntri;
        *nbuffer+=nbuf;
        mcompbuffer+=nbuf;
        mfullbuffer+=nfull;
        time_max=time;
      }
    }
    if(istep%isoframestep!=0||(settmin_i==1&&time<tmin_i)||(settmax_i==1&&time>tmax_i)||skip_frame==1){
    }
    else{
      (*nisosteps)++;
      if(mfullbuffer>(*maxfullbuffer))*maxfullbuffer=mfullbuffer;
      if(mcompbuffer>(*maxcompbuffer))*maxcompbuffer=mcompbuffer;
    }
    istep++;
  }

#ifdef EGZ
  *isostreamptr=EGZ_FOPEN(isofile,"rb",0,2);
#else
  *isostreamptr=EGZ_FOPEN(isofile,"rb");
#endif
  if(*levelsptr==NULL){
    if(NewMemory((void **)levelsptr,*nisolevels*sizeof(float))==0)return;
  }
  else{
    if(ResizeMemory((void **)levelsptr,*nisolevels*sizeof(float))==0)return;
  }
  EGZ_FREAD(&one,4,1,*isostreamptr);
  EGZ_FREAD(&nnn,4,1,*isostreamptr);
  EGZ_FREAD(*levelsptr,4,(unsigned int)(*nisolevels),*isostreamptr);
  EGZ_FREAD(nnormaltable,4,1,*isostreamptr);
  NewMemory((void **)normaltable,3*(*nnormaltable)*sizeof(short));
  EGZ_FREAD(*normaltable,2,3*(*nnormaltable),*isostreamptr);
}

/* ------------------ getisolevels ------------------------ */

void getisolevels(const char *isofile, int dataflag, float **levelsptr, float ***colorlevelsptr, int *nisolevels){
  int one;
  int version;
  int len[3],labellengths=0;
  int nlevels;
  EGZ_FILE *isostreamptr;
  int i;
  float **colorlevels=NULL;

#ifdef EGZ
  isostreamptr=EGZ_FOPEN(isofile,"rb",0,2);
#else
  isostreamptr=EGZ_FOPEN(isofile,"rb");
#endif

  EGZ_FREAD(&one,4,1,isostreamptr);
  if(dataflag!=0){
    EGZ_FREAD(&version,4,1,isostreamptr);
  }
  else{
    version=1;
  }
  EGZ_FREAD(len,4,3,isostreamptr);
  labellengths=len[0]+len[1]+len[2];
  EGZ_FSEEK(isostreamptr,labellengths+4,SEEK_CUR);
  EGZ_FREAD(&nlevels,4,1,isostreamptr);
  *nisolevels=nlevels;
  FREEMEMORY(*levelsptr);
  NewMemory((void **)levelsptr,nlevels*sizeof(float));
  EGZ_FREAD(*levelsptr,4,(unsigned int)(nlevels),isostreamptr);
  EGZ_FCLOSE(isostreamptr);
  NewMemory((void **)&colorlevels,nlevels*sizeof(float *));
  for(i=0;i<nlevels;i++){
    float *colorlevel;
    
    colorlevels[i]=NULL;
  }
  *colorlevelsptr=colorlevels;
  
}

/* ------------------ getisosizes ------------------------ */

void getisosizes(const char *isofile, int dataflag, EGZ_FILE **isostreamptr, int *nvertices, int *ntriangles, 
                 float **levelsptr, int *nisolevels, int *nisosteps, int isoframestep, 
                 float *tmin, float *tmax, int endian){
	int len[3],labellengths=0;
	int nlevels, n;
  int nvertices_i, ntriangles_i;
	int i;
	float time, time_max;
	int beg;
	int version;
  int one;
  int skip;
  float ttmin, ttmax;

#ifdef EGZ
  *isostreamptr=EGZ_FOPEN(isofile,"rb",0,2);
#else
  *isostreamptr=EGZ_FOPEN(isofile,"rb");
#endif

  *tmin=1000000000.;
  *tmax=-1000000000.;
  EGZ_FREAD(&one,4,1,*isostreamptr);
  if(dataflag!=0){
    EGZ_FREAD(&version,4,1,*isostreamptr);
  }
  else{
    version=1;
  }
  EGZ_FREAD(len,4,3,*isostreamptr);
  labellengths=len[0]+len[1]+len[2];
  EGZ_FSEEK(*isostreamptr,labellengths+4,SEEK_CUR);
  EGZ_FREAD(&nlevels,4,1,*isostreamptr);
  *nisolevels=nlevels;
  if(*levelsptr==NULL){
    if(NewMemory((void **)levelsptr,*nisolevels*sizeof(float))==0)return;
  }
  else{
    if(ResizeMemory((void **)levelsptr,*nisolevels*sizeof(float))==0)return;
  }
  EGZ_FREAD(*levelsptr,4,(unsigned int)(nlevels),*isostreamptr);
  *nisosteps=0; *nvertices=0; *ntriangles=0;
  beg=EGZ_FTELL(*isostreamptr);
  i=0;
  time_max=-1000000.0;
  for(;;){
    int skip_frame;

    EGZ_FREAD(&time,4,1,*isostreamptr);
    skip_frame=1;
    if(time>time_max){
      skip_frame=0;
      time_max=time;
    }
    if(EGZ_FEOF(*isostreamptr)!=0)break;
	  for(n=0;n<nlevels;n++){
      EGZ_FREAD(&nvertices_i,4,1,*isostreamptr);
      if(EGZ_FEOF(*isostreamptr)!=0)break;
	    EGZ_FREAD(&ntriangles_i,4,1,*isostreamptr);
      if(EGZ_FEOF(*isostreamptr)!=0)break;
      skip=0;
      if(nvertices_i>0){
        skip += 6*nvertices_i;
        EGZ_FSEEK(*isostreamptr,skip,SEEK_CUR);
        skip=0;
        if(dataflag==1){
          EGZ_FREAD(&ttmin,4,1,*isostreamptr);
          if(ttmin<*tmin)*tmin=ttmin;
          EGZ_FREAD(&ttmax,4,1,*isostreamptr);
          if(ttmax>*tmax)*tmax=ttmax;
          skip += 2*nvertices_i;
        }
      }
	    if(nvertices_i<256){                 /* number of triangles */
	      skip+=ntriangles_i;
      }
      else if(nvertices_i>=256&&nvertices_i<65536){
	      skip+=ntriangles_i*2;
      }
	    else{
	     skip+=ntriangles_i*4;
      }
      EGZ_FSEEK(*isostreamptr,skip,SEEK_CUR);
    }
    if(skip_frame==1)continue;
    i++;
    if(i%isoframestep!=0)continue;
    if((settmin_i==1&&time<tmin_i))continue;
    if((settmax_i==1&&time>tmax_i))continue;

    *nvertices += nvertices_i;
	*ntriangles += ntriangles_i;
	*nisosteps += 1;
  }
  EGZ_FSEEK(*isostreamptr,beg,SEEK_SET);
  if(dataflag==1&&axissmooth==1){
    smoothlabel(tmin,tmax,nrgb);
  }
}

/* ------------------ readiso ------------------------ */

void readiso(const char *file, int ifile, int flag, int *errorcode){
  extern int isoframestep;
  int itime,ilevel,itri,ivert,n,ii,iitime;
  isosurface *asurface;
  int nisopoints, nisotriangles;

  float time, time_max;
  EGZ_FILE *isostream;

  int blocknumber;
  int error;
  int skip;
  float factor, offset[3];
  float *iso_colors;
  int n_iso_colors;
  
  mesh *meshi;
  iso *ib;

  int local_starttime=0, local_stoptime=0;
  FILE_SIZE file_size=0;
  int local_starttime0=0, local_stoptime0=0;  
  float delta_time, delta_time0;

  local_starttime0 = glutGet(GLUT_ELAPSED_TIME);
  
  ib = isoinfo+ifile;
  if(ib->loaded==0&&flag==UNLOAD)return;
  blocknumber=ib->blocknumber;
  ib->isoupdate_timestep=-1;
  meshi = meshinfo+blocknumber;
  unloadiso(meshi);
  unload_iso_trans();
  ib->loaded=0;
  ib->display=0;
  plotstate=getplotstate(DYNAMIC_PLOTS);
  updatetimes();
  *errorcode = 0;

#ifdef _DEBUG
  if(flag==UNLOAD){
    printf("After iso unload: ");
    PrintAllMemoryInfo;
  }
#endif
  update_isotype();
  if(flag==UNLOAD){
    updatemenu=1;
    {
      iso *isoi;

      loaded_isomesh=get_loaded_isomesh();
      update_iso_showlevels();
    }
    return;
  }
  meshi->isofilenum=ifile;
  highlight_mesh = blocknumber;
  
  factor = (meshi->xyzmaxdiff/xyzmaxdiff)/65535.0;
  offset[0]=(meshi->xbar0-xbar0)/xyzmaxdiff;
  offset[1]=(meshi->ybar0-ybar0)/xyzmaxdiff;
  offset[2]=(meshi->zbar0-zbar0)/xyzmaxdiff;

  if(iso_ambient_ini==NULL||n_iso_ambient_ini==0){
    iso_colors=iso_ambient;
    n_iso_colors=n_iso_ambient;
  }
  else{
    iso_colors=iso_ambient_ini;
    n_iso_colors=n_iso_ambient_ini;
  }

  getisosizes(file, ib->dataflag, &isostream, &nisopoints, &nisotriangles, 
    &meshi->isolevels, &meshi->nisolevels, &meshi->nisosteps, isoframestep, 
    &ib->tmin, &ib->tmax, endian_data);

  file_size=get_filesize(file);

  if(meshi->isolevels==NULL){
    readiso("",ifile,UNLOAD,&error);
    *errorcode=1;
    return;
  }               
  ASSERT(meshi->isotimes==NULL);
  if(NewMemory((void **)&meshi->isotimes,sizeof(float)*meshi->nisosteps)==0){
    readiso("",ifile,UNLOAD,&error);
    *errorcode=1;
    return;
  }
  ASSERT(meshi->showlevels==NULL);
  if(NewMemory((void **)&meshi->showlevels,sizeof(int)*meshi->nisolevels)==0){
    *errorcode=1;
    readiso("",ifile,UNLOAD,&error);
    return;
  }
  for(ilevel=0;ilevel<meshi->nisolevels;ilevel++){
    meshi->showlevels[ilevel]=1;
  }
  isomin=meshi->isolevels[0];
  isomax=meshi->isolevels[0];
  meshi->isomin_index=0;
  meshi->isomax_index=0;
  for(ilevel=1;ilevel<meshi->nisolevels;ilevel++){
    if(meshi->isolevels[ilevel]<isomin){
      isomin=meshi->isolevels[ilevel];
      meshi->isomin_index=ilevel;
    }
    if(meshi->isolevels[ilevel]>isomax){
      isomax=meshi->isolevels[ilevel];
      meshi->isomax_index=ilevel;
    }
  }
  ASSERT(meshi->animatedsurfaces==NULL);
  if(NewMemory((void **)&meshi->animatedsurfaces,meshi->nisolevels*meshi->nisosteps*sizeof(isosurface))==0){
    *errorcode=1;
    readiso("",ifile,UNLOAD,&error);
    return;
  }
  if(ResizeMemory((void **)&meshi->isotimes,sizeof(float)*meshi->nisosteps)==0){
    *errorcode=1;
    readiso("",ifile,UNLOAD,&error);
    return;
  }

  asurface=meshi->animatedsurfaces;
  iitime=0;
  itime=0;
  time_max = -1000000.0;
  local_starttime = glutGet(GLUT_ELAPSED_TIME);
  for(;;){
    int skip_frame;
    int ntri_total;

    skip_frame=0;
    iitime++;

    EGZ_FREAD(&time,4,1,isostream);
    if(EGZ_FEOF(isostream)!=0)break;
    skip_frame=1;
    if(time>time_max){
      skip_frame=0;
      time_max=time;
    }
    meshi->isotimes[itime]=time;
    if(iitime%isoframestep!=0||(settmin_i==1&&time<tmin_i)||(settmax_i==1&&time>tmax_i)||skip_frame==1){
    }
    else{
      printf("isosurface time=%f\n",time);
    }
    ntri_total=0;
    for(ilevel=0;ilevel<meshi->nisolevels;ilevel++){
      int nvertices_i, ntriangles_i;
          
      asurface->dataflag=ib->dataflag;
        
      EGZ_FREAD(&nvertices_i,4,1,isostream);
      if(EGZ_FEOF(isostream)!=0)break;
      EGZ_FREAD(&ntriangles_i,4,1,isostream);
      if(EGZ_FEOF(isostream)!=0)break;
      asurface->niso_triangles=ntriangles_i/3;
      asurface->niso_vertices=nvertices_i;
        
      if(iitime%isoframestep!=0||(settmin_i==1&&time<tmin_i)||(settmax_i==1&&time>tmax_i)||skip_frame==1){
        skip=0;
        if(nvertices_i<=0||ntriangles_i<=0)continue;
        skip += (6*nvertices_i);
        if(ib->dataflag==1)skip += (8 + 2*nvertices_i);
  	    if(nvertices_i<256){
  	      skip += (ntriangles_i);
        }
        else if(nvertices_i>=256&&nvertices_i<65536){
  	      skip += (ntriangles_i*2);
        }
	      else{
  	      skip += (ntriangles_i*4);
        }
        EGZ_FSEEK(isostream,skip,SEEK_CUR);
        continue;
      }
      
      asurface->iso_triangles=NULL;
      asurface->iso_vertices=NULL; 
      if(nvertices_i>0){
        unsigned short *verti;
        unsigned short *vertices_i;
          
        NewMemory((void **)&asurface->iso_vertices,nvertices_i*sizeof(isovert));
        NewMemory((void **)&vertices_i,3*nvertices_i*sizeof(unsigned short));
        verti = vertices_i;
        EGZ_FREAD(vertices_i,2,(unsigned int)(3*nvertices_i),isostream);
        for(ivert=0;ivert<nvertices_i;ivert++){
          isovert *isoverti;
          float *xyz,*vertnorm;
            
          isoverti = asurface->iso_vertices+ivert;
          xyz = isoverti->xyz;
          xyz[0]=offset[0]+factor*(*verti++); 
          xyz[1]=offset[1]+factor*(*verti++); 
          xyz[2]=offset[2]+factor*(*verti++); 
          isoverti->flag=0;
            
          vertnorm=isoverti->norm;
          vertnorm[0]=0.0;
          vertnorm[1]=0.0;
          vertnorm[2]=0.0;
          if(ilevel==0&&strcmp(ib->surface_label.shortlabel,"hrrpuv")==0){
            isoverti->color=hrrpuv_iso_color;
          }
          else{
            isoverti->color=iso_colors+4*ilevel;
          }
        }
        FREEMEMORY(vertices_i);

        if(ib->dataflag==1){
          unsigned short *tvertices_i;
          float tcolor0, tcolorfactor, tcolor;
          
          EGZ_FREAD(&asurface->tmin,4,1,isostream);
          EGZ_FREAD(&asurface->tmax,4,1,isostream);
          NewMemory((void **)&tvertices_i,nvertices_i*sizeof(unsigned short));
          EGZ_FREAD(tvertices_i,2,(unsigned int)nvertices_i,isostream);
          if(ib->tmax>ib->tmin){
            tcolor0 = (asurface->tmin-ib->tmin)/(ib->tmax-ib->tmin);
            tcolorfactor = (asurface->tmax-asurface->tmin)/65535.;
            tcolorfactor /= (ib->tmax-ib->tmin);
          }
          else{
            tcolor0=0.5;
            tcolorfactor=0.0;
          }
          for(ivert=0;ivert<nvertices_i;ivert++){
            isovert *isoverti;
            unsigned char colorindex;
                          
            isoverti = asurface->iso_vertices+ivert;
            tcolor = tcolor0 + tvertices_i[ivert]*tcolorfactor;
            if(tcolor<0.0)tcolor=0.0;
            if(tcolor>1.0)tcolor=1.0;
            colorindex = (unsigned char)(tcolor*255);
            isoverti->color = rgb_iso+4*colorindex;
            isoverti->texturecolor=tcolor;
          }
          FREEMEMORY(tvertices_i);
        }
      }
      if(EGZ_FEOF(isostream)!=0)break;
      if(ntriangles_i>0){
        unsigned char *triangles1_i;
        unsigned short *triangles2_i;
        int *triangles_i;
          
        NewMemory((void **)&triangles_i,ntriangles_i*sizeof(int));
        if(nvertices_i<256&&nvertices_i>0){
          NewMemory((void **)&triangles1_i,ntriangles_i*sizeof(unsigned char));
          EGZ_FREAD(triangles1_i,1,(unsigned int)ntriangles_i,isostream);
          for(itri=0;itri<ntriangles_i;itri++){
            triangles_i[itri]=triangles1_i[itri];
          }
          FREEMEMORY(triangles1_i);
        }
        else if(nvertices_i>=256&&nvertices_i<65536){
          NewMemory((void **)&triangles2_i,ntriangles_i*sizeof(unsigned short));
          EGZ_FREAD(triangles2_i,2,(unsigned int)ntriangles_i,isostream);
          for(itri=0;itri<ntriangles_i;itri++){
            triangles_i[itri]=triangles2_i[itri];
          }
          FREEMEMORY(triangles2_i);
        }
        else{
          EGZ_FREAD(triangles_i,4,(unsigned int)ntriangles_i,isostream);
        } 
        NewMemory((void **)&asurface->iso_triangles,(ntriangles_i/3)*sizeof(isotri));
        for(itri=0;itri<ntriangles_i/3;itri++){
          isotri *isotrii;
              
          isotrii=asurface->iso_triangles+itri;
          isotrii->v1=asurface->iso_vertices+triangles_i[3*itri];
          isotrii->v2=asurface->iso_vertices+triangles_i[3*itri+1];
          isotrii->v3=asurface->iso_vertices+triangles_i[3*itri+2];
          isotrii->xyzmid[0]=(isotrii->v1->xyz[0]+isotrii->v2->xyz[0]+isotrii->v3->xyz[0])/3.0;
          isotrii->xyzmid[1]=(isotrii->v1->xyz[1]+isotrii->v2->xyz[1]+isotrii->v3->xyz[1])/3.0;
          isotrii->xyzmid[2]=(isotrii->v1->xyz[2]+isotrii->v2->xyz[2]+isotrii->v3->xyz[2])/3.0;
          isotrii->distance=-1.0;
          if(ilevel==0&&strcmp(ib->surface_label.shortlabel,"hrrpuv")==0){
            ib->colorlevels[ilevel]=hrrpuv_iso_color;
          }
          else{
            ib->colorlevels[ilevel]=iso_colors+4*ilevel;
          }
          isotrii->color=ib->colorlevels+ilevel;
          if(ib->dataflag==0){
            isotrii->v1->color=*isotrii->color;
            isotrii->v2->color=*isotrii->color;
            isotrii->v3->color=*isotrii->color;
          }
        }
        FREEMEMORY(triangles_i);
      }
        
      if(EGZ_FEOF(isostream)!=0)break;

      for(itri=0;itri<ntriangles_i/3;itri++){
        isotri *isotrii;
        float *v1, *v2, *v3;
        float *trinorm, *vertnorm;
        float area;
        float out[3];
                    
        isotrii = asurface->iso_triangles+itri;
        v1=isotrii->v1->xyz;
        v2=isotrii->v2->xyz;
        v3=isotrii->v3->xyz;
        calcNormal2f(v1,v2,v3,out,&area);
          
        trinorm=isotrii->norm;
        trinorm[0]=out[0];
        trinorm[1]=out[1];
        trinorm[2]=out[2];
        
        vertnorm = isotrii->v1->norm;
        vertnorm[0] += out[0]*area;
        vertnorm[1] += out[1]*area;
        vertnorm[2] += out[2]*area;
        
        vertnorm = isotrii->v2->norm;
        vertnorm[0] += out[0]*area;
        vertnorm[1] += out[1]*area;
        vertnorm[2] += out[2]*area;
          
        vertnorm = isotrii->v3->norm;
        vertnorm[0] += out[0]*area;
        vertnorm[1] += out[1]*area;
        vertnorm[2] += out[2]*area;
      }
      for(ivert=0;ivert<nvertices_i;ivert++){
        ReduceToUnit(asurface->iso_vertices[ivert].norm);
      }
      ntri_total+=asurface->niso_triangles;
      asurface++;
    }
   
    if(skip_frame==1||iitime%isoframestep!=0||(settmin_i==1&&time<tmin_i)||(settmax_i==1&&time>tmax_i)){
     // if(skip_frame==1)jj--;
    }
    else{
      itime++;
      if(itime>=meshi->nisosteps)break;
    }
  }

  local_stoptime = glutGet(GLUT_ELAPSED_TIME);
  delta_time = (local_stoptime-local_starttime)/1000.0;
  EGZ_FCLOSE(isostream);
  if(*errorcode!=0){
    unloadiso(meshi);
    readiso("",ifile,UNLOAD,&error);
    return;
    }


  ib->loaded=1;
  ib->display=1;
  loaded_isomesh=get_loaded_isomesh();
  update_iso_showlevels();
  ReadIsoFile=1;
  plotstate=getplotstate(DYNAMIC_PLOTS);
  updatemenu=1;
  iisotype=getisotype(ib);

  if(ib->dataflag==1){
    iisottype = getisottype(ib);
    sync_isobounds(iisottype);
    setisolabels(ib->tmin, ib->tmax, ib, errorcode);
    CheckMemory;
  }

  updatetimes();
#ifdef _DEBUG
  printf("After iso load: ");
  PrintMemoryInfo;
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

  glutPostRedisplay();
  CheckMemory;
}

/* ------------------ unloadiso_iso_trans ------------------------ */

void unload_iso_trans(void){
  if(iso_trans_list!=NULL){
    int i;

    for(i=0;i<niso_timesteps;i++){
      FREEMEMORY(iso_trans_list[i]);
    }
    FREEMEMORY(niso_trans_list);
    FREEMEMORY(iso_trans_list);
  }
  if(iso_opaques_list!=NULL){
      int i;

      for(i=0;i<niso_timesteps;i++){
        FREEMEMORY(iso_opaques_list[i]);
      }
    FREEMEMORY(niso_opaques_list);
    FREEMEMORY(iso_opaques_list);
  }

  niso_trans=0;
  niso_opaques=0;
}

/* ------------------ unloadiso ------------------------ */

void unloadiso(mesh *meshi){
  isosurface *asurface;
  int n,j;
  iso *ib;
  int nloaded=0;
  int i;
  mesh *meshi2;

  if(meshi->isofilenum==-1)return;
  ib = isoinfo + meshi->isofilenum;
  if(meshi->nisosteps>0&&meshi->nisolevels>0){
    for(i=0;i<meshi->nisosteps*meshi->nisolevels;i++){
      asurface=meshi->animatedsurfaces+i;
      FREEMEMORY(asurface->iso_triangles);
      FREEMEMORY(asurface->iso_vertices);
    }
    CheckMemoryOff;
    FREEMEMORY(meshi->isotimes);
    FREEMEMORY(meshi->animatedsurfaces);
    FREEMEMORY(meshi->showlevels);
  }
  meshi->nisosteps=0;
  FREEMEMORY(ib->normaltable);

  unload_iso_trans();

  ib->loaded=0;
  ib->display=0;
  plotstate=getplotstate(DYNAMIC_PLOTS);
  meshi->isofilenum=-1;
  for(i=0;i<nmeshes;i++){
    meshi2 = meshinfo+i;
    if(meshi2->isofilenum!=-1)nloaded++;
  }
  if(nloaded==0){
    ReadIsoFile=0;
  }

  updatetimes();
  updatemenu=1;
  IDLE();

  return;
}

/* ------------------ drawiso ------------------------ */

void drawiso(int tranflag){
  int i;
  isosurface *asurface;
  float *iso_colors;
  int n_iso_colors;
  int *showlevels, nisolevels;
  iso *isoi=NULL;
  float iso_color_tmp[4];
  float *iso_color_ptr;
  int iso_lighting;
  mesh *meshi;

  meshi = loaded_isomesh;


  CheckMemory;
  if(tranflag==DRAW_TRANSPARENT&&visAIso!=1)return;
  if(meshi->isofilenum>=0){
    isoi = isoinfo + meshi->isofilenum;
  }

  iso_lighting=1;

  showlevels=meshi->showlevels;
  nisolevels=meshi->nisolevels;

  if(iso_ambient_ini==NULL||n_iso_ambient_ini==0){
    iso_colors=iso_ambient;
    n_iso_colors=n_iso_ambient;
  }
  else{
    iso_colors=iso_ambient_ini;
    n_iso_colors=n_iso_ambient_ini;
  }

  if(visAIso==1){
    isotri **iso_list_start;
    int niso_list_start;
    float *colorptr=NULL;
    float *colorptr_old=NULL;

    asurface = meshi->animatedsurfaces + meshi->iiso*meshi->nisolevels;
    if(cullfaces==1)glDisable(GL_CULL_FACE);

    iso_specular[3] = 1.0;
    if(tranflag==DRAW_TRANSPARENT)transparenton();

    if(usetexturebar==1&&isoi->dataflag==1){
      glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
      glEnable(GL_TEXTURE_1D);
      glBindTexture(GL_TEXTURE_1D,texture_iso_colorbar_id);
    }

    glPushAttrib(GL_LIGHTING_BIT);
    if(iso_lighting==1){
      glEnable(GL_LIGHTING);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,iso_specular);
      glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,iso_shininess);
      glEnable(GL_COLOR_MATERIAL);
    }
    glBegin(GL_TRIANGLES);

    if(tranflag==DRAW_TRANSPARENT){
      iso_list_start=iso_trans;
      niso_list_start=niso_trans;
    }
    else{
      iso_list_start=iso_opaques;
      niso_list_start=niso_opaques;
    }
    CheckMemory;
    if(usetexturebar==1&&isoi->dataflag==1){
      for(i=0;i<niso_list_start;i++){
        isotri *tri;
        float *colorptr;
        isovert *v1, *v2, *v3;
        
        tri=iso_list_start[i];

        v1 = tri->v1;
        v2 = tri->v2;
        v3 = tri->v3;

        glTexCoord1f(v1->texturecolor);
        glNormal3fv(v1->norm);
        glVertex3fv(v1->xyz);
        
        glTexCoord1f(v2->texturecolor);
        glNormal3fv(v2->norm);
        glVertex3fv(v2->xyz);
        
        glTexCoord1f(v3->texturecolor);
        glNormal3fv(v3->norm);
        glVertex3fv(v3->xyz);
      }
    }
    else{
      for(i=0;i<niso_list_start;i++){
        isotri *tri;
        float *colorptr;
        isovert *v1, *v2, *v3;
        
        tri=iso_list_start[i];

        v1 = tri->v1;
        v2 = tri->v2;
        v3 = tri->v3;

        glColor4fv(v1->color);
        glNormal3fv(v1->norm);
        glVertex3fv(v1->xyz);
        
        glColor4fv(v2->color);
        glNormal3fv(v2->norm);
        glVertex3fv(v2->xyz);
        
        glColor4fv(v3->color);
        glNormal3fv(v3->norm);
        glVertex3fv(v3->xyz);
      }
    }
    glEnd();

    glPopAttrib();
    if(usetexturebar==1&&isoi->dataflag==1)glDisable(GL_TEXTURE_1D);


    if(tranflag==DRAW_TRANSPARENT)transparentoff();
    if(cullfaces==1)glEnable(GL_CULL_FACE);
    CheckMemory;
  }

  if(visAIso==2){
    asurface = meshi->animatedsurfaces + meshi->iiso*meshi->nisolevels;

    glPushAttrib(GL_LIGHTING_BIT);
    antialias(1);
    glLineWidth(isolinewidth);
    glBegin(GL_LINES);
    for(i=0;i<niso_trans;i++){
      isotri *tri;
      float *xyz1, *xyz2, *xyz3;
        
      tri=iso_trans[i];

      xyz1 = tri->v1->xyz;
      xyz2 = tri->v2->xyz;
      xyz3 = tri->v3->xyz;

      glColor3fv(*tri->color);

      glVertex3fv(xyz1);
      glVertex3fv(xyz2);
        
      glVertex3fv(xyz2);
      glVertex3fv(xyz3);
        
      glVertex3fv(xyz3);
      glVertex3fv(xyz1);
    }
    for(i=0;i<niso_opaques;i++){
      isotri *tri;
      float *xyz1, *xyz2, *xyz3;
        
      tri=iso_opaques[i];

      xyz1 = tri->v1->xyz;
      xyz2 = tri->v2->xyz;
      xyz3 = tri->v3->xyz;

      glColor3fv(*tri->color);

      glVertex3fv(xyz1);
      glVertex3fv(xyz2);
        
      glVertex3fv(xyz2);
      glVertex3fv(xyz3);
        
      glVertex3fv(xyz3);
      glVertex3fv(xyz1);
    }
    glEnd();
    antialias(0);
    glPopAttrib();
  }

  if(visAIso==3){
    asurface = meshi->animatedsurfaces + meshi->iiso*meshi->nisolevels;

    antialias(1);
    glPointSize(isopointsize);
    asurface--;
    glBegin(GL_POINTS);
    for(i=0;i<niso_trans;i++){
      isotri *tri;
      float *xyz1, *xyz2, *xyz3;
        
      tri=iso_trans[i];

      xyz1 = tri->v1->xyz;
      xyz2 = tri->v2->xyz;
      xyz3 = tri->v3->xyz;
      
      glColor3fv(*tri->color);

      glVertex3fv(xyz1);
      glVertex3fv(xyz2);
      glVertex3fv(xyz3);
    }
    for(i=0;i<niso_opaques;i++){
      isotri *tri;
      float *xyz1, *xyz2, *xyz3;
        
      tri=iso_opaques[i];

      xyz1 = tri->v1->xyz;
      xyz2 = tri->v2->xyz;
      xyz3 = tri->v3->xyz;
      
      glColor3fv(*tri->color);

      glVertex3fv(xyz1);
      glVertex3fv(xyz2);
      glVertex3fv(xyz3);
    }
    glEnd();
    antialias(0);
  }
}

/* ------------------ drawstaticiso ------------------------ */

void drawstaticiso(const isosurface *asurface,int surfacetype, 
                   int smoothnorm, int trans_flag, int data_type, 
                   float line_width){
  int j,k;
  float vv1[3],vv2[3],vv3[3];
  float vv1n[3],vv2n[3],vv3n[3];
  short *norm=NULL;
  unsigned short *v1, *v2, *v3;
  unsigned short *vertices_i=NULL;
  int *triangles_i=NULL;
  int nvertices;
  int i1, i2, i3;
  short *norm1,*norm2,*norm3,*vertexnorm=NULL;
  int ntriangles;
  float xyzmin[3], xyzmaxdiff_local;
  int drawing_transparent, drawing_blockage_transparent, drawing_vent_transparent;
  int drawing_smooth;
  int transparenton_flag=0;

  get_drawing_parms(&drawing_smooth, &drawing_transparent, &drawing_blockage_transparent, &drawing_vent_transparent);

  xyzmin[0] = asurface->xmin;
  xyzmin[1] = asurface->ymin;
  xyzmin[2] = asurface->zmin;
  xyzmaxdiff_local = asurface->xyzmaxdiff;

  nvertices=asurface->nvertices;
  ntriangles=asurface->ntriangles/3;
  if(ntriangles==0)return;
  if(surfacetype==1||surfacetype==-1){
    float rgbtemp[4];
    float *col;

    col = asurface->color;
    if(setbw==1){
      rgbtemp[0]=0.299*col[0]+0.587*col[1]+0.114*col[2];
      rgbtemp[1]=rgbtemp[0];
      rgbtemp[2]=rgbtemp[0];
    }
    else{
      rgbtemp[0]=col[0];
      rgbtemp[1]=col[1];
      rgbtemp[2]=col[2];
    }

    if(smooth_block_solid==0){
      rgbtemp[3]=asurface->color[3];
    }
    else{
      rgbtemp[3]=1.0;
    }
    if(data_type!=0){
      if(rgbtemp[3]<1.0&&trans_flag!=DRAW_TRANSPARENT)return;
      if(rgbtemp[3]>=1.0&&trans_flag==DRAW_TRANSPARENT)return;
    }
    if(
      trans_flag==DRAW_TRANSPARENT&&
      (
      (data_type==0&&use_transparency_data==1)||
      (data_type==1&&drawing_blockage_transparent==1)
      )
      ){
        if(rgbtemp[3]<0.99){
          drawing_transparent=1;
          drawing_blockage_transparent=1;
          transparenton_flag=1;
          transparenton();
        }
    }
    iso_specular[3] = 1.0;
    if(asurface->cullfaces==1)glDisable(GL_CULL_FACE);
    glPushAttrib(GL_LIGHTING_BIT);
    if(surfacetype==1){
      glEnable(GL_LIGHTING);
      glEnable(GL_COLOR_MATERIAL);
    }
    glBegin(GL_TRIANGLES);
    if(surfacetype==1){
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,asurface->color);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,iso_specular);
    }

    if(transparenton_flag==1){
	  glColor4fv(rgbtemp);
    }
    else{
	  glColor3fv(rgbtemp);
    }
    vertices_i=asurface->vertices;
    triangles_i=asurface->triangles;
    norm=asurface->norm;
    vertexnorm=asurface->vertexnorm;
    for(j=0;j<ntriangles;j++){
      i1=3*triangles_i[3*j];
      i2=3*triangles_i[3*j+1];
      i3=3*triangles_i[3*j+2];
      v1=vertices_i+i1;
      v2=vertices_i+i2;
      v3=vertices_i+i3;
      for(k=0;k<3;k++){
        vv1[k]=xyzmin[k]+xyzmaxdiff_local*v1[k]/65535.;
        vv2[k]=xyzmin[k]+xyzmaxdiff_local*v2[k]/65535.;
        vv3[k]=xyzmin[k]+xyzmaxdiff_local*v3[k]/65535.;
      }
		  if(smoothnorm==1){
        norm1 = vertexnorm+i1;
	      norm2 = vertexnorm+i2;
	      norm3 = vertexnorm+i3;
	      glNormal3sv(norm1);
        glVertex3fv(vv1);
	      glNormal3sv(norm2);
        glVertex3fv(vv2);
	      glNormal3sv(norm3);
        glVertex3fv(vv3);
      }
	    else{
        glNormal3sv(norm);
        glVertex3fv(vv1);
        glVertex3fv(vv2);
        glVertex3fv(vv3);
        norm += 3;
      }
    }
    glEnd();
    if(asurface->cullfaces==1)glEnable(GL_CULL_FACE);
    if(surfacetype==1){
      glDisable(GL_COLOR_MATERIAL);
      glDisable(GL_LIGHTING);
    }

    glPopAttrib();
    if(transparenton_flag==1)transparentoff();  
  }

  if(surfacetype==2){
    glPushMatrix();
    antialias(1);
    glLineWidth(line_width);
    glBegin(GL_LINES);
    glColor3fv(asurface->color);
    vertices_i=asurface->vertices;
    triangles_i=asurface->triangles;
    for(j=0;j<ntriangles;j++){
      i1=3*triangles_i[3*j];
      i2=3*triangles_i[3*j+1];
      i3=3*triangles_i[3*j+2];
      v1=vertices_i+i1;
      v2=vertices_i+i2;
      v3=vertices_i+i3;
      for(k=0;k<3;k++){
        vv1[k]=xyzmin[k]+xyzmaxdiff_local*v1[k]/65535.;
        vv2[k]=xyzmin[k]+xyzmaxdiff_local*v2[k]/65535.;
        vv3[k]=xyzmin[k]+xyzmaxdiff_local*v3[k]/65535.;
      }
      glVertex3fv(vv1);
      glVertex3fv(vv2);
      glVertex3fv(vv2);
      glVertex3fv(vv3);
      glVertex3fv(vv3);
      glVertex3fv(vv1);
    }
    glEnd();
    antialias(0);
    glPopMatrix();
  }

  if(surfacetype==3){
    glPushMatrix();
    antialias(1);
    glPointSize(plot3dpointsize);
    glBegin(GL_POINTS);
    glColor3fv(asurface->color);
    nvertices=asurface->nvertices;
    ntriangles=asurface->ntriangles/3;
    vertices_i=asurface->vertices;
    triangles_i=asurface->triangles;
    for(j=0;j<nvertices;j++){
      v1=vertices_i+3*j;
      for(k=0;k<3;k++){
        vv1[k]=xyzmin[k]+xyzmaxdiff_local*v1[k]/65535.;
      }

      glVertex3fv(vv1);
    }
    glEnd();
    antialias(0);
    glPopMatrix();
  }

  if(showisonormals==1){

    glPushMatrix();
    antialias(1);
    glLineWidth(line_width);
    glBegin(GL_LINES);
    glColor3f((float)1.,(float)1.,(float)1.);
    for(j=0;j<ntriangles;j++){
      i1=3*triangles_i[3*j];
      i2=3*triangles_i[3*j+1];
      i3=3*triangles_i[3*j+2];
      v1=vertices_i+i1;
      v2=vertices_i+i2;
      v3=vertices_i+i3;
      for(k=0;k<3;k++){
        vv1[k]=xyzmin[k]+xyzmaxdiff_local*v1[k]/65535.;
        vv2[k]=xyzmin[k]+xyzmaxdiff_local*v2[k]/65535.;
        vv3[k]=xyzmin[k]+xyzmaxdiff_local*v3[k]/65535.;
      }

	  	if(isonormtype==1){
        norm1 = vertexnorm+i1;
		    norm2 = vertexnorm+i2;
		    norm3 = vertexnorm+i3;
        for(k=0;k<3;k++){
          vv1n[k]=vv1[k]+norm1[k]/(8.*32768.)/4.0;
          vv2n[k]=vv2[k]+norm2[k]/(8.*32768.)/4.0;
          vv3n[k]=vv3[k]+norm3[k]/(8.*32768.)/4.0;
        }

        glVertex3fv(vv1);
        glVertex3fv(vv1n);
        glVertex3fv(vv2);
        glVertex3fv(vv2n);
        glVertex3fv(vv3);
        glVertex3fv(vv3n);
      }
		  else{
        for(k=0;k<3;k++){
          vv1n[k]=vv1[k]+norm[k]/(8.*32768.)/4.0;
          vv2n[k]=vv2[k]+norm[k]/(8.*32768.)/4.0;
          vv3n[k]=vv3[k]+norm[k]/(8.*32768.)/4.0;
        }

        glVertex3fv(vv1);
        glVertex3fv(vv1n);
        glVertex3fv(vv2);
        glVertex3fv(vv2n);
        glVertex3fv(vv3);
        glVertex3fv(vv3n);
        norm += 3;
      }
    }
    glEnd();
    antialias(0);
    glPopMatrix();
  }
}

/* ------------------ updateslicetypes ------------------------ */

void updateisotypes(void){
  int i;
  iso *isoi;

  nisotypes = 0;
  for(i=0;i<nisoinfo;i++){
    isoi = isoinfo+i;
    if(getisoindex(isoi)==-1)isotypes[nisotypes++]=i;
  }
  for(i=0;i<nisoinfo;i++){
    isoi = isoinfo+i;
    isoi->type=getisotype(isoi);
  }
}

/* ------------------ getisoindex ------------------------ */

int getisoindex(const iso *isoi){
  iso *isoi2;
  int j;

  for(j=0;j<nisotypes;j++){
    isoi2 = isoinfo+isotypes[j];
    if(strcmp(isoi->surface_label.longlabel,isoi2->surface_label.longlabel)==0)return isotypes[j];
  }
  return -1;
}

/* ------------------ getisotype ------------------------ */

int getisotype(const iso *isoi){
  iso *isoi2;
  int j;

  for(j=0;j<nisotypes;j++){
    isoi2 = isoinfo+isotypes[j];

    if(strcmp(isoi->surface_label.longlabel,isoi2->surface_label.longlabel)==0)return j;
  }
  return -1;
}

/* ------------------ getisottype ------------------------ */

int getisottype(const iso *isoi){
  iso *isoi2;
  int j;
  int jj;

  if(isoi->dataflag==0)return -1;
  jj = 0;
  for(j=0;j<nisoinfo;j++){
    isoi2 = isoinfo+j;

    if(isoi2->dataflag==0)continue;
    if(isoi2->firstshort==0)continue;
    if(strcmp(isoi->color_label.longlabel,isoi2->color_label.longlabel)==0)return jj;
    jj++;
  }
  return -1;
}

/* ------------------ update_isotype ------------------------ */

void update_isotype(void){
  int i;
  iso *isoi;


  for(i=0;i<nisoinfo;i++){
    isoi = isoinfo + i;
    if(isoi->loaded==0)continue;
    if(isoi->display==1&&isoi->type==iisotype)return;
  }

  for(i=0;i<nisoinfo;i++){
    isoi = isoinfo + i;
    if(isoi->loaded==0)continue;
    if(isoi->display==1){
      iisotype = getisoindex(isoi);
      return;
    }
  }

  iisotype = -1;
  return;
    
}

/* ------------------ isocompare ------------------------ */

int isocompare( const void *arg1, const void *arg2 ){
  iso *isoi, *isoj;

  isoi = isoinfo + *(int *)arg1;
  isoj = isoinfo + *(int *)arg2;

  if(strcmp(isoi->surface_label.longlabel,isoj->surface_label.longlabel)<0)return -1;
  if(strcmp(isoi->surface_label.longlabel,isoj->surface_label.longlabel)>0)return 1;
  if(isoi->blocknumber<isoj->blocknumber)return -1;
  if(isoi->blocknumber>isoj->blocknumber)return 1;
  return 0;
}

/* ------------------ updateisomenulabels ------------------------ */

void updateisomenulabels(void){
  int i;
  iso *isoi;
  char label[128];

  if(nisoinfo>0){
    FREEMEMORY(isoorderindex);
    NewMemory((void **)&isoorderindex,sizeof(int)*nisoinfo);
    for(i=0;i<nisoinfo;i++){
      isoorderindex[i]=i;
    }
    qsort( (int *)isoorderindex, (size_t)nisoinfo, sizeof(int), isocompare );

    for(i=0;i<nisoinfo;i++){
      isoi = isoinfo + i;
      STRCPY(isoi->menulabel,isoi->surface_label.longlabel);
      if(nmeshes>1){
	      mesh *isomesh;

		    isomesh = meshinfo + isoi->blocknumber;
        sprintf(label,"%s",isomesh->label);
        STRCAT(isoi->menulabel,", ");
        STRCAT(isoi->menulabel,label);
      }
      if(showfiles==1){
        STRCAT(isoi->menulabel,", ");
        STRCAT(isoi->menulabel,isoi->file);
      }
    } 
  }
}

/* ------------------ update_iso_showlevels ------------------------ */

void update_iso_showlevels(void){
  int nisolevels;
  int *showlevels;
  int i, j;
  mesh *meshi;

  if(loaded_isomesh==NULL)return;

  nisolevels=loaded_isomesh->nisolevels;
  showlevels=loaded_isomesh->showlevels;

  for(j=0;j<nmeshes;j++){
    meshi = meshinfo+j;
    if(meshi->isofilenum==-1)continue;
    for(i=0;i<nisolevels;i++){
      if(i<meshi->nisolevels)meshi->showlevels[i]=showlevels[i];
    }
  }
}

/* ------------------ uncompress_isodataframe ------------------------ */
#ifdef USE_ZLIB
void uncompress_isodataframe(isosurface *asurface_in, isosurface *asurface_out, int n, mesh *meshi){
  uLong countin;
  uLongf countout;
  unsigned char *compressed_data, *full_data;
  int ilevel, ivert, itri;
  float factor, offset[3];
  iso *ib;
  float *iso_colors;
  int n_iso_colors;
  int nverts, ntris;
  int *triangles_i;

  factor = (meshi->xyzmaxdiff/xyzmaxdiff)/65535.0;
  offset[0]=(meshi->xbar0-xbar0)/xyzmaxdiff;
  offset[1]=(meshi->ybar0-ybar0)/xyzmaxdiff;
  offset[2]=(meshi->zbar0-zbar0)/xyzmaxdiff;

  if(iso_ambient_ini==NULL||n_iso_ambient_ini==0){
    iso_colors=iso_ambient;
    n_iso_colors=n_iso_ambient;
  }
  else{
    iso_colors=iso_ambient_ini;
    n_iso_colors=n_iso_ambient_ini;
  }

  ib = isoinfo + meshi->isofilenum;

  for(ilevel=0;ilevel<n;ilevel++){
    unsigned short *verti;
    unsigned char *normi;
    int niso_vertices;
    int niso_triangles;
    
    compressed_data = asurface_in->comp_bufferframe;
    countin = asurface_in->ncomp_bufferframe;
    countout = asurface_in->nfull_bufferframe;
    full_data = asurface_in->full_bufferframe;

    uncompress(full_data,&countout,compressed_data,countin);

    niso_vertices = asurface_in->niso_vertices;
    niso_triangles = asurface_in->niso_triangles;
    if(niso_vertices==0)return;


    verti = (unsigned short *)full_data;
    for(ivert=0;ivert<niso_vertices;ivert++){
      isovert *isoverti;
      float *xyz,*vertnorm;
            
      isoverti = asurface_out->iso_vertices+ivert;
      xyz = isoverti->xyz;
      xyz[0]=offset[0]+factor*(*verti++); 
      xyz[1]=offset[1]+factor*(*verti++); 
      xyz[2]=offset[2]+factor*(*verti++); 
      isoverti->flag=0;
            
      vertnorm=isoverti->norm;
      vertnorm[0]=0.0;
      vertnorm[1]=0.0;
      vertnorm[2]=0.0;
      if(ilevel==0&&strcmp(ib->surface_label.shortlabel,"hrrpuv")==0){
        isoverti->color=hrrpuv_iso_color;
      }
      else{
        isoverti->color=iso_colors+4*ilevel;
      }
    }

    normi = (unsigned char *)full_data + 6*niso_vertices;
    for(ivert=0;ivert<niso_vertices;ivert++){
      isovert *isoverti;
      float *norm, *spherenorm;
      int index;

      isoverti = asurface_out->iso_vertices+ivert;
      norm=isoverti->norm;
      index = *normi++;
      spherenorm=sphereinfo->normals + 3*index;
      norm[0]=spherenorm[0];
      norm[1]=spherenorm[1];
      norm[2]=spherenorm[2];
    }

    NewMemory((void **)&triangles_i,3*niso_triangles*sizeof(int));
    if(niso_vertices<256){
      unsigned char *index;

      index = (unsigned char *)(full_data + 7*niso_vertices);
      for(itri=0;itri<niso_triangles;itri++){
        triangles_i[3*itri+0] = *index++;
        triangles_i[3*itri+1] = *index++;
        triangles_i[3*itri+2] = *index++;
      }
    }
    else if(niso_vertices>=256&&niso_vertices<65536){
      unsigned short *index;

      index = (unsigned short *)(full_data + 7*niso_vertices);
      for(itri=0;itri<niso_triangles;itri++){
        triangles_i[3*itri+0] = *index++;
        triangles_i[3*itri+1] = *index++;
        triangles_i[3*itri+2] = *index++;
      }
    }
    else if(niso_vertices>=65536){
      int *index;

      index = (int *)(full_data + 7*niso_vertices);
      for(itri=0;itri<niso_triangles;itri++){
        triangles_i[3*itri+0] = *index++;
        triangles_i[3*itri+1] = *index++;
        triangles_i[3*itri+2] = *index++;
      }
    }

    for(itri=0;itri<niso_triangles;itri++){
      isotri *isotrii;
              
      isotrii=asurface_out->iso_triangles+itri;
      isotrii->v1=asurface_out->iso_vertices+triangles_i[3*itri];
      isotrii->v2=asurface_out->iso_vertices+triangles_i[3*itri+1];
      isotrii->v3=asurface_out->iso_vertices+triangles_i[3*itri+2];
      isotrii->xyzmid[0]=(isotrii->v1->xyz[0]+isotrii->v2->xyz[0]+isotrii->v3->xyz[0])/3.0;
      isotrii->xyzmid[1]=(isotrii->v1->xyz[1]+isotrii->v2->xyz[1]+isotrii->v3->xyz[1])/3.0;
      isotrii->xyzmid[2]=(isotrii->v1->xyz[2]+isotrii->v2->xyz[2]+isotrii->v3->xyz[2])/3.0;
      isotrii->distance=-1.0;
      if(ilevel==0&&strcmp(ib->surface_label.shortlabel,"hrrpuv")==0){
        ib->colorlevels[ilevel]=hrrpuv_iso_color;
      }
      else{
        ib->colorlevels[ilevel]=iso_colors+4*ilevel;
      }
      isotrii->color=ib->colorlevels+ilevel;
      if(ib->dataflag==0){
        isotrii->v1->color=*isotrii->color;
        isotrii->v2->color=*isotrii->color;
        isotrii->v3->color=*isotrii->color;
      }
    }
    FREEMEMORY(triangles_i);
    if(countout!= asurface_in->nfull_bufferframe){
      printf("problems with decompressing iso data\n");
      printf(" %i %i\n",asurface_in->nfull_bufferframe,(int)countout);
    }
    asurface_in++;
    asurface_out++;
  }

}
#endif

/* ------------------ setisolabels ------------------------ */

void setisolabels(float smin, float smax, 
                    iso *sd, int *errorcode){
  char *scale;
  int isotype;
  databounds *sb;

  isotype=getisottype(sd);
  sb = isobounds + isotype;
  sb->label=&(sd->color_label);


  *errorcode=0;
  printf("setting up iso labels \n");
  scale=sb->scale;
  getIsoLabels(smin,smax,nrgb,
                sb->colorlabels,&scale,sb->levels256);
}

/* ------------------ sync_isobounds ------------------------ */

void sync_isobounds(int isottype){
  int i,j,ii,kk,ncount;
  isosurface *asurface;
  int firsttime=1;
  float tmin, tmax;

  // find number of iso-surfaces with values 

  ncount=0;
  for(i=0;i<nisoinfo;i++){
    iso *isoi;

    isoi = isoinfo + i;
    if(isoi->loaded==0||isoi->type!=iisotype||isoi->dataflag==0)continue;
    if(iisottype!=getisottype(isoi))continue;
    ncount++;
  }
  if(ncount<=1)return;

  // find min and max bounds for valued iso-surfaces

  for(i=0;i<nisoinfo;i++){
    iso *isoi;

    isoi = isoinfo + i;
    if(isoi->loaded==0||isoi->type!=iisotype||isoi->dataflag==0)continue;
    if(iisottype!=getisottype(isoi))continue;
    if(firsttime==1){
      firsttime=0;
      tmin=isoi->tmin;
      tmax=isoi->tmax;
    }
    else{
      if(tmin<isoi->tmin)isoi->tmin=tmin;
      if(tmax>isoi->tmax)isoi->tmax=tmax;
    }
  }

  // set min and max bounds for valued iso-surfaces

  for(i=0;i<nisoinfo;i++){
    iso *isoi;

    isoi = isoinfo + i;
    if(isoi->loaded==0||isoi->type!=iisotype||isoi->dataflag==0)continue;
    if(iisottype!=getisottype(isoi))continue;
    isoi->tmin=tmin;
    isoi->tmax=tmax;
  }

  // rescale all data

  for(i=0;i<nisoinfo;i++){
    iso *isoi;
    mesh *meshi;

    isoi = isoinfo + i;
    if(isoi->loaded==0||isoi->type!=iisotype||isoi->dataflag==0)continue;
    if(iisottype!=getisottype(isoi))continue;
    
    meshi = meshinfo + isoi->blocknumber;
    asurface=meshi->animatedsurfaces;

    for(ii=0;ii<meshi->nisosteps;ii++){
      for(j=0;j<meshi->nisolevels;j++){
        float tcolor, tcolor0, tcolorfactor;

        if(isoi->tmax>isoi->tmin){
          tcolor0 = (asurface->tmin-isoi->tmin)/(isoi->tmax-isoi->tmin);
          tcolorfactor = (asurface->tmax-asurface->tmin)/65535.;
          tcolorfactor /= (isoi->tmax-isoi->tmin);
        }
        else{
          tcolor0=0.5;
          tcolorfactor=0.0;
        }
        for(kk=0;kk<asurface->nvertices;kk++){
          tcolor = tcolor0 + asurface->tvertices[kk]*tcolorfactor;
          if(tcolor<0.0)tcolor=0.0;
          if(tcolor>1.0)tcolor=1.0;
          asurface->color8[kk] = (unsigned char)(tcolor*255);
        }
        asurface++;
      }
    }
  }
}

/* ------------------ compareisonodes ------------------------ */

int compare_iso_triangles( const void *arg1, const void *arg2 ){
  isotri *trii, *trij;

  trii = *(isotri **)arg1;
  trij = *(isotri **)arg2;

  if(trii->distance<trij->distance)return  1;
  if(trii->distance>trij->distance)return -1;
  return 0;
}

/* ------------------ sort_triangles ------------------------ */

void sort_iso_triangles(float *mm){
  int itri;
  int newflag;
  int dosort=0;

  if(niso_trans==0)return;
  newflag=1-iso_trans[0]->v1->flag;
  for(itri=0;itri<niso_trans;itri++){
    isotri *tri;
    float xyzeye[3];
    float *xyz;
    isovert *v1, *v2, *v3;

    tri = iso_trans[itri];
    v1 = tri->v1;
    v2 = tri->v2;
    v3 = tri->v3;
    if(v1->flag!=newflag){
      v1->flag=newflag;
      xyz = v1->xyz;
      xyzeye[0] = mm[0]*xyz[0] + mm[4]*xyz[1] +  mm[8]*xyz[2] + mm[12];
      xyzeye[1] = mm[1]*xyz[0] + mm[5]*xyz[1] +  mm[9]*xyz[2] + mm[13];
      xyzeye[2] = mm[2]*xyz[0] + mm[6]*xyz[1] + mm[10]*xyz[2] + mm[14];
      xyzeye[0]/=mscale[0];
      xyzeye[1]/=mscale[1];
      xyzeye[2]/=mscale[2];
      v1->distance=xyzeye[0]*xyzeye[0]+xyzeye[1]*xyzeye[1]+xyzeye[2]*xyzeye[2];
    }
    if(v2->flag!=newflag){
      v2->flag=newflag;
      xyz = v2->xyz;
      xyzeye[0] = mm[0]*xyz[0] + mm[4]*xyz[1] +  mm[8]*xyz[2] + mm[12];
      xyzeye[1] = mm[1]*xyz[0] + mm[5]*xyz[1] +  mm[9]*xyz[2] + mm[13];
      xyzeye[2] = mm[2]*xyz[0] + mm[6]*xyz[1] + mm[10]*xyz[2] + mm[14];
      xyzeye[0]/=mscale[0];
      xyzeye[1]/=mscale[1];
      xyzeye[2]/=mscale[2];
      v2->distance=xyzeye[0]*xyzeye[0]+xyzeye[1]*xyzeye[1]+xyzeye[2]*xyzeye[2];
    }
    if(v3->flag!=newflag){
      v3->flag=newflag;
      xyz = v3->xyz;
      xyzeye[0] = mm[0]*xyz[0] + mm[4]*xyz[1] +  mm[8]*xyz[2] + mm[12];
      xyzeye[1] = mm[1]*xyz[0] + mm[5]*xyz[1] +  mm[9]*xyz[2] + mm[13];
      xyzeye[2] = mm[2]*xyz[0] + mm[6]*xyz[1] + mm[10]*xyz[2] + mm[14];
      xyzeye[0]/=mscale[0];
      xyzeye[1]/=mscale[1];
      xyzeye[2]/=mscale[2];
      v3->distance=xyzeye[0]*xyzeye[0]+xyzeye[1]*xyzeye[1]+xyzeye[2]*xyzeye[2];
    }
    tri->distance=(v1->distance+v2->distance+v3->distance);
    if(itri>0&&dosort==0&&tri->distance>iso_trans[itri-1]->distance)dosort==1;
  }
  if(dosort==1)qsort((isotri **)iso_trans,(size_t)niso_trans,sizeof(isotri **),compare_iso_triangles);
}

/* ------------------ update_isotri_list ------------------------ */

void Update_Isotris(int flag){
  int ilev,itri;
  isosurface *asurfi;
  isotri **iso_trans_tmp,**iso_opaques_tmp;
  int *showlevels;
  mesh *meshi;
  float *colorptr;
  isosurface *asurface;
  int ntris;
  iso *loaded_iso;

  if(loaded_isomesh==NULL||loaded_isomesh->isofilenum==-1)return;
  loaded_iso=isoinfo + loaded_isomesh->isofilenum;
 
  if(iso_trans_list==NULL||iso_opaques_list==NULL){
    int iitime;

    niso_timesteps=loaded_isomesh->nisosteps;
    if(iso_trans_list==NULL){
      int i;

      NewMemory((void **)&niso_trans_list,niso_timesteps*sizeof(int));
      NewMemory((void **)&iso_trans_list,niso_timesteps*sizeof(isotri **));
      for(i=0;i<niso_timesteps;i++){
        iso_trans_list[i]=NULL;
      }
    }
    if(iso_opaques_list==NULL){
      int i;

      NewMemory((void **)&niso_opaques_list,niso_timesteps*sizeof(int));
      NewMemory((void **)&iso_opaques_list,niso_timesteps*sizeof(isotri **));
      for(i=0;i<niso_timesteps;i++){
        iso_opaques_list[i]=NULL;
      }
    }
    for(iitime=0;iitime<niso_timesteps;iitime++){
      int i;

      ntris=0;
      for(i=0;i<nisoinfo;i++){
        iso *isoi;
        int ilev;
    
        isoi = isoinfo+i;
        if(isoi->loaded==0||isoi->display==0)continue;

        meshi = meshinfo + isoi->blocknumber;
        asurface = meshi->animatedsurfaces + iitime*meshi->nisolevels;
        for(ilev=0;ilev<meshi->nisolevels;ilev++){
          asurfi = asurface + ilev;
          ntris+=asurfi->niso_triangles;
        }
      }
      if(ntris>0){
        NewMemory((void **)&iso_trans,ntris*sizeof(isotri *));
        iso_trans_list[iitime]=iso_trans;
        NewMemory((void **)&iso_opaques,ntris*sizeof(isotri *));
        iso_opaques_list[iitime]=iso_opaques;
      }
    }
    flag=1;
  }
  if(flag==1){
    int iitime;

    for(iitime=0;iitime<niso_timesteps;iitime++){
      niso_trans_list[iitime]=-1;
      niso_opaques_list[iitime]=-1;
    }
  }

  iso_trans=iso_trans_list[loaded_isomesh->iiso];
  iso_opaques=iso_opaques_list[loaded_isomesh->iiso];
  niso_trans=niso_trans_list[loaded_isomesh->iiso];
  niso_opaques=niso_opaques_list[loaded_isomesh->iiso];

  if(niso_trans==-1||niso_opaques==-1){
    int i;

    flag=1;
    iso_trans_tmp=iso_trans;
    iso_opaques_tmp=iso_opaques;
    niso_trans=0;
    niso_opaques=0;
    for(i=0;i<nisoinfo;i++){
      iso *isoi;
    
      isoi = isoinfo+i;
      if(isoi->loaded==0||isoi->display==0)continue;
  
      CheckMemory;
      meshi = meshinfo + isoi->blocknumber;
      asurface = meshi->animatedsurfaces + meshi->iiso*meshi->nisolevels;
      showlevels=meshi->showlevels;
  
      if(transparent_state==ALL_TRANSPARENT){
        for(ilev=0;ilev<meshi->nisolevels;ilev++){
          if(showlevels[ilev]==0)continue;
          asurfi = asurface + ilev;
          if(asurfi->iso_triangles>0){
            niso_trans += asurfi->niso_triangles;
            for(itri=0;itri<asurfi->niso_triangles;itri++){
              *iso_trans_tmp++=asurfi->iso_triangles+itri;
            }
            colorptr=isoi->colorlevels[ilev];
            colorptr[3]=transparentlevel;
          }
        }
      }
      else if(transparent_state==MIN_SOLID){
        for(ilev=0;ilev<1;ilev++){
          if(showlevels[ilev]==0)continue;
          asurfi = asurface + ilev;
          if(asurfi->niso_triangles>0){
            niso_opaques += asurfi->niso_triangles;
            for(itri=0;itri<asurfi->niso_triangles;itri++){
              *iso_opaques_tmp++=asurfi->iso_triangles+itri;
            }
            colorptr=isoi->colorlevels[ilev];
            colorptr[3]=1.0;
          }
        }
        for(ilev=1;ilev<meshi->nisolevels;ilev++){
          if(showlevels[ilev]==0)continue;
          asurfi = asurface + ilev;
          if(asurfi->niso_triangles>0){
            niso_trans += asurfi->niso_triangles;
            for(itri=0;itri<asurfi->niso_triangles;itri++){
              *iso_trans_tmp++=asurfi->iso_triangles+itri;
            }
            colorptr=isoi->colorlevels[ilev];
            colorptr[3]=transparentlevel;
          }
        }
      }
      else if(transparent_state==MAX_SOLID){
        for(ilev=0;ilev<meshi->nisolevels-1;ilev++){
          if(showlevels[ilev]==0)continue;
          asurfi = asurface + ilev;
          if(asurfi->niso_triangles>0){
            niso_trans += asurfi->niso_triangles;
            for(itri=0;itri<asurfi->niso_triangles;itri++){
              *iso_trans_tmp++=asurfi->iso_triangles+itri;
            } 
            colorptr=isoi->colorlevels[ilev];
            colorptr[3]=transparentlevel;
          }
        }
        for(ilev=meshi->nisolevels-1;ilev<meshi->nisolevels;ilev++){
          if(showlevels[ilev]==0)continue;
          asurfi = asurface + ilev;
          if(asurfi->niso_triangles>0){
            niso_opaques += asurfi->niso_triangles;
            for(itri=0;itri<asurfi->niso_triangles;itri++){
              *iso_opaques_tmp++=asurfi->iso_triangles+itri;
            }
            colorptr=isoi->colorlevels[ilev];
            colorptr[3]=1.0;
          }
        }
      }
      else if(transparent_state==ALL_SOLID){
        for(ilev=0;ilev<meshi->nisolevels;ilev++){
          CheckMemory;
          if(showlevels[ilev]==0)continue;
          asurfi = asurface + ilev;
          if(asurfi->niso_triangles>0){
            niso_opaques += asurfi->niso_triangles;
            for(itri=0;itri<asurfi->niso_triangles;itri++){
              *iso_opaques_tmp++=asurfi->iso_triangles+itri;
            }
            colorptr=isoi->colorlevels[ilev];
            colorptr[3]=1.0;
          }
        }
      }
    }
  }

  if(sort_transparency==1&&niso_trans>0){
    sort_iso_triangles(modelview_scratch);
  }
  niso_trans_list[loaded_isomesh->iiso]=niso_trans;
  niso_opaques_list[loaded_isomesh->iiso]=niso_opaques;

  CheckMemory;
}

/* ------------------ get_loaded_isomesh ------------------------ */

mesh *get_loaded_isomesh(void){
  mesh *return_mesh;
  int i,nsteps=-1;

  if(isoinfo==NULL)return NULL;
  return_mesh=NULL;
  for(i=0;i<nisoinfo;i++){
    mesh *mesh2;
    iso *isoi;

    isoi = isoinfo + i;
    if(isoi->loaded==0)continue;
    mesh2 = meshinfo + isoi->blocknumber;
    if(nsteps==-1||mesh2->nisosteps<nsteps){
      return_mesh = mesh2;
      nsteps=mesh2->nisosteps;
    }
  }
  return return_mesh;
}
