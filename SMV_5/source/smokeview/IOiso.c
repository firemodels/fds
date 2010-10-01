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
#include "egz_stdio.h"
#include "MALLOC.h"
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

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

/* ------------------ getcisolevels ------------------------ */

void getcisolevels(const char *isofile, float **levelsptr, int *nisolevels){
  EGZ_FILE *isostream;
  char buffer[1024];
  int one;

  *nisolevels=0;
#ifdef EGZ
  isostream=EGZ_FOPEN(isofile,"rb",0,2);
#else
  isostream=EGZ_FOPEN(isofile,"rb");
#endif
  if(isostream==NULL)return;
  EGZ_FREAD(&one,4,1,isostream);
  EGZ_FREAD(nisolevels,4,1,isostream);
  FREEMEMORY(*levelsptr);
  NewMemory((void **)levelsptr,*nisolevels*sizeof(float));
  EGZ_FREAD(*levelsptr,4,(unsigned int)(*nisolevels),isostream);
  EGZ_FCLOSE(isostream);
}

/* ------------------ getisolevels ------------------------ */

void getisolevels(const char *isofile, int dataflag, float **levelsptr, int *nisolevels){
  int one;
  int version;
  int len[3],labellengths=0;
  int nlevels;
  EGZ_FILE *isostreamptr;

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
  int i,j,n,ii,jj;
  float out[3];
  isosurface *asurface;
  int nisopoints, nisotriangles;
  int maxfullbuffer,maxcompbuffer;
  unsigned char *comp_buffer;
  unsigned char *full_buffer;

  float time, time_max;
  EGZ_FILE *isostream;
  float *xyznorm=NULL;

  unsigned short *vertices_i=NULL, *tvertices_i=NULL;
  int *triangles_i; 
  unsigned char *triangles1_i; 
  unsigned short *triangles2_i;
  unsigned char *color8=NULL;
  int i1, i2, i3;
  int ntriangles_i, nvertices_i;
  unsigned short *v1, *v2, *v3;
  short *norm,*vertexnorm;
  float area;
  float isomin, isomax;
  int blocknumber;
  int error;
  int skip;
  float tcolor, tcolor0, tcolorfactor;

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
  meshi = meshinfo+blocknumber;
  unloadiso(meshi);
  FREEMEMORY(meshi->isotimes);
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

      loaded_isomesh=NULL;
      for(i=0;i<niso_files;i++){
        isoi=isoinfo+i;
        if(isoi->loaded==1){
          loaded_isomesh=meshinfo+isoi->blocknumber;
          break;
        }
      }
      update_iso_showlevels();
    }
    return;
  }
  meshi->isofilenum=ifile;
  highlight_mesh = blocknumber;
   
  if(ib->compression_type==1){
    int nbuffer;

    getisoheader(ib->comp_file,&isostream,ib->size_file,&nisopoints, &nisotriangles,&nbuffer,&maxfullbuffer, &maxcompbuffer,
      &meshi->isolevels, &meshi->nisolevels, &meshi->nisosteps, isoframestep, &ib->normaltable, &ib->nnormaltable);
    getfile_size(ib->comp_file,&file_size);
    if(nbuffer>0){
      NewMemory((void **)&ib->comp_buffer,nbuffer);
      maxfullbuffer=1.01*maxfullbuffer+600;
      maxcompbuffer=1.01*maxcompbuffer+600;
      NewMemory((void **)&ib->full_bufferframe,maxfullbuffer);
      NewMemory((void **)&ib->comp_bufferframe,maxcompbuffer);
    }
  }
  else{
    getisosizes(file, ib->dataflag, &isostream, 
      &nisopoints, &nisotriangles, 
      &meshi->isolevels, &meshi->nisolevels, &meshi->nisosteps, isoframestep, 
      &ib->tmin, &ib->tmax,
      endian_data);

    getfile_size(file,&file_size);
  }
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
  for(i=0;i<meshi->nisolevels;i++){
    meshi->showlevels[i]=1;
  }
  isomin=meshi->isolevels[0];
  isomax=meshi->isolevels[0];
  meshi->isomin_index=0;
  meshi->isomax_index=0;
  for(i=1;i<meshi->nisolevels;i++){
    if(meshi->isolevels[i]<isomin){
      isomin=meshi->isolevels[i];
      meshi->isomin_index=i;
    }
    if(meshi->isolevels[i]>isomax){
      isomax=meshi->isolevels[i];
      meshi->isomax_index=i;
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

  comp_buffer=ib->comp_buffer;
  asurface=meshi->animatedsurfaces;
  jj=0;
  i=0;
  time_max = -1000000.0;
  local_starttime = glutGet(GLUT_ELAPSED_TIME);
  for(;;){
    int skip_frame;

    skip_frame=0;
    jj++;
  // 1
  // nisolevels
  // levels(0), ..., levels(nisolevels-1)
  // npoints
  // xyz(0), ..., xyz(npoints-1)

  // time
  // nvertices
  // ntriangles
  // nbuffer
  // buffer(0), ..., buffer(nbuffer-1)

  // ***** 
    if(ib->compression_type==1){
      int nvertices, ntriangles, nbuffer, nfull;
      int ibuffer[4];

      full_buffer=ib->full_bufferframe;
      EGZ_FREAD(&time,4,1,isostream);
      if(EGZ_FEOF(isostream)!=0)break;
      meshi->isotimes[i]=time;
      printf("isosurface time=%.2f\n",time);
      for(j=0;j<meshi->nisolevels;j++){
        EGZ_FREAD(ibuffer,4,4,isostream);
        nvertices=ibuffer[0];
        ntriangles=ibuffer[1];
        nfull=ibuffer[2];
        nbuffer=ibuffer[3];
        asurface->compression_type=ib->compression_type;
        asurface->comp_bufferframe=NULL;
        asurface->full_bufferframe=NULL;
        asurface->ncomp_bufferframe=0;
        asurface->nfull_bufferframe=0;
        if(nbuffer>0){
          asurface->comp_bufferframe=comp_buffer;
          asurface->ncomp_bufferframe=nbuffer;

          asurface->full_bufferframe=full_buffer;
          asurface->nfull_bufferframe=nfull;

          EGZ_FREAD(comp_buffer,1,nbuffer,isostream);

          comp_buffer+=nbuffer;
          full_buffer+=nfull;
        }
        if(settmin_i==1&&time<tmin_i)continue;
        if(settmax_i==1&&time>tmax_i)continue;


        asurface->dataflag=ib->dataflag;
	      asurface->nvertices=nvertices;
	      asurface->ntriangles=ntriangles;
	      asurface->triangles=NULL;
        asurface->triangles1=NULL;
        asurface->triangles2=NULL;
	      asurface->vertices=NULL;
        asurface->s_norm=NULL;
        
        asurface->vertices=(unsigned short *)asurface->full_bufferframe;
        asurface->s_norm=(unsigned char *)(asurface->full_bufferframe+6*nvertices);
        if(nvertices>0&&nvertices<256){
          asurface->triangles1=(unsigned char *)(asurface->full_bufferframe+7*nvertices);
        }
        else if(nvertices>=256&&nvertices<65536){
          asurface->triangles2=(unsigned short *)(asurface->full_bufferframe+7*nvertices);
        }
        else if(nvertices>=65536){
          asurface->triangles=(int *)(asurface->full_bufferframe+7*nvertices);
        }
        asurface->tvertices=NULL;
        asurface->color8=NULL;
        asurface->norm=NULL;
	      asurface->vertexnorm=NULL;
        asurface++;
      }
    }
    else{
      EGZ_FREAD(&time,4,1,isostream);
      if(EGZ_FEOF(isostream)!=0)break;
      skip_frame=1;
      if(time>time_max){
        skip_frame=0;
        time_max=time;
      }
      meshi->isotimes[i]=time;
      if(jj%isoframestep!=0||(settmin_i==1&&time<tmin_i)||(settmax_i==1&&time>tmax_i)||skip_frame==1){
      }
      else{
        printf("isosurface time=%f\n",time);
      }
      for(j=0;j<meshi->nisolevels;j++){
        asurface->dataflag=ib->dataflag;
        EGZ_FREAD(&nvertices_i,4,1,isostream);
        if(EGZ_FEOF(isostream)!=0)break;
        EGZ_FREAD(&ntriangles_i,4,1,isostream);
        if(EGZ_FEOF(isostream)!=0)break;
        if(jj%isoframestep!=0||(settmin_i==1&&time<tmin_i)||(settmax_i==1&&time>tmax_i)||skip_frame==1){
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
        triangles_i=NULL;
        triangles1_i=NULL;
        triangles2_i=NULL;
	      vertices_i=NULL;
        tvertices_i=NULL;
        color8=NULL;
        vertexnorm=NULL;
        norm=NULL;
        FREEMEMORY(xyznorm);
       
#define FREELOCAL_ISO FREEMEMORY(xyznorm);FREEMEMORY(vertexnorm);FREEMEMORY(vertices_i);FREEMEMORY(tvertices_i);\
                      FREEMEMORY(triangles1_i);FREEMEMORY(triangles2_i);FREEMEMORY(triangles_i);\
                      FREEMEMORY(norm);FREEMEMORY(color8)

        if(nvertices_i>0){
          if( NewMemory((void **)&xyznorm,3*nvertices_i*sizeof(float))==0 ){
            FREELOCAL_ISO;
            break;
          }
          if(NewMemory((void **)&vertexnorm,3*nvertices_i*sizeof(short))==0){
            FREELOCAL_ISO;
            break;
          }
          if(NewMemory((void **)&vertices_i,3*nvertices_i*sizeof(unsigned short))==0){
            FREELOCAL_ISO;
            break;
          }
          if(ib->dataflag==1){
            if(NewMemory((void **)&tvertices_i,nvertices_i*sizeof(unsigned short))==0){
              FREELOCAL_ISO;
              break;
            }
            if(NewMemory((void **)&color8,nvertices_i*sizeof(unsigned char))==0){
              FREELOCAL_ISO;
              break;
            }
          }
        }
        for(ii=0;ii<3*nvertices_i;ii++){
          xyznorm[ii]=0.0;
        }
        if(nvertices_i>0){
          EGZ_FREAD(vertices_i,2,(unsigned int)(3*nvertices_i),isostream);
        }
        if(ib->dataflag==1&&nvertices_i>0){
          EGZ_FREAD(&asurface->tmin,4,1,isostream);
          EGZ_FREAD(&asurface->tmax,4,1,isostream);
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
          for(ii=0;ii<nvertices_i;ii++){
            tcolor = tcolor0 + tvertices_i[ii]*tcolorfactor;
            if(tcolor<0.0)tcolor=0.0;
            if(tcolor>1.0)tcolor=1.0;
            color8[ii] = (unsigned char)(tcolor*255);
          }
        }
	      if(EGZ_FEOF(isostream)!=0)break;
        if(nvertices_i<256&&nvertices_i>0){
          if(ntriangles_i>0){
            if(NewMemory((void **)&triangles1_i,ntriangles_i*sizeof(unsigned char))==0){
              FREELOCAL_ISO;
              break;
            }
  	    	  EGZ_FREAD(triangles1_i,1,(unsigned int)ntriangles_i,isostream);
          }
        }
        else if(nvertices_i>=256&&nvertices_i<65536){
          if(ntriangles_i>0){
            if(NewMemory((void **)&triangles2_i,ntriangles_i*sizeof(unsigned short))==0){
              FREELOCAL_ISO;
              break;
            }
  		      EGZ_FREAD(triangles2_i,2,(unsigned int)ntriangles_i,isostream);
          }
        }
        else{
          if(ntriangles_i>0){
            if(NewMemory((void **)&triangles_i,ntriangles_i*sizeof(int))==0){
              FREELOCAL_ISO;
              break;
            }
  		      EGZ_FREAD(triangles_i,4,(unsigned int)ntriangles_i,isostream);
          } 
        }
	      if(EGZ_FEOF(isostream)!=0)break;

        if(ntriangles_i>0){
          if(NewMemory((void **)&norm,ntriangles_i*sizeof(short))==0){
            FREELOCAL_ISO;
            break;
          }
        }

        for(n=0;n<ntriangles_i/3;n++){
          if(nvertices_i<256){
            i1=3*triangles1_i[3*n];
            i2=3*triangles1_i[3*n+1];
            i3=3*triangles1_i[3*n+2];
          }
          else if(nvertices_i>=256&&nvertices_i<65536){
            i1=3*triangles2_i[3*n];
            i2=3*triangles2_i[3*n+1];
            i3=3*triangles2_i[3*n+2];
          }
          else{
            i1=3*triangles_i[3*n];
            i2=3*triangles_i[3*n+1];
            i3=3*triangles_i[3*n+2];
          }
          v1=vertices_i+i1;
          v2=vertices_i+i2;
          v3=vertices_i+i3;
          calcNormal2(v1,v2,v3,out,&area);
          norm[3*n  ]=(short)(out[0]*32767);
          norm[3*n+1]=(short)(out[1]*32767);
          norm[3*n+2]=(short)(out[2]*32767);
    	    xyznorm[i1  ] += out[0]*area;
		      xyznorm[i1+1] += out[1]*area;
		      xyznorm[i1+2] += out[2]*area;
		      xyznorm[i2  ] += out[0]*area;
		      xyznorm[i2+1] += out[1]*area;
		      xyznorm[i2+2] += out[2]*area;
		      xyznorm[i3  ] += out[0]*area;
		      xyznorm[i3+1] += out[1]*area;
		      xyznorm[i3+2] += out[2]*area;
        }
	      for(n=0;n<nvertices_i;n++){
		      ReduceToUnit(xyznorm+3*n);
          vertexnorm[3*n  ]=(short)(xyznorm[3*n  ]*32767);
          vertexnorm[3*n+1]=(short)(xyznorm[3*n+1]*32767);
          vertexnorm[3*n+2]=(short)(xyznorm[3*n+2]*32767);
        }
	      FREEMEMORY(xyznorm);
	      asurface->nvertices=nvertices_i;
	      asurface->ntriangles=ntriangles_i;
	      asurface->triangles=triangles_i;
        asurface->triangles1=triangles1_i;
        asurface->triangles2=triangles2_i;
	      asurface->vertices=vertices_i;
        if(ib->dataflag==1){
          asurface->tvertices=tvertices_i;
          asurface->color8=color8;
        }
        asurface->norm=norm;
	      asurface->vertexnorm=vertexnorm;

        asurface++;
      }
    }
    if(ib->compression_type==1){
      i++;
      if(i>=meshi->nisosteps)break;
      continue;
    }
    if(skip_frame==1||jj%isoframestep!=0||(settmin_i==1&&time<tmin_i)||(settmax_i==1&&time>tmax_i)){
     // if(skip_frame==1)jj--;
    }
    else{
      i++;
      if(i>=meshi->nisosteps)break;
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
  loaded_isomesh=meshinfo+ib->blocknumber;
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
  GLUTPOSTREDISPLAY
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
    asurface=meshi->animatedsurfaces;
    CheckMemoryOff;
    if(ib->compression_type==0){
      for(n=0;n<meshi->nisosteps;n++){
        for(j=0;j<meshi->nisolevels;j++){
          FREEMEMORY(asurface->triangles);
          FREEMEMORY(asurface->triangles1);
          FREEMEMORY(asurface->triangles2);
          FREEMEMORY(asurface->vertices);
          FREEMEMORY(asurface->norm);
          FREEMEMORY(asurface->vertexnorm);
          asurface++;
        }
      }
    }
    CheckMemoryOn;
    FREEMEMORY(meshi->isotimes);
    FREEMEMORY(meshi->animatedsurfaces);
    FREEMEMORY(meshi->showlevels);
  }
  meshi->nisosteps=0;
  FREEMEMORY(ib->comp_bufferframe);
  FREEMEMORY(ib->full_bufferframe);
  FREEMEMORY(ib->comp_buffer);
  FREEMEMORY(ib->normaltable);
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

void drawiso(const mesh *meshi,int tranflag){
  int i, j,k;
  float vv1[3],vv2[3],vv3[3];
  float vv1n[3],vv2n[3],vv3n[3];
  isosurface *asurface;
  short *norm;
  unsigned short *v1, *v2, *v3, tval1=0, tval2=0, tval3=0;
  unsigned short *vertices_i=NULL,*tvertices_i=NULL;
  int *triangles_i;
  unsigned short *triangles2_i;
  unsigned char *triangles1_i;
  unsigned char *color8;
  int nvertices;
  int i1, i2, i3;
  short *norm1,*norm2,*norm3,*vertexnorm;
  float *iso_colors;
  int n_iso_colors;
  int icolor;
  int ntriangles;
  int *showlevels, nisolevels;
  int isomin_index,isomax_index;
  float factor, offset[3];
  iso *isoi=NULL;
  float iso_color_tmp[4];
  float *iso_color_ptr;
  int iso_lighting;

  if(meshi->isofilenum>=0){
    isoi = isoinfo + meshi->isofilenum;
  }
  if(isoi->dataflag==1){
    iso_lighting=0;
  }
  else{
    iso_lighting=1;
  }

  showlevels=meshi->showlevels;
  nisolevels=meshi->nisolevels;
  isomin_index=meshi->isomin_index;
  isomax_index=meshi->isomax_index;
  factor = (meshi->xyzmaxdiff/xyzmaxdiff)/65535.0;
  offset[0]=(meshi->xbar0-xbar0)/xyzmaxdiff;
  offset[1]=(meshi->ybar0-ybar0)/xyzmaxdiff;
  offset[2]=(meshi->zbar0-zbar0)/xyzmaxdiff;


  if(tranflag==DRAW_TRANSPARENT&&visAIso!=1)return;

  if(iso_ambient_ini==NULL||n_iso_ambient_ini==0){
    iso_colors=iso_ambient;
    n_iso_colors=n_iso_ambient;
  }
  else{
    iso_colors=iso_ambient_ini;
    n_iso_colors=n_iso_ambient_ini;
  }

  if(visAIso==1){

    asurface = meshi->animatedsurfaces + meshi->iiso*meshi->nisolevels - 1;
    if(cullfaces==1)glDisable(GL_CULL_FACE);

    iso_specular[3] = 1.0;
    if(tranflag==DRAW_TRANSPARENT)transparenton();

    glPushAttrib(GL_LIGHTING_BIT);
    if(iso_lighting==1){
      glEnable(GL_LIGHTING);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,iso_specular);
      glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,iso_shininess);
    }
    glBegin(GL_TRIANGLES);

    for(i=0;i<nisolevels;i++){
      asurface++;
      icolor=i;
      if(icolor>n_iso_colors-1)icolor=n_iso_colors-1;
      if(showlevels[i]==0)continue;
      if(tranflag==DRAW_TRANSPARENT){
        if(transparent_state==ALL_TRANSPARENT){
        }
        else if(transparent_state==MIN_SOLID){
          if(i==isomin_index)continue;
        }
        else if(transparent_state==MAX_SOLID){
          if(i==isomax_index)continue;
        }
        else if(transparent_state==ALL_SOLID){
          continue;
        }
      }
      else if(tranflag==DRAW_SOLID){
        if(transparent_state==ALL_TRANSPARENT){
          continue;
        }
        else if(transparent_state==MIN_SOLID){
          if(i!=isomin_index)continue;
        }
        else if(transparent_state==MAX_SOLID){
          if(i!=isomax_index)continue;
        }
        else if(transparent_state==ALL_SOLID){
        }
      }
      nvertices=asurface->nvertices;
      ntriangles=asurface->ntriangles/3;
      if(ntriangles==0)continue;
      if(i==0&&strcmp(isoi->surface_label.shortlabel,"hrrpuv")==0){
        iso_color_ptr=hrrpuv_iso_color;
      }
      else{
        iso_color_ptr=iso_colors+4*icolor;
      }
      if(setbw!=0){
        float greylevel;

        greylevel=color2bw(iso_color_ptr);
        iso_color_tmp[0]=greylevel;
        iso_color_tmp[1]=greylevel;
        iso_color_tmp[2]=greylevel;
        iso_color_ptr=iso_color_tmp;
      }
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,iso_color_ptr);
      vertices_i=asurface->vertices;
      if(asurface->dataflag==1){
        tvertices_i=asurface->tvertices;
        color8 = asurface->color8;
      }
      triangles_i=asurface->triangles;
      triangles1_i=asurface->triangles1;
      triangles2_i=asurface->triangles2;
      norm=asurface->norm;
  	  vertexnorm=asurface->vertexnorm;
      for(j=0;j<ntriangles;j++){
        if(nvertices<256){
          i1=3*triangles1_i[3*j];
          i2=3*triangles1_i[3*j+1];
          i3=3*triangles1_i[3*j+2];
        }
        else if(nvertices>=256&&nvertices<65536){
          i1=3*triangles2_i[3*j];
          i2=3*triangles2_i[3*j+1];
          i3=3*triangles2_i[3*j+2];
        }
        else{
          i1=3*triangles_i[3*j];
          i2=3*triangles_i[3*j+1];
          i3=3*triangles_i[3*j+2];
        }
        v1=vertices_i+i1;
        v2=vertices_i+i2;
        v3=vertices_i+i3;
        for(k=0;k<3;k++){
          vv1[k]=offset[k]+factor*v1[k];
          vv2[k]=offset[k]+factor*v2[k];
          vv3[k]=offset[k]+factor*v3[k];
        }
        if(asurface->dataflag==1){
          tval1=color8[i1/3];
          tval2=color8[i2/3];
          tval3=color8[i3/3];
        }
	    	if(isonormtype==1
          ||(isonormtype==0&&isoi->compression_type==1)
          ){
          if(isoi->compression_type==1){
            norm1 = isoi->normaltable + 3*asurface->s_norm[i1/3];
            norm2 = isoi->normaltable + 3*asurface->s_norm[i2/3];
            norm3 = isoi->normaltable + 3*asurface->s_norm[i3/3];
          }
          else{
            norm1 = vertexnorm+i1;
	  	      norm2 = vertexnorm+i2;
		        norm3 = vertexnorm+i3;
          }
          if(asurface->dataflag==1){
            glNormal3sv(norm1);
            glColor4fv(rgb_full[tval1]);
            glVertex3fv(vv1);

  		      glNormal3sv(norm2);
            glColor4fv(rgb_full[tval2]);
            glVertex3fv(vv2);

   		      glNormal3sv(norm3);
            glColor4fv(rgb_full[tval3]);
            glVertex3fv(vv3);
          }
          else{
            glNormal3sv(norm1);
            glVertex3fv(vv1);

  		      glNormal3sv(norm2);
            glVertex3fv(vv2);

   		      glNormal3sv(norm3);
            glVertex3fv(vv3);
          }
        }
		    else{
          if(asurface->dataflag==1){
            glNormal3sv(norm);
            glColor4fv(rgb_full[tval1]);
            glVertex3fv(vv1);

            glColor4fv(rgb_full[tval2]);
            glVertex3fv(vv2);

            glColor4fv(rgb_full[tval3]);
            glVertex3fv(vv3);
          }
          else{
            glNormal3sv(norm);
            glVertex3fv(vv1);
            glVertex3fv(vv2);
            glVertex3fv(vv3);
          }
          norm += 3;
        }
      }
    }
    glEnd();
    if(asurface->dataflag==1)glDisable(GL_COLOR_MATERIAL);

    glPopAttrib();

    if(tranflag==DRAW_TRANSPARENT)transparentoff();
    if(cullfaces==1)glEnable(GL_CULL_FACE);
  }

  if(visAIso==2){
    asurface = meshi->animatedsurfaces + meshi->iiso*meshi->nisolevels-1;

    glPushAttrib(GL_LIGHTING_BIT);
    antialias(1);
    glLineWidth(isolinewidth);
    glBegin(GL_LINES);
    for(i=0;i<nisolevels;i++){
      asurface++;
      if(asurface->dataflag==1)tvertices_i=asurface->tvertices;
      if(showlevels[i]==0)continue;
      icolor=i;
      if(icolor>n_iso_colors-1)icolor=n_iso_colors-1;
      if(setbw==0){
        iso_color_ptr = iso_colors+4*icolor;
      }
      else{
        float greylevel;

        greylevel=color2bw(iso_colors+4*icolor);
        iso_color_tmp[0]=greylevel;
        iso_color_tmp[1]=greylevel;
        iso_color_tmp[2]=greylevel;
        iso_color_ptr=iso_color_tmp;
      }
      glColor3fv(iso_color_ptr);
      nvertices=asurface->nvertices;
      ntriangles=asurface->ntriangles/3;
      vertices_i=asurface->vertices;
      triangles_i=asurface->triangles;
      triangles1_i=asurface->triangles1;
      triangles2_i=asurface->triangles2;
      for(j=0;j<ntriangles;j++){
        if(nvertices<256){
          i1=3*triangles1_i[3*j];
          i2=3*triangles1_i[3*j+1];
          i3=3*triangles1_i[3*j+2];
        }
        else if(nvertices>=256&&nvertices<65536){
          i1=3*triangles2_i[3*j];
          i2=3*triangles2_i[3*j+1];
          i3=3*triangles2_i[3*j+2];
        }
        else{
          i1=3*triangles_i[3*j];
          i2=3*triangles_i[3*j+1];
          i3=3*triangles_i[3*j+2];
        }
        v1=vertices_i+i1;
        v2=vertices_i+i2;
        v3=vertices_i+i3;
        for(k=0;k<3;k++){
          vv1[k]=offset[k]+factor*v1[k];
          vv2[k]=offset[k]+factor*v2[k];
          vv3[k]=offset[k]+factor*v3[k];
        }
        if(asurface->dataflag==1){
          tval1=tvertices_i[i1/3]/256;
          if(tval1>255)tval1=255;
          tval2=tvertices_i[i2/3]/256;
          if(tval2>255)tval2=255;
          tval3=tvertices_i[i3/3]/256;
          if(tval3>255)tval3=255;
        }
        if(asurface->dataflag==1){
          glColor4fv(rgb_full[tval1]);
          glVertex3fv(vv1);
          glColor4fv(rgb_full[tval2]);
          glVertex3fv(vv2);
          glVertex3fv(vv2);
          glColor4fv(rgb_full[tval3]);
          glVertex3fv(vv3);
          glVertex3fv(vv3);
          glColor4fv(rgb_full[tval1]);
          glVertex3fv(vv1);
        }
        else{
          glVertex3fv(vv1);
          glVertex3fv(vv2);
          glVertex3fv(vv2);
          glVertex3fv(vv3);
          glVertex3fv(vv3);
          glVertex3fv(vv1);
        }
      }
    }
    glEnd();
    antialias(0);
    glPopAttrib();

  }

  if(showisonormals==1){
    asurface = meshi->animatedsurfaces + meshi->iiso*meshi->nisolevels - 1;

    antialias(1);
    glLineWidth(isolinewidth);
    glBegin(GL_LINES);
    for(i=0;i<nisolevels;i++){
      asurface++;
      if(showlevels[i]==0)continue;
      glColor3f((float)1.,(float)1.,(float)1.);
      nvertices=asurface->nvertices;
      ntriangles=asurface->ntriangles/3;
      vertices_i=asurface->vertices;
      if(asurface->dataflag==1)tvertices_i=asurface->tvertices;
      triangles_i=asurface->triangles;
      triangles1_i=asurface->triangles1;
      triangles2_i=asurface->triangles2;
      norm=asurface->norm;
  	  vertexnorm=asurface->vertexnorm;
      for(j=0;j<ntriangles;j++){
        if(nvertices<256){
          i1=3*triangles1_i[3*j];
          i2=3*triangles1_i[3*j+1];
          i3=3*triangles1_i[3*j+2];
        }
        else if(nvertices>=256&&nvertices<65536){
          i1=3*triangles2_i[3*j];
          i2=3*triangles2_i[3*j+1];
          i3=3*triangles2_i[3*j+2];
        }
        else{
          i1=3*triangles_i[3*j];
          i2=3*triangles_i[3*j+1];
          i3=3*triangles_i[3*j+2];
        }
        v1=vertices_i+i1;
        v2=vertices_i+i2;
        v3=vertices_i+i3;
        for(k=0;k<3;k++){
          vv1[k]=offset[k]+factor*v1[k];
          vv2[k]=offset[k]+factor*v2[k];
          vv3[k]=offset[k]+factor*v3[k];
        }
	    	if(isonormtype==1
          ||(isonormtype==0&&isoi->compression_type==1)
          ){
          if(isoi->compression_type==1){
            norm1 = isoi->normaltable + 3*asurface->s_norm[i1/3];
            norm2 = isoi->normaltable + 3*asurface->s_norm[i2/3];
            norm3 = isoi->normaltable + 3*asurface->s_norm[i3/3];
          }
          else{
            norm1 = vertexnorm+i1;
	  	      norm2 = vertexnorm+i2;
		        norm3 = vertexnorm+i3;
          }
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
    }
    glEnd();
    antialias(0);
  }

  if(visAIso==3){
    asurface = meshi->animatedsurfaces + meshi->iiso*meshi->nisolevels - 1;

    antialias(1);
    glPointSize(isopointsize);
    glBegin(GL_POINTS);
    for(i=0;i<nisolevels;i++){
      asurface++;
      if(showlevels[i]==0)continue;
      icolor=i;
      if(icolor>n_iso_colors-1)icolor=n_iso_colors-1;
      if(setbw==0){
        iso_color_ptr = iso_colors+4*icolor;
      }
      else{
        float greylevel;

        greylevel=color2bw(iso_colors+4*icolor);
        iso_color_tmp[0]=greylevel;
        iso_color_tmp[1]=greylevel;
        iso_color_tmp[2]=greylevel;
        iso_color_ptr=iso_color_tmp;
      }
      glColor3fv(iso_color_ptr);
      nvertices=asurface->nvertices;
      ntriangles=asurface->ntriangles/3;
      vertices_i=asurface->vertices;
      if(asurface->dataflag==1)tvertices_i=asurface->tvertices;
      triangles_i=asurface->triangles;
      triangles1_i=asurface->triangles1;
      triangles2_i=asurface->triangles2;
      for(j=0;j<nvertices;j++){
        v1=vertices_i+3*j;
        for(k=0;k<3;k++){
          vv1[k]=offset[k]+factor*v1[k];
        }
        if(asurface->dataflag==1){
          tval1=tvertices_i[j]/256;
          if(tval1>255)tval1=255;
          glColor4fv(rgb_full[tval1]);
          glVertex3fv(vv1);
        }
        else{
          glVertex3fv(vv1);
        }
      }
      asurface++;
    }
    glEnd();
    antialias(0);
  }
}

/* ------------------ drawtiso ------------------------ */

void drawtiso(const mesh *meshi,int tranflag){
  int i, j,k;
  float vv1[3],vv2[3],vv3[3];
  float vv1n[3],vv2n[3],vv3n[3];
  isosurface *asurface;
  short *norm;
  unsigned short *v1, *v2, *v3, tval1=0, tval2=0, tval3=0;
  float ttval1, ttval2, ttval3;
  unsigned short *vertices_i=NULL,*tvertices_i=NULL;
  int *triangles_i;
  unsigned short *triangles2_i;
  unsigned char *triangles1_i;
  unsigned char *color8;
  int nvertices;
  int i1, i2, i3;
  short *norm1,*norm2,*norm3,*vertexnorm;
  float *iso_colors;
  int n_iso_colors;
  int icolor;
  int ntriangles;
  int *showlevels, nisolevels;
  int isomin_index,isomax_index;
  float factor, offset[3];
  iso *isoi=NULL;
  float iso_color_tmp[4];
  float *iso_color_ptr;
  int iso_lighting;

  if(meshi->isofilenum>=0){
    isoi = isoinfo + meshi->isofilenum;
  }
  if(isoi->dataflag==1){
    iso_lighting=0;
  }
  else{
    iso_lighting=1;
  }

  showlevels=meshi->showlevels;
  nisolevels=meshi->nisolevels;
  isomin_index=meshi->isomin_index;
  isomax_index=meshi->isomax_index;
  factor = (meshi->xyzmaxdiff/xyzmaxdiff)/65535.0;
  offset[0]=(meshi->xbar0-xbar0)/xyzmaxdiff;
  offset[1]=(meshi->ybar0-ybar0)/xyzmaxdiff;
  offset[2]=(meshi->zbar0-zbar0)/xyzmaxdiff;


  if(tranflag==DRAW_TRANSPARENT&&visAIso!=1)return;

  if(iso_ambient_ini==NULL||n_iso_ambient_ini==0){
    iso_colors=iso_ambient;
    n_iso_colors=n_iso_ambient;
  }
  else{
    iso_colors=iso_ambient_ini;
    n_iso_colors=n_iso_ambient_ini;
  }

  if(visAIso==1){


    asurface = meshi->animatedsurfaces + meshi->iiso*meshi->nisolevels - 1;
    if(cullfaces==1)glDisable(GL_CULL_FACE);

    iso_specular[3] = 1.0;
    if(tranflag==DRAW_TRANSPARENT)transparenton();
    glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
    glEnable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D,texture_iso_colorbar_id);

    glPushAttrib(GL_LIGHTING_BIT);
    if(iso_lighting==1){
      glEnable(GL_LIGHTING);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,iso_specular);
      glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,iso_shininess);
    }
    glBegin(GL_TRIANGLES);

    for(i=0;i<nisolevels;i++){
      asurface++;
      icolor=i;
      if(icolor>n_iso_colors-1)icolor=n_iso_colors-1;
      if(showlevels[i]==0)continue;
      if(tranflag==DRAW_TRANSPARENT){
        if(transparent_state==ALL_TRANSPARENT){
        }
        else if(transparent_state==MIN_SOLID){
          if(i==isomin_index)continue;
        }
        else if(transparent_state==MAX_SOLID){
          if(i==isomax_index)continue;
        }
        else if(transparent_state==ALL_SOLID){
          continue;
        }
      }
      else if(tranflag==DRAW_SOLID){
        if(transparent_state==ALL_TRANSPARENT){
          continue;
        }
        else if(transparent_state==MIN_SOLID){
          if(i!=isomin_index)continue;
        }
        else if(transparent_state==MAX_SOLID){
          if(i!=isomax_index)continue;
        }
        else if(transparent_state==ALL_SOLID){
        }
      }
      nvertices=asurface->nvertices;
      ntriangles=asurface->ntriangles/3;
      if(ntriangles==0)continue;
      if(setbw==0){
        iso_color_ptr = iso_colors+4*icolor;
      }
      else{
        float greylevel;

        greylevel=color2bw(iso_colors+4*icolor);
        iso_color_tmp[0]=greylevel;
        iso_color_tmp[1]=greylevel;
        iso_color_tmp[2]=greylevel;
        iso_color_ptr=iso_color_tmp;
      }
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,iso_color_ptr);
      vertices_i=asurface->vertices;
      if(asurface->dataflag==1){
        tvertices_i=asurface->tvertices;
        color8 = asurface->color8;
      }
      triangles_i=asurface->triangles;
      triangles1_i=asurface->triangles1;
      triangles2_i=asurface->triangles2;
      norm=asurface->norm;
  	  vertexnorm=asurface->vertexnorm;
      for(j=0;j<ntriangles;j++){
        if(nvertices<256){
          i1=3*triangles1_i[3*j];
          i2=3*triangles1_i[3*j+1];
          i3=3*triangles1_i[3*j+2];
        }
        else if(nvertices>=256&&nvertices<65536){
          i1=3*triangles2_i[3*j];
          i2=3*triangles2_i[3*j+1];
          i3=3*triangles2_i[3*j+2];
        }
        else{
          i1=3*triangles_i[3*j];
          i2=3*triangles_i[3*j+1];
          i3=3*triangles_i[3*j+2];
        }
        v1=vertices_i+i1;
        v2=vertices_i+i2;
        v3=vertices_i+i3;
        for(k=0;k<3;k++){
          vv1[k]=offset[k]+factor*v1[k];
          vv2[k]=offset[k]+factor*v2[k];
          vv3[k]=offset[k]+factor*v3[k];
        }
        if(asurface->dataflag==1){
          tval1=color8[i1/3];
          tval2=color8[i2/3];
          tval3=color8[i3/3];
          ttval1=tval1/255.0;
          ttval2=tval2/255.0;
          ttval3=tval3/255.0;
        }
	    	if(isonormtype==1
          ||(isonormtype==0&&isoi->compression_type==1)
          ){
          if(isoi->compression_type==1){
            norm1 = isoi->normaltable + 3*asurface->s_norm[i1/3];
            norm2 = isoi->normaltable + 3*asurface->s_norm[i2/3];
            norm3 = isoi->normaltable + 3*asurface->s_norm[i3/3];
          }
          else{
            norm1 = vertexnorm+i1;
	  	      norm2 = vertexnorm+i2;
		        norm3 = vertexnorm+i3;
          }
          if(asurface->dataflag==1){
            glNormal3sv(norm1);
            //glColor4fv(rgb_full[tval1]);
            glTexCoord1f(ttval1);
            glVertex3fv(vv1);

  		      glNormal3sv(norm2);
            //glColor4fv(rgb_full[tval2]);
            glTexCoord1f(ttval2);
            glVertex3fv(vv2);

   		      glNormal3sv(norm3);
            //glColor4fv(rgb_full[tval3]);
            glTexCoord1f(ttval3);
            glVertex3fv(vv3);
          }
          else{
            glNormal3sv(norm1);
            glVertex3fv(vv1);

  		      glNormal3sv(norm2);
            glVertex3fv(vv2);

   		      glNormal3sv(norm3);
            glVertex3fv(vv3);
          }
        }
		    else{
          if(asurface->dataflag==1){
            glNormal3sv(norm);
            //glColor4fv(rgb_full[tval1]);
            glTexCoord1f(ttval1);
            glVertex3fv(vv1);

            //glColor4fv(rgb_full[tval2]);
            glTexCoord1f(ttval2);
            glVertex3fv(vv2);

            //glColor4fv(rgb_full[tval3]);
            glTexCoord1f(ttval3);
            glVertex3fv(vv3);
          }
          else{
            glNormal3sv(norm);
            glVertex3fv(vv1);
            glVertex3fv(vv2);
            glVertex3fv(vv3);
          }
          norm += 3;
        }
      }
    }
    glEnd();
    if(asurface->dataflag==1)glDisable(GL_COLOR_MATERIAL);

    glPopAttrib();

    if(tranflag==DRAW_TRANSPARENT)transparentoff();
    glDisable(GL_TEXTURE_1D);
    if(cullfaces==1)glEnable(GL_CULL_FACE);
  }

  if(visAIso==2){
    asurface = meshi->animatedsurfaces + meshi->iiso*meshi->nisolevels-1;

    glPushAttrib(GL_LIGHTING_BIT);
    antialias(1);
    glLineWidth(isolinewidth);
    glBegin(GL_LINES);
    for(i=0;i<nisolevels;i++){
      asurface++;
      if(asurface->dataflag==1)tvertices_i=asurface->tvertices;
      if(showlevels[i]==0)continue;
      icolor=i;
      if(icolor>n_iso_colors-1)icolor=n_iso_colors-1;
      if(setbw==0){
        iso_color_ptr = iso_colors+4*icolor;
      }
      else{
        float greylevel;

        greylevel=color2bw(iso_colors+4*icolor);
        iso_color_tmp[0]=greylevel;
        iso_color_tmp[1]=greylevel;
        iso_color_tmp[2]=greylevel;
        iso_color_ptr=iso_color_tmp;
      }
      glColor3fv(iso_color_ptr);
      nvertices=asurface->nvertices;
      ntriangles=asurface->ntriangles/3;
      vertices_i=asurface->vertices;
      triangles_i=asurface->triangles;
      triangles1_i=asurface->triangles1;
      triangles2_i=asurface->triangles2;
      for(j=0;j<ntriangles;j++){
        if(nvertices<256){
          i1=3*triangles1_i[3*j];
          i2=3*triangles1_i[3*j+1];
          i3=3*triangles1_i[3*j+2];
        }
        else if(nvertices>=256&&nvertices<65536){
          i1=3*triangles2_i[3*j];
          i2=3*triangles2_i[3*j+1];
          i3=3*triangles2_i[3*j+2];
        }
        else{
          i1=3*triangles_i[3*j];
          i2=3*triangles_i[3*j+1];
          i3=3*triangles_i[3*j+2];
        }
        v1=vertices_i+i1;
        v2=vertices_i+i2;
        v3=vertices_i+i3;
        for(k=0;k<3;k++){
          vv1[k]=offset[k]+factor*v1[k];
          vv2[k]=offset[k]+factor*v2[k];
          vv3[k]=offset[k]+factor*v3[k];
        }
        if(asurface->dataflag==1){
          tval1=tvertices_i[i1/3]/256;
          if(tval1>255)tval1=255;
          tval2=tvertices_i[i2/3]/256;
          if(tval2>255)tval2=255;
          tval3=tvertices_i[i3/3]/256;
          if(tval3>255)tval3=255;
        }
        if(asurface->dataflag==1){
          glColor4fv(rgb_full[tval1]);
          glVertex3fv(vv1);
          glColor4fv(rgb_full[tval2]);
          glVertex3fv(vv2);
          glVertex3fv(vv2);
          glColor4fv(rgb_full[tval3]);
          glVertex3fv(vv3);
          glVertex3fv(vv3);
          glColor4fv(rgb_full[tval1]);
          glVertex3fv(vv1);
        }
        else{
          glVertex3fv(vv1);
          glVertex3fv(vv2);
          glVertex3fv(vv2);
          glVertex3fv(vv3);
          glVertex3fv(vv3);
          glVertex3fv(vv1);
        }
      }
    }
    glEnd();
    antialias(0);
    glPopAttrib();

  }

  if(showisonormals==1){
    asurface = meshi->animatedsurfaces + meshi->iiso*meshi->nisolevels - 1;

    antialias(1);
    glLineWidth(isolinewidth);
    glBegin(GL_LINES);
    for(i=0;i<nisolevels;i++){
      asurface++;
      if(showlevels[i]==0)continue;
      glColor3f((float)1.,(float)1.,(float)1.);
      nvertices=asurface->nvertices;
      ntriangles=asurface->ntriangles/3;
      vertices_i=asurface->vertices;
      if(asurface->dataflag==1)tvertices_i=asurface->tvertices;
      triangles_i=asurface->triangles;
      triangles1_i=asurface->triangles1;
      triangles2_i=asurface->triangles2;
      norm=asurface->norm;
  	  vertexnorm=asurface->vertexnorm;
      for(j=0;j<ntriangles;j++){
        if(nvertices<256){
          i1=3*triangles1_i[3*j];
          i2=3*triangles1_i[3*j+1];
          i3=3*triangles1_i[3*j+2];
        }
        else if(nvertices>=256&&nvertices<65536){
          i1=3*triangles2_i[3*j];
          i2=3*triangles2_i[3*j+1];
          i3=3*triangles2_i[3*j+2];
        }
        else{
          i1=3*triangles_i[3*j];
          i2=3*triangles_i[3*j+1];
          i3=3*triangles_i[3*j+2];
        }
        v1=vertices_i+i1;
        v2=vertices_i+i2;
        v3=vertices_i+i3;
        for(k=0;k<3;k++){
          vv1[k]=offset[k]+factor*v1[k];
          vv2[k]=offset[k]+factor*v2[k];
          vv3[k]=offset[k]+factor*v3[k];
        }
	    	if(isonormtype==1
          ||(isonormtype==0&&isoi->compression_type==1)
          ){
          if(isoi->compression_type==1){
            norm1 = isoi->normaltable + 3*asurface->s_norm[i1/3];
            norm2 = isoi->normaltable + 3*asurface->s_norm[i2/3];
            norm3 = isoi->normaltable + 3*asurface->s_norm[i3/3];
          }
          else{
            norm1 = vertexnorm+i1;
	  	      norm2 = vertexnorm+i2;
		        norm3 = vertexnorm+i3;
          }
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
    }
    glEnd();
    antialias(0);
  }

  if(visAIso==3){
    asurface = meshi->animatedsurfaces + meshi->iiso*meshi->nisolevels - 1;

    antialias(1);
    glPointSize(plot3dpointsize);
    glBegin(GL_POINTS);
    for(i=0;i<nisolevels;i++){
      asurface++;
      if(showlevels[i]==0)continue;
      icolor=i;
      if(icolor>n_iso_colors-1)icolor=n_iso_colors-1;
      if(setbw==0){
        iso_color_ptr = iso_colors+4*icolor;
      }
      else{
        float greylevel;

        greylevel=color2bw(iso_colors+4*icolor);
        iso_color_tmp[0]=greylevel;
        iso_color_tmp[1]=greylevel;
        iso_color_tmp[2]=greylevel;
        iso_color_ptr=iso_color_tmp;
      }
      glColor3fv(iso_color_ptr);
      nvertices=asurface->nvertices;
      ntriangles=asurface->ntriangles/3;
      vertices_i=asurface->vertices;
      if(asurface->dataflag==1)tvertices_i=asurface->tvertices;
      triangles_i=asurface->triangles;
      triangles1_i=asurface->triangles1;
      triangles2_i=asurface->triangles2;
      for(j=0;j<nvertices;j++){
        v1=vertices_i+3*j;
        for(k=0;k<3;k++){
          vv1[k]=offset[k]+factor*v1[k];
        }
        if(asurface->dataflag==1){
          tval1=tvertices_i[j]/256;
          if(tval1>255)tval1=255;
          glColor4fv(rgb_full[tval1]);
          glVertex3fv(vv1);
        }
        else{
          glVertex3fv(vv1);
        }
      }
      asurface++;
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
  static int jjj=0;
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
      (data_type==0&&transparentflag==1)||
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
  for(i=0;i<niso_files;i++){
    isoi = isoinfo+i;
    if(getisoindex(isoi)==-1)isotypes[nisotypes++]=i;
  }
  for(i=0;i<niso_files;i++){
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
  for(j=0;j<niso_files;j++){
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


  for(i=0;i<niso_files;i++){
    isoi = isoinfo + i;
    if(isoi->loaded==0)continue;
    if(isoi->display==1&&isoi->type==iisotype)return;
  }

  for(i=0;i<niso_files;i++){
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

  if(niso_files>0){
    FREEMEMORY(isoorderindex);
    NewMemory((void **)&isoorderindex,sizeof(int)*niso_files);
    for(i=0;i<niso_files;i++){
      isoorderindex[i]=i;
    }
    qsort( (int *)isoorderindex, (size_t)niso_files, sizeof(int), isocompare );

    for(i=0;i<niso_files;i++){
      isoi = isoinfo + i;
      STRCPY(isoi->menulabel,isoi->surface_label.longlabel);
      if(nmeshes>1){
        sprintf(label,"Mesh %i",1+isoi->blocknumber);
        STRCAT(isoi->menulabel,", ");
        STRCAT(isoi->menulabel,label);
      }
      if(showfiles==1){
        STRCAT(isoi->menulabel,", ");
        STRCAT(isoi->menulabel,isoi->file);
      }
      if(isoi->compression_type==1){
        STRCAT(isoi->menulabel," (ZLIB)");
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
void uncompress_isodataframe(isosurface *asurface, int n){
  uLong countin;
  uLongf countout;
  unsigned char *compressed_data, *full_data;
  int i;

  for(i=0;i<n;i++){
    asurface++;
    compressed_data = asurface->comp_bufferframe;
    countin = asurface->ncomp_bufferframe;
    countout = asurface->nfull_bufferframe;
    full_data = asurface->full_bufferframe;
    uncompress(full_data,&countout,compressed_data,countin);
    if(countout!= asurface->nfull_bufferframe){
      printf("problems with decompressing iso data\n");
      printf(" %i %i\n",asurface->nfull_bufferframe,(int)countout);
    }
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
  for(i=0;i<niso_files;i++){
    iso *isoi;

    isoi = isoinfo + i;
    if(isoi->loaded==0||isoi->type!=iisotype||isoi->dataflag==0)continue;
    if(iisottype!=getisottype(isoi))continue;
    ncount++;
  }
  if(ncount<=1)return;

  // find min and max bounds for valued iso-surfaces

  for(i=0;i<niso_files;i++){
    iso *isoi;
    mesh *meshi;

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

  for(i=0;i<niso_files;i++){
    iso *isoi;
    mesh *meshi;

    isoi = isoinfo + i;
    if(isoi->loaded==0||isoi->type!=iisotype||isoi->dataflag==0)continue;
    if(iisottype!=getisottype(isoi))continue;
    isoi->tmin=tmin;
    isoi->tmax=tmax;
  }

  // rescale all data

  for(i=0;i<niso_files;i++){
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
