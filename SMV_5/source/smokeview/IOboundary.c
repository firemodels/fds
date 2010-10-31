// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>  
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include "flowfiles.h"
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "MALLOC.h"
#include "egz_stdio.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"
#include "smv_endian.h"
#include "smokeviewdefs.h"

// svn revision character string
char IOboundary_revision[]="$Revision$";



int getpatchfacedir(mesh *gb, int i1, int i2, int j1, int j2, int k1, int k2, 
                    int *blockonpatch, mesh **meshonpatch);
int getpatchface2dir(mesh *gb, int i1, int i2, int j1, int j2, int k1, int k2, int patchdir,
                    int *blockonpatch, mesh **meshonpatch);
int getpatchindex(const patch *patchi);

/* ------------------ readpatch ------------------------ */

void readpatch(int ifile, int flag, int *errorcode){
  int error;
  FILE_SIZE lenfile;
  int patchfilenum;
  float *xyzpatchcopy;
  float *xyzpatch_ignitecopy;
  int *patchblankcopy;
  char patchlonglabel[31], patchshortlabel[31], patchunit[31];
  int mxpatch_frames;
  int i1, i2, j1, j2, k1, k2, n;
  int i, j, k;
  int headersize, framesize;
  int statfile;
  STRUCTSTAT statbuffer;
  int nbb;
  int ibartemp,jbartemp,kbartemp;
  float *xplttemp,*yplttemp,*zplttemp;
  int blocknumber;
  patch *pi,*patchbase,*patchi;
  mesh *meshi;
  float patchmin_global, patchmax_global;
  int local_first,nsize,iblock;
  int ext_wall;

  blockagedata *bc;
  int nn;
  int filenum;
  char *patchscale;
  int npatchloaded;
  int ncompressed_buffer;
  char *file;
  int local_starttime=0, local_stoptime=0;
  FILE_SIZE file_size=0;
  int local_starttime0=0, local_stoptime0=0;  
  float delta_time, delta_time0;
  int file_unit;

  if(patchinfo[ifile].loaded==0&&flag==UNLOAD)return;

  local_starttime0 = glutGet(GLUT_ELAPSED_TIME);
  local_first=1;
  CheckMemory;
  patchfilenum=ifile;
  pi = patchinfo + ifile;
  file = pi->file;
  blocknumber = pi->blocknumber;
  highlight_mesh = blocknumber;
  meshi = meshinfo+blocknumber;
  update_current_mesh(meshi);
  filenum = meshi->patchfilenum;
  if(filenum>=0&&filenum<npatch_files){
    patchinfo[filenum].loaded=0;
    patchinfo[filenum].display=0;
  }

  meshi->patchfilenum=ifile;
  pi->display=0;
  plotstate=getplotstate(DYNAMIC_PLOTS);

  nbb = meshi->nbptrs;
  if(nbb==0)nbb=1;
  updatefaces=1;
  *errorcode=0;
  FREEMEMORY(meshi->blockonpatch);
  FREEMEMORY(meshi->meshonpatch);
  FREEMEMORY(meshi->patchdir);
  FREEMEMORY(meshi->patch_surfindex);
  FREEMEMORY(meshi->pi1);
  FREEMEMORY(meshi->pi2);
  FREEMEMORY(meshi->pj1);
  FREEMEMORY(meshi->pj2);
  FREEMEMORY(meshi->pk1);
  FREEMEMORY(meshi->pk2);
  FREEMEMORY(meshi->patchrow);
  FREEMEMORY(meshi->patchcol);
  FREEMEMORY(meshi->blockstart);
  FREEMEMORY(meshi->zipoffset);
  FREEMEMORY(meshi->zipsize);
  FREEMEMORY(meshi->patchtype);
  FREEMEMORY(meshi->visPatches);
  FREEMEMORY(meshi->xyzpatch);
  FREEMEMORY(meshi->xyzpatch_threshold);
  FREEMEMORY(meshi->ipqq);
  FREEMEMORY(meshi->ipqq_zlib);
  FREEMEMORY(meshi->ipqqi_zlib);
  FREEMEMORY(meshi->pqq);
  FREEMEMORY(meshi->thresholdtime);
  FREEMEMORY(meshi->patchtimes);
  FREEMEMORY(meshi->patchblank);

  if(meshi->patch_contours!=NULL){  
    ASSERT(meshi->npatches>0&&meshi->mxpatch_frames>0);
    for(i=0;i<meshi->npatches*meshi->mxpatch_frames;i++){
      if(meshi->patch_contours[i]!=NULL){
        freecontour(meshi->patch_contours[i]);
      }
    }
    FREEMEMORY(meshi->patch_contours);
  }

  if(flag==UNLOAD){
    update_patchtype();
    update_unit_defs();
    updatetimes();
    meshi->npatches=0;
    {
      int enableflag=1;

      for(i=0;i<npatch_files;i++){
        patch *patchii;

        patchii = patchinfo + i;
        if(patchii->loaded==1&&patchii->compression_type==1){
          enableflag=0;
          break;
        }
      }
      if(enableflag==1)enable_boundary_glui();
    }
    updatemenu=1;
#ifdef _DEBUG
    printf("After boundary file unload: ");
    PrintMemoryInfo;
#endif
    return;
  }
  if(ifile>=0&&ifile<npatch_files){
    global2localpatchbounds(patchinfo[ifile].label.shortlabel);
  }

  if(colorlabelpatch!=NULL){
    for(n=0;n<MAXRGB;n++){
      FREEMEMORY(colorlabelpatch[n]);
    }
    FREEMEMORY(colorlabelpatch);
  }
  patchinfo[ifile].extreme_max=0;
  patchinfo[ifile].extreme_min=0;
  ibartemp=meshi->ibar;
  jbartemp=meshi->jbar;
  kbartemp=meshi->kbar;
  xplttemp=meshi->xplt;
  yplttemp=meshi->yplt;
  zplttemp=meshi->zplt;
  do_threshold=0;
  
  if(activate_threshold==1){
    if(
      strncmp(patchinfo[ifile].label.shortlabel,"TEMP",4) == 0||
      strncmp(patchinfo[ifile].label.shortlabel,"temp",4) == 0
      ){
      do_threshold=1;
    }
  }

  update_patch_hist(pi);

  lenfile = strlen(file);
  file_unit=15;
  FORTget_file_unit(&file_unit,&file_unit);
  if(patchinfo[ifile].compression_type==0){
    FILE_SIZE labellen=LABELLEN;

    FORTgetpatchsizes1(&file_unit,file,patchlonglabel,patchshortlabel,patchunit,&endian,&meshi->npatches,&headersize,&error,
                       lenfile,labellen,labellen,labellen);
    if(error!=0){
      readpatch(ifile,UNLOAD,&error);
      *errorcode=1;
      return;
    }
  }
  else{
    meshi->npatches=0;
    getpatchheader(file,&meshi->npatches,&patchmin,&patchmax);
  }
  if(meshi->npatches>0){
    if(
       NewMemory((void **)&meshi->meshonpatch,sizeof(mesh *)*meshi->npatches)==0||
       NewMemory((void **)&meshi->blockonpatch,sizeof(int)*meshi->npatches)==0||
       NewMemory((void **)&meshi->patchdir    ,sizeof(int)*meshi->npatches)==0||
       NewMemory((void **)&meshi->patch_surfindex   ,sizeof(int)*meshi->npatches)==0||
       NewMemory((void **)&meshi->pi1         ,sizeof(int)*meshi->npatches)==0||
       NewMemory((void **)&meshi->pi2         ,sizeof(int)*meshi->npatches)==0||
       NewMemory((void **)&meshi->pj1         ,sizeof(int)*meshi->npatches)==0||
       NewMemory((void **)&meshi->pj2         ,sizeof(int)*meshi->npatches)==0||
       NewMemory((void **)&meshi->pk1         ,sizeof(int)*meshi->npatches)==0||
       NewMemory((void **)&meshi->pk2         ,sizeof(int)*meshi->npatches)==0||
       NewMemory((void **)&meshi->patchtype   ,sizeof(int)*meshi->npatches)==0||
       NewMemory((void **)&meshi->visPatches  ,sizeof(int)*meshi->npatches)==0||
       NewMemory((void **)&meshi->patchrow    ,sizeof(int)*meshi->npatches)==0||
       NewMemory((void **)&meshi->patchcol    ,sizeof(int)*meshi->npatches)==0||
       NewMemory((void **)&meshi->blockstart  ,sizeof(int)*(1+meshi->npatches))==0){
      *errorcode=1;
      if(patchinfo[ifile].compression_type==0){
        FORTclosefortranfile(&file_unit);
      }
      readpatch(ifile,UNLOAD,&error);
      return;
    }
  }

  if(patchinfo[ifile].compression_type==0){
    FORTgetpatchsizes2(&file_unit,&patchinfo[ifile].version,
      &meshi->npatches,&meshi->npatchsize,
      meshi->pi1,meshi->pi2,meshi->pj1,meshi->pj2,meshi->pk1,meshi->pk2,meshi->patchdir,
      &headersize,&framesize);

    // loadpatchbysteps
    //  0 - load entire uncompressed data set
    //  1 - load uncompressed data set one frame at a time
    //  2 - load compressed data set

    loadpatchbysteps=0;
    if(flag==LOAD){
      mxpatch_frames = mxframes+51;
      statfile=STAT(file,&statbuffer);
      if(statfile==0&&framesize!=0){
        {
          int file_frames;
          file_frames=(statbuffer.st_size-headersize)/framesize+51;
          if(file_frames<mxpatch_frames)mxpatch_frames=file_frames;
        }
      }
      meshi->mxpatch_frames=mxpatch_frames;


  /*
  If the min and max boundary file values are specified then we don't have
  to read in the whole file to determine the bounds.  In this case, memory is allocated
  one time step at a time rather than for all time steps.
  */

      if(statfile==0&&setpatchmin==1&&setpatchmax==1)loadpatchbysteps=1;
    }
  }
  else{
    int nnsize=0;

    getpatchheader2(file,
      &patchinfo[ifile].version,
      meshi->pi1,meshi->pi2,
      meshi->pj1,meshi->pj2,
      meshi->pk1,meshi->pk2,
      meshi->patchdir);
    for(i=0;i<meshi->npatches;i++){
      int ii1, ii2, jj1, jj2, kk1, kk2;

      ii1=meshi->pi1[i];
      ii2=meshi->pi2[i];
      jj1=meshi->pj1[i];
      jj2=meshi->pj2[i];
      kk1=meshi->pk1[i];
      kk2=meshi->pk2[i];
      nnsize += (ii2+1-ii1)*(jj2+1-jj1)*(kk2+1-kk1);
    }
    meshi->npatchsize=nnsize;
    loadpatchbysteps=2;
  }

  if(meshi->npatchsize>0){
    if(
       NewMemory((void **)&meshi->xyzpatch,3*sizeof(float)*meshi->npatchsize)==0||
       NewMemory((void **)&meshi->xyzpatch_threshold,3*sizeof(float)*meshi->npatchsize)==0||
       NewMemory((void **)&meshi->thresholdtime,sizeof(float)*meshi->npatchsize)==0||
       NewMemory((void **)&meshi->patchblank,meshi->npatchsize*sizeof(int))==0
       ){
      *errorcode=1;
      patchinfo[ifile].loaded=0;
      patchinfo[ifile].display=0;
      if(patchinfo[ifile].compression_type==0){
        FORTclosefortranfile(&file_unit);
      }
      readpatch(ifile,UNLOAD,&error);
      return;
    }
  }
  for(n=0;n<meshi->npatchsize;n++){
    meshi->patchblank[n]=1;
  }
  xyzpatchcopy = meshi->xyzpatch;
  xyzpatch_ignitecopy = meshi->xyzpatch_threshold;
  patchblankcopy = meshi->patchblank;
  meshi->blockstart[0]=0;
  for(n=0;n<meshi->nbptrs;n++){
    bc=meshi->blockageinfoptrs[n];
    for(j=0;j<6;j++){
      bc->patchvis[j]=1;
    }
  }

  for(n=0;n<meshi->npatches;n++){
    float dxx, dyy, dzz, ig_factor;
    float dx_factor, dy_factor, dz_factor;

    i1=meshi->pi1[n]; i2=meshi->pi2[n];
    j1=meshi->pj1[n]; j2=meshi->pj2[n];
    k1=meshi->pk1[n]; k2=meshi->pk2[n];
    if(pi->version==0){
      meshi->patchdir[n]=getpatchfacedir(meshi,i1,i2,j1,j2,k1,k2,
        meshi->blockonpatch+n,meshi->meshonpatch+n);
    }
    else{
      int patchdir;

      patchdir=meshi->patchdir[n];
      getpatchface2dir(meshi,i1,i2,j1,j2,k1,k2,patchdir,
        meshi->blockonpatch+n,meshi->meshonpatch+n);
      if(meshi->patchdir[n]==2||meshi->patchdir[n]==-2){
        meshi->patchdir[n]=-meshi->patchdir[n];
      }
    }
    meshi->patch_surfindex[n]=0;
    dxx = 0.0;
    dyy = 0.0;
    dzz = 0.0;
    ig_factor=0.03;
    
    switch (meshi->patchdir[n]){
    case -1:
      meshi->patch_surfindex[n]=0;
      dxx = -meshi->xplt[1]*ig_factor;
      break;
    case 1:
      meshi->patch_surfindex[n]=1;
      dxx = meshi->xplt[1]*ig_factor;
      break;
    case -2:
      meshi->patch_surfindex[n]=2;
      dyy = meshi->yplt[1]*ig_factor;
      break;
    case 2:
      meshi->patch_surfindex[n]=3;
      dyy = -meshi->yplt[1]*ig_factor;
      break;
    case -3:
      meshi->patch_surfindex[n]=4;
      dzz = -meshi->zplt[1]*ig_factor;
      break;
    case 3:
      meshi->patch_surfindex[n]=5;
      dzz = meshi->zplt[1]*ig_factor;
      break;
    default:
      ASSERT(FFALSE);
    }

    
    meshi->patchtype[n]=INTERIORwall;
    if(i1==i2){
      meshi->patchcol[n] = j2 + 1 - j1;
      meshi->patchrow[n] = k2 + 1 - k1;

      ext_wall=0;
      if(j1==0&&j2==jbartemp&&k1==0&&k2==kbartemp){
        if(i1==0||i2==ibartemp){
          ext_wall=1;
          if(i1==0)meshi->patchtype[n]=LEFTwall;
          if(i2==ibartemp)meshi->patchtype[n]=RIGHTwall;
        }
      }
      if(ext_wall==0){
        // an internal wall so set blank to 1 then zero out where there are vents
        for(k=k1;k<=k2;k++){
          if(k==k1){
            dz_factor=-meshi->zplt[1]*ig_factor;
          }
          else if(k==k2){
            dz_factor=meshi->zplt[1]*ig_factor;
          }
          else{
            dz_factor=0.0;
          }
          for(j=j1;j<=j2;j++){
            if(j==j1){
              dy_factor=-meshi->yplt[1]*ig_factor;
            }
            else if(j==j2){
              dy_factor=meshi->yplt[1]*ig_factor;
            }
            else{
              dy_factor=0.0;
            }
            *xyzpatchcopy++ = xplttemp[i1];
            *xyzpatchcopy++ = yplttemp[j];
            *xyzpatchcopy++ = zplttemp[k];
            *xyzpatch_ignitecopy++ = xplttemp[i1]+dxx;
            *xyzpatch_ignitecopy++ = yplttemp[j]+dy_factor;
            *xyzpatch_ignitecopy++ = zplttemp[k]+dz_factor;
            *patchblankcopy++ = nodeinvent(meshi,i1,j,k,1);
          }
        }
      }
      else{
        int ii;

        // an external wall so set blank to 0 then set to one where there are dummy vents
        ii=0;
        for(k=k1;k<=k2;k++){
          if(k==k1){
            dz_factor=-meshi->zplt[1]*ig_factor;
          }
          else if(k==k2){
            dz_factor=meshi->zplt[1]*ig_factor;
          }
          else{
            dz_factor=0.0;
          }
          for(j=j1;j<=j2;j++){
            if(j==j1){
              dy_factor=-meshi->yplt[1]*ig_factor;
            }
            else if(j==j2){
              dy_factor=meshi->yplt[1]*ig_factor;
            }
            else{
              dy_factor=0.0;
            }
            *xyzpatchcopy++ = xplttemp[i1];
            *xyzpatchcopy++ = yplttemp[j];
            *xyzpatchcopy++ = zplttemp[k];
            *xyzpatch_ignitecopy++ = xplttemp[i1]+dxx;
            *xyzpatch_ignitecopy++ = yplttemp[j]+dy_factor;
            *xyzpatch_ignitecopy++ = zplttemp[k]+dz_factor;
            patchblankcopy[ii++]=0;
          }
        }
        nodein_extvent(
          n,
          patchblankcopy,meshi,i1,i2,j1,j2,k1,k2);
        patchblankcopy += (k2+1-k1)*(j2+1-j1);
      }
    }
    else if(j1==j2){
      meshi->patchcol[n] = i2 + 1 - i1;
      meshi->patchrow[n] = k2 + 1 - k1;

      ext_wall=0;
      if(i1==0&&i2==ibartemp&&k1==0&&k2==kbartemp){
        if(j1==0||j2==jbartemp){
          ext_wall=1;
          if(j1==0)meshi->patchtype[n]=FRONTwall;
          if(j2==jbartemp)meshi->patchtype[n]=BACKwall;
        }
      }
      if(ext_wall==0){
        for(k=k1;k<=k2;k++){
          if(k==k1){
            dz_factor=-meshi->zplt[1]*ig_factor;
          }
          else if(k==k2){
            dz_factor=meshi->zplt[1]*ig_factor;
          }
          else{
            dz_factor=0.0;
          }
          for(i=i1;i<=i2;i++){
            if(i==i1){
              dx_factor=-meshi->xplt[1]*ig_factor;
            }
            else if(i==i2){
              dx_factor=meshi->xplt[1]*ig_factor;
            }
            else{
              dx_factor=0.0;
            }
            *xyzpatchcopy++ = xplttemp[i];
            *xyzpatchcopy++ = yplttemp[j1];
            *xyzpatchcopy++ = zplttemp[k];
            *xyzpatch_ignitecopy++ = xplttemp[i]+dx_factor;
            *xyzpatch_ignitecopy++ = yplttemp[j1]+dyy;
            *xyzpatch_ignitecopy++ = zplttemp[k]+dz_factor;
            *patchblankcopy++ = nodeinvent(meshi,i,j1,k,2);
          }
        }
      }
      else{
        int ii;

        // an external wall so set blank to 0 then zero out where there are vents
        ii=0;
        for(k=k1;k<=k2;k++){
          if(k==k1){
            dz_factor=-meshi->zplt[1]*ig_factor;
          }
          else if(k==k2){
            dz_factor=meshi->zplt[1]*ig_factor;
          }
          else{
            dz_factor=0.0;
          }
          for(i=i1;i<=i2;i++){
            if(i==i1){
              dx_factor=-meshi->xplt[1]*ig_factor;
            }
            else if(i==i2){
              dx_factor=meshi->xplt[1]*ig_factor;
            }
            else{
              dx_factor=0.0;
            }
            *xyzpatchcopy++ = xplttemp[i];
            *xyzpatchcopy++ = yplttemp[j1];
            *xyzpatchcopy++ = zplttemp[k];
            *xyzpatch_ignitecopy++ = xplttemp[i]+dx_factor;
            *xyzpatch_ignitecopy++ = yplttemp[j1]+dyy;
            *xyzpatch_ignitecopy++ = zplttemp[k]+dz_factor;
            patchblankcopy[ii++]=0;
          }
        }
        nodein_extvent(
          n,
          patchblankcopy,meshi,i1,i2,j1,j2,k1,k2);
        patchblankcopy += (k2+1-k1)*(i2+1-i1);
      }
    }
    else if(k1==k2){
      meshi->patchcol[n] = i2 + 1 - i1;
      meshi->patchrow[n] = j2 + 1 - j1;

      ext_wall=0;
      if(i1==0&&i2==ibartemp&&j1==0&&j2==jbartemp){
        if(k1==0||k2==kbartemp){
          ext_wall=1;
          if(k1==0)meshi->patchtype[n]=DOWNwall;
          if(k2==kbartemp)meshi->patchtype[n]=UPwall;
        }
      }
      if(ext_wall==0){
        for(j=j1;j<=j2;j++){
          if(j==j1){
            dy_factor=-meshi->yplt[1]*ig_factor;
          }
          else if(j==j2){
            dy_factor=meshi->yplt[1]*ig_factor;
          }
          else{
            dy_factor=0.0;
          }
          for(i=i1;i<=i2;i++){
            if(i==i1){
              dx_factor=-meshi->xplt[1]*ig_factor;
            }
            else if(i==i2){
              dx_factor=meshi->xplt[1]*ig_factor;
            }
            else{
              dx_factor=0.0;
            }
            *xyzpatchcopy++ = xplttemp[i];
            *xyzpatchcopy++ = yplttemp[j];
            *xyzpatchcopy++ = zplttemp[k1];
            *xyzpatch_ignitecopy++ = xplttemp[i]+dx_factor;
            *xyzpatch_ignitecopy++ = yplttemp[j]+dy_factor;
            *xyzpatch_ignitecopy++ = zplttemp[k1]+dzz;
            *patchblankcopy++ = nodeinvent(meshi,i,j,k1,3);
          }
        }
      }
      else{
        int ii;

      // an external wall so set blank to 0 then zero out where there are vents

        ii=0;
        for(j=j1;j<=j2;j++){
          if(j==j1){
            dy_factor=-meshi->yplt[1]*ig_factor;
          }
          else if(j==j2){
            dy_factor=meshi->yplt[1]*ig_factor;
          }
          else{
            dy_factor=0.0;
          }
          for(i=i1;i<=i2;i++){
            if(i==i1){
              dx_factor=-meshi->xplt[1]*ig_factor;
            }
            else if(i==i2){
              dx_factor=meshi->xplt[1]*ig_factor;
            }
            else{
              dx_factor=0.0;
            }
            *xyzpatchcopy++ = xplttemp[i];
            *xyzpatchcopy++ = yplttemp[j];
            *xyzpatchcopy++ = zplttemp[k1];
            *xyzpatch_ignitecopy++ = xplttemp[i]+dx_factor;
            *xyzpatch_ignitecopy++ = yplttemp[j]+dy_factor;
            *xyzpatch_ignitecopy++ = zplttemp[k1]+dzz;
            patchblankcopy[ii++]=0;
          }
        }
        nodein_extvent(
          n,
          patchblankcopy,meshi,i1,i2,j1,j2,k1,k2);
        patchblankcopy += (j2+1-j1)*(i2+1-i1);
      }
    }
    meshi->blockstart[n+1]=meshi->blockstart[n]+meshi->patchrow[n]*meshi->patchcol[n];
    meshi->visPatches[n]=visPatchType[meshi->patchtype[n]];
  }

  meshi->patchfacevis2=0;
  for(n=0;n<meshi->nbptrs;n++){
    bc=meshi->blockageinfoptrs[n];
    bc->patchvis[6]=0;
    for(j=0;j<6;j++){
      bc->patchvis[6]+=bc->patchvis[j];
    }
    if(bc->patchvis[6]!=0){
      bc->patchvis[6]=1;
      meshi->patchfacevis2=1;
    }
  }

  meshi->pqq = NULL;
  switch (loadpatchbysteps){
  case 0:
    while(meshi->pqq==NULL&&mxpatch_frames>100){
      mxpatch_frames-=50;
      meshi->mxpatch_frames=mxpatch_frames;
      NewMemory((void **)&meshi->pqq,sizeof(float)*mxpatch_frames*meshi->npatchsize);
    }
    if(meshi->pqq==NULL){
      NewMemory((void **)&meshi->pqq,sizeof(float)*mxpatch_frames*meshi->npatchsize);
    }
    break;
  case 1:
    npqq = meshi->npatchsize*meshi->mxpatch_frames;
    if(
      NewMemory((void **)&meshi->pqq,sizeof(float)*meshi->npatchsize)==0||
      NewMemory((void **)&meshi->ipqq,sizeof(unsigned char)*npqq)==0){
      *errorcode=1;
      FORTclosefortranfile(&file_unit);
      readpatch(ifile,UNLOAD,&error);
      return;
    }
    break;
  case 2:
    getpatchsizeinfo(patchinfo+ifile, &mxpatch_frames, &ncompressed_buffer);
    NewMemory((void **)&meshi->ipqq_zlib,sizeof(unsigned char)*ncompressed_buffer);
    NewMemory((void **)&meshi->ipqqi_zlib,sizeof(unsigned char)*meshi->npatchsize);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }

  NewMemory((void **)&meshi->patchtimes,sizeof(float)*mxpatch_frames);
  NewMemory((void **)&meshi->zipoffset,sizeof(unsigned int)*mxpatch_frames);
  NewMemory((void **)&meshi->zipsize,sizeof(unsigned int)*mxpatch_frames);
  if(meshi->patchtimes==NULL){
    *errorcode=1;
    if(patchinfo[ifile].compression_type!=2){
      FORTclosefortranfile(&file_unit);
    }
    readpatch(ifile,UNLOAD,&error);
    return;
  }
  if(loadpatchbysteps==2){
    getpatchdata_zlib(patchinfo+ifile,meshi->ipqq_zlib,ncompressed_buffer,
      meshi->patchtimes,meshi->zipoffset,meshi->zipsize,mxpatch_frames);
    meshi->npatch_frames=mxpatch_frames;
  }
  else{
    if(meshi->pqq==NULL){
      *errorcode=1;
      FORTclosefortranfile(&file_unit);
      readpatch(ifile,UNLOAD,&error);
      return;
    }
    meshi->npatch_frames=0;
  }

  file_size=get_filesize(file);
  local_starttime = glutGet(GLUT_ELAPSED_TIME);
  for (i=0;i<mxpatch_frames;){
    if(loadpatchbysteps==1){
      meshi->pqqi = meshi->pqq;
      meshi->ipqqi = meshi->ipqq + i*meshi->npatchsize;
    }
    else if(loadpatchbysteps==0){
      meshi->pqqi = meshi->pqq + i*meshi->npatchsize;
    }
    meshi->patchtimesi = meshi->patchtimes + i;

    error=0;
    if(loadpatchbysteps==0||loadpatchbysteps==1){
      for (n=0;n<boundframestep;n++){
        if(error==0){
          int npqqi;

          FORTgetpatchdata(&file_unit,&meshi->npatches,
          meshi->pi1,meshi->pi2,
          meshi->pj1,meshi->pj2,
          meshi->pk1,meshi->pk2,
          meshi->patchtimesi,meshi->pqqi,&npqqi,
          &error);
        }
      }
    }
    if(do_threshold==1){
      if(local_first==1){
        nn=0;
        for(n=0;n<meshi->npatchsize;n++){
          meshi->thresholdtime[n]=-1.0;
        }
        local_first=0;
      }

      /* create char contour plot for each patch corresponding to a blockage */

      {
        float chartemp_ign;
        surface *surfi;
        int surf_index;

        nn=0;
#ifdef USE_ZLIB
        if(loadpatchbysteps==2)uncompress_patchdataframe(meshi,i);
#endif
        for(n=0;n<meshi->npatches;n++){
          mesh *meshblock;
          float dval;

          iblock=meshi->blockonpatch[n];
          meshblock = meshi->meshonpatch[n];
          nsize=meshi->patchrow[n]*meshi->patchcol[n];
          ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
          if(iblock!=-1&&meshblock!=NULL){
            switch (loadpatchbysteps){
            case 0:
            case 1:
              for(j=0;j<nsize;j++){
                if(meshi->thresholdtime[nn+j]<0.0&&meshi->pqqi[nn+j]>=temp_threshold){
                  meshi->thresholdtime[nn+j]=meshi->patchtimes[i];
                }
              }
              break;
            case 2:
              dval = (patchmax-patchmin)/255.0;
              for(j=0;j<nsize;j++){
                float val;
                int ival;

                ival = meshi->ipqqi_zlib[nn+j];
                val = patchmin + dval*ival;
                if(meshi->thresholdtime[nn+j]<0.0&&val>=temp_threshold){
                  meshi->thresholdtime[nn+j]=meshi->patchtimes[i];
                }
              }
              break;
            default:
              ASSERT(FFALSE);
              break;
            }
          }
          nn+=nsize;
        }
      }
    }
    if(loadpatchbysteps==1){
      getBoundaryColors2(
        meshi->pqqi, meshi->npatchsize, meshi->ipqqi, 
                 setpatchmin,&patchmin, setpatchmax,&patchmax, 
                 &patchmin_global, &patchmax_global,
                 nrgb_full,
                 &patchinfo[ifile].extreme_min,&patchinfo[ifile].extreme_max);
    }
    CheckMemory;
    if(error!=0)break;
    if(settmax_b!=0&&*meshi->patchtimesi>tmax_b)break;

    if(loadpatchbysteps==0||loadpatchbysteps==1){
      if(!(settmin_b!=0&&*meshi->patchtimesi<tmin_b)){
        printf("boundary time=%.2f\n",*meshi->patchtimesi);

        meshi->npatch_frames++;
        if(meshi->npatch_frames + 1 > mxpatch_frames){
          printf("reallocating memory\n");
          mxpatch_frames = meshi->npatch_frames + 50; /* this + 50 must match - 50 below */
          meshi->mxpatch_frames=mxpatch_frames;
          if(
            ResizeMemory((void **)&meshi->pqq,           mxpatch_frames*meshi->npatchsize*sizeof(float))==0||
            ResizeMemory((void **)&meshi->patchtimes,    mxpatch_frames*sizeof(float))==0
           ){
            *errorcode=1;
            readpatch(ifile,UNLOAD,&error);
            FORTclosefortranfile(&file_unit);
            return;
          }
        }
        i++;
      }
    }
    if(loadpatchbysteps==2){
      i++;
    }
  }
  local_stoptime = glutGet(GLUT_ELAPSED_TIME);
  delta_time = (local_stoptime-local_starttime)/1000.0;
  
  /* convert patch values into integers pointing to an rgb color table */

  if(loadpatchbysteps==0){
    npqq = meshi->npatch_frames*meshi->npatchsize;
    if(npqq==0||NewMemory((void **)&meshi->ipqq,sizeof(unsigned char)*npqq)==0){
      *errorcode=1;
      FORTclosefortranfile(&file_unit);
      readpatch(ifile,UNLOAD,&error);
      return;
    }
  }

  printf("computing boundary color levels \n");
  if(NewMemory((void **)&colorlabelpatch,MAXRGB*sizeof(char *))==0){
    *errorcode=1;
    if(loadpatchbysteps!=2){
      FORTclosefortranfile(&file_unit);
    }
    readpatch(ifile,UNLOAD,&error);
    return;
  }
  for(n=0;n<MAXRGB;n++){
    colorlabelpatch[n]=NULL;
  }
  for(n=0;n<nrgb;n++){
    if(NewMemory((void **)&colorlabelpatch[n],11)==0){
      *errorcode=1;
      if(loadpatchbysteps!=2){
        FORTclosefortranfile(&file_unit);
      }
      readpatch(ifile,UNLOAD,&error);
      return;
    }
  }
  patchi = patchinfo + ifile;
  patchscale = patchi->scale;
  patchbase = patchinfo + getpatchindex(patchi);
  npatchloaded=0;
  patchinfo[ifile].loaded=1;
  ipatchtype=getpatchtype(patchinfo+ifile);
  for(i=0;i<npatch_files;i++){
    if(patchinfo[i].loaded==1&&patchinfo[i].type==ipatchtype)npatchloaded++;
  }
  switch(loadpatchbysteps){
  case 0:
    getBoundaryColors3(patchinfo + ifile,meshi->pqq, npqq, meshi->ipqq, 
      setpatchmin,&patchmin, setpatchmax,&patchmax, 
      &patchmin_global, &patchmax_global,
      nrgb, colorlabelpatch,patchscale,boundarylevels256,
      &patchinfo[ifile].extreme_min,&patchinfo[ifile].extreme_max);
    break;
  case 1:
    getBoundaryLabels(
      patchmin, patchmax, 
      colorlabelpatch,patchscale,boundarylevels256,nrgb);
    break;
  case 2:
    getBoundaryLabels(
      patchmin, patchmax, 
      colorlabelpatch,patchscale,boundarylevels256,nrgb);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  strcpy(patchbase->scale,patchi->scale);
  if(do_threshold==1){
    meshi->surface_tempmax=patchmax_global;
    meshi->surface_tempmin=patchmin_global;
  }

  local2globalpatchbounds(patchinfo[ifile].label.shortlabel);
  updatepatchlistindex(patchfilenum);

  FREEMEMORY(meshi->pqq);
  patchinfo[ifile].loaded=1;
  patchinfo[ifile].display=1;
  ipatchtype=getpatchtype(patchinfo+ifile);
  showexterior=1-showexterior;
  allexterior = 1-allexterior;
  ShowPatchMenu(EXTERIORwallmenu);
  plotstate=getplotstate(DYNAMIC_PLOTS);
  if(patchinfo[ifile].compression_type==1)disable_boundary_glui();
  updatetimes();
  update_unit_defs();
  updatechopcolors();
#ifdef _DEBUG
  printf("After boundary file load: ");
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

/* ------------------ nodeinblockage ------------------------ */

int nodeinblockage(const mesh *meshnode, int i,int j,int k, int *imesh, int *iblockage){
  int ii,jj;
  mesh *meshii;
  blockagedata *bc;
  float xn, yn, zn;
  float xm_min, xm_max;
  float ym_min, ym_max;
  float zm_min, zm_max;
  float xb_min, xb_max;
  float yb_min, yb_max;
  float zb_min, zb_max;
  float *xplt, *yplt, *zplt;


  xn = meshnode->xplt[i];
  yn = meshnode->yplt[j];
  zn = meshnode->zplt[k];

  *imesh=-1;

  for(ii=0;ii<nmeshes;ii++){
    meshii = meshinfo + ii;
    if(meshnode==meshii)continue;

    xplt = meshii->xplt;
    yplt = meshii->yplt;
    zplt = meshii->zplt;

    xm_min=xplt[0];
    xm_max=meshii->xbar;
    ym_min=yplt[0];
    ym_max=meshii->ybar;
    zm_min=zplt[0];
    zm_max=meshii->zbar;
    if(xn<xm_min||xn>xm_max)continue;
    if(yn<ym_min||yn>ym_max)continue;
    if(zn<zm_min||zn>zm_max)continue;


    for(jj=0;jj<meshii->nbptrs;jj++){
      bc = meshii->blockageinfoptrs[jj];
      if(bc->hole==1)continue;
      /*
      xb_min=bc->xmin;
      xb_max=bc->xmax;
      yb_min=bc->ymin;
      yb_max=bc->ymax;
      zb_min=bc->zmin;
      zb_max=bc->zmax;
      */
      xb_min=xplt[bc->ijk[0]];
      xb_max=xplt[bc->ijk[1]];
      yb_min=yplt[bc->ijk[2]];
      yb_max=yplt[bc->ijk[3]];
      zb_min=zplt[bc->ijk[4]];
      zb_max=zplt[bc->ijk[5]];
      if(xb_min<=xn&&xn<=xb_max&&
         yb_min<=yn&&yn<=yb_max&&
         zb_min<=zn&&zn<=zb_max){
        *imesh=ii;
        *iblockage=jj;
        return 1;
      }
    }
  }
  return 0;
}

/* ------------------ nodeinvent ------------------------ */

int nodeinvent(const mesh *meshi, int i,int j,int k, int dir){
  int ii;
  ventdata *vi;
  vi=meshi->ventinfo;
  for(ii=0;ii<meshi->nvents;ii++){
    vi = meshi->ventinfo+ii;
    if(vi->hideboundary==1){
      switch (dir){
      case 1:
        if(vi->imin==i&&i==vi->imax&&
           vi->jmin <j&&j <vi->jmax&&
           vi->kmin <k&&k <vi->kmax){
          return 0;
        }
        break;
      case 2:
        if(vi->jmin==j&&j==vi->jmax&&
           vi->imin <i&&i <vi->imax&&
           vi->kmin <k&&k <vi->kmax){
          return 0;
        }
        break;
      case 3:
        if(vi->kmin==k&&k==vi->kmax&&
           vi->imin <i&&i <vi->imax&&
           vi->jmin <j&&j <vi->jmax){
          return 0;
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
      }
    }
  }
  return 1;
}

/* ------------------ nodeinvent ------------------------ */

void nodein_extvent(
                    int ipatch, 
                    int *patchblank, const mesh *meshi, 
                    int i1, int i2, 
                    int j1, int j2, 
                    int k1, int k2){
  int ii,iii;
  ventdata *vi;
  int imin, jmin, kmin, imax, jmax, kmax;
  int i, j, k;
  int dir;
  int imesh, iblockage;
  mesh *meshblock;

  dir=0;
  if(i1==i2)dir=1;
  if(j1==j2)dir=2;
  if(k1==k2)dir=3;

  for(ii=0;ii<meshi->nvents;ii++){
    vi = meshi->ventinfo+ii;
    if(vi->hideboundary==1)continue;
    if(vi->dir2!=dir)continue;
    switch (dir){
    case 1:
      if(vi->imin!=i1)continue;
      if(vi->jmax<j1)continue;
      if(vi->jmin>j2)continue;
      if(vi->kmax<k1)continue;
      if(vi->kmin>k2)continue;
      jmin=vi->jmin;
      jmax=vi->jmax;
      if(jmin<j1)jmin=j1;
      if(jmax>j2)jmax=j2;
      kmin=vi->kmin;
      kmax=vi->kmax;
      if(kmin<k1)kmin=k1;
      if(kmax>k2)kmax=k2;
      for(k=kmin;k<=kmax;k++){
        for(j=jmin;j<=jmax;j++){
          iii=(k-k1)*(j2+1-j1) + (j-j1);
          patchblank[iii]=1;
        }
      }
      break;
    case 2:
      if(vi->jmin!=j1)continue;
      if(vi->imax<i1)continue;
      if(vi->imin>i2)continue;
      if(vi->kmax<k1)continue;
      if(vi->kmin>k2)continue;
      imin=vi->imin;
      imax=vi->imax;
      if(imin<i1)imin=i1;
      if(imax>i2)imax=i2;
      kmin=vi->kmin;
      kmax=vi->kmax;
      if(kmin<k1)kmin=k1;
      if(kmax>k2)kmax=k2;
      for(k=kmin;k<=kmax;k++){
        for(i=imin;i<=imax;i++){
          iii=(k-k1)*(i2+1-i1) + (i-i1);
          patchblank[iii]=1;
        }
      }
      break;
    case 3:
      if(vi->kmin!=k1)continue;
      if(vi->imax<i1)continue;
      if(vi->imin>i2)continue;
      if(vi->jmax<j1)continue;
      if(vi->jmin>j2)continue;
      imin=vi->imin;
      imax=vi->imax;
      if(imin<i1)imin=i1;
      if(imax>i2)imax=i2;
      jmin=vi->jmin;
      jmax=vi->jmax;
      if(jmin<j1)jmin=j1;
      if(jmax>j2)jmax=j2;
      for(j=jmin;j<=jmax;j++){
        for(i=imin;i<=imax;i++){
          iii=(j-j1)*(i2+1-i1) + (i-i1);
          patchblank[iii]=1;
        }
      }
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
  }
  switch (dir){
  case 1:
    for(k=k1;k<=k2;k++){
      for(j=j1;j<=j2;j++){
        iii=(k-k1)*(j2+1-j1) + (j-j1);
        if(patchblank[iii]==1)continue;
        patchblank[iii] = nodeinblockage(meshi,i1,j,k,&imesh,&iblockage);
      }
    }
    break;
  case 2:
    for(k=k1;k<=k2;k++){
      for(i=i1;i<=i2;i++){
        iii=(k-k1)*(i2+1-i1) + (i-i1);
        if(patchblank[iii]==1)continue;
        patchblank[iii]=nodeinblockage(meshi,i,j1,k,&imesh,&iblockage);
        if(imesh!=-1){
          meshblock = meshinfo+imesh;
          ASSERT(iblockage>=0&&iblockage<meshblock->nbptrs);
          meshi->blockonpatch[ipatch]=iblockage;
          meshi->meshonpatch[ipatch]=meshblock;
        }
      }
    }
    break;
  case 3:
    for(j=j1;j<=j2;j++){
      for(i=i1;i<=i2;i++){
        iii=(j-j1)*(i2+1-i1) + (i-i1);
        if(patchblank[iii]==1)continue;
        patchblank[iii]=nodeinblockage(meshi,i,j,k1,&imesh,&iblockage);
      }
    }
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
}

/* ------------------ ispatchtype ------------------------ */

int ispatchtype(int type){
  int n,i;
  mesh *meshi;

  for(i=0;i<nmeshes;i++){
    meshi=meshinfo+i;
    for(n=0;n<meshi->npatches;n++){
      if(meshi->patchtype[n]==type)return 1;
    }
  }
  return 0;
}

/* ------------------ local2globalpatchbounds ------------------------ */

void local2globalpatchbounds(const char *key){
  int i;
  
  for(i=0;i<npatch_files;i++){
    if(strcmp(patchinfo[i].label.shortlabel,key)==0){
      patchinfo[i].valmin=patchmin;
      patchinfo[i].valmax=patchmax;
      patchinfo[i].setvalmin=setpatchmin;
      patchinfo[i].setvalmax=setpatchmax;

      patchinfo[i].chopmin=patchchopmin;
      patchinfo[i].chopmax=patchchopmax;
      patchinfo[i].setchopmin=setpatchchopmin;
      patchinfo[i].setchopmax=setpatchchopmax;

    }
  }
}

/* ------------------ global2localpatchbounds ------------------------ */

void global2localpatchbounds(const char *key){
  int i;
  
  for(i=0;i<npatch_files;i++){
    if(strcmp(patchinfo[i].label.shortlabel,key)==0){
      patchmin=patchinfo[i].valmin;
      patchmax=patchinfo[i].valmax;
      setpatchmin=patchinfo[i].setvalmin;
      setpatchmax=patchinfo[i].setvalmax;

      patchchopmin=patchinfo[i].chopmin;
      patchchopmax=patchinfo[i].chopmax;
      setpatchchopmin=patchinfo[i].setchopmin;
      setpatchchopmax=patchinfo[i].setchopmax;
      update_hidepatchsurface();
      
      local2globalpatchbounds(key);
      return;
    }
  }
}

/* ------------------ drawpatch_texture ------------------------ */

void drawpatch_texture(const mesh *meshi){
  float r11, r12, r21, r22;
  int n,nn;
  int nrow, ncol, irow, icol;
  float *xyzp1, *xyzp2;
  unsigned char *ipq1, *ipq2;
  unsigned char *ipqqcopy;
  float *xyzpatchcopy;
  int *patchblankcopy, *ip1, *ip2;
  float *patchtimes;
  int *visPatches;
  float *xyzpatch;
  int *patchdir, *patchrow, *patchcol;
  int *blockstart;
  int *patchblank;
  unsigned char *ipqqi;
  int iblock;
  blockagedata *bc;
  patch *patchi;
  mesh *meshblock;
  float clear_color[4]={1.0,1.0,1.0,1.0};
  float dboundx,dboundy,dboundz;
  float *xplt, *yplt, *zplt;

  if(vis_threshold==1&&vis_onlythreshold==1&&do_threshold==1)return;

  if(hidepatchsurface==0){
    xplt=meshi->xplt;
    yplt=meshi->yplt;
    zplt=meshi->zplt;

    dboundx = (xplt[1]-xplt[0])/10.0;
    dboundy = (yplt[1]-yplt[0])/10.0;
    dboundz = (zplt[1]-zplt[0])/10.0;
  }

  patchtimes=meshi->patchtimes;
  visPatches=meshi->visPatches;
  xyzpatch=meshi->xyzpatch;
  patchdir=meshi->patchdir;
  patchrow=meshi->patchrow;
  patchcol=meshi->patchcol;
  blockstart=meshi->blockstart;
  patchblank=meshi->patchblank;
  patchi=patchinfo+meshi->patchfilenum;
  switch(patchi->compression_type){
  case 0:
    ASSERT(meshi->ipqqi!=NULL);
    ipqqi=meshi->ipqqi;
    break;
  case 1:
    ASSERT(meshi->ipqqi_zlib!=NULL);
    ipqqi=meshi->ipqqi_zlib;
    break;
  default:
    ASSERT(0);
  }
  patchi = patchinfo + meshi->patchfilenum;

  if(patchtimes[0]>times[itime]||patchi->display==0)return;
  if(cullfaces==1)glDisable(GL_CULL_FACE);

  /* if a contour boundary does not match a blockage face then draw "both sides" of boundary */

  if((use_transparency_data==1&&showcolorbarlines==1)||setpatchchopmin==1||setpatchchopmax==1)transparenton();
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
  glEnable(GL_TEXTURE_1D);
  glBindTexture(GL_TEXTURE_1D,texture_patch_colorbar_id);

  nn =0;
  glBegin(GL_TRIANGLES);
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(visPatches[n]==1&&patchdir[n]==0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      ipqqcopy = ipqqi + blockstart[n];

      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        ipq1 = ipqqcopy + irow*ncol;
        ipq2 = ipq1 + ncol;
        ip2 = ip1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            r11 = (float)(*ipq1)/255.0;
            r12 = (float)(*(ipq1+1))/255.0;
            r21 = (float)(*ipq2)/255.0;
            r22 = (float)(*(ipq2+1))/255.0;
            if(abs(r11-r22)<abs(r12-r21)){
              glTexCoord1f(r11);glVertex3fv(xyzp1);
              glTexCoord1f(r12);glVertex3fv(xyzp1+3);
              glTexCoord1f(r22);glVertex3fv(xyzp2+3);
              glTexCoord1f(r11);glVertex3fv(xyzp1);
              glTexCoord1f(r22);glVertex3fv(xyzp2+3);
              glTexCoord1f(r21);glVertex3fv(xyzp2);

            }
            else{
              glTexCoord1f(r11);glVertex3fv(xyzp1);
              glTexCoord1f(r12);glVertex3fv(xyzp1+3);
              glTexCoord1f(r21);glVertex3fv(xyzp2);
              glTexCoord1f(r12);glVertex3fv(xyzp1+3);
              glTexCoord1f(r22);glVertex3fv(xyzp2+3);
              glTexCoord1f(r21);glVertex3fv(xyzp2);
            }
          }
          ipq1++; ipq2++; ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
    }
    nn += patchrow[n]*patchcol[n];
  }
  glEnd();
  if(cullfaces==1)glEnable(GL_CULL_FACE);

  /* if a contour boundary DOES match a blockage face then draw "one sides" of boundary */

  nn=0;
  if(hidepatchsurface==1){
    glBegin(GL_TRIANGLES);
  }
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock=meshi->meshonpatch[n];
    if(iblock!=-1){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(meshi->visPatches[n]==1&&meshi->patchdir[n]>0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      ipqqcopy = ipqqi + blockstart[n];
      if(hidepatchsurface==0){
        glPushMatrix();
        switch (meshi->patchdir[n]){
          case 1:
            glTranslatef(dboundx,0.0,0.0);
            break;
          case 2:
            glTranslatef(0.0,-dboundy,0.0);
            break;
          case 3:
            glTranslatef(0.00,0.0,dboundz);
            break;
        }
        glBegin(GL_TRIANGLES);
      }
      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ipq1 = ipqqcopy + irow*ncol;
        ip1 = patchblankcopy + irow*ncol;

        xyzp2 = xyzp1 + 3*ncol;
        ipq2 = ipq1 + ncol;
        ip2 = ip1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            r11 = (float)(*ipq1)/255.0;
            r12 = (float)(*(ipq1+1))/255.0;
            r21 = (float)(*ipq2)/255.0;
            r22 = (float)(*(ipq2+1))/255.0;
            if(abs(*ipq1-*(ipq2+1))<abs(*(ipq1+1)-*ipq2)){
              glTexCoord1f(r11);glVertex3fv(xyzp1);
              glTexCoord1f(r12);glVertex3fv(xyzp1+3);
              glTexCoord1f(r22);glVertex3fv(xyzp2+3);
              glTexCoord1f(r11);glVertex3fv(xyzp1);
              glTexCoord1f(r22);glVertex3fv(xyzp2+3);
              glTexCoord1f(r21);glVertex3fv(xyzp2);
            }
            else{
              glTexCoord1f(r11);glVertex3fv(xyzp1);
              glTexCoord1f(r12);glVertex3fv(xyzp1+3);
              glTexCoord1f(r21);glVertex3fv(xyzp2);
              glTexCoord1f(r12);glVertex3fv(xyzp1+3);
              glTexCoord1f(r22);glVertex3fv(xyzp2+3);
              glTexCoord1f(r21);glVertex3fv(xyzp2);
            }
          }
          ipq1++; ipq2++; ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
      if(hidepatchsurface==0){
        glEnd();
        glPopMatrix();
      }
    }
    nn += patchrow[n]*patchcol[n];
  }

  /* if a contour boundary DOES match a blockage face then draw "one sides" of boundary */
  nn=0;
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(visPatches[n]==1&&patchdir[n]<0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      ipqqcopy = ipqqi + blockstart[n];
      if(hidepatchsurface==0){
        glPushMatrix();
        switch (meshi->patchdir[n]){
          case -1:
            glTranslatef(-dboundx,0.0,0.0);
            break;
          case -2:
            glTranslatef(0.0,dboundy,0.0);
            break;
          case -3:
            glTranslatef(0.0,0.0,-dboundz);
            break;
        }
        glBegin(GL_TRIANGLES);
      }
      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        ipq1 = ipqqcopy + irow*ncol;
        ipq2 = ipq1 + ncol;
        ip2 = ip1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            r11 = (float)(*ipq1)/255.0;
            r12 = (float)(*(ipq1+1))/255.0;
            r21 = (float)(*ipq2)/255.0;
            r22 = (float)(*(ipq2+1))/255.0;
            if(abs(*ipq1-*(ipq2+1))<abs(*(ipq1+1)-*ipq2)){
              glTexCoord1f(r11);glVertex3fv(xyzp1);
              glTexCoord1f(r22);glVertex3fv(xyzp2+3);
              glTexCoord1f(r12);glVertex3fv(xyzp1+3);
              glTexCoord1f(r11);glVertex3fv(xyzp1);
              glTexCoord1f(r21);glVertex3fv(xyzp2);
              glTexCoord1f(r22);glVertex3fv(xyzp2+3);
            }
            else{
              glTexCoord1f(r11);glVertex3fv(xyzp1);
              glTexCoord1f(r21);glVertex3fv(xyzp2);
              glTexCoord1f(r12);glVertex3fv(xyzp1+3);
              glTexCoord1f(r12);glVertex3fv(xyzp1+3);
              glTexCoord1f(r21);glVertex3fv(xyzp2);
              glTexCoord1f(r22);glVertex3fv(xyzp2+3);
            }
          }
          ipq1++; ipq2++; ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
      if(hidepatchsurface==0){
        glEnd();
        glPopMatrix();
      }
    }
    nn += patchrow[n]*patchcol[n];
  }
  if(hidepatchsurface==1){
    glEnd();
  }
  glDisable(GL_TEXTURE_1D);
  if((use_transparency_data==1&&showcolorbarlines==1)||setpatchchopmin==1||setpatchchopmax==1)transparentoff();
}

/* ------------------ drawpatch_texture_threshold ------------------------ */

void drawpatch_texture_threshold(const mesh *meshi){
  float r11, r12, r21, r22;
  int n,nn,nn1,nn2;
  int nrow, ncol, irow, icol;
  float *xyzp1, *xyzp2;
  unsigned char *ipq1, *ipq2;
  unsigned char *ipqqcopy;
  float *xyzpatchcopy;
  int *patchblankcopy, *ip1, *ip2;
  float *patchtimes;
  int *visPatches;
  float *xyzpatch;
  int *patchdir, *patchrow, *patchcol;
  int *blockstart;
  int *patchblank;
  unsigned char *ipqqi;
  int iblock;
  blockagedata *bc;
  patch *patchi;
  float *color11, *color12, *color21, *color22;
  mesh *meshblock;
  float burn_color[4]={0.0,0.0,0.0,1.0};
  float clear_color[4]={1.0,1.0,1.0,1.0};

  if(vis_threshold==1&&vis_onlythreshold==1&&do_threshold==1)return;

  patchtimes=meshi->patchtimes;
  visPatches=meshi->visPatches;
  xyzpatch=meshi->xyzpatch;
  patchdir=meshi->patchdir;
  patchrow=meshi->patchrow;
  patchcol=meshi->patchcol;
  blockstart=meshi->blockstart;
  patchblank=meshi->patchblank;
  patchi=patchinfo+meshi->patchfilenum;
  switch(patchi->compression_type){
  case 0:
    ASSERT(meshi->ipqqi!=NULL);
    ipqqi=meshi->ipqqi;
    break;
  case 1:
    ASSERT(meshi->ipqqi_zlib!=NULL);
    ipqqi=meshi->ipqqi_zlib;
    break;
  default:
    ASSERT(0);
  }
  patchi = patchinfo + meshi->patchfilenum;

  if(patchtimes[0]>times[itime]||patchi->display==0)return;
  if(cullfaces==1)glDisable(GL_CULL_FACE);

  /* if a contour boundary does not match a blockage face then draw "both sides" of boundary */

  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
  glEnable(GL_TEXTURE_1D);
  glBindTexture(GL_TEXTURE_1D,texture_colorbar_id);

  nn =0;
  glBegin(GL_TRIANGLES);
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(visPatches[n]==1&&patchdir[n]==0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      ipqqcopy = ipqqi + blockstart[n];

      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        ipq1 = ipqqcopy + irow*ncol;
        ipq2 = ipq1 + ncol;
        ip2 = ip1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            r11 = (float)(*ipq1)/255.0;
            r12 = (float)(*(ipq1+1))/255.0;
            r21 = (float)(*ipq2)/255.0;
            r22 = (float)(*(ipq2+1))/255.0;
            color11=clear_color;
            color12=clear_color;
            color21=clear_color;
            color22=clear_color;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol  ])color11=burn_color;
            if(meshi->thresholdtime[nn1+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol+1])color12=burn_color;
            if(meshi->thresholdtime[nn2+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol  ])color21=burn_color;
            if(meshi->thresholdtime[nn2+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol+1])color22=burn_color;
            if(color11==color12&&color11==color21&&color11==color22){
              glColor4fv(color11);
              glTexCoord1f(r11);glVertex3fv(xyzp1);
              glTexCoord1f(r12);glVertex3fv(xyzp1+3);
              glTexCoord1f(r22);glVertex3fv(xyzp2+3);
              glTexCoord1f(r11);glVertex3fv(xyzp1);
              glTexCoord1f(r22);glVertex3fv(xyzp2+3);
              glTexCoord1f(r21);glVertex3fv(xyzp2);
            }
            else{
             if(abs(r11-r22)<abs(r12-r21)){
               glTexCoord1f(r11);glColor4fv(color11);glVertex3fv(xyzp1);
               glTexCoord1f(r12);glColor4fv(color12);glVertex3fv(xyzp1+3);
               glTexCoord1f(r22);glColor4fv(color22);glVertex3fv(xyzp2+3);
               glTexCoord1f(r11);glColor4fv(color11);glVertex3fv(xyzp1);
               glTexCoord1f(r22);glColor4fv(color22);glVertex3fv(xyzp2+3);
               glTexCoord1f(r21);glColor4fv(color21);glVertex3fv(xyzp2);
             }
             else{
               glTexCoord1f(r11);glColor4fv(color11);glVertex3fv(xyzp1);
               glTexCoord1f(r12);glColor4fv(color12);glVertex3fv(xyzp1+3);
               glTexCoord1f(r21);glColor4fv(color21);glVertex3fv(xyzp2);
               glTexCoord1f(r12);glColor4fv(color12);glVertex3fv(xyzp1+3);
               glTexCoord1f(r22);glColor4fv(color22);glVertex3fv(xyzp2+3);
               glTexCoord1f(r21);glColor4fv(color21);glVertex3fv(xyzp2);
             }
            }
          }
          ipq1++; ipq2++; ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
    }
    nn += patchrow[n]*patchcol[n];
  }
  glEnd();
  if(cullfaces==1)glEnable(GL_CULL_FACE);

  /* if a contour boundary DOES match a blockage face then draw "one sides" of boundary */

  nn=0;
  glBegin(GL_TRIANGLES);
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock=meshi->meshonpatch[n];
    if(iblock!=-1){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(meshi->visPatches[n]==1&&meshi->patchdir[n]>0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      ipqqcopy = ipqqi + blockstart[n];
      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ipq1 = ipqqcopy + irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;

        xyzp2 = xyzp1 + 3*ncol;
        ipq2 = ipq1 + ncol;
        ip2 = ip1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            r11 = (float)(*ipq1)/255.0;
            r12 = (float)(*(ipq1+1))/255.0;
            r21 = (float)(*ipq2)/255.0;
            r22 = (float)(*(ipq2+1))/255.0;
            color11=clear_color;
            color12=clear_color;
            color21=clear_color;
            color22=clear_color;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol  ])color11=burn_color;
            if(meshi->thresholdtime[nn1+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol+1])color12=burn_color;
            if(meshi->thresholdtime[nn2+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol  ])color21=burn_color;
            if(meshi->thresholdtime[nn2+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol+1])color22=burn_color;
            if(color11==color12&&color11==color21&&color11==color22){
              glColor4fv(color11);
              glTexCoord1f(r11);glVertex3fv(xyzp1);
              glTexCoord1f(r12);glVertex3fv(xyzp1+3);
              glTexCoord1f(r22);glVertex3fv(xyzp2+3);
              glTexCoord1f(r11);glVertex3fv(xyzp1);
              glTexCoord1f(r22);glVertex3fv(xyzp2+3);
              glTexCoord1f(r21);glVertex3fv(xyzp2);
            }
            else{
              if(abs(*ipq1-*(ipq2+1))<abs(*(ipq1+1)-*ipq2)){
                glTexCoord1f(r11);glColor4fv(color11);glVertex3fv(xyzp1);
                glTexCoord1f(r12);glColor4fv(color12);glVertex3fv(xyzp1+3);
                glTexCoord1f(r22);glColor4fv(color22);glVertex3fv(xyzp2+3);
                glTexCoord1f(r11);glColor4fv(color11);glVertex3fv(xyzp1);
                glTexCoord1f(r22);glColor4fv(color22);glVertex3fv(xyzp2+3);
                glTexCoord1f(r21);glColor4fv(color21);glVertex3fv(xyzp2);
              }
              else{
                glTexCoord1f(r11);glColor4fv(color11);glVertex3fv(xyzp1);
                glTexCoord1f(r12);glColor4fv(color12);glVertex3fv(xyzp1+3);
                glTexCoord1f(r21);glColor4fv(color21);glVertex3fv(xyzp2);
                glTexCoord1f(r12);glColor4fv(color12);glVertex3fv(xyzp1+3);
                glTexCoord1f(r22);glColor4fv(color22);glVertex3fv(xyzp2+3);
                glTexCoord1f(r21);glColor4fv(color21);glVertex3fv(xyzp2);
              }
            }
          }
          ipq1++; ipq2++; ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
    }
    nn += patchrow[n]*patchcol[n];
  }

  /* if a contour boundary DOES match a blockage face then draw "one sides" of boundary */
  nn=0;
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(visPatches[n]==1&&patchdir[n]<0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      ipqqcopy = ipqqi + blockstart[n];
      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        ipq1 = ipqqcopy + irow*ncol;
        ipq2 = ipq1 + ncol;
        ip2 = ip1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            r11 = (float)(*ipq1)/255.0;
            r12 = (float)(*(ipq1+1))/255.0;
            r21 = (float)(*ipq2)/255.0;
            r22 = (float)(*(ipq2+1))/255.0;
            color11=clear_color;
            color12=clear_color;
            color21=clear_color;
            color22=clear_color;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol  ])color11=burn_color;
            if(meshi->thresholdtime[nn1+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol+1])color12=burn_color;
            if(meshi->thresholdtime[nn2+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol  ])color21=burn_color;
            if(meshi->thresholdtime[nn2+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol+1])color22=burn_color;
            if(color11==color12&&color11==color21&&color11==color22){
              glColor4fv(color11);
              glTexCoord1f(r11);glVertex3fv(xyzp1);
              glTexCoord1f(r22);glVertex3fv(xyzp2+3);
              glTexCoord1f(r12);glVertex3fv(xyzp1+3);
              glTexCoord1f(r11);glVertex3fv(xyzp1);
              glTexCoord1f(r21);glVertex3fv(xyzp2);
              glTexCoord1f(r22);glVertex3fv(xyzp2+3);
            }
            else{
              if(abs(*ipq1-*(ipq2+1))<abs(*(ipq1+1)-*ipq2)){
                glTexCoord1f(r11);glColor4fv(color11);glVertex3fv(xyzp1);
                glTexCoord1f(r22);glColor4fv(color22);glVertex3fv(xyzp2+3);
                glTexCoord1f(r12);glColor4fv(color12);glVertex3fv(xyzp1+3);
                glTexCoord1f(r11);glColor4fv(color11);glVertex3fv(xyzp1);
                glTexCoord1f(r21);glColor4fv(color21);glVertex3fv(xyzp2);
                glTexCoord1f(r22);glColor4fv(color22);glVertex3fv(xyzp2+3);
              }
              else{
                glTexCoord1f(r11);glColor4fv(color11);glVertex3fv(xyzp1);
                glTexCoord1f(r21);glColor4fv(color21);glVertex3fv(xyzp2);
                glTexCoord1f(r12);glColor4fv(color12);glVertex3fv(xyzp1+3);
                glTexCoord1f(r12);glColor4fv(color12);glVertex3fv(xyzp1+3);
                glTexCoord1f(r21);glColor4fv(color21);glVertex3fv(xyzp2);
                glTexCoord1f(r22);glColor4fv(color22);glVertex3fv(xyzp2+3);
              }
            }
          }
          ipq1++; ipq2++; ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
    }
    nn += patchrow[n]*patchcol[n];
  }
  glEnd();
  glDisable(GL_TEXTURE_1D);
}

/* ------------------ drawpatch_threshold_cellcenter ------------------------ */

void drawpatch_threshold_cellcenter(const mesh *meshi){
  float r11, r12, r21, r22;
  int n,nn,nn1;
  int nrow, ncol, irow, icol;
  float *xyzp1, *xyzp2;
  unsigned char *ipq1, *ipq2;
  unsigned char *ipqqcopy;
  float *xyzpatchcopy;
  int *patchblankcopy, *ip1, *ip2;
  float *patchtimes;
  int *visPatches;
  float *xyzpatch;
  int *patchdir, *patchrow, *patchcol;
  int *blockstart;
  int *patchblank;
  unsigned char *ipqqi;
  int iblock;
  blockagedata *bc;
  patch *patchi;
  float *color11;
  mesh *meshblock;
  float burn_color[4]={0.0,0.0,0.0,1.0};
  float clear_color[4]={1.0,1.0,1.0,1.0};

  if(vis_threshold==1&&vis_onlythreshold==1&&do_threshold==1)return;

  patchtimes=meshi->patchtimes;
  visPatches=meshi->visPatches;
  xyzpatch=meshi->xyzpatch;
  patchdir=meshi->patchdir;
  patchrow=meshi->patchrow;
  patchcol=meshi->patchcol;
  blockstart=meshi->blockstart;
  patchblank=meshi->patchblank;
  patchi=patchinfo+meshi->patchfilenum;
  switch(patchi->compression_type){
  case 0:
    ASSERT(meshi->ipqqi!=NULL);
    ipqqi=meshi->ipqqi;
    break;
  case 1:
    ASSERT(meshi->ipqqi_zlib!=NULL);
    ipqqi=meshi->ipqqi_zlib;
    break;
  default:
    ASSERT(0);
  }
  patchi = patchinfo + meshi->patchfilenum;

  if(patchtimes[0]>times[itime]||patchi->display==0)return;
  if(cullfaces==1)glDisable(GL_CULL_FACE);

  /* if a contour boundary does not match a blockage face then draw "both sides" of boundary */

  nn =0;
  glBegin(GL_TRIANGLES);
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(visPatches[n]==1&&patchdir[n]==0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      ipqqcopy = ipqqi + blockstart[n];

      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        ipq1 = ipqqcopy + irow*ncol;
        ipq2 = ipq1 + ncol;
        ip2 = ip1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            color11=clear_color;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol  ])color11=burn_color;

            glColor4fv(color11);
            glVertex3fv(xyzp1);
            glVertex3fv(xyzp1+3);
            glVertex3fv(xyzp2+3);

            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2+3);
            glVertex3fv(xyzp2);
          }
          ipq1++; ipq2++; ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
    }
    nn += patchrow[n]*patchcol[n];
  }
  glEnd();
  if(cullfaces==1)glEnable(GL_CULL_FACE);

  /* if a contour boundary DOES match a blockage face then draw "one sides" of boundary */

  nn=0;
  glBegin(GL_TRIANGLES);
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock=meshi->meshonpatch[n];
    if(iblock!=-1){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(meshi->visPatches[n]==1&&meshi->patchdir[n]>0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      ipqqcopy = ipqqi + blockstart[n];
      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ipq1 = ipqqcopy + irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;

        xyzp2 = xyzp1 + 3*ncol;
        ipq2 = ipq1 + ncol;
        ip2 = ip1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            color11=clear_color;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol  ])color11=burn_color;
           
            glColor4fv(color11);
            glVertex3fv(xyzp1);
            glVertex3fv(xyzp1+3);
            glVertex3fv(xyzp2+3);

            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2+3);
            glVertex3fv(xyzp2);
          }
          ipq1++; ipq2++; ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
    }
    nn += patchrow[n]*patchcol[n];
  }

  /* if a contour boundary DOES match a blockage face then draw "one sides" of boundary */
  nn=0;
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(visPatches[n]==1&&patchdir[n]<0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      ipqqcopy = ipqqi + blockstart[n];
      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        ipq1 = ipqqcopy + irow*ncol;
        ipq2 = ipq1 + ncol;
        ip2 = ip1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            color11=clear_color;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol  ])color11=burn_color;
          
            glColor4fv(color11);
            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2+3);
            glVertex3fv(xyzp1+3);
            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2);
            glVertex3fv(xyzp2+3);
          }
          ipq1++; ipq2++; ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
    }
    nn += patchrow[n]*patchcol[n];
  }
  glEnd();
}

/* ------------------ drawpatch ------------------------ */

void drawpatch(const mesh *meshi){
  int n,nn,nn1,nn2;
  int nrow, ncol, irow, icol;
  float *xyzp1, *xyzp2;
  unsigned char *ipq1, *ipq2;
  unsigned char *ipqqcopy;
  float *xyzpatchcopy;
  int *patchblankcopy, *ip1, *ip2;
  float *patchtimes;
  int *visPatches;
  float *xyzpatch;
  int *patchdir, *patchrow, *patchcol;
  int *blockstart;
  int *patchblank;
  unsigned char *ipqqi;
  int iblock;
  blockagedata *bc;
  patch *patchi;
  float *color11, *color12, *color21, *color22;
  mesh *meshblock;
  float dboundx,dboundy,dboundz;
  float *xplt, *yplt, *zplt;

  if(vis_threshold==1&&vis_onlythreshold==1&&do_threshold==1)return;

  if(hidepatchsurface==0){
    xplt=meshi->xplt;
    yplt=meshi->yplt;
    zplt=meshi->zplt;

    dboundx = (xplt[1]-xplt[0])/10.0;
    dboundy = (yplt[1]-yplt[0])/10.0;
    dboundz = (zplt[1]-zplt[0])/10.0;
  }

  patchtimes=meshi->patchtimes;
  visPatches=meshi->visPatches;
  xyzpatch=meshi->xyzpatch;
  patchdir=meshi->patchdir;
  patchrow=meshi->patchrow;
  patchcol=meshi->patchcol;
  blockstart=meshi->blockstart;
  patchblank=meshi->patchblank;
  patchi=patchinfo+meshi->patchfilenum;
  switch(patchi->compression_type){
  case 0:
    ASSERT(meshi->ipqqi!=NULL);
    ipqqi=meshi->ipqqi;
    break;
  case 1:
    ASSERT(meshi->ipqqi_zlib!=NULL);
    ipqqi=meshi->ipqqi_zlib;
    break;
  default:
    ASSERT(0);
  }
  patchi = patchinfo + meshi->patchfilenum;

  if(patchtimes[0]>times[itime]||patchi->display==0)return;
  if(cullfaces==1)glDisable(GL_CULL_FACE);

  /* if a contour boundary does not match a blockage face then draw "both sides" of boundary */

  if((use_transparency_data==1&&showcolorbarlines==1)||setpatchchopmin==1||setpatchchopmax==1)transparenton();
  nn =0;
  glBegin(GL_TRIANGLES);
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(visPatches[n]==1&&patchdir[n]==0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      ipqqcopy = ipqqi + blockstart[n];
      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        ipq1 = ipqqcopy + irow*ncol;
        ipq2 = ipq1 + ncol;
        ip2 = ip1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            color11 = rgb_patch+4*(*ipq1);
            color12 = rgb_patch+4*(*(ipq1+1));
            color21 = rgb_patch+4*(*ipq2);
            color22 = rgb_patch+4*(*(ipq2+1));
            if(vis_threshold==1&&vis_onlythreshold==0&&do_threshold==1){
              if(meshi->thresholdtime[nn1+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol  ])color11=&char_color[0];
              if(meshi->thresholdtime[nn1+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol+1])color12=&char_color[0];
              if(meshi->thresholdtime[nn2+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol  ])color21=&char_color[0];
              if(meshi->thresholdtime[nn2+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol+1])color22=&char_color[0];
            }
            if(abs(*ipq1-*(ipq2+1))<abs(*(ipq1+1)-*ipq2)){
              glColor4fv(color11); 
              glVertex3fv(xyzp1);

              glColor4fv(color12); 
              glVertex3fv(xyzp1+3);

              glColor4fv(color22); 
              glVertex3fv(xyzp2+3);

              glColor4fv(color11); 
              glVertex3fv(xyzp1);

              glColor4fv(color22); 
              glVertex3fv(xyzp2+3);

              glColor4fv(color21); 
              glVertex3fv(xyzp2);
            }
            else{
              glColor4fv(color11); 
              glVertex3fv(xyzp1);

              glColor4fv(color12); 
              glVertex3fv(xyzp1+3);

              glColor4fv(color21); 
              glVertex3fv(xyzp2);

              glColor4fv(color12); 
              glVertex3fv(xyzp1+3);

              glColor4fv(color22); 
              glVertex3fv(xyzp2+3);

              glColor4fv(color21); 
              glVertex3fv(xyzp2);
            }
          }
          ipq1++; ipq2++; ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
    }
    nn += patchrow[n]*patchcol[n];
  }
  glEnd();
  if(cullfaces==1)glEnable(GL_CULL_FACE);

  /* if a contour boundary DOES match a blockage face then draw "one sides" of boundary */

  nn=0;
  if(hidepatchsurface==1){
    glBegin(GL_TRIANGLES);
  }
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock=meshi->meshonpatch[n];
    if(iblock!=-1){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(meshi->visPatches[n]==1&&meshi->patchdir[n]>0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      ipqqcopy = ipqqi + blockstart[n];
      if(hidepatchsurface==0){
        glPushMatrix();
        switch (meshi->patchdir[n]){
          case 1:
            glTranslatef(dboundx,0.0,0.0);
            break;
          case 2:
            glTranslatef(0.0,-dboundy,0.0);
            break;
          case 3:
            glTranslatef(0.00,0.0,dboundz);
            break;
        }
        glBegin(GL_TRIANGLES);
      }
      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ipq1 = ipqqcopy + irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;

        xyzp2 = xyzp1 + 3*ncol;
        ipq2 = ipq1 + ncol;
        ip2 = ip1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            color11 = rgb_patch+4*(*ipq1);
            color12 = rgb_patch+4*(*(ipq1+1));
            color21 = rgb_patch+4*(*ipq2);
            color22 = rgb_patch+4*(*(ipq2+1));
            if(vis_threshold==1&&vis_onlythreshold==0&&do_threshold==1){
              if(meshi->thresholdtime[nn1+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol  ])color11=&char_color[0];
              if(meshi->thresholdtime[nn1+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol+1])color12=&char_color[0];
              if(meshi->thresholdtime[nn2+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol  ])color21=&char_color[0];
              if(meshi->thresholdtime[nn2+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol+1])color22=&char_color[0];
            }
            if(abs(*ipq1-*(ipq2+1))<abs(*(ipq1+1)-*ipq2)){
              glColor4fv(color11); 
              glVertex3fv(xyzp1);

              glColor4fv(color12); 
              glVertex3fv(xyzp1+3);

              glColor4fv(color22); 
              glVertex3fv(xyzp2+3);

              glColor4fv(color11); 
              glVertex3fv(xyzp1);

              glColor4fv(color22); 
              glVertex3fv(xyzp2+3);

              glColor4fv(color21); 
              glVertex3fv(xyzp2);
            }
            else{
              glColor4fv(color11); 
              glVertex3fv(xyzp1);

              glColor4fv(color12); 
              glVertex3fv(xyzp1+3);

              glColor4fv(color21); 
              glVertex3fv(xyzp2);

              glColor4fv(color12); 
              glVertex3fv(xyzp1+3);

              glColor4fv(color22); 
              glVertex3fv(xyzp2+3);

              glColor4fv(color21); 
              glVertex3fv(xyzp2);
            }
          }
          ipq1++; ipq2++; ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
      if(hidepatchsurface==0){
        glEnd();
        glPopMatrix();
      }
    }
    nn += patchrow[n]*patchcol[n];
  }

  /* if a contour boundary DOES match a blockage face then draw "one sides" of boundary */
  nn=0;
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(visPatches[n]==1&&patchdir[n]<0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      ipqqcopy = ipqqi + blockstart[n];
      if(hidepatchsurface==0){
        glPushMatrix();
        switch (meshi->patchdir[n]){
          case -1:
            glTranslatef(-dboundx,0.0,0.0);
            break;
          case -2:
            glTranslatef(0.0,dboundy,0.0);
            break;
          case -3:
            glTranslatef(0.00,0.0,-dboundz);
            break;
        }
        glBegin(GL_TRIANGLES);
      }
      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        ipq1 = ipqqcopy + irow*ncol;
        ipq2 = ipq1 + ncol;
        ip2 = ip1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            color11 = rgb_patch+4*(*ipq1);
            color12 = rgb_patch+4*(*(ipq1+1));
            color21 = rgb_patch+4*(*ipq2);
            color22 = rgb_patch+4*(*(ipq2+1));
            if(vis_threshold==1&&vis_onlythreshold==0&&do_threshold==1){
              if(meshi->thresholdtime[nn1+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol  ])color11=&char_color[0];
              if(meshi->thresholdtime[nn1+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol+1])color12=&char_color[0];
              if(meshi->thresholdtime[nn2+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol  ])color21=&char_color[0];
              if(meshi->thresholdtime[nn2+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol+1])color22=&char_color[0];
            }
            if(abs(*ipq1-*(ipq2+1))<abs(*(ipq1+1)-*ipq2)){
              glColor4fv(color11); 
              glVertex3fv(xyzp1);

              glColor4fv(color22); 
              glVertex3fv(xyzp2+3);

              glColor4fv(color12); 
              glVertex3fv(xyzp1+3);

              glColor4fv(color11); 
              glVertex3fv(xyzp1);

              glColor4fv(color21); 
              glVertex3fv(xyzp2);

              glColor4fv(color22); 
              glVertex3fv(xyzp2+3);
            }
            else{
              glColor4fv(color11); 
              glVertex3fv(xyzp1);

              glColor4fv(color21); 
              glVertex3fv(xyzp2);

              glColor4fv(color12); 
              glVertex3fv(xyzp1+3);

              glColor4fv(color12); 
              glVertex3fv(xyzp1+3);

              glColor4fv(color21); 
              glVertex3fv(xyzp2);

              glColor4fv(color22); 
              glVertex3fv(xyzp2+3);
            }
          }
          ipq1++; ipq2++; ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
      if(hidepatchsurface==0){
        glEnd();
        glPopMatrix();
      }
    }
    nn += patchrow[n]*patchcol[n];
  }
  if(hidepatchsurface==1){
    glEnd();
  }
  if((use_transparency_data==1&&showcolorbarlines==1)||setpatchchopmin==1||setpatchchopmax==1)transparentoff();
}

/* ------------------ drawpatch_cellcenter ------------------------ */

void drawpatch_cellcenter(const mesh *meshi){
  int n,nn,nn1;
  int nrow, ncol, irow, icol;
  float *xyzp1, *xyzp2;
  unsigned char *ipq1, *ipq2;
  unsigned char *ipqqcopy;
  float *xyzpatchcopy;
  int *patchblankcopy, *ip1, *ip2;
  float *patchtimes;
  int *visPatches;
  float *xyzpatch;
  int *patchdir, *patchrow, *patchcol;
  int *blockstart;
  int *patchblank;
  unsigned char *ipqqi;
  int iblock;
  blockagedata *bc;
  patch *patchi;
  float *color11, *color12, *color21, *color22;
  mesh *meshblock;

  float dboundx,dboundy,dboundz;
  float *xplt, *yplt, *zplt;

  if(vis_threshold==1&&vis_onlythreshold==1&&do_threshold==1)return;

  if(hidepatchsurface==0){
    xplt=meshi->xplt;
    yplt=meshi->yplt;
    zplt=meshi->zplt;

    dboundx = (xplt[1]-xplt[0])/10.0;
    dboundy = (yplt[1]-yplt[0])/10.0;
    dboundz = (zplt[1]-zplt[0])/10.0;
  }

  patchtimes=meshi->patchtimes;
  visPatches=meshi->visPatches;
  xyzpatch=meshi->xyzpatch;
  patchdir=meshi->patchdir;
  patchrow=meshi->patchrow;
  patchcol=meshi->patchcol;
  blockstart=meshi->blockstart;
  patchblank=meshi->patchblank;
  patchi=patchinfo+meshi->patchfilenum;
  switch(patchi->compression_type){
  case 0:
    ASSERT(meshi->ipqqi!=NULL);
    ipqqi=meshi->ipqqi;
    break;
  case 1:
    ASSERT(meshi->ipqqi_zlib!=NULL);
    ipqqi=meshi->ipqqi_zlib;
    break;
  default:
    ASSERT(0);
  }
  patchi = patchinfo + meshi->patchfilenum;

  if(patchtimes[0]>times[itime]||patchi->display==0)return;
  if(cullfaces==1)glDisable(GL_CULL_FACE);

  /* if a contour boundary does not match a blockage face then draw "both sides" of boundary */

  if((use_transparency_data==1&&showcolorbarlines==1)||setpatchchopmin==1||setpatchchopmax==1)transparenton();
  nn =0;
  glBegin(GL_TRIANGLES);
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(visPatches[n]==1&&patchdir[n]==0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      ipqqcopy = ipqqi + blockstart[n];
      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        ipq1 = ipqqcopy + irow*ncol;
        ipq2 = ipq1 + ncol;
        ip2 = ip1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            color11 = rgb_patch+4*(*ipq1);
            if(vis_threshold==1&&vis_onlythreshold==0&&do_threshold==1){
              if(meshi->thresholdtime[nn1+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol  ])color11=&char_color[0];
            }
            glColor4fv(color11); 
            glVertex3fv(xyzp1);
            glVertex3fv(xyzp1+3);
            glVertex3fv(xyzp2+3);

            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2+3);
            glVertex3fv(xyzp2);
          }
          ipq1++; ipq2++; ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
    }
    nn += patchrow[n]*patchcol[n];
  }
  glEnd();
  if(cullfaces==1)glEnable(GL_CULL_FACE);

  /* if a contour boundary DOES match a blockage face then draw "one sides" of boundary */

  nn=0;
  if(hidepatchsurface==1){
    glBegin(GL_TRIANGLES);
  }
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock=meshi->meshonpatch[n];
    if(iblock!=-1){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(meshi->visPatches[n]==1&&meshi->patchdir[n]>0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      ipqqcopy = ipqqi + blockstart[n];
      if(hidepatchsurface==0){
        glPushMatrix();
        switch (meshi->patchdir[n]){
          case 1:
            glTranslatef(dboundx,0.0,0.0);
            break;
          case 2:
            glTranslatef(0.0,-dboundy,0.0);
            break;
          case 3:
            glTranslatef(0.0,0.0,dboundz);
            break;
        }
        glBegin(GL_TRIANGLES);
      }
      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ipq1 = ipqqcopy + irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;

        xyzp2 = xyzp1 + 3*ncol;
        ipq2 = ipq1 + ncol;
        ip2 = ip1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            color11 = rgb_patch+4*(*ipq1);
            if(vis_threshold==1&&vis_onlythreshold==0&&do_threshold==1){
              if(meshi->thresholdtime[nn1+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol  ])color11=&char_color[0];
            }
            glColor4fv(color11); 
            glVertex3fv(xyzp1);
            glVertex3fv(xyzp1+3);
            glVertex3fv(xyzp2+3);

            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2+3);
            glVertex3fv(xyzp2);
          }
          ipq1++; ipq2++; ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
      if(hidepatchsurface==0){
        glEnd();
        glPopMatrix();
      }
    }
    nn += patchrow[n]*patchcol[n];
  }

  /* if a contour boundary DOES match a blockage face then draw "one sides" of boundary */
  nn=0;
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(visPatches[n]==1&&patchdir[n]<0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      ipqqcopy = ipqqi + blockstart[n];
      if(hidepatchsurface==0){
        glPushMatrix();
        switch (meshi->patchdir[n]){
          case -1:
            glTranslatef(-dboundx,0.0,0.0);
            break;
          case -2:
            glTranslatef(0.0,dboundy,0.0);
            break;
          case -3:
            glTranslatef(0.0,0.0,-dboundz);
            break;
        }
        glBegin(GL_TRIANGLES);
      }
      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        ipq1 = ipqqcopy + irow*ncol;
        ipq2 = ipq1 + ncol;
        ip2 = ip1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            color11 = rgb_patch+4*(*ipq1);
            if(vis_threshold==1&&vis_onlythreshold==0&&do_threshold==1){
              if(meshi->thresholdtime[nn1+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol  ])color11=&char_color[0];
            }
            glColor4fv(color11); 
            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2+3);
            glVertex3fv(xyzp1+3);

            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2);
            glVertex3fv(xyzp2+3);
          }
          ipq1++; ipq2++; ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
      if(hidepatchsurface==0){
        glEnd();
        glPopMatrix();
      }
    }
    nn += patchrow[n]*patchcol[n];
  }
  if(hidepatchsurface==1){
    glEnd();
  }
  if((use_transparency_data==1&&showcolorbarlines==1)||setpatchchopmin==1||setpatchchopmax==1)transparentoff();
}

/* ------------------ drawolythreshold ------------------------ */

void drawonlythreshold(const mesh *meshi){
  int n,nn,nn1,nn2;
  int nrow, ncol, irow, icol;
  float *xyzp1, *xyzp2;
  float *xyzpatchcopy;
  int *patchblankcopy, *ip1, *ip2;
  float *patchtimes;
  int *visPatches;
  float *xyzpatch;
  int *patchdir, *patchrow, *patchcol;
  int *blockstart;
  int *patchblank;
  int iblock;
  blockagedata *bc;
  patch *patchi;
  float *color11, *color12, *color21, *color22;
  float *color_black;
  mesh *meshblock;

  if(vis_threshold==0||vis_onlythreshold==0||do_threshold==0)return;

  patchtimes=meshi->patchtimes;
  visPatches=meshi->visPatches;
  xyzpatch=meshi->xyzpatch_threshold;
  patchdir=meshi->patchdir;
  patchrow=meshi->patchrow;
  patchcol=meshi->patchcol;
  blockstart=meshi->blockstart;
  patchblank=meshi->patchblank;
  patchi=patchinfo+meshi->patchfilenum;
  switch(patchi->compression_type){
  case 0:
    ASSERT(meshi->ipqqi!=NULL);
    break;
  case 1:
    ASSERT(meshi->ipqqi_zlib!=NULL);
    break;
  default:
    ASSERT(0);
  }
  patchi = patchinfo + meshi->patchfilenum;

  if(patchtimes[0]>times[itime]||patchi->display==0)return;
  if(cullfaces==1)glDisable(GL_CULL_FACE);

  /* if a contour boundary does not match a blockage face then draw "both sides" of boundary */

  nn =0;
  color_black=&char_color[0];
  glBegin(GL_TRIANGLES);
  glColor4fv(color_black); 
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(visPatches[n]==1&&patchdir[n]==0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        ip2 = ip1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            int nnulls;

            nnulls=4;
            color11 = NULL;
            color12 = NULL;
            color21 = NULL;
            color22 = NULL;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol  ]){
              color11=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn1+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol+1]){
              color12=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn2+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol  ]){
              color21=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn2+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol+1]){
              color22=&char_color[0];
              nnulls--;
            }

            if(nnulls==0){
              glVertex3fv(xyzp1);
              glVertex3fv(xyzp1+3);
              glVertex3fv(xyzp2+3);

              glVertex3fv(xyzp1);
              glVertex3fv(xyzp2+3);
              glVertex3fv(xyzp2);
            }
            if(nnulls==1){
              if(color11!=NULL)glVertex3fv(xyzp1);
              if(color12!=NULL)glVertex3fv(xyzp1+3);
              if(color22!=NULL)glVertex3fv(xyzp2+3);
              if(color21!=NULL)glVertex3fv(xyzp2);
            }
          }
          ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
    }
    nn += patchrow[n]*patchcol[n];
  }
  glEnd();
  if(cullfaces==1)glEnable(GL_CULL_FACE);

  /* if a contour boundary DOES match a blockage face then draw "one sides" of boundary */

  glBegin(GL_TRIANGLES);
  glColor4fv(color_black); 
  nn=0;
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock=meshi->meshonpatch[n];
    if(iblock!=-1){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(meshi->visPatches[n]==1&&meshi->patchdir[n]>0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;

        xyzp2 = xyzp1 + 3*ncol;
        ip2 = ip1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            int nnulls;

            color11 = NULL;
            color12 = NULL;
            color21 = NULL;
            color22 = NULL;
            nnulls=4;

            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol  ]){
              color11=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn1+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol+1]){
              color12=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn2+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol  ]){
              color21=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn2+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol+1]){
              color22=&char_color[0];
              nnulls--;
            }

            if(nnulls==0){
              glVertex3fv(xyzp1);
              glVertex3fv(xyzp1+3);
              glVertex3fv(xyzp2+3);

              glVertex3fv(xyzp1);
              glVertex3fv(xyzp2+3);
              glVertex3fv(xyzp2);
            }
            if(nnulls==1){
              if(color11!=NULL)glVertex3fv(xyzp1);
              if(color12!=NULL)glVertex3fv(xyzp1+3);
              if(color22!=NULL)glVertex3fv(xyzp2+3);
              if(color21!=NULL)glVertex3fv(xyzp2);
            }
          }
          ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
    }
    nn += patchrow[n]*patchcol[n];
  }

  /* if a contour boundary DOES match a blockage face then draw "one sides" of boundary */
  nn=0;
  for(n=0;n<meshi->npatches;n++){
    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itime]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(visPatches[n]==1&&patchdir[n]<0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      for(irow=0;irow<nrow-1;irow++){
        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        ip1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        ip2 = ip1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*ip1==1&&*ip2==1&&*(ip1+1)==1&&*(ip2+1)==1){
            int nnulls;

            color11 = NULL;
            color12 = NULL;
            color21 = NULL;
            color22 = NULL;
            nnulls=4;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol  ]){
              color11=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn1+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn1+icol+1]){
              color12=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn2+icol  ]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol  ]){
              color21=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn2+icol+1]>=0.0&&times[itime]>meshi->thresholdtime[nn2+icol+1]){
              color22=&char_color[0];
              nnulls--;
            }
            if(nnulls==0){
              glVertex3fv(xyzp1);
              glVertex3fv(xyzp2+3);
              glVertex3fv(xyzp1+3);

              glVertex3fv(xyzp1);
              glVertex3fv(xyzp2);
              glVertex3fv(xyzp2+3);
            }
            if(nnulls==1){
              if(color21!=NULL)glVertex3fv(xyzp2);
              if(color22!=NULL)glVertex3fv(xyzp2+3);
              if(color12!=NULL)glVertex3fv(xyzp1+3);
              if(color11!=NULL)glVertex3fv(xyzp1);
            }
          }
          ip1++; ip2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
    }
    nn += patchrow[n]*patchcol[n];
  }
  glEnd();
}

/* ------------------ getpatchfacedir ------------------------ */

int getpatchfacedir(mesh *meshi, int i1, int i2, int j1, int j2, int k1, int k2, 
                    int *blockonpatch, mesh **meshonpatch){
  int i;
  blockagedata *bc;
  mesh *meshblock;
  int imesh;

  *meshonpatch=NULL;
  if(i1==i2){
    for(i=0;i<meshi->nbptrs;i++){
      bc=meshi->blockageinfoptrs[i];
      if(j1==bc->ijk[JMIN]&&j2==bc->ijk[JMAX]&&
         k1==bc->ijk[KMIN]&&k2==bc->ijk[KMAX]){
        if(i1==bc->ijk[IMIN]){
          bc->patchvis[0]=0;
          *blockonpatch=i;
          *meshonpatch=meshi;
        }
        if(i1==bc->ijk[IMAX]){
          bc->patchvis[1]=0;
          *blockonpatch=i;
          *meshonpatch=meshi;
        }
        if(i1==bc->ijk[IMIN])return(-1);
        if(i1==bc->ijk[IMAX])return(1);
      }
    }
    for(imesh=0;imesh<nmeshes;imesh++){
      meshblock = meshinfo + imesh;
      if(meshblock==meshi)continue;
      for(i=0;i<meshblock->nbptrs;i++){
        bc=meshblock->blockageinfoptrs[i];
      }
    }
  }
  else if(j1==j2){
    for(i=0;i<meshi->nbptrs;i++){
      bc=meshi->blockageinfoptrs[i];
      if(i1==bc->ijk[IMIN]&&i2==bc->ijk[IMAX]&&
         k1==bc->ijk[KMIN]&&k2==bc->ijk[KMAX]){
        if(j1==bc->ijk[JMIN]){
          bc->patchvis[2]=0;
          *blockonpatch=i;
          *meshonpatch=meshi;
        }
        if(j1==bc->ijk[JMAX]){
          bc->patchvis[3]=0;
          *blockonpatch=i;
          *meshonpatch=meshi;
        }
        if(j1==bc->ijk[JMIN])return(2);
        if(j1==bc->ijk[JMAX])return(-2);
      }
    }
  }
  else if(k1==k2){
    for(i=0;i<meshi->nbptrs;i++){
      bc=meshi->blockageinfoptrs[i];
      if(i1==bc->ijk[IMIN]&&i2==bc->ijk[IMAX]&&
         j1==bc->ijk[JMIN]&&j2==bc->ijk[JMAX]){
        if(k1==bc->ijk[KMIN]){
          bc->patchvis[4]=0;
          *blockonpatch=i;
          *meshonpatch=meshi;
        }
        if(k1==bc->ijk[KMAX]){
          bc->patchvis[5]=0;
          *blockonpatch=i;
          *meshonpatch=meshi;
        }
        if(k1==bc->ijk[KMIN])return(-3);
        if(k1==bc->ijk[KMAX])return(3);
      }
    }
  }
  *blockonpatch=-1;
  *meshonpatch=NULL;
  if(i1==i2){
	  if(i1==0&&j1==0&&j2==meshi->jbar&&k1==0&&k2==meshi->kbar){
		  return(1);
	  }
	  if(i1==meshi->ibar&&j1==0&&j2==meshi->jbar&&k1==0&&k2==meshi->kbar){
		  return(-1);
	  }
  }
  else if(j1==j2){
	  if(j1==0&&i1==0&&i2==meshi->ibar&&k1==0&&k2==meshi->kbar){
		  return(-1);
	  }
	  if(j1==meshi->jbar&&i1==0&&i2==meshi->ibar&&k1==0&&k2==meshi->kbar){
		  return(1);
	  }
  }
  else if(k1==k2){
	  if(k1==0&&j1==0&&j2==meshi->jbar&&i1==0&&i2==meshi->ibar){
		  return(1);
	  }
	  if(k1==meshi->kbar&&j1==0&&j2==meshi->jbar&&i1==0&&i2==meshi->ibar){
		  return(-1);
	  }
  }
  return(0);
}


/* ------------------ getpatchfacedir ------------------------ */

int getpatchface2dir(mesh *meshi, int i1, int i2, int j1, int j2, int k1, int k2, int patchdir,
                    int *blockonpatch, mesh **meshonpatch){
  int i;
  blockagedata *bc;
  mesh *meshblock;
  int imesh;

  *meshonpatch=NULL;
  if(i1==i2){
    for(i=0;i<meshi->nbptrs;i++){
      bc=meshi->blockageinfoptrs[i];
      if(j1==bc->ijk[JMIN]&&j2==bc->ijk[JMAX]&&
         k1==bc->ijk[KMIN]&&k2==bc->ijk[KMAX]){
        if(i1==bc->ijk[IMIN]&&patchdir==-1){
          bc->patchvis[0]=0;
          *blockonpatch=i;
          *meshonpatch=meshi;
          return(-1);
        }
        if(i1==bc->ijk[IMAX]&&patchdir==1){
          bc->patchvis[1]=0;
          *blockonpatch=i;
          *meshonpatch=meshi;
          return(1);
        }
      }
    }
    for(imesh=0;imesh<nmeshes;imesh++){
      meshblock = meshinfo + imesh;
      if(meshblock==meshi)continue;
      for(i=0;i<meshblock->nbptrs;i++){
        bc=meshblock->blockageinfoptrs[i];
      }
    }
  }
  else if(j1==j2){
    for(i=0;i<meshi->nbptrs;i++){
      bc=meshi->blockageinfoptrs[i];
      if(i1==bc->ijk[IMIN]&&i2==bc->ijk[IMAX]&&
         k1==bc->ijk[KMIN]&&k2==bc->ijk[KMAX]){
        if(j1==bc->ijk[JMIN]&&patchdir==-2){
          bc->patchvis[2]=0;
          *blockonpatch=i;
          *meshonpatch=meshi;
          return(2);
        }
        if(j1==bc->ijk[JMAX]&&patchdir==2){
          bc->patchvis[3]=0;
          *blockonpatch=i;
          *meshonpatch=meshi;
          return(-2);
        }
      }
    }
  }
  else if(k1==k2){
    for(i=0;i<meshi->nbptrs;i++){
      bc=meshi->blockageinfoptrs[i];
      if(i1==bc->ijk[IMIN]&&i2==bc->ijk[IMAX]&&
         j1==bc->ijk[JMIN]&&j2==bc->ijk[JMAX]){
        if(k1==bc->ijk[KMIN]&&patchdir==-3){
          bc->patchvis[4]=0;
          *blockonpatch=i;
          *meshonpatch=meshi;
          return(-3);
        }
        if(k1==bc->ijk[KMAX]&&patchdir==3){
          bc->patchvis[5]=0;
          *blockonpatch=i;
          *meshonpatch=meshi;
          return(3);
        }
      }
    }
  }
  *blockonpatch=-1;
  *meshonpatch=NULL;
  if(i1==i2){
	  if(i1==0&&j1==0&&j2==meshi->jbar&&k1==0&&k2==meshi->kbar){
		  return(1);
	  }
	  if(i1==meshi->ibar&&j1==0&&j2==meshi->jbar&&k1==0&&k2==meshi->kbar){
		  return(-1);
	  }
  }
  else if(j1==j2){
	  if(j1==0&&i1==0&&i2==meshi->ibar&&k1==0&&k2==meshi->kbar){
		  return(-1);
	  }
	  if(j1==meshi->jbar&&i1==0&&i2==meshi->ibar&&k1==0&&k2==meshi->kbar){
		  return(1);
	  }
  }
  else if(k1==k2){
	  if(k1==0&&j1==0&&j2==meshi->jbar&&i1==0&&i2==meshi->ibar){
		  return(1);
	  }
	  if(k1==meshi->kbar&&j1==0&&j2==meshi->jbar&&i1==0&&i2==meshi->ibar){
		  return(-1);
	  }
  }
  return(0);
}

/* ------------------ updatepatchtypes ------------------------ */

void updatepatchtypes(void){
  int i;
  patch *patchi;

  npatchtypes = 0;
  for(i=0;i<npatch_files;i++){
    patchi = patchinfo+i;
    if(getpatchindex(patchi)==-1)patchtypes[npatchtypes++]=i;
  }
  for(i=0;i<npatch_files;i++){
    patchi = patchinfo+i;
    patchi->type=getpatchtype(patchi);
  }
}

/* ------------------ getpatchindex ------------------------ */

int getpatchindex(const patch *patchi){
  patch *patchi2;
  int j;

  for(j=0;j<npatchtypes;j++){
    patchi2 = patchinfo+patchtypes[j];
    if(strcmp(patchi->label.shortlabel,patchi2->label.shortlabel)==0)return patchtypes[j];
  }
  return -1;
}

/* ------------------ getpatchtype ------------------------ */

int getpatchtype(const patch *patchi){
  patch *patchi2;
  int j;

  for(j=0;j<npatchtypes;j++){
    patchi2 = patchinfo+patchtypes[j];
    if(strcmp(patchi->label.shortlabel,patchi2->label.shortlabel)==0)return j;
  }
  return -1;
}

/* ------------------ update_patchtype ------------------------ */

void update_patchtype(){
  int i;
  patch *patchi;

  for(i=0;i<npatch_files;i++){
    patchi = patchinfo + i;
    if(patchi->loaded==1&&patchi->display==1&&patchi->type==ipatchtype)return;
  }

  for(i=0;i<npatch_files;i++){
    patchi = patchinfo + i;
    if(patchi->loaded==1&&patchi->display==1){
      ipatchtype = getpatchindex(patchi);
      return;
    }
  }
  /*
  for(i=0;i<npatch_files;i++){
    patchi = patchinfo + i;
    if(patchi->loaded==1){
      ipatchtype = getpatchindex(patchi);
      patchi->display=1;
    }
  }
  */

  ipatchtype = -1;
  return;
    
}

/* ------------------ patchcompare ------------------------ */

int patchcompare( const void *arg1, const void *arg2 ){
  patch *patchi, *patchj;

  patchi = patchinfo + *(int *)arg1;
  patchj = patchinfo + *(int *)arg2;

  if(strcmp(patchi->label.longlabel,patchj->label.longlabel)<0)return -1;
  if(strcmp(patchi->label.longlabel,patchj->label.longlabel)>0)return 1;
  if(patchi->blocknumber<patchj->blocknumber)return -1;
  if(patchi->blocknumber>patchj->blocknumber)return 1;
  return 0;
}

/* ------------------ updatepatchmenulabels ------------------------ */

void updatepatchmenulabels(void){
  int i;
  patch *patchi;
  char label[128];
  STRUCTSTAT statbuffer;

  if(npatch_files>0){
    FREEMEMORY(patchorderindex);
    NewMemory((void **)&patchorderindex,sizeof(int)*npatch_files);
    for(i=0;i<npatch_files;i++){
      patchorderindex[i]=i;
    }
    qsort( (int *)patchorderindex, (size_t)npatch_files, sizeof(int), patchcompare );

    for(i=0;i<npatch_files;i++){
      patchi = patchinfo + i;
      STRCPY(patchi->menulabel,patchi->label.longlabel);
      if(nmeshes>1){
        sprintf(label,"Mesh %i",1+patchi->blocknumber);
        STRCAT(patchi->menulabel,", ");
        STRCAT(patchi->menulabel,label);
      }
      if(STAT(patchi->comp_file,&statbuffer)==0){
        patchi->file=patchi->comp_file;
        patchi->compression_type=1;
      }
      else{
        patchi->file=patchi->reg_file;
        patchi->compression_type=0;
      }
      if(showfiles==1){
        STRCAT(patchi->menulabel,", ");
        STRCAT(patchi->menulabel,patchi->file);
      }
      if(patchi->compression_type==1){
        STRCAT(patchi->menulabel," (ZLIB)");
      }
    } 
  }


}

/* ------------------ getpatchheader ------------------------ */

void getpatchheader(char *file,int *npatches, float *ppatchmin, float *ppatchmax){
  EGZ_FILE *stream;
  float minmax[2];

  stream=EGZ_FOPEN(file,"rb",0,2);

  if(stream==NULL)return;

  // endian
  // completion (0/1)
  // fileversion (compressed format)
  // version  (bndf version)
  // global min max (used to perform conversion)
  // local min max  (min max found for this file)
  // npatch
  // i1,i2,j1,j2,k1,k2,idir,dummy,dummy (npatch times)
  // time
  // compressed size of frame
  // compressed buffer

  EGZ_FSEEK(stream,16,SEEK_CUR);

  EGZ_FREAD(minmax,4,2,stream);
  EGZ_FSEEK(stream,8,SEEK_CUR);
  EGZ_FREAD(npatches,4,1,stream);
  *ppatchmin=minmax[0];
  *ppatchmax=minmax[1];
  EGZ_FCLOSE(stream);
}

/* ------------------ getpatchheader2 ------------------------ */

void getpatchheader2(char *file,
      int *version,
      int *pi1,int *pi2,
      int *pj1,int *pj2,
      int *pk1,int *pk2,
      int *patchdir){
  int i;
  int buffer[9];
  EGZ_FILE *stream;
  int npatches;

  // endian
  // completion (0/1)
  // fileversion (compressed format)
  // version  (bndf version)
  // global min max (used to perform conversion)
  // local min max  (min max found for this file)
  // npatch
  // i1,i2,j1,j2,k1,k2,idir,dummy,dummy (npatch times)
  // time
  // compressed size of frame
  // compressed buffer

  stream=EGZ_FOPEN(file,"rb",0,2);

  if(stream==NULL)return;

  EGZ_FSEEK(stream,12,SEEK_CUR);
  EGZ_FREAD(version,4,1,stream);
  EGZ_FSEEK(stream,16,SEEK_CUR);
  EGZ_FREAD(&npatches,4,1,stream);
  for(i=0;i<npatches;i++){
    buffer[6]=0;
    if(*version==0){
      EGZ_FREAD(buffer,4,6,stream);
    }
    else{
      EGZ_FREAD(buffer,4,9,stream);
    }
    pi1[i]=buffer[0];
    pi2[i]=buffer[1];
    pj1[i]=buffer[2];
    pj2[i]=buffer[3];
    pk1[i]=buffer[4];
    pk2[i]=buffer[5];
    patchdir[i]=buffer[6];
  }
  EGZ_FCLOSE(stream);
}

/* ------------------ getpatchsizeinfo ------------------------ */

void getpatchsizeinfo(patch *patchi, int *nframes, int *buffersize){
  FILE *streamsize;
  EGZ_FILE *stream;
  int nf=0, bsize=0;
  float local_time;
  char buffer[255];
  int full_size, compressed_size;
  int npatches;
  int buff[9];
  int version;
  int size;
  size_t return_code;
  int i;
  char sizefile[1024];
  int local_count;
  float time_max;

  // endian
  // completion (0/1)
  // fileversion (compressed format)
  // version  (bndf version)
  // global min max (used to perform conversion)
  // local min max  (min max found for this file)
  // npatch
  // i1,i2,j1,j2,k1,k2,idir,dummy,dummy (npatch times)
  // time
  // compressed size of frame
  // compressed buffer

  strcpy(sizefile,patchi->size_file);
  strcat(sizefile,".szz");
  streamsize=fopen(sizefile,"r");
  if(streamsize==NULL){
    *nframes=0;
    *buffersize=0;

    strcpy(sizefile,patchi->size_file);
    strcat(sizefile,".sz");
    streamsize=fopen(sizefile,"r");

    stream=EGZ_FOPEN(patchi->file,"rb",0,2);
    if(stream==NULL)return;

    streamsize=fopen(sizefile,"w");
    if(streamsize==NULL){
      EGZ_FCLOSE(stream);
      return;
    }

    EGZ_FSEEK(stream,12,SEEK_CUR);
    EGZ_FREAD(&version,4,1,stream);
    EGZ_FSEEK(stream,16,SEEK_CUR);
    EGZ_FREAD(&npatches,4,1,stream);
    size=0;
    return_code=0;
    for(i=0;i<npatches;i++){
      if(version==0){
        return_code=EGZ_FREAD(buff,4,6,stream);
      }
      else{
        return_code=EGZ_FREAD(buff,4,9,stream);
      }
      if(return_code==0)break;
      size+=(buff[1]+1-buff[0])*(buff[3]+1-buff[2])*(buff[5]+1-buff[4]);
    }
    if(return_code==0){
      EGZ_FCLOSE(stream);
      fclose(streamsize);
      return;
    }
    for(;;){
      return_code=EGZ_FREAD(&local_time,4,1,stream);
      if(return_code==0)break;
      return_code=EGZ_FREAD(&compressed_size,4,1,stream);
      if(return_code==0)break;
      return_code=EGZ_FSEEK(stream,compressed_size,SEEK_CUR);
      if(return_code!=0)break;
      fprintf(streamsize,"%f %i %i\n",local_time,size,compressed_size);
    }
    EGZ_FCLOSE(stream);
    fclose(streamsize);
    streamsize=fopen(sizefile,"r");
    if(streamsize==NULL)return;
  }

  local_count=-1;
  time_max=-1000000.0;
  for(;;){
    int frame_skip;

    if(fgets(buffer,255,streamsize)==NULL)break;

    sscanf(buffer,"%f %i %i",&local_time,&full_size,&compressed_size);
    frame_skip=1;
    if(local_time>time_max){
      time_max=local_time;
      frame_skip=0;
    }
    if(frame_skip==1)continue;

    local_count++;
    if(local_count%boundframestep!=0)continue;

    nf++;
    bsize+=compressed_size;
  }
  *nframes=nf;
  *buffersize=bsize;
  fclose(streamsize);
}

/* ------------------ getpatchdata_zlib ------------------------ */

void getpatchdata_zlib(patch *patchi,unsigned char *data,int ndata, 
                       float *local_times, unsigned int *zipoffset, unsigned int *zipsize, int ntimes){
  EGZ_FILE *stream;
  float local_time;
  unsigned int compressed_size;
  int npatches;
  int version;
  int return_code;
  int i;
  int local_skip;
  unsigned char *datacopy;
  unsigned int offset;
  int local_count;
  float time_max;

  // endian
  // completion (0/1)
  // fileversion (compressed format)
  // version  (bndf version)
  // global min max (used to perform conversion)
  // local min max  (min max found for this file)
  // npatch
  // i1,i2,j1,j2,k1,k2,idir,dummy,dummy (npatch times)
  // time
  // compressed size of frame
  // compressed buffer

  stream=EGZ_FOPEN(patchi->file,"rb",0,2);
  if(stream==NULL)return;

  EGZ_FSEEK(stream,12,SEEK_CUR);
  EGZ_FREAD(&version,4,1,stream);
  EGZ_FSEEK(stream,16,SEEK_CUR);
  EGZ_FREAD(&npatches,4,1,stream);
  if(version==0){
    local_skip = 6*npatches*4;
  }
  else{
    local_skip = 9*npatches*4;
  }
  return_code=EGZ_FSEEK(stream,local_skip,SEEK_CUR);
  if(return_code!=0){
    EGZ_FCLOSE(stream);
    return;
  }
  datacopy=data;
  offset=0;
  local_count=-1;
  i=-1;
  time_max=-1000000.0;
  CheckMemory;
  for(;;){
    int skip_frame;

    if(EGZ_FREAD(&local_time,4,1,stream)==0)break;
    skip_frame=1;
    if(local_time>time_max){
      time_max=local_time;
      skip_frame=0;
      local_count++;
    }
    if(EGZ_FREAD(&compressed_size,4,1,stream)==0)break;
    if(skip_frame==0&&local_count%boundframestep==0){
      if(EGZ_FREAD(datacopy,1,compressed_size,stream)==0)break;
    }
    else{
      EGZ_FSEEK(stream,compressed_size,SEEK_CUR);
    }
               
    if(skip_frame==1||local_count%boundframestep!=0)continue;
    i++;
    if(i>=ntimes)break;
    printf("boundary time=%.2f\n",local_time);
    ASSERT(i<ntimes);
    local_times[i]=local_time;
    zipoffset[i]=offset;
    zipsize[i]=compressed_size;
    datacopy+=compressed_size;
    offset+=compressed_size;
  }
  EGZ_FCLOSE(stream);
}

/* ------------------ uncompress_patchdataframe ------------------------ */

#ifdef USE_ZLIB
void uncompress_patchdataframe(mesh *meshi,int local_iframe){
  unsigned int countin;
  uLongf countout;
  unsigned char *compressed_data;

  compressed_data = meshi->ipqq_zlib+meshi->zipoffset[local_iframe];
  countin = meshi->zipsize[local_iframe];
  countout=meshi->npatchsize;

  uncompress(meshi->ipqqi_zlib,&countout,compressed_data,countin);

}
#endif

/* ------------------ update_hidepatchsurface ------------------------ */

void update_hidepatchsurface(void){
  int hidepatchsurface_old;
      
  hidepatchsurface_old=hidepatchsurface;
  if(setpatchchopmin==1||setpatchchopmax==1){
    hidepatchsurface=0;
  }
  else{
    hidepatchsurface=1;
  }
  if(hidepatchsurface_old!=hidepatchsurface)updatefacelists=1;
}


/* ------------------ update_patch_hist ------------------------ */

void Update_All_Patch_Bounds(void){
  int i;
  int total=0;

  for(i=0;i<npatch_files;i++){
    patch *patchi;

    patchi = patchinfo + i;
    total+=update_patch_hist(patchi);
    update_patch_bounds(patchi);
  }
  if(total==0){
    printf("Boundary file bounds already computed.\n");
  }
  else{
    printf("Bounds for %i boundary files computed\n",total);
  }
}

/* ------------------ update_patch_hist ------------------------ */

int update_patch_hist(patch *patchj){
  int i;
  int endiandata;
  int first=1;
  int sum=0;

  if(patchj->setvalmax==SET_MAX&&patchj->setvalmin==SET_MIN)return 0;
  endiandata=getendian();

  for(i=0;i<npatch_files;i++){
    int endian, npatches, error, boundaryunitnumber;
    patch *patchi;
    int unit1;
    FILE_SIZE lenfile;
    int error1;
    int *pi1, *pi2, *pj1, *pj2, *pk1, *pk2, *patchdir, *patchsize;
    float patchtime1, *patchframe;
    int patchframesize;
    int j;
    mesh *meshi;

    patchi = patchinfo + i;
    if(patchi->inuse_getbounds==1||patchi->histogram->complete==1||patchi->type!=patchj->type||patchi->cellcenter!=patchj->cellcenter||patchi->bounds.defined==1)continue;

    patchi->inuse_getbounds=1;
    meshi = meshinfo + patchi->blocknumber;

    if(first==1){
      printf("Determining %s percentile and global data bounds\n",patchi->label.longlabel);
      first=0;
    }
    printf("  Examining %s\n",patchi->file);
    sum++;
    lenfile=strlen(patchi->file);

    FORTget_file_unit(&unit1,&patchi->unit_start);
    FORTgetboundaryheader1(patchi->file,&unit1,&endian, &npatches, &error, lenfile);
    if(npatches==0){
      FORTclosefortranfile(&unit1);
      continue;
    }

    NewMemory((void **)&pi1,npatches*sizeof(int));
    NewMemory((void **)&pi2,npatches*sizeof(int));
    NewMemory((void **)&pj1,npatches*sizeof(int));
    NewMemory((void **)&pj2,npatches*sizeof(int));
    NewMemory((void **)&pk1,npatches*sizeof(int));
    NewMemory((void **)&pk2,npatches*sizeof(int));
    NewMemory((void **)&patchdir,npatches*sizeof(int));
    NewMemory((void **)&patchsize,npatches*sizeof(int));

    FORTgetboundaryheader2(&unit1, &patchi->version, &npatches, pi1, pi2, pj1, pj2, pk1, pk2, patchdir);

    patchframesize=0;
    for(j=0;j<npatches;j++){
      int patchsize;

      patchsize = (pi2[j]+1-pi1[j]);
      patchsize *= (pj2[j]+1-pj1[j]);
      patchsize *= (pk2[j]+1-pk1[j]);
      patchframesize+=patchsize;
    }
    
    FORTget_file_unit(&unit1,&patchi->unit_start);
    FORTopenboundary(patchi->file,&unit1,&endiandata,&patchi->version,&error1,lenfile);

    NewMemory((void **)&patchframe,patchframesize*sizeof(float));
    init_histogram(patchi->histogram);
    error1=0;
    while(error1==0){
      int ndummy;

      FORTgetpatchdata(&unit1, &npatches, 
        pi1, pi2, pj1, pj2, pk1, pk2, &patchtime1, patchframe, &ndummy,&error1);
      update_histogram(patchframe,patchframesize,patchi->histogram);
    }
    FORTclosefortranfile(&unit1);
    FREEMEMORY(patchframe);
    FREEMEMORY(pi1);
    FREEMEMORY(pi2);
    FREEMEMORY(pj1);
    FREEMEMORY(pj2);
    FREEMEMORY(pk1);
    FREEMEMORY(pk2);
    FREEMEMORY(patchdir);
    FREEMEMORY(patchsize);
    complete_histogram(patchi->histogram);
  }
  return sum;
}


