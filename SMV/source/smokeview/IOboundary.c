#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include GLUT_H

#include "smv_endian.h"
#include "update.h"
#include "smokeviewvars.h"
#include "compress.h"

/* ------------------ output_Patchdata ------------------------ */

void output_Patchdata(char *csvfile, char *patchfile, meshdata *meshi){
  int iframe;
  float *vals;
  float *xplt, *yplt, *zplt;
  FILE *csvstream=NULL;

  csvstream=fopen(csvfile,"w");
  if(csvstream==NULL)return;
  fprintf(csvstream,"%s\n",patchfile);
  fprintf(csvstream,"time interval:,%f,%f\n",patchout_tmin,patchout_tmax);
  fprintf(csvstream,"region:,%f,%f,%f,%f,%f,%f\n\n",patchout_xmin,patchout_xmax,patchout_ymin,patchout_ymax,patchout_zmin,patchout_zmax);

  vals=meshi->patchval;
  xplt = meshi->xplt_orig;
  yplt = meshi->yplt_orig;
  zplt = meshi->zplt_orig;

  for(iframe=0;iframe<meshi->mxpatch_frames;iframe++){
    int ipatch;
    float pt;

    pt = meshi->patch_times[iframe];
    if(patchout_tmin<patchout_tmax&&(pt<patchout_tmin||pt>patchout_tmax)){
      vals+=meshi->npatchsize;
      continue;
    }

    for(ipatch=0;ipatch<meshi->npatches;ipatch++){
      int i1, i2, j1, j2, k1, k2;
      int imin, imax, jmin, jmax, kmin, kmax;
      int i, j, k;

      i1 = meshi->pi1[ipatch];
      i2 = meshi->pi2[ipatch];
      j1 = meshi->pj1[ipatch];
      j2 = meshi->pj2[ipatch];
      k1 = meshi->pk1[ipatch];
      k2 = meshi->pk2[ipatch];
      if(patchout_xmin<patchout_xmax&&(patchout_xmax<xplt[i1]||patchout_xmin>xplt[i2]))continue;
      if(patchout_ymin<patchout_ymax&&(patchout_ymax<yplt[j1]||patchout_ymin>yplt[j2]))continue;
      if(patchout_zmin<patchout_zmax&&(patchout_zmax<zplt[k1]||patchout_zmin>zplt[k2]))continue;
      imin=i1;
      imax=i2;
      jmin=j1;
      jmax=j2;
      kmin=k1;
      kmax=k2;
      for(i=0;i<meshi->ibar;i++){
        if(xplt[i]<=patchout_xmin&&patchout_xmin<=xplt[i+1])imin=i;
        if(xplt[i]<=patchout_xmax&&patchout_xmax<=xplt[i+1])imax=i;
      }
      for(j=0;j<meshi->jbar;j++){
        if(yplt[j]<=patchout_ymin&&patchout_ymin<=yplt[j+1])jmin=j;
        if(yplt[j]<=patchout_ymax&&patchout_ymax<=yplt[j+1])jmax=j;
      }
      for(k=0;k<meshi->kbar;k++){
        if(zplt[k]<=patchout_zmin&&patchout_zmin<=zplt[k+1])kmin=k;
        if(zplt[k]<=patchout_zmax&&patchout_zmax<=zplt[k+1])kmax=k;
      }

      fprintf(csvstream,"time:,%f,patch %i, of, %i\n",meshi->patch_times[iframe],ipatch+1,meshi->npatches);
      fprintf(csvstream,"region:,%i,%i,%i,%i,%i,%i\n",i1,i2,j1,j2,k1,k2);
      fprintf(csvstream,",%f,%f,%f,%f,%f,%f\n\n",xplt[i1],xplt[i2],yplt[j1],yplt[j2],zplt[k1],zplt[k2]);
      if(i1==i2){
        for(k=k1;k<=k2;k++){
          int out;

          if(k==k1){
            fprintf(csvstream,"Z\\Y,");
            for(j=jmin;j<=jmax;j++){
              fprintf(csvstream,"%f,",yplt[j]);
            }
            fprintf(csvstream,"\n");
          }
          if(k>=kmin&&k<=kmax){
            fprintf(csvstream,"%f,",zplt[k]);
          }

          out=0;
          for(j=j1;j<=j2;j++){
            if(k>=kmin&&k<=kmax&&j>=jmin&&j<=jmax){
              fprintf(csvstream,"%f,",*vals);
              out=1;
            }
            vals++;
          }
          if(out==1)fprintf(csvstream,"\n");
        }
      }
      else if(j1==j2){
        for(k=k1;k<=k2;k++){
          int out;

          if(k==k1){
            fprintf(csvstream,"Z\\X,");
            for(i=imin;i<=imax;i++){
              fprintf(csvstream,"%f,",xplt[i]);
            }
            fprintf(csvstream,"\n");
          }
          if(k>=kmin&&k<=kmax){
            fprintf(csvstream,"%f,",zplt[k]);
          }

          out=0;
          for(i=i1;i<=i2;i++){
            if(k>=kmin&&k<=kmax&&i>=imin&&i<=imax){
              fprintf(csvstream,"%f,",*vals);
              out=1;
            }
            vals++;
          }
          if(out==1)fprintf(csvstream,"\n");
        }
      }
      else{
        for(j=j1;j<j2;j++){
          int out;

          if(j==j1){
            fprintf(csvstream,"Y\\X,");
            for(i=imin;i<=imax;i++){
              fprintf(csvstream,"%f,",xplt[i]);
            }
            fprintf(csvstream,"\n");
          }
          if(j>=jmin&&j<=jmax){
            fprintf(csvstream,"%f,",yplt[j]);
          }

          out=0;
          for(i=i1;i<=i2;i++){
            if(i>=imin&&i<=imax&&j>=jmin&&j<=jmax){
              fprintf(csvstream,"%f,",*vals);
              out=1;
            }
            vals++;
          }
          if(out==1)fprintf(csvstream,"\n");
        }
      }
    }
  }
  fclose(csvstream);

}

/* ------------------ getpatchfacedir ------------------------ */

int getpatchfacedir(meshdata *meshi, int i1, int i2, int j1, int j2, int k1, int k2,
  int *blockonpatch, meshdata **meshonpatch){
  int i;

  *meshonpatch = NULL;
  if(i1 == i2){
    for(i = 0; i < meshi->nbptrs; i++){
      blockagedata *bc;

      bc = meshi->blockageinfoptrs[i];
      if(j1 == bc->ijk[JMIN] && j2 == bc->ijk[JMAX] &&
        k1 == bc->ijk[KMIN] && k2 == bc->ijk[KMAX]){
        if(i1 == bc->ijk[IMIN]){
          bc->patchvis[0] = 0;
          *blockonpatch = i;
          *meshonpatch = meshi;
        }
        if(i1 == bc->ijk[IMAX]){
          bc->patchvis[1] = 0;
          *blockonpatch = i;
          *meshonpatch = meshi;
        }
        if(i1 == bc->ijk[IMIN])return(-1);
        if(i1 == bc->ijk[IMAX])return(1);
      }
    }
  }
  else if(j1 == j2){
    for(i = 0; i < meshi->nbptrs; i++){
      blockagedata *bc;

      bc = meshi->blockageinfoptrs[i];
      if(i1 == bc->ijk[IMIN] && i2 == bc->ijk[IMAX] &&
        k1 == bc->ijk[KMIN] && k2 == bc->ijk[KMAX]){
        if(j1 == bc->ijk[JMIN]){
          bc->patchvis[2] = 0;
          *blockonpatch = i;
          *meshonpatch = meshi;
        }
        if(j1 == bc->ijk[JMAX]){
          bc->patchvis[3] = 0;
          *blockonpatch = i;
          *meshonpatch = meshi;
        }
        if(j1 == bc->ijk[JMIN])return(2);
        if(j1 == bc->ijk[JMAX])return(-2);
      }
    }
  }
  else if(k1 == k2){
    for(i = 0; i < meshi->nbptrs; i++){
      blockagedata *bc;

      bc = meshi->blockageinfoptrs[i];
      if(i1 == bc->ijk[IMIN] && i2 == bc->ijk[IMAX] &&
        j1 == bc->ijk[JMIN] && j2 == bc->ijk[JMAX]){
        if(k1 == bc->ijk[KMIN]){
          bc->patchvis[4] = 0;
          *blockonpatch = i;
          *meshonpatch = meshi;
        }
        if(k1 == bc->ijk[KMAX]){
          bc->patchvis[5] = 0;
          *blockonpatch = i;
          *meshonpatch = meshi;
        }
        if(k1 == bc->ijk[KMIN])return(-3);
        if(k1 == bc->ijk[KMAX])return(3);
      }
    }
  }
  *blockonpatch = -1;
  *meshonpatch = NULL;
  if(i1 == i2){
    if(i1 == 0 && j1 == 0 && j2 == meshi->jbar&&k1 == 0 && k2 == meshi->kbar){
      return(1);
    }
    if(i1 == meshi->ibar&&j1 == 0 && j2 == meshi->jbar&&k1 == 0 && k2 == meshi->kbar){
      return(-1);
    }
  }
  else if(j1 == j2){
    if(j1 == 0 && i1 == 0 && i2 == meshi->ibar&&k1 == 0 && k2 == meshi->kbar){
      return(-1);
    }
    if(j1 == meshi->jbar&&i1 == 0 && i2 == meshi->ibar&&k1 == 0 && k2 == meshi->kbar){
      return(1);
    }
  }
  else if(k1 == k2){
    if(k1 == 0 && j1 == 0 && j2 == meshi->jbar&&i1 == 0 && i2 == meshi->ibar){
      return(1);
    }
    if(k1 == meshi->kbar&&j1 == 0 && j2 == meshi->jbar&&i1 == 0 && i2 == meshi->ibar){
      return(-1);
    }
  }
  return(0);
}

/* ------------------ getpatchface2dir ------------------------ */

int getpatchface2dir(meshdata *meshi, int i1, int i2, int j1, int j2, int k1, int k2, int patchdir,
  int *blockonpatch, meshdata **meshonpatch){
  int i;
  blockagedata *bc;

  *meshonpatch = NULL;
  if(i1 == i2){
    for(i = 0; i < meshi->nbptrs; i++){
      bc = meshi->blockageinfoptrs[i];
      if(j1 == bc->ijk[JMIN] && j2 == bc->ijk[JMAX] &&
        k1 == bc->ijk[KMIN] && k2 == bc->ijk[KMAX]){
        if(i1 == bc->ijk[IMIN] && patchdir == XDIRNEG){
          bc->patchvis[0] = 0;
          *blockonpatch = i;
          *meshonpatch = meshi;
          return(-1);
        }
        if(i1 == bc->ijk[IMAX] && patchdir == XDIR){
          bc->patchvis[1] = 0;
          *blockonpatch = i;
          *meshonpatch = meshi;
          return(1);
        }
      }
    }
  }
  else if(j1 == j2){
    for(i = 0; i < meshi->nbptrs; i++){
      bc = meshi->blockageinfoptrs[i];
      if(i1 == bc->ijk[IMIN] && i2 == bc->ijk[IMAX] &&
        k1 == bc->ijk[KMIN] && k2 == bc->ijk[KMAX]){
        if(j1 == bc->ijk[JMIN] && patchdir == YDIRNEG){
          bc->patchvis[2] = 0;
          *blockonpatch = i;
          *meshonpatch = meshi;
          return(2);
        }
        if(j1 == bc->ijk[JMAX] && patchdir == YDIR){
          bc->patchvis[3] = 0;
          *blockonpatch = i;
          *meshonpatch = meshi;
          return(-2);
        }
      }
    }
  }
  else if(k1 == k2){
    for(i = 0; i < meshi->nbptrs; i++){
      bc = meshi->blockageinfoptrs[i];
      if(i1 == bc->ijk[IMIN] && i2 == bc->ijk[IMAX] &&
        j1 == bc->ijk[JMIN] && j2 == bc->ijk[JMAX]){
        if(k1 == bc->ijk[KMIN] && patchdir == ZDIRNEG){
          bc->patchvis[4] = 0;
          *blockonpatch = i;
          *meshonpatch = meshi;
          return(-3);
        }
        if(k1 == bc->ijk[KMAX] && patchdir == ZDIR){
          bc->patchvis[5] = 0;
          *blockonpatch = i;
          *meshonpatch = meshi;
          return(3);
        }
      }
    }
  }
  *blockonpatch = -1;
  *meshonpatch = NULL;
  if(i1 == i2){
    if(i1 == 0 && j1 == 0 && j2 == meshi->jbar&&k1 == 0 && k2 == meshi->kbar){
      return(1);
    }
    if(i1 == meshi->ibar&&j1 == 0 && j2 == meshi->jbar&&k1 == 0 && k2 == meshi->kbar){
      return(-1);
    }
  }
  else if(j1 == j2){
    if(j1 == 0 && i1 == 0 && i2 == meshi->ibar&&k1 == 0 && k2 == meshi->kbar){
      return(-1);
    }
    if(j1 == meshi->jbar&&i1 == 0 && i2 == meshi->ibar&&k1 == 0 && k2 == meshi->kbar){
      return(1);
    }
  }
  else if(k1 == k2){
    if(k1 == 0 && j1 == 0 && j2 == meshi->jbar&&i1 == 0 && i2 == meshi->ibar){
      return(1);
    }
    if(k1 == meshi->kbar&&j1 == 0 && j2 == meshi->jbar&&i1 == 0 && i2 == meshi->ibar){
      return(-1);
    }
  }
  return(0);
}

/* ------------------ getpatchindex ------------------------ */

int getpatchindex(const patchdata *patchi){
  int j;

  for(j = 0; j < npatchtypes; j++){
    patchdata *patchi2;

    patchi2 = patchinfo + patchtypes[j];
    if(strcmp(patchi->label.shortlabel, patchi2->label.shortlabel) == 0)return patchtypes[j];
  }
  return -1;
}

/* ------------------ readpatch_bndf ------------------------ */

void readpatch_bndf(int ifile, int flag, int *errorcode){
  int error;
  FILE_SIZE lenfile;
  int patchfilenum;
  float *xyzpatchcopy;
  float *xyzpatch_ignitecopy;
  int *patchblankcopy;
  int mxpatch_frames;
  int n;
  int ii;
  int headersize, framesize;
  int statfile;
  STRUCTSTAT statbuffer;
  int nbb;
  int ibartemp,jbartemp,kbartemp;
  float *xplttemp,*yplttemp,*zplttemp;
  int blocknumber;
  patchdata *patchi,*patchbase;
  meshdata *meshi;
  float patchmin_global, patchmax_global;
  int local_first,nsize,iblock;
  int npatchvals;
  char patchcsvfile[1024];

  int nn;
  int filenum;
  char *patchscale;
  int ncompressed_buffer;
  char *file;
  int local_starttime=0, local_stoptime=0;
  FILE_SIZE file_size=0;
  int local_starttime0=0, local_stoptime0=0;
  float delta_time, delta_time0;
  int file_unit;
  int wallcenter=0;

  patchi = patchinfo + ifile;
  if(patchi->loaded==0&&flag==UNLOAD)return;
  if(strcmp(patchi->label.shortlabel,"wc")==0)wallcenter=1;

  if(output_patchdata==1){
      sprintf(patchcsvfile,"%s_bndf_%04i.csv",fdsprefix,ifile);
  }

  local_first=1;
  CheckMemory;
  patchfilenum=ifile;
  file = patchi->file;
  blocknumber = patchi->blocknumber;
  highlight_mesh = blocknumber;
  meshi = meshinfo+blocknumber;
  update_current_mesh(meshi);
  filenum = meshi->patchfilenum;
  if(filenum>=0&&filenum<npatchinfo){
    patchi->loaded=0;
    patchi->display=0;
  }

  meshi->patchfilenum=ifile;
  patchi->display=0;
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
  FREEMEMORY(meshi->patchventcolors);
  FREEMEMORY(meshi->cpatchval);
  FREEMEMORY(meshi->cpatchval_zlib);
  FREEMEMORY(meshi->cpatchval_iframe_zlib);
  FREEMEMORY(meshi->patchval);
  FREEMEMORY(meshi->thresholdtime);
  FREEMEMORY(meshi->patch_times);
  FREEMEMORY(meshi->patchblank);

  if(meshi->patch_contours!=NULL){
    int i;

    ASSERT(meshi->npatches>0&&meshi->mxpatch_frames>0);
    for(i=0;i<meshi->npatches*meshi->mxpatch_frames;i++){
      if(meshi->patch_contours[i]!=NULL){
        freecontour(meshi->patch_contours[i]);
      }
    }
    FREEMEMORY(meshi->patch_contours);
  }

  if(flag==UNLOAD){
    int enableflag=1;
    int i;

    update_patchtype();
    update_unit_defs();
    Update_Times();
    meshi->npatches=0;
    for(i=0;i<npatchinfo;i++){
      patchdata *patchii;

      patchii = patchinfo + i;
      if(patchii->loaded==1&&patchii->compression_type==COMPRESSED_ZLIB){
        enableflag=0;
        break;
      }
    }
    if(enableflag==1)enable_boundary_glui();
    updatemenu=1;
#ifdef pp_MEMPRINT
    PRINTF("After boundary file unload: \n");
    PrintMemoryInfo;
#endif
    return;
  }
  if(ifile>=0&&ifile<npatchinfo){
    global2localpatchbounds(patchi->label.shortlabel);
  }

  if(colorlabelpatch!=NULL){
    for(n=0;n<MAXRGB;n++){
      FREEMEMORY(colorlabelpatch[n]);
    }
    FREEMEMORY(colorlabelpatch);
  }
  patchi->extreme_max=0;
  patchi->extreme_min=0;
  ibartemp=meshi->ibar;
  jbartemp=meshi->jbar;
  kbartemp=meshi->kbar;
  xplttemp=meshi->xplt;
  yplttemp=meshi->yplt;
  zplttemp=meshi->zplt;
  do_threshold=0;

  if(activate_threshold==1){
    if(
      strncmp(patchi->label.shortlabel,"TEMP",4) == 0||
      strncmp(patchi->label.shortlabel,"temp",4) == 0
      ){
      do_threshold=1;
    }
  }

  update_patch_hist(patchi);

  lenfile = strlen(file);
  file_unit=15;
  FORTget_file_unit(&file_unit,&file_unit);
  if(patchi->compression_type==UNCOMPRESSED){
    FILE_SIZE labellen=LABELLEN;
    char patchlonglabel[31], patchshortlabel[31], patchunit[31];

    FORTgetpatchsizes1(&file_unit,file,patchlonglabel,patchshortlabel,patchunit,&meshi->npatches,&headersize,&error,
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
    patchmin_global = patchmin;
    patchmax_global = patchmax;
  }
  if(meshi->npatches>0){
    if(
       NewMemory((void **)&meshi->meshonpatch,sizeof(meshdata *)*meshi->npatches)==0||
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
      if(patchi->compression_type==UNCOMPRESSED){
        FORTclosefortranfile(&file_unit);
      }
      readpatch(ifile,UNLOAD,&error);
      return;
    }
  }

  if(patchi->compression_type==UNCOMPRESSED){
    FORTgetpatchsizes2(&file_unit,&patchi->version,
      &meshi->npatches,&meshi->npatchsize,
      meshi->pi1,meshi->pi2,meshi->pj1,meshi->pj2,meshi->pk1,meshi->pk2,meshi->patchdir,
      &headersize,&framesize);

    // loadpatchbysteps
    //  0 - load entire uncompressed data set
    //  1 - load uncompressed data set one frame at a time
    //  2 - load compressed data set

    loadpatchbysteps=UNCOMPRESSED_ALLFRAMES;
    if(flag==LOAD){
      mxpatch_frames = MAXFRAMES+51;
      statfile=STAT(file,&statbuffer);
      if(statfile==0&&framesize!=0){
        int file_frames;

        file_frames=(statbuffer.st_size-headersize)/framesize+51;
        if(file_frames<mxpatch_frames)mxpatch_frames=file_frames;
      }
      meshi->mxpatch_frames=mxpatch_frames;


  /*
  If the min and max boundary file values are specified then we don't have
  to read in the whole file to determine the bounds.  In this case, memory is allocated
  one time step at a time rather than for all time steps.
  */

      if(statfile==0&&setpatchmin==1&&setpatchmax==1&&cache_boundarydata==0)loadpatchbysteps=UNCOMPRESSED_BYFRAME;
    }
  }
  else{
    int nnsize=0;
    int i;

    getpatchheader2(file,
      &patchi->version,
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
    loadpatchbysteps=COMPRESSED_ALLFRAMES;
  }

  if(meshi->npatchsize>0){
    if(
       NewMemory((void **)&meshi->xyzpatch,3*sizeof(float)*meshi->npatchsize)==0||
       NewMemory((void **)&meshi->xyzpatch_threshold,3*sizeof(float)*meshi->npatchsize)==0||
       NewMemory((void **)&meshi->thresholdtime,sizeof(float)*meshi->npatchsize)==0||
       NewMemory((void **)&meshi->patchblank,meshi->npatchsize*sizeof(int))==0
       ){
      *errorcode=1;
      patchi->loaded=0;
      patchi->display=0;
      if(patchi->compression_type==UNCOMPRESSED){
        FORTclosefortranfile(&file_unit);
      }
      readpatch(ifile,UNLOAD,&error);
      return;
    }
  }
  for(n=0;n<meshi->npatchsize;n++){
    meshi->patchblank[n]=GAS;
  }
  xyzpatchcopy = meshi->xyzpatch;
  xyzpatch_ignitecopy = meshi->xyzpatch_threshold;
  patchblankcopy = meshi->patchblank;
  meshi->blockstart[0]=0;
  for(n=0;n<meshi->nbptrs;n++){
    blockagedata *bc;
    int j;

    bc=meshi->blockageinfoptrs[n];
    for(j=0;j<6;j++){
      bc->patchvis[j]=1;
    }
  }

  for(n=0;n<meshi->npatches;n++){
    float dxx, dyy, dzz, ig_factor;
    float dx_factor, dy_factor, dz_factor;
    int i1, i2, j1, j2, k1, k2;


    i1=meshi->pi1[n];
    i2=meshi->pi2[n];
    j1=meshi->pj1[n];
    j2=meshi->pj2[n];
    k1=meshi->pk1[n];
    k2=meshi->pk2[n];
    if(patchi->version==0){
      meshi->patchdir[n]=getpatchfacedir(meshi,i1,i2,j1,j2,k1,k2,
        meshi->blockonpatch+n,meshi->meshonpatch+n);
    }
    else{
      int patchdir;

      patchdir=meshi->patchdir[n];
      getpatchface2dir(meshi,i1,i2,j1,j2,k1,k2,patchdir,
        meshi->blockonpatch+n,meshi->meshonpatch+n);
      if(meshi->patchdir[n]==YDIR||meshi->patchdir[n]==YDIRNEG){
        meshi->patchdir[n]=-meshi->patchdir[n];
      }
    }
    meshi->patch_surfindex[n]=0;
    dxx = 0.0;
    dyy = 0.0;
    dzz = 0.0;
    ig_factor=0.03;

    switch(meshi->patchdir[n]){
    case XDIRNEG:
      meshi->patch_surfindex[n]=0;
      dxx = -meshi->xplt[1]*ig_factor;
      break;
    case XDIR:
      meshi->patch_surfindex[n]=1;
      dxx = meshi->xplt[1]*ig_factor;
      break;
    case YDIRNEG:
      meshi->patch_surfindex[n]=2;
      dyy = meshi->yplt[1]*ig_factor;
      break;
    case YDIR:
      meshi->patch_surfindex[n]=3;
      dyy = -meshi->yplt[1]*ig_factor;
      break;
    case ZDIRNEG:
      meshi->patch_surfindex[n]=4;
      dzz = -meshi->zplt[1]*ig_factor;
      break;
    case ZDIR:
      meshi->patch_surfindex[n]=5;
      dzz = meshi->zplt[1]*ig_factor;
      break;
    default:
      ASSERT(FFALSE);
    }


    meshi->patchtype[n]=INTERIORwall;
    if(i1==i2){
      int ext_wall;

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
        int k;

        // an internal wall so set blank to 1 then zero out where there are vents
        for(k=k1;k<=k2;k++){
          int j;

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
            *patchblankcopy++ = nodeinvent(meshi,i1,j,k,1,wallcenter);
          }
        }
      }
      else{
        int iii;
        int k;

        // an external wall so set blank to 0 then set to one where there are dummy vents
        iii=0;
        for(k=k1;k<=k2;k++){
          int j;

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
            patchblankcopy[iii++]=SOLID;
          }
        }
        nodein_extvent(n,patchblankcopy,meshi,i1,i2,j1,j2,k1,k2,wallcenter);
        patchblankcopy += (k2+1-k1)*(j2+1-j1);
      }
    }
    else if(j1==j2){
      int ext_wall;

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
        int k;

        for(k=k1;k<=k2;k++){
          int i;

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
            *patchblankcopy++ = nodeinvent(meshi,i,j1,k,2,wallcenter);
          }
        }
      }
      else{
        int iii;
        int k;

        // an external wall so set blank to 0 then zero out where there are vents
        iii=0;
        for(k=k1;k<=k2;k++){
          int i;

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
            patchblankcopy[iii++]=SOLID;
          }
        }
        nodein_extvent(n,patchblankcopy,meshi,i1,i2,j1,j2,k1,k2,wallcenter);
        patchblankcopy += (k2+1-k1)*(i2+1-i1);
      }
    }
    else if(k1==k2){
      int ext_wall;

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
        int j;

        for(j=j1;j<=j2;j++){
          int i;

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
            *patchblankcopy++ = nodeinvent(meshi,i,j,k1,3,wallcenter);
          }
        }
      }
      else{
        int iii;
        int j;

      // an external wall so set blank to 0 then zero out where there are vents

        iii=0;
        for(j=j1;j<=j2;j++){
          int i;

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
            patchblankcopy[iii++]=SOLID;
          }
        }
        nodein_extvent(n,patchblankcopy,meshi,i1,i2,j1,j2,k1,k2,wallcenter);
        patchblankcopy += (j2+1-j1)*(i2+1-i1);
      }
    }
    meshi->blockstart[n+1]=meshi->blockstart[n]+meshi->patchrow[n]*meshi->patchcol[n];
    meshi->visPatches[n]=visPatchType[meshi->patchtype[n]];
  }

  for(n=0;n<meshi->nbptrs;n++){
    blockagedata *bc;
    int j;

    bc=meshi->blockageinfoptrs[n];
    bc->patchvis[6]=0;
    for(j=0;j<6;j++){
      bc->patchvis[6]+=bc->patchvis[j];
    }
    if(bc->patchvis[6]!=0){
      bc->patchvis[6]=1;
    }
  }

  meshi->patchval = NULL;
  switch(loadpatchbysteps){
  case UNCOMPRESSED_ALLFRAMES:
    while(meshi->patchval==NULL&&mxpatch_frames>100){
      mxpatch_frames-=50;
      meshi->mxpatch_frames=mxpatch_frames;
      NewMemory((void **)&meshi->patchval,sizeof(float)*mxpatch_frames*meshi->npatchsize);
    }
    if(meshi->patchval==NULL){
      NewMemory((void **)&meshi->patchval,sizeof(float)*mxpatch_frames*meshi->npatchsize);
    }
    break;
  case UNCOMPRESSED_BYFRAME:
    npatchvals = meshi->npatchsize*meshi->mxpatch_frames;
    if(
      NewMemory((void **)&meshi->patchval,sizeof(float)*meshi->npatchsize)==0||
      NewMemory((void **)&meshi->cpatchval,sizeof(unsigned char)*npatchvals)==0){
      *errorcode=1;
      FORTclosefortranfile(&file_unit);
      readpatch(ifile,UNLOAD,&error);
      return;
    }
    break;
  case COMPRESSED_ALLFRAMES:
    getpatchsizeinfo(patchi, &mxpatch_frames, &ncompressed_buffer);
    NewMemory((void **)&meshi->cpatchval_zlib,sizeof(unsigned char)*ncompressed_buffer);
    NewMemory((void **)&meshi->cpatchval_iframe_zlib,sizeof(unsigned char)*meshi->npatchsize);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }

  NewMemory((void **)&meshi->patch_times,sizeof(float)*mxpatch_frames);
  NewMemory((void **)&meshi->zipoffset,sizeof(unsigned int)*mxpatch_frames);
  NewMemory((void **)&meshi->zipsize,sizeof(unsigned int)*mxpatch_frames);
  if(meshi->patch_times==NULL){
    *errorcode=1;
    FORTclosefortranfile(&file_unit);
    readpatch(ifile,UNLOAD,&error);
    return;
  }
  if(loadpatchbysteps==COMPRESSED_ALLFRAMES){
    getpatchdata_zlib(patchi,meshi->cpatchval_zlib,ncompressed_buffer,
      meshi->patch_times,meshi->zipoffset,meshi->zipsize,mxpatch_frames);
    meshi->npatch_times=mxpatch_frames;
  }
  else{
    if(meshi->patchval==NULL){
      *errorcode=1;
      FORTclosefortranfile(&file_unit);
      readpatch(ifile,UNLOAD,&error);
      return;
    }
    meshi->npatch_times=0;
  }

  file_size=get_filesize(file);
  local_starttime = glutGet(GLUT_ELAPSED_TIME);
  for (ii=0;ii<mxpatch_frames;){
    if(loadpatchbysteps==UNCOMPRESSED_BYFRAME){
      meshi->patchval_iframe = meshi->patchval;
      meshi->cpatchval_iframe = meshi->cpatchval + ii*meshi->npatchsize;
    }
    else if(loadpatchbysteps==UNCOMPRESSED_ALLFRAMES){
      meshi->patchval_iframe = meshi->patchval + ii*meshi->npatchsize;
    }
    meshi->patch_timesi = meshi->patch_times + ii;

    error=0;
    if(loadpatchbysteps==UNCOMPRESSED_ALLFRAMES||loadpatchbysteps==UNCOMPRESSED_BYFRAME){
      for (n=0;n<boundframestep;n++){
        if(error==0){
          int npatchval_iframe;

          FORTgetpatchdata(&file_unit,&meshi->npatches,
          meshi->pi1,meshi->pi2,
          meshi->pj1,meshi->pj2,
          meshi->pk1,meshi->pk2,
          meshi->patch_timesi,meshi->patchval_iframe,&npatchval_iframe,&error);
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

        nn=0;
        if(loadpatchbysteps==COMPRESSED_ALLFRAMES)uncompress_patchdataframe(meshi,ii);
        for(n=0;n<meshi->npatches;n++){
          meshdata *meshblock;
          float dval;
          int j;

          iblock=meshi->blockonpatch[n];
          meshblock = meshi->meshonpatch[n];
          nsize=meshi->patchrow[n]*meshi->patchcol[n];
          ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
          if(iblock!=-1&&meshblock!=NULL){
            switch(loadpatchbysteps){
            case UNCOMPRESSED_ALLFRAMES:
            case UNCOMPRESSED_BYFRAME:
              for(j=0;j<nsize;j++){
                if(meshi->thresholdtime[nn+j]<0.0&&meshi->patchval_iframe[nn+j]>=temp_threshold){
                  meshi->thresholdtime[nn+j]=meshi->patch_times[ii];
                }
              }
              break;
            case COMPRESSED_ALLFRAMES:
              dval = (patchmax-patchmin)/255.0;
              for(j=0;j<nsize;j++){
                float val;
                int ival;

                ival = meshi->cpatchval_iframe_zlib[nn+j];
                val = patchmin + dval*ival;
                if(meshi->thresholdtime[nn+j]<0.0&&val>=temp_threshold){
                  meshi->thresholdtime[nn+j]=meshi->patch_times[ii];
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
    if(loadpatchbysteps==UNCOMPRESSED_BYFRAME){
      getBoundaryColors2(
        meshi->patchval_iframe, meshi->npatchsize, meshi->cpatchval_iframe,
                 setpatchmin,&patchmin, setpatchmax,&patchmax,
                 &patchmin_global, &patchmax_global,
                 nrgb_full,
                 &patchi->extreme_min,&patchi->extreme_max);
    }
    CheckMemory;
    if(error!=0)break;
    if(settmax_b!=0&&*meshi->patch_timesi>tmax_b)break;

    switch(loadpatchbysteps){
      case UNCOMPRESSED_ALLFRAMES:
      case UNCOMPRESSED_BYFRAME:
        if(!(settmin_b!=0&&*meshi->patch_timesi<tmin_b)){
          PRINTF("boundary time=%.2f\n",*meshi->patch_timesi);

          meshi->npatch_times++;
          if(meshi->npatch_times + 1 > mxpatch_frames){
            PRINTF("reallocating memory\n");
            mxpatch_frames = meshi->npatch_times + 50; /* this + 50 must match - 50 below */
            meshi->mxpatch_frames=mxpatch_frames;
            if(
              ResizeMemory((void **)&meshi->patchval,           mxpatch_frames*meshi->npatchsize*sizeof(float))==0||
              ResizeMemory((void **)&meshi->patch_times,    mxpatch_frames*sizeof(float))==0
             ){
              *errorcode=1;
              readpatch(ifile,UNLOAD,&error);
              FORTclosefortranfile(&file_unit);
              return;
            }
          }
          ii++;
        }
        break;
      case COMPRESSED_ALLFRAMES:
        ii++;
        break;
    }
  }
  local_stoptime = glutGet(GLUT_ELAPSED_TIME);
  delta_time = (local_stoptime-local_starttime)/1000.0;

  /* convert patch values into integers pointing to an rgb color table */

  if(loadpatchbysteps==UNCOMPRESSED_ALLFRAMES){
    if(output_patchdata==1){
      output_Patchdata(patchcsvfile,patchi->file,meshi);
    }
    npatchvals = meshi->npatch_times*meshi->npatchsize;
    if(npatchvals==0||NewMemory((void **)&meshi->cpatchval,sizeof(unsigned char)*npatchvals)==0){
      *errorcode=1;
      FORTclosefortranfile(&file_unit);
      readpatch(ifile,UNLOAD,&error);
      return;
    }
  }

  PRINTF("computing boundary color levels \n");
  if(NewMemory((void **)&colorlabelpatch,MAXRGB*sizeof(char *))==0){
    *errorcode=1;
    if(loadpatchbysteps!=COMPRESSED_ALLFRAMES){
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
      if(loadpatchbysteps!=COMPRESSED_ALLFRAMES){
        FORTclosefortranfile(&file_unit);
      }
      readpatch(ifile,UNLOAD,&error);
      return;
    }
  }
  patchscale = patchi->scale;
  patchbase = patchinfo + getpatchindex(patchi);
  patchi->loaded=1;
  ipatchtype=getpatchtype(patchi);
  switch(loadpatchbysteps){
  case UNCOMPRESSED_ALLFRAMES:
    getBoundaryColors3(patchi,meshi->patchval, npatchvals, meshi->cpatchval,
      setpatchmin,&patchmin, setpatchmax,&patchmax,
      &patchmin_global, &patchmax_global,
      nrgb, colorlabelpatch,patchscale,boundarylevels256,
      &patchi->extreme_min,&patchi->extreme_max);
    break;
  case UNCOMPRESSED_BYFRAME:
    getBoundaryLabels(
      patchmin, patchmax,
      colorlabelpatch,patchscale,boundarylevels256,nrgb);
    break;
  case COMPRESSED_ALLFRAMES:
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

  local2globalpatchbounds(patchi->label.shortlabel);
  updatepatchlistindex(patchfilenum);

  if(wallcenter==1){
    int i;

    init_vent_colors();
    FREEMEMORY(meshi->patchventcolors);
    NewMemory((void **)&meshi->patchventcolors,sizeof(float *)*npatchvals);
    for(i=0;i<npatchvals;i++){
      int vent_index;

      vent_index = CLAMP(meshi->patchval[i]+0.1,0,nventcolors-1);
      meshi->patchventcolors[i]=ventcolors[vent_index];
    }
  }
  if(cache_boundarydata==0){
    FREEMEMORY(meshi->patchval);
  }
  patchi->loaded=1;
  patchi->display=1;
  ipatchtype=getpatchtype(patchi);
  showexterior=1-showexterior;
  allexterior = 1-allexterior;
  ShowPatchMenu(EXTERIORwallmenu);
  plotstate=getplotstate(DYNAMIC_PLOTS);
  if(patchi->compression_type==COMPRESSED_ZLIB)disable_boundary_glui();
  Update_Times();
  update_unit_defs();
  updatechopcolors();
#ifdef pp_MEMPRINT
  PRINTF("After boundary file load: \n");
  PrintMemoryInfo;
#endif
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

  glutPostRedisplay();
}

/* ------------------ readpatch ------------------------ */

void readpatch(int ifile, int load_flag, int *errorcode){
  patchdata *patchi;

  patchi = patchinfo + ifile;
  if(patchi->filetype==PATCH_GEOMETRY){
    ASSERT(ifile>=0&&ifile<ngeominfo);
    read_geomdata(ifile,load_flag,errorcode);
  }
  else{
    ASSERT(ifile>=0&&ifile<npatchinfo);
    readpatch_bndf(ifile,load_flag,errorcode);
  }
}

/* ------------------ nodeinblockage ------------------------ */

int nodeinblockage(const meshdata *meshnode, int i,int j,int k, int *imesh, int *iblockage){
  int ii;
  float xn, yn, zn;

  xn = meshnode->xplt[i];
  yn = meshnode->yplt[j];
  zn = meshnode->zplt[k];

  *imesh=-1;

  for(ii=0;ii<nmeshes;ii++){
    int jj;
    meshdata *meshii;
    blockagedata *bc;
    float xm_min, xm_max;
    float ym_min, ym_max;
    float zm_min, zm_max;
    float xb_min, xb_max;
    float yb_min, yb_max;
    float zb_min, zb_max;
    float *xplt, *yplt, *zplt;

    meshii = meshinfo + ii;
    if(meshnode==meshii)continue;

    xplt = meshii->xplt;
    yplt = meshii->yplt;
    zplt = meshii->zplt;

    xm_min=xplt[0];
    xm_max=meshii->xyz_bar[XXX];
    ym_min=yplt[0];
    ym_max=meshii->xyz_bar[YYY];
    zm_min=zplt[0];
    zm_max=meshii->xyz_bar[ZZZ];
    if(xn<xm_min||xn>xm_max)continue;
    if(yn<ym_min||yn>ym_max)continue;
    if(zn<zm_min||zn>zm_max)continue;


    for(jj=0;jj<meshii->nbptrs;jj++){
      bc = meshii->blockageinfoptrs[jj];
      if(bc->hole==1)continue;
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

int nodeinvent(const meshdata *meshi, int i,int j,int k, int dir,int option){
  int ii;

  if(option==1)return 1;
  for(ii=0;ii<meshi->nvents;ii++){
    ventdata *vi;

    vi = meshi->ventinfo+ii;
    if(vi->hideboundary==1){
      switch(dir){
      case XDIR:
        if(vi->imin==i&&i==vi->imax&&
           vi->jmin <j&&j <vi->jmax&&
           vi->kmin <k&&k <vi->kmax){
          return 0;
        }
        break;
      case YDIR:
        if(vi->jmin==j&&j==vi->jmax&&
           vi->imin <i&&i <vi->imax&&
           vi->kmin <k&&k <vi->kmax){
          return 0;
        }
        break;
      case ZDIR:
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

/* ------------------ nodein_extvent ------------------------ */

void nodein_extvent(int ipatch, int *patchblank, const meshdata *meshi,
                    int i1, int i2, int j1, int j2, int k1, int k2, int option){
  int ii, dir=0;

  if(i1==i2)dir=1;
  if(j1==j2)dir=2;
  if(k1==k2)dir=3;

  for(ii=0;ii<meshi->nvents;ii++){
    ventdata *vi;
    int imin, jmin, kmin, imax, jmax, kmax;

    vi = meshi->ventinfo+ii;
    if(vi->hideboundary==1&&option==0)continue;
    if(vi->dir2!=dir)continue;
    switch(dir){
      int i, j, k;

    case XDIR:
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
          int iii;

          iii=(k-k1)*(j2+1-j1) + (j-j1);
          patchblank[iii]=GAS;
        }
      }
      break;
    case YDIR:
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
          int iii;

          iii=(k-k1)*(i2+1-i1) + (i-i1);
          patchblank[iii]=GAS;
        }
      }
      break;
    case ZDIR:
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
          int iii;

          iii=(j-j1)*(i2+1-i1) + (i-i1);
          patchblank[iii]=GAS;
        }
      }
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
  }
  switch(dir){
    int i, j, k;

  case XDIR:
    for(k=k1;k<=k2;k++){
      for(j=j1;j<=j2;j++){
        int iii,imesh,iblockage;

        iii=(k-k1)*(j2+1-j1) + (j-j1);
        if(patchblank[iii]==GAS)continue;
        patchblank[iii] = nodeinblockage(meshi,i1,j,k,&imesh,&iblockage);
      }
    }
    break;
  case YDIR:
    for(k=k1;k<=k2;k++){
      for(i=i1;i<=i2;i++){
        int iii,imesh,iblockage;
        meshdata *meshblock;

        iii=(k-k1)*(i2+1-i1) + (i-i1);
        if(patchblank[iii]==GAS)continue;
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
  case ZDIR:
    for(j=j1;j<=j2;j++){
      for(i=i1;i<=i2;i++){
        int iii,imesh,iblockage;

        iii=(j-j1)*(i2+1-i1) + (i-i1);
        if(patchblank[iii]==GAS)continue;
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
  int i;

  for(i=0;i<nmeshes;i++){
    meshdata *meshi;
    int n;

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

  for(i=0;i<npatchinfo;i++){
    patchdata *patchi;

    patchi = patchinfo + i;
    if(strcmp(patchi->label.shortlabel,key)==0){
      patchi->valmin=patchmin;
      patchi->valmax=patchmax;
      patchi->setvalmin=setpatchmin;
      patchi->setvalmax=setpatchmax;

      patchi->chopmin=patchchopmin;
      patchi->chopmax=patchchopmax;
      patchi->setchopmin=setpatchchopmin;
      patchi->setchopmax=setpatchchopmax;
    }
  }
}

/* ------------------ global2localpatchbounds ------------------------ */

void global2localpatchbounds(const char *key){
  int i;

  for(i=0;i<npatchinfo;i++){
    patchdata *patchi;

    patchi = patchinfo + i;
    if(strcmp(patchi->label.shortlabel,key)==0){
      patchmin=patchi->valmin;
      patchmax=patchi->valmax;
      setpatchmin=patchi->setvalmin;
      setpatchmax=patchi->setvalmax;

      patchchopmin=patchi->chopmin;
      patchchopmax=patchi->chopmax;
      setpatchchopmin=patchi->setchopmin;
      setpatchchopmax=patchi->setchopmax;
      patchmin_unit = (unsigned char *)patchi->label.unit;
      patchmax_unit = patchmin_unit;
      update_glui_patch_units();
      update_hidepatchsurface();

      local2globalpatchbounds(key);
      return;
    }
  }
}

/* ------------------ drawpatch_texture ------------------------ */

void drawpatch_texture(const meshdata *meshi){
  float r11, r12, r21, r22;
  int n;
  int nrow, ncol, irow, icol;
  unsigned char *cpatchval_iframe_copy;
  float *xyzpatchcopy;
  int *patchblankcopy;
  float *patch_times;
  int *visPatches;
  float *xyzpatch;
  int *patchdir, *patchrow, *patchcol, *patchtype;
  int *blockstart;
  int *patchblank;
  unsigned char *cpatchval_iframe;
  int iblock;
  blockagedata *bc;
  patchdata *patchi;
  meshdata *meshblock;
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

  patch_times=meshi->patch_times;
  visPatches=meshi->visPatches;
  xyzpatch=meshi->xyzpatch;
  patchdir=meshi->patchdir;
  patchtype=meshi->patchtype;
  patchrow=meshi->patchrow;
  patchcol=meshi->patchcol;
  blockstart=meshi->blockstart;
  patchblank=meshi->patchblank;
  patchi=patchinfo+meshi->patchfilenum;
  switch(patchi->compression_type){
  case UNCOMPRESSED:
    ASSERT(meshi->cpatchval_iframe!=NULL);
    cpatchval_iframe=meshi->cpatchval_iframe;
    break;
  case COMPRESSED_ZLIB:
    ASSERT(meshi->cpatchval_iframe_zlib!=NULL);
    cpatchval_iframe=meshi->cpatchval_iframe_zlib;
    break;
  default:
    ASSERT(FFALSE);
  }
  patchi = patchinfo + meshi->patchfilenum;

  if(patch_times[0]>global_times[itimes]||patchi->display==0)return;
  if(cullfaces==1)glDisable(GL_CULL_FACE);

  /* if a contour boundary does not match a blockage face then draw "both sides" of boundary */

  if((use_transparency_data==1&&contour_type==LINE_CONTOURS)||setpatchchopmin==1||setpatchchopmax==1)transparenton();
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
  glEnable(GL_TEXTURE_1D);
  glBindTexture(GL_TEXTURE_1D,texture_patch_colorbar_id);

  glBegin(GL_TRIANGLES);
  for(n=0;n<meshi->npatches;n++){
    int drawit;

    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
        continue;
      }
    }
    drawit=0;
    if(visPatches[n]==1&&patchdir[n]==0)drawit=1;
    if(patchtype[n]!=INTERIORwall&&showpatch_both==1)drawit=1;
    if(drawit==1){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      cpatchval_iframe_copy = cpatchval_iframe + blockstart[n];

      for(irow=0;irow<nrow-1;irow++){
        unsigned char *cpatchval1, *cpatchval2;
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        patchblank1 = patchblankcopy + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        cpatchval1 = cpatchval_iframe_copy + irow*ncol;
        cpatchval2 = cpatchval1 + ncol;
        patchblank2 = patchblank1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            r11 = (float)(*cpatchval1)/255.0;
            r12 = (float)(*(cpatchval1+1))/255.0;
            r21 = (float)(*cpatchval2)/255.0;
            r22 = (float)(*(cpatchval2+1))/255.0;
            if(ABS(r11-r22)<ABS(r12-r21)){
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
          cpatchval1++; cpatchval2++; patchblank1++; patchblank2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
    }
  }
  glEnd();
  if(cullfaces==1)glEnable(GL_CULL_FACE);

  /* if a contour boundary DOES match a blockage face then draw "one sides" of boundary */

  if(hidepatchsurface==1){
    glBegin(GL_TRIANGLES);
  }
  for(n=0;n<meshi->npatches;n++){
    int drawit;

    iblock = meshi->blockonpatch[n];
    meshblock=meshi->meshonpatch[n];
    if(iblock!=-1){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
        continue;
      }
    }
    drawit=0;
    if(meshi->visPatches[n]==1&&meshi->patchdir[n]>0){
      if(patchtype[n]==INTERIORwall||showpatch_both==0){
        drawit=1;
      }
    }
    if(drawit==1){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      cpatchval_iframe_copy = cpatchval_iframe + blockstart[n];
      if(hidepatchsurface==0){
        glPushMatrix();
        switch(meshi->patchdir[n]){
          case XDIR:
            glTranslatef(dboundx,0.0,0.0);
            break;
          case YDIR:
            glTranslatef(0.0,-dboundy,0.0);
            break;
          case ZDIR:
            glTranslatef(0.00,0.0,dboundz);
            break;
          default:
            ASSERT(FFALSE);
            break;
        }
        glBegin(GL_TRIANGLES);
      }
      for(irow=0;irow<nrow-1;irow++){
        unsigned char *cpatchval1, *cpatchval2;
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        cpatchval1 = cpatchval_iframe_copy + irow*ncol;
        patchblank1 = patchblankcopy + irow*ncol;

        xyzp2 = xyzp1 + 3*ncol;
        cpatchval2 = cpatchval1 + ncol;
        patchblank2 = patchblank1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            r11 = (float)(*cpatchval1)/255.0;
            r12 = (float)(*(cpatchval1+1))/255.0;
            r21 = (float)(*cpatchval2)/255.0;
            r22 = (float)(*(cpatchval2+1))/255.0;
            if(ABS(*cpatchval1-*(cpatchval2+1))<ABS(*(cpatchval1+1)-*cpatchval2)){
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
          cpatchval1++; cpatchval2++; patchblank1++; patchblank2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
      if(hidepatchsurface==0){
        glEnd();
        glPopMatrix();
      }
    }
  }

  /* if a contour boundary DOES match a blockage face then draw "one sides" of boundary */
  for(n=0;n<meshi->npatches;n++){
    int drawit;

    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
        continue;
      }
    }
    drawit=0;
    if(visPatches[n]==1&&patchdir[n]<0){
      if(patchtype[n]==INTERIORwall||showpatch_both==0){
        drawit=1;
      }
    }
    if(drawit==1){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      cpatchval_iframe_copy = cpatchval_iframe + blockstart[n];
      if(hidepatchsurface==0){
        glPushMatrix();
        switch(meshi->patchdir[n]){
          case XDIRNEG:
            glTranslatef(-dboundx,0.0,0.0);
            break;
          case YDIRNEG:
            glTranslatef(0.0,dboundy,0.0);
            break;
          case ZDIRNEG:
            glTranslatef(0.0,0.0,-dboundz);
            break;
          default:
            ASSERT(FFALSE);
            break;
        }
        glBegin(GL_TRIANGLES);
      }
      for(irow=0;irow<nrow-1;irow++){
        unsigned char *cpatchval1, *cpatchval2;
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        patchblank1 = patchblankcopy + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        cpatchval1 = cpatchval_iframe_copy + irow*ncol;
        cpatchval2 = cpatchval1 + ncol;
        patchblank2 = patchblank1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            r11 = (float)(*cpatchval1)/255.0;
            r12 = (float)(*(cpatchval1+1))/255.0;
            r21 = (float)(*cpatchval2)/255.0;
            r22 = (float)(*(cpatchval2+1))/255.0;
            if(ABS(*cpatchval1-*(cpatchval2+1))<ABS(*(cpatchval1+1)-*cpatchval2)){
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
          cpatchval1++; cpatchval2++; patchblank1++; patchblank2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
      if(hidepatchsurface==0){
        glEnd();
        glPopMatrix();
      }
    }
  }
  if(hidepatchsurface==1){
    glEnd();
  }
  glDisable(GL_TEXTURE_1D);
  if((use_transparency_data==1&&contour_type==LINE_CONTOURS)||setpatchchopmin==1||setpatchchopmax==1)transparentoff();
}

/* ------------------ drawpatch_texture_threshold ------------------------ */

void drawpatch_texture_threshold(const meshdata *meshi){
  float r11, r12, r21, r22;
  int n,nn,nn1,nn2;
  int nrow, ncol, irow, icol;
  unsigned char *cpatchval1, *cpatchval2;
  unsigned char *cpatchval_iframe_copy;
  float *xyzpatchcopy;
  int *patchblankcopy;
  float *patch_times;
  int *visPatches;
  float *xyzpatch;
  int *patchdir, *patchrow, *patchcol;
  int *blockstart;
  int *patchblank;
  unsigned char *cpatchval_iframe;
  int iblock;
  blockagedata *bc;
  patchdata *patchi;
  float *color11, *color12, *color21, *color22;
  meshdata *meshblock;
  float burn_color[4]={0.0,0.0,0.0,1.0};
  float clear_color[4]={1.0,1.0,1.0,1.0};

  if(vis_threshold==1&&vis_onlythreshold==1&&do_threshold==1)return;

  patch_times=meshi->patch_times;
  visPatches=meshi->visPatches;
  xyzpatch=meshi->xyzpatch;
  patchdir=meshi->patchdir;
  patchrow=meshi->patchrow;
  patchcol=meshi->patchcol;
  blockstart=meshi->blockstart;
  patchblank=meshi->patchblank;
  patchi=patchinfo+meshi->patchfilenum;
  switch(patchi->compression_type){
  case UNCOMPRESSED:
    ASSERT(meshi->cpatchval_iframe!=NULL);
    cpatchval_iframe=meshi->cpatchval_iframe;
    break;
  case COMPRESSED_ZLIB:
    ASSERT(meshi->cpatchval_iframe_zlib!=NULL);
    cpatchval_iframe=meshi->cpatchval_iframe_zlib;
    break;
  default:
    ASSERT(FFALSE);
  }
  patchi = patchinfo + meshi->patchfilenum;

  if(patch_times[0]>global_times[itimes]||patchi->display==0)return;
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
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(visPatches[n]==1&&patchdir[n]==0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      cpatchval_iframe_copy = cpatchval_iframe + blockstart[n];

      for(irow=0;irow<nrow-1;irow++){
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        patchblank1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        cpatchval1 = cpatchval_iframe_copy + irow*ncol;
        cpatchval2 = cpatchval1 + ncol;
        patchblank2 = patchblank1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            r11 = (float)(*cpatchval1)/255.0;
            r12 = (float)(*(cpatchval1+1))/255.0;
            r21 = (float)(*cpatchval2)/255.0;
            r22 = (float)(*(cpatchval2+1))/255.0;
            color11=clear_color;
            color12=clear_color;
            color21=clear_color;
            color22=clear_color;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol  ])color11=burn_color;
            if(meshi->thresholdtime[nn1+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol+1])color12=burn_color;
            if(meshi->thresholdtime[nn2+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol  ])color21=burn_color;
            if(meshi->thresholdtime[nn2+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol+1])color22=burn_color;
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
             if(ABS(r11-r22)<ABS(r12-r21)){
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
          cpatchval1++; cpatchval2++; patchblank1++; patchblank2++;
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
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(meshi->visPatches[n]==1&&meshi->patchdir[n]>0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      cpatchval_iframe_copy = cpatchval_iframe + blockstart[n];
      for(irow=0;irow<nrow-1;irow++){
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        cpatchval1 = cpatchval_iframe_copy + irow*ncol;
        patchblank1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;

        xyzp2 = xyzp1 + 3*ncol;
        cpatchval2 = cpatchval1 + ncol;
        patchblank2 = patchblank1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            r11 = (float)(*cpatchval1)/255.0;
            r12 = (float)(*(cpatchval1+1))/255.0;
            r21 = (float)(*cpatchval2)/255.0;
            r22 = (float)(*(cpatchval2+1))/255.0;
            color11=clear_color;
            color12=clear_color;
            color21=clear_color;
            color22=clear_color;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol  ])color11=burn_color;
            if(meshi->thresholdtime[nn1+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol+1])color12=burn_color;
            if(meshi->thresholdtime[nn2+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol  ])color21=burn_color;
            if(meshi->thresholdtime[nn2+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol+1])color22=burn_color;
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
              if(ABS(*cpatchval1-*(cpatchval2+1))<ABS(*(cpatchval1+1)-*cpatchval2)){
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
          cpatchval1++; cpatchval2++; patchblank1++; patchblank2++;
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
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    if(visPatches[n]==1&&patchdir[n]<0){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      cpatchval_iframe_copy = cpatchval_iframe + blockstart[n];
      for(irow=0;irow<nrow-1;irow++){
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        patchblank1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        cpatchval1 = cpatchval_iframe_copy + irow*ncol;
        cpatchval2 = cpatchval1 + ncol;
        patchblank2 = patchblank1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            r11 = (float)(*cpatchval1)/255.0;
            r12 = (float)(*(cpatchval1+1))/255.0;
            r21 = (float)(*cpatchval2)/255.0;
            r22 = (float)(*(cpatchval2+1))/255.0;
            color11=clear_color;
            color12=clear_color;
            color21=clear_color;
            color22=clear_color;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol  ])color11=burn_color;
            if(meshi->thresholdtime[nn1+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol+1])color12=burn_color;
            if(meshi->thresholdtime[nn2+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol  ])color21=burn_color;
            if(meshi->thresholdtime[nn2+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol+1])color22=burn_color;
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
              if(ABS(*cpatchval1-*(cpatchval2+1))<ABS(*(cpatchval1+1)-*cpatchval2)){
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
          cpatchval1++; cpatchval2++; patchblank1++; patchblank2++;
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

void drawpatch_threshold_cellcenter(const meshdata *meshi){
  int n,nn,nn1;
  int nrow, ncol, irow, icol;
  float *xyzpatchcopy;
  int *patchblankcopy;
  float *patch_times;
  int *visPatches;
  float *xyzpatch;
  int *patchdir, *patchrow, *patchcol;
  int *blockstart;
  int *patchblank;
  int iblock;
  blockagedata *bc;
  patchdata *patchi;
  float *color11;
  meshdata *meshblock;
  float burn_color[4]={0.0,0.0,0.0,1.0};
  float clear_color[4]={1.0,1.0,1.0,1.0};

  if(vis_threshold==1&&vis_onlythreshold==1&&do_threshold==1)return;

  patch_times=meshi->patch_times;
  visPatches=meshi->visPatches;
  xyzpatch=meshi->xyzpatch;
  patchdir=meshi->patchdir;
  patchrow=meshi->patchrow;
  patchcol=meshi->patchcol;
  blockstart=meshi->blockstart;
  patchblank=meshi->patchblank;
  patchi=patchinfo+meshi->patchfilenum;
  switch(patchi->compression_type){
  case UNCOMPRESSED:
    ASSERT(meshi->cpatchval_iframe!=NULL);
    break;
  case COMPRESSED_ZLIB:
    ASSERT(meshi->cpatchval_iframe_zlib!=NULL);
    break;
  default:
    ASSERT(FFALSE);
  }
  patchi = patchinfo + meshi->patchfilenum;

  if(patch_times[0]>global_times[itimes]||patchi->display==0)return;
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
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
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
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        patchblank1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        patchblank2 = patchblank1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            color11=clear_color;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol  ])color11=burn_color;

            glColor4fv(color11);
            glVertex3fv(xyzp1);
            glVertex3fv(xyzp1+3);
            glVertex3fv(xyzp2+3);

            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2+3);
            glVertex3fv(xyzp2);
          }
          patchblank1++; patchblank2++;
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
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
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
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        patchblank1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;

        xyzp2 = xyzp1 + 3*ncol;
        patchblank2 = patchblank1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            color11=clear_color;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol  ])color11=burn_color;

            glColor4fv(color11);
            glVertex3fv(xyzp1);
            glVertex3fv(xyzp1+3);
            glVertex3fv(xyzp2+3);

            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2+3);
            glVertex3fv(xyzp2);
          }
          patchblank1++; patchblank2++;
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
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
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
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        patchblank1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        patchblank2 = patchblank1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            color11=clear_color;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol  ])color11=burn_color;

            glColor4fv(color11);
            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2+3);
            glVertex3fv(xyzp1+3);
            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2);
            glVertex3fv(xyzp2+3);
          }
          patchblank1++; patchblank2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
    }
    nn += patchrow[n]*patchcol[n];
  }
  glEnd();
}

/* ------------------ draw_patchframe ------------------------ */

void draw_patchframe(int flag){
  meshdata *meshi;
  int i;

  for(i=0;i<npatchinfo;i++){
    patchdata *patchi;

    patchi = patchinfo + i;
    if(patchi->filetype == PATCH_GEOMETRY && patchi->loaded == 1 && patchi->display == 1){
      if(flag == DRAW_OPAQUE){
        if(patchi->slice == 0){
          draw_geomdata(flag, patchi, GEOM_STATIC);
          draw_geomdata(flag, patchi, GEOM_DYNAMIC);
        }
      }
      else{
        if(patchi->slice == 1){
          draw_geomdata(flag, patchi, GEOM_STATIC);
          draw_geomdata(flag, patchi, GEOM_DYNAMIC);
        }
      }
    }
  }
  if(flag == DRAW_TRANSPARENT)return;
  for(i=0;i<nmeshes;i++){
    meshi=meshinfo+i;
    if(meshi->npatches>0){
      int filenum;

      filenum=meshi->patchfilenum;
      if(filenum!=-1){
        patchdata *patchi;

        patchi = patchinfo + filenum;
        if(patchi->loaded==0||patchi->display==0||patchi->type!=ipatchtype)continue;
        if(usetexturebar!=0){
          if(vis_threshold==1&&do_threshold==1){
            if(patchi->filetype==PATCH_CELL_CENTER){
              drawpatch_threshold_cellcenter(meshi);
            }
            else if(patchi->filetype==PATCH_NODE_CENTER){
              drawpatch_texture_threshold(meshi);
            }
          }
          else{
            if(patchi->filetype==PATCH_CELL_CENTER){
              drawpatch_cellcenter(meshi);
            }
            else if(patchi->filetype==PATCH_NODE_CENTER){
              drawpatch_texture(meshi);
            }
          }
        }
        else{
          if(patchi->filetype==PATCH_CELL_CENTER){
            drawpatch_cellcenter(meshi);
          }
          else if(patchi->filetype==PATCH_NODE_CENTER){
            drawpatch(meshi);
          }
        }
        if(vis_threshold==1&&vis_onlythreshold==1&&do_threshold==1)drawonlythreshold(meshi);
      }
    }
  }
}

/* ------------------ drawpatch ------------------------ */

void drawpatch(const meshdata *meshi){
  int n,nn,nn1,nn2;
  int nrow, ncol, irow, icol;
  unsigned char *cpatchval1, *cpatchval2;
  unsigned char *cpatchval_iframe_copy;
  float *xyzpatchcopy;
  int *patchblankcopy;
  float *patch_times;
  int *visPatches;
  float *xyzpatch;
  int *patchdir, *patchrow, *patchcol, *patchtype;
  int *blockstart;
  int *patchblank;
  unsigned char *cpatchval_iframe;
  int iblock;
  blockagedata *bc;
  patchdata *patchi;
  float *color11, *color12, *color21, *color22;
  meshdata *meshblock;
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

  patch_times=meshi->patch_times;
  visPatches=meshi->visPatches;
  xyzpatch=meshi->xyzpatch;
  patchdir=meshi->patchdir;
  patchtype=meshi->patchtype;
  patchrow=meshi->patchrow;
  patchcol=meshi->patchcol;
  blockstart=meshi->blockstart;
  patchblank=meshi->patchblank;
  patchi=patchinfo+meshi->patchfilenum;
  switch(patchi->compression_type){
  case UNCOMPRESSED:
    ASSERT(meshi->cpatchval_iframe!=NULL);
    cpatchval_iframe=meshi->cpatchval_iframe;
    break;
  case COMPRESSED_ZLIB:
    ASSERT(meshi->cpatchval_iframe_zlib!=NULL);
    cpatchval_iframe=meshi->cpatchval_iframe_zlib;
    break;
  default:
    ASSERT(FFALSE);
  }
  patchi = patchinfo + meshi->patchfilenum;

  if(patch_times[0]>global_times[itimes]||patchi->display==0)return;
  if(cullfaces==1)glDisable(GL_CULL_FACE);

  /* if a contour boundary does not match a blockage face then draw "both sides" of boundary */

  if((use_transparency_data==1&&contour_type==LINE_CONTOURS)||setpatchchopmin==1||setpatchchopmax==1)transparenton();
  nn =0;
  glBegin(GL_TRIANGLES);
  for(n=0;n<meshi->npatches;n++){
    int drawit;

    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    drawit=0;
    if(visPatches[n]==1&&patchdir[n]==0)drawit=1;
    if(patchtype[n]!=INTERIORwall&&showpatch_both==1)drawit=1;
    if(drawit==1){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      cpatchval_iframe_copy = cpatchval_iframe + blockstart[n];
      for(irow=0;irow<nrow-1;irow++){
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        patchblank1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        cpatchval1 = cpatchval_iframe_copy + irow*ncol;
        cpatchval2 = cpatchval1 + ncol;
        patchblank2 = patchblank1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            color11 = rgb_patch+4*(*cpatchval1);
            color12 = rgb_patch+4*(*(cpatchval1+1));
            color21 = rgb_patch+4*(*cpatchval2);
            color22 = rgb_patch+4*(*(cpatchval2+1));
            if(vis_threshold==1&&vis_onlythreshold==0&&do_threshold==1){
              if(meshi->thresholdtime[nn1+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol  ])color11=&char_color[0];
              if(meshi->thresholdtime[nn1+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol+1])color12=&char_color[0];
              if(meshi->thresholdtime[nn2+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol  ])color21=&char_color[0];
              if(meshi->thresholdtime[nn2+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol+1])color22=&char_color[0];
            }
            if(ABS(*cpatchval1-*(cpatchval2+1))<ABS(*(cpatchval1+1)-*cpatchval2)){
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
          cpatchval1++; cpatchval2++; patchblank1++; patchblank2++;
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
    int drawit;

    iblock = meshi->blockonpatch[n];
    meshblock=meshi->meshonpatch[n];
    if(iblock!=-1){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    drawit=0;
    if(meshi->visPatches[n]==1&&meshi->patchdir[n]>0){
      if(patchtype[n]==INTERIORwall||showpatch_both==0){
        drawit=1;
      }
    }
    if(drawit==1){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      cpatchval_iframe_copy = cpatchval_iframe + blockstart[n];
      if(hidepatchsurface==0){
        glPushMatrix();
        switch(meshi->patchdir[n]){
          case XDIR:
            glTranslatef(dboundx,0.0,0.0);
            break;
          case YDIR:
            glTranslatef(0.0,-dboundy,0.0);
            break;
          case ZDIR:
            glTranslatef(0.00,0.0,dboundz);
            break;
          default:
            ASSERT(FFALSE);
            break;
        }
        glBegin(GL_TRIANGLES);
      }
      for(irow=0;irow<nrow-1;irow++){
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        cpatchval1 = cpatchval_iframe_copy + irow*ncol;
        patchblank1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;

        xyzp2 = xyzp1 + 3*ncol;
        cpatchval2 = cpatchval1 + ncol;
        patchblank2 = patchblank1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            color11 = rgb_patch+4*(*cpatchval1);
            color12 = rgb_patch+4*(*(cpatchval1+1));
            color21 = rgb_patch+4*(*cpatchval2);
            color22 = rgb_patch+4*(*(cpatchval2+1));
            if(vis_threshold==1&&vis_onlythreshold==0&&do_threshold==1){
              if(meshi->thresholdtime[nn1+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol  ])color11=&char_color[0];
              if(meshi->thresholdtime[nn1+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol+1])color12=&char_color[0];
              if(meshi->thresholdtime[nn2+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol  ])color21=&char_color[0];
              if(meshi->thresholdtime[nn2+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol+1])color22=&char_color[0];
            }
            if(ABS(*cpatchval1-*(cpatchval2+1))<ABS(*(cpatchval1+1)-*cpatchval2)){
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
          cpatchval1++; cpatchval2++; patchblank1++; patchblank2++;
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
    int drawit;

    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    drawit=0;
    if(visPatches[n]==1&&patchdir[n]<0){
      if(patchtype[n]==INTERIORwall||showpatch_both==0){
        drawit=1;
      }
    }
    if(drawit==1){
      nrow=patchrow[n];
      ncol=patchcol[n];
      xyzpatchcopy = xyzpatch + 3*blockstart[n];
      patchblankcopy = patchblank + blockstart[n];
      cpatchval_iframe_copy = cpatchval_iframe + blockstart[n];
      if(hidepatchsurface==0){
        glPushMatrix();
        switch(meshi->patchdir[n]){
          case XDIRNEG:
            glTranslatef(-dboundx,0.0,0.0);
            break;
          case YDIRNEG:
            glTranslatef(0.0,dboundy,0.0);
            break;
          case ZDIRNEG:
            glTranslatef(0.00,0.0,-dboundz);
            break;
          default:
            ASSERT(FFALSE);
            break;
        }
        glBegin(GL_TRIANGLES);
      }
      for(irow=0;irow<nrow-1;irow++){
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        patchblank1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        cpatchval1 = cpatchval_iframe_copy + irow*ncol;
        cpatchval2 = cpatchval1 + ncol;
        patchblank2 = patchblank1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            color11 = rgb_patch+4*(*cpatchval1);
            color12 = rgb_patch+4*(*(cpatchval1+1));
            color21 = rgb_patch+4*(*cpatchval2);
            color22 = rgb_patch+4*(*(cpatchval2+1));
            if(vis_threshold==1&&vis_onlythreshold==0&&do_threshold==1){
              if(meshi->thresholdtime[nn1+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol  ])color11=&char_color[0];
              if(meshi->thresholdtime[nn1+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol+1])color12=&char_color[0];
              if(meshi->thresholdtime[nn2+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol  ])color21=&char_color[0];
              if(meshi->thresholdtime[nn2+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol+1])color22=&char_color[0];
            }
            if(ABS(*cpatchval1-*(cpatchval2+1))<ABS(*(cpatchval1+1)-*cpatchval2)){
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
          cpatchval1++; cpatchval2++; patchblank1++; patchblank2++;
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
  if((use_transparency_data==1&&contour_type==LINE_CONTOURS)||setpatchchopmin==1||setpatchchopmax==1)transparentoff();
}

/* ------------------ drawpatch_cellcenter ------------------------ */

void drawpatch_cellcenter(const meshdata *meshi){
  int n,nn,nn1;
  int nrow, ncol, irow, icol;
  unsigned char *cpatchval1;
  unsigned char *cpatchval_iframe_copy;
  float *patch_times;
  int *visPatches;
  int *patchdir, *patchrow, *patchcol, *patchtype;
  int *blockstart;
  unsigned char *cpatchval_iframe;
  int iblock;
  blockagedata *bc;
  patchdata *patchi;
  float *color11;
  meshdata *meshblock;

  float dboundx,dboundy,dboundz;
  float *xplt, *yplt, *zplt;
  float **patchventcolors;

  if(vis_threshold==1&&vis_onlythreshold==1&&do_threshold==1)return;

  if(hidepatchsurface==0){
    xplt=meshi->xplt;
    yplt=meshi->yplt;
    zplt=meshi->zplt;

    dboundx = (xplt[1]-xplt[0])/10.0;
    dboundy = (yplt[1]-yplt[0])/10.0;
    dboundz = (zplt[1]-zplt[0])/10.0;
  }

  patch_times=meshi->patch_times;
  visPatches=meshi->visPatches;
  patchdir=meshi->patchdir;
  patchtype=meshi->patchtype;
  patchrow=meshi->patchrow;
  patchcol=meshi->patchcol;
  blockstart=meshi->blockstart;
  patchventcolors=meshi->patchventcolors;
  patchi=patchinfo+meshi->patchfilenum;

  switch(patchi->compression_type){
  case UNCOMPRESSED:
    ASSERT(meshi->cpatchval_iframe!=NULL);
    cpatchval_iframe=meshi->cpatchval_iframe;
    break;
  case COMPRESSED_ZLIB:
    ASSERT(meshi->cpatchval_iframe_zlib!=NULL);
    cpatchval_iframe=meshi->cpatchval_iframe_zlib;
    break;
  default:
    ASSERT(FFALSE);
  }
  patchi = patchinfo + meshi->patchfilenum;

  if(patch_times[0]>global_times[itimes]||patchi->display==0)return;
  if(cullfaces==1)glDisable(GL_CULL_FACE);

  /* if a contour boundary does not match a blockage face then draw "both sides" of boundary */

  if((use_transparency_data==1&&contour_type==LINE_CONTOURS)||setpatchchopmin==1||setpatchchopmax==1)transparenton();
  nn =0;
  glBegin(GL_TRIANGLES);
  for(n=0;n<meshi->npatches;n++){
    int drawit;

    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    drawit=0;
    if(visPatches[n]==1&&patchdir[n]==0)drawit=1;
    if(patchtype[n]!=INTERIORwall&&showpatch_both==1)drawit=1;
    if(drawit==1){
      nrow=patchrow[n];
      ncol=patchcol[n];
      cpatchval_iframe_copy = cpatchval_iframe + blockstart[n];
      for(irow=0;irow<nrow-1;irow++){
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = meshi->xyzpatch + 3*blockstart[n] + 3*irow*ncol;
        patchblank1 = meshi->patchblank + blockstart[n] + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        cpatchval1 = cpatchval_iframe_copy + irow*ncol;
        patchblank2 = patchblank1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            if(patchventcolors==NULL){
              color11 = rgb_patch+4*(*cpatchval1);
              if(vis_threshold==1&&vis_onlythreshold==0&&do_threshold==1){
                if(meshi->thresholdtime[nn1+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol  ])color11=&char_color[0];
              }
            }
            else{
              color11=patchventcolors[(cpatchval1-cpatchval_iframe)];
            }
            glColor4fv(color11);
            glVertex3fv(xyzp1);
            glVertex3fv(xyzp1+3);
            glVertex3fv(xyzp2+3);

            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2+3);
            glVertex3fv(xyzp2);
          }
          cpatchval1++;
          patchblank1++;
          patchblank2++;
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
    int drawit;

    iblock = meshi->blockonpatch[n];
    meshblock=meshi->meshonpatch[n];
    if(iblock!=-1){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    drawit=0;
    if(meshi->visPatches[n]==1&&meshi->patchdir[n]>0){
      if(patchtype[n]==INTERIORwall||showpatch_both==0){
        drawit=1;
      }
    }
    if(drawit==1){
      nrow=patchrow[n];
      ncol=patchcol[n];
      cpatchval_iframe_copy = cpatchval_iframe + blockstart[n];
      if(hidepatchsurface==0){
        glPushMatrix();
        switch(meshi->patchdir[n]){
          case XDIR:
            glTranslatef(dboundx,0.0,0.0);
            break;
          case YDIR:
            glTranslatef(0.0,-dboundy,0.0);
            break;
          case ZDIR:
            glTranslatef(0.0,0.0,dboundz);
            break;
          default:
            ASSERT(FFALSE);
            break;
        }
        glBegin(GL_TRIANGLES);
      }
      for(irow=0;irow<nrow-1;irow++){
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = meshi->xyzpatch + 3*blockstart[n] + 3*irow*ncol;
        cpatchval1 = cpatchval_iframe_copy + irow*ncol;
        patchblank1 = meshi->patchblank + blockstart[n] + irow*ncol;
        nn1 = nn + irow*ncol;

        xyzp2 = xyzp1 + 3*ncol;
        patchblank2 = patchblank1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            if(patchventcolors==NULL){
              color11 = rgb_patch+4*(*cpatchval1);
              if(vis_threshold==1&&vis_onlythreshold==0&&do_threshold==1){
                if(meshi->thresholdtime[nn1+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol  ])color11=&char_color[0];
              }
            }
            else{
              color11=patchventcolors[(cpatchval1-cpatchval_iframe)];
            }
            glColor4fv(color11);
            glVertex3fv(xyzp1);
            glVertex3fv(xyzp1+3);
            glVertex3fv(xyzp2+3);

            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2+3);
            glVertex3fv(xyzp2);
          }
          cpatchval1++;
          patchblank1++;
          patchblank2++;
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
    int drawit;

    iblock = meshi->blockonpatch[n];
    meshblock = meshi->meshonpatch[n];
    ASSERT((iblock!=-1&&meshblock!=NULL)||(iblock==-1&&meshblock==NULL));
    if(iblock!=-1&&meshblock!=NULL){
      bc=meshblock->blockageinfoptrs[iblock];
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
        nn += patchrow[n]*patchcol[n];
        continue;
      }
    }
    drawit=0;
    if(visPatches[n]==1&&patchdir[n]<0){
      if(patchtype[n]==INTERIORwall||showpatch_both==0){
        drawit=1;
      }
    }
    if(drawit==1){
      nrow=patchrow[n];
      ncol=patchcol[n];
      cpatchval_iframe_copy = cpatchval_iframe + blockstart[n];
      if(hidepatchsurface==0){
        glPushMatrix();
        switch(meshi->patchdir[n]){
          case XDIRNEG:
            glTranslatef(-dboundx,0.0,0.0);
            break;
          case YDIRNEG:
            glTranslatef(0.0,dboundy,0.0);
            break;
          case ZDIRNEG:
            glTranslatef(0.0,0.0,-dboundz);
            break;
          default:
            ASSERT(FFALSE);
            break;
        }
        glBegin(GL_TRIANGLES);
      }
      for(irow=0;irow<nrow-1;irow++){
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = meshi->xyzpatch + 3*blockstart[n] + 3*irow*ncol;
        patchblank1 = meshi->patchblank + blockstart[n] + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        cpatchval1 = cpatchval_iframe_copy + irow*ncol;
        patchblank2 = patchblank1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            if(patchventcolors==NULL){
              color11 = rgb_patch+4*(*cpatchval1);
              if(vis_threshold==1&&vis_onlythreshold==0&&do_threshold==1){
                if(meshi->thresholdtime[nn1+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol  ])color11=&char_color[0];
              }
            }
            else{
              color11=patchventcolors[(cpatchval1-cpatchval_iframe)];
            }
            glColor4fv(color11);
            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2+3);
            glVertex3fv(xyzp1+3);

            glVertex3fv(xyzp1);
            glVertex3fv(xyzp2);
            glVertex3fv(xyzp2+3);
          }
          cpatchval1++;
          patchblank1++;
          patchblank2++;
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
  if((use_transparency_data==1&&contour_type==LINE_CONTOURS)||setpatchchopmin==1||setpatchchopmax==1)transparentoff();
}

/* ------------------ drawolythreshold ------------------------ */

void drawonlythreshold(const meshdata *meshi){
  int n,nn,nn1,nn2;
  int nrow, ncol, irow, icol;
  float *xyzpatchcopy;
  int *patchblankcopy;
  float *patch_times;
  int *visPatches;
  float *xyzpatch;
  int *patchdir, *patchrow, *patchcol;
  int *blockstart;
  int *patchblank;
  int iblock;
  blockagedata *bc;
  patchdata *patchi;
  float *color11, *color12, *color21, *color22;
  float *color_black;
  meshdata *meshblock;

  if(vis_threshold==0||vis_onlythreshold==0||do_threshold==0)return;

  patch_times=meshi->patch_times;
  visPatches=meshi->visPatches;
  xyzpatch=meshi->xyzpatch_threshold;
  patchdir=meshi->patchdir;
  patchrow=meshi->patchrow;
  patchcol=meshi->patchcol;
  blockstart=meshi->blockstart;
  patchblank=meshi->patchblank;
  patchi=patchinfo+meshi->patchfilenum;
  switch(patchi->compression_type){
  case UNCOMPRESSED:
    ASSERT(meshi->cpatchval_iframe!=NULL);
    break;
  case COMPRESSED_ZLIB:
    ASSERT(meshi->cpatchval_iframe_zlib!=NULL);
    break;
  default:
    ASSERT(FFALSE);
  }
  patchi = patchinfo + meshi->patchfilenum;

  if(patch_times[0]>global_times[itimes]||patchi->display==0)return;
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
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
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
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        patchblank1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        patchblank2 = patchblank1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            int nnulls;

            nnulls=4;
            color11 = NULL;
            color12 = NULL;
            color21 = NULL;
            color22 = NULL;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol  ]){
              color11=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn1+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol+1]){
              color12=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn2+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol  ]){
              color21=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn2+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol+1]){
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
          patchblank1++; patchblank2++;
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
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
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
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        patchblank1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;

        xyzp2 = xyzp1 + 3*ncol;
        patchblank2 = patchblank1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            int nnulls;

            color11 = NULL;
            color12 = NULL;
            color21 = NULL;
            color22 = NULL;
            nnulls=4;

            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol  ]){
              color11=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn1+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol+1]){
              color12=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn2+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol  ]){
              color21=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn2+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol+1]){
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
          patchblank1++; patchblank2++;
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
      if(bc->showtimelist!=NULL&&bc->showtimelist[itimes]==0){
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
        int *patchblank1, *patchblank2;
        float *xyzp1, *xyzp2;

        xyzp1 = xyzpatchcopy + 3*irow*ncol;
        patchblank1 = patchblankcopy + irow*ncol;
        nn1 = nn + irow*ncol;
        xyzp2 = xyzp1 + 3*ncol;
        patchblank2 = patchblank1 + ncol;
        nn2 = nn1 + ncol;

        for(icol=0;icol<ncol-1;icol++){
          if(*patchblank1==GAS&&*patchblank2==GAS&&*(patchblank1+1)==GAS&&*(patchblank2+1)==GAS){
            int nnulls;

            color11 = NULL;
            color12 = NULL;
            color21 = NULL;
            color22 = NULL;
            nnulls=4;
            if(meshi->thresholdtime[nn1+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol  ]){
              color11=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn1+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn1+icol+1]){
              color12=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn2+icol  ]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol  ]){
              color21=&char_color[0];
              nnulls--;
            }
            if(meshi->thresholdtime[nn2+icol+1]>=0.0&&global_times[itimes]>meshi->thresholdtime[nn2+icol+1]){
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
          patchblank1++; patchblank2++;
          xyzp1+=3;
          xyzp2+=3;
        }
      }
    }
    nn += patchrow[n]*patchcol[n];
  }
  glEnd();
}

/* ------------------ updatepatchtypes ------------------------ */

void updatepatchtypes(void){
  int i;
  patchdata *patchi;

  npatchtypes = 0;
  for(i=0;i<npatchinfo;i++){
    patchi = patchinfo+i;
    if(getpatchindex(patchi)==-1)patchtypes[npatchtypes++]=i;
  }
  for(i=0;i<npatchinfo;i++){
    patchi = patchinfo+i;
    patchi->type=getpatchtype(patchi);
  }
}

/* ------------------ getpatchtype ------------------------ */

int getpatchtype(const patchdata *patchi){
  int j;

  for(j=0;j<npatchtypes;j++){
    patchdata *patchi2;

    patchi2 = patchinfo+patchtypes[j];
    if(strcmp(patchi->label.shortlabel,patchi2->label.shortlabel)==0)return j;
  }
  return -1;
}

/* ------------------ update_patchtype ------------------------ */

void update_patchtype(){
  int i;

  for(i=0;i<npatchinfo;i++){
    patchdata *patchi;

    patchi = patchinfo + i;
    if(patchi->loaded==1&&patchi->display==1&&patchi->type==ipatchtype)return;
  }

  for(i=0;i<npatchinfo;i++){
    patchdata *patchi;

    patchi = patchinfo + i;
    if(patchi->loaded==1&&patchi->display==1){
      ipatchtype = getpatchindex(patchi);
      return;
    }
  }
  ipatchtype = -1;
  return;

}

/* ------------------ patchcompare ------------------------ */

int patchcompare( const void *arg1, const void *arg2 ){
  patchdata *patchi, *patchj;

  patchi = patchinfo + *(int *)arg1;
  patchj = patchinfo + *(int *)arg2;

  if(strcmp(patchi->label.longlabel,patchj->label.longlabel)<0)return -1;
  if(strcmp(patchi->label.longlabel,patchj->label.longlabel)>0)return 1;
  if(patchi->blocknumber<patchj->blocknumber)return -1;
  if(patchi->blocknumber>patchj->blocknumber)return 1;
  return 0;
}

/* ------------------ update_patch_menulabels ------------------------ */

void update_patch_menulabels(void){
  int i;
  patchdata *patchi;
  char label[128];
  STRUCTSTAT statbuffer;

  if(npatchinfo>0){
    FREEMEMORY(patchorderindex);
    NewMemory((void **)&patchorderindex,sizeof(int)*npatchinfo);
    for(i=0;i<npatchinfo;i++){
      patchorderindex[i]=i;
    }
    qsort( (int *)patchorderindex, (size_t)npatchinfo, sizeof(int), patchcompare );

    for(i=0;i<npatchinfo;i++){
      patchi = patchinfo + i;
      STRCPY(patchi->menulabel,patchi->label.longlabel);
      if(nmeshes>1){
        meshdata *patchmesh;

        patchmesh = meshinfo + patchi->blocknumber;
        sprintf(label,"%s",patchmesh->label);
        STRCAT(patchi->menulabel,", ");
        STRCAT(patchi->menulabel,label);
      }
      if(patchi->filetype == PATCH_GEOMETRY){
        if(strlen(patchi->gslicedir) != 0){
          STRCAT(patchi->menulabel, ", ");
          STRCAT(patchi->menulabel, patchi->gslicedir);
        }
      }
      if(STAT(patchi->comp_file,&statbuffer)==0){
        patchi->file=patchi->comp_file;
        patchi->compression_type=COMPRESSED_ZLIB;
      }
      else{
        patchi->file=patchi->reg_file;
        patchi->compression_type=UNCOMPRESSED;
      }
      if(showfiles==1){
        STRCAT(patchi->menulabel,", ");
        STRCAT(patchi->menulabel,patchi->file);
      }
      if(patchi->compression_type==COMPRESSED_ZLIB){
        STRCAT(patchi->menulabel," (ZLIB)");
      }
    }
  }


}

/* ------------------ getpatchheader ------------------------ */

void getpatchheader(char *file,int *npatches, float *ppatchmin, float *ppatchmax){
  FILE *stream;
  float minmax[2];

  stream=fopen(file,"rb");

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

  FSEEK(stream,16,SEEK_CUR);

  fread(minmax,4,2,stream);
  FSEEK(stream,8,SEEK_CUR);
  fread(npatches,4,1,stream);
  *ppatchmin=minmax[0];
  *ppatchmax=minmax[1];
  fclose(stream);
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
  FILE *stream;
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

  stream=fopen(file,"rb");

  if(stream==NULL)return;

  FSEEK(stream,12,SEEK_CUR);
  fread(version,4,1,stream);
  FSEEK(stream,16,SEEK_CUR);
  fread(&npatches,4,1,stream);
  for(i=0;i<npatches;i++){
    buffer[6]=0;
    if(*version==0){
      fread(buffer,4,6,stream);
    }
    else{
      fread(buffer,4,9,stream);
    }
    pi1[i]=buffer[0];
    pi2[i]=buffer[1];
    pj1[i]=buffer[2];
    pj2[i]=buffer[3];
    pk1[i]=buffer[4];
    pk2[i]=buffer[5];
    patchdir[i]=buffer[6];
  }
  fclose(stream);
}

/* ------------------ getpatchsizeinfo ------------------------ */

void getpatchsizeinfo(patchdata *patchi, int *nframes, int *buffersize){
  FILE *streamsize;
  FILE *stream;
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

    stream=fopen(patchi->file,"rb");
    if(stream==NULL)return;

    streamsize=fopen(sizefile,"w");
    if(streamsize==NULL){
      fclose(stream);
      return;
    }

    FSEEK(stream,12,SEEK_CUR);
    fread(&version,4,1,stream);
    FSEEK(stream,16,SEEK_CUR);
    fread(&npatches,4,1,stream);
    size=0;
    return_code=0;
    for(i=0;i<npatches;i++){
      if(version==0){
        return_code=fread(buff,4,6,stream);
      }
      else{
        return_code=fread(buff,4,9,stream);
      }
      if(return_code==0)break;
      size+=(buff[1]+1-buff[0])*(buff[3]+1-buff[2])*(buff[5]+1-buff[4]);
    }
    if(return_code==0){
      fclose(stream);
      fclose(streamsize);
      return;
    }
    for(;;){
      return_code=fread(&local_time,4,1,stream);
      if(return_code==0)break;
      return_code=fread(&compressed_size,4,1,stream);
      if(return_code==0)break;
      return_code=FSEEK(stream,compressed_size,SEEK_CUR);
      if(return_code!=0)break;
      fprintf(streamsize,"%f %i %i\n",local_time,size,compressed_size);
    }
    fclose(stream);
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

void getpatchdata_zlib(patchdata *patchi,unsigned char *data,int ndata,
                       float *local_times, unsigned int *zipoffset, unsigned int *zipsize, int ntimes_local){
  FILE *stream;
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

  stream=fopen(patchi->file,"rb");
  if(stream==NULL)return;

  FSEEK(stream,12,SEEK_CUR);
  fread(&version,4,1,stream);
  FSEEK(stream,16,SEEK_CUR);
  fread(&npatches,4,1,stream);
  if(version==0){
    local_skip = 6*npatches*4;
  }
  else{
    local_skip = 9*npatches*4;
  }
  return_code=FSEEK(stream,local_skip,SEEK_CUR);
  if(return_code!=0){
    fclose(stream);
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

    if(fread(&local_time,4,1,stream)==0)break;
    skip_frame=1;
    if(local_time>time_max){
      time_max=local_time;
      skip_frame=0;
      local_count++;
    }
    if(fread(&compressed_size,4,1,stream)==0)break;
    if(skip_frame==0&&local_count%boundframestep==0){
      if(fread(datacopy,1,compressed_size,stream)==0)break;
    }
    else{
      FSEEK(stream,compressed_size,SEEK_CUR);
    }

    if(skip_frame==1||local_count%boundframestep!=0)continue;
    i++;
    if(i>=ntimes_local)break;
    PRINTF("boundary time=%.2f\n",local_time);
    ASSERT(i<ntimes_local);
    local_times[i]=local_time;
    zipoffset[i]=offset;
    zipsize[i]=compressed_size;
    datacopy+=compressed_size;
    offset+=compressed_size;
  }
  fclose(stream);
}

/* ------------------ uncompress_patchdataframe ------------------------ */

void uncompress_patchdataframe(meshdata *meshi,int local_iframe){
  unsigned int countin;
  uLongf countout;
  unsigned char *compressed_data;

  compressed_data = meshi->cpatchval_zlib+meshi->zipoffset[local_iframe];
  countin = meshi->zipsize[local_iframe];
  countout=meshi->npatchsize;

  uncompress_zlib(meshi->cpatchval_iframe_zlib,&countout,compressed_data,countin);

}

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


/* ------------------ Update_All_Patch_Bounds ------------------------ */

void Update_All_Patch_Bounds_st(void){
  int i;
  int total=0;

  LOCK_COMPRESS;
  for(i=0;i<npatchinfo;i++){
    patchdata *patchi;

    patchi = patchinfo + i;
    total+=update_patch_hist(patchi);
    update_patch_bounds(patchi);
  }
  if(total==0){
    PRINTF("Boundary file bounds already computed.\n");
  }
  else{
    PRINTF("Bounds for %i boundary files computed\n",total);
  }
  UNLOCK_COMPRESS;
}

/* ------------------ update_patch_hist ------------------------ */

int update_patch_hist(patchdata *patchj){
  int i;
  int first=1;
  int sum=0;

  if(patchj->setvalmax==SET_MAX&&patchj->setvalmin==SET_MIN)return 0;

  for(i=0;i<npatchinfo;i++){
    int npatches, error;
    patchdata *patchi;
    int unit1;
    FILE_SIZE lenfile;
    int error1;
    int *pi1, *pi2, *pj1, *pj2, *pk1, *pk2, *patchdir, *patchsize;
    float patchtime1, *patchframe;
    int patchframesize;
    int j;
    time_t modtime;

    patchi = patchinfo + i;
    if(patchi->type!=patchj->type||patchi->filetype!=patchj->filetype||patchi->filetype==PATCH_GEOMETRY)continue;
    modtime=file_modtime(patchi->file);
    if(modtime>patchi->modtime){
      patchi->modtime=modtime;
      patchi->inuse_getbounds=0;
      patchi->histogram->complete=0;
      patchi->bounds.defined=0;
    }
    if(patchi->inuse_getbounds==1||patchi->histogram->complete==1||patchi->bounds.defined==1)continue;

    patchi->inuse_getbounds=1;

    if(first==1){
      PRINTF("Determining %s percentile and global data bounds\n",patchi->label.longlabel);
      first=0;
    }
    PRINTF("  Examining %s\n",patchi->file);
    sum++;
    lenfile=strlen(patchi->file);

    FORTget_file_unit(&unit1,&patchi->unit_start);
    FORTgetboundaryheader1(patchi->file,&unit1, &npatches, &error, lenfile);
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
      int npatchsize;

      npatchsize = (pi2[j]+1-pi1[j]);
      npatchsize *= (pj2[j]+1-pj1[j]);
      npatchsize *= (pk2[j]+1-pk1[j]);
      patchframesize+=npatchsize;
    }

    FORTget_file_unit(&unit1,&patchi->unit_start);
    FORTopenboundary(patchi->file,&unit1,&patchi->version,&error1,lenfile);

    NewMemory((void **)&patchframe,patchframesize*sizeof(float));
    reset_histogram(patchi->histogram);
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


