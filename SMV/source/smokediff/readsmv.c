#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "svdiff.h"
#include "MALLOC.h"

/* ------------------ readsmv ------------------------ */

int readsmv(FILE *streamsmv, FILE *stream_out, casedata *smvcase){

  int igrid,ipdim;
  int islice,iplot3d,iboundary;
  char buffer[255];
  meshdata *meshinfo=NULL;
  slice *sliceinfo=NULL;
  boundary *boundaryinfo=NULL;
  plot3d *plot3dinfo=NULL;
  int nmeshes, nsliceinfo, nplot3dinfo, nboundary_files;
  int itrnx, itrny, itrnz;

  igrid=0;
  ipdim=0;
  nmeshes=0;
  nsliceinfo=0;
  nplot3dinfo=0;
  nboundary_files=0;
  itrnx=0;
  itrny=0;
  itrnz=0;

  while(!feof(streamsmv)){
    if(fgets(buffer,255,streamsmv)==NULL)break;
    CheckMemory;
    if(strncmp(buffer," ",1)==0)continue;

    if(
      match(buffer,"SLCF") == 1||
      match(buffer,"SLCC") == 1||
      match(buffer, "SLCD") == 1 ||
      match(buffer, "SLFL") == 1 ||
      match(buffer,"SLCT") == 1
      ){
      nsliceinfo++;
      continue;
    }
    if(match(buffer,"BNDF") == 1||
       match(buffer,"BNDC") == 1
       ){
      nboundary_files++;
      continue;
    }
    if(
      match(buffer,"PL3D") == 1){
      nplot3dinfo++;
      continue;
    }
    if(match(buffer,"GRID") == 1){
      nmeshes++;
      continue;
    }
    if(match(buffer,"PDIM") == 1){
      ipdim++;
      continue;
    }
  }

  if(nmeshes!=ipdim){
    fprintf(stderr,"*** Error (fatal): number of GRID statements (%i) not equal to\n",nmeshes);
    fprintf(stderr,"                 number of PDIM statements (%i)\n",ipdim);
    exit(0);
  }

  // allocate memory for mesh info

  if(nmeshes>0&&nmeshes==ipdim){
    NewMemory((void **)&meshinfo,nmeshes*sizeof(meshdata));
  }
  smvcase->meshinfo = meshinfo;
  smvcase->nmeshes = nmeshes;

  // allocate memory for slice file info

  if(nsliceinfo>0){
    slice *slicei;
    int i;

    NewMemory((void **)&sliceinfo,nsliceinfo*sizeof(slice));
    for(i=0;i<nsliceinfo;i++){
      slicei = sliceinfo + i;
      slicei->file=NULL;
    }
  }
  smvcase->sliceinfo=sliceinfo;
  smvcase->nsliceinfo = nsliceinfo;

  // allocate memory for boundary file info

  if(nboundary_files>0){
    boundary *boundaryi;
    int i;

    NewMemory((void **)&boundaryinfo,nboundary_files*sizeof(boundary));
    for(i=0;i<nboundary_files;i++){
      boundaryi = boundaryinfo + i;
      boundaryi->file=NULL;
    }
  }
  smvcase->boundaryinfo = boundaryinfo;
  smvcase->nboundary_files = nboundary_files;

  // allocate memory for plot3d file info

  if(nplot3dinfo>0){
    NewMemory((void **)&plot3dinfo,nplot3dinfo*sizeof(plot3d));
  }
  smvcase->nplot3dinfo = nplot3dinfo;
  smvcase->plot3dinfo = plot3dinfo;

  islice=0;
  iplot3d=0;
  iboundary=0;
  ipdim=0;
  igrid=0;
  rewind(streamsmv);
  while(!feof(streamsmv)){
    if(fgets(buffer,255,streamsmv)==NULL)break;
    CheckMemory;
    if(stream_out==NULL&&strncmp(buffer," ",1)==0)continue;

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ GRID ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"GRID") == 1){
      meshdata *meshi;
      float *xp, *yp, *zp;
      int ibar, jbar, kbar;

      meshi=meshinfo+igrid;
      igrid++;
      fgets(buffer,255,streamsmv);
      sscanf(buffer,"%i %i %i",&ibar,&jbar,&kbar);
      NewMemory((void **)&xp,sizeof(float)*(ibar+1));
      NewMemory((void **)&yp,sizeof(float)*(jbar+1));
      NewMemory((void **)&zp,sizeof(float)*(kbar+1));
      meshi->ibar=ibar;
      meshi->jbar=jbar;
      meshi->kbar=kbar;
      meshi->xplt=xp;
      meshi->yplt=yp;
      meshi->zplt=zp;

      if(stream_out!=NULL){
        trim_back(buffer);
        fprintf(stream_out,"GRID\n%s\n",buffer);
      }

      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ PDIM ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"PDIM") == 1){
      meshdata *meshi;

      meshi=meshinfo+ipdim;
      ipdim++;
      fgets(buffer,255,streamsmv);
      sscanf(buffer,"%f %f %f %f %f %f",&meshi->xbar0,&meshi->xbar,&meshi->ybar0,&meshi->ybar,&meshi->zbar0,&meshi->zbar);
      if(stream_out!=NULL){
        trim_back(buffer);
        fprintf(stream_out,"PDIM\n%s\n",buffer);
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ TRNX ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"TRNX")==1){
      float *xpltcopy, *xplt;
      int ibar, idummy, nn;
      meshdata *meshi;

      if(stream_out!=NULL){
        trim_back(buffer);
        fprintf(stream_out,"%s\n",buffer);
      }
      itrnx++;
      meshi = meshinfo + itrnx - 1;
      xpltcopy=meshi->xplt;
      xplt = meshi->xplt;

      ibar=meshi->ibar;
      fgets(buffer,255,streamsmv);
      if(stream_out!=NULL){
        trim_back(buffer);
        fprintf(stream_out,"%s\n",buffer);
      }
      sscanf(buffer,"%i ",&idummy);
      for(nn=0;nn<idummy;nn++){
        fgets(buffer,255,streamsmv);
        if(stream_out!=NULL){
          trim_back(buffer);
          fprintf(stream_out,"%s\n",buffer);
        }
      }
      for(nn=0;nn<=ibar;nn++){
        fgets(buffer,255,streamsmv);
        if(stream_out!=NULL){
          trim_back(buffer);
          fprintf(stream_out,"%s\n",buffer);
        }
        sscanf(buffer,"%i %f",&idummy,xpltcopy);
        xpltcopy++;
      }
      meshi->dx=xplt[1]-xplt[0];

      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ TRNY ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"TRNY")==1){
      float *ypltcopy, *yplt;
      int jbar, idummy, nn;
      meshdata *meshi;

      if(stream_out!=NULL){
        trim_back(buffer);
        fprintf(stream_out,"%s\n",buffer);
      }
      itrny++;
      meshi = meshinfo + itrny - 1;
      yplt = meshi->yplt;
      ypltcopy=meshi->yplt;
      jbar=meshi->jbar;
      fgets(buffer,255,streamsmv);
      if(stream_out!=NULL){
        trim_back(buffer);
        fprintf(stream_out,"%s\n",buffer);
      }
      sscanf(buffer,"%i ",&idummy);
      for(nn=0;nn<idummy;nn++){
        fgets(buffer,255,streamsmv);
        if(stream_out!=NULL){
          trim_back(buffer);
          fprintf(stream_out,"%s\n",buffer);
        }
      }
      for(nn=0;nn<=jbar;nn++){
        fgets(buffer,255,streamsmv);
        if(stream_out!=NULL){
          trim_back(buffer);
          fprintf(stream_out,"%s\n",buffer);
        }
        sscanf(buffer,"%i %f",&idummy,ypltcopy);
        ypltcopy++;
      }
      meshi->dy=yplt[1]-yplt[0];
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ TRNZ ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"TRNZ")==1){
      float *zpltcopy,*zplt;
      int kbar, idummy, nn;
      meshdata *meshi;

      if(stream_out!=NULL){
        trim_back(buffer);
        fprintf(stream_out,"%s\n",buffer);
      }
      itrnz++;
      meshi = meshinfo + itrnz - 1;
      zplt = meshi->zplt;
      zpltcopy=meshi->zplt;
      kbar=meshi->kbar;
      fgets(buffer,255,streamsmv);
      if(stream_out!=NULL){
        trim_back(buffer);
        fprintf(stream_out,"%s\n",buffer);
      }
      sscanf(buffer,"%i ",&idummy);
      for(nn=0;nn<idummy;nn++){
        fgets(buffer,255,streamsmv);
        if(stream_out!=NULL){
          trim_back(buffer);
          fprintf(stream_out,"%s\n",buffer);
        }
      }
      for(nn=0;nn<=kbar;nn++){
        fgets(buffer,255,streamsmv);
        if(stream_out!=NULL){
          trim_back(buffer);
          fprintf(stream_out,"%s\n",buffer);
        }
        sscanf(buffer,"%i %f",&idummy,zpltcopy);
        zpltcopy++;
      }
      meshi->dz = zplt[1]-zplt[0];
      continue;
    }
    if(match(buffer,"ENDF") == 1){
      char endian_filename[1024];
      FILE *ENDIANfile;
      int endian=0, endian_native, endian_data, len;

      if(fgets(buffer,255,streamsmv)==NULL)break;
      len=strlen(buffer);
      buffer[len-1]='\0';
      trim_back(buffer);
      fullfile(endian_filename,smvcase->dir,buffer);
      ENDIANfile = fopen(endian_filename,"rb");
      if(ENDIANfile!=NULL){
        endian_native = getendian();
        FSEEK(ENDIANfile,4,SEEK_SET);
        fread(&endian_data,4,1,ENDIANfile);
        fclose(ENDIANfile);
        endian=endian_native;
        if(endian_data!=1)endian=1-endian_native;
        smvcase->endian=endian;
      }
      if(stream_out!=NULL){
        int lenout;

        make_outfile(endian_filename, NULL, buffer, ".end");
        fprintf(stream_out,"ENDF\n %s\n",endian_filename);
        make_outfile(endian_filename, destdir, buffer, ".end");
        lenout=strlen(endian_filename);
        FORTendianout(endian_filename,lenout);
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ PL3D ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"PL3D") == 1){
      meshdata *plot3dmesh;
      plot3d *plot3di;
      float time_local;
      int meshnumber=1;
      char full_file[1024];
      FILE_SIZE filesize;

      if(strlen(buffer)>4){
        sscanf(buffer+4,"%f %i",&time_local,&meshnumber);
      }

      plot3dmesh = meshinfo + meshnumber - 1;

      plot3di = plot3dinfo + iplot3d;
      plot3di->plot3dmesh=plot3dmesh;
      plot3di->time=time_local;
      trim_back(buffer);
      strcpy(plot3di->keyword,buffer);

      fgets(buffer,255,streamsmv);
      fullfile(full_file,smvcase->dir,buffer);
      if(getfileinfo(full_file,NULL,&filesize)==0){
        int i;

        NewMemory((void **)&plot3di->file,(unsigned int)(strlen(full_file)+1));
        for(i = 0; i < 5; i++){
          NewMemory((void **)&plot3di->histogram[i], sizeof(histogramdata));
          init_histogram(plot3di->histogram[i],NHIST_BUCKETS);
        }

        CheckMemory;
        strcpy(plot3di->file,trim_front(buffer));
        CheckMemory;
        if(readlabels(plot3di->labels+0,streamsmv)==2)break;
        if(readlabels(plot3di->labels+1,streamsmv)==2)break;
        if(readlabels(plot3di->labels+2,streamsmv)==2)break;
        if(readlabels(plot3di->labels+3,streamsmv)==2)break;
        if(readlabels(plot3di->labels+4,streamsmv)==2)break;

        CheckMemory;

        iplot3d++;
      }
      else{
        if(display_warnings==1)fprintf(stderr,"*** Warning: the file, %s, does not exist.\n",full_file);
        CheckMemory;
        if(readlabels(plot3di->labels+0,streamsmv)==2)break;
        if(readlabels(plot3di->labels+1,streamsmv)==2)break;
        if(readlabels(plot3di->labels+2,streamsmv)==2)break;
        if(readlabels(plot3di->labels+3,streamsmv)==2)break;
        if(readlabels(plot3di->labels+4,streamsmv)==2)break;
        nplot3dinfo--;
        smvcase->nplot3dinfo=nplot3dinfo;
      }
      continue;
    }

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SLCF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(
      match(buffer,"SLCF") == 1||
      match(buffer,"SLCC") == 1||
      match(buffer, "SLCD") == 1 ||
      match(buffer, "SLFL") == 1 ||
      match(buffer,"SLCT") == 1)
    {
      int version_local=0;
      int len;
      FILE_SIZE filesize;
      slice *slicei;
      int meshnumber=0;
      meshdata *slicemesh;
      char full_file[1024];

      len=strlen(buffer);
      if(len>4){
        sscanf(buffer+4,"%i %i",&meshnumber,&version_local);
      }

      slicei = sliceinfo + islice;

      slicemesh = smvcase->meshinfo+meshnumber-1;
      slicei->slicemesh = slicemesh;
      trim_back(buffer);

      strcpy(slicei->keyword,buffer);

      if(match(buffer,"SLCF") == 1){
        slicei->slicetype=1;
      }
      if(match(buffer,"SLCC") == 1||match(buffer, "SLCD") == 1){
          slicei->slicetype = 2;
      }
      if(match(buffer,"SLFL") == 1){
        slicei->slicetype=3;
      }
      if(match(buffer,"SLCT") == 1){
        slicei->slicetype=4;
      }

      slicei->version=version_local;

      if(fgets(buffer,255,streamsmv)==NULL)break;
      trim_back(buffer);
      if(strlen(buffer)==0)break;
      fullfile(full_file,smvcase->dir,buffer);
      if(getfileinfo(full_file,NULL,&filesize)==0){
        int is1=-1, is2=-1, js1=-1, js2=-1, ks1=-1, ks2=-1;
        int ni, nj, nk;
        int error, lenfile;
        int endian;
        float *xplt, *yplt, *zplt;

        endian=getendian();
        NewMemory((void **)&slicei->file,(unsigned int)(strlen(full_file)+1));
        NewMemory((void **)&slicei->histogram,sizeof(histogramdata));
        init_histogram(slicei->histogram,NHIST_BUCKETS);
        STRCPY(slicei->file, trim_front(buffer));
        if(readlabels(&slicei->label,streamsmv)==2){
          fprintf(stderr,"*** Warning: problem reading SLCF entry\n");
          break;
        }
        slicei->filesize=filesize;
        lenfile=strlen(full_file);
        FORTgetsliceparms(full_file,&is1,&is2,&js1,&js2,&ks1,&ks2,&ni,&nj,&nk,&slicei->volslice,&error,lenfile);
        slicei->is1=is1;
        slicei->is2=is2;
        slicei->js1=js1;
        slicei->js2=js2;
        slicei->ks1=ks1;
        slicei->ks2=ks2;
        xplt = slicemesh->xplt;
        yplt = slicemesh->yplt;
        zplt = slicemesh->zplt;
        slicei->xmin = xplt[is1];
        slicei->xmax = xplt[is2];
        slicei->ymin = yplt[js1];
        slicei->ymax = yplt[js2];
        slicei->zmin = zplt[ks1];
        slicei->zmax = zplt[ks2];
        slicei->slice2=NULL;

        islice++;
      }
      else{
        if(display_warnings==1)fprintf(stderr,"*** Warning: the file, %s, does not exist.\n",buffer);
        if(readlabels(&sliceinfo[islice].label,streamsmv)==2)break;
        nsliceinfo--;
        smvcase->nsliceinfo=nsliceinfo;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ BNDF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"BNDF") == 1||
       match(buffer,"BNDC") == 1
       ){
      int version_local=0;
      int len;
      FILE_SIZE filesize;
      boundary *boundaryi;
      int meshnumber=0;
      meshdata *boundarymesh;
      char full_file[1024];

      len=strlen(buffer);
      if(len>4){
        sscanf(buffer+4,"%i %i",&meshnumber,&version_local);
      }

      boundaryi = boundaryinfo + iboundary;

      boundarymesh = smvcase->meshinfo+meshnumber-1;
      boundaryi->boundarymesh = boundarymesh;
      trim_back(buffer);

      strcpy(boundaryi->keyword,buffer);

      boundaryi->version=version_local;

      if(match(buffer,"BNDF") == 1){
        boundaryi->boundarytype=1;
      }
      if(match(buffer,"BNDC") == 1){
        boundaryi->boundarytype=2;
      }

      if(fgets(buffer,255,streamsmv)==NULL)break;
      trim_back(buffer);
      if(strlen(buffer)==0)break;
      fullfile(full_file,smvcase->dir,buffer);
      if(getfileinfo(full_file,NULL,&filesize)==0){
        int lenfile, endian, npatches, error, boundaryunitnumber;

        NewMemory((void **)&boundaryi->file,(unsigned int)(strlen(full_file)+1));
        NewMemory((void **)&boundaryi->histogram,sizeof(histogramdata));
        init_histogram(boundaryi->histogram,NHIST_BUCKETS);
        STRCPY(boundaryi->file, trim_front(buffer));
        if(readlabels(&boundaryi->label,streamsmv)==2){
          fprintf(stderr,"*** Warning: problem reading BNDF entry\n");
          break;
        }
        boundaryi->filesize=filesize;
        lenfile=strlen(full_file);
        endian=getendian();
        boundaryunitnumber=15;
        FORTgetboundaryheader1(full_file,&boundaryunitnumber, &npatches, &error, lenfile);
        if(npatches>0){
          int *pi1, *pi2, *pj1, *pj2, *pk1, *pk2, *patchdir, *patch2index, *patchsize, *qoffset;
          int i;

          NewMemory((void **)&pi1,npatches*sizeof(int));
          NewMemory((void **)&pi2,npatches*sizeof(int));
          NewMemory((void **)&pj1,npatches*sizeof(int));
          NewMemory((void **)&pj2,npatches*sizeof(int));
          NewMemory((void **)&pk1,npatches*sizeof(int));
          NewMemory((void **)&pk2,npatches*sizeof(int));
          NewMemory((void **)&patchdir,npatches*sizeof(int));
          NewMemory((void **)&patch2index,npatches*sizeof(int));
          NewMemory((void **)&patchsize,npatches*sizeof(int));
          NewMemory((void **)&qoffset,npatches*sizeof(int));
          boundaryi->pi1=pi1;
          boundaryi->pi2=pi2;
          boundaryi->pj1=pj1;
          boundaryi->pj2=pj2;
          boundaryi->pk1=pk1;
          boundaryi->pk2=pk2;
          boundaryi->patchdir=patchdir;
          boundaryi->patch2index=patch2index;
          boundaryi->npatches=npatches;
          boundaryi->patchsize=patchsize;
          boundaryi->qoffset=qoffset;
          FORTgetboundaryheader2(&boundaryunitnumber, &version_local, &npatches, pi1, pi2, pj1, pj2, pk1, pk2, patchdir);
          for(i=0;i<npatches;i++){
            boundaryi->patchsize[i] = (pi2[i]+1-pi1[i])*(pj2[i]+1-pj1[i])*(pk2[i]+1-pk1[i]);
            if(i==0){
              boundaryi->qoffset[i]=0;
            }
            else{
              boundaryi->qoffset[i]=boundaryi->qoffset[i-1]+boundaryi->patchsize[i-1];
            }
          }
          CheckMemory;
       }

        boundaryi->boundary2=NULL;

        iboundary++;
      }
      else{
        fprintf(stderr,"*** Warning: the file, %s, does not exist.\n",buffer);
        if(readlabels(&boundaryinfo[iboundary].label,streamsmv)==2)break;
        nboundary_files--;
        smvcase->nboundary_files=nboundary_files;
      }
      continue;
    }

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ vis keywords not differenced++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

    // skip over the following keywords

    if(
       match(buffer,"ISOF") == 1||
       match(buffer,"ISOG") == 1||
       match(buffer,"TISOF")==1||
       match(buffer,"SMOKE3D")==1||
       match(buffer,"SMOKF3D")==1||
       match(buffer,"VSMOKF3D")==1||
       match(buffer,"PART")==1||
       match(buffer,"EVAC")==1||
       match(buffer,"PRT5")==1||
       match(buffer,"EVA5")==1
       ){
      char comm[1024];

      strcpy(comm,buffer);
      fgets(buffer,255,streamsmv);
      if(match(comm,"PRT5")==1||match(comm,"EVA5")==1){
        int i, nlines;

        fgets(buffer,255,streamsmv);
        sscanf(buffer,"%i",&nlines);
        for(i=0;i<nlines;i++){
          fgets(buffer,255,streamsmv);
        }
      }
      else{
        fgets(buffer,255,streamsmv);
        fgets(buffer,255,streamsmv);
        fgets(buffer,255,streamsmv);
      }
      if(match(comm,"TISOF")==1){
        int i;

        for(i=0;i<3;i++){
          fgets(buffer,255,streamsmv);
        }
      }
      continue;
    }
    if(stream_out!=NULL){
      trim_back(buffer);
      fprintf(stream_out,"%s\n",buffer);
    }
    continue;
  }
  if(stream_out!=NULL){
    fprintf(stream_out,"SMOKEDIFF\n");
  }
  return 0;
}

