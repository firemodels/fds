#include "options.h"
#include "zlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "svzip.h"
#include "MALLOC.h"

int ReadSMV(char *smvfile){

  FILE *streamsmv;
  int ioffset;
  int unit_start=15;
  int igrid,ipdim;
  int ipatch,ipatch_seq;
  int iplot3d, iplot3d_seq;
  int ismoke3d, ismoke3d_seq;
  int islice, islice_seq;
#ifdef pp_PART
  int ipart_seq;
#endif
#define BUFFERSIZE 255
  char buffer[BUFFERSIZE];

  igrid=0;
  ipdim=0;
  streamsmv=fopen(smvfile,"r");
  if(streamsmv==NULL){
    PRINTF("The file: %s could not be opened\n",smvfile);
    return 1;
  }

  while(!feof(streamsmv)){
    if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
    CheckMemory;
    if(strncmp(buffer," ",1)==0)continue;

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SMOKE3D ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(Match(buffer,"SMOKE3D") == 1){
      nsmoke3dinfo++;
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ BNDF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(Match(buffer,"BNDF") == 1){
      npatchinfo++;
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ PL3D ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(Match(buffer,"PL3D") == 1){
      nplot3dinfo++;
      continue;
    }
#ifdef pp_PART
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ CLASS_OF_PARTICLES++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(Match(buffer,"CLASS_OF_PARTICLES") == 1){
      int i,nclasses;

      npartclassinfo++;
      fgets(buffer,BUFFERSIZE,streamsmv);
      fgets(buffer,BUFFERSIZE,streamsmv);
      fgets(buffer,BUFFERSIZE,streamsmv);
      sscanf(buffer,"%i",&nclasses);
      if(nclasses>0)maxpart5propinfo+=nclasses;
      for(i=0;i<nclasses;i++){
        fgets(buffer,BUFFERSIZE,streamsmv);
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ PART ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(Match(buffer,"PRT5") == 1){
      npartinfo++;
      continue;
    }
#endif
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SLCF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(
      Match(buffer,"SLCF") == 1||
      Match(buffer,"SLCC") == 1||
      Match(buffer, "SLCD") == 1 ||
      Match(buffer, "SLFL") == 1 ||
      Match(buffer,"SLCT") == 1
      ){
      nsliceinfo++;
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ GRID ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(Match(buffer,"GRID") == 1){
      nmeshes++;
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ PDIM ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(Match(buffer,"PDIM") == 1){
      ipdim++;
      continue;
    }
  }

  // allocate memory for smoke3d file info

  if(nsmoke3dinfo>0){
    smoke3d *smoke3di;
    int i;

    NewMemory((void **)&smoke3dinfo,nsmoke3dinfo*sizeof(smoke3d));
    for(i=0;i<nsmoke3dinfo;i++){
      smoke3di = smoke3dinfo + i;
      smoke3di->file=NULL;
      smoke3di->filebase=NULL;
    }
  }

  // allocate memory for boundary file info

  if(npatchinfo>0){
    patch *patchi;
    int i;

    NewMemory((void **)&patchinfo,npatchinfo*sizeof(patch));
    for(i=0;i<npatchinfo;i++){
      patchi = patchinfo + i;
      patchi->file=NULL;
      patchi->filebase=NULL;
      patchi->setvalmin=0;
      patchi->setvalmax=0;
      patchi->valmin=0.0;
      patchi->valmax=1.0;
    }
  }

  // allocate memory for plot3d file info

  if(nplot3dinfo>0){
    int i;

    NewMemory((void **)&plot3dinfo,nplot3dinfo*sizeof(plot3d));
    for(i=0;i<nplot3dinfo;i++){
      int j;
      plot3d *plot3di;

      plot3di = plot3dinfo + i;
      plot3di->file=NULL;
      plot3di->filebase=NULL;
      for(j=0;j<5;j++){
        plot3di->bounds[j].setvalmin=0;
        plot3di->bounds[j].setvalmax=0;
        plot3di->bounds[j].valmin=0.0;
        plot3di->bounds[j].valmax=1.0;
      }
    }
  }

  if(nmeshes>0&&nmeshes==ipdim){
    NewMemory((void **)&meshinfo,nmeshes*sizeof(meshdata));
  }
  else{
  }
  // allocate memory for slice file info

  if(nsliceinfo>0){
    slice *slicei;
    int i;

    NewMemory((void **)&sliceinfo,nsliceinfo*sizeof(slice));
    for(i=0;i<nsliceinfo;i++){
      slicei = sliceinfo + i;
      slicei->file=NULL;
      slicei->filebase=NULL;

      slicei->setvalmin=0;
      slicei->setvalmax=0;
      slicei->valmax=1.0;
      slicei->valmin=0.0;

      slicei->setchopvalmin=0;
      slicei->setchopvalmax=0;
      slicei->chopvalmax=1.0;
      slicei->chopvalmin=0.0;

      slicei->doit=1;
    }
  }

  // allocate memory for particle file info

#ifdef pp_PART
  if(npartinfo>0){
    part *parti;
    int i;

    NewMemory((void **)&partinfo,npartinfo*sizeof(part));
    for(i=0;i<npartinfo;i++){
      parti = partinfo + i;
      parti->file=NULL;
      parti->filebase=NULL;
      parti->setvalmin=0;
      parti->setvalmax=0;
      parti->valmax=1.0;
      parti->valmin=0.0;
    }
  }
  if(npartclassinfo>0){
    NewMemory((void **)&partclassinfo,npartclassinfo*sizeof(partclassdata));
  }
  if(maxpart5propinfo>0){
    NewMemory((void **)&part5propinfo,maxpart5propinfo*sizeof(partpropdata));
  }
#endif

  // read in smv file a second time_local to compress files

  ioffset=0;
  ipatch=0;
  ipatch_seq=0;
#ifdef pp_PART
  ipart_seq=0;
  npartclassinfo=0;
  npart5propinfo=0;
  npartinfo=0;
#endif
  iplot3d=0;
  iplot3d_seq=0;
  islice=0;
  islice_seq=0;
  ipdim=0;
  igrid=0;
  ismoke3d=0;
  ismoke3d_seq=0;
  rewind(streamsmv);
#ifndef pp_THREAD
  if(GLOBcleanfiles==0)PRINTF("Compressing .bf, .iso, .s3d, and .sf data files referenced in %s\n\n",smvfile);
#endif
  if(GLOBcleanfiles==1){
    PRINTF("Removing compressed .bf, .iso, .s3d and .sf data files referenced in %s\n",smvfile);
    PRINTF("   (Each removal occurs only if the corresponding uncompressed file exists)\n\n");
  }
  while(!feof(streamsmv)){
    patch *patchi;

    if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
    CheckMemory;
    if(strncmp(buffer," ",1)==0)continue;
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ ENDF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(Match(buffer,"ENDF") == 1){
      FILE *endianstream;
      int one;

      GLOBendf=1;
      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      TrimBack(buffer);
      strcpy(GLOBendianfilebase,buffer);
      FREEMEMORY(GLOBendianfile);
      if(GLOBsourcedir==NULL){
        NewMemory((void **)&GLOBendianfile,strlen(buffer)+1);
        strcpy(GLOBendianfile,buffer);
      }
      else{
        int lendir=0;

        if(GLOBsourcedir!=NULL)lendir=strlen(GLOBsourcedir);
        NewMemory((void **)&GLOBendianfile,strlen(buffer)+lendir+1);
        strcpy(GLOBendianfile,GLOBsourcedir);
        strcat(GLOBendianfile,buffer);
      }
      endianstream=fopen(GLOBendianfile,"rb");
      if(endianstream!=NULL){
        FSEEK(endianstream,4,SEEK_CUR);
        fread(&one,4,1,endianstream);
        if(one==1){
          endianswitch=0;
        }
        else{
          endianswitch=1;
        }
        fclose(endianstream);
      }
      continue;
    }

    if(Match(buffer,"GRID") == 1){
      meshdata *meshi;

      meshi=meshinfo+igrid;
      igrid++;
      fgets(buffer,BUFFERSIZE,streamsmv);
      sscanf(buffer,"%i %i %i",&meshi->ibar,&meshi->jbar,&meshi->kbar);
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ PDIM ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(Match(buffer,"PDIM") == 1){
      meshdata *meshi;

      meshi=meshinfo+ipdim;
      ipdim++;
      fgets(buffer,BUFFERSIZE,streamsmv);
      sscanf(buffer,"%f %f %f %f %f %f",&meshi->xbar0,&meshi->xbar,&meshi->ybar0,&meshi->ybar,&meshi->zbar0,&meshi->zbar);
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SYST ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(Match(buffer,"SYST") == 1){
      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      GLOBsyst=1;
      TrimBack(buffer);
      if(Match(buffer,"SGI") == 1||Match(buffer,"AIX")==1){
        if(getendian()==0){
          endianswitch=1;
        }
        else{
          endianswitch=0;
        }
      }
      else{
        if(getendian()==0){
          endianswitch=0;
        }
        else{
          endianswitch=1;
        }
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ OFFSET ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(Match(buffer,"OFFSET") == 1){
      float dummy;

      ioffset++;
      fgets(buffer,255,streamsmv);
      sscanf(buffer,"%f %f %f",&dummy,&dummy,&dummy);
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SMOKE3D ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(Match(buffer,"SMOKE3D") == 1){
      smoke3d *smoke3di;
      FILE_SIZE filesize;
      int filelen;
      char *buffer2;

      smoke3di = smoke3dinfo + ismoke3d;
      smoke3di->unit_start=unit_start++;
      ismoke3d_seq++;
      smoke3di->seq_id = ismoke3d_seq;
      smoke3di->autozip = 0;
      smoke3di->inuse=0;
      smoke3di->compressed=0;
      smoke3di->smokemesh=meshinfo + ioffset - 1;

      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      TrimBack(buffer);
      buffer2=TrimFront(buffer);
      filelen=strlen(buffer2);
      if(GLOBsourcedir!=NULL){
        filelen+=strlen(GLOBsourcedir)+1;
      }
      if(filelen<=0)break;
      if(getfileinfo(buffer2,GLOBsourcedir,&filesize)==0){
        NewMemory((void **)&smoke3di->file,(unsigned int)(filelen+1));
        NewMemory((void **)&smoke3di->filebase,(unsigned int)(filelen+1));
        STRCPY(smoke3di->filebase,buffer2);
        if(GLOBsourcedir!=NULL){
          STRCPY(smoke3di->file,GLOBsourcedir);
          STRCAT(smoke3di->file,buffer2);
        }
        else{
          STRCPY(smoke3di->file,buffer2);
        }
        smoke3di->filesize=filesize;
        if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
        buffer2 = TrimFront(buffer);
        TrimBack(buffer2);
        if(strcmp(buffer2,"HRRPUV")==0){
          smoke3di->is_soot=0;
        }
        else{
          smoke3di->is_soot=1;
        }
        ismoke3d++;
      }
      else{
        fprintf(stderr,"*** Warning: the file, %s, does not exist.\n",buffer);
        nsmoke3dinfo--;
      }
      continue;
    }
#ifdef pp_PART
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++ CLASS_OF_PARTICLES +++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

    if(Match(buffer,"CLASS_OF_PARTICLES") == 1){
      partclassdata *partclassi;
      int j;
      char *percen;

      partclassi = partclassinfo + npartclassinfo;

      fgets(buffer,BUFFERSIZE,streamsmv);
      percen=strchr(buffer,'%');
      if(percen!=NULL)percen=0;
      TrimBack(buffer);
      NewMemory((void **)&partclassi->name,strlen(buffer)+1);
      strcpy(partclassi->name,buffer);

      fgets(buffer,BUFFERSIZE,streamsmv);

      fgets(buffer,BUFFERSIZE,streamsmv);
      sscanf(buffer,"%i",&partclassi->ntypes);
      if(partclassi->ntypes>0){
        NewMemory((void **)&partclassi->labels,partclassi->ntypes*sizeof(flowlabels));
        for(j=0;j<partclassi->ntypes;j++){
          flowlabels *labelj;
          partpropdata *part5propi;

          labelj = partclassi->labels+j;
          labelj->longlabel=NULL;
          labelj->shortlabel=NULL;
          labelj->unit=NULL;
          ReadLabels(labelj,streamsmv);
          part5propi=getpartprop(labelj->shortlabel);
          if(part5propi==NULL){
            part5propi = part5propinfo + npart5propinfo;
            part5propi->label.longlabel=labelj->longlabel;
            part5propi->label.shortlabel=labelj->shortlabel;
            part5propi->label.unit=labelj->unit;
            part5propi->setvalmin=0;
            part5propi->valmin=1.0;
            part5propi->setvalmax=0;
            part5propi->valmax=0.0;
            npart5propinfo++;
          }
        }
      }
      npartclassinfo++;
      continue;
    }

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ PART ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(Match(buffer,"PRT5") == 1){
      int version_local=0,meshindex=1;
      int i;
      char *buffer2;
      int len;
      part *parti;
      FILE_SIZE filesize;

      len=strlen(buffer);
      if(len>4){
        buffer2=buffer+4;
        sscanf(buffer2,"%i %i",&meshindex,&version_local);
        meshindex--;

      }

      parti = partinfo + npartinfo;
      parti->unit_start=unit_start++;
      parti->partmesh = meshinfo + meshindex;
      ipart_seq++;
      parti->seq_id = ipart_seq;
      parti->autozip = 0;
      parti->inuse=0;
      parti->compressed=0;
      parti->compressed2=0;
      parti->inuse_part2iso=0;

      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      TrimBack(buffer);
      buffer2=TrimFront(buffer);
      if(strlen(buffer2)==0)break;
      if(getfileinfo(buffer2,GLOBsourcedir,&filesize)==0){
        int lendir=0;

        if(GLOBsourcedir!=NULL)lendir=strlen(GLOBsourcedir);
        NewMemory((void **)&parti->file,(unsigned int)(strlen(buffer2)+lendir+1));
        NewMemory((void **)&parti->filebase,(unsigned int)(strlen(buffer2)+1));
        STRCPY(parti->filebase,buffer2);
        if(GLOBsourcedir!=NULL){
          STRCPY(parti->file,GLOBsourcedir);
          STRCAT(parti->file,buffer2);
        }
        else{
          STRCPY(parti->file,buffer2);
        }
        parti->filesize=filesize;
        npartinfo++;
      }
      else{
        fprintf(stderr,"*** Warning: the file, %s, does not exist.\n",buffer2);
      }
      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      sscanf(buffer,"%i",&parti->nclasses);
      if(parti->nclasses>0){
        NewMemory((void **)&parti->classptr,parti->nclasses*sizeof(partclassdata *));
      }
      else{
        parti->nclasses=0;
      }
      for(i=0;i<parti->nclasses;i++){
        int classindex;

        if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
        sscanf(buffer,"%i",&classindex);
        parti->classptr[i]=partclassinfo + classindex - 1;
      }
      continue;
    }
#endif
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ BNDF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(Match(buffer,"BNDF") == 1){
      int version_local=0,dummy;
      char *buffer2;
      int len;
      FILE_SIZE filesize;

      len=strlen(buffer);
      if(len>4){
        buffer2=buffer+4;
        sscanf(buffer2,"%i %i",&dummy,&version_local);
      }

      patchi = patchinfo + ipatch;
      patchi->unit_start = unit_start++;
      ipatch_seq++;
      patchi->seq_id = ipatch;
      patchi->autozip = 0;
      patchi->inuse=0;
      patchi->compressed=0;
      patchi->version=version_local;

      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      TrimBack(buffer);
      buffer2=TrimFront(buffer);
      if(strlen(buffer2)==0)break;
      if(getfileinfo(buffer2,GLOBsourcedir,&filesize)==0){
        int lendir=0;

        if(GLOBsourcedir!=NULL)lendir=strlen(GLOBsourcedir);
        NewMemory((void **)&patchi->file,(unsigned int)(strlen(buffer2)+lendir+1));
        NewMemory((void **)&patchi->filebase,(unsigned int)(strlen(buffer2)+1));
        STRCPY(patchi->filebase,buffer2);
        if(GLOBsourcedir!=NULL){
          STRCPY(patchi->file,GLOBsourcedir);
          STRCAT(patchi->file,buffer2);
        }
        else{
          STRCPY(patchi->file,buffer2);
        }
        if(ReadLabels(&patchi->label,streamsmv)==2){
          fprintf(stderr,"*** Warning: problem reading BNDF entry\n");
          break;
        }
        patchi->filesize=filesize;
        if(GLOBget_boundary_bounds==1){
          int npatches, error, boundaryunitnumber;
          FILE_SIZE lenfile;

          NewMemory((void **)&patchi->histogram,sizeof(histogramdata));
          patchi->histogram->buckets = NULL;
          patchi->histogram->buckets_2d = NULL;
          lenfile = strlen(patchi->file);
          boundaryunitnumber=15;
          FORTgetboundaryheader1(patchi->file,&boundaryunitnumber, &npatches, &error, lenfile);
          if(npatches>0){
            int *pi1, *pi2, *pj1, *pj2, *pk1, *pk2, *patchdir, *patchsize;
            int i;

            NewMemory((void **)&pi1,npatches*sizeof(int));
            NewMemory((void **)&pi2,npatches*sizeof(int));
            NewMemory((void **)&pj1,npatches*sizeof(int));
            NewMemory((void **)&pj2,npatches*sizeof(int));
            NewMemory((void **)&pk1,npatches*sizeof(int));
            NewMemory((void **)&pk2,npatches*sizeof(int));
            NewMemory((void **)&patchdir,npatches*sizeof(int));
            NewMemory((void **)&patchsize,npatches*sizeof(int));
            patchi->pi1=pi1;
            patchi->pi2=pi2;
            patchi->pj1=pj1;
            patchi->pj2=pj2;
            patchi->pk1=pk1;
            patchi->pk2=pk2;
            patchi->patchdir=patchdir;
            patchi->npatches=npatches;
            patchi->patchsize=patchsize;
            FORTgetboundaryheader2(&boundaryunitnumber, &version_local, &npatches, pi1, pi2, pj1, pj2, pk1, pk2, patchdir);
            for(i=0;i<npatches;i++){
              patchi->patchsize[i] = (pi2[i]+1-pi1[i])*(pj2[i]+1-pj1[i])*(pk2[i]+1-pk1[i]);
            }
            CheckMemory;
          }
        }
        ipatch++;
      }
      else{
        fprintf(stderr,"*** Warning: the file, %s, does not exist.\n",buffer);
        if(ReadLabels(&patchi->label,streamsmv)==2)break;
        npatchinfo--;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SLCF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(
      Match(buffer,"SLCF") == 1||
      Match(buffer,"SLCC") == 1||
      Match(buffer, "SLCD") == 1 ||
      Match(buffer, "SLFL") == 1 ||
      Match(buffer,"SLCT") == 1)
    {
      int version_local=0,dummy;
      char *buffer2;
      int len;
      FILE_SIZE filesize;
      slice *slicei;
      int blocknumber;

      len=strlen(buffer);
      if(len>4){
        buffer2=buffer+4;
        sscanf(buffer2,"%i %i",&dummy,&version_local);
      }

      if(nmeshes>1){
        blocknumber=ioffset-1;
      }
      else{
        blocknumber=0;
      }
      if(len>5){
        buffer2=buffer+4;
        sscanf(buffer2,"%i",&blocknumber);
        blocknumber--;
      }

      islice_seq++;
      slicei = sliceinfo + islice;
      slicei->blocknumber=blocknumber;
      slicei->unit_start=unit_start++;
      slicei->version=version_local;
      slicei->seq_id = islice_seq;
      slicei->autozip = 0;
      slicei->inuse=0;
      slicei->involuse=0;
      slicei->compressed=0;
      slicei->vol_compressed=0;

      if(GLOBget_slice_bounds==1){
        NewMemory((void **)&slicei->histogram,sizeof(histogramdata));
        slicei->histogram->buckets = NULL;
        slicei->histogram->buckets_2d = NULL;
      }

      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      TrimBack(buffer);
      buffer2=TrimFront(buffer);
      if(strlen(buffer2)==0)break;
      if(getfileinfo(buffer2,GLOBsourcedir,&filesize)==0){
        int lendir=0;

        if(GLOBsourcedir!=NULL)lendir=strlen(GLOBsourcedir);
        NewMemory((void **)&slicei->file,(unsigned int)(strlen(buffer2)+lendir+1));
        NewMemory((void **)&slicei->filebase,(unsigned int)(strlen(buffer2)+1));
        STRCPY(slicei->filebase,buffer2);
        if(GLOBsourcedir!=NULL){
          STRCPY(slicei->file,GLOBsourcedir);
          STRCAT(slicei->file,buffer2);
        }
        else{
          STRCPY(slicei->file,buffer2);
        }
        if(ReadLabels(&slicei->label,streamsmv)==2){
          fprintf(stderr,"*** Warning: problem reading SLCF entry\n");
          break;
        }
        slicei->filesize=filesize;
        islice++;
      }
      else{
        fprintf(stderr,"*** Warning: the file, %s, does not exist.\n",buffer2);
        if(ReadLabels(&sliceinfo[islice].label,streamsmv)==2)break;
        nsliceinfo--;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ PL3D ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(Match(buffer,"PL3D") == 1){
      int version_local=0;
      char *buffer2;
      FILE_SIZE filesize;
      plot3d *plot3di;
      int blocknumber;
      float time_local;
      int blocktemp;

      if(nmeshes>1){
        blocknumber=igrid-1;
      }
      else{
        blocknumber=0;
      }
      if(strlen(buffer)>5){
        buffer2 = buffer+5;
        blocktemp=1;
        sscanf(buffer2,"%s %f %i",buffer2,&time_local,&blocktemp);
        if(blocktemp>0&&blocktemp<=nmeshes)blocknumber = blocktemp-1;
      }
      else{
        time_local=-1.0;
      }

      plot3di = plot3dinfo + iplot3d;
      plot3di->unit_start=unit_start++;
      iplot3d_seq++;
      plot3di->seq_id = iplot3d;
      plot3di->autozip = 0;
      plot3di->version=version_local;
      plot3di->plot3d_mesh=meshinfo + blocknumber;
      plot3di->time=time_local;
      plot3di->inuse=0;
      plot3di->compressed=0;

      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      TrimBack(buffer);
      buffer2=TrimFront(buffer);
      if(strlen(buffer2)==0)break;
      if(getfileinfo(buffer2,GLOBsourcedir,&filesize)==0){
        int lendir=0;

        if(GLOBsourcedir!=NULL)lendir=strlen(GLOBsourcedir);
        NewMemory((void **)&plot3di->file,(unsigned int)(strlen(buffer2)+lendir+1));
        NewMemory((void **)&plot3di->filebase,(unsigned int)(strlen(buffer2)+1));
        STRCPY(plot3di->filebase,buffer2);
        if(GLOBsourcedir!=NULL){
          STRCPY(plot3di->file,GLOBsourcedir);
        }
        else{
          STRCPY(plot3di->file,"");
        }
        STRCAT(plot3di->file,buffer2);
        if(ReadLabels(&plot3di->labels[0],streamsmv)==2||
           ReadLabels(&plot3di->labels[1],streamsmv)==2||
           ReadLabels(&plot3di->labels[2],streamsmv)==2||
           ReadLabels(&plot3di->labels[3],streamsmv)==2||
           ReadLabels(&plot3di->labels[4],streamsmv)==2){
          fprintf(stderr,"*** Warning: problem reading PL3D entry\n");
          break;
        }
        plot3di->filesize=filesize;
        iplot3d++;
      }
      else{
        fprintf(stderr,"*** Warning: the file, %s, does not exist.\n",buffer);
        if(ReadLabels(&plot3dinfo[iplot3d].labels[0],streamsmv)==2)break;
        if(ReadLabels(&plot3dinfo[iplot3d].labels[1],streamsmv)==2)break;
        if(ReadLabels(&plot3dinfo[iplot3d].labels[2],streamsmv)==2)break;
        if(ReadLabels(&plot3dinfo[iplot3d].labels[3],streamsmv)==2)break;
        if(ReadLabels(&plot3dinfo[iplot3d].labels[4],streamsmv)==2)break;
        nplot3dinfo--;
      }
      continue;
    }
  }
  {
    int i;

    for(i=0;i<nmeshes;i++){
      meshdata *meshi;
      int ii, jj, kk;
      float *xplt, *yplt, *zplt;
      float *xpltcell, *ypltcell, *zpltcell;

      meshi = meshinfo + i;
      meshi->dx = (meshi->xbar-meshi->xbar0)/meshi->ibar;
      meshi->dy = (meshi->ybar-meshi->ybar0)/meshi->jbar;
      meshi->dz = (meshi->zbar-meshi->zbar0)/meshi->kbar;
      meshi->dxx = (meshi->xbar-meshi->xbar0)/65535;
      meshi->dyy = (meshi->ybar-meshi->ybar0)/65535;
      meshi->dzz = (meshi->zbar-meshi->zbar0)/65535;

      meshi->xplt=NULL;
      NewMemory((void **)&meshi->xplt,(meshi->ibar+1)*sizeof(float));
      xplt = meshi->xplt;
      for(ii=0;ii<meshi->ibar;ii++){
        xplt[ii] = meshi->xbar0 + ii*meshi->dx;
      }
      xplt[meshi->ibar]=meshi->xbar;
      CheckMemory;

      meshi->xpltcell=NULL;
      NewMemory((void **)&meshi->xpltcell,(meshi->ibar+2)*sizeof(float));
      xpltcell = meshi->xpltcell+1;
      for(ii=-1;ii<meshi->ibar+1;ii++){
        xpltcell[ii] = meshi->xbar0 + (ii+0.5)*meshi->dx;
      }
      CheckMemory;

      meshi->yplt=NULL;
      NewMemory((void **)&meshi->yplt,(meshi->jbar+1)*sizeof(float));
      yplt = meshi->yplt;
      for(jj=0;jj<meshi->jbar;jj++){
        yplt[jj] = meshi->ybar0 + jj*meshi->dy;
      }
      yplt[meshi->jbar]=meshi->ybar;
      CheckMemory;

      meshi->ypltcell=NULL;
      NewMemory((void **)&meshi->ypltcell,(meshi->jbar+2)*sizeof(float));
      ypltcell = meshi->ypltcell+1;
      for(ii=-1;ii<meshi->jbar+1;ii++){
        ypltcell[ii] = meshi->ybar0 + (ii+0.5)*meshi->dy;
      }
      CheckMemory;

      meshi->zplt=NULL;
      NewMemory((void **)&meshi->zplt,(meshi->kbar+1)*sizeof(float));
      zplt = meshi->zplt;
      for(kk=0;kk<meshi->kbar;kk++){
        zplt[kk] = meshi->zbar0 + kk*meshi->dz;
      }
      zplt[meshi->kbar]=meshi->zbar;
      CheckMemory;

      meshi->zpltcell=NULL;
      NewMemory((void **)&meshi->zpltcell,(meshi->kbar+2)*sizeof(float));
      zpltcell = meshi->zpltcell+1;
      for(ii=-1;ii<meshi->kbar+1;ii++){
        zpltcell[ii] = meshi->zbar0 + (ii+0.5)*meshi->dz;
      }
      CheckMemory;

    }
  }
  init_volrender();
  return 0;
}


/* ------------------ ReadINI ------------------------ */

void ReadINI(char *casenameini){
  char *smoketemp;
  char globalini[256],smokeviewini[256];
  char *globaliniptr, *smokeviewiniptr;

  smoketemp = getenv("SMOKEVIEWINI");
  if(smoketemp==NULL)smoketemp=getenv("smokeviewini");
  if(smoketemp==NULL)smoketemp=getenv("svini");
  if(smoketemp==NULL)smoketemp=getenv("SVINI");
  if(smoketemp!=NULL){
    strcpy(globalini,smoketemp);
    strcat(globalini,dirseparator);
    strcat(globalini,"smokeview.ini");
  }
  globaliniptr=globalini;
  smokeviewiniptr=smokeviewini;
  strcpy(smokeviewini,"smokeview.ini");
  if(smoketemp!=NULL)ReadINI2(globaliniptr);
  ReadINI2(smokeviewiniptr);
  ReadINI2(casenameini);
}

/* ------------------ ReadINI2 ------------------------ */

void ReadINI2(char *inifile){
  char buffer[255],buffer2[255];
  char *type_buffer;
  FILE *stream;
  patch *patchi;

  stream=fopen(inifile,"r");
  if(stream==NULL)return;

  while(!feof(stream)){
    if(fgets(buffer,BUFFERSIZE,stream)==NULL)break;

    if(Match(buffer,"V_SLICE")==1){
      int setslicemin, setslicemax;
      float slicemin, slicemax;
      slice *slicei;

      fgets(buffer,BUFFERSIZE,stream);
      strcpy(buffer2,"");
      sscanf(buffer,"%i %f %i %f %s",&setslicemin,&slicemin,&setslicemax,&slicemax,buffer2);
      type_buffer=TrimFront(buffer2);
      TrimBack(type_buffer);
      slicei=getslice(type_buffer);
      if(slicei!=NULL){
        slicei->setvalmax=setslicemax;
        slicei->setvalmin=setslicemin;
        slicei->valmax=slicemax;
        slicei->valmin=slicemin;
      }
      continue;
    }
    if(Match(buffer,"V_PLOT3D")==1){
      int nplot3d_vars;
      plot3d *plot3di;
      int i;

      if(plot3dinfo==NULL)continue;
      plot3di = plot3dinfo;

      fgets(buffer,BUFFERSIZE,stream);
      nplot3d_vars=5;
      sscanf(buffer,"%i",&nplot3d_vars);
      if(nplot3d_vars<0)nplot3d_vars=0;
      if(nplot3d_vars>5)nplot3d_vars=5;

      for(i=0;i<nplot3d_vars;i++){
        int iplot3d;
        int setvalmin, setvalmax;
        float valmin, valmax;

        fgets(buffer,BUFFERSIZE,stream);
        sscanf(buffer,"%i %i %f %i %f",&iplot3d,&setvalmin,&valmin,&setvalmax,&valmax);
        iplot3d--;
        if(iplot3d>=0&&iplot3d<5){
          plot3di->bounds[iplot3d].setvalmin=setvalmin;
          plot3di->bounds[iplot3d].setvalmax=setvalmax;
          plot3di->bounds[iplot3d].valmin=valmin;
          plot3di->bounds[iplot3d].valmax=valmax;
        }
      }
      continue;
    }
    if(Match(buffer,"C_SLICE")==1){
      int setchopslicemin, setchopslicemax;
      float chopslicemin, chopslicemax;
      slice *slicei;

      fgets(buffer,BUFFERSIZE,stream);
      strcpy(buffer2,"");
      sscanf(buffer,"%i %f %i %f %s",&setchopslicemin,&chopslicemin,&setchopslicemax,&chopslicemax,buffer2);
      type_buffer=TrimFront(buffer2);
      TrimBack(type_buffer);
      slicei=getslice(type_buffer);
      if(slicei!=NULL){
        slicei->setchopvalmax=setchopslicemax;
        slicei->setchopvalmin=setchopslicemin;
        slicei->chopvalmax=chopslicemax;
        slicei->chopvalmin=chopslicemin;
      }
      continue;
    }
    if(Match(buffer,"V_BOUNDARY")==1){
      int setpatchmin, setpatchmax;
      float patchmin, patchmax;

      fgets(buffer,BUFFERSIZE,stream);
      strcpy(buffer2,"");
      sscanf(buffer,"%i %f %i %f %s",&setpatchmin,&patchmin,&setpatchmax,&patchmax,buffer2);
      type_buffer=TrimFront(buffer2);
      TrimBack(type_buffer);
      patchi=getpatch(type_buffer);
      if(patchi!=NULL){
        patchi->setvalmax=setpatchmax;
        patchi->setvalmin=setpatchmin;
        patchi->valmax=patchmax;
        patchi->valmin=patchmin;
      }
      continue;
    }
    if(GLOBframeskip<1&&Match(buffer,"SLICEZIPSTEP")==1){
	    fgets(buffer,BUFFERSIZE,stream);
	    sscanf(buffer,"%i",&GLOBslicezipstep);
	    if(GLOBslicezipstep<1)GLOBslicezipstep=1;
      continue;
    }
    if(GLOBframeskip<1&&Match(buffer,"SMOKE3DZIPSTEP")==1){
	    fgets(buffer,BUFFERSIZE,stream);
	    sscanf(buffer,"%i",&GLOBsmoke3dzipstep);
	    if(GLOBsmoke3dzipstep<1)GLOBsmoke3dzipstep=1;
      continue;
    }
    if(GLOBframeskip<1&&Match(buffer,"BOUNDZIPSTEP")==1){
	    fgets(buffer,BUFFERSIZE,stream);
	    sscanf(buffer,"%i",&GLOBboundzipstep);
	    if(GLOBboundzipstep<1)GLOBboundzipstep=1;
      continue;
    }

#ifdef pp_PART
    if(Match(buffer,"V_PARTICLES")==1){
      int setpartmin, setpartmax;
      float partmin, partmax;
      partpropdata *partpropi;

      fgets(buffer,BUFFERSIZE,stream);
      strcpy(buffer2,"");
      sscanf(buffer,"%i %f %i %f %s",&setpartmin,&partmin,&setpartmax,&partmax,buffer2);
      type_buffer=TrimFront(buffer2);
      TrimBack(type_buffer);
      partpropi=getpartprop(type_buffer);
      if(partpropi!=NULL){
        partpropi->setvalmax=setpartmax;
        partpropi->setvalmin=setpartmin;
        partpropi->valmax=partmax;
        partpropi->valmin=partmin;
      }
      continue;
    }
#endif
    if(Match(buffer,"SLICEAUTO")==1){
      int nslice_auto=0;
      int i;
      int seq_id;

      fgets(buffer,BUFFERSIZE,stream);
      sscanf(buffer,"%i",&nslice_auto);
      for(i=0;i<nslice_auto;i++){
        fgets(buffer,BUFFERSIZE,stream);
        sscanf(buffer,"%i",&seq_id);
        get_startup_slice(seq_id);
      }
      continue;
    }
    if(Match(buffer,"S3DAUTO")==1){
      int n3dsmokes=0;
      int i;
      int seq_id;

      fgets(buffer,BUFFERSIZE,stream);
      sscanf(buffer,"%i",&n3dsmokes);
      for(i=0;i<n3dsmokes;i++){
        fgets(buffer,BUFFERSIZE,stream);
        sscanf(buffer,"%i",&seq_id);
        get_startup_smoke(seq_id);
      }
      continue;
    }
    if(Match(buffer,"PATCHAUTO")==1){
      int n3dsmokes=0;
      int i;
      int seq_id;

      fgets(buffer,BUFFERSIZE,stream);
      sscanf(buffer,"%i",&n3dsmokes);
      for(i=0;i<n3dsmokes;i++){
        fgets(buffer,BUFFERSIZE,stream);
        sscanf(buffer,"%i",&seq_id);
        get_startup_patch(seq_id);
      }
      continue;
    }
  }
  fclose(stream);
  return;

}

 /* ------------------ get_startup_patch ------------------------ */

  void get_startup_patch(int seq_id){
    int i;
    for(i=0;i<npatchinfo;i++){
      patch *patchi;

      patchi = patchinfo + i;
      if(patchi->seq_id==seq_id){
        patchi->autozip=1;
        return;
      }
    }
  }

 /* ------------------ get_startup_smoke3d ------------------------ */

  void get_startup_smoke(int seq_id){
    int i;
    for(i=0;i<nsmoke3dinfo;i++){
      smoke3d *smoke3di;

      smoke3di = smoke3dinfo + i;

      if(smoke3di->seq_id==seq_id){
        smoke3di->autozip=1;
        return;
      }
    }
  }

 /* ------------------ get_startup_slice ------------------------ */

  void get_startup_slice(int seq_id){
    int i;
    for(i=0;i<nsliceinfo;i++){
      slice *slicei;

      slicei = sliceinfo + i;

      if(slicei->seq_id==seq_id){
        slicei->autozip=1;
        return;
      }
    }
  }


/* ------------------ init_volrender ------------------------ */

void init_volrender(void){
  int i;

  nvolrenderinfo=0;
  for(i=0;i<nmeshes;i++){
    meshdata *meshi;
    volrenderdata *vr;

    meshi = meshinfo + i;
    vr = &(meshi->volrenderinfo);
    vr->rendermesh=meshi;
    vr->fire=NULL;
    vr->smoke=NULL;
  }
  for(i=0;i<nsliceinfo;i++){
    slice *slicei;
    char *shortlabel;
    int blocknumber;
    meshdata *meshi;
    volrenderdata *vr;
    int ni, nj, nk;

    slicei = sliceinfo + i;
    slicei->isvolslice=0;
    slicei->voltype=0;
    blocknumber = slicei->blocknumber;
    if(blocknumber<0||blocknumber>=nmeshes)continue;
    meshi = meshinfo + blocknumber;
    getsliceparms_c(slicei->file,&ni,&nj,&nk);

    if(ni!=meshi->ibar+1||nj!=meshi->jbar+1||nk!=meshi->kbar+1)continue;
    vr = &(meshi->volrenderinfo);
    shortlabel = slicei->label.shortlabel;

    if(STRCMP(shortlabel,"temp")==0){
      vr->fire=slicei;
     continue;
    }
    if(STRCMP(shortlabel,"rho_Soot")==0){
      vr->smoke=slicei;
      continue;
    }
  }
  nvolrenderinfo=0;
  for(i=0;i<nmeshes;i++){
    meshdata *meshi;
    volrenderdata *vr;

    meshi = meshinfo + i;
    vr = &(meshi->volrenderinfo);
    if(vr->smoke!=NULL){
      nvolrenderinfo++;
      vr->smoke->isvolslice=1;
      vr->smoke->voltype=1;
      if(vr->fire!=NULL){
        vr->fire->isvolslice=1;
        vr->fire->voltype=2;
      }
    }
  }
}
