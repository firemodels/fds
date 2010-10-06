// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include "zlib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "egz_stdio.h"
#include "svzip.h"
#include "MALLOC.h"

// svn revision character string
char readfiles_revision[]="$Revision$";

int readsmv(char *smvfile){
  
  FILE *streamsmv;
  int iiso,igrid,ipdim, iiso_seq;
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
    printf("The file: %s could not be opened\n",smvfile);
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
    if(match(buffer,"SMOKE3D",7) == 1){
      nsmoke3d_files++;
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ BNDF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"BNDF",4) == 1){
      npatch_files++;
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ PL3D ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"PL3D",4) == 1){
      nplot3d_files++;
      continue;
    }
#ifdef pp_PART
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ CLASS_OF_PARTICLES++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"CLASS_OF_PARTICLES",18) == 1){
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
    if(match(buffer,"PRT5",4) == 1){
      npart_files++;
      continue;
    }
#endif
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SLCF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(
      match(buffer,"SLCF",4) == 1||
      match(buffer,"SLCC",4) == 1||
      match(buffer,"SLFL",4) == 1||
      match(buffer,"SLCT",4) == 1
      ){
      nslice_files++;
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ ISOF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"ISOF",4) == 1){
      niso_files++;
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ GRID ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"GRID",4) == 1){
      nmeshes++;
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ PDIM ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"PDIM",4) == 1){
      ipdim++;
      continue;
    }
  }

  // allocate memory for smoke3d file info

  if(nsmoke3d_files>0){
    smoke3d *smoke3di;
    int i;

    NewMemory((void **)&smoke3dinfo,nsmoke3d_files*sizeof(smoke3d));
    for(i=0;i<nsmoke3d_files;i++){
      smoke3di = smoke3dinfo + i;
      smoke3di->file=NULL;
      smoke3di->filebase=NULL;
    }
  }

  // allocate memory for boundary file info

  if(npatch_files>0){
    patch *patchi;
    int i;

    NewMemory((void **)&patchinfo,npatch_files*sizeof(patch));
    for(i=0;i<npatch_files;i++){
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

  if(nplot3d_files>0){
    int i;

    NewMemory((void **)&plot3dinfo,nplot3d_files*sizeof(plot3d));
    for(i=0;i<nplot3d_files;i++){
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
    NewMemory((void **)&meshinfo,nmeshes*sizeof(mesh));
    doiso=1;
  }
  else{
    doiso=0;
  }
  // allocate memory for slice file info

  if(nslice_files>0){
    slice *slicei;
    int i;

    NewMemory((void **)&sliceinfo,nslice_files*sizeof(slice));
    for(i=0;i<nslice_files;i++){
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

  // allocate memory for isosurface file info

  if(niso_files>0){
    iso *isoi;
    int i;

    NewMemory((void **)&isoinfo,niso_files*sizeof(iso));
    for(i=0;i<niso_files;i++){
      isoi = isoinfo + i;
      isoi->file=NULL;
      isoi->filebase=NULL;
      isoi->isolevels=NULL;
    }
  }

  // allocate memory for particle file info

#ifdef pp_PART
  if(npart_files>0){
    part *parti;
    int i;

    NewMemory((void **)&partinfo,npart_files*sizeof(part));
    for(i=0;i<npart_files;i++){
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
    NewMemory((void **)&partclassinfo,npartclassinfo*sizeof(part5class));
  }
  if(maxpart5propinfo>0){
    NewMemory((void **)&part5propinfo,maxpart5propinfo*sizeof(part5prop));
  }
#endif
  
  // read in smv file a second time to compress files

  ipatch=0;
  ipatch_seq=0;
#ifdef pp_PART
  ipart_seq=0;
  npartclassinfo=0;
  npart5propinfo=0;
  npart_files=0;
#endif
  iplot3d=0;
  iplot3d_seq=0;
  islice=0;
  islice_seq=0;
  iiso=0;
  iiso_seq=0;
  ipdim=0;
  igrid=0;
  ismoke3d=0;
  ismoke3d_seq=0;
  rewind(streamsmv);
  if(cleanfiles==0)printf("Compressing .bf, .iso, .s3d, and .sf data files referenced in %s\n",smvfile);
  if(cleanfiles==1){
    printf("Removing compressed .bf, .iso, .s3d and .sf data files referenced in %s\n",smvfile);
    printf("   (Each removal occurs only if the corresponding uncompressed file exists)\n\n");
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
    if(match(buffer,"ENDF",4) == 1){
      FILE *endianstream;
      int one;

      endf=1;
      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      trim(buffer);
      strcpy(endianfilebase,buffer);
      FREEMEMORY(endianfile);
      if(sourcedir==NULL){
        NewMemory((void **)&endianfile,strlen(buffer)+1);
        strcpy(endianfile,buffer);
      }
      else{
        NewMemory((void **)&endianfile,strlen(buffer)+strlen(sourcedir)+1);
        strcpy(endianfile,sourcedir);
        strcat(endianfile,buffer);
      }
      endianstream=fopen(endianfile,"rb");
      if(endianstream!=NULL){
        fseek(endianstream,4,SEEK_CUR);
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

    if(match(buffer,"GRID",4) == 1){
      mesh *meshi;

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
    if(match(buffer,"PDIM",4) == 1){
      mesh *meshi;

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
    if(match(buffer,"SYST",4) == 1){
      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      syst=1;
      trim(buffer);
      if(match(buffer,"SGI",3) == 1||match(buffer,"AIX",3)==1){
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
    ++++++++++++++++++++++ SMOKE3D ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"SMOKE3D",7) == 1){
      smoke3d *smoke3di;
      int filesize;
      int filelen;
      char *buffer2;

      smoke3di = smoke3dinfo + ismoke3d;
      ismoke3d_seq++;
      smoke3di->seq_id = ismoke3d_seq;
      smoke3di->autozip = 0;
      smoke3di->inuse=0;

      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      trim(buffer);
      buffer2=trim_front(buffer);
      filelen=strlen(buffer2);
      if(sourcedir!=NULL){
        filelen+=lensourcedir+1;
      }
      if(filelen<=0)break;
      if(getfileinfo(buffer2,sourcedir,&filesize)==0){
        NewMemory((void **)&smoke3di->file,(unsigned int)(filelen+lensourcedir+1));
        NewMemory((void **)&smoke3di->filebase,(unsigned int)(filelen+1));
        STRCPY(smoke3di->filebase,buffer2);
        if(sourcedir!=NULL){
          STRCPY(smoke3di->file,sourcedir);
          STRCAT(smoke3di->file,buffer2);
        }
        else{
          STRCPY(smoke3di->file,buffer2);
        }
        smoke3di->filesize=filesize;
        ismoke3d++;
      }
      else{
        printf("*** Warning: the file, %s, does not exist.\n",buffer);
        nsmoke3d_files--;
      }
      continue;
    }
#ifdef pp_PART
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++ CLASS_OF_PARTICLES +++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

    if(match(buffer,"CLASS_OF_PARTICLES",18) == 1){
      part5class *partclassi;
      int j;
      char *percen;

      partclassi = partclassinfo + npartclassinfo;

      fgets(buffer,BUFFERSIZE,streamsmv);
      percen=strchr(buffer,'%');
      if(percen!=NULL)percen=0;
      trim(buffer);
      NewMemory((void **)&partclassi->name,strlen(buffer)+1);
      strcpy(partclassi->name,buffer);

      fgets(buffer,BUFFERSIZE,streamsmv);

      fgets(buffer,BUFFERSIZE,streamsmv);
      sscanf(buffer,"%i",&partclassi->ntypes);
      if(partclassi->ntypes>0){
        NewMemory((void **)&partclassi->labels,partclassi->ntypes*sizeof(flowlabels));
        for(j=0;j<partclassi->ntypes;j++){
          flowlabels *labelj;
          part5prop *part5propi;

          labelj = partclassi->labels+j;
          labelj->longlabel=NULL;
          labelj->shortlabel=NULL;
          labelj->unit=NULL;
          readlabels(labelj,streamsmv);
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
    if(match(buffer,"PRT5",4) == 1){
      int version=0,meshindex=1;
      int i;
      char *buffer2;
      int len;
      part *parti;
      int filesize;

      len=strlen(buffer);
      if(len>4){
        buffer2=buffer+4;
        sscanf(buffer2,"%i %i",&meshindex,&version);
        meshindex--;

      }

      parti = partinfo + npart_files;
      parti->partmesh = meshinfo + meshindex;
      ipart_seq++;
      parti->seq_id = ipart_seq;
      parti->autozip = 0;
      parti->inuse=0;
      parti->inuse_part2iso=0;

      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      trim(buffer);
      buffer2=trim_front(buffer);
      if(strlen(buffer2)<=0)break;
      if(getfileinfo(buffer2,sourcedir,&filesize)==0){
        NewMemory((void **)&parti->file,(unsigned int)(strlen(buffer2)+lensourcedir+1));
        NewMemory((void **)&parti->filebase,(unsigned int)(strlen(buffer2)+1));
        STRCPY(parti->filebase,buffer2);
        if(sourcedir!=NULL){
          STRCPY(parti->file,sourcedir);
          STRCAT(parti->file,buffer2);
        }
        else{
          STRCPY(parti->file,buffer2);
        }
        parti->filesize=filesize;
        npart_files++;
      }
      else{
        printf("*** Warning: the file, %s, does not exist.\n",buffer2);
      }
      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      sscanf(buffer,"%i",&parti->nclasses);
      if(parti->nclasses>0){
        NewMemory((void **)&parti->classptr,parti->nclasses*sizeof(part5class *));
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
    if(match(buffer,"BNDF",4) == 1){
      int version=0,dummy;
      char *buffer2;
      int len;
      int filesize;

      len=strlen(buffer);
      if(len>4){
        buffer2=buffer+4;
        sscanf(buffer2,"%i %i",&dummy,&version);
      }

      patchi = patchinfo + ipatch;
      ipatch_seq++;
      patchi->seq_id = ipatch;
      patchi->autozip = 0;
      patchi->inuse=0;
      patchi->version=version;

      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      trim(buffer);
      buffer2=trim_front(buffer);
      if(strlen(buffer2)<=0)break;
      if(getfileinfo(buffer2,sourcedir,&filesize)==0){
        NewMemory((void **)&patchi->file,(unsigned int)(strlen(buffer2)+lensourcedir+1));
        NewMemory((void **)&patchi->filebase,(unsigned int)(strlen(buffer2)+1));
        STRCPY(patchi->filebase,buffer2);
        if(sourcedir!=NULL){
          STRCPY(patchi->file,sourcedir);
          STRCAT(patchi->file,buffer2);
        }
        else{
          STRCPY(patchi->file,buffer2);
        }
        if(readlabels(&patchi->label,streamsmv)==2){
          printf("*** Warning: problem reading BNDF entry\n");
          break;
        }
        patchi->filesize=filesize;
        if(get_boundary_bounds==1){
          int endian, npatches, error, boundaryunitnumber;
          FILE_SIZE lenfile;

          NewMemory((void **)&patchi->histogram,sizeof(histogramdata));
          lenfile=strlen(patchi->file);
          endian=getendian();
          boundaryunitnumber=15;
          FORTgetboundaryheader1(patchi->file,&boundaryunitnumber,&endian, &npatches, &error, lenfile);
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
            FORTgetboundaryheader2(&boundaryunitnumber, &version, &npatches, pi1, pi2, pj1, pj2, pk1, pk2, patchdir);
            for(i=0;i<npatches;i++){
              patchi->patchsize[i] = (pi2[i]+1-pi1[i])*(pj2[i]+1-pj1[i])*(pk2[i]+1-pk1[i]);
            }
            CheckMemory;
          }
        }
        ipatch++;
      }
      else{
        printf("*** Warning: the file, %s, does not exist.\n",buffer);
        if(readlabels(&patchinfo[ipatch].label,streamsmv)==2)break;
        npatch_files--;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SLCF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(
      match(buffer,"SLCF",4) == 1||
      match(buffer,"SLCC",4) == 1||
      match(buffer,"SLFL",4) == 1||
      match(buffer,"SLCT",4) == 1)
    {
      int version=0,dummy;
      char *buffer2;
      int len;
      int filesize;
      slice *slicei;
      char buffer_rle[255];


      len=strlen(buffer);
      if(len>4){
        buffer2=buffer+4;
        sscanf(buffer2,"%i %i",&dummy,&version);
      }

      islice_seq++;
      slicei = sliceinfo + islice;
      slicei->version=version;
      slicei->seq_id = islice_seq;
      slicei->autozip = 0;
      slicei->inuse=0;
        
      if(get_slice_bounds==1){
        NewMemory((void **)&slicei->histogram,sizeof(histogramdata));
      }

      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      trim(buffer);
      buffer2=trim_front(buffer);
      if(strlen(buffer2)<=0)break;
      strcpy(buffer_rle,buffer2);
      strcat(buffer_rle,".rle");
      slicei->rle=-1;
      if(getfileinfo(buffer_rle,sourcedir,&filesize)==0){
        NewMemory((void **)&slicei->file,(unsigned int)(strlen(buffer_rle)+lensourcedir+1));
        NewMemory((void **)&slicei->filebase,(unsigned int)(strlen(buffer_rle)+1));
        STRCPY(slicei->filebase,buffer_rle);
        slicei->rle=1;
        if(sourcedir!=NULL){
          STRCPY(slicei->file,sourcedir);
          STRCAT(slicei->file,buffer_rle);
        }
        else{
          STRCPY(slicei->file,buffer_rle);
        }
        if(readlabels(&slicei->label,streamsmv)==2){
          printf("*** Warning: problem reading SLCF entry\n");
          break;
        }
        slicei->filesize=filesize;
        islice++;
      }
      if(slicei->rle==-1&&getfileinfo(buffer2,sourcedir,&filesize)==0){
        NewMemory((void **)&slicei->file,(unsigned int)(strlen(buffer2)+lensourcedir+1));
        NewMemory((void **)&slicei->filebase,(unsigned int)(strlen(buffer2)+1));
        STRCPY(slicei->filebase,buffer2);
        slicei->rle=0;
        if(sourcedir!=NULL){
          STRCPY(slicei->file,sourcedir);
          STRCAT(slicei->file,buffer2);
        }
        else{
          STRCPY(slicei->file,buffer2);
        }
        if(readlabels(&slicei->label,streamsmv)==2){
          printf("*** Warning: problem reading SLCF entry\n");
          break;
        }
        slicei->filesize=filesize;
        islice++;
      }
      if(slicei->rle==-1){
        printf("*** Warning: the file, %s, does not exist.\n",buffer2);
        if(readlabels(&sliceinfo[islice].label,streamsmv)==2)break;
        nslice_files--;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ PL3D ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"PL3D",4) == 1){
      int version=0;
      char *buffer2;
      int filesize;
      plot3d *plot3di;
      int blocknumber;
      float time;
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
        sscanf(buffer2,"%s %f %i",buffer2,&time,&blocktemp);
        if(blocktemp>0&&blocktemp<=nmeshes)blocknumber = blocktemp-1;
      }
      else{
        time=-1.0;
      }

      plot3di = plot3dinfo + iplot3d;
      iplot3d_seq++;
      plot3di->seq_id = iplot3d;
      plot3di->autozip = 0;
      plot3di->version=version;
      plot3di->plot3d_mesh=meshinfo + blocknumber;
      plot3di->time=time;
      plot3di->inuse=0;

      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      trim(buffer);
      buffer2=trim_front(buffer);
      if(strlen(buffer2)<=0)break;
      if(getfileinfo(buffer2,sourcedir,&filesize)==0){
        NewMemory((void **)&plot3di->file,(unsigned int)(strlen(buffer2)+lensourcedir+1));
        NewMemory((void **)&plot3di->filebase,(unsigned int)(strlen(buffer2)+1));
        STRCPY(plot3di->filebase,buffer2);
        if(sourcedir!=NULL){
          STRCPY(plot3di->file,sourcedir);
        }
        else{
          STRCPY(plot3di->file,"");
        }
        STRCAT(plot3di->file,buffer2);
        if(readlabels(&plot3di->labels[0],streamsmv)==2||
           readlabels(&plot3di->labels[1],streamsmv)==2||
           readlabels(&plot3di->labels[2],streamsmv)==2||
           readlabels(&plot3di->labels[3],streamsmv)==2||
           readlabels(&plot3di->labels[4],streamsmv)==2){
          printf("*** Warning: problem reading PL3D entry\n");
          break;
        }
        plot3di->filesize=filesize;
        iplot3d++;
      }
      else{
        printf("*** Warning: the file, %s, does not exist.\n",buffer);
        if(readlabels(&plot3dinfo[iplot3d].labels[0],streamsmv)==2)break;
        if(readlabels(&plot3dinfo[iplot3d].labels[1],streamsmv)==2)break;
        if(readlabels(&plot3dinfo[iplot3d].labels[2],streamsmv)==2)break;
        if(readlabels(&plot3dinfo[iplot3d].labels[3],streamsmv)==2)break;
        if(readlabels(&plot3dinfo[iplot3d].labels[4],streamsmv)==2)break;
        nplot3d_files--;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ ISOF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"ISOF",4) == 1){

      int version=0;
      int blocknumber=0;
      char *buffer2;
      int len;
      int filesize;
      iso *isoi;

      CheckMemory;
      trim(buffer);
      len=strlen(buffer);
      if(len>4){
        buffer2=buffer+4;
        sscanf(buffer2,"%i %i",&blocknumber,&version);
        blocknumber--;
      }

      CheckMemory;
      isoi = isoinfo + iiso;
      iiso_seq++;
      isoi->seq_id = iiso_seq;
      isoi->autozip = 0;
      isoi->version=version;
      isoi->dataflag=0;
      isoi->blocknumber=blocknumber;
      isoi->file=NULL;
      isoi->filebase=NULL;
      isoi->inuse=0;

      if(fgets(buffer,BUFFERSIZE,streamsmv)==NULL)break;
      trim(buffer);
      buffer2=trim_front(buffer);
      if(strlen(buffer2)<=0)break;
      if(getfileinfo(buffer2,sourcedir,&filesize)==0){
        int filelen;

        filelen = strlen(buffer2)+lensourcedir+1;
        NewMemory((void **)&isoi->file,filelen);
        NewMemory((void **)&isoi->filebase,strlen(buffer2)+1);
        STRCPY(isoi->filebase,buffer2);
        if(sourcedir!=NULL){
          STRCPY(isoi->file,sourcedir);
          STRCAT(isoi->file,buffer2);
        }
        else{
          STRCPY(isoi->file,buffer2);
        }
        if(readlabels(&isoi->label,streamsmv)==2){
          printf("*** Warning: problem reading SLCF entry\n");
          break;
        }
        isoi->filesize=filesize;
        iiso++;
      }
      else{
        printf("*** Warning: the file, %s, does not exist.\n",buffer);
        if(readlabels(&isoinfo[iiso].label,streamsmv)==2)break;
        niso_files--;
      }
      continue;
    }

  }
  {
    int i;

    for(i=0;i<nmeshes;i++){
      mesh *meshi;
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
  return 0;
}


/* ------------------ readini ------------------------ */

void readini(char *casenameini){
  char *smoketemp;
  char globalini[256],smokeviewini[256];
  char *globaliniptr, *smokeviewiniptr;
#ifdef pp_cvf
char dirseparator[]="\\";
#else
char dirseparator[]="/";
#endif

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
  if(smoketemp!=NULL)readini2(globaliniptr);
  readini2(smokeviewiniptr);
  readini2(casenameini);
}

/* ------------------ readini2 ------------------------ */

void readini2(char *inifile){
  char buffer[255],buffer2[255];
  char *type_buffer;
  FILE *stream;
  patch *patchi;

  stream=fopen(inifile,"r");
  if(stream==NULL)return;

  while(!feof(stream)){
    if(fgets(buffer,BUFFERSIZE,stream)==NULL)break;

    if(match(buffer,"V_SLICE",7)==1){
      int setslicemin, setslicemax;
      float slicemin, slicemax;
      slice *slicei;

      fgets(buffer,BUFFERSIZE,stream);
      strcpy(buffer2,"");
      sscanf(buffer,"%i %f %i %f %s",&setslicemin,&slicemin,&setslicemax,&slicemax,buffer2);
      type_buffer=trim_front(buffer2);
      trim(type_buffer);
      slicei=getslice(type_buffer);
      if(slicei!=NULL){
        slicei->setvalmax=setslicemax;
        slicei->setvalmin=setslicemin;
        slicei->valmax=slicemax;
        slicei->valmin=slicemin;
      }
      continue;
    }
    if(match(buffer,"V_PLOT3D",8)==1){
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
    if(match(buffer,"C_SLICE",7)==1){
      int setchopslicemin, setchopslicemax;
      float chopslicemin, chopslicemax;
      slice *slicei;

      fgets(buffer,BUFFERSIZE,stream);
      strcpy(buffer2,"");
      sscanf(buffer,"%i %f %i %f %s",&setchopslicemin,&chopslicemin,&setchopslicemax,&chopslicemax,buffer2);
      type_buffer=trim_front(buffer2);
      trim(type_buffer);
      slicei=getslice(type_buffer);
      if(slicei!=NULL){
        slicei->setchopvalmax=setchopslicemax;
        slicei->setchopvalmin=setchopslicemin;
        slicei->chopvalmax=chopslicemax;
        slicei->chopvalmin=chopslicemin;
      }
      continue;
    }
    if(match(buffer,"V_BOUNDARY",8)==1){
      int setpatchmin, setpatchmax;
      float patchmin, patchmax;

      fgets(buffer,BUFFERSIZE,stream);
      strcpy(buffer2,"");
      sscanf(buffer,"%i %f %i %f %s",&setpatchmin,&patchmin,&setpatchmax,&patchmax,buffer2);
      type_buffer=trim_front(buffer2);
      trim(type_buffer);
      patchi=getpatch(type_buffer);
      if(patchi!=NULL){
        patchi->setvalmax=setpatchmax;
        patchi->setvalmin=setpatchmin;
        patchi->valmax=patchmax;
        patchi->valmin=patchmin;
      }
      continue;
    }
    if(frameskip<1&&match(buffer,"SLICEZIPSTEP",12)==1){
	    fgets(buffer,BUFFERSIZE,stream);
	    sscanf(buffer,"%i",&slicezipstep);
	    if(slicezipstep<1)slicezipstep=1;
      continue;
    }
    if(frameskip<1&&match(buffer,"ISOZIPSTEP",10)==1){
	    fgets(buffer,BUFFERSIZE,stream);
	    sscanf(buffer,"%i",&isozipstep);
	    if(isozipstep<1)isozipstep=1;
      continue;
    }
    if(frameskip<1&&match(buffer,"SMOKE3DZIPSTEP",14)==1){
	    fgets(buffer,BUFFERSIZE,stream);
	    sscanf(buffer,"%i",&smoke3dzipstep);
	    if(smoke3dzipstep<1)smoke3dzipstep=1;
      continue;
    }
    if(frameskip<1&&match(buffer,"BOUNDZIPSTEP",12)==1){
	    fgets(buffer,BUFFERSIZE,stream);
	    sscanf(buffer,"%i",&boundzipstep);
	    if(boundzipstep<1)boundzipstep=1;
      continue;
    }

#ifdef pp_PART
    if(match(buffer,"V_PARTICLES",11)==1){
      int setpartmin, setpartmax;
      float partmin, partmax;
      part5prop *partpropi;

      fgets(buffer,BUFFERSIZE,stream);
      strcpy(buffer2,"");
      sscanf(buffer,"%i %f %i %f %s",&setpartmin,&partmin,&setpartmax,&partmax,buffer2);
      type_buffer=trim_front(buffer2);
      trim(type_buffer);
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
    if(match(buffer,"SLICEAUTO",9)==1){
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
    if(match(buffer,"ISOAUTO",7)==1){
      int n3dsmokes=0;
      int i;
      int seq_id;

      fgets(buffer,BUFFERSIZE,stream);
      sscanf(buffer,"%i",&n3dsmokes);
      for(i=0;i<n3dsmokes;i++){
        fgets(buffer,BUFFERSIZE,stream);
        sscanf(buffer,"%i",&seq_id);
        get_startup_iso(seq_id);
      }
      continue;
    }
    if(match(buffer,"S3DAUTO",7)==1){
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
    if(match(buffer,"PATCHAUTO",9)==1){
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
    for(i=0;i<npatch_files;i++){
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
    for(i=0;i<nsmoke3d_files;i++){
      smoke3d *smoke3di;

      smoke3di = smoke3dinfo + i;

      if(smoke3di->seq_id==seq_id){
        smoke3di->autozip=1;
        return;
      }
    }
  }


 /* ------------------ get_startup_iso ------------------------ */

  void get_startup_iso(int seq_id){
    int i;
    for(i=0;i<niso_files;i++){
      iso *isoi;

      isoi = isoinfo + i;

      if(isoi->seq_id==seq_id){
        isoi->autozip=1;
        return;
      }
    }
  }

 /* ------------------ get_startup_slice ------------------------ */

  void get_startup_slice(int seq_id){
    int i;
    for(i=0;i<nslice_files;i++){
      slice *slicei;

      slicei = sliceinfo + i;

      if(slicei->seq_id==seq_id){
        slicei->autozip=1;
        return;
      }
    }
  }

