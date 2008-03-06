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
  int ismoke3d, ismoke3d_seq;
  int islice, islice_seq;
#ifdef pp_PART
  int ipart,ipart_seq;
#endif
  char buffer[255];

  igrid=0;
  ipdim=0;
  streamsmv=fopen(smvfile,"r");
  if(streamsmv==NULL){
    printf("The file: %s could not be opened\n",smvfile);
    return 1;
  }

  while(!feof(streamsmv)){
    if(fgets(buffer,255,streamsmv)==NULL)break;
    CheckMemory;
    if(strncmp(buffer," ",1)==0)continue;

    if(match(buffer,"SMOKE3D",7) == 1){
      nsmoke3d_files++;
      continue;
    }
    if(match(buffer,"BNDF",4) == 1){
      npatch_files++;
      continue;
    }
#ifdef pp_PART
    if(match(buffer,"PART",4) == 1){
      npart_files++;
      continue;
    }
#endif
    if(match(buffer,"SLCF",4) == 1){
      nslice_files++;
      continue;
    }
    if(match(buffer,"ISOF",4) == 1){
      niso_files++;
      continue;
    }
    if(match(buffer,"GRID",4) == 1){
      nmeshes++;
      continue;
    }
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
#ifdef pp_LIGHT
      smoke3di->light_q_rect=NULL;
      smoke3di->smoke_mesh=NULL;
#endif
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
      patchi->valmax=1.0;
      patchi->valmin=0.0;
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
      slicei->doit=1;
    }
  }


  
  // allocate memory for slice file info

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
#endif
  
  // read in smv file a second time to compress files

  ipatch=0;
  ipatch_seq=0;
#ifdef pp_PART
  ipart=0;
  ipart_seq=0;
#endif
  islice=0;
  islice_seq=0;
  iiso=0;
  iiso_seq=0;
  ipdim=0;
  igrid=0;
  ismoke3d=0;
  ismoke3d_seq=0;
  rewind(streamsmv);
  if(cleanfiles==0)printf("Compressing .bf, .s3d, and .sf data files referenced in %s\n",smvfile);
  if(cleanfiles==1){
    printf("Removing compressed .bf, .s3d and .sf data files found in %s\n",smvfile);
    printf("   (Each removal occurs only if the corresponding uncompressed file exists)\n\n");
  }
  while(!feof(streamsmv)){
    patch *patchi;

    if(fgets(buffer,255,streamsmv)==NULL)break;
    CheckMemory;
    if(strncmp(buffer," ",1)==0)continue;

    if(match(buffer,"ENDF",4) == 1){
      FILE *endianstream;
      int one;

      endf=1;
      if(fgets(buffer,255,streamsmv)==NULL)break;
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

    if(doiso==1&&match(buffer,"GRID",4) == 1){
      mesh *meshi;

      meshi=meshinfo+igrid;
      igrid++;
      fgets(buffer,255,streamsmv);
      sscanf(buffer,"%i %i %i",&meshi->ibar,&meshi->jbar,&meshi->kbar);
      continue;
    }
    if(doiso==1&&match(buffer,"PDIM",4) == 1){
      mesh *meshi;

      meshi=meshinfo+ipdim;
      ipdim++;
      fgets(buffer,255,streamsmv);
      sscanf(buffer,"%f %f %f %f %f %f",&meshi->xbar0,&meshi->xbar,&meshi->ybar0,&meshi->ybar,&meshi->zbar0,&meshi->zbar);
      continue;
    }

    if(match(buffer,"SYST",4) == 1){
      if(fgets(buffer,255,streamsmv)==NULL)break;
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

    if(match(buffer,"SMOKE3D",7) == 1){
      smoke3d *smoke3di;
      int filesize;
      int filelen;
      int blocknumber;

      smoke3di = smoke3dinfo + ismoke3d;
#ifdef pp_LIGHT
      if(strlen(buffer)>8){
        blocknumber=1;
        sscanf(buffer+8,"%i",&blocknumber);
        blocknumber--;
        smoke3di->smoke_mesh=meshinfo + blocknumber;
      }
#endif
      ismoke3d_seq++;
      smoke3di->seq_id = ismoke3d_seq;
      smoke3di->autozip = 0;
      if(fgets(buffer,255,streamsmv)==NULL)break;
      trim(buffer);
      filelen=strlen(buffer);
      if(sourcedir!=NULL){
        filelen+=lensourcedir+1;
      }
      if(filelen<=0)break;
      if(getfileinfo(buffer,sourcedir,&filesize)==0){
        NewMemory((void **)&smoke3di->file,(unsigned int)(filelen+lensourcedir+1));
        NewMemory((void **)&smoke3di->filebase,(unsigned int)(filelen+1));
        STRCPY(smoke3di->filebase,buffer);
        if(sourcedir!=NULL){
          STRCPY(smoke3di->file,sourcedir);
          STRCAT(smoke3di->file,buffer);
        }
        else{
          STRCPY(smoke3di->file,buffer);
        }
        smoke3di->filesize=filesize;
#ifdef pp_LIGHT
        if(readlabels(&smoke3di->label,streamsmv)!=2){
          ismoke3d++;
        }
#else
        ismoke3d++;
#endif
      }
      else{
        printf("*** Warning: the file, %s, does not exist.\n",buffer);
        nsmoke3d_files--;
      }
#ifdef pp_LIGHT
      {
        flowlabels *label;

        label=&smoke3di->label;
        smoke3di->type=0;
        if(label->shortlabel!=NULL){
          if(label->shortlabel!=NULL&&strncmp(smoke3di->label.shortlabel,"soot",4)==0){
            smoke3di->type=1;
          }
          else if(strncmp(smoke3di->label.shortlabel,"hrrpuv",6)==0){
            smoke3di->type=2;
          }
          else if(strncmp(smoke3di->label.shortlabel,"water",5)==0){
            smoke3di->type=3;
          }
          else{
            smoke3di->type=1;
          }
        }
      }
#endif
      continue;
    }


#ifdef pp_PART
    if(match(buffer,"PART",4) == 1){
      int version=0,dummy;
      char *buffer2;
      int len;
      part *parti;
      int filesize;

      len=strlen(buffer);
      if(len>4){
        buffer2=buffer+4;
        sscanf(buffer2,"%i %i",&dummy,&version);
      }

      parti = partinfo + ipart;
      ipart_seq++;
      parti->seq_id = ipart_seq;
      parti->autozip = 0;

      if(fgets(buffer,255,streamsmv)==NULL)break;
      trim(buffer);
      if(strlen(buffer)<=0)break;
      if(getfileinfo(buffer,sourcedir,&filesize)==0){
        NewMemory((void **)&parti->file,(unsigned int)(strlen(buffer)+lensourcedir+1));
        NewMemory((void **)&parti->filebase,(unsigned int)(strlen(buffer)+1));
        STRCPY(parti->filebase,buffer);
        if(sourcedir!=NULL){
          STRCPY(parti->file,sourcedir);
          STRCAT(parti->file,buffer);
        }
        else{
          STRCPY(parti->file,buffer);
        }
        if(readlabels(&parti->label,streamsmv)==2){
          printf("*** Warning: problem reading BNDF entry\n");
          break;
        }
        parti->filesize=filesize;
        ipart++;
      }
      else{
        printf("*** Warning: the file, %s, does not exist.\n",buffer);
        if(readlabels(&partinfo[ipart].label,streamsmv)==2)break;
        npart_files--;
      }
      continue;
    }
#endif

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
      patchi->version=version;

      if(fgets(buffer,255,streamsmv)==NULL)break;
      trim(buffer);
      if(strlen(buffer)<=0)break;
      if(getfileinfo(buffer,sourcedir,&filesize)==0){
        NewMemory((void **)&patchi->file,(unsigned int)(strlen(buffer)+lensourcedir+1));
        NewMemory((void **)&patchi->filebase,(unsigned int)(strlen(buffer)+1));
        STRCPY(patchi->filebase,buffer);
        if(sourcedir!=NULL){
          STRCPY(patchi->file,sourcedir);
          STRCAT(patchi->file,buffer);
        }
        else{
          STRCPY(patchi->file,buffer);
        }
        if(readlabels(&patchi->label,streamsmv)==2){
          printf("*** Warning: problem reading BNDF entry\n");
          break;
        }
        patchi->filesize=filesize;
        ipatch++;
      }
      else{
        printf("*** Warning: the file, %s, does not exist.\n",buffer);
        if(readlabels(&patchinfo[ipatch].label,streamsmv)==2)break;
        npatch_files--;
      }
      continue;
    }
    if(match(buffer,"SLCF",4) == 1){
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

      if(fgets(buffer,255,streamsmv)==NULL)break;
      trim(buffer);
      if(strlen(buffer)<=0)break;
      strcpy(buffer_rle,buffer);
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
      if(slicei->rle==-1&&getfileinfo(buffer,sourcedir,&filesize)==0){
        NewMemory((void **)&slicei->file,(unsigned int)(strlen(buffer)+lensourcedir+1));
        NewMemory((void **)&slicei->filebase,(unsigned int)(strlen(buffer)+1));
        STRCPY(slicei->filebase,buffer);
        slicei->rle=0;
        if(sourcedir!=NULL){
          STRCPY(slicei->file,sourcedir);
          STRCAT(slicei->file,buffer);
        }
        else{
          STRCPY(slicei->file,buffer);
        }
        if(readlabels(&slicei->label,streamsmv)==2){
          printf("*** Warning: problem reading SLCF entry\n");
          break;
        }
        slicei->filesize=filesize;
        islice++;
      }
      if(slicei->rle==-1){
        printf("*** Warning: the file, %s, does not exist.\n",buffer);
        if(readlabels(&sliceinfo[islice].label,streamsmv)==2)break;
        nslice_files--;
      }
      continue;
    }
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

      if(fgets(buffer,255,streamsmv)==NULL)break;
      trim(buffer);
      if(strlen(buffer)<=0)break;
      if(getfileinfo(buffer,sourcedir,&filesize)==0){
        int filelen;

        filelen = strlen(buffer)+lensourcedir+1;
        NewMemory((void **)&isoi->file,filelen);
        NewMemory((void **)&isoi->filebase,strlen(buffer)+1);
        STRCPY(isoi->filebase,buffer);
        if(sourcedir!=NULL){
          STRCPY(isoi->file,sourcedir);
          STRCAT(isoi->file,buffer);
        }
        else{
          STRCPY(isoi->file,buffer);
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

#ifdef pp_LIGHT
    if(meshinfo!=NULL)light_delta=meshinfo->dx;
#endif
    for(i=0;i<nmeshes;i++){
      mesh *meshi;
      int ii, jj, kk;
      float *xplt, *yplt, *zplt;
#ifdef pp_LIGHT
      float dx, dy, dz;
#endif

      meshi = meshinfo + i;
      meshi->dx = (meshi->xbar-meshi->xbar0)/meshi->ibar;
      meshi->dy = (meshi->ybar-meshi->ybar0)/meshi->jbar;
      meshi->dz = (meshi->zbar-meshi->zbar0)/meshi->kbar;
#ifdef pp_LIGHT
      dx = meshi->xbar - meshi->xbar0;
      dy = meshi->ybar - meshi->ybar0;
      dz = meshi->zbar - meshi->zbar0;
      meshi->dxyzmax = sqrt(dx*dx+dy*dy+dz*dz);
      if(meshi->dx<light_delta)light_delta=meshi->dx;
      if(meshi->dy<light_delta)light_delta=meshi->dy;
      if(meshi->dz<light_delta)light_delta=meshi->dz;
#endif
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

      meshi->yplt=NULL;
      NewMemory((void **)&meshi->yplt,(meshi->jbar+1)*sizeof(float));
      yplt = meshi->yplt;
      for(jj=0;jj<meshi->jbar;jj++){
        yplt[jj] = meshi->ybar0 + jj*meshi->dy;
      }
      yplt[meshi->jbar]=meshi->ybar;
      CheckMemory;

      meshi->zplt=NULL;
      NewMemory((void **)&meshi->zplt,(meshi->kbar+1)*sizeof(float));
      zplt = meshi->zplt;
      for(kk=0;kk<meshi->kbar;kk++){
        zplt[kk] = meshi->zbar0 + kk*meshi->dz;
      }
      zplt[meshi->kbar]=meshi->zbar;
      CheckMemory;

#ifdef pp_LIGHT
      meshi->photon_cell=NULL;
      if(make_lighting_file==1){
        int nni, nnj, nnk;

        nni = meshi->ibar;
        nnj = meshi->jbar;
        nnk = meshi->kbar;

        NewMemory((void **)&meshi->photon_cell,nni*nnj*nnk*sizeof(float));
      }
#endif

    }
  }
  return 0;
}


/* ------------------ readini ------------------------ */

void readini(char *casenameini){
  char *smoketemp;
  char globalini[256];
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
  
  if(globalini!=NULL)readini2(globalini);
  readini2("smokeview.ini");
  readini2(casenameini);
#ifdef pp_LIGHT
  {
    char lightini[256];

    strcpy(lightini,casenameini);
    lightini[strlen(lightini)-4]=0;
    strcat(lightini,".lit");
    readini2(lightini);
  }
#endif
}

/* ------------------ readini2 ------------------------ */

void readini2(char *inifile){
  char buffer[255],buffer2[255];
  char *type_buffer;
  FILE *stream;
  patch *patchi;

  stream=fopen(inifile,"r");
  if(stream==NULL)return;

#ifdef pp_LIGHT
  FREEMEMORY(lightinfo);
  nlightinfo=0;

  // pass 1

  while(!feof(stream)){
    if(fgets(buffer,255,stream)==NULL)break;

    if(match(buffer,"L_POINT",7)==1){
      fgets(buffer,255,stream);
      nlightinfo++;
      continue;
    }
    if(match(buffer,"L_MOVEPOINT",7)==1){
      fgets(buffer,255,stream);
      nlightinfo++;
      continue;
    }
    if(match(buffer,"L_LINE",7)==1){
      fgets(buffer,255,stream);
      nlightinfo++;
      continue;
    }
    if(match(buffer,"L_REGION",7)==1){
      fgets(buffer,255,stream);
      nlightinfo++;
      continue;
    }
  }
  if(nlightinfo>0){
    NewMemory((void **)&lightinfo,nlightinfo*sizeof(lightdata));
    NewMemory((void **)&light_cdf,(nlightinfo+1)*sizeof(float));
  }
  nlightinfo=0;
  rewind(stream);

  // pass 2

#endif
  while(!feof(stream)){
    if(fgets(buffer,255,stream)==NULL)break;

#ifdef pp_LIGHT
    if(match(buffer,"L_PHOTONS",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&nphotons);
      if(nphotons<1)nphotons=NPHOTONS;
      continue;
    }
    if(match(buffer,"L_MINMAX",7)==1){
      float l_min=light_min, l_max=light_max;

      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f",&l_min,&l_max);
      if(l_min<=0.0||l_max<=0.0){
        printf("*** error: light flux bounds must be positive\n");
        printf("           using default values of:\n");
        printf("  light_min=%f light_max=%f\n",light_min,light_max);
        continue;
      }
      if(l_min>l_max){
        printf("*** error: light min must be smaller than light max\n");
        printf("           using default values of:\n");
        printf("  light_min=%f light_max=%f\n",light_min,light_max);
        continue;
      }
      light_min=l_min;
      light_max=l_max;
      continue;
    }
    if(match(buffer,"L_POINT",7)==1||match(buffer,"L_MOVEPOINT",11)==1){
      lightdata *lighti;

      lighti = lightinfo + nlightinfo;
      lighti->type=0;
      if(match(buffer,"L_MOVEPOINT",11)==1)lighti->move=1;
      lighti->dir=0;
      if(lighti->move==0){
        float *xyz;
  
        xyz = lighti->xyz1;
        fgets(buffer,255,stream);
        sscanf(buffer,"%f %f %f %f",xyz,xyz+1,xyz+2,&lighti->q);
      }
      else{
        float *xyz, *xyz2;

        xyz = lighti->xyz1;
        xyz2 = lighti->xyz2;
        fgets(buffer,255,stream);
        sscanf(buffer,"%f %f %f %f",&lighti->t1,  xyz,  xyz+1, xyz+2);
        fgets(buffer,255,stream);
        sscanf(buffer,"%f %f %f %f",&lighti->t2, xyz2, xyz2+1,xyz2+2);
        fgets(buffer,255,stream);
        sscanf(buffer,"%f",&lighti->q);
      }
      nlightinfo++;
      continue;
    }
    if(match(buffer,"L_LINE",7)==1){
      lightdata *lighti;
      float *xyz1, *xyz2;

      lighti = lightinfo + nlightinfo;
      lighti->type=1;
      lighti->dir=0;
      xyz1 = lighti->xyz1;
      xyz2 = lighti->xyz2;
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f %f %f %f %f",xyz1,xyz1+1,xyz1+2,xyz2,xyz2+1,xyz2+2,&lighti->q);
      nlightinfo++;
      continue;
    }
    if(match(buffer,"L_REGION",7)==1){
      lightdata *lighti;
      float *xyz1, *xyz2;
      int dir=0;

      lighti = lightinfo + nlightinfo;
      lighti->type=2;
      xyz1 = lighti->xyz1;
      xyz2 = lighti->xyz2;
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f %f %f %f %f %i",xyz1,xyz1+1,xyz1+2,xyz2,xyz2+1,xyz2+2,&lighti->q,&dir);
      if(abs(dir)<0||abs(dir)>3)dir=0;
      lighti->dir=dir;
      nlightinfo++;
      continue;
    }
    if(match(buffer,"L_DELTA",7)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f",&light_delta);
      continue;
    }
#endif
    if(match(buffer,"V_SLICE",7)==1){
      int setslicemin, setslicemax;
      float slicemin, slicemax;
      slice *slicei;

      fgets(buffer,255,stream);
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
    if(match(buffer,"V_BOUNDARY",8)==1){
      int setpatchmin, setpatchmax;
      float patchmin, patchmax;

      fgets(buffer,255,stream);
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
    if(match(buffer,"SLICEZIPSTEP",12)==1){
	    fgets(buffer,255,stream);
	    sscanf(buffer,"%i",&slicezipstep);
	    if(slicezipstep<1)slicezipstep=1;
      continue;
    }
    if(match(buffer,"ISOZIPSTEP",10)==1){
	    fgets(buffer,255,stream);
	    sscanf(buffer,"%i",&isozipstep);
	    if(isozipstep<1)isozipstep=1;
      continue;
    }
    if(match(buffer,"SMOKE3DZIPSTEP",14)==1){
	    fgets(buffer,255,stream);
	    sscanf(buffer,"%i",&smoke3dzipstep);
	    if(smoke3dzipstep<1)smoke3dzipstep=1;
      continue;
    }
    if(match(buffer,"BOUNDZIPSTEP",12)==1){
	    fgets(buffer,255,stream);
	    sscanf(buffer,"%i",&boundzipstep);
	    if(boundzipstep<1)boundzipstep=1;
      continue;
    }

#ifdef pp_PART
    if(match(buffer,"V_PARTICLES",11)==1){
      int setpartmin, setpartmax;
      float partmin, partmax;
      part *parti;

      fgets(buffer,255,stream);
      strcpy(buffer2,"");
      sscanf(buffer,"%i %f %i %f %s",&setpartmin,&partmin,&setpartmax,&partmax,buffer2);
      type_buffer=trim_front(buffer2);
      trim(type_buffer);
      parti=getpart(type_buffer);
      if(parti!=NULL){
        parti->setvalmax=setpartmax;
        parti->setvalmin=setpartmin;
        parti->valmax=partmax;
        parti->valmin=partmin;
      }
      continue;
    }
#endif
    if(match(buffer,"SLICEAUTO",9)==1){
      int nslice_auto=0;
      int i;
      int seq_id;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&nslice_auto);
      for(i=0;i<nslice_auto;i++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&seq_id);
        get_startup_slice(seq_id);
      }
      continue;
    }
    if(match(buffer,"ISOAUTO",7)==1){
      int n3dsmokes=0;
      int i;
      int seq_id;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&n3dsmokes);
      for(i=0;i<n3dsmokes;i++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&seq_id);
        get_startup_iso(seq_id);
      }
      continue;
    }
    if(match(buffer,"S3DAUTO",7)==1){
      int n3dsmokes=0;
      int i;
      int seq_id;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&n3dsmokes);
      for(i=0;i<n3dsmokes;i++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&seq_id);
        get_startup_smoke(seq_id);
      }
      continue;
    }
    if(match(buffer,"PATCHAUTO",9)==1){
      int n3dsmokes=0;
      int i;
      int seq_id;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&n3dsmokes);
      for(i=0;i<n3dsmokes;i++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&seq_id);
        get_startup_patch(seq_id);
      }
      continue;
    }
  }
#ifdef pp_LIGHT
  {
    int i;

    if(nlightinfo>0){
      light_cdf[0]=0.0;
      for(i=1;i<nlightinfo+1;i++){
        lightdata *lighti;

        lighti = lightinfo + i - 1;
        light_cdf[i]=light_cdf[i-1]+lighti->q;
      }
      for(i=0;i<nlightinfo+1;i++){
        light_cdf[i]/=light_cdf[nlightinfo];
      }
    }
  }
#endif
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

