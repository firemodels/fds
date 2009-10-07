// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "svdiff.h"
#include "MALLOC.h"

// svn revision character string
char readsmv_revision[]="$Revision$";

/* ------------------ readsmv ------------------------ */

int readsmv(FILE *streamsmv, FILE *stream_out, casedata *smvcase){
  
  int igrid,ipdim;
  int islice,iplot3d;
  char buffer[255];
  mesh *meshinfo;
  slice *sliceinfo;
  plot3d *plot3dinfo;
  int nmeshes, nslice_files, nplot3d_files;
  int itrnx, itrny, itrnz;

  igrid=0;
  ipdim=0;
  nmeshes=0;
  nslice_files=0;
  nplot3d_files=0;
  itrnx=0;
  itrny=0;
  itrnz=0;

  while(!feof(streamsmv)){
    if(fgets(buffer,255,streamsmv)==NULL)break;
    CheckMemory;
    if(strncmp(buffer," ",1)==0)continue;

    if(
      match(buffer,"SLCF",4) == 1||
      match(buffer,"SLCC",4) == 1||
      match(buffer,"SLFL",4) == 1||
      match(buffer,"SLCT",4) == 1
      ){
      nslice_files++;
      continue;
    }
    if(
      match(buffer,"PL3D",4) == 1){
      nplot3d_files++;
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

  if(nmeshes!=ipdim){
    printf("*** fatal error: number of GRID statements (%i) not equal to\n",nmeshes);
    printf("                 number of PDIM statements (%i)\n",ipdim);
    exit(0);
  }

  // allocate memory for mesh info

  if(nmeshes>0&&nmeshes==ipdim){
    NewMemory((void **)&meshinfo,nmeshes*sizeof(mesh));
    smvcase->meshinfo=meshinfo;
    smvcase->nmeshes=nmeshes;
  }

  // allocate memory for slice file info

  if(nslice_files>0){
    slice *slicei;
    int i;

    NewMemory((void **)&sliceinfo,nslice_files*sizeof(slice));
    smvcase->nslice_files=nslice_files;
    smvcase->sliceinfo=sliceinfo;
    for(i=0;i<nslice_files;i++){
      slicei = sliceinfo + i;
      slicei->file=NULL;
    }
  }

  // allocate memory for plot3d file info

  if(nplot3d_files>0){
    NewMemory((void **)&plot3dinfo,nplot3d_files*sizeof(plot3d));
    smvcase->nplot3d_files=nplot3d_files;
    smvcase->plot3dinfo=plot3dinfo;
  }

  islice=0;
  iplot3d=0;
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
    if(match(buffer,"GRID",4) == 1){
      mesh *meshi;
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
        trim(buffer);
        fprintf(stream_out,"GRID\n%s\n",buffer);
      }
      
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
      fgets(buffer,255,streamsmv);
      sscanf(buffer,"%f %f %f %f %f %f",&meshi->xbar0,&meshi->xbar,&meshi->ybar0,&meshi->ybar,&meshi->zbar0,&meshi->zbar);
      if(stream_out!=NULL){
        trim(buffer);
        fprintf(stream_out,"PDIM\n%s\n",buffer);
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ TRNX ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"TRNX",4)==1){
      float *xpltcopy, *xplt;
      int ibar, idummy, nn;
      mesh *meshi;

      if(stream_out!=NULL){
        trim(buffer);
        fprintf(stream_out,"%s\n",buffer);
      }
      itrnx++;
      meshi = meshinfo + itrnx - 1;
      xpltcopy=meshi->xplt;
      xplt = meshi->xplt;

      ibar=meshi->ibar;
      fgets(buffer,255,streamsmv);
      if(stream_out!=NULL){
        trim(buffer);
        fprintf(stream_out,"%s\n",buffer);
      }
      sscanf(buffer,"%i ",&idummy);
      for(nn=0;nn<idummy;nn++){
        fgets(buffer,255,streamsmv);
        if(stream_out!=NULL){
          trim(buffer);
          fprintf(stream_out,"%s\n",buffer);
        }
      }
      for(nn=0;nn<=ibar;nn++){
        fgets(buffer,255,streamsmv);
        if(stream_out!=NULL){
          trim(buffer);
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
    if(match(buffer,"TRNY",4)==1){
      float *ypltcopy, *yplt;
      int jbar, idummy, nn;
      mesh *meshi;

      if(stream_out!=NULL){
        trim(buffer);
        fprintf(stream_out,"%s\n",buffer);
      }
      itrny++;
      meshi = meshinfo + itrny - 1;
      yplt = meshi->yplt;
      ypltcopy=meshi->yplt;
      jbar=meshi->jbar;
      fgets(buffer,255,streamsmv);
      if(stream_out!=NULL){
        trim(buffer);
        fprintf(stream_out,"%s\n",buffer);
      }
      sscanf(buffer,"%i ",&idummy);
      for(nn=0;nn<idummy;nn++){
        fgets(buffer,255,streamsmv);
        if(stream_out!=NULL){
          trim(buffer);
          fprintf(stream_out,"%s\n",buffer);
        }
      }
      for(nn=0;nn<=jbar;nn++){
        fgets(buffer,255,streamsmv);
        if(stream_out!=NULL){
          trim(buffer);
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
    if(match(buffer,"TRNZ",4)==1){
      float *zpltcopy,*zplt;
      int kbar, idummy, nn;
      mesh *meshi;

      if(stream_out!=NULL){
        trim(buffer);
        fprintf(stream_out,"%s\n",buffer);
      }
      itrnz++;
      meshi = meshinfo + itrnz - 1;
      zplt = meshi->zplt;
      zpltcopy=meshi->zplt;
      kbar=meshi->kbar;
      fgets(buffer,255,streamsmv);
      if(stream_out!=NULL){
        trim(buffer);
        fprintf(stream_out,"%s\n",buffer);
      }
      sscanf(buffer,"%i ",&idummy);
      for(nn=0;nn<idummy;nn++){
        fgets(buffer,255,streamsmv);
        if(stream_out!=NULL){
          trim(buffer);
          fprintf(stream_out,"%s\n",buffer);
        }
      }
      for(nn=0;nn<=kbar;nn++){
        fgets(buffer,255,streamsmv);
        if(stream_out!=NULL){
          trim(buffer);
          fprintf(stream_out,"%s\n",buffer);
        }
        sscanf(buffer,"%i %f",&idummy,zpltcopy);
        zpltcopy++;
      }
      meshi->dz = zplt[1]-zplt[0];
      continue;
    }
    if(match(buffer,"ENDF",4) == 1){
      char endianfilename[1024];
      FILE *ENDIANfile;
      int endian=0, endian_native, endian_data, len;

      if(fgets(buffer,255,streamsmv)==NULL)break;
      len=strlen(buffer);
      buffer[len-1]='\0';
      trim(buffer);
      fullfile(endianfilename,smvcase->dir,buffer);
      ENDIANfile = fopen(endianfilename,"rb");
      if(ENDIANfile!=NULL){
        endian_native = getendian();
        fseek(ENDIANfile,4,SEEK_SET);
        fread(&endian_data,4,1,ENDIANfile);
        fclose(ENDIANfile);
        endian=endian_native;
        if(endian_data!=1)endian=1-endian_native;
        smvcase->endian=endian;
      }
      if(stream_out!=NULL){
        make_outfile(endianfilename, NULL, buffer, ".end");
        fprintf(stream_out,"ENDF\n %s\n",buffer);
        make_outfile(endianfilename, destdir, buffer, ".end");
        ENDIANfile=fopen(endianfilename,"wb");
        if(ENDIANfile!=NULL){
          int one=1;

          fwrite(&one,4,1,ENDIANfile);
          fclose(ENDIANfile);
        }
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ PL3D ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"PL3D",4) == 1){
      mesh *plot3dmesh;
      plot3d *plot3di;
      float time;
      int meshnumber=1;

      if(strlen(buffer)>4){
        sscanf(buffer+4,"%f %i",&time,&meshnumber);
      }

      plot3dmesh = meshinfo + meshnumber - 1;

      plot3di = plot3dinfo + iplot3d;
      plot3di->plot3dmesh=plot3dmesh;
      plot3di->time=time;
      trim(buffer);
      strcpy(plot3di->keyword,buffer);

      fgets(buffer,255,streamsmv);
      trim(buffer);
      NewMemory((void **)&plot3di->file,(unsigned int)(strlen(buffer)+1));
      CheckMemory;
      strcpy(plot3di->file,buffer);
      CheckMemory;
      if(readlabels(plot3di->labels+0,streamsmv)==2)break;
      if(readlabels(plot3di->labels+1,streamsmv)==2)break;
      if(readlabels(plot3di->labels+2,streamsmv)==2)break;
      if(readlabels(plot3di->labels+3,streamsmv)==2)break;
      if(readlabels(plot3di->labels+4,streamsmv)==2)break;

      CheckMemory;
      
      iplot3d++;
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
      int version=0;
      int len;
      int filesize;
      slice *slicei;
      int meshnumber=0;
      mesh *slicemesh;
      float *xplt, *yplt, *zplt;
      char full_file[1024];

      len=strlen(buffer);
      if(len>4){
        sscanf(buffer+4,"%i %i",&meshnumber,&version);
      }

      slicei = sliceinfo + islice;

      slicemesh = smvcase->meshinfo+meshnumber-1;
      slicei->slicemesh = slicemesh;
      trim(buffer);

      strcpy(slicei->keyword,buffer);

      if(match(buffer,"SLCF",4) == 1){
        slicei->slicetype=1;
      }
      if(match(buffer,"SLCC",4) == 1){
        slicei->slicetype=2;
      }
      if(match(buffer,"SLFL",4) == 1){
        slicei->slicetype=3;
      }
      if(match(buffer,"SLCT",4) == 1){
        slicei->slicetype=4;
      }
      xplt = slicemesh->xplt;
      yplt = slicemesh->yplt;
      zplt = slicemesh->zplt;

      slicei->version=version;

      if(fgets(buffer,255,streamsmv)==NULL)break;
      trim(buffer);
      if(strlen(buffer)<=0)break;
      fullfile(full_file,smvcase->dir,buffer);
      if(getfileinfo(full_file,NULL,&filesize)==0){
        int is1, is2, js1, js2, ks1, ks2;
        int error, lenfile;
        int endian;

        NewMemory((void **)&slicei->file,(unsigned int)(strlen(buffer)+1));
        STRCPY(slicei->file,buffer);
        if(readlabels(&slicei->label,streamsmv)==2){
          printf("*** Warning: problem reading SLCF entry\n");
          break;
        }
        slicei->filesize=filesize;
        lenfile=strlen(full_file);
        FORTgetsliceparms(full_file,&endian,&is1,&is2,&js1,&js2,&ks1,&ks2,&slicei->volslice,&error,lenfile);
        slicei->is1=is1;
        slicei->is2=is2;
        slicei->js1=js1;
        slicei->js2=js2;
        slicei->ks1=ks1;
        slicei->ks2=ks2;
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
        printf("*** Warning: the file, %s, does not exist.\n",buffer);
        if(readlabels(&sliceinfo[islice].label,streamsmv)==2)break;
        nslice_files--;
        smvcase->nslice_files=nslice_files;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ vis keywords not differenced++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

    // skip over the following keywords

    if(match(buffer,"ISOF",4) == 1||
       match(buffer,"TISOF",5)==1||
       match(buffer,"SMOKE3D",7)==1||
       match(buffer,"BNDF",4)==1||
       match(buffer,"BNDC",4)==1||
       match(buffer,"PART",4)==1||
       match(buffer,"EVAC",4)==1||
       match(buffer,"PRT5",4)==1||
       match(buffer,"EVA5",4)==1
       ){
      char comm[1024];

      strcpy(comm,buffer);
      fgets(buffer,255,streamsmv);
      if(match(comm,"PRT5",4)==1||match(comm,"EVA5",4)==1){
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
      if(match(comm,"TISOF",5)==1){
        int i;

        for(i=0;i<3;i++){
          fgets(buffer,255,streamsmv);
        }
      }
      continue;
    }
    if(stream_out!=NULL){
      trim(buffer);
      fprintf(stream_out,"%s\n",buffer);
    }
    continue;
  }
  return 0;
}

