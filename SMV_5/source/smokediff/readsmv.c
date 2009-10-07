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
  int islice;
  char buffer[255];
  mesh *meshinfo;
  slice *sliceinfo;
  int nmeshes, nslice_files;
  int itrnx, itrny, itrnz;

  igrid=0;
  ipdim=0;
  nmeshes=0;
  nslice_files=0;
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
  islice=0;
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
       match(buffer,"PL3D",4)==1||
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
      if(match(comm,"PL3D",4)==1){
        int i;

        for(i=0;i<4*3;i++){
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

/* ------------------ setup_slice ------------------------ */

void setup_slice(FILE *stream_out){
  casedata *case1, *case2;
  int i;

  case1 = caseinfo;
  case2 = caseinfo + 1;

  for(i=0;i<case1->nslice_files;i++){
    slice *slicei;

    slicei = case1->sliceinfo + i;
    slicei->slice2 = getslice(slicei,case2);
    if(slicei->slice2!=NULL&&stream_out!=NULL){
      char outfile[1024];

      fprintf(stream_out,"%s\n",slicei->keyword);
      make_outfile(outfile,NULL,slicei->file,".sf");
      fprintf(stream_out,"%s\n",outfile);
      fprintf(stream_out,"%s\n",slicei->label.longlabel);
      fprintf(stream_out,"%s\n",slicei->label.shortlabel);
      fprintf(stream_out,"%s\n",slicei->label.unit);
    }
  }
}

/* ------------------ getslice ------------------------ */

slice *getslice(slice *slicein, casedata *case2){
  int i;
  float dx, dy, dz;

  dx = slicein->slicemesh->dx/2.0;
  dy = slicein->slicemesh->dy/2.0;
  dz = slicein->slicemesh->dz/2.0;

  for(i=0;i<case2->nslice_files;i++){
    slice *sliceout;

    sliceout = case2->sliceinfo + i;
    if(slicein->slicetype!=sliceout->slicetype)continue;
    if(strcmp(slicein->label.longlabel,sliceout->label.longlabel)!=0)continue;
    if(fabs(slicein->xmin-sliceout->xmin)>dx)continue;
    if(fabs(slicein->xmax-sliceout->xmax)>dx)continue;
    if(fabs(slicein->ymin-sliceout->ymin)>dy)continue;
    if(fabs(slicein->ymax-sliceout->ymax)>dy)continue;
    if(fabs(slicein->zmin-sliceout->zmin)>dz)continue;
    if(fabs(slicein->zmax-sliceout->zmax)>dz)continue;
    return sliceout;
  }
  return NULL;
}

/* ------------------ getslice ------------------------ */

void diff_slices(void){
  int j;

  for(j=0;j<caseinfo->nslice_files;j++){
    char *file1, *file2;
    char fullfile1[1024], fullfile2[1024], outfile[1024];
    slice *slicei, *slice1, *slice2;
    FILE *stream;
    int unit1, unit2, unit3;
    FILE_SIZE len1,len2;
    int is1a, is2a, js1a, js2a, ks1a, ks2a;
    int is1b, is2b, js1b, js2b, ks1b, ks2b;
    int error1,error2;
    float time1, *qframe1;
    int nqframe1;
    float time2, *qframe2;
    int nqframe2;
    float *qframeout;
    int i;
    int len;

    slicei = caseinfo->sliceinfo+j;
    slice1 = slicei;
    if(slicei->slice2==NULL)continue;
    slice2 = slicei->slice2;
    file1 = slicei->file;
    file2 = slicei->slice2->file;
    fullfile(fullfile1,sourcedir1,file1);
    fullfile(fullfile2,sourcedir2,file2);

    stream=fopen(fullfile1,"r");
    if(stream==NULL)continue;
    fclose(stream);

    stream=fopen(fullfile2,"r");
    if(stream==NULL)continue;
    fclose(stream);

    make_outfile(outfile,destdir,file1,".sf");
    if(strlen(outfile)==0)continue;
    stream=fopen(outfile,"w");
    if(stream==NULL)continue;
    fclose(stream);

    unit1=11;
    len1=strlen(fullfile1);
    FORTopenslice(fullfile1,&unit1,&caseinfo->endian,&is1a,&is2a,&js1a,&js2a,&ks1a,&ks2a,&error1,len1);
    unit2=12;
    len2=strlen(fullfile2);
    FORTopenslice(fullfile2,&unit2,&caseinfo->endian,&is1b,&is2b,&js1b,&js2b,&ks1b,&ks2b,&error2,len2);
    if(is1a!=is1b||js1a!=js1b||ks1a!=ks1b||
       is2a!=is2b||js2a!=js2b||ks2a!=ks2b||
       error1!=0||error2!=0){
      FORTclosefortranfile(&unit1);
      FORTclosefortranfile(&unit2);
      if(error1!=0||error2!=0){
        if(error1=0)printf("*** problem opening %s\n",fullfile1);
        if(error2=0)printf("*** problem opening %s\n",fullfile2);
      }
      if(is1a!=is1b||js1a!=js1b||ks1a!=ks1b||
         is2a!=is2b||js2a!=js2b||ks2a!=ks2b){
        printf("*** integer slice bounds do not match for\n",fullfile1);
        printf("    %i %i %i %i %i %i\n",is1a, is2a, js1a, js2a, ks1a, ks2a);
        printf("    %i %i %i %i %i %i\n",is1b, is2b, js1b, js2b, ks1b, ks2b);
        printf(" %f %f %f %f %f %f\n",slice1->xmin,slice1->xmax,slice1->ymin,slice1->ymax,slice1->zmin,slice1->zmax);
        printf(" %f %f %f %f %f %f\n",slice2->xmin,slice2->xmax,slice2->ymin,slice2->ymax,slice2->zmin,slice2->zmax);
      }
      continue;
    }

    nqframe1 = (is2a+1-is1a)*(js2a+1-js1a)*(ks2a+1-ks1a);
    NewMemory((void **)&qframe1,nqframe1*sizeof(float));
    NewMemory((void **)&qframeout,nqframe1*sizeof(float));
    nqframe2 = (is2b+1-is1b)*(js2b+1-js1b)*(ks2b+1-ks1b);
    NewMemory((void **)&qframe2,nqframe2*sizeof(float));

    len=strlen(outfile);
    unit3=13;
    FORToutsliceheader(outfile,&unit3,&is1a,&is2a,&js1a,&js2a,&ks1a,&ks2a,&error1,len);
    if(error1!=0){
      FORTclosefortranfile(&unit1);
      FORTclosefortranfile(&unit2);
      printf("*** problem writing out header for %s\n",fullfile1);
      continue;
    }
    printf("differencing %s and %s\n",fullfile1,fullfile2);
    for(;;){
      FORTgetsliceframe(&unit1,&is1a,&is2a,&js1a,&js2a,&ks1a,&ks2a,&time1,qframe1,&error1);
      FORTgetsliceframe(&unit2,&is1b,&is2b,&js1b,&js2b,&ks1b,&ks2b,&time2,qframe2,&error2);
      if(error1!=0||error2!=0)break;
      for(i=0;i<nqframe1;i++){
        qframeout[i]=qframe2[i]-qframe1[i];
      }
      FORToutsliceframe(&unit3,&is1a,&is2a,&js1a,&js2a,&ks1a,&ks2a,&time1,qframeout,&error1);
    }


    FORTclosefortranfile(&unit1);
    FORTclosefortranfile(&unit2);
    FORTclosefortranfile(&unit3);
    FREEMEMORY(qframe1);
    FREEMEMORY(qframe2);
    FREEMEMORY(qframeout);
  }
}

/* ------------------ make_fileout ------------------------ */

void make_outfile(char *outfile, char *destdir, char *file1, char *ext){
  char filecopy[1024], *file1_noext;

  strcpy(filecopy,file1);
  file1_noext=strstr(filecopy,ext);
  strcpy(outfile,"");
  if(file1_noext==NULL)return;
  file1_noext[0]='\0';
  if(destdir!=NULL){
    strcpy(outfile,destdir);
  }
  strcat(outfile,filecopy);
  strcat(outfile,"_diff");
  strcat(outfile,ext);
}