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
char CNViso_revision[]="$Revision$";

void getisosizes(EGZ_FILE *isostreamptr,float **levelsptr, int *nisolevels, int *sizebefore);
void endian_switch(void *val, int nval);
int tri_compare( const void *arg1, const void *arg2 );
int tri_compare1( const void *arg1, const void *arg2 );
int tri_compare2( const void *arg1, const void *arg2 );
int tri_compare4( const void *arg1, const void *arg2 );

/* ------------------ convert_iso ------------------------ */

int convert_iso(iso *isoi){

int sortflag;
int *triangles_i;
unsigned char *triangles1_i;
unsigned short *triangles2_i;
int *sortindex=NULL;
int nsort;
int nsortMAX=-1;
unsigned char  *triangle1_copy=NULL;
unsigned short *triangle2_copy=NULL;
int   *triangle_copy=NULL;

  char isofile_svz[1024], isosizefile_svz[1024];
  EGZ_FILE *ISOFILE;
  FILE *isostream, *isosizestream;

  int nvertices_i, ntriangles_i;
  unsigned char *vertexnorm;
  unsigned short *vertices_i;
  int iframe, j, jj, skip;
  float time;
  int sizebefore, sizeafter;
  int nvertices_iMAX=-1, ntriangles_iMAX=-1, ntriangles1_iMAX=-1, ntriangles2_iMAX=-1;
  int nbufferMAX=-1;
  int nbuffer;
  unsigned char *compressed_buffer=NULL, *full_buffer=NULL;
  char xxx[2];
  long data_loc;
  int percent_done;
  int percent_next=10;
  vert *vertinfo;
  int one=1;
  float time_max;
//  float *isolevels;

  strcpy(xxx,"X");

#ifdef EGZ
  ISOFILE=EGZ_FOPEN(isoi->file,"rb",0,2);
#else
  ISOFILE=EGZ_FOPEN(isoi->file,"rb");
#endif

  if(ISOFILE==NULL){
    printf("  %s could not be opened\n",isoi->file);
    return 0;
  }

  // set up iso compressed file

  if(destdir!=NULL){
    strcpy(isofile_svz,destdir);
    strcat(isofile_svz,isoi->filebase);
  }
  else{
    strcpy(isofile_svz,isoi->file);
  }
  strcat(isofile_svz,".svz");

  if(destdir!=NULL){
    strcpy(isosizefile_svz,destdir);
    strcat(isosizefile_svz,isoi->filebase);
  }
  else{
    strcpy(isosizefile_svz,isoi->file);
  }
  strcat(isosizefile_svz,".sz");

  // clean iso zip files

  if(cleanfiles==1){
    isostream=fopen(isofile_svz,"rb");
    if(isostream!=NULL){
      fclose(isostream);
      printf("  Removing %s.\n",isofile_svz);
      UNLINK(isofile_svz);
      LOCK_COMPRESS;
      filesremoved++;
      UNLOCK_COMPRESS;
    }
    isosizestream=fopen(isosizefile_svz,"rb");
    if(isosizestream!=NULL){
      fclose(isosizestream);
      printf("  Removing %s.\n",isosizefile_svz);
      UNLINK(isosizefile_svz);
      LOCK_COMPRESS;
      filesremoved++;
      UNLOCK_COMPRESS;
    }
    return 0;
  }

  if(overwrite_iso==0){
    isostream=fopen(isofile_svz,"rb");
    if(isostream!=NULL){
      fclose(isostream);
      printf("  %s exists.\n",isofile_svz);
      printf("     Use the -f option to overwrite smokezip compressed files\n");
      return 0;
    }
  }

  isostream=fopen(isofile_svz,"wb");
  isosizestream=fopen(isosizefile_svz,"w");
  if(isostream==NULL||isosizestream==NULL){
    if(isostream==NULL){
      printf("  %s could not be opened for writing\n",isofile_svz);
    }
    if(isosizestream==NULL){
      printf("  %s could not be opened for writing\n",isosizefile_svz);
    }
    if(isostream!=NULL)fclose(isostream);
    if(isosizestream!=NULL)fclose(isosizestream);
    EGZ_FCLOSE(ISOFILE);
    return 0;
  }

  // read data
  sizebefore=0;
  sizeafter=0;

  getisosizes(ISOFILE, &isoi->isolevels, &isoi->nisolevels, &sizebefore);
  fprintf(isosizestream,"%i\n",isoi->nisolevels);
//  isolevels=isoi->isolevels;

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

  if(cleanfiles==0)printf("Compressing iso-surface file %s\n",isoi->file);

  fwrite(&one,4,1,isostream);
  fwrite(&isoi->nisolevels,4,1,isostream);
  fwrite(isoi->isolevels,4,isoi->nisolevels,isostream);
  fwrite(&sphereinfo.npoints,4,1,isostream);
  fwrite(sphereinfo.snormals,2,3*sphereinfo.npoints,isostream);

  sizeafter+=(12+4*isoi->nisolevels+6*sphereinfo.npoints);

  triangles_i=NULL;
  triangles1_i=NULL;
  triangles2_i=NULL;
	vertices_i=NULL;
  vertexnorm=NULL;
  vertinfo=NULL;

  jj=-1;
  iframe=0;
  time_max=-1000000.0;
  for(;;){
    int memory_fail=0;

    jj++;
    EGZ_FREAD(&time,4,1,ISOFILE);
    if(EGZ_FEOF(ISOFILE)!=0)goto wrapup;
    sizebefore+=4;

    if(jj%isozipstep==0){
      fprintf(isosizestream," %f\n",time);
      fwrite(&time,4,1,isostream);
    }
    sizeafter+=4;
    for(j=0;j<isoi->nisolevels;j++){
      EGZ_FREAD(&nvertices_i,4,1,ISOFILE);
	    if(EGZ_FEOF(ISOFILE)!=0)goto wrapup;
      sizebefore+=4;
	    EGZ_FREAD(&ntriangles_i,4,1,ISOFILE);
	    if(EGZ_FEOF(ISOFILE)!=0)goto wrapup;

      sizebefore+=4;
      if(jj%isozipstep!=0||time<time_max){
        if(jj%isozipstep==0&&time<time_max)jj--;
        skip=0;
        if(nvertices_i<=0||ntriangles_i<=0)continue;
        skip += (6*nvertices_i);
        if(isoi->dataflag==1)skip += (8 + 2*nvertices_i);
  	    if(nvertices_i<256){
  	      skip += (ntriangles_i);
        }
        else if(nvertices_i>=256&&nvertices_i<65536){
  	      skip += (ntriangles_i*2);
        }
	      else{
  	      skip += (ntriangles_i*4);
        }
        EGZ_FSEEK(ISOFILE,skip,SEEK_CUR);
        continue;
      }
      time_max=time;

#define FREELOCAL_ISO FREEMEMORY(vertexnorm);FREEMEMORY(vertices_i);FREEMEMORY(vertinfo);\
                      FREEMEMORY(triangles1_i);FREEMEMORY(triangles2_i);FREEMEMORY(triangles_i);;

      // allocate memory for and read in vertices

      if(nvertices_i>0){
        if(nvertices_i>nvertices_iMAX){
          FREEMEMORY(vertexnorm);
          FREEMEMORY(vertices_i);
          if(NewMemory((void **)&vertexnorm,nvertices_i*sizeof(unsigned char))==0){
            FREELOCAL_ISO;
            break;
          }
          if(NewMemory((void **)&vertices_i,3*nvertices_i*sizeof(unsigned short))==0){
            FREELOCAL_ISO;
            memory_fail=1;
            break;
          }
          if(NewMemory((void **)&vertinfo,nvertices_i*sizeof(vert))==0){
            FREELOCAL_ISO;
            memory_fail=1;
            break;
          }
          nvertices_iMAX=nvertices_i;
        }
        EGZ_FREAD(vertices_i,2,3*nvertices_i,ISOFILE);
  	    if(EGZ_FEOF(ISOFILE)!=0)goto wrapup;
        sizebefore+=2*3*nvertices_i;
//        sizebefore+=2*3*nvertices_i;  //  not in .iso file but needed in smokeview for uncompressed isosurrface display
//        sizebefore+=6*ntriangles_i;   // not in iso file but also needed in smokeview
      }

      // read in triangle data

      // case 1: less than 256 vertices so need only unsigned char's

      if(nvertices_i<256&&nvertices_i>0){
        if(ntriangles_i>0){
          if(ntriangles_i>ntriangles1_iMAX){
            FREEMEMORY(triangles1_i);
            if(
              NewMemory((void **)&triangles1_i,ntriangles_i*sizeof(unsigned char))==0){
              FREELOCAL_ISO;
              memory_fail=1;
              break;
            }
            ntriangles1_iMAX=ntriangles_i;
          }
  	    	EGZ_FREAD(triangles1_i,1,ntriangles_i,ISOFILE);
    	    if(EGZ_FEOF(ISOFILE)!=0)goto wrapup;
          sizebefore+=ntriangles_i;
        }
      }

      // case 2: 256<=nvertices_i<65536 so need unsigned short's

      else if(nvertices_i>=256&&nvertices_i<65536){
        if(ntriangles_i>0){
          if(ntriangles_i>ntriangles2_iMAX){
            FREEMEMORY(triangles2_i);
            if(
              NewMemory((void **)&triangles2_i,ntriangles_i*sizeof(unsigned short))==0){
              FREELOCAL_ISO;
              memory_fail=1;
              break;
            }
            ntriangles2_iMAX=ntriangles_i;
          }
  		    EGZ_FREAD(triangles2_i,2,ntriangles_i,ISOFILE);
     	    if(EGZ_FEOF(ISOFILE)!=0)goto wrapup;
          sizebefore+=2*ntriangles_i;
        }
      }

      // case 3: more than 65535 vertices so need int's

      else{
        if(ntriangles_i>0){
          if(ntriangles_i>ntriangles_iMAX){
            FREEMEMORY(triangles_i);
            if(NewMemory((void **)&triangles_i,ntriangles_i*sizeof(int))==0){
              FREELOCAL_ISO;
              memory_fail=1;
              break;
            }
            ntriangles_iMAX=ntriangles_i;
          }
  		    EGZ_FREAD(triangles_i,4,ntriangles_i,ISOFILE);
    	    if(EGZ_FEOF(ISOFILE)!=0)goto wrapup;
          sizebefore+=(4*ntriangles_i);
        }
      }

      // COMPRESS DATA
      // data:
      // time, nvertices_i, ntriangles_i, vertices_i, triangles_i
      //
      fwrite(&nvertices_i,4,1,isostream);
      fwrite(&ntriangles_i,4,1,isostream);
      sizeafter+=8;

      nbuffer=0;
      if(ntriangles_i<=0||nvertices_i<=0){
        int nfull=0;

        fwrite(&nfull,4,1,isostream);
        fwrite(&nbuffer,4,1,isostream);
        fprintf(isosizestream," %i %i %i %i\n",nvertices_i,ntriangles_i,0,0);
      }
      if(ntriangles_i>0&&nvertices_i>0){
        int i;
#ifdef WIN32
        uLongf ncompressed;
//        int ncompressed;
#else
        int ncompressed;
#endif
        int nfull;

        // allocate buffer needed for compression

        nbuffer=7*nvertices_i;
        if(ntriangles_i<256){
          nbuffer+=ntriangles_i;
        }
        else if(ntriangles_i>=256&&ntriangles_i<65536){
          nbuffer+=2*ntriangles_i;
        }
        else if(ntriangles_i>=65536){
          nbuffer+=4*ntriangles_i;
        }
        nfull=nbuffer;
        nbuffer = 1.01*nbuffer+600;
        ncompressed=nbuffer;
        if(nbuffer>nbufferMAX){
          FREEMEMORY(compressed_buffer);
          FREEMEMORY(full_buffer);
          NewMemory((void **)&full_buffer,nbuffer);
          NewMemory((void **)&compressed_buffer,nbuffer);
          nbufferMAX=nbuffer;
        }

        // compute normals

        for(i=0;i<nvertices_i;i++){
          vert *verti;
          float *xyz;

          verti = vertinfo + i;
          xyz = verti->normal;
          xyz[0]=0.0;
          xyz[1]=0.0;
          xyz[2]=0.0;
        }

        for(i=0;i<ntriangles_i/3;i++){
          unsigned short *v1, *v2, *v3;
          float normal[3];
          float area;
          int ind[3];
          int jjj;

          if(nvertices_i<256){
            ind[0]=triangles1_i[3*i];
            ind[1]=triangles1_i[3*i+1];
            ind[2]=triangles1_i[3*i+2];
          }
          else if(nvertices_i>=256&&nvertices_i<65536){
            ind[0]=triangles2_i[3*i];
            ind[1]=triangles2_i[3*i+1];
            ind[2]=triangles2_i[3*i+2];
          }
          else{
            ind[0]=triangles_i[3*i];
            ind[1]=triangles_i[3*i+1];
            ind[2]=triangles_i[3*i+2];
          }
          v1 = vertices_i + 3*ind[0];
          v2 = vertices_i + 3*ind[1];
          v3 = vertices_i + 3*ind[2];
          Normal(v1,v2,v3,normal,&area);
          for(jjj=0;jjj<3;jjj++){
            vert *vertj;

            vertj = vertinfo + ind[jjj];
            vertj->normal[0]+=area*normal[0];
            vertj->normal[1]+=area*normal[1];
            vertj->normal[2]+=area*normal[2];
          }
        }
        for(i=0;i<nvertices_i;i++){
          vert *verti;

          verti = vertinfo + i;
          vertexnorm[i]=getnormalindex(&sphereinfo,verti->normal);
        }

// void Normal(unsigned short *v1, unsigned short *v2, unsigned short *v3, float *normal){


        // reduce precision of vertices to improve compression - but don't reduce too much

        for(i=0;i<3*nvertices_i;i++){
          vertices_i[i]/=32;
          vertices_i[i]*=32;
        }
        
        // reorder vertices putting vertex with smallest index first
        //   but keep cyclical order in tact

        for(i=0;i<ntriangles_i/3;i++){
          unsigned char *tri1,dumm1[3];
          unsigned short *tri2,dumm2[3];
          int *tri4,dumm4[3];

          if(nvertices_i<256){
            tri1 = triangles1_i + 3*i;
            if(tri1[0]<=tri1[1]&&tri1[0]<=tri1[2])continue;
            if(tri1[1]<=tri1[0]&&tri1[1]<=tri1[2]){
              dumm1[0]=tri1[1];
              dumm1[1]=tri1[2];
              dumm1[2]=tri1[0];
              tri1[0]=dumm1[0];
              tri1[1]=dumm1[1];
              tri1[2]=dumm1[2];
              continue;
            }
            dumm1[0]=tri1[2];
            dumm1[1]=tri1[0];
            dumm1[2]=tri1[1];
            tri1[0]=dumm1[0];
            tri1[1]=dumm1[1];
            tri1[2]=dumm1[2];
            continue;
          }
          else if(nvertices_i>=256&&nvertices_i<65536){
            tri2 = triangles2_i + 3*i;
            if(tri2[0]<=tri2[1]&&tri2[0]<=tri2[2])continue;
            if(tri2[1]<=tri2[0]&&tri2[1]<=tri2[2]){
              dumm2[0]=tri2[1];
              dumm2[1]=tri2[2];
              dumm2[2]=tri2[0];
              tri2[0]=dumm2[0];
              tri2[1]=dumm2[1];
              tri2[2]=dumm2[2];
              continue;
            }
            dumm2[0]=tri2[2];
            dumm2[1]=tri2[0];
            dumm2[2]=tri2[1];
            tri2[0]=dumm2[0];
            tri2[1]=dumm2[1];
            tri2[2]=dumm2[2];
            continue;
          }
          else if(nvertices_i>=65536){
            tri4 = triangles_i + 3*i;
            if(tri4[0]<=tri4[1]&&tri4[0]<=tri4[2])continue;
            if(tri4[1]<=tri4[0]&&tri4[1]<=tri4[2]){
              dumm4[0]=tri4[1];
              dumm4[1]=tri4[2];
              dumm4[2]=tri4[0];
              tri4[0]=dumm4[0];
              tri4[1]=dumm4[1];
              tri4[2]=dumm4[2];
              continue;
            }
            dumm4[0]=tri4[2];
            dumm4[1]=tri4[0];
            dumm4[2]=tri4[1];
            tri4[0]=dumm4[0];
            tri4[1]=dumm4[1];
            tri4[2]=dumm4[2];
            continue;
          }
        }

        {
          nsort=ntriangles_i/3;
          if(nsort>nsortMAX){
            FREEMEMORY(sortindex);
            NewMemory((void **)&sortindex,nsort*sizeof(int));
            FREEMEMORY(triangle1_copy);
            NewMemory((void **)&triangle1_copy,ntriangles_i);
            FREEMEMORY(triangle2_copy);
            NewMemory((void **)&triangle2_copy,ntriangles_i*sizeof(unsigned short));
            FREEMEMORY(triangle_copy);
            NewMemory((void **)&triangle_copy,ntriangles_i*sizeof(int));
            nsortMAX=nsort;
          }
          for(i=0;i<nsort;i++){
            sortindex[i]=i;
          }
          if(nvertices_i<256){
            sortflag=1;
          }
          else if(nvertices_i>=256&&nvertices_i<65536){
            sortflag=2;
          }
          else {
            sortflag=4;
          }
          switch (sortflag){
          case 1:
            qsort((unsigned char *)triangles1_i,(size_t)ntriangles_i/3,3*sizeof(unsigned char),tri_compare1);
            for(i=0;i<ntriangles_i/3;i++){
              triangle1_copy[3*i+0]=triangles1_i[3*i+0];
              triangle1_copy[3*i+1]=triangles1_i[3*i+1];
              triangle1_copy[3*i+2]=triangles1_i[3*i+2];
            }
            break;
          case 2:
            qsort((unsigned short *)triangles2_i,(size_t)ntriangles_i/3,3*sizeof(unsigned short),tri_compare2);
            for(i=0;i<ntriangles_i/3;i++){
              triangle2_copy[3*i+0]=triangles2_i[3*i+0];
              triangle2_copy[3*i+1]=triangles2_i[3*i+1];
              triangle2_copy[3*i+2]=triangles2_i[3*i+2];
            }
            break;
          case 4:
            qsort((int *)triangles_i,(size_t)ntriangles_i/3,3*sizeof(int),tri_compare4);
            for(i=0;i<ntriangles_i/3;i++){
              triangle_copy[3*i+0]=triangles_i[3*i+0];
              triangle_copy[3*i+1]=triangles_i[3*i+1];
              triangle_copy[3*i+2]=triangles_i[3*i+2];
            }
            break;
          }
        }

        // copy filtered, sorted data into buffer to be compressed

        memcpy(full_buffer,vertices_i,6*nvertices_i);
        memcpy(full_buffer+6*nvertices_i,vertexnorm,nvertices_i);
        if(nvertices_i<256){
          memcpy(full_buffer+7*nvertices_i,triangle1_copy,ntriangles_i);
        }
        else if(nvertices_i>=256&&nvertices_i<65536){
          memcpy(full_buffer+7*nvertices_i,triangle2_copy,2*ntriangles_i);
        }
        else if(nvertices_i>=65536){
          memcpy(full_buffer+7*nvertices_i,triangle_copy,4*ntriangles_i);
        }

        // compress data
        compress(compressed_buffer, &ncompressed, full_buffer, nfull);

        // write out compressed data

        fwrite(&nfull,4,1,isostream);
        fwrite(&ncompressed,4,1,isostream);
        fwrite(compressed_buffer,1,ncompressed,isostream);
        sizeafter+=(8+ncompressed);
        fprintf(isosizestream," %i %i %i %i\n",nvertices_i,ntriangles_i,nfull,(int)ncompressed);


        data_loc=EGZ_FTELL(ISOFILE);
        percent_done=100.0*(float)data_loc/(float)isoi->filesize;
        if(percent_done>percent_next){
          printf(" %i%s",percent_next,pp);
          LOCK_COMPRESS;
          fflush(stdout);
          UNLOCK_COMPRESS;
          percent_next+=10;
        }
      }
    }
    LOCK_ISOS;
    if(memory_fail==1){
      printf("memory allocation error\n");
      UNLOCK_ISOS;
      break;
    }
    iframe++;
    UNLOCK_ISOS;
  }
wrapup:
  FREELOCAL_ISO;
  printf(" 100%s completed\n",pp);
  EGZ_FCLOSE(ISOFILE);
  fclose(isostream);
  fclose(isosizestream);
  {
    char before_label[256],after_label[256];
    getfilesizelabel(sizebefore,before_label);
    getfilesizelabel(sizeafter,after_label);
//    printf("    records=%i, ",count);
    printf("Sizes: original=%s, ",before_label);
    printf("compressed=%s (%4.1f%s reduction)\n",after_label,(float)sizebefore/(float)sizeafter,xxx);
  }
  return 0;
}


//  int *triangles_i;
//  unsigned char *triangles1_i;
//  short *triangles2_i;

/* ------------------ compare ------------------------ */
/*
int tri_compare( const void *arg1, const void *arg2 ){
  int i, j;
  unsigned char x1[3], y1[3];
  unsigned short x2[3], y2[3];
  int x4[3], y4[3];

  i=*(int *)arg1;
  j=*(int *)arg2;
  switch (sortflag){
  case 1:
    x1[0] = triangles1_i[3*i];
    y1[0] = triangles1_i[3*j];
    if(x1[0]<y1[0])return -1;
    if(x1[0]>y1[0])return 1;

    x1[1] = triangles1_i[3*i+1];
    y1[1] = triangles1_i[3*j+1];
    if(x1[1]<y1[1])return -1;
    if(x1[1]>y1[1])return 1;

    x1[2] = triangles1_i[3*i+2];
    y1[2] = triangles1_i[3*j+2];
    if(x1[2]<y1[2])return -1;
    if(x1[2]>y1[2])return 1;

    return 0;
    break;
  case 2:
    x2[0] = triangles2_i[3*i];
    y2[0] = triangles2_i[3*j];
    if(x2[0]<y2[0])return -1;
    if(x2[0]>y2[0])return 1;

    x2[1] = triangles2_i[3*i+1];
    y2[1] = triangles2_i[3*j+1];
    if(x2[1]<y2[1])return -1;
    if(x2[1]>y2[1])return 1;

    x2[2] = triangles2_i[3*i+2];
    y2[2] = triangles2_i[3*j+2];
    if(x2[2]<y2[2])return -1;
    if(x2[2]>y2[2])return 1;

    return 0;
    break;
  case 4:
    x4[0] = triangles_i[3*i];
    y4[0] = triangles_i[3*j];
    if(x4[0]<y4[0])return -1;
    if(x4[0]>y4[0])return 1;

    x4[1] = triangles_i[3*i+1];
    y4[1] = triangles_i[3*j+1];
    if(x4[1]<y4[1])return -1;
    if(x4[1]>y4[1])return 1;

    x4[2] = triangles_i[3*i+2];
    y4[2] = triangles_i[3*j+2];
    if(x4[2]<y4[2])return -1;
    if(x4[2]>y4[2])return 1;

    return 0;
    break;
  }
  return 0;
}
*/
/* ------------------ tri_compare4 ------------------------ */

int tri_compare4( const void *arg1, const void *arg2 ){
  int *x1, *y1;

  x1=(int *)arg1;
  y1=(int *)arg2;

  if(x1[0]<y1[0])return -1;
  if(x1[0]>y1[0])return 1;

  if(x1[1]<y1[1])return -1;
  if(x1[1]>y1[1])return 1;

  if(x1[2]<y1[2])return -1;
  if(x1[2]>y1[2])return 1;

  return 0;
}

/* ------------------ tri_compare2 ------------------------ */

int tri_compare2( const void *arg1, const void *arg2 ){
  unsigned short *x1, *y1;

  x1=(unsigned short *)arg1;
  y1=(unsigned short *)arg2;

  if(x1[0]<y1[0])return -1;
  if(x1[0]>y1[0])return 1;

  if(x1[1]<y1[1])return -1;
  if(x1[1]>y1[1])return 1;

  if(x1[2]<y1[2])return -1;
  if(x1[2]>y1[2])return 1;

  return 0;
}

/* ------------------ tri_compare1 ------------------------ */

int tri_compare1( const void *arg1, const void *arg2 ){
  unsigned char *x1, *y1;

  x1=(unsigned char *)arg1;
  y1=(unsigned char *)arg2;

  if(x1[0]<y1[0])return -1;
  if(x1[0]>y1[0])return 1;

  if(x1[1]<y1[1])return -1;
  if(x1[1]>y1[1])return 1;

  if(x1[2]<y1[2])return -1;
  if(x1[2]>y1[2])return 1;

  return 0;
}

/* ------------------ convert_iso ------------------------ */

void *compress_isos(void *arg){
  int i;
  iso *isoi;

  if(niso_files>0){
    LOCK_ISOS;
    if(first_initsphere==1){
      first_initsphere=0;
      initspherepoints(&sphereinfo,14);
    }
    UNLOCK_ISOS;
  // convert and compress files
    printf("\n");
    for(i=0;i<niso_files;i++){
      isoi = isoinfo + i;
      if(autozip==1&&isoi->autozip==0)continue;
      LOCK_ISOS;
      if(isoi->inuse==1){
        UNLOCK_ISOS;
        continue;
      }
      isoi->inuse=1;
      UNLOCK_ISOS;
      convert_iso(isoi);
    }
  }
  return NULL;
}

/* ------------------ getisosizes ------------------------ */

void getisosizes(EGZ_FILE *ISOFILE,float **levelsptr, int *nisolevels, int *sizebefore){
  int len[3],labellengths=0;
  int nlevels;
  int one;

  *sizebefore=0;
  EGZ_FREAD(&one,4,1,ISOFILE);
  *sizebefore+=4;
  EGZ_FREAD(len,4,3,ISOFILE);
  *sizebefore+=12;
  labellengths=len[0]+len[1]+len[2];
  EGZ_FSEEK(ISOFILE,labellengths+4,SEEK_CUR);
  *sizebefore+=(labellengths+4);;
  EGZ_FREAD(&nlevels,4,1,ISOFILE);
  *sizebefore+=4;
  *nisolevels=nlevels;
  if(*levelsptr==NULL){
    if(NewMemory((void **)levelsptr,*nisolevels*sizeof(float))==0)return;
  }
  else{
    if(ResizeMemory((void **)levelsptr,*nisolevels*sizeof(float))==0)return;
  }
  EGZ_FREAD(*levelsptr,4,nlevels,ISOFILE);
  *sizebefore+=(4*nlevels);
}
