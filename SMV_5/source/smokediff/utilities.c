// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "svdiff.h"
#include "svn_revision.h"
#include "MALLOC.h"

// svn revision character string
char utilities_revision[]="$Revision$";

/* ------------------ match ------------------------ */

int match(const char *buffer, const char *key, unsigned int lenkey){
  if(strncmp(buffer,key,lenkey) == 0)return(1);
  return(0);
}


/* ------------------ match ------------------------ */

int mesh_match(mesh *mesh1, mesh *mesh2){
  int ibar, jbar, kbar;

  if(mesh1->ibar!=mesh2->ibar)return 0;
  if(mesh1->jbar!=mesh2->jbar)return 0;
  if(mesh1->kbar!=mesh2->kbar)return 0;
  ibar=mesh1->ibar;
  jbar=mesh1->jbar;
  kbar=mesh1->kbar;
  if(fabs(mesh1->xplt[0]-mesh2->xplt[0])>mesh1->dx)return 0;
  if(fabs(mesh1->yplt[0]-mesh2->yplt[0])>mesh1->dy)return 0;
  if(fabs(mesh1->zplt[0]-mesh2->zplt[0])>mesh1->dz)return 0;
  if(fabs(mesh1->xplt[ibar]-mesh2->xplt[ibar])>mesh1->dx)return 0;
  if(fabs(mesh1->yplt[jbar]-mesh2->yplt[jbar])>mesh1->dy)return 0;
  if(fabs(mesh1->zplt[kbar]-mesh2->zplt[kbar])>mesh1->dz)return 0;
  return 1;
}

/* ------------------ trim_front ------------------------ */

char *trim_front(char *line){
  char *blank=" ";
  const char *c;
  unsigned int i,len;

  c = line;
  len=strlen(line);
  for(i=0;i<len;i++){
    if(strncmp(c++,blank,1)!=0)return line+i;
  }
  return line;
}

/* ------------------ trim ------------------------ */

void trim(char *line){
  char *blank=" ";
  const char *c;
  const char *lf="\n";
  unsigned int len;
  unsigned int i;
  len = strlen(line);
  c = line+len-1;
  for(i=0; i<len; i++){
    if(strncmp(c,blank,1)!=0&&strncmp(c,lf,1)!=0){
      c++; 
      line[c-line]='\0';
      return;
    }
    c--;
  }
  *line='\0';
}

/* ------------------ readlabels ------------------------ */

int readlabels(flowlabels *flowlabel, FILE *stream){
  char buffer[255];
  unsigned int len;

  flowlabel->longlabel=NULL;
  flowlabel->shortlabel=NULL;
  flowlabel->unit=NULL;

  if(fgets(buffer,255,stream)==NULL)return 2;

  len=strlen(buffer);
  buffer[len-1]='\0';
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void **)&flowlabel->longlabel,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->longlabel,buffer);


  if(fgets(buffer,255,stream)==NULL)return 2;
  len=strlen(buffer);
  buffer[len-1]='\0';
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void **)&flowlabel->shortlabel,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->shortlabel,buffer);

  if(fgets(buffer,255,stream)==NULL)return 2;
  len=strlen(buffer);
  buffer[len-1]='\0';
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void *)&flowlabel->unit,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->unit,buffer);
  return 0;
}

  /* ------------------ getfileinfo ------------------------ */

int getfileinfo(char *filename, char *source_dir, int *filesize){
  STRUCTSTAT statbuffer;
  int statfile;
  char buffer[1024];

  if(source_dir==NULL){
    strcpy(buffer,filename);
  }
  else{
    strcpy(buffer,source_dir);
    strcat(buffer,filename);
  }
  if(filesize!=NULL)*filesize=0;
  statfile=STAT(buffer,&statbuffer);
  if(statfile!=0)return statfile;
  if(filesize!=NULL)*filesize=statbuffer.st_size;
  return statfile;
}

/* ------------------ version ------------------------ */

void version(void){
    char smv_version[100];
    int svn_num;

    getSMDiffversion(smv_version);  // get Smokeview version (ie 5.x.z)
    svn_num=getmaxrevision();    // get svn revision number
    printf("\n");
    printf("Smokediff\n\n");
    printf("Version: %s\n",smv_version);
    printf("SVN Revision Number: %i\n",svn_num);
    printf("Compile Date: %s\n",__DATE__);
#ifdef X64
    printf("Platform: WIN64\n");
#endif
#ifdef WIN32
#ifndef X64
    printf("Platform: WIN32\n");
#endif
#endif
#ifdef pp_OSX
    printf("Platform: OS X\n");
#endif
#ifdef pp_LINUX
    printf("Platform: LINUX\n");
#endif

}

/* ------------------ getmaxrev ------------------------ */

#define MAXREV(cval) max_revision=imax(getrevision(cval),max_revision)
int getmaxrevision(void){
  int max_revision=0;

  MAXREV(main_revision);
  MAXREV(utilities_revision);
  MAXREV(readsmv_revision);
  MAXREV(dmalloc_revision);
  MAXREV(assert_revision);
  MAXREV(IOdslice_revision);
  MAXREV(IOdboundary_revision);
  MAXREV(IOdplot_revision);
  return max_revision;
}

/* ------------------ getSMVversion ------------------------ */

void getSMDiffversion(char *SMDiffversion){
  strcpy(SMDiffversion,SMDiffVERSION);
}

/* ------------------ imax ------------------------ */

int imax(int a, int b){
  if(a>b){
    return a;
  }
  else{
    return b;
  }
}

/* ------------------ getrevision ------------------------ */

int getrevision(char *svn){
  char svn_string[256];
  char *svn_ptr;
  int return_val;

  svn_ptr=svn_string;
  svn=strchr(svn,':');
  if(svn==NULL||strlen(svn)<=4)return 0;
  
  svn++;
  strcpy(svn_ptr,svn);
  svn_ptr=trim_front(svn_ptr);
  svn_ptr[strlen(svn_ptr)-1]=0;
  trim(svn_ptr);
  sscanf(svn_ptr,"%i",&return_val);
  return return_val;
}

/* ------------------ fullfile ------------------------ */

void fullfile(char *fileout, char *dir, char *file){
  strcpy(fileout,"");
  if(dir!=NULL)strcat(fileout,dir);
  strcat(fileout,file);
}

/* ------------------ getendian ------------------------ */

int getendian(void){
  short val;
  char *cval;
  val=1;
  cval = (char *)&val+1;
  return (int)(*cval);
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

/* ------------------ get_hist_val ------------------------ */

float get_hist_val(bucketdata *bucket, float cdf){
  int cutoff, count;
  int i;
  float returnval;

  cutoff = cdf*bucket->ntotal;
  count=0;
  for(i=0;i<NBUCKETS;i++){
    count+=bucket->buckets[i];
    if(count>cutoff){
      returnval = bucket->valmin + (float)(i+0.5)*(bucket->valmax-bucket->valmin)/(float)NBUCKETS;
      return returnval;
    }
  }
  return bucket->valmax;
}

/* ------------------ init_buckets ------------------------ */

void init_buckets(bucketdata *bucket){
  int i;

  for(i=0;i<NBUCKETS;i++){
    bucket->buckets[i]=0;
  }
  bucket->ntotal=0;
  bucket->minmax_defined=0;
}

/* ------------------ update_buckets ------------------------ */

void update_buckets(float *vals, int nvals, bucketdata *bucket){
  int i;
  float bucket_copy[NBUCKETS];
  float valmin_vals, valmax_vals;
  float valmin_new, valmax_new;
  int remap;

  // find min and max of vals array

  valmin_vals = vals[0];
  valmax_vals = valmin_vals;
  for(i=1;i<nvals;i++){
    if(vals[i]<valmin_vals)valmin_vals=vals[i];
    if(vals[i]>valmax_vals)valmax_vals=vals[i];
  }

  remap=0;
  if(bucket->minmax_defined==0){
    bucket->minmax_defined=1;
    valmin_new=valmin_vals;
    valmax_new=valmax_vals;
  }
  else{

    // update valmin and valmax

    valmin_new=bucket->valmin;
    valmax_new=bucket->valmax;
    if(valmin_vals<valmin_new){
      valmin_new=valmin_vals;
      remap=1;
    }
    if(valmax_vals>valmax_new){
      valmax_new=valmin_vals;
      remap=1;
    }
  }
  if(remap==1){

    // map old bucket array using new min/max bounds

    for(i=0;i<NBUCKETS;i++){
      bucket_copy[i]=0;
    }
    if(valmin_new==valmax_new){
      for(i=0;i<NBUCKETS;i++){
        if(bucket->buckets[i]==0)continue;
        bucket_copy[0]+=bucket->buckets[i];
      }
    }
    else{
      float dbucket;

      dbucket = (bucket->valmax-bucket->valmin)/(float)NBUCKETS;
      for(i=0;i<NBUCKETS;i++){
        float val;
        int ival;

        if(bucket->buckets[i]==0)continue;
        val = bucket->valmin + (float)(i+0.5)*dbucket;
        ival = NBUCKETS*(val-valmin_new)/(valmax_new-valmin_new);
        if(ival<0)ival=0;
        if(ival>NBUCKETS-1)ival=NBUCKETS-1;
        bucket_copy[ival]+=bucket->buckets[i];
      }
    }
    for(i=0;i<NBUCKETS;i++){
      bucket->buckets[i]=bucket_copy[i];
    }
  }

  // update bucket with vals passed in

  if(valmax_new==valmin_new){
    bucket->buckets[0]+=nvals;
  }
  else{
    for(i=0;i<nvals;i++){
      int ival;

      ival = NBUCKETS*(vals[i]-valmin_new)/(valmax_new-valmin_new);
      if(ival<0)ival=0;
      if(ival>NBUCKETS-1)ival=NBUCKETS-1;
      bucket->buckets[ival]++;
    }
  }
  bucket->valmin=valmin_new;
  bucket->valmax=valmax_new;
  bucket->ntotal+=nvals;
}


