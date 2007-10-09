// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "zlib.h"
#include "egz_stdio.h"
#include "svzip.h"
#include "MALLOC.h"
#include "svn_revision.h"

// svn revision character string
char utilities_revision[]="$Revision$";

#define MARK 255

#define FORTREAD(read) fseek(BOUNDARYFILE,4,SEEK_CUR);returncode=read;fseek(BOUNDARYFILE,4,SEEK_CUR);

/* ------------------ match ------------------------ */

int match(const char *buffer, const char *key, unsigned int lenkey){
  if(strncmp(buffer,key,lenkey) == 0)return(1);
  return(0);
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

/* ------------------ irle ------------------------ */

unsigned int irle(unsigned char *buffer_in, int nchars_in, unsigned char *buffer_out){
  int nrepeats,nn;
  unsigned char thischar, *buffer_in_end;

  nn=0;
  buffer_in_end  = buffer_in  + nchars_in;

  while(buffer_in<buffer_in_end){
    if(*buffer_in==MARK){
      if(buffer_in+2>=buffer_in_end)break;
      buffer_in++;
      thischar=*buffer_in++;
      nrepeats=*buffer_in++;
      nn+=nrepeats;
      memset(buffer_out,thischar,nrepeats);
      buffer_out+=nrepeats;
    }
    else{
      *buffer_out++=*buffer_in++;
      nn++;
    }


  }
  return nn;
}

/* ------------------ readlabels ------------------------ */

int readlabels(flowlabels *flowlabel, FILE *stream){
  char buffer[255];
  unsigned int len;

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

  /* ------------------ smoothlabel ------------------------ */

void smoothlabel(float *a, float *b, int n){
  float delta, factor, logdelta;
  int ndigits;

  delta = (*b-*a)/(n-2);
  if(delta==0.0)return;
  logdelta = log10((double)delta);
  ndigits=logdelta-1;
  if(logdelta<=1)ndigits--;
  factor = 5*pow(10,ndigits);
  delta = (int)(delta/factor + 0.5f)*factor;

  *a = factor*(int)(*a/factor+0.5f);
  *b = *a + (n-2)*delta;

}

  /* ------------------ trimzeros ------------------------ */

void trimzeros(char *line){
  unsigned int i;
  unsigned int len;
  char *c;

  len = strlen(line);
  c = line + len-1;
  for(i=len-1;i>0;i--){
    if(*c=='0'){
      c--;
      if(*c=='.'){
        line[i+1]='\0';
        return;
      }
      continue;
    }
    line[i+1]='\0';
    return;
  }
  line[0]='\0';
}

  /* ------------------ getfileinfo ------------------------ */

int getfileinfo(char *filename, char *source_dir, int *filesize){
  struct stat statbuffer;
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
  statfile=stat(buffer,&statbuffer);
  if(statfile!=0)return statfile;
  if(filesize!=NULL)*filesize=statbuffer.st_size;
  return statfile;
}

  /* ------------------ getfilesizelabel ------------------------ */

void getfilesizelabel(int size, char *sizelabel){
  int leftsize,rightsize;

#define sizeGB   1000000000
#define size100MB 100000000
#define size10MB   10000000
#define sizeMB     1000000
#define size100KB    100000
#define size10KB      10000

  if(size>=sizeGB){
    size/=size10MB;
    leftsize=size/100;
    rightsize=size-100*leftsize;
    sprintf(sizelabel,"%i.%02i GB",leftsize,rightsize);
  }
  else if(size>=size100MB&&size<sizeGB){
    size/=sizeMB;
    leftsize=size;
    sprintf(sizelabel,"%i MB",leftsize);
  }
  else if(size>=size10MB&&size<size100MB){
    size/=size100KB;
    leftsize=size/10;
    rightsize=size-10*leftsize;
    sprintf(sizelabel,"%i.%i MB",leftsize,rightsize);
  }
  else if(size>=sizeMB&&size<size10MB){
    size/=size10KB;
    leftsize=size/100;
    rightsize=size-100*leftsize;
    sprintf(sizelabel,"%i.%02i MB",leftsize,rightsize);
  }
  else if(size>=size100KB&&size<sizeMB){
    size/=1000;
    leftsize=size;
    sprintf(sizelabel,"%i KB",leftsize);
  }
  else{
    size/=10;
    leftsize=size/100;
    rightsize=size-100*leftsize;
    sprintf(sizelabel,"%i.%02i KB",leftsize,rightsize);
  }
}

/* ------------------ calcNormal2 ------------------------ */

void Normal(unsigned short *v1, unsigned short *v2, unsigned short *v3, float *normal, float *area){
  float u[3], v[3];

  float norm2;

  u[0]=v2[0]-v1[0];
  u[1]=v2[1]-v1[1];
  u[2]=v2[2]-v1[2];

  v[0]=v3[0]-v1[0];
  v[1]=v3[1]-v1[1];
  v[2]=v3[2]-v1[2];

  normal[0] = u[1]*v[2] - u[2]*v[1];
  normal[1] = u[2]*v[0] - u[0]*v[2];
  normal[2] = u[0]*v[1] - u[1]*v[0];

  norm2 = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  *area = norm2/2.0;

  normal[0]/=norm2;
  normal[1]/=norm2;
  normal[2]/=norm2;


}

/* ------------------ version ------------------------ */

void version(void){
    char smv_version[100];
    int svn_num;

    getSMZversion(smv_version);  // get Smokeview version (ie 5.x.z)
    svn_num=getmaxrevision();    // get svn revision number
    printf("\n");
    printf("Smokezip\n\n");
    printf("Version: %s\n",smv_version);
    printf("SVN Revision Number: %i\n",svn_num);
    printf("Compile Date: %s\n",__DATE__);
#ifdef WIN32
    printf("Platform: WIN32\n");
#endif
#ifdef pp_OSX
    printf("Platform: OS X\n");
#endif
#ifdef pp_LINUX
    printf("Platform: LINUX\n");
#endif

}

int imax(int a, int b){
  if(a>b){
    return a;
  }
  else{
    return b;
  }
}

#define MAXREV(cval) max_revision=imax(getrevision(cval),max_revision)
int getmaxrevision(void){
  int max_revision=0;

  MAXREV(assert_revision);
  MAXREV(CNV3dsmoke_revision);
  MAXREV(CNVboundary_revision);
  MAXREV(CNViso_revision);
  MAXREV(CNVpart_revision);
  MAXREV(CNVslice_revision);
  MAXREV(csphere_revision);
  MAXREV(dmalloc_revision);
  MAXREV(egz_stdio_revision);
  MAXREV(endian_revision);
  MAXREV(main_revision);
  MAXREV(readfiles_revision);
  MAXREV(stats_revision);
  MAXREV(utilities_revision);
  return max_revision;
}


/* ------------------ getSMVversion ------------------------ */

void getSMZversion(char *SMZversion){
  strcpy(SMZversion,"1.1.1");
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
