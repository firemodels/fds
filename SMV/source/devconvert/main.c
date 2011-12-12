// $Date: 2011-03-24 16:22:30 -0400 (Thu, 24 Mar 2011) $ 
// $Revision: 7970 $
// $Author: gforney $

#define INMAIN
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "svn_revision.h"
#include "datadefs.h"
#include "string_util.h"
#include "MALLOC.h"

//dummy change to bump version number to 0.9

// svn revision character string
char main_revision[]="$Revision: 7970 $";

void usage(char *prog);
void version(char *prog);
void getPROGversion(char *PROGversion);
int getmaxrevision(void);
int getrevision(char *svn);

/* ------------------ gettokrns ------------------------ */

int gettokens(char *tokens, char **tokenptrs){
  int ntokenptrs;
  char *token;
  int i;

  ntokenptrs=0;
  token=strtok(tokens,",");
  while(token!=NULL){
    tokenptrs[ntokenptrs++]=token;
    token=strtok(NULL,",");
  }
  for(i=0;i<ntokenptrs;i++){
    trim(tokenptrs[i]);
    tokenptrs[i]=trim_front(tokenptrs[i]);
  }
  return ntokenptrs;
}

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char *prog;
  char *arg,*csv;
  char file_in[256],file_out[256];
  FILE *stream_in, *stream_out;
  int in_header;
  int buffer_len, nrows, ncols;
  char *buffer,*labels,**labelptrs;
  char *datalabels,**datalabelptrs;
  int nlabelptrs,ndatalabelptrs;
  int *transfer,ntransfer,itransfer;
  int i;

  prog=argv[0];

  if(argc==1){
   version("devconvert");
   return 1;
  }

  arg=argv[1];
  csv=strstr(arg,".csv");
  if(csv!=NULL)*csv=0;
  strcpy(file_in,arg);
  strcat(file_in,".csv");
  strcpy(file_out,arg);
  strcat(file_out,"_exp.csv");

  stream_in=fopen(file_in,"r");
  if(stream_in==NULL){
    printf("***error: The file %s could not be opened for input\n",file_in);
    return 1;
  }

  stream_out=fopen(file_out,"w");
  if(stream_out==NULL){
    printf("***error: The file %s could not be opened for output\n",file_out);
    return 1;
  }

  initMM();

  buffer_len=getrowcols(stream_in, &nrows, &ncols);
  buffer_len+=10;
  
  NewMemory((void **)&buffer,buffer_len);
  NewMemory((void **)&labels,buffer_len);
  NewMemory((void **)&labelptrs,buffer_len*sizeof(char *));
  NewMemory((void **)&datalabels,buffer_len);
  NewMemory((void **)&datalabelptrs,buffer_len*sizeof(char *));
  NewMemory((void **)&transfer,buffer_len*sizeof(int));

  if(fgets(labels,buffer_len,stream_in)==NULL){
    printf("***error: The file %s is empty\n",file_in);
    return 1;
  }
  if(fgets(labels,buffer_len,stream_in)==NULL){
    printf("***error: The file %s is empty\n",file_in);
    return 1;
  }
  
  nlabelptrs=gettokens(labels,labelptrs);
  ntransfer=0;
  for(i=0;i<nlabelptrs;i++){
    char *token;

    token=labelptrs[i];
    if(strcmp(token,"time")==0){
      transfer[i]=1;
      ntransfer++;
    }
    else if(strncmp(token,"ws",2)==0){
      transfer[i]=2;
      ntransfer++;
    }
    else if(strncmp(token,"wd",2)==0){
      transfer[i]=3;
      ntransfer++;
    }
    else{
      transfer[i]=0;
    }
  }
  itransfer=0;
  for(i=0;i<nlabelptrs;i++){
    if(transfer[i]!=0){
      fprintf(stream_out,"%s",labelptrs[i]);
      itransfer++;
      if(itransfer!=ntransfer)fprintf(stream_out,",");
    }
  }
  fprintf(stream_out,"\n");
  while(!feof(stream_in)){

    CheckMemory;
    if(fgets(datalabels,buffer_len,stream_in)==NULL)break;
    ndatalabelptrs=gettokens(datalabels,datalabelptrs);
    itransfer=0;
    for(i=0;i<ndatalabelptrs;i++){
      char *token;

      if(transfer[i]==0)continue;
      token=datalabelptrs[i];
      if(transfer[i]==1){
        if(strchr(token,':')!=NULL){
          char *hour, *min, *sec;
          int time=0;

          hour=strtok(token,":");
          min=strtok(NULL,":");
          sec=strtok(NULL,":");
          if(hour!=NULL)time+=3600*atoi(hour);
          if(min!=NULL)time+=60*atoi(min);
          if(sec!=NULL)time+=atoi(sec);
          fprintf(stream_out,"%i",time);
        }
        else{
          fprintf(stream_out,"%s",token);
        }
      }
      if(transfer[i]==2){
        if(strcmp(token,"99.99")==0){
          fprintf(stream_out,"NULL");
        }
        else{
          fprintf(stream_out,"%s",token);
        }
      }
      if(transfer[i]==3){
        if(strcmp(token,"9999")==0){
          fprintf(stream_out,"NULL");
        }
        else{
          fprintf(stream_out,"%s",token);
        }
      }
      itransfer++;
      if(itransfer!=ntransfer)fprintf(stream_out,",");
    }
    fprintf(stream_out,"\n");
  }

  return 0;
}

/* ------------------ usage ------------------------ */

void usage(char *prog){
  char prog_version[100];
  int svn_num;

  getPROGversion(prog_version);  // get version (ie 5.x.z)
  svn_num=getmaxrevision();    // get svn revision number

  printf("\n");
  printf("devconvert %s(%i) - %s\n",prog_version,svn_num,__DATE__);
  printf("  Convert a sodar spreadsheet data file for use by Smokeview:\n\n");
  printf("  %s",prog);
  printf(" prog datafile.csv\n\n");

  printf("where\n\n");

  printf("  datafile.csv  - spreadsheet file to be converted\n");
}

/* ------------------ version ------------------------ */

void version(char *prog){
    char version[100];
    int svn_num;

    getPROGversion(version);  // get Smokeview version (ie 5.x.z)
    svn_num=getmaxrevision();    // get svn revision number
    printf("\n");
    printf("%s\n\n",prog);
    printf("Version: %s\n",version);
    printf("SVN Revision Number: %i\n",svn_num);
    printf("Compile Date: %s\n",__DATE__);
}

/* ------------------ getPROGversion ------------------------ */

void getPROGversion(char *PROGversion){
  strcpy(PROGversion,PROGVERSION);
}

/* ------------------ getmaxrev ------------------------ */

#define MAXREV(cval) max_revision=MAX(getrevision(cval),max_revision)
int getmaxrevision(void){
  int max_revision=0;

  MAXREV(main_revision);
  return max_revision;
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



