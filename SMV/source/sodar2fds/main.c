// $Date$ 
// $Revision$
// $Author$

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
char main_revision[]="$Revision$";

void usage(char *prog);
void version(char *prog);
void getPROGversion(char *PROGversion);
int getmaxrevision(void);

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
  char *arg,*csv,*argin=NULL,*argout=NULL;
  char file_in[256],file_out[256];
  FILE *stream_in, *stream_out;
  int buffer_len, nrows, ncols;
  char *buffer,*labels,**labelptrs;
  char *datalabels,**datalabelptrs;
  int nlabelptrs,ndatalabelptrs;
  int *transfer,ntransfer,itransfer;
  float *zdev;
  float xyzoffset[3]={0.0,0.0,0.0};
  int i;
  char prefix[256],percen[2];
  int useprefix=0;
  char coffset[255];

  strcpy(percen,"%");
  strcpy(prefix,"");

  prog=argv[0];

  if(argc==1){
   version("devconvert");
   return 1;
  }

  for(i=1;i<argc;i++){
    int lenarg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(strcmp(arg,"-offset")==0){
      i++;
      if(i>=argc)continue;
      arg=argv[i];
      if(strlen(arg)>1){
        sscanf(arg,"%f %f %f",xyzoffset,xyzoffset+1,xyzoffset+2);
      }
      continue;
    }
    else if(strcmp(arg,"-prefix")==0){
      useprefix=1;
      i++;
      if(i>=argc)continue;
      arg=argv[i];
      strcpy(prefix,arg);
      strcat(prefix,"_");
      continue;
    }
    if(argin==NULL){
      argin=arg;
      continue;
    }
    if(argout==NULL){
      argout=arg;
      continue;
    }
  }

  if(argin==NULL){
    printf("***error: An input file was not specified\n");
    return 1;
  }
  csv=strstr(argin,".csv");
  if(csv!=NULL)*csv=0;
  strcpy(file_in,argin);
  strcat(file_in,".csv");
  if(argout==NULL){
    strcpy(file_out,argin);
    strcat(file_out,"_exp.csv");
  }
  else{
    strcpy(file_out,argout);
  }

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

  initMALLOC();

  buffer_len=getrowcols(stream_in, &nrows, &ncols);
  buffer_len+=10;
  
  NewMemory((void **)&buffer,buffer_len);
  NewMemory((void **)&labels,buffer_len);
  NewMemory((void **)&labelptrs,buffer_len*sizeof(char *));
  NewMemory((void **)&datalabels,buffer_len);
  NewMemory((void **)&datalabelptrs,buffer_len*sizeof(char *));
  NewMemory((void **)&transfer,buffer_len*sizeof(int));
  NewMemory((void **)&zdev,buffer_len*sizeof(float));

  if(fgets(labels,buffer_len,stream_in)==NULL){
    printf("***error: The file %s is empty\n",file_in);
    return 1;
  }
  while(strncmp(labels,"Sodar",5)==0){
    if(fgets(labels,buffer_len,stream_in)==NULL){
      printf("***error: The file %s is empty\n",file_in);
      return 1;
    }
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
      zdev[i]=(float)atof(token+2);
      ntransfer++;
    }
    else if(strncmp(token,"wd",2)==0){
      zdev[i]=(float)atof(token+2);
      transfer[i]=3;
      ntransfer++;
    }
    else if(strncmp(token,"sds",3)==0){
      zdev[i]=(float)atof(token+3);
      transfer[i]=4;
      ntransfer++;
    }
    else if(strncmp(token,"sdwd",4)==0){
      zdev[i]=(float)atof(token+4);
      transfer[i]=5;
      ntransfer++;
    }
    else{
      transfer[i]=0;
    }
  }
  fprintf(stream_out,"//HEADER\n");
  for(i=0;i<nlabelptrs;i++){
    char token2[256];

    strcpy(token2,"");
    if(useprefix==1)strcat(token2,prefix);
    strcat(token2,labelptrs[i]);
    sprintf(coffset,"%f %f %f",xyzoffset[0],xyzoffset[1],xyzoffset[2]+zdev[i]);
    trimmzeros(coffset);
    if(transfer[i]==2){
      fprintf(stream_out,"DEVICE\n");
      fprintf(stream_out," %s %s VELOCITY %s sensor\n",token2,percen,percen);
      fprintf(stream_out," %s\n",coffset);
    }
    else if(transfer[i]==3){
      fprintf(stream_out,"DEVICE\n");
      fprintf(stream_out," %s %s ANGLE %s sensor\n",token2,percen,percen);
      fprintf(stream_out," %s\n",coffset);
    }
    else if(transfer[i]==4){
      fprintf(stream_out,"DEVICE\n");
      fprintf(stream_out," %s %s SD_VELOCITY %s sensor\n",token2,percen,percen);
      fprintf(stream_out," %s\n",coffset);
    }
    else if(transfer[i]==5){
      fprintf(stream_out,"DEVICE\n");
      fprintf(stream_out," %s %s SD_ANGLE %s sensor\n",token2,percen,percen);
      fprintf(stream_out," %s\n",coffset);
    }
  }
  fprintf(stream_out,"//DATA\n");
  itransfer=0;
  for(i=0;i<nlabelptrs;i++){
    if(transfer[i]==1){
      fprintf(stream_out,"s");
      itransfer++;
      if(itransfer!=ntransfer)fprintf(stream_out,",");
    }
    if(transfer[i]==2){
      fprintf(stream_out,"m/s");
      itransfer++;
      if(itransfer!=ntransfer)fprintf(stream_out,",");
    }
    if(transfer[i]==3){
      fprintf(stream_out,"deg");
      itransfer++;
      if(itransfer!=ntransfer)fprintf(stream_out,",");
    }
    if(transfer[i]==4){
      fprintf(stream_out,"m/s");
      itransfer++;
      if(itransfer!=ntransfer)fprintf(stream_out,",");
    }
    if(transfer[i]==5){
      fprintf(stream_out,"deg");
      itransfer++;
      if(itransfer!=ntransfer)fprintf(stream_out,",");
    }
  }
  fprintf(stream_out,"\n");
  itransfer=0;
  for(i=0;i<nlabelptrs;i++){
    char token2[256];

    strcpy(token2,"");
    if(useprefix==1)strcat(token2,prefix);
    strcat(token2,labelptrs[i]);
    if(transfer[i]!=0){
      if(transfer[i]==1){
        fprintf(stream_out,"%s",labelptrs[i]);
      }
      else{
        fprintf(stream_out,"%s",token2);
      }
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
          int time_local=0;

          hour=strtok(token,":");
          min=strtok(NULL,":");
          sec=strtok(NULL,":");
          if(hour!=NULL)time_local+=3600*atoi(hour);
          if(min!=NULL)time_local+=60*atoi(min);
          if(sec!=NULL)time_local+=atoi(sec);
          fprintf(stream_out,"%i",time_local);
        }
        else{
          fprintf(stream_out,"%s",token);
        }
      }
      if(transfer[i]==2||transfer[i]==4){
        if(strcmp(token,"99.99")==0){
          fprintf(stream_out,"NULL");
        }
        else{
          fprintf(stream_out,"%s",token);
        }
      }
      if(transfer[i]==3||transfer[i]==5){
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
  printf("sodar2fds %s(%i) - %s\n",prog_version,svn_num,__DATE__);
  printf("  Convert a sodar spreadsheet data file for use by Smokeview:\n\n");
  printf("  %s",prog);
  printf(" prog [-prefix label] [-offset x y z] datafile\n\n");

  printf("where\n\n");

  printf("  -prefix label  - prefix column headers with label\n");
  printf("  -offset x y z  - offset sensor locations by (x,y,z)\n");
  printf("  datafile.csv   - spreadsheet file to be converted\n");
}

/* ------------------ version ------------------------ */

void version(char *prog){
    char version_local[100];
    int svn_num;

    getPROGversion(version_local);  // get Smokeview version (ie 5.x.z)
    svn_num=getmaxrevision();    // get svn revision number
    printf("\n");
    printf("%s\n\n",prog);
    printf("Version: %s\n",version_local);
    printf("SVN Revision Number: %i\n",svn_num);
    printf("Compile Date: %s\n",__DATE__);
}

/* ------------------ getPROGversion ------------------------ */

void getPROGversion(char *PROGversion){
  strcpy(PROGversion,PROGVERSION);
}

/* ------------------ getmaxrev ------------------------ */

#define MAXREV(cval) rev=getrevision(cval);max_revision=MAX(rev,max_revision)
int getmaxrevision(void){
  int max_revision=0,rev;

  MAXREV(main_revision);
  return max_revision;
}
