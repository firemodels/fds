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

//dummy change to bump version number to  0.9

// svn revision character string
char main_revision[]="$Revision$";

void usage(char *prog);
void version(char *prog);
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

/* ------------------ daytime2sec ------------------------ */

int daytime2sec(char *tokenorig){          
  char token[256];
  char *hour=NULL, *min=NULL, *sec=NULL;
  char *month=NULL, *day=NULL, *year=NULL;
  int imonth, iday, iyear=2000, ileap;
  int time_local;
  int days_local;
  int month2days[]={0,31,59,90,120,151,181,212,243,273,304,334};
#define SECS_IN_DAY 86400
  char *slash1=NULL, *slash2=NULL, *colen1=NULL, *colen2=NULL;

  strcpy(token,tokenorig); 

  slash1=strchr(token,'/');
  if(slash1!=NULL)slash2=strchr(slash1+1,'/');
  colen1=strchr(token,':');
  if(colen1!=NULL)colen2=strchr(colen1+1,':');

  if(slash1==NULL){
    hour=token;
  }
  else if(slash1!=NULL&&slash2==NULL){
    char *dayend;

    month=token;
    day=slash1+1;
    dayend=strchr(day,' ');
    if(dayend!=NULL)*dayend=0;
    hour=dayend+1;
  }
  else{
    char *yearend;

    month=token;
    day=slash1+1;
    year=slash2+1;
    yearend=strchr(year,' ');
    if(yearend!=NULL)*yearend=0;
    hour=yearend+1;
  }
  if(colen1!=NULL){
    min=colen1+1;
    if(colen2==NULL){
      char *minend;

      minend=strchr(min,' ');
      if(minend!=NULL)*minend=0;
    }
    else{
      char *secend;

      sec=colen2+1;
      secend=strchr(sec,' ');
      if(secend!=NULL)*secend=0;
    }
  }
  if(slash1!=NULL)*slash1=0;
  if(slash2!=NULL)*slash2=0;
  if(colen1!=NULL)*colen1=0;
  if(colen2!=NULL)*colen2=0;


  days_local=0;
  time_local=0;
  if(month!=NULL){
    iyear = atoi(year)-2000;
    imonth = atoi(month);
    iday = atoi(day);
    ileap = iyear/4 + 1;
    if(iyear%4==0&&imonth<3)ileap--;
    days_local += iyear*365;
    days_local += month2days[imonth-1];
    days_local += iday - 1 +ileap;
    time_local += SECS_IN_DAY*days_local;
  }
  if(hour!=NULL)time_local+=3600*atoi(hour);
  if(min!=NULL)time_local+=60*atoi(min);
  if(sec!=NULL)time_local+=atoi(sec);
  return time_local;
}

/* ------------------ diffdate ------------------------ */

int diffdate(char *token, char *tokenbase){
  int difft;

  difft = daytime2sec(token) - daytime2sec(tokenbase);
  return difft;
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
  int is_sodar_file=1;
  char tokenbase[256], *tokenbaseptr=NULL;
  char *datelabelptr=NULL, datelabel[256];

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
    if(strcmp(arg,"-wv")==0){
      is_sodar_file=0;
      continue;
    }
    if(strcmp(arg,"-offset")==0){
      i++;
      if(i>=argc)continue;
      arg=argv[i];
      if(strlen(arg)>1){
        sscanf(arg,"%f %f %f",xyzoffset,xyzoffset+1,xyzoffset+2);
      }
      continue;
    }
    else if(strcmp(arg,"-date")==0){
      datelabelptr=datelabel;
      i++;
      if(i>=argc)continue;
      arg=argv[i];
      strcpy(datelabel,arg);
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
  if(is_sodar_file==1){
    while(strncmp(labels,"Sodar",5)==0){
      if(fgets(labels,buffer_len,stream_in)==NULL){
        printf("***error: The file %s is empty\n",file_in);
        return 1;
      }
    }
  }
  else{
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
    else if(strcmp(token,"TIMESTAMP")==0){
      transfer[i]=1;
      ntransfer++;
    }
    else if(strncmp(token,"ws",2)==0){
      transfer[i]=2;
      zdev[i]=(float)atof(token+2);
      ntransfer++;
    }
    else if(strncmp(token,"WS",2)==0){
      transfer[i]=2;
      zdev[i]=0.0;
      ntransfer++;
    }
    else if(strncmp(token,"wd",2)==0){
      zdev[i]=(float)atof(token+2);
      transfer[i]=3;
      ntransfer++;
    }
    else if(strcmp(token,"WindDir")==0){
      zdev[i]=0.0;
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
          int time_local=0;

          if(tokenbaseptr==NULL){
            tokenbaseptr=tokenbase;
            strcpy(tokenbase,token);
            time_local=0;
          }
          else{
            time_local = diffdate(token,tokenbaseptr);
          }
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
  printf("wind2fds %s(%i) - %s\n",prog_version,svn_num,__DATE__);
  printf("  Convert spreadheets containing wind data to files compatible with Smokeview:\n\n");
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

/* ------------------ getmaxrev ------------------------ */

int getmaxrevision(void){
  int max_revision=0,rev;

  MAXREV(main_revision);
  return max_revision;
}
