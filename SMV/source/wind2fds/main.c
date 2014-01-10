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
  char *c_dateptr=NULL, c_date[256];
  char *c_mindateptr=NULL, c_mindate[256];
  char *c_maxdateptr=NULL, c_maxdate[256];
  int have_mintime=0, have_maxtime=0;
  unsigned int i_mindate, i_maxdate;
  int lendate=0;

  strcpy(percen,"%");
  strcpy(prefix,"");

  prog=argv[0];

  if(argc==1){
   version("wind2fds");
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
      c_dateptr=c_date;
      i++;
      if(i>=argc)continue;
      arg=argv[i];
      strcpy(c_date,arg);
      lendate=strlen(c_date);
      continue;
    }
    else if(strcmp(arg,"-mindate")==0){
      i++;
      if(i>argc)continue;
      arg=argv[i];
      c_mindateptr=c_mindate;
      strcpy(c_mindateptr,arg);
      if(strchr(c_mindateptr,':')!=NULL){
        have_mintime=1;
        i_mindate=date2sec(c_mindateptr);
      }
      else{
        i_mindate=date2day(c_mindateptr);
      }
      continue;
    }
    else if(strcmp(arg,"-maxdate")==0){
      i++;
      if(i>argc)continue;
      arg=argv[i];
      c_maxdateptr=c_maxdate;
      strcpy(c_maxdateptr,arg);
      if(strchr(c_maxdateptr,':')!=NULL){
        have_maxtime=1;
        i_maxdate=date2sec(c_maxdateptr);
      }
      else{
        i_maxdate=date2day(c_maxdateptr);
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
    else if(strcmp(arg,"-h")==0||strcmp(arg,"-help")==0){
      usage(prog);
      exit(0);
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
    fprintf(stderr,"*** Error: An input file was not specified\n");
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
    fprintf(stderr,"*** Error: The file %s could not be opened for input\n",file_in);
    return 1;
  }

  stream_out=fopen(file_out,"w");
  if(stream_out==NULL){
    fprintf(stderr,"*** Error: The file %s could not be opened for output\n",file_out);
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
    fprintf(stderr,"*** Error: The file %s is empty\n",file_in);
    return 1;
  }
  if(is_sodar_file==1){
    while(strncmp(labels,"Sodar",5)==0){
      if(fgets(labels,buffer_len,stream_in)==NULL){
        fprintf(stderr,"*** Error: The file %s is empty\n",file_in);
        return 1;
      }
    }
  }
  else{
    if(fgets(labels,buffer_len,stream_in)==NULL){
      fprintf(stderr,"*** Error: The file %s is empty\n",file_in);
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
  if(c_dateptr!=NULL){
    fprintf(stream_out,"  //  date: %s\n",c_dateptr);
  }
  if(c_mindateptr!=NULL){
    fprintf(stream_out,"  //  mindate: %s\n",c_mindateptr);
  }
  if(c_maxdateptr!=NULL){
    fprintf(stream_out,"  //  maxdate: %s\n",c_maxdateptr);
  }
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
      fprintf(stream_out,"s,s");
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
        fprintf(stream_out,"%s,time_orig",labelptrs[i]);
      }
      else{
        fprintf(stream_out,"%s",token2);
      }
      itransfer++;
      if(itransfer!=ntransfer)fprintf(stream_out,",");
    }
  }
  if(is_sodar_file==0){
    if(fgets(labels,buffer_len,stream_in)==NULL){
      fprintf(stderr,"*** Error: The file %s is empty\n",file_in);
      return 1;
    }
    if(fgets(labels,buffer_len,stream_in)==NULL){
      fprintf(stderr,"*** Error: The file %s is empty\n",file_in);
      return 1;
    }
  }
  fprintf(stream_out,"\n");
  while(!feof(stream_in)){
    int skip_time;

    CheckMemory;
    if(fgets(datalabels,buffer_len,stream_in)==NULL)break;
    ndatalabelptrs=gettokens(datalabels,datalabelptrs);
    itransfer=0;
    skip_time=0;
    for(i=0;i<ndatalabelptrs;i++){
      char *token;

      if(transfer[i]==0)continue;
      token=datalabelptrs[i];
      if(transfer[i]==1){
        if(strchr(token,':')!=NULL){
          unsigned int time_local=0;

          if(c_dateptr!=NULL&&strncmp(c_dateptr,token,lendate)!=0){
            skip_time=1;
            break;
          }
          if(c_mindateptr!=NULL){
            if(have_mintime==0&&date2day(token)<i_mindate||have_mintime==1&&date2sec(token)<i_mindate){
              skip_time=1;
              break;
            }
          }
          if(c_maxdateptr!=NULL){
            if(have_maxtime==0&&date2day(token)>i_maxdate||have_maxtime==1&&date2sec(token)>i_maxdate){
              skip_time=1;
              break;
            }
          }
          if(tokenbaseptr==NULL){
            tokenbaseptr=tokenbase;
            strcpy(tokenbase,token);
            time_local=0;
          }
          else{
            time_local = diffdate(token,tokenbaseptr);
          }
          fprintf(stream_out,"%i,%s",time_local,token);
        }
        else{
          fprintf(stream_out,"%s",token);
        }
      }
      if(transfer[i]==2||transfer[i]==4){
        if(strcmp(token,"99.99")==0||strcmp(token,"NAN")==0){
          fprintf(stream_out,"NULL");
        }
        else{
          fprintf(stream_out,"%s",token);
        }
      }
      if(transfer[i]==3||transfer[i]==5){
        if(strcmp(token,"9999")==0||strcmp(token,"NAN")==0){
          fprintf(stream_out,"NULL");
        }
        else{
          fprintf(stream_out,"%s",token);
        }
      }
      itransfer++;
      if(itransfer!=ntransfer)fprintf(stream_out,",");
    }
    if(skip_time==1)continue;
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

  printf("  -h             - displays this message\n");
  printf("  -prefix label  - prefix column headers with label\n");
  printf("  -offset x y z  - offset sensor locations by (x,y,z)\n");
  printf("  -wv            - converting a non-sodar file\n");
  printf("  -date mm/dd/yyyy - only convert data recorded on the mm/dd/yyyy\n");
  printf("  -mindate \"mm/dd/yyyy [hh:mm:ss]\" - ignore data recorded before specified date\n");
  printf("  -maxdate \"mm/dd/yyyy [hh:mm:ss]\" - ignore data recorded after specified date\n");
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
