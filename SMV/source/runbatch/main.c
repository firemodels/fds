#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void usage(void);
int batchfile_exists(char *filename);

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  int i;
  char *arg;

  if(argc==1){
    usage();
    return 1;
  }
  for(i=1;i<argc;i++){
    arg=argv[i];
    if(arg[0]=='-'){
      if(strlen(arg)>1){
        switch(arg[1]){
          case 'h':
              usage();
              return 1;
              break;
          default:
            printf("Unknown option: %s\n",arg);
            usage();
            return 1;
        }
      }
    }
    else{
      break;
    }
  }
  if(batchfile_exists(arg)==1){
    system(arg);
  }
  else{
    if(arg!=NULL&&strlen(arg)>0){
      fprintf(stderr,"ERROR: the batch file, %s, does not exist.\n",arg);
    }
    else{
      fprintf(stderr,"ERROR: batch file does not exist.\n");
    }
  }
}

/* ------------------ usage ------------------------ */

void usage(void){
  printf("\n");
  printf("runbatch batchfile\n");
  printf("  Runs a batch file as an executable\n\n");
}

/* ------------------ batchfile_exists ------------------------ */

int batchfile_exists(char *filename){
  char batchfilename[1024];

  if(filename==NULL)return 0;
  strcpy(batchfilename,filename);
  if(strstr(filename,".bat")==NULL)strcat(batchfilename,".bat");
  if(_access(batchfilename,0)==-1){
    return 0;
  }
  else{
    return 1;
  }
}
