// $Date$ 
// $Revision$
// $Author$

#define INMAIN
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <process.h>
#include <windows.h>
#include "svn_revision.h"
#include "background.h"

//dummy change to bump version number to 0.9

// svn revision character string
char main_revision[]="$Revision$";
void run_command(void);

void usage(char *prog);

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char *prog;
  int i;
  int argstart=-1;
  float delay_time=0.0;
  int itime;
  char *arg;
  char *command;
  char *command_arg;

  prog=argv[0];

  if(argc==1){
    version();
    return 1;
  }

  for(i=1;i<argc;i++){
    int lenarg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'){
      if(lenarg>1){
        switch(arg[1]){
	        case 'd':
            i++;
            if(i<argc){
              arg=argv[i];
              sscanf(arg,"%f",&delay_time);
              if(delay_time<0.0)delay_time=0.0;
            }
		        break;
          case 'h':
            usage(prog);
            return 1;
            break;
          default:
            printf("Unknown option: %s\n",arg);
            usage(prog);
            return 1;
            break;
	      }
      }
    }
    else{
      argstart=i;
      break;

    }
  }

  if(argstart<0)return 0;

  itime = delay_time*1000;
  command=argv[argstart];
  if(delay_time>0.0){
    printf("In %f seconds, executing command: %s\n",delay_time,command);
  }
  else{
    printf("Executing command: %s\n",command);
  }

  Sleep(itime);
  printf("before spawn\n");

  _spawnvp(_P_NOWAIT,command, argv+argstart);


  printf("background exiting\n");

  return 0;
}

void usage(char *prog){
  char prog_version[100];
  int svn_num;

  getPROGversion(prog_version);  // get version (ie 5.x.z)
  svn_num=getmaxrevision();    // get svn revision number

  printf("\n");
  printf("  background %s(%i) - %s\n\n",prog_version,svn_num,__DATE__);
  printf("  Runs a windows job in the background\n\n");
  printf("  %s",prog);
  printf(" [-d delay time (s) ] prog prog_arguments\n\n");

  printf("  -d dtime - wait dtime seconds before running program prog in the background\n");
}
