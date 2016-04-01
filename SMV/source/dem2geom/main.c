#define INMAIN
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "string_util.h"


/* ------------------ usage ------------------------ */

void usage(char *prog){
 char githash[256];
 char gitdate[256];

  getGitInfo(githash,gitdate);    // get githash

  fprintf(stderr, "\n%s (%s) %s\n", prog, githash, __DATE__);
}

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char buffer[1024],*buffptr;
  int i;
  char *filein=NULL,*fileout=NULL,*prog;
  FILE *streamin=NULL,*streamout=NULL;

  set_stdout(stdout);
  buffptr=buffer;
  prog=argv[0];
  for(i=1;i<argc;i++){
    int lenarg;
    char *arg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      switch(arg[1]){
      case 'h':
        usage(prog);
        exit(1);
        break;
      case 'v':
        version("dem2geom");
        exit(1);
        break;
      default:
        usage(prog);
        exit(1);
        break;
      }
    }
  }
  if(filein==NULL||fileout==NULL){
    usage(prog);
    exit(1);
  }
}
