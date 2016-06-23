#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ply.h"

/* ------------------ usage ------------------------ */

void usage(void){
  printf(" ply2fds filename.ply\n\n");
  printf(" -h - display this message\n");
  printf(" -i - output plyfile info\n");
}

/* ------------------ get_ply_info ------------------------ */

void get_ply_info(char *plyfile){
}

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  int i;
  char *plyfile=NULL;
  int GETINFO = 0;

  if(argc==1){
    usage();
    return 1;
  }

  for(i = 1; i<argc; i++){
    int lenarg;
    char *arg;

    arg = argv[i];
    lenarg = strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      switch(arg[1]){
      case 'h':
        usage();
        return 1;
        break;
      case 'i':
        GETINFO = 1;
        break;
      default:
        usage();
        return 1;
        break;
      }
    }
    else{
      if(plyfile==NULL){
        plyfile = argv[i];
      }
    }
  }
  if(plyfile==NULL){
    usage();
    return 1;
  }
  if(GETINFO==1){
    get_ply_info(plyfile);
    return 0;
  }
}
