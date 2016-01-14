#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------ usage ------------------------ */

void usage(void){
  printf(" ply2fds filename.ply\n\n");
  printf(" -h - display this message\n");
}

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  int i;
  char *filebase=NULL;

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
      default:
        usage();
        return 1;
        break;
      }
    }
    else{
      if(filebase==NULL){
        filebase = argv[i];
      }
    }
  }
  if(filebase==NULL){
    usage();
    return 1;
  }

}
