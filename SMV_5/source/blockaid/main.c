// $Date$ 
// $Revision$
// $Author$

#define INMAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "svn_revision.h"
#include "blockaid.h"

// svn revision character string
char main_revision[]="$Revision$";

/* ------------------ main ------------------------ */

int main(int argc, char **argv){

  int i;
  char *arg, *prog;
  char *fdsfile=NULL;

  prog=argv[0];
  fdsfile=NULL;
  if(argc==1){
    version();
    return 1;
  }

  for(i=1;i<argc;i++){
    size_t lenarg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      switch(arg[1]){
      case 'b':
        break;
      default:
        usage(prog);
        return 1;
      }
    }
    else{
      if(fdsfile==NULL){
        fdsfile=argv[i];
      }
    }
  }

  // construct smv filename
  
  if(fdsfile==NULL){
    usage(prog);
    return 1;
  }

  // make sure smv file name exists

  if(getfileinfo(fdsfile,NULL,NULL)!=0){
    printf("file: %s does not exist\n",fdsfile);
    return 1;
  }
  startup();
  readfds(fdsfile);

}

void startup(void){
  blockaid_first=&ba_first;
  blockaid_last=&ba_last;
  blockaid_first->prev=NULL;
  blockaid_first->next=blockaid_last;
  blockaid_last->prev=blockaid_first;
  blockaid_last->next=NULL;
  nblockaid=0;
}