// $Date: 2012-08-20 19:39:01 -0400 (Mon, 20 Aug 2012) $ 
// $Revision: 12156 $
// $Author: koverholt $

#define INMAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "svn_revision.h"
#include "blockaid.h"
#include "MALLOC.h"

// svn revision character string
char main_revision[]="$Revision: 12156 $";

/* ------------------ main ------------------------ */

int main(int argc, char **argv){

  int i;
  char *arg, *prog;
  char *in_file_base=NULL;

#ifdef WIN32
  strcpy(dirsep,"\\");
#else
  strcpy(dirsep,"/");
#endif
  libdir=NULL;

  prog=argv[0];
  in_file_base=NULL;
  if(argc==1){
    version();
    return 1;
  }

  force_write = 0;

  for(i=1;i<argc;i++){
    size_t lenarg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      switch(arg[1]){
      case 'f':
        force_write=1;
        break;
      case 'h':
        usage();
        return 1;
        break;
      case 'l':
        if(i+1<argc){
          int lenarg;

          lenarg=strlen(argv[i+1]);
          NewMemory((void **)&libdir,lenarg+2);
          strcpy(libdir,argv[i+1]);
          if(libdir[lenarg-1]!=dirsep[0]){
            strcat(libdir,dirsep);
          }
          i++;
        }
        break;

      default:
        usage();
        return 1;
      }
    }
    else{
      if(in_file_base==NULL){
        in_file_base=argv[i];
      }
    }
  }

  // construct input and output filenames
  
  if(in_file_base==NULL){
    usage();
    return 1;
  }

  // make sure smv file name exists


  startup();
  readfds(in_file_base);

}

void startup(void){
  nkeyvalstack=0;
  blockaid_first=&ba_first;
  blockaid_last=&ba_last;
  blockaid_first->prev=NULL;
  blockaid_first->next=blockaid_last;
  blockaid_first->keyword_list=NULL;
  blockaid_first->val_list=NULL;
  blockaid_first->nkeywordlist=0;
  blockaid_last->prev=blockaid_first;
  blockaid_last->next=NULL;
  nblockaid=0;
  NewMemory((void **)&grouplist,MAXRECURSE*sizeof(blockaiddata *));
  NewMemory((void **)&offset_rotate,4*MAXRECURSE*sizeof(float));
}