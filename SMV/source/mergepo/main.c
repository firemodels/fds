// $Date: 2011-05-20 00:21:00 -0400 (Fri, 20 May 2011) $ 
// $Revision: 8337 $
// $Author: gforney $

#define INMAIN

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MALLOC.h"
#include "translate.h"
#include "string_util.h"
     
void usage(char *prog);

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  int i;
  char *arg,*prog, *file1=NULL, *file_template=NULL;

  initMM();
  prog=argv[0];
  for(i=1;i<argc;i++){
    int lenarg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      switch(arg[1]){
      case 'a':
        break;
      default:
        usage(prog);
        return 1;
      }
    }
    else{
      if(file1==NULL){
        file1=argv[i];
      }
      else{
        file_template=argv[i];
      }
    }
  }
  if(file1!=NULL&&file_template!=NULL){
    int return1, return2;
    trdata *trinfo1, *trinfo_template;
    int ntrinfo1, ntrinfo_template;

    return1 = parse_lang(file1, &trinfo1, &ntrinfo1);
    return2 = parse_lang(file_template, &trinfo_template, &ntrinfo_template);
    if(return1==1&&return2==1){
      int i;

      for(i=0;i<ntrinfo_template;i++){
        trdata *tri, *trval;

        tri = trinfo_template + i;
        trval = bsearch(tri,trinfo1,ntrinfo1,sizeof(trdata),compare_trdata);
        if(trval!=NULL){
          tri->value=tri->value;
        }
        else{
          tri->value=NULL;
        }
      }
      qsort(trinfo_template,ntrinfo_template,sizeof(trdata),compare_trdata2);
      for(i=0;i<ntrinfo_template;i++){
        trdata *tri;

        tri = trinfo_template + i;
        printf("\n");
        printf("msgid \"%s\"\n",tri->key);
        if(tri->value!=NULL){
          printf("msgstr \"%s\"\n",tri->value);
        }
        else{
          printf("msgstr \"\"\n");
        }
      }
      
    }

  }
}

/* ------------------ usage ------------------------ */

void usage(char *prog){
}

