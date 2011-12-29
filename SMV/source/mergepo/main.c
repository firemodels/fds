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
  char *arg,*prog, *file_lang=NULL, *file_template=NULL;

  initMALLOC();
  prog=argv[0];
  for(i=1;i<argc;i++){
    int lenarg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      switch(arg[1]){
      default:
        usage(prog);
        return 1;
      }
    }
    else{
      if(file_lang==NULL){
        file_lang=argv[i];
      }
      else{
        file_template=argv[i];
      }
    }
  }
  if(file_lang==NULL||file_template==NULL){
    usage(prog);
    return 1;
  }
  if(file_lang!=NULL&&file_template!=NULL){
    int return1, return2;
    trdata *trinfo_lang, *trinfo_template;
    int ntrinfo_lang, ntrinfo_template;

    return1 = parse_lang(file_lang, &trinfo_lang, &ntrinfo_lang);
    return2 = parse_lang(file_template, &trinfo_template, &ntrinfo_template);
    if(return1==1&&return2==1){
      int i;

      printf("\n");
      for(i=0;i<ntrinfo_template;i++){
        trdata *tr_template_i, *tr_otherlang;

        tr_template_i = trinfo_template + i;
        // foreach string in template look for a translation and store it in lang
        tr_otherlang = bsearch(tr_template_i,trinfo_lang,ntrinfo_lang,sizeof(trdata),compare_trdata);
        if(tr_otherlang!=NULL){
          tr_template_i->value=tr_otherlang->value;
        }
        else{
          tr_template_i->value=NULL;
        }
        printf("msgid \"%s\"\n",tr_template_i->key);
        if(tr_template_i->value!=NULL){
          printf("msgstr \"%s\"\n",tr_template_i->value);
        }
        else{
          printf("msgstr \"\"\n");
        }
        if(i!=ntrinfo_template-1)printf("\n");
      }
    }
  }
}

/* ------------------ usage ------------------------ */

void usage(char *prog){
  printf("%s smokeview_template.po smokeview_xx.po\n",prog);
  printf("Merge two .po files, typically smokeview_template.po with smokeview_xx.po\n");
  printf("where xx is the language being translated.\n");
  printf("\n");
  printf("The updated .po files is output to stdout\n");
}

