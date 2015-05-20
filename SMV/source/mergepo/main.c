#define INMAIN

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "MALLOC.h"
#include "translate.h"
     
void usage(char *prog);

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  int i;
  char *arg,*prog, *file_lang=NULL, *file_template=NULL;
  int add_comments=0;

  initMALLOC();
  prog=argv[0];
  for(i=1;i<argc;i++){
    int lenarg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      switch(arg[1]){
      case 'c':
        add_comments=1;
        break;
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
  if(add_comments==1){
    printf("// Translate the phrase after each msgid and place\n");
    printf("// it on the following line after msgstr, as in:\n");
    printf("// msgid \"original phrase\"\n");
    printf("// msgstr \"translated phrase\"\n");
    printf("// \n");
    printf("// If you have trouble with some of the terms feel free\n");
    printf("// to ask.  It is also OK to leave troublesome terms\n");
    printf("// untranslated, they will be output in English\n");
    printf("// \n");
    printf("// Suggested translation priorities:\n");
    printf("// 1.  translate terms in menus, e.g. Load/Unload, \n");
    printf("//     Show/Hide etc.\n");
    printf("// 2.  translate terms in dialog boxes.\n");
    printf("// 3.  Now go through the following list and translate terms.\n");
    printf("//     that are left.  \n");
    printf("// \n");
  }
  if(file_lang!=NULL&&file_template!=NULL){
    int return1, return2;
    trdata *trinfo_lang, *trinfo_template;
    int ntrinfo_lang, ntrinfo_template;

    return1 = parse_lang(file_lang, &trinfo_lang, &ntrinfo_lang);
    return2 = parse_lang(file_template, &trinfo_template, &ntrinfo_template);
    if(return1==1&&return2==1){
      printf("\n");
      for(i=0;i<ntrinfo_template;i++){
        trdata *tri, *tr_otherlangi;

        tri = trinfo_template + i;
        // foreach string in template look for a translation and store it in lang
        tr_otherlangi = bsearch(tri,trinfo_lang,ntrinfo_lang,sizeof(trdata),compare_trdata);
        if(tr_otherlangi!=NULL){
          tri->value=tr_otherlangi->value;
        }
        else{
          tri->value=NULL;
        }
        printf("msgid \"%s\"\n",tri->key);
        if(tri->value!=NULL){
          printf("msgstr \"%s\"\n",tri->value);
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

