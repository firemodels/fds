// $Date$ 
// $Revision$
// $Author$

#define IN_TRANSLATE
#include "pragmas.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "translate.h"
// svn revision character string
char translate_revision[]="$Revision$";

char *translate2(char *string){
  return string;
}

/* ------------------ translate ------------------------ */

char *translate(char *string){
  char *valin,*valout;
  int i;

  if(tr_lang!=0){
    int len;

    len = strlen(string);
    for(i=0;i<len;i++){
      valin=string++;
      if(*valin!=' '&&*valin!='*')break;
      tr_string[i]=*valin;
    }
    valout=translate2(valin);
    strcpy(tr_string+i,valout);
    return tr_string;
  }
  return string;
}
