// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "string_util.h"
     
int add_msgstring=0;
void usage(char *prog);
void trim(char *line);

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char buffer[1024];
  int i;
  char *arg,*prog;

  prog=argv[0];
  for(i=1;i<argc;i++){
    int lenarg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      switch(arg[1]){
      case 'a':
        add_msgstring=1;
        break;
      default:
        usage(prog);
        return 1;
      }
    }
    else{
      usage(prog);
      return 1;
    }
  }
  if(add_msgstring==0){
    while(!feof(stdin)){
      char *beg,*end, *beg2;

      fgets(buffer,sizeof(buffer),stdin);
      beg=strstr(buffer,"_(\"");
      if(beg==NULL)continue;
      beg+=2;
      for(beg2=beg+1;beg2<buffer+sizeof(buffer);beg2++){
        if(*beg2==' '||*beg2=='*'||*beg2=='#')continue;
        beg=beg2-1;
        *beg='"';
        break;
      }
      for(end=beg+1;end<buffer+sizeof(buffer);end++){
        if(*end=='\"')break;
      }
      *(end+1)=0;
      printf("msgid %s\n",beg);
    }
  }
  else{
    while(!feof(stdin)){
      fgets(buffer,sizeof(buffer),stdin);
      trim(buffer);
      printf("\n");
      printf("%s\n",buffer);
      printf("msgstr \"\"\n");
    }
  }
}

/* ------------------ usage ------------------------ */

void usage(char *prog){
  printf("%s [-a] < stdin > stdout \n",prog);
  printf("Create a .po file by parsing a collection of .c/.h/.cpp files\n");
  printf("looking for strings of the form _(\"....\") , outputting each\n");
  printf("string found as \n");
  printf("MSGID \".....\"\n");
  printf("If the -a option is  used then the string\n");
  printf("MSGSTR \"\"\n");
  printf("is also output\n");
}


