// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>
#include <string.h>
#include "string_util.h"
     
int add_msgstring=0;
void usage(char *prog);

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char buffer[1024];
  int ii;
  char *arg,*prog;
  FILE *stream;

  //stream=fopen("menus.c","r");
  stream=stdin;
  prog=argv[0];
  for(ii=1;ii<argc;ii++){
    int lenarg;

    arg=argv[ii];
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
    while(!feof(stream)){
      char *beg,*end, *beg2;

      fgets(buffer,sizeof(buffer),stream);
      beg=strstr(buffer,"_(\"");
      if(beg==NULL)continue;
      beg+=2;
      for(beg2=beg+1;beg2<buffer+sizeof(buffer);beg2++){
        char c;
   
        c = *beg2;
        if((c<'a'||c>'z')&&(c<'A'||c>'Z'))continue;
        beg=beg2-1;
        *beg='"';
        break;
      }
      end=strstr(beg+1,"\"");
      if(end!=NULL){
        int i,len;

        end[1]=0;
        len=strlen(beg);
        for(i=len-2;i>=0;i--){
          char c;
   
          c = beg[i];
          if((c<'a'||c>'z')&&(c<'A'||c>'Z'))continue;
          beg[i+1]='"';
          beg[i+2]=0;
          break;
        }
        printf("msgid %s\n",beg);
      }
    }
  }
  else{
    while(!feof(stream)){
      char *beg;

      fgets(buffer,sizeof(buffer),stream);
      trim(buffer);
      beg=strstr(buffer,"msgid");
      if(beg!=NULL&&beg==buffer){
        printf("%s\n",buffer);
        printf("msgstr \"\"\n");
        printf("\n");
      }
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


