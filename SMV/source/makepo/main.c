// $Date$ 
// $Revision$
// $Author$

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

       
void usage(char *prog);

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
  while(!feof(stdin)){
    char *beg,*end;

    fgets(buffer,sizeof(buffer),stdin);
    beg=strstr(buffer,"_(\"");
    if(beg==NULL)continue;
    beg+=2;
    for(end=beg+1;end<buffer+sizeof(buffer);end++){
      if(*end=='\"')break;
    }
    *(end+1)=0;
    printf("msgid %s\n",beg);
  }
}

/* ------------------ usage ------------------------ */

void usage(char *prog){
}

