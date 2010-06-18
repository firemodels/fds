// $Date$ 
// $Revision$
// $Author$

// convert the FDS cases Linux/OSX csh script to an equivalent Windows bat version 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void trim(char *line);
void usage(char *prog);

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char buffer[1024],*buffptr;
  int i;
  char *filein=NULL,*fileout=NULL,*prog;
  FILE *streamin=NULL,*streamout=NULL;

  buffptr=buffer;
  prog=argv[0];
  for(i=1;i<argc;i++){
    int lenarg;
    char *arg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      switch(arg[1]){
      case 'h':
        usage(prog);
        exit(1);
        break;
      default:
        usage(prog);
        exit(1);
        break;
      }
    }
    else{
      if(filein==NULL){
        filein=arg;
        continue;
      }
      if(fileout==NULL){
        fileout=arg;
        continue;
      }
    }
  }
  if(filein==NULL||fileout==NULL){
    usage(prog);
    exit(1);
  }
  streamin=fopen(filein,"r");
  streamout=fopen(fileout,"w");

  if(streamin==NULL||streamout==NULL){
    if(streamin==NULL){
      fprintf(stderr,"unable to open %s for input\n",filein);
    }
    if(streamout==NULL){
      fprintf(stderr,"unable to open %s for output\n",fileout);
    }
    exit(1);
  }
  while(1){
    if(fgets(buffer,1024,streamin)==NULL)break;
    trim(buffer);
    if(strlen(buffer)==0||buffer[0]=='#')continue;
    if(strncmp(buffer,"$RUNFDS",7)==0){
      fprintf(streamout,"%s %s\n","%RUNFDS%",buffptr+8);
    }
    else{
      fprintf(streamout,"%s\n",buffer);
    }
  }
}

/* ------------------ usage ------------------------ */

void usage(char *prog){
  fprintf(stderr,"%s file_in file_out\n\n",prog);
  fprintf(stderr,"  convert the csh script file file_in to a windows bat file equivalent by\n");
  fprintf(stderr,"  ignoring lines beginning with the csh comment character # and converting \n");
  fprintf(stderr,"  the csh symbol $RUNFDS to the DOS symbol %RUNFDS% \n");
}

/* ------------------ trim ------------------------ */

void trim(char *line){
  char *blank=" ";
  const char *c;
  const char *lf="\n";
  unsigned int len;
  unsigned int i;
  len = strlen(line);
  c = line+len-1;
  for(i=0; i<len; i++){
    if(strncmp(c,blank,1)!=0&&strncmp(c,lf,1)!=0){
      c++; 
      line[c-line]='\0';
      return;
    }
    c--;
  }
  *line='\0';
}

