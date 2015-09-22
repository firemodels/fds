// convert the Linux/OSX script containing a list FDS cases 
// to an equivalent Windows bat version 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

void trim(char *line);
void usage(char *prog);
char *trim_front(char *line);

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char buffer[1024],*buffptr;
  char buffer2[1024];
  int i;
  char *filein=NULL,*fileout=NULL,*prog;
  FILE *streamin=NULL,*streamout=NULL;
  int lendata;

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
  fprintf(streamout,"@echo off\n");
  for(;;){
    if(fgets(buffer,1024,streamin)==NULL)break;
    trim(buffer);
    if(strlen(buffer)==0){
      fprintf(streamout,"\n");
      continue;
    }
    if(buffer[0]=='#'){
      if(strlen(buffer)>0){
        fprintf(streamout,":: %s\n",buffer+1);
      }
      else{
        fprintf(streamout,"::\n");
      }
      continue;
    }
    if(buffer[0]=='$'){
      char *comm_beg, *comm_end, *data;
      char *casename;
      int j;
      char *datato, *datafrom;
      
      comm_beg=buffer+1;
      comm_end=strchr(buffer,' ');
      data = comm_end+1;
      *comm_end=0;

      trim(data);
      fprintf(streamout,"%s%s%s %s\n","%",comm_beg,"%",data);
      continue;

    }
    fprintf(streamout,"%s\n",buffer);
  }
}

/* ------------------ usage ------------------------ */

void usage(char *prog){
  char RUNFDS[32];

  strcpy(RUNFDS,"%RUNFDS%");
  fprintf(stderr,"%s file_in file_out\n\n",prog);
  fprintf(stderr,"  convert the Linux/OSX script file file_in to a windows bat file equivalent by\n");
  fprintf(stderr,"  ignoring lines beginning with the comment character # and converting \n");
  fprintf(stderr,"  the Linux/OSX symbol $RUNFDS to the DOS symbol %s\n",RUNFDS);
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

/* ------------------ trim_front ------------------------ */

char *trim_front(char *line){

// returns a pointer to the first non-blank character in the character string line

  char *c;

  for(c=line;c<=line+strlen(line)-1;c++){
    if(!isspace(*c))return c;
  }
  return line;
}
