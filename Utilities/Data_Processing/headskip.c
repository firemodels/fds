// $Date: 2008-12-12 16:54:31 -0500 (Fri, 12 Dec 2008) $ 
// $Revision: 2842 $
// $Author: gforney $

// copy all but the first n from input to output 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// svn revision character string
char main_revision[]="$Revision: 2842 $";
void trim(char *line);
void usage(char *prog);

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char buffer[1024];
  int i,skip=0;
  char *filebase=NULL,*prog;
  FILE *stream=NULL;

  prog=argv[0];
  for(i=1;i<argc;i++){
    int lenarg;
    char *arg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      switch(arg[1]){
      case 's':
        i++;
        if(i<argc){
          sscanf(argv[i],"%i",&skip);
          if(skip<0)skip=0;
        }
        break;
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
      if(filebase==NULL){
        filebase=argv[i];
      }
    }
  }
  if(filebase==NULL){
    stream=stdin;
  }
  else{
    stream=fopen(filebase,"r");
  }
  if(stream==NULL){
    if(filebase!=NULL){
      fprintf(stderr,"unable to open %s for input\n",filebase);
    }
    else{
      fprintf(stderr,"unable to read from stdin\n");
    }
    exit(1);
  }
  for(i=0;i<skip;i++){
    if(fgets(buffer,1024,stdin)==NULL)exit;
  }
  while(1){
    if(fgets(buffer,1024,stdin)==NULL)break;
    trim(buffer);
    printf("%s\n",buffer);
  }
}

/* ------------------ usage ------------------------ */

void usage(char *prog){
  fprintf(stderr,"%s [-s n] < file_in > file_out\n\n",prog);
  fprintf(stderr,"  copy file_in to file_out while optionally skipping\n");
  fprintf(stderr,"  the first n lines of file_in\n");
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

