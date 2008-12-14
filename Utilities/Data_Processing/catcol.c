// $Date: 2008-12-12 16:54:31 -0500 (Fri, 12 Dec 2008) $ 
// $Revision: 2842 $
// $Author: gforney $

// colcat merges 3 files together.  The result of the merge
//      is written to stdout (so should be re-direected to a file)

//      The merging takes place to the right of a file
//      file rather than at the bottom.  This allows 3
//      separate spread sheets to be combined so that the
//      data from each may be plotted on one  python script.

//      The intent is that this program be made more general
//      to handle any reasonable number of files.  Now it solves
//      a particular problem for the MFRI tower data (plotting
//      3 HRR cases at once)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// svn revision character string
char main_revision[]="$Revision: 2842 $";
void trim(char *line);

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char buffer[1024];
  char **filelist;
  FILE **streams;
  int nfiles;
  int i;

  nfiles = argc-1;
  if(nfiles<1){
    fprintf(stderr,"%s file 1, file 2, ..., file n\n\n",argv[0]);
    fprintf(stderr,"  Concatenates files 2 to the right of 1 then 3 to the\n");
    fprintf(stderr,"  right of 2 and so on assuming that each file is\n");
    fprintf(stderr,"  comma deliminated. The result of the concatenation\n");
    fprintf(stderr,"  is output to stdout\n");
    exit(0);
  }
  
  if(nfiles>10){
    fprintf(stderr,"  Concatenation of more than 10 files not supported\n");
    exit(1);
  }
  streams=malloc((argc-1)*sizeof(FILE *));
  for(i=0;i<nfiles;i++){
    streams[i]=fopen(argv[i+1],"r");
    if(streams[i]==NULL){
      fprintf(stderr,"  The file %s could not be opened for input.\n",argv[i]);
    }
  }

  while(1){
//  as soon as one file ends the merging stops
   
    for(i=0;i<nfiles;i++){ 
      if(fgets(buffer,1024,streams[i])==NULL){
        if(i!=0)printf("\n");
        exit(1);
      }
      trim(buffer);
      printf("%s",buffer);
      if(i==nfiles-1){
        printf("\n");
      }
      else{
        printf(",");
      }
    }
  }
  exit(0);
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

