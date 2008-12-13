// $Date: 2008-12-12 16:54:31 -0500 (Fri, 12 Dec 2008) $ 
// $Revision: 2842 $
// $Author: gforney $

// colcat merges 3 files together.  The results of the merge
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
  char buffer1[1024];
  char buffer2[1024];
  char buffer3[1024];
  FILE *stream1, *stream2, *stream3;
  char *file1, *file2, *file3;

  if(argc<4){
    fprintf(stderr,"aborting: 3 input files required (only %i found)\n",argc-1);
    exit(1);
  }
  file1 = argv[1];
  file2 = argv[2];
  file3 = argv[3];
  stream1=fopen(file1,"r");
  stream2=fopen(file2,"r");
  stream3=fopen(file3,"r");
  if(stream1==NULL||stream2==NULL||stream3==NULL){
    fprintf(stderr,"aborting ...\n");
    if(stream1==NULL){
      fprintf(stderr,"  The file %s could not be opened for input.\n",file1);
    }
    if(stream2==NULL){
      fprintf(stderr,"  The file %s could not be opened for input.\n",file2);
    }
    if(stream3==NULL){
      fprintf(stderr,"  The file %s could not be opened for input.\n",file3);
    }
    exit(1);
  }
  while(1){
//  as soon as one file ends the merging stops
    if(fgets(buffer1,1024,stream1)==NULL)break;
    if(fgets(buffer2,1024,stream2)==NULL)break;
    if(fgets(buffer3,1024,stream3)==NULL)break;
    trim(buffer1);
    trim(buffer2);
    trim(buffer3);
    printf("%s,%s,%s\n",buffer1,buffer2,buffer3);
  }
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

