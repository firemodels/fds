// $Date: 2008-12-12 16:54:31 -0500 (Fri, 12 Dec 2008) $ 
// $Revision: 2842 $
// $Author: gforney $

// copy all but the first line from input to output 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// svn revision character string
char main_revision[]="$Revision: 2842 $";
void trim(char *line);

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char buffer[1024];

  if(fgets(buffer,1024,stdin)==NULL)exit;
  while(1){
    if(fgets(buffer,1024,stdin)==NULL)break;
    trim(buffer);
    printf("%s\n",buffer);
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

