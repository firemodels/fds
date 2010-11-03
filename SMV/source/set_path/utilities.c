// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include "svn_revision.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "MALLOC.h"
#include "ASSERT.h"

char utilities_revision[]="$Revision$";

/* ------------------ imax ------------------------ */

int imax(int a, int b){
  if(a>b){
    return a;
  }
  else{
    return b;
  }
}

/* ------------------ trim ------------------------ */

void trim(char *line){
  size_t len, i;
  unsigned char *c;
  
  len = strlen(line);
  c = (unsigned char *)(line+len-1);
  for(i=0; i<len; i++){
    if(isspace(*c)==0){
      c++; 
      line[c-(unsigned char *)line]='\0';
      return;
    }
    c--;
  }
  *line='\0';
}

/* ------------------ trim_front ------------------------ */

char *trim_front(char *line){
  unsigned char *c;
  size_t i,len;

  c = (unsigned char *)line;
  len=strlen(line);
  for(i=0;i<len;i++){
    if(isspace(*c)==0)return line+i;
    c++;
  }
  return line+len;
}

/* ------------------ STRCMP ------------------------ */

int STRCMP(const char *s1, const char *s2){
  while (toupper(*s1) == toupper(*s2++)){
		if (*s1++ == 0)return (0);
  }
	return (toupper(*(const unsigned char *)s1) - toupper(*(const unsigned char *)(s2 - 1)));
}

/* ------------------ STRSTR ------------------------ */

char *STRSTR(char *c, const char *key){
  char C[10000],*result, KEY[10000];
  size_t i, len_c,len_key;

  if(c==NULL||key==NULL)return NULL;
  len_c=strlen(c);
  len_key=strlen(key);
  if(len_c<1||len_key<1)return NULL;

  for(i=0;i<len_c;i++){
    C[i]=toupper(c[i]);
  }
  C[len_c]='\0';
  for(i=0;i<len_key;i++){
    KEY[i]=toupper(key[i]);
  }
  KEY[len_key]='\0';
  result = strstr(C,KEY);
  if(result==NULL){
    return NULL;
  }
  else{
    return result + (c - C); 
  }
}

/* ------------------ getrevision ------------------------ */

int getrevision(char *svn){
  char svn_string[256];
  char *svn_ptr;
  int return_val;

  svn_ptr=svn_string;
  svn=strchr(svn,':');
  if(svn==NULL||strlen(svn)<=4)return 0;
  
  svn++;
  strcpy(svn_ptr,svn);
  svn_ptr=trim_front(svn_ptr);
  svn_ptr[strlen(svn_ptr)-1]=0;
  trim(svn_ptr);
  sscanf(svn_ptr,"%i",&return_val);
  return return_val;
}

/* ------------------ getmaxrev ------------------------ */

#define MAXREV(cval) max_revision=imax(getrevision(cval),max_revision)
int getmaxrevision(void){
  int max_revision=0;

  MAXREV(main_revision);
  MAXREV(utilities_revision);
  return max_revision;
}
