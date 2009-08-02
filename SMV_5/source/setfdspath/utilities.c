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
  char *blank=" ";
  const char *c;
  const char *lf="\n", *cr="\r";
  size_t len, i;
  
  len = strlen(line);
  c = line+len-1;
  for(i=0; i<len; i++){
    if(strncmp(c,blank,1)!=0&&strncmp(c,lf,1)!=0&&strncmp(c,cr,1)!=0){
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
  char *blank=" ";
  char *tab="\t";
  const char *c;
  size_t i,len;

  c = line;
  len=strlen(line);
  for(i=0;i<len;i++){
    if(strncmp(c,blank,1)!=0&&strncmp(c,tab,1)!=0)return line+i;
    c++;
  }
  return line;
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
  char *C,*CCOPY,*CC,*cc,*result;
  char *KEY,*KEYCOPY,*KEY2;
  size_t i, len,len2;
  int diff;

  if(c==NULL||key==NULL)return NULL;
  len=strlen(c);
  len2=strlen(key);
  if(len<1||len2<1)return NULL;
  if(NewMemory((void **)&C,(unsigned int)(len+1))==0)return NULL;
  CC=C;
  cc=c;
  CCOPY=C;
  if(NewMemory((void **)&KEY,(unsigned int)(len2+1))==0){
    FreeMemory(C);
    return NULL;
  }
  KEY2=KEY;
  KEYCOPY=KEY;
  for(i=0;i<len;i++){
    *CC++=(char)toupper(*cc++);
  }
  for(i=0;i<len2;i++){
    *KEY2++=(char)toupper(*key++);
  }
  *CC='\0';
  *KEY2='\0';
  result = strstr(C,KEY);
  if(result!=NULL)diff = result - C;
  FREEMEMORY(CCOPY);
  FREEMEMORY(KEYCOPY);
  if(result==NULL)return NULL;
  return c + diff;

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
  MAXREV(assert_revision);
  MAXREV(dmalloc_revision);
  MAXREV(utilities_revision);
  return max_revision;
}
