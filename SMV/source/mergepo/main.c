// $Date: 2011-05-20 00:21:00 -0400 (Fri, 20 May 2011) $ 
// $Revision: 8337 $
// $Author: gforney $

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MALLOC.h"

     
typedef struct {
  char *key, *value;
} trdata;

void usage(char *prog);
void trim(char *line);
int compare_trdata( const void *arg1, const void *arg2 );
char *trim_front(char *line);
char *getstring(char *buffer);
int parse_lang(char *file, trdata **trinfoptr, int *ntrinfoptr);

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  int i;
  char *arg,*prog, *file1=NULL, *file_template=NULL;

  initMM();
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
      if(file1==NULL){
        file1=argv[i];
      }
      else{
        file_template=argv[i];
      }
    }
  }
  if(file1!=NULL&&file_template!=NULL){
    int return1, return2;
    trdata *trinfo1, *trinfo_template;
    int ntrinfo1, ntrinfo_template;

    return1 = parse_lang(file1, &trinfo1, &ntrinfo1);
    return2 = parse_lang(file_template, &trinfo_template, &ntrinfo_template);
    if(return1==1&&return2==1){
      int i;

      for(i=0;i<ntrinfo_template;i++){
        trdata *tri, *trval;

        tri = trinfo_template + i;
        trval = bsearch(tri,trinfo1,ntrinfo1,sizeof(trdata),compare_trdata);
        printf("\n");
        printf("msgid \"%s\"\n",tri->key);
        if(trval!=NULL&&trval->value!=NULL){
          printf("msgstr \"%s\"\n",trval->value);
        }
        else{
          printf("msgstr \"\"\n");
        }
      }
    }

  }
}

/* ------------------ trim_front ------------------------ */

char *trim_front(char *line){
  char *blank=" ";
  const char *c;
  size_t i,len;

  c = line;
  len=strlen(line);
  for(i=0;i<len;i++){
    if(strncmp(c++,blank,1)!=0)return line+i;
  }
  return line;
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


/* ------------------ isocompare ------------------------ */

int compare_trdata( const void *arg1, const void *arg2 ){
  trdata *tri, *trj;

  tri = (trdata *)arg1;
  trj = (trdata *)arg2;

  return strcmp(tri->key,trj->key);
}

/* ------------------ parse_lang ------------------------ */

int parse_lang(char *file, trdata **trinfoptr, int *ntrinfoptr){
  FILE *stream;
  char buffer[1000];
  char *buf;
  char *key, *value;
  trdata *trinfo;
  int ntrinfo;


  ntrinfo=0;
  stream=fopen(file,"r");
  if(stream==NULL)return 0;

  while(!feof(stream)){
    if(fgets(buffer,1000,stream)==NULL)break;
    buf=trim_front(buffer);
    trim(buf);
    key = strstr(buf,"msgid");
    if(key!=NULL&&key==buf)ntrinfo++;
  }
  if(ntrinfo==0)return 0;

  NewMemory((void **)&trinfo,sizeof(trdata)*ntrinfo);

  ntrinfo=0;
  rewind(stream);
  while(!feof(stream)){
    trdata *tri;
    char *keybuf,*valbuf;

    if(fgets(buffer,1000,stream)==NULL)break;
    buf=trim_front(buffer);
    trim(buf);
    key = strstr(buf,"msgid");
    if(key==NULL||key!=buf)continue;
    ntrinfo++;
    tri = trinfo + ntrinfo - 1;
    key = getstring(key+5);
    if(key==NULL){
      tri->key=NULL;
    }
    else{
      NewMemory((void **)&keybuf,strlen(key)+1);
      strcpy(keybuf,key);
      tri->key=keybuf;
    }

    if(fgets(buffer,1000,stream)==NULL)break;
    buf=trim_front(buffer);
    trim(buf);
    value = getstring(buf+6);
    if(value==NULL){
      tri->value=value;
    }
    else{
      NewMemory((void **)&valbuf,strlen(value)+1);
      strcpy(valbuf,value);
      tri->value=valbuf;
    }
  }
  fclose(stream);

  qsort(trinfo,ntrinfo,sizeof(trdata),compare_trdata);
  *ntrinfoptr=ntrinfo;
  *trinfoptr=trinfo;
  return 1;
}

/* ------------------ getstring ------------------------ */

char *getstring(char *buffer){
  char *begin=NULL;
  int len,i;
// look for a quoted string in buffer -  "la;kjsf;lkjaf"
  
  len=strlen(buffer);

  for(i=0;i<len;i++){
    if(buffer[i]=='"'){
      begin=buffer+i+1;
      break;
    }
  }
  if(begin==NULL)return NULL;
  for(i=begin-buffer+1;i<len;i++){
    if(buffer[i]=='"'){
      buffer[i]=0;
      break;
    }
  }
  if(i==len)buffer[i-1]=0;
  for(i=0;i<strlen(begin);i++){
    if(begin[i]!=' ')return begin;
  }
  return NULL;
}

/* ------------------ usage ------------------------ */

void usage(char *prog){
}

