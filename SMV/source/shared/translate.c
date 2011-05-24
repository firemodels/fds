// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#define IN_TRANSLATE
#include "pragmas.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include "MALLOC.h"
#include "translate.h"
#include "string_util.h"

// svn revision character string
char translate_revision[]="$Revision$";
extern int tr_english;

/* ------------------ isocompare ------------------------ */

int compare_trdata( const void *arg1, const void *arg2 ){
  trdata *tri, *trj;

  tri = (trdata *)arg1;
  trj = (trdata *)arg2;

  return strcmp(tri->key,trj->key);
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

/* ------------------ parse_lang ------------------------ */

int parse_lang(char *file){
  FILE *stream;
  char buffer[1000];
  char *buf;
  char *key, *value;


  ntrinfo=0;
  FREEMEMORY(trinfo);
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
  return 1;
}

/* ------------------ init_translate ------------------------ */

void init_translate(char *smokeviewbindir, char *tr_name){
  char *LANG;

  if(tr_name!=NULL){
    LANG=tr_name;
  }
  else{
    LANG=getenv("LANG");
  }

  tr_lang=0;
  if(LANG!=NULL&&strncmp(LANG,"en",2)!=0){
    FILE *stream;
    int lensmokebindir;
    char lang[256];

    lang[0]=tolower(LANG[0]);
    lang[1]=tolower(LANG[1]);
    lang[2]=0;
    lensmokebindir=strlen(smokeviewbindir);
    NewMemory((void **)&smokeview_lang,(unsigned int)(lensmokebindir+15+1));
    STRCPY(smokeview_lang,smokeviewbindir);
    STRCAT(smokeview_lang,"smokeview_");
    STRCAT(smokeview_lang,lang);
    STRCAT(smokeview_lang,".po");
    
    stream=fopen(smokeview_lang,"r");
    if(stream!=NULL){
      fclose(stream);
      tr_lang=1;
    }
    tr_lang=parse_lang(smokeview_lang);
  }
}

/* ------------------ translate2 ------------------------ */

char *translate2(char *string){
  char *value=NULL;
  trdata *trval;
  trdata trkey;

  trkey.key=string;

  trval = bsearch(&trkey,trinfo,ntrinfo,sizeof(trdata),compare_trdata);
  if(trval!=NULL)value=trval->value;
  return value;
}

/* ------------------ translate ------------------------ */

char *translate(char *string){
  char *valin,*valout;
  int i;

  if(tr_lang!=0&&tr_english==0){
    int len;

    len = strlen(string);
    for(i=0;i<len;i++){
      valin=string+i;
      if(*valin!=' '&&*valin!='*')break;
      tr_string[i]=*valin;
    }
    valout=translate2(valin);
    if(valout==NULL)return string;
    strcpy(tr_string+i,valout);
    return tr_string;
  }
  return string;
}
