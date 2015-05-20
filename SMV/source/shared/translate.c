#include "options.h"
#define IN_TRANSLATE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "MALLOC.h"
#include "translate.h"

/* ------------------ compare_trdata ------------------------ */

int compare_trdata( const void *arg1, const void *arg2 ){
  trdata *tri, *trj;
  int compval;

  tri = (trdata *)arg1; 
  trj = (trdata *)arg2;

  compval = STRCMP(tri->key,trj->key);
  if(compval!=0)return compval;
  return strcmp(tri->key,trj->key);
}

/* ------------------ parse_lang ------------------------ */

int parse_lang(char *file, trdata **trinfoptr, int *ntrinfoptr){
  /*! \fn int parse_lang(char *file, trdata **trinfoptr, int *ntrinfoptr)
      \brief read a po file and put english/foreign language string pairs into
             the trinfo data structure
  */
  FILE *stream;
  trdata *trinfo_local;
  int ntrinfo_local;

  ntrinfo_local=*ntrinfoptr;
  trinfo_local=*trinfoptr;
  if(ntrinfo_local>0){
    int i;

    for(i=0;i<ntrinfo_local;i++){
      trdata *tri;
  
      tri = trinfo_local + i;
      FREEMEMORY(tri->key);
      FREEMEMORY(tri->value);
    }
    FREEMEMORY(trinfo_local);
  }

  ntrinfo_local=0;
  if(file==NULL)return 0;
  stream=fopen(file,"r");
  if(stream==NULL)return 0;

  while(!feof(stream)){
    char buffer[1000];
    char *buf;
    char *key;

    if(fgets(buffer,1000,stream)==NULL)break;
    buf=trim_front(buffer);
    trim(buf);
    if(strlen(buf)>=2&&strncmp(buf,"//",2)==0)continue;
    key = strstr(buf,"msgid");
    if(key!=NULL&&key==buf)ntrinfo_local++;
  }
  if(ntrinfo_local==0){
    fclose(stream);
    return 0;
  }

  NewMemory((void **)&trinfo_local,sizeof(trdata)*ntrinfo_local);
  *trinfoptr=trinfo_local;
  *ntrinfoptr=ntrinfo_local;

  ntrinfo_local=0;
  rewind(stream);
  while(!feof(stream)){
    char buffer[1000];
    trdata *tri;
    char *keybuf,*valbuf;
    char *buf;
    char *key, *value;
    int doit;

    if(fgets(buffer,1000,stream)==NULL)break;
    buf=trim_front(buffer);
    trim(buf);
    if(strlen(buf)>=2&&strncmp(buf,"//",2)==0)continue;
    key = strstr(buf,"msgid");
    if(key==NULL||key!=buf)continue;
    ntrinfo_local++;
    tri = trinfo_local + ntrinfo_local - 1;
    key = getstring(key+5);
    if(key==NULL){
      tri->key=NULL;
    }
    else{
      NewMemory((void **)&keybuf,strlen(key)+1);
      strcpy(keybuf,key);
      tri->key=keybuf;
    }

    for(doit=1;doit==1;){
      doit=0;
      if(fgets(buffer,1000,stream)==NULL)break;
      buf=trim_front(buffer);
      trim(buf);
      if(strlen(buf)>=2&&strncmp(buf,"//",2)==0)doit=1;
    }
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

  qsort(trinfo_local,ntrinfo_local,sizeof(trdata),compare_trdata);
  return 1;
}

/* ------------------ init_translate ------------------------ */

void init_translate(char *bindir, char *tr_name){
  /*! \fn void init_translate(char *bindir, char *tr_name)
      \brief initialize po language translation data structures
  */
  char *LANG;

  if(tr_name!=NULL){
    LANG=tr_name;
  }
  else{
    LANG=getenv("LANG");
  }

  tr_otherlang=0;
  if(LANG!=NULL&&strncmp(LANG,"en",2)!=0){
    FILE *stream;
    int lensmokebindir;
    char lang[256];

    lang[0]=tolower(LANG[0]);
    lang[1]=tolower(LANG[1]);
    lang[2]=0;
    lensmokebindir=strlen(bindir);
    FREEMEMORY(smokeview_lang);
    NewMemory((void **)&smokeview_lang,(unsigned int)(lensmokebindir+15+1));
    STRCPY(smokeview_lang,bindir);
    STRCAT(smokeview_lang,"smokeview_");
    STRCAT(smokeview_lang,lang);
    STRCAT(smokeview_lang,".po");
    
    stream=fopen(smokeview_lang,"r");
    if(stream!=NULL){
      fclose(stream);
      tr_otherlang=1;
    }
    tr_otherlang=parse_lang(smokeview_lang,&trinfo,&ntrinfo);
    if(tr_otherlang==1){
      PRINTF("Using translation file: %s",smokeview_lang);
      PRINTF("\n");
    }
    else{
      fprintf(stderr,"*** Warning: Failed to parse translation file: %s",smokeview_lang);
      fprintf(stderr,"\n");
      fprintf(stderr,"%s\n","Menus will be in English");
    }
  }
}

/* ------------------ translate ------------------------ */

char *translate(char *string){
  /*! \fn char *translate(char *string)
      \brief return the translation of string, return string if translation not found
  */
  int i, len, nchars_before=0, nchars_after=0;
  unsigned int nchars_in=0;
  char *string_before, *string_in, *string_out, *string_after;

  if(tr_otherlang==0)return string;


  len = strlen(string);

  // find leading non-alpha characters

  string_in = string;
  string_before = tr_string_before;
  for(i = 0; i<len; i++){
    char C,D;
    char c;

    c=string[i];
    C=toupper(c);
    D=toupper(string[i+1]);
    if((C>='A'&&C<='Z')||( (c=='1'||c=='2'||c=='3') && D=='D' ) ){
      nchars_before=i;
      string_in=string+i;
      if(nchars_before>0){
        string_before=tr_string_before;
        string_before[nchars_before]=0;
      }
      else{
        string_before=NULL;
      }
      break;
    }
    tr_string_before[i]=c;
  }

  // find trailing non-alpha characters

  string_after = string+len;
  for(i=len-1;i>=nchars_before;i--){
    char c;

    c=string[i];
    if((c>='a'&&c<='z')||(c>='A'&&c<='Z')){
      nchars_after=len-1-i;
      string_after=string+i+1;
      if(nchars_after>0){
        memcpy(tr_string_after,string_after,nchars_after);
        tr_string_after[nchars_after]=0;
        string_after=tr_string_after;
      }
      else{
        string_after=NULL;
      }
      break;
    }
  }
  nchars_in=len-nchars_before-nchars_after;
  if(nchars_in==0)return string;
  memcpy(tr_string_in,string_in,nchars_in);
  tr_string_in[nchars_in]=0;
  string_in=tr_string_in;

  // translate string_in

  {
    trdata *tr_out, tr_in;

    tr_in.key=string_in;

    tr_out = bsearch(&tr_in,trinfo,ntrinfo,sizeof(trdata),compare_trdata);
    if(tr_out==NULL||tr_out->value==NULL)return string;
    string_out=tr_out->value;
  }

  strcpy(tr_string,"");
  if(string_before!=NULL)strcat(tr_string,string_before);
  strcat(tr_string,string_out);
  if(string_after!=NULL)strcat(tr_string,string_after);
  return tr_string;
}
