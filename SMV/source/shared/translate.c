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
  *trinfoptr=trinfo;
  *ntrinfoptr=ntrinfo;

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
      printf(_("Using translation file: %s\n"),smokeview_lang);
    }
    else{
      printf(_("Failed to parse translation file: %s\n"),smokeview_lang);
      printf(_("Menus will be in English\n"));
    }
  }
}

/* ------------------ translate ------------------------ */

char *translate(char *string){
  /*! \fn char *translate(char *string)
      \brief return the translation of string, return string if translation not found
  */
  char c;
  int i, len, nchars_before=0, nchars_in=0, nchars_after=0;
  char *string_before, *string_in, *string_out, *string_after;

  if(tr_otherlang==0)return string;


  len = strlen(string);

  // find leading non-alpha characters

  for(i=0;i<len;i++){
    c=string[i];
    if((c>='a'&&c<='z')||(c>='A'&&c<='Z')){
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

  for(i=len-1;i>=nchars_before;i--){
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
