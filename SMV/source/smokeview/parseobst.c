// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include "flowfiles.h"
#include "MALLOC.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"

// svn revision character string
char parseobst_revision[]="$Revision$";

/* ------------------ getlabels ------------------------ */

void getlabels(const char *filein){

  FILE *stream_in;
  char buffer[1000];
  int fdsobstcount=0;
  int i,j;
  char *obstlabel;
  mesh *meshi;
  blockagedata *bc;
  int id;
  size_t lenlabel;
  char **obstlabels=NULL;
  int nobstlabels=0;

  if(filein==NULL)return;
  stream_in = fopen(filein,"r");
  if(stream_in==NULL)return;

  while(!feof(stream_in)){
    if(fgets(buffer,1000,stream_in)==NULL)break;

    if(STRSTR(buffer,"&OBST")==NULL)continue;
    fdsobstcount++;
  }
  nobstlabels=fdsobstcount;
  if(nobstlabels>0)NewMemory((void **)&obstlabels,nobstlabels*sizeof(char *));
  for(i=0;i<nobstlabels;i++){
    obstlabels[i]=NULL;
  }
  rewind(stream_in);
  fdsobstcount=0;
  while(!feof(stream_in)){
    if(fgets(buffer,1000,stream_in)==NULL)break;

    if(STRSTR(buffer,"&OBST")==NULL)continue;
    fdsobstcount++;
    while((obstlabel=strstr(buffer,"/"))==NULL){
      fgets(buffer,1000,stream_in);
    }
    obstlabel++;
    lenlabel=strlen(obstlabel);
    obstlabel=trim_front(obstlabel);
    trim(obstlabel);
    lenlabel=strlen(obstlabel);
    if(lenlabel>0){
      NewMemory((void **)&obstlabels[fdsobstcount-1],(unsigned int)(lenlabel+1));
      strcpy(obstlabels[fdsobstcount-1],obstlabel);
    }
  }
  fclose(stream_in);

  for(i=0;i<nmeshes;i++){
    meshi = meshinfo + i;
    for(j=0;j<meshi->nbptrs;j++){
      bc = meshi->blockageinfoptrs[j];
      id = bc->id-1;
      if(id>=0&&id<nobstlabels){
        if(obstlabels[id]!=NULL){
          lenlabel=strlen(obstlabels[id]);
          ResizeMemory((void **)&bc->label,(unsigned int)(lenlabel+1));
          strcpy(bc->label,obstlabels[id]);
        }
      }
    }
  }
  for(i=0;i<nobstlabels;i++){
    FREEMEMORY(obstlabels[i]);
  }
  FREEMEMORY(obstlabels);
}
