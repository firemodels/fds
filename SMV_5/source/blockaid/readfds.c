// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "svn_revision.h"
#include "blockaid.h"

// svn revision character string
char readfds_revision[]="$Revision$";

/* ------------------ readfds ------------------------ */

int readfds(char *fdsfile){
  
  FILE *streamfds;
#define LENBUFFER 10000
  char buffer[LENBUFFER];
  int in_assembly=0;

  streamfds=fopen(fdsfile,"r");
  if(streamfds==NULL){
    printf("The file: %s could not be opened\n",fdsfile);
    return 1;
  }

// pass 1

  while(!feof(streamfds)){
    blockaiddata *assembly;

    if(get_fds_line(streamfds, buffer, LENBUFFER)==-1)break;
    trim(buffer);

    if(match(buffer,"&BGRP",5)==1){
      in_assembly=1;
      assembly=create_assembly(buffer); // BASM ID='....' ORIG=x,y,z
      continue;
    }
    if(match(buffer,"&EGRP",5)==1){
      in_assembly=0;
      continue;
    }
    if(in_assembly==0)continue;

    // using info in current buffer add to assembly data structures
    update_assembly(assembly,buffer);

  }

  // pass 2

  rewind(streamfds);
  in_assembly=0;
  while(!feof(streamfds)){
    if(get_fds_line(streamfds, buffer, LENBUFFER)==-1)break;
    trim(buffer);

    if(match(buffer,"&BGRP",5)==1){
      in_assembly=1;
      continue;
    }
    if(match(buffer,"&EGRP",5)==1){
      in_assembly=0;
      continue;
    }
    if(in_assembly==0&&match(buffer,"&GRP",4)==1){
      in_assembly=2;
    }
    switch (in_assembly){
      case 0:  // regular line, output it
        printf("%s\n",buffer);
        break;
      case 1:  // inside an assembly defn, skip
        break;
      case 2:  // assembly line, apply translation and rotation then output assembly lines
        in_assembly=0;
        expand_assembly(buffer,1);
        break;
    }
  }
  return 0;
}

/* ------------------ get_fds_line ------------------------ */

int get_fds_line(FILE *stream, char *fdsbuffer, unsigned int len_fdsbuffer){
  int copyback=0;
  size_t lenbuffer2;
  char buffer[LENBUFFER], buffer2[LENBUFFER];
  int is_command=0;

  copyback=0;
  if(fgets(buffer,LENBUFFER,stream)==NULL||strlen(buffer)>len_fdsbuffer)return -1;
  strcpy(buffer2,buffer);
  lenbuffer2=0;
  if(buffer[0]=='&')is_command=1;
  while(is_command==1&&strstr(buffer,"/")==NULL){
    if(fgets(buffer,LENBUFFER,stream)==NULL)return -1;
    lenbuffer2+=strlen(buffer);
    if(lenbuffer2>len_fdsbuffer||lenbuffer2>LENBUFFER)return -1;
    strcat(buffer2,buffer);
    copyback=1;
  }
  if(copyback==1){
    strcpy(fdsbuffer,buffer2);
  }
  else{
    strcpy(fdsbuffer,buffer);
  }
  return (int)strlen(fdsbuffer);
}

/* ------------------ expand_assembly ------------------------ */

void expand_assembly(char *buffer, int first_time){
  char buffer2[1000],buffer3[1000];
  float offset[3], rotate;
  float offset2[3], rotate2;
  char *id,*id2;
  blockaiddata *assem,*assemi;
  fdsdata *thisline;
  char charxb[32], charoffset[32], charrotate[32];
  int i;

  offset[0]=0.0;
  offset[1]=0.0;
  offset[2]=0.0;
  get_irvals(buffer, "OFFSET", 3, NULL, offset, NULL, NULL);

  rotate=0.0;
  get_irvals(buffer, "ROTATE", 1, NULL, &rotate, NULL, NULL);

  id=getkeyid(buffer,"GRP_ID");
  assem=get_assembly(id);

  if(first_time==1){
    for(assemi=blockaid_first->next;assemi->next!=NULL;assemi=assemi->next){
      assemi->in_use=0;
    }
    assem->in_use=1;
    printf("\n MAJOR GROUP: %s offset=%f,%f,%f rotate=%f\n\n",assem->id,offset[0],offset[1],offset[2],rotate);
  }
  else{
    //if(assem->in_use==1){
    //  printf(" *** fatal error: A &GRP may not contained within another &GRP with the same ID\n");
    //  //exit(1);
   // }
    printf("\n    MINOR GROUP: %s offset=%f,%f,%f rotate=%f\n\n",assem->id,offset[0],offset[1],offset[2],rotate);
    assem->in_use=1;
  }

  for(thisline=assem->first_line->next;thisline->next!=NULL;thisline=thisline->next){
    float xb[6];

    if(thisline->line_after!=NULL&&thisline->line_before!=NULL){
      if(thisline->type==1){
        printf("%s ",thisline->line_before);
        for(i=0;i<6;i++){
          xb[i]=thisline->xb[i];
        }
        rotatexy(xb,xb+2,assem->orig,rotate);
        rotatexy(xb+1,xb+3,assem->orig,rotate);
        for(i=0;i<6;i++){
          xb[i]+=offset[i/2];
        }
        for(i=0;i<6;i++){
          sprintf(charxb,"%f",xb[i]);
          trimzeros(charxb);
          if(i==5){
            printf("%s",charxb);
          }
          else{
            printf("%s,",charxb);
          }
        }
        printf("%s\n",thisline->line_after);
      }
      else if(thisline->type==2){
        strcpy(buffer2,thisline->line);
        offset2[0]=0.0;
        offset2[1]=0.0;
        offset2[2]=0.0;
        get_irvals(buffer2, "OFFSET", 3, NULL, offset2, NULL, NULL);
        offset2[0]+=offset[0];
        offset2[1]+=offset[1];
        offset2[2]+=offset[2];

        rotate2=0.0;
        get_irvals(buffer2, "ROTATE", 1, NULL, &rotate2, NULL, NULL);
        rotate2+=rotate;

        id2=getkeyid(buffer2,"GRP_ID");

        strcpy(buffer3,"&GRP GRP_ID='");
        strcat(buffer3,id2);
        strcat(buffer3,"' OFFSET=");

        sprintf(charoffset,"%f",offset2[0]);
        trimzeros(charoffset);
        strcat(buffer3,charoffset);
        strcat(buffer3,", ");

        sprintf(charoffset,"%f",offset2[1]);
        trimzeros(charoffset);
        strcat(buffer3,charoffset);
        strcat(buffer3,", ");

        sprintf(charoffset,"%f",offset2[2]);
        trimzeros(charoffset);
        strcat(buffer3,charoffset);
        strcat(buffer3,", ROTATE=");

        sprintf(charrotate,"%f",rotate2);
        trimzeros(charrotate);
        strcat(buffer3,charrotate);
        strcat(buffer3," /");


        expand_assembly(buffer3,0);
      }
    }
  }

}

/* ------------------ create_assembly ------------------------ */

blockaiddata *create_assembly(char *buffer){
  blockaiddata *blockaidi, *bprev, *bnext;
  float *orig;
  char *id;
  fdsdata *first_line, *last_line;

  blockaidi=malloc(sizeof(blockaiddata));
  bprev=blockaid_first;
  bnext=blockaid_first->next;
  blockaidi->prev=bprev;
  blockaidi->next=bnext;
  bprev->next=blockaidi;
  bnext->prev=blockaidi;

  orig=blockaidi->orig;
  orig[0]=0.0;
  orig[1]=0.0;
  orig[2]=0.0;

  get_irvals(buffer, "ORIG", 3, NULL, orig, NULL, NULL);
  id=getkeyid(buffer,"GRP_ID");
  if(id!=NULL){
    blockaidi->id=malloc(strlen(id)+1);
    strcpy(blockaidi->id,id);
  }
/*
typedef struct _fdsdata {
  char *line;
  struct _fdsdata *prev, *next;
} fdsdata;
*/
  blockaidi->first_line=&blockaidi->f_line;
  blockaidi->last_line=&blockaidi->l_line;
  first_line=blockaidi->first_line;
  last_line=blockaidi->last_line;
  first_line->line=NULL;
  first_line->next=last_line;
  first_line->prev=NULL;
  last_line->line=NULL;
  last_line->next=NULL;
  last_line->prev=first_line;


  return blockaidi;
}

/* ------------------ update_assembly ------------------------ */

/*
typedef struct _fdsdata {
  char *line;
  struct _fdsdata *prev, *next;
} fdsdata;

typedef struct _blockaiddata {
  char *id;
  float orig[3];
  struct _fdsdata *first_line, *last_line;
  struct _fdsdata f_line, l_line;
  struct _blockaiddata *prev, *next;
} blockaiddata;
*/

void update_assembly(blockaiddata *assembly,char *buffer){
  fdsdata *prev, *next, *thisfds;
  size_t len;

  next=assembly->last_line;
  prev=next->prev;

  thisfds=malloc(sizeof(fdsdata));
  thisfds->line=NULL;
  if(match(buffer,"&OBST",5)==1||
    match(buffer,"&HOLE",5)==1||
    match(buffer,"&VENT",5)==1){
    thisfds->type=1;
  }
  else if(match(buffer,"&GRP",4)==1){
    thisfds->type=2;
  }
  else{
    thisfds->type=0;
  }
  if(buffer!=NULL){
    len=strlen(buffer);
    thisfds->line=malloc(len+1);
    thisfds->linecopy=malloc(len+1);
    strcpy(thisfds->line,buffer);
    strcpy(thisfds->linecopy,buffer);
    if(thisfds->type==1){
      get_irvals(buffer, "XB", 6, NULL, thisfds->xb,&thisfds->ibeg,&thisfds->iend);
    }
    else if(thisfds->type==2){
      get_irvals(buffer, "OFFSET", 3, NULL, thisfds->xb,&thisfds->ibeg,&thisfds->iend);
    }

    if(thisfds->type!=0&&(thisfds->ibeg>=0&&thisfds->iend>=0)){
      thisfds->linecopy[thisfds->ibeg]=0;
      thisfds->line_before=thisfds->linecopy;
      thisfds->line_after=thisfds->linecopy+ thisfds->iend;
    }
    else{
      thisfds->line_before=NULL;
      thisfds->line_after=NULL;
    }


  }
  prev->next=thisfds;
  thisfds->prev=prev;
  thisfds->next=next;
  next->prev=thisfds;
}

void remove_assembly(blockaiddata *assemb){
}

/* ------------------ get_assembly ------------------------ */

blockaiddata *get_assembly(char *id){
  blockaiddata *assm;

  for(assm=blockaid_first->next;assm->next!=NULL;assm=assm->next){
    if(strcmp(assm->id,id)==0)return assm;
  }
  return NULL;
}

/* ------------------ init_assemdata ------------------------ */

void init_assemdata(char *id, float *orig, blockaiddata *prev, blockaiddata *next){
  blockaiddata *newassem;
  fdsdata *fl, *ll;

  newassem = malloc(sizeof(blockaiddata));
  strcpy(newassem->id,id);
  newassem->orig[0]=orig[0];
  newassem->orig[1]=orig[1];
  newassem->orig[2]=orig[2];
  newassem->prev=prev;
  newassem->next=next;
  prev->next=newassem;
  next->prev=newassem;
  newassem->first_line=malloc(sizeof(fdsdata));
  newassem->last_line=malloc(sizeof(fdsdata));
  fl = newassem->first_line;
  ll = newassem->last_line;

  fl->line=NULL;
  fl->prev=NULL;
  fl->next=ll;

  ll->line=NULL;
  ll->prev=fl;
  ll->next=NULL;
}
