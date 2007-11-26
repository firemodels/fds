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
#include "MALLOC.h"

// svn revision character string
void init_boundbox0(void);
char readfds_revision[]="$Revision$";
int compare( const void *arg1, const void *arg2 );
int read_pass1(char *fdsfile, int recurse_level);

/* ------------------ readfds ------------------------ */

int readfds(char *fdsfile){
  
  FILE *streamfds;
#define LENBUFFER 10000
  char buffer[LENBUFFER];
  int in_assembly=0;

  if(read_pass1(fdsfile,0)!=0)return 1;

  init_bb();

  // pass 2

  streamfds=fopen(fdsfile,"r");
  if(streamfds==NULL){
    printf("The file: %s could not be opened\n",fdsfile);
    return 1;
  }
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
    if(match(buffer,"&INCL",5)==1){
      in_assembly=3;
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
        expand_assembly(buffer,0);
        break;
      case 3:  // &INCL line, don't write it out
        in_assembly=0;
        break;
    }
  }
  return 0;
}

/* ------------------ read_pass1 ------------------------ */

int read_pass1(char *fdsfile, int recurse_level){
  FILE *streamfds;
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
      assembly=create_assembly(buffer); // &BASM ID='....' ORIG=x,y,z
      continue;
    }
    if(match(buffer,"&EGRP",5)==1){ // &EASM /
      in_assembly=0;
      continue;
    }
    if(match(buffer,"&INCL",5)==1){ // &INCL FILE='kdkdkkdk' /
      char *file;

      file=getkeyid(buffer,"FILE");
      if(file!=NULL&&recurse_level<10){
        read_pass1(file,recurse_level+1);
      }
      continue;
    }
    if(in_assembly==0)continue;
   // using info in current buffer add to assembly data structures
    update_assembly(assembly,buffer);
  }
  fclose(streamfds);
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

/* ------------------ get_boundbox ------------------------ */

void get_boundbox(blockaiddata *assem,int recurse_level){
  fdsdata *thisline;
  int i;
  float *bb_min, *bb_max, *bb_dxyz;

  if(assem->bb_box_defined==1)return;

  if(recurse_level>MAXRECURSE){
    printf(" *** Fatal error:  recursion level must be less than %i\n",MAXRECURSE);
    return;
  }

  bb_min=assem->bb_min;
  bb_max=assem->bb_max;
  bb_dxyz=assem->bb_dxyz;

  for(thisline=assem->first_line->next;thisline->next!=NULL;thisline=thisline->next){
    float xb[6];

    if(thisline->line_after!=NULL&&thisline->line_before!=NULL){
      if(thisline->type==1&&thisline->is_obst==1){
        for(i=0;i<6;i++){
          xb[i]=thisline->xb[i]-assem->xyz0[i/2];
        }
      }
      else if(thisline->type==2){
        char linebuffer[1024];
        char *id2;
        blockaiddata *assem2;        
        float offset2[3], rotate2;
        int j;

        strcpy(linebuffer,thisline->line);

        offset2[0]=0.0;
        offset2[1]=0.0;
        offset2[2]=0.0;
        get_irvals(linebuffer, "XYZ", 3, NULL, offset2, NULL, NULL);

        rotate2=0.0;
        get_irvals(linebuffer, "ROTATE", 1, NULL, &rotate2, NULL, NULL);

        id2=getkeyid(linebuffer,"GRP_ID");
        assem2=get_assembly(id2);

        get_boundbox(assem2,recurse_level+1);

        xb[0]=assem2->bb_min[0];
        xb[1]=assem2->bb_max[0];
        xb[2]=assem2->bb_min[1];
        xb[3]=assem2->bb_max[1];
        xb[4]=assem2->bb_min[2];
        xb[5]=assem2->bb_max[2];

        rotatexy(xb,xb+2,  rotate2, assem2->bb_dxyz);
        rotatexy(xb+1,xb+3,rotate2, assem2->bb_dxyz);
        reorder(xb);
        reorder(xb+2);
        reorder(xb+4);
        for(j=0;j<6;j++){
          xb[j]+=offset2[j/2];
        }
      }
      else{
        continue;
      }

      if(xb[0]<bb_min[0])bb_min[0]=xb[0];
      if(xb[1]>bb_max[0])bb_max[0]=xb[1];

      if(xb[2]<bb_min[1])bb_min[1]=xb[2];
      if(xb[3]>bb_max[1])bb_max[1]=xb[3];

      if(xb[4]<bb_min[2])bb_min[2]=xb[4];
      if(xb[5]>bb_max[2])bb_max[2]=xb[5];

      bb_dxyz[0]=bb_max[0]-bb_min[0];
      bb_dxyz[1]=bb_max[1]-bb_min[1];
      bb_dxyz[2]=bb_max[2]-bb_min[2];
    }
  }
}

/* ------------------ init_boundbox0 ------------------------ */

void init_boundbox0(void){
  fdsdata *thisline;
  blockaiddata *assem;

  for(assem=blockaid_first->next;assem->next!=NULL;assem=assem->next){
    float xb[6];
    int do_this_assm;
    float *bb_min, *bb_max, *bb_dxyz;

    if(assem->bb_box_defined==1)continue;

    bb_min=assem->bb_min;
    bb_max=assem->bb_max;
    bb_dxyz=assem->bb_dxyz;

    do_this_assm=1;
    for(thisline=assem->first_line->next;thisline->next!=NULL;thisline=thisline->next){

      if(thisline->type==1&&thisline->is_obst==1){
      }
      else if(thisline->type==2){
        do_this_assm=0;
        break;
      }
    }
    if(do_this_assm==0)continue;

    for(thisline=assem->first_line->next;thisline->next!=NULL;thisline=thisline->next){
      int i;

      if(thisline->type!=1||thisline->is_obst!=1)continue;

      for(i=0;i<6;i++){
        xb[i]=thisline->xb[i]-assem->xyz0[i/2];
      }

      if(xb[0]<bb_min[0])bb_min[0]=xb[0];
      if(xb[1]>bb_max[0])bb_max[0]=xb[1];

      if(xb[2]<bb_min[1])bb_min[1]=xb[2];
      if(xb[3]>bb_max[1])bb_max[1]=xb[3];

      if(xb[4]<bb_min[2])bb_min[2]=xb[4];
      if(xb[5]>bb_max[2])bb_max[2]=xb[5];
    }
    bb_dxyz[0]=bb_max[0]-bb_min[0];
    bb_dxyz[1]=bb_max[1]-bb_min[1];
    bb_dxyz[2]=bb_max[2]-bb_min[2];
    assem->bb_box_defined=1;
  }
}
/* ------------------ expand_assembly ------------------------ */

void expand_assembly(char *buffer, int recurse_level){
  float offset[3], rotate;
  char *id;
  char blank[100];
  blockaiddata *assem;
  fdsdata *thisline;
  char charxb[32];
  int i,j;

  if(recurse_level>MAXRECURSE){
    printf(" *** Fatal error:  recursion level must be less than %i\n",MAXRECURSE);
    return;
  }

  offset[0]=0.0;
  offset[1]=0.0;
  offset[2]=0.0;
  get_irvals(buffer, "XYZ", 3, NULL, offset, NULL, NULL);

  rotate=0.0;
  get_irvals(buffer, "ROTATE", 1, NULL, &rotate, NULL, NULL);

  id=getkeyid(buffer,"GRP_ID");
  assem=get_assembly(id);
  if(assem==NULL){
    printf(" **** warning ****\n");
    printf("      The blockage assembly, %s, is not defined\n\n",id);
    return;
  }

  assemblylist[recurse_level]=assem;
  offset_rotate[4*recurse_level  ] =offset[0];
  offset_rotate[4*recurse_level+1]=offset[1];
  offset_rotate[4*recurse_level+2]=offset[2];
  offset_rotate[4*recurse_level+3]=rotate ;


  if(recurse_level==0){
    printf("\n MAJOR GROUP: %s offset=%f,%f,%f rotate=%f\n",
      assem->id,offset[0],offset[1],offset[2],rotate);
  }
  else{
    int lenspace;

    for(i=0;i<recurse_level;i++){
      if(strcmp(id,assemblylist[i]->id)==0){
        printf(" **** warning ****\n");
        printf("      Block defintions with Id's: \n");
        for(j=0;j<recurse_level;j++){
          printf(" %s,",assemblylist[j]->id);
        }
        printf(" %s\n",assemblylist[recurse_level]->id);
        printf(" are defined circularly.  Their expansion is halted\n");
        printf(" **** warning ****\n");
        return;
      }
    }
    
    lenspace = recurse_level;
    if(lenspace>4)lenspace=4;
    strcpy(blank,"");
    for(j=0;j<lenspace;j++){
      strcat(blank,"  ");
    }
    printf("\n%s MINOR GROUP: %s offset=%f,%f,%f rotate=%f\n",
      blank,assem->id,offset[0],offset[1],offset[2],rotate);
  }

  for(thisline=assem->first_line->next;thisline->next!=NULL;thisline=thisline->next){
    float xb[6];

    if(thisline->line_after!=NULL&&thisline->line_before!=NULL){
      if(thisline->type==1){
        float *xyz, *rotate;

        if(thisline->blockaid!=NULL&&thisline->blockaid->nkeywordlist>0){
          char line_before2[10000];

          subst_keys(line_before2,thisline,thisline->line_before);
          printf("%s ",line_before2);
        }
        else{
          printf("%s ",thisline->line_before);
        }
        for(i=0;i<6;i++){
          xb[i]=thisline->xb[i]-assem->xyz0[i/2];
        }

        for(i=recurse_level;i>=0;i--){
          blockaiddata *assemi;

          assemi = assemblylist[i];
          xyz = offset_rotate+4*i;
          rotate = offset_rotate+4*i+3;
          rotatexy(xb,xb+2,  rotate[0],assemi->bb_dxyz);
          rotatexy(xb+1,xb+3,rotate[0],assemi->bb_dxyz);
          for(j=0;j<6;j++){
            xb[j]+=xyz[j/2];
          }
          reorder(xb);
          reorder(xb+2);
          reorder(xb+4);

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
        if(thisline->blockaid!=NULL&&thisline->blockaid->nkeywordlist>0){
          char line_after2[10000];

          subst_keys(line_after2,thisline,thisline->line_after);
          printf("%s ",line_after2);
        }
        else{
          printf("%s\n",thisline->line_after);
        }

      }
      else if(thisline->type==2){
        char linebuffer[1024];

        strcpy(linebuffer,thisline->line);
        expand_assembly(linebuffer,recurse_level+1);
      }
    }
  }

}

/* ------------------ reorder ------------------------ */

void reorder(float *xy){
  float xymin, xymax;

  xymin = xy[0];
  if(xy[1]<xymin)xymin=xy[1];

  xymax = xy[0];
  if(xy[1]>xymax)xymax=xy[1];

  xy[0]=xymin;
  xy[1]=xymax;
}

/* ------------------ subst ------------------------ */

void subst(char *line, char *keybeg,char *keyend,char *val){
  char lineend[10000];
  size_t len_from, len_to;

  if(keybeg==NULL||keyend==NULL)return;

  len_from=keyend-keybeg;
  if(len_from<1)return;

  len_to=strlen(val);
  if(len_to<0)return;

  strcpy(lineend,keyend);
  if(len_to>0)strncpy(keybeg,val,len_to);
  strcat(keybeg+len_to,lineend);

}

/* ------------------ get_keyend ------------------------ */

char *get_keyend(char *keybeg){
  char *key;

  key=keybeg;
  for(key=keybeg;;key++){
    if(*key=='\0')return NULL;
    if(*key=' '||*key==','||*key=='/')return key;
    if(key-keybeg>1000)return NULL;
  }
}

/* ------------------ get_val ------------------------ */

char *get_val(char *keybeg, char *keyend){
  char key2[100],*key;
  char *c2;

  if(keybeg==NULL||keyend==NULL)return NULL;

  c2=key2;
  for(key=keybeg;key<keyend;key++){
    *c2++=*key;
  }
  *c2='\0';
  //  get value associated with key2

  return NULL;
}

/* ------------------ subst_keys ------------------------ */

void subst_keys(char *line_after,fdsdata *thisline,char *line_before){
  char line[10000];
  char *keybeg, *keyend, *val;

  strcpy(line,line_before);
  keybeg=strstr(line,"#");
  keyend=get_keyend(keybeg);

  val=get_val(keybeg,keyend);
  while(val!=NULL){
    subst(line,keybeg,keyend,val);

    keybeg=strstr(line,"#");
    keyend=get_keyend(keybeg);

    val=get_val(keybeg,keyend);
  }
  strcpy(line_after,line);
  return;
}
  

/* ------------------ get_keywords ------------------------ */

void get_keywords(blockaiddata *blockaidi, char *buffer){
  char buffer2[1024], *buffptr;
  char *key, *endkey, *val, *endval;
  int nkeys;

  strcpy(buffer2,buffer);
  buffptr=buffer2;
  nkeys=0;
  for(key=strstr(buffptr,"#");key!=NULL;){
    nkeys++;
    buffptr=key+1;
    key=strstr(buffptr,"#");
  }
  blockaidi->nkeywordlist=0;
  if(nkeys==0)return;

  NewMemory((void **)&blockaidi->keyword_list,nkeys*sizeof(char *));
  NewMemory((void **)&blockaidi->val_list,nkeys*sizeof(char *));
  nkeys=0;
  buffptr=buffer2;

  for(key=strstr(buffptr,"#");key!=NULL;){
    size_t lenkey, lenval;
    char *ckey, *cval;

    endkey=strstr(key+1,"=");
    if(endkey==NULL)break;

    val=strstr(endkey+1,"'");
    if(val==NULL)break;
    val++;

    endval=strstr(val,"'");
    if(endval==NULL)break;
    *endkey='\0';
    *endval='\0';
    lenkey=strlen(key);
    lenval=strlen(val);
    if(lenkey<1||lenval<1)continue;

    NewMemory((void **)&ckey,lenkey+1);
    NewMemory((void **)&cval,lenval+1);
    strcpy(ckey,key);
    strcpy(cval,val);
    blockaidi->keyword_list[nkeys]=ckey;
    blockaidi->val_list[nkeys]=cval;
    nkeys++;
    buffptr=endval+1;
    key=strstr(buffptr,"#");
  }
  blockaidi->nkeywordlist=nkeys;
  return;
}

/* ------------------ create_assembly ------------------------ */

blockaiddata *create_assembly(char *buffer){
  blockaiddata *blockaidi, *bprev, *bnext;
  float *orig;
  char *id;
  fdsdata *first_line, *last_line;

  NewMemory((void **)&blockaidi,sizeof(blockaiddata));
  bprev=blockaid_first;
  bnext=blockaid_first->next;
  blockaidi->prev=bprev;
  blockaidi->next=bnext;
  bprev->next=blockaidi;
  bnext->prev=blockaidi;
  get_keywords(blockaidi,buffer);

  orig=blockaidi->xyz0;

  if(get_irvals(buffer, "ORIG", 3, NULL, orig, NULL, NULL)==3){
    // first 3 positions of orig are defined in get_irvals
    orig[3]=1.0;
  }
  else{
    orig[0]=0.0;
    orig[1]=0.0;
    orig[2]=0.0;
    orig[3]=1.0;
  }
  id=getkeyid(buffer,"GRP_ID");
  if(id!=NULL){
    NewMemory((void **)&blockaidi->id,strlen(id)+1);
    strcpy(blockaidi->id,id);
  }
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

void update_assembly(blockaiddata *assembly,char *buffer){
  fdsdata *prev, *next, *thisfds;
  size_t len;
  int is_obst=0;

  next=assembly->last_line;
  prev=next->prev;

  NewMemory((void **)&thisfds,sizeof(fdsdata));

  thisfds->line=NULL;
  thisfds->is_obst=0;
  thisfds->blockaid=NULL;
  if(match(buffer,"&OBST",5)==1||
    match(buffer,"&HOLE",5)==1||
    match(buffer,"&VENT",5)==1){
    thisfds->type=1;
    if(match(buffer,"&OBST",5)==1){
      is_obst=1;
      thisfds->is_obst=1;
    }
  }
  else if(match(buffer,"&GRP",4)==1){
    thisfds->type=2;
    if(strstr(buffer,"#")!=NULL){
      NewMemory((void **)&thisfds->blockaid,sizeof(blockaiddata ));
      if(thisfds->blockaid!=NULL){
        get_keywords(thisfds->blockaid, buffer);
      }
    }
    else{
      thisfds->blockaid=NULL;
    }
  }
  else{
    thisfds->type=0;
  }
  if(buffer!=NULL){
    len=strlen(buffer);
    NewMemory((void **)&thisfds->line,len+1);
    NewMemory((void **)&thisfds->linecopy,len+1);
    strcpy(thisfds->line,buffer);
    strcpy(thisfds->linecopy,buffer);
    if(thisfds->type==1){
      float *orig,oorig[3], *xb;
        
      get_irvals(buffer, "XB", 6, NULL, thisfds->xb,&thisfds->ibeg,&thisfds->iend);
      if(is_obst==1){

        xb = thisfds->xb;
        if(assembly->xyz0[3]<0.0){
          float *orig;

          orig = assembly->xyz0;
          if(xb[0]<orig[0])orig[0]=xb[0];
          if(xb[2]<orig[1])orig[1]=xb[2];
          if(xb[4]<orig[2])orig[2]=xb[4];
        }
        else{
          oorig[0]=0.0;
          oorig[1]=0.0;
          oorig[2]=0.0;
          orig=oorig;
        }
      }
    }
    else if(thisfds->type==2){
      get_irvals(buffer, "XYZ", 3, NULL, thisfds->xb,&thisfds->ibeg,&thisfds->iend);
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

/* ------------------ remove_assemblyb ------------------------ */

void remove_assembly(blockaiddata *assemb){
}

/* ------------------ init_bb ------------------------ */

void init_bb(void){
  blockaiddata *assm;

  for(assm=blockaid_first->next;assm->next!=NULL;assm=assm->next){
    assm->bb_min[0]=MAXPOS;
    assm->bb_max[0]=MINPOS;

    assm->bb_min[1]=MAXPOS;
    assm->bb_max[1]=MINPOS;

    assm->bb_min[2]=MAXPOS;
    assm->bb_max[2]=MINPOS;

    assm->bb_box_defined=0;
  }
  init_boundbox0();
  for(assm=blockaid_first->next;assm->next!=NULL;assm=assm->next){
    float *x1, *x2, *dx;

    get_boundbox(assm,0);

    x1 = assm->bb_min;
    x2 = assm->bb_max;
    dx = assm->bb_dxyz;

    printf("assm=%s\n",assm->id);
    printf("  bounds=(%f,%f), (%f,%f), (%f,%f)\n",x1[0],x2[0],x1[1],x2[1],x1[2],x2[2]);
  }
}

/* ------------------ get_assembly ------------------------ */

blockaiddata *get_assembly(char *id){
  blockaiddata *assm;

  if(id==NULL)return NULL;
  for(assm=blockaid_first->next;assm->next!=NULL;assm=assm->next){
    if(strcmp(assm->id,id)==0)return assm;
  }
  return NULL;
}


/* ------------------ init_assemdata ------------------------ */

void init_assemdata(char *id, float *orig, blockaiddata *prev, blockaiddata *next){
  blockaiddata *newassem;
  fdsdata *fl, *ll;

  NewMemory((void **)&newassem,sizeof(blockaiddata));
  strcpy(newassem->id,id);
  newassem->xyz0[0]=orig[0];
  newassem->xyz0[1]=orig[1];
  newassem->xyz0[2]=orig[2];
  newassem->prev=prev;
  newassem->next=next;

  prev->next=newassem;
  next->prev=newassem;
  NewMemory((void **)&newassem->first_line,sizeof(fdsdata));

  NewMemory((void **)&newassem->last_line,sizeof(fdsdata));
  fl = newassem->first_line;
  ll = newassem->last_line;

  fl->line=NULL;
  fl->prev=NULL;
  fl->next=ll;

  ll->line=NULL;
  ll->prev=fl;
  ll->next=NULL;
}
