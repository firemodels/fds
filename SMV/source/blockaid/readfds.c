// $Date: 2012-09-13 12:38:42 -0400 (Thu, 13 Sep 2012) $ 
// $Revision: 12611 $
// $Author: gforney $

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
char readfds_revision[]="$Revision: 12611 $";
int compare( const void *arg1, const void *arg2 );
int read_pass1(char *fdsfile, int recurse_level);
void expand_shell(FILE *stream_out, char *buffer);

/* ------------------ readfds ------------------------ */

int readfds(char *in_file_base){
  
  char in_file[256], out_file[56];
  char file[256];
  char *ext;

  FILE *stream_in, *stream_out;
#define LENBUFFER 10000
  char buffer[LENBUFFER], repbuffer[LENBUFFER];
  int irep, nrep;
  int in_group=0;

  strcpy(file,in_file_base);
  trim(file);
  ext=strrchr(file,'.');
  if(ext==NULL||ext==file){
    ext=NULL;
  }
  else{
    *ext='\0';
  }

  if(ext==NULL){
    strcpy(in_file,file);
    strcat(in_file,".fof");
  }
  else{
    strcpy(in_file,in_file_base);
  }

  strcpy(out_file,file);
  strcat(out_file,".fds");

  if(strcmp(in_file,out_file)==0){
    printf("ERROR: The input and output files: %s are the same.\n",in_file);
    return 1;
  }

  if(getfileinfo(in_file,NULL,NULL)!=0){
    printf("ERROR: The input file file %s does not exist\n",in_file);
    return 1;
  }

  stream_in=fopen(in_file,"r");
  if(stream_in==NULL){
    printf("ERROR: The input file %s does not exist or can not be opened for input\n",in_file);
    return 1;
  }
  fclose(stream_in);

  stream_out=fopen(out_file,"r");
  if(stream_out!=NULL){
    if(force_write!=1){
      printf("ERROR: The output file %s exists, use the -f option to force an overwrite.\n",out_file);
      fclose(stream_out);
      return 1;
    }
    fclose(stream_out);
  }

  if(read_pass1(in_file,0)!=0)return 1;

  init_bb();

  // pass 2

  stream_in=fopen(in_file,"r");
  if(stream_in==NULL){
    printf("ERROR: The input file %s could not be opened.\n",in_file);
    return 1;
  }
  stream_out=fopen(out_file,"w");
  if(stream_out==NULL){
    printf("ERROR: The ouput file %s could not be opened for output.\n",out_file);
    return 1;
  }

  in_group=0;
  irep=0;
  nrep=0;
  while((nrep>0&&irep<nrep)||!feof(stream_in)){
    if(get_fds_line(stream_in, buffer, repbuffer, LENBUFFER, &irep, &nrep)==-1)break;
    if(nrep!=0&&irep>=nrep){
      irep=0;
      nrep=0;
    }
    trim(buffer);

    if(match(buffer,"&BGRP",5)==1){
      in_group=1;
      continue;
    }
    if(match(buffer,"&EGRP",5)==1){
      in_group=0;
      continue;
    }
    if(match(buffer,"&INCL",5)==1){
      in_group=3;
      continue;
    }
    if(in_group==0&&match(buffer,"&GRP",4)==1){
      in_group=2;
    }
    if(in_group==0&&match(buffer,"&SHELL",6)==1){
      in_group=4;
    }
    switch (in_group){
      case 0:  // regular line, output it
        fprintf(stream_out,"%s\n",buffer);
        break;
      case 1:  // inside an group defn, skip
        break;
      case 2:  // group line, apply translation and rotation then output group lines
        in_group=0;
        get_keywords(blockaid_first, buffer);

        nkeyvalstack=1;
        keyvalstack[nkeyvalstack-1].keyword_list=blockaid_first->keyword_list;
        keyvalstack[nkeyvalstack-1].val_list=blockaid_first->val_list;
        keyvalstack[nkeyvalstack-1].nkeywordlist=blockaid_first->nkeywordlist;

        expand_group(stream_out,buffer,0);
        break;
      case 3:  // &INCL line, don't write it out
        in_group=0;
        break;
      case 4:  // &SHELL
        expand_shell(stream_out,buffer);
        in_group=0;
        break;
    }
  }
  return 0;
}

/* ------------------ read_pass1 ------------------------ */

int read_pass1(char *in_file, int recurse_level){
  FILE *stream_in;
  char buffer[LENBUFFER], repbuffer[LENBUFFER];
  int irep, nrep;
  int in_group=0;

  stream_in=fopen(in_file,"r");
  if(stream_in==NULL){
    if(libdir!=NULL){
      char libfile[1024];

      strcpy(libfile,libdir);
      strcat(libfile,in_file);
      stream_in=fopen(libfile,"r");
      if(stream_in==NULL){
        printf("ERROR: The input files %s or %s could not be opened for input.\n",in_file,libfile);
        return 1;
      }
    }
    else{
      printf("ERROR: The input file %s could not be opened for input.\n",in_file);
      return 1;
    }
  }

// pass 1
  irep=0;
  nrep=0;
  while((nrep>0&&irep<nrep)||!feof(stream_in)){
    blockaiddata *group;

    if(get_fds_line(stream_in, buffer, repbuffer, LENBUFFER, &irep, &nrep)==-1)break;
    if(nrep!=0&&irep>=nrep){
      irep=0;
      nrep=0;
    }
    trim(buffer);

    if(match(buffer,"&BGRP",5)==1){
      in_group=1;
      group=create_group(buffer); // &BGRP ID='....' ORIGIN=x,y,z
      if(group==NULL){
        printf("FATAL ERROR: failed to create group for the following:\n");
        printf("%s\n",buffer);
        exit(1);
      }
      continue;
    }
    if(match(buffer,"&EGRP",5)==1){ // &EGRP /
      in_group=0;
      continue;
    }
    if(match(buffer,"&INCL",5)==1){ // &INCL FILE='kdkdkkdk' /
      char *file;

      file=get_keyid(buffer,"FILE");
      if(file!=NULL){
        if(recurse_level<10){
          read_pass1(file,recurse_level+1);
        }
        else{
          printf("ERROR: include file %s nested to deeply\n",file);
        }
      }
      continue;
    }
    if(in_group==0)continue;
   // using info in current buffer add to group data structures
    update_group(group,buffer);
  }
  fclose(stream_in);
  return 0;
}

/* ------------------ get_fds_line ------------------------ */

int get_fds_line(FILE *stream, 
                 char *fdsbuffer, char *repbuffer, unsigned int len_fdsbuffer, 
                 int *irep, int *nrep){
  int copyback=0;
  size_t lenbuffer2;
  char buffer[LENBUFFER], buffer2[LENBUFFER], buffer3[LENBUFFER];
  int is_command=0;

  if(*nrep==0){
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
    if(match(fdsbuffer,"&OBST",5)==1){
      int repeat;
      int ibeg, iend;
      float dxyz[6];
      char *crepeat, *cdxyz;

      strcpy(buffer3,fdsbuffer);
      strcpy(repbuffer,buffer3);
      if(get_irvals(buffer3, "COPIES", 1, &repeat, NULL,&ibeg,&iend)==1&&repeat>0){

        crepeat=strstr(buffer3,"COPIES");
        if(crepeat!=NULL)ibeg=crepeat-buffer3;
        subst_string(buffer3,ibeg,iend,NULL); // remove REPEAT sub-string

        if(get_irvals(buffer3, "DXYZ", 6, NULL, dxyz, &ibeg,&iend)==6){

          cdxyz=strstr(buffer3,"DXYZ");
          if(cdxyz!=NULL)ibeg=cdxyz-buffer3;
          subst_string(buffer3,ibeg,iend,NULL); // remove DXYZ sub-string

          *nrep=repeat;
          *irep=0;
          strcpy(fdsbuffer,buffer3);
        }
      }
    }
    return (int)strlen(fdsbuffer);
  }
  else{
    int repeat;
    int ibeg, iend;
    float xb[6],dxyz[6];
    int i;
    char xbstring[256], xb2string[256], *cbeg;

    (*irep)++;
    strcpy(buffer3,repbuffer);
    if(get_irvals(buffer3, "COPIES", 1, &repeat, NULL,&ibeg,&iend)==1&&repeat>0){
      cbeg=strstr(buffer3,"COPIES");
      if(cbeg!=NULL)ibeg=cbeg-buffer3;
      subst_string(buffer3,ibeg,iend,NULL);  // remove REPEAT sub-string
      if(get_irvals(buffer3, "DXYZ", 6, NULL, dxyz, &ibeg,&iend)==6){
        cbeg=strstr(buffer3,"DXYZ");
        if(cbeg!=NULL)ibeg=cbeg-buffer3;
        subst_string(buffer3,ibeg,iend,NULL); // remove DXYZ sub-string
        if(get_irvals(buffer3, "XB", 6, NULL, xb, &ibeg,&iend)==6){
          cbeg=strstr(buffer3,"XB");
          if(cbeg!=NULL)ibeg=cbeg-buffer3;
          for(i=0;i<6;i++){
            xb[i]+=(*irep-1)*dxyz[i];
          }
          float2string(xb,6,xbstring);
          strcpy(xb2string,"XB=");
          strcat(xb2string,xbstring);
          strcat(xb2string," ");
          subst_string(buffer3,ibeg,iend,xb2string);
          strcpy(fdsbuffer,buffer3);
          return (int)strlen(fdsbuffer);
        }
      }
    }
  }
  *nrep=0;
  *irep=0;
  return 0;
}

/* ------------------ expand_shell ------------------------ */

#define OBST_SHELL(SIDE) \
  for(i=0;i<6;i++){\
    xb2[i]=xb[i];\
  }\
  SIDE;\
  strcpy(obst_string,buffer);\
  float2string(xb2,6,xbstring);\
  subst_string(obst_string,ibeg,iend,xbstring);\
  fprintf(stream_out,"%s\n",obst_string)

void expand_shell(FILE *stream_out, char *buffer){
  float xb[6], xb2[6], delta;
  char xbstring[MAXLINE];
  char obst_string[MAXLINE];
  char *delta_beg;
  int i;
  int ibeg, iend;
  int have_delta=1;
  float fullblock=0;
  char *sides_beg;
  int sides[6];

  if(get_irvals(buffer, "DELTA", 1, NULL, &delta,&ibeg, &iend)!=1){
    delta=0.0;
    have_delta=0;
  }
  delta_beg=strstr(buffer,"DELTA");
  if(delta_beg!=NULL){
    ibeg=(int)(delta_beg-buffer);
    subst_string(buffer,ibeg,iend,NULL);  // remove DELTA substring
  }

  for(i=0;i<6;i++){
    sides[i]=1;
  }
  get_irvals(buffer, "SIDES", 6, sides, NULL,&ibeg, &iend);
  for(i=0;i<6;i++){
    if(sides[i]!=0)sides[i]=1;
  }
  sides_beg=strstr(buffer,"SIDES");
  if(sides_beg!=NULL){
    ibeg=(int)(sides_beg-buffer);
    subst_string(buffer,ibeg,iend,NULL);  // remove SIDES substring
  }

  subst_string(buffer,0,5,"&OBST"); // replace initial &SHELL with &OBST

  if(get_irvals(buffer, "XB", 6, NULL, xb, &ibeg, &iend)!=6)return;

  if(have_delta==0||delta<0.000001)return;

  if(delta>(xb[1]-xb[0])/2.0)fullblock=1;
  if(delta>(xb[3]-xb[2])/2.0)fullblock=1;
  if(delta>(xb[5]-xb[4])/2.0)fullblock=1;

  if(fullblock==0){
    if(sides[0]==1){
      OBST_SHELL(xb2[1]=xb[0]+delta;xb2[2]=xb[2]+delta;xb2[3]=xb[3]-delta;xb2[4]=xb[4]+delta;xb2[5]=xb[5]-delta); // left
    }
    if(sides[1]==1){
      OBST_SHELL(xb2[0]=xb[1]-delta;xb2[2]=xb[2]+delta;xb2[3]=xb[3]-delta;xb2[4]=xb[4]+delta;xb2[5]=xb[5]-delta); // right
    }
    if(sides[2]==1){
      OBST_SHELL(xb2[3]=xb[2]+delta;xb2[4]=xb[4]+delta;xb2[5]=xb[5]-delta); // front
    }
    if(sides[3]==1){
      OBST_SHELL(xb2[2]=xb[3]-delta;xb2[4]=xb[4]+delta;xb2[5]=xb[5]-delta); // back
    }
    if(sides[4]==1){
      OBST_SHELL(xb2[5]=xb[4]+delta); // bottom
    }
    if(sides[5]==1){
      OBST_SHELL(xb2[4]=xb[5]-delta); // top
    }
  }
  else{
    OBST_SHELL(xb2[1]=xb[1]);
  }
}

/* ------------------ get_boundbox ------------------------ */

int get_boundbox(blockaiddata *group,int recurse_level){
  fdsdata *thisline;
  int i;
  float *bb_min, *bb_max, *bb_dxyz;

  group->nloaded++;
  if(group->loadonce==1&&group->nloaded>1)return 0;
  if(group->bb_box_defined==1)return 1;
  

  if(recurse_level>MAXRECURSE){
    printf(" *** Fatal error:  recursion level must be less than %i\n",MAXRECURSE);
    return 1;
  }

  bb_min=group->bb_min;
  bb_max=group->bb_max;
  bb_dxyz=group->bb_dxyz;

  for(thisline=group->first_line->next;thisline->next!=NULL;thisline=thisline->next){
    float xb[6];

    if(thisline->line_after!=NULL&&thisline->line_before!=NULL){
      if(thisline->type==1&&(thisline->is_obst==1||thisline->is_shell==1)){
        for(i=0;i<6;i++){
          xb[i]=thisline->xb[i]-group->xyz0[i/2];
        }
      }
      else if(thisline->type==2){
        char linebuffer[MAXLINE];
        char *id2;
        blockaiddata *group2;        
        float offset2[3], rotate2;
        int j;

        strcpy(linebuffer,thisline->line);

        offset2[0]=0.0;
        offset2[1]=0.0;
        offset2[2]=0.0;
        get_irvals(linebuffer, "XYZ", 3, NULL, offset2, NULL, NULL);

        rotate2=0.0;
        get_irvals(linebuffer, "ROTATE", 1, NULL, &rotate2, NULL, NULL);

        id2=get_keyid(linebuffer,"ID");
        if(id2==NULL){
          printf("ERROR: The keyword, ID, has not been defined properly on the following line:\n");
          printf("%s\n",linebuffer);
          return 1;
        }
        group2=get_group(id2);
        if(group2==NULL){
          printf("ERROR: The group, %s, referenced on the following line has not been defined.\n",id2);
          printf("%s\n",linebuffer);
          return 1;
        }

        if(get_boundbox(group2,recurse_level+1)==0)continue;;

        xb[0]=group2->bb_min[0];
        xb[1]=group2->bb_max[0];
        xb[2]=group2->bb_min[1];
        xb[3]=group2->bb_max[1];
        xb[4]=group2->bb_min[2];
        xb[5]=group2->bb_max[2];

        rotatexy(xb,xb+2,  rotate2, group2->bb_dxyz);
        rotatexy(xb+1,xb+3,rotate2, group2->bb_dxyz);
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
  return 1;
}

/* ------------------ init_boundbox0 ------------------------ */

void init_boundbox0(void){
  fdsdata *thisline;
  blockaiddata *group;

  for(group=blockaid_first->next;group->next!=NULL;group=group->next){
    float xb[6];
    int do_this_assm;
    float *bb_min, *bb_max, *bb_dxyz;

    if(group->bb_box_defined==1)continue;

    bb_min=group->bb_min;
    bb_max=group->bb_max;
    bb_dxyz=group->bb_dxyz;

    do_this_assm=1;
    for(thisline=group->first_line->next;thisline->next!=NULL;thisline=thisline->next){

      if(thisline->type==1&&(thisline->is_obst==1||thisline->is_shell==1)){
      }
      else if(thisline->type==2){
        do_this_assm=0;
        break;
      }
    }
    if(do_this_assm==0)continue;

    for(thisline=group->first_line->next;thisline->next!=NULL;thisline=thisline->next){
      int i;

      if(thisline->type!=1||(thisline->is_obst!=1&&thisline->is_shell!=1))continue;

      for(i=0;i<6;i++){
        xb[i]=thisline->xb[i]-group->xyz0[i/2];
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
    group->bb_box_defined=1;
  }
}
/* ------------------ expand_group ------------------------ */

void expand_group(FILE *stream_out, char *buffer, int recurse_level){
  float offset[3], rotate;
  char *id;
  blockaiddata *group;
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

  id=get_keyid(buffer,"ID");
  if(id==NULL){
    printf("ERROR: The keyword, ID, is not defined properly on the following line:\n");
    printf("%s\n",buffer);
    return;
  }

  group=get_group(id);
  if(group==NULL){
    printf("ERROR: The group, %s, referenced on the following line has not been defined.\n",id);
    printf("%s\n",buffer);
    return;
  }

  group->nloaded++;
  if(group->loadonce==1&&group->nloaded>1)return;

  grouplist[recurse_level]=group;

  offset_rotate[4*recurse_level  ] =offset[0];
  offset_rotate[4*recurse_level+1]=offset[1];
  offset_rotate[4*recurse_level+2]=offset[2];
  offset_rotate[4*recurse_level+3]=rotate ;


  for(i=0;i<recurse_level;i++){
    if(strcmp(id,grouplist[i]->id)==0){
      printf(" **** warning ****\n");
      printf("      Block defintions with Id's: \n");
      for(j=0;j<recurse_level;j++){
        printf(" %s,",grouplist[j]->id);
      }
      printf(" %s\n",grouplist[recurse_level]->id);
      printf(" are defined circularly.  Their expansion is halted\n");
      printf(" **** warning ****\n");
      return;
    }
  }

  fprintf(stream_out,"\n *** GROUP: %s ",group->id);
  if(group->type==1){
    fprintf(stream_out,"offset=%f,%f,%f rotate=%f\n",
      offset[0],offset[1],offset[2],rotate);
  }
  fprintf(stream_out,"\n");

  for(thisline=group->first_line->next;thisline->next!=NULL;thisline=thisline->next){
    float xb[6];

    if(thisline->line_after==NULL||thisline->line_before==NULL){
      if(thisline->type==1&&thisline->line!=NULL){
        char *block_key;
        char line2[MAXLINE];

        nkeyvalstack++;
        keyvalstack[nkeyvalstack-1].keyword_list=group->keyword_list;
        keyvalstack[nkeyvalstack-1].val_list=group->val_list;
        keyvalstack[nkeyvalstack-1].nkeywordlist=group->nkeywordlist;

        block_key=strstr(thisline->line,"#");
        if(block_key!=NULL){
          char line_before2[MAXLINE];

          strcpy(line_before2,thisline->line);
          subst_keys(line_before2,recurse_level);
          strcpy(line2,line_before2);
        }
        else{
          strcpy(line2,thisline->line);
        }
        fprintf(stream_out,"%s\n",line2);
      }
      if(thisline->type==2&&thisline->line!=NULL){
        char linebuffer[MAXLINE];

        strcpy(linebuffer,thisline->line);
        expand_group(stream_out,linebuffer,recurse_level+1);      
      }
    }

    if(thisline->line_after!=NULL&&thisline->line_before!=NULL){
      if(thisline->type==1){
        float *xyz, *rotate;
        char *block_key;
        char line2[MAXLINE];

        nkeyvalstack++;
        keyvalstack[nkeyvalstack-1].keyword_list=group->keyword_list;
        keyvalstack[nkeyvalstack-1].val_list=group->val_list;
        keyvalstack[nkeyvalstack-1].nkeywordlist=group->nkeywordlist;

        block_key=strstr(thisline->line_before,"#");
        if(block_key!=NULL){
          char line_before2[10000];

          strcpy(line_before2,thisline->line_before);
          subst_keys(line_before2,recurse_level);
          strcpy(line2,line_before2);
        }
        else{
          strcpy(line2,thisline->line_before);
        }
        for(i=0;i<6;i++){
          xb[i]=thisline->xb[i]-group->xyz0[i/2];
        }

        for(i=recurse_level;i>=0;i--){
          blockaiddata *groupi;

          groupi = grouplist[i];
          xyz = offset_rotate+4*i;
          rotate = offset_rotate+4*i+3;
          rotatexy(xb,xb+2,  rotate[0],groupi->bb_dxyz);
          rotatexy(xb+1,xb+3,rotate[0],groupi->bb_dxyz);
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
            strcat(line2,charxb);
          }
          else{
            strcat(line2,charxb);
            strcat(line2,",");
          }
        }
        block_key=strstr(thisline->line_after,"#");
        if(block_key!=NULL){
          char line_after2[10000];

          strcpy(line_after2,thisline->line_after);
          subst_keys(line_after2,recurse_level);
          strcat(line2,line_after2);
        }
        else{
          strcat(line2,thisline->line_after);
        }
        if(thisline->is_shell==0){
          fprintf(stream_out,"%s\n",line2);
        }
        else{
          expand_shell(stream_out,line2);
        }
        nkeyvalstack--;

      }
      else if(thisline->type==2){
        char linebuffer[MAXLINE];
        blockaiddata *blockaidi;

        nkeyvalstack++;
        blockaidi = thisline->blockaid;
        if(blockaidi!=NULL){
          keyvalstack[nkeyvalstack-1].keyword_list=blockaidi->keyword_list;
          keyvalstack[nkeyvalstack-1].val_list=blockaidi->val_list;
          keyvalstack[nkeyvalstack-1].nkeywordlist=blockaidi->nkeywordlist;
        }
        else{
          keyvalstack[nkeyvalstack-1].keyword_list=NULL;
          keyvalstack[nkeyvalstack-1].val_list=NULL;
          keyvalstack[nkeyvalstack-1].nkeywordlist=0;
        }

        nkeyvalstack++;
        keyvalstack[nkeyvalstack-1].keyword_list=group->keyword_list;
        keyvalstack[nkeyvalstack-1].val_list=group->val_list;
        keyvalstack[nkeyvalstack-1].nkeywordlist=group->nkeywordlist;

        strcpy(linebuffer,thisline->line);
        expand_group(stream_out,linebuffer,recurse_level+1);
        nkeyvalstack-=2;
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

void subst(char *keybeg,char *keyend,char *val){
  char lineend[10000];
  size_t len_from, len_to;
  size_t i;
  char *c;
  char *valcopy;

  if(keybeg==NULL||keyend==NULL)return;

  len_from=keyend-keybeg;
  if(len_from<1)return;

  len_to=strlen(val);
  if(len_to<0)return;

  strcpy(lineend,keyend);
  c=keybeg;
  valcopy=val;
  for(i=0;i<len_to;i++){
    *c++=*valcopy++;
  }
  for(i=0;i<strlen(lineend);i++){
    *c++=lineend[i];
  }
  *c='\0';

}

/* ------------------ get_val ------------------------ */

char *get_val(char *key, int recurse_level){
  int i,j;

  if(key==NULL||strlen(key)==0)return NULL;

  for(i=0;i<nkeyvalstack;i++){
    char **keylist, **vallist;

    keylist=keyvalstack[i].keyword_list;
    if(keylist==NULL)continue;
    vallist=keyvalstack[i].val_list;
    for(j=0;j<keyvalstack[i].nkeywordlist;j++){
      if(keylist[j]==NULL)continue;
      if(strcmp(keylist[j],key)==0){
        return vallist[j];
      }
    }
  }
  return NULL;
}

/* ------------------ get_key ------------------------ */

void get_key(char *line, char *key, char **keybeg, char **keyend){
  char *keyptr, *lineptr;

  *keybeg=strstr(line,"#");
  if(*keybeg==NULL)return;
  keyptr=key;
  lineptr=*keybeg;
  for(;;){
    *keyptr=*lineptr;
    if(*lineptr=='\0'||*lineptr==' '||*lineptr==','||*lineptr=='/'||*lineptr=='\''){
      *keyptr='\0';
      *keyend=lineptr;
      trim(key);
      return;
    }
    lineptr++;
    keyptr++;
  }
}
/* ------------------ subst_keys ------------------------ */

void subst_keys(char *line_in, int recurse_level){
  char line[10000],key[100];
  char *keybeg, *keyend, *val;

  strcpy(line,line_in);
  get_key(line,key, &keybeg, &keyend);
  if(keybeg!=NULL){
    val=get_val(key,recurse_level);
  }
  else{
    val=NULL;
  }
  while(val!=NULL){
    subst(keybeg,keyend,val);

    get_key(line,key, &keybeg, &keyend);
    if(keybeg!=NULL){
      val=get_val(key,recurse_level);
    }
    else{
      val=NULL;
    }
  }
  strcpy(line_in ,line);
  get_key(line,key, &keybeg, &keyend);
  if(keybeg!=NULL){
    *keyend='\0';
    trim(keybeg);
    printf("WARNING: No value found for the key %s\n",keybeg);
  }
  return;
}
  

/* ------------------ get_keywords ------------------------ */

void get_keywords(blockaiddata *blockaidi, char *buffer){
  char buffer2[MAXLINE], *buffptr;
  char *key, *endkey, *val, *endval;
  int nkeys;
  int i;

  for(i=0;i<blockaidi->nkeywordlist;i++){
    char *key, *val;

    key=blockaidi->keyword_list[i];
    val=blockaidi->val_list[i];
    FREEMEMORY(key);
    FREEMEMORY(val);
  }
  FREEMEMORY(blockaidi->keyword_list);
  FREEMEMORY(blockaidi->val_list);
  blockaidi->nkeywordlist=0;

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
    char *ckey, *cval, *endline;

    endkey=strstr(key+1,"=");
    if(endkey==NULL)break;

    val=strstr(endkey+1,"'");
    if(val==NULL)break;
    val++;
    
    endline=strstr(val,"\n");
    endval=strstr(val,"'");
    if(endval==NULL)break;
    *endkey='\0';
    *endval='\0';
    trim(key);
    trim(val);
    lenkey=strlen(key);
    lenval=strlen(val);
    if(endline!=NULL&&endval>endline){
      fprintf(stderr,"**********\n");
      printf("ERROR: A keyword value can't be split across a line\n");
      printf("    keyword: %s\n",key);
      printf("    value: %s\n",val);
      fprintf(stderr,"**********\n");
      break;
    }
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
  for(i=0;i<nkeys;i++){
    char *key, *val;

    key=blockaidi->keyword_list[i];
    val=blockaidi->val_list[i];
#ifdef _DEBUG
    printf("key=%s val=%s\n",key,val);
#endif
  }
  return;
}

/* ------------------ create_group ------------------------ */

blockaiddata *create_group(char *buffer){
  blockaiddata *blockaidi, *bprev, *bnext;
  float *orig;
  char *id;
  fdsdata *first_line, *last_line;

  NewMemory((void **)&blockaidi,sizeof(blockaiddata));
  blockaidi->nloaded=0;
  if(strstr(buffer,"LOADONCE")==NULL){
    blockaidi->loadonce=0;
  }
  else{
    blockaidi->loadonce=1;
  }
  blockaidi->type=0;
  blockaidi->keyword_list=NULL;
  blockaidi->val_list=NULL;
  blockaidi->nkeywordlist=0;
  bprev=blockaid_first;
  bnext=blockaid_first->next;
  blockaidi->prev=bprev;
  blockaidi->next=bnext;
  bprev->next=blockaidi;
  bnext->prev=blockaidi;
  get_keywords(blockaidi,buffer);

  orig=blockaidi->xyz0;

  if(get_irvals(buffer, "ORIGIN", 3, NULL, orig, NULL, NULL)==3){
    // first 3 positions of orig are defined in get_irvals
    orig[3]=1.0;
  }
  else{
    orig[0]=0.0;
    orig[1]=0.0;
    orig[2]=0.0;
    orig[3]=1.0;
  }
  id=get_keyid(buffer,"ID");
  if(id!=NULL){
    NewMemory((void **)&blockaidi->id,strlen(id)+1);
    strcpy(blockaidi->id,id);
  }
  else{
    printf("ERROR: The keyword, ID, is not defined properly on the following line:\n");
    printf("%s\n",buffer);
    return NULL;
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

/* ------------------ update_group ------------------------ */

void update_group(blockaiddata *group,char *buffer){
  fdsdata *prev, *next, *thisfds;
  size_t len;
  int is_obst=0;
  int is_shell=0;

  next=group->last_line;
  prev=next->prev;

  NewMemory((void **)&thisfds,sizeof(fdsdata));

  thisfds->line=NULL;
  thisfds->is_obst=0;
  thisfds->is_shell=0;
  thisfds->blockaid=NULL;
  if(match(buffer,"&GRP",4)==1){
    blockaiddata *blockaid;

    thisfds->type=2;
    if(strstr(buffer,"#")!=NULL){

      NewMemory((void **)&blockaid,sizeof(blockaiddata ));
      blockaid->nloaded=0;
      blockaid->loadonce=0;
      blockaid->keyword_list=NULL;
      blockaid->val_list=NULL;
      blockaid->nkeywordlist=0;
      if(blockaid!=NULL){
        get_keywords(blockaid, buffer);
      }
    }
    else{
      blockaid=NULL;
    }
    thisfds->blockaid=blockaid;
  }
  else if(match(buffer,"&BGRP",5)==1||match(buffer,"&EGRP",5)==1){
    printf("ERROR: &BGRP and &EGRP cannot be nested\n");
    thisfds->type=0;
  }
  else{
    thisfds->type=1;
    if(
      match(buffer,"&OBST",5)==1||
      match(buffer,"&VENT",5)==1||
      match(buffer,"&SHELL",6)==1||
      match(buffer,"&HOLE",5)==1
      ){
        group->type=1;
    }
    if(match(buffer,"&OBST",5)==1){
      is_obst=1;
      thisfds->is_obst=1;
    }
    if(match(buffer,"&SHELL",6)==1){
      is_shell=1;
      thisfds->is_shell=1;
    }
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
      if(is_obst==1||is_shell==1){

        xb = thisfds->xb;
        if(group->xyz0[3]<0.0){
          float *orig;

          orig = group->xyz0;
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
        if(is_shell==1){
          get_irvals(buffer, "DELTA", 1, NULL, &thisfds->delta,NULL,NULL);
        }
      }
    }
    else if(thisfds->type==2){
      get_irvals(buffer, "XYZ", 3, NULL, thisfds->xb,&thisfds->ibeg,&thisfds->iend);
    }

    thisfds->line_before=NULL;
    thisfds->line_after=NULL;
    if(thisfds->type!=0){
      if(thisfds->ibeg>=0&&thisfds->iend>=0){
        thisfds->linecopy[thisfds->ibeg]=0;
        thisfds->line_before=thisfds->linecopy;
        thisfds->line_after=thisfds->linecopy+ thisfds->iend;
      }
    }


  }
  prev->next=thisfds;
  thisfds->prev=prev;
  thisfds->next=next;
  next->prev=thisfds;
}

/* ------------------ remove_groupb ------------------------ */

void remove_group(blockaiddata *groupb){
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
//#ifdef _DEBUG
  for(assm=blockaid_first->next;assm->next!=NULL;assm=assm->next){
    float *x1, *x2, *dx;

    get_boundbox(assm,0);

    x1 = assm->bb_min;
    x2 = assm->bb_max;
    dx = assm->bb_dxyz;

#ifdef _DEBUG
    printf("assm=%s\n",assm->id);
    printf("  bounds=(%f,%f), (%f,%f), (%f,%f)\n",x1[0],x2[0],x1[1],x2[1],x1[2],x2[2]);
#endif
  }
  for(assm=blockaid_first->next;assm->next!=NULL;assm=assm->next){
    assm->nloaded=0;
  }
//#endif
}

/* ------------------ get_group ------------------------ */

blockaiddata *get_group(char *id){
  blockaiddata *assm;

  if(id==NULL)return NULL;
  for(assm=blockaid_first->next;assm->next!=NULL;assm=assm->next){
    if(strcmp(assm->id,id)==0)return assm;
  }
  return NULL;
}


/* ------------------ init_groupdata ------------------------ */

void init_groupdata(char *id, float *orig, blockaiddata *prev, blockaiddata *next){
  blockaiddata *newgroup;
  fdsdata *fl, *ll;

  NewMemory((void **)&newgroup,sizeof(blockaiddata));
  strcpy(newgroup->id,id);
  newgroup->loadonce=0;
  newgroup->nloaded=0;
  newgroup->xyz0[0]=orig[0];
  newgroup->xyz0[1]=orig[1];
  newgroup->xyz0[2]=orig[2];
  newgroup->prev=prev;
  newgroup->next=next;
  newgroup->keyword_list=NULL;
  newgroup->val_list=NULL;
  newgroup->nkeywordlist=0;


  prev->next=newgroup;
  next->prev=newgroup;
  NewMemory((void **)&newgroup->first_line,sizeof(fdsdata));

  NewMemory((void **)&newgroup->last_line,sizeof(fdsdata));
  fl = newgroup->first_line;
  ll = newgroup->last_line;

  fl->line=NULL;
  fl->prev=NULL;
  fl->next=ll;

  ll->line=NULL;
  ll->prev=fl;
  ll->next=NULL;
}
