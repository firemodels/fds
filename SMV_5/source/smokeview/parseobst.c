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

void headsubst(char *line, FILE *stream_out);
void parseobst_xb(char *line, int *i1, int *i2);
void parseobst_surf(char *line, int *i1, int *i2);
void xb_obstsubst(char *oldobst, blockagedata *bc);
void output_new_obst(FILE *stream_out, FILE *stream_out2);
void output_obst(char *buffer,int *smvobstcount, int fdsobstcount, FILE *stream_out, FILE *stream_out2);
void obst_valsout(char *buffer,float xb1,float xb2,float yb1,float yb2,float zb1,float zb2);

/* ------------------ outputchangedblockages ------------------------ */

void outputchangedblockages(void){
  FILE *streamout;

  if(fds_filein!=NULL&&fds_fileout!=NULL&&fds_fileout2){
    convert_fdsfile(fds_filein, fds_fileout, fds_fileout2);
  }
  if(fds_filein==NULL){
    streamout=NULL;
    if(fds_fileout!=NULL)streamout=fopen(fds_fileout,"w");
    if(streamout==NULL)streamout=stdout;
    if(streamout!=stdout)fclose(streamout);
    printf("\n*** Warning: The FDS input file name was either not specified \n");
    printf(" in the smokeview (.smv) file or does not exist. Add the keyword/input \n");
    printf(" file name pair:\n\nINPF\ncasename.data\n\n");
    printf(" where casename.data is the name of the FDS input file\n\n");
  }


}

/* ------------------ convert_fdsfile ------------------------ */

void convert_fdsfile(const char *filein, const char *fileout, const char *fileout2){

  FILE *stream_in, *stream_out, *stream_out2;
  char *buffer2=NULL,buffer[1000];
  size_t sizebuffer2;
  size_t lenbuffer2,lenbuffer;
  int fdsobstcount=0,smvobstcount=0;
  int copyback;
  int checkhead=1;

  if(filein==NULL||fileout==NULL||fileout2==NULL)return;
  stream_in = fopen(filein,"r");
  if(stream_in==NULL)return;

  stream_out = fopen(fileout,"w");
  if(stream_out==NULL){
    fclose(stream_in);
    return;
  }

  stream_out2 = fopen(fileout2,"w");
  if(stream_out2==NULL){
    fclose(stream_in);
    fclose(stream_out);
    return;
  }

  sizebuffer2=1000;
  if(NewMemory((void **)&buffer2,(unsigned int)(sizebuffer2+1))==0){
    fclose(stream_in);
    fclose(stream_out);
    fclose(stream_out2);
    return;
  }

  while(!feof(stream_in)){
    if(fgets(buffer,1000,stream_in)==NULL)break;
//    trim(buffer);

    /* output non &OBST lines without change */
    if(checkhead==1&&STRSTR(buffer,"&HEAD")!=NULL){
      checkhead=0;
      headsubst(buffer,stream_out);
      continue;
    }

    if(STRSTR(buffer,"&OBST")==NULL){
      fputs(buffer,stream_out);
      continue;
    }
    fdsobstcount++;

    if(fdsobstcount>0&&fdsobstcount<=nchanged_idlist){
      if(changed_idlist[fdsobstcount]==0){
        fputs(buffer,stream_out);
        continue;
      }
    }


    /* concatenate multi-line &OBST lines into 1 */

    strcpy(buffer2,buffer);
    copyback=0;
    while(strstr(buffer,"/")==NULL){
      fgets(buffer,1000,stream_in);
      lenbuffer=strlen(buffer); 
      lenbuffer2=strlen(buffer2);
      if(lenbuffer2+lenbuffer+2>sizebuffer2){
        sizebuffer2 = lenbuffer2+lenbuffer+2+1000;
        ResizeMemory((void **)&buffer2,(unsigned int)sizebuffer2);
      }
      strcat(buffer2,buffer);
      copyback=1;
    }
    if(copyback==1)strcpy(buffer,buffer2);
    output_obst(buffer,&smvobstcount,fdsobstcount,stream_out,stream_out2);
  }
  smvobstcount=0;
  if(fdsobstcount==0)output_obst(buffer,&smvobstcount,0,stream_out,stream_out2);
  output_new_obst(stream_out, stream_out2);


  FREEMEMORY(buffer2);
  fclose(stream_in);
  fclose(stream_out);
  fclose(stream_out2);

}

/* ------------------ output_new_obst ------------------------ */

void output_new_obst(FILE *stream_out, FILE *stream_out2){
  int i;
  float xb1, xb2, yb1, yb2, zb1, zb2;
  float *xplt, *yplt, *zplt;
  blockagedata *bc;
  mesh *meshi;
  char buffer2[1000];
  char buffer[1000];


  for(i=0;i<nselectblocks;i++){
    bc = selectblockinfo[sortedblocklist[i]];
    if(bc->id==-1){
      meshi = meshinfo + bc->meshindex;
      xplt = meshi->xplt_orig;
      yplt = meshi->yplt_orig;
      zplt = meshi->zplt_orig;
      xb1 = xplt[bc->ijk[IMIN]];
      xb2 = xplt[bc->ijk[IMAX]];
      yb1 = yplt[bc->ijk[JMIN]];
      yb2 = yplt[bc->ijk[JMAX]];
      zb1 = zplt[bc->ijk[KMIN]];
      zb2 = zplt[bc->ijk[KMAX]];

      strcpy(buffer,"&OBST XB=");
      obst_valsout(buffer2,xb1,xb2,yb1,yb2,zb1,zb2);
      strcat(buffer,buffer2);
      strcat(buffer,"/");
      if(bc->label!=NULL){
        strcat(buffer,"  ");
        strcat(buffer,bc->label);
      }
      fprintf(stream_out,"%s\n",buffer);
      fprintf(stream_out2,"%s\n",buffer);
      continue;
    }
  }
}

/* ------------------ output_obst ------------------------ */

void output_obst(char *buffer,int *smvobstcount, int fdsobstcount, FILE *stream_out, FILE *stream_out2){
  int i;
  mesh *meshi;
  float *xplt, *yplt, *zplt;
  blockagedata *bc;
  float xb1, xb2, yb1, yb2, zb1, zb2;
  char *slashptr;
  size_t len;
  char buffer2[1000];

  for(i=*smvobstcount;i<nselectblocks;i++){
    bc = selectblockinfo[sortedblocklist[i]];
    *smvobstcount=i;
    if(bc->id==-1)continue;
    if(bc->id==-1){
      meshi = meshinfo + bc->meshindex;
      xplt = meshi->xplt_orig;
      yplt = meshi->yplt_orig;
      zplt = meshi->zplt_orig;
      xb1 = xplt[bc->ijk[IMIN]];
      xb2 = xplt[bc->ijk[IMAX]];
      yb1 = yplt[bc->ijk[JMIN]];
      yb2 = yplt[bc->ijk[JMAX]];
      zb1 = zplt[bc->ijk[KMIN]];
      zb2 = zplt[bc->ijk[KMAX]];

  //    sprintf(buffer,"&OBST XB=%f, %f, %f, %f, %f, %f/",xb1,xb2,yb1,yb2,zb1,zb2);
      strcpy(buffer,"&OBST XB=");
      obst_valsout(buffer2,xb1,xb2,yb1,yb2,zb1,zb2);
      strcat(buffer,buffer2);
      strcat(buffer,"/");
      if(bc->label!=NULL){
        strcat(buffer,"  ");
        strcat(buffer,bc->label);
      }
      fprintf(stream_out,"%s\n",buffer);
      fprintf(stream_out2,"%s\n",buffer);
      continue;
    }
    if(bc->id>fdsobstcount)break;
    if(bc->id==fdsobstcount){
      if(bc->changed==1||bc->changed_surface==1||
        (bc->id>0&&bc->id<nchanged_idlist&&changed_idlist[bc->id]==1)){
        xb_obstsubst(buffer,bc);
      }
      slashptr = strstr(buffer,"/");
      if(slashptr!=NULL){
        slashptr[1]='\0';
        strcat(buffer,"  ");
        if(bc->label!=NULL){
          len=strlen(bc->label);
          if(len>0&&bc->label[0]!='*'){
            strcat(buffer,bc->label);
          }
        }
        strcat(buffer,"\n");
      }
      if(bc->del==0&&bc->hidden==0)fputs(buffer,stream_out);

    }
  }
}

/* ------------------ getkeyparam ------------------------ */

char *getkeyparam(char *source, int *n, const char *key){
  size_t len;
  char *keyptr,*s1,*e1;

  *n=0;
  if(source==NULL||key==NULL)return NULL;
  len=strlen(key);
  if(len==0)return NULL;
  keyptr=STRSTR(source,key);
  if(keyptr==NULL)return NULL;
  s1=strstr(keyptr,"'");
  if(s1==NULL){
    s1=strstr(keyptr,"\"");
    if(s1==NULL)return NULL;
    s1++;
    e1=strstr(s1,"\"");
    if(e1!=NULL){
      e1[0]='\0';
    }
    *n=strlen(s1);
    return s1;
  }
  s1++;
  e1=strstr(s1,"'");
  if(e1!=NULL){
    e1[0]='\0';
  }
  *n=strlen(s1);
  return s1;
}

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

/* ------------------ xb_obstsubst ------------------------ */

void xb_obstsubst(char *source,blockagedata *bc){
  int i1, i2;
  float xb1, xb2, yb1, yb2, zb1, zb2;
  float *xplt, *yplt, *zplt;
  mesh *meshi;
  char buffer[1000];
  char buffer2[1000];
  char buffer3[1000];
  size_t len;
  char *slashptr;
  char *surflabel;

  if(source==NULL)return;
  slashptr=strstr(source,"/");
  if(slashptr!=NULL)slashptr[1]='\0';
  i1=-1;i2=-1;

  parseobst_xb(source,&i1,&i2);
  if(i1>=0&&i2>=0){
    strncpy(buffer,source,(unsigned int)i1);
    buffer[i1]='\0';
    meshi=meshinfo + bc->meshindex;
    xplt = meshi->xplt_orig;
    yplt = meshi->yplt_orig;
    zplt = meshi->zplt_orig;
    xb1 = xplt[bc->ijk[IMIN]];
    xb2 = xplt[bc->ijk[IMAX]];
    yb1 = yplt[bc->ijk[JMIN]];
    yb2 = yplt[bc->ijk[JMAX]];
    zb1 = zplt[bc->ijk[KMIN]];
    zb2 = zplt[bc->ijk[KMAX]];
    strcpy(buffer2,"XB=");
    obst_valsout(buffer3,xb1,xb2,yb1,yb2,zb1,zb2);
    strcat(buffer2,buffer3);
//    sprintf(buffer2," XB=%f, %f, %f, %f, %f, %f ",xb1,xb2,yb1,yb2,zb1,zb2);
    len = strlen(buffer2);
    if(len+strlen(source+i2+1)<1000)strcpy(buffer2+len,source+i2+1);
    strcat(buffer,buffer2);
    strcpy(source,buffer);
  }

  i1=-1;i2=-1;
  if(bc->changed_surface==1){
    parseobst_surf(source,&i1,&i2);
    if(i1<0&&i2<0){
      slashptr=strstr(source,"/");
      if(slashptr!=NULL){
        i1 = slashptr-source;
        i2=i1-1;
      }
    }
    if(i1>=0&&i2>=0){
      strncpy(buffer,source,(unsigned int)i1);
      buffer[i1]='\0';
      switch (bc->walltype){
       case WALL_1:
         surflabel=bc->surf[0]->surfacelabel;
         if(strlen(surfacedefaultlabel)==0||strcmp(surfacedefaultlabel,surflabel)!=0){
          sprintf(buffer2," SURF_ID='%s'",bc->surf[0]->surfacelabel);
         }
         else{
           buffer2[0]='\0';
         }
        break;

       case WALL_3:
         sprintf(buffer2," SURF_IDS='%s','%s','%s'",
         bc->surf[UP_Z]->surfacelabel,
         bc->surf[UP_Y]->surfacelabel,
         bc->surf[DOWN_Z]->surfacelabel
         );
       break;

       case WALL_6:
       sprintf(buffer2," SURF_ID6='%s','%s','%s','%s','%s','%s'",
         bc->surf[DOWN_X]->surfacelabel,
         bc->surf[UP_X]->surfacelabel,
         bc->surf[DOWN_Y]->surfacelabel,
         bc->surf[UP_Y]->surfacelabel,
         bc->surf[DOWN_Z]->surfacelabel,
         bc->surf[UP_Z]->surfacelabel
         );
       break;
       default:
         ASSERT(0);
         break;
      }
      len=strlen(buffer2);
      if(len+strlen(source+i2+1)<1000)strcpy(buffer2+len,source+i2+1);
      strcat(buffer,buffer2);
      strcpy(source,buffer);
    }
  }

}

/* ------------------ obst_valsout ------------------------ */

void obst_valsout(char *buffer,float xb1,float xb2,float yb1,float yb2,float zb1,float zb2){
  char buffer2[1000];
  float vals[6];
  int i;

  vals[0]=xb1;
  vals[1]=xb2;
  vals[2]=yb1;
  vals[3]=yb2;
  vals[4]=zb1;
  vals[5]=zb2;


  strcpy(buffer,"");

  for(i=0;i<6;i++){
    sprintf(buffer2,"%f",vals[i]);
    trimzeros(buffer2);
    strcat(buffer,buffer2);
    if(i!=5)strcat(buffer,",");
  }
}

/* ------------------ parseobst_xb ------------------------ */

void parseobst_xb(char *line, int *i1, int *i2){
  char *l1, *l2, *xbptr, *c;
  char line2[1000];
  size_t len;
  int mode,token;
  unsigned int i;
#define INBLANK 0
#define INTOKEN 1

  len = strlen(line);
  *i1 = -1;
  *i2 = -1;
  if(len<1)return;

  l1=line;
  l2=line2;
  for(i=0;i<len;i++){
    *l2++=(char)toupper(*l1++);
  }
  *l2='\0';
  if(strstr(line2,"&OBST")==NULL)return;
  xbptr = strstr(line2,"XB");
  if(xbptr==NULL)return;
  *i1 = xbptr - line2;
  mode = INBLANK;
  c = xbptr+2;
  token=0;
  for(i=*i1+2;i<len;i++){
    switch (mode) {
    case INBLANK:
      if(*c==' '||*c==','||*c=='\n')break;
      if(token==0&&*c=='=')break;
      token++;
      mode = INTOKEN;
      break;
    case INTOKEN:
      if(*c==' '||*c==','||*c=='/'){
        if(token==6){
          *i2=i-1;
          return;
        }
        if(*c=='/'){
          *i2=-1;
          return;
        }
        mode=INBLANK;
      }
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
    c++;
  }
  *i1=-1;
  return;
}

/* ------------------ parseobst_surf ------------------------ */

void parseobst_surf(char *line, int *i1, int *i2){
  char *l1, *l2, *surfptr,*surfptr1=NULL,*surfptr3,*surfptr6, *c;
  char line2[1000];
  size_t len,i;
  int nquote=0;
  int lensurf;
  int quotes_expected;

  len = strlen(line);
  *i1 = -1;
  *i2 = -1;
  if(len<1)return;

  l1=line;
  l2=line2;
  for(i=0;i<len;i++){
    *l2++=(char)toupper(*l1++);
  }
  *l2='\0';
  if(strstr(line2,"&OBST")==NULL)return;

  surfptr1 = strstr(line2,"SURF_ID");
  if(surfptr1==NULL)return;
  surfptr3 = strstr(line2,"SURF_IDS");
  surfptr6 = strstr(line2,"SURF_ID6");

  if(surfptr3==NULL&&surfptr6==NULL){
    lensurf=7;
    quotes_expected=2;
    surfptr=surfptr1;
  }
  else if(surfptr3!=NULL){
    lensurf=8;
    quotes_expected=6;
    surfptr = surfptr3;
  }
  else if(surfptr6!=NULL){
    surfptr = surfptr6;
    lensurf=8;
    quotes_expected=12;
  }

  *i1 = surfptr - line2;

  c = surfptr+lensurf;
  for(i=*i1+lensurf;i<len;i++){
    if(*c=='\''||*c=='"')nquote++;
    if(nquote==quotes_expected){
      *i2=i;
      return;
    }
    c++;
  }
  *i1=-1;
  return;
}

/* ------------------ headsubst ------------------------ */

void headsubst(char *line, FILE *stream_out){
  char *l1, *l2, *chidptr,*c;
  char line2[1000];
  int lenhead;
  size_t len,i;
  int nquote=0;
  int quotes_expected;
  int i1=-1,i2=-1;
  char *ext;

  len = strlen(line);
  if(len<1)return;

  l1=line;
  l2=line2;
  for(i=0;i<len;i++){
    *l2++=(char)toupper(*l1++);
  }
  *l2='\0';
  if(strstr(line2,"&HEAD")==NULL){
    fputs(line,stream_out);
    return;
  }

  chidptr = strstr(line2,"CHID");
  if(chidptr==NULL){
    fputs(line,stream_out);
    return;
  }

  lenhead=4;
  quotes_expected=2;

  i1 = chidptr - line2;

  c = chidptr+lenhead;
  for(i=i1+lenhead;i<len;i++){
    if(*c=='\''||*c=='"')nquote++;
    if(nquote==quotes_expected){
      i2=i;
      break;
    }
    c++;
  }
  if(i2!=-1){
    for(c=line;c<line+i1;c++){
      putc(*c,stream_out);
    }
    strcpy(line2,fds_fileout);
    ext=strrchr(line2,'.');
    if(ext!=NULL)line2[ext-line2]='\0';
    fprintf(stream_out,"CHID=\'%s\'",line2);
    for(c=line+i2+1;c<line+len;c++){
      putc(*c,stream_out);
    }
  }
  else{
    fputs(line,stream_out);
  }
  return;
}
