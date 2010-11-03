// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>  
#include <stdlib.h>
#include "flowfiles.h"
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "MALLOC.h"
#include "smokeviewvars.h"

// svn revision character string
char IOtarget_revision[]="$Revision$";

/* ------------------ readtarget ------------------------ */

void readtarget(const char *file, int ifile, int flag, int *errorcode){
  FILE *stream;
  float ttargmin, ttargmax;
  char buffer[255];
  int n,colorindex;
  float time;
  int exitloop = 0;
  float *t, *x, *y, *z, *x2, *y2, *z2;
  float xt, yt, zt;
  float xt2, yt2, zt2;
  int i;
  float r, g, b;
  int nsteps;

  CheckMemory;
  *errorcode=0;
  if(targinfo[ifile].type==2){
    readtarget2(file,ifile,flag,errorcode);
    return;
  }
  if(target_positions != NULL){
    for(n=0;n<ntargets;n++){
      FREEMEMORY(target_positions[n].x);
      FREEMEMORY(target_positions[n].y);
      FREEMEMORY(target_positions[n].z);
      FREEMEMORY(target_positions[n].x2);
      FREEMEMORY(target_positions[n].y2);
      FREEMEMORY(target_positions[n].z2);
      FREEMEMORY(target_positions[n].t);
    }
    FREEMEMORY(target_positions);
  }
  FREEMEMORY(targtimes);
  targfilenum = ifile;
  if(flag == UNLOAD){
    targinfo[ifile].loaded=0;
    targinfo[ifile].display=0;
    plotstate=getplotstate(DYNAMIC_PLOTS);
    visTarg = 0;
    ReadTargFile=0;
    ntargets=0;
    updatetimes();
    updatemenu=1;
    return;
  }

  if( (stream=fopen(file,"r"))==NULL)return;

  fgets(buffer,255,stream);
  sscanf(buffer,"%i",&ntargets);
  FREEMEMORY(target_positions);
  NewMemory((void **)&target_positions,ntargets*sizeof(targpos));
  NewMemory((void **)&targtimes,NTARGTIMES*sizeof(float));
  CheckMemory;

  /* read target colors */

  for(n=0;n<ntargets;n++){
    fgets(buffer,255,stream);
    sscanf(buffer,"%i %f %f %f",&colorindex,&r,&g,&b);
    sscanf(buffer,"%i ",&colorindex);
    if(colorindex<0){
      if(r<0.0)r=0.0;
      if(r>1.0)r=1.0;
      if(g<0.0)g=0.0;
      if(g>1.0)g=1.0;
      if(b<0.0)b=0.0;
      if(b>1.0)b=1.0;
      target_positions[n].rgb[0]=r;
      target_positions[n].rgb[1]=g;
      target_positions[n].rgb[2]=b;
    }
    else{
      if(colorindex>=0&&colorindex<nrgb2){
        target_positions[n].rgb[0]=rgb2[colorindex][0];
        target_positions[n].rgb[1]=rgb2[colorindex][1];
        target_positions[n].rgb[2]=rgb2[colorindex][2];
      }
      else{
        target_positions[n].rgb[0]=rgb2[1][0];
        target_positions[n].rgb[1]=rgb2[1][1];
        target_positions[n].rgb[2]=rgb2[1][2];
      }
    }
  }

  for(n=0;n<ntargets;n++){
    if(fgets(buffer,255,stream)==NULL)break;
    sscanf(buffer,"%i ",&nsteps);
    if(nsteps<=0)break;
    target_positions[n].nsteps=nsteps;
    t=NULL;
    x=NULL;
    y=NULL;
    z=NULL;
    x2=NULL;
    y2=NULL;
    z2=NULL;
    CheckMemory;
    NewMemory((void **)&t,nsteps*sizeof(float));
    NewMemory((void **)&x,nsteps*sizeof(float));
    NewMemory((void **)&y,nsteps*sizeof(float));
    NewMemory((void **)&z,nsteps*sizeof(float));
    NewMemory((void **)&x2,nsteps*sizeof(float));
    NewMemory((void **)&y2,nsteps*sizeof(float));
    NewMemory((void **)&z2,nsteps*sizeof(float));
    target_positions[n].t=t;
    target_positions[n].x=x;
    target_positions[n].y=y;
    target_positions[n].z=z;
    target_positions[n].x2=x2;
    target_positions[n].y2=y2;
    target_positions[n].z2=z2;
    for(i=0;i<nsteps;i++){
      if(fgets(buffer,255,stream)!=NULL){
        sscanf(buffer,"%f %f %f %f %f %f %f",&time,&xt,&yt,&zt,&xt2,&yt2,&zt2);
        target_positions[n].t[i]=time;
        target_positions[n].x[i]=xt;
        target_positions[n].y[i]=yt;
        target_positions[n].z[i]=zt;
        target_positions[n].x2[i]=xt2;
        target_positions[n].y2[i]=yt2;
        target_positions[n].z2[i]=zt2;
      }
      else{exitloop=1;}
    }
    if(exitloop==1)break;
  }
  if(ntargets>0){
    ttargmin=1000000.;
    ttargmax=-1000000.;
    for(n=0;n<ntargets;n++){
      nsteps = target_positions[n].nsteps;
      if(nsteps>0){
        if(target_positions[n].t[nsteps-1]<ttargmin)ttargmin=target_positions[n].t[0];
        if(target_positions[n].t[nsteps-1]>ttargmax)ttargmax=target_positions[n].t[nsteps-1];
      }
    }
    if(ttargmax>ttargmin){
      for(n=0;n<NTARGTIMES;n++){
        targtimes[n]=ttargmin+n*(ttargmax-ttargmin)/(NTARGTIMES-1);
      }
    }
  }
  printf("min=%f max=%f\n",ttargmin,ttargmax);
  ReadTargFile=1;
  visTarg=1;
  targinfo[ifile].loaded=1;
  targinfo[ifile].display=1;
  plotstate=getplotstate(DYNAMIC_PLOTS);
  updatetimes();
  fclose(stream);
  IDLE();

}

/* ------------------ readtarget2 ------------------------ */

void readtarget2(const char *file, int ifile, int flag, int *errorcode){
  FILE *stream;
  char buffer[255];
  int n;
  float time;
  float *t, *x, *y, *z;
  float xt, yt, zt;
  int i;
  int nsteps;
  float vals_local[6];
  float *vals,valmin,valmax,val;
  unsigned char *color;
  int doit;

  CheckMemory;
  *errorcode=0;
  if(target_positions != NULL){
    for(n=0;n<ntargets;n++){
      FREEMEMORY(target_positions[n].x);
      FREEMEMORY(target_positions[n].y);
      FREEMEMORY(target_positions[n].z);
      FREEMEMORY(target_positions[n].t);
      FREEMEMORY(target_positions[n].vals);
    }
    FREEMEMORY(target_positions);
  }
  FREEMEMORY(targtimes);
  targfilenum = ifile;
  if(flag == UNLOAD){
    targinfo[ifile].loaded=0;
    targinfo[ifile].display=0;
    plotstate=getplotstate(DYNAMIC_PLOTS);
    visTarg = 0;
    ReadTargFile=0;
    ntargets=0;
    updatetimes();
    updatemenu=1;
    return;
  }

  if( (stream=fopen(file,"r"))==NULL)return;

  fgets(buffer,255,stream);
  fgets(buffer,255,stream);
  fgets(buffer,255,stream);
  ntargtimes=0;
  while(fgets(buffer,255,stream)!=NULL){
    ntargtimes++;
  }
  rewind(stream);

  fgets(buffer,255,stream);
  fgets(buffer,255,stream);
  fgets(buffer,255,stream);
  ntargets=1;
  FREEMEMORY(target_positions);
  NewMemory((void **)&target_positions,ntargets*sizeof(targpos));
  NewMemory((void **)&targtimes,ntargtimes*sizeof(float));
  CheckMemory;

  n=0;
  doit=1;
  if(fgets(buffer,255,stream)==NULL)doit=0;
  if(doit==1){
    sscanf(buffer,"%i ",&nsteps);
    if(nsteps<=0)doit=0;
  }
  if(doit==0){
    fclose(stream);
    readtarget2("", ifile, UNLOAD, errorcode);
  }
  target_positions[n].nsteps=nsteps;
  t=NULL;
  x=NULL;
  y=NULL;
  z=NULL;
  CheckMemory;
  NewMemory((void **)&t,ntargtimes*sizeof(float));
  NewMemory((void **)&x,ntargtimes*sizeof(float));
  NewMemory((void **)&y,ntargtimes*sizeof(float));
  NewMemory((void **)&z,ntargtimes*sizeof(float));
  NewMemory((void **)&vals,6*ntargtimes*sizeof(float));
  NewMemory((void **)&color,6*ntargtimes*sizeof(unsigned char));
  target_positions[n].t=t;
  target_positions[n].x=x;
  target_positions[n].y=y;
  target_positions[n].z=z;
  target_positions[n].color=color;
  target_positions[n].vals=vals;
  for(i=0;i<ntargtimes;i++){
    if(fgets(buffer,255,stream)!=NULL){
      sscanf(buffer,"%f %f %f %f %f %f %f %f %f %f",
        &time,&xt,&yt,&zt,
        vals_local,vals_local+1,vals_local+2,vals_local+3,vals_local+4,vals_local+5
        );
      target_positions[n].t[i]=time;
      targtimes[i]=time;
      target_positions[n].x[i]=xt;
      target_positions[n].y[i]=yt;
      target_positions[n].z[i]=zt;
      target_positions[n].vals[6*i+0]=vals_local[0];
      target_positions[n].vals[6*i+1]=vals_local[1];
      target_positions[n].vals[6*i+2]=vals_local[2];
      target_positions[n].vals[6*i+3]=vals_local[3];
      target_positions[n].vals[6*i+4]=vals_local[4];
      target_positions[n].vals[6*i+5]=vals_local[5];
    }
    else{
      break;
    }
  }

  vals = target_positions[n].vals;
  valmin=vals[0];
  valmax=valmin;
  for(i=1;i<6*ntargtimes;i++){
    if(vals[i]<valmin)valmin=vals[i];
    if(vals[i]>valmax)valmax=vals[i];
  }
  printf("data bounds: min=%f max=%f\n",valmin,valmax);
  if(settargetmin==1)valmin=targetmin;
  if(settargetmax==1)valmax=targetmax;
  printf("applied bounds: min=%f max=%f\n",valmin,valmax);

  target_positions[n].valmin=valmin;
  target_positions[n].valmax=valmax;
  for(i=0;i<6*ntargtimes;i++){
    if(valmax>valmin){
      val = (vals[i]-valmin)/(valmax-valmin);
    }
    else{
      val = 0.5;
    }
    if(val<0.0)val=0.0;
    if(val>1.0)val=1.0;
    color[i]=(unsigned char)(val*255);
  }
  ReadTargFile=1;
  visTarg=1;
  targinfo[ifile].loaded=1;
  targinfo[ifile].display=1;
  plotstate=getplotstate(DYNAMIC_PLOTS);
  updatetimes();
  fclose(stream);
  IDLE();

}

/* ------------------ gettargetposition ------------------------ */

int gettargetposition(int itarget, float time, float *x, float *y, float *z){
  targpos *tp;
  int i,nsteps;
  float factor;
  float denom;

  if(itarget<0||itarget>=ntargets)return INVISIBLE;
  tp = target_positions+itarget;
  nsteps = tp->nsteps;
  if(nsteps<=1)return INVISIBLE;
  if(time<tp->t[0]||time>tp->t[nsteps-1])return INVISIBLE;
  for(i=0;i<nsteps-1;i++){
    if(tp->t[i]<=time&&time<=tp->t[i+1]){
      denom=tp->t[i+1]-tp->t[i];
      if(denom!=0.0){factor=(time-tp->t[i])/denom;}
      else{factor=1.0;}
      *x=(1.0-factor)*tp->x[i]+factor*tp->x[i+1];
      *y=(1.0-factor)*tp->y[i]+factor*tp->y[i+1];
      *z=(1.0-factor)*tp->z[i]+factor*tp->z[i+1];
      *x = (*x-xbar0)/xyzmaxdiff;
      *y = (*y-ybar0)/xyzmaxdiff;
      *z = (*z-zbar0)/xyzmaxdiff;
      return VISIBLE;
    }
  }
  return INVISIBLE;
}


/* ----------------------- drawTargets ----------------------------- */

void drawTargets(void){
  float time_val;
  float xtarget, ytarget, ztarget;
  int i, j;

  time_val=0.0;
  if(times!=NULL)time_val = times[itime];
  for(j=0;j<ntarg_files;j++){
    if(targinfo[j].loaded==0||targinfo[j].display==0)continue;
    switch (targinfo[j].type){
     case 1:
     for(i=0;i<ntargets;i++){
       if(gettargetposition(i,time_val,&xtarget,&ytarget,&ztarget)==1){
         glPointSize((float)4.0);
         glBegin(GL_POINTS);
         glColor3fv(target_positions[i].rgb);
         glVertex3f(xtarget,ytarget,ztarget);
         glEnd();
       }
     }
     break;
     case 2:
       for(i=0;i<ntargets;i++){
         unsigned char *color;

         j = targtimeslist[itime];
         color = target_positions[i].color;
         xtarget=(target_positions[i].x[j]-xbar0)/xyzmaxdiff;;
         ytarget=(target_positions[i].y[j]-ybar0)/xyzmaxdiff;
         ztarget=(target_positions[i].z[j]-zbar0)/xyzmaxdiff;
         glPointSize((float)4.0);
         glBegin(GL_QUADS);
#define DELTA 0.1f
#define RIGHT 0
#define LEFT  1
#define BACK 2
#define FRONT 3
#define TOP 4
#define BOTTOM 5
         /* left */
         glColor3fv(rgb_full[color[6*j+LEFT]]);
         glVertex3f(xtarget,ytarget,      ztarget      );
         glVertex3f(xtarget,ytarget,      ztarget+DELTA);
         glVertex3f(xtarget,ytarget+DELTA,ztarget+DELTA);
         glVertex3f(xtarget,ytarget+DELTA,ztarget      );

         /* right */
         glColor3fv(rgb_full[color[6*j+RIGHT]]);
         glVertex3f(xtarget+DELTA,ytarget,      ztarget      );
         glVertex3f(xtarget+DELTA,ytarget+DELTA,ztarget      );
         glVertex3f(xtarget+DELTA,ytarget+DELTA,ztarget+DELTA);
         glVertex3f(xtarget+DELTA,ytarget,      ztarget+DELTA);

         /* front */
         glColor3fv(rgb_full[color[6*j+FRONT]]);
         glVertex3f(xtarget,      ytarget,ztarget      );
         glVertex3f(xtarget+DELTA,ytarget,ztarget      );
         glVertex3f(xtarget+DELTA,ytarget,ztarget+DELTA);
         glVertex3f(xtarget,      ytarget,ztarget+DELTA);

         /* back */
         glColor3fv(rgb_full[color[6*j+BACK]]);
         glVertex3f(xtarget,      ytarget+DELTA,ztarget      );
         glVertex3f(xtarget,      ytarget+DELTA,ztarget+DELTA);
         glVertex3f(xtarget+DELTA,ytarget+DELTA,ztarget+DELTA);
         glVertex3f(xtarget+DELTA,ytarget+DELTA,ztarget      );

         /* bottom */
         glColor3fv(rgb_full[color[6*j+BOTTOM]]);
         glVertex3f(xtarget,      ytarget,      ztarget);
         glVertex3f(xtarget,      ytarget+DELTA,ztarget);
         glVertex3f(xtarget+DELTA,ytarget+DELTA,ztarget);
         glVertex3f(xtarget+DELTA,ytarget,      ztarget);

         /* top */
         glColor3fv(rgb_full[color[6*j+TOP]]);
         glVertex3f(xtarget,      ytarget,      ztarget+DELTA);
         glVertex3f(xtarget+DELTA,ytarget,      ztarget+DELTA);
         glVertex3f(xtarget+DELTA,ytarget+DELTA,ztarget+DELTA);
         glVertex3f(xtarget,      ytarget+DELTA,ztarget+DELTA);



         glEnd();
       }
       break;
     default:
       ASSERT(FFALSE);
       break;
    }
  }
}



