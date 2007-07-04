#include "options.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "flowfiles.h"
#include "MALLOC.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"


/* ------------------ setColorbarClipPlanes ------------------------ */

void setColorbarClipPlanes(int flag){
  static GLdouble clipplane_x[4], clipplane_y[4], clipplane_z[4];
  static GLdouble clipplane_X[4], clipplane_Y[4], clipplane_Z[4];

  if(flag==1){
      clipplane_x[0]=1.0;
      clipplane_x[1]=0.0;
      clipplane_x[2]=0.0;
      clipplane_x[3]=-2.0;
      glClipPlane(GL_CLIP_PLANE0,clipplane_x);
      glEnable(GL_CLIP_PLANE0);

      clipplane_X[0]=-1.0;
      clipplane_X[1]=0.0;
      clipplane_X[2]=0.0;
      clipplane_X[3]=2.0;
      glClipPlane(GL_CLIP_PLANE3,clipplane_X);
      glEnable(GL_CLIP_PLANE3);

      clipplane_y[0]=0.0;
      clipplane_y[1]=1.0;
      clipplane_y[2]=0.0;
      clipplane_y[3]=-2.0;
      glClipPlane(GL_CLIP_PLANE1,clipplane_y);
      glEnable(GL_CLIP_PLANE1);

      clipplane_Y[0]=0.0;
      clipplane_Y[1]=-1.0;
      clipplane_Y[2]=0.0;
      clipplane_Y[3]=2.0;
      glClipPlane(GL_CLIP_PLANE4,clipplane_Y);
      glEnable(GL_CLIP_PLANE4);

      clipplane_z[0]=0.0;
      clipplane_z[1]=0.0;
      clipplane_z[2]=1.0;
      clipplane_z[3]=-2.0;
      glClipPlane(GL_CLIP_PLANE2,clipplane_z);
      glEnable(GL_CLIP_PLANE2);

      clipplane_Z[0]=0.0;
      clipplane_Z[1]=0.0;
      clipplane_Z[2]=-1.0;
      clipplane_Z[3]=2.0;
      glClipPlane(GL_CLIP_PLANE5,clipplane_Z);
      glEnable(GL_CLIP_PLANE5);
  }
  else{
    glDisable(GL_CLIP_PLANE0);
    glDisable(GL_CLIP_PLANE1);
    glDisable(GL_CLIP_PLANE2);
    glDisable(GL_CLIP_PLANE3);
    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);
  }
}

/* ------------------ addcolorbar ------------------------ */

void addcolorbar(int icolorbar){
  colorbardata *cb_to, *cb_from;
  int i;

  ncolorbars++;
  CheckMemory;
  ResizeMemory((void **)&colorbarinfo,ncolorbars*sizeof(colorbardata));
  cb_from = colorbarinfo + icolorbar;
  CheckMemory;

      // new colorbar

  cb_to=colorbarinfo+ncolorbars-1;
  NewMemory((void **)&cb_to->rgb,3*MAXRGB*sizeof(float));
  strcpy(cb_to->label,"Copy of ");
  strcat(cb_to->label,cb_from->label);
  cb_to->label_ptr=cb_to->label;
  cb_to->npoints=cb_from->npoints;
  cb_to->pointindex=cb_from->pointindex;
  NewMemory((void **)&cb_to->flegs,cb_to->npoints*sizeof(float));
  NewMemory((void **)&cb_to->c_vals,cb_to->npoints*sizeof(float));
  NewMemory((void **)&cb_to->rgbnodes,6*cb_to->npoints*sizeof(float));
  NewMemory((void **)&cb_to->jumpflag,cb_to->npoints*sizeof(int));

  for(i=0;i<6*cb_to->npoints;i++){
    cb_to->rgbnodes[i]=cb_from->rgbnodes[i];
  }
  for(i=0;i<cb_to->npoints;i++){
    cb_to->c_vals[i]=cb_from->c_vals[i];
    cb_to->jumpflag[i]=cb_from->jumpflag[i];
    cb_to->flegs[i]=cb_from->flegs[i];
  }

  remapcolorbar(cb_to);

}

/* ------------------ drawcolorbarpath ------------------------ */

void drawcolorbarpath(void){
  int i;
  float *rrgb,*rrgb2;
  colorbardata *cbi;
  float *rgbleft, *rgbright;

  cbi = colorbarinfo + colorbartype;

  glPointSize(5.0);
  glBegin(GL_POINTS);
  for(i=0;i<255;i++){
    rrgb=cbi->rgb+3*i;
    glColor3fv(rrgb);
    glVertex3fv(rrgb);
  }
  glEnd();

  glPointSize(10.0);
  glBegin(GL_POINTS);
  for(i=0;i<cbi->npoints;i++){
    rrgb=cbi->rgbnodes+6*i;
    rrgb2=rrgb+3;
    glColor3fv(rrgb);
    glVertex3fv(rrgb);
    /*
    if(i!=cbi->npoints-1){
      glColor3fv(rrgb2);
      glVertex3fv(rrgb2);
    }
    */
  }
#define PLEFT -0.01
#define PRIGHT 1.01

#define PLEFT2 -0.1
#define PRIGHT2 1.1

//  glColor3fv(foregroundcolor);
//  glVertex3f(PLEFT,PLEFT,PLEFT);
//  glVertex3f(PLEFT,PRIGHT,PLEFT);
//  glVertex3f(PRIGHT,PLEFT,PLEFT);
//  glVertex3f(PRIGHT,PRIGHT,PLEFT);
//  glVertex3f(PLEFT,PLEFT,PRIGHT);
//  glVertex3f(PLEFT,PRIGHT,PRIGHT);
//  glVertex3f(PRIGHT,PLEFT,PRIGHT);
//  glVertex3f(PRIGHT,PRIGHT,PRIGHT);

  glEnd();
  // draw rgb color axese

  glLineWidth(5.0);
  glBegin(GL_LINES);
  glColor3f(1.0,0.0,0.0);
  glVertex3f( PLEFT2,PLEFT2,PLEFT2);
  glVertex3f(PRIGHT2,PLEFT2,PLEFT2);

  glColor3f(0.0,1.0,0.0);
  glVertex3f(PLEFT2, PLEFT2,PLEFT2);
  glVertex3f(PLEFT2,PRIGHT2,PLEFT2);

  glColor3f(0.0,0.0,1.0);
  glVertex3f(PLEFT2,PLEFT2, PLEFT2);
  glVertex3f(PLEFT2,PLEFT2,PRIGHT2);

  glEnd();

  if(colorbarpoint>=0&&colorbarpoint<cbi->npoints){
    rgbleft = cbi->rgbnodes+6*colorbarpoint;
    rgbright = rgbleft - 3;

    glPointSize(20.0);
    glBegin(GL_POINTS);
    glColor3fv(rgbleft);
    glVertex3fv(rgbleft);
    if(cbi->jumpflag[colorbarpoint]==1){
      glColor3fv(rgbright);
      glVertex3fv(rgbright);
    }
    glEnd();
    if(cbi->jumpflag[colorbarpoint]==1){
      output3Text(foregroundcolor,  rgbleft[0]+2*PLEFT, rgbleft[1]+2*PLEFT, rgbleft[2]+2*PLEFT, "Right");
      output3Text(foregroundcolor, rgbright[0]+2*PLEFT,rgbright[1]+2*PLEFT,rgbright[2]+2*PLEFT,"Left");
    }

  }


  {
    float zbot;

    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      rrgb=cbi->rgb+3*i;
      glColor3fv(rrgb);
      zbot=(float)i/255.0;
      glVertex3f(1.1,0.0,zbot);
      glVertex3f(1.3,0.0,zbot);
    }
    glEnd();
    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      rrgb=cbi->rgb+3*i;
      glColor3fv(rrgb);
      zbot=(float)i/255.0;
      glVertex3f(1.3,0.0,zbot);
      glVertex3f(1.1,0.0,zbot);
    }
    glEnd();
    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      rrgb=cbi->rgb+3*i;
      glColor3fv(rrgb);
      zbot=(float)i/255.0;
      glVertex3f(1.2,-0.1,zbot);
      glVertex3f(1.2, 0.1,zbot);
    }
    glEnd();
    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      rrgb=cbi->rgb+3*i;
      glColor3fv(rrgb);
      zbot=(float)i/255.0;
      glVertex3f(1.2, 0.1,zbot);
      glVertex3f(1.2,-0.1,zbot);
    }
    glEnd();
  }
}

/* ------------------ freecolorbars ------------------------ */

void freecolorbars(void){
  int i;

  if(ncolorbars>0&&colorbarinfo!=NULL){
    {
      colorbardata *cbi;

      for(i=0;i<ncolorbars;i++){
        cbi=colorbarinfo+i;
        freecolorbar(cbi);
      }
    }
    FREEMEMORY(colorbarinfo);
    ncolorbars=0;
  }
  if(ncolorbarsini>0){
    {
      colorbardata *cbi;

      for(i=0;i<ncolorbarsini;i++){
        cbi=colorbariniinfo+i;
        freecolorbar(cbi);
      }
    }
    FREEMEMORY(colorbariniinfo);
    ncolorbarsini=0;
  }
}

/* ------------------ remapcolorbar ------------------------ */

void remapcolorbar(colorbardata *cbi){
  int i;
  float sum;

  sum = 0.0;
  for(i=0;i<cbi->npoints-1;i++){
    sum+=cbi->flegs[i];
  }
  for(i=1;i<cbi->npoints;i++){
    float *rgbleft, *rgbright;
    float *rrr;

    rrr = cbi->rgbnodes+6*i;
    if(cbi->jumpflag[i]==0){
       rgbleft = cbi->rgbnodes+6*i;
      rgbright = cbi->rgbnodes+6*i - 3;

      rgbright[0] = rgbleft[0];
      rgbright[1] = rgbleft[1];
      rgbright[2] = rgbleft[2];

    }
  }

  if(sum==0.0)sum=1.0;
  for(i=0;i<cbi->npoints-1;i++){
    cbi->flegs[i]/=sum;
  }
  {
    int npoints=0,npoints_seg;

    for(i=0;i<cbi->npoints-1;i++){
      if(i==cbi->npoints-2){
        npoints_seg=256-npoints;
      }
      else{
        npoints_seg=256*cbi->flegs[i];
      }
      npoints += npoints_seg;
      interpcolor(cbi->rgbnodes+6*i,cbi->rgbnodes+6*i+3,cbi->rgb+3*(npoints-npoints_seg),npoints_seg,cbi->jumpflag[i]);
    }
  }
}

/* ------------------ freecolorbar ------------------------ */

void freecolorbar(colorbardata *cbi){
  if(cbi==NULL)return;
  FREEMEMORY(cbi->rgb);
  FREEMEMORY(cbi->rgbnodes);
  FREEMEMORY(cbi->flegs);
  FREEMEMORY(cbi->c_vals);
  FREEMEMORY(cbi->jumpflag);
  cbi->npoints=0;
}

/* ------------------ interpcolor ------------------------ */

void interpcolor(float *col1, float *col2,float *rrgb,int npoints_seg, int jumpflag){

  float dr, dg, db;
  int i;
  int nn;


  // if two adjacent segments are continuous then don't include the last point

  if(jumpflag==0){
    nn = npoints_seg-1;
  }
  else{
    nn = npoints_seg;
  }
  if(nn<1)nn=1;

  dr=(col2[0]-col1[0])/nn;
  dg=(col2[1]-col1[1])/nn;
  db=(col2[2]-col1[2])/nn;

  for(i=0;i<npoints_seg;i++){
    nn=3*i;
    rrgb[nn]   = col1[0] + i*dr;
    rrgb[nn+1] = col1[1] + i*dg;
    rrgb[nn+2] = col1[2] + i*db;
  }

}

/* ------------------ initdefaultcolorbars ------------------------ */

void initdefaultcolorbars(void){
  int i;

  if(colorbarinfo==NULL){
    if(ncolorbars==0)ncolorbars=ndefaultcolorbars;
    NewMemory((void **)&colorbarinfo,ncolorbars*sizeof(colorbardata));

    {
      colorbardata *cbi;
      int ii;

      for(i=0;i<ncolorbars;i++){
        cbi = colorbarinfo + i;
        cbi->rgb=NULL;
        cbi->flegs=NULL;
        cbi->c_vals=NULL;
        cbi->rgbnodes=NULL;
        cbi->jumpflag=NULL;
        cbi->npoints=0;
        cbi->pointindex=0;
      }

      // original colorbar

      cbi=colorbarinfo;
      NewMemory((void **)&cbi->rgb,3*MAXRGB*sizeof(float));
      strcpy(cbi->label,"Original");
      cbi->label_ptr=cbi->label;
      cbi->npoints=nrgb;
      NewMemory((void **)&cbi->flegs,cbi->npoints*sizeof(float));
      NewMemory((void **)&cbi->c_vals,cbi->npoints*sizeof(float));
      NewMemory((void **)&cbi->rgbnodes,6*cbi->npoints*sizeof(float));
      NewMemory((void **)&cbi->jumpflag,cbi->npoints*sizeof(float));
      for(i=0;i<cbi->npoints;i++){
        ii = 6*i;
        cbi->rgbnodes[ii]  =rgb[i][0];
        cbi->rgbnodes[ii+1]=rgb[i][1];
        cbi->rgbnodes[ii+2]=rgb[i][2];
        if(i!=cbi->npoints-1){
          cbi->rgbnodes[ii+3]=rgb[i+1][0];
          cbi->rgbnodes[ii+4]=rgb[i+1][1];
          cbi->rgbnodes[ii+5]=rgb[i+1][2];
          cbi->flegs[i]=1.0/(float)(cbi->npoints-1);
        }
        cbi->jumpflag[i]=0;
      }

      // rainbow colorbar

      cbi=colorbarinfo+1;
      NewMemory((void **)&cbi->rgb,3*MAXRGB*sizeof(float));
      strcpy(cbi->label,"Rainbow");
      cbi->label_ptr=cbi->label;
      cbi->npoints=5;
      NewMemory((void **)&cbi->flegs,cbi->npoints*sizeof(float));
      NewMemory((void **)&cbi->c_vals,cbi->npoints*sizeof(float));
      NewMemory((void **)&cbi->rgbnodes,6*cbi->npoints*sizeof(float));
      NewMemory((void **)&cbi->jumpflag,cbi->npoints*sizeof(float));

      cbi->rgbnodes[0]=0.0;
      cbi->rgbnodes[1]=0.0;
      cbi->rgbnodes[2]=1.0;

      cbi->rgbnodes[3]=0.0;
      cbi->rgbnodes[4]=1.0;
      cbi->rgbnodes[5]=1.0;

      cbi->rgbnodes[6]=0.0;
      cbi->rgbnodes[7]=1.0;
      cbi->rgbnodes[8]=1.0;

      cbi->rgbnodes[9]=0.0;
      cbi->rgbnodes[10]=1.0;
      cbi->rgbnodes[11]=0.0;

      cbi->rgbnodes[12]=0.0;
      cbi->rgbnodes[13]=1.0;
      cbi->rgbnodes[14]=0.0;

      cbi->rgbnodes[15]=1.0;
      cbi->rgbnodes[16]=1.0;
      cbi->rgbnodes[17]=0.0;

      cbi->rgbnodes[18]=1.0;
      cbi->rgbnodes[19]=1.0;
      cbi->rgbnodes[20]=0.0;

      cbi->rgbnodes[21]=1.0;
      cbi->rgbnodes[22]=0.0;
      cbi->rgbnodes[23]=0.0;

      cbi->rgbnodes[24]=1.0;
      cbi->rgbnodes[25]=0.0;
      cbi->rgbnodes[26]=0.0;

      cbi->rgbnodes[27]=1.0;
      cbi->rgbnodes[28]=0.0;
      cbi->rgbnodes[29]=0.0;


      cbi->jumpflag[0]=0;
      cbi->jumpflag[1]=0;
      cbi->jumpflag[2]=0;
      cbi->jumpflag[3]=0;

      {
        int nlegs;

        nlegs = cbi->npoints-1;
        cbi->flegs[0]=1.0/(float)nlegs;
        cbi->flegs[1]=1.0/(float)nlegs;
        cbi->flegs[2]=1.0/(float)nlegs;
        cbi->flegs[3]=1.0/(float)nlegs;
      }

      // b&w colorbar

      cbi=colorbarinfo+2;
      NewMemory((void **)&cbi->rgb,3*MAXRGB*sizeof(float));
      strcpy(cbi->label,"Black and White");
      cbi->label_ptr=cbi->label;

      cbi->npoints=3;
      NewMemory((void **)&cbi->flegs,cbi->npoints*sizeof(float));
      NewMemory((void **)&cbi->c_vals,cbi->npoints*sizeof(float));
      NewMemory((void **)&cbi->rgbnodes,6*cbi->npoints*sizeof(float));
      NewMemory((void **)&cbi->jumpflag,cbi->npoints*sizeof(float));

      cbi->rgbnodes[0]=1.0;
      cbi->rgbnodes[1]=1.0;
      cbi->rgbnodes[2]=1.0;

      cbi->rgbnodes[3]=0.5;
      cbi->rgbnodes[4]=0.5;
      cbi->rgbnodes[5]=0.5;

      cbi->rgbnodes[6]=0.5;
      cbi->rgbnodes[7]=0.5;
      cbi->rgbnodes[8]=0.5;

      cbi->rgbnodes[9] =0.0;
      cbi->rgbnodes[10]=0.0;
      cbi->rgbnodes[11]=0.0;

      cbi->rgbnodes[12] =0.0;
      cbi->rgbnodes[13]=0.0;
      cbi->rgbnodes[14]=0.0;

      cbi->rgbnodes[15] =0.0;
      cbi->rgbnodes[16]=0.0;
      cbi->rgbnodes[17]=0.0;

      cbi->jumpflag[0]=0;
      cbi->jumpflag[1]=0;

      {
        int nlegs;

        nlegs = cbi->npoints-1;
        cbi->flegs[0]=1.0/(float)nlegs;
        cbi->flegs[1]=1.0/(float)nlegs;
      }

    }
  }
  remapcolorbar(colorbarinfo);
  remapcolorbar(colorbarinfo+1);
  remapcolorbar(colorbarinfo+2);

}

/* ------------------ initcolorbar ------------------------ */

void initcolorbars(void){

  FREEMEMORY(colorbarinfo);

//  NewMemory((void **)&colorbarinfo,ncolorbars*sizeof(colorbardata));

  initdefaultcolorbars();

  // add code here to init colorbars read in from .ini file

}

/* ------------------ ResizeColorbar ------------------------ */

void ResizeColorbar(colorbardata *cbi, int n){
  cbi->npoints=n;
  ResizeMemory((void **)&cbi->flegs,n*sizeof(float));
  ResizeMemory((void **)&cbi->c_vals,n*sizeof(float));
  ResizeMemory((void **)&cbi->rgbnodes,6*n*sizeof(float));
  ResizeMemory((void **)&cbi->jumpflag,n*sizeof(int));
}
