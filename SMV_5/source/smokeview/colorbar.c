// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "flowfiles.h"
#include "MALLOC.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

char colorbar_revision[]="$Revision$";

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

  ncolorbars++;
  CheckMemory;
  ResizeMemory((void **)&colorbarinfo,ncolorbars*sizeof(colorbardata));
  cb_from = colorbarinfo + icolorbar;
  CheckMemory;

      // new colorbar

  cb_to=colorbarinfo+ncolorbars-1;

  memcpy(cb_to,cb_from,sizeof(colorbardata));
  strcpy(cb_to->label,"Copy of ");
  strcat(cb_to->label,cb_from->label);
  cb_to->label_ptr=cb_to->label;

  remapcolorbar(cb_to);

}

/* ------------------ drawcolorbarpath ------------------------ */

void drawcolorbarpath(void){
  int i;
  float *rrgb;
  colorbardata *cbi;
  float *rgbleft, *rgbright;

  cbi = colorbarinfo + colorbartype;
  glPointSize(5.0);
  glBegin(GL_POINTS);
  for(i=0;i<255;i++){
    rrgb=cbi->colorbar+3*i;
    glColor3fv(rrgb);
    glVertex3fv(rrgb);
  }
  glEnd();

  glPointSize(10.0);
  glBegin(GL_POINTS);
  for(i=0;i<cbi->nlegs;i++){
    rrgb=cbi->leg_rgb+6*i;
    glColor3fv(rrgb);
    glVertex3fv(rrgb);
  }
#define PLEFT -0.01
#define PRIGHT 1.01

#define PLEFT2 -0.1
#define PRIGHT2 1.1

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

  if(colorbarpoint>=0&&colorbarpoint<cbi->nlegs){
    rgbleft = cbi->leg_rgb+6*colorbarpoint;
    rgbright = rgbleft - 3;

    glPointSize(20.0);
    glBegin(GL_POINTS);
    glColor3fv(rgbleft);
    glVertex3fv(rgbleft);
    if(cbi->splitflag[colorbarpoint]==1){
      glColor3fv(rgbright);
      glVertex3fv(rgbright);
    }
    glEnd();
    if(cbi->splitflag[colorbarpoint]==1){
      output3Text(foregroundcolor,  rgbleft[0]+2*PLEFT, rgbleft[1]+2*PLEFT, rgbleft[2]+2*PLEFT, "Right");
      output3Text(foregroundcolor, rgbright[0]+2*PLEFT,rgbright[1]+2*PLEFT,rgbright[2]+2*PLEFT,"Left");
    }

  }


  {
    float zbot;
    float dzpoint;


    glPointSize(10.0);
    glBegin(GL_POINTS);
    for(i=0;i<cbi->nlegs;i++){
      rrgb = cbi->colorbar+3*cbi->colorbar_index[i];
      dzpoint = (float)cbi->colorbar_index[i]/255.0;
      glColor3fv(rrgb);
      glVertex3f(1.5,0.0,dzpoint);
    }
    glEnd();

    for(i=0;i<cbi->nlegs;i++){
      int ii;
      char cbuff[1024];

      ii = cbi->colorbar_index[i];
      dzpoint = (float)cbi->colorbar_index[i]/255.0;
      sprintf(cbuff,"%i",ii);
      output3Text(foregroundcolor, 1.55,0.0,dzpoint,cbuff);
    }
    if(colorbarpoint>=0&&colorbarpoint<cbi->nlegs){
      glPointSize(20.0);
      glBegin(GL_POINTS);
      rrgb = cbi->colorbar+3*cbi->colorbar_index[colorbarpoint];
      dzpoint = (float)cbi->colorbar_index[colorbarpoint]/255.0;
      glColor3fv(rrgb);
      glVertex3f(1.5,0.0,dzpoint);
      glEnd();
    }

    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      rrgb=cbi->colorbar+3*i;
      glColor3fv(rrgb);
      zbot=(float)i/255.0;
      glVertex3f(1.1,0.0,zbot);
      glVertex3f(1.3,0.0,zbot);
    }
    glEnd();

    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      rrgb=cbi->colorbar+3*i;
      glColor3fv(rrgb);
      zbot=(float)i/255.0;
      glVertex3f(1.3,0.0,zbot);
      glVertex3f(1.1,0.0,zbot);
    }
    glEnd();
    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      rrgb=cbi->colorbar+3*i;
      glColor3fv(rrgb);
      zbot=(float)i/255.0;
      glVertex3f(1.2,-0.1,zbot);
      glVertex3f(1.2, 0.1,zbot);
    }
    glEnd();
    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      rrgb=cbi->colorbar+3*i;
      glColor3fv(rrgb);
      zbot=(float)i/255.0;
      glVertex3f(1.2, 0.1,zbot);
      glVertex3f(1.2,-0.1,zbot);
    }
    glEnd();
  }
}


/* ------------------ adjust_colorbar_splits ------------------------ */

void adjust_colorbar_splits(colorbardata *cbi){
  int i,j;
  int nsplits;

  nsplits=0;
  for(i=0;i<cbi->nlegs;i++){
    if(i==0||i==cbi->nlegs-1||cbi->splitflag[i]==1){
      cbi->splits[nsplits++]=i;
    }
  }
  cbi->nsplits=nsplits;

  for(i=0;i<cbi->nsplits-1;i++){
    int i1, i2;
    unsigned char *ci;

    i1 = cbi->splits[i];
    i2 = cbi->splits[i+1];
    ci = cbi->colorbar_index;
    for(j=i1+1;j<i2;j++){
      ci[j] = ((float)(i2-j)*ci[i1] + (float)(j-i1)*ci[i2])/(float)(i2-i1);
    }
  }
}
/* ------------------ remapcolorbar ------------------------ */

void remapcolorbar(colorbardata *cbi){
  int i;

  CheckMemory;
  for(i=1;i<cbi->nlegs;i++){
    float *rgbleft, *rgbright;

    if(cbi->splitflag[i]==0){
      rgbleft = cbi->leg_rgb+6*i;
      rgbright = cbi->leg_rgb+6*i - 3;

      rgbright[0] = rgbleft[0];
      rgbright[1] = rgbleft[1];
      rgbright[2] = rgbleft[2];
    }
  }
  CheckMemory;

  adjust_colorbar_splits(cbi);

  {
    int nlegs=0,nlegs_seg;

    for(i=0;i<cbi->nlegs-1;i++){
      if(i==cbi->nlegs-2){
        nlegs_seg=256-nlegs;
      }
      else{
        float frac;

        frac=(float)(cbi->colorbar_index[i+1]-cbi->colorbar_index[i])/255.0;
        nlegs_seg=256*frac;
      }
      nlegs += nlegs_seg;
      interpcolor(cbi->leg_rgb+6*i,cbi->leg_rgb+6*i+3,cbi->colorbar+3*(nlegs-nlegs_seg),nlegs_seg,cbi->splitflag[i]);
      CheckMemory;
    }
  }
  CheckMemory;
}

/* ------------------ interpcolor ------------------------ */

void interpcolor(float *col1, float *col2,float *rrgb,int nlegs_seg, int splitflag){

  float dr, dg, db;
  int i;
  int nn;


  // if two adjacent segments are continuous then don't include the last point

  if(splitflag==0){
    nn = nlegs_seg-1;
  }
  else{
    nn = nlegs_seg;
  }
  if(nn<1)nn=1;

  dr=(col2[0]-col1[0])/nn;
  dg=(col2[1]-col1[1])/nn;
  db=(col2[2]-col1[2])/nn;
  CheckMemory;
  for(i=0;i<nlegs_seg;i++){
    nn=3*i;
    rrgb[nn]   = col1[0] + i*dr;
    rrgb[nn+1] = col1[1] + i*dg;
    rrgb[nn+2] = col1[2] + i*db;
  }

}

/* ------------------ initdefaultcolorbars ------------------------ */

void initdefaultcolorbars(void){
  int i;
  float dval;
  int nlegs;
  colorbardata *cbi;
  int ii;

  FREEMEMORY(colorbarinfo);
  ncolorbars=ndefaultcolorbars;
  NewMemory((void **)&colorbarinfo,ncolorbars*sizeof(colorbardata));

  for(i=0;i<ncolorbars;i++){
    cbi = colorbarinfo + i;
    cbi->nlegs=0;
    cbi->legindex=0;
  }
      // rainbow colorbar

  cbi=colorbarinfo;
  strcpy(cbi->label,"Rainbow");
  cbi->label_ptr=cbi->label;
  cbi->nlegs=5;
  cbi->valmin=0.0;
  cbi->valmax=100.0;

  cbi->leg_rgb[0]=0.0;
  cbi->leg_rgb[1]=0.0;
  cbi->leg_rgb[2]=1.0;

  cbi->leg_rgb[3]=0.0;
  cbi->leg_rgb[4]=1.0;
  cbi->leg_rgb[5]=1.0;

  cbi->leg_rgb[6]=0.0;
  cbi->leg_rgb[7]=1.0;
  cbi->leg_rgb[8]=1.0;

  cbi->leg_rgb[9]=0.0;
  cbi->leg_rgb[10]=1.0;
  cbi->leg_rgb[11]=0.0;

  cbi->leg_rgb[12]=0.0;
  cbi->leg_rgb[13]=1.0;
  cbi->leg_rgb[14]=0.0;

  cbi->leg_rgb[15]=1.0;
  cbi->leg_rgb[16]=1.0;
  cbi->leg_rgb[17]=0.0;

  cbi->leg_rgb[18]=1.0;
  cbi->leg_rgb[19]=1.0;
  cbi->leg_rgb[20]=0.0;

  cbi->leg_rgb[21]=1.0;
  cbi->leg_rgb[22]=0.0;
  cbi->leg_rgb[23]=0.0;

  cbi->leg_rgb[24]=1.0;
  cbi->leg_rgb[25]=0.0;
  cbi->leg_rgb[26]=0.0;

  cbi->leg_rgb[27]=1.0;
  cbi->leg_rgb[28]=0.0;
  cbi->leg_rgb[29]=0.0;

  cbi->splitflag[0]=0;
  cbi->splitflag[1]=0;
  cbi->splitflag[2]=0;
  cbi->splitflag[3]=0;
  cbi->splitflag[4]=0;

  nlegs = cbi->nlegs-1;
  cbi->colorbar_index[0]=0;
  cbi->colorbar_index[1]=64;
  cbi->colorbar_index[2]=128;
  cbi->colorbar_index[3]=192;
  cbi->colorbar_index[4]=255;
  dval = (cb_valmax-cb_valmin)/(float)nlegs;
  cbi->legvals[0]=cb_valmin;
  cbi->legvals[1]=cb_valmin+dval;
  cbi->legvals[2]=cb_valmin+2*dval;
  cbi->legvals[3]=cb_valmin+3*dval;
  cbi->legvals[4]=cb_valmax;

      // b&w colorbar

  cbi=colorbarinfo+1;
  strcpy(cbi->label,"Black and White");
  cbi->label_ptr=cbi->label;
  cbi->valmin=0.0;
  cbi->valmax=100.0;

  cbi->nlegs=3;

  cbi->leg_rgb[0]=1.0;
  cbi->leg_rgb[1]=1.0;
  cbi->leg_rgb[2]=1.0;

  cbi->leg_rgb[3]=0.5;
  cbi->leg_rgb[4]=0.5;
  cbi->leg_rgb[5]=0.5;

  cbi->leg_rgb[6]=0.5;
  cbi->leg_rgb[7]=0.5;
  cbi->leg_rgb[8]=0.5;

  cbi->leg_rgb[9] =0.0;
  cbi->leg_rgb[10]=0.0;
  cbi->leg_rgb[11]=0.0;

  cbi->leg_rgb[12] =0.0;
  cbi->leg_rgb[13]=0.0;
  cbi->leg_rgb[14]=0.0;

  cbi->leg_rgb[15] =0.0;
  cbi->leg_rgb[16]=0.0;
  cbi->leg_rgb[17]=0.0;

  cbi->splitflag[0]=0;
  cbi->splitflag[1]=0;
  cbi->splitflag[2]=0;

  nlegs = cbi->nlegs-1;
  cbi->colorbar_index[0]=0;
  cbi->colorbar_index[1]=128;
  cbi->colorbar_index[2]=255;
  cbi->legvals[0]=cb_valmin;
  cbi->legvals[1]=cb_valmax;

  remapcolorbar(colorbarinfo);
  remapcolorbar(colorbarinfo+1);
}


