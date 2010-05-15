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
  unsigned char *rrgb;
  colorbardata *cbi;
  unsigned char *rgbleft;

  cbi = colorbarinfo + colorbartype;
  glPointSize(5.0);
  glBegin(GL_POINTS);
  for(i=0;i<255;i++){
    float *rgb;

    rgb=cbi->colorbar+3*i;
    glColor3fv(rgb);
    glVertex3fv(rgb);
  }
  glEnd();

  glPointSize(10.0);
  glBegin(GL_POINTS);
  for(i=0;i<cbi->nnodes;i++){
    rrgb=cbi->rgb_node+3*i;
    glColor3ubv(rrgb);
    glVertex3f(rrgb[0]/255.0,rrgb[1]/255.0,rrgb[2]/255.0);
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

  if(colorbarpoint>=0&&colorbarpoint<cbi->nnodes){
    rgbleft = cbi->rgb_node+3*colorbarpoint;

    glPointSize(20.0);
    glBegin(GL_POINTS);
    glColor3ubv(rgbleft);
    glVertex3f(rgbleft[0]/255.0,rgbleft[1]/255.0,rgbleft[2]/255.0);
    glEnd();
  }


  {
    float zbot;
    float dzpoint;

    glPointSize(10.0);
    glBegin(GL_POINTS);
    for(i=0;i<cbi->nnodes;i++){
      float *rgb;

      rgb = cbi->colorbar+3*cbi->index_node[i];
      dzpoint = (float)cbi->index_node[i]/255.0;
      glColor3fv(rgb);
      glVertex3f(1.5,0.0,dzpoint);
    }
    glEnd();

    for(i=0;i<cbi->nnodes;i++){
      int ii;
      char cbuff[1024];

      ii = cbi->index_node[i];
      dzpoint = (float)cbi->index_node[i]/255.0;
      sprintf(cbuff,"%i",ii);
      output3Text(foregroundcolor, 1.55,0.0,dzpoint,cbuff);
    }
    if(colorbarpoint>=0&&colorbarpoint<cbi->nnodes){
      float *rgb;

      glPointSize(20.0);
      glBegin(GL_POINTS);
      rgb = cbi->colorbar+3*cbi->index_node[colorbarpoint];
      dzpoint = (float)cbi->index_node[colorbarpoint]/255.0;
      glColor3fv(rgb);
      glVertex3f(1.5,0.0,dzpoint);
      glEnd();
    }

    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      float *rgb;

      rgb=cbi->colorbar+3*i;
      glColor3fv(rgb);
      zbot=(float)i/255.0;
      glVertex3f(1.1,0.0,zbot);
      glVertex3f(1.3,0.0,zbot);
    }
    glEnd();

    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      float *rgb;

      rgb=cbi->colorbar+3*i;
      glColor3fv(rgb);
      zbot=(float)i/255.0;
      glVertex3f(1.3,0.0,zbot);
      glVertex3f(1.1,0.0,zbot);
    }
    glEnd();
    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      float *rgb;

      rgb=cbi->colorbar+3*i;
      glColor3fv(rgb);
      zbot=(float)i/255.0;
      glVertex3f(1.2,-0.1,zbot);
      glVertex3f(1.2, 0.1,zbot);
    }
    glEnd();
    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      float *rgb;

      rgb=cbi->colorbar+3*i;
      glColor3fv(rgb);
      zbot=(float)i/255.0;
      glVertex3f(1.2, 0.1,zbot);
      glVertex3f(1.2,-0.1,zbot);
    }
    glEnd();
  }
}

/* ------------------ remapcolorbar ------------------------ */

void remapcolorbar(colorbardata *cbi){
  int i,j,i1,i2;
  float factor,*colorbar;
  unsigned char *rgb_node;
  unsigned char *alpha;

  CheckMemory;
  colorbar=cbi->colorbar;
  rgb_node=cbi->rgb_node;
  alpha=cbi->alpha;
  for(i=0;i<cbi->index_node[0];i++){
    colorbar[0+3*i]=rgb_node[0]/255.0;
    colorbar[1+3*i]=rgb_node[1]/255.0;
    colorbar[2+3*i]=rgb_node[2]/255.0;
    if(rgb_node[0]==0&&rgb_node[1]==1&&rgb_node[2]==2){
      alpha[i]=0;
    }
    else{
      alpha[i]=255;
    }
  }
  for(i=0;i<cbi->nnodes-1;i++){
    i1 = cbi->index_node[i];
    i2 = cbi->index_node[i+1];
    if(i2==i1)continue;
    rgb_node = cbi->rgb_node+3*i;
    for(j=i1;j<i2;j++){
      factor = (float)(j-i1)/(float)(i2-i1);
      colorbar[0+3*j]=((1.0-factor)*rgb_node[0]+factor*rgb_node[3])/255.0;
      colorbar[1+3*j]=((1.0-factor)*rgb_node[1]+factor*rgb_node[4])/255.0;
      colorbar[2+3*j]=((1.0-factor)*rgb_node[2]+factor*rgb_node[5])/255.0;
      if(
        rgb_node[0]==0&&rgb_node[1]==1&&rgb_node[2]==2&&
        rgb_node[3]==0&&rgb_node[4]==1&&rgb_node[5]==2
        ){
        alpha[j]=0;
      }
      else{
        alpha[j]=255;
      }
    }
  }
  rgb_node = cbi->rgb_node+3*(cbi->nnodes-1);
  for(i=cbi->index_node[cbi->nnodes-1];i<256;i++){
    colorbar[0+3*i]=rgb_node[0]/255.0;
    colorbar[1+3*i]=rgb_node[1]/255.0;
    colorbar[2+3*i]=rgb_node[2]/255.0;
    if(rgb_node[0]==0&&rgb_node[1]==1&&rgb_node[2]==2){
      alpha[i]=0;
    }
    else{
      alpha[i]=255;
    }
  }
  if(show_extremedata==1){
    colorbar[0]=rgb_below_min[0];
    colorbar[1]=rgb_below_min[1];
    colorbar[2]=rgb_below_min[2];
    colorbar[0+3*255]=rgb_above_max[0];
    colorbar[1+3*255]=rgb_above_max[1];
    colorbar[2+3*255]=rgb_above_max[2];
  }
  CheckMemory;
}

/* ------------------ initdefaultcolorbars ------------------------ */

void initdefaultcolorbars(void){
  int i;
  float dval;
  int nlegs;
  colorbardata *cbi;
  int ii;

  ndefaultcolorbars=7;
  
  FREEMEMORY(colorbarinfo);
  ncolorbars=ndefaultcolorbars;
  NewMemory((void **)&colorbarinfo,ncolorbars*sizeof(colorbardata));

  // rainbow colorbar

  cbi=colorbarinfo;


  strcpy(cbi->label,"Rainbow");
  cbi->label_ptr=cbi->label;
  cbi->nnodes=5;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=0;
  cbi->rgb_node[1]=0;
  cbi->rgb_node[2]=255;

  cbi->index_node[1]=64;
  cbi->rgb_node[3]=0;
  cbi->rgb_node[4]=255;
  cbi->rgb_node[5]=255;

  cbi->index_node[2]=128;
  cbi->rgb_node[6]=0;
  cbi->rgb_node[7]=255;
  cbi->rgb_node[8]=0;

  cbi->index_node[3]=192;
  cbi->rgb_node[9]=255;
  cbi->rgb_node[10]=255;
  cbi->rgb_node[11]=0;

  cbi->index_node[4]=255;
  cbi->rgb_node[12]=255;
  cbi->rgb_node[13]=0;
  cbi->rgb_node[14]=0;
  cbi++;

  // jet colorbar

  strcpy(cbi->label,"Rainbow 2");
  cbi->label_ptr=cbi->label;
  cbi->nnodes=6;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=0;
  cbi->rgb_node[1]=0;
  cbi->rgb_node[2]=143;

  cbi->index_node[1]=32;
  cbi->rgb_node[3]=0;
  cbi->rgb_node[4]=0;
  cbi->rgb_node[5]=255;

  cbi->index_node[2]=96;
  cbi->rgb_node[6]=0;
  cbi->rgb_node[7]=255;
  cbi->rgb_node[8]=255;

  cbi->index_node[3]=160;
  cbi->rgb_node[9]=255;
  cbi->rgb_node[10]=255;
  cbi->rgb_node[11]=0;

  cbi->index_node[4]=234;
  cbi->rgb_node[12]=255;
  cbi->rgb_node[13]=0;
  cbi->rgb_node[14]=0;

  cbi->index_node[5]=255;
  cbi->rgb_node[15]=128;
  cbi->rgb_node[16]=0;
  cbi->rgb_node[17]=0;
  cbi++;
  
  // yellow/red

  strcpy(cbi->label,"yellow->red");
  cbi->label_ptr=cbi->label;
  cbi->nnodes=2;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=255;
  cbi->rgb_node[1]=255;
  cbi->rgb_node[2]=0;

  cbi->index_node[1]=255;
  cbi->rgb_node[3]=255;
  cbi->rgb_node[4]=0;
  cbi->rgb_node[5]=0;
  cbi++;
  
  // blue/green/red

  strcpy(cbi->label,"blue->green->red");
  cbi->label_ptr=cbi->label;
  cbi->nnodes=3;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=0;
  cbi->rgb_node[1]=0;
  cbi->rgb_node[2]=255;

  cbi->index_node[1]=128;
  cbi->rgb_node[3]=0;
  cbi->rgb_node[4]=255;
  cbi->rgb_node[5]=0;

  cbi->index_node[2]=255;
  cbi->rgb_node[6]=255;
  cbi->rgb_node[7]=0;
  cbi->rgb_node[8]=0;
  cbi++;

  // blue->red split

  strcpy(cbi->label,"blue->red split");
  cbi->label_ptr=cbi->label;

  cbi->nnodes=4;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=0;
  cbi->rgb_node[1]=0;
  cbi->rgb_node[2]=255;

  cbi->index_node[1]=128;
  cbi->rgb_node[3]=0;
  cbi->rgb_node[4]=255;
  cbi->rgb_node[5]=255;

  cbi->index_node[2]=128;
  cbi->rgb_node[6]=255;
  cbi->rgb_node[7]=255;
  cbi->rgb_node[8]=0;

  cbi->index_node[3]=255;
  cbi->rgb_node[9]=255;
  cbi->rgb_node[10]=0;
  cbi->rgb_node[11]=0;
  cbi++;

  // blue->red split

  strcpy(cbi->label,"red/black");
  cbi->label_ptr=cbi->label;

  cbi->nnodes=6;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=255;
  cbi->rgb_node[1]=255;
  cbi->rgb_node[2]=255;

  cbi->index_node[1]=120;
  cbi->rgb_node[3]=255;
  cbi->rgb_node[4]=255;
  cbi->rgb_node[5]=255;

  cbi->index_node[2]=120;
  cbi->rgb_node[6]=255;
  cbi->rgb_node[7]=0;
  cbi->rgb_node[8]=0;

  cbi->index_node[3]=136;
  cbi->rgb_node[9]=255;
  cbi->rgb_node[10]=0;
  cbi->rgb_node[11]=0;

  cbi->index_node[4]=136;
  cbi->rgb_node[9]=64;
  cbi->rgb_node[10]=64;
  cbi->rgb_node[11]=64;

  cbi->index_node[5]=255;
  cbi->rgb_node[9]=64;
  cbi->rgb_node[10]=64;
  cbi->rgb_node[11]=64;
  cbi++;

  // b&w colorbar

  bw_colorbar_index = cbi - colorbarinfo;
  strcpy(cbi->label,"white->black");
  cbi->label_ptr=cbi->label;

  cbi->nnodes=2;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=255;
  cbi->rgb_node[1]=255;
  cbi->rgb_node[2]=255;

  cbi->index_node[1]=255;
  cbi->rgb_node[3] =0;
  cbi->rgb_node[4]=0;
  cbi->rgb_node[5]=0;
  cbi++;

// construct colormaps from color node info

  rgb_below_min[0]=255-64;
  rgb_below_min[1]=255-64;
  rgb_below_min[2]=255-64;
  rgb_above_max[0]=0;
  rgb_above_max[1]=0;
  rgb_above_max[2]=0;

  for(i=0;i<ndefaultcolorbars;i++){
    unsigned char *rgb_min, *rgb_max;

    cbi = colorbarinfo + i;

    remapcolorbar(cbi);
    update_colorbar_splits(cbi);
  }
}

/* ------------------ update_colorbar_splits ------------------------ */

void update_colorbar_splits(colorbardata *cbi){
  int i;

  cbi->nsplits=0;
  for(i=1;i<cbi->nnodes;i++){
    if(cbi->index_node[i]==cbi->index_node[i-1]){
      cbi->splits[cbi->nsplits]=i;
      cbi->nsplits++;
    }
  }
}


