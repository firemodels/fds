// $Date$ 
// $Revision$
// $Author$

// svn revision character string
char colorbar_revision[]="$Revision$";

#include "options.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "string_util.h"
#include "smokeviewvars.h"

#define DYFONT (-0.5)

/* ------------------ UpdateTimeLabels ------------------------ */

void UpdateTimeLabels(void){
  float time0;
  int hour, min, sec,sec10;
 
  time0 = timeoffset;
  if(global_times!=NULL)time0 = timeoffset + global_times[itimes];
  if(vishmsTimelabel==1){
    hour = time0/3600;
    min = (int)(time0/60.0 - 60.0*hour);
    sec10 = (int)(10*(time0 -  60.0*min - 3600.0*hour));
    sec = sec10/10;
    sec10 = sec10 - 10*sec;
    sprintf(timelabel,"  %i:%.2i:%.2i.%i",hour,min,sec,sec10);
  }
  else{
    float dt;
    char timeval[30], *timevalptr;

    if(nglobal_times>1){
      dt=global_times[1]-global_times[0];
    }
    else{
      dt=0.0;
    }
    if(dt<0.0)dt=-dt;
    timevalptr=time2timelabel(time0,dt,timeval);
    strcpy(timelabel,"Time: ");
    strcat(timelabel,timevalptr);
  }
  sprintf(framelabel,"Frame: %i",itimes);
  if(hrrinfo!=NULL&&hrrinfo->display==1&&hrrinfo->loaded==1){
    float hrr;

    hrr = hrrinfo->hrrval[hrrinfo->itime];
    if(hrr<1.0){
      sprintf(hrrinfo->hrrlabel,"HRR: %4.1f",hrr*1000.0);
    }
    else if(hrr>1000.0){
      sprintf(hrrinfo->hrrlabel,"HRR: %4.1f MW",hrr/1000.0);
    }
    else{
      sprintf(hrrinfo->hrrlabel,"HRR: %4.1f kW",hrr);
    }
  }
}

/* ------------------ drawTimeBar ------------------------ */

void drawTimeBar(void){
  float xleft=.175f, xright=1.0f, ybot=0.10f, ytop=.35f, xxright;

  glDisable(GL_LIGHTING);
  xleft = xtimeleft;
  if(fontindex==LARGE_FONT)xleft=xtimeleft+0.11;
  xright = xtimeright;

  glLineWidth(linewidth);
  glBegin(GL_LINE_LOOP);
  glColor4fv(timebarcolor);
  glVertex2f(xleft,ybot);
  glVertex2f(xright,ybot);
  glVertex2f(xright,ytop);
  glVertex2f(xleft,ytop);
  glEnd();

  if(nglobal_times != 1){
    xxright = xleft + (float)itimes*(xright-xleft)/(nglobal_times-1);
  }
  else{
    xxright=xright;
  }
  glBegin(GL_POLYGON);
  glColor4fv(timebarcolor);
  glVertex2f(xleft,ybot);
  glVertex2f(xxright,ybot);
  glVertex2f(xxright,ytop);
  glVertex2f(xleft,ytop);
  glEnd();

}

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
    float *rgbi;

    rgbi=cbi->colorbar+3*i;
    glColor3fv(rgbi);
    glVertex3fv(rgbi);
  }
  glEnd();

  glPointSize(10.0);
  glBegin(GL_POINTS);
  for(i=0;i<cbi->nnodes;i++){
    rrgb=cbi->rgb_node+3*i;
    glColor3ubv(rrgb);
    glVertex3f(rrgb[0]/255.0,rrgb[1]/255.0,rrgb[2]/255.0);
  }
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
      float *rgbi;

      rgbi = cbi->colorbar+3*cbi->index_node[i];
      dzpoint = (float)cbi->index_node[i]/255.0;
      glColor3fv(rgbi);
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
      float *rgbi;

      glPointSize(20.0);
      glBegin(GL_POINTS);
      rgbi = cbi->colorbar+3*cbi->index_node[colorbarpoint];
      dzpoint = (float)cbi->index_node[colorbarpoint]/255.0;
      glColor3fv(rgbi);
      glVertex3f(1.5,0.0,dzpoint);
      glEnd();
    }

    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      float *rgbi;

      rgbi=cbi->colorbar+3*i;
      glColor3fv(rgbi);
      zbot=(float)i/255.0;
      glVertex3f(1.1,0.0,zbot);
      glVertex3f(1.3,0.0,zbot);
    }
    glEnd();

    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      float *rgbi;

      rgbi=cbi->colorbar+3*i;
      glColor3fv(rgbi);
      zbot=(float)i/255.0;
      glVertex3f(1.3,0.0,zbot);
      glVertex3f(1.1,0.0,zbot);
    }
    glEnd();
    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      float *rgbi;

      rgbi=cbi->colorbar+3*i;
      glColor3fv(rgbi);
      zbot=(float)i/255.0;
      glVertex3f(1.2,-0.1,zbot);
      glVertex3f(1.2, 0.1,zbot);
    }
    glEnd();
    glBegin(GL_QUAD_STRIP);
    for(i=0;i<256;i++){
      float *rgbi;

      rgbi=cbi->colorbar+3*i;
      glColor3fv(rgbi);
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
    colorbar[3*i]=rgb_node[0]/255.0;
    colorbar[1+3*i]=rgb_node[1]/255.0;
    colorbar[2+3*i]=rgb_node[2]/255.0;
    if(
      (rgb_node[0]==0&&rgb_node[1]==1&&rgb_node[2]==2)||
      (rgb_node[0]==253&&rgb_node[1]==254&&rgb_node[2]==255)
      ){
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
      colorbar[3*j]=MIX(factor,rgb_node[3],rgb_node[0])/255.0;
      colorbar[1+3*j]=MIX(factor,rgb_node[4],rgb_node[1])/255.0;
      colorbar[2+3*j]=MIX(factor,rgb_node[5],rgb_node[2])/255.0;
      if(
        (rgb_node[0]==0&&rgb_node[1]==1&&rgb_node[2]==2&&
        rgb_node[3]==0&&rgb_node[4]==1&&rgb_node[5]==2)||
        (rgb_node[0]==253&&rgb_node[1]==254&&rgb_node[2]==255&&
         rgb_node[3]==253&&rgb_node[4]==254&&rgb_node[5]==255)
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
    colorbar[3*i]=rgb_node[0]/255.0;
    colorbar[1+3*i]=rgb_node[1]/255.0;
    colorbar[2+3*i]=rgb_node[2]/255.0;
    if(
      (rgb_node[0]==0&&rgb_node[1]==1&&rgb_node[2]==2)||
      (rgb_node[0]==253&&rgb_node[1]==254&&rgb_node[2]==255)
      )
    {
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
    colorbar[3*255]=rgb_above_max[0];
    colorbar[1+3*255]=rgb_above_max[1];
    colorbar[2+3*255]=rgb_above_max[2];
  }
  CheckMemory;
}

/* ------------------ initdefaultcolorbars ------------------------ */

void initdefaultcolorbars(void){
  int i;
  colorbardata *cbi;

  ndefaultcolorbars=13;
  
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


  // blue->yellow->orange split

  strcpy(cbi->label,"blue->yellow->orange split");
  cbi->label_ptr=cbi->label;

  cbi->nnodes=6;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=0;
  cbi->rgb_node[1]=0;
  cbi->rgb_node[2]=255;

  cbi->index_node[1]=26;
  cbi->rgb_node[3]=0;
  cbi->rgb_node[4]=0;
  cbi->rgb_node[5]=255;

  cbi->index_node[2]=26;
  cbi->rgb_node[6]=255;
  cbi->rgb_node[7]=255;
  cbi->rgb_node[8]=0;

  cbi->index_node[3]=85;
  cbi->rgb_node[9]=255;
  cbi->rgb_node[10]=255;
  cbi->rgb_node[11]=0;

  cbi->index_node[4]=85;
  cbi->rgb_node[12]=255;
  cbi->rgb_node[13]=155;
  cbi->rgb_node[14]=0;

  cbi->index_node[5]=255;
  cbi->rgb_node[15]=255;
  cbi->rgb_node[16]=155;
  cbi->rgb_node[17]=0;
  cbi++;

  // fire (original)

  fire_colorbar_index=cbi-colorbarinfo;
  fire_colorbar=cbi;
  strcpy(cbi->label,"fire (original)");
  cbi->label_ptr=cbi->label;

  cbi->nnodes=4;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=0;
  cbi->rgb_node[1]=0;
  cbi->rgb_node[2]=0;

  cbi->index_node[1]=127;
  cbi->rgb_node[3]=0;
  cbi->rgb_node[4]=0;
  cbi->rgb_node[5]=0;

  cbi->index_node[2]=128;
  cbi->rgb_node[6]=255;
  cbi->rgb_node[7]=128;
  cbi->rgb_node[8]=0;

  cbi->index_node[3]=255;
  cbi->rgb_node[9]=255;
  cbi->rgb_node[10]=128;
  cbi->rgb_node[11]=0;
  cbi++;

  // fire (black->orange)

  fire_colorbar_index=cbi-colorbarinfo;
  fire_colorbar=cbi;
  strcpy(cbi->label,"fire (black->orange)");
  cbi->label_ptr=cbi->label;

  cbi->nnodes=11;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=0;
  cbi->rgb_node[1]=0;
  cbi->rgb_node[2]=0;

  cbi->index_node[1]=128;
  cbi->rgb_node[3]=0;
  cbi->rgb_node[4]=0;
  cbi->rgb_node[5]=0;

  cbi->index_node[2]=129;
  cbi->rgb_node[6]=255;
  cbi->rgb_node[7]=51;
  cbi->rgb_node[8]=0;

  cbi->index_node[3]=143;
  cbi->rgb_node[9]=255;
  cbi->rgb_node[10]=69;
  cbi->rgb_node[11]=0;

  cbi->index_node[4]=159;
  cbi->rgb_node[12]=255;
  cbi->rgb_node[13]=82;
  cbi->rgb_node[14]=0;

  cbi->index_node[5]=175;
  cbi->rgb_node[15]=255;
  cbi->rgb_node[16]=93;
  cbi->rgb_node[17]=0;

  cbi->index_node[6]=191;
  cbi->rgb_node[18]=255;
  cbi->rgb_node[19]=102;
  cbi->rgb_node[20]=0;

  cbi->index_node[7]=207;
  cbi->rgb_node[21]=255;
  cbi->rgb_node[22]=111;
  cbi->rgb_node[23]=0;

  cbi->index_node[8]=223;
  cbi->rgb_node[24]=255;
  cbi->rgb_node[25]=118;
  cbi->rgb_node[26]=0;

  cbi->index_node[9]=239;
  cbi->rgb_node[27]=255;
  cbi->rgb_node[28]=124;
  cbi->rgb_node[29]=0;

  cbi->index_node[10]=255;
  cbi->rgb_node[30]=255;
  cbi->rgb_node[31]=130;
  cbi->rgb_node[32]=0;
  cbi++;

  // fire (new)

  fire_colorbar_index=cbi-colorbarinfo;
  fire_colorbar=cbi;
  strcpy(cbi->label,"fire (new)");
  cbi->label_ptr=cbi->label;

  cbi->nnodes=7;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=0;
  cbi->rgb_node[1]=0;
  cbi->rgb_node[2]=0;

  cbi->index_node[1]=128;
  cbi->rgb_node[3]=0;
  cbi->rgb_node[4]=0;
  cbi->rgb_node[5]=0;

  cbi->index_node[2]=129;
  cbi->rgb_node[6]=255;
  cbi->rgb_node[7]=0;
  cbi->rgb_node[8]=0;

  cbi->index_node[3]=205;
  cbi->rgb_node[9]=255;
  cbi->rgb_node[10]=69;
  cbi->rgb_node[11]=0;

  cbi->index_node[4]=222;
  cbi->rgb_node[12]=255;
  cbi->rgb_node[13]=126;
  cbi->rgb_node[14]=0;

  cbi->index_node[5]=239;
  cbi->rgb_node[15]=255;
  cbi->rgb_node[16]=141;
  cbi->rgb_node[17]=0;

  cbi->index_node[6]=255;
  cbi->rgb_node[15]=255;
  cbi->rgb_node[16]=165;
  cbi->rgb_node[17]=0;
  cbi++;
  
  // fire (custom)

  fire_custom_colorbar=cbi;
  strcpy(cbi->label,"fire (custom)");
  cbi->label_ptr=cbi->label;

  cbi->nnodes=4;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=smoke_shade*255;
  cbi->rgb_node[1]=smoke_shade*255;
  cbi->rgb_node[2]=smoke_shade*255;

  cbi->index_node[1]=126;
  cbi->rgb_node[3]=smoke_shade*255;
  cbi->rgb_node[4]=smoke_shade*255;
  cbi->rgb_node[5]=smoke_shade*255;

  cbi->index_node[2]=131;
  cbi->rgb_node[6]=fire_red;
  cbi->rgb_node[7]=fire_green;
  cbi->rgb_node[8]=fire_blue;

  cbi->index_node[3]=255;
  cbi->rgb_node[9]=fire_red;
  cbi->rgb_node[10]=fire_green;
  cbi->rgb_node[11]=fire_blue;

  cbi++;

  // fire line (level set)

  levelset_colorbar=cbi-colorbarinfo;
  strcpy(cbi->label,"fire line (level set)");
  cbi->label_ptr=cbi->label;

  cbi->nnodes=6;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=64;
  cbi->rgb_node[1]=64;
  cbi->rgb_node[2]=64;

  cbi->index_node[1]=120;
  cbi->rgb_node[3]=64;
  cbi->rgb_node[4]=64;
  cbi->rgb_node[5]=64;

  cbi->index_node[2]=120;
  cbi->rgb_node[6]=255;
  cbi->rgb_node[7]=0;
  cbi->rgb_node[8]=0;

  cbi->index_node[3]=136;
  cbi->rgb_node[9]=255;
  cbi->rgb_node[10]=0;
  cbi->rgb_node[11]=0;

  cbi->index_node[4]=136;
  cbi->rgb_node[12]=0;
  cbi->rgb_node[13]=1;
  cbi->rgb_node[14]=2;

  cbi->index_node[5]=255;
  cbi->rgb_node[15]=0;
  cbi->rgb_node[16]=1;
  cbi->rgb_node[17]=2;
  cbi++;


  // fire line (wall thickness)

  wallthickness_colorbar=cbi-colorbarinfo;
  strcpy(cbi->label,"fire line (wall thickness)");
  cbi->label_ptr=cbi->label;

  cbi->nnodes=4;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=0;
  cbi->rgb_node[1]=0;
  cbi->rgb_node[2]=0;

  cbi->index_node[1]=32;
  cbi->rgb_node[3]=0;
  cbi->rgb_node[4]=0;
  cbi->rgb_node[5]=0;

  cbi->index_node[2]=32;
  cbi->rgb_node[6]=253;
  cbi->rgb_node[7]=254;
  cbi->rgb_node[8]=255;

  cbi->index_node[3]=255;
  cbi->rgb_node[9]=253;
  cbi->rgb_node[10]=254;
  cbi->rgb_node[11]=255;

  cbi++;

  // black->white

  bw_colorbar_index = cbi - colorbarinfo;
  strcpy(cbi->label,"black->white");
  cbi->label_ptr=cbi->label;

  cbi->nnodes=2;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=0;
  cbi->rgb_node[1]=0;
  cbi->rgb_node[2]=0;

  cbi->index_node[1]=255;
  cbi->rgb_node[3] =255;
  cbi->rgb_node[4]=255;
  cbi->rgb_node[5]=255;
  cbi++;

// construct colormaps from color node info

  rgb_below_min[0]=255-64;
  rgb_below_min[1]=255-64;
  rgb_below_min[2]=255-64;
  rgb_above_max[0]=0;
  rgb_above_max[1]=0;
  rgb_above_max[2]=0;

  for(i=0;i<ndefaultcolorbars;i++){
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

/* ------------------ get_label_position ------------------------ */

int get_label_position(float position, float dyfont, float barbot){
  float diff_position, min_diff;
  int iposition;
  int i;

  iposition=-1;
  min_diff = 1000.0;
  for (i=0; i<nrgb-1; i++){
    float vert_position;

    vert_position = (float)(i)*(float)(nrgb+DYFONT)/(float)(nrgb-2) + barbot-dyfont/2.0;
    diff_position = position - vert_position;
    if(diff_position<0.0)diff_position=-diff_position;
    if(diff_position<min_diff){
      iposition=i;
      min_diff=diff_position;
    }
  }
  if(min_diff>0.2)iposition=-1;
  return iposition;
}

/* ------------------ drawColorBars ------------------------ */

void drawColorBars(void){
  float dyfont;
  int dyscreen;
  int temp;

  int i,i3;
  float labeltop;
  float barleft;
  float dbar=.20f, dtext=.15f, dspace=.07f, space;
  float barbot=-0.5f;
  float right[3],left[3],bottom[4];
  int ileft=0;
  int leftzone, leftsmoke, leftslice, leftpatch, leftiso;
  int lastzone, lastsmoke, lastslice, lastpatch, lastiso;
  float yy,yy2;
  int iposition;
  float position;
  float tttval, tttmin, tttmax;
  float val;
  char slicelabel[256], p3dlabel[256], boundarylabel[256], partlabel[256], zonelabel[256],isolabel[256];
  char unitlabel[256];
  float *p3lev;
  databounds *sb;
  patchdata *patchi;

  int sliceunitclass,sliceunittype;
  int sliceflag=0;
  int isoflag=0;
  float *slicefactor=NULL;
  float slicefactor2[2];
  float slicerange,isorange;
  float *isofactor=NULL;

  int plot3dunitclass, plot3dunittype;
  int plot3dflag=0;
  float *plot3dfactor=NULL;
  float plot3dfactor2[2];
  float plot3drange;

  int patchunitclass, patchunittype;
  int zoneunitclass, zoneunittype;
  int patchflag=0;
  int zoneflag=0;
  float *patchfactor=NULL;
  float *zonefactor=NULL;
  float patchrange=0.0;
  float zonerange;

  int partunitclass, partunittype;
  int partflag=0;
  float *partfactor=NULL;
  float partrange=0.0;

  char partunitlabel2[256], partshortlabel2[256];

  GLfloat *color1, *color2;

  color1=&(foregroundcolor[0]);
  color2=&(redcolor[0]);

  barleft = barright-dbar;
  right[0] = barleft;
  space = dtext + dspace;
  if(fontindex==LARGE_FONT)space *= 1.5;
  left[0] = barleft - space;
  if(left[0]<0.0f)left[0]=0.0;

  right[1]=left[0];
  left[1] = left[0] - space;
  if(left[1]<0.0f)left[1]=0.0;

  right[2]=left[1];
  left[2] = left[1]-space;
  if(left[2]<0.0f)left[2]=0.0;

  labeltop=nrgb+1.0;

  temp = (int)(1.2f*dwinH);
  dyscreen=screenHeight-temp-fontHoffset-2*titlesafe_offset;

  switch (fontindex){

  case SMALL_FONT:
    dyfont = (float)(small_font_height+1)*(float)(nrgb+1)/(float)dyscreen;
    break;
  case LARGE_FONT:
    dyfont = (float)(large_font_height+1)*(float)(nrgb+1)/(float)dyscreen;
    break;
  default:
    ASSERT(FFALSE);
    break;
  }

  bottom[0]=labeltop-dyfont;
  bottom[1]=labeltop-2*dyfont;
  bottom[2]=labeltop-3*dyfont;
  bottom[3]=labeltop-4*dyfont;
  if(showiso_colorbar==1||showevac_colorbar==1||
    (showsmoke==1&&parttype!=0)||showslice==1||
    (showvslice==1&&vslicecolorbarflag==1)||
    showpatch==1||
    (showzone==1&&sethazardcolor==0)||
    showplot3d==1){
    SNIFF_ERRORS("before colorbar");
    CheckMemory;
    if(showplot3d==1&&contour_type==STEPPED_CONTOURS){
      float top, tophat;

      // draw plot3d colorbars

      top = barbot+nrgb+DYFONT;
      tophat = top - (top-barbot)/(nrgb-2);
      glBegin(GL_QUADS);
      for (i = 0; i < nrgb-2; i++){
        float *rgb_plot3d_local;

        rgb_plot3d_local = rgb_plot3d_contour[i];

        yy = MIX2(i,nrgb-3,tophat,barbot);
        yy2 = MIX2(i+1,nrgb-3,tophat,barbot);

        if(rgb_plot3d_local[3]!=0.0){
          glColor4fv(rgb_plot3d_local);
          glVertex2f(barleft, yy); 
          glVertex2f(barright,yy);
       
          glVertex2f(barright,yy2);
          glVertex2f(barleft, yy2);
        }
      }
      glEnd();
      if(show_extremedata==1){
        float barmid;
        float *rgb_plot3d_local;

        rgb_plot3d_local = rgb_plot3d_contour[nrgb-2];
        barmid = (barleft+barright)/2.0;
        i=-1;
        yy = MIX2(i+0.5,nrgb-3,tophat,barbot);
        yy2 = MIX2(i+1,nrgb-3,tophat,barbot);


        if(show_extreme_below==1||show_extreme_above==1)glEnable(GL_POLYGON_SMOOTH);

        if(show_extreme_below==1&&rgb_plot3d_local[3]!=0.0){     
          glBegin(GL_TRIANGLES);
          glColor4fv(rgb_plot3d_local);

          glVertex2f(barleft,yy2);
          glVertex2f(barmid,yy);
          glVertex2f(barright,yy2);
          glEnd();
        }

        i=nrgb-2;
        yy = MIX2(i,nrgb-3,tophat,barbot);
        yy2 = MIX2(i+0.5,nrgb-3,tophat,barbot);

        rgb_plot3d_local = rgb_plot3d_contour[nrgb-1];
        if(show_extreme_above==1&&rgb_plot3d_local[3]!=0.0){
          glBegin(GL_TRIANGLES);
          glColor4fv(rgb_plot3d_local);
          glVertex2f(barleft, yy); 
          glVertex2f(barright,yy);
       
          glVertex2f(barmid, yy2);
          glEnd();
        }
        if(show_extreme_below==1||show_extreme_above==1)glDisable(GL_POLYGON_SMOOTH);
      }
    }
    else{
      float top, tophat;

      // draw all other colorbars
      top = barbot+nrgb+DYFONT;
      tophat = top - (top-barbot)/255.0;
      glBegin(GL_QUADS);
      for (i = 0; i < nrgb_full-1; i++){
        float *rgb_cb,*rgb_cb2;

        rgb_cb=rgb_full[i];

        yy = MIX2(i,255,tophat,barbot);
        yy2 = MIX2(i+1,255,tophat,barbot);
        i3=i+1;
        if(i==nrgb_full-2)i3=i;
        rgb_cb2=rgb_full[i3];

        if(rgb_cb[3]!=0.0&&rgb_cb2[3]!=0.0){
          glColor4fv(rgb_cb); 
          glVertex2f(barleft, yy);
          glVertex2f(barright,yy);

          glColor4fv(rgb_cb2);
          glVertex2f(barright,yy2);
          glVertex2f(barleft,yy2);
        }
      }
      glEnd();
    }
    if(show_extremedata==1){
      float top, tophat, barmid;

      barmid=(barleft+barright)/2.0;
      top = barbot+nrgb+DYFONT;
      tophat = top - (top-barbot)/(nrgb-2);
      i=-1;
      yy = MIX2(i+0.5,nrgb-3,tophat,barbot);
      yy2 = MIX2(i+1,nrgb-3,tophat,barbot);

      if(show_extreme_below==1||show_extreme_above==1)glEnable(GL_POLYGON_SMOOTH);

      if(show_extreme_below==1){
        glBegin(GL_TRIANGLES);
        glColor4fv(rgb_full[0]);

        glVertex2f(barleft, yy2);
        glVertex2f(barmid, yy); 
        glVertex2f(barright,yy2);
        glEnd();
      }

      i=nrgb-2;
      yy = MIX2(i,nrgb-3,tophat,barbot)    +(barbot-tophat)/255.0;
      yy2 = MIX2(i+0.5,nrgb-3,tophat,barbot)+(barbot-tophat)/255.0;

      if(show_extreme_above==1){
        glBegin(GL_TRIANGLES);
        glColor4fv(rgb_full[nrgb_full-1]);
        glVertex2f(barleft, yy); 
        glVertex2f(barright,yy);
        glVertex2f(barmid, yy2);
        glEnd();
      }
      if(show_extreme_below==1||show_extreme_above==1)glDisable(GL_POLYGON_SMOOTH);
    }
  }
  lastsmoke=0;
  lastslice=0;
  lastpatch=0;
  lastzone=0;
  lastiso=0;
  leftsmoke=0;
  leftslice=0;
  leftpatch=0;
  leftzone=0;
  leftiso=0;
  ileft=0;
  if(showiso_colorbar==1){
    leftiso=ileft;
    lastsmoke=0;lastslice=0;lastpatch=0;lastzone=0;lastiso=1;
    ileft++;
  }
  if(showevac_colorbar==1||showsmoke==1){
    if(parttype!=0){
      leftsmoke=ileft;
      lastsmoke=1;lastslice=0;lastpatch=0;lastzone=0;lastiso=0;
      ileft++;
    }
  }
  if(showslice==1||(showvslice==1&&vslicecolorbarflag==1)){
    leftslice=ileft;
    lastsmoke=0;lastslice=1;lastpatch=0;lastzone=0;lastiso=0;
    ileft++;
  }
  if(showpatch==1){
    leftpatch=ileft;
    lastsmoke=0;lastslice=0;lastpatch=1;lastzone=0;lastiso=0;
  }

  strcpy(partshortlabel2,"");
  strcpy(partunitlabel2,"");
  if(showevac_colorbar==1||showsmoke==1){
    if(parttype!=0){
      if(showsmoke==1&&showevac==0)outputBarText(right[leftsmoke],bottom[0],color1,"Part");
      if(showevac==1)outputBarText(right[leftsmoke],bottom[0],color1,"Human");
    }
    if(parttype==-1){
      strcpy(partshortlabel2,"temp");
      strcpy(partunitlabel2,"C");
    }
    else if(parttype==-2){
      strcpy(partshortlabel2,"HRRPUV");
      strcpy(partunitlabel2,"kW/m^3");
    }
    else{
      if(partshortlabel!=NULL)strcpy(partshortlabel2,partshortlabel);
      if(partunitlabel!=NULL)strcpy(partunitlabel2,partunitlabel);
    }
    if(parttype!=0){
      getunitinfo(partunitlabel2,&partunitclass,&partunittype);
      if(partunitclass>=0&&partunitclass<nunitclasses){
        if(partunittype>=0){
          partflag=1;
          partfactor=unitclasses[partunitclass].units[partunittype].scale;
          strcpy(partunitlabel,unitclasses[partunitclass].units[partunittype].unit);
        }
      }
      outputBarText(right[leftsmoke],bottom[1],color1,partshortlabel);
      outputBarText(right[leftsmoke],bottom[2],color1,partunitlabel);
      outputBarText(right[leftsmoke],bottom[3],color1,partscale);
    }
  }
  if(showslice==1||(showvslice==1&&vslicecolorbarflag==1)){
    sb = slicebounds + islicetype;
    strcpy(unitlabel,sb->label->unit);
    getunitinfo(sb->label->unit,&sliceunitclass,&sliceunittype);
    if(sliceunitclass>=0&&sliceunitclass<nunitclasses){
      if(sliceunittype>0){
        sliceflag=1;
        slicefactor=unitclasses[sliceunitclass].units[sliceunittype].scale;
        strcpy(unitlabel,unitclasses[sliceunitclass].units[sliceunittype].unit);
      }
    }
    outputBarText(right[leftslice],bottom[0],color1,"Slice");
    outputBarText(right[leftslice],bottom[1],color1,sb->label->shortlabel);
    outputBarText(right[leftslice],bottom[2],color1,unitlabel);
    if(strcmp(unitlabel,"ppm")==0&&slicefactor!=NULL){
      slicefactor2[0]=*slicefactor*sb->fscale;
      slicefactor2[1]=0.0;
      slicefactor=slicefactor2;
    }
    else{
      outputBarText(right[leftslice],bottom[3],color1,sb->scale);
    }
  }
  if(showiso_colorbar==1){
    sb = isobounds + iisottype;
    strcpy(unitlabel,sb->label->unit);
    /*
    getunitinfo(sb->label->unit,&isounitclass,&isounittype);
    if(isounitclass>=0&&isounitclass<nunitclasses){
      if(isounittype>0){
        isoflag=1;
        isofactor=unitclasses[isounitclass].units[isounittype].scale;
        strcpy(unitlabel,unitclasses[isounitclass].units[isounittype].unit);
      }
    }
    */
    outputBarText(right[leftiso],bottom[0],color1,"Iso");
    outputBarText(right[leftiso],bottom[1],color1,sb->label->shortlabel);
    outputBarText(right[leftiso],bottom[2],color1,unitlabel);
    outputBarText(right[leftiso],bottom[3],color1,sb->scale);
  }
  if(showpatch==1){
    patchi = patchinfo + patchtypes[ipatchtype];
    strcpy(unitlabel,patchi->label.unit);
    getunitinfo(patchi->label.unit,&patchunitclass,&patchunittype);
    if(patchunitclass>=0&&patchunitclass<nunitclasses){
      if(patchunittype>0){
        patchflag=1;
        patchfactor=unitclasses[patchunitclass].units[patchunittype].scale;
        strcpy(unitlabel,unitclasses[patchunitclass].units[patchunittype].unit);
      }
    }
    outputBarText(right[leftpatch],bottom[0],color1,"Bndry");
    outputBarText(right[leftpatch],bottom[1],color1,patchi->label.shortlabel);
    outputBarText(right[leftpatch],bottom[2],color1,unitlabel);
    outputBarText(right[leftpatch],bottom[3],color1,patchi->scale);
  }
  if(showplot3d==1){
    char *p3label;
    char *up3label;

    up3label = plot3dinfo[0].label[plotn-1].unit;
    strcpy(unitlabel,up3label);
    getunitinfo(up3label,&plot3dunitclass,&plot3dunittype);
    if(plot3dunitclass>=0&&plot3dunitclass<nunitclasses){
      if(plot3dunittype>0){
        plot3dflag=1;
        plot3dfactor=unitclasses[plot3dunitclass].units[plot3dunittype].scale;
        strcpy(unitlabel,unitclasses[plot3dunitclass].units[plot3dunittype].unit);
      }
    }
    p3label = plot3dinfo[0].label[plotn-1].shortlabel;
    outputBarText(right[0],bottom[0],color1,"Plot3d");
    outputBarText(right[0],bottom[1],color1,p3label);
    outputBarText(right[0],bottom[2],color1,unitlabel);
    if(strcmp(unitlabel,"ppm")==0&&plot3dfactor!=NULL){
      plot3dfactor2[0]=*plot3dfactor*fscalep3[plotn-1];
      plot3dfactor2[1]=0.0;
      plot3dfactor=plot3dfactor2;
    }
    else{
      outputBarText(right[0],bottom[3],color1,scalep3[plotn-1]);
    }

  }
  if(showzone==1&&sethazardcolor==0){
    strcpy(unitlabel,"C");
    getunitinfo(unitlabel,&zoneunitclass,&zoneunittype);
    if(zoneunitclass>=0&&zoneunitclass<nunitclasses){
      if(zoneunittype>0){
        zoneflag=1;
        zonefactor=unitclasses[zoneunitclass].units[zoneunittype].scale;
        strcpy(unitlabel,unitclasses[zoneunitclass].units[zoneunittype].unit);
      }
    }
    outputBarText(right[leftzone],bottom[0],color1,"Zone");
    outputBarText(right[leftzone],bottom[1],color1,"Temp");
    outputBarText(right[leftzone],bottom[2],color1,unitlabel);
    outputBarText(right[leftzone],bottom[3],color1,zonescale);

  }
  if(showiso_colorbar==1){
    sb = isobounds + iisottype;
    tttmin = sb->levels256[0];
    tttmax = sb->levels256[255];
    isorange=tttmax-tttmin;
    iposition=-1;
    if(global_changecolorindex!=-1){
      char isocolorlabel[256];
      char *isocolorlabel_ptr=NULL;

      tttval = sb->levels256[valindex];
      num2string(isolabel,tttval,isorange);
      isocolorlabel_ptr=isolabel;
      if(isoflag==1){
        scalefloat2string(tttval,isocolorlabel, isofactor, isorange);
        isocolorlabel_ptr=isocolorlabel;
      }
      position = (float)global_changecolorindex/255.0*(float)(nrgb+DYFONT)+barbot-dyfont/2.0;
      iposition = get_label_position(position,dyfont,barbot);
      outputBarText(right[leftiso],position,color2,isocolorlabel_ptr);
    }
    for (i=0; i<nrgb-1; i++){
      float vert_position;
      char isocolorlabel[256];
      char *isocolorlabel_ptr=NULL;

      vert_position = (float)(i)*(float)(nrgb+DYFONT)/(float)(nrgb-2) + barbot-dyfont/2.0;
      if(iposition==i)continue;
      isocolorlabel_ptr=&(sb->colorlabels[i+1][0]);
      if(isoflag==1){
        val = tttmin + i*isorange/(nrgb-2);
        scalefloat2string(val,isocolorlabel, isofactor, isorange);
        isocolorlabel_ptr=isocolorlabel;
      }
      outputBarText(right[leftiso],vert_position,color1,isocolorlabel_ptr);
    }
  }
  if(showevac_colorbar==1||(showsmoke==1&&parttype!=0)){
    float *partlevels256_ptr;

    partlevels256_ptr=partlevels256;
    if(prop_index>=0&&prop_index<npart5prop){
      partlevels256_ptr=part5propinfo[prop_index].ppartlevels256;
    }

    iposition=-1;
    tttmin = partlevels256_ptr[0];
    tttmax = partlevels256_ptr[255];
    partrange = tttmax - tttmin;
    if(global_changecolorindex!=-1){
     char partcolorlabel[256];
     char *partcolorlabel_ptr=NULL;

      tttval = partlevels256_ptr[valindex];
      num2string(partlabel,tttval,partrange);
      partcolorlabel_ptr=partlabel;
      if(partflag==1){
        scalefloat2string(tttval,partcolorlabel, partfactor, partrange);
        partcolorlabel_ptr=partcolorlabel;
      }
      position = (float)global_changecolorindex/255.0*(float)(nrgb+DYFONT)+barbot-dyfont/2.0;
      iposition = get_label_position(position,dyfont,barbot);
      outputBarText(right[leftsmoke],position,color2,partcolorlabel_ptr);
    }
    for (i=0; i<nrgb-1; i++){
      float vert_position;
      char partcolorlabel[256],partcolorlabel2[256];
      char *partcolorlabel_ptr=NULL;

      vert_position = (float)(i)*(float)(nrgb+DYFONT)/(float)(nrgb-2) + barbot-dyfont/2.0;
      if(iposition==i)continue;
      if(prop_index>=0&&prop_index<npart5prop){
        partcolorlabel_ptr=&part5propinfo[prop_index].partlabels[i+1][0];
      }
      else{
        if(colorlabelpart!=NULL){
          partcolorlabel_ptr=&colorlabelpart[i+1][0];
        }
        else{
          partcolorlabel_ptr=NULL;
        }
      }
      if(partflag==1){
        val = tttmin + i*partrange/(nrgb-2);
        scalefloat2string(val,partcolorlabel, partfactor, partrange);
        scalestring(partcolorlabel_ptr,partcolorlabel, partfactor, partrange);
        partcolorlabel_ptr=partcolorlabel;
      }
      if(lastsmoke==1&&showcolorbar_gtge==1){
        if(show_extreme_below==1&&i==0){
          strcpy(partcolorlabel2,"<=");
          strcat(partcolorlabel2,partcolorlabel_ptr);
          partcolorlabel_ptr=partcolorlabel2;
        }
        if(show_extreme_above==1&&i==nrgb-2){
          strcpy(partcolorlabel2,">=");
          strcat(partcolorlabel2,partcolorlabel_ptr);
          partcolorlabel_ptr=partcolorlabel2;
        }
      }
      outputBarText(right[leftsmoke],vert_position,color1,partcolorlabel_ptr);
    }
  }
  if(showslice==1||(showvslice==1&&vslicecolorbarflag==1)){
    sb=slicebounds+islicetype;
    tttmin = sb->levels256[0];
    tttmax = sb->levels256[255];
    slicerange=tttmax-tttmin;
    iposition=-1;
    if(global_changecolorindex!=-1){
      char slicecolorlabel[256];
      char *slicecolorlabel_ptr=NULL;

      tttval = sb->levels256[valindex];
      num2string(slicelabel,tttval,slicerange);
      slicecolorlabel_ptr=slicelabel;
      if(sliceflag==1){
        scalefloat2string(tttval,slicecolorlabel, slicefactor, slicerange);
        slicecolorlabel_ptr=slicecolorlabel;
      }
      position = (float)global_changecolorindex/255.0*(float)(nrgb+DYFONT)+barbot-dyfont/2.0;
      iposition = get_label_position(position,dyfont,barbot);
      outputBarText(right[leftslice],position,color2,slicecolorlabel_ptr);
    }
    for (i=0; i<nrgb-1; i++){
      float vert_position;
      char slicecolorlabel[256],slicecolorlabel2[256];
      char *slicecolorlabel_ptr=NULL;

      vert_position = (float)(i)*(float)(nrgb+DYFONT)/(float)(nrgb-2) + barbot-dyfont/2.0;
      if(iposition==i)continue;
      slicecolorlabel_ptr=&(sb->colorlabels[i+1][0]);
      if(sliceflag==1){
        val = tttmin + i*slicerange/(nrgb-2);
        scalefloat2string(val,slicecolorlabel, slicefactor, slicerange);
        slicecolorlabel_ptr=slicecolorlabel;
      }
      if(lastslice==1&&showcolorbar_gtge==1){
        if(show_extreme_below==1&&i==0){
          strcpy(slicecolorlabel2,"<=");
          strcat(slicecolorlabel2,slicecolorlabel_ptr);
          slicecolorlabel_ptr=slicecolorlabel2;
        }
        if(show_extreme_above==1&&i==nrgb-2){
          strcpy(slicecolorlabel2,">=");
          strcat(slicecolorlabel2,slicecolorlabel_ptr);
          slicecolorlabel_ptr=slicecolorlabel2;
        }
      }
      outputBarText(right[leftslice],vert_position,color1,slicecolorlabel_ptr);
    }
  }
  if(showpatch==1){
    iposition=-1;
    tttmin = boundarylevels256[0];
    tttmax = boundarylevels256[255];
    patchrange=tttmax-tttmin;
    if(global_changecolorindex!=-1){
      char patchcolorlabel[256];
      char *patchcolorlabel_ptr=NULL;

      // draw boundary file value selected with mouse
      tttval = boundarylevels256[valindex];
      num2string(boundarylabel,tttval,tttmax-tttmin);
      position = (float)global_changecolorindex/255.0*(float)(nrgb+DYFONT)+barbot-dyfont/2.0;
      iposition = get_label_position(position,dyfont,barbot);
      patchcolorlabel_ptr=&(boundarylabel[0]);
      if(patchflag==1){
        scalefloat2string(tttval,patchcolorlabel, patchfactor, patchrange);
        patchcolorlabel_ptr=patchcolorlabel;
      }
      outputBarText(right[leftpatch],position,color2,patchcolorlabel_ptr);
    }
    for (i=0; i<nrgb-1; i++){
      char patchcolorlabel[256],patchcolorlabel2[256];
      char *patchcolorlabel_ptr=NULL;
      float vert_position;

      vert_position = (float)(i)*(float)(nrgb+DYFONT)/(float)(nrgb-2) + barbot-dyfont/2.0;

      if(iposition==i)continue;
      patchcolorlabel_ptr=&colorlabelpatch[i+1][0];
      if(patchflag==1){
        val = tttmin + i*patchrange/(nrgb-2);
        scalefloat2string(val,patchcolorlabel, patchfactor, patchrange);
        patchcolorlabel_ptr=patchcolorlabel;
      }
      if(lastpatch==1&&showcolorbar_gtge==1){
        if(show_extreme_below==1&&i==0){
          strcpy(patchcolorlabel2,"<=");
          strcat(patchcolorlabel2,patchcolorlabel_ptr);
          patchcolorlabel_ptr=patchcolorlabel2;
        }
        if(show_extreme_above==1&&i==nrgb-2){
          strcpy(patchcolorlabel2,">=");
          strcat(patchcolorlabel2,patchcolorlabel_ptr);
          patchcolorlabel_ptr=patchcolorlabel2;
        }
      }
      outputBarText(right[leftpatch],vert_position,color1,patchcolorlabel_ptr);
    }
  }
  if(showzone==1&&sethazardcolor==0){
    iposition=-1;
    tttmin = zonelevels256[0];
    tttmax = zonelevels256[255];
    zonerange=tttmax-tttmin;
    if(global_changecolorindex!=-1){
      char zonecolorlabel[256];
      char *zonecolorlabel_ptr=NULL;

      tttval = zonelevels256[valindex];
      num2string(zonelabel,tttval,tttmax-tttmin);
      position = (float)global_changecolorindex/255.0*(float)(nrgb+DYFONT)+barbot-dyfont/2.0;
      iposition = get_label_position(position,dyfont,barbot);
      zonecolorlabel_ptr=&(zonelabel[0]);
      if(zoneflag==1){
        scalefloat2string(tttval,zonecolorlabel, zonefactor, zonerange);
        zonecolorlabel_ptr=zonecolorlabel;
      }
      outputBarText(right[leftzone],position,color2,zonecolorlabel_ptr);
    }
    for (i=0; i<nrgb-1; i++){
      float vert_position;
      char zonecolorlabel[256];
      char *zonecolorlabel_ptr=NULL;

      vert_position = (float)(i)*(float)(nrgb+DYFONT)/(float)(nrgb-2) + barbot-dyfont/2.0;
      if(iposition==i)continue;
      zonecolorlabel_ptr=&colorlabelzone[i+1][0];
      if(zoneflag==1){
        val = tttmin + (i-1)*zonerange/(nrgb-2);
        scalefloat2string(val,zonecolorlabel, zonefactor, zonerange);
        zonecolorlabel_ptr=zonecolorlabel;
      }
      outputBarText(right[leftzone],vert_position,color1,zonecolorlabel_ptr);
    }
  }
  if(showplot3d==1){
    iposition=-1;
    p3lev = p3levels256[plotn-1];
    tttmin = p3lev[0];
    tttmax = p3lev[255];
    plot3drange = tttmax - tttmin;
    if(global_changecolorindex!=-1){
      char plot3dcolorlabel[256];
      char *plot3dcolorlabel_ptr=NULL;

      tttval = p3lev[valindex];
      num2string(p3dlabel,tttval,tttmax-tttmin);
      plot3dcolorlabel_ptr = p3dlabel;
      if(plot3dflag==1){
        scalefloat2string(tttval,plot3dcolorlabel, plot3dfactor, plot3drange);
        plot3dcolorlabel_ptr=plot3dcolorlabel;
      }
      if(contour_type==SHADED_CONTOURS||contour_type==LINE_CONTOURS){
        position = (float)global_changecolorindex/255.0*(float)(nrgb+DYFONT)+barbot-dyfont/2.0;
      }
      else{
         yy = MIX2(i,nrgb,nrgb+DYFONT+barbot,barbot);
        yy2 = MIX2(i+1,nrgb,nrgb+DYFONT+barbot,barbot);
        position = (yy+yy2)/2.0;
      }
      position = (float)global_changecolorindex/255.0*(float)(nrgb+DYFONT)+barbot-dyfont/2.0;
      iposition = get_label_position(position,dyfont,barbot);
      outputBarText(right[0],position,color2,plot3dcolorlabel_ptr);
    }
    if(visiso==0){
      float ddbar,bar0,vert_position;

      ddbar = (float)(nrgb+DYFONT)/(float)(nrgb-2);
      bar0 = barbot-dyfont/2.0;

      for (i=0; i<nrgb-1; i++){
        char plot3dcolorlabel[256];
        char plot3dcolorlabel2[256];
        char *plot3dcolorlabel_ptr=NULL;

        vert_position = (float)(i)*ddbar+bar0;
        if(iposition==i)continue;
        plot3dcolorlabel_ptr=&colorlabelp3[plotn-1][i][0];
        if(plot3dflag==1){
          val = tttmin + i*plot3drange/(nrgb-2);
          scalefloat2string(val,plot3dcolorlabel, plot3dfactor, plot3drange);
          plot3dcolorlabel_ptr=plot3dcolorlabel;
        }
        if(showcolorbar_gtge==1){
          if((show_extremedata==0||show_extreme_below==0)&&i==0){
            strcpy(plot3dcolorlabel2,"<=");
            strcat(plot3dcolorlabel2,plot3dcolorlabel_ptr);
            plot3dcolorlabel_ptr=plot3dcolorlabel2;
          }  
          if((show_extremedata==0||show_extreme_above==0)&&i==nrgb-2){
            strcpy(plot3dcolorlabel2,">=");
            strcat(plot3dcolorlabel2,plot3dcolorlabel_ptr);
            plot3dcolorlabel_ptr=plot3dcolorlabel2;
          }
        }
        outputBarText(right[0],vert_position,color1,plot3dcolorlabel_ptr);
      }
    }
    else{
      float ddbar,bar0,vert_position;

      ddbar = (float)(nrgb+DYFONT)/(float)(nrgb-2);
      bar0 = barbot-dyfont/2.0;

      for (i=0; i<nrgb-2; i++){
        char plot3dcolorlabel[256];
        char *plot3dcolorlabel_ptr=NULL;

        vert_position = (float)(i+0.5)*ddbar + bar0;

        if(iposition==i)continue;
        plot3dcolorlabel_ptr=&colorlabeliso[plotn-1][i][0];
        if(plot3dflag==1){
          val = tttmin + (i-1)*plot3drange/(nrgb-2);
          scalefloat2string(val,plot3dcolorlabel, plot3dfactor, plot3drange);
          plot3dcolorlabel_ptr=plot3dcolorlabel;
        }
        if(isolevelindex==i||isolevelindex2==i){
          outputBarText(right[0],vert_position,color2,plot3dcolorlabel_ptr);
        }
        else{
          outputBarText(right[0],vert_position,color1,plot3dcolorlabel_ptr);
        }
      }
    }

  }
}



