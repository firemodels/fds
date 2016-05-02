#include "options.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include GLUT_H

#include "smokeviewvars.h"

/* ------------------ UpdateTimeLabels ------------------------ */

void UpdateTimeLabels(void){
  float time0;

  time0 = timeoffset;
  if(global_times!=NULL)time0 = timeoffset + global_times[itimes];
  if(vishmsTimelabel==1){
    int hour, min, sec,sec10;

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
      sprintf(hrrinfo->hrrlabel,"HRR: %4.1f W",hrr*1000.0);
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

void drawTimeBar(float xleft, float xright, float ybot, float ytop){
  float xxright;

  if(xright<=xleft)return;
  glDisable(GL_LIGHTING);

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

/* ------------------ newcolorbar ------------------------ */

colorbardata *newcolorbar(char *name, unsigned char *table, int ntable){
  colorbardata *newcolorbar;
  int i;
  unsigned char *rgb_node;

  ncolorbars++;
  CheckMemory;
  ResizeMemory((void **)&colorbarinfo, ncolorbars*sizeof(colorbardata));
  CheckMemory;

  // new colorbar

  newcolorbar = colorbarinfo+ncolorbars-1;

  strcpy(newcolorbar->label, name);
  newcolorbar->label_ptr = newcolorbar->label;
  newcolorbar->nnodes = ntable;
  newcolorbar->nodehilight = 0;
  rgb_node = newcolorbar->rgb_node;
  for(i = 0; i<ntable; i++){
    int ii;

    ii = i*255/(ntable-1);
    newcolorbar->index_node[i]=ii;
    *rgb_node++ = *table++;
    *rgb_node++ = *table++;
    *rgb_node++ = *table++;
  }

  remapcolorbar(newcolorbar);
  return newcolorbar;
}

/* ------------------ addcolorbar ------------------------ */

void addcolorbar(int icolorbar){
  colorbardata *cb_to, *cb_from;

  ncolorbars++;
  CheckMemory;
  ResizeMemory((void **)&colorbarinfo,ncolorbars*sizeof(colorbardata));
  UpdateCurrentColorbar(colorbarinfo + colorbartype);

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
  colorbardata *cbi;
  int ncolors;

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
    unsigned char *rrgb;

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
    unsigned char *rgbleft;

    rgbleft = cbi->rgb_node+3*colorbarpoint;

    glPointSize(20.0);
    glBegin(GL_POINTS);
    glColor3ubv(rgbleft);
    glVertex3f(rgbleft[0]/255.0,rgbleft[1]/255.0,rgbleft[2]/255.0);
    glEnd();
  }

  {
    float xdenorm, ydenorm, zdenorm;

    glPointSize(10.0);
    glBegin(GL_POINTS);
    for(i=0;i<cbi->nnodes;i++){
      float *rgbi;
      float dzpoint;

      rgbi = cbi->colorbar+3*cbi->index_node[i];
      dzpoint = (float)cbi->index_node[i]/255.0;
      glColor3fv(rgbi);
      glVertex3f(1.5,0.0,dzpoint);
    }
    glEnd();

    xdenorm = DENORMALIZE_X(1.55);
    ydenorm = DENORMALIZE_Y(0.0);
    if(fontindex==SCALED_FONT)scale_3dfont();
    glPushMatrix();
    glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
    glTranslatef(-xbar0,-ybar0,-zbar0);
    for(i=0;i<cbi->nnodes;i++){
      char cbuff[1024];
      float dzpoint;

      dzpoint = (float)cbi->index_node[i]/255.0;
      zdenorm = DENORMALIZE_Z(dzpoint);
      sprintf(cbuff,"%i",(int)cbi->index_node[i]);
      output3Text(foregroundcolor, xdenorm,ydenorm,zdenorm, cbuff);
    }
    glPopMatrix();
    glLineWidth(5.0);
    if(colorbarpoint>=0&&colorbarpoint<cbi->nnodes){
      float *rgbi;
      float dzpoint;

      glPointSize(20.0);
      glBegin(GL_POINTS);
      rgbi = cbi->colorbar+3*cbi->index_node[colorbarpoint];
      dzpoint = (float)cbi->index_node[colorbarpoint]/255.0;
      glColor3fv(rgbi);
      glVertex3f(1.5,0.0,dzpoint);
      glEnd();
    }
    if(show_firecolormap==1){
      char vvlabel[255];
      float vval_min, vval_cutoff, vval_max;

      if(smoke_render_option==RENDER_SLICE){
        vval_min=global_hrrpuv_min;
        vval_cutoff=global_hrrpuv_cutoff;
        vval_max=global_hrrpuv_max;
      }
      else{
        vval_min=temperature_min;
        vval_cutoff=temperature_cutoff;
        vval_max=temperature_max;
      }
      sprintf(vvlabel,"%4.0f",vval_min);
      output3Text(foregroundcolor, 1.0,0.0,0.0,vvlabel);

      sprintf(vvlabel,"%4.0f",vval_cutoff);
      output3Text(foregroundcolor, 1.0,0.0,(vval_cutoff-vval_min)/(vval_max-vval_min),vvlabel);

      sprintf(vvlabel,"%4.0f",vval_max);
      output3Text(foregroundcolor, 1.0,0.0,1.0,vvlabel);
    }

    if(show_firecolormap==1){
      ncolors=MAXSMOKERGB-1;
    }
    else{
      ncolors=MAXRGB-1;
    }
    glBegin(GL_TRIANGLES);
    for(i=1;i<ncolors;i++){
      float *rgbi;
      float zbot, ztop;

      if(show_firecolormap==1){
        rgbi=rgb_volsmokecolormap+4*i;
      }
      else{
        rgbi=cbi->colorbar+3*i;
      }
      glColor3fv(rgbi);
      zbot=(float)i/(float)ncolors;
      ztop=(float)(i+1)/(float)ncolors;

      glVertex3f(1.1,0.0,zbot);
      glVertex3f(1.3,0.0,zbot);
      glVertex3f(1.3,0.0,ztop);

      glVertex3f(1.1,0.0,zbot);
      glVertex3f(1.3,0.0,ztop);
      glVertex3f(1.3,0.0,zbot);

      glVertex3f(1.1,0.0,zbot);
      glVertex3f(1.3,0.0,ztop);
      glVertex3f(1.1,0.0,ztop);

      glVertex3f(1.1,0.0,zbot);
      glVertex3f(1.1,0.0,ztop);
      glVertex3f(1.3,0.0,ztop);

      glVertex3f(1.2,-0.1,zbot);
      glVertex3f(1.2, 0.1,zbot);
      glVertex3f(1.2, 0.1,ztop);

      glVertex3f(1.2,-0.1,zbot);
      glVertex3f(1.2, 0.1,ztop);
      glVertex3f(1.2, 0.1,zbot);

      glVertex3f(1.2,-0.1,zbot);
      glVertex3f(1.2, 0.1,ztop);
      glVertex3f(1.2,-0.1,ztop);

      glVertex3f(1.2,-0.1,zbot);
      glVertex3f(1.2,-0.1,ztop);
      glVertex3f(1.2, 0.1,ztop);
    }
    glEnd();
  }
}

/* ------------------ getcolorbar ------------------------ */

colorbardata *getcolorbar(char *label){
  int i;

  for(i=0;i<ncolorbars;i++){
    colorbardata *cb;

    cb = colorbarinfo + i;
    if(strcmp(cb->label,label)==0)return cb;
  }
  return NULL;
}

/* ------------------ UpdateCurrentColorbar ------------------------ */
#define FILEUPDATE 6
void UpdateCurrentColorbar(colorbardata *cb){
  int jj=0,fed_loaded=0;

  current_colorbar = cb;
  if(current_colorbar != NULL&&strcmp(current_colorbar->label, "FED") == 0){
    is_fed_colorbar = 1;
  }
  else{
    is_fed_colorbar = 0;
  }
  for(jj=0;jj<nslice_loaded;jj++){
    slicedata *slicej;
    int j;

    j = slice_loaded_list[jj];
    slicej = sliceinfo + j;
    if(slicej->display==0)continue;
    if(slicej->is_fed==1){
      fed_loaded=1;
      break;
    }
  }
  if(is_fed_colorbar==1&&fed_loaded==1)Slice_CB(FILEUPDATE);
}

/* ------------------ remapcolorbar ------------------------ */

void remapcolorbar(colorbardata *cbi){
  int i;
  float *colorbar;
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
    int i1,i2,j;

    i1 = cbi->index_node[i];
    i2 = cbi->index_node[i+1];
    if(i2==i1)continue;
    rgb_node = cbi->rgb_node+3*i;
    for(j=i1;j<i2;j++){
      float factor;

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
  if(show_extreme_mindata==1){
    colorbar[0]=rgb_below_min[0];
    colorbar[1]=rgb_below_min[1];
    colorbar[2]=rgb_below_min[2];
  }
  if(show_extreme_maxdata==1){
    colorbar[0+3*255]=rgb_above_max[0];
    colorbar[1+3*255]=rgb_above_max[1];
    colorbar[2+3*255]=rgb_above_max[2];
  }
  CheckMemory;
}

/* ------------------ remap_colorbartype ------------------------ */

void remap_colorbartype(int cb_oldtype, char *cb_newname){
  switch(cb_oldtype){
    case 0:
      strcpy(cb_newname,"Rainbow");
      break;
    case 1:
      strcpy(cb_newname,"Rainbow 2");
      break;
    case 2:
      strcpy(cb_newname,"yellow->red");
      break;
    case 3:
      strcpy(cb_newname,"blue->green->red");
      break;
    case 4:
      strcpy(cb_newname,"blue->red split");
      break;
    case 5:
      strcpy(cb_newname,"FED");
      break;
    case 6:
      //strcpy(cb_newname,"fire (original)");
      strcpy(cb_newname,"fire 2");
      break;
    case 7:
     // strcpy(cb_newname,"fire (black->orange)");
      strcpy(cb_newname,"fire 2");
      break;
    case 8:
      //strcpy(cb_newname,"fire (new)");
      strcpy(cb_newname,"fire 2");
      break;
    case 9:
      //strcpy(cb_newname,"fire (new2)");
      strcpy(cb_newname,"fire 2");
      break;
    case 10:
      //strcpy(cb_newname,"fire (custom)");
      strcpy(cb_newname,"fire 2");
      break;
    case 11:
      strcpy(cb_newname,"fire line (level set)");
      break;
    case 12:
      strcpy(cb_newname,"fire line (wall thickness)");
      break;
    case 13:
      strcpy(cb_newname,"black->white");
      break;
    default:
#define NCOLORBARS_PREV 14
      if(cb_oldtype>=NCOLORBARS_PREV){
        cb_oldtype -= (NCOLORBARS_PREV-ndefaultcolorbars);
      }
      if(cb_oldtype>=0&&cb_oldtype<ncolorbars){
        colorbardata *cb;

        cb = colorbarinfo + cb_oldtype;
        strcpy(cb_newname,cb->label);
      }
      else{
        strcpy(cb_newname,"Rainbow");
      }
      break;
  }
}

/* ------------------ initdefaultcolorbars ------------------------ */

void initdefaultcolorbars(void){
  int i;
  colorbardata *cbi;

  ndefaultcolorbars=11;

  FREEMEMORY(colorbarinfo);
  ncolorbars=ndefaultcolorbars;
  NewMemory((void **)&colorbarinfo,ncolorbars*sizeof(colorbardata));
  UpdateCurrentColorbar(colorbarinfo + colorbartype);


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

  // Rainbow 2 colorbar

  strcpy(cbi->label,"Rainbow 2");
  cbi->label_ptr=cbi->label;
  cbi->nnodes=4;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=0;
  cbi->rgb_node[1]=0;
  cbi->rgb_node[2]=255;

  cbi->index_node[1]=85;
  cbi->rgb_node[3]=0;
  cbi->rgb_node[4]=255;
  cbi->rgb_node[5]=255;

  cbi->index_node[2]=170;
  cbi->rgb_node[6]=255;
  cbi->rgb_node[7]=255;
  cbi->rgb_node[8]=0;

  cbi->index_node[3]=255;
  cbi->rgb_node[9]=255;
  cbi->rgb_node[10]=0;
  cbi->rgb_node[11]=0;
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

  // FED

  strcpy(cbi->label,"FED");
  cbi->label_ptr=cbi->label;

  cbi->nnodes=6;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=96;
  cbi->rgb_node[1]=96;
  cbi->rgb_node[2]=255;

  cbi->index_node[1]=26; // 0.295276,0.307087
  cbi->rgb_node[3]=96;
  cbi->rgb_node[4]=96;
  cbi->rgb_node[5]=255;

  cbi->index_node[2]=26;
  cbi->rgb_node[6]=255;
  cbi->rgb_node[7]=255;
  cbi->rgb_node[8]=0;

  cbi->index_node[3]=85; // 0.992126,1.003937
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
  strcpy(cbi->label,"fire");
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

  // fire 2

  fire_colorbar_index=cbi-colorbarinfo;
  fire_colorbar=cbi;
  strcpy(cbi->label,"fire 2");
  cbi->label_ptr=cbi->label;

  cbi->nnodes=10;
  cbi->nodehilight=0;

  cbi->index_node[0]=0;
  cbi->rgb_node[0]=0;
  cbi->rgb_node[1]=0;
  cbi->rgb_node[2]=0;

  cbi->index_node[1]=127;
  cbi->rgb_node[3]=38;
  cbi->rgb_node[4]=0;
  cbi->rgb_node[5]=0;

  cbi->index_node[2]=128;
  cbi->rgb_node[6]=219;
  cbi->rgb_node[7]=68;
  cbi->rgb_node[8]=21;

  cbi->index_node[3]=160;
  cbi->rgb_node[9]=255;
  cbi->rgb_node[10]=125;
  cbi->rgb_node[11]=36;

  cbi->index_node[4]=183;
  cbi->rgb_node[12]=255;
  cbi->rgb_node[13]=157;
  cbi->rgb_node[14]=52;

  cbi->index_node[5]=198;
  cbi->rgb_node[15]=255;
  cbi->rgb_node[16]=170;
  cbi->rgb_node[17]=63;

  cbi->index_node[6]=214;
  cbi->rgb_node[18]=255;
  cbi->rgb_node[19]=198;
  cbi->rgb_node[20]=93;

  cbi->index_node[7]=229;
  cbi->rgb_node[21]=255;
  cbi->rgb_node[22]=208;
  cbi->rgb_node[23]=109;

  cbi->index_node[8]=244;
  cbi->rgb_node[24]=255;
  cbi->rgb_node[25]=234;
  cbi->rgb_node[26]=161;

  cbi->index_node[9]=255;
  cbi->rgb_node[27]=255;
  cbi->rgb_node[28]=255;
  cbi->rgb_node[29]=238;
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

// construct colormaps from color node info

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

/* ------------------ drawColorBars ------------------------ */

void drawColorBars(void){
  int ilabel=0;

  int i,i3;
  int ileft=0;
  int leftzone, leftsmoke, leftslice, leftpatch, leftiso;
  int iposition;

  int sliceflag=0;
  int isoflag=0;
  float *slicefactor=NULL;
  float slicefactor2[2];
  float *isofactor=NULL;

  int plot3dflag=0;
  float *plot3dfactor=NULL;
  float plot3dfactor2[2];
  float plot3drange;

  int patchflag=0;
  int zoneflag=0;
  float *patchfactor=NULL;
  float *zonefactor=NULL;
  float patchrange=0.0;
  float zonerange;

  int partflag=0;
  float *partfactor=NULL;
  float partrange=0.0;

  int fed_slice=0;

  GLfloat *foreground_color, *red_color;

  foreground_color=&(foregroundcolor[0]);
  red_color=&(redcolor[0]);

  if(showiso_colorbar==1||showevac_colorbar==1||
    (showsmoke==1&&parttype!=0)||showslice==1||
    (showvslice==1&&vslicecolorbarflag==1)||
    (showpatch==1&&wc_flag==0)||
    (showzone==1&&zonecolortype==ZONETEMP_COLOR)||
    showplot3d==1){

    SNIFF_ERRORS("before colorbar");
    CheckMemory;
    if(showslice==1||(showvslice==1&&vslicecolorbarflag==1)){
      boundsdata *sb;

      sb = slicebounds + islicetype;

      if(strcmp(sb->label->shortlabel,"FED")==0){
          if(current_colorbar!=NULL){
            strcpy(default_fed_colorbar,current_colorbar->label);
            if(strcmp(current_colorbar->label,"FED")==0){
              fed_slice=1;
              if(strcmp(sb->colorlabels[1],"0.00")!=0||strcmp(sb->colorlabels[nrgb-1],"3.00")!=0)fed_slice=0;
            }
        }
      }
    }

    // -------------- draw plot3d colorbars ------------

    if(showplot3d==1&&contour_type==STEPPED_CONTOURS){
      glBegin(GL_QUADS);
      for (i = 0; i < nrgb-2; i++){
        float *rgb_plot3d_local;
        float ybot, ytop;

        rgb_plot3d_local = rgb_plot3d_contour[i];

        ybot = MIX2(i,nrgb-3,colorbar_top_pos,colorbar_down_pos);
        ytop = MIX2(i+1,nrgb-3,colorbar_top_pos,colorbar_down_pos);

        if(rgb_plot3d_local[3]!=0.0){
          glColor4fv(rgb_plot3d_local);
          glVertex2f((float)colorbar_left_pos, ybot);
          glVertex2f(colorbar_right_pos,ybot);

          glVertex2f(colorbar_right_pos,ytop);
          glVertex2f(colorbar_left_pos, ytop);
        }
      }
      glEnd();
      if(show_extreme_mindata==1||show_extreme_maxdata==1){
        float barmid;
        float *rgb_plot3d_local;
        float ybot, ytop;

        rgb_plot3d_local = rgb_plot3d_contour[nrgb-2];
        barmid = (colorbar_left_pos+colorbar_right_pos)/2.0;
        i=-1;
        ytop = MIX2(i+0.5,nrgb-3,colorbar_top_pos,colorbar_down_pos);
        ybot = MIX2(i+1,nrgb-3,colorbar_top_pos,colorbar_down_pos);

        if(have_extreme_mindata==1||have_extreme_maxdata==1)glEnable(GL_POLYGON_SMOOTH);

        if(show_extreme_mindata==1&&have_extreme_mindata==1&&rgb_plot3d_local[3]!=0.0){
          glBegin(GL_TRIANGLES);
          glColor4fv(rgb_plot3d_local);

          glVertex2f(colorbar_left_pos,ybot);
          glVertex2f(barmid,ytop);
          glVertex2f(colorbar_right_pos,ybot);
          glEnd();
        }

        i=nrgb-2;
        ybot = MIX2(i,nrgb-3,colorbar_top_pos,colorbar_down_pos);
        ytop = MIX2(i+0.5,nrgb-3,colorbar_top_pos,colorbar_down_pos);

        rgb_plot3d_local = rgb_plot3d_contour[nrgb-1];
        if(show_extreme_maxdata==1&&have_extreme_maxdata==1&&rgb_plot3d_local[3]!=0.0){
          glBegin(GL_TRIANGLES);
          glColor4fv(rgb_plot3d_local);
          glVertex2f(colorbar_left_pos, ybot);
          glVertex2f(colorbar_right_pos,ybot);
          glVertex2f(barmid, ytop);
          glEnd();
        }
        if(have_extreme_mindata==1||have_extreme_maxdata==1)glDisable(GL_POLYGON_SMOOTH);
      }
    }
    else{

      // -------------- draw all other colorbars ------------

      if(show_fed_area==1&&fed_slice==1&&fed_areas!=NULL){
        char area_label[256];
        char percen[]="%";
        float yy;

        glPushMatrix();
        glTranslatef(
          colorbar_left_pos,
          0.0,
          0.0);
        sprintf(area_label,"%i%s",fed_areas[0],percen);
        yy = MIX2(0.15,3.0,colorbar_top_pos,colorbar_down_pos)-VP_colorbar.text_height/2;
        outputBarText(0.0,yy,foreground_color,area_label);

        sprintf(area_label,"%i%s",fed_areas[1],percen);
        yy = MIX2(0.65,3.0,colorbar_top_pos,colorbar_down_pos)-VP_colorbar.text_height/2;
        outputBarText(0.0,yy,foreground_color,area_label);

        sprintf(area_label,"%i%s",fed_areas[2],percen);
        yy = MIX2(2.0,3.0,colorbar_top_pos,colorbar_down_pos)-VP_colorbar.text_height/2;
        outputBarText(0.0,yy,foreground_color,area_label);

        sprintf(area_label,"%i%s",fed_areas[3],percen);
        yy = MIX2(3.0,3.0,colorbar_top_pos,colorbar_down_pos)-VP_colorbar.text_height/2;
        outputBarText(0.0,yy+10,foreground_color,area_label);
        glPopMatrix();
      }

      glBegin(GL_QUADS);
      for (i = 0; i < nrgb_full-1; i++){
        float *rgb_cb,*rgb_cb2;
        float yy, yy2;

        rgb_cb=rgb_full[i];

        yy = MIX2(i,255,colorbar_top_pos,colorbar_down_pos);
        yy2 = MIX2(i+1,255,colorbar_top_pos,colorbar_down_pos);
        i3=i+1;
        if(i==nrgb_full-2)i3=i;
        rgb_cb2=rgb_full[i3];

        if(rgb_cb[3]!=0.0&&rgb_cb2[3]!=0.0){
          glColor4fv(rgb_cb);
          glVertex2f(colorbar_left_pos, yy);
          glVertex2f(colorbar_right_pos,yy);

          glColor4fv(rgb_cb2);
          glVertex2f(colorbar_right_pos,yy2);
          glVertex2f(colorbar_left_pos,yy2);
        }
      }
      glEnd();
    }
    if(show_extreme_mindata==1||show_extreme_maxdata==1){
      float barmid;

      barmid=(colorbar_right_pos+colorbar_left_pos)/2.0;

      if(have_extreme_mindata==1||have_extreme_maxdata==1)glEnable(GL_POLYGON_SMOOTH);

      if(show_extreme_mindata==1&&have_extreme_mindata==1){
        glBegin(GL_TRIANGLES);
        glColor4fv(rgb_full[0]);

        glVertex2f( colorbar_left_pos, colorbar_down_pos);
        glVertex2f(            barmid, colorbar_down_pos-0.866*colorbar_delta);
        glVertex2f(colorbar_right_pos, colorbar_down_pos);
        glEnd();
      }

      if(show_extreme_maxdata==1&&have_extreme_maxdata==1){
        glBegin(GL_TRIANGLES);
        glColor4fv(rgb_full[nrgb_full-1]);
        glVertex2f(colorbar_right_pos, colorbar_top_pos);
        glVertex2f(            barmid, colorbar_top_pos+0.866*colorbar_delta);
        glVertex2f( colorbar_left_pos, colorbar_top_pos);
        glEnd();
      }
      if(have_extreme_mindata==1||have_extreme_maxdata==1)glDisable(GL_POLYGON_SMOOTH);
    }
  }

  // -------------- compute columns where left labels will occur ------------

  leftsmoke=0;
  leftslice=0;
  leftpatch=0;
  leftiso=0;
  ileft=0;
  if(showiso_colorbar==1){
    leftiso=ileft;
    ileft++;
  }
  if(showevac_colorbar==1||showsmoke==1){
    if(parttype!=0){
      leftsmoke=ileft;
      ileft++;
    }
  }
  if(showslice==1||(showvslice==1&&vslicecolorbarflag==1)){
    leftslice=ileft;
    ileft++;
  }
  if(showpatch==1&&wc_flag==0){
    leftpatch=ileft;
  }
  leftzone = ileft;

  // -------------- particle file top labels ------------

  if(showevac_colorbar==1||showsmoke==1){
    char partunitlabel2[256], partshortlabel2[256];

    strcpy(partshortlabel2,"");
    strcpy(partunitlabel2,"");

    glPushMatrix();
    glTranslatef(
      colorbar_left_pos-colorbar_label_width,
      colorbar_top_pos+v_space+colorbar_delta,
      0.0);
    if(parttype!=0){
      if(showsmoke==1&&showevac==0)outputBarText(0.0,3*(VP_colorbar.text_height+v_space),foreground_color,"Part");
      if(showevac==1)outputBarText(0.0,3*(VP_colorbar.text_height+v_space),foreground_color,"Human");
    }
    if(parttype==-1){
      strcpy(partshortlabel2,"temp");
      strcpy(partunitlabel2,degC);
    }
    else if(parttype==-2){
      strcpy(partshortlabel2,"HRRPUV");
      strcpy(partunitlabel2,"kW/m3");
    }
    else{
      if(partshortlabel!=NULL)strcpy(partshortlabel2,partshortlabel);
      if(partunitlabel!=NULL)strcpy(partunitlabel2,partunitlabel);
    }
    if(parttype!=0){
      int partunitclass, partunittype;

      getunitinfo(partunitlabel2,&partunitclass,&partunittype);
      if(partunitclass>=0&&partunitclass<nunitclasses){
        if(partunittype>=0){
          partflag=1;
          partfactor=unitclasses[partunitclass].units[partunittype].scale;
          strcpy(partunitlabel,unitclasses[partunitclass].units[partunittype].unit);
        }
      }
      outputBarText(0.0,2*(VP_colorbar.text_height+v_space),foreground_color,partshortlabel);
      outputBarText(0.0,(VP_colorbar.text_height+v_space),foreground_color,partunitlabel);
      outputBarText(0.0,0.0,foreground_color,partscale);
    }
    glPopMatrix();
  }

  // -------------- slice file top labels ------------

  if(showslice==1||(showvslice==1&&vslicecolorbarflag==1)){
    char unitlabel[256];
    int sliceunitclass,sliceunittype;
    boundsdata *sb;

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
    glPushMatrix();
    glTranslatef(
      colorbar_left_pos-colorbar_label_width,
      colorbar_top_pos+v_space+colorbar_delta,
      0.0);
    glTranslatef(-leftslice*(colorbar_label_width+h_space),0.0,0.0);
    outputBarText(0.0,3*(VP_colorbar.text_height+v_space),foreground_color,"Slice");
    outputBarText(0.0,2*(VP_colorbar.text_height+v_space),foreground_color,sb->label->shortlabel);
    outputBarText(0.0,  (VP_colorbar.text_height+v_space),foreground_color,unitlabel);
    if(strcmp(unitlabel,"ppm")==0&&slicefactor!=NULL){
      slicefactor2[0]=*slicefactor*sb->fscale;
      slicefactor2[1]=0.0;
      slicefactor=slicefactor2;
    }
    else{
      outputBarText(0.0,0.0,foreground_color,sb->scale);
    }
    glPopMatrix();
    ilabel++;
  }

  // -------------- isosurface top labels ------------

  if(showiso_colorbar==1){
    char unitlabel[256];
    boundsdata *sb;

    sb = isobounds + iisottype;
    strcpy(unitlabel,sb->label->unit);
    glPushMatrix();
    glTranslatef(
      colorbar_left_pos-colorbar_label_width,
      colorbar_top_pos+v_space+colorbar_delta,
      0.0);
    glTranslatef(-leftiso*(colorbar_label_width+h_space),0.0,0.0);
    outputBarText(0.0,3*(VP_colorbar.text_height+v_space),foreground_color,"Iso");
    outputBarText(0.0,2*(VP_colorbar.text_height+v_space),foreground_color,sb->label->shortlabel);
    outputBarText(0.0,(VP_colorbar.text_height+v_space),foreground_color,unitlabel);
    outputBarText(0.0,0.0,foreground_color,sb->scale);
    glPopMatrix();
  }

  // -------------- boundary file top labels ------------

  if(showpatch==1&&wc_flag==0){
    char unitlabel[256];
    patchdata *patchi;
    int patchunitclass, patchunittype;

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
    glPushMatrix();
    glTranslatef(
      colorbar_left_pos-colorbar_label_width,
      colorbar_top_pos+v_space+colorbar_delta,
      0.0);
    glTranslatef(-leftpatch*(colorbar_label_width+h_space),0.0,0.0);
    outputBarText(0.0,3*(VP_colorbar.text_height+v_space),foreground_color,"Bndry");
    outputBarText(0.0,2*(VP_colorbar.text_height+v_space),foreground_color,patchi->label.shortlabel);
    outputBarText(0.0,  (VP_colorbar.text_height+v_space),foreground_color,unitlabel);
    outputBarText(0.0,0.0,foreground_color,patchi->scale);
    glPopMatrix();
  }

  // -------------- plot3d top labels ------------

  if(showplot3d==1){
    char *p3label;
    char *up3label;
    char unitlabel[256];
    int plot3dunitclass, plot3dunittype;

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
    glPushMatrix();
    glTranslatef(
      colorbar_left_pos-colorbar_label_width,
      colorbar_top_pos+v_space+colorbar_delta,
      0.0);
    outputBarText(0.0,3*(VP_colorbar.text_height+v_space),foreground_color,"Plot3d");
    outputBarText(0.0,2*(VP_colorbar.text_height+v_space),foreground_color,p3label);
    outputBarText(0.0, (VP_colorbar.text_height+v_space),foreground_color,unitlabel);
    if(strcmp(unitlabel,"ppm")==0&&plot3dfactor!=NULL){
      plot3dfactor2[0]=*plot3dfactor*fscalep3[plotn-1];
      plot3dfactor2[1]=0.0;
      plot3dfactor=plot3dfactor2;
    }
    else{
      outputBarText(0.0,0.0,foreground_color,scalep3[plotn-1]);
    }
    glPopMatrix();
  }
  if(showzone==1&&zonecolortype==ZONETEMP_COLOR){
    char unitlabel[256];
    int zoneunitclass, zoneunittype;

    strcpy(unitlabel,degC);
    getunitinfo(unitlabel,&zoneunitclass,&zoneunittype);
    if(zoneunitclass>=0&&zoneunitclass<nunitclasses){
      if(zoneunittype>0){
        zoneflag=1;
        zonefactor=unitclasses[zoneunitclass].units[zoneunittype].scale;
        strcpy(unitlabel,unitclasses[zoneunitclass].units[zoneunittype].unit);
      }
    }
    glPushMatrix();
    glTranslatef(colorbar_left_pos-colorbar_label_width,colorbar_top_pos+v_space+colorbar_delta,0.0);
    glTranslatef(-leftzone*(colorbar_label_width+h_space),0.0,0.0);
    outputBarText(0.0,3*(VP_colorbar.text_height+v_space),foreground_color,"Zone");
    outputBarText(0.0,2*(VP_colorbar.text_height+v_space),foreground_color,"Temp");
    outputBarText(0.0,(VP_colorbar.text_height+v_space),foreground_color,unitlabel);
    outputBarText(0.0,0.0,foreground_color,zonescale);
    glPopMatrix();
    SNIFF_ERRORS("After ZONE labels");
  }

  // -------------- isosurface left labels ------------

  if(showiso_colorbar==1){
    float tttval, tttmin, tttmax;
    boundsdata *sb;
    float isorange;

    sb = isobounds + iisottype;
    tttmin = sb->levels256[0];
    tttmax = sb->levels256[255];
    isorange=tttmax-tttmin;
    iposition=-1;
    glPushMatrix();
    glTranslatef(colorbar_left_pos-colorbar_label_width,-VP_colorbar.text_height/2.0,0.0);
    glTranslatef(-leftiso*(colorbar_label_width+h_space),0.0,0.0);
    if(global_colorbar_index!=-1){
      char isocolorlabel[256],isolabel[256];
      char *isocolorlabel_ptr=NULL;
      float vert_position;

      tttval = sb->levels256[valindex];
      num2string(isolabel,tttval,isorange);
      isocolorlabel_ptr=isolabel;
      if(isoflag==1){
        scalefloat2string(tttval,isocolorlabel, isofactor, isorange);
        isocolorlabel_ptr=isocolorlabel;
      }
      vert_position = MIX2(global_colorbar_index,255,colorbar_top_pos,colorbar_down_pos);
      iposition = MIX2(global_colorbar_index,255,nrgb-1,0);
      outputBarText(0.0,vert_position,red_color,isocolorlabel_ptr);
    }
    for (i=0; i<nrgb-1; i++){
      float vert_position;
      char isocolorlabel[256];
      char *isocolorlabel_ptr=NULL;

      vert_position = MIX2(i,nrgb-2,colorbar_top_pos,colorbar_down_pos);
      if(iposition==i)continue;
      isocolorlabel_ptr=&(sb->colorlabels[i+1][0]);
      if(isoflag==1){
        float val;

        val = tttmin + i*isorange/(nrgb-2);
        scalefloat2string(val,isocolorlabel, isofactor, isorange);
        isocolorlabel_ptr=isocolorlabel;
      }
      outputBarText(0.0,vert_position,foreground_color,isocolorlabel_ptr);
    }
    glPopMatrix();
  }

  // -------------- particle left labels ------------

  if(showevac_colorbar==1||(showsmoke==1&&parttype!=0)){
    float *partlevels256_ptr;
    float tttval, tttmin, tttmax;

    partlevels256_ptr=partlevels256;
    if(prop_index>=0&&prop_index<npart5prop){
      partlevels256_ptr=part5propinfo[prop_index].ppartlevels256;
    }

    iposition=-1;
    tttmin = partlevels256_ptr[0];
    tttmax = partlevels256_ptr[255];
    partrange = tttmax - tttmin;
    glPushMatrix();
    glTranslatef(colorbar_left_pos-colorbar_label_width,-VP_colorbar.text_height/2.0,0.0);
    glTranslatef(-leftsmoke*(colorbar_label_width+h_space),0.0,0.0);
    if(global_colorbar_index!=-1){
      char partcolorlabel[256],*partcolorlabel_ptr=NULL,partlabel[256];
      float vert_position;

      tttval = partlevels256_ptr[valindex];
      num2string(partlabel,tttval,partrange);
      partcolorlabel_ptr=partlabel;
      if(partflag==1){
        scalefloat2string(tttval,partcolorlabel, partfactor, partrange);
        partcolorlabel_ptr=partcolorlabel;
      }
      vert_position = MIX2(global_colorbar_index,255,colorbar_top_pos,colorbar_down_pos);
      iposition = MIX2(global_colorbar_index,255,nrgb-1,0);
      outputBarText(0.0,vert_position,red_color,partcolorlabel_ptr);
    }
    for (i=0; i<nrgb-1; i++){
      float vert_position;
      char partcolorlabel[256];
      char *partcolorlabel_ptr=NULL;

      vert_position = MIX2(i,nrgb-2,colorbar_top_pos,colorbar_down_pos);
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
        float val;

        val = tttmin + i*partrange/(nrgb-2);
        scalefloat2string(val,partcolorlabel, partfactor, partrange);
        scalestring(partcolorlabel_ptr,partcolorlabel, partfactor, partrange);
        partcolorlabel_ptr=partcolorlabel;
      }
      outputBarText(0.0,vert_position,foreground_color,partcolorlabel_ptr);
    }
    glPopMatrix();
  }

  // -------------- slice left labels ------------

  if(showslice==1||(showvslice==1&&vslicecolorbarflag==1)){
    float tttval, tttmin, tttmax;
    boundsdata *sb;
    float slicerange;

    sb=slicebounds+islicetype;
    tttmin = sb->levels256[0];
    tttmax = sb->levels256[255];
    slicerange=tttmax-tttmin;
    iposition=-1;
    glPushMatrix();
    glTranslatef(colorbar_left_pos-colorbar_label_width,-VP_colorbar.text_height/2.0,0.0);
    glTranslatef(-leftslice*(colorbar_label_width+h_space),0.0,0.0);
    if(global_colorbar_index!=-1){
      char slicelabel[256], slicecolorlabel[256];
      char *slicecolorlabel_ptr=NULL;
      float vert_position;

      tttval = sb->levels256[valindex];
      num2string(slicelabel,tttval,slicerange);
      slicecolorlabel_ptr=slicelabel;
      if(sliceflag==1){
        scalefloat2string(tttval,slicecolorlabel, slicefactor, slicerange);
        slicecolorlabel_ptr=slicecolorlabel;
      }
      vert_position = MIX2(global_colorbar_index,255,colorbar_top_pos,colorbar_down_pos);
      iposition = MIX2(global_colorbar_index,255,nrgb-1,0);
      outputBarText(0.0,vert_position,red_color,slicecolorlabel_ptr);
    }
    if(fed_slice==1){
      for (i=0; i<nrgb-1; i++){
        float vert_position;

        vert_position = MIX2(0.0,3.0,colorbar_top_pos,colorbar_down_pos);
        outputBarText(0.0,vert_position,foreground_color,"0.00");

        vert_position = MIX2(0.3,3.0,colorbar_top_pos,colorbar_down_pos);
        outputBarText(0.0,vert_position,foreground_color,"0.30");

        vert_position = MIX2(1.0,3.0,colorbar_top_pos,colorbar_down_pos);
        outputBarText(0.0,vert_position,foreground_color,"1.00");

        vert_position = MIX2(3.0,3.0,colorbar_top_pos,colorbar_down_pos);
        outputBarText(0.0,vert_position,foreground_color,"3.00");
      }
    }
    else{
      for (i=0; i<nrgb-1; i++){
        float vert_position;
        char slicecolorlabel[256];
        char *slicecolorlabel_ptr=NULL;

        vert_position = MIX2(i,nrgb-2,colorbar_top_pos,colorbar_down_pos);
        if(iposition==i)continue;
        slicecolorlabel_ptr=&(sb->colorlabels[i+1][0]);
        if(sliceflag==1){
          float val;

          val = tttmin + i*slicerange/(nrgb-2);
          scalefloat2string(val,slicecolorlabel, slicefactor, slicerange);
          slicecolorlabel_ptr=slicecolorlabel;
        }
        outputBarText(0.0,vert_position,foreground_color,slicecolorlabel_ptr);
      }
    }
    glPopMatrix();
  }

  // -------------- boundary left labels ------------

  if(showpatch==1&&wc_flag==0){
    float tttval, tttmin, tttmax;

    iposition=-1;
    tttmin = boundarylevels256[0];
    tttmax = boundarylevels256[255];
    patchrange=tttmax-tttmin;
    glPushMatrix();
    glTranslatef(colorbar_left_pos-colorbar_label_width,-VP_colorbar.text_height/2.0,0.0);
    glTranslatef(-leftpatch*(colorbar_label_width+h_space),0.0,0.0);
    if(global_colorbar_index!=-1){
      char patchcolorlabel[256],boundarylabel[256],*patchcolorlabel_ptr=NULL;
      float vert_position;

      // draw boundary file value selected with mouse
      tttval = boundarylevels256[valindex];
      num2string(boundarylabel,tttval,tttmax-tttmin);
      patchcolorlabel_ptr=&(boundarylabel[0]);
      if(patchflag==1){
        scalefloat2string(tttval,patchcolorlabel, patchfactor, patchrange);
        patchcolorlabel_ptr=patchcolorlabel;
      }
      vert_position = MIX2(global_colorbar_index,255,colorbar_top_pos,colorbar_down_pos);
      iposition = MIX2(global_colorbar_index,255,nrgb-1,0);
      outputBarText(0.0,vert_position,red_color,patchcolorlabel_ptr);
    }
    for (i=0; i<nrgb-1; i++){
      char patchcolorlabel[256];
      char *patchcolorlabel_ptr=NULL;
      float vert_position;

      vert_position = MIX2(i,nrgb-2,colorbar_top_pos,colorbar_down_pos);

      if(iposition==i)continue;
      patchcolorlabel_ptr=&colorlabelpatch[i+1][0];
      if(patchflag==1){
        float val;

        val = tttmin + i*patchrange/(nrgb-2);
        scalefloat2string(val,patchcolorlabel, patchfactor, patchrange);
        patchcolorlabel_ptr=patchcolorlabel;
      }
      outputBarText(0.0,vert_position,foreground_color,patchcolorlabel_ptr);
    }
    glPopMatrix();
  }

  // -------------- zone left labels ------------

  if(showzone==1&&zonecolortype==ZONETEMP_COLOR){
    float tttval, tttmin, tttmax;

    iposition=-1;
    tttmin = zonelevels256[0];
    tttmax = zonelevels256[255];
    zonerange=tttmax-tttmin;
    glPushMatrix();
    glTranslatef(colorbar_left_pos-colorbar_label_width,-VP_colorbar.text_height/2.0,0.0);
    glTranslatef(-leftzone*(colorbar_label_width+h_space),0.0,0.0);
    if(global_colorbar_index!=-1){
      char zonecolorlabel[256],*zonecolorlabel_ptr=NULL,zonelabel[256];
      float vert_position;

      tttval = zonelevels256[valindex];
      num2string(zonelabel,tttval,tttmax-tttmin);
      zonecolorlabel_ptr=&(zonelabel[0]);
      if(zoneflag==1){
        scalefloat2string(tttval,zonecolorlabel, zonefactor, zonerange);
        zonecolorlabel_ptr=zonecolorlabel;
      }
      vert_position = MIX2(global_colorbar_index,255,colorbar_top_pos,colorbar_down_pos);
      iposition = MIX2(global_colorbar_index,255,nrgb-1,0);
      outputBarText(0.0,vert_position,red_color,zonecolorlabel_ptr);
    }
    for (i=0; i<nrgb-1; i++){
      float vert_position;
      char zonecolorlabel[256];
      char *zonecolorlabel_ptr=NULL;

      vert_position = MIX2(i,nrgb-2,colorbar_top_pos,colorbar_down_pos);
      if(iposition==i)continue;
      zonecolorlabel_ptr=&colorlabelzone[i+1][0];
      if(zoneflag==1){
        float val;

        val = tttmin + (i-1)*zonerange/(nrgb-2);
        scalefloat2string(val,zonecolorlabel, zonefactor, zonerange);
        zonecolorlabel_ptr=zonecolorlabel;
      }
      outputBarText(0.0,vert_position,foreground_color,zonecolorlabel_ptr);
    }
    SNIFF_ERRORS("after zone left labels");
    glPopMatrix();
  }

  // -------------- plot3d left labels ------------

  if(showplot3d==1){
    float *p3lev;
    float tttval, tttmin, tttmax;

    iposition=-1;
    p3lev = p3levels256[plotn-1];
    tttmin = p3lev[0];
    tttmax = p3lev[255];
    plot3drange = tttmax - tttmin;
    glPushMatrix();
    glTranslatef(colorbar_left_pos-colorbar_label_width,-VP_colorbar.text_height/2.0,0.0);
    if(global_colorbar_index!=-1){
      char plot3dcolorlabel[256], p3dlabel[256], *plot3dcolorlabel_ptr=NULL;
      float vert_position;

      tttval = p3lev[valindex];
      num2string(p3dlabel,tttval,tttmax-tttmin);
      plot3dcolorlabel_ptr = p3dlabel;
      if(plot3dflag==1){
        scalefloat2string(tttval,plot3dcolorlabel, plot3dfactor, plot3drange);
        plot3dcolorlabel_ptr=plot3dcolorlabel;
      }
      vert_position = MIX2(global_colorbar_index,255,colorbar_top_pos,colorbar_down_pos);
      iposition = MIX2(global_colorbar_index,255,nrgb-1,0);
      outputBarText(0.0,vert_position,red_color,plot3dcolorlabel_ptr);
    }
    if(visiso==0){
      float vert_position;

      for (i=0; i<nrgb-1; i++){
        char plot3dcolorlabel[256];
        char *plot3dcolorlabel_ptr=NULL;

        vert_position = MIX2(i,nrgb-2,colorbar_top_pos,colorbar_down_pos);
        if(iposition==i)continue;
        plot3dcolorlabel_ptr=&colorlabelp3[plotn-1][i][0];
        if(plot3dflag==1){
          float val;

          val = tttmin + i*plot3drange/(nrgb-2);
          scalefloat2string(val,plot3dcolorlabel, plot3dfactor, plot3drange);
          plot3dcolorlabel_ptr=plot3dcolorlabel;
        }
        outputBarText(0.0,vert_position,foreground_color,plot3dcolorlabel_ptr);
      }
    }
    else{
      float vert_position;

      for (i=0; i<nrgb-2; i++){
        char plot3dcolorlabel[256];
        char *plot3dcolorlabel_ptr=NULL;

        vert_position = MIX2(i,nrgb-2,colorbar_top_pos,colorbar_down_pos);

        if(iposition==i)continue;
        plot3dcolorlabel_ptr=&colorlabeliso[plotn-1][i][0];
        if(plot3dflag==1){
          float val;

          val = tttmin + (i-1)*plot3drange/(nrgb-2);
          scalefloat2string(val,plot3dcolorlabel, plot3dfactor, plot3drange);
          plot3dcolorlabel_ptr=plot3dcolorlabel;
        }
        if(isolevelindex==i||isolevelindex2==i){
          outputBarText(0.0,vert_position,red_color,plot3dcolorlabel_ptr);
        }
        else{
          outputBarText(0.0,vert_position,foreground_color,plot3dcolorlabel_ptr);
        }
      }
    }
    glPopMatrix();
  }
}



