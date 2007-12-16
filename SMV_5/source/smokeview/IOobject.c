// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>  
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "flowfiles.h"
#include "MALLOC.h"
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

// svn revision character string
char IOobject_revision[]="$Revision$";

#define SV_TRANSLATE  100
#define SV_ROTATEX    101
#define SV_ROTATEY    102
#define SV_ROTATEZ    103
#define SV_SCALEXYZ   104
#define SV_SCALE      105
#ifdef pp_AVATAR
#define SV_GETUSERVALS    106
#endif

#define SV_TRANSLATE_NUMARGS  3
#define SV_ROTATEX_NUMARGS    1
#define SV_ROTATEY_NUMARGS    1
#define SV_ROTATEZ_NUMARGS    1
#define SV_SCALEXYZ_NUMARGS   3
#define SV_SCALE_NUMARGS      1
#ifdef pp_AVATAR
#define SV_GETUSERVALS_NUMARGS    1
#endif

#define SV_DRAWCUBE      200
#define SV_DRAWSPHERE    201
#define SV_DRAWDISK      202
#define SV_DRAWLINE      203
#define SV_DRAWCIRCLE    204
#define SV_DRAWTRUNCCONE 205
#define SV_DRAWNOTCHPLATE 206
#define SV_DRAWRING      207
#define SV_DRAWCONE      208
#define SV_DRAWHEXDISK   209
#define SV_DRAWPOLYDISK  210
#define SV_DRAWPOINT     211


#define SV_DRAWCUBE_NUMARGS      1
#define SV_DRAWSPHERE_NUMARGS    1
#define SV_DRAWDISK_NUMARGS      2
#define SV_DRAWLINE_NUMARGS      6
#define SV_DRAWCIRCLE_NUMARGS    1
#define SV_DRAWTRUNCCONE_NUMARGS 3
#define SV_DRAWNOTCHPLATE_NUMARGS 4
#define SV_DRAWRING_NUMARGS      3
#define SV_DRAWCONE_NUMARGS      2
#define SV_DRAWHEXDISK_NUMARGS   2
#define SV_DRAWPOLYDISK_NUMARGS   3
#define SV_DRAWPOINT_NUMARGS     0

#define SV_PUSH       300
#define SV_POP        301
#define SV_SETCOLOR   302
#define SV_SETBW      303
#define SV_SETLINEWIDTH 304
#define SV_SETPOINTSIZE 305

#define SV_NO_OP      999

#define SV_PUSH_NUMARGS       0
#define SV_POP_NUMARGS        0
#define SV_SETCOLOR_NUMARGS   3
#define SV_SETBW_NUMARGS      1
#define SV_SETLINEWIDTH_NUMARGS 1
#define SV_SETPOINTSIZE_NUMARGS 1

#define SV_ERR -1

void reporterror(char *buffer, char *token, int numargs_found, int numargs_expected);

void drawcone(float d1, float height, float *rgbcolor);
void drawtrunccone(float d1, float d2, float height, float *rgbcolor);
void drawline(float *xyz1, float *xyz2, float *rgbcolor);
void drawcircle(float diameter, float *rgbcolor);
void drawpoint(float *rgbcolor);
void drawsphere(float diameter, float *rgbcolor);
void drawcube(float size, float *rgbcolor);
void drawdisk(float diameter, float height, float *rgbcolor);
void drawhexdisk(float diameter, float height, float *rgbcolor);
void drawpolydisk(int nsides, float diameter, float height, float *rgbcolor);
void drawring(float d_inner, float d_outer, float height, float *rgbcolor);
void drawnotchplate(float diameter, float height, float notchheight, float direction, float *rgbcolor);
void draw_SVOBJECT(sv_object *object, int iframe);
sv_object *get_object(char *label);
void free_object(sv_object *object);
void remove_comment(char *buffer);
void freecircle(void);
void initcircle(unsigned int npoints);
void getargsops(char *buffer,float **args,int *nargs, int **ops, int *nops);

static float *xcirc=NULL, *ycirc=NULL;
static int ncirc;
static float specular[4]={0.4,0.4,0.4,1.0};

/* ----------------------- draw_devices ----------------------------- */

void draw_devices(void){
  device *devicei;
  int i;
  float *xyz;

  glPushMatrix();
  glPushAttrib(GL_POINT_BIT|GL_LINE_BIT);
  glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);
  glTranslatef(-xbar0,-ybar0,-zbar0);
  for(i=0;i<ndeviceinfo;i++){
    devicei = deviceinfo + i;

    //if((visSensor==1&&ntc_total>0)||
    //   (visSprink==1&&nspr_total>0)||
    //   (visHeat==1&&nheat_total>0)){

    if(devicei->object->visible==0)continue;
    if(isZoneFireModel==1&&strcmp(devicei->object->label,"target")==0&&visSensor==0)continue;
    xyz = devicei->xyz;
    glPushMatrix();
    glTranslatef(xyz[0],xyz[1],xyz[2]);

    if(devicei->angle_az!=0.0){
      glRotatef(devicei->angle_az,0.0,0.0,1.0);
    }
    if(devicei->angle_elev!=0.0){
      glRotatef(-devicei->angle_elev,0.0,1.0,0.0);
    }

    if(showtime==1&&itime>=0&&itime<ntimes&&
      times!=NULL&&devicei->act_time>=0.0&&times[itime]>=devicei->act_time&&
      devicei->object->nframes>1){
      draw_SVOBJECT(devicei->object,1);
    }
    else{
      draw_SVOBJECT(devicei->object,0);
    }
    glPopMatrix();
  }

  glPopAttrib();
  glPopMatrix();
  drawTargetNorm();
}

/* ----------------------- draw_devices ----------------------------- */

void drawTargetNorm(void){
  int i;
  device *devicei;
  float *xyz, *xyznorm;

  if(isZoneFireModel==1&&hasSensorNorm==1&&visSensor==1&&visSensorNorm==1){
    glPushMatrix();
    glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);
    glBegin(GL_LINES);
    glColor4fv(sensornormcolor);

    for(i=0;i<ndeviceinfo;i++){
      float xyz2[3];

      devicei = deviceinfo + i;

      if(devicei->object->visible==0)continue;
      if(strcmp(devicei->object->label,"sensor")==0&&visSensor==0)continue;
      if(isZoneFireModel==1&&strcmp(devicei->object->label,"target")==0&&visSensor==0)continue;
      xyz = devicei->xyz;
      xyznorm = devicei->xyznorm;
      glVertex3fv(xyz);
      xyz2[0]=xyz[0]+devicenorm_length*xyznorm[0];
      xyz2[1]=xyz[1]+devicenorm_length*xyznorm[1];
      xyz2[2]=xyz[2]+devicenorm_length*xyznorm[2];
      glVertex3fv(xyz2);
    }
    glEnd();
    glPopMatrix();
  }
}

/* ----------------------- draw_SMVOBJECT ----------------------------- */

void draw_SVOBJECT(sv_object *object, int iframe){
  sv_object_frame *framei;
  int *op;
  float *arg;
  int iarg,iop;
  float *rgbptr;
  float rgbcolor[4];

  framei=object->obj_frames[iframe];
  ASSERT(framei->error==0||framei->error==1);
  if(framei->error!=0)framei=error_frame;

  rgbcolor[0]=1.0;
  rgbcolor[1]=0.0;
  rgbcolor[2]=0.0;
  rgbcolor[3]=1.0;
  rgbptr=rgbcolor;
  glPushMatrix();
  iarg = 0;
  iop = 0;

  if(framei->display_list_ID==-1){
    int displaylist_id;

    displaylist_id = glGenLists(1);
    if(displaylist_id!=0){
      framei->display_list_ID=displaylist_id;
      glNewList(displaylist_id,GL_COMPILE_AND_EXECUTE);
    }

	glEnable(GL_LIGHTING);

	glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,specular);

    glEnable(GL_COLOR_MATERIAL);

    while(iop<framei->nops){
      arg = framei->args + iarg;
      op = framei->ops + iop;
      switch (*op){
#ifdef pp_AVATAR
      case SV_GETUSERVALS:
        if(iarg+SV_GETUSERVALS_NUMARGS<=framei->nargs){
          int iiarg;

          iiarg=arg[0]+0.5;
          topvalstack-=iiarg;
          if(topvalstack>=0&&iarg+1+iiarg<framei->nargs){
            int i;

            for(i=0;i<iiarg;i++){
              arg[1+i]=valstack[topvalstack+i];
            }
          }
        }
        iarg+=1;
        break;
#endif
      case SV_TRANSLATE:
        if(iarg+SV_TRANSLATE_NUMARGS<=framei->nargs)glTranslatef(arg[0],arg[1],arg[2]);
        iarg+=3;
        break;
      case SV_ROTATEX:
        if(iarg+SV_ROTATEX_NUMARGS<=framei->nargs)glRotatef(arg[0],1.0,0.0,0.0);
        iarg++;
        break;
      case SV_ROTATEY:
        if(iarg+SV_ROTATEY_NUMARGS<=framei->nargs)glRotatef(arg[0],0.0,1.0,0.0);
        iarg++;
        break;
      case SV_ROTATEZ:
        if(iarg+SV_ROTATEZ_NUMARGS<=framei->nargs)glRotatef(arg[0],0.0,0.0,1.0);
        iarg++;
        break;
      case SV_SCALEXYZ:
        if(iarg+SV_SCALEXYZ_NUMARGS<=framei->nargs)glScalef(arg[0],arg[1],arg[2]);
        iarg+=3;
        break;
      case SV_SCALE:
        if(iarg+SV_SCALE_NUMARGS<=framei->nargs)glScalef(arg[0],arg[1],arg[2]);
        iarg++;
        break;
      case SV_DRAWCUBE:
        if(iarg+SV_DRAWCUBE_NUMARGS<=framei->nargs)drawcube(arg[0],rgbptr);
        rgbptr=NULL;
        iarg++;
        break;
      case SV_DRAWDISK:
        if(iarg+SV_DRAWDISK_NUMARGS<=framei->nargs)drawdisk(arg[0],arg[1], rgbptr);
        rgbptr=NULL;
        iarg+=2;
        break;
      case SV_DRAWHEXDISK:
        if(iarg+SV_DRAWHEXDISK_NUMARGS<=framei->nargs)drawhexdisk(arg[0],arg[1], rgbptr);
        rgbptr=NULL;
        iarg+=2;
        break;
      case SV_DRAWPOLYDISK:
        if(iarg+SV_DRAWPOLYDISK_NUMARGS<=framei->nargs){
          int nsides;
  
          nsides = arg[0]+0.5;
          drawpolydisk(nsides, arg[1],arg[2], rgbptr);
          rgbptr=NULL;
        }
        iarg+=3;
        break;
      case SV_DRAWRING:
        if(iarg+SV_DRAWRING_NUMARGS<=framei->nargs)drawring(arg[0],arg[1], arg[2], rgbptr);
        rgbptr=NULL;
        iarg+=3;
        break;
      case SV_DRAWNOTCHPLATE:
        if(iarg+SV_DRAWNOTCHPLATE_NUMARGS<=framei->nargs)drawnotchplate(arg[0],arg[1], arg[2], arg[3], rgbptr);
        rgbptr=NULL;
        iarg+=4;
        break;
      case SV_DRAWTRUNCCONE:
        if(iarg+SV_DRAWTRUNCCONE_NUMARGS<=framei->nargs)drawtrunccone(arg[0],arg[1],arg[2], rgbptr);
        rgbptr=NULL;
        iarg+=3;
        break;
      case SV_DRAWCONE:
        if(iarg+SV_DRAWCONE_NUMARGS<=framei->nargs)drawcone(arg[0],arg[1], rgbptr);
        rgbptr=NULL;
        iarg+=2;
        break;
      case SV_DRAWSPHERE:
        if(iarg+SV_DRAWSPHERE_NUMARGS<=framei->nargs)drawsphere(arg[0],rgbptr);
        rgbptr=NULL;
        iarg++;
        break;
      case SV_DRAWCIRCLE:
        if(iarg+SV_DRAWCIRCLE_NUMARGS<=framei->nargs)drawcircle(arg[0],rgbptr);
        rgbptr=NULL;
        iarg++;
        break;
      case SV_DRAWPOINT:
        if(iarg+SV_DRAWPOINT_NUMARGS<=framei->nargs)drawpoint(rgbptr);
        rgbptr=NULL;
        break;
      case SV_SETCOLOR:
        if(iarg+SV_SETCOLOR_NUMARGS<=framei->nargs){
          rgbcolor[0]=arg[0];
          rgbcolor[1]=arg[1];
          rgbcolor[2]=arg[2];
          rgbcolor[3]=1.0;
          rgbptr=rgbcolor;
        }
        iarg+=3;
        break;
      case SV_SETLINEWIDTH:
        if(iarg+SV_SETLINEWIDTH_NUMARGS<=framei->nargs){
          glLineWidth(arg[0]);
          iarg++;
        }
        break;
      case SV_SETPOINTSIZE:
        if(iarg+SV_SETPOINTSIZE_NUMARGS<=framei->nargs){
          glPointSize(arg[0]);
          iarg++;
        }
        break;
      case SV_SETBW:
        if(iarg+SV_SETBW_NUMARGS<=framei->nargs){
          rgbcolor[0]=arg[0];
          rgbcolor[1]=arg[0];
          rgbcolor[2]=arg[0];
          rgbcolor[3]=1.0;
          rgbptr=rgbcolor;
        }
        iarg+=1;
        break;
      case SV_DRAWLINE:
        if(iarg+SV_DRAWLINE_NUMARGS<=framei->nargs)drawline(arg,arg+3,rgbptr);
        rgbptr=NULL;
        iarg+=6;
        break;
      case SV_PUSH:
        glPushMatrix();
        break;
      case SV_POP:
        glPopMatrix();
        break;
      case SV_NO_OP:
        break;
      case SV_ERR:
        break;
      default:
        ASSERT(FFALSE);
        break;
      }
      iop++;
    }
    if(displaylist_id!=0){
      glEndList();
    }

    glDisable(GL_COLOR_MATERIAL);
	glDisable(GL_LIGHTING);
  }
  else{
    glCallList(framei->display_list_ID);
  }
  glPopMatrix();

}

/* ----------------------- drawline ----------------------------- */

void drawline(float *xyz1, float *xyz2, float *rgbcolor){
  glBegin(GL_LINES);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);
  glVertex3fv(xyz1);
  glVertex3fv(xyz2);
  glEnd();
}


/* ----------------------- drawsphere ----------------------------- */

void drawsphere(float diameter, float *rgbcolor){
  int i,j;
  float *radsphere, *zsphere;

  if(ncirc==0)initcircle(12);
  radsphere=ycirc;
  zsphere=xcirc;
  glPushMatrix();
  glScalef(diameter/2.0,diameter/2.0,diameter/2.0);

  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);
  for(i=0;i<ncirc/2+1;i++){
    for(j=0;j<ncirc;j++){
      glNormal3f(radsphere[  i]*xcirc[  j],radsphere[  i]*ycirc[  j],zsphere[  i]);
      glVertex3f(radsphere[  i]*xcirc[  j],radsphere[  i]*ycirc[  j],zsphere[  i]);

      glNormal3f(radsphere[i+1]*xcirc[  j],radsphere[i+1]*ycirc[  j],zsphere[i+1]);
      glVertex3f(radsphere[i+1]*xcirc[  j],radsphere[i+1]*ycirc[  j],zsphere[i+1]);

      glNormal3f(radsphere[i+1]*xcirc[j+1],radsphere[i+1]*ycirc[j+1],zsphere[i+1]);
      glVertex3f(radsphere[i+1]*xcirc[j+1],radsphere[i+1]*ycirc[j+1],zsphere[i+1]);

      glNormal3f(radsphere[  i]*xcirc[j+1],radsphere[  i]*ycirc[j+1],zsphere[  i]);
      glVertex3f(radsphere[  i]*xcirc[j+1],radsphere[  i]*ycirc[j+1],zsphere[  i]);
    }
  }
  glEnd();
  glPopMatrix();
}

/* ----------------------- drawpoint ----------------------------- */

void drawpoint(float *rgbcolor){
  glBegin(GL_POINTS);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);
  glVertex3f(0.0,0.0,0.0);
  glEnd();
}

/* ----------------------- drawcircle ----------------------------- */

void drawcircle(float diameter,float *rgbcolor){
  int i;

  if(ncirc==0)initcircle(12);
  glBegin(GL_LINE_LOOP);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);
  for(i=0;i<ncirc;i++){
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0);
  }
  glEnd();
}

/* ----------------------- drawcube ----------------------------- */

void drawcube(float size, float *rgbcolor){
  float s2;

  s2 = size/2.0;

  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);


  glNormal3f(0.0,0.0,-1.0);
  glVertex3f(-s2,-s2,-s2);  // 1
  glVertex3f(-s2, s2,-s2);  // 4
  glVertex3f( s2, s2,-s2);  // 3
  glVertex3f( s2,-s2,-s2);  // 2

  glNormal3f(0.0,0.0,1.0);
  glVertex3f(-s2,-s2, s2);  // 5
  glVertex3f( s2,-s2, s2);  // 6
  glVertex3f( s2, s2, s2);  // 7
  glVertex3f(-s2, s2, s2);  // 8

  glNormal3f(0.0,-1.0,0.0);
  glVertex3f(-s2,-s2,-s2);  // 1
  glVertex3f( s2,-s2,-s2);  // 2
  glVertex3f( s2,-s2, s2);  // 6
  glVertex3f(-s2,-s2, s2);  // 5
                    
  glNormal3f(0.0,1.0,0.0);
  glVertex3f( s2, s2,-s2);  // 3
  glVertex3f(-s2, s2,-s2);  // 4
  glVertex3f(-s2, s2, s2);  // 8
  glVertex3f( s2, s2, s2);  // 7

  glNormal3f(-1.0,0.0,0.0);
  glVertex3f(-s2,-s2,-s2);  // 1
  glVertex3f(-s2,-s2, s2);  // 5
  glVertex3f(-s2, s2, s2);  // 8
  glVertex3f(-s2, s2,-s2);  // 4
                     
  glNormal3f(1.0,0.0,0.0);
  glVertex3f( s2,-s2,-s2);  // 2
  glVertex3f( s2, s2,-s2);  // 3
  glVertex3f( s2, s2, s2);  // 7
  glVertex3f( s2,-s2, s2);  // 6
  glEnd();

}



/* ----------------------- drawring ----------------------------- */

void drawring(float diam_inner, float diam_outer, float height, float *rgbcolor){
  int i;

  if(ncirc==0)initcircle(12);
  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);

  for(i=0;i<ncirc;i++){
    glNormal3f(xcirc[i],ycirc[i],0.0);
    glVertex3f(diam_outer*xcirc[  i]/2.0,diam_outer*ycirc[  i]/2.0,0.0); // 1

    glNormal3f(xcirc[i+1],ycirc[i+1],0.0);
    glVertex3f(diam_outer*xcirc[i+1]/2.0,diam_outer*ycirc[i+1]/2.0,0.0); // 2

    glNormal3f(xcirc[i+1],ycirc[i+1],0.0);
    glVertex3f(diam_outer*xcirc[i+1]/2.0,diam_outer*ycirc[i+1]/2.0, height); // 3

    glNormal3f(xcirc[i],ycirc[i],0.0);
    glVertex3f(diam_outer*xcirc[  i]/2.0,diam_outer*ycirc[  i]/2.0, height); // 4

    glNormal3f(xcirc[i],ycirc[i],0.0);
    glVertex3f(diam_inner*xcirc[  i]/2.0,diam_inner*ycirc[  i]/2.0,0.0); // 1

    glNormal3f(xcirc[i],ycirc[i],0.0);
    glVertex3f(diam_inner*xcirc[  i]/2.0,diam_inner*ycirc[  i]/2.0, height); // 4

    glNormal3f(xcirc[i+1],ycirc[i+1],0.0);
    glVertex3f(diam_inner*xcirc[i+1]/2.0,diam_inner*ycirc[i+1]/2.0, height); // 3

    glNormal3f(xcirc[i+1],ycirc[i+1],0.0);
    glVertex3f(diam_inner*xcirc[i+1]/2.0,diam_inner*ycirc[i+1]/2.0,0.0); // 2
  }
  glNormal3f(0.0,0.0,1.0);
  for(i=0;i<ncirc;i++){
    glVertex3f(diam_outer*xcirc[  i]/2.0,diam_outer*ycirc[  i]/2.0,height);
    glVertex3f(diam_outer*xcirc[i+1]/2.0,diam_outer*ycirc[i+1]/2.0,height);
    glVertex3f(diam_inner*xcirc[i+1]/2.0,diam_inner*ycirc[i+1]/2.0,height);
    glVertex3f(diam_inner*xcirc[  i]/2.0,diam_inner*ycirc[  i]/2.0,height);
  }
  glNormal3f(0.0,0.0,-1.0);
  for(i=0;i<ncirc;i++){
    glVertex3f(diam_outer*xcirc[  i]/2.0,diam_outer*ycirc[  i]/2.0,0.0);
    glVertex3f(diam_inner*xcirc[  i]/2.0,diam_inner*ycirc[  i]/2.0,0.0);
    glVertex3f(diam_inner*xcirc[i+1]/2.0,diam_inner*ycirc[i+1]/2.0,0.0);
    glVertex3f(diam_outer*xcirc[i+1]/2.0,diam_outer*ycirc[i+1]/2.0,0.0);
  }
  glEnd();

}

/* ----------------------- drawdisk ----------------------------- */

void drawdisk(float diameter, float height, float *rgbcolor){
  int i;

  if(ncirc==0)initcircle(12);
  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);

  for(i=0;i<ncirc;i++){
    glNormal3f(xcirc[i],ycirc[i],0.0);
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1

    glNormal3f(xcirc[i+1],ycirc[i+1],0.0);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,0.0); // 2

    glNormal3f(xcirc[i+1],ycirc[i+1],0.0);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height); // 3

    glNormal3f(xcirc[i],ycirc[i],0.0);
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height); // 4
  }
  glEnd();

  glBegin(GL_TRIANGLES);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);

  glNormal3f(0.0,0.0,-1.0);
  for(i=0;i<ncirc;i++){
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0);
    glVertex3f(                    0.0,                    0.0,0.0);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,0.0);
  }
  glNormal3f(0.0,0.0,1.0);
  for(i=0;i<ncirc;i++){
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height);
    glVertex3f(                    0.0,                    0.0, height);
  }
  glEnd();

}


/* ----------------------- drawhexdisk ----------------------------- */

void drawpolydisk(int nsides, float diameter, float height, float *rgbcolor){
  int i;

  float x[33], y[33], xnorm[32], ynorm[32];
  float radius;
  float factor,factor2,pi;

  if(nsides>32)nsides=32;
  if(nsides<3)nsides=3;

  pi = 4.0*atan(1.0);
  factor=2.0*pi/nsides;
  factor2 = factor/2.0;

  for(i=0;i<nsides;i++){
    x[i]=cos(i*factor);
    y[i]=sin(i*factor);
    xnorm[i] = cos(factor2+i*factor);
    ynorm[i] = sin(factor2+i*factor);
  }
  x[nsides] = x[0];
  y[nsides] = y[0];

  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);

  radius = diameter/2.0;

  for(i=0;i<nsides;i++){
    glNormal3f(xnorm[i], ynorm[i], 0.0);
    glVertex3f(radius*x[  i],radius*y[  i],0.0); // 1

    glVertex3f(radius*x[i+1],radius*y[i+1],0.0); // 2

    glVertex3f(radius*x[i+1],radius*y[i+1], height); // 3

    glVertex3f(radius*x[  i],radius*y[  i], height); // 4
  }
  glEnd();

  glBegin(GL_TRIANGLES);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);

  glNormal3f(0.0,0.0,-1.0);
  for(i=0;i<nsides;i++){
    glVertex3f(radius*x[  i],radius*y[  i],0.0);
    glVertex3f(            0.0,            0.0,0.0);
    glVertex3f(radius*x[i+1],radius*y[i+1],0.0);
  }
  glNormal3f(0.0,0.0,1.0);
  for(i=0;i<nsides;i++){
    glVertex3f(radius*x[  i],radius*y[  i], height);
    glVertex3f(radius*x[i+1],radius*y[i+1], height);
    glVertex3f(            0.0,            0.0, height);
  }
  glEnd();
}

/* ----------------------- drawhexdisk ----------------------------- */

void drawhexdisk(float diameter, float height, float *rgbcolor){
  int i;

  float x[7]={0.866,0.0,-0.866,-0.866,0.0 ,0.866,0.866};
  float y[7]={0.5,  1.0, 0.5,  -0.5, -1.0,-0.5,  0.5};
  float xnorm[6]={0.500, -0.500, -1.0,-0.500,  0.500, 1.0};
  float ynorm[6]={0.866,  0.866,  0.0,-0.866, -0.866, 0.0};
  float radius;

  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);

  radius = diameter/2.0;

  for(i=0;i<6;i++){
    glNormal3f(xnorm[i], ynorm[i], 0.0);
    glVertex3f(radius*x[  i],radius*y[  i],0.0); // 1

    glVertex3f(radius*x[i+1],radius*y[i+1],0.0); // 2

    glVertex3f(radius*x[i+1],radius*y[i+1], height); // 3

    glVertex3f(radius*x[  i],radius*y[  i], height); // 4
  }
  glEnd();

  glBegin(GL_TRIANGLES);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);

  glNormal3f(0.0,0.0,-1.0);
  for(i=0;i<6;i++){
    glVertex3f(radius*x[  i],radius*y[  i],0.0);
    glVertex3f(            0.0,            0.0,0.0);
    glVertex3f(radius*x[i+1],radius*y[i+1],0.0);
  }
  glNormal3f(0.0,0.0,1.0);
  for(i=0;i<6;i++){
    glVertex3f(radius*x[  i],radius*y[  i], height);
    glVertex3f(radius*x[i+1],radius*y[i+1], height);
    glVertex3f(            0.0,            0.0, height);
  }
  glEnd();
}


/* ----------------------- drawnotchplate ----------------------------- */

void drawnotchplate(float diameter, float height, float notchheight, float direction, float *rgbcolor){
  int i;
  float diameter2;

  diameter2 = diameter + notchheight;

  if(ncirc==0)initcircle(12);
  if(cullfaces==1)glDisable(GL_CULL_FACE);


  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);

  for(i=0;i<ncirc;i++){
    float xmid, ymid;

    xmid = (xcirc[i]+xcirc[i+1])/2.0;
    ymid = (ycirc[i]+ycirc[i+1])/2.0;

    glNormal3f(xcirc[i],ycirc[i],0.0);
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1

    glNormal3f(xcirc[i+1],ycirc[i+1],0.0);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,0.0); // 2

    glNormal3f(xcirc[i+1],ycirc[i+1],0.0);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height); // 3

    glNormal3f(xcirc[i],ycirc[i],0.0);
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height); // 4

    // draw notch

    if(direction<0.0){
      glNormal3f(xcirc[i],ycirc[i],0.0);
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1

      glNormal3f(xcirc[i],ycirc[i],0.0);
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, 0.0-notchheight); // 4

      glNormal3f(xmid,ymid,0.0);
      glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0, 0.0-notchheight); // 3

      glNormal3f(xmid,ymid,0.0);
      glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0,0.0); // 2
    }
    else{
      // top plate
      glNormal3f(0.0,0.0,1.0);
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,height); // 1t

      glNormal3f(0.0,0.0,1.0);
      glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, height); // 4t

      glNormal3f(0.0,0.0,1.0);
      glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, height); // 3t

      glNormal3f(0.0,0.0,1.0);
      glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0,height); // 2t

      // bottom plate

      glNormal3f(0.0,0.0,-1.0);
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1b

      glNormal3f(0.0,0.0,-1.0);
      glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0,0.0); // 2b

      glNormal3f(0.0,0.0,-1.0);
      glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, 0.0); // 3b

      glNormal3f(0.0,0.0,-1.0);
      glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, 0.0); // 4b

      // front plate

      glNormal3f(xcirc[i],ycirc[i],0.0);
      glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, height); // 4t-1

      glNormal3f(xcirc[i],ycirc[i],0.0);
      glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, 0.0); // 4b-4

      glNormal3f(xmid,ymid,0.0);
      glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, 0.0); // 3b-3

      glNormal3f(xmid,ymid,0.0);
      glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, height); // 3t-2

      // left plate

      glNormal3f(-ycirc[i],xcirc[i],0.0);
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,height); // 1t

      glNormal3f(-ycirc[i],xcirc[i],0.0);
      glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, height); // 4t

      glNormal3f(-ycirc[i],xcirc[i],0.0);
      glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, 0.0); // 4b

      glNormal3f(-ycirc[i],xcirc[i],0.0);
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1b

      // right plate

      glNormal3f(ymid,-xmid,0.0);
      glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0,height); // 1t

      glNormal3f(ymid,-xmid,0.0);
      glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0,0.0); // 1b

      glNormal3f(ymid,-xmid,0.0);
      glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, 0.0); // 4b

      glNormal3f(ymid,-xmid,0.0);
      glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, height); // 4t
    }

  }
  glEnd();

  glBegin(GL_TRIANGLES);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);

  glNormal3f(0.0,0.0,-1.0);
  for(i=0;i<ncirc;i++){
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0);
    glVertex3f(                    0.0,                    0.0,0.0);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,0.0);
  }
  glNormal3f(0.0,0.0,1.0);
  for(i=0;i<ncirc;i++){
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height);
    glVertex3f(                    0.0,                    0.0, height);
  }
  glEnd();
  if(cullfaces==1)glEnable(GL_CULL_FACE);
}

/* ----------------------- drawcone ----------------------------- */

void drawcone(float d1, float height, float *rgbcolor){
  int i;
  float dz;

  if(ncirc==0)initcircle(12);
  if(height<=0.0)height=0.0001;

  dz = d1/height;

  glBegin(GL_TRIANGLES);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);

  for(i=0;i<ncirc;i++){
    glNormal3f(xcirc[i],ycirc[i],dz);
    glVertex3f(d1*xcirc[  i]/2.0,d1*ycirc[  i]/2.0,0.0); // 1

    glNormal3f(xcirc[i+1],ycirc[i+1],dz);
    glVertex3f(d1*xcirc[i+1]/2.0,d1*ycirc[i+1]/2.0,0.0); // 2

    glNormal3f(xcirc[i],ycirc[i],dz);
    glVertex3f(0.0,0.0, height); // 3
  }
  glNormal3f(0.0,0.0,-1.0);
  for(i=0;i<ncirc;i++){
    glVertex3f(d1*xcirc[  i]/2.0,d1*ycirc[  i]/2.0,0.0);
    glVertex3f(                    0.0,                    0.0,0.0);
    glVertex3f(d1*xcirc[i+1]/2.0,d1*ycirc[i+1]/2.0,0.0);
  }
  glEnd();
}


/* ----------------------- drawtrunccone ----------------------------- */

void drawtrunccone(float d1, float d2, float height, float *rgbcolor){
  int i;
  float dz;

  if(ncirc==0)initcircle(12);
  if(height<=0.0)height=0.0001;

  dz = -(d2-d1)/height;

  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);

  for(i=0;i<ncirc;i++){
    glNormal3f(xcirc[i],ycirc[i],dz);
    glVertex3f(d1*xcirc[  i]/2.0,d1*ycirc[  i]/2.0,0.0); // 1

    glNormal3f(xcirc[i+1],ycirc[i+1],dz);
    glVertex3f(d1*xcirc[i+1]/2.0,d1*ycirc[i+1]/2.0,0.0); // 2

    glNormal3f(xcirc[i+1],ycirc[i+1],dz);
    glVertex3f(d2*xcirc[i+1]/2.0,d2*ycirc[i+1]/2.0, height); // 3

    glNormal3f(xcirc[i],ycirc[i],dz);
    glVertex3f(d2*xcirc[  i]/2.0,d2*ycirc[  i]/2.0, height); // 4
  }
  glEnd();

  glBegin(GL_TRIANGLES);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);

  glNormal3f(0.0,0.0,-1.0);
  for(i=0;i<ncirc;i++){
    glVertex3f(d1*xcirc[  i]/2.0,d1*ycirc[  i]/2.0,0.0);
    glVertex3f(                    0.0,                    0.0,0.0);
    glVertex3f(d1*xcirc[i+1]/2.0,d1*ycirc[i+1]/2.0,0.0);
  }
  glNormal3f(0.0,0.0,1.0);
  for(i=0;i<ncirc;i++){
    glVertex3f(d2*xcirc[  i]/2.0,d2*ycirc[  i]/2.0, height);
    glVertex3f(d2*xcirc[i+1]/2.0,d2*ycirc[i+1]/2.0, height);
    glVertex3f(                    0.0,                    0.0, height);
  }
  glEnd();

}

/* ----------------------- get_SMVOBJECT_type ----------------------------- */

sv_object *get_SVOBJECT_type(char *label){
  int i;
  sv_object *objecti;

  label = trim_front(label);
  for(i=0;i<ndevice_defs;i++){
    objecti = device_defs[i];
    if(strcmp(label,objecti->label)==0){
      objecti->used=1;
      return objecti;
    }
  }
  return NULL;
}

/* ----------------------- initcircle ----------------------------- */

void initcircle(unsigned int npoints){
  float drad, pi;
  int i;

  if(ncirc!=0)freecircle();
  if(npoints<2)return;
  ncirc=npoints;
  NewMemory( (void **)&xcirc,(ncirc+1)*sizeof(float));
  NewMemory( (void **)&ycirc,(ncirc+1)*sizeof(float));
  pi=4.0*atan(1.0);
  drad=2.0*pi/(float)ncirc;
  for(i=0;i<ncirc;i++){
    xcirc[i] = cos(i*drad);
    ycirc[i] = sin(i*drad);
  }
  xcirc[ncirc]=xcirc[0];
  ycirc[ncirc]=ycirc[0];

}

/* ----------------------- freecircle ----------------------------- */

void freecircle(void){
  FREEMEMORY(xcirc);
  FREEMEMORY(ycirc);
  ncirc=0;
}

/* ----------------------- init_SMVOBJECT1 ----------------------------- */

sv_object *init_SVOBJECT1(char *label, char *commands, int visible){
  sv_object *object;
  sv_object_frame *framei;

  NewMemory( (void **)&object,sizeof(sv_object));
  object->used=0;
  object->visible=visible;
  strcpy(object->label,label);
  object->nframes=1;
  object->obj_frames=NULL;
  NewMemory((void **)&framei,sizeof(sv_object_frame));
  NewMemory((void **)&object->obj_frames,object->nframes*sizeof(sv_object_frame *));

  object->obj_frames[0]=framei;
  getargsops(commands,&framei->args,&framei->nargs,&framei->ops,&framei->nops);
  framei->display_list_ID=-1;
  framei->error=0;

  return object;
}

/* ----------------------- make_error_frame ----------------------------- */

void make_error_frame(void){
  char buffer[256];

  error_frame=NULL;
  NewMemory((void **)&error_frame,sizeof(sv_object_frame));
  strcpy(buffer,"1.0 0.0 0.0 setcolor 0.1 drawsphere");
  getargsops(buffer,&error_frame->args,&error_frame->nargs,&error_frame->ops,&error_frame->nops);
  error_frame->display_list_ID=-1;
}

/* ----------------------- init_SMVOBJECT2 ----------------------------- */

sv_object *init_SVOBJECT2(char *label, char *commandsoff, char *commandson, int visible){
  sv_object *object;
  int i;

  NewMemory( (void **)&object,sizeof(sv_object));
  object->used=0;
  object->visible=visible;
  strcpy(object->label,label);
  object->nframes=2;
  object->obj_frames=NULL;
  NewMemory((void **)&object->obj_frames,object->nframes*sizeof(sv_object_frame *));

  for(i=0;i<object->nframes;i++){
    sv_object_frame *framei;

    if(i==0){
      NewMemory((void **)&framei,sizeof(sv_object_frame));
      object->obj_frames[0]=framei;
      framei->error=0;
      getargsops(commandsoff,&framei->args,&framei->nargs,&framei->ops,&framei->nops);
      framei->display_list_ID=-1;
    }
    else{
      NewMemory((void **)&framei,sizeof(sv_object_frame));
      object->obj_frames[1]=framei;
      getargsops(commandson,&framei->args,&framei->nargs,&framei->ops,&framei->nops);
      framei->error=0;
      framei->display_list_ID=-1;
    }
    

  }

  return object;
}


/* ----------------------- getargsops ----------------------------- */

void getargsops(char *buffer,float **args,int *nargs, int **ops, int *nops){
  char *token;
  char buffer2[256];
  char buffer_save[256];
  char *bufptr;
  float *local_args;
  int *local_ops;
  int iop;
  int local_nops, local_nargs;
  int numargs;
//  int error_code;

  strcpy(buffer_save,buffer);

  strcpy(buffer2,buffer);
  bufptr=buffer2;

  local_nargs=0;
  local_nops=0;
  token=strtok(buffer," ");
  while(token!=NULL){
    char c;

    c=token[0];
    if(c>='a'&&c<='z'||c>='A'&&c<='Z'){
      local_nops++;
    }
    else{
      local_nargs++;
    }
    token=strtok(NULL," ");
  }
  local_args=NULL;
  if(local_nargs>0)NewMemory((void  **)&local_args,local_nargs*sizeof(float));
  *args=local_args;
  *nargs=local_nargs;

  local_ops=NULL;
  if(local_nops>0)NewMemory((void **)&local_ops,local_nops*sizeof(int));
  *ops=local_ops;
  *nops=local_nops;
  if(local_nops==0)return;

  token=strtok(bufptr," ");
  numargs = 0;

  while(token!=NULL){
    char c;

    c=token[0];
    if(c>='a'&&c<='z'||c>='A'&&c<='Z'){
      if(strcmp(token,"translate")==0){
        iop=SV_TRANSLATE;
        reporterror(buffer_save,token,numargs,SV_TRANSLATE_NUMARGS);
      }
      else if(strcmp(token,"rotatex")==0){
        iop=SV_ROTATEX;
        reporterror(buffer_save,token,numargs,SV_ROTATEX_NUMARGS);
      }
      else if(strcmp(token,"rotatey")==0){
        iop=SV_ROTATEY;
        reporterror(buffer_save,token,numargs,SV_ROTATEY_NUMARGS);
      }
      else if(strcmp(token,"rotatez")==0){
        iop=SV_ROTATEZ;
        reporterror(buffer_save,token,numargs,SV_ROTATEZ_NUMARGS);
      }
      else if(strcmp(token,"scalexyz")==0){
        iop=SV_SCALEXYZ;
        reporterror(buffer_save,token,numargs,SV_SCALEXYZ_NUMARGS);
      }
      else if(strcmp(token,"scale")==0&&strcmp(token,"scalexyz")!=0){
        iop=SV_SCALE;
        reporterror(buffer_save,token,numargs,SV_SCALE_NUMARGS);
      }
      else if(strcmp(token,"drawcube")==0){
        iop=SV_DRAWCUBE;
        reporterror(buffer_save,token,numargs,SV_DRAWCUBE_NUMARGS);
      }
      else if(strcmp(token,"drawdisk")==0){
        iop=SV_DRAWDISK;
        reporterror(buffer_save,token,numargs,SV_DRAWDISK_NUMARGS);
      }
      else if(strcmp(token,"drawhexdisk")==0){
        iop=SV_DRAWHEXDISK;
        reporterror(buffer_save,token,numargs,SV_DRAWHEXDISK_NUMARGS);
      }
      else if(strcmp(token,"drawpolydisk")==0){
        iop=SV_DRAWPOLYDISK;
        reporterror(buffer_save,token,numargs,SV_DRAWPOLYDISK_NUMARGS);
      }
      else if(strcmp(token,"drawring")==0){
        iop=SV_DRAWRING;
        reporterror(buffer_save,token,numargs,SV_DRAWRING_NUMARGS);
      }
      else if(strcmp(token,"drawnotchplate")==0){
        iop=SV_DRAWNOTCHPLATE;
        reporterror(buffer_save,token,numargs,SV_DRAWNOTCHPLATE_NUMARGS);
      }
      else if(strcmp(token,"drawtrunccone")==0){
        iop=SV_DRAWTRUNCCONE;
        reporterror(buffer_save,token,numargs,SV_DRAWTRUNCCONE_NUMARGS);
      }
      else if(strcmp(token,"drawcone")==0){
        iop=SV_DRAWCONE;
        reporterror(buffer_save,token,numargs,SV_DRAWCONE_NUMARGS);
      }
      else if(strcmp(token,"drawsphere")==0){
        iop=SV_DRAWSPHERE;
        reporterror(buffer_save,token,numargs,SV_DRAWSPHERE_NUMARGS);
      }
      else if(strcmp(token,"drawline")==0){
        iop=SV_DRAWLINE;
        reporterror(buffer_save,token,numargs,SV_DRAWLINE_NUMARGS);
      }
      else if(strcmp(token,"drawpoint")==0){
        iop=SV_DRAWPOINT;
        reporterror(buffer_save,token,numargs,SV_DRAWPOINT_NUMARGS);
      }
      else if(strcmp(token,"drawcircle")==0){
        iop=SV_DRAWCIRCLE;
        reporterror(buffer_save,token,numargs,SV_DRAWCIRCLE_NUMARGS);
      }
      else if(strcmp(token,"setcolor")==0){
        iop=SV_SETCOLOR;
        reporterror(buffer_save,token,numargs,SV_SETCOLOR_NUMARGS);
      }
      else if(strcmp(token,"setlinewidth")==0){
        iop=SV_SETLINEWIDTH;
        reporterror(buffer_save,token,numargs,SV_SETLINEWIDTH_NUMARGS);
      }
      else if(strcmp(token,"setpointsize")==0){
        iop=SV_SETPOINTSIZE;
        reporterror(buffer_save,token,numargs,SV_SETPOINTSIZE_NUMARGS);
      }
      else if(strcmp(token,"setbw")==0){
        iop=SV_SETBW;
        reporterror(buffer_save,token,numargs,SV_SETBW_NUMARGS);
      }
      else if(strcmp(token,"push")==0){
        iop=SV_PUSH;
        reporterror(buffer_save,token,numargs,SV_PUSH_NUMARGS);
      }
      else if(strcmp(token,"pop")==0){
        iop=SV_POP;
        reporterror(buffer_save,token,numargs,SV_POP_NUMARGS);
      }
#ifdef pp_AVATAR
      else if(strcmp(token,"getuservals")==0){
        iop=SV_GETUSERVALS;
        reporterror(buffer_save,token,numargs,SV_GETUSERVALS_NUMARGS);
      }
#endif
      else{
        iop=SV_ERR;
      }
      *local_ops++=iop;
      numargs=0;
    }
    else{
      sscanf(token,"%f",local_args);
      local_args++;
      numargs++;
    }
    token=strtok(NULL," ");
  }
}


/* ----------------------- read_device_defs ----------------------------- */

int read_device_defs(char *file){
  FILE *stream;
  char buffer[256], *trim_buffer;
  sv_object *temp_object, *prev_object, *next_object, *current_object;
  sv_object_frame *current_frame;
  int firstdef;
  float *arglist;
  int *oplist, nargs, nops;
  sv_object *object_start, *objecti;
  size_t lenbuffer;
  int ndevices=0;

  stream=fopen(file,"r");
  if(stream==NULL)return 0;
  printf("Reading device definitions from: %s\n",file);

  firstdef=-1;
  while(!feof(stream)){
    CheckMemory;
    if(fgets(buffer,255,stream)==NULL)break;
    remove_comment(buffer);
    trim(buffer);
    trim_buffer=trim_front(buffer);
    lenbuffer=strlen(buffer);
    if(lenbuffer<1)continue;


#ifdef pp_AVATAR
    if(match(buffer,"DEVICEDEF",9) == 1||
       match(buffer,"AVATARDEF",9) == 1
      ){
        int is_avatar=0;
#else
    if(match(buffer,"DEVICEDEF",9) == 1){
#endif
      char *label;

      sv_object_frame *first_frame, *last_frame;

#ifdef pp_AVATAR
      if(match(buffer,"AVATARDEF",9) == 1){
        is_avatar=1;
      }  
#endif
      ndevices++;
      if(fgets(buffer,255,stream)==NULL)break;
      remove_comment(buffer);
      trim(buffer);
      label = trim_front(buffer);
      temp_object=get_object(label);
      if(temp_object!=NULL){
        free_object(temp_object);
      }
  
      NewMemory((void **)&current_object,sizeof(sv_object));

      strcpy(current_object->label,label);
      prev_object = device_def_last.prev;
      next_object = &device_def_last;

      prev_object->next=current_object;
      next_object->prev=current_object;

      current_object->next=next_object;
      current_object->prev=prev_object;
      current_object->visible=1;

      current_object->nframes=0;

      first_frame = &current_object->first_frame;
      last_frame = &current_object->last_frame;
#ifdef pp_AVATAR
      current_object->type=is_avatar;
#endif

      first_frame->next=last_frame;
      first_frame->prev=NULL;
      last_frame->prev=first_frame;
      last_frame->next=NULL;

      firstdef=1;
      continue;
    }
    if(match(trim_buffer,"NEWFRAME",8) == 1||firstdef==1){
      sv_object_frame *prev_frame,*next_frame;

      if(firstdef==-1)continue;
      NewMemory((void **)&current_frame,sizeof(sv_object_frame));

      next_frame=&current_object->last_frame;
      prev_frame=next_frame->prev;

      next_frame->prev=current_frame;
      prev_frame->next=current_frame;

      current_frame->next=next_frame;
      current_frame->prev=prev_frame;

      current_frame->args=NULL;
      current_frame->ops=NULL;
      current_frame->nargs=0;
      current_frame->nops=0;
      current_frame->display_list_ID=-1;

      current_object->nframes++;

      firstdef=0;
      if(match(trim_buffer,"NEWFRAME",8)==1)continue;
    }
    getargsops(buffer,&arglist, &nargs, &oplist, &nops);
    if(nargs>0){
      if(current_frame->nargs==0){
        current_frame->args=arglist;
        current_frame->nargs=nargs;
      }
      else{
        int i, ncurrent;

        ncurrent = current_frame->nargs;
        ResizeMemory((void **)&current_frame->args,(ncurrent+nargs)*sizeof(float));
        for(i=0;i<nargs;i++){
          current_frame->args[ncurrent+i] = arglist[i];
        }
        current_frame->nargs=ncurrent+nargs;
        FREEMEMORY(arglist);
      }
    }
    if(nops>0){
      if(current_frame->nops==0){
        current_frame->ops=oplist;
        current_frame->nops=nops;
      }
      else{
        int i, ncurrent;

        ncurrent = current_frame->nops;
        ResizeMemory((void **)&current_frame->ops,(ncurrent+nops)*sizeof(int));
        for(i=0;i<nops;i++){
          current_frame->ops[ncurrent+i] = oplist[i];
        }
        current_frame->nops=ncurrent+nops;
        FREEMEMORY(oplist);
      }
    }

  }

  fclose(stream);

  object_start = device_def_first.next;
  objecti = object_start;
  ndevice_defs=0;
  for(;objecti->next!=NULL;){
    CheckMemory;
    ndevice_defs++;
    objecti->obj_frames=NULL;
    if(objecti->nframes>0){
      NewMemory((void **)&objecti->obj_frames,objecti->nframes*sizeof(sv_object_frame *));
    }
    objecti=objecti->next;
  }
  FREEMEMORY(device_defs);
  if(ndevice_defs>0){
    int i,j;

    NewMemory((void **)&device_defs,ndevice_defs*sizeof(sv_object *));

    object_start = device_def_first.next;
    objecti = object_start;
    i=0;
    for(;objecti->next!=NULL;){
      sv_object_frame *frame_start, *framei;

      CheckMemory;
      device_defs[i]=objecti;
      i++;
      frame_start = objecti->first_frame.next;
      framei = frame_start;
      j=0;
      for(;framei->next!=NULL;){
        int iop, npushpop=0;

        CheckMemory;
        objecti->obj_frames[j]=framei;
        framei->error=0;
        for(iop=0;iop<framei->nops;iop++){
          int op;

          op = framei->ops[iop];
          if(op==SV_PUSH){
            npushpop++;
          }
          else if(op==SV_POP){
            npushpop--;
            if(npushpop<0){
              npushpop=0;
              framei->ops[iop]=SV_NO_OP;
            }
          }
        }
        if(npushpop>0){
          framei->error=1;
        }
        framei=framei->next;
        j++;
      }
      objecti=objecti->next;
    }
  }
  return ndevices;
}

/* ----------------------- reporterror ----------------------------- */

void reporterror(char *buffer, char *token, int numargs_found, int numargs_expected){
  if(numargs_found==numargs_expected)return;
  printf("*** Error:  %i arguments were found for the token \"%s\" while parsing: \"%s\".  %i arguments were expected.\n",numargs_found,token,buffer,numargs_expected);
}

/* ----------------------- get_device_label ----------------------------- */

char *get_device_label(char *buffer){
  char *label_present;

  label_present=strstr(buffer,"#");
  if(label_present==NULL)return NULL;
  if(strlen(label_present)<=1){
    label_present[0]=0;
    return NULL;
  }
  label_present[0]=0;
  label_present++;
  label_present=trim_front(label_present);
  trim(label_present);
  if(strlen(label_present)==0)return NULL;
  return label_present;
}

/* ----------------------- get_object ----------------------------- */

sv_object *get_object(char *label){
  sv_object *objecti,*object_start;

  object_start = device_def_first.next;
  objecti = object_start;
  for(;objecti->next!=NULL;objecti=objecti->next){
    if(strcmp(objecti->label,label)==0)return objecti;
  }
  return NULL;
}

/* ----------------------- freeall_objects ----------------------------- */

void freeall_objects(void){
  sv_object *object;

  for(;;){
    object = device_def_last.prev;
    if(object->prev==NULL)break;
    free_object(object);
  }
}

/* ----------------------- free_object ----------------------------- */

void free_object(sv_object *object){
  sv_object *prev, *next;
  sv_object_frame *framei, *frame_start;

  prev = object->prev;
  next = object->next;

  prev->next=next;
  next->prev=prev;

  frame_start = &object->first_frame;
  framei = frame_start->next;
  for(;framei->next!=NULL;){
    sv_object_frame *next_frame;

    next_frame=framei->next;
    if(framei->nargs>0)FREEMEMORY(framei->args);
    if(framei->nops>0)FREEMEMORY(framei->ops);
    FREEMEMORY(framei);
    framei=next_frame;
  }
  FREEMEMORY(object);
}

/* ----------------------- remove_comment ----------------------------- */

void remove_comment(char *buffer){
  char *comment;
// test
  comment = strstr(buffer,"//");
  if(comment!=NULL)comment[0]=0;
  return;
}

/* ----------------------- init_device_defs ----------------------------- */

void init_device_defs(void){
#ifndef pp_AVATAR
  if(ndeviceinfo>0){
#endif
    char com_buffer[1024];
    char com_buffer2[1024];


    {
      char objectfile[1024];

      svofile_exists = 0;

      if(smvprogdir!=NULL){
        strcpy(objectfile,smvprogdir);
        strcat(objectfile,"devices.svo");
        read_device_defs(objectfile);
      }

      strcpy(objectfile,"devices.svo");
      read_device_defs(objectfile);

      strcpy(objectfile,fdsprefix);
      strcat(objectfile,".svo");
      read_device_defs(objectfile);

#ifdef pp_AVATAR
      init_avatar();
#endif
    }

    if(isZoneFireModel==1){
      strcpy(com_buffer,"1.0 1.0 0.0 setcolor 0.02 0.05 drawdisk");
      device_defs_backup[0] = init_SVOBJECT1("target", com_buffer,1);
    }
    else{
      strcpy(com_buffer,"1.0 1.0 0.0 setcolor 0.038 drawcube");
      device_defs_backup[0] = init_SVOBJECT1("sensor", com_buffer,1);
    }
    strcpy(com_buffer, "0.0 1.0 0.0 setcolor 0.038 drawcube");
    strcpy(com_buffer2,"1.0 0.0 0.0 setcolor 0.038 drawcube");
    device_defs_backup[1] = init_SVOBJECT2("heat_detector", com_buffer, com_buffer2,1);

    strcpy(com_buffer, "0.0 1.0 0.0 setcolor 0.038 drawcube");
    strcpy(com_buffer2,"1.0 0.0 0.0 setcolor 0.038 drawcube");
    device_defs_backup[2] = init_SVOBJECT2("sprinkler_upright", com_buffer, com_buffer2,1);

    strcpy(com_buffer, "0.5 0.5 0.5 setcolor 0.2 0.05 drawdisk");
    strcpy(com_buffer2,"1.0 0.0 0.0 setcolor 0.2 0.05 drawdisk");
    device_defs_backup[3] = init_SVOBJECT2("smoke_detector", com_buffer, com_buffer2,1);

    make_error_frame();
    
    if(ndevice_defs==0){

      ndevice_defs=4;
      FREEMEMORY(device_defs);
      NewMemory((void **)&device_defs,4*sizeof(sv_object *));
      device_defs[0] = device_defs_backup[0];
      device_defs[1] = device_defs_backup[1];
      device_defs[2] = device_defs_backup[2];
      device_defs[3] = device_defs_backup[3];
    }
    //if(isZoneFireModel==1){
    //  int ii;
    //  for(ii=0;ii<ndevice_defs;ii++){
    //   if(strcmp(device_defs[ii]->label,"sensor")==0)strcpy(device_defs[ii]->label,"target");
    //  }
    //}
#ifndef pp_AVATAR
  }
  else{
    ndevice_defs=0;
  }
#endif
}
#ifdef pp_AVATAR
  void init_avatar(void){
    int iavatar_types;
    sv_object *objecti,*object_start;
    char com_buffer[1024];
    
    object_start = device_def_first.next;
    navatar_types=2;
    for(objecti = object_start;objecti->next!=NULL;objecti=objecti->next){
      if(objecti->type==1)navatar_types++;
    }
    NewMemory((void **)&avatar_types,navatar_types*sizeof(sv_object *));

//    strcpy(com_buffer,"1.0 0.0 0.0 setcolor 0.0 0.0 0.0 1.0 0.0 0.0 drawline 0.0 0.0 0.0 0.0 0.0 1.0 drawline");
    strcpy(com_buffer,"1.0 0.0 0.0 setcolor 0.03 0.1 drawdisk 0.0 0.0 1.0 setcolor 90.0 rotatey 0.03 0.2 drawdisk");
    avatar_defs_backup[0] = init_SVOBJECT1("Avatar_1", com_buffer,1);
    avatar_defs_backup[0]->type=1;

//    strcpy(com_buffer,"0.0 0.0 1.0 setcolor 0.1 drawcircle");
    strcpy(com_buffer,"1.0 1.0 0.0 setcolor 0.02 0.05 drawdisk");
    avatar_defs_backup[1] = init_SVOBJECT1("Avatar_2", com_buffer,1);
    avatar_defs_backup[1]->type=1;

    avatar_types[0]=avatar_defs_backup[0];
    avatar_types[1]=avatar_defs_backup[1];

    iavatar_types=2;
    for(objecti = object_start;objecti->next!=NULL;objecti=objecti->next){
      if(objecti->type==0)continue;
      avatar_types[iavatar_types++]=objecti;
    }
    iavatar_types=0;
  }
#endif
