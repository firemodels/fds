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

#define CIRCLE_SEGS 12

#define SV_TRANSLATE  100
#define SV_ROTATEX    101
#define SV_ROTATEY    102
#define SV_ROTATEZ    103
#define SV_SCALEXYZ   104
#define SV_SCALE      105
#define SV_OFFSETX 108
#define SV_OFFSETY 109
#define SV_OFFSETZ 110
#define SV_MULTIADDT 112
#define SV_CLIP 113
#define SV_MIRRORCLIP 114
#define SV_PERIODICCLIP 115

#define SV_TRANSLATE_NUMARGS  3
#define SV_ROTATEX_NUMARGS    1
#define SV_ROTATEY_NUMARGS    1
#define SV_ROTATEZ_NUMARGS    1
#define SV_SCALEXYZ_NUMARGS   3
#define SV_SCALE_NUMARGS      1
#define SV_OFFSETX_NUMARGS 1
#define SV_OFFSETY_NUMARGS 1
#define SV_OFFSETZ_NUMARGS 1
#define SV_MULTIADDT_NUMARGS 3
#define SV_CLIP_NUMARGS 4
#define SV_MIRRORCLIP_NUMARGS 4
#define SV_PERIODICCLIP_NUMARGS 4

#define SV_TRANSLATE_NUMOUTARGS  0
#define SV_ROTATEX_NUMOUTARGS    0
#define SV_ROTATEY_NUMOUTARGS    0
#define SV_ROTATEZ_NUMOUTARGS    0
#define SV_SCALEXYZ_NUMOUTARGS   0
#define SV_SCALE_NUMOUTARGS      0
#define SV_OFFSETX_NUMOUTARGS 0
#define SV_OFFSETY_NUMOUTARGS 0
#define SV_OFFSETZ_NUMOUTARGS 0
#define SV_MULTIADDT_NUMOUTARGS 1
#define SV_CLIP_NUMOUTARGS 1
#define SV_MIRRORCLIP_NUMOUTARGS 1
#define SV_PERIODICCLIP_NUMOUTARGS 1

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
#define SV_DRAWARC       212
#define SV_DRAWCDISK     213
#define SV_DRAWTSPHERE   214
#define SV_DRAWARCDISK   215

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
#define SV_DRAWARC_NUMARGS       2
#define SV_DRAWCDISK_NUMARGS     2
#define SV_DRAWTSPHERE_NUMARGS   2
#define SV_DRAWARCDISK_NUMARGS   3

#define SV_DRAWCUBE_NUMOUTARGS      0
#define SV_DRAWSPHERE_NUMOUTARGS    0
#define SV_DRAWDISK_NUMOUTARGS      0
#define SV_DRAWLINE_NUMOUTARGS      0
#define SV_DRAWCIRCLE_NUMOUTARGS    0
#define SV_DRAWTRUNCCONE_NUMOUTARGS 0
#define SV_DRAWNOTCHPLATE_NUMOUTARGS 0
#define SV_DRAWRING_NUMOUTARGS      0
#define SV_DRAWCONE_NUMOUTARGS      0
#define SV_DRAWHEXDISK_NUMOUTARGS   0
#define SV_DRAWPOLYDISK_NUMOUTARGS   0
#define SV_DRAWPOINT_NUMOUTARGS     0
#define SV_DRAWARC_NUMOUTARGS       0
#define SV_DRAWCDISK_NUMOUTARGS     0
#define SV_DRAWTSPHERE_NUMOUTARGS   0
#define SV_DRAWARCDISK_NUMOUTARGS   0

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

#define SV_PUSH_NUMOUTARGS       0
#define SV_POP_NUMOUTARGS        0
#define SV_SETCOLOR_NUMOUTARGS   0
#define SV_SETBW_NUMOUTARGS      0
#define SV_SETLINEWIDTH_NUMOUTARGS 0
#define SV_SETPOINTSIZE_NUMOUTARGS 0

#define SV_ERR -1

#define NLAT device_sphere_segments
#define NLONG (2*device_sphere_segments)

#define TOKEN_FLOAT 0
#define TOKEN_COMMAND 1
#define TOKEN_GETVAL 2

char *parse_device_frame(char *buffer, FILE *stream, int *eof, sv_object_frame *frame);
void reporterror(char *buffer, char *token, int numargs_found, int numargs_expected);
float get_point2box_dist(float boxmin[3], float boxmax[3], float p1[3], float p2[3]);

void drawcone(float d1, float height, unsigned char *rgbcolor);
void drawtrunccone(float d1, float d2, float height, unsigned char *rgbcolor);
void drawline(float *xyz1, float *xyz2, unsigned char *rgbcolor);
void drawarc(float angle, float diameter, unsigned char *rgbcolor);
void drawcircle(float diameter, unsigned char *rgbcolor);
void drawpoint(unsigned char *rgbcolor);
void drawsphere(float diameter, unsigned char *rgbcolor);
void drawtsphere(float texture_index, float diameter, unsigned char *rgbcolor);
void drawcube(float size, unsigned char *rgbcolor);
void drawcdisk(float diameter, float height, unsigned char *rgbcolor);
void drawdisk(float diameter, float height, unsigned char *rgbcolor);
void drawarcdisk(float angle, float diameter, float height, unsigned char *rgbcolor);
void drawhexdisk(float diameter, float height, unsigned char *rgbcolor);
void drawpolydisk(int nsides, float diameter, float height, unsigned char *rgbcolor);
void drawring(float d_inner, float d_outer, float height, unsigned char *rgbcolor);
void drawnotchplate(float diameter, float height, float notchheight, float direction, unsigned char *rgbcolor);
void draw_SVOBJECT(sv_object *object, int iframe, float *valstack, int nvalstack);
sv_object *get_object(char *label);
void free_object(sv_object *object);
void remove_comment(char *buffer);
void freecircle(void);
void initcircle(unsigned int npoints);

static float *xcirc=NULL, *ycirc=NULL;
static int ncirc;
static float *cos_long=NULL, *sin_long=NULL, *cos_lat=NULL, *sin_lat=NULL;
static int nlat, nlong;
static float specular[4]={0.4,0.4,0.4,1.0};
unsigned char *rgbimage=NULL;
int rgbsize=0;

/* ------------------ get_mesh ------------------------ */

 mesh *get_mesh(float xyz[3]){
  int i;
  mesh *meshi;
  int ibar, jbar, kbar;
  float xmin, ymin, zmin;
  float xmax, ymax, zmax;

  for(i=0;i<nmeshes;i++){
    meshi = meshinfo + i;

    ibar=meshi->ibar;
    jbar=meshi->jbar;
    kbar=meshi->kbar;
    xmin=meshi->xplt_orig[0];
    ymin=meshi->yplt_orig[0];
    zmin=meshi->zplt_orig[0];
    xmax=meshi->xplt_orig[ibar];
    ymax=meshi->yplt_orig[jbar];
    zmax=meshi->zplt_orig[kbar];

    if(xmin<=xyz[0]&&xyz[0]<=xmax&&
       ymin<=xyz[1]&&xyz[1]<=ymax&&
       zmin<=xyz[2]&&xyz[2]<=zmax){
         return meshi;
    }
  }
  return NULL;
}

/* ------------------ get_world_eyepos ------------------------ */

void get_world_eyepos(float *mm, float eyepos[3]){
    /*
      ( m0 m4 m8  m12 ) (x)    (0)
      ( m1 m5 m9  m13 ) (y)    (0)
      ( m2 m6 m10 m14 ) (z)  = (0)
      ( m3 m7 m11 m15 ) (1)    (1)

       ( m0 m4  m8 )      (m12)
   Q=  ( m1 m5  m9 )  u = (m13)
       ( m2 m6 m10 )      (m14)
      
      (Q   u) (x)     (0)      
      (v^T 1) (y)   = (1)
       
      m3=m7=m11=0, v^T=0, y=1   Qx+u=0 => x=-Q^Tu
    */

  eyepos[0] = -(mm[0]*mm[12]+mm[1]*mm[13]+ mm[2]*mm[14])/mscale[0];
  eyepos[1] = -(mm[4]*mm[12]+mm[5]*mm[13]+ mm[6]*mm[14])/mscale[1];
  eyepos[2] = -(mm[8]*mm[12]+mm[9]*mm[13]+mm[10]*mm[14])/mscale[2];
  eyepos[0] = xbar0 + xyzmaxdiff*eyepos[0];
  eyepos[1] = ybar0 + xyzmaxdiff*eyepos[1];
  eyepos[2] = zbar0 + xyzmaxdiff*eyepos[2];
}

/* ----------------------- getsmokesensors ----------------------------- */

void getsmokesensors(void){
  int doit, i;
//  int index;
  int width, height;
  unsigned char rgbval[3];

  width = screenWidth;
  height = screenHeight;

  doit=0;
  for(i=0;i<ndeviceinfo;i++){
    device *devicei;
    char *label;

    devicei = deviceinfo + i;
    label = devicei->object->label;
    if(STRCMP(label,"smokesensor")!=0)continue;
    doit=1;
  }
  if(doit==0)return;

  if(rgbimage==NULL||rgbsize!=width*height){
    if(rgbimage!=NULL){
      FREEMEMORY(rgbimage);
    }
    rgbsize=width*height;
    NewMemory( (void **)&rgbimage,3*rgbsize*sizeof(GLubyte));
  }
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(0,0,width,height, GL_RGB, GL_UNSIGNED_BYTE, rgbimage);

  for(i=0;i<ndeviceinfo;i++){
    device *devicei;
    char *label;
    int row, col;
    int index,val;

    devicei = deviceinfo + i;
    label = devicei->object->label;

    
    if(STRCMP(label,"smokesensor")!=0)continue;

    col = devicei->screenijk[0];
    row = devicei->screenijk[1];

    if(col<0||col>width-1||row<0||row>height-1){
      val=-1;
    }
    else{
      index=row*width+col;
      val=rgbimage[3*index];
    }
    devicei->visval=val;
  }

}

/* ----------------------- getdevice_screencoords ----------------------------- */

void getdevice_screencoords(void){
  double mv_setup[16], projection_setup[16];
  GLint viewport_setup[4];
//  double d_ijk[3];
  int i;
  int doit;

  doit=0;
  for(i=0;i<ndeviceinfo;i++){
    device *devicei;
    char *label;

    devicei = deviceinfo + i;
    label = devicei->object->label;
    if(STRCMP(label,"smokesensor")!=0)continue;
    doit=1;
  }
  if(doit==0)return;

  glGetDoublev(GL_MODELVIEW_MATRIX,mv_setup);
  glGetDoublev(GL_PROJECTION_MATRIX,projection_setup);
  glGetIntegerv(GL_VIEWPORT, viewport_setup);
  for(i=0;i<ndeviceinfo;i++){
    float *xyz;
    double d_ijk[3];
    device *devicei;
    int *ijk;
    char *label;
    mesh *device_mesh;

    devicei = deviceinfo + i;
    label = devicei->object->label;
    
    if(STRCMP(label,"smokesensor")!=0)continue;
    xyz = devicei->xyz;
    device_mesh = devicei->device_mesh;
    devicei->eyedist = get_point2box_dist(device_mesh->boxmin,device_mesh->boxmax,xyz,world_eyepos);
    ijk = devicei->screenijk;
    gluProject(xyz[0],xyz[1],xyz[2],mv_setup,projection_setup,viewport_setup,d_ijk,d_ijk+1,d_ijk+2);
    ijk[0] = d_ijk[0];
    ijk[1] = d_ijk[1];
    ijk[2] = d_ijk[2];
  }
}

/* ----------------------- draw_devices_val ----------------------------- */

void draw_devices_val(void){
  device *devicei;
  int i;
  float *xyz, *xyznorm;
  float white[3]={1.0,1.0,1.0};
  float black[3]={0.0,0.0,0.0};
  int doit=0;

  glPushMatrix();
  glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);
  glTranslatef(-xbar0,-ybar0,-zbar0);
  if(active_smokesensors==1&&show_smokesensors!=0){
    getdevice_screencoords();
  }
  for(i=0;i<ndeviceinfo;i++){
    devicei = deviceinfo + i;

    if(devicei->object->visible==0)continue;
    xyz = devicei->xyz;
    xyznorm = devicei->xyznorm;
    if(active_smokesensors==1&&show_smokesensors!=0&&STRCMP(devicei->object->label,"smokesensor")==0){
      char label[256];
      float val;
      int ival;

      if(doit==0){
        glBlendFunc(GL_ONE,GL_ZERO);
        doit=1;
      }
      switch (show_smokesensors){
        case 1:
          sprintf(label,"%i",devicei->visval);
          break;
        case 2:
          val = devicei->visval/255.0;
          sprintf(label,"%.2f",val);
          trimzeros(label);
          break;
        case 3:
        case 4:
          ival = devicei->visval;
          if(ival==255){
            strcpy(label,"Inf");
          }
          else{
            float light_extinct;

            val = ival/255.0;
            light_extinct = -log(val)/devicei->eyedist;
            val = smoke3d_cvis/light_extinct;
            if(val<10.0){
              sprintf(label,"%.1f",val);
            }
            else{
              sprintf(label,"%.0f",val);
            }
            trimzeros(label);
          }
          break;
      }
      if(devicei->visval>128){
        output3Text(black,xyz[0]+0.2*xyznorm[0],xyz[1]+0.2*xyznorm[1],xyz[2]+0.2*xyznorm[2],label);
      }
      else{
        output3Text(white,xyz[0]+0.2*xyznorm[0],xyz[1]+0.2*xyznorm[1],xyz[2]+0.2*xyznorm[2],label);
      }
    }
  }
  if(doit==1){
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  }
  glPopMatrix();
}

  /* ----------------------- draw_devices ----------------------------- */

void draw_devices(void){
  device *devicei;
  int i;
  float *xyz;

  if(select_device==0||show_mode!=SELECT){
    for(i=0;i<ndeviceinfo;i++){
      devicei = deviceinfo + i;
      if(devicei->object->visible==0)continue;
      if(devicei->plane_surface!=NULL){
        int j;

        for(j=0;j<nmeshes;j++){
          drawstaticiso(devicei->plane_surface[j],-1,0,2,0,devicei->line_width);
          drawstaticiso(devicei->plane_surface[j],2,0,2,0,devicei->line_width);
        }
        continue;
      }
    }
  }

  glPushMatrix();
  glPushAttrib(GL_POINT_BIT|GL_LINE_BIT);
  glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);
  glTranslatef(-xbar0,-ybar0,-zbar0);
  for(i=0;i<ndeviceinfo;i++){
    int tagval;
    int save_use_displaylist;

    devicei = deviceinfo + i;

    if(devicei->object->visible==0)continue;
    if(devicei->plane_surface!=NULL){
      continue;
    }
    if(isZoneFireModel==1&&STRCMP(devicei->object->label,"target")==0&&visSensor==0)continue;

    save_use_displaylist=devicei->object->use_displaylist;
    tagval=i+1;
    if(select_device==1&&show_mode==SELECT){

      select_device_color[0]=tagval>>(ngreenbits+nbluebits);
      select_device_color[1]=tagval>>nbluebits;
      select_device_color[2]=tagval&rgbmask[nbluebits-1];
      select_device_color_ptr=select_device_color;
      devicei->object->use_displaylist=0;
    }
    else{
      if(selected_device_tag>0&&select_device==1&&selected_device_tag==tagval){
        select_device_color_ptr=select_device_color;
        select_device_color[0]=255;
        select_device_color[1]=0;
        select_device_color[2]=0;
        devicei->object->use_displaylist=0;
      }
      else{
        select_device_color_ptr=NULL;
      }
    }

    xyz = devicei->xyz;
    glPushMatrix();
    glTranslatef(xyz[0],xyz[1],xyz[2]);

    if(active_smokesensors==1&&show_smokesensors!=0&&STRCMP(devicei->object->label,"smokesensor")==0){
      float *xyznorm;

      xyznorm = devicei->xyznorm;
      xyznorm[0]=world_eyepos[0]-devicei->xyz[0];
      xyznorm[1]=world_eyepos[1]-devicei->xyz[1];
      xyznorm[2]=world_eyepos[2]-devicei->xyz[2];

      get_elevaz(xyznorm,&devicei->angle_elev,&devicei->angle_az);
    }
    if(devicei->angle_az!=0.0){
      glRotatef(devicei->angle_az,0.0,0.0,1.0);
    }
    if(devicei->angle_elev!=0.0){
      glRotatef(-devicei->angle_elev,0.0,1.0,0.0);
    }
    if(sensorrelsize!=1.0){
      glScalef(sensorrelsize,sensorrelsize,sensorrelsize);
    }
    if(devicei->nparams>0){
      int j;

      for(j=0;j<devicei->nparams;j++){
        valstack[j]=devicei->params[j];
      }
    }
    nvalstack=devicei->nparams;
    ntexturestack=0;
    if(devicei->ntextures>0){
      if(devicei->textureinfo!=NULL){
        int j;

        for(j=0;j<devicei->ntextures;j++){
          valstack[nvalstack+0+j]=devicei->params[nvalstack+j];
        }
      }
      ntexturestack=devicei->ntextures;
    }
    if(showtime==1&&itime>=0&&itime<ntimes&&devicei->showstatelist!=NULL){
      int state;

      state=devicei->showstatelist[itime];
      draw_SVOBJECT(devicei->object,state,valstack,nvalstack);
    }
    else{
      draw_SVOBJECT(devicei->object,devicei->state0,valstack,nvalstack);
    }
    devicei->object->use_displaylist=save_use_displaylist;
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
      if(STRCMP(devicei->object->label,"sensor")==0&&visSensor==0)continue;
      if(isZoneFireModel==1&&STRCMP(devicei->object->label,"target")==0&&visSensor==0)continue;
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

void draw_SVOBJECT(sv_object *object, int iframe, float *valstack, int nvalstack){
  sv_object_frame *framei;
  int *op;
  unsigned char *rgbptr;
  unsigned char rgbcolor[4];
  int displaylist_id=0;
  int ii;

  if(iframe>object->nframes-1||iframe<0)iframe=0;
  framei=object->obj_frames[iframe];
  ASSERT(framei->error==0||framei->error==1);
  if(framei->error!=0)framei=error_frame;

  rgbcolor[0]=255;
  rgbcolor[1]=0;
  rgbcolor[2]=0;
  rgbcolor[3]=255;
  rgbptr=rgbcolor;
  glPushMatrix();

  if(valstack!=NULL&&nvalstack>0){
    int i;

    for(i=0;i<nvalstack;i++){
      tokendata *toki;

      toki = framei->tokens + i;
      if(i>framei->ntokens-1)break;
      toki->var=valstack[i];
    }
  }

  if(framei->display_list_ID!=-1&&object->use_displaylist==1){
    if(framei->use_bw==setbw){
      glCallList(framei->display_list_ID);
      glPopMatrix();
      return;
    }
    else{
      framei->use_bw=setbw;
      glDeleteLists(framei->display_list_ID,1);
      framei->display_list_ID=-1;
    }
  }

  if(object->use_displaylist==1){   
    displaylist_id = glGenLists(1);
    if(displaylist_id!=0){
      framei->display_list_ID=displaylist_id;
      glNewList(displaylist_id,GL_COMPILE_AND_EXECUTE);
    }
  }

  if(select_device_color_ptr==NULL){
    glEnable(GL_LIGHTING);

    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,specular);

    glEnable(GL_COLOR_MATERIAL);
  }

  for(ii=0;ii<framei->ncommands;ii++){
    tokendata *toki,*tok1,*tok2,*tok3,*tok4;
#define NARGVAL 6
    float arg[NARGVAL], *argptr;
    int j;

    if(select_device_color_ptr==NULL){
      rgbptr=rgbcolor;
    }
    else{
      rgbptr=select_device_color_ptr;
    }
    toki = framei->command_list[ii];
    for(j=0;j<toki->nvars;j++){
      tokendata *tokj;
      
      tokj = toki - toki->nvars + j;
      arg[j] = *(tokj->varptr);
    }
    if(toki->nvars>0){
      argptr=(toki-1)->varptr;
    }

    switch (toki->command){
    case SV_MULTIADDT:
      {
        float val1, val2, val_result;
        float time_val=0.0;

        val1=arg[0];
        val2=arg[1];

        if(ntimes>0){
          time_val=times[itime];
        }

        val_result=val1*time_val+val2;

        *argptr=val_result;
      }
      break;
    case SV_CLIP:
      {
        int argval, argmin, argmax, stackskip;
        float val, valmin, valmax;

        val=arg[0];

        valmin=arg[1];

        valmax=arg[2];

        if(val<valmin)val=valmin;
        if(val>valmax)val=valmax;
        
        *argptr=val;
      }
      break;
    case SV_MIRRORCLIP:
      {
        float val, valmin, valmax;
        float val2, valmax2;
        float val_result;

        val=arg[0];

        valmin=arg[1];

        valmax=arg[2];

        val2=val-valmin;
        valmax2=valmax-valmin;

        val2=fmod(val2,2.0*valmax2);
        if(val2>valmax2)val2=2.0*valmax2-val2;

        val_result = val2 + valmin;

        *argptr=val_result;
      }
      break;
    case SV_PERIODICCLIP:
      {
        float val, valmin, valmax;
        float val2, valmax2;
        float val_result;

        val=arg[0];

        valmin=arg[1];

        valmax=arg[2];

        val2=val-valmin;
        valmax2=valmax-valmin;

        val2=fmod(val2,valmax2);

        val_result = val2 + valmin;

        *argptr=val_result;
      }
      break;
    case SV_TRANSLATE:
      glTranslatef(arg[0],arg[1],arg[2]);
      break;
    case SV_OFFSETX:
      glTranslatef(arg[0],0.0,0.0);
      break;
    case SV_OFFSETY:
      glTranslatef(0.0,arg[0],0.0);
      break;
    case SV_OFFSETZ:
      glTranslatef(0.0,0.0,arg[0]);
      break;
    case SV_ROTATEX:
      glRotatef(arg[0],1.0,0.0,0.0);
      break;
    case SV_ROTATEY:
      glRotatef(arg[0],0.0,1.0,0.0);
      break;
    case SV_ROTATEZ:
      glRotatef(arg[0],0.0,0.0,1.0);
      break;
    case SV_SCALEXYZ:
      glScalef(arg[0],arg[1],arg[2]);
      break;
    case SV_SCALE:
      glScalef(arg[0],arg[1],arg[2]);
      break;
    case SV_DRAWCUBE:
      drawcube(arg[0],rgbptr);
      rgbptr=NULL;
      break;
    case SV_DRAWDISK:
      drawdisk(arg[0],arg[1], rgbptr);
      rgbptr=NULL;
      break;
    case SV_DRAWARCDISK:
      drawarcdisk(arg[0],arg[1], arg[2], rgbptr);
      rgbptr=NULL;
      break;
    case SV_DRAWCDISK:
      drawcdisk(arg[0],arg[1], rgbptr);
      rgbptr=NULL;
      break;
    case SV_DRAWHEXDISK:
      drawhexdisk(arg[0],arg[1], rgbptr);
      rgbptr=NULL;
      break;
    case SV_DRAWPOLYDISK:
      {
        int nsides;
  
        nsides = arg[0]+0.5;
        drawpolydisk(nsides, arg[1],arg[2], rgbptr);
        rgbptr=NULL;
      }
      break;
    case SV_DRAWRING:
      drawring(arg[0],arg[1], arg[2], rgbptr);
      rgbptr=NULL;
      break;
    case SV_DRAWNOTCHPLATE:
      drawnotchplate(arg[0],arg[1], arg[2], arg[3], rgbptr);
      rgbptr=NULL;
      break;
    case SV_DRAWTRUNCCONE:
      drawtrunccone(arg[0],arg[1],arg[2], rgbptr);
      rgbptr=NULL;
      break;
    case SV_DRAWCONE:
      drawcone(arg[0],arg[1], rgbptr);
      rgbptr=NULL;
      break;
    case SV_DRAWTSPHERE:
      {
        drawtsphere(arg[0],arg[1],rgbptr);
      }
      rgbptr=NULL;
      break;
    case SV_DRAWSPHERE:
      drawsphere(arg[0],rgbptr);
      rgbptr=NULL;
      break;
    case SV_DRAWCIRCLE:
      drawcircle(arg[0],rgbptr);
      rgbptr=NULL;
      break;
    case SV_DRAWARC:
      drawarc(arg[0],arg[1],rgbptr);
      rgbptr=NULL;
      break;
    case SV_DRAWPOINT:
      drawpoint(rgbptr);
      rgbptr=NULL;
      break;
    case SV_SETCOLOR:
      {
        if(setbw==1){
          float grey;

          grey = color2bw(arg);
          rgbcolor[0]=255*grey;
          rgbcolor[1]=255*grey;
          rgbcolor[2]=255*grey;
          rgbcolor[3]=255;
        }
        else{
          rgbcolor[0]=255*arg[0];
          rgbcolor[1]=255*arg[1];
          rgbcolor[2]=255*arg[2];
          rgbcolor[3]=255;
        }
        if(select_device_color_ptr==NULL){
          rgbptr=rgbcolor;
        }
        else{
          rgbptr=select_device_color_ptr;
        }
      }
      break;
    case SV_SETLINEWIDTH:
      {
        glLineWidth(arg[0]);
      }
      break;
    case SV_SETPOINTSIZE:
      {
        glPointSize(arg[0]);
      }
      break;
    case SV_SETBW:
      {
        rgbcolor[0]=255*arg[0];
        rgbcolor[1]=255*arg[0];
        rgbcolor[2]=255*arg[0];
        rgbcolor[3]=255;
        if(select_device_color_ptr==NULL){
          rgbptr=rgbcolor;
        }
        else{
          rgbptr=select_device_color_ptr;
        }
      }
      break;
    case SV_DRAWLINE:
      drawline(arg,arg+3,rgbptr);
      rgbptr=NULL;
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
  }
  if(select_device_color==NULL){
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
  }

  if(object->use_displaylist==1&&displaylist_id!=0){
    glEndList();
  }

  glPopMatrix();

}

/* ----------------------- drawline ----------------------------- */

void drawline(float *xyz1, float *xyz2, unsigned char *rgbcolor){
  glBegin(GL_LINES);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
  glVertex3fv(xyz1);
  glVertex3fv(xyz2);
  glEnd();
}

/* ----------------------- drawsphere ----------------------------- */

void drawtsphere(float texture_index_ptr,float diameter, unsigned char *rgbcolor){
  texture *texti;
  float latitude, longitude;
  int texture_index;

  if(texture_index_ptr>=0.0&&texture_index_ptr<ntexturestack){
    texture_index=valstack[nvalstack+(int)(texture_index_ptr+0.5)]+0.5;
  }
  else{
    texture_index=-1.0;
  }
  if(texture_index<0||texture_index>ntextures-1){
    texti=NULL;
  }
  else{
    int itext;

    texti = textureinfo + texture_index;
    if(texti->loaded==0||texti->display==0)texti=NULL;
  }
  if(texti==NULL){
    drawsphere(diameter,rgbcolor);
  }
  else{
    int i,j;

    glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
    glEnable(GL_TEXTURE_2D);

    glBindTexture(GL_TEXTURE_2D,texti->name);

    glPushMatrix();
    glScalef(diameter/2.0,diameter/2.0,diameter/2.0);
    if(cos_lat==NULL)initspheresegs(NLAT,NLONG);
    glBegin(GL_QUADS);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
    for(j=0;j<NLAT;j++){
      float ti,tip1;
      float tj,tjp1;

      tj = 1.0-(float)j/NLAT;
      tjp1 = 1.0-(float)(j+1)/NLAT;
      for(i=0;i<NLONG;i++){
        float x, y, z;

        ti = 1.0-(float)i/NLONG;
        tip1 = 1.0-(float)(i+1)/NLONG;

        x = cos_long[i]*cos_lat[j];
        y = sin_long[i]*cos_lat[j];
        z = sin_lat[j];

        glNormal3f(x,y,z);
        glTexCoord2f(ti,tj);
        glVertex3f(x,y,z);

        x = cos_long[i+1]*cos_lat[j];
        y = sin_long[i+1]*cos_lat[j];
        z = sin_lat[j];
        glNormal3f(x,y,z);
        glTexCoord2f(tip1,tj);
        glVertex3f(x,y,z);

        x = cos_long[i+1]*cos_lat[j+1];
        y = sin_long[i+1]*cos_lat[j+1];
        z = sin_lat[j+1];
        glNormal3f(x,y,z);
        glTexCoord2f(tip1,tjp1);
        glVertex3f(x,y,z);

        x = cos_long[i]*cos_lat[j+1];
        y = sin_long[i]*cos_lat[j+1];
        z = sin_lat[j+1];

        glNormal3f(x,y,z);
        glTexCoord2f(ti,tjp1);
        glVertex3f(x,y,z);
      }
    }
    glEnd();
    glPopMatrix();
    if(texti!=NULL){
      glDisable(GL_TEXTURE_2D);
    }
  }
}

/* ----------------------- drawsphere ----------------------------- */

void drawsphere(float diameter, unsigned char *rgbcolor){
  int i,j;

  if(cos_lat==NULL)initspheresegs(NLAT,NLONG);

  glPushMatrix();
  glScalef(diameter/2.0,diameter/2.0,diameter/2.0);

  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
  for(j=0;j<NLAT;j++){
    for(i=0;i<NLONG;i++){
      float x, y, z;

      x = cos_long[i]*cos_lat[j];
      y = sin_long[i]*cos_lat[j];
      z = sin_lat[j];

      glNormal3f(x,y,z);
      glVertex3f(x,y,z);

      x = cos_long[i+1]*cos_lat[j];
      y = sin_long[i+1]*cos_lat[j];
      z = sin_lat[j];
      glNormal3f(x,y,z);
      glVertex3f(x,y,z);

      x = cos_long[i+1]*cos_lat[j+1];
      y = sin_long[i+1]*cos_lat[j+1];
      z = sin_lat[j+1];
      glNormal3f(x,y,z);
      glVertex3f(x,y,z);

      x = cos_long[i]*cos_lat[j+1];
      y = sin_long[i]*cos_lat[j+1];
      z = sin_lat[j+1];

      glNormal3f(x,y,z);
      glVertex3f(x,y,z);
    }
  }
  glEnd();
  glPopMatrix();
}

/* ----------------------- drawpoint ----------------------------- */

void drawpoint(unsigned char *rgbcolor){
  glBegin(GL_POINTS);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
  glVertex3f(0.0,0.0,0.0);
  glEnd();
}

/* ----------------------- drawcircle ----------------------------- */

void drawcircle(float diameter,unsigned char *rgbcolor){
  int i;

  if(ncirc==0)initcircle(CIRCLE_SEGS);
  glBegin(GL_LINE_LOOP);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
  for(i=0;i<ncirc;i++){
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0);
  }
  glEnd();
}

/* ----------------------- drawarc ----------------------------- */

void drawarc(float angle, float diameter,unsigned char *rgbcolor){
  int i, iarc;

  if(ncirc==0)initcircle(CIRCLE_SEGS);
  iarc = CIRCLE_SEGS*(angle+180.0/CIRCLE_SEGS)/360.0;
  if(iarc<2)iarc=2;
  if(iarc>CIRCLE_SEGS)iarc=CIRCLE_SEGS;
  glBegin(GL_LINE_LOOP);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
  for(i=0;i<iarc;i++){
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0);
  }
  glEnd();
}

/* ----------------------- drawcube ----------------------------- */

void drawcube(float size, unsigned char *rgbcolor){
  float s2;

  s2 = size/2.0;

  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);


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

void drawring(float diam_inner, float diam_outer, float height, unsigned char *rgbcolor){
  int i;

  if(ncirc==0)initcircle(CIRCLE_SEGS);
  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

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

void drawdisk(float diameter, float height, unsigned char *rgbcolor){
  int i;

  if(ncirc==0)initcircle(CIRCLE_SEGS);
  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

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
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

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

/* ----------------------- drawarcdisk ----------------------------- */

void drawarcdisk(float angle, float diameter, float height, unsigned char *rgbcolor){
  int i, iarc;

  if(ncirc==0)initcircle(CIRCLE_SEGS);

  iarc = CIRCLE_SEGS*angle/360.0 + 0.5;
  if(iarc<2)iarc=2;
  if(iarc>CIRCLE_SEGS)iarc=CIRCLE_SEGS;

  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

  for(i=0;i<iarc;i++){
    glNormal3f(xcirc[i],ycirc[i],0.0);
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1

    glNormal3f(xcirc[i+1],ycirc[i+1],0.0);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,0.0); // 2

    glNormal3f(xcirc[i+1],ycirc[i+1],0.0);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height); // 3

    glNormal3f(xcirc[i],ycirc[i],0.0);
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height); // 4
  }
  glNormal3f(0.0,-1.0,0.0);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(diameter*xcirc[  0]/2.0,diameter*ycirc[  0]/2.0,0.0); // 1
  glVertex3f(diameter*xcirc[  0]/2.0,diameter*ycirc[  0]/2.0,height); // 1
  glVertex3f(0.0,0.0,height);

  glNormal3f(-ycirc[iarc-1],xcirc[iarc-1],0.0);
  glVertex3f(0.0,0.0,0.0);
  glVertex3f(diameter*xcirc[  iarc-1]/2.0,diameter*ycirc[  iarc-1]/2.0,0.0); // 1
  glVertex3f(diameter*xcirc[  iarc-1]/2.0,diameter*ycirc[  iarc-1]/2.0,height); // 1
  glVertex3f(0.0,0.0,height);
  glEnd();

  glBegin(GL_TRIANGLES);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

  glNormal3f(0.0,0.0,-1.0);
  for(i=0;i<iarc;i++){
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0);
    glVertex3f(                    0.0,                    0.0,0.0);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,0.0);
  }
  glNormal3f(0.0,0.0,1.0);
  for(i=0;i<iarc;i++){
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height);
    glVertex3f(                    0.0,                    0.0, height);
  }
  glEnd();

}

/* ----------------------- drawcdisk ----------------------------- */

void drawcdisk(float diameter, float height, unsigned char *rgbcolor){
  int i;

  if(ncirc==0)initcircle(CIRCLE_SEGS);

  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

  for(i=0;i<ncirc;i++){
    glNormal3f(xcirc[i],ycirc[i],0.0);
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,-height/2.00); // 1

    glNormal3f(xcirc[i+1],ycirc[i+1],0.0);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,-height/2.0); // 2

    glNormal3f(xcirc[i+1],ycirc[i+1],0.0);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height/2.0); // 3

    glNormal3f(xcirc[i],ycirc[i],0.0);
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height/2.0); // 4
  }
  glEnd();

  glBegin(GL_TRIANGLES);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

  glNormal3f(0.0,0.0,-1.0);
  for(i=0;i<ncirc;i++){
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,-height/2.0);
    glVertex3f(                    0.0,                    0.0,-height/2.0);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,-height/2.0);
  }
  glNormal3f(0.0,0.0,1.0);
  for(i=0;i<ncirc;i++){
    glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height/2.0);
    glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height/2.0);
    glVertex3f(                    0.0,                    0.0, height/2.0);
  }
  glEnd();

}

/* ----------------------- drawhexdisk ----------------------------- */

void drawpolydisk(int nsides, float diameter, float height, unsigned char *rgbcolor){
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
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

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
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

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

void drawhexdisk(float diameter, float height, unsigned char *rgbcolor){
  int i;

  float x[7]={0.866,0.0,-0.866,-0.866,0.0 ,0.866,0.866};
  float y[7]={0.5,  1.0, 0.5,  -0.5, -1.0,-0.5,  0.5};
  float xnorm[6]={0.500, -0.500, -1.0,-0.500,  0.500, 1.0};
  float ynorm[6]={0.866,  0.866,  0.0,-0.866, -0.866, 0.0};
  float radius;

  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

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
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

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

void drawnotchplate(float diameter, float height, float notchheight, float direction, unsigned char *rgbcolor){
  int i;
  float diameter2;

  diameter2 = diameter + notchheight;

  if(ncirc==0)initcircle(CIRCLE_SEGS);
  if(cullfaces==1)glDisable(GL_CULL_FACE);


  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

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
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

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

void drawcone(float d1, float height, unsigned char *rgbcolor){
  int i;
  float factor, denom, rad;
  float hdr;

  if(ncirc==0)initcircle(CIRCLE_SEGS);
  if(height<=0.0)height=0.0001;


  rad = d1/2.0;
  hdr = height/rad;
  denom = 1.0/sqrt(1.0+hdr*hdr);
  factor = hdr*denom;

  glBegin(GL_TRIANGLES);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

  for(i=0;i<ncirc;i++){
    glNormal3f(factor*xcirc[i],factor*ycirc[i],denom);
    glVertex3f(rad*xcirc[  i],rad*ycirc[  i],0.0); // 1

    glNormal3f(factor*xcirc[i+1],factor*ycirc[i+1],denom);
    glVertex3f(rad*xcirc[i+1],rad*ycirc[i+1],0.0); // 2

    glNormal3f(factor*xcirc[i],factor*ycirc[i],denom);
    glVertex3f(0.0,0.0, height); // 3
  }
  glNormal3f(0.0,0.0,-1.0);
  for(i=0;i<ncirc;i++){
    glVertex3f(rad*xcirc[  i],rad*ycirc[  i],0.0);
    glVertex3f(                    0.0,                    0.0,0.0);
    glVertex3f(rad*xcirc[i+1],rad*ycirc[i+1],0.0);
  }
  glEnd();
}

/* ----------------------- drawtrunccone ----------------------------- */

void drawtrunccone(float d1, float d2, float height, unsigned char *rgbcolor){
  int i;
  float dz;

  if(ncirc==0)initcircle(CIRCLE_SEGS);
  if(height<=0.0)height=0.0001;

  dz = -(d2-d1)/height;

  glBegin(GL_QUADS);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

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
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

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
    if(STRCMP(label,objecti->label)==0){
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

/* ----------------------- initspheresegs ----------------------------- */

void initspheresegs(int nlat, int nlong){
  float dlat, dlong, pi;
  int i;

  FREEMEMORY(cos_lat);
  FREEMEMORY(sin_lat);
  FREEMEMORY(cos_long);
  FREEMEMORY(sin_long);
  NewMemory( (void **)&cos_lat,(nlat+1)*sizeof(float));
  NewMemory( (void **)&sin_lat,(nlat+1)*sizeof(float));
  NewMemory( (void **)&cos_long,(nlong+1)*sizeof(float));
  NewMemory( (void **)&sin_long,(nlong+1)*sizeof(float));
  pi=4.0*atan(1.0);

  dlat=pi/(float)nlat;
  for(i=0;i<=nlat;i++){
    float angle;

    angle = -pi/2.0 + i*dlat;
    cos_lat[i] = cos(angle);
    sin_lat[i] = sin(angle);
  }

  dlong=2.0*pi/(float)nlong;
  for(i=0;i<nlong;i++){
    float angle;

    angle = i*dlong;
    cos_long[i] = cos(angle);
    sin_long[i] = sin(angle);
  }
  cos_long[nlong]=cos_long[0];
  sin_long[nlong]=sin_long[0];
}

/* ----------------------- freecircle ----------------------------- */

void freecircle(void){
  FREEMEMORY(xcirc);
  FREEMEMORY(ycirc);
  ncirc=0;
}

/* ----------------------- init_SVOBJECT1 ----------------------------- */

sv_object *init_SVOBJECT1(char *label, char *commands, int visible){
  sv_object *object;
  sv_object_frame *framei;
  int eof;

  NewMemory( (void **)&object,sizeof(sv_object));
  object->use_displaylist=1;
  object->select_mode=0;
  object->used=0;
  object->visible=visible;
  strcpy(object->label,label);
  object->nframes=1;
  object->obj_frames=NULL;
  NewMemory((void **)&framei,sizeof(sv_object_frame));
  NewMemory((void **)&object->obj_frames,object->nframes*sizeof(sv_object_frame *));

  object->obj_frames[0]=framei;
  framei->device=object;
  parse_device_frame(commands, NULL, &eof, framei);
  framei->display_list_ID=-1;
  framei->error=0;
  framei->use_bw=setbw;

  return object;
}

/* ----------------------- make_error_frame ----------------------------- */

void make_error_frame(void){
  char buffer[256];
  int eof;

  error_frame=NULL;
  NewMemory((void **)&error_frame,sizeof(sv_object_frame));
  strcpy(buffer,"1.0 0.0 0.0 setcolor 0.1 drawsphere");
  error_frame->device=NULL;
  parse_device_frame(buffer, NULL, &eof, error_frame);
  error_frame->display_list_ID=-1;
}

/* ----------------------- init_SMVOBJECT2 ----------------------------- */

sv_object *init_SVOBJECT2(char *label, char *commandsoff, char *commandson, int visible){
  sv_object *object;
  int i;

  NewMemory( (void **)&object,sizeof(sv_object));
  object->use_displaylist=1;
  object->select_mode=0;
  object->used=0;
  object->visible=visible;
  strcpy(object->label,label);
  object->nframes=2;
  object->obj_frames=NULL;
  NewMemory((void **)&object->obj_frames,object->nframes*sizeof(sv_object_frame *));

  for(i=0;i<object->nframes;i++){
    sv_object_frame *framei;
    int eof;

    if(i==0){
      NewMemory((void **)&framei,sizeof(sv_object_frame));
      object->obj_frames[0]=framei;
      framei->error=0;
      framei->device=object;
      parse_device_frame(commandsoff, NULL, &eof, framei);
      framei->display_list_ID=-1;
      framei->use_bw=setbw;
    }
    else{
      NewMemory((void **)&framei,sizeof(sv_object_frame));
      object->obj_frames[1]=framei;
      framei->device=object;
      parse_device_frame(commandson, NULL, &eof, framei);
      framei->error=0;
      framei->display_list_ID=-1;
      framei->use_bw=setbw;
    }
  }
  return object;
}

/* ----------------------- gettoken ----------------------------- */

void get_token_id(char *token, int *opptr, int *num_opptr, int *num_outopptr, int *use_displaylist){

  int op, num_op, num_outop;
  
  *use_displaylist=0;
  
  if(STRCMP(token,"translate")==0){
    op=SV_TRANSLATE;
    num_op=SV_TRANSLATE_NUMARGS;
    num_outop=SV_TRANSLATE_NUMOUTARGS;
  }
  else if(STRCMP(token,"offsetx")==0){
    op=SV_OFFSETX;
    num_op=SV_OFFSETX_NUMARGS;
    num_outop=SV_OFFSETX_NUMOUTARGS;
  }
  else if(STRCMP(token,"offsety")==0){
    op=SV_OFFSETY;
    num_op=SV_OFFSETY_NUMARGS;
    num_outop=SV_OFFSETY_NUMOUTARGS;
  }
  else if(STRCMP(token,"offsetz")==0){
    op=SV_OFFSETZ;
    num_op=SV_OFFSETZ_NUMARGS;
    num_outop=SV_OFFSETZ_NUMOUTARGS;
  }
  else if(STRCMP(token,"rotatex")==0){
    op=SV_ROTATEX;
    num_op=SV_ROTATEX_NUMARGS;
    num_outop=SV_ROTATEX_NUMOUTARGS;
  }
  else if(STRCMP(token,"rotatey")==0){
    op=SV_ROTATEY;
    num_op=SV_ROTATEY_NUMARGS;
    num_outop=SV_ROTATEY_NUMOUTARGS;
  }
  else if(STRCMP(token,"rotatez")==0){
    op=SV_ROTATEZ;
    num_op=SV_ROTATEZ_NUMARGS;
    num_outop=SV_ROTATEZ_NUMOUTARGS;
  }
  else if(STRCMP(token,"scalexyz")==0){
    op=SV_SCALEXYZ;
    num_op=SV_SCALEXYZ_NUMARGS;
    num_outop=SV_SCALEXYZ_NUMOUTARGS;
  }
  else if(STRCMP(token,"scale")==0&&STRCMP(token,"scalexyz")!=0){
    op=SV_SCALE;
    num_op=SV_SCALE_NUMARGS;
    num_outop=SV_SCALE_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawcube")==0){
    op=SV_DRAWCUBE;
    num_op=SV_DRAWCUBE_NUMARGS;
    num_outop=SV_DRAWCUBE_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawdisk")==0){
    op=SV_DRAWDISK;
    num_op=SV_DRAWDISK_NUMARGS;
    num_outop=SV_DRAWDISK_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawcdisk")==0){
    op=SV_DRAWCDISK;
    num_op=SV_DRAWCDISK_NUMARGS;
    num_outop=SV_DRAWCDISK_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawhexdisk")==0){
    op=SV_DRAWHEXDISK;
    num_op=SV_DRAWHEXDISK_NUMARGS;
    num_outop=SV_DRAWHEXDISK_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawpolydisk")==0){
    op=SV_DRAWPOLYDISK;
    num_op=SV_DRAWPOLYDISK_NUMARGS;
    num_outop=SV_DRAWPOLYDISK_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawring")==0){
    op=SV_DRAWRING;
    num_op=SV_DRAWRING_NUMARGS;
    num_outop=SV_DRAWRING_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawnotchplate")==0){
    op=SV_DRAWNOTCHPLATE;
    num_op=SV_DRAWNOTCHPLATE_NUMARGS;
    num_outop=SV_DRAWNOTCHPLATE_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawtrunccone")==0){
    op=SV_DRAWTRUNCCONE;
    num_op=SV_DRAWTRUNCCONE_NUMARGS;
    num_outop=SV_DRAWTRUNCCONE_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawcone")==0){
    op=SV_DRAWCONE;
    num_op=SV_DRAWCONE_NUMARGS;
    num_outop=SV_DRAWCONE_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawtsphere")==0){
    op=SV_DRAWTSPHERE;
    *use_displaylist=0;
    num_op=SV_DRAWTSPHERE_NUMARGS;
    num_outop=SV_DRAWTSPHERE_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawsphere")==0){
    op=SV_DRAWSPHERE;
    num_op=SV_DRAWSPHERE_NUMARGS;
    num_outop=SV_DRAWSPHERE_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawline")==0){
    op=SV_DRAWLINE;
    num_op=SV_DRAWLINE_NUMARGS;
    num_outop=SV_DRAWLINE_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawpoint")==0){
    op=SV_DRAWPOINT;
    num_op=SV_DRAWPOINT_NUMARGS;
    num_outop=SV_DRAWPOINT_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawcircle")==0){
    op=SV_DRAWCIRCLE;
    num_op=SV_DRAWCIRCLE_NUMARGS;
    num_outop=SV_DRAWCIRCLE_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawarc")==0){
    op=SV_DRAWARC;
    num_op=SV_DRAWARC_NUMARGS;
    num_outop=SV_DRAWARC_NUMOUTARGS;
  }
  else if(STRCMP(token,"setcolor")==0){
    op=SV_SETCOLOR;
    num_op=SV_SETCOLOR_NUMARGS;
    num_outop=SV_SETCOLOR_NUMOUTARGS;
  }
  else if(STRCMP(token,"setlinewidth")==0){
    op=SV_SETLINEWIDTH;
    num_op=SV_SETLINEWIDTH_NUMARGS;
    num_outop=SV_SETLINEWIDTH_NUMOUTARGS;
  }
  else if(STRCMP(token,"setpointsize")==0){
    op=SV_SETPOINTSIZE;
    num_op=SV_SETPOINTSIZE_NUMARGS;
    num_outop=SV_SETPOINTSIZE_NUMOUTARGS;
  }
  else if(STRCMP(token,"setbw")==0){
    op=SV_SETBW;
    num_op=SV_SETBW_NUMARGS;
    num_outop=SV_SETBW_NUMOUTARGS;
  }
  else if(STRCMP(token,"push")==0){
    op=SV_PUSH;
    num_op=SV_PUSH_NUMARGS;
    num_outop=SV_PUSH_NUMOUTARGS;
  }
  else if(STRCMP(token,"pop")==0){
    op=SV_POP;
    num_op=SV_POP_NUMARGS;
    num_outop=SV_POP_NUMOUTARGS;
  }
  else if(STRCMP(token,"multiaddt")==0){
    op=SV_MULTIADDT;
    *use_displaylist=0;
    num_op=SV_MULTIADDT_NUMARGS;
    num_outop=SV_MULTIADDT_NUMOUTARGS;
  }
  else if(STRCMP(token,"clip")==0){
    op=SV_CLIP;
    *use_displaylist=0;
    num_op=SV_CLIP_NUMARGS;
    num_outop=SV_CLIP_NUMOUTARGS;
  }
  else if(STRCMP(token,"mirrorclip")==0){
    op=SV_MIRRORCLIP;
    *use_displaylist=0;
    num_op=SV_MIRRORCLIP_NUMARGS;
    num_outop=SV_MIRRORCLIP_NUMOUTARGS;
  }
  else if(STRCMP(token,"periodicclip")==0){
    op=SV_PERIODICCLIP;
    *use_displaylist=0;
    num_op=SV_PERIODICCLIP_NUMARGS;
    num_outop=SV_PERIODICCLIP_NUMOUTARGS;
  }
  else{
    op=SV_ERR;
    num_op=0;
 }
 *opptr=op;
 *num_opptr=num_op;
 *num_outopptr=num_outop;
}

/* ----------------------- get_token_loc ----------------------------- */

int get_token_loc(char *var,sv_object_frame *frame){
  int i;
  int return_val;

  for(i=0;i<frame->nsymbols;i++){
    int ii;
    tokendata *toki;
    char *token_var;

    ii = frame->symbols[i];
    toki = frame->tokens+ii;
    token_var = toki->token+1;
    if(strcmp(var,token_var)==0)return ii;
  }
  return -1;
}

/* ----------------------- parse_frame ----------------------------- */

char *parse_device_frame(char *buffer, FILE *stream, int *eof, sv_object_frame *frame){
  char  object_buffer[100000];
  char object_buffer2[100000];
  int ntokens;
  char *token;
  char *buffer_ptr=NULL,*buffer2;
  int i;
  int nsymbols,ncommands;

  *eof = 0;

  // concatentate frame

  trim(buffer);
  strcpy(object_buffer,buffer);
  while(stream!=NULL&&!feof(stream)){
    if(fgets(buffer,255,stream)==NULL){
      *eof=1;
      break;
    }
    remove_comment(buffer);
    trim(buffer);
    buffer2=trim_front(buffer);
    if(match(buffer2,"DEVICEDEF",9) == 1||
       match(buffer2,"AVATARDEF",9) == 1||
       match(buffer2,"NEWFRAME",8) == 1){
         buffer_ptr=buffer2;
         break;
    }
    strcat(object_buffer," ");
    strcat(object_buffer,buffer2);
  }
  strcpy(object_buffer2,object_buffer);

// count tokens

  ntokens=0;
  strcpy(object_buffer2,object_buffer);
  token=strtok(object_buffer2," ");
  while(token!=NULL){
    ntokens++;
    token=strtok(NULL," ");
  }
  frame->ntokens=ntokens;
  NewMemory((void **)&frame->tokens,ntokens*sizeof(tokendata));
  if(ntokens>0)NewMemory((void **)&frame->symbols,ntokens*sizeof(int));

  // count symbols and commands, zero out access counter

  strcpy(object_buffer2,object_buffer);
  token=strtok(object_buffer2," ");
  nsymbols=0;
  ncommands=0;
  for(i=0;i<ntokens;i++){
    tokendata *toki;
    char c;

    toki = frame->tokens + i;
    toki->token=token;
    toki->reads=0;

    c = token[0];

    if(c==':'){
      frame->symbols[nsymbols++]=i;
    }
    if(c>='a'&&c<='z'||c>='A'&&c<='Z')ncommands++;

    token=strtok(NULL," ");
  }
  frame->nsymbols=nsymbols;
  if(ntokens>0)NewMemory((void **)&frame->command_list,ntokens*sizeof(tokendata *));
  frame->ncommands=ncommands;

  // fill in token data structure

  nsymbols=0;
  ncommands=0;
  for(i=0;i<ntokens;i++){
    tokendata *toki, *first_token=NULL;
    char c;

    toki = frame->tokens + i;

    c = toki->token[0];
    if(first_token==NULL&&c!=':')first_token=toki;
    if(c>='a'&&c<='z'||c>='A'&&c<='Z'){
      int use_displaylist;
      int nargs_actual, noutargs_actual;
      tokendata *this_token, *last_token;

      toki->type=TOKEN_COMMAND;
      get_token_id(toki->token, &toki->command, &toki->nvars, &toki->noutvars, &use_displaylist);
      frame->command_list[ncommands]=toki;
      if(frame->device!=NULL)frame->device->use_displaylist=use_displaylist;
      if(ncommands>0){
        this_token=toki;
        last_token=frame->command_list[ncommands-1];
        nargs_actual = this_token-last_token - 1;
      }
      else{
        nargs_actual = toki-first_token;
        nargs_actual = toki->nvars;
      }
      if(nargs_actual!=toki->nvars){
        printf("*** error: The command %s in device %s has %i arguments, %i were expected\n",
          toki->token,frame->device->label,nargs_actual,toki->nvars);
      }
      if(nargs_actual==toki->nvars){
        int ii;

        noutargs_actual=0;
        for(ii=0;ii<nargs_actual;ii++){
          tokendata *tokii;

          tokii=toki-1-ii;
          if(tokii<frame->tokens)break;
          c=tokii->token[0];
          if(c!=':')continue;
          noutargs_actual++;
        }
        if(noutargs_actual!=toki->noutvars){
          printf("*** error: The command %s in device %s has %i output arguments, %i were expected\n",
            toki->token,frame->device->label,noutargs_actual,toki->noutvars);
        }
      }
      ncommands++;
    }
    else if(c=='$'){
      char vartoken[255];
      char *vartokenptr;
      tokendata *tokdest;

      toki->loc=get_token_loc(toki->token+1,frame);
      if(toki->loc>=0){
        tokdest = frame->tokens+toki->loc;
        toki->varptr=&tokdest->var;
        tokdest->reads++;
      }
      else{
        toki->varptr=NULL;
        printf("*** error: The label %s in device %s is not defined\n",toki->token,frame->device->label);
      }

      toki->type=TOKEN_GETVAL;
    }
    else if(c==':'){
      toki->type=TOKEN_FLOAT;
      toki->varptr=&toki->var;
      nsymbols++;
    }
    else{
      toki->type=TOKEN_FLOAT;
      sscanf(toki->token,"%f%",&toki->var);
      toki->varptr=&toki->var;
    }
  }
  for(i=0;i<ntokens;i++){
    tokendata *toki;
    char c;

    toki = frame->tokens + i;
    c=toki->token[0];
    if(c!=':')continue;
    if(toki->reads==0){
   //   printf("*** warning: token %s in device %s was not used\n",
   //     toki->token,frame->device->label);
    }
  }

  return buffer_ptr;
}

/* ----------------------- read_device_defs ----------------------------- */

int read_device_defs(char *file){
  FILE *stream;
  char buffer[256], *trim_buffer;
  char *buffer_ptr;
  sv_object *temp_object, *prev_object, *next_object, *current_object;
  sv_object_frame *current_frame;
  int firstdef;
  float *arglist;
  int *oplist, nargs, nops;
  sv_object *object_start, *objecti;
  size_t lenbuffer;
  int ndevices=0;
  int eof=0;

  stream=fopen(file,"r");
  if(stream==NULL)return 0;
  printf("      Reading device definitions from: %s\n",file);

  firstdef=-1;
  buffer_ptr=NULL;
  while(!feof(stream)){
    CheckMemory;
    if(buffer_ptr==NULL){
      if(eof==1||fgets(buffer,255,stream)==NULL)break;
      buffer_ptr=buffer;
    }
    remove_comment(buffer_ptr);
    trim(buffer_ptr);
    trim_buffer=trim_front(buffer_ptr);
    lenbuffer=strlen(buffer_ptr);
    if(lenbuffer<1){
      buffer_ptr=NULL;
      continue;
    }


    if(match(buffer_ptr,"DEVICEDEF",9) == 1||
       match(buffer_ptr,"AVATARDEF",9) == 1
      ){
        int is_avatar=0;
      char *label;

      sv_object_frame *first_frame, *last_frame;

      if(match(buffer_ptr,"AVATARDEF",9) == 1){
        is_avatar=1;
      }  
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
      current_object->use_displaylist=1;
      current_object->select_mode=0;
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
      current_object->type=is_avatar;

      first_frame->next=last_frame;
      first_frame->prev=NULL;
      last_frame->prev=first_frame;
      last_frame->next=NULL;

      firstdef=1;
      buffer_ptr=NULL;
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

      current_frame->display_list_ID=-1;
      current_frame->use_bw=setbw;
      current_frame->device=current_object;

      current_object->nframes++;

      firstdef=0;
      if(match(trim_buffer,"NEWFRAME",8)==1){
        buffer_ptr=NULL;
        continue;
      }
    }
    buffer_ptr=parse_device_frame(buffer,stream,&eof,current_frame);
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
        int iop, npushpop=0, ii;

        CheckMemory;
        objecti->obj_frames[j]=framei;
        framei->error=0;
        for(ii=0;ii<framei->ncommands;ii++){
          tokendata *command;
          char *c;
          int op;


          command = framei->command_list[ii];

          op = command->command;
          if(op==SV_PUSH){
            npushpop++;
          }
          else if(op==SV_POP){
            npushpop--;
            if(npushpop<0){
              npushpop=0;
              command->command=SV_NO_OP;
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
#ifdef pp_MESSAGE
  {
    char message[256];
    sprintf(message,"***Error: %i arguments were found (%i expected) for the token \"%s\" while parsing: \"%s\".\n",
      numargs_found,numargs_expected,token,buffer);
    error_message(message);
  }
#else
    printf("***Error: %i arguments were found (%i expected) for the token \"%s\" while parsing: \"%s\".\n",
      numargs_found,numargs_expected,token,buffer);
#endif
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
    if(STRCMP(objecti->label,label)==0)return objecti;
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
    if(framei->nsymbols>0)FREEMEMORY(framei->symbols);
    if(framei->ntokens>0)FREEMEMORY(framei->tokens);
    FREEMEMORY(framei);
    framei=next_frame;
  }
  FREEMEMORY(object);
}

/* ----------------------- remove_comment ----------------------------- */

void remove_comment(char *buffer){
  char *comment;
  comment = strstr(buffer,"//");
  if(comment!=NULL)comment[0]=0;
  return;
}

/* ----------------------- update_device_defs ----------------------------- */

void update_device_objects(void){
  int i;

  for(i=0;i<ndeviceinfo;i++){
    device *devicei;

    devicei = deviceinfo + i;

    devicei->object = get_SVOBJECT_type(devicei->label);
    if(devicei->object==NULL){
      devicei->object = device_defs_backup[0];
    }
    if(devicei->ntextures>0){
      int j;
      for(j=0;j<devicei->ntextures;j++){
        texture *texti;
        int texture_index;

        texti = devicei->textureinfo+j;
        texture_index = texti - textureinfo;
        if(texture_index<0||texture_index>ntextures-1)texture_index=-1;
        devicei->params[devicei->nparams+j]=texture_index;
      }
    }
  }

}

/* ----------------------- init_device_defs ----------------------------- */

void init_device_defs(void){
    char com_buffer[1024];
    char com_buffer2[1024];


    {
      char objectfile[1024];

      svofile_exists = 0;

      if(smvprogdir!=NULL){
        strcpy(objectfile,smvprogdir);
        strcat(objectfile,"devices.svo");
        read_device_defs(objectfile);
       // strcat(objectfile,"test.svo");
       // read_device_defs2(objectfile);
      }

      strcpy(objectfile,"devices.svo");
      read_device_defs(objectfile);

      strcpy(objectfile,fdsprefix);
      strcat(objectfile,".svo");
      read_device_defs(objectfile);

      init_avatar();
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
}

/* ----------------------- init_avatar ----------------------------- */

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

  strcpy(com_buffer,"0.0 0.0 1.0 translate 1.0 0.0 0.0 setcolor 0.03 0.1 drawdisk 0.0 0.0 1.0 setcolor 90.0 rotatey 0.03 0.2 drawdisk");
  avatar_defs_backup[0] = init_SVOBJECT1("Avatar_1", com_buffer,1);
  avatar_defs_backup[0]->type=1;

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

/* ----------------------- dist ----------------------------- */

float dist(float p1[3], float p2[3]){
  float dx, dy, dz;

  dx = p1[0] - p2[0];
  dy = p1[1] - p2[1];
  dz = p1[2] - p2[2];
  return sqrt(dx*dx+dy*dy+dz*dz);
}

/* ----------------------- init_avatar ----------------------------- */

float get_point2box_dist(float boxmin[3], float boxmax[3], float p1[3], float p2orig[3]){
  int i;
  float tt;
  int doit=0;
  float dx, dy, dz;
  float xx, yy, zz;
  float p2[3];

  // box - xmin, ymin, zmin, xmax, ymax, zmax

  // if p1 is outside of box then return dist(p1,p2)

  for(i=0;i<3;i++){
    if(p1[i]<boxmin[i])return dist(p1,p2orig);
    if(p1[i]>boxmax[i])return dist(p1,p2orig);
    p2[i]=p2orig[i];
  }

  // if p1 and p2 are both inside box then return dist(p1,p2)

  for(i=0;i<3;i++){
    if(p2[i]<boxmin[i]){
      doit=1;
      break;
    }
    if(p2[i]>boxmax[i]){
      doit=1;
      break;
    }
  }
  if(doit==0)return dist(p1,p2);

  dx = p2[0]-p1[0];
  dy = p2[1]-p1[1];
  dz = p2[2]-p1[2];

  if(p1[0]>=boxmin[0]&&boxmin[0]>=p2[0]){
    if(dx!=0.0){
      tt=(boxmin[0]-p1[0])/dx;
      xx = boxmin[0];
      yy = p1[1] + tt*dy;
      zz = p1[2] + tt*dz;
      if(boxmin[1]<=yy&&yy<=boxmax[1]&&boxmin[2]<=zz&&zz<=boxmax[2]){
        p2[0]=xx;
        p2[1]=yy;
        p2[2]=zz;
        return dist(p1,p2);
      }
    }
  }
  if(p1[0]<=boxmax[0]&&boxmax[0]<=p2[0]){
    if(dx!=0.0){
      tt=(boxmax[0]-p1[0])/dx;
      xx = boxmax[0];
      yy = p1[1] + tt*dy;
      zz = p1[2] + tt*dz;
      if(boxmin[1]<=yy&&yy<=boxmax[1]&&boxmin[2]<=zz&&zz<=boxmax[2]){
        p2[0]=xx;
        p2[1]=yy;
        p2[2]=zz;
        return dist(p1,p2);
      }
    }
  }
  if(p1[1]>=boxmin[1]&&boxmin[1]>=p2[1]){
    if(dy!=0.0){
      tt=(boxmin[1]-p1[1])/dy;
      xx = p1[0] + tt*dx;
      yy = boxmin[1];
      zz = p1[2] + tt*dz;
      if(boxmin[0]<=xx&&xx<=boxmax[0]&&boxmin[2]<=zz&&zz<=boxmax[2]){
        p2[0]=xx;
        p2[1]=yy;
        p2[2]=zz;
        return dist(p1,p2);
      }
    }
  }
  if(p1[1]<=boxmax[1]&&boxmax[1]<=p2[1]){
    if(dy!=0.0){
      tt=(boxmax[1]-p1[1])/dy;
      xx = p1[0] + tt*dx;
      yy = boxmax[1];
      zz = p1[2] + tt*dz;
      if(boxmin[0]<=xx&&xx<=boxmax[0]&&boxmin[2]<=zz&&zz<=boxmax[2]){
        p2[0]=xx;
        p2[1]=yy;
        p2[2]=zz;
        return dist(p1,p2);
      }
    }
  }
  if(p1[2]>=boxmin[2]&&boxmin[2]>=p2[2]){
    if(dz!=0.0){
      tt=(boxmin[2]-p1[2])/dz;
      xx = p1[0] + tt*dx;
      yy = p1[1] + tt*dy;
      zz = boxmin[2];
      if(boxmin[0]<=xx&&xx<=boxmax[0]&&boxmin[1]<=yy&&yy<=boxmax[1]){
        p2[0]=xx;
        p2[1]=yy;
        p2[2]=zz;
        return dist(p1,p2);
      }
    }
  }
  if(p1[2]<=boxmax[2]&&boxmax[2]<=p2[2]){
    if(dz!=0.0){
      tt=(boxmax[2]-p1[2])/dz;
      xx = p1[0] + tt*dx;
      yy = p1[1] + tt*dy;
      zz = boxmin[2];
      if(boxmin[0]<=xx&&xx<=boxmax[0]&&boxmin[1]<=yy&&yy<=boxmax[1]){
        p2[0]=xx;
        p2[1]=yy;
        p2[2]=zz;
        return dist(p1,p2);
      }
    }
  }
  ASSERT(FALSE);
  return dist(p1,p2);
}

/* ----------------------- dist2plane ------------------------ */

float dist2plane(float x, float y, float z, float xyzp[3], float xyzpn[3]){
  float return_val;
  float xyz[3];
  int i;

  xyz[0]=x;
  xyz[1]=y;
  xyz[2]=z;
  return_val=0.0;
  for(i=0;i<3;i++){
    return_val+=(xyz[i]-xyzp[i])*xyzpn[i];
  }
  return return_val;
}

/* ----------------------- init_device_plane ------------------------ */

void init_device_plane(device *devicei){
  int colorindex;
  int i;
  float level=0.0;
  float xx[2], yy[2], zz[2];

/* stuff min and max grid data into a more convenient form 
  assuming the following grid numbering scheme

       5-------6
     / |      /| 
   /   |     / | 
  4 -------7   |
  |    |   |   |  
  Z    1---|---2
  |  Y     |  /
  |/       |/
  0--X-----3     

  */
  if(devicei->plane_surface==NULL)return;
  if(devicei->color==NULL){
    float rgbcolor[4];

    rgbcolor[0]=1.0;
    rgbcolor[1]=0.0;
    rgbcolor[2]=0.0;
    rgbcolor[3]=1.0;
    devicei->color=getcolorptr(rgbcolor);
  }
  colorindex=0;
  for(i=0;i<nmeshes;i++){
    int j;
    mesh *meshi;
    float xvert[12], yvert[12], zvert[12];
    int triangles[18];
    int nvert, ntriangles;
    int nodeindexes[8], closestnodes[18];
    float vals[8];

    InitIsosurface(devicei->plane_surface[i],level,devicei->color,colorindex);
    devicei->plane_surface[i]->cullfaces=1;

    meshi = meshinfo + i;

    xx[0]=meshi->xbar0;
    xx[1]=xbar0+xyzmaxdiff*meshi->xbar;
    yy[0]=meshi->ybar0;
    yy[1]=ybar0+xyzmaxdiff*meshi->ybar;
    zz[0]=meshi->zbar0;
    zz[1]=zbar0+xyzmaxdiff*meshi->zbar;
    for(j=0;j<8;j++){
      nodeindexes[j]=j;
    }
    vals[0]=dist2plane(xx[0],yy[0],zz[0],devicei->xyz,devicei->xyznorm);
    vals[1]=dist2plane(xx[0],yy[1],zz[0],devicei->xyz,devicei->xyznorm);
    vals[2]=dist2plane(xx[1],yy[1],zz[0],devicei->xyz,devicei->xyznorm);
    vals[3]=dist2plane(xx[1],yy[0],zz[0],devicei->xyz,devicei->xyznorm);
    vals[4]=dist2plane(xx[0],yy[0],zz[1],devicei->xyz,devicei->xyznorm);
    vals[5]=dist2plane(xx[0],yy[1],zz[1],devicei->xyz,devicei->xyznorm);
    vals[6]=dist2plane(xx[1],yy[1],zz[1],devicei->xyz,devicei->xyznorm);
    vals[7]=dist2plane(xx[1],yy[0],zz[1],devicei->xyz,devicei->xyznorm);

    xx[0]=(meshi->xbar0-xbar0)/xyzmaxdiff;
    xx[1]=meshi->xbar;
    yy[0]=(meshi->ybar0-ybar0)/xyzmaxdiff;
    yy[1]=meshi->ybar;
    zz[0]=(meshi->zbar0-zbar0)/xyzmaxdiff;
    zz[1]=meshi->zbar;

    GetIsobox(xx, yy, zz, vals, NULL, nodeindexes, level,
              xvert, yvert, zvert, NULL, closestnodes, &nvert, triangles, &ntriangles);

    UpdateIsosurface(devicei->plane_surface[i], xvert, yvert, zvert, NULL,
                     closestnodes, nvert, triangles, ntriangles);
    GetNormalSurface(devicei->plane_surface[i]);
    CompressIsosurface(devicei->plane_surface[i],1,
          xbar0,2*xbar,ybar0,2*ybar,zbar0,zbar);
    SmoothIsoSurface(devicei->plane_surface[i]);
  }

}

/* ----------------------- init_device ----------------------------- */

void init_device(device *devicei, float *xyz, float *xyzn, int state0, int nparams, float *params, char *labelptr){
  float norm;
  int i;

  devicei->labelptr=devicei->label;
  devicei->color=NULL;
  devicei->line_width=1.0;
  if(labelptr!=NULL){
    strcpy(devicei->label,labelptr);
  }
  if(STRCMP(devicei->object->label,"plane")==0){
    float color[4];

    NewMemory( (void **)&devicei->plane_surface,nmeshes*sizeof(isosurface *));
    for(i=0;i<nmeshes;i++){
      NewMemory( (void **)&devicei->plane_surface[i],sizeof(isosurface));
    }
    if(nparams>=3){
      color[0]=params[0];
      color[1]=params[1];
      color[2]=params[2];
      color[3]=1.0;
      devicei->color=getcolorptr(color);
    }
    if(nparams>=4){
      devicei->line_width=params[3];
    }
  }
  else{
    devicei->plane_surface=NULL;
  }
  if(xyz!=NULL){
    devicei->xyz[0]=xyz[0];
    devicei->xyz[1]=xyz[1];
    devicei->xyz[2]=xyz[2];
  }
  norm = sqrt(xyzn[0]*xyzn[0]+xyzn[1]*xyzn[1]+xyzn[2]*xyzn[2]);
  if(norm!=0.0){
    devicei->xyznorm[0]=xyzn[0]/norm;
    devicei->xyznorm[1]=xyzn[1]/norm;
    devicei->xyznorm[2]=xyzn[2]/norm;
  }
  else{
    devicei->xyznorm[0]=0.0;
    devicei->xyznorm[1]=0.0;
    devicei->xyznorm[2]=1.0;
  }
  devicei->nstate_changes=0;
  devicei->istate_changes=0;
  devicei->act_times=NULL;
  devicei->state_values=NULL;
  devicei->showstatelist=NULL;
  devicei->act_time=-1.0;
  devicei->device_mesh=NULL;
  devicei->state0=state0;
  devicei->nparams=nparams;
  devicei->params=params;
  if(nparams>0&&params!=NULL){
    for(i=0;i<nparams;i++){
      devicei->params[i]=params[i];
    }
  }
}
