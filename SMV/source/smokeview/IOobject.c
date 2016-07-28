#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "datadefs.h"
#include "smokeviewvars.h"

#define VEL_INVALID 0
#define VEL_CARTESIAN 1
#define VEL_POLAR 2

#define VECTOR_LINE 0
#define VECTOR_ARROW 1
#define VECTOR_OBJECT 2
#define VECTOR_PROFILE 3

#define CIRCLE_SEGS 12

#define SV_TRANSLATE  100
#define SV_ROTATEX    101
#define SV_ROTATEY    102
#define SV_ROTATEZ    103
#define SV_SCALEXYZ   104
#define SV_SCALE      105
#define SV_SCALEAUTO  106
#define SV_SCALEGRID  107
#define SV_OFFSETX 108
#define SV_OFFSETY 109
#define SV_OFFSETZ 110
#define SV_MULTIADDT 112
#define SV_CLIP 113
#define SV_MIRRORCLIP 114
#define SV_PERIODICCLIP 115
#define SV_ADD 116
#define SV_SUB 117
#define SV_MULT 118
#define SV_DIV 119
#define SV_GETT 120
#define SV_IF 121
#define SV_ELSE 122
#define SV_ENDIF 123
#define SV_GT 124
#define SV_GE 125
#define SV_LT 126
#define SV_LE 127
#define SV_AND 128
#define SV_OR 129
#define SV_ABS 130
#define SV_EQ 131
#define SV_ROTATEXYZ 132
#define SV_GTRANSLATE  133
#define SV_ROTATEAXIS 134
#define SV_ROTATEEYE 135
#define SV_INCLUDE 136
#define SV_INCLUDEF 137
#define SV_RANDXY 138
#define SV_RANDXZ 139
#define SV_RANDYZ 140
#define SV_RANDXYZ 141
#define SV_ORIENX 142
#define SV_ORIENY 143
#define SV_ORIENZ 144
#define SV_CLIPX 146
#define SV_CLIPY 147
#define SV_CLIPZ 148
#define SV_CLIPOFF 149

#define SV_TRANSLATE_NUMARGS  3
#define SV_ROTATEX_NUMARGS    1
#define SV_ROTATEY_NUMARGS    1
#define SV_ROTATEZ_NUMARGS    1
#define SV_SCALEXYZ_NUMARGS   3
#define SV_SCALEAUTO_NUMARGS  1
#define SV_SCALEGRID_NUMARGS  1
#define SV_SCALE_NUMARGS      1
#define SV_OFFSETX_NUMARGS 1
#define SV_OFFSETY_NUMARGS 1
#define SV_OFFSETZ_NUMARGS 1
#define SV_MULTIADDT_NUMARGS 3
#define SV_CLIP_NUMARGS 4
#define SV_MIRRORCLIP_NUMARGS 4
#define SV_PERIODICCLIP_NUMARGS 4
#define SV_ADD_NUMARGS 3
#define SV_SUB_NUMARGS 3
#define SV_MULT_NUMARGS 3
#define SV_DIV_NUMARGS 3
#define SV_GETT_NUMARGS 1
#define SV_IF_NUMARGS 1
#define SV_ELSE_NUMARGS 0
#define SV_ENDIF_NUMARGS 0
#define SV_GT_NUMARGS 3
#define SV_GE_NUMARGS 3
#define SV_LT_NUMARGS 3
#define SV_LE_NUMARGS 3
#define SV_AND_NUMARGS 3
#define SV_OR_NUMARGS 3
#define SV_ABS_NUMARGS 2
#define SV_EQ_NUMARGS 2
#define SV_ROTATEXYZ_NUMARGS 3
#define SV_GTRANSLATE_NUMARGS  3
#define SV_ROTATEAXIS_NUMARGS 4
#define SV_ROTATEEYE_NUMARGS 0
#define SV_INCLUDE_NUMARGS 1
#define SV_INCLUDEF_NUMARGS 2
#define SV_RANDXY_NUMARGS 1
#define SV_RANDXZ_NUMARGS 1
#define SV_RANDYZ_NUMARGS 1
#define SV_RANDXYZ_NUMARGS 1
#define SV_ORIENX_NUMARGS 3
#define SV_ORIENY_NUMARGS 3
#define SV_ORIENZ_NUMARGS 3
#define SV_CLIPX_NUMARGS 4
#define SV_CLIPY_NUMARGS 4
#define SV_CLIPZ_NUMARGS 4
#define SV_CLIPOFF_NUMARGS 0

#define SV_TRANSLATE_NUMOUTARGS  0
#define SV_ROTATEX_NUMOUTARGS    0
#define SV_ROTATEY_NUMOUTARGS    0
#define SV_ROTATEZ_NUMOUTARGS    0
#define SV_SCALEXYZ_NUMOUTARGS   0
#define SV_SCALEAUTO_NUMOUTARGS  0
#define SV_SCALEGRID_NUMOUTARGS  0
#define SV_SCALE_NUMOUTARGS      0
#define SV_OFFSETX_NUMOUTARGS 0
#define SV_OFFSETY_NUMOUTARGS 0
#define SV_OFFSETZ_NUMOUTARGS 0
#define SV_MULTIADDT_NUMOUTARGS 1
#define SV_CLIP_NUMOUTARGS 1
#define SV_MIRRORCLIP_NUMOUTARGS 1
#define SV_PERIODICCLIP_NUMOUTARGS 1
#define SV_ADD_NUMOUTARGS 1
#define SV_SUB_NUMOUTARGS 1
#define SV_MULT_NUMOUTARGS 1
#define SV_DIV_NUMOUTARGS 1
#define SV_GETT_NUMOUTARGS 1
#define SV_IF_NUMOUTARGS 0
#define SV_ELSE_NUMOUTARGS 0
#define SV_ENDIF_NUMOUTARGS 0
#define SV_GT_NUMOUTARGS 1
#define SV_GE_NUMOUTARGS 1
#define SV_LT_NUMOUTARGS 1
#define SV_LE_NUMOUTARGS 1
#define SV_AND_NUMOUTARGS 1
#define SV_OR_NUMOUTARGS 1
#define SV_ABS_NUMOUTARGS 1
#define SV_EQ_NUMOUTARGS 0
#define SV_ROTATEXYZ_NUMOUTARGS 0
#define SV_GTRANSLATE_NUMOUTARGS  0
#define SV_ROTATEAXIS_NUMOUTARGS 0
#define SV_ROTATEEYE_NUMOUTARGS 0
#define SV_INCLUDE_NUMOUTARGS 0
#define SV_INCLUDEF_NUMOUTARGS 0
#define SV_RANDXY_NUMOUTARGS 0
#define SV_RANDXZ_NUMOUTARGS 0
#define SV_RANDYZ_NUMOUTARGS 0
#define SV_RANDXYZ_NUMOUTARGS 0
#define SV_ORIENX_NUMOUTARGS 0
#define SV_ORIENY_NUMOUTARGS 0
#define SV_ORIENZ_NUMOUTARGS 0
#define SV_CLIPX_NUMOUTARGS 0
#define SV_CLIPY_NUMOUTARGS 0
#define SV_CLIPZ_NUMOUTARGS 0
#define SV_CLIPOFF_NUMOUTARGS 0

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
#define SV_DRAWSQUARE    216
#define SV_DRAWVENT      217
#define SV_DRAWCUBEC     218
#define SV_DRAWHSPHERE   219
#define SV_DRAWTRIBLOCK   220
#define SV_DRAWFILLEDCIRCLE    221

#define SV_DRAWCUBE_NUMARGS      1
#define SV_DRAWSPHERE_NUMARGS    1
#define SV_DRAWDISK_NUMARGS      2
#define SV_DRAWLINE_NUMARGS      6
#define SV_DRAWCIRCLE_NUMARGS    1
#define SV_DRAWFILLEDCIRCLE_NUMARGS    1
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
#define SV_DRAWSQUARE_NUMARGS 1
#define SV_DRAWVENT_NUMARGS 2
#define SV_DRAWCUBEC_NUMARGS      1
#define SV_DRAWHSPHERE_NUMARGS    1
#define SV_DRAWTRIBLOCK_NUMARGS    2

#define SV_DRAWCUBE_NUMOUTARGS      0
#define SV_DRAWSPHERE_NUMOUTARGS    0
#define SV_DRAWDISK_NUMOUTARGS      0
#define SV_DRAWLINE_NUMOUTARGS      0
#define SV_DRAWCIRCLE_NUMOUTARGS    0
#define SV_DRAWFILLEDCIRCLE_NUMOUTARGS    0
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
#define SV_DRAWSQUARE_NUMOUTARGS 0
#define SV_DRAWVENT_NUMOUTARGS 0
#define SV_DRAWCUBEC_NUMOUTARGS      0
#define SV_DRAWHSPHERE_NUMOUTARGS    0
#define SV_DRAWTRIBLOCK_NUMOUTARGS    0

#define SV_PUSH       300
#define SV_POP        301
#define SV_SETRGB   302
#define SV_SETRGBVAL   303
#define SV_SETBW      304
#define SV_SETLINEWIDTH 305
#define SV_SETPOINTSIZE 306
#define SV_SETCOLOR   307
#define SV_GETTEXTUREINDEX 308

#define SV_NO_OP      999

#define SV_PUSH_NUMARGS       0
#define SV_POP_NUMARGS        0
#define SV_SETRGB_NUMARGS   3
#define SV_SETRGBVAL_NUMARGS   3
#define SV_SETBW_NUMARGS      1
#define SV_SETLINEWIDTH_NUMARGS 1
#define SV_SETPOINTSIZE_NUMARGS 1
#define SV_SETCOLOR_NUMARGS   1
#define SV_GETTEXTUREINDEX_NUMARGS 2
#define SV_NO_OP_NUMARGS 0

#define SV_PUSH_NUMOUTARGS       0
#define SV_POP_NUMOUTARGS        0
#define SV_SETRGB_NUMOUTARGS   0
#define SV_SETRGBVAL_NUMOUTARGS   0
#define SV_SETBW_NUMOUTARGS      0
#define SV_SETLINEWIDTH_NUMOUTARGS 0
#define SV_SETPOINTSIZE_NUMOUTARGS 0
#define SV_SETCOLOR_NUMOUTARGS   0
#define SV_GETTEXTUREINDEX_NUMOUTARGS 1
#define SV_NO_OP_NUMOUTARGS 0


#define SV_ERR -1

#define NLAT device_sphere_segments
#define NLONG (2*device_sphere_segments)

#define TOKEN_FLOAT 0
#define TOKEN_COMMAND 1
#define TOKEN_GETVAL 2
#define TOKEN_STRING 3
#define TOKEN_TEXTURE 4

char *parse_device_frame(char *buffer, FILE *stream, int *eof, sv_object_frame *frame);
void reporterror(char *buffer, char *token, int numargs_found, int numargs_expected);
float get_point2box_dist(float boxmin[3], float boxmax[3], float p1[3], float p2[3]);

void drawsphereseg(float anglemin, float anglemax, float rmin, float rmax);
void rotateeye(void);
void rotateaxis(float angle, float ax, float ay, float az);
void rotatexyz(float x, float y, float z);
void drawcone(float d1, float height, unsigned char *rgbcolor);
void drawtrunccone(float d1, float d2, float height, unsigned char *rgbcolor);
void drawline(float *xyz1, float *xyz2, unsigned char *rgbcolor);
void drawarc(float angle, float diameter, unsigned char *rgbcolor);
void drawpoint(unsigned char *rgbcolor);
void drawsphere(float diameter, unsigned char *rgbcolor);
void drawhsphere(float diameter, unsigned char *rgbcolor);
void drawtriblock(float size, float height, unsigned char *rgbcolor);
void drawtsphere(int texture_index, float diameter, unsigned char *rgbcolor);
void drawcube(float size, unsigned char *rgbcolor);
void drawsquare(float size, unsigned char *rgbcolor);
void drawvent(float width, float height, unsigned char *rgbcolor);
void drawcdisk(float diameter, float height, unsigned char *rgbcolor);
void drawdisk(float diameter, float height, unsigned char *rgbcolor);
void drawarcdisk(float angle, float diameter, float height, unsigned char *rgbcolor);
void drawhexdisk(float diameter, float height, unsigned char *rgbcolor);
void drawpolydisk(int nsides, float diameter, float height, unsigned char *rgbcolor);
void drawring(float d_inner, float d_outer, float height, unsigned char *rgbcolor);
void drawnotchplate(float diameter, float height, float notchheight, float direction, unsigned char *rgbcolor);
void draw_SVOBJECT(sv_object *object, int frame_index_local, propdata *prop, int recurse_level, float *valrgb, int vis_override);
void free_object(sv_object *object);
void freecircle(circdata *circinfo);

static float *cos_long=NULL, *sin_long=NULL, *cos_lat=NULL, *sin_lat=NULL;
static float specular[4]={0.4,0.4,0.4,1.0};
unsigned char *rgbimage=NULL;
int rgbsize=0;

/* ------------------ get_world_eyepos ------------------------ */

void get_world_eyepos(float *mm, float user_eyepos[3],float scaled_eyepos_local[3]){
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

  scaled_eyepos_local[0] = -(mm[0]*mm[12]+mm[1]*mm[13]+ mm[2]*mm[14])/mscale[0];
  scaled_eyepos_local[1] = -(mm[4]*mm[12]+mm[5]*mm[13]+ mm[6]*mm[14])/mscale[1];
  scaled_eyepos_local[2] = -(mm[8]*mm[12]+mm[9]*mm[13]+mm[10]*mm[14])/mscale[2];
  DENORMALIZE_XYZ(user_eyepos,scaled_eyepos);
}

/* ----------------------- getsmokesensors ----------------------------- */

void getsmokesensors(void){
  int doit, i;
  int width, height;

  width = screenWidth;
  height = screenHeight;

  doit=0;
  for(i=0;i<ndeviceinfo;i++){
    devicedata *devicei;
    char *label;

    devicei = deviceinfo + i;
    label = devicei->object->label;
    if(STRCMP(label,"smokesensor")!=0)continue;
    doit=1;
    break;
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
    devicedata *devicei;
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
    devicedata *devicei;
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
    devicedata *devicei;
    int *ijk;
    char *label;
    meshdata *device_mesh;

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
  devicedata *devicei;
  int i;
  float *xyz, *xyznorm;
  float white[3]={1.0,1.0,1.0};
  float black[3]={0.0,0.0,0.0};
  int doit=0;

  if(fontindex==SCALED_FONT)ScaleFont3D();
  glPushMatrix();
  glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
  glTranslatef(-xbar0,-ybar0,-zbar0);
  if(active_smokesensors==1&&show_smokesensors!=SMOKESENSORS_HIDDEN){
    getdevice_screencoords();
  }
  for(i=0;i<ndeviceinfo;i++){
    devicei = deviceinfo + i;

    if(devicei->object->visible==0)continue;
    xyz = devicei->xyz;
    xyznorm = devicei->xyznorm;
    if(active_smokesensors==1&&show_smokesensors!=SMOKESENSORS_HIDDEN&&STRCMP(devicei->object->label,"smokesensor")==0){
      char label[256];
      float val;
      int ival;

      if(doit==0){
        glBlendFunc(GL_ONE,GL_ZERO);
        doit=1;
      }
      switch(show_smokesensors){
        case SMOKESENSORS_0255:
          sprintf(label,"%i",devicei->visval);
          break;
        case SMOKESENSORS_01:
          val = devicei->visval/255.0;
          sprintf(label,"%.2f",val);
          trimzeros(label);
          break;
        case SMOKESENSORS_SCALED:
        case SMOKESENSORS_0INF:
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
        default:
          ASSERT(FFALSE);
          break;
      }
      if(devicei->visval>128){
        Output3Text(black,xyz[0]+0.2*xyznorm[0],xyz[1]+0.2*xyznorm[1],xyz[2]+0.2*xyznorm[2],label);
      }
      else{
        Output3Text(white,xyz[0]+0.2*xyznorm[0],xyz[1]+0.2*xyznorm[1],xyz[2]+0.2*xyznorm[2],label);
      }
    }
  }
  if(doit==1){
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  }
  glPopMatrix();
}

/* ----------------------- get_vdevice_vel ----------------------------- */

void get_vdevice_vel(float time_local, vdevicedata *vdevicei, float *vel, float *angle_local, float *dvel, float *dangle, int *velocity_type){
  float uvel=0.0, vvel=0.0, wvel=0.0;
  devicedata *udev, *vdev, *wdev;
  int validu=1,validv=1,validw=1;

  udev = vdevicei->udev;
  vdev = vdevicei->vdev;
  wdev = vdevicei->wdev;

  *velocity_type=VEL_INVALID;
  if(udev!=NULL){
    uvel=get_device_val(time_local,udev,&validu);
  }
  if(vdev!=NULL){
    vvel=get_device_val(time_local,vdev,&validv);
  }
  if(wdev!=NULL){
    wvel=get_device_val(time_local,wdev,&validw);
  }
  if(validu==1&&validv==1&&validw==1){
    vel[0]=uvel;
    vel[1]=vvel;
    vel[2]=wvel;
    *velocity_type=VEL_CARTESIAN;
  }
  if(vdevicei->veldev!=NULL&&vdevicei->angledev!=NULL){
    float  velocity,  ang;
    float dvelocity=0.0, dang=0.0;
    int  valid_velocity=0,  valid_angle=0;
    int dvalid_velocity=0, dvalid_angle=0;

    velocity=get_device_val(time_local,vdevicei->veldev,&valid_velocity);
    if(vdevicei->sd_veldev!=NULL){
      dvelocity=get_device_val(time_local,vdevicei->sd_veldev,&dvalid_velocity);
      if(dvalid_velocity==0)dvelocity=0.0;
    }
    ang=get_device_val(time_local,vdevicei->angledev,&valid_angle);
    if(vdevicei->sd_angledev!=NULL){
      dang=get_device_val(time_local,vdevicei->sd_angledev,&dvalid_angle);
      if(dvalid_angle==0)dang=0.0;
    }
    if(valid_velocity==1&&valid_angle==1){
      vel[0]=velocity;
      dvel[0]=dvelocity;
      angle_local[0]=ang;
      dangle[0]=dang;
      *velocity_type=VEL_POLAR;
    }
  }
}

#define IN_INTERVAL(IVAL) \
  if(time_local>=times_local[(IVAL)]&&time_local<=times_local[(IVAL)+1]){\
    if(time_local-times_local[(IVAL)]<times_local[(IVAL)+1]-time_local){\
      devicei->val=devicei->vals[(IVAL)];\
      *valid=devicei->valids[(IVAL)];\
    }\
    else{\
      devicei->val=devicei->vals[(IVAL)+1];\
      *valid=devicei->valids[(IVAL)+1];\
    }\
    devicei->ival=(IVAL);\
    return devicei->val;\
  }

/* ----------------------- get_devices_val ----------------------------- */

float get_device_val(float time_local, devicedata *devicei, int *valid){
  int nvals;
  int ival;
  float *times_local;

  nvals = devicei->nvals;
  ival = devicei->ival;
  times_local = devicei->times;

  if(nvals==0||times_local==NULL){
    *valid=0;
    return 0.0;
  }
  IN_INTERVAL(ival);
  if(ival<nvals-1){
    IN_INTERVAL(ival+1);
  }
  if(ival>0){
    IN_INTERVAL(ival-1);
  }

  if(time_local<=times_local[0]){
    devicei->val=devicei->vals[0];
    devicei->ival=0;
    *valid=devicei->valids[0];
  }
  else if(time_local>=times_local[nvals-1]){
    devicei->val=devicei->vals[nvals-1];
    devicei->ival=nvals-2;
    *valid=devicei->valids[nvals-1];
  }
  else{
    int low, mid, high;

    low = 0;
    high = nvals-1;

    while (high-low>1){
      mid = (low+high)/2;
      if(time_local>times_local[mid]){
        low=mid;
      }
      else{
        high=mid;
      }
    }
    devicei->ival=low;
    devicei->val=devicei->vals[low];
    *valid=devicei->valids[low];
  }

  return devicei->val;
}

/* ----------------------- get_device_color ----------------------------- */

unsigned char *get_device_color(devicedata *devicei, unsigned char *colorval,float valmin, float valmax){
  float val;
  int valid,colorindex;
  float *rgb_local;

  if(devicei==NULL||valmax<=valmin)return NULL;
  val=get_device_val(global_times[itimes],devicei,&valid);
  if(valid!=1)return NULL;
  val = (val-valmin)/(valmax-valmin);
  colorindex=CLAMP(255*val,1,254);
  rgb_local=current_colorbar->colorbar+3*colorindex;
  colorval[0]=255*rgb_local[0];
  colorval[1]=255*rgb_local[1];
  colorval[2]=255*rgb_local[2];
  return colorval;
}

/* ----------------------- Output_Device_Val ----------------------------- */

void Output_Device_Val(devicedata *devicei){
  char label[1000];
  float val;
  int valid;

  if(fontindex==SCALED_FONT)ScaleFont3D();
  val=get_device_val(global_times[itimes],devicei,&valid);
  if(valid==1){
    f_units *unitclass;
    char *unit;
    char valuelabel[100];

    unitclass = GetUnitClass(devicei->unit);
    unit=devicei->unit;
    if(unitclass!=NULL){
      char *unit_type;
      f_unit *funit;

      funit = unitclass->units+unitclass->unit_index;
      unit = funit->unit;
      unit_type = unitclass->unitclass;
      val = GetUnitVal(unit_type, val);
    }
    strcpy(label, "");
    sprintf(valuelabel, "%.1f",val);
    if(showdevicetype == 1){
      strcat(label, devicei->quantity);
      strcat(label, " ");
    }
    strcat(label, valuelabel);
    strcat(label, " ");
    if(showdeviceunit == 1)strcat(label, unit);
    Output3Text(foregroundcolor,0.0,0.0,0.0,label);
  }
  else{
    sprintf(label,"not available");
    Output3Text(foregroundcolor,0.0,0.0,0.0,label);
  }
}

/* ----------------------- draw_pilot ----------------------------- */

#ifdef pp_PILOT
void draw_pilot1(void){
  int i;

  if(showtime == 1 && itimes >= 0 && itimes<nglobal_times&&vispilot == 1 && nvdeviceinfo>0){
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, block_ambient2);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glEnable(GL_COLOR_MATERIAL);

    glPushMatrix();
    glScalef(SCALE2SMV(1.0), SCALE2SMV(1.0), SCALE2SMV(1.0));
    glTranslatef(-xbar0, -ybar0, -zbar0);
    glColor3fv(foregroundcolor);
    glLineWidth(vectorlinewidth);
    for(i = 0; i < nvdeviceinfo; i++){
      vdevicedata *vdevi;
      float *xyz;
      int k;
      pilotdata *piloti;
      float dangle;

      vdevi = vdeviceinfo + i;
      if(vdevi->unique == 0)continue;
      xyz = vdevi->valdev->xyz;

      piloti = &(vdevi->pilotinfo);

      dangle = 360.0 / (float)piloti->nbuckets;
      glBegin(GL_LINES);
      for(k = 0; k < piloti->nbuckets; k++){
        float angle, cosang, sinang;

        angle = (float)k*dangle;
        cosang = cos(DEG2RAD*angle);
        sinang = sin(DEG2RAD*angle);
        glVertex3f(xyz[0], xyz[1], xyz[2]);
        glVertex3f(xyz[0] - SCALE2FDS(piloti->fraction[k])*cosang, xyz[1] - SCALE2FDS(piloti->fraction[k])*sinang, xyz[2]);
      }
      glEnd();
    }
    glPopMatrix();
    glDisable(GL_LIGHTING);
  }
}

/* ----------------------- draw_pilot2 ----------------------------- */

void draw_pilot2(void){
  int i;

  if(showtime == 1 && itimes >= 0 && itimes<nglobal_times&&vispilot == 1 && nvdeviceinfo>0){
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, block_ambient2);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular);
    glEnable(GL_COLOR_MATERIAL);

    glPushMatrix();
    glScalef(SCALE2SMV(1.0), SCALE2SMV(1.0), SCALE2SMV(1.0));
    glTranslatef(-xbar0, -ybar0, -zbar0);
    glColor3fv(foregroundcolor);
    glLineWidth(vectorlinewidth);
    for(i = 0; i < nvdeviceinfo; i++){
      vdevicedata *vdevi;
      float *xyz;
      int k;
      pilotdata *piloti;
      float dangle;

      vdevi = vdeviceinfo + i;
      if(vdevi->unique == 0)continue;
      xyz = vdevi->valdev->xyz;

      piloti = &(vdevi->pilotinfo);

      dangle = 360.0 / (float)piloti->nbuckets;
      glBegin(GL_LINES);
      for(k = 0; k < piloti->nbuckets; k++){
        float angle, cosang, sinang;
        int kk;

        kk = k;
        angle = (float)kk*dangle;
        cosang = cos(DEG2RAD*angle);
        sinang = sin(DEG2RAD*angle);
        glVertex3f(xyz[0] - SCALE2FDS(piloti->fraction[kk])*cosang, xyz[1] - SCALE2FDS(piloti->fraction[kk])*sinang, xyz[2]);

        kk = k+1;
        if(kk == piloti->nbuckets - 1)kk = 0;
        angle = (float)kk*dangle;
        cosang = cos(DEG2RAD*angle);
        sinang = sin(DEG2RAD*angle);
        glVertex3f(xyz[0] - SCALE2FDS(piloti->fraction[kk])*cosang, xyz[1] - SCALE2FDS(piloti->fraction[kk])*sinang, xyz[2]);
      }
      glEnd();
    }
    glPopMatrix();
    glDisable(GL_LIGHTING);
  }
}

/* ----------------------- draw_pilot ----------------------------- */

void draw_pilot(void){
  switch(pilot_viewtype){
  case 0:
    draw_pilot1();
    break;
  case 1:
    draw_pilot2();
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
}
#endif

/* ----------------------- draw_devices ----------------------------- */

void draw_devices(void){
  int drawobjects_as_vectors;
  int ii;

  if(select_device==0||show_mode!=SELECTOBJECT){
    int i;

    for(i=0;i<ndeviceinfo;i++){
      devicedata *devicei;

      devicei = deviceinfo + i;
      if(devicei->object->visible==0)continue;
      if(devicei->in_zone_csv==1)continue;
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
  drawobjects_as_vectors=0;
  if(showtime==1&&itimes>=0&&itimes<nglobal_times&&showvdeviceval==1&&nvdeviceinfo>0){
    unsigned char arrow_color[4];
    float arrow_color_float[4];
    int j;

    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,specular);
    glEnable(GL_COLOR_MATERIAL);

    glPushMatrix();
    glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
    glTranslatef(-xbar0,-ybar0,-zbar0);
    glPointSize(vectorpointsize);
    arrow_color[0]=255*foregroundcolor[0];
    arrow_color[1]=255*foregroundcolor[1];
    arrow_color[2]=255*foregroundcolor[2];
    arrow_color[3]=255;
    glColor3ubv(arrow_color);
    for(j = 0; j < ntreedeviceinfo; j++){
      treedevicedata *treei;
      int i,first;

      treei = treedeviceinfo + j;
      if(vectortype==VECTOR_PROFILE&&treei->n<mintreesize)continue;
      first = 1;
      for(i = treei->first; i <= treei->last; i++){
        vdevicedata *vdevi;
        devicedata *devicei;
        float vel[3], angle, dvel, dangle;
        float *xyz;
        int velocity_type;
        vdevicesortdata *vdevsorti;

        vdevsorti = vdevices_sorted + i;
        if(vectortype == VECTOR_PROFILE){
          if(vdevsorti->dir == XDIR && vis_xtree == 0)continue;
          if(vdevsorti->dir == YDIR && vis_ytree == 0)continue;
          if(vdevsorti->dir == ZDIR && vis_ztree == 0)continue;
        }
        else{
          if(vdevsorti->dir != ZDIR)continue;
        }

        vdevi = vdevsorti->vdeviceinfo;
        devicei = vdevi->colordev;
        if(devicei == NULL)continue;
        if(vdevi->unique == 0)continue;
        xyz = vdevi->valdev->xyz;
        get_vdevice_vel(global_times[itimes], vdevi, vel, &angle, &dvel, &dangle, &velocity_type);
        if(colordeviceval == 1){
          int type, vistype = 0;

          type = devicei->type2;
          if(type >= 0 && type < ndevicetypes)vistype = devicetypes[type]->type2vis;
          if(vistype == 1){
            unsigned char color[4], *colorptr;

            colorptr = get_device_color(devicei, color, device_valmin, device_valmax);
            if(colorptr != NULL){
              arrow_color[0] = colorptr[0];
              arrow_color[1] = colorptr[1];
              arrow_color[2] = colorptr[2];
              arrow_color[3] = 255;
            }
            else{
              arrow_color[0] = 255 * foregroundcolor[0];
              arrow_color[1] = 255 * foregroundcolor[1];
              arrow_color[2] = 255 * foregroundcolor[2];
              arrow_color[3] = 255;
            }
            glColor3ubv(arrow_color);
          }
        }
        arrow_color_float[0] = (float)arrow_color[0] / 255.0;
        arrow_color_float[1] = (float)arrow_color[1] / 255.0;
        arrow_color_float[2] = (float)arrow_color[2] / 255.0;
        arrow_color_float[3] = (float)arrow_color[3] / 255.0;
        if(velocity_type == VEL_CARTESIAN){
          float xyz1_old[3], xyz2_old[3];
          float xyz1_new[3], xyz2_new[3];
          unsigned char arrow_color_old[4], *arrow_color_new;
          float dxyz[3], vec0[3] = {0.0, 0.0, 0.0}, zvec[3] = {0.0, 0.0, 1.0};
          float axis[3], speed;
          int state = 0;
          int jj;

          for(jj = 0; jj < 3; jj++){
            dxyz[jj] = 0.5*SCALE2FDS(vel[jj]) / max_dev_vel;
          }
          speed = sqrt(dxyz[0] * dxyz[0] + dxyz[1] * dxyz[1] + dxyz[2] * dxyz[2]);

          switch(vectortype){
          case VECTOR_PROFILE:
            xyz2_new[0] = xyz[0] + dxyz[0];
            xyz2_new[1] = xyz[1] + dxyz[1];
            xyz2_new[2] = xyz[2] + dxyz[2];

            xyz1_new[0] = xyz[0];
            xyz1_new[1] = xyz[1];
            xyz1_new[2] = xyz[2];

            arrow_color_new = arrow_color;
            if(first==1){
              first = 0;
              memcpy(xyz1_old, xyz1_new, 3 * sizeof(float));
              memcpy(xyz2_old, xyz2_new, 3 * sizeof(float));
              memcpy(arrow_color_old, arrow_color_new, 4 * sizeof(unsigned char));
              continue;
            }

            //  draw triangles for following rectangle

            //   xyz1_new---------xyz2_new
            //      |       /       |
            //   xyz1_old---------xyz2_old
            glBegin(GL_TRIANGLES);
            glColor3ubv(arrow_color_old);
            glVertex3fv(xyz1_old);
            glVertex3fv(xyz2_old);
            glColor3ubv(arrow_color_new);
            glVertex3fv(xyz2_new);

            glColor3ubv(arrow_color_old);
            glVertex3fv(xyz1_old);
            glColor3ubv(arrow_color_new);
            glVertex3fv(xyz2_new);
            glColor3ubv(arrow_color_old);
            glVertex3fv(xyz2_old);

            glColor3ubv(arrow_color_old);
            glVertex3fv(xyz1_old);
            glColor3ubv(arrow_color_new);
            glVertex3fv(xyz2_new);
            glVertex3fv(xyz1_new);

            glColor3ubv(arrow_color_old);
            glVertex3fv(xyz1_old);
            glColor3ubv(arrow_color_new);
            glVertex3fv(xyz1_new);
            glVertex3fv(xyz2_new);

            glEnd();

            memcpy(xyz1_old, xyz1_new, 3 * sizeof(float));
            memcpy(xyz2_old, xyz2_new, 3 * sizeof(float));
            memcpy(arrow_color_old, arrow_color_new, 4 * sizeof(unsigned char));
            break;
          case VECTOR_LINE:
            rotateu2v(zvec, dxyz, axis, &angle);
            glPushMatrix();
            glTranslatef(xyz[0], xyz[1], xyz[2]);
            glRotatef(RAD2DEG*angle, axis[0], axis[1], axis[2]);
            glScalef(1.0, 1.0, vector_baselength*speed);
            glColor3ubv(arrow_color);
            glBegin(GL_LINES);
            glVertex3fv(vec0);
            glVertex3fv(zvec);
            glEnd();
            glBegin(GL_POINTS);
            glVertex3fv(zvec);
            glEnd();
            glPopMatrix();
            break;
          case VECTOR_ARROW:
            rotateu2v(zvec, dxyz, axis, &angle);
            glPushMatrix();
            glTranslatef(xyz[0], xyz[1], xyz[2]);
            glRotatef(RAD2DEG*angle, axis[0], axis[1], axis[2]);
            glPushMatrix();
            glScalef(1.0, 1.0, speed*vector_baselength);
            drawdisk(vector_basediameter, 1.0, arrow_color);
            glPopMatrix();
            glTranslatef(0.0, 0.0, speed*vector_baselength);
            drawcone(vector_headdiameter, vector_headlength, arrow_color);
            glPopMatrix();
            break;
          case VECTOR_OBJECT:
            rotateu2v(zvec, dxyz, axis, &angle);
            glPushMatrix();
            glTranslatef(xyz[0], xyz[1], xyz[2]);
            glRotatef(RAD2DEG*angle, axis[0], axis[1], axis[2]);
            drawobjects_as_vectors = 1;
            glScalef(sensorrelsize*vector_baselength, sensorrelsize*vector_baselength, sensorrelsize*vector_baselength);
            draw_SVOBJECT(devicei->object, state, devicei->prop, 0, arrow_color_float, 1);
            glPopMatrix();
            break;
          default:
            ASSERT(FFALSE);
            break;
          }
        }
        if(velocity_type == VEL_POLAR){
          float vv;
          float xyz1_old[3], xyz2_old[3];
          float xyz1_new[3], xyz2_new[3];
          unsigned char arrow_color_old[4], *arrow_color_new;
          float anglemin, anglemax, rmin, rmax;


          switch(vectortype){
          case VECTOR_PROFILE:
//            cos(-alpha)  -sin(-alpha)  0
//  rot(z) =  sin(-alpha)   cos(-alpha)  0
//                0            0       1

//            1   0        0             1 0  0
//  rot(x) =  0 cos(90) -sin(90)    =    1 0 -1
//            0 sin(90)  cos(90)         0 1  0

//                   cos(alpha)  0  -sin(alpha)
// rot(z)*rot(x) =  -sin(alpha)  0  -cos(alpha)
//                        0      1        0

// rot(z)*rot(x)*(0,0,vv) = (-sin(alpha)*vv,-cos(alpha)*vv,0)

            vv = SCALE2FDS(vel[0]) / max_dev_vel;
            xyz2_new[0] = xyz[0] - sin(angle*DEG2RAD)*vv;
            xyz2_new[1] = xyz[1] - cos(angle*DEG2RAD)*vv;
            xyz2_new[2] = xyz[2];

            xyz1_new[0] = xyz[0];
            xyz1_new[1] = xyz[1];
            xyz1_new[2] = xyz[2];

            arrow_color_new = arrow_color;
            if(i == treei->first){
              xyz1_old[0]=xyz1_new[0];
              xyz1_old[1]=xyz1_new[1];
              xyz1_old[2]=xyz1_new[2];

              xyz2_old[0]=xyz2_new[0];
              xyz2_old[1]=xyz2_new[1];
              xyz2_old[2]=xyz2_new[2];

              arrow_color_old[0] = arrow_color_new[0];
              arrow_color_old[1] = arrow_color_new[1];
              arrow_color_old[2] = arrow_color_new[2];
              arrow_color_old[3] = arrow_color_new[3];
              continue;
            }

            //  draw triangles for following rectangle

            //   xyz1_new---------xyz2_new
            //      |       /       |
            //   xyz1_old---------xyz2_old
            glBegin(GL_TRIANGLES);
            glColor3ubv(arrow_color_old);
            glVertex3fv(xyz1_old);
            glVertex3fv(xyz2_old);
            glColor3ubv(arrow_color_new);
            glVertex3fv(xyz2_new);

            glColor3ubv(arrow_color_old);
            glVertex3fv(xyz1_old);
            glColor3ubv(arrow_color_new);
            glVertex3fv(xyz2_new);
            glColor3ubv(arrow_color_old);
            glVertex3fv(xyz2_old);

            glColor3ubv(arrow_color_old);
            glVertex3fv(xyz1_old);
            glColor3ubv(arrow_color_new);
            glVertex3fv(xyz2_new);
            glVertex3fv(xyz1_new);

            glColor3ubv(arrow_color_old);
            glVertex3fv(xyz1_old);
            glColor3ubv(arrow_color_new);
            glVertex3fv(xyz1_new);
            glVertex3fv(xyz2_new);

            glEnd();

            xyz1_old[0]=xyz1_new[0];
            xyz1_old[1]=xyz1_new[1];
            xyz1_old[2]=xyz1_new[2];

            xyz2_old[0]=xyz2_new[0];
            xyz2_old[1]=xyz2_new[1];
            xyz2_old[2]=xyz2_new[2];

            arrow_color_old[0]=arrow_color[0];
            arrow_color_old[1]=arrow_color[1];
            arrow_color_old[2]=arrow_color[2];
            arrow_color_old[3]=arrow_color[3];
            break;
          case VECTOR_LINE:
            vv = SCALE2FDS(vel[0]) / max_dev_vel;
            glPushMatrix();
            glTranslatef(xyz[0], xyz[1], xyz[2]);
            glRotatef(-angle, 0.0, 0.0, 1.0);
            glRotatef(90.0, 1.0, 0.0, 0.0);
            glColor3ubv(arrow_color);
            glBegin(GL_LINES);
            glVertex3f(0.0, 0.0, 0.0);
            glVertex3f(0.0, 0.0, vv);
            glEnd();
            glPopMatrix();
            break;
          case VECTOR_ARROW:
            vv = SCALE2FDS(vel[0]) / max_dev_vel;
            glPushMatrix();
            glTranslatef(xyz[0], xyz[1], xyz[2]);
            glRotatef(-angle, 0.0, 0.0, 1.0);
            glRotatef(90.0, 1.0, 0.0, 0.0);
            glColor3ubv(arrow_color);

            glPushMatrix();
            glScalef(1.0, 1.0, vv*vector_baselength);
            drawdisk(vector_basediameter*xyzmaxdiff / 10.0, 1.0, arrow_color);
            glPopMatrix();
            glTranslatef(0.0, 0.0, vv*vector_baselength);
            drawcone(vector_headdiameter*xyzmaxdiff / 10.0, vector_headlength*xyzmaxdiff / 10.0, arrow_color);
            glPopMatrix();
            break;
          case VECTOR_OBJECT:
            vv = SCALE2FDS(vel[0]) / max_dev_vel;
            glPushMatrix();
            glTranslatef(xyz[0], xyz[1], xyz[2]);
            glRotatef(-angle, 0.0, 0.0, 1.0);
            glRotatef(90.0, 1.0, 0.0, 0.0);
            glPushMatrix();
            glScalef(1.0, 1.0, vv*vector_baselength);
            drawdisk(vector_basediameter*xyzmaxdiff / 10.0, 1.0, arrow_color);
            glPopMatrix();

            anglemin = -dangle*DEG2RAD;
            anglemax = -dangle*DEG2RAD;
            rmin = MAX(vv - dvel, 0.0);
            rmax = vv + dvel;
            drawsphereseg(anglemin, anglemax, rmin, rmax);
            glPopMatrix();
            break;
          default:
            ASSERT(FFALSE);
            break;
          }
        }
      }
    }
    glPopMatrix();
    glDisable(GL_LIGHTING);
  }

  glPushMatrix();
  glPushAttrib(GL_POINT_BIT|GL_LINE_BIT);
  glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
  glTranslatef(-xbar0,-ybar0,-zbar0);
  for(ii=0;ii<ndeviceinfo;ii++){
    devicedata *devicei;
    int tagval;
    int save_use_displaylist;
    propdata *prop;
    int j;
    float dpsi;
    float *xyz;

    devicei = deviceinfo + ii;
    prop=devicei->prop;

    if(devicei->object->visible==0||(devicei->prop!=NULL&&devicei->prop->smv_object->visible==0))continue;
    if(devicei->plane_surface!=NULL)continue;
    if(isZoneFireModel==1&&STRCMP(devicei->object->label,"target")==0&&visSensor==0)continue;
    if(devicei->in_zone_csv==1)continue;
    if(isZoneFireModel==1&&STRCMP(devicei->label,"TIME")==0)continue;
    save_use_displaylist=devicei->object->use_displaylist;
    tagval=ii+1;
    if(select_device==1&&show_mode==SELECTOBJECT){

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

    if(show_device_orientation==1){
      float *xyznorm;

      xyznorm=devicei->xyznorm;
      glPushMatrix();
      glScalef(orientation_scale/5.0,orientation_scale/5.0,orientation_scale/5.0);
      glBegin(GL_LINES);
      glVertex3f(0.0,0.0,0.0);
      glVertex3f(xyznorm[0],xyznorm[1],xyznorm[2]);
      glEnd();
      glPopMatrix();
    }
    dpsi=0.0;
    if((active_smokesensors==1&&show_smokesensors!=SMOKESENSORS_HIDDEN&&STRCMP(devicei->object->label,"smokesensor")==0)||
       STRCMP(devicei->object->label,"thermocouple")==0
      ){
      float *xyznorm;

      xyznorm = devicei->xyznorm;
      xyznorm[0]=world_eyepos[0]-devicei->xyz[0];
      xyznorm[1]=world_eyepos[1]-devicei->xyz[1];
      xyznorm[2]=world_eyepos[2]-devicei->xyz[2];

      GetElevAz(xyznorm,&devicei->dtheta,devicei->rotate_axis, &dpsi);
    }
    {
      float *axis;

      axis = devicei->rotate_axis;
      // the statement below causes problems in objects.svo definitions
      glRotatef(devicei->dtheta,axis[0],axis[1],axis[2]);
      glRotatef(-dpsi,0.0,0.0,1.0);
    }
    if(sensorrelsize!=1.0){
      glScalef(sensorrelsize,sensorrelsize,sensorrelsize);
    }
    prop=devicei->prop;
    if(prop!=NULL){
      prop->rotate_axis=devicei->rotate_axis;
      prop->rotate_angle=devicei->dtheta;
    }
    if(devicei->nparams>0&&prop!=NULL){
      prop->nvars_indep=devicei->nparams;
      if(prop->fvals==NULL){
        prop->nvars_indep=devicei->nparams;
        NewMemory((void **)&prop->fvals,prop->nvars_indep*sizeof(float));
      }
      if(prop->vars_indep_index==NULL){
        prop->nvars_indep=devicei->nparams;
        NewMemory((void **)&prop->vars_indep_index,prop->nvars_indep*sizeof(int));
      }
      for(j=0;j<devicei->nparams;j++){
        prop->fvals[j]=devicei->params[j];
        prop->vars_indep_index[j]=j;
      }
    }
    if(showtime==1&&itimes>=0&&itimes<nglobal_times&&showdeviceval==1&&ndevicetypes>0){
      int type,vistype=0;

      type=devicei->type2;
      if(type>=0&&type<ndevicetypes)vistype=devicetypes[type]->type2vis;
      if(vistype==1){
        Output_Device_Val(devicei);
      }
    }
    if(drawobjects_as_vectors==0){
      if(showtime==1&&itimes>=0&&itimes<nglobal_times){
        int state;

        if(devicei->showstatelist==NULL){
          state=devicei->state0;
        }
        else{
          state=devicei->showstatelist[itimes];
        }
        if(colordeviceval==1){
          int type,vistype=0;

          type=devicei->type2;
          if(type>=0&&type<ndevicetypes)vistype=devicetypes[type]->type2vis;
          if(vistype==1)draw_SVOBJECT(devicei->object,state,prop,0,NULL,0);
        }
        else{
          draw_SVOBJECT(devicei->object,state,prop,0,NULL,0);
        }
      }
      else{
        draw_SVOBJECT(devicei->object,devicei->state0,prop,0,NULL,0);
      }
    }
    if(devicei->nparams>0&&prop!=NULL){
      prop->nvars_indep=0;
    }
    devicei->object->use_displaylist=save_use_displaylist;
    glPopMatrix();
  }

  glPopAttrib();
  glPopMatrix();
  drawTargetNorm();
}

/* ----------------------- drawTargetNorm ----------------------------- */

void drawTargetNorm(void){
  int i;
  devicedata *devicei;
  float *xyz, *xyznorm;

  if(isZoneFireModel==1&&hasSensorNorm==1&&visSensor==1&&visSensorNorm==1){
    glPushMatrix();
    glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
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

/* ----------------------- draw_SVOBJECT ----------------------------- */

void draw_SVOBJECT(sv_object *object_dev, int iframe_local, propdata *prop, int recurse_level,float *valrgb,int vis_override){
  sv_object_frame *framei,*frame0;
  tokendata *toknext;
  unsigned char *rgbptr_local;
  unsigned char rgbcolor[4];
  int displaylist_id=0;
  int ii;
  sv_object *object;
  int use_material;

  if(prop!=NULL){
    object=prop->smv_object;
  }
  else{
    object=object_dev;
  }
  if(object->visible==0&&vis_override==0)return;
  if(object == missing_device&&show_missing_objects == 0)return;
  if(iframe_local>object->nframes-1||iframe_local<0)iframe_local=0;
  framei=object->obj_frames[iframe_local];
  frame0=object->obj_frames[0];

  ASSERT(framei->error==0||framei->error==1);

  if(framei->error==1){
    object=error_device;
    framei=error_device->obj_frames[0];
    prop=NULL;
  }

  rgbcolor[0]=255;
  rgbcolor[1]=0;
  rgbcolor[2]=0;
  rgbcolor[3]=255;
  rgbptr_local=rgbcolor;
  glPushMatrix();

// copy in default values ( :var=value in objects.svo file )

  for(ii=0;ii<framei->ntokens;ii++){
    tokendata *toki;

    toki = framei->tokens+ii;
    if(toki->is_label==1){
      toki->var=toki->default_val;
    }
    if(toki->is_texturefile==1){
      strcpy(toki->string,toki->default_string);
    }
  }

  // copy values

  if(prop!=NULL){
    int i;

    // copy time dependent evac data

    if(prop->draw_evac==1&&frame0->nevac_tokens>0){
      tokendata *tok00;

      tok00 = frame0->tokens;
      for(i=0;i<NEVAC_TOKENS;i++){
        tokendata *toki,*tok0;
        int itok;

        tok0 = frame0->evac_tokens[i];
        if(tok0==NULL)continue;
        itok = tok0 - tok00;
        toki = framei->tokens + itok;
        toki->var=tok0->evac_var;
      }
    }

    // copy time dependent data using variables from the class_of_... lines

    if(prop->nvars_dep>0){
      for(i=0;i<prop->nvars_dep;i++){
        tokendata *toki;
        int index;

        index = prop->vars_dep_index[i];
        if(index<0||index>framei->ntokens-1)continue;
        toki = framei->tokens + index;
        toki->var=prop->fvars_dep[i];
      }
    }

    // copy static data from PROP line

    for(i=0;i<prop->nvars_indep;i++){
      tokendata *toki;
      int index;

      index = prop->vars_indep_index[i];
      if(index<0||index>framei->ntokens-1)continue;
      toki = framei->tokens + index;
      toki->var=prop->fvals[i];
      if(prop->svals!=NULL&&prop->svals[i]!=NULL&&strlen(prop->svals[i])>0){
        strcpy(toki->string,prop->svals[i]);
      }
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

  use_material=0;
  if(select_device_color_ptr==NULL&&recurse_level==0){
    glEnable(GL_LIGHTING);

    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,specular);

    glEnable(GL_COLOR_MATERIAL);
    use_material=1;
  }
  toknext=NULL;
  for(ii=0;;ii++){
    tokendata *toki;
#define NARGVAL 6
    float arg[NARGVAL], *argptr;
    int j;

    if(ii==0){
      toki=framei->command_list[0];
    }
    else{
      toki=toknext;
    }
    if(toki==NULL)break;
    toknext=toki->next;

    if(select_device_color_ptr==NULL){
      rgbptr_local=rgbcolor;
    }
    else{
      rgbptr_local=select_device_color_ptr;
    }
    for(j=0;j<toki->nvars;j++){
      tokendata *tokj;

      tokj = toki - toki->nvars + j;
      arg[j] = *(tokj->varptr);
    }
    if(toki->nvars>0){
      argptr=(toki-1)->varptr;
    }

    switch(toki->command){
    case SV_ADD:
      {
        float val1, val2, val_result;

        val1=arg[0];
        val2=arg[1];


        val_result=val1+val2;

        *argptr=val_result;
      }
      break;
    case SV_ORIENX:
      if(arg[2]<10.0){
        float u[3]={1.0,0.0,0.0}, axis[3], angle;

        rotateu2v(u, arg, axis, &angle);
        glRotatef(RAD2DEG*angle,axis[0],axis[1],axis[2]);

      }
      break;
    case SV_ORIENY:
      if(arg[2]<10.0){
        float u[3]={0.0,1.0,0.0}, axis[3], angle;

        rotateu2v(u, arg, axis, &angle);
        glRotatef(RAD2DEG*angle,axis[0],axis[1],axis[2]);

      }
      break;
    case SV_ORIENZ:
      if(arg[2]<10.0){
        float u[3]={0.0,0.0,1.0}, axis[3], angle;

        rotateu2v(u, arg, axis, &angle);
        glRotatef(RAD2DEG*angle,axis[0],axis[1],axis[2]);

      }
      break;
    case SV_RANDXY:
      if(ABS(arg[0]-1.0)<0.01){
        float random_angle=0.0;

        random_angle=rand_ab(prop->tag_number,0.0,360.0);
        glRotatef(random_angle,0.0,0.0,1.0);
      }
      break;
    case SV_RANDXZ:
      if(ABS(arg[0]-1.0)<0.01){
        float random_angle=0.0;

        random_angle=rand_ab(prop->tag_number,0.0,360.0);
        glRotatef(random_angle,0.0,1.0,0.0);
      }
      break;
    case SV_RANDYZ:
      if(ABS(arg[0]-1.0)<0.01){
        float random_angle=0.0;

        random_angle=rand_ab(prop->tag_number,0.0,360.0);
        glRotatef(random_angle,1.0,0.0,0.0);
      }
      break;
    case SV_RANDXYZ:
      if(ABS(arg[0]-1.0)<0.01){
        float zz, tt, rr, xx, yy, olddir[3], newdir[3], axis[3], angle;

//    Choose z uniformly distributed in [-1,1].
//    Choose t uniformly distributed on [0, 2*pi).
//    Let r = sqrt(1-z^2).
//    Let x = r * cos(t).
//    Let y = r * sin(t).

        zz = rand_ab(2*prop->tag_number-1,-1.0,1.0);
        tt = rand_ab(2*prop->tag_number,0.0,2.0*PI);
        rr = sqrt(ABS(1.0-zz*zz));
        xx = rr*cos(tt);
        yy = rr*sin(tt);
        olddir[0]=1.0;
        olddir[1]=0.0;
        olddir[2]=0.0;
        newdir[0]=xx;
        newdir[1]=yy;
        newdir[2]=zz;
        rotateu2v(olddir, newdir, axis, &angle);
        glRotatef(RAD2DEG*angle,axis[0],axis[1],axis[2]);
      }
      break;
   	case SV_INCLUDE:
	  case SV_INCLUDEF:
	    {
        sv_object *included_object;
        int iframe_local2;
	      char *object_name;

	      if(toki->included_object==NULL){
	        if(toki->command==SV_INCLUDEF){
	          iframe_local2=arg[0];
		      }
	        else{
            iframe_local2=0;
		      }
          object_name = (toki-1)->string;
          included_object = get_SVOBJECT_type(object_name,missing_device);
	        toki->included_frame=iframe_local2;
	        toki->included_object=included_object;
        }
	      else{
          iframe_local2=toki->included_frame;
          included_object = toki->included_object;
        }
        draw_SVOBJECT(included_object, iframe_local2, NULL, recurse_level+1,NULL,0);
      }
      break;
    case SV_ABS:
      if(arg[0]<0.0){
        *argptr=-arg[0];
      }
      else{
        *argptr=arg[0];
      }
      break;
    case SV_SUB:
      {
        float val1, val2, val_result;

        val1=arg[0];
        val2=arg[1];


        val_result=val1-val2;

        *argptr=val_result;
      }
      break;
    case SV_MULT:
      {
        float val1, val2, val_result;

        val1=arg[0];
        val2=arg[1];


        val_result=val1*val2;

        *argptr=val_result;
      }
      break;
    case SV_DIV:
      {
        float val1, val2, val_result;

        val1=arg[0];
        val2=arg[1];


        if(val2==0.0){
          val_result=0.0;
        }
        else{
          val_result=val1/val2;
        }

        *argptr=val_result;
      }
      break;
    case SV_GETT:
      {
        float val_result;
        float time_val=0.0;

        if(nglobal_times>0){
          time_val=global_times[itimes];
        }

        val_result=time_val;

        *argptr=val_result;
      }
      break;
    case SV_MULTIADDT:
      {
        float val1, val2, val_result;
        float time_val=0.0;

        val1=arg[0];
        val2=arg[1];

        if(nglobal_times>0){
          time_val=global_times[itimes];
        }

        val_result=val1*time_val+val2;

        *argptr=val_result;
      }
      break;
    case SV_CLIP:
      {
        float val, valmin, valmax;

        val=arg[0];
        valmin=arg[1];
        valmax=arg[2];

        *argptr=CLAMP(val,valmin,valmax);
      }
      break;
    case SV_CLIPX:
      {
        clipdata objclip,*ci;

        ci=&objclip;
        ci->clip_xmin=arg[0];
        ci->xmin=arg[1];
        ci->clip_xmax=arg[2];
        ci->xmax=arg[3];
        ci->clip_ymin=-1;
        ci->clip_ymax=-1;
        ci->clip_zmin=-1;
        ci->clip_zmax=-1;
        setClipPlanes(ci,CLIP_ON);
      }
      break;
    case SV_CLIPY:
      {
        clipdata objclip,*ci;

        ci=&objclip;
        ci->clip_ymin=arg[0];
        ci->ymin=arg[1];
        ci->clip_ymax=arg[2];
        ci->ymax=arg[3];
        ci->clip_xmin=-1;
        ci->clip_xmax=-1;
        ci->clip_zmin=-1;
        ci->clip_zmax=-1;
        setClipPlanes(ci,CLIP_ON);
      }
      break;
    case SV_CLIPZ:
      {
        clipdata objclip,*ci;

        ci=&objclip;
        ci->clip_zmin=arg[0];
        ci->zmin=arg[1];
        ci->clip_zmax=arg[2];
        ci->zmax=arg[3];
        ci->clip_xmin=-1;
        ci->clip_xmin=-1;
        ci->clip_ymin=-1;
        ci->clip_ymax=-1;
        setClipPlanes(ci,CLIP_ON);
      }
      break;
    case SV_CLIPOFF:
      setClipPlanes(NULL,CLIP_OFF);
      break;
    case SV_MIRRORCLIP:
      {
        float val, valmin, valmax;
        float val_rel, valmax_rel;

        val=arg[0];
        valmin=arg[1];
        valmax=arg[2];

        valmax_rel=valmax-valmin;
        val_rel=fmod(val-valmin,2.0*valmax_rel);
        if(val_rel<0.0)val_rel+=2.0*valmax_rel;
        if(val_rel>valmax_rel)val_rel=2.0*valmax_rel-val_rel;

        *argptr = val_rel + valmin;
      }
      break;
    case SV_PERIODICCLIP:
      {
        float val, valmin, valmax;
        float val_rel, valmax_rel;

        val=arg[0];
        valmin=arg[1];
        valmax=arg[2];

        val_rel=val-valmin;
        valmax_rel=valmax-valmin;

        val_rel=fmod(val_rel,valmax_rel);
        if(val_rel<0.0)val+=valmax_rel;

        *argptr = val_rel + valmin;
      }
      break;
    case SV_GTRANSLATE:
      if(prop!=NULL){
        float *axis;

        axis = prop->rotate_axis;
        glRotatef(-prop->rotate_angle,axis[0],axis[1],axis[2]);
        glTranslatef(arg[0],arg[1],arg[2]);
        glRotatef(prop->rotate_angle,axis[0],axis[1],axis[2]);
      }
      else{
        glTranslatef(arg[0],arg[1],arg[2]);
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
    case SV_IF:
      if(ABS(arg[0])<=0.001){
        toknext=toki->elsenext;
      }
      break;
    case SV_ELSE:
    case SV_ENDIF:
      break;
    case SV_AND:
      if(ABS(arg[0])>=0.001&&ABS(arg[1])>=0.001){
        *argptr=1.0;
      }
      else{
        *argptr=0.0;
      }
      break;
    case SV_OR:
      if(ABS(arg[0])>=0.001||ABS(arg[1])>=0.001){
        *argptr=1.0;
      }
      else{
        *argptr=0.0;
      }
      break;
    case SV_EQ:
      {
        float *to_ptr, *from_ptr;

        to_ptr=(toki-2)->varptr;
        from_ptr=(toki-1)->varptr;
        *to_ptr=*from_ptr;
      }
      break;
    case SV_GT:
      if(arg[0]>arg[1]){
        *argptr=1.0;
      }
      else{
        *argptr=0.0;
      }
      break;
    case SV_GE:
      if(arg[0]>=arg[1]){
        *argptr=1.0;
      }
      else{
        *argptr=0.0;
      }
      break;
    case SV_LT:
      if(arg[0]<arg[1]){
        *argptr=1.0;
      }
      else{
        *argptr=0.0;
      }
      break;
    case SV_LE:
      if(arg[0]<=arg[1]){
        *argptr=1.0;
      }
      else{
        *argptr=0.0;
      }
      break;
    case SV_ROTATEXYZ:
      rotatexyz(arg[0],arg[1],arg[2]);
      break;
    case SV_ROTATEAXIS:
      rotateaxis(arg[0],arg[1],arg[2],arg[3]);
      break;
    case SV_ROTATEEYE:
      rotateeye();
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
    case SV_SCALEAUTO:
      glScalef(SCALE2FDS(arg[0]),SCALE2FDS(arg[0]),SCALE2FDS(arg[0]));
      break;
    case SV_SCALEGRID:
      glScalef(arg[0]*min_gridcell_size,arg[0]*min_gridcell_size,arg[0]*min_gridcell_size);
      break;
    case SV_SCALE:
      glScalef(arg[0],arg[1],arg[2]);
      break;
    case SV_DRAWCUBE:
      drawcube(arg[0],rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWCUBEC:
      drawcubec(arg[0],rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWSQUARE:
      drawsquare(arg[0],rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWVENT:
      drawvent(arg[0],arg[1],rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWDISK:
      drawdisk(arg[0],arg[1], rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWARCDISK:
      drawarcdisk(arg[0],arg[1], arg[2], rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWCDISK:
      drawcdisk(arg[0],arg[1], rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWHEXDISK:
      drawhexdisk(arg[0],arg[1], rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWPOLYDISK:
      {
        int nsides;

        nsides = arg[0]+0.5;
        drawpolydisk(nsides, arg[1],arg[2], rgbptr_local);
        rgbptr_local=NULL;
      }
      break;
    case SV_DRAWRING:
      drawring(arg[0],arg[1], arg[2], rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWNOTCHPLATE:
      drawnotchplate(arg[0],arg[1], arg[2], arg[3], rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWTRUNCCONE:
      drawtrunccone(arg[0],arg[1],arg[2], rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWCONE:
      drawcone(arg[0],arg[1], rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWTSPHERE:
      {
        int texture_index;

        texture_index = arg[0]+0.5;
        drawtsphere(texture_index,arg[1],rgbptr_local);
      }
      rgbptr_local=NULL;
      break;
    case SV_DRAWSPHERE:
      drawsphere(arg[0],rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWHSPHERE:
      drawhsphere(arg[0],rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWTRIBLOCK:
      drawtriblock(arg[0],arg[1],rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWCIRCLE:
      drawcircle(arg[0],rgbptr_local, &object_circ);
      rgbptr_local=NULL;
      break;
    case SV_DRAWFILLEDCIRCLE:
      drawfilledcircle(arg[0],rgbptr_local, &object_circ);
      rgbptr_local=NULL;
      break;
    case SV_DRAWARC:
      drawarc(arg[0],arg[1],rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_DRAWPOINT:
      drawpoint(rgbptr_local);
      rgbptr_local=NULL;
      break;
    case SV_GETTEXTUREINDEX:
      {
        char *texturefile;
        int i;
        int textureindex=0;

        texturefile = (toki-2)->stringptr;

        for(i=0;i<ndevice_texture_list;i++){
          if(strcmp(device_texture_list[i],texturefile)==0){
            textureindex=i;
            break;
          }
        }
        *argptr=textureindex;
      }
      break;
    case SV_SETCOLOR:
    case SV_SETRGB:
    case SV_SETRGBVAL:
      if(toki->command==SV_SETCOLOR){
        FILE_SIZE lenstring;
        int iarg[3];
        char *stringptr;

        stringptr = (toki-1)->string;

        lenstring=(FILE_SIZE)strlen(stringptr);
        FORTcolor2rgb(iarg,stringptr,lenstring);
        arg[0]=iarg[0];
        arg[1]=iarg[1];
        arg[2]=iarg[2];
      }
      if(valrgb!=NULL&&toki->command==SV_SETRGBVAL){
        arg[0]=valrgb[0];
        arg[1]=valrgb[1];
        arg[2]=valrgb[2];
      }
      if(setbw==1){
        float grey;

        grey = color2bw(arg);
        rgbcolor[0]=grey;
        rgbcolor[1]=grey;
        rgbcolor[2]=grey;
        rgbcolor[3]=255;
      }
      else{
        rgbcolor[0]=arg[0];
        rgbcolor[1]=arg[1];
        rgbcolor[2]=arg[2];
        rgbcolor[3]=255;
      }
      if(select_device_color_ptr==NULL){
        rgbptr_local=rgbcolor;
      }
      else{
        rgbptr_local=select_device_color_ptr;
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
          rgbptr_local=rgbcolor;
        }
        else{
          rgbptr_local=select_device_color_ptr;
        }
      }
      break;
    case SV_DRAWLINE:
      drawline(arg,arg+3,rgbptr_local);
      rgbptr_local=NULL;
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
  if(use_material==1&&recurse_level==0){
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

/* ----------------------- drawtsphere ----------------------------- */

void drawtsphere(int texture_index,float diameter, unsigned char *rgbcolor){
  texturedata *texti;

  if(texture_index<0||texture_index>ntextures-1){
    texti=NULL;
  }
  else{
    texti = textureinfo + texture_index;
    if(texti->loaded==0||texti->display==0)texti=NULL;
  }
  if(texti!=NULL&&object_outlines==0){
    int i,j;

    glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
    glEnable(GL_TEXTURE_2D);

    glBindTexture(GL_TEXTURE_2D,texti->name);

    glPushMatrix();
    glScalef(diameter/2.0,diameter/2.0,diameter/2.0);
    if(cos_lat==NULL)Init_Sphere(NLAT,NLONG);
    glBegin(GL_QUADS);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
    for(j=0;j<NLAT;j++){
      float ti,tip1;
      float tj,tjp1;

      tj = 1.0-(float)j/NLAT;
      tjp1 = 1.0-(float)(j+1)/NLAT;
      for(i=0;i<NLONG;i++){
        float x, y, z;

        ti = 1.0-(float)i/(float)NLONG;
        tip1 = 1.0-(float)(i+1)/(float)NLONG;

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
    glDisable(GL_TEXTURE_2D);
  }
  else{
    drawsphere(diameter,rgbcolor);
  }
}

/* ----------------------- drawsphereseg ----------------------------- */

void drawsphereseg(float anglemin, float anglemax, float rmin, float rmax){
  int i, j;
  float ai, aip1, aj, ajp1;
  float danglei,danglej;
  float cosi, cosip1, sini, sinip1;
  float cosj, cosjp1, sinj, sinjp1;
  float colorin[4]={1.0,0.0,0.0,1.0};
  float colorout[4]={0.0,1.0,0.0,1.0};
  float coloredge[4]={0.0,0.0,1.0,1.0};
  float colori[4]={1.0,0.0,0.0,1.0};
  float colorip1[4]={1.0,0.0,0.0,1.0};

  anglemax=0.0;

  danglei = (anglemax-anglemin)/(float)NLAT;
  danglej = 2.0*PI/(float)NLONG;

  if(object_outlines==0){
  glColor4fv(coloredge);
  glBegin(GL_QUADS);
  ai = anglemin;
  cosi = cos(ai);
  sini = sin(ai);
  for(j=0;j<NLONG;j++){
    aj = j*danglej;
    ajp1 = (j+1)*danglej;
    cosj = cos(aj);
    cosjp1 = cos(ajp1);
    sinj = sin(aj);
    sinjp1 = sin(ajp1);
    glNormal3f(-cosj,-sinj,0.0);
    glVertex3f(rmin*sini*cosj,rmin*sini*sinj,cosi*rmin);

    glNormal3f(-cosjp1,-sinjp1,0.0);
    glVertex3f(rmin*sini*cosjp1,rmin*sini*sinjp1,cosi*rmin);

    glNormal3f(-cosjp1,-sinjp1,0.0);
    glVertex3f(rmax*sini*cosjp1,rmax*sini*sinjp1,cosi*rmax);

    glNormal3f(-cosj,-sinj,0.0);
    glVertex3f(rmax*sini*cosj,rmax*sini*sinj,cosi*rmax);
  }

  memcpy(colori,colorin,4*sizeof(float));
  memcpy(colorip1,colorin,4*sizeof(float));
  for(i=0;i<NLAT;i++){
    ai = anglemin + i*danglei;
    aip1 = anglemin + (i+1)*danglei;
    cosi = cos(ai);
    cosip1 = cos(aip1);
    sini = sin(ai);
    sinip1 = sin(aip1);
    colori[1]=0.6*(float)i/(float)NLAT;
    colorip1[1]=0.6*(float)(i+1)/(float)NLAT;
    glColor4fv(colori);
    for(j=0;j<NLONG;j++){
      aj = j*danglej;
      ajp1 = (j+1)*danglej;
      cosj = cos(aj);
      cosjp1 = cos(ajp1);
      sinj = sin(aj);
      sinjp1 = sin(ajp1);
      glColor4fv(colori);
      glNormal3f(-sini*cosj,-sini*sinj,-cosi);
      glVertex3f(rmin*sini*cosj,rmin*sini*sinj,cosi*rmin);

      glColor4fv(colorip1);
      glNormal3f(-sinip1*cosj,-sinip1*sinj,-cosip1);
      glVertex3f(rmin*sinip1*cosj,rmin*sinip1*sinj,cosip1*rmin);

      glColor4fv(colorip1);
      glNormal3f(-sinip1*cosjp1,-sinip1*sinjp1,-cosip1);
      glVertex3f(rmin*sinip1*cosjp1,rmin*sinip1*sinjp1,cosip1*rmin);

      glColor4fv(colori);
      glNormal3f(-sini*cosjp1,-sini*sinjp1,-cosi);
      glVertex3f(rmin*sini*cosjp1,rmin*sini*sinjp1,cosi*rmin);
    }
  }

  memcpy(colori,colorout,4*sizeof(float));
  memcpy(colorip1,colorout,4*sizeof(float));
  for(i=0;i<NLAT;i++){
    ai = anglemin + i*danglei;
    aip1 = anglemin + (i+1)*danglei;
    cosi = cos(ai);
    cosip1 = cos(aip1);
    sini = sin(ai);
    sinip1 = sin(aip1);
    colori[2]=0.6*(float)i/(float)NLAT;
    colorip1[2]=0.6*(float)(i+1)/(float)NLAT;
    for(j=0;j<NLONG;j++){
      aj = j*danglej;
      ajp1 = (j+1)*danglej;
      cosj = cos(aj);
      cosjp1 = cos(ajp1);
      sinj = sin(aj);
      sinjp1 = sin(ajp1);

      glColor4fv(colori);
      glNormal3f(sini*cosj,sini*sinj,cosi);
      glVertex3f(rmax*sini*cosj,rmax*sini*sinj,cosi*rmax);

      glNormal3f(sini*cosjp1,sini*sinjp1,cosi);
      glVertex3f(rmax*sini*cosjp1,rmax*sini*sinjp1,cosi*rmax);

      glColor4fv(colorip1);
      glNormal3f(sinip1*cosjp1,sinip1*sinjp1,cosip1);
      glVertex3f(rmax*sinip1*cosjp1,rmax*sinip1*sinjp1,cosip1*rmax);

      glNormal3f(sinip1*cosj,sinip1*sinj,cosip1);
      glVertex3f(rmax*sinip1*cosj,rmax*sinip1*sinj,cosip1*rmax);
    }
  }
  glEnd();
  }
  else{
  ai = anglemin;
  glColor4fv(coloredge);
  glBegin(GL_LINES);
  cosi = cos(ai);
  sini = sin(ai);
  for(j=0;j<NLONG;j++){
    aj = j*danglej;
    ajp1 = (j+1)*danglej;
    cosj = cos(aj);
    cosjp1 = cos(ajp1);
    sinj = sin(aj);
    sinjp1 = sin(ajp1);
    glVertex3f(rmin*sini*cosj,rmin*sini*sinj,cosi*rmin);

    glVertex3f(rmax*sini*cosj,rmax*sini*sinj,cosi*rmax);
  }

  memcpy(colori,colorin,4*sizeof(float));
  memcpy(colorip1,colorin,4*sizeof(float));
  for(i=0;i<NLAT;i++){
    ai = anglemin + i*danglei;
    aip1 = anglemin + (i+1)*danglei;
    cosi = cos(ai);
    cosip1 = cos(aip1);
    sini = sin(ai);
    sinip1 = sin(aip1);
    colori[1]=0.6*(float)i/(float)NLAT;
    colorip1[1]=0.6*(float)(i+1)/(float)NLAT;
    glColor4fv(colori);
    for(j=0;j<NLONG;j++){
      aj = j*danglej;
      ajp1 = (j+1)*danglej;
      cosj = cos(aj);
      cosjp1 = cos(ajp1);
      sinj = sin(aj);
      sinjp1 = sin(ajp1);
      glColor4fv(colori);
      glVertex3f(rmin*sini*cosj,rmin*sini*sinj,cosi*rmin);

      glColor4fv(colorip1);
      glVertex3f(rmin*sinip1*cosj,rmin*sinip1*sinj,cosip1*rmin);

      glColor4fv(colorip1);
      glVertex3f(rmin*sinip1*cosjp1,rmin*sinip1*sinjp1,cosip1*rmin);

      glColor4fv(colori);
      glVertex3f(rmin*sini*cosjp1,rmin*sini*sinjp1,cosi*rmin);
    }
  }

  memcpy(colori,colorout,4*sizeof(float));
  memcpy(colorip1,colorout,4*sizeof(float));
  for(i=0;i<NLAT;i++){
    ai = anglemin + i*danglei;
    aip1 = anglemin + (i+1)*danglei;
    cosi = cos(ai);
    cosip1 = cos(aip1);
    sini = sin(ai);
    sinip1 = sin(aip1);
    colori[2]=0.6*(float)i/(float)NLAT;
    colorip1[2]=0.6*(float)(i+1)/(float)NLAT;
    for(j=0;j<NLONG;j++){
      aj = j*danglej;
      ajp1 = (j+1)*danglej;
      cosj = cos(aj);
      cosjp1 = cos(ajp1);
      sinj = sin(aj);
      sinjp1 = sin(ajp1);

      glColor4fv(colori);
      glVertex3f(rmax*sini*cosj,rmax*sini*sinj,cosi*rmax);

      glVertex3f(rmax*sini*cosjp1,rmax*sini*sinjp1,cosi*rmax);

      glColor4fv(colorip1);
      glVertex3f(rmax*sinip1*cosjp1,rmax*sinip1*sinjp1,cosip1*rmax);

      glVertex3f(rmax*sinip1*cosj,rmax*sinip1*sinj,cosip1*rmax);
    }
  }
  glEnd();
  }
}

/* ----------------------- drawsphere ----------------------------- */

void drawsphere(float diameter, unsigned char *rgbcolor){
  int i,j;

  if(cos_lat==NULL)Init_Sphere(NLAT,NLONG);

  glPushMatrix();
  glScalef(diameter/2.0,diameter/2.0,diameter/2.0);

  if(object_outlines==0){
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
  }
  else{
    glBegin(GL_LINE_LOOP);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
    for(j=0;j<NLAT;j++){
      for(i=0;i<NLONG;i++){
        float x, y, z;

        x = cos_long[i]*cos_lat[j];
        y = sin_long[i]*cos_lat[j];
        z = sin_lat[j];
        glVertex3f(x,y,z);

        x = cos_long[i+1]*cos_lat[j];
        y = sin_long[i+1]*cos_lat[j];
        z = sin_lat[j];
        glVertex3f(x,y,z);

        x = cos_long[i+1]*cos_lat[j+1];
        y = sin_long[i+1]*cos_lat[j+1];
        z = sin_lat[j+1];
        glVertex3f(x,y,z);

        x = cos_long[i]*cos_lat[j+1];
        y = sin_long[i]*cos_lat[j+1];
        z = sin_lat[j+1];
        glVertex3f(x,y,z);
      }
    }
    glEnd();
  }
  glPopMatrix();
}

/* ----------------------- drawhsphere ----------------------------- */

void drawhsphere(float diameter, unsigned char *rgbcolor){
  int i,j;

  if(cos_lat==NULL)Init_Sphere(NLAT,NLONG);

  glPushMatrix();
  glScalef(diameter/2.0,diameter/2.0,diameter/2.0);

  if(object_outlines==0){
    glBegin(GL_QUADS);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
    for(j=NLAT/2;j<NLAT;j++){
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
    glBegin(GL_TRIANGLES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
    for(i=0;i<NLONG;i++){
      float x, y, z;

      x = cos_long[i+1];
      y = sin_long[i+1];
      z = 0.0;

      glNormal3f(0.0,0.0,-1.0);
      glVertex3f(x,y,z);

      x = cos_long[i];
      y = sin_long[i];
      z = 0.0;
      glVertex3f(x,y,z);

      x = 0.0;
      y = 0.0;
      z = 0.0;
      glVertex3f(x,y,z);
    }
    glEnd();
  }
  else{
    glBegin(GL_LINE_LOOP);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
    for(j=NLAT/2;j<NLAT;j++){
      for(i=0;i<NLONG;i++){
        float x, y, z;

        x = cos_long[i]*cos_lat[j];
        y = sin_long[i]*cos_lat[j];
        z = sin_lat[j];
        glVertex3f(x,y,z);

        x = cos_long[i+1]*cos_lat[j];
        y = sin_long[i+1]*cos_lat[j];
        z = sin_lat[j];
        glVertex3f(x,y,z);

        x = cos_long[i+1]*cos_lat[j+1];
        y = sin_long[i+1]*cos_lat[j+1];
        z = sin_lat[j+1];
        glVertex3f(x,y,z);

        x = cos_long[i]*cos_lat[j+1];
        y = sin_long[i]*cos_lat[j+1];
        z = sin_lat[j+1];
        glVertex3f(x,y,z);
      }
    }
    glEnd();
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
    for(i=0;i<NLONG;i++){
      float x, y, z;

      x = cos_long[i];
      y = sin_long[i];
      z = 0.0;
      glVertex3f(x,y,z);

      x = 0.0;
      y = 0.0;
      z = 0.0;
      glVertex3f(x,y,z);
    }
    glEnd();
  }
  glPopMatrix();
}

/* ----------------------- drawpoint ----------------------------- */

void drawpoint(unsigned char *rgbcolor){
  glBegin(GL_POINTS);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
  glVertex3f(0.0,0.0,0.0);
  glEnd();
}


/* ----------------------- drawrectangle ----------------------------- */

void drawrectangle(float width,float height, unsigned char *rgbcolor){
  glBegin(GL_LINE_LOOP);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(width, 0.0, 0.0);
  glVertex3f(width, height, 0.0);
  glVertex3f(0.0, height, 0.0);
  glVertex3f(0.0, 0.0, 0.0);

  glEnd();
}

/* ----------------------- drawfilledrectangle ----------------------------- */

void drawfilledrectangle(float width,float height, unsigned char *rgbcolor){
  if(object_outlines==0){
    glBegin(GL_TRIANGLES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(width, 0.0, 0.0);
    glVertex3f(width, height, 0.0);

    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(width, height, 0.0);
    glVertex3f(0.0, height, 0.0);

    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(width, height, 0.0);
    glVertex3f(width, 0.0, 0.0);

    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, height, 0.0);
    glVertex3f(width, height, 0.0);
    glEnd();
  }
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(width, 0.0, 0.0);

    glVertex3f(width, 0.0, 0.0);
    glVertex3f(width, height, 0.0);

    glVertex3f(width, height, 0.0);
    glVertex3f(0.0, 0.0, 0.0);

    glVertex3f(0.0, 0.0, 0.0);
    glVertex3f(0.0, height, 0.0);

    glVertex3f(0.0, height, 0.0);
    glVertex3f(width, height, 0.0);
    glEnd();
  }
}

/* ----------------------- drawfilledcircle ----------------------------- */

void drawfilledcircle(float diameter,unsigned char *rgbcolor, circdata *circinfo){
  int i;
  int ncirc;
  float *xcirc, *ycirc;

  if(circinfo->ncirc==0)Init_Circle(CIRCLE_SEGS,circinfo);
  ncirc = circinfo->ncirc;
  xcirc = circinfo->xcirc;
  ycirc = circinfo->ycirc;

  if(object_outlines==0){
    glBegin(GL_TRIANGLES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
    for(i=0;i<ncirc;i++){
      int i2;

      i2 = (i+1)%ncirc;
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0);
      glVertex3f(0.0,0.0,0.0);
      glVertex3f(diameter*xcirc[  i2]/2.0,diameter*ycirc[  i2]/2.0,0.0);
    }
    for(i=0;i<ncirc;i++){
      int i2;

      i2 = (i+1)%ncirc;
      glVertex3f(diameter*xcirc[  i2]/2.0,diameter*ycirc[  i2]/2.0,0.0);
      glVertex3f(0.0,0.0,0.0);
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0);
    }
    glEnd();
  }
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
    for(i=0;i<ncirc;i++){
      int i2;

      i2 = (i+1)%ncirc;
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0);
      glVertex3f(0.0,0.0,0.0);
      glVertex3f(diameter*xcirc[  i2]/2.0,diameter*ycirc[  i2]/2.0,0.0);
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0);
    }
    glEnd();
  }
}

/* ----------------------- drawcircle ----------------------------- */

void drawcircle(float diameter,unsigned char *rgbcolor, circdata *circinfo){
  int i;
  int ncirc;
  float *xcirc, *ycirc;

  if(circinfo->ncirc==0)Init_Circle(CIRCLE_SEGS,circinfo);
  ncirc = circinfo->ncirc;
  xcirc = circinfo->xcirc;
  ycirc = circinfo->ycirc;

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
  float *xcirc, *ycirc;

  if(object_circ.ncirc==0)Init_Circle(CIRCLE_SEGS,&object_circ);
  xcirc = object_circ.xcirc;
  ycirc = object_circ.ycirc;

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

  if(object_outlines==0){
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
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);


    glVertex3f(-s2,-s2,-s2);
    glVertex3f(-s2,-s2, s2);
    glVertex3f( s2,-s2,-s2);
    glVertex3f( s2,-s2, s2);
    glVertex3f(-s2, s2,-s2);
    glVertex3f(-s2, s2, s2);
    glVertex3f( s2, s2,-s2);
    glVertex3f( s2, s2, s2);

    glVertex3f(-s2,-s2,-s2);
    glVertex3f(-s2, s2,-s2);
    glVertex3f( s2,-s2,-s2);
    glVertex3f( s2, s2,-s2);
    glVertex3f(-s2,-s2, s2);
    glVertex3f(-s2, s2, s2);
    glVertex3f( s2,-s2, s2);
    glVertex3f( s2, s2, s2);

    glVertex3f(-s2,-s2,-s2);
    glVertex3f( s2,-s2,-s2);
    glVertex3f(-s2,-s2, s2);
    glVertex3f( s2,-s2, s2);
    glVertex3f(-s2, s2,-s2);
    glVertex3f( s2, s2,-s2);
    glVertex3f(-s2, s2, s2);
    glVertex3f( s2, s2, s2);
    glEnd();
  }

}

/* ----------------------- drawcubec_outline ----------------------------- */

void drawcubec_outline(float size, unsigned char *rgbcolor){
  float s1,s2;

  s2 = size;
  s1 = 0.0;

    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    glVertex3f(s1,s1,s1);
    glVertex3f(s1,s1,s2);
    glVertex3f(s1,s2,s1);
    glVertex3f(s1,s2,s2);
    glVertex3f(s2,s1,s1);
    glVertex3f(s2,s1,s2);
    glVertex3f(s2,s2,s1);
    glVertex3f(s2,s2,s2);

    glVertex3f(s1,s1,s1);
    glVertex3f(s1,s2,s1);
    glVertex3f(s1,s1,s2);
    glVertex3f(s1,s2,s2);
    glVertex3f(s2,s1,s1);
    glVertex3f(s2,s2,s1);
    glVertex3f(s2,s1,s2);
    glVertex3f(s2,s2,s2);

    glVertex3f(s1,s1,s1);
    glVertex3f(s2,s1,s1);
    glVertex3f(s1,s1,s2);
    glVertex3f(s2,s1,s2);
    glVertex3f(s1,s2,s1);
    glVertex3f(s2,s2,s1);
    glVertex3f(s1,s2,s2);
    glVertex3f(s2,s2,s2);
    glEnd();
}

/* ----------------------- drawbox_outline ----------------------------- */

void drawbox_outline(float x1, float x2, float y1, float y2, float z1, float z2, float *rgbcolor){

  glBegin(GL_LINES);
  if(rgbcolor!=NULL)glColor3fv(rgbcolor);

  glVertex3f(x1,y1,z1);
  glVertex3f(x1,y1,z2);
  glVertex3f(x2,y1,z1);
  glVertex3f(x2,y1,z2);
  glVertex3f(x2,y2,z1);
  glVertex3f(x2,y2,z2);
  glVertex3f(x1,y2,z1);
  glVertex3f(x1,y2,z2);

  glVertex3f(x1,y1,z1);
  glVertex3f(x1,y2,z1);
  glVertex3f(x2,y1,z1);
  glVertex3f(x2,y2,z1);
  glVertex3f(x2,y1,z2);
  glVertex3f(x2,y2,z2);
  glVertex3f(x1,y1,z2);
  glVertex3f(x1,y2,z2);

  glVertex3f(x1,y1,z1);
  glVertex3f(x2,y1,z1);
  glVertex3f(x1,y2,z1);
  glVertex3f(x2,y2,z1);
  glVertex3f(x1,y2,z2);
  glVertex3f(x2,y2,z2);
  glVertex3f(x1,y1,z2);
  glVertex3f(x2,y1,z2);
  glEnd();
}

/* ----------------------- drawcubec ----------------------------- */

void drawcubec(float size, unsigned char *rgbcolor){
  float s1,s2;

  s2 = size;
  s1 = 0.0;

  if(object_outlines==0){
    glBegin(GL_QUADS);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    glNormal3f(0.0,0.0,-1.0);
    glVertex3f(s1,s1,s1);  // 1
    glVertex3f(s1, s2,s1);  // 4
    glVertex3f( s2, s2,s1);  // 3
    glVertex3f( s2,s1,s1);  // 2

    glNormal3f(0.0,0.0,1.0);
    glVertex3f(s1,s1, s2);  // 5
    glVertex3f( s2,s1, s2);  // 6
    glVertex3f( s2, s2, s2);  // 7
    glVertex3f(s1, s2, s2);  // 8

    glNormal3f(0.0,-1.0,0.0);
    glVertex3f(s1,s1,s1);  // 1
    glVertex3f( s2,s1,s1);  // 2
    glVertex3f( s2,s1, s2);  // 6
    glVertex3f(s1,s1, s2);  // 5

    glNormal3f(0.0,1.0,0.0);
    glVertex3f( s2, s2,s1);  // 3
    glVertex3f(s1, s2,s1);  // 4
    glVertex3f(s1, s2, s2);  // 8
    glVertex3f( s2, s2, s2);  // 7

    glNormal3f(-1.0,0.0,0.0);
    glVertex3f(s1,s1,s1);  // 1
    glVertex3f(s1,s1, s2);  // 5
    glVertex3f(s1, s2, s2);  // 8
    glVertex3f(s1, s2,s1);  // 4

    glNormal3f(1.0,0.0,0.0);
    glVertex3f( s2,s1,s1);  // 2
    glVertex3f( s2, s2,s1);  // 3
    glVertex3f( s2, s2, s2);  // 7
    glVertex3f( s2,s1, s2);  // 6
    glEnd();
  }
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    glVertex3f(s1,s1,s1);
    glVertex3f(s1,s1,s2);
    glVertex3f(s1,s2,s1);
    glVertex3f(s1,s2,s2);
    glVertex3f(s2,s1,s1);
    glVertex3f(s2,s1,s2);
    glVertex3f(s2,s2,s1);
    glVertex3f(s2,s2,s2);

    glVertex3f(s1,s1,s1);
    glVertex3f(s1,s2,s1);
    glVertex3f(s1,s1,s2);
    glVertex3f(s1,s2,s2);
    glVertex3f(s2,s1,s1);
    glVertex3f(s2,s2,s1);
    glVertex3f(s2,s1,s2);
    glVertex3f(s2,s2,s2);

    glVertex3f(s1,s1,s1);
    glVertex3f(s2,s1,s1);
    glVertex3f(s1,s1,s2);
    glVertex3f(s2,s1,s2);
    glVertex3f(s1,s2,s1);
    glVertex3f(s2,s2,s1);
    glVertex3f(s1,s2,s2);
    glVertex3f(s2,s2,s2);
    glEnd();
  }

}

/* ----------------------- drawtriblock ----------------------------- */

void drawtriblock(float s, float h, unsigned char *rgbcolor){
  float sd2;
  float ny=0.0, nz=1.0, denom;

  sd2 = s/2.0;
  if(object_outlines==0){
    denom = sqrt(h*h+sd2*sd2);
    if(denom>0.0){
      ny = h/denom;
      nz = sd2/denom;
    }
    glBegin(GL_TRIANGLES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    glNormal3f( 0.0, 0.0,-1.0);
    glVertex3f( 0.0, 0.0, 0.0);  // 1
    glVertex3f( 0.0,   s, 0.0);  // 4
    glVertex3f(   s,   s, 0.0);  // 3
    glVertex3f( 0.0, 0.0, 0.0);  // 1
    glVertex3f(   s,   s, 0.0);  // 3
    glVertex3f(   s, 0.0, 0.0);  // 2

    glNormal3f(-1.0, 0.0, 0.0);
    glVertex3f( 0.0, 0.0, 0.0);  // 1
    glVertex3f( 0.0, sd2,   h);  // 6
    glVertex3f( 0.0,   s, 0.0);  // 4

    glNormal3f( 1.0, 0.0, 0.0);
    glVertex3f(   s, 0.0, 0.0);  // 2
    glVertex3f(   s,   s, 0.0);  // 3
    glVertex3f(   s, sd2,   h);  // 5

    glNormal3f( 0.0, -ny,  nz);
    glVertex3f( 0.0, 0.0, 0.0);  // 1
    glVertex3f(   s, sd2,   h);  // 5
    glVertex3f( 0.0, sd2,   h);  // 6
    glVertex3f( 0.0, 0.0, 0.0);  // 1
    glVertex3f(   s, 0.0, 0.0);  // 2
    glVertex3f(   s, sd2,   h);  // 5

    glNormal3f( 0.0,  ny,  nz);
    glVertex3f(   s,   s, 0.0);  // 3
    glVertex3f( 0.0, sd2,   h);  // 6
    glVertex3f(   s, sd2,   h);  // 5
    glVertex3f(   s,   s, 0.0);  // 3
    glVertex3f( 0.0,   s, 0.0);  // 4
    glVertex3f( 0.0, sd2,   h);  // 6

    glEnd();
  }
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    glVertex3f( 0.0, 0.0, 0.0);  // 1
    glVertex3f(   s, 0.0, 0.0);  // 2
    glVertex3f(   s, 0.0, 0.0);  // 2
    glVertex3f(   s,   s, 0.0);  // 3
    glVertex3f(   s,   s, 0.0);  // 3
    glVertex3f( 0.0,   s, 0.0);  // 4
    glVertex3f( 0.0,   s, 0.0);  // 4
    glVertex3f( 0.0, 0.0, 0.0);  // 1

    glVertex3f( 0.0, 0.0, 0.0);  // 1
    glVertex3f( 0.0, sd2,   h);  // 6
    glVertex3f( 0.0, sd2,   h);  // 6
    glVertex3f( 0.0,   s, 0.0);  // 4

    glVertex3f(   s, 0.0, 0.0);  // 2
    glVertex3f(   s, sd2,   h);  // 5
    glVertex3f(   s, sd2,   h);  // 5
    glVertex3f(   s,   s, 0.0);  // 3

    glVertex3f( 0.0, sd2,   h);  // 6
    glVertex3f(   s, sd2,   h);  // 5

    glEnd();
  }
}

/* ----------------------- drawvent ----------------------------- */

void drawvent(float width, float height, unsigned char *rgbcolor){
  float wd2, hd2, dw, dh;
  int i;
  float dslot;
#define NSLOTS 10
#define FACTOR 2

  dslot = height/(NSLOTS+NSLOTS-1+FACTOR+FACTOR);
  wd2 = width/2.0;
  hd2 = height/2.0;
  dw = dslot*FACTOR;
  dh = dw;

  if(object_outlines==0){
    glBegin(GL_QUADS);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    glNormal3f(0.0,0.0,1.0);

    glVertex3f(-wd2,   -hd2,0.0);
    glVertex3f(-wd2+dw,-hd2,0.0);
    glVertex3f(-wd2+dw, hd2,0.0);
    glVertex3f(-wd2,    hd2,0.0);

    glVertex3f(wd2-dw,-hd2,0.0);
    glVertex3f(wd2,   -hd2,0.0);
    glVertex3f(wd2,    hd2,0.0);
    glVertex3f(wd2-dw, hd2,0.0);

    glVertex3f(-wd2+dw,-hd2,   0.0);
    glVertex3f( wd2-dw,-hd2,   0.0);
    glVertex3f( wd2-dw,-hd2+dh,0.0);
    glVertex3f(-wd2+dw,-hd2+dh,0.0);

    glVertex3f(-wd2+dw,hd2-dh,   0.0);
    glVertex3f( wd2-dw,hd2-dh,   0.0);
    glVertex3f( wd2-dw,hd2   ,0.0);
    glVertex3f(-wd2+dw,hd2   ,0.0);

    glNormal3f(0.0,0.0,-1.0);
    glVertex3f(-wd2,    hd2,0.0);
    glVertex3f(-wd2+dw, hd2,0.0);
    glVertex3f(-wd2+dw,-hd2,0.0);
    glVertex3f(-wd2,   -hd2,0.0);

    glVertex3f(wd2-dw, hd2,0.0);
    glVertex3f(wd2,    hd2,0.0);
    glVertex3f(wd2,   -hd2,0.0);
    glVertex3f(wd2-dw,-hd2,0.0);

    glVertex3f(-wd2+dw,-hd2+dh,0.0);
    glVertex3f( wd2-dw,-hd2+dh,0.0);
    glVertex3f( wd2-dw,-hd2,   0.0);
    glVertex3f(-wd2+dw,-hd2,   0.0);

    glVertex3f(-wd2+dw,hd2   ,0.0);
    glVertex3f( wd2-dw,hd2   ,0.0);
    glVertex3f( wd2-dw,hd2-dh,   0.0);
    glVertex3f(-wd2+dw,hd2-dh,   0.0);

    glNormal3f(0.0,0.0,1.0);
    for(i=0;i<NSLOTS;i++){
      float yy, yy2;

      yy = -hd2+(2*i+FACTOR+1)*dslot;
      yy2 = yy + dslot;
      glVertex3f(-wd2,yy,0.0);
      glVertex3f( wd2,yy,0.0);
      glVertex3f( wd2,yy2,0.0);
      glVertex3f(-wd2,yy2,0.0);
    }
    glNormal3f(0.0,0.0,-1.0);
    for(i=0;i<NSLOTS;i++){
      float yy, yy2;

      yy = -hd2+(2*i+FACTOR+1)*dslot;
      yy2 = yy + dslot;
      glVertex3f(-wd2,yy2,0.0);
      glVertex3f( wd2,yy2,0.0);
      glVertex3f( wd2,yy,0.0);
      glVertex3f(-wd2,yy,0.0);
    }
    glEnd();
  }
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    glVertex3f(-wd2,   -hd2,0.0);
    glVertex3f(-wd2+dw,-hd2,0.0);

    glVertex3f(-wd2+dw,-hd2,0.0);
    glVertex3f(-wd2+dw, hd2,0.0);

    glVertex3f(-wd2+dw, hd2,0.0);
    glVertex3f(-wd2,    hd2,0.0);

    glVertex3f(-wd2,    hd2,0.0);
    glVertex3f(-wd2,   -hd2,0.0);

    glVertex3f(wd2-dw,-hd2,0.0);
    glVertex3f(wd2,   -hd2,0.0);

    glVertex3f(wd2,   -hd2,0.0);
    glVertex3f(wd2,    hd2,0.0);

    glVertex3f(wd2,    hd2,0.0);
    glVertex3f(wd2-dw, hd2,0.0);

    glVertex3f(wd2-dw, hd2,0.0);
    glVertex3f(wd2-dw,-hd2,0.0);

    glVertex3f(-wd2+dw,-hd2,   0.0);
    glVertex3f( wd2-dw,-hd2,   0.0);

    glVertex3f( wd2-dw,-hd2,   0.0);
    glVertex3f( wd2-dw,-hd2+dh,0.0);

    glVertex3f( wd2-dw,-hd2+dh,0.0);
    glVertex3f(-wd2+dw,-hd2+dh,0.0);

    glVertex3f(-wd2+dw,-hd2+dh,0.0);
    glVertex3f(-wd2+dw,-hd2,   0.0);

    glVertex3f(-wd2+dw,hd2-dh,   0.0);
    glVertex3f( wd2-dw,hd2-dh,   0.0);

    glVertex3f( wd2-dw,hd2-dh,   0.0);
    glVertex3f( wd2-dw,hd2   ,0.0);

    glVertex3f( wd2-dw,hd2   ,0.0);
    glVertex3f(-wd2+dw,hd2   ,0.0);

    glVertex3f(-wd2+dw,hd2   ,0.0);
    glVertex3f(-wd2+dw,hd2-dh,   0.0);

    glVertex3f(-wd2,    hd2,0.0);
    glVertex3f(-wd2+dw, hd2,0.0);

    glVertex3f(-wd2+dw, hd2,0.0);
    glVertex3f(-wd2+dw,-hd2,0.0);

    glVertex3f(-wd2+dw,-hd2,0.0);
    glVertex3f(-wd2,   -hd2,0.0);

    glVertex3f(-wd2,   -hd2,0.0);
    glVertex3f(-wd2,    hd2,0.0);

    glVertex3f(wd2-dw, hd2,0.0);
    glVertex3f(wd2,    hd2,0.0);

    glVertex3f(wd2,    hd2,0.0);
    glVertex3f(wd2,   -hd2,0.0);

    glVertex3f(wd2,   -hd2,0.0);
    glVertex3f(wd2-dw,-hd2,0.0);

    glVertex3f(wd2-dw,-hd2,0.0);
    glVertex3f(wd2-dw, hd2,0.0);

    glVertex3f(-wd2+dw,-hd2+dh,0.0);
    glVertex3f( wd2-dw,-hd2+dh,0.0);

    glVertex3f( wd2-dw,-hd2+dh,0.0);
    glVertex3f( wd2-dw,-hd2,   0.0);

    glVertex3f( wd2-dw,-hd2,   0.0);
    glVertex3f(-wd2+dw,-hd2,   0.0);

    glVertex3f(-wd2+dw,-hd2,   0.0);
    glVertex3f(-wd2+dw,-hd2+dh,0.0);

    glVertex3f(-wd2+dw,hd2   ,0.0);
    glVertex3f( wd2-dw,hd2   ,0.0);

    glVertex3f( wd2-dw,hd2   ,0.0);
    glVertex3f( wd2-dw,hd2-dh,   0.0);

    glVertex3f( wd2-dw,hd2-dh,   0.0);
    glVertex3f(-wd2+dw,hd2-dh,   0.0);

    glVertex3f(-wd2+dw,hd2-dh,   0.0);
    glVertex3f(-wd2+dw,hd2   ,0.0);

    for(i=0;i<NSLOTS-1;i++){
      float yy, yy2;

      yy = -hd2+(2*i+FACTOR+1)*dslot;
      yy2 = yy + dslot;
      glVertex3f(-wd2+dw,yy2,0.0);
      glVertex3f( wd2-dw,yy2,0.0);

      glVertex3f( wd2-dw,yy2,0.0);
      glVertex3f( wd2-dw,yy,0.0);

      glVertex3f( wd2-dw,yy,0.0);
      glVertex3f(-wd2+dw,yy,0.0);

      glVertex3f(-wd2+dw,yy,0.0);
      glVertex3f(-wd2+dw,yy2,0.0);
    }
    glEnd();
  }

}

/* ----------------------- drawsquare ----------------------------- */

void drawsquare(float size, unsigned char *rgbcolor){
  float s2;

  s2 = size/2.0;

  if(object_outlines==0){
    glBegin(GL_QUADS);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    glNormal3f(0.0,0.0,-1.0);
    glVertex3f(-s2,-s2,0.0);  // 1
    glVertex3f(-s2, s2,0.0);  // 4
    glVertex3f( s2, s2,0.0);  // 3
    glVertex3f( s2,-s2,0.0);  // 2

    glNormal3f(0.0,0.0,1.0);
    glVertex3f( s2,-s2,0.0);  // 2
    glVertex3f( s2, s2,0.0);  // 3
    glVertex3f(-s2, s2,0.0);  // 4
    glVertex3f(-s2,-s2,0.0);  // 1
    glEnd();
  }
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    glVertex3f(-s2,-s2,0.0);  // 1
    glVertex3f(-s2, s2,0.0);  // 4

    glVertex3f(-s2, s2,0.0);  // 4
    glVertex3f( s2, s2,0.0);  // 3

    glVertex3f( s2, s2,0.0);  // 3
    glVertex3f( s2,-s2,0.0);  // 2

    glVertex3f( s2,-s2,0.0);  // 2
    glVertex3f(-s2,-s2,0.0);  // 1
    glEnd();
  }

}

/* ----------------------- drawring ----------------------------- */

void drawring(float diam_inner, float diam_outer, float height, unsigned char *rgbcolor){
  int i;
  int ncirc;
  float *xcirc, *ycirc;

  if(object_circ.ncirc==0)Init_Circle(CIRCLE_SEGS,&object_circ);
  ncirc = object_circ.ncirc;
  xcirc = object_circ.xcirc;
  ycirc = object_circ.ycirc;

  if(object_outlines==0){
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
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    for(i=0;i<ncirc;i++){
      glVertex3f(diam_outer*xcirc[  i]/2.0,diam_outer*ycirc[  i]/2.0,0.0); // 1
      glVertex3f(diam_outer*xcirc[i+1]/2.0,diam_outer*ycirc[i+1]/2.0,0.0); // 2

      glVertex3f(diam_outer*xcirc[i+1]/2.0,diam_outer*ycirc[i+1]/2.0,0.0); // 2
      glVertex3f(diam_outer*xcirc[i+1]/2.0,diam_outer*ycirc[i+1]/2.0, height); // 3

      glVertex3f(diam_outer*xcirc[i+1]/2.0,diam_outer*ycirc[i+1]/2.0, height); // 3
      glVertex3f(diam_outer*xcirc[  i]/2.0,diam_outer*ycirc[  i]/2.0, height); // 4

      glVertex3f(diam_outer*xcirc[  i]/2.0,diam_outer*ycirc[  i]/2.0, height); // 4
      glVertex3f(diam_outer*xcirc[  i]/2.0,diam_outer*ycirc[  i]/2.0,0.0); // 1

      glVertex3f(diam_inner*xcirc[  i]/2.0,diam_inner*ycirc[  i]/2.0,0.0); // 1
      glVertex3f(diam_inner*xcirc[  i]/2.0,diam_inner*ycirc[  i]/2.0, height); // 4

      glVertex3f(diam_inner*xcirc[  i]/2.0,diam_inner*ycirc[  i]/2.0, height); // 4
      glVertex3f(diam_inner*xcirc[i+1]/2.0,diam_inner*ycirc[i+1]/2.0, height); // 3

      glVertex3f(diam_inner*xcirc[i+1]/2.0,diam_inner*ycirc[i+1]/2.0, height); // 3
      glVertex3f(diam_inner*xcirc[i+1]/2.0,diam_inner*ycirc[i+1]/2.0,0.0); // 2

      glVertex3f(diam_inner*xcirc[i+1]/2.0,diam_inner*ycirc[i+1]/2.0,0.0); // 2
      glVertex3f(diam_inner*xcirc[  i]/2.0,diam_inner*ycirc[  i]/2.0,0.0); // 1
    }
    for(i=0;i<ncirc;i++){
      glVertex3f(diam_outer*xcirc[  i]/2.0,diam_outer*ycirc[  i]/2.0,height);
      glVertex3f(diam_outer*xcirc[i+1]/2.0,diam_outer*ycirc[i+1]/2.0,height);

      glVertex3f(diam_outer*xcirc[i+1]/2.0,diam_outer*ycirc[i+1]/2.0,height);
      glVertex3f(diam_inner*xcirc[i+1]/2.0,diam_inner*ycirc[i+1]/2.0,height);

      glVertex3f(diam_inner*xcirc[i+1]/2.0,diam_inner*ycirc[i+1]/2.0,height);
      glVertex3f(diam_inner*xcirc[  i]/2.0,diam_inner*ycirc[  i]/2.0,height);

      glVertex3f(diam_inner*xcirc[  i]/2.0,diam_inner*ycirc[  i]/2.0,height);
      glVertex3f(diam_outer*xcirc[  i]/2.0,diam_outer*ycirc[  i]/2.0,height);
    }
    for(i=0;i<ncirc;i++){
      glVertex3f(diam_outer*xcirc[  i]/2.0,diam_outer*ycirc[  i]/2.0,0.0);
      glVertex3f(diam_inner*xcirc[  i]/2.0,diam_inner*ycirc[  i]/2.0,0.0);

      glVertex3f(diam_inner*xcirc[  i]/2.0,diam_inner*ycirc[  i]/2.0,0.0);
      glVertex3f(diam_inner*xcirc[i+1]/2.0,diam_inner*ycirc[i+1]/2.0,0.0);

      glVertex3f(diam_inner*xcirc[i+1]/2.0,diam_inner*ycirc[i+1]/2.0,0.0);
      glVertex3f(diam_outer*xcirc[i+1]/2.0,diam_outer*ycirc[i+1]/2.0,0.0);

      glVertex3f(diam_outer*xcirc[i+1]/2.0,diam_outer*ycirc[i+1]/2.0,0.0);
      glVertex3f(diam_outer*xcirc[  i]/2.0,diam_outer*ycirc[  i]/2.0,0.0);
    }
    glEnd();
  }

}


/* ----------------------- rotateeye ----------------------------- */

void rotateeye(void){
  rotatexyz(partfacedir[0],partfacedir[1],partfacedir[2]);
}

/* ----------------------- rotateaxis ----------------------------- */

void rotateaxis(float angle, float ax, float ay, float az){
  glRotatef(angle,ax,ay,az);
}

/* ----------------------- rotatexyz ----------------------------- */

void rotatexyz(float x, float y, float z){
  float angle;
  float normxy,normxyz;

  normxy=x*x+y*y;
  normxy=sqrt(normxy);
  if(normxy<0.00001)return;
  normxyz=x*x+y*y+z*z;
  normxyz=sqrt(normxyz);
  if(normxyz<0.00001)return;
  angle=RAD2DEG*acos(z/normxyz);
  glRotatef(angle,-y/normxy,x/normxy,0.0);
}

/* ----------------------- drawdisk ----------------------------- */

void drawdisk(float diameter, float height, unsigned char *rgbcolor){
  int i;
  int ncirc;
  float *xcirc, *ycirc;

  if(object_circ.ncirc==0)Init_Circle(CIRCLE_SEGS,&object_circ);
  ncirc = object_circ.ncirc;
  xcirc = object_circ.xcirc;
  ycirc = object_circ.ycirc;

  if(object_outlines==0){
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
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    for(i=0;i<ncirc;i++){
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1
      glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,0.0); // 2

      glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,0.0); // 2
      glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height); // 3

      glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height); // 3
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height); // 4

      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height); // 4
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1
    }
    glEnd();

    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    for(i=0;i<ncirc;i++){
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0);
	  glVertex3f(                    0.0,                    0.0,0.0);

	  glVertex3f(                    0.0,                    0.0,0.0);
	  glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,0.0);

	  glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,0.0);
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0);
    }
    for(i=0;i<ncirc;i++){
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height);
	  glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height);

	  glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height);
	  glVertex3f(                    0.0,                    0.0, height);

	  glVertex3f(                    0.0,                    0.0, height);
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height);
    }
    glEnd();
  }
}

/* ----------------------- drawarcdisk ----------------------------- */

void drawarcdisk(float angle, float diameter, float height, unsigned char *rgbcolor){
  int i, iarc;

  if(cos_lat==NULL)Init_Sphere(NLAT,NLONG);

  iarc = NLONG*angle/360.0 + 0.5;
  if(iarc<2)iarc=2;
  if(iarc>NLONG)iarc=NLONG;

  if(object_outlines==0){
    glBegin(GL_QUADS);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    for(i=0;i<iarc;i++){
      glNormal3f(cos_long[i],sin_long[i],0.0);
      glVertex3f(diameter*cos_long[  i]/2.0,diameter*sin_long[  i]/2.0,0.0); // 1

      glNormal3f(cos_long[i+1],sin_long[i+1],0.0);
      glVertex3f(diameter*cos_long[i+1]/2.0,diameter*sin_long[i+1]/2.0,0.0); // 2

      glNormal3f(cos_long[i+1],sin_long[i+1],0.0);
      glVertex3f(diameter*cos_long[i+1]/2.0,diameter*sin_long[i+1]/2.0, height); // 3

      glNormal3f(cos_long[i],sin_long[i],0.0);
      glVertex3f(diameter*cos_long[  i]/2.0,diameter*sin_long[  i]/2.0, height); // 4
    }

    glNormal3f(0.0,-1.0,0.0);
    glVertex3f(0.0,0.0,0.0);
    glVertex3f(diameter*cos_long[  0]/2.0,diameter*sin_long[  0]/2.0,0.0); // 1
    glVertex3f(diameter*cos_long[  0]/2.0,diameter*sin_long[  0]/2.0,height); // 1
    glVertex3f(0.0,0.0,height);

    glNormal3f(sin_long[iarc-1],-cos_long[iarc-1],0.0);
    glVertex3f(0.0,0.0,height);
    glVertex3f(diameter*cos_long[  iarc]/2.0,diameter*sin_long[  iarc]/2.0,height); // 1
    glVertex3f(diameter*cos_long[  iarc]/2.0,diameter*sin_long[  iarc]/2.0,0.0); // 1
    glVertex3f(0.0,0.0,0.0);
    glEnd();

    glBegin(GL_TRIANGLES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    glNormal3f(0.0,0.0,-1.0);
    for(i=0;i<iarc;i++){
      glVertex3f(diameter*cos_long[  i]/2.0,diameter*sin_long[  i]/2.0,0.0);
      glVertex3f(                    0.0,                    0.0,0.0);
      glVertex3f(diameter*cos_long[i+1]/2.0,diameter*sin_long[i+1]/2.0,0.0);
    }
    glNormal3f(0.0,0.0,1.0);
    for(i=0;i<iarc;i++){
      glVertex3f(diameter*cos_long[  i]/2.0,diameter*sin_long[  i]/2.0, height);
      glVertex3f(diameter*cos_long[i+1]/2.0,diameter*sin_long[i+1]/2.0, height);
      glVertex3f(                    0.0,                    0.0, height);
    }
    glEnd();
  }
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    for(i=0;i<iarc;i++){
      glVertex3f(diameter*cos_long[  i]/2.0,diameter*sin_long[  i]/2.0,0.0); // 1
      glVertex3f(diameter*cos_long[i+1]/2.0,diameter*sin_long[i+1]/2.0,0.0); // 2

      glVertex3f(diameter*cos_long[i+1]/2.0,diameter*sin_long[i+1]/2.0,0.0); // 2
      glVertex3f(diameter*cos_long[i+1]/2.0,diameter*sin_long[i+1]/2.0, height); // 3

      glVertex3f(diameter*cos_long[i+1]/2.0,diameter*sin_long[i+1]/2.0, height); // 3
      glVertex3f(diameter*cos_long[  i]/2.0,diameter*sin_long[  i]/2.0, height); // 4

      glVertex3f(diameter*cos_long[  i]/2.0,diameter*sin_long[  i]/2.0, height); // 4
      glVertex3f(diameter*cos_long[  i]/2.0,diameter*sin_long[  i]/2.0,0.0); // 1
    }

    glVertex3f(0.0,0.0,0.0);
    glVertex3f(diameter*cos_long[  0]/2.0,diameter*sin_long[  0]/2.0,0.0); // 1

    glVertex3f(diameter*cos_long[  0]/2.0,diameter*sin_long[  0]/2.0,0.0); // 1
    glVertex3f(diameter*cos_long[  0]/2.0,diameter*sin_long[  0]/2.0,height); // 1

    glVertex3f(diameter*cos_long[  0]/2.0,diameter*sin_long[  0]/2.0,height); // 1
    glVertex3f(0.0,0.0,height);

    glVertex3f(0.0,0.0,height);
    glVertex3f(0.0,0.0,0.0);

    glVertex3f(0.0,0.0,height);
    glVertex3f(diameter*cos_long[  iarc]/2.0,diameter*sin_long[  iarc]/2.0,height); // 1

    glVertex3f(diameter*cos_long[  iarc]/2.0,diameter*sin_long[  iarc]/2.0,height); // 1
    glVertex3f(diameter*cos_long[  iarc]/2.0,diameter*sin_long[  iarc]/2.0,0.0); // 1

    glVertex3f(diameter*cos_long[  iarc]/2.0,diameter*sin_long[  iarc]/2.0,0.0); // 1
    glVertex3f(0.0,0.0,0.0);

    glVertex3f(0.0,0.0,0.0);
    glVertex3f(0.0,0.0,height);
    glEnd();

    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    for(i=0;i<iarc;i++){
      glVertex3f(diameter*cos_long[  i]/2.0,diameter*sin_long[  i]/2.0,0.0);
      glVertex3f(                    0.0,                    0.0,0.0);

      glVertex3f(                    0.0,                    0.0,0.0);
      glVertex3f(diameter*cos_long[i+1]/2.0,diameter*sin_long[i+1]/2.0,0.0);

      glVertex3f(diameter*cos_long[i+1]/2.0,diameter*sin_long[i+1]/2.0,0.0);
      glVertex3f(diameter*cos_long[  i]/2.0,diameter*sin_long[  i]/2.0,0.0);
    }
    for(i=0;i<iarc;i++){
      glVertex3f(diameter*cos_long[  i]/2.0,diameter*sin_long[  i]/2.0, height);
      glVertex3f(diameter*cos_long[i+1]/2.0,diameter*sin_long[i+1]/2.0, height);

      glVertex3f(diameter*cos_long[i+1]/2.0,diameter*sin_long[i+1]/2.0, height);
      glVertex3f(                    0.0,                    0.0, height);

      glVertex3f(                    0.0,                    0.0, height);
      glVertex3f(diameter*cos_long[  i]/2.0,diameter*sin_long[  i]/2.0, height);
    }
    glEnd();
  }

}

/* ----------------------- drawcdisk ----------------------------- */

void drawcdisk(float diameter, float height, unsigned char *rgbcolor){
  int i;
  int ncirc;
  float *xcirc, *ycirc;

  if(object_circ.ncirc==0)Init_Circle(CIRCLE_SEGS,&object_circ);
  ncirc = object_circ.ncirc;
  xcirc = object_circ.xcirc;
  ycirc = object_circ.ycirc;

  if(object_outlines==0){
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
  }
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    for(i=0;i<ncirc;i++){
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,-height/2.00); // 1
	  glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,-height/2.0); // 2

	  glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,-height/2.0); // 2
	  glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height/2.0); // 3

	  glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height/2.0); // 3
	  glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height/2.0); // 4

	  glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height/2.0); // 4
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,-height/2.00); // 1
    }
    glEnd();
  }

  if(object_outlines==0){
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
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    for(i=0;i<ncirc;i++){
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,-height/2.0);
	  glVertex3f(                    0.0,                    0.0,-height/2.0);

	  glVertex3f(                    0.0,                    0.0,-height/2.0);
	  glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,-height/2.0);

	  glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,-height/2.0);
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,-height/2.0);
    }
    for(i=0;i<ncirc;i++){
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height/2.0);
	  glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height/2.0);

	  glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height/2.0);
	  glVertex3f(                    0.0,                    0.0, height/2.0);

	  glVertex3f(                    0.0,                    0.0, height/2.0);
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height/2.0);
    }
    glEnd();
  }
}

/* ----------------------- drawpolydisk ----------------------------- */

void drawpolydisk(int nsides, float diameter, float height, unsigned char *rgbcolor){
  int i;

  float x[33], y[33], xnorm[32], ynorm[32];
  float radius;
  float factor,factor2;

  if(nsides>32)nsides=32;
  if(nsides<3)nsides=3;

  factor=2.0*PI/nsides;
  factor2 = factor/2.0;

  for(i=0;i<nsides;i++){
    x[i]=cos(i*factor);
    y[i]=sin(i*factor);
    xnorm[i] = cos(factor2+i*factor);
    ynorm[i] = sin(factor2+i*factor);
  }
  x[nsides] = x[0];
  y[nsides] = y[0];

  if(object_outlines==0){
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
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    radius = diameter/2.0;

    for(i=0;i<nsides;i++){
      glVertex3f(radius*x[  i],radius*y[  i],0.0); // 1
      glVertex3f(radius*x[i+1],radius*y[i+1],0.0); // 2

      glVertex3f(radius*x[i+1],radius*y[i+1],0.0); // 2
      glVertex3f(radius*x[i+1],radius*y[i+1], height); // 3

      glVertex3f(radius*x[i+1],radius*y[i+1], height); // 3
      glVertex3f(radius*x[  i],radius*y[  i], height); // 4

      glVertex3f(radius*x[  i],radius*y[  i], height); // 4
      glVertex3f(radius*x[  i],radius*y[  i],0.0); // 1
    }
    glEnd();

    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    for(i=0;i<nsides;i++){
      glVertex3f(radius*x[  i],radius*y[  i],0.0);
      glVertex3f(            0.0,            0.0,0.0);

      glVertex3f(            0.0,            0.0,0.0);
      glVertex3f(radius*x[i+1],radius*y[i+1],0.0);

      glVertex3f(radius*x[i+1],radius*y[i+1],0.0);
      glVertex3f(radius*x[  i],radius*y[  i],0.0);
    }
    for(i=0;i<nsides;i++){
      glVertex3f(radius*x[  i],radius*y[  i], height);
      glVertex3f(radius*x[i+1],radius*y[i+1], height);

      glVertex3f(radius*x[i+1],radius*y[i+1], height);
      glVertex3f(            0.0,            0.0, height);

      glVertex3f(            0.0,            0.0, height);
      glVertex3f(radius*x[  i],radius*y[  i], height);
    }
    glEnd();
  }
}

/* ----------------------- drawhexdisk ----------------------------- */

void drawhexdisk(float diameter, float height, unsigned char *rgbcolor){
  int i;

  float x[7]={0.866,0.0,-0.866,-0.866,0.0 ,0.866,0.866};
  float y[7]={0.5,  1.0, 0.5,  -0.5, -1.0,-0.5,  0.5};
  float xnorm[6]={0.500, -0.500, -1.0,-0.500,  0.500, 1.0};
  float ynorm[6]={0.866,  0.866,  0.0,-0.866, -0.866, 0.0};
  float radius;

  if(object_outlines==0){
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
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    radius = diameter/2.0;

    for(i=0;i<6;i++){
      glVertex3f(radius*x[  i],radius*y[  i],0.0); // 1
      glVertex3f(radius*x[i+1],radius*y[i+1],0.0); // 2

      glVertex3f(radius*x[i+1],radius*y[i+1],0.0); // 2
      glVertex3f(radius*x[i+1],radius*y[i+1], height); // 3

      glVertex3f(radius*x[i+1],radius*y[i+1], height); // 3
      glVertex3f(radius*x[  i],radius*y[  i], height); // 4

      glVertex3f(radius*x[  i],radius*y[  i], height); // 4
      glVertex3f(radius*x[  i],radius*y[  i],0.0); // 1
    }
    glEnd();

    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    for(i=0;i<6;i++){
      glVertex3f(radius*x[  i],radius*y[  i],0.0);
      glVertex3f(            0.0,            0.0,0.0);

      glVertex3f(            0.0,            0.0,0.0);
      glVertex3f(radius*x[i+1],radius*y[i+1],0.0);

      glVertex3f(radius*x[i+1],radius*y[i+1],0.0);
      glVertex3f(radius*x[  i],radius*y[  i],0.0);
    }
    for(i=0;i<6;i++){
      glVertex3f(radius*x[  i],radius*y[  i], height);
      glVertex3f(radius*x[i+1],radius*y[i+1], height);

      glVertex3f(radius*x[i+1],radius*y[i+1], height);
      glVertex3f(            0.0,            0.0, height);

      glVertex3f(            0.0,            0.0, height);
      glVertex3f(radius*x[  i],radius*y[  i], height);
    }
    glEnd();
  }
}

/* ----------------------- drawnotchplate ----------------------------- */

void drawnotchplate(float diameter, float height, float notchheight, float direction, unsigned char *rgbcolor){
  int i;
  float diameter2;

  int ncirc;
  float *xcirc, *ycirc;

  if(object_circ.ncirc==0)Init_Circle(CIRCLE_SEGS,&object_circ);
  ncirc = object_circ.ncirc;
  xcirc = object_circ.xcirc;
  ycirc = object_circ.ycirc;

  diameter2 = diameter + notchheight;

  if(object_outlines==0){
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
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    for(i=0;i<ncirc;i++){
      float xmid, ymid;

      xmid = (xcirc[i]+xcirc[i+1])/2.0;
      ymid = (ycirc[i]+ycirc[i+1])/2.0;

      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1
      glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,0.0); // 2

      glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,0.0); // 2
      glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height); // 3

      glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height); // 3
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height); // 4

      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height); // 4
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1

    // draw notch

      if(direction<0.0){
        glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1
        glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, 0.0-notchheight); // 4

        glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, 0.0-notchheight); // 4
        glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0, 0.0-notchheight); // 3

        glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0, 0.0-notchheight); // 3
        glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0,0.0); // 2

        glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0,0.0); // 2
        glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1
      }
      else{
      // top plate
        glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,height); // 1t
        glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, height); // 4t

        glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, height); // 4t
        glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, height); // 3t

        glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, height); // 3t
        glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0,height); // 2t

        glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0,height); // 2t
        glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,height); // 1t

      // bottom plate

        glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1b
        glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0,0.0); // 2b

        glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0,0.0); // 2b
        glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, 0.0); // 3b

        glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, 0.0); // 3b
        glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, 0.0); // 4b

        glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, 0.0); // 4b
        glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1b

      // front plate

        glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, height); // 4t-1
        glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, 0.0); // 4b-4

        glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, 0.0); // 4b-4
        glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, 0.0); // 3b-3

        glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, 0.0); // 3b-3
        glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, height); // 3t-2

        glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, height); // 3t-2
        glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, height); // 4t-1

      // left plate

        glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,height); // 1t
        glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, height); // 4t

        glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, height); // 4t
        glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, 0.0); // 4b

        glVertex3f(diameter2*xcirc[  i]/2.0,diameter2*ycirc[  i]/2.0, 0.0); // 4b
        glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1b

        glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0); // 1b
        glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,height); // 1t

      // right plate

        glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0,height); // 1t
        glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0,0.0); // 1b

        glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0,0.0); // 1b
        glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, 0.0); // 4b

        glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, 0.0); // 4b
        glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, height); // 4t

        glVertex3f(diameter2*xmid/2.0,diameter2*ymid/2.0, height); // 4t
        glVertex3f(diameter*xmid/2.0,diameter*ymid/2.0,height); // 1t
      }
    }
    glEnd();

    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    for(i=0;i<ncirc;i++){
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0);
      glVertex3f(                    0.0,                    0.0,0.0);

      glVertex3f(                    0.0,                    0.0,0.0);
      glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,0.0);

      glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0,0.0);
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0,0.0);
    }
    for(i=0;i<ncirc;i++){
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height);
      glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height);

      glVertex3f(diameter*xcirc[i+1]/2.0,diameter*ycirc[i+1]/2.0, height);
      glVertex3f(                    0.0,                    0.0, height);

      glVertex3f(                    0.0,                    0.0, height);
      glVertex3f(diameter*xcirc[  i]/2.0,diameter*ycirc[  i]/2.0, height);
    }
    glEnd();
  }
}

/* ----------------------- drawcone ----------------------------- */

void drawcone(float d1, float height, unsigned char *rgbcolor){
  int i;
  float factor, denom, rad;
  float hdr;
  int ncirc;
  float *xcirc, *ycirc;

  if(object_circ.ncirc==0)Init_Circle(CIRCLE_SEGS,&object_circ);
  ncirc = object_circ.ncirc;
  xcirc = object_circ.xcirc;
  ycirc = object_circ.ycirc;
  if(height<=0.0)height=0.0001;

  rad = d1/2.0;
  hdr = height/rad;
  denom = 1.0/sqrt(1.0+hdr*hdr);
  factor = hdr*denom;

  if(object_outlines==0){
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
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    for(i=0;i<ncirc;i++){
      glVertex3f(rad*xcirc[  i],rad*ycirc[  i],0.0); // 1
      glVertex3f(rad*xcirc[i+1],rad*ycirc[i+1],0.0); // 2

      glVertex3f(rad*xcirc[i+1],rad*ycirc[i+1],0.0); // 2
      glVertex3f(0.0,0.0, height); // 3

      glVertex3f(0.0,0.0, height); // 3
      glVertex3f(rad*xcirc[  i],rad*ycirc[  i],0.0); // 1
    }
    for(i=0;i<ncirc;i++){
      glVertex3f(rad*xcirc[  i],rad*ycirc[  i],0.0);
      glVertex3f(                    0.0,                    0.0,0.0);

      glVertex3f(                    0.0,                    0.0,0.0);
      glVertex3f(rad*xcirc[i+1],rad*ycirc[i+1],0.0);

      glVertex3f(rad*xcirc[i+1],rad*ycirc[i+1],0.0);
      glVertex3f(rad*xcirc[  i],rad*ycirc[  i],0.0);
    }
    glEnd();
  }
}

/* ----------------------- drawtrunccone ----------------------------- */

void drawtrunccone(float d1, float d2, float height, unsigned char *rgbcolor){
  int i;
  float dz;
  int ncirc;
  float *xcirc, *ycirc;

  if(object_circ.ncirc==0)Init_Circle(CIRCLE_SEGS,&object_circ);
  ncirc = object_circ.ncirc;
  xcirc = object_circ.xcirc;
  ycirc = object_circ.ycirc;

  if(height<=0.0)height=0.0001;
  dz = -(d2-d1)/height;
  if(object_outlines==0){
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
  else{
    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    for(i=0;i<ncirc;i++){
      glVertex3f(d1*xcirc[  i]/2.0,d1*ycirc[  i]/2.0,0.0); // 1
      glVertex3f(d1*xcirc[i+1]/2.0,d1*ycirc[i+1]/2.0,0.0); // 2

      glVertex3f(d1*xcirc[i+1]/2.0,d1*ycirc[i+1]/2.0,0.0); // 2
      glVertex3f(d2*xcirc[i+1]/2.0,d2*ycirc[i+1]/2.0, height); // 3

      glVertex3f(d2*xcirc[i+1]/2.0,d2*ycirc[i+1]/2.0, height); // 3
      glVertex3f(d2*xcirc[  i]/2.0,d2*ycirc[  i]/2.0, height); // 4

      glVertex3f(d2*xcirc[  i]/2.0,d2*ycirc[  i]/2.0, height); // 4
      glVertex3f(d1*xcirc[  i]/2.0,d1*ycirc[  i]/2.0,0.0); // 1
    }
    glEnd();

    glBegin(GL_LINES);
    if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

    for(i=0;i<ncirc;i++){
      glVertex3f(d1*xcirc[  i]/2.0,d1*ycirc[  i]/2.0,0.0);
      glVertex3f(                    0.0,                    0.0,0.0);

      glVertex3f(                    0.0,                    0.0,0.0);
      glVertex3f(d1*xcirc[i+1]/2.0,d1*ycirc[i+1]/2.0,0.0);

      glVertex3f(d1*xcirc[i+1]/2.0,d1*ycirc[i+1]/2.0,0.0);
      glVertex3f(d1*xcirc[  i]/2.0,d1*ycirc[  i]/2.0,0.0);

    }
    for(i=0;i<ncirc;i++){
      glVertex3f(d2*xcirc[  i]/2.0,d2*ycirc[  i]/2.0, height);
      glVertex3f(d2*xcirc[i+1]/2.0,d2*ycirc[i+1]/2.0, height);

      glVertex3f(d2*xcirc[i+1]/2.0,d2*ycirc[i+1]/2.0, height);
      glVertex3f(                    0.0,                    0.0, height);

      glVertex3f(                    0.0,                    0.0, height);
      glVertex3f(d2*xcirc[  i]/2.0,d2*ycirc[  i]/2.0, height);
    }
    glEnd();
  }
}

/* ----------------------- get_SMVOBJECT_type ----------------------------- */

sv_object *get_SVOBJECT_type(char *olabel,sv_object *default_object){
  int i;
  sv_object *objecti;
  char label[256],*labelptr;

  if(olabel==NULL)return default_object;
  strcpy(label,olabel);
  labelptr=label;
  trim_back(label);
  labelptr = trim_front(label);
  if(strlen(labelptr)==0)return default_object;
  for(i=0;i<nobject_defs;i++){
    objecti = object_defs[i];
    if(STRCMP(labelptr,objecti->label)==0){
      objecti->used=1;
      return objecti;
    }
  }
  return default_object;
}

/* ----------------------- get_SMVOBJECT_type2 ----------------------------- */

sv_object *get_SVOBJECT_type2(char *olabel,sv_object *default_object){
  sv_object *object_start, *objecti;
  char label[256],*labelptr;

  if(olabel==NULL)return default_object;
  strcpy(label,olabel);
  labelptr=label;
  trim_back(label);
  labelptr = trim_front(label);
  if(strlen(labelptr)==0)return default_object;
  object_start = object_def_first.next;
  objecti = object_start;
  for(;objecti->next!=NULL;){
    if(STRCMP(labelptr,objecti->label)==0){
      objecti->used=1;
      return objecti;
    }
    objecti=objecti->next;
  }
  return default_object;
}

/* ----------------------- Init_Circle ----------------------------- */

void Init_Circle(unsigned int npoints, circdata *circinfo){
  float drad;
  int i;
  float *xcirc, *ycirc;
  int ncirc;

  if(circinfo->ncirc!=0)freecircle(circinfo);
  if(npoints<2)return;
  ncirc=npoints;
  NewMemory( (void **)&xcirc,(ncirc+1)*sizeof(float));
  NewMemory( (void **)&ycirc,(ncirc+1)*sizeof(float));
  drad=2.0*PI/(float)ncirc;


  for(i=0;i<ncirc;i++){
    xcirc[i] = cos(i*drad);
    ycirc[i] = sin(i*drad);
  }
  xcirc[ncirc]=xcirc[0];
  ycirc[ncirc]=ycirc[0];

  circinfo->xcirc=xcirc;
  circinfo->ycirc=ycirc;
  circinfo->ncirc=npoints;
}

/* ----------------------- Init_Sphere ----------------------------- */

void Init_Sphere(int nlat, int nlong){
  float dlat, dlong;
  int i;

  FREEMEMORY(cos_lat);
  FREEMEMORY(sin_lat);
  FREEMEMORY(cos_long);
  FREEMEMORY(sin_long);
  NewMemory( (void **)&cos_lat,(nlat+1)*sizeof(float));
  NewMemory( (void **)&sin_lat,(nlat+1)*sizeof(float));
  NewMemory( (void **)&cos_long,(nlong+1)*sizeof(float));
  NewMemory( (void **)&sin_long,(nlong+1)*sizeof(float));

  dlat=PI/(float)nlat;
  for(i=0;i<=nlat;i++){
    float angle;

    angle = -PI/2.0 + i*dlat;
    cos_lat[i] = cos(angle);
    sin_lat[i] = sin(angle);
  }

  dlong=2.0*PI/(float)nlong;
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

void freecircle(circdata *circinfo){
  FREEMEMORY(circinfo->xcirc);
  FREEMEMORY(circinfo->ycirc);
  circinfo->ncirc=0;
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
  framei->ntextures=0;
  framei->nevac_tokens=0;

  return object;
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
      framei->ntextures=0;
      framei->nevac_tokens=0;
    }
    else{
      NewMemory((void **)&framei,sizeof(sv_object_frame));
      object->obj_frames[1]=framei;
      framei->device=object;
      parse_device_frame(commandson, NULL, &eof, framei);
      framei->error=0;
      framei->display_list_ID=-1;
      framei->use_bw=setbw;
      framei->ntextures=0;
      framei->nevac_tokens=0;
    }
  }
  return object;
}

/* ----------------------- get_token_id ----------------------------- */

int get_token_id(char *token, int *opptr, int *num_opptr, int *num_outopptr, int *use_displaylist){

  int op, num_op, num_outop;
  int return_val;

  *use_displaylist=0;

  return_val=0;
  if(STRCMP(token,"translate")==0){
    op=SV_TRANSLATE;
    num_op=SV_TRANSLATE_NUMARGS;
    num_outop=SV_TRANSLATE_NUMOUTARGS;
  }
  else if(STRCMP(token,"gtranslate")==0){
    op=SV_GTRANSLATE;
    num_op=SV_GTRANSLATE_NUMARGS;
    num_outop=SV_GTRANSLATE_NUMOUTARGS;
  }
  else if(STRCMP(token,"no_op")==0){
    op=SV_NO_OP;
    num_op=SV_NO_OP_NUMARGS;
    num_outop=SV_NO_OP_NUMOUTARGS;
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
  else if(STRCMP(token,"rotatexyz")==0){
    op=SV_ROTATEXYZ;
    num_op=SV_ROTATEXYZ_NUMARGS;
    num_outop=SV_ROTATEXYZ_NUMOUTARGS;
  }
  else if(STRCMP(token,"rotateaxis")==0){
    op=SV_ROTATEAXIS;
    num_op=SV_ROTATEAXIS_NUMARGS;
    num_outop=SV_ROTATEAXIS_NUMOUTARGS;
  }
  else if(STRCMP(token,"rotateeye")==0){
    op=SV_ROTATEEYE;
    num_op=SV_ROTATEEYE_NUMARGS;
    num_outop=SV_ROTATEEYE_NUMOUTARGS;
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
  else if(STRCMP(token,"if")==0){
    op=SV_IF;
    num_op=SV_IF_NUMARGS;
    num_outop=SV_IF_NUMOUTARGS;
  }
  else if(STRCMP(token,"else")==0){
    op=SV_ELSE;
    num_op=SV_ELSE_NUMARGS;
    num_outop=SV_ELSE_NUMOUTARGS;
  }
  else if(STRCMP(token,"endif")==0){
    op=SV_ENDIF;
    num_op=SV_ENDIF_NUMARGS;
    num_outop=SV_ENDIF_NUMOUTARGS;
  }
  else if(STRCMP(token,"LT")==0){
    op=SV_LT;
    num_op=SV_LT_NUMARGS;
    num_outop=SV_LT_NUMOUTARGS;
  }
  else if(STRCMP(token,"LE")==0){
    op=SV_LE;
    num_op=SV_LE_NUMARGS;
    num_outop=SV_LE_NUMOUTARGS;
  }
  else if(STRCMP(token,"GT")==0){
    op=SV_GT;
    num_op=SV_GT_NUMARGS;
    num_outop=SV_GT_NUMOUTARGS;
  }
  else if(STRCMP(token,"EQ")==0){
    op=SV_EQ;
    num_op=SV_EQ_NUMARGS;
    num_outop=SV_EQ_NUMOUTARGS;
  }
  else if(STRCMP(token,"GE")==0){
    op=SV_GE;
    num_op=SV_GE_NUMARGS;
    num_outop=SV_GE_NUMOUTARGS;
  }
  else if(STRCMP(token,"AND")==0){
    op=SV_AND;
    num_op=SV_AND_NUMARGS;
    num_outop=SV_AND_NUMOUTARGS;
  }
  else if(STRCMP(token,"OR")==0){
    op=SV_OR;
    num_op=SV_OR_NUMARGS;
    num_outop=SV_OR_NUMOUTARGS;
  }
  else if(STRCMP(token,"scalexyz")==0){
    op=SV_SCALEXYZ;
    num_op=SV_SCALEXYZ_NUMARGS;
    num_outop=SV_SCALEXYZ_NUMOUTARGS;
  }
  else if(STRCMP(token,"scaleauto")==0){
    op=SV_SCALEAUTO;
    num_op=SV_SCALEAUTO_NUMARGS;
    num_outop=SV_SCALEAUTO_NUMOUTARGS;
  }
  else if(STRCMP(token,"scalegrid")==0){
    op=SV_SCALEGRID;
    num_op=SV_SCALEGRID_NUMARGS;
    num_outop=SV_SCALEGRID_NUMOUTARGS;
  }
  else if(STRCMP(token,"scale")==0&&STRCMP(token,"scalexyz")!=0&&STRCMP(token,"scaleauto")!=0&&STRCMP(token,"scalegrid")!=0){
    op=SV_SCALE;
    num_op=SV_SCALE_NUMARGS;
    num_outop=SV_SCALE_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawcube")==0){
    op=SV_DRAWCUBE;
    num_op=SV_DRAWCUBE_NUMARGS;
    num_outop=SV_DRAWCUBE_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawcubec")==0){
    op=SV_DRAWCUBEC;
    num_op=SV_DRAWCUBEC_NUMARGS;
    num_outop=SV_DRAWCUBEC_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawvent")==0){
    op=SV_DRAWVENT;
    num_op=SV_DRAWVENT_NUMARGS;
    num_outop=SV_DRAWVENT_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawsquare")==0){
    op=SV_DRAWSQUARE;
    num_op=SV_DRAWSQUARE_NUMARGS;
    num_outop=SV_DRAWSQUARE_NUMOUTARGS;
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
  else if(STRCMP(token,"drawhsphere")==0){
    op=SV_DRAWHSPHERE;
    num_op=SV_DRAWHSPHERE_NUMARGS;
    num_outop=SV_DRAWHSPHERE_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawtriblock")==0){
    op=SV_DRAWTRIBLOCK;
    num_op=SV_DRAWTRIBLOCK_NUMARGS;
    num_outop=SV_DRAWTRIBLOCK_NUMOUTARGS;
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
  else if(STRCMP(token,"drawfilledcircle")==0){
    op=SV_DRAWFILLEDCIRCLE;
    num_op=SV_DRAWFILLEDCIRCLE_NUMARGS;
    num_outop=SV_DRAWFILLEDCIRCLE_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawarc")==0){
    op=SV_DRAWARC;
    num_op=SV_DRAWARC_NUMARGS;
    num_outop=SV_DRAWARC_NUMOUTARGS;
  }
  else if(STRCMP(token,"drawarcdisk")==0){
    op=SV_DRAWARCDISK;
    num_op=SV_DRAWARCDISK_NUMARGS;
    num_outop=SV_DRAWARCDISK_NUMOUTARGS;
  }
  else if(STRCMP(token,"setcolor")==0){
    op=SV_SETCOLOR;
    num_op=SV_SETCOLOR_NUMARGS;
    num_outop=SV_SETCOLOR_NUMOUTARGS;
  }
  else if(STRCMP(token,"gettextureindex")==0){
    op=SV_GETTEXTUREINDEX;
    num_op=SV_GETTEXTUREINDEX_NUMARGS;
    num_outop=SV_GETTEXTUREINDEX_NUMOUTARGS;
  }
  else if(STRCMP(token,"setrgb")==0){
    op=SV_SETRGB;
    num_op=SV_SETRGB_NUMARGS;
    num_outop=SV_SETRGB_NUMOUTARGS;
  }
  else if(STRCMP(token,"setrgbval")==0){
    op=SV_SETRGBVAL;
    num_op=SV_SETRGBVAL_NUMARGS;
    num_outop=SV_SETRGBVAL_NUMOUTARGS;
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
  else if(STRCMP(token,"include")==0){
    op=SV_INCLUDE;
    num_op=SV_INCLUDE_NUMARGS;
    num_outop=SV_INCLUDE_NUMOUTARGS;
  }
  else if(STRCMP(token,"orienx")==0){
    op=SV_ORIENX;
    num_op=SV_ORIENX_NUMARGS;
    num_outop=SV_ORIENX_NUMOUTARGS;
  }
  else if(STRCMP(token,"orieny")==0){
    op=SV_ORIENY;
    num_op=SV_ORIENY_NUMARGS;
    num_outop=SV_ORIENY_NUMOUTARGS;
  }
  else if(STRCMP(token,"orienz")==0){
    op=SV_ORIENZ;
    num_op=SV_ORIENZ_NUMARGS;
    num_outop=SV_ORIENZ_NUMOUTARGS;
  }
  else if(STRCMP(token,"randxy")==0){
    op=SV_RANDXY;
    num_op=SV_RANDXY_NUMARGS;
    num_outop=SV_RANDXY_NUMOUTARGS;
  }
  else if(STRCMP(token,"randxz")==0){
    op=SV_RANDXZ;
    num_op=SV_RANDXZ_NUMARGS;
    num_outop=SV_RANDXZ_NUMOUTARGS;
  }
  else if(STRCMP(token,"randyz")==0){
    op=SV_RANDYZ;
    num_op=SV_RANDYZ_NUMARGS;
    num_outop=SV_RANDYZ_NUMOUTARGS;
  }
  else if(STRCMP(token,"randxyz")==0){
    op=SV_RANDXYZ;
    num_op=SV_RANDXYZ_NUMARGS;
    num_outop=SV_RANDXYZ_NUMOUTARGS;
  }
  else if(STRCMP(token,"includef")==0){
    op=SV_INCLUDEF;
    num_op=SV_INCLUDEF_NUMARGS;
    num_outop=SV_INCLUDEF_NUMOUTARGS;
  }
  else if(STRCMP(token,"abs")==0){
    op=SV_ABS;
    num_op=SV_ABS_NUMARGS;
    num_outop=SV_ABS_NUMOUTARGS;
  }
  else if(STRCMP(token,"add")==0){
    op=SV_ADD;
    num_op=SV_ADD_NUMARGS;
    num_outop=SV_ADD_NUMOUTARGS;
  }
  else if(STRCMP(token,"sub")==0){
    op=SV_SUB;
    num_op=SV_SUB_NUMARGS;
    num_outop=SV_SUB_NUMOUTARGS;
  }
  else if(STRCMP(token,"mult")==0){
    op=SV_MULT;
    num_op=SV_MULT_NUMARGS;
    num_outop=SV_MULT_NUMOUTARGS;
  }
  else if(STRCMP(token,"div")==0){
    op=SV_DIV;
    num_op=SV_DIV_NUMARGS;
    num_outop=SV_DIV_NUMOUTARGS;
  }
  else if(STRCMP(token,"gett")==0){
    op=SV_GETT;
    num_op=SV_GETT_NUMARGS;
    num_outop=SV_GETT_NUMOUTARGS;
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
  else if(STRCMP(token,"clipx")==0){
    op=SV_CLIP;
    *use_displaylist=0;
    num_op=SV_CLIP_NUMARGS;
    num_outop=SV_CLIP_NUMOUTARGS;
  }
  else if(STRCMP(token,"clipy")==0){
    op=SV_CLIP;
    *use_displaylist=0;
    num_op=SV_CLIP_NUMARGS;
    num_outop=SV_CLIP_NUMOUTARGS;
  }
  else if(STRCMP(token,"clipz")==0){
    op=SV_CLIP;
    *use_displaylist=0;
    num_op=SV_CLIP_NUMARGS;
    num_outop=SV_CLIP_NUMOUTARGS;
  }
  else if(STRCMP(token,"clipoff")==0){
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
    num_outop=0;
    return_val=1;
  }
  *opptr=op;
  *num_opptr=num_op;
  *num_outopptr=num_outop;
  return return_val;
}

/* ----------------------- get_token_loc ----------------------------- */

int get_token_loc(char *var,sv_object_frame *frame){
  int i;

  for(i=0;i<frame->nsymbols;i++){
    int ii;
    tokendata *toki;
    char *token_var;

    ii = frame->symbols[i];
    toki = frame->tokens+ii;
    token_var = toki->tokenlabel+1;
    if(STRCMP(var,token_var)==0)return ii;
  }
  return -1;
}

/* ----------------------- get_token_ptr ----------------------------- */

tokendata *get_token_ptr(char *var,sv_object_frame *frame){
  int i;

  for(i=0;i<frame->nsymbols;i++){
    int ii;
    tokendata *toki;
    char *token_var;

    ii = frame->symbols[i];
    toki = frame->tokens+ii;
    token_var = toki->tokenlabel+1;
    if(STRCMP(var,token_var)==0)return toki;
  }
  return NULL;
}

/* ----------------------- parse_device_frame ----------------------------- */

char *parse_device_frame(char *buffer, FILE *stream, int *eof, sv_object_frame *frame){
#define BUFFER_SIZE 10000

  char  object_buffer[10*BUFFER_SIZE];
  int ntokens;
  char *token,*tokens[BUFFER_SIZE];
  char *buffer_ptr=NULL,*buffer2;
  int i;
  int nsymbols,ncommands;
  int ntextures_local=0;
  int last_command_index=0;

  *eof = 0;

  frame->error=0;
  trim_back(buffer);
  strcpy(object_buffer,buffer);
  while(stream!=NULL&&!feof(stream)){
    if(fgets(buffer,255,stream)==NULL){
      *eof=1;
      break;
    }
    buffer2=remove_comment(buffer);
    if(match(buffer2,"OBJECTDEF") == 1||
       match(buffer2,"AVATARDEF") == 1||
       match(buffer2,"NEWFRAME") == 1){
         buffer_ptr=buffer2;
         break;
    }
    strcat(object_buffer," ");
    strcat(object_buffer,buffer2);
  }
  parse_object_string(object_buffer,tokens,&ntokens);
  frame->ntokens=ntokens;
  if(ntokens>0){
    NewMemory((void **)&frame->tokens,ntokens*sizeof(tokendata));
    NewMemory((void **)&frame->symbols,ntokens*sizeof(int));
    NewMemory((void **)&frame->command_list,ntokens*sizeof(tokendata *));
  }

  // count symbols and commands, zero out access counter

  nsymbols=0;
  ncommands=0;
  for(i=0;i<ntokens;i++){
    tokendata *toki;
    char c;
    token=tokens[i];
    toki = frame->tokens + i;
    strcpy(toki->token,token);
    strcpy(toki->tokenlabel,token);
    strcpy(toki->tokenfulllabel,token);
    toki->reads=0;

    c = token[0];

    if(c==':'){
      frame->symbols[nsymbols++]=i;
    }
    if(c>='a'&&c<='z'||c>='A'&&c<='Z')ncommands++;
  }
  frame->nsymbols=nsymbols;
  frame->ncommands=ncommands;

  // fill in token data structure

  nsymbols=0;
  ncommands=0;
  for(i=0;i<ntokens;i++){
    tokendata *toki, *first_token=NULL;
    char c;

    toki = frame->tokens + i;

    c = toki->token[0];
    toki->is_label=0;
    toki->is_string=0;
    toki->is_texturefile=0;
    toki->next=NULL;
    if(first_token==NULL&&c!=':')first_token=toki;
    if(c>='a'&&c<='z'||c>='A'&&c<='Z'){
      int use_displaylist;
      int nargs_actual, noutargs_actual;
      tokendata *this_token, *last_token;
      int error_code;

      toki->type=TOKEN_COMMAND;
      error_code=get_token_id(toki->token, &toki->command, &toki->nvars, &toki->noutvars, &use_displaylist);
	    toki->included_frame=0;
	    toki->included_object=NULL;
      if(error_code==1){
        frame->error=1;
        fprintf(stderr,"*** Error: unable to identify the command, %s, while parsing:\n\n",toki->token);
        fprintf(stderr,"      %s\n\n",object_buffer);
      }
      frame->command_list[ncommands]=toki;
      if(frame->device!=NULL)frame->device->use_displaylist=use_displaylist;
      if(ncommands>0){
        this_token=toki;
        last_token=frame->command_list[ncommands-1];
        last_token->next=this_token;
        this_token->next=NULL;
        nargs_actual = i - last_command_index - 1;
      }
      else{
        nargs_actual = toki-first_token;
        nargs_actual = toki->nvars;
      }
      if(nargs_actual!=toki->nvars){
        frame->error=1;
        if(toki->nvars==1){
          fprintf(stderr,"*** Error: The command %s in device %s has %i arguments, %i was expected\n",
            toki->token,frame->device->label,nargs_actual,toki->nvars);
        }
        else{
          fprintf(stderr,"*** Error: The command %s in device %s has %i arguments, %i were expected\n",
            toki->token,frame->device->label,nargs_actual,toki->nvars);
        }
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
          if(toki->noutvars==1){
            fprintf(stderr,"*** Error: The command %s in device %s has %i output arguments, %i was expected\n",
              toki->token,frame->device->label,noutargs_actual,toki->noutvars);
            }
          else{
            fprintf(stderr,"*** Error: The command %s in device %s has %i output arguments, %i were expected\n",
              toki->token,frame->device->label,noutargs_actual,toki->noutvars);
          }
        }
      }
      ncommands++;
      last_command_index=i;
    }
    else if(c=='$'){
      tokendata *tokdest;

      toki->loc=get_token_loc(toki->token+1,frame);
      if(toki->loc>=0){
        tokdest = frame->tokens+toki->loc;
        toki->varptr=&tokdest->var;
        toki->stringptr=tokdest->string;
        tokdest->reads++;
      }
      else{
        frame->error=1;
        toki->varptr=NULL;
        toki->stringptr=NULL;
        fprintf(stderr,"*** Error: The label %s in device %s is not defined\n",toki->token,frame->device->label);
      }

      toki->type=TOKEN_GETVAL;
    }
    else if(c==':'){
      char *var, *val, *equal;
      char bufcopy[1024];

      strcpy(bufcopy,toki->token);
      var=strtok(bufcopy,"=");
      toki->default_val=0.0;
      strcpy(toki->default_string,"");
      if(var!=NULL){
        val=strtok(NULL,"=");
        if(val!=NULL){
          char *quoted_string;

          quoted_string=strstr(val,"\"");
          if(quoted_string!=NULL){
            int len;

            toki->is_string=1;
            quoted_string++;
            len=strlen(quoted_string);
            if(quoted_string[len-1]=='"')quoted_string[len-1]=' ';
            trim_back(quoted_string);
            quoted_string=trim_front(quoted_string);
            strcpy(toki->default_string,quoted_string);
            quoted_string=strstr(quoted_string,"t%");
            if(quoted_string!=NULL){
              quoted_string+=2;
              quoted_string=trim_front(quoted_string);
              strcpy(toki->default_string,quoted_string);
              toki->is_texturefile=1;
            }
            toki->default_val=0.0;
          }
          else{
            strcpy(toki->tokenlabel,var);
            sscanf(val,"%f",&toki->default_val);
          }
        }
      }
      equal=strchr(toki->token,'=');
      if(equal!=NULL)*equal=0;
      toki->type=TOKEN_FLOAT;
      toki->varptr=&toki->var;
      toki->stringptr=toki->string;
      toki->is_label=1;
      nsymbols++;
    }
    else if(c=='"'){
      char string_copy[256], *sptr;
      int lenstr;
      char *texturefile;

      toki->type=TOKEN_STRING;
      toki->var=0.0;
      toki->varptr=&toki->var;
      toki->stringptr=toki->string;
      sptr=string_copy;
      strcpy(sptr,toki->token);
      sptr++;
      texturefile=strstr(sptr,"t%");
      if(texturefile!=NULL){
        sptr=texturefile+2;
        toki->type=TOKEN_TEXTURE;
        ntextures_local++;
      }
      lenstr=strlen(sptr);
      if(sptr[lenstr-1]=='"')sptr[lenstr-1]=' ';
      trim_back(sptr);
      sptr=trim_front(sptr);
      strcpy(toki->string,sptr);
    }
    else{
      toki->type=TOKEN_FLOAT;
      sscanf(toki->token,"%f",&toki->var);
      toki->varptr=&toki->var;
    }
  }
  frame->ntextures=ntextures_local;
  for(i=0;i<ntokens;i++){
    tokendata *toki;
    char c;

    toki = frame->tokens + i;
    c=toki->token[0];
    if(c!=':')continue;
#ifdef _DEBUG
    if(toki->reads==0){
      fprintf(stderr,"*** Warning: token %s in device %s was not used\n",
        toki->token,frame->device->label);
    }
#endif
  }

  // define data structures for conditional tokens

  for(i=0;i<ncommands;i++){
    tokendata *toki;

    toki = frame->command_list[i];
    switch(toki->command){
      int j,if_level;

      case SV_IF:
        if_level=0;
        for(j=i+1;j<ncommands;j++){
          tokendata *tokj;

          tokj = frame->command_list[j];
          if(tokj->command==SV_IF){
            if_level++;
            continue;
          }
          if(if_level>0&&tokj->command==SV_ENDIF){
            if_level--;
            continue;
          }
          if(if_level==0&&(tokj->command==SV_ELSE||tokj->command==SV_ENDIF)){
            toki->elsenext=frame->command_list[j+1];
            break;
          }
        }
        break;
      case SV_ELSE:
        if_level=0;
        for(j=i+1;j<ncommands;j++){
          tokendata *tokj;

          tokj = frame->command_list[j];
          if(tokj->command==SV_IF){
            if_level++;
            continue;
          }
          if(if_level>0&&tokj->command==SV_ENDIF){
            if_level--;
            continue;
          }
          if(if_level==0&&tokj->command==SV_ENDIF){
            toki->next=frame->command_list[j+1];
            break;
          }
        }
        break;
      default:
        break;
    }
  }
  return buffer_ptr;
}

/* ----------------------- getdevice ----------------------------- */

devicedata *getdevice(char *label,int index){
  int i;

  if(strlen(label)>=4&&strncmp(label,"null",4)==0&&index>=0&&index<ndeviceinfo){
    return deviceinfo + index;
  }
  for(i=0;i<ndeviceinfo;i++){
    devicedata *devicei;

    devicei = deviceinfo + i;
    if(STRCMP(devicei->label,label)==0)return devicei;
  }
  return NULL;
}

/* ----------------------- rewind_device_file ----------------------------- */

void rewind_device_file(FILE *stream){
#define BUFFER_LEN 255
  char buffer[BUFFER_LEN],*comma;
  int found_data=0,buffer_len=BUFFER_LEN;

  fgets(buffer,buffer_len,stream);
  comma=strchr(buffer,',');
  if(comma!=NULL)*comma=0;
  trim_back(buffer);
  if(strcmp(buffer,"//HEADER")!=0){
    rewind(stream);
    return;
  }
  while(!feof(stream)){
    fgets(buffer,buffer_len,stream);
    comma=strchr(buffer,',');
    if(comma!=NULL)*comma=0;
    trim_back(buffer);
    if(strcmp(buffer,"//DATA")==0){
      found_data=1;
      break;
    }
  }
  if(found_data==0){
    fprintf(stderr,"*** Warning //DATA keyword not found in spreadsheet file\n");
  }
  return;
}

/* ----------------------- get_ndevices ----------------------------- */

int get_ndevices(char *file){
  FILE *stream;
  char buffer[BUFFER_LEN],*comma;
  int buffer_len=BUFFER_LEN,nd=0;

  if(file==NULL)return 0;
  stream=fopen(file,"r");
  if(stream==NULL)return 0;
  fgets(buffer,buffer_len,stream);
  comma=strchr(buffer,',');
  if(comma!=NULL)*comma=0;
  trim_back(buffer);
  if(strcmp(buffer,"//HEADER")!=0){
    fclose(stream);
    return 0;
  }

  while(!feof(stream)){
    fgets(buffer,buffer_len,stream);
    comma=strchr(buffer,',');
    if(comma!=NULL)*comma=0;
    trim_back(buffer);
    if(strcmp(buffer,"//DATA")==0){
      break;
    }
    if(strcmp(buffer,"DEVICE")==0){
      nd++;
    }
  }
  fclose(stream);
  return nd;
}

#define EPSDEV 0.01

/* ------------------ comparev2devices ------------------------ */

int comparev2devices(const void *arg1, const void *arg2){
  vdevicesortdata *vdevi, *vdevj;
  float *xyzi, *xyzj;
  int diri, dirj;

  vdevi = (vdevicesortdata *)arg1;
  vdevj = (vdevicesortdata *)arg2;
  diri = vdevi->dir;
  dirj = vdevj->dir;
  xyzi = vdevi->vdeviceinfo->valdev->xyz;
  xyzj = vdevj->vdeviceinfo->valdev->xyz;
  if(diri - dirj < 0)return -1;
  if(diri - dirj > 0)return 1;
  switch(diri){
  case XDIR:
    if(xyzi[1] - xyzj[1]<-EPSDEV)return -1;
    if(xyzi[1] - xyzj[1]>EPSDEV)return 1;
    if(xyzi[2] - xyzj[2]<-EPSDEV)return -1;
    if(xyzi[2] - xyzj[2]>+EPSDEV)return 1;
    break;
  case YDIR:
    if(xyzi[0] - xyzj[0]<-EPSDEV)return -1;
    if(xyzi[0] - xyzj[0]>EPSDEV)return 1;
    if(xyzi[2] - xyzj[2]<-EPSDEV)return -1;
    if(xyzi[2] - xyzj[2]>+EPSDEV)return 1;
    break;
  case ZDIR:
    if(xyzi[0] - xyzj[0]<-EPSDEV)return -1;
    if(xyzi[0] - xyzj[0]>EPSDEV)return 1;
    if(xyzi[1] - xyzj[1]<-EPSDEV)return -1;
    if(xyzi[1] - xyzj[1]>+EPSDEV)return 1;
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  return 0;
}

/* ------------------ comparev3devices ------------------------ */

int comparev3devices( const void *arg1, const void *arg2 ){
  vdevicesortdata *vdevi, *vdevj;
  float *xyzi, *xyzj;
  int diri, dirj;

  vdevi = (vdevicesortdata *)arg1;
  vdevj = (vdevicesortdata *)arg2;
  diri = vdevi->dir;
  dirj = vdevj->dir;
  xyzi = vdevi->vdeviceinfo->valdev->xyz;
  xyzj = vdevj->vdeviceinfo->valdev->xyz;
  if(diri - dirj < 0)return -1;
  if(diri - dirj > 0)return 1;
  switch(diri){
  case XDIR:
    if(xyzi[1]-xyzj[1]<-EPSDEV)return -1;
    if(xyzi[1]-xyzj[1]>+EPSDEV)return 1;
    if(xyzi[2]-xyzj[2]<-EPSDEV)return -1;
    if(xyzi[2]-xyzj[2]>+EPSDEV)return 1;
    if(xyzi[0]-xyzj[0]<-EPSDEV)return -1;
    if(xyzi[0]-xyzj[0]>EPSDEV)return 1;
    break;
  case YDIR:
    if(xyzi[0]-xyzj[0]<-EPSDEV)return -1;
    if(xyzi[0]-xyzj[0]>EPSDEV)return 1;
    if(xyzi[2]-xyzj[2]<-EPSDEV)return -1;
    if(xyzi[2]-xyzj[2]>+EPSDEV)return 1;
    if(xyzi[1]-xyzj[1]<-EPSDEV)return -1;
    if(xyzi[1]-xyzj[1]>+EPSDEV)return 1;
    break;
  case ZDIR:
    if(xyzi[0]-xyzj[0]<-EPSDEV)return -1;
    if(xyzi[0]-xyzj[0]>EPSDEV)return 1;
    if(xyzi[1]-xyzj[1]<-EPSDEV)return -1;
    if(xyzi[1]-xyzj[1]>+EPSDEV)return 1;
    if(xyzi[2]-xyzj[2]<-EPSDEV)return -1;
    if(xyzi[2]-xyzj[2]>+EPSDEV)return 1;
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  return 0;
}

/* ----------------------- setup_tree_devices ----------------------------- */

void setup_tree_devices(void){
  int i;
  treedevicedata *treei;

  if(nvdeviceinfo==0)return;
  if(ntreedeviceinfo>0){
    FREEMEMORY(treedeviceinfo);
    ntreedeviceinfo=0;
  }

  qsort((vdevicedata **)vdevices_sorted,3*(size_t)nvdeviceinfo,sizeof(vdevicesortdata),comparev3devices);

  ntreedeviceinfo = 1;
  for(i = 1; i < 3*nvdeviceinfo; i++){
    if(comparev2devices(vdevices_sorted+i, vdevices_sorted+i-1) != 0)ntreedeviceinfo++;
  }

  NewMemory((void **)&treedeviceinfo,ntreedeviceinfo*sizeof(treedevicedata));

  ntreedeviceinfo = 1;
  treei = treedeviceinfo;
  treei->first = 0;
  for(i = 1; i < 3*nvdeviceinfo; i++){
    if(comparev2devices(vdevices_sorted + i, vdevices_sorted + i - 1) != 0){
      treei->last = i-1;
      treei = treedeviceinfo + ntreedeviceinfo;
      treei->first = i;
      ntreedeviceinfo++;
    }
  }
  treei->last = 3*nvdeviceinfo - 1;

  max_device_tree=0;
  for(i = 0; i < ntreedeviceinfo; i++){
    int j, n;

    treei = treedeviceinfo + i;
    n = 0;
    for(j = treei->first; j <= treei->last; j++){
      vdevicedata *vdevi;
      vdevicesortdata *vdevsorti;

      vdevsorti = vdevices_sorted + j;
      vdevi = vdevsorti->vdeviceinfo;
      if(vdevi->unique != 0)n++;
    }
    treei->n = n;
    max_device_tree=MAX(max_device_tree,n);
  }
}

/* ----------------------- setup_zone_devs ----------------------------- */

void setup_zone_devs(void){
  int i;

  show_missing_objects = 0;
  for(i=0;i<nzoneinfo;i++){
    FILE *stream;
    char *file;
    int nrows, ncols, buffer_len,ntokens;
    char *buffer=NULL,**devclabels=NULL;
    zonedata *zonei;
    int j;

    zonei = zoneinfo + i;
    if(zonei->csv!=1)continue;
    file = zonei->file;

    stream=fopen(file,"r");
    if(stream==NULL)continue;
    buffer_len=getrowcols(stream,&nrows,&ncols);
    if(nrows<=0||ncols<=0||buffer_len<=0){
      fclose(stream);
      continue;
    }
    buffer_len+=10;
    rewind(stream);

    NewMemory((void **)&buffer,buffer_len);
    NewMemory((void **)&devclabels,ncols*sizeof(char *));
    fgets(buffer,buffer_len,stream);
    fgets(buffer,buffer_len,stream);
    parsecsv(buffer,devclabels,ncols,&ntokens);
    for(j=0;j<ntokens;j++){
      devicedata *devi;

      trim_back(devclabels[j]);
      devclabels[j]=trim_front(devclabels[j]);
      devi = getdevice(devclabels[j],-1);
      if(devi!=NULL)devi->in_zone_csv=1;
    }
    FREEMEMORY(devclabels);
    FREEMEMORY(buffer);
    fclose(stream);
  }
}

/* ----------------------- read_device_data ----------------------------- */

void read_device_data(char *file, int filetype, int loadstatus){
  FILE *stream;
  int nrows, ncols;
  int irow;
  float *vals=NULL;
  int *valids=NULL;
  int i;
  char *buffer, *buffer2;
  char **devcunits=NULL, **devclabels=NULL;
  devicedata **devices=NULL;
  int ntokens;
  int buffer_len;
  float *times_local;

// unload data

  if(loadstatus==UNLOAD){
    for(i=0;i<ndeviceinfo;i++){
      devicedata *devicei;

      devicei = deviceinfo + i;
      if(devicei->filetype!=filetype)continue;
      FREEMEMORY(devicei->vals);
      FREEMEMORY(devicei->valids);
    }
    for(i=0;i<ndeviceinfo;i++){
      devicedata *devicei;
      int j;

      devicei = deviceinfo + i;
      if(devicei->filetype!=filetype||devicei->times==NULL)continue;
      times_local = devicei->times;
      FREEMEMORY(devicei->times);
      for(j=i+1;j<ndeviceinfo;j++){
        devicedata *devicej;

        devicej = deviceinfo + j;
        if(devicej->filetype!=filetype)continue;
        if(times_local==devicej->times)devicej->times=NULL;
      }
    }
    return;
  }

  // find number of rows and columns

  stream=fopen(file,"r");
  if(stream==NULL)return;
  rewind_device_file(stream);
  buffer_len=getrowcols(stream,&nrows,&ncols);
  if(nrows<=0||ncols<=0||buffer_len<=0){
    fclose(stream);
    return;
  }
  buffer_len+=10;
  rewind_device_file(stream);

  NewMemory((void **)&buffer,buffer_len);
  NewMemory((void **)&buffer2,buffer_len);
  NewMemory((void **)&vals,ncols*sizeof(float));
  NewMemory((void **)&valids,ncols*sizeof(int));
  NewMemory((void **)&devcunits,ncols*sizeof(char *));
  NewMemory((void **)&devclabels,ncols*sizeof(char *));
  NewMemory((void **)&devices,ncols*sizeof(devicedata *));

  for(i = 0; i<ncols; i++){
    devices[i] = NULL;
  }

  fgets(buffer,buffer_len,stream);
  parsecsv(buffer,devcunits,ncols,&ntokens);
  for(i=0;i<ntokens;i++){
    trim_back(devcunits[i]);
    devcunits[i]=trim_front(devcunits[i]);
  }

  fgets(buffer2,buffer_len,stream);
  parsecsv(buffer2,devclabels,ncols,&ntokens);
  for(i=0;i<ntokens;i++){
    trim_back(devclabels[i]);
    devclabels[i]=trim_front(devclabels[i]);
  }

  NewMemory((void **)&times_local,nrows*sizeof(float));
  for(i=1;i<ntokens;i++){
    devicedata *devicei;

    devicei = getdevice(devclabels[i],i-1);
    devices[i]=devicei;
#ifdef _DEBUG
    if(devicei==NULL){
      fprintf(stderr,"*** Error: spreadsheet entry: %s is not present in %s\n",devclabels[i],smv_filename);
      continue;
    }
#endif
    if(devicei==NULL)continue;
    devicei->filetype=filetype;
    if(filetype==CSV_FDS)devicei->in_devc_csv=1;
    NewMemory((void **)&devicei->vals,nrows*sizeof(float));
    NewMemory((void **)&devicei->valids,nrows*sizeof(int));
    devicei->times=times_local;
#ifdef pp_DEG
    if(strcmp(devcunits[i],"C")==0){
      strcpy(devicei->unit,degC);
    }
    else{
      strcpy(devicei->unit,devcunits[i]);
    }
#else
    strcpy(devicei->unit,devcunits[i]);
#endif
    devicei->nvals=nrows-2;
  }

  for(irow=2;irow<nrows;irow++){
    int icol=0;

    fgets(buffer,buffer_len,stream);
    fparsecsv(buffer,vals,valids,ncols,&ntokens);
    times_local[irow-2]=vals[icol];
    for(icol=1;icol<ncols;icol++){
      devicedata *devicei;

      devicei = devices[icol];
      if(devicei==NULL)continue;
      devicei->vals[irow-2]=vals[icol];
      devicei->valids[irow-2]=valids[icol];
    }
  }
  FREEMEMORY(buffer);
  FREEMEMORY(buffer2);
  fclose(stream);

  FREEMEMORY(vals);
  FREEMEMORY(valids);
  FREEMEMORY(devcunits);
  FREEMEMORY(devclabels)
  FREEMEMORY(devices);
}

/* ----------------------- get_device ----------------------------- */

devicedata *get_device(float *xyzval, char *device_label, int device_type){
  int j;

  for(j=0;j<ndeviceinfo;j++){
    devicedata *devj;
    float *xyz;

    devj = deviceinfo + j;
    if(devj->filetype!=device_type)continue;
    xyz = devj->xyz;
    if(strcmp(devj->quantity,device_label)!=0)continue;
    if(ABS(xyz[0]-xyzval[0])>EPSDEV)continue;
    if(ABS(xyz[1]-xyzval[1])>EPSDEV)continue;
    if(ABS(xyz[2]-xyzval[2])>EPSDEV)continue;
    return devj;
  }
  return NULL;
}

/* ----------------------- get_vdevice ----------------------------- */

vdevicedata *get_vdevice(float *xyzval){
  int j;

  for(j=0;j<nvdeviceinfo;j++){
    vdevicedata *vdevj;
    float *xyzj;

    vdevj = vdeviceinfo + j;

    xyzj = vdevj->valdev->xyz;
    if(ABS(xyzval[0]-xyzj[0])>EPSDEV)continue;
    if(ABS(xyzval[1]-xyzj[1])>EPSDEV)continue;
    if(ABS(xyzval[2]-xyzj[2])>EPSDEV)continue;
    return vdevj;
  }
  return NULL;
}

/* ----------------------- update_colordevs ----------------------------- */

void update_colordevs(void){
  int i;
  devicedata *colordev;

  colordev = devicetypes[devicetypes_index];

  for(i=0;i<nvdeviceinfo;i++){
    vdevicedata *vdevi;

    vdevi = vdeviceinfo + i;
    vdevi->colordev=NULL;
  }
  for(i=0;i<ndeviceinfo;i++){
    devicedata *devi;
    vdevicedata *vdevi;
    devi = deviceinfo + i;
    vdevi = devi->vdevice;
    if(vdevi==NULL)continue;
    if(strcmp(colordev->quantity,devi->quantity)==0){
      vdevi->colordev=devi;
    }
  }
}

/* ----------------------- is_dup_device_label ----------------------------- */
#define BEFORE 0
#define AFTER 1

int is_dup_device_label(int index, int direction){
  int i,i1,i2;
  devicedata *dev_index;

  if(direction==BEFORE){
    i1=0;
    i2=index;
  }
  else{
    i1=index+1;
    i2=ndeviceinfo;
  }
  dev_index = deviceinfo + index;
  if(index<0||index>=ndeviceinfo||dev_index->label==NULL||STRCMP(dev_index->label,"null")==0||dev_index->in_devc_csv==0)return 0;

  for(i=i1;i<i2;i++){
    devicedata *devi;

    devi = deviceinfo + i;
    if(devi->label==NULL||STRCMP(devi->label,"null")==0)continue;
    if(STRCMP(dev_index->label,devi->label)==0)return 1;
  }
  return 0;
}

/* ----------------------- setup_pilot_data ----------------------------- */

#ifdef pp_PILOT
#ifdef pp_WINDROSE
void setup_pilot_data(int nbuckets, int nr, int ntheta, int flag){
#else
void setup_pilot_data(int nbuckets){
#endif
  int i;
  float dangle;

  dangle = 360.0 / (float)nbuckets;
  for(i = 0; i < nvdeviceinfo; i++){
    vdevicedata *vdevicei;
    devicedata *udev, *vdev, *wdev;
    devicedata *angledev, *veldev;
    int j, ibucket;
	pilotdata *piloti;

    vdevicei = vdeviceinfo + i;
    udev = vdevicei->udev;
    vdev = vdevicei->vdev;
    wdev = vdevicei->wdev;
    angledev = vdevicei->angledev;
    veldev = vdevicei->veldev;

    piloti = &(vdevicei->pilotinfo);
    {
       float *vel, *fraction;

       vel = piloti->vel;
       fraction = piloti->fraction;
       vel = piloti->vel;
       FREEMEMORY(fraction);
       FREEMEMORY(vel);
       NewMemory((void **)&fraction, nbuckets*sizeof(float));
       NewMemory((void **)&vel, nbuckets*sizeof(float));
       piloti->vel = vel;
       piloti->fraction = fraction;
       piloti->nbuckets = nbuckets;
    }

    for(j = 0; j < nbuckets; j++){
      piloti->fraction[j] = 0.0;
      piloti->vel[j] = 0.0;
    }
    piloti->total = 0;
    if(udev != NULL&&vdev != NULL){
      int nvals;

      nvals = MIN(udev->nvals, vdev->nvals);
      if(wdev!=NULL)nvals = MIN(nvals, wdev->nvals);
      for(j = 0; j<nvals; j++){
        float uval, vval, wval = 0.0, vel, veluv, angle;

        uval = udev->vals[j];
        vval = vdev->vals[j];
        if(wdev != NULL)wval = wdev->vals[j];
        vel = sqrt(uval*uval + vval*vval + wval*wval);
        veluv = sqrt(uval*uval + vval*vval);
        if(veluv>0.0){
          angle = fmod(180.0 + atan2(vval, uval)*RAD2DEG + dangle/2.0, 360.0);
          ibucket = CLAMP(angle / dangle, 0, nbuckets-1);
          piloti->fraction[ibucket]++;
          piloti->vel[ibucket] += vel;
        }
      }
#ifdef pp_WINDROSE
      {
        float rmin, rmax;
        histogramdata *histogram;

        histogram = &(piloti->histogram);
        if(flag != FIRST_TIME){
          free_histogram2d(histogram);
        }
        init_histogram2d(histogram, nr, ntheta);
        get_2dminmax(udev->vals, vdev->vals, nvals, &rmin, &rmax, HIST_COMPUTE_BOUNDS);
        copy_uvdata2histogram(udev->vals,vdev->vals,nvals,rmin,rmax,histogram);
      }
#endif
    }
    else if(angledev != NULL&&veldev != NULL){
      int nvals;

      nvals = MIN(angledev->nvals, veldev->nvals);
      for(j = 0; j < nvals; j++){
        float vel, angle;

        angle = angledev->vals[j];
        vel = veldev->vals[j];
        angle = fmod(angle + dangle/2.0, 360.0);
        ibucket = CLAMP(angle / dangle, 0, nbuckets-1);
        piloti->fraction[ibucket]++;
        piloti->vel[ibucket] += vel;
      }
#ifdef pp_WINDROSE
      {
        float rmin, rmax;
        histogramdata *histogram;

        histogram = &(piloti->histogram);
        if(flag != FIRST_TIME){
          free_histogram2d(histogram);
        }
        init_histogram2d(histogram, nr, ntheta);
        get_polarminmax(veldev->vals, nvals, &rmin, &rmax, HIST_COMPUTE_BOUNDS);
        copy_polardata2histogram(veldev->vals,angledev->vals,nvals,rmin,rmax,histogram);
      }
#endif
    }
    else{
      continue;
    }
    for(j = 0; j<nbuckets; j++){
      piloti->total += piloti->fraction[j];
      if(piloti->fraction[j]>0.0){
        piloti->vel[j] /= piloti->fraction[j];
      }
    }
    if(piloti->total > 0){
      for(j = 0; j < nbuckets; j++){
        piloti->fraction[j] /= piloti->total;
      }
    }
  }
}
#endif

/* ----------------------- setup_device_data ----------------------------- */

void setup_device_data(void){
  float *vals=NULL;
  int *valids=NULL;
  int i;
  char **devcunits=NULL, **devclabels=NULL;
  int is_dup;

  if(ndeviceinfo==0)return;
  FREEMEMORY(vdeviceinfo);
  NewMemory((void **)&vdeviceinfo,ndeviceinfo*sizeof(vdevicedata));
  FREEMEMORY(vdevices_sorted);
  NewMemory((void **)&vdevices_sorted,3*ndeviceinfo*sizeof(vdevicesortdata));
  nvdeviceinfo=0;
  for(i=0;i<ndeviceinfo;i++){
    vdevicedata *vdevi;
    devicedata *devi,*devj;
    float *xyzval;

    devi = deviceinfo + i;
    xyzval=devi->xyz;
    devi->vdevice=NULL;

    vdevi = vdeviceinfo + nvdeviceinfo;
    vdevi->valdev = devi;
    vdevi->udev=NULL;
    vdevi->vdev=NULL;
    vdevi->wdev=NULL;
    vdevi->angledev=NULL;
    vdevi->veldev=NULL;
    vdevi->sd_angledev=NULL;
    vdevi->sd_veldev=NULL;
    vdevi->colordev=NULL;
#ifdef pp_PILOT
    vdevi->pilotinfo.vel=NULL;
    vdevi->pilotinfo.fraction=NULL;
    vdevi->pilotinfo.nbuckets=0;
#endif

    devj = get_device(xyzval,"VELOCITY",CSV_EXP);
    if(devj!=NULL){
      vdevi->veldev=devj;
      vdevi->filetype=CSV_EXP;
    }

    devj = get_device(xyzval,"SD_VELOCITY",CSV_EXP);
    if(devj!=NULL){
      vdevi->sd_veldev=devj;
      vdevi->filetype=CSV_EXP;
    }

    devj = get_device(xyzval,"ANGLE",CSV_EXP);
    if(devj!=NULL){
      vdevi->angledev=devj;
      vdevi->filetype=CSV_EXP;
    }

    devj = get_device(xyzval,"SD_ANGLE",CSV_EXP);
    if(devj!=NULL){
      vdevi->sd_angledev=devj;
      vdevi->filetype=CSV_EXP;
    }

    devj = get_device(xyzval,"U-VELOCITY",CSV_FDS);
    if(devj!=NULL){
      vdevi->udev=devj;
      vdevi->filetype=CSV_FDS;
    }

    devj = get_device(xyzval,"V-VELOCITY",CSV_FDS);
    if(devj!=NULL){
      vdevi->vdev=devj;
      vdevi->filetype=CSV_FDS;
    }

    devj = get_device(xyzval,"W-VELOCITY",CSV_FDS);
    if(devj!=NULL){
      vdevi->wdev=devj;
      vdevi->filetype=CSV_FDS;
    }

    if(vdevi->udev!=NULL||vdevi->vdev!=NULL||vdevi->wdev!=NULL||
      vdevi->angledev!=NULL||vdevi->veldev!=NULL){
      vdevi->unique=1;
      nvdeviceinfo++;
    }
  }

  // look for duplicate device labels

  is_dup=0;
  for(i=0;i<ndeviceinfo;i++){
    devicedata *devi;

    devi = deviceinfo + i;
    if(devi->label==NULL||STRCMP(devi->label,"null")==0)continue;
    if(is_dup_device_label(i,AFTER)==1){
      is_dup=1;
      break;
    }
  }
  if(is_dup==1){
    int ii;

    fprintf(stderr,"*** Warning: Duplicate device labels: ");
    for(ii=0;ii<ndeviceinfo;ii++){
      devicedata *devi;

      devi = deviceinfo + ii;
      if(devi->label==NULL||STRCMP(devi->label,"null")==0)continue;
      if(is_dup_device_label(ii,BEFORE)==0&&is_dup_device_label(ii,AFTER)==1){
        fprintf(stderr," %s,",devi->label);
      }
    }
    fprintf(stderr," found in %s\n",fds_filein);
  }
  for(i=0;i<nvdeviceinfo;i++){
    vdevicedata *vdevi;
    int j;
    float *xyzi;

    vdevi = vdeviceinfo + i;
    xyzi = vdevi->valdev->xyz;
    for(j=i+1;j<nvdeviceinfo;j++){
      vdevicedata *vdevj;
      float *xyzj;

      vdevj = vdeviceinfo + j;
      if(vdevj->unique==0)continue;
      xyzj = vdevj->valdev->xyz;
      if(ABS(xyzi[0]-xyzj[0])>EPSDEV)continue;
      if(ABS(xyzi[1]-xyzj[1])>EPSDEV)continue;
      if(ABS(xyzi[2]-xyzj[2])>EPSDEV)continue;
      vdevj->unique=0;
    }
  }
  max_dev_vel=-1.0;
  for(i=0;i<nvdeviceinfo;i++){
    vdevicedata *vdevi;
    devicedata *devval;
    int j;

    vdevi = vdeviceinfo + i;
    if(vdevi->unique==0)continue;
    devval = vdevi->valdev;
    if(vdevi->udev!=NULL)vdevi->udev->vdevice=vdevi;
    if(vdevi->vdev!=NULL)vdevi->vdev->vdevice=vdevi;
    if(vdevi->wdev!=NULL)vdevi->wdev->vdevice=vdevi;
    if(vdevi->filetype==CSV_FDS){
      devicedata *udev,*vdev,*wdev;

      udev=vdevi->udev;
      vdev=vdevi->vdev;
      wdev=vdevi->wdev;
      for(j=0;j<devval->nvals;j++){
        float uvel=0.0, vvel=0.0, wvel=0.0;
        float speed;

        if(udev!=NULL)uvel=udev->vals[j];
        if(vdev!=NULL)vvel=vdev->vals[j];
        if(wdev!=NULL)wvel=wdev->vals[j];
        speed = sqrt(uvel*uvel+vvel*vvel+wvel*wvel);
        if(speed>max_dev_vel)max_dev_vel=speed;
      }
    }
    if(vdevi->filetype==CSV_EXP){
      devicedata *veldev;

      veldev=vdevi->veldev;
      if(veldev!=NULL){
        for(j=0;j<devval->nvals;j++){
          if(veldev->valids[j]==1){
            float speed;

            speed=veldev->vals[j];
            if(speed>max_dev_vel)max_dev_vel=speed;
          }
        }
      }
    }
  }

  // find devices linked with each vdevice

  if(ndeviceinfo>0){
    for(i=0;i<ndeviceinfo;i++){
      devicedata *devi;
      float *xyzi;
      vdevicedata *vdevj;

      devi = deviceinfo + i;
      if(devi->vdevice!=NULL)continue;
      xyzi = devi->xyz;
      vdevj = get_vdevice(xyzi);
      if(vdevj!=NULL)devi->vdevice=vdevj;
    }
  }

  if(ndeviceinfo>0){
    ndevicetypes=0;
    FREEMEMORY(devicetypes);
    NewMemory((void **)&devicetypes,ndeviceinfo*sizeof(devicedata *));
    for(i=0;i<ndeviceinfo;i++){
      devicedata *devi;

      devi = deviceinfo + i;
      devi->type2=-1;
    }
    for(i=0;i<ndeviceinfo;i++){
      int j;
      devicedata *devi;

      devi = deviceinfo + i;
      if(devi->type2>=0||strlen(devi->quantity)==0)continue;
      devi->type2=ndevicetypes;
      devi->type2vis=0;
      devicetypes[ndevicetypes++]=devi;
      for(j=i+1;j<ndeviceinfo;j++){
        devicedata *devj;

        devj = deviceinfo + j;
        if(devj->type2<0&&strcmp(devi->quantity,devj->quantity)==0){
          devj->type2=devi->type2;
        }
      }
    }
    if(ndevicetypes>0)devicetypes[0]->type2vis=1;
  }
  for(i=0;i<nvdeviceinfo;i++){
    vdevicesortdata *vdevsorti;

    vdevsorti = vdevices_sorted + i;
    vdevsorti->vdeviceinfo = vdeviceinfo + i;
    vdevsorti->dir = XDIR;

    vdevsorti = vdevices_sorted + nvdeviceinfo + i;
    vdevsorti->vdeviceinfo = vdeviceinfo + i;
    vdevsorti->dir = YDIR;

    vdevsorti = vdevices_sorted + 2*nvdeviceinfo + i;
    vdevsorti->vdeviceinfo = vdeviceinfo + i;
    vdevsorti->dir = ZDIR;
  }

  setup_tree_devices();
  update_colordevs();

  // convert velocities to pilot chart format
#ifdef pp_PILOT
#ifdef pp_WINDROSE
  setup_pilot_data(npilot_buckets,npilot_nr,npilot_ntheta,FIRST_TIME);
#else
  setup_pilot_data(npilot_buckets);
#endif
#endif

  FREEMEMORY(vals);
  FREEMEMORY(valids);
  FREEMEMORY(devcunits);
  FREEMEMORY(devclabels)
}

/* ----------------------- read_object_defs ----------------------------- */

int read_object_defs(char *file){
  FILE *stream;
  char buffer[256], *trim_buffer;
  char *buffer_ptr;
  sv_object *temp_object, *prev_object, *next_object, *current_object;
  sv_object_frame *current_frame;
  int firstdef;
  sv_object *object_start, *objecti;
  size_t lenbuffer;
  int ndevices=0;
  int eof=0;

 // freeall_objects();

  stream=fopen(file,"r");
  if(stream==NULL)return 0;
  PRINTF("Processing object file:  %s\n",file);

  firstdef=-1;
  buffer_ptr=NULL;
  while(!feof(stream)){
    CheckMemory;
    if(buffer_ptr==NULL){
      if(eof==1||fgets(buffer,255,stream)==NULL)break;
      buffer_ptr=buffer;
    }
    trim_buffer=remove_comment(buffer_ptr);
    lenbuffer=strlen(buffer_ptr);
    if(lenbuffer<1){
      buffer_ptr=NULL;
      continue;
    }


    if(match(buffer_ptr,"OBJECTDEF") == 1||
       match(buffer_ptr,"AVATARDEF") == 1
      ){
        int object_type=IS_NOT_AVATAR;
      char *label;

      sv_object_frame *first_frame, *last_frame;

      if(match(buffer_ptr,"AVATARDEF") == 1){
        object_type=IS_AVATAR;
      }
      ndevices++;
      if(fgets(buffer,255,stream)==NULL)break;
      label=remove_comment(buffer);
      temp_object=get_object(label);
      if(temp_object!=NULL){
        free_object(temp_object);
      }

      NewMemory((void **)&current_object,sizeof(sv_object));
      current_object->used=0;
      current_object->use_displaylist=1;
      current_object->select_mode=0;
      strcpy(current_object->label,label);
      prev_object = object_def_last.prev;
      next_object = &object_def_last;

      prev_object->next=current_object;
      next_object->prev=current_object;

      current_object->next=next_object;
      current_object->prev=prev_object;
      current_object->visible=1;

      current_object->nframes=0;

      first_frame = &current_object->first_frame;
      last_frame = &current_object->last_frame;
      current_object->type=object_type;

      first_frame->next=last_frame;
      first_frame->prev=NULL;
      last_frame->prev=first_frame;
      last_frame->next=NULL;

      firstdef=1;
      buffer_ptr=NULL;
      continue;
    }
    if(match(trim_buffer,"NEWFRAME") == 1||firstdef==1){
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
      current_frame->nevac_tokens=0;

      current_object->nframes++;

      firstdef=0;
      if(match(trim_buffer,"NEWFRAME")==1){
        buffer_ptr=NULL;
        continue;
      }
    }
    buffer_ptr=parse_device_frame(buffer,stream,&eof,current_frame);
  }
  fclose(stream);

  object_start = object_def_first.next;
  objecti = object_start;
  nobject_defs=0;
  for(;objecti->next!=NULL;){
    CheckMemory;
    nobject_defs++;
    objecti->obj_frames=NULL;
    if(objecti->nframes>0){
      NewMemory((void **)&objecti->obj_frames,objecti->nframes*sizeof(sv_object_frame *));
    }
    objecti=objecti->next;
  }
  FREEMEMORY(object_defs);
  if(nobject_defs>0){
    int i,j;

    NewMemory((void **)&object_defs,nobject_defs*sizeof(sv_object *));

    object_start = object_def_first.next;
    objecti = object_start;
    i=0;
    for(;objecti->next!=NULL;){
      sv_object_frame *frame_start, *framei;

      CheckMemory;
      object_defs[i]=objecti;
      i++;
      frame_start = objecti->first_frame.next;
      framei = frame_start;
      j=0;
      for(;framei->next!=NULL;){
        int npushpop=0, ii;

        CheckMemory;
        objecti->obj_frames[j]=framei;
        for(ii=0;ii<framei->ncommands;ii++){
          tokendata *command;
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
          fprintf(stderr,"*** Error: The number of push and pop commands are not equal\n");
          framei->error=1;
        }
        framei=framei->next;
        j++;
      }
      objecti=objecti->next;
    }
  }
  PRINTF("Object file processing complete\n\n");
  return ndevices;
}

/* ----------------------- reporterror ----------------------------- */

void reporterror(char *buffer, char *token, int numargs_found, int numargs_expected){
  if(numargs_found==numargs_expected)return;
  fprintf(stderr,"*** Error: %i arguments were found (%i expected) for the token, %s, while parsing: %s\n",numargs_found,numargs_expected,token,buffer);
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
  trim_back(label_present);
  if(strlen(label_present)==0)return NULL;
  return label_present;
}

/* ----------------------- get_object ----------------------------- */

sv_object *get_object(char *label){
  sv_object *objecti,*object_start;

  object_start = object_def_first.next;
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
    object = object_def_last.prev;
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
    if(framei->nsymbols>0){
      FREEMEMORY(framei->symbols);
    }
    if(framei->ntokens>0){
      FREEMEMORY(framei->tokens);
    }
    FREEMEMORY(framei);
    framei=next_frame;
  }
  FREEMEMORY(object);
}

/* ----------------------- update_device_textures ----------------------------- */

void update_device_textures(void){

  // create a list of device textures

  int i;

  for(i=0;i<ndeviceinfo;i++){
    devicedata *devicei;

    devicei = deviceinfo + i;

    if(devicei->object==NULL){
      devicei->object = get_SVOBJECT_type(devicei->label,missing_device);
    }
  }

  device_texture_list=NULL;
  ndevice_texture_list=0;

  // count device textures

  for(i=0;i<ndeviceinfo;i++){
    devicedata *devicei;
    sv_object *object;
    int j;

    devicei = deviceinfo + i;
    object = devicei->object;
    for(j=0;j<object->nframes;j++){
      sv_object_frame *frame;

      frame = object->obj_frames[j];
      ndevice_texture_list+=frame->ntextures;
    }
  }
  for(i=0;i<npropinfo;i++){
    propdata *propi;

    propi = propinfo + i;

    ndevice_texture_list += propi->ntextures;
  }

  // allocate data structures and fill in list

  if(ndevice_texture_list>0){
    NewMemory((void **)&device_texture_list,ndevice_texture_list*sizeof(char *));
    NewMemory((void **)&device_texture_list_index,ndevice_texture_list*sizeof(int));
    ndevice_texture_list=0;
    for(i=0;i<ndeviceinfo;i++){
      devicedata *devicei;
      sv_object *object;
      int j;

      devicei = deviceinfo + i;
      object = devicei->object;
      for(j=0;j<object->nframes;j++){
        sv_object_frame *frame;
        int k;

        frame = object->obj_frames[j];
        if(frame->ntextures==0)continue;
        for(k=0;k<frame->ntokens;k++){
          tokendata *toki;
          int kk;
          int dup;

          toki = frame->tokens + k;
          if(toki->type!=TOKEN_TEXTURE)continue;
          dup=0;
          for(kk=0;kk<ndevice_texture_list;kk++){
            if(strcmp(device_texture_list[kk],toki->string)==0){
              dup=1;
              break;
            }
          }
          if(dup==0)device_texture_list[ndevice_texture_list++]=toki->string;
        }
      }
    }
    for(i=0;i<npropinfo;i++){
      propdata *propi;
      int j;

      propi = propinfo + i;
      if(propi->ntextures==0)continue;
      for(j=0;j<propi->ntextures;j++){
        int dup;
        char *texturefile;
        int kk;

        texturefile=propi->texturefiles[j];
        dup=0;
        for(kk=0;kk<ndevice_texture_list;kk++){
          if(strcmp(device_texture_list[kk],texturefile)==0){
            dup=1;
            break;
          }
        }
        if(dup==0)device_texture_list[ndevice_texture_list++]=texturefile;
      }
    }
  }
}

/* ----------------------- init_object_defs ----------------------------- */

void init_object_defs(void){
  char com_buffer[1024];
  char com_buffer2[1024];
  char objectfile[1024];
  int i;

  svofile_exists = 0;

  if(smokeview_bindir!=NULL){
    strcpy(objectfile,smokeview_bindir);
    strcat(objectfile,"objects.svo");
    read_object_defs(objectfile);
  }

  strcpy(objectfile,"objects.svo");
  read_object_defs(objectfile);

  strcpy(objectfile,fdsprefix);
  strcat(objectfile,".svo");
  read_object_defs(objectfile);

  init_avatar();

  if(isZoneFireModel==1){
    strcpy(com_buffer,"255 255 0 setrgb 0.02 0.05 drawdisk");
    target_object_backup = init_SVOBJECT1("target", com_buffer,1);
  }
  else{
    strcpy(com_buffer,"255 255 0 setrgb 0.038 drawcube");
    target_object_backup = init_SVOBJECT1("sensor", com_buffer,1);
  }

  strcpy(com_buffer,"255 255 0 setrgb 0.038 drawcube");
  thcp_object_backup = init_SVOBJECT1("thcp", com_buffer,1);

  strcpy(com_buffer, "0 255 0 setrgb 0.038 drawcube");
  strcpy(com_buffer2,"255 0 0 setrgb 0.038 drawcube");
  heat_detector_object_backup = init_SVOBJECT2("heat_detector", com_buffer, com_buffer2,1);

  strcpy(com_buffer, "0 255 0 setrgb 0.038 drawcube");
  strcpy(com_buffer2,"255 0 0 setrgb 0.038 drawcube");
  sprinkler_upright_object_backup = init_SVOBJECT2("sprinkler_upright", com_buffer, com_buffer2,1);

  strcpy(com_buffer, "127 127 127 setrgb 0.2 0.05 drawdisk");
  strcpy(com_buffer2,"255 0 0 setrgb 0.2 0.05 drawdisk");
  smoke_detector_object_backup = init_SVOBJECT2("smoke_detector", com_buffer, com_buffer2,1);

  strcpy(com_buffer, "255 0 0 setrgb push 45.0 rotatey -0.1 offsetz 0.05 0.2 drawdisk pop push -45.0 rotatey -0.1 offsetz 0.05 0.2 drawdisk pop");
  error_device = init_SVOBJECT1("error_device", com_buffer,1);

  if(missing_device==NULL){
    strcpy(com_buffer, "0 0 255 setrgb push 45.0 rotatey -0.1 offsetz 0.05 0.2 drawdisk pop push -45.0 rotatey -0.1 offsetz 0.05 0.2 drawdisk pop");
    missing_device = init_SVOBJECT1("missing_device", com_buffer,1);
  }

  if(nobject_defs==0){
    nobject_defs=4;
    FREEMEMORY(object_defs);
    NewMemory((void **)&object_defs,4*sizeof(sv_object *));
    object_defs[0] = target_object_backup;
    object_defs[1] = heat_detector_object_backup;
    object_defs[2] = sprinkler_upright_object_backup;
    object_defs[3] = smoke_detector_object_backup;
  }

  for(i=0;i<navatar_types;i++){
    sv_object_frame *obj_frame;
    int n;
    tokendata **evac_tokens,*evac_token;

    CheckMemory;

    obj_frame=avatar_types[i]->obj_frames[0];
    evac_tokens = obj_frame->evac_tokens;

    n=0;

    evac_token=get_token_ptr("W",obj_frame);
    evac_tokens[n++]=evac_token;

    evac_token=get_token_ptr("D",obj_frame);
    evac_tokens[n++]=evac_token;

    evac_token=get_token_ptr("H1",obj_frame);
    evac_tokens[n++]=evac_token;

    evac_token=get_token_ptr("SX",obj_frame);
    evac_tokens[n++]=evac_token;

    evac_token=get_token_ptr("SY",obj_frame);
    evac_tokens[n++]=evac_token;

    evac_token=get_token_ptr("SZ",obj_frame);
    evac_tokens[n++]=evac_token;

    evac_token=get_token_ptr("R",obj_frame);
    evac_tokens[n++]=evac_token;

    evac_token=get_token_ptr("G",obj_frame);
    evac_tokens[n++]=evac_token;

    evac_token=get_token_ptr("B",obj_frame);
    evac_tokens[n++]=evac_token;

    evac_token=get_token_ptr("HX",obj_frame);
    evac_tokens[n++]=evac_token;

    evac_token=get_token_ptr("HY",obj_frame);
    evac_tokens[n++]=evac_token;

    evac_token=get_token_ptr("HZ",obj_frame);
    evac_tokens[n++]=evac_token;
  }
}

/* ----------------------- update_object_used ----------------------------- */

void update_object_used(void){
  int i;

  for(i = 0; i<nobject_defs; i++){
    sv_object *obj_typei;

    obj_typei = object_defs[i];
    obj_typei->used_by_device = 0;
  }
  for(i = 0; i<ndeviceinfo; i++){
    devicedata *devicei;
    propdata *propi;
    int jj;

    devicei = deviceinfo+i;
    propi = devicei->prop;
    if(propi==NULL)continue;
    for(jj = 0; jj<propi->nsmokeview_ids; jj++){
      sv_object *objectj;

      objectj = propi->smv_objects[jj];
      objectj->used_by_device = 1;
    }
  }
  for(i = 0; i<npart5prop; i++){
    partpropdata *partpropi;
    int j;

    partpropi = part5propinfo+i;
    for(j = 0; j<npartclassinfo; j++){
      partclassdata *partclassj;
      propdata *propi;
      int jj;

      if(partpropi->class_present[j]==0)continue;
      partclassj = partclassinfo+j;
      propi = partclassj->prop;
      if(propi==NULL)continue;
      for(jj = 0; jj<propi->nsmokeview_ids; jj++){
        sv_object *objectj;

        objectj = propi->smv_objects[jj];
        objectj->used_by_device = 1;
      }
    }
  }
}

/* ----------------------- init_avatar ----------------------------- */

void init_avatar(void){
  int iavatar_types_local;
  sv_object *objecti,*object_start;
  char com_buffer[1024];
  char labels[1024];

  strcpy(labels,":DUM1 :DUM2 :DUM3 :W :D :H1 :SX :SY :SZ :R :G :B :HX :HY :HZ ");

  object_start = object_def_first.next;
  navatar_types=2;
  for(objecti = object_start;objecti->next!=NULL;objecti=objecti->next){
    if(objecti->type==IS_AVATAR)navatar_types++;
  }
  NewMemory((void **)&avatar_types,navatar_types*sizeof(sv_object *));

  strcpy(com_buffer,labels);
  strcat(com_buffer,"0.0 0.0 1.0 translate 255 0 0 setrgb 0.03 0.1 drawdisk 0 0 255 setrgb 90.0 rotatey 0.03 0.2 drawdisk");
  avatar_defs_backup[0] = init_SVOBJECT1("Avatar_1", com_buffer,1);
  avatar_defs_backup[0]->type=IS_AVATAR;

  strcpy(com_buffer,labels);
  strcat(com_buffer,"255 255 0 setrgb 0.02 0.05 drawdisk");
  avatar_defs_backup[1] = init_SVOBJECT1("Avatar_2", com_buffer,1);
  avatar_defs_backup[1]->type=IS_AVATAR;

  avatar_types[0]=avatar_defs_backup[0];
  avatar_types[1]=avatar_defs_backup[1];

  iavatar_types_local=2;
  for(objecti = object_start;objecti->next!=NULL;objecti=objecti->next){
    if(objecti->type==IS_NOT_AVATAR)continue;
    avatar_types[iavatar_types_local++]=objecti;
  }
  iavatar_types_local=0;
}

/* ----------------------- dist ----------------------------- */

float dist(float p1[3], float p2[3]){
  float dx, dy, dz;

  dx = p1[0] - p2[0];
  dy = p1[1] - p2[1];
  dz = p1[2] - p2[2];
  return sqrt(dx*dx+dy*dy+dz*dz);
}

/* ----------------------- get_point2box_dist ----------------------------- */

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
  ASSERT(FFALSE);
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

void init_device_plane(devicedata *devicei){
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
    meshdata *meshi;
    float xvert[12], yvert[12], zvert[12];
    int triangles[18];
    int nvert, ntriangles;
    int nodeindexes[8], closestnodes[18];
    float vals[8];

    InitIsosurface(devicei->plane_surface[i],level,devicei->color,colorindex);
    devicei->plane_surface[i]->cullfaces=1;

    meshi = meshinfo + i;

    xx[0]=meshi->xyz_bar0[XXX];
    xx[1]=DENORMALIZE_X(meshi->xyz_bar[XXX]);

    yy[0]=meshi->xyz_bar0[YYY];
    yy[1]=DENORMALIZE_Y(meshi->xyz_bar[YYY]);

    zz[0]=meshi->xyz_bar0[ZZZ];
    zz[1]=DENORMALIZE_Z(meshi->xyz_bar[ZZZ]);

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

    xx[0]=NORMALIZE_X(meshi->xyz_bar0[XXX]);
    yy[0]=NORMALIZE_Y(meshi->xyz_bar0[YYY]);
    zz[0]=NORMALIZE_Z(meshi->xyz_bar0[ZZZ]);
    xx[1]=meshi->xyz_bar[XXX];
    yy[1]=meshi->xyz_bar[YYY];
    zz[1]=meshi->xyz_bar[ZZZ];

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

void init_device(devicedata *devicei, float *xyz, float *xyzn, int state0, int nparams, float *params, char *labelptr){
  float norm;
  int i;

  devicei->nvals=0;
  devicei->filetype=-1;
  devicei->in_zone_csv=0;
  devicei->in_devc_csv=0;
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
  devicei->times=NULL;
  devicei->vals=NULL;
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
  devicei->ival=0;
  if(nparams>0&&params!=NULL){
    for(i=0;i<nparams;i++){
      devicei->params[i]=params[i];
    }
  }
}

/* ----------------------- get_indep_var_indices ----------------------------- */

void get_indep_var_indices(sv_object *smv_object,
        char **var_indep_strings, int nvars_indep, int *index){

  int i;
  sv_object_frame *obj_frame;

  obj_frame=smv_object->obj_frames[0];

  for(i=0;i<nvars_indep;i++){
    char *var;

    var = var_indep_strings[i];
    index[i]=get_token_loc(var,obj_frame);
  }
}

/* ----------------------- get_evac_indices ----------------------------- */

void get_evac_indices(sv_object *smv_object,int *evac_index,int *nevac_index){

  int n;

  sv_object_frame *obj_frame;

  obj_frame=smv_object->obj_frames[0];

  n=0;

  evac_index[n++]=get_token_loc("W",obj_frame);
  evac_index[n++]=get_token_loc("D",obj_frame);
  evac_index[n++]=get_token_loc("H1",obj_frame);
  evac_index[n++]=get_token_loc("SX",obj_frame);
  evac_index[n++]=get_token_loc("SY",obj_frame);
  evac_index[n++]=get_token_loc("SZ",obj_frame);
  evac_index[n++]=get_token_loc("R",obj_frame);
  evac_index[n++]=get_token_loc("G",obj_frame);
  evac_index[n++]=get_token_loc("B",obj_frame);
  evac_index[n++]=get_token_loc("HX",obj_frame);
  evac_index[n++]=get_token_loc("HY",obj_frame);
  evac_index[n++]=get_token_loc("HZ",obj_frame);

  *nevac_index=n;
}

/* ----------------------- update_partclass_depend ----------------------------- */

void update_partclass_depend(partclassdata *partclassi){
  int i;

  if(partclassi->prop!=NULL){
    sv_object_frame *obj_frame;
    int nvar;

    if(partclassi->kind==HUMANS){
      partclassi->prop->draw_evac=1;
    }
    else{
      partclassi->prop->draw_evac=0;
    }
    obj_frame=partclassi->prop->smv_object->obj_frames[0];
    for(i=0;i<partclassi->nvars_dep-3;i++){
      char *var;

      var=partclassi->vars_dep[i];
      partclassi->vars_dep_index[i]=get_token_loc(var,obj_frame);
    }
    nvar = partclassi->nvars_dep;
    partclassi->vars_dep_index[nvar-3]=get_token_loc("R",obj_frame);
    partclassi->vars_dep_index[nvar-2]=get_token_loc("G",obj_frame);
    partclassi->vars_dep_index[nvar-1]=get_token_loc("B",obj_frame);
  }
}

/* ------------------ normalize ------------------------ */

void normalize(float *xyz, int n){
  float norm,norm2;
  int i;

  norm2 = 0.0;

  for(i=0;i<n;i++){
    norm2 += xyz[i]*xyz[i];
  }
  norm = sqrt(norm2);
  if(norm<0.00001){
    for(i=0;i<n-1;i++){
      xyz[i]=0.0;
    }
    xyz[n-1]=1.0;
  }
  else{
    for(i=0;i<n;i++){
      xyz[i]/=norm;
    }
  }
}

/* ------------------ parse_object_string ------------------------ */

void parse_object_string(char *string,char **tokens, int *ntokens){
  int i, len, in_quote, in_token, last_in_token,ntok2=0;
  char *c;
  char *tokens_head[BUFFER_SIZE], *tokens_tail[BUFFER_SIZE];
  int in_head=1,nhead=0,ntail=0;

  c=string;
  in_quote=0;
  in_token=0;
  last_in_token=0;
  len=strlen(string);
  for(i=0;i<=len;i++){
    switch(*c){
      case '"':
        in_quote=1-in_quote;
        in_token=1;
        break;
      case ' ':
        if(in_quote==0){
          in_token=0;
        }
        break;
      case 0:
        in_token=0;
        break;
      default:
        in_token=1;
        break;
    }
    if(in_token>last_in_token){
      if(in_head==1&&c[0]==':'){
        tokens_head[nhead++]=c;
      }
      else{
        tokens_tail[ntail++]=c;
        in_head=0;
      }
    }
    if(in_token<last_in_token){
      char *tok;
      int in_head2;

      *c=0;
      if(ntail>0)tok = tokens_tail[ntail-1];
      if(ntail>0&&(strcmp(tok,"include")==0||strcmp(tok,"includef")==0)){
        int j;
        sv_object *included_object;
        int iframe_local;
	      char *object_name;
	      int nparms;
	      sv_object_frame *frame;
        int len2;

        object_name=tokens_tail[ntail-2];
        if(object_name[0]=='"')object_name++;
        len2=strlen(object_name);
        if(object_name[len2-1]=='"')object_name[len2-1]=0;

        if(missing_device==NULL){
          char com_buffer[1024];

          strcpy(com_buffer, "0 0 255 setrgb push 45.0 rotatey -0.1 offsetz 0.05 0.2 drawdisk pop push -45.0 rotatey -0.1 offsetz 0.05 0.2 drawdisk pop");
          missing_device = init_SVOBJECT1("missing_device", com_buffer,1);
        }

        included_object = get_SVOBJECT_type2(object_name,missing_device);

        if(strcmp(tok,"includef")==0&&included_object!=missing_device&&ntail>2){
          char *iframe_label;

          iframe_label=tokens_tail[ntail-3];
          sscanf(iframe_label,"%i",&iframe_local);
          if(iframe_local<0)iframe_local=0;
          if(iframe_local>included_object->nframes-1)iframe_local=included_object->nframes-1;
          nparms=3;
        }
        else{
          iframe_local=0;
          nparms=2;
        }
        ntail-=nparms;
        for(j=0,frame=included_object->first_frame.next;frame->next!=NULL;j++,frame=frame->next){
          if(j==iframe_local)break;
        }
        in_head2=1;
        for(j=0;j<frame->ntokens;j++){
          char *cc;

          cc = frame->tokens[j].tokenlabel;
          if(in_head2==1&&cc[0]==':'){
            tokens_head[nhead++]=frame->tokens[j].tokenfulllabel;
          }
          else{
            in_head2=0;
            tokens_tail[ntail++]=cc;
          }
        }
      }
    }
    last_in_token=in_token;
    c++;
  }
  ntok2=0;
  for(i=0;i<nhead;i++){
    tokens[ntok2++]=tokens_head[i];
  }
  for(i=0;i<ntail;i++){
    tokens[ntok2++]=tokens_tail[i];
  }
  *ntokens=ntok2;
}


