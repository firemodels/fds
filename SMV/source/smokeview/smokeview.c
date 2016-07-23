#include "options.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "glew.h"
#include GLUT_H

#include "smokeviewvars.h"

#ifdef WIN32
#include <direct.h>
#endif

/* ------------------ _Sniff_Errors ------------------------ */

void _Sniff_Errors(char *whereat){
  int error;

  while((error=glGetError())!=GL_NO_ERROR){
    char *glu_error;

    glu_error=(char *)gluErrorString((unsigned int)error);
    fprintf(stderr,"*** Error: OpenGL error:%s, where:%s %i\n",
      glu_error,whereat,snifferrornumber);
      snifferrornumber++;
  }
}

/* ------------------ updateLights ------------------------ */

void updateLights(float *pos1, float *pos2){
  int i;
  GLfloat ambientlight2[4], diffuselight2[4];
  int lightCount;
  float div;

  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, lightmodel_localviewer == 0? GL_FALSE : GL_TRUE);
  glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, lightmodel_separatespecularcolor == 0? GL_SINGLE_COLOR : GL_SEPARATE_SPECULAR_COLOR);

  lightCount = 0;
  if(light_enabled0){
    ++lightCount;
  }
  if(light_enabled1){
    ++lightCount;
  }

  div = lightCount > 0? 1.0f/(float)lightCount : 1.0f;
  for(i=0;i<3;i++){
    ambientlight2[i]=ambientlight[i]*div;
    diffuselight2[i]=diffuselight[i]*div;
  }
  ambientlight2[3]=1.0;
  diffuselight2[3]=1.0;
  if(light_enabled0){
    glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuselight2);
    glLightfv(GL_LIGHT0,GL_AMBIENT,ambientlight2);
    if(pos1!=NULL)glLightfv(GL_LIGHT0,GL_POSITION,pos1);
    glEnable(GL_LIGHT0);
  }
  else{
    glDisable(GL_LIGHT0);
  }

  if(light_enabled1){
    glLightfv(GL_LIGHT1,GL_DIFFUSE,diffuselight2);
    glLightfv(GL_LIGHT1,GL_AMBIENT,ambientlight2);
    if(pos2!=NULL)glLightfv(GL_LIGHT1,GL_POSITION,pos2);
    glEnable(GL_LIGHT1);
  }
  else{
    glDisable(GL_LIGHT1);
  }

  UpdateLIGHTS=0;
}

/* ------------------ antialias ------------------------ */

void antialias(int flag){
  if(antialiasflag==1){
    if(flag==1){
      glEnable(GL_LINE_SMOOTH);
      glEnable(GL_BLEND);
      glEnable(GL_POINT_SMOOTH);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glHint(GL_LINE_SMOOTH_HINT,GL_DONT_CARE);
    }
    if(flag==0){
      glDisable(GL_LINE_SMOOTH);
      glDisable(GL_POINT_SMOOTH);
      glDisable(GL_BLEND);
    }
  }
}

/* ------------------ transparenton ------------------------ */

void transparenton(void){
  glEnable(GL_BLEND);
  glDepthMask(GL_FALSE);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
}

/* ------------------ transparentoff ------------------------ */

void transparentoff(void){
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
}


/* ------------------ SetViewPoint ------------------------ */

void SetViewPoint(int option){
  in_external=0;
  switch(option){
    int rotation_type_save;
    int projection_type_save;

  case RESTORE_EXTERIOR_VIEW_ZOOM:
    break;
  case RESTORE_EXTERIOR_VIEW:
    in_external=1;
    rotation_type_save = camera_current->rotation_type;
    projection_type_save = camera_current->projection_type;
    CopyCamera(camera_current,camera_external);
    camera_current->rotation_type=rotation_type_save;
    camera_current->projection_type=projection_type_save;
    if(camera_current->projection_type==1){
      camera_current->eye[1]=camera_current->isometric_y;
    }
    break;
  case RESTORE_INTERIOR_VIEW:
    rotation_type_save = camera_current->rotation_type;
    projection_type_save = camera_current->projection_type;
    CopyCamera(camera_current,camera_internal);
    camera_current->rotation_type=rotation_type_save;
    camera_current->projection_type=projection_type_save;
    break;
  case RESTORE_SAVED_VIEW:
    CopyCamera(camera_current,camera_save);
    break;
  case 3:
  case 4:
    ASSERT(FFALSE);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  if(rotation_type==ROTATION_3AXIS){
    float azimuth, elevation,axis[3];
    float quat_temp[4];
    float x, y, z;

    azimuth = camera_current->az_elev[0]*DEG2RAD;
    elevation = camera_current->az_elev[1]*DEG2RAD;

    x = cos(azimuth);
    y = sin(azimuth);
    z = cos(elevation);

    axis[0]=0.0;
    axis[1]=0.0;
    axis[2]=1.0;

    angleaxis2quat(azimuth,axis,quat_temp);

    axis[0]=x;
    axis[1]=y;
    axis[2]=0.0;

    angleaxis2quat(acos(z),axis,quat_general);

    mult_quat(quat_temp,quat_general,quat_general);

    quat2rot(quat_general,quat_rotation);
  }
  if(option==RESTORE_EXTERIOR_VIEW_ZOOM)camera_current->zoom=zooms[zoomindex];
  zoom=camera_current->zoom;
  update_glui_zoom();
}

/* ------------------ init_volrender_script ------------------------ */

void init_volrender_script(char *prefix, char *tour_label, int startframe, int skipframe){
  scriptfiledata *sfd;
  FILE *script_stream;

  if(volrender_scriptname==NULL){
    int len;

    len = strlen(fdsprefix)+strlen("_volrender.ssf")+1;
    NewMemory((void **)&volrender_scriptname,(unsigned int)(len));
    STRCPY(volrender_scriptname,fdsprefix);
    STRCAT(volrender_scriptname,"_volrender.ssf");
  }

  sfd = insert_scriptfile(volrender_scriptname);
  if(sfd!=NULL)default_script=sfd;
  script_stream=fopen(volrender_scriptname,"w");
  if(script_stream!=NULL){
    fprintf(script_stream,"RENDERDIR\n");
    fprintf(script_stream," .\n");
    if(tour_label!=NULL&&strcmp(tour_label,"Manual")!=0){
      fprintf(script_stream,"LOADTOUR\n");
      fprintf(script_stream," %s\n",tour_label);
    }
    fprintf(script_stream,"VOLSMOKERENDERALL\n");
    fprintf(script_stream," %i %i\n",skipframe,startframe);
    fprintf(script_stream," %s\n",prefix);
    runscript=1;
    fclose(script_stream);
  }
}

/* ------------------ display_version_info ------------------------ */

void display_version_info(char *progname){
  PRINTversion(progname);
  if(fds_version!=NULL){
    PRINTF("FDS Build: %s\n",fds_githash);
  }
  if(smokeviewpath!=NULL){
    PRINTF("Smokeview path: %s\n",smokeviewpath);
  }
  if(smokezippath!=NULL){
    PRINTF("Smokezip path: %s\n",smokezippath);
  }
  if(texturedir!=NULL){
    PRINTF("Texture directory path: %s\n",texturedir);
  }
}

