// $Date: 2015-03-18 14:50:48 -0400 (Wed, 18 Mar 2015) $ 
// $Revision: 22025 $
// $Author: gforney $

#define CPP
#include "options.h"

// svn revision character string
extern "C" char glui_clip_revision[];
char glui_clip_revision[]="$Revision: 22025 $";

#include <stdio.h>
#include <string.h>
#include GLUT_H
#include <math.h>

#include "smokeviewvars.h"

GLUI *glui_clip=NULL;

GLUI_RadioGroup *radio_clip=NULL;

GLUI_Spinner *SPINNER_clip_xmax=NULL, *SPINNER_clip_xmin=NULL;
GLUI_Spinner *SPINNER_clip_ymax=NULL, *SPINNER_clip_ymin=NULL;
GLUI_Spinner *SPINNER_clip_zmax=NULL, *SPINNER_clip_zmin=NULL;

GLUI_Checkbox *CHECKBOX_clip_xmin=NULL, *CHECKBOX_clip_xmax=NULL;
GLUI_Checkbox *CHECKBOX_clip_ymin=NULL, *CHECKBOX_clip_ymax=NULL;
GLUI_Checkbox *CHECKBOX_clip_zmin=NULL, *CHECKBOX_clip_zmax=NULL;

GLUI_Panel *PANEL_clip_lower=NULL, *PANEL_clip_upper=NULL, *PANEL_clip=NULL,*panel_wrapup=NULL;
GLUI_Panel *PANEL_clipx=NULL, *PANEL_clipX=NULL;
GLUI_Panel *PANEL_clipy=NULL, *PANEL_clipY=NULL;
GLUI_Panel *PANEL_clipz=NULL, *PANEL_clipZ=NULL;
GLUI_Panel *PANEL_blockageview=NULL;

GLUI_Listbox *LIST_mesh=NULL;

GLUI_RadioButton *RADIOBUTTON_clip_1a=NULL;
GLUI_RadioButton *RADIOBUTTON_clip_1b=NULL;
GLUI_RadioButton *RADIOBUTTON_clip_1c=NULL;

GLUI_Button *BUTTON_clip_1=NULL;
GLUI_Button *BUTTON_clip_2=NULL;

#define CLIP_xlower 0
#define CLIP_ylower 1
#define CLIP_zlower 2
#define CLIP_xupper 3
#define CLIP_yupper 4
#define CLIP_zupper 5

#define CLIP_all 12

#define SPINNER_xlower 13
#define SPINNER_ylower 14
#define SPINNER_zlower 15
#define SPINNER_xupper 16
#define SPINNER_yupper 17
#define SPINNER_zupper 18

#define INI_VALS -1
#define DEFAULT_VALS -2

#define CLIP_CLOSE 99
#define SAVE_SETTINGS 98
#define CLIP_MESH 80



void CLIP_CB(int var);
void set_clip_controls(int val);

/* ------------------ glui_clip_setup ------------------------ */

extern "C" void glui_clip_setup(int main_window){
  int i;

  update_glui_clip=0;
  if(glui_clip!=NULL){
    glui_clip->close();
    glui_clip=NULL;
  }
  glui_clip = GLUI_Master.create_glui("Clipping",0,0,0);
  glui_clip->hide();

  PANEL_clip = glui_clip->add_panel("",GLUI_PANEL_NONE);
  PANEL_clip_lower = glui_clip->add_panel_to_panel(PANEL_clip,_("Clip Lower"));
  PANEL_clipx = glui_clip->add_panel_to_panel(PANEL_clip_lower,"X",GLUI_PANEL_NONE);
  SPINNER_clip_xmin=glui_clip->add_spinner_to_panel(PANEL_clipx,"X",GLUI_SPINNER_FLOAT,&clipinfo.xmin,SPINNER_xlower,CLIP_CB);
  SPINNER_clip_xmin->set_float_limits(xclip_min,xclip_max,GLUI_LIMIT_CLAMP);
  glui_clip->add_column_to_panel(PANEL_clipx,false);
  CHECKBOX_clip_xmin=glui_clip->add_checkbox_to_panel(PANEL_clipx,"",&clipinfo.clip_xmin,CLIP_xlower,CLIP_CB);

  PANEL_clipy = glui_clip->add_panel_to_panel(PANEL_clip_lower,"Y",GLUI_PANEL_NONE);
  SPINNER_clip_ymin=glui_clip->add_spinner_to_panel(PANEL_clipy,"Y",GLUI_SPINNER_FLOAT,&clipinfo.ymin,SPINNER_ylower,CLIP_CB);
  SPINNER_clip_ymin->set_float_limits(yclip_min,yclip_max,GLUI_LIMIT_CLAMP);
  glui_clip->add_column_to_panel(PANEL_clipy,false);
  CHECKBOX_clip_ymin=glui_clip->add_checkbox_to_panel(PANEL_clipy,"",&clipinfo.clip_ymin,CLIP_ylower,CLIP_CB);

  PANEL_clipz = glui_clip->add_panel_to_panel(PANEL_clip_lower,"Z",GLUI_PANEL_NONE);
  SPINNER_clip_zmin=glui_clip->add_spinner_to_panel(PANEL_clipz,"Z",GLUI_SPINNER_FLOAT,&clipinfo.zmin,SPINNER_zlower,CLIP_CB);
  SPINNER_clip_zmin->set_float_limits(zclip_min,zclip_max,GLUI_LIMIT_CLAMP);
  glui_clip->add_column_to_panel(PANEL_clipz,false);
  CHECKBOX_clip_zmin=glui_clip->add_checkbox_to_panel(PANEL_clipz,"",&clipinfo.clip_zmin,CLIP_zlower,CLIP_CB);

  radio_clip = glui_clip->add_radiogroup_to_panel(PANEL_clip,&clip_mode,CLIP_all,CLIP_CB);
  RADIOBUTTON_clip_1a=glui_clip->add_radiobutton_to_group(radio_clip,_("Clipping disabled"));
  RADIOBUTTON_clip_1b=glui_clip->add_radiobutton_to_group(radio_clip,_("Clip blockages and data"));
  RADIOBUTTON_clip_1c=glui_clip->add_radiobutton_to_group(radio_clip,_("Clip blockages"));
  RADIOBUTTON_clip_1c=glui_clip->add_radiobutton_to_group(radio_clip,_("Clip data"));

  glui_clip->add_column_to_panel(PANEL_clip,false);

  PANEL_clip_upper = glui_clip->add_panel_to_panel(PANEL_clip,_("Clip upper"));

  PANEL_clipX = glui_clip->add_panel_to_panel(PANEL_clip_upper,"X",GLUI_PANEL_NONE);
  SPINNER_clip_xmax=glui_clip->add_spinner_to_panel(PANEL_clipX,"X",GLUI_SPINNER_FLOAT,&clipinfo.xmax,SPINNER_xupper,CLIP_CB);
  SPINNER_clip_xmax->set_float_limits(xclip_min,xclip_max,GLUI_LIMIT_CLAMP);
  glui_clip->add_column_to_panel(PANEL_clipX,false);
  CHECKBOX_clip_xmax=glui_clip->add_checkbox_to_panel(PANEL_clipX,"",&clipinfo.clip_xmax,CLIP_xupper,CLIP_CB);

  PANEL_clipY = glui_clip->add_panel_to_panel(PANEL_clip_upper,"Y",GLUI_PANEL_NONE);
  SPINNER_clip_ymax=glui_clip->add_spinner_to_panel(PANEL_clipY,"Y",GLUI_SPINNER_FLOAT,&clipinfo.ymax,SPINNER_yupper,CLIP_CB);
  SPINNER_clip_ymax->set_float_limits(yclip_min,yclip_max,GLUI_LIMIT_CLAMP);
  glui_clip->add_column_to_panel(PANEL_clipY,false);
  CHECKBOX_clip_ymax=glui_clip->add_checkbox_to_panel(PANEL_clipY,"",&clipinfo.clip_ymax,CLIP_yupper,CLIP_CB);

  PANEL_clipZ = glui_clip->add_panel_to_panel(PANEL_clip_upper,"Z",GLUI_PANEL_NONE);
  SPINNER_clip_zmax=glui_clip->add_spinner_to_panel(PANEL_clipZ,"Z",GLUI_SPINNER_FLOAT,&clipinfo.zmax,SPINNER_zupper,CLIP_CB);
  SPINNER_clip_zmax->set_float_limits(zclip_min,zclip_max,GLUI_LIMIT_CLAMP);
  glui_clip->add_column_to_panel(PANEL_clipZ,false);
  CHECKBOX_clip_zmax=glui_clip->add_checkbox_to_panel(PANEL_clipZ,"",&clipinfo.clip_zmax,CLIP_zupper,CLIP_CB);

  PANEL_blockageview = glui_clip->add_rollout_to_panel(PANEL_clip,"Hide blockages",false);
  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi = meshinfo + i;
    glui_clip->add_checkbox_to_panel(PANEL_blockageview,meshi->label,&meshi->blockvis);
  }

  panel_wrapup = glui_clip->add_panel_to_panel(PANEL_clip,"",GLUI_PANEL_NONE);

  glui_clip->add_column_to_panel(panel_wrapup,false);

  BUTTON_clip_1=glui_clip->add_button_to_panel(panel_wrapup,_("Save settings"),SAVE_SETTINGS,CLIP_CB);

  glui_clip->add_column_to_panel(panel_wrapup,false);

  BUTTON_clip_2=glui_clip->add_button_to_panel(panel_wrapup,_("Close"),CLIP_CLOSE,CLIP_CB);

  if(updateclipvals==1){
    set_clip_controls(INI_VALS);  // clip vals from ini file
  }
  else{
    if(clip_mesh==0){
      set_clip_controls(DEFAULT_VALS);  // clip vals from global scene
    }
    else{
      set_clip_controls(clip_mesh);  // clip vals from mesh clip_mesh
    }
  }

  glui_clip->set_main_gfx_window( main_window );
}

/* ------------------ hide_glui_clip ------------------------ */

extern "C" void hide_glui_clip(void){
  if(glui_clip!=NULL)glui_clip->hide();
  updatemenu=1;
}

/* ------------------ show_glui_clip ------------------------ */

extern "C" void show_glui_clip(void){
  if(glui_clip!=NULL)glui_clip->show();
}


/* ------------------ Update_Glui_Clip ------------------------ */

extern "C" void Update_Glui_Clip(void){
  if(CHECKBOX_clip_xmin!=NULL&&CHECKBOX_clip_ymin!=NULL&&CHECKBOX_clip_zmin!=NULL&&
     CHECKBOX_clip_xmax!=NULL&&CHECKBOX_clip_ymax!=NULL&&CHECKBOX_clip_zmax!=NULL){

    CHECKBOX_clip_xmin->set_int_val(clipinfo.clip_xmin);
    CHECKBOX_clip_ymin->set_int_val(clipinfo.clip_ymin);
    CHECKBOX_clip_zmin->set_int_val(clipinfo.clip_zmin);
    CHECKBOX_clip_xmax->set_int_val(clipinfo.clip_xmax);
    CHECKBOX_clip_ymax->set_int_val(clipinfo.clip_ymax);
    CHECKBOX_clip_zmax->set_int_val(clipinfo.clip_zmax);
    if(radio_clip!=NULL)radio_clip->set_int_val(clip_mode);
    CLIP_CB(CLIP_all);
  }
}

/* ------------------ CLIP_CB ------------------------ */

void CLIP_CB(int var){
  int i;

  glutPostRedisplay();
  switch(var){
  case CLIP_MESH:
    if(clip_mesh==0){
      set_clip_controls(DEFAULT_VALS);
    }
    else{
      set_clip_controls(clip_mesh);
    }
    break;
  case SAVE_SETTINGS:
    writeini(LOCAL_INI,NULL);
    break;
  case CLIP_CLOSE:
    hide_glui_clip();
    break;
  case CLIP_xlower:
    if(clipinfo.clip_xmin==0)SPINNER_clip_xmin->disable();
    if(clipinfo.clip_xmin==1)SPINNER_clip_xmin->enable();
    updatefacelists=1;
    break;
  case CLIP_ylower:
    if(clipinfo.clip_ymin==0)SPINNER_clip_ymin->disable();
    if(clipinfo.clip_ymin==1)SPINNER_clip_ymin->enable();
    updatefacelists=1;
    break;
  case CLIP_zlower:
    if(clipinfo.clip_zmin==0)SPINNER_clip_zmin->disable();
    if(clipinfo.clip_zmin==1)SPINNER_clip_zmin->enable();
    updatefacelists=1;
    break;
  case CLIP_xupper:
    if(clipinfo.clip_xmax==0)SPINNER_clip_xmax->disable();
    if(clipinfo.clip_xmax==1)SPINNER_clip_xmax->enable();
    updatefacelists=1;
    break;
  case CLIP_yupper:
    if(clipinfo.clip_ymax==0)SPINNER_clip_ymax->disable();
    if(clipinfo.clip_ymax==1)SPINNER_clip_ymax->enable();
    updatefacelists=1;
    break;
  case CLIP_zupper:
    if(clipinfo.clip_zmax==0)SPINNER_clip_zmax->disable();
    if(clipinfo.clip_zmax==1)SPINNER_clip_zmax->enable();
    updatefacelists=1;
    break;
  case CLIP_all:
    updatefacelists=1;
    update_clipplanes();
    if(clip_mode!=CLIP_OFF){
      for(i=0;i<6;i++){
        CLIP_CB(i);
      }
      CHECKBOX_clip_xmin->enable();
      CHECKBOX_clip_ymin->enable();
      CHECKBOX_clip_zmin->enable();
      CHECKBOX_clip_xmax->enable();
      CHECKBOX_clip_ymax->enable();
      CHECKBOX_clip_zmax->enable();
    }
    else{
      SPINNER_clip_xmin->disable();
      SPINNER_clip_ymin->disable();
      SPINNER_clip_zmin->disable();
      SPINNER_clip_xmax->disable();
      SPINNER_clip_ymax->disable();
      SPINNER_clip_zmax->disable();

      CHECKBOX_clip_xmin->disable();
      CHECKBOX_clip_ymin->disable();
      CHECKBOX_clip_zmin->disable();
      CHECKBOX_clip_xmax->disable();
      CHECKBOX_clip_ymax->disable();
      CHECKBOX_clip_zmax->disable();
    }
    break;
  case SPINNER_xlower:
    SPINNER_clip_xmax->set_float_limits(clipinfo.xmin,xclip_max,GLUI_LIMIT_CLAMP);
    updatefacelists=1;
    break;
  case SPINNER_xupper:
    SPINNER_clip_xmin->set_float_limits(xclip_min,clipinfo.xmax,GLUI_LIMIT_CLAMP);
    updatefacelists=1;
    break;
  case SPINNER_ylower:
    SPINNER_clip_ymax->set_float_limits(clipinfo.ymin,yclip_max,GLUI_LIMIT_CLAMP);
    updatefacelists=1;
    break;
  case SPINNER_yupper:
    SPINNER_clip_ymin->set_float_limits(yclip_min,clipinfo.ymax,GLUI_LIMIT_CLAMP);
    updatefacelists=1;
    break;
  case SPINNER_zlower:
    SPINNER_clip_zmax->set_float_limits(clipinfo.zmin,zclip_max,GLUI_LIMIT_CLAMP);
    updatefacelists=1;
    break;
  case SPINNER_zupper:
    SPINNER_clip_zmin->set_float_limits(zclip_min,clipinfo.zmax,GLUI_LIMIT_CLAMP);
    updatefacelists=1;
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  if(var>=CLIP_xlower&&var<=CLIP_zupper){
    clip2cam(camera_current);
  }
}

/* ------------------ update_clip_all ------------------------ */

extern "C" void update_clip_all(void){
  CLIP_CB(CLIP_all);
  radio_clip->set_int_val(clip_mode);
}

/* ------------------ set_clip_controls ------------------------ */

void set_clip_controls(int val){
  int i;

  for(i=0;i<6;i++){
    CLIP_CB(i);
  }
  if(val==DEFAULT_VALS){
    clipinfo.xmin = xclip_min;
    clipinfo.ymin = yclip_min;
    clipinfo.zmin = zclip_min;
    clipinfo.xmax = xclip_max;
    clipinfo.ymax = yclip_max;
    clipinfo.zmax = zclip_max;
  }
  if(val>=1&&val<=nmeshes){
    mesh *meshi;
    float *xplt, *yplt, *zplt;

    float dxclip, dyclip, dzclip;

    dxclip = (xbarORIG-xbar0ORIG)/1000.0;
    dyclip = (ybarORIG-ybar0ORIG)/1000.0;
    dzclip = (zbarORIG-zbar0ORIG)/1000.0;

    meshi = meshinfo + val - 1;

    xplt = meshi->xplt_orig;
    yplt = meshi->yplt_orig;
    zplt = meshi->zplt_orig;

    clipinfo.xmin = xplt[0] - dxclip;
    clipinfo.ymin = yplt[0] - dyclip;
    clipinfo.zmin = zplt[0] - dzclip;
    clipinfo.xmax = xplt[meshi->ibar] + dxclip;
    clipinfo.ymax = yplt[meshi->jbar] + dyclip;
    clipinfo.zmax = zplt[meshi->kbar] + dzclip;
  }
  SPINNER_clip_xmin->set_float_val(clipinfo.xmin);
  SPINNER_clip_ymin->set_float_val(clipinfo.ymin);
  SPINNER_clip_zmin->set_float_val(clipinfo.zmin);
  SPINNER_clip_xmax->set_float_val(clipinfo.xmax);
  SPINNER_clip_ymax->set_float_val(clipinfo.ymax);
  SPINNER_clip_zmax->set_float_val(clipinfo.zmax);
}

