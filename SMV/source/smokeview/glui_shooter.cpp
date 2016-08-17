#define CPP
#include "options.h"

#include <stdio.h>
#include <string.h>
#include GLUT_H
#include <math.h>

#include "update.h"
#include "smokeviewvars.h"

GLUI *glui_shooter=NULL;

GLUI_Rollout *ROLLOUT_shooter_frame=NULL;
GLUI_Rollout *ROLLOUT_shooter_misc=NULL;
GLUI_Rollout *ROLLOUT_shooter_velocity = NULL;

GLUI_Panel *PANEL_shooter_frameA=NULL;
GLUI_Panel *PANEL_shooter_frameB=NULL;
GLUI_Panel *PANEL_shooter_frameC=NULL;
GLUI_Panel *PANEL_shooter_frameD=NULL;
GLUI_Panel *PANEL_shooter_frameE=NULL;
GLUI_Panel *PANEL_shooter_frameF=NULL;
GLUI_Panel *PANEL_shooter_frameG=NULL;
GLUI_Panel *PANEL_shooter_frameH=NULL;
GLUI_Panel *PANEL_shooter_win = NULL;

GLUI_RadioGroup *RADIO_shooter_vel_type=NULL;

GLUI_RadioButton *RADIOBUTTON_plot3dtype=NULL;
GLUI_RadioButton *RADIOBUTTON_shooter_1=NULL;

GLUI_Spinner *SPINNER_shooter_x=NULL;
GLUI_Spinner *SPINNER_shooter_y=NULL;
GLUI_Spinner *SPINNER_shooter_z=NULL;
GLUI_Spinner *SPINNER_shooter_u=NULL;
GLUI_Spinner *SPINNER_shooter_v=NULL;
GLUI_Spinner *SPINNER_shooter_w=NULL;
GLUI_Spinner *SPINNER_shooter_dx=NULL;
GLUI_Spinner *SPINNER_shooter_dy=NULL;
GLUI_Spinner *SPINNER_shooter_dz=NULL;
GLUI_Spinner *SPINNER_shooter_nparts=NULL;
GLUI_Spinner *SPINNER_shooter_fps=NULL;
GLUI_Spinner *SPINNER_shooter_veldir=NULL;
GLUI_Spinner *SPINNER_shooter_v_inf=NULL;
GLUI_Spinner *SPINNER_shooter_u0=NULL;
GLUI_Spinner *SPINNER_shooter_z0=NULL;
GLUI_Spinner *SPINNER_shooter_p=NULL;
GLUI_Spinner *SPINNER_shooter_duration=NULL;
GLUI_Spinner *SPINNER_shooter_history=NULL;
GLUI_Spinner *SPINNER_shooter_1=NULL;

GLUI_Button *BUTTON_shooter_loadplot3d=NULL;
GLUI_Button *BUTTON_shooter_1=NULL;
GLUI_Button *BUTTON_shooter_2=NULL;
GLUI_Button *BUTTON_shooter_3=NULL;

GLUI_Checkbox *CHECKBOX_shooter_1=NULL;
GLUI_Checkbox *CHECKBOX_shooter_2=NULL;
GLUI_Checkbox *CHECKBOX_shooter_3=NULL;

GLUI_Listbox *LIST_shooter_times=NULL;


#define SHOOTER_VEL_TYPE 101
#define SHOOTER_APPLY 102
#define SHOOTER_DURATION 103
#define SHOOTER_FPS 104
#define SHOOTER_NPARTS 105
#define SHOOTER_HISTORY 106
#define SHOOTER_XYZ 107
#define SHOOTER_DXYZ 108
#define SHOOTER_VEL 109
#define SHOOTER_SHOW 110
#define SHOOTER_FIRSTFRAME 111
#define SHOOTER_TIME 112
#define SHOOTER_LOADPLOT3D 113
#define SHOOTER_UVW 114
#define SHOOTER_TERMINAL_VEL 115

#define SAVE_SETTINGS 900
#define SHOOTER_CLOSE 901

void SHOOTER_CB(int var);
extern "C" int allocate_shooter(void);
extern "C" void init_shooter_data(void);
extern "C" void Plot3DListMenu(int val);

#define START_SHOOTER_ROLLOUT 0
#define BACKGROUND_SHOOTER_ROLLOUT 1
#define MISC_SHOOTER_ROLLOUT 2

procdata shooterprocinfo[3];
int nshooterprocinfo = 0;

/* ------------------ Shooter_Rollout_CB ------------------------ */

void Shooter_Rollout_CB(int var){
  toggle_rollout(shooterprocinfo, nshooterprocinfo, var);
}

/* ------------------ hide_glui_shooter ------------------------ */

extern "C" void hide_glui_shooter(void){
  if(glui_shooter!=NULL){
    glui_shooter->hide();
    updatemenu=1;
  }
}

/* ------------------ show_glui_shooter ------------------------ */

extern "C" void show_glui_shooter(void){
  if(glui_shooter!=NULL){
    glui_shooter->show();
    updatemenu=1;
  }
}

/* ------------------ glui_shooter_setup ------------------------ */

extern "C" void glui_shooter_setup(int main_window){

  update_glui_shooter=0;
  if(glui_shooter!=NULL){
    glui_shooter->close();
    glui_shooter=NULL;
  }
  glui_shooter = GLUI_Master.create_glui(_d("Particle tracking"),0,0,0 );
  glui_shooter->hide();

  ROLLOUT_shooter_frame = glui_shooter->add_rollout(_d("Starting locations/velocities"), true, START_SHOOTER_ROLLOUT, Shooter_Rollout_CB);
  ADDPROCINFO(shooterprocinfo, nshooterprocinfo, ROLLOUT_shooter_frame, START_SHOOTER_ROLLOUT);

  PANEL_shooter_frameE=glui_shooter->add_panel_to_panel(ROLLOUT_shooter_frame,_d("Positions"),false);

  PANEL_shooter_frameA=glui_shooter->add_panel_to_panel(PANEL_shooter_frameE,_d("Center"));
  glui_shooter->add_column_to_panel(PANEL_shooter_frameE,false);
  PANEL_shooter_frameB=glui_shooter->add_panel_to_panel(PANEL_shooter_frameE,_d("Size"));

  SPINNER_shooter_x=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameA,"x",GLUI_SPINNER_FLOAT,shooter_xyz,SHOOTER_XYZ,SHOOTER_CB);
  SPINNER_shooter_x->set_float_limits(xbar0,xbarORIG);

  SPINNER_shooter_y=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameA,"y",GLUI_SPINNER_FLOAT,shooter_xyz+1,SHOOTER_XYZ,SHOOTER_CB);
  SPINNER_shooter_y->set_float_limits(ybar0,ybarORIG);

  SPINNER_shooter_z=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameA,"z",GLUI_SPINNER_FLOAT,shooter_xyz+2,SHOOTER_XYZ,SHOOTER_CB);
  SPINNER_shooter_z->set_float_limits(zbar0,zbarORIG);

  SPINNER_shooter_dx=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameB,"dx",GLUI_SPINNER_FLOAT,shooter_dxyz,SHOOTER_DXYZ,SHOOTER_CB);
  SPINNER_shooter_dx->set_float_limits(0.0,xbarORIG-xbar0);

  SPINNER_shooter_dy=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameB,"dy",GLUI_SPINNER_FLOAT,shooter_dxyz+1,SHOOTER_DXYZ,SHOOTER_CB);
  SPINNER_shooter_dy->set_float_limits(0.0,ybarORIG-ybar0);

  SPINNER_shooter_dz=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameB,"dz",GLUI_SPINNER_FLOAT,shooter_dxyz+2,SHOOTER_DXYZ,SHOOTER_CB);
  SPINNER_shooter_dz->set_float_limits(0.0,zbarORIG-zbar0);

  PANEL_shooter_frameF=glui_shooter->add_panel_to_panel(ROLLOUT_shooter_frame,_d("Velocities"));
  SPINNER_shooter_u=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameF,"u",GLUI_SPINNER_FLOAT,shooter_uvw,SHOOTER_UVW,SHOOTER_CB);
  SPINNER_shooter_v=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameF,"v",GLUI_SPINNER_FLOAT,shooter_uvw+1,SHOOTER_UVW,SHOOTER_CB);
  SPINNER_shooter_w=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameF,"w",GLUI_SPINNER_FLOAT,shooter_uvw+2,SHOOTER_UVW,SHOOTER_CB);

  ROLLOUT_shooter_velocity = glui_shooter->add_rollout(_d("Background velocity field"), false, BACKGROUND_SHOOTER_ROLLOUT, Shooter_Rollout_CB);
  ADDPROCINFO(shooterprocinfo, nshooterprocinfo, ROLLOUT_shooter_velocity,BACKGROUND_SHOOTER_ROLLOUT);

  RADIO_shooter_vel_type=glui_shooter->add_radiogroup_to_panel(ROLLOUT_shooter_velocity,&shooter_vel_type,SHOOTER_VEL_TYPE,SHOOTER_CB);
  RADIOBUTTON_plot3dtype=glui_shooter->add_radiobutton_to_group(RADIO_shooter_vel_type,"PLOT3D");
  RADIOBUTTON_shooter_1=glui_shooter->add_radiobutton_to_group(RADIO_shooter_vel_type,_d("Power law"));
  if(nplot3dtimelist>0&&plot3dtimelist!=NULL){
  }
  else{
    shooter_vel_type=1;
    RADIO_shooter_vel_type->set_int_val(shooter_vel_type);
    RADIOBUTTON_plot3dtype->disable();
  }

  if(nplot3dtimelist>0&&plot3dtimelist!=NULL){
    int i;

    PANEL_shooter_frameC=glui_shooter->add_panel_to_panel(ROLLOUT_shooter_velocity,"PLOT3D");
    BUTTON_shooter_loadplot3d=glui_shooter->add_button_to_panel(PANEL_shooter_frameC,_d("Load"),SHOOTER_LOADPLOT3D,SHOOTER_CB);
    LIST_shooter_times = glui_shooter->add_listbox_to_panel(PANEL_shooter_frameC,_d("Time:"),&shooter_itime,SHOOTER_TIME,SHOOTER_CB);
    for(i=0;i<nplot3dtimelist;i++){
      char label[255];

      sprintf(label,"%f",plot3dtimelist[i]);
      TrimZeros(label);
      LIST_shooter_times->add_item(i,label);
    }
  }

  PANEL_shooter_frameD=glui_shooter->add_panel_to_panel(ROLLOUT_shooter_velocity,_d("Power law"));
  SPINNER_shooter_u0=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameD,_d("reference velocity, U0 (m/s)"),GLUI_SPINNER_FLOAT,&shooter_u0,SHOOTER_VEL,SHOOTER_CB);
  SPINNER_shooter_z0=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameD,_d("reference elevation, Z0 (m)"),GLUI_SPINNER_FLOAT,&shooter_z0,SHOOTER_VEL,SHOOTER_CB);
  SPINNER_shooter_p=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameD,_d("decay, p"),GLUI_SPINNER_FLOAT,&shooter_p,SHOOTER_VEL,SHOOTER_CB);
  SPINNER_shooter_veldir=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameD,_d("velocity direction (deg)"),GLUI_SPINNER_FLOAT,&shooter_veldir,SHOOTER_VEL,SHOOTER_CB);
  SPINNER_shooter_veldir->set_float_limits(-180.0,180.0);

  ROLLOUT_shooter_misc = glui_shooter->add_rollout("Misc", false, MISC_SHOOTER_ROLLOUT, Shooter_Rollout_CB);
  ADDPROCINFO(shooterprocinfo, nshooterprocinfo, ROLLOUT_shooter_misc, MISC_SHOOTER_ROLLOUT);

  PANEL_shooter_frameG = glui_shooter->add_panel_to_panel(ROLLOUT_shooter_misc, "", false);
  glui_shooter->add_column_to_panel(ROLLOUT_shooter_misc,false);
  PANEL_shooter_frameH=glui_shooter->add_panel_to_panel(ROLLOUT_shooter_misc,"",false);

  CHECKBOX_shooter_1=glui_shooter->add_checkbox_to_panel(PANEL_shooter_frameG,_d("Show particles"),&visShooter,SHOOTER_SHOW,SHOOTER_CB);
  CHECKBOX_shooter_2=glui_shooter->add_checkbox_to_panel(PANEL_shooter_frameG,_d("Update continuously"),&shooter_cont_update);
  CHECKBOX_shooter_3=glui_shooter->add_checkbox_to_panel(PANEL_shooter_frameG,_d("Show only first frame"),&shooter_firstframe,SHOOTER_FIRSTFRAME,SHOOTER_CB);
  SPINNER_shooter_v_inf=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameG,_d("terminal velocity"),GLUI_SPINNER_FLOAT,&shooter_v_inf,SHOOTER_TERMINAL_VEL,SHOOTER_CB);
  BUTTON_shooter_1=glui_shooter->add_button_to_panel(PANEL_shooter_frameG,_d("Compute tracks"),SHOOTER_APPLY,SHOOTER_CB);

  SPINNER_shooter_1=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameH,_d("Particle size"),GLUI_SPINNER_FLOAT,&shooterpointsize);
  SPINNER_shooter_nparts=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameH,_d("number of particles"),GLUI_SPINNER_INT,&shooter_nparts,SHOOTER_NPARTS,SHOOTER_CB);

  SPINNER_shooter_fps=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameH,_d("frames per second"),GLUI_SPINNER_INT,&shooter_fps,SHOOTER_FPS,SHOOTER_CB);

  SPINNER_shooter_duration=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameH,_d("duration (s)"),GLUI_SPINNER_FLOAT,&shooter_duration,SHOOTER_DURATION,SHOOTER_CB);
 // SPINNER_shooter_history=glui_shooter->add_spinner_to_panel(PANEL_shooter_frameH,"history (s)",
 //   GLUI_SPINNER_FLOAT,&shooter_history,SHOOTER_HISTORY,SHOOTER_CB);
 // SPINNER_shooter_history->disable();

  SHOOTER_CB(SHOOTER_NPARTS);
  SHOOTER_CB(SHOOTER_FPS);
  SHOOTER_CB(SHOOTER_DURATION);

  PANEL_shooter_win=glui_shooter->add_panel("",GLUI_PANEL_NONE);

  BUTTON_shooter_2=glui_shooter->add_button_to_panel(PANEL_shooter_win,_d("Save settings"),SAVE_SETTINGS,SHOOTER_CB);
  glui_shooter->add_column_to_panel(PANEL_shooter_win,false);
  BUTTON_shooter_3=glui_shooter->add_button_to_panel(PANEL_shooter_win,_d("Close"),SHOOTER_CLOSE,SHOOTER_CB);

  SHOOTER_CB(SHOOTER_VEL_TYPE);
  SHOOTER_CB(SHOOTER_VEL);

  glui_shooter->set_main_gfx_window( main_window );
}

/* ------------------ SHOOTER_CB ------------------------ */

void SHOOTER_CB(int var){
  float ang;
  if(shooter_firstframe==1){
    ResetItimes0();
  }
  switch(var){
    case SHOOTER_LOADPLOT3D:
      PRINTF("Loading PLOT3D data at time: %f\n",plot3dtimelist[shooter_itime]);
      Plot3DListMenu(shooter_itime);
      SHOOTER_CB(SHOOTER_APPLY);
      break;
    case SHOOTER_TIME:
      break;
    case SHOOTER_FIRSTFRAME:
      break;
    case SHOOTER_SHOW:
      plotstate=GetPlotState(DYNAMIC_PLOTS);
      UpdateTimes();
      break;
    case SHOOTER_TERMINAL_VEL:
      if(shooter_v_inf<0.0){
        shooter_v_inf=0.0;
        SPINNER_shooter_v_inf->set_float_val(shooter_v_inf);
      }
      if(shooter_cont_update==1){
        SHOOTER_CB(SHOOTER_APPLY);
      }
      break;
    case SHOOTER_VEL:
      ang = DEG2RAD*shooter_veldir;
      shooter_velz=0.0;
      shooter_velx = shooter_u0*cos(ang);
      shooter_vely = shooter_u0*sin(ang);
      if(shooter_cont_update==1){
        SHOOTER_CB(SHOOTER_APPLY);
      }
      break;
    case SHOOTER_UVW:
    case SHOOTER_XYZ:
      if(shooter_active==1){
        init_shooter_data();
      }
      if(shooter_cont_update==1){
        SHOOTER_CB(SHOOTER_APPLY);
      }
      break;
    case SHOOTER_DXYZ:
      if(shooter_active==1){
        init_shooter_data();
      }
      if(shooter_cont_update==1){
        SHOOTER_CB(SHOOTER_APPLY);
      }
      break;
    case SHOOTER_NPARTS:
      if(shooter_nparts<1&&SPINNER_shooter_nparts!=NULL){
        shooter_nparts=1;
        SPINNER_shooter_nparts->set_int_val(shooter_nparts);
      }
      if(shooter_cont_update==1){
        SHOOTER_CB(SHOOTER_APPLY);
      }
      break;
    case SHOOTER_FPS:
      if(shooter_fps<1&&SPINNER_shooter_fps!=NULL){
        shooter_fps=1;
        SPINNER_shooter_fps->set_int_val(shooter_fps);
      }
      if(shooter_active==1){
        init_shooter_data();
      }
      if(shooter_cont_update==1){
        SHOOTER_CB(SHOOTER_APPLY);
      }
      break;
    case SHOOTER_HISTORY:
      break;
    case SHOOTER_DURATION:
      if(shooter_duration<1.0&&SPINNER_shooter_duration!=NULL){
        shooter_duration=1.0;
        SPINNER_shooter_duration->set_float_val(shooter_duration);
      }
      if(shooter_cont_update==1){
        SHOOTER_CB(SHOOTER_APPLY);
      }
      break;
    case SHOOTER_APPLY:
      nshooter_frames=(int)(shooter_duration*shooter_fps);
      max_shooter_points=nshooter_frames*shooter_nparts;

      if(allocate_shooter()==0){
        solve_shooter_data();
        plotstate=GetPlotState(DYNAMIC_PLOTS);
        UpdateTimes();
      }
      break;
    case SHOOTER_VEL_TYPE:
      if(shooter_vel_type==1){
        SPINNER_shooter_u0->enable();
        SPINNER_shooter_z0->enable();
        SPINNER_shooter_p->enable();
        SPINNER_shooter_veldir->enable();
        if(LIST_shooter_times!=NULL)LIST_shooter_times->disable();
        if(BUTTON_shooter_loadplot3d!=NULL)BUTTON_shooter_loadplot3d->disable();
      }
      else{
        SPINNER_shooter_u0->disable();
        SPINNER_shooter_z0->disable();
        SPINNER_shooter_p->disable();
        SPINNER_shooter_veldir->disable();
        if(LIST_shooter_times!=NULL)LIST_shooter_times->enable();
        if(BUTTON_shooter_loadplot3d!=NULL)BUTTON_shooter_loadplot3d->enable();
      }
      if(shooter_cont_update==1){
        SHOOTER_CB(SHOOTER_APPLY);
      }
      break;
    case SAVE_SETTINGS:
      WriteINI(LOCAL_INI,NULL);
      break;
    case SHOOTER_CLOSE:
      hide_glui_shooter();
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
}
