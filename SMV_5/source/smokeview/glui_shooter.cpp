// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#ifdef pp_SHOOTER
#include <string.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <math.h>
#include "glui.h"
#include "flowfiles.h"
#define CPP
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

GLUI *glui_shooter=NULL;


// svn revision character string
extern "C" char glui_shooter_revision[]="$Revision$";
// $Date$ $Author$

GLUI_Panel *panel_shooter_frame=NULL;
GLUI_Panel *panel_shooter_frameA=NULL;
GLUI_Panel *panel_shooter_frameB=NULL;
GLUI_Panel *panel_shooter_velocity=NULL;
GLUI_Panel *panel_shooter_win=NULL;
GLUI_RadioGroup *RADIO_shooter_vel_type=NULL;
GLUI_Spinner *SPINNER_shooter_x=NULL;
GLUI_Spinner *SPINNER_shooter_y=NULL;
GLUI_Spinner *SPINNER_shooter_z=NULL;
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

#define SHOOTER_VEL_TYPE 101
#define SHOOTER_APPLY 102
#define SHOOTER_DURATION 103
#define SHOOTER_FPS 104
#define SHOOTER_NPARTS 105
#define SHOOTER_HISTORY 106
#define SHOOTER_XYZ 107
#define SHOOTER_DXYZ 108
#define SHOOTER_VEL 109

#define SAVE_SETTINGS 900
#define SHOOTER_CLOSE 901

void SHOOTER_CB(int var);
extern "C" int allocate_shooter(void);
extern "C" void init_shooter_data(void);

/* ------------------ hide_glui_shooter ------------------------ */

extern "C" void hide_shooter(void){
  if(glui_shooter!=NULL){
    glui_shooter->hide();
    showshooter=0;
    updatemenu=1;
  }
}

/* ------------------ show_shooter ------------------------ */

extern "C" void show_shooter(void){
  if(glui_shooter!=NULL){
    glui_shooter->show();
    showshooter=1;
    updatemenu=1;
  }
}

/* ------------------ glui_shooter_setup ------------------------ */

extern "C" void glui_shooter_setup(int main_window){  

  glui_shooter = GLUI_Master.create_glui( "Particle Shooting",0,0,0 );
  if(showshooter==0)glui_shooter->hide();


  panel_shooter_frame=glui_shooter->add_panel("Initial Frame");

  panel_shooter_frameA=glui_shooter->add_panel_to_panel(panel_shooter_frame,"Center");
  glui_shooter->add_column_to_panel(panel_shooter_frame,false);
  panel_shooter_frameB=glui_shooter->add_panel_to_panel(panel_shooter_frame,"Size");

  SPINNER_shooter_x=glui_shooter->add_spinner_to_panel(panel_shooter_frameA,"x",
    GLUI_SPINNER_FLOAT,shooter_xyz,SHOOTER_XYZ,SHOOTER_CB);
  SPINNER_shooter_x->set_float_limits(xbar0,xbarORIG);

  SPINNER_shooter_y=glui_shooter->add_spinner_to_panel(panel_shooter_frameA,"y",
    GLUI_SPINNER_FLOAT,shooter_xyz+1,SHOOTER_XYZ,SHOOTER_CB);
  SPINNER_shooter_y->set_float_limits(ybar0,ybarORIG);

  SPINNER_shooter_z=glui_shooter->add_spinner_to_panel(panel_shooter_frameA,"z",
    GLUI_SPINNER_FLOAT,shooter_xyz+2,SHOOTER_XYZ,SHOOTER_CB);
  SPINNER_shooter_z->set_float_limits(zbar0,zbarORIG);

  SPINNER_shooter_dx=glui_shooter->add_spinner_to_panel(panel_shooter_frameB,"dx",
    GLUI_SPINNER_FLOAT,shooter_dxyz,SHOOTER_DXYZ,SHOOTER_CB);
  SPINNER_shooter_dx->set_float_limits(0.0,xbarORIG-xbar0);

  SPINNER_shooter_dy=glui_shooter->add_spinner_to_panel(panel_shooter_frameB,"dy",
    GLUI_SPINNER_FLOAT,shooter_dxyz+1,SHOOTER_DXYZ,SHOOTER_CB);
  SPINNER_shooter_dy->set_float_limits(0.0,ybarORIG-ybar0);

  SPINNER_shooter_dz=glui_shooter->add_spinner_to_panel(panel_shooter_frameB,"dz",
    GLUI_SPINNER_FLOAT,shooter_dxyz+2,SHOOTER_DXYZ,SHOOTER_CB);
  SPINNER_shooter_dz->set_float_limits(0.0,zbarORIG-zbar0);

  glui_shooter->add_spinner_to_panel(panel_shooter_frameA,"size",GLUI_SPINNER_FLOAT,&shooterpointsize);
  SPINNER_shooter_nparts=glui_shooter->add_spinner_to_panel(panel_shooter_frame,"number of particles",
    GLUI_SPINNER_INT,&shooter_nparts,SHOOTER_NPARTS,SHOOTER_CB);

  SPINNER_shooter_fps=glui_shooter->add_spinner_to_panel(panel_shooter_frame,"frames per second",
    GLUI_SPINNER_INT,&shooter_fps,SHOOTER_FPS,SHOOTER_CB);

  SPINNER_shooter_duration=glui_shooter->add_spinner_to_panel(panel_shooter_frame,"duration (s)",
    GLUI_SPINNER_FLOAT,&shooter_duration,SHOOTER_DURATION,SHOOTER_CB);
  SPINNER_shooter_history=glui_shooter->add_spinner_to_panel(panel_shooter_frame,"history (s)",
    GLUI_SPINNER_FLOAT,&shooter_history,SHOOTER_HISTORY,SHOOTER_CB);
  glui_shooter->add_checkbox_to_panel(panel_shooter_frame,"Show particles",&show_shooter_points);

  SHOOTER_CB(SHOOTER_NPARTS);
  SHOOTER_CB(SHOOTER_FPS);
  SHOOTER_CB(SHOOTER_DURATION);

  panel_shooter_velocity=glui_shooter->add_panel("Initial Velocity");
  
  RADIO_shooter_vel_type=glui_shooter->add_radiogroup_to_panel(panel_shooter_velocity,&shooter_vel_type,
    SHOOTER_VEL_TYPE,SHOOTER_CB);
  glui_shooter->add_radiobutton_to_group(RADIO_shooter_vel_type,"Use PLOT3D velocity data");
  glui_shooter->add_radiobutton_to_group(RADIO_shooter_vel_type,"Use specified Profile");

  SPINNER_shooter_v_inf=glui_shooter->add_spinner_to_panel(panel_shooter_velocity,"terminal velocity, v_inf (m/s)",
    GLUI_SPINNER_FLOAT,&shooter_v_inf,SHOOTER_VEL,SHOOTER_CB);
  SPINNER_shooter_u0=glui_shooter->add_spinner_to_panel(panel_shooter_velocity,"reference velocity, U0 (m/s)",
    GLUI_SPINNER_FLOAT,&shooter_u0,SHOOTER_VEL,SHOOTER_CB);
  SPINNER_shooter_z0=glui_shooter->add_spinner_to_panel(panel_shooter_velocity,"reference elevation, Z0 (m)",
    GLUI_SPINNER_FLOAT,&shooter_z0,SHOOTER_VEL,SHOOTER_CB);
  SPINNER_shooter_p=glui_shooter->add_spinner_to_panel(panel_shooter_velocity,"decay, p",
    GLUI_SPINNER_FLOAT,&shooter_p,SHOOTER_VEL,SHOOTER_CB);
  SPINNER_shooter_veldir=glui_shooter->add_spinner_to_panel(panel_shooter_velocity,"velocity direction (deg)",
    GLUI_SPINNER_FLOAT,&shooter_veldir,SHOOTER_VEL,SHOOTER_CB);
  SPINNER_shooter_veldir->set_float_limits(-180.0,180.0);

  panel_shooter_win=glui_shooter->add_panel("",GLUI_PANEL_NONE);

  glui_shooter->add_button_to_panel(panel_shooter_win,"Apply",SHOOTER_APPLY,SHOOTER_CB);
  glui_shooter->add_column_to_panel(panel_shooter_win,false);
  glui_shooter->add_button_to_panel(panel_shooter_win,"Save Settings",SAVE_SETTINGS,SHOOTER_CB);
  glui_shooter->add_column_to_panel(panel_shooter_win,false);
  glui_shooter->add_button_to_panel(panel_shooter_win,"Close",SHOOTER_CLOSE,SHOOTER_CB);

  SHOOTER_CB(SHOOTER_VEL_TYPE);
  SHOOTER_CB(SHOOTER_VEL);

  glui_shooter->set_main_gfx_window( main_window );
}

/* ------------------ SHOOTER_CB ------------------------ */

void SHOOTER_CB(int var){
  float pi,ang;

  switch (var){
    case SHOOTER_VEL:
      pi = 4.0*atan(1.0);
      ang = 2.0*pi*shooter_veldir/360.0;
      shooter_velz=0.0;
      shooter_velx = shooter_u0*cos(ang);
      shooter_vely = shooter_u0*sin(ang);
      break;
    case SHOOTER_XYZ:
      if(shooter_active==1){
        init_shooter_data();
      }
      break;
    case SHOOTER_DXYZ:
      if(shooter_active==1){
        init_shooter_data();
      }
      break;
    case SHOOTER_NPARTS:
      if(shooter_nparts<1&&SPINNER_shooter_nparts!=NULL){
        shooter_nparts=1;
        SPINNER_shooter_nparts->set_int_val(shooter_nparts);
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
      break;
    case SHOOTER_HISTORY:
      break;
    case SHOOTER_DURATION:
      if(shooter_duration<1.0&&SPINNER_shooter_duration!=NULL){
        shooter_duration=1.0;
        SPINNER_shooter_duration->set_float_val(shooter_duration);
      }
      break;
    case SHOOTER_APPLY:
      max_shooter_frames=shooter_duration*shooter_fps;
      max_shooter_frames*=(shooter_history*shooter_fps);
      max_shooter_points=max_shooter_frames*shooter_nparts;

      if(allocate_shooter()==0){
        init_shooter_data();
      }
      break;
    case SHOOTER_VEL_TYPE:
      if(shooter_vel_type==1){
        SPINNER_shooter_v_inf->enable();
        SPINNER_shooter_u0->enable();
        SPINNER_shooter_z0->enable();
        SPINNER_shooter_p->enable();
        SPINNER_shooter_veldir->enable();
      }
      else{
        SPINNER_shooter_v_inf->disable();
        SPINNER_shooter_u0->disable();
        SPINNER_shooter_z0->disable();
        SPINNER_shooter_p->disable();
        SPINNER_shooter_veldir->disable();
     }
      break;
    case SAVE_SETTINGS:
      writeini(LOCAL_INI);
      break;
    case SHOOTER_CLOSE:
      hide_shooter();
      break;
  }
}

#endif
