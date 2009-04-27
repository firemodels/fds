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
GLUI_Spinner *SPINNER_shooter_velmag=NULL;

#define SHOOTER_VEL_TYPE 101

#define SAVE_SETTINGS 900
#define SHOOTER_CLOSE 901

void SHOOTER_CB(int var);

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

  SPINNER_shooter_x=glui_shooter->add_spinner_to_panel(panel_shooter_frameA,"x",GLUI_SPINNER_FLOAT,shooter_xyz);
  SPINNER_shooter_y=glui_shooter->add_spinner_to_panel(panel_shooter_frameA,"y",GLUI_SPINNER_FLOAT,shooter_xyz+1);
  SPINNER_shooter_z=glui_shooter->add_spinner_to_panel(panel_shooter_frameA,"z",GLUI_SPINNER_FLOAT,shooter_xyz+2);
  SPINNER_shooter_dx=glui_shooter->add_spinner_to_panel(panel_shooter_frameB,"dx",GLUI_SPINNER_FLOAT,shooter_dxyz);
  SPINNER_shooter_dy=glui_shooter->add_spinner_to_panel(panel_shooter_frameB,"dy",GLUI_SPINNER_FLOAT,shooter_dxyz+1);
  SPINNER_shooter_dz=glui_shooter->add_spinner_to_panel(panel_shooter_frameB,"dz",GLUI_SPINNER_FLOAT,shooter_dxyz+2);

  SPINNER_shooter_nparts=glui_shooter->add_spinner_to_panel(panel_shooter_frame,"number of particles",
    GLUI_SPINNER_INT,&shooter_nparts);
  SPINNER_shooter_nparts->set_w(1200);
  SPINNER_shooter_nparts->set_int_limits(1,100);
  SPINNER_shooter_fps=glui_shooter->add_spinner_to_panel(panel_shooter_frame,"frames per second",
    GLUI_SPINNER_INT,&shooter_fps);
  SPINNER_shooter_fps->set_int_limits(1,100);

  panel_shooter_velocity=glui_shooter->add_panel("Initial Velocity");
  
  RADIO_shooter_vel_type=glui_shooter->add_radiogroup_to_panel(panel_shooter_velocity,&shooter_vel_type,
    SHOOTER_VEL_TYPE,SHOOTER_CB);
  glui_shooter->add_radiobutton_to_group(RADIO_shooter_vel_type,"Use PLOT3D velocity data");
  glui_shooter->add_radiobutton_to_group(RADIO_shooter_vel_type,"Use specified Profile");

  SPINNER_shooter_velmag=glui_shooter->add_spinner_to_panel(panel_shooter_velocity,"velocity magnitude (m/s)",
    GLUI_SPINNER_FLOAT,&shooter_velmag);
  SPINNER_shooter_veldir=glui_shooter->add_spinner_to_panel(panel_shooter_velocity,"velocity direction (deg)",
    GLUI_SPINNER_FLOAT,&shooter_veldir);
  SPINNER_shooter_veldir->set_float_limits(-180.0,180.0);

  panel_shooter_win=glui_shooter->add_panel("",GLUI_PANEL_NONE);

  glui_shooter->add_button_to_panel(panel_shooter_win,"Save Settings",SAVE_SETTINGS,SHOOTER_CB);
  glui_shooter->add_column_to_panel(panel_shooter_win,false);
  glui_shooter->add_button_to_panel(panel_shooter_win,"Close",SHOOTER_CLOSE,SHOOTER_CB);

  SHOOTER_CB(SHOOTER_VEL_TYPE);

  glui_shooter->set_main_gfx_window( main_window );
}

/* ------------------ SHOOTER_CB ------------------------ */

void SHOOTER_CB(int var){
  switch (var){
    case SHOOTER_VEL_TYPE:
      if(shooter_vel_type==1){
        SPINNER_shooter_velmag->enable();
        SPINNER_shooter_veldir->enable();
      }
      else{
        SPINNER_shooter_velmag->disable();
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
