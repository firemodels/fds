// $Date$ 
// $Revision$
// $Author$

#include "options.h"
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

// svn revision character string
extern "C" char glui_stereo_revision[]="$Revision$";

GLUI_Panel *panel_stereo_method=NULL;
GLUI_RadioGroup *RADIO_showstereo=NULL;
GLUI_RadioButton *RADIO_seq=NULL;
GLUI *glui_stereo=NULL;
GLUI_Spinner *SPINNER_zero_parallax=NULL, *SPINNER_stereo_offset=NULL;


#define STEREO_CLOSE 0
#define STEREO_RESET 2
#define STEREO_FRAME 3
#define SAVE_SETTINGS 999

void STEREO_CB(int var);

extern "C" void update_glui_stereo(void){
  int showstereoOLD;

  if(RADIO_showstereo!=NULL){
    RADIO_showstereo->set_int_val(showstereo);
  }
  showstereoOLD=showstereo;
  if(showstereoOLD==3&&showstereo!=3){
    if(setbw!=setbwSAVE){
      setbw=1-setbwSAVE;
      ShadeMenu(2);
    }
  }
  else if(showstereoOLD!=3&&showstereo==3){
    setbwSAVE=setbw;
    if(setbw==0){
      setbwSAVE=setbw;
      ShadeMenu(2);
    }
  }
}

/* ------------------ glui_stereo_setup ------------------------ */

extern "C" void glui_stereo_setup(int main_window){
  if(glui_stereo!=NULL)glui_stereo->close();
  glui_stereo = GLUI_Master.create_glui("stereo",0,0,0);
  if(showgluistereo==0)glui_stereo->hide();
  
  panel_stereo_method = glui_stereo->add_panel("Stereo Method");
  RADIO_showstereo = glui_stereo->add_radiogroup_to_panel(panel_stereo_method,&showstereo);
  glui_stereo->add_radiobutton_to_group(RADIO_showstereo,"Off");
  RADIO_seq=glui_stereo->add_radiobutton_to_group(RADIO_showstereo,"Sucessive frames");
  if(videoSTEREO==0)RADIO_seq->disable();
  glui_stereo->add_radiobutton_to_group(RADIO_showstereo,"Left/Right");
  glui_stereo->add_radiobutton_to_group(RADIO_showstereo,"Red/Blue");

  fzero*=xyzmaxdiff;
  SPINNER_zero_parallax=glui_stereo->add_spinner("Distance to zero parallax plane (m)",GLUI_SPINNER_FLOAT,&fzero);
  //SPINNER_zero_parallax->set_float_limits(0.1*xyzmaxdiff,2.0*xyzmaxdiff,GLUI_LIMIT_CLAMP);

  update_glui_stereo();

  glui_stereo->add_button("Reset",STEREO_RESET,STEREO_CB);
  glui_stereo->add_button("Save Settings",SAVE_SETTINGS,STEREO_CB);
  glui_stereo->add_button("Close",STEREO_CLOSE,STEREO_CB);
  
  glui_stereo->set_main_gfx_window( main_window );
}

/* ------------------ hide_glui_stereo ------------------------ */

extern "C" void hide_glui_stereo(void){
  if(glui_stereo!=NULL)glui_stereo->hide();
  showgluistereo=0;
  updatemenu=1;
}

/* ------------------ show_glui_stereo ------------------------ */

extern "C" void show_glui_stereo(void){
  if(glui_stereo!=NULL)glui_stereo->show();
}

/* ------------------ STEREO_CB ------------------------ */

void STEREO_CB(int var){

  switch (var){
  case STEREO_RESET:
    SPINNER_zero_parallax->set_float_val(0.25*xyzmaxdiff);
    break;
  case STEREO_CLOSE:
    hide_glui_stereo();
    break;
  case SAVE_SETTINGS:
    writeini(LOCAL_INI);
    break;
  }
}

