#define CPP
#include "options.h"

#include <stdio.h>
#include <string.h>
#include GLUT_H
#include <math.h>

#include "smokeviewvars.h"

GLUI *glui_stereo=NULL;

GLUI_Panel *PANEL_stereo_method=NULL;

GLUI_RadioGroup *RADIO_stereotype=NULL;
GLUI_RadioGroup *RADIO_stereotype_frame=NULL;

GLUI_RadioButton *RADIOBUTTON_seq=NULL;
GLUI_RadioButton *RADIOBUTTON_1=NULL;
GLUI_RadioButton *RADIOBUTTON_2=NULL;
GLUI_RadioButton *RADIOBUTTON_3=NULL;
GLUI_RadioButton *RADIOBUTTON_4=NULL;
GLUI_RadioButton *RADIOBUTTON_5=NULL;

GLUI_Spinner *SPINNER_zero_parallax=NULL, *SPINNER_right_green2=NULL, *SPINNER_right_blue2=NULL;

GLUI_Button *BUTTON_stereo_1=NULL;
GLUI_Button *BUTTON_stereo_2=NULL;
GLUI_Button *BUTTON_stereo_3=NULL;

#define STEREO_CLOSE 0
#define STEREO_RESET 2
#define STEREO_SHOW 4
#define STEREO_GREEN 5
#define STEREO_BLUE 6
#define SAVE_SETTINGS 999

void STEREO_CB(int var);

/* ------------------ Update_Glui_Stereo ------------------------ */

extern "C" void Update_Glui_Stereo(void){
  if(RADIO_stereotype!=NULL){
    RADIO_stereotype->set_int_val(stereotype);
  }
  if(stereotypeOLD==3&&stereotype!=3){
    if(setbw!=setbwSAVE){
      setbw=1-setbwSAVE;
      ColorbarMenu(COLORBAR_TOGGLE_BW);
    }
  }
  else if(stereotypeOLD!=3&&stereotype==3){
    setbwSAVE=setbw;
    if(setbw==0){
      setbwSAVE=setbw;
      ColorbarMenu(COLORBAR_TOGGLE_BW);
    }
  }
}

/* ------------------ glui_stereo_setup ------------------------ */

extern "C" void glui_stereo_setup(int main_window){
  update_glui_stereo=0;
  if(glui_stereo!=NULL){
    glui_stereo->close();
    glui_stereo=NULL;
  }
  if(glui_stereo!=NULL)glui_stereo->close();
  glui_stereo = GLUI_Master.create_glui("Stereo",0,0,0);
  glui_stereo->hide();

  PANEL_stereo_method = glui_stereo->add_panel(_d("Stereo Method"));
  RADIO_stereotype = glui_stereo->add_radiogroup_to_panel(PANEL_stereo_method,&stereotype,STEREO_SHOW,STEREO_CB);
  RADIOBUTTON_1=glui_stereo->add_radiobutton_to_group(RADIO_stereotype,_d("Off"));
  RADIOBUTTON_seq=glui_stereo->add_radiobutton_to_group(RADIO_stereotype,_d("Successive frames"));
  if(videoSTEREO==0)RADIOBUTTON_seq->disable();
  RADIOBUTTON_2=glui_stereo->add_radiobutton_to_group(RADIO_stereotype,_d("Left/Right"));
  RADIOBUTTON_3=glui_stereo->add_radiobutton_to_group(RADIO_stereotype,_d("Red/Blue"));
  RADIOBUTTON_4=glui_stereo->add_radiobutton_to_group(RADIO_stereotype,_d("Red/Cyan"));
  RADIOBUTTON_5=glui_stereo->add_radiobutton_to_group(RADIO_stereotype,_d("Custom Red/Custom Blue"));
  SPINNER_right_green2=glui_stereo->add_spinner_to_panel(PANEL_stereo_method,_d("green"),GLUI_SPINNER_FLOAT,&right_green,STEREO_GREEN,STEREO_CB);
  SPINNER_right_blue2= glui_stereo->add_spinner_to_panel(PANEL_stereo_method, _d("blue"),GLUI_SPINNER_FLOAT,&right_blue,STEREO_BLUE,STEREO_CB);

  SPINNER_right_green2->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);
  SPINNER_right_blue2->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);

  fzero=SCALE2FDS(fzero);
  SPINNER_zero_parallax=glui_stereo->add_spinner(_d("Distance to zero parallax plane (m)"),GLUI_SPINNER_FLOAT,&fzero);
  glui_stereo->add_checkbox("Show stereo parallax",&show_parallax);
  RADIO_stereotype_frame = glui_stereo->add_radiogroup(&stereotype_frame);
  glui_stereo->add_radiobutton_to_group(RADIO_stereotype_frame,_d("Left eye"));
  glui_stereo->add_radiobutton_to_group(RADIO_stereotype_frame,_d("Right eye"));
  glui_stereo->add_radiobutton_to_group(RADIO_stereotype_frame,_d("Both eyes"));
  //SPINNER_zero_parallax->set_float_limits(0.1*xyzmaxdiff,2.0*xyzmaxdiff,GLUI_LIMIT_CLAMP);
  STEREO_CB(STEREO_SHOW);
  Update_Glui_Stereo();

  BUTTON_stereo_1=glui_stereo->add_button(_d("Reset"),STEREO_RESET,STEREO_CB);
  BUTTON_stereo_2=glui_stereo->add_button(_d("Save settings"),SAVE_SETTINGS,STEREO_CB);
  BUTTON_stereo_3=glui_stereo->add_button(_d("Close"),STEREO_CLOSE,STEREO_CB);

  glui_stereo->set_main_gfx_window( main_window );
}

/* ------------------ hide_glui_stereo ------------------------ */

extern "C" void hide_glui_stereo(void){
  if(glui_stereo!=NULL)glui_stereo->hide();
  updatemenu=1;
}

/* ------------------ show_glui_stereo ------------------------ */

extern "C" void show_glui_stereo(void){
  if(glui_stereo!=NULL)glui_stereo->show();
}

/* ------------------ STEREO_CB ------------------------ */

void STEREO_CB(int var){

  switch(var){
  case STEREO_GREEN:
   // right_blue=1.0-right_green;
   // SPINNER_right_blue2->set_float_val(right_blue);
    break;
  case STEREO_BLUE:
   // right_green=1.0-right_blue;
   // SPINNER_right_green2->set_float_val(right_green);
    break;
  case STEREO_SHOW:
    if(stereotypeOLD!=stereotype){
      Update_Glui_Stereo();
      stereotypeOLD=stereotype;
    }
    if(stereotype==STEREO_CUSTOM){
      SPINNER_right_blue2->enable();
      SPINNER_right_green2->enable();
    }
    else{
      SPINNER_right_blue2->disable();
      SPINNER_right_green2->disable();
    }
    break;
  case STEREO_RESET:
    SPINNER_zero_parallax->set_float_val(SCALE2FDS(0.25));
    break;
  case STEREO_CLOSE:
    hide_glui_stereo();
    break;
  case SAVE_SETTINGS:
    writeini(LOCAL_INI,NULL);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
}

