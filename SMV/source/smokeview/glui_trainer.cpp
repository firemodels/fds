// $Date$ 
// $Revision$
// $Author$

#define CPP
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
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

// svn revision character string
extern "C" char glui_trainer_revision[]="$Revision$";

void ROTATE_CB(int var);
void TRAINER_CB(int var);

#define LOAD_SMOKE 100
#define LOAD_TEMP 101
#define LOAD_OXY 102

#define MOVETYPE 200
#define TRAINERPATH 201
#define TOGGLE_VIEW 203
#define TRAINER_PAUSE 212
#define TRAINER_LEFTRIGHT_INOUT 221
#define TRAINER_UPDOWN 223
#define TRAINER_AZ_ELEV 231

#define TRAINERVIEWPOINTS 300
#define TRAINEROUTLINE 301
#define MENU_OUTLINEVIEW -104

GLUI *glui_trainer=NULL;
GLUI_Checkbox *CHECKBOX_pause=NULL;
GLUI_Checkbox *CHECKBOX_outline=NULL;
GLUI_Listbox *LIST_trainerpath=NULL,*LIST_viewpoint=NULL;
GLUI_Panel *panel_smokeview=NULL;
GLUI_Panel *panel_explore=NULL;
GLUI_Panel *panel_manual=NULL;//, *panel_automatic=NULL;
GLUI_Panel *panel_move=NULL,*panel_rotate2=NULL;
GLUI_Button *BUTTON_smoke3d=NULL, *BUTTON_temp=NULL, *BUTTON_oxy=NULL, *BUTTON_unload=NULL;
GLUI_Button *BUTTON_toggle_view=NULL;
GLUI_Translation *TRANS_updown=NULL,*TRANS_leftright_inout=NULL;
GLUI_Translation *TRANS_az_elev=NULL;
GLUI_StaticText *STATIC_alert=NULL;
GLUI *glui_alert=NULL;

extern "C" void update_glui_viewlist(void){
  if(trainer_viewpoints!=-1){
    LIST_viewpoint->set_int_val(-1);
  }
}

extern "C" void show_load_alert(void){
  showalert=1;
  if(glui_alert!=NULL)glui_alert->show();
}
extern "C" void hide_load_alert(void){
  showalert=0;
  if(glui_alert!=NULL)glui_alert->hide();
}
/* ------------------ glui_tour_setup ------------------------ */

extern "C" void glui_alert_setup(int main_window){
  if(glui_alert!=NULL)glui_alert->close();
  glui_alert = GLUI_Master.create_glui("",0,screenWidth/2,screenHeight/2);
  glui_alert->hide();
  STATIC_alert = glui_alert->add_statictext("Loading smoke and fire data");
}

/* ------------------ hide_glui_trainer ------------------------ */

extern "C" void hide_trainer(void){
  if(glui_trainer!=NULL){
    glui_trainer->hide();
    showtrainer=0;
    updatemenu=1;
  }
}

/* ------------------ show_trainer ------------------------ */

extern "C" void show_trainer(void){
  if(glui_trainer!=NULL){
    glui_trainer->show();
    showtrainer=1;
    updatemenu=1;
  }
}

/* ------------------ update_trainer_outline ------------------------ */

extern "C" void update_trainer_outline(void){
  if(visBlocks==visBLOCKOutline){
    trainer_outline=0;
  }
  else{
    trainer_outline=1;
  }
  if(CHECKBOX_outline!=NULL)CHECKBOX_outline->set_int_val(trainer_outline);
}

/* ------------------ update_trainer_moves ------------------------ */

extern "C" void update_trainer_moves(void){
  float *eye_xyz;
  float *az, *elev;

  az = camera_current->angle_zx;
  elev = camera_current->angle_zx+1;

  eye_xyz = camera_current->eye;

  if(TRANS_leftright_inout!=NULL){
    TRANS_leftright_inout->set_x(eye_xyz[0]);
    TRANS_leftright_inout->set_y(eye_xyz[1]);
    TRANS_leftright_inout->set_speed(1.0/(float)screenWidth);
  }

  if(TRANS_updown!=NULL){
    TRANS_updown->set_x(eye_xyz[2]);
    TRANS_updown->set_speed(1.0/(float)screenHeight);
  }

  if(TRANS_az_elev!=NULL){
    TRANS_az_elev->set_x(*az);
    TRANS_az_elev->set_speed(180.0/(float)screenHeight);
    TRANS_az_elev->set_y(*elev);
  }

}

/* ------------------ glui_trainer_setup ------------------------ */

extern "C" void glui_trainer_setup(int main_window){

  if(glui_trainer!=NULL)glui_trainer->close();
  glui_trainer = GLUI_Master.create_glui("Demonstrator",0,screenWidth+12,0);
  if(showgluitrainer==0)glui_trainer->hide();
  
  glui_trainer->set_main_gfx_window( main_window );
  panel_smokeview = glui_trainer->add_panel("Data");
  BUTTON_smoke3d = glui_trainer->add_button_to_panel(panel_smokeview,"Smoke/Fire",LOAD_SMOKE,TRAINER_CB);
  if(AnySmoke(NULL)==0)BUTTON_smoke3d->disable();
  BUTTON_temp = glui_trainer->add_button_to_panel(panel_smokeview,"Temperature",LOAD_TEMP,TRAINER_CB);
  if(AnySlices("TEMPERATURE")==0)BUTTON_temp->disable();
  BUTTON_oxy = glui_trainer->add_button_to_panel(panel_smokeview,"Oxygen",LOAD_OXY,TRAINER_CB);
  if(AnySlices("oxygen")==0&&AnySlices("oxygen VOLUME FRACTION")==0){
    BUTTON_oxy->disable();
  }
  
  panel_explore = glui_trainer->add_panel("Explore",true);

  trainer_path=-1;
  LIST_trainerpath = glui_trainer->add_listbox_to_panel(panel_explore,"Path:",&trainer_path,TRAINERPATH,TRAINER_CB);
  {
    int i;
    LIST_trainerpath->add_item(-1,"Manual");
    LIST_trainerpath->add_item(-2,"-");
    for(i=0;i<ntours;i++){
      tourdata *touri;

      touri = tourinfo + i;
      LIST_trainerpath->add_item(i,touri->menulabel);
    }
  }

  LIST_viewpoint = glui_trainer->add_listbox_to_panel(panel_explore,"Viewpoint:",&trainer_viewpoints,TRAINERVIEWPOINTS,TRAINER_CB);
  {
    camera *ca;

    for(ca=camera_list_first.next;ca->next!=NULL;ca=ca->next){
      char line[256];

      if(strcmp(ca->name,"internal")==0)continue;
      if(strcmp(ca->name,"external")==0){
        strcpy(line,"Default");
        LIST_viewpoint->add_item(ca->view_id,line);
        LIST_viewpoint->add_item(-1,"-");
      }
      else{
        strcpy(line,ca->name);
        LIST_viewpoint->add_item(ca->view_id,line);
      }
      if(ca->view_id>=1)ntrainer_viewpoints++;
    }
  }
  BUTTON_toggle_view = glui_trainer->add_button_to_panel(panel_explore,"Toggle View",TOGGLE_VIEW,TRAINER_CB);
  if(ntrainer_viewpoints<=2)BUTTON_toggle_view->disable();

  CHECKBOX_outline = glui_trainer->add_checkbox_to_panel(panel_explore,"Show walls",&trainer_outline,TRAINEROUTLINE,TRAINER_CB);
  CHECKBOX_pause = glui_trainer->add_checkbox_to_panel(panel_explore,"Pause",&trainer_pause,TRAINER_PAUSE,TRAINER_CB);

  update_trainer_outline();
  panel_move = glui_trainer->add_panel_to_panel(panel_explore,"Move",false);
  TRANS_leftright_inout = glui_trainer->add_translation_to_panel(panel_move,"Horizontal",
    GLUI_TRANSLATION_XY,trainer_xzy,TRAINER_LEFTRIGHT_INOUT,ROTATE_CB);
  glui_trainer->add_column_to_panel(panel_move,false);

  TRANS_updown = glui_trainer->add_translation_to_panel(panel_move,"Vertical",
    GLUI_TRANSLATION_Y,trainer_xzy+2,TRAINER_UPDOWN,ROTATE_CB);
  glui_trainer->add_column_to_panel(panel_move,false);

  TRANS_az_elev = glui_trainer->add_translation_to_panel(panel_move,"Rotate",
    GLUI_TRANSLATION_XY,trainer_ab,TRAINER_AZ_ELEV,ROTATE_CB);

  update_trainer_moves();

  TRAINER_CB(MOVETYPE);
  TRAINER_CB(TRAINERVIEWPOINTS);

}
#define ANGLE_LEFT -1
#define ANGLE_RIGHT 1
#define GO_FORWARD 0
#define GO_BACKWARD 2
#define MOVE_UP 3
#define MOVE_DOWN 4

/* ------------------ ROTATE_CB ------------------------ */

void ROTATE_CB(int var){

  float *eye_xyz, *az, *elev;

  eye_xyz = camera_current->eye;
  az = camera_current->angle_zx;
  elev = camera_current->angle_zx+1;


  if(eyeview!=0){
    eyeview=0;
    handle_eyeview(0);
    ResetView(RESTORE_EXTERIOR_VIEW);
  }

  if(trainer_viewpoints!=-1){
    LIST_viewpoint->set_int_val(-1);
  }
  if(trainer_path!=-1){
    LIST_trainerpath->set_int_val(-1);
    TRAINER_CB(TRAINERPATH);
  }
  switch (var){
  case TRAINER_AZ_ELEV:
    *az = TRANS_az_elev->get_x();
    *elev = -TRANS_az_elev->get_y();
    break;
  case TRAINER_LEFTRIGHT_INOUT:
    eye_xyz[0]=TRANS_leftright_inout->get_x();
    eye_xyz[1]=TRANS_leftright_inout->get_y();
    break;
  case TRAINER_UPDOWN:
    eye_xyz[2]=TRANS_updown->get_y();
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  camera_current->dirty=1;
}

/* ------------------ TRAINER_CB ------------------------ */

void TRAINER_CB(int var){

  switch (var){
    int i;

  case TRAINER_PAUSE:
    stept=trainer_pause;
    keyboard('t',0,0);
    break;
  case TOGGLE_VIEW:
    if(ntrainer_viewpoints<=0)break;
    trainer_viewpoints--;
    if(trainer_viewpoints<2){
      trainer_viewpoints=ntrainer_viewpoints;
    }
    LIST_viewpoint->set_int_val(trainer_viewpoints);
    TRAINER_CB(TRAINERVIEWPOINTS);
  break;
  case TRAINERVIEWPOINTS:
    if(trainer_viewpoints!=-1){
      ResetMenu(trainer_viewpoints);
    }
    if(trainer_path!=-1){
      int viewpoint_save;

      viewpoint_save=trainer_viewpoints;
      LIST_trainerpath->set_int_val(-1);
      TRAINER_CB(TRAINERPATH);
      if(viewpoint_save!=-11){
        LIST_viewpoint->set_int_val(viewpoint_save);
        ResetMenu(viewpoint_save);
      }
    }
    break;
  case TRAINEROUTLINE:
    if(trainer_outline==0){
      visBlocks=visBLOCKAsInput;
    }
    else{
      visBlocks=visBLOCKOutline;
    }
    ResetMenu(MENU_OUTLINEVIEW);
    break;
  case TRAINERPATH:
    TRAINER_CB(MOVETYPE);
    if(trainer_viewpoints!=1){
      LIST_viewpoint->set_int_val(-1);
    }
    switch (trainer_path){
    case -1:
      if(trainer_path_old!=-1){
        trainer_pause=0;
        CHECKBOX_pause->set_int_val(trainer_pause);
        TourMenu(-2);
        eyeview=0;
        handle_eyeview(0);
        from_glui_trainer=1;
        trainee_location=0;
      }
      TRAINER_CB(MOVETYPE);
      break;
    case -2:
      break;
    default:
      if(eyeview!=1){
        eyeview=1;
        handle_eyeview(0);
      }
      for(i=0;i<ntours;i++){
        tourdata *touri;

        touri = tourinfo + i;
        touri->display=0;
      }
      viewtourfrompath=1;
      TourMenu(trainer_path);
    }
    trainer_path_old=trainer_path;
    break;
  case MOVETYPE:
    eyeview=0;
    handle_eyeview(0);
    ResetView(RESTORE_EXTERIOR_VIEW);
    break;
  case LOAD_SMOKE:
    TrainerViewMenu(1);
    break;
  case LOAD_TEMP:
// kind of a hack, having to put in code seg twice, but this is required to get data chopping to work
    if(slicebounds!=NULL&&islicetype!=-1){
      if(setslicechopmin==1||setslicechopmax==1){
       setslicebounds(islicetype);
      }
    }
    TrainerViewMenu(2);
    updatechopcolors();
    if(slicebounds!=NULL&&islicetype!=-1){
      if(setslicechopmin==1||setslicechopmax==1){
       setslicebounds(islicetype);
      }
    }
    TrainerViewMenu(2);
    updatechopcolors();
    colorbarflip=1;
    ColorBarMenu(-2);
    break;
  case LOAD_OXY:
// kind of a hack, having to put in code seg twice, but this is required to get data chopping to work
    if(slicebounds!=NULL&&islicetype!=-1){
      if(setslicechopmin==1||setslicechopmax==1){
       setslicebounds(islicetype);
      }
    }
    TrainerViewMenu(3);
    updatechopcolors();
    if(slicebounds!=NULL&&islicetype!=-1){
      if(setslicechopmin==1||setslicechopmax==1){
       setslicebounds(islicetype);
      }
    }
    TrainerViewMenu(3);
    updatechopcolors();
    colorbarflip=0;
    ColorBarMenu(-2);
    break;
  }
}

