#define CPP
#include "options.h"

#include <stdio.h>
#include <string.h>
#include GLUT_H
#include <math.h>

#include "smokeviewvars.h"

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

GLUI_Panel *PANEL_smokeview=NULL;
GLUI_Panel *PANEL_explore=NULL;
GLUI_Panel *PANEL_manual=NULL;//, *PANEL_automatic=NULL;
GLUI_Panel *PANEL_move=NULL,*PANEL_rotate2=NULL;

GLUI_Button *BUTTON_smoke3d=NULL, *BUTTON_temp=NULL, *BUTTON_oxy=NULL, *BUTTON_unload=NULL;
GLUI_Button *BUTTON_toggle_view=NULL;

GLUI_Translation *TRANSLATE_updown=NULL,*TRANSLATE_leftright_inout=NULL;
GLUI_Translation *TRANSLATE_az_elev=NULL;

GLUI_StaticText *STATIC_alert=NULL;

GLUI *glui_alert=NULL;

/* ------------------ update_glui_viewlist ------------------------ */

extern "C" void update_glui_viewlist(void){
  if(trainer_viewpoints!=-1){
    LIST_viewpoint->set_int_val(-1);
  }
}

/* ------------------ show_glui_alert ------------------------ */

extern "C" void show_glui_alert(void){
  if(glui_alert!=NULL)glui_alert->show();
}

/* ------------------ hide_glui_alert ------------------------ */

extern "C" void hide_glui_alert(void){
  if(glui_alert!=NULL)glui_alert->hide();
}

/* ------------------ glui_alert_setup ------------------------ */

extern "C" void glui_alert_setup(int main_window){
  update_glui_alert=0;
  if(glui_alert!=NULL){
    glui_alert->close();
    glui_alert=NULL;
  }
  glui_alert = GLUI_Master.create_glui("",0,screenWidth/2,screenHeight/2);
  glui_alert->hide();
  STATIC_alert = glui_alert->add_statictext(_("Loading smoke and fire data"));
}

/* ------------------ hide_glui_trainer ------------------------ */

extern "C" void hide_glui_trainer(void){
  if(glui_trainer!=NULL){
    glui_trainer->hide();
    showtrainer_dialog=0;
    updatemenu=1;
  }
}

/* ------------------ show_glui_trainer ------------------------ */

extern "C" void show_glui_trainer(void){
  if(glui_trainer!=NULL){
    glui_trainer->show();
    showtrainer_dialog=1;
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

  az = camera_current->az_elev;
  elev = camera_current->az_elev+1;

  eye_xyz = camera_current->eye;

  if(TRANSLATE_leftright_inout!=NULL){
    TRANSLATE_leftright_inout->set_x(eye_xyz[0]);
    TRANSLATE_leftright_inout->set_y(eye_xyz[1]);
    TRANSLATE_leftright_inout->set_speed(1.0/(float)screenWidth);
  }

  if(TRANSLATE_updown!=NULL){
    TRANSLATE_updown->set_x(eye_xyz[2]);
    TRANSLATE_updown->set_speed(1.0/(float)screenHeight);
  }

  if(TRANSLATE_az_elev!=NULL){
    TRANSLATE_az_elev->set_x(*az);
    TRANSLATE_az_elev->set_speed(180.0/(float)screenHeight);
    TRANSLATE_az_elev->set_y(*elev);
  }

}

/* ------------------ glui_trainer_setup ------------------------ */

extern "C" void glui_trainer_setup(int main_window){

  update_glui_trainer=0;
  if(glui_trainer!=NULL){
    glui_trainer->close();
    glui_trainer=NULL;
  }
  if(glui_trainer!=NULL)glui_trainer->close();
  glui_trainer = GLUI_Master.create_glui(_("Demonstrator"),0,screenWidth+12,0);
  if(showgluitrainer==0)glui_trainer->hide();
  
  glui_trainer->set_main_gfx_window( main_window );
  PANEL_smokeview = glui_trainer->add_panel(_("Data"));
  BUTTON_smoke3d = glui_trainer->add_button_to_panel(PANEL_smokeview,_("Smoke/Fire"),LOAD_SMOKE,TRAINER_CB);
  if(AnySmoke(NULL)==0)BUTTON_smoke3d->disable();
  BUTTON_temp = glui_trainer->add_button_to_panel(PANEL_smokeview,_("Temperature"),LOAD_TEMP,TRAINER_CB);
  if(AnySlices("TEMPERATURE")==0)BUTTON_temp->disable();
  BUTTON_oxy = glui_trainer->add_button_to_panel(PANEL_smokeview,_("Oxygen"),LOAD_OXY,TRAINER_CB);
  if(AnySlices("oxygen")==0&&AnySlices(_("oxygen VOLUME FRACTION"))==0){
    BUTTON_oxy->disable();
  }
  
  PANEL_explore = glui_trainer->add_panel("Explore",true);

  trainer_path=-1;
  LIST_trainerpath = glui_trainer->add_listbox_to_panel(PANEL_explore,_("Path:"),&trainer_path,TRAINERPATH,TRAINER_CB);
  {
    int i;
    LIST_trainerpath->add_item(-1,_("Manual"));
    LIST_trainerpath->add_item(-2,"-");
    for(i=0;i<ntours;i++){
      tourdata *touri;

      touri = tourinfo + i;
      LIST_trainerpath->add_item(i,touri->menulabel);
    }
  }

  LIST_viewpoint = glui_trainer->add_listbox_to_panel(PANEL_explore,_("Viewpoint:"),&trainer_viewpoints,TRAINERVIEWPOINTS,TRAINER_CB);
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
  BUTTON_toggle_view = glui_trainer->add_button_to_panel(PANEL_explore,_("Toggle View"),TOGGLE_VIEW,TRAINER_CB);
  if(ntrainer_viewpoints<=2)BUTTON_toggle_view->disable();

  CHECKBOX_outline = glui_trainer->add_checkbox_to_panel(PANEL_explore,_("Show walls"),&trainer_outline,TRAINEROUTLINE,TRAINER_CB);
  CHECKBOX_pause = glui_trainer->add_checkbox_to_panel(PANEL_explore,_("Pause"),&trainer_pause,TRAINER_PAUSE,TRAINER_CB);

  update_trainer_outline();
  PANEL_move = glui_trainer->add_panel_to_panel(PANEL_explore,_("Move"),false);
  TRANSLATE_leftright_inout = glui_trainer->add_translation_to_panel(PANEL_move,_("Horizontal"),
    GLUI_TRANSLATION_XY,trainer_xzy,TRAINER_LEFTRIGHT_INOUT,ROTATE_CB);
  glui_trainer->add_column_to_panel(PANEL_move,false);

  TRANSLATE_updown = glui_trainer->add_translation_to_panel(PANEL_move,_("Vertical"),
    GLUI_TRANSLATION_Y,trainer_xzy+2,TRAINER_UPDOWN,ROTATE_CB);
  glui_trainer->add_column_to_panel(PANEL_move,false);

  TRANSLATE_az_elev = glui_trainer->add_translation_to_panel(PANEL_move,_("Rotate"),
    GLUI_TRANSLATION_XY,trainer_ab,TRAINER_AZ_ELEV,ROTATE_CB);

  update_trainer_moves();

  TRAINER_CB(MOVETYPE);
  TRAINER_CB(TRAINERVIEWPOINTS);

}

/* ------------------ ROTATE_CB ------------------------ */

void ROTATE_CB(int var){

  float *eye_xyz, *az, *elev;

  eye_xyz = camera_current->eye;
  az = camera_current->az_elev;
  elev = camera_current->az_elev+1;


  if(rotation_type!=ROTATION_2AXIS){
    rotation_type=ROTATION_2AXIS;
    handle_rotation_type(ROTATION_2AXIS);
    ResetView(RESTORE_EXTERIOR_VIEW);
  }

  if(trainer_viewpoints!=-1){
    LIST_viewpoint->set_int_val(-1);
  }
  if(trainer_path!=-1){
    LIST_trainerpath->set_int_val(-1);
    TRAINER_CB(TRAINERPATH);
  }
  switch(var){
  case TRAINER_AZ_ELEV:
    *az = TRANSLATE_az_elev->get_x();
    *elev = -TRANSLATE_az_elev->get_y();
    break;
  case TRAINER_LEFTRIGHT_INOUT:
    eye_xyz[0]=TRANSLATE_leftright_inout->get_x();
    eye_xyz[1]=TRANSLATE_leftright_inout->get_y();
    break;
  case TRAINER_UPDOWN:
    eye_xyz[2]=TRANSLATE_updown->get_y();
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  camera_current->dirty=1;
}

/* ------------------ TRAINER_CB ------------------------ */

void TRAINER_CB(int var){

  switch(var){
    int i;

  case TRAINER_PAUSE:
    stept=trainer_pause;
    keyboard('t',FROM_SMOKEVIEW);
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
    switch(trainer_path){
    case -1:
      if(trainer_path_old!=-1){
        trainer_pause=0;
        CHECKBOX_pause->set_int_val(trainer_pause);
        TourMenu(MENU_TOUR_MANUAL);
        rotation_type=ROTATION_2AXIS;
        handle_rotation_type(ROTATION_2AXIS);
        from_glui_trainer=1;
        trainee_location=0;
      }
      TRAINER_CB(MOVETYPE);
      break;
    case -2:
      break;
    default:
      if(rotation_type!=EYE_CENTERED){
        rotation_type=EYE_CENTERED;
        handle_rotation_type(ROTATION_2AXIS);
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
    rotation_type=ROTATION_2AXIS;
    handle_rotation_type(ROTATION_2AXIS);
    ResetView(RESTORE_EXTERIOR_VIEW);
    break;
  case LOAD_SMOKE:
    TrainerViewMenu(MENU_TRAINER_smoke);
    break;
  case LOAD_TEMP:
// kind of a hack, having to put in code seg twice, but this is required to get data chopping to work
    if(slicebounds!=NULL&&islicetype!=-1){
      if(setslicechopmin==1||setslicechopmax==1){
       setslicebounds(islicetype);
      }
    }
    TrainerViewMenu(MENU_TRAINER_temp);
    updatechopcolors();
    if(slicebounds!=NULL&&islicetype!=-1){
      if(setslicechopmin==1||setslicechopmax==1){
       setslicebounds(islicetype);
      }
    }
    updatechopcolors();
    colorbarflip=1;
    ColorBarMenu(COLORBAR_FLIP);
    break;
  case LOAD_OXY:
// kind of a hack, having to put in code seg twice, but this is required to get data chopping to work
    if(slicebounds!=NULL&&islicetype!=-1){
      if(setslicechopmin==1||setslicechopmax==1){
       setslicebounds(islicetype);
      }
    }
    TrainerViewMenu(MENU_TRAINER_oxy);
    updatechopcolors();
    if(slicebounds!=NULL&&islicetype!=-1){
      if(setslicechopmin==1||setslicechopmax==1){
       setslicebounds(islicetype);
      }
    }
    updatechopcolors();
    colorbarflip=0;
    ColorBarMenu(COLORBAR_FLIP);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
}

