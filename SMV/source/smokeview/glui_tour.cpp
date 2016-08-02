#define CPP
#include "options.h"

#include <stdio.h>
#include <string.h>
#include GLUT_H
#include <math.h>

#include "update.h"
#include "smokeviewvars.h"

static int viewtype1=REL_VIEW;
static int viewtype2=REL_VIEW;
static float tour_ttt, tour_az_path=0.0, tour_tension=0.0;
static float tour_view_xyz[3]={0.0,0.0,0.0}, tour_elev_path=0.0;
static int tour_hide=0;
static float tour_zoom=1.0;
static char tour_label[sizeof(GLUI_String)];

int nexttour(void);
int prevtour(void);
void TOUR_CB(int var);

GLUI *glui_tour=NULL;

GLUI_Rollout *ROLLOUT_avatar=NULL;
GLUI_Rollout *ROLLOUT_tour = NULL;
GLUI_Rollout *ROLLOUT_keyframe = NULL;
GLUI_Rollout *ROLLOUT_settings = NULL;

GLUI_Panel *PANEL_settingskeyframe=NULL;
GLUI_Panel *PANEL_tension=NULL;
GLUI_Panel *PANEL_path=NULL;
GLUI_Panel *PANEL_tour1=NULL;
GLUI_Panel *PANEL_tour2=NULL;
GLUI_Panel *PANEL_tour3=NULL;
GLUI_Panel *PANEL_tour4=NULL;
GLUI_Panel *PANEL_close_tour=NULL;
GLUI_Panel *PANEL_view=NULL;
GLUI_Panel *PANEL_view2=NULL;
GLUI_Panel *PANEL_view3=NULL;
GLUI_Panel *PANEL_pos=NULL;
GLUI_Panel *PANEL_pos2=NULL;
GLUI_Panel *PANEL_pos3=NULL;
GLUI_Panel *PANEL_view_xyz=NULL;
GLUI_Panel *PANEL_view_angle=NULL;

GLUI_Checkbox *CHECKBOX_snap=NULL,*CHECKBOX_view=NULL,*CHECKBOX_showtourroute=NULL,*CHECKBOX_constantvel=NULL;
GLUI_Checkbox *CHECKBOX_showtour_locus=NULL,*CHECKBOX_showintermediate=NULL;
GLUI_Checkbox *CHECKBOX_globaltension_flag=NULL, *CHECKBOX_tourhide=NULL;
GLUI_Checkbox *CHECKBOX_view1=NULL;
GLUI_Checkbox *CHECKBOX_view2=NULL;

GLUI_Spinner *SPINNER_t=NULL,*SPINNER_x=NULL, *SPINNER_y=NULL,*SPINNER_z=NULL;
GLUI_Spinner *SPINNER_viewx=NULL, *SPINNER_viewy=NULL,*SPINNER_viewz=NULL;
GLUI_Spinner *SPINNER_az_path=NULL;
GLUI_Spinner *SPINNER_globaltourtension=NULL;
GLUI_Spinner *SPINNER_tourtension=NULL;
GLUI_Spinner *SPINNER_tourzoom=NULL;
GLUI_Spinner *SPINNER_elev_path=NULL;

GLUI_Button *BUTTON_next_tour=NULL,*BUTTON_prev_tour=NULL;

GLUI_EditText *EDIT_label=NULL;

GLUI_Listbox *LISTBOX_tour=NULL;
GLUI_Listbox *LISTBOX_avatar=NULL;

#define TOUR_CLOSE 99
#define KEYFRAME_tXYZ 1
#define KEYFRAME_INSERT 2
#define KEYFRAME_DELETE 3
#define KEYFRAME_PREVIOUS 4
#define KEYFRAME_NEXT 5
#define SAVE_SETTINGS 7
#define SHOWTOURROUTE 8
#define VIEWTOURFROMPATH 9
#define CONSTANTTOURVEL 11
#define GLOBAL_TENSION 15
#define GLOBAL_TENSIONFLAG 12
#define TOUR_INSERT 16
#define TOUR_PREVIOUS 17
#define TOUR_NEXT 18
#define TOUR_LABEL 19
#define TOUR_HIDE 20
#define KEYFRAME_viewXYZ 22
#define VIEWSNAP 23
#define TOUR_LIST 24
#define TOUR_AVATAR 31
#define VIEW1 26
#define VIEW2 30
#define VIEW_times 27
#define TOUR_UPDATELABEL 28
#define TOUR_USECURRENT 29

#define TOURMENU(f) callfrom_tourglui=1;TourMenu(f);callfrom_tourglui=0;

#define TOURS_TOURS_ROLLOUT 0
#define SETTINGS_TOURS_ROLLOUT 1
#define KEYFRAME_TOURS_ROLLOUT 2

procdata toursprocinfo[3];
int ntoursprocinfo = 0;

/* ------------------ is_tour_open ------------------------ */

int is_tour_open(void){
  return 0;
}

/* ------------------ Tours_Rollout_CB ------------------------ */

void Tours_Rollout_CB(int var){
  toggle_rollout(toursprocinfo, ntoursprocinfo, var);
}


/* ------------------ update_tour_state ------------------------ */

extern "C" void update_tour_state(void){
  TOUR_CB(SHOWTOURROUTE);
  TOUR_CB(VIEWTOURFROMPATH);
  TOUR_CB(VIEWSNAP);
  TOUR_CB(GLOBAL_TENSION);
  TOUR_CB(GLOBAL_TENSIONFLAG);
}

/* ------------------ add_delete_keyframe ------------------------ */

void add_delete_keyframe(int flag){
  if(flag==ADD_KEYFRAME)TOUR_CB(KEYFRAME_INSERT);
  if(flag==DELETE_KEYFRAME)TOUR_CB(KEYFRAME_DELETE);
}


/* ------------------ update_edit_tour ------------------------ */

extern "C" void update_edit_tour(void){
  TOUR_CB(SHOWTOURROUTE);
}

/* ------------------ update_tour_parms ------------------------ */

extern "C" void update_tour_parms(void){
  TOUR_CB(KEYFRAME_tXYZ);
}

/* ------------------ add_new_tour ------------------------ */

extern "C" void add_new_tour(void){
  TOUR_CB(TOUR_INSERT);
}

/* ------------------ glui_tour_setup ------------------------ */

extern "C" void glui_tour_setup(int main_window){

  int i;

  update_glui_tour=0;
  if(glui_tour!=NULL){
    glui_tour->close();
    glui_tour=NULL;
  }
  glui_tour = GLUI_Master.create_glui(_d("Tours"),0,0,0);
  glui_tour->hide();


  ROLLOUT_tour = glui_tour->add_rollout("Tours",true,TOURS_TOURS_ROLLOUT, Tours_Rollout_CB);
  ADDPROCINFO(toursprocinfo, ntoursprocinfo, ROLLOUT_tour, TOURS_TOURS_ROLLOUT);


  PANEL_tour1 = glui_tour->add_panel_to_panel(ROLLOUT_tour,"",GLUI_PANEL_NONE);

  BUTTON_next_tour=glui_tour->add_button_to_panel(PANEL_tour1,_d("Next"),TOUR_NEXT,TOUR_CB);
  glui_tour->add_column_to_panel(PANEL_tour1);
  BUTTON_prev_tour=glui_tour->add_button_to_panel(PANEL_tour1,_d("Previous"),TOUR_PREVIOUS,TOUR_CB);

  PANEL_tour4 = glui_tour->add_panel_to_panel(ROLLOUT_tour,"",GLUI_PANEL_NONE);

  if(ntours>0){
    selectedtour_index = TOURINDEX_MANUAL;
    selectedtour_index_old = TOURINDEX_MANUAL;
    LISTBOX_tour=glui_tour->add_listbox_to_panel(PANEL_tour4,"",&selectedtour_index,TOUR_LIST,TOUR_CB);

    LISTBOX_tour->add_item(TOURINDEX_MANUAL, "Manual");
    LISTBOX_tour->add_item(-999,"-");
    for(i=0;i<ntours;i++){
      tourdata *touri;

      touri = tourinfo + i;
      LISTBOX_tour->add_item(i,touri->label);
    }
    LISTBOX_tour->set_int_val(selectedtour_index);
    glui_tour->add_column_to_panel(PANEL_tour1,false);
  }
  glui_tour->add_column_to_panel(PANEL_tour4,false);
  glui_tour->add_button_to_panel(PANEL_tour4,_d("New"),TOUR_INSERT,TOUR_CB);

  ROLLOUT_avatar = glui_tour->add_rollout_to_panel(ROLLOUT_tour,"Label/Avatars",false);
  PANEL_tour3 = glui_tour->add_panel_to_panel(ROLLOUT_avatar,"",GLUI_PANEL_NONE);
  EDIT_label=glui_tour->add_edittext_to_panel(PANEL_tour3,"Label:",GLUI_EDITTEXT_TEXT,tour_label,TOUR_LABEL,TOUR_CB);
  glui_tour->add_column_to_panel(PANEL_tour3,false);
  glui_tour->add_button_to_panel(PANEL_tour3,_d("Update"),TOUR_UPDATELABEL,TOUR_CB);
  EDIT_label->set_w(240);

  PANEL_tour2 = glui_tour->add_panel_to_panel(ROLLOUT_avatar,"",GLUI_PANEL_NONE);
  if(navatar_types>0){
    LISTBOX_avatar=glui_tour->add_listbox_to_panel(PANEL_tour2,_d("Avatar:"),&glui_avatar_index,TOUR_AVATAR,TOUR_CB);
    for(i=0;i<navatar_types;i++){
      LISTBOX_avatar->add_item(i,avatar_types[i]->label);
    }
    if(tourlocus_type==0){
      glui_avatar_index=-1;
    }
    else if(tourlocus_type==1){
      glui_avatar_index=-2;
    }
    else{
      glui_avatar_index=iavatar_types;
    }
    LISTBOX_avatar->set_int_val(glui_avatar_index);
    glui_tour->add_column_to_panel(PANEL_tour2,false);
    CHECKBOX_showtour_locus=glui_tour->add_checkbox_to_panel(PANEL_tour2,_d("Show avatar"),&show_tourlocus);
  }

  ROLLOUT_settings = glui_tour->add_rollout(_d("Settings"),true,SETTINGS_TOURS_ROLLOUT, Tours_Rollout_CB);
  ADDPROCINFO(toursprocinfo, ntoursprocinfo, ROLLOUT_settings, SETTINGS_TOURS_ROLLOUT);

  CHECKBOX_showtourroute=glui_tour->add_checkbox_to_panel(ROLLOUT_settings,_d("Edit tour"),&edittour,SHOWTOURROUTE,TOUR_CB);
  CHECKBOX_view=glui_tour->add_checkbox_to_panel(ROLLOUT_settings,_d("View from tour path"),&viewtourfrompath,VIEWTOURFROMPATH,TOUR_CB);
  CHECKBOX_snap=glui_tour->add_checkbox_to_panel(ROLLOUT_settings,_d("View from selected keyframe"),&keyframe_snap,VIEWSNAP,TOUR_CB);
  CHECKBOX_constantvel=glui_tour->add_checkbox_to_panel(ROLLOUT_settings,_d("Constant speed"),&tour_constant_vel,CONSTANTTOURVEL,TOUR_CB);
  CHECKBOX_showintermediate=glui_tour->add_checkbox_to_panel(ROLLOUT_settings,_d("Show intermediate path nodes"),&show_path_knots);
#ifdef _DEBUG
  glui_tour->add_checkbox_to_panel(ROLLOUT_settings,_d("Antialias tour path line"),&tour_antialias);
#endif

  PANEL_path = glui_tour->add_panel_to_panel(ROLLOUT_settings,_d("Duration"),true);

  glui_tour->add_spinner_to_panel(PANEL_path,_d("start time"),GLUI_SPINNER_FLOAT,&view_tstart,VIEW_times,TOUR_CB);
  glui_tour->add_spinner_to_panel(PANEL_path,_d("stop time:"),GLUI_SPINNER_FLOAT,&view_tstop, VIEW_times,TOUR_CB);
  glui_tour->add_spinner_to_panel(PANEL_path,_d("points"),    GLUI_SPINNER_INT,&view_ntimes,  VIEW_times,TOUR_CB);

  PANEL_tension = glui_tour->add_panel_to_panel(ROLLOUT_settings,_d("Tension"),true);
  CHECKBOX_globaltension_flag=glui_tour->add_checkbox_to_panel(PANEL_tension,_d("Global"),&tour_global_tension_flag,GLOBAL_TENSIONFLAG,TOUR_CB);
  SPINNER_globaltourtension=glui_tour->add_spinner_to_panel(PANEL_tension,_d("All keyframes"),GLUI_SPINNER_FLOAT,&tour_global_tension,GLOBAL_TENSION,TOUR_CB);
  SPINNER_tourtension=glui_tour->add_spinner_to_panel(PANEL_tension,_d("Selected keyframe"),GLUI_SPINNER_FLOAT,&tour_tension,KEYFRAME_tXYZ,TOUR_CB);

  SPINNER_globaltourtension->set_float_limits(-1.0,1.0,GLUI_LIMIT_CLAMP);
  SPINNER_tourtension->set_float_limits(-1.0,1.0,GLUI_LIMIT_CLAMP);

  ROLLOUT_keyframe = glui_tour->add_rollout("Keyframe",true,KEYFRAME_TOURS_ROLLOUT, Tours_Rollout_CB);
  ADDPROCINFO(toursprocinfo, ntoursprocinfo, ROLLOUT_keyframe, KEYFRAME_TOURS_ROLLOUT);


  PANEL_pos = glui_tour->add_panel_to_panel(ROLLOUT_keyframe,"",GLUI_PANEL_NONE);

  PANEL_pos3 = glui_tour->add_panel_to_panel(PANEL_pos,"",GLUI_PANEL_NONE);
  glui_tour->add_button_to_panel(PANEL_pos3,_d("Next"),    KEYFRAME_NEXT,TOUR_CB);
  glui_tour->add_button_to_panel(PANEL_pos3,_d("Previous"),KEYFRAME_PREVIOUS,TOUR_CB);
  glui_tour->add_button_to_panel(PANEL_pos3,_d("Add"),KEYFRAME_INSERT,TOUR_CB);
  glui_tour->add_button_to_panel(PANEL_pos3,_d("Delete"),KEYFRAME_DELETE,TOUR_CB);

  glui_tour->add_column_to_panel(PANEL_pos,false);

  PANEL_pos2 = glui_tour->add_panel_to_panel(PANEL_pos,"Position/time");
  SPINNER_t=glui_tour->add_spinner_to_panel(PANEL_pos2,"t:",GLUI_SPINNER_FLOAT,&tour_ttt,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_x=glui_tour->add_spinner_to_panel(PANEL_pos2,"X:",GLUI_SPINNER_FLOAT,tour_xyz,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_y=glui_tour->add_spinner_to_panel(PANEL_pos2,"Y:",GLUI_SPINNER_FLOAT,tour_xyz+1,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_z=glui_tour->add_spinner_to_panel(PANEL_pos2,"Z:",GLUI_SPINNER_FLOAT,tour_xyz+2,KEYFRAME_tXYZ,TOUR_CB);

  PANEL_view = glui_tour->add_panel_to_panel(ROLLOUT_keyframe,"View direction",true);
  PANEL_view3 = glui_tour->add_panel_to_panel(PANEL_view,"",GLUI_PANEL_NONE);
  CHECKBOX_view1=glui_tour->add_checkbox_to_panel(PANEL_view3,"Absolute",&viewtype1,VIEW1,TOUR_CB);
  glui_tour->add_column_to_panel(PANEL_view3,false);
  CHECKBOX_view2=glui_tour->add_checkbox_to_panel(PANEL_view3,"Relative to path",&viewtype2,VIEW2,TOUR_CB);

  PANEL_view2 = glui_tour->add_panel_to_panel(PANEL_view,"",GLUI_PANEL_NONE);
  PANEL_view_xyz = glui_tour->add_panel_to_panel(PANEL_view2,"",GLUI_PANEL_NONE);
  glui_tour->add_column_to_panel(PANEL_view2,false);
  PANEL_view_angle = glui_tour->add_panel_to_panel(PANEL_view2,"",GLUI_PANEL_NONE);


  SPINNER_viewx=glui_tour->add_spinner_to_panel(PANEL_view_xyz,"X",GLUI_SPINNER_FLOAT,tour_view_xyz,KEYFRAME_viewXYZ,TOUR_CB);
  SPINNER_viewy=glui_tour->add_spinner_to_panel(PANEL_view_xyz,"Y",GLUI_SPINNER_FLOAT,tour_view_xyz+1,KEYFRAME_viewXYZ,TOUR_CB);
  SPINNER_viewz=glui_tour->add_spinner_to_panel(PANEL_view_xyz,"Z",GLUI_SPINNER_FLOAT,tour_view_xyz+2,KEYFRAME_viewXYZ,TOUR_CB);
  SPINNER_az_path=glui_tour->add_spinner_to_panel(PANEL_view_angle,_d("Azimuth:"),GLUI_SPINNER_FLOAT,&tour_az_path,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_elev_path=glui_tour->add_spinner_to_panel(PANEL_view_angle,_d("Elevation:"),GLUI_SPINNER_FLOAT,&tour_elev_path,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_elev_path->set_float_limits(-90.0,90.0);
  SPINNER_tourzoom=glui_tour->add_spinner_to_panel(PANEL_view_angle,_d("Zoom:"),GLUI_SPINNER_FLOAT,&tour_zoom,KEYFRAME_tXYZ,TOUR_CB);

  PANEL_close_tour = glui_tour->add_panel("",false);
  glui_tour->add_button_to_panel(PANEL_close_tour,_d("Save settings"),SAVE_SETTINGS,TOUR_CB);
  glui_tour->add_column_to_panel(PANEL_close_tour,false);
  glui_tour->add_button_to_panel(PANEL_close_tour,"Close",TOUR_CLOSE,TOUR_CB);

  ROLLOUT_keyframe->close();
  ROLLOUT_settings->close();

  glui_tour->set_main_gfx_window( main_window );

  TOUR_CB(VIEW1);
  TOUR_CB(GLOBAL_TENSIONFLAG);

  SPINNER_az_path->set_float_limits(-180.0,180.0);
  update_tourcontrols();
  update_tourlist=1;
}

/* ------------------ Update_Tourlist(void) ------------------------ */

extern "C" void Update_Tourlist(void){

  update_tourlist=0;
  TOUR_CB(TOUR_LIST);
}

/* ------------------ hide_glui_tour ------------------------ */

extern "C" void hide_glui_tour(void){
  if(glui_tour!=NULL)glui_tour->hide();
  showtour_dialog=0;
  updatemenu=1;
}

/* ------------------ show_glui_tour ------------------------ */

extern "C" void show_glui_tour(void){
  showtour_dialog=1;
  if(glui_tour!=NULL)glui_tour->show();
  updatemenu=1;
}

/* ------------------ trim_val ------------------------ */

extern "C" float trim_val(float val){
  if(val<0.000001&&val>-0.000001){
    return 0.0;
  }
  else{
    return val;
  }
}

/* ------------------ update_glui_keyframe ------------------------ */

extern "C" void update_glui_keyframe(void){
  SPINNER_x->set_float_val(tour_xyz[0]);
  SPINNER_y->set_float_val(tour_xyz[1]);
  SPINNER_z->set_float_val(tour_xyz[2]);
}

/* ------------------ set_glui_keyframe ------------------------ */

extern "C" void set_glui_keyframe(void){
  tourdata *ti;
  float *eye,*xyz_view;

  if(selected_frame==NULL)return;

  ti = selected_tour;
  if(ti==NULL)return;

  tour_hide=1-ti->display;
  if(selected_tour!=NULL)strcpy(tour_label,selected_tour->label);
  glui_avatar_index=ti->glui_avatar_index;
  TOUR_CB(TOUR_AVATAR);
  LISTBOX_avatar->set_int_val(glui_avatar_index);
  eye = selected_frame->nodeval.eye;
  xyz_view = selected_frame->nodeval.xyz_view_abs;

  tour_ttt = selected_frame->disp_time;
  tour_xyz[0] = trim_val(DENORMALIZE_X(eye[0]));
  tour_xyz[1] = trim_val(DENORMALIZE_Y(eye[1]));
  tour_xyz[2] = trim_val(DENORMALIZE_Z(eye[2]));
  tour_view_xyz[0] = trim_val(DENORMALIZE_X(xyz_view[0]));
  tour_view_xyz[1] = trim_val(DENORMALIZE_Y(xyz_view[1]));
  tour_view_xyz[2] = trim_val(DENORMALIZE_Z(xyz_view[2]));
  tour_az_path = selected_frame->az_path;
  tour_continuity=selected_frame->continuity;
  tour_bias=selected_frame->bias;
  viewtype1=selected_frame->viewtype;
  viewtype2=1-viewtype1;
  tour_tension=selected_frame->tension;
  tour_zoom=selected_frame->nodeval.zoom;
  tour_elev_path=selected_frame->nodeval.elev_path;

  tour_global_tension_flag=selected_tour->global_tension_flag;
  tour_global_tension=selected_tour->global_tension;


  if(SPINNER_t==NULL)return;

  CHECKBOX_globaltension_flag->set_int_val(tour_global_tension_flag);
  SPINNER_globaltourtension->set_float_val(tour_global_tension);
  if(tour_global_tension_flag==1){
    SPINNER_globaltourtension->enable();
    SPINNER_tourtension->disable();
  }
  else{
    SPINNER_globaltourtension->disable();
    SPINNER_tourtension->enable();
  }

  {
    float time_temp;

    time_temp=tour_ttt;
    SPINNER_t->set_float_limits(selected_frame->prev->disp_time,selected_frame->next->disp_time);
    tour_ttt=time_temp;
    SPINNER_t->set_float_val(tour_ttt);
  }

  if(edittour==1){
    if(tour_constant_vel==0){
      SPINNER_t->enable();
    }
    else{
      SPINNER_t->disable();
    }
  }

  SPINNER_x->set_float_val(tour_xyz[0]);
  SPINNER_y->set_float_val(tour_xyz[1]);
  SPINNER_z->set_float_val(tour_xyz[2]);
  SPINNER_tourtension->set_float_val(tour_tension);
  SPINNER_tourzoom->set_float_val(tour_zoom);
  SPINNER_viewx->set_float_val(tour_view_xyz[0]);
  SPINNER_viewy->set_float_val(tour_view_xyz[1]);
  SPINNER_viewz->set_float_val(tour_view_xyz[2]);
  SPINNER_az_path->set_float_val(tour_az_path);
  SPINNER_elev_path->set_float_val(tour_elev_path);
  if(CHECKBOX_tourhide!=NULL)CHECKBOX_tourhide->set_int_val(tour_hide);
  EDIT_label->set_text(tour_label);

  if(edittour==1){
    if(viewtype1==ABS_VIEW){
      SPINNER_az_path->disable();
      SPINNER_elev_path->disable();
      SPINNER_viewx->enable();
      SPINNER_viewy->enable();
      SPINNER_viewz->enable();
    }
    else{
      SPINNER_az_path->enable();
      SPINNER_elev_path->enable();
      SPINNER_viewx->disable();
      SPINNER_viewy->disable();
      SPINNER_viewz->disable();
    }
  }
  else{
    SPINNER_az_path->disable();
    SPINNER_elev_path->disable();
    SPINNER_viewx->disable();
    SPINNER_viewy->disable();
    SPINNER_viewz->disable();
  }
  CHECKBOX_view1->set_int_val(viewtype1);
  CHECKBOX_view2->set_int_val(viewtype2);
}

extern "C" void update_tourindex(void){
  update_selectedtour_index=0;
  selectedtour_index=selectedtour_index_ini;
  TOUR_CB(TOUR_LIST);
}

/* ------------------ TOUR_CB ------------------------ */

void TOUR_CB(int var){
  keyframe *thiskey,*nextkey,*newframe;
  keyframe *lastkey;
  tourdata *thistour=NULL;
  float *xyz_view,*eye;

  float key_xyz[3];
  float key_params[3];
  float key_time_in, key_az_path, key_view[3], key_zoom;
  float key_elev_path, key_bank;

  if(ntours==0&&var!=TOUR_INSERT&&var!=TOUR_CLOSE&&var!=SAVE_SETTINGS){
    return;
  }
  //if(selected_frame==NULL&&tourinfo!=NULL){
  //  selected_frame=tourinfo[0].first_frame.next;
  //  selected_tour=tourinfo;
  //  if(selected_frame->next==NULL)selected_frame=NULL;
  //  set_glui_keyframe();
  //}
  if(selected_frame!=NULL){
    thistour=selected_tour;
  }

  switch(var){
  case TOUR_USECURRENT:
    break;
  case TOUR_NEXT:
    if(nexttour()==1){
      selected_tour->display=0;
      TOURMENU(selectedtour_index);
      set_glui_keyframe();
    }
    break;
  case TOUR_PREVIOUS:
    if(prevtour()==1){
      selected_tour->display=0;
      TOURMENU(selectedtour_index);
      set_glui_keyframe();
    }
    break;
  case TOUR_CLOSE:
    hide_glui_tour();
    break;
  case SAVE_SETTINGS:
    WriteINI(LOCAL_INI,NULL);
    break;
  case SHOWTOURROUTE:
    edittour = 1 - edittour;
    TOURMENU(MENU_TOUR_SHOWDIALOG);
    update_tourcontrols();
    TOUR_CB(VIEW1);
    updatemenu=0;
    break;
  case VIEWSNAP:
    if(viewtourfrompath==1&&keyframe_snap==1){
      viewtourfrompath=0;
      CHECKBOX_view->set_int_val(viewtourfrompath);
    }
    TOUR_CB(VIEWTOURFROMPATH);
    updatemenu=0;
    break;
  case VIEWTOURFROMPATH:
    if(viewtourfrompath==1&&keyframe_snap==1){
      keyframe_snap=0;
      CHECKBOX_snap->set_int_val(keyframe_snap);
    }
    viewtourfrompath = 1 - viewtourfrompath;
    TOURMENU(MENU_TOUR_VIEWFROMROUTE);
    break;
  case VIEW2:
    viewtype1=1-viewtype2;
    viewtype2=1-viewtype1;
    TOUR_CB(VIEW1);
    CHECKBOX_view1->set_int_val(viewtype1);
    break;
  case VIEW1:
    viewtype2 = 1 - viewtype1;
    CHECKBOX_view2->set_int_val(viewtype2);
    if(viewtype1==1&&edittour==1){
      SPINNER_az_path->disable();
      SPINNER_elev_path->disable();
      SPINNER_viewx->enable();
      SPINNER_viewy->enable();
      SPINNER_viewz->enable();
      if(selected_frame!=NULL){
        xyzview2azelev(selected_frame,NULL,NULL);
        SPINNER_az_path->set_float_val(tour_az_path);
        SPINNER_elev_path->set_float_val(tour_elev_path);
      }
    }
    else if(viewtype1==REL_VIEW&&edittour==1){
      SPINNER_az_path->enable();
      SPINNER_elev_path->enable();
      SPINNER_viewx->disable();
      SPINNER_viewy->disable();
      SPINNER_viewz->disable();
    }
    else if(edittour==0){
      SPINNER_az_path->disable();
      SPINNER_elev_path->disable();
      SPINNER_viewx->disable();
      SPINNER_viewy->disable();
      SPINNER_viewz->disable();
    }
    if(selected_frame!=NULL){
      selected_frame->viewtype=viewtype1;
    }
    createtourpaths();
    break;
  case VIEW_times:
    ReallocTourMemory();
    createtourpaths();
    UpdateTimes();
    break;
  case KEYFRAME_viewXYZ:
    if(selected_frame!=NULL){
      if(selected_tour-tourinfo==0)dirtycircletour=1;
      selected_tour->startup=0;
      xyz_view = selected_frame->nodeval.xyz_view_abs;
      NORMALIZE_XYZ(xyz_view,tour_view_xyz);

      xyzview2azelev(selected_frame,&tour_az_path,&tour_elev_path);
      SPINNER_az_path->set_float_val(tour_az_path);
      SPINNER_elev_path->set_float_val(tour_elev_path);

      createtourpaths();
      selected_frame->selected=1;
    }
    break;
  case KEYFRAME_tXYZ:
    if(selected_frame!=NULL){
      if(selected_tour-tourinfo==0)dirtycircletour=1;
      selected_tour->startup=0;
      eye = selected_frame->nodeval.eye;
      xyz_view = selected_frame->nodeval.xyz_view_abs;

      if(tour_constant_vel==0){
        selected_frame->noncon_time=tour_ttt;
        selected_frame->disp_time=tour_ttt;
      }
      NORMALIZE_XYZ(eye,tour_xyz);
      if(viewtype1==REL_VIEW){
        tour_az_path = SPINNER_az_path->get_float_val();
      }
      selected_frame->az_path=tour_az_path;
      selected_frame->nodeval.elev_path=tour_elev_path;

      selected_frame->tension=tour_tension;
      selected_frame->bias=tour_bias;
      selected_frame->continuity=tour_continuity;
      selected_frame->viewtype=viewtype1;
      selected_frame->nodeval.zoom=tour_zoom;
      NORMALIZE_XYZ(xyz_view,tour_view_xyz);
      createtourpaths();
      selected_frame->selected=1;
      if(viewtype1==ABS_VIEW){
        TOUR_CB(KEYFRAME_viewXYZ);
      }
    }
    break;
  case GLOBAL_TENSIONFLAG:
    if(selected_tour!=NULL&&SPINNER_globaltourtension!=NULL){
      selected_tour->global_tension_flag=tour_global_tension_flag;
      if(tour_global_tension_flag==1){
        SPINNER_globaltourtension->enable();
        SPINNER_tourtension->disable();
      }
      else{
        SPINNER_globaltourtension->disable();
        SPINNER_tourtension->enable();
      }
      createtourpaths();
    }
    break;
  case GLOBAL_TENSION:
    if(selected_tour!=NULL){
      selected_tour->global_tension=tour_global_tension;
      createtourpaths();
    }
    break;
  case KEYFRAME_NEXT:
    if(selected_frame==NULL&&tourinfo!=NULL){
      selected_frame=&(tourinfo[0].first_frame);
      selected_tour=tourinfo;
    }
    if(selected_frame!=NULL){
      thistour=selected_tour;
      if(selected_frame->next!=&thistour->last_frame){
        new_select(selected_frame->next);
      }
      else{
        new_select(thistour->first_frame.next);
      }
    }
    set_glui_keyframe();
    break;
  case KEYFRAME_PREVIOUS:
    if(selected_frame==NULL&&tourinfo!=NULL){
      selected_frame=&(tourinfo[0].last_frame);
      selected_tour=tourinfo;
    }
    if(selected_frame!=NULL){
      thistour=selected_tour;
      selected_tour=thistour;
      if(selected_frame->prev!=&thistour->first_frame){
        new_select(selected_frame->prev);
      }
      else{
        new_select(thistour->last_frame.prev);
      }
    }
    set_glui_keyframe();
    break;
  case CONSTANTTOURVEL:
    update_tourcontrols();
    createtourpaths();
    UpdateTimes();
    set_glui_keyframe();
    break;
  case KEYFRAME_INSERT:
    if(selected_frame!=NULL){
      thistour=selected_tour;
      thiskey=selected_frame;
      nextkey=thiskey->next;
      if(nextkey==&thistour->last_frame){
        lastkey=thiskey->prev;
        key_xyz[0]=DENORMALIZE_X(2*thiskey->nodeval.eye[0]-lastkey->nodeval.eye[0]);
        key_xyz[1]=DENORMALIZE_Y(2*thiskey->nodeval.eye[1]-lastkey->nodeval.eye[1]);
        key_xyz[2]=DENORMALIZE_Z(2*thiskey->nodeval.eye[2]-lastkey->nodeval.eye[2]);
        key_az_path = (2*thiskey->az_path-lastkey->az_path);
        key_elev_path=(2*thiskey->nodeval.elev_path-lastkey->nodeval.elev_path);
        key_time_in = thiskey->noncon_time;
        thiskey->noncon_time=(thiskey->noncon_time+lastkey->noncon_time)/2.0;
        key_params[0]=(2*thiskey->bias-lastkey->bias);
        key_params[1]=(2*thiskey->continuity-lastkey->continuity);
        key_params[2]=(2*thiskey->tension-lastkey->tension);
        key_view[0]=DENORMALIZE_X(2*thiskey->nodeval.xyz_view_abs[0]-lastkey->nodeval.xyz_view_abs[0]);
        key_view[1]=DENORMALIZE_Y(2*thiskey->nodeval.xyz_view_abs[1]-lastkey->nodeval.xyz_view_abs[1]);
        key_view[2]=DENORMALIZE_Z(2*thiskey->nodeval.xyz_view_abs[2]-lastkey->nodeval.xyz_view_abs[2]);
        key_zoom = (2*thiskey->nodeval.zoom - lastkey->nodeval.zoom);
        key_bank = (2*thiskey->bank - lastkey->bank);
        viewtype1=thiskey->viewtype;
        viewtype2=1-viewtype1;
      }
      else{
        key_xyz[0]=DENORMALIZE_X((thiskey->nodeval.eye[0]+nextkey->nodeval.eye[0])/2.0);
        key_xyz[1]=DENORMALIZE_Y((thiskey->nodeval.eye[1]+nextkey->nodeval.eye[1])/2.0);
        key_xyz[2]=DENORMALIZE_Z((thiskey->nodeval.eye[2]+nextkey->nodeval.eye[2])/2.0);
        key_az_path = (thiskey->az_path+nextkey->az_path)/2.0;
        key_elev_path=(thiskey->nodeval.elev_path+nextkey->nodeval.elev_path)/2.0;
        key_time_in = (thiskey->noncon_time+nextkey->noncon_time)/2.0;
        key_params[0]=(thiskey->bias+nextkey->bias)/2.0;
        key_params[1]=(thiskey->continuity+nextkey->continuity)/2.0;
        key_params[2]=(thiskey->tension+nextkey->tension)/2.0;
        key_view[0]=DENORMALIZE_X((thiskey->nodeval.xyz_view_abs[0]+nextkey->nodeval.xyz_view_abs[0])/2.0);
        key_view[1]=DENORMALIZE_Y((thiskey->nodeval.xyz_view_abs[1]+nextkey->nodeval.xyz_view_abs[1])/2.0);
        key_view[2]=DENORMALIZE_Z((thiskey->nodeval.xyz_view_abs[2]+nextkey->nodeval.xyz_view_abs[2])/2.0);
        key_zoom = (thiskey->nodeval.zoom + nextkey->nodeval.zoom)/2.0;
        key_bank = (thiskey->bank + nextkey->bank)/2.0;
        if(thiskey->viewtype==REL_VIEW&&nextkey->viewtype==REL_VIEW){
          viewtype1=REL_VIEW;
        }
        else{
          viewtype1=ABS_VIEW;
          if(thiskey->viewtype==ABS_VIEW){
            DENORMALIZE_XYZ(key_view,thiskey->nodeval.xyz_view_abs);
            key_elev_path = thiskey->nodeval.elev_path;
          }
          if(thiskey->viewtype==REL_VIEW&&nextkey->viewtype==ABS_VIEW){
            DENORMALIZE_XYZ(key_view,nextkey->nodeval.xyz_view_abs);
            key_elev_path = nextkey->nodeval.elev_path;
          }
        }
        viewtype2=1-viewtype1;
      }
      newframe=add_frame(selected_frame,key_time_in,key_xyz,key_az_path,key_elev_path,key_bank,key_params,viewtype1,key_zoom,key_view);
      createtourpaths();
      new_select(newframe);
      set_glui_keyframe();
    }
    break;
  case KEYFRAME_DELETE:
    if(selected_frame!=NULL){
      selected_frame=delete_frame(selected_frame);
      if(selected_frame!=NULL){
        selected_frame->selected=1;
        createtourpaths();
      }
      else{
        if(thistour!=NULL)delete_tour(thistour-tourinfo);
      }
    }
    break;
  case TOUR_AVATAR:
    if(selected_tour->glui_avatar_index!=glui_avatar_index){
      selected_tour->glui_avatar_index=glui_avatar_index;
// hack to make touring avatar show up
//      avatar_types[glui_avatar_index]->visible=1;
      updatemenu=1;
    }
    if(glui_avatar_index==-1){
      tourlocus_type=0;
    }
    else if(glui_avatar_index==-2){
      tourlocus_type=1;
    }
    else{
      tourlocus_type=2;
      iavatar_types=glui_avatar_index;
    }
    break;
  case TOUR_UPDATELABEL:
    // supposed to fall through to TOUR_LIST
  case TOUR_LIST:
    if(selectedtour_index==-999){
      selectedtour_index=selectedtour_index_old;
      if(selectedtour_index==-999)selectedtour_index = TOURINDEX_MANUAL;
      TOUR_CB(TOUR_LIST);
      return;
    }
    switch(selectedtour_index){
    case TOURINDEX_ALL:
      TOURMENU(MENU_TOUR_SHOWALL); // show all tours
      set_glui_keyframe();
      break;
    case TOURINDEX_MANUAL:
      edittour=0;
      TOURMENU(MENU_TOUR_CLEARALL);  // reset tour vis to ini values
      break;
    case TOURINDEX_DEFAULT:
      TOURMENU(MENU_TOUR_DEFAULT);  // default tour
      break;
    default:
      selected_tour=tourinfo + selectedtour_index;
      selected_frame=selected_tour->first_frame.next;
      selected_tour->display=0;
      TOURMENU(selectedtour_index);
      set_glui_keyframe();
      break;
    }
    delete_tourlist();
    create_tourlist();
    updateviewtour();
    update_tourcontrols();
    selectedtour_index_old=selectedtour_index;
    break;
  case TOUR_INSERT:
    if(CHECKBOX_showtourroute!=NULL&&edittour==0)CHECKBOX_showtourroute->set_int_val(1);
    thistour=add_tour(NULL);
    selected_frame=thistour->first_frame.next;
    selected_tour=thistour;
    selectedtour_index = thistour - tourinfo;
    selectedtour_index_old=selectedtour_index;
    set_glui_keyframe();
    createtourpaths();
    updateviewtour();
    update_tourcontrols();
    selected_tour->display=0;
    TOURMENU(selectedtour_index);
    updatemenu=1;
    break;
  case TOUR_LABEL:
    if(thistour!=NULL){
      strcpy(thistour->label,tour_label);
      set_glui_keyframe();
      if(LISTBOX_tour!=NULL){
        LISTBOX_tour->delete_item(thistour-tourinfo);
        LISTBOX_tour->add_item(thistour-tourinfo,thistour->label);
      }
      update_tour_menulabels();
      updatemenu=1;
    }
    break;
  case TOUR_HIDE:
    if(thistour!=NULL){
      if(tour_hide==1){
        thistour->display=1;
        TOURMENU(thistour-tourinfo);
        nexttour();
        set_glui_keyframe();
        thistour->display=0;
      }
      else{
        thistour->display=1;
      }
      updatemenu=1;
      delete_tourlist();
      create_tourlist();
      updateviewtour();
      update_tourcontrols();
    }
    break;
  default:
    ASSERT(FFALSE);
  }
}


/* ------------------ delete_tourlsit ------------------------ */

extern "C" void delete_tourlist(void){
  int i;

  if(LISTBOX_tour==NULL)return;
  for(i=0;i<ntours;i++){
    LISTBOX_tour->delete_item(i);
  }
  delete_vol_tourlist(); //xx comment this line if smokebot fails with seg fault
}

/* ------------------ create_tourlist ------------------------ */

extern "C" void create_tourlist(void){
  int i;

  if(LISTBOX_tour==NULL)return;
  for(i=0;i<ntours;i++){
    tourdata *touri;
    char label[1000];

    touri = tourinfo + i;
    strcpy(label,"");
    if(i==selectedtour_index)strcat(label,"*");
    if(touri->label!=NULL&&strlen(touri->label)>0)strcat(label,touri->label);
    if(strlen(label)>0){
      LISTBOX_tour->add_item(i,label);
    }
    else{
      LISTBOX_tour->add_item(i,"error");
    }
  }
  if(selectedtour_index>=-1&&selectedtour_index<ntours)LISTBOX_tour->set_int_val(selectedtour_index);

 create_vol_tourlist(); //xx comment this line if smokebot fails with seg fault
}

/* ------------------ nexttour ------------------------ */

int nexttour(void){
  int i;

  i = selectedtour_index + 1;
  if(i>ntours-1)i=0;
  if(i>=0&&i<ntours){
    selectedtour_index=i;
    selected_tour = tourinfo + i;
    selected_frame = selected_tour->first_frame.next;
    return 1;
  }
  return 0;
}

/* ------------------ prevtour ------------------------ */

int prevtour(void){
  int i;

  i=selectedtour_index-1;
  if(i<0)i=ntours-1;
  if(i>=0&&i<ntours){
    selectedtour_index=i;
    selected_tour = tourinfo + i;
    selected_frame = selected_tour->first_frame.next;
    return 1;
  }
  return 0;
}

/* ------------------ update_tourcontrols ------------------------ */

extern "C" void update_tourcontrols(void){

  if(BUTTON_next_tour==NULL)return;
  if(BUTTON_prev_tour==NULL)return;
  if(ROLLOUT_keyframe==NULL)return;
  if(SPINNER_t==NULL)return;
  if(CHECKBOX_showtourroute!=NULL)CHECKBOX_showtourroute->set_int_val(edittour);
  if(CHECKBOX_view!=NULL)CHECKBOX_view->set_int_val(viewtourfrompath);
  if(ntours>1){
    BUTTON_next_tour->enable();
    BUTTON_prev_tour->enable();
  }
  else{
    BUTTON_next_tour->disable();
    BUTTON_prev_tour->disable();
  }
  if(ntours>0&&edittour==1){
    if(SPINNER_t!=NULL)SPINNER_t->enable();
    if(ROLLOUT_keyframe!=NULL&&ROLLOUT_keyframe->enabled==0)ROLLOUT_keyframe->enable();
    if(SPINNER_az_path!=NULL)SPINNER_az_path->enable();
    if(SPINNER_elev_path!=NULL)SPINNER_elev_path->enable();
  }
  else{
    if(SPINNER_t!=NULL)SPINNER_t->disable();
    if(ROLLOUT_keyframe!=NULL&&ROLLOUT_keyframe->enabled==1)ROLLOUT_keyframe->disable();
    if(SPINNER_az_path!=NULL)SPINNER_az_path->disable();
    if(SPINNER_elev_path!=NULL)SPINNER_elev_path->disable();
  }

  if(CHECKBOX_tourhide!=NULL){
    if(viewanytours>0&&edittour==1){
      CHECKBOX_tourhide->enable();
    }
    else{
      CHECKBOX_tourhide->disable();
    }
  }
  if(selected_tour!=NULL){
    selectedtour_index = selected_tour-tourinfo;
    LISTBOX_tour->set_int_val(selectedtour_index);
    LISTBOX_avatar->enable();
    CHECKBOX_showtour_locus->enable();
  }
  else{
    selectedtour_index = TOURINDEX_MANUAL;
    LISTBOX_tour->set_int_val(selectedtour_index);
    LISTBOX_avatar->disable();
    CHECKBOX_showtour_locus->disable();
  }
  if(edittour==1){
    if(tour_constant_vel==0){
      SPINNER_t->enable();
    }
    else{
      SPINNER_t->disable();
    }
  }

}

/* ------------------ update_globaltension ------------------------ */

extern "C" void update_globaltension(void){
  TOUR_CB(GLOBAL_TENSIONFLAG);
}
