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
extern "C" char glui_tour_revision[]="$Revision$";

static int viewtype=0;
static float tour_x=0.0, tour_y=0.0, tour_z=0.0, tour_ttt, tour_az_path=0.0, tour_tension=0.0;
#ifdef pp_TOUR
static float tour_az_scene=0.0,tour_elev_scene=0.0;
static int viewtype1=0,viewtype2=1, viewtype3=0;
#endif
static float tour_viewx=0.0, tour_viewy=0.0, tour_viewz=0.0, tour_elev_path=0.0;
static int tour_hide=0;
static int tour_global_tension_flag=1;
static float tour_global_tension=0.0;
static float tour_zoom=1.0;
static char tour_label[sizeof(GLUI_String)];

int nexttour(void);
int prevtour(void);
void TOUR_CB(int var);
void setviewcontrols(void);

GLUI_Spinner *SPINNER_t=NULL,*SPINNER_x=NULL, *SPINNER_y=NULL,*SPINNER_z=NULL;
GLUI_Spinner *SPINNER_viewx=NULL, *SPINNER_viewy=NULL,*SPINNER_viewz=NULL;
GLUI_Spinner *SPINNER_az_path=NULL;
#ifdef pp_TOUR
GLUI_Spinner *SPINNER_az_scene=NULL;
GLUI_Spinner *SPINNER_elev_scene=NULL;
#endif

GLUI *glui_advancedtour=NULL, *glui_tour=NULL;
GLUI_Rollout *panel_tour=NULL;
#ifdef pp_TOUR
GLUI_Rollout *panel_view1=NULL;
GLUI_Rollout *panel_view2=NULL;
GLUI_Rollout *panel_view3=NULL;
#endif
GLUI_Panel *panel_path=NULL;
GLUI_Panel *panel_advancedkeyframe=NULL;
GLUI_Panel *panel_keyframe=NULL;
GLUI_Panel *panel_tour1=NULL;
GLUI_Panel *panel_tour2=NULL;
GLUI_Panel *panel_tour4=NULL;
GLUI_Panel *panel_view=NULL;
GLUI_Panel *panel_advanced=NULL;
GLUI_Panel *panel_spline=NULL;
GLUI_Panel *panel_settings=NULL;
GLUI_Panel *panel_movedir=NULL;
GLUI_Checkbox *snap_checkbox=NULL,*view_checkbox=NULL,*showtourroute_checkbox=NULL,*CHECKBOX_constantvel=NULL;
GLUI_Checkbox *CHECKBOXshowtourlocus=NULL,*CHECKBOXshowintermediate=NULL;
GLUI_Checkbox *CHECKBOXglobaltension_flag=NULL, *CHECKBOXtourhide=NULL;
GLUI_Checkbox *CHECKBOXview1=NULL;
#ifdef pp_TOUR
GLUI_Checkbox *CHECKBOXview2=NULL;
GLUI_Checkbox *CHECKBOXview3=NULL;
#endif
GLUI_Checkbox *CHECKBOX_usecurrent=NULL;
GLUI_Spinner *SPINNER_globaltourtension=NULL;
GLUI_Spinner *SPINNER_tourtension=NULL;
GLUI_Spinner *SPINNER_tourbias=NULL;
GLUI_Spinner *SPINNER_tourcontinuity=NULL;
//GLUI_Spinner *SPINNER_tourzoom=NULL;
GLUI_Spinner *SPINNER_elev_path=NULL;
GLUI_Button *BUTTON_settings=NULL,*BUTTONnext_tour=NULL,*BUTTONprev_tour=NULL;
GLUI_EditText *EDITlabel=NULL;
GLUI_Listbox *LISTBOX_tour=NULL;

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
#define ADVANCEDTOUR_CLOSE 101
#define ADVANCEDTOUR_OPEN 102
#define TOUR_HIDE 20
#define KEYFRAME_viewXYZ 22
#define VIEWSNAP 23
#define TOUR_LIST 24
#define VIEW1 26
#ifdef pp_TOUR
#define VIEW2 30
#define VIEW3 31
#endif
#define VIEW_times 27
#define TOUR_UPDATELABEL 28
#define TOUR_USECURRENT 29

#define TOURMENU(f) callfrom_tourglui=1;TourMenu(f);callfrom_tourglui=0;

extern "C" void add_new_tour(void){
  TOUR_CB(TOUR_INSERT);
}

extern "C" void glui_tour_setup(int main_window){

  tourdata *touri;
  int i;

  if(glui_tour!=NULL)glui_tour->close();
  glui_tour = GLUI_Master.create_glui("Edit Tours",0,0,0);
  if(showgluitour==0)glui_tour->hide();


  panel_settings = glui_tour->add_panel("Tour Settings");
  view_checkbox=glui_tour->add_checkbox_to_panel(panel_settings,"View From Tour Path",&viewtourfrompath,VIEWTOURFROMPATH,TOUR_CB);
  snap_checkbox=glui_tour->add_checkbox_to_panel(panel_settings,"View From Selected Keyframe",&keyframe_snap,VIEWSNAP,TOUR_CB);
  glui_tour->add_separator_to_panel(panel_settings);
  CHECKBOX_constantvel=glui_tour->add_checkbox_to_panel(panel_settings,"Constant Speed",&tour_constant_vel,CONSTANTTOURVEL,TOUR_CB);
  showtourroute_checkbox=glui_tour->add_checkbox_to_panel(panel_settings,"Edit Tour Path",&edittour,SHOWTOURROUTE,TOUR_CB);
  CHECKBOXshowtourlocus=glui_tour->add_checkbox_to_panel(panel_settings,"Show Avatar",&show_tourlocus);
  CHECKBOXshowintermediate=glui_tour->add_checkbox_to_panel(panel_settings,"Show Intermediate Path Nodes",&show_path_knots);

  glui_tour->add_separator_to_panel(panel_settings);

  glui_tour->add_button_to_panel(panel_settings,"New Tour",TOUR_INSERT,TOUR_CB);
  BUTTONnext_tour=glui_tour->add_button_to_panel(panel_settings,"Next Tour",TOUR_NEXT,TOUR_CB);
  BUTTONprev_tour=glui_tour->add_button_to_panel(panel_settings,"Previous Tour",TOUR_PREVIOUS,TOUR_CB);
  if(ntours>0){
    selectedtour_index=-1;
    selectedtour_index_old=-1;
    LISTBOX_tour=glui_tour->add_listbox_to_panel(panel_settings,"Select Tour:",&selectedtour_index,TOUR_LIST,TOUR_CB);

    LISTBOX_tour->add_item(-1,"Manual");
    LISTBOX_tour->add_item(-999,"-");
    for(i=0;i<ntours;i++){
      touri = tourinfo + i;
      LISTBOX_tour->add_item(i,touri->label);
    }
    LISTBOX_tour->set_int_val(selectedtour_index);
  }
  glui_tour->add_separator_to_panel(panel_settings);

  glui_tour->add_button_to_panel(panel_settings,"Update Tour Label",TOUR_UPDATELABEL,TOUR_CB);
  EDITlabel=glui_tour->add_edittext_to_panel(panel_settings,"Tour Label",GLUI_EDITTEXT_TEXT,tour_label,TOUR_LABEL,TOUR_CB);
  EDITlabel->set_w(240);

  panel_keyframe = glui_tour->add_panel("Keyframe");
  
#ifdef pp_TOUR
  CHECKBOX_usecurrent=glui_tour->add_checkbox_to_panel(panel_keyframe,"Add keyframe using current position",
    &tour_usecurrent,TOUR_USECURRENT,TOUR_CB);
#endif
  panel_tour1 = glui_tour->add_panel_to_panel(panel_keyframe,"",GLUI_PANEL_NONE);
  glui_tour->add_button_to_panel(panel_tour1,"Next",    KEYFRAME_NEXT,TOUR_CB);
  glui_tour->add_button_to_panel(panel_tour1,"Previous",KEYFRAME_PREVIOUS,TOUR_CB);
  SPINNER_t=glui_tour->add_spinner_to_panel(panel_tour1,"t:",GLUI_SPINNER_FLOAT,&tour_ttt,KEYFRAME_tXYZ,TOUR_CB);
  glui_tour->add_column_to_panel(panel_tour1,false);
  glui_tour->add_button_to_panel(panel_tour1,"Add",KEYFRAME_INSERT,TOUR_CB);
  glui_tour->add_button_to_panel(panel_tour1,"Delete",KEYFRAME_DELETE,TOUR_CB);
#ifndef pp_TOUR
  CHECKBOXview1=glui_tour->add_checkbox_to_panel(panel_tour1,"X,Y,Z View",&viewtype,VIEW1,TOUR_CB);
#endif

  panel_movedir = glui_tour->add_panel_to_panel(panel_keyframe,"",GLUI_PANEL_NONE);
  glui_tour->add_statictext_to_panel(panel_movedir,"Position");
  SPINNER_x=glui_tour->add_spinner_to_panel(panel_movedir,"X:",GLUI_SPINNER_FLOAT,&tour_x,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_y=glui_tour->add_spinner_to_panel(panel_movedir,"Y:",GLUI_SPINNER_FLOAT,&tour_y,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_z=glui_tour->add_spinner_to_panel(panel_movedir,"Z:",GLUI_SPINNER_FLOAT,&tour_z,KEYFRAME_tXYZ,TOUR_CB);

  glui_tour->add_column_to_panel(panel_movedir,false);
  glui_tour->add_statictext_to_panel(panel_movedir,"View Direction");
#ifdef pp_TOUR

  CHECKBOXview1=glui_tour->add_checkbox_to_panel(panel_movedir,"Position",&viewtype1,VIEW1,TOUR_CB);
  CHECKBOXview2=glui_tour->add_checkbox_to_panel(panel_movedir,"Path relative",&viewtype2,VIEW2,TOUR_CB);
  CHECKBOXview3=glui_tour->add_checkbox_to_panel(panel_movedir,"Scene relative",&viewtype3,VIEW3,TOUR_CB);

  panel_view1 = glui_tour->add_rollout_to_panel(panel_movedir,"Position",false);
  SPINNER_viewx=glui_tour->add_spinner_to_panel(panel_view1,"X",GLUI_SPINNER_FLOAT,&tour_viewx,KEYFRAME_viewXYZ,TOUR_CB);
  SPINNER_viewy=glui_tour->add_spinner_to_panel(panel_view1,"Y",GLUI_SPINNER_FLOAT,&tour_viewy,KEYFRAME_viewXYZ,TOUR_CB);
  SPINNER_viewz=glui_tour->add_spinner_to_panel(panel_view1,"Z",GLUI_SPINNER_FLOAT,&tour_viewz,KEYFRAME_viewXYZ,TOUR_CB);

  panel_view2 = glui_tour->add_rollout_to_panel(panel_movedir,"Path relative",false);
  SPINNER_az_path=glui_tour->add_spinner_to_panel(panel_view2,"Azimuth:",GLUI_SPINNER_FLOAT,&tour_az_path,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_elev_path=glui_tour->add_spinner_to_panel(panel_view2,"Elevation:",GLUI_SPINNER_FLOAT,&tour_elev_path,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_elev_path->set_float_limits(-90.0,90.0);

  panel_view3 = glui_tour->add_rollout_to_panel(panel_movedir,"Scene relative",false);
  SPINNER_az_scene=glui_tour->add_spinner_to_panel(panel_view3,"Azimuth:",GLUI_SPINNER_FLOAT,&tour_az_scene,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_elev_scene=glui_tour->add_spinner_to_panel(panel_view3,"Elevation:",GLUI_SPINNER_FLOAT,&tour_elev_scene,KEYFRAME_tXYZ,TOUR_CB);
#else
  SPINNER_viewx=glui_tour->add_spinner_to_panel(panel_movedir,"X",GLUI_SPINNER_FLOAT,&tour_viewx,KEYFRAME_viewXYZ,TOUR_CB);
  SPINNER_viewy=glui_tour->add_spinner_to_panel(panel_movedir,"Y",GLUI_SPINNER_FLOAT,&tour_viewy,KEYFRAME_viewXYZ,TOUR_CB);
  SPINNER_viewz=glui_tour->add_spinner_to_panel(panel_movedir,"Z",GLUI_SPINNER_FLOAT,&tour_viewz,KEYFRAME_viewXYZ,TOUR_CB);
  SPINNER_az_path=glui_tour->add_spinner_to_panel(panel_movedir,"Azimuth:",GLUI_SPINNER_FLOAT,&tour_az_path,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_elev_path=glui_tour->add_spinner_to_panel(panel_movedir,"Elevation:",GLUI_SPINNER_FLOAT,&tour_elev_path,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_elev_path->set_float_limits(-90.0,90.0);
#endif

  
 // if(tour_constant_vel==0){
 //   SPINNER_t->enable();
 // }
 // else{
 //   SPINNER_t->disable();
 // }

 
  glui_tour->add_separator();

  panel_tour2 = glui_tour->add_panel("",false);
  glui_tour->add_button_to_panel(panel_tour2,"Advanced Settings",ADVANCEDTOUR_OPEN,TOUR_CB);
  glui_tour->add_column_to_panel(panel_tour2,false);
  BUTTON_settings=glui_tour->add_button_to_panel(panel_tour2,"Save Settings",SAVE_SETTINGS,TOUR_CB);

  glui_tour->add_button("Close",TOUR_CLOSE,TOUR_CB);

  glui_tour->set_main_gfx_window( main_window );
  update_tourcontrols();
}

/* ------------------ glui_advancedtour_setup ------------------------ */

extern "C" void glui_advancedtour_setup(int main_window){

  if(glui_advancedtour!=NULL)glui_advancedtour->close();
  glui_advancedtour = GLUI_Master.create_glui("Advanced Settings",0,0,0);
  glui_advancedtour->hide();



  panel_advanced = glui_advancedtour->add_panel("Advanced Settings",false);

  panel_path = glui_advancedtour->add_panel_to_panel(panel_advanced,"Duration/Points (all tours)",true);

  glui_advancedtour->add_spinner_to_panel(panel_path,"start time",GLUI_SPINNER_FLOAT,&view_tstart,VIEW_times,TOUR_CB);
  glui_advancedtour->add_spinner_to_panel(panel_path,"stop time:",GLUI_SPINNER_FLOAT,&view_tstop, VIEW_times,TOUR_CB);
  glui_advancedtour->add_spinner_to_panel(panel_path,"points",    GLUI_SPINNER_INT,&view_ntimes,  VIEW_times,TOUR_CB);

  panel_advancedkeyframe = glui_advancedtour->add_panel_to_panel(panel_advanced,"Edit Keyframe's View/Tension",true);

  glui_advancedtour->add_button_to_panel(panel_advancedkeyframe,"Next",    KEYFRAME_NEXT,TOUR_CB);
  glui_advancedtour->add_button_to_panel(panel_advancedkeyframe,"Previous",KEYFRAME_PREVIOUS,TOUR_CB);

  panel_spline = glui_advancedtour->add_panel_to_panel(panel_advancedkeyframe,"Tension",true);
  CHECKBOXglobaltension_flag=glui_advancedtour->add_checkbox_to_panel(panel_spline,"Global",&tour_global_tension_flag,GLOBAL_TENSIONFLAG,TOUR_CB);
   SPINNER_globaltourtension=glui_advancedtour->add_spinner_to_panel(panel_spline,"All keyframes",
    GLUI_SPINNER_FLOAT,&tour_global_tension,GLOBAL_TENSION,TOUR_CB);
  SPINNER_globaltourtension->set_float_limits(-1.0,1.0,GLUI_LIMIT_CLAMP);

  SPINNER_tourtension=glui_advancedtour->add_spinner_to_panel(panel_spline,"Selected keyframeframe",
    GLUI_SPINNER_FLOAT,&tour_tension,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_tourtension->set_float_limits(-1.0,1.0,GLUI_LIMIT_CLAMP);

  TOUR_CB(VIEW1);
  TOUR_CB(GLOBAL_TENSIONFLAG);

  SPINNER_az_path->set_float_limits(-180.0,180.0);
  
  glui_advancedtour->add_button("Close",ADVANCEDTOUR_CLOSE,TOUR_CB);

  glui_advancedtour->set_main_gfx_window( main_window );
  update_tourcontrols();
  update_tourlist=1;
 // TOUR_CB(TOUR_LIST);
}

extern "C" void Update_Tourlist(void){

  update_tourlist=0;
  TOUR_CB(TOUR_LIST);
}
/* ------------------ hide_glui_advancedtour ------------------------ */

extern "C" void hide_glui_advancedtour(void){
  if(glui_advancedtour!=NULL)glui_advancedtour->hide();
}

/* ------------------ hide_glui_tour ------------------------ */

extern "C" void hide_glui_tour(void){
  if(glui_tour!=NULL)glui_tour->hide();
  hide_glui_advancedtour();
  showgluitour=0;
  updatemenu=1;
}

/* ------------------ show_glui_tour ------------------------ */

extern "C" void show_glui_tour(void){
  if(glui_tour!=NULL)glui_tour->show();
}

/* ------------------ show_glui_advancedtour ------------------------ */

extern "C" void show_glui_advancedtour(void){
  if(glui_tour!=NULL)glui_advancedtour->show();
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

/* ------------------ set_glui_keyframe ------------------------ */

extern "C" void set_glui_keyframe(){
  tourdata *ti;
  float *eye,*aview;

  if(selected_frame==NULL)return;

  ti = selected_tour;
  if(ti!=NULL){
    tour_hide=1-ti->display;
  }
  if(selected_tour!=NULL)strcpy(tour_label,selected_tour->label);

  eye = selected_frame->nodeval.eye;
  aview = selected_frame->nodeval.aview;

  tour_ttt = selected_frame->disp_time;
  tour_x = trim_val(xbar0 + eye[0]*xyzmaxdiff);
  tour_y = trim_val(ybar0 + eye[1]*xyzmaxdiff);
  tour_z = trim_val(zbar0 + eye[2]*xyzmaxdiff);
  tour_viewx = trim_val(xbar0 + xyzmaxdiff*aview[0]);
  tour_viewy = trim_val(ybar0 + xyzmaxdiff*aview[1]);
  tour_viewz = trim_val(zbar0 + xyzmaxdiff*aview[2]);
  tour_az_path = selected_frame->az_path;
#ifdef pp_TOUR
  tour_az_scene = selected_frame->az_scene;
  tour_elev_scene = selected_frame->nodeval.elev_scene;
#endif
  tour_continuity=selected_frame->continuity;
  tour_bias=selected_frame->bias;
#ifdef pp_TOUR
  viewtype1=selected_frame->viewtype1;
  viewtype2=selected_frame->viewtype2;
  viewtype3=selected_frame->viewtype3;
  setviewcontrols();
#else
  viewtype=selected_frame->viewtype;
#endif
  tour_tension=selected_frame->tension;
  tour_zoom=selected_frame->nodeval.zoom;
  tour_elev_path=selected_frame->nodeval.elev_path;
 
  tour_global_tension_flag=selected_tour->global_tension_flag;
  tour_global_tension=selected_tour->global_tension;


  if(SPINNER_t==NULL)return;

  CHECKBOXglobaltension_flag->set_int_val(tour_global_tension_flag);
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

  SPINNER_x->set_float_val(tour_x);
  SPINNER_y->set_float_val(tour_y);
  SPINNER_z->set_float_val(tour_z);
  if(SPINNER_tourbias!=NULL)SPINNER_tourbias->set_float_val(tour_bias);
  if(SPINNER_tourcontinuity!=NULL)SPINNER_tourcontinuity->set_float_val(tour_continuity);
  SPINNER_tourtension->set_float_val(tour_tension);
//    SPINNER_tourzoom->set_float_val(tour_zoom);
  SPINNER_viewx->set_float_val(tour_viewx);
  SPINNER_viewy->set_float_val(tour_viewy);
  SPINNER_viewz->set_float_val(tour_viewz);
  SPINNER_az_path->set_float_val(tour_az_path);
  SPINNER_elev_path->set_float_val(tour_elev_path);
#ifdef pp_TOUR
  SPINNER_az_scene->set_float_val(tour_az_scene);
  SPINNER_elev_scene->set_float_val(tour_elev_scene);
#endif
  if(ti!=NULL&&CHECKBOXtourhide!=NULL)CHECKBOXtourhide->set_int_val(tour_hide);
  EDITlabel->set_text(tour_label);

  if(edittour==1){
#ifdef pp_TOUR
    setviewcontrols();
#else
    if(viewtype==1){
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
#endif
  }
#ifdef pp_TOUR
  CHECKBOXview1->set_int_val(viewtype1);
#else
  CHECKBOXview1->set_int_val(viewtype);
#endif
}

extern "C" void update_tourindex(void){
  update_selectedtour_index=0;
  selectedtour_index=selectedtour_index_ini;
  TOUR_CB(TOUR_LIST);
}

/* ------------------ TOUR_CB ------------------------ */

void TOUR_CB(int var){
  keyframe *thiskey,*nextkey,*newframe;
#ifndef pp_TOUR
  keyframe *lastkey;
  float dummy;
#endif
  tourdata *thistour=NULL;
  float *aview,*eye;

  float key_xyz[3];
  float key_params[3];
  float key_time_in, key_az_path, key_view[3], key_zoom;
  float key_elev_path, key_bank;
#ifdef pp_TOUR
  float key_az_scene, key_elev_scene;
#endif

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

  switch (var){
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
  case ADVANCEDTOUR_OPEN:
    show_glui_advancedtour();
    break;
  case ADVANCEDTOUR_CLOSE:
    hide_glui_advancedtour();
    break;
  case TOUR_CLOSE:
    hide_glui_tour();
    TOUR_CB(ADVANCEDTOUR_CLOSE);
    break;
  case SAVE_SETTINGS:
    writeini(LOCAL_INI);
    break;
  case SHOWTOURROUTE:
    edittour = 1 - edittour;
    TOURMENU(-4);
    /*
    if(panel_keyframe!=NULL){
      if(edittour==1){
        panel_keyframe->enable();
      }
      else{
        panel_keyframe->disable();
      }
    }
    if(panel_advanced!=NULL){
      if(edittour==1){
        panel_advanced->enable();
      }
      else{
        panel_advanced->disable();
      }
    }
    */
    update_tourcontrols();
    TOUR_CB(VIEW1);
    break;
  case VIEWSNAP:
    if(viewtourfrompath==1&&keyframe_snap==1){
      viewtourfrompath=0;
      view_checkbox->set_int_val(viewtourfrompath);
    }
    TOUR_CB(VIEWTOURFROMPATH);
    break;
  case VIEWTOURFROMPATH:
    if(viewtourfrompath==1&&keyframe_snap==1){
      keyframe_snap=0;
      snap_checkbox->set_int_val(keyframe_snap);
    }
    viewtourfrompath = 1 - viewtourfrompath;
    TOURMENU(-5);
    break;
  case VIEW1:
#ifdef pp_TOUR
    if(viewtype1==1){
      viewtype2=0;
      viewtype3=0;
    }
    else{
      viewtype1=0;
      viewtype2=1;
      viewtype3=0;
    }
    setviewcontrols();
#else
    if(viewtype==1){
      SPINNER_az_path->disable();
      SPINNER_elev_path->disable();
      SPINNER_viewx->enable();
      SPINNER_viewy->enable();
      SPINNER_viewz->enable();
      if(selected_frame!=NULL){
        adjustviewangle(selected_frame,&dummy,&dummy);
        SPINNER_az_path->set_float_val(tour_az_path);
        SPINNER_elev_path->set_float_val(tour_elev_path);
      }
    }
    else if(viewtype==0){
      SPINNER_az_path->enable();
      SPINNER_elev_path->enable();
      SPINNER_viewx->disable();
      SPINNER_viewy->disable();
      SPINNER_viewz->disable();
    }
#endif
    if(selected_frame!=NULL){
#ifdef pp_TOUR
      selected_frame->viewtype1=viewtype1;
      selected_frame->viewtype2=viewtype2;
      selected_frame->viewtype3=viewtype3;
      setviewcontrols();
#endif
      selected_frame->viewtype=viewtype;
    }
    createtourpaths();
    break;
#ifdef pp_TOUR
  case VIEW2:
    if(viewtype2==1){
      viewtype1=0;
      viewtype3=0;
    }
    else{
      viewtype1=0;
      viewtype2=0;
      viewtype3=1;
    }
    if(selected_frame!=NULL){
      selected_frame->viewtype1=viewtype1;
      selected_frame->viewtype2=viewtype2;
      selected_frame->viewtype3=viewtype3;
    }
    setviewcontrols();
    selected_frame->viewtype=viewtype;
    createtourpaths();
    break;
  case VIEW3:
    if(viewtype3==1){
      viewtype1=0;
      viewtype2=0;
    }
    else{
      viewtype1=0;
      viewtype2=1;
      viewtype3=0;
    }
    if(selected_frame!=NULL){
      selected_frame->viewtype1=viewtype1;
      selected_frame->viewtype2=viewtype2;
      selected_frame->viewtype3=viewtype3;
    }
    setviewcontrols();
    selected_frame->viewtype=viewtype;
    createtourpaths();
    break;
#endif
  case VIEW_times:
    ReallocTourMemory();
    createtourpaths();
    updatetimes();
    break;
  case KEYFRAME_viewXYZ:
    if(selected_frame!=NULL){
      if(selected_tour-tourinfo==0)dirtycircletour=1;
      selected_tour->startup=0;
      aview = selected_frame->nodeval.aview;
      aview[0]=(tour_viewx-xbar0)/xyzmaxdiff;
      aview[1]=(tour_viewy-ybar0)/xyzmaxdiff;
      aview[2]=(tour_viewz-zbar0)/xyzmaxdiff;

#ifdef pp_TOUR
      setviewcontrols();
//      adjustviewangle(selected_frame,viewtype);
//      SPINNER_az_path->set_float_val(tour_az_path);
//      SPINNER_elev_path->set_float_val(tour_elev_path);
#else
      adjustviewangle(selected_frame,&tour_az_path,&tour_elev_path);
      SPINNER_az_path->set_float_val(tour_az_path);
      SPINNER_elev_path->set_float_val(tour_elev_path);
#endif

      createtourpaths();
      selected_frame->selected=1;
    }
    break;
  case KEYFRAME_tXYZ:
    if(selected_frame!=NULL){
      if(selected_tour-tourinfo==0)dirtycircletour=1;
      selected_tour->startup=0;
      eye = selected_frame->nodeval.eye;
      aview = selected_frame->nodeval.aview;

      /*
      {
        int change_time;

        change_time=0;
        if(tour_ttt<selected_frame->prev->display_time){
          change_time=1;
          tour_ttt=selected_frame->prev->display_time;
        }
        if(tour_ttt>selected_frame->next->display_time){
          change_time=1;
          tour_ttt=selected_frame->next->display_time;
        }
        if(change_time==1){
          SPINNER_t->set_float_val(tour_ttt);
        }
      }

        */

      if(tour_constant_vel==0){
        selected_frame->noncon_time=tour_ttt;
        selected_frame->disp_time=tour_ttt;
      }
      eye[0]=(tour_x-xbar0)/xyzmaxdiff;
      eye[1]=(tour_y-ybar0)/xyzmaxdiff;
      eye[2]=(tour_z-zbar0)/xyzmaxdiff;
      if(viewtype==0){
        tour_az_path = SPINNER_az_path->get_float_val();
      }
      selected_frame->az_path=tour_az_path;
      selected_frame->nodeval.elev_path=tour_elev_path;
#ifdef pp_TOUR
      selected_frame->az_scene=tour_az_scene;
      selected_frame->nodeval.elev_scene=tour_elev_scene;
#endif

      selected_frame->tension=tour_tension;
      selected_frame->bias=tour_bias;
      selected_frame->continuity=tour_continuity;
      selected_frame->viewtype=viewtype;
      selected_frame->nodeval.zoom=tour_zoom;
      aview[0]=(tour_viewx-xbar0)/xyzmaxdiff;
      aview[1]=(tour_viewy-ybar0)/xyzmaxdiff;
      aview[2]=(tour_viewz-zbar0)/xyzmaxdiff;
      createtourpaths();
      selected_frame->selected=1;
      if(viewtype==1){
        TOUR_CB(KEYFRAME_viewXYZ);
      }
#ifdef pp_TOUR
      setviewcontrols();
//      selected_frame->az_scene=tour_az_scene;
//      selected_frame->nodeval.elev_scene=tour_elev_scene;
#endif
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
    updatetimes();
    set_glui_keyframe();
    break;
#ifdef pp_TOUR
  case KEYFRAME_INSERT:
    thistour=selected_tour;
    if(times!=NULL&&tour_usecurrent==1){
      keyframe *keyj;
      keyframe *first_frame,*last_frame;
      float time0;

      time0 = timeoffset + times[itime];
      first_frame = (thistour->first_frame).next;
      last_frame = (thistour->last_frame).prev;
      if(first_frame->next!=NULL&&last_frame!=NULL){
        thiskey=first_frame;
        if(time0>=first_frame->noncon_time&&time0<last_frame->noncon_time){
          for(keyj=(thistour->first_frame).next;keyj->next!=NULL;keyj=keyj->next){
            if(keyj->noncon_time<=time0&&time0<keyj->next->noncon_time){
              thiskey=keyj;
              break;
            }
          }
        }
        else if(time0>last_frame->noncon_time){
          thiskey=last_frame;
        }
      }
      else{
        thiskey=selected_frame;
      }
    }
    else{
      thiskey=selected_frame;
    }
    if(thiskey==NULL)return;
    if(thiskey->next!=NULL){
      nextkey=thiskey->next;
      if(thiskey->viewtype!=nextkey->viewtype)nextkey=thiskey;
    }
    else{
      nextkey=thiskey;
    }
    if(times!=NULL&&tour_usecurrent==1){
      float *xyz, *angles;
      float time0;
      float az2, el2;

      time0 = timeoffset + times[itime];

      xyz = camera_current->eye;
      angles = camera_current->angle_zx;

      key_xyz[0]=xbar0 + xyzmaxdiff*xyz[0];
      key_xyz[1]=ybar0 + xyzmaxdiff*xyz[1];
      key_xyz[2]=zbar0 + xyzmaxdiff*xyz[2];
      key_az_scene = camera_current->direction_angle + camera_current->view_angle;
      key_elev_scene=camera_current->elevation_angle;
      az2 = PI*key_az_scene/180.0;
      el2 = PI*key_elev_scene/180.0;
      key_view[0]=key_xyz[0] + sin(az2)*cos(el2);
      key_view[1]=key_xyz[1] + cos(az2)*cos(el2);
      key_view[2]=key_xyz[2] + sin(el2);
      key_time_in = time0;
      viewtype=2;
      printf(" tour time=%f xyz=%f %f %f angles=%f %f\n",key_time_in,key_xyz[0],key_xyz[1],key_xyz[2],
        key_az_scene,key_elev_scene);
    }
    else{
      key_xyz[0]=xbar0 + xyzmaxdiff*(thiskey->nodeval.eye[0]+nextkey->nodeval.eye[0])/2.0;
      key_xyz[1]=ybar0 + xyzmaxdiff*(thiskey->nodeval.eye[1]+nextkey->nodeval.eye[1])/2.0;
      key_xyz[2]=zbar0 + xyzmaxdiff*(thiskey->nodeval.eye[2]+nextkey->nodeval.eye[2])/2.0;
      key_az_path = (thiskey->az_path+nextkey->az_path)/2.0;
      key_elev_path=(thiskey->nodeval.elev_path+nextkey->nodeval.elev_path)/2.0;
      key_az_scene = (thiskey->az_scene+nextkey->az_scene)/2.0;
      key_elev_scene=(thiskey->nodeval.elev_scene+nextkey->nodeval.elev_scene)/2.0;
      key_time_in = (thiskey->noncon_time+nextkey->noncon_time)/2.0;
      key_view[0]=xbar0 + xyzmaxdiff*(thiskey->nodeval.aview[0]+nextkey->nodeval.aview[0])/2.0;
      key_view[1]=ybar0 + xyzmaxdiff*(thiskey->nodeval.aview[1]+nextkey->nodeval.aview[1])/2.0;
      key_view[2]=zbar0 + xyzmaxdiff*(thiskey->nodeval.aview[2]+nextkey->nodeval.aview[2])/2.0;
      viewtype=thiskey->viewtype;
    }

    key_params[0]=(thiskey->bias+nextkey->bias)/2.0;
    key_params[1]=(thiskey->continuity+nextkey->continuity)/2.0;
    key_params[2]=(thiskey->tension+nextkey->tension)/2.0;
    key_zoom = (thiskey->nodeval.zoom + nextkey->nodeval.zoom)/2.0;
    key_bank = (thiskey->bank + nextkey->bank)/2.0;

    newframe=add_frame(selected_frame,key_time_in,key_xyz,key_az_path,key_az_scene,key_elev_path,key_elev_scene,key_bank,
    key_params,viewtype,key_zoom,key_view);
    new_select(newframe);
    set_glui_keyframe();
    adjustviewangle(newframe,2);
    createtourpaths();

    break;
#else
  case KEYFRAME_INSERT:
    if(selected_frame!=NULL){
      thistour=selected_tour;
      thiskey=selected_frame;
      nextkey=thiskey->next;
      if(nextkey==&thistour->last_frame){
        lastkey=thiskey->prev;
        key_xyz[0]=xbar0 + xyzmaxdiff*(2*thiskey->nodeval.eye[0]-lastkey->nodeval.eye[0]);
        key_xyz[1]=ybar0 + xyzmaxdiff*(2*thiskey->nodeval.eye[1]-lastkey->nodeval.eye[1]);
        key_xyz[2]=zbar0 + xyzmaxdiff*(2*thiskey->nodeval.eye[2]-lastkey->nodeval.eye[2]);
        key_az_path = (2*thiskey->az_path-lastkey->az_path);
        key_elev_path=(2*thiskey->nodeval.elev_path-lastkey->nodeval.elev_path);
        key_time_in = thiskey->noncon_time;
        thiskey->noncon_time=(thiskey->noncon_time+lastkey->noncon_time)/2.0;
        key_params[0]=(2*thiskey->bias-lastkey->bias);
        key_params[1]=(2*thiskey->continuity-lastkey->continuity);
        key_params[2]=(2*thiskey->tension-lastkey->tension);
        key_view[0]=xbar0 + xyzmaxdiff*(2*thiskey->nodeval.aview[0]-lastkey->nodeval.aview[0]);
        key_view[1]=ybar0 + xyzmaxdiff*(2*thiskey->nodeval.aview[1]-lastkey->nodeval.aview[1]);
        key_view[2]=zbar0 + xyzmaxdiff*(2*thiskey->nodeval.aview[2]-lastkey->nodeval.aview[2]);
        key_zoom = (2*thiskey->nodeval.zoom - lastkey->nodeval.zoom);
        key_bank = (2*thiskey->bank - lastkey->bank);
        viewtype=thiskey->viewtype;
      }
      else{
        key_xyz[0]=xbar0 + xyzmaxdiff*(thiskey->nodeval.eye[0]+nextkey->nodeval.eye[0])/2.0;
        key_xyz[1]=ybar0 + xyzmaxdiff*(thiskey->nodeval.eye[1]+nextkey->nodeval.eye[1])/2.0;
        key_xyz[2]=zbar0 + xyzmaxdiff*(thiskey->nodeval.eye[2]+nextkey->nodeval.eye[2])/2.0;
        key_az_path = (thiskey->az_path+nextkey->az_path)/2.0;
        key_elev_path=(thiskey->nodeval.elev_path+nextkey->nodeval.elev_path)/2.0;
        key_time_in = (thiskey->noncon_time+nextkey->noncon_time)/2.0;
        key_params[0]=(thiskey->bias+nextkey->bias)/2.0;
        key_params[1]=(thiskey->continuity+nextkey->continuity)/2.0;
        key_params[2]=(thiskey->tension+nextkey->tension)/2.0;
        key_view[0]=xbar0 + xyzmaxdiff*(thiskey->nodeval.aview[0]+nextkey->nodeval.aview[0])/2.0;
        key_view[1]=ybar0 + xyzmaxdiff*(thiskey->nodeval.aview[1]+nextkey->nodeval.aview[1])/2.0;
        key_view[2]=zbar0 + xyzmaxdiff*(thiskey->nodeval.aview[2]+nextkey->nodeval.aview[2])/2.0;
        key_zoom = (thiskey->nodeval.zoom + nextkey->nodeval.zoom)/2.0;
        key_bank = (thiskey->bank + nextkey->bank)/2.0;
        if(thiskey->viewtype==0&&nextkey->viewtype==0){
          viewtype=0;
        }
        else{
          viewtype=1;
          if(thiskey->viewtype==1){
            key_view[0]=xbar0 + xyzmaxdiff*thiskey->nodeval.aview[0];
            key_view[1]=ybar0 + xyzmaxdiff*thiskey->nodeval.aview[1];
            key_view[2]=zbar0 + xyzmaxdiff*thiskey->nodeval.aview[2];
            key_elev_path = thiskey->nodeval.elev_path;
          }
          if(thiskey->viewtype==0&&nextkey->viewtype==1){
            key_view[0]=xbar0 + xyzmaxdiff*nextkey->nodeval.aview[0];
            key_view[1]=ybar0 + xyzmaxdiff*nextkey->nodeval.aview[1];
            key_view[2]=zbar0 + xyzmaxdiff*nextkey->nodeval.aview[2];
            key_elev_path = nextkey->nodeval.elev_path;
          }
        }
      }
      newframe=add_frame(selected_frame,key_time_in,key_xyz,key_az_path,key_elev_path,key_bank,
      key_params,viewtype,key_zoom,key_view);
      createtourpaths();
      new_select(newframe);
      set_glui_keyframe();
    }
    break;
#endif
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
  case TOUR_UPDATELABEL:
    // supposed to fall through to TOUR_LIST
  case TOUR_LIST:
    if(selectedtour_index==-999){
      selectedtour_index=selectedtour_index_old;
      if(selectedtour_index==-999)selectedtour_index=-1;
      TOUR_CB(TOUR_LIST);
      return;
    }
    switch (selectedtour_index){
    case -3:
      TOURMENU(-3); // show all tours
      set_glui_keyframe();
      break;
    case -1:
      edittour=0;
      TOURMENU(-2);  // hide all tours
      break;
    case -4:
      TOURMENU(-1);  // default tour
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
    if(showtourroute_checkbox!=NULL&&edittour==0)showtourroute_checkbox->set_int_val(1);
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
      updatetourmenulabels();
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

#ifdef pp_TOUR
/* ------------------ setviewcontrols ------------------------ */

void setviewcontrols(void){
  if(viewtype1!=CHECKBOXview1->get_int_val())CHECKBOXview1->set_int_val(viewtype1);
  if(viewtype2!=CHECKBOXview2->get_int_val())CHECKBOXview2->set_int_val(viewtype2);
  if(viewtype3!=CHECKBOXview3->get_int_val())CHECKBOXview3->set_int_val(viewtype3);

  if(viewtype1==1)viewtype=1;
  if(viewtype2==1)viewtype=0;
  if(viewtype3==1)viewtype=2;
  switch (viewtype){
    case 0:
//      panel_view1->close();
      if(panel_view2->is_open==0)panel_view2->open();
//      panel_view3->close();
      break;
    case 1:
      if(panel_view1->is_open==0)panel_view1->open();
//      panel_view2->close();
//      panel_view3->close();
      break;
    case 2:
//      panel_view1->close();
//      panel_view2->close();
      if(panel_view3->is_open==0)panel_view3->open();
      break;
  }
  if(selected_frame!=NULL){
    keyframe *sf;

    sf=selected_frame;
    adjustviewangle(sf,viewtype);
    SPINNER_az_path->set_float_val(sf->az_path);
    SPINNER_elev_path->set_float_val(sf->nodeval.elev_path);

    SPINNER_az_scene->set_float_val(sf->az_scene);
    SPINNER_elev_scene->set_float_val(sf->nodeval.elev_scene);

    SPINNER_viewx->set_float_val(xbar0+xyzmaxdiff*sf->nodeval.aview[0]);
    SPINNER_viewy->set_float_val(ybar0+xyzmaxdiff*sf->nodeval.aview[1]);
    SPINNER_viewz->set_float_val(zbar0+xyzmaxdiff*sf->nodeval.aview[2]);
  }
}
#endif

/* ------------------ delete_tourlsit ------------------------ */

void delete_tourlist(void){
  int i;
  if(LISTBOX_tour==NULL)return;
  for(i=0;i<ntours;i++){
    LISTBOX_tour->delete_item(i);
  }
}

/* ------------------ create_tourlist ------------------------ */

void create_tourlist(void){
  int i;
  tourdata *touri;
  char label[1000];
  
  if(LISTBOX_tour==NULL)return;
  for(i=0;i<ntours;i++){
    touri = tourinfo + i;
    if(i==selectedtour_index){
      strcpy(label,"*");
      strcat(label,touri->label);
      LISTBOX_tour->add_item(i,label);
    }
    else{
      LISTBOX_tour->add_item(i,touri->label);
    }
  }
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

void update_tourcontrols(void){

  if(BUTTONnext_tour==NULL)return;
  if(BUTTONprev_tour==NULL)return;
  if(panel_keyframe==NULL)return;
  if(panel_advanced==NULL)return;
  if(SPINNER_t==NULL)return;
  if(showtourroute_checkbox!=NULL)showtourroute_checkbox->set_int_val(edittour);
  if(view_checkbox!=NULL)view_checkbox->set_int_val(viewtourfrompath);
  if(ntours>1){
    BUTTONnext_tour->enable();
    BUTTONprev_tour->enable();
  }
  else{
    BUTTONnext_tour->disable();
    BUTTONprev_tour->disable();
  }
  if(ntours>0&&edittour==1){
    if(SPINNER_t!=NULL)SPINNER_t->enable();
    if(panel_keyframe!=NULL&&panel_keyframe->enabled==0)panel_keyframe->enable();
    if(panel_advanced!=NULL&&panel_advanced->enabled==0)panel_advanced->enable();
  }
  else{
    if(SPINNER_t!=NULL)SPINNER_t->disable();
    if(panel_keyframe!=NULL&&panel_keyframe->enabled==1)panel_keyframe->disable();
    if(panel_advanced!=NULL&&panel_advanced->enabled==1)panel_advanced->disable();
  }

  if(CHECKBOXtourhide!=NULL){
    if(viewanytours>0&&edittour==1){
      CHECKBOXtourhide->enable();
    }
    else{
      CHECKBOXtourhide->disable();
    }
  }
  if(selected_tour!=NULL){
    selectedtour_index = selected_tour-tourinfo;
    LISTBOX_tour->set_int_val(selectedtour_index);
  }
  else{
    selectedtour_index = -1;
    LISTBOX_tour->set_int_val(selectedtour_index);
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

void update_globaltension(void){
  TOUR_CB(GLOBAL_TENSIONFLAG);
}
