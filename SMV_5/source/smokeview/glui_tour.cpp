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
extern "C" char glui_tour_revision[]="$Revision$";

static int viewtype=0;
static float tour_x=0.0, tour_y=0.0, tour_z=0.0, tour_ttt, tour_az_path=0.0, tour_tension=0.0;
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

GLUI *glui_advancedtour=NULL, *glui_tour=NULL;
GLUI_Rollout *panel_tour=NULL;
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
GLUI_Checkbox *CHECKBOX_usecurrent=NULL;
GLUI_Spinner *SPINNER_globaltourtension=NULL;
GLUI_Spinner *SPINNER_tourtension=NULL;
GLUI_Spinner *SPINNER_tourbias=NULL;
GLUI_Spinner *SPINNER_tourcontinuity=NULL;
GLUI_Spinner *SPINNER_tourzoom=NULL;
GLUI_Spinner *SPINNER_elev_path=NULL;
GLUI_Button *BUTTON_settings=NULL,*BUTTONnext_tour=NULL,*BUTTONprev_tour=NULL;
GLUI_EditText *EDITlabel=NULL;
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
#define ADVANCEDTOUR_CLOSE 101
#define ADVANCEDTOUR_OPEN 102
#define TOUR_HIDE 20
#define KEYFRAME_viewXYZ 22
#define VIEWSNAP 23
#define TOUR_LIST 24
#define TOUR_AVATAR 31
#define VIEW1 26
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
  if(navatar_types>0){
    LISTBOX_avatar=glui_tour->add_listbox_to_panel(panel_settings,"Avatar:",&glui_avatar_index,TOUR_AVATAR,TOUR_CB);

    for(i=0;i<navatar_types;i++){
      touri = tourinfo + i;
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
  }

  panel_keyframe = glui_tour->add_panel("Keyframe");
  
  panel_tour1 = glui_tour->add_panel_to_panel(panel_keyframe,"",GLUI_PANEL_NONE);
  glui_tour->add_button_to_panel(panel_tour1,"Next",    KEYFRAME_NEXT,TOUR_CB);
  glui_tour->add_button_to_panel(panel_tour1,"Previous",KEYFRAME_PREVIOUS,TOUR_CB);
  SPINNER_t=glui_tour->add_spinner_to_panel(panel_tour1,"t:",GLUI_SPINNER_FLOAT,&tour_ttt,KEYFRAME_tXYZ,TOUR_CB);
  glui_tour->add_column_to_panel(panel_tour1,false);
  glui_tour->add_button_to_panel(panel_tour1,"Add",KEYFRAME_INSERT,TOUR_CB);
  glui_tour->add_button_to_panel(panel_tour1,"Delete",KEYFRAME_DELETE,TOUR_CB);
  CHECKBOXview1=glui_tour->add_checkbox_to_panel(panel_tour1,"X,Y,Z View",&viewtype,VIEW1,TOUR_CB);

  panel_movedir = glui_tour->add_panel_to_panel(panel_keyframe,"",GLUI_PANEL_NONE);
  glui_tour->add_statictext_to_panel(panel_movedir,"Position");
  SPINNER_x=glui_tour->add_spinner_to_panel(panel_movedir,"X:",GLUI_SPINNER_FLOAT,&tour_x,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_y=glui_tour->add_spinner_to_panel(panel_movedir,"Y:",GLUI_SPINNER_FLOAT,&tour_y,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_z=glui_tour->add_spinner_to_panel(panel_movedir,"Z:",GLUI_SPINNER_FLOAT,&tour_z,KEYFRAME_tXYZ,TOUR_CB);

  glui_tour->add_column_to_panel(panel_movedir,false);
  glui_tour->add_statictext_to_panel(panel_movedir,"View Direction");
  SPINNER_viewx=glui_tour->add_spinner_to_panel(panel_movedir,"X",GLUI_SPINNER_FLOAT,&tour_viewx,KEYFRAME_viewXYZ,TOUR_CB);
  SPINNER_viewy=glui_tour->add_spinner_to_panel(panel_movedir,"Y",GLUI_SPINNER_FLOAT,&tour_viewy,KEYFRAME_viewXYZ,TOUR_CB);
  SPINNER_viewz=glui_tour->add_spinner_to_panel(panel_movedir,"Z",GLUI_SPINNER_FLOAT,&tour_viewz,KEYFRAME_viewXYZ,TOUR_CB);
  SPINNER_az_path=glui_tour->add_spinner_to_panel(panel_movedir,"Azimuth:",GLUI_SPINNER_FLOAT,&tour_az_path,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_elev_path=glui_tour->add_spinner_to_panel(panel_movedir,"Elevation:",GLUI_SPINNER_FLOAT,&tour_elev_path,KEYFRAME_tXYZ,TOUR_CB);
  SPINNER_elev_path->set_float_limits(-90.0,90.0);
  SPINNER_tourzoom=glui_tour->add_spinner_to_panel(panel_movedir,"Zoom:",GLUI_SPINNER_FLOAT,&tour_zoom,KEYFRAME_tXYZ,TOUR_CB);
  
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
  glui_avatar_index=ti->glui_avatar_index;
  TOUR_CB(TOUR_AVATAR);
  LISTBOX_avatar->set_int_val(glui_avatar_index);
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
  tour_continuity=selected_frame->continuity;
  tour_bias=selected_frame->bias;
  viewtype=selected_frame->viewtype;
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
  SPINNER_tourzoom->set_float_val(tour_zoom);
  SPINNER_viewx->set_float_val(tour_viewx);
  SPINNER_viewy->set_float_val(tour_viewy);
  SPINNER_viewz->set_float_val(tour_viewz);
  SPINNER_az_path->set_float_val(tour_az_path);
  SPINNER_elev_path->set_float_val(tour_elev_path);
  if(ti!=NULL&&CHECKBOXtourhide!=NULL)CHECKBOXtourhide->set_int_val(tour_hide);
  EDITlabel->set_text(tour_label);

  if(edittour==1){
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
  }
  CHECKBOXview1->set_int_val(viewtype);
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
  float dummy;
  tourdata *thistour=NULL;
  float *aview,*eye;

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
    update_tourcontrols();
    TOUR_CB(VIEW1);
    updatemenu=0;
    break;
  case VIEWSNAP:
    if(viewtourfrompath==1&&keyframe_snap==1){
      viewtourfrompath=0;
      view_checkbox->set_int_val(viewtourfrompath);
    }
    TOUR_CB(VIEWTOURFROMPATH);
    updatemenu=0;
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
    if(selected_frame!=NULL){
      selected_frame->viewtype=viewtype;
    }
    createtourpaths();
    break;
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

      adjustviewangle(selected_frame,&tour_az_path,&tour_elev_path);
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
      TOURMENU(-13);  // reset tour vis to ini values
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
