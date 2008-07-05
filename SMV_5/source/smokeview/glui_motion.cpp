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
extern "C" char glui_motion_revision[]="$Revision$";

#define TRANSLATE_XY 101
#define ROTATE_ZX 102
#define GLUI_Z 2
#define MESH_LIST 4
#define EYE_ROTATE 5
#define EYE_ROTATE_90 6
#define EYELEVEL 7
#define FLOORLEVEL 8
#define SPEED 9
#define EYEVIEW_LEVEL 10
#define EYE_X 21
#define EYE_Y 22
#define EYE_Z 23
#define CRAWL 0
#define WALK 1
#define RUN 2
#define PROJECTION 24

#define LABEL_VIEW 4

#define LIST_VIEW 5
#define ADD_LIST_VIEW 6
#define DELETE_LIST_VIEW 7
#define RESTORE_LIST_VIEW 8
#define SAVE_LIST_VIEW 9
#define STARTUP 10
#define CYCLEVIEWS 11
#define ZOOM 12
#define APERTURE 15
#define CURSOR 13
#define SAVE_SETTINGS 14
#define WINDOW_RESIZE 16
#define WINDOWSIZE_LIST 17
#define SNAPVIEW 21

void EYEVIEW_CB(int var);
void BUTTON_hide2_CB(int var);
void BUTTON_Reset_CB(int var);
void TRANSLATE_CB(int var);

GLUI_Listbox *meshlist1=NULL;
GLUI *glui_motion=NULL;
GLUI_Panel *panel_rotatebuttons=NULL, *panel_translate=NULL,*panel_close=NULL;
GLUI_Panel *panel_rotate=NULL;
GLUI_Panel *panel_speed=NULL;
GLUI_Panel *panel_height=NULL;
GLUI_Rollout *panel_motion=NULL;
GLUI_Panel *panel_translate2=NULL,*panel_translate3=NULL;
GLUI_Rollout *panel_projection=NULL;
GLUI_Panel *panel_anglebuttons=NULL;
GLUI_RadioGroup *projection_radio=NULL,*eyeview_radio=NULL;
GLUI_RadioGroup *eyelevel_radio=NULL;
GLUI_Translation *rotate_zx=NULL,*eyerotate_z=NULL;
GLUI_Translation *translate_z=NULL,*translate_xy=NULL;
GLUI_Checkbox *blockpath_checkbox=NULL,*cursor_checkbox=NULL;
GLUI_Button *eyerotate90_z=NULL,*eyelevel=NULL, *floorlevel=NULL, *reset_saved_view=NULL;
GLUI_Button *restore_view=NULL,*save_view=NULL,*add_view=NULL,*delete_view=NULL;
GLUI_Button *startup_button=NULL,*cycle_views_button=NULL;
GLUI_Rollout *reset_panel=NULL;
GLUI_EditText *edit_view_label=NULL;
GLUI_Listbox *view_lists=NULL;
GLUI_Listbox *LIST_windowsize=NULL;
GLUI_Spinner *SPINNER_zoom=NULL,*SPINNER_aperture=NULL;
GLUI_Spinner *SPINNER_speed_crawl=NULL, *SPINNER_speed_walk=NULL;
GLUI_Spinner *SPINNER_xx=NULL, *SPINNER_yy=NULL, *SPINNER_zz=NULL;
GLUI_Spinner *SPINNER_window_width=NULL, *SPINNER_window_height=NULL;
GLUI_Button *window_update=NULL ;
GLUI_Button *button_snap=NULL;

void enable_disable_views(void);

extern "C" void gluiIdle(void){
  GLUI_Master.set_glutIdleFunc(Idle);
}
extern "C" void gluiIdleNULL(void){
  GLUI_Master.set_glutIdleFunc(NULL);
}


/* ------------------ reset_glui_view ------------------------ */

extern "C" void reset_glui_view(int ival){
  view_lists->set_int_val(ival);
  selected_view=ival;
  BUTTON_Reset_CB(RESTORE_LIST_VIEW);
  enable_disable_views();
}

/* ------------------ glui_motion_setup ------------------------ */

extern "C" void enable_reset_saved_view(void){
  if(reset_saved_view!=NULL)reset_saved_view->enable();
}

/* ------------------ update_glui_speed ------------------------ */


extern "C" void update_glui_speed(void){
  TRANSLATE_CB(EYEVIEW_LEVEL);
  if(eyelevel_radio!=NULL)eyelevel_radio->set_int_val(eyeview_level);
}

/* ------------------ update_glui_zoom ------------------------ */

extern "C" void update_glui_zoom(void){
  if(SPINNER_zoom!=NULL)SPINNER_zoom->set_float_val(zoom);
  aperture_glui=zoom2aperture(zoom);
  if(SPINNER_aperture!=NULL)SPINNER_aperture->set_float_val(aperture_glui);
}

/* ------------------ update_camera_label ------------------------ */

extern "C" void update_camera_label(void){
  edit_view_label->set_text(camera_label);
}

/* ------------------ update_cursor_checkbox ------------------------ */

extern "C" void update_cursor_checkbox(void){
  cursor_checkbox->set_int_val(cursorPlot3D);
}

/* ------------------ update_view_list ------------------------ */

extern "C" void update_view_gluilist(void){
  camera *ca;
  for(ca=camera_list_first.next;ca->next!=NULL;ca=ca->next){
    view_lists->add_item(ca->view_id,ca->name);
  }
  view_lists->set_int_val(startup_view_ini);
  selected_view=startup_view_ini;
  enable_disable_views();
  BUTTON_Reset_CB(RESTORE_LIST_VIEW);

}

/* ------------------ glui_motion_setup ------------------------ */

extern "C" void glui_motion_setup(int main_window){
  int i;
#define TRANSLATE_SPEED 0.005
  mesh *meshi;
  int *rotation_index;

  float *eye_xyz;
  if(camera_label!=NULL){
    free(camera_label);
    camera_label=NULL;
  }
  camera_label=(char *)malloc(sizeof(GLUI_String));
  strcpy(camera_label,"current");

  eye_xyz=camera_current->eye;

  if(glui_motion!=NULL)glui_motion->close();
  glui_motion = GLUI_Master.create_glui("Motion/View",0,0,0);
  if(showmotion==0)glui_motion->hide();

  panel_motion = glui_motion->add_rollout("Motion");

  panel_translate2 = glui_motion->add_panel_to_panel(panel_motion,"",GLUI_PANEL_NONE);
  d_eye_xyz[0]=0.0;
  d_eye_xyz[1]=0.0;
  d_eye_xyz[2]=0.0;
  dsave_eye_xyz[0]=0.0;
  dsave_eye_xyz[1]=0.0;
  dsave_eye_xyz[2]=0.0;

  translate_xy=glui_motion->add_translation_to_panel(panel_translate2,"Horizontal",GLUI_TRANSLATION_XY,d_eye_xyz,TRANSLATE_XY,TRANSLATE_CB);
  translate_xy->set_speed(TRANSLATE_SPEED);

  glui_motion->add_column_to_panel(panel_translate2,false);

  translate_z=glui_motion->add_translation_to_panel(panel_translate2,"Vertical",GLUI_TRANSLATION_Y,eye_xyz+2,GLUI_Z,TRANSLATE_CB);
  translate_z->set_speed(TRANSLATE_SPEED);

//  panel_rotate = glui_motion->add_panel_to_panel(panel_motion,"Rotate");
  panel_rotatebuttons = glui_motion->add_panel_to_panel(panel_motion,"",GLUI_PANEL_NONE);

  rotate_zx=glui_motion->add_translation_to_panel(panel_rotatebuttons,"Rotate",GLUI_TRANSLATION_XY,motion_ab,ROTATE_ZX,TRANSLATE_CB);
  glui_motion->add_column_to_panel(panel_rotatebuttons,false);

  eyerotate_z=glui_motion->add_translation_to_panel(panel_rotatebuttons,"first person",GLUI_TRANSLATION_X,
    motion_dir,EYE_ROTATE,TRANSLATE_CB);
  eyerotate_z->set_speed(180.0/(float)screenWidth);
  eyerotate_z->disable();
 
  eyeview_radio=glui_motion->add_radiogroup_to_panel(panel_motion,&eyeview,0,EYEVIEW_CB);
  glui_motion->add_radiobutton_to_group(eyeview_radio,"general rotations");
  glui_motion->add_radiobutton_to_group(eyeview_radio,"first person movement");
  glui_motion->add_radiobutton_to_group(eyeview_radio,"level rotations");

  rotation_index=&camera_current->rotation_index;
  *rotation_index=nmeshes;
  rotation_index_OLD=nmeshes;
  if(nmeshes>1){
    meshlist1 = glui_motion->add_listbox_to_panel(panel_motion,"Rotate about:",rotation_index,MESH_LIST,TRANSLATE_CB);
    for(i=0;i<nmeshes;i++){
      meshi = meshinfo + i;
      meshlist1->add_item(i,meshi->label);
    }
    meshlist1->add_item(nmeshes,"world center");
    meshlist1->set_int_val(*rotation_index);
  }
  panel_anglebuttons = glui_motion->add_panel_to_panel(panel_motion,"",GLUI_PANEL_NONE);
  eyerotate90_z=glui_motion->add_button_to_panel(panel_anglebuttons,"90 deg",EYE_ROTATE_90,TRANSLATE_CB);
  eyerotate90_z->disable();
  eyerotate90_z->set_alignment(GLUI_ALIGN_LEFT);
//  glui_motion->add_column_to_panel(panel_anglebuttons,false);
  button_snap=glui_motion->add_button_to_panel(panel_anglebuttons,"Snap",SNAPVIEW,TRANSLATE_CB);

  //glui_motion->add_column(false);

  panel_projection = glui_motion->add_rollout("View",false);
  projection_radio=glui_motion->add_radiogroup_to_panel(panel_projection,&projection_type,PROJECTION,TRANSLATE_CB);
  glui_motion->add_radiobutton_to_group(projection_radio,"Perspective");
  glui_motion->add_radiobutton_to_group(projection_radio,"Size Preserving");
  SPINNER_zoom=glui_motion->add_spinner_to_panel(panel_projection,"zoom",GLUI_SPINNER_FLOAT,&zoom,
    ZOOM,TRANSLATE_CB);
  SPINNER_zoom->set_float_limits(0.10,10.0,GLUI_LIMIT_CLAMP);
  aperture_glui=zoom2aperture(zoom);
  SPINNER_aperture=glui_motion->add_spinner_to_panel(panel_projection,"aperture",GLUI_SPINNER_FLOAT,&aperture_glui,
    APERTURE,TRANSLATE_CB);
  glui_motion->add_separator_to_panel(panel_projection);
  glui_motion->add_statictext_to_panel(panel_projection,"Window Size");
  
  LIST_windowsize = glui_motion->add_listbox_to_panel(panel_projection,"",&windowsize_pointer,WINDOWSIZE_LIST,TRANSLATE_CB);
  LIST_windowsize->add_item(0,"Custom");
  LIST_windowsize->add_item(1,"-");
  LIST_windowsize->add_item(2, "320x240");
  LIST_windowsize->add_item(3, "640x480");
  LIST_windowsize->add_item(7, "720x480");
  if(max_screenWidth>=800&&max_screenHeight>=480)LIST_windowsize->add_item(4, "800x640");
  if(max_screenWidth>=1024&&max_screenHeight>=768)  LIST_windowsize->add_item(5,"1024x768");
  if(max_screenWidth>=1280&&max_screenHeight>=1024)  LIST_windowsize->add_item(6,"1280x1024");
  update_windowsizelist();

  SPINNER_window_width = glui_motion->add_spinner_to_panel(panel_projection,"width",GLUI_SPINNER_INT,&glui_screenWidth);
  SPINNER_window_width->set_int_limits(100,max_screenWidth);
  SPINNER_window_height = glui_motion->add_spinner_to_panel(panel_projection,"height",GLUI_SPINNER_INT,&glui_screenHeight);
  SPINNER_window_height->set_int_limits(100,max_screenHeight);
  window_update=glui_motion->add_button_to_panel(panel_projection,"Apply",WINDOW_RESIZE,TRANSLATE_CB);

  reset_panel = glui_motion->add_rollout("Save/Restore Views",false);

  restore_view=glui_motion->add_button_to_panel(reset_panel,"Restore",RESTORE_LIST_VIEW,BUTTON_Reset_CB);
  save_view=glui_motion->add_button_to_panel(reset_panel,"Replace",SAVE_LIST_VIEW,BUTTON_Reset_CB);
  add_view=glui_motion->add_button_to_panel(reset_panel,"Add",ADD_LIST_VIEW,BUTTON_Reset_CB);
  delete_view=glui_motion->add_button_to_panel(reset_panel,"Delete",DELETE_LIST_VIEW,BUTTON_Reset_CB);

  glui_motion->add_column_to_panel(reset_panel,false);

  view_lists = glui_motion->add_listbox_to_panel(reset_panel,"Select",&i_view_list,LIST_VIEW,BUTTON_Reset_CB);
  edit_view_label=glui_motion->add_edittext_to_panel(reset_panel,"Edit",GLUI_EDITTEXT_TEXT,camera_label,LABEL_VIEW,BUTTON_Reset_CB);
  startup_button=glui_motion->add_button_to_panel(reset_panel,"view at startup",STARTUP,BUTTON_Reset_CB);
  cycle_views_button=glui_motion->add_button_to_panel(reset_panel,"cycle user views",CYCLEVIEWS,BUTTON_Reset_CB);

  cursor_checkbox=glui_motion->add_checkbox("Map cursor keys for Plot3D use",&cursorPlot3D,CURSOR,TRANSLATE_CB);
  panel_close = glui_motion->add_panel("",GLUI_PANEL_NONE);

  glui_motion->add_button_to_panel(panel_close,"Save Settings",SAVE_SETTINGS,BUTTON_hide2_CB);

  glui_motion->add_column_to_panel(panel_close,false);

  glui_motion->add_button_to_panel(panel_close,"Close",1,BUTTON_hide2_CB);

  showhide_translate(eyeview);
  glui_motion->set_main_gfx_window( main_window );
}

/* ------------------ enable_disable_views ------------------------ */

void enable_disable_views(void){
  int ival;
  ival=view_lists->get_int_val();
  selected_view=ival;
  camera *cex;

  cex=&camera_list_first;
  cex=cex->next;
  cex=cex->next;
  cex=cex->next;
  if(cex->next==NULL){
    cycle_views_button->disable();
  }
  else{
    cycle_views_button->enable();
  }

  switch (ival){
  case 0:
  case 1:
    save_view->disable();
    delete_view->disable();
   // edit_view_label->disable();
    break;
  default:
    edit_view_label->enable();
    if(restore_view!=NULL)restore_view->enable();
    save_view->enable();
    add_view->enable();
    delete_view->enable();

    break;
  }
}

/* ------------------ update_windowsizelist ------------------------ */

extern "C" void update_windowsizelist(void){
  windowsize_pointer=0;
  glui_screenWidth=screenWidth;
  glui_screenHeight=screenHeight;
  if(SPINNER_window_width!=NULL)SPINNER_window_width->set_int_val(screenWidth);
  if(SPINNER_window_height!=NULL)SPINNER_window_height->set_int_val(screenHeight);

  if(screenWidth==320&&screenHeight==240){
    windowsize_pointer=2;
  }
  if(screenWidth==720&&screenHeight==480){
    windowsize_pointer=7;
  }
  if(screenWidth==640&&screenHeight==480){
    windowsize_pointer=3;
  }
  if(screenWidth==800&&screenHeight==640){
    windowsize_pointer=4;
  }
  if(screenWidth==1024&&screenHeight==768){
    windowsize_pointer=5;
  }
  if(screenWidth==1280&&screenHeight==1024){
    windowsize_pointer=6;
  }
  if(LIST_windowsize!=NULL)LIST_windowsize->set_int_val(windowsize_pointer);
}


/* ------------------ update_blockpath ------------------------ */

extern "C" void update_blockpath(void){
  if(blockpath_checkbox!=NULL)blockpath_checkbox->set_int_val(pass_through);
}


/* ------------------ update_translate ------------------------ */

extern "C" void update_translate(void){
  float *eye_xyz,*angle_zx;

  eye_xyz = camera_current->eye;
  angle_zx = camera_current->angle_zx;

  d_eye_xyz[0]=eye_xyz[0]-eye_xyz0[0];
  d_eye_xyz[1]=eye_xyz[1]-eye_xyz0[1];
  d_eye_xyz[2]=eye_xyz[2]-eye_xyz0[2];

  translate_xy->set_x(d_eye_xyz[0]);
  if(eyeview==WORLD_CENTERED_LEVEL){
    d_eye_xyz[1]=0.0;
  }
  translate_xy->set_y(d_eye_xyz[1]);
  translate_z->set_y(eye_xyz[2]);
  rotate_zx->set_x(angle_zx[0]);
  rotate_zx->set_y(angle_zx[1]);
  eyerotate_z->set_x(camera_current->direction_angle);
}

/* ------------------ update_projection_type ------------------------ */

extern "C" void update_projection_type(void){
  if(projection_radio!=NULL)projection_radio->set_int_val(projection_type);
  if(projection_type==1){
    if(SPINNER_zoom!=NULL)    SPINNER_zoom->disable();
    if(SPINNER_aperture!=NULL)SPINNER_aperture->disable();
  }
  else{
    if(SPINNER_zoom!=NULL)    SPINNER_zoom->enable();
    if(SPINNER_aperture!=NULL)SPINNER_aperture->enable();
  }
}

extern "C" void update_eyerotate(void){
  eyerotate_z->set_x(camera_current->direction_angle);
}

/* ------------------ showhide_translate ------------------------ */

extern "C" void showhide_translate(int var){
  float *eye_xyz;

  eye_xyz = camera_current->eye;

  eye_xyz0[0]=eye_xyz[0];
  eye_xyz0[1]=eye_xyz[1];
  eye_xyz0[2]=eye_xyz[2];
  d_eye_xyz[0]=0.0;
  d_eye_xyz[1]=0.0;
  switch (var){
  case WORLD_CENTERED:
    if(panel_translate!=NULL)panel_translate->enable();
    if(rotate_zx!=NULL)rotate_zx->enable();
    if(eyerotate_z!=NULL)eyerotate_z->disable();
    if(eyerotate90_z!=NULL)eyerotate90_z->disable();
    if(blockpath_checkbox!=NULL)blockpath_checkbox->disable();
    if(eyelevel_radio!=NULL)eyelevel_radio->disable();
    if(panel_speed!=NULL)panel_speed->disable();
    if(panel_height!=NULL)panel_height->disable();
    if(eyeview_radio!=NULL)eyeview_radio->set_int_val(eyeview);
    if(eyelevel!=NULL)eyelevel->disable();
    if(floorlevel!=NULL)floorlevel->disable();
    if(meshlist1!=NULL)meshlist1->enable();
    if(button_snap!=NULL)button_snap->enable();
    break;
  case EYE_CENTERED:
    if(panel_translate!=NULL)panel_translate->enable();
    if(rotate_zx!=NULL)rotate_zx->disable();
    if(eyerotate_z!=NULL)eyerotate_z->enable();
    if(eyerotate90_z!=NULL)eyerotate90_z->enable();
    if(blockpath_checkbox!=NULL)blockpath_checkbox->enable();
    if(panel_speed!=NULL)panel_speed->enable();
    if(panel_height!=NULL)panel_height->enable();
    if(eyeview_radio!=NULL)eyeview_radio->set_int_val(eyeview);
    if(eyelevel!=NULL)eyelevel->enable();
    if(floorlevel!=NULL)floorlevel->enable();
    if(meshlist1!=NULL)meshlist1->disable();
    if(button_snap!=NULL)button_snap->disable();
    break;
  case WORLD_CENTERED_LEVEL:
    if(panel_translate!=NULL)panel_translate->enable();
    if(rotate_zx!=NULL)rotate_zx->enable();
    if(eyerotate_z!=NULL)eyerotate_z->disable();
    if(eyerotate90_z!=NULL)eyerotate90_z->disable();
    if(blockpath_checkbox!=NULL)blockpath_checkbox->disable();
    if(eyelevel_radio!=NULL)eyelevel_radio->disable();
    if(panel_speed!=NULL)panel_speed->disable();
    if(panel_height!=NULL)panel_height->disable();
    if(eyeview_radio!=NULL)eyeview_radio->set_int_val(eyeview);
    if(eyelevel!=NULL)eyelevel->disable();
    if(floorlevel!=NULL)floorlevel->disable();
    if(meshlist1!=NULL)meshlist1->enable();
    if(button_snap!=NULL)button_snap->enable();
    break;
  default:
    ASSERT(FFALSE);
  }

}

/* ------------------ TRANSLATE_CB ------------------------ */

void TRANSLATE_CB(int var){
  float dx, dy;
  float dx2, dy2;
  float *eye_xyz;
  float *direction_angle;
  float *cos_direction_angle, *sin_direction_angle;
  int *rotation_index;

  if(var==CURSOR){
    updatemenu=1;
    return;
  }
  eye_xyz = camera_current->eye;
  direction_angle=&camera_current->direction_angle;
  cos_direction_angle=&camera_current->cos_direction_angle;
  sin_direction_angle=&camera_current->sin_direction_angle;
  rotation_index = &camera_current->rotation_index;
  if(selected_view!=-999){
    selected_view=-999;
    updatemenu=1;
  }

  switch (var){

    case EYE_ROTATE:
      *direction_angle=motion_dir[0];
      *cos_direction_angle = cos(PI*(*direction_angle)/180.0);
      *sin_direction_angle = sin(PI*(*direction_angle)/180.0);
      if(glui_move_mode!=EYE_ROTATE){
        eye_xyz0[0]=eye_xyz[0];
        eye_xyz0[1]=eye_xyz[1];
        update_translate();
      }
      glui_move_mode=EYE_ROTATE;
      return;
      break;
    case EYE_ROTATE_90:
      {
        float diffangle;
        int intangle;

        intangle = (int)((*direction_angle+45)/90)*90.0;
        diffangle = *direction_angle-intangle;
        if(diffangle<0.0)diffangle = -diffangle;
        if(diffangle>1.0){
          *direction_angle=intangle;
        }
        else{
          *direction_angle+=90.0;
        }
      }
      if(*direction_angle>=360.0)*direction_angle-=360.0;
      TRANSLATE_CB(EYE_ROTATE);
      glui_move_mode=EYE_ROTATE_90;
      return;
    case ROTATE_ZX:
      {
        float *angle_zx;
        angle_zx = camera_current->angle_zx;
        angle_zx[0] = rotate_zx->get_x();
        angle_zx[1] = -rotate_zx->get_y();
      }
      break;
    case WINDOWSIZE_LIST:
      switch (windowsize_pointer){
        case 2:
          glui_screenWidth=320;
          glui_screenHeight=240;
          break;
        case 3:
          glui_screenWidth=640;
          glui_screenHeight=480;
          break;
        case 7:
          glui_screenWidth=720;
          glui_screenHeight=480;
          break;
        case 4:
          glui_screenWidth=800;
          glui_screenHeight=640;
          break;
        case 5:
          glui_screenWidth=1024;
          glui_screenHeight=768;
          break;
        case 6:
          glui_screenWidth=1280;
          glui_screenHeight=1024;
          break;
      }
      if(windowsize_pointer>=2){
        SPINNER_window_width->set_int_val(glui_screenWidth);
        SPINNER_window_height->set_int_val(glui_screenHeight);
        screenWidth=glui_screenWidth;
        screenHeight=glui_screenHeight;
        ResizeWindow(screenWidth,screenHeight);
      }
      break;
    case SNAPVIEW:
      snap_view_angles();
      break;
    case WINDOW_RESIZE:
      screenWidth=glui_screenWidth;
      screenHeight=glui_screenHeight;
      update_windowsizelist();
      ResizeWindow(screenWidth,screenHeight);
      break;

    case PROJECTION:
      if(projection_type==0){
        projection_type=1;
      }
      else{
        projection_type=0;
      }
      ZoomMenu(-2);
      camera_current->projection_type=projection_type;
      return;
      break;
    case EYELEVEL:
      desired_view_height=1.5;
      break;
    case FLOORLEVEL:
      desired_view_height=0.6;
      break;
    case SPEED:
      switch (eyeview_level){
      case CRAWL:
        setspeed(speed_crawl);
        break;
      case WALK:
        setspeed(speed_walk);
        break;
      }
      break;
    case EYEVIEW_LEVEL:
      switch (eyeview_level){
      case 0:
        desired_view_height=0.6;
        setspeed(speed_crawl);
        break;
      case 1:
        desired_view_height=1.5;
        setspeed(speed_walk);
        break;
      case 2:
        desired_view_height=1.5;
        setspeed(speed_walk);
        break;
      }
      break;
    case MESH_LIST:
      if(*rotation_index>=0&&*rotation_index<nmeshes){
        update_current_mesh(meshinfo + (*rotation_index));
        update_rotation_index(*rotation_index);
      }
      else{
        update_current_mesh(meshinfo);
        update_rotation_index(nmeshes);
      }
      return;
      break;
    case ZOOM:
      zoomindex=-1;
      camera_current->zoom=zoom;
      aperture_glui=zoom2aperture(zoom);
      if(SPINNER_aperture!=NULL)SPINNER_aperture->set_float_val(aperture_glui);
      break;
    case APERTURE:
      zoom=aperture2zoom(aperture_glui);
      if(zoom<0.1||zoom>10.0){
        if(zoom<0.1)zoom=0.1;
        if(zoom>10.0)zoom=10.0;
        aperture_glui=zoom2aperture(zoom);
        if(SPINNER_aperture!=NULL)SPINNER_aperture->set_float_val(aperture_glui);
      }
      zoomindex=-1;
      camera_current->zoom=zoom;
      if(SPINNER_zoom!=NULL)SPINNER_zoom->set_float_val(zoom);
      break;
//    default:
//      ASSERT(FFALSE);
//      break;
  }
  
  dx = d_eye_xyz[0];
  dy = d_eye_xyz[1];
  if(var==EYE_ROTATE){
    dy=motion_dir[1]*TRANSLATE_SPEED*(float)screenWidth/1800.0;
  }
  if(eyeview==EYE_CENTERED){
    *cos_direction_angle = cos(PI*(*direction_angle)/180.0);
    *sin_direction_angle = sin(PI*(*direction_angle)/180.0);
    dx2 = *cos_direction_angle*dx + *sin_direction_angle*dy;
    dy2 = -(*sin_direction_angle)*dx + (*cos_direction_angle)*dy;
    dx = dx2;
    dy = dy2;
  }

  if(glui_move_mode==EYE_ROTATE){
    getnewpos(eye_xyz,dx,dy,0.0,1.0);
    eye_xyz0[0]=eye_xyz[0];
    eye_xyz0[1]=eye_xyz[1];
    eye_xyz0[2]=eye_xyz[2];
  }
  else{
    eye_xyz[0] = eye_xyz0[0] + dx;
    eye_xyz[1] = eye_xyz0[1] + dy;
  }

  switch (var){
    case EYE_ROTATE:
    case TRANSLATE_XY:
      if(glui_move_mode==EYE_ROTATE){
        eye_xyz0[0]=eye_xyz[0];
        eye_xyz0[1]=eye_xyz[1];
        update_translate();
      }
      if(translate_xy!=NULL){
        translate_xy->set_x(d_eye_xyz[0]);
        translate_xy->set_y(d_eye_xyz[1]);
      }
      glui_move_mode=TRANSLATE_XY;
      break;
    case GLUI_Z:
      if(glui_move_mode==EYE_ROTATE){
        eye_xyz0[0]=eye_xyz[0];
        eye_xyz0[1]=eye_xyz[1];
        update_translate();
      }
      glui_move_mode=GLUI_Z;
      break;
    case APERTURE:
    case ZOOM:
    case FLOORLEVEL:
    case EYELEVEL:
    case SPEED:
    case EYEVIEW_LEVEL:
    case PROJECTION:
    case WINDOW_RESIZE:
    case WINDOWSIZE_LIST:
    case SNAPVIEW:
    case ROTATE_ZX:
      break;
    default:
      ASSERT(FFALSE);
      break;
  }


}

/* ------------------ update_meshlist1 ------------------------ */

extern "C" void update_meshlist1(int val){
  if(meshlist1==NULL)return;
  meshlist1->set_int_val(val);
  if(val>=0&&val<nmeshes){
    eyeview_radio->set_int_val(0);
    handle_eyeview(0);
  }
}

/* ------------------ hide_glui_motion ------------------------ */

extern "C" void hide_glui_motion(void){
  if(glui_motion!=NULL)glui_motion->hide();
  showmotion=0;
}

/* ------------------ show_glui_motion_setup ------------------------ */

extern "C" void show_glui_motion(void){
  if(glui_motion!=NULL)glui_motion->show();
}

/* ------------------ BUTTON_hide2_CB ------------------------ */

void BUTTON_hide2_CB(int var){
  switch (var){
  case 1:
    if(glui_motion!=NULL)glui_motion->hide();
    showmotion=0;
    updatemenu=1;
    break;
  case SAVE_SETTINGS:
    updatemenu=1;
    writeini(LOCAL_INI);
    break;
  }
}

/* ------------------ BUTTON_Reset_CB ------------------------ */

void BUTTON_Reset_CB(int var){
  int ival;
  int eyeview_save;
  camera *cam1,*cex,*ca;
  char *label;
  camera *prev, *next;
  int view_id;

  switch (var){
  case RESTORE_EXTERIOR_VIEW:
  case RESTORE_INTERIOR_VIEW:
  case RESTORE_SAVED_VIEW:
    ResetView(var);
    break;
  case SAVE_VIEW:
    strcpy(camera_current->name,camera_label);
    reset_saved_view->enable();
    ViewpointMenu(SAVE_VIEW);
    break;
  case LABEL_VIEW:
    updatemenu=1;
    break;
  case SAVE_LIST_VIEW:
    ival=view_lists->get_int_val();
    selected_view=ival;
    label=edit_view_label->get_text();
    cex=&camera_list_first;
    cex=cex->next;
    cex=cex->next;
    for(ca=cex;ca->next!=NULL;ca=ca->next){
      if(ca->view_id==ival)break;
    }
    if(ival==ca->view_id){
      cam1=ca;
    }
    else{
      return;
    }
    prev=ca->prev;
    next=ca->next;
    view_id=ca->view_id;
    copy_camera(ca,camera_current);
    ca->prev=prev;
    ca->next=next;
    ca->view_id=view_id;

    view_lists->delete_item(ival);
    view_lists->add_item(ival,label);
    strcpy(ca->name,label);

    break;
  case ADD_LIST_VIEW:
    add_list_view(NULL);
    break;
  case DELETE_LIST_VIEW:
    ival=view_lists->get_int_val();
    label=edit_view_label->get_text();
    cex=&camera_list_first;
    cex=cex->next;
    cex=cex->next;
    for(ca=cex;ca->next!=NULL;ca=ca->next){
      if(ca->view_id==ival)break;
    }
    if(ival==ca->view_id){
      cam1=ca;
    }
    else{
      return;
    }
    view_lists->delete_item(ival);
    prev=cam1->prev;
    delete_camera(cam1);
    view_lists->set_int_val(prev->view_id);
    selected_view=prev->view_id;
    BUTTON_Reset_CB(RESTORE_LIST_VIEW);
    enable_disable_views();
    break;
  case RESTORE_LIST_VIEW:
    ival=view_lists->get_int_val();
    selected_view=ival;
    for(ca=camera_list_first.next;ca->next!=NULL;ca=ca->next){
      if(ca->view_id==ival)break;
    }

   eyeview_save = ca->eyeview;
   copy_camera(camera_current,ca);
   camera_current->eyeview=eyeview_save;
   edit_view_label->set_text(ca->name);
   break;
  case LIST_VIEW:
    BUTTON_Reset_CB(RESTORE_LIST_VIEW);
    enable_disable_views();
    break;
  case STARTUP:
    startup_view_ini=view_lists->get_int_val();
    selected_view=startup_view_ini;
    writeini(LOCAL_INI);
    break;
  case CYCLEVIEWS:
    ival=view_lists->get_int_val();
    selected_view=ival;
    cex=&camera_list_first;
    cex=cex->next;
    cex=cex->next;
    switch (ival){
    case 0:
    case 1:
      cex=cex->next;
      if(cex->next==NULL)return;
      ival=cex->view_id;
      break;
    default:
      for(ca=cex;ca->next!=NULL;ca=ca->next){
        if(ca->view_id==ival)break;
      }
      cex=ca->next;
      if(cex->next==NULL){
        cex=&camera_list_first;
        cex=cex->next;
        cex=cex->next;
        cex=cex->next;
        if(cex->next==NULL)return;
        ival=cex->view_id;
      }
      else{
        ival=cex->view_id;
      }
      break;
    }
    view_lists->set_int_val(ival);
    selected_view=ival;
    BUTTON_Reset_CB(RESTORE_LIST_VIEW);
    break;
    default:
      ASSERT(FFALSE);
      break;
  }
}

void set_startup_view(void){
  BUTTON_Reset_CB(STARTUP);
}

void add_list_view(char *label_in){
  int ival;
  char *label;
  camera *cam1,*cam2,*cex,*ca;

  ival=view_lists->get_int_val();
  selected_view=ival;
  label=label_in;
  if(label==NULL)label=edit_view_label->get_text();
  cex=&camera_list_first;
  cex=cex->next;
  cex=cex->next;
  for(ca=cex;ca->next!=NULL;ca=ca->next){
    if(ca->view_id==ival)break;
  }
  if(ival==ca->view_id){
    cam1=ca;
  }
  else{
    cam1=cex;
  }
  cam2 = insert_camera(cam1,camera_current,label);
  if(cam2!=NULL){
    view_lists->add_item(cam2->view_id,cam2->name);
    view_lists->set_int_val(cam2->view_id);
    selected_view=cam2->view_id;
  }
  enable_disable_views();
}

/* ------------------ EYEVIEW_CB ------------------------ */

void EYEVIEW_CB(int var){
  handle_eyeview(0);
}

