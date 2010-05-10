// $Date: 2008-01-21 21:39:14 -0500 (Mon, 21 Jan 2008) $ 
// $Revision: 1222 $
// $Author: gforney $

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
extern "C" char glui_wui_revision[]="$Revision: 1222 $";
GLUI_Panel *panel_terrain=NULL;

#define TERRAIN_COLORS 35
#define TERRAIN_VERT 34
#define WUI_CLOSE 99
#define SAVE_SETTINGS 98
#define TERRAIN_TYPE 36
#define TERRAIN_MIN 37
#define TERRAIN_MAX 38
#define TERRAIN_FIRE_LINE_UPDATE 39

#define SETVALMIN 1
#define SETVALMAX 2
#define VALMIN 3
#define VALMAX 4
#define FILETYPEINDEX 5
#define CHOPVALMIN 13
#define CHOPVALMAX 14
#define SETCHOPMINVAL 15
#define SETCHOPMAXVAL 16
#define CHOPUPDATE 17
#define FILEUPDATE 6

GLUI *glui_wui=NULL;
GLUI_Panel *panel_wui=NULL;
GLUI_Panel *panel_fire_line=NULL;
GLUI_Panel *panel_terrain_hidden1=NULL;
GLUI_Panel *panel_terrain_color=NULL;
GLUI_Panel *panel_terrain_type=NULL;
GLUI_RadioGroup *RADIO_terrain_type=NULL;
GLUI_RadioButton *RADIO_texture=NULL;
GLUI_Spinner *SPINNER_red_min=NULL;
GLUI_Spinner *SPINNER_green_min=NULL;
GLUI_Spinner *SPINNER_blue_min=NULL;
GLUI_Spinner *SPINNER_red_max=NULL;
GLUI_Spinner *SPINNER_green_max=NULL;
GLUI_Spinner *SPINNER_blue_max=NULL;
GLUI_Spinner *SPINNER_vertical_factor=NULL;
GLUI_Spinner *SPINNER_fire_line_min=NULL;
GLUI_Spinner *SPINNER_fire_line_max=NULL;

extern "C" void WUI_CB(int var);

/* ------------------ update_glui_wui ------------------------ */

extern "C" void update_glui_wui(void){
  if(RADIO_terrain_type!=NULL){
    RADIO_terrain_type->set_int_val(visTerrainType);
  }
}

/* ------------------ glui_wui_setup ------------------------ */

extern "C" void glui_wui_setup(int main_window){

  if(glui_wui!=NULL)glui_wui->close();
  glui_wui = GLUI_Master.create_glui("Terrain",0,0,0);
  if(showwui==0)glui_wui->hide();

  panel_terrain = glui_wui->add_panel("",GLUI_PANEL_NONE);

  if(nterraininfo>0){
    panel_terrain_color=glui_wui->add_panel_to_panel(panel_terrain,"3D surface color");
//  glui_labels->add_statictext_to_panel(panel_tick2,"                    x");
//  glui_labels->add_column_to_panel(panel_tick2,false);

    glui_wui->add_statictext_to_panel(panel_terrain_color,"                  min");
    SPINNER_red_min=glui_wui->add_spinner_to_panel(panel_terrain_color,"red",GLUI_SPINNER_INT,
      terrain_rgba_zmin,TERRAIN_COLORS,WUI_CB);
    SPINNER_green_min=glui_wui->add_spinner_to_panel(panel_terrain_color,"green",GLUI_SPINNER_INT,
      terrain_rgba_zmin+1,TERRAIN_COLORS,WUI_CB);
    SPINNER_blue_min=glui_wui->add_spinner_to_panel(panel_terrain_color,"blue",GLUI_SPINNER_INT,
      terrain_rgba_zmin+2,TERRAIN_COLORS,WUI_CB);
    glui_wui->add_column_to_panel(panel_terrain_color,false);

    glui_wui->add_statictext_to_panel(panel_terrain_color,"                  max");
    SPINNER_red_max=glui_wui->add_spinner_to_panel(panel_terrain_color,"",GLUI_SPINNER_INT,
      terrain_rgba_zmax,TERRAIN_COLORS,WUI_CB);
    SPINNER_green_max=glui_wui->add_spinner_to_panel(panel_terrain_color,"",GLUI_SPINNER_INT,
      terrain_rgba_zmax+1,TERRAIN_COLORS,WUI_CB);
    SPINNER_blue_max=glui_wui->add_spinner_to_panel(panel_terrain_color,"",GLUI_SPINNER_INT,
      terrain_rgba_zmax+2,TERRAIN_COLORS,WUI_CB);

      SPINNER_red_min->set_int_limits(0,255,GLUI_LIMIT_CLAMP);
    SPINNER_green_min->set_int_limits(0,255,GLUI_LIMIT_CLAMP);
     SPINNER_blue_min->set_int_limits(0,255,GLUI_LIMIT_CLAMP);
      SPINNER_red_max->set_int_limits(0,255,GLUI_LIMIT_CLAMP);
    SPINNER_green_max->set_int_limits(0,255,GLUI_LIMIT_CLAMP);
     SPINNER_blue_max->set_int_limits(0,255,GLUI_LIMIT_CLAMP);

    panel_terrain_hidden1=glui_wui->add_panel_to_panel(panel_terrain,"",GLUI_PANEL_NONE);
    panel_terrain_type=glui_wui->add_panel_to_panel(panel_terrain_hidden1,"Surface view");
    RADIO_terrain_type=glui_wui->add_radiogroup_to_panel(panel_terrain_type,&visTerrainType,
      TERRAIN_TYPE,WUI_CB);
    glui_wui->add_radiobutton_to_group(RADIO_terrain_type,"3D surface");
    glui_wui->add_radiobutton_to_group(RADIO_terrain_type,"2D stepped");
    glui_wui->add_radiobutton_to_group(RADIO_terrain_type,"2D lines");
    RADIO_texture=glui_wui->add_radiobutton_to_group(RADIO_terrain_type,"Image");
    glui_wui->add_radiobutton_to_group(RADIO_terrain_type,"Hidden");

    if(terrain_texture==NULL||terrain_texture->loaded==0){
      RADIO_texture->disable();
    }

    glui_wui->add_column_to_panel(panel_terrain_hidden1,false);

    panel_fire_line=glui_wui->add_panel_to_panel(panel_terrain_hidden1,"Fire line");
    SPINNER_fire_line_min=glui_wui->add_spinner_to_panel(panel_fire_line,
      "chop below:",GLUI_SPINNER_FLOAT,&fire_line_min,TERRAIN_MIN,WUI_CB);
    SPINNER_fire_line_max=glui_wui->add_spinner_to_panel(panel_fire_line,
      "color red above:",GLUI_SPINNER_FLOAT,&fire_line_max,TERRAIN_MAX,WUI_CB);
    glui_wui->add_button_to_panel(panel_fire_line,"Update",TERRAIN_FIRE_LINE_UPDATE,WUI_CB);

    SPINNER_vertical_factor=glui_wui->add_spinner_to_panel(panel_terrain_hidden1,
      "vertical exaggeration",GLUI_SPINNER_FLOAT,&vertical_factor,TERRAIN_VERT,WUI_CB);
     SPINNER_vertical_factor->set_float_limits(0.25,4.0,GLUI_LIMIT_CLAMP);

    glui_wui->add_button("Save Settings",SAVE_SETTINGS,WUI_CB);

    glui_wui->add_button("Close",WUI_CLOSE,WUI_CB);

  }

  glui_wui->set_main_gfx_window( main_window );
}

/* ------------------ hide_glui_wui ------------------------ */

extern "C" void hide_glui_wui(void){
  if(glui_wui!=NULL)glui_wui->hide();
  showwui=0;
  updatemenu=1;
}

/* ------------------ show_glui_wui ------------------------ */

extern "C" void show_glui_wui(void){
  if(glui_wui!=NULL)glui_wui->show();
}

/* ------------------ WUI_CB ------------------------ */

extern "C" void WUI_CB(int var){
  int i;
  int fire_line_type;

  switch (var){
    case TERRAIN_MIN:
    case TERRAIN_MAX:
      break;
    case TERRAIN_FIRE_LINE_UPDATE:
      fire_line_type=getslicetype_fromlabel("Fire line");
      if(fire_line_type<0)break;
      list_slice_index=fire_line_type;
      updateslicelist(list_slice_index);


      Slice_CB(FILETYPEINDEX);

      setslicemin=1;
      setslicemax=1;
      slicemin=20.0;
      slicemax=fire_line_max;
      setslicechopmin=1;
      setslicechopmax=0;
      slicechopmin=fire_line_min;

      Slice_CB(SETVALMIN);
      Slice_CB(SETVALMAX);
      Slice_CB(VALMIN);
      Slice_CB(VALMAX);

      Slice_CB(SETCHOPMINVAL);
      Slice_CB(SETCHOPMAXVAL);
      Slice_CB(CHOPVALMIN);
      Slice_CB(CHOPVALMAX);

      Slice_CB(FILEUPDATE);
      Slice_CB(CHOPUPDATE);
      break;
    case TERRAIN_COLORS:
      update_terrain_colors();
      break;
    case TERRAIN_VERT:
      update_terrain(0,vertical_factor);
      break;
    case TERRAIN_TYPE:
      if(visTerrainType==0){
        planar_terrain_slice=0;
      }
      else{
        planar_terrain_slice=1;
      }
      updatemenu=1;
      break;
  case SAVE_SETTINGS:
    writeini(LOCAL_INI);
    break;
  case WUI_CLOSE:
    hide_glui_wui();
    break;

  default:
    ASSERT(0);
    break;
  }
}
