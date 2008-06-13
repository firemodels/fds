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
GLUI_Panel *panel_terrain_color=NULL;

#define TERRAIN_COLORS 35

GLUI *glui_wui=NULL;
GLUI_Panel *panel_wui=NULL;
GLUI_Spinner *SPINNER_red_min=NULL;
GLUI_Spinner *SPINNER_green_min=NULL;
GLUI_Spinner *SPINNER_blue_min=NULL;
GLUI_Spinner *SPINNER_red_max=NULL;
GLUI_Spinner *SPINNER_green_max=NULL;
GLUI_Spinner *SPINNER_blue_max=NULL;
GLUI_Spinner *SPINNER_vertical_factor=NULL;

void WUI_CB(int var);
/* ------------------ glui_wui_setup ------------------------ */

extern "C" void glui_wui_setup(int main_window){

  if(glui_wui!=NULL)glui_wui->close();
  glui_wui = GLUI_Master.create_glui("WUI Display Properties",0,0,0);
  if(showwui==0)glui_wui->hide();

  panel_wui = glui_wui->add_panel("",GLUI_PANEL_NONE);

  if(nterraininfo>0){
    panel_terrain_color = glui_wui->add_panel("Terrain Colors");
    SPINNER_red_min=glui_wui->add_spinner_to_panel(panel_terrain_color,"min red",GLUI_SPINNER_FLOAT,
      terrain_rgba_zmin,TERRAIN_COLORS,WUI_CB);
    SPINNER_green_min=glui_wui->add_spinner_to_panel(panel_terrain_color,"min green",GLUI_SPINNER_FLOAT,
      terrain_rgba_zmin+1,TERRAIN_COLORS,WUI_CB);
    SPINNER_blue_min=glui_wui->add_spinner_to_panel(panel_terrain_color,"min blue",GLUI_SPINNER_FLOAT,
      terrain_rgba_zmin+2,TERRAIN_COLORS,WUI_CB);
    SPINNER_red_max=glui_wui->add_spinner_to_panel(panel_terrain_color,"max red",GLUI_SPINNER_FLOAT,
      terrain_rgba_zmax,TERRAIN_COLORS,WUI_CB);
    SPINNER_green_max=glui_wui->add_spinner_to_panel(panel_terrain_color,"max green",GLUI_SPINNER_FLOAT,
      terrain_rgba_zmax+1,TERRAIN_COLORS,WUI_CB);
    SPINNER_blue_max=glui_wui->add_spinner_to_panel(panel_terrain_color,"max blue",GLUI_SPINNER_FLOAT,
      terrain_rgba_zmax+2,TERRAIN_COLORS,WUI_CB);
      SPINNER_red_min->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);
    SPINNER_green_min->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);
     SPINNER_blue_min->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);
      SPINNER_red_max->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);
    SPINNER_green_max->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);
     SPINNER_blue_max->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);

    SPINNER_vertical_factor=glui_wui->add_spinner_to_panel(panel_terrain_color,"vertical factor",GLUI_SPINNER_FLOAT,
      &vertical_factor,TERRAIN_COLORS,WUI_CB);
     SPINNER_vertical_factor->set_float_limits(0.25,4.0,GLUI_LIMIT_CLAMP);
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

void WUI_CB(int var){
  int i;

  switch (var){
    case TERRAIN_COLORS:
      update_terrain(0,vertical_factor);
      break;
  default:
    ASSERT(0);
    break;
  }
}
