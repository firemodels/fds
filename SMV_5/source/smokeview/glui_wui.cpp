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
GLUI *glui_wui=NULL;
GLUI_Panel *panel_wui=NULL;

/* ------------------ glui_wui_setup ------------------------ */

extern "C" void glui_wui_setup(int main_window){

  if(glui_wui!=NULL)glui_wui->close();
  glui_wui = GLUI_Master.create_glui("WUI",0,0,0);
  if(showwui==0)glui_wui->hide();

  panel_wui = glui_wui->add_panel("",GLUI_PANEL_NONE);

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

/* ------------------ CLIP_CB ------------------------ */

void WUI_CB(int var){
  int i;

  switch (var){
  default:
    ASSERT(0);
    break;
  }
}
