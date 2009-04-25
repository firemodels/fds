// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#ifdef pp_SHOOTER
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

GLUI *glui_shooter=NULL;


// svn revision character string
extern "C" char glui_shooter_revision[]="$Revision$";
// $Date$ $Author$

#define SAVE_SETTINGS 900
#define CLOSE_SHOOTER 901

void SHOOTER_CB(int var);

/* ------------------ hide_glui_shooter ------------------------ */

extern "C" void hide_shooter(void){
  if(glui_shooter!=NULL){
    glui_shooter->hide();
    showshooter=0;
    updatemenu=1;
  }
}

/* ------------------ show_shooter ------------------------ */

extern "C" void show_shooter(void){
  if(glui_shooter!=NULL){
    glui_shooter->show();
    showshooter=1;
    updatemenu=1;
  }
}

/* ------------------ glui_shooter_setup ------------------------ */

extern "C" void glui_shooter_setup(int main_window){  

  glui_shooter = GLUI_Master.create_glui( "Particle Shooting",0,0,0 );
  if(showshooter==0)glui_shooter->hide();


  glui_shooter->add_button("Save Settings",SAVE_SETTINGS,SHOOTER_CB);

  glui_shooter->add_button("Close",CLOSE_SHOOTER,SHOOTER_CB);

  glui_shooter->set_main_gfx_window( main_window );
}

/* ------------------ SHOOTER_CB ------------------------ */

void SHOOTER_CB(int var){
}

#endif
