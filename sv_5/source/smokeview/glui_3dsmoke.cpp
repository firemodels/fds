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
#include "MALLOC.h"
#define CPP
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

// svn revision character string
extern "C" char glui_3dsmoke_revision[]="$Revision: 614 $";

extern "C" void update2_glui_smoke3dframestep(void);
extern GLUI_Rollout *panel_smoke3d;
extern GLUI *glui_bounds;

#define IDLE() Idle();

void SMOKE_3D_CB(int var);
void update_alpha(void);

#define SMOKE_3D_CLOSE 0
#define FIRE_RED 1
#define FIRE_GREEN 2
#define FIRE_BLUE 3
#define FIRE_HALFDEPTH 4
#define FIRE_ALPHA 5
#define FIRE_CUTOFF 6
#define SMOKE_SHADE 7
#define SMOKE_THICK 8
#define SAVE_SETTINGS 9
#define FRAMELOADING 10
#define SMOKETEST 11

GLUI *glui_3dsmoke=NULL;
GLUI_RadioGroup *alphagroup=NULL,*skipframes;
GLUI_Checkbox *CHECKBOX_smokecullflag=NULL;
GLUI_Checkbox *CHECKBOX_smokedrawtest=NULL;
GLUI_Checkbox *CHECKBOX_smokedrawtest2=NULL;
GLUI_Checkbox *CHECKBOX_smoke3d_external=NULL;
GLUI_Checkbox *CHECKBOX_zlib=NULL;
GLUI_Spinner *SPINNER_smokedrawtest_nummin=NULL;
GLUI_Spinner *SPINNER_smoke3dframestep=NULL;
GLUI_Spinner *SPINNER_smokedrawtest_nummax=NULL;
GLUI_Spinner *SPINNER_smoke3d_thick=NULL;
GLUI_Spinner *SPINNER_smoke3d_smoke_shade=NULL;
GLUI_Spinner *SPINNER_smoke3d_fire_red=NULL;
GLUI_Spinner *SPINNER_smoke3d_fire_green=NULL;
GLUI_Spinner *SPINNER_smoke3d_fire_blue=NULL;
GLUI_Spinner *SPINNER_smoke3d_hrrpuv_cutoff=NULL;
GLUI_Spinner *SPINNER_smoke3d_fire_halfdepth=NULL;
GLUI_Spinner **SPINNER_smoke3d_hrrpuv_cutoffptr=NULL;
GLUI_Panel *panel_buttons,*panel1=NULL,*panel2=NULL,*panel3=NULL,*panel4=NULL,*panel5=NULL,*panel6=NULL,*panel7=NULL;
GLUI_Spinner *SPINNER_extinct=NULL;
GLUI_Spinner *SPINNER_smokedens=NULL;
GLUI_Spinner *SPINNER_pathlength=NULL;
GLUI_StaticText *TEXT_smokealpha=NULL;
GLUI_StaticText *TEXT_smokedepth=NULL;
GLUI_Checkbox *CHECKBOX_show_smoketest=NULL;

extern "C" void update_smoke3dflags(void){
  alphagroup->set_int_val(adjustalphaflag);
  CHECKBOX_smokecullflag->set_int_val(smokecullflag);
  if(CHECKBOX_smokedrawtest!=NULL)CHECKBOX_smokedrawtest->set_int_val(smokedrawtest);
  if(CHECKBOX_smokedrawtest2!=NULL)CHECKBOX_smokedrawtest2->set_int_val(smokedrawtest2);
  skipframes->set_int_val(smokeskipm1);

}

extern "C" void glui_3dsmoke_setup2(int main_window){

  int i;

  
  if(nsmoke3d<=0)return;
  SPINNER_smoke3d_hrrpuv_cutoffptr=(GLUI_Spinner **)malloc(nmeshes*sizeof(GLUI_Spinner *));
  
  glui_3dsmoke=glui_bounds;
  panel4 = glui_3dsmoke->add_panel_to_panel(panel_smoke3d,"Absorption Adjustments");
  panel4->set_alignment(GLUI_ALIGN_LEFT);
  alphagroup = glui_3dsmoke->add_radiogroup_to_panel(panel4,&adjustalphaflag);
  glui_3dsmoke->add_radiobutton_to_group(alphagroup,"none");
  glui_3dsmoke->add_radiobutton_to_group(alphagroup,"adjust off-center");
  glui_3dsmoke->add_radiobutton_to_group(alphagroup,"adjust off-center & zero at boundaries");
  glui_3dsmoke->add_radiobutton_to_group(alphagroup,"zero at boundaries");

  if(smoketest==1){
    panel7 = glui_3dsmoke->add_panel_to_panel(panel_smoke3d,"Test Smoke");
    panel7->set_alignment(GLUI_ALIGN_LEFT);
    CHECKBOX_show_smoketest=glui_3dsmoke->add_checkbox_to_panel(panel7,"Show Test Smoke",&show_smoketest);
    SPINNER_extinct=glui_3dsmoke->add_spinner_to_panel(panel7,"Mass Extinction Coeff (m2/g)",GLUI_SPINNER_FLOAT,&smoke_extinct,SMOKETEST,SMOKE_3D_CB);
    SPINNER_extinct->set_float_limits(0.0,10.0);
    SPINNER_smokedens=glui_3dsmoke->add_spinner_to_panel(panel7,"Smoke Density (g/m3)",GLUI_SPINNER_FLOAT,&smoke_dens,SMOKETEST,SMOKE_3D_CB);
    SPINNER_smokedens->set_float_limits(0.0,1.0);
    SPINNER_pathlength=glui_3dsmoke->add_spinner_to_panel(panel7,"Path Length (m)",GLUI_SPINNER_FLOAT,&smoke_pathlength,SMOKETEST,SMOKE_3D_CB);
    SPINNER_pathlength->set_float_limits(0.0,20.0);
    TEXT_smokealpha=glui_3dsmoke->add_statictext_to_panel(panel7,"Alpha");
    TEXT_smokedepth=glui_3dsmoke->add_statictext_to_panel(panel7,"Depth");
    update_alpha();
  }




  panel1 = glui_3dsmoke->add_panel_to_panel(panel_smoke3d,"",GLUI_PANEL_NONE);

  panel2 = glui_3dsmoke->add_panel_to_panel(panel1,"Fire");
  panel2->set_alignment(GLUI_ALIGN_LEFT);
  SPINNER_smoke3d_fire_red=glui_3dsmoke->add_spinner_to_panel(panel2,"red",GLUI_SPINNER_INT,&fire_red,FIRE_RED,SMOKE_3D_CB);
  SPINNER_smoke3d_fire_red->set_int_limits(0,255);
  SPINNER_smoke3d_fire_green=glui_3dsmoke->add_spinner_to_panel(panel2,"green",GLUI_SPINNER_INT,&fire_green,FIRE_GREEN,SMOKE_3D_CB);
  SPINNER_smoke3d_fire_green->set_int_limits(0,255);
  SPINNER_smoke3d_fire_blue=glui_3dsmoke->add_spinner_to_panel(panel2,"blue",GLUI_SPINNER_INT,&fire_blue,FIRE_BLUE,SMOKE_3D_CB);
  SPINNER_smoke3d_fire_blue->set_int_limits(0,255);
  SPINNER_smoke3d_fire_halfdepth=glui_3dsmoke->add_spinner_to_panel(panel2,"50% flame depth",GLUI_SPINNER_FLOAT,&fire_halfdepth,FIRE_HALFDEPTH,SMOKE_3D_CB);
  SPINNER_smoke3d_fire_halfdepth->set_float_limits(0.0,10.0);


  {
    mesh *meshi;

    for(i=0;i<selected_case->nmeshes;i++){
      meshi = selected_case->meshinfo + i;
      SPINNER_smoke3d_hrrpuv_cutoffptr[i]=glui_3dsmoke->add_spinner_to_panel
        (panel2,meshi->label,GLUI_SPINNER_FLOAT,&meshi->hrrpuv_cutoff,FIRE_CUTOFF,SMOKE_3D_CB);
      SPINNER_smoke3d_hrrpuv_cutoffptr[i]->set_float_limits(0.0,1200);
    }
  }

 glui_3dsmoke->add_column_to_panel(panel1,false);

  panel6 = glui_3dsmoke->add_panel_to_panel(panel1,"Smoke");
  panel6->set_alignment(GLUI_ALIGN_LEFT);
  SPINNER_smoke3d_smoke_shade=glui_3dsmoke->add_spinner_to_panel(panel6,"Gray Level",GLUI_SPINNER_INT,&smoke_shade,SMOKE_SHADE,SMOKE_3D_CB);
  SPINNER_smoke3d_smoke_shade->set_int_limits(0,255);
  SPINNER_smoke3d_thick=glui_3dsmoke->add_spinner_to_panel(panel6,"Thickness",GLUI_SPINNER_INT,&smoke3d_thick,SMOKE_THICK,SMOKE_3D_CB);
  SPINNER_smoke3d_thick->set_int_limits(0,7);
#ifdef _DEBUG
  CHECKBOX_smoke3d_external=glui_3dsmoke->add_checkbox_to_panel(panel6,"View smoke externally (ONLY)",&smoke3d_external);
#endif



  panel5 = glui_3dsmoke->add_panel_to_panel(panel1,"Slices");
  panel5->set_alignment(GLUI_ALIGN_LEFT);
  CHECKBOX_smokecullflag=glui_3dsmoke->add_checkbox_to_panel(panel5,"Cull hidden slices",&smokecullflag);
#ifdef _DEBUG
  CHECKBOX_smokedrawtest=glui_3dsmoke->add_checkbox_to_panel(panel5,"Show Only Back Slices",&smokedrawtest);
  CHECKBOX_smokedrawtest2=glui_3dsmoke->add_checkbox_to_panel(panel5,"Show Only Dir X Slices",&smokedrawtest2);


  SPINNER_smokedrawtest_nummin=glui_3dsmoke->add_spinner_to_panel(panel5,"Back Slice",GLUI_SPINNER_INT,&smokedrawtest_nummin);
  SPINNER_smokedrawtest_nummin->set_int_limits(1,ijkbarmax);

  SPINNER_smokedrawtest_nummax=glui_3dsmoke->add_spinner_to_panel(panel5,"Front Slice",GLUI_SPINNER_INT,&smokedrawtest_nummax);
  SPINNER_smokedrawtest_nummax->set_int_limits(1,ijkbarmax);
#endif
  skipframes = glui_3dsmoke->add_radiogroup_to_panel(panel5,&smokeskipm1);
  glui_3dsmoke->add_radiobutton_to_group(skipframes,"display all");
  glui_3dsmoke->add_radiobutton_to_group(skipframes,"   ... every 2nd");
  glui_3dsmoke->add_radiobutton_to_group(skipframes,"   ... every 3rd");




//  panel_buttons = glui_3dsmoke->add_panel("",false);

//  SPINNER_smoke3dframestep=glui_3dsmoke->add_spinner_to_panel(panel_buttons,"Frame Skip",GLUI_SPINNER_INT,&smoke3dframeskip,
//    FRAMELOADING,SMOKE_3D_CB);
//  SPINNER_smoke3dframestep->set_int_limits(0,100);

//  glui_3dsmoke->add_column_to_panel(panel_buttons,false);

//  glui_3dsmoke->add_button_to_panel(panel_buttons,"Save Settings",SAVE_SETTINGS,SMOKE_3D_CB);

//  glui_3dsmoke->add_column_to_panel(panel_buttons,false);

//  glui_3dsmoke->add_button_to_panel(panel_buttons,"Close",SMOKE_3D_CLOSE,SMOKE_3D_CB);
  
//  glui_3dsmoke->set_main_gfx_window( main_window );
}

/* ------------------ glui_3dsmoke_setup ------------------------ */

extern "C" void glui_3dsmoke_setup(int main_window){

  int i;

  SPINNER_smoke3d_hrrpuv_cutoffptr=(GLUI_Spinner **)malloc(nmeshes*sizeof(GLUI_Spinner *));



  if(glui_3dsmoke!=NULL)glui_3dsmoke->close();
  glui_3dsmoke = GLUI_Master.create_glui("3D Smoke",0,0,0);
  if(showglui3dsmoke==0)glui_3dsmoke->hide();

  panel4 = glui_3dsmoke->add_panel("Absorption Adjustments");
  panel4->set_alignment(GLUI_ALIGN_LEFT);
  alphagroup = glui_3dsmoke->add_radiogroup_to_panel(panel4,&adjustalphaflag);
  glui_3dsmoke->add_radiobutton_to_group(alphagroup,"none");
  glui_3dsmoke->add_radiobutton_to_group(alphagroup,"adjust off-center");
  glui_3dsmoke->add_radiobutton_to_group(alphagroup,"adjust off-center & zero at boundaries");

  panel1 = glui_3dsmoke->add_panel("",GLUI_PANEL_NONE);

  panel2 = glui_3dsmoke->add_panel_to_panel(panel1,"Fire");
  panel2->set_alignment(GLUI_ALIGN_LEFT);
  SPINNER_smoke3d_fire_red=glui_3dsmoke->add_spinner_to_panel(panel2,"red",GLUI_SPINNER_INT,&fire_red,FIRE_RED,SMOKE_3D_CB);
  SPINNER_smoke3d_fire_red->set_int_limits(0,255);
  SPINNER_smoke3d_fire_green=glui_3dsmoke->add_spinner_to_panel(panel2,"green",GLUI_SPINNER_INT,&fire_green,FIRE_GREEN,SMOKE_3D_CB);
  SPINNER_smoke3d_fire_green->set_int_limits(0,255);
  SPINNER_smoke3d_fire_blue=glui_3dsmoke->add_spinner_to_panel(panel2,"blue",GLUI_SPINNER_INT,&fire_blue,FIRE_BLUE,SMOKE_3D_CB);
  SPINNER_smoke3d_fire_blue->set_int_limits(0,255);
  SPINNER_smoke3d_fire_halfdepth=glui_3dsmoke->add_spinner_to_panel(panel2,"50% flame depth",GLUI_SPINNER_FLOAT,&fire_halfdepth,FIRE_HALFDEPTH,SMOKE_3D_CB);
  SPINNER_smoke3d_fire_halfdepth->set_float_limits(0.0,10.0);


  {
    mesh *meshi;

    for(i=0;i<selected_case->nmeshes;i++){
      meshi = selected_case->meshinfo + i;
      SPINNER_smoke3d_hrrpuv_cutoffptr[i]=glui_3dsmoke->add_spinner_to_panel
        (panel2,meshi->label,GLUI_SPINNER_FLOAT,&meshi->hrrpuv_cutoff,FIRE_CUTOFF,SMOKE_3D_CB);
      SPINNER_smoke3d_hrrpuv_cutoffptr[i]->set_float_limits(0.0,1200);
    }
  }

 glui_3dsmoke->add_column_to_panel(panel1,false);

  panel6 = glui_3dsmoke->add_panel_to_panel(panel1,"Smoke");
  panel6->set_alignment(GLUI_ALIGN_LEFT);
  SPINNER_smoke3d_smoke_shade=glui_3dsmoke->add_spinner_to_panel(panel6,"Gray Level",GLUI_SPINNER_INT,&smoke_shade,SMOKE_SHADE,SMOKE_3D_CB);
  SPINNER_smoke3d_smoke_shade->set_int_limits(0,255);
  SPINNER_smoke3d_thick=glui_3dsmoke->add_spinner_to_panel(panel6,"Thickness",GLUI_SPINNER_INT,&smoke3d_thick,SMOKE_THICK,SMOKE_3D_CB);
  SPINNER_smoke3d_thick->set_int_limits(0,7);
#ifdef _DEBUG
  CHECKBOX_smoke3d_external=glui_3dsmoke->add_checkbox_to_panel(panel6,"View smoke externally (ONLY)",&smoke3d_external);
#endif



  panel5 = glui_3dsmoke->add_panel_to_panel(panel1,"Slices");
  panel5->set_alignment(GLUI_ALIGN_LEFT);
  CHECKBOX_smokecullflag=glui_3dsmoke->add_checkbox_to_panel(panel5,"Cull hidden slices",&smokecullflag);
#ifdef _DEBUG
  CHECKBOX_smokedrawtest=glui_3dsmoke->add_checkbox_to_panel(panel5,"Show Only Back Slices",&smokedrawtest);
  CHECKBOX_smokedrawtest2=glui_3dsmoke->add_checkbox_to_panel(panel5,"Show Only Dir X Slices",&smokedrawtest2);


  SPINNER_smokedrawtest_nummin=glui_3dsmoke->add_spinner_to_panel(panel5,"Back Slice",GLUI_SPINNER_INT,&smokedrawtest_nummin);
  SPINNER_smokedrawtest_nummin->set_int_limits(1,ijkbarmax);

  SPINNER_smokedrawtest_nummax=glui_3dsmoke->add_spinner_to_panel(panel5,"Front Slice",GLUI_SPINNER_INT,&smokedrawtest_nummax);
  SPINNER_smokedrawtest_nummax->set_int_limits(1,ijkbarmax);
#endif
  skipframes = glui_3dsmoke->add_radiogroup_to_panel(panel5,&smokeskipm1);
  glui_3dsmoke->add_radiobutton_to_group(skipframes,"display all");
  glui_3dsmoke->add_radiobutton_to_group(skipframes,"   ... every 2nd");
  glui_3dsmoke->add_radiobutton_to_group(skipframes,"   ... every 3rd");




  panel_buttons = glui_3dsmoke->add_panel("",false);

  SPINNER_smoke3dframestep=glui_3dsmoke->add_spinner_to_panel(panel_buttons,"Frame Skip",GLUI_SPINNER_INT,&smoke3dframeskip,
    FRAMELOADING,SMOKE_3D_CB);
  SPINNER_smoke3dframestep->set_int_limits(0,100);

  glui_3dsmoke->add_column_to_panel(panel_buttons,false);

  glui_3dsmoke->add_button_to_panel(panel_buttons,"Save Settings",SAVE_SETTINGS,SMOKE_3D_CB);

  glui_3dsmoke->add_column_to_panel(panel_buttons,false);

  glui_3dsmoke->add_button_to_panel(panel_buttons,"Close",SMOKE_3D_CLOSE,SMOKE_3D_CB);
  
  glui_3dsmoke->set_main_gfx_window( main_window );
}

/* ------------------ hide_glui_stereo ------------------------ */

extern "C" void hide_glui_3dsmoke(void){
  if(glui_3dsmoke!=NULL)glui_3dsmoke->hide();
  showglui3dsmoke=0;
  updatemenu=1;
}

/* ------------------ show_glui_3dsmoke ------------------------ */

extern "C" void show_glui_3dsmoke(void){
  if(glui_3dsmoke!=NULL)glui_3dsmoke->show();
}


/* ------------------ show_glui_3dsmoke ------------------------ */

extern "C" void update_glui_smoke3dframestep(void){
//  SPINNER_smoke3dframestep->set_int_val(smoke3dframeskip);
}

void update_alpha(void){
  char label[100];
  char label1[100],label2[100],label3[100];
  float depth;
  float factor;

  factor = 1.0 - exp(-smoke_extinct*smoke_dens*smoke_pathlength);
  smoke_alpha = 255*factor;
  if(smoke_alpha<0)smoke_alpha=0;
  if(smoke_alpha>255)smoke_alpha=255;
  sprintf(label1,"%f",smoke_extinct);
  trimzeros(label1);
  sprintf(label2,"%f",smoke_dens);
  trimzeros(label2);
  sprintf(label3,"%f",smoke_pathlength);
  trimzeros(label3);
  sprintf(label,"alpha=%i=255*(1.0-exp(-%s*%s*%s))",smoke_alpha,label1,label2,label3);
  if(panel7!=NULL){
    TEXT_smokealpha->set_text(label);
  }
  
  if(smoke_extinct!=0.0&&smoke_dens!=0){
    depth=0.693147/(smoke_extinct*smoke_dens);
    sprintf(label,"50%s smoke depth=%f","%",depth);
  }
  else{
    sprintf(label,"50%s smoke depth=***","%");
  }
  if(panel7!=NULL){
    TEXT_smokedepth->set_text(label);
  }


}

/* ------------------ 3dsmoke_CB ------------------------ */


void SMOKE_3D_CB(int var){
  switch (var){
  case SMOKETEST:
    update_alpha();
    break;
  case FRAMELOADING:
    smoke3dframestep = smoke3dframeskip+1;
    update2_glui_smoke3dframestep();
    updatemenu=1;
    break;
  case SAVE_SETTINGS:
    writeini(LOCAL_INI);
    break;
  case SMOKE_3D_CLOSE:
    hide_glui_3dsmoke();
    break;
  case FIRE_RED:
  case FIRE_GREEN:
  case FIRE_BLUE:
  case FIRE_HALFDEPTH:
  case FIRE_CUTOFF:
  case FIRE_ALPHA:
  case SMOKE_SHADE:
  case SMOKE_THICK:
  glutPostRedisplay();
  force_redisplay=1;
  IDLE();

    break;
  default:
#ifdef _DEBUG
    abort();
#endif
    break;
  }
}

