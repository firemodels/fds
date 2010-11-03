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
#include "MALLOC.h"
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

// svn revision character string
extern "C" char glui_3dsmoke_revision[]="$Revision$";

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
#define GLOBAL_FIRE_CUTOFF 15
#define SMOKE_SHADE 7
#ifdef pp_GPU
#define SMOKE_RTHICK 8
#else
#define SMOKE_THICK 8
#endif
#define SAVE_SETTINGS 9
#define FRAMELOADING 10
#define SMOKETEST 11
#ifdef pp_CULL
#define CULL_SMOKE 12
#define CULL_PORTSIZE 14
#endif
#ifdef pp_GPU
#define GPU_SMOKE 13
#endif

#ifdef pp_CULL
GLUI_Spinner *SPINNER_cull_portsize=NULL;
GLUI_Checkbox *CHECKBOX_show_cullports=NULL;
#endif
GLUI_Checkbox *CHECKBOX_smokecullflag=NULL;
GLUI *glui_3dsmoke=NULL;
GLUI_RadioGroup *alphagroup=NULL,*skipframes,*radio_smokesensors=NULL;
GLUI_Spinner *SPINNER_cvis=NULL;
GLUI_Checkbox *CHECKBOX_test_smokesensors=NULL;
#ifdef pp_GPU
GLUI_Checkbox *CHECKBOX_smokeGPU=NULL;
#endif
GLUI_Checkbox *CHECKBOX_smokedrawtest=NULL;
GLUI_Checkbox *CHECKBOX_smokedrawtest2=NULL;
GLUI_Checkbox *CHECKBOX_smoke3d_external=NULL;
GLUI_Checkbox *CHECKBOX_zlib=NULL;
GLUI_Spinner *SPINNER_smokedrawtest_nummin=NULL;
GLUI_Spinner *SPINNER_smoke3dframestep=NULL;
GLUI_Spinner *SPINNER_smokedrawtest_nummax=NULL;
#ifdef pp_GPU
GLUI_Spinner *SPINNER_smoke3d_rthick=NULL;
#else
GLUI_Spinner *SPINNER_smoke3d_thick=NULL;
#endif
GLUI_Spinner *SPINNER_smoke3d_smoke_shade=NULL;
GLUI_Spinner *SPINNER_smoke3d_fire_red=NULL;
GLUI_Spinner *SPINNER_smoke3d_fire_green=NULL;
GLUI_Spinner *SPINNER_smoke3d_fire_blue=NULL;
GLUI_Spinner *SPINNER_smoke3d_hrrpuv_cutoff=NULL;
GLUI_Spinner *SPINNER_smoke3d_fire_halfdepth=NULL;
GLUI_Spinner **SPINNER_smoke3d_hrrpuv_cutoffptr=NULL;
GLUI_Panel *panel_overall=NULL,*panel_fire=NULL,*panel_absorption=NULL,*panel_smokesensor=NULL;
GLUI_Panel *panel_slices=NULL,*panel_smoke=NULL,*panel_testsmoke=NULL;
GLUI_Spinner *SPINNER_extinct=NULL;
GLUI_Spinner *SPINNER_smokedens=NULL;
GLUI_Spinner *SPINNER_pathlength=NULL;
GLUI_StaticText *TEXT_smokealpha=NULL;
GLUI_StaticText *TEXT_smokedepth=NULL;
GLUI_Checkbox *CHECKBOX_show_smoketest=NULL;

/* ------------------ update_smoke3dflags ------------------------ */

extern "C" void update_smoke3dflags(void){
  alphagroup->set_int_val(adjustalphaflag);
#ifdef pp_GPU
  if(CHECKBOX_smokeGPU!=NULL)CHECKBOX_smokeGPU->set_int_val(usegpu);
#endif
#ifdef pp_CULL
  CHECKBOX_smokecullflag->set_int_val(cullsmoke);
#else
  CHECKBOX_smokecullflag->set_int_val(smokecullflag);
#endif
  if(CHECKBOX_smokedrawtest!=NULL)CHECKBOX_smokedrawtest->set_int_val(smokedrawtest);
  if(CHECKBOX_smokedrawtest2!=NULL)CHECKBOX_smokedrawtest2->set_int_val(smokedrawtest2);
  skipframes->set_int_val(smokeskipm1);
#ifdef pp_GPU
  SMOKE_3D_CB(GPU_SMOKE);
#endif
#ifdef pp_CULL
  SMOKE_3D_CB(CULL_SMOKE);
#endif
  GLUTPOSTREDISPLAY
}

/* ------------------ glui_3dsmoke_setup ------------------------ */

extern "C" void glui_3dsmoke_setup(int main_window){

  int i;

  
  if(nsmoke3d_files<=0)return;
  SPINNER_smoke3d_hrrpuv_cutoffptr=(GLUI_Spinner **)malloc((nmeshes+1)*sizeof(GLUI_Spinner *));
  
  glui_3dsmoke=glui_bounds;
  panel_absorption = glui_3dsmoke->add_panel_to_panel(panel_smoke3d,"Absorption Adjustments");
  panel_absorption->set_alignment(GLUI_ALIGN_LEFT);
  alphagroup = glui_3dsmoke->add_radiogroup_to_panel(panel_absorption,&adjustalphaflag);
  glui_3dsmoke->add_radiobutton_to_group(alphagroup,"none");
  glui_3dsmoke->add_radiobutton_to_group(alphagroup,"adjust off-center");
  glui_3dsmoke->add_radiobutton_to_group(alphagroup,"adjust off-center & zero at boundaries");
  glui_3dsmoke->add_radiobutton_to_group(alphagroup,"zero at boundaries");

  if(smoketest==1){
    panel_testsmoke = glui_3dsmoke->add_panel_to_panel(panel_smoke3d,"Test Smoke");
    panel_testsmoke->set_alignment(GLUI_ALIGN_LEFT);
    CHECKBOX_show_smoketest=glui_3dsmoke->add_checkbox_to_panel(panel_testsmoke,"Show Test Smoke",&show_smoketest);
    SPINNER_extinct=glui_3dsmoke->add_spinner_to_panel(panel_testsmoke,"Mass Extinction Coeff (m2/g)",GLUI_SPINNER_FLOAT,&smoke_extinct,SMOKETEST,SMOKE_3D_CB);
    SPINNER_extinct->set_float_limits(0.0,10.0);
    SPINNER_smokedens=glui_3dsmoke->add_spinner_to_panel(panel_testsmoke,"Smoke Density (g/m3)",GLUI_SPINNER_FLOAT,&smoke_dens,SMOKETEST,SMOKE_3D_CB);
    SPINNER_smokedens->set_float_limits(0.0,1.0);
    SPINNER_pathlength=glui_3dsmoke->add_spinner_to_panel(panel_testsmoke,"Path Length (m)",GLUI_SPINNER_FLOAT,&smoke_pathlength,SMOKETEST,SMOKE_3D_CB);
    SPINNER_pathlength->set_float_limits(0.0,20.0);
    TEXT_smokealpha=glui_3dsmoke->add_statictext_to_panel(panel_testsmoke,"Alpha");
    TEXT_smokedepth=glui_3dsmoke->add_statictext_to_panel(panel_testsmoke,"Depth");
    update_alpha();
  }




  panel_overall = glui_3dsmoke->add_panel_to_panel(panel_smoke3d,"",GLUI_PANEL_NONE);
  if(active_smokesensors==1){
    panel_smokesensor = glui_3dsmoke->add_panel_to_panel(panel_overall,"Visibility");
    radio_smokesensors = glui_3dsmoke->add_radiogroup_to_panel(panel_smokesensor,&show_smokesensors);
    glui_3dsmoke->add_radiobutton_to_group(radio_smokesensors,"Hidden");
    glui_3dsmoke->add_radiobutton_to_group(radio_smokesensors,"Grey (0-255)");
    glui_3dsmoke->add_radiobutton_to_group(radio_smokesensors,"I/I0 (0.0-1.0)");
    glui_3dsmoke->add_radiobutton_to_group(radio_smokesensors,"Scaled optical depth (SCD)");
    glui_3dsmoke->add_statictext_to_panel(panel_smokesensor,"SCD=C/K=C*L/Ln(I/I0) (0-Inf)");
    SPINNER_cvis=glui_3dsmoke->add_spinner_to_panel(panel_smokesensor,"C",GLUI_SPINNER_FLOAT,&smoke3d_cvis);
    SPINNER_cvis->set_float_limits(1.0,20.0);
#ifdef _DEBUG
    CHECKBOX_test_smokesensors=glui_3dsmoke->add_checkbox_to_panel(panel_smokesensor,"Test visibility sensor",&test_smokesensors);
#endif
  }

  panel_fire = glui_3dsmoke->add_panel_to_panel(panel_overall,"Fire");
  panel_fire->set_alignment(GLUI_ALIGN_LEFT);
  SPINNER_smoke3d_fire_red=glui_3dsmoke->add_spinner_to_panel(panel_fire,"red",GLUI_SPINNER_INT,&fire_red,FIRE_RED,SMOKE_3D_CB);
  SPINNER_smoke3d_fire_red->set_int_limits(0,255);
  SPINNER_smoke3d_fire_green=glui_3dsmoke->add_spinner_to_panel(panel_fire,"green",GLUI_SPINNER_INT,&fire_green,FIRE_GREEN,SMOKE_3D_CB);
  SPINNER_smoke3d_fire_green->set_int_limits(0,255);
  SPINNER_smoke3d_fire_blue=glui_3dsmoke->add_spinner_to_panel(panel_fire,"blue",GLUI_SPINNER_INT,&fire_blue,FIRE_BLUE,SMOKE_3D_CB);
  SPINNER_smoke3d_fire_blue->set_int_limits(0,255);
  SPINNER_smoke3d_fire_halfdepth=glui_3dsmoke->add_spinner_to_panel(panel_fire,"50% flame depth (m)",GLUI_SPINNER_FLOAT,&fire_halfdepth,FIRE_HALFDEPTH,SMOKE_3D_CB);
  SPINNER_smoke3d_fire_halfdepth->set_float_limits(0.0,10.0);


  {
    mesh *meshi;
    
#define HRRPUV_CUTOFF_MAX (hrrpuv_max_smv-0.01)

    glui_3dsmoke->add_statictext_to_panel(panel_fire,"HRRPUV cutoff (kW/m3):");

    if(nmeshes>1){
      SPINNER_smoke3d_hrrpuv_cutoffptr[nmeshes]=glui_3dsmoke->add_spinner_to_panel
        (panel_fire,"All meshes",GLUI_SPINNER_FLOAT,&global_hrrpuv_cutoff,GLOBAL_FIRE_CUTOFF,SMOKE_3D_CB);
      SPINNER_smoke3d_hrrpuv_cutoffptr[nmeshes]->set_float_limits(0.0,HRRPUV_CUTOFF_MAX);
    }
    else{
      SPINNER_smoke3d_hrrpuv_cutoffptr[nmeshes]=NULL;
    }
    for(i=0;i<nmeshes;i++){
      meshi = meshinfo + i;
      if(meshi->hrrpuv_cutoff>HRRPUV_CUTOFF_MAX)meshi->hrrpuv_cutoff=HRRPUV_CUTOFF_MAX;
      SPINNER_smoke3d_hrrpuv_cutoffptr[i]=glui_3dsmoke->add_spinner_to_panel
        (panel_fire,meshi->label,GLUI_SPINNER_FLOAT,&meshi->hrrpuv_cutoff,FIRE_CUTOFF,SMOKE_3D_CB);
      SPINNER_smoke3d_hrrpuv_cutoffptr[i]->set_float_limits(0.0,HRRPUV_CUTOFF_MAX);
    }
  }

  glui_3dsmoke->add_column_to_panel(panel_overall,false);

  panel_smoke = glui_3dsmoke->add_panel_to_panel(panel_overall,"Smoke");
  panel_smoke->set_alignment(GLUI_ALIGN_LEFT);
  SPINNER_smoke3d_smoke_shade=glui_3dsmoke->add_spinner_to_panel(panel_smoke,"Grey Level",GLUI_SPINNER_INT,&smoke_shade,SMOKE_SHADE,SMOKE_3D_CB);
  SPINNER_smoke3d_smoke_shade->set_int_limits(0,255);
#ifdef pp_GPU
  SPINNER_smoke3d_rthick=glui_3dsmoke->add_spinner_to_panel(panel_smoke,"Thickness",
    GLUI_SPINNER_FLOAT,&smoke3d_rthick,SMOKE_RTHICK,SMOKE_3D_CB);
  SPINNER_smoke3d_rthick->set_float_limits(1.0,255.0);
  smoke3d_thick = log_base2(smoke3d_rthick);
#else
  SPINNER_smoke3d_thick=glui_3dsmoke->add_spinner_to_panel(panel_smoke,"Thickness",
    GLUI_SPINNER_INT,&smoke3d_thick,SMOKE_THICK,SMOKE_3D_CB);
  SPINNER_smoke3d_thick->set_int_limits(0,7);
#endif
#ifdef _DEBUG
  CHECKBOX_smoke3d_external=glui_3dsmoke->add_checkbox_to_panel(panel_smoke,"View smoke externally (ONLY)",&smoke3d_external);
#endif



  panel_slices = glui_3dsmoke->add_panel_to_panel(panel_overall,"Slices");
  panel_slices->set_alignment(GLUI_ALIGN_LEFT);
#ifdef pp_GPU
  CHECKBOX_smokeGPU=glui_3dsmoke->add_checkbox_to_panel(panel_slices,"Use GPU",&usegpu,GPU_SMOKE,SMOKE_3D_CB);
  if(gpuactive==0){
    usegpu=0;
    CHECKBOX_smokeGPU->disable();
  }
#endif
#ifdef pp_CULL
  CHECKBOX_smokecullflag=glui_3dsmoke->add_checkbox_to_panel(panel_slices,"Cull hidden slices",&cullsmoke,CULL_SMOKE,SMOKE_3D_CB);
  if(cullactive==0){
    cullsmoke=0;
    CHECKBOX_smokecullflag->disable();
  }
  CHECKBOX_show_cullports=glui_3dsmoke->add_checkbox_to_panel(panel_slices,"Show cull ports",&show_cullports);
  SPINNER_cull_portsize=glui_3dsmoke->add_spinner_to_panel(panel_slices,"Cull port size",GLUI_SPINNER_INT,&cull_portsize,CULL_PORTSIZE,SMOKE_3D_CB);
  {
    int ijk_max=0;
    for(i=0;i<nmeshes;i++){
      mesh *meshi;

      meshi = meshinfo + i;
      if(ijk_max<meshi->ibar+1)ijk_max=meshi->ibar+1;
      if(ijk_max<meshi->jbar+1)ijk_max=meshi->jbar+1;
      if(ijk_max<meshi->kbar+1)ijk_max=meshi->kbar+1;
    }
    SPINNER_cull_portsize->set_int_limits(3,ijk_max);
  }
#else
  CHECKBOX_smokecullflag=glui_3dsmoke->add_checkbox_to_panel(panel_slices,"Cull hidden slices",&smokecullflag);
#endif
#ifdef _DEBUG
  CHECKBOX_smokedrawtest=glui_3dsmoke->add_checkbox_to_panel(panel_slices,"Show Only Back Slices",&smokedrawtest);
  CHECKBOX_smokedrawtest2=glui_3dsmoke->add_checkbox_to_panel(panel_slices,"Show Only Dir X Slices",&smokedrawtest2);


  SPINNER_smokedrawtest_nummin=glui_3dsmoke->add_spinner_to_panel(panel_slices,"Back Slice",GLUI_SPINNER_INT,&smokedrawtest_nummin);
  SPINNER_smokedrawtest_nummin->set_int_limits(1,ijkbarmax);

  SPINNER_smokedrawtest_nummax=glui_3dsmoke->add_spinner_to_panel(panel_slices,"Front Slice",GLUI_SPINNER_INT,&smokedrawtest_nummax);
  SPINNER_smokedrawtest_nummax->set_int_limits(1,ijkbarmax);
#endif
  skipframes = glui_3dsmoke->add_radiogroup_to_panel(panel_slices,&smokeskipm1);
  glui_3dsmoke->add_radiobutton_to_group(skipframes,"display all");
  glui_3dsmoke->add_radiobutton_to_group(skipframes,"   ... every 2nd");
  glui_3dsmoke->add_radiobutton_to_group(skipframes,"   ... every 3rd");

#ifdef pp_GPU
  SMOKE_3D_CB(GPU_SMOKE);
#endif
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
  if(panel_testsmoke!=NULL){
    TEXT_smokealpha->set_text(label);
  }
  
  if(smoke_extinct!=0.0&&smoke_dens!=0){
    depth=0.693147/(smoke_extinct*smoke_dens);
    sprintf(label,"50%s smoke depth=%f","%",depth);
  }
  else{
    sprintf(label,"50%s smoke depth=***","%");
  }
  if(panel_testsmoke!=NULL){
    TEXT_smokedepth->set_text(label);
  }


}

/* ------------------ 3dsmoke_CB ------------------------ */


void SMOKE_3D_CB(int var){
  int i;

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
  case GLOBAL_FIRE_CUTOFF:
    for(i=0;i<nmeshes;i++){
      SPINNER_smoke3d_hrrpuv_cutoffptr[i]->set_float_val(global_hrrpuv_cutoff);
    }
    GLUTPOSTREDISPLAY
    force_redisplay=1;
    IDLE();
    break;
  case FIRE_RED:
  case FIRE_GREEN:
  case FIRE_BLUE:
  case FIRE_HALFDEPTH:
  case FIRE_CUTOFF:
  case FIRE_ALPHA:
  case SMOKE_SHADE:
    GLUTPOSTREDISPLAY
    force_redisplay=1;
    IDLE();
    break;
#ifdef pp_GPU
  case SMOKE_RTHICK:
  
    smoke3d_thick = log_base2(smoke3d_rthick);
    GLUTPOSTREDISPLAY
    force_redisplay=1;
    IDLE();
    break;
#else
  case SMOKE_THICK:
    GLUTPOSTREDISPLAY
    force_redisplay=1;
    IDLE();
    break;
#endif
#ifdef pp_CULL
  case CULL_PORTSIZE:
    initcull(cullsmoke);
    break;
  case CULL_SMOKE:
    initcull(cullsmoke);
    break;
#endif
#ifdef pp_GPU
  case GPU_SMOKE:
    if(usegpu==1){
      skipframes->set_int_val(0);
      skipframes->disable();
#ifdef pp_CULL
      if(cullactive==1){
        CHECKBOX_smokecullflag->enable();
      }
      SPINNER_cull_portsize->enable();
      CHECKBOX_show_cullports->enable();
#endif
    }
    else{
      skipframes->enable();
#ifdef pp_CULL
      CHECKBOX_smokecullflag->disable();
      SPINNER_cull_portsize->disable();
      CHECKBOX_show_cullports->disable();
#endif
    }
    break;
#endif
  default:
#ifdef _DEBUG
    abort();
#endif
    break;
  }
}
