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
#include "translate.h"
#include "MALLOC.h"
#include "datadefs.h"

// svn revision character string
extern "C" char glui_3dsmoke_revision[]="$Revision$";

extern "C" void SmokeColorBarMenu(int fire_colorbar_index);
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
#define SMOKE_COLORBAR_LIST 16
#define USE_FIRESMOKEMAP 17
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
#define VOL_SMOKE 13
#define VOL_NGRID 18

GLUI_Listbox *LISTBOX_smoke_colorbar=NULL;
#ifdef pp_CULL
GLUI_Spinner *SPINNER_cull_portsize=NULL;
GLUI_Checkbox *CHECKBOX_show_cullports=NULL;
#endif
GLUI_Checkbox *CHECKBOX_use_firesmokemap=NULL;
GLUI_Checkbox *CHECKBOX_smokecullflag=NULL;
GLUI *glui_3dsmoke=NULL;
GLUI_RadioGroup *alphagroup=NULL,*skipframes,*radio_smokesensors=NULL;
GLUI_Spinner *SPINNER_cvis=NULL;
GLUI_Checkbox *CHECKBOX_test_smokesensors=NULL;
GLUI_Checkbox *CHECKBOX_smokeGPU=NULL,*CHECKBOX_usevolrender=NULL;
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
GLUI_Panel *panel_overall=NULL,*panel_colormap=NULL,*panel_absorption=NULL,*panel_smokesensor=NULL;
GLUI_Panel *panel_slices=NULL,*panel_testsmoke=NULL;
GLUI_Spinner *SPINNER_extinct=NULL;
GLUI_Spinner *SPINNER_smokedens=NULL;
GLUI_Spinner *SPINNER_pathlength=NULL;
GLUI_StaticText *TEXT_smokealpha=NULL,*STATIC_hrrpuvcolor=NULL;
GLUI_StaticText *TEXT_smokedepth=NULL;
GLUI_Checkbox *CHECKBOX_show_smoketest=NULL;

/* ------------------ update_smoke3dflags ------------------------ */

extern "C" void update_smoke3dflags(void){
  alphagroup->set_int_val(adjustalphaflag);
#ifdef pp_GPU
  if(CHECKBOX_smokeGPU!=NULL)CHECKBOX_smokeGPU->set_int_val(usegpu);
#endif
  if(CHECKBOX_usevolrender!=NULL)CHECKBOX_usevolrender->set_int_val(usevolrender);
#ifdef pp_CULL
  CHECKBOX_smokecullflag->set_int_val(cullsmoke);
#else
  CHECKBOX_smokecullflag->set_int_val(smokecullflag);
#endif
  if(CHECKBOX_smokedrawtest!=NULL)CHECKBOX_smokedrawtest->set_int_val(smokedrawtest);
  if(CHECKBOX_smokedrawtest2!=NULL)CHECKBOX_smokedrawtest2->set_int_val(smokedrawtest2);
  skipframes->set_int_val(smokeskipm1);
  SMOKE_3D_CB(VOL_SMOKE);
#ifdef pp_CULL
  SMOKE_3D_CB(CULL_SMOKE);
#endif
  glutPostRedisplay();
}

/* ------------------ glui_3dsmoke_setup ------------------------ */

extern "C" void glui_3dsmoke_setup(int main_window){

  int i;

  
  if(nsmoke3dinfo<=0)return;
  NewMemory((void **)&SPINNER_smoke3d_hrrpuv_cutoffptr,(nmeshes+1)*sizeof(GLUI_Spinner *));
  
  glui_3dsmoke=glui_bounds;

  if(smoketest==1){
    panel_testsmoke = glui_3dsmoke->add_panel_to_panel(panel_smoke3d,_("Test smoke"));
    panel_testsmoke->set_alignment(GLUI_ALIGN_LEFT);
    CHECKBOX_show_smoketest=glui_3dsmoke->add_checkbox_to_panel(panel_testsmoke,_("Show test smoke"),&show_smoketest);
    SPINNER_extinct=glui_3dsmoke->add_spinner_to_panel(panel_testsmoke,_("Mass extinction coeff (m2/g)"),GLUI_SPINNER_FLOAT,&smoke_extinct,SMOKETEST,SMOKE_3D_CB);
    SPINNER_extinct->set_float_limits(0.0,10.0);
    SPINNER_smokedens=glui_3dsmoke->add_spinner_to_panel(panel_testsmoke,_("Smoke density (g/m3)"),GLUI_SPINNER_FLOAT,&smoke_dens,SMOKETEST,SMOKE_3D_CB);
    SPINNER_smokedens->set_float_limits(0.0,1.0);
    SPINNER_pathlength=glui_3dsmoke->add_spinner_to_panel(panel_testsmoke,_("Path length (m)"),GLUI_SPINNER_FLOAT,&smoke_pathlength,SMOKETEST,SMOKE_3D_CB);
    SPINNER_pathlength->set_float_limits(0.0,20.0);
    TEXT_smokealpha=glui_3dsmoke->add_statictext_to_panel(panel_testsmoke,_("Alpha"));
    TEXT_smokedepth=glui_3dsmoke->add_statictext_to_panel(panel_testsmoke,_("Depth"));
    update_alpha();
  }




  panel_overall = glui_3dsmoke->add_panel_to_panel(panel_smoke3d,"",GLUI_PANEL_NONE);
  if(active_smokesensors==1){
    panel_smokesensor = glui_3dsmoke->add_panel_to_panel(panel_overall,_("Visibility"));
    radio_smokesensors = glui_3dsmoke->add_radiogroup_to_panel(panel_smokesensor,&show_smokesensors);
    glui_3dsmoke->add_radiobutton_to_group(radio_smokesensors,_("Hidden"));
    glui_3dsmoke->add_radiobutton_to_group(radio_smokesensors,_("Grey (0-255)"));
    glui_3dsmoke->add_radiobutton_to_group(radio_smokesensors,"I/I0 (0.0-1.0)");
    glui_3dsmoke->add_radiobutton_to_group(radio_smokesensors,_("Scaled optical depth (SCD)"));
    glui_3dsmoke->add_statictext_to_panel(panel_smokesensor,"SCD=C/K=C*L/Ln(I/I0) (0-Inf)");
    SPINNER_cvis=glui_3dsmoke->add_spinner_to_panel(panel_smokesensor,"C",GLUI_SPINNER_FLOAT,&smoke3d_cvis);
    SPINNER_cvis->set_float_limits(1.0,20.0);
#ifdef _DEBUG
    CHECKBOX_test_smokesensors=glui_3dsmoke->add_checkbox_to_panel(panel_smokesensor,"Test visibility sensor",&test_smokesensors);
#endif
  }

  panel_colormap = glui_3dsmoke->add_panel_to_panel(panel_overall,_("Fire/Smoke display"));

  CHECKBOX_use_firesmokemap=glui_3dsmoke->add_checkbox_to_panel(panel_colormap,_("Use colormap"),&use_firesmokemap,USE_FIRESMOKEMAP,SMOKE_3D_CB);
  if(ncolorbars>0){
    LISTBOX_smoke_colorbar=glui_3dsmoke->add_listbox_to_panel(panel_colormap,_("colormap:"),&fire_colorbar_index,SMOKE_COLORBAR_LIST,SMOKE_3D_CB);

    for(i=0;i<ncolorbars;i++){
      colorbardata *cbi;

      cbi = colorbarinfo + i;
      cbi->label_ptr=cbi->label;
      LISTBOX_smoke_colorbar->add_item(i,cbi->label_ptr);
    }
    LISTBOX_smoke_colorbar->set_int_val(fire_colorbar_index);
  }

  STATIC_hrrpuvcolor=glui_3dsmoke->add_statictext_to_panel(panel_colormap,_("fire color:"));
  SPINNER_smoke3d_fire_red=glui_3dsmoke->add_spinner_to_panel(panel_colormap,_("red"),GLUI_SPINNER_INT,&fire_red,FIRE_RED,SMOKE_3D_CB);
  SPINNER_smoke3d_fire_red->set_int_limits(0,255);
  SPINNER_smoke3d_fire_green=glui_3dsmoke->add_spinner_to_panel(panel_colormap,_("green"),GLUI_SPINNER_INT,&fire_green,FIRE_GREEN,SMOKE_3D_CB);
  SPINNER_smoke3d_fire_green->set_int_limits(0,255);
  SPINNER_smoke3d_fire_blue=glui_3dsmoke->add_spinner_to_panel(panel_colormap,_("blue"),GLUI_SPINNER_INT,&fire_blue,FIRE_BLUE,SMOKE_3D_CB);
  SPINNER_smoke3d_fire_blue->set_int_limits(0,255);

  SPINNER_smoke3d_smoke_shade=glui_3dsmoke->add_spinner_to_panel(panel_colormap,_("smoke albedo"),GLUI_SPINNER_FLOAT,&smoke_shade,SMOKE_SHADE,SMOKE_3D_CB);
  SPINNER_smoke3d_smoke_shade->set_float_limits(0.0,1.0);

  glui_3dsmoke->add_separator_to_panel(panel_colormap);


  {
    mesh *meshi;
    
#define HRRPUV_CUTOFF_MAX (hrrpuv_max_smv-0.01)

    glui_3dsmoke->add_statictext_to_panel(panel_colormap,_("HRRPUV cutoff (kW/m3):"));

//    SPINNER_smoke3d_fire_halfdepth=glui_3dsmoke->add_spinner_to_panel(panel_colormap,"50% fire depth (m)",GLUI_SPINNER_FLOAT,&fire_halfdepth,FIRE_HALFDEPTH,SMOKE_3D_CB);
//    SPINNER_smoke3d_fire_halfdepth->set_float_limits(0.0,10.0);
    if(nmeshes>1){
      SPINNER_smoke3d_hrrpuv_cutoffptr[nmeshes]=glui_3dsmoke->add_spinner_to_panel
        (panel_colormap,_("All meshes"),GLUI_SPINNER_FLOAT,&global_hrrpuv_cutoff,GLOBAL_FIRE_CUTOFF,SMOKE_3D_CB);
      SPINNER_smoke3d_hrrpuv_cutoffptr[nmeshes]->set_float_limits(0.0,HRRPUV_CUTOFF_MAX);
    }
    else{
      SPINNER_smoke3d_hrrpuv_cutoffptr[nmeshes]=NULL;
    }
    for(i=0;i<nmeshes;i++){
      meshi = meshinfo + i;
      if(meshi->hrrpuv_cutoff>HRRPUV_CUTOFF_MAX)meshi->hrrpuv_cutoff=HRRPUV_CUTOFF_MAX;
      SPINNER_smoke3d_hrrpuv_cutoffptr[i]=glui_3dsmoke->add_spinner_to_panel
        (panel_colormap,meshi->label,GLUI_SPINNER_FLOAT,&meshi->hrrpuv_cutoff,FIRE_CUTOFF,SMOKE_3D_CB);
      SPINNER_smoke3d_hrrpuv_cutoffptr[i]->set_float_limits(0.0,HRRPUV_CUTOFF_MAX);
    }
  }

  #ifdef pp_GPU
  SPINNER_smoke3d_rthick=glui_3dsmoke->add_spinner_to_panel(panel_colormap,_("Thickness"),
    GLUI_SPINNER_FLOAT,&smoke3d_rthick,SMOKE_RTHICK,SMOKE_3D_CB);
  SPINNER_smoke3d_rthick->set_float_limits(1.0,255.0);
  smoke3d_thick = log_base2(smoke3d_rthick);
#else
  SPINNER_smoke3d_thick=glui_3dsmoke->add_spinner_to_panel(panel_colormap,"Thickness",
    GLUI_SPINNER_INT,&smoke3d_thick,SMOKE_THICK,SMOKE_3D_CB);
  SPINNER_smoke3d_thick->set_int_limits(0,7);
#endif
#ifdef _DEBUG
  CHECKBOX_smoke3d_external=glui_3dsmoke->add_checkbox_to_panel(panel_colormap,"View smoke externally (ONLY)",&smoke3d_external);
#endif
  SMOKE_3D_CB(USE_FIRESMOKEMAP);

  glui_3dsmoke->add_column_to_panel(panel_overall,false);

  panel_absorption = glui_3dsmoke->add_panel_to_panel(panel_overall,_("Absorption adjustments"));
  panel_absorption->set_alignment(GLUI_ALIGN_LEFT);
  alphagroup = glui_3dsmoke->add_radiogroup_to_panel(panel_absorption,&adjustalphaflag);
  glui_3dsmoke->add_radiobutton_to_group(alphagroup,_("none"));
  glui_3dsmoke->add_radiobutton_to_group(alphagroup,_("adjust off-center"));
  glui_3dsmoke->add_radiobutton_to_group(alphagroup,_("zero at boundaries"));
  glui_3dsmoke->add_radiobutton_to_group(alphagroup,_("both"));

  panel_slices = glui_3dsmoke->add_panel_to_panel(panel_overall,"Slices");
  panel_slices->set_alignment(GLUI_ALIGN_LEFT);
#ifdef pp_GPU
  CHECKBOX_smokeGPU=glui_3dsmoke->add_checkbox_to_panel(panel_slices,_("Use GPU"),&usegpu,VOL_SMOKE,SMOKE_3D_CB);
#endif
  CHECKBOX_usevolrender=glui_3dsmoke->add_checkbox_to_panel(panel_slices,_("Use full volume rendering"),&usevolrender,VOL_SMOKE,SMOKE_3D_CB);

  glui_3dsmoke->add_checkbox_to_panel(panel_slices,"block smoke",&block_volsmoke);
  glui_3dsmoke->add_checkbox_to_panel(panel_slices,"debug",&smoke3dVoldebug);
 
#ifdef pp_GPU
  if(gpuactive==0){
    usegpu=0;
    CHECKBOX_smokeGPU->disable();
  }
#endif
#ifdef pp_CULL
  CHECKBOX_smokecullflag=glui_3dsmoke->add_checkbox_to_panel(panel_slices,_("Cull hidden slices"),&cullsmoke,CULL_SMOKE,SMOKE_3D_CB);
  if(cullactive==0){
    cullsmoke=0;
    CHECKBOX_smokecullflag->disable();
  }
  CHECKBOX_show_cullports=glui_3dsmoke->add_checkbox_to_panel(panel_slices,_("Show cull ports"),&show_cullports);
  SPINNER_cull_portsize=glui_3dsmoke->add_spinner_to_panel(panel_slices,_("Cull port size"),GLUI_SPINNER_INT,&cull_portsize,CULL_PORTSIZE,SMOKE_3D_CB);
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
  CHECKBOX_smokedrawtest=glui_3dsmoke->add_checkbox_to_panel(panel_slices,"Show only back slices",&smokedrawtest);
  CHECKBOX_smokedrawtest2=glui_3dsmoke->add_checkbox_to_panel(panel_slices,"Show only X slices",&smokedrawtest2);


  SPINNER_smokedrawtest_nummin=glui_3dsmoke->add_spinner_to_panel(panel_slices,"Back slice",GLUI_SPINNER_INT,&smokedrawtest_nummin);
  SPINNER_smokedrawtest_nummin->set_int_limits(1,ijkbarmax);

  SPINNER_smokedrawtest_nummax=glui_3dsmoke->add_spinner_to_panel(panel_slices,"Front slice",GLUI_SPINNER_INT,&smokedrawtest_nummax);
  SPINNER_smokedrawtest_nummax->set_int_limits(1,ijkbarmax);
#endif
  skipframes = glui_3dsmoke->add_radiogroup_to_panel(panel_slices,&smokeskipm1);
  glui_3dsmoke->add_radiobutton_to_group(skipframes,_("display all"));
  glui_3dsmoke->add_radiobutton_to_group(skipframes,_("   ... every 2nd"));
  glui_3dsmoke->add_radiobutton_to_group(skipframes,_("   ... every 3rd"));

#ifdef pp_GPU
  SMOKE_3D_CB(VOL_SMOKE);
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
  smoke_alpha=CLAMP(smoke_alpha,0,255);
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
  case USE_FIRESMOKEMAP:
    if(use_firesmokemap==1){
      LISTBOX_smoke_colorbar->enable();
      SPINNER_smoke3d_fire_red->disable();
      SPINNER_smoke3d_fire_green->disable();
      SPINNER_smoke3d_fire_blue->disable();
      SPINNER_smoke3d_smoke_shade->disable();
      STATIC_hrrpuvcolor->disable();
      if(fire_colorbar_index_save!=-1){
        SmokeColorBarMenu(fire_colorbar_index_save);
      }
      else{
        SmokeColorBarMenu(fire_colorbar_index);
      }
    }
    else{
      LISTBOX_smoke_colorbar->disable();
      SPINNER_smoke3d_fire_red->enable();
      SPINNER_smoke3d_fire_green->enable();
      SPINNER_smoke3d_fire_blue->enable();
      SPINNER_smoke3d_smoke_shade->enable();
      STATIC_hrrpuvcolor->enable();
      fire_colorbar_index_save=fire_colorbar_index;
      SmokeColorBarMenu(fire_custom_colorbar-colorbarinfo);
    }
    if(LISTBOX_smoke_colorbar->get_int_val()!=fire_colorbar_index){
      LISTBOX_smoke_colorbar->set_int_val(fire_colorbar_index);
    }
    break;
  case SMOKE_COLORBAR_LIST:
    SmokeColorBarMenu(fire_colorbar_index);
    updatemenu=1;
    break;
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
    glutPostRedisplay();
    force_redisplay=1;
    IDLE();
    break;
  case FIRE_RED:
  case FIRE_GREEN:
  case FIRE_BLUE:
  case SMOKE_SHADE:
    glutPostRedisplay();
    if(fire_custom_colorbar!=NULL){
      unsigned char *rgb_node;

      rgb_node=fire_custom_colorbar->rgb_node;
      rgb_node[0]=smoke_shade*255;
      rgb_node[1]=smoke_shade*255;
      rgb_node[2]=smoke_shade*255;
      rgb_node[3]=smoke_shade*255;
      rgb_node[4]=smoke_shade*255;
      rgb_node[5]=smoke_shade*255;
      rgb_node[6]=fire_red;
      rgb_node[7]=fire_green;
      rgb_node[8]=fire_blue;
      rgb_node[9]=fire_red;
      rgb_node[10]=fire_green;
      rgb_node[11]=fire_blue;
      remapcolorbar(fire_custom_colorbar);
      updatecolors(-1);
    }
    force_redisplay=1;
    IDLE();
    break;
  case FIRE_HALFDEPTH:
  case FIRE_CUTOFF:
  case FIRE_ALPHA:
    glutPostRedisplay();
    force_redisplay=1;
    IDLE();
    break;
#ifdef pp_GPU
  case SMOKE_RTHICK:
  
    smoke3d_thick = log_base2(smoke3d_rthick);
    glutPostRedisplay();
    force_redisplay=1;
    IDLE();
    break;
#else
  case SMOKE_THICK:
    glutPostRedisplay();
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
  case VOL_NGRID:
    glutPostRedisplay();
    break;
  case VOL_SMOKE:
#ifdef pp_GPU
    if(usegpu==1&&usevolrender==0){
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
      if(usevolrender==1){
        skipframes->disable();
      }
      else{
        skipframes->enable();
      }
#ifdef pp_CULL
      CHECKBOX_smokecullflag->disable();
      SPINNER_cull_portsize->disable();
      CHECKBOX_show_cullports->disable();
#endif
    }
    break;
#else
    if(usevolrender==1){
      skipframes->disable();
    }
    else{
      skipframes->enable();
    }
#endif
    break;
  default:
#ifdef _DEBUG
    abort();
#endif
    break;
  }
}
