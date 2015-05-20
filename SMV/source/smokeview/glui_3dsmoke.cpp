#define CPP
#include "options.h"

#include <stdio.h>
#include <string.h>
#include GLUT_H
#include <math.h>

#include "smokeviewvars.h"

extern GLUI_Rollout *ROLLOUT_smoke3d;
extern GLUI *glui_bounds;

extern "C" void Smoke3d_CB(int var);
extern "C" void init_volrender_surface(int firstcall);

#define FIRE_RED 1
#define FIRE_GREEN 2
#define FIRE_BLUE 3
#define FIRE_HALFDEPTH 4
#define GLOBAL_FIRE_CUTOFF 15
#define SMOKE_SHADE 7
#define SMOKE_COLORBAR_LIST 16
#define FIRECOLORMAP_TYPE 17
#define SHOW_FIRECOLORMAP 25
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
#define SMOKE_OPTIONS 19
#define LOAD_COMPRESSED_DATA 20
#define TEMP_MIN 21
#define TEMP_CUTOFF 22
#define TEMP_MAX 23
#define COMBINE_MESHES 24
#define NONGPU_VOL_FACTOR 26
#define GPU_VOL_FACTOR 27
#define GENERATE_IMAGES 28
#define START_FRAME 29
#define SKIP_FRAME 30
#define CANCEL_GENERATE_IMAGES 31
#define VOL_TOUR_LIST 46
#define VOL_PREFIX 47
#define VOL_UNLOAD_ALL 48

// two defines below are also defined elsewhere

#define SCRIPT_CANCEL_NOW 45
#define TOUR_LIST 24

void Script_CB(int var);
void TOUR_CB(int var);

GLUI *glui_3dsmoke=NULL;

GLUI_EditText *EDIT_vol_prefix=NULL;

GLUI_Listbox *LISTBOX_VOL_tour=NULL;

GLUI_Button *BUTTON_volunload=NULL;
GLUI_Button *BUTTON_startrender=NULL;
GLUI_Button *BUTTON_cancelrender=NULL;

GLUI_Listbox *LISTBOX_smoke_colorbar=NULL;

GLUI_RadioGroup *RADIO_alpha=NULL;
GLUI_RadioGroup *RADIO_render=NULL;
GLUI_RadioGroup *RADIO_skipframes=NULL;
GLUI_RadioGroup *RADIO_smokesensors=NULL;
GLUI_RadioGroup *RADIO_loadvol=NULL;
GLUI_RadioGroup *RADIO_use_colormap=NULL;

GLUI_RadioButton *RADIOBUTTON_direct=NULL,*RADIOBUTTON_constraint=NULL, *RADIOBUTTON_noconstraint=NULL;

GLUI_Spinner *SPINNER_startframe=NULL;
GLUI_Spinner *SPINNER_skipframe=NULL;
#ifdef pp_CULL
GLUI_Spinner *SPINNER_cull_portsize=NULL;
#endif
GLUI_Spinner *SPINNER_hrrpuv_cutoff=NULL;
GLUI_Spinner *SPINNER_nongpu_vol_factor=NULL;
GLUI_Spinner *SPINNER_gpu_vol_factor=NULL;

GLUI_Spinner *SPINNER_temperature_min=NULL;
GLUI_Spinner *SPINNER_temperature_cutoff=NULL;
GLUI_Spinner *SPINNER_temperature_max=NULL;
GLUI_Spinner *SPINNER_fire_opacity_factor=NULL;
GLUI_Spinner *SPINNER_mass_extinct=NULL;
GLUI_Spinner *SPINNER_cvis=NULL;
GLUI_Spinner *SPINNER_smokedrawtest_nummin=NULL;
GLUI_Spinner *SPINNER_smokedrawtest_nummax=NULL;
#ifdef pp_GPU
GLUI_Spinner *SPINNER_smoke3d_rthick=NULL;
#else
GLUI_Spinner *SPINNER_smoke3d_thick=NULL;
#endif
GLUI_Spinner *SPINNER_smoke3d_smoke_albedo=NULL;
GLUI_Spinner *SPINNER_smoke3d_fire_red=NULL;
GLUI_Spinner *SPINNER_smoke3d_fire_green=NULL;
GLUI_Spinner *SPINNER_smoke3d_fire_blue=NULL;
GLUI_Spinner *SPINNER_smoke3d_fire_halfdepth=NULL;
GLUI_Spinner *SPINNER_extinct=NULL;
GLUI_Spinner *SPINNER_smokedens=NULL;
GLUI_Spinner *SPINNER_pathlength=NULL;

GLUI_Checkbox *CHECKBOX_combine_meshes=NULL;
#ifdef pp_CULL
GLUI_Checkbox *CHECKBOX_show_cullports=NULL;
#endif
GLUI_Checkbox *CHECKBOX_usevolrender=NULL;
GLUI_Checkbox *CHECKBOX_compress_volsmoke=NULL;
GLUI_Checkbox *CHECKBOX_smokecullflag=NULL;
GLUI_Checkbox *CHECKBOX_test_smokesensors=NULL;
GLUI_Checkbox *CHECKBOX_smokeGPU=NULL;
GLUI_Checkbox *CHECKBOX_smokedrawtest=NULL;
GLUI_Checkbox *CHECKBOX_smokedrawtest2=NULL;
GLUI_Checkbox *CHECKBOX_zlib=NULL;
GLUI_Checkbox **CHECKBOX_meshvisptr=NULL;
GLUI_Checkbox *CHECKBOX_meshvis=NULL;
GLUI_Checkbox *CHECKBOX_show_smoketest=NULL;

GLUI_Panel *PANEL_overall=NULL;
GLUI_Panel *PANEL_colormap2=NULL;
GLUI_Panel *PANEL_colormap3a=NULL;
GLUI_Panel *PANEL_colormap3b=NULL;
GLUI_Panel *PANEL_colormap=NULL;
GLUI_Panel *PANEL_absorption=NULL,*PANEL_smokesensor=NULL;
GLUI_Panel *PANEL_testsmoke=NULL;

GLUI_Rollout *ROLLOUT_colormap3 = NULL;
GLUI_Rollout *PANEL_generate_images = NULL;
GLUI_Rollout *PANEL_colormap2a=NULL;
GLUI_Rollout *PANEL_colormap2b=NULL;
GLUI_Rollout *PANEL_meshvis=NULL;
GLUI_Rollout *ROLLOUT_slices=NULL;
GLUI_Rollout *ROLLOUT_volume=NULL;

GLUI_StaticText *TEXT_smokealpha=NULL;
GLUI_StaticText *TEXT_smokedepth=NULL;


extern "C" void UnLoadVolSmoke3DMenu(int var);

#define VOLRENDER_ROLLOUT 0
#define SLICERENDER_ROLLOUT 1

procdata smokeprocinfo[2];
int nsmokeprocinfo = 0;

/* ------------------ Smoke_Rollout_CB ------------------------ */

void Smoke_Rollout_CB(int var){
  if(var == VOLRENDER_ROLLOUT)smoke_render_option = 1;
    if(var==SLICERENDER_ROLLOUT)smoke_render_option=0;
    Smoke3d_CB(SMOKE_OPTIONS);
    toggle_rollout(smokeprocinfo, nsmokeprocinfo, var);
}

/* ------------------ delete_vol_tourlsit ------------------------ */

extern "C" void delete_vol_tourlist(void){
  int i;

  if(LISTBOX_VOL_tour==NULL)return;
  for(i=0;i<ntours;i++){
    LISTBOX_VOL_tour->delete_item(i);
  }
}

/* ------------------ create_vol_tourlist ------------------------ */

extern "C" void create_vol_tourlist(void){
  int i;
  
  if(LISTBOX_VOL_tour==NULL)return;
  for(i=0;i<ntours;i++){
    tourdata *touri;
    char label[1000];

    touri = tourinfo + i;
    strcpy(label,"");
    if(i==selectedtour_index)strcat(label,"*");
    if(touri->label!=NULL&&strlen(touri->label)>0)strcat(label,touri->label);
    if(strlen(label)>0){
      LISTBOX_VOL_tour->add_item(i,label);
    }
    else{
      LISTBOX_VOL_tour->add_item(i,"error");
    }
  }
  if(selectedtour_index>=-1&&selectedtour_index<ntours){
    LISTBOX_VOL_tour->set_int_val(selectedtour_index);
  }
}

/* ------------------ update_combine_meshes ------------------------ */

extern "C" void update_combine_meshes(void){
  if(CHECKBOX_combine_meshes!=NULL)CHECKBOX_combine_meshes->set_int_val(combine_meshes);
}

/* ------------------ update_gpu ------------------------ */

extern "C" void update_gpu(void){
#ifdef pp_GPU
  if(nsmoke3dinfo>0&&CHECKBOX_smokeGPU!=NULL){
    CHECKBOX_smokeGPU->set_int_val(usegpu);  
  }
#endif
}

/* ------------------ update_smoke3dflags ------------------------ */

extern "C" void update_smoke3dflags(void){
  RADIO_alpha->set_int_val(adjustalphaflag);
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
  RADIO_skipframes->set_int_val(smokeskipm1);
  Smoke3d_CB(VOL_SMOKE);
#ifdef pp_CULL
  Smoke3d_CB(CULL_SMOKE);
#endif
  glutPostRedisplay();
}

/* ------------------ Update_Smoke_Type ------------------------ */

void Update_Smoke_Type(void){
  switch(smoke_render_option){
    case RENDER_SLICE:
      if(ROLLOUT_slices!=NULL)ROLLOUT_slices->open();
      if(ROLLOUT_volume!=NULL)ROLLOUT_volume->close();
      break;
    case RENDER_VOLUME:
      if(ROLLOUT_slices!=NULL)ROLLOUT_slices->close();
      if(ROLLOUT_volume!=NULL)ROLLOUT_volume->open();
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
}

/* ------------------ update_apha ------------------------ */

void update_alpha(void){
  char label[100];
  char label1[100],label2[100],label3[100];
  float depth;
  float factor;

  factor = 1.0 - exp(-smoke_extinct*smoke_dens*smoke_pathlength);
  smoke_alpha = (int)(255*factor);
  if(smoke_alpha<0)smoke_alpha=0;
  if(smoke_alpha>255)smoke_alpha=255;
  sprintf(label1,"%f",smoke_extinct);
  trimzeros(label1);
  sprintf(label2,"%f",smoke_dens);
  trimzeros(label2);
  sprintf(label3,"%f",smoke_pathlength);
  trimzeros(label3);
  sprintf(label,"alpha=%i=255*(1.0-exp(-%s*%s*%s))",smoke_alpha,label1,label2,label3);
  if(PANEL_testsmoke!=NULL){
    TEXT_smokealpha->set_text(label);
  }
  
  if(smoke_extinct!=0.0&&smoke_dens!=0){
    depth=0.693147/(smoke_extinct*smoke_dens);
    sprintf(label,"50%s smoke depth=%f","%",depth);
  }
  else{
    sprintf(label,"50%s smoke depth=***","%");
  }
  if(PANEL_testsmoke!=NULL){
    TEXT_smokedepth->set_text(label);
  }
}

/* ------------------ glui_3dsmoke_setup ------------------------ */

extern "C" void glui_3dsmoke_setup(int main_window){

  int i;

  
  if(nsmoke3dinfo<=0&&nvolrenderinfo<=0)return;
  if(CHECKBOX_meshvisptr!=NULL)FREEMEMORY(CHECKBOX_meshvisptr);
  NewMemory((void **)&CHECKBOX_meshvisptr,nmeshes*sizeof(GLUI_Checkbox *));
  
  glui_3dsmoke=glui_bounds;

  if(smoketest==1){
    PANEL_testsmoke = glui_3dsmoke->add_panel_to_panel(ROLLOUT_smoke3d,_("Test smoke"));
    PANEL_testsmoke->set_alignment(GLUI_ALIGN_LEFT);
    CHECKBOX_show_smoketest=glui_3dsmoke->add_checkbox_to_panel(PANEL_testsmoke,_("Show test smoke"),&show_smoketest);
    SPINNER_extinct=glui_3dsmoke->add_spinner_to_panel(PANEL_testsmoke,_("Mass extinction coeff (m2/g)"),GLUI_SPINNER_FLOAT,&smoke_extinct,SMOKETEST,Smoke3d_CB);
    SPINNER_extinct->set_float_limits(0.0,10.0);
    SPINNER_smokedens=glui_3dsmoke->add_spinner_to_panel(PANEL_testsmoke,_("Smoke density (g/m3)"),GLUI_SPINNER_FLOAT,&smoke_dens,SMOKETEST,Smoke3d_CB);
    SPINNER_smokedens->set_float_limits(0.0,1.0);
    SPINNER_pathlength=glui_3dsmoke->add_spinner_to_panel(PANEL_testsmoke,_("Path length (m)"),GLUI_SPINNER_FLOAT,&smoke_pathlength,SMOKETEST,Smoke3d_CB);
    SPINNER_pathlength->set_float_limits(0.0,20.0);
    TEXT_smokealpha=glui_3dsmoke->add_statictext_to_panel(PANEL_testsmoke,_("Alpha"));
    TEXT_smokedepth=glui_3dsmoke->add_statictext_to_panel(PANEL_testsmoke,_("Depth"));
    update_alpha();
  }

  PANEL_overall = glui_3dsmoke->add_panel_to_panel(ROLLOUT_smoke3d,"",GLUI_PANEL_NONE);
#ifdef pp_GPU
  CHECKBOX_smokeGPU=glui_3dsmoke->add_checkbox_to_panel(PANEL_overall,_("Use GPU"),&usegpu,VOL_SMOKE,Smoke3d_CB);
#endif

  if(active_smokesensors==1){
    PANEL_smokesensor = glui_3dsmoke->add_panel_to_panel(PANEL_overall,_("Visibility"));
    RADIO_smokesensors = glui_3dsmoke->add_radiogroup_to_panel(PANEL_smokesensor,&show_smokesensors);
    glui_3dsmoke->add_radiobutton_to_group(RADIO_smokesensors,_("Hidden"));
    glui_3dsmoke->add_radiobutton_to_group(RADIO_smokesensors,_("Grey (0-255)"));
    glui_3dsmoke->add_radiobutton_to_group(RADIO_smokesensors,"I/I0 (0.0-1.0)");
    glui_3dsmoke->add_radiobutton_to_group(RADIO_smokesensors,_("Scaled optical depth (SCD)"));
    glui_3dsmoke->add_statictext_to_panel(PANEL_smokesensor,"SCD=C/K=C*L/Ln(I/I0) (0-Inf)");
    SPINNER_cvis=glui_3dsmoke->add_spinner_to_panel(PANEL_smokesensor,"C",GLUI_SPINNER_FLOAT,&smoke3d_cvis);
    SPINNER_cvis->set_float_limits(1.0,20.0);
#ifdef _DEBUG
    CHECKBOX_test_smokesensors=glui_3dsmoke->add_checkbox_to_panel(PANEL_smokesensor,"Test visibility sensor",&test_smokesensors);
#endif
  }

  PANEL_colormap = glui_3dsmoke->add_panel_to_panel(PANEL_overall,_("Color"));

  glui_3dsmoke->add_checkbox_to_panel(PANEL_colormap,"Show colormap",&show_firecolormap,SHOW_FIRECOLORMAP,Smoke3d_CB);
  RADIO_use_colormap = glui_3dsmoke->add_radiogroup_to_panel(PANEL_colormap,&firecolormap_type,FIRECOLORMAP_TYPE,Smoke3d_CB);
  RADIOBUTTON_direct=glui_3dsmoke->add_radiobutton_to_group(RADIO_use_colormap,"Use specified color, opacity, albedo");
  RADIOBUTTON_constraint=glui_3dsmoke->add_radiobutton_to_group(RADIO_use_colormap,"Use colormap with constraints");
  RADIOBUTTON_noconstraint=glui_3dsmoke->add_radiobutton_to_group(RADIO_use_colormap,"Use colormap without constraints");

  ROLLOUT_colormap3 = glui_3dsmoke->add_rollout_to_panel(PANEL_colormap,"fire color/opacity, smoke albedo",false);

  PANEL_colormap3a = glui_3dsmoke->add_panel_to_panel(ROLLOUT_colormap3,"",GLUI_PANEL_NONE);
  SPINNER_smoke3d_fire_red=glui_3dsmoke->add_spinner_to_panel(PANEL_colormap3a,_("red"),GLUI_SPINNER_INT,&fire_red,FIRE_RED,Smoke3d_CB);
  SPINNER_smoke3d_fire_red->set_int_limits(0,255);
  SPINNER_smoke3d_fire_green=glui_3dsmoke->add_spinner_to_panel(PANEL_colormap3a,_("green"),GLUI_SPINNER_INT,&fire_green,FIRE_GREEN,Smoke3d_CB);
  SPINNER_smoke3d_fire_green->set_int_limits(0,255);
  SPINNER_smoke3d_fire_blue=glui_3dsmoke->add_spinner_to_panel(PANEL_colormap3a,_("blue"),GLUI_SPINNER_INT,&fire_blue,FIRE_BLUE,Smoke3d_CB);
  SPINNER_smoke3d_fire_blue->set_int_limits(0,255);

  PANEL_colormap3b = glui_3dsmoke->add_panel_to_panel(ROLLOUT_colormap3,"",GLUI_PANEL_NONE);
  SPINNER_smoke3d_fire_halfdepth=glui_3dsmoke->add_spinner_to_panel(PANEL_colormap3b,_("fire half depth (m)"),GLUI_SPINNER_FLOAT,&fire_halfdepth,FIRE_HALFDEPTH,Smoke3d_CB);
  SPINNER_smoke3d_fire_halfdepth->set_float_limits(0.0,10.0);
  SPINNER_smoke3d_smoke_albedo=glui_3dsmoke->add_spinner_to_panel(PANEL_colormap3b,_("smoke albedo"),GLUI_SPINNER_FLOAT,&smoke_albedo,SMOKE_SHADE,Smoke3d_CB);
  SPINNER_smoke3d_smoke_albedo->set_float_limits(0.0,1.0);

  if(ncolorbars>0){
    LISTBOX_smoke_colorbar=glui_3dsmoke->add_listbox_to_panel(ROLLOUT_colormap3,_("colormap:"),&fire_colorbar_index,SMOKE_COLORBAR_LIST,Smoke3d_CB);

    for(i=0;i<ncolorbars;i++){
      colorbardata *cbi;

      cbi = colorbarinfo + i;
      cbi->label_ptr=cbi->label;
      LISTBOX_smoke_colorbar->add_item(i,cbi->label_ptr);
    }
    LISTBOX_smoke_colorbar->set_int_val(fire_colorbar_index);
  }

  PANEL_colormap2 = glui_3dsmoke->add_panel_to_panel(PANEL_colormap,"",GLUI_PANEL_NONE);

#define HRRPUV_CUTOFF_MAX (hrrpuv_max_smv-0.01)

  PANEL_colormap2b = glui_3dsmoke->add_rollout_to_panel(PANEL_colormap2,"HRRPUV (kW/m3)");
  SPINNER_hrrpuv_cutoff=glui_3dsmoke->add_spinner_to_panel(PANEL_colormap2b,_("cutoff"),GLUI_SPINNER_FLOAT,&global_hrrpuv_cutoff,GLOBAL_FIRE_CUTOFF,Smoke3d_CB);

  SPINNER_hrrpuv_cutoff->set_float_limits(0.0,HRRPUV_CUTOFF_MAX);

  {
    char label[256];

    strcpy(label,"Temperature (");
    strcat(label,degC);
    strcat(label,")");
    PANEL_colormap2a = glui_3dsmoke->add_rollout_to_panel(PANEL_colormap2,label);
  }
  SPINNER_temperature_min=glui_3dsmoke->add_spinner_to_panel(PANEL_colormap2a,_("min"),GLUI_SPINNER_FLOAT,
    &temperature_min,TEMP_MIN,Smoke3d_CB);
  SPINNER_temperature_cutoff=glui_3dsmoke->add_spinner_to_panel(PANEL_colormap2a,_("cutoff"),GLUI_SPINNER_FLOAT,
    &temperature_cutoff,TEMP_CUTOFF,Smoke3d_CB);
  SPINNER_temperature_max=glui_3dsmoke->add_spinner_to_panel(PANEL_colormap2a,_("max"),GLUI_SPINNER_FLOAT,
    &temperature_max,TEMP_MAX,Smoke3d_CB);

  Smoke3d_CB(TEMP_MIN);
  Smoke3d_CB(TEMP_CUTOFF);
  Smoke3d_CB(TEMP_MAX);

  if(nsmoke3dinfo>0){
    PANEL_meshvis = glui_3dsmoke->add_rollout_to_panel(PANEL_overall,"Mesh Visibility",false);
    for(i=0;i<nmeshes;i++){
      mesh *meshi;

      meshi = meshinfo + i;
      glui_3dsmoke->add_checkbox_to_panel(PANEL_meshvis,meshi->label,meshvisptr+i);
    }
  }

  Smoke3d_CB(FIRECOLORMAP_TYPE);

  glui_3dsmoke->add_column_to_panel(PANEL_overall,false);

  if(nsmoke3dinfo<=0||nvolrenderinfo<=0){
    smoke_render_option=RENDER_SLICE;
    if(nsmoke3dinfo>0)smoke_render_option=RENDER_SLICE;
    if(nvolrenderinfo>0)smoke_render_option=RENDER_VOLUME;
  }

  // slice render dialog
  
  if(nsmoke3dinfo>0){
    ROLLOUT_slices = glui_3dsmoke->add_rollout_to_panel(PANEL_overall,_("Slice render settings"),false, SLICERENDER_ROLLOUT, Smoke_Rollout_CB);
    ADDPROCINFO(smokeprocinfo, nsmokeprocinfo, ROLLOUT_slices, SLICERENDER_ROLLOUT);
    ROLLOUT_slices->set_alignment(GLUI_ALIGN_LEFT);
 
#ifdef pp_GPU
    if(gpuactive==0){
      usegpu=0;
      CHECKBOX_smokeGPU->disable();
    }
#endif
#ifdef pp_CULL
    CHECKBOX_smokecullflag=glui_3dsmoke->add_checkbox_to_panel(ROLLOUT_slices,_("Cull hidden slices"),&cullsmoke,CULL_SMOKE,Smoke3d_CB);
    if(cullactive==0){
      cullsmoke=0;
      CHECKBOX_smokecullflag->disable();
    }
    CHECKBOX_show_cullports=glui_3dsmoke->add_checkbox_to_panel(ROLLOUT_slices,_("Show cull ports"),&show_cullports);
    SPINNER_cull_portsize=glui_3dsmoke->add_spinner_to_panel(ROLLOUT_slices,_("Cull port size"),GLUI_SPINNER_INT,&cull_portsize,CULL_PORTSIZE,Smoke3d_CB);
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
    CHECKBOX_smokecullflag=glui_3dsmoke->add_checkbox_to_panel(ROLLOUT_slices,"Cull hidden slices",&smokecullflag);
#endif
#ifdef _DEBUG
    CHECKBOX_smokedrawtest=glui_3dsmoke->add_checkbox_to_panel(ROLLOUT_slices,"Show only back slices",&smokedrawtest);
    CHECKBOX_smokedrawtest2=glui_3dsmoke->add_checkbox_to_panel(ROLLOUT_slices,"Show only X slices",&smokedrawtest2);

    SPINNER_smokedrawtest_nummin=glui_3dsmoke->add_spinner_to_panel(ROLLOUT_slices,"Back slice",GLUI_SPINNER_INT,&smokedrawtest_nummin);
    SPINNER_smokedrawtest_nummin->set_int_limits(1,ijkbarmax);

    SPINNER_smokedrawtest_nummax=glui_3dsmoke->add_spinner_to_panel(ROLLOUT_slices,"Front slice",GLUI_SPINNER_INT,&smokedrawtest_nummax);
    SPINNER_smokedrawtest_nummax->set_int_limits(1,ijkbarmax);
#endif
    RADIO_skipframes = glui_3dsmoke->add_radiogroup_to_panel(ROLLOUT_slices,&smokeskipm1);
    glui_3dsmoke->add_radiobutton_to_group(RADIO_skipframes,_("Display all"));
    glui_3dsmoke->add_radiobutton_to_group(RADIO_skipframes,_("   ... Every 2nd"));
    glui_3dsmoke->add_radiobutton_to_group(RADIO_skipframes,_("   ... Every 3rd"));
#ifdef pp_GPU
    SPINNER_smoke3d_rthick=glui_3dsmoke->add_spinner_to_panel(ROLLOUT_slices,_("Thickness"),
      GLUI_SPINNER_FLOAT,&smoke3d_rthick,SMOKE_RTHICK,Smoke3d_CB);
    SPINNER_smoke3d_rthick->set_float_limits(1.0,255.0);
    smoke3d_thick = log_base2(smoke3d_rthick);
#else
    SPINNER_smoke3d_thick=glui_3dsmoke->add_spinner_to_panel(PANEL_colormap,"Thickness",
    GLUI_SPINNER_INT,&smoke3d_thick,SMOKE_THICK,Smoke3d_CB);
    SPINNER_smoke3d_thick->set_int_limits(0,7);
#endif

    PANEL_absorption = glui_3dsmoke->add_panel_to_panel(ROLLOUT_slices,_("Absorption adjustments"));
    PANEL_absorption->set_alignment(GLUI_ALIGN_LEFT);
    RADIO_alpha = glui_3dsmoke->add_radiogroup_to_panel(PANEL_absorption,&adjustalphaflag);
    glui_3dsmoke->add_radiobutton_to_group(RADIO_alpha,_("None"));
    glui_3dsmoke->add_radiobutton_to_group(RADIO_alpha,_("adjust off-center"));
    glui_3dsmoke->add_radiobutton_to_group(RADIO_alpha,_("zero at boundaries"));
    glui_3dsmoke->add_radiobutton_to_group(RADIO_alpha,_("both"));
  }


  // volume render dialog

  if(nvolrenderinfo > 0){
    ROLLOUT_volume = glui_3dsmoke->add_rollout_to_panel(PANEL_overall, _("Volume render settings"), false, VOLRENDER_ROLLOUT, Smoke_Rollout_CB);
    ADDPROCINFO(smokeprocinfo, nsmokeprocinfo, ROLLOUT_volume, VOLRENDER_ROLLOUT);

    if(have_volcompressed == 1){
      RADIO_loadvol = glui_3dsmoke->add_radiogroup_to_panel(ROLLOUT_volume, &glui_load_volcompressed, LOAD_COMPRESSED_DATA, Smoke3d_CB);
      glui_3dsmoke->add_radiobutton_to_group(RADIO_loadvol, _("Load full data"));
      glui_3dsmoke->add_radiobutton_to_group(RADIO_loadvol, _("Load compressed data"));
    }
    glui_3dsmoke->add_checkbox_to_panel(ROLLOUT_volume, _("Load data in background"), &use_multi_threading);
    CHECKBOX_compress_volsmoke = glui_3dsmoke->add_checkbox_to_panel(ROLLOUT_volume, _("Compress data while loading"), &glui_compress_volsmoke);
    if(have_volcompressed == 1){
      Smoke3d_CB(LOAD_COMPRESSED_DATA);
    }
    glui_3dsmoke->add_checkbox_to_panel(ROLLOUT_volume, _("Display data as b/w"), &volbw);
    glui_3dsmoke->add_checkbox_to_panel(ROLLOUT_volume, _("Show data while moving scene"), &show_volsmoke_moving);
    glui_3dsmoke->add_checkbox_to_panel(ROLLOUT_volume, _("Load data only at render times"), &load_at_rendertimes);

    SPINNER_fire_opacity_factor = glui_3dsmoke->add_spinner_to_panel(ROLLOUT_volume, _("Fire opacity multiplier"), GLUI_SPINNER_FLOAT, &fire_opacity_factor);
    SPINNER_fire_opacity_factor->set_float_limits(1.0, 10.0);
    SPINNER_mass_extinct = glui_3dsmoke->add_spinner_to_panel(ROLLOUT_volume, _("Mass extinction coeff"), GLUI_SPINNER_FLOAT, &mass_extinct);
    SPINNER_mass_extinct->set_float_limits(100.0, 100000.0);
    glui_3dsmoke->add_checkbox_to_panel(ROLLOUT_volume, _("Freeze"), &freeze_volsmoke);
#ifdef _DEBUG
    CHECKBOX_usevolrender = glui_3dsmoke->add_checkbox_to_panel(ROLLOUT_volume, _("Show"), &usevolrender, VOL_SMOKE, Smoke3d_CB);
    glui_3dsmoke->add_checkbox_to_panel(ROLLOUT_volume, "block smoke", &block_volsmoke);
    glui_3dsmoke->add_checkbox_to_panel(ROLLOUT_volume, "debug", &smoke3dVoldebug);
#endif
    CHECKBOX_combine_meshes = glui_3dsmoke->add_checkbox_to_panel(ROLLOUT_volume, _("Combine meshes"), &combine_meshes, COMBINE_MESHES, Smoke3d_CB);
    SPINNER_nongpu_vol_factor = glui_3dsmoke->add_spinner_to_panel(ROLLOUT_volume, _("non-gpu grid multiplier"), GLUI_SPINNER_FLOAT, &nongpu_vol_factor, NONGPU_VOL_FACTOR, Smoke3d_CB);
    SPINNER_nongpu_vol_factor->set_float_limits(1.0, 10.0);
    SPINNER_gpu_vol_factor = glui_3dsmoke->add_spinner_to_panel(ROLLOUT_volume, _("gpu grid multiplier"), GLUI_SPINNER_FLOAT, &gpu_vol_factor, GPU_VOL_FACTOR, Smoke3d_CB);
    SPINNER_gpu_vol_factor->set_float_limits(1.0, 10.0);

    PANEL_generate_images = glui_3dsmoke->add_rollout_to_panel(ROLLOUT_volume, _("Generate images"), false);

    SPINNER_startframe = glui_3dsmoke->add_spinner_to_panel(PANEL_generate_images, _("start frame"), GLUI_SPINNER_INT, &startframe0, START_FRAME, Smoke3d_CB);
    SPINNER_skipframe = glui_3dsmoke->add_spinner_to_panel(PANEL_generate_images, _("skip frame"), GLUI_SPINNER_INT, &skipframe0, SKIP_FRAME, Smoke3d_CB);
    Smoke3d_CB(START_FRAME);
    Smoke3d_CB(SKIP_FRAME);
    if(ntours > 0){
      selectedtour_index = -1;
      selectedtour_index_old = -1;
      LISTBOX_VOL_tour = glui_3dsmoke->add_listbox_to_panel(PANEL_generate_images, "Tour:", &selectedtour_index, VOL_TOUR_LIST, Smoke3d_CB);

      LISTBOX_VOL_tour->add_item(-1, "Manual");
      LISTBOX_VOL_tour->add_item(-999, "-");
      for(i = 0; i < ntours; i++){
        tourdata *touri;

        touri = tourinfo + i;
        LISTBOX_VOL_tour->add_item(i, touri->label);
      }
      LISTBOX_VOL_tour->set_int_val(selectedtour_index);
    }

    strcpy(vol_prefix, fdsprefix);
    EDIT_vol_prefix = glui_3dsmoke->add_edittext_to_panel(PANEL_generate_images, "image prefix:", GLUI_EDITTEXT_TEXT, vol_prefix, VOL_PREFIX, Smoke3d_CB);
    EDIT_vol_prefix->set_w(200);

    BUTTON_startrender = glui_3dsmoke->add_button_to_panel(PANEL_generate_images, _("Generate images"), GENERATE_IMAGES, Smoke3d_CB);
    BUTTON_cancelrender = glui_3dsmoke->add_button_to_panel(PANEL_generate_images, _("Cancel"), CANCEL_GENERATE_IMAGES, Smoke3d_CB);
    BUTTON_volunload = glui_3dsmoke->add_button_to_panel(PANEL_generate_images, _("Unload"), VOL_UNLOAD_ALL, Smoke3d_CB);
  }

  Update_Smoke_Type();

#ifdef pp_GPU
  Smoke3d_CB(VOL_SMOKE);
#endif
  Smoke3d_CB(SMOKE_OPTIONS);
}

/* ------------------ 3dsmoke_CB ------------------------ */

extern "C" void Smoke3d_CB(int var){
  int i;
  char *tour_label;
  char *vol_prefixptr;

  updatemenu=1;
  switch(var){
  float temp_min, temp_max;

  case VOL_UNLOAD_ALL:
    UnLoadVolSmoke3DMenu(-1);
    break;
  case VOL_PREFIX:
    break;
  case VOL_TOUR_LIST:
    TOUR_CB(TOUR_LIST);
    break;
  case START_FRAME:
    if(startframe0<0){
      startframe0=0;
      SPINNER_startframe->set_int_val(startframe0);
    }
    break;
  case SKIP_FRAME:
    if(skipframe0<1){
      skipframe0=1;
      SPINNER_skipframe->set_int_val(skipframe0);
    }
    break;
  case CANCEL_GENERATE_IMAGES:
    Script_CB(SCRIPT_CANCEL_NOW);
    break;
  case GENERATE_IMAGES:
    if(selected_tour==NULL){
      tour_label=NULL;
    }
    else{
      tour_label=selected_tour->label;
    }
    trim(vol_prefix);
    vol_prefixptr=trim_front(vol_prefix);
    if(strlen(vol_prefixptr)==0)vol_prefixptr=fdsprefix;
    init_volrender_script(vol_prefixptr, tour_label, startframe0, skipframe0);
    break;
  case NONGPU_VOL_FACTOR:
    init_volrender_surface(NOT_FIRSTCALL);
    break;
  case GPU_VOL_FACTOR:
    break;
  case COMBINE_MESHES:
    define_volsmoke_textures();
    break;
  case SHOW_FIRECOLORMAP:
    Update_Smokecolormap(smoke_render_option);
    if(show_firecolormap==1){
      show_glui_colorbar();
    }
    else{
      hide_glui_colorbar();
    }
    break;
  case TEMP_MIN:
    temp_min = 20.0;
    temp_max = (float)(10.0*(int)(temperature_cutoff/10.0)-10.0);
    SPINNER_temperature_min->set_float_limits(temp_min,temp_max);
    Update_Smokecolormap(smoke_render_option);
    break;
  case TEMP_CUTOFF:
    temp_min = (float)(10*(int)(temperature_min/10.0) + 10.0);
    temp_max = (float)(10*(int)(temperature_max/10.0) - 10.0);
    SPINNER_temperature_cutoff->set_float_limits(temp_min,temp_max);
    Update_Smokecolormap(smoke_render_option);
    break;
  case TEMP_MAX:
    temp_min = (float)(10*(int)(temperature_cutoff/10.0)+10.0);
    temp_max = 1800.0;
    SPINNER_temperature_max->set_float_limits(temp_min,temp_max);
    Update_Smokecolormap(smoke_render_option);
    break;
  case LOAD_COMPRESSED_DATA:
    if(load_volcompressed==1){
      CHECKBOX_compress_volsmoke->disable();
    }
    else{
      CHECKBOX_compress_volsmoke->enable();
    }
    break;
  case SMOKE_OPTIONS:
    if(firecolormap_type!=FIRECOLORMAP_NOCONSTRAINT&&smoke_render_option==RENDER_SLICE){
      PANEL_colormap2b->enable();
      SPINNER_hrrpuv_cutoff->enable();
      PANEL_colormap2b->open();
    }
    else{
      PANEL_colormap2b->disable();
      PANEL_colormap2b->close();
    }
    if(smoke_render_option==RENDER_SLICE){
      if(PANEL_colormap2a!=NULL){
        PANEL_colormap2a->disable();
        PANEL_colormap2a->close();
      }
      if(PANEL_absorption!=NULL)PANEL_absorption->enable();
      firecolormap_type=firecolormap_type_save;
      RADIO_use_colormap->set_int_val(firecolormap_type);
      RADIOBUTTON_direct->enable();
      SPINNER_smoke3d_fire_halfdepth->enable();
    }
    else{
      if(PANEL_colormap2a!=NULL){
        PANEL_colormap2a->enable();
        PANEL_colormap2a->open();
      }
      if(PANEL_absorption!=NULL)PANEL_absorption->disable();
      firecolormap_type_save=firecolormap_type;
      firecolormap_type=FIRECOLORMAP_CONSTRAINT;
      RADIO_use_colormap->set_int_val(firecolormap_type);
      RADIOBUTTON_direct->disable();
      SPINNER_smoke3d_fire_halfdepth->disable();
    }
    Smoke3d_CB(FIRECOLORMAP_TYPE);
    Update_Smoke_Type();
    break;
  case FIRECOLORMAP_TYPE:
    if(firecolormap_type==FIRECOLORMAP_CONSTRAINT&&smoke_render_option==RENDER_VOLUME){
      PANEL_colormap2a->open();
      PANEL_colormap2a->enable();
    }
    else{
      PANEL_colormap2a->close();
      PANEL_colormap2a->disable();
    }
    if(firecolormap_type!=FIRECOLORMAP_NOCONSTRAINT&&smoke_render_option==RENDER_SLICE){
      PANEL_colormap2b->enable();
      SPINNER_hrrpuv_cutoff->enable();
      PANEL_colormap2b->open();
    }
    else{
      PANEL_colormap2b->disable();
      PANEL_colormap2b->close();
    }
    if(firecolormap_type!=FIRECOLORMAP_DIRECT){
      LISTBOX_smoke_colorbar->enable();
      ROLLOUT_colormap3->disable();
      ROLLOUT_colormap3->close();
      if(fire_colorbar_index_save!=-1){
        SmokeColorBarMenu(fire_colorbar_index_save);
      }
      else{
        SmokeColorBarMenu(fire_colorbar_index);
      }
    }
    else{
      LISTBOX_smoke_colorbar->disable();
      ROLLOUT_colormap3->enable();
      ROLLOUT_colormap3->open();
      SPINNER_smoke3d_smoke_albedo->enable();
      SPINNER_smoke3d_fire_red->enable();
      SPINNER_smoke3d_fire_green->enable();
      SPINNER_smoke3d_fire_blue->enable();
      SPINNER_smoke3d_fire_halfdepth->enable();

      fire_colorbar_index_save=fire_colorbar_index;
    }
    if(LISTBOX_smoke_colorbar->get_int_val()!=fire_colorbar_index){
      LISTBOX_smoke_colorbar->set_int_val(fire_colorbar_index);
    }
    Update_Smokecolormap(smoke_render_option);
    break;
  case SMOKE_COLORBAR_LIST:
    SmokeColorBarMenu(fire_colorbar_index);
    Update_Smokecolormap(smoke_render_option);
    updatemenu=1;
    break;
  case SMOKETEST:
    update_alpha();
    break;
  case FRAMELOADING:
    smoke3dframestep = smoke3dframeskip+1;
    updatemenu=1;
    break;
  case SAVE_SETTINGS:
    writeini(LOCAL_INI,NULL);
    break;
  case GLOBAL_FIRE_CUTOFF:
    glutPostRedisplay();
    force_redisplay=1;
    Idle_CB();
    Update_Smokecolormap(smoke_render_option);
    break;
  case FIRE_RED:
  case FIRE_GREEN:
  case FIRE_BLUE:
  case SMOKE_SHADE:
    glutPostRedisplay();
    force_redisplay=1;
    Update_Smokecolormap(smoke_render_option);
    Idle_CB();
    break;
  case FIRE_HALFDEPTH:
    for(i=0;i<nmeshes;i++){
      mesh *meshi;

      meshi = meshinfo + i;
      meshi->update_firehalfdepth=1;
    }
    glutPostRedisplay();
    force_redisplay=1;
    Update_Smokecolormap(smoke_render_option);
    Idle_CB();
   break;  
#ifdef pp_GPU
  case SMOKE_RTHICK:
  
    smoke3d_thick = log_base2(smoke3d_rthick);
    glutPostRedisplay();
    force_redisplay=1;
    Idle_CB();
    break;
#else
  case SMOKE_THICK:
    glutPostRedisplay();
    force_redisplay=1;
    Idle_CB();
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
    {
      volrenderdata *vr;

      vr = &meshinfo->volrenderinfo;
      if(vr!=NULL&&vr->smokeslice!=NULL&&vr->smokeslice->slicetype==SLICE_CENTER){
        if(usegpu==1&&combine_meshes==1){
          combine_meshes=0;
          update_combine_meshes();
          Smoke3d_CB(COMBINE_MESHES);
        }
        if(usegpu==0&&combine_meshes==0){
          combine_meshes=1;
          update_combine_meshes();
          Smoke3d_CB(COMBINE_MESHES);
        }
      }
    }
    if(smoke_render_option==RENDER_SLICE){
#ifdef pp_GPU
      if(usegpu==1){
        RADIO_skipframes->set_int_val(0);
        RADIO_skipframes->disable();
#ifdef pp_CULL
        if(cullactive==1){
          CHECKBOX_smokecullflag->enable();
        }
        SPINNER_cull_portsize->enable();
        CHECKBOX_show_cullports->enable();
#endif
      }
      else{
        RADIO_skipframes->enable();
#ifdef pp_CULL
        CHECKBOX_smokecullflag->disable();
        SPINNER_cull_portsize->disable();
        CHECKBOX_show_cullports->disable();
#endif
      }
#else
      RADIO_skipframes->enable();
#endif
    }
    break;
  default:
#ifdef _DEBUG
    abort();
#endif
    break;
  }
}
