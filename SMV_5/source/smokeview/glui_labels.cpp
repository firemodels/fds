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

// svn revision character string2
extern "C" char glui_labels_revision[]="$Revision$";

extern "C" void Labels_CB(int value);
int nevacloaded,nplot3dloaded,nsmoke3dloaded,nisoloaded,nsliceloaded,nvsliceloaded,npartloaded,npatchloaded;

GLUI *glui_labels=NULL;

GLUI_Spinner *SPINNER_tick_xmin=NULL;
GLUI_Spinner *SPINNER_tick_ymin=NULL;
GLUI_Spinner *SPINNER_tick_zmin=NULL;
GLUI_Spinner *SPINNER_tick_xmax=NULL;
GLUI_Spinner *SPINNER_tick_ymax=NULL;
GLUI_Spinner *SPINNER_tick_zmax=NULL;

GLUI_Rollout *panel_user_tick=NULL;
GLUI_Rollout *panel_label1=NULL;
GLUI_Panel *panel_label2=NULL;
GLUI_Panel *panel_tick1;
GLUI_Panel *panel_tick1a;
GLUI_Panel *panel_tick1b;
GLUI_Panel *panel_tick2;
GLUI_Spinner *SPINNER_sensorrelsize=NULL;
GLUI_Spinner *SPINNER_tick_x0=NULL;
GLUI_Spinner *SPINNER_tick_y0=NULL;
GLUI_Spinner *SPINNER_tick_z0=NULL;
GLUI_Spinner *SPINNER_tick_dx0=NULL;
GLUI_Spinner *SPINNER_tick_dy0=NULL;
GLUI_Spinner *SPINNER_tick_dz0=NULL;
GLUI_Panel *panel_transparency=NULL;
GLUI_Spinner *SPINNER_labels_transparency_face=NULL;
GLUI_Checkbox *CHECKBOX_labels_colorbar=NULL;
GLUI_Checkbox *CHECKBOX_labels_timebar=NULL;
GLUI_Checkbox *CHECKBOX_labels_ticks=NULL;
GLUI_Checkbox *CHECKBOX_labels_title=NULL;
GLUI_Checkbox *CHECKBOX_labels_axis=NULL;
GLUI_Checkbox *CHECKBOX_labels_hms=NULL;
GLUI_Spinner *SPINNER_subtick=NULL;

GLUI_Checkbox *CHECKBOX_labels_framerate=NULL;
GLUI_Checkbox *CHECKBOX_labels_timelabel=NULL;
GLUI_Checkbox *CHECKBOX_labels_framelabel=NULL;
GLUI_Checkbox *CHECKBOX_labels_hrrlabel=NULL;
GLUI_Checkbox *CHECKBOX_labels_hrrcutoff=NULL;
GLUI_Checkbox *CHECKBOX_labels_availmemory=NULL;
GLUI_Checkbox *CHECKBOX_labels_labels=NULL;
GLUI_Checkbox *CHECKBOX_labels_gridloc=NULL;
GLUI_Checkbox *CHECKBOX_labels_average=NULL;
GLUI_Checkbox *CHECKBOX_vis_user_ticks=NULL;
GLUI_Checkbox *CHECKBOX_user_ticks_show_x=NULL;
GLUI_Checkbox *CHECKBOX_user_ticks_show_y=NULL;
GLUI_Checkbox *CHECKBOX_user_ticks_show_z=NULL;
GLUI_Checkbox *CHECKBOX_tick_auto=NULL;

GLUI_Checkbox *CHECKBOX_labels_flip=NULL;
GLUI_Checkbox *CHECKBOX_labels_shade=NULL;
GLUI_Checkbox *CHECKBOX_labels_transparent_override=NULL;

GLUI_RadioGroup *RADIO_fontsize=NULL,*RADIO_showhide=NULL;
GLUI_Button *Button_EVAC=NULL;
GLUI_Button *Button_PART=NULL;
GLUI_Button *Button_SLICE=NULL;
GLUI_Button *Button_VSLICE=NULL;
GLUI_Button *Button_PLOT3D=NULL;
GLUI_Button *Button_3DSMOKE=NULL;
GLUI_Button *Button_BOUNDARY=NULL;
GLUI_Button *Button_ISO=NULL;
GLUI_Button *Button_BENCHMARK=NULL;

GLUI_Panel *panel_showhide=NULL;

#define LABELS_label 0
#define FRAME_label 21
#define HRR_label 22
#define HRRPUVCUTOFF_label 23
#define LABELS_showall 1
#define LABELS_hideall 2
#define LABELS_close 3
#define LABELS_flip 4
#define LABELS_shade 5
#define LABELS_transparent 6
#define LABELS_fontsize 7
#define LABELS_ticks 8
#define LABELS_sensorsize 20

#define LABELS_particleshow    10
#define LABELS_sliceshow       11
#define LABELS_vsliceshow      12
#define LABELS_boundaryshow    13
#define LABELS_3dsmokeshow     14
#define LABELS_isosurfaceshow  15
#define LABELS_evacshow 19
#define LABELS_PLOT3D 16
#define LABELS_BENCHMARK 17
#define LABELS_HMS 18
#define SAVE_SETTINGS 99

/* ------------------ glui_labels_setup ------------------------ */

extern "C" void glui_labels_setup(int main_window){

  if(glui_labels!=NULL)glui_labels->close();
  glui_labels = GLUI_Master.create_glui("Display",0,0,0);
  if(showlabels==0)glui_labels->hide();

  panel_label1 = glui_labels->add_rollout("General Settings",true);
  CHECKBOX_labels_colorbar=glui_labels->add_checkbox_to_panel(panel_label1,"Color Bar",&visColorLabels,LABELS_label,Labels_CB);
  CHECKBOX_labels_timebar=glui_labels->add_checkbox_to_panel(panel_label1,"Time Bar",&visTimeLabels,LABELS_label,Labels_CB);
  CHECKBOX_labels_timelabel=glui_labels->add_checkbox_to_panel(panel_label1,"Time Label",&visTimelabel,LABELS_label,Labels_CB);
  CHECKBOX_labels_framelabel=glui_labels->add_checkbox_to_panel(panel_label1,"Frame Label",&visFramelabel,FRAME_label,Labels_CB);
  CHECKBOX_labels_hrrlabel=glui_labels->add_checkbox_to_panel(panel_label1,"HRR Label",&visHRRlabel,HRR_label,Labels_CB);
  CHECKBOX_labels_hrrcutoff=glui_labels->add_checkbox_to_panel(panel_label1,"HRRPUV cutoff",&show_hrrcutoff,HRRPUVCUTOFF_label,Labels_CB);
  CHECKBOX_labels_ticks=glui_labels->add_checkbox_to_panel(panel_label1,"FDS Ticks",&visTicks,LABELS_label,Labels_CB);
  if(ntotal_blockages>0||isZoneFireModel==0){
    CHECKBOX_labels_gridloc=glui_labels->add_checkbox_to_panel(panel_label1,"Grid Loc",&visgridloc,LABELS_label,Labels_CB);
  }
  if(nslice_files>0)CHECKBOX_labels_average=glui_labels->add_checkbox_to_panel(panel_label1,"Average",&vis_slice_average,LABELS_label,Labels_CB);
  glui_labels->add_column_to_panel(panel_label1,false);


  CHECKBOX_labels_title=glui_labels->add_checkbox_to_panel(panel_label1,"Title",&visTitle0,LABELS_label,Labels_CB);
  CHECKBOX_labels_axis=glui_labels->add_checkbox_to_panel(panel_label1,"Axis",&visaxislabels,LABELS_label,Labels_CB);
  CHECKBOX_labels_framerate=glui_labels->add_checkbox_to_panel(panel_label1,"Frame Rate",&visFramerate,LABELS_label,Labels_CB);
#ifdef pp_memstatus
  CHECKBOX_labels_availmemory=glui_labels->add_checkbox_to_panel(panel_label1,"Memory Load",&visAvailmemory,LABELS_label,Labels_CB);
#endif
  CHECKBOX_labels_labels=glui_labels->add_checkbox_to_panel(panel_label1,"Text labels",&visLabels,LABELS_label,Labels_CB);
  glui_labels->add_button_to_panel(panel_label1,"Show All",LABELS_showall,Labels_CB);
  glui_labels->add_button_to_panel(panel_label1,"Hide All",LABELS_hideall,Labels_CB);


  glui_labels->add_column_to_panel(panel_label1,true);

  CHECKBOX_labels_flip=glui_labels->add_checkbox_to_panel(panel_label1,"Flip Background",&background_flip,LABELS_flip,Labels_CB);
  CHECKBOX_labels_shade=glui_labels->add_checkbox_to_panel(panel_label1,"Shades of Grey",&setbw,LABELS_shade,Labels_CB);

  if(nface_transparent>0){
    panel_transparency = glui_labels->add_panel_to_panel(panel_label1,"Geometry Transparency");
    CHECKBOX_labels_transparent_override=glui_labels->add_checkbox_to_panel(panel_transparency,"Use level:",&use_transparency_geom,LABELS_transparent,Labels_CB);
    SPINNER_labels_transparency_face=glui_labels->add_spinner_to_panel(panel_transparency,"",GLUI_SPINNER_FLOAT,&transparency_geom,LABELS_transparent,Labels_CB);
    SPINNER_labels_transparency_face->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);
    Labels_CB(LABELS_transparent);
  }


  CHECKBOX_labels_hms=glui_labels->add_checkbox_to_panel(panel_label1,"hms time label",&vishmsTimelabel,LABELS_HMS,Labels_CB);

  RADIO_fontsize = glui_labels->add_radiogroup_to_panel(panel_label1,&fontindex,LABELS_fontsize,Labels_CB);
  glui_labels->add_radiobutton_to_group(RADIO_fontsize,"small font");
  glui_labels->add_radiobutton_to_group(RADIO_fontsize,"large font");

  SPINNER_sensorrelsize=glui_labels->add_spinner_to_panel(panel_label1,"Sensor Scaling",GLUI_SPINNER_FLOAT,&sensorrelsize,LABELS_sensorsize,Labels_CB);

  panel_user_tick = glui_labels->add_rollout("User Tick Settings",false);


  panel_tick1 = glui_labels->add_panel_to_panel(panel_user_tick,"Display",true);
  panel_tick1a = glui_labels->add_panel_to_panel(panel_tick1,"",false);

  CHECKBOX_vis_user_ticks=glui_labels->add_checkbox_to_panel(panel_tick1a,"Show User Ticks",&vis_user_ticks);
  glui_labels->add_column_to_panel(panel_tick1a,false);
  SPINNER_subtick=glui_labels->add_spinner_to_panel(panel_tick1a,"sub-intervals",GLUI_SPINNER_INT,&user_tick_sub); 
  SPINNER_subtick->set_int_limits(1,10,GLUI_LIMIT_CLAMP);

  panel_tick1b = glui_labels->add_panel_to_panel(panel_tick1,"",false);
  CHECKBOX_tick_auto=glui_labels->add_checkbox_to_panel(panel_tick1b,"Auto place (2D)",&auto_user_tick_placement,LABELS_ticks,Labels_CB);
  glui_labels->add_column_to_panel(panel_tick1b,false);
  CHECKBOX_user_ticks_show_x=glui_labels->add_checkbox_to_panel(panel_tick1b,"x",&user_tick_show_x);
  glui_labels->add_column_to_panel(panel_tick1b,false);
  CHECKBOX_user_ticks_show_y=glui_labels->add_checkbox_to_panel(panel_tick1b,"y",&user_tick_show_y);
  glui_labels->add_column_to_panel(panel_tick1b,false);
  CHECKBOX_user_ticks_show_z=glui_labels->add_checkbox_to_panel(panel_tick1b,"z",&user_tick_show_z);
  Labels_CB(LABELS_ticks);

  panel_tick2 = glui_labels->add_panel_to_panel(panel_user_tick,"Parameters",true);
  glui_labels->add_statictext_to_panel(panel_tick2,"                    x");
  SPINNER_tick_x0=glui_labels->add_spinner_to_panel(panel_tick2,"origin",GLUI_SPINNER_FLOAT,user_tick_origin);
  SPINNER_tick_xmin=glui_labels->add_spinner_to_panel(panel_tick2,"min",GLUI_SPINNER_FLOAT,user_tick_min);
  SPINNER_tick_xmax=glui_labels->add_spinner_to_panel(panel_tick2,"max",GLUI_SPINNER_FLOAT,user_tick_max);
  SPINNER_tick_dx0=glui_labels->add_spinner_to_panel(panel_tick2,"step",GLUI_SPINNER_FLOAT,user_tick_step);

  glui_labels->add_column_to_panel(panel_tick2,false);

  glui_labels->add_statictext_to_panel(panel_tick2,"                    y");
  SPINNER_tick_y0=glui_labels->add_spinner_to_panel(panel_tick2,"",GLUI_SPINNER_FLOAT,user_tick_origin+1);
  SPINNER_tick_ymin=glui_labels->add_spinner_to_panel(panel_tick2,"",GLUI_SPINNER_FLOAT,user_tick_min+1);
  SPINNER_tick_ymax=glui_labels->add_spinner_to_panel(panel_tick2,"",GLUI_SPINNER_FLOAT,user_tick_max+1);
  SPINNER_tick_dy0=glui_labels->add_spinner_to_panel(panel_tick2,"",GLUI_SPINNER_FLOAT,user_tick_step+1);

  glui_labels->add_column_to_panel(panel_tick2,false);

  glui_labels->add_statictext_to_panel(panel_tick2,"                    z");
  SPINNER_tick_z0=glui_labels->add_spinner_to_panel(panel_tick2,"",GLUI_SPINNER_FLOAT,user_tick_origin+2);
  SPINNER_tick_zmin=glui_labels->add_spinner_to_panel(panel_tick2,"",GLUI_SPINNER_FLOAT,user_tick_min+2);
  SPINNER_tick_zmax=glui_labels->add_spinner_to_panel(panel_tick2,"",GLUI_SPINNER_FLOAT,user_tick_max+2);
  SPINNER_tick_dz0=glui_labels->add_spinner_to_panel(panel_tick2,"",GLUI_SPINNER_FLOAT,user_tick_step+2);
  
  if((npart_files>0)||nslice_files>0||nvslice>0||niso_files>0||npatch_files||nsmoke3d_files>0||nplot3d_files>0){
    panel_showhide = glui_labels->add_rollout("Show/Hide Loaded Files",false);

    RADIO_showhide = glui_labels->add_radiogroup_to_panel(panel_showhide,&showhide_option);
    glui_labels->add_radiobutton_to_group(RADIO_showhide,"Show");
    glui_labels->add_radiobutton_to_group(RADIO_showhide,"Show Only");
    glui_labels->add_radiobutton_to_group(RADIO_showhide,"Hide");

    glui_labels->add_column_to_panel(panel_showhide,false);

    if(nevac>0){}
    if(npart_files>0&&nevac!=npart_files)Button_PART=glui_labels->add_button_to_panel(panel_showhide,"Particle",LABELS_particleshow,Labels_CB);
    if(nevac>0)Button_EVAC=glui_labels->add_button_to_panel(panel_showhide,"Evacuation",LABELS_evacshow,Labels_CB);
    if(nslice_files>0)Button_SLICE=glui_labels->add_button_to_panel(panel_showhide,"Slice",LABELS_sliceshow,Labels_CB);
    if(nvslice>0)Button_VSLICE=glui_labels->add_button_to_panel(panel_showhide,"Vector",LABELS_vsliceshow,Labels_CB);
    if(niso_files>0)Button_ISO=glui_labels->add_button_to_panel(panel_showhide,"Isosurface",LABELS_isosurfaceshow,Labels_CB);
    if(npatch_files>0)Button_BOUNDARY=glui_labels->add_button_to_panel(panel_showhide,"Boundary",LABELS_boundaryshow,Labels_CB);
    if(nsmoke3d_files>0)Button_3DSMOKE=glui_labels->add_button_to_panel(panel_showhide,"3D Smoke",LABELS_3dsmokeshow,Labels_CB);
    if(nplot3d_files>0)Button_PLOT3D=glui_labels->add_button_to_panel(panel_showhide,"Plot3D",LABELS_PLOT3D,Labels_CB);

    update_showhidebuttons();
  }

  panel_label2 = glui_labels->add_panel("",false);
  Button_BENCHMARK=glui_labels->add_button_to_panel(panel_label2,"Benchmark",LABELS_BENCHMARK,Labels_CB);
  glui_labels->add_column_to_panel(panel_label2,false);

  glui_labels->add_button_to_panel(panel_label2,"Save Settings",SAVE_SETTINGS,Labels_CB);
  glui_labels->add_column_to_panel(panel_label2,false);

  glui_labels->add_button_to_panel(panel_label2,"Close",LABELS_close,Labels_CB);

  glui_labels->set_main_gfx_window( main_window );
}

/* ------------------ update_fileload  ------------------------ */

extern "C" void update_fileload(void){
  int i;
  particle *parti;
  slice *slicei;
  iso *isoi;
  patch *patchi;
  smoke3d *smoke3di;
  plot3d *plot3di;
  vslice *vslicei;

  npartloaded=0;
  nevacloaded=0;
  for(i=0;i<npart_files;i++){
    parti = partinfo+i;
    if(parti->loaded==1&&parti->evac==0){
      npartloaded++;
    }
    if(parti->loaded==1&&parti->evac==1){
      nevacloaded++;
    }
  }

  nsliceloaded=0;
  for(i=0;i<nslice_files;i++){
    slicei = sliceinfo+i;
    if(slicei->loaded==1){
      nsliceloaded++;
    }
  }

  nvsliceloaded=0;
  for(i=0;i<nvslice;i++){
    vslicei = vsliceinfo+i;
    if(vslicei->loaded==1){
      nvsliceloaded++;
    }
  }

  nisoloaded=0;
  for(i=0;i<niso_files;i++){
    isoi = isoinfo+i;
    if(isoi->loaded==1){
      nisoloaded++;
    }
  }

  npatchloaded=0;
  for(i=0;i<npatch_files;i++){
    patchi = patchinfo+i;
    if(patchi->loaded==1){
      npatchloaded++;
    }
  }

  nsmoke3dloaded=0;
  for(i=0;i<nsmoke3d_files;i++){
    smoke3di = smoke3dinfo+i;
    if(smoke3di->loaded==1){
      nsmoke3dloaded++;
    }
  }

  nplot3dloaded=0;
  for(i=0;i<nplot3d_files;i++){
    plot3di = plot3dinfo+i;
    if(plot3di->loaded==1){
      nplot3dloaded++;
    }
  }
}

/* ------------------ update_showhidebuttons ------------------------ */

extern "C" void update_showhidebuttons(void){

  update_fileload();
  if(Button_PART!=NULL){
    if(npartloaded==0){
      Button_PART->disable();
    }
    else{
      Button_PART->enable();
    }
  }
  
  if(Button_SLICE!=NULL){
    if(nsliceloaded==0){
      Button_SLICE->disable();
    }
    else{
      Button_SLICE->enable();
    }
  }

  if(Button_VSLICE!=NULL){
    if(nvsliceloaded==0){
      Button_VSLICE->disable();
    }
    else{
      Button_VSLICE->enable();
    }
  }

  if(Button_ISO!=NULL){
    if(nisoloaded==0){
      Button_ISO->disable();
    }
    else{
      Button_ISO->enable();
    }
  }

  if(Button_BOUNDARY!=NULL){
    if(npatchloaded==0){
      Button_BOUNDARY->disable();
    }
    else{
      Button_BOUNDARY->enable();
    }
  }

  if(Button_3DSMOKE!=NULL){
    if(nsmoke3dloaded==0){
      Button_3DSMOKE->disable();
    }
    else{
      Button_3DSMOKE->enable();
    }
  }

  if(Button_PLOT3D!=NULL){
    if(nplot3dloaded==0){
      Button_PLOT3D->disable();
    }
    else{
      Button_PLOT3D->enable();
    }
  }

  if(nplot3dloaded==0&&nsmoke3dloaded==0&&nisoloaded==0&&nsliceloaded==0&&npartloaded==0&&npatchloaded==0){
    if(RADIO_showhide!=NULL)RADIO_showhide->disable();
  }
  else{
    if(RADIO_showhide!=NULL)RADIO_showhide->enable();
  }
}

/* ------------------ hide_glui_labels ------------------------ */

extern "C" void hide_glui_labels(void){
  if(glui_labels!=NULL)glui_labels->hide();
  showlabels=0;
  updatemenu=1;
}

/* ------------------ show_glui_labels ------------------------ */

extern "C" void show_glui_labels(void){
  if(glui_labels!=NULL)glui_labels->show();
}

/* ------------------ Labels_CB ------------------------ */

void Labels_CB(int var){
  updatemenu=1;
  switch (var){
  case LABELS_sensorsize:
    if(sensorrelsize<sensorrelsizeMIN){
      sensorrelsize=sensorrelsizeMIN;
      if(SPINNER_sensorrelsize!=NULL){
        SPINNER_sensorrelsize->set_float_val(sensorrelsize);
      }
    }
    break;
  case SAVE_SETTINGS:
    writeini(LOCAL_INI);
    break;
  case LABELS_HMS:
    break;
  case LABELS_BENCHMARK:
    if(showtime==1){
      ResetMenu(1);        // change to external view
      ResetMenu(-103);     // reset time to 0
      TourMenu(-2);        // hide all tours
      FrameRateMenu(1000); // set frame rate to unlimited
      visFramerate=0;
      LabelMenu(3);        // show frame rate label
      benchmark_flag=1;
      printf("\n*** benchmarking started\n");
      printf(  "    SINGLE rather than DOUBLE buffering is used during the benchmark\n");
      printf(  "    for more accurate results.  The \"flashing\" that occurs is normal.\n");
    }
    else{
      printf("\n*** Error: a file needs to be loaded before benchmarking\n");
      printf("           can be performed\n");
    }
    break;
  case LABELS_showall:
    LabelMenu(4);
    break;
  case LABELS_hideall:
    LabelMenu(5);
    break;
  case LABELS_flip:
    background_flip = 1 - background_flip;
    ShadeMenu(1);
    break;
  case LABELS_shade:
    setbw = 1 - setbw;
    ShadeMenu(2);
    break;
  case LABELS_transparent:
    break;
  case LABELS_close:
    hide_glui_labels();
    break;
  case LABELS_fontsize:
    FontMenu(fontindex);
    break;
  case LABELS_evacshow:
    switch (showhide_option){
    case 0:
      EvacShowMenu(3);
      break;
    case 1:
      EvacShowMenu(3);
      if(npartloaded!=0)ParticleShowMenu(HIDEALL_PARTICLE);
      if(nsmoke3dloaded!=0)Smoke3DShowMenu(HIDEALL_SMOKE3D);
      if(nisoloaded!=0)IsoShowMenu(HIDEALL_ISO);
      if(nsliceloaded!=0)ShowHideSliceMenu(HIDEALL_SLICE);
      if(nvsliceloaded!=0)ShowVSliceMenu(HIDEALL_VSLICE);
      if(npatchloaded!=0)ShowPatchMenu(HIDEALL_BOUNDARY);
      break;
    case 2:
      EvacShowMenu(4);
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
    break;
    break;
  case  LABELS_particleshow:
    switch (showhide_option){
    case 0:
      ParticleShowMenu(3);
      break;
    case 1:
      ParticleShowMenu(3);
      if(nevacloaded!=0)EvacShowMenu(HIDEALL_PARTICLE);
      if(nsmoke3dloaded!=0)Smoke3DShowMenu(HIDEALL_SMOKE3D);
      if(nisoloaded!=0)IsoShowMenu(HIDEALL_ISO);
      if(nsliceloaded!=0)ShowHideSliceMenu(HIDEALL_SLICE);
      if(nvsliceloaded!=0)ShowVSliceMenu(HIDEALL_VSLICE);
      if(npatchloaded!=0)ShowPatchMenu(HIDEALL_BOUNDARY);
      break;
    case 2:
      ParticleShowMenu(4);
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
    break;
  case  LABELS_sliceshow:
    switch (showhide_option){
    case 0:
      ShowHideSliceMenu(SHOWALL_SLICE);
      break;
    case 1:
      ShowHideSliceMenu(SHOWALL_SLICE);
      if(nevacloaded!=0)EvacShowMenu(HIDEALL_PARTICLE);
      if(nvsliceloaded!=0)ShowVSliceMenu(HIDEALL_VSLICE);
      if(npatchloaded!=0)ShowPatchMenu(HIDEALL_BOUNDARY);
      if(nsmoke3dloaded!=0)Smoke3DShowMenu(HIDEALL_SMOKE3D);
      if(nisoloaded!=0)IsoShowMenu(HIDEALL_ISO);
      if(npartloaded!=0)ParticleShowMenu(HIDEALL_PARTICLE);
      break;
    case 2:
      ShowHideSliceMenu(HIDEALL_SLICE);
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
    break;
  case  LABELS_vsliceshow:
    switch (showhide_option){
    case 0:
      ShowVSliceMenu(SHOWALL_VSLICE);
      break;
    case 1:
      ShowVSliceMenu(SHOWALL_VSLICE);
      if(nevacloaded!=0)EvacShowMenu(HIDEALL_PARTICLE);
      if(npatchloaded!=0)ShowPatchMenu(HIDEALL_BOUNDARY);
      if(nsmoke3dloaded!=0)Smoke3DShowMenu(HIDEALL_SMOKE3D);
      if(nisoloaded!=0)IsoShowMenu(HIDEALL_ISO);
      if(npartloaded!=0)ParticleShowMenu(HIDEALL_PARTICLE);
      if(nsliceloaded!=0)ShowHideSliceMenu(HIDEALL_SLICE);
      break;
    case 2:
      ShowHideSliceMenu(HIDEALL_SLICE);
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
    break;
  case  LABELS_boundaryshow:
    switch (showhide_option){
    case 0:
      ShowPatchMenu(SHOWALL_BOUNDARY);
      break;
    case 1:
      ShowPatchMenu(SHOWALL_BOUNDARY);
      if(nevacloaded!=0)EvacShowMenu(HIDEALL_PARTICLE);
      if(nsmoke3dloaded!=0)Smoke3DShowMenu(HIDEALL_SMOKE3D);
      if(npartloaded!=0)ParticleShowMenu(HIDEALL_PARTICLE);
      if(nvsliceloaded!=0)ShowVSliceMenu(HIDEALL_VSLICE);
      if(nsliceloaded!=0)ShowHideSliceMenu(HIDEALL_SLICE);
      if(nisoloaded!=0)IsoShowMenu(HIDEALL_ISO);
      break;
    case 2:
      ShowPatchMenu(HIDEALL_BOUNDARY);
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
    break;
  case  LABELS_3dsmokeshow:
    switch (showhide_option){
    case 0:
      Smoke3DShowMenu(SHOWALL_SMOKE3D);
      break;
    case 1:
      Smoke3DShowMenu(SHOWALL_SMOKE3D);
      if(nevacloaded!=0)EvacShowMenu(HIDEALL_PARTICLE);
      if(npatchloaded!=0)ShowPatchMenu(HIDEALL_BOUNDARY);
      if(npartloaded!=0)ParticleShowMenu(HIDEALL_PARTICLE);
      if(nvsliceloaded!=0)ShowVSliceMenu(HIDEALL_VSLICE);
      if(nsliceloaded!=0)ShowHideSliceMenu(HIDEALL_SLICE);
      if(nisoloaded!=0)IsoShowMenu(HIDEALL_ISO);
      break;
    case 2:
      Smoke3DShowMenu(HIDEALL_SMOKE3D);
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
    break;
  case  LABELS_isosurfaceshow:
    switch (showhide_option){
    case 0:
      IsoShowMenu(10001);
      break;
    case 1:
      IsoShowMenu(SHOWALL_ISO);
      if(nevacloaded!=0)EvacShowMenu(HIDEALL_PARTICLE);
      if(nsmoke3dloaded!=0)Smoke3DShowMenu(HIDEALL_SMOKE3D);
      if(npatchloaded!=0)ShowPatchMenu(HIDEALL_BOUNDARY);
      if(npartloaded!=0)ParticleShowMenu(HIDEALL_PARTICLE);
      if(nvsliceloaded!=0)ShowVSliceMenu(HIDEALL_VSLICE);
      if(nsliceloaded!=0)ShowHideSliceMenu(HIDEALL_SLICE);
      break;
    case 2:
      IsoShowMenu(HIDEALL_ISO);
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
    break;
  case  LABELS_PLOT3D:
    switch (showhide_option){
    case 0:
    case 1:
      Plot3DShowMenu(SHOWALL_PLOT3D);
      break;
    case 2:
      Plot3DShowMenu(HIDEALL_PLOT3D);
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
    break;
  case LABELS_label:
    break;
  case FRAME_label:
    visFramelabel=1-visFramelabel;
    LabelMenu(9);
    break;
  case HRR_label:
    visHRRlabel=1-visHRRlabel;
    LabelMenu(16);
    break;
  case HRRPUVCUTOFF_label:
    break;
  case LABELS_ticks:
    if(auto_user_tick_placement==1){
      CHECKBOX_user_ticks_show_x->disable();
      CHECKBOX_user_ticks_show_y->disable();
      CHECKBOX_user_ticks_show_z->disable();
    }
    else{
      CHECKBOX_user_ticks_show_x->enable();
      CHECKBOX_user_ticks_show_y->enable();
      CHECKBOX_user_ticks_show_z->enable();
    }
    break;
  default:
    ASSERT(FFALSE);
  }
}


/* ------------------ set_labels_controls ------------------------ */

extern "C" void set_labels_controls(){

  if(CHECKBOX_vis_user_ticks!=NULL)CHECKBOX_vis_user_ticks->set_int_val(vis_user_ticks);
  if(CHECKBOX_labels_hrrlabel!=NULL)CHECKBOX_labels_hrrlabel->set_int_val(visHRRlabel);
  if(CHECKBOX_labels_hrrcutoff!=NULL)CHECKBOX_labels_hrrcutoff->set_int_val(show_hrrcutoff);
  if(CHECKBOX_labels_title!=NULL)CHECKBOX_labels_title->set_int_val(visTitle0);
  if(CHECKBOX_labels_colorbar!=NULL)CHECKBOX_labels_colorbar->set_int_val(visColorLabels);
  if(CHECKBOX_labels_timebar!=NULL)CHECKBOX_labels_timebar->set_int_val(visTimeLabels);
  if(CHECKBOX_labels_timelabel!=NULL)CHECKBOX_labels_timelabel->set_int_val(visTimeLabels);
  if(CHECKBOX_labels_framelabel!=NULL)CHECKBOX_labels_framelabel->set_int_val(visFramelabel);
  if(CHECKBOX_labels_ticks!=NULL)CHECKBOX_labels_ticks->set_int_val(visTicks);
  if(CHECKBOX_labels_axis!=NULL)CHECKBOX_labels_axis->set_int_val(visaxislabels);
  if(CHECKBOX_labels_framerate!=NULL)CHECKBOX_labels_framerate->set_int_val(visFramerate);
  if(CHECKBOX_labels_average!=NULL)CHECKBOX_labels_average->set_int_val(vis_slice_average);
#ifdef pp_memstatus
  if(CHECKBOX_labels_availmemory!=NULL)CHECKBOX_labels_availmemory->set_int_val(visAvailmemory);
#endif
  if(CHECKBOX_labels_labels!=NULL)CHECKBOX_labels_labels->set_int_val(visLabels);

  if(CHECKBOX_labels_flip!=NULL)CHECKBOX_labels_flip->set_int_val(background_flip);
  if(CHECKBOX_labels_shade!=NULL)CHECKBOX_labels_shade->set_int_val(setbw);
  if(RADIO_fontsize!=NULL)RADIO_fontsize->set_int_val(fontindex);
  if(CHECKBOX_labels_hms!=NULL)CHECKBOX_labels_hms->set_int_val(vishmsTimelabel);
  if(CHECKBOX_labels_gridloc!=NULL)CHECKBOX_labels_gridloc->set_int_val(visgridloc);

}

