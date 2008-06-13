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

// svn revision character string2
extern "C" char glui_labels_revision[]="$Revision$";

extern "C" void Labels_CB(int value);
int nevacloaded,nplot3dloaded,nsmoke3dloaded,nisoloaded,nsliceloaded,nvsliceloaded,npartloaded,npatchloaded;

GLUI *glui_labels=NULL;

GLUI_Panel *panel_transparency=NULL;
GLUI_Spinner *SPINNER_labels_transparency=NULL;
GLUI_Checkbox *CHECKBOX_labels_colorbar=NULL;
GLUI_Checkbox *CHECKBOX_labels_timebar=NULL;
GLUI_Checkbox *CHECKBOX_labels_ticks=NULL;
GLUI_Checkbox *CHECKBOX_labels_title=NULL;
GLUI_Checkbox *CHECKBOX_labels_axis=NULL;
GLUI_Checkbox *CHECKBOX_labels_hms=NULL;

GLUI_Checkbox *CHECKBOX_labels_framerate=NULL;
GLUI_Checkbox *CHECKBOX_labels_timelabel=NULL;
GLUI_Checkbox *CHECKBOX_labels_framelabel=NULL;
GLUI_Checkbox *CHECKBOX_labels_hrrlabel=NULL;
GLUI_Checkbox *CHECKBOX_labels_hrrcutoff=NULL;
GLUI_Checkbox *CHECKBOX_labels_availmemory=NULL;
GLUI_Checkbox *CHECKBOX_labels_labels=NULL;
GLUI_Checkbox *CHECKBOX_labels_gridloc=NULL;
GLUI_Checkbox *CHECKBOX_labels_average=NULL;

GLUI_Checkbox *CHECKBOX_labels_flip=NULL;
GLUI_Checkbox *CHECKBOX_labels_shade=NULL;
GLUI_Checkbox *CHECKBOX_labels_transparent=NULL;
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

  CHECKBOX_labels_colorbar=glui_labels->add_checkbox("Color Bar",&visColorLabels,LABELS_label,Labels_CB);
  CHECKBOX_labels_timebar=glui_labels->add_checkbox("Time Bar",&visTimeLabels,LABELS_label,Labels_CB);
  CHECKBOX_labels_timelabel=glui_labels->add_checkbox("Time Label",&visTimelabel,LABELS_label,Labels_CB);
  CHECKBOX_labels_framelabel=glui_labels->add_checkbox("Frame Label",&visFramelabel,FRAME_label,Labels_CB);
  CHECKBOX_labels_hrrlabel=glui_labels->add_checkbox("HRR Label",&visHRRlabel,HRR_label,Labels_CB);
  CHECKBOX_labels_hrrcutoff=glui_labels->add_checkbox("HRRPUV cutoff",&show_hrrcutoff,HRRPUVCUTOFF_label,Labels_CB);
  CHECKBOX_labels_ticks=glui_labels->add_checkbox("Ticks",&visTicks,LABELS_label,Labels_CB);
  if(ntotal_blockages>0||isZoneFireModel==0){
    CHECKBOX_labels_gridloc=glui_labels->add_checkbox("Grid Loc",&visgridloc,LABELS_label,Labels_CB);
  }
  if(nslice>0)CHECKBOX_labels_average=glui_labels->add_checkbox("Average",&vis_slice_average,LABELS_label,Labels_CB);
  glui_labels->add_column(false);


  CHECKBOX_labels_title=glui_labels->add_checkbox("Title",&visTitle0,LABELS_label,Labels_CB);
  CHECKBOX_labels_axis=glui_labels->add_checkbox("Axis",&visaxislabels,LABELS_label,Labels_CB);
  CHECKBOX_labels_framerate=glui_labels->add_checkbox("Frame Rate",&visFramerate,LABELS_label,Labels_CB);
#ifdef pp_memstatus
  CHECKBOX_labels_availmemory=glui_labels->add_checkbox("Memory Load",&visAvailmemory,LABELS_label,Labels_CB);
#endif
  CHECKBOX_labels_labels=glui_labels->add_checkbox("Text labels",&visLabels,LABELS_label,Labels_CB);
  glui_labels->add_button("Show All",LABELS_showall,Labels_CB);
  glui_labels->add_button("Hide All",LABELS_hideall,Labels_CB);


  glui_labels->add_column(true);

  CHECKBOX_labels_flip=glui_labels->add_checkbox("Flip Background",&flip,LABELS_flip,Labels_CB);
  CHECKBOX_labels_shade=glui_labels->add_checkbox("Color",&setbw,LABELS_shade,Labels_CB);

  if(nface_transparent>0){
    panel_transparency = glui_labels->add_panel("Transparency");
    CHECKBOX_labels_transparent=glui_labels->add_checkbox_to_panel(panel_transparency,"Activate",&transparentflag,LABELS_transparent,Labels_CB);
    CHECKBOX_labels_transparent_override=glui_labels->add_checkbox_to_panel(panel_transparency,"Override",&transparency_override,LABELS_transparent,Labels_CB);
    SPINNER_labels_transparency=glui_labels->add_spinner_to_panel(panel_transparency,"Override value",GLUI_SPINNER_FLOAT,&transparency_level,LABELS_transparent,Labels_CB);
    SPINNER_labels_transparency->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);
  }
  else{
    CHECKBOX_labels_transparent=glui_labels->add_checkbox("Transparency",&transparentflag,LABELS_transparent,Labels_CB);
  }
  Labels_CB(LABELS_transparent);

  CHECKBOX_labels_hms=glui_labels->add_checkbox("hms time label",&vishmsTimelabel,LABELS_HMS,Labels_CB);


  RADIO_fontsize = glui_labels->add_radiogroup(&fontindex,LABELS_fontsize,Labels_CB);
  glui_labels->add_radiobutton_to_group(RADIO_fontsize,"small font");
  glui_labels->add_radiobutton_to_group(RADIO_fontsize,"large font");

  if((npartinfo>0)||nslice>0||nvslice>0||niso>0||npatch_files||nsmoke3d>0||nplot3d>0){
    glui_labels->add_column(true);

    panel_showhide = glui_labels->add_panel("Show/Hide Loaded Files");

    RADIO_showhide = glui_labels->add_radiogroup_to_panel(panel_showhide,&showhide_option);
    glui_labels->add_radiobutton_to_group(RADIO_showhide,"Show");
    glui_labels->add_radiobutton_to_group(RADIO_showhide,"Show Only");
    glui_labels->add_radiobutton_to_group(RADIO_showhide,"Hide");

    glui_labels->add_column_to_panel(panel_showhide,false);

    if(nevac>0){}
    if(npartinfo>0&&nevac!=npartinfo)Button_PART=glui_labels->add_button_to_panel(panel_showhide,"Particle",LABELS_particleshow,Labels_CB);
    if(nevac>0)Button_EVAC=glui_labels->add_button_to_panel(panel_showhide,"Evacuation",LABELS_evacshow,Labels_CB);
    if(nslice>0)Button_SLICE=glui_labels->add_button_to_panel(panel_showhide,"Slice",LABELS_sliceshow,Labels_CB);
    if(nvslice>0)Button_VSLICE=glui_labels->add_button_to_panel(panel_showhide,"Vector",LABELS_vsliceshow,Labels_CB);
    if(niso>0)Button_ISO=glui_labels->add_button_to_panel(panel_showhide,"Isosurface",LABELS_isosurfaceshow,Labels_CB);
    if(npatch_files>0)Button_BOUNDARY=glui_labels->add_button_to_panel(panel_showhide,"Boundary",LABELS_boundaryshow,Labels_CB);
    if(nsmoke3d>0)Button_3DSMOKE=glui_labels->add_button_to_panel(panel_showhide,"3D Smoke",LABELS_3dsmokeshow,Labels_CB);
    if(nplot3d>0)Button_PLOT3D=glui_labels->add_button_to_panel(panel_showhide,"Plot3D",LABELS_PLOT3D,Labels_CB);

    update_showhidebuttons();
  }

  glui_labels->add_column(false);

  Button_BENCHMARK=glui_labels->add_button("Benchmark",LABELS_BENCHMARK,Labels_CB);

  glui_labels->add_button("Save Settings",SAVE_SETTINGS,Labels_CB);

  glui_labels->add_button("Close",LABELS_close,Labels_CB);

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
  for(i=0;i<npartinfo;i++){
    parti = partinfo+i;
    if(parti->loaded==1&&parti->evac==0){
      npartloaded++;
    }
    if(parti->loaded==1&&parti->evac==1){
      nevacloaded++;
    }
  }

  nsliceloaded=0;
  for(i=0;i<nslice;i++){
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
  for(i=0;i<niso;i++){
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
  for(i=0;i<nsmoke3d;i++){
    smoke3di = smoke3dinfo+i;
    if(smoke3di->loaded==1){
      nsmoke3dloaded++;
    }
  }

  nplot3dloaded=0;
  for(i=0;i<nplot3d;i++){
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
    flip = 1 - flip;
    ShadeMenu(1);
    break;
  case LABELS_shade:
    setbw = 1 - setbw;
    ShadeMenu(2);
    break;
  case LABELS_transparent:
    if(CHECKBOX_labels_transparent_override!=NULL){
      if(transparentflag==1)CHECKBOX_labels_transparent_override->enable();
      if(transparentflag==0)CHECKBOX_labels_transparent_override->disable();
    }
    if(SPINNER_labels_transparency!=NULL){
      if(transparentflag==1&&transparency_override==1){
        SPINNER_labels_transparency->enable();
      }
      else{
        SPINNER_labels_transparency->disable();
      }
    }
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
  default:
    ASSERT(FFALSE);
  }
}


/* ------------------ set_labels_controls ------------------------ */

extern "C" void set_labels_controls(){

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

  if(CHECKBOX_labels_flip!=NULL)CHECKBOX_labels_flip->set_int_val(flip);
  if(CHECKBOX_labels_transparent!=NULL)CHECKBOX_labels_transparent->set_int_val(transparentflag);
  if(CHECKBOX_labels_shade!=NULL)CHECKBOX_labels_shade->set_int_val(setbw);
  if(RADIO_fontsize!=NULL)RADIO_fontsize->set_int_val(fontindex);
  if(CHECKBOX_labels_hms!=NULL)CHECKBOX_labels_hms->set_int_val(vishmsTimelabel);
  if(CHECKBOX_labels_gridloc!=NULL)CHECKBOX_labels_gridloc->set_int_val(visgridloc);

}

