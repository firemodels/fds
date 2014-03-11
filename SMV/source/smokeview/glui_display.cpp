// $Date$ 
// $Revision$
// $Author$

#define CPP
#include "options.h"

// svn revision character string2
extern "C" char glui_labels_revision[];
char glui_labels_revision[]="$Revision$";

extern "C" void FileShow_CB(int var);
extern "C" void ShowHideMenu(int val);
extern "C" void colorbar_global2local(void);
extern "C" void Volume_CB(int var);

void Text_Labels_CB(int var);

#include <stdio.h>
#include <string.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <math.h>

#include "smokeviewvars.h"

int nevacloaded,nplot3dloaded,nsmoke3dloaded,nisoloaded,nsliceloaded,nvsliceloaded,npartloaded,npatchloaded;

GLUI *glui_labels=NULL;

GLUI_EditText *EDIT_LB_label_string=NULL;

GLUI_Spinner *SPINNER_colorband=NULL;
GLUI_Spinner *SPINNER_labels_transparency_data=NULL;
#ifdef pp_BETA
GLUI_Spinner *SPINNER_cullgeom_portsize=NULL;
#endif

GLUI_Listbox *LIST_colorbar2=NULL;
GLUI_Listbox *LIST_LB_labels=NULL;

GLUI_Spinner *SPINNER_down_red=NULL,*SPINNER_down_green=NULL,*SPINNER_down_blue=NULL;
GLUI_Spinner *SPINNER_up_red=NULL,*SPINNER_up_green=NULL,*SPINNER_up_blue=NULL;
GLUI_Spinner *SPINNER_LB_time_start=NULL;
GLUI_Spinner *SPINNER_LB_time_stop=NULL;
GLUI_Spinner *SPINNER_LB_red=NULL;
GLUI_Spinner *SPINNER_LB_green=NULL;
GLUI_Spinner *SPINNER_LB_blue=NULL;
GLUI_Spinner *SPINNER_LB_x=NULL;
GLUI_Spinner *SPINNER_LB_y=NULL;
GLUI_Spinner *SPINNER_LB_z=NULL;
GLUI_Spinner *SPINNER_tick_xmin=NULL;
GLUI_Spinner *SPINNER_tick_ymin=NULL;
GLUI_Spinner *SPINNER_tick_zmin=NULL;
GLUI_Spinner *SPINNER_tick_xmax=NULL;
GLUI_Spinner *SPINNER_tick_ymax=NULL;
GLUI_Spinner *SPINNER_tick_zmax=NULL;
GLUI_Spinner *SPINNER_gridlinewidth=NULL;
GLUI_Spinner *SPINNER_linewidth=NULL;
GLUI_Spinner *SPINNER_tick_x0=NULL;
GLUI_Spinner *SPINNER_tick_y0=NULL;
GLUI_Spinner *SPINNER_tick_z0=NULL;
GLUI_Spinner *SPINNER_tick_dx0=NULL;
GLUI_Spinner *SPINNER_tick_dy0=NULL;
GLUI_Spinner *SPINNER_tick_dz0=NULL;
GLUI_Spinner *SPINNER_labels_transparency_face=NULL;
GLUI_Spinner *SPINNER_subtick=NULL;
GLUI_Spinner *SPINNER_scaled_font2d_height=NULL;
GLUI_Spinner *SPINNER_scaled_font3d_height=NULL;
GLUI_Spinner *SPINNER_scaled_font2d_height2width=NULL;
GLUI_Spinner *SPINNER_scaled_font3d_height2width=NULL;
GLUI_Spinner *SPINNER_scaled_font3d_thickness=NULL;
GLUI_Spinner *SPINNER_scaled_font2d_thickness=NULL;

GLUI_Checkbox *CHECKBOX_labels_meshlabel=NULL;
GLUI_Checkbox *CHECKBOX_labels_version=NULL;
GLUI_Checkbox *CHECKBOX_vis_user_ticks=NULL;
GLUI_Checkbox *CHECKBOX_vis_user_ticks2=NULL;
GLUI_Checkbox *CHECKBOX_show_extreme_mindata=NULL;
GLUI_Checkbox *CHECKBOX_show_extreme_maxdata=NULL;
GLUI_Checkbox *CHECKBOX_colorbarflip=NULL;
#ifdef pp_BETA
GLUI_Checkbox *CHECKBOX_cullgeom=NULL;
#endif
GLUI_Checkbox *CHECKBOX_axislabels_smooth=NULL, *CHECKBOX_transparentflag=NULL;
GLUI_Checkbox *CHECKBOX_LB_visLabels=NULL;
GLUI_Checkbox *CHECKBOX_LB_label_use_foreground=NULL;
GLUI_Checkbox *CHECKBOX_LB_label_show_always=NULL;
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
GLUI_Checkbox *CHECKBOX_user_ticks_show_x=NULL;
GLUI_Checkbox *CHECKBOX_user_ticks_show_y=NULL;
GLUI_Checkbox *CHECKBOX_user_ticks_show_z=NULL;
GLUI_Checkbox *CHECKBOX_tick_auto=NULL;
GLUI_Checkbox *CHECKBOX_label_1=NULL;
GLUI_Checkbox *CHECKBOX_label_2=NULL;
GLUI_Checkbox *CHECKBOX_label_3=NULL;
GLUI_Checkbox *CHECKBOX_labels_flip=NULL;
GLUI_Checkbox *CHECKBOX_labels_shade=NULL;
GLUI_Checkbox *CHECKBOX_labels_transparent_override=NULL;

GLUI_Rollout *ROLLOUT_scene=NULL;
GLUI_Rollout *ROLLOUT_font=NULL;
GLUI_Rollout *ROLLOUT_user_labels=NULL;
GLUI_Rollout *ROLLOUT_user_tick=NULL;
GLUI_Rollout *ROLLOUT_label1=NULL;

GLUI_Panel *PANEL_extreme=NULL,*PANEL_cb8=NULL,*PANEL_cb7=NULL;
GLUI_Panel *PANEL_extreme_min=NULL, *PANEL_extreme_max=NULL;
GLUI_Panel *PANEL_extreme2=NULL;
GLUI_Panel *PANEL_cb11=NULL;
GLUI_Panel *PANEL_contours=NULL;
GLUI_Panel *PANEL_gen1=NULL, *PANEL_gen2=NULL, *PANEL_gen3=NULL;
GLUI_Panel *PANEL_LB_panel1=NULL, *PANEL_LB_panel2=NULL, *PANEL_LB_panel3=NULL;
GLUI_Panel *PANEL_LB_panel4=NULL, *PANEL_LB_panel5=NULL, *PANEL_LB_panel6=NULL;
GLUI_Panel *PANEL_LB_color=NULL, *PANEL_LB_time=NULL;
GLUI_Panel *PANEL_LB_position=NULL;
GLUI_Panel *PANEL_label2=NULL;
GLUI_Panel *PANEL_tick1;
GLUI_Panel *PANEL_tick1a;
GLUI_Panel *PANEL_tick1b;
GLUI_Panel *PANEL_tick2;
GLUI_Panel *PANEL_transparency=NULL;
GLUI_Panel *PANEL_showhide=NULL;
GLUI_Panel *PANEL_font2d=NULL;
GLUI_Panel *PANEL_font3d=NULL;

GLUI_RadioGroup *RADIO2_plot3d_display=NULL;
GLUI_RadioGroup *RADIO_fontsize=NULL,*RADIO_showhide=NULL;
GLUI_RadioButton *RADIOBUTTON_label_1a=NULL;
GLUI_RadioButton *RADIOBUTTON_label_1b=NULL;
GLUI_RadioButton *RADIOBUTTON_label_1c=NULL;

GLUI_Button *BUTTON_LB_label_previous=NULL;
GLUI_Button *BUTTON_LB_label_next=NULL;
GLUI_Button *BUTTON_LB_label_update=NULL;
GLUI_Button *BUTTON_LB_label_add=NULL;
GLUI_Button *BUTTON_LB_label_delete=NULL;
GLUI_Button *BUTTON_LB_label_set=NULL;
GLUI_Button *BUTTON_EVAC=NULL;
GLUI_Button *BUTTON_PART=NULL;
GLUI_Button *BUTTON_SLICE=NULL;
GLUI_Button *BUTTON_VSLICE=NULL;
GLUI_Button *BUTTON_PLOT3D=NULL;
GLUI_Button *BUTTON_3DSMOKE=NULL;
GLUI_Button *BUTTON_BOUNDARY=NULL;
GLUI_Button *BUTTON_ISO=NULL;
GLUI_Button *BUTTON_label_1=NULL;
GLUI_Button *BUTTON_label_2=NULL;
GLUI_Button *BUTTON_label_3=NULL;
GLUI_Button *BUTTON_label_4=NULL;

#ifdef pp_GEOMTEST
GLUI_Rollout *ROLLOUT_geomtest=NULL;
GLUI_Panel *PANEL_geom1=NULL;
GLUI_Panel *PANEL_geom1a=NULL;
GLUI_Panel *PANEL_geom1b=NULL;
GLUI_Panel *PANEL_geom1c=NULL;
GLUI_Panel *PANEL_geom1d=NULL;
GLUI_Panel *PANEL_geom2=NULL;
GLUI_Panel *PANEL_geom2a=NULL;
GLUI_Panel *PANEL_geom2b=NULL;
GLUI_Panel *PANEL_geom2c=NULL;
GLUI_Spinner *SPINNER_box_bounds[6];
GLUI_Spinner *SPINNER_box_translate[3];
GLUI_Spinner *SPINNER_tetra_vertices[12];
#endif

#define VOL_BOXTRANSLATE 0

#define COLORBAR_EXTREME_RGB 15
#define COLORBAR_EXTREME 16
#define FLIP 19

#define LB_LIST 0
#define LB_ADD 1
#define LB_DELETE 2
#define LB_RGB 3
#define LB_XYZ 4
#define LB_STARTSTOP 5
#define LB_SHOWALWAYS 6
#define LB_FOREGROUND 7
#define LB_UPDATE 8
#define LB_PREVIOUS 9
#define LB_NEXT 10
#define LB_VISLABELS 11

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
#define LABELS_drawface 24
#define LABELS_hide_overlaps 25
#define LABELS_version 26
#define LABELS_meshlabel 27
#define LABELS_usertick 28
#define LABELS_usertick2 29

#define FILESHOW_particle    10
#define FILESHOW_slice       11
#define FILESHOW_vslice      12
#define FILESHOW_boundary    13
#define FILESHOW_3dsmoke     14
#define FILESHOW_isosurface  15
#define FILESHOW_evac 19
#define FILESHOW_plot3d 16
#define LABELS_HMS 18
#define SAVE_SETTINGS 99

#define COLORBAR_SMOOTH 113
#define COLORBAND 115

#define COLORBAR_LIST2 112
#define DATA_transparent 26
#define TRANSPARENTLEVEL 110

#define UPDATEPLOT 10
extern "C" void PLOT3D_CB(int var);
extern "C" void Extreme_CB(int var);


int cb_up_rgb[3],cb_down_rgb[3];


/* ------------------ glui_labels_rename ------------------------ */

extern "C" void update_glui_label_text(void){
  if(LABEL_Get_Nuserlabels()>0){
    labeldata *gl;

    gl=&LABEL_local;

    LIST_LB_labels->set_int_val(gl->glui_id);
    EDIT_LB_label_string->set_text(gl->name);
    SPINNER_LB_x->set_float_val(gl->xyz[0]);
    SPINNER_LB_y->set_float_val(gl->xyz[1]);
    SPINNER_LB_z->set_float_val(gl->xyz[2]);
    SPINNER_LB_time_start->set_float_val(gl->tstart_stop[0]);
    SPINNER_LB_time_stop->set_float_val(gl->tstart_stop[1]);
    CHECKBOX_LB_label_show_always->set_int_val(gl->show_always);

    SPINNER_LB_red->set_int_val(gl->rgb[0]);
    SPINNER_LB_green->set_int_val(gl->rgb[1]);
    SPINNER_LB_blue->set_int_val(gl->rgb[2]);
    CHECKBOX_LB_label_use_foreground->set_int_val(gl->useforegroundcolor);

    CHECKBOX_LB_visLabels->enable();
    LIST_LB_labels->enable();
    EDIT_LB_label_string->enable();
    SPINNER_LB_x->enable();
    SPINNER_LB_y->enable();
    SPINNER_LB_z->enable();
    SPINNER_LB_time_start->enable();
    SPINNER_LB_time_stop->enable();
    SPINNER_LB_red->enable();
    SPINNER_LB_green->enable();
    SPINNER_LB_blue->enable();
    CHECKBOX_LB_label_show_always->enable();
  }
  else{
    CHECKBOX_LB_visLabels->disable();
    LIST_LB_labels->disable();
    EDIT_LB_label_string->disable();
    SPINNER_LB_x->disable();
    SPINNER_LB_y->disable();
    SPINNER_LB_z->disable();
    SPINNER_LB_time_start->disable();
    SPINNER_LB_time_stop->disable();
    SPINNER_LB_red->disable();
    SPINNER_LB_green->disable();
    SPINNER_LB_blue->disable();
    CHECKBOX_LB_label_show_always->disable();
  }
}

  /* ------------------ glui_labels_rename ------------------------ */

extern "C" void glui_update_fontindex(void){
  if(RADIO_fontsize!=NULL){
    if(fontindex==SCALED_FONT){
      SPINNER_scaled_font2d_height->enable();
      SPINNER_scaled_font3d_height->enable();
      SPINNER_scaled_font2d_height2width->enable();
      SPINNER_scaled_font3d_height2width->enable();
      SPINNER_scaled_font2d_thickness->enable();
      SPINNER_scaled_font3d_thickness->enable();
    }
    else{
      SPINNER_scaled_font2d_height->disable();
      SPINNER_scaled_font3d_height->disable();
      SPINNER_scaled_font2d_height2width->disable();
      SPINNER_scaled_font3d_height2width->disable();
      SPINNER_scaled_font2d_thickness->disable();
      SPINNER_scaled_font3d_thickness->disable();
    }
  }
}
/* ------------------ glui_labels_rename ------------------------ */

extern "C" void glui_labels_rename(void){

  ROLLOUT_label1->set_name(_("General Settings"));
  CHECKBOX_labels_colorbar->set_name(_("Colorbar"));
  CHECKBOX_labels_timebar->set_name(_("Time bar"));
  CHECKBOX_labels_timelabel->set_name(_("Time label"));
  CHECKBOX_labels_framelabel->set_name(_("Frame label"));
  CHECKBOX_labels_hrrlabel->set_name(_("HRR label"));
  CHECKBOX_labels_hrrcutoff->set_name(_("HRRPUV cutoff"));
  CHECKBOX_labels_ticks->set_name(_("FDS ticks"));
  if(ntotal_blockages>0||isZoneFireModel==0){
    CHECKBOX_labels_gridloc->set_name(_("Grid loc"));
  }
  if(nsliceinfo>0)CHECKBOX_labels_average->set_name(_("Average"));

  CHECKBOX_labels_title->set_name(_("Title"));
  CHECKBOX_labels_axis->set_name(_("Axis"));
  CHECKBOX_labels_framerate->set_name(_("Frame rate"));
#ifdef pp_memstatus
  CHECKBOX_labels_availmemory->set_name(_("Memory load"));
#endif
  CHECKBOX_labels_labels->set_name(_("Text labels"));
  CHECKBOX_labels_meshlabel->set_name(_("Mesh label"));
  
  CHECKBOX_label_1->set_name(_("Fast blockage drawing"));
  CHECKBOX_label_2->set_name(_("Sort transparent faces"));
  CHECKBOX_label_3->set_name(_("Hide overlaps"));
  BUTTON_label_1->set_name(_("Show all"));
  BUTTON_label_2->set_name(_("Hide all"));

  CHECKBOX_labels_flip->set_name(_("Flip background"));
  CHECKBOX_labels_shade->set_name(_("Shades of grey"));

  if(nface_transparent>0){
    PANEL_transparency->set_name(_("Geometry transparency"));
    CHECKBOX_labels_transparent_override->set_name(_("Use level:"));
  }


  CHECKBOX_labels_hms->set_name(_("hms time label"));

  RADIOBUTTON_label_1a->set_name(_("small"));
  RADIOBUTTON_label_1b->set_name(_("large"));
  RADIOBUTTON_label_1c->set_name(_("scaled"));


  ROLLOUT_user_tick->set_name("User tick settings");

  PANEL_tick1->set_name(_("Display settings"));

  CHECKBOX_vis_user_ticks->set_name(_("User ticks"));
  SPINNER_subtick->set_name(_("sub-intervals")); 
  CHECKBOX_tick_auto->set_name(_("Auto place (2D)"));

  PANEL_tick2->set_name(_("Parameters"));
  SPINNER_tick_x0->set_name(_("origin"));
  SPINNER_tick_xmin->set_name(_("Min"));
  SPINNER_tick_xmax->set_name(_("Max"));
  SPINNER_tick_dx0->set_name(_("Step"));

  if((npartinfo>0)||nsliceinfo>0||nvsliceinfo>0||nisoinfo>0||npatchinfo||nsmoke3dinfo>0||nplot3dinfo>0){
    PANEL_showhide->set_name(_("Show/Hide Loaded Files"));
  }

//    glui_labels->add_radiobutton_to_group(RADIO_showhide,_("Show"));
//    glui_labels->add_radiobutton_to_group(RADIO_showhide,_("Show Only"));
//    glui_labels->add_radiobutton_to_group(RADIO_showhide,_("Hide"));

  BUTTON_label_3->set_name(_("Save settings"));
  BUTTON_label_4->set_name(_("Close"));
}


/* ------------------ glui_labels_setup ------------------------ */

extern "C" void glui_labels_setup(int main_window){
  labeldata *gl;
  int init_status=0;

  update_glui_labels=0;
  if(glui_labels!=NULL){
    glui_labels->close();
    glui_labels=NULL;
  }
  glui_labels = GLUI_Master.create_glui("Display settings",0,0,0);
  if(showdisplay_dialog==0)glui_labels->hide();

  // -------------- General Settings -------------------

  ROLLOUT_label1 = glui_labels->add_rollout(_("General"),true);
  PANEL_gen1=glui_labels->add_panel_to_panel(ROLLOUT_label1,"",GLUI_PANEL_NONE);
  CHECKBOX_labels_colorbar=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("Colorbar"),&visColorbarLabels,LABELS_label,Labels_CB);
  CHECKBOX_labels_timebar=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("Time bar"),&visTimeLabels,LABELS_label,Labels_CB);
  CHECKBOX_labels_timelabel=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("Time label"),&visTimelabel,LABELS_label,Labels_CB);
  CHECKBOX_labels_framelabel=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("Frame label"),&visFramelabel,FRAME_label,Labels_CB);
  CHECKBOX_labels_hrrlabel=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("HRR label"),&visHRRlabel,HRR_label,Labels_CB);
  CHECKBOX_labels_hrrcutoff=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("HRRPUV cutoff"),&show_hrrcutoff,HRRPUVCUTOFF_label,Labels_CB);
  CHECKBOX_labels_ticks=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("FDS ticks"),&visTicks,LABELS_label,Labels_CB);
  if(ntickinfo>0){
    CHECKBOX_labels_ticks->enable();
  }
  else{
    CHECKBOX_labels_ticks->disable();
    visTicks=0;
    CHECKBOX_labels_ticks->set_int_val(visTicks);
  }
  CHECKBOX_vis_user_ticks2=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("User ticks"),&vis_user_ticks,LABELS_usertick2,Labels_CB);
  CHECKBOX_labels_version=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("Version info"),&gversion,LABELS_version,Labels_CB);

  glui_labels->add_column_to_panel(PANEL_gen1,false);
  if(ntotal_blockages>0||isZoneFireModel==0){
    CHECKBOX_labels_gridloc=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("Grid loc"),&visgridloc,LABELS_label,Labels_CB);
  }
  if(nsliceinfo>0)CHECKBOX_labels_average=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("Average"),&vis_slice_average,LABELS_label,Labels_CB);

  CHECKBOX_labels_title=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("Title"),&visTitle,LABELS_label,Labels_CB);
  CHECKBOX_labels_axis=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("Axis"),&visaxislabels,LABELS_label,Labels_CB);
  CHECKBOX_labels_framerate=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("Frame rate"),&visFramerate,LABELS_label,Labels_CB);
#ifdef pp_memstatus
  CHECKBOX_labels_availmemory=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("Memory load"),&visAvailmemory,LABELS_label,Labels_CB);
#endif
  CHECKBOX_labels_labels=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("Text labels"),&visLabels,LABELS_label,Labels_CB);
  CHECKBOX_labels_meshlabel=glui_labels->add_checkbox_to_panel(PANEL_gen1,_("Mesh label"),&visBlocklabel,LABELS_meshlabel,Labels_CB);

  PANEL_gen2=glui_labels->add_panel_to_panel(ROLLOUT_label1,"",GLUI_PANEL_NONE);

  BUTTON_label_1=glui_labels->add_button_to_panel(PANEL_gen2,_("Show all"),LABELS_showall,Labels_CB);
  glui_labels->add_column_to_panel(PANEL_gen2,false);
  BUTTON_label_2=glui_labels->add_button_to_panel(PANEL_gen2,_("Hide all"),LABELS_hideall,Labels_CB);

  glui_labels->add_separator_to_panel(ROLLOUT_label1);

  PANEL_gen3=glui_labels->add_panel_to_panel(ROLLOUT_label1,"",GLUI_PANEL_NONE);

  CHECKBOX_labels_flip=glui_labels->add_checkbox_to_panel(PANEL_gen3,_("Flip background"),&background_flip,LABELS_flip,Labels_CB);
  CHECKBOX_labels_hms=glui_labels->add_checkbox_to_panel(PANEL_gen3,_("hms time label"),&vishmsTimelabel,LABELS_HMS,Labels_CB);
  SPINNER_linewidth=glui_labels->add_spinner_to_panel(PANEL_gen3,"blockage line width",GLUI_SPINNER_FLOAT,&linewidth);
  SPINNER_linewidth->set_float_limits(1.0,10.0,GLUI_LIMIT_CLAMP);
  SPINNER_gridlinewidth=glui_labels->add_spinner_to_panel(PANEL_gen3,"grid line width",GLUI_SPINNER_FLOAT,&gridlinewidth);
  SPINNER_gridlinewidth->set_float_limits(1.0,10.0,GLUI_LIMIT_CLAMP);

  glui_labels->add_column_to_panel(PANEL_gen3,false);

  CHECKBOX_label_1=glui_labels->add_checkbox_to_panel(PANEL_gen3,_("Fast blockage drawing"),&use_new_drawface,LABELS_drawface,Labels_CB);
  CHECKBOX_label_2=glui_labels->add_checkbox_to_panel(PANEL_gen3,_("Sort transparent faces"),&sort_transparent_faces,LABELS_drawface,Labels_CB);
  CHECKBOX_label_3=glui_labels->add_checkbox_to_panel(PANEL_gen3,_("Hide overlaps"),&hide_overlaps,LABELS_hide_overlaps,Labels_CB);
 
  if(nface_transparent>0){
    glui_labels->add_column_to_panel(PANEL_gen1,true);
    PANEL_transparency = glui_labels->add_panel_to_panel(PANEL_gen3,_("Geometry transparency"));
    CHECKBOX_labels_transparent_override=glui_labels->add_checkbox_to_panel(PANEL_transparency,_("Use level:"),&use_transparency_geom,LABELS_transparent,Labels_CB);
    SPINNER_labels_transparency_face=glui_labels->add_spinner_to_panel(PANEL_transparency,"",GLUI_SPINNER_FLOAT,&transparency_geom,LABELS_transparent,Labels_CB);
    SPINNER_labels_transparency_face->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);
    Labels_CB(LABELS_transparent);
  }

  // -------------- Data coloring -------------------

  ROLLOUT_scene = glui_labels->add_rollout("Data coloring",false);

  if(ncolorbars>0){
    int i;

    selectedcolorbar_index2=-1;
    LIST_colorbar2=glui_labels->add_listbox_to_panel(ROLLOUT_scene,_("Colorbar:"),&selectedcolorbar_index2,COLORBAR_LIST2,Slice_CB);

    for(i=0;i<ncolorbars;i++){
      colorbardata *cbi;

      cbi = colorbarinfo + i;
      cbi->label_ptr=cbi->label;
      LIST_colorbar2->add_item(i,cbi->label_ptr);
    }
    LIST_colorbar2->set_int_val(colorbartype);
  }

  PANEL_cb11=glui_labels->add_panel_to_panel(ROLLOUT_scene,"",GLUI_PANEL_NONE);

  PANEL_contours = glui_labels->add_panel_to_panel(PANEL_cb11,_("Colorbar type:"));
  RADIO2_plot3d_display=glui_labels->add_radiogroup_to_panel(PANEL_contours,&contour_type,UPDATEPLOT,PLOT3D_CB);
  glui_labels->add_radiobutton_to_group(RADIO2_plot3d_display,_("Continuous"));
  glui_labels->add_radiobutton_to_group(RADIO2_plot3d_display,_("Stepped"));
  glui_labels->add_radiobutton_to_group(RADIO2_plot3d_display,_("Line"));

  SPINNER_colorband=glui_labels->add_spinner_to_panel(PANEL_cb11,"Selection width:",GLUI_SPINNER_INT,&colorband,COLORBAND,Slice_CB);
  SPINNER_colorband->set_int_limits(1,10);

  glui_labels->add_column_to_panel(PANEL_cb11,false);

  CHECKBOX_labels_shade=glui_labels->add_checkbox_to_panel(PANEL_cb11,_("color -> grey"),&setbw,LABELS_shade,Labels_CB);
  CHECKBOX_colorbarflip=glui_labels->add_checkbox_to_panel(PANEL_cb11,_("flip colorbar"),&colorbarflip,FLIP,Labels_CB);
  CHECKBOX_axislabels_smooth=glui_labels->add_checkbox_to_panel(PANEL_cb11,_("Smooth colorbar labels"),&axislabels_smooth,COLORBAR_SMOOTH,Slice_CB);
  CHECKBOX_transparentflag=glui_labels->add_checkbox_to_panel(PANEL_cb11,_("Use transparency"),
    &use_transparency_data,DATA_transparent,Slice_CB);
  SPINNER_labels_transparency_data=glui_labels->add_spinner_to_panel(PANEL_cb11,_("transparency level"),
    GLUI_SPINNER_FLOAT,&transparent_level,TRANSPARENTLEVEL,Slice_CB);
  SPINNER_labels_transparency_data->set_w(0);
  SPINNER_labels_transparency_data->set_float_limits(0.0,1.0,GLUI_LIMIT_CLAMP);


  PANEL_extreme = glui_labels->add_panel_to_panel(ROLLOUT_scene,"",GLUI_PANEL_NONE);

  PANEL_extreme2 = glui_labels->add_panel_to_panel(PANEL_extreme,"Highlight extreme data");

  PANEL_extreme_min = glui_labels->add_panel_to_panel(PANEL_extreme2,"",GLUI_PANEL_NONE);
  CHECKBOX_show_extreme_mindata=glui_labels->add_checkbox_to_panel(PANEL_extreme_min,_("Color below min"),&show_extreme_mindata,COLORBAR_EXTREME,Extreme_CB);
  
  SPINNER_down_red=  glui_labels->add_spinner_to_panel(PANEL_extreme_min,_("red"),  GLUI_SPINNER_INT,cb_down_rgb,COLORBAR_EXTREME_RGB,Extreme_CB);
  SPINNER_down_green=glui_labels->add_spinner_to_panel(PANEL_extreme_min,_("green"),GLUI_SPINNER_INT,cb_down_rgb+1,COLORBAR_EXTREME_RGB,Extreme_CB);
  SPINNER_down_blue= glui_labels->add_spinner_to_panel(PANEL_extreme_min,_("blue"), GLUI_SPINNER_INT,cb_down_rgb+2,COLORBAR_EXTREME_RGB,Extreme_CB);
  SPINNER_down_red->set_int_limits(0,255);
  SPINNER_down_green->set_int_limits(0,255);
  SPINNER_down_blue->set_int_limits(0,255);

  glui_labels->add_column_to_panel(PANEL_extreme2,false);

  PANEL_extreme_max = glui_labels->add_panel_to_panel(PANEL_extreme2,"",GLUI_PANEL_NONE);

  CHECKBOX_show_extreme_maxdata=glui_labels->add_checkbox_to_panel(PANEL_extreme_max,_("Color above max"),&show_extreme_maxdata,COLORBAR_EXTREME,Extreme_CB);

  SPINNER_up_red=  glui_labels->add_spinner_to_panel(PANEL_extreme_max,_("red"),  GLUI_SPINNER_INT,cb_up_rgb,COLORBAR_EXTREME_RGB,Extreme_CB);
  SPINNER_up_green=glui_labels->add_spinner_to_panel(PANEL_extreme_max,_("green"),GLUI_SPINNER_INT,cb_up_rgb+1,COLORBAR_EXTREME_RGB,Extreme_CB);
  SPINNER_up_blue= glui_labels->add_spinner_to_panel(PANEL_extreme_max,_("blue"), GLUI_SPINNER_INT,cb_up_rgb+2,COLORBAR_EXTREME_RGB,Extreme_CB);
  SPINNER_up_red->set_int_limits(0,255);
  SPINNER_up_green->set_int_limits(0,255);
  SPINNER_up_blue->set_int_limits(0,255);
  colorbar_global2local();

  // -------------- Show/Hide Loaded files -------------------

  if((npartinfo>0)||nsliceinfo>0||nvsliceinfo>0||nisoinfo>0||npatchinfo||nsmoke3dinfo>0||nplot3dinfo>0){
    PANEL_showhide = glui_labels->add_rollout("Show/Hide Loaded Files",false);

    RADIO_showhide = glui_labels->add_radiogroup_to_panel(PANEL_showhide,&showhide_option);
    glui_labels->add_radiobutton_to_group(RADIO_showhide,_("Show"));
    glui_labels->add_radiobutton_to_group(RADIO_showhide,_("Show Only"));
    glui_labels->add_radiobutton_to_group(RADIO_showhide,_("Hide"));

    glui_labels->add_column_to_panel(PANEL_showhide,false);

    if(nevac>0){}
    if(npartinfo>0&&nevac!=npartinfo)BUTTON_PART=glui_labels->add_button_to_panel(PANEL_showhide,"Particle",FILESHOW_particle,FileShow_CB);
    if(nevac>0)BUTTON_EVAC=glui_labels->add_button_to_panel(PANEL_showhide,"Evacuation",FILESHOW_evac,FileShow_CB);
    if(nsliceinfo>0)BUTTON_SLICE=glui_labels->add_button_to_panel(PANEL_showhide,"Slice",FILESHOW_slice,FileShow_CB);
    if(nvsliceinfo>0)BUTTON_VSLICE=glui_labels->add_button_to_panel(PANEL_showhide,"Vector",FILESHOW_vslice,FileShow_CB);
    if(nisoinfo>0)BUTTON_ISO=glui_labels->add_button_to_panel(PANEL_showhide,"Isosurface",FILESHOW_isosurface,FileShow_CB);
    if(npatchinfo>0)BUTTON_BOUNDARY=glui_labels->add_button_to_panel(PANEL_showhide,"Boundary",FILESHOW_boundary,FileShow_CB);
    if(nsmoke3dinfo>0)BUTTON_3DSMOKE=glui_labels->add_button_to_panel(PANEL_showhide,"3D smoke",FILESHOW_3dsmoke,FileShow_CB);
    if(nplot3dinfo>0)BUTTON_PLOT3D=glui_labels->add_button_to_panel(PANEL_showhide,"Plot3D",FILESHOW_plot3d,FileShow_CB);

    update_showhidebuttons();
  }

  // -------------- Fonts -------------------

  ROLLOUT_font = glui_labels->add_rollout("Fonts",false);
  RADIO_fontsize = glui_labels->add_radiogroup_to_panel(ROLLOUT_font,&fontindex,LABELS_fontsize,Labels_CB);
  RADIOBUTTON_label_1a=glui_labels->add_radiobutton_to_group(RADIO_fontsize,_("small"));
  RADIOBUTTON_label_1b=glui_labels->add_radiobutton_to_group(RADIO_fontsize,_("large"));
  RADIOBUTTON_label_1c=glui_labels->add_radiobutton_to_group(RADIO_fontsize,_("scaled"));

  PANEL_font2d = glui_labels->add_panel_to_panel(ROLLOUT_font,"labels",true);
  SPINNER_scaled_font2d_height=glui_labels->add_spinner_to_panel(PANEL_font2d,"height:",GLUI_SPINNER_INT,&scaled_font2d_height);
  SPINNER_scaled_font2d_height2width=glui_labels->add_spinner_to_panel(PANEL_font2d,"height/width",GLUI_SPINNER_FLOAT,&scaled_font2d_height2width);
  SPINNER_scaled_font2d_height2width->set_float_limits(0.5,1.5);
  SPINNER_scaled_font2d_thickness=glui_labels->add_spinner_to_panel(PANEL_font2d,"thickness:",GLUI_SPINNER_INT,&scaled_font2d_thickness);
  SPINNER_scaled_font2d_thickness->set_int_limits(1,10);

  PANEL_font3d = glui_labels->add_panel_to_panel(ROLLOUT_font,"scene",true);
  SPINNER_scaled_font3d_height=glui_labels->add_spinner_to_panel(PANEL_font3d,"height:",GLUI_SPINNER_INT,&scaled_font3d_height);
  SPINNER_scaled_font3d_height2width=glui_labels->add_spinner_to_panel(PANEL_font3d,"height/width:",GLUI_SPINNER_FLOAT,&scaled_font3d_height2width);
  SPINNER_scaled_font3d_height2width->set_float_limits(0.5,1.5);
  SPINNER_scaled_font3d_thickness=glui_labels->add_spinner_to_panel(PANEL_font3d,"thickness:",GLUI_SPINNER_INT,&scaled_font3d_thickness);
  SPINNER_scaled_font3d_thickness->set_int_limits(1,10);
  glui_update_fontindex();

  // -------------- User tick settings -------------------

  ROLLOUT_user_tick = glui_labels->add_rollout("User ticks",false);


  PANEL_tick1 = glui_labels->add_panel_to_panel(ROLLOUT_user_tick,_("Display"),true);
  PANEL_tick1a = glui_labels->add_panel_to_panel(PANEL_tick1,"",false);

  CHECKBOX_vis_user_ticks=glui_labels->add_checkbox_to_panel(PANEL_tick1a,_("Show user ticks"),&vis_user_ticks,LABELS_usertick,Labels_CB);
  glui_labels->add_column_to_panel(PANEL_tick1a,false);
  SPINNER_subtick=glui_labels->add_spinner_to_panel(PANEL_tick1a,_("sub-intervals"),GLUI_SPINNER_INT,&user_tick_sub); 
  SPINNER_subtick->set_int_limits(1,10,GLUI_LIMIT_CLAMP);

  PANEL_tick1b = glui_labels->add_panel_to_panel(PANEL_tick1,"",false);
  CHECKBOX_tick_auto=glui_labels->add_checkbox_to_panel(PANEL_tick1b,_("Auto place (2D)"),&auto_user_tick_placement,LABELS_ticks,Labels_CB);
  glui_labels->add_column_to_panel(PANEL_tick1b,false);
  CHECKBOX_user_ticks_show_x=glui_labels->add_checkbox_to_panel(PANEL_tick1b,"x",&user_tick_show_x);
  glui_labels->add_column_to_panel(PANEL_tick1b,false);
  CHECKBOX_user_ticks_show_y=glui_labels->add_checkbox_to_panel(PANEL_tick1b,"y",&user_tick_show_y);
  glui_labels->add_column_to_panel(PANEL_tick1b,false);
  CHECKBOX_user_ticks_show_z=glui_labels->add_checkbox_to_panel(PANEL_tick1b,"z",&user_tick_show_z);
  Labels_CB(LABELS_ticks);

  PANEL_tick2 = glui_labels->add_panel_to_panel(ROLLOUT_user_tick,_("Parameters"),true);
  glui_labels->add_statictext_to_panel(PANEL_tick2,"                    x");
  SPINNER_tick_x0=glui_labels->add_spinner_to_panel(PANEL_tick2,_("origin"),GLUI_SPINNER_FLOAT,user_tick_origin);
  SPINNER_tick_xmin=glui_labels->add_spinner_to_panel(PANEL_tick2,_("Min"),GLUI_SPINNER_FLOAT,user_tick_min);
  SPINNER_tick_xmax=glui_labels->add_spinner_to_panel(PANEL_tick2,_("Max"),GLUI_SPINNER_FLOAT,user_tick_max);
  SPINNER_tick_dx0=glui_labels->add_spinner_to_panel(PANEL_tick2,_("Step"),GLUI_SPINNER_FLOAT,user_tick_step);

  glui_labels->add_column_to_panel(PANEL_tick2,false);

  glui_labels->add_statictext_to_panel(PANEL_tick2,"                    y");
  SPINNER_tick_y0=glui_labels->add_spinner_to_panel(PANEL_tick2,"",GLUI_SPINNER_FLOAT,user_tick_origin+1);
  SPINNER_tick_ymin=glui_labels->add_spinner_to_panel(PANEL_tick2,"",GLUI_SPINNER_FLOAT,user_tick_min+1);
  SPINNER_tick_ymax=glui_labels->add_spinner_to_panel(PANEL_tick2,"",GLUI_SPINNER_FLOAT,user_tick_max+1);
  SPINNER_tick_dy0=glui_labels->add_spinner_to_panel(PANEL_tick2,"",GLUI_SPINNER_FLOAT,user_tick_step+1);

  glui_labels->add_column_to_panel(PANEL_tick2,false);

  glui_labels->add_statictext_to_panel(PANEL_tick2,"                    z");
  SPINNER_tick_z0=glui_labels->add_spinner_to_panel(PANEL_tick2,"",GLUI_SPINNER_FLOAT,user_tick_origin+2);
  SPINNER_tick_zmin=glui_labels->add_spinner_to_panel(PANEL_tick2,"",GLUI_SPINNER_FLOAT,user_tick_min+2);
  SPINNER_tick_zmax=glui_labels->add_spinner_to_panel(PANEL_tick2,"",GLUI_SPINNER_FLOAT,user_tick_max+2);
  SPINNER_tick_dz0=glui_labels->add_spinner_to_panel(PANEL_tick2,"",GLUI_SPINNER_FLOAT,user_tick_step+2);

  // -------------- User labels -------------------

  gl=&LABEL_local;
  init_status=LABEL_Init(gl);
  ROLLOUT_user_labels = glui_labels->add_rollout("Labels",false);

  PANEL_LB_panel1 = glui_labels->add_panel_to_panel(ROLLOUT_user_labels,"",GLUI_PANEL_NONE);


  PANEL_LB_panel3 = glui_labels->add_panel_to_panel(ROLLOUT_user_labels,"Labels");
 
  CHECKBOX_LB_visLabels=glui_labels->add_checkbox_to_panel(PANEL_LB_panel3,"Show labels",&visLabels,LB_VISLABELS,Text_Labels_CB);

  PANEL_LB_panel4 = glui_labels->add_panel_to_panel(PANEL_LB_panel3,"",GLUI_PANEL_NONE);
  BUTTON_LB_label_add=glui_labels->add_button_to_panel(PANEL_LB_panel4,"Add",LB_ADD,Text_Labels_CB);
  glui_labels->add_column_to_panel(PANEL_LB_panel4,false);
  BUTTON_LB_label_delete=glui_labels->add_button_to_panel(PANEL_LB_panel4,"Delete",LB_DELETE,Text_Labels_CB);

  LIST_LB_labels=glui_labels->add_listbox_to_panel(PANEL_LB_panel3,"Select",&label_list_index,LB_LIST,Text_Labels_CB);
  {
    labeldata *thislabel;
    int count=0;

    for(thislabel=label_first_ptr->next;thislabel->next!=NULL;thislabel=thislabel->next){
      if(thislabel->labeltype==TYPE_SMV){
        thislabel->glui_id=-1;
        continue;
      }
      thislabel->glui_id=count;
      LIST_LB_labels->add_item(count++,thislabel->name);
    }
  }
  PANEL_LB_panel2 = glui_labels->add_panel_to_panel(PANEL_LB_panel3,"",GLUI_PANEL_NONE);
  EDIT_LB_label_string=glui_labels->add_edittext_to_panel(PANEL_LB_panel2,"Edit:",GLUI_EDITTEXT_TEXT,gl->name,LB_UPDATE,Text_Labels_CB);
  glui_labels->add_column_to_panel(PANEL_LB_panel2,false);
  BUTTON_LB_label_update=glui_labels->add_button_to_panel(PANEL_LB_panel2,"Update",LB_UPDATE,Text_Labels_CB);

  PANEL_LB_panel6 = glui_labels->add_panel_to_panel(PANEL_LB_panel3,"",GLUI_PANEL_NONE);
  BUTTON_LB_label_previous=glui_labels->add_button_to_panel(PANEL_LB_panel6,"Previous",LB_PREVIOUS,Text_Labels_CB);
  glui_labels->add_column_to_panel(PANEL_LB_panel6,false);
  BUTTON_LB_label_next=glui_labels->add_button_to_panel(PANEL_LB_panel6,"Next",LB_NEXT,Text_Labels_CB);

  PANEL_LB_panel5 = glui_labels->add_panel_to_panel(ROLLOUT_user_labels,"",GLUI_PANEL_NONE);
  PANEL_LB_position=glui_labels->add_panel_to_panel(PANEL_LB_panel5,"position");
  SPINNER_LB_x=glui_labels->add_spinner_to_panel(PANEL_LB_position,"x",GLUI_SPINNER_FLOAT,gl->xyz,LB_XYZ,Text_Labels_CB);
  SPINNER_LB_y=glui_labels->add_spinner_to_panel(PANEL_LB_position,"y",GLUI_SPINNER_FLOAT,gl->xyz+1,LB_XYZ,Text_Labels_CB);
  SPINNER_LB_z=glui_labels->add_spinner_to_panel(PANEL_LB_position,"z",GLUI_SPINNER_FLOAT,gl->xyz+2,LB_XYZ,Text_Labels_CB);
  {
    float xmin, ymin, zmin, xmax, ymax, zmax;

    xmin = xbar0ORIG - 0.25*(xbarORIG-xbar0ORIG);
    xmax = xbarORIG + 0.25*(xbarORIG-xbar0ORIG);
    ymin = ybar0ORIG - 0.25*(ybarORIG-ybar0ORIG);
    ymax = ybarORIG + 0.25*(ybarORIG-ybar0ORIG);
    zmin = zbar0ORIG - 0.25*(zbarORIG-zbar0ORIG);
    zmax = zbarORIG + 0.25*(zbarORIG-zbar0ORIG);
    SPINNER_LB_x->set_float_limits(xmin,xmax);
    SPINNER_LB_y->set_float_limits(ymin,ymax);
    SPINNER_LB_z->set_float_limits(zmin,zmax);

  }

  glui_labels->add_column_to_panel(PANEL_LB_panel5,false);
  PANEL_LB_time=glui_labels->add_panel_to_panel(PANEL_LB_panel5,"time");
  SPINNER_LB_time_start=glui_labels->add_spinner_to_panel(PANEL_LB_time,"start",GLUI_SPINNER_FLOAT,gl->tstart_stop,LB_STARTSTOP,Text_Labels_CB);
  SPINNER_LB_time_stop=glui_labels->add_spinner_to_panel(PANEL_LB_time,"stop",GLUI_SPINNER_FLOAT,gl->tstart_stop+1,LB_STARTSTOP,Text_Labels_CB);
  CHECKBOX_LB_label_show_always=glui_labels->add_checkbox_to_panel(PANEL_LB_time,"Show always",&gl->show_always,LB_SHOWALWAYS,Text_Labels_CB);

  PANEL_LB_color=glui_labels->add_panel_to_panel(ROLLOUT_user_labels,"color");
  SPINNER_LB_red=glui_labels->add_spinner_to_panel(PANEL_LB_color,"red",GLUI_SPINNER_INT,gl->rgb,LB_RGB,Text_Labels_CB);
  SPINNER_LB_green=glui_labels->add_spinner_to_panel(PANEL_LB_color,"green",GLUI_SPINNER_INT,gl->rgb+1,LB_RGB,Text_Labels_CB);
  SPINNER_LB_blue=glui_labels->add_spinner_to_panel(PANEL_LB_color,"blue",GLUI_SPINNER_INT,gl->rgb+2,LB_RGB,Text_Labels_CB);
  SPINNER_LB_red->set_int_limits(0,255);
  SPINNER_LB_green->set_int_limits(0,255);
  SPINNER_LB_blue->set_int_limits(0,255);
  CHECKBOX_LB_label_use_foreground=glui_labels->add_checkbox_to_panel(PANEL_LB_color,"Use foreground color",&gl->useforegroundcolor,LB_FOREGROUND,Text_Labels_CB);
  Text_Labels_CB(LB_LIST);

  // -------------- Cube/Tetra intersection test -------------------
#ifdef pp_GEOMTEST
  ROLLOUT_geomtest = glui_labels->add_rollout("Cube/Tetra intersection test",false);
  glui_labels->add_checkbox_to_panel(ROLLOUT_geomtest,"show",&show_geomtest);
  PANEL_geom1=glui_labels->add_panel_to_panel(ROLLOUT_geomtest,"box bounds");

  PANEL_geom1d=glui_labels->add_panel_to_panel(PANEL_geom1,"",GLUI_PANEL_NONE);
  PANEL_geom1a=glui_labels->add_panel_to_panel(PANEL_geom1d,"",GLUI_PANEL_NONE);
  glui_labels->add_column_to_panel(PANEL_geom1d,false);
  PANEL_geom1b=glui_labels->add_panel_to_panel(PANEL_geom1d,"",GLUI_PANEL_NONE);

  PANEL_geom1c=glui_labels->add_panel_to_panel(PANEL_geom1,"",GLUI_PANEL_NONE);

  SPINNER_box_bounds[0]=glui_labels->add_spinner_to_panel(PANEL_geom1a,"xmin",GLUI_SPINNER_FLOAT,box_bounds2,VOL_BOXTRANSLATE,Volume_CB);
  SPINNER_box_bounds[2]=glui_labels->add_spinner_to_panel(PANEL_geom1a,"ymin",GLUI_SPINNER_FLOAT,box_bounds2+2,VOL_BOXTRANSLATE,Volume_CB);
  SPINNER_box_bounds[4]=glui_labels->add_spinner_to_panel(PANEL_geom1a,"zmin",GLUI_SPINNER_FLOAT,box_bounds2+4,VOL_BOXTRANSLATE,Volume_CB);
  SPINNER_box_bounds[1]=glui_labels->add_spinner_to_panel(PANEL_geom1b,"xmax",GLUI_SPINNER_FLOAT,box_bounds2+1,VOL_BOXTRANSLATE,Volume_CB);
  SPINNER_box_bounds[3]=glui_labels->add_spinner_to_panel(PANEL_geom1b,"ymax",GLUI_SPINNER_FLOAT,box_bounds2+3,VOL_BOXTRANSLATE,Volume_CB);
  SPINNER_box_bounds[5]=glui_labels->add_spinner_to_panel(PANEL_geom1b,"zmax",GLUI_SPINNER_FLOAT,box_bounds2+5,VOL_BOXTRANSLATE,Volume_CB);

  SPINNER_box_translate[0]=glui_labels->add_spinner_to_panel(PANEL_geom1c,"translate: x",GLUI_SPINNER_FLOAT,box_translate,VOL_BOXTRANSLATE,Volume_CB);
  glui_labels->add_column_to_panel(PANEL_geom1c,false);
  SPINNER_box_translate[1]=glui_labels->add_spinner_to_panel(PANEL_geom1c,"y",GLUI_SPINNER_FLOAT,box_translate+1,VOL_BOXTRANSLATE,Volume_CB);
  glui_labels->add_column_to_panel(PANEL_geom1c,false);
  SPINNER_box_translate[2]=glui_labels->add_spinner_to_panel(PANEL_geom1c,"z",GLUI_SPINNER_FLOAT,box_translate+2,VOL_BOXTRANSLATE,Volume_CB);
  Volume_CB(VOL_BOXTRANSLATE);

  PANEL_geom2=glui_labels->add_panel_to_panel(ROLLOUT_geomtest,"tetrahedron vertices");
  PANEL_geom2a=glui_labels->add_panel_to_panel(PANEL_geom2,"",GLUI_PANEL_NONE);
  glui_labels->add_column_to_panel(PANEL_geom2,false);
  PANEL_geom2b=glui_labels->add_panel_to_panel(PANEL_geom2,"",GLUI_PANEL_NONE);
  glui_labels->add_column_to_panel(PANEL_geom2,false);
  PANEL_geom2c=glui_labels->add_panel_to_panel(PANEL_geom2,"",GLUI_PANEL_NONE);

  SPINNER_tetra_vertices[0]=glui_labels->add_spinner_to_panel(PANEL_geom2a,"v1 x:",GLUI_SPINNER_FLOAT,tetra_vertices);
  SPINNER_tetra_vertices[3]=glui_labels->add_spinner_to_panel(PANEL_geom2a,"v2 x:",GLUI_SPINNER_FLOAT,tetra_vertices+3);
  SPINNER_tetra_vertices[6]=glui_labels->add_spinner_to_panel(PANEL_geom2a,"v3 x:",GLUI_SPINNER_FLOAT,tetra_vertices+6);
  SPINNER_tetra_vertices[9]=glui_labels->add_spinner_to_panel(PANEL_geom2a,"v4 x:",GLUI_SPINNER_FLOAT,tetra_vertices+9);

  SPINNER_tetra_vertices[1]=glui_labels->add_spinner_to_panel(PANEL_geom2b,"y:",GLUI_SPINNER_FLOAT,tetra_vertices+1);
  SPINNER_tetra_vertices[4]=glui_labels->add_spinner_to_panel(PANEL_geom2b,"y:",GLUI_SPINNER_FLOAT,tetra_vertices+4);
  SPINNER_tetra_vertices[7]=glui_labels->add_spinner_to_panel(PANEL_geom2b,"y:",GLUI_SPINNER_FLOAT,tetra_vertices+7);
  SPINNER_tetra_vertices[10]=glui_labels->add_spinner_to_panel(PANEL_geom2b,"y:",GLUI_SPINNER_FLOAT,tetra_vertices+10);

  SPINNER_tetra_vertices[2]=glui_labels->add_spinner_to_panel(PANEL_geom2c,"z:",GLUI_SPINNER_FLOAT,tetra_vertices+2);
  SPINNER_tetra_vertices[5]=glui_labels->add_spinner_to_panel(PANEL_geom2c,"z:",GLUI_SPINNER_FLOAT,tetra_vertices+5);
  SPINNER_tetra_vertices[8]=glui_labels->add_spinner_to_panel(PANEL_geom2c,"z:",GLUI_SPINNER_FLOAT,tetra_vertices+8);
  SPINNER_tetra_vertices[11]=glui_labels->add_spinner_to_panel(PANEL_geom2c,"z:",GLUI_SPINNER_FLOAT,tetra_vertices+11);

  glui_labels->add_checkbox_to_panel(ROLLOUT_geomtest,"show intersection",&show_intersection);
  {
    int i;
    
    for(i=0;i<10;i++){
      char label[100];

      sprintf(label,"face %i",i);
      glui_labels->add_checkbox_to_panel(ROLLOUT_geomtest,label,tetrabox_vis+i);
    }
  }

#endif

  // -------------- 

  PANEL_label2 = glui_labels->add_panel("",false);
  glui_labels->add_column_to_panel(PANEL_label2,false);

  BUTTON_label_3=glui_labels->add_button_to_panel(PANEL_label2,_("Save settings"),SAVE_SETTINGS,Labels_CB);
  glui_labels->add_column_to_panel(PANEL_label2,false);

  BUTTON_label_4=glui_labels->add_button_to_panel(PANEL_label2,_("Close"),LABELS_close,Labels_CB);

  glui_labels->set_main_gfx_window( main_window );
}

/* ------------------ update_fileload  ------------------------ */

extern "C" void update_fileload(void){
  int i;
  partdata *parti;
  slicedata *slicei;
  isodata *isoi;
  patchdata *patchi;
  smoke3ddata *smoke3di;
  plot3ddata *plot3di;
  vslicedata *vslicei;

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
  for(i=0;i<nsliceinfo;i++){
    slicei = sliceinfo+i;
    if(slicei->loaded==1){
      nsliceloaded++;
    }
  }

  nvsliceloaded=0;
  for(i=0;i<nvsliceinfo;i++){
    vslicei = vsliceinfo+i;
    if(vslicei->loaded==1){
      nvsliceloaded++;
    }
  }

  nisoloaded=0;
  for(i=0;i<nisoinfo;i++){
    isoi = isoinfo+i;
    if(isoi->loaded==1){
      nisoloaded++;
    }
  }

  npatchloaded=0;
  for(i=0;i<npatchinfo;i++){
    patchi = patchinfo+i;
    if(patchi->loaded==1){
      npatchloaded++;
    }
  }

  nsmoke3dloaded=0;
  for(i=0;i<nsmoke3dinfo;i++){
    smoke3di = smoke3dinfo+i;
    if(smoke3di->loaded==1){
      nsmoke3dloaded++;
    }
  }

  nplot3dloaded=0;
  for(i=0;i<nplot3dinfo;i++){
    plot3di = plot3dinfo+i;
    if(plot3di->loaded==1){
      nplot3dloaded++;
    }
  }
}

/* ------------------ update_showhidebuttons ------------------------ */

extern "C" void update_showhidebuttons(void){

  update_fileload();
  if(CHECKBOX_label_3!=NULL){
    CHECKBOX_label_3->set_int_val(hide_overlaps);
  }
  if(BUTTON_PART!=NULL){
    if(npartloaded==0){
      BUTTON_PART->disable();
    }
    else{
      BUTTON_PART->enable();
    }
  }
  
  if(BUTTON_SLICE!=NULL){
    if(nsliceloaded==0){
      BUTTON_SLICE->disable();
    }
    else{
      BUTTON_SLICE->enable();
    }
  }

  if(BUTTON_VSLICE!=NULL){
    if(nvsliceloaded==0){
      BUTTON_VSLICE->disable();
    }
    else{
      BUTTON_VSLICE->enable();
    }
  }

  if(BUTTON_ISO!=NULL){
    if(nisoloaded==0){
      BUTTON_ISO->disable();
    }
    else{
      BUTTON_ISO->enable();
    }
  }

  if(BUTTON_BOUNDARY!=NULL){
    if(npatchloaded==0){
      BUTTON_BOUNDARY->disable();
    }
    else{
      BUTTON_BOUNDARY->enable();
    }
  }

  if(BUTTON_3DSMOKE!=NULL){
    if(nsmoke3dloaded==0){
      BUTTON_3DSMOKE->disable();
    }
    else{
      BUTTON_3DSMOKE->enable();
    }
  }

  if(BUTTON_PLOT3D!=NULL){
    if(nplot3dloaded==0){
      BUTTON_PLOT3D->disable();
    }
    else{
      BUTTON_PLOT3D->enable();
    }
  }

  if(nplot3dloaded==0&&nsmoke3dloaded==0&&nisoloaded==0&&nsliceloaded==0&&npartloaded==0&&npatchloaded==0){
    if(RADIO_showhide!=NULL)RADIO_showhide->disable();
  }
  else{
    if(RADIO_showhide!=NULL)RADIO_showhide->enable();
  }
}

/* ------------------ hide_glui_display ------------------------ */

extern "C" void hide_glui_display(void){
  if(glui_labels!=NULL)glui_labels->hide();
  showdisplay_dialog_save=showdisplay_dialog;
  showdisplay_dialog=0;
  updatemenu=1;
}

/* ------------------ show_glui_display ------------------------ */

extern "C" void show_glui_display(void){
  showdisplay_dialog=1;
  if(glui_labels!=NULL)glui_labels->show();
}

/* ------------------ Text_labels_CB ------------------------ */

void Text_Labels_CB(int var){
  labeldata *thislabel,*gl,*new_label;
  int count;
  char name[300];
  int len;

  len=sizeof(GLUI_String);

  gl=&LABEL_local;
  switch (var){
    case LB_VISLABELS:
      updatemenu=1;
      break;
    case LB_UPDATE:
      for(thislabel=label_first_ptr->next;thislabel->next!=NULL;thislabel=thislabel->next){
        if(thislabel->glui_id<0)continue;
        LIST_LB_labels->delete_item(thislabel->glui_id);
      }
      strcpy(LABEL_global_ptr->name,gl->name);
      //LABEL_resort(LABEL_global_ptr);

      count=0;
      for(thislabel=label_first_ptr->next;thislabel->next!=NULL;thislabel=thislabel->next){
        if(thislabel->labeltype==TYPE_SMV)continue;
        thislabel->glui_id=count;
        LIST_LB_labels->add_item(count++,thislabel->name);
      }
      break;
    case LB_STARTSTOP:
      memcpy(LABEL_global_ptr->tstart_stop,gl->tstart_stop,2*sizeof(float));
      break;
    case LB_SHOWALWAYS:
      memcpy(&LABEL_global_ptr->show_always,&gl->show_always,sizeof(int));
      break;
    case LB_FOREGROUND:
      memcpy(&LABEL_global_ptr->useforegroundcolor,&gl->useforegroundcolor,sizeof(int));
      break;
    case LB_PREVIOUS:
      new_label=LABEL_get(LIST_LB_labels->curr_text);
      new_label=LABEL_Previous(new_label);
      if(new_label==NULL)break;
      LABEL_global_ptr=new_label;
      if(new_label!=NULL){
        LABEL_copy(gl,new_label);
        update_glui_label_text();
      }
      break;
    case LB_NEXT:
      new_label=LABEL_get(LIST_LB_labels->curr_text);
      new_label=LABEL_Next(new_label);
      if(new_label==NULL)break;
      LABEL_global_ptr=new_label;
      if(new_label!=NULL){
        LABEL_copy(gl,new_label);
        update_glui_label_text();
      }
      break;
    case LB_LIST:
      new_label=LABEL_get(LIST_LB_labels->curr_text);
      LABEL_global_ptr=new_label;
      if(new_label!=NULL){
        LABEL_copy(gl,new_label);
      }
      update_glui_label_text();
      break;
    case LB_ADD:
      updatemenu=1;
      if(LABEL_Get_Nuserlabels()>0){
        strcpy(name,"copy of ");
        strcat(name,gl->name);
        strcpy(gl->name,name);
      }
      else{
        gl=&LABEL_default;
      }
      gl->labeltype=TYPE_INI;
      for(thislabel=label_first_ptr->next;thislabel->next!=NULL;thislabel=thislabel->next){
        if(thislabel->glui_id<0)continue;
        LIST_LB_labels->delete_item(thislabel->glui_id);
      }
      LABEL_insert(gl);
      count=0;
      for(thislabel=label_first_ptr->next;thislabel->next!=NULL;thislabel=thislabel->next){
        if(thislabel->labeltype==TYPE_SMV)continue;
        thislabel->glui_id=count;
        LIST_LB_labels->add_item(count++,thislabel->name);
      }
      Text_Labels_CB(LB_LIST);
      break;
    case LB_DELETE:
      strcpy(name,LIST_LB_labels->curr_text);
      for(thislabel=label_first_ptr->next;thislabel->next!=NULL;thislabel=thislabel->next){
        if(thislabel->glui_id<0)continue;
        LIST_LB_labels->delete_item(thislabel->glui_id);
      }
      thislabel=LABEL_get(name);
      if(thislabel!=NULL){
        LABEL_delete(thislabel);
      }
      count=0;
      for(thislabel=label_first_ptr->next;thislabel->next!=NULL;thislabel=thislabel->next){
        if(thislabel->labeltype==TYPE_SMV)continue;
        thislabel->glui_id=count;
        LIST_LB_labels->add_item(count++,thislabel->name);
      }
      Text_Labels_CB(LB_LIST);
      break;
    case LB_RGB:
      gl->frgb[0]=gl->rgb[0]/255.0;
      gl->frgb[1]=gl->rgb[1]/255.0;
      gl->frgb[2]=gl->rgb[2]/255.0;
      memcpy(LABEL_global_ptr->frgb,gl->frgb,3*sizeof(float));
      memcpy(LABEL_global_ptr->rgb,gl->rgb,3*sizeof(int));
      break;
    case LB_XYZ:
      memcpy(LABEL_global_ptr->xyz,gl->xyz,3*sizeof(float));
      break;
    default:
      ASSERT(0);
      break;
  }
}

/* ------------------ Volume_CB ------------------------ */
#ifdef pp_GEOMTEST
extern "C" void Volume_CB(int var){
  switch (var){
    case VOL_BOXTRANSLATE:
      box_bounds[0]=box_bounds2[0]+box_translate[0];
      box_bounds[1]=box_bounds2[1]+box_translate[0];
      box_bounds[2]=box_bounds2[2]+box_translate[1];
      box_bounds[3]=box_bounds2[3]+box_translate[1];
      box_bounds[4]=box_bounds2[4]+box_translate[2];
      box_bounds[5]=box_bounds2[5]+box_translate[2];
      break;
    default:
      ASSERT(0);
      break;
  }
}
#endif

/* ------------------ Labels_CB ------------------------ */

extern "C" void Labels_CB(int var){
  updatemenu=1;
  switch (var){
    case FLIP:
      colorbarflip = 1 - colorbarflip;
      ColorBarMenu(COLORBAR_FLIP);
      break;
  case LABELS_hide_overlaps:
    updatefacelists=1;
    updatehiddenfaces=1;
    UpdateHiddenFaces();
    glutPostRedisplay();
    break;
#ifdef pp_BETA
  case LABELS_drawface:
    /*
    if(use_new_drawface==1){
      CHECKBOX_cullgeom->enable();
      if(cullgeom==1){
        SPINNER_cullgeom_portsize->enable();
      }
      else{
        SPINNER_cullgeom_portsize->disable();
      }
    }
    else{
      CHECKBOX_cullgeom->disable();
      SPINNER_cullgeom_portsize->disable();
    }
    update_initcullgeom=1;
    set_cull_vis();
    */
    updatefacelists=1;
    break;
#endif

  case LABELS_usertick:
    CHECKBOX_vis_user_ticks2->set_int_val(vis_user_ticks);
    break;
  case LABELS_usertick2:
    CHECKBOX_vis_user_ticks->set_int_val(vis_user_ticks);
    break;
  case LABELS_version:
    break;
  case LABELS_meshlabel:
    break;
  case SAVE_SETTINGS:
    writeini(LOCAL_INI,NULL);
    break;
  case LABELS_HMS:
    break;
  case LABELS_showall:
    LabelMenu(4);
    break;
  case LABELS_hideall:
    LabelMenu(5);
    break;
  case LABELS_flip:
    background_flip = 1 - background_flip;
    ShowHideMenu(15);
    break;
  case LABELS_shade:
    setbw = 1 - setbw;
    ColorBarMenu(COLORBAR_TOGGLE_BW);
    break;
  case LABELS_transparent:
    break;
  case LABELS_close:
    hide_glui_display();
    break;
  case LABELS_fontsize:
    FontMenu(fontindex);
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

  if(CHECKBOX_LB_visLabels!=NULL)CHECKBOX_LB_visLabels->set_int_val(visLabels);
  if(CHECKBOX_vis_user_ticks!=NULL)CHECKBOX_vis_user_ticks->set_int_val(vis_user_ticks);
  if(CHECKBOX_labels_hrrlabel!=NULL)CHECKBOX_labels_hrrlabel->set_int_val(visHRRlabel);
  if(CHECKBOX_labels_hrrcutoff!=NULL)CHECKBOX_labels_hrrcutoff->set_int_val(show_hrrcutoff);
  if(CHECKBOX_labels_title!=NULL)CHECKBOX_labels_title->set_int_val(visTitle);
  if(CHECKBOX_labels_colorbar!=NULL)CHECKBOX_labels_colorbar->set_int_val(visColorbarLabels);
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
  if(CHECKBOX_labels_version!=NULL)CHECKBOX_labels_version->set_int_val(gversion);
  if(CHECKBOX_labels_meshlabel!=NULL)CHECKBOX_labels_meshlabel->set_int_val(visBlocklabel);
  if(CHECKBOX_vis_user_ticks2!=NULL)CHECKBOX_vis_user_ticks2->set_int_val(vis_user_ticks);
}


/* ------------------ update_colorbarflip ------------------------ */

extern "C" void update_colorbarflip(void){
  CHECKBOX_colorbarflip->set_int_val(colorbarflip);
}

/* ------------------ update_colorbar_list2 ------------------------ */

extern "C" void update_colorbar_list2(void){
  LIST_colorbar2->set_int_val(selectedcolorbar_index2);
}

/* ------------------ add_colorbar_list2 ------------------------ */

extern "C" void add_colorbar_list2(int index, char *label){
  LIST_colorbar2->add_item(index,label);
}

/* ------------------ set_colorbar_list_index ------------------------ */

extern "C" void set_colorbar_list_index(int val){
  if(LIST_colorbar2!=NULL)LIST_colorbar2->set_int_val(val);
}

/* ------------------ get_colorbar_list_index ------------------------ */

extern "C" int get_colorbar_list_index(void){
  return LIST_colorbar2->get_int_val();
}

/* ------------------ update_axislabels_smooth ------------------------ */

extern "C" void update_axislabels_smooth(void){
  CHECKBOX_axislabels_smooth->set_int_val(axislabels_smooth);
}

/* ------------------ transparency ------------------------ */

extern "C" void update_transparency(void){
  CHECKBOX_transparentflag->set_int_val(use_transparency_data);
}

/* ------------------ FileShow_CB ------------------------ */

extern "C" void FileShow_CB(int var){
  updatemenu=1;
  switch (var){
    case  FILESHOW_plot3d:
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
    case FILESHOW_evac:
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
    case  FILESHOW_particle:
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
    case  FILESHOW_slice:
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
    case  FILESHOW_vslice:
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
    case  FILESHOW_boundary:
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
    case  FILESHOW_3dsmoke:
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
    case  FILESHOW_isosurface:
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
    default:
      break;
  }
}

/* ------------------ Extreme_CB ------------------------ */

extern "C" void Extreme_CB(int var){
  colorbardata *cbi;
  unsigned char *rgb_nodes;
  int i;

  switch (var){
    case COLORBAR_EXTREME:
      if(show_extreme_mindata==1){
        if(SPINNER_down_red!=NULL)SPINNER_down_red->enable();
        if(SPINNER_down_green!=NULL)SPINNER_down_green->enable();
        if(SPINNER_down_blue!=NULL)SPINNER_down_blue->enable();
      }
      else{
        if(SPINNER_down_red!=NULL)SPINNER_down_red->disable();
        if(SPINNER_down_green!=NULL)SPINNER_down_green->disable();
        if(SPINNER_down_blue!=NULL)SPINNER_down_blue->disable();
      }
      if(show_extreme_maxdata==1){
        if(SPINNER_up_red!=NULL)SPINNER_up_red->enable();
        if(SPINNER_up_green!=NULL)SPINNER_up_green->enable();
        if(SPINNER_up_blue!=NULL)SPINNER_up_blue->enable();
      }
      else{
        if(SPINNER_up_red!=NULL)SPINNER_up_red->disable();
        if(SPINNER_up_green!=NULL)SPINNER_up_green->disable();
        if(SPINNER_up_blue!=NULL)SPINNER_up_blue->disable();
      }
      if(colorbartype<0||colorbartype>=ncolorbars)return;
      cbi = colorbarinfo + colorbartype;
      remapcolorbar(cbi);
      UpdateRGBColors(COLORBAR_INDEX_NONE);
      updatemenu=1;
      break;
    case COLORBAR_EXTREME_RGB:
      if(colorbartype<0||colorbartype>=ncolorbars)return;
      cbi = colorbarinfo + colorbartype;

      rgb_nodes=rgb_above_max;
      for(i=0;i<3;i++){
        rgb_nodes[i]=cb_up_rgb[i];
      }
      rgb_nodes=rgb_below_min;
      for(i=0;i<3;i++){
        rgb_nodes[i]=cb_down_rgb[i];
      }
      remapcolorbar(cbi);
      UpdateRGBColors(COLORBAR_INDEX_NONE);
      break;
    default:
      break;
  }
}

/* ------------------ update_extreme_vals ------------------------ */

extern "C" void update_extreme_vals(void){
  unsigned char *rgb_local;

  rgb_local = rgb_below_min;
  if(SPINNER_down_red!=NULL)SPINNER_down_red->set_int_val(  (int)(rgb_local[0]));
  if(SPINNER_down_green!=NULL)SPINNER_down_green->set_int_val(  (int)(rgb_local[1]));
  if(SPINNER_down_blue!=NULL)SPINNER_down_blue->set_int_val(  (int)(rgb_local[2]));

  rgb_local = rgb_above_max;
  if(SPINNER_up_red!=NULL)SPINNER_up_red->set_int_val(  (int)(rgb_local[0]));
  if(SPINNER_up_green!=NULL)SPINNER_up_green->set_int_val(  (int)(rgb_local[1]));
  if(SPINNER_up_blue!=NULL)SPINNER_up_blue->set_int_val(  (int)(rgb_local[2]));
}

/* ------------------ update_camera_label ------------------------ */

extern "C" void update_extreme(void){
  if(CHECKBOX_show_extreme_mindata!=NULL){
    CHECKBOX_show_extreme_mindata->set_int_val(show_extreme_mindata);
  }
  if(CHECKBOX_show_extreme_maxdata!=NULL){
    CHECKBOX_show_extreme_maxdata->set_int_val(show_extreme_maxdata);
  }
  Extreme_CB(COLORBAR_EXTREME);
}
