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

// svn revision character string
extern "C" char glui_bounds_revision[]="$Revision$";
// $Date$ $Author$

extern "C" void update_slice_contours(int slice_type_index, float line_min, float line_max, int nline_values);
extern "C" int file_exist(char *file);
extern "C" void keyboard(unsigned char key, int x, int y);
extern "C" void glui_script_disable(void);
extern "C" void glui_script_enable();
extern "C" void writeini(int var);
extern "C" void ScriptMenu(int var);
extern "C" void update_glui_smoke3dframestep(void);
extern "C" void ParticlePropShowMenu(int value);
extern "C" void colorbar_global2local(void);
extern "C" void ReloadMenu(int var);

void COLORBAR_CB(int var);
void SETslicemax(int setslicemax, float slicemax);
void SETslicemin(int setslicemin, float slicemin);
void BUTTON_hide_CB(int var);
void PART_CB(int var), Bound_CB(int var), Slice_CB(int var), PLOT3D_CB(int var), Iso_CB(int var), Smoke3D_CB(int var);
void Time_CB(int var);
void Script_CB(int var);
void boundmenu(GLUI_Rollout **rollout_bound, GLUI_Rollout **rollout_chop, GLUI_Panel *panel, char *button_title,
          GLUI_EditText **con_min,GLUI_EditText **con_max,
          GLUI_RadioGroup **con_setmin,GLUI_RadioGroup **con_setmax,
          GLUI_Checkbox **con_setchopmin, GLUI_Checkbox **con_setchopmax,
          GLUI_EditText **con_chopmin, GLUI_EditText **con_chopmax,
          int *setminval, int *setmaxval,
          float *minval, float *maxval,
          int *setchopminval, int *setchopmaxval,
          float *chopminval, float *chopmaxval,
          int updatebounds,
          int truncatebounds,
          GLUI_Update_CB FILE_CB);

GLUI_Rollout *rollout_slice_bound=NULL;
GLUI_Rollout *rollout_slice_chop=NULL;

#define TRUNCATEMIN 8
#define TRUNCATEMAX 9
#define SETVALMIN 1
#define SETVALMAX 2
#define VALMIN 3
#define VALMAX 4
#define FILETYPEINDEX 5
#define FILEUPDATE 6
#define FILERELOAD 7
#define UPDATEPLOT 10
#define PLOTISO 11
#define SHOWCHAR 12
#define CHOPVALMIN 13
#define CHOPVALMAX 14
#define SETCHOPMINVAL 15
#define SETCHOPMAXVAL 16
#define CHOPUPDATE 17
#define FRAMELOADING 18
#define STREAKLENGTH 19
#define TRACERS 21
#define PLOTISOTYPE 22
#ifdef pp_TRANSFORM
#define TRANSFORM_SLICE 23
#define RESET_SLICE 24
#endif
#define AVERAGE_DATA 201
#define TURB_DATA 202
#define UNLOAD_QDATA 203
#define SET_TIME 204
#define TBOUNDS 205
#define TBOUNDS_USE 206
#define RELOAD_DATA 207

#define SCRIPT_START 31
#define SCRIPT_STOP 32
#define SCRIPT_LIST 33
#define SCRIPT_SAVEINI 34
#define SCRIPT_EDIT_INI 35
#define SCRIPT_SETSUFFIX 36
#define SCRIPT_RUNSCRIPT 37
#define SCRIPT_LOADINI 38
#define SCRIPT_LISTINI 39
#define SCRIPT_EDIT_RENDERDIR 40
#define SCRIPT_RENDER 41
#define SCRIPT_RENDER_SUFFIX 42
#define SCRIPT_RENDER_DIR 43

#define SAVE_SETTINGS 99
#define CLOSE 98
#define COMPRESS_FILES 97
#define OVERWRITE 96
#define COMPRESS_AUTOLOADED 91
#define ERASE 95
#define STARTUP 94
#define SAVE_FILE_LIST 93
#define LOAD_FILES 92
#define COLORBAR_EXTREME2 109
#define TRANSPARENTLEVEL 110
#define COLORBAR_LINES 111
#define COLORBAR_LIST2 112
#define COLORBAR_SMOOTH 113

#define UPDATE_VECTOR 101

#define TRUNCATEBOUNDS 1
#define DONT_TRUNCATEBOUNDS 0
#define UPDATEBOUNDS 1
#define DONT_UPDATEBOUNDS 0

#ifdef pp_SLICECONTOURS
#define LINE_CONTOUR_VALUE 301
#define UPDATE_LINE_CONTOUR_VALUE 302
#endif

GLUI_Button *BUTTON_compress=NULL;
#ifdef pp_TRANSFORM
GLUI_Panel *panel_slice_transform=NULL;
#endif
GLUI_Panel *panel_pan1=NULL;
GLUI_Panel *panel_pan2=NULL;
GLUI_Panel *panel_pan3=NULL;
GLUI_Panel *panel_contours=NULL;
GLUI_Panel *panel_isosurface=NULL;
GLUI_Panel *panel_vector=NULL;
GLUI_Panel *panel_slice_vector=NULL;
GLUI_Spinner *SPINNER_transparentlevel=NULL;
#ifdef pp_SLICECONTOURS
GLUI_Panel *panel_line_contour=NULL;
GLUI_Spinner *SPINNER_line_contour_num=NULL;
GLUI_Spinner *SPINNER_line_contour_min=NULL;
GLUI_Spinner *SPINNER_line_contour_max=NULL;
GLUI_Button *BUTTON_update_line_contour=NULL;
#endif
GLUI_Spinner *SPINNER_timebounds=NULL;
GLUI_Panel *panel_script1=NULL;
GLUI_Panel *panel_script1a=NULL;
GLUI_Panel *panel_script1b=NULL;
GLUI_Panel *panel_script1c=NULL;
GLUI_Panel *panel_script2=NULL;
GLUI_Panel *panel_script2a=NULL;
GLUI_Panel *panel_script2b=NULL;
GLUI_Panel *panel_script3=NULL;
GLUI_Listbox *LIST_scriptlist=NULL;
GLUI_Listbox *LISTBOX_colorbar2=NULL;
GLUI_Listbox *LIST_ini_list=NULL;
GLUI_Button *BUTTON_ini_load=NULL;
GLUI_Rollout *rollout_SCRIPT=NULL;
GLUI_Button *BUTTON_script_stop=NULL;
GLUI_Button *BUTTON_script_start=NULL;
GLUI_Button *BUTTON_script_saveini=NULL;
GLUI_Button *BUTTON_script_render=NULL;
GLUI_EditText *EDIT_ini=NULL;
GLUI_EditText *EDIT_renderdir=NULL;
GLUI_EditText *EDIT_rendersuffix=NULL;
GLUI_Button *BUTTON_script_setsuffix=NULL;
GLUI_Button *BUTTON_script_runscript=NULL;
GLUI_Rollout *rollout_BOUNDARY=NULL;
GLUI_Rollout *rollout_AUTOLOAD=NULL;
GLUI_Rollout *rollout_compress=NULL;
GLUI *glui_bounds=NULL;
GLUI_Rollout *panel_plot3d=NULL,*panel_evac=NULL,*panel_part=NULL,*panel_slice=NULL,*panel_bound=NULL,*panel_iso=NULL,*panel_smoke3d=NULL;
GLUI_Rollout *panel_time=NULL,*panel_colorbar=NULL;

GLUI_RadioGroup *bf_rlist=NULL, *p3rlist=NULL,*slice_rlist=NULL;
GLUI_RadioGroup *part5_rlist=NULL;
GLUI_RadioGroup *iso_isotype=NULL;
GLUI_RadioGroup *plot3d_isotype=NULL;
GLUI_RadioButton *plot3d_iso_hidden=NULL;
GLUI_RadioGroup *plot3d_display=NULL;
GLUI_EditText *con_slice_min=NULL, *con_slice_max=NULL;
GLUI_EditText *con_slice_chopmin=NULL, *con_slice_chopmax=NULL;
GLUI_EditText *con_part_chopmin=NULL, *con_part_chopmax=NULL;
GLUI_RadioGroup *con_slice_setmin=NULL, *con_slice_setmax=NULL;
GLUI_Checkbox *showchar_checkbox=NULL, *showonlychar_checkbox;
GLUI_Checkbox *CHECKBOX_showcolorbarlines=NULL;
GLUI_Checkbox *CHECKBOX_axissmooth=NULL;
GLUI_Checkbox *CHECKBOX_extreme2=NULL;
GLUI_Checkbox *startup_checkbox=NULL;
GLUI_Checkbox *check_overwrite_all=NULL;
GLUI_Checkbox *check_compress_autoloaded=NULL;
GLUI_Checkbox *check_erase_all=NULL;
GLUI_Checkbox *check_multi_task=NULL;
GLUI_Checkbox *con_slice_setchopmin=NULL;
GLUI_Checkbox *con_slice_setchopmax=NULL;
GLUI_Checkbox *con_p3_setchopmin=NULL, *con_p3_setchopmax=NULL;
GLUI_Checkbox *con_bf_setchopmin=NULL, *con_bf_setchopmax=NULL;
GLUI_Checkbox *con_part_setchopmin=NULL, *con_part_setchopmax=NULL;
GLUI_Checkbox *showtracer_checkbox=NULL;
GLUI_Checkbox *CHECKBOX_cellcenter_slice_interp=NULL;
GLUI_Checkbox *CHECKBOX_skip_subslice=NULL;
GLUI_Checkbox *CHECKBOX_turb_slice=NULL;
GLUI_Checkbox *CHECKBOX_average_slice=NULL;
GLUI_Checkbox *CHECKBOX_unload_qdata=NULL;
GLUI_Checkbox *CHECKBOX_use_tload_begin=NULL;
GLUI_Checkbox *CHECKBOX_use_tload_end=NULL;
GLUI_Checkbox *CHECKBOX_use_tload_skip=NULL;
GLUI_Spinner *SPINNER_tload_begin=NULL;
GLUI_Spinner *SPINNER_tload_end=NULL;
GLUI_Spinner *SPINNER_tload_skip=NULL;
GLUI_Panel *panel_time2=NULL;
GLUI_Panel *panel_time1a=NULL;
GLUI_Panel *panel_time2a=NULL;
GLUI_Panel *panel_time2b=NULL;
GLUI_Panel *panel_time2c=NULL;
GLUI_Button *BUTTON_RELOAD=NULL;
GLUI_Button *BUTTON_SETTIME=NULL;

GLUI_Spinner *SPINNER_plot3d_vectorpointsize=NULL,*SPINNER_plot3d_vectorlinewidth=NULL;
GLUI_Spinner *SPINNER_sliceaverage=NULL;
GLUI_Spinner *SPINNER_partpointstep=NULL;
GLUI_Spinner *SPINNER_smoke3dzipstep=NULL;
GLUI_Spinner *SPINNER_slicezipstep=NULL;
GLUI_Spinner *SPINNER_isozipstep=NULL;
GLUI_Spinner *SPINNER_boundzipstep=NULL;
GLUI_Spinner *SPINNER_partstreaklength=NULL;
GLUI_Spinner *SPINNER_partpointsize=NULL;
GLUI_Spinner *SPINNER_isopointsize=NULL;
GLUI_Spinner *SPINNER_isolinewidth=NULL;
GLUI_Spinner *SPINNER_plot3dpointsize=NULL;
GLUI_Spinner *SPINNER_plot3dlinewidth=NULL;
GLUI_Spinner *SPINNER_streaklinewidth=NULL;
GLUI_Spinner *SPINNER_vectorpointsize=NULL;
GLUI_Spinner *SPINNER_vectorlinewidth=NULL;



GLUI_EditText *con_bf_min=NULL, *con_bf_max=NULL;
GLUI_EditText *con_bf_chopmin=NULL, *con_bf_chopmax=NULL;
GLUI_RadioGroup *con_bf_setmin=NULL, *con_bf_setmax=NULL;


GLUI_EditText *con_part_min=NULL, *con_part_max=NULL;
GLUI_RadioGroup *con_part_setmin=NULL, *con_part_setmax=NULL;

GLUI_EditText *con_p3_min=NULL, *con_p3_max=NULL;
GLUI_EditText *con_p3_chopmin=NULL, *con_p3_chopmax=NULL;
GLUI_RadioGroup *con_p3_setmin=NULL, *con_p3_setmax=NULL;

/* ------------------ update_colorbar_smooth ------------------------ */

extern "C" void update_colorbar_smooth(void){
  CHECKBOX_axissmooth->set_int_val(axissmooth);
}

/* ------------------ update_colorbar_list2 ------------------------ */

extern "C" void update_colorbar_list2(void){
  LISTBOX_colorbar2->set_int_val(selectedcolorbar_index2);
}

/* ------------------ update_update_extreme2 ------------------------ */

extern "C" void update_extreme2(void){
  if(CHECKBOX_extreme2!=NULL)CHECKBOX_extreme2->set_int_val(show_extremedata);
}

/* ------------------ update_glui_plot3d ------------------------ */

extern "C" void update_glui_plot3d(void){
  PLOT3D_CB(UNLOAD_QDATA);
}

/* ------------------ PART_CB_INIT ------------------------ */

extern "C" void PART_CB_INIT(void){
  PART_CB(FILETYPEINDEX);
}

/* ------------------ glui_bounds_setup ------------------------ */

extern "C" void glui_bounds_setup(int main_window){  
  int i;
  int firstcount=0;
  int nradio;

  overwrite_all=0;
  glui_bounds = GLUI_Master.create_glui( "File/Bounds Settings",0,0,0 );
  if(showbounds==0)glui_bounds->hide();

  /*  3d smoke   */

  if(nsmoke3d_files>0){
    panel_smoke3d = glui_bounds->add_rollout("3D Smoke",false);
  }

  /*  Boundary File Bounds   */

  if(npatch_files>0){
    glui_active=1;
    panel_bound = glui_bounds->add_rollout("Boundary",false);

    bf_rlist = glui_bounds->add_radiogroup_to_panel(panel_bound,&list_patch_index,FILETYPEINDEX,Bound_CB);
    nradio=0;
    for(i=0;i<npatch_files;i++){
      if(patchinfo[i].firstshort==1)nradio++;
    }
    if(nradio>1){
      for(i=0;i<npatch_files;i++){
        if(patchinfo[i].firstshort==1)glui_bounds->add_radiobutton_to_group(bf_rlist,patchinfo[i].label.shortlabel);
      }
      if(activate_threshold==1){
        glui_bounds->add_separator_to_panel(panel_bound);
        showchar_checkbox=glui_bounds->add_checkbox_to_panel(panel_bound,"Show temp threshold",
          &vis_threshold,SHOWCHAR,Bound_CB);
        showonlychar_checkbox=glui_bounds->add_checkbox_to_panel(panel_bound,"Show only temp threshold",
          &vis_onlythreshold,SHOWCHAR,Bound_CB);
        glui_bounds->add_spinner_to_panel(panel_bound,"temp threshold (C)",GLUI_SPINNER_FLOAT,
          &temp_threshold);
        Bound_CB(SHOWCHAR);
      }
      glui_bounds->add_column_to_panel(panel_bound,false);
    }
    else{
      if(activate_threshold==1){
        showchar_checkbox=glui_bounds->add_checkbox_to_panel(panel_bound,"Show temp threshold",&vis_threshold,SHOWCHAR,Bound_CB);
        showonlychar_checkbox=glui_bounds->add_checkbox_to_panel(panel_bound,"Show only threshold",&vis_onlythreshold,SHOWCHAR,Bound_CB);
        glui_bounds->add_spinner_to_panel(panel_bound,"temp threshold (C)",GLUI_SPINNER_FLOAT,&temp_threshold);
        Bound_CB(SHOWCHAR);
      }
    }

    boundmenu(&rollout_BOUNDARY,NULL,panel_bound,"Reload Boundary File(s)",
      &con_bf_min,&con_bf_max,&con_bf_setmin,&con_bf_setmax,
      NULL,NULL,NULL,NULL,
      &setpatchmin,&setpatchmax,&patchmin,&patchmax,
      NULL,NULL,NULL,NULL,
      DONT_UPDATEBOUNDS,DONT_TRUNCATEBOUNDS,
      Bound_CB);

  }

  /*  Iso File Load Bounds   */

  if(niso_files>0){
    panel_iso = glui_bounds->add_rollout("Isosurface",false);
    
    SPINNER_isopointsize=glui_bounds->add_spinner_to_panel(panel_iso,"Point size",GLUI_SPINNER_FLOAT,
      &isopointsize);
    SPINNER_isopointsize->set_float_limits(1.0,10.0);

    SPINNER_isolinewidth=glui_bounds->add_spinner_to_panel(panel_iso,"Line width",GLUI_SPINNER_FLOAT,
      &isolinewidth);
    SPINNER_isolinewidth->set_float_limits(1.0,10.0);

    iso_isotype = glui_bounds->add_radiogroup_to_panel(panel_iso,&visAIso);

    glui_bounds->add_radiobutton_to_group(iso_isotype,"Hidden");
    glui_bounds->add_radiobutton_to_group(iso_isotype,"Solid");
    glui_bounds->add_radiobutton_to_group(iso_isotype,"Outline");
    glui_bounds->add_radiobutton_to_group(iso_isotype,"Points");
    
  }

  /* Particle File Bounds  */

//SVEXTERN part5prop *part5propinfo;
//SVEXTERN int npart5prop;

  if(npart_files>0&&nevac!=npart_files){
    glui_active=1;
    panel_part = glui_bounds->add_rollout("Particle",false);
    if(npart5prop>0){
      ipart5prop=0;
      ipart5prop_old=0;
      part5_rlist = glui_bounds->add_radiogroup_to_panel(panel_part,&ipart5prop,FILETYPEINDEX,PART_CB);

      for(i=0;i<npart5prop;i++){
        part5prop *partpropi;

        partpropi = part5propinfo + i;
        glui_bounds->add_radiobutton_to_group(part5_rlist,partpropi->label->shortlabel);
      }

      glui_bounds->add_column_to_panel(panel_part,false);

      {
        part5prop *propi;

        propi = part5propinfo;
        setpartmin=propi->setvalmin;
        setpartmax=propi->setvalmax;
        setpartchopmin=propi->setchopmin;
        setpartchopmax=propi->setchopmax;
      }

      PART_CB(FILETYPEINDEX);
      PART_CB(SETVALMIN);
      PART_CB(SETVALMAX);
    }


    {
      char boundmenulabel[100];
      GLUI_EditText **editcon;

#ifdef _DEBUG
      editcon=&con_part_chopmin;
      editcon=NULL;
#else
      editcon=NULL;
#endif

      strcpy(boundmenulabel,"Reload Particle File(s)");
      boundmenu(NULL,NULL,panel_part,boundmenulabel,
        &con_part_min,&con_part_max,&con_part_setmin,&con_part_setmax,
        &con_part_setchopmin, &con_part_setchopmax,
        editcon, &con_part_chopmax,
        &setpartmin,&setpartmax,&partmin,&partmax,
        &setpartchopmin,&setpartchopmax,&partchopmin,&partchopmax,
        DONT_UPDATEBOUNDS,DONT_TRUNCATEBOUNDS,
        PART_CB);
        if(partinfo!=NULL&&partinfo->version==0){
          SPINNER_partpointstep=glui_bounds->add_spinner_to_panel(panel_part,"Point Skip",GLUI_SPINNER_INT,
            &partpointskip,FRAMELOADING,PART_CB);
          SPINNER_partpointstep->set_int_limits(0,100);
        }
        SPINNER_partpointsize=glui_bounds->add_spinner_to_panel(panel_part,"Particle size",GLUI_SPINNER_FLOAT,
          &partpointsize);
        SPINNER_partpointsize->set_float_limits(1.0,10.0);
        SPINNER_streaklinewidth=glui_bounds->add_spinner_to_panel(panel_part,"Streak line width",GLUI_SPINNER_FLOAT,
          &streaklinewidth);
        SPINNER_streaklinewidth->set_float_limits(1.0,10.0);

        SPINNER_partstreaklength=glui_bounds->add_spinner_to_panel(panel_part,"Streak Length (s)",GLUI_SPINNER_FLOAT,
          &float_streak5value,STREAKLENGTH,PART_CB);
        SPINNER_partstreaklength->set_float_limits(0.0,10.0);
        
        showtracer_checkbox=glui_bounds->add_checkbox_to_panel(panel_part,"Always show tracers",&show_tracers_always,
          TRACERS,PART_CB);
    }
  }

  if(nevac>0){
    glui_active=1;
    panel_evac = glui_bounds->add_rollout("Evacuation",false);
    glui_bounds->add_checkbox_to_panel(panel_evac,"Select Avatar",&select_avatar);
   // glui_bounds->add_checkbox_to_panel(panel_evac,"View from selected Avatar",&view_from_selected_avatar);
  }

  /* Plot3D file bounds */

  if(nplot3d_files>0){
    glui_active=1;
    panel_plot3d = glui_bounds->add_rollout("Plot3D",false);

    p3rlist = glui_bounds->add_radiogroup_to_panel(panel_plot3d,&list_p3_index,FILETYPEINDEX,PLOT3D_CB);
    for(i=0;i<mxplot3dvars;i++){
      glui_bounds->add_radiobutton_to_group(p3rlist,plot3dinfo[0].label[i].shortlabel);
    }
    CHECKBOX_unload_qdata=glui_bounds->add_checkbox_to_panel(panel_plot3d,"Free memory after coloring data",
      &unload_qdata,UNLOAD_QDATA,PLOT3D_CB);

    panel_pan3 = glui_bounds->add_panel_to_panel(panel_plot3d,"",GLUI_PANEL_NONE);
    panel_contours = glui_bounds->add_panel_to_panel(panel_pan3,"Contours");
    plot3d_display=glui_bounds->add_radiogroup_to_panel(panel_contours,&p3cont2d,UPDATEPLOT,PLOT3D_CB);
    glui_bounds->add_radiobutton_to_group(plot3d_display,"Shaded Contours");
    glui_bounds->add_radiobutton_to_group(plot3d_display,"Stepped Contours");
    glui_bounds->add_radiobutton_to_group(plot3d_display,"Line Contours");
    
    //glui_bounds->add_column_to_panel(panel_pan3);
    panel_vector = glui_bounds->add_panel_to_panel(panel_pan3,"Vector");
    glui_bounds->add_checkbox_to_panel(panel_vector,"Show Vectors",&visVector,UPDATEPLOT,PLOT3D_CB);
    SPINNER_plot3d_vectorpointsize=glui_bounds->add_spinner_to_panel(panel_vector,"Point size",GLUI_SPINNER_FLOAT,&vectorpointsize,UPDATE_VECTOR,PLOT3D_CB);
    SPINNER_plot3d_vectorpointsize->set_float_limits(1.0,10.0);
    SPINNER_plot3d_vectorlinewidth=glui_bounds->add_spinner_to_panel(panel_vector,"Line width",GLUI_SPINNER_FLOAT,&vectorlinewidth,UPDATE_VECTOR,PLOT3D_CB);
    SPINNER_plot3d_vectorlinewidth->set_float_limits(1.0,10.0);

    panel_isosurface = glui_bounds->add_panel_to_panel(panel_plot3d,"Isosurface");
    panel_pan1 = glui_bounds->add_panel_to_panel(panel_isosurface,"",GLUI_PANEL_NONE);

    glui_bounds->add_checkbox_to_panel(panel_pan1,"Show Isosurface",&visiso,PLOTISO,PLOT3D_CB);
    SPINNER_plot3dpointsize=glui_bounds->add_spinner_to_panel(panel_pan1,"Point size",GLUI_SPINNER_FLOAT,
      &plot3dpointsize);
    SPINNER_plot3dpointsize->set_float_limits(1.0,10.0);

    SPINNER_plot3dlinewidth=glui_bounds->add_spinner_to_panel(panel_pan1,"Line width",GLUI_SPINNER_FLOAT,
      &plot3dlinewidth);
    SPINNER_plot3dlinewidth->set_float_limits(1.0,10.0);
//    glui_bounds->add_column_to_panel(panel_isosurface);
    panel_pan2 = glui_bounds->add_panel_to_panel(panel_isosurface,"",GLUI_PANEL_NONE);
    plot3d_isotype=glui_bounds->add_radiogroup_to_panel(panel_pan2,&p3dsurfacetype,PLOTISOTYPE,PLOT3D_CB);
    plot3d_iso_hidden=glui_bounds->add_radiobutton_to_group(plot3d_isotype,"Hidden");
    glui_bounds->add_radiobutton_to_group(plot3d_isotype,"Solid");
    glui_bounds->add_radiobutton_to_group(plot3d_isotype,"Outline");
    glui_bounds->add_radiobutton_to_group(plot3d_isotype,"Points");
    plot3d_iso_hidden->disable();


    p3min_temp=p3min[0];
    p3max_temp=p3max[0];
    p3chopmin_temp=p3chopmin[0];
    p3chopmax_temp=p3chopmax[0];
    glui_bounds->add_column_to_panel(panel_plot3d,false);
    boundmenu(NULL,NULL,panel_plot3d,"Reload Plot3D File(s)",
      &con_p3_min,&con_p3_max,&con_p3_setmin,&con_p3_setmax,
      &con_p3_setchopmin, &con_p3_setchopmax,
      &con_p3_chopmin, &con_p3_chopmax,
      &setp3min_temp,&setp3max_temp,&p3min_temp,&p3max_temp,
      &setp3chopmin_temp, &setp3chopmax_temp,&p3chopmin_temp,&p3chopmax_temp,
      DONT_UPDATEBOUNDS,TRUNCATEBOUNDS,
      PLOT3D_CB);
    PLOT3D_CB(UNLOAD_QDATA);
  }

  /*  Slice File Bounds   */

  if(nslice_files>0){
    int index;

    glui_active=1;
    panel_slice = glui_bounds->add_rollout("Slice",false);

    slice_rlist = glui_bounds->add_radiogroup_to_panel(panel_slice,&list_slice_index,FILETYPEINDEX,Slice_CB);

    index=0;
    for(i=0;i<nslice_files;i++){
      if(sliceinfo[i].firstshort==1){
        GLUI_RadioButton *RADIO_slicetype;

        RADIO_slicetype=glui_bounds->add_radiobutton_to_group(slice_rlist,sliceinfo[i].label.shortlabel);
        if(strcmp(sliceinfo[i].label.shortlabel,"Fire line")==0){
          RADIO_slicetype->disable();
          fire_line_index=index;
        }
        index++;
        firstcount++;
      }
    }

    glui_bounds->add_column_to_panel(panel_slice,false);

    boundmenu(&rollout_slice_bound,&rollout_slice_chop,panel_slice,"Reload Slice File(s)",
      &con_slice_min,&con_slice_max,&con_slice_setmin,&con_slice_setmax,
      &con_slice_setchopmin, &con_slice_setchopmax,
      &con_slice_chopmin, &con_slice_chopmax,
      &setslicemin,&setslicemax,&slicemin,&slicemax,
      &setslicechopmin, &setslicechopmax,
      &slicechopmin, &slicechopmax,
      UPDATEBOUNDS,DONT_TRUNCATEBOUNDS,
      Slice_CB);
#ifdef pp_TRANSFORM
    panel_slice_transform=glui_bounds->add_rollout_to_panel(panel_slice,"Transform data",false);
    glui_bounds->add_button_to_panel(panel_slice_transform,"Transform data",TRANSFORM_SLICE,Slice_CB);
    glui_bounds->add_button_to_panel(panel_slice_transform,"Reset data",RESET_SLICE,Slice_CB);
#endif
    SPINNER_transparentlevel=glui_bounds->add_spinner_to_panel(panel_slice,"Transparent level",GLUI_SPINNER_FLOAT,&transparentlevel,TRANSPARENTLEVEL,Slice_CB);
    SPINNER_transparentlevel->set_float_limits(0.0,1.0);
    CHECKBOX_average_slice=glui_bounds->add_checkbox_to_panel(panel_slice,"Averaged Slice Data",&slice_average_flag,AVERAGE_DATA,Slice_CB);
    //CHECKBOX_turb_slice=glui_bounds->add_checkbox_to_panel(panel_slice,"Turbulence Resolution",&slice_turbprop_flag,TURB_DATA,Slice_CB);
    SPINNER_sliceaverage=glui_bounds->add_spinner_to_panel(panel_slice,"Time Interval",GLUI_SPINNER_FLOAT,&slice_average_interval);
    {
      float tttmax=120.0;

      if(view_tstop>tttmax)tttmax=view_tstop;
      SPINNER_sliceaverage->set_float_limits(0.0,tttmax);
    }

    panel_slice_vector=glui_bounds->add_panel_to_panel(panel_slice,"Vector");
    SPINNER_vectorpointsize=glui_bounds->add_spinner_to_panel(panel_slice_vector,"Point size",GLUI_SPINNER_FLOAT,
      &vectorpointsize,UPDATE_VECTOR,Slice_CB);
    SPINNER_vectorpointsize->set_float_limits(1.0,10.0);
    SPINNER_vectorlinewidth=glui_bounds->add_spinner_to_panel(panel_slice_vector,"Line width",GLUI_SPINNER_FLOAT,
      &vectorlinewidth,UPDATE_VECTOR,Slice_CB);
    SPINNER_vectorlinewidth->set_float_limits(1.0,10.0);

#ifdef pp_SLICECONTOURS

    panel_line_contour = glui_bounds->add_panel_to_panel(panel_slice,"Line Contours");
    slice_line_contour_min=0.0;
    slice_line_contour_max=1.0;
    SPINNER_line_contour_min=glui_bounds->add_spinner_to_panel(panel_line_contour,"Min",GLUI_SPINNER_FLOAT,
      &slice_line_contour_min,LINE_CONTOUR_VALUE,Slice_CB);
    SPINNER_line_contour_max=glui_bounds->add_spinner_to_panel(panel_line_contour,"Max",GLUI_SPINNER_FLOAT,
      &slice_line_contour_max,LINE_CONTOUR_VALUE,Slice_CB);
    slice_line_contour_num=1;
    SPINNER_line_contour_num=glui_bounds->add_spinner_to_panel(panel_line_contour,"Number of contours",GLUI_SPINNER_INT,
      &slice_line_contour_num,LINE_CONTOUR_VALUE,Slice_CB);
    BUTTON_update_line_contour=glui_bounds->add_button_to_panel(panel_line_contour,"Update Contours",UPDATE_LINE_CONTOUR_VALUE,Slice_CB);
    glui_bounds->add_checkbox_to_panel(panel_line_contour,"Show contours",&vis_slice_contours);
#endif

    if(n_embedded_meshes>0){
      CHECKBOX_skip_subslice=glui_bounds->add_checkbox_to_panel(panel_slice,"Skip coarse sub-slice",&skip_slice_in_embedded_mesh);
    }
    if(cellcenter_slice_active==1){
      CHECKBOX_cellcenter_slice_interp = glui_bounds->add_checkbox_to_panel(panel_slice,"Interpolate cell centered slices",
        &cellcenter_interp);
    }
    glui_bounds->add_checkbox_to_panel(panel_slice,"Output data to file",&output_slicedata);
    Slice_CB(FILETYPEINDEX);
  }

  glui_bounds->add_separator();
  panel_time = glui_bounds->add_rollout("Time",false);
  panel_time1a = glui_bounds->add_panel_to_panel(panel_time,"",false);
  SPINNER_timebounds=glui_bounds->add_spinner_to_panel(panel_time1a,"time:",GLUI_SPINNER_FLOAT,&glui_time);
  glui_bounds->add_column_to_panel(panel_time1a,false);
  SPINNER_timebounds->set_float_limits(0.0,3600.0*24);
  BUTTON_SETTIME=glui_bounds->add_button_to_panel(panel_time1a,"Set",SET_TIME,Time_CB);
  
  panel_time2 = glui_bounds->add_panel_to_panel(panel_time,"Data Loading",true);

  panel_time2a = glui_bounds->add_panel_to_panel(panel_time2,"",false);
  SPINNER_tload_begin=glui_bounds->add_spinner_to_panel(panel_time2a,"tmin",GLUI_SPINNER_FLOAT,&tload_begin,TBOUNDS,Time_CB);
  glui_bounds->add_column_to_panel(panel_time2a,false);
  CHECKBOX_use_tload_begin=glui_bounds->add_checkbox_to_panel(panel_time2a,"",&use_tload_begin,TBOUNDS_USE,Time_CB);

  panel_time2b = glui_bounds->add_panel_to_panel(panel_time2,"",false);
  SPINNER_tload_end=glui_bounds->add_spinner_to_panel(panel_time2b,"tmax",GLUI_SPINNER_FLOAT,&tload_end,TBOUNDS,Time_CB);
  glui_bounds->add_column_to_panel(panel_time2b,false);
  CHECKBOX_use_tload_end=glui_bounds->add_checkbox_to_panel(panel_time2b,"",&use_tload_end,TBOUNDS_USE,Time_CB);

  panel_time2c = glui_bounds->add_panel_to_panel(panel_time2,"",false);
  SPINNER_tload_skip=glui_bounds->add_spinner_to_panel(panel_time2c,"frame skip",GLUI_SPINNER_INT,&tload_skip,TBOUNDS,Time_CB);
  glui_bounds->add_column_to_panel(panel_time2c,false);
  CHECKBOX_use_tload_skip=glui_bounds->add_checkbox_to_panel(panel_time2c,"",&use_tload_skip,TBOUNDS_USE,Time_CB);
  SPINNER_tload_skip->set_int_limits(0,100);
  
  BUTTON_RELOAD=glui_bounds->add_button_to_panel(panel_time2,"Reload",RELOAD_DATA,Time_CB);

  Time_CB(TBOUNDS_USE);

  panel_colorbar = glui_bounds->add_rollout("Data Coloring",false);
  if(ncolorbars>0){
    selectedcolorbar_index2=-1;
    LISTBOX_colorbar2=glui_bounds->add_listbox_to_panel(panel_colorbar,"Colorbar:",&selectedcolorbar_index2,COLORBAR_LIST2,Slice_CB);

    for(i=0;i<ncolorbars;i++){
      colorbardata *cbi;

      cbi = colorbarinfo + i;
      cbi->label_ptr=cbi->label;
      LISTBOX_colorbar2->add_item(i,cbi->label_ptr);
    }
    LISTBOX_colorbar2->set_int_val(colorbartype);
  }
  CHECKBOX_axissmooth=glui_bounds->add_checkbox_to_panel(panel_colorbar,"Smooth colorbar labels",&axissmooth,COLORBAR_SMOOTH,Slice_CB);
  CHECKBOX_showcolorbarlines=glui_bounds->add_checkbox_to_panel(panel_colorbar,"Lines",&showcolorbarlines,COLORBAR_LINES,Slice_CB);
  CHECKBOX_extreme2=glui_bounds->add_checkbox_to_panel(panel_colorbar,"Highlight extreme data",&show_extremedata,
    COLORBAR_EXTREME2,Slice_CB);

#ifdef pp_COMPRESS
  if(smokezippath!=NULL&&(npatch_files>0||nsmoke3d_files>0||nslice_files>0)){
    glui_bounds->add_separator();
    rollout_compress=glui_bounds->add_rollout("Compress Files (Smokezip)",false);
    check_erase_all=glui_bounds->add_checkbox_to_panel(rollout_compress,"Erase compressed files",
      &erase_all,ERASE,Bound_CB);
    check_overwrite_all=glui_bounds->add_checkbox_to_panel(rollout_compress,"Overwrite compressed files",
      &overwrite_all,OVERWRITE,Bound_CB);
    check_compress_autoloaded=glui_bounds->add_checkbox_to_panel(rollout_compress,"Compress only autoloaded files",
      &compress_autoloaded,COMPRESS_AUTOLOADED,Bound_CB);
#ifdef _DEBUG
    check_multi_task=glui_bounds->add_checkbox_to_panel(rollout_compress,"Multi-thread Compression",
      &mt_compress);
#endif
    if(nslice_files>0){
      SPINNER_slicezipstep=glui_bounds->add_spinner_to_panel(rollout_compress,"Slice Frame Skip",GLUI_SPINNER_INT,&slicezipskip,
        FRAMELOADING,Slice_CB);
      SPINNER_slicezipstep->set_int_limits(0,100);
    }
    if(niso_files>0){
      SPINNER_isozipstep=glui_bounds->add_spinner_to_panel(rollout_compress,"ISO Frame Skip",GLUI_SPINNER_INT,&isozipskip,
        FRAMELOADING,Iso_CB);
      SPINNER_isozipstep->set_int_limits(0,100);
    }
    if(nsmoke3d_files>0){
      SPINNER_smoke3dzipstep=glui_bounds->add_spinner_to_panel(rollout_compress,"3D Smoke Frame Skip",GLUI_SPINNER_INT,&smoke3dzipskip,
        FRAMELOADING,Smoke3D_CB);
      SPINNER_smoke3dzipstep->set_int_limits(0,100);
    }
    if(npatch_files>0){
      SPINNER_boundzipstep=glui_bounds->add_spinner_to_panel(rollout_compress,"Boundary File Frame Skip",
        GLUI_SPINNER_INT,&boundzipskip,FRAMELOADING,Bound_CB);
      SPINNER_boundzipstep->set_int_limits(0,100);
    }
    BUTTON_compress=glui_bounds->add_button_to_panel(rollout_compress,"Run smokezip",COMPRESS_FILES,Bound_CB);
  }
#endif

  rollout_SCRIPT = glui_bounds->add_rollout("Scripts/Config",false);

  panel_script1 = glui_bounds->add_panel_to_panel(rollout_SCRIPT,"Script files",true);

  panel_script1a = glui_bounds->add_panel_to_panel(panel_script1,"",false);
  BUTTON_script_start=glui_bounds->add_button_to_panel(panel_script1a,"Start Recording",SCRIPT_START,Script_CB);
  glui_bounds->add_column_to_panel(panel_script1a,false);
  BUTTON_script_stop=glui_bounds->add_button_to_panel(panel_script1a,"Stop Recording",SCRIPT_STOP,Script_CB);
  BUTTON_script_stop->disable();

  panel_script1b = glui_bounds->add_panel_to_panel(panel_script1,"",false);
  BUTTON_script_runscript=glui_bounds->add_button_to_panel(panel_script1b,"Run script",SCRIPT_RUNSCRIPT,Script_CB);
  glui_bounds->add_column_to_panel(panel_script1b,false);
  LIST_scriptlist = glui_bounds->add_listbox_to_panel(panel_script1b,"Select:",&script_index,SCRIPT_LIST,Script_CB);
    {
      scriptfiledata *scriptfile;

        for(scriptfile=first_scriptfile.next;scriptfile->next!=NULL;scriptfile=scriptfile->next){
          char *file;
          int len;

          file=scriptfile->file;
          if(file==NULL)continue;
          if(file_exist(file)==0)continue;
          len = strlen(file);
          if(len<=0)continue;

          LIST_scriptlist->add_item(scriptfile->id,file);
        }
        Script_CB(SCRIPT_LIST);
    }

  panel_script3 = glui_bounds->add_panel_to_panel(rollout_SCRIPT,"Render images",true);
  EDIT_renderdir=glui_bounds->add_edittext_to_panel(panel_script3,"render directory:",
    GLUI_EDITTEXT_TEXT,script_renderdir,SCRIPT_RENDER_DIR,Script_CB);
  EDIT_renderdir->set_w(260);
  panel_script1c = glui_bounds->add_panel_to_panel(panel_script3,"",false);
  BUTTON_script_render=glui_bounds->add_button_to_panel(panel_script1c,"Render",SCRIPT_RENDER,Script_CB);
  glui_bounds->add_column_to_panel(panel_script1c,false);
  EDIT_rendersuffix=glui_bounds->add_edittext_to_panel(panel_script1c,"suffix:",
    GLUI_EDITTEXT_TEXT,script_renderfilesuffix,SCRIPT_RENDER_SUFFIX,Script_CB);
  EDIT_rendersuffix->set_w(130);
  Script_CB(SCRIPT_RENDER_SUFFIX);

  panel_script2 = glui_bounds->add_panel_to_panel(rollout_SCRIPT,"Config files",true);

  panel_script2a = glui_bounds->add_panel_to_panel(panel_script2,"",false);
  EDIT_ini=glui_bounds->add_edittext_to_panel(panel_script2a,"suffix:",GLUI_EDITTEXT_TEXT,script_inifile_suffix,SCRIPT_EDIT_INI,Script_CB);
  glui_bounds->add_column_to_panel(panel_script2a,false);
  BUTTON_script_setsuffix=glui_bounds->add_button_to_panel(panel_script2a,"Set",SCRIPT_SETSUFFIX,Script_CB);
  glui_bounds->add_column_to_panel(panel_script2a,false);
  BUTTON_script_saveini=glui_bounds->add_button_to_panel(panel_script2a,"Save: ",SCRIPT_SAVEINI,Script_CB);
  Script_CB(SCRIPT_EDIT_INI);

  panel_script2b = glui_bounds->add_panel_to_panel(panel_script2,"",false);
  ini_index=-2;
  LIST_ini_list = glui_bounds->add_listbox_to_panel(panel_script2b,"Select:",&ini_index,SCRIPT_LISTINI,Script_CB);
  if(file_exist(caseinifilename)==1){
    ini_index=-1;
    LIST_ini_list->add_item(-1,caseinifilename);
  }
  {
    inifiledata *inifile;

    for(inifile=first_inifile.next;inifile->next!=NULL;inifile=inifile->next){
      if(inifile->file!=NULL&&file_exist(inifile->file)==1){
        if(ini_index==-2)ini_index=inifile->id;
        LIST_ini_list->add_item(inifile->id,inifile->file);
      }
    }
  }
  glui_bounds->add_column_to_panel(panel_script2b,false);
  BUTTON_ini_load=glui_bounds->add_button_to_panel(panel_script2b,"Load",SCRIPT_LOADINI,Script_CB);

  rollout_AUTOLOAD = glui_bounds->add_rollout("Auto Load",false);
  glui_bounds->add_checkbox_to_panel(rollout_AUTOLOAD,"Auto Load at Startup",
    &loadfiles_at_startup,STARTUP,Bound_CB);
  glui_bounds->add_button_to_panel(rollout_AUTOLOAD,"Save Auto Load File List",SAVE_FILE_LIST,Bound_CB);
  glui_bounds->add_button_to_panel(rollout_AUTOLOAD,"Auto Load Now",LOAD_FILES,Bound_CB);

  glui_bounds->add_button("Save Settings",SAVE_SETTINGS,BUTTON_hide_CB);

  glui_bounds->add_button("Close",CLOSE,BUTTON_hide_CB);

  glui_bounds->set_main_gfx_window( main_window );
}

/* ------------------ compress_onoff ------------------------ */

extern "C" void compress_onoff(int flag){
  if(flag==0){
    if(BUTTON_compress!=NULL)BUTTON_compress->disable();
    if(check_overwrite_all!=NULL)check_overwrite_all->disable();
    if(check_erase_all!=NULL)check_erase_all->disable();
    if(check_multi_task!=NULL)check_multi_task->disable();
  }
  else{
    if(BUTTON_compress!=NULL)BUTTON_compress->enable();
    if(check_overwrite_all!=NULL)check_overwrite_all->enable();
    if(check_erase_all!=NULL)check_erase_all->enable();
    if(check_multi_task!=NULL)check_multi_task->enable();
  }
}

/* ------------------ boundmenu ------------------------ */

void boundmenu(GLUI_Rollout **bound_rollout,GLUI_Rollout **chop_rollout, GLUI_Panel *panel, char *button_title,
          GLUI_EditText **con_min,GLUI_EditText **con_max,
          GLUI_RadioGroup **con_setmin,GLUI_RadioGroup **con_setmax,
          GLUI_Checkbox **con_setchopmin, GLUI_Checkbox **con_setchopmax,
          GLUI_EditText **con_chopmin, GLUI_EditText **con_chopmax,

          int *setminval, int *setmaxval,
          float *minval, float *maxval,
          int *setchopminval, int *setchopmaxval,
          float *chopminval, float *chopmaxval,
          int updatebounds,
          int truncatebounds,
          GLUI_Update_CB FILE_CB){

  GLUI_Panel *panel_a,*panel_b,*panel_c;
  GLUI_Rollout *panel_e=NULL,*panel_g=NULL;
  GLUI_Panel *panel_f=NULL;

  panel_g = glui_bounds->add_rollout_to_panel(panel,"Bound Data",false);
  if(bound_rollout!=NULL)*bound_rollout=panel_g;
  panel_a = glui_bounds->add_panel_to_panel(panel_g,"",GLUI_PANEL_NONE);

  *con_min = glui_bounds->add_edittext_to_panel(panel_a,"",GLUI_EDITTEXT_FLOAT,minval,VALMIN,FILE_CB);
  if(*setminval==0){
    (*con_min)->disable();
  }
  glui_bounds->add_column_to_panel(panel_a,false);
  *con_setmin = glui_bounds->add_radiogroup_to_panel(panel_a,setminval,SETVALMIN,FILE_CB);
  glui_bounds->add_radiobutton_to_group(*con_setmin,"percentile min");
  glui_bounds->add_radiobutton_to_group(*con_setmin,"set min");
  glui_bounds->add_radiobutton_to_group(*con_setmin,"global min");
  panel_b = glui_bounds->add_panel_to_panel(panel_g,"",GLUI_PANEL_NONE);

  *con_max = glui_bounds->add_edittext_to_panel(panel_b,"",GLUI_EDITTEXT_FLOAT,maxval,VALMAX,FILE_CB);
  if(*setminval==0){
    (*con_max)->disable();
  }
  glui_bounds->add_column_to_panel(panel_b,false);
    *con_setmax = glui_bounds->add_radiogroup_to_panel(panel_b,setmaxval,SETVALMAX,FILE_CB);
    glui_bounds->add_radiobutton_to_group(*con_setmax,"percentile max");
    glui_bounds->add_radiobutton_to_group(*con_setmax,"set max");
    glui_bounds->add_radiobutton_to_group(*con_setmax,"global max");

  panel_c = glui_bounds->add_panel_to_panel(panel_g,"",GLUI_PANEL_NONE);
  
  if(updatebounds==UPDATEBOUNDS){
    glui_bounds->add_button_to_panel(panel_c,"Update",FILEUPDATE,FILE_CB);
  }
  else{
    glui_bounds->add_button_to_panel(panel_c,button_title,FILERELOAD,FILE_CB);
  }           

  if(con_chopmin!=NULL&&con_chopmax!=NULL&&con_setchopmin!=NULL&&con_setchopmax!=NULL){
    panel_e = glui_bounds->add_rollout_to_panel(panel,"Truncate Data",false);
    if(chop_rollout!=NULL)*chop_rollout=panel_e;
    panel_f = glui_bounds->add_panel_to_panel(panel_e,"",GLUI_PANEL_NONE);
    *con_chopmin = glui_bounds->add_edittext_to_panel(panel_f,"",GLUI_EDITTEXT_FLOAT,chopminval,CHOPVALMIN,FILE_CB);
    *con_chopmax = glui_bounds->add_edittext_to_panel(panel_f,"",GLUI_EDITTEXT_FLOAT,chopmaxval,CHOPVALMAX,FILE_CB);
    glui_bounds->add_column_to_panel(panel_f,false);
    *con_setchopmin=glui_bounds->add_checkbox_to_panel(panel_f,"Below",setchopminval,SETCHOPMINVAL,FILE_CB);
    *con_setchopmax=glui_bounds->add_checkbox_to_panel(panel_f,"Above",setchopmaxval,SETCHOPMAXVAL,FILE_CB);
//  FILE_CB(SETCHOPMINVAL);
//ERROR    FILE_CB(SETCHOPMAXVAL);
    glui_bounds->add_button_to_panel(panel_e,"Update",CHOPUPDATE,FILE_CB);
  }

//  panel_d = glui_bounds->add_panel_to_panel(panel,"",GLUI_PANEL_NONE);

}

/* ------------------ PLOT3D_CB ------------------------ */

void PLOT3D_CB(int var){
  int i;

  switch (var){
  case UNLOAD_QDATA:
    if(unload_qdata==1){
     panel_isosurface->disable();
    }
    else{
      int enable_isosurface;
   
      enable_isosurface=1;
      for(i=0;i<nmeshes;i++){
        mesh *meshi;
        plot3d *plot3di;

        meshi = meshinfo + i;
        if(meshi->plot3dfilenum==-1)continue;
        plot3di = plot3dinfo + meshi->plot3dfilenum;
        if(plot3di->loaded==0||plot3di->display==0)continue;
        if(meshi->qdata==NULL){
          enable_isosurface=0;
          break;
        }
      }
      if(enable_isosurface==1)panel_isosurface->enable();
      if(enable_isosurface==0)panel_isosurface->disable();
    }
    break;
  case UPDATE_VECTOR:
    if(SPINNER_vectorpointsize!=NULL&&SPINNER_vectorlinewidth!=NULL){
      SPINNER_vectorpointsize->set_float_val(vectorpointsize);
      SPINNER_vectorlinewidth->set_float_val(vectorlinewidth);
    }
    break;
  case CHOPUPDATE:
    updatechopcolors();
    break;
  case SETCHOPMINVAL:
    updatechopcolors();
    switch (setp3chopmin_temp){
      case 0:
      con_p3_chopmin->disable();
      break;
      case 1:
      con_p3_chopmin->enable();
      break;
    }
  case SETCHOPMAXVAL:
    updatechopcolors();
    switch (setp3chopmax_temp){
      case 0:
      con_p3_chopmax->disable();
      break;
      case 1:
      con_p3_chopmax->enable();
      break;
    }
    break;
  case CHOPVALMIN:
//      &setp3chopmin_temp, &setp3chopmax_temp,&p3chopmin_temp,&p3chopmax_temp,
    p3chopmin[list_p3_index]=p3chopmin_temp;
    setp3chopmin[list_p3_index]=setp3chopmin_temp;

    updatechopcolors();
    break;
  case CHOPVALMAX:
    p3chopmax[list_p3_index]=p3chopmax_temp;
    setp3chopmax[list_p3_index]=setp3chopmax_temp;

    updatechopcolors();
    break;
  case PLOTISO:
    visiso = 1 - visiso;
    handleiso();
    GLUTPOSTREDISPLAY
    break;
  case PLOTISOTYPE:
    updatemenu=1;
    break;
  case UPDATEPLOT:
    if(p3cont2d==LINE_CONTOURS){
      showcolorbarlines=1;
      p3cont2d=SHADED_CONTOURS;
    }
    else{
      showcolorbarlines=0;
    }
    updatecolors(-1);
    updatemenu=1;
    GLUTPOSTREDISPLAY
    break;
  case FILETYPEINDEX:
   p3min[list_p3_index_old]=p3min_temp;
   p3max[list_p3_index_old]=p3max_temp;
   setp3min[list_p3_index_old]=setp3min_temp;
   setp3max[list_p3_index_old]=setp3max_temp;
   p3chopmin[list_p3_index_old]=p3chopmin_temp;
   p3chopmax[list_p3_index_old]=p3chopmax_temp;
   setp3chopmin[list_p3_index_old]=setp3chopmin_temp;
   setp3chopmax[list_p3_index_old]=setp3chopmax_temp;


   p3min_temp=p3min[list_p3_index];
   p3max_temp=p3max[list_p3_index];
   setp3min_temp=setp3min[list_p3_index];
   setp3max_temp=setp3max[list_p3_index];
   p3chopmin_temp=p3chopmin[list_p3_index];
   p3chopmax_temp=p3chopmax[list_p3_index];
   setp3chopmin_temp=setp3chopmin[list_p3_index];
   setp3chopmax_temp=setp3chopmax[list_p3_index];


   con_p3_min->set_float_val(p3min_temp);
   con_p3_max->set_float_val(p3max_temp);
   con_p3_chopmin->set_float_val(p3chopmin_temp);
   con_p3_chopmax->set_float_val(p3chopmax_temp);


   list_p3_index_old=list_p3_index;
   con_p3_setmin->set_int_val(setp3min_temp);
   con_p3_setmax->set_int_val(setp3max_temp);
   con_p3_setchopmin->set_int_val(setp3chopmin_temp);
   con_p3_setchopmax->set_int_val(setp3chopmax_temp);
   PLOT3D_CB(SETCHOPMINVAL);
   PLOT3D_CB(SETCHOPMAXVAL);
   PLOT3D_CB(SETVALMIN);
   PLOT3D_CB(SETVALMAX);
   break;
  case SETVALMIN:
   switch (setp3min_temp){
    case PERCENTILE_MIN:
    case GLOBAL_MIN:
      con_p3_min->disable();
      break;
    case SET_MIN:
    case CHOP_MIN:
      con_p3_min->enable();
      break;
   }
   break;
  case SETVALMAX:
     switch (setp3max_temp){
      case PERCENTILE_MIN:
      case GLOBAL_MIN:
        con_p3_max->disable();
        break;
      case SET_MIN:
      case CHOP_MIN:
        con_p3_max->enable();
        break;
     }
   break;
  case FILEUPDATE:
   p3min[list_p3_index] =    p3min_temp;
   p3max[list_p3_index] =    p3max_temp;
   setp3min[list_p3_index] = setp3min_temp;
   setp3max[list_p3_index] = setp3max_temp;
   break;
  case FILERELOAD:
   PLOT3D_CB(FILEUPDATE);
   for(i=0;i<nplot3d_files;i++){
     if(plot3dinfo[i].loaded==0)continue;
     LoadPlot3dMenu(i);
   }
   updateglui();
   break;
  }
}

/* ------------------ updatetracers ------------------------ */

extern "C" void updatetracers(void){
  if(showtracer_checkbox==NULL)return;
  showtracer_checkbox->set_int_val(show_tracers_always);
}


/* ------------------ update_glui_isotype ------------------------ */

extern "C" void update_glui_isotype(void){
  iso_isotype->set_int_val(visAIso);
}


/* ------------------ update_glui_isotype ------------------------ */

extern "C" void update_glui_plot3dtype(void){
  plot3d_isotype->set_int_val(p3dsurfacetype);
}
/* ------------------ updatechar ------------------------ */

extern "C" void updatechar(void){
  if(showchar_checkbox==NULL)return;
  if(canshow_threshold==1){
    showchar_checkbox->enable();
  }
  else{
    showchar_checkbox->disable();
  }
  showchar_checkbox->set_int_val(vis_threshold);
  Bound_CB(SHOWCHAR);
}

/* ------------------ updateplot3dlistindex ------------------------ */

extern "C" void updateplot3dlistindex(void){
  int i;

  if(p3rlist==NULL)return;
  i = p3rlist->get_int_val();
  if(i!=plotn-1){
    p3min[i]=p3min_temp;
    p3max[i]=p3max_temp;
    setp3min[i]=setp3min_temp;
    setp3max[i]=setp3max_temp;

    p3chopmin[i]=p3chopmin_temp;
    p3chopmax[i]=p3chopmax_temp;
    setp3chopmin[i]=setp3chopmin_temp;
    setp3chopmax[i]=setp3chopmax_temp;
  }

  i=plotn-1;
  list_p3_index_old=i;
  if(i<0)i=0;
  if(i>mxplot3dvars-1)i=mxplot3dvars-1;
  p3rlist->set_int_val(i);
  p3min_temp = p3min[i];
  p3max_temp = p3max[i];
  setp3min_temp = setp3min[i];
  setp3max_temp = setp3max[i];

  p3chopmin_temp = p3chopmin[i];
  p3chopmax_temp = p3chopmax[i];
  setp3chopmin_temp = setp3chopmin[i];
  setp3chopmax_temp = setp3chopmax[i];
  if(nplot3d_files>0){
    PLOT3D_CB(SETVALMIN);
    PLOT3D_CB(SETVALMAX);
    PLOT3D_CB(SETCHOPMINVAL);
    PLOT3D_CB(SETCHOPMAXVAL);
  }
  updatechopcolors();
  updateglui();
}

void Iso_CB(int var){
  switch (var) {
  case FRAMELOADING:
    isoframestep=isoframeskip+1;
    isozipstep=isozipskip+1;
    updatemenu=1;
    break;
  }
}
void Smoke3D_CB(int var){
  switch (var) {
  case FRAMELOADING:
    smoke3dframestep=smoke3dframeskip+1;
    smoke3dzipstep=smoke3dzipskip+1;
    update_glui_smoke3dframestep();
    updatemenu=1;
    break;
  }
}

/* ------------------ add_scriptlist ------------------------ */

extern "C" void add_scriptlist(char *file, int id){
  if(file!=NULL&&strlen(file)>0&&LIST_scriptlist!=NULL){
    LIST_scriptlist->add_item(id,file);
  }
}

/* ------------------ glui_script_enable ------------------------ */

  void glui_script_enable(void){
    BUTTON_script_start->enable();
    BUTTON_script_stop->enable();
    BUTTON_script_runscript->enable();
    LIST_scriptlist->enable();
    BUTTON_script_saveini->enable();
    BUTTON_script_setsuffix->enable();
    EDIT_ini->enable();
  }

/* ------------------ glui_script_disable ------------------------ */

  void glui_script_disable(void){
    BUTTON_script_start->disable();
    BUTTON_script_stop->disable();
    BUTTON_script_runscript->disable();
    LIST_scriptlist->disable();
    BUTTON_script_saveini->disable();
    BUTTON_script_setsuffix->disable();
    EDIT_ini->disable();
  }

  /* ------------------ Script_CB ------------------------ */

  void Script_CB(int var){
    char label[1024];
    char *name;
    int id;
    int len,i;
    int set_renderlabel;

    switch (var){
    case SCRIPT_RENDER_DIR:
      strcpy(label,script_renderdir);
      trim(label);
      name = trim_front(label);
      set_renderlabel=0;
      if(name!=NULL&&strlen(name)!=strlen(script_renderdir)){
        strcpy(script_renderdir,name);
        set_renderlabel=1;
      }
      name=script_renderdir;
      len = strlen(script_renderdir);
      if(len==0)break;
      for(i=0;i<len;i++){
#ifdef WIN32
        if(name[i]=='/'){
          set_renderlabel=1;
          name[i]='\\';
        }
#else
        if(name[i]=='\\'){
          set_renderlabel=1;
          name[i]='/';
        }
#endif
      }
#ifdef WIN32
      if(name[len-1]!='\\'){
        set_renderlabel=1;
        strcat(name,dirseparator);        
      }
#else
      if(name[len-1]!='/'){
        set_renderlabel=1;
        strcat(name,dirseparator);        
      }
#endif
      if(set_renderlabel==1){
        EDIT_renderdir->set_text(script_renderdir);
      }
      break;
    case SCRIPT_RENDER:
      keyboard('r',0,0);
      break;
    case SCRIPT_RENDER_SUFFIX:
    {
      char *suffix;

      trim(script_renderfilesuffix);
      suffix = trim_front(script_renderfilesuffix);
      strcpy(script_renderfile,"");
      if(strlen(suffix)>0){
        strcpy(script_renderfile,fdsprefix);
        strcat(script_renderfile,"_");
        strcat(script_renderfile,suffix);
        strcpy(label,"Render: ");
        strcat(label,script_renderfile);
      }
      else{
        strcpy(label,"Render");
      }
      BUTTON_script_render->set_name(label);
    }
      break;
    case SCRIPT_START:
      BUTTON_script_start->disable();
      BUTTON_script_stop->enable();
      BUTTON_script_runscript->disable();
      EDIT_renderdir->disable();

      ScriptMenu(START_RECORDING_SCRIPT);
      break;
    case SCRIPT_STOP:
      ScriptMenu(STOP_RECORDING_SCRIPT);
      BUTTON_script_start->enable();
      BUTTON_script_stop->disable();
      BUTTON_script_runscript->enable();
      EDIT_renderdir->enable();
      break;
    case SCRIPT_RUNSCRIPT:
      name = 5+BUTTON_script_runscript->name;
      printf("running script: %s\n",name);
      ScriptMenu(LIST_scriptlist->get_int_val());
      break;
    case SCRIPT_LIST:
      id = LIST_scriptlist->get_int_val();
      name = get_scriptfilename(id);
      if(name!=NULL&&strlen(name)>0){
        strcpy(label,"Run: ");
        strcat(label,name);
        BUTTON_script_runscript->set_name(label);
      }
      break;
    case SCRIPT_SAVEINI:
      name = 5+BUTTON_script_saveini->name;
      if(strlen(name)>0){
        inifiledata *inifile;

        strcpy(scriptinifilename,name);
        inifile=insert_inifile(name);
        writeini(SCRIPT_INI);
        if(inifile!=NULL&&LIST_ini_list!=NULL){
          LIST_ini_list->add_item(inifile->id,inifile->file);
        }
      }
      break;
    case SCRIPT_LOADINI:
      {
        char *inifilename;

        id = LIST_ini_list->get_int_val();
        if(id==-1){
          readini(0);
          if(scriptoutstream!=NULL){
            fprintf(scriptoutstream,"LOADINIFILE\n");
            fprintf(scriptoutstream," %s\n",caseinifilename);
          }
        }
        else if(id>=0){
          inifilename = get_inifilename(id);
          if(inifilename==NULL||strlen(inifilename)<=0)break;
          scriptinifilename2=scriptinifilename;
          strcpy(scriptinifilename,inifilename);
          windowresized=0;
          readini(2);
          if(scriptoutstream!=NULL){
            fprintf(scriptoutstream,"LOADINIFILE\n");
            fprintf(scriptoutstream," %s\n",scriptinifilename);
          }
          scriptinifilename2=NULL;
        }
      }
      break;
    case SCRIPT_EDIT_INI:
      strcpy(label,"Save ");
      strcat(label,fdsprefix);
      trim(script_inifile_suffix);
      if(strlen(script_inifile_suffix)>0){
        strcat(label,"_");
        strcat(label,script_inifile_suffix);
      }
      strcat(label,".ini");
      BUTTON_script_saveini->set_name(label);
      break;
    case SCRIPT_SETSUFFIX:
      break;
    }
  }

/* ------------------ Bound_CB ------------------------ */

  void Bound_CB(int var){
  patch *pi;
  int i;

  switch (var) {
  case FRAMELOADING:
    boundframestep=boundframeskip+1;
    boundzipstep=boundzipskip+1;
    updatemenu=1;
    break;
  case SETCHOPMINVAL:
    switch (setpatchchopmin){
      case 0:
        con_bf_chopmin->disable();
      break;
      case 1:
        con_bf_chopmin->enable();
      break;
    }
  case SETCHOPMAXVAL:
    switch (setpatchchopmax){
      case 0:
        con_bf_chopmax->disable();
      break;
      case 1:
        con_bf_chopmax->enable();
      break;
    }
    break;
  case SHOWCHAR:
    if(showchar_checkbox!=NULL&&showonlychar_checkbox!=NULL){
      if(vis_threshold==1){
        showonlychar_checkbox->enable();
      }
      else{
        showonlychar_checkbox->disable();
      }
    }
    updatemenu=1;
    updatefacelists=1;
    break;
  case FILETYPEINDEX:
    local2globalpatchbounds(patchlabellist[list_patch_index_old]);
    global2localpatchbounds(patchlabellist[list_patch_index]);

    con_bf_min->set_float_val(patchmin);
    con_bf_max->set_float_val(patchmax);

    Bound_CB(SETVALMIN);
    Bound_CB(SETVALMAX);
    if(con_bf_setmin!=NULL)con_bf_setmin->set_int_val(setpatchmin);
    if(con_bf_setmax!=NULL)con_bf_setmax->set_int_val(setpatchmax);

    list_patch_index_old = list_patch_index;
    break;
  case SETVALMIN:
    switch (setpatchmin){
    case PERCENTILE_MIN:
    case GLOBAL_MIN:
      con_bf_min->disable();
      break;
    case SET_MIN:
      con_bf_min->enable();
      break;
    }
    Bound_CB(FILEUPDATE);
    break;
  case SETVALMAX:
    switch (setpatchmax){
    case PERCENTILE_MAX:
    case GLOBAL_MAX:
      con_bf_max->disable();
      break;
    case SET_MAX:
      con_bf_max->enable();
      break;
    }
    Bound_CB(FILEUPDATE);
  break;
  case FILEUPDATE:
    for(i=0;i<npatch_files;i++){
      if(strcmp(patchinfo[i].label.shortlabel,patchlabellist[list_patch_index])==0){
        patchinfo[i].valmin=patchmin;
        patchinfo[i].valmax=patchmax;
        patchinfo[i].setvalmin=setpatchmin;
        patchinfo[i].setvalmax=setpatchmax;
      }
    }
    break;
  case FILERELOAD:
    Bound_CB(FILEUPDATE);
    for(i=0;i<npatch_files;i++){
      pi=patchinfo+i;
      if(pi->loaded==0)continue;
      LoadPatchMenu(i);
    }
    con_bf_min->set_float_val(patchmin);
    con_bf_max->set_float_val(patchmax);
    break;
  case COMPRESS_FILES:
    compress_svzip();
    break;
  case COMPRESS_AUTOLOADED:
    updatemenu=1;
    break;
  case OVERWRITE:
    if(overwrite_all==1){
      check_erase_all->set_int_val(0);
    }
    updatemenu=1;
    break;
  case ERASE:
    if(erase_all==1){
      check_overwrite_all->set_int_val(0);
    }
    updatemenu=1;
    break;
  case STARTUP:
    BUTTON_hide_CB(SAVE_SETTINGS);
    break;
  case SAVE_FILE_LIST:
    set_3dsmoke_startup();
    BUTTON_hide_CB(SAVE_SETTINGS);
    break;
  case LOAD_FILES:
      load_startup_smoke();
    break;
  }
}
/* ------------------ updatepatchlistindex ------------------------ */

extern "C" void updatepatchlistindex(int patchfilenum){
  int i;
  for(i=0;i<npatch2;i++){
    if(strcmp(patchlabellist[i],patchinfo[patchfilenum].label.shortlabel)==0){
      bf_rlist->set_int_val(i);
      list_patch_index_old=list_patch_index;
      global2localpatchbounds(patchlabellist[i]);
      con_bf_setmin->set_int_val(setpatchmin);
      con_bf_setmax->set_int_val(setpatchmax);
      con_bf_min->set_float_val(patchmin);
      con_bf_max->set_float_val(patchmax);
      if(setpatchmin==SET_MIN){
        con_bf_min->enable();
      }
      else{
        con_bf_min->disable();
      }
      if(setpatchmax==SET_MAX){
        con_bf_max->enable();
      }
      else{
        con_bf_max->disable();
      }
      return;
    }
  }
}

/* ------------------ update_glui_streakvalue ------------------------ */

extern "C" void update_glui_cellcenter_interp(void){
  if(CHECKBOX_cellcenter_slice_interp!=NULL)CHECKBOX_cellcenter_slice_interp->set_int_val(cellcenter_interp);
  if(CHECKBOX_skip_subslice!=NULL)CHECKBOX_skip_subslice->set_int_val(skip_slice_in_embedded_mesh);
}

/* ------------------ update_glui_streakvalue ------------------------ */

extern "C" void update_glui_streakvalue(float rvalue){
  float_streak5value=rvalue;
  if(SPINNER_partstreaklength!=NULL)SPINNER_partstreaklength->set_float_val(rvalue);
}

/* ------------------ PART_CB ------------------------ */

void PART_CB(int var){
  part5prop *propi, *propi0;

  propi = part5propinfo + ipart5prop;
  propi0 = part5propinfo + ipart5prop_old;
  switch (var){
  case VALMIN:
    if(setpartmax==SET_MAX)propi->user_max=partmax;
    break;
  case VALMAX:
    if(setpartmin==SET_MIN)propi->user_min=partmin;
    break;
  case FILETYPEINDEX:

    // save data from controls

    propi0->setvalmin=setpartmin;
    propi0->setvalmax=setpartmax;
    if(setpartmin==SET_MIN){
      propi0->user_min=partmin;
    }
    if(setpartmax==SET_MAX){
      propi0->user_max=partmax;
    }
    propi0->setchopmin=setpartchopmin;
    propi0->setchopmax=setpartchopmax;
    propi0->chopmin=partchopmin;
    propi0->chopmax=partchopmax;

    // copy data to controls

    setpartmin=propi->setvalmin;
    setpartmax=propi->setvalmax;
    PART_CB(SETVALMIN);
    PART_CB(SETVALMAX);

    setpartchopmin=propi->setchopmin;
    setpartchopmax=propi->setchopmax;
    partchopmin=propi->chopmin;
    partchopmax=propi->chopmax;

    // update controls

    if(con_part_setmin!=NULL)con_part_setmin->set_int_val(setpartmin);
    if(con_part_setmax!=NULL)con_part_setmax->set_int_val(setpartmax);


    if(con_part_chopmin!=NULL)con_part_chopmin->set_float_val(partchopmin);
    if(con_part_chopmax!=NULL)con_part_chopmax->set_float_val(partchopmax);
    if(con_part_setchopmin!=NULL)con_part_setchopmin->set_int_val(setpartchopmin);
    if(con_part_setchopmax!=NULL)con_part_setchopmax->set_int_val(setpartchopmax);

    ipart5prop_old = ipart5prop;

    if(con_part_setchopmin!=NULL)PART_CB(SETCHOPMINVAL);
    if(con_part_setchopmax!=NULL)PART_CB(SETCHOPMAXVAL);

    break;
  case STREAKLENGTH:
    update_streakvalue(float_streak5value);
    if(float_streak5value==0.0){
      streak5show=0;
    }
    else{
      streak5show=1;
    }
    updatemenu=1;
    break;
  case TRACERS:
    updatemenu=1;
    break;
  case FRAMELOADING:
    partframestep=partframeskip+1;
    partpointstep=partpointskip+1;
    evacframestep=evacframeskip+1;
    evacframestep=evacframeskip+1;
    if(partpointstep!=partpointstep_old){
      update_all_partvis2();
      partpointstep_old=partpointstep;
    }

    updatemenu=1;
    break;
  case CHOPUPDATE:
    updatechopcolors();
    break;
  case SETCHOPMINVAL:
    updatechopcolors();
    switch (setpartchopmin){
      case 0:
      con_part_chopmin->disable();
      break;
      case 1:
      con_part_chopmin->enable();
      break;
    }
  case SETCHOPMAXVAL:
    updatechopcolors();
    switch (setpartchopmax){
      case 0:
      con_part_chopmax->disable();
      break;
      case 1:
      con_part_chopmax->enable();
      break;
    }
    break;
  case CHOPVALMIN:
    if(con_part_min!=NULL)con_part_min->set_float_val(partchopmin);
//    SETslicemin(setslicemin,slicemin);
    updatechopcolors();
    break;
  case CHOPVALMAX:
    if(con_part_max!=NULL)con_part_max->set_float_val(partchopmax);
//    SETslicemax(setslicemax,slicemax);
    updatechopcolors();
    break;
  case SETVALMIN:
    if(setpartmin_old==SET_MIN){
      if(propi0!=NULL)propi0->user_min=partmin;
    }
    setpartmin_old=setpartmin;
    switch (setpartmin){
    case PERCENTILE_MIN:
      if(propi!=NULL)partmin=propi->percentile_min;
      if(con_part_min!=NULL)con_part_min->disable();
      break;
    case GLOBAL_MIN:
      if(propi!=NULL)partmin=propi->global_min;
      if(con_part_min!=NULL)con_part_min->disable();
      break;
    case SET_MIN:
      if(propi!=NULL)partmin=propi->user_min;
      if(con_part_min!=NULL)con_part_min->enable();
      break;
    }
    if(propi!=NULL)propi->valmin=partmin;
    if(con_part_min!=NULL)con_part_min->set_float_val(partmin);
   break;
  case SETVALMAX:
    if(setpartmax_old==SET_MAX){
      if(propi0!=NULL)propi0->user_max=partmax;
    }
    setpartmax_old=setpartmax;
    switch (setpartmax){
    case PERCENTILE_MAX:
      if(propi!=NULL)partmax=propi->percentile_max;
      if(con_part_max!=NULL)con_part_max->disable();
      break;
    case GLOBAL_MAX:
      if(propi!=NULL)partmax=propi->global_max;
      if(con_part_max!=NULL)con_part_max->disable();
      break;
    case SET_MAX:
      if(propi!=NULL)partmax=propi->user_max;
      if(con_part_max!=NULL)con_part_max->enable();
      break;
    }
    if(propi!=NULL)propi->valmax=partmax;
    if(con_part_max!=NULL)con_part_max->set_float_val(partmax);
   break;
  case FILERELOAD:
    {
      int prop_index_SAVE;

     prop_index_SAVE=prop_index;
     PART_CB(FILETYPEINDEX);
     if(con_part_min!=NULL&&setpartmin==SET_MIN)PART_CB(SETVALMIN);
     if(con_part_max!=NULL&&setpartmax==SET_MAX)PART_CB(SETVALMAX);
     ParticleMenu(0);
     updateglui();
     ParticlePropShowMenu(prop_index_SAVE);
    }
   break;
  }
}

/* ------------------ Time_CB ------------------------ */

void Time_CB(int var){

  updatemenu=1;
  switch (var){
  case SET_TIME:
    settimeval(glui_time);
    break;
  case TBOUNDS:
    if(use_tload_begin==1||use_tload_end==1||use_tload_skip==1){
      update_tbounds();
    }
    break;
  case TBOUNDS_USE:
    if(use_tload_begin==1){
      SPINNER_tload_begin->enable();
    }
    else{
      SPINNER_tload_begin->disable();
    }
    if(use_tload_end==1){
      SPINNER_tload_end->enable();
    }
    else{
      SPINNER_tload_end->disable();
    }
    if(use_tload_skip==1){
      SPINNER_tload_skip->enable();
    }
    else{
      SPINNER_tload_skip->disable();
    }
    update_tbounds();
    break;
  case RELOAD_DATA:
    ReloadMenu(0);
    break;
  }
}

/* ------------------ SLICE_CB ------------------------ */

extern "C" void Slice_CB(int var){
  int error,i;
  int ii;
  slice *sd;

  updatemenu=1;
  if(var==COLORBAR_EXTREME2){
    update_extreme();
    return;
  }
  if(var==COLORBAR_LIST2){
      selectedcolorbar_index=LISTBOX_colorbar2->get_int_val();
      update_colorbar_list();
      ColorBarMenu(selectedcolorbar_index);
      colorbar_global2local();
  }
  if(var==COLORBAR_LINES){
    if(showcolorbarlines==1)p3cont2d=SHADED_CONTOURS;
    update_plot3d_display();
    updatecolors(-1);
    return;
  }
  if(var==COLORBAR_SMOOTH){
    updatemenu=1;
    return;
  }
  switch (var){
    case TRANSPARENTLEVEL:
      updatecolors(-1);
      break;
#ifdef pp_SLICECONTOURS
    case LINE_CONTOUR_VALUE:
      if(slice_line_contour_num<1){
        slice_line_contour_num=1;
        SPINNER_line_contour_num->set_int_val(slice_line_contour_num);
      }
      if(slice_line_contour_num==1&&slice_line_contour_min!=slice_line_contour_max){
        slice_line_contour_max=slice_line_contour_min;
        SPINNER_line_contour_max->set_float_val(slice_line_contour_max);
      }
      slicebounds[list_slice_index].line_contour_min=slice_line_contour_min;
      slicebounds[list_slice_index].line_contour_max=slice_line_contour_max;
      slicebounds[list_slice_index].line_contour_num=slice_line_contour_num;
      break;
    case UPDATE_LINE_CONTOUR_VALUE:
      update_slice_contours(list_slice_index,slice_line_contour_min, slice_line_contour_max,slice_line_contour_num);
      break;
#endif
  case AVERAGE_DATA:
    if(slice_average_flag==1&&slice_turbprop_flag==1){
      slice_turbprop_flag=0;
      if(CHECKBOX_turb_slice!=NULL)CHECKBOX_turb_slice->set_int_val(slice_turbprop_flag);
    }
    break;
  case TURB_DATA:
    if(slice_average_flag==1&&slice_turbprop_flag==1){
      slice_average_flag=0;
      CHECKBOX_average_slice->set_int_val(slice_average_flag);
    }
    break;
  case UPDATE_VECTOR:
    if(SPINNER_plot3d_vectorpointsize!=NULL&&SPINNER_plot3d_vectorlinewidth!=NULL){
      SPINNER_plot3d_vectorpointsize->set_float_val(vectorpointsize);
      SPINNER_plot3d_vectorlinewidth->set_float_val(vectorlinewidth);
    }
    break;
  case FRAMELOADING:
    sliceframestep=sliceframeskip+1;
    slicezipstep=slicezipskip+1;
    updatemenu=1;
    break;
  case CHOPUPDATE:
    updatechopcolors();
    break;
  case SETCHOPMINVAL:
    updatechopcolors();
    SETslicemin(setslicemin,slicemin);
    switch (setslicechopmin){
      case 0:
      con_slice_chopmin->disable();
      break;
      case 1:
      con_slice_chopmin->enable();
      break;
    }
  case SETCHOPMAXVAL:
    updatechopcolors();
    SETslicemax(setslicemax,slicemax);
    switch (setslicechopmax){
      case 0:
      con_slice_chopmax->disable();
      break;
      case 1:
      con_slice_chopmax->enable();
      break;
    }
    break;
  case CHOPVALMIN:
    ASSERT(con_slice_min!=NULL);
    con_slice_min->set_float_val(slicemin);
    SETslicemin(setslicemin,slicemin);
    updatechopcolors();
    break;
  case CHOPVALMAX:
    ASSERT(con_slice_max!=NULL);
    con_slice_max->set_float_val(slicemax);
    SETslicemax(setslicemax,slicemax);
    updatechopcolors();
    break;
  case SETVALMIN:
    ASSERT(con_slice_min!=NULL);
    switch (setslicemin){
    case PERCENTILE_MIN:
    case GLOBAL_MIN:
      con_slice_min->disable();
      break;
    case SET_MIN:
      con_slice_min->enable();
      break;
    }
    if(con_slice_setmin!=NULL)con_slice_setmin->set_int_val(setslicemin);
    SETslicemin(setslicemin,slicemin);
    break;
  case SETVALMAX:
    ASSERT(con_slice_max!=NULL);
    switch (setslicemax){
    case PERCENTILE_MAX:
    case GLOBAL_MAX:
      con_slice_max->disable();
      break;
    case SET_MAX:
      con_slice_max->enable();
      break;
    }
    if(con_slice_setmax!=NULL)con_slice_setmax->set_int_val(setslicemax);
    SETslicemax(setslicemax,slicemax);
    break;
  case VALMIN:
    ASSERT(con_slice_min!=NULL);
    con_slice_min->set_float_val(slicemin);
    SETslicemin(setslicemin,slicemin);
    break;
  case VALMAX:
    ASSERT(con_slice_max!=NULL);
    con_slice_max->set_float_val(slicemax);
    SETslicemax(setslicemax,slicemax);
    break;
  case FILETYPEINDEX:
    if(slice_bounds_dialog==1&&list_slice_index==fire_line_index){
      slice_bounds_dialog=0;
      if(rollout_slice_bound!=NULL){
        rollout_slice_bound->close();
        rollout_slice_bound->disable();
      }
      if(rollout_slice_chop!=NULL){
        rollout_slice_chop->close();
        rollout_slice_chop->disable();
      }
    }
    if(slice_bounds_dialog==0&&list_slice_index!=fire_line_index){
      slice_bounds_dialog=1;
      if(rollout_slice_bound!=NULL){
        rollout_slice_bound->enable();
      }
      if(rollout_slice_chop!=NULL){
        rollout_slice_chop->enable();
      }
    }
    setslicebounds(list_slice_index);
    ASSERT(con_slice_min!=NULL);
    con_slice_min->set_float_val(slicemin);
    switch (setslicemin){
    case PERCENTILE_MIN:
    case GLOBAL_MIN:
      con_slice_min->disable();
      break;
    case SET_MIN:
      con_slice_min->enable();
      break;
    }
    ASSERT(con_slice_max!=NULL);
    con_slice_max->set_float_val(slicemax);
    switch (setslicemax){
    case PERCENTILE_MIN:
    case GLOBAL_MAX:
      con_slice_max->disable();
      break;
    case SET_MAX:
      con_slice_max->enable();
      break;
    }

    con_slice_setmin->set_int_val(setslicemin);
    con_slice_setmax->set_int_val(setslicemax);

    con_slice_chopmin->set_float_val(slicechopmin);
    con_slice_chopmax->set_float_val(slicechopmax);
    con_slice_setchopmin->set_int_val(setslicechopmin);
    con_slice_setchopmax->set_int_val(setslicechopmax);

    if(setslicechopmin==1){
      con_slice_chopmin->enable();
    }
    else{
      con_slice_chopmin->disable();
    }
    if(setslicechopmax==1){
      con_slice_chopmax->enable();
    }
    else{
      con_slice_chopmax->disable();
    }
#ifdef pp_SLICECONTOURS
    SPINNER_line_contour_min->set_float_val(slice_line_contour_min);
    SPINNER_line_contour_max->set_float_val(slice_line_contour_max);
    SPINNER_line_contour_num->set_int_val(slice_line_contour_num);
#endif
    break;
  case FILEUPDATE:
    for(ii=0;ii<nslice_loaded;ii++){
      i = slice_loaded_list[ii];
      sd = sliceinfo + i;
      if(sd->type!=islicetype)continue;
      readslice("",i, RESETBOUNDS, &error);
    }
    break;
  case FILERELOAD:
    Slice_CB(FILEUPDATE);
    if(slicefilenum>=0&&slicefilenum<nslice_files){
      LoadSliceMenu(slicefilenum);
    }
    else{
      LoadSliceMenu(0);
    }
    updateglui();
    break;
  }
}

/* ------------------ SETslicemin ------------------------ */

void SETslicemin(int setslicemin, float slicemin){
  slicebounds[list_slice_index].setvalmin=setslicemin;
  slicebounds[list_slice_index].valmin=slicemin;
  slicebounds[list_slice_index].setchopmin=setslicechopmin;
  slicebounds[list_slice_index].chopmin=slicechopmin;
}

/* ------------------ SETslicemax ------------------------ */

void SETslicemax(int setslicemax, float slicemax){
  slicebounds[list_slice_index].setvalmax=setslicemax;
  slicebounds[list_slice_index].valmax=slicemax;
  slicebounds[list_slice_index].setchopmax=setslicechopmax;
  slicebounds[list_slice_index].chopmax=slicechopmax;
}

/* ------------------ updateslicelist ------------------------ */

extern "C" void updateslicelist(int index){
  slice_rlist->set_int_val(index);
}

/* ------------------ updateslicelistindex ------------------------ */

extern "C" void updateslicelistindex(int sfn){
  int i;
  int slicefiletype;
  slice *sd;
  if(sfn<0){
    updateslicefilenum();
    sfn=slicefilenum;
  }
  sd = sliceinfo+sfn;
  slicefiletype=getsliceindex(sd);
  if(slicefiletype>=0&&slicefiletype<nslice2){
    i = slicefiletype;
    slice_rlist->set_int_val(i);
    setslicebounds(i);
    list_slice_index=i;
    Slice_CB(SETVALMIN);
    Slice_CB(SETVALMAX);
    Slice_CB(VALMIN);
    Slice_CB(VALMAX);
    Slice_CB(SETCHOPMINVAL);
    Slice_CB(SETCHOPMAXVAL);
    Slice_CB(CHOPVALMIN);
    Slice_CB(CHOPVALMAX);
  }
}

/* ------------------ updateglui ------------------------ */

void updateglui(void){
  GLUI_Master.sync_live_all();
}

/* ------------------ BUTTON_hide_CB ------------------------ */

void BUTTON_hide_CB(int var){
  switch (var){
  case CLOSE:
    glui_bounds->hide();
    showbounds=0;
    updatemenu=1;
    break;
  case SAVE_SETTINGS:
    writeini(LOCAL_INI);
    break;
  case COMPRESS_FILES:
    printf("compressing\n");
    break;
  }
}

/* ------------------ show_glui ------------------------ */

extern "C" void show_glui_bounds(void){
  int islice, ipatch;
  if(nslice_files>0){
    islice=slice_rlist->get_int_val();
    setslicebounds(islice);
    Slice_CB(SETVALMIN);
    Slice_CB(SETVALMAX);
    Slice_CB(VALMIN);
    Slice_CB(VALMAX);
  }
  if(npatch_files>0){
    ipatch=bf_rlist->get_int_val();
    global2localpatchbounds(patchlabellist[ipatch]);
    Bound_CB(SETVALMIN);
    Bound_CB(SETVALMAX);
  }
  if(npart_files>0&&npart_files!=nevac){
    PART_CB(SETVALMIN);
    PART_CB(SETVALMAX);
  }
  if(nplot3d_files>0){
    PLOT3D_CB(SETVALMIN);
    PLOT3D_CB(SETVALMAX);
  }

  if(nslice_files>0||npatch_files>0)updateglui();

  updatechar();
  glui_bounds->show();
}

/* ------------------ enable_boundary_glui ------------------------ */

extern "C" void enable_boundary_glui(void){
  rollout_BOUNDARY->enable();
}

/* ------------------ disable_boundary_glui ------------------------ */

extern "C" void disable_boundary_glui(void){
  rollout_BOUNDARY->disable();
}

/* ------------------ update_overwrite ------------------------ */

extern "C" void update_overwrite(void){
  if(check_overwrite_all!=NULL)check_overwrite_all->set_int_val(overwrite_all);
  if(check_compress_autoloaded!=NULL)check_compress_autoloaded->set_int_val(compress_autoloaded);
}

/* ------------------ hide_glui_bounds ------------------------ */

extern "C" void hide_glui_bounds(void){
  if(glui_bounds!=NULL)glui_bounds->hide();
  showbounds=0;
}

/* ------------------ update_plot3d_display ------------------------ */

extern "C" void update_plot3d_display(void){
  if(plot3d_display!=NULL)plot3d_display->set_int_val(p3cont2d);
}

/* ------------------ update2_glui_smoke3dframestep ------------------------ */

extern "C" void update2_glui_smoke3dframestep(void){
}

/* ------------------ open_smokeplane ------------------------ */

extern "C" void open_smokepanel(void){
  if(panel_smoke3d!=NULL){
    panel_smoke3d->open();
  }
}

/* ------------------ open_smokezippanel ------------------------ */

extern "C" void open_smokezippanel(void){
  if(rollout_compress!=NULL){
    rollout_compress->open();
  }
}

/* ------------------ updateGluiTimeBounds ------------------------ */

extern "C" void updateGluiTimeBounds(float time_min, float time_max){
  if(SPINNER_timebounds!=NULL){
    SPINNER_timebounds->set_float_limits(time_min,time_max);
  }
}

/* ------------------ update_tbounds ------------------------ */

void update_tbounds(void){
  settmin_p=use_tload_begin;
  settmax_p=use_tload_end;
  tmin_p=tload_begin;
  tmax_p=tload_end;

  settmin_s=use_tload_begin;
  settmax_s=use_tload_end;
  tmin_s=tload_begin;
  tmax_s=tload_end;

  settmin_i=use_tload_begin;
  settmax_i=use_tload_end;
  tmin_i=tload_begin;
  tmax_i=tload_end;

  settmin_s=use_tload_begin;
  settmax_s=use_tload_end;
  tmin_s=tload_begin;
  tmax_s=tload_end;

  settmin_b=use_tload_begin;
  settmax_b=use_tload_end;
  tmin_b=tload_begin;
  tmax_b=tload_end;

  if(use_tload_skip==1){
    smoke3dframeskip=tload_skip;
    boundframeskip=tload_skip;
    isoframeskip=tload_skip;
    partframeskip=tload_skip;
    evacframeskip=tload_skip;
    sliceframeskip=tload_skip;
  }
  else{
    smoke3dframeskip=0;
    boundframeskip=0;
    isoframeskip=0;
    partframeskip=0;
    evacframeskip=0;
    sliceframeskip=0;
  }

  Smoke3D_CB(FRAMELOADING);
  Bound_CB(FRAMELOADING);
  Iso_CB(FRAMELOADING);
  PART_CB(FRAMELOADING);
  Slice_CB(FRAMELOADING);
}

