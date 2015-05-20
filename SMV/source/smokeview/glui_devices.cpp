#define CPP
#include "options.h"

#include <stdio.h>
#include <string.h>
#include GLUT_H
#include <math.h>

#include "smokeviewvars.h"

#define DEVICE_sensorsize 20
#define SHOWDEVICEVALS 26
#define COLORDEVICEVALS 27
#define DEVICE_devicetypes 28
#define SAVE_SETTINGS 99
#define DEVICE_close 3
#define DEVICE_show_orientation 4

#define OPEN_UP 0
#define OPEN_DOWN 1
#define OPEN_FILEINDEX 2
#define OPEN_OPEN 3
#define OPEN_CANCEL 4
#define OPEN_FILTER 5
#define OPEN_APPLY_FILTER 6
#define OPEN_UPDATE_LIST 7

class CGluiOpen {
};

void Open_CB(int var);

int gluiopen_file_index=0;
int gluiopen_nfilelist=0;
char gluiopen_path_dir[1024];
filelistdata *gluiopen_filelist;
char gluiopen_filter[sizeof(GLUI_String)];
char gluiopen_filter2[sizeof(GLUI_String)];

GLUI *glui_device=NULL;

GLUI_Button *BUTTON_open_down=NULL ;
GLUI_Button *BUTTON_device_1=NULL;
GLUI_Button *BUTTON_device_2=NULL;

GLUI_Panel *PANEL_objects=NULL;
GLUI_Panel *PANEL_velocityvectors=NULL;
GLUI_Panel *PANEL_vectors=NULL;
GLUI_Panel *PANEL_arrow_base=NULL;
GLUI_Panel *PANEL_arrow_height=NULL;
GLUI_Panel *PANEL_devicevalues=NULL;
GLUI_Panel *PANEL_smvobjects=NULL;
GLUI_Panel *PANEL_devicevis=NULL;
GLUI_Panel *PANEL_label3=NULL;
GLUI_Rollout *ROLLOUT_arrow_dimensions=NULL;
GLUI_Panel *PANEL_vector_type=NULL;

GLUI_Listbox *LIST_open=NULL;

GLUI_EditText *EDIT_filter=NULL;

GLUI_Spinner *SPINNER_sensorrelsize=NULL;
GLUI_Spinner *SPINNER_orientation_scale=NULL;

GLUI_RadioGroup *RADIO_devicetypes=NULL;
GLUI_RadioGroup *RADIO_vectortype=NULL;

GLUI_Checkbox *CHECKBOX_device_1=NULL;
GLUI_Checkbox *CHECKBOX_device_2=NULL;
GLUI_Checkbox *CHECKBOX_device_3=NULL;
GLUI_Checkbox *CHECKBOX_device_4=NULL;
GLUI_Checkbox *CHECKBOX_device_orientation=NULL;

void Device_CB(int var);

/* ------------------ update_device_size ------------------------ */

extern "C" void update_device_size(void){
  if(sensorrelsize<sensorrelsizeMIN)sensorrelsize = sensorrelsizeMIN;
  if(SPINNER_sensorrelsize!=NULL)SPINNER_sensorrelsize->set_float_val(sensorrelsize);
}

/* ------------------ update_device_orientation ------------------------ */

extern "C" void update_device_orientation(void){
  if(CHECKBOX_device_orientation!=NULL)CHECKBOX_device_orientation->set_int_val(show_device_orientation);
}

/* ------------------ update_glui_devices ------------------------ */

extern "C" void update_glui_devices(void){
  Device_CB(SHOWDEVICEVALS);
  Device_CB(COLORDEVICEVALS);
  Device_CB(DEVICE_devicetypes);
}

/* ------------------ glui_device_setup ------------------------ */

extern "C" void glui_device_setup(int main_window){

  update_glui_device=0;
  if(glui_device!=NULL){
    glui_device->close();
    glui_device=NULL;
  }
  glui_device = GLUI_Master.create_glui("Device/Objects",0,0,0);
  glui_device->hide();

  if(ndeviceinfo>0){
    int i;

    PANEL_objects = glui_device->add_panel("Devices/Objects",false);

    PANEL_smvobjects = glui_device->add_panel_to_panel(PANEL_objects,"Objects",true);
    SPINNER_sensorrelsize=glui_device->add_spinner_to_panel(PANEL_smvobjects,_("Scale"),GLUI_SPINNER_FLOAT,&sensorrelsize,DEVICE_sensorsize,Device_CB);
    CHECKBOX_device_3=glui_device->add_checkbox_to_panel(PANEL_smvobjects,_("Outline"),&object_outlines);
    CHECKBOX_device_orientation=glui_device->add_checkbox_to_panel(PANEL_smvobjects,_("Orientation"),&show_device_orientation,DEVICE_show_orientation,Device_CB);
    SPINNER_orientation_scale=glui_device->add_spinner_to_panel(PANEL_smvobjects,_("Orientation scale"),GLUI_SPINNER_FLOAT,&orientation_scale);
    SPINNER_orientation_scale->set_float_limits(0.1,10.0);

    if(get_num_activedevices()>0||isZoneFireModel==1){
      PANEL_velocityvectors = glui_device->add_panel_to_panel(PANEL_objects, "Flow Vectors", true);
      if(nvdeviceinfo==0)PANEL_velocityvectors->disable();
      CHECKBOX_device_1=glui_device->add_checkbox_to_panel(PANEL_velocityvectors,_("Show"),&showvdeviceval);
      PANEL_vector_type=glui_device->add_panel_to_panel(PANEL_velocityvectors,"type",true);
      RADIO_vectortype=glui_device->add_radiogroup_to_panel(PANEL_vector_type,&vectortype);
      glui_device->add_radiobutton_to_group(RADIO_vectortype,"line");
      glui_device->add_radiobutton_to_group(RADIO_vectortype,"arrow");
      glui_device->add_radiobutton_to_group(RADIO_vectortype,"object");
      ROLLOUT_arrow_dimensions=glui_device->add_rollout_to_panel(PANEL_velocityvectors,"dimensions",false);
      PANEL_arrow_base=glui_device->add_panel_to_panel(ROLLOUT_arrow_dimensions,"base",true);
      glui_device->add_spinner_to_panel(PANEL_arrow_base,_("height"),GLUI_SPINNER_FLOAT,&vector_baseheight);
      glui_device->add_spinner_to_panel(PANEL_arrow_base,_("diameter"),GLUI_SPINNER_FLOAT,&vector_basediameter);
      PANEL_arrow_height=glui_device->add_panel_to_panel(ROLLOUT_arrow_dimensions,"height",true);
      glui_device->add_spinner_to_panel(PANEL_arrow_height,_("height"),GLUI_SPINNER_FLOAT,&vector_headheight);
      glui_device->add_spinner_to_panel(PANEL_arrow_height,_("diameter"),GLUI_SPINNER_FLOAT,&vector_headdiameter);
#ifdef pp_PILOT
      glui_device->add_checkbox_to_panel(PANEL_velocityvectors,_("Pilot view"),&vispilot);
#endif

      PANEL_devicevalues = glui_device->add_panel_to_panel(PANEL_objects,"Device values",true);

      CHECKBOX_device_2=glui_device->add_checkbox_to_panel(PANEL_devicevalues,_("Show values"),&showdeviceval,SHOWDEVICEVALS,Device_CB);
      CHECKBOX_device_4=glui_device->add_checkbox_to_panel(PANEL_devicevalues,_("Color"),&colordeviceval,COLORDEVICEVALS,Device_CB);
      glui_device->add_spinner_to_panel(PANEL_devicevalues,"min",GLUI_SPINNER_FLOAT,&device_valmin);
      glui_device->add_spinner_to_panel(PANEL_devicevalues,"max",GLUI_SPINNER_FLOAT,&device_valmax);

      PANEL_devicevis=glui_device->add_panel_to_panel(PANEL_devicevalues,"",false);
      devicetypes_index=CLAMP(devicetypes_index,0,ndevicetypes-1);
      RADIO_devicetypes=glui_device->add_radiogroup_to_panel(PANEL_devicevis,&devicetypes_index,DEVICE_devicetypes,Device_CB);
      for(i=0;i<ndevicetypes;i++){
        glui_device->add_radiobutton_to_group(RADIO_devicetypes,devicetypes[i]->quantity);
      }

      update_glui_devices();
    }
  }

  PANEL_label3 = glui_device->add_panel("",false);
  glui_device->add_column_to_panel(PANEL_label3,false);

  BUTTON_device_1=glui_device->add_button_to_panel(PANEL_label3,_("Save settings"),SAVE_SETTINGS,Device_CB);
  glui_device->add_column_to_panel(PANEL_label3,false);

  BUTTON_device_2=glui_device->add_button_to_panel(PANEL_label3,_("Close"),DEVICE_close,Device_CB);

  glui_device->set_main_gfx_window( main_window );
}

/* ------------------ hide_glui_device ------------------------ */

extern "C" void hide_glui_device(void){
  if(glui_device!=NULL)glui_device->hide();
  updatemenu=1;
}

/* ------------------ show_glui_device ------------------------ */

extern "C" void show_glui_device(void){
  if(glui_device!=NULL)glui_device->show();
}

/* ------------------ Device_CB ------------------------ */

void Open_CB(int var){
  int i;
  filelistdata *filei;
  char *open_filter_ptr;


  switch(var){
    case OPEN_UP:
      strcat(gluiopen_path_dir,dirseparator);
      strcat(gluiopen_path_dir,"..");
      Open_CB(OPEN_UPDATE_LIST);
      break;
    case OPEN_DOWN:
      if(gluiopen_filelist==NULL)break;
      filei = gluiopen_filelist + gluiopen_file_index;
      if(filei->type==1){
        strcat(gluiopen_path_dir,dirseparator);
        strcat(gluiopen_path_dir,filei->file);
        Open_CB(OPEN_UPDATE_LIST);
      }
      break;
    case OPEN_FILEINDEX:
      PRINTF("in OPEN_FILEINDEX\n");
      if(gluiopen_filelist==NULL)break;
      filei = gluiopen_filelist + gluiopen_file_index;
      if(filei->type==1){
        BUTTON_open_down->enable();
      }
      else{
        BUTTON_open_down->disable();
      }
      break;
    case OPEN_OPEN:
      if(gluiopen_filelist==NULL)break;
      filei = gluiopen_filelist + gluiopen_file_index;
      if(filei->type==1){
        strcat(gluiopen_path_dir,dirseparator);
        strcat(gluiopen_path_dir,filei->file);
        Open_CB(OPEN_UPDATE_LIST);
      }
      else{
        PRINTF("opening file: %s\n",filei->file);
      }
      break;
    case OPEN_CANCEL:
      break;
    case OPEN_FILTER:
      break;
    case OPEN_APPLY_FILTER:
      strcpy(gluiopen_filter2,gluiopen_filter);
      trim(gluiopen_filter2);
      open_filter_ptr = trim_front(gluiopen_filter2);
      EDIT_filter->set_text(open_filter_ptr);
      Open_CB(OPEN_UPDATE_LIST);
      break;
    case OPEN_UPDATE_LIST:
      LIST_open->delete_item("");
      for(i=0;i<gluiopen_nfilelist;i++){
        char label[1024];

        strcpy(label,"");
        if(gluiopen_filelist[i].type==1){
          strcat(label,"> ");
        }
        strcat(label,gluiopen_filelist[i].file);
        LIST_open->delete_item(label);
      }
      free_filelist(gluiopen_filelist,&gluiopen_nfilelist);
      gluiopen_nfilelist=get_nfilelist(gluiopen_path_dir,gluiopen_filter);
      if(gluiopen_nfilelist==0){
        LIST_open->add_item(0,"");
      }
      get_filelist(gluiopen_path_dir, gluiopen_filter,gluiopen_nfilelist,&gluiopen_filelist);
      if(gluiopen_nfilelist>0&&gluiopen_filelist[0].type==1){
        BUTTON_open_down->enable();
      }
      else{
        BUTTON_open_down->disable();
      }
      for(i=0;i<gluiopen_nfilelist;i++){
        char label[1024];

        strcpy(label,"");
        if(gluiopen_filelist[i].type==1){
          strcat(label,"> ");
        }
        strcat(label,gluiopen_filelist[i].file);
        LIST_open->add_item(i,label);
      }
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
}

/* ------------------ Device_CB ------------------------ */

void Device_CB(int var){
  int i;

  updatemenu=1;
  switch(var){
  case DEVICE_show_orientation:
    updatemenu=1;
    break;
  case DEVICE_devicetypes:
    for(i=0;i<ndevicetypes;i++){
      devicetypes[i]->type2vis=0;
    }
    if(ndevicetypes>0){
      devicetypes[devicetypes_index]->type2vis=1;
      update_colordevs();
    }
    break;
  case SHOWDEVICEVALS:
  case COLORDEVICEVALS:
    if(PANEL_devicevis!=NULL){
      if(colordeviceval==1||showdeviceval==1){
        PANEL_devicevis->enable();
      }
      else{
        PANEL_devicevis->disable();
      }
    }

    break;
  case DEVICE_sensorsize:
    if(sensorrelsize<sensorrelsizeMIN){
      sensorrelsize=sensorrelsizeMIN;
      if(SPINNER_sensorrelsize!=NULL){
        SPINNER_sensorrelsize->set_float_val(sensorrelsize);
      }
    }
    break;
  case SAVE_SETTINGS:
    writeini(LOCAL_INI,NULL);
    break;
  case DEVICE_close:
    hide_glui_device();
    break;
  default:
    ASSERT(FFALSE);
  }
}
