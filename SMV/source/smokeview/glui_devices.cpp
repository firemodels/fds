// $Date$ 
// $Revision$
// $Author$

// svn revision character string2
extern "C" char glui_device_revision[]="$Revision$";

#define CPP
#include "options.h"
#include <string.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <math.h>

#include "smokeviewvars.h"

#define DEVICE_sensorsize 20
#define SHOWDEVICEVALS 26
#define DEVICE_devicetypes 27
#define SAVE_SETTINGS 99
#define DEVICE_close 3

#ifdef pp_OPEN
#define OPEN_UP 0
#define OPEN_FILEINDEX 1
#define OPEN_OPEN 2
#define OPEN_CANCEL 3
#define OPEN_FILTER 4
#define OPEN_APPLY_FILTER 5
#define OPEN_UPDATE_LIST 6

void Open_CB(int var);

int gluiopen_file_index=0;
int gluiopen_nfilelist=0;
char gluiopen_path_dir[1024];
filelistdata *gluiopen_filelist;
char gluiopen_filter[sizeof(GLUI_String)];
char gluiopen_filter2[sizeof(GLUI_String)];

GLUI_Panel *gluiopen_panel_open=NULL;
GLUI_Panel *gluiopen_panel_open2=NULL;
GLUI_Panel *gluiopen_panel_open3=NULL;
GLUI_Listbox *gluiopen_LISTBOX_open=NULL;
GLUI_EditText *gluiopen_EDIT_filter=NULL;
#endif


GLUI *glui_device=NULL;
GLUI_Panel *panel_objects=NULL;
GLUI_Spinner *SPINNER_sensorrelsize=NULL;
GLUI_Panel *panel_devicevis=NULL;
GLUI_RadioGroup *RADIO_devicetypes=NULL;
GLUI_Panel *panel_label3=NULL;

void Device_CB(int var);


/* ------------------ glui_device_setup ------------------------ */

extern "C" void glui_device_setup(int main_window){

  if(glui_device!=NULL)glui_device->close();
  glui_device = GLUI_Master.create_glui("Device/Ojbects",0,0,0);
  if(showdevice_dialog==0)glui_device->hide();


  if(ndeviceinfo>0){
    int i;

    panel_objects = glui_device->add_panel("Devices/Objects",false);
    SPINNER_sensorrelsize=glui_device->add_spinner_to_panel(panel_objects,_("Scaling"),GLUI_SPINNER_FLOAT,&sensorrelsize,DEVICE_sensorsize,Device_CB);
    if(ndevicetypes>0){
      glui_device->add_checkbox_to_panel(panel_objects,_("Show velocity vectors"),&showvdeviceval);
      glui_device->add_checkbox_to_panel(panel_objects,_("Show values"),&showdeviceval,SHOWDEVICEVALS,Device_CB);
      glui_device->add_checkbox_to_panel(panel_objects,_("Outline"),&object_outlines);
      panel_devicevis=glui_device->add_panel_to_panel(panel_objects,"",false);
      RADIO_devicetypes=glui_device->add_radiogroup_to_panel(panel_devicevis,&devicetypes_index,DEVICE_devicetypes,Device_CB);
      for(i=0;i<ndevicetypes;i++){
        glui_device->add_radiobutton_to_group(RADIO_devicetypes,devicetypes[i]->quantity);
      }
      Device_CB(SHOWDEVICEVALS);
      Device_CB(DEVICE_devicetypes);
    }

  }

  panel_label3 = glui_device->add_panel("",false);
  glui_device->add_column_to_panel(panel_label3,false);

  glui_device->add_button_to_panel(panel_label3,_("Save settings"),SAVE_SETTINGS,Device_CB);
  glui_device->add_column_to_panel(panel_label3,false);

  glui_device->add_button_to_panel(panel_label3,_("Close"),DEVICE_close,Device_CB);

#ifdef pp_OPEN
  strcpy(gluiopen_filter,"*.csv");
  gluiopen_panel_open = glui_device->add_panel("Open",true);
  glui_device->add_button_to_panel(gluiopen_panel_open,_("Up"),OPEN_UP,Open_CB);
  gluiopen_file_index=0;
  gluiopen_LISTBOX_open=glui_device->add_listbox_to_panel(gluiopen_panel_open,"",&gluiopen_file_index,OPEN_FILEINDEX,Open_CB);
  strcpy(gluiopen_path_dir,".");
  Open_CB(OPEN_UPDATE_LIST);
  gluiopen_panel_open2 = glui_device->add_panel_to_panel(gluiopen_panel_open,"",false);
  gluiopen_EDIT_filter=glui_device->add_edittext_to_panel(gluiopen_panel_open2,"filter:",GLUI_EDITTEXT_TEXT,gluiopen_filter,OPEN_FILTER,Open_CB);
  glui_device->add_column_to_panel(gluiopen_panel_open2);
  glui_device->add_button_to_panel(gluiopen_panel_open2,_("Apply Filter"),OPEN_APPLY_FILTER,Open_CB);

  gluiopen_panel_open3 = glui_device->add_panel_to_panel(gluiopen_panel_open,"",false);
  glui_device->add_button_to_panel(gluiopen_panel_open3,_("Open"),OPEN_OPEN,Open_CB);
  glui_device->add_column_to_panel(gluiopen_panel_open3);
  glui_device->add_button_to_panel(gluiopen_panel_open3,_("Cancel"),OPEN_CANCEL,Open_CB);

#endif

  glui_device->set_main_gfx_window( main_window );
}

/* ------------------ hide_glui_device ------------------------ */

extern "C" void hide_glui_device(void){
  if(glui_device!=NULL)glui_device->hide();
  showdevice_dialog=0;
  updatemenu=1;
}

/* ------------------ show_glui_device ------------------------ */

extern "C" void show_glui_device(void){
  if(glui_device!=NULL)glui_device->show();
}

#ifdef pp_OPEN

/* ------------------ Device_CB ------------------------ */

void Open_CB(int var){
  int i;
  filelistdata *filei;
  char *open_filter_ptr;


  switch (var){
    case OPEN_UP:
      strcat(gluiopen_path_dir,dirseparator);
      strcat(gluiopen_path_dir,"..");
      Open_CB(OPEN_UPDATE_LIST);
      break;
    case OPEN_FILEINDEX:
      if(gluiopen_filelist==NULL)break;
      filei = gluiopen_filelist + gluiopen_file_index;
      if(filei->type==1){
        strcat(gluiopen_path_dir,dirseparator);
        strcat(gluiopen_path_dir,filei->file);
        Open_CB(OPEN_UPDATE_LIST);
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
        printf("opening file: %s\n",filei->file);
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
      gluiopen_EDIT_filter->set_text(open_filter_ptr);
      Open_CB(OPEN_UPDATE_LIST);
      break;
    case OPEN_UPDATE_LIST:
      gluiopen_LISTBOX_open->delete_item("");
      for(i=0;i<gluiopen_nfilelist;i++){
        char label[1024];

        strcpy(label,"");
        if(gluiopen_filelist[i].type==1){
          strcat(label,"> ");
        }
        strcat(label,gluiopen_filelist[i].file);
        gluiopen_LISTBOX_open->delete_item(label);
      }
      free_filelist(gluiopen_filelist,&gluiopen_nfilelist);
      gluiopen_nfilelist=get_nfilelist(gluiopen_path_dir,gluiopen_filter);
      if(gluiopen_nfilelist==0){
        gluiopen_LISTBOX_open->add_item(0,"");
      }
      get_filelist(gluiopen_path_dir, gluiopen_filter,gluiopen_nfilelist,&gluiopen_filelist);
      for(i=0;i<gluiopen_nfilelist;i++){
        char label[1024];

        strcpy(label,"");
        if(gluiopen_filelist[i].type==1){
          strcat(label,"> ");
        }
        strcat(label,gluiopen_filelist[i].file);
        gluiopen_LISTBOX_open->add_item(i,label);
      }
      break;
    default:
      ASSERT(0);
      break;
  }
}
#endif

/* ------------------ Device_CB ------------------------ */

void Device_CB(int var){
  int i;

  updatemenu=1;
  switch (var){
  case DEVICE_devicetypes:
    for(i=0;i<ndevicetypes;i++){
      devicetypes[i]->type2vis=0;
    }
    devicetypes[devicetypes_index]->type2vis=1;
    break;
  case SHOWDEVICEVALS:
    if(panel_devicevis!=NULL){
      if(showdeviceval==1){
        panel_devicevis->enable();
      }
      else{
        panel_devicevis->disable();
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
    writeini(LOCAL_INI);
    break;
  case DEVICE_close:
    hide_glui_device();
    break;
  default:
    ASSERT(FFALSE);
  }
}
