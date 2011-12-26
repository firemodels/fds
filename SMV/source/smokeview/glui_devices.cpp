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
#include "smokeviewvars.h"
#include "translate.h"

// svn revision character string2
extern "C" char glui_device_revision[]="$Revision$";

#define DEVICE_sensorsize 20
#define SHOWDEVICEVALS 26
#define DEVICE_devicetypes 27
#define SAVE_SETTINGS 99
#define DEVICE_close 3


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
