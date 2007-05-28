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

GLUI_Listbox *LISTBOX_colorbar=NULL;
GLUI *glui_colorbar=NULL;
GLUI_Spinner *SPINNER_left_red=NULL;
GLUI_Spinner *SPINNER_left_green=NULL;
GLUI_Spinner *SPINNER_left_blue=NULL;
GLUI_Spinner *SPINNER_right_red=NULL;
GLUI_Spinner *SPINNER_right_green=NULL;
GLUI_Spinner *SPINNER_right_blue=NULL;
GLUI_Button *BUTTON_next=NULL,*BUTTON_prev=NULL;

int selectedcolorbar_index;
float cb_rgb[6];

#define COLORBAR_LIST 0
#define COLORBAR_CLOSE 1
#define COLORBAR_RGB 2
#define COLORBAR_NEXT 3
#define COLORBAR_PREV 4

extern "C" void colorbar_global2local(void);

void COLORBAR_CB(int var);

/* ------------------ hide_glui_colorbar ------------------------ */

void hide_glui_colorbar(void){
  if(glui_colorbar!=NULL)glui_colorbar->hide();
  showcolorbar=0;
  updatemenu=1;
}

/* ------------------ show_glui_colorbar ------------------------ */

void show_glui_colorbar(void){
  if(glui_colorbar!=NULL)glui_colorbar->show();
}

/* ------------------ glui_colorbar_setup ------------------------ */

extern "C" void glui_colorbar_setup(int main_window){
  colorbardata *cbi;
  int i;

  if(glui_colorbar!=NULL)glui_colorbar->close();
  glui_colorbar = GLUI_Master.create_glui("Colorbar Editor",0,0,0);
  if(showcolorbar==0)glui_colorbar->hide();

  if(ncolorbars>0){
    selectedcolorbar_index=-1;
    LISTBOX_colorbar=glui_colorbar->add_listbox("Select Colorbar:",&selectedcolorbar_index,COLORBAR_LIST,COLORBAR_CB);

    for(i=0;i<ncolorbars;i++){
      cbi = colorbarinfo + i;
      LISTBOX_colorbar->add_item(i,cbi->label);
    }
    LISTBOX_colorbar->set_int_val(colorbartype);
  }
  SPINNER_left_red=  glui_colorbar->add_spinner("left red",  GLUI_SPINNER_FLOAT,cb_rgb,  COLORBAR_RGB,COLORBAR_CB);
  SPINNER_left_green=glui_colorbar->add_spinner("left green",GLUI_SPINNER_FLOAT,cb_rgb+1,COLORBAR_RGB,COLORBAR_CB);
  SPINNER_left_blue= glui_colorbar->add_spinner("left blue", GLUI_SPINNER_FLOAT,cb_rgb+2,COLORBAR_RGB,COLORBAR_CB);
  SPINNER_right_red=  glui_colorbar->add_spinner("right red",  GLUI_SPINNER_FLOAT,cb_rgb+3,COLORBAR_RGB,COLORBAR_CB);
  SPINNER_right_green=glui_colorbar->add_spinner("right green",GLUI_SPINNER_FLOAT,cb_rgb+4,COLORBAR_RGB,COLORBAR_CB);
  SPINNER_right_blue= glui_colorbar->add_spinner("right blue", GLUI_SPINNER_FLOAT,cb_rgb+5,COLORBAR_RGB,COLORBAR_CB);
  SPINNER_left_red->set_float_limits(0.0,1.0);
  SPINNER_left_green->set_float_limits(0.0,1.0);
  SPINNER_left_blue->set_float_limits(0.0,1.0);
  SPINNER_right_red->set_float_limits(0.0,1.0);
  SPINNER_right_green->set_float_limits(0.0,1.0);
  SPINNER_right_blue->set_float_limits(0.0,1.0);
  colorbar_global2local();

  BUTTON_next=glui_colorbar->add_button("Next",COLORBAR_NEXT,COLORBAR_CB);
  BUTTON_prev=glui_colorbar->add_button("Prev",COLORBAR_PREV,COLORBAR_CB);
  glui_colorbar->add_button("Close",COLORBAR_CLOSE,COLORBAR_CB);

  glui_colorbar->set_main_gfx_window( main_window );
}

/* ------------------ COLORBAR_CB ------------------------ */

void COLORBAR_CB(int var){
  colorbardata *cbi;
  float *rgb_nodes;
  int i;

  switch (var){
  case COLORBAR_RGB:
    if(colorbartype<0||colorbartype>=ncolorbars)return;
    cbi = colorbarinfo + colorbartype;
    if(colorbarframe<0||colorbarframe>=cbi->nlegs)return;
    rgb_nodes=cbi->rgbnodes+6*colorbarframe;
    for(i=0;i<6;i++){
      rgb_nodes[i]=cb_rgb[i];
    }
    initcolorbars();
    break;
  case COLORBAR_LIST:
    {
      int icolorbar;


      icolorbar=LISTBOX_colorbar->get_int_val();
      ColorBarMenu(icolorbar);
      colorbar_global2local();
    }
    break;
  case COLORBAR_CLOSE:
    glui_colorbar->hide();
    break;
  case COLORBAR_NEXT:
  case COLORBAR_PREV:
    if(colorbartype<0||colorbartype>=ncolorbars)return;
    cbi = colorbarinfo + colorbartype;
    if(var==COLORBAR_NEXT){
      colorbarframe++;
      if(colorbarframe>cbi->nlegs-1)colorbarframe=0;
    }
    else if(var==COLORBAR_PREV){
      colorbarframe--;
      if(colorbarframe<0)colorbarframe=cbi->nlegs-1;
    }
    cbi->legindex=colorbarframe;
    colorbar_global2local();
    break;
  }
}

/* ------------------ colorbar_global2local ------------------------ */

extern "C" void colorbar_global2local(void){
  colorbardata *cbi;
  int ii;
  float *rgb;
  int icolorbar;

  if(colorbartype>=0&&colorbartype<ncolorbars){
    cbi = colorbarinfo + colorbartype;
    if(cbi->rgbnodes==NULL)return;
    colorbarframe=cbi->legindex;
    ii = 6*colorbarframe;
    rgb = cbi->rgbnodes+ii;

    SPINNER_left_red->set_float_val(   rgb[0]);
    SPINNER_left_green->set_float_val( rgb[1]);
    SPINNER_left_blue->set_float_val(  rgb[2]);
    icolorbar=LISTBOX_colorbar->get_int_val();
    if(icolorbar>=ndefaultcolorbars){
      SPINNER_left_red->enable();
      SPINNER_left_green->enable();
      SPINNER_left_blue->enable();
      if(cbi->contflag[colorbarframe]==0||colorbarframe==cbi->nlegs-1){
        SPINNER_right_red->enable();
        SPINNER_right_green->enable();
        SPINNER_right_blue->enable();
      }
      else{
        SPINNER_right_red->disable();
        SPINNER_right_green->disable();
        SPINNER_right_blue->disable();
      }
    }
    else{
      SPINNER_left_red->disable();
      SPINNER_left_green->disable();
      SPINNER_left_blue->disable();
      SPINNER_right_red->disable();
      SPINNER_right_green->disable();
      SPINNER_right_blue->disable();
    }
    SPINNER_right_red->set_float_val(  rgb[3]);
    SPINNER_right_green->set_float_val(rgb[4]);
    SPINNER_right_blue->set_float_val( rgb[5]);
  }
}
