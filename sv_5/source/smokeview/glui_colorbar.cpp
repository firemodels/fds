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
GLUI_Button *BUTTON_new=NULL;
GLUI_Button *BUTTON_addpoint=NULL;
GLUI_Button *BUTTON_deletepoint=NULL;
GLUI_Button *BUTTON_savesettings=NULL;
GLUI_Checkbox *CHECKBOX_jump=NULL;
GLUI_Spinner *SPINNER_valmin=NULL;
GLUI_Spinner *SPINNER_valmax=NULL;
GLUI_Spinner *SPINNER_val=NULL;

int selectedcolorbar_index;
float cb_rgb[6];
int cb_jump;
float cb_valmin, cb_valmax, cb_val;

#define COLORBAR_LIST 0
#define COLORBAR_CLOSE 1
#define COLORBAR_RGB 2
#define COLORBAR_NEXT 3
#define COLORBAR_PREV 4
#define COLORBAR_NEW 5
#define COLORBAR_VALMINMAX 6
#define COLORBAR_ADDPOINT 7
#define COLORBAR_DELETEPOINT 8
#define COLORBAR_SAVE 9

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
  CHECKBOX_jump = glui_colorbar->add_checkbox("jump",&cb_jump,COLORBAR_RGB,COLORBAR_CB);
  SPINNER_valmin=  glui_colorbar->add_spinner("min",  GLUI_SPINNER_FLOAT,&cb_valmin,  COLORBAR_VALMINMAX,COLORBAR_CB);
  SPINNER_valmax=  glui_colorbar->add_spinner("max",  GLUI_SPINNER_FLOAT,&cb_valmax,  COLORBAR_VALMINMAX,COLORBAR_CB);
  SPINNER_val=  glui_colorbar->add_spinner("val",  GLUI_SPINNER_FLOAT,&cb_val,  COLORBAR_VALMINMAX,COLORBAR_CB);
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

  BUTTON_addpoint=glui_colorbar->add_button("Add point",COLORBAR_ADDPOINT,COLORBAR_CB);
  BUTTON_deletepoint=glui_colorbar->add_button("Delete point",COLORBAR_DELETEPOINT,COLORBAR_CB);

  BUTTON_next=glui_colorbar->add_button("Next",COLORBAR_NEXT,COLORBAR_CB);
  BUTTON_prev=glui_colorbar->add_button("Prev",COLORBAR_PREV,COLORBAR_CB);
  BUTTON_new=glui_colorbar->add_button("New Colorbar",COLORBAR_NEW,COLORBAR_CB);

  colorbar_global2local();

  glui_colorbar->add_button("Save Settings",COLORBAR_SAVE,COLORBAR_CB);
  glui_colorbar->add_button("Close",COLORBAR_CLOSE,COLORBAR_CB);

  glui_colorbar->set_main_gfx_window( main_window );
}

/* ------------------ COLORBAR_CB ------------------------ */

void COLORBAR_CB(int var){
  colorbardata *cbi;
  float *rgb_nodes;
  int i;
  int npoints;
  float fleg;
  float *rgb1, *rgb1a;
  float *rgb2, *rgb2a;
  float *rgb0, *rgb0a;

  switch (var){
  case COLORBAR_SAVE:
    updatemenu=1;
    writeini(LOCAL_INI);
    break;
  case COLORBAR_ADDPOINT:
    if(colorbartype<ndefaultcolorbars||colorbartype>=ncolorbars)return;
    cbi = colorbarinfo + colorbartype;
    if(colorbarpoint<=0||colorbarpoint>cbi->npoints-1)return;

    npoints = cbi->npoints+1;
    ResizeColorbar(cbi,npoints);
    for(i=cbi->npoints-1;i>=colorbarpoint+1;i--){

      cbi->c_vals[i]=cbi->c_vals[i-1];
      cbi->jumpflag[i]=0;
      rgb2 = cbi->rgbnodes+6*i;
      rgb1 = rgb2-6;
      rgb2a = rgb2+3;
      rgb1a = rgb1+3;
      rgb2[0] =rgb1[0];
      rgb2[1] =rgb1[1];
      rgb2[2] =rgb1[2];
      rgb2a[0]=rgb1a[0];
      rgb2a[1]=rgb1a[1];
      rgb2a[2]=rgb1a[2];
    }
    for(i=0;i<cbi->npoints;i++){
      cbi->flegs[i]=1.0/(float)(cbi->npoints-1);
    }

    cbi->c_vals[colorbarpoint]=(cbi->c_vals[colorbarpoint]+cbi->c_vals[colorbarpoint-1])/2.0;
    cbi->jumpflag[colorbarpoint]=0;
    rgb1 = cbi->rgbnodes+6*colorbarpoint;
    rgb2 = rgb1 + 6;
    rgb1a = rgb1 + 3;
    rgb2a = rgb2 + 3;
    rgb0 = rgb1 - 6;
    rgb0a = rgb0 + 3;

    rgb1[0] = (rgb1[0] + rgb0[0])/2.0;
    rgb1[1] = (rgb1[1] + rgb0[1])/2.0;
    rgb1[2] = (rgb1[2] + rgb0[2])/2.0;

    rgb0a[0] = rgb1[0];
    rgb0a[1] = rgb1[1];
    rgb0a[2] = rgb1[2];

    rgb1a[0] = rgb2[0];
    rgb1a[1] = rgb2[1];
    rgb1a[2] = rgb2[2];

    colorbar_global2local();
    remapcolorbar(cbi);
    updatecolors(-1);

    if(colorbarpoint==cbi->npoints)colorbarpoint=cbi->npoints-1;
    break;
  case COLORBAR_DELETEPOINT:
    if(colorbartype<ndefaultcolorbars||colorbartype>=ncolorbars)return;
    cbi = colorbarinfo + colorbartype;
    if(colorbarpoint<0||colorbarpoint>cbi->npoints-1)return;

    if(cbi->npoints<=2)return;
    for(i=colorbarpoint+1;i<cbi->npoints;i++){
      float *rgb1, *rgb1a;
      float *rgb2, *rgb2a;

      cbi->c_vals[i-1]=cbi->c_vals[i];
      cbi->jumpflag[i-1]=cbi->jumpflag[i];
      rgb2 = cbi->rgbnodes+6*i;
      rgb1 = rgb2-6;
      rgb2a = rgb2+3;
      rgb1a = rgb1+3;
      rgb1[0]=rgb2[0];
      rgb1[1]=rgb2[1];
      rgb1[2]=rgb2[2];
      rgb1a[0]=rgb2a[0];
      rgb1a[1]=rgb2a[1];
      rgb1a[2]=rgb2a[2];
    }
    cbi->npoints--;
    for(i=0;i<cbi->npoints;i++){
      cbi->flegs[i]=1.0/(float)(cbi->npoints-1);
    }
    remapcolorbar(cbi);
    updatecolors(-1);
    if(colorbarpoint==cbi->npoints)colorbarpoint=cbi->npoints-1;
    break;
  case COLORBAR_RGB:
    if(colorbartype<0||colorbartype>=ncolorbars)return;
    cbi = colorbarinfo + colorbartype;
    if(colorbarpoint<0||colorbarpoint>cbi->npoints-1)return;
    if(cbi->jumpflag[colorbarpoint]!=cb_jump){
      cbi->jumpflag[colorbarpoint]=cb_jump;
      colorbar_global2local();
    }
    rgb_nodes=cbi->rgbnodes+6*colorbarpoint;
    for(i=0;i<3;i++){
      rgb_nodes[i]=cb_rgb[i];
    }
    if(cbi->jumpflag[colorbarpoint]==1){
      for(i=3;i<6;i++){
        rgb_nodes[i-6]=cb_rgb[i];
      }
    }
   // initcolorbars();
    remapcolorbar(cbi);
    updatecolors(-1);
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
    viscolorbarpath=0;
    glui_colorbar->hide();
    break;
  case COLORBAR_NEXT:
  case COLORBAR_PREV:
    if(colorbartype<0||colorbartype>=ncolorbars)return;
    cbi = colorbarinfo + colorbartype;
    if(var==COLORBAR_NEXT){
      colorbarpoint++;
      if(colorbarpoint>cbi->npoints-1)colorbarpoint=0;
    }
    else if(var==COLORBAR_PREV){
      colorbarpoint--;
      if(colorbarpoint<0)colorbarpoint=cbi->npoints-1;
    }
    cbi->pointindex=colorbarpoint;
    colorbar_global2local();
    break;
  case COLORBAR_NEW:
    addcolorbar("test");
    LISTBOX_colorbar->add_item(ncolorbars-1,"test");
    LISTBOX_colorbar->set_int_val(ncolorbars-1);
    colorbartype=ncolorbars-1;
    COLORBAR_CB(COLORBAR_LIST);
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
    colorbarpoint=cbi->pointindex;
    cb_jump = cbi->jumpflag[colorbarpoint];
    CHECKBOX_jump->set_int_val(cb_jump);

    ii = 6*colorbarpoint;
    rgb = cbi->rgbnodes+ii;

    SPINNER_left_red->set_float_val(   rgb[0]);
    SPINNER_left_green->set_float_val( rgb[1]);
    SPINNER_left_blue->set_float_val(  rgb[2]);
    icolorbar=LISTBOX_colorbar->get_int_val();
    if(icolorbar>=ndefaultcolorbars){
      SPINNER_left_red->enable();
      SPINNER_left_green->enable();
      SPINNER_left_blue->enable();
      if(cbi->jumpflag[colorbarpoint]==1){
        SPINNER_right_red->enable();
        SPINNER_right_green->enable();
        SPINNER_right_blue->enable();
      }
      else{
        SPINNER_right_red->disable();
        SPINNER_right_green->disable();
        SPINNER_right_blue->disable();
      }
     BUTTON_next->enable();
     BUTTON_prev->enable();
     BUTTON_addpoint->enable();
     BUTTON_deletepoint->enable();
     if(colorbarpoint==0){
       CHECKBOX_jump->disable();
     }
     else{
       CHECKBOX_jump->enable();
     }
     SPINNER_val->enable();
    }
    else{
     BUTTON_next->disable();
     BUTTON_prev->disable();
     BUTTON_addpoint->disable();
     BUTTON_deletepoint->disable();
     CHECKBOX_jump->disable();
     SPINNER_val->disable();

      SPINNER_left_red->disable();
      SPINNER_left_green->disable();
      SPINNER_left_blue->disable();
      SPINNER_right_red->disable();
      SPINNER_right_green->disable();
      SPINNER_right_blue->disable();
    }
    if(colorbarpoint>0){
      SPINNER_right_red->set_float_val(  rgb[3-6]);
      SPINNER_right_green->set_float_val(rgb[4-6]);
      SPINNER_right_blue->set_float_val( rgb[5-6]);
    }
  }
}
