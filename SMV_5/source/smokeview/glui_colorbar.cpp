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
extern "C" char glui_colorbar_revision[]="$Revision$";

//
// setColorbarClipPlanes(1);

GLUI_Panel *panel_cb1=NULL;
GLUI_Panel *panel_cb2=NULL;
GLUI_Panel *panel_cb2L=NULL;
GLUI_Panel *panel_cb2R=NULL;
GLUI_Panel *panel_cb2R2=NULL;
GLUI_Panel *panel_cb3=NULL;
GLUI_Panel *panel_cb4=NULL;
GLUI_Panel *panel_cb4L=NULL;
GLUI_Panel *panel_cb4R=NULL;
GLUI_Panel *panel_point=NULL;
GLUI_Panel *panel_cb5=NULL;
GLUI_Panel *panel_cb5a=NULL;
GLUI_Panel *panel_cb5b=NULL;
GLUI_Panel *panel_cb6=NULL;
GLUI_Panel *panel_cb7=NULL;
GLUI_Panel *panel_cb8=NULL;

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
GLUI_Button *BUTTON_delete=NULL;
GLUI_Button *BUTTON_addpoint=NULL;
GLUI_Button *BUTTON_deletepoint=NULL;
GLUI_Button *BUTTON_savesettings=NULL;
GLUI_Button *BUTTON_update=NULL;
GLUI_Checkbox *CHECKBOX_jump=NULL;
GLUI_Spinner *SPINNER_valmin=NULL;
GLUI_Spinner *SPINNER_valmax=NULL;
GLUI_Spinner *SPINNER_val=NULL;
GLUI_Spinner *SPINNER_colorindex=NULL;
GLUI_Checkbox *CHECKBOX_hidesv=NULL;
GLUI_EditText *EDITTEXT_colorbar_label=NULL;
GLUI_StaticText *STATICTEXT_left=NULL, *STATICTEXT_right=NULL;


int selectedcolorbar_index;
int cb_rgb[6];
int cb_jump;

#define COLORBAR_LIST 0
#define COLORBAR_CLOSE 1
#define COLORBAR_RGB 2
#define COLORBAR_NEXT 3
#define COLORBAR_PREV 4
#define COLORBAR_NEW 5
#define COLORBAR_VAL 6
#define COLORBAR_ADDPOINT 7
#define COLORBAR_DELETEPOINT 8
#define COLORBAR_SAVE 9
#define COLORBAR_LABEL 10
#define COLORBAR_UPDATE 11
#define COLORBAR_COLORINDEX 12
#define COLORBAR_MINMAX 13
#define COLORBAR_DELETE 14

extern "C" void colorbar_global2local(void);

void COLORBAR_CB(int var);

/* ------------------ update_camera_label ------------------------ */

extern "C" void update_colorbar_label(void){
  EDITTEXT_colorbar_label->set_text(colorbar_label);
}

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

  cb_valmin=0.0;
  cb_valmax=100.0;
  cb_val=50.0;
  cb_colorindex=128;

  if(colorbar_label!=NULL){
    free(colorbar_label);
    colorbar_label=NULL;
  }
  colorbar_label=(char *)malloc(sizeof(GLUI_String));
  strcpy(colorbar_label,"New Colorbar");

  if(glui_colorbar!=NULL)glui_colorbar->close();
  glui_colorbar = GLUI_Master.create_glui("Colorbar Editor",0,0,0);
  if(showcolorbar==0)glui_colorbar->hide();

  panel_cb2R2 = glui_colorbar->add_panel("",GLUI_PANEL_NONE);
  BUTTON_new=glui_colorbar->add_button_to_panel(panel_cb2R2,"New Colorbar",COLORBAR_NEW,COLORBAR_CB);
  glui_colorbar->add_column_to_panel(panel_cb2R2,false);
  BUTTON_delete=glui_colorbar->add_button_to_panel(panel_cb2R2,"Delete Colorbar",COLORBAR_DELETE,COLORBAR_CB);
  glui_colorbar->add_column_to_panel(panel_cb2R2,false);
  cb_hidesv=1;
  CHECKBOX_hidesv = glui_colorbar->add_checkbox_to_panel(panel_cb2R2,"Hide Scene",&cb_hidesv);

  panel_cb1 = glui_colorbar->add_panel("Colorbar");
  if(ncolorbars>0){
    selectedcolorbar_index=-1;
    LISTBOX_colorbar=glui_colorbar->add_listbox_to_panel(panel_cb1,"",&selectedcolorbar_index,COLORBAR_LIST,COLORBAR_CB);

    for(i=0;i<ncolorbars;i++){
      cbi = colorbarinfo + i;
      cbi->label_ptr=cbi->label;
      LISTBOX_colorbar->add_item(i,cbi->label_ptr);
    }
    LISTBOX_colorbar->set_int_val(colorbartype);
  }
  EDITTEXT_colorbar_label=glui_colorbar->add_edittext_to_panel(panel_cb1,"Label",GLUI_EDITTEXT_TEXT,colorbar_label,COLORBAR_LABEL,COLORBAR_CB);  
  BUTTON_update=glui_colorbar->add_button_to_panel(panel_cb1,"Update label",COLORBAR_UPDATE,COLORBAR_CB);
  glui_colorbar->add_column_to_panel(panel_cb1,false);

  SPINNER_valmin=  glui_colorbar->add_spinner_to_panel(panel_cb1,"Min value",  GLUI_SPINNER_FLOAT,&cb_valmin,  COLORBAR_MINMAX,COLORBAR_CB);
  SPINNER_valmax=  glui_colorbar->add_spinner_to_panel(panel_cb1,"Max value",  GLUI_SPINNER_FLOAT,&cb_valmax,  COLORBAR_MINMAX,COLORBAR_CB);

  panel_point = glui_colorbar->add_panel("Colorbar Node");
  
  panel_cb5 = glui_colorbar->add_panel_to_panel(panel_point,"",GLUI_PANEL_NONE);

  BUTTON_prev=glui_colorbar->add_button_to_panel(panel_cb5,"Previous",COLORBAR_PREV,COLORBAR_CB);
  BUTTON_deletepoint=glui_colorbar->add_button_to_panel(panel_cb5,"Delete",COLORBAR_DELETEPOINT,COLORBAR_CB);
  
  glui_colorbar->add_column_to_panel(panel_cb5,false);

  BUTTON_next=glui_colorbar->add_button_to_panel(panel_cb5,"Next",COLORBAR_NEXT,COLORBAR_CB);
  BUTTON_addpoint=glui_colorbar->add_button_to_panel(panel_cb5,"Insert",COLORBAR_ADDPOINT,COLORBAR_CB);

  panel_cb4 = glui_colorbar->add_panel_to_panel(panel_point,"",GLUI_PANEL_NONE);
  STATICTEXT_left=glui_colorbar->add_statictext_to_panel(panel_cb4,"");
  SPINNER_right_red=  glui_colorbar->add_spinner_to_panel(panel_cb4,"red",  GLUI_SPINNER_INT,cb_rgb,COLORBAR_RGB,COLORBAR_CB);
  SPINNER_right_green=glui_colorbar->add_spinner_to_panel(panel_cb4,"green",GLUI_SPINNER_INT,cb_rgb+1,COLORBAR_RGB,COLORBAR_CB);
  SPINNER_right_blue= glui_colorbar->add_spinner_to_panel(panel_cb4,"blue", GLUI_SPINNER_INT,cb_rgb+2,COLORBAR_RGB,COLORBAR_CB);

  glui_colorbar->add_column_to_panel(panel_cb4,false);
  STATICTEXT_right=glui_colorbar->add_statictext_to_panel(panel_cb4,"Color");
  SPINNER_left_red=  glui_colorbar->add_spinner_to_panel(panel_cb4,"",  GLUI_SPINNER_INT,cb_rgb+3,  COLORBAR_RGB,COLORBAR_CB);
  SPINNER_left_green=glui_colorbar->add_spinner_to_panel(panel_cb4,"",GLUI_SPINNER_INT,cb_rgb+4,COLORBAR_RGB,COLORBAR_CB);
  SPINNER_left_blue= glui_colorbar->add_spinner_to_panel(panel_cb4,"", GLUI_SPINNER_INT,cb_rgb+5,COLORBAR_RGB,COLORBAR_CB);


  panel_cb7 = glui_colorbar->add_panel_to_panel(panel_point,"",GLUI_PANEL_NONE);
  CHECKBOX_jump = glui_colorbar->add_checkbox_to_panel(panel_cb7,"Split Colorbar",&cb_jump,COLORBAR_RGB,COLORBAR_CB);
  glui_colorbar->add_column_to_panel(panel_cb7,false);
  SPINNER_val=  glui_colorbar->add_spinner_to_panel(panel_cb7,"Value",  GLUI_SPINNER_FLOAT,&cb_val,  COLORBAR_VAL,COLORBAR_CB);
  SPINNER_colorindex=  glui_colorbar->add_spinner_to_panel(panel_cb7,"Colorbar index",  GLUI_SPINNER_INT,&cb_colorindex,  COLORBAR_COLORINDEX,COLORBAR_CB);
  SPINNER_colorindex->set_int_limits(0,255);

  SPINNER_left_red->set_int_limits(0,255);
  SPINNER_left_green->set_int_limits(0,255);
  SPINNER_left_blue->set_int_limits(0,255);
  SPINNER_right_red->set_int_limits(0,255);
  SPINNER_right_green->set_int_limits(0,255);
  SPINNER_right_blue->set_int_limits(0,255);

  colorbar_global2local();

  panel_cb8 = glui_colorbar->add_panel("",GLUI_PANEL_NONE);
  glui_colorbar->add_button_to_panel(panel_cb8,"Save Settings",COLORBAR_SAVE,COLORBAR_CB);
  glui_colorbar->add_column_to_panel(panel_cb8,false);
  glui_colorbar->add_button_to_panel(panel_cb8,"Close",COLORBAR_CLOSE,COLORBAR_CB);

  glui_colorbar->set_main_gfx_window( main_window );
}

/* ------------------ COLORBAR_CB ------------------------ */

void COLORBAR_CB(int var){
  colorbardata *cbi;
  float *rgb_nodes;
  int i;
  float *rgb1, *rgb2;
  float *rgb1a, *rgb2a;
  float *rgb0, *rgb0a;

  switch (var){
  case COLORBAR_COLORINDEX:
    if(colorbartype>=ndefaultcolorbars&&colorbartype<ncolorbars){
      cbi = colorbarinfo + colorbartype;
      
//    SPINNER_val->set_float_val(cbi->legvals[colorbarpoint]);
//    SPINNER_colorindex->set_int_val(cbi->colorbar_index[colorbarpoint]);
      cb_val=cb_valmin + (cb_valmax-cb_valmin)*((float)cb_colorindex/255.0);
      cbi->legvals[colorbarpoint]=cb_val;
      cbi->colorbar_index[colorbarpoint]=cb_colorindex;

      colorbar_global2local();
      remapcolorbar(cbi);
      updatecolors(-1);
    }
    break;
  case COLORBAR_MINMAX:
    if(colorbartype>=ndefaultcolorbars&&colorbartype<ncolorbars){
      cbi = colorbarinfo + colorbartype;
      for(i=0;i<cbi->nlegs-1;i++){
        cbi->legvals[i]=cb_valmin + (float)i*(cb_valmax-cb_valmin)/(cbi->nlegs-1);
        cbi->colorbar_index[i]=255*(cbi->legvals[i]-cb_valmin)/(cb_valmax-cb_valmin);
      }
      cbi->legvals[cbi->nlegs-1]=cb_valmax;
      cbi->colorbar_index[cbi->nlegs-1]=255;
      cbi->valmin=cb_valmin;
      cbi->valmax=cb_valmax;

      SPINNER_val->set_float_limits(cb_valmin,cb_valmax);
      colorbar_global2local();
      remapcolorbar(cbi);
      updatecolors(-1);
    }
    break;
  case COLORBAR_VAL:
    if(colorbartype>=ndefaultcolorbars&&colorbartype<ncolorbars){
      cbi = colorbarinfo + colorbartype;
      cb_colorindex=255*(cb_val-cb_valmin)/(cb_valmax-cb_valmin);
      if(cb_colorindex<0)cb_colorindex=0;
      if(cb_colorindex>255)cb_colorindex=255;
      cbi->legvals[colorbarpoint]=cb_val;
      cbi->colorbar_index[colorbarpoint]=cb_colorindex;
      colorbar_global2local();
      remapcolorbar(cbi);
      updatecolors(-1);
    }
    break;
  case COLORBAR_UPDATE:
    COLORBAR_CB(COLORBAR_LABEL);
    break;
  case COLORBAR_LABEL:
    if(colorbartype>=ndefaultcolorbars&&colorbartype<ncolorbars){
      cbi = colorbarinfo + colorbartype;
      strcpy(cbi->label,colorbar_label);
      LISTBOX_colorbar->delete_item(colorbartype);
      LISTBOX_colorbar->add_item(colorbartype,colorbar_label);
      LISTBOX_colorbar->set_int_val(0);
      LISTBOX_colorbar->set_int_val(colorbartype);
      updatemenu=1;
    }
    break;
  case COLORBAR_SAVE:
    updatemenu=1;
    writeini(LOCAL_INI);
    break;
  case COLORBAR_ADDPOINT:
    if(colorbartype<ndefaultcolorbars||colorbartype>=ncolorbars)return;
    cbi = colorbarinfo + colorbartype;
    if(colorbarpoint<=0||colorbarpoint>cbi->nlegs-1)return;

    cbi->nlegs++;
    for(i=cbi->nlegs-1;i>=colorbarpoint+1;i--){
      cbi->legvals[i]=cbi->legvals[i-1];
      cbi->colorbar_index[i]=cbi->colorbar_index[i-1];
      cbi->splitflag[i]=cbi->splitflag[i-1];
      rgb2 = cbi->leg_rgb+6*i;
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
    cbi->colorbar_index[colorbarpoint]=(cbi->colorbar_index[colorbarpoint]+cbi->colorbar_index[colorbarpoint-1])/2;
    cbi->legvals[colorbarpoint]=(cbi->legvals[colorbarpoint]+cbi->legvals[colorbarpoint-1])/2.0;

    cbi->splitflag[colorbarpoint]=0;
    rgb1 = cbi->leg_rgb+6*colorbarpoint;
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

    if(colorbarpoint==cbi->nlegs)colorbarpoint=cbi->nlegs-1;
    break;
  case COLORBAR_DELETEPOINT:
    if(colorbartype<ndefaultcolorbars||colorbartype>=ncolorbars)return;
    cbi = colorbarinfo + colorbartype;
    if(colorbarpoint<0||colorbarpoint>cbi->nlegs-1)return;

    if(cbi->nlegs<=2)return;
    for(i=colorbarpoint+1;i<cbi->nlegs;i++){
      cbi->legvals[i-1]=cbi->legvals[i];
      cbi->splitflag[i-1]=cbi->splitflag[i];
      cbi->colorbar_index[i-1]=cbi->colorbar_index[i];
      rgb2 = cbi->leg_rgb+6*i;
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
    cbi->nlegs--;
    remapcolorbar(cbi);
    updatecolors(-1);
    if(colorbarpoint==cbi->nlegs)colorbarpoint=cbi->nlegs-1;
    break;
  case COLORBAR_RGB:
    if(colorbartype<0||colorbartype>=ncolorbars)return;
    cbi = colorbarinfo + colorbartype;
    if(colorbarpoint<0||colorbarpoint>cbi->nlegs-1)return;
    if(cbi->splitflag[colorbarpoint]!=cb_jump){
      cbi->splitflag[colorbarpoint]=cb_jump;
      colorbar_global2local();
    }
    rgb_nodes=cbi->leg_rgb+6*colorbarpoint;

    for(i=0;i<3;i++){
      rgb_nodes[i]=(float)cb_rgb[i+3]/255.0;
    }
    if(cbi->splitflag[colorbarpoint]==1){
      for(i=3;i<6;i++){
        rgb_nodes[i-6]=(float)cb_rgb[i-3]/255.0;
      }
    }
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
      if(colorbarpoint>cbi->nlegs-1)colorbarpoint=0;
    }
    else if(var==COLORBAR_PREV){
      colorbarpoint--;
      if(colorbarpoint<0)colorbarpoint=cbi->nlegs-1;
    }
    cbi->legindex=colorbarpoint;
    colorbar_global2local();
    break;
  case COLORBAR_NEW:
    if(colorbartype<0||colorbartype>=ncolorbars)return;
    addcolorbar(colorbartype);
    colorbartype=ncolorbars-1;
    cbi = colorbarinfo + colorbartype;  //addcolorbar resizes (and possibly moves) colorbarinfo
    LISTBOX_colorbar->add_item(colorbartype,cbi->label);
    LISTBOX_colorbar->set_int_val(colorbartype);
    COLORBAR_CB(COLORBAR_LIST);
    break;
  case COLORBAR_DELETE:
    if(colorbartype>=ndefaultcolorbars&&colorbartype<ncolorbars){
      colorbardata *cb_from, *cb_to;

      for(i=colorbartype;i<ncolorbars-1;i++){
        cb_to = colorbarinfo + i;
        cb_from = cb_to + 1;
        memcpy(cb_to,cb_from,sizeof(colorbardata));
        cb_to->label_ptr=cb_to->label;
      }
      for(i=colorbartype;i<ncolorbars;i++){
        LISTBOX_colorbar->delete_item(i);
      }
      ncolorbars--;
      for(i=colorbartype;i<ncolorbars;i++){
        cbi = colorbarinfo + i;
        LISTBOX_colorbar->add_item(i,cbi->label_ptr);
      }
      if(colorbartype==ncolorbars)colorbartype--;
      LISTBOX_colorbar->set_int_val(0);
      COLORBAR_CB(COLORBAR_LIST);
    }
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
    colorbarpoint=cbi->legindex;
    
    cb_jump = cbi->splitflag[colorbarpoint];
    CHECKBOX_jump->set_int_val(cb_jump);
    
    SPINNER_val->set_float_val(cbi->legvals[colorbarpoint]);
    SPINNER_colorindex->set_int_val(cbi->colorbar_index[colorbarpoint]);
    SPINNER_valmin->set_float_val(  cbi->valmin);
    SPINNER_valmax->set_float_val(  cbi->valmax);

    ii = 6*colorbarpoint;
    rgb = cbi->leg_rgb+ii;

    BUTTON_next->enable();
    BUTTON_prev->enable();

    strcpy(colorbar_label,cbi->label);
    EDITTEXT_colorbar_label->set_text(colorbar_label);
    SPINNER_left_red->set_int_val(   (int)(255*rgb[0]));
    SPINNER_left_green->set_int_val( (int)(255*rgb[1]));
    SPINNER_left_blue->set_int_val(  (int)(255*rgb[2]));
    icolorbar=LISTBOX_colorbar->get_int_val();
    if(icolorbar>=ndefaultcolorbars){
      BUTTON_delete->enable();
      SPINNER_valmin->enable();
      SPINNER_valmax->enable();
      EDITTEXT_colorbar_label->enable();
      SPINNER_left_red->enable();
      SPINNER_left_green->enable();
      SPINNER_left_blue->enable();
      if(cbi->splitflag[colorbarpoint]==1){
        SPINNER_right_red->enable();
        SPINNER_right_green->enable();
        SPINNER_right_blue->enable();
      }
      else{
        SPINNER_right_red->disable();
        SPINNER_right_green->disable();
        SPINNER_right_blue->disable();
      }
      BUTTON_addpoint->enable();
      BUTTON_deletepoint->enable();
      if(colorbarpoint==0||colorbarpoint==cbi->nlegs-1){
        CHECKBOX_jump->disable();
        SPINNER_val->disable();
        SPINNER_colorindex->disable();
      }
      else{
        CHECKBOX_jump->enable();
        SPINNER_val->enable();
        SPINNER_colorindex->enable();
      }
    }
    else{
      BUTTON_delete->disable();
      SPINNER_valmin->disable();
      SPINNER_valmax->disable();
      EDITTEXT_colorbar_label->disable();
      BUTTON_addpoint->disable();
      BUTTON_deletepoint->disable();
      CHECKBOX_jump->disable();
      SPINNER_val->disable();
      SPINNER_colorindex->disable();

      SPINNER_left_red->disable();
      SPINNER_left_green->disable();
      SPINNER_left_blue->disable();
      SPINNER_right_red->disable();
      SPINNER_right_green->disable();
      SPINNER_right_blue->disable();
    }
    if(colorbarpoint>0){
      SPINNER_right_red->set_float_val(  (int)(255*rgb[3-6]));
      SPINNER_right_green->set_float_val((int)(255*rgb[4-6]));
      SPINNER_right_blue->set_float_val( (int)(255*rgb[5-6]));
    }
  }
  if(CHECKBOX_jump->get_int_val()==1){
    STATICTEXT_left->set_text("Color left of split");
    STATICTEXT_right->set_text("Color right of split");
  }
  else{
    STATICTEXT_left->set_text("");
    STATICTEXT_right->set_text("Color");
  }
}
