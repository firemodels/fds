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
GLUI_Panel *panel_extreme=NULL;
GLUI_Panel *panel_cb5=NULL;
GLUI_Panel *panel_cb5a=NULL;
GLUI_Panel *panel_cb5b=NULL;
GLUI_Panel *panel_cb6=NULL;
GLUI_Panel *panel_cb7=NULL;
GLUI_Panel *panel_cb8=NULL;
GLUI_Panel *panel_cb9=NULL;

GLUI_Listbox *LISTBOX_colorbar=NULL;
GLUI *glui_colorbar=NULL;
GLUI_Spinner *SPINNER_down_red=NULL;
GLUI_Spinner *SPINNER_down_green=NULL;
GLUI_Spinner *SPINNER_down_blue=NULL;
GLUI_Spinner *SPINNER_up_red=NULL;
GLUI_Spinner *SPINNER_up_green=NULL;
GLUI_Spinner *SPINNER_up_blue=NULL;
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
GLUI_Checkbox *CHECKBOX_usebounds=NULL;
GLUI_Spinner *SPINNER_valmin=NULL;
GLUI_Spinner *SPINNER_valmax=NULL;
GLUI_Spinner *SPINNER_val=NULL;
GLUI_Spinner *SPINNER_colorindex=NULL;
GLUI_Checkbox *CHECKBOX_hidesv=NULL;
GLUI_EditText *EDITTEXT_colorbar_label=NULL;
GLUI_StaticText *STATICTEXT_left=NULL, *STATICTEXT_right=NULL;


int selectedcolorbar_index;
int cb_rgb[3],cb_up_rgb[3],cb_down_rgb[3];
int cb_usecolorbar_extreme;

#define COLORBAR_LIST 0
#define COLORBAR_CLOSE 1
#define COLORBAR_RGB 2
#define COLORBAR_NEXT 3
#define COLORBAR_PREV 4
#define COLORBAR_NEW 5
#define COLORBAR_ADDPOINT 7
#define COLORBAR_DELETEPOINT 8
#define COLORBAR_SAVE 9
#define COLORBAR_LABEL 10
#define COLORBAR_UPDATE 11
#define COLORBAR_COLORINDEX 12
#define COLORBAR_DELETE 14
#define COLORBAR_EXTREME_RGB 15
#define COLORBAR_EXTREME 16

extern "C" void colorbar_global2local(void);

void COLORBAR_CB(int var);


/* ------------------ update_camera_label ------------------------ */

extern "C" void update_extreme(void){
  CHECKBOX_usebounds->set_int_val(show_extremedata);
}


/* ------------------ update_camera_label ------------------------ */

extern "C" void update_colorbar_type(void){
  LISTBOX_colorbar->set_int_val(colorbartype);
}

/* ------------------ update_camera_label ------------------------ */

extern "C" void update_colorbar_label(void){
  EDITTEXT_colorbar_label->set_text(colorbar_label);
}

/* ------------------ hide_glui_colorbar ------------------------ */

void hide_glui_colorbar(void){
  showcolorbar=0;
  if(glui_colorbar!=NULL){
    Reshape(screenWidth,screenHeight);
    ResetView(RESTORE_EXTERIOR_VIEW);
    glui_colorbar->hide();
  }
  updatemenu=1;
}

/* ------------------ show_glui_colorbar ------------------------ */

void show_glui_colorbar(void){
// show colorbar dialog box and redefine initial view point
  if(glui_colorbar!=NULL){
    Reshape(screenWidth,screenHeight);
    ResetView(RESTORE_EXTERIOR_VIEW);
    glui_colorbar->show();
  }
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

  panel_point = glui_colorbar->add_panel("Node");
  
  panel_cb5 = glui_colorbar->add_panel_to_panel(panel_point,"",GLUI_PANEL_NONE);

  BUTTON_prev=glui_colorbar->add_button_to_panel(panel_cb5,"Previous",COLORBAR_PREV,COLORBAR_CB);
  BUTTON_deletepoint=glui_colorbar->add_button_to_panel(panel_cb5,"Delete",COLORBAR_DELETEPOINT,COLORBAR_CB);
  
  glui_colorbar->add_column_to_panel(panel_cb5,false);

  BUTTON_next=glui_colorbar->add_button_to_panel(panel_cb5,"Next",COLORBAR_NEXT,COLORBAR_CB);
  BUTTON_addpoint=glui_colorbar->add_button_to_panel(panel_cb5,"Insert",COLORBAR_ADDPOINT,COLORBAR_CB);

  panel_cb4 = glui_colorbar->add_panel_to_panel(panel_point,"",GLUI_PANEL_NONE);
  SPINNER_colorindex=  glui_colorbar->add_spinner_to_panel(panel_cb4,"node index",  GLUI_SPINNER_INT,&cb_colorindex,  COLORBAR_COLORINDEX,COLORBAR_CB);
  SPINNER_colorindex->set_int_limits(0,255);
  SPINNER_right_red=  glui_colorbar->add_spinner_to_panel(panel_cb4,"red",  GLUI_SPINNER_INT,cb_rgb,COLORBAR_RGB,COLORBAR_CB);
  SPINNER_right_green=glui_colorbar->add_spinner_to_panel(panel_cb4,"green",GLUI_SPINNER_INT,cb_rgb+1,COLORBAR_RGB,COLORBAR_CB);
  SPINNER_right_blue= glui_colorbar->add_spinner_to_panel(panel_cb4,"blue", GLUI_SPINNER_INT,cb_rgb+2,COLORBAR_RGB,COLORBAR_CB);

  SPINNER_right_red->set_int_limits(0,255);
  SPINNER_right_green->set_int_limits(0,255);
  SPINNER_right_blue->set_int_limits(0,255);

  panel_extreme = glui_colorbar->add_panel("");

  CHECKBOX_usebounds=glui_colorbar->add_checkbox_to_panel(panel_extreme,"Highlight extreme data",&show_extremedata,
    COLORBAR_EXTREME,COLORBAR_CB);
  panel_cb9 = glui_colorbar->add_panel_to_panel(panel_extreme,"",GLUI_PANEL_NONE);
  panel_cb8 = glui_colorbar->add_panel_to_panel(panel_cb9,"Below specified min");
  SPINNER_down_red=  glui_colorbar->add_spinner_to_panel(panel_cb8,"red",  GLUI_SPINNER_INT,cb_down_rgb,COLORBAR_EXTREME_RGB,COLORBAR_CB);
  SPINNER_down_green=glui_colorbar->add_spinner_to_panel(panel_cb8,"green",GLUI_SPINNER_INT,cb_down_rgb+1,COLORBAR_EXTREME_RGB,COLORBAR_CB);
  SPINNER_down_blue= glui_colorbar->add_spinner_to_panel(panel_cb8,"blue", GLUI_SPINNER_INT,cb_down_rgb+2,COLORBAR_EXTREME_RGB,COLORBAR_CB);

  glui_colorbar->add_column_to_panel(panel_cb9);

  panel_cb7 = glui_colorbar->add_panel_to_panel(panel_cb9,"Above specified max");
  SPINNER_up_red=  glui_colorbar->add_spinner_to_panel(panel_cb7,"red",  GLUI_SPINNER_INT,cb_up_rgb,COLORBAR_EXTREME_RGB,COLORBAR_CB);
  SPINNER_up_green=glui_colorbar->add_spinner_to_panel(panel_cb7,"green",GLUI_SPINNER_INT,cb_up_rgb+1,COLORBAR_EXTREME_RGB,COLORBAR_CB);
  SPINNER_up_blue= glui_colorbar->add_spinner_to_panel(panel_cb7,"blue", GLUI_SPINNER_INT,cb_up_rgb+2,COLORBAR_EXTREME_RGB,COLORBAR_CB);

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
  unsigned char *rgb_nodes;
  int i;

  switch (var){
  case COLORBAR_COLORINDEX:
    if(colorbartype>=ndefaultcolorbars&&colorbartype<ncolorbars){
      cbi = colorbarinfo + colorbartype;
      
      cbi->index_node[colorbarpoint]=cb_colorindex;

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
      char *clabel;

      cbi = colorbarinfo + colorbartype;
      clabel=EDITTEXT_colorbar_label->get_text();
      strcpy(cbi->label,clabel);
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
    if(colorbarpoint<=0||colorbarpoint>cbi->nnodes-1)return;

    cbi->nnodes++;
    if(colorbarpoint<1)colorbarpoint=1;
    cbi->nodehilight=colorbarpoint;
    for(i=cbi->nnodes-1;i>=colorbarpoint+1;i--){
      unsigned char *rgb1, *rgb2;

      rgb2 = cbi->rgb_node+3*i;
      rgb1 = rgb2-3;
      rgb2[0] =rgb1[0];
      rgb2[1] =rgb1[1];
      rgb2[2] =rgb1[2];
      cbi->index_node[i]=cbi->index_node[i-1];
    }
    {
      unsigned char *rnew, *rbef, *raft;
      unsigned char *inew, *ibef, *iaft;

      rnew = cbi->rgb_node+3*colorbarpoint;
      rbef = rnew-3;
      raft = rnew+3;
      rnew[0]=((float)rbef[0]+(float)raft[0])/2.0;
      rnew[1]=((float)rbef[1]+(float)raft[1])/2.0;
      rnew[2]=((float)rbef[2]+(float)raft[2])/2.0;

      inew = cbi->index_node+colorbarpoint;
      ibef = inew-1;
      iaft = inew+1;
      *inew = (*ibef+*iaft)/2;
    }

    colorbar_global2local();
    remapcolorbar(cbi);
    updatecolors(-1);

    if(colorbarpoint==cbi->nnodes)colorbarpoint=cbi->nnodes-1;
    break;
  case COLORBAR_DELETEPOINT:
    if(colorbartype<ndefaultcolorbars||colorbartype>=ncolorbars)return;
    cbi = colorbarinfo + colorbartype;
    if(colorbarpoint<0||colorbarpoint>cbi->nnodes-1)return;

    if(cbi->nnodes<2)return;
    for(i=colorbarpoint+1;i<cbi->nnodes;i++){
      unsigned char *rgb1, *rgb2;

      cbi->index_node[i-1]=cbi->index_node[i];
      rgb2 = cbi->rgb_node+3*i;
      rgb1 = rgb2-3;
      rgb1[0]=rgb2[0];
      rgb1[1]=rgb2[1];
      rgb1[2]=rgb2[2];
    }
    cbi->nnodes--;
    remapcolorbar(cbi);
    updatecolors(-1);
    if(colorbarpoint==cbi->nnodes)colorbarpoint=cbi->nnodes-1;
    break;
  case COLORBAR_EXTREME:
    if(show_extremedata==1){
      SPINNER_down_red->enable();
      SPINNER_down_green->enable();
      SPINNER_down_blue->enable();
      SPINNER_up_red->enable();
      SPINNER_up_green->enable();
      SPINNER_up_blue->enable();
    }
    else{
      SPINNER_down_red->disable();
      SPINNER_down_green->disable();
      SPINNER_down_blue->disable();
      SPINNER_up_red->disable();
      SPINNER_up_green->disable();
      SPINNER_up_blue->disable();
    }
    if(colorbartype<0||colorbartype>=ncolorbars)return;
    cbi = colorbarinfo + colorbartype;
    remapcolorbar(cbi);
    updatecolors(-1);
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
    updatecolors(-1);
    break;
  case COLORBAR_RGB:
    if(colorbartype<0||colorbartype>=ncolorbars)return;
    cbi = colorbarinfo + colorbartype;
    if(colorbarpoint<0||colorbarpoint>cbi->nnodes-1)return;
    rgb_nodes=cbi->rgb_node+3*colorbarpoint;

    for(i=0;i<3;i++){
      rgb_nodes[i]=cb_rgb[i];
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
      if(colorbarpoint>cbi->nnodes-1)colorbarpoint=0;
    }
    else if(var==COLORBAR_PREV){
      colorbarpoint--;
      if(colorbarpoint<0)colorbarpoint=cbi->nnodes-1;
    }
    cbi->nodehilight=colorbarpoint;
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
  unsigned char *rgb;
  int icolorbar;

  if(colorbartype<0||colorbartype>=ncolorbars)return;

  cbi = colorbarinfo + colorbartype;
  colorbarpoint=cbi->nodehilight;
    
  SPINNER_colorindex->set_int_val(cbi->index_node[colorbarpoint]);

  BUTTON_next->enable();
  BUTTON_prev->enable();

  strcpy(colorbar_label,cbi->label);
  EDITTEXT_colorbar_label->set_text(colorbar_label);
  icolorbar=LISTBOX_colorbar->get_int_val();
  if(icolorbar>=ndefaultcolorbars){
    BUTTON_delete->enable();
    EDITTEXT_colorbar_label->enable();
    SPINNER_right_red->enable();
    SPINNER_right_green->enable();
    SPINNER_right_blue->enable();
    BUTTON_addpoint->enable();
    BUTTON_deletepoint->enable();
    SPINNER_colorindex->enable();
  }
  else{
    BUTTON_delete->disable();
    EDITTEXT_colorbar_label->disable();
    SPINNER_right_red->disable();
    SPINNER_right_green->disable();
    SPINNER_right_blue->disable();
    BUTTON_addpoint->disable();
    BUTTON_deletepoint->disable();
    SPINNER_colorindex->disable();
  }
  rgb = cbi->rgb_node+3*colorbarpoint;
  SPINNER_right_red->set_int_val(  (int)(rgb[0]));
  SPINNER_right_green->set_int_val((int)(rgb[1]));
  SPINNER_right_blue->set_int_val( (int)(rgb[2]));

  rgb = rgb_below_min;
  SPINNER_down_red->set_int_val(  (int)(rgb[0]));
  SPINNER_down_green->set_int_val(  (int)(rgb[1]));
  SPINNER_down_blue->set_int_val(  (int)(rgb[2]));

  rgb = rgb_above_max;
  SPINNER_up_red->set_int_val(  (int)(rgb[0]));
  SPINNER_up_green->set_int_val(  (int)(rgb[1]));
  SPINNER_up_blue->set_int_val(  (int)(rgb[2]));

  COLORBAR_CB(COLORBAR_EXTREME);


}
