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
#include "MALLOC.h"
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"
#include "translate.h"
#include "MALLOC.h"

// svn revision character string
extern "C" char glui_vol3dsmoke_revision[]="$Revision$";

extern "C" void SMOKE_3D_CB(int var);
extern "C" void update_gpu(void);
extern GLUI_Rollout *panel_volsmoke3d;
void VOLSMOKE_3D_CB(int var);
GLUI_Panel *panel_volslices=NULL, *panel_volcolormap=NULL, *panel_voloverall=NULL;
GLUI_Spinner *SPINNERvol_smoke3d_smoke_shade=NULL;
GLUI_Listbox *LISTBOXvol_smoke_colorbar=NULL;
GLUI *glui_vol3dsmoke=NULL;
GLUI_Checkbox *CHECKBOXvol_usevolrender=NULL;
GLUI_Checkbox *CHECKBOXvol_use_firesmokemap=NULL;
GLUI_Checkbox *CHECKBOXvol_smokeGPU=NULL;

#define SMOKE_COLORBAR_LIST 16
#define USE_FIRESMOKEMAP 17
#define SMOKE_SHADE 7
#define VOL_SMOKE 13

extern GLUI *glui_bounds;

/* ------------------ update_volgpu ------------------------ */

extern "C" void update_volgpu(void){
#ifdef pp_GPU
  if(nvolrenderinfo>0&&CHECKBOXvol_smokeGPU!=NULL){
    CHECKBOXvol_smokeGPU->set_int_val(usegpu);  
  }
#endif
}

/* ------------------ hide_glui_3dsmoke ------------------------ */

extern "C" void hide_glui_vol3dsmoke(void){
  if(glui_vol3dsmoke!=NULL)glui_vol3dsmoke->hide();
  showgluivol3dsmoke=0;
  updatemenu=1;
}

/* ------------------ show_glui_vol3dsmoke ------------------------ */

extern "C" void show_glui_vol3dsmoke(void){
  if(glui_vol3dsmoke!=NULL)glui_vol3dsmoke->show();
}

/* ------------------ update_volsmoke3dflags ------------------------ */

extern "C" void update_volsmoke3dflags(void){
  if(CHECKBOXvol_usevolrender!=NULL)CHECKBOXvol_usevolrender->set_int_val(usevolrender);
  glutPostRedisplay();
}

/* ------------------ glui_3dsmoke_setup ------------------------ */

extern "C" void glui_vol3dsmoke_setup(int main_window){

  int i;

  if(nvolrenderinfo<=0)return;
  
  glui_vol3dsmoke=glui_bounds;

  panel_voloverall = glui_vol3dsmoke->add_panel_to_panel(panel_volsmoke3d,"",GLUI_PANEL_NONE);
#ifdef pp_GPU
  CHECKBOXvol_smokeGPU=glui_vol3dsmoke->add_checkbox_to_panel(panel_voloverall,_("Use GPU"),&usegpu,VOL_SMOKE,VOLSMOKE_3D_CB);
#endif
  if(ncolorbars>0){
    LISTBOXvol_smoke_colorbar=glui_vol3dsmoke->add_listbox_to_panel(panel_voloverall,_("colormap:"),&fire_colorbar_index,SMOKE_COLORBAR_LIST,VOLSMOKE_3D_CB);

    for(i=0;i<ncolorbars;i++){
      colorbardata *cbi;

      cbi = colorbarinfo + i;
      cbi->label_ptr=cbi->label;
      LISTBOXvol_smoke_colorbar->add_item(i,cbi->label_ptr);
    }
    LISTBOXvol_smoke_colorbar->set_int_val(fire_colorbar_index);
  }

  SPINNERvol_smoke3d_smoke_shade=glui_vol3dsmoke->add_spinner_to_panel(panel_voloverall,_("smoke albedo"),GLUI_SPINNER_FLOAT,&smoke_shade,SMOKE_SHADE,VOLSMOKE_3D_CB);
  SPINNERvol_smoke3d_smoke_shade->set_float_limits(0.0,1.0);


#ifdef _DEBUG
  CHECKBOXvol_usevolrender=glui_vol3dsmoke->add_checkbox_to_panel(panel_voloverall,_("Use full volume rendering"),&usevolrender,VOL_SMOKE,VOLSMOKE_3D_CB);
#endif
  glui_vol3dsmoke->add_checkbox_to_panel(panel_voloverall,_("Load data in background"),&use_multi_threading);
  glui_vol3dsmoke->add_checkbox_to_panel(panel_voloverall,_("Load data only at render times"),&load_at_rendertimes);
#ifdef _DEBUG
  glui_vol3dsmoke->add_checkbox_to_panel(panel_voloverall,"block smoke",&block_volsmoke);
  glui_vol3dsmoke->add_checkbox_to_panel(panel_voloverall,"debug",&smoke3dVoldebug);
#endif

}

/* ------------------ VOLSMOKE_3D_CB ------------------------ */

void VOLSMOKE_3D_CB(int var){
  switch (var){
    case USE_FIRESMOKEMAP:
      break;
  case VOL_SMOKE:
#ifdef pp_GPU
    SMOKE_3D_CB(VOL_SMOKE);
    update_gpu();
    updatemenu=1;
#endif
    break;
  }
}

