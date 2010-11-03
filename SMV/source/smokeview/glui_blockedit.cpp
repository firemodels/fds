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
#include "glui.h"
#include "flowfiles.h"
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

// svn revision character string
extern "C" char glui_blockedit_revision[]="$Revision$";

#define XMIN_SPIN 20
#define YMIN_SPIN 21
#define ZMIN_SPIN 22
#define XMAX_SPIN 23
#define YMAX_SPIN 24
#define ZMAX_SPIN 25
#define UPDATE_LIST 31
#define RADIO_WALL 32

GLUI_StaticText *statictext_blockage_index=NULL;
GLUI_StaticText *statictext_mesh_index=NULL;
GLUI_StaticText *statictext_label=NULL;

GLUI *glui_edit=NULL;
GLUI_Panel *panel_obj_select=NULL,*panel_surface=NULL;
GLUI_Panel *panel_obj_stretch2=NULL,*panel_obj_stretch3=NULL, *panel_obj_stretch4=NULL;
GLUI_EditText *edittext_xmin=NULL, *edittext_ymin=NULL, *edittext_zmin=NULL;
GLUI_EditText *edittext_xmax=NULL, *edittext_ymax=NULL, *edittext_zmax=NULL;
GLUI_Listbox *surfacelists[7]={NULL,NULL,NULL,NULL,NULL,NULL,NULL};
GLUI_Checkbox *blockage_checkbox=NULL;

extern "C" void OBJECT_CB(int var);

void BUTTON_hide3_CB(int var);

char a_updatelabel[1000];
char *updatelabel=NULL;

/* ------------------ glui_edit_setup ------------------------ */

extern "C" void glui_edit_setup(int main_window){
  int ibar,jbar,kbar;
  float *xplt_orig, *yplt_orig, *zplt_orig;
  surface *surfi;
  char *surfacelabel;
  int i;
  mesh *meshi;

  ibar=current_mesh->ibar;
  jbar=current_mesh->jbar;
  kbar=current_mesh->kbar;
  xplt_orig=current_mesh->xplt_orig;
  yplt_orig=current_mesh->yplt_orig;
  zplt_orig=current_mesh->zplt_orig;

  if(glui_edit!=NULL)glui_edit->close();
  glui_edit = GLUI_Master.create_glui("Blockage Info",0,0,0);
  if(showedit==0)glui_edit->hide();

  panel_obj_select = glui_edit->add_panel("SURFs");

  panel_surface=glui_edit->add_panel_to_panel(panel_obj_select,"",GLUI_PANEL_NONE);

  glui_edit->add_column_to_panel(panel_surface,false);

  if(nsurfaces>0){
    int not_used=0;

    glui_edit->add_statictext_to_panel(panel_surface,"");

    surfacelists[NOT_USED] = glui_edit->add_listbox_to_panel(panel_surface,"Unused SURFs",
      surface_indices+NOT_USED,UPDATE_LIST,OBJECT_CB);
    surfacelists[NOT_USED]->set_w(260);
    for(i=0;i<nsurfaces;i++){
      surfi = surfaceinfo + sorted_surfidlist[i];
      if(surfi->used_by_obst==1)continue;
      if(surfi->used_by_vent==1)continue;
      if(surfi->obst_surface==0)continue;
      not_used++;
    }
    if(not_used>0){
      for(i=0;i<nsurfaces;i++){
        surfi = surfaceinfo + sorted_surfidlist[i];
        if(surfi->used_by_obst==1)continue;
        if(surfi->used_by_vent==1)continue;
        if(surfi->obst_surface==0)continue;
        surfacelabel = surfi->surfacelabel;
        surfacelists[NOT_USED]->add_item(i,surfacelabel);
        surface_indices[NOT_USED]=i;
        surface_indices_bak[NOT_USED]=surface_indices[NOT_USED];
      }
    }
    else{
      surfacelists[NOT_USED]->add_item(0,"none");
      surface_indices[NOT_USED]=0;
      surface_indices_bak[NOT_USED]=surface_indices[NOT_USED];
    }
    surfacelists[NOT_USED]->set_int_val(surface_indices[NOT_USED]);

    surfacelists[DOWN_X] = glui_edit->add_listbox_to_panel(panel_surface,"Left",surface_indices+DOWN_X,UPDATE_LIST,OBJECT_CB);
    surfacelists[DOWN_X]->set_w(260);
    for(i=0;i<nsurfaces;i++){
      surfi = surfaceinfo + sorted_surfidlist[i];
      if(surfi->used_by_obst!=1)continue;
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      surfacelists[DOWN_X]->add_item(i,surfacelabel);
    }

    surfacelists[UP_X] = glui_edit->add_listbox_to_panel(panel_surface,"Right",surface_indices+UP_X,UPDATE_LIST,OBJECT_CB);
    surfacelists[UP_X]->set_w(260);
    for(i=0;i<nsurfaces;i++){
      surfi = surfaceinfo + sorted_surfidlist[i];
      if(surfi->used_by_obst!=1)continue;
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      surfacelists[UP_X]->add_item(i,surfacelabel);
    }

    surfacelists[DOWN_Y] = glui_edit->add_listbox_to_panel(panel_surface,"Front",surface_indices+DOWN_Y,UPDATE_LIST,OBJECT_CB);
    surfacelists[DOWN_Y]->set_w(260);
    for(i=0;i<nsurfaces;i++){
      surfi = surfaceinfo + sorted_surfidlist[i];
      if(surfi->used_by_obst!=1)continue;
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      surfacelists[DOWN_Y]->add_item(i,surfacelabel);
    }

    surfacelists[UP_Y] = glui_edit->add_listbox_to_panel(panel_surface,"Back",surface_indices+UP_Y,UPDATE_LIST,OBJECT_CB);
    surfacelists[UP_Y]->set_w(260);
    for(i=0;i<nsurfaces;i++){
      surfi = surfaceinfo + sorted_surfidlist[i];
      if(surfi->used_by_obst!=1)continue;
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      surfacelists[UP_Y]->add_item(i,surfacelabel);
    }

    surfacelists[DOWN_Z] = glui_edit->add_listbox_to_panel(panel_surface,"Down",surface_indices+DOWN_Z,UPDATE_LIST,OBJECT_CB);
    surfacelists[DOWN_Z]->set_w(260);
    for(i=0;i<nsurfaces;i++){
      surfi = surfaceinfo + sorted_surfidlist[i];
      if(surfi->used_by_obst!=1)continue;
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      surfacelists[DOWN_Z]->add_item(i,surfacelabel);
    }

    surfacelists[UP_Z] = glui_edit->add_listbox_to_panel(panel_surface,"Up",surface_indices+UP_Z,UPDATE_LIST,OBJECT_CB);
    surfacelists[UP_Z]->set_w(260);
    for(i=0;i<nsurfaces;i++){
      surfi = surfaceinfo + sorted_surfidlist[i];
      if(surfi->used_by_obst!=1)continue;
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      surfacelists[UP_Z]->add_item(i,surfacelabel);
    }

    OBJECT_CB(RADIO_WALL);
    for(i=0;i<6;i++){
      surfacelists[i]->disable();
    }
  }

  glui_edit->add_column(false);

  panel_obj_stretch2 = glui_edit->add_panel("Coordinates");

  blockage_checkbox=glui_edit->add_checkbox_to_panel(panel_obj_stretch2,"Dimensions snapped to grid",&blockage_snapped,
    BLOCKAGE_AS_INPUT,OBJECT_CB);
  panel_obj_stretch3 = glui_edit->add_panel_to_panel(panel_obj_stretch2,"",GLUI_PANEL_NONE);
  edittext_xmin=glui_edit->add_edittext_to_panel(panel_obj_stretch3,"x",GLUI_EDITTEXT_FLOAT,&glui_block_xmin,XMIN_SPIN,OBJECT_CB);
  edittext_ymin=glui_edit->add_edittext_to_panel(panel_obj_stretch3,"y",GLUI_EDITTEXT_FLOAT,&glui_block_ymin,YMIN_SPIN,OBJECT_CB);
  edittext_zmin=glui_edit->add_edittext_to_panel(panel_obj_stretch3,"z",GLUI_EDITTEXT_FLOAT,&glui_block_zmin,ZMIN_SPIN,OBJECT_CB);

  glui_edit->add_column_to_panel(panel_obj_stretch3,false);
  edittext_xmax=glui_edit->add_edittext_to_panel(panel_obj_stretch3,"",GLUI_EDITTEXT_FLOAT,&glui_block_xmax,XMAX_SPIN,OBJECT_CB);
  edittext_ymax=glui_edit->add_edittext_to_panel(panel_obj_stretch3,"",GLUI_EDITTEXT_FLOAT,&glui_block_ymax,YMAX_SPIN,OBJECT_CB);
  edittext_zmax=glui_edit->add_edittext_to_panel(panel_obj_stretch3,"",GLUI_EDITTEXT_FLOAT,&glui_block_zmax,ZMAX_SPIN,OBJECT_CB);

  edittext_xmin->disable();
  edittext_ymin->disable();
  edittext_zmin->disable();

  edittext_xmax->disable();
  edittext_ymax->disable();
  edittext_zmax->disable();
  OBJECT_CB(BLOCKAGE_AS_INPUT);

  edittext_xmin->set_float_limits(xplt_orig[0],xplt_orig[ibar],GLUI_LIMIT_CLAMP);
  edittext_xmax->set_float_limits(xplt_orig[0],xplt_orig[ibar],GLUI_LIMIT_CLAMP);
  edittext_ymin->set_float_limits(yplt_orig[0],yplt_orig[jbar],GLUI_LIMIT_CLAMP);
  edittext_ymax->set_float_limits(yplt_orig[0],yplt_orig[jbar],GLUI_LIMIT_CLAMP);
  edittext_zmin->set_float_limits(zplt_orig[0],zplt_orig[kbar],GLUI_LIMIT_CLAMP);
  edittext_zmax->set_float_limits(zplt_orig[0],zplt_orig[kbar],GLUI_LIMIT_CLAMP);

  panel_obj_stretch4=glui_edit->add_panel("",GLUI_PANEL_NONE);

  if(nmeshes>1){
    char meshlabel[255];

    strcpy(meshlabel,"Mesh label:");
    strcat(meshlabel,meshinfo->label);
    statictext_mesh_index=glui_edit->add_statictext_to_panel(panel_obj_stretch4,meshlabel);
  }
  statictext_blockage_index=glui_edit->add_statictext_to_panel(panel_obj_stretch4,"&OBST number: ");
  statictext_label=glui_edit->add_statictext_to_panel(panel_obj_stretch4,"&OBST label:");
  glui_edit->add_separator_to_panel(panel_obj_stretch4);
  glui_edit->add_button_to_panel(panel_obj_stretch4,"Close",CLOSE_WINDOW,BUTTON_hide3_CB);

  glui_edit->set_main_gfx_window( main_window );
}

/* ------------------ hide_glui_edit ------------------------ */

extern "C" void hide_glui_edit(void){
  blockageSelect=0;
  if(glui_edit!=NULL)glui_edit->hide();
  showedit=0;
  updatemenu=1;
  editwindow_status=CLOSE_WINDOW;
}

/* ------------------ show_glui_edit ------------------------ */

extern "C" void show_glui_edit(void){
  blockageSelect=1;
  update_blockvals(0);
  glui_edit->show();
}

extern "C" void DialogMenu(int value);


/* ------------------ BUTTON_hide3_CB ------------------------ */

void BUTTON_hide3_CB(int var){
  blockagedata *bc;
  blockagedata *bchighlight_save;
  int i,j,k;
  mesh *meshi;
  switch (var){
  case CLOSE_WINDOW: 
    DialogMenu(16);
    smooth_blockages();
    break;
  case UPDATE_WINDOW:
    blockages_dirty=0;
    break;
  case CANCEL_WINDOW:
    break;
  default:
    ASSERT(FFALSE);
    break;
  }

}

/* ------------------ update_blockvals ------------------------ */

extern "C" void update_blockvals(int flag){
  float xmin, xmax, ymin, ymax, zmin, zmax;
  int imin, jmin, kmin;
  char *label;
  int i;
  int temp;
  float *xplt_orig, *yplt_orig, *zplt_orig;
  int ibar, jbar, kbar;

  get_blockvals(&xmin,&xmax,&ymin,&ymax,&zmin,&zmax,&imin,&jmin,&kmin);

  xplt_orig = current_mesh->xplt_orig;
  yplt_orig = current_mesh->yplt_orig;
  zplt_orig = current_mesh->zplt_orig;
  ibar = current_mesh->ibar;
  jbar = current_mesh->jbar;
  kbar = current_mesh->kbar;

  edittext_xmin->set_float_limits(xplt_orig[0],xplt_orig[ibar],GLUI_LIMIT_CLAMP);
  edittext_xmax->set_float_limits(xplt_orig[0],xplt_orig[ibar],GLUI_LIMIT_CLAMP);
  edittext_ymin->set_float_limits(yplt_orig[0],yplt_orig[jbar],GLUI_LIMIT_CLAMP);
  edittext_ymax->set_float_limits(yplt_orig[0],yplt_orig[jbar],GLUI_LIMIT_CLAMP);
  edittext_zmin->set_float_limits(zplt_orig[0],zplt_orig[kbar],GLUI_LIMIT_CLAMP);
  edittext_zmax->set_float_limits(zplt_orig[0],zplt_orig[kbar],GLUI_LIMIT_CLAMP);

  edittext_xmin->set_float_val(xmin);
  edittext_xmax->set_float_val(xmax);
  edittext_ymin->set_float_val(ymin);
  edittext_ymax->set_float_val(ymax);
  edittext_zmin->set_float_val(zmin);
  edittext_zmax->set_float_val(zmax);
  if(bchighlight!=NULL&&nsurfaces>0){
    wall_case=bchighlight->walltype;
#ifdef pp_WALLGROUP
    wallgroup->set_int_val(wall_case);
#endif
    OBJECT_CB(RADIO_WALL);
  }

  if(flag==1){
    if(bchighlight!=NULL){
      char dialog_label[255];
      mesh *blockmesh;

      if(nmeshes>1){
        blockmesh = meshinfo + bchighlight->meshindex;
        sprintf(dialog_label,"Mesh label: %s",blockmesh->label);
        statictext_mesh_index->set_text(dialog_label);
      }
      sprintf(dialog_label,"&OBST index: %i",bchighlight->id);
      statictext_blockage_index->set_text(dialog_label);
      strcpy(dialog_label,"&OBST label: ");
      strcat(dialog_label,bchighlight->label);
      statictext_label->set_text(dialog_label);

      switch (wall_case){
      case WALL_1:
        temp=bchighlight->surf_index[UP_Z];
        for(i=0;i<6;i++){
          bchighlight->surf_index[i]=temp;
        }
        break;
      case WALL_3:
        temp=bchighlight->surf_index[UP_Y];
        bchighlight->surf_index[DOWN_X]=temp;
        bchighlight->surf_index[DOWN_Y]=temp;
        bchighlight->surf_index[UP_X]=temp;
        break;
      case WALL_6:
        break;
      default:
        ASSERT(FFALSE);
        break;
      }

      if(nsurfaces>0){
        for(i=0;i<6;i++){
          surface_indices[i] = inv_sorted_surfidlist[bchighlight->surf_index[i]];
          surface_indices_bak[i] = inv_sorted_surfidlist[bchighlight->surf_index[i]];
          surfacelists[i]->set_int_val(surface_indices[i]);
        }
      }
    }
    else{
      if(nsurfaces>0){
        for(i=0;i<6;i++){
          surface_indices[i]=inv_sorted_surfidlist[0];
          surface_indices_bak[i]=inv_sorted_surfidlist[0];
          surfacelists[i]->set_int_val(surface_indices[i]);
        }
      }
    }
  }
}

/* ------------------ OBJECT_CB ------------------------ */

void OBJECT_CB(int var){
  int i,temp;
  switch (var){
    case UPDATE_LIST:
      if(surface_indices[NOT_USED]!=surface_indices_bak[NOT_USED]){
        surface_indices[NOT_USED]=surface_indices_bak[NOT_USED];
        surfacelists[NOT_USED]->set_int_val(surface_indices[NOT_USED]);
      }
      switch (wall_case){
      case WALL_1:
        temp=surface_indices_bak[UP_Z];
        if(nsurfaces>0){
          for(i=0;i<6;i++){
            surface_indices[i]=temp;
            surfacelists[i]->set_int_val(temp);
          }
        }
        break;
      case WALL_3:
        if(nsurfaces>0){
          for(i=0;i<6;i++){
            temp=surface_indices_bak[i];
            surface_indices[i]=temp;
            surfacelists[i]->set_int_val(temp);
          }
        }
        break;
      case WALL_6:
        if(nsurfaces>0){
          for(i=0;i<6;i++){
            temp=surface_indices_bak[i];
            surface_indices[i]=temp;
            surfacelists[i]->set_int_val(temp);
          }
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
      }

      if(bchighlight!=NULL){
        for(i=0;i<6;i++){
          bchighlight->surf[i]=surfaceinfo+sorted_surfidlist[surface_indices_bak[i]];
          bchighlight->surf_index[i]=sorted_surfidlist[surface_indices_bak[i]];
        }
        bchighlight->changed_surface=1;
        if(bchighlight->id>0&&bchighlight->id<=nchanged_idlist){
          changed_idlist[bchighlight->id]=1;
        }
        blockages_dirty=1;
        updateusetextures();
        update_faces();
      }
      break;
    case RADIO_WALL:
      if(nsurfaces==0)break;
      if(bchighlight!=NULL){
        bchighlight->walltype=wall_case;
      }
      switch (wall_case){
      case WALL_6:
        for(i=0;i<6;i++){
          surfacelists[i]->enable();
        }
        surfacelists[DOWN_Z]->set_name("z lower face");
        surfacelists[UP_Z]->set_name("z upper face");
        surfacelists[DOWN_Y]->set_name("y lower face");
        surfacelists[UP_Y]->set_name("y upper face");
        surfacelists[DOWN_X]->set_name("x lower face");
        surfacelists[UP_X]->set_name("x upper face");
        break;
      case WALL_3:
        for(i=0;i<6;i++){
          surfacelists[i]->disable();
        }
        surfacelists[DOWN_Z]->enable();
        surfacelists[UP_Z]->enable();
        surfacelists[UP_Y]->enable();

        surfacelists[DOWN_Z]->set_name("z lower face");
        surfacelists[UP_Z]->set_name("z upper face");
        surfacelists[UP_Y]->set_name("side faces");
        surfacelists[DOWN_Y]->set_name("");
        surfacelists[DOWN_X]->set_name("");
        surfacelists[UP_X]->set_name("");

        break;
      case WALL_1:
        for(i=0;i<6;i++){
          surfacelists[i]->disable();
        }
        surfacelists[UP_Z]->enable();
        surfacelists[UP_Z]->set_name("All faces");

        surfacelists[DOWN_Z]->set_name("");
        surfacelists[DOWN_Y]->set_name("");
        surfacelists[UP_Y]->set_name("");
        surfacelists[DOWN_X]->set_name("");
        surfacelists[UP_X]->set_name("");
        break;
      default:
        ASSERT(FFALSE);
        break;
      }
      OBJECT_CB(UPDATE_LIST);
      break;
      case BLOCKAGE_AS_INPUT2:
        blockage_snapped=1-blockage_as_input;
        blockage_checkbox->set_int_val(blockage_snapped);
      case BLOCKAGE_AS_INPUT:
        blockage_as_input=1-blockage_snapped;
        if(blockage_as_input==1){
          blocklocation=BLOCKlocation_exact;
        }
        else{
          blocklocation=BLOCKlocation_grid;
        }
        update_blockvals(0);
        break;
    default:
      ASSERT(FFALSE);
      break;
  }
}
