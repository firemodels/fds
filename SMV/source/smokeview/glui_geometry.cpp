// $Date$ 
// $Revision$
// $Author$

#define CPP
#include "options.h"

// svn revision character string
extern "C" char glui_blockedit_revision[];
char glui_blockedit_revision[]="$Revision$";


#include <stdio.h>
#include <string.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "smokeviewvars.h"

#define XMIN_SPIN 20
#define YMIN_SPIN 21
#define ZMIN_SPIN 22
#define XMAX_SPIN 23
#define YMAX_SPIN 24
#define ZMAX_SPIN 25
#define UPDATE_LIST 31
#define RADIO_WALL 32
#define SAVE_SETTINGS 33

GLUI *glui_geometry=NULL;

GLUI_Button *BUTTON_blockage_1=NULL;

GLUI_Checkbox *CHECKBOX_blockage=NULL;

GLUI_EditText *EDIT_xmin=NULL, *EDIT_ymin=NULL, *EDIT_zmin=NULL;
GLUI_EditText *EDIT_xmax=NULL, *EDIT_ymax=NULL, *EDIT_zmax=NULL;

GLUI_Listbox *LIST_surface[7]={NULL,NULL,NULL,NULL,NULL,NULL,NULL};

GLUI_Panel *PANEL_obj_select=NULL,*PANEL_surface=NULL;
GLUI_Panel *PANEL_obj_stretch2=NULL,*PANEL_obj_stretch3=NULL, *PANEL_obj_stretch4=NULL;

GLUI_Rollout *ROLLOUT_blockedit=NULL;
GLUI_Rollout *ROLLOUT_geometry=NULL;

GLUI_Spinner *SPINNER_face_factor=NULL;

GLUI_StaticText *STATIC_blockage_index=NULL;
GLUI_StaticText *STATIC_mesh_index=NULL;
GLUI_StaticText *STATIC_label=NULL;

void Blockedit_DLG_CB(int var);

char a_updatelabel[1000];
char *updatelabel=NULL;

/* ------------------ glui_geometry_setup ------------------------ */

extern "C" void glui_geometry_setup(int main_window){
  int ibar,jbar,kbar;
  float *xplt_orig, *yplt_orig, *zplt_orig;
  surfdata *surfi;
  char *surfacelabel;
  int i;

  ibar=current_mesh->ibar;
  jbar=current_mesh->jbar;
  kbar=current_mesh->kbar;
  xplt_orig=current_mesh->xplt_orig;
  yplt_orig=current_mesh->yplt_orig;
  zplt_orig=current_mesh->zplt_orig;

  update_glui_geometry=0;
  if(glui_geometry!=NULL){
    glui_geometry->close();
    glui_geometry=NULL;
  }
  glui_geometry = GLUI_Master.create_glui("Geometry",0,0,0);
  if(showedit_dialog==0)glui_geometry->hide();

  ROLLOUT_blockedit = glui_geometry->add_rollout("Structured",false);
  PANEL_obj_select = glui_geometry->add_panel_to_panel(ROLLOUT_blockedit,"SURFs");

  PANEL_surface=glui_geometry->add_panel_to_panel(PANEL_obj_select,"",GLUI_PANEL_NONE);

  glui_geometry->add_column_to_panel(PANEL_surface,false);

  if(nsurfinfo>0){
    int not_used=0;

    glui_geometry->add_statictext_to_panel(PANEL_surface,"");

    LIST_surface[NOT_USED] = glui_geometry->add_listbox_to_panel(PANEL_surface,_("Unused SURFs"),surface_indices+NOT_USED,UPDATE_LIST,OBJECT_CB);
    LIST_surface[NOT_USED]->set_w(260);
    for(i=0;i<nsurfinfo;i++){
      surfi = surfinfo + sorted_surfidlist[i];
      if(surfi->used_by_obst==1)continue;
      if(surfi->used_by_vent==1)continue;
      if(surfi->obst_surface==0)continue;
      not_used++;
    }
    if(not_used>0){
      for(i=0;i<nsurfinfo;i++){
        surfi = surfinfo + sorted_surfidlist[i];
        if(surfi->used_by_obst==1)continue;
        if(surfi->used_by_vent==1)continue;
        if(surfi->obst_surface==0)continue;
        surfacelabel = surfi->surfacelabel;
        LIST_surface[NOT_USED]->add_item(i,surfacelabel);
        surface_indices[NOT_USED]=i;
        surface_indices_bak[NOT_USED]=surface_indices[NOT_USED];
      }
    }
    else{
      LIST_surface[NOT_USED]->add_item(0,_("None"));
      surface_indices[NOT_USED]=0;
      surface_indices_bak[NOT_USED]=surface_indices[NOT_USED];
    }
    LIST_surface[NOT_USED]->set_int_val(surface_indices[NOT_USED]);

    LIST_surface[DOWN_X] = glui_geometry->add_listbox_to_panel(PANEL_surface,_("Left"),surface_indices+DOWN_X,UPDATE_LIST,OBJECT_CB);
    LIST_surface[DOWN_X]->set_w(260);
    for(i=0;i<nsurfinfo;i++){
      surfi = surfinfo + sorted_surfidlist[i];
      if(surfi->used_by_obst!=1)continue;
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      LIST_surface[DOWN_X]->add_item(i,surfacelabel);
    }

    LIST_surface[UP_X] = glui_geometry->add_listbox_to_panel(PANEL_surface,_("Right"),surface_indices+UP_X,UPDATE_LIST,OBJECT_CB);
    LIST_surface[UP_X]->set_w(260);
    for(i=0;i<nsurfinfo;i++){
      surfi = surfinfo + sorted_surfidlist[i];
      if(surfi->used_by_obst!=1)continue;
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      LIST_surface[UP_X]->add_item(i,surfacelabel);
    }

    LIST_surface[DOWN_Y] = glui_geometry->add_listbox_to_panel(PANEL_surface,_("Front"),surface_indices+DOWN_Y,UPDATE_LIST,OBJECT_CB);
    LIST_surface[DOWN_Y]->set_w(260);
    for(i=0;i<nsurfinfo;i++){
      surfi = surfinfo + sorted_surfidlist[i];
      if(surfi->used_by_obst!=1)continue;
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      LIST_surface[DOWN_Y]->add_item(i,surfacelabel);
    }

    LIST_surface[UP_Y] = glui_geometry->add_listbox_to_panel(PANEL_surface,_("Back"),surface_indices+UP_Y,UPDATE_LIST,OBJECT_CB);
    LIST_surface[UP_Y]->set_w(260);
    for(i=0;i<nsurfinfo;i++){
      surfi = surfinfo + sorted_surfidlist[i];
      if(surfi->used_by_obst!=1)continue;
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      LIST_surface[UP_Y]->add_item(i,surfacelabel);
    }

    LIST_surface[DOWN_Z] = glui_geometry->add_listbox_to_panel(PANEL_surface,_("Down"),surface_indices+DOWN_Z,UPDATE_LIST,OBJECT_CB);
    LIST_surface[DOWN_Z]->set_w(260);
    for(i=0;i<nsurfinfo;i++){
      surfi = surfinfo + sorted_surfidlist[i];
      if(surfi->used_by_obst!=1)continue;
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      LIST_surface[DOWN_Z]->add_item(i,surfacelabel);
    }

    LIST_surface[UP_Z] = glui_geometry->add_listbox_to_panel(PANEL_surface,_("Up"),surface_indices+UP_Z,UPDATE_LIST,OBJECT_CB);
    LIST_surface[UP_Z]->set_w(260);
    for(i=0;i<nsurfinfo;i++){
      surfi = surfinfo + sorted_surfidlist[i];
      if(surfi->used_by_obst!=1)continue;
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      LIST_surface[UP_Z]->add_item(i,surfacelabel);
    }

    OBJECT_CB(RADIO_WALL);
    for(i=0;i<6;i++){
      LIST_surface[i]->disable();
    }
  }

  glui_geometry->add_column_to_panel(ROLLOUT_blockedit,false);

  PANEL_obj_stretch4=glui_geometry->add_panel_to_panel(ROLLOUT_blockedit,"",GLUI_PANEL_NONE);

  {
    char meshlabel[255];

    strcpy(meshlabel,_("Mesh:"));
    strcat(meshlabel,meshinfo->label);
    STATIC_mesh_index=glui_geometry->add_statictext_to_panel(PANEL_obj_stretch4,meshlabel);

  }
  STATIC_blockage_index=glui_geometry->add_statictext_to_panel(PANEL_obj_stretch4,"&OBST number: ");
  STATIC_label=glui_geometry->add_statictext_to_panel(PANEL_obj_stretch4,"&OBST label:");

  PANEL_obj_stretch2 = glui_geometry->add_panel_to_panel(ROLLOUT_blockedit,"Coordinates");

  CHECKBOX_blockage=glui_geometry->add_checkbox_to_panel(PANEL_obj_stretch2,_("Dimensions snapped to grid"),&blockage_snapped,
    BLOCKAGE_AS_INPUT,OBJECT_CB);
  PANEL_obj_stretch3 = glui_geometry->add_panel_to_panel(PANEL_obj_stretch2,"",GLUI_PANEL_NONE);
  EDIT_xmin=glui_geometry->add_edittext_to_panel(PANEL_obj_stretch3,"x",GLUI_EDITTEXT_FLOAT,&glui_block_xmin,XMIN_SPIN,OBJECT_CB);
  EDIT_ymin=glui_geometry->add_edittext_to_panel(PANEL_obj_stretch3,"y",GLUI_EDITTEXT_FLOAT,&glui_block_ymin,YMIN_SPIN,OBJECT_CB);
  EDIT_zmin=glui_geometry->add_edittext_to_panel(PANEL_obj_stretch3,"z",GLUI_EDITTEXT_FLOAT,&glui_block_zmin,ZMIN_SPIN,OBJECT_CB);

  glui_geometry->add_column_to_panel(PANEL_obj_stretch3,false);
  EDIT_xmax=glui_geometry->add_edittext_to_panel(PANEL_obj_stretch3,"",GLUI_EDITTEXT_FLOAT,&glui_block_xmax,XMAX_SPIN,OBJECT_CB);
  EDIT_ymax=glui_geometry->add_edittext_to_panel(PANEL_obj_stretch3,"",GLUI_EDITTEXT_FLOAT,&glui_block_ymax,YMAX_SPIN,OBJECT_CB);
  EDIT_zmax=glui_geometry->add_edittext_to_panel(PANEL_obj_stretch3,"",GLUI_EDITTEXT_FLOAT,&glui_block_zmax,ZMAX_SPIN,OBJECT_CB);

  EDIT_xmin->disable();
  EDIT_ymin->disable();
  EDIT_zmin->disable();

  EDIT_xmax->disable();
  EDIT_ymax->disable();
  EDIT_zmax->disable();
  OBJECT_CB(BLOCKAGE_AS_INPUT);

  EDIT_xmin->set_float_limits(xplt_orig[0],xplt_orig[ibar],GLUI_LIMIT_CLAMP);
  EDIT_xmax->set_float_limits(xplt_orig[0],xplt_orig[ibar],GLUI_LIMIT_CLAMP);
  EDIT_ymin->set_float_limits(yplt_orig[0],yplt_orig[jbar],GLUI_LIMIT_CLAMP);
  EDIT_ymax->set_float_limits(yplt_orig[0],yplt_orig[jbar],GLUI_LIMIT_CLAMP);
  EDIT_zmin->set_float_limits(zplt_orig[0],zplt_orig[kbar],GLUI_LIMIT_CLAMP);
  EDIT_zmax->set_float_limits(zplt_orig[0],zplt_orig[kbar],GLUI_LIMIT_CLAMP);

  ROLLOUT_geometry = glui_geometry->add_rollout("Unstructured",false);
  SPINNER_face_factor=glui_geometry->add_spinner_to_panel(ROLLOUT_geometry,"face factor",GLUI_SPINNER_FLOAT,&face_factor);
  SPINNER_face_factor->set_float_limits(0.0,0.5);

  glui_geometry->add_separator();
  glui_geometry->add_button(_("Save settings"),SAVE_SETTINGS,Blockedit_DLG_CB);
  BUTTON_blockage_1=glui_geometry->add_button(_("Close"),CLOSE_WINDOW,Blockedit_DLG_CB);


  glui_geometry->set_main_gfx_window( main_window );
}

/* ------------------ hide_glui_geometry ------------------------ */

extern "C" void hide_glui_geometry(void){
  blockageSelect=0;
  if(glui_geometry!=NULL)glui_geometry->hide();
  showedit_dialog_save=showedit_dialog;
  showedit_dialog=0;
  updatemenu=1;
  editwindow_status=CLOSE_WINDOW;
}

/* ------------------ show_glui_geometry ------------------------ */

extern "C" void show_glui_geometry(void){
  showedit_dialog=1;
  blockageSelect=1;
  update_blockvals(NOT_SELECT_BLOCKS);
  glui_geometry->show();
}

/* ------------------ Blockedit_DLG_CB ------------------------ */

void Blockedit_DLG_CB(int var){
  switch (var){
  case SAVE_SETTINGS:
    updatemenu=1;
    writeini(LOCAL_INI,NULL);
    break;
  case CLOSE_WINDOW: 
    DialogMenu(DIALOG_VIEW_GEOM);
    smooth_blockages();
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

  EDIT_xmin->set_float_limits(xplt_orig[0],xplt_orig[ibar],GLUI_LIMIT_CLAMP);
  EDIT_xmax->set_float_limits(xplt_orig[0],xplt_orig[ibar],GLUI_LIMIT_CLAMP);
  EDIT_ymin->set_float_limits(yplt_orig[0],yplt_orig[jbar],GLUI_LIMIT_CLAMP);
  EDIT_ymax->set_float_limits(yplt_orig[0],yplt_orig[jbar],GLUI_LIMIT_CLAMP);
  EDIT_zmin->set_float_limits(zplt_orig[0],zplt_orig[kbar],GLUI_LIMIT_CLAMP);
  EDIT_zmax->set_float_limits(zplt_orig[0],zplt_orig[kbar],GLUI_LIMIT_CLAMP);

  EDIT_xmin->set_float_val(xmin);
  EDIT_xmax->set_float_val(xmax);
  EDIT_ymin->set_float_val(ymin);
  EDIT_ymax->set_float_val(ymax);
  EDIT_zmin->set_float_val(zmin);
  EDIT_zmax->set_float_val(zmax);
  if(bchighlight!=NULL&&nsurfinfo>0){
    wall_case=bchighlight->walltype;
    OBJECT_CB(RADIO_WALL);
  }

  if(flag==SELECT_BLOCKS){
    if(bchighlight!=NULL){
      char dialog_label[255];
      mesh *blockmesh;

      if(nmeshes>1){
        blockmesh = meshinfo + bchighlight->meshindex;
        sprintf(dialog_label,"Mesh label: %s",blockmesh->label);
        STATIC_mesh_index->set_text(dialog_label);
      }
      sprintf(dialog_label,"&OBST index: %i",bchighlight->blockage_id);
      STATIC_blockage_index->set_text(dialog_label);
      strcpy(dialog_label,"&OBST label: ");
      strcat(dialog_label,bchighlight->label);
      STATIC_label->set_text(dialog_label);

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

      if(nsurfinfo>0){
        for(i=0;i<6;i++){
          surface_indices[i] = inv_sorted_surfidlist[bchighlight->surf_index[i]];
          surface_indices_bak[i] = inv_sorted_surfidlist[bchighlight->surf_index[i]];
          LIST_surface[i]->set_int_val(surface_indices[i]);
        }
      }
    }
    else{
      if(nsurfinfo>0){
        for(i=0;i<6;i++){
          surface_indices[i]=inv_sorted_surfidlist[0];
          surface_indices_bak[i]=inv_sorted_surfidlist[0];
          LIST_surface[i]->set_int_val(surface_indices[i]);
        }
      }
    }
  }
}

/* ------------------ OBJECT_CB ------------------------ */

extern "C" void OBJECT_CB(int var){
  int i,temp;
  switch (var){
    case UPDATE_LIST:
      if(surface_indices[NOT_USED]!=surface_indices_bak[NOT_USED]){
        surface_indices[NOT_USED]=surface_indices_bak[NOT_USED];
        LIST_surface[NOT_USED]->set_int_val(surface_indices[NOT_USED]);
      }
      switch (wall_case){
      case WALL_1:
        temp=surface_indices_bak[UP_Z];
        if(nsurfinfo>0){
          for(i=0;i<6;i++){
            surface_indices[i]=temp;
            LIST_surface[i]->set_int_val(temp);
          }
        }
        break;
      case WALL_3:
        if(nsurfinfo>0){
          for(i=0;i<6;i++){
            temp=surface_indices_bak[i];
            surface_indices[i]=temp;
            LIST_surface[i]->set_int_val(temp);
          }
        }
        break;
      case WALL_6:
        if(nsurfinfo>0){
          for(i=0;i<6;i++){
            temp=surface_indices_bak[i];
            surface_indices[i]=temp;
            LIST_surface[i]->set_int_val(temp);
          }
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
      }

      if(bchighlight!=NULL){
        for(i=0;i<6;i++){
          bchighlight->surf[i]=surfinfo+sorted_surfidlist[surface_indices_bak[i]];
          bchighlight->surf_index[i]=sorted_surfidlist[surface_indices_bak[i]];
        }
        bchighlight->changed_surface=1;
        if(bchighlight->blockage_id>0&&bchighlight->blockage_id<=nchanged_idlist){
          changed_idlist[bchighlight->blockage_id]=1;
        }
        blockages_dirty=1;
        updateusetextures();
        update_faces();
      }
      break;
    case RADIO_WALL:
      if(nsurfinfo==0)break;
      if(bchighlight!=NULL){
        bchighlight->walltype=wall_case;
      }
      switch (wall_case){
      case WALL_6:
        for(i=0;i<6;i++){
          LIST_surface[i]->enable();
        }
        LIST_surface[DOWN_Z]->set_name("z lower face");
        LIST_surface[UP_Z]->set_name("z upper face");
        LIST_surface[DOWN_Y]->set_name("y lower face");
        LIST_surface[UP_Y]->set_name("y upper face");
        LIST_surface[DOWN_X]->set_name("x lower face");
        LIST_surface[UP_X]->set_name("x upper face");
        break;
      case WALL_3:
        for(i=0;i<6;i++){
          LIST_surface[i]->disable();
        }
        LIST_surface[DOWN_Z]->enable();
        LIST_surface[UP_Z]->enable();
        LIST_surface[UP_Y]->enable();

        LIST_surface[DOWN_Z]->set_name("z lower face");
        LIST_surface[UP_Z]->set_name("z upper face");
        LIST_surface[UP_Y]->set_name("side faces");
        LIST_surface[DOWN_Y]->set_name("");
        LIST_surface[DOWN_X]->set_name("");
        LIST_surface[UP_X]->set_name("");

        break;
      case WALL_1:
        for(i=0;i<6;i++){
          LIST_surface[i]->disable();
        }
        LIST_surface[UP_Z]->enable();
        LIST_surface[UP_Z]->set_name("All faces");

        LIST_surface[DOWN_Z]->set_name("");
        LIST_surface[DOWN_Y]->set_name("");
        LIST_surface[UP_Y]->set_name("");
        LIST_surface[DOWN_X]->set_name("");
        LIST_surface[UP_X]->set_name("");
        break;
      default:
        ASSERT(FFALSE);
        break;
      }
      OBJECT_CB(UPDATE_LIST);
      break;
      case BLOCKAGE_AS_INPUT2:
      case BLOCKAGE_AS_INPUT:
        if(var==BLOCKAGE_AS_INPUT2){
          blockage_snapped=1-blockage_as_input;
          CHECKBOX_blockage->set_int_val(blockage_snapped);
        }
        blockage_as_input=1-blockage_snapped;
        if(blockage_as_input==1){
          blocklocation=BLOCKlocation_exact;
        }
        else{
          blocklocation=BLOCKlocation_grid;
        }
        update_blockvals(NOT_SELECT_BLOCKS);
        break;
    default:
      ASSERT(FFALSE);
      break;
  }
}
