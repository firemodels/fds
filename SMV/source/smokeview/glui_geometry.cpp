#define CPP
#include "options.h"

extern "C" void Volume_CB(int var);

#include <stdio.h>
#include <string.h>
#include GLUT_H

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
#define VISAXISLABELS 34
#define SHOW_TETRA 35

GLUI_Panel *PANEL_geom_surface=NULL;
GLUI_Panel *PANEL_geom_interior=NULL;
GLUI_Checkbox *CHECKBOX_surface_solid=NULL, *CHECKBOX_surface_outline;
GLUI_Checkbox *CHECKBOX_interior_solid=NULL, *CHECKBOX_interior_outline;

GLUI_Rollout *ROLLOUT_geomtest=NULL;
GLUI_Panel *PANEL_geom1=NULL;
GLUI_Panel *PANEL_geom1a=NULL;
GLUI_Panel *PANEL_geom1b=NULL;
GLUI_Panel *PANEL_geom1c=NULL;
GLUI_Panel *PANEL_geom1d=NULL;
GLUI_Panel *PANEL_geom2=NULL;
GLUI_Panel *PANEL_geom2a=NULL;
GLUI_Panel *PANEL_geom2b=NULL;
GLUI_Panel *PANEL_geom2c=NULL;
GLUI_Panel *PANEL_geom3a=NULL;
GLUI_Panel *PANEL_geom3b=NULL;
GLUI_Panel *PANEL_geom3c=NULL;
GLUI_Panel *PANEL_geom3ab=NULL;
GLUI_Panel *PANEL_geom3abc=NULL;
GLUI_Spinner *SPINNER_box_bounds[6];
GLUI_Spinner *SPINNER_box_translate[3];
GLUI_Spinner *SPINNER_tetra_vertices[12];
GLUI_Checkbox *CHECKBOX_tetrabox_showhide[10];
GLUI_Checkbox *CHECKBOX_visaxislabels;

#define VOL_BOXTRANSLATE 0
#define VOL_TETRA 1
#define UPDATE_VOLBOX_CONTROLS 2
#define VOL_SHOWHIDE 3


GLUI *glui_geometry=NULL;

GLUI_Button *BUTTON_blockage_1=NULL;

GLUI_Checkbox *CHECKBOX_blockage=NULL;

GLUI_EditText *EDIT_xmin=NULL, *EDIT_ymin=NULL, *EDIT_zmin=NULL;
GLUI_EditText *EDIT_xmax=NULL, *EDIT_ymax=NULL, *EDIT_zmax=NULL;

GLUI_Listbox *LIST_surface[7]={NULL,NULL,NULL,NULL,NULL,NULL,NULL};

GLUI_Panel *PANEL_obj_select=NULL,*PANEL_surface=NULL,*PANEL_interior=NULL,*PANEL_geom_showhide;
GLUI_Panel *PANEL_obj_stretch2=NULL,*PANEL_obj_stretch3=NULL, *PANEL_obj_stretch4=NULL;

GLUI_Rollout *ROLLOUT_structured=NULL;
GLUI_Rollout *ROLLOUT_unstructured=NULL;

GLUI_Spinner *SPINNER_face_factor=NULL;
GLUI_Spinner *SPINNER_tetra_line_thickness=NULL;
GLUI_Spinner *SPINNER_tetra_point_size = NULL;

GLUI_StaticText *STATIC_blockage_index=NULL;
GLUI_StaticText *STATIC_mesh_index=NULL;
GLUI_StaticText *STATIC_label=NULL;

void Blockedit_DLG_CB(int var);

char a_updatelabel[1000];
char *updatelabel=NULL;

/* ------------------ update_axislabels ------------------------ */

extern "C" void update_visaxislabels(void){
  if(CHECKBOX_visaxislabels!=NULL)CHECKBOX_visaxislabels->set_int_val(visaxislabels);
}

/* ------------------ update_geometry_controls ------------------------ */

extern "C" void update_geometry_controls(void){
  if(CHECKBOX_surface_solid!=NULL)CHECKBOX_surface_solid->set_int_val(showtrisurface);
  if(CHECKBOX_surface_outline!=NULL)CHECKBOX_surface_outline->set_int_val(showtrioutline);
  if(CHECKBOX_interior_solid!=NULL)CHECKBOX_surface_solid->set_int_val(show_geometry_interior_solid);
  if(CHECKBOX_interior_outline!=NULL)CHECKBOX_surface_outline->set_int_val(show_geometry_interior_outline);
}

/* ------------------ get_geom_dialog_state ------------------------ */

extern "C" void get_geom_dialog_state(void){
  if(ROLLOUT_structured!=NULL){
    if(ROLLOUT_structured->is_open){
      structured_isopen=1;
    }
    else{
      structured_isopen=0;
    }
  }
  if(ROLLOUT_unstructured!=NULL){
    if(ROLLOUT_unstructured->is_open){
      unstructured_isopen=1;
    }
    else{
      unstructured_isopen=0;
    }
  }
}

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

  ROLLOUT_structured = glui_geometry->add_rollout("Structured",false);
  if(structured_isopen==1)ROLLOUT_structured->open();
  PANEL_obj_select = glui_geometry->add_panel_to_panel(ROLLOUT_structured,"SURFs");

  PANEL_surface=glui_geometry->add_panel_to_panel(PANEL_obj_select,"",GLUI_PANEL_NONE);

  glui_geometry->add_column_to_panel(PANEL_surface,false);

  if(nsurfinfo>0){
    glui_geometry->add_statictext_to_panel(PANEL_surface,"");

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

  glui_geometry->add_column_to_panel(ROLLOUT_structured,false);

  PANEL_obj_stretch4=glui_geometry->add_panel_to_panel(ROLLOUT_structured,"",GLUI_PANEL_NONE);

  {
    char meshlabel[255];

    strcpy(meshlabel,_("Mesh:"));
    strcat(meshlabel,meshinfo->label);
    STATIC_mesh_index=glui_geometry->add_statictext_to_panel(PANEL_obj_stretch4,meshlabel);

  }
  STATIC_blockage_index=glui_geometry->add_statictext_to_panel(PANEL_obj_stretch4,"&OBST number: ");
  STATIC_label=glui_geometry->add_statictext_to_panel(PANEL_obj_stretch4,"&OBST label:");

  PANEL_obj_stretch2 = glui_geometry->add_panel_to_panel(ROLLOUT_structured,"Coordinates");

  CHECKBOX_blockage=glui_geometry->add_checkbox_to_panel(PANEL_obj_stretch2,_("Dimensions snapped to grid"),&blockage_snapped,
    BLOCKAGE_AS_INPUT,OBJECT_CB);
  CHECKBOX_visaxislabels=glui_geometry->add_checkbox_to_panel(PANEL_obj_stretch2,_("Show axis labels"),&visaxislabels,VISAXISLABELS,OBJECT_CB);
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

#ifdef pp_GEOMTEST
  ROLLOUT_unstructured = glui_geometry->add_rollout("Unstructured",false);
  if(unstructured_isopen==1)ROLLOUT_unstructured->open();
  SPINNER_face_factor=glui_geometry->add_spinner_to_panel(ROLLOUT_unstructured,"face factor",GLUI_SPINNER_FLOAT,&face_factor);
  SPINNER_face_factor->set_float_limits(0.0,0.5);

  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi = meshinfo + i;
    if(meshi->ncutcells>0){
      glui_geometry->add_checkbox_to_panel(ROLLOUT_unstructured,_("Show cutcells"),&show_cutcells);
      break;
    }
  }
  PANEL_geom_showhide = glui_geometry->add_panel_to_panel(ROLLOUT_unstructured,"",GLUI_PANEL_NONE);
  PANEL_surface = glui_geometry->add_panel_to_panel(PANEL_geom_showhide,"surface");
  CHECKBOX_surface_solid=glui_geometry->add_checkbox_to_panel(PANEL_surface,"solid",&showtrisurface,VOL_SHOWHIDE,Volume_CB);
  CHECKBOX_surface_outline=glui_geometry->add_checkbox_to_panel(PANEL_surface,"outline",&showtrioutline,VOL_SHOWHIDE,Volume_CB);

  glui_geometry->add_column_to_panel(PANEL_geom_showhide,false);
  PANEL_interior = glui_geometry->add_panel_to_panel(PANEL_geom_showhide,"interior");
  CHECKBOX_interior_solid=glui_geometry->add_checkbox_to_panel(PANEL_interior,"solid",&show_geometry_interior_solid,VOL_SHOWHIDE,Volume_CB);
  CHECKBOX_interior_outline=glui_geometry->add_checkbox_to_panel(PANEL_interior,"outline",&show_geometry_interior_outline,VOL_SHOWHIDE,Volume_CB);

  // -------------- Cube/Tetra intersection test -------------------

  ROLLOUT_geomtest = glui_geometry->add_rollout_to_panel(ROLLOUT_unstructured,"Cube/Tetra intersection test",false);
  glui_geometry->add_checkbox_to_panel(ROLLOUT_geomtest, "show intersection region", &show_geomtest, SHOW_TETRA, Volume_CB);
  glui_geometry->add_checkbox_to_panel(ROLLOUT_geomtest,"show area labels", &show_tetratest_labels);
  SPINNER_tetra_line_thickness=glui_geometry->add_spinner_to_panel(ROLLOUT_geomtest,"line thickness",GLUI_SPINNER_FLOAT,&tetra_line_thickness);
//  SPINNER_tetra_line_thickness->set_float_limits(1.0, 10.0);
  SPINNER_tetra_point_size = glui_geometry->add_spinner_to_panel(ROLLOUT_geomtest, "point size", GLUI_SPINNER_FLOAT, &tetra_point_size);
//  SPINNER_tetra_point_size->set_float_limits(1.0, 20.0);

  PANEL_geom1 = glui_geometry->add_panel_to_panel(ROLLOUT_geomtest,"box bounding planes");

  PANEL_geom1d=glui_geometry->add_panel_to_panel(PANEL_geom1,"",GLUI_PANEL_NONE);
  PANEL_geom1a=glui_geometry->add_panel_to_panel(PANEL_geom1d,"",GLUI_PANEL_NONE);
  glui_geometry->add_column_to_panel(PANEL_geom1d,false);
  PANEL_geom1b=glui_geometry->add_panel_to_panel(PANEL_geom1d,"",GLUI_PANEL_NONE);

  PANEL_geom1c=glui_geometry->add_panel_to_panel(PANEL_geom1,"",GLUI_PANEL_NONE);

  SPINNER_box_bounds[0]=glui_geometry->add_spinner_to_panel(PANEL_geom1a,"xmin",GLUI_SPINNER_FLOAT,box_bounds2,VOL_BOXTRANSLATE,Volume_CB);
  SPINNER_box_bounds[2]=glui_geometry->add_spinner_to_panel(PANEL_geom1a,"ymin",GLUI_SPINNER_FLOAT,box_bounds2+2,VOL_BOXTRANSLATE,Volume_CB);
  SPINNER_box_bounds[4]=glui_geometry->add_spinner_to_panel(PANEL_geom1a,"zmin",GLUI_SPINNER_FLOAT,box_bounds2+4,VOL_BOXTRANSLATE,Volume_CB);
  SPINNER_box_bounds[1]=glui_geometry->add_spinner_to_panel(PANEL_geom1b,"xmax",GLUI_SPINNER_FLOAT,box_bounds2+1,VOL_BOXTRANSLATE,Volume_CB);
  SPINNER_box_bounds[3]=glui_geometry->add_spinner_to_panel(PANEL_geom1b,"ymax",GLUI_SPINNER_FLOAT,box_bounds2+3,VOL_BOXTRANSLATE,Volume_CB);
  SPINNER_box_bounds[5]=glui_geometry->add_spinner_to_panel(PANEL_geom1b,"zmax",GLUI_SPINNER_FLOAT,box_bounds2+5,VOL_BOXTRANSLATE,Volume_CB);

  SPINNER_box_translate[0]=glui_geometry->add_spinner_to_panel(PANEL_geom1c,"translate: x",GLUI_SPINNER_FLOAT,box_translate,VOL_BOXTRANSLATE,Volume_CB);
  glui_geometry->add_column_to_panel(PANEL_geom1c,false);
  SPINNER_box_translate[1]=glui_geometry->add_spinner_to_panel(PANEL_geom1c,"y",GLUI_SPINNER_FLOAT,box_translate+1,VOL_BOXTRANSLATE,Volume_CB);
  glui_geometry->add_column_to_panel(PANEL_geom1c,false);
  SPINNER_box_translate[2]=glui_geometry->add_spinner_to_panel(PANEL_geom1c,"z",GLUI_SPINNER_FLOAT,box_translate+2,VOL_BOXTRANSLATE,Volume_CB);
  Volume_CB(VOL_BOXTRANSLATE);
  PANEL_geom2=glui_geometry->add_panel_to_panel(ROLLOUT_geomtest,"tetrahedron vertices");
  PANEL_geom2a=glui_geometry->add_panel_to_panel(PANEL_geom2,"",GLUI_PANEL_NONE);
  glui_geometry->add_column_to_panel(PANEL_geom2,false);
  PANEL_geom2b=glui_geometry->add_panel_to_panel(PANEL_geom2,"",GLUI_PANEL_NONE);
  glui_geometry->add_column_to_panel(PANEL_geom2,false);
  PANEL_geom2c=glui_geometry->add_panel_to_panel(PANEL_geom2,"",GLUI_PANEL_NONE);

  SPINNER_tetra_vertices[0]=glui_geometry->add_spinner_to_panel(PANEL_geom2a,"v1 x:",GLUI_SPINNER_FLOAT,tetra_vertices,VOL_TETRA,Volume_CB);
  SPINNER_tetra_vertices[3]=glui_geometry->add_spinner_to_panel(PANEL_geom2a,"v2 x:",GLUI_SPINNER_FLOAT,tetra_vertices+3,VOL_TETRA,Volume_CB);
  SPINNER_tetra_vertices[6]=glui_geometry->add_spinner_to_panel(PANEL_geom2a,"v3 x:",GLUI_SPINNER_FLOAT,tetra_vertices+6,VOL_TETRA,Volume_CB);
  SPINNER_tetra_vertices[9]=glui_geometry->add_spinner_to_panel(PANEL_geom2a,"v4 x:",GLUI_SPINNER_FLOAT,tetra_vertices+9,VOL_TETRA,Volume_CB);

  SPINNER_tetra_vertices[1]=glui_geometry->add_spinner_to_panel(PANEL_geom2b,"y:",GLUI_SPINNER_FLOAT,tetra_vertices+1,VOL_TETRA,Volume_CB);
  SPINNER_tetra_vertices[4]=glui_geometry->add_spinner_to_panel(PANEL_geom2b,"y:",GLUI_SPINNER_FLOAT,tetra_vertices+4,VOL_TETRA,Volume_CB);
  SPINNER_tetra_vertices[7]=glui_geometry->add_spinner_to_panel(PANEL_geom2b,"y:",GLUI_SPINNER_FLOAT,tetra_vertices+7,VOL_TETRA,Volume_CB);
  SPINNER_tetra_vertices[10]=glui_geometry->add_spinner_to_panel(PANEL_geom2b,"y:",GLUI_SPINNER_FLOAT,tetra_vertices+10,VOL_TETRA,Volume_CB);

  SPINNER_tetra_vertices[2]=glui_geometry->add_spinner_to_panel(PANEL_geom2c,"z:",GLUI_SPINNER_FLOAT,tetra_vertices+2,VOL_TETRA,Volume_CB);
  SPINNER_tetra_vertices[5]=glui_geometry->add_spinner_to_panel(PANEL_geom2c,"z:",GLUI_SPINNER_FLOAT,tetra_vertices+5,VOL_TETRA,Volume_CB);
  SPINNER_tetra_vertices[8]=glui_geometry->add_spinner_to_panel(PANEL_geom2c,"z:",GLUI_SPINNER_FLOAT,tetra_vertices+8,VOL_TETRA,Volume_CB);
  SPINNER_tetra_vertices[11]=glui_geometry->add_spinner_to_panel(PANEL_geom2c,"z:",GLUI_SPINNER_FLOAT,tetra_vertices+11,VOL_TETRA,Volume_CB);

  PANEL_geom3abc=glui_geometry->add_panel_to_panel(ROLLOUT_geomtest,"box/tetrahedron faces",GLUI_PANEL_NONE);

  glui_geometry->add_checkbox_to_panel(ROLLOUT_geomtest,"tetra test",&show_test_in_tetra);
  glui_geometry->add_spinner_to_panel(ROLLOUT_geomtest,"tetra x:",GLUI_SPINNER_FLOAT,tetra_xyz);
  glui_geometry->add_spinner_to_panel(ROLLOUT_geomtest,"tetra y:",GLUI_SPINNER_FLOAT,tetra_xyz+1);
  glui_geometry->add_spinner_to_panel(ROLLOUT_geomtest,"tetra z:",GLUI_SPINNER_FLOAT,tetra_xyz+2);

  PANEL_geom3ab=glui_geometry->add_panel_to_panel(PANEL_geom3abc,"box");
  PANEL_geom3a=glui_geometry->add_panel_to_panel(PANEL_geom3ab,"",GLUI_PANEL_NONE);
  glui_geometry->add_column_to_panel(PANEL_geom3ab,false);
  PANEL_geom3b=glui_geometry->add_panel_to_panel(PANEL_geom3ab,"",GLUI_PANEL_NONE);


  glui_geometry->add_column_to_panel(PANEL_geom3abc,false);
  PANEL_geom3c=glui_geometry->add_panel_to_panel(PANEL_geom3abc,"tetrahedron");

  CHECKBOX_tetrabox_showhide[0]=glui_geometry->add_checkbox_to_panel(PANEL_geom3a,"xmin",tetrabox_vis+0);
  CHECKBOX_tetrabox_showhide[1]=glui_geometry->add_checkbox_to_panel(PANEL_geom3b,"xmax",tetrabox_vis+1);
  CHECKBOX_tetrabox_showhide[2]=glui_geometry->add_checkbox_to_panel(PANEL_geom3a,"ymin",tetrabox_vis+2);
  CHECKBOX_tetrabox_showhide[3]=glui_geometry->add_checkbox_to_panel(PANEL_geom3b,"ymax",tetrabox_vis+3);
  CHECKBOX_tetrabox_showhide[4]=glui_geometry->add_checkbox_to_panel(PANEL_geom3a,"zmin",tetrabox_vis+4);
  CHECKBOX_tetrabox_showhide[5]=glui_geometry->add_checkbox_to_panel(PANEL_geom3b,"zmax",tetrabox_vis+5);
  CHECKBOX_tetrabox_showhide[8]=glui_geometry->add_checkbox_to_panel(PANEL_geom3c,"v1 v3 v4",tetrabox_vis+6);
  CHECKBOX_tetrabox_showhide[7]=glui_geometry->add_checkbox_to_panel(PANEL_geom3c,"v2 v3 v4",tetrabox_vis+7);
  CHECKBOX_tetrabox_showhide[6]=glui_geometry->add_checkbox_to_panel(PANEL_geom3c,"v1 v2 v4",tetrabox_vis+8);
  CHECKBOX_tetrabox_showhide[9]=glui_geometry->add_checkbox_to_panel(PANEL_geom3c,"v1 v2 v3",tetrabox_vis+9);
#endif

  glui_geometry->add_separator();
  glui_geometry->add_button(_("Save settings"),SAVE_SETTINGS,Blockedit_DLG_CB);
  BUTTON_blockage_1=glui_geometry->add_button(_("Close"),CLOSE_WINDOW,Blockedit_DLG_CB);


  glui_geometry->set_main_gfx_window( main_window );
}

/* ------------------ Volume_CB ------------------------ */

extern "C" void Volume_CB(int var){
  int i;
  switch(var){
    case SHOW_TETRA:
      if(show_geomtest==1){
        BlockageMenu(visBLOCKHide);
        VentMenu(HIDE_ALL_VENTS);
      }
      break;
    case VOL_BOXTRANSLATE:
      box_bounds[0]=box_bounds2[0]+box_translate[0];
      box_bounds[1]=box_bounds2[1]+box_translate[0];
      box_bounds[2]=box_bounds2[2]+box_translate[1];
      box_bounds[3]=box_bounds2[3]+box_translate[1];
      box_bounds[4]=box_bounds2[4]+box_translate[2];
      box_bounds[5]=box_bounds2[5]+box_translate[2];
      update_volbox_controls=1;
      break;
    case VOL_TETRA:
      update_volbox_controls=1;
      break;
    case UPDATE_VOLBOX_CONTROLS:
      update_volbox_controls=0;
      for(i=0;i<10;i++){
        if(face_vis[i]!=face_vis_old[i]){
          //if(face_vis[i]==1){
          //  CHECKBOX_tetrabox_showhide[i]->enable();
          // }
          // else{
          //   CHECKBOX_tetrabox_showhide[i]->disable();
          // }
          face_vis_old[i]=face_vis[i];
        }
      }
      break;
    case VOL_SHOWHIDE:
      updatemenu=1;
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
}

/* ------------------ hide_glui_geometry ------------------------ */

extern "C" void hide_glui_geometry(void){
  blockageSelect=0;
  if(glui_geometry!=NULL)glui_geometry->hide();
  showedit_dialog=0;
  updatemenu=1;
  editwindow_status=CLOSE_WINDOW;
}

/* ------------------ show_glui_geometry ------------------------ */

extern "C" void show_glui_geometry(void){
  showedit_dialog=1;
  blockageSelect=1;
  Update_Blockvals(NOT_SELECT_BLOCKS);
  if(glui_geometry!=NULL)glui_geometry->show();
}

/* ------------------ Blockedit_DLG_CB ------------------------ */

void Blockedit_DLG_CB(int var){
  switch(var){
  case SAVE_SETTINGS:
    updatemenu=1;
    writeini(LOCAL_INI,NULL);
    break;
  case CLOSE_WINDOW: 
    DialogMenu(DIALOG_GEOMETRY);
    smooth_blockages();
    break;
  default:
    ASSERT(FFALSE);
    break;
  }

}

/* ------------------ Update_Blockvals ------------------------ */

extern "C" void Update_Blockvals(int flag){
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

      switch(wall_case){
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
  switch(var){
    case VISAXISLABELS:
      updatemenu=1;
      break;
    case UPDATE_LIST:
      switch(wall_case){
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
        UpdateFaces();
      }
      break;
    case RADIO_WALL:
      if(nsurfinfo==0)break;
      if(bchighlight!=NULL){
        bchighlight->walltype=wall_case;
      }
      switch(wall_case){
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
        if(blocklocation!=BLOCKlocation_cad){
          if(blockage_as_input==1){
            blocklocation=BLOCKlocation_exact;
          }
          else{
            blocklocation=BLOCKlocation_grid;
          }
        }
        Update_Blockvals(NOT_SELECT_BLOCKS);
        break;
    default:
      ASSERT(FFALSE);
      break;
  }
}
