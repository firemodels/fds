#include "options.h"
#include <string.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "glui.h"
#include "flowfiles.h"
#define CPP
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

// svn revision character string
extern "C" char glui_blockedit_revision[]="$Revision$";

#define NEW_OBJECT 1
#define DELETE_OBJECT 4
#define UNDELETE_OBJECT 5
#define RADIO_XYZ 6
#define RADIO_OBJECT_TYPE 7
#define XMIN_SPIN 20
#define YMIN_SPIN 21
#define ZMIN_SPIN 22
#define XMAX_SPIN 23
#define YMAX_SPIN 24
#define ZMAX_SPIN 25
#define UPDATE_LIST 31
#define RADIO_WALL 32
#define SET_LABEL 33
#define MESH_LIST 34

GLUI *glui_edit=NULL;
GLUI_Panel *panel_label=NULL,*panel_obj_select=NULL,*panel_surface=NULL;
GLUI_Panel *panel_obj_new=NULL;
GLUI_RadioGroup *group2=NULL;
GLUI_RadioGroup *wallgroup=NULL;
GLUI_Panel *panel_obj_stretch=NULL, *panel_obj_stretch2=NULL,*panel_obj_stretch3=NULL, *panel_obj_stretch4=NULL;
GLUI_EditText *edittext_xmin=NULL, *edittext_ymin=NULL, *edittext_zmin=NULL;
GLUI_EditText *edittext_xmax=NULL, *edittext_ymax=NULL, *edittext_zmax=NULL;
GLUI_EditText *edittextlabel=NULL;
GLUI_Listbox *surfacelists[6]={NULL,NULL,NULL,NULL,NULL,NULL};
GLUI_Listbox *meshlist=NULL;
GLUI_Spinner *SPINNER_move=NULL, *SPINNER_stretch_white=NULL,*SPINNER_stretch_black;
GLUI_Button *Update_Button=NULL;
GLUI_Button *Cancel_Button=NULL;
GLUI_Button *Undo_Button=NULL;
GLUI_Button *UndoAll_Button=NULL;
GLUI_Checkbox *blockage_checkbox=NULL;

extern "C" void OBJECT_CB(int var);

void OBJECT_MOVE_STRETCH_CB(int var);
void BUTTON_hide3_CB(int var);

char a_updatelabel[1000];
char *updatelabel=NULL;


extern "C" void update_highlight_mesh(void){
  if(meshlist!=NULL)meshlist->set_int_val(highlight_mesh);
}

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
  glui_edit = GLUI_Master.create_glui("Edit",0,0,0);
  if(showedit==0)glui_edit->hide();

  panel_obj_select = glui_edit->add_panel("Select");


  panel_obj_new=glui_edit->add_panel_to_panel(panel_obj_select,"",GLUI_PANEL_NONE);
  if(nmeshes>1){
    meshlist = glui_edit->add_listbox_to_panel(panel_obj_new,"Mesh",&highlight_mesh,MESH_LIST,OBJECT_CB);
    for(i=0;i<selected_case->nmeshes;i++){
      meshi = selected_case->meshinfo + i;
      meshlist->add_item(i,meshi->label);
    }
  }

  glui_edit->add_button_to_panel(panel_obj_new,"New",NEW_OBJECT,OBJECT_CB);
  glui_edit->add_button_to_panel(panel_obj_new,"Delete",DELETE_OBJECT,OBJECT_CB);
  glui_edit->add_button_to_panel(panel_obj_new,"UnDelete",UNDELETE_OBJECT,OBJECT_CB);

  panel_label=glui_edit->add_panel_to_panel(panel_obj_select,"",GLUI_PANEL_NONE);

  edittextlabel = glui_edit->add_edittext_to_panel(panel_label,"Label",GLUI_EDITTEXT_TEXT,NULL,SET_LABEL,OBJECT_CB);
  edittextlabel->set_w(260);

  panel_surface=glui_edit->add_panel_to_panel(panel_obj_select,"",GLUI_PANEL_NONE);

  if(nsurfaces>0){
    wallgroup = glui_edit->add_radiogroup_to_panel(panel_surface,&wall_case,RADIO_WALL,OBJECT_CB);
    glui_edit->add_radiobutton_to_group(wallgroup,"1 Wall");
    glui_edit->add_radiobutton_to_group(wallgroup,"3 Wall");
    glui_edit->add_radiobutton_to_group(wallgroup,"6 Wall");
  }
  glui_edit->add_checkbox_to_panel(panel_surface,"Color",&visNormalEditColors);


  glui_edit->add_column_to_panel(panel_surface,false);

  if(nsurfaces>0){
    glui_edit->add_statictext_to_panel(panel_surface,"Surface Type:");
    surfacelists[DOWN_X] = glui_edit->add_listbox_to_panel(panel_surface,"Left",surface_indices+DOWN_X,UPDATE_LIST,OBJECT_CB);
    surfacelists[DOWN_X]->set_w(260);
    for(i=0;i<nsurfaces;i++){
      surfi = surfaceinfo + sorted_surfidlist[i];
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      surfacelists[DOWN_X]->add_item(i,surfacelabel);
    }

    surfacelists[UP_X] = glui_edit->add_listbox_to_panel(panel_surface,"Right",surface_indices+UP_X,UPDATE_LIST,OBJECT_CB);
    surfacelists[UP_X]->set_w(260);
    for(i=0;i<nsurfaces;i++){
      surfi = surfaceinfo + sorted_surfidlist[i];
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      surfacelists[UP_X]->add_item(i,surfacelabel);
    }

    surfacelists[DOWN_Y] = glui_edit->add_listbox_to_panel(panel_surface,"Front",surface_indices+DOWN_Y,UPDATE_LIST,OBJECT_CB);
    surfacelists[DOWN_Y]->set_w(260);
    for(i=0;i<nsurfaces;i++){
      surfi = surfaceinfo + sorted_surfidlist[i];
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      surfacelists[DOWN_Y]->add_item(i,surfacelabel);
    }

    surfacelists[UP_Y] = glui_edit->add_listbox_to_panel(panel_surface,"Back",surface_indices+UP_Y,UPDATE_LIST,OBJECT_CB);
    surfacelists[UP_Y]->set_w(260);
    for(i=0;i<nsurfaces;i++){
      surfi = surfaceinfo + sorted_surfidlist[i];
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      surfacelists[UP_Y]->add_item(i,surfacelabel);
    }

    surfacelists[DOWN_Z] = glui_edit->add_listbox_to_panel(panel_surface,"Down",surface_indices+DOWN_Z,UPDATE_LIST,OBJECT_CB);
    surfacelists[DOWN_Z]->set_w(260);
    for(i=0;i<nsurfaces;i++){
      surfi = surfaceinfo + sorted_surfidlist[i];
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      surfacelists[DOWN_Z]->add_item(i,surfacelabel);
    }

    surfacelists[UP_Z] = glui_edit->add_listbox_to_panel(panel_surface,"Up",surface_indices+UP_Z,UPDATE_LIST,OBJECT_CB);
    surfacelists[UP_Z]->set_w(260);
    for(i=0;i<nsurfaces;i++){
      surfi = surfaceinfo + sorted_surfidlist[i];
      if(surfi->obst_surface==0)continue;
      surfacelabel = surfi->surfacelabel;
      surfacelists[UP_Z]->add_item(i,surfacelabel);
    }

    OBJECT_CB(RADIO_WALL);
  }

  glui_edit->add_column(false);

  panel_obj_stretch2 = glui_edit->add_panel("Move/Stretch");

  panel_obj_stretch = glui_edit->add_panel_to_panel(panel_obj_stretch2,"",GLUI_PANEL_NONE);
  group2 = glui_edit->add_radiogroup_to_panel(panel_obj_stretch,&xyz_dir,RADIO_XYZ,OBJECT_CB);
  glui_edit->add_radiobutton_to_group(group2,"X");
  glui_edit->add_radiobutton_to_group(group2,"Y");
  glui_edit->add_radiobutton_to_group(group2,"Z");
  glui_edit->add_column_to_panel(panel_obj_stretch,false);
  SPINNER_stretch_black=glui_edit->add_spinner_to_panel(panel_obj_stretch,"Stretch black",GLUI_SPINNER_INT,&stretch_var_black,STRETCH_BLACK,OBJECT_MOVE_STRETCH_CB);
  SPINNER_stretch_black->set_int_limits(0,ibar,GLUI_LIMIT_CLAMP);
  SPINNER_stretch_white=glui_edit->add_spinner_to_panel(panel_obj_stretch,"Stretch white",GLUI_SPINNER_INT,&stretch_var_white,STRETCH_WHITE,OBJECT_MOVE_STRETCH_CB);
  SPINNER_stretch_white->set_int_limits(0,ibar,GLUI_LIMIT_CLAMP);
  SPINNER_move=glui_edit->add_spinner_to_panel(panel_obj_stretch,"Move",GLUI_SPINNER_INT,&move_var,MOVE,OBJECT_MOVE_STRETCH_CB);
  SPINNER_move->set_int_limits(0,ibar,GLUI_LIMIT_CLAMP);

  panel_obj_stretch3 = glui_edit->add_panel_to_panel(panel_obj_stretch2,"",GLUI_PANEL_NONE);
  edittext_xmin=glui_edit->add_edittext_to_panel(panel_obj_stretch3,"x",GLUI_EDITTEXT_FLOAT,&glui_block_xmin,XMIN_SPIN,OBJECT_CB);
  edittext_ymin=glui_edit->add_edittext_to_panel(panel_obj_stretch3,"y",GLUI_EDITTEXT_FLOAT,&glui_block_ymin,YMIN_SPIN,OBJECT_CB);
  edittext_zmin=glui_edit->add_edittext_to_panel(panel_obj_stretch3,"z",GLUI_EDITTEXT_FLOAT,&glui_block_zmin,ZMIN_SPIN,OBJECT_CB);

  glui_edit->add_column_to_panel(panel_obj_stretch3,false);
  edittext_xmax=glui_edit->add_edittext_to_panel(panel_obj_stretch3,"",GLUI_EDITTEXT_FLOAT,&glui_block_xmax,XMAX_SPIN,OBJECT_CB);
  edittext_ymax=glui_edit->add_edittext_to_panel(panel_obj_stretch3,"",GLUI_EDITTEXT_FLOAT,&glui_block_ymax,YMAX_SPIN,OBJECT_CB);
  edittext_zmax=glui_edit->add_edittext_to_panel(panel_obj_stretch3,"",GLUI_EDITTEXT_FLOAT,&glui_block_zmax,ZMAX_SPIN,OBJECT_CB);

  blockage_checkbox=glui_edit->add_checkbox_to_panel(panel_obj_stretch3,"Blockages as Input",&blockage_as_input,BLOCKAGE_AS_INPUT,OBJECT_CB);

  edittext_xmin->set_float_limits(xplt_orig[0],xplt_orig[ibar],GLUI_LIMIT_CLAMP);
  edittext_xmax->set_float_limits(xplt_orig[0],xplt_orig[ibar],GLUI_LIMIT_CLAMP);
  edittext_ymin->set_float_limits(yplt_orig[0],yplt_orig[jbar],GLUI_LIMIT_CLAMP);
  edittext_ymax->set_float_limits(yplt_orig[0],yplt_orig[jbar],GLUI_LIMIT_CLAMP);
  edittext_zmin->set_float_limits(zplt_orig[0],zplt_orig[kbar],GLUI_LIMIT_CLAMP);
  edittext_zmax->set_float_limits(zplt_orig[0],zplt_orig[kbar],GLUI_LIMIT_CLAMP);

  panel_obj_stretch4=glui_edit->add_panel("",GLUI_PANEL_NONE);
  Undo_Button=glui_edit->add_button_to_panel(panel_obj_stretch4,"Undo Selected Blockage",UNDO_BLOCKAGE,BUTTON_hide3_CB);
  UndoAll_Button=glui_edit->add_button_to_panel(panel_obj_stretch4,"Undo All Blockages",UNDO_ALL_BLOCKAGES,BUTTON_hide3_CB);

  glui_edit->add_separator_to_panel(panel_obj_stretch4);
  if(fds_fileout!=NULL){
    updatelabel=a_updatelabel;
    strcpy(updatelabel,"Update ");
    strcat(updatelabel,fds_fileout);
    Update_Button=glui_edit->add_button_to_panel(panel_obj_stretch4,updatelabel,UPDATE_WINDOW,BUTTON_hide3_CB);
  }
  if(fds_fileout==NULL){
    Update_Button=glui_edit->add_button_to_panel(panel_obj_stretch4,"(FDS Data File Unknown)",UPDATE_WINDOW,BUTTON_hide3_CB);
    Update_Button->disable();
  }

  glui_edit->add_button_to_panel(panel_obj_stretch4,"Close",CLOSE_WINDOW,BUTTON_hide3_CB);

  glui_edit->set_main_gfx_window( main_window );
}

/* ------------------ hide_glui_edit ------------------------ */

extern "C" void hide_glui_edit(void){
  blockageSelect=0;
  if(glui_edit!=NULL)glui_edit->hide();
  showedit=0;
  updatemenu=1;
  outputchangedblockages();
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
    if(blockages_dirty==1)outputchangedblockages();
    DialogMenu(16);
    smooth_blockages();
    break;
  case UPDATE_WINDOW:
    outputchangedblockages();
    blockages_dirty=0;
    break;
  case CANCEL_WINDOW:
    break;
  case UNDO_BLOCKAGE:
    if(bchighlight!=NULL){
      bchighlight->xmin=bchighlight->xyzORIG[0];
      bchighlight->xmax=bchighlight->xyzORIG[1];
      bchighlight->ymin=bchighlight->xyzORIG[2];
      bchighlight->ymax=bchighlight->xyzORIG[3];
      bchighlight->zmin=bchighlight->xyzORIG[4];
      bchighlight->zmax=bchighlight->xyzORIG[5];
      bchighlight->ijk[0]=bchighlight->ijkORIG[0];
      bchighlight->ijk[1]=bchighlight->ijkORIG[1];
      bchighlight->ijk[2]=bchighlight->ijkORIG[2];
      bchighlight->ijk[3]=bchighlight->ijkORIG[3];
      bchighlight->ijk[4]=bchighlight->ijkORIG[4];
      bchighlight->ijk[5]=bchighlight->ijkORIG[5];
      for(i=0;i<6;i++){
        bchighlight->surf[i]=bchighlight->surfORIG[i];
        bchighlight->surf_index[i]=bchighlight->surf_indexORIG[i];
      }
      for(i=0;i<6;i++){
        surface_indices[i]=inv_sorted_surfidlist[bchighlight->surf_index[i]];
      }
      bchighlight->walltype=bchighlight->walltypeORIG;
      updateusetextures();
      update_blockvals(1);
     // OBJECT_CB(RADIO_WALL);
      update_faces();
    }
    break;
  case UNDO_ALL_BLOCKAGES:
    bchighlight_save=bchighlight;
    for(i=0;i<selected_case->nmeshes;i++){
      meshi=selected_case->meshinfo+i;
      for(j=0;j<meshi->nbptrs;j++){
        bc=meshi->blockageinfoptrs[j];
        if(bc==NULL)continue;
        if(bc->changed==0&&bc->changed_surface==0)continue;
        bchighlight=bc;
        bc->xmin=bc->xyzORIG[0];
        bc->xmax=bc->xyzORIG[1];
        bc->ymin=bc->xyzORIG[2];
        bc->ymax=bc->xyzORIG[3];
        bc->zmin=bc->xyzORIG[4];
        bc->zmax=bc->xyzORIG[5];
        bc->ijk[0]=bc->ijkORIG[0];
        bc->ijk[1]=bc->ijkORIG[1];
        bc->ijk[2]=bc->ijkORIG[2];
        bc->ijk[3]=bc->ijkORIG[3];
        bc->ijk[4]=bc->ijkORIG[4];
        bc->ijk[5]=bc->ijkORIG[5];
        for(k=0;k<6;k++){
          bc->surf[k]=bc->surfORIG[k];
          bc->surf_index[k]=bc->surf_indexORIG[k];
        }
        for(i=0;i<6;i++){
          surface_indices[i]=inv_sorted_surfidlist[bc->surf_index[i]];
        }
        bc->walltype=bc->walltypeORIG;
        update_blockvals(1);
        OBJECT_CB(RADIO_WALL);
      }
    }
    bchighlight=bchighlight_save;
    updateusetextures();
    update_faces();
        
    break;
  default:
    ASSERT(FFALSE);
    break;
  }

}

/* ------------------ update_movestretch ------------------------ */

extern "C" void update_movestretch(blockagedata *bc){

  int ibar,jbar,kbar;
  int imin, imax, jmin, jmax, kmin, kmax;

  ibar=current_mesh->ibar;
  jbar=current_mesh->jbar;
  kbar=current_mesh->kbar;

  if(bc==NULL){
    imin=0;
    imax=1;
    jmin=0;
    jmax=1;
    kmin=0;
    kmax=1;
  }
  else{
    imin=bc->ijk[IMIN];
    imax=bc->ijk[IMAX];
    jmin=bc->ijk[JMIN];
    jmax=bc->ijk[JMAX];
    kmin=bc->ijk[KMIN];
    kmax=bc->ijk[KMAX];
  }

  switch (xyz_dir){
  case XDIR:
    SPINNER_move->set_int_limits(0,ibar,GLUI_LIMIT_CLAMP);
    SPINNER_move->set_int_val(imin);
    SPINNER_stretch_black->set_int_limits(0,ibar,GLUI_LIMIT_CLAMP);
    SPINNER_stretch_black->set_int_val(imin);
    SPINNER_stretch_white->set_int_limits(0,ibar,GLUI_LIMIT_CLAMP);
    SPINNER_stretch_white->set_int_val(imax);
    break;
  case YDIR:
    SPINNER_move->set_int_limits(0,jbar,GLUI_LIMIT_CLAMP);
    SPINNER_move->set_int_val(jmin);
    SPINNER_stretch_black->set_int_limits(0,jbar,GLUI_LIMIT_CLAMP);
    SPINNER_stretch_black->set_int_val(jmin);
    SPINNER_stretch_white->set_int_limits(0,jbar,GLUI_LIMIT_CLAMP);
    SPINNER_stretch_white->set_int_val(jmax);
    break;
  case ZDIR:
    SPINNER_move->set_int_limits(0,kbar,GLUI_LIMIT_CLAMP);
    SPINNER_move->set_int_val(kmin);
    SPINNER_stretch_black->set_int_limits(0,ibar,GLUI_LIMIT_CLAMP);
    SPINNER_stretch_black->set_int_val(kmin);
    SPINNER_stretch_white->set_int_limits(0,kbar,GLUI_LIMIT_CLAMP);
    SPINNER_stretch_white->set_int_val(kmax);
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
    wallgroup->set_int_val(wall_case);
    OBJECT_CB(RADIO_WALL);
  }

  if(flag==1){
    if(bchighlight_old!=NULL){
      label=edittextlabel->get_text();
      obstlabelcopy(&bchighlight_old->label,label);
    }
    if(bchighlight!=NULL){
      edittextlabel->set_text(bchighlight->label);

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
          surfacelists[i]->set_int_val(surface_indices[i]);
        }
      }
    }
    else{
      edittextlabel->set_text("");
      if(nsurfaces>0){
        for(i=0;i<6;i++){
          surface_indices[i]=inv_sorted_surfidlist[0];
          surfacelists[i]->set_int_val(surface_indices[i]);
        }
      }
    }
  }
  update_movestretch(bchighlight);
}

/* ------------------ OBJECT_MOVE_STRETCH_CB ------------------------ */

void OBJECT_MOVE_STRETCH_CB(int var){
  switch (var){
  case MOVE:
    which_face=0;
    moveiblockage(move_var);
    break;
  case STRETCH_BLACK:
    which_face=-1;
    stretchiblockage(stretch_var_black);
    break;
  case STRETCH_WHITE:
    which_face=1;
    stretchiblockage(stretch_var_white);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  update_blockvals(0);
  update_faces();
}

/* ------------------ update_xyzdir ------------------------ */

extern "C" void update_xyzdir(int dir){
  xyz_blockage_dir=dir;
  switch (dir) {
  case 0:
    group2->set_int_val(dir);
    break;
  case 1:
    group2->set_int_val(dir);
    break;
  case 2:
    group2->set_int_val(dir);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
}

/* ------------------ OBJECT_CB ------------------------ */

void OBJECT_CB(int var){
  int i,temp;
  switch (var){
    case MESH_LIST:
      update_current_mesh(selected_case->meshinfo + highlight_mesh);
      update_rotation_index(highlight_mesh);
    break;
    case SET_LABEL:
    if(bchighlight!=NULL){
      obstlabelcopy(&bchighlight->label,edittextlabel->get_text());
    }
    break;
    case NEW_OBJECT:
      addnewobject();
      break;
    case DELETE_OBJECT:
      deleteobject();
      break;
    case UNDELETE_OBJECT:
      undeleteobject();
      break;
    case RADIO_XYZ:
      xyz_blockage_dir=xyz_dir;
      update_movestretch(bchighlight);
      update_xyzdir(xyz_dir);
      break;
    case RADIO_OBJECT_TYPE:
      xyz_dir=xyz_blockage_dir;
      group2->set_int_val(xyz_dir);
      break;
    case XMIN_SPIN:
    case XMAX_SPIN:
    case YMIN_SPIN:
    case YMAX_SPIN:
    case ZMIN_SPIN:
    case ZMAX_SPIN:
      movefblockage(
        &glui_block_xmin, &glui_block_xmax,
        &glui_block_ymin, &glui_block_ymax,
        &glui_block_zmin, &glui_block_zmax);
      update_blockvals(0);
      break;

    case UPDATE_LIST:
      switch (wall_case){
      case WALL_1:
        temp=surface_indices[UP_Z];
        if(nsurfaces>0){
          for(i=0;i<6;i++){
            if(i==UP_Z)continue;
            surface_indices[i]=temp;
            surfacelists[i]->set_int_val(temp);
          }
        }
        break;
      case WALL_3:
        temp=surface_indices[UP_Y];
        if(nsurfaces>0){
          for(i=0;i<6;i++){
            if(i==UP_Z||i==DOWN_Z||i==UP_Y)continue;
            surface_indices[i]=temp;
            surfacelists[i]->set_int_val(temp);
          }
        }
        break;
      case WALL_6:
        break;
      default:
        ASSERT(FFALSE);
        break;
      }

      if(bchighlight!=NULL){
        for(i=0;i<6;i++){
          bchighlight->surf[i]=surfaceinfo+sorted_surfidlist[surface_indices[i]];
          bchighlight->surf_index[i]=sorted_surfidlist[surface_indices[i]];
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
        surfacelists[DOWN_Z]->set_name("zmin");
        surfacelists[UP_Z]->set_name("zmax");
        surfacelists[DOWN_Y]->set_name("ymin");
        surfacelists[UP_Y]->set_name("ymax");
        surfacelists[DOWN_X]->set_name("xmin");
        surfacelists[UP_X]->set_name("xmax");
        break;
      case WALL_3:
        for(i=0;i<6;i++){
          surfacelists[i]->disable();
        }
        surfacelists[DOWN_Z]->enable();
        surfacelists[UP_Z]->enable();
        surfacelists[UP_Y]->enable();

        surfacelists[DOWN_Z]->set_name("zmin");
        surfacelists[UP_Z]->set_name("zmax");
        surfacelists[UP_Y]->set_name("vertical");

        surfacelists[DOWN_Y]->set_name("");
        surfacelists[DOWN_X]->set_name("");
        surfacelists[UP_X]->set_name("");
        break;
      case WALL_1:
        for(i=0;i<6;i++){
          surfacelists[i]->disable();
        }
        surfacelists[UP_Z]->enable();
        surfacelists[UP_Z]->set_name("All");

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
        blockage_checkbox->set_int_val(blockage_as_input);
      case BLOCKAGE_AS_INPUT:
        if(blockage_as_input==1){
          edittext_xmin->disable();
          edittext_ymin->disable();
          edittext_zmin->disable();

          edittext_xmax->disable();
          edittext_ymax->disable();
          edittext_zmax->disable();
          SPINNER_stretch_black->disable();
          SPINNER_stretch_white->disable();
          SPINNER_move->disable();
          blocklocation=BLOCKlocation_exact;
        }
        else{
          edittext_xmin->enable();
          edittext_ymin->enable();
          edittext_zmin->enable();

          edittext_xmax->enable();
          edittext_ymax->enable();
          edittext_zmax->enable();
          SPINNER_stretch_black->enable();
          SPINNER_stretch_white->enable();
          SPINNER_move->enable();
          blocklocation=BLOCKlocation_grid;
        }
        break;
    default:
      ASSERT(FFALSE);
      break;
  }
}
