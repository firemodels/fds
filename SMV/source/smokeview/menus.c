#include "options.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include GLUT_H

#include "string_util.h"
#include "update.h"
#include "smokeviewvars.h"
#include "IOvolsmoke.h"

#ifdef WIN32
#include <direct.h>
#endif

#define MENU_SMOKE3D_IBLANK -2

#define MENU_KEEP_ALL -2
#define MENU_KEEP_FINE -3
#define MENU_KEEP_COARSE -4

#define MENU_SLICECOLORDEFER -5

#define MENU_OPTION_TRAINERMENU 2

#define MENU_UPDATEBOUNDS -3

#define MENU_DUMMY3 -2

#define MENU_ERASECOMPRESS 1
#define MENU_OVERWRITECOMPRESS 2
#define MENU_COMPRESSNOW 3
#define MENU_COMPRESSAUTOLOAD 4

#define MENU_TRAINER_CLEAR 998
#define MENU_MAIN_QUIT 3

#define MENU_READCASEINI -1
#define MENU_READINI 1
#define MENU_WRITEINI 2
#define MENU_WRITECASEINI 3
#define MENU_READSVO 4

#define MENU_DUMMY2 -1

#define MENU_PLOT3D_DUMMY 997
#define MENU_PLOT3D_Z 1
#define MENU_PLOT3D_Y 2
#define MENU_PLOT3D_X 3
#define MENU_PLOT3D_CONT 4
#define MENU_PLOT3D_SHOWALL 5
#define MENU_PLOT3D_HIDEALL 6

#define MENU_MAIN_TRAINERTOGGLE 997

#define MENU_UNLOADSMOKE3D_UNLOADALLSOOT -1
#define MENU_UNLOADSMOKE3D_UNLOADALLFIRE -2
#define MENU_UNLOADSMOKE3D_UNLOADALLWATER -3

#define MENU_UNLOADTERRAIN_UNLOADALL -10
#define MENU_UNLOADTERRAIN_DUMMY -1

#define MENU_LOADTERRAIN_LOADALL -9
#define MENU_LOADTERRAIN_UNLOAD -10
#define MENU_LOADTERRAIN_DUMMY -1

#define MENU_LOADVSLICE_SHOWALL -20

#define MENU_EVAC_ALLMESHES -11
#define MENU_EVAC_UNLOADALL -1
#define MENU_EVAC_DUMMY -2

#define MENU_PARTICLE_UNLOAD -1
#define MENU_PARTICLE_DUMMY -2
#define MENU_PARTICLE_ALLMESHES -11

#define MENU_UNLOADEVAC_UNLOADALL -1

#define MENU_UNLOADPARTICLE_UNLOADALL -1

#define MENU_AVATAR_DEFINED -1

#define MENU_PARTSHOW_PARTICLES 1
#define MENU_PARTSHOW_DROPLETS 2
#define MENU_PARTSHOW_SHOWALL 3
#define MENU_PARTSHOW_HIDEALL 4
#define MENU_PARTSHOW_STATIC 5

#define MENU_PROP_DUMMY -1
#define MENU_PROP_SHOWALL -2
#define MENU_PROP_HIDEALL -3
#define MENU_PROP_HIDEPART -4
#define MENU_PROP_HIDEAVATAR -5
#define MENU_PROP_TRACERS -6

#define MENU_STREAK_HIDE -2
#define MENU_STREAK_HEAD -3

#define MENU_VECTOR_SHOW -2

#define MENU_SURFACE_SMOOTH 0
#define MENU_SURFACE_FACET 1
#define MENU_SURFACE_OUTLINE 2
#define MENU_SURFACE_POINTS 3

#define MENU_ISOSHOW_SHOWALL 99
#define MENU_ISOSHOW_HIDEALL 98
#define MENU_ISOSHOW_ALLSOLID 94
#define MENU_ISOSHOW_ALLTRANSPARENT 95
#define MENU_ISOSHOW_MINSOLID 96
#define MENU_ISOSHOW_MAXSOLID 97

#define MENU_ISOSHOW_SOLID 1
#define MENU_ISOSHOW_OUTLINE 2
#define MENU_ISOSHOW_POINTS 3
#define MENU_ISOSHOW_SMOOTH 4
#define MENU_ISOSHOW_NORMALS 5

#define MENU_ZONE_2DTEMP 6
#define MENU_ZONE_2DHAZARD 5
#define MENU_ZONE_3DSMOKE 7
#define MENU_ZONE_HORIZONTAL 1
#define MENU_ZONE_VERTICAL 2
#define MENU_ZONE_LAYERHIDE 4
#define MENU_ZONE_VENTS 14
#define MENU_ZONE_HVENTS 15
#define MENU_ZONE_VVENTS 16
#define MENU_ZONE_MVENTS 17
#define MENU_ZONE_FIRES 18
#define MENU_ZONE_VENT_SLAB 19
#define MENU_ZONE_VENT_PROFILE 20

#define MENU_SHOWSLICE_INBLOCKAGE -11
#define MENU_SHOWSLICE_SLICEANDVECTORS -15
#define MENU_SHOWSLICE_TERRAIN -13
#define MENU_SHOWSLICE_OFFSET -12
#define MENU_SHOWSLICE_FEDAREA -14

#define MENU_VENT_OPEN 14
#define MENU_VENT_EXTERIOR 16
#define MENU_VENT_OTHER 21
#define MENU_VENT_OUTLINE 15
#define MENU_VENT_TWOINTERIOR 18
#define MENU_VENT_TWOEXTERIOR 19
#define MENU_VENT_TRANSPARENT 20
#define MENU_VENT_CIRCLE 23
#define MENU_VENT_RECTANGLE 24
#define MENU_VENT_CIRCLEHIDE 25
#define MENU_VENT_CIRCLEOUTLINE 26

#define MENU_TIMEVIEW -103
#define SAVE_VIEWPOINT -101
#define MENU_STARTUPVIEW -102
#define MENU_OUTLINEVIEW -104
#define MENU_SIZEPRESERVING -105
#define MENU_DUMMY -999

#define MENU_SHOWHIDE_EVAC 13
#define MENU_SHOWHIDE_PRINT 16
#define MENU_SHOWHIDE_PARTICLES 1
#define MENU_SHOWHIDE_SENSOR 9
#define MENU_SHOWHIDE_SENSOR_NORM 14
#define MENU_SHOWHIDE_OFFSET 12

#define MENU_UNITS_RESET -1
#define MENU_UNITS_SHOWALL -3
#define MENU_UNITS_HMS -2

#define MENU_HELP_ISSUES -1
#define MENU_HELP_DOWNLOADS -2
#define MENU_HELP_DOCUMENTATION -3
#define MENU_HELP_FDSWEB -4

#define GRID_yz 1
#define GRID_xz 2
#define GRID_xy 3
#define GRID_showall 4
#define GRID_hideall 5
#define GRID_grid 7
#define GRID_probe 8

#define OBJECT_SHOWALL -1
#define OBJECT_HIDEALL -2
#define OBJECT_SELECT -3
#define OBJECT_OUTLINE -4
#define OBJECT_ORIENTATION -5
#define OBJECT_MISSING -6

void AddScriptList(char *file, int id);
void UpdateGluiRender(void);
void PropMenu(int value);
void UnLoadVolSmoke3DMenu(int value);
void LoadVolSmoke3DMenu(int value);
void UpdateScriptStep(void);
#define ISO_COLORS 4
void IsoCB(int var);
void UpdateSliceDups(void);
void UpdateVSliceDups(void);
void UnloadVSliceMenu(int value);


#ifdef WIN32

/* ------------------ OpenSMVFile ------------------------ */

void OpenSMVFile(char *filebuffer,int filebufferlength,int *openfile){
  char stringFilter[]="Smokeview Files (*.smv)\0*.smv\0\0\0";
  char strTitle[80]="Select Smokeview Case";
  int buffersize;
  char smv_directory[1024];
  OPENFILENAME fileinfo;

  *openfile=0;
  buffersize=sizeof(OPENFILENAME);

  STRCPY(filebuffer,"");
  fileinfo.lStructSize=(unsigned long)buffersize;
  fileinfo.hwndOwner=NULL;
  fileinfo.lpstrFilter=stringFilter;
  fileinfo.lpstrCustomFilter=NULL;
  fileinfo.lpstrFile=filebuffer;
  fileinfo.nMaxFile=(unsigned long)filebufferlength;
  fileinfo.lpstrFileTitle=NULL;
  fileinfo.nMaxFileTitle=80;
  fileinfo.lpstrInitialDir=NULL;
  fileinfo.lpstrTitle=strTitle;
  fileinfo.Flags=0;
  fileinfo.lpstrDefExt=NULL;

  if(GetOpenFileName(&fileinfo)){
    STRCPY(smv_directory,"");
    strncat(smv_directory,filebuffer,fileinfo.nFileOffset);
    if( _chdir( smv_directory )   ){
      PRINTF( "Unable to locate the directory: %s\n", smv_directory );
    }
    else{
      *openfile=1;
    }
  }
}
#endif

/* ------------------ HideAllSmoke ------------------------ */

void HideAllSmoke(void){
  int i;
  for(i = 0; i < nsmoke3dinfo; i++){
    smoke3ddata *smoke3di;

    smoke3di = smoke3dinfo + i;
    if(smoke3di->loaded == 1)smoke3di->display = 0;
  }
  for(i = 0; i < nisoinfo; i++){
    isodata *isoi;

    isoi = isoinfo + i;
    if(isoi->loaded == 1)isoi->display = 0;
  }
}

/* ------------------ HideAllSlices ------------------------ */

void HideAllSlices(void){
  int i;

  glutSetCursor(GLUT_CURSOR_WAIT);
  for(i = 0; i < nsliceinfo; i++){
    sliceinfo[i].display = 0;
  }
  updatemenu = 1;
  glutPostRedisplay();
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ ShowAllSmoke ------------------------ */

void ShowAllSmoke(void){
  int i;
  for(i = 0; i < nsmoke3dinfo; i++){
    smoke3ddata *smoke3di;

    smoke3di = smoke3dinfo + i;
    if(smoke3di->loaded == 1)smoke3di->display = 1;
  }
  for(i = 0; i < nisoinfo; i++){
    isodata *isoi;

    isoi = isoinfo + i;
    if(isoi->loaded == 1)isoi->display = 1;
  }
}

/* ------------------ ShowMultiSliceMenu ------------------------ */

void ShowMultiSliceMenu(int value){
  multislicedata *mslicei;
  slicedata *sd;
  int mdisplay;
  int i;

  updatemenu = 1;
  glutPostRedisplay();
  switch(value){
  case SHOW_ALL:
  case HIDE_ALL:
    ShowHideSliceMenu(value);
    return;
  case MENU_SHOWSLICE_INBLOCKAGE:
    show_slice_in_obst = 1 - show_slice_in_obst;
	update_show_slice_in_obst();
    break;
  case -12:
    offset_slice = 1 - offset_slice;
    break;
  case -14:
    show_fed_area = 1 - show_fed_area;
    break;
  default:
    mslicei = multisliceinfo + value;
    mdisplay = 0;
    if(islicetype == mslicei->type){
      if(plotstate != DYNAMIC_PLOTS){
        plotstate = DYNAMIC_PLOTS;
        mdisplay = 1;
      }
      else{
        mdisplay = 1 - mslicei->display;
      }
    }
    else{
      plotstate = DYNAMIC_PLOTS;
      islicetype = mslicei->type;
      mdisplay = 1;
    }
    for(i = 0; i < mslicei->nslices; i++){
      sd = sliceinfo + mslicei->islices[i];
      if(sd->loaded == 0)continue;
      sd->display = mdisplay;
    }
    break;
  }
  updateslicefilenum();
  plotstate = GetPlotState(DYNAMIC_PLOTS);

  updateglui();
  updateslicelistindex(slicefilenum);
  UpdateShow();
}

/* ------------------ ShowAllSlices ------------------------ */

void ShowAllSlices(char *type1, char *type2){
  int i;

  glutSetCursor(GLUT_CURSOR_WAIT);
  if(trainer_showall_mslice == 1){
    for(i = 0; i < nsliceinfo; i++){
      sliceinfo[i].display = 0;
      if(sliceinfo[i].loaded == 0)continue;
      if(
        type1 != NULL&&STRCMP(sliceinfo[i].label.longlabel, type1) == 0 ||
        type2 != NULL&&STRCMP(sliceinfo[i].label.longlabel, type2) == 0
        ){
        sliceinfo[i].display = 1;
        islicetype = sliceinfo[i].type;
      }
    }
  }
  else{
    int msliceindex;

    if(trainerload == 2){
      if(trainerload == trainerload_old){
        trainer_temp_index++;
        if(trainer_temp_index > trainer_temp_n - 1){
          trainer_temp_index = 0;
        }
      }
      msliceindex = trainer_temp_indexes[trainer_temp_index];
    }
    else{
      if(trainerload == trainerload_old){
        trainer_oxy_index++;
        if(trainer_oxy_index > trainer_oxy_n - 1){
          trainer_oxy_index = 0;
        }
      }
      msliceindex = trainer_oxy_indexes[trainer_oxy_index];
    }
    ShowMultiSliceMenu(HIDE_ALL);
    ShowMultiSliceMenu(msliceindex);
  }
  updatemenu = 1;
  glutPostRedisplay();
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ TrainerViewMenu ------------------------ */

void TrainerViewMenu(int value){
  switch(value){
  case MENU_TRAINER_smoke:   // realistic
    HideAllSlices();
    trainerload=1;
    ShowAllSmoke();
    trainerload_old=1;
    break;
  case MENU_TRAINER_temp:  // temperature slices
    HideAllSmoke();
    trainerload=2;
    ShowAllSlices("TEMPERATURE",NULL);
    trainerload_old=2;
    break;
  case MENU_TRAINER_oxy:  //  oxygen slices
    HideAllSmoke();
    trainerload=3;
    ShowAllSlices("OXYGEN","OXYGEN VOLUME FRACTION");
    trainerload_old=3;
    break;
  case MENU_TRAINER_CLEAR: // unload
    LoadUnloadMenu(UNLOADALL);
    trainerload=0;
    trainerload_old=0;
    break;
  default:
    ASSERT(FFALSE);
  }
  updatemenu=1;
  glutPostRedisplay();
}

/* ------------------ MainMenu ------------------------ */

void MainMenu(int value){

  if(value==MENU_MAIN_QUIT){
    if(scriptoutstream!=NULL){
      ScriptMenu(SCRIPT_STOP_RECORDING);
    }
    exit(0);
  }
  if(value==MENU_MAIN_TRAINERTOGGLE){
    trainer_mode=1-trainer_mode;
  }
  updatemenu=1;
  glutPostRedisplay();
}

/* ------------------ StaticVariableMenu ------------------------ */

void StaticVariableMenu(int value){

  plotn=value;
  plotstate=STATIC_PLOTS;
  visGrid=0;
  if(visiso==1){
    updateshowstep(1,ISO);
  }
  updatesurface();
  if(visx_all==1){
    updateshowstep(1,XDIR);
  }
  if(visy_all==1){
    updateshowstep(1,YDIR);
  }
  if(visz_all==1){
    updateshowstep(1,ZDIR);
  }
  if(visx_all==0&&visy_all==0&&visz_all==0){
    updateshowstep(1,YDIR);
  }
  updateallplotslices();
  updatemenu=1;
  glutPostRedisplay();
  updateplot3dlistindex();
}

/* ------------------ IsoVariableMenu ------------------------ */

void IsoVariableMenu(int value){
  if(ReadPlot3dFile==1){
    plotn=value;
    if(visx_all==1){
      updateshowstep(1,XDIR);
    }
    if(visy_all==1){
      updateshowstep(1,YDIR);
    }
    if(visz_all==1){
      updateshowstep(1,ZDIR);
    }
    updateshowstep(1,ISO);
    updatesurface();
    plotstate=STATIC_PLOTS;
    updateplotslice(XDIR);
    updateplotslice(YDIR);
    updateplotslice(ZDIR);
    updatemenu=1;
    glutPostRedisplay();
    updateplot3dlistindex();
  }
}

/* ------------------ LabelMenu ------------------------ */

void LabelMenu(int value){
  updatemenu=1;
  if(value == MENU_DUMMY)return;
  glutPostRedisplay();
  switch(value){
  case MENU_LABEL_colorbar:
    visColorbar=1-visColorbar;
    break;
  case MENU_LABEL_timebar:
    visTimebar=1-visTimebar;
    break;
  case MENU_LABEL_title:
    visTitle=1-visTitle;
    break;
  case MENU_LABEL_framerate:
    visFramerate = 1 - visFramerate;
    break;
   case MENU_LABEL_ShowAll:
    visUSERticks=1;
    visColorbar=1;
    visTimebar=1;
    visTitle=1;
    visFramerate=1;
#ifdef pp_memstatus
    visAvailmemory=1;
#endif
    visaxislabels=1;
    visTimelabel=1;
    visFramelabel=1;
    visLabels=1;
    visMeshlabel=1;
    vis_slice_average=1;
    if(ntickinfo>0)visFDSticks=1;
    visgridloc=1;
    visHRRlabel=1;
    show_hrrcutoff=1;
    visFramelabel=1;
	  if(hrrinfo != NULL&&hrrinfo->display != 1)UpdateHrrinfo(1);
    gversion=1;
    break;
   case MENU_LABEL_HideAll:
    visUSERticks=0;
    visColorbar=0;
    visTimebar=0;
    visTitle=0;
    visFramerate=0;
    visaxislabels=0;
    visLabels=0;
    visTimelabel=0;
    visFramelabel=0;
    visMeshlabel=0;
    visHRRlabel=0;
    show_hrrcutoff=0;
	if (hrrinfo != NULL&&hrrinfo->display != 0)UpdateHrrinfo(0);
    if(ntickinfo>0)visFDSticks=0;
    visgridloc=0;
    vis_slice_average=0;
    gversion=0;
#ifdef pp_memstatus
    visAvailmemory=0;
#endif
    break;
   case MENU_LABEL_northangle:
     vis_northangle = 1-vis_northangle;
     break;
   case MENU_LABEL_axis:
    visaxislabels = 1 - visaxislabels;
    update_visaxislabels();
    break;
   case MENU_LABEL_textlabels:
     visLabels = 1 - visLabels;
     break;
   case MENU_LABEL_timelabel:
     visTimelabel=1-visTimelabel;
     if(visTimelabel==1)visTimebar=1;
     break;
   case MENU_LABEL_framelabel:
     visFramelabel=1-visFramelabel;
     if(visFramelabel==1)visTimebar=1;
     if(visFramelabel==1){
       visHRRlabel=0;
	   UpdateHrrinfo(visHRRlabel);
     }
     break;
   case MENU_LABEL_meshlabel:
     visMeshlabel=1-visMeshlabel;
     break;
#ifdef pp_memstatus
   case MENU_LABEL_memload:
     visAvailmemory = 1 - visAvailmemory;
     break;
#endif
#ifdef pp_MEMDEBUG
   case MENU_LABEL_memusage:
     visUsagememory = 1 - visUsagememory;
#ifdef pp_memstatus
     if(visUsagememory==1)visAvailmemory=0;
#endif
     break;
#endif
   case MENU_LABEL_fdsticks:
     visFDSticks=1-visFDSticks;
     break;
   case MENU_LABEL_hmslabel:
     vishmsTimelabel = 1 - vishmsTimelabel;
     break;
   case MENU_LABEL_grid:
     visgridloc = 1 - visgridloc;
     break;
   case MENU_LABEL_sliceaverage:
     vis_slice_average = 1 - vis_slice_average;
     break;
   case MENU_LABEL_hrr:
     visHRRlabel=1-visHRRlabel;
     if(visHRRlabel==1)visTimebar=1;
	 UpdateHrrinfo(visHRRlabel);
     break;
   case MENU_LABEL_hrrcutoff:
     show_hrrcutoff=1-show_hrrcutoff;
     break;
   case MENU_LABEL_userticks:
     visUSERticks = 1 - visUSERticks;
     break;
   case MENU_LABEL_gversion:
     gversion=1-gversion;
     break;
   default:
     ASSERT(FFALSE);
     break;
  }
  set_labels_controls();
}

/* ------------------ SmokeColorbarMenu ------------------------ */

void SmokeColorbarMenu(int value){
  if(value==MENU_DUMMY)return;
  updatemenu=1;

  value = CLAMP(value, 0, ncolorbars - 1);
  fire_colorbar_index=value;
  fire_colorbar = colorbarinfo + value;
  UpdateRGBColors(COLORBAR_INDEX_NONE);
  if(FlowDir>0){
    keyboard('-',FROM_SMOKEVIEW);
    keyboard(' ',FROM_SMOKEVIEW);
  }
  else{
    keyboard(' ',FROM_SMOKEVIEW);
    keyboard('-',FROM_SMOKEVIEW);
  }
  glutPostRedisplay();
}

/* ------------------ ColorbarMenu ------------------------ */

void ColorbarMenu(int value){
  if(value==MENU_DUMMY)return;
  updatemenu=1;
  glutPostRedisplay();
  if(value<0){
    switch(value){
    case COLORBAR_FLIP:
      colorbarflip=1-colorbarflip;
      update_colorbarflip();
      break;
    case COLORBAR_RESET:
      show_extreme_mindata=0;
      show_extreme_maxdata=0;
      colorbarflip=0;
      contour_type=SHADED_CONTOURS;
      setbw=0;
      update_extreme();
      UpdateRGBColors(COLORBAR_INDEX_NONE);
      break;
    case COLORBAR_HIGHLIGHT_BELOW:
      show_extreme_mindata=1-show_extreme_mindata;
      update_extreme();
      UpdateRGBColors(COLORBAR_INDEX_NONE);
      break;
    case COLORBAR_HIGHLIGHT_ABOVE:
      show_extreme_maxdata=1-show_extreme_maxdata;
      update_extreme();
      UpdateRGBColors(COLORBAR_INDEX_NONE);
      break;
    case COLORBAR_TOGGLE_BW_DATA:
      setbwdata = 1 - setbwdata;
      if(setbwdata==1){
        colorbartype_save=colorbartype;
        ColorbarMenu(bw_colorbar_index);
      }
      else{
        ColorbarMenu(colorbartype_save);
      }
      IsoCB(ISO_COLORS);
      break;
    case COLORBAR_TOGGLE_BW:
      setbw=1-setbw;
      InitRGB();
      set_labels_controls();
      break;
   case COLORBAR_TRANSPARENT:
     use_transparency_data=1-use_transparency_data;
     UpdateRGBColors(COLORBAR_INDEX_NONE);
     set_labels_controls();
     update_transparency();
     break;
   case COLORBAR_CONTINUOUS:
     contour_type=SHADED_CONTOURS;
     UpdateRGBColors(COLORBAR_INDEX_NONE);
     break;
   case COLORBAR_STEPPED:
     contour_type=STEPPED_CONTOURS;
     UpdateRGBColors(COLORBAR_INDEX_NONE);
     break;
   case COLORBAR_LINES:
     contour_type=LINE_CONTOURS;
     UpdateRGBColors(COLORBAR_INDEX_NONE);
     break;
   default:
     ASSERT(FFALSE);
     break;
   }
  }
  if(value>=0){
    colorbartype=value;
    selectedcolorbar_index2=colorbartype;
    UpdateCurrentColorbar(colorbarinfo+colorbartype);
    update_colorbar_type();
    update_colorbar_list2();
    if(colorbartype == bw_colorbar_index){
      setbwdata = 1;
    }
    else{
      setbwdata = 0;
    }
    IsoCB(ISO_COLORS);
    set_labels_controls();
  }
  if(value>-10){
    UpdateRGBColors(COLORBAR_INDEX_NONE);
  }
}

/* ------------------ Smoke3DShowMenu ------------------------ */

void Smoke3DShowMenu(int value){
  smoke3ddata *smoke3di;
  int i;

  updatemenu=1;
  glutPostRedisplay();
  if(value<0){
    switch(value){
    case SHOW_ALL:
      plotstate=DYNAMIC_PLOTS;
      for(i=0;i<nsmoke3dinfo;i++){
        smoke3di = smoke3dinfo + i;
        if(smoke3di->loaded==1)smoke3di->display=1;
      }
      break;
    case HIDE_ALL:
      for(i=0;i<nsmoke3dinfo;i++){
        smoke3di = smoke3dinfo + i;
        if(smoke3di->loaded==1)smoke3di->display=0;
      }
      break;
    case HAVE_LIGHT:
      show_smoke_lighting = 1 - show_smoke_lighting;
      break;
    default:
      ASSERT(FFALSE);
    }
  }
  else{
    smoke3di = smoke3dinfo + value;
    if(plotstate!=DYNAMIC_PLOTS){
      plotstate=DYNAMIC_PLOTS;
      smoke3di->display=1;
    }
    else{
      smoke3di->display = 1 - smoke3di->display;
    }
  }
}

/* ------------------ IsoShowMenu ------------------------ */

void IsoShowMenu(int value){
  int i;
  int nisolevels, *showlevels;
  isodata *isoi;

  nisolevels=loaded_isomesh->nisolevels;
  showlevels=loaded_isomesh->showlevels;

  switch(value){
  case  MENU_ISOSHOW_SMOOTH:
    smooth_iso_normal=1-smooth_iso_normal;
    break;
  case MENU_ISOSHOW_NORMALS:
    show_iso_normal = 1 - show_iso_normal;
    break;
  case MENU_ISOSHOW_SOLID:
  case MENU_ISOSHOW_OUTLINE:
  case MENU_ISOSHOW_POINTS:
    if(value == MENU_ISOSHOW_SOLID)show_iso_solid=1-show_iso_solid;
    if(value == MENU_ISOSHOW_OUTLINE)show_iso_outline = 1 - show_iso_outline;
    if(value == MENU_ISOSHOW_POINTS)show_iso_verts = 1 - show_iso_verts;
    visAIso=show_iso_solid*1+show_iso_outline*2+show_iso_verts*4;
    if(visAIso!=0){
      plotstate=DYNAMIC_PLOTS;
    }
    update_glui_isotype();
    break;
   case MENU_ISOSHOW_ALLSOLID:
    transparent_state=ALL_SOLID;
    if(loaded_isomesh==NULL)break;
    for(i=0;i<loaded_isomesh->nisolevels;i++){
      surfdata *surfi;

      surfi = surfinfo + nsurfinfo + 1 + i;
      surfi->transparent_level=1.0;
    }
    use_transparency_data=0;
    break;
   case MENU_ISOSHOW_ALLTRANSPARENT:
    transparent_state=ALL_TRANSPARENT;
    if(loaded_isomesh==NULL)break;
    for(i=0;i<loaded_isomesh->nisolevels;i++){
      surfdata *surfi;

      surfi = surfinfo + nsurfinfo + 1 + i;
      surfi->transparent_level=transparent_level;
    }
    use_transparency_data=1;
    break;
   case MENU_ISOSHOW_MINSOLID:
    transparent_state=MIN_SOLID;
    if(loaded_isomesh==NULL)break;
    for(i=0;i<loaded_isomesh->nisolevels;i++){
      surfdata *surfi;

      surfi = surfinfo + nsurfinfo + 1 + i;
      surfi->transparent_level=transparent_level;
    }
    surfinfo[nsurfinfo+1].transparent_level=1.0;
    use_transparency_data=1;
    break;
   case MENU_ISOSHOW_MAXSOLID:
    transparent_state=MAX_SOLID;
    if(loaded_isomesh==NULL)break;
    for(i=0;i<loaded_isomesh->nisolevels;i++){
      surfdata *surfi;

      surfi = surfinfo + nsurfinfo + 1 + i;
      surfi->transparent_level=transparent_level;
    }
    use_transparency_data=1;
    surfinfo[nsurfinfo+1+loaded_isomesh->nisolevels-1].transparent_level=1.0;
    break;
   case MENU_ISOSHOW_HIDEALL:
    show_iso_solid=0;
    show_iso_outline=0;
    show_iso_verts=0;
    visAIso=show_iso_solid*1+show_iso_outline*2+show_iso_verts*4;
    for(i=0;i<nisolevels;i++){
      showlevels[i]=0;
    }
    break;
   case MENU_ISOSHOW_SHOWALL:
    show_iso_solid=1;
    show_iso_outline=0;
    show_iso_verts=0;
    visAIso=show_iso_solid*1+show_iso_outline*2+show_iso_verts*4;
    for(i=0;i<nisolevels;i++){
      showlevels[i]=1;
    }
    break;
   default:
    if(value>99&&value<999&&value-100<nisolevels){
     showlevels[value-100] = 1 - showlevels[value-100];
    }
    else if(value>=1000&&value<=10000){      // we can only have 9900 isosurface files
     isoi = isoinfo + value - 1000;          // hope that is enough!!
     if(plotstate!=DYNAMIC_PLOTS){
       plotstate=DYNAMIC_PLOTS;
       isoi->display=1;
       iisotype=isoi->type;
     }
     else{
       if(isoi->type==iisotype){
         isoi->display = 1 - isoi->display;
         update_isotype();
       }
       else{
         isoi->display=1;
         iisotype=isoi->type;
       }
     }
     UpdateShow();
    }
    else if(value>=SHOWALL_ISO){
      if(value==SHOWALL_ISO){
        plotstate=DYNAMIC_PLOTS;
        for(i=0;i<nisoinfo;i++){
          isoinfo[i].display=1;
        }
      }
      else if(value==HIDEALL_ISO){
        for(i=0;i<nisoinfo;i++){
          isoinfo[i].display=0;
        }
      }
     UpdateShow();
    }
  }
  update_iso_showlevels();
  Update_Isotris(1);

  updatemenu=1;
  glutPostRedisplay();
}

/* ------------------ ShowVSliceMenu ------------------------ */

void ShowVSliceMenu(int value){
  int i;
  vslicedata *vd;

  if(value == MENU_DUMMY)return;
  updatemenu = 1;
  glutPostRedisplay();
  if(value==SHOW_ALL){
    for(i=0;i<nvsliceinfo;i++){
      vd = vsliceinfo+i;
      if(vd->loaded==0)continue;
      vd->display=1;
    }
    show_all_slices=1;
    UpdateTimes();
    return;
  }
  if(value==HIDE_ALL){
    for(i=0;i<nvsliceinfo;i++){
      vd = vsliceinfo+i;
      if(vd->loaded==0)continue;
      vd->display=0;
    }
    show_all_slices=0;
    UpdateTimes();
    return;
  }
  if(value == MENU_SHOWSLICE_INBLOCKAGE){
    show_slice_in_obst=1-show_slice_in_obst;
	update_show_slice_in_obst();
    return;
  }
  if(value == MENU_SHOWSLICE_OFFSET){
    offset_slice=1-offset_slice;
    return;
  }
  if(value == MENU_SHOWSLICE_SLICEANDVECTORS){
    show_slices_and_vectors=1-show_slices_and_vectors;
    return;
  }
  vd = vsliceinfo + value;
  if(islicetype==sliceinfo[vd->ival].type){
    if(plotstate!=DYNAMIC_PLOTS){
      plotstate=DYNAMIC_PLOTS;
      vd->display=1;
    }
    else{
      vd->display = 1 - vd->display;
    }
    if(vd->iu!=-1){
      slicedata *sd;

      sd = sliceinfo+vd->iu;
      sd->vloaded=vd->display;
    }
    if(vd->iv!=-1){
      slicedata *sd;

      sd = sliceinfo+vd->iv;
      sd->vloaded=vd->display;
    }
    if(vd->iw!=-1){
      slicedata *sd;

      sd = sliceinfo+vd->iw;
      sd->vloaded=vd->display;
    }
    if(vd->ival!=-1){
      slicedata *sd;

      sd = sliceinfo+vd->ival;
      sd->vloaded=vd->display;
    }
  }
  else{
    islicetype = sliceinfo[vd->ival].type;
    vd->display=1;
  }
  plotstate=GetPlotState(DYNAMIC_PLOTS);
  UpdateShow();
}

/* ------------------ ShowHideSliceMenu ------------------------ */

void ShowHideSliceMenu(int value){
  int i;

  if(value == MENU_DUMMY)return;
  updatemenu=1;
  glutPostRedisplay();
  if(value<0){
    switch(value){
    case SHOW_ALL:
      for(i=0;i<nsliceinfo;i++){
        sliceinfo[i].display=1;
      }
      show_all_slices=1;
      break;
    case HIDE_ALL:
      for(i=0;i<nsliceinfo;i++){
        sliceinfo[i].display=0;
      }
      show_all_slices=0;
      break;
    case MENU_SHOWSLICE_INBLOCKAGE:
      show_slice_in_obst=1-show_slice_in_obst;
      update_show_slice_in_obst();
      break;
    case MENU_SHOWSLICE_OFFSET:
      offset_slice=1-offset_slice;
      break;
    case MENU_SHOWSLICE_TERRAIN:
      planar_terrain_slice=1-planar_terrain_slice;
      break;
    case MENU_SHOWSLICE_FEDAREA:
      show_fed_area=1-show_fed_area;
      break;
    case MENU_SHOWSLICE_SLICEANDVECTORS:
      show_slices_and_vectors=1-show_slices_and_vectors;
      return;
    default:
      ASSERT(FFALSE);
    }
  }
  else{
    slicedata *sd;

    sd = sliceinfo + value;
    if(islicetype==sd->type){
      if(plotstate!=DYNAMIC_PLOTS){
        plotstate=DYNAMIC_PLOTS;
        sd->display=1;
      }
      else{
        sd->display = 1 - sd->display;
      }
    }
    else{
      plotstate=DYNAMIC_PLOTS;
      islicetype = sd->type;
      sd->display=1;
    }
  }
  updateslicefilenum();
  plotstate=GetPlotState(DYNAMIC_PLOTS);

  updateglui();
  updateslicelistindex(slicefilenum);
  UpdateShow();
}

/* ------------------ ShowHideMenu ------------------------ */

void ShowHideMenu(int value){
  if(value==MENU_DUMMY)return;
  updatemenu=1;
  glutPostRedisplay();
  switch(value){
#ifdef pp_MEMPRINT
  case MENU_SHOWHIDE_PRINT:
    PrintMemoryInfo;
    break;
#endif
  case MENU_SHOWHIDE_FLIP:
   background_flip = 1-background_flip;
   UpdateRGBColors(COLORBAR_INDEX_NONE);
   set_labels_controls();
   break;
  case MENU_SHOWHIDE_EVAC:
    if(plotstate==DYNAMIC_PLOTS){
      visEvac=1-visEvac;
    }
    else{
      plotstate=DYNAMIC_PLOTS;
      visEvac=1;
    }
    UpdateTimes();
    break;
  case MENU_SHOWHIDE_PARTICLES:
    if(plotstate==DYNAMIC_PLOTS){
      visParticles=1-visParticles;
    }
    else{
      plotstate=DYNAMIC_PLOTS;
      visParticles=1;
    }
    UpdateTimes();
    break;
  case MENU_SHOWHIDE_SENSOR:
    visSensor=1-visSensor;
    break;
  case MENU_SHOWHIDE_SENSOR_NORM:
    visSensorNorm=1-visSensorNorm;
    break;
  case MENU_SHOWHIDE_OFFSET:
    if(titlesafe_offset==0){
      titlesafe_offset=titlesafe_offsetBASE;
    }
    else{
      titlesafe_offset=0;
    }
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
}

/* ------------------ ViewpointMenu ------------------------ */

void ViewpointMenu(int value){
  if(value==MENU_DUMMY)return;
  updatemenu=1;
  glutPostRedisplay();
  switch(value){
  case TOGGLE_TITLE_SAFE:
    if(titlesafe_offset==0){
      titlesafe_offset=titlesafe_offsetBASE;
    }
    else{
      titlesafe_offset=0;
    }
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
}

/* ------------------ DialogMenu ------------------------ */

void DialogMenu(int value){
  glutPostRedisplay();
  switch(value){
  case DIALOG_SHOOTER:
    show_glui_shooter();
    break;
  case DIALOG_TRAINER:
    show_glui_trainer();
    break;
  case DIALOG_DEVICE:
    show_glui_device();
    break;
  case DIALOG_3DSMOKE:
  case DIALOG_AUTOLOAD:
  case DIALOG_BOUNDS:
  case DIALOG_CONFIG:
  case DIALOG_SCRIPT:
  case DIALOG_SHOWFILES:
  case DIALOG_SMOKEZIP:
  case DIALOG_TIME:
    show_glui_bounds(value);
    break;
  case DIALOG_MOTION:
  case DIALOG_RENDER:
  case DIALOG_MOVIE:
  case DIALOG_SCALING:
  case DIALOG_VIEW:
  case DIALOG_WINDOW:
    show_glui_motion(value);
    break;
  case DIALOG_TICKS:
  case DIALOG_COLORING:
  case DIALOG_FONTS:
  case DIALOG_LABELS:
  case DIALOG_DISPLAY:
    show_glui_display(value);
    break;
  case DIALOG_TOUR:
   show_glui_tour();
   break;
  case DIALOG_CLIP:
    show_glui_clip();
    break;
  case DIALOG_STEREO:
    show_glui_stereo();
    break;
  case DIALOG_WUI:
    show_glui_wui();
    break;
  case DIALOG_COLORBAR:
    showcolorbar_dialog=1-showcolorbar_dialog;
    if(showcolorbar_dialog==1){
      show_glui_colorbar();
    }
    if(showcolorbar_dialog==0){
      hide_glui_colorbar();
    }
    break;
  case DIALOG_GEOMETRY:
    showedit_dialog=1-showedit_dialog;
    if(showedit_dialog==1){
      if(fds_filein!=NULL&&updategetobstlabels==1){
        CheckMemoryOff;
        GetObstLabels(fds_filein);
        CheckMemoryOn;
        updategetobstlabels=0;
      }
      visBlocksSave=visBlocks;
      show_glui_geometry();
      visBlocks=visBLOCKNormal;
    }
    if(showedit_dialog==0){
      hide_glui_geometry();
      visBlocks=visBlocksSave;
    }
    update_trainer_outline();

    break;
  case DIALOG_HIDEALL:
    showcolorbar_dialog = 0;
    hide_glui_shooter();
    hide_glui_display();
    hide_glui_bounds();
    hide_glui_motion();
    hide_glui_tour();
    hide_glui_clip();
    hide_glui_wui();
    hide_glui_stereo();
    hide_glui_colorbar();
    if(showedit_dialog==1)DialogMenu(DIALOG_GEOMETRY);
    hide_glui_trainer();
    hide_glui_device();
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  updatemenu=1;
}

/* ------------------ ZoomMenu ------------------------ */

void ZoomMenu(int value){
  if(value==MENU_DUMMY)return;
  updatemenu=1;
  if(opengldefined==1){
    glutPostRedisplay();
  }
  zoomindex=value;
  if(zoomindex==-1){
    if(zoom<zooms[0]){
      zoom=zooms[0];
      zoomindex=0;
    }
    if(zoom>zooms[4]){
      zoom=zooms[4];
      zoomindex=4;
    }
    if(projection_type!=0){
      camera_current->projection_type=projection_type;
      SetViewPoint(RESTORE_EXTERIOR_VIEW);
      update_projection_type();
    }
  }
  else if(zoomindex==UPDATE_PROJECTION){
    camera_current->projection_type=projection_type;
    update_projection_type();
    if(projection_type==0){
      UpdateCameraYpos(camera_current);
    }
    else{
      camera_current->eye[1]=camera_current->isometric_y;
    }

  }
  else{
    if(zoomindex<0)zoomindex=2;
    if(zoomindex>4)zoomindex=2;
    zoom=zooms[zoomindex];
    if(projection_type!=0){
      SetViewPoint(RESTORE_EXTERIOR_VIEW_ZOOM);
      camera_current->projection_type=projection_type;
      update_projection_type();
    }
  }
  camera_current->zoom=zoom;
  update_glui_zoom();
}

/* ------------------ ApertureMenu ------------------------ */

void ApertureMenu(int value){
  updatemenu=1;
  if(opengldefined==1){
    glutPostRedisplay();
  }
  apertureindex = CLAMP(value, 0, 4);
  aperture=apertures[apertureindex];
}

/* ------------------ FontMenu ------------------------ */

void FontMenu(int value){
  updatemenu=1;
  if(opengldefined==1){
    glutPostRedisplay();
  }
  switch(value){
  case SMALL_FONT:
    fontindex=SMALL_FONT;
    large_font=GLUT_BITMAP_HELVETICA_12;
    small_font=GLUT_BITMAP_HELVETICA_10;
    large_font_height=12;
    small_font_height=10;
    break;
  case LARGE_FONT:
    fontindex=LARGE_FONT;
    large_font=GLUT_BITMAP_HELVETICA_18;
    small_font=GLUT_BITMAP_HELVETICA_18;
    large_font_height=18;
    small_font_height=18;
    break;
  case SCALED_FONT:
    fontindex=SCALED_FONT;
    break;
  default:
    ASSERT(FFALSE);
  }
  glui_update_fontindex();
  set_labels_controls();
}

/* ------------------ UnitsMenu ------------------------ */

void UnitsMenu(int value){
  int unitclass, unit_index;
  int i;

  unitclass = value/1000;
  unit_index = value - unitclass*1000;
  unitclasses[unitclass].unit_index=unit_index;
  if(value==MENU_UNITS_RESET){
    for(i=0;i<nunitclasses;i++){
      unitclasses[i].unit_index=0;
    }
  }
  else if(value==MENU_UNITS_HMS){
    vishmsTimelabel = 1 - vishmsTimelabel;
    set_labels_controls();

  }
  else if(value==MENU_UNITS_SHOWALL){
    show_all_units = 1 - show_all_units;
  }
  updatemenu=1;
  glutPostRedisplay();
}

/* ------------------ OptionMenu ------------------------ */

void OptionMenu(int value){
  if(value==MENU_DUMMY)return;
  updatemenu=1;
  glutPostRedisplay();
  if(value == MENU_OPTION_TRAINERMENU){
    trainer_mode=1;
    if(showtrainer_dialog==0){
      show_glui_trainer();
    }
    FontMenu(LARGE_FONT);
  }
}

/* ------------------ GetNextViewLabel ------------------------ */

void GetNextViewLabel(char *label){
  cameradata *ca;
  int i;

  for(i=1;;i++){
    char view[256],*newview;

    sprintf(view,"view %i",i);
    newview=NULL;
    for(ca=camera_list_first.next;ca->next!=NULL;ca=ca->next){
      if(strcmp(view,ca->name)==0){
        newview=ca->name;
        break;
      }
    }
    if(newview==NULL){
      strcpy(label,view);
      return;
    }
  }
}

/* ------------------ ResetMenu ------------------------ */

void ResetMenu(int value){
  char view_label[256];

  if(value==MENU_DUMMY)return;
  switch(value){
  case MENU_SIZEPRESERVING:
    projection_type = 1 - projection_type;
    Motion_CB(PROJECTION);
    break;
  case MENU_OUTLINEVIEW:
    if(visBlocks==visBLOCKOutline){
      BlockageMenu(visBLOCKAsInput);
    }
    else{
      BlockageMenu(visBLOCKOutline);
    }
    break;
  case MENU_TIMEVIEW:
    UpdateTimes();
    break;
  case SAVE_VIEWPOINT:
    GetNextViewLabel(view_label);
    add_list_view(view_label);
    break;
  case MENU_STARTUPVIEW:
    if(selected_view==MENU_DUMMY)ResetMenu(SAVE_VIEWPOINT);
    set_startup_view();
    break;
  default:
    ASSERT(value>=0);
    if(value<100000){
      reset_glui_view(value);
      if(scriptoutstream!=NULL){
        fprintf(scriptoutstream,"SETVIEWPOINT\n");
        fprintf(scriptoutstream," %s\n",camera_current->name);
      }
    }
    break;
  }
  //updatezoommenu=1; // updating zoom causes a bug when restoring views from the menu
                      // kept commented code in for future reference
  updatemenu=1;
  glutPostRedisplay();
}

/* ------------------ RenderState ------------------------ */

void RenderState(int onoff){
  rendering_status=onoff;
  if(onoff==RENDER_ON){
    update_screeninfo = 1;
    saveW=screenWidth;
    saveH=screenHeight;
    if(renderW==0||renderH==0){
      ResizeWindow(screenWidth,screenHeight);
    }
    else{
      if(renderW>max_screenWidth){
        ResizeWindow(max_screenWidth,max_screenHeight);
      }
      else{
        ResizeWindow(renderW,renderH);
      }
    }
  }
  else{
    Enable360Zoom();
    render_mode = RENDER_XYSINGLE;
    setScreenSize(&saveW,&saveH);
    ResizeWindow(screenWidth,screenHeight);
  }
}

/* ------------------ RenderMenu ------------------------ */

void RenderMenu(int value){
  slicedata *sd;
  int i,n;
  meshdata *meshi;

  updatemenu=1;
  if(value>=11000)return;
  if(opengldefined==1){
    glutPostRedisplay();
  }
  if(value>=10000&&value<=10005){
    nrender_rows=value-10000;
    update_nrender_rows();
    return;
  }
  switch(value){
  case RenderCustom:
    render_window_size = value;
    renderW = script_render_width;
    renderH = script_render_height;
    render_size_index = value;
    break;
  case Render320:
    render_window_size=value;
    renderW=320;
    renderH=240;
    render_size_index=value;
    break;
  case Render640:
    render_window_size=value;
    renderW=640;
    renderH=480;
    render_size_index=value;
    break;
  case RenderWindow:
    render_window_size=value;
    renderW=0;
    renderH=0;
    render_size_index=value;
    break;
  case RENDER_CURRENT_SINGLE:
    render_from_menu=1;
    keyboard('r',FROM_SMOKEVIEW);
     break;
  case RENDER_CURRENT_360:
    LabelMenu(MENU_LABEL_HideAll);
    GetViewportInfo();
    RenderMenu(RENDER_CURRENT_SINGLE);
    render_from_menu = 1;
    keyboard('R', FROM_SMOKEVIEW);
    break;
  case RENDER_CURRENT_MULTIPLE:
    if(nrender_rows==1)RenderMenu(RENDER_CURRENT_SINGLE);
    render_from_menu=1;
    keyboard('R',FROM_SMOKEVIEW);
    break;
  case RenderCancel:
    RenderState(RENDER_OFF);
    break;
  case RenderLABELframenumber:
    render_label_type=RENDER_LABEL_FRAMENUM;
    update_glui_filelabel(render_label_type);
    break;
  case RenderLABELtime:
    render_label_type=RENDER_LABEL_TIME;
    update_glui_filelabel(render_label_type);
    break;
  case RenderPNG:
     render_filetype=PNG;
     updatemenu=1;
     break;
  case RenderJPEG:
     render_filetype=JPEG;
     updatemenu=1;
     break;
  default:
    if(RenderTime==0&&touring==0)return;
    if(touring==1){
      rendertourcount=0;
    }
    if(stept==0){
      keyboard('t',FROM_SMOKEVIEW);
    }
    RenderState(RENDER_ON);
    ResetItimes0();
    for(i=0;i<nsliceinfo;i++){
      sd=sliceinfo+i;
      sd->itime=0;
    }
    frame_index=first_frame_index;
    for(i=0;i<nmeshes;i++){
      meshi=meshinfo+i;
      meshi->patch_itime=0;
    }
    UpdateTimeLabels();
    RenderSkip=value;
    FlowDir=1;
    for(n=0;n<nglobal_times;n++){
      render_frame[n]=0;
    }
    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"RENDERALL\n");
      fprintf(scriptoutstream," %i\n",RenderSkip);
      fprintf(scriptoutstream,"\n");
    }
    render_times = RENDER_ALLTIMES;
    break;
  }
  update_nrender_rows();
}

/* ------------------ EvacShowMenu ------------------------ */

void EvacShowMenu(int value){
  partdata *parti;
  int i;

  if(nevac==0)return;
  if(value==MENU_DUMMY)return;
  if(value<0){
    value = -value;
    value--;
    parti = partinfo + value;
    parti->display = 1 - parti->display;
    updatemenu=1;
    glutPostRedisplay();
    plotstate=GetPlotState(DYNAMIC_PLOTS);
    return;
  }
  if(plotstate==DYNAMIC_PLOTS){
    switch(value){
    case 3:
      visEvac=1;
      for(i=0;i<npartinfo;i++){
        parti = partinfo + i;
        if(parti->loaded==0||parti->evac==0)continue;
        parti->display=1;
      }
      break;
    case 4:
      visEvac=0;
      for(i=0;i<npartinfo;i++){
        parti = partinfo + i;
        if(parti->loaded==0||parti->evac==0)continue;
        parti->display=0;
      }
      break;
      default:
        ASSERT(FFALSE);
        break;
    }
  }
  else{
    switch(value){
    case 3:
      visEvac=1;
      for(i=0;i<npartinfo;i++){
        parti = partinfo + i;
        if(parti->loaded==0||parti->evac==0)continue;
        parti->display=1;
      }
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
  }
  updatemenu=1;
  plotstate=GetPlotState(DYNAMIC_PLOTS);
  glutPostRedisplay();

}

/* ------------------ ParticleShowMenu ------------------------ */

void ParticleShowMenu(int value){
  partdata *parti;
  int i;

  if(npartinfo==0)return;
  if(value==MENU_DUMMY)return;
  if(value<0){
    value = -value;
    value--;
    parti = partinfo + value;
    parti->display = 1 - parti->display;
    updatemenu=1;
    glutPostRedisplay();
    plotstate=GetPlotState(DYNAMIC_PLOTS);
    return;
  }
  if(plotstate==DYNAMIC_PLOTS){
    switch(value){
      case MENU_PARTSHOW_PARTICLES:
        if(visSmokePart==2){
          visSmokePart=0;
        }
        else{
          visSmokePart=2;
        }
        break;
      case MENU_PARTSHOW_DROPLETS:
        visSprinkPart = 1 - visSprinkPart;
        break;
      case MENU_PARTSHOW_SHOWALL:
        visSprinkPart=1;
        visSmokePart=2;
        for(i=0;i<npartinfo;i++){
          parti = partinfo + i;
          if(parti->loaded==0||parti->evac==1)continue;
          parti->display=1;
        }
        break;
      case MENU_PARTSHOW_STATIC:
        break;
      case MENU_PARTSHOW_HIDEALL:
        visSprinkPart=0;
        visSmokePart=0;
        for(i=0;i<npartinfo;i++){
          parti = partinfo + i;
          if(parti->loaded==0||parti->evac==1)continue;
          parti->display=0;
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
    /*
    for(i=0;i<npartinfo;i++){
      parti = partinfo + i;
      if(parti->loaded==0||parti->evac==1)continue;
      parti->display_droplet=0;
      parti->display_smoke=0;
      if(visSmokePart!=0)parti->display_smoke=1;
      if(visSprinkPart==1)parti->display_droplet=1;
    }
    */
    if(visSprinkPart==1||visSmokePart!=0){
      visParticles=1;
    }
    else{
      visParticles=0;
    }
  }
  else{
    switch(value){
      case 1:
        visSmokePart = 2;
        break;
      case 2:
        visSprinkPart = 1;
        break;
      case 3:
        visSprinkPart=1;
        visSmokePart=2;
        for(i=0;i<npartinfo;i++){
          parti = partinfo + i;
          if(parti->loaded==0)continue;
          parti->display=1;
        }
        break;
      case 5:
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
    /*
    for(i=0;i<npartinfo;i++){
      parti = partinfo + i;
      if(parti->loaded==0||parti->evac==1)continue;
      parti->display_droplet=0;
      parti->display_smoke=0;
      if(visSmokePart!=0)parti->display_smoke=1;
      if(visSprinkPart==1)parti->display_droplet=1;
    }
    */
    if(visSmokePart!=0||visSprinkPart==1){
      visParticles=1;
    }
  }
  updatemenu=1;
  plotstate=GetPlotState(DYNAMIC_PLOTS);
  glutPostRedisplay();
}

#define MENU_FRAMERATE_Realtime 2001
#define MENU_FRAMERATE_2xRealtime 2002
#define MENU_FRAMERATE_4xRealtime 2004

/* ------------------ FrameRateMenu ------------------------ */

void FrameRateMenu(int value){
  updateUpdateFrameRateMenu=0;
  realtime_flag=0;
  frameinterval=1;
  if(value > 0){
    switch(value){
    case MENU_FRAMERATE_Realtime:
      if(nglobal_times>0){
        if(global_times!=NULL)frameinterval=1000.*(global_times[nglobal_times-1]-global_times[0])/nglobal_times;
      }
      realtime_flag=1;
      break;
    case MENU_FRAMERATE_2xRealtime:
      if(global_times!=NULL)frameinterval=1000.*(global_times[nglobal_times-1]-global_times[0])/nglobal_times;
      frameinterval /= 2.0;
      realtime_flag=2;
      break;
    case MENU_FRAMERATE_4xRealtime:
      if(global_times!=NULL)frameinterval=1000.*(global_times[nglobal_times-1]-global_times[0])/nglobal_times;
      frameinterval /= 4.0;
      realtime_flag=4;
      break;
    default:
      frameinterval = 1000./value;
      if(frameinterval<1.0){frameinterval = 0.0;}
      break;
    }
    if(global_times==NULL&&realtime_flag!=0)updateUpdateFrameRateMenu=1;
  }
  else{
    keyboard('t',FROM_SMOKEVIEW);
    RenderState(RENDER_OFF);
    FlowDir=1;
  }
  frameratevalue=value;
  updatemenu=1;
  if(opengldefined==1){
    glutPostRedisplay();
  }
  reset_gltime();
}

/* ------------------ IsoSurfaceTypeMenu ------------------------ */

void IsoSurfaceTypeMenu(int value){
  if(ReadPlot3dFile==1){
    switch(value){
    case MENU_SURFACE_SMOOTH:
      p3dsurfacesmooth=1;
      p3dsurfacetype=SURFACE_SOLID;
      break;
    case MENU_SURFACE_FACET:
      p3dsurfacesmooth=0;
      p3dsurfacetype=SURFACE_SOLID;
      break;
    case MENU_SURFACE_OUTLINE:
      p3dsurfacetype=SURFACE_OUTLINE;
      break;
    case MENU_SURFACE_POINTS:
      p3dsurfacetype=SURFACE_POINTS;
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
    update_glui_plot3dtype();
    updatemenu=1;
    glutPostRedisplay();
  }
}

/* ------------------ IsoSurfaceMenu ------------------------ */

void IsoSurfaceMenu(int value){
  if(ReadPlot3dFile==1){
    updatemenu=1;
    glutPostRedisplay();
    if(value==1){
      updateshowstep(0,ISO);
    }
    if(value==2){
      p3dsurfacesmooth = 1 - p3dsurfacesmooth;
    }
  }
}

/* ------------------ LevelMenu ------------------------ */

void LevelMenu(int value){
  if(ReadPlot3dFile==1){
    plotiso[plotn-1]=value;
    updateshowstep(1,ISO);
    updatesurface();
    updatemenu=1;
    glutPostRedisplay();
  }
}

/* ------------------ HelpMenu ------------------------ */

void HelpMenu(int value){
  switch(value){
    case MENU_HELP_ISSUES:
#ifdef pp_OSX
      system("open https://github.com/firemodels/fds-smv/issues");
#endif
#ifdef WIN32
      ShellExecute(NULL, "open", "https://github.com/firemodels/fds-smv/issues", NULL, NULL, SW_SHOWNORMAL);
#endif
      break;
    case MENU_HELP_DOWNLOADS:
#ifdef pp_OSX
      system("open https://pages.nist.gov/fds-smv/downloads.html");
#endif
#ifdef WIN32
      ShellExecute(NULL, "open", "https://pages.nist.gov/fds-smv/downloads.html", NULL, NULL, SW_SHOWNORMAL);
#endif
      break;
    case MENU_HELP_DOCUMENTATION:
#ifdef pp_OSX
      system("open https://pages.nist.gov/fds-smv/");
#endif
#ifdef WIN32
      ShellExecute(NULL, "open", "https://pages.nist.gov/fds-smv/", NULL, NULL, SW_SHOWNORMAL);
#endif
      break;
    case MENU_HELP_FDSWEB:
#ifdef pp_OSX
      system("open https://pages.nist.gov/fds-smv/");
#endif
#ifdef WIN32
      ShellExecute(NULL, "open", "https://pages.nist.gov/fds-smv/", NULL, NULL, SW_SHOWNORMAL);
#endif
      break;
    case MENU_DUMMY:
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
}

/* ------------------ VectorSkipMenu ------------------------ */

void VectorSkipMenu(int value){
  if(value==-1)return; /* dummy label in menu */
  if(value==MENU_VECTOR_SHOW){       /* toggle vector visibility */
    visVector=1-visVector;
    if(vectorspresent==0)visVector=0;
    updatemenu=1;
    glutPostRedisplay();
    return;
  }
  vectorskip=value;
  visVector=1;
  updatemenu=1;
  glutPostRedisplay();
}

/* ------------------ TextureShowMenu ------------------------ */

void TextureShowMenu(int value){
  texturedata *texti;
  int i;
  int texturedisplay=0;
  int texture_flag=0;

  updatefacelists=1;
  if(value>=0){
    texti = textureinfo + value;
    texti->display = 1-texti->display;
    if(texti->display==1)texturedisplay=1;
    for(i=0;i<ntextures;i++){
      texti = textureinfo + i;
      if(texti->loaded==0||texti->used==0)continue;
      if(texti->display==0){
        showall_textures=0;
        texture_flag=1;
        break;
      }
    }
    if(texture_flag==0)showall_textures=1;
  }
  else{
    switch(value){
    case MENU_TEXTURE_SHOWALL:
      for(i=0;i<ntextures;i++){
        texti = textureinfo + i;
        if(texti->loaded==0||texti->used==0)continue;
        texti->display=1;
        texturedisplay=1;
      }
      showall_textures=1;
      break;
    case MENU_TEXTURE_HIDEALL:
      for(i=0;i<ntextures;i++){
        texti = textureinfo + i;
        if(texti->loaded==0||texti->used==0)continue;
        texti->display=0;
      }
      showall_textures=0;
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
  }
  visGeomTextures=0;
  if(texturedisplay==1){
    for(i=0;i<ngeominfo;i++){
      geomdata *geomi;
      surfdata *surf;
      texturedata *textii=NULL;

      geomi = geominfo + i;
      surf = geomi->surf;
      if(surf!=NULL)textii=surf->textureinfo;
      if(textii!=NULL&&textii->display==1){
        visGeomTextures=1;
        break;
      }
    }
    if(value!=visBLOCKOutline&&value!=visBLOCKSolidOutline&&value!=visBLOCKHide){
      BlockageMenu(visBLOCKAsInput);
    }
  }
  updatemenu=1;
  glutPostRedisplay();

}

/* ------------------ Plot3DShowMenu ------------------------ */

void Plot3DShowMenu(int value){
  int i;

  if(value==MENU_PLOT3D_DUMMY)return;
  switch(value){
  case MENU_PLOT3D_Z:
      visz_all=1-visz_all;
      break;
  case MENU_PLOT3D_Y:
      visy_all=1-visy_all;
      break;
  case MENU_PLOT3D_X:
      visx_all=1-visx_all;
      break;
  case MENU_PLOT3D_CONT:
      switch(contour_type){
        case SHADED_CONTOURS:
          contour_type=STEPPED_CONTOURS;
          break;
        case STEPPED_CONTOURS:
        case LINE_CONTOURS:
          contour_type=SHADED_CONTOURS;
          break;
        default:
          ASSERT(FFALSE);
          break;
      }
      break;
  case MENU_PLOT3D_SHOWALL:
      visx_all=1;
      visy_all=1;
      visz_all=1;
      break;
  case MENU_PLOT3D_HIDEALL:
      visx_all=0;
      visy_all=0;
      visz_all=0;
      plotstate=DYNAMIC_PLOTS;
      break;
   case HIDEALL_PLOT3D:
     for(i=0;i<nplot3dinfo;i++){
       if(plot3dinfo[i].loaded==1)plot3dinfo[i].display=0;
     }
     break;
   case SHOWALL_PLOT3D:
     for(i=0;i<nplot3dinfo;i++){
       if(plot3dinfo[i].loaded==1)plot3dinfo[i].display=1;
     }
     break;
   default:
     if(value>=1000){
       if(plotstate==STATIC_PLOTS){
         plot3dinfo[value-1000].display=1-plot3dinfo[value-1000].display;
       }
       else{
         plot3dinfo[value-1000].display=1;
       }
     }
     break;
  }
  plotstate=GetPlotState(STATIC_PLOTS);
  if(plotstate==STATIC_PLOTS&&visiso==1){
    updatesurface();
  }
  updatemenu=1;
  glutPostRedisplay();
}


/* ------------------ GridSliceMenu ------------------------ */

void GridSliceMenu(int value){
  switch(value){
  case GRID_xy:
    visz_all=1-visz_all;
    if(visz_all==1&&visGrid==0)visGrid=1;
    break;
  case GRID_xz:
    visy_all=1-visy_all;
    if(visy_all==1&&visGrid==0)visGrid=1;
    break;
  case GRID_yz:
    visx_all=1-visx_all;
    if(visx_all==1&&visGrid==0)visGrid=1;
    break;
  case GRID_showall:
    visx_all=1;
    visy_all=1;
    visz_all=1;
    visGrid=1;
    break;
  case GRID_hideall:
    visx_all=0;
    visy_all=0;
    visz_all=0;
    break;
  case MENU_DUMMY:
    break;
  case GRID_grid:
    switch(visGrid){
      case GridProbe:
        visGrid=noGridProbe;
        break;
      case GridnoProbe:
        visGrid=noGridnoProbe;
        break;
      case noGridProbe:
        visGrid=GridProbe;
        break;
      case noGridnoProbe:
        visGrid=GridnoProbe;
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
    break;
  case GRID_probe:
    switch(visGrid){
      case GridProbe:
        visGrid=GridnoProbe;
        break;
      case GridnoProbe:
        visGrid=GridProbe;
        break;
      case noGridProbe:
        visGrid=noGridnoProbe;
        break;
      case noGridnoProbe:
        visGrid=noGridProbe;
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  updatemenu=1;
  glutPostRedisplay();
}

#ifdef pp_COMPRESS

/* ------------------ CompressMenu ------------------------ */

void CompressMenu(int value){
  if(value==MENU_DUMMY)return;
  switch(value){
  case MENU_ERASECOMPRESS:
    erase_all=1;
    overwrite_all=0;
    update_overwrite();
    compress_svzip();
    break;
  case MENU_OVERWRITECOMPRESS:
    erase_all=0;
    overwrite_all=1-overwrite_all;
    update_overwrite();
    break;
  case MENU_COMPRESSNOW:
    erase_all=0;
    compress_svzip();
    break;
  case MENU_COMPRESSAUTOLOAD:
    compress_autoloaded=1-compress_autoloaded;
    update_overwrite();
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  updatemenu=1;
}
#endif


/* ------------------ IniSubMenu ------------------------ */

void IniSubMenu(int value){
  if(value==MENU_READCASEINI){
    ReadINI(NULL);
  }
  else{
    char *ini_filename;
    char *script_filename2;

    ini_filename = get_inifilename(value);
    if(ini_filename==NULL||strlen(ini_filename)==0)return;
    script_filename2=script_filename;
    strcpy(script_filename,ini_filename);
    windowresized=0;
    ReadINI(script_filename2);
  }
}

/* ------------------ SmokeviewIniMenu ------------------------ */

void SmokeviewIniMenu(int value){
  switch(value){
  case MENU_READINI:
    ReadINI(NULL);
    UpdateRGBColors(COLORBAR_INDEX_NONE);
    break;
  case MENU_WRITEINI:
    WriteINI(GLOBAL_INI,NULL);
    break;
  case MENU_WRITECASEINI:
    WriteINI(LOCAL_INI,NULL);
    break;
  case MENU_READSVO:
    init_object_defs();
    break;
  case MENU_DUMMY:
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  updatemenu=1;
  glutPostRedisplay();
}

/* ------------------ PeriodicReloads ------------------------ */

void PeriodicReloads(int value){
  if(periodic_reloads!=0){
    LoadUnloadMenu(RELOADALL);
    glutTimerFunc((unsigned int)value,PeriodicReloads,value);
  }
}


/* ------------------ ScriptMenu2 ------------------------ */

void ScriptMenu2(int value){
  script_step=1;
  UpdateScriptStep();
  ScriptMenu(value);
}

/* ------------------ ScriptMenu ------------------------ */

void ScriptMenu(int value){
  int error_code;
  scriptfiledata *scriptfile;
  char newscriptfilename[1024];

  if(value==MENU_DUMMY)return;
  updatemenu=1;
  glutPostRedisplay();
  switch(value){
    case SCRIPT_FILE_LOADING:
      defer_file_loading = 1 - defer_file_loading;
      update_defer();
      break;
    case SCRIPT_STEP:
      script_step=1-script_step;
      break;
    case SCRIPT_CANCEL:
      current_script_command=NULL;
      runscript=0;
      first_frame_index=0;
      skip_render_frames=0;
      script_startframe=-1;
      script_skipframe=-1;
      script_step=0;
      glui_script_enable();
      break;
    case SCRIPT_CONTINUE:
      script_step=0;
      break;
    case SCRIPT_START_RECORDING2:
      defer_file_loading = 1;
      update_defer();
      ScriptMenu(SCRIPT_START_RECORDING);
      break;
    case SCRIPT_START_RECORDING:
      update_script_start();
      get_newscriptfilename(newscriptfilename);
      script_recording = insert_scriptfile(newscriptfilename);
      scriptoutstream=fopen(newscriptfilename,"w");
      if(scriptoutstream!=NULL){
        PRINTF("Script recorder on\n");
        script_recording->recording=1;
        {
          char *renderdir;

          TrimBack(script_renderdir);
          renderdir = TrimFront(script_renderdir);
          if(strlen(renderdir)>0&&strcmp(renderdir,".")!=0){
            fprintf(scriptoutstream,"RENDERDIR\n");
            fprintf(scriptoutstream," %s\n",renderdir);
          }
          else{
            fprintf(scriptoutstream,"RENDERDIR\n");
            fprintf(scriptoutstream," .\n");
          }
        }
        fprintf(scriptoutstream,"XSCENECLIP\n");
        fprintf(scriptoutstream," %i %f %i %f\n",clipinfo.clip_xmin,clipinfo.xmin,clipinfo.clip_xmax,clipinfo.xmax);
        fprintf(scriptoutstream,"YSCENECLIP\n");
        fprintf(scriptoutstream," %i %f %i %f\n",clipinfo.clip_ymin,clipinfo.ymin,clipinfo.clip_ymax,clipinfo.ymax);
        fprintf(scriptoutstream,"ZSCENECLIP\n");
        fprintf(scriptoutstream," %i %f %i %f\n",clipinfo.clip_zmin,clipinfo.zmin,clipinfo.clip_zmax,clipinfo.zmax);
        fprintf(scriptoutstream,"SCENECLIP\n");
        fprintf(scriptoutstream," %i\n",clip_mode);
      }
      else{
        script_recording->recording=0;
        script_recording=NULL;
        fprintf(stderr,"*** Error: The script file %s could not be opened for writing.",newscriptfilename);
      }
      break;
    case SCRIPT_STOP_RECORDING:
      if(script_recording!=NULL){
        script_recording->recording=0;
        AddScriptList(script_recording->file,script_recording->id);
        script_recording=NULL;
      }
      if(scriptoutstream!=NULL){
        fclose(scriptoutstream);
        scriptoutstream=NULL;
        PRINTF("Script recorder off\n");
      }
      update_script_stop();
      break;
    default:
      for(scriptfile=first_scriptfile.next;scriptfile->next!=NULL;scriptfile=scriptfile->next){
        char *file;

        file=scriptfile->file;
        if(file==NULL)continue;
        if(scriptfile->id!=value)continue;
        error_code=compile_script(file);
        if(error_code==0){
      //    ReadINI(NULL);
          start_script();
        }
        else{
          fprintf(stderr,"*** Error (fatal): unable to open script file");
          if(file!=NULL)fprintf(stderr,": %s",file);
          fprintf(stderr,"\n");
          if(from_commandline==1)exit(1);
        }
        break;
      }
      break;
  }
}

/* ------------------ ScriptMenu ------------------------ */
#ifdef pp_LUA
// prototype here rather than create a header file for lua_api
int load_script(char *filename);
void LuaScriptMenu(int value){
  luascriptfiledata *luascriptfile;

  if(value==MENU_DUMMY)return;
  updatemenu=1;
  glutPostRedisplay();
  switch(value){
    // case SCRIPT_FILE_LOADING:
    //   defer_file_loading = 1 - defer_file_loading;
    //   update_defer();
    //   break;
    default:
      for(luascriptfile=first_luascriptfile.next;luascriptfile->next!=NULL;luascriptfile=luascriptfile->next){
        char *file;

        file=luascriptfile->file;
        if(file==NULL)continue;
        if(luascriptfile->id!=value)continue;
        // set the runluascript variable to true, this must be done before
        // calling loadscript, as it both uses and modifies this variable
        runluascript = 1;
        // load the script
        load_script(luascriptfile->file);
        // let the display callback do its work, i.e. just return to the main
        // loop, DisplayCB will run through the script.
        break;
      }
      break;
  }
}
#endif

/* ------------------ ReLoadMenu ------------------------ */

void ReloadMenu(int value){
  int msecs;

  updatemenu=1;
  periodic_value=value;
  switch(value){
  case STOP_RENDERING:
    periodic_reloads=0;
    break;
  case RELOAD_NOW:
    LoadUnloadMenu(RELOADALL);
    break;
  default:
    periodic_reloads=1;
    msecs = value*60*1000;
    glutTimerFunc((unsigned int)msecs,PeriodicReloads,msecs);
    break;
  }
}


/* ------------------ AboutMenu ------------------------ */

void AboutMenu(int value){
}

/* ------------------ LoadUnloadMenu ------------------------ */

void LoadUnloadMenu(int value){
  int errorcode;
  int i;
  int ii;

  if(value==MENU_DUMMY)return;
  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value==UNLOADALL){
   // leaving code here commented in case I later decide to unload terrain files
   // for(i=0;i<nterraininfo;i++){
   //   readterrain("",i,UNLOAD,&errorcode);
   // }
    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"UNLOADALL\n");
    }
    if(hrr_csv_filename!=NULL){
      ReadHRR(UNLOAD, &errorcode);
    }
    if(nvolrenderinfo>0){
      LoadVolSmoke3DMenu(UNLOAD_ALL);
    }
    for(i = 0; i < nsliceinfo; i++){
      slicedata *slicei;

      slicei = sliceinfo + i;
      readslice(slicei->file, i, UNLOAD, DEFER_SLICECOLOR,&errorcode);
    }
    for(i = 0; i<nplot3dinfo; i++){
      readplot3d("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<npatchinfo;i++){
      readpatch(i,UNLOAD,&errorcode);
    }
    for(i=0;i<npartinfo;i++){
      readpart("",i,UNLOAD,PARTDATA,&errorcode);
    }
    for(i=0;i<nisoinfo;i++){
      readiso("",i,UNLOAD,NULL,&errorcode);
    }
    for(i=0;i<nzoneinfo;i++){
      readzone(i,UNLOAD,&errorcode);
    }
    for(i=0;i<nsmoke3dinfo;i++){
      readsmoke3d(i,UNLOAD,&errorcode);
    }
    if(nvolrenderinfo>0){
      UnLoadVolSmoke3DMenu(UNLOAD_ALL);
    }
    updatemenu=1;
    glutPostRedisplay();
  }
  if(value==RELOADALL){
    int last_slice_loaded;

    LOCK_COMPRESS
    ReadSMVDynamic(smv_filename);
    if(hrr_csv_filename!=NULL){
      ReadHRR(LOAD, &errorcode);
    }
    islicetype_save=islicetype;
    for(i=0;i<nsliceinfo;i++){
      sliceinfo[i].reload=1;
    }
    for(i=0;i<nterraininfo;i++){
      if(terraininfo[i].loaded==1){
        readterrain(terraininfo[i].file,i,LOAD,&errorcode);
      }
    }
    for(i=0;i<nvsliceinfo;i++){
      if(vsliceinfo[i].loaded==1){
        readvslice(i,LOAD,&errorcode);
      }
    }
    if(nslice_loaded>1)last_slice_loaded = slice_loaded_list[nslice_loaded-1];
    for(ii = nslice_loaded - 1; ii>=0; ii--){
      slicedata *slicei;


      i = slice_loaded_list[ii];
      slicei = sliceinfo + i;
      if(slicei->reload == 1){
        last_slice_loaded = i;
        break;
      }
    }
    for(ii = 0; ii<nslice_loaded; ii++){
      slicedata *slicei;


      i = slice_loaded_list[ii];
      slicei = sliceinfo + i;
      if(slicei->reload==1){
        int set_slicecolor;

        set_slicecolor = DEFER_SLICECOLOR;
        if(i == last_slice_loaded)set_slicecolor = SET_SLICECOLOR;
        readslice(slicei->file,i,LOAD,set_slicecolor,&errorcode);
      }
    }
    islicetype=islicetype_save;
    for(i=0;i<nplot3dinfo;i++){
      if(plot3dinfo[i].loaded==1){
        readplot3d(plot3dinfo[i].file,i,LOAD,&errorcode);
      }
    }
    for(ii=0;ii<npatch_loaded;ii++){
      i = patch_loaded_list[ii];
      readpatch(i,LOAD,&errorcode);
    }
    for(i=0;i<nsmoke3dinfo;i++){
      if(smoke3dinfo[i].loaded==1){
        readsmoke3d(i,LOAD,&errorcode);
      }
    }
    for(i=0;i<npartinfo;i++){
      if(partinfo[i].loaded==1){
        partinfo[i].reload=1;
        readpart(partinfo[i].file,i,UNLOAD,PARTDATA,&errorcode);
      }
      else{
        partinfo[i].reload=0;
      }
    }
    npartframes_max=get_min_partframes();
    for(i=0;i<npartinfo;i++){
      if(partinfo[i].reload==1){
        readpart(partinfo[i].file, i, UNLOAD, PARTDATA,&errorcode);
      }
    }
    for(i=0;i<npartinfo;i++){
      if(partinfo[i].reload==1){
        readpart(partinfo[i].file, i, LOAD, PARTDATA,&errorcode);
      }
    }
    update_readiso_geom_wrapup = UPDATE_ISO_START_ALL;
    for(i = 0; i<nisoinfo; i++){
      isodata *isoi;

      isoi = isoinfo + i;
      if(isoi->loaded==0)continue;
      readiso(isoi->file,i,LOAD,NULL,&errorcode);
    }
    if(update_readiso_geom_wrapup == UPDATE_ISO_ALL_NOW)readiso_geom_wrapup();
    update_readiso_geom_wrapup = UPDATE_ISO_OFF;
    UNLOCK_COMPRESS
  //  plotstate=DYNAMIC_PLOTS;
  //  visParticles=1;
    updatemenu=1;
    glutPostRedisplay();
  }
  if(value==SHOWFILES){
    glutPostRedisplay();
    showfiles=1-showfiles;
    updatemenu=1;
    update_slice_menulabels();
    update_vslice_menulabels();
    update_smoke3d_menulabels();
    update_patch_menulabels();
    update_iso_menulabels();
    update_part_menulabels();
    update_tour_menulabels();
    update_plot3d_menulabels();
  }
  if(value==REDIRECT){
    updatemenu=1;
    glutPostRedisplay();
    redirect=1-redirect;
    if(LOG_FILENAME!=NULL){
      fclose(LOG_FILENAME);
      LOG_FILENAME=NULL;
    }
    if(redirect==1){
      LOG_FILENAME=fopen(log_filename,"w");
      if(LOG_FILENAME==NULL)redirect=0;
    }
    if(redirect==1){
      set_stdout(LOG_FILENAME);
    }
    else{
      set_stdout(stdout);
    }
  }
  glutSetCursor(GLUT_CURSOR_RIGHT_ARROW);
}

void AvatarEvacMenu(int value){
  if(value==MENU_DUMMY)return;
  iavatar_evac=value;
  updatemenu=1;
}

/* ------------------ TourMenu ------------------------ */

void TourMenu(int value){
  tourdata *touri;
  int i;

  if(value==MENU_DUMMY)return;
  touring=0;
  updatemenu=1;
  glutPostRedisplay();
  switch(value){
  case MENU_TOUR_EDIT:
    DialogMenu(DIALOG_TOUR);
    break;
  case MENU_TOUR_NEW:
    add_new_tour();
    DialogMenu(DIALOG_TOUR);
    break;
  case MENU_TOUR_CLEARALL:
    for(i=0;i<ntours;i++){  // clear all tours
      touri = tourinfo + i;
      touri->display=touri->display2;
    }
    if(viewtourfrompath==1){
      SetViewPoint(RESTORE_EXTERIOR_VIEW);
    }
    from_glui_trainer=0;
    for(i=0;i<ntours;i++){
      touri = tourinfo + i;
      if(touri->display==1){
        selected_tour=touri;
        break;
      }
    }
    selected_tour=NULL;
    break;
  case MENU_TOUR_MANUAL:
    for(i=0;i<ntours;i++){  // clear all tours
      touri = tourinfo + i;
      touri->display=0;
    }
    if(viewtourfrompath==1){
      SetViewPoint(RESTORE_EXTERIOR_VIEW);
    }
    from_glui_trainer=0;
    selected_tour=NULL;
    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"UNLOADTOUR\n");
    }
    DialogMenu(DIALOG_TOUR);
    break;
  case MENU_TOUR_SHOWDIALOG:
    edittour=1-edittour;
    if(edittour==1&&showtour_dialog==0){
      show_glui_tour();
    }
    break;
  case MENU_TOUR_SHOWALL:               // show all tours
    for(i=0;i<ntours;i++){
      touri = tourinfo + i;
      touri->display=1;
    }
    plotstate=GetPlotState(DYNAMIC_PLOTS);
    break;
  case MENU_TOUR_VIEWFROMROUTE:               // view from route
    viewtourfrompath = 1 - viewtourfrompath;
    if(viewtourfrompath==0)SetViewPoint(RESTORE_EXTERIOR_VIEW);
    break;
  case MENU_TOUR_DEFAULT:
    for(i=0;i<ntours;i++){
      touri = tourinfo + i;
      touri->display=0;
    }
    SetViewPoint(RESTORE_EXTERIOR_VIEW);
    defaulttour();
    break;
  default:
    if(value<-22){
      tourlocus_type=2;
      iavatar_types=(-value-23);
      if(selectedtour_index>=0&&selectedtour_index<ntours){
        tourinfo[selectedtour_index].glui_avatar_index=iavatar_types;
      }
    }

    //  show one tour

    if(value>=0&&value<ntours){
      int j;

      touri = tourinfo + value;
      touri->display = 1 - touri->display;
      if(touri->display==1){
        selectedtour_index=value;
        selected_frame=touri->first_frame.next;
        selected_tour=touri;
      }
      else{
        for(j=0;j<ntours;j++){
          tourdata *tourj;

          tourj = tourinfo + j;
          if(touri==tourj||tourj->display==0)continue;
          selectedtour_index=j;
          selected_frame=tourj->first_frame.next;
          selected_tour=tourj;
          break;
        }
      }
    }
    break;
  }
  updateviewtour();
  delete_tourlist();
  create_tourlist();
  update_tourcontrols();
  plotstate=GetPlotState(DYNAMIC_PLOTS);
  if(value!=-5&&value!=-4)UpdateTimes();
  callfrom_tourglui=0;

}

/* ------------------ SetTour ------------------------ */

void SetTour(tourdata *thetour){
  int tournumber;

  if(thetour==NULL)return;
  tournumber = thetour - tourinfo;
  TourMenu(tournumber);
}

/* ------------------ EvacMenu ------------------------ */

void EvacMenu(int value){
  int errorcode;

  if(value==MENU_EVAC_DUMMY)return;
  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value==MENU_EVAC_ALLMESHES){
    int i;

    for(i=0;i<npartinfo;i++){
      partdata *parti;

      parti=partinfo + i;
      if(parti->evac==0)continue;
      readpart(parti->file, i, UNLOAD, PARTDATA,&errorcode);
    }
    npartframes_max=get_min_partframes();
    for(i=0;i<npartinfo;i++){
      partdata *parti;

      parti=partinfo + i;
      if(parti->evac==0)continue;
      ReadEvacFile=1;
      readpart(parti->file, i, LOAD, PARTDATA,&errorcode);
      if(scriptoutstream!=NULL){
        fprintf(scriptoutstream,"LOADFILE\n");
        fprintf(scriptoutstream," %s\n",parti->file);
      }
    }
    force_redisplay=1;
    UpdateFrameNumber(0);
  }
  if(value>=0){
    ReadEvacFile=1;
    npartframes_max=get_min_partframes();
    readpart(partinfo[value].file, value, LOAD, PARTDATA,&errorcode);
    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"LOADFILE\n");
      fprintf(scriptoutstream," %s\n",partinfo[value].file);
    }
  }
  else if(value==MENU_EVAC_UNLOADALL){
    int i;

    for(i=0;i<npartinfo;i++){
      if(partinfo[i].evac==0)continue;
      readpart("", i, UNLOAD, PARTDATA,&errorcode);
    }
  }
  updatemenu=1;
  glutPostRedisplay();
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ UpdateStreakValue ------------------------ */

void UpdateStreakValue(float value){
  partdata *parti=NULL;
  int i;

  streak_index=-1;
  for(i=0;i<nstreak_rvalue;i++){
    if(ABS(value-streak_rvalue[i])<0.01){
      streak_index=i;
      float_streak5value=streak_rvalue[i];
      break;
    }
  }
  for(i=0;i<npartinfo;i++){
    parti = partinfo + i;
    if(parti->loaded==1)break;
  }
  if(parti!=NULL&&parti->loaded==1&&parti->ntimes>1){
    for(i=0;i<parti->ntimes-1;i++){
      if(parti->times[i]<=value&&value<parti->times[i+1]){
        streak5step=i;
        break;
      }
    }
  }
}
/* ------------------ ParticleStreakShowMenu ------------------------ */

void ParticleStreakShowMenu(int value){
  float rvalue;

  if(value==-1)return;
  if(value==MENU_STREAK_HIDE){
    streak5show=0;
    streak5step=0;
  }
  else if(value==MENU_STREAK_HEAD){
    showstreakhead=1-showstreakhead;
  }
  else{
    streak5show=1;
    streak5step=0;
    rvalue=streak_rvalue[value];
    UpdateStreakValue(rvalue-0.001);
    update_glui_streakvalue(rvalue);

  }
  updatemenu=1;
  glutPostRedisplay();
}

/* ------------------ Particle5ShowMenu ------------------------ */

void Particle5ShowMenu(int value){
}

/* ------------------ ParticlePropShowMenu ------------------------ */

void ParticlePropShowMenu(int value){
  partpropdata *propi;

  int propvalue;

  propvalue = (-value)/10000-1;
  value = -((-value)%10000);

  if(value==MENU_PROP_DUMMY)return;

  if(value>=0){
    int iprop;
    int i;

    part5show=1;
    parttype=0;
    iprop = value;
    for(i=0;i<npart5prop;i++){
      propi = part5propinfo + i;
      propi->display=0;
    }

    propi = part5propinfo + iprop;
    last_prop_display=iprop;
    propi->display=1;
    part5colorindex=iprop;

    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"PARTCLASSCOLOR\n");
      fprintf(scriptoutstream," %s\n",propi->label->longlabel);
    }
    current_property = propi;
    if(iprop!=0){
      parttype=1;
    }
    prop_index = iprop;
    partshortlabel=propi->label->shortlabel;
    partunitlabel=propi->label->unit;
    partscale=propi->scale;
  }
  else if(value==MENU_PROP_SHOWALL){
    if(current_property!=NULL){
      unsigned char *vis;
      int i;

      vis = current_property->class_vis;
      for(i=0;i< npartclassinfo;i++){
        vis[i]=1;
      }
    }
  }
  else if(value==MENU_PROP_HIDEALL){
    if(current_property!=NULL){
      unsigned char *vis;
      int i;

      vis = current_property->class_vis;
      for(i=0;i< npartclassinfo;i++){
        vis[i]=0;
      }
    }

  }
  else if(value==MENU_PROP_HIDEPART){
    int i;
    int unhide=1;

    for(i=0;i<npart5prop;i++){
      propi = part5propinfo + i;
      if(propi->particle_property==1){
        if(propi->display==1)unhide=0;
        propi->display=0;
      }
    }
    part5show=0;
    parttype=0;
    if(unhide==1&&last_prop_display>=0){
      ParticlePropShowMenu(last_prop_display);
    }
  }
  else if(value==MENU_PROP_HIDEAVATAR){
    int i;

    for(i=0;i<npart5prop;i++){
      propi = part5propinfo + i;
      if(propi->human_property==1){
        propi->display=0;
      }
    }
    part5show=0;
    parttype=0;
  }
  else if(value==MENU_PROP_TRACERS){
    show_tracers_always=1-show_tracers_always;
    updatetracers();
  }
  else{
    int iclass;
    int vistype;

    iclass =  (-value - 10)/5;
    vistype = (-value - 10)%5;
    if(vistype==0){
      if(current_property!=NULL){
        unsigned char *vis;

        vis = current_property->class_vis;
        vis[iclass] = 1 - vis[iclass];
        if(scriptoutstream!=NULL){
          fprintf(scriptoutstream,"PARTCLASSTYPE\n");
          fprintf(scriptoutstream," %s\n",current_property->label->longlabel);
        }
      }
    }
    else{
      partclassdata *partclassj;

      partclassj = partclassinfo + iclass;
      partclassj->vis_type=vistype;
      PropMenu(propvalue);
    }

  }
  updatemenu=1;
  glutPostRedisplay();
}

/* ------------------ LoadParticleMenu ------------------------ */

void LoadParticleMenu(int value){
  int errorcode,i;
  partdata *parti;

  get_allpart_histogram();
  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value>=0){
    char  *partfile;

    ReadPartFile=1;
    partfile = partinfo[value].file;
    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"LOADFILE\n");
      fprintf(scriptoutstream," %s\n",partfile);
    }
    npartframes_max=get_min_partframes();
    readpart(partfile, value, LOAD, PARTDATA,&errorcode);
  }
  else{
    if(value==-1){
      for(i=0;i<npartinfo;i++){
        if(partinfo[i].evac==1)continue;
        readpart("", i, UNLOAD, PARTDATA,&errorcode);
      }
    }
    else{
      if(scriptoutstream!=NULL){
        fprintf(scriptoutstream,"LOADPARTICLES\n");
      }
      npartframes_max=get_min_partframes();
      if(value==PARTFILE_LOADALL){
        for(i = 0; i<npartinfo; i++){
          parti = partinfo+i;
          if(parti->evac==1)continue;
          readpart(parti->file, i, UNLOAD, PARTDATA, &errorcode);
        }
      }
      for(i=0;i<npartinfo;i++){
        parti = partinfo + i;
        if(parti->evac==1)continue;
        if(parti->loaded==0&&value==PARTFILE_RELOADALL)continue;
        readpart(parti->file, i, LOAD, PARTDATA,&errorcode);
      }
      force_redisplay=1;
      UpdateFrameNumber(0);
    }
  }
  updatemenu=1;
  glutPostRedisplay();
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ ZoneMenu ------------------------ */

void ZoneMenu(int value){
  int i,errorcode;
  if(value>=0){
    if(scriptoutstream!=NULL){
      zonedata *zonei;

      zonei = zoneinfo + value;
      fprintf(scriptoutstream,"LOADFILE\n");
      fprintf(scriptoutstream," %s\n",zonei->file);
    }
    readzone(value,LOAD,&errorcode);
  }
  else{
    for(i=0;i<nzoneinfo;i++){
      readzone(i,UNLOAD,&errorcode);
    }
  }
  updatemenu=1;
  glutPostRedisplay();
}

/* ------------------ UnloadVSliceMenu ------------------------ */

void UnloadVSliceMenu(int value){
  int errorcode,i;

  updatemenu=1;
  glutPostRedisplay();
  if(value>=0){
    readvslice(value,UNLOAD,&errorcode);
  }
  else if(value==UNLOAD_ALL){
    for(i=0;i<nvsliceinfo;i++){
      readvslice(i,UNLOAD,&errorcode);
    }
  }
  else if(value==-2){
    int unload_index;

    unload_index=last_vslice_loadstack();
    if(unload_index>=0&&unload_index<nvsliceinfo){
      readvslice(unload_index,UNLOAD,&errorcode);
    }
  }
}

/* ------------------ UnloadPatchMenu ------------------------ */

void UnloadPatchMenu(int value){
  int errorcode,i;

  updatemenu=1;
  glutPostRedisplay();
  if(value>=0){
    readpatch(value,UNLOAD,&errorcode);
  }
  else{
    for(i=0;i<npatchinfo;i++){
      readpatch(i,UNLOAD,&errorcode);
    }
  }
}

/* ------------------ UnloadIsoMenu ------------------------ */

void UnloadIsoMenu(int value){
  int errorcode,i;

  updatemenu=1;
  glutPostRedisplay();
  if(value>=0){
    readiso("",value,UNLOAD,NULL,&errorcode);
  }
  else{
    for(i=0;i<nisoinfo;i++){
      readiso("",i,UNLOAD,NULL,&errorcode);
    }
  }
}

/* ------------------ UnloadPlot3dMenu ------------------------ */

void UnloadPlot3dMenu(int value){
  int errorcode,i;

  updatemenu=1;
  glutPostRedisplay();
  if(value>=0){
    readplot3d("",value,UNLOAD,&errorcode);
  }
  else{
    for(i=0;i<nplot3dinfo;i++){
      readplot3d("",i,UNLOAD,&errorcode);
    }
  }
}

/* ------------------ UnloadEvacMenu ------------------------ */

void UnloadEvacMenu(int value){
  int errorcode,i;

  updatemenu=1;
  glutPostRedisplay();
  if(value>=0){
    readpart("", value, UNLOAD, PARTDATA,&errorcode);
  }
  else{
    for(i=0;i<npartinfo;i++){
      if(partinfo[i].evac==0)continue;
      readpart("", i, UNLOAD, PARTDATA,&errorcode);
    }
  }
}

/* ------------------ UnloadPartMenu ------------------------ */

void UnloadPartMenu(int value){
  int errorcode,i;

  updatemenu=1;
  glutPostRedisplay();
  if(value>=0){
    readpart("", value, UNLOAD, PARTDATA,&errorcode);
  }
  else{
    for(i=0;i<npartinfo;i++){
      if(partinfo[i].evac==1)continue;
      readpart("", i, UNLOAD, PARTDATA,&errorcode);
    }
  }
}

/* ------------------ LoadVSliceMenu ------------------------ */

void LoadVSliceMenu(int value){
  int errorcode;
  int i;

  if(value==MENU_DUMMY)return;
  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value==UNLOAD_ALL){
    for(i=0;i<nvsliceinfo;i++){
      readvslice(i,UNLOAD,&errorcode);
    }
    return;
  }
  else if(value==MENU_LOADVSLICE_SHOWALL){
    showallslicevectors=1-showallslicevectors;
    updatemenu=1;
    glutPostRedisplay();
  }
  else if(value>=0){
    vslicedata *vslicei;
    slicedata *slicei;

    readvslice(value, LOAD, &errorcode);
    vslicei = vsliceinfo + value;
    slicei = vslicei->val;
    if(script_multivslice==0&&slicei!=NULL&&scriptoutstream!=NULL){
      fprintf(scriptoutstream,"LOADVSLICEM\n");
      fprintf(scriptoutstream," %s\n", slicei->label.longlabel);
      fprintf(scriptoutstream," %i %f\n", slicei->idir, slicei->position_orig);
      fprintf(scriptoutstream," %i\n", slicei->blocknumber+1);
    }
  }
  else{
    int submenutype;
    char *submenulabel;
    vslicedata *vslicei;
    slicedata *slicei;
    int dir;

    value = -(1000 + value);
    submenutype=value/4;
    dir=value%4;
    submenutype=subvslice_menuindex[submenutype];
    vslicei = vsliceinfo + submenutype;
    slicei = sliceinfo + vslicei->ival;
    submenulabel = slicei->label.longlabel;
    for(i=0;i<nvsliceinfo;i++){
      char *longlabel;

      vslicei = vsliceinfo + i;
      slicei=sliceinfo + vslicei->ival;
      longlabel = slicei->label.longlabel;
      if(strcmp(longlabel,submenulabel)!=0)continue;
      if(dir!=0&&dir!=slicei->idir)continue;
      readvslice(i,LOAD,&errorcode);
    }
  }
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ UnloadSliceMenu ------------------------ */

void UnloadSliceMenu(int value){
  int errorcode,i;

  updatemenu=1;
  glutPostRedisplay();
  if(value>=0){
    readslice("",value,UNLOAD,SET_SLICECOLOR,&errorcode);
  }
  else{
    if(value==UNLOAD_ALL){
      for(i=0;i<nsliceinfo;i++){
        readslice("",i,UNLOAD,DEFER_SLICECOLOR,&errorcode);
      }
    }
    else if(value==UNLOAD_LAST){
      int unload_index;

      unload_index=last_slice_loadstack();
      if(unload_index>=0&&unload_index<nsliceinfo){
        readslice("",unload_index,UNLOAD,SET_SLICECOLOR,&errorcode);
      }
    }
  }
}

/* ------------------ UnLoadMultiVSliceMenu ------------------------ */

void UnloadMultiVSliceMenu(int value){
  int i;
  multivslicedata *mvslicei;

  if(value>=0){
    mvslicei = multivsliceinfo + value;
    for(i=0;i<mvslicei->nvslices;i++){
      UnloadSliceMenu(mvslicei->ivslices[i]);
    }
  }
  else{
    UnloadSliceMenu(UNLOAD_ALL);
  }
}

/* ------------------ UnLoadMultiSliceMenu ------------------------ */

void UnloadMultiSliceMenu(int value){
  int i;
  multislicedata *mslicei;

  if(value>=0){
    mslicei = multisliceinfo + value;
    for(i=0;i<mslicei->nslices;i++){
      UnloadSliceMenu(mslicei->islices[i]);
    }
  }
  else{
    UnloadSliceMenu(UNLOAD_ALL);
  }
}

/* ------------------ ShowVolSmoke3DMenu ------------------------ */

void ShowVolSmoke3DMenu(int value){
  int i;

  if(value==MENU_DUMMY)return;
  updatemenu=1;
  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value>=0){
    meshdata *meshi;
    volrenderdata *vr;

    meshi = meshinfo + value;
    vr = &(meshi->volrenderinfo);
    if(vr->fireslice!=NULL||vr->smokeslice!=NULL){
      if(vr->loaded==1){
        vr->display=1-vr->display;
        PRINTF("%s vis state:%i\n",meshi->label,vr->display);
      }
    }
  }
  else if(value==HIDE_ALL){  // hide all
    for(i=0;i<nmeshes;i++){
      meshdata *meshi;
      volrenderdata *vr;

      meshi = meshinfo + i;
      vr = &(meshi->volrenderinfo);
      if(vr->fireslice==NULL||vr->smokeslice==NULL)continue;
      if(vr->loaded==1){
        vr->display=0;
        PRINTF("%s vis state:%i\n",meshi->label,vr->display);
      }
    }
  }
  else if(value==SHOW_ALL){  // show all
    for(i=0;i<nmeshes;i++){
      meshdata *meshi;
      volrenderdata *vr;

      meshi = meshinfo + i;
      vr = &(meshi->volrenderinfo);
      if(vr->fireslice==NULL||vr->smokeslice==NULL)continue;
      if(vr->loaded==1){
        vr->display=1;
        PRINTF("%s vis state:%i\n",meshi->label,vr->display);
      }
    }
  }
  updatemenu=1;
  glutPostRedisplay();
}

/* ------------------ UnLoadVolSmoke3DMenu ------------------------ */

void UnLoadVolSmoke3DMenu(int value){
  int i;

  if(value==MENU_DUMMY)return;
  read_vol_mesh=VOL_UNLOAD;
  updatemenu=1;
  if(value<0){
    if(value==UNLOAD_ALL){
      for(i=0;i<nmeshes;i++){
        meshdata *meshi;
        volrenderdata *vr;

        meshi = meshinfo + i;
        vr = &(meshi->volrenderinfo);
        if(vr->fireslice==NULL||vr->smokeslice==NULL)continue;
        if(vr->loaded==1){
          unload_volsmoke_allframes(vr);
        }
      }
    }
  }
  else{
    meshdata *meshi;
    volrenderdata *vr;
    slicedata *fireslice, *smokeslice;

    meshi = meshinfo + value;
    vr = &(meshi->volrenderinfo);
    fireslice = vr->fireslice;
    smokeslice = vr->smokeslice;
    if(fireslice!=NULL||smokeslice!=NULL){
      unload_volsmoke_allframes(vr);
    }
  }
  updatemenu=1;
  read_vol_mesh=VOL_READNONE;
  glutPostRedisplay();
}

/* ------------------ LoadVolSmoke3DMenu ------------------------ */

void LoadVolSmoke3DMenu(int value){
  if(value==MENU_DUMMY)return;
  updatemenu=1;
  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value>=0){
    meshdata *meshi;
    volrenderdata *vr;
    slicedata *fireslice, *smokeslice;

    update_smokecolorbar=1;
    meshi = meshinfo + value;
    vr = &(meshi->volrenderinfo);
    fireslice = vr->fireslice;
    smokeslice = vr->smokeslice;
    if(smokeslice!=NULL&&fireslice!=NULL){
      if(scriptoutstream!=NULL){
        fprintf(scriptoutstream,"LOADVOLSMOKE\n");
        fprintf(scriptoutstream," %i\n",value);
      }
      if(read_vol_mesh==VOL_READNONE){
        read_vol_mesh=value;
        read_volsmoke_allframes_allmeshes();
      }
      else{
        fprintf(stderr,"*** Warning: 3D smoke is currently being loaded\n");
        fprintf(stderr,"   Load data when this is complete.\n");
      }
    }
  }
  else if(value==UNLOAD_ALL){  // unload all
    if(read_vol_mesh==VOL_READNONE){
      UnLoadVolSmoke3DMenu(value);
    }
      else{
        if(read_vol_mesh==VOL_UNLOAD){
          fprintf(stderr,"*** Warning: data is currently being unloaded\n");
        }
        else{
          fprintf(stderr,"*** Warning: data is currently being loaded\n");
        }
        fprintf(stderr,"    Continue when this is complete.\n");
      }
  }
  else if(value==LOAD_ALL){  // load all
    update_smokecolorbar=1;
    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"LOADVOLSMOKE\n");
      fprintf(scriptoutstream," -1\n");
    }
    if(read_vol_mesh==VOL_READNONE){
      read_vol_mesh=VOL_READALL;
      read_volsmoke_allframes_allmeshes();
    }
    else{
      if(read_vol_mesh==VOL_UNLOAD){
        fprintf(stderr,"*** Warning: data is currently being unloaded\n");
      }
      else{
        fprintf(stderr,"*** Warning: data is currently being loaded\n");
      }
      fprintf(stderr,"    Continue when this is complete.\n");
    }
  }
  updatemenu=1;
  Idle_CB();
  glutPostRedisplay();
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ UnLoadSmoke3DMenu ------------------------ */

void UnLoadSmoke3DMenu(int value){
  int errorcode;
  int i;
  smoke3ddata *smoke3di;

  if(value==MENU_DUMMY)return;
  updatemenu=1;
  if(value<0){
    value= -value;
    for(i=0;i<nsmoke3dinfo;i++){
      smoke3di = smoke3dinfo + i;
      if(smoke3di->loaded==1&&smoke3di->type==value){
        readsmoke3d(i,UNLOAD,&errorcode);
      }
    }
  }
  else{
    readsmoke3d(value,UNLOAD,&errorcode);
  }
}

/* ------------------ LoadSmoke3DMenu ------------------------ */

void LoadSmoke3DMenu(int value){
  int i,errorcode;
  smoke3ddata *smoke3di, *smoke3dj;

  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value>=0){
    if(scriptoutstream!=NULL){
      char *file;

      file = smoke3dinfo[value].file;
      fprintf(scriptoutstream,"LOADFILE\n");
      fprintf(scriptoutstream," %s\n",file);
    }
    if(scriptoutstream==NULL||defer_file_loading==0){
      readsmoke3d(value,LOAD,&errorcode);
    }
  }
  else if(value==UNLOAD_ALL){
    for(i=0;i<nsmoke3dinfo;i++){
      readsmoke3d(i,UNLOAD,&errorcode);
    }
  }
  else if(value==MENU_SMOKE3D_IBLANK){
    update_makeiblank_smoke3d = 1;
  }
  else if(value==-9){
    if(scriptoutstream==NULL||defer_file_loading==0){
      for(i=0;i<nsmoke3dinfo;i++){
        smoke3di = smoke3dinfo + i;
        if(smoke3di->loaded==1)continue;
        readsmoke3d(i,LOAD,&errorcode);
      }
    }
    ASSERT(FFALSE); // check to see if this code segment is used
  }
  else if(value<=-10){
    value = -(value + 10);
    smoke3dj = smoke3dinfo + value;
    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"LOAD3DSMOKE\n");
      fprintf(scriptoutstream," %s\n",smoke3dj->label.longlabel);
    }
    if(scriptoutstream==NULL||defer_file_loading==0){
      for(i=0;i<nsmoke3dinfo;i++){
        smoke3di = smoke3dinfo + i;
        if(strcmp(smoke3di->label.shortlabel,smoke3dj->label.shortlabel)==0){
          readsmoke3d(i,LOAD,&errorcode);
        }
      }
    }
  }
  updatemenu=1;
  glutPostRedisplay();
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ AnySmoke ------------------------ */

int AnySmoke(char *type){

  if(nsmoke3dinfo>0)return 1;
  return 0;
}

/* ------------------ AnySlices ------------------------ */

int AnySlices(char *type){
  int i;

  for(i=0;i<nsliceinfo;i++){
    if(STRCMP(sliceinfo[i].label.longlabel,type)==0)return 1;
  }
  return 0;
}

/* ------------------ UnLoadTerrainMenu ------------------------ */

void UnloadTerrainMenu(int value){
  int i;
  int errorcode;

  if(value >= 0 && value < nterraininfo){
    readterrain("", value, UNLOAD, &errorcode);
  }
  else if(value == MENU_UNLOADTERRAIN_UNLOADALL){
    for(i = 0; i < nterraininfo; i++){
      UnloadTerrainMenu(i);
    }
  }
  updatemenu = 1;
  glutPostRedisplay();

}

/* ------------------ LoadTerrainMenu ------------------------ */

void LoadTerrainMenu(int value){
  int i;
  int errorcode;

  if(value>=0&&value<nterraininfo){
    terraindata *terri;

    terri = terraininfo + value;
    readterrain(terri->file,value,LOAD,&errorcode);
  }
  else if(value==MENU_LOADTERRAIN_UNLOAD){
    UnloadTerrainMenu(value);
  }
  else if(value==MENU_LOADTERRAIN_LOADALL){
    for(i=0;i<nterraininfo;i++){
      LoadTerrainMenu(i);
    }
  }
  updatemenu=1;
  glutPostRedisplay();
}

/* ------------------ DefineAllFEDs ------------------------ */

void DefineAllFEDs(void){
  int i;

  compute_fed=0;
  for(i=nsliceinfo-nfedinfo;i<nsliceinfo;i++){
    LoadSliceMenu(i);
    UnloadSliceMenu(i);
  }
  exit(0);
}

/* ------------------ LoadSlicei ------------------------ */

void LoadSlicei(int set_slicecolor, int value){
  slicedata *slicei;
  int errorcode;

  slicei = sliceinfo + value;
  slicei->loading=1;
  if(script_multislice == 0 && scriptoutstream != NULL){
    fprintf(scriptoutstream, "LOADSLICEM\n");
    fprintf(scriptoutstream, " %s\n", slicei->label.longlabel);
    fprintf(scriptoutstream, " %i %f\n", slicei->idir, slicei->position_orig);
    fprintf(scriptoutstream, " %i\n", slicei->blocknumber + 1);
  }
  if(scriptoutstream == NULL || defer_file_loading == 0){
    if(value < nsliceinfo - nfedinfo){
      readslice(slicei->file, value, LOAD, set_slicecolor, &errorcode);
    }
    else{
      readfed(value, LOAD, FED_SLICE, &errorcode);
    }
  }
  slicei->loading=0;
}

/* ------------------ LoadSliceMenu ------------------------ */

void LoadSliceMenu(int value){
  int errorcode,i;

  if(value==MENU_DUMMY)return;
  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value>=0){
    LoadSlicei(SET_SLICECOLOR,value);
  }
  else{
    if(value==UNLOAD_ALL){
      for(i=0;i<nsliceinfo;i++){
        readslice("",i,UNLOAD,DEFER_SLICECOLOR,&errorcode);
      }
    }
    else if(value==MENU_SHOWSLICE_INBLOCKAGE){
      show_slice_in_obst=1-show_slice_in_obst;
      update_show_slice_in_obst();
    }
    else{
      int submenutype;
      char *submenulabel;
      slicedata *slicei;
      int dir;
      int last_slice;

      value = -(1000 + value);
      submenutype=value/4;
      dir=value%4;
      submenutype=subslice_menuindex[submenutype];
      slicei = sliceinfo + submenutype;
      submenulabel = slicei->label.longlabel;
      last_slice = nsliceinfo - 1;
      for(i = nsliceinfo-1; i>=0; i--){
        char *longlabel;

        slicei = sliceinfo + i;
        longlabel = slicei->label.longlabel;
        if(strcmp(longlabel, submenulabel) != 0)continue;
        if(dir != 0 && dir != slicei->idir)continue;
        last_slice = i;
        break;
      }
      for(i = 0; i<nsliceinfo; i++){
        char *longlabel;
        int set_slicecolor;

        slicei = sliceinfo + i;
        longlabel = slicei->label.longlabel;
        if(strcmp(longlabel,submenulabel)!=0)continue;
        if(dir!=0&&dir!=slicei->idir)continue;
        set_slicecolor = DEFER_SLICECOLOR;
        if(i == last_slice)set_slicecolor = SET_SLICECOLOR;
        readslice(slicei->file,i,LOAD,set_slicecolor,&errorcode);
      }
    }
  }
  updatemenu=1;
  glutPostRedisplay();
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ LoadVMultiSliceMenu ------------------------ */

void LoadMultiVSliceMenu(int value){
  int i;
  multivslicedata *mvslicei;

  if(value==MENU_DUMMY)return;
  if(value>=0){
    mvslicei = multivsliceinfo + value;
    if(scriptoutstream!=NULL){
      if(mvslicei->nvslices>0){
        slicedata *slicei;

        slicei = sliceinfo + mvslicei->ivslices[0];
        fprintf(scriptoutstream,"LOADVSLICE\n");
        fprintf(scriptoutstream," %s\n",slicei->label.longlabel);
        fprintf(scriptoutstream," %i %f\n",slicei->idir,slicei->position_orig);
        script_multivslice=1;
      }
    }
    if(scriptoutstream==NULL||defer_file_loading==0){
      for(i=0;i<mvslicei->nvslices;i++){
        vslicedata *vslicei;

        vslicei = vsliceinfo + mvslicei->ivslices[i];
        if(vslicei->skip==0&&vslicei->loaded==0)LoadVSliceMenu(mvslicei->ivslices[i]);
        if(vslicei->skip==1&&vslicei->loaded==1)UnloadVSliceMenu(mvslicei->ivslices[i]);
      }
    }
    script_multivslice=0;
  }
  else{
    switch(value){
      case -20:
        showallslicevectors=1-showallslicevectors;
        updatemenu=1;
        glutPostRedisplay();
        return;
      case UNLOAD_ALL:
        LoadVSliceMenu(UNLOAD_ALL);
        break;

#ifdef pp_SLICEDUP
      case MENU_KEEP_ALL:
      if(vectorslicedup_option!=SLICEDUP_KEEPALL){
        vectorslicedup_option = SLICEDUP_KEEPALL;
        updatemenu = 1;
        glutPostRedisplay();
        UpdateVSliceDups();
        update_slicedup_dialog();
      }
      break;

      case  MENU_KEEP_COARSE:
      if(vectorslicedup_option!=SLICEDUP_KEEPCOARSE){
        vectorslicedup_option = SLICEDUP_KEEPCOARSE;
        updatemenu = 1;
        glutPostRedisplay();
        UpdateVSliceDups();
        update_slicedup_dialog();
      }
      break;

      case MENU_KEEP_FINE:
      if(vectorslicedup_option!=SLICEDUP_KEEPFINE){
        vectorslicedup_option = SLICEDUP_KEEPFINE;
        updatemenu = 1;
        glutPostRedisplay();
        UpdateVSliceDups();
        update_slicedup_dialog();
      }
      break;
#endif
      default:
        ASSERT(FFALSE);
        break;
    }
  }
}

/* ------------------ LoadAllMSlices ------------------------ */

void LoadAllMSlices(int last_slice, multislicedata *mslicei){
  int i;

  for(i = 0; i < mslicei->nslices; i++){
    slicedata *slicei;
    int set_slicecolor;

    slicei = sliceinfo + mslicei->islices[i];
    set_slicecolor = DEFER_SLICECOLOR;
    if(last_slice == i)set_slicecolor = SET_SLICECOLOR;
    if(slicei->skip == 0 && slicei->loaded == 0){
      LoadSlicei(set_slicecolor,mslicei->islices[i]);
    }
  }
}

/* ------------------ LoadMultiSliceMenu ------------------------ */

void LoadMultiSliceMenu(int value){
  int i;
  multislicedata *mslicei;

  if(value==MENU_DUMMY)return;
  if(value>=0){
    mslicei = multisliceinfo + value;
    if(scriptoutstream!=NULL){
      if(mslicei->nslices>0){
        slicedata *slicei;

        slicei = sliceinfo + mslicei->islices[0];
        fprintf(scriptoutstream,"LOADSLICE\n");
        fprintf(scriptoutstream," %s\n",slicei->label.longlabel);
        fprintf(scriptoutstream," %i %f\n",slicei->idir,slicei->position_orig);
        script_multislice=1;
      }
    }
    if(scriptoutstream==NULL||defer_file_loading==0){
      int last_slice;

      last_slice = mslicei->nslices - 1;
      for(i = mslicei->nslices-1; i >=0; i--){
        slicedata *slicei;

        slicei = sliceinfo + mslicei->islices[i];
        if(slicei->skip == 0 && slicei->loaded == 0){
          last_slice = i;
          break;
        }
      }
      for(i = 0; i < mslicei->nslices; i++){
        slicedata *slicei;

        slicei = sliceinfo + mslicei->islices[i];
        if(slicei->skip == 1 && slicei->loaded == 1)UnloadSliceMenu(mslicei->islices[i]);
      }
      LoadAllMSlices(last_slice,mslicei);
      if(mslicei->nslices>0&&mslicei->islices[0]>=nsliceinfo-nfedinfo){
        output_mfed_csv(mslicei);
      }
    }
    script_multislice=0;
  }
  else{
    switch(value){
      case UNLOAD_ALL:
      LoadSliceMenu(UNLOAD_ALL);
      break;
#ifdef pp_SLICEDUP
      case MENU_KEEP_ALL:
      if(slicedup_option!=SLICEDUP_KEEPALL){
        slicedup_option = SLICEDUP_KEEPALL;
        updatemenu = 1;
        glutPostRedisplay();
        UpdateSliceDups();
        update_slicedup_dialog();
      }
      break;

      case  MENU_KEEP_COARSE:
      if(slicedup_option!=SLICEDUP_KEEPCOARSE){
        slicedup_option = SLICEDUP_KEEPCOARSE;
        updatemenu = 1;
        glutPostRedisplay();
        UpdateSliceDups();
        update_slicedup_dialog();
      }
      break;

      case MENU_KEEP_FINE:
      if(slicedup_option!=SLICEDUP_KEEPFINE){
        slicedup_option = SLICEDUP_KEEPFINE;
        updatemenu = 1;
        glutPostRedisplay();
        UpdateSliceDups();
        update_slicedup_dialog();
      }
      break;
#endif
      case MENU_SLICECOLORDEFER:
        use_set_slicecolor = 1 - use_set_slicecolor;
        updatemenu = 1;
        break;
      default:
      ASSERT(FFALSE);
      break;
    }
  }
}

/* ------------------ Plot3DListMenu ------------------------ */

void Plot3DListMenu(int value){
  int i;
  plot3ddata *plot3di;

  if(value<0||value>=nplot3dtimelist)return;
  LoadPlot3dMenu(UNLOAD_ALL);
  if(scriptoutstream!=NULL){
    fprintf(scriptoutstream,"LOADPLOT3D\n");
    fprintf(scriptoutstream," %f\n",plot3dtimelist[value]);
  }
  for(i=0;i<nplot3dinfo;i++){
    plot3di = plot3dinfo + i;
    if(ABS(plot3di->time-plot3dtimelist[value])<0.5){
      LoadPlot3dMenu(i);
    }
  }
}

/* ------------------ UpdateMenu ------------------------ */

void UpdateMenu(void){
  updatemenu=1;
  glutPostRedisplay();
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ LoadPlot3DMenu ------------------------ */

void LoadPlot3dMenu(int value){
  int errorcode;
  int i;

  if(value==MENU_PLOT3D_DUMMY)return;
  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value>=0){
    char *plot3dfile;

    ReadPlot3dFile=1;
    plot3dfile = plot3dinfo[value].file;
    if(scriptoutstream!=NULL&&loadplot3dall==0){
      fprintf(scriptoutstream,"LOADPLOT3D\n");
      fprintf(scriptoutstream," %i %f\n",
        plot3dinfo[value].blocknumber+1,plot3dinfo[value].time);
    }
    if(scriptoutstream==NULL||defer_file_loading==0){
      readplot3d(plot3dfile,value,LOAD,&errorcode);
    }
  }
  else if(value==UNLOAD_ALL){
    for(i=0;i<nplot3dinfo;i++){
      readplot3d("",i,UNLOAD,&errorcode);
    }
  }
  else{
    value+=100000;
    loadplot3dall=1;
    Plot3DListMenu(value);
    loadplot3dall=0;
  }
  updatemenu=1;
  glutPostRedisplay();
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ LoadIsoi ------------------------ */

void LoadIsoi(int value){
  char *file;
  isodata *isoi;
  int errorcode;

  ReadIsoFile=1;
  isoi = isoinfo + value;
  file=isoi->file;
  isoi->loading=1;
  if(script_iso==0&&scriptoutstream!=NULL){
    fprintf(scriptoutstream,"LOADISOM\n");
    fprintf(scriptoutstream, " %s\n", isoi->surface_label.longlabel);
    fprintf(scriptoutstream, " %i\n", isoi->blocknumber+1);
  }
  if(scriptoutstream==NULL||defer_file_loading==0){
    readiso(file,value,LOAD,NULL,&errorcode);
    if(update_readiso_geom_wrapup == UPDATE_ISO_ONE_NOW)readiso_geom_wrapup();
  }
  isoi->loading=0;
}

  /* ------------------ LoadAllIsos ------------------------ */

void LoadAllIsos(int iso_type){
  int i;

  for(i = 0; i < nisoinfo; i++){
    isodata *isoi;

    isoi = isoinfo + i;
    if(iso_type == isoi->type)LoadIsoi(i);
  }
}

/* ------------------ LoadIsoMenu ------------------------ */

void LoadIsoMenu(int value){
  int errorcode;
  int i;
  int ii;

  if(value==MENU_DUMMY3)return;
  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value>=0){
    LoadIsoi(value);
  }
  if(value==-1){
    for(i=0;i<nisoinfo;i++){
      isodata *isoi;

      isoi = isoinfo + i;
      if(isoi->loaded==1)readiso("",i,UNLOAD,NULL,&errorcode);
    }
  }
  if(value<=-10){
    isodata *isoi;

    ii = -(value + 10);
    isoi = isoinfo + ii;
    if(scriptoutstream!=NULL){
      script_iso=1;
      fprintf(scriptoutstream,"LOADISO\n");
      fprintf(scriptoutstream," %s\n",isoi->surface_label.longlabel);
    }
    if(scriptoutstream==NULL||defer_file_loading==0){
      update_readiso_geom_wrapup = UPDATE_ISO_START_ALL;
      LoadAllIsos(isoi->type);
      if(update_readiso_geom_wrapup == UPDATE_ISO_ALL_NOW)readiso_geom_wrapup();
      update_readiso_geom_wrapup = UPDATE_ISO_OFF;
    }
    script_iso=0;
  }
  updatemenu=1;
  glutPostRedisplay();
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ LoadPatchMenu ------------------------ */

void LoadPatchMenu(int value){
  int errorcode;
  int i,ii;
  int patchtypenew;

  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value>=0){
    patchtypenew=GetPatchType(patchinfo+value);
    if(patchtypenew!=-1){
      for(ii=0;ii<npatch_loaded;ii++){
        patchdata *patchi;

        i = patch_loaded_list[ii];
        patchi = patchinfo + i;
        if(patchi->type!=patchtypenew)readpatch(i,UNLOAD,&errorcode);
      }
    }
    if(scriptoutstream!=NULL){
      patchdata *patchi;

      patchi = patchinfo + value;
      fprintf(scriptoutstream,"// LOADFILE\n");
      fprintf(scriptoutstream,"//  %s\n",patchi->file);
      fprintf(scriptoutstream, "LOADBOUNDARYM\n");
      fprintf(scriptoutstream, " %s\n", patchi->label.longlabel);
      fprintf(scriptoutstream, " %i\n", patchi->blocknumber+1);
    }
    if(scriptoutstream==NULL||defer_file_loading==0){
      LOCK_COMPRESS
      readpatch(value,LOAD,&errorcode);
      UNLOCK_COMPRESS
    }
  }
  else if(value<=-10){
    patchdata *patchj;

    value = -(value + 10);
    patchj = patchinfo + value;
    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"LOADBOUNDARY\n");
      fprintf(scriptoutstream," %s\n",patchj->label.longlabel);
    }
    if(scriptoutstream==NULL||defer_file_loading==0){
      for(i=0;i<npatchinfo;i++){
        patchdata *patchi;

        patchi = patchinfo + i;
        if(strcmp(patchi->label.longlabel,patchj->label.longlabel)==0&&patchi->filetype==patchj->filetype){
          LOCK_COMPRESS
          readpatch(i,LOAD,&errorcode);
          UNLOCK_COMPRESS
        }
      }
    }
    force_redisplay=1;
    UpdateFrameNumber(0);
  }
  else if(value==MENU_UPDATEBOUNDS){
    Update_All_Patch_Bounds();
  }
  else{
    for(i=0;i<npatchinfo;i++){
      readpatch(i,UNLOAD,&errorcode);
    }
  }
  updatemenu=1;
  glutPostRedisplay();
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ ShowPatchMenu ------------------------ */

void ShowPatchMenu(int value){
  updatemenu=1;
  updatefacelists=1;
  glutPostRedisplay();
  if(value>=1000){
    patchdata *patchi;

    patchi=patchinfo+value-1000;
    if(patchi->type==ipatchtype){
      patchi->display=1-patchi->display;
    }
    else{
      patchi->display=1;
      ipatchtype=GetPatchType(patchi);
    }
    UpdatePatchType();
  }
  if(value==SHOW_CHAR){
    vis_threshold = 1 - vis_threshold;
    updatechar();
  }
  if(value==SHOWALL_BOUNDARY){
    int ii;

    for(ii=0;ii<npatch_loaded;ii++){
      patchdata *patchi;
      int i;

      i = patch_loaded_list[ii];
      patchi = patchinfo + i;
      patchi->display=1;
    }
  }
  if(value==HIDEALL_BOUNDARY){
    int ii;

    for(ii=0;ii<npatch_loaded;ii++){
      patchdata *patchi;
      int i;

      i = patch_loaded_list[ii];
      patchi = patchinfo + i;
      patchi->display=0;
    }
  }
  if(value<0){
    if(value==EXTERIORwallmenu){
      int i,n,val;

      allexterior = 1-allexterior;
  	  showexterior=1-showexterior;
      val = allexterior;
      for(n=0;n<current_mesh->npatches;n++){
        if(current_mesh->patchtype[n]!=INTERIORwall){
          current_mesh->visPatches[n]=val;
        }
      }
      for(i=1;i<7;i++){
        visPatchType[i]=val;
      }
    }
    else if(value==INTERIORwallmenu){
      int n,val;

      allinterior = 1 - allinterior;
      val = allinterior;
      visPatchType[INTERIORwall]=val;
      for(n=0;n<current_mesh->npatches;n++){
        if(current_mesh->patchtype[n]==INTERIORwall){
          current_mesh->visPatches[n]=val;
        }
      }
    }
    else if(value == SOLIDpatchmenu){
      show_patch_solid = 1 - show_patch_solid;
    }
    else if(value == OUTLINEpatchmenu){
      show_patch_outline = 1 - show_patch_outline;
    }
    else if(value == POINTSpatchmenu){
      show_patch_verts = 1 - show_patch_verts;
    }
    else if(value==INSOLIDpatchmenu){
      show_patch_insolid = 1-show_patch_insolid;
    }
    else if(value==INGASpatchmenu){
      show_patch_ingas = 1-show_patch_ingas;
    }
    else if(value == INCUTCELLpatchmenu){
      show_patch_incutcell = 1 - show_patch_incutcell;
    }
    else if(value != DUMMYwallmenu){
      int n;

      value = -(value+2); /* map xxxwallmenu to xxxwall */
      for(n=0;n<current_mesh->npatches;n++){
        if(current_mesh->patchtype[n]==value){
          current_mesh->visPatches[n] = 1 - current_mesh->visPatches[n];
          visPatchType[value]=current_mesh->visPatches[n];
        }
      }
    }
  }
  plotstate=GetPlotState(DYNAMIC_PLOTS);
}

/* ------------------ VentMenu ------------------------ */

void VentMenu(int value){
  if(value==-1)return;
  switch(value){
  case SHOW_ALL_VENTS: // show all vents
    visVents=1;
    visOpenVents=1;
    visDummyVents=1;
    visOtherVents=1;
    visCircularVents=VENT_CIRCLE;
    break;
  case MENU_VENT_OPEN:
    visOpenVents=1-visOpenVents;
    break;
  case MENU_VENT_OUTLINE:
    visOpenVentsAsOutline = 1 - visOpenVentsAsOutline;
    break;
  case MENU_VENT_EXTERIOR:
    visDummyVents = 1 - visDummyVents;
    break;
  case MENU_VENT_TWOINTERIOR:
     show_bothsides_int=1-show_bothsides_int;
     updatefaces=1;
     break;
  case MENU_VENT_TWOEXTERIOR:
     show_bothsides_ext = 1 - show_bothsides_ext;
     updatefaces=1;
     break;
  case MENU_VENT_TRANSPARENT:
     show_transparent_vents=1-show_transparent_vents;
     updatefaces=1;
     break;
  case MENU_VENT_OTHER:
     visOtherVents=1-visOtherVents;
     break;
   case HIDE_ALL_VENTS: // Hide all vents
     visVents=0;
     visOpenVents=0;
     visDummyVents=0;
     visOtherVents=0;
     visCircularVents=VENT_HIDE;
     break;
   case MENU_VENT_CIRCLE:
     visCircularVents=VENT_CIRCLE;
     break;
   case MENU_VENT_RECTANGLE:
     visCircularVents=VENT_RECTANGLE;
     break;
   case MENU_VENT_CIRCLEHIDE:
     visCircularVents=VENT_HIDE;
     break;
   case MENU_VENT_CIRCLEOUTLINE:
     circle_outline=1-circle_outline;
     break;
  default:
    ASSERT(FFALSE);
    break;
  }
  updatefacelists=1;
  updatemenu=1;
  glutPostRedisplay();
}
#define GEOMETRY_SOLID 0
#define GEOMETRY_OUTLINE 1
#define GEOMETRY_SOLIDOUTLINE 2
#define GEOMETRY_INTERIOR_SOLID 9
#define GEOMETRY_INTERIOR_OUTLINE 12
#define GEOMETRY_HIDE 7
#define GEOMETRY_TETRA_HIDE 11
#define GEOMETRY_SHOWNORMAL 3
#define GEOMETRY_SORTFACES 6
#define GEOMETRY_SMOOTHNORMAL 4
#define GEOMETRY_SHOWDIAGNOSTICS 13
#define GEOMETRY_HILIGHTSKINNY 5
#define GEOMETRY_HIDEALL 8

/* ------------------ ImmersedMenu ------------------------ */

void ImmersedMenu(int value){
  updatemenu=1;
  switch(value){

    case GEOMETRY_INTERIOR_SOLID:
      show_volumes_solid=1-show_volumes_solid;
      break;
    case GEOMETRY_INTERIOR_OUTLINE:
      show_volumes_outline=1-show_volumes_outline;
      break;
    case GEOMETRY_TETRA_HIDE:
      show_volumes_solid=0;
      show_volumes_outline=0;
      break;
    case GEOMETRY_SOLIDOUTLINE:
      if(show_faces_solid==1&&show_faces_outline==1){
        show_faces_solid=1;
        show_faces_outline=0;
      }
      else{
        show_faces_solid=1;
        show_faces_outline=1;
      }
      break;
    case GEOMETRY_SOLID:
      if(show_faces_solid==1&&show_faces_outline==1){
        show_faces_solid=1;
        show_faces_outline=0;
      }
      else if(show_faces_solid==1&&show_faces_outline==0){
        show_faces_solid=0;
        show_faces_outline=1;
      }
      else if(show_faces_solid==0&&show_faces_outline==1){
        show_faces_solid=1;
        show_faces_outline=0;
      }
      else{
        show_faces_solid=1;
        show_faces_outline=0;
      }
      break;
    case GEOMETRY_OUTLINE:
      if(show_faces_solid==1&&show_faces_outline==1){
        show_faces_solid=0;
        show_faces_outline=1;
      }
      else if(show_faces_solid==1&&show_faces_outline==0){
        show_faces_solid=0;
        show_faces_outline=1;
      }
      else if(show_faces_solid==0&&show_faces_outline==1){
        show_faces_solid=1;
        show_faces_outline=0;
      }
      else{
        show_faces_solid=0;
        show_faces_outline=1;
      }
      break;
    case GEOMETRY_SHOWNORMAL:
      show_geom_normal=1-show_geom_normal;
      break;
    case GEOMETRY_SMOOTHNORMAL:
      smooth_geom_normal=1-smooth_geom_normal;
      break;
    case GEOMETRY_HILIGHTSKINNY:
      hilight_skinny = 1 - hilight_skinny;
      break;
    case GEOMETRY_SORTFACES:
      sort_geometry=1-sort_geometry;
      break;
    case GEOMETRY_SHOWDIAGNOSTICS:
      show_geometry_diagnostics = 1 - show_geometry_diagnostics;
      break;
    case GEOMETRY_HIDE:
      show_faces_solid=0;
      show_faces_outline=0;
      break;
    case GEOMETRY_HIDEALL:
      ImmersedMenu(GEOMETRY_HIDE);
      ImmersedMenu(GEOMETRY_TETRA_HIDE);
      show_geom_normal = 0;
      break;
    case MENU_DUMMY:
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
  update_geometry_controls();

  glutPostRedisplay();
}

/* ------------------ BlockageMenu ------------------------ */

void BlockageMenu(int value){
  int change_state=0;

  if(solid_state<0)solid_state=visBlocks;
  if(outline_state<0)outline_state=OUTLINE_NONE;
  switch(value){
    case visBLOCKOutlineColor:
      outline_color_flag = 1 - outline_color_flag;
      updatefaces=1;
      break;
    case visBLOCKOnlyOutline:
      if(outline_state!=OUTLINE_ONLY){
        outline_state=OUTLINE_ONLY;
      }
      else{
        outline_state=OUTLINE_NONE;
      }
      change_state=1;
      break;
    case visCADOpaque:
      viscadopaque = 1 - viscadopaque;
      break;
    case visBLOCKAddOutline:
      if(outline_state!=OUTLINE_ADDED){
        outline_state=OUTLINE_ADDED;
        if(solid_state==visBLOCKHide)solid_state=visBLOCKAsInput;
      }
      else{
        outline_state=OUTLINE_NONE;
      }
      change_state=1;
      break;
    case visBLOCKAsInput:
      solid_state=visBLOCKAsInput;
      if(outline_state==OUTLINE_ONLY)outline_state=OUTLINE_ADDED;
      change_state=1;
      break;
    case visBLOCKNormal:
      solid_state=visBLOCKNormal;
      if(outline_state==OUTLINE_ONLY)outline_state=OUTLINE_ADDED;
      change_state=1;
      break;
    case visBLOCKHide:
      outline_state=OUTLINE_NONE;
      solid_state=visBLOCKHide;
      change_state=1;
      break;
    case BLOCKlocation_grid:
    case BLOCKlocation_exact:
    case BLOCKlocation_cad:
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
  if(change_state==1){
    switch(outline_state){
      case OUTLINE_NONE:
        value=solid_state;
        break;
      case OUTLINE_ONLY:
        value=visBLOCKOutline;
        break;
      case OUTLINE_ADDED:
        switch(solid_state){
          case visBLOCKAsInput:
            value=visBLOCKAsInputOutline;
          break;
          case visBLOCKNormal:
            value=visBLOCKSolidOutline;
            break;
          case BLOCKAGE_HIDDEN:
            value=visBLOCKOutline;
            break;
          default:
            ASSERT(FFALSE);
            break;
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
  }

  switch(value){
   case visCADOpaque:
   case visBLOCKOutlineColor:
     break;
   case visBLOCKAsInputOutline:
   case visBLOCKAsInput:
     visBlocks=value;
     update_trainer_outline();
     break;
   case visBLOCKNormal:
   case visBLOCKOutline:
   case visBLOCKHide:
   case visBLOCKSolidOutline:
     visBlocks=value;
     if(value==visBLOCKSolidOutline||visBLOCKold==visBLOCKSolidOutline)updatefaces=1;
     update_trainer_outline();
     break;
   case BLOCKlocation_grid:
   case BLOCKlocation_exact:
   case BLOCKlocation_cad:
     blocklocation=value;
     break;
   case BLOCKtexture_cad:
     visCadTextures=1-visCadTextures;
     break;
   case visBLOCKTransparent:
     visTransparentBlockage=1-visTransparentBlockage;
     break;
   default:
     if(value<0){
       value=-value-1;
       if(value>=0&&value<=npropinfo-1){
         propdata *propi;

         propi = propinfo + value;
         propi->blockvis=1-propi->blockvis;
       }
     }
     else{
       ASSERT(FFALSE);
     }
     break;
  }
  visBLOCKold=value;
  updatemenu=1;
 // updatefaces=1;
  updatefacelists=1;
  updatehiddenfaces=1;
  glutPostRedisplay();
}

/* ------------------ RotateTypeMenu ------------------------ */

void RotateTypeMenu(int value){
  if(value==MENU_DUMMY)return;
  rotation_type = value;
  update_rotation_type(rotation_type);
  rotation_type_CB(rotation_type);
  updatemenu=1;
  glutPostRedisplay();
}

/* ------------------ TitleMenu ------------------------ */

void TitleMenu(int value){
  updatemenu=1;
  glutPostRedisplay();
  switch(value){
  case 0:
    visTitle = 1 - visTitle;
    break;
  case 3:
    visTitle=1;
    break;
  case 4:
    visTitle=0;
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
}

/* ------------------ PropMenu ------------------------ */

void PropMenu(int value){
  int iprop, iobject;

  // value = iobject*npropinfo + iprop

  iprop = value%npropinfo;
  iobject = value/npropinfo;
  if(iprop>=0&&iprop<npropinfo){
    propdata *propi;

    propi = propinfo + iprop;
    if(iobject>=0&&iobject<propi->nsmokeview_ids){
      int i;

      propi->smokeview_id=propi->smokeview_ids[iobject];
      propi->smv_object=propi->smv_objects[iobject];
      updatemenu=1;
      get_indep_var_indices(propi->smv_object,
        propi->vars_indep,propi->nvars_indep,
        propi->vars_indep_index);

      for(i=0;i<npartclassinfo;i++){
        partclassdata *partclassi;

        partclassi = partclassinfo + i;
        update_partclass_depend(partclassi);

      }

      glutPostRedisplay();
    }
  }
}

/* ------------------ ShowObjectsMenu ------------------------ */

void ShowObjectsMenu(int value){
  sv_object *objecti;
  int i;

  if(value>=0&&value<nobject_defs){
    objecti = object_defs[value];
    objecti->visible = 1 - objecti->visible;
  }
  else if(value == OBJECT_MISSING){
    updatemenu = 1;
    show_missing_objects = 1 - show_missing_objects;
  }
  else if(value==OBJECT_SHOWALL){
    for(i=0;i<nobject_defs;i++){
      objecti = object_defs[i];
      objecti->visible=1;
    }
  }
  else if(value==OBJECT_HIDEALL){
    for(i=0;i<nobject_defs;i++){
      objecti = object_defs[i];
      objecti->visible=0;
    }
  }
  else if(value==OBJECT_SELECT){
    select_device=1-select_device;
  }
  else if(value==OBJECT_OUTLINE){
    object_outlines=1-object_outlines;
  }
  else if(value==OBJECT_ORIENTATION){
    show_device_orientation=1-show_device_orientation;
    update_device_orientation();
  }
  else if(value==MENU_DUMMY){
  }
  else{
    device_sphere_segments=ABS(value);
    Init_Sphere(device_sphere_segments,2*device_sphere_segments);
    Init_Circle(2*device_sphere_segments,&object_circ);
  }
  updatemenu=1;
  glutPostRedisplay();
}

/* ------------------ ZoneShowMenu ------------------------ */

void ZoneShowMenu(int value){
  switch(value){
  case MENU_DUMMY:
    return;
  case MENU_ZONE_HORIZONTAL:
    if(zonecolortype==ZONESMOKE_COLOR)zonecolortype = ZONETEMP_COLOR;
    visVZone=0;
    visHZone=1;
    visZone=1;
    break;
  case MENU_ZONE_VERTICAL:
    if(zonecolortype==ZONESMOKE_COLOR)zonecolortype = ZONETEMP_COLOR;
    visVZone = 1;
    visHZone=0;
    visZone=1;
    break;
  case MENU_ZONE_LAYERHIDE:
    visVZone=0;
    visHZone=0;
    visSZone=0;
    break;
  case MENU_ZONE_2DHAZARD:
    zonecolortype=ZONEHAZARD_COLOR;
    visSZone=0;
    if(visVZone==0&&visHZone==0)visVZone=1;
    visZone=1;
    break;
  case MENU_ZONE_2DTEMP:
    zonecolortype=ZONETEMP_COLOR;
    visSZone=0;
    if(visVZone==0&&visHZone==0)visVZone=1;
    visZone=1;
    break;
  case MENU_ZONE_3DSMOKE:
    zonecolortype=ZONESMOKE_COLOR;
    visSZone=1;
    visZone=1;
    break;
  case MENU_ZONE_VENTS:
    visVentFlow=1-visVentFlow;
    if(visVentFlow==1){
      visVentHFlow=1;
      visVentVFlow=1;
      visVentMFlow=1;
    }
    else{
      visVentHFlow=0;
      visVentVFlow=0;
      visVentMFlow=0;
    }
    break;
#define VISVENTFLOW     if((nzhvents>0&&visVentHFlow==1)||(nzvvents>0&&visVentVFlow==1)||(nzmvents>0&&visVentMFlow==1)){\
      visVentFlow = 1;\
    }\
    else{\
      visVentFlow = 0;\
    }
  case MENU_ZONE_HVENTS:
    visVentHFlow = 1-visVentHFlow;
    VISVENTFLOW;
    break;
  case MENU_ZONE_VVENTS:
    visVentVFlow = 1-visVentVFlow;
    VISVENTFLOW;
    break;
  case MENU_ZONE_MVENTS:
    visVentMFlow = 1-visVentMFlow;
    VISVENTFLOW;
    break;
  case MENU_ZONE_VENT_SLAB:
    visventslab = 1 - visventslab;
    if(visventslab==1)visventprofile=0;
    if(visventprofile==1||visventslab==1){
      visVentHFlow=1;
    }
    else{
      visVentHFlow=0;
    }
    VISVENTFLOW;
    break;
  case MENU_ZONE_VENT_PROFILE:
    visventprofile = 1 - visventprofile;
    if(visventprofile==1)visventslab=0;
    if(visventprofile==1||visventslab==1){
      visVentHFlow=1;
    }
    else{
      visVentHFlow=0;
    }
    VISVENTFLOW;
    break;
  case MENU_ZONE_FIRES:
    viszonefire=1-viszonefire;
    break;
  default:
    ASSERT(FFALSE);
  }
  updatemenu=1;
  glutPostRedisplay();
}

#define GEOM_Vents 15
#define GEOM_Outline 3
#define GEOM_TriangleCount 14
#define GEOM_ShowAll 11
#define GEOM_HideAll 13

/* ------------------ GeometryMenu ------------------------ */

void GeometryMenu(int value){

  switch(value){
  case GEOM_TriangleCount:
    show_triangle_count=1-show_triangle_count;
    break;
  case GEOM_Outline:
    if(isZoneFireModel==0)visFrame=1-visFrame;
    break;
  case 5:
    visFloor=1-visFloor;
    break;
  case 6:
    visWalls=1-visWalls;
    break;
  case 7:
    visCeiling=1-visCeiling;
    break;

  case 17+TERRAIN_3D:
  case 17+TERRAIN_2D_STEPPED:
  case 17+TERRAIN_2D_LINE:
  case 17+TERRAIN_3D_MAP:
  case 17+TERRAIN_HIDDEN:
    if(value==17+TERRAIN_HIDDEN){
      BlockageMenu(visBlocksSave);
      if(visOtherVents!=visOtherVentsSAVE){
        visOtherVents=visOtherVentsSAVE;
      }
    }
    else{
      if(visOtherVents!=0){
        visOtherVentsSAVE=visOtherVents;
        visOtherVents=0;
      }
    }
    visTerrainType=value-17;
    if(visTerrainType==TERRAIN_3D){
      planar_terrain_slice=0;
    }
    else{
      planar_terrain_slice=1;
    }
    Update_Glui_Wui();
    break;
  case GEOM_ShowAll:
    if(isZoneFireModel)visFrame=1;
    /*
    visFloor=1;
    visWalls=1;
    visCeiling=1;
    */
    visVents=1;
    BlockageMenu(visBLOCKAsInput);
    break;
  case GEOM_HideAll:
    visFrame=0;
    visFloor=0;
    visWalls=0;
    visCeiling=0;
    visVents=0;
    visGrid=0;
    BlockageMenu(visBLOCKHide);
    break;
  case GEOM_Vents:
    visVents=1-visVents;
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  updatefacelists=1;
  updatemenu=1;
  glutPostRedisplay();
}

/* ------------------ GetNumActiveDevices ------------------------ */

int GetNumActiveDevices(void){
  int num_activedevices = 0;

  if(nobject_defs > 0){
    int i;

    for(i = 0; i < nobject_defs; i++){
      sv_object *obj_typei;

      obj_typei = object_defs[i];
      if(obj_typei->used_by_device == 1)num_activedevices++;
    }
  }
  return num_activedevices;
}

/* ------------------ GetNTotalVents ------------------------ */

int GetNTotalVents(void){
  int ntotal_vents = 0;
  int i;

  for(i = 0; i < nmeshes; i++){
    meshdata *meshi;

    meshi = meshinfo + i;
    ntotal_vents += meshi->nvents;
  }
  return ntotal_vents;
}

/* ------------------ IsPatchType ------------------------ */

int IsPatchType(int type){
  int i;

  for(i = 0; i < nmeshes; i++){
    meshdata *meshi;
    int n;

    meshi = meshinfo + i;
    for(n = 0; n < meshi->npatches; n++){
      if(meshi->patchtype[n] == type)return 1;
    }
  }
  return 0;
}

/* ------------------ InitMenus ------------------------ */

void InitMenus(int unload){
  int i;
  int nsmoke3dloaded,nvolsmoke3dloaded;
  int nsliceloaded,nvsliceloaded2,nvsliceloaded,nmultisliceloaded;
  int nvslice0, nvslice1, nvslice2,nvsliceloaded0,nvsliceloaded1;
  int npartloaded,npart5loaded,nevacloaded;
  int npatchloaded;
  int nplot3dloaded;
  int nisoloaded;
  int showhide_data = 0;

static int filesdialogmenu = 0, viewdialogmenu = 0, datadialogmenu = 0, windowdialogmenu=0;
static int labelmenu=0, colorbarmenu=0, colorbarsmenu=0, colorbarshademenu, smokecolorbarmenu=0, showhidemenu=0;
static int optionmenu=0, rotatetypemenu=0;
static int resetmenu=0, frameratemenu=0, rendermenu=0, smokeviewinimenu=0, inisubmenu=0, resolutionmultipliermenu=0;
static int startrenderingmenu=0;
#ifdef pp_COMPRESS
static int compressmenu=0;
#endif
static int showhideslicemenu=0,showvslicemenu=0;
static int plot3dshowmenu=0, staticvariablemenu=0, helpmenu=0, webhelpmenu=0, keyboardhelpmenu=0, mousehelpmenu=0;
static int vectorskipmenu=0,unitsmenu=0;
static int isosurfacemenu=0, isovariablemenu=0, levelmenu=0;
static int fontmenu=0, aperturemenu=0,dialogmenu=0,zoommenu=0;
static int gridslicemenu=0, blockagemenu=0, immersedmenu=0, immersedinteriormenu=0, immersedsurfacemenu=0, loadpatchmenu=0, ventmenu=0, circularventmenu=0;
static int loadisomenu=0, isosurfacetypemenu=0;
static int geometrymenu=0, loadunloadmenu=0, reloadmenu=0, aboutmenu=0, disclaimermenu=0, terrain_showmenu=0;
static int scriptmenu=0;
static int scriptlistmenu=0,scriptsteplistmenu=0,scriptrecordmenu=0;
#ifdef pp_LUA
static int luascriptmenu=0;
static int luascriptlistmenu=0;
#endif
static int loadplot3dmenu=0, unloadvslicemenu=0, unloadslicemenu=0;
static int loadterrainmenu=0, unloadterrainmenu=0;
static int loadsmoke3dmenu=0,loadsmoke3dsootmenu=0,loadsmoke3dhrrmenu=0,loadsmoke3dwatermenu=0;
static int loadvolsmoke3dmenu=0,showvolsmoke3dmenu=0;
static int unloadsmoke3dmenu=0,unloadvolsmoke3dmenu=0;
static int unloadevacmenu=0, unloadpartmenu=0, loadslicemenu=0, loadmultislicemenu=0;
static int *loadsubvslicemenu=NULL, nloadsubvslicemenu=0;
static int *loadsubslicemenu=NULL, nloadsubslicemenu=0, iloadsubslicemenu=0;
static int *loadsubmslicemenu=NULL, nloadsubmslicemenu=0;
static int *loadsubmvslicemenu=NULL, nloadsubmvslicemenu=0;
static int *loadsubplot3dmenu=NULL, nloadsubplot3dmenu=0;
static int loadmultivslicemenu=0, unloadmultivslicemenu=0;
static int unloadmultislicemenu=0, vslicemenu=0, staticslicemenu=0;
static int evacmenu=0, particlemenu=0, particlesubmenu=0, showpatchmenu=0, zonemenu=0, isoshowmenu=0, isoshowsubmenu=0, isolevelmenu=0, smoke3dshowmenu=0;
static int particlepropshowmenu=0,humanpropshowmenu=0;
static int *particlepropshowsubmenu=NULL;
static int particlestreakshowmenu=0;
static int tourmenu=0;
static int avatartourmenu=0,avatarevacmenu=0;
static int trainerviewmenu=0,mainmenu=0,zoneshowmenu=0,particleshowmenu=0,evacshowmenu=0;
static int showobjectsmenu=0,spheresegmentmenu=0,propmenu=0;
static int unloadplot3dmenu=0, unloadpatchmenu=0, unloadisomenu=0;
static int showmultislicemenu=0;
static int textureshowmenu=0;
#ifdef _DEBUG
static int menu_count=0;
static int in_menu=0;
#endif

updatemenu=0;
#ifdef _DEBUG
  PRINTF("Updating Menus %i In menu %i\n",menu_count++,in_menu);
  in_menu=1;
#endif
  update_showhidebuttons();
  glutPostRedisplay();

  nsliceloaded=0;
  for(i=0;i<nsliceinfo;i++){
    slicedata *sd;

    sd = sliceinfo + i;
    if(sd->loaded==1)nsliceloaded++;
  }

  nsmoke3dloaded=0;
  for(i=0;i<nsmoke3dinfo;i++){
    smoke3ddata *smoke3di;

    smoke3di=smoke3dinfo+i;
    if(smoke3di->loaded==1)nsmoke3dloaded++;
  }

  nvolsmoke3dloaded=0;
  for(i=0;i<nmeshes;i++){
    meshdata *meshi;
    volrenderdata *vr;

    meshi = meshinfo + i;
    vr = &(meshi->volrenderinfo);
    if(vr->fireslice==NULL||vr->smokeslice==NULL)continue;
    if(vr->loaded==1){
      nvolsmoke3dloaded++;
    }
  }

  nvsliceloaded=0;
  for(i=0;i<nvsliceinfo;i++){
    vslicedata *vd;

    vd = vsliceinfo + i;
    if(vd->loaded==1)nvsliceloaded++;
  }

  for(i=0;i<nmultisliceinfo;i++){
    multislicedata *mslicei;
    int j;

    mslicei = multisliceinfo + i;
    mslicei->loaded=0;
    mslicei->display=0;
    for(j=0;j<mslicei->nslices;j++){
      slicedata *sd;

      sd = sliceinfo + mslicei->islices[j];
      if(sd->loaded==1)mslicei->loaded++;
      if(sd->display==1)mslicei->display++;
    }
    if(mslicei->loaded>0&&mslicei->loaded<mslicei->nslices){
      mslicei->loaded=-1;
    }
    else if(mslicei->loaded==mslicei->nslices){
      mslicei->loaded=1;
    }
    if(mslicei->display>0&&mslicei->display<mslicei->nslices){
      mslicei->display=-1;
    }
    else if(mslicei->display==mslicei->nslices){
      mslicei->display=1;
    }
  }
  for(i=0;i<nmultivsliceinfo;i++){
    multivslicedata *mvslicei;
    int j;

    mvslicei = multivsliceinfo + i;
    mvslicei->loaded=0;
    mvslicei->display=0;
    for(j=0;j<mvslicei->nvslices;j++){
      vslicedata *vd;

      vd = vsliceinfo + mvslicei->ivslices[j];
      if(vd->loaded==1)mvslicei->loaded++;
      if(vd->display==1)mvslicei->display++;
    }
    if(mvslicei->loaded>0&&mvslicei->loaded<mvslicei->nvslices){
      mvslicei->loaded=-1;
    }
    else if(mvslicei->loaded==mvslicei->nvslices){
      mvslicei->loaded=1;
    }
    if(mvslicei->display>0&&mvslicei->display<mvslicei->nvslices){
      mvslicei->display=-1;
    }
    else if(mvslicei->display==mvslicei->nvslices){
      mvslicei->display=1;
    }
  }


  npart5loaded=0;
  npartloaded=0;
  nevacloaded=0;
  for(i=0;i<npartinfo;i++){
    partdata *parti;

    parti = partinfo+i;
    if(parti->loaded==1&&parti->evac==0)npartloaded++;
    if(parti->loaded==1&&parti->evac==1)nevacloaded++;
    if(parti->loaded==1){
      npart5loaded++;
    }
  }

  nplot3dloaded=0;
  for(i=0;i<nplot3dinfo;i++){
    plot3ddata *plot3di;

    plot3di = plot3dinfo + i;
    if(plot3di->loaded==1)nplot3dloaded++;
  }

  nisoloaded=0;
  for(i=0;i<nisoinfo;i++){
    isodata *isoi;

    isoi = isoinfo + i;
    if(isoi->loaded==1)nisoloaded++;
  }

  npatchloaded=0;
  for(i=0;i<npatchinfo;i++){
    patchdata *patchi;

    patchi = patchinfo + i;
    if(patchi->loaded==1)npatchloaded++;
  }

#define CREATEMENU(menu,Menu) menu=glutCreateMenu(Menu);\
  if(nmenus<10000){\
    strcpy(menuinfo[nmenus].label,#Menu);\
    menuinfo[nmenus++].menuvar=menu;\
  }

  {
    for(i=0;i<nmenus;i++){
      menudata *menui;

      menui = menuinfo + i;

      if(menui->menuvar>0){
        glutDestroyMenu(menui->menuvar);
      }
    }
    nmenus=0;
  }
  if(nloadsubslicemenu>0){
    FREEMEMORY(loadsubslicemenu);
  }
  if(nloadsubmslicemenu>0){
    FREEMEMORY(loadsubmslicemenu);
  }
  if(nloadsubvslicemenu>0){
    FREEMEMORY(loadsubvslicemenu);
  }
  FREEMEMORY(particlepropshowsubmenu);
  if(nloadsubmvslicemenu>0){
    FREEMEMORY(loadsubmvslicemenu);
  }
  if(nloadsubplot3dmenu>0){
    FREEMEMORY(loadsubplot3dmenu);
  }
  if(unload==UNLOAD)return;



/* --------------------------------patch menu -------------------------- */
  if(npatchinfo>0){
    int npatchslice = 0;

    CREATEMENU(showpatchmenu,ShowPatchMenu);
    npatchloaded=0;
    {
      int local_do_threshold=0;
      char menulabel[1024];
      int ii;

      for(ii=0;ii<npatchinfo;ii++){
        patchdata *patchi;

        i = patchorderindex[ii];
        patchi = patchinfo+i;
        if(patchi->loaded==0)continue;
        npatchloaded++;
        if(patchi->slice == 1)npatchslice++;
        if(patchi->display==1&&patchi->type==ipatchtype){
          STRCPY(menulabel,"*");
          STRCAT(menulabel,patchi->menulabel);
        }
        else{
          STRCPY(menulabel,patchi->menulabel);
        }
        glutAddMenuEntry(menulabel,1000+i);
        if(activate_threshold==1){
          if(
            strncmp(patchi->label.shortlabel,"TEMP",4) == 0||
            strncmp(patchi->label.shortlabel,"temp",4) == 0
            ){
            local_do_threshold=1;
          }
        }
      }
      if(npatchslice>0){
        glutAddMenuEntry("Geometry slice data", DUMMYwallmenu);
        if(show_patch_solid==1){
          glutAddMenuEntry("  *solid", SOLIDpatchmenu);
        }
        else{
          glutAddMenuEntry("  solid", SOLIDpatchmenu);
        }
        if(show_patch_outline==1){
          glutAddMenuEntry("  *outline", OUTLINEpatchmenu);
        }
        else{
          glutAddMenuEntry("  outline", OUTLINEpatchmenu);
        }
        if(show_patch_verts==1){
          glutAddMenuEntry("  *points", POINTSpatchmenu);
        }
        else{
          glutAddMenuEntry("  points", POINTSpatchmenu);
        }
        glutAddMenuEntry("-", DUMMYwallmenu);
        if(show_patch_insolid==1){
          glutAddMenuEntry("  *in solid", INSOLIDpatchmenu);
        }
        else{
          glutAddMenuEntry("  in solid", INSOLIDpatchmenu);
        }
      }
      if(show_patch_ingas==1){
        glutAddMenuEntry("  *in gas", INGASpatchmenu);
      }
      else{
        glutAddMenuEntry("  in gas", INGASpatchmenu);
      }
      if(show_patch_incutcell == 1){
        glutAddMenuEntry("  *in cutcell", INCUTCELLpatchmenu);
      }
      else{
        glutAddMenuEntry("  in cutcell", INCUTCELLpatchmenu);
      }
      if(activate_threshold == 1 && local_do_threshold == 1){
        glutAddMenuEntry("-",DUMMYwallmenu);
        if(vis_threshold==1)glutAddMenuEntry("*char",SHOW_CHAR);
        if(vis_threshold==0)glutAddMenuEntry("char",SHOW_CHAR);
      }

    }
    if(npatchloaded>0){
      glutAddMenuEntry("-",DUMMYwallmenu);
    }
    if(npatchloaded>1){
      glutAddMenuEntry(_("Show all boundary files"),SHOWALL_BOUNDARY);
      glutAddMenuEntry(_("Hide all boundary files"),HIDEALL_BOUNDARY);
      glutAddMenuEntry("-",DUMMYwallmenu);
    }
    if(showexterior==1){
      glutAddMenuEntry(_("*Exterior"),EXTERIORwallmenu);
    }
    if(showexterior==0){
      glutAddMenuEntry(_("Exterior"),EXTERIORwallmenu);
    }
    if(visPatchType[INTERIORwall]==1){
      glutAddMenuEntry(_("*Interior"),INTERIORwallmenu);
    }
    if(visPatchType[INTERIORwall]==0){
      glutAddMenuEntry(_("Interior"),INTERIORwallmenu);
    }
    if(IsPatchType(FRONTwall)==1&&visPatchType[FRONTwall]==1){
      glutAddMenuEntry(_("*Front"),FRONTwallmenu);
    }
    if(IsPatchType(FRONTwall)==1&&visPatchType[FRONTwall]==0){
      glutAddMenuEntry(_("Front"),FRONTwallmenu);
    }
    if(IsPatchType(BACKwall)==1&& visPatchType[BACKwall]==1){
      glutAddMenuEntry(_("*Back"),BACKwallmenu);
    }
    if(IsPatchType(BACKwall)==1&& visPatchType[BACKwall]==0){
      glutAddMenuEntry(_("Back"),BACKwallmenu);
    }
    if(IsPatchType(LEFTwall)==1&& visPatchType[LEFTwall]==1){
      glutAddMenuEntry(_("*Left"),LEFTwallmenu);
    }
    if(IsPatchType(LEFTwall)==1&& visPatchType[LEFTwall]==0){
      glutAddMenuEntry(_("Left"),LEFTwallmenu);
    }
    if(IsPatchType(RIGHTwall)==1&&visPatchType[RIGHTwall]==1){
      glutAddMenuEntry(_("*Right"),RIGHTwallmenu);
    }
    if(IsPatchType(RIGHTwall)==1&&visPatchType[RIGHTwall]==0){
      glutAddMenuEntry(_("Right"),RIGHTwallmenu);
    }
    if(IsPatchType(UPwall)==1&&   visPatchType[UPwall]==1){
      glutAddMenuEntry(_("*Up"),UPwallmenu);
    }
    if(IsPatchType(UPwall)==1&&   visPatchType[UPwall]==0){
      glutAddMenuEntry(_("Up"),UPwallmenu);
    }
    if(IsPatchType(DOWNwall)==1&& visPatchType[DOWNwall]==1){
      glutAddMenuEntry(_("*Down"),DOWNwallmenu);
    }
    if(IsPatchType(DOWNwall)==1&& visPatchType[DOWNwall]==0){
      glutAddMenuEntry(_("Down"),DOWNwallmenu);
    }
  }

/* --------------------------------surface menu -------------------------- */

  CREATEMENU(immersedsurfacemenu,ImmersedMenu);
  if(show_faces_solid==1&&show_faces_outline==1){
    glutAddMenuEntry(_("*Solid and outline"),GEOMETRY_SOLIDOUTLINE);
  }
  else{
    glutAddMenuEntry(_("Solid and outline"),GEOMETRY_SOLIDOUTLINE);
  }
  if(show_faces_solid==1&&show_faces_outline==0){
    glutAddMenuEntry(_("*Solid only"),GEOMETRY_SOLID);
  }
  else{
    glutAddMenuEntry(_("Solid only"),GEOMETRY_SOLID);
  }
  if(show_faces_outline==1&&show_faces_solid==0){
    glutAddMenuEntry(_("*Outline only"),GEOMETRY_OUTLINE);
  }
  else{
    glutAddMenuEntry(_("Outline only"),GEOMETRY_OUTLINE);
  }
  if(show_faces_solid == 0 && show_faces_outline == 0){
    glutAddMenuEntry(_("*Hide"),GEOMETRY_HIDE);
  }
  else{
    glutAddMenuEntry(_("Hide"),GEOMETRY_HIDE);
  }

/* --------------------------------interior geometry menu -------------------------- */

  CREATEMENU(immersedinteriormenu,ImmersedMenu);
  if(have_volume==1){
    if(show_volumes_solid==1)glutAddMenuEntry(_("*Solid"),GEOMETRY_INTERIOR_SOLID);
    if(show_volumes_solid==0)glutAddMenuEntry(_("Solid"),GEOMETRY_INTERIOR_SOLID);
    if(show_volumes_outline==1)glutAddMenuEntry(_("*Outline"),GEOMETRY_INTERIOR_OUTLINE);
    if(show_volumes_outline==0)glutAddMenuEntry(_("Outline"),GEOMETRY_INTERIOR_OUTLINE);
    if(show_volumes_outline == 0 && show_volumes_solid == 0){
      glutAddMenuEntry(_("*Hide"),GEOMETRY_TETRA_HIDE);
    }
    else{
      glutAddMenuEntry(_("Hide"),GEOMETRY_TETRA_HIDE);
    }
  }

/* --------------------------------surface geometry menu -------------------------- */

  CREATEMENU(immersedmenu,ImmersedMenu);
  glutAddSubMenu(_("Faces"),immersedsurfacemenu);
  if(have_volume==1){
    glutAddSubMenu(_("Volumes"),immersedinteriormenu);
  }
  if(sort_geometry==1){
    glutAddMenuEntry(_("*Sort faces"), GEOMETRY_SORTFACES);
  }
  else{
    glutAddMenuEntry(_("Sort faces"), GEOMETRY_SORTFACES);
  }
  if(show_geom_normal == 1){
    glutAddMenuEntry(_("*Show normal"), GEOMETRY_SHOWNORMAL);
  }
  else{
    glutAddMenuEntry(_("Show normal"), GEOMETRY_SHOWNORMAL);
  }
  if(smooth_geom_normal==1){
    glutAddMenuEntry(_("*Smooth normal"), GEOMETRY_SMOOTHNORMAL);
  }
  else{
    glutAddMenuEntry(_("Smooth normal"), GEOMETRY_SMOOTHNORMAL);
  }
  if(ngeomdiaginfo>0){
    if(show_geometry_diagnostics == 1){
      glutAddMenuEntry(_("*Show geometry diagnostics"), GEOMETRY_SHOWDIAGNOSTICS);
    }
    else{
      glutAddMenuEntry(_("Show geometry diagnostics"), GEOMETRY_SHOWDIAGNOSTICS);
    }
  }
  if(hilight_skinny == 1){
    glutAddMenuEntry(_("*Hilight skinny triangles"), GEOMETRY_HILIGHTSKINNY);
  }
  else{
    glutAddMenuEntry(_("Hilight skinny triangles"), GEOMETRY_HILIGHTSKINNY);
  }
  if(show_faces_solid == 0 && show_faces_outline == 0 && show_volumes_solid == 0){
    glutAddMenuEntry(_("*Hide all"), GEOMETRY_HIDEALL);
  }
  else{
    glutAddMenuEntry(_("Hide all"), GEOMETRY_HIDEALL);
  }

/* --------------------------------blockage menu -------------------------- */

  CREATEMENU(blockagemenu,BlockageMenu);
  glutAddMenuEntry(_("View Method:"),MENU_DUMMY);
  if(visBlocks==visBLOCKAsInput||visBlocks==visBLOCKAsInputOutline){
    glutAddMenuEntry(_("   *Defined in input file"),visBLOCKAsInput);
  }
   else{
    glutAddMenuEntry(_("   Defined in input file"),visBLOCKAsInput);
  }
  if(visBlocks==visBLOCKNormal||visBlocks==visBLOCKSolidOutline){
    glutAddMenuEntry(_("   *Solid"),visBLOCKNormal);
    if(ntransparentblocks>0){
      if(visTransparentBlockage==1){
         glutAddMenuEntry(_("      *Transparent"),visBLOCKTransparent);
      }
      else{
         glutAddMenuEntry(_("      Transparent"),visBLOCKTransparent);
      }
    }
  }
  else{
    glutAddMenuEntry(_("   Solid"),visBLOCKNormal);
  }
  if(outline_state==OUTLINE_ONLY){
    glutAddMenuEntry(_("   *Outline Only"),visBLOCKOnlyOutline);
  }
  else{
    glutAddMenuEntry(_("   Outline Only"),visBLOCKOnlyOutline);
  }
  if(outline_state==OUTLINE_ADDED){
    glutAddMenuEntry(_("   *Outline Added"),visBLOCKAddOutline);
  }
  else{
    glutAddMenuEntry(_("   Outline Added"),visBLOCKAddOutline);
  }
  if(ncadgeom>0){
    if(viscadopaque==1){
      glutAddMenuEntry(_("   *Cad surface drawn opaque"),visCADOpaque);
    }
    else{
      glutAddMenuEntry(_("   Cad surface drawn opaque"),visCADOpaque);
    }
  }
  if(visBlocks==visBLOCKHide){
    glutAddMenuEntry(_("   *Hidden"),visBLOCKHide);
  }
  else{
    glutAddMenuEntry(_("   Hidden"),visBLOCKHide);
  }
  glutAddMenuEntry("-",MENU_DUMMY);
  glutAddMenuEntry(_(" Outline color:"),MENU_DUMMY);
  if(outline_color_flag==1){
    glutAddMenuEntry(_("   use blockage"),visBLOCKOutlineColor);
    glutAddMenuEntry(_("   *use foreground"),visBLOCKOutlineColor);
  }
  else{
    glutAddMenuEntry(_("   *use blockage"),visBLOCKOutlineColor);
    glutAddMenuEntry(_("   use foreground"),visBLOCKOutlineColor);
  }
  {
    int nblockprop=0;

    for(i=0;i<npropinfo;i++){
      propdata *propi;

      propi = propinfo + i;
      if(propi->inblockage==1)nblockprop++;
    }
    if(nblockprop>0){
      char propmenulabel[255];

      glutAddMenuEntry("-",MENU_DUMMY);
      glutAddMenuEntry(_("Show/Hide blockage types:"),MENU_DUMMY);
      for(i=0;i<npropinfo;i++){
        propdata *propi;

        propi = propinfo + i;
        if(propi->inblockage==1){
          strcpy(propmenulabel,"    ");
          if(propi->blockvis==1)strcat(propmenulabel,"*");
          strcat(propmenulabel,propi->label);
          glutAddMenuEntry(propmenulabel,-i-1);
        }
      }
    }
  }
  glutAddMenuEntry("-",MENU_DUMMY);
  glutAddMenuEntry(_("Locations:"),MENU_DUMMY);
  if(blocklocation==BLOCKlocation_grid){
    glutAddMenuEntry(_("   *Actual"),BLOCKlocation_grid);
  }
  else{
    glutAddMenuEntry(_("   Actual"),BLOCKlocation_grid);
  }
  if(blocklocation==BLOCKlocation_exact){
    glutAddMenuEntry(_("   *Requested"),BLOCKlocation_exact);
  }
  else{
    glutAddMenuEntry(_("   Requested"),BLOCKlocation_exact);
  }
  if(ncadgeom>0){
    if(blocklocation==BLOCKlocation_cad){
      glutAddMenuEntry(_("   *Cad"),BLOCKlocation_cad);
    }
    else{
      glutAddMenuEntry(_("   Cad"),BLOCKlocation_cad);
    }
    {
      cadgeomdata *cd;
      cadlookdata *cdi;
      int showtexturemenu;

      showtexturemenu=0;
      for(i=0;i<ncadgeom;i++){
        int j;

        cd = cadgeominfo + i;
        for(j=0;j<cd->ncadlookinfo;j++){
          cdi = cd->cadlookinfo+j;
          if(cdi->textureinfo.loaded==1){
            showtexturemenu=1;
            break;
          }
        }
        if(showtexturemenu==1)break;
      }
      if(showtexturemenu==1){
        if(visCadTextures==1){
          glutAddMenuEntry(_(" *Show CAD textures"),BLOCKtexture_cad);
        }
        else{
          glutAddMenuEntry(_(" Show CAD textures"),BLOCKtexture_cad);
        }
      }
    }
  }


/* --------------------------------level menu -------------------------- */

  if(nplot3dinfo>0){
    CREATEMENU(levelmenu,LevelMenu);
    for(i=1;i<nrgb-1;i++){
      if(colorlabeliso!=NULL){
        char *colorlabel;
        char levellabel2[256];

        colorlabel=&colorlabeliso[plotn-1][nrgb-2-i][0];
        strcpy(levellabel2,"");
        if(plotiso[plotn-1]==nrgb-2-i&&visiso==1){
          strcat(levellabel2,"*");
        }
        strcat(levellabel2,colorlabel);
        glutAddMenuEntry(levellabel2,nrgb-2-i);
      }
      else{
        char chari[4];

        if(plotiso[plotn]==i&&visiso==1){
          sprintf(chari,"*%i",i+1);
        }
        else{
          sprintf(chari,"%i",i+1);
        }
        glutAddMenuEntry(chari,i+1);
      }
    }
  }

/* --------------------------------static variable menu -------------------------- */

  if(nplot3dinfo>0){
    int n;

    CREATEMENU(staticvariablemenu,StaticVariableMenu);
    for(n=0;n<numplot3dvars;n++){
      char *p3label;

      p3label = plot3dinfo[0].label[n].shortlabel;
      if(plotn-1==n){
        char menulabel[1024];

        STRCPY(menulabel,"*");
        STRCAT(menulabel,p3label);
        glutAddMenuEntry(menulabel,n+1);
      }
      else{
        glutAddMenuEntry(p3label,n+1);
      }
    }
  }

/* --------------------------------iso variable menu -------------------------- */

  if(nplot3dinfo>0){
    int n;

    CREATEMENU(isovariablemenu,IsoVariableMenu);
    for(n=0;n<numplot3dvars;n++){
      char *p3label;

      p3label = plot3dinfo[0].label[n].shortlabel;
      if(plotn-1==n&&visiso==1){
        char menulabel[1024];

        STRCPY(menulabel,"*");
        STRCAT(menulabel,p3label);
        glutAddMenuEntry(menulabel,n+1);
      }
      else{
        glutAddMenuEntry(p3label,n+1);
      }
    }
  }

/* --------------------------------iso surface menu -------------------------- */
  if(nplot3dinfo>0){
    CREATEMENU(isosurfacetypemenu,IsoSurfaceTypeMenu);
    if(p3dsurfacesmooth==1&&p3dsurfacetype==SURFACE_SOLID){
      glutAddMenuEntry(_("*Smooth"),MENU_SURFACE_SMOOTH);
    }
     else{
       glutAddMenuEntry(_("Smooth"),MENU_SURFACE_SMOOTH);
     }
     if(p3dsurfacesmooth==0&&p3dsurfacetype==SURFACE_SOLID){
       glutAddMenuEntry(_("*Facets"),MENU_SURFACE_FACET);
     }
    else{
      glutAddMenuEntry(_("Facets"),MENU_SURFACE_FACET);
    }
    if(p3dsurfacetype==SURFACE_OUTLINE)glutAddMenuEntry(_("*Triangles"),SURFACE_OUTLINE);
    if(p3dsurfacetype!=SURFACE_OUTLINE)glutAddMenuEntry(_("Triangles"),SURFACE_OUTLINE);
    if(p3dsurfacetype == SURFACE_POINTS)glutAddMenuEntry(_("*Points"), SURFACE_POINTS);
    if(p3dsurfacetype != SURFACE_POINTS)glutAddMenuEntry(_("Points"), SURFACE_POINTS);

    CREATEMENU(isosurfacemenu,IsoSurfaceMenu);
    glutAddSubMenu(_("Solution variable"),isovariablemenu);
    glutAddSubMenu(_("Solution value"),levelmenu);
    glutAddSubMenu(_("Surface type"),isosurfacetypemenu);
    glutAddMenuEntry(_("Hide"),1);
  }

/* --------------------------------vector skip menu -------------------------- */

  if(nplot3dinfo>0){
    CREATEMENU(vectorskipmenu,VectorSkipMenu);
    if(visVector==1)glutAddMenuEntry(_("*Show"),MENU_VECTOR_SHOW);
    if(visVector!=1)glutAddMenuEntry(_("Show"),MENU_VECTOR_SHOW);
    glutAddMenuEntry(_("Frequency:"),-1);
    if(vectorskip==1)glutAddMenuEntry(_("*All"),1);
    if(vectorskip!=1)glutAddMenuEntry(_("All"),1);
    if(vectorskip==2)glutAddMenuEntry(_("*Every 2nd"),2);
    if(vectorskip!=2)glutAddMenuEntry(_("Every 2nd"),2);
    if(vectorskip==3)glutAddMenuEntry(_("*Every 3rd"),3);
    if(vectorskip!=3)glutAddMenuEntry(_("Every 3rd"),3);
    if(vectorskip==4)glutAddMenuEntry(_("*Every 4th"),4);
    if(vectorskip!=4)glutAddMenuEntry(_("Every 4th"),4);
  }

  if(ntextures_loaded_used>0){
    int ntextures_used;

    CREATEMENU(textureshowmenu,TextureShowMenu);
    ntextures_used=0;
    for(i=0;i<ntextures;i++){
      texturedata *texti;
      char menulabel[1024];

      texti = textureinfo + i;
      if(texti->loaded==0||texti->used==0)continue;
      ntextures_used++;
      if(texti->display==1){
        STRCPY(menulabel,"*");
        STRCAT(menulabel,texti->file);
        glutAddMenuEntry(menulabel,i);
      }
      else{
        STRCPY(menulabel,texti->file);
        glutAddMenuEntry(menulabel,i);
      }
    }
    if(ntextures_used>1){
      glutAddMenuEntry("-",MENU_DUMMY);
      glutAddMenuEntry(_("Show all"),MENU_TEXTURE_SHOWALL);
      glutAddMenuEntry(_("Hide all"),MENU_TEXTURE_HIDEALL);
    }
  }

/* --------------------------------Plot3d Show menu -------------------------- */
  if(nplot3dinfo>0){
    CREATEMENU(staticslicemenu,Plot3DShowMenu);
    glutAddSubMenu(_("Solution variable"),staticvariablemenu);
    if(visz_all==1)glutAddMenuEntry(_("*xy plane"), MENU_PLOT3D_Z);
    if(visz_all==0)glutAddMenuEntry(_("xy plane"), MENU_PLOT3D_Z);
    if(visy_all==1)glutAddMenuEntry(_("*xz plane"), MENU_PLOT3D_Y);
    if(visy_all==0)glutAddMenuEntry(_("xz plane"), MENU_PLOT3D_Y);
    if(visx_all==1)glutAddMenuEntry(_("*yz plane"), MENU_PLOT3D_X);
    if(visx_all==0)glutAddMenuEntry(_("yz plane"), MENU_PLOT3D_X);
    if(vectorspresent==1)glutAddSubMenu(_("Flow vectors"),vectorskipmenu);
    if(contour_type==SHADED_CONTOURS){
      glutAddMenuEntry(_("*Continuous contours"), MENU_PLOT3D_CONT);
    }
    if(contour_type!=SHADED_CONTOURS){
      glutAddMenuEntry(_("Continuous contours"), MENU_PLOT3D_CONT);
    }
    glutAddMenuEntry(_("Show all planes"), MENU_PLOT3D_SHOWALL);
    glutAddMenuEntry(_("Hide all planes"), MENU_PLOT3D_HIDEALL);

    CREATEMENU(plot3dshowmenu,Plot3DShowMenu);
    if(nplot3dloaded>0){
      int ii;

      for(ii=0;ii<nplot3dinfo;ii++){
        plot3ddata *plot3di;
        char menulabel[1024];

        i=plot3dorderindex[ii];
        plot3di = plot3dinfo + i;
        if(ii==0){
          glutAddMenuEntry(plot3di->longlabel, MENU_PLOT3D_DUMMY);
        }
        else{
          if(strcmp(plot3di->longlabel, plot3dinfo[plot3dorderindex[ii - 1]].longlabel) != 0){
            glutAddMenuEntry(plot3di->longlabel, MENU_PLOT3D_DUMMY);
          }
        }
        if(plot3di->loaded==0)continue;
        if(plotstate==STATIC_PLOTS&&plot3di->display==1){
          STRCPY(menulabel,"*");
          STRCAT(menulabel,plot3di->menulabel);
        }
        else{
          STRCPY(menulabel,plot3di->menulabel);
        }
        glutAddMenuEntry(menulabel,1000+i);
      }
      if(nplot3dloaded>1){
        glutAddMenuEntry("-",MENU_PLOT3D_DUMMY);
        glutAddMenuEntry(_("Show all PLOT3D files"),SHOWALL_PLOT3D);
        glutAddMenuEntry(_("Hide all PLOT3D files"),HIDEALL_PLOT3D);
      }
      glutAddMenuEntry("-",MENU_PLOT3D_DUMMY);
    }
    glutAddSubMenu(_("2D contours"),staticslicemenu);
    if(cache_qdata==1){
      glutAddSubMenu(_("3D contours"),isosurfacemenu);
    }

  }

/* --------------------------------grid slice menu -------------------------- */

  CREATEMENU(gridslicemenu,GridSliceMenu);
  if(visGrid==GridnoProbe||visGrid==GridProbe){
    glutAddMenuEntry(_("*show grid"),GRID_grid);
  }
  else{
    glutAddMenuEntry(_("show grid"),GRID_grid);
  }
  if(visGrid==GridProbe||visGrid==noGridProbe){
    glutAddMenuEntry(_("*show grid location"),GRID_probe);
  }
  else{
    glutAddMenuEntry(_("show grid location"),GRID_probe);
  }
  glutAddMenuEntry("-",MENU_DUMMY);
  if(visz_all==1){
    glutAddMenuEntry(_("*xy plane"),GRID_xy);
  }
  else{
    glutAddMenuEntry(_("xy plane"),GRID_xy);
  }
  if(visy_all==1){
    glutAddMenuEntry(_("*xz plane"),GRID_xz);
  }
  else{
    glutAddMenuEntry(_("xz plane"),GRID_xz);
  }
  if(visx_all==1){
    glutAddMenuEntry(_("*yz plane"),GRID_yz);
  }
  else{
    glutAddMenuEntry(_("yz plane"),GRID_yz);
  }
  if(visx_all==0||visy_all==0||visz_all==0){
    glutAddMenuEntry(_("Show all planes"),GRID_showall);
  }
  if(visx_all==1||visy_all==1||visz_all==1){
    glutAddMenuEntry(_("Hide all planes"),GRID_hideall);
  }

  CREATEMENU(circularventmenu,VentMenu);
  if(visCircularVents==VENT_CIRCLE){
    glutAddMenuEntry(_("*As circle"), MENU_VENT_CIRCLE);
    glutAddMenuEntry(_("As rectangle"), MENU_VENT_RECTANGLE);
    glutAddMenuEntry(_("Hide"), MENU_VENT_CIRCLEHIDE);
  }
  if(visCircularVents==VENT_RECTANGLE){
    glutAddMenuEntry(_("As circle"), MENU_VENT_CIRCLE);
    glutAddMenuEntry(_("*As rectangle"), MENU_VENT_RECTANGLE);
    glutAddMenuEntry(_("Hide"), MENU_VENT_CIRCLEHIDE);
  }
  if(visCircularVents==VENT_HIDE){
    glutAddMenuEntry(_("As circle"), MENU_VENT_CIRCLE);
    glutAddMenuEntry(_("As rectangle"), MENU_VENT_RECTANGLE);
    glutAddMenuEntry(_("*Hide"), MENU_VENT_CIRCLEHIDE);
  }
  glutAddMenuEntry("-",MENU_DUMMY2);
  if(circle_outline == 1)glutAddMenuEntry("*Outline", MENU_VENT_CIRCLEOUTLINE);
  if(circle_outline == 0)glutAddMenuEntry("Outline", MENU_VENT_CIRCLEOUTLINE);

/* --------------------------------vent menu -------------------------- */

  CREATEMENU(ventmenu,VentMenu);
  if(GetNTotalVents()>0){
    if(nopenvents>0){
      if(visOpenVents == 1)glutAddMenuEntry(_("*Open"), MENU_VENT_OPEN);
      if(visOpenVents == 0)glutAddMenuEntry(_("Open"), MENU_VENT_OPEN);
    }
    if(ndummyvents>0){
      if(visDummyVents == 1)glutAddMenuEntry(_("*Exterior"), MENU_VENT_EXTERIOR);
      if(visDummyVents == 0)glutAddMenuEntry(_("Exterior"), MENU_VENT_EXTERIOR);
    }
    if(ncvents>0){
      if(visCircularVents!=VENT_HIDE)glutAddSubMenu(_("*Circular"),circularventmenu);
      if(visCircularVents==VENT_HIDE)glutAddSubMenu(_("Circular"),circularventmenu);
    }
    if(GetNTotalVents()>nopenvents+ndummyvents){
      if(visOtherVents == 1)glutAddMenuEntry(_("*Other"), MENU_VENT_OTHER);
      if(visOtherVents == 0)glutAddMenuEntry(_("Other"), MENU_VENT_OTHER);
    }
    if(visOpenVents==1&&visDummyVents==1&&visOtherVents==1){
      glutAddMenuEntry(_("*Show all"),SHOW_ALL_VENTS);
    }
    else{
      glutAddMenuEntry(_("Show all"),SHOW_ALL_VENTS);
    }
    if(visOpenVents==0&&visDummyVents==0&&visOtherVents==0){
      glutAddMenuEntry(_("*Hide all"),HIDE_ALL_VENTS);
    }
    else{
      glutAddMenuEntry(_("Hide all"),HIDE_ALL_VENTS);
    }
    glutAddMenuEntry("-",MENU_DUMMY2);
    if(nopenvents_nonoutline>0){
      if(visOpenVentsAsOutline == 1)glutAddMenuEntry(_("*Open vents as outlines"), MENU_VENT_OUTLINE);
      if(visOpenVentsAsOutline == 0)glutAddMenuEntry(_("Open vents as outlines"), MENU_VENT_OUTLINE);
    }
    if(have_vents_int==1){
      if(show_bothsides_int == 1)glutAddMenuEntry(_("*Two sided (interior)"), MENU_VENT_TWOINTERIOR);
      if(show_bothsides_int == 0)glutAddMenuEntry(_("Two sided (interior)"), MENU_VENT_TWOINTERIOR);
    }
    if(show_bothsides_ext == 1)glutAddMenuEntry(_("*Two sided (exterior)"), MENU_VENT_TWOEXTERIOR);
    if(show_bothsides_ext == 0)glutAddMenuEntry(_("Two sided (exterior)"), MENU_VENT_TWOEXTERIOR);
    if(nvent_transparent>0){
      if(show_transparent_vents == 1)glutAddMenuEntry(_("*Transparent"), MENU_VENT_TRANSPARENT);
      if(show_transparent_vents == 0)glutAddMenuEntry(_("Transparent"), MENU_VENT_TRANSPARENT);
    }
  }

/* --------------------------------terrain_showmenu -------------------------- */

  CREATEMENU(terrain_showmenu,GeometryMenu);
  if(visTerrainType==TERRAIN_3D)glutAddMenuEntry(_("*3D surface"),17+TERRAIN_3D);
  if(visTerrainType!=TERRAIN_3D)glutAddMenuEntry(_("3D surface"),17+TERRAIN_3D);
  if(visTerrainType==TERRAIN_2D_STEPPED)glutAddMenuEntry(_("*2D stepped"),17+TERRAIN_2D_STEPPED);
  if(visTerrainType!=TERRAIN_2D_STEPPED)glutAddMenuEntry(_("2D stepped"),17+TERRAIN_2D_STEPPED);
  if(visTerrainType==TERRAIN_2D_LINE)glutAddMenuEntry(_("*2D lines"),17+TERRAIN_2D_LINE);
  if(visTerrainType!=TERRAIN_2D_LINE)glutAddMenuEntry(_("2D lines"),17+TERRAIN_2D_LINE);
  if(terrain_texture!=NULL&&terrain_texture->loaded==1){
    if(visTerrainType==TERRAIN_3D_MAP)glutAddMenuEntry(_("*Image"),17+TERRAIN_3D_MAP);
    if(visTerrainType!=TERRAIN_3D_MAP)glutAddMenuEntry(_("Image"),17+TERRAIN_3D_MAP);
  }
  if(visTerrainType==TERRAIN_HIDDEN)glutAddMenuEntry(_("*Hidden"),17+TERRAIN_HIDDEN);
  if(visTerrainType!=TERRAIN_HIDDEN)glutAddMenuEntry(_("Hidden"),17+TERRAIN_HIDDEN);

  if(nobject_defs>0){
    int multiprop;

    multiprop=0;
    for(i=0;i<npropinfo;i++){
      propdata *propi;

      propi = propinfo + i;
      if(propi->nsmokeview_ids>1)multiprop=1;
    }
    if(multiprop==1){
      for(i=0;i<npropinfo;i++){
        propdata *propi;

        propi = propinfo + i;
        CREATEMENU(propi->menu_id,PropMenu);
        if(propi->nsmokeview_ids>1){
          int jj;
          char menulabel[1024];

          for(jj=0;jj<propi->nsmokeview_ids;jj++){
            strcpy(menulabel,"");
            if(propi->smokeview_ids[jj]==propi->smokeview_id){
              strcat(menulabel,"*");
            }
            strcat(menulabel,propi->smokeview_ids[jj]);
            glutAddMenuEntry(menulabel,jj*npropinfo+i);
          }
        }
      }
      CREATEMENU(propmenu,PropMenu);
      for(i=0;i<npropinfo;i++){
        propdata *propi;

        propi = propinfo + i;
        if(propi->nsmokeview_ids>1){
          glutAddSubMenu(propi->label,propi->menu_id);
        }
      }
    }

    CREATEMENU(spheresegmentmenu,ShowObjectsMenu);
    if(device_sphere_segments==6){
      glutAddMenuEntry("   *6",-6);
    }
    else{
      glutAddMenuEntry("   6",-6);
    }
    if(device_sphere_segments==12){
      glutAddMenuEntry("   *12",-12);
    }
    else{
      glutAddMenuEntry("   12",-12);
    }
    if(device_sphere_segments==24){
      glutAddMenuEntry("   *24",-24);
    }
    else{
      glutAddMenuEntry("   24",-24);
    }
    if(device_sphere_segments==48){
      glutAddMenuEntry("   *48",-48);
    }
    else{
      glutAddMenuEntry("   48",-48);
    }
  }

  if(nobject_defs>0){
    CREATEMENU(showobjectsmenu,ShowObjectsMenu);
    for(i=0;i<nobject_defs;i++){
      sv_object *obj_typei;

      obj_typei = object_defs[i];
      if(obj_typei->used_by_device==1){
        char obj_menu[256];

        strcpy(obj_menu,"");
        if(obj_typei->visible==1){
          strcat(obj_menu,"*");
        }
        strcat(obj_menu,obj_typei->label);
        glutAddMenuEntry(obj_menu,i);
      }
    }
    if(have_missing_objects == 1&&isZoneFireModel==0){
      glutAddMenuEntry("-", MENU_DUMMY);
      if(show_missing_objects==1)glutAddMenuEntry("*undefined",OBJECT_MISSING);
      if(show_missing_objects == 0)glutAddMenuEntry("undefined",OBJECT_MISSING);
    }
    glutAddMenuEntry("-",MENU_DUMMY);
    if(ndeviceinfo>0){
      if(select_device==1){
        glutAddMenuEntry(_("*Select"),OBJECT_SELECT);
      }
      else{
        glutAddMenuEntry(_("Select"),OBJECT_SELECT);
      }
    }
    if(object_outlines==0)glutAddMenuEntry(_("Outline"),OBJECT_OUTLINE);
    if(object_outlines==1)glutAddMenuEntry(_("*Outline"),OBJECT_OUTLINE);
    glutAddMenuEntry(_("Show all"),OBJECT_SHOWALL);
    glutAddMenuEntry(_("Hide all"),OBJECT_HIDEALL);
    if(show_device_orientation==1){
      glutAddMenuEntry(_("*Show orientation"),OBJECT_ORIENTATION);
    }
    else{
      glutAddMenuEntry(_("Show orientation"),OBJECT_ORIENTATION);
    }
    glutAddMenuEntry("-",MENU_DUMMY);
    glutAddSubMenu(_("Segments"),spheresegmentmenu);

  }

/* --------------------------------geometry menu -------------------------- */

  CREATEMENU(geometrymenu,GeometryMenu);
  if(ntotal_blockages>0)glutAddSubMenu(_("Obstacles"),blockagemenu);
  if(ngeominfo>0)glutAddSubMenu(_("Immersed"),immersedmenu);
  if(GetNumActiveDevices()>0||ncvents>0){
    glutAddSubMenu(_("Objects"),showobjectsmenu);
  }
  if(nterraininfo>0){
    glutAddSubMenu(_("Terrain"),terrain_showmenu);
  }
  if(GetNTotalVents()>0)glutAddSubMenu(_("Surfaces"), ventmenu);
  if(nzvents > 0){
    if(visVents == 1){
      glutAddMenuEntry(_("*Vents"), GEOM_Vents);
    }
    else{
      glutAddMenuEntry(_("Vents"), GEOM_Vents);
    }
  }
  if(ntotal_blockages>0 || isZoneFireModel == 1){
    glutAddSubMenu(_("Grid"),gridslicemenu);
  }
  if(isZoneFireModel==0){
    if(visFrame==1)glutAddMenuEntry(_("*Outline"), GEOM_Outline);
    if(visFrame==0)glutAddMenuEntry(_("Outline"), GEOM_Outline);
  }
  else{
    visFrame=0;
  }
#ifdef _DEBUG
  if(show_triangle_count==1)glutAddMenuEntry(_("*Triangle count"), GEOM_TriangleCount);
  if(show_triangle_count==0)glutAddMenuEntry(_("Triangle count"), GEOM_TriangleCount);
#endif
  glutAddMenuEntry(_("Show all"), GEOM_ShowAll);
  glutAddMenuEntry(_("Hide all"), GEOM_HideAll);

/* --------------------------------label menu -------------------------- */

  CREATEMENU(labelmenu, LabelMenu);
  if(visColorbar==1)glutAddMenuEntry(_("*Colorbar"),MENU_LABEL_colorbar);
  if(visColorbar==0)glutAddMenuEntry(_("Colorbar"),MENU_LABEL_colorbar);
  if(visTimebar==1)glutAddMenuEntry(_("*Time bar"),MENU_LABEL_timebar);
  if(visTimebar==0)glutAddMenuEntry(_("Time bar"),MENU_LABEL_timebar);
  if(visTitle == 1)glutAddMenuEntry(_("*Title"), MENU_LABEL_title);
  if(visTitle == 0)glutAddMenuEntry(_("Title"), MENU_LABEL_title);
  glutAddMenuEntry("-", MENU_DUMMY);

  if(visaxislabels == 1)glutAddMenuEntry(_("*Axis"), MENU_LABEL_axis);
  if(visaxislabels == 0)glutAddMenuEntry(_("Axis"), MENU_LABEL_axis);
  if(have_northangle==1){
    if(vis_northangle==1)glutAddMenuEntry(_("*North"), MENU_LABEL_northangle);
    if(vis_northangle==0)glutAddMenuEntry(_("North"), MENU_LABEL_northangle);
  }
  if(ntickinfo>0){
    if(visFDSticks == 0)glutAddMenuEntry(_("FDS generated ticks"), MENU_LABEL_fdsticks);
    if(visFDSticks == 1)glutAddMenuEntry(_("*FDS generated ticks"), MENU_LABEL_fdsticks);
  }
  if(visFramelabel == 1)glutAddMenuEntry(_("*Frame"), MENU_LABEL_framelabel);
  if(visFramelabel == 0)glutAddMenuEntry(_("Frame"), MENU_LABEL_framelabel);
  if(visFramerate == 1)glutAddMenuEntry(_("*Frame rate"), MENU_LABEL_framerate);
  if(visFramerate == 0)glutAddMenuEntry(_("Frame rate"), MENU_LABEL_framerate);
  if(visgridloc == 1)glutAddMenuEntry(_("*Grid locations"), MENU_LABEL_grid);
  if(visgridloc == 0)glutAddMenuEntry(_("Grid locations"), MENU_LABEL_grid);
  if(show_hrrcutoff_active == 1){
    if(show_hrrcutoff == 1)glutAddMenuEntry(_("*HRRPUV cutoff"), MENU_LABEL_hrrcutoff);
    if(show_hrrcutoff == 0)glutAddMenuEntry(_("HRRPUV cutoff"), MENU_LABEL_hrrcutoff);
  }
  if(hrrinfo != NULL){
    if(visHRRlabel == 1)glutAddMenuEntry(_("*HRR"), MENU_LABEL_hrr);
    if(visHRRlabel == 0)glutAddMenuEntry(_("HRR"), MENU_LABEL_hrr);
  }
#ifdef pp_memstatus
  if(visAvailmemory == 1)glutAddMenuEntry(_("*Memory load"), MENU_LABEL_memload);
  if(visAvailmemory == 0)glutAddMenuEntry(_("Memory load"), MENU_LABEL_memload);
#endif
#ifdef pp_MEMDEBUG
  if(visUsagememory == 1)glutAddMenuEntry(_("*Memory usage"), MENU_LABEL_memusage);
  if(visUsagememory == 0)glutAddMenuEntry(_("Memory usage"), MENU_LABEL_memusage);
#endif
  if(visMeshlabel == 1)glutAddMenuEntry(_("*Mesh"), MENU_LABEL_meshlabel);
  if(visMeshlabel == 0)glutAddMenuEntry(_("Mesh"), MENU_LABEL_meshlabel);
  if(vis_slice_average == 1)glutAddMenuEntry(_("*Slice average"), MENU_LABEL_sliceaverage);
  if(vis_slice_average == 0)glutAddMenuEntry(_("Slice average"), MENU_LABEL_sliceaverage);
  if(LabelGetNUserLabels() > 0){
    if(visLabels == 1)glutAddMenuEntry(_("*Text labels"), MENU_LABEL_textlabels);
    if(visLabels == 0)glutAddMenuEntry(_("Text labels"), MENU_LABEL_textlabels);
  }
  if(visTimelabel == 1)glutAddMenuEntry(_("*Time"), MENU_LABEL_timelabel);
  if(visTimelabel == 0)glutAddMenuEntry(_("Time"), MENU_LABEL_timelabel);
  if(visUSERticks == 1)glutAddMenuEntry(_("*User settable ticks"), MENU_LABEL_userticks);
  if(visUSERticks == 0)glutAddMenuEntry(_("User settable ticks"), MENU_LABEL_userticks);
  if(gversion == 1)glutAddMenuEntry(_("*Version info"), MENU_LABEL_gversion);
  if(gversion == 0)glutAddMenuEntry(_("Version info"), MENU_LABEL_gversion);

  glutAddMenuEntry("-", MENU_DUMMY);
  glutAddMenuEntry(_("Show all"), MENU_LABEL_ShowAll);
  glutAddMenuEntry(_("Hide all"), MENU_LABEL_HideAll);

/* --------------------------------rotate type menu -------------------------- */

  CREATEMENU(rotatetypemenu,RotateTypeMenu);
  glutAddMenuEntry("Scene centered:",MENU_DUMMY);
  switch(rotation_type){
  case EYE_CENTERED:
    glutAddMenuEntry("  2 axis",ROTATION_2AXIS);
    glutAddMenuEntry("  Level (1 axis)",ROTATION_1AXIS);
    glutAddMenuEntry("  3 axis",ROTATION_3AXIS);
    glutAddMenuEntry("*Eye centered",EYE_CENTERED);
    break;
  case ROTATION_2AXIS:
    glutAddMenuEntry("  *2 axis",ROTATION_2AXIS);
    glutAddMenuEntry("  Level (1 axis)",ROTATION_1AXIS);
    glutAddMenuEntry("  3 axis",ROTATION_3AXIS);
    glutAddMenuEntry("Eye centered",EYE_CENTERED);
    break;
  case ROTATION_1AXIS:
    glutAddMenuEntry("  2 axis",ROTATION_2AXIS);
    glutAddMenuEntry("  *Level (1 axis)",ROTATION_1AXIS);
    glutAddMenuEntry("  3 axis",ROTATION_3AXIS);
    glutAddMenuEntry("Eye centered",EYE_CENTERED);
    break;
  case ROTATION_3AXIS:
    glutAddMenuEntry("  2 axis",ROTATION_2AXIS);
    glutAddMenuEntry("  Level (1 axis)",ROTATION_1AXIS);
    glutAddMenuEntry("  *3 axis",ROTATION_3AXIS);
    glutAddMenuEntry("Eye centered",EYE_CENTERED);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }

/* --------------------------------zone show menu -------------------------- */

  if(nzoneinfo>0&&(ReadZoneFile==1||nzvents>0)){
    CREATEMENU(zoneshowmenu,ZoneShowMenu);
    glutAddMenuEntry(_("Layers"),MENU_DUMMY);
    glutAddMenuEntry(_("   Representation:"),MENU_DUMMY);
    if(visZone==1&&zonecolortype==ZONETEMP_COLOR){
      glutAddMenuEntry(_("      *Temperature"), MENU_ZONE_2DTEMP);
#ifdef pp_HAZARD
      glutAddMenuEntry(_("      Hazard"), MENU_ZONE_2DHAZARD);
#endif
      glutAddMenuEntry(_("      Smoke"), 7);
    }
    else if(visZone==1&&zonecolortype==ZONEHAZARD_COLOR){
      glutAddMenuEntry(_("      Temperature"), MENU_ZONE_2DTEMP);
#ifdef pp_HAZARD
      glutAddMenuEntry(_("      *Hazard"), MENU_ZONE_2DHAZARD);
#endif
      glutAddMenuEntry(_("      Smoke"), MENU_ZONE_3DSMOKE);
    }
    else{
      glutAddMenuEntry(_("      Temperature"), MENU_ZONE_2DTEMP);
#ifdef pp_HAZARD
      glutAddMenuEntry(_("      Hazard"), MENU_ZONE_2DHAZARD);
#endif
      glutAddMenuEntry(_("      *Smoke"), MENU_ZONE_3DSMOKE);
    }
    glutAddMenuEntry(_("   Orientation:"), MENU_DUMMY);
    if(visZone==1){
      if(zonecolortype==ZONESMOKE_COLOR){
        glutAddMenuEntry(_("      Horizontal"), MENU_ZONE_HORIZONTAL);
        glutAddMenuEntry(_("      Vertical"), MENU_ZONE_VERTICAL);
      }
      else{
        if(visHZone==1)glutAddMenuEntry(_("      *Horizontal"), MENU_ZONE_HORIZONTAL);
        if(visHZone==0)glutAddMenuEntry(_("      Horizontal"), MENU_ZONE_HORIZONTAL);
        if(visVZone==1)glutAddMenuEntry(_("      *Vertical"), MENU_ZONE_VERTICAL);
        if(visVZone==0)glutAddMenuEntry(_("      Vertical"), MENU_ZONE_VERTICAL);
      }
    }
    else{
      glutAddMenuEntry(_("      Horizontal"), MENU_ZONE_HORIZONTAL);
      glutAddMenuEntry(_("      Vertical"), MENU_ZONE_VERTICAL);
    }
    if(visZone==0){
      glutAddMenuEntry(_("   *Hide"), MENU_ZONE_LAYERHIDE);
    }
    else{
      glutAddMenuEntry(_("   Hide"), MENU_ZONE_LAYERHIDE);
    }
    if(nzvents>0){
      if(visVentFlow==1){
        glutAddMenuEntry(_("*Vent flow"), MENU_ZONE_VENTS);
      }
      else{
        glutAddMenuEntry(_("Vent flow"), MENU_ZONE_VENTS);
      }
      if(nzhvents>0){
        if(visVentHFlow==1){
          glutAddMenuEntry(_("   *Horizontal"), MENU_ZONE_HVENTS);
        }
        else{
          glutAddMenuEntry(_("   Horizontal"), MENU_ZONE_HVENTS);
        }
        if(have_ventslab_flow==0){
          if(visventprofile==1){
            glutAddMenuEntry(_("      *velocity profile"), MENU_ZONE_VENT_PROFILE);
          }
          else{
            glutAddMenuEntry(_("      velocity profile"), MENU_ZONE_VENT_PROFILE);
          }
        }
        else{
          if(visventslab==1){
            glutAddMenuEntry(_("      *mass slab"), MENU_ZONE_VENT_SLAB);
          }
          else{
            glutAddMenuEntry(_("      mass slab"), MENU_ZONE_VENT_SLAB);
          }
          if(visventprofile==1){
            glutAddMenuEntry(_("      *velocity profile"), MENU_ZONE_VENT_PROFILE);
          }
          else{
            glutAddMenuEntry(_("      velocity profile"), MENU_ZONE_VENT_PROFILE);
          }
        }
      }
      if(nzvvents>0){
        if(visVentVFlow==1){
          glutAddMenuEntry(_("   *Vertical"), MENU_ZONE_VVENTS);
        }
        else{
          glutAddMenuEntry(_("   Vertical"), MENU_ZONE_VVENTS);
        }
      }
      if(nzmvents>0){
        if(visVentMFlow==1){
          glutAddMenuEntry(_("   *Mechancial"), MENU_ZONE_MVENTS);
        }
        else{
          glutAddMenuEntry(_("   Mechancial"), MENU_ZONE_MVENTS);
        }
      }
    }
    if(nfires>0){
      if(viszonefire==1){
        glutAddMenuEntry(_("*Fires"), MENU_ZONE_FIRES);
      }
      else{
        glutAddMenuEntry(_("Fires"), MENU_ZONE_FIRES);
      }
    }
  }

  /* --------------------------------particle class show menu -------------------------- */

  if(npartclassinfo>0){
    int ntypes;

    CREATEMENU(particlestreakshowmenu,ParticleStreakShowMenu);
    {
      int iii;
      char streaklabel[1024];

      streak_rvalue[nstreak_rvalue-1]=tmax_part;
      for(iii=0;iii<nstreak_rvalue;iii++){
        if(iii==streak_index){
          sprintf(streaklabel,"*%f",streak_rvalue[iii]);
        }
        else{
          sprintf(streaklabel,"%f",streak_rvalue[iii]);
        }
        TrimZeros(streaklabel);
        strcat(streaklabel," s");
        glutAddMenuEntry(streaklabel,iii);
      }
    }
    glutAddMenuEntry("-",MENU_DUMMY2);
    if(showstreakhead==1){
      glutAddMenuEntry(_("*Particle head"),MENU_STREAK_HEAD);
    }
    else{
      glutAddMenuEntry(_("Particle head"),MENU_STREAK_HEAD);
    }
    glutAddMenuEntry(_("Hide"),MENU_STREAK_HIDE);

// allocate memory for particle property sub-menus

    if(npart5prop*npartclassinfo>0){
      NewMemory((void **)&particlepropshowsubmenu,npart5prop*npartclassinfo*sizeof(int));
    }

      ntypes=0;
      for(i=0;i<npart5prop;i++){
        partpropdata *propi;
        int j;

        propi = part5propinfo + i;
        if(propi->display==0)continue;
        for(j=0;j<npartclassinfo;j++){
          partclassdata *partclassj;
          char menulabel[1024];

          if(propi->class_present[j]==0)continue;
          partclassj = partclassinfo + j;
          CREATEMENU(particlepropshowsubmenu[ntypes],ParticlePropShowMenu);
          ntypes++;
          if(propi->class_vis[j]==1){
            strcpy(menulabel,"  *");
          }
          else{
            strcpy(menulabel,"  ");
          }
          strcat(menulabel,partclassj->name);
          if(partclassj->col_diameter>=0||partclassj->col_length>=0||partclassj->device_name!=NULL||
             (partclassj->prop!=NULL&&partclassj->prop->smokeview_id!=NULL)||
             (partclassj->col_u_vel>=0&&partclassj->col_v_vel>=0&&partclassj->col_w_vel>=0)
            ){
            if(propi->class_vis[j]==1){
              strcpy(menulabel,_("using:"));
            }
            else{
              strcpy(menulabel,_("using:"));
            }
            glutAddMenuEntry(menulabel,-10-5*j);
            if(partclassj->kind!=HUMANS){
              if(partclassj->vis_type==PART_POINTS){
                glutAddMenuEntry(_("    *points"),-10-5*j-PART_POINTS);
              }
              else{
                glutAddMenuEntry(_("    points"),-10-5*j-PART_POINTS);
              }
              if(partclassj->col_diameter>=0||partclassj->device_name!=NULL){
                if(partclassj->vis_type==PART_SPHERES){
                  glutAddMenuEntry(_("    *spheres"),-10-5*j-PART_SPHERES);
                }
                else{
                  glutAddMenuEntry(_("    spheres"),-10-5*j-PART_SPHERES);
                }
              }
              if(partclassj->col_length>=0||partclassj->device_name!=NULL||
                (partclassj->col_u_vel>=0&&partclassj->col_v_vel>=0&&partclassj->col_w_vel>=0)){
                if(partclassj->vis_type==PART_LINES){
                  glutAddMenuEntry(_("    *Lines"),-10-5*j-PART_LINES);
                }
                else{
                  glutAddMenuEntry(_("    Lines"),-10-5*j-PART_LINES);
                }
              }
            }
            if(partclassj->smv_device!=NULL&&partclassj->device_name!=NULL||
              (partclassj->prop!=NULL&&partclassj->prop->smokeview_id!=NULL)
              ){
              if(partclassj->device_name!=NULL){
                strcpy(menulabel,"    ");
                if(partclassj->vis_type==PART_SMV_DEVICE){
                  strcat(menulabel,"*");
                }
                strcat(menulabel,partclassj->device_name);
                glutAddMenuEntry(menulabel,-10-5*j-PART_SMV_DEVICE);
              }
              else if(partclassj->prop!=NULL&&partclassj->prop->smokeview_id!=NULL){
                int iii;
                propdata *propclass;

              // value = iobject*npropinfo + iprop
                propclass=partclassj->prop;
                for(iii=0;iii<propclass->nsmokeview_ids;iii++){
                  int propvalue, showvalue, menuvalue;

                  propvalue = iii*npropinfo + (propclass-propinfo);
                  showvalue = -10-5*j-PART_SMV_DEVICE;
                  menuvalue = (-1-propvalue)*10000 + showvalue;
                  // propvalue = (-menuvalue)/10000-1;
                  // showvalue = -((-menuvalue)%10000)
                  strcpy(menulabel,"    ");
                  if(partclassj->vis_type==PART_SMV_DEVICE&&propclass->smokeview_ids[iii]==propclass->smokeview_id){
                    strcat(menulabel,"*");
                  }
                  strcat(menulabel,propclass->smokeview_ids[iii]);
                  glutAddMenuEntry(menulabel,menuvalue);
                }
              }
            }
          }
          else{
            glutAddMenuEntry(menulabel,-10-5*j);
          }
        }
      }

    CREATEMENU(particlepropshowmenu,ParticlePropShowMenu);
    if(npart5prop>=0){
      glutAddMenuEntry(_("Color with:"),MENU_PROP_DUMMY);
      for(i=0;i<npart5prop;i++){
        partpropdata *propi;
        char menulabel[1024];

        propi = part5propinfo + i;
        if(propi->particle_property==0)continue;
        if(propi->display==1){
          strcpy(menulabel,"  *");
        }
        else{
          strcpy(menulabel,"  ");
        }
        strcat(menulabel,propi->label->longlabel);
        glutAddMenuEntry(menulabel,i);
      }

      if(part5show==0)glutAddMenuEntry(_("  *Hide"), MENU_PROP_HIDEPART);
      if(part5show==1)glutAddMenuEntry(_("  Hide"), MENU_PROP_HIDEPART);
      glutAddMenuEntry("-",MENU_PROP_DUMMY);

      glutAddMenuEntry(_("Draw"),MENU_PROP_DUMMY);
      ntypes=0;
      for(i=0;i<npart5prop;i++){
        partpropdata *propi;
        int j;

        propi = part5propinfo + i;
        if(propi->display==0)continue;
        for(j=0;j<npartclassinfo;j++){
          partclassdata *partclassj;
          char menulabel[1024];

          if(propi->class_present[j]==0)continue;
          partclassj = partclassinfo + j;
          if(partclassj->kind==HUMANS)continue;
          ntypes++;
          if(propi->class_vis[j]==1){
            strcpy(menulabel,"  *");
          }
          else{
            strcpy(menulabel,"  ");
          }
          strcat(menulabel,partclassj->name);
          glutAddSubMenu(menulabel,particlepropshowsubmenu[ntypes-1]);
        }
      }

      if(ntypes>1){
        glutAddMenuEntry(_("  Show all"),MENU_PROP_SHOWALL);
        glutAddMenuEntry(_("  Hide all"),MENU_PROP_HIDEALL);
      }
      glutAddMenuEntry("-",MENU_PROP_DUMMY);
      if(streak5show==1){
        glutAddSubMenu(_("*Streaks"),particlestreakshowmenu);
      }
      else{
        glutAddSubMenu(_("Streaks"),particlestreakshowmenu);
      }
      glutAddMenuEntry("-",MENU_PROP_DUMMY);
      if(show_tracers_always==0)glutAddMenuEntry(_("Show tracers always"),MENU_PROP_TRACERS);
      if(show_tracers_always==1)glutAddMenuEntry(_("*Show tracers always"), MENU_PROP_TRACERS);
    }

    CREATEMENU(humanpropshowmenu,ParticlePropShowMenu);
    if(npart5prop>=0){
      glutAddMenuEntry(_("Color with:"),MENU_PROP_DUMMY);
      for(i=0;i<npart5prop;i++){
        partpropdata *propi;
        char menulabel[1024];

        propi = part5propinfo + i;
        if(propi->human_property==0)continue;
        if(propi->display==1){
          strcpy(menulabel,"  *");
        }
        else{
          strcpy(menulabel,"  ");
        }
        strcat(menulabel,propi->label->longlabel);
        glutAddMenuEntry(menulabel,i);
      }

      if(part5show==0)glutAddMenuEntry(_("  *Hide"),MENU_PROP_HIDEAVATAR);
      if(part5show==1)glutAddMenuEntry(_("  Hide"), MENU_PROP_HIDEAVATAR);
      glutAddMenuEntry("-",MENU_PROP_DUMMY);
      glutAddMenuEntry(_("Draw"),MENU_PROP_DUMMY);
      ntypes=0;
      for(i=0;i<npart5prop;i++){
        partpropdata *propi;
        int j;

        propi = part5propinfo + i;
        if(propi->display==0)continue;
        for(j=0;j<npartclassinfo;j++){
          partclassdata *partclassj;
          char menulabel[1024];

          if(propi->class_present[j]==0)continue;
          partclassj = partclassinfo + j;
          if(partclassj->kind!=HUMANS)continue;
          ntypes++;
          if(propi->class_vis[j]==1){
            strcpy(menulabel,"  *");
          }
          else{
            strcpy(menulabel,"  ");
          }
          strcat(menulabel,partclassj->name);
          glutAddSubMenu(menulabel,particlepropshowsubmenu[ntypes-1]);
        }
        //break;
      }
      if(ntypes>1){
        glutAddMenuEntry(_("  Show all"),MENU_PROP_SHOWALL);
        glutAddMenuEntry(_("  Hide all"),MENU_PROP_HIDEALL);
      }
      glutAddMenuEntry("-",MENU_PROP_DUMMY);
      if(streak5show==1){
        glutAddSubMenu(_("  *Streaks"),particlestreakshowmenu);
      }
      else{
        glutAddSubMenu(_("  Streaks"),particlestreakshowmenu);
      }
    }
  }

/* --------------------------------particle show menu -------------------------- */

  if(npartinfo>0&&nevac!=npartinfo){
    int ii;
    int showall;

    CREATEMENU(particleshowmenu,ParticleShowMenu);
    for(ii=0;ii<npartinfo;ii++){
      partdata *parti;
      char menulabel[1024];

      i = partorderindex[ii];
      parti = partinfo + i;
      if(parti->loaded==0)continue;
      if(parti->evac==1)continue;
      STRCPY(menulabel,"");
      if(parti->display==1)STRCAT(menulabel,"*");
      STRCAT(menulabel,parti->menulabel);
      glutAddMenuEntry(menulabel,-1-i);
    }

    glutAddMenuEntry("-",MENU_DUMMY);
    if(plotstate==DYNAMIC_PLOTS&&visSmokePart!=0){
      if(visSmokePart==2)glutAddMenuEntry(_("*Particles"),MENU_PARTSHOW_PARTICLES);
      if(visSmokePart==1)glutAddMenuEntry(_("#Particles"), MENU_PARTSHOW_PARTICLES);
    }
    else{
      glutAddMenuEntry(_("Particles"), MENU_PARTSHOW_PARTICLES);
    }
    if(havesprinkpart==1){
      if(plotstate==DYNAMIC_PLOTS&&visSprinkPart==1){
        glutAddMenuEntry(_("*Droplets"), MENU_PARTSHOW_DROPLETS);
      }
      else{
        glutAddMenuEntry(_("Droplets"), MENU_PARTSHOW_DROPLETS);
      }
    }
    showall=0;
    if(plotstate==DYNAMIC_PLOTS){
      if(visSprinkPart==1&&visSmokePart!=0)showall=1;
    }
    glutAddMenuEntry("-",MENU_DUMMY);
    if(showall==1){
      glutAddMenuEntry(_("*Show all"), MENU_PARTSHOW_SHOWALL);
    }
    else{
      glutAddMenuEntry(_("Show all"), MENU_PARTSHOW_SHOWALL);
    }
    if(plotstate==DYNAMIC_PLOTS){
      int hideall;

      hideall=1;
      if(visSmokePart!=0)hideall=0;
      if(havesprinkpart==1&&visSprinkPart==1)hideall=0;
      if(hideall==1){
        glutAddMenuEntry(_("*Hide all"), MENU_PARTSHOW_HIDEALL);
      }
      else{
        glutAddMenuEntry(_("Hide all"), MENU_PARTSHOW_HIDEALL);
      }
    }
  }

/* --------------------------------Evac show menu -------------------------- */

  if(nevac>0){
    int ii;

    CREATEMENU(evacshowmenu,EvacShowMenu);
    for(ii=0;ii<npartinfo;ii++){
      partdata *parti;
      char menulabel[1024];

      i = partorderindex[ii];
      parti = partinfo + i;
      if(parti->loaded==0)continue;
      if(parti->evac==0)continue;
      STRCPY(menulabel,"");
      if(parti->display==1)STRCAT(menulabel,"*");
      STRCAT(menulabel,parti->menulabel);
      glutAddMenuEntry(menulabel,-1-i);
    }
    glutAddMenuEntry("-",MENU_DUMMY);
    glutAddMenuEntry(_("Show all"), MENU_PARTSHOW_SHOWALL);
    if(plotstate==DYNAMIC_PLOTS){
      glutAddMenuEntry(_("Hide all"), MENU_PARTSHOW_HIDEALL);
    }
  }

/* -------------------------------- colorbarmenu -------------------------- */

  if(nsmoke3dinfo>0&&Read3DSmoke3DFile==1||nvolrenderinfo>0){
    colorbardata *cbi;
    char ccolorbarmenu[256];

    CREATEMENU(smokecolorbarmenu,SmokeColorbarMenu);

    glutAddMenuEntry(_("Smoke map:"),MENU_DUMMY);
    for(i=0;i<ncolorbars;i++){
      cbi = colorbarinfo + i;

      strcpy(ccolorbarmenu,"  ");
      if(fire_colorbar_index==i){
        strcat(ccolorbarmenu,"*");
        strcat(ccolorbarmenu,cbi->label);
      }
      else{
        strcat(ccolorbarmenu,cbi->label);
      }
      glutAddMenuEntry(ccolorbarmenu,i);
    }
  }


  /* --------------------------------smoke3d showmenu -------------------------- */
  if(nsmoke3dinfo>0&&Read3DSmoke3DFile==1){
    {
      if(nsmoke3dloaded>0){
        CREATEMENU(smoke3dshowmenu,Smoke3DShowMenu);
        for(i=0;i<nsmoke3dinfo;i++){
          smoke3ddata *smoke3di;
          char menulabel[1024];

          smoke3di = smoke3dinfo + i;
          if(smoke3di->loaded==0)continue;
          strcpy(menulabel,"");
          if(smoke3di->display==1)strcpy(menulabel,"*");
          strcat(menulabel,smoke3di->menulabel);
          glutAddMenuEntry(menulabel,i);
        }
        glutAddSubMenu(_("Smoke color map"),smokecolorbarmenu);
        if(have_lighting==1){
          if(show_smoke_lighting==1)glutAddMenuEntry(_("*Light smoke"),HAVE_LIGHT);
          if(show_smoke_lighting==0)glutAddMenuEntry(_("Light smoke"),HAVE_LIGHT);
        }
        if(nsmoke3dinfo>1){
          glutAddMenuEntry("-",MENU_DUMMY);
          glutAddMenuEntry(_("Show all"),SHOW_ALL);
          glutAddMenuEntry(_("Hide all"),HIDE_ALL);
        }
      }
    }
  }

/* --------------------------------iso level menu -------------------------- */

  if(loaded_isomesh!=NULL&&nisoinfo>0&&ReadIsoFile==1){
    CREATEMENU(isolevelmenu,IsoShowMenu);
    if(loaded_isomesh->nisolevels>0&&loaded_isomesh->showlevels!=NULL){
      int showflag,hideflag;
      showflag=1;
      hideflag=1;
      for(i=0;i<loaded_isomesh->nisolevels;i++){
        char levellabel[1024];

        if(loaded_isomesh->showlevels[i]==1){
          sprintf(levellabel,"*%f ",loaded_isomesh->isolevels[i]);
          hideflag=0;
        }
        else{
          showflag=0;
          sprintf(levellabel,"%f ",loaded_isomesh->isolevels[i]);
        }
        if(loaded_isomesh->isofilenum!=-1){
          STRCAT(levellabel,isoinfo[loaded_isomesh->isofilenum].surface_label.unit);
        }
        else{
          STRCAT(levellabel,"");
        }
        glutAddMenuEntry(levellabel,100+i);

      }
      if(showflag == 1)glutAddMenuEntry(_("*Show all levels"), MENU_ISOSHOW_SHOWALL);
      if(showflag == 0)glutAddMenuEntry(_("Show all levels"), MENU_ISOSHOW_SHOWALL);
      if(hideflag == 1)glutAddMenuEntry(_("*Hide all levels"), MENU_ISOSHOW_HIDEALL);
      if(hideflag == 0)glutAddMenuEntry(_("Hide all levels"), MENU_ISOSHOW_HIDEALL);
      glutAddMenuEntry("---",93);
      if(transparent_state == ALL_SOLID)glutAddMenuEntry(_("*All levels solid"), MENU_ISOSHOW_ALLSOLID);
      if(transparent_state != ALL_SOLID)glutAddMenuEntry(_("All levels solid"), MENU_ISOSHOW_ALLSOLID);
      if(transparent_state == ALL_TRANSPARENT)glutAddMenuEntry(_("*All levels transparent"), MENU_ISOSHOW_ALLTRANSPARENT);
      if(transparent_state != ALL_TRANSPARENT)glutAddMenuEntry(_("All levels transparent"), MENU_ISOSHOW_ALLTRANSPARENT);
      if(transparent_state == MIN_SOLID)glutAddMenuEntry(_("*Minimum level solid"), MENU_ISOSHOW_MINSOLID);
      if(transparent_state != MIN_SOLID)glutAddMenuEntry(_("Minimum level solid"), MENU_ISOSHOW_MINSOLID);
      if(transparent_state == MAX_SOLID)glutAddMenuEntry(_("*Maximum level solid"), MENU_ISOSHOW_MAXSOLID);
      if(transparent_state != MAX_SOLID)glutAddMenuEntry(_("Maximum level solid"), MENU_ISOSHOW_MAXSOLID);
    }
    else{
      glutAddMenuEntry(_("Show"), MENU_ISOSHOW_SHOWALL);
      if(visAIso == 0)glutAddMenuEntry(_("*Hide"), MENU_ISOSHOW_HIDEALL);
      if(visAIso != 0)glutAddMenuEntry(_("Hide"), MENU_ISOSHOW_HIDEALL);
    }

/* --------------------------------iso show menu -------------------------- */

    if(nisoinfo>0&&ReadIsoFile==1){
      meshdata *hmesh;
      isodata *iso2;
      int ii;

      CREATEMENU(isoshowsubmenu,IsoShowMenu);
      iso2=NULL;
      for(ii=0;ii<nisoinfo;ii++){
        isodata *isoi;
        char menulabel[1024];

        i = isoorderindex[ii];
        isoi = isoinfo + i;
        if(isoi->loaded==0)continue;
        if(iso2==NULL&&isoi->type==iisotype)iso2=isoi;
        if(plotstate==DYNAMIC_PLOTS&&isoi->display==1&&isoi->type==iisotype){
          iso2=isoi;
          STRCPY(menulabel,"*");
          STRCAT(menulabel,isoi->menulabel);
        }
        else{
          STRCPY(menulabel,isoi->menulabel);
        }
        glutAddMenuEntry(menulabel,1000+i);
      }
      CREATEMENU(isoshowmenu, IsoShowMenu);
      if(iso2!=NULL){
        char menulabel[1024];

        glutAddSubMenu(iso2->surface_label.longlabel, isoshowsubmenu);
        STRCPY(menulabel, _("Show all"));
        STRCAT(menulabel," ");
        STRCAT(menulabel,iso2->surface_label.longlabel);
        STRCAT(menulabel," ");
        STRCAT(menulabel,_("isosurfaces"));
        glutAddMenuEntry(menulabel,SHOWALL_ISO);
        STRCPY(menulabel,_("Hide all"));
        STRCAT(menulabel," ");
        STRCAT(menulabel,iso2->surface_label.longlabel);
        STRCAT(menulabel," ");
        STRCAT(menulabel,_("isosurfaces"));
        glutAddMenuEntry(menulabel,HIDEALL_ISO);
        glutAddMenuEntry("-", MENU_DUMMY);
      }
      if((visAIso & 1) == 1)glutAddMenuEntry(_("*Solid"), MENU_ISOSHOW_SOLID);
      if((visAIso & 1) != 1)glutAddMenuEntry(_("Solid"), MENU_ISOSHOW_SOLID);
      if((visAIso & 2) == 2)glutAddMenuEntry(_("*Outline"), MENU_ISOSHOW_OUTLINE);
      if((visAIso & 2) != 2)glutAddMenuEntry(_("Outline"), MENU_ISOSHOW_OUTLINE);
      if((visAIso & 4) == 4)glutAddMenuEntry(_("*Points"), MENU_ISOSHOW_POINTS);
      if((visAIso & 4) != 4)glutAddMenuEntry(_("Points"), MENU_ISOSHOW_POINTS);
      hmesh=meshinfo+highlight_mesh;
      if(hmesh->isofilenum!=-1){
        char levellabel[1024];

        STRCPY(levellabel,isoinfo[hmesh->isofilenum].surface_label.shortlabel);
        STRCAT(levellabel," ");
        STRCAT(levellabel,_("Levels"));
        glutAddSubMenu(levellabel,isolevelmenu);
      }
      if(niso_compressed==0){
        if(smooth_iso_normal == 1)glutAddMenuEntry(_("*Smooth"), MENU_ISOSHOW_SMOOTH);
        if(smooth_iso_normal == 0)glutAddMenuEntry(_("Smooth"), MENU_ISOSHOW_SMOOTH);
      }
      if(show_iso_normal == 1)glutAddMenuEntry(_("*Show normals"), MENU_ISOSHOW_NORMALS);
      if(show_iso_normal == 0)glutAddMenuEntry(_("Show normals"), MENU_ISOSHOW_NORMALS);
    }
  }

/* -------------------------------- colorbarmenu -------------------------- */

  CREATEMENU(colorbarshademenu,ColorbarMenu);
  if(contour_type==SHADED_CONTOURS){
    glutAddMenuEntry("*Continuous",COLORBAR_CONTINUOUS);
    glutAddMenuEntry("Stepped",COLORBAR_STEPPED);
    glutAddMenuEntry("Lines",COLORBAR_LINES);
  }
  else if(contour_type==STEPPED_CONTOURS){
    glutAddMenuEntry("Continuous",COLORBAR_CONTINUOUS);
    glutAddMenuEntry("*Stepped",COLORBAR_STEPPED);
    glutAddMenuEntry("Lines",COLORBAR_LINES);
  }else if(contour_type==LINE_CONTOURS){
    glutAddMenuEntry("Continuous",COLORBAR_CONTINUOUS);
    glutAddMenuEntry("Stepped",COLORBAR_STEPPED);
    glutAddMenuEntry("*Lines",COLORBAR_LINES);
  }
  glutAddMenuEntry("-",MENU_DUMMY);
  if(show_extreme_maxdata == 1){
    glutAddMenuEntry(_("  *Highlight data above specified max"), COLORBAR_HIGHLIGHT_ABOVE);
  }
  else{
    glutAddMenuEntry(_("  Highlight data above specified max"), COLORBAR_HIGHLIGHT_ABOVE);
  }
  if(show_extreme_mindata == 1){
    glutAddMenuEntry(_("  *Highlight data below specified min"), COLORBAR_HIGHLIGHT_BELOW);
  }
  else{
    glutAddMenuEntry(_("  Highlight data below specified min"), COLORBAR_HIGHLIGHT_BELOW);
  }
  if(colorbarflip == 1){
    glutAddMenuEntry(_("  *Flip"), COLORBAR_FLIP);
  }
  else{
    glutAddMenuEntry(_("  Flip"), COLORBAR_FLIP);
  }

  CREATEMENU(colorbarsmenu,ColorbarMenu);
  {
    colorbardata *cbi;
    char ccolorbarmenu[256];

    for(i=0;i<ncolorbars;i++){
      cbi = colorbarinfo + i;

      strcpy(ccolorbarmenu,"  ");
      if(colorbartype==i){
        strcat(ccolorbarmenu,"*");
        strcat(ccolorbarmenu,cbi->label);
      }
      else{
        strcat(ccolorbarmenu,cbi->label);
      }
      glutAddMenuEntry(ccolorbarmenu,i);
    }
  }

/* -------------------------------- colorbarmenu -------------------------- */

  CREATEMENU(colorbarmenu,ColorbarMenu);
  glutAddSubMenu(_("Colorbar"),colorbarsmenu);
  glutAddSubMenu(_("Colorbar type"), colorbarshademenu);
  if(use_transparency_data==1){
    glutAddMenuEntry(_("  *Transparent (data)"),COLORBAR_TRANSPARENT);
  }
  else{
    glutAddMenuEntry(_("  Transparent (data)"),COLORBAR_TRANSPARENT);
  }
  if(setbwdata == 1)glutAddMenuEntry("*Black/White (data)", COLORBAR_TOGGLE_BW_DATA);
  if(setbwdata == 0)glutAddMenuEntry("Black/White  (data)", COLORBAR_TOGGLE_BW_DATA);
  if(setbw == 1)glutAddMenuEntry("*Black/White (geometry)", COLORBAR_TOGGLE_BW);
  if(setbw == 0)glutAddMenuEntry("Black/White (geometry)", COLORBAR_TOGGLE_BW);
  glutAddMenuEntry(_("  Reset"), COLORBAR_RESET);

/* --------------------------------showVslice menu -------------------------- */
  if(nvsliceloaded==0){
    vd_shown=NULL;
  }
  if(nvsliceinfo>0&&nvsliceloaded>0){
    CREATEMENU(showvslicemenu,ShowVSliceMenu);
    nvsliceloaded0=0;
    for(i=0;i<nvsliceinfo;i++){
      vslicedata *vd;
      slicedata *sd;
      char menulabel[1024];

      vd = vsliceinfo + i;
      if(vd->loaded==0)continue;
      sd = sliceinfo + vd->ival;
      nvsliceloaded0++;
      STRCPY(menulabel,"");
      if(plotstate==DYNAMIC_PLOTS&&sd->type==islicetype&&vd->display==1){
        vd_shown=vd;
        STRCAT(menulabel,"*");
      }
      STRCAT(menulabel,sd->menulabel2);
      if(sd->slicelabel!=NULL){
        STRCAT(menulabel," - ");
        STRCAT(menulabel,sd->slicelabel);
      }
      glutAddMenuEntry(menulabel,i);
    }
    if(show_slice_in_obst == 1)glutAddMenuEntry(_("*Show vector slice in blockage"), MENU_SHOWSLICE_INBLOCKAGE);
    if(show_slice_in_obst == 0)glutAddMenuEntry(_("Show vector slice in blockage"), MENU_SHOWSLICE_INBLOCKAGE);
    if(show_slices_and_vectors == 1)glutAddMenuEntry(_("*Show slices and vectors"), MENU_SHOWSLICE_SLICEANDVECTORS);
    if(show_slices_and_vectors == 0)glutAddMenuEntry(_("Show slices and vectors"), MENU_SHOWSLICE_SLICEANDVECTORS);
    if(offset_slice == 1)glutAddMenuEntry(_("*Offset vector slice"), MENU_SHOWSLICE_OFFSET);
    if(offset_slice == 0)glutAddMenuEntry(_("Offset vector slice"), MENU_SHOWSLICE_OFFSET);
    if(vd_shown!=NULL&&nvsliceloaded0!=0){
      char menulabel[1024];

      glutAddMenuEntry("",MENU_DUMMY);
      STRCPY(menulabel,_("Show all"));
      STRCAT(menulabel," ");
      STRCAT(menulabel,sliceinfo[vd_shown->ival].label.longlabel);
      STRCAT(menulabel," ");
      STRCAT(menulabel,_("vector slices"));
      glutAddMenuEntry(menulabel,SHOW_ALL);
      STRCPY(menulabel,_("Hide all"));
      STRCAT(menulabel," ");
      STRCAT(menulabel,sliceinfo[vd_shown->ival].label.longlabel);
      STRCAT(menulabel," ");
      STRCAT(menulabel,_("vector slices"));
      glutAddMenuEntry(menulabel,HIDE_ALL);
    }
  }
  if(nsliceinfo>0&&nmultisliceinfo<nsliceinfo){
    CREATEMENU(showmultislicemenu,ShowMultiSliceMenu);
    for(i=0;i<nmultisliceinfo;i++){
      slicedata *sd;
      char menulabel[1024];
      multislicedata *mslicei;

      mslicei = multisliceinfo + i;
      if(mslicei->loaded==0)continue;
      sd = sliceinfo + mslicei->islices[0];
      STRCPY(menulabel,"");
      if(plotstate==DYNAMIC_PLOTS&&mslicei->display!=0&&mslicei->type==islicetype){
        if(mslicei->display==1){
          STRCAT(menulabel,"*");
        }
        else if(mslicei->display==-1){
          STRCAT(menulabel,"#");
        }
      }
      STRCAT(menulabel,mslicei->menulabel2);
      if(sd->slicelabel!=NULL){
        STRCAT(menulabel," - ");
        STRCAT(menulabel,sd->slicelabel);
      }
      glutAddMenuEntry(menulabel,i);
    }
    if(show_slice_in_obst == 1)glutAddMenuEntry(_("*Show multi slice in blockage"), MENU_SHOWSLICE_INBLOCKAGE);
    if(show_slice_in_obst == 0)glutAddMenuEntry(_("Show multi slice in blockage"), MENU_SHOWSLICE_INBLOCKAGE);
    if(offset_slice == 1)glutAddMenuEntry(_("*Offset slice"), MENU_SHOWSLICE_OFFSET);
    if(offset_slice == 0)glutAddMenuEntry(_("Offset slice"), MENU_SHOWSLICE_OFFSET);
    if(nfedinfo>0){
      if(show_fed_area == 1)glutAddMenuEntry("*Show FED areas", MENU_SHOWSLICE_FEDAREA);
      if(show_fed_area == 0)glutAddMenuEntry("Show FED areas", MENU_SHOWSLICE_FEDAREA);
    }
  }

/* --------------------------------showslice menu -------------------------- */
  if(nsliceloaded==0){
    sd_shown=NULL;
  }
  if(nsliceinfo>0&&nsliceloaded>0){
    int ii;

    CREATEMENU(showhideslicemenu,ShowHideSliceMenu);
    for(ii=0;ii<nslice_loaded;ii++){
      slicedata *sd;
      char menulabel[1024];

      i = slice_loaded_list[ii];
      sd = sliceinfo + i;
      if(sd_shown==NULL&&sd->type==islicetype){
        sd_shown = sd;
      }
      STRCPY(menulabel,"");
      if(plotstate==DYNAMIC_PLOTS&&sd->display==1&&sd->type==islicetype){
        sd_shown=sd;
        STRCAT(menulabel,"*");
      }
      STRCAT(menulabel,sd->menulabel2);
      if(sd->slicelabel!=NULL){
        STRCAT(menulabel," - ");
        STRCAT(menulabel,sd->slicelabel);
      }
      glutAddMenuEntry(menulabel,i);
    }
    glutAddMenuEntry("-",MENU_DUMMY);
    if(have_terrain_slice() == 1){
      if(planar_terrain_slice == 1)glutAddMenuEntry(_("*Planar terrain slice"), MENU_SHOWSLICE_TERRAIN);
      if(planar_terrain_slice == 0)glutAddMenuEntry(_("Planar terrain slice"), MENU_SHOWSLICE_TERRAIN);
    }
    if(show_slice_in_obst == 1)glutAddMenuEntry(_("*Show slice in blockage"), MENU_SHOWSLICE_INBLOCKAGE);
    if(show_slice_in_obst == 0)glutAddMenuEntry(_("Show slice in blockage"), MENU_SHOWSLICE_INBLOCKAGE);
    if(offset_slice == 1)glutAddMenuEntry(_("*Offset slice"), MENU_SHOWSLICE_OFFSET);
    if(offset_slice == 0)glutAddMenuEntry(_("Offset slice"), MENU_SHOWSLICE_OFFSET);
    if(show_slices_and_vectors == 1)glutAddMenuEntry(_("*Show slices and vectors"), MENU_SHOWSLICE_SLICEANDVECTORS);
    if(show_slices_and_vectors == 0)glutAddMenuEntry(_("Show slices and vectors"), MENU_SHOWSLICE_SLICEANDVECTORS);
    if(nfedinfo>0){
      if(show_fed_area == 1)glutAddMenuEntry("*Show FED areas", MENU_SHOWSLICE_FEDAREA);
      if(show_fed_area == 0)glutAddMenuEntry("Show FED areas", MENU_SHOWSLICE_FEDAREA);
    }
    if(nsliceloaded>0&&sd_shown!=NULL){
      char menulabel[1024];

      glutAddMenuEntry("-",MENU_DUMMY);
      STRCPY(menulabel,_("Show all"));
      STRCAT(menulabel," ");
      STRCAT(menulabel,sd_shown->label.longlabel);
      STRCAT(menulabel," ");
      STRCAT(menulabel,_("slices"));
      glutAddMenuEntry(menulabel,SHOW_ALL);
      STRCPY(menulabel,_("Hide all"));
      STRCAT(menulabel," ");
      STRCAT(menulabel,sd_shown->label.longlabel);
      STRCAT(menulabel," ");
      STRCAT(menulabel,_("slices"));
      glutAddMenuEntry(menulabel,HIDE_ALL);
    }
  }

/* -------------------------------- avatartour menu -------------------------- */

  CREATEMENU(avatarevacmenu,AvatarEvacMenu);
  if(navatar_types>0){
    if(iavatar_evac==-1){
      glutAddMenuEntry(_("*Defined in evac file"),MENU_AVATAR_DEFINED);
    }
    else{
      glutAddMenuEntry(_("Defined in evac file"),MENU_AVATAR_DEFINED);
    }
    glutAddMenuEntry("-",MENU_DUMMY);
    for(i=0;i<navatar_types;i++){
      char menulabel[1024];

      strcpy(menulabel,"");
      if(iavatar_evac==i){
        strcat(menulabel,"*");
      }
      strcat(menulabel,avatar_types[i]->label);
      glutAddMenuEntry(menulabel,i);
    }
  }
  CREATEMENU(avatartourmenu,TourMenu);
  if(navatar_types>0){
    if(selectedtour_index>=0&&selectedtour_index<ntours){
      tourdata *touri;
      char menulabel[1024];

      touri = tourinfo + selectedtour_index;
      strcpy(menulabel,"For ");
      strcat(menulabel,touri->label);
      glutAddMenuEntry(menulabel,MENU_DUMMY);
      glutAddMenuEntry("-",MENU_DUMMY);
    }

    for(i=0;i<navatar_types;i++){
      char menulabel[1024];

      strcpy(menulabel,"");
      if(tourlocus_type==2&&iavatar_types==i){
        strcat(menulabel,"*");
      }
      strcat(menulabel,avatar_types[i]->label);
      glutAddMenuEntry(menulabel,-23-i);
    }
  }

    /* --------------------------------tour menu -------------------------- */

  CREATEMENU(tourmenu,TourMenu);

  glutAddMenuEntry(_("New..."),MENU_TOUR_NEW);
  if(ntours>0){
    if(showtour_dialog==1){
      glutAddMenuEntry(_("*Edit..."),MENU_TOUR_EDIT);
    }
    else{
      glutAddMenuEntry(_("Edit..."),MENU_TOUR_EDIT);
    }
    if(ntours>0)glutAddMenuEntry("-",MENU_DUMMY);
    for(i=0;i<ntours;i++){
      tourdata *touri;
      int glui_avatar_index_local;
      char menulabel[1024];

      touri = tourinfo + i;
      if(touri->display==1){
        STRCPY(menulabel,"");
        if(selectedtour_index==i){
          STRCAT(menulabel,"@");
        }
        STRCAT(menulabel,"*");
        STRCAT(menulabel,touri->menulabel);
      }
      else{
        STRCPY(menulabel,touri->menulabel);
      }
      glui_avatar_index_local = touri->glui_avatar_index;
      if(glui_avatar_index_local>=0&&glui_avatar_index_local<navatar_types){
        sv_object *avatari;

        avatari=avatar_types[glui_avatar_index_local];
        strcat(menulabel,"(");
        strcat(menulabel,avatari->label);
        strcat(menulabel,")");
      }
      glutAddMenuEntry(menulabel,i);
    }
    if(selectedtour_index>=0){
      char menulabel[1024];

      strcpy(menulabel,"");
      if(viewtourfrompath==1)strcat(menulabel,"*");
      strcat(menulabel,"View from ");
      strcat(menulabel,tourinfo[selectedtour_index].label);
      glutAddMenuEntry(menulabel,MENU_TOUR_VIEWFROMROUTE);
    }
    glutAddMenuEntry("-",MENU_DUMMY);
    glutAddMenuEntry(_("Show all"),MENU_TOUR_SHOWALL);
    glutAddMenuEntry(_("Hide all"),MENU_TOUR_MANUAL);
  }

 /* --------------------------------Show Volume smoke menu -------------------------- */

  if(nvolsmoke3dloaded>0){
    CREATEMENU(showvolsmoke3dmenu,ShowVolSmoke3DMenu);
    if(nvolsmoke3dloaded>1){
      char vlabel[256];

      strcpy(vlabel,_("3D smoke (volume rendered)"));
      strcat(vlabel,_(" - Show all"));
      glutAddMenuEntry(vlabel,SHOW_ALL);

      strcpy(vlabel,_("3D smoke (volume rendered)"));
      strcat(vlabel,_(" - Hide all"));
      glutAddMenuEntry(vlabel,HIDE_ALL);
      glutAddMenuEntry("-",MENU_DUMMY);
    }
    for(i=0;i<nmeshes;i++){
      meshdata *meshi;
      volrenderdata *vr;
      char menulabel[1024];

      meshi = meshinfo + i;
      vr = &(meshi->volrenderinfo);
      if(vr->fireslice==NULL||vr->smokeslice==NULL)continue;
      if(vr->loaded==0)continue;
      strcpy(menulabel,"");
      if(vr->display==1)strcat(menulabel,"*");
      strcat(menulabel,meshi->label);
      glutAddMenuEntry(menulabel,i);
    }
    glutAddMenuEntry("-",MENU_DUMMY);
    glutAddSubMenu(_("Smoke color map"),smokecolorbarmenu);
  }

  CREATEMENU(aperturemenu,ApertureMenu);
  if(apertureindex==0)glutAddMenuEntry("*30",0);
  if(apertureindex!=0)glutAddMenuEntry("30",0);
  if(apertureindex==1)glutAddMenuEntry("*45",1);
  if(apertureindex!=1)glutAddMenuEntry("45",1);
  if(apertureindex==2)glutAddMenuEntry("*60",2);
  if(apertureindex!=2)glutAddMenuEntry("60",2);
  if(apertureindex==3)glutAddMenuEntry("*75",3);
  if(apertureindex!=3)glutAddMenuEntry("75",3);
  if(apertureindex==4)glutAddMenuEntry("*90",4);
  if(apertureindex!=4)glutAddMenuEntry("90",4);

  CREATEMENU(zoommenu,ZoomMenu);
  if(zoomindex==0)glutAddMenuEntry("*0.25",0);
  if(zoomindex!=0)glutAddMenuEntry("0.25",0);
  if(zoomindex==1)glutAddMenuEntry("*0.50",1);
  if(zoomindex!=1)glutAddMenuEntry("0.50",1);
  if(zoomindex==2)glutAddMenuEntry("*1.0",2);
  if(zoomindex!=2)glutAddMenuEntry("1.0",2);
  if(zoomindex==3)glutAddMenuEntry("*2.0",3);
  if(zoomindex!=3)glutAddMenuEntry("2.0",3);
  if(zoomindex==4)glutAddMenuEntry("*4.0",4);
  if(zoomindex!=4)glutAddMenuEntry("4.0",4);

  /* --------------------------------reset menu -------------------------- */

  CREATEMENU(resetmenu,ResetMenu);
  {
    char line[256];
    cameradata *ca;

    if(trainer_mode==1){
      if(visBlocks==visBLOCKOutline){
        glutAddMenuEntry(_("*Outline"),MENU_OUTLINEVIEW);
      }
      else{
        glutAddMenuEntry(_("Outline"),MENU_OUTLINEVIEW);
      }
      glutAddMenuEntry("-",MENU_DUMMY);
    }
    if(trainer_mode==0){
      glutAddMenuEntry(_("Save"),SAVE_VIEWPOINT);
      glutAddMenuEntry(_("Set as Startup"),MENU_STARTUPVIEW);
      glutAddSubMenu(_("Zoom"),zoommenu);
      if(projection_type==1)glutAddMenuEntry(_("Switch to perspective view       ALT v"),MENU_SIZEPRESERVING);
      if(projection_type==0)glutAddMenuEntry(_("Switch to size preserving view   ALT v"),MENU_SIZEPRESERVING);
      glutAddMenuEntry("-",MENU_DUMMY);
    }
    for(ca=camera_list_first.next;ca->next!=NULL;ca=ca->next){
      if(trainer_mode==1&&strcmp(ca->name,_("internal"))==0)continue;
      strcpy(line,"");
      if(ca->view_id==selected_view){
        strcat(line,"*");
      }
      if(trainer_mode==1&&strcmp(ca->name,_("external"))==0){
        strcat(line,_("Outside"));
      }
      else{
        strcat(line,ca->name);
      }
      glutAddMenuEntry(line,ca->view_id);
    }
  }
  if(trainer_mode==0&&showtime==1){
    glutAddMenuEntry("-",MENU_DUMMY);
    glutAddMenuEntry(_("Time"),MENU_TIMEVIEW);
  }

  /* --------------------------------showhide menu -------------------------- */

  showhide_data = 0;
  CREATEMENU(showhidemenu,ShowHideMenu);
  glutAddSubMenu(_("Color"), colorbarmenu);
  glutAddSubMenu(_("Geometry"),geometrymenu);
  glutAddSubMenu(_("Labels"),labelmenu);
  glutAddSubMenu(_("Viewpoints"), resetmenu);
  glutAddMenuEntry("-", MENU_DUMMY);
  if(Read3DSmoke3DFile==1&&nsmoke3dloaded>0){
    showhide_data = 1;
    glutAddSubMenu(_("3D smoke"), smoke3dshowmenu);
  }
  if(nvolsmoke3dloaded>0){
    char vlabel[256];

    showhide_data = 1;
    strcpy(vlabel, _("3D smoke (volume rendered)"));
    glutAddSubMenu(vlabel, showvolsmoke3dmenu);
  }
  if(npatchloaded>0){
    showhide_data = 1;
    glutAddSubMenu(_("Boundaries"), showpatchmenu);
  }
  {
    int human_present=0;
    int particle_present=0;

    if(npart5loaded>0){
      int ii;

      showhide_data = 1;
      for(ii = 0; ii<npartinfo; ii++){
        partdata *parti;

        parti = partinfo + ii;
        if(parti->loaded==0)continue;
        if(parti->evac==1)human_present=1;
        if(parti->evac==0)particle_present=1;
      }
      if(particle_present==1){
        glutAddSubMenu(_("Particles"),particlepropshowmenu);
      }
      if(human_present==1){
        glutAddSubMenu(_("Humans"),humanpropshowmenu);
      }
    }
  }

  if(ReadIsoFile==1){
    int niso_loaded=0;

    for(i=0;i<nisoinfo;i++){
      isodata *isoi;

      isoi = isoinfo + i;
      if(isoi->loaded==1)niso_loaded++;
    }

    if(niso_loaded>1){
     glutAddSubMenu(_("Isosurfaces"),isoshowmenu);
    }
    else{
     glutAddSubMenu(_("Isosurfaces"),isoshowmenu);
    }
    showhide_data = 1;
  }
  if(ReadPlot3dFile==1){
    showhide_data = 1;
    glutAddSubMenu(_("Plot3d"), plot3dshowmenu);
  }

  nvslice0=0, nvslice1=0, nvslice2=0;
  nvsliceloaded0=0, nvsliceloaded1=0;
  nvsliceloaded2=0;
  for(i=0;i<nvsliceinfo;i++){
    vslicedata *vd;

    vd = vsliceinfo+i;
    switch(vd->vec_type){
    case 0:
      nvslice0++;
      if(vd->loaded==1)nvsliceloaded0++;
      break;
    case 1:
      nvslice1++;
      if(vd->loaded==1)nvsliceloaded1++;
      break;
    case 2:
      nvslice2++;
      if(vd->loaded==1)nvsliceloaded1++;
      break;
     default:
      ASSERT(FFALSE);
      break;
    }
  }
  if(nvsliceloaded0+nvsliceloaded1+nvsliceloaded2>0){
    showhide_data = 1;
    glutAddSubMenu(_("Slices (vector)"), showvslicemenu);
  }
  if(nsliceloaded>0){
    showhide_data = 1;
    glutAddSubMenu(_("Slices"), showhideslicemenu);
    if(nmultisliceinfo<nsliceinfo){
      glutAddSubMenu(_("Multi-Slices"),showmultislicemenu);
    }
  }

  if(nzoneinfo>0&&(ReadZoneFile==1||nzvents>0)){
    showhide_data = 1;
    glutAddSubMenu(_("Zone"), zoneshowmenu);
  }
  if(nobject_defs>0){
    int num_activedevices=0;

    for(i = 0; i<nobject_defs; i++){
      sv_object *obj_typei;

      obj_typei = object_defs[i];
      if(obj_typei->used==1)num_activedevices++;
    }

    if(num_activedevices>0){
      /*
      if(isZoneFireModel==0||(isZoneFireModel==1&&num_activedevices>1)){
        glutAddSubMenu(_("Objects"),showobjectsmenu);
      }
      */
    }
    else{
      for(i=0;i<nobject_defs;i++){
        sv_object *obj_typei;

        obj_typei = object_defs[i];
        if(obj_typei->used==1){
          char obj_menu[256];

          strcpy(obj_menu,"");
          if(obj_typei->visible==1){
            strcat(obj_menu,"*");
          }
          strcat(obj_menu,obj_typei->label);
          showhide_data = 1;
          glutAddMenuEntry(obj_menu, i);
        }
      }
    }
  }
  if(ntc_total>0){
    showhide_data = 1;
    if(isZoneFireModel==1){
      if(visSensor==1)glutAddMenuEntry(_("*Targets"), MENU_SHOWHIDE_SENSOR);
      if(visSensor==0)glutAddMenuEntry(_("Targets"), MENU_SHOWHIDE_SENSOR);
      if(hasSensorNorm==1){
        if(visSensorNorm==1)glutAddMenuEntry(_("*Target orientation"), MENU_SHOWHIDE_SENSOR_NORM);
        if(visSensorNorm==0)glutAddMenuEntry(_("Target orientation"), MENU_SHOWHIDE_SENSOR_NORM);
      }
    }
    else{
      if(visSensor==1)glutAddMenuEntry(_("*Thermocouples"), MENU_SHOWHIDE_SENSOR);
      if(visSensor==0)glutAddMenuEntry(_("Thermocouples"), MENU_SHOWHIDE_SENSOR);
      if(hasSensorNorm==1){
        if(visSensorNorm==1)glutAddMenuEntry(_("*Thermocouple norms"), MENU_SHOWHIDE_SENSOR_NORM);
        if(visSensorNorm==0)glutAddMenuEntry(_("Thermocouple norms"), MENU_SHOWHIDE_SENSOR_NORM);
      }
    }
  }
  if(showhide_data==1)glutAddMenuEntry("-", MENU_DUMMY);
  if(background_flip==1){
    glutAddMenuEntry(_("*Flip background"), MENU_SHOWHIDE_FLIP);
  }
  else{
    glutAddMenuEntry(_("Flip background"), MENU_SHOWHIDE_FLIP);
  }
  if(titlesafe_offset==0)glutAddMenuEntry(_("Offset window"), MENU_SHOWHIDE_OFFSET);
  if(titlesafe_offset!=0)glutAddMenuEntry(_("*Offset window"),MENU_SHOWHIDE_OFFSET);
  if(ntextures_loaded_used>0){
    glutAddSubMenu(_("Textures"),textureshowmenu);
  }
#ifdef pp_MEMPRINT
  glutAddMenuEntry("Show Memory block info",MENU_SHOWHIDE_PRINT);
#endif


/* --------------------------------frame rate menu -------------------------- */

  CREATEMENU(frameratemenu,FrameRateMenu);
  if(frameratevalue==1)glutAddMenuEntry("*1 FPS",1);
  if(frameratevalue!=1)glutAddMenuEntry("1 FPS",1);
  if(frameratevalue==2)glutAddMenuEntry("*2 FPS",2);
  if(frameratevalue!=2)glutAddMenuEntry("2 FPS",2);
  if(frameratevalue==4)glutAddMenuEntry("*4 FPS",4);
  if(frameratevalue!=4)glutAddMenuEntry("4 FPS",4);
  if(frameratevalue==8)glutAddMenuEntry("*8 FPS",8);
  if(frameratevalue!=8)glutAddMenuEntry("8 FPS",8);
  if(frameratevalue==10)glutAddMenuEntry("*10 FPS",10);
  if(frameratevalue!=10)glutAddMenuEntry("10 FPS",10);
  if(frameratevalue==15)glutAddMenuEntry("*15 FPS",15);
  if(frameratevalue!=15)glutAddMenuEntry("15 FPS",15);
  if(frameratevalue==30)glutAddMenuEntry("*30 FPS",30);
  if(frameratevalue!=30)glutAddMenuEntry("30 FPS",30);
  if(frameratevalue==2001)glutAddMenuEntry(_("*Real time"),MENU_FRAMERATE_Realtime);
  if(frameratevalue!=2001)glutAddMenuEntry(_("Real time"), MENU_FRAMERATE_Realtime);
  if(frameratevalue==2002)glutAddMenuEntry(_("*2 x Real time"), MENU_FRAMERATE_2xRealtime);
  if(frameratevalue!=2002)glutAddMenuEntry(_("2 x Real time"), MENU_FRAMERATE_2xRealtime);
  if(frameratevalue==2004)glutAddMenuEntry(_("*4 x Real time"), MENU_FRAMERATE_4xRealtime);
  if(frameratevalue!=2004)glutAddMenuEntry(_("4 x Real time"), MENU_FRAMERATE_4xRealtime);
  if(frameratevalue!=1000)glutAddMenuEntry(_("Unlimited"),1000);
  if(frameratevalue==1000)glutAddMenuEntry(_("*Unlimited"),1000);
  if(frameratevalue<0){
    glutAddMenuEntry(_("*Step"),-1);
  }
  else{
    glutAddMenuEntry(_("Step"),-1);
  }

/* --------------------------------render menu -------------------------- */
  {
    char renderwindow[1024];
    char renderwindow2[1024];
    char renderwindow3[1024];
    char rendertemp[1024];
    int render_current=0;

    strcpy(renderwindow,"  ");
    if(renderW==320)strcat(renderwindow,"*");
    strcat(renderwindow,"320x240");

    strcpy(renderwindow2,"  ");
    if(renderW==640)strcat(renderwindow2,"*");
    strcat(renderwindow2,"640x480");

    sprintf(rendertemp,"%i%s%i (current)",screenWidth,"x",screenHeight);
    strcpy(renderwindow3,"  ");
    if(renderW!=320&&renderW!=640&&renderW!=2*screenWidth){
      render_current=1;
      strcat(renderwindow3,"*");
    }
    strcat(renderwindow3,rendertemp);

    CREATEMENU(startrenderingmenu,RenderMenu);
    glutAddMenuEntry(_("  One frame (single part)"),RENDER_CURRENT_SINGLE);
    if(render_current==1){
      char menulabel[1024];

      sprintf(menulabel,"  One frame (%i x %i parts)",nrender_rows,nrender_rows);
      glutAddMenuEntry(menulabel,RENDER_CURRENT_MULTIPLE);
    }
    if(RenderTime==1||touring==1){
      glutAddMenuEntry(_("  All frames"),1);
      glutAddMenuEntry(_("  Every 2nd frame"),2);
      glutAddMenuEntry(_("  Every 3rd frame"),3);
      glutAddMenuEntry(_("  Every 4th frame"),4);
      glutAddMenuEntry(_("  Every 5th frame"),5);
      glutAddMenuEntry(_("  Every 10th frame"),10);
      glutAddMenuEntry(_("  Every 20th frame"),20);
      glutAddMenuEntry(_("  Cancel"),RenderCancel);
    }

    CREATEMENU(resolutionmultipliermenu,RenderMenu);
    if(nrender_rows==2){
      glutAddMenuEntry("  *2",HIDEALL_ISO);
    }
    else{
      glutAddMenuEntry("  2",HIDEALL_ISO);
    }
    if(nrender_rows==3){
      glutAddMenuEntry("  *3",10003);
    }
    else{
      glutAddMenuEntry("  3",10003);
    }
    if(nrender_rows==4){
      glutAddMenuEntry("  *4",10004);
    }
    else{
      glutAddMenuEntry("  4",10004);
    }
    if(nrender_rows==5){
      glutAddMenuEntry("  *5",10005);
    }
    else{
      glutAddMenuEntry("  5",10005);
    }
    if(nrender_rows==6)glutAddMenuEntry("  *6",10005);
    if(nrender_rows==7)glutAddMenuEntry("  *7",10005);
    if(nrender_rows==8)glutAddMenuEntry("  *8",10005);
    if(nrender_rows==9)glutAddMenuEntry("  *9",10005);
    if(nrender_rows==10)glutAddMenuEntry("  *10",10005);

    CREATEMENU(rendermenu,RenderMenu);
    glutAddMenuEntry(_("Resolution:"),11000);
    glutAddMenuEntry(renderwindow,Render320);
    glutAddMenuEntry(renderwindow2,Render640);
    glutAddMenuEntry(renderwindow3,RenderWindow);

    if(render_current==1)glutAddSubMenu(_("Resolution multiplier"),resolutionmultipliermenu);

    glutAddMenuEntry(_("Type:"),11000);
    if(render_filetype==PNG){
      glutAddMenuEntry("  *PNG",RenderPNG);
      glutAddMenuEntry("  JPEG",RenderJPEG);
    }
    if(render_filetype==JPEG){
      glutAddMenuEntry("  PNG",RenderPNG);
      glutAddMenuEntry("  *JPEG",RenderJPEG);
    }
    if(render_filetype==IMAGE_NONE){
      glutAddMenuEntry("  PNG",RenderPNG);
      glutAddMenuEntry("  JPEG",RenderJPEG);
    }

    glutAddMenuEntry(_("File suffix:"),11000);
    if(render_label_type==RENDER_LABEL_FRAMENUM){
      glutAddMenuEntry(_("  *Frame number"),RenderLABELframenumber);
      glutAddMenuEntry(_("  Time"),RenderLABELtime);
    }
    if(render_label_type==RENDER_LABEL_TIME){
      glutAddMenuEntry(_("  Frame number"),RenderLABELframenumber);
      glutAddMenuEntry(_("  *Time"),RenderLABELtime);
    }

    glutAddSubMenu(_("Start rendering:"),startrenderingmenu);
    UpdateGluiRender();
  }

  /* --------------------------------filesdialog menu -------------------------- */

  CREATEMENU(filesdialogmenu, DialogMenu);
  glutAddMenuEntry(_("Auto load data files..."), DIALOG_AUTOLOAD);
#ifdef pp_COMPRESS
  if(smokezippath!=NULL&&(npatchinfo>0||nsmoke3dinfo>0||nsliceinfo>0)){
    glutAddMenuEntry(_("Compress data files...  ALT z"), DIALOG_SMOKEZIP);
  }
#endif
  glutAddMenuEntry(_("Save/load configuration files..."), DIALOG_CONFIG);
  glutAddMenuEntry(_("Render images..."), DIALOG_RENDER);
  if(have_ffmpeg==1){
    glutAddMenuEntry(_("Make movies..."), DIALOG_MOVIE);
  }
  glutAddMenuEntry(_("Record/run scripts..."), DIALOG_SCRIPT);

  /* --------------------------------viewdialog menu -------------------------- */

  CREATEMENU(viewdialogmenu, DialogMenu);
  glutAddMenuEntry(_("Create/edit tours...  ALT t"), DIALOG_TOUR);
  glutAddMenuEntry(_("Edit colorbar...  ALT C"), DIALOG_COLORBAR);
  if(isZoneFireModel==0){
    glutAddMenuEntry(_("Examine geometry...  ALT e"), DIALOG_GEOMETRY);
  }
  glutAddMenuEntry(_("Stereo settings..."), DIALOG_STEREO);
  if(trainer_active==1){
    glutAddMenuEntry(_("Trainer..."), DIALOG_TRAINER);
  }

  /* --------------------------------datadialog menu -------------------------- */

  CREATEMENU(datadialogmenu, DialogMenu);
  glutAddMenuEntry(_("Coloring..."), DIALOG_COLORING);
  if(ndeviceinfo>0&&GetNumActiveDevices()>0){
    glutAddMenuEntry(_("Devices/Objects..."), DIALOG_DEVICE);
  }
  glutAddMenuEntry(_("Show/Hide..."), DIALOG_SHOWFILES);
  glutAddMenuEntry(_("Particle tracking..."), DIALOG_SHOOTER);
  glutAddMenuEntry(_("Time bounds..."), DIALOG_TIME);
  if(nterraininfo>0){
    glutAddMenuEntry(_("WUI display... ALT w..."), DIALOG_WUI);
  }

  /* --------------------------------window menu -------------------------- */

  CREATEMENU(windowdialogmenu, DialogMenu);
  glutAddMenuEntry(_("Fonts..."), DIALOG_FONTS);
  glutAddMenuEntry(_("User ticks..."), DIALOG_TICKS);
  glutAddMenuEntry(_("Labels..."), DIALOG_LABELS);
  glutAddMenuEntry(_("Properties..."), DIALOG_WINDOW);
  glutAddMenuEntry(_("Scaling..."), DIALOG_SCALING);

  /* --------------------------------dialog menu -------------------------- */

  CREATEMENU(dialogmenu,DialogMenu);

  glutAddMenuEntry(_("Clip scene...  ALT c"), DIALOG_CLIP);
  glutAddMenuEntry(_("Data bounds... ALT b"), DIALOG_BOUNDS);
  glutAddMenuEntry(_("Display...  ALT d"), DIALOG_DISPLAY);
  glutAddMenuEntry(_("Motion...  ALT m"),DIALOG_MOTION);
  glutAddMenuEntry(_("Viewpoints... ALT g"),DIALOG_VIEW);

  glutAddMenuEntry("-",MENU_DUMMY2);

  glutAddSubMenu(_("Data"), datadialogmenu);
  glutAddSubMenu(_("Files"), filesdialogmenu);
  glutAddSubMenu(_("View"), viewdialogmenu);
  glutAddSubMenu(_("Window"), windowdialogmenu);

  glutAddMenuEntry("-",MENU_DUMMY2);
  glutAddMenuEntry(_("Close all dialogs  ALT x"),DIALOG_HIDEALL);

  /* -------------------------------- font menu -------------------------- */

  if(showfontmenu==1){
    CREATEMENU(fontmenu,FontMenu);
    switch(fontindex){
    case SMALL_FONT:
      glutAddMenuEntry(_("*Normal"),SMALL_FONT);
      glutAddMenuEntry(_("Large"),LARGE_FONT);
#ifdef pp_BETA
      glutAddMenuEntry(_("Scaled"),SCALED_FONT);
#endif
      break;
    case LARGE_FONT:
      glutAddMenuEntry(_("Normal"),SMALL_FONT);
      glutAddMenuEntry(_("*Large"),LARGE_FONT);
#ifdef pp_BETA
      glutAddMenuEntry(_("Scaled"),SCALED_FONT);
#endif
      break;
#ifdef pp_BETA
    case SCALED_FONT:
      glutAddMenuEntry(_("Normal"),SMALL_FONT);
      glutAddMenuEntry(_("Large"),LARGE_FONT);
      glutAddMenuEntry(_("*Scaled"),SCALED_FONT);
#endif
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
  }

  /* -------------------------------- units menu -------------------------- */

  if(nunitclasses>0){
    for(i=0;i<nunitclasses;i++){
      f_units *uci;
      int j;

      uci = unitclasses + i;
      CheckMemory;

      CREATEMENU(uci->submenuid,UnitsMenu);

      for(j=0;j<uci->nunits;j++){
        char menulabel[1024];

        if(uci->unit_index==j){
          strcpy(menulabel,"*");
          strcat(menulabel,uci->units[j].unit);
        }
        else{
          strcpy(menulabel,uci->units[j].unit);
        }
        if(smokediff==1&&uci->diff_index==j&&uci->units[j].rel_defined==1){
          strcat(menulabel," rel to ");
          strcat(menulabel,uci->units[j].rel_val);
          strcat(menulabel," ");
          strcat(menulabel,uci->units[0].unit);
        }
        glutAddMenuEntry(menulabel,1000*i+j);
      }
    }
    CREATEMENU(unitsmenu,UnitsMenu);
    for(i=0;i<nunitclasses;i++){
      f_units *uci;

      uci = unitclasses + i;
      if(uci->visible==1||show_all_units==1||strcmp(uci->unitclass,"Distance")==0){
        glutAddSubMenu(uci->unitclass,uci->submenuid);
      }
    }
    if(vishmsTimelabel==0)glutAddMenuEntry(_("time (h:m:s)"), MENU_UNITS_HMS);
    if(vishmsTimelabel==1)glutAddMenuEntry(_("*time (h:m:s)"), MENU_UNITS_HMS);
#ifdef pp_BETA
    if(show_all_units==1)glutAddMenuEntry("*show all units", MENU_UNITS_SHOWALL);
    if(show_all_units==0)glutAddMenuEntry("show all units", MENU_UNITS_SHOWALL);
#endif
    glutAddMenuEntry(_("Reset"), MENU_UNITS_RESET);
  }

/* --------------------------------option menu -------------------------- */

  CREATEMENU(optionmenu,OptionMenu);
  if(nunitclasses>0)glutAddSubMenu(_("Display Units"),unitsmenu);
  glutAddSubMenu(_("Rotation type"),rotatetypemenu);
  glutAddSubMenu(_("Max frame rate"),frameratemenu);
  glutAddSubMenu(_("Render"),rendermenu);
  glutAddSubMenu(_("Tours"),tourmenu);
  if(showfontmenu==1)glutAddSubMenu(_("Font"),fontmenu);
  if(trainer_active==1)glutAddMenuEntry(_("Trainer menu"),MENU_OPTION_TRAINERMENU);

/* -------------------------------- about menu -------------------------- */

  CREATEMENU(disclaimermenu,AboutMenu);
  glutAddMenuEntry("The US Department of Commerce makes no warranty, expressed or",MENU_DUMMY);
  glutAddMenuEntry("implied, to users of Smokeview, and accepts no responsibility",MENU_DUMMY);
  glutAddMenuEntry("for its use. Users of Smokeview assume sole responsibility under",MENU_DUMMY);
  glutAddMenuEntry("Federal law for determining the appropriateness of its use in any",MENU_DUMMY);
  glutAddMenuEntry("particular application; for any conclusions drawn from the results",MENU_DUMMY);
  glutAddMenuEntry("of its use; and for any actions taken or not taken as a result of",MENU_DUMMY);
  glutAddMenuEntry("analysis performed using this tools.",MENU_DUMMY);
  glutAddMenuEntry("",MENU_DUMMY);
  glutAddMenuEntry("Smokeview and the companion program FDS is intended for use only",MENU_DUMMY);
  glutAddMenuEntry("by those competent in the fields of fluid dynamics, thermodynamics,",MENU_DUMMY);
  glutAddMenuEntry("combustion, and heat transfer, and is intended only to supplement",MENU_DUMMY);
  glutAddMenuEntry("the informed judgment of the qualified user. These software packages",MENU_DUMMY);
  glutAddMenuEntry("may or may not have predictive capability when applied to a specific",MENU_DUMMY);
  glutAddMenuEntry("set of factual circumstances.  Lack of accurate predictions could lead",MENU_DUMMY);
  glutAddMenuEntry("to erroneous conclusions with regard to fire safety.  All results",MENU_DUMMY);
  glutAddMenuEntry("should be evaluated by an informed user.",1);

/* -------------------------------- about menu -------------------------- */

  CREATEMENU(aboutmenu,AboutMenu);
  glutAddMenuEntry(release_title,1);
  {
#ifdef pp_GPU
    char version_label[256];
#endif
    char menulabel[1024];

    sprintf(menulabel,"  Smokeview (64 bit) build: %s",smv_githash);
    glutAddMenuEntry(menulabel,1);
    if(fds_version!=NULL){
      sprintf(menulabel, "  FDS version: %s", fds_version);
      glutAddMenuEntry(menulabel, 1);
    }
    if(fds_githash!=NULL){
      sprintf(menulabel,"  FDS build: %s",fds_githash);
      glutAddMenuEntry(menulabel,1);
    }
#ifdef pp_GPU
    strcpy(version_label,_("  OpenGL version:"));
    strcat(version_label," ");
    strcat(version_label,(char *)glGetString(GL_VERSION));
    glutAddMenuEntry(version_label,1);
    if(gpuactive==1){
      if(usegpu==1){
        strcpy(menulabel,_("  GPU activated. (Press G to deactivate)"));
      }
      else{
        strcpy(menulabel,_("  GPU available but not in use. (Press G to activate)"));
      }
    }
    else{
      strcpy(menulabel,_("  GPU not available"));
    }
    glutAddMenuEntry(menulabel,1);
#endif
#ifdef pp_CULL
    if(cullactive==1&&gpuactive==1){
      if(cullsmoke==1&&usegpu==1){
        strcpy(menulabel,_("  Smoke culling activated. (Press C to deactivate)"));
      }
      else{
        strcpy(menulabel,_("  Smoke culling available but not in use. ( To activate: "));
        strcat(menulabel,_(" Press"));
        if(usegpu==0 && cullsmoke==1)strcat(menulabel," G.)");
        if(usegpu==1 && cullsmoke==0)strcat(menulabel," C.)");
        if(usegpu==0 && cullsmoke==0)strcat(menulabel," G then C.)");
      }
    }
    else{
      strcpy(menulabel,_("  Smoke culling not available"));
    }
    glutAddMenuEntry(menulabel,1);
#endif
    glutAddSubMenu(_("Disclaimer"),disclaimermenu);
  }

  /* --------------------------------web help menu -------------------------- */

  CREATEMENU(webhelpmenu,HelpMenu);

#ifdef WIN32
  glutAddMenuEntry(_("Obtain Documentation"), MENU_HELP_DOCUMENTATION);
  glutAddMenuEntry(_("Report problems"), MENU_HELP_ISSUES);
  glutAddMenuEntry(_("Download software updates"), MENU_HELP_DOWNLOADS);
  glutAddMenuEntry(_("FDS/Smokeview website"), MENU_HELP_FDSWEB);
#endif
#ifdef pp_OSX
  glutAddMenuEntry(_("Obtain Documentation"),MENU_HELP_DOCUMENTATION);
  glutAddMenuEntry(_("Report problems"),MENU_HELP_ISSUES);
  glutAddMenuEntry(_("Download software updates"),MENU_HELP_DOWNLOADS);
  glutAddMenuEntry(_("FDS/Smokeview website"),MENU_HELP_FDSWEB);
#endif
#ifndef WIN32
#ifndef PP_OSX
  glutAddMenuEntry(_("Download documentation at  http://fire.nist.gov/fds/documentation.html"),MENU_DUMMY);
  glutAddMenuEntry(_("Report a problem at http://code.google.com/p/fds-smv/issues/"),MENU_DUMMY);
  glutAddMenuEntry(_("Check for updates at http://code.google.com/p/fds-smv/downloads/"),MENU_DUMMY);
  glutAddMenuEntry(_("FDS/Smokeview website: http://fire.nist.gov/fds"),MENU_DUMMY);
#endif
#endif

  /* --------------------------------keyboard help menu -------------------------- */

  CREATEMENU(keyboardhelpmenu,HelpMenu);
  if(plotstate==DYNAMIC_PLOTS){
    glutAddMenuEntry(_("Animation"),MENU_DUMMY);
    glutAddMenuEntry(_("  t: set/unset single time step mode"), MENU_DUMMY);
    glutAddMenuEntry(_("  0: reset animation to the initial time"), MENU_DUMMY);
    glutAddMenuEntry(_("  T: toggle method for interpolating data color"), MENU_DUMMY);
    if(cellcenter_slice_active==1){
      glutAddMenuEntry(_("     (also, toggles cell center display on/off)"), MENU_DUMMY);
      glutAddMenuEntry(_("  @: display FDS values in cell centered slices"), MENU_DUMMY);
    }
    glutAddMenuEntry(_("  u: reload files"), MENU_DUMMY);
    glutAddMenuEntry(_("  H: toggle  visibility of slice and vector slice files"), MENU_DUMMY);
    glutAddMenuEntry(_("  L: unload last slice file loaded"), MENU_DUMMY);
    glutAddMenuEntry(_("  1-9: number of frames to skip"), MENU_DUMMY);
  }
  if(rotation_type==EYE_CENTERED){
    glutAddMenuEntry(_("Motion"), MENU_DUMMY);
    glutAddMenuEntry(_("   left/right cursor: rotate left/right"), MENU_DUMMY);
    glutAddMenuEntry(_("      up/down cursor: move forward/backward"), MENU_DUMMY);
    glutAddMenuEntry(_(" CTRL:up/down cursor: move forward/backward 5 times slower"), MENU_DUMMY);
    glutAddMenuEntry(_(" SHIFT: left/right cursor: rotate 90 degrees"), MENU_DUMMY);
    glutAddMenuEntry(_("    ALT:left/right cursor: slide left/right"), MENU_DUMMY);
    glutAddMenuEntry(_("    ALT:   up/down cursor: slide up/down"), MENU_DUMMY);
    glutAddMenuEntry(_("     INSERT/HOME/PageUP: tilt down/reset/tilt up"), MENU_DUMMY);
  }
  if(plotstate==STATIC_PLOTS){
    glutAddMenuEntry(_("Plot3D"), MENU_DUMMY);
    glutAddMenuEntry(_("  x,y,z: toggle contour plot visibility along x, y and z axis"), MENU_DUMMY);
    glutAddMenuEntry(_("  p: increment plot3d variable"), MENU_DUMMY);
    glutAddMenuEntry(_("  P: toggle cursor key mappings"), MENU_DUMMY);
    glutAddMenuEntry(_("  v: toggle flow vector visiblity"), MENU_DUMMY);
    glutAddMenuEntry(_("  a/ALT a: increase/decrease flow vector length by 1.5"), MENU_DUMMY);
    glutAddMenuEntry(_("  s: change interval between adjacent vectors"), MENU_DUMMY);
    glutAddMenuEntry(_("  c: toggle between continuous and 2D stepped contours"), MENU_DUMMY);
    glutAddMenuEntry(_("  i: toggle iso-surface visibility"), MENU_DUMMY);
  }
  glutAddMenuEntry(_("Misc"), MENU_DUMMY);
  glutAddMenuEntry(_("  r: render the current scene to an image file"), MENU_DUMMY);
  glutAddMenuEntry(_("  R:   (same as r but at twice the resolution)"), MENU_DUMMY);
  if(ntotal_blockages>0||isZoneFireModel==1){
    glutAddMenuEntry(_("  g: toggle grid visibility"), MENU_DUMMY);
  }
  glutAddMenuEntry(_("  e: toggle between view rotation types: scene centered 2 axis, 1 axis, 3 axis and eye centered"), MENU_DUMMY);
  glutAddMenuEntry(_("  q: display blockage locations as specified by user or by FDS"), MENU_DUMMY);
  if(ntotal_blockages>0){
    glutAddMenuEntry(_("  O: toggle blockage view (normal <--> outline)"), MENU_DUMMY);
    glutAddMenuEntry(_("  ALT o: cycle between all blockage view types"), MENU_DUMMY);
  }
  if(ndeviceinfo>0&&GetNumActiveDevices()>0){
    glutAddMenuEntry("  j/ALT j: increase/decrease object size", MENU_DUMMY);
  }
  glutAddMenuEntry("  ALT r: toggle research mode (global min/max for coloring data, turn off axis label smoothing)", MENU_DUMMY);
  glutAddMenuEntry(_("  W: toggle clipping - use Options/Clip menu to specify clipping planes"), MENU_DUMMY);
  glutAddMenuEntry(_("  -/space bar: decrement/increment time step, 2D contour planes, 3D contour levels"), MENU_DUMMY);
  glutAddMenuEntry("", MENU_DUMMY);
  glutAddMenuEntry(_("  ALT v: toggle projection  method (between perspective and size preserving)"), MENU_DUMMY);
  if(n_embedded_meshes>0){
    glutAddMenuEntry(_("  ALT u: toggle coarse slice display in embedded mesh"), MENU_DUMMY);
  }
  if(cellcenter_slice_active==1){
    glutAddMenuEntry(_("  ALT y: if current slice is cell centered, toggle interpolation on/off"), MENU_DUMMY);
  }
  if(caseini_filename!=NULL&&strlen(caseini_filename)>0){
    char inilabel[512];

    sprintf(inilabel,_("  #: save settings to %s"),caseini_filename);
    glutAddMenuEntry(inilabel,MENU_DUMMY);
  }
  else{
    glutAddMenuEntry(_("  #: save settings (create casename.ini file)"), MENU_DUMMY);
  }
  glutAddMenuEntry(_("  !: snap scene to closest 45 degree orientation"), MENU_DUMMY);
  glutAddMenuEntry(_("  ~: level the scene"),2);
  glutAddMenuEntry(_("  &: toggle line anti-aliasing (draw lines smoothly)"), MENU_DUMMY);

  /* --------------------------------mouse help menu -------------------------- */

  CREATEMENU(mousehelpmenu,HelpMenu);
  switch(rotation_type){
    case ROTATION_2AXIS:
      glutAddMenuEntry(_("horizontal/vertical: rotate about z, x axis"), MENU_DUMMY);
      break;
    case ROTATION_1AXIS:
      glutAddMenuEntry(_("horizontal: rotate about z axis"), MENU_DUMMY);
      break;
    case ROTATION_3AXIS:
      glutAddMenuEntry(_("horizontal/vertical: rotate about z, x axis (click near scene center)"), MENU_DUMMY);
      glutAddMenuEntry(_("clock/counter clockwise: rotate about y axis (click near scene edge)"), MENU_DUMMY);
      break;
    case EYE_CENTERED:
      glutAddMenuEntry(_("horizontal/vertical: rotate about user location"), MENU_DUMMY);
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
  switch(rotation_type){
    case EYE_CENTERED:
      break;
    case ROTATION_2AXIS:
    case ROTATION_1AXIS:
    case ROTATION_3AXIS:
      glutAddMenuEntry(_("CTRL horizontal/vertical: translate along x, y axis"), MENU_DUMMY);
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
  glutAddMenuEntry(_("ALT vertical: translate along z axis"), MENU_DUMMY);
  if(SHOW_gslice_data==1){
    glutAddMenuEntry(_("double-click: rotate/translate 3D node-centered slice"), MENU_DUMMY);
  }

  /* --------------------------------help menu -------------------------- */

  CREATEMENU(helpmenu,HelpMenu);
  glutAddSubMenu(_("Web"),webhelpmenu);
  glutAddSubMenu(_("Shortcuts"),keyboardhelpmenu);
  glutAddSubMenu(_("Mouse"),mousehelpmenu);
  glutAddSubMenu(_("About"),aboutmenu);

  /* --------------------------------particle menu -------------------------- */

  if(npartinfo>0&&nevac!=npartinfo){
    int ii;

    CREATEMENU(unloadpartmenu,UnloadPartMenu);
    for(ii=0;ii<npartinfo;ii++){
      partdata *parti;
      char menulabel[1024];

      i = partorderindex[ii];
      parti = partinfo + i;
      if(parti->loaded==1&&parti->evac==0){
        STRCPY(menulabel,parti->menulabel);
        glutAddMenuEntry(menulabel,i);
      }
    }

    glutAddMenuEntry(_("Unload all"),MENU_UNLOADPARTICLE_UNLOADALL);

    if(nmeshes==1){
      CREATEMENU(particlemenu,LoadParticleMenu);
    }
    else{
      CREATEMENU(particlesubmenu,LoadParticleMenu);
    }
    for(ii=0;ii<npartinfo;ii++){
      char menulabel[1024];

      i = partorderindex[ii];
      if(partinfo[i].evac==1)continue;
      if(partinfo[i].loaded==1){
        STRCPY(menulabel,"*");
        STRCAT(menulabel,partinfo[i].menulabel);
      }
      else{
        STRCPY(menulabel,partinfo[i].menulabel);
      }
      glutAddMenuEntry(menulabel,i);
    }
    if(nmeshes>1){
      char menulabel[1024];

      CREATEMENU(particlemenu,LoadParticleMenu);
      if(npartinfo > 0){
        strcpy(menulabel, _("Particles"));
        strcat(menulabel, " - ");
        strcat(menulabel, _("All meshes"));
        glutAddMenuEntry(menulabel, MENU_PARTICLE_ALLMESHES);
        strcpy(menulabel, _("Particles"));
        strcat(menulabel, " - ");
        strcat(menulabel, _("Single mesh"));
        glutAddSubMenu(menulabel, particlesubmenu);
        glutAddMenuEntry("-", MENU_PARTICLE_DUMMY);
      }
    }

    if(npartloaded<=1){
      glutAddMenuEntry(_("Unload"),MENU_PARTICLE_UNLOAD);
    }
     else{
       glutAddSubMenu(_("Unload"),unloadpartmenu);
     }
  }

  if(nevac>0){
    int ii;

    CREATEMENU(unloadevacmenu,UnloadEvacMenu);
    for(ii=0;ii<npartinfo;ii++){
      partdata *parti;
      char menulabel[1024];

      i = partorderindex[ii];
      parti = partinfo + i;
      if(parti->loaded==1&&parti->evac==0){
        STRCPY(menulabel,parti->menulabel);
        glutAddMenuEntry(menulabel,i);
      }
    }
    glutAddMenuEntry(_("Unload all"),MENU_UNLOADEVAC_UNLOADALL);

    CREATEMENU(evacmenu,EvacMenu);
    {
      int nevacs=0,nevacloaded2=0;
      char menulabel[1024];

      for(ii=0;ii<npartinfo;ii++){
        partdata *parti;

        parti = partinfo + ii;
        if(parti->evac==1){
          if(parti->loaded==1)nevacloaded2++;
          nevacs++;
        }
      }

      if(nevacs>1){
        strcpy(menulabel,_("Humans - all meshes"));
        glutAddMenuEntry(menulabel,MENU_EVAC_ALLMESHES);
        glutAddMenuEntry("-",MENU_EVAC_DUMMY);
      }
      for(ii=0;ii<npartinfo;ii++){
        i = partorderindex[ii];
        if(partinfo[i].evac==0)continue;
        if(partinfo[i].loaded==1){
          STRCPY(menulabel,"*");
          STRCAT(menulabel,partinfo[i].menulabel);
        }
        else{
          STRCPY(menulabel,partinfo[i].menulabel);
        }
        glutAddMenuEntry(menulabel,i);
      }
      if(nevacloaded2<=1){
        glutAddMenuEntry(_("Unload"),MENU_EVAC_UNLOADALL);
      }
       else{
         glutAddSubMenu(_("Unload"),unloadevacmenu);
       }
  }
    }

/* --------------------------------unload and load vslice menus -------------------------- */

  if(nvsliceinfo>0){
    int ii;

    if(nmultivsliceinfo<nvsliceinfo){
      CREATEMENU(unloadmultivslicemenu,UnloadMultiVSliceMenu);
      for(i=0;i<nmultivsliceinfo;i++){
        multivslicedata *mvslicei;

        mvslicei = multivsliceinfo + i;
        if(mvslicei->loaded!=0){
          glutAddMenuEntry(mvslicei->menulabel2,i);
        }
      }
      glutAddMenuEntry(_("Unload all"),UNLOAD_ALL);

      nloadsubmvslicemenu=1;
      for(i=1;i<nmultivsliceinfo;i++){
        vslicedata *vi, *vim1;
        slicedata *si, *sim1;

        vi = vsliceinfo + (multivsliceinfo+i)->ivslices[0];
        vim1 = vsliceinfo + (multivsliceinfo+i-1)->ivslices[0];
        si = sliceinfo + vi->ival;
        sim1 = sliceinfo + vim1->ival;
        if(strcmp(si->label.longlabel,sim1->label.longlabel)!=0){
          nloadsubmvslicemenu++;
        }
      }
      NewMemory((void **)&loadsubmvslicemenu,nloadsubmvslicemenu*sizeof(int));
      for(i=0;i<nloadsubmvslicemenu;i++){
        loadsubmvslicemenu[i]=0;
      }

      nmultisliceloaded=0;
      nloadsubmvslicemenu=0;
      for(i=0;i<nmultivsliceinfo;i++){
        vslicedata *vi, *vim1;
        slicedata *si, *sim1;
        char menulabel[1024];
        multivslicedata *mvslicei;

        mvslicei = multivsliceinfo + i;

        vi = vsliceinfo + mvslicei->ivslices[0];
        si = sliceinfo + vi->ival;
        if(i>0){
          vim1 = vsliceinfo + (multivsliceinfo+i-1)->ivslices[0];
          sim1 = sliceinfo + vim1->ival;
        }
        if(i==0||(i!=0&&strcmp(si->label.longlabel,sim1->label.longlabel)!=0)){
          CREATEMENU(loadsubmvslicemenu[nloadsubmvslicemenu],LoadMultiVSliceMenu);
          nloadsubmvslicemenu++;
        }

        STRCPY(menulabel,"");
        if(mvslicei->loaded==1){
          STRCAT(menulabel,"*");
          nmultisliceloaded++;
        }
        else if(mvslicei->loaded==-1){
          STRCAT(menulabel,"#");
        }
        else{
        }
        STRCAT(menulabel,mvslicei->menulabel);
        if(si->slicelabel!=NULL){
          STRCAT(menulabel," - ");
          STRCAT(menulabel,si->slicelabel);
        }
        glutAddMenuEntry(menulabel,i);
      }

      nloadsubmvslicemenu=0;
      CREATEMENU(loadmultivslicemenu,LoadMultiVSliceMenu);
      for(i=0;i<nmultivsliceinfo;i++){
        vslicedata *vi, *vim1;
        slicedata *si, *sim1;

        vi = vsliceinfo + (multivsliceinfo+i)->ivslices[0];
        si = sliceinfo + vi->ival;
		    if(i>0){
		      vim1 = vsliceinfo + (multivsliceinfo+i-1)->ivslices[0];
          sim1 = sliceinfo + vim1->ival;
	    	}
        if(i==0||(i>0&&strcmp(si->label.longlabel,sim1->label.longlabel)!=0)){
          if(si->vec_comp==0||showallslicevectors==1){
            char mlabel[1024], mlabel2[1024];

            STRCPY(mlabel,si->label.longlabel);
            if(i==0&&si->mesh_type>0||(i>0&&si->mesh_type!=sim1->mesh_type)){
              sprintf(mlabel2,"*** Evac type %i meshes ***",si->mesh_type);
              glutAddMenuEntry(mlabel2,MENU_DUMMY);
            }
            glutAddSubMenu(mlabel,loadsubmvslicemenu[nloadsubmvslicemenu]);
          }
          nloadsubmvslicemenu++;
        }
      }
      if(nmultivsliceinfo>0)glutAddMenuEntry("-",MENU_DUMMY);
#ifdef pp_SLICEDUP
      if(nslicedups > 0){
        glutAddMenuEntry("Duplicate vector slices", MENU_DUMMY);
        if(vectorslicedup_option == SLICEDUP_KEEPALL){
          glutAddMenuEntry("  *keep all", MENU_KEEP_ALL);
        }
        else{
          glutAddMenuEntry("  keep all", MENU_KEEP_ALL);
        }
        if(vectorslicedup_option == SLICEDUP_KEEPFINE){
          glutAddMenuEntry("  *keep fine", MENU_KEEP_FINE);
        }
        else{
          glutAddMenuEntry("  keep fine", MENU_KEEP_FINE);
        }
        if(vectorslicedup_option == SLICEDUP_KEEPCOARSE){
          glutAddMenuEntry("  *keep coarse", MENU_KEEP_COARSE);
        }
        else{
          glutAddMenuEntry("  keep coarse", MENU_KEEP_COARSE);
        }
        glutAddMenuEntry("-", MENU_DUMMY);
      }
#endif

      if(showallslicevectors == 0)glutAddMenuEntry(_("Show all vector slice menu entries"), MENU_LOADVSLICE_SHOWALL);
      if(showallslicevectors == 1)glutAddMenuEntry(_("*Show all vector slice menu entries"), MENU_LOADVSLICE_SHOWALL);
      if(nmultisliceloaded>1){
        glutAddSubMenu(_("Unload"),unloadmultivslicemenu);
      }
      else{
        glutAddMenuEntry(_("Unload"),UNLOAD_ALL);
      }
    }

    CREATEMENU(unloadvslicemenu,UnloadVSliceMenu);
    for(ii=0;ii<nvsliceinfo;ii++){
      vslicedata *vd;

      i = vsliceorderindex[ii];
      vd = vsliceinfo + i;
      if(vd->loaded==0)continue;
      glutAddMenuEntry(vd->menulabel2,i);
    }
    glutAddMenuEntry("-",MENU_DUMMY);
    //glutAddMenuEntry("Unload last",-2);
    glutAddMenuEntry(_("Unload all"),UNLOAD_ALL);

    if(nvslice0>0){
      vslicedata *vd, *vdim1,*vdip1;
      if(nvsliceinfo>0){
        nloadsubvslicemenu=1;
        for(ii=1;ii<nvsliceinfo;ii++){
          slicedata *sd, *sdm1;

          i=vsliceorderindex[ii];
          vd = vsliceinfo + i;
          sd = sliceinfo + vd->ival;
		      if(ii>0){
            vdim1 = vsliceinfo + vsliceorderindex[ii-1];
            sdm1 = sliceinfo + vdim1->ival;
		      }
          if(ii==0||strcmp(sd->label.longlabel,sdm1->label.longlabel)!=0){
            nloadsubvslicemenu++;
          }
        }
        NewMemory((void **)&loadsubvslicemenu,nloadsubvslicemenu*sizeof(int));
        for(i=0;i<nloadsubvslicemenu;i++){
          loadsubvslicemenu[i]=0;
        }
        nloadsubvslicemenu=0;
        for(ii=0;ii<nvsliceinfo;ii++){
          slicedata *sd, *sdm1, *sdp1;
          char menulabel[1024];

          i=vsliceorderindex[ii];
          vd = vsliceinfo + i;
          sd = sliceinfo + vd->ival;

          if(ii!=0){
            vdim1 = vsliceinfo + vsliceorderindex[ii-1];
            sdm1 = sliceinfo + vdim1->ival;
          }
          if(ii!=nvsliceinfo-1){
            vdip1 = vsliceinfo + vsliceorderindex[ii+1];
            sdp1 = sliceinfo + vdip1->ival;
          }

          if(ii==0||strcmp(sd->label.longlabel,sdm1->label.longlabel)!=0){
            CREATEMENU(loadsubvslicemenu[nloadsubvslicemenu],LoadVSliceMenu);
          }
          if(vd->loaded==1){
            STRCPY(menulabel,"*");
            STRCAT(menulabel,sd->menulabel);
          }
          else{
            STRCPY(menulabel,sd->menulabel);
          }
          if(sd->vec_comp==0||showallslicevectors==1)glutAddMenuEntry(menulabel,i);
          if(ii==nvsliceinfo-1||strcmp(sd->label.longlabel,sdp1->label.longlabel)!=0){
            subvslice_menuindex[nloadsubvslicemenu]=vsliceorderindex[ii];
            if(sd->ndirxyz[1]+sd->ndirxyz[2]+sd->ndirxyz[3]>1){
              glutAddMenuEntry("-",MENU_DUMMY);
            }
            if(sd->ndirxyz[1]>1){
              glutAddMenuEntry(_("Load All x"),-1000-4*nloadsubvslicemenu-1);
            }
            if(sd->ndirxyz[2]>1){
              glutAddMenuEntry(_("Load All y"),-1000-4*nloadsubvslicemenu-2);
            }
            if(sd->ndirxyz[3]>1){
              glutAddMenuEntry(_("Load All z"),-1000-4*nloadsubvslicemenu-3);
            }
            if(sd->ndirxyz[1]+sd->ndirxyz[2]+sd->ndirxyz[3]>1){
              glutAddMenuEntry(_("Load All"),-1000-4*nloadsubvslicemenu);
            }
          }
          if(ii==0||strcmp(sd->label.longlabel,sdm1->label.longlabel)!=0){
            nloadsubvslicemenu++;
          }
        }
        CREATEMENU(vslicemenu,LoadVSliceMenu);
        nloadsubvslicemenu=0;
        for(ii=0;ii<nvsliceinfo;ii++){
          slicedata *sd, *sdm1;

          i=vsliceorderindex[ii];
          vd = vsliceinfo + i;
          sd = sliceinfo + vd->ival;
		  if(ii>0){
            vdim1 = vsliceinfo + vsliceorderindex[ii-1];
            sdm1 = sliceinfo + vdim1->ival;
		  }
          if(ii==0||strcmp(sd->label.longlabel,sdm1->label.longlabel)!=0){
            if(sd->vec_comp==0||showallslicevectors==1){
              char mlabel[1024], mlabel2[1024];

              STRCPY(mlabel,sd->label.longlabel);
              if(ii==0&&sd->mesh_type>0||(ii>0&&sd->mesh_type!=sdm1->mesh_type)){
                sprintf(mlabel2,"*** Evac type %i meshdata ***",sd->mesh_type);
                glutAddMenuEntry(mlabel2,MENU_DUMMY);
              }
              glutAddSubMenu(mlabel,loadsubvslicemenu[nloadsubvslicemenu]);
            }
            nloadsubvslicemenu++;
          }
        }
      }
    }
    if(nvsliceinfo>0)glutAddMenuEntry("-",MENU_DUMMY);
    if(showallslicevectors==0)glutAddMenuEntry(_("Show all vector slice menu entries"), MENU_LOADVSLICE_SHOWALL);
    if(showallslicevectors==1)glutAddMenuEntry(_("*Show all vector slice menu entries"), MENU_LOADVSLICE_SHOWALL);
    if(nvsliceloaded>1){
      glutAddSubMenu(_("Unload"),unloadvslicemenu);
    }
    else{
     glutAddMenuEntry(_("Unload"),UNLOAD_ALL);
    }
  }

  /* --------------------------------unload and load slice menus -------------------------- */

  if(nterraininfo>0){
    int nterrainloaded=0;

    CREATEMENU(unloadterrainmenu,UnloadTerrainMenu);
    for(i=0;i<nterraininfo;i++){
      terraindata *terri;

      terri = terraininfo + i;
      if(terri->loaded==1){
        nterrainloaded++;
        glutAddMenuEntry(terri->file,i);
      }
    }
    if(nterrainloaded>1){
        glutAddMenuEntry("-",MENU_UNLOADTERRAIN_DUMMY);
        glutAddMenuEntry(_("Unload all"), MENU_UNLOADTERRAIN_UNLOADALL);
    }
    CREATEMENU(loadterrainmenu,LoadTerrainMenu);
    if(nterraininfo>1){
      glutAddMenuEntry(_("All terrains"), MENU_LOADTERRAIN_LOADALL);
      glutAddMenuEntry("-",MENU_LOADTERRAIN_DUMMY);
    }
    /*
    leaving code commented in case I later decide to load/unload terrain files
    for(i=0;i<nterraininfo;i++){
      char menulabel[256];

      terraindata *terri;

      terri = terraininfo + i;
      strcpy(menulabel,"");
      if(terri->loaded==1)strcat(menulabel,"*");
      strcat(menulabel,terri->file);
      glutAddMenuEntry(menulabel,i);
    }
    */
    if(nterrainloaded==1){
      glutAddMenuEntry("-",MENU_LOADTERRAIN_DUMMY);
      glutAddMenuEntry(_("Unload terrain"), MENU_LOADTERRAIN_UNLOAD);
    }
    else if(nterrainloaded>1){
      glutAddMenuEntry("-",MENU_LOADTERRAIN_DUMMY);
      glutAddSubMenu(_("Unload terrain"),unloadterrainmenu);
    }
  }
    if(nsliceinfo>0){

      if(nmultisliceinfo<nsliceinfo){
        CREATEMENU(unloadmultislicemenu,UnloadMultiSliceMenu);
        nmultisliceloaded=0;
        for(i=0;i<nmultisliceinfo;i++){
          multislicedata *mslicei;

          mslicei = multisliceinfo + i;
          if(mslicei->loaded!=0){
            glutAddMenuEntry(mslicei->menulabel2,i);
          }
        }
        glutAddMenuEntry(_("Unload all"),UNLOAD_ALL);

        nloadsubmslicemenu=1;
        for(i=1;i<nmultisliceinfo;i++){
          slicedata *sd, *sdim1;

          sd = sliceinfo+(multisliceinfo + i)->islices[0];
          sdim1 = sliceinfo+(multisliceinfo + i-1)->islices[0];
          if(strcmp(sd->label.longlabel,sdim1->label.longlabel)!=0)nloadsubmslicemenu++;
        }
        NewMemory((void **)&loadsubmslicemenu,nloadsubmslicemenu*sizeof(int));
        for(i=0;i<nloadsubmslicemenu;i++){
          loadsubmslicemenu[i]=0;
        }
        nloadsubmslicemenu=0;
        for(i=0;i<nmultisliceinfo;i++){
          slicedata *sd, *sdim1;
          char menulabel[1024];
          multislicedata *mslicei;

          sd = sliceinfo+(multisliceinfo + i)->islices[0];
          if(i>0)sdim1 = sliceinfo+(multisliceinfo + i-1)->islices[0];
          mslicei = multisliceinfo + i;
          if(i==0||strcmp(sd->label.longlabel,sdim1->label.longlabel)!=0){
            CREATEMENU(loadsubmslicemenu[nloadsubmslicemenu],LoadMultiSliceMenu);
            nloadsubmslicemenu++;
          }
          STRCPY(menulabel,"");
          if(mslicei->loaded==1){
            STRCAT(menulabel,"*");
            nmultisliceloaded++;
          }
          else if(mslicei->loaded==-1){
            STRCAT(menulabel,"#");
            nmultisliceloaded++;
          }
          STRCAT(menulabel,mslicei->menulabel);
          if(sd->slicelabel!=NULL){
            STRCAT(menulabel," - ");
            STRCAT(menulabel,sd->slicelabel);
          }
          glutAddMenuEntry(menulabel,i);
        }
        CREATEMENU(loadmultislicemenu,LoadMultiSliceMenu);
        nloadsubmslicemenu=0;
        for(i=0;i<nmultisliceinfo;i++){
          slicedata *sd, *sdim1;

          sd = sliceinfo+(multisliceinfo + i)->islices[0];
          if(i>0)sdim1 = sliceinfo+(multisliceinfo + i-1)->islices[0];

          if(i==0||strcmp(sd->label.longlabel,sdim1->label.longlabel)!=0){
            char mlabel[1024], mlabel2[1024];

            STRCPY(mlabel,sd->label.longlabel);
            if(i==0&&sd->mesh_type>0||(i>0&&sd->mesh_type!=sdim1->mesh_type)){
              sprintf(mlabel2,"*** Evac type %i meshes ***",sd->mesh_type);
              if(sd->slicetype==SLICE_CELL_CENTER){
                flowlabels *label;

                label = &sd->label;
                if(strcmp(label->shortlabel,"U-VEL")==0||strcmp(label->shortlabel,"V-VEL")==0||strcmp(label->shortlabel,"W-VEL")==0){
                  continue;
                }
              }
              glutAddMenuEntry(mlabel2,MENU_DUMMY);
            }
            if(sd->slicetype==SLICE_CELL_CENTER){
              flowlabels *label;

              label = &sd->label;
              if(strcmp(label->shortlabel,"U-VEL")==0||strcmp(label->shortlabel,"V-VEL")==0||strcmp(label->shortlabel,"W-VEL")==0){
                continue;
              }
            }
            glutAddSubMenu(mlabel,             loadsubmslicemenu[nloadsubmslicemenu]);
            nloadsubmslicemenu++;
          }
        }
        if(nmultisliceinfo>0)glutAddMenuEntry("-",MENU_DUMMY);
#ifdef pp_SLICEDUP
        if(nslicedups > 0){
          glutAddMenuEntry("Duplicate slices", MENU_DUMMY);
          if(slicedup_option == SLICEDUP_KEEPALL){
            glutAddMenuEntry("  *keep all", MENU_KEEP_ALL);
          }
          else{
            glutAddMenuEntry("  keep all", MENU_KEEP_ALL);
          }
          if(slicedup_option == SLICEDUP_KEEPFINE){
            glutAddMenuEntry("  *keep fine", MENU_KEEP_FINE);
          }
          else{
            glutAddMenuEntry("  keep fine", MENU_KEEP_FINE);
          }
          if(slicedup_option == SLICEDUP_KEEPCOARSE){
            glutAddMenuEntry("  *keep coarse", MENU_KEEP_COARSE);
          }
          else{
            glutAddMenuEntry("  keep coarse", MENU_KEEP_COARSE);
          }
          glutAddMenuEntry("-", MENU_DUMMY);
        }
#endif
        if(use_set_slicecolor==1){
          glutAddMenuEntry("  *defer slice coloring", MENU_SLICECOLORDEFER);
        }
        else{
          glutAddMenuEntry("  defer slice coloring", MENU_SLICECOLORDEFER);
        }
        if(nmultisliceloaded>1){
          glutAddSubMenu(_("Unload"),unloadmultislicemenu);
        }
        else{
          glutAddMenuEntry(_("Unload"),UNLOAD_ALL);
        }

      }
      CREATEMENU(unloadslicemenu,UnloadSliceMenu);
      for(i=0;i<nsliceinfo;i++){
        slicedata *sd;
        char menulabel[1024];

        sd = sliceinfo + sliceorderindex[i];
        if(sd->loaded==1){
          STRCPY(menulabel,sd->menulabel2);
          glutAddMenuEntry(menulabel,sliceorderindex[i]);
        }
      }
      glutAddMenuEntry("-",MENU_DUMMY);
      glutAddMenuEntry(_("Unload last"),UNLOAD_LAST);
      glutAddMenuEntry(_("Unload all"),UNLOAD_ALL);

//*** this is where I would put the "sub-slice" menus ordered by type
      nloadsubslicemenu=1;
      for(i=1;i<nsliceinfo;i++){
        slicedata *sd,*sdim1;

        sd = sliceinfo + sliceorderindex[i];
        sdim1 = sliceinfo + sliceorderindex[i-1];
        if(strcmp(sd->label.longlabel,sdim1->label.longlabel)!=0)nloadsubslicemenu++;
      }
      NewMemory((void **)&loadsubslicemenu,nloadsubslicemenu*sizeof(int));
      for(i=0;i<nloadsubslicemenu;i++){
        loadsubslicemenu[i]=0;
      }
      iloadsubslicemenu=0;
      for(i=0;i<nsliceinfo;i++){
        slicedata *sd,*sdim1,*sdip1;
        char menulabel[1024];

        if(i!=0){
          sdim1 = sliceinfo + sliceorderindex[i-1];
        }
        sd = sliceinfo + sliceorderindex[i];
        if(i!=nsliceinfo-1){
          sdip1 = sliceinfo + sliceorderindex[i+1];
        }
        if(i==0||strcmp(sd->label.longlabel,sdim1->label.longlabel)!=0){
          CREATEMENU(loadsubslicemenu[iloadsubslicemenu],LoadSliceMenu);
        }
        STRCPY(menulabel,"");
        if(sd->loaded==1){
          STRCAT(menulabel,"*");
        }
        STRCAT(menulabel,sd->menulabel);
        if(sd->slicelabel!=NULL){
          STRCAT(menulabel," - ");
          STRCAT(menulabel,sd->slicelabel);
        }
        if(sd->menu_show==1){
          glutAddMenuEntry(menulabel,sliceorderindex[i]);
        }
        if(i==nsliceinfo-1||strcmp(sd->label.longlabel,sdip1->label.longlabel)!=0){
          subslice_menuindex[iloadsubslicemenu]=sliceorderindex[i];
  		    if(sd->ndirxyz[1]+sd->ndirxyz[2]+sd->ndirxyz[3]>1){
            glutAddMenuEntry("-",MENU_DUMMY);
		      }
		      if(sd->ndirxyz[1]>1){
            glutAddMenuEntry(_("Load All x"),-1000-4*iloadsubslicemenu-1);
		      }
		      if(sd->ndirxyz[2]>1){
            glutAddMenuEntry(_("Load All y"),-1000-4*iloadsubslicemenu-2);
		      }
		      if(sd->ndirxyz[3]>1){
            glutAddMenuEntry(_("Load All z"),-1000-4*iloadsubslicemenu-3);
		      }
		      if(sd->ndirxyz[1]+sd->ndirxyz[2]+sd->ndirxyz[3]>1){
            glutAddMenuEntry(_("Load All"),-1000-4*iloadsubslicemenu);
		      }
        }
        if(i==0||strcmp(sd->label.longlabel,sdim1->label.longlabel)!=0){
          iloadsubslicemenu++;
        }
      }
      CREATEMENU(loadslicemenu,LoadSliceMenu);
      iloadsubslicemenu=0;
      for(i=0;i<nsliceinfo;i++){
        slicedata *sd,*sdim1;

        sd = sliceinfo + sliceorderindex[i];
        if(i>0)sdim1 = sliceinfo + sliceorderindex[i-1];
        if(i==0||strcmp(sd->label.longlabel,sdim1->label.longlabel)!=0){
          char mlabel[1024],mlabel2[1024];;

          STRCPY(mlabel,sd->label.longlabel);
          if(i==0&&sd->mesh_type>0||(i>0&&sd->mesh_type!=sdim1->mesh_type)){
            if(sd->menu_show==1){
              sprintf(mlabel2,"*** Evac type %i meshdata ***",sd->mesh_type);
              glutAddMenuEntry(mlabel2,MENU_DUMMY);
            }
          }
          if(sd->menu_show==1)glutAddSubMenu(mlabel,loadsubslicemenu[iloadsubslicemenu]);
          iloadsubslicemenu++;
        }
      }
      glutAddMenuEntry("-", MENU_DUMMY);
      if(show_slice_in_obst == 1)glutAddMenuEntry("*Show slice in blockage", MENU_SHOWSLICE_INBLOCKAGE);
      if(show_slice_in_obst == 0)glutAddMenuEntry("Show slice in blockage", MENU_SHOWSLICE_INBLOCKAGE);
      glutAddMenuEntry("-",MENU_DUMMY);
      if(nsliceloaded>1){
        glutAddSubMenu(_("Unload"),unloadslicemenu);
      }
      else{
        glutAddMenuEntry(_("Unload"),UNLOAD_ALL);
      }
    }

/* --------------------------------unload and load 3d vol smoke menus -------------------------- */

    if(nvolsmoke3dloaded>0){
      CREATEMENU(unloadvolsmoke3dmenu,UnLoadVolSmoke3DMenu);
      if(nvolsmoke3dloaded>1){
        char vlabel[256];

        strcpy(vlabel,_("3D smoke (volume rendered)"));
        strcat(vlabel,_(" - All meshes"));
        glutAddMenuEntry(vlabel,UNLOAD_ALL);
      }
      for(i=0;i<nmeshes;i++){
        meshdata *meshi;
        volrenderdata *vr;

        meshi = meshinfo + i;
        vr = &(meshi->volrenderinfo);
        if(vr->fireslice==NULL||vr->smokeslice==NULL)continue;
        if(vr->loaded==0)continue;
        glutAddMenuEntry(meshi->label,i);
      }
    }
    if(nvolrenderinfo>0){
      CREATEMENU(loadvolsmoke3dmenu,LoadVolSmoke3DMenu);
      if(nvolrenderinfo>1){
        char vlabel[256];

        strcpy(vlabel,_("3D smoke (volume rendered)"));
        strcat(vlabel,_(" - All meshes"));
        glutAddMenuEntry(vlabel,LOAD_ALL);
        glutAddMenuEntry("-",MENU_DUMMY);
      }
      for(i=0;i<nmeshes;i++){
        meshdata *meshi;
        volrenderdata *vr;
        char menulabel[1024];

        meshi = meshinfo + i;
        vr = &(meshi->volrenderinfo);
        if(vr->fireslice==NULL||vr->smokeslice==NULL)continue;
        strcpy(menulabel,"");
        if(vr->loaded==1)strcat(menulabel,"*");
        strcat(menulabel,meshi->label);
        glutAddMenuEntry(menulabel,i);
      }
      if(nvolsmoke3dloaded==1)glutAddMenuEntry(_("Unload"),UNLOAD_ALL);
      if(nvolsmoke3dloaded>1)glutAddSubMenu(_("Unload"),unloadvolsmoke3dmenu);
    }

    /* --------------------------------unload and load 3d smoke menus -------------------------- */

    {
      smoke3ddata *smoke3di;
      if(nsmoke3dloaded>0){
        CREATEMENU(unloadsmoke3dmenu,UnLoadSmoke3DMenu);
        {
          int nsootloaded=0,nhrrloaded=0,nwaterloaded=0;

          for(i=0;i<nsmoke3dinfo;i++){
            smoke3di=smoke3dinfo + i;
            if(smoke3di->loaded==0)continue;
            switch(smoke3di->type){
            case SOOT:
              nsootloaded++;
              break;
            case FIRE:
              nhrrloaded++;
              break;
            case WATER:
              nwaterloaded++;
              break;
            default:
              ASSERT(FFALSE);
              break;
            }
          }
          if(nsootloaded>1) glutAddMenuEntry(_("SOOT MASS FRACTION - all meshes"), MENU_UNLOADSMOKE3D_UNLOADALLSOOT);
          if(nhrrloaded>1)  glutAddMenuEntry(_("HRRPUV - all meshes"), MENU_UNLOADSMOKE3D_UNLOADALLFIRE);
          if(nwaterloaded>1)glutAddMenuEntry(_("water - all meshes"), MENU_UNLOADSMOKE3D_UNLOADALLWATER);
          if(nsootloaded>1||nhrrloaded>1||nwaterloaded>1)glutAddMenuEntry("-",MENU_DUMMY);
        }
        for(i=0;i<nsmoke3dinfo;i++){
          smoke3di=smoke3dinfo+i;
          if(smoke3di->loaded==0)continue;
          glutAddMenuEntry(smoke3di->menulabel,i);
        }
      }
    }
    {
      smoke3ddata *smoke3di;
      int n_soot_menu=0, n_hrr_menu=0, n_water_menu=0;

      if(nsmoke3dinfo>0){
        char menulabel[1024];

        if(nmeshes==1){
          CREATEMENU(loadsmoke3dmenu,LoadSmoke3DMenu);
        }
        if(nmeshes>1){
          CREATEMENU(loadsmoke3dsootmenu,LoadSmoke3DMenu);
        }
        for(i=0;i<nsmoke3dinfo;i++){
          smoke3di = smoke3dinfo + i;
          if(smoke3di->type!=SOOT)continue;
          n_soot_menu++;
          strcpy(menulabel,"");
          if(smoke3di->loaded==1){
            strcat(menulabel,"*");
          }
          strcat(menulabel,smoke3di->menulabel);
          glutAddMenuEntry(menulabel,i);
        }
        if(nmeshes>1){
          CREATEMENU(loadsmoke3dhrrmenu,LoadSmoke3DMenu);
        }
        for(i=0;i<nsmoke3dinfo;i++){
          smoke3di = smoke3dinfo + i;
          if(smoke3di->type!=FIRE)continue;
          n_hrr_menu++;
          strcpy(menulabel,"");
          if(smoke3di->loaded==1){
            strcat(menulabel,"*");
          }
          strcat(menulabel,smoke3di->menulabel);
          glutAddMenuEntry(menulabel,i);
        }
        if(nmeshes>1){
          CREATEMENU(loadsmoke3dwatermenu,LoadSmoke3DMenu);
        }
        for(i=0;i<nsmoke3dinfo;i++){
          smoke3di = smoke3dinfo + i;
          if(smoke3di->type!=WATER)continue;
          n_water_menu++;
          strcpy(menulabel,"");
          if(smoke3di->loaded==1){
            strcat(menulabel,"*");
          }
          strcat(menulabel,smoke3di->menulabel);
          glutAddMenuEntry(menulabel,i);
        }
        if(nmeshes>1){
          CREATEMENU(loadsmoke3dmenu,LoadSmoke3DMenu);
        }
        {
          int useitem;
          smoke3ddata *smoke3dj;

          if(nmeshes>1){
            for(i=0;i<nsmoke3dinfo;i++){
              int j;

              useitem=i;
              smoke3di=smoke3dinfo + i;
              for(j=0;j<i;j++){
                smoke3dj = smoke3dinfo + j;
                if(strcmp(smoke3di->label.longlabel,smoke3dj->label.longlabel)==0){
                  useitem=-1;
                  break;
                }
              }
              if(useitem!=-1){
                strcpy(menulabel,smoke3di->label.longlabel);
                strcat(menulabel," - ");
                strcat(menulabel,_("All meshes"));
                glutAddMenuEntry(menulabel,-useitem-10);
              }
            }
            glutAddMenuEntry("-",MENU_DUMMY3);
          }
          if(nmeshes>1){
            if(n_soot_menu>0)glutAddSubMenu(_("SOOT MASS FRACTION - single mesh"),loadsmoke3dsootmenu);
            if(n_hrr_menu>0)glutAddSubMenu(_("HRRPUV - single mesh"),loadsmoke3dhrrmenu);
            if(n_water_menu>0)glutAddSubMenu(_("Water - single mesh"),loadsmoke3dwatermenu);
          }
        }
        if(use_iblank==0){
          glutAddMenuEntry("-", MENU_DUMMY3);
          glutAddMenuEntry(_("Initialize smoke blockage info"), MENU_SMOKE3D_IBLANK);
          if(nsmoke3dloaded>=1)glutAddMenuEntry("-", MENU_DUMMY3);
        }
        if(nsmoke3dloaded==1)glutAddMenuEntry(_("Unload"),UNLOAD_ALL);
        if(nsmoke3dloaded>1)glutAddSubMenu(_("Unload"),unloadsmoke3dmenu);
      }
    }

/* --------------------------------plot3d menu -------------------------- */

    if(nplot3dinfo>0){
      plot3ddata *plot3dim1, *plot3di;
      char menulabel[1024];
      int ii;

      CREATEMENU(unloadplot3dmenu,UnloadPlot3dMenu);
      for(ii=0;ii<nplot3dinfo;ii++){
        i=plot3dorderindex[ii];
        plot3di = plot3dinfo + i;
        if(ii==0){
          strcpy(menulabel,plot3di->longlabel);
          glutAddMenuEntry(menulabel,MENU_PLOT3D_DUMMY);
        }
        if(ii!=0&&strcmp(plot3di->longlabel,plot3dinfo[plot3dorderindex[ii-1]].longlabel)!=0){
          glutAddMenuEntry(plot3di->longlabel,MENU_PLOT3D_DUMMY);
        }
        if(plot3di->loaded==0)continue;
        STRCPY(menulabel,plot3dinfo[i].menulabel);
        glutAddMenuEntry(menulabel,i);
      }
      glutAddMenuEntry("Unload all",UNLOAD_ALL);



      nloadsubplot3dmenu=1;
      for(ii=1;ii<nplot3dinfo;ii++){
        int im1;

        i = plot3dorderindex[ii];
        im1 = plot3dorderindex[ii-1];
        plot3di = plot3dinfo + i;
        plot3dim1 = plot3dinfo + im1;
        if(ABS(plot3di->time-plot3dim1->time)>0.1)nloadsubplot3dmenu++;
      }
      NewMemory((void **)&loadsubplot3dmenu,nloadsubplot3dmenu*sizeof(int));
      for(i=0;i<nloadsubplot3dmenu;i++){
        loadsubplot3dmenu[i]=0;
      }

      nloadsubplot3dmenu=0;
      i = plot3dorderindex[0];
      plot3di = plot3dinfo + i;
      CREATEMENU(loadsubplot3dmenu[nloadsubplot3dmenu],LoadPlot3dMenu);
      strcpy(menulabel,"");
      if(plot3di->loaded==1){
        strcat(menulabel,"*");
      }
      strcat(menulabel,plot3di->menulabel);
      glutAddMenuEntry(menulabel,i);
      nloadsubplot3dmenu++;

      for(ii=1;ii<nplot3dinfo;ii++){
        int im1;

        i = plot3dorderindex[ii];
        im1 = plot3dorderindex[ii-1];
        plot3di = plot3dinfo + i;
        plot3dim1 = plot3dinfo + im1;
        if(ABS(plot3di->time-plot3dim1->time)>0.1){
          if(nmeshes>1)glutAddMenuEntry("  All meshes",-100000+nloadsubplot3dmenu-1);
          CREATEMENU(loadsubplot3dmenu[nloadsubplot3dmenu],LoadPlot3dMenu);
          nloadsubplot3dmenu++;
        }
        strcpy(menulabel,"");
        if(plot3di->loaded==1){
          strcat(menulabel,"*");
        }
        strcat(menulabel,plot3di->menulabel);
        glutAddMenuEntry(menulabel,i);
      }
      if(nmeshes>1)glutAddMenuEntry(_("  All meshes"),-100000+nloadsubplot3dmenu-1);

      nloadsubplot3dmenu=0;
      CREATEMENU(loadplot3dmenu,LoadPlot3dMenu);
      for(ii=0;ii<nplot3dinfo;ii++){
        int im1;

        i = plot3dorderindex[ii];
        plot3di = plot3dinfo + i;
        if(ii==0){
          strcpy(menulabel,plot3di->longlabel);
          glutAddMenuEntry(menulabel,MENU_PLOT3D_DUMMY);
          sprintf(menulabel,"  %f",plot3di->time);
          TrimZeros(menulabel);
          strcat(menulabel," s");
          if(nmeshes>1){
            glutAddSubMenu(menulabel,loadsubplot3dmenu[nloadsubplot3dmenu]);
          }
          else{
            strcpy(menulabel,"  ");
            if(plot3di->loaded==1){
              strcat(menulabel,"*");
            }
            strcat(menulabel,plot3di->menulabel);
            glutAddMenuEntry(menulabel,i);
          }
          nloadsubplot3dmenu++;
        }
        if(ii!=0){
          i = plot3dorderindex[ii];
          im1 = plot3dorderindex[ii-1];
          plot3di = plot3dinfo + i;
          plot3dim1 = plot3dinfo + im1;
          if(strcmp(plot3di->longlabel,plot3dim1->longlabel)!=0){
            glutAddMenuEntry(plot3di->longlabel,MENU_PLOT3D_DUMMY);
          }
          if(ABS(plot3di->time-plot3dim1->time)>0.1){
            sprintf(menulabel,"  %f",plot3di->time);
            TrimZeros(menulabel);
            strcat(menulabel," s");
            if(nmeshes>1){
              glutAddSubMenu(menulabel,loadsubplot3dmenu[nloadsubplot3dmenu]);
            }
            else{
              strcpy(menulabel,"  ");
              if(plot3di->loaded==1){
                strcat(menulabel,"*");
              }
              strcat(menulabel,plot3di->menulabel);
              glutAddMenuEntry(menulabel,i);
            }
            nloadsubplot3dmenu++;
          }
        }
      }
      if(nplot3dloaded>1){
        glutAddSubMenu(_("Unload"),unloadplot3dmenu);
      }
      else{
       glutAddMenuEntry(_("Unload"),UNLOAD_ALL);
      }
    }

/* --------------------------------load patch menu -------------------------- */

    if(npatchinfo>0){
      int ii;

      nloadpatchsubmenus=0;
      CREATEMENU(unloadpatchmenu,UnloadPatchMenu);
      for(ii=0;ii<npatchinfo;ii++){
        patchdata *patchi;
        char menulabel[1024];

        i = patchorderindex[ii];
        patchi = patchinfo + i;
        if(patchi->loaded==0)continue;
        STRCPY(menulabel,patchi->menulabel);
        glutAddMenuEntry(menulabel,i);
      }
      glutAddMenuEntry(_("Unload all"),UNLOAD_ALL);

      if(nmeshes>1&&loadpatchsubmenus==NULL){
        NewMemory((void **)&loadpatchsubmenus,npatchinfo*sizeof(int));
      }

      if(nmeshes>1){
        CREATEMENU(loadpatchsubmenus[nloadpatchsubmenus],LoadPatchMenu);
        nloadpatchsubmenus++;
      }
      else{
        CREATEMENU(loadpatchmenu,LoadPatchMenu);
      }

      for(ii=0;ii<npatchinfo;ii++){
        patchdata *patchim1, *patchi;
        char menulabel[1024];

        i = patchorderindex[ii];
        patchi = patchinfo + i;
        if(ii>0){
          patchim1 = patchinfo + patchorderindex[ii-1];
          if(nmeshes>1&&strcmp(patchim1->label.longlabel,patchi->label.longlabel)!=0){
            CREATEMENU(loadpatchsubmenus[nloadpatchsubmenus],LoadPatchMenu);
            nloadpatchsubmenus++;
          }
        }

        if(patchi->loaded==1){
          STRCPY(menulabel,"*");
          STRCAT(menulabel,patchi->menulabel);
        }
        else{
          STRCPY(menulabel,patchi->menulabel);
        }
        glutAddMenuEntry(menulabel,i);
      }

      if(nmeshes>1){
        CREATEMENU(loadpatchmenu,LoadPatchMenu);
      }

      {
        int useitem;
        patchdata *patchi, *patchj;

        if(nmeshes>1){
          char menulabel[1024];

          for(i=0;i<npatchinfo;i++){
            int j;

            useitem=i;
            patchi = patchinfo + i;
            for(j=0;j<i;j++){
              patchj = patchinfo + j;
              if(strcmp(patchi->label.longlabel,patchj->label.longlabel)==0){
                useitem=-1;
                break;
              }
            }
            if(useitem!=-1){
              strcpy(menulabel,patchi->label.longlabel);
              strcat(menulabel," - ");
              strcat(menulabel,_("All meshes"));
              glutAddMenuEntry(menulabel,-useitem-10);
            }
          }
          glutAddMenuEntry("-",MENU_DUMMY3);
          for(ii=0;ii<npatchinfo;ii++){
            patchdata *patch1, *patch2;

            i = patchorderindex[ii];
            patch2 = patchinfo + i;
            if(ii==0){
              nloadpatchsubmenus=0;
              strcpy(menulabel,patch2->label.longlabel);
              strcat(menulabel," - ");
              strcat(menulabel,_("Single mesh"));
              glutAddSubMenu(menulabel,loadpatchsubmenus[nloadpatchsubmenus]);
              nloadpatchsubmenus++;
            }
            else{
              patch1 = patchinfo + patchorderindex[ii-1];
              if(strcmp(patch1->label.longlabel,patch2->label.longlabel)!=0){
                strcpy(menulabel,patch2->label.longlabel);
                strcat(menulabel," - Single mesh");
                glutAddSubMenu(menulabel,loadpatchsubmenus[nloadpatchsubmenus]);
                nloadpatchsubmenus++;
              }
            }
          }
        }
      }
      glutAddMenuEntry("-",MENU_DUMMY3);
      glutAddMenuEntry(_("Update bounds"),MENU_UPDATEBOUNDS);
      if(npatchloaded>1){
        glutAddSubMenu(_("Unload"),unloadpatchmenu);
      }
      else{
       glutAddMenuEntry(_("Unload"),UNLOAD_ALL);
      }
    }

/* --------------------------------load iso menu -------------------------- */

    if(nisoinfo>0){
      int ii;

      CREATEMENU(unloadisomenu,UnloadIsoMenu);
      for(ii=0;ii<nisoinfo;ii++){
        isodata *isoi;
        char menulabel[1024];

        i = isoorderindex[ii];
        isoi = isoinfo + i;
        if(isoi->loaded==0)continue;
        STRCPY(menulabel,isoi->menulabel);
        glutAddMenuEntry(menulabel,i);
      }
      glutAddMenuEntry("Unload all",UNLOAD_ALL);

      if(nisoinfo>0){
        if(isosubmenus==NULL){
          NewMemory((void **)&isosubmenus,nisoinfo*sizeof(int));
        }
        nisosubmenus=0;

        CREATEMENU(isosubmenus[nisosubmenus],LoadIsoMenu);
        nisosubmenus++;
      }

      if(nmeshes==1){
        CREATEMENU(loadisomenu,LoadIsoMenu);
      }
      for(ii=0;ii<nisoinfo;ii++){
        isodata *iso1, *iso2;
        char menulabel[1024];

        i = isoorderindex[ii];
        if(ii>0){
          iso1 = isoinfo + isoorderindex[ii-1];
          iso2 = isoinfo + isoorderindex[ii];
          if(nmeshes>1&&strcmp(iso1->surface_label.longlabel,iso2->surface_label.longlabel)!=0){
            CREATEMENU(isosubmenus[nisosubmenus],LoadIsoMenu);
            nisosubmenus++;
          }
        }
        if(isoinfo[i].loaded==1){
          STRCPY(menulabel,"*");
          STRCAT(menulabel,isoinfo[i].menulabel);
        }
        else{
          STRCPY(menulabel,isoinfo[i].menulabel);
        }
        glutAddMenuEntry(menulabel,i);
      }

      {
        int useitem;
        isodata *isoi, *isoj;

        if(nmeshes>1){
         CREATEMENU(loadisomenu,LoadIsoMenu);
          for(i=0;i<nisoinfo;i++){
            int j;

            useitem=i;
            isoi = isoinfo + i;
            for(j=0;j<i;j++){
              isoj = isoinfo + j;
              if(strcmp(isoi->surface_label.longlabel,isoj->surface_label.longlabel)==0){
                useitem=-1;
                break;
              }
            }
            if(useitem!=-1){
              char menulabel[1024];

              strcpy(menulabel,isoi->surface_label.longlabel);
              strcat(menulabel," - ");
              strcat(menulabel,_("All meshes"));
              glutAddMenuEntry(menulabel,-useitem-10);
            }
          }
          glutAddMenuEntry("-",MENU_DUMMY3);

          for(ii=0;ii<nisoinfo;ii++){
            isodata *iso1, *iso2;
            char menulabel[1024];

            i = isoorderindex[ii];
            iso1 = isoinfo + i;
            if(ii==0){
              nisosubmenus=0;
              strcpy(menulabel,iso1->surface_label.longlabel);
              strcat(menulabel," - Single mesh");
              glutAddSubMenu(menulabel,isosubmenus[nisosubmenus]);
              nisosubmenus++;
            }
            else{
              iso2 = isoinfo + isoorderindex[ii-1];
              if(strcmp(iso1->surface_label.longlabel,iso2->surface_label.longlabel)!=0){
                strcpy(menulabel,iso1->surface_label.longlabel);
                strcat(menulabel," - Single mesh");
                glutAddSubMenu(menulabel,isosubmenus[nisosubmenus]);
                nisosubmenus++;
              }
            }
          }

       }
     }

      if(nisoloaded>1){
        glutAddSubMenu(_("Unload"),unloadisomenu);
      }
      else{
       glutAddMenuEntry(_("Unload"),UNLOAD_ALL);
      }
    }

/* --------------------------------zone menu -------------------------- */

    if(nzoneinfo>0){
      CREATEMENU(zonemenu,ZoneMenu);
      for(i=0;i<nzoneinfo;i++){
        zonedata *zonei;
        char menulabel[1024];
        int n;

        zonei = zoneinfo + i;
        if(zonefilenum==i){
          STRCPY(menulabel,"*");
          STRCAT(menulabel,zonei->file);
        }
        else{STRCPY(menulabel,zonei->file);}
        STRCAT(menulabel,", ");
        for(n=0;n<3;n++){
          STRCAT(menulabel,zonei->label[n].shortlabel);
          STRCAT(menulabel,", ");
        }
        STRCAT(menulabel,zonei->label[3].shortlabel);
        glutAddMenuEntry(menulabel,i);
      }
      glutAddMenuEntry("Unload",UNLOAD_ALL);

    }
/* -------------------------------- compress menu -------------------------- */

#ifdef pp_COMPRESS
    if(smokezippath != NULL && (npatchinfo > 0 || nsmoke3dinfo > 0 || nsliceinfo > 0)){
    CREATEMENU(compressmenu,CompressMenu);
    glutAddMenuEntry(_("Compression options"),MENU_DUMMY);  // -c
    if(overwrite_all==1){
      glutAddMenuEntry(_("  *Overwrite compressed files"),MENU_OVERWRITECOMPRESS);  // -f
    }
    else{
      glutAddMenuEntry(_("  Overwrite compressed files"),MENU_OVERWRITECOMPRESS);  // -f
    }
    if(compress_autoloaded==1){
      glutAddMenuEntry(_("  *Compress only autoloaded files"),MENU_COMPRESSAUTOLOAD);  // -f
    }
    else{
      glutAddMenuEntry(_("  Compress only autoloaded files"),MENU_COMPRESSAUTOLOAD);  // -f
    }
    glutAddMenuEntry("-",MENU_DUMMY);  // -c
    glutAddMenuEntry(_("Compress now"),MENU_COMPRESSNOW);
    glutAddMenuEntry(_("Erase compressed files"),MENU_ERASECOMPRESS);  // -c
  }
#endif


/* --------------------------------inisub menu -------------------------- */
  {
    int n_inifiles;
    inifiledata *inifile;

    n_inifiles=0;
    for(inifile=first_inifile.next;inifile->next!=NULL;inifile=inifile->next){
      if(inifile->file!=NULL&&file_exists(inifile->file)==1){
        n_inifiles++;
      }
    }
    if(n_inifiles>0){
      CREATEMENU(inisubmenu,IniSubMenu);
      if(caseini_filename!=NULL&&file_exists(caseini_filename)==1){
        glutAddMenuEntry(caseini_filename,MENU_READCASEINI);
      }
      for(inifile=first_inifile.next;inifile->next!=NULL;inifile=inifile->next){
        if(inifile->file!=NULL&&file_exists(inifile->file)==1){
          glutAddMenuEntry(inifile->file,inifile->id);
        }
      }
    }
  }

/* --------------------------------smokeviewini menu -------------------------- */

    CREATEMENU(smokeviewinimenu,SmokeviewIniMenu);
   {
    inifiledata *inifile;
    int n_inifiles;

    n_inifiles=0;
    for(inifile=first_inifile.next;inifile->next!=NULL;inifile=inifile->next){
      if(inifile->file!=NULL&&file_exists(inifile->file)==1){
        n_inifiles++;
      }
    }
    if( n_inifiles>0||file_exists(INIfile)==1||file_exists(caseini_filename)==1||file_exists(smokeviewini)==1){
      if(n_inifiles==0){
        glutAddMenuEntry(_("Read ini files"),MENU_READINI);
      }
      else{
        glutAddSubMenu(_("Read ini files"),inisubmenu);
      }
    }
  }

    glutAddMenuEntry(_(WRITEINIfile),MENU_WRITEINI);

    {
      char caselabel[255];

      STRCPY(caselabel,_("Write"));
      STRCAT(caselabel," ");
      STRCAT(caselabel,caseini_filename);

      glutAddMenuEntry(caselabel,MENU_WRITECASEINI);
    }

    if(ndeviceinfo>0){
      glutAddMenuEntry("-",MENU_DUMMY);
      glutAddMenuEntry(_("Read .svo files"),MENU_READSVO);
    }

    CREATEMENU(reloadmenu,ReloadMenu);
    glutAddMenuEntry(_("Reload Now"),RELOAD_NOW);
    if(periodic_value==1)glutAddMenuEntry(_("*Reload every 1 minute"),1);
    if(periodic_value!=1)glutAddMenuEntry(_("Reload every 1 minute"),1);
    if(periodic_value==5)glutAddMenuEntry(_("*Reload every 5 minutes"),5);
    if(periodic_value!=5)glutAddMenuEntry(_("Reload every 5 minutes"),5);
    if(periodic_value==10)glutAddMenuEntry(_("*Reload every 10 minutes"),10);
    if(periodic_value!=10)glutAddMenuEntry(_("Reload every 10 minutes"),10);
    glutAddMenuEntry(_("Stop Rendering"),STOP_RENDERING);


    {
      int nscripts;

      nscripts=0;
      if(script_recording==NULL){
        scriptfiledata *scriptfile;
        STRUCTSTAT statbuffer;

        for(scriptfile=first_scriptfile.next;scriptfile->next!=NULL;scriptfile=scriptfile->next){
          char *file;
          int len;

          file=scriptfile->file;
          if(file==NULL)continue;
          len = strlen(file);
          if(len<=0)continue;
          if(STAT(file,&statbuffer)!=0)continue;

          nscripts++;
        }

        if(nscripts>0){
          CREATEMENU(scriptlistmenu,ScriptMenu);
          for(scriptfile=first_scriptfile.next;scriptfile->next!=NULL;scriptfile=scriptfile->next){
            char *file;
            int len;
            char menulabel[1024];

            file=scriptfile->file;
            if(file==NULL)continue;
            len = strlen(file);
            if(len<=0)continue;
            if(STAT(file,&statbuffer)!=0)continue;

            strcpy(menulabel,"  ");
            strcat(menulabel,file);
            glutAddMenuEntry(menulabel,scriptfile->id);
          }
          CREATEMENU(scriptsteplistmenu,ScriptMenu2);
          for(scriptfile=first_scriptfile.next;scriptfile->next!=NULL;scriptfile=scriptfile->next){
            char *file;
            int len;
            char menulabel[1024];

            file=scriptfile->file;
            if(file==NULL)continue;
            len = strlen(file);
            if(len<=0)continue;
            if(STAT(file,&statbuffer)!=0)continue;

            strcpy(menulabel,"  ");
            strcat(menulabel,file);
            glutAddMenuEntry(menulabel,scriptfile->id);
          }
        }
      }
#ifdef pp_LUA
    {
      int nluascripts;

      nluascripts=0;
      // For Lua, the list of scripts is simply a list of filenames in the
      // directory with the right extension.
      luascriptfiledata *luascriptfile;
      STRUCTSTAT luastatbuffer;

      for(luascriptfile=first_luascriptfile.next;luascriptfile->next!=NULL;luascriptfile=luascriptfile->next){
        char *file;
        int len;

        file=luascriptfile->file;
        if(file==NULL)continue;
        len = strlen(file);
        if(len<=0)continue;
        if(STAT(file,&luastatbuffer)!=0)continue;

        nluascripts++;
      }
      if(nluascripts>0){
        CREATEMENU(luascriptlistmenu,LuaScriptMenu);
        for(luascriptfile=first_luascriptfile.next;luascriptfile->next!=NULL;luascriptfile=luascriptfile->next){
          char *file;
          int len;
          char menulabel[1024];

          file=luascriptfile->file;
          if(file==NULL)continue;
          len = strlen(file);
          if(len<=0)continue;
          if(STAT(file,&luastatbuffer)!=0)continue;

          strcpy(menulabel,"  ");
          strcat(menulabel,file);
          glutAddMenuEntry(menulabel,luascriptfile->id);
        }
      }
      CREATEMENU(luascriptmenu,LuaScriptMenu);
      if(nluascripts>0){
        glutAddSubMenu(_("Run"),luascriptlistmenu);
      }
    }
#endif

      CREATEMENU(scriptrecordmenu,ScriptMenu);
      if(script_recording==NULL){
        glutAddMenuEntry(_("Start"),SCRIPT_START_RECORDING);
        glutAddMenuEntry(_("Start (disable file loading)"),SCRIPT_START_RECORDING2);
      }
      glutAddMenuEntry(_("Stop"),SCRIPT_STOP_RECORDING);

      CREATEMENU(scriptmenu,ScriptMenu);
      if(nscripts>0){
        glutAddSubMenu(_("Run"),scriptlistmenu);
        glutAddSubMenu(_("Step (using ^)"),scriptsteplistmenu);
        if(script_step==1)glutAddMenuEntry(_("Continue"),SCRIPT_CONTINUE);
        if(script_step==1||current_script_command!=NULL)glutAddMenuEntry(_("Cancel"),SCRIPT_CANCEL);
      }
      glutAddSubMenu(_("Record"),scriptrecordmenu);
    }

  /* --------------------------------loadunload menu -------------------------- */
    {
      char loadmenulabel[100];
      char steplabel[100];

      CREATEMENU(loadunloadmenu,LoadUnloadMenu);
      strcpy(steplabel,_("error: steplabel not defined"));
      if(nsmoke3dinfo>0){
        strcpy(loadmenulabel,_("3D smoke"));
        if(smoke3dframeskip>0){
          sprintf(steplabel,"/Skip %i",smoke3dframeskip);
          strcat(loadmenulabel,steplabel);
        }
        glutAddSubMenu(loadmenulabel,loadsmoke3dmenu);
      }
      if(nvolrenderinfo>0&&smokediff==0){
        char vlabel[256];

        strcpy(vlabel,_("3D smoke (volume rendered)"));
        glutAddSubMenu(vlabel,loadvolsmoke3dmenu);
      }
      if(manual_terrain==1&&nterraininfo>0){
        glutAddSubMenu(_("Terrain"),loadterrainmenu);
      }
      if(nsliceinfo>0&&nmultisliceinfo<nsliceinfo){
        strcpy(loadmenulabel,_("Multi-Slices"));
        if(sliceframeskip>0){
          sprintf(steplabel,"/Skip %i",sliceframeskip);
          strcat(loadmenulabel,steplabel);
        }
        glutAddSubMenu(loadmenulabel,loadmultislicemenu);
      }
      if(nvsliceinfo>0&&nmultivsliceinfo<nvsliceinfo){
        strcpy(loadmenulabel,_("Multi-Vector Slices"));
        if(sliceframeskip>0){
          sprintf(steplabel,"/Skip %i",sliceframeskip);
          strcat(loadmenulabel,steplabel);
        }
        glutAddSubMenu(loadmenulabel,loadmultivslicemenu);
      }
      if(nsliceinfo>0){
        strcpy(loadmenulabel,"Slice file");
        if(sliceframeskip>0){
          sprintf(steplabel,"/Skip %i",sliceframeskip);
          strcat(loadmenulabel,steplabel);
        }
        glutAddSubMenu(loadmenulabel,loadslicemenu);
      }
      if(nvsliceinfo>0){
        strcpy(loadmenulabel,_("Vector slices"));
        if(sliceframestep>1){
          sprintf(steplabel,"/Skip %i",sliceframeskip);
          strcat(loadmenulabel,steplabel);
        }
        glutAddSubMenu(loadmenulabel,vslicemenu);
      }
      if(nisoinfo>0){
        strcpy(loadmenulabel,"Isosurface file");
        if(isoframeskip_global>0){
          sprintf(steplabel,"/Skip %i",isoframeskip_global);
          strcat(loadmenulabel,steplabel);
        }
        glutAddSubMenu(loadmenulabel,loadisomenu);
      }
      if(npatchinfo>0){
        strcpy(loadmenulabel,"Boundary file");
        if(boundframeskip>0){
          sprintf(steplabel,"/Skip %i",boundframeskip);
          strcat(loadmenulabel,steplabel);
        }
        glutAddSubMenu(loadmenulabel,loadpatchmenu);
      }
      if(npartinfo>0){
        if(nevac!=npartinfo){
          strcpy(loadmenulabel,"Particle file");
          if(partframeskip>0){
            sprintf(steplabel,"/Skip Frame %i",partframeskip);
            strcat(loadmenulabel,steplabel);
          }
          glutAddSubMenu(loadmenulabel,particlemenu);
        }
        if(nevac>0){
          strcpy(loadmenulabel,_("Evacuation"));
          if(partframeskip>0){
            sprintf(steplabel,"/Skip Frame %i",partframeskip);
            strcat(loadmenulabel,steplabel);
          }
          glutAddSubMenu(loadmenulabel,evacmenu);
        }
      }
      if(nplot3dinfo>0)glutAddSubMenu("Plot3d file",loadplot3dmenu);
      if(nzoneinfo>0){
        strcpy(loadmenulabel,"Zone fire file");
        glutAddSubMenu(loadmenulabel,zonemenu);
      }
      if(glui_active==1){
        glutAddMenuEntry("-",MENU_DUMMY);
      }
      glutAddSubMenu(_("Configuration files"),smokeviewinimenu);
      glutAddSubMenu(_("Scripts"),scriptmenu);
#ifdef pp_LUA
      glutAddSubMenu(_("Lua Scripts"),luascriptmenu);
#endif
#ifdef pp_COMPRESS
      if(smokezippath!=NULL&&(npatchinfo>0||nsmoke3dinfo>0||nsliceinfo>0)){
        glutAddSubMenu(_("Compression"),compressmenu);
      }
#endif
      if(showfiles==1)glutAddMenuEntry(_("*Show file names"),SHOWFILES);
      if(showfiles==0)glutAddMenuEntry(_("Show file names"),SHOWFILES);

      {
        char menulabel[1024];

        strcpy(menulabel,"");
        if(redirect==1)strcat(menulabel,"*");
        strcat(menulabel,"Redirect messages to ");
        strcat(menulabel,log_filename);
        glutAddMenuEntry(menulabel,REDIRECT);
      }

      glutAddSubMenu(_("Reload"),reloadmenu);
      glutAddMenuEntry(_("Unload all"),UNLOADALL);
    }

/* --------------------------------main menu -------------------------- */
    if(trainer_mode==1){
      CREATEMENU(trainerviewmenu,TrainerViewMenu);
      if(AnySmoke(NULL)==1){
        if(trainerload==1)glutAddMenuEntry(_("*Realistic"),MENU_TRAINER_smoke);
        if(trainerload!=1)glutAddMenuEntry(_("Realistic"),MENU_TRAINER_smoke);
      }
      if(AnySlices("TEMPERATURE")==1){
        if(trainerload==2)glutAddMenuEntry(_("*Temperature"),MENU_TRAINER_temp);
        if(trainerload!=2)glutAddMenuEntry(_("Temperature"),MENU_TRAINER_temp);
      }
      if(AnySlices("oxygen")==1||
         AnySlices("oxygen VOLUME FRACTION")==1){
        if(trainerload==3)glutAddMenuEntry(_("*Oxygen"),MENU_TRAINER_oxy);
        if(trainerload!=3)glutAddMenuEntry(_("Oxygen"),MENU_TRAINER_oxy);
      }
      glutAddMenuEntry(_("Clear"),MENU_TRAINER_CLEAR);
    }

    CREATEMENU(mainmenu,MainMenu);
    if(trainer_mode==0){
      glutAddSubMenu(_("Load/Unload"),loadunloadmenu);
      glutAddSubMenu(_("Show/Hide"),showhidemenu);
      glutAddSubMenu(_("Options"),optionmenu);
      glutAddSubMenu(_("Dialogs"),dialogmenu);
      glutAddSubMenu(_("Help"),helpmenu);
      glutAddMenuEntry(_("Quit"),MENU_MAIN_QUIT);
    }
    if(trainer_active==1){
      if(trainer_mode==1){
        glutAddMenuEntry(_("Smokeview menus"),MENU_MAIN_TRAINERTOGGLE);
      }
      else{
        glutAddMenuEntry(_("Trainer menus"),MENU_MAIN_TRAINERTOGGLE);
      }
    }
    updatemenu=0;
#ifdef _DEBUG
  in_menu=0;
  PRINTF("nmenus=%i\n",nmenus);
#endif

}

/* ------------------ MenuStatus ------------------------ */

void MenuStatus_CB(int status, int x, int y){
  float *eye_xyz;
  menustatus=status;
  /* keep scene from "bouncing" around when leaving a menu */
  start_xyz0[0]=x;
  start_xyz0[1]=y;
  /*touring=0;*/
  mouse_down_xy0[0]=x; mouse_down_xy0[1]=y;
  eye_xyz = camera_current->eye;
  eye_xyz0[0]=eye_xyz[0];
  eye_xyz0[1]=eye_xyz[1];
  eye_xyz0[2]=eye_xyz[2];
}
