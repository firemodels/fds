// $Date$ 
// $Revision$
// $Author$

#define MENU_TIMEVIEW -103
#define MENU_SAVEVIEW -101
#define MENU_STARTUPVIEW -102
#define MENU_OUTLINEVIEW -104
#define MENU_SIZEPRESERVING -105
#define MENU_DUMMY -999

// svn revision character string
char menu_revision[]="$Revision$";
void ScriptMenu(int var);
void add_scriptlist(char *file, int id);
void update_glui_render();
void PropMenu(int value);

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
      printf( "Unable to locate the directory: %s\n", smv_directory );
    }
    else{
      *openfile=1;
    }
  }
}
#endif

/* ------------------ TrainerViewMenu ------------------------ */

void TrainerViewMenu(int value){
  switch (value) {
  case 1:   // realistic
    HideAllSlices();
    ShowAllSmoke();
    trainerload=1;
    break;
  case 2:  // temperature slices
//    LoadSmoke3DMenu(-1);
    HideAllSmoke();
    ShowAllSlices("TEMPERATURE",NULL);
    trainerload=2;
    break;
  case 3:  //  oxygen slices
//    LoadSmoke3DMenu(-1);
    HideAllSmoke();
    ShowAllSlices("OXYGEN","OXYGEN VOLUME FRACTION");
    trainerload=3;
    break;
  case 998: // unload
//    LoadSmoke3DMenu(-1);
//    UnloadSliceMenu(-1);
    LoadUnloadMenu(UNLOADALL);
    trainerload=0;
    break;
  default:
    ASSERT(FFALSE);
  }
  updatemenu=1;
  GLUTPOSTREDISPLAY
}

/* ------------------ MainMenu ------------------------ */

void MainMenu(int value){

  if(value==3){
    if(scriptoutstream!=NULL){
      ScriptMenu(STOP_RECORDING_SCRIPT);
    }
    exit(0);
  }
  if(value==1){
    defaulttour();
  }
  if(value==997){
    trainer_mode=1-trainer_mode;
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY  
}

/* ------------------ sv_MainMenu ------------------------ */

void svWINAPI sv_MainMenu(int value){
  MainMenu(value);
}

/* ------------------ StaticVariableMenu ------------------------ */

void StaticVariableMenu(int value){
  mesh *meshi;
  meshi=current_mesh;
  plotn=value;
  plotstate=STATIC_PLOTS;
  visGrid=0;
  if(visiso==1){
    updateshowstep(1,ISO);
  }
  updatesurface();
  if(meshi->visx==1){
    updateshowstep(1,DIRX);
  }
  if(meshi->visy==1){
    updateshowstep(1,DIRY);
  }
  if(meshi->visz==1){
    updateshowstep(1,DIRZ);
  }
  if(meshi->visx==0&&meshi->visy==0&&meshi->visz==0){
    updateshowstep(1,DIRY);
  }
  updateallplotslices();
  updatemenu=1;  
  GLUTPOSTREDISPLAY  
  updateplot3dlistindex();
}

/* ------------------ IsoVariableMenu ------------------------ */

void IsoVariableMenu(int value){
  mesh *meshi;
  meshi=current_mesh;
  if(ReadPlot3dFile==1){
    plotn=value;
    if(meshi->visx==1){updateshowstep(1,DIRX);}
    if(meshi->visy==1){updateshowstep(1,DIRY);}
    if(meshi->visz==1){updateshowstep(1,DIRZ);}
    updateshowstep(1,ISO);
    updatesurface();
    plotstate=STATIC_PLOTS;
    updateplotslice(1);
    updateplotslice(2);
    updateplotslice(3);
    updatemenu=1;  
    GLUTPOSTREDISPLAY
    updateplot3dlistindex();
  }
}

/* ------------------ LabelMenu ------------------------ */

void LabelMenu(int value){
  updatemenu=1;  
  GLUTPOSTREDISPLAY
  ASSERTFLAG(visColorLabels);
  ASSERTFLAG(visTimeLabels);
  ASSERTFLAG(visTitle0);
  ASSERTFLAG(visFramerate);
  ASSERTFLAG(visaxislabels);
  ASSERTFLAG(vis_slice_average);
  switch (value){
   case 0:
    visColorLabels=1-visColorLabels;
    break;
   case 1:
    visTimeLabels=1-visTimeLabels;
    break;
   case 2:
    visTitle0=1-visTitle0;
    break;
   case 3:
    visFramerate = 1 - visFramerate;
    break;
   case 4:
    vis_user_ticks=1;
    visColorLabels=1;
    visTimeLabels=1;
    visTitle0=1;
    visTitle1=1;
    visTitle2=1;
    visFramerate=1;
#ifdef pp_memstatus
    visAvailmemory=1;
#endif
    visaxislabels=1;
    visTimelabel=1;
    visFramelabel=1;
    visLabels=1;
    visBlocklabel=1;
    vis_slice_average=1;
    if(nticks>0)visTicks=1;
    visgridloc=1;
    visHRRlabel=1;
    show_hrrcutoff=1;
    visFramelabel=1;
    if(hrrinfo!=NULL){
      if(hrrinfo->display==0){
        hrrinfo->display=1;
        updatetimes();
      }
      hrrinfo->display=1;
    }
    break;
   case 5:
    vis_user_ticks=0;
    visColorLabels=0;
    visTimeLabels=0;
    visTitle0=0;
    visTitle1=0;
    visTitle2=0;
    visFramerate=0;
    visaxislabels=0;
    visLabels=0;
    visTimelabel=0;
    visFramelabel=0;
    visBlocklabel=0;
    visHRRlabel=0;
    show_hrrcutoff=0;
    if(hrrinfo!=NULL){
      if(hrrinfo->display==1){
        hrrinfo->display=0;
        updatetimes();
      }
      hrrinfo->display=0;
    }
    if(nticks>0)visTicks=0;
    visgridloc=0;
    vis_slice_average=0;
#ifdef pp_memstatus
    visAvailmemory=0;
#endif
    break;
   case 6:
    visaxislabels = 1 - visaxislabels;
    break;
   case 7:
     visLabels = 1 - visLabels;
     break;
   case 8:
     visTimelabel=1-visTimelabel;
     if(visTimelabel==1)visTimeLabels=1;
     break;
   case 9:
     visFramelabel=1-visFramelabel;
     if(visFramelabel==1)visTimeLabels=1;
     if(visFramelabel==1){
       visHRRlabel=0;
       if(hrrinfo!=NULL){
         hrrinfo->display=visHRRlabel;
         updatetimes();
       }
     }
     break;
   case 10:
     visBlocklabel=1-visBlocklabel;
     break;
#ifdef pp_memstatus
   case 11:
     visAvailmemory = 1 - visAvailmemory;
     break;
#endif
   case 12:
     visTicks=1-visTicks;
     break;
   case 13:
     vishmsTimelabel = 1 - vishmsTimelabel;
     break;
   case 14:
     visgridloc = 1 - visgridloc;
     break;
   case 15:
     vis_slice_average = 1 - vis_slice_average;
     break;
   case 16:
     visHRRlabel=1-visHRRlabel;
     if(visHRRlabel==1)visTimeLabels=1;
     if(hrrinfo!=NULL){
       hrrinfo->display=visHRRlabel;
       updatetimes();
     }
     break;
   case 17:
     show_hrrcutoff=1-show_hrrcutoff;
     break;
   case 18:
     vis_user_ticks = 1 - vis_user_ticks;
     break;
   default:
     ASSERT(FFALSE);
     break;
  }
  updateshowtitles();
  set_labels_controls();

}

/* ------------------ updateshowtitles ------------------------ */

void updateshowtitles(void){
  ntitles=0;
  if(visTitle0==1)ntitles++;
  if(strlen(TITLE1)!=0&&visTitle1==1){
    ntitles++;
    showtitle1=1;
  }
  if(strlen(TITLE2)!=0&&visTitle2==1){
    ntitles++;
    showtitle2=1;
  }
  visTitle=0;
  if(visTitle0==1||showtitle1==1||showtitle2==1)visTitle=1;
}

/* ------------------ LightingMenu ------------------------ */

void LightingMenu(int value){
    if(value<0)return;
    ASSERTFLAG(visLIGHT0);
    ASSERTFLAG(visLIGHT1);
    switch (value){
      case 1: 
        visLIGHT0 = 1 - visLIGHT0; 
        break;
      case 2: 
        visLIGHT1 = 1 - visLIGHT1; 
        break;
      case 3: 
        visLIGHT1 = 1 - visLIGHT1; 
        visLIGHT0 = 1 - visLIGHT0; 
        break;
      case 4: 
        visLIGHT0=1; 
        visLIGHT1=1; 
        break;
      case 5: 
        visLIGHT0=0; 
        visLIGHT1=0; 
        break;
#ifdef pp_LIGHT
      case 6:
        show_smokelighting=1-show_smokelighting;
        update_showlight();
        break;
#endif
      default:
        ASSERT(FFALSE);
        break;
    }
    UpdateLIGHTS=1;
    updatemenu=1;  
    GLUTPOSTREDISPLAY
}

/* ------------------ ColorBarMenu ------------------------ */

void ColorBarMenu(int value){
  if(value==-999)return;
  updatemenu=1;
  GLUTPOSTREDISPLAY
  if(value<0){
    switch (value){
    case -2:
      colorbarflip=1-colorbarflip;
      break;
    case -3:
      colorbarcycle++;
      if(colorbarcycle>=nrgb)colorbarcycle=0;
      break;
    case -4:
      colorbarcycle=0;
      setbw=0;
      break;
    case -5:
      viscolorbarpath=1-viscolorbarpath;
      break;
    case -6:
      if(colorband==1){
        colorband=5;
      }
      else{
        colorband=1;
      }
      break;
    case -7:
      show_extremedata=1-show_extremedata;
      update_extreme();
      updatecolors(-1);
      break;
    case -11:
     background_flip = 1-background_flip;
     updatecolors(-1);
     set_labels_controls();
     break;
    case -12:
     setbw=1-setbw;
     if(setbw==1){
       colorbartype_save=colorbartype;
       ColorBarMenu(bw_colorbar_index);
     }
     else{
       colorbartype=colorbartype_save;
       ColorBarMenu(colorbartype);
     }
     updatecolors(-1);
     set_labels_controls();
     break;
   case -13:
     transparentflag=1-transparentflag;
     updatecolors(-1);
     set_labels_controls();
     break;
   case -14:
     colorbarflip=1-colorbarflip;updatecolors(-1);
     break;
   case -15:
     colorbarcycle++;
     if(colorbarcycle>=nrgb)colorbarcycle=0;
     updatecolors(-1);
     break;
   case -16:
     colorbarcycle=0;
     background_flip=0;
     setbw=0;
     updatecolors(-1);
     break;
   default:
     ASSERT(FFALSE);
     break;
   }
  }
  if(value>=0){
    colorbartype=value;
    current_colorbar=colorbarinfo+colorbartype;
    update_colorbar_type();
    if(value==bw_colorbar_index){
      setbw=1;
    }
    else{
      setbw=0;
    }
  }
  if(value>-10){
    updatecolors(-1);
  }
}

/* ------------------ ShadeMenu ------------------------ */

void ShadeMenu(int value){
  ColorBarMenu(-10-value);
}
/* ------------------ Smoke3DShowMenu ------------------------ */

void Smoke3DShowMenu(int value){
  smoke3d *smoke3di;
  int i;

  updatemenu=1;  
  GLUTPOSTREDISPLAY
  if(value<0){
    switch (value){
    case SHOW_ALL:
      plotstate=DYNAMIC_PLOTS;
      for(i=0;i<nsmoke3d_files;i++){
        smoke3di = smoke3dinfo + i;
        if(smoke3di->loaded==1)smoke3di->display=1;
      }
      break;
    case HIDE_ALL:
      for(i=0;i<nsmoke3d_files;i++){
        smoke3di = smoke3dinfo + i;
        if(smoke3di->loaded==1)smoke3di->display=0;
      }
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
  iso *isoi;

  nisolevels=loaded_isomesh->nisolevels;
  showlevels=loaded_isomesh->showlevels;

  ASSERTFLAG(isonormtype);
  ASSERTFLAG(showisonormals);
  switch (value){
   case  4:
    isonormtype = 1 - isonormtype;
    break;
   case 5:
    showisonormals = 1 - showisonormals;
    break;
   case 1:
   case 2:
   case 3:
    visAIso=value;
    if(visAIso!=0){
      plotstate=DYNAMIC_PLOTS;
    }
    update_glui_isotype();
    break;
   case 94:
    transparent_state=ALL_SOLID;
    break;
   case 95:
    transparent_state=ALL_TRANSPARENT;
    break;
   case 96:
    transparent_state=MIN_SOLID;
    break;
   case 97:
    transparent_state=MAX_SOLID;
    break;
   case 98:
    for(i=0;i<nisolevels;i++){
      showlevels[i]=0;
    }
    break;
   case 99:
    for(i=0;i<nisolevels;i++){
      showlevels[i]=1;
    }
    break;
   default:
    if(value>99&&value<999&&value-100<nisolevels){
     ASSERTFLAG(showlevels[value-100]);
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
     updateShow();
    }
    else if(value>=10001){
      if(value==10001){
        plotstate=DYNAMIC_PLOTS;
        for(i=0;i<niso_files;i++){
          isoinfo[i].display=1;
        }
      }
      else if(value==10002){
        for(i=0;i<niso_files;i++){
          isoinfo[i].display=0;
        }
      }
     updateShow();
    }
  }
  
  update_iso_showlevels();

  updatemenu=1;  
  GLUTPOSTREDISPLAY
}

/* ------------------ ShowVSliceMenu ------------------------ */

void ShowVSliceMenu(int value){
  vslice *vd;
  slice *sd;
  int i;
  updatemenu=1;  
  GLUTPOSTREDISPLAY
  if(value==SHOW_ALL){
    for(i=0;i<nvslice;i++){
      vd = vsliceinfo + i;
      if(vd->loaded==0)continue;
      vd->display=1;
    }
    show_all_slices=1;
    updatetimes();
    return;
  }
  if(value==HIDE_ALL){
    for(i=0;i<nvslice;i++){
      vd = vsliceinfo + i;
      if(vd->loaded==0)continue;
      vd->display=0;
    }
    show_all_slices=0;
    updatetimes();
    return;
  }
  if(value==-11){
    show_slice_in_obst=1-show_slice_in_obst;
    return;
  }
  if(value==-12){
    offset_slice=1-offset_slice;
    return;
  }
  vd = vsliceinfo + value;
  ASSERTFLAG(vd->display);
  if(islicetype==sliceinfo[vd->ival].type){
    if(plotstate!=DYNAMIC_PLOTS){
      plotstate=DYNAMIC_PLOTS;
      vd->display=1;
    }
    else{
      vd->display = 1 - vd->display;
    }
    if(vd->iu!=-1){
      sd=sliceinfo+vd->iu;
      sd->vloaded=vd->display;
    }
    if(vd->iv!=-1){
      sd=sliceinfo+vd->iv;
      sd->vloaded=vd->display;
    }
    if(vd->iw!=-1){
      sd=sliceinfo+vd->iw;
      sd->vloaded=vd->display;
    }
    if(vd->ival!=-1){
      sd=sliceinfo+vd->ival;
      sd->vloaded=vd->display;
    }
  }
  else{
    islicetype = sliceinfo[vd->ival].type;
    sd=sliceinfo+vd->ival;
    vd->display=1;
  }
  plotstate=getplotstate(DYNAMIC_PLOTS);
  updateShow();
}

/* ------------------ ShowHideSliceMenu ------------------------ */

void ShowHideSliceMenu(int value){
  int i;

  updatemenu=1;  
  GLUTPOSTREDISPLAY
  if(value<0){
    switch (value){
    case SHOW_ALL:
      for(i=0;i<nslice_files;i++){
        sliceinfo[i].display=1;
      }
      show_all_slices=1;
      break;
    case HIDE_ALL:
      for(i=0;i<nslice_files;i++){
        sliceinfo[i].display=0;
      }
      show_all_slices=0;
      break;
    case -11:
      show_slice_in_obst=1-show_slice_in_obst;
      break;
    case -12:
      offset_slice=1-offset_slice;
      break;
    case -13:
      planar_terrain_slice=1-planar_terrain_slice;
      break;
    default:
      ASSERT(FFALSE);
    }
  }
  else{
    slice *sd;

    sd = sliceinfo + value;
    ASSERTFLAG(sd->display);
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
  plotstate=getplotstate(DYNAMIC_PLOTS);

  updateglui();
  updateslicelistindex(slicefilenum);
  updateShow();
}

/* ------------------ ShowMultiSliceMenu ------------------------ */

void ShowMultiSliceMenu(int value){
  multislice *mslicei;
  slice *sd;
  int mdisplay;
  int i;

  updatemenu=1;  
  GLUTPOSTREDISPLAY
  switch (value){
  case SHOW_ALL:
  case HIDE_ALL:
    ShowHideSliceMenu(value);
    return;
  case -11:
    show_slice_in_obst=1-show_slice_in_obst;
    break;
  case -12:
    offset_slice=1-offset_slice;
    break;
  default:
    mslicei = multisliceinfo + value;
    mdisplay=0;
    if(islicetype==mslicei->type){
      if(plotstate!=DYNAMIC_PLOTS){
        plotstate=DYNAMIC_PLOTS;
        mdisplay=1;
      }
      else{
        mdisplay = 1 - mslicei->display;
      }
    }
    else{
      plotstate=DYNAMIC_PLOTS;
      islicetype=mslicei->type;
      mdisplay=1;
    }
    for(i=0;i<mslicei->nslices;i++){
      sd = sliceinfo + mslicei->islices[i];
      if(sd->loaded==0)continue;
      sd->display=mdisplay;
    }
    break;
  }
  updateslicefilenum();
  plotstate=getplotstate(DYNAMIC_PLOTS);

  updateglui();
  updateslicelistindex(slicefilenum);
  updateShow();
}

/* ------------------ ShowHideMenu ------------------------ */

void ShowHideMenu(int value){
  updatemenu=1;  
  GLUTPOSTREDISPLAY
  switch (value){
  case 13:
    if(plotstate==DYNAMIC_PLOTS){
      visEvac=1-visEvac;
    }
    else{
      plotstate=DYNAMIC_PLOTS;
      visEvac=1;
    }
    updatetimes();
    break;
  case 1:
    if(plotstate==DYNAMIC_PLOTS){
      visSmoke=1-visSmoke;
    }
    else{
      plotstate=DYNAMIC_PLOTS;
      visSmoke=1;
    }
    updatetimes();
    break;
  case 2:
    visTarg=1-visTarg;
    break;
  case 9:
    visSensor=1-visSensor;
    break;
  case 14:
    visSensorNorm=1-visSensorNorm;
    break;
  case 10:
    visHeat=1-visHeat;
    break;
  case 11:
    visSprink=1-visSprink;
    break;
  case 12:
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
  if(value==999)return;
  updatemenu=1;  
  GLUTPOSTREDISPLAY
  switch (value){
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

  switch (value){
#ifdef pp_SHOOTER
  case 27:
    showshooterDLG=1-showshooterDLG;
    if(showshooterDLG==1)show_shooter();
    if(showshooterDLG==0)hide_shooter();
    break;
#endif
  case 25:
    showtrainer=1-showtrainer;
    if(showtrainer==1)show_trainer();
    if(showtrainer==0)hide_trainer();
    break;
  case 22:
    showlabels=1-showlabels;
    if(showlabels==1)show_glui_labels();
    if(showlabels==0)hide_glui_labels();
    break;
  case 14:
  case 20:
  case 24:
    showbounds=1-showbounds;
    showglui3dsmoke=1-showglui3dsmoke;
    if(showbounds==1)show_glui_bounds();
    if(showbounds==0)hide_glui_bounds();
    if(value==20&&showbounds==1){
      open_smokepanel();
    }
    if(value==24&&showbounds==1){
      open_smokezippanel();
    }
    break;
  case 15:
    showmotion=1-showmotion;
    if(showmotion==1)show_glui_motion();
    if(showmotion==0)hide_glui_motion();
    break;
  case 21:
   showgluitour=1-showgluitour;
   if(showgluitour==1)show_glui_tour();
   if(showgluitour==0)hide_glui_tour();
   break;
  case 18:
    showclip=1-showclip;
    if(showclip==1)show_glui_clip();
    if(showclip==0)hide_glui_clip();
    break;
  case 19:
    showgluistereo=1-showgluistereo;
    if(showgluistereo==1)show_glui_stereo();
    if(showgluistereo==0)hide_glui_stereo();
    break;
  case 26:
    showwui=1-showwui;
    if(showwui==1)show_glui_wui();
    if(showwui==0)hide_glui_wui();
    break;
  case 23:
    showcolorbar=1-showcolorbar;
    if(showcolorbar==1){
      viscolorbarpath=1;
      show_glui_colorbar();
    }
    if(showcolorbar==0){
      viscolorbarpath=0;
      hide_glui_colorbar();
    }
    break;
  case 16:
    showedit=1-showedit;
    if(showedit==1){
      if(fds_filein!=NULL&&updategetlabels==1){
        CheckMemoryOff;
        getlabels(fds_filein);
        CheckMemoryOn;
        updategetlabels=0;
      }
      visBlocksSave=visBlocks;
      show_glui_edit();
      visBlocks=visBLOCKNormal;
    }
    if(showedit==0){
      hide_glui_edit();
      visBlocks=visBlocksSave;
    }
    update_trainer_outline();

    break;
  case -2:
    showlabels=0;
    hide_glui_labels();
    showbounds=0;
    showglui3dsmoke=0;
    hide_glui_bounds();
    showmotion=0;
    hide_glui_motion();
    showgluitour=0;
    hide_glui_tour();
    showclip=0;
    hide_glui_clip();
    showwui=0;
    hide_glui_wui();
    showgluistereo=0;
    hide_glui_stereo();
    showcolorbar=0;
    hide_glui_colorbar();
    if(showedit==1)DialogMenu(16);
    showtrainer=0;
    hide_trainer();
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  updatemenu=1;
}

/* ------------------ ZoomMenu ------------------------ */

void ZoomMenu(int value){
  if(value==999)return;
  updatemenu=1;  
  if(opengldefined==1){
    GLUTPOSTREDISPLAY
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
      ResetView(RESTORE_EXTERIOR_VIEW);
      update_projection_type();
    }
  }
  else if(zoomindex==-2){
    camera_current->projection_type=projection_type;
    update_projection_type();
    if(projection_type==0){
      update_camera_ypos(camera_current);
      //ResetView(RESTORE_EXTERIOR_VIEW);
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
      ResetView(RESTORE_EXTERIOR_VIEW_ZOOM);
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
    GLUTPOSTREDISPLAY
  }
  apertureindex=value;
  if(apertureindex<0)apertureindex=0;
  if(apertureindex>4)apertureindex=4;
  aperture=apertures[apertureindex];
}

/* ------------------ FontMenu ------------------------ */

void FontMenu(int value){
  updatemenu=1;  
  if(opengldefined==1){
    GLUTPOSTREDISPLAY
  }
  switch (value){
  case SMALL_FONT:
    fontindex=SMALL_FONT;
    large_font=GLUT_BITMAP_HELVETICA_12;
    small_font=GLUT_BITMAP_HELVETICA_10;
    large_font_height=12;
    small_font_height=10;
    fontWoffset=0;
    fontHoffset=0;
    dwinH=dwinHbase;
    dwinH=1.2*dwinHbase;
    break;
  case LARGE_FONT:
    fontindex=LARGE_FONT;
    large_font=GLUT_BITMAP_HELVETICA_18;
    small_font=GLUT_BITMAP_HELVETICA_18;
    large_font_height=18;
    small_font_height=18;
    fontWoffset=0;
    fontHoffset=0;
    dwinH=1.2*dwinHbase;
    break;
  case LARGE_FONT_SAFE:
    fontindex=LARGE_FONT_SAFE;
    large_font=GLUT_BITMAP_HELVETICA_18;
    small_font=GLUT_BITMAP_HELVETICA_18;
    large_font_height=18;
    small_font_height=18;
    fontWoffset=50;
    fontHoffset=50;
    break;
  default:
    ASSERT(FFALSE);
  }
  set_labels_controls();
}

/* ------------------ UnitsMenu ------------------------ */

void UnitsMenu(int value){
  int unitclass, unittype;
  int i;

  unitclass = value/1000;
  unittype = value - unitclass*1000;
  unitclasses[unitclass].active=unittype;
  if(value==-1){
    for(i=0;i<nunitclasses;i++){
      unitclasses[i].active=0;
    }
  }
  else if(value==-2){
    vishmsTimelabel = 1 - vishmsTimelabel;
    set_labels_controls();

  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY
}

/* ------------------ OptionMenu ------------------------ */

void OptionMenu(int value){
  if(value==999)return;
  updatemenu=1;  
  GLUTPOSTREDISPLAY
  if(value==1){
    Labels_CB(17); // run the benchmark
  }
  if(value==2){
    trainer_mode=1;
    if(showtrainer==0){
      showtrainer=1;
      show_trainer();
    }
    FontMenu(1);
    updatemenu=1;
  }
}

/* ------------------ ResetMenu ------------------------ */

void ResetMenu(int value){
  char line[256];

  if(value==MENU_DUMMY)return;
  switch (value){
  case MENU_SIZEPRESERVING:
    projection_type = 1 - projection_type;
    TRANSLATE_CB(PROJECTION);
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
    updatetimes();
    break;
  case MENU_SAVEVIEW:
    menu_view_number++;
    sprintf(line,"view %i",menu_view_number);
    add_list_view(line);
    break;
  case MENU_STARTUPVIEW:
    if(selected_view==-999)ResetMenu(MENU_SAVEVIEW);
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
  GLUTPOSTREDISPLAY
}

/* ------------------ RenderState ------------------------ */

void RenderState(int onoff){
  if(onoff==1){
    saveW=screenWidth;
    saveH=screenHeight;
    RenderGif=1;
    if(renderW==0||renderH==0){
      ResizeWindow(screenWidth,screenHeight);
    }
    else{
      if(renderW>max_screenWidth){
        render_double=render_double_state;
        ResizeWindow(max_screenWidth,max_screenHeight);
      }
      else{
        ResizeWindow(renderW,renderH);
      }
    }
  }
  else{
    render_double=0;
    render_double_state=0;
    RenderGif=0;
    screenWidth=saveW;
    screenHeight=saveH;
    ResizeWindow(screenWidth,screenHeight);
  }
}

/* ------------------ RenderMenu ------------------------ */

void RenderMenu(int value){
  slice *sd;
  int i,n;
  mesh *meshi;

  updatemenu=1;
  if(value>=10000)return;
  if(opengldefined==1){
    GLUTPOSTREDISPLAY
  }
  switch (value){
  case Render320:
    render_option=value;
    render_double_state=0;
    renderW=320;
    renderH=240;
    break;
  case Render640:
    render_option=value;
    render_double_state=0;
    renderW=640;
    renderH=480;
    break;
  case Render2Window:
    render_option=value;
    renderW=2*screenWidth;
    renderH=2*screenHeight;
    render_double_state=2;
    break;
  case RenderWindow:
    render_option=value;
    render_double_state=0;
    renderW=0;
    renderH=0;
    break;
  case RenderFull:
    render_option=value;
    render_double_state=0;
    renderW=max_screenWidth;
    renderH=max_screenHeight;
    break;
  case Render2Full:
    render_option=value;
    renderW=2*max_screenWidth;
    renderH=2*max_screenHeight;
    render_double_state=2;
    break;
  case Render4Full:
    render_option=value;
    renderW=4*max_screenWidth;
    renderH=4*max_screenHeight;
    render_double_state=4;
    break;
  case RenderCancel:
    render_double=0;
    render_double_state=0;
    RenderState(0);
    break;
  case RenderOnce:
    if(render_double_state!=0){
      render_double=render_double_state;
      keyboard('R',0,0);
    }
    else{
      keyboard('r',0,0);
    }
     break;
  case RenderPNG:
     renderfiletype=0;
     updatemenu=1;  
     break;
#ifdef pp_JPEG
  case RenderJPEG:
     renderfiletype=1;
     updatemenu=1;  
     break;
#endif
#ifdef pp_GDGIF
  case RenderGIF:
     renderfiletype=2;
     updatemenu=1;  
     break;
#endif
  default:
    if(RenderTime==0&&touring==0)return;
    if(touring==1){
      rendertourcount=0;
      tourangle=0.0;
    }
    RenderState(1);
    itime=0;
    for(i=0;i<nslice_files;i++){
      sd=sliceinfo+i;
      sd->islice=0;
    }
    iframe=iframebeg;
    for(i=0;i<nmeshes;i++){
      meshi=meshinfo+i;
      meshi->ipatch=0;
    }
    UpdateTimeLabels();
    RenderSkip=value;
    FlowDir=1;
    for(n=0;n<ntimes;n++){
      render_frame[n]=0;
    }
    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"RENDERALL\n");
      fprintf(scriptoutstream," %i\n",RenderSkip);
      fprintf(scriptoutstream,"\n");
    }
    break;
  }
}

/* ------------------ EvacShowMenu ------------------------ */

void EvacShowMenu(int value){
  particle *parti;
  int i;

  if(nevac==0)return;
  if(value==999)return;
  ASSERTFLAG(visEvac);
  if(value<0){
    value = -value;
    value--;
    parti = partinfo + value;
    parti->display = 1 - parti->display;
    updatemenu=1;  
    GLUTPOSTREDISPLAY
    plotstate=getplotstate(DYNAMIC_PLOTS);
    return;
  }
  if(plotstate==DYNAMIC_PLOTS){
    switch (value){
    case 3:
      visEvac=1;
      for(i=0;i<npart_files;i++){
        parti = partinfo + i;
        if(parti->loaded==0||parti->evac==0)continue;
        parti->display=1;
      }
      break;
    case 4:
      visEvac=0;
      for(i=0;i<npart_files;i++){
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
    switch (value){
    case 3:
      visEvac=1;
      for(i=0;i<npart_files;i++){
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
  plotstate=getplotstate(DYNAMIC_PLOTS);
  GLUTPOSTREDISPLAY

}

/* ------------------ ParticleShowMenu ------------------------ */

void ParticleShowMenu(int value){
  particle *parti;
  int i;

  if(npart_files==0)return;
  if(value==999)return;
  ASSERTFLAG(visSmoke);
  ASSERTFLAG(visEvac);
  ASSERTFLAG(visSprinkPart);
  ASSERTFLAG(visStaticSmoke);
  ASSERTFLAG(visSprinkPart);
  if(value<0){
    value = -value;
    value--;
    parti = partinfo + value;
    parti->display = 1 - parti->display;
    updatemenu=1;  
    GLUTPOSTREDISPLAY
    plotstate=getplotstate(DYNAMIC_PLOTS);
    return;
  }
  if(plotstate==DYNAMIC_PLOTS){
    switch (value){
      case 1:
        if(visSmokePart==2){
          visSmokePart=0;
        }
        else{
          visSmokePart=2;
        }
        break;
      case 2: 
        visSprinkPart = 1 - visSprinkPart; 
        break;
      case 3: 
        visSprinkPart=1; 
        visSmokePart=2; 
        visStaticSmoke=1; 
        for(i=0;i<npart_files;i++){
          parti = partinfo + i;
          if(parti->loaded==0||parti->evac==1)continue;
          parti->display=1;
        }
        break;
      case 5: 
        visStaticSmoke = 1 - visStaticSmoke; 
        break;
      case 4: 
        visSprinkPart=0; 
        visSmokePart=0; 
        visStaticSmoke=0;
        for(i=0;i<npart_files;i++){
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
    for(i=0;i<npart_files;i++){
      parti = partinfo + i;
      if(parti->loaded==0||parti->evac==1)continue;
      parti->display_droplet=0;
      parti->display_smoke=0;
      if(visSmokePart!=0)parti->display_smoke=1;
      if(visSprinkPart==1)parti->display_droplet=1;
    }
    */
    if(visSprinkPart==1||visSmokePart!=0){
      visSmoke=1;
    }
    else{
      visSmoke=0;
    }
  }
  else{
  //  visSmokePart=0; 
  //  visSprinkPart=0;
    switch (value){
      case 1: 
        visSmokePart = 2; 
        break;
      case 2: 
        visSprinkPart = 1; 
        break;
      case 3: 
        visSprinkPart=1; 
        visSmokePart=2; 
        visStaticSmoke=1; 
        for(i=0;i<npart_files;i++){
          parti = partinfo + i;
          if(parti->loaded==0)continue;
          parti->display=1;
        }
        break;
      case 5: 
        visStaticSmoke=1; 
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
    /*
    for(i=0;i<npart_files;i++){
      parti = partinfo + i;
      if(parti->loaded==0||parti->evac==1)continue;
      parti->display_droplet=0;
      parti->display_smoke=0;
      if(visSmokePart!=0)parti->display_smoke=1;
      if(visSprinkPart==1)parti->display_droplet=1;
    }
    */
    if(visSmokePart!=0||visSprinkPart==1){
      visSmoke=1;
    }
  }
  updatemenu=1;  
  plotstate=getplotstate(DYNAMIC_PLOTS);
  GLUTPOSTREDISPLAY
}

/* ------------------ FrameRateMenu ------------------------ */

//void keyboard(unsigned char key, int x, int y);

void FrameRateMenu(int value){
  updateUpdateFrameRateMenu=0;
  realtime_flag=0;
  frameinterval=1;
  if(value > 0){
    switch (value){
    case 2001:
      if(ntimes>0){
        if(times!=NULL)frameinterval=1000.*(times[ntimes-1]-times[0])/ntimes;
      }
      realtime_flag=1;
      break;
    case 2002:
      if(times!=NULL)frameinterval=1000.*(times[ntimes-1]-times[0])/ntimes;
      frameinterval /= 2.0;
      realtime_flag=2;
      break;
    case 2004:
      if(times!=NULL)frameinterval=1000.*(times[ntimes-1]-times[0])/ntimes;
      frameinterval /= 4.0;
      realtime_flag=4;
      break;
    default:
      frameinterval = 1000./value;
      if(frameinterval<1.0){frameinterval = 0.0;}
      break;
    }
    if(times==NULL&&realtime_flag!=0)updateUpdateFrameRateMenu=1;
  }
  else{
    keyboard('t',0,0);
    RenderState(0);
    FlowDir=1;
  }
  frameratevalue=value;
  updatemenu=1;  
  if(opengldefined==1){
    GLUTPOSTREDISPLAY
  }
  reset_gltime();
}

/* ------------------ sv_FrameRateMenu ------------------------ */

void svWINAPI sv_FrameRateMenu(int value){
  FrameRateMenu(value);
}

/* ------------------ IsoSurfaceTypeMenu ------------------------ */

void IsoSurfaceTypeMenu(int value){
  if(ReadPlot3dFile==1){
    switch (value){
    case 0:
      p3dsurfacesmooth=1;
      p3dsurfacetype=1;
      break;
    case 1:
      p3dsurfacesmooth=0;
      p3dsurfacetype=1;
      break;
    case 2:
      p3dsurfacetype=2;
      break;
    case 3:
      p3dsurfacetype=3;
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
    update_glui_plot3dtype();
    updatemenu=1;  
    GLUTPOSTREDISPLAY
  }
}

/* ------------------ IsoBlockMenu ------------------------ */

void IsoBlockMenu(int value){
  if(ReadPlot3dFile==1){
    isooffset=value;
    if(isooffset<1)isooffset=offsetmax;
    if(isooffset>offsetmax)isooffset=1;
    updatesurface();
    updateshowstep(1,ISO);
    updatemenu=1;  
    GLUTPOSTREDISPLAY
  }
}

/* ------------------ IsoSurfaceMenu ------------------------ */

void IsoSurfaceMenu(int value){
  if(ReadPlot3dFile==1){
    updatemenu=1;  
    GLUTPOSTREDISPLAY
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
    GLUTPOSTREDISPLAY
  }
}

/* ------------------ HelpMenu ------------------------ */

void HelpMenu(int value){
  switch (value){
    case -1:
#ifdef pp_OSX
      system("open http://code.google.com/p/fds-smv/issues/");
#endif
#ifdef WIN32
      ShellExecute(NULL, "open", "http://code.google.com/p/fds-smv/issues/", NULL, NULL, SW_SHOWNORMAL);
#endif
      break;
    case -2:
#ifdef pp_OSX
      system("open http://code.google.com/p/fds-smv/downloads/");
#endif
#ifdef WIN32
      ShellExecute(NULL, "open", "http://code.google.com/p/fds-smv/downloads/", NULL, NULL, SW_SHOWNORMAL);
#endif
      break;
    case -3:
#ifdef pp_OSX
      system("open http://fire.nist.gov/fds/documentation.html");
#endif
#ifdef WIN32
      ShellExecute(NULL, "open", "http://fire.nist.gov/fds/documentation.html", NULL, NULL, SW_SHOWNORMAL);
#endif
      break;
    case -4:
#ifdef pp_OSX
      system("open http://fire.nist.gov/fds/");
#endif
#ifdef WIN32
      ShellExecute(NULL, "open", "http://fire.nist.gov/fds/", NULL, NULL, SW_SHOWNORMAL);
#endif
      break;
  }
}


/* ------------------ VectorSkipMenu ------------------------ */

void VectorSkipMenu(int value){
  if(value==-1)return; /* dummy label in menu */
  if(value==-2){       /* toggle vector visibility */
    visVector=1-visVector;
    if(vectorspresent==0)visVector=0;
    updatemenu=1;  
    GLUTPOSTREDISPLAY
    return;
  }
  vectorskip=value;
  visVector=1;
  updatemenu=1;  
  GLUTPOSTREDISPLAY
}

/* ------------------ TextureShowMenu ------------------------ */

void TextureShowMenu(int value){
  texture *texti;
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
    switch (value){
    case -1:
      for(i=0;i<ntextures;i++){
        texti = textureinfo + i;
        if(texti->loaded==0||texti->used==0)continue;
        texti->display=1;
        texturedisplay=1;
      }
      showall_textures=1;
      break;
    case -2:
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
  if(texturedisplay==1){
    if(value!=visBLOCKOutline&&value!=visBLOCKHide){
      BlockageMenu(visBLOCKAsInput);
    }
  }
  updatemenu=1;
  GLUTPOSTREDISPLAY

}

/* ------------------ Plot3DShowMenu ------------------------ */

void Plot3DShowMenu(int value){
  mesh *meshi;
  int i;

  meshi=current_mesh;
  switch (value){
   case 1:
    updateshowstep(1-meshi->visz,DIRZ);
    break;
   case 2:
    updateshowstep(1-meshi->visy,DIRY);
    break;
   case 3:
    updateshowstep(1-meshi->visx,DIRX);
    break;
   case 4:
    switch (p3cont2d){
     case 0:
      p3cont2d=1;
      break;
     case 1:
     case 2:
      p3cont2d=0;
      break;
     default:
       ASSERT(FFALSE);
       break;
    }
    break;
   case 5:
    meshi->visx=1;
    meshi->visy=1;
    meshi->visz=1;
    break;
   case 6:
    meshi->visx=0;
    meshi->visy=0;
    meshi->visz=0;
    plotstate=DYNAMIC_PLOTS;
    break;
   case 7:
    visVector=1-visVector;
    if(vectorspresent==0)visVector=0;
    break;
   case HIDEALL_PLOT3D:
     for(i=0;i<nplot3d_files;i++){
       if(plot3dinfo[i].loaded==1)plot3dinfo[i].display=0;
     }

     break;
   case SHOWALL_PLOT3D:
     for(i=0;i<nplot3d_files;i++){
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
  plotstate=getplotstate(STATIC_PLOTS);
  if(plotstate==STATIC_PLOTS&&visiso==1){
    updatesurface();
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY
}


/* ------------------ GridSliceMenu ------------------------ */

void GridSliceMenu(int value){
  mesh *meshi;
  int justTurnedOn=0;

  meshi=current_mesh;
  if(visGrid==0){
    visGrid=1;
    justTurnedOn=1;
  }
  switch (value){
  case 1:
    if(justTurnedOn==1&&meshi->visz==1)meshi->visz=0;
    updateshowstep(1-meshi->visz,DIRZ);
    if(meshi->visz==1)visGrid=1;
    break;
  case 2:
    if(justTurnedOn==1&&meshi->visy==1)meshi->visy=0;
    updateshowstep(1-meshi->visy,DIRY);
    if(meshi->visy==1)visGrid=1;
    break;
  case 3:
    if(justTurnedOn==1&&meshi->visx==1)meshi->visx=0;
    updateshowstep(1-meshi->visx,DIRX);
    if(meshi->visx==1)visGrid=1;
    break;
  case 4:
    meshi->visx=1;
    meshi->visy=1;
    meshi->visz=1;
    visGrid=1;
    break;
  case 5:
    meshi->visx=0;
    meshi->visy=0;
    meshi->visz=0;
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  if(meshi->visx==0&&meshi->visy==0&&meshi->visz==0){
    visGrid=0;
  }
  togglegridstate(visGrid);
  updatemenu=1;  
  GLUTPOSTREDISPLAY
}

#ifdef pp_COMPRESS
/* ------------------ CompressMenu ------------------------ */

void CompressMenu(int value){
  if(value==999)return;
  switch (value){
  case 1:
    erase_all=1;
    overwrite_all=0;
    update_overwrite();
    compress_svzip();
    break;
  case 2:
    erase_all=0;
    overwrite_all=1-overwrite_all;
    update_overwrite();
    break;
  case 3:
    erase_all=0;
    compress_svzip();
    break;
  case 4:
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
  if(value==-1){
    readini(0);
  }
  else{
    char *inifilename;

    inifilename = get_inifilename(value);
    if(inifilename==NULL||strlen(inifilename)<=0)return;
    scriptinifilename2=scriptinifilename;
    strcpy(scriptinifilename,inifilename);
    windowresized=0;
    readini(2);
    scriptinifilename2=NULL;
  }
}

/* ------------------ SmokeviewiniMenu ------------------------ */

void SmokeviewiniMenu(int value){
  switch (value){
  case 1:
    readini(0);
    updatecolors(-1);
    updateshowtitles();
    break;
  case 2:
    writeini(GLOBAL_INI);
    break;
  case 3:
    writeini(LOCAL_INI);
    break;
  case 4:
    init_object_defs();
    break;
  case 999:
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY
}


/* ------------------ PeriodicReloads ------------------------ */

void PeriodicReloads(int value){
  if(periodic_reloads!=0){
    LoadUnloadMenu(RELOADALL);
    glutTimerFunc((unsigned int)value,PeriodicReloads,value);
  }
}

/* ------------------ ScriptMenu ------------------------ */

void ScriptMenu(int value){
  int error_code;
  scriptfiledata *scriptfile;
  char newscriptfilename[1024];

  if(value==999)return;
  updatemenu=1;
  GLUTPOSTREDISPLAY
  switch (value){
    case START_RECORDING_SCRIPT:
      get_newscriptfilename(newscriptfilename);
      script_recording = insert_scriptfile(newscriptfilename);
      scriptoutstream=fopen(newscriptfilename,"w");
      if(scriptoutstream!=NULL){
        printf("Script recorder on\n");
        script_recording->recording=1;
        {
          char *renderdir;

          trim(script_renderdir);
          renderdir = trim_front(script_renderdir);
          if(strlen(renderdir)>0&&strcmp(renderdir,".")!=0){
            fprintf(scriptoutstream,"RENDERDIR\n");
            fprintf(scriptoutstream," %s\n",renderdir);
          }
          else{
            fprintf(scriptoutstream,"RENDERDIR\n");
            fprintf(scriptoutstream," .\n");
          }
        }
      }
      else{
        script_recording->recording=0;
        script_recording=NULL;
#ifdef pp_MESSAGE
        {
          char message[256];
          sprintf(message,"*** warning: The script file %s could not be opened for writing\n",newscriptfilename);
          warning_message(message);
        }
#else
        printf("*** warning: The script file %s could not be opened for writing\n",newscriptfilename);
#endif
      }
      break;
    case STOP_RECORDING_SCRIPT:
      if(script_recording!=NULL){
        script_recording->recording=0;
        add_scriptlist(script_recording->file,script_recording->id);
        script_recording=NULL;
      }
      if(scriptoutstream!=NULL){
        fclose(scriptoutstream);
        scriptoutstream=NULL;
        writeini(SCRIPT_INI);
        printf("Script recorder off\n");
      }
      break;
    default:
      for(scriptfile=first_scriptfile.next;scriptfile->next!=NULL;scriptfile=scriptfile->next){
        char *file;

        file=scriptfile->file;
        if(file==NULL)continue;
        if(scriptfile->id!=value)continue;
        error_code=compile_script(file);
        if(error_code==0){
          readini(1);
          start_script();
        }
        break;
      }
      break;
  }
}

/* ------------------ ReLoadMenu ------------------------ */

void ReloadMenu(int value){
  int msecs;

  updatemenu=1;
  periodic_value=value;
  switch (value){
  case -1:
    periodic_reloads=0;
    break;
  case 0:
    LoadUnloadMenu(RELOADALL);
    break;
  default:
    periodic_reloads=1;
    msecs = value*60*1000;
    glutTimerFunc((unsigned int)msecs,PeriodicReloads,msecs);
    break;
  }
}


/* ------------------ DisclaimerMenu ------------------------ */

void DisclaimerMenu(int value){
}

/* ------------------ LoadUnloadMenu ------------------------ */

void LoadUnloadMenu(int value){
  int errorcode;
  size_t len;
  int i;
  int ii;

  if(value==999)return;
  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value==UNLOADALL){
   // leaving code here commented in case I later decide to unload terrain files
   // for(i=0;i<nterraininfo;i++){
   //   readterrain("",i,UNLOAD,&errorcode);
   // }
    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"UNLOADALL\n");
    }
    if(hrrfilename!=NULL){
      readhrr(UNLOAD, &errorcode);
    }
    for(i=0;i<nslice_files;i++){
      readslice("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<nplot3d_files;i++){
      readplot3d("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<npatch_files;i++){
      readpatch(i,UNLOAD,&errorcode);
    }
    for(i=0;i<npart_files;i++){
      readpart("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<niso_files;i++){
      readiso("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<nzone;i++){
      readzone("",i,UNLOAD,&errorcode);
    }
    for(i=0;i<nsmoke3d_files;i++){
      readsmoke3d(i,UNLOAD,&errorcode);
    }
    //plotstate=DYNAMIC_PLOTS;
    //visSmoke=1;
    updatemenu=1;  
    GLUTPOSTREDISPLAY
  }
  if(value==RELOADALL){
    LOCK_COMPRESS
    if(hrrfilename!=NULL){
      readhrr(LOAD, &errorcode);
    }
    islicetype_save=islicetype;
    for(i=0;i<nslice_files;i++){
      sliceinfo[i].reload=1;
    }
    for(i=0;i<nterraininfo;i++){
      if(terraininfo[i].loaded==1){
        readterrain(terraininfo[i].file,i,LOAD,&errorcode);
      }
    }
    for(i=0;i<nvslice;i++){
      if(vsliceinfo[i].loaded==1){
        readvslice(i,LOAD,&errorcode);
      }
    }
    for(ii=0;ii<nslice_loaded;ii++){
      i = slice_loaded_list[ii];
      if(sliceinfo[i].reload==1){
        readslice(sliceinfo[i].file,i,LOAD,&errorcode);
      }
    }
    islicetype=islicetype_save;
    for(i=0;i<nplot3d_files;i++){
      if(plot3dinfo[i].loaded==1){
        readplot3d(plot3dinfo[i].file,i,LOAD,&errorcode);
      }
    }
    for(ii=0;ii<npatch_loaded;ii++){
      i = patch_loaded_list[ii];
      readpatch(i,LOAD,&errorcode);
    }
    for(i=0;i<nsmoke3d_files;i++){
      if(smoke3dinfo[i].loaded==1){
        readsmoke3d(i,LOAD,&errorcode);
      }
    }
    for(i=0;i<npart_files;i++){
      if(partinfo[i].loaded==1){
        readpart(partinfo[i].file,i,LOAD,&errorcode);
      }
    }
    for(i=0;i<niso_files;i++){
      iso *isoi;

      isoi = isoinfo + i;
      if(isoi->loaded==0)continue;
      readiso(isoi->file,i,LOAD,&errorcode);
    }
    UNLOCK_COMPRESS
  //  plotstate=DYNAMIC_PLOTS;
  //  visSmoke=1;
    updatemenu=1;  
    GLUTPOSTREDISPLAY
  }
  if(value==2&&smvfilename!=NULL){
#ifdef _DEBUG
    CheckMemory;
#endif
    sv_unload();
    FREEMEMORY(fdsprefix2);
    len=0;
    if(fdsprefix!=NULL)len=strlen(fdsprefix);
    if(len!=0&&NewMemory((void **)&fdsprefix2,(unsigned int)(len+1))!=0){
      STRCPY(fdsprefix2,fdsprefix);
      sv_startup(fdsprefix2,0);
    }
  }
#ifdef pp_OPEN
  if(value==3){
    OpenSMVFile(openfilebuffer,1024,&openfileflag);
    if(openfileflag==1){
      updateOpenSMVFile=1;
    }
    else{
      updateOpenSMVFile=0;
    }
    if(updateOpenSMVFile==1){
      DialogMenu(-2); // hide all glui menus
      sv_unload();
      PrintAllMemoryInfo;
      sv_init0();
      sv_startup(openfilebuffer,0);
      updateOpenSMVFile=0;
    }
  }
#endif
  if(value==SHOWFILES){
    GLUTPOSTREDISPLAY  
    showfiles=1-showfiles;
    updatemenu=1;
    updateslicemenulabels();
    updatevslicemenulabels();
   // updatesmokemenulabels();
    updatesmoke3dmenulabels();
    updatepatchmenulabels();
    updateisomenulabels();
    updatepartmenulabels();
    updatetourmenulabels();
    updateplot3dmenulabels();
  }
  glutSetCursor(GLUT_CURSOR_RIGHT_ARROW);
}


/* ------------------ ShowTourMenu ------------------------ */

void ShowTourMenu(int value){
}


/* ------------------ AvatarTourMenu ------------------------ */

void AvatarTourMenu(int value){
}

void AvatarEvacMenu(int value){
  if(value==-999)return;
  iavatar_evac=value;
  updatemenu=1;
}

/* ------------------ TourMenu ------------------------ */

void TourMenu(int value){
  tourdata *touri;
  int i;

  if(value==-999)return;
  touring=0;
  updatemenu=1;
  GLUTPOSTREDISPLAY
  switch (value){
  case -12:
    add_new_tour();
    if(showgluitour==0){
      showgluitour=1;
      show_glui_tour();
    }
    break;
  case -13:               // clear all tours
    for(i=0;i<ntours;i++){
      touri = tourinfo + i;
      touri->display=touri->display2;
    }
    if(viewtourfrompath==1
//      &&from_glui_trainer==0
      ){
      ResetView(RESTORE_EXTERIOR_VIEW);
    }
    from_glui_trainer=0;
//    viewtourfrompath=0;
    for(i=0;i<ntours;i++){
      touri = tourinfo + i;
      if(touri->display==1){
        selected_tour=touri;
        break;
      }
    }
    selected_tour=NULL;
    break;
  case -2:               // clear all tours
    for(i=0;i<ntours;i++){
      touri = tourinfo + i;
      touri->display=0;
    }
    if(viewtourfrompath==1
//      &&from_glui_trainer==0
      ){
      ResetView(RESTORE_EXTERIOR_VIEW);
    }
    from_glui_trainer=0;
//    viewtourfrompath=0;
    selected_tour=NULL;
    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"UNLOADTOUR\n");
    }
    break;
  case -4:
    edittour=1-edittour;
    if(edittour==1&&showgluitour==0){
      showgluitour=1;
      show_glui_tour();
    }
    break;
  case -3:               // show all tours
    for(i=0;i<ntours;i++){
      touri = tourinfo + i;
      touri->display=1;
    }
    plotstate=getplotstate(DYNAMIC_PLOTS);
    break;
  case -5:               // view from route
    viewtourfrompath = 1 - viewtourfrompath;
    if(viewtourfrompath==0)ResetView(RESTORE_EXTERIOR_VIEW);
    break;
  case -6:
    tour_constant_vel=1-tour_constant_vel;
    createtourpaths();
    updatetimes();
    break;
  case -1:
    for(i=0;i<ntours;i++){
      touri = tourinfo + i;
      touri->display=0;
    }
    ResetView(RESTORE_EXTERIOR_VIEW);
    defaulttour();
    break;
  case -8:  // crawl
    eyeview_level=0;
    pass_through=0;
    update_glui_speed();
    if(trainer_mode==1){
      eyeview=EYE_CENTERED;
      handle_eyeview(0);
    }
    break;
  case -9:  // walk
    eyeview_level=1;
    pass_through=0;
    update_glui_speed();
    if(trainer_mode==1){
      eyeview=EYE_CENTERED;
      handle_eyeview(0);
    }
    break;
  case -10: // overview
    eyeview_level=2;
    pass_through=1;
    update_glui_speed();
    if(trainer_mode==1){
      eyeview=EYE_CENTERED;
      handle_eyeview(0);
    }
    break;
  case -11: // bird's eye
    break;
  case -21:
    tourlocus_type=0;
    break;
  case -22:
    tourlocus_type=1;
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
  plotstate=getplotstate(DYNAMIC_PLOTS);
  if(value!=-5&&value!=-4)updatetimes();
  callfrom_tourglui=0;

}

/* ------------------ SetTour ------------------------ */

void SetTour(tourdata *thetour){
  int tournumber;

  if(thetour==NULL)return;
  tournumber = thetour - tourinfo;
  TourMenu(tournumber);
}

/* ------------------ targetMenu ------------------------ */

void TargetMenu(int value){
  int errorcode,i;
  if(value>=0){
    readtarget(targinfo[value].file,value,LOAD,&errorcode);
  }
  else{
    if(value==-1){
      for(i=0;i<ntarg_files;i++){
        readtarget("",i,UNLOAD,&errorcode);
      }
    }
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY
}

/* ------------------ EvacMenu ------------------------ */

void EvacMenu(int value){
  int errorcode;
  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value==-11){
    int i;

    for(i=0;i<npart_files;i++){
      particle *parti;

      parti=partinfo + i;
      if(parti->evac==0)continue;
      ReadEvacFile=1;
      readpart(parti->file,i,LOAD,&errorcode);
      if(scriptoutstream!=NULL){
        fprintf(scriptoutstream,"LOADFILE\n");
        fprintf(scriptoutstream," %s\n",parti->file);
      }
    }
    force_redisplay=1;
    update_framenumber(0);
  }
  if(value>=0){
    ReadEvacFile=1;
    readpart(partinfo[value].file,value,LOAD,&errorcode);
    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"LOADFILE\n");
      fprintf(scriptoutstream," %s\n",partinfo[value].file);
    }
  }
  else{
    if(value==-1){
      int i;

      for(i=0;i<npart_files;i++){
        if(partinfo[i].evac==0)continue;
        readpart("",i,UNLOAD,&errorcode);
      }
    }
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ update_streakvalue ------------------------ */

void update_streakvalue(float value){
  particle *parti=NULL;
  int i;

  streak_index=-1;
  for(i=0;i<nstreak_value;i++){
    float absdiff;

    absdiff = value-streak_rvalue[i];
    if(absdiff<0.0)absdiff=-absdiff;

    if(absdiff<0.01){
      streak_index=i;
      float_streak5value=streak_rvalue[i];
      break;
    }
  }
  for(i=0;i<npart_files;i++){
    parti = partinfo + i;
    if(parti->loaded==1)break;
  }
  if(parti!=NULL&&parti->nframes>1){
    for(i=0;i<parti->nframes-1;i++){
      if(parti->ptimes[i]<=value&&value<parti->ptimes[i+1]){
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
  if(value==-2){
    streak5show=0;
    streak5step=0;
  }
  else if(value==-3){
    showstreakhead=1-showstreakhead;
  }
  else{
    streak5show=1;
    streak5step=0;
    rvalue=streak_rvalue[value];
    update_streakvalue(rvalue);
    update_glui_streakvalue(rvalue);

  }
  updatemenu=1;
  GLUTPOSTREDISPLAY
}

/* ------------------ Particle5ShowMenu ------------------------ */

void Particle5ShowMenu(int value){
}

/* ------------------ ParticlePropShowMenu ------------------------ */

void ParticlePropShowMenu(int value){
  part5prop *propi;

  int propvalue, showvalue;

  propvalue = (-value)/10000-1;
  value = -((-value)%10000);

  if(value==-1)return;

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
  else if(value==-2){
    if(current_property!=NULL){
      unsigned char *vis;
      int i;

      vis = current_property->class_vis;
      for(i=0;i< npartclassinfo;i++){
        vis[i]=1;
      }
    }
  }
  else if(value==-3){
    if(current_property!=NULL){
      unsigned char *vis;
      int i;

      vis = current_property->class_vis;
      for(i=0;i< npartclassinfo;i++){
        vis[i]=0;
      }
    }

  }
  else if(value==-4){
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
  else if(value==-5){
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
  else if(value==-6){
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
      part5class *partclassj;

      partclassj = partclassinfo + iclass;
      partclassj->vis_type=vistype;
      PropMenu(propvalue);
    }

  }
  updatemenu=1;
  GLUTPOSTREDISPLAY
}

/* ------------------ ParticleMenu ------------------------ */

void ParticleMenu(int value){
  int errorcode,i;
  int whichpart;
  particle *parti, *partj;

  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value>=0){
    char  *partfile;

    ReadPartFile=1;
    partfile = partinfo[value].file;
    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"LOADFILE\n");
      fprintf(scriptoutstream," %s\n",partfile);
    }
    readpart(partfile,value,LOAD,&errorcode);
  }
  else{
    if(value==-1){
      for(i=0;i<npart_files;i++){
        if(partinfo[i].evac==1)continue;
        readpart("",i,UNLOAD,&errorcode);
      }
    }
    else{
      ReadPartFile=1;
      whichpart=-(10+value);
      partj = partinfo + whichpart;
      if(scriptoutstream!=NULL){
        fprintf(scriptoutstream,"LOADPARTICLES\n");
      }
      for(i=0;i<npart_files;i++){
        parti = partinfo + i;
        if(parti->evac==1)continue;
        if(parti->version==1||strcmp(parti->label.longlabel,partj->label.longlabel)==0){
          readpart(parti->file,i,LOAD,&errorcode);
        }
      }
      force_redisplay=1;
      update_framenumber(0);
    }
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ ZoneMenu ------------------------ */

void ZoneMenu(int value){
  int i,errorcode;
  if(value>=0){
    if(scriptoutstream!=NULL){
      char *file;

      file = zoneinfo[value].file;
      fprintf(scriptoutstream,"LOADFILE\n");
      fprintf(scriptoutstream," %s\n",file);
    }
    readzone(zoneinfo[value].file,value,LOAD,&errorcode);
  }
  else{
    for(i=0;i<nzone;i++){
      readzone("",i,UNLOAD,&errorcode);
    }
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY
}

/* ------------------ UnloadVSliceMenu ------------------------ */

void UnloadVSliceMenu(int value){
  int errorcode,i;

  updatemenu=1;  
  GLUTPOSTREDISPLAY
  if(value>=0){
    readvslice(value,UNLOAD,&errorcode);
  }
  else{
    for(i=0;i<nvslice;i++){
      readvslice(i,UNLOAD,&errorcode);
    }
  }
}

/* ------------------ UnloadPatchMenu ------------------------ */

void UnloadPatchMenu(int value){
  int errorcode,i;

  updatemenu=1;  
  GLUTPOSTREDISPLAY
  if(value>=0){
    readpatch(value,UNLOAD,&errorcode);
  }
  else{
    for(i=0;i<npatch_files;i++){
      readpatch(i,UNLOAD,&errorcode);
    }
  }
}

/* ------------------ UnloadIsoMenu ------------------------ */

void UnloadIsoMenu(int value){
  int errorcode,i;

  updatemenu=1;  
  GLUTPOSTREDISPLAY
  if(value>=0){
    readiso("",value,UNLOAD,&errorcode);
  }
  else{
    for(i=0;i<niso_files;i++){
      readiso("",i,UNLOAD,&errorcode);
    }
  }
}

/* ------------------ UnloadPlot3dMenu ------------------------ */

void UnloadPlot3dMenu(int value){
  int errorcode,i;

  updatemenu=1;  
  GLUTPOSTREDISPLAY
  if(value>=0){
    readplot3d("",value,UNLOAD,&errorcode);
  }
  else{
    for(i=0;i<nplot3d_files;i++){
      readplot3d("",i,UNLOAD,&errorcode);
    }
  }
}

/* ------------------ UnloadEvacMenu ------------------------ */

void UnloadEvacMenu(int value){
  int errorcode,i;

  updatemenu=1;  
  GLUTPOSTREDISPLAY
  if(value>=0){
    readpart("",value,UNLOAD,&errorcode);
  }
  else{
    for(i=0;i<npart_files;i++){
      if(partinfo[i].evac==0)continue;
      readpart("",i,UNLOAD,&errorcode);
    }
  }
}

/* ------------------ UnloadPartMenu ------------------------ */

void UnloadPartMenu(int value){
  int errorcode,i;

  updatemenu=1;  
  GLUTPOSTREDISPLAY
  if(value>=0){
    readpart("",value,UNLOAD,&errorcode);
  }
  else{
    for(i=0;i<npart_files;i++){
      if(partinfo[i].evac==1)continue;
      readpart("",i,UNLOAD,&errorcode);
    }
  }
}

/* ------------------ LoadVSliceMenu ------------------------ */

void LoadVSliceMenu(int value){
  int errorcode;
  int i;

  if(value==-999)return;
  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value==-1){
    for(i=0;i<nvslice;i++){
      readvslice(i,UNLOAD,&errorcode);
    }
    return;
  }
  if(value==-20){
    showallslicevectors=1-showallslicevectors;
    updatemenu=1;  
    GLUTPOSTREDISPLAY
  }
  if(value>=0){
    vslice *vslicei;
    slice *slicei;

    readvslice(value, LOAD, &errorcode);
    vslicei = vsliceinfo + value;
    slicei = vslicei->val;
    if(script_multivslice==0&&slicei!=NULL&&scriptoutstream!=NULL){
      char *file;

      file=slicei->file;
      fprintf(scriptoutstream,"LOADVFILE\n");
      fprintf(scriptoutstream," %s\n",file);
    }
  }
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ UnloadSliceMenu ------------------------ */

void UnloadSliceMenu(int value){
  int errorcode,i;

  updatemenu=1;  
  GLUTPOSTREDISPLAY
  if(value>=0){
    readslice("",value,UNLOAD,&errorcode);
  }
  else{
    for(i=0;i<nslice_files;i++){
      readslice("",i,UNLOAD,&errorcode);
    }
  }
}

/* ------------------ UnLoadMultiVSliceMenu ------------------------ */

void UnloadMultiVSliceMenu(int value){
  int i;
  multivslice *mvslicei;

  if(value>=0){
    mvslicei = multivsliceinfo + value;
    for(i=0;i<mvslicei->nvslices;i++){
      UnloadSliceMenu(mvslicei->ivslices[i]);
    }
  }
  else{
    UnloadSliceMenu(-1);
  }
}


/* ------------------ UnLoadMultiSliceMenu ------------------------ */

void UnloadMultiSliceMenu(int value){
  int i;
  multislice *mslicei;

  if(value>=0){
    mslicei = multisliceinfo + value;
    for(i=0;i<mslicei->nslices;i++){
      UnloadSliceMenu(mslicei->islices[i]);
    }
  }
  else{
    UnloadSliceMenu(-1);
  }
}

/* ------------------ UnLoadSmoke3DMenu ------------------------ */

void UnLoadSmoke3DMenu(int value){
  int errorcode;
  int i;
  smoke3d *smoke3di;

  if(value==999)return;
  updatemenu=1;
  if(value<0){
    value= -value;
    for(i=0;i<nsmoke3d_files;i++){
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
  smoke3d *smoke3di, *smoke3dj;

  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value>=0){
    if(scriptoutstream!=NULL){
      char *file;

      file = smoke3dinfo[value].file;
      fprintf(scriptoutstream,"LOADFILE\n");
      fprintf(scriptoutstream," %s\n",file);
    }
  readsmoke3d(value,LOAD,&errorcode);

  }
  else if(value==-1){
    for(i=0;i<nsmoke3d_files;i++){
      readsmoke3d(i,UNLOAD,&errorcode);
    }
  }
  if(value==-9){
    lock_allsmoke=1;
    for(i=0;i<nsmoke3d_files;i++){
      smoke3di = smoke3dinfo + i;
      if(smoke3di->loaded==1)continue;
      readsmoke3d(i,LOAD,&errorcode);
    }
    lock_allsmoke=0;
  }
  if(value<=-10){
    value = -(value + 10);
    smoke3dj = smoke3dinfo + value;
    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"LOAD3DSMOKE\n");
      fprintf(scriptoutstream," %s\n",smoke3dj->label.longlabel);
    }
    for(i=0;i<nsmoke3d_files;i++){
      smoke3di = smoke3dinfo + i;
      if(strcmp(smoke3di->label.shortlabel,smoke3dj->label.shortlabel)==0){
        readsmoke3d(i,LOAD,&errorcode);
      }
    }
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ AnySmoke ------------------------ */

int AnySmoke(char *type){

  if(nsmoke3d_files>0)return 1;
  return 0;
}

/* ------------------ AnySlices ------------------------ */

int AnySlices(char *type){
  int i;

  for(i=0;i<nslice_files;i++){
    if(STRCMP(sliceinfo[i].label.longlabel,type)==0)return 1;
  }
  return 0;
}

/* ------------------ LoadAllSlices ------------------------ */

void HideAllSlices(void){
  int i;

  glutSetCursor(GLUT_CURSOR_WAIT);
  for(i=0;i<nslice_files;i++){
    sliceinfo[i].display=0;
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ ShowAllSmoke ------------------------ */

void ShowAllSmoke(void){
  int i;
  for(i=0;i<nsmoke3d_files;i++){
    smoke3d *smoke3di;

    smoke3di = smoke3dinfo + i;
    if(smoke3di->loaded==1)smoke3di->display=1;
  }
  for(i=0;i<niso_files;i++){
    iso *isoi;

    isoi = isoinfo + i;
    if(isoi->loaded==1)isoi->display=1;
  }
}

/* ------------------ HideAllSmoke ------------------------ */

void HideAllSmoke(void){
  int i;
  for(i=0;i<nsmoke3d_files;i++){
    smoke3d *smoke3di;

    smoke3di = smoke3dinfo + i;
    if(smoke3di->loaded==1)smoke3di->display=0;
  }
  for(i=0;i<niso_files;i++){
    iso *isoi;

    isoi = isoinfo + i;
    if(isoi->loaded==1)isoi->display=0;
  }
}

/* ------------------ ShowAllSlices ------------------------ */

void ShowAllSlices(char *type1, char *type2){
  int i;

  glutSetCursor(GLUT_CURSOR_WAIT);
  for(i=0;i<nslice_files;i++){
    sliceinfo[i].display=0;
    if(sliceinfo[i].loaded==0)continue;
    if(
      type1!=NULL&&STRCMP(sliceinfo[i].label.longlabel,type1)==0||
      type2!=NULL&&STRCMP(sliceinfo[i].label.longlabel,type2)==0
      ){
      sliceinfo[i].display=1;
      islicetype=sliceinfo[i].type;
    }
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

void UnloadTerrainMenu(int value);
void LoadTerrainMenu(int value);

/* ------------------ LoadTerrainMenu ------------------------ */

void LoadTerrainMenu(int value){
  int i;
  int errorcode;

  if(value>=0&&value<nterraininfo){
    terraindata *terri;

    terri = terraininfo + value;
    readterrain(terri->file,value,LOAD,&errorcode);
  }
  else if(value==-10){
    UnloadTerrainMenu(value);
  }
  else if(value==-9){
    for(i=0;i<nterraininfo;i++){
      LoadTerrainMenu(i);
    }
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY
}

/* ------------------ UnLoadTerrainMenu ------------------------ */

void UnloadTerrainMenu(int value){
  int i;
  int errorcode;

  if(value>=0&&value<nterraininfo){
    readterrain("",value,UNLOAD,&errorcode);
  }
  else if(value==-10){
    for(i=0;i<nterraininfo;i++){
      UnloadTerrainMenu(i);
    }
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY

}

/* ------------------ LoadSliceMenu ------------------------ */

void LoadSliceMenu(int value){
  int errorcode,i;

  if(value==-999)return;
  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value>=0){
    char *file;

    file = sliceinfo[value].file;
    if(script_multislice==0&&scriptoutstream!=NULL){
      fprintf(scriptoutstream,"LOADFILE\n");
      fprintf(scriptoutstream," %s\n",file);
    }
    readslice(file,value,LOAD,&errorcode);
  }
  else{
    for(i=0;i<nslice_files;i++){
      readslice("",i,UNLOAD,&errorcode);
    }
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ LoadVMultiSliceMenu ------------------------ */

void LoadMultiVSliceMenu(int value){
  int i;
  multivslice *mvslicei;

  if(value==-999)return;
  if(value==-20){
    showallslicevectors=1-showallslicevectors;
    updatemenu=1;  
    GLUTPOSTREDISPLAY
    return;
  }
  if(value>=0){
    mvslicei = multivsliceinfo + value;
    if(scriptoutstream!=NULL){
      if(mvslicei->nvslices>0){
        slice *slicei;

        slicei = sliceinfo + mvslicei->ivslices[0];
        fprintf(scriptoutstream,"LOADVSLICE\n");
        fprintf(scriptoutstream," %s\n",slicei->label.longlabel);
        fprintf(scriptoutstream," %i %f\n",slicei->idir,slicei->position);
        script_multivslice=1;
      }
    }
    for(i=0;i<mvslicei->nvslices;i++){
      LoadVSliceMenu(mvslicei->ivslices[i]);
    }
    script_multivslice=0;
  }
  else{
    LoadVSliceMenu(-1);
  }
}


/* ------------------ LoadMultiSliceMenu ------------------------ */

void LoadMultiSliceMenu(int value){
  int i;
  multislice *mslicei;

  if(value==-999)return;
  if(value>=0){
    mslicei = multisliceinfo + value;
    if(scriptoutstream!=NULL){
      if(mslicei->nslices>0){
        slice *slicei;

        slicei = sliceinfo + mslicei->islices[0];
        fprintf(scriptoutstream,"LOADSLICE\n");
        fprintf(scriptoutstream," %s\n",slicei->label.longlabel);
        fprintf(scriptoutstream," %i %f\n",slicei->idir,slicei->position);
        script_multislice=1;
      }
    }
    for(i=0;i<mslicei->nslices;i++){
      LoadSliceMenu(mslicei->islices[i]);
    } 
    script_multislice=0;
  }
  else{
    LoadSliceMenu(-1);
  }
}

/* ------------------ Plot3DListMenu ------------------------ */

void Plot3DListMenu(int value){
  int i;
  plot3d *plot3di;

  if(value<0||value>=nplot3dtimelist)return;
  LoadPlot3dMenu(-1);
  if(scriptoutstream!=NULL){
    fprintf(scriptoutstream,"LOADPLOT3D\n");
    fprintf(scriptoutstream," %f\n",plot3dtimelist[value]);
  }
  for(i=0;i<nplot3d_files;i++){
    plot3di = plot3dinfo + i;
    if(fabs(plot3di->time-plot3dtimelist[value])<0.5){
      LoadPlot3dMenu(i);
    }
  }
}

/* ------------------ update_menu ------------------------ */

void update_menu(void){
  updatemenu=1; 
  GLUTPOSTREDISPLAY
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ LoadPlot3DMenu ------------------------ */

void LoadPlot3dMenu(int value){
  int errorcode;
  int i;

  if(value==997)return;
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
    readplot3d(plot3dfile,value,LOAD,&errorcode);
  }
  else if(value==-1){
    for(i=0;i<nplot3d_files;i++){
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
  GLUTPOSTREDISPLAY
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ LoadIsoMenu ------------------------ */

void LoadIsoMenu(int value){
  int errorcode;
  int i;
  int ii;
  iso *isoii, *isoi;

  if(value==-2)return;
  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value>=0){
    char *file;

    ReadIsoFile=1;
    file=isoinfo[value].file;
    if(script_iso==0&&scriptoutstream!=NULL){
      fprintf(scriptoutstream,"LOADFILE\n");
      fprintf(scriptoutstream," %s\n",file);
    }
    readiso(file,value,LOAD,&errorcode);
  }
  if(value==-1){
    for(i=0;i<niso_files;i++){
      isoii = isoinfo + i;
      if(isoii->loaded==1)readiso("",i,UNLOAD,&errorcode);
    }
  }
  if(value<=-10){
    ii = -(value + 10);
    isoii = isoinfo + ii;
    if(scriptoutstream!=NULL){
      script_iso=1;
      fprintf(scriptoutstream,"LOADISO\n");
      fprintf(scriptoutstream," %s\n",isoii->surface_label.longlabel);
    }
    for(i=0;i<niso_files;i++){
      isoi = isoinfo + i;
      if(isoii->type!=isoi->type)continue;
      LoadIsoMenu(i);
    }
    script_iso=0;
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ LoadPatchMenu ------------------------ */

void LoadPatchMenu(int value){
  int errorcode;
  int i,ii;
  int patchtypenew;
  patch *patchi;

  glutSetCursor(GLUT_CURSOR_WAIT);
  if(value>=0){
    patchtypenew=getpatchtype(patchinfo+value);
    if(patchtypenew!=-1){
      for(ii=0;ii<npatch_loaded;ii++){
        i = patch_loaded_list[ii];
        patchi = patchinfo + i;
        if(patchi->type!=patchtypenew)readpatch(i,UNLOAD,&errorcode);
      }
    }
    if(scriptoutstream!=NULL){
      char *file;

      file = patchinfo[value].file;
      fprintf(scriptoutstream,"LOADFILE\n");
      fprintf(scriptoutstream," %s\n",file);
    }
    LOCK_COMPRESS
    readpatch(value,LOAD,&errorcode);
    UNLOCK_COMPRESS
  }
  else if(value<=-10){
    patch *patchj;

    value = -(value + 10);
    patchj = patchinfo + value;
    if(scriptoutstream!=NULL){
      fprintf(scriptoutstream,"LOADBOUNDARY\n");
      fprintf(scriptoutstream," %s\n",patchj->label.longlabel);
    }
    for(i=0;i<npatch_files;i++){
      patchi = patchinfo + i;
      if(strcmp(patchi->label.shortlabel,patchj->label.shortlabel)==0){
        LOCK_COMPRESS
        readpatch(i,LOAD,&errorcode);
        UNLOCK_COMPRESS
      }
    }
    force_redisplay=1;
    update_framenumber(0);
  }
  else{
    for(i=0;i<npatch_files;i++){
      readpatch(i,UNLOAD,&errorcode);
    }
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY
  glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
}

/* ------------------ ShowPatchMenu ------------------------ */

void ShowPatchMenu(int value){
  int val=0,n,i,ii;
//  int sum=0;
  mesh *meshi;
  patch *patchi;
  meshi=current_mesh;

  updatemenu=1;  
  updatefacelists=1;
  GLUTPOSTREDISPLAY
  if(value>=1000){
    patchi=patchinfo+value-1000;
    if(patchi->type==ipatchtype){
      patchi->display=1-patchi->display;
    }
    else{
      patchi->display=1;
      ipatchtype=getpatchtype(patchi);
    }
    update_patchtype();
  }
  if(value==SHOW_CHAR){
    vis_threshold = 1 - vis_threshold;
    updatechar();
  }
  if(value==SHOWALL_BOUNDARY){
    for(ii=0;ii<npatch_loaded;ii++){
      i = patch_loaded_list[ii];
      patchi = patchinfo + i;
      patchi->display=1;
    }
  }
  if(value==HIDEALL_BOUNDARY){
    for(ii=0;ii<npatch_loaded;ii++){
      i = patch_loaded_list[ii];
      patchi = patchinfo + i;
      patchi->display=0;
    }
  }
  if(value<0){
    if(value==EXTERIORwallmenu){
      allexterior = 1-allexterior;
	  showexterior=1-showexterior;
//      for(i=1;i<7;i++){
//        sum+=visPatchType[i];
//      }
      val = allexterior;
      for(n=0;n<meshi->npatches;n++){
        if(meshi->patchtype[n]!=INTERIORwall){
          meshi->visPatches[n]=val;
        }
      }
      for(i=1;i<7;i++){
        visPatchType[i]=val;
      }
    }
    else if(value==INTERIORwallmenu){
      ASSERTFLAG(allinterior);
      allinterior = 1 - allinterior;
      val = allinterior;
      visPatchType[INTERIORwall]=val;
      for(n=0;n<meshi->npatches;n++){
        if(meshi->patchtype[n]==INTERIORwall){
          meshi->visPatches[n]=val;
        }
      }
    }
    else if(value!=DUMMYwallmenu){
      value = -(value+2); /* map xxxwallmenu to xxxwall */
      for(n=0;n<meshi->npatches;n++){
        if(meshi->patchtype[n]==value){
          ASSERTFLAG(meshi->visPatches[n]);
          meshi->visPatches[n] = 1 - meshi->visPatches[n];
          visPatchType[value]=meshi->visPatches[n];
        }
      }
    }
  }
  plotstate=getplotstate(DYNAMIC_PLOTS);
}


/* ------------------ VentMenu ------------------------ */

void VentMenu(int value){

  if(value==-1)return;
  switch (value){
  case 10:
    visVents=1-visVents;
    break;
  case 14:
    visOpenVents=1-visOpenVents;
    break;
  case 15:
    visOpenVentsAsOutline = 1 - visOpenVentsAsOutline;
    break;
  case 16:
    visDummyVents = 1 - visDummyVents;
    break;
   case 18:
     show_bothsides_int=1-show_bothsides_int;
     updatefaces=1;
     break;
   case 19:
     show_bothsides_ext = 1 - show_bothsides_ext;
     updatefaces=1;
     break;
   case 20:
     show_transparent_vents=1-show_transparent_vents;
     updatefaces=1;
     break;
  default:
    ASSERT(FFALSE);
    break;
  }
  updatefacelists=1;
  updatemenu=1;  
  GLUTPOSTREDISPLAY
}

/* ------------------ BlockageMenu ------------------------ */

//     visBlocks (visBLOCKNormal, visBLOCKAsInput )
//     visSmoothAsNormal
//     visTransparentBlockage

void BlockageMenu(int value){
  switch (value){
//   case visBLOCKFacet:
//   case visBLOCKSmooth:
   case visBLOCKAsInput:
     visBlocks=value;
     update_trainer_outline();
     
     break;
   case visSmoothBLOCKSolid:
     smooth_block_solid = 1 - smooth_block_solid;
     break;
   case visBLOCKNormal:
   case visBLOCKOutline:
   case visBLOCKHide:
     visBlocks=value;
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
   case visBLOCKSmoothAsNormal:
     visSmoothAsNormal = 1 - visSmoothAsNormal;
     break;
   case visBLOCKTransparent:
     visTransparentBlockage=1-visTransparentBlockage;
     break;
   case SMOOTH_BLOCKAGES:
     menusmooth=1;
     updatesmoothblocks=1;
     break;
   case SMOOTH_ATSTART:
     sb_atstart=1-sb_atstart;
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
  updatemenu=1;  
 // updatefaces=1;
  updatefacelists=1;
  updatehiddenfaces=1;
  GLUTPOSTREDISPLAY
}

/* ------------------ sv_BlockageMenu ------------------------ */

void svWINAPI sv_BlockageMenu(int value){
  BlockageMenu(value);
}


/* ------------------ RotateTypeMenu ------------------------ */

void RotateTypeMenu(int value){
  eyeview = value;
  updatemenu=1;  
  GLUTPOSTREDISPLAY
}

/* ------------------ TitleMenu ------------------------ */

void TitleMenu(int value){
  updatemenu=1;  
  GLUTPOSTREDISPLAY
  ASSERTFLAG(visTitle0);
  ASSERTFLAG(visTitle1);
  ASSERTFLAG(visTitle2);
  switch (value){
  case 0:
    visTitle0 = 1 - visTitle0;
    break;
  case 1:
    visTitle1 = 1 - visTitle1;
    break;
  case 2:
    visTitle2 = 1 - visTitle2;
    break;
  case 3:
    visTitle0=1;
    visTitle1=1;
    visTitle2=1;
    break;
  case 4:
    visTitle0=0;
    visTitle1=0;
    visTitle2=0;
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
        part5class *partclassi;

        partclassi = partclassinfo + i;
        update_partclass_depend(partclassi);

      }

      GLUTPOSTREDISPLAY
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
  else if(value==-1){
    for(i=0;i<nobject_defs;i++){
      objecti = object_defs[i];
      objecti->visible=1;
    }
  }
  else if(value==-2){
    for(i=0;i<nobject_defs;i++){
      objecti = object_defs[i];
      objecti->visible=0;
    }
  }
  else if(value==-3){
    select_device=1-select_device;
  }
  else if(value==-999){
  }
  else{
    int oldval;

    oldval=device_sphere_segments;
    value=-value;
    device_sphere_segments=value;
    if(oldval!=device_sphere_segments){
      initspheresegs(device_sphere_segments,2*device_sphere_segments);
    }
  }
  updatemenu=1;
  GLUTPOSTREDISPLAY
}

/* ------------------ ZoneShowMenu ------------------------ */

void ZoneShowMenu(int value){
  switch (value){
  case 999:
    return;
    break;
  case 1:
    visVZone=0;
    visHZone=1;
    break;
  case 2:
    visVZone=1;
    visHZone=0;
    break;
  case 3:
    visVZone=1; 
    visHZone=1;
    break;
  case 4:
    visVZone=0; 
    visHZone=0;
    break;
  case 5:
  case 6:
    sethazardcolor=1-sethazardcolor;
    break;
  case 11:         //solid
    visVentSolid=1;
    visVentLines=0;
    visVents=1;
    break;
  case 12:         //lines
    visVentLines=1;
    visVentSolid=0;
    visVents=1;
    break;
  case 13:         //hide
    visVentLines=0;
    visVentSolid=0;
    break;
  case 14:
    visVents=1-visVents;
    if(visVents==0){
      visVentLines=0;
      visVentSolid=0;
    }
    if(visVents==1){
      visVentSolid=1;
      visVentLines=0;
    }
    break;
  default:
    ASSERT(FFALSE);
  }
  updatemenu=1;  
  GLUTPOSTREDISPLAY
}

/* ------------------ GeometryMenu ------------------------ */

void GeometryMenu(int value){

  switch (value){
  case 3:
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
    //shaded 17 0
    //stepped 18 1
    //line    19 2
    //texture 20 3
    //hidden 21 4
  case 17:
  case 18:
  case 19:
  case 20:
  case 21:
    if(value==21){
      BlockageMenu(visBlocksSave);
    }
    else{
      visBlocksSave=visBlocks;
      BlockageMenu(visBLOCKHide);
    }
    visTerrainType=value-17;
    if(visTerrainType==0){
      planar_terrain_slice=0;
    }
    else{
      planar_terrain_slice=1;
    }
    update_glui_wui();
    break;
  case 11:
    if(isZoneFireModel)visFrame=1;
    /*
    visFloor=1;
    visWalls=1;
    visCeiling=1;
    */
    visVents=1;
    BlockageMenu(visBLOCKAsInput);
    break;
  case 13:
    visFrame=0;
    visFloor=0;
    visWalls=0;
    visCeiling=0;
    visVents=0;
    visGrid=0;
    BlockageMenu(visBLOCKHide);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  updatefacelists=1;
  updatemenu=1;  
  GLUTPOSTREDISPLAY
}

/* ------------------ MENU_vslice ------------------------ */

void MENU_vslice(int vec_type){
  int ii,i;
  slice *sd;
  vslice *vd;
  char check[]="*";
  char *menulabel;
  char a_menulabel[1000];
  menulabel=a_menulabel;

  for(ii=0;ii<nvslice;ii++){
    i=vsliceorderindex[ii];
    vd = vsliceinfo + i;
    if(vd->vec_type!=vec_type)continue;
    sd = sliceinfo + vd->ival;
    if(vd->loaded==1){
      STRCPY(menulabel,check);
      STRCAT(menulabel,sd->menulabel);
    }
    else{
      STRCPY(menulabel,sd->menulabel);
    }
    if(sd->vec_comp==0||showallslicevectors==1)glutAddMenuEntry(menulabel,i);
  }
}

/* ------------------ InitMenus ------------------------ */

void InitMenus(int unload){
  int showflag, hideflag;
  int n,i,ii;
  FILE *stream;
  FILE *stream2;
  FILE *stream3;
  int readinifile=0;
  char caselabel[255],chari[4];
  int showall,hideall;
  char levellabel[1024];
  int nsliceloaded;
  int nevacloaded;
  int nsmoke3dloaded;
  int nvsliceloaded2,nvsliceloaded;
  int nmultisliceloaded;
  int npartloaded;
  int npart5loaded;
  int npart4loaded;
  int npatchloaded;
  int nplot3dloaded;
  int nisoloaded;
  int nvslice0, nvslice1, nvslice2;
  int nvsliceloaded0;
  int nvsliceloaded1;
  char check[]="*";
  multislice *mslicei;
  multivslice *mvslicei;
  mesh *cmesh;
  char a_menulabel[1000];
  char *menulabel;
  int j;
  int ntextures_used;
  int multiprop;

static int titlemenu=0, labelmenu=0, colorbarmenu=0, lightingmenu=0, showhidemenu=0;
static int optionmenu=0, rotatetypemenu=0;
static int resetmenu=0, frameratemenu=0, rendermenu=0, smokeviewinimenu=0, inisubmenu=0;
#ifdef pp_COMPRESS
static int compressmenu=0;
#endif
static int showhideslicemenu=0,showvslicemenu=0;
static int plot3dshowmenu=0, staticvariablemenu=0, helpmenu=0;
static int vectorskipmenu=0,unitsmenu=0;
static int isosurfacemenu=0, isovariablemenu=0, levelmenu=0;
static int isoblockmenu=0, fontmenu=0, aperturemenu=0,dialogmenu=0,zoommenu=0;
static int gridslicemenu=0, blockagemenu=0, loadpatchmenu=0, ventmenu=0;
static int loadisomenu=0, isosurfacetypemenu=0;
static int geometrymenu=0, loadunloadmenu=0, reloadmenu=0, disclaimermenu=0, terrain_showmenu=0;
static int scriptmenu=0;
static int loadplot3dmenu=0, unloadvslicemenu=0, unloadslicemenu=0;
static int loadterrainmenu=0, unloadterrainmenu=0;
static int loadsmoke3dmenu=0,loadsmoke3dsootmenu=0,loadsmoke3dhrrmenu=0,loadsmoke3dwatermenu=0;
static int unloadsmoke3dmenu=0;
static int unloadevacmenu=0, unloadpartmenu=0, loadslicemenu=0, loadmultislicemenu=0;
static int *loadsubvslicemenu=NULL, nloadsubvslicemenu=0;
static int *loadsubslicemenu=NULL, nloadsubslicemenu=0, iloadsubslicemenu=0;
static int *loadsubmslicemenu=NULL, nloadsubmslicemenu=0;
static int *loadsubmvslicemenu=NULL, nloadsubmvslicemenu=0;
static int *loadsubplot3dmenu=NULL, nloadsubplot3dmenu=0;
static int loadmultivslicemenu=0, unloadmultivslicemenu=0;
static int unloadmultislicemenu=0, vslicemenu=0, staticslicemenu=0;
static int evacmenu=0, particlemenu=0, particlesubmenu=0, showpatchmenu=0, zonemenu=0, isoshowmenu=0, isolevelmenu=0, smoke3dshowmenu=0;
static int particle5showmenu=0;
static int particlepropshowmenu=0,humanpropshowmenu=0;
static int *particlepropshowsubmenu=NULL;
static int particlestreakshowmenu=0;
static int tourmenu=0;
static int showtourmenu=0, avatartourmenu=0,avatarevacmenu=0;
static int trainerviewmenu=0,mainmenu=0,zoneshowmenu=0,particleshowmenu=0,evacshowmenu=0,targetmenu=0;
static int showobjectsmenu=0,spheresegmentmenu=0,propmenu=0;
static int unloadplot3dmenu=0, unloadpatchmenu=0, unloadisomenu=0;
static int showmultislicemenu=0;
static int textureshowmenu=0;
#ifdef _DEBUG
static int menu_count=0;
static int in_menu=0;
#endif


#ifdef _DEBUG
  printf("Updating Menus %i In menu %i\n",menu_count++,in_menu);
  in_menu=1;
#endif
  update_showhidebuttons();
  GLUTPOSTREDISPLAY
  cmesh=current_mesh;

  menulabel = a_menulabel;
  nsmoke3dloaded=0;
  nsliceloaded=0;
  nvsliceloaded=0;
  for(i=0;i<nslice_files;i++){
    slice *sd;

    sd = sliceinfo + i;
    if(sd->loaded==1)nsliceloaded++;
  }
  for(i=0;i<nsmoke3d_files;i++){
    smoke3d *smoke3di;

    smoke3di=smoke3dinfo+i;
    if(smoke3di->loaded==1)nsmoke3dloaded++;
  }
  for(i=0;i<nvslice;i++){
    vslice *vd;

    vd = vsliceinfo + i;
    if(vd->loaded==1)nvsliceloaded++;
  }

  for(i=0;i<nmultislices;i++){
    mslicei = multisliceinfo + i;
    mslicei->loaded=0;
    mslicei->display=0;
    for(j=0;j<mslicei->nslices;j++){
      slice *sd;

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
  for(i=0;i<nmultivslices;i++){
    mvslicei = multivsliceinfo + i;
    mvslicei->loaded=0;
    mvslicei->display=0;
    for(j=0;j<mvslicei->nvslices;j++){
      vslice *vd;

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
  npart4loaded=0;
  npartloaded=0;
  nevacloaded=0;
  for(i=0;i<npart_files;i++){
    particle *parti;

    parti = partinfo+i;
    if(parti->loaded==1&&parti->evac==0)npartloaded++;
    if(parti->loaded==1&&parti->evac==1)nevacloaded++;
    if(parti->loaded==1){
      if(parti->version==1)npart5loaded++;
      if(parti->version==0)npart4loaded++;
    }
  }

  nplot3dloaded=0;
  for(i=0;i<nplot3d_files;i++){
    plot3d *plot3di;

    plot3di = plot3dinfo + i;
    if(plot3di->loaded==1)nplot3dloaded++;
  }

  nisoloaded=0;
  for(i=0;i<niso_files;i++){
    iso *isoi;

    isoi = isoinfo + i;
    if(isoi->loaded==1)nisoloaded++;
  }

  npatchloaded=0;
  for(i=0;i<npatch_files;i++){
    patch *patchi;

    patchi = patchinfo + i;
    if(patchi->loaded==1)npatchloaded++;
  }

#define DESTROYMENU(f) 
#define CREATEMENU(menu,Menu) menu=glutCreateMenu(Menu);\
  if(nmenus<10000){\
    strcpy(menuinfo[nmenus].label,#Menu);\
    menuinfo[nmenus++].menuvar=menu;\
  }

  {
    for(i=0;i<nmenus;i++){
      int menuvar;
      menudata *menui;

      menui = menuinfo + i;

      menuvar = menui->menuvar;
      if(menuvar!=0){
        glutDestroyMenu(menuvar);
        menuvar=0;
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
  if(npatch_files>0){
    CREATEMENU(showpatchmenu,ShowPatchMenu);
    npatchloaded=0;
    {
      int local_do_threshold=0;

      for(ii=0;ii<npatch_files;ii++){
        patch *patchi;

        i = patchorderindex[ii];
        patchi = patchinfo+i;
        if(patchi->loaded==0)continue;
        npatchloaded++;
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
      if(activate_threshold==1&&local_do_threshold==1){
        glutAddMenuEntry("-",DUMMYwallmenu);
        if(vis_threshold==1)glutAddMenuEntry("*char",SHOW_CHAR);
        if(vis_threshold==0)glutAddMenuEntry("char",SHOW_CHAR);
      }

    }
    if(npatchloaded>0){
      glutAddMenuEntry("-",DUMMYwallmenu);
    }
    if(npatchloaded>1){
      glutAddMenuEntry("Show All Boundary Files",SHOWALL_BOUNDARY);
      glutAddMenuEntry("Hide All Boundary Files",HIDEALL_BOUNDARY);
      glutAddMenuEntry("-",DUMMYwallmenu);
    }
    if(showexterior==1){
      glutAddMenuEntry("*Exterior",EXTERIORwallmenu);
    }
    if(showexterior==0){
      glutAddMenuEntry("Exterior",EXTERIORwallmenu);
    }
    if(visPatchType[INTERIORwall]==1){
      glutAddMenuEntry("*Interior",INTERIORwallmenu);
    }
    if(visPatchType[INTERIORwall]==0){
      glutAddMenuEntry("Interior",INTERIORwallmenu);
    }
    if(ispatchtype(FRONTwall)==1&&visPatchType[FRONTwall]==1){
      glutAddMenuEntry("*Front",FRONTwallmenu);
    }
    if(ispatchtype(FRONTwall)==1&&visPatchType[FRONTwall]==0){
      glutAddMenuEntry("Front",FRONTwallmenu);
    }
    if(ispatchtype(BACKwall)==1&& visPatchType[BACKwall]==1){
      glutAddMenuEntry("*Back",BACKwallmenu);
    }
    if(ispatchtype(BACKwall)==1&& visPatchType[BACKwall]==0){
      glutAddMenuEntry("Back",BACKwallmenu);
    }
    if(ispatchtype(LEFTwall)==1&& visPatchType[LEFTwall]==1){ 
      glutAddMenuEntry("*Left",LEFTwallmenu);
    }
    if(ispatchtype(LEFTwall)==1&& visPatchType[LEFTwall]==0){
      glutAddMenuEntry("Left",LEFTwallmenu);
    }
    if(ispatchtype(RIGHTwall)==1&&visPatchType[RIGHTwall]==1){
      glutAddMenuEntry("*Right",RIGHTwallmenu);
    }
    if(ispatchtype(RIGHTwall)==1&&visPatchType[RIGHTwall]==0){
      glutAddMenuEntry("Right",RIGHTwallmenu);
    }
    if(ispatchtype(UPwall)==1&&   visPatchType[UPwall]==1){
      glutAddMenuEntry("*Up",UPwallmenu);
    }
    if(ispatchtype(UPwall)==1&&   visPatchType[UPwall]==0){
      glutAddMenuEntry("Up",UPwallmenu);
    }
    if(ispatchtype(DOWNwall)==1&& visPatchType[DOWNwall]==1){
      glutAddMenuEntry("*Down",DOWNwallmenu);
    }
    if(ispatchtype(DOWNwall)==1&& visPatchType[DOWNwall]==0){
      glutAddMenuEntry("Down",DOWNwallmenu);
    }
  }

/* --------------------------------blockage menu -------------------------- */
  CREATEMENU(blockagemenu,BlockageMenu);
  if(use_menusmooth==1){
    if(sb_atstart==1){
      glutAddMenuEntry("*Smooth blockages at startup",SMOOTH_ATSTART);
    }
    else{
      glutAddMenuEntry("Smooth blockages at startup",SMOOTH_ATSTART);
    }
    glutAddMenuEntry("Smooth Blockages Now",SMOOTH_BLOCKAGES);
  }
  glutAddMenuEntry("View Method:",999);
  if(visBlocks==visBLOCKAsInput){
    glutAddMenuEntry("   *Defined In Input File",visBLOCKAsInput);
  }
   else{
    glutAddMenuEntry("   Defined In Input File",visBLOCKAsInput);
  }
  if(use_menusmooth==1){
    if(smooth_block_solid==1){
      glutAddMenuEntry("       *Smooth blockages drawn opaque",visSmoothBLOCKSolid);
    }
    else{
      glutAddMenuEntry("       Smooth blockages drawn opaque",visSmoothBLOCKSolid);
    }
  }
  if(visBlocks==visBLOCKNormal){
    glutAddMenuEntry("   *Solid",visBLOCKNormal);
    if(nsmoothblocks>0){
      if(visSmoothAsNormal==1){
         glutAddMenuEntry("      Smooth",visBLOCKSmoothAsNormal);
      }
      else{
         glutAddMenuEntry("      *Smooth",visBLOCKSmoothAsNormal);
      }
    }
    if(ntransparentblocks>0){
      if(visTransparentBlockage==1){
         glutAddMenuEntry("      *Transparent",visBLOCKTransparent);
      }
      else{
         glutAddMenuEntry("      Transparent",visBLOCKTransparent);
      }
    }
  }
  else{
    glutAddMenuEntry("   Solid",visBLOCKNormal);
  }
  if(visBlocks==visBLOCKOutline){
    glutAddMenuEntry("   *Outline",visBLOCKOutline);
  }
   else{
     glutAddMenuEntry("   Outline",visBLOCKOutline);
   }
  if(visBlocks==visBLOCKHide){
    glutAddMenuEntry("   *Hidden",visBLOCKHide);
  }
  else{
    glutAddMenuEntry("   Hidden",visBLOCKHide);
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

      glutAddMenuEntry("-",999);
      glutAddMenuEntry("Show/Hide Blockage Types:",999);
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
  glutAddMenuEntry("-",999);
  glutAddMenuEntry("Locations:",999);
  if(blocklocation==BLOCKlocation_grid){
    glutAddMenuEntry("   *Actual",BLOCKlocation_grid);
  }
  else{
    glutAddMenuEntry("   Actual",BLOCKlocation_grid);
  }
  if(blocklocation==BLOCKlocation_exact){
    glutAddMenuEntry("   *Requested",BLOCKlocation_exact);
  }
  else{
    glutAddMenuEntry("   Requested",BLOCKlocation_exact);
  }
  if(ncadgeom>0){
    if(blocklocation==BLOCKlocation_cad){
      glutAddMenuEntry("   *Cad",BLOCKlocation_cad);
    }
    else{
      glutAddMenuEntry("   Cad",BLOCKlocation_cad);
    }
    {
      cadgeom *cd;
      cadlook *cdi;
      int showtexturemenu;

      showtexturemenu=0;
      for(i=0;i<ncadgeom;i++){
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
          glutAddMenuEntry(" *Show Cad Textures",BLOCKtexture_cad);
        }
        else{
          glutAddMenuEntry(" Show Cad Textures",BLOCKtexture_cad);
        }
      }
    }
  }


/* --------------------------------level menu -------------------------- */

  if(nplot3d_files>0){
    CREATEMENU(levelmenu,LevelMenu);
    for(i=1;i<nrgb-1;i++){
      if(colorlabeliso!=NULL){
        char *colorlabel;
        char levellabel[256];

        colorlabel=&colorlabeliso[plotn-1][nrgb-2-i][0];
        strcpy(levellabel,"");
        if(plotiso[plotn-1]==nrgb-2-i&&visiso==1){
          strcat(levellabel,"*");
        }
        strcat(levellabel,colorlabel);
        glutAddMenuEntry(levellabel,nrgb-2-i);
      }
      else{
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

  if(nplot3d_files>0){
    CREATEMENU(staticvariablemenu,StaticVariableMenu);
    for(n=0;n<numplot3dvars;n++){
      if(plotn-1==n){
        STRCPY(menulabel,"*");
        STRCAT(menulabel,shortp3label[n]);
        glutAddMenuEntry(menulabel,n+1);
      }
      else{
        glutAddMenuEntry(shortp3label[n],n+1);
      }
    }
  }

/* --------------------------------iso block menu -------------------------- */

  if(nplot3d_files>0){
    CREATEMENU(isoblockmenu,IsoBlockMenu);
    for(i=1;i<=offsetmax;i++){
      if(isooffset==i){
        sprintf(chari,"*%i",i);
      }
      else{
        sprintf(chari,"%i",i);
      }
      glutAddMenuEntry(chari,i);
    }
  }

/* --------------------------------iso variable menu -------------------------- */

  if(nplot3d_files>0){
    CREATEMENU(isovariablemenu,IsoVariableMenu);
    for(n=0;n<numplot3dvars;n++){
      if(plotn-1==n&&visiso==1){
        STRCPY(menulabel,"*");
        STRCAT(menulabel,shortp3label[n]);
        glutAddMenuEntry(menulabel,n+1);
      }
      else{
        glutAddMenuEntry(shortp3label[n],n+1);
      }
    }
  }

/* --------------------------------iso surface menu -------------------------- */
  if(nplot3d_files>0){
    CREATEMENU(isosurfacetypemenu,IsoSurfaceTypeMenu);
    if(p3dsurfacesmooth==1&&p3dsurfacetype==1){
      glutAddMenuEntry("*Smooth",0);
    }
     else{
       glutAddMenuEntry("Smooth",0);
     }
     if(p3dsurfacesmooth==0&&p3dsurfacetype==1){
       glutAddMenuEntry("*Facets",1);
     }
    else{
      glutAddMenuEntry("Facets",1);
    }
    if(p3dsurfacetype==2)glutAddMenuEntry("*Triangles",2);
    if(p3dsurfacetype!=2)glutAddMenuEntry("Triangles",2);
    if(p3dsurfacetype==3)glutAddMenuEntry("*Points",3);
    if(p3dsurfacetype!=3)glutAddMenuEntry("Points",3);

    CREATEMENU(isosurfacemenu,IsoSurfaceMenu);
    glutAddSubMenu("Solution Variable",isovariablemenu);
    glutAddSubMenu("Solution Value",levelmenu);
    glutAddSubMenu("Block Size",isoblockmenu);
    glutAddSubMenu("Surface Type",isosurfacetypemenu);
    glutAddMenuEntry("Hide",1);
  }

/* --------------------------------vector skip menu -------------------------- */

  if(nplot3d_files>0){
    CREATEMENU(vectorskipmenu,VectorSkipMenu);
    if(visVector==1)glutAddMenuEntry("*Show",-2);
    if(visVector!=1)glutAddMenuEntry("Show",-2);
    glutAddMenuEntry("Frequency:",-1);
    if(vectorskip==1)glutAddMenuEntry("*All",1);
    if(vectorskip!=1)glutAddMenuEntry("All",1);
    if(vectorskip==2)glutAddMenuEntry("*Every 2nd",2);
    if(vectorskip!=2)glutAddMenuEntry("Every 2nd",2);
    if(vectorskip==3)glutAddMenuEntry("*Every 3rd",3);
    if(vectorskip!=3)glutAddMenuEntry("Every 3rd",3);
    if(vectorskip==4)glutAddMenuEntry("*Every 4th",4);
    if(vectorskip!=4)glutAddMenuEntry("Every 4th",4);
  }

  if(ntextures_loaded_used>0){
    CREATEMENU(textureshowmenu,TextureShowMenu);
    ntextures_used=0;
    for(i=0;i<ntextures;i++){
      texture *texti;

      texti = textureinfo + i;
      if(texti->loaded==0||texti->used==0)continue;
      ntextures_used++;
      if(texti->display==1){
        STRCPY(menulabel,check);
        STRCAT(menulabel,texti->file);  
        glutAddMenuEntry(menulabel,i);
      }
      else{
        STRCPY(menulabel,texti->file);
        glutAddMenuEntry(menulabel,i);
      }
    }
    if(ntextures_used>1){
      glutAddMenuEntry("-",-999);
      glutAddMenuEntry("Show All",-1);
      glutAddMenuEntry("Hide All",-2);
    }
  }

/* --------------------------------static slice menu -------------------------- */
  if(nplot3d_files>0){
    CREATEMENU(staticslicemenu,Plot3DShowMenu);
    glutAddSubMenu("Solution Variable",staticvariablemenu);
    if(cmesh->visz==1)glutAddMenuEntry("*xy plane",1);
    if(cmesh->visz==0)glutAddMenuEntry("xy plane",1);
    if(cmesh->visy==1)glutAddMenuEntry("*xz  plane",2);
    if(cmesh->visy==0)glutAddMenuEntry("xz  plane",2);
    if(cmesh->visx==1)glutAddMenuEntry("*yz  plane",3);
    if(cmesh->visx==0)glutAddMenuEntry("yz  plane",3);
    if(vectorspresent==1)glutAddSubMenu("Flow vectors",vectorskipmenu);
    if(p3cont2d==SHADED_CONTOURS){
      glutAddMenuEntry("*Continuous Contours",4);
    }
    if(p3cont2d!=SHADED_CONTOURS){
      glutAddMenuEntry("Continuous Contours",4);
    }
    glutAddMenuEntry("Show all planes in active mesh",5);
    glutAddMenuEntry("Hide all planes in active mesh",6);

    CREATEMENU(plot3dshowmenu,Plot3DShowMenu);
    if(nplot3dloaded>0){
      for(ii=0;ii<nplot3d_files;ii++){
        plot3d *plot3di;

        i=plot3dorderindex[ii];
        plot3di = plot3dinfo + i;
        if(ii==0){
          glutAddMenuEntry(plot3di->longlabel,997);
        }
        if(ii!=0&&strcmp(plot3di->longlabel,plot3dinfo[plot3dorderindex[ii-1]].longlabel)!=0){
          glutAddMenuEntry(plot3di->longlabel,997);
        }
        if(plot3di->loaded==0)continue;
        if(plotstate==STATIC_PLOTS&&plot3di->display==1){
          STRCPY(menulabel,check);
          STRCAT(menulabel,plot3di->menulabel);  
        }
        else{
          STRCPY(menulabel,plot3di->menulabel);
        }
        glutAddMenuEntry(menulabel,1000+i);
      }
      if(nplot3dloaded>1){
        glutAddMenuEntry("-",997);
        glutAddMenuEntry("Show All PLOT3D Files",SHOWALL_PLOT3D);
        glutAddMenuEntry("Hide All PLOT3D Files",HIDEALL_PLOT3D);
      }
      glutAddMenuEntry("-",997);
    }
    glutAddSubMenu("2D Contours",staticslicemenu);
    if(unload_qdata==0){
      glutAddSubMenu("3D Contours",isosurfacemenu);
    }

  }

/* --------------------------------grid slice menu -------------------------- */

  CREATEMENU(gridslicemenu,GridSliceMenu);
  if(cmesh->visz==1){
    if(visGrid==1){
      glutAddMenuEntry("*xy plane",1);
    }
    else{
      glutAddMenuEntry("*xy plane (grid turned off)",1);
    }
  }
  else{
    glutAddMenuEntry("xy plane",1);
  }
  if(cmesh->visy==1){
    if(visGrid==1){
      glutAddMenuEntry("*xz plane",2);
    }
    else{
      glutAddMenuEntry("*xz plane (grid turned off)",2);
    }
  }
  else{
    glutAddMenuEntry("xz plane",2);
  }
  if(cmesh->visx==1){
    if(visGrid==1){
      glutAddMenuEntry("*yz plane",3);
    }
    else{
      glutAddMenuEntry("*yz plane (grid turned off)",3);
    }
  }
  else{
    glutAddMenuEntry("yz plane",3);
  }
  glutAddMenuEntry("Show All",4);
  if(visGrid==1)glutAddMenuEntry("Hide All",5);

/* --------------------------------vent menu -------------------------- */

  CREATEMENU(ventmenu,VentMenu);
  {
    int ntotal_vents=0;

    for(i=0;i<nmeshes;i++){
      mesh *meshi;

      meshi=meshinfo+i;
      ntotal_vents+=meshi->nvents;
    }
    if(ntotal_vents>0){
      if(visVents==1)glutAddMenuEntry("*All",10);
      if(visVents==0)glutAddMenuEntry("All",10);
      if(nopenvents>0){
        if(visOpenVents==1)glutAddMenuEntry("  *Open Vents",14);
        if(visOpenVents==0)glutAddMenuEntry("  Open Vents",14);
      }
      if(ndummyvents>0){
        if(visDummyVents==1)glutAddMenuEntry("  *Dummy Vents",16);
        if(visDummyVents==0)glutAddMenuEntry("  Dummy Vents",16);
      }
      glutAddMenuEntry("-",-1);
      if(nopenvents_nonoutline>0){
        if(visOpenVentsAsOutline==1)glutAddMenuEntry("*Open Vents As Outlines",15);
        if(visOpenVentsAsOutline==0)glutAddMenuEntry("Open Vents As Outlines",15);
      }
      if(have_vents_int==1){
        if(show_bothsides_int==1)glutAddMenuEntry("*Two sided (interior)",18);
        if(show_bothsides_int==0)glutAddMenuEntry("Two sided (interior)",18);
      }
      if(show_bothsides_ext==1)glutAddMenuEntry("*Two sided (exterior)",19);
      if(show_bothsides_ext==0)glutAddMenuEntry("Two sided (exterior)",19);
      if(nvent_transparent>0){
        if(show_transparent_vents==1)glutAddMenuEntry("*Transparent",20);
        if(show_transparent_vents==0)glutAddMenuEntry("Transparent",20);
      }
    }
  }
  
/* --------------------------------terrain_showmenu -------------------------- */

  CREATEMENU(terrain_showmenu,GeometryMenu);

  if(visTerrainType==0)glutAddMenuEntry("*3D surface",17);
  if(visTerrainType!=0)glutAddMenuEntry("3D surface",17);
  if(visTerrainType==1)glutAddMenuEntry("*2D stepped",18);
  if(visTerrainType!=1)glutAddMenuEntry("2D stepped",18);
  if(visTerrainType==2)glutAddMenuEntry("*2D lines",19);
  if(visTerrainType!=2)glutAddMenuEntry("2D lines",19);
  if(terrain_texture!=NULL&&terrain_texture->loaded==1){
    if(visTerrainType==3)glutAddMenuEntry("*Image",20);
    if(visTerrainType!=3)glutAddMenuEntry("Image",20);
  }
  if(visTerrainType==4)glutAddMenuEntry("*Hidden",21);
  if(visTerrainType!=4)glutAddMenuEntry("Hidden",21);
    
  if(nobject_defs>0){
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
  }

  if(nobject_defs>0){
    CREATEMENU(showobjectsmenu,ShowObjectsMenu);
    for(i=0;i<nobject_defs;i++){
      sv_object *obj_typei;

      obj_typei = object_defs[i];
      obj_typei->used_by_device=0;
    }
    for(i=0;i<ndeviceinfo;i++){
      device *devicei;
      propdata *propi;
      int jj;

      devicei = deviceinfo + i;
      propi = devicei->prop;
      for(jj=0;jj<propi->nsmokeview_ids;jj++){
        sv_object *objectj;

        objectj = propi->smv_objects[jj];
        objectj->used_by_device=1;
      }
    }
    for(i=0;i<npart5prop;i++){
      part5prop *partpropi;

      partpropi = part5propinfo + i;
      for(j=0;j<npartclassinfo;j++){
        part5class *partclassj;
        propdata *propi;
        int jj;

        if(partpropi->class_present[j]==0)continue;
        partclassj = partclassinfo + j;
        propi = partclassj->prop;
        if(propi==NULL)continue;
        for(jj=0;jj<propi->nsmokeview_ids;jj++){
          sv_object *objectj;

          objectj = propi->smv_objects[jj];
          objectj->used_by_device=1;
        }
      }
    }

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
    glutAddMenuEntry("-",-999);
    if(ndeviceinfo>0){
      if(select_device==1){
        glutAddMenuEntry("*Select",-3);
      }
      else{
        glutAddMenuEntry("Select",-3);
      }
    }
    glutAddMenuEntry("Show All",-1);
    glutAddMenuEntry("Hide All",-2);
    glutAddMenuEntry("-",-999);
    glutAddSubMenu("Segments",spheresegmentmenu);

  }

/* --------------------------------geometry menu -------------------------- */

  CREATEMENU(geometrymenu,GeometryMenu);
  if(showedit==0&&ntotal_blockages>0)glutAddSubMenu("Blockages",blockagemenu);
  if(nobject_defs>0){
    int num_activedevices=0;

    for(i=0;i<nobject_defs;i++){
      sv_object *obj_typei;

      obj_typei = object_defs[i];
      if(obj_typei->used==1)num_activedevices++;
    }
    if(isZoneFireModel==0||(isZoneFireModel==1&&num_activedevices>1)){
      glutAddSubMenu("Objects",showobjectsmenu);
    }
  }
    //shaded 17 0
    //stepped 18 1
    //line    19 2
    //texture 20 3
    //hidden 21 4
  if(nterraininfo>0){
    glutAddSubMenu("Terrain",terrain_showmenu);
  }
 glutAddSubMenu("Vents",ventmenu);
 if(ntotal_blockages>0||isZoneFireModel==0){
    glutAddSubMenu("Grid",gridslicemenu);
  }
  if(isZoneFireModel==0){
    if(visFrame==1)glutAddMenuEntry("*Outline",3);
    if(visFrame==0)glutAddMenuEntry("Outline",3);
  }
  else{
    visFrame=0;
  }
  glutAddMenuEntry("Show All",11);
  glutAddMenuEntry("Hide All",13);

/* --------------------------------title menu -------------------------- */

  ntitles=1;
  if(strlen(TITLE1)!=0)ntitles++;
  if(strlen(TITLE2)!=0)ntitles++;
  if(ntitles>1){
    CREATEMENU(titlemenu,TitleMenu);
    if(visTitle0==1)glutAddMenuEntry("*Default",0);
    if(visTitle0==0)glutAddMenuEntry("Default",0);
    if(strlen(TITLE1)!=0){
      if(visTitle1==1){glutAddMenuEntry("*User Title 1",1);}
      if(visTitle1==0){glutAddMenuEntry("User Title 1",1);}
    }
    if(strlen(TITLE2)!=0){
      if(visTitle2==1){glutAddMenuEntry("*User Title 2",2);}
      if(visTitle2==0){glutAddMenuEntry("User Title 2",2);}
    }
    glutAddMenuEntry("Show All",3);
    glutAddMenuEntry("Hide All",4);
  }

/* --------------------------------label menu -------------------------- */

  CREATEMENU(labelmenu,LabelMenu);
  if(visColorLabels==1)glutAddMenuEntry("*Color Bars",0);
  if(visColorLabels==0)glutAddMenuEntry("Color Bars",0);
  if(visTimeLabels==1)glutAddMenuEntry("*Time Bars",1);
  if(visTimeLabels==0)glutAddMenuEntry("Time Bars",1);
  if(nticks>0){
    if(visTicks==0)glutAddMenuEntry("FDS Generated Ticks",12);
    if(visTicks==1)glutAddMenuEntry("*FDS Generated Ticks",12);
  }
  if(vis_user_ticks==1)glutAddMenuEntry("*User Settable Ticks",18);
  if(vis_user_ticks==0)glutAddMenuEntry("User Settable Ticks",18);

  if(ntitles>1){
    glutAddSubMenu("Titles",titlemenu);
  }
  else{
    if(visTitle==1)glutAddMenuEntry("*Title",2);
    if(visTitle==0)glutAddMenuEntry("Title",2);
  }
  if(visaxislabels==1)glutAddMenuEntry("*Axis",6);
  if(visaxislabels==0)glutAddMenuEntry("Axis",6);
  if(visFramerate==1)glutAddMenuEntry("*Frame Rate",3);
  if(visFramerate==0)glutAddMenuEntry("Frame Rate",3);
  if(visTimelabel==1)glutAddMenuEntry("*Time label",8);
  if(visTimelabel==0)glutAddMenuEntry("Time label",8);
  if(visFramelabel==1)glutAddMenuEntry("*Frame label",9);
  if(visFramelabel==0)glutAddMenuEntry("Frame label",9);
  if(hrrinfo!=NULL){
    if(visHRRlabel==1)glutAddMenuEntry("*HRR label",16);
    if(visHRRlabel==0)glutAddMenuEntry("HRR label",16);
  }
  if(show_hrrcutoff_active==1){
    if(show_hrrcutoff==1)glutAddMenuEntry("*HRRPUV cutoff",17);
    if(show_hrrcutoff==0)glutAddMenuEntry("HRRPUV cutoff",17);
  }
  if(vis_slice_average==1)glutAddMenuEntry("*Slice Average",15);
  if(vis_slice_average==0)glutAddMenuEntry("Slice Average",15);
  if(nmeshes>1){
    if(visBlocklabel==1)glutAddMenuEntry("*Mesh label",10);
    if(visBlocklabel==0)glutAddMenuEntry("Mesh label",10);
  }
#ifdef pp_memstatus
  if(visAvailmemory==1)glutAddMenuEntry("*Memory Load",11);
  if(visAvailmemory==0)glutAddMenuEntry("Memory Load",11);
#endif
  if(nlabels>0){
    if(visLabels==1)glutAddMenuEntry("*Text Labels",7);
    if(visLabels==0)glutAddMenuEntry("Text Labels",7);
  }
  if(ntotal_blockages>0||isZoneFireModel==0){
    if(visgridloc==1)glutAddMenuEntry("*Grid Locations",14);
    if(visgridloc==0)glutAddMenuEntry("Grid Locations",14);
  }

  glutAddMenuEntry("Show All",4);
  glutAddMenuEntry("Hide All",5);

/* --------------------------------rotate type menu -------------------------- */

  CREATEMENU(rotatetypemenu,RotateTypeMenu);
  switch (eyeview){
  case EYE_CENTERED:
    glutAddMenuEntry("World Centered",0);
    glutAddMenuEntry("*Eye Centered",1);
    glutAddMenuEntry("World Centered, level rotation",2);
  break;
  case WORLD_CENTERED:
    glutAddMenuEntry("*World Centered",0);
    glutAddMenuEntry("Eye Centered",1);
    glutAddMenuEntry("World Centered, level rotation",2);
    break;
  case WORLD_CENTERED_LEVEL:
    glutAddMenuEntry("World Centered",0);
    glutAddMenuEntry("Eye Centered",1);
    glutAddMenuEntry("*World Centered, level rotation",2);
  break;
  default:
    ASSERT(FFALSE);
    break;
  }

/* --------------------------------zone show menu -------------------------- */

  if(nzone>0){
    CREATEMENU(zoneshowmenu,ZoneShowMenu);
    glutAddMenuEntry("Layers",999);
    if(ReadZoneFile==1){
      if(sethazardcolor==1){
        glutAddMenuEntry("  Colors: *Hazard",5);
        glutAddMenuEntry("  Colors: Temperature",6);
      }
      else{
        glutAddMenuEntry("  Colors: Hazard",5);
        glutAddMenuEntry("  Colors: *Temperature",6);
      }
      if(visHZone==1)glutAddMenuEntry("  *Horizontal",1);
      if(visHZone==0)glutAddMenuEntry("  Horizontal",1);
      if(visVZone==1)glutAddMenuEntry("  *Vertical",2);
      if(visVZone==0)glutAddMenuEntry("  Vertical",2);
      if(visHZone==0&&visVZone==0){
        glutAddMenuEntry("  *Hide",4);
      }
      else{
        glutAddMenuEntry("  Hide",4);
      }
    }
    if(nzvents>0){
      if(visVents==1){
        glutAddMenuEntry("*Vents",14);
      }
      else{
        glutAddMenuEntry("Vents",14);
      }
      if(ReadZoneFile==1){
        if(visVentSolid==1){
          glutAddMenuEntry("   *Solid",11);
        }
        else{
          glutAddMenuEntry("   Solid",11);
        }
        if(visVentLines==1){
          glutAddMenuEntry("   *Lines",12);
        }
        else{
          glutAddMenuEntry("   Lines",12);
        }
        if(visVentSolid==0&&visVentLines==0){
          glutAddMenuEntry("   *Hide",13);
        }
        else{
          glutAddMenuEntry("   Hide",13);
        }
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

      for(iii=0;iii<nstreak_value;iii++){
        if(iii==streak_index){
          sprintf(streaklabel,"*%f",streak_rvalue[iii]);
        }
        else{
          sprintf(streaklabel,"%f",streak_rvalue[iii]);
        }
        trimzeros(streaklabel);
        strcat(streaklabel," s");
        glutAddMenuEntry(streaklabel,iii);
      }
    }
    glutAddMenuEntry("-",-1);
    if(showstreakhead==1){
      glutAddMenuEntry("*Particle Head",-3);
    }
    else{
      glutAddMenuEntry("Particle Head",-3);
    }
    glutAddMenuEntry("Hide",-2);

// allocate memory for particle property sub-menus

    if(npart5prop*npartclassinfo>0){
      NewMemory((void **)&particlepropshowsubmenu,npart5prop*npartclassinfo*sizeof(int));
    }

      ntypes=0;
      for(i=0;i<npart5prop;i++){
        part5prop *propi;

        propi = part5propinfo + i;
        if(propi->display==0)continue;
        for(j=0;j<npartclassinfo;j++){
          part5class *partclassj;

          if(propi->class_present[j]==0)continue;
          partclassj = partclassinfo + j;
          if(partclassj->kind==HUMANS)continue;
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
              strcpy(menulabel,"using:");
            }
            else{
              strcpy(menulabel,"using:");
            }
            glutAddMenuEntry(menulabel,-10-5*j);
            if(partclassj->vis_type==PART_POINTS){
              glutAddMenuEntry("    *points",-10-5*j-PART_POINTS);
            }
            else{
              glutAddMenuEntry("    points",-10-5*j-PART_POINTS);
            }
            if(partclassj->col_diameter>=0||partclassj->device_name!=NULL){
              if(partclassj->vis_type==PART_SPHERES){
                glutAddMenuEntry("    *spheres",-10-5*j-PART_SPHERES);
              }
              else{
                glutAddMenuEntry("    spheres",-10-5*j-PART_SPHERES);
              }
            }
            if(partclassj->col_length>=0||partclassj->device_name!=NULL||
              (partclassj->col_u_vel>=0&&partclassj->col_v_vel>=0&&partclassj->col_w_vel>=0)){
              if(partclassj->vis_type==PART_LINES){
                glutAddMenuEntry("    *lines",-10-5*j-PART_LINES);
              }
              else{
                glutAddMenuEntry("    lines",-10-5*j-PART_LINES);
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
      glutAddMenuEntry("Color with:",-1);
      for(i=0;i<npart5prop;i++){
        part5prop *propi;

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
    
      if(part5show==0)glutAddMenuEntry("  *Hide",-4);
      if(part5show==1)glutAddMenuEntry("  Hide",-4);
      glutAddMenuEntry("-",-1);

      glutAddMenuEntry("Draw",-1);
      ntypes=0;
      for(i=0;i<npart5prop;i++){
        part5prop *propi;

        propi = part5propinfo + i;
        if(propi->display==0)continue;
        for(j=0;j<npartclassinfo;j++){
          part5class *partclassj;

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
        glutAddMenuEntry("  Show All",-2);
        glutAddMenuEntry("  Hide All",-3);
      }
      glutAddMenuEntry("-",-1);
      if(streak5show==1){
        glutAddSubMenu("*Streaks",particlestreakshowmenu);
      }
      else{
        glutAddSubMenu("Streaks",particlestreakshowmenu);
      }
      glutAddMenuEntry("-",-1);
      if(show_tracers_always==0)glutAddMenuEntry("Show tracers always",-6);
      if(show_tracers_always==1)glutAddMenuEntry("*Show tracers always",-6);
    }

    CREATEMENU(humanpropshowmenu,ParticlePropShowMenu);
    if(npart5prop>=0){
      glutAddMenuEntry("Color:",-1);
      for(i=0;i<npart5prop;i++){
        part5prop *propi;

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
    
      if(part5show==0)glutAddMenuEntry("  *Hide",-5);
      if(part5show==1)glutAddMenuEntry("  Hide",-5);
      glutAddMenuEntry("-",-1);
      glutAddMenuEntry("Type:",-1);
      ntypes=0;
      for(i=0;i<npart5prop;i++){
        part5prop *propi;

        propi = part5propinfo + i;
        if(propi->display==0)continue;
        for(j=0;j<npartclassinfo;j++){
          part5class *partclassj;

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
          glutAddMenuEntry(menulabel,-10-5*j);
        }
        break;
      }
      if(streak5show==1){
        glutAddSubMenu("  *Streaks",particlestreakshowmenu);
      }
      else{
        glutAddSubMenu("  Streaks",particlestreakshowmenu);
      }
      if(ntypes>1){
        glutAddMenuEntry("  Show All",-2);
        glutAddMenuEntry("  Hide All",-3);
      }
    }
  }

/* --------------------------------particle show menu -------------------------- */

  if(npart_files>0&&nevac!=npart_files){
    CREATEMENU(particleshowmenu,ParticleShowMenu);
    for(ii=0;ii<npart_files;ii++){
      particle *parti;

      i = partorderindex[ii];
      parti = partinfo + i;
      if(parti->loaded==0)continue;
      if(parti->evac==1)continue;
      STRCPY(menulabel,"");
      if(parti->display==1)STRCAT(menulabel,"*");
      STRCAT(menulabel,parti->menulabel);
      glutAddMenuEntry(menulabel,-1-i);
    }
    glutAddMenuEntry("-",999);
    if(plotstate==DYNAMIC_PLOTS&&visSmokePart!=0){
      if(visSmokePart==2)glutAddMenuEntry("*Particles",1);
      if(visSmokePart==1)glutAddMenuEntry("#Particles",1);
    }
    else{
      glutAddMenuEntry("Particles",1);
    }
    if(staticframe0==1){
      if(visStaticSmoke==1){
        glutAddMenuEntry("*Particles (static)",5);
      }
      else{
        glutAddMenuEntry("Particles (static)",5);
      }
    }
    if(havesprinkpart==1){
      if(plotstate==DYNAMIC_PLOTS&&visSprinkPart==1){
        glutAddMenuEntry("*Droplets",2);
      }
      else{
        glutAddMenuEntry("Droplets",2);
      }
    }
    showall=0;
    if(plotstate==DYNAMIC_PLOTS){
      if(visSprinkPart==1&&visSmokePart!=0)showall=1;
      if(staticframe0==1&&visStaticSmoke==0)showall=0;
    }
    glutAddMenuEntry("-",999);
    if(showall==1){
      glutAddMenuEntry("*Show All",3);
    }
    else{
      glutAddMenuEntry("Show All",3);
    }
    if(plotstate==DYNAMIC_PLOTS){
      hideall=1;
      if(visSmokePart!=0)hideall=0;
      if(havesprinkpart==1&&visSprinkPart==1)hideall=0;
      if(staticframe0==1&&visStaticSmoke==1)hideall=0;
      if(hideall==1){
        glutAddMenuEntry("*Hide All",4);
      }
      else{
        glutAddMenuEntry("Hide All",4);
      }
    }
  }


/* --------------------------------Evac show menu -------------------------- */

  if(nevac>0){
    CREATEMENU(evacshowmenu,EvacShowMenu);
    for(ii=0;ii<npart_files;ii++){
      particle *parti;

      i = partorderindex[ii];
      parti = partinfo + i;
      if(parti->loaded==0)continue;
      if(parti->evac==0)continue;
      STRCPY(menulabel,"");
      if(parti->display==1)STRCAT(menulabel,"*");
      STRCAT(menulabel,parti->menulabel);
      glutAddMenuEntry(menulabel,-1-i);
    }
    glutAddMenuEntry("-",999);
    glutAddMenuEntry("Show All",3);
    if(plotstate==DYNAMIC_PLOTS){
      glutAddMenuEntry("Hide All",4);
    }
  }

/* --------------------------------lighting menu -------------------------- */

  if(showlightmenu==1){
    CREATEMENU(lightingmenu,LightingMenu);
    glutAddMenuEntry("OpenGL Lighting",-1);
    if(visLIGHT0==1)glutAddMenuEntry("  *Right",1);
    if(visLIGHT0==0)glutAddMenuEntry("  Right",1);
    if(visLIGHT1==1)glutAddMenuEntry("  *Left",2);
    if(visLIGHT1==0)glutAddMenuEntry("  Left",2);
    glutAddMenuEntry("  Flip",3);
    glutAddMenuEntry("  All",4);
    glutAddMenuEntry("  None",5);
#ifdef pp_LIGHT
    if(show_smokelighting==1)glutAddMenuEntry("*Smoke Lighting",6);
    if(show_smokelighting==0)glutAddMenuEntry("Smoke Lighting",6);
#endif
  }

  /* --------------------------------smoke3d showmenu -------------------------- */
  if(nsmoke3d_files>0&&Read3DSmoke3DFile==1){
    {
      if(nsmoke3dloaded>0){
        CREATEMENU(smoke3dshowmenu,Smoke3DShowMenu);
        for(i=0;i<nsmoke3d_files;i++){
          smoke3d *smoke3di;

          smoke3di = smoke3dinfo + i;
          if(smoke3di->loaded==0)continue;
          strcpy(menulabel,"");
          if(smoke3di->display==1)strcpy(menulabel,"*");
          strcat(menulabel,smoke3di->menulabel);
          glutAddMenuEntry(menulabel,i);
        }
        if(nsmoke3d_files>1){
          glutAddMenuEntry("-",-999);
          glutAddMenuEntry("Show All",SHOW_ALL);
          glutAddMenuEntry("Hide All",HIDE_ALL);
        }
      }
    }
  }

/* --------------------------------iso level menu -------------------------- */

  if(loaded_isomesh!=NULL&&niso_files>0&&ReadIsoFile==1){
    CREATEMENU(isolevelmenu,IsoShowMenu);
    if(loaded_isomesh->nisolevels>0&&loaded_isomesh->showlevels!=NULL){
      showflag=1;
      hideflag=1;
      for(i=0;i<loaded_isomesh->nisolevels;i++){
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
      if(showflag==1)glutAddMenuEntry("*Show All Levels",99);
      if(showflag==0)glutAddMenuEntry("Show All Levels",99);
      if(hideflag==1)glutAddMenuEntry("*Hide All Levels",98);
      if(hideflag==0)glutAddMenuEntry("Hide All Levels",98);
      glutAddMenuEntry("---",93);
      if(transparent_state==ALL_SOLID)glutAddMenuEntry("*All Levels SOLID",94);
      if(transparent_state!=ALL_SOLID)glutAddMenuEntry("All Levels SOLID",94);
      if(transparent_state==ALL_TRANSPARENT)glutAddMenuEntry("*All Levels Transparent",95);
      if(transparent_state!=ALL_TRANSPARENT)glutAddMenuEntry("All Levels Transparent",95);
      if(transparent_state==MIN_SOLID)glutAddMenuEntry("*Minimum Level Solid",96);
      if(transparent_state!=MIN_SOLID)glutAddMenuEntry("Minimum Level Solid",96);
      if(transparent_state==MAX_SOLID)glutAddMenuEntry("*Maximum Level Solid",97);
      if(transparent_state!=MAX_SOLID)glutAddMenuEntry("Maximum Level Solid",97);
    }
    else{
      glutAddMenuEntry("Show",99);
      if(visAIso==0)glutAddMenuEntry("*Hide",98);
      if(visAIso!=0)glutAddMenuEntry("Hide",98);
    }

/* --------------------------------iso show menu -------------------------- */

    if(niso_files>0&&ReadIsoFile==1){
      mesh *hmesh;
      iso *iso2;

      CREATEMENU(isoshowmenu,IsoShowMenu);
      iso2=NULL;
      for(ii=0;ii<niso_files;ii++){
        iso *isoi;

        i = isoorderindex[ii];
        isoi = isoinfo + i;
        if(isoi->loaded==0)continue;
        if(iso2==NULL&&isoi->type==iisotype){
          iso2=isoi;
        }
        if(plotstate==DYNAMIC_PLOTS&&isoi->display==1&&isoi->type==iisotype){
          iso2=isoi;
          STRCPY(menulabel,check);
          STRCAT(menulabel,isoi->menulabel);  
        }
        else{
          STRCPY(menulabel,isoi->menulabel);
        }
        glutAddMenuEntry(menulabel,1000+i);
      }
      if(iso2!=NULL){
        glutAddMenuEntry("-",999);
        STRCPY(menulabel,"Show All ");
        STRCAT(menulabel,iso2->surface_label.longlabel);
        STRCAT(menulabel," isosurfaces");
        glutAddMenuEntry(menulabel,10001);
        STRCPY(menulabel,"Hide All ");
        STRCAT(menulabel,iso2->surface_label.longlabel);
        STRCAT(menulabel," isosurfaces");
        glutAddMenuEntry(menulabel,10002);
      }
      glutAddMenuEntry("-",999);
      if(visAIso==1)glutAddMenuEntry("*Solid",1);
      if(visAIso!=1)glutAddMenuEntry("Solid",1);
      if(visAIso==2)glutAddMenuEntry("*Outline",2);
      if(visAIso!=2)glutAddMenuEntry("Outline",2);
      if(visAIso!=3)glutAddMenuEntry("Points",3);
      if(visAIso==3)glutAddMenuEntry("*Points",3);
      hmesh=meshinfo+highlight_mesh;
      if(hmesh->isofilenum!=-1){
        STRCPY(levellabel,isoinfo[hmesh->isofilenum].surface_label.shortlabel);
        STRCAT(levellabel," Levels");
        glutAddSubMenu(levellabel,isolevelmenu);
      }
      if(niso_compressed==0){
        if(isonormtype==1)glutAddMenuEntry("*Smooth",4);
        if(isonormtype==0)glutAddMenuEntry("Smooth",4);
      }
      if(showisonormals==1)glutAddMenuEntry("*Show normals",5);
      if(showisonormals==0)glutAddMenuEntry("Show normals",5);
    }
  }

/* -------------------------------- colorbarmenu -------------------------- */

  CREATEMENU(colorbarmenu,ColorBarMenu);
  {
    colorbardata *cbi;
    char ccolorbarmenu[256];

    glutAddMenuEntry("Colorbar:",-999);
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
  glutAddMenuEntry("-",-999);
  glutAddMenuEntry("Modify colorbar:",-999);
  if(show_extremedata==1){
    glutAddMenuEntry("  *Highlight extreme data",-7);
  }
  else{
    glutAddMenuEntry("  Highlight extreme data",-7);
  }
  if(colorband==1){
    glutAddMenuEntry("  *Narrow colorband",-6);
  }
  else{
    glutAddMenuEntry("  *Wide colorband",-6);
  }
  if(colorbarflip==1){
    glutAddMenuEntry("  *Flip colorbar",-2);
  }
  else{
    glutAddMenuEntry("  Flip colorbar",-2);
  }
  glutAddMenuEntry("  Cycle colorbar",-3);
  glutAddMenuEntry("  Reset",-4);
  glutAddMenuEntry("-",-999);
  glutAddMenuEntry("Shades:",-999);
  if(background_flip==1){
    glutAddMenuEntry("  *Flip background",-11);
  }
  else{
    glutAddMenuEntry("  Flip background",-11);
  }
  if(setbw==0){
    glutAddMenuEntry("  *Color/BW",-12);
  }
  else{
    glutAddMenuEntry("  Color/*BW",-12);
  }
  if(transparentflag==1){
    glutAddMenuEntry("  *Transparent (data)",-13);
  }
  else{
    glutAddMenuEntry("  Transparent (data)",-13);
  }

/* --------------------------------showVslice menu -------------------------- */
  if(nvsliceloaded==0){
    vd_shown=NULL;
  }
  if(nvslice>0&&nvsliceloaded>0){
    CREATEMENU(showvslicemenu,ShowVSliceMenu);
    nvsliceloaded0=0;
    for(i=0;i<nvslice;i++){
      vslice *vd;

      vd = vsliceinfo + i;
      if(vd->loaded==0)continue;
      nvsliceloaded0++;
      STRCPY(menulabel,"");
      if(plotstate==DYNAMIC_PLOTS&&sliceinfo[vd->ival].type==islicetype&&vd->display==1){
        vd_shown=vd;
        STRCPY(menulabel,"*");
      }
      STRCAT(menulabel,sliceinfo[vd->ival].menulabel2);
      glutAddMenuEntry(menulabel,i);
    }
    if(show_slice_in_obst==1)glutAddMenuEntry("*Show vector slice in blockage",-11);
    if(show_slice_in_obst==0)glutAddMenuEntry("Show vector slice in blockage",-11);
    if(offset_slice==1)glutAddMenuEntry("*Offset vector slice",-12);
    if(offset_slice==0)glutAddMenuEntry("Offset vector slice",-12);
    if(vd_shown!=NULL&&nvsliceloaded0!=0){
      glutAddMenuEntry("",-10);
      STRCPY(menulabel,"Show All ");
      STRCAT(menulabel,sliceinfo[vd_shown->ival].label.longlabel);
      STRCAT(menulabel," vector slices");
      glutAddMenuEntry(menulabel,SHOW_ALL);
      STRCPY(menulabel,"Hide All ");
      STRCAT(menulabel,sliceinfo[vd_shown->ival].label.longlabel);
      STRCAT(menulabel," vector slices");
      glutAddMenuEntry(menulabel,HIDE_ALL);
    }
  }
  if(nslice_files>0&&nmultislices<nslice_files){
    CREATEMENU(showmultislicemenu,ShowMultiSliceMenu);
    for(i=0;i<nmultislices;i++){
      mslicei = multisliceinfo + i;
      if(mslicei->loaded==0)continue;
      if(plotstate==DYNAMIC_PLOTS&&mslicei->display!=0&&mslicei->type==islicetype){
        if(mslicei->display==1){
          STRCPY(menulabel,"*");
          STRCAT(menulabel,mslicei->menulabel2);
        }
        else if(mslicei->display==-1){
          STRCPY(menulabel,"#");
          STRCAT(menulabel,mslicei->menulabel2);
        }
      }
      else{
        STRCPY(menulabel,mslicei->menulabel2);
      }
      glutAddMenuEntry(menulabel,i);
    }
    if(show_slice_in_obst==1)glutAddMenuEntry("*Show multi slice in blockage",-11);
    if(show_slice_in_obst==0)glutAddMenuEntry("Show multi slice in blockage",-11);
    if(offset_slice==1)glutAddMenuEntry("*Offset slice",-12);
    if(offset_slice==0)glutAddMenuEntry("Offset slice",-12);
  }

/* --------------------------------showslice menu -------------------------- */
  if(nsliceloaded==0){
    sd_shown=NULL;
  }
  if(nslice_files>0&&nsliceloaded>0){
    CREATEMENU(showhideslicemenu,ShowHideSliceMenu);
    for(ii=0;ii<nslice_loaded;ii++){
      slice *sd;

      i = slice_loaded_list[ii];
      sd = sliceinfo + i;
      if(sd_shown==NULL&&sd->type==islicetype){
        sd_shown = sd;
      }
      if(plotstate==DYNAMIC_PLOTS&&sd->display==1&&sd->type==islicetype){
        sd_shown=sd;
        STRCPY(menulabel,check);
        STRCAT(menulabel,sd->menulabel2);  
      }
      else{
        STRCPY(menulabel,sd->menulabel2);  
      }
      glutAddMenuEntry(menulabel,i);
    }
    glutAddMenuEntry("-",-10);
    if(have_terrain_slice()==1){
      if(planar_terrain_slice==1)glutAddMenuEntry("*Planar terrain slice",-13);
      if(planar_terrain_slice==0)glutAddMenuEntry("Planar terrain slice",-13);
    }
    if(show_slice_in_obst==1)glutAddMenuEntry("*Show slice in blockage",-11);
    if(show_slice_in_obst==0)glutAddMenuEntry("Show slice in blockage",-11);
    if(offset_slice==1)glutAddMenuEntry("*Offset slice",-12);
    if(offset_slice==0)glutAddMenuEntry("Offset slice",-12);
    if(nsliceloaded>0&&sd_shown!=NULL){
      glutAddMenuEntry("-",-10);
      STRCPY(menulabel,"Show All ");
      STRCAT(menulabel,sd_shown->label.longlabel);
      STRCAT(menulabel," slices");
      glutAddMenuEntry(menulabel,SHOW_ALL);
      STRCPY(menulabel,"Hide All ");
      STRCAT(menulabel,sd_shown->label.longlabel);
      STRCAT(menulabel," slices");
      glutAddMenuEntry(menulabel,HIDE_ALL);
    }
  }

/* -------------------------------- avatartour menu -------------------------- */

  CREATEMENU(avatarevacmenu,AvatarEvacMenu);
  if(navatar_types>0){
    if(iavatar_evac==-1){
      glutAddMenuEntry("*Defined in Evac File",-1);
    }
    else{
      glutAddMenuEntry("Defined in Evac File",-1);
    }
    glutAddMenuEntry("-",-999);
    for(i=0;i<navatar_types;i++){
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

      touri = tourinfo + selectedtour_index;
      strcpy(menulabel,"For ");
      strcat(menulabel,touri->label);
      glutAddMenuEntry(menulabel,-999);
      glutAddMenuEntry("-",-999);
    }

    for(i=0;i<navatar_types;i++){
      strcpy(menulabel,"");
      if(tourlocus_type==2&&iavatar_types==i){
        strcat(menulabel,"*");
      }
      strcat(menulabel,avatar_types[i]->label);
      glutAddMenuEntry(menulabel,-23-i);
    }
  }

  /* -------------------------------- showtour menu -------------------------- */

  CREATEMENU(showtourmenu,TourMenu);

    for(i=0;i<ntours;i++){
      tourdata *touri;
      int glui_avatar_index;

      touri = tourinfo + i;
      if(touri->display==1){
        STRCPY(menulabel,"");
        if(selectedtour_index==i){
          STRCAT(menulabel,"@");
        }
        STRCAT(menulabel,check);
        STRCAT(menulabel,touri->menulabel);  
      }
      else{
        STRCPY(menulabel,touri->menulabel);
      }
      glui_avatar_index = touri->glui_avatar_index;
      if(glui_avatar_index>=0&&glui_avatar_index<navatar_types){
        sv_object *avatari;

        avatari=avatar_types[glui_avatar_index];
        strcat(menulabel,"(");
        strcat(menulabel,avatari->label);
        strcat(menulabel,")");
      }
      glutAddMenuEntry(menulabel,i);
    }
    glutAddMenuEntry("-",-999);
    glutAddMenuEntry("Show All Tours",-3);
    glutAddMenuEntry("Hide All Tours",-2);
    if(selectedtour_index>=0){
      strcpy(menulabel,"");
      if(viewtourfrompath==1)strcat(menulabel,"*");
      strcat(menulabel,"View from ");
      strcat(menulabel,tourinfo[selectedtour_index].label);
      glutAddMenuEntry(menulabel,-5);
    }

    /* --------------------------------tour menu -------------------------- */

    CREATEMENU(tourmenu,TourMenu);
      
    glutAddMenuEntry("New...",-12);
    if(ntours>0){
      if(trainer_mode==0){
        glutAddMenuEntry("Manual",-2);
      }
      else{
        glutAddMenuEntry("Crawl",-8);
        glutAddMenuEntry("Walk",-9);
        glutAddMenuEntry("Overview",-10);
      }
      if(trainer_mode==0||(trainer_mode==1&&ntours>0))glutAddMenuEntry("-",-999);
    }
  if(ntours>0){
    glutAddSubMenu("Show/Hide",showtourmenu);
  }


/* --------------------------------showhide menu -------------------------- */

  CREATEMENU(showhidemenu,ShowHideMenu);
  if(ntotal_blockages>0||isZoneFireModel==0){
    glutAddSubMenu("Geometry",geometrymenu);
  }
  glutAddSubMenu("Labels",labelmenu);
  {
    int human_present=0;
    int particle_present=0;

    if(npart5loaded>0){
      for(ii=0;ii<npart_files;ii++){
        particle *parti;

        parti = partinfo + ii;
        if(parti->loaded==0)continue;
        if(parti->evac==1)human_present=1;
        if(parti->evac==0)particle_present=1;
      }
      if(particle_present==1){
        glutAddSubMenu("Particles",particlepropshowmenu);
      }
      if(human_present==1){
        glutAddSubMenu("Humans",humanpropshowmenu);
      }
    }
  }
  if(nevac>0&&navatar_types>0){
    glutAddSubMenu("Use Avatar:",avatarevacmenu);
  }

  if(npart4loaded>0){
    if(havesprinkpart!=0||staticframe0!=0||npartloaded>1){
      glutAddSubMenu("Particles",particleshowmenu);
    }
    else{
      if(ReadPartFile==1&&showsmoke==1)glutAddMenuEntry("*Particles",1);
      if(ReadPartFile==1&&showsmoke==0)glutAddMenuEntry("Particles",1);
    }
  }
  if(partinfo!=NULL&&partinfo[0].version==1){
  }
  else{
    if(nevacloaded>1){
      glutAddSubMenu("Evacuation",evacshowmenu);
    }
    else{
      if(ReadEvacFile==1&&showevac==1)glutAddMenuEntry("*Evacuation",13);
      if(ReadEvacFile==1&&showevac==0)glutAddMenuEntry("Evacuation",13);
    }
  }
  if(ReadIsoFile==1){
    int niso_loaded=0;
    for(i=0;i<niso_files;i++){
      iso *isoi;

      isoi = isoinfo + i;
      if(isoi->loaded==1)niso_loaded++;
    }

    if(niso_loaded>1){
     glutAddSubMenu("Animated Surfaces",isoshowmenu);
    }
    else{
     glutAddSubMenu("Animated Surface",isoshowmenu);
    }
  }
  if(Read3DSmoke3DFile==1&&nsmoke3dloaded>0){
      glutAddSubMenu("3D Smoke",smoke3dshowmenu);
  }

  nvslice0=0, nvslice1=0, nvslice2=0;
  nvsliceloaded0=0, nvsliceloaded1=0;
  nvsliceloaded2=0;
  for(i=0;i<nvslice;i++){
    vslice *vd;

    vd = vsliceinfo+i;
    switch (vd->vec_type){
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
    glutAddSubMenu("Animated Vector Slices",showvslicemenu);
  }

  if(nsliceloaded>0){
    glutAddSubMenu("Animated Slices",showhideslicemenu);
    if(nmultislices<nslice_files){
      glutAddSubMenu("Animated Multi-Slices",showmultislicemenu);
    }
  }

  if(ReadPlot3dFile==1){
    glutAddSubMenu("Plot3d",plot3dshowmenu);
  }
  if(npatchloaded>0){
    glutAddSubMenu("Boundaries",showpatchmenu);
  }
  if(nzone>0&&(ReadZoneFile==1||nzvents>0))glutAddSubMenu("Zone",zoneshowmenu);
  if(ReadTargFile==1){
    if(showtarget==1)glutAddMenuEntry("*Targets",2);
    if(showtarget==0)glutAddMenuEntry("Targets",2);
  }
  if(nobject_defs>0){
    int num_activedevices=0;

    if(multiprop==1){
     // glutAddSubMenu("Properties",propmenu);
    }

    for(i=0;i<nobject_defs;i++){
      sv_object *obj_typei;

      obj_typei = object_defs[i];
      if(obj_typei->used==1)num_activedevices++;
    }

    if(num_activedevices>0){
      /*
      if(isZoneFireModel==0||(isZoneFireModel==1&&num_activedevices>1)){
        glutAddSubMenu("Objects",showobjectsmenu);
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
          glutAddMenuEntry(obj_menu,i);
        }
      }
    }
  }
  if(ntc_total>0){
    if(isZoneFireModel==1){
      if(visSensor==1)glutAddMenuEntry("*Targets",9);
      if(visSensor==0)glutAddMenuEntry("Targets",9);
      if(hasSensorNorm==1){
        if(visSensorNorm==1)glutAddMenuEntry("*Target Orientation",14);
        if(visSensorNorm==0)glutAddMenuEntry("Target Orientation",14);
      }
    }
    else{
      if(visSensor==1)glutAddMenuEntry("*Thermocouples",9);
      if(visSensor==0)glutAddMenuEntry("Thermocouples",9);
      if(hasSensorNorm==1){
        if(visSensorNorm==1)glutAddMenuEntry("*Thermocouple Norms",14);
        if(visSensorNorm==0)glutAddMenuEntry("Thermocouple Norms",14);
      }
    }
  }
  if(titlesafe_offset==0)glutAddMenuEntry("Offset Window",12);
  if(titlesafe_offset!=0)glutAddMenuEntry("*Offset Window",12);
  if(ntextures_loaded_used>0){
    glutAddSubMenu("Textures",textureshowmenu);
  }

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
  if(frameratevalue==2001)glutAddMenuEntry("*Real time",2001);
  if(frameratevalue!=2001)glutAddMenuEntry("Real time",2001);
  if(frameratevalue==2002)glutAddMenuEntry("*2 x Real time",2002);
  if(frameratevalue!=2002)glutAddMenuEntry("2 x Real time",2002);
  if(frameratevalue==2004)glutAddMenuEntry("*4 x Real time",2004);
  if(frameratevalue!=2004)glutAddMenuEntry("4 x Real time",2004);
  if(frameratevalue!=1000)glutAddMenuEntry("Unlimited",1000);
  if(frameratevalue==1000)glutAddMenuEntry("*Unlimited",1000);
  if(frameratevalue<0){glutAddMenuEntry("*Step",-1);}
   else{glutAddMenuEntry("Step",-1);}

/* --------------------------------render menu -------------------------- */
   {
     char renderwindow[1024];
     char renderwindow2[1024];
     char renderwindow3[1024];
     char renderwindow4[1024];


     sprintf(renderwindow,"%i%s%i(current)",screenWidth,"x",screenHeight);
     strcpy(renderwindow2,"*");
     strcat(renderwindow2,renderwindow);
     sprintf(renderwindow3,"%i%s%i(2*current)",2*screenWidth,"x",2*screenHeight);
     strcpy(renderwindow4,"*");
     strcat(renderwindow4,renderwindow3);


  CREATEMENU(rendermenu,RenderMenu);
  glutAddMenuEntry("SIZE",10000);
  if(renderW==320){
    glutAddMenuEntry("*320x240",Render320);
    glutAddMenuEntry("640x480",Render640);
    glutAddMenuEntry(renderwindow,RenderWindow);
    glutAddMenuEntry(renderwindow3,Render2Window);
  }
  else if(renderW==640){
    glutAddMenuEntry("320x240",Render320);
    glutAddMenuEntry("*640x480",Render640);
    glutAddMenuEntry(renderwindow,RenderWindow);
    glutAddMenuEntry(renderwindow3,Render2Window);
  }
  else if(renderW==2*screenWidth){
    glutAddMenuEntry("320x240",Render320);
    glutAddMenuEntry("640x480",Render640);
    glutAddMenuEntry(renderwindow,RenderWindow);
    glutAddMenuEntry(renderwindow4,Render2Window);
  }
  else{
    glutAddMenuEntry("320x240",Render320);
    glutAddMenuEntry("640x480",Render640);
    glutAddMenuEntry(renderwindow2,RenderWindow);
    glutAddMenuEntry(renderwindow3,Render2Window);
  }
  glutAddMenuEntry("-",10000);
  glutAddMenuEntry("TYPE",10000);
  if(renderfiletype==0){
    glutAddMenuEntry("*PNG",RenderPNG);
#ifdef pp_JPEG
    glutAddMenuEntry("JPEG",RenderJPEG);
#endif
#ifdef pp_GDGIF
    glutAddMenuEntry("GIF",RenderGIF);
#endif
  }
  if(renderfiletype==1){
    glutAddMenuEntry("PNG",RenderPNG);
#ifdef pp_JPEG
    glutAddMenuEntry("*JPEG",RenderJPEG);
#endif
#ifdef pp_GDGIF
    glutAddMenuEntry("GIF",RenderGIF);
#endif
  }
  if(renderfiletype==2){
    glutAddMenuEntry("PNG",RenderPNG);
#ifdef pp_JPEG
    glutAddMenuEntry("JPEG",RenderJPEG);
#endif
#ifdef pp_GDGIF
    glutAddMenuEntry("*GIF",RenderGIF);
#endif
  }
  glutAddMenuEntry("-",10000);
  glutAddMenuEntry("NUMBER",10000);
  glutAddMenuEntry("One Frame",RenderOnce);
  update_glui_render();
  if(RenderTime==1||touring==1){
   // if(render_double_state==0){
    glutAddMenuEntry("All Frames",1);
    glutAddMenuEntry("Every 2nd Frame",2);
    glutAddMenuEntry("Every 3rd Frame",3);
    glutAddMenuEntry("Every 4th Frame",4);
    glutAddMenuEntry("Every 5th Frame",5);
    glutAddMenuEntry("Every 10th Frame",10);
    glutAddMenuEntry("Every 20th Frame",20);
    glutAddMenuEntry("Cancel",RenderCancel);
   // }
  }
   }

   /* --------------------------------viewpoint menu -------------------------- */


  CREATEMENU(dialogmenu,DialogMenu);
  if(nterraininfo>0){
    if(showwui==1)glutAddMenuEntry("*WUI Display... ALT+w",26);
    if(showwui==0)glutAddMenuEntry("WUI Display... ALT+w...",26);
  }
  if(showclip==1)glutAddMenuEntry("*Clip Geometry...  ALT+c",18);
  if(showclip==0)glutAddMenuEntry("Clip Geometry...  ALT+c",18);
#ifdef pp_COMPRESS
  if(smokezippath!=NULL&&(npatch_files>0||nsmoke3d_files>0||nslice_files>0)){
    if(showbounds==1)glutAddMenuEntry("*Compression/Smokezip...  ALT+z",24);
    if(showbounds==0)glutAddMenuEntry("Compression/Smokezip...  ALT+z",24);
  }
#endif
  if(showcolorbar==1)glutAddMenuEntry("*Customize Colorbar...  ALT+c",23);
  if(showcolorbar==0)glutAddMenuEntry("Customize Colorbar...  ALT+c",23);
  if(showlabels==1)glutAddMenuEntry("*Display...  ALT+d",22);
  if(showlabels==0)glutAddMenuEntry("Display...  ALT+d",22);
  if(isZoneFireModel==0){
    if(showedit==1)glutAddMenuEntry("*Examine Blockages...  ALT+e",16);
    if(showedit==0)glutAddMenuEntry("Examine Blockages...  ALT+e",16);
  }
  if(showbounds==1)glutAddMenuEntry("*File/Bound/Script Settings...  ALT+f",14);
  if(showbounds==0)glutAddMenuEntry("File/Bound/Script Settings...  ALT+f",14);
  if(showmotion==1)glutAddMenuEntry("*Motion/View/Render...  ALT+m",15);
  if(showmotion==0)glutAddMenuEntry("Motion/View/Render...  ALT+m",15);
#ifdef pp_SHOOTER
  if(showshooterDLG==1)glutAddMenuEntry("*Particle Tracking...",27);
  if(showshooterDLG==0)glutAddMenuEntry("Particle Tracking...",27);
#endif
  if(nsmoke3d_files>0){
    if(showbounds==1)glutAddMenuEntry("*Smoke3D Parameters...  ALT+s",20);
    if(showbounds==0)glutAddMenuEntry("Smoke3D Parameters...  ALT+s",20);
  }
  if(showgluistereo==1)glutAddMenuEntry("*Stereo Parameters...",19);
  if(showgluistereo==0)glutAddMenuEntry("Stereo Parameters...",19);
  if(showgluitour==1)glutAddMenuEntry("*Tours...  ALT+t",21);
  if(showgluitour==0)glutAddMenuEntry("Tours...  ALT+t",21);
  if(trainer_active==1){
    if(showtrainer==1)glutAddMenuEntry("*Trainer...",25);
    if(showtrainer==0)glutAddMenuEntry("Trainer...",25);
  }
  glutAddMenuEntry("-",-1);
  glutAddMenuEntry("Close All Dialogs  ALT+x",-2);

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

  /* -------------------------------- font menu -------------------------- */

  if(showfontmenu==1){
    CREATEMENU(fontmenu,FontMenu);
    switch (fontindex){
    case SMALL_FONT:
      glutAddMenuEntry("*Normal",0);
      glutAddMenuEntry("Large",1);
      break;
    case LARGE_FONT:
      glutAddMenuEntry("Normal",0);
      glutAddMenuEntry("*Large",1);
      break;
    case LARGE_FONT_SAFE:
      glutAddMenuEntry("Normal",0);
      glutAddMenuEntry("Large",1);
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
  }
  if(nunitclasses>0){
    for(i=0;i<nunitclasses;i++){
      f_units *uci;

      uci = unitclasses + i;
      CheckMemory;

      CREATEMENU(uci->submenuid,UnitsMenu);

      for(j=0;j<uci->nunits;j++){
        if(uci->active==j){
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
      if(uci->visible==0)continue;
      glutAddSubMenu(uci->unitclass,uci->submenuid);
    }
    if(vishmsTimelabel==0)glutAddMenuEntry("time (h:m:s)",-2);
    if(vishmsTimelabel==1)glutAddMenuEntry("*time (h:m:s)",-2);
    glutAddMenuEntry("Reset",-1);
  }

/* --------------------------------option menu -------------------------- */

  CREATEMENU(optionmenu,OptionMenu);
  glutAddSubMenu("Colorbars/Shades",colorbarmenu);
  if(nunitclasses>0)glutAddSubMenu("Units",unitsmenu);
#ifdef pp_SHOWLIGHT
  if(showlightmenu==1)glutAddSubMenu("Lighting",lightingmenu);
#endif
  glutAddSubMenu("Rotation",rotatetypemenu);
  glutAddSubMenu("Max Frame Rate",frameratemenu);
#ifdef pp_RENDER
  glutAddSubMenu("Render",rendermenu);
#endif
  if(showfontmenu==1)glutAddSubMenu("Font Size",fontmenu);
//  glutAddSubMenu("Aperture",aperturemenu);
#ifdef pp_ALPHA
  glutAddMenuEntry("Benchmark",1);
#endif
  if(trainer_active==1)glutAddMenuEntry("Trainer Menu",2);
  
  /* --------------------------------reset menu -------------------------- */

  CREATEMENU(resetmenu,ResetMenu);
  {
    char line[256];
    camera *ca;

    if(trainer_mode==1){
      if(visBlocks==visBLOCKOutline){
        glutAddMenuEntry("*Outline",MENU_OUTLINEVIEW);
      }
      else{
        glutAddMenuEntry("Outline",MENU_OUTLINEVIEW);
      }
      glutAddMenuEntry("-",MENU_DUMMY);
    }
    if(trainer_mode==0){
      glutAddMenuEntry("Save",MENU_SAVEVIEW);
      glutAddMenuEntry("Set as Startup",MENU_STARTUPVIEW);
      glutAddSubMenu("Zoom",zoommenu);
      if(projection_type==1)glutAddMenuEntry("Switch to perspective view       ALT+v",MENU_SIZEPRESERVING);
      if(projection_type==0)glutAddMenuEntry("Switch to size preserving view   ALT+v",MENU_SIZEPRESERVING);
      glutAddMenuEntry("-",MENU_DUMMY);
    }
    for(ca=camera_list_first.next;ca->next!=NULL;ca=ca->next){
      if(trainer_mode==1&&strcmp(ca->name,"internal")==0)continue;
      strcpy(line,"");
      if(ca->view_id==selected_view){
        strcat(line,"*");
      }
      if(trainer_mode==1&&strcmp(ca->name,"external")==0){
        strcat(line,"Outside");
      }
      else{
        strcat(line,ca->name);
      }
      glutAddMenuEntry(line,ca->view_id);
    }
  }
  if(trainer_mode==0&&showtime==1){
    glutAddMenuEntry("-",MENU_DUMMY);
    glutAddMenuEntry("Time",MENU_TIMEVIEW);
  }

/* -------------------------------- disclaimer menu -------------------------- */

  CREATEMENU(disclaimermenu,DisclaimerMenu);
  glutAddMenuEntry(TITLERELEASE,1);
  {
#ifdef pp_GPU
    char version_label[256];
#endif

#ifdef BIT64
    sprintf(menulabel,"  Smokeview (64 bit) revision:%i",revision_smv);
#else
    sprintf(menulabel,"  Smokeview (32 bit) revision:%i",revision_smv);
#endif
    glutAddMenuEntry(menulabel,1);
    if(revision_fds>0){
      sprintf(menulabel,"  FDS revision:%i",revision_fds);
      glutAddMenuEntry(menulabel,1);
    }
#ifdef pp_GPU
    strcpy(version_label,"  OpenGL Version: "); 
    strcat(version_label,(char *)glGetString(GL_VERSION));
    glutAddMenuEntry(version_label,1);
    if(gpuactive==1){
      if(usegpu==1){
        strcpy(menulabel,"  GPU activated. (Press G to deactivate)");
      }
      else{
        strcpy(menulabel,"  GPU available but not in use. (Press G to activate)");
      }
    }
    else{
      strcpy(menulabel,"  GPU not available");
    }
    glutAddMenuEntry(menulabel,1);
#endif
#ifdef pp_CULL
    if(cullactive==1&&gpuactive==1){
      if(cullsmoke==1&&usegpu==1){
        strcpy(menulabel,"  Smoke culling activated. (Press C to deactivate)");
      }
      else{
        strcpy(menulabel,"  Smoke culling available but not in use. ( To activate: ");
        if(usegpu==0 && cullsmoke==1)strcat(menulabel," Press G.)");
        if(usegpu==1 && cullsmoke==0)strcat(menulabel," Press C.)");
        if(usegpu==0 && cullsmoke==0)strcat(menulabel," Press G then C.)");
      }
    }
    else{
      strcpy(menulabel,"  Smoke culling not available");
    }
    glutAddMenuEntry(menulabel,1);
#endif
  }
  glutAddMenuEntry("-",1);
  glutAddMenuEntry("The US Department of Commerce makes no warranty, expressed or",1);
  glutAddMenuEntry("implied, to users of Smokeview, and accepts no responsibility",1);
  glutAddMenuEntry("for its use. Users of Smokeview assume sole responsibility under",1);
  glutAddMenuEntry("Federal law for determining the appropriateness of its use in any",1);
  glutAddMenuEntry("particular application; for any conclusions drawn from the results",1); 
  glutAddMenuEntry("of its use; and for any actions taken or not taken as a result of",1); 
  glutAddMenuEntry("analysis performed using this tools.",1);
  glutAddMenuEntry("",1);
  glutAddMenuEntry("Smokeview and the companion program FDS is intended for use only",1);
  glutAddMenuEntry("by those competent in the fields of fluid dynamics, thermodynamics,",1);
  glutAddMenuEntry("combustion, and heat transfer, and is intended only to supplement",1);
  glutAddMenuEntry("the informed judgment of the qualified user. These software packages",1);
  glutAddMenuEntry("may or may not have predictive capability when applied to a specific",1);
  glutAddMenuEntry("set of factual circumstances.  Lack of accurate predictions could lead",1);
  glutAddMenuEntry("to erroneous conclusions with regard to fire safety.  All results",1);
  glutAddMenuEntry("should be evaluated by an informed user.",1);

/* --------------------------------help menu -------------------------- */

  {
    int displayblank=0;
  CREATEMENU(helpmenu,HelpMenu);
#ifdef WIN32
  glutAddMenuEntry("Web Links",1);
  glutAddMenuEntry("  Documentation",-3);
  glutAddMenuEntry("  Issue Tracker",-1);
  glutAddMenuEntry("  Updates",-2);
  glutAddMenuEntry("  FDS/Smokeview Website",-4);
  glutAddMenuEntry("",1);
#endif
#ifdef pp_OSX
  glutAddMenuEntry("Web Links",1);
  glutAddMenuEntry("  Documentation",-3);
  glutAddMenuEntry("  Issue Tracker",-1);
  glutAddMenuEntry("  Updates",-2);
  glutAddMenuEntry("  FDS/Smokeview Website",-4);
  glutAddMenuEntry("",1);
#endif
#ifndef WIN32
#ifndef PP_OSX
  glutAddMenuEntry("Web Resources",1);
  glutAddMenuEntry("*** Download documentation at  http://fire.nist.gov/fds/documentation.html",1);
  glutAddMenuEntry("*** Report a problem at http://code.google.com/p/fds-smv/issues/",1);
  glutAddMenuEntry("*** Check for updates at http://code.google.com/p/fds-smv/downloads/",1);
  glutAddMenuEntry("*** FDS/Smokeview website: http://fire.nist.gov/fds",1);
  glutAddMenuEntry("-",1);
#endif
#endif
  if(plotstate==DYNAMIC_PLOTS){
    glutAddMenuEntry("Animation keyboard commands",1);
    glutAddMenuEntry("  t: set/unset single time step mode",6);
    glutAddMenuEntry("  o: reset animation to the initial time",6);
    glutAddMenuEntry("  T: toggle texture method for drawing slice and boundary colors",6);
    glutAddMenuEntry("  u: reload files",6);
    glutAddMenuEntry("1-9: number of frames to skip",6);
    displayblank=1;
  }
  if(plotstate==STATIC_PLOTS){
    if(displayblank==1)glutAddMenuEntry("",1);
    glutAddMenuEntry("Plot3D keyboard commands",1);
    glutAddMenuEntry("  x,y,z: toggle contour plot visibility along x, y and z axis",3);
    glutAddMenuEntry("  p: increment plot3d variable",2);
    glutAddMenuEntry("  P: toggle cursor key mappings",2);
    glutAddMenuEntry("  v: toggle flow vector visiblity",3);
    glutAddMenuEntry("  a: change flow vector lengths",3);
    glutAddMenuEntry("  s: change interval between adjacent vectors",3);
    glutAddMenuEntry("  c: toggle between continuous and 2D stepped contours",3);
    glutAddMenuEntry("  i: toggle iso-surface visibility",2);
    displayblank=1;
  }
  if(displayblank==1)glutAddMenuEntry("",1);
  glutAddMenuEntry("Misc keyboard commands",1);
  glutAddMenuEntry("  r: render the current scene to an image file",7);
  glutAddMenuEntry("  R:   (same as r but at twice the resolution)",7);
  if(ntotal_blockages>0||isZoneFireModel==0){
    glutAddMenuEntry("  g: toggle grid visibility",2);
  }
  glutAddMenuEntry("  e: toggle between eye, world and world/level rotation motion",7);
  glutAddMenuEntry("  w: toggle clipping - use Options/Clip menu to specify clipping planes",7);
  glutAddMenuEntry("  -: decrement time step, 2D contour planes, 3D contour levels",2);
  glutAddMenuEntry("  space bar: increment time step, 2D contour planes, 3D contour levels",2);
  glutAddMenuEntry("",1);
  glutAddMenuEntry("  SHIFT+t: toggle slice file color interpolation methods",2);
  if(cellcenter_slice_active==1){
    glutAddMenuEntry("           (also, toggles cell center display on/off)",2);
  }
  glutAddMenuEntry("  ALT+v: toggle projection  method (between perspective and size preserving)",2);
  if(n_embedded_meshes>0){
    glutAddMenuEntry("  ALT+u: toggle coarse slice display in embedded mesh",2);
  }
  if(cellcenter_slice_active==1){
    glutAddMenuEntry("  ALT+y: if current slice is cell centered, toggle interpolation on/off",2);
  }
  if(caseinifilename!=NULL&&strlen(caseinifilename)>0){
    char inilabel[512];

    sprintf(inilabel,"  #: save settings to %s",caseinifilename);
    glutAddMenuEntry(inilabel,2);
  }
  else{
    glutAddMenuEntry("  #: save settings (create casename.ini file)",2);
  }
  glutAddMenuEntry("  !: snap scene's view angles",2);
  glutAddMenuEntry("  &: toggle line anti-aliasing (draw lines smoothly)",2);
  glutAddMenuEntry("",1);
  glutAddMenuEntry("Mouse Motion",1);
  glutAddMenuEntry("      : rotate around z, x axis",1);
  glutAddMenuEntry("  CTRL: translate along x, y axis",1);
  glutAddMenuEntry("   ALT: translate along z axis",1);
  if(eyeview==EYE_CENTERED){
    glutAddMenuEntry("",1);
    glutAddMenuEntry("Keyboard Motion",1);
    glutAddMenuEntry("   left/right cursor: rotate left/right",1);
    glutAddMenuEntry("      up/down cursor: move forward/backward",1);
    glutAddMenuEntry(" CTRL:up/down cursor: move forward/backward 5 times slower",1);
    glutAddMenuEntry(" SHIFT: left/right cursor: rotate 90 degrees",1);
    glutAddMenuEntry("    ALT:left/right cursor: slide left/right",1);
    glutAddMenuEntry("    ALT:   up/down cursor: slide up/down",1);
    glutAddMenuEntry("     INSERT/HOME/PageUP: tilt down/reset/tilt up",1);
  }
  }

  /* -------------------------------- target menu -------------------------- */

  if(ntarg_files>0){
    CREATEMENU(targetmenu,TargetMenu);
    for(i=0;i<ntarg_files;i++){
      if(targfilenum==i){
        STRCPY(menulabel,check);
        STRCAT(menulabel,targinfo[i].file);  
      }
      else{STRCPY(menulabel,targinfo[i].file);}
      glutAddMenuEntry(menulabel,i);
    }
    glutAddMenuEntry("Unload",-1);
    CheckMemory;
  }

  /* --------------------------------particle menu -------------------------- */

  if(npart_files>0&&nevac!=npart_files){
    CREATEMENU(unloadpartmenu,UnloadPartMenu);
    for(ii=0;ii<npart_files;ii++){
      particle *parti;

      i = partorderindex[ii];
      parti = partinfo + i;
      if(parti->loaded==1&&parti->evac==0){
        STRCPY(menulabel,parti->menulabel);  
        glutAddMenuEntry(menulabel,i);
      }
    }
    glutAddMenuEntry("Unload All",-1);

    if(nmeshes==1){
      CREATEMENU(particlemenu,ParticleMenu);
    }
    else{
      CREATEMENU(particlesubmenu,ParticleMenu);
    }
    for(ii=0;ii<npart_files;ii++){
      i = partorderindex[ii];
      if(partinfo[i].evac==1)continue;
      if(partinfo[i].loaded==1){
        STRCPY(menulabel,check);
        STRCAT(menulabel,partinfo[i].menulabel);  
      }
      else{
        STRCPY(menulabel,partinfo[i].menulabel);
      }
      glutAddMenuEntry(menulabel,i);
    }
    {
      int useitem;
      int atleastone=0;
      particle *parti, *partj;

      if(nmeshes>1){
        CREATEMENU(particlemenu,ParticleMenu);
        if(npart_files>0){
          if(partinfo->version==1){
            strcpy(menulabel,"Particles - All meshes");
            glutAddMenuEntry(menulabel,-11);
            glutAddSubMenu("Particles - Single mesh",particlesubmenu);
            glutAddMenuEntry("-",-2);
          }
          else{
            for(i=0;i<npart_files;i++){
              useitem=i;
              parti = partinfo + i;
              if(parti->evac==1)continue;
              for(j=0;j<i;j++){
                partj = partinfo + j;
                if(partj->evac==1)continue;
                if(strcmp(parti->label.longlabel,partj->label.longlabel)==0){
                  useitem=-1;
                  break;
                }
              }
              if(useitem!=-1){
                atleastone=1;
                strcpy(menulabel,parti->label.longlabel);
                strcat(menulabel," - All meshes");
                glutAddMenuEntry(menulabel,-useitem-10);
              }
            }
            glutAddSubMenu("Particles - Single mesh",particlesubmenu);
            if(atleastone==1)glutAddMenuEntry("-",-2);
          }
        }
      }
    }

    if(npartloaded<=1){
      glutAddMenuEntry("Unload",-1);
    }
     else{
       glutAddSubMenu("Unload",unloadpartmenu);
     }
  }

  if(nevac>0){
    CREATEMENU(unloadevacmenu,UnloadEvacMenu);
    for(ii=0;ii<npart_files;ii++){
      particle *parti;

      i = partorderindex[ii];
      parti = partinfo + i;
      if(parti->loaded==1&&parti->evac==0){
        STRCPY(menulabel,parti->menulabel);  
        glutAddMenuEntry(menulabel,i);
      }
    }
    glutAddMenuEntry("Unload All",-1);

    CREATEMENU(evacmenu,EvacMenu);
    {
      int nevacs=0,nevacloaded2=0;

      for(ii=0;ii<npart_files;ii++){
        particle *parti;

        parti = partinfo + ii;
        if(parti->evac==1){
          if(parti->loaded==1)nevacloaded2++;
          nevacs++;
        }
      }
      if(nevacs>1){
        strcpy(menulabel,"Humans - All meshes");
        glutAddMenuEntry(menulabel,-11);
        glutAddMenuEntry("-",-2);
      }
      for(ii=0;ii<npart_files;ii++){
        i = partorderindex[ii];
        if(partinfo[i].evac==0)continue;
        if(partinfo[i].loaded==1){
          STRCPY(menulabel,check);
          STRCAT(menulabel,partinfo[i].menulabel);  
        }
        else{
          STRCPY(menulabel,partinfo[i].menulabel);
        }
        glutAddMenuEntry(menulabel,i);
      }
      if(nevacloaded2<=1){
        glutAddMenuEntry("Unload",-1);
      }
       else{
         glutAddSubMenu("Unload",unloadevacmenu);
       }
  }
    }

/* --------------------------------unload and load vslice menus -------------------------- */

  if(nvslice>0){

    if(nmultivslices<nvslice){
      CREATEMENU(unloadmultivslicemenu,UnloadMultiVSliceMenu);
      for(i=0;i<nmultivslices;i++){
        mvslicei = multivsliceinfo + i;
        if(mvslicei->loaded!=0){
          glutAddMenuEntry(mvslicei->menulabel2,i);
        }
      }
      glutAddMenuEntry("Unload All",-1);

      nloadsubmvslicemenu=1;
      for(i=1;i<nmultivslices;i++){
        vslice *vi, *vim1;
        slice *si, *sim1;

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
      for(i=0;i<nmultivslices;i++){
        vslice *vi, *vim1;
        slice *si, *sim1;

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

        if(mvslicei->loaded==1){
          STRCPY(menulabel,"*");
          STRCAT(menulabel,mvslicei->menulabel);
          nmultisliceloaded++;
        }
        else if(mvslicei->loaded==-1){
          STRCPY(menulabel,"#");
          STRCAT(menulabel,mvslicei->menulabel);
        }
        else{
          STRCPY(menulabel,mvslicei->menulabel);
        }
        glutAddMenuEntry(menulabel,i);
      }

      nloadsubmvslicemenu=0;
      CREATEMENU(loadmultivslicemenu,LoadMultiVSliceMenu);
      for(i=0;i<nmultivslices;i++){
        vslice *vi, *vim1;
        slice *si, *sim1;

        vi = vsliceinfo + (multivsliceinfo+i)->ivslices[0];
        if(i>0)vim1 = vsliceinfo + (multivsliceinfo+i-1)->ivslices[0];
        si = sliceinfo + vi->ival;
        if(i>0)sim1 = sliceinfo + vim1->ival;
        if(i==0||(i>0&&strcmp(si->label.longlabel,sim1->label.longlabel)!=0)){
          if(si->vec_comp==0||showallslicevectors==1){
            char mlabel[1024], mlabel2[1024];

            STRCPY(mlabel,si->label.longlabel);
            if(i>0&&si->mesh_type!=sim1->mesh_type){
              sprintf(mlabel2,"*** Evac type %i meshes ***",si->mesh_type);
              glutAddMenuEntry(mlabel2,-999);
            }
            glutAddSubMenu(mlabel,loadsubmvslicemenu[nloadsubmvslicemenu]);
          }
          nloadsubmvslicemenu++;
        }
      }
      if(nmultivslices>0)glutAddMenuEntry("-",-999);
      if(showallslicevectors==0)glutAddMenuEntry("Show All Vector Slice Entries",-20);
      if(showallslicevectors==1)glutAddMenuEntry("*Show All Vector Slice Entries",-20);
      if(nmultisliceloaded>1){
        glutAddSubMenu("Unload",unloadmultivslicemenu);
      }
      else{
        glutAddMenuEntry("Unload",-1);
      }
    }

    CREATEMENU(unloadvslicemenu,UnloadVSliceMenu);
    for(ii=0;ii<nvslice;ii++){
      vslice *vd;

      i = vsliceorderindex[ii];
      vd = vsliceinfo + i;
      if(vd->loaded==0)continue;
      glutAddMenuEntry(vd->menulabel2,i);
    }
    glutAddMenuEntry("Unload All",-1);

    if(nvslice0>0){
      vslice *vd, *vdim1;
      if(nvslice>0){
        nloadsubvslicemenu=1;
        for(ii=1;ii<nvslice;ii++){
          slice *sd, *sdm1;

          i=vsliceorderindex[ii];
          vd = vsliceinfo + i;
          vdim1 = vsliceinfo + vsliceorderindex[ii-1];
          sd = sliceinfo + vd->ival;
          sdm1 = sliceinfo + vdim1->ival;
          if(strcmp(sd->label.longlabel,sdm1->label.longlabel)!=0){
            nloadsubvslicemenu++;
          }
        }
        NewMemory((void **)&loadsubvslicemenu,nloadsubvslicemenu*sizeof(int));
        for(i=0;i<nloadsubvslicemenu;i++){
          loadsubvslicemenu[i]=0;
        }
        nloadsubvslicemenu=0;
        for(ii=0;ii<nvslice;ii++){
          slice *sd, *sdm1;

          i=vsliceorderindex[ii];
          vd = vsliceinfo + i;
          if(ii>0)vdim1 = vsliceinfo + vsliceorderindex[ii-1];
          sd = sliceinfo + vd->ival;
          if(ii>0)sdm1 = sliceinfo + vdim1->ival;
          if(ii==0||strcmp(sd->label.longlabel,sdm1->label.longlabel)!=0){
            CREATEMENU(loadsubvslicemenu[nloadsubvslicemenu],LoadVSliceMenu);
            nloadsubvslicemenu++;
          }
          sd = sliceinfo + vd->ival;
          if(vd->loaded==1){
            STRCPY(menulabel,check);
            STRCAT(menulabel,sd->menulabel);
          }
          else{
            STRCPY(menulabel,sd->menulabel);
          }
          if(sd->vec_comp==0||showallslicevectors==1)glutAddMenuEntry(menulabel,i);
        }
        CREATEMENU(vslicemenu,LoadVSliceMenu);
        nloadsubvslicemenu=0;
        for(ii=0;ii<nvslice;ii++){
          slice *sd, *sdm1;

          i=vsliceorderindex[ii];
          vd = vsliceinfo + i;
          if(ii>0)vdim1 = vsliceinfo + vsliceorderindex[ii-1];
          sd = sliceinfo + vd->ival;
          if(ii>0)sdm1 = sliceinfo + vdim1->ival;
          if(ii==0||strcmp(sd->label.longlabel,sdm1->label.longlabel)!=0){
            if(sd->vec_comp==0||showallslicevectors==1){
              char mlabel[1024], mlabel2[1024];

              STRCPY(mlabel,sd->label.longlabel);
              if(ii>0&&sd->mesh_type!=sdm1->mesh_type){
                sprintf(mlabel2,"*** Evac type %i mesh ***",sd->mesh_type);
                glutAddMenuEntry(mlabel2,-999);
              }
              glutAddSubMenu(mlabel,loadsubvslicemenu[nloadsubvslicemenu]);
            }
            nloadsubvslicemenu++;
          }
        }
      }
    } 
    if(nvslice>0)glutAddMenuEntry("-",-999);
    if(showallslicevectors==0)glutAddMenuEntry("Show All Vector Slice Entries",-20);
    if(showallslicevectors==1)glutAddMenuEntry("*Show All Vector Slice Entries",-20);
    if(nvsliceloaded>1){
      glutAddSubMenu("Unload",unloadvslicemenu);
    }
    else{
     glutAddMenuEntry("Unload",-1);
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
        glutAddMenuEntry("-",-1);
        glutAddMenuEntry("Unload All",-10);
    }
    CREATEMENU(loadterrainmenu,LoadTerrainMenu);
    if(nterraininfo>1){
      glutAddMenuEntry("All Terrains",-9);
      glutAddMenuEntry("-",-1);
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
      glutAddMenuEntry("-",-1);
      glutAddMenuEntry("Unload Terrain",-10);
    }
    else if(nterrainloaded>1){
      glutAddMenuEntry("-",-1);
      glutAddSubMenu("Unload Terrain",unloadterrainmenu);
    }
  }
    if(nslice_files>0){

      if(nmultislices<nslice_files){
        CREATEMENU(unloadmultislicemenu,UnloadMultiSliceMenu);
        nmultisliceloaded=0;
        for(i=0;i<nmultislices;i++){
          mslicei = multisliceinfo + i;
          if(mslicei->loaded!=0){
            glutAddMenuEntry(mslicei->menulabel2,i);
          }
        }
        glutAddMenuEntry("Unload All",-1);

        nloadsubmslicemenu=1;
        for(i=1;i<nmultislices;i++){
          slice *sd, *sdim1;

          sd = sliceinfo+(multisliceinfo + i)->islices[0];
          sdim1 = sliceinfo+(multisliceinfo + i-1)->islices[0];
          if(strcmp(sd->label.longlabel,sdim1->label.longlabel)!=0)nloadsubmslicemenu++;
        }
        NewMemory((void **)&loadsubmslicemenu,nloadsubmslicemenu*sizeof(int));
        for(i=0;i<nloadsubmslicemenu;i++){
          loadsubmslicemenu[i]=0;
        }
        nloadsubmslicemenu=0;
        for(i=0;i<nmultislices;i++){
          slice *sd, *sdim1;

          sd = sliceinfo+(multisliceinfo + i)->islices[0];
          if(i>0)sdim1 = sliceinfo+(multisliceinfo + i-1)->islices[0];
          mslicei = multisliceinfo + i;
          if(i==0||strcmp(sd->label.longlabel,sdim1->label.longlabel)!=0){
            CREATEMENU(loadsubmslicemenu[nloadsubmslicemenu],LoadMultiSliceMenu);
            nloadsubmslicemenu++;
          }
          if(mslicei->loaded==1){
            STRCPY(menulabel,"*");
            STRCAT(menulabel,mslicei->menulabel);
            nmultisliceloaded++;
          }
          else if(mslicei->loaded==-1){
            STRCPY(menulabel,"#");
            STRCAT(menulabel,mslicei->menulabel);
            nmultisliceloaded++;
          }
          else{
            STRCPY(menulabel,mslicei->menulabel);
          }
          glutAddMenuEntry(menulabel,i);
        }
        CREATEMENU(loadmultislicemenu,LoadMultiSliceMenu);
        nloadsubmslicemenu=0;
        for(i=0;i<nmultislices;i++){
          slice *sd, *sdim1;

          sd = sliceinfo+(multisliceinfo + i)->islices[0];
          if(i>0)sdim1 = sliceinfo+(multisliceinfo + i-1)->islices[0];

          if(i==0||strcmp(sd->label.longlabel,sdim1->label.longlabel)!=0){
            char mlabel[1024], mlabel2[1024];

            STRCPY(mlabel,sd->label.longlabel);
            if(i>0&&sd->mesh_type!=sdim1->mesh_type){
              sprintf(mlabel2,"*** Evac type %i meshes ***",sd->mesh_type);
              glutAddMenuEntry(mlabel2,-999);
            }
            glutAddSubMenu(mlabel,             loadsubmslicemenu[nloadsubmslicemenu]);
            nloadsubmslicemenu++;
          }
        }
        if(nmultislices>0)glutAddMenuEntry("-",-999);
        if(nmultisliceloaded>1){
          glutAddSubMenu("Unload",unloadmultislicemenu);
        }
        else{
          glutAddMenuEntry("Unload",-1);
        }

      }
      CREATEMENU(unloadslicemenu,UnloadSliceMenu);
      for(i=0;i<nslice_files;i++){
        slice *sd;

        sd = sliceinfo + sliceorderindex[i];
        if(sd->loaded==1){
          STRCPY(menulabel,sd->menulabel2);  
          glutAddMenuEntry(menulabel,sliceorderindex[i]);
        }
      }
      glutAddMenuEntry("Unload All",-1);

//*** this is where I would put the "sub-slice" menus ordered by type
      nloadsubslicemenu=1;
      for(i=1;i<nslice_files;i++){
        slice *sd,*sdim1;

        sd = sliceinfo + sliceorderindex[i];
        sdim1 = sliceinfo + sliceorderindex[i-1];
        if(strcmp(sd->label.longlabel,sdim1->label.longlabel)!=0)nloadsubslicemenu++;
      }
      NewMemory((void **)&loadsubslicemenu,nloadsubslicemenu*sizeof(int));
      for(i=0;i<nloadsubslicemenu;i++){
        loadsubslicemenu[i]=0;
      }
      iloadsubslicemenu=0;
      for(i=0;i<nslice_files;i++){
        slice *sd,*sdim1;

        sd = sliceinfo + sliceorderindex[i];
        sdim1 = sliceinfo + sliceorderindex[i-1];
        if(i==0||strcmp(sd->label.longlabel,sdim1->label.longlabel)!=0){
          CREATEMENU(loadsubslicemenu[iloadsubslicemenu],LoadSliceMenu);
          iloadsubslicemenu++;
        }
        if(sd->loaded==1){
          STRCPY(menulabel,check);
          STRCAT(menulabel,sd->menulabel);
        }
        else{
          STRCPY(menulabel,sd->menulabel);
        }
        glutAddMenuEntry(menulabel,sliceorderindex[i]);
      }
      CREATEMENU(loadslicemenu,LoadSliceMenu);
      iloadsubslicemenu=0;
      for(i=0;i<nslice_files;i++){
        slice *sd,*sdim1;

        sd = sliceinfo + sliceorderindex[i];
        sdim1 = sliceinfo + sliceorderindex[i-1];
        if(i==0||strcmp(sd->label.longlabel,sdim1->label.longlabel)!=0){
          char mlabel[1024],mlabel2[1024];;

          STRCPY(mlabel,sd->label.longlabel);
          if(i>0&&sd->mesh_type!=sdim1->mesh_type){
            sprintf(mlabel2,"*** Evac type %i mesh ***",sd->mesh_type);
            glutAddMenuEntry(mlabel2,-999);
          }
          glutAddSubMenu(mlabel,loadsubslicemenu[iloadsubslicemenu]);
          iloadsubslicemenu++;
        }
      }
      glutAddMenuEntry("-",-999);
      if(nsliceloaded>1){
        glutAddSubMenu("Unload",unloadslicemenu);
      }
      else{
        glutAddMenuEntry("Unload",-1);
      }
    }

/* --------------------------------unload and load 3d smoke menus -------------------------- */

    {
      smoke3d *smoke3di;
      if(nsmoke3dloaded>0){
        CREATEMENU(unloadsmoke3dmenu,UnLoadSmoke3DMenu);
        {
          int nsootloaded=0,nhrrloaded=0,nwaterloaded=0;

          for(i=0;i<nsmoke3d_files;i++){
            smoke3di=smoke3dinfo + i;
            if(smoke3di->loaded==0)continue;
            switch (smoke3di->type){
            case 1:
              nsootloaded++;
              break;
            case 2:
              nhrrloaded++;
              break;
            case 3:
              nwaterloaded++;
              break;
            default:
              ASSERT(FFALSE);
              break;
            }
          }
          if(nsootloaded>1) glutAddMenuEntry("soot density - All meshes",-1);
          if(nhrrloaded>1)  glutAddMenuEntry("HRRPUV - All meshes",-2);
          if(nwaterloaded>1)glutAddMenuEntry("water - All meshes",-3);
          if(nsootloaded>1||nhrrloaded>1||nwaterloaded>1)glutAddMenuEntry("-",999);
        }
        for(i=0;i<nsmoke3d_files;i++){
          smoke3di=smoke3dinfo+i;
          if(smoke3di->loaded==0)continue;
          glutAddMenuEntry(smoke3di->menulabel,i);
        }
      }
    }
    {
      smoke3d *smoke3di;
      int n_soot_menu=0, n_hrr_menu=0, n_water_menu=0;

      if(nsmoke3d_files>0){
        if(nmeshes==1){
          CREATEMENU(loadsmoke3dmenu,LoadSmoke3DMenu);
        }
        if(nmeshes>1){
          CREATEMENU(loadsmoke3dsootmenu,LoadSmoke3DMenu);
        }
        for(i=0;i<nsmoke3d_files;i++){
          smoke3di = smoke3dinfo + i;
          if(smoke3di->type!=1)continue;
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
        for(i=0;i<nsmoke3d_files;i++){
          smoke3di = smoke3dinfo + i;
          if(smoke3di->type!=2)continue;
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
        for(i=0;i<nsmoke3d_files;i++){
          smoke3di = smoke3dinfo + i;
          if(smoke3di->type!=3)continue;
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
          smoke3d *smoke3dj;

          if(nmeshes>1){
            for(i=0;i<nsmoke3d_files;i++){
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
                strcat(menulabel," - All meshes");
                glutAddMenuEntry(menulabel,-useitem-10);
              }
            }
            glutAddMenuEntry("-",-2);
          }
          if(nmeshes>1){
            if(n_soot_menu>0)glutAddSubMenu("soot MASS FRACTION - Single mesh",loadsmoke3dsootmenu);
            if(n_hrr_menu>0)glutAddSubMenu("HRRPUV - Single mesh",loadsmoke3dhrrmenu);
            if(n_water_menu>0)glutAddSubMenu("Water - Single mesh",loadsmoke3dwatermenu);
          }
        }
        if(nsmoke3dloaded==1)glutAddMenuEntry("Unload",-1);
        if(nsmoke3dloaded>1)glutAddSubMenu("Unload",unloadsmoke3dmenu);
      }
    }

/* --------------------------------plot3d menu -------------------------- */

    if(nplot3d_files>0){
      plot3d *plot3dim1, *plot3di;

      CREATEMENU(unloadplot3dmenu,UnloadPlot3dMenu);
      for(ii=0;ii<nplot3d_files;ii++){
        i=plot3dorderindex[ii];
        plot3di = plot3dinfo + i;
        if(ii==0){
          strcpy(menulabel,plot3di->longlabel);
          glutAddMenuEntry(menulabel,997);
        }
        if(ii!=0&&strcmp(plot3di->longlabel,plot3dinfo[plot3dorderindex[ii-1]].longlabel)!=0){
          glutAddMenuEntry(plot3di->longlabel,997);
        }
        if(plot3di->loaded==0)continue;
        STRCPY(menulabel,plot3dinfo[i].menulabel);  
        glutAddMenuEntry(menulabel,i);
      }
      glutAddMenuEntry("Unload All",-1);


      
      nloadsubplot3dmenu=1;
      for(ii=1;ii<nplot3d_files;ii++){
        int im1;

        i = plot3dorderindex[ii];
        im1 = plot3dorderindex[ii-1];
        plot3di = plot3dinfo + i;
        plot3dim1 = plot3dinfo + im1;
        if(fabs(plot3di->time-plot3dim1->time)>0.1)nloadsubplot3dmenu++;
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
        strcat(menulabel,check);
      }
      strcat(menulabel,plot3di->menulabel);
      glutAddMenuEntry(menulabel,i);
      nloadsubplot3dmenu++;

      for(ii=1;ii<nplot3d_files;ii++){
        int im1;

        i = plot3dorderindex[ii];
        im1 = plot3dorderindex[ii-1];
        plot3di = plot3dinfo + i;
        plot3dim1 = plot3dinfo + im1;
        if(fabs(plot3di->time-plot3dim1->time)>0.1){
          if(nmeshes>1)glutAddMenuEntry("  All Meshes",-100000+nloadsubplot3dmenu-1);
          CREATEMENU(loadsubplot3dmenu[nloadsubplot3dmenu],LoadPlot3dMenu);
          nloadsubplot3dmenu++;
        }
        strcpy(menulabel,"");
        if(plot3di->loaded==1){
          strcat(menulabel,check);
        }
        strcat(menulabel,plot3di->menulabel);
        glutAddMenuEntry(menulabel,i);
      }
      if(nmeshes>1)glutAddMenuEntry("  All Meshes",-100000+nloadsubplot3dmenu-1);

      nloadsubplot3dmenu=0;
      CREATEMENU(loadplot3dmenu,LoadPlot3dMenu);
      for(ii=0;ii<nplot3d_files;ii++){
        int im1;

        i = plot3dorderindex[ii];
        plot3di = plot3dinfo + i;
        if(ii==0){
          strcpy(menulabel,plot3di->longlabel);
          glutAddMenuEntry(menulabel,997);
          sprintf(menulabel,"  %f",plot3di->time);
          trimzeros(menulabel);
          strcat(menulabel," s");
          if(nmeshes>1){
            glutAddSubMenu(menulabel,loadsubplot3dmenu[nloadsubplot3dmenu]);
          }
          else{
            strcpy(menulabel,"  ");
            if(plot3di->loaded==1){
              strcat(menulabel,check);
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
            glutAddMenuEntry(plot3di->longlabel,997);
          }
          if(fabs(plot3di->time-plot3dim1->time)>0.1){
            sprintf(menulabel,"  %f",plot3di->time);
            trimzeros(menulabel);
            strcat(menulabel," s");
            if(nmeshes>1){
              glutAddSubMenu(menulabel,loadsubplot3dmenu[nloadsubplot3dmenu]);
            }
            else{
              strcpy(menulabel,"  ");
              if(plot3di->loaded==1){
                strcat(menulabel,check);
              }
              strcat(menulabel,plot3di->menulabel);
              glutAddMenuEntry(menulabel,i);
            }
            nloadsubplot3dmenu++;
          }
        }
      }
      if(nplot3dloaded>1){
        glutAddSubMenu("Unload",unloadplot3dmenu);
      }
      else{
       glutAddMenuEntry("Unload",-1);
      }
    }

/* --------------------------------load patch menu -------------------------- */

    if(npatch_files>0){
      nloadpatchsubmenus=0;
      CREATEMENU(unloadpatchmenu,UnloadPatchMenu);
      for(ii=0;ii<npatch_files;ii++){
        patch *patchi;

        i = patchorderindex[ii];
        patchi = patchinfo + i;
        if(patchi->loaded==0)continue;
        STRCPY(menulabel,patchinfo[i].menulabel);
        glutAddMenuEntry(menulabel,i);
      }
      glutAddMenuEntry("Unload All",-1);

      if(nmeshes>1&&loadpatchsubmenus==NULL){
        NewMemory((void **)&loadpatchsubmenus,npatch_files*sizeof(int));
      }

      if(nmeshes>1){
        CREATEMENU(loadpatchsubmenus[nloadpatchsubmenus],LoadPatchMenu);
        nloadpatchsubmenus++;
      }
      else{
        CREATEMENU(loadpatchmenu,LoadPatchMenu);
      }

      for(ii=0;ii<npatch_files;ii++){
        patch *patch1, *patch2;

        i = patchorderindex[ii];
        patch2 = patchinfo + i;
        if(ii>0){
          patch1 = patchinfo + patchorderindex[ii-1];
          if(nmeshes>1&&strcmp(patch1->label.longlabel,patch2->label.longlabel)!=0){
            CREATEMENU(loadpatchsubmenus[nloadpatchsubmenus],LoadPatchMenu);
            nloadpatchsubmenus++;
          }
        }

        if(patchinfo[i].loaded==1){
          STRCPY(menulabel,check);
          STRCAT(menulabel,patchinfo[i].menulabel);  
        }
        else{
          STRCPY(menulabel,patchinfo[i].menulabel);
        }
        glutAddMenuEntry(menulabel,i);
      }

      if(nmeshes>1){
        CREATEMENU(loadpatchmenu,LoadPatchMenu);
      }

      {
        int useitem;
        patch *patchi, *patchj;

        if(nmeshes>1){
          for(i=0;i<npatch_files;i++){
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
              strcat(menulabel," - All meshes");
              glutAddMenuEntry(menulabel,-useitem-10);
            }
          }
          glutAddMenuEntry("-",-2);

          for(ii=0;ii<npatch_files;ii++){
            patch *patch1, *patch2;

            i = patchorderindex[ii];
            patch2 = patchinfo + i;
            if(ii==0){
              nloadpatchsubmenus=0;
              strcpy(menulabel,patch2->label.longlabel);
              strcat(menulabel," - Single mesh");
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

      if(npatchloaded>1){
        glutAddSubMenu("Unload",unloadpatchmenu);
      }
      else{
       glutAddMenuEntry("Unload",-1);
      }
    }

/* --------------------------------load iso menu -------------------------- */

    if(niso_files>0){
      CREATEMENU(unloadisomenu,UnloadIsoMenu);
      for(ii=0;ii<niso_files;ii++){
        iso *isoi;

        i = isoorderindex[ii];
        isoi = isoinfo + i;
        if(isoi->loaded==0)continue;
        STRCPY(menulabel,isoinfo[i].menulabel);  
        glutAddMenuEntry(menulabel,i);
      }
      glutAddMenuEntry("Unload All",-1);

      if(niso_files>0){
        if(isosubmenus==NULL){
          NewMemory((void **)&isosubmenus,niso_files*sizeof(int));
        }
        nisosubmenus=0;

        CREATEMENU(isosubmenus[nisosubmenus],LoadIsoMenu);
        nisosubmenus++;
      }

      if(nmeshes==1){
        CREATEMENU(loadisomenu,LoadIsoMenu);
      }
      for(ii=0;ii<niso_files;ii++){
        iso *iso1, *iso2;

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
          STRCPY(menulabel,check);
          STRCAT(menulabel,isoinfo[i].menulabel);  
        }
        else{
          STRCPY(menulabel,isoinfo[i].menulabel);
        }
        glutAddMenuEntry(menulabel,i);
      }

      {
        int useitem;
        iso *isoi, *isoj;

        if(nmeshes>1){
         CREATEMENU(loadisomenu,LoadIsoMenu);
          for(i=0;i<niso_files;i++){
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
              strcpy(menulabel,isoi->surface_label.longlabel);
              strcat(menulabel," - All meshes");
              glutAddMenuEntry(menulabel,-useitem-10);
            }
          }
          glutAddMenuEntry("-",-2);

          for(ii=0;ii<niso_files;ii++){
            iso *iso1, *iso2;

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
        glutAddSubMenu("Unload",unloadisomenu);
      }
      else{
       glutAddMenuEntry("Unload",-1);
      }
    }

/* --------------------------------zone menu -------------------------- */

    if(nzone>0){
      CREATEMENU(zonemenu,ZoneMenu);
      for(i=0;i<nzone;i++){
        if(zonefilenum==i){
          STRCPY(menulabel,check);
          STRCAT(menulabel,zoneinfo[i].file);  
        }
        else{STRCPY(menulabel,zoneinfo[i].file);}
        STRCAT(menulabel,", ");
        for(n=0;n<3;n++){
          STRCAT(menulabel,zoneinfo[i].label[n].shortlabel);
          STRCAT(menulabel,", ");
        }
        STRCAT(menulabel,zoneinfo[i].label[3].shortlabel);
        glutAddMenuEntry(menulabel,i);
      }
      glutAddMenuEntry("Unload",-1);

    }
/* -------------------------------- compress menu -------------------------- */

#ifdef pp_COMPRESS
  if(smokezippath!=NULL&&(npatch_files>0||nsmoke3d_files>0||nslice_files>0)){
    CREATEMENU(compressmenu,CompressMenu);
    glutAddMenuEntry("Compression Options",999);  // -c
    if(overwrite_all==1){
      glutAddMenuEntry("  *Overwrite compressed files",2);  // -f
    }
    else{
      glutAddMenuEntry("  Overwrite compressed files",2);  // -f
    }
    if(compress_autoloaded==1){
      glutAddMenuEntry("  *Compress only autoloaded files",4);  // -f
    }
    else{
      glutAddMenuEntry("  Compress only autoloaded files",4);  // -f
    }
    glutAddMenuEntry("-",999);  // -c
    glutAddMenuEntry("Compress now",3);
    glutAddMenuEntry("Erase compressed files",1);  // -c
  }
#endif


/* --------------------------------inisub menu -------------------------- */
  {
    int n_inifiles;
    inifiledata *inifile;

    n_inifiles=0;
    for(inifile=first_inifile.next;inifile->next!=NULL;inifile=inifile->next){
      if(inifile->file!=NULL&&file_exist(inifile->file)==1){
        n_inifiles++;
      }
    }
    if(n_inifiles>0){
      CREATEMENU(inisubmenu,IniSubMenu);
      if(caseinifilename!=NULL&&file_exist(caseinifilename)==1){
        glutAddMenuEntry(caseinifilename,-1);
      }
      for(inifile=first_inifile.next;inifile->next!=NULL;inifile=inifile->next){
        if(inifile->file!=NULL&&file_exist(inifile->file)==1){
          glutAddMenuEntry(inifile->file,inifile->id);
        }
      }
    }
  }

/* --------------------------------smokeviewini menu -------------------------- */

    CREATEMENU(smokeviewinimenu,SmokeviewiniMenu);


   {
    inifiledata *inifile;
    int n_inifiles;

    n_inifiles=0;
    for(inifile=first_inifile.next;inifile->next!=NULL;inifile=inifile->next){
      if(inifile->file!=NULL&&file_exist(inifile->file)==1){
        n_inifiles++;
      }
    }
    if( n_inifiles>0||file_exist(INIfile)==1||file_exist(caseinifilename)==1||file_exist(smokeviewini)==1){
      if(n_inifiles==0){
        glutAddMenuEntry("Read ini files",1);
      }
      else{
        glutAddSubMenu("Read ini files",inisubmenu);
      }
    }
  }

    glutAddMenuEntry(WRITEINIfile,2);

    STRCPY(caselabel,"Write ");
    STRCAT(caselabel,caseinifilename);

    glutAddMenuEntry(caselabel,3);

    if(ndeviceinfo>0){
      glutAddMenuEntry("-",999);
      glutAddMenuEntry("Read .svo files",4);
    }

    CREATEMENU(reloadmenu,ReloadMenu);
    glutAddMenuEntry("Reload Now",0);
    if(periodic_value==1)glutAddMenuEntry("*Reload every 1 min",1);
    if(periodic_value!=1)glutAddMenuEntry("Reload every 1 min",1);
    if(periodic_value==5)glutAddMenuEntry("*Reload every 5 min",5);
    if(periodic_value!=5)glutAddMenuEntry("Reload every 5 min",5);
    if(periodic_value==10)glutAddMenuEntry("*Reload every 10 min",10);
    if(periodic_value!=10)glutAddMenuEntry("Reload every 10 min",10);
    glutAddMenuEntry("Stop Rendering",-1);

    CREATEMENU(scriptmenu,ScriptMenu);

    if(script_recording==NULL){
      scriptfiledata *scriptfile;
      STRUCTSTAT statbuffer;
      int nscripts=0;

      nscripts=0;
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
        glutAddMenuEntry("Run Script:",-999);
        for(scriptfile=first_scriptfile.next;scriptfile->next!=NULL;scriptfile=scriptfile->next){
          char *file;
          int len;

          file=scriptfile->file;
          if(file==NULL)continue;
          len = strlen(file);
          if(len<=0)continue;
          if(STAT(file,&statbuffer)!=0)continue;

          strcpy(menulabel,"  ");
          strcat(menulabel,file);
          glutAddMenuEntry(menulabel,scriptfile->id);
        }
        glutAddMenuEntry("-",-999);
      }

    }
    glutAddMenuEntry("Create Script:",-999);
    if(script_recording==NULL)glutAddMenuEntry("  Start Recording",START_RECORDING_SCRIPT);
    glutAddMenuEntry("  Stop Recording",STOP_RECORDING_SCRIPT);

  /* --------------------------------loadunload menu -------------------------- */
    {
      char loadmenulabel[100];
      char steplabel[100];

      CREATEMENU(loadunloadmenu,LoadUnloadMenu);
#ifdef pp_OPEN
      glutAddMenuEntry("Open Smokeview (.smv) file",3);
#endif
      strcpy(steplabel,"error: steplabel not defined");
      if(nsmoke3d_files>0){
        strcpy(loadmenulabel,"3D Smoke");
        if(smoke3dframeskip>0){
          sprintf(steplabel,"/Skip %i",smoke3dframeskip);
          strcat(loadmenulabel,steplabel);
        }
        glutAddSubMenu(loadmenulabel,loadsmoke3dmenu);
      }
      if(manual_terrain==1&&nterraininfo>0){
        glutAddSubMenu("Terrain",loadterrainmenu);
      }
      if(nslice_files>0&&nmultislices<nslice_files){
        strcpy(loadmenulabel,"Multi-Slices");
        if(sliceframeskip>0){
          sprintf(steplabel,"/Skip %i",sliceframeskip);
          strcat(loadmenulabel,steplabel);
        }
        glutAddSubMenu(loadmenulabel,loadmultislicemenu);
      }
      if(nvslice>0&&nmultivslices<nvslice){
        strcpy(loadmenulabel,"Multi-Vector Slices");
        if(sliceframeskip>0){
          sprintf(steplabel,"/Skip %i",sliceframeskip);
          strcat(loadmenulabel,steplabel);
        }
        glutAddSubMenu(loadmenulabel,loadmultivslicemenu);
      }
      if(nslice_files>0){
        strcpy(loadmenulabel,"Slice File");
        if(sliceframeskip>0){
          sprintf(steplabel,"/Skip %i",sliceframeskip);
          strcat(loadmenulabel,steplabel);
        }
        glutAddSubMenu(loadmenulabel,loadslicemenu);
      }
      if(nvslice>0){
        strcpy(loadmenulabel,"Vector slices");
        if(sliceframestep>1){
          sprintf(steplabel,"/Skip %i",sliceframeskip);
          strcat(loadmenulabel,steplabel);
        }
        glutAddSubMenu(loadmenulabel,vslicemenu);
      }
      if(niso_files>0){
        strcpy(loadmenulabel,"Isosurface File");
        if(isoframeskip>0){
          sprintf(steplabel,"/Skip %i",isoframeskip);
          strcat(loadmenulabel,steplabel);
        }
        glutAddSubMenu(loadmenulabel,loadisomenu);
      }
      if(npatch_files>0){
        strcpy(loadmenulabel,"Boundary File");
        if(boundframeskip>0){
          sprintf(steplabel,"/Skip %i",boundframeskip);
          strcat(loadmenulabel,steplabel);
        }
        glutAddSubMenu(loadmenulabel,loadpatchmenu);
      }
      if(npart_files>0){
        if(nevac!=npart_files){
          strcpy(loadmenulabel,"Particle File");
          if(partframeskip>0||partpointskip>0){
            if(partframeskip>0&&partpointskip>0){
              sprintf(steplabel,"/Skip Frame %i, Point %i",partframeskip,partpointskip);
            }
            else if(partframeskip<=0&&partpointskip>0){
              sprintf(steplabel,"/Skip Point %i",partpointskip);
            }
            else if(partframeskip>0&&partpointskip<=0){
              sprintf(steplabel,"/Skip Frame %i",partframeskip);
            }
            strcat(loadmenulabel,steplabel);
          }
          glutAddSubMenu(loadmenulabel,particlemenu);
        }
        if(nevac>0){
          strcpy(loadmenulabel,"Evacuation");
          if(partframeskip>0||partpointskip>0){
            if(partframeskip>0&&partpointskip>0){
              sprintf(steplabel,"/Skip Frame %i, Point %i",partframeskip,partpointskip);
            }
            else if(partframeskip<=0&&partpointskip>0){
              sprintf(steplabel,"/Skip Point %i",partpointskip);
            }
            else if(partframeskip>0&&partpointskip<=0){
              sprintf(steplabel,"/Skip Frame %i",partframeskip);
            }
            strcat(loadmenulabel,steplabel);
          }
          glutAddSubMenu(loadmenulabel,evacmenu);
        }
      }
      if(nplot3d_files>0)glutAddSubMenu("Plot3d File",loadplot3dmenu);
      if(ntarg_files>0){
        strcpy(loadmenulabel,"Target File");
        glutAddSubMenu(loadmenulabel,targetmenu);
      }
      if(nzone>0){
        strcpy(loadmenulabel,"Zone Fire File");
        glutAddSubMenu(loadmenulabel,zonemenu);
      }
#ifdef pp_OPEN
        glutAddMenuEntry("-",999);
#else
      if(glui_active==1){      
        glutAddMenuEntry("-",999);
      }
#endif
      glutAddSubMenu("Configuration files",smokeviewinimenu);
      glutAddSubMenu("Script Options",scriptmenu);
#ifdef pp_COMPRESS
      if(smokezippath!=NULL&&(npatch_files>0||nsmoke3d_files>0||nslice_files>0)){
        glutAddSubMenu("Compression",compressmenu);
      }
#endif
#ifdef _DEBUG
      glutAddMenuEntry(smvmenufile,2);
#endif
      if(showfiles==1)glutAddMenuEntry("*Show File Names",SHOWFILES);
      if(showfiles==0)glutAddMenuEntry("Show File Names",SHOWFILES);
      glutAddSubMenu("Reload",reloadmenu);
      glutAddMenuEntry("Unload All",UNLOADALL);
    }

/* --------------------------------main menu -------------------------- */
    if(trainer_mode==1){
      CREATEMENU(trainerviewmenu,TrainerViewMenu);
      if(AnySmoke(NULL)==1){
        if(trainerload==1)glutAddMenuEntry("*Realistic",1);
        if(trainerload!=1)glutAddMenuEntry("Realistic",1);
      }
      if(AnySlices("TEMPERATURE")==1){
        if(trainerload==2)glutAddMenuEntry("*Temperature",2);
        if(trainerload!=2)glutAddMenuEntry("Temperature",2);
      }
      if(AnySlices("oxygen")==1||
         AnySlices("oxygen VOLUME FRACTION")==1){
        if(trainerload==3)glutAddMenuEntry("*Oxygen",3);
        if(trainerload!=3)glutAddMenuEntry("Oxygen",3);
      }
      glutAddMenuEntry("Clear",998);
    }

    CREATEMENU(mainmenu,MainMenu);
    if(trainer_mode==0){
      glutAddSubMenu("Load/Unload",loadunloadmenu);
      glutAddSubMenu("Show/Hide",showhidemenu);
      glutAddSubMenu("Options",optionmenu);
      glutAddSubMenu("Dialogs",dialogmenu);
      glutAddSubMenu("Tours",tourmenu);
      glutAddSubMenu("View",resetmenu);
      glutAddSubMenu("Help",helpmenu);
    }
    if(trainer_active==1){
      if(trainer_mode==1){
        glutAddMenuEntry("Smokeview Menus",997);
      }
      else{
        glutAddMenuEntry("Trainer Menus",997);
      }
    }
  glutAddMenuEntry("Quit",3);
  glutAddSubMenu("About Smokeview...",disclaimermenu);
  updatemenu=0;
#ifdef _DEBUG
  in_menu=0;
#endif

}

/* ------------------ MenuStatus ------------------------ */

void MenuStatus(int status, int x, int y){
  float *eye_xyz;
  menustatus=status;
  /* keep scene from "bouncing" around when leaving a menu */
  start_xyz0[0]=x;
  start_xyz0[1]=y;
  /*touring=0;*/
  xm0=x; ym0=y;
  eye_xyz = camera_current->eye;
  eye_xyz0[0]=eye_xyz[0];
  eye_xyz0[1]=eye_xyz[1];
  eye_xyz0[2]=eye_xyz[2];
}
