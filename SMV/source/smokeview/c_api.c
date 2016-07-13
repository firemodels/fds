
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <math.h>
#include GLUT_H

#include "update.h"
#include "smokeviewvars.h"
#include "IOvolsmoke.h"

#include "c_api.h"
#include "gd.h"

#ifdef WIN32
#define snprintf _snprintf
#endif

// function prototypes for functions drawn from other areas of smokeview
// from startup.c
void readboundini(void);
void init_lang(void);
void Init(void);
// from menus.c
void update_menu(void);
void LoadVolSmoke3DMenu(int value);
void UnLoadVolSmoke3DMenu(int value);

int set_slice_bound_min(const char *slice_type, int set, float value) {
  for (int i = 0; i < nslice2; i++) {
      printf("setting %s min bound ", slice_type);
      if(set) {printf("ON");} else {printf("OFF");}
      printf(" with value of %f\n", value);
      if(!strcmp(slice_type, slicebounds[i].datalabel)) {
          slicebounds[i].setvalmin=set;
          slicebounds[i].valmin=value;
      }
  }
  updateslicebounds();
  Slice_CB(6); // TODO: remove constant
}

float get_slice_bound_min(const char *slice_type) {
  float min, max;
  for (int i = 0; i < nslice2; i++) {
      if(!strcmp(slice_type, slicebounds[i].datalabel)) {
          min=slicebounds[i].valmin;
          // max=slicebounds[i].valmax;
      }
  }
  return min;
}

float get_slice_bound_max(const char *slice_type) {
  float min, max;
  for (int i = 0; i < nslice2; i++) {
      if(!strcmp(slice_type, slicebounds[i].datalabel)) {
          // min=slicebounds[i].valmin;
          max=slicebounds[i].valmax;
      }
  }
  return max;
}

int set_slice_bound_max(const char *slice_type, int set, float value) {
  for (int i = 0; i < nslice2; i++) {
      printf("setting %s max bound ", slice_type);
      if(set) {printf("ON");} else {printf("OFF");}
      printf(" with value of %f\n", value);
      if(!strcmp(slice_type, slicebounds[i].datalabel)) {
          slicebounds[i].setvalmax=set;
          slicebounds[i].valmax=value;
      }
  }
  updateslicebounds();
  Slice_CB(6); // TODO: remove constant
}

/* ------------------ loadsmvall ------------------------ */
// Loads a SMV file into smokeview. This should be completely independent from
// the setup of smokeview, and should be able to be run multiple times (or not
// at all). This is based on setup_case from startup.c. Currently it features
// initialisation of some GUI elements that have yet to be factored out.
int loadsmvall(const char *input_filename) {
  printf("about to load %s\n", input_filename);
  int return_code;
  // fdsprefix and input_filename_ext are global and defined in smokeviewvars.h
  // TODO: move these into the model information namespace
  parse_smv_filepath(input_filename, fdsprefix, input_filename_ext);
  printf("fdsprefix: %s\ninput_filename_ext: %s\n", fdsprefix,
         input_filename_ext);
  return_code=loadsmv(fdsprefix, input_filename_ext);
  if(return_code==0&&update_bounds==1)return_code=Update_Bounds();
  if(return_code!=0)return 1;
  // if(convert_ini==1){
    // readini(ini_from);
  // }
}

// This function takes a filepath to an smv file an finds the casename
// and the extension, which are returned in the 2nd and 3rd arguments (the
// 2nd and 3rd aguments a pre-existing strings).
int parse_smv_filepath(const char *smv_filepath, char *fdsprefix,
                       char *input_filename_ext) {
  int len_casename;
  strcpy(input_filename_ext,"");
  len_casename = (int) strlen(smv_filepath);
  if(len_casename>4){
    char *c_ext;

    c_ext=strrchr(smv_filepath,'.');
    if(c_ext!=NULL){
      STRCPY(input_filename_ext,c_ext);
      to_lower(input_filename_ext);

      if(c_ext!=NULL&&
        (strcmp(input_filename_ext,".smv")==0||
         strcmp(input_filename_ext,".svd")==0||
         strcmp(input_filename_ext,".smt")==0)
         ){
        // c_ext[0]=0;
        STRCPY(fdsprefix,smv_filepath);
        fdsprefix[strlen(fdsprefix)-4] = 0;
        strcpy(movie_name, fdsprefix);
        strcpy(render_file_base, fdsprefix);
        FREEMEMORY(trainer_filename);
        NewMemory((void **)&trainer_filename,(unsigned int)(len_casename+7));
        STRCPY(trainer_filename,smv_filepath);
        STRCAT(trainer_filename,".svd");
        FREEMEMORY(test_filename);
        NewMemory((void **)&test_filename,(unsigned int)(len_casename+7));
        STRCPY(test_filename,smv_filepath);
        STRCAT(test_filename,".smt");
      }
    }
  }
  return 0;
}

/* ------------------ loadsmv ------------------------ */
int loadsmv(char *input_filename, char *input_filename_ext){
  int return_code;
  char *input_file;
  if(input_filename==NULL){
      printf("ERROR: input_filename is NULL\n");
  }
  if(input_filename_ext==NULL){
      printf("ERROR: input_filename_ext is NULL\n");
  }
  fflush(stdout);
  printf("loading smv: %s\n", input_filename);
  printf("ext: %s\n", input_filename_ext);
  fflush(stdout);
  /*
  warning: the following line was commented out!! (perhaps it broke something)
     this line is necessary in order to define smv_filename and trainer_filename
  */
 // parse_commandlines(argc, argv);
  return_code=-1;
  if(strcmp(input_filename_ext,".svd")==0||demo_option==1){
    trainer_mode=1;
    trainer_active=1;
    if(strcmp(input_filename_ext,".svd")==0){
      input_file=trainer_filename;
    }
    else if(strcmp(input_filename_ext,".smt")==0){
      input_file=test_filename;
    }
    else{
      input_file=smv_filename;
    }
    return_code=readsmv(input_file,iso_filename);
    if(return_code==0){
      show_glui_trainer();
      show_glui_alert();
    }
  }
  else{
    input_file=smv_filename;
    return_code=readsmv(input_file,iso_filename);
  }
  switch(return_code){
    case 1:
      fprintf(stderr,"*** Error: Smokeview file, %s, not found\n",input_file);
      return 1;
    case 2:
      fprintf(stderr,"*** Error: problem reading Smokeview file, %s\n",
              input_file);
      return 2;
    case 0:
      readsmv_dynamic(input_file);
      break;
    default:
      ASSERT(FFALSE);
  }

  /* initialize units */

  InitUnits();
  init_unit_defs();
  set_unit_vis();

  CheckMemory;
  readini(NULL);
  readboundini();
  if(use_graphics==0)return 0;
#ifdef pp_LANG
  init_lang();
#endif

  if(ntours==0)setup_tour();
  glui_colorbar_setup(mainwindow_id);
  glui_motion_setup(mainwindow_id);
  glui_bounds_setup(mainwindow_id);
  glui_shooter_setup(mainwindow_id);
  glui_geometry_setup(mainwindow_id);
  glui_clip_setup(mainwindow_id);
  // glui_console_setup(mainwindow_id);
  glui_wui_setup(mainwindow_id);
  glui_labels_setup(mainwindow_id);
  glui_device_setup(mainwindow_id);
  glui_tour_setup(mainwindow_id);
  glui_alert_setup(mainwindow_id);
  glui_stereo_setup(mainwindow_id);
  glui_3dsmoke_setup(mainwindow_id);

  if(UpdateLIGHTS==1)updateLights(light_position0,light_position1);

  glutReshapeWindow(screenWidth,screenHeight);

  glutSetWindow(mainwindow_id);
  glutShowWindow();
  glutSetWindowTitle(fdsprefix);
  Init();
  glui_trainer_setup(mainwindow_id);
  glutDetachMenu(GLUT_RIGHT_BUTTON);
  InitMenus(LOAD);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
  if(trainer_mode==1){
    show_glui_trainer();
    show_glui_alert();
  }
  return 0;
}

/* ------------------ loadfile ------------------------ */

int loadfile(const char *filename) {
  int errorcode;

  FREEMEMORY(loaded_file);
  PRINTF("loading file %s\n\n",filename);
  if(filename!=NULL&&strlen(filename)>0){
    NewMemory((void **)&loaded_file,strlen(filename)+1);
    strcpy(loaded_file,filename);
  }
  for (int i=0;i<nsliceinfo;i++){
    slicedata *sd;

    sd = sliceinfo + i;
    if(strcmp(sd->file,filename)==0){
      if(i<nsliceinfo-nfedinfo){
        readslice(sd->file,i,LOAD,SET_SLICECOLOR,&errorcode);
      }
      else{
        readfed(i,LOAD,FED_SLICE,&errorcode);
      }
      return errorcode;
    }
  }
  for (int i=0;i<npatchinfo;i++){
    patchdata *patchi;

    patchi = patchinfo + i;
    if(strcmp(patchi->file,filename)==0){
      readpatch(i,LOAD,&errorcode);
      return errorcode;
    }
  }
  npartframes_max=get_min_partframes();
  for (int i=0;i<npartinfo;i++){
    partdata *parti;

    parti = partinfo + i;
    if(strcmp(parti->file,filename)==0){
      readpart(parti->file,i,LOAD,PARTDATA,&errorcode);
      return errorcode;
    }
  }
  for (int i=0;i<nisoinfo;i++){
    isodata *isoi;

    isoi = isoinfo + i;
    if(strcmp(isoi->file,filename)==0){
      readiso(isoi->file,i,LOAD,NULL,&errorcode);
      return errorcode;
    }
  }
  for (int i=0;i<nsmoke3dinfo;i++){
    smoke3ddata *smoke3di;

    smoke3di = smoke3dinfo + i;
    if(strcmp(smoke3di->file,filename)==0){
      PRINTF("loading smoke3d file:  %s\n", filename);
      readsmoke3d(i,LOAD,&errorcode);
      PRINTF("loading complete:  %d\n", errorcode);
      return errorcode;
    }
  }
  for (int i=0;i<nzoneinfo;i++){
    zonedata *zonei;

    zonei = zoneinfo + i;
    if(strcmp(zonei->file,filename)==0){
      readzone(i,LOAD,&errorcode);
      return errorcode;
    }
  }
  for (int i=0;i<nplot3dinfo;i++){
    plot3ddata *plot3di;

    plot3di = plot3dinfo + i;
    if(strcmp(plot3di->file,filename)==0){
      ReadPlot3dFile=1;
      readplot3d(plot3di->file,i,LOAD,&errorcode);
      update_menu();
      return errorcode;
    }
  }

  fprintf(stderr,"*** Error: file %s failed to load\n",filename);
  return 1;
}

/* ------------------ loadinifile ------------------------ */

void loadinifile(const char *filepath){
  PRINTF("loading ini file %s\n\n",filepath);
  windowresized=0;
  char f[1048];
  strcpy(f,filepath);
  readini(f);
}

/* ------------------ loadvfile ------------------------ */

int loadvfile(const char *filepath){
  FREEMEMORY(loaded_file);
  PRINTF("loading vector slice file %s\n\n",filepath);
  for (int i=0;i<nvsliceinfo;i++){
    slicedata *val;
    vslicedata *vslicei;

    vslicei = vsliceinfo + i;
    val = sliceinfo + vslicei->ival;
    if(val==NULL)continue;
    if(strcmp(val->reg_file,filepath)==0){
      LoadVSliceMenu(i);
      if(filepath!=NULL&&strlen(filepath)>0){
        NewMemory((void **)&loaded_file,strlen(filepath)+1);
        strcpy(loaded_file,filepath);
      }
      return 0;
    }
  }
  fprintf(stderr,"*** Error: Vector slice file %s was not loaded\n",
          filepath);
  return 1;
}

/* ------------------ loadboundaryfile ------------------------ */

void loadboundaryfile(const char *filepath){
  int errorcode;
  int count=0;

  FREEMEMORY(loaded_file);
  PRINTF("loading boundary files of type: %s\n\n",filepath);

  for (int i=0;i<npatchinfo;i++){
    patchdata *patchi;

    patchi = patchinfo + i;
    if(strcmp(patchi->label.longlabel,filepath)==0){
      LOCK_COMPRESS
      readpatch(i,LOAD,&errorcode);
      if(filepath!=NULL&&strlen(filepath)>0){
        FREEMEMORY(loaded_file);
        NewMemory((void **)&loaded_file,strlen(filepath)+1);
        strcpy(loaded_file,filepath);
      }
      count++;
      UNLOCK_COMPRESS
    }
  }
  if(count==0)fprintf(stderr,"*** Error: Boundary files of type %s failed to"
                             "load\n",filepath);
  force_redisplay=1;
  updatemenu=1;
  Update_Framenumber(0);

}

/* ------------------ label ------------------------ */

void label(const char *label) {
    PRINTF("*******************************\n");
    PRINTF("*** %s ***\n",label);
    PRINTF("*******************************\n");
}

/*
  Specify offset clip in pixels.
*/
void renderclip(int flag, int left, int right, int bottom, int top) {
    clip_rendered_scene = flag;
    render_clip_left    = left;
    render_clip_right   = right;
    render_clip_bottom  = bottom;
    render_clip_top     = top;
}

/* ------------------ render ------------------------ */

int render(const char *filename) {
    //runluascript=0;
	Display_CB();
    //runluascript=1;
    //strcpy(render_file_base,filename);
    printf("basename(c): %s\n", filename);
	return RenderFrameLua(VIEW_CENTER, filename);
}


// construct filepath for image to be renderd
char* form_filename(int view_mode, char *renderfile_name, char *renderfile_dir,
                   char *renderfile_path, int woffset, int hoffset, int screenH,
                   const char *basename) {

    // char renderfile_ext[4]; // does not include the '.'
    char suffix[20];
    char* renderfile_ext;
    char* view_suffix;

    // determine the extension to be used, and set renderfile_ext to it
    switch(renderfiletype) {
        case 0:
            renderfile_ext = ext_png;
            break;
        case 1:
            renderfile_ext = ext_jpg;
            break;
        default:
            renderfiletype = 2;
            renderfile_ext = ext_png;
            break;
    }

    // if the basename has not been specified, use a predefined method to
    // determine the filename
    if(basename == NULL) {
        printf("basename is null\n");


        switch(view_mode) {
            case VIEW_LEFT:
                if(showstereo==STEREO_LR){
                    view_suffix = "_L";
                }
                break;
            case VIEW_RIGHT:
                if(showstereo==STEREO_LR){
                    view_suffix = "_R";
                }
                break;
            case VIEW_CENTER:
                view_suffix = "";
                break;
            default:
                ASSERT(FFALSE);
                break;
        }

        if(can_write_to_dir(renderfile_dir)==0){
            printf("Creating directory: %s\n", renderfile_dir);

// #if defined(WIN32)
//             CreateDirectory (renderfile_dir, NULL);
// #elif defined(_MINGW32_)
//             CreateDirectory (renderfile_dir, NULL);
// #else
//             mkdir(renderfile_dir, 0700);
// #endif
            // TODO: ensure this can be made cross-platform
            if (strlen(renderfile_dir)>0) {
                printf("making dir: %s", renderfile_dir);
#ifdef MINGW
                mkdir(renderfile_dir);
#elif defined(pp_LINUX)
                mkdir(renderfile_dir, 0755);
#endif
#ifdef pp_OSX
                mkdir(renderfile_dir, 0755);
#endif
            }
        }
        if(showstereo==STEREO_LR
            &&(view_mode==VIEW_LEFT||view_mode==VIEW_RIGHT)){
          hoffset=screenHeight/4;
          screenH = screenHeight/2;
          if(view_mode==VIEW_RIGHT)woffset=screenWidth;
        }

        snprintf(renderfile_name, 1024,
                  "%s%s%s",
                  chidfilebase, view_suffix, renderfile_ext);
        printf("chidfilebase is: %s\n", chidfilebase);
        printf("fdsprefix is: %s\n", fdsprefix);
        printf("directory is: %s\n", renderfile_dir);
        printf("filename is: %s\n", renderfile_name);
    } else {
        snprintf(renderfile_name, 1024, "%s%s", basename, renderfile_ext);
    }
    return renderfile_name;
}

// This is function fulfills the exact same purpose as the original RenderFrame
// function, except that it takes a second argument, basename. This could be
// be used as a drop in replacement as long as all existing calls are modified
// to use basename = NULL.
  /* ------------------ RenderFrameLua ------------------------ */
// The second argument to RenderFrameLua is the name that should be given to the
// rendered file. If basename == NULL, then a default filename is formed based
// on the chosen frame and rendering options.
int RenderFrameLua(int view_mode, const char *basename) {
  char renderfile_name[1024]; // the name the file (including extension)
  char renderfile_dir[1024]; // the directory into which the image will be
                            // rendered
  char renderfile_path[2048]; // the full path of the rendered image
  int woffset=0,hoffset=0;
  int screenH;
  int return_code;

  if(script_dir_path != NULL){
      strcpy(renderfile_dir, script_dir_path);
  } else {
      strcpy(renderfile_dir, ".");
  }

#ifdef WIN32
  // reset display idle timer to prevent screen saver from activating
  SetThreadExecutionState(ES_DISPLAY_REQUIRED);
#endif

  screenH = screenHeight;
  // we should not be rendering under these conditions
  if(view_mode==VIEW_LEFT&&showstereo==STEREO_RB)return 0;
  // construct filename for image to be rendered
  form_filename(view_mode, renderfile_name, renderfile_dir, renderfile_path,
                woffset, hoffset, screenH, basename);

  printf("renderfile_name: %s\n", renderfile_name);
  // render image
  return_code = SVimage2file(renderfile_dir,renderfile_name,renderfiletype,
                             woffset,screenWidth,hoffset,screenH);
  if(RenderTime==1&&output_slicedata==1){
    output_Slicedata();
  }
  return return_code;
}

int RenderFrameLuaVar(int view_mode, gdImagePtr *RENDERimage) {
  char renderfile_name[1024]; // the name the file (including extension)
  char renderfile_dir[1024]; // the directory into which the image will be
                            // rendered
  char renderfile_path[2048]; // the full path of the rendered image
  int woffset=0,hoffset=0;
  int screenH;
  int return_code;

#ifdef WIN32
  // reset display idle timer to prevent screen saver from activating
  SetThreadExecutionState(ES_DISPLAY_REQUIRED);
#endif

  screenH = screenHeight;
  // we should not be rendering under these conditions
  if(view_mode==VIEW_LEFT&&showstereo==STEREO_RB)return 0;
  // render image
  return_code = SVimage2var(renderfiletype,
                             woffset,screenWidth,hoffset,screenH, RENDERimage);
  if(RenderTime==1&&output_slicedata==1){
    output_Slicedata();
  }
  return return_code;
}

/* ------------------ settourkeyframe ------------------------ */

void settourkeyframe(float keyframe_time) {
  keyframe *keyj,*minkey=NULL;
  tourdata *touri;
  float minkeytime=1000000000.0;

  if(selected_tour==NULL)return;
  touri=selected_tour;
  for(keyj=(touri->first_frame).next;keyj->next!=NULL;keyj=keyj->next){
    float diff_time;

    if(keyj==(touri->first_frame).next){
      minkey=keyj;
      minkeytime=ABS(keyframe_time-keyj->nodeval.time);
      continue;
    }
    diff_time=ABS(keyframe_time-keyj->nodeval.time);
    if(diff_time<minkeytime){
      minkey=keyj;
      minkeytime=diff_time;
    }
  }
  if(minkey!=NULL){
    new_select(minkey);
    set_glui_keyframe();
    update_tourcontrols();
  }
}

/* ------------------ gsliceview ------------------------ */

void gsliceview(int data, int show_triangles, int show_triangulation,
                int show_normal) {
  vis_gslice_data = data;
  show_gslice_triangles = show_triangles;
  show_gslice_triangulation = show_triangulation;
  show_gslice_normal = show_normal;
  update_gslice = 1;
}

/* ------------------ gslicepos ------------------------ */

void gslicepos(float x, float y, float z) {
  gslice_xyz[0] = x;
  gslice_xyz[1] = y;
  gslice_xyz[2] = z;
  update_gslice = 1;
}

/* ------------------ gsliceorien ------------------------ */

void gsliceorien(float az, float elev) {
  gslice_normal_azelev[0] = az;
  gslice_normal_azelev[1] = elev;
  update_gslice = 1;
}

/* ------------------ settourview ------------------------ */

void settourview(int edittourArg, int mode, int show_tourlocusArg,
                 float tour_global_tensionArg) {
  edittour=edittourArg;
  show_tourlocus=show_tourlocusArg;
  tour_global_tension_flag=1;
  tour_global_tension=tour_global_tensionArg;
  switch(mode){
    case 0:
      viewtourfrompath=0;
      keyframe_snap=0;
      break;
    case 1:
      viewtourfrompath=1;
      keyframe_snap=0;
      break;
    case 2:
      viewtourfrompath=0;
      keyframe_snap=1;
      break;
    default:
      viewtourfrompath=0;
      keyframe_snap=0;
      break;
  }
  update_tour_state();
}

int getframe() {
    int framenumber = itimes;
    return framenumber;
}

float gettime() {
    return global_times[itimes];
}
/* ------------------ settime ------------------------ */

int settime(float timeval) {
  int i,imin;
  float valmin;

  PRINTF("setting time to %f\n\n",timeval);
  if(global_times!=NULL&&nglobal_times>0){
    if(timeval<global_times[0])timeval=global_times[0];
    if(timeval>global_times[nglobal_times-1]-0.0001){
      float dt;

      dt = timeval-global_times[nglobal_times-1]-0.0001;
#ifdef pp_SETTIME
      if(nglobal_times>1&&dt>global_times[1]-global_times[0]){
        fprintf(stderr,"*** Error: data not available at time requested\n");
        fprintf(stderr,"           time: %f s, min time: %f, max time: %f s\n",
          timeval,global_times[0],global_times[nglobal_times-1]);
        if(loaded_file!=NULL)fprintf(stderr,"           loaded file: %s\n",
                                     loaded_file);
        if(script_labelstring!=NULL)fprintf(stderr,"                 "
                                            "label: %s\n",script_labelstring);
      }
#endif
      timeval=global_times[nglobal_times-1]-0.0001;
    }
    valmin=ABS(global_times[0]-timeval);
    imin=0;
    for(i=1;i<nglobal_times-1;i++){
      float val;

      val = ABS(global_times[i]-timeval);
      if(val<valmin){
        valmin=val;
        imin=i;
      }
    }
    itimes=imin;
    script_itime=imin;
    stept=0;
    force_redisplay=1;
    Update_Framenumber(0);
    UpdateTimeLabels();
    PRINTF("script: time set to %f s (i.e. frame number: %d)\n\n",
         global_times[itimes], itimes);
    return 0;
  } else {
    PRINTF("no data loaded, time could not be set\n");
    return 1;
  }

}

void set_slice_in_obst(int setting) {
	show_slice_in_obst = setting;
	if(show_slice_in_obst==0)PRINTF("Not showing slices witin blockages.\n");
  if(show_slice_in_obst==1)PRINTF("Showing slices within blockages.\n");
  // updateslicefilenum();
  // plotstate=getplotstate(DYNAMIC_PLOTS);
	//
  // updateglui();
  // updateslicelistindex(slicefilenum);
  // Update_Show();
}

int get_slice_in_obst() {
	return show_slice_in_obst;
}

// colorbar visibility
void set_colorbar_visibility(int setting) {
  visColorbar = setting;
  if(visColorbar==0)PRINTF("Colorbar hidden\n");
  if(visColorbar==1)PRINTF("Colorbar visible\n");
}

int get_colorbar_visibility() {
  return visColorbar;
}

void toggle_colorbar_visibility() {
  visColorbar = 1 - visColorbar;
  if(visColorbar==0)PRINTF("Colorbar hidden\n");
  if(visColorbar==1)PRINTF("Colorbar visible\n");
}

// timebar visibility
void set_timebar_visibility(int setting) {
  visTimebar = setting;
  if(visTimebar==0)PRINTF("Time bar hidden\n");
  if(visTimebar==1)PRINTF("Time bar visible\n");
}

int get_timebar_visibility() {
  return visTimebar;
}

void toggle_timebar_visibility() {
  visTimebar = 1 - visTimebar;
  if(visTimebar==0)PRINTF("Time bar hidden\n");
  if(visTimebar==1)PRINTF("Time bar visible\n");
}

// title visibility
void set_title_visibility(int setting) {
  visTitle = setting;
  if(visTitle==0)PRINTF("Title hidden\n");
  if(visTitle==1)PRINTF("Title visible\n");
}

int get_title_visibility() {
  return visTitle;
}

void toggle_title_visibility() {
  visTitle = 1 - visTitle;
  if(visTitle==0)PRINTF("Title hidden\n");
  if(visTitle==1)PRINTF("Title visible\n");
}

// axis visibility
void set_axis_visibility(int setting) {
  visaxislabels = setting;
  update_visaxislabels();
  if(visaxislabels==0)PRINTF("Axis labels hidden\n");
  if(visaxislabels==1)PRINTF("Axis labels visible\n");
}

int get_axis_visibility() {
  return visaxislabels;
}

void toggle_axis_visibility() {
  visaxislabels = 1 - visaxislabels;
  update_visaxislabels();
  if(visaxislabels==0)PRINTF("Axis labels hidden\n");
  if(visaxislabels==1)PRINTF("Axis labels visible\n");
}

// framelabel visibility
void set_framelabel_visibility(int setting) {
  visFramelabel = setting;
  // The frame label should not be shown without the timebar
  // so show timebar if necessary.
  if(visFramelabel==1)visTimebar=1;
  if(visFramelabel==1){
    visHRRlabel=0;
    if(hrrinfo!=NULL){
      hrrinfo->display=visHRRlabel;
      Update_Times();
    }
  }
  if(visFramelabel==0)PRINTF("Frame label hidden\n");
  if(visFramelabel==1)PRINTF("Frame label visible\n");
}

int get_framelabel_visibility() {
  return visFramelabel;
}

void toggle_framelabel_visibility() {
  visFramelabel = 1 - visFramelabel;
  // The frame label should not be shown without the timebar
  // so show timebar if necessary.
  if(visFramelabel==1)visTimebar=1;
    if(visFramelabel==1){
      visHRRlabel=0;
      if(hrrinfo!=NULL){
        hrrinfo->display=visHRRlabel;
        Update_Times();
      }
  }
  if(visFramelabel==0)PRINTF("Frame label hidden\n");
  if(visFramelabel==1)PRINTF("Frame label visible\n");
}

// framerate visibility
void set_framerate_visibility(int setting) {
  visFramerate = setting;
  if(visFramerate==0)PRINTF("Frame rate hidden\n");
  if(visFramerate==1)PRINTF("Frame rate visible\n");
}

int get_framerate_visibility() {
  return visFramerate;
}

void toggle_framerate_visibility() {
  visFramerate = 1 - visFramerate;
  if(visFramerate==0)PRINTF("Frame rate hidden\n");
  if(visFramerate==1)PRINTF("Frame rate visible\n");
}

// grid locations visibility
void set_gridloc_visibility(int setting) {
  visgridloc = setting;
  if(visgridloc==0)PRINTF("Grid locations hidden\n");
  if(visgridloc==1)PRINTF("Grid locations visible\n");
}

int get_gridloc_visibility() {
  return visgridloc;
}

void toggle_gridloc_visibility() {
  visgridloc = 1 - visgridloc;
  if(visgridloc==0)PRINTF("Grid locations hidden\n");
  if(visgridloc==1)PRINTF("Grid locations visible\n");
}

// HRRPUV cutoff visibility
void set_hrrcutoff_visibility(int setting) {
  show_hrrcutoff = setting;
  if(show_hrrcutoff==0)PRINTF("HRR cutoff hidden\n");
  if(show_hrrcutoff==1)PRINTF("HRR cutoff visible\n");
}

int get_hrrcutoff_visibility() {
  return show_hrrcutoff;
}

void toggle_hrrcutoff_visibility() {
  show_hrrcutoff = 1 - show_hrrcutoff;
  if(show_hrrcutoff==0)PRINTF("HRR cutoff hidden\n");
  if(show_hrrcutoff==1)PRINTF("HRR cutoff visible\n");
}

// HRR label
void set_hrrlabel_visibility(int setting) {
  visHRRlabel = setting;
  if (hrrinfo != NULL&&hrrinfo->display != 0)Update_hrrinfo(0);
  if(show_hrrcutoff==0)PRINTF("HRR label hidden\n");
  if(show_hrrcutoff==1)PRINTF("HRR label visible\n");
}

int get_hrrlabel_visibility() {
  return visHRRlabel;
}

void toggle_hrrlabel_visibility() {
  visHRRlabel = 1 - visHRRlabel;
  if (hrrinfo != NULL&&hrrinfo->display != 0)Update_hrrinfo(0);
  if(show_hrrcutoff==0)PRINTF("HRR label hidden\n");
  if(show_hrrcutoff==1)PRINTF("HRR label visible\n");
}

// memory load
#ifdef pp_memstatus
void set_memload_visibility(int setting) {
  visAvailmemory = setting;
  if(visAvailmemory==0)PRINTF("Memory load hidden\n");
  if(visAvailmemory==1)PRINTF("Memory load visible\n");
}

int get_memload_visibility() {
  return visAvailmemory;
}

void toggle_memload_visibility() {
  visAvailmemory = 1 - visAvailmemory;
  if(visAvailmemory==0)PRINTF("Memory load hidden\n");
  if(visAvailmemory==1)PRINTF("Memory load visible\n");
}
#endif

// mesh label
void set_meshlabel_visibility(int setting) {
  visMeshlabel = setting;
  if(visMeshlabel==0)PRINTF("Mesh label hidden\n");
  if(visMeshlabel==1)PRINTF("Mesh label visible\n");
}

int get_meshlabel_visibility() {
  return visMeshlabel;
}

void toggle_meshlabel_visibility() {
  visMeshlabel = 1 - visMeshlabel;
  if(visMeshlabel==0)PRINTF("Mesh label hidden\n");
  if(visMeshlabel==1)PRINTF("Mesh label visible\n");
}

// slice average
void set_slice_average_visibility(int setting) {
  vis_slice_average = setting;
  if(vis_slice_average==0)PRINTF("Slice average hidden\n");
  if(vis_slice_average==1)PRINTF("Slice average visible\n");
}

int get_slice_average_visibility() {
  return vis_slice_average;
}

void toggle_slice_average_visibility() {
  vis_slice_average = 1 - vis_slice_average;
  if(vis_slice_average==0)PRINTF("Slice average hidden\n");
  if(vis_slice_average==1)PRINTF("Slice average visible\n");
}

// time
void set_time_visibility(int setting) {
  visTimelabel = setting;
  if(visTimelabel==0)PRINTF("Time label hidden\n");
  if(visTimelabel==1)PRINTF("Time label visible\n");
}

int get_time_visibility() {
  return visTimelabel;
}

void toggle_time_visibility() {
  visTimelabel = 1 - visTimelabel;
  if(visTimelabel==0)PRINTF("Time label hidden\n");
  if(visTimelabel==1)PRINTF("Time label visible\n");
}

// user settable ticks
void set_user_ticks_visibility(int setting) {
  visUSERticks = setting;
  if(visUSERticks==0)PRINTF("User settable ticks hidden\n");
  if(visUSERticks==1)PRINTF("User settable ticks visible\n");
}

int get_user_ticks_visibility() {
  return visUSERticks;
}

void toggle_user_ticks_visibility() {
  visUSERticks = 1 - visUSERticks;
  if(visUSERticks==0)PRINTF("User settable ticks hidden\n");
  if(visUSERticks==1)PRINTF("User settable ticks visible\n");
}

//version info
void set_version_info_visibility(int setting) {
  gversion = setting;
  if(gversion==0)PRINTF("Version info hidden\n");
  if(gversion==1)PRINTF("Version info visible\n");
}

int get_version_info_visibility() {
  return gversion;
}

void toggle_version_info_visibility() {
  gversion = 1 - gversion;
  if(gversion==0)PRINTF("Version info hidden\n");
  if(gversion==1)PRINTF("Version info visible\n");
}

void set_all_label_visibility(int setting) {
  set_colorbar_visibility(setting);
  set_timebar_visibility(setting);
  set_title_visibility(setting);
  set_axis_visibility(setting);
  set_framelabel_visibility(setting);
  set_framerate_visibility(setting);
  set_gridloc_visibility(setting);
  set_hrrcutoff_visibility(setting);
  set_hrrlabel_visibility(setting);
  #ifdef pp_memstatus
  set_memload_visibility(setting);
  #endif
  set_meshlabel_visibility(setting);
  set_slice_average_visibility(setting);
  set_time_visibility(setting);
  set_user_ticks_visibility(setting);
  set_version_info_visibility(setting);
}

// Display Units
// time
void set_timehms(int setting) {
  vishmsTimelabel = 1 - vishmsTimelabel;
  set_labels_controls();
  if(vishmsTimelabel==0)PRINTF("Time label in h:m:s\n");
  if(vishmsTimelabel==1)PRINTF("Time label in s\n");
}

int get_timehms() {
  return vishmsTimelabel;
}

void toggle_timehms() {
  vishmsTimelabel = 1 - vishmsTimelabel;
  set_labels_controls();
  if(vishmsTimelabel==0)PRINTF("Time label in h:m:s\n");
  if(vishmsTimelabel==1)PRINTF("Time label in s\n");
}

// arbitrary
void set_units(int unitclass, int unit_index) {
  unitclasses[unitclass].unit_index=unit_index;
  updatemenu=1;
  glutPostRedisplay();
}

void set_units_default() {
  for (int i=0;i<nunitclasses;i++){
      unitclasses[i].unit_index=0;
    }
  updatemenu=1;
  glutPostRedisplay();
}

void set_unitclass_default(int unitclass) {
  unitclasses[unitclass].unit_index=0;
  updatemenu=1;
  glutPostRedisplay();
}

// Show/Hide Geometry
// Obstacles
// View Method

/** Set the method for viewing blockages.
 * The 'setting' integer is used as following:
 * 0 - Defined in input file
 * 1 - Solid
 * 2 - Outine only
 * 3 - Outline added
 * 4 - Hidden
 */
int blockage_view_method(int setting) {
  int value;
  switch(setting) {
    case 0:
      value=visBLOCKAsInput;
      break;
    case 1:
      value=visBLOCKNormal;
      break;
    case 2:
      value=visBLOCKOutline;
      break;
    case 3:
      value=visBLOCKAddOutline;
      break;
    case 4:
      value=visBLOCKHide;
      break;
    default:
      return 1;
      break;
  }
  // TODO
  // The below is the menu code verbatim. Simplify to contain only the
  // necessary code.
  BlockageMenu(value);
  return 0;
}

/** Set the color to be used when drawing blockage outlines.
 * The 'setting' integer is used as following:
 * 0 - Use blockage
 * 1 - Use foreground
 */
int blockage_outline_color(int setting) {
  switch(setting) {
    case 0:
      outline_color_flag = 0;
      updatefaces=1;
      break;
    case 1:
      outline_color_flag = 1;
      updatefaces=1;
      break;
    default:
      return 1;
      break;
  }
  return 0;
}

/** Determine how the blockages should be displayed.
 * The 'setting' integer is used as following:
 * 0 - grid - Snapped to the grid as used by FDS.
 * 1 - exact - As specified.
 * 2 - cad - Using CAD geometry.
 * This is used for the BLOCKLOCATION .ini option.
 */
int blockage_locations(int setting) {
  switch(setting) {
    case 0:
      blocklocation=BLOCKlocation_grid;
      break;
    case 1:
      blocklocation=BLOCKlocation_exact;
    case 2:
      blocklocation=BLOCKlocation_cad;
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
  return 0;
}

void setframe(int framenumber) {
  itimes=framenumber;
  script_itime=itimes;
  stept=0;
  force_redisplay=1;
  Update_Framenumber(0);
  UpdateTimeLabels();
}

/* ------------------ loadvolsmoke ------------------------ */
/*
  Load files needed to view volume rendered smoke. One may either load files for
  all meshes or for one particular mesh. Use meshnumber = -1 for all meshes.
*/
void loadvolsmoke(int meshnumber) {
  int imesh;

  imesh = meshnumber;
  if(imesh==-1){
    read_vol_mesh=VOL_READALL;
    read_volsmoke_allframes_allmeshes2(NULL);
  }
  else if(imesh>=0&&imesh<nmeshes){
    meshdata *meshi;
    volrenderdata *vr;

    meshi = meshinfo + imesh;
    vr = &meshi->volrenderinfo;
    read_volsmoke_allframes(vr);
  }
}

/* ------------------ loadvolsmokeframe ------------------------ */
/*
  As with loadvolsmoke, but for a single frame indicated my framenumber. Flag is
  set to 1 when calling from a script. Reason unkown.
*/
void loadvolsmokeframe(int meshnumber, int framenumber, int flag) {
  int framenum, index;
  int first = 1;
  int i;

  index = meshnumber;
  framenum = framenumber;
  if(index > nmeshes - 1)index = -1;
  for (i = 0; i < nmeshes; i++){
    if(index == i || index < 0){
      meshdata *meshi;
      volrenderdata *vr;

      meshi = meshinfo + i;
      vr = &meshi->volrenderinfo;
      free_volsmoke_frame(vr, framenum);
      read_volsmoke_frame(vr, framenum, &first);
      if(vr->times_defined == 0){
        vr->times_defined = 1;
        get_volsmoke_all_times(vr);
      }
      vr->loaded = 1;
      vr->display = 1;
    }
  }
  plotstate = getplotstate(DYNAMIC_PLOTS);
  stept = 1;
  Update_Times();
  force_redisplay = 1;
  Update_Framenumber(framenum);
  i = framenum;
  itimes = i;
  script_itime = i;
  stept = 1;
  force_redisplay = 1;
  Update_Framenumber(0);
  UpdateTimeLabels();
  // TODO: replace with a call to render()
  keyboard('r', FROM_SMOKEVIEW);
  RenderOnceNow = 0;
  if(flag == 1)script_render = 1;// called when only rendering a single frame
}

/* ------------------ load3dsmoke ------------------------ */

void load3dsmoke(const char *smoke_type){
  int errorcode;
  int count=0;

  FREEMEMORY(loaded_file);
  PRINTF("script: loading smoke3d files of type: %s\n\n",smoke_type);

  for (int i=0;i<nsmoke3dinfo;i++){
    smoke3ddata *smoke3di;

    smoke3di = smoke3dinfo + i;
    if(match_upper(smoke3di->label.longlabel,smoke_type)==1){
      readsmoke3d(i,LOAD,&errorcode);
      if(smoke_type!=NULL&&strlen(smoke_type)>0){
        FREEMEMORY(loaded_file);
        NewMemory((void **)&loaded_file,strlen(smoke_type)+1);
        strcpy(loaded_file,smoke_type);
      }
      count++;
    }
  }
  if(count==0)fprintf(stderr,"*** Error: Smoke3d files of type %s failed to "
                      "load\n",smoke_type);
  force_redisplay=1;
  updatemenu=1;

}

void rendertype(const char *type) {
    if(STRCMP(type, "JPG")==0){
      update_render_type(JPEG);
    }
    else{
      update_render_type(PNG);
    }
}

int get_rendertype() {
    return renderfiletype;
}

void set_movietype(const char *type) {
    if(STRCMP(type, "WMV") == 0){
      update_movie_type(WMV);
    }
    if(STRCMP(type, "MP4") == 0){
      update_movie_type(MP4);
    }
    else{
      update_movie_type(AVI);
    }
}

int get_movietype() {
    return moviefiletype;
}

void makemovie(const char *name, const char *base, float framerate) {
    strcpy(movie_name, name);
    strcpy(render_file_base, base);
    movie_framerate=framerate;
    Render_CB(MAKE_MOVIE);
}

/* ------------------ script_loadtour ------------------------ */

int loadtour(const char *tourname) {
  int count=0;
  int errorcode = 0;
  PRINTF("loading tour %s\n\n", tourname);

  for (int i=0;i<ntours;i++){
    tourdata *touri;

    touri = tourinfo + i;
    if(strcmp(touri->label,tourname)==0){
      TourMenu(i);
      viewtourfrompath=0;
      TourMenu(MENU_TOUR_VIEWFROMROUTE);
      count++;
      break;
    }
  }

  if(count==0) {
      fprintf(stderr,"*** Error: The tour %s failed to load\n",
                      tourname);
      errorcode = 1;
  }
  force_redisplay=1;
  updatemenu=1;
  return errorcode;
}

/* ------------------ loadparticles ------------------------ */

void loadparticles(const char *name){
  int errorcode;
  int count=0;

  FREEMEMORY(loaded_file);

  PRINTF("script: loading particles files\n\n");

  npartframes_max=get_min_partframes();
  for (int i=0;i<npartinfo;i++){
    partdata *parti;

    parti = partinfo + i;
    if(parti->evac==1)continue;
    readpart(parti->file,i,UNLOAD,PARTDATA,&errorcode);
    count++;
  }
  for (int i=0;i<npartinfo;i++){
    partdata *parti;

    parti = partinfo + i;
    if(parti->evac==1)continue;
      readpart(parti->file,i,LOAD,PARTDATA,&errorcode);
      if(name!=NULL&&strlen(name)>0){
        FREEMEMORY(loaded_file);
        NewMemory((void **)&loaded_file,strlen(name)+1);
        strcpy(loaded_file,name);
      }
      count++;
  }
  if(count==0)fprintf(stderr,"*** Error: Particles files failed to load\n");
  force_redisplay=1;
  Update_Framenumber(0);
  updatemenu=1;
}
/* ------------------ partclasscolor ------------------------ */

void partclasscolor(const char *color){
  int count=0;

  for (int i=0;i<npart5prop;i++){
	partpropdata *propi;

	propi = part5propinfo + i;
	if(propi->particle_property==0)continue;
	if(strcmp(propi->label->longlabel,color)==0){
	  ParticlePropShowMenu(i);
	  count++;
	}
  }
  if(count==0)fprintf(stderr,
      "*** Error: particle class color: %s failed to be set\n",color);
}

/* ------------------ partclasstype ------------------------ */

void partclasstype(const char *part_type){
  int count=0;

  for (int i=0;i<npart5prop;i++){
    partpropdata *propi;
    int j;

    propi = part5propinfo + i;
    if(propi->display==0)continue;
    for (int j=0;j<npartclassinfo;j++){
      partclassdata *partclassj;

      if(propi->class_present[j]==0)continue;
      partclassj = partclassinfo + j;
      if(partclassj->kind==HUMANS)continue;
      if(strcmp(partclassj->name,part_type)==0){
        ParticlePropShowMenu(-10-j);
        count++;
      }
    }
  }
  if(count==0)fprintf(stderr,"*** Error: particle class type %s failed to be "
                             "set\n",part_type);
}

/* ------------------ script_plot3dprops ------------------------ */

void plot3dprops(int variable_index, int showvector, int vector_length_index,
                 int display_type, float vector_length) {
  int p_index;

  p_index = variable_index;
  if(p_index<1)p_index=1;
  if(p_index>5)p_index=5;

  visVector = showvector;
  if(visVector!=1)visVector=0;

  plotn = p_index;
  if(plotn<1){
    plotn=numplot3dvars;
  }
  if(plotn>numplot3dvars){
    plotn=1;
  }
  updateallplotslices();
  if(visiso==1)updatesurface();
  updateplot3dlistindex();

  vecfactor=1.0;
  if(vector_length>=0.0)vecfactor=vector_length;
  update_vector_widgets();

  PRINTF("vecfactor=%f\n",vecfactor);

  contour_type=CLAMP(display_type,0,2);
  update_plot3d_display();

  if(visVector==1&&ReadPlot3dFile==1){
    meshdata *gbsave,*gbi;

    gbsave=current_mesh;
    for (int i=0;i<nmeshes;i++){
      gbi = meshinfo + i;
      if(gbi->plot3dfilenum==-1)continue;
      update_current_mesh(gbi);
      updateplotslice(XDIR);
      updateplotslice(YDIR);
      updateplotslice(ZDIR);
    }
    update_current_mesh(gbsave);
  }
}

// TODO: this function uses 5 int script values, but the documentation only
// lists three. Find out what is going on here.
// /* ------------------ script_showplot3ddata ------------------------ */
//
// void script_showplot3ddata(int meshnumber, int plane_orientation, int display,
//                            float position) {
//   meshdata *meshi;
//   int imesh, dir, showhide;
//   float val;
//   int isolevel;
//
//   imesh = scripti->ival-1;
//   if(imesh<0||imesh>nmeshes-1)return;
//
//   meshi = meshinfo + imesh;
//   update_current_mesh(meshi);
//
//   dir = scripti->ival2;
//   if(dir<1)dir=1;
//   if(dir>4)dir=4;
//
//   plotn=scripti->ival3;
//
//   showhide = scripti->ival4;
//   val = scripti->fval;
//
//   switch(dir){
//     case 1:
//       visx_all=showhide;
//       iplotx_all=get_index(val,1,plotx_all,nplotx_all);
//       next_xindex(1,0);
//       next_xindex(-1,0);
//       break;
//     case 2:
//       visy_all=showhide;
//       iploty_all=get_index(val,2,ploty_all,nploty_all);
//       next_yindex(1,0);
//       next_yindex(-1,0);
//       break;
//     case 3:
//       visz_all=showhide;
//       iplotz_all=get_index(val,3,plotz_all,nplotz_all);
//       next_zindex(1,0);
//       next_zindex(-1,0);
//       break;
//     case 4:
//       isolevel=scripti->ival5;
//       plotiso[plotn-1]=isolevel;
//       updateshowstep(showhide,ISO);
//       updatesurface();
//       updatemenu=1;
//       break;
//     default:
//       ASSERT(FFALSE);
//       break;
//   }
//   updateplotslice(dir);
//
// }


/* ------------------ loadplot3d ------------------------ */

void loadplot3d(int meshnumber, float time_local){
  int blocknum;
  int count=0;

  time_local = time_local;
  blocknum = meshnumber-1;

  for(int i=0;i<nplot3dinfo;i++){
    plot3ddata *plot3di;

    plot3di = plot3dinfo + i;
    if(plot3di->blocknumber==blocknum&&ABS(plot3di->time-time_local)<0.5){
      count++;
      LoadPlot3dMenu(i);
    }
  }
  UpdateRGBColors(COLORBAR_INDEX_NONE);
  set_labels_controls();
  if(count==0)fprintf(stderr,"*** Error: Plot3d file failed to load\n");

  //update_menu();
}

/* ------------------ loadiso ------------------------ */

void loadiso(const char *type) {
  int count=0;

  FREEMEMORY(loaded_file);
  PRINTF("loading isosurface files of type: %s\n\n",type);

  update_readiso_geom_wrapup = UPDATE_ISO_START_ALL;
  for (int i = 0; i<nisoinfo; i++){
    int errorcode;
    isodata *isoi;

    isoi = isoinfo + i;
    if(STRCMP(isoi->surface_label.longlabel,type)==0){
      readiso(isoi->file,i,LOAD,NULL,&errorcode);
      if(type != NULL&&strlen(type)>0){
        FREEMEMORY(loaded_file);
        NewMemory((void **)&loaded_file,strlen(type)+1);
        strcpy(loaded_file,type);
      }
      count++;
    }
  }
  if(update_readiso_geom_wrapup == UPDATE_ISO_ALL_NOW)readiso_geom_wrapup();
  update_readiso_geom_wrapup = UPDATE_ISO_OFF;
  if(count == 0)fprintf(stderr, "*** Error: Isosurface files of type %s failed "
                                "to load\n", type);
  force_redisplay=1;
  updatemenu=1;
}

/* ------------------ loadslice ------------------------ */

void loadslice(const char *type, int axis, float distance){
  int count=0;

  PRINTF("loading slice files of type: %s\n\n",type);

  for (int i=0;i<nmultisliceinfo;i++){
    multislicedata *mslicei;
    slicedata *slicei;
    float delta_orig;

    mslicei = multisliceinfo + i;
    if(mslicei->nslices<=0)continue;
    slicei = sliceinfo + mslicei->islices[0];
    printf("slice name: %s, axis: %d, position: %f\n", slicei->label.longlabel,
           slicei->idir, slicei->position_orig);
    if(match_upper(slicei->label.longlabel,type)==0)continue;
    if(slicei->idir!=axis)continue;
    delta_orig = slicei->position_orig - distance;
    if(delta_orig<0.0)delta_orig = -delta_orig;
    if(delta_orig>slicei->delta_orig)continue;

    for (int j=0;j<mslicei->nslices;j++){
      LoadSliceMenu(mslicei->islices[j]);
      count++;
    }
    break;
  }
  if(count==0)fprintf(stderr,"*** Error: Slice files of type %s failed to "
                             "load\n",type);
}

/* ------------------ loadvslice ------------------------ */

void loadvslice(const char *type, int axis, float distance){
  float delta_orig;
  int count=0;

  PRINTF("loading vector slice files of type: %s\n\n",type);

  for (int i=0;i<nmultivsliceinfo;i++){
    multivslicedata *mvslicei;
    slicedata *slicei;

    mvslicei = multivsliceinfo + i;
    if(mvslicei->nvslices<=0)continue;
    slicei = sliceinfo + mvslicei->ivslices[0];
    if(match_upper(slicei->label.longlabel,type)==0)continue;
    if(slicei->idir!=axis)continue;
    delta_orig = slicei->position_orig - distance;
    if(delta_orig<0.0)delta_orig = -delta_orig;
    if(delta_orig>slicei->delta_orig)continue;

    for (int j=0;j<mvslicei->nvslices;j++){
      LoadVSliceMenu(mvslicei->ivslices[j]);
      count++;
    }
    break;
  }
  if(count==0)fprintf(stderr,"*** Error: Vector slice files of type %s failed "
                             "to load\n",type);
}

/* ------------------ unloadslice ------------------------ */

void unloadslice(int value){
  int errorcode;

  updatemenu=1;
  glutPostRedisplay();
  if(value>=0){
    readslice("",value,UNLOAD,SET_SLICECOLOR,&errorcode);
  }
  else{
    if(value==UNLOAD_ALL){
      for (int i=0;i<nsliceinfo;i++){
        readslice("",i,UNLOAD,SET_SLICECOLOR,&errorcode);
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

/* ------------------ unloadall ------------------------ */

void unloadall() {
  int errorcode;
  int ii;

   // leaving code here commented in case I later decide to unload terrain files
   // for(i=0;i<nterraininfo;i++){
   //   readterrain("",i,UNLOAD,&errorcode);
   // }
#ifdef pp_DEVICE
    if(devc_csv_filename!=NULL){
      read_device_data(devc_csv_filename,CSV_FDS,UNLOAD);
    }
    if(exp_csv_filename!=NULL){
      read_device_data(exp_csv_filename,CSV_EXP,UNLOAD);
    }
#endif
    if(nvolrenderinfo>0){
      LoadVolSmoke3DMenu(UNLOAD_ALL);
    }
    for (int i=0;i<nsliceinfo;i++){
      readslice("",i,UNLOAD,SET_SLICECOLOR,&errorcode);
    }
    for (int i=0;i<nplot3dinfo;i++){
      readplot3d("",i,UNLOAD,&errorcode);
    }
    for (int i=0;i<npatchinfo;i++){
      readpatch(i,UNLOAD,&errorcode);
    }
    for (int i=0;i<npartinfo;i++){
      readpart("",i,UNLOAD,PARTDATA,&errorcode);
    }
    for (int i=0;i<nisoinfo;i++){
      readiso("",i,UNLOAD,NULL,&errorcode);
    }
    for (int i=0;i<nzoneinfo;i++){
      readzone(i,UNLOAD,&errorcode);
    }
    for (int i=0;i<nsmoke3dinfo;i++){
      readsmoke3d(i,UNLOAD,&errorcode);
    }
    if(nvolrenderinfo>0){
      UnLoadVolSmoke3DMenu(UNLOAD_ALL);
    }
}

void unloadtour() {
    TourMenu(MENU_TOUR_MANUAL);
}

/* ------------------ exit_smokeview ------------------------ */

void exit_smokeview() {
	PRINTF("exiting...\n");
	exit(EXIT_SUCCESS);
}

/* ------------------ setviewpoint ------------------------ */

int setviewpoint(const char *viewpoint){
  cameradata *ca;
  int count=0;
  int errorcode = 0;
  PRINTF("setting viewpoint to %s\n\n",viewpoint);
  for(ca=camera_list_first.next;ca->next!=NULL;ca=ca->next){
    if(strcmp(viewpoint,ca->name)==0){
      ResetMenu(ca->view_id);
      count++;
      break;
    }
  }
  if(count==0){
      errorcode = 1;
      fprintf(stderr,"*** Error: The viewpoint %s was not found\n",viewpoint);
  }
  fprintf(stderr, "Viewpoint set to %s\n", camera_current->name);
  return errorcode;
}

int get_clipping_mode() {
  return clip_mode;
}

void set_clipping_mode(int mode) {
    clip_mode=mode;
    updatefacelists=1;
}

void set_sceneclip_x(int clipMin, float min, int clipMax, float max) {
    clipinfo.clip_xmin=clipMin;
    clipinfo.xmin = min;

    clipinfo.clip_xmax=clipMax;
    clipinfo.xmax = max;
    updatefacelists=1;
}

void set_sceneclip_x_min(int flag, float value) {
    clipinfo.clip_xmin = flag;
    clipinfo.xmin = value;
    updatefacelists = 1;
}

void set_sceneclip_x_max(int flag, float value) {
    clipinfo.clip_xmax = flag;
    clipinfo.xmax = value;
    updatefacelists = 1;
}

void set_sceneclip_y(int clipMin, float min, int clipMax, float max) {
    clipinfo.clip_ymin=clipMin;
    clipinfo.ymin = min;

    clipinfo.clip_ymax=clipMax;
    clipinfo.ymax = max;
    updatefacelists=1;
}

void set_sceneclip_y_min(int flag, float value) {
    clipinfo.clip_ymin = flag;
    clipinfo.ymin = value;
    updatefacelists = 1;
}

void set_sceneclip_y_max(int flag, float value) {
    clipinfo.clip_ymax = flag;
    clipinfo.ymax = value;
    updatefacelists = 1;
}

void set_sceneclip_z(int clipMin, float min, int clipMax, float max) {
    clipinfo.clip_zmin=clipMin;
    clipinfo.zmin = min;

    clipinfo.clip_zmax=clipMax;
    clipinfo.zmax = max;
    updatefacelists=1;
    update_clip_all();
}

void set_sceneclip_z_min(int flag, float value) {
    clipinfo.clip_zmin = flag;
    clipinfo.zmin = value;
    updatefacelists = 1;
}

void set_sceneclip_z_max(int flag, float value) {
    clipinfo.clip_zmax = flag;
    clipinfo.zmax = value;
    updatefacelists = 1;
}

/* ------------------ setrenderdir ------------------------ */

int setrenderdir(const char *dir) {
  printf("c_api: setting renderdir to: %s\n", dir);
  // TODO: as lua gives us consts, but most smv code uses non-const, we
  // must make a non-const copy
  int l = strlen(dir);
  char *dir_path_temp = malloc(l+1);
  strncpy(dir_path_temp,dir,l+1);
	if(dir!=NULL&&strlen(dir_path_temp)>0){
      if(can_write_to_dir(dir_path_temp)==0){
        fprintf(stderr,"*** Error: Cannot write to the RENDERDIR "
                "directory: %s\n",dir_path_temp);
        return 1;
      } else {
        free(script_dir_path);
        script_dir_path=dir_path_temp;
        PRINTF("c_api: renderdir set to: %s\n", script_dir_path);
        return 0;
      }
   	} else {
      // TODO: why would we ever want to set the render directory to NULL
      script_dir_path=NULL;
      FREEMEMORY(dir_path_temp);
      return 1;
    }

}

/* ------------------ setcolorbarindex ------------------------ */
void setcolorbarindex(int chosen_index) {
	UpdateRGBColors(chosen_index);
}

/* ------------------ setcolorbarindex ------------------------ */
int getcolorbarindex() {
	return global_colorbar_index;
}

/* ------------------ setwindowsize ------------------------ */
void setwindowsize(int width, int height) {
  printf("Setting window size to %dx%d\n", width, height);
  glutReshapeWindow(width,height);
  ResizeWindow(width, height);
  Reshape_CB(width, height);

}

/* ------------------ setgridvisibility ------------------------ */

void setgridvisibility(int selection) {
	visGrid = selection;
	// selection may be one of:
	// - noGridnoProbe
	// - GridnoProbe
	// - GridProbe
	// - noGridProbe
	if(visGrid==GridProbe||visGrid==noGridProbe)visgridloc=1;
}

/* ------------------ setgridparms ------------------------ */

void setgridparms(int x_vis, int y_vis, int z_vis,
                  int x_plot, int y_plot, int z_plot) {
	visx_all = x_vis;
	visy_all = y_vis;
	visz_all = z_vis;

	iplotx_all = x_plot;
	iploty_all = y_plot;
	iplotz_all = z_plot;

	if(iplotx_all>nplotx_all-1)iplotx_all=0;
	if(iploty_all>nploty_all-1)iploty_all=0;
	if(iplotz_all>nplotz_all-1)iplotz_all=0;
}

/** Set the direction of the colorbar.
 * Is used for the .ini option FLIP. If the values is TRUE, the colorbar
 * runs in the opposite direction to that specified.
 */
void setcolorbarflip(int flip) {
	colorbarflip = flip;
	update_colorbarflip();
  UpdateRGBColors(COLORBAR_INDEX_NONE);
}

/** Get whether the direction of the colorbar is flipped.
 */
int getcolorbarflip(int flip) {
    return colorbarflip;
}

// Camera API
// These function live-modify the current view by modifying "camera_current".
void camera_set_rotation_type(int rotation_type) {
  camera_current->rotation_type = rotation_type;
}

int camera_get_rotation_type() {
  return camera_current->rotation_type;
}

// TODO: How does the rotation index work.
void camera_set_rotation_index(int rotation_index) {
  camera_current->rotation_index = rotation_index;
}

int camera_get_rotation_index() {
  return camera_current->rotation_index;
}

void camera_set_viewdir(float xcen, float ycen, float zcen) {
  printf("c_api: Setting viewDir to %f %f %f\n", xcen, ycen, zcen);
  camera_current->xcen = xcen;
  camera_current->ycen = ycen;
  camera_current->zcen = zcen;
  // camera_set_xcen(xcen);
  // camera_set_ycen(ycen);
  // camera_set_zcen(zcen);
}

// xcen
void camera_set_xcen(float xcen) {
  printf("Setting xcen to %f\n", xcen);
  camera_current->xcen = xcen;
}
float camera_get_xcen() {
  return camera_current->xcen;
}

// ycen
void camera_set_ycen(float ycen) {
  printf("Setting ycen to %f\n", ycen);
  camera_current->ycen = ycen;
}
float camera_get_ycen() {
  return camera_current->ycen;
}

// zcen
void camera_set_zcen(float zcen) {
  printf("Setting zcen to %f\n", zcen);
  camera_current->zcen = zcen;
}
float camera_get_zcen() {
  return camera_current->zcen;
}

// eyex
void camera_mod_eyex(float delta) {
  camera_current->eye[0] = camera_current->eye[0] + delta;
}
void camera_set_eyex(float eyex) {
  camera_current->eye[0] = eyex;
}

float camera_get_eyex() {
  return camera_current->eye[0];
}

// eyey
void camera_mod_eyey(float delta) {
  camera_current->eye[1] = camera_current->eye[1] + delta;
}
void camera_set_eyey(float eyey) {
  camera_current->eye[1] = eyey;
}
float camera_get_eyey() {
  return camera_current->eye[1];
}

// eyez
void camera_mod_eyez(float delta) {
  camera_current->eye[2] = camera_current->eye[2] + delta;
}
void camera_set_eyez(float eyez) {
  camera_current->eye[2] = eyez;
}
float camera_get_eyez() {
  return camera_current->eye[2];
}

// azimuth
void camera_mod_az(float delta) {
  camera_current->az_elev[0] = camera_current->az_elev[0] + delta;
}
void camera_set_az(float az) {
  camera_current->az_elev[0] = az;
}
float camera_get_az() {
  return camera_current->az_elev[0];
}

// elevation
void camera_mod_elev(float delta) {
  camera_current->az_elev[1] = camera_current->az_elev[1] + delta;
}
void camera_set_elev(float elev) {
  camera_current->az_elev[1] = elev;
}
float camera_get_elev() {
  return camera_current->az_elev[1];
}

// projection_type
void camera_toggle_projection_type() {
  camera_current->projection_type = 1 - camera_current->projection_type;
}
int camera_set_projection_type(int projection_type) {
  camera_current->projection_type = projection_type;
  // 1 is orthogonal
  // 0 is perspective
  return 0;
}
int camera_get_projection_type() {
  return camera_current->projection_type;
}


// .ini config options
   // *** COLOR/LIGHTING ***
int set_ambientlight(float r, float g, float b) {
  ambientlight[0] = r;
  ambientlight[1] = g;
  ambientlight[2] = b;
  return 0;
} // AMBIENTLIGHT

int set_backgroundcolor(float r, float g, float b) {
  backgroundbasecolor[0] = r;
  backgroundbasecolor[1] = g;
  backgroundbasecolor[2] = b;
  return 0;
} // BACKGROUNDCOLOR

int set_blockcolor(float r, float g, float b) {
  block_ambient2[0] = r;
  block_ambient2[1] = g;
  block_ambient2[2] = b;
  return 0;
} // BLOCKCOLOR

int set_blockshininess(float v) {
  block_shininess = v;
  return 0;
} // BLOCKSHININESS

int set_blockspecular(float r, float g, float b) {
  block_specular2[0] = r;
  block_specular2[1] = g;
  block_specular2[2] = b;
  return 0;
} // BLOCKSPECULAR

int set_boundcolor(float r, float g, float b) {
  boundcolor[0] = r;
  boundcolor[1] = g;
  boundcolor[2] = b;
  return 0;
} // BOUNDCOLOR

int set_colorbar_textureflag(int v)  {
  usetexturebar = v;
  return 0;
}
int get_colorbar_textureflag() {
  return usetexturebar;
}

int set_colorbar_colors(int ncolors, float colors[][3]) {
  float *rgb_ini_copy;
  float *rgb_ini_copy_p;
  CheckMemory;
  if(NewMemory((void **)&rgb_ini_copy, 4 * ncolors*sizeof(float)) == 0)return 2;
  rgb_ini_copy_p = rgb_ini_copy;
  for (int i = 0; i < ncolors; i++) {
    float *r = rgb_ini_copy_p;
    float *g = rgb_ini_copy_p + 1;
    float *b = rgb_ini_copy_p + 2;
    *r = colors[i][0];
    *g = colors[i][1];
    *b = colors[i][2];
    rgb_ini_copy_p +=3;
  }

  FREEMEMORY(rgb_ini);
  rgb_ini = rgb_ini_copy;
  nrgb_ini = ncolors;
  initrgb();
  return 0;
}

int set_color2bar_colors(int ncolors, float colors[][3]) {
  float *rgb_ini_copy;
  float *rgb_ini_copy_p;
  CheckMemory;
  if(NewMemory((void **)&rgb_ini_copy, 4 * ncolors*sizeof(float)) == 0)return 2;;
  rgb_ini_copy_p = rgb_ini_copy;
  for (int i = 0; i < ncolors; i++) {
    float *r = rgb_ini_copy_p;
    float *g = rgb_ini_copy_p + 1;
    float *b = rgb_ini_copy_p + 2;
    *r = colors[i][0];
    *g = colors[i][1];
    *b = colors[i][2];
    rgb_ini_copy_p +=3;
  }

  FREEMEMORY(rgb2_ini);
  rgb2_ini = rgb_ini_copy;
  nrgb2_ini = ncolors;
  return 0;
}

int set_diffuselight(float r, float g, float b) {
  diffuselight[0] = r;
  diffuselight[1] = g;
  diffuselight[2] = b;
  return 0;
}// DIFFUSELIGHT

int set_directioncolor(float r, float g, float b) {
  direction_color[0] = r;
  direction_color[1] = g;
  direction_color[2] = b;
  return 0;
} // DIRECTIONCOLOR

int set_flip(int v) {
  background_flip = v;
  return 0;
} // FLIP

int set_foregroundcolor(float r, float g, float b) {
  foregroundbasecolor[0] = r;
  foregroundbasecolor[1] = g;
  foregroundbasecolor[2] = b;
  return 0;
} // FOREGROUNDCOLOR

int set_heatoffcolor(float r, float g, float b) {
  heatoffcolor[0] = r;
  heatoffcolor[1] = g;
  heatoffcolor[2] = b;
  return 0;
} // HEATOFFCOLOR

int set_heatoncolor(float r, float g, float b) {
  heatoncolor[0] = r;
  heatoncolor[1] = g;
  heatoncolor[2] = b;
  return 0;
} // HEATONCOLOR

int set_isocolors(float shininess, float default_opaqueness, float specular[3],
                  int nlevels, float colors[][4]) {

  iso_shininess = shininess;
  iso_transparency = default_opaqueness;
  iso_specular[0] = specular[0];
  iso_specular[1] = specular[1];
  iso_specular[2] = specular[2];

  for (int i = 0; i < nlevels; i++) {
    iso_colors[0] = CLAMP(colors[i][0], 0.0, 1.0);
    iso_colors[1] = CLAMP(colors[i][1], 0.0, 1.0);
    iso_colors[2] = CLAMP(colors[i][2], 0.0, 1.0);
    iso_colors[3] = CLAMP(colors[i][3], 0.0, 1.0);
  }
  update_isocolors();
  update_iso_colorlevel();
  return 0;
} // ISOCOLORS

int set_colortable(int ncolors, int colors[][4], char **names) {
  int nctableinfo;
  int i;
  colortabledata *ctableinfo = NULL;
  nctableinfo = ncolors;
  nctableinfo = MAX(nctableinfo, 0);
  if(nctableinfo>0){
    NewMemory((void **)&ctableinfo, nctableinfo*sizeof(colortabledata));
    for(i = 0; i<nctableinfo; i++){
      colortabledata *rgbi;
      rgbi = ctableinfo + i;
      // TODO: This sets the default alpha value to 255, as per the
      // original readsmv.c function, but is defunct in this context
      // as this value is required by the function prototype.
      // color[i][3] = 255;
      strcpy(rgbi->label, names[i]);
      rgbi->color[0] = CLAMP(colors[i][0], 0, 255);
      rgbi->color[1] = CLAMP(colors[i][1], 0, 255);
      rgbi->color[2] = CLAMP(colors[i][2], 0, 255);
      rgbi->color[3] = CLAMP(colors[i][3], 0, 255);
    }
    UpdateColorTable(ctableinfo, nctableinfo);
    FREEMEMORY(ctableinfo);
  }
  return 0;
} // COLORTABLE

int set_light0(int setting) {
  light_enabled0 = setting;
  UpdateLIGHTS = 1;
  return 0;
} // LIGHT0

int set_light1(int setting) {
  light_enabled1 = setting;
  UpdateLIGHTS = 1;
  return 0;
} // LIGHT1

int set_lightpos0(float a, float b, float c, float d) {
  light_position0[0] = a;
  light_position0[1] = a;
  light_position0[2] = a;
  light_position0[3] = a;
  return 0;
} // LIGHTPOS0

int set_lightpos1(float a, float b, float c, float d) {
  light_position1[0] = a;
  light_position1[1] = a;
  light_position1[2] = a;
  light_position1[3] = a;
  return 0;
} // LIGHTPOS1

int set_lightmodellocalviewer(int setting) {
  lightmodel_localviewer = setting == 0 ? GL_FALSE : GL_TRUE;
  UpdateLIGHTS = 1;
  return 0;
} // LIGHTMODELLOCALVIEWER

int set_lightmodelseparatespecularcolor(int setting) {
  lightmodel_separatespecularcolor = setting;
  UpdateLIGHTS = 1;
  return 0;
} // LIGHTMODELSEPARATESPECULARCOLOR

int set_sensorcolor(float r, float g, float b) {
  sensorcolor[0] = r;
  sensorcolor[1] = g;
  sensorcolor[2] = b;
  return 0;
} // SENSORCOLOR

int set_sensornormcolor(float r, float g, float b) {
  sensornormcolor[0] = r;
  sensornormcolor[1] = g;
  sensornormcolor[2] = b;
  return 0;
} // SENSORNORMCOLOR

int set_bw(int geo_setting, int data_setting) {
  setbw = geo_setting;
  setbwdata = data_setting;
} // SETBW

int set_sprinkleroffcolor(float r, float g, float b) {
  sprinkoffcolor[0] = r;
  sprinkoffcolor[1] = g;
  sprinkoffcolor[2] = b;
  return 0;
} // SPRINKOFFCOLOR

int set_sprinkleroncolor(float r, float g, float b) {
  sprinkoncolor[0] = r;
  sprinkoncolor[1] = g;
  sprinkoncolor[2] = b;
  return 0;
} // SPRINKONCOLOR

int set_staticpartcolor(float r, float g, float b) {
  static_color[0] = r;
  static_color[1] = g;
  static_color[2] = b;
  return 0;
} // STATICPARTCOLOR

int set_timebarcolor(float r, float g, float b) {
  timebarcolor[0] = r;
  timebarcolor[1] = g;
  timebarcolor[2] = b;
  return 0;
} // TIMEBARCOLOR

int set_ventcolor(float r, float g, float b) {
  ventcolor[0] = r;
  ventcolor[1] = g;
  ventcolor[2] = b;
  return 0;
} // VENTCOLOR


// --    *** SIZES/OFFSETS ***
int set_gridlinewidth(float v) {
  gridlinewidth = v;
  return 0;
} // GRIDLINEWIDTH

int set_isolinewidth(float v) {
  isolinewidth = v;
  return 0;
} // ISOLINEWIDTH

int set_isopointsize(float v) {
  isopointsize = v;
  return 0;
} // ISOPOINTSIZE

int set_linewidth(float v) {
  linewidth = v;
  return 0;
} // LINEWIDTH

int set_partpointsize(float v) {
  partpointsize = v;
  return 0;
} // PARTPOINTSIZE

int set_plot3dlinewidth(float v){
  plot3dlinewidth = v;
  return 0;
} // PLOT3DLINEWIDTH

int set_plot3dpointsize(float v) {
  plot3dpointsize = v;
  return 0;
} // PLOT3DPOINTSIZE

int set_sensorabssize(float v) {
  sensorabssize = v;
  return 0;
} // SENSORABSSIZE

int set_sensorrelsize(float v) {
  sensorrelsize = v;
  return 0;
} // SENSORRELSIZE

int set_sliceoffset(float v) {
  sliceoffset_factor = v;
  return 0;
} // SLICEOFFSET

int set_smoothlines(int v) {
  antialiasflag = v;
  return 0;
} // SMOOTHLINES

int set_spheresegs(int v) {
  device_sphere_segments = v;
  return 0;
} // SPHERESEGS

int set_sprinklerabssize(float v) {
  sprinklerabssize = v;
  return 0;
} // SPRINKLERABSSIZE

int set_streaklinewidth(float v) {
  streaklinewidth = v;
  return 0;
} // STREAKLINEWIDTH

int set_ticklinewidth(float v) {
  ticklinewidth = v;
  return 0;
} // TICKLINEWIDTH

int set_usenewdrawface(int v) {
  use_new_drawface = v;
  return 0;
} // USENEWDRAWFACE

int set_veccontours(int v) {
  show_slices_and_vectors = v;
  return 0;
} // VECCONTOURS

int set_veclength(int a, float b, float c) {
  int dummy1 = a; // TODO: what is the point of this value
  vecfactor = b;
  int dummy2 = c; // TODO: what is the point of this value
  return 0;
} // VECLENGTH

int set_vectorlinewidth(float a, float b) {
  vectorlinewidth = a;
  slice_line_contour_width = b;
  return 0;
} // VECTORLINEWIDTH

int set_vectorpointsize(float v) {
  vectorpointsize = v;
  return 0;
} // VECTORPOINTSIZE

int set_ventlinewidth(float v) {
  ventlinewidth = v;
  return 0;
} // VENTLINEWIDTH

int set_ventoffset(float v) {
  ventoffset_factor = v;
  return 0;
} // VENTOFFSET

int set_windowoffset(int v) {
  titlesafe_offsetBASE = v;
  return 0;
} // WINDOWOFFSET

int set_windowwidth(int v) {
  screenWidth = v;
  return 0;
} // WINDOWWIDTH

int set_windowheight(int v) {
  screenHeight = v;
  return 0;
} // WINDOWHEIGHT


// --  *** DATA LOADING ***
int set_boundzipstep(int v) {
  boundzipstep = v;
  return 0;
} // BOUNDZIPSTEP

int set_fed(int v) {
  regenerate_fed = v;
  return 0;
} // FED

int set_fedcolorbar(const char *name) {
  if(strlen(name)>0) {
    strcpy(default_fed_colorbar, name);
    return 0;
  } else {
    return 1;
  }
} // FEDCOLORBAR

int set_isozipstep(int v) {
  isozipstep = v;
  return 0;
} // ISOZIPSTEP

int set_nopart(int v) {
  nopart = v;
  return 0;
} // NOPART

// int set_partpointstep(int v) {
//   partpointstep = v;
//   return 0;
// } // PARTPOINTSTEP

int set_showfedarea(int v) {
  show_fed_area = v;
  return 0;
} // SHOWFEDAREA

int set_sliceaverage(int flag, float interval, int vis) {
  slice_average_flag = flag;
  slice_average_interval = interval;
  vis_slice_average = vis;
  return 0;
} // SLICEAVERAGE

int set_slicedataout(int v) {
  output_slicedata = v;
  return 0;
} // SLICEDATAOUT

int set_slicezipstep(int v) {
  slicezipstep = v;
  return 0;
} // SLICEZIPSTEP

int set_smoke3dzipstep(int v) {
  smoke3dzipstep = v;
  return 0;
} // SMOKE3DZIPSTEP

int set_userrotate(int index, int show_center, float x, float y, float z) {
  glui_rotation_index = index;
  slice_average_interval = show_center;
  xcenCUSTOM = x;
  ycenCUSTOM = y;
  zcenCUSTOM = z;
  return 0;
} // USER_ROTATE


// --  *** VIEW PARAMETERS ***
int set_aperture(int v) {
  apertureindex = v;
  return 0;
} // APERTURE

int set_axissmooth(int v) {
  axislabels_smooth = v;
  return 0;
} // AXISSMOOTH

// provided above
int set_blocklocation(int v) {
  blocklocation = v;
  return 0;
} // BLOCKLOCATION

int set_boundarytwoside(int v) {
  showpatch_both = v;
  return 0;
} // BOUNDARYTWOSIDE

int set_clip(float v_near, float v_far) {
  nearclip = v_near;
  farclip = v_far;
  return 0;
} // CLIP

int set_contourtype(int v) {
  contour_type = v;
  return 0;
} // CONTOURTYPE

int set_cullfaces(int v) {
  cullfaces = v;
  return 0;
} // CULLFACES

int set_texturelighting(int v) {
  enable_texture_lighting = v;
  return 0;
} // ENABLETEXTURELIGHTING

int set_eyeview(int v) {
  rotation_type = v;
  return 0;
} // EYEVIEW

int set_eyex(float v) {
  eyexfactor = v;
  return 0;
} // EYEX

int set_eyey(float v) {
  eyeyfactor = v;
  return 0;
} // EYEY

int set_eyez(float v) {
  eyezfactor = v;
  return 0;
} // EYEZ

int set_fontsize(int v) {
  fontindex = v;
  return 0;
} // FONTSIZE

int set_frameratevalue(int v) {
  frameratevalue = v;
  return 0;
} // FRAMERATEVALUE

int set_geomdiags(int structured, int unstructured, int diagnostics) {
  structured_isopen = structured;
  unstructured = unstructured_isopen;
  show_geometry_diagnostics = diagnostics;
  return 0;
} // GEOMDIAGS

// int set_geomshow(int )
// GEOMSHOW
int set_showfaces_interior(int v) {
  show_faces_interior = v;
  return 0;
}
int set_showfaces_exterior(int v) {
  show_faces_exterior = v;
  return 0;
}
int set_showfaces_solid(int v) {
  frameratevalue = v;
  return 0;
}
int set_showfaces_outline(int v) {
  show_faces_outline = v;
  return 0;
}
int set_smoothgeomnormal(int v) {
  smooth_geom_normal = v;
  return 0;
}
int set_showvolumes_interior(int v) {
  show_volumes_interior = v;
  return 0;
}
int set_showvolumes_exterior(int v) {
  show_volumes_exterior = v;
  return 0;
}
int set_showvolumes_solid(int v) {
  show_volumes_solid = v;
  return 0;
}
int set_showvolumes_outline(int v) {
  show_volumes_outline = v;
  return 0;
}
int set_geomvertexag(int v) {
  geom_vert_exag = v;
  return 0;
}
int set_geommaxangle(int v) {
  geom_max_angle = v;
  return 0;
}

int set_gversion(int v) {
  gversion = v;
  return 0;
} // GVERSION

int set_isotran2(int v) {
  transparent_state = v;
  return 0;
} // ISOTRAN2

int set_meshvis(int n, int vals[]) {
  meshdata *meshi;
  for (int i = 0; i < n; i++) {
    if(i>nmeshes - 1)break;
    meshi = meshinfo + i;
    meshi->blockvis = vals[i];
    ONEORZERO(meshi->blockvis);
  }
  return 0;
} // MESHVIS

int set_meshoffset(int meshnum, int value) {
  if(meshnum >= 0 && meshnum<nmeshes){
    meshdata *meshi;

    meshi = meshinfo + meshnum;
    meshi->mesh_offset_ptr = meshi->mesh_offset;
    return 0;
  }
  return 1;
} // MESHOFFSET

int set_northangle(int vis, float x, float y, float z) {
  vis_northangle = vis;
  northangle_position[0] = x;
  northangle_position[1] = y;
  northangle_position[2] = z;
  return 0;
} // NORTHANGLE

int set_offsetslice(int v) {
  offset_slice = v;
  return 0;
} // OFFSETSLICE

int set_outlinemode(int a, int b) {
  highlight_flag = a;
  outline_color_flag = b;
  return 0;
} // OUTLINEMODE

int set_p3dsurfacetype(int v) {
  p3dsurfacetype = v;
  return 0;
} // P3DSURFACETYPE

int set_p3dsurfacesmooth(int v) {
  p3dsurfacesmooth = v;
  return 0;
} // P3DSURFACESMOOTH

int set_projection(int v) {
  projection_type = v;
  return 0;
} // PROJECTION

int set_scaledfont(int height2d, float height2dwidth, int thickness2d,
                   int height3d, float height3dwidth, int thickness3d) {
  scaled_font2d_height = scaled_font2d_height;
  scaled_font2d_height2width = scaled_font2d_height2width;
  thickness2d = scaled_font2d_thickness;
  scaled_font3d_height = scaled_font3d_height;
  scaled_font3d_height2width = scaled_font3d_height2width;
  thickness3d = scaled_font3d_thickness;
  return 0;
} // SCALEDFONT

int set_showalltextures(int v) {
  showall_textures = v;
  return 0;
} // SHOWALLTEXTURES

int set_showaxislabels(int v) {
  visaxislabels = v;
  return 0;
} // SHOWAXISLABELS TODO: duplicate

int set_showblocklabel(int v) {
  visMeshlabel = v;
  return 0;
} // SHOWBLOCKLABEL

int set_showblocks(int v) {
  visBlocks = v;
  return 0;
} // SHOWBLOCKS

int set_showcadandgrid(int v) {
  show_cad_and_grid = v;
  return 0;
} // SHOWCADANDGRID

int set_showcadopaque(int v) {
  viscadopaque = v;
  return 0;
} // SHOWCADOPAQUE

int set_showceiling(int v) {
  visCeiling = v;
  return 0;
} // SHOWCEILING

int set_showcolorbars(int v) {
  visColorbar = v;
  return 0;
} // SHOWCOLORBARS

int set_showcvents(int a, int b) {
  visCircularVents = a;
  circle_outline = b;
  return 0;
} // SHOWCVENTS

int set_showdummyvents(int v) {
  visDummyVents = v;
  return 0;
} // SHOWDUMMYVENTS

int set_showevacslices(int a, int b, int c) {
  show_evac_slices = a;
  constant_evac_coloring = b;
  show_evac_colorbar = c;
  return 0;
} // SHOWEVACSLICES

int set_showfloor(int v) {
  visFloor = v;
  return 0;
} // SHOWFLOOR

int set_showframe(int v) {
  visFrame = v;
  return 0;
} // SHOWFRAME

int set_showframelabel(int v) {
  visFramelabel = v;
  return 0;
} // SHOWFRAMELABEL

int set_showframerate(int v) {
  visFramerate = v;
  return 0;
} // SHOWFRAMERATE

int set_showgrid(int v) {
  visGrid = v;
  return 0;
} // SHOWGRID

int set_showgridloc(int v) {
  visgridloc = v;
  return 0;
} // SHOWGRIDLOC

int set_showhmstimelabel(int v) {
  vishmsTimelabel = v;
  return 0;
} // SHOWHMSTIMELABEL

int set_showhrrcutoff(int v) {
  visHRRlabel = v;
  return 0;
} // SHOWHRRCUTOFF

int set_showiso(int v) {
  visAIso = v;
  return 0;
} // SHOWISO

int set_showisonormals(int v) {
  show_iso_normal = v;
  return 0;
} // SHOWISONORMALS

int set_showlabels(int v) {
  visLabels = v;
  return 0;
} // SHOWLABELS

#ifdef pp_memstatus
int set_showmemload(int v) {
  visAvailmemory = v;
  return 0;
} // SHOWMEMLOAD
#endif

// int set_shownormalwhensmooth(int v); // SHOWNORMALWHENSMOOTH
int set_showopenvents(int a, int b) {
  visOpenVents = a;
  visOpenVentsAsOutline = b;
  return 0;
} // SHOWOPENVENTS

int set_showothervents(int v) {
  visOtherVents = v;
  return 0;
} // SHOWOTHERVENTS

int set_showsensors(int a, int b) {
  visSensor = a;
  visSensorNorm = b;
  return 0;
} // SHOWSENSORS

int set_showsliceinobst(int v) {
  show_slice_in_obst = v;
  return 0;
} // SHOWSLICEINOBST

int set_showsmokepart(int v) {
  visSmokePart = v;
  return 0;
} // SHOWSMOKEPART

int set_showsprinkpart(int v) {
  visSprinkPart = v;
  return 0;
} // SHOWSPRINKPART

int set_showstreak(int show, int step, int showhead, int index) {
  streak5show = show;
  streak5step = step;
  showstreakhead = showhead;
  streak_index = index;
  return 0;
} // SHOWSTREAK

int set_showterrain(int v){
  visTerrainType = v;
  return 0;
} // SHOWTERRAIN

int set_showtetras(int a, int b){
  show_volumes_solid = a;
  show_volumes_outline = b;
  return 0;
} // SHOWTETRAS

int set_showthreshold(int a, int b, float c){
  vis_threshold = a;
  vis_onlythreshold = b;
  temp_threshold = c;
  return 0;
} // SHOWTHRESHOLD

int set_showticks(int v){
  visFDSticks = v;
  return 0;
} // SHOWTICKS

int set_showtimebar(int v){
  visTimebar = v;
  return 0;
} // SHOWTIMEBAR

int set_showtimelabel(int v){
  visTimelabel = v;
  return 0;
} // SHOWTIMELABEL

int set_showtitle(int v){
  visTitle = v;
  return 0;
} // SHOWTITLE

int set_showtracersalways(int v){
  show_tracers_always = v;
  return 0;
} // SHOWTRACERSALWAYS

int set_showtriangles(int a, int b, int c, int d, int e, int f){
  show_iso_solid = a;
  show_iso_outline = b;
  show_iso_verts = c;
  show_iso_normal = d;
  smooth_iso_normal = e;
  return 0;
} // SHOWTRIANGLES

int set_showtransparent(int v){
  visTransparentBlockage = v;
  return 0;
} // SHOWTRANSPARENT

int set_showtransparentvents(int v){
  show_transparent_vents = v;
  return 0;
} // SHOWTRANSPARENTVENTS

int set_showtrianglecount(int v){
  show_triangle_count = v;
  return 0;
} // SHOWTRIANGLECOUNT

int set_showventflow(int a, int b, int c, int d, int e){
  visVentHFlow = a;
  visventslab = b;
  visventprofile = c;
  visVentVFlow = d;
  visVentMFlow = e;
  return 0;
} // SHOWVENTFLOW

int set_showvents(int v) {
  visVents = v;
  return 0;
} // SHOWVENTS

int set_showwalls(int v) {
  visWalls = v;
  return 0;
} // SHOWWALLS

int set_skipembedslice(int v) {
  skip_slice_in_embedded_mesh = v;
  return 0;
} // SKIPEMBEDSLICE

#ifdef pp_SLICEUP
int set_slicedup(int a, int b){
  slicedup_option = a;
  vectorslicedup_option = b;
  return 0;
} // SLICEDUP
#endif

int set_smokesensors(int a, int b){
  show_smokesensors = a;
  test_smokesensors = b;
  return 0;
} // SMOKESENSORS

// int set_smoothblocksolid(int v); // SMOOTHBLOCKSOLID
#ifdef pp_LANG
int set_startuplang(const char *lang) {
  char *bufptr;

  strncpy(startup_lang_code, lang, 2);
  startup_lang_code[2] = '\0';
  if(strcmp(startup_lang_code, "en") != 0){
    show_lang_menu = 1;
  }
  if(tr_name == NULL){
    int langlen;

    langlen = strlen(bufptr);
    NewMemory((void **)&tr_name, langlen + 48 + 1);
    strcpy(tr_name, bufptr);
  }
  return 0;
} // STARTUPLANG
#endif

int set_stereo(int v) {
  showstereo = v;
  return 0;
} // STEREO

int set_surfinc(int v) {
  surfincrement = v;
  return 0;
} // SURFINC

int set_terrainparams(int r_min, int g_min, int b_min,
                      int r_max, int g_max, int b_max, int v) {
  terrain_rgba_zmin[0] = r_min;
  terrain_rgba_zmin[1] = g_min;
  terrain_rgba_zmin[2] = b_min;
  terrain_rgba_zmin[0] = r_max;
  terrain_rgba_zmin[1] = g_max;
  terrain_rgba_zmin[2] = b_max;
  vertical_factor = v;
  return 0;
} // TERRAINPARMS

int set_titlesafe(int v) {
  titlesafe_offset = v;
  return 0;
} // TITLESAFE

int set_trainermode(int v) {
  trainer_mode = v;
  return 0;
} // TRAINERMODE

int set_trainerview(int v) {
  trainerview = v;
  return 0;
} // TRAINERVIEW

int set_transparent(int use_flag, float level) {
  use_transparency_data = use_flag;
  transparent_level = level;
} // TRANSPARENT

int set_treeparms(int minsize, int visx, int visy, int visz) {
  mintreesize = minsize;
  vis_xtree = visx;
  vis_ytree = visy;
  vis_ztree = visz;
  return 0;
} // TREEPARMS

int set_twosidedvents(int internal, int external) {
  show_bothsides_int = internal;
  show_bothsides_ext = external;
  return 0;
} // TWOSIDEDVENTS

int set_vectorskip(int v) {
  vectorskip = v;
  return 0;
} // VECTORSKIP

int set_volsmoke(int a, int b, int c, int d, int e,
                 float f, float g, float h, float i,
                 float j, float k, float l) {
  glui_compress_volsmoke = a;
  use_multi_threading = b;
  load_at_rendertimes = c;
  volbw = d;
  show_volsmoke_moving = e;
  temperature_min = f;
  temperature_cutoff = g;
  temperature_max = h;
  fire_opacity_factor = i;
  mass_extinct = j;
  gpu_vol_factor = k;
  nongpu_vol_factor = l;
  return 0;
} // VOLSMOKE

int set_zoom(int a, float b) {
  zoomindex = a;
  zoom = b;
  return 0;
} // ZOOM

// *** MISC ***
int set_cellcentertext(int v) {
  cell_center_text = v;
  return 0;
} // CELLCENTERTEXT

int set_inputfile(const char *filename) {
  size_t len;
  len = strlen(filename);

  FREEMEMORY(INI_fds_filein);
  if(NewMemory((void **)&INI_fds_filein, (unsigned int)(len + 1)) == 0)return 2;
  STRCPY(INI_fds_filein, filename);
  return 0;
} // INPUT_FILE

int set_labelstartupview(const char *startupview) {
  strcpy(label_startup_view, startupview);
  update_startup_view = 1;
  return 0;
} // LABELSTARTUPVIEW

int set_pixelskip(int v) {
  pixel_skip = v;
  return 0;
} // PIXELSKIP

int set_renderclip(int use_flag, int left, int right, int bottom, int top) {
  clip_rendered_scene = use_flag;
  render_clip_left = left;
  render_clip_right = right;
  render_clip_bottom = bottom;
  render_clip_top = top;
  return 0;
} // RENDERCLIP

int set_renderfilelabel(int v) {
  renderfilelabel = v;
  return 0;
} // RENDERFILELABEL

int set_renderfiletype(int render, int movie) {
  renderfiletype = render;
  moviefiletype = movie;
  return 0;
} // RENDERFILETYPE

int set_skybox() {
  // skyboxdata *skyi;

  // free_skybox();
  // nskyboxinfo = 1;
  // NewMemory((void **)&skyboxinfo, nskyboxinfo*sizeof(skyboxdata));
  // skyi = skyboxinfo;

  // for(i = 0; i<6; i++){
  //   fgets(buffer, 255, stream);
  //   loadskytexture(buffer, skyi->face + i);
  // }
  ASSERT(FFALSE);
  return 0;
} // SKYBOX TODO

int set_renderoption(int opt, int rows) {
  render_option = opt;
  nrender_rows = rows;
  return 0;
} // RENDEROPTION

int set_unitclasses(int n, int indices[]) {
  int i;
  for(i = 0; i<n; i++){
    if(i>nunitclasses - 1)continue;
    unitclasses[i].unit_index = indices[i];
  }
  return 0;
} // UNITCLASSES

int set_zaxisangles(float a, float b, float c) {
  zaxis_angles[0] = a;
  zaxis_angles[1] = b;
  zaxis_angles[2] = c;
  return 0;
}

// *** 3D SMOKE INFO ***
int set_adjustalpha(int v) {
  adjustalphaflag = v;
  return 0;
} // ADJUSTALPHA

int set_colorbartype(int type, const char *label) {
  update_colorbartype = 1;
  colorbartype = type;
  strcpy(colorbarname, label);
  return 0;
} // COLORBARTYPE

int set_extremecolors(int rmin, int gmin, int bmin,
                      int rmax, int gmax, int bmax) {
  rgb_below_min[0] = CLAMP(rmin, 0, 255);
  rgb_below_min[1] = CLAMP(gmin, 0, 255);
  rgb_below_min[2] = CLAMP(bmin, 0, 255);

  rgb_above_max[0] = CLAMP(rmax, 0, 255);
  rgb_above_max[1] = CLAMP(gmax, 0, 255);
  rgb_above_max[2] = CLAMP(bmax, 0, 255);
  return 0;
} // EXTREMECOLORS

int set_firecolor(int r, int g, int b) {
  fire_red = r;
  fire_green = g;
  fire_blue = b;
  return 0;
} // FIRECOLOR

int set_firecolormap(int type, int index) {
  firecolormap_type = type;
  fire_colorbar_index = index;
  return 0;
} // FIRECOLORMAP

int set_firedepth(float v) {
  fire_halfdepth = v;
  return 0;
} // FIREDEPTH

// int set_gcolorbar(int ncolorbarini, ) {
//   colorbardata *cbi;
//   int r1, g1, b1;
//   int n;

//   initdefaultcolorbars();

//   ncolorbars = ndefaultcolorbars + ncolorbarini;
//   if(ncolorbarini>0)ResizeMemory((void **)&colorbarinfo, ncolorbars*sizeof(colorbardata));

//   for(n = ndefaultcolorbars; n<ncolorbars; n++){
//     char *cb_buffptr;

//     cbi = colorbarinfo + n;
//     fgets(buffer, 255, stream);
//     trim_back(buffer);
//     cb_buffptr = trim_front(buffer);
//     strcpy(cbi->label, cb_buffptr);

//     fgets(buffer, 255, stream);
//     sscanf(buffer, "%i %i", &cbi->nnodes, &cbi->nodehilight);
//     if(cbi->nnodes<0)cbi->nnodes = 0;
//     if(cbi->nodehilight<0 || cbi->nodehilight >= cbi->nnodes){
//       cbi->nodehilight = 0;
//     }

//     cbi->label_ptr = cbi->label;
//     for(i = 0; i<cbi->nnodes; i++){
//       int icbar;
//       int nn;

//       fgets(buffer, 255, stream);
//       r1 = -1; g1 = -1; b1 = -1;
//       sscanf(buffer, "%i %i %i %i", &icbar, &r1, &g1, &b1);
//       cbi->index_node[i] = icbar;
//       nn = 3 * i;
//       cbi->rgb_node[nn] = r1;
//       cbi->rgb_node[nn + 1] = g1;
//       cbi->rgb_node[nn + 2] = b1;
//     }
//     remapcolorbar(cbi);
//   }
//   return 0;
// } // GCOLORBAR

int set_showextremedata(int show_extremedata, int below, int above) {
  // int below = -1, above = -1, show_extremedata;
  if(below == -1 && above == -1){
    if(below == -1)below = 0;
    if(below != 0)below = 1;
    if(above == -1)above = 0;
    if(above != 0)above = 1;
  }
  else{
    if(show_extremedata != 1)show_extremedata = 0;
    if(show_extremedata == 1){
      below = 1;
      above = 1;
    }
    else{
      below = 0;
      above = 0;
    }
  }
  show_extreme_mindata = below;
  show_extreme_maxdata = above;
  return 0;
} // SHOWEXTREMEDATA

int set_smokecolor(int r, int g, int b) {
  smoke_red = r;
  smoke_green = g;
  smoke_blue = b;
  return 0;
} // SMOKECOLOR

int set_smokecull(int v) {
#ifdef pp_CULL
  if(gpuactive == 1){
    cullsmoke = v;
    if(cullsmoke != 0)cullsmoke = 1;
  }
  else{
    cullsmoke = 0;
  }
#else
  smokecullflag = v;
#endif
  return 0;
} // SMOKECULL

int set_smokeskip(int v) {
  smokeskipm1 = v;
  return 0;
} // SMOKESKIP

int set_smokealbedo(float v) {
  smoke_albedo = v;
  return 0;
} // SMOKEALBEDO

#ifdef pp_GPU
int set_smokerthick(float v) {
  smoke3d_rthick = v;
  smoke3d_rthick = CLAMP(smoke3d_rthick, 1.0, 255.0);
  smoke3d_thick = log_base2(smoke3d_rthick);
  return 0;
} // SMOKERTHICK
#endif

int set_smokethick(float v) {
  smoke3d_thick = v;
  return 0;
} // SMOKETHICK

#ifdef pp_GPU
int set_usegpu(int v) {
  usegpu = v;
  return 0;
}
#endif

// *** ZONE FIRE PARAMETRES ***
int set_showhazardcolors(int v) {
  zonecolortype = v;
  return 0;
} // SHOWHAZARDCOLORS

int set_showhzone(int v) {
  visHZone = v;
  return 0;
} // SHOWHZONE

int set_showszone(int v) {
  visSZone = v;
  return 0;
} // SHOWSZONE

int set_showvzone(int v) {
  visVZone = v;
  return 0;
} // SHOWVZONE

int set_showzonefire(int v) {
  viszonefire = v;
  return 0;
} // SHOWZONEFIRE

// *** TOUR INFO ***
int set_showpathnodes(int v) {
  show_path_knots = v;
  return 0;
} // SHOWPATHNODES

int set_showtourroute(int v) {
  edittour = v;
  return 0;
} // SHOWTOURROUTE

// TOURCOLORS
int set_tourcolors_selectedpathline(float r, float g, float b) {
  tourcol_selectedpathline[0] = r;
  tourcol_selectedpathline[1] = g;
  tourcol_selectedpathline[2] = b;
  return 0;
}
int set_tourcolors_selectedpathlineknots(float r, float g, float b) {
  tourcol_selectedpathlineknots[0] = r;
  tourcol_selectedpathlineknots[1] = g;
  tourcol_selectedpathlineknots[2] = b;
  return 0;
}
int set_tourcolors_selectedknot(float r, float g, float b) {
  tourcol_selectedknot[0] = r;
  tourcol_selectedknot[1] = g;
  tourcol_selectedknot[2] = b;
  return 0;
}
int set_tourcolors_pathline(float r, float g, float b) {
  tourcol_pathline[0] = r;
  tourcol_pathline[1] = g;
  tourcol_pathline[2] = b;
  return 0;
}
int set_tourcolors_pathknots(float r, float g, float b) {
  tourcol_pathknots[0] = r;
  tourcol_pathknots[1] = g;
  tourcol_pathknots[2] = b;
  return 0;
}
int set_tourcolors_text(float r, float g, float b) {
  tourcol_text[0] = r;
  tourcol_text[1] = g;
  tourcol_text[2] = b;
  return 0;
}
int set_tourcolors_avatar(float r, float g, float b) {
  tourcol_avatar[0] = r;
  tourcol_avatar[1] = g;
  tourcol_avatar[2] = b;
  return 0;
}

int set_tourconstantvel(int v) {
  tour_constant_vel = v;
  return 0;
} // TOURCONSTANTVEL

int set_viewalltours(int v) {
  viewalltours = v;
  return 0;
} // VIEWALLTOURS

int set_viewtimes(float start, float stop, int ntimes) {
  view_tstart = start;
  view_tstop = stop;
  view_ntimes = ntimes;
  return 0;
} // VIEWTIMES

int set_viewtourfrompath(int v) {
  viewtourfrompath = v;
  return 0;
} // VIEWTOURFROMPATH

int set_avatarevac(int v) {
  iavatar_evac = v;
  return 0;
} // AVATAREVAC

int set_geometrytest(int a, int b, float c, float d, int vals[],
                     float b1Vals[], float b2Vals[], float b3Vals[]) {
  int *v;
  int ii;
  geomtest_option = a;
  show_tetratest_labels = b;
  tetra_line_thickness = c;
  tetra_point_size = d;
  v = tetrabox_vis;
  ONEORZERO(show_tetratest_labels);
  for(ii = 0; ii<10; ii++){
    v[ii] = vals[ii];
    ONEORZERO(v[ii]);
  }
  for(ii = 0; ii<6; ii++){
    box_bounds2[ii] = b1Vals[ii];
  }
  for(ii = 0; ii<12; ii++){
     tetra_vertices[ii] = b2Vals[ii];
  }
  for(ii = 0; ii<3; ii++){
    box_translate[ii] = b3Vals[ii];
  }
  return 0;
} //  GEOMETRYTEST

int set_devicevectordimensions(float baselength, float basediameter,
                               float headlength, float headdiameter) {
  vector_baselength = baselength;
  vector_basediameter = basediameter;
  vector_headlength = headlength;
  vector_headdiameter = headdiameter;
  return 0;
} // DEVICEVECTORDIMENSIONS

int set_devicebounds(float min, float max) {
  device_valmin = min;
  device_valmax = max;
  return 0;
} // DEVICEBOUNDS

int set_deviceorientation(int a, float b) {
  show_device_orientation = a;
  orientation_scale = b;
  show_device_orientation = CLAMP(show_device_orientation, 0, 1);
  orientation_scale = CLAMP(orientation_scale, 0.1, 10.0);
} // DEVICEORIENTATION

int set_gridparms(int vx, int vy, int vz, int px, int py, int pz) {
  visx_all = vx;
  visy_all = vy;
  visz_all = vz;

  iplotx_all = px;
  iploty_all = py;
  iplotz_all = pz;

  if(iplotx_all>nplotx_all - 1)iplotx_all = 0;
  if(iploty_all>nploty_all - 1)iploty_all = 0;
  if(iplotz_all>nplotz_all - 1)iplotz_all = 0;

  return 0;
} // GRIDPARMS

int set_gsliceparms(int vis_data, int vis_triangles, int vis_triangulation,
                    int vis_normal, float xyz[], float azelev[]) {
  vis_gslice_data = vis_data;
  show_gslice_triangles = vis_triangles;
  show_gslice_triangulation = vis_triangulation;
  show_gslice_normal = vis_normal;
  ONEORZERO(vis_gslice_data);
  ONEORZERO(show_gslice_triangles);
  ONEORZERO(show_gslice_triangulation);
  ONEORZERO(show_gslice_normal);

  gslice_xyz[0] = xyz[0];
  gslice_xyz[1] = xyz[1];
  gslice_xyz[2] = xyz[2];

  gslice_normal_azelev[0] = azelev[0];
  gslice_normal_azelev[1] = azelev[1];

  update_gslice = 1;

  return 0;
} // GSLICEPARMS

int set_loadfilesatstartup(int v) {
  loadfiles_at_startup = v;
  return 0;
} // LOADFILESATSTARTUP
int set_mscale(float a, float b, float c) {
  mscale[0] = a;
  mscale[1] = b;
  mscale[2] = c;
  return 0;
} // MSCALE

int set_sliceauto(int n, int vals[]) {
  int i;
  int n3dsmokes = 0;
  int seq_id;
  n3dsmokes = n; // TODO: is n3dsmokes the right variable.
  // TODO: this discards  the values. Verify.
  for(i = 0; i<n3dsmokes; i++){
    seq_id = vals[i];
    get_startup_slice(seq_id);
  }
  update_load_Files = 1;
  return 0;
} // SLICEAUTO

int set_msliceauto(int n, int vals[]) {
  int i;
  int n3dsmokes = 0;
  int seq_id;
  n3dsmokes = n; // TODO: is n3dsmokes the right variable
  for(i = 0; i<n3dsmokes; i++){
    seq_id = vals[i];

    if(seq_id >= 0 && seq_id<nmultisliceinfo){
      multislicedata *mslicei;

      mslicei = multisliceinfo + seq_id;
      mslicei->autoload = 1;
    }
  }
  update_load_Files = 1;
  return 0;
} // MSLICEAUTO

int set_compressauto(int v) {
  compress_autoloaded = v;
  return 0;
} // COMPRESSAUTO

// int set_part5propdisp(int vals[]) {
//   char *token;

//   for(i = 0; i<npart5prop; i++){
//     partpropdata *propi;
//     int j;

//     propi = part5propinfo + i;
//     fgets(buffer, 255, stream);

//     trim_back(buffer);
//     token = strtok(buffer, " ");
//     j = 0;
//     while(token != NULL&&j<npartclassinfo){
//       int visval;

//       sscanf(token, "%i", &visval);
//       propi->class_vis[j] = visval;
//       token = strtok(NULL, " ");
//       j++;
//     }
//   }
//   CheckMemory;
//   continue;
// } // PART5PROPDISP

// int set_part5color(int n, int vals[]) {
//   int i;
//   for(i = 0; i<npart5prop; i++){
//     partpropdata *propi;

//     propi = part5propinfo + i;
//     propi->display = 0;
//   }
//   part5colorindex = 0;
//   i = n;
//   if(i >= 0 && i<npart5prop){
//     partpropdata *propi;

//     part5colorindex = i;
//     propi = part5propinfo + i;
//     propi->display = 1;
//   }
//   continue;
//   return 0;
// } // PART5COLOR

int set_propindex(int nvals, int vals[][2]) {
  int i;
  for(i = 0; i<nvals; i++){
    propdata *propi;
    int ind, val;
    ind = vals[0][2];
    val = vals[1][2];
    if(ind<0 || ind>npropinfo - 1)return 0;
    propi = propinfo + ind;
    if(val<0 || val>propi->nsmokeview_ids - 1)return 0;
    propi->smokeview_id = propi->smokeview_ids[val];
    propi->smv_object = propi->smv_objects[val];
  }
  for(i = 0; i<npartclassinfo; i++){
    partclassdata *partclassi;

    partclassi = partclassinfo + i;
    update_partclass_depend(partclassi);

  }
  return 0;
} // PROPINDEX

int set_shooter(float xyz[], float dxyz[], float uvw[],
                float velmag, float veldir, float pointsize,
                int fps, int vel_type, int nparts, int vis, int cont_update,
                float duration, float v_inf) {
  shooter_xyz[0] = xyz[0];
  shooter_xyz[1] = xyz[1];
  shooter_xyz[2] = xyz[2];

  shooter_dxyz[0] = dxyz[0];
  shooter_dxyz[1] = dxyz[1];
  shooter_dxyz[2] = dxyz[2];

  shooter_uvw[0] = uvw[0];
  shooter_uvw[1] = uvw[1];
  shooter_uvw[2] = uvw[2];

  shooter_velmag = velmag;
  shooter_veldir = veldir;
  shooterpointsize = pointsize;

  shooter_fps = fps;
  shooter_vel_type = vel_type;
  shooter_nparts = nparts;
  visShooter = vis;
  shooter_cont_update = cont_update;

  shooter_duration = duration;
  shooter_v_inf = v_inf;

  return 0;
} // SHOOTER

int set_showdevices(int ndevices_ini, const char **names) {
  sv_object *obj_typei;
  char *dev_label;
  int i;
  char tempname[255]; // temporary buffer to convert from const string

  for(i = 0; i<nobject_defs; i++){
    obj_typei = object_defs[i];
    obj_typei->visible = 0;
  }
  for(i = 0; i<ndevices_ini; i++){
    strncpy(tempname, names[i], 255 - 1); // use temp buffer
    obj_typei = get_object(tempname);
    // obj_typei = get_object(names[i]);
    if(obj_typei != NULL){
      obj_typei->visible = 1;
    }
  }
  return 0;
} // SHOWDEVICES

int set_showdevicevals(int vshowdeviceval, int vshowvdeviceval,
    int vdevicetypes_index, int vcolordeviceval, int vvectortype,
    int vvispilot, int vshowdevicetype, int vshowdeviceunit) {
  showdeviceval = vshowdeviceval;
  showvdeviceval = vshowvdeviceval;
  devicetypes_index = vdevicetypes_index;
  colordeviceval = vcolordeviceval;
  vectortype = vvectortype;
  vispilot = vvispilot;
  showdevicetype = vshowdevicetype;
  showdeviceunit = vshowdeviceunit;
  devicetypes_index = CLAMP(vdevicetypes_index, 0, ndevicetypes - 1);
  update_glui_devices();
} // SHOWDEVICEVALS

int set_showmissingobjects(int v) {
  show_missing_objects = v;
  ONEORZERO(show_missing_objects);
  return 0;
} // SHOWMISSINGOBJECTS

int set_tourindex(int v) {
  selectedtour_index_ini = v;
  if(selectedtour_index_ini < 0)selectedtour_index_ini = -1;
  update_selectedtour_index = 1;
  return 0;
} // TOURINDEX

int set_userticks(int vis, int auto_place, int sub, float origin[],
                  float min[], float max[], float step[],
                  int show_x, int show_y, int show_z) {
  visUSERticks = vis;
  auto_user_tick_placement = auto_place;
  user_tick_sub = sub;

  user_tick_origin[0] = origin[0];
  user_tick_origin[1] = origin[1];
  user_tick_origin[2] = origin[2];

  user_tick_min[0] = min[0];
  user_tick_min[1] = min[1];
  user_tick_min[2] = min[2];

  user_tick_max[0] = max[0];
  user_tick_max[1] = max[1];
  user_tick_max[2] = max[2];

  user_tick_step[0] = step[0];
  user_tick_step[1] = step[1];
  user_tick_step[2] = step[2];

  user_tick_show_x = show_x;
  user_tick_show_y = show_y;
  user_tick_show_z = show_z;

  return 0;
} // USERTICKS

int set_c_particles(int minFlag, float minValue, int maxFlag, float maxValue,
                    const char *label) {
  if (label == NULL) {
   label = "";
  }
  int l = strlen(label);
  char label_copy[l+1];
  // convert to mutable string (mainly to avoid discard const warnings)
  strcpy(label_copy, label);
  if(npart5prop>0){
    int label_index = 0;
    if(strlen(label)>0)label_index = get_partprop_index_s(label_copy);
    if(label_index >= 0 && label_index<npart5prop){
      partpropdata *propi;

      propi = part5propinfo + label_index;
      propi->setchopmin = minFlag;
      propi->setchopmax = maxFlag;
      propi->chopmin = minValue;
      propi->chopmax = maxValue;
    }
  }
  return 0;
}

int set_c_slice(int minFlag, float minValue, int maxFlag, float maxValue,
                    const char *label) {
  int i;
  // if there is a label, use it
  if(strcmp(label, "") != 0){
    for(i = 0; i<nslice2; i++){
      if(strcmp(slicebounds[i].datalabel, label) != 0)continue;
      slicebounds[i].setchopmin = minFlag;
      slicebounds[i].setchopmax = maxFlag;
      slicebounds[i].chopmin = minValue;
      slicebounds[i].chopmax = maxValue;
      break;
    }
  // if there is no label apply values to all slice types
  } else{
    for(i = 0; i<nslice2; i++){
      slicebounds[i].setchopmin = minFlag;
      slicebounds[i].setchopmax = maxFlag;
      slicebounds[i].chopmin = minValue;
      slicebounds[i].chopmax = maxValue;
    }
  }
  return 0;
}

int set_cache_boundarydata(int setting) {
  cache_boundarydata = setting;
  return 0;
} // CACHE_BOUNDARYDATA

int set_cache_qdata(int setting) {
  cache_qdata = setting;
  return 0;
} // CACHE_QDATA

int set_percentilelevel(int setting) {
  percentile_level = setting;
  if(percentile_level<0.0)percentile_level = 0.01;
  if(percentile_level>0.5)percentile_level = 0.01;
  return 0;
} // PERCENTILELEVEL

int set_timeoffset(int setting) {
  timeoffset = setting;
  return 0;
} // TIMEOFFSET

int set_patchdataout(int outputFlag, float tmin, float tmax, float xmin,
                     float xmax, float ymin, float ymax, float zmin,
                     float zmax) {
  output_patchdata = outputFlag;
  patchout_tmin = tmin;
  patchout_tmax = tmax;
  patchout_xmin = xmin;
  patchout_xmax = xmax;
  patchout_ymin = ymin;
  patchout_ymax = ymax;
  patchout_zmin = zmin;
  patchout_zmax = zmax;
  ONEORZERO(output_patchdata);
  return 0;
} // PATCHDATAOUT

int set_c_plot3d(int n3d, int minFlags[], int minVals[], int maxFlags[],
                 int maxVals[]) {
  if(n3d>mxplot3dvars)n3d = mxplot3dvars;
  for (int i = 0; i<n3d; i++){
    setp3chopmin[i] = minFlags[i];
    setp3chopmax[i] = maxFlags[i];
    p3chopmin[i] = minVals[i];
    p3chopmax[i] = maxVals[i];
  }
  return 0;
} // C_PLOT3D

int set_v_plot3d(int n3d, int minFlags[], int minVals[], int maxFlags[],
                 int maxVals[]) {
  if(n3d>mxplot3dvars)n3d = mxplot3dvars;
  for (int i = 0; i<n3d; i++){
    setp3min[i] = minFlags[i];
    setp3max[i] = maxFlags[i];
    p3min[i] = minVals[i];
    p3max[i] = maxVals[i];
  }
  return 0;
} // V_PLOT3D

int set_tload(int beginFlag, float beginVal, int endFlag, int endVal,
              int skipFlag, int skipVal) {
  use_tload_begin = beginFlag;
  tload_begin = beginVal;
  use_tload_end = endFlag;
  tload_end = endVal;
  use_tload_skip = skipFlag;
  tload_skip = skipVal;
  return 0;
} // TLOAD

int set_v5_particles(int minFlag, float minValue, int maxFlag, float maxValue,
                    const char *label) {
  if (label == NULL) {
   label = "";
  }
  int l = strlen(label);
  char label_copy[l+1];
  // convert to mutable string (mainly to avoid discard const warnings)
  strcpy(label_copy, label);
  if(npart5prop>0){
    int label_index = 0;

    if(strlen(label)>0)label_index = get_partprop_index_s(label_copy);
    if(label_index >= 0 && label_index<npart5prop){
      partpropdata *propi;

      propi = part5propinfo + label_index;
      propi->setvalmin = minFlag;
      propi->setvalmax = maxFlag;
      propi->valmin = minValue;
      propi->valmax = maxValue;
      switch(minFlag){
        case PERCENTILE_MIN:
          propi->percentile_min = minValue;
          break;
        case GLOBAL_MIN:
          propi->global_min = minValue;
          break;
        case SET_MIN:
          propi->user_min = minValue;
          break;
        default:
          ASSERT(FFALSE);
          break;
      }
      switch(maxFlag){
        case PERCENTILE_MAX:
          propi->percentile_max = maxValue;
          break;
        case GLOBAL_MAX:
          propi->global_max = maxValue;
          break;
        case SET_MAX:
          propi->user_max = maxValue;
          break;
        default:
          ASSERT(FFALSE);
          break;
      }
    }
  }
  return 0;
}

int set_v_particles(int minFlag, float minValue, int maxFlag, float maxValue) {
  setpartmin = minFlag;
  partmin = minValue;
  setpartmax = maxFlag;
  partmax = maxValue;
  return 0;
} // V_PARTICLES

int set_v_target(int minFlag, float minValue, int maxFlag, float maxValue) {
  settargetmin = minFlag;
  targetmin = minValue;
  settargetmax = maxFlag;
  targetmax = maxValue;
  return 0;
} // V_TARGET

int set_v_slice(int minFlag, float minValue, int maxFlag, float maxValue,
                    const char *label, float lineMin, float lineMax,
                    int lineNum) {
  // if there is a label to apply, use it
  if(strcmp(label, "") != 0){
    for (int i = 0; i<nslice2; i++){
      if(strcmp(slicebounds[i].datalabel, label) != 0)continue;
      slicebounds[i].setvalmin = minFlag;
      slicebounds[i].setvalmax = maxFlag;
      slicebounds[i].valmin = minValue;
      slicebounds[i].valmax = maxValue;

      slicebounds[i].line_contour_min = lineMin;
      slicebounds[i].line_contour_max = lineMax;
      slicebounds[i].line_contour_num = lineNum;
      break;
    }
  // if there is no label apply values to all slice types
  } else{
    for (int i = 0; i<nslice2; i++){
      slicebounds[i].setvalmin = minFlag;
      slicebounds[i].setvalmax = maxFlag;
      slicebounds[i].valmin = minValue;
      slicebounds[i].valmax = maxValue;

      slicebounds[i].line_contour_min = lineMin;
      slicebounds[i].line_contour_max = lineMax;
      slicebounds[i].line_contour_num = lineNum;
    }
  }
  return 0;
} // V_SLICE

int show_smoke3d_showall() {
  smoke3ddata *smoke3di;

  updatemenu=1;
  glutPostRedisplay();
  plotstate=DYNAMIC_PLOTS;
  for (int i=0;i<nsmoke3dinfo;i++){
    smoke3di = smoke3dinfo + i;
    if(smoke3di->loaded==1)smoke3di->display=1;
  }
  glutPostRedisplay();
  Update_Show();
  return 0;
}

int show_smoke3d_hideall() {
  smoke3ddata *smoke3di;

  updatemenu=1;
  glutPostRedisplay();
  for (int i=0;i<nsmoke3dinfo;i++){
    smoke3di = smoke3dinfo + i;
    if(smoke3di->loaded==1)smoke3di->display=0;
  }
  Update_Show();
  return 0;
}

int show_slices_showall(){
  updatemenu=1;
  glutPostRedisplay();
  for (int i=0;i<nsliceinfo;i++){
    sliceinfo[i].display=1;
  }
  show_all_slices=1;
  updateslicefilenum();
  plotstate=getplotstate(DYNAMIC_PLOTS);

  updateglui();
  updateslicelistindex(slicefilenum);
  Update_Show();
  glutPostRedisplay();
  return 0;
}

int show_slices_hideall(){
  updatemenu=1;
  glutPostRedisplay();
  for (int i=0;i<nsliceinfo;i++){
    sliceinfo[i].display=0;
  }
  show_all_slices=0;
  updateslicefilenum();
  plotstate=getplotstate(DYNAMIC_PLOTS);

  updateglui();
  updateslicelistindex(slicefilenum);
  Update_Show();
  return 0;
}
