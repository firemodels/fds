#include "options.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include "smokeviewvars.h"

#include GLUT_H
#include "gd.h"

#define RENDER_START 3
//#define pp_RENDER360NEW

/* ------------------ does_movie_exist ------------------------ */

int does_movie_exist(char *mov_name, char *moviefile){
  char *movie;

  if(mov_name == NULL || strlen(mov_name) < 1)return 0;
  trim_back(mov_name);
  movie = trim_front(mov_name);
  strcpy(moviefile, movie);
  strcat(moviefile, ".mp4");
  if(file_exists(moviefile) == 1)return 1;
  return 0;
}

/* ------------------ PlayMovie ------------------------ */

void PlayMovie(void){
  char command_line[1024], moviefile_path[1024];

  if(play_movie_now==0)return;
  if(file_exists(get_moviefile_path(moviefile_path)) == 1){
    strcpy(command_line, "ffplay ");
    strcat(command_line,moviefile_path);
    psystem(command_line);
  }
  else{
    PRINTF("*** Error: the movie file, %s, does not exist\n", moviefile_path);
  }
}

/* ------------------ get_moviefile_path ------------------------ */

char *get_moviefile_path(char *moviefile_path){
  char moviefile[1024], *movie;

  trim_back(movie_name);
  movie = trim_front(movie_name);
  strcpy(moviefile, movie);
  strcat(moviefile, movie_ext);
  strcpy(moviefile_path, "");

  // if a script is running  prepend path defined in script

  if(script_dir_path != NULL&&strlen(script_dir_path) > 0){
    strcat(moviefile_path, script_dir_path);
    strcat(moviefile_path, dirseparator);
  }
  strcat(moviefile_path, moviefile);
  return moviefile_path;
}

/* ------------------ MakeMovie ------------------------ */

void MakeMovie(void){
  char command_line[1024];
  char frame0[1024];
  char moviefile_path[1024],overwrite_flag[10],image_ext[10], movie_frames[1024];
  int make_movie_now=1;

// wait to make movie until after images are rendered

  if(render_state == ON)return;

  if(renderfiletype==JPEG){
    strcpy(image_ext, ".jpg");
  }
  else{
    strcpy(image_ext, ".png");
  }

// if the first frame doesn't exist then generate images

  strcpy(frame0, render_file_base);
  strcat(frame0, "_0001");
  strcat(frame0, image_ext);
  if(runscript==0&&file_exists(frame0)==0){
    Render_CB(RENDER_START);
    return;
  }

// construct full pathname of movie

  get_moviefile_path(moviefile_path);

// add -y option if overwriting movie file

  if(overwrite_movie == 1){
    strcpy(overwrite_flag, "-y ");
  }
  else{
    strcpy(overwrite_flag, "");
    if(file_exists(moviefile_path) == 1&&script_dir_path==NULL){
       PRINTF("*** Warning: The movie file %s exists.  Set movie overwrite checkbox in movie dialog box.\n",moviefile_path);
       make_movie_now=0;
    }
  }


  if(make_movie_now==1){
// construct name of frames used to make movie

    strcpy(movie_frames, render_file_base);
    strcat(movie_frames,"_%04d");
    strcat(movie_frames, image_ext);

  // form command line for making movie

    sprintf(command_line, "ffmpeg %s -r %i -i ", overwrite_flag,movie_framerate);
    strcat(command_line, movie_frames);
    strcat(command_line, " ");
    strcat(command_line, moviefile_path);

// make movie

    system(command_line);
  }

// enable movie making button

  enable_disable_makemovie(ON);
  enable_disable_playmovie();

  update_makemovie = 0;
}

/* ------------------ Render ------------------------ */

void Render(int view_mode){
  if(current_script_command!=NULL&&(current_script_command->command==SCRIPT_VOLSMOKERENDERALL||current_script_command->command==SCRIPT_ISORENDERALL)){
    if( (render_frame[itimes]>0&&showstereo==STEREO_NONE)||(render_frame[itimes]>1&&showstereo!=STEREO_NONE) ){
      if(itimes==0){
        current_script_command->remove_frame=itimes;
        current_script_command->exit=1;
        stept=0;
        return;
      }
    }
  //  render_frame[itimes]++; //xxx check whether this is needed
    if( (render_frame[itimes]>0&&showstereo==STEREO_NONE)||(render_frame[itimes]>1&&showstereo!=STEREO_NONE) ){
      current_script_command->remove_frame=itimes;
    }
  }
  if(RenderOnceNow==0&&RenderOnceNowR==0&&RenderOnceNowL==0&&render_state==RENDER_ON&&render_multi==0){
    if(plotstate==DYNAMIC_PLOTS && nglobal_times>0){
     if(itimes>=0&&itimes<nglobal_times&&
       ((render_frame[itimes] == 0&&showstereo==STEREO_NONE)||(render_frame[itimes]<2&&showstereo!=STEREO_NONE))
       ){
       render_frame[itimes]++;
       RenderFrame(view_mode);
     }
     else{
       ASSERT(RenderSkip>0);
       RenderState(RENDER_OFF);
       RenderSkip=1;
     }
    }
    if(touring == 1 && nglobal_times == 0){
      if(rendertourcount % RenderSkip == 0)RenderFrame(view_mode);
      rendertourcount++;
      if(nglobal_times>0)tourangle_global += (float)(2.0*PI/((float)nglobal_times/(float)RenderSkip));
      if(nglobal_times==0)tourangle_global += (float)(2.0*PI/((float)maxtourframes/(float)RenderSkip));
      if(tourangle_global>2.0*PI){
        RenderState(RENDER_OFF);
        RenderSkip=1;
        tourangle_global=0.0;
      }
    }
  }

  if(render_multi==0){
    SNIFF_ERRORS("after render");
  }

  if(RenderOnceNow==1||RenderOnceNowL==1||RenderOnceNowR==1){
    if(render_multi==0)RenderFrame(view_mode);
    RenderOnceNow=0;
    if(view_mode==VIEW_LEFT)RenderOnceNowL=0;
    if(view_mode==VIEW_RIGHT)RenderOnceNowR=0;
    if(RenderOnceNowR==0&&RenderOnceNowL==0&&render_multi==0){
      RenderState(RENDER_OFF);
      RenderSkip=1;
    }
  }
  if(script_render==1){
    script_render=0;
    RenderState(RENDER_OFF);
  }
}

  /* ------------------ RenderFrame ------------------------ */

void RenderFrame(int view_mode){
  char renderfile_name[1024], renderfile_dir[1024], renderfile_full[1024], renderfile_suffix[1024], *renderfile_ext;
  char *renderfile_dir_ptr=NULL;
  int use_scriptfile;
  int woffset=0,hoffset=0;
  int screenH;

#ifdef WIN32
  SetThreadExecutionState(ES_DISPLAY_REQUIRED); // reset display idle timer to prevent screen saver from activating
#endif

  screenH = screenHeight;
  if(view_mode==VIEW_LEFT&&showstereo==STEREO_RB)return;

// construct filename for image to be rendered

  strcpy(renderfile_dir,"");
  strcpy(renderfile_suffix,"");
  use_scriptfile=0;

  // filename base

  if(current_script_command==NULL){
    strcpy(renderfile_name,render_file_base);
  }
  else{
    char suffix[20];

    if(
      (current_script_command->command==SCRIPT_RENDERONCE||
       current_script_command->command==SCRIPT_RENDERALL||
       current_script_command->command==SCRIPT_VOLSMOKERENDERALL||
       current_script_command->command==SCRIPT_ISORENDERALL
       )&&
       current_script_command->cval2!=NULL
       ){
        strcpy(renderfile_name,current_script_command->cval2);
        use_scriptfile=1;
    }
    else{
      strcpy(renderfile_name,fdsprefix);
    }
    if(script_dir_path!=NULL&&strlen(script_dir_path)>0){
      if(strlen(script_dir_path)==2&&script_dir_path[0]=='.'&&script_dir_path[1]==dirseparator[0]){
      }
      else{
        strcpy(renderfile_dir,script_dir_path);
        renderfile_dir_ptr=renderfile_dir;
      }
    }
    strcpy(suffix,"");
    switch(view_mode){
    case VIEW_LEFT:
      if(showstereo==STEREO_LR){
        strcat(suffix,"_L");
      }
      break;
    case VIEW_RIGHT:
      if(showstereo==STEREO_LR){
        strcat(suffix,"_R");
      }
      break;
    case VIEW_CENTER:
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
    strcat(renderfile_suffix,suffix);
  }

  // directory

  if(can_write_to_dir(renderfile_dir)==0){
    if(can_write_to_dir(smokeviewtempdir)==1){
      strcpy(renderfile_dir,smokeviewtempdir);
    }
    else{
      fprintf(stderr,"*** Error: unable to output render file\n");
      return;
    }
  }

  // filename suffix

  if(use_scriptfile==0||
    (current_script_command!=NULL&&
    (current_script_command->command==SCRIPT_RENDERALL||
     current_script_command->command==SCRIPT_VOLSMOKERENDERALL||
     current_script_command->command==SCRIPT_ISORENDERALL
     ))){
    int image_num;
    char suffix[20];

    strcpy(renderfile_suffix,"_");
    if(RenderTime==0){
      image_num=seqnum;
    }
    else{
      if(skip_render_frames==1){
        image_num=itimes;
      }
      else{
        image_num=itimes/RenderSkip;
      }
    }
    if(renderfilelabel==0||RenderTime==0){
      float time_local;
      int code;

      if(RenderTime==0){
        sprintf(suffix,"s%04i",image_num);
      }
      else{
        sprintf(suffix,"%04i",image_num);
      }
      code = getplot3dtime(&time_local);
      if(code==1&&renderfilelabel==1){
        char timelabel_local[20], *timelabelptr, dt=1.0;

        timelabelptr = time2timelabel(time_local,dt,timelabel_local);
        strcat(suffix,"_");
        strcat(suffix,timelabelptr);
        strcat(suffix,"s");
      }
    }
    else{
      float time_local;
      char timelabel_local[20], *timelabelptr;
      float dt;

      time_local = global_times[itimes];
      dt = global_times[1]-global_times[0];
      if(dt<0.0)dt=-dt;
      timelabelptr = time2timelabel(time_local,dt,timelabel_local);
      strcpy(suffix,timelabelptr);
      strcat(suffix,"s");
    }
    switch(view_mode){
    case VIEW_CENTER:
      if(RenderTime==0)seqnum++;
      break;
    case VIEW_LEFT:
        if(showstereo==STEREO_LR){
          strcat(suffix,"_L");
        }
      break;
    case VIEW_RIGHT:
      if(showstereo==STEREO_NONE||showstereo==STEREO_TIME||showstereo==STEREO_LR){
        strcat(suffix,"_R");
      }
      if(RenderTime==0)seqnum++;
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
    strcat(renderfile_suffix,suffix);
  }

  if(showstereo==STEREO_LR&&(view_mode==VIEW_LEFT||view_mode==VIEW_RIGHT)){
    hoffset=screenHeight/4;
    screenH = screenHeight/2;
    if(view_mode==VIEW_RIGHT)woffset=screenWidth;
  }

  // filename extension

  switch(renderfiletype){
  case 0:
    renderfile_ext=ext_png;
    break;
  case 1:
    renderfile_ext=ext_jpg;
    break;
  default:
    renderfiletype=2;
    renderfile_ext=ext_png;
    break;
  }

  // form full filename from parts

  strcpy(renderfile_full,renderfile_name);
  if(strlen(renderfile_suffix)>0)strcat(renderfile_full,renderfile_suffix);
  strcat(renderfile_full,renderfile_ext);

  // render image

  SVimage2file(renderfile_dir_ptr,renderfile_full,renderfiletype,woffset,screenWidth,hoffset,screenH);
  if(RenderTime==1&&output_slicedata==1){
    output_Slicedata();
  }
}

/* ------------------ getscreenbuffer --------- */

GLubyte *getscreenbuffer(void){

  GLubyte *OpenGLimage=NULL;

  int x=0, y=0;

  NewMemory((void **)&OpenGLimage,screenWidth * screenHeight * sizeof(GLubyte) * 3);

  if(OpenGLimage==NULL)return NULL;

  glPixelStorei(GL_PACK_ALIGNMENT, 1);

  /* get the image from the OpenGL frame buffer */

  glReadPixels(x, y, screenWidth, screenHeight, GL_RGB, GL_UNSIGNED_BYTE, OpenGLimage);

  return OpenGLimage;

}

/* ------------------ mergescreenbuffers ------------------------ */

int mergescreenbuffers(int nscreen_rows, GLubyte **screenbuffers){

  char *renderfile,renderfile_base[1024],renderfile2[1024];
  char *ext;
  FILE *RENDERfile=NULL;
  GLubyte *p;
  gdImagePtr RENDERimage;
  unsigned int r, g, b;
  int i,j,rgb_local;
  int nscreen_cols;
  int irow;

  nscreen_cols=nscreen_rows;
  switch(renderfiletype){
  case PNG:
    ext=ext_png;
    break;
  case JPEG:
    ext=ext_jpg;
    break;
  default:
    renderfiletype=0;
    ext=ext_png;
    break;
  }

  if(scriptoutstream!=NULL&&current_script_command!=NULL&&current_script_command->cval2!=NULL){
    strcpy(renderfile2,current_script_command->cval2);
  }
  else{
    strcpy(renderfile2,fdsprefix);
  }
  if(RenderTime==1){
      sprintf(renderfile_base,"%s_%04i",renderfile2,itimes/RenderSkip);
  }
  if(RenderTime==0){
      sprintf(renderfile_base,"%s_s%04i",renderfile2,seqnum);
      seqnum++;
  }
  strcat(renderfile_base,ext);
  renderfile=get_filename(smokeviewtempdir,renderfile_base,tempdir_flag);
  if(renderfile==NULL){
    fprintf(stderr,"*** Error: unable to write to %s",renderfile_base);
    return 1;
  }
  RENDERfile = fopen(renderfile, "wb");
  PRINTF("Rendering to: %s .",renderfile);
  if(RENDERfile == NULL){
    fprintf(stderr,"*** Error: unable to write to %s",renderfile);
    FREEMEMORY(renderfile);
    return 1;
  }
  RENDERimage = gdImageCreateTrueColor(nscreen_cols*screenWidth,nscreen_rows*screenHeight);

  p = *screenbuffers++;
  for(irow=0;irow<nscreen_rows;irow++){
    int icol;

    for(icol=0;icol<nscreen_cols;icol++){
      for (i = (nscreen_rows-irow)*screenHeight-1 ; i>=(nscreen_rows-irow-1)*screenHeight; i--){
        for(j=icol*screenWidth;j<(icol+1)*screenWidth;j++){
          r=*p++; g=*p++; b=*p++;
          rgb_local = (r<<16)|(g<<8)|b;
          gdImageSetPixel(RENDERimage,j,i,rgb_local);

        }
      }
      p=*screenbuffers++;
    }
  }

  /* output the image */

  switch(renderfiletype){
  case PNG:
    gdImagePng(RENDERimage,RENDERfile);
    break;
  case JPEG:
    gdImageJpeg(RENDERimage,RENDERfile,-1);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  fclose(RENDERfile);

  /* free up memory used by both OpenGL and GIF images */

  gdImageDestroy(RENDERimage);
  FREEMEMORY(renderfile);
  PRINTF(" Completed\n");
  return 0;
}

/* ------------------ getscreenmap360 ------------------------ */

unsigned int getscreenmap360(float *xyz) {
  int ibuff, i, j;

  for (ibuff = 0; ibuff < nscreeninfo; ibuff++){
    screendata *screeni;
    float *view, *up, *right, t;
    float A, B;

    screeni = screeninfo + ibuff;
    view = screeni->view;
    up = screeni->up;
    right = screeni->right;

    t = DOT3(xyz,view);
//    if (t < 0.0)continue;

    A = DOT3(xyz, right)/t;
    B = DOT3(xyz, up)/t;
    
    if(ABS(A)<=screeni->width/2.0 && ABS(B)<=screeni->height/2.0){
      int ix, iy, index;
      unsigned int return_val;

      ix = CLAMP(screeni->nwidth*(screeni->width / 2.0 + A) / screeni->width,    0, screeni->nwidth - 1);
      iy = CLAMP(screeni->nheight*(screeni->height / 2.0 + B) / screeni->height, 0, screeni->nheight - 1);
      index = iy*screeni->nwidth + ix;
      return_val = ((ibuff+1) << 24) |  index;
      return return_val;
    }
  }
  return 0;
}

/* ------------------ setup_screeninfo ------------------------ */

void setup_screeninfo(void){
  int ibuf;

  nscreeninfo = 26;
  NewMemory((void **)&screeninfo, nscreeninfo * sizeof(screendata));

  for(ibuf = 0; ibuf < nscreeninfo; ibuf++){
    screendata *screeni;
    float azimuth, elevation;
    float *right, *view, *up;
    float sina, cosa;
    float cose, sine;
    float aspect_ratio;
    int j, *map;
    float aperture_width, aperture_height;

    aperture_width = 45.0;
    screeni = screeninfo + ibuf;
    screeni->map = NULL;
    screeni->nwidth=VP_scene.width;
    screeni->nheight=VP_scene.height;
    aspect_ratio = (float)screeni->nwidth/(float)screeni->nheight;
    screeni->width=2.0*tan(DEG2RAD*aperture_width/2.0);
    screeni->height = screeni->width / aspect_ratio;
    aperture_height = 2.0*RAD2DEG*atan(screeni->height / 2.0);

    if(ibuf == 0){
      azimuth = 0.0;
      elevation = -90;
    }
    else if(ibuf >= 1 && ibuf < 9){
      int ii;

      ii = ibuf - 1;
      azimuth = ii*45.0;
      elevation = -90.0+aperture_height;
    }
    else if(ibuf >= 9 && ibuf < 17){
      int ii;

      ii = ibuf - 9;
      azimuth = ii*45.0;
      elevation = -90.0 + 2.0*aperture_height;
    }
    else if(ibuf >= 17 && ibuf < 25){
      int ii;

      ii = ibuf - 17;
      azimuth = ii*45.0;
      elevation = 45.0;
      elevation = -90.0 + 3.0*aperture_height;
    }
    else if(ibuf == 25){
      azimuth = 0.0;
      elevation = -90.0 + 4.0*aperture_height;
    }

    cosa = cos(DEG2RAD*azimuth);
    sina = sin(DEG2RAD*azimuth);
    cose = cos(DEG2RAD*elevation);
    sine = sin(DEG2RAD*elevation);

    view = screeni->view;
    view[0] = sina*cose;
    view[1] = cosa*cose;
    view[2] = sine;

    up = screeni->up;
    if(ABS(sine) < 0.001){
      up[0] = 0.0;
      up[1] = 0.0;
      up[2] = 1.0;
    }
    else {
      float denom;

      denom = sqrt(1.0 + cose*cose/(sine*sine));
      up[0] = -sina / denom;
      up[1] = -cosa / denom;
      up[2] = (cose / sine) / denom;
    }

    right = screeni->right;
    if(sine < 0.0){
      right[0] = -cosa;
      right[1] = sina;
      right[2] = 0.0;
    }
    else {
      right[0] = cosa;
      right[1] = -sina;
      right[2] = 0.0;
    }

    NewMemory((void **)&map, screeni->nheight*screeni->nwidth * sizeof(int));
    screeni->map = map;

    for(j = 0; j < screeni->nheight; j++){
      int i;
      float facty;

      facty = screeni->height*(float)(j - screeni->nheight / 2)/ (float)screeni->nheight;

      for(i = 0; i < screeni->nwidth; i++){
        float factx, xyz[3], xyznorm;
        float elev, az;
        int ijk, ijk360;
        int i360, j360;

        factx = screeni->width*(float)(i - screeni->nwidth / 2) / (float)screeni->nwidth;

        xyz[0] = view[0] + factx*right[0] + facty*up[0];
        xyz[1] = view[1] + factx*right[1] + facty*up[1];
        xyz[2] = view[2] + factx*right[2] + facty*up[2];
        xyznorm = sqrt(xyz[0] * xyz[0] + xyz[1] * xyz[1] + xyz[2] * xyz[2]);

        elev = CLAMP(asin(xyz[2] / xyznorm)*RAD2DEG, -90.0, 90.0);
        az = CLAMP(atan2(xyz[0], xyz[1])*RAD2DEG,-180.0,180.0);

        ijk = j*screeni->nwidth + i;
        j360 = CLAMP(nheight360-1-nheight360*(elev + 90.0) / 180.0,0,nheight360-1);
        i360 = CLAMP(nwidth360*(az + 180.0) / 360.0,0,nwidth360-1);
        ijk360 = j360*nwidth360 + i360;
        map[ijk] = ijk360;
      }
    }
  }

  NewMemory((void **)&screenmap360, nwidth360*nheight360 * sizeof(unsigned int));
  {
    int i,j;
    float *cosa, *sina, *cose, *sine;

    NewMemory((void **)&sina, nwidth360 * sizeof(float));
    NewMemory((void **)&cosa, nwidth360 * sizeof(float));
    NewMemory((void **)&sine, nheight360 * sizeof(float));
    NewMemory((void **)&cose, nheight360 * sizeof(float));

    for (i = 0; i < nwidth360; i++){
      float alpha;

      alpha = -180.0 + (float)i*360.0 / (float)nwidth360;
      sina[i] = sin(DEG2RAD*alpha);
      cosa[i] = cos(DEG2RAD*alpha);
    }
    for (i = 0; i < nheight360; i++){
      float eps;

      eps = -90.0 + (float)i*180.0 / (float)nheight360;
      sine[i] = sin(DEG2RAD*eps);
      cose[i] = cos(DEG2RAD*eps);
    }
    for (j = 0; j < nheight360; j++){
      for (i = 0; i < nwidth360; i++){
        float xyz[3];
        int buff;

        xyz[0] = sina[i] * cose[j];
        xyz[1] = cosa[i] * cose[j];
        xyz[2] = sine[j];
        screenmap360[j*nwidth360 + i] = getscreenmap360(xyz);
      }
    }
    FREEMEMORY(sina);
    FREEMEMORY(cosa);
    FREEMEMORY(sine);
    FREEMEMORY(cose);
  }
}

/* ------------------ mergescreenbuffers360 ------------------------ */

int mergescreenbuffers360(void){

  char *renderfile, renderfile_base[1024], renderfile2[1024];
  char *ext;
  FILE *RENDERfile = NULL;
  gdImagePtr RENDERimage;
#ifndef pp_RENDER360NEW
  int ibuff;
#endif
  int i, j, ijk360;
  int *screenbuffer;

  switch (renderfiletype){
  case PNG:
    ext = ext_png;
    break;
  case JPEG:
    ext = ext_jpg;
    break;
  default:
    renderfiletype = 0;
    ext = ext_png;
    break;
  }

  if(scriptoutstream != NULL&&current_script_command != NULL&&current_script_command->cval2 != NULL){
    strcpy(renderfile2, current_script_command->cval2);
  }
  else {
    strcpy(renderfile2, fdsprefix);
  }
  if(RenderTime == 1){
    sprintf(renderfile_base, "%s_%04i", renderfile2, itimes / RenderSkip);
  }
  if(RenderTime == 0){
    sprintf(renderfile_base, "%s_s%04i", renderfile2, seqnum);
    seqnum++;
  }
  strcat(renderfile_base, ext);
  renderfile = get_filename(smokeviewtempdir, renderfile_base, tempdir_flag);
  if(renderfile == NULL){
    fprintf(stderr, "*** Error: unable to write to %s", renderfile_base);
    return 1;
  }
  RENDERfile = fopen(renderfile, "wb");
  PRINTF("Rendering to: %s .", renderfile);
  if(RENDERfile == NULL){
    fprintf(stderr, "*** Error: unable to write to %s", renderfile);
    FREEMEMORY(renderfile);
    return 1;
  }
  RENDERimage = gdImageCreateTrueColor(nwidth360, nheight360);
  NewMemory((void **)&screenbuffer,nwidth360*nheight360 * sizeof(int));

  for(i=0;i<nwidth360*nheight360;i++){
    screenbuffer[i]=0;
  }

#ifdef pp_RENDER360NEW
  ijk360 = 0;
  for(j=0;j<nheight360;j++){
    for(i=0;i<nwidth360;i++){
      GLubyte *p;
      int ibuff, ijk, rgb_local;
      screendata *screeni;
      unsigned int r, g, b;

      ibuff = screenmap360[ijk360] >> 24;
      if(ibuff == 0)continue;
      ibuff--;
      screeni = screeninfo + ibuff;
      p = screeni->screenbuffer;

      ijk = screenmap360[ijk360] & 0xffffff;
      r=p[3*ijk];
      g=p[3*ijk+1];
      b=p[3*ijk+2];
      rgb_local = (r<<16)|(g<<8)|b;
      screenbuffer[ijk360]=rgb_local;
      ijk360++;
    }
  }
#else
  for(ibuff=0;ibuff<nscreeninfo;ibuff++){
    GLubyte *p;
    screendata *screeni;
    int *map;

#ifdef pp_RENDER360_DEBUG
    if(screenvis[ibuff] == 0)continue;
#endif
    screeni = screeninfo + ibuff;
    map = screeni->map;

    p = screeni->screenbuffer;

    for(j=0;j<screenHeight;j++){
      for(i=0;i<screenWidth;i++){
        unsigned int r, g, b;
        int rgb_local;
        int ijk;

        ijk = j*screenWidth + i;
        r=p[3*ijk];
        g=p[3*ijk+1];
        b=p[3*ijk+2];
        rgb_local = (r<<16)|(g<<8)|b;
        ijk360 = map[ijk];
        screenbuffer[ijk360]=rgb_local;
      }
    }
  }
#endif

  ijk360 = 0;
  for(j=0;j<nheight360;j++){
    for(i=0;i<nwidth360;i++){
      gdImageSetPixel(RENDERimage, i, j, screenbuffer[ijk360++]);
    }
  }

  /* output the image */

  switch (renderfiletype){
  case PNG:
    gdImagePng(RENDERimage, RENDERfile);
    break;
  case JPEG:
    gdImageJpeg(RENDERimage, RENDERfile, -1);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  fclose(RENDERfile);

  /* free up memory used by both OpenGL and GIF images */

  gdImageDestroy(RENDERimage);
  FREEMEMORY(renderfile);
  FREEMEMORY(screenbuffer);
  PRINTF(" Completed\n");
  return 0;
}

/* ------------------ SVimage2file ------------------------ */

int SVimage2file(char *directory, char *RENDERfilename, int rendertype, int woffset, int width, int hoffset, int height){

  FILE *RENDERfile;
  char *renderfile=NULL;
  GLubyte *OpenGLimage, *p;
  gdImagePtr RENDERimage;
  unsigned int r, g, b;
  int i,j,rgb_local;
  int width_beg, width_end, height_beg, height_end;
  int width2, height2;

  width_beg=woffset;
  width_end=width+woffset;
  height_beg=hoffset;
  height_end=hoffset+height;
  if(clip_rendered_scene==1){
    width_beg+=render_clip_left;
    width_end-=render_clip_right;
    height_beg+=render_clip_bottom;
    height_end-=render_clip_top;
  }
  width2 = width_end-width_beg;
  height2 = height_end-height_beg;

  if(directory==NULL){
    renderfile=get_filename(smokeviewtempdir,RENDERfilename,tempdir_flag);
  }
  else{
    renderfile=get_filename(directory,RENDERfilename,1);
  }
  if(renderfile == NULL){
    fprintf(stderr,"*** Error: Unable to write to %s\n",RENDERfilename);
    return 1;
  }
  RENDERfile = fopen(renderfile, "wb");
  if(RENDERfile == NULL){
    fprintf(stderr,"*** Error: Unable to write to %s\n",renderfile);
    return 1;
  }
  NewMemory((void **)&OpenGLimage,width2 * height2 * sizeof(GLubyte) * 3);
  if(OpenGLimage == NULL){
    fprintf(stderr,"*** Error allocating memory buffer for render file:%s\n",renderfile);
    fclose(RENDERfile);
    return 1;
  }
  PRINTF("Rendering to: %s .",renderfile);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);

  /* get the image from the OpenGL frame buffer */

  glReadPixels(width_beg, height_beg, width2, height2, GL_RGB, GL_UNSIGNED_BYTE, OpenGLimage);

  /* copy the image from OpenGL memory to GIF memory */

  p = OpenGLimage;

  RENDERimage = gdImageCreateTrueColor(width2,height2);

  for (i = height2-1 ; i>=0; i--){
    for(j=0;j<width2;j++){
      r=*p++; g=*p++; b=*p++;
      rgb_local = (r<<16)|(g<<8)|b;
      gdImageSetPixel(RENDERimage,j,i,rgb_local);
    }
  }

  if(test_smokesensors==1&&active_smokesensors==1&&show_smokesensors!=SMOKESENSORS_HIDDEN){
    int idev;

    for(idev=0;idev<ndeviceinfo;idev++){
      devicedata *devicei;
      int idev_col, idev_row;
      int col_offset, row_offset;
      unsigned int red=255<<16;

      devicei = deviceinfo + idev;

      if(devicei->object->visible==0)continue;
      if(strcmp(devicei->object->label,"smokesensor")!=0)continue;
      idev_row = devicei->screenijk[0];
      idev_col = devicei->screenijk[1];
      for(col_offset=-3;col_offset<4;col_offset++){
        for(row_offset=-3;row_offset<4;row_offset++){
          int irow, icol;

          irow = idev_row+row_offset;
          if(irow<0)irow=0;
          if(irow>width-1)irow=width-1;

          icol = height - 1 - (idev_col+col_offset);
          if(icol<0)icol=0;
          if(icol>height-1)icol=height-1;

          gdImageSetPixel(RENDERimage,irow,icol,red);
        }
      }
    }
  }

  /* output the gif image */

  switch(rendertype){
  case PNG:
    gdImagePng(RENDERimage,RENDERfile);
    break;
  case JPEG:
    gdImageJpeg(RENDERimage,RENDERfile,-1);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }

  /* free up memory used by both OpenGL and GIF images */

  fclose(RENDERfile);
  FREEMEMORY(renderfile);

  gdImageDestroy(RENDERimage);
  FREEMEMORY(OpenGLimage);
  PRINTF(" Completed.\n");
  return 0;
}

/* ------------------ SVimage2var ------------------------ */
#ifdef pp_LUA
int SVimage2var(int rendertype,
    int woffset, int width, int hoffset, int height, gdImagePtr *RENDERimage) {

  GLubyte *OpenGLimage, *p;
  unsigned int r, g, b;
  int i,j,rgb_local;
  int width_beg, width_end, height_beg, height_end;
  int width2, height2;

  width_beg=woffset;
  width_end=width+woffset;
  height_beg=hoffset;
  height_end=hoffset+height;
  if(clip_rendered_scene==1){
    width_beg+=render_clip_left;
    width_end-=render_clip_right;
    height_beg+=render_clip_bottom;
    height_end-=render_clip_top;
  }
  width2 = width_end-width_beg;
  height2 = height_end-height_beg;

  NewMemory((void **)&OpenGLimage,width2 * height2 * sizeof(GLubyte) * 3);
  if(OpenGLimage == NULL){
    fprintf(stderr,"*** Error allocating memory buffer for render var\n");
    return 1;
  }
  glPixelStorei(GL_PACK_ALIGNMENT, 1);

  /* get the image from the OpenGL frame buffer */

  glReadPixels(width_beg, height_beg, width2, height2, GL_RGB, GL_UNSIGNED_BYTE, OpenGLimage);

  /* copy the image from OpenGL memory to GIF memory */

  p = OpenGLimage;

  *RENDERimage = gdImageCreateTrueColor(width2,height2);

  for (i = height2-1 ; i>=0; i--){
    for(j=0;j<width2;j++){
      r=*p++; g=*p++; b=*p++;
      rgb_local = (r<<16)|(g<<8)|b;
      gdImageSetPixel(*RENDERimage,j,i,rgb_local);
    }
  }
  if(test_smokesensors==1&&active_smokesensors==1&&show_smokesensors!=SMOKESENSORS_HIDDEN){
    int idev;

    for(idev=0;idev<ndeviceinfo;idev++){
      devicedata *devicei;
      int idev_col, idev_row;
      int col_offset, row_offset;
      unsigned int red=255<<16;

      devicei = deviceinfo + idev;

      if(devicei->object->visible==0)continue;
      if(strcmp(devicei->object->label,"smokesensor")!=0)continue;
      idev_row = devicei->screenijk[0];
      idev_col = devicei->screenijk[1];
      for(col_offset=-3;col_offset<4;col_offset++){
        for(row_offset=-3;row_offset<4;row_offset++){
          int irow, icol;

          irow = idev_row+row_offset;
          if(irow<0)irow=0;
          if(irow>width-1)irow=width-1;

          icol = height - 1 - (idev_col+col_offset);
          if(icol<0)icol=0;
          if(icol>height-1)icol=height-1;

          gdImageSetPixel(*RENDERimage,irow,icol,red);
        }
      }
    }
  }

  /* free up memory used by OpenGL image */
  FREEMEMORY(OpenGLimage);
  PRINTF(" Completed.\n");
  return 0;
}
#endif

/* ------------------ readpicture ------------------------ */

unsigned char *readpicture(char *filename, int *width, int *height, int printflag){
  char *ext;
  unsigned char *returncode;
  char *filebuffer=NULL;
  int allocated;
  STRUCTSTAT statbuffer;

  if(filename==NULL)return NULL;
  if(STAT(filename,&statbuffer)==0){
    filebuffer=filename;
    allocated=0;
  }
  else{
    size_t lenbuffer;

    if(texturedir==NULL){
      if(printflag==1){
        fprintf(stderr,"*** Error: texture file: %s unavailable\n",filename);
      }
      return NULL;
    }
    else{
      FILE *stream;

      lenbuffer=strlen(filename)+strlen(texturedir)+1;
      NewMemory((void **)&filebuffer,(unsigned int)(lenbuffer+1));
      allocated=1;
      strcpy(filebuffer,texturedir);
      strcat(filebuffer,dirseparator);
      strcat(filebuffer,filename);
      stream=fopen(filebuffer,"rb");
      if(stream==NULL){
        if(printflag==1){
          fprintf(stderr,"*** Error: texture file: %s unavailable\n",filebuffer);
        }
        FREEMEMORY(filebuffer);
        return NULL;
      }
      else{
        fclose(stream);
      }
    }
  }


  if(printflag==1)PRINTF("Loading texture:%s ",filebuffer);
  ext = filebuffer + strlen(filebuffer) - 4;
  if(strncmp(ext,".jpg",4)==0||strncmp(ext,".JPG",4)==0){
    returncode = readjpeg(filebuffer,width,height,pixel_skip);
  }
  else if(strncmp(ext,".png",4)==0||strncmp(ext,".PNG",4)==0){
    returncode = readpng(filebuffer,width,height);
  }
  else{
    if(allocated==1){
      FREEMEMORY(filebuffer);
    }
    return NULL;
  }
  if(allocated==1){
    FREEMEMORY(filebuffer);
  }
  if(printflag==1){
    if(returncode!=NULL){
      PRINTF(" - completed\n");
    }
    else{
      PRINTF(" - failed\n");
      fprintf(stderr,"*** Error: attempt to input %s failed\n",filename);
    }
  }
  return returncode;

}

/* ------------------ readjpeg ------------------------ */

unsigned char *readjpeg(const char *filename,int *width, int *height, int skip_local){

  FILE *file;
  gdImagePtr image;
  unsigned char *dataptr,*dptr;
  int i,j;
  unsigned int intrgb;
  int WIDTH, HEIGHT;
  int NEWWIDTH, NEWHEIGHT;
  int jump;

  file = fopen(filename, "rb");
  if(file == NULL)return NULL;
  image = gdImageCreateFromJpeg(file);
  fclose(file);
  if(image==NULL)return NULL;
  WIDTH=gdImageSX(image);
  HEIGHT=gdImageSY(image);
  if(skip_local<0)skip_local=0;
  jump = skip_local + 1;
  NEWWIDTH=WIDTH/jump;
  if(WIDTH%jump!=0)NEWWIDTH++;
  NEWHEIGHT=HEIGHT/jump;
  if(HEIGHT%jump!=0)NEWHEIGHT++;
  *width=NEWWIDTH;
  *height=NEWHEIGHT;
  if( NewMemory((void **)&dataptr,(unsigned int)(4*NEWWIDTH*NEWHEIGHT) )==0){
    gdImageDestroy(image);
    return NULL;
  }
  dptr=dataptr;
  for (i = 0; i<HEIGHT; i+=jump){
    for(j=0;j<WIDTH;j+=jump){
      intrgb=(unsigned int)gdImageGetPixel(image,j,(unsigned int)(HEIGHT-(1+i)));
      *dptr++ = (intrgb>>16)&255;
      *dptr++ = (intrgb>>8)&255;
      *dptr++ = intrgb&255;
      *dptr++=0xff;
    }
  }
  gdImageDestroy(image);
  return dataptr;

}

/* ------------------ readpng ------------------------ */

unsigned char *readpng(const char *filename,int *width, int *height){

  FILE *file;
  gdImagePtr image;
  unsigned char *dataptr,*dptr;
  int i,j;
  unsigned int intrgb;

  file = fopen(filename, "rb");
  if(file == NULL)return NULL;
  image = gdImageCreateFromPng(file);
  fclose(file);
  *width=gdImageSX(image);
  *height=gdImageSY(image);
  if( NewMemory((void **)&dataptr,(unsigned int)(4*(*width)*(*height)) )==0){
    gdImageDestroy(image);
    return NULL;
  }
  dptr=dataptr;
  for (i = 0; i<*height; i++){
    for(j=0;j<*width;j++){
      intrgb=(unsigned int)gdImageGetPixel(image,j,(unsigned int)(*height-(1+i)));
      *dptr++ = (intrgb>>16)&255;
      *dptr++ = (intrgb>>8)&255;
      *dptr++ = intrgb&255;
      *dptr++ = 0xff;
    }
  }
  gdImageDestroy(image);
  return dataptr;

}
