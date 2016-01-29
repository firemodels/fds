#include "options.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "smokeviewvars.h"

#include GLUT_H
#include "gd.h"

#if defined(WIN32)
#include <windows.h>
#elif defined(_MINGW32_)
#include <sys/stat.h>
#else
#include <sys/stat.h>
#endif

#define RENDER_TYPE 0
#define RENDER_START 3

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
       RenderFrame(view_mode, NULL);
     }
     else{
       ASSERT(RenderSkip>0);
       RenderState(RENDER_OFF);
       RenderSkip=1;
     }
    }
    if(touring == 1 && nglobal_times == 0){
      if(rendertourcount % RenderSkip == 0)RenderFrame(view_mode, NULL);
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
    if(render_multi==0)RenderFrame(view_mode, NULL); 
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

char* form_filename(int view_mode, char *renderfile_name, char *renderfile_dir,
                   char *renderfile_path, int woffset, int hoffset, int screenH,
                   char *basename);
  /* ------------------ RenderFrame ------------------------ */
// The second argument to RenderFrame is the name that should be given to the
// rendered file. If basename == NULL, then a default filename is formed based
// on the chosen frame and rendering options.
void RenderFrame(int view_mode, char *basename){
  char renderfile_name[1024]; // the name the file (including extension)
  char renderfile_dir[1024]; // the directory into which the image will be rendered
  char renderfile_path[2048]; // the full path of the rendered image
  int woffset=0,hoffset=0;
  int screenH;

  if(script_dir_path != NULL){
      strcpy(renderfile_dir, script_dir_path);
  } else {
      strcpy(renderfile_dir, ".");
  }

#ifdef WIN32
  SetThreadExecutionState(ES_DISPLAY_REQUIRED); // reset display idle timer to prevent screen saver from activating
#endif

  screenH = screenHeight;
  if(view_mode==VIEW_LEFT&&showstereo==STEREO_RB)return;
  // construct filename for image to be rendered
  form_filename(view_mode, renderfile_name, renderfile_dir, renderfile_path,
                woffset, hoffset, screenH, basename);
                
  printf("renderfile_name: %s\n", renderfile_name);
  // render image
  SVimage2file(renderfile_dir,renderfile_name,renderfiletype,woffset,screenWidth,hoffset,screenH);
  if(RenderTime==1&&output_slicedata==1){
    output_Slicedata();
  }
}

// construct filepath for image to be renderd
char* form_filename(int view_mode, char *renderfile_name, char *renderfile_dir,
                   char *renderfile_path, int woffset, int hoffset, int screenH,
                   char *basename) {
    
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
                mkdir(renderfile_dir);
            }
        }
        if(showstereo==STEREO_LR&&(view_mode==VIEW_LEFT||view_mode==VIEW_RIGHT)){
          hoffset=screenHeight/4;
          screenH = screenHeight/2;
          if(view_mode==VIEW_RIGHT)woffset=screenWidth;
        }

        snprintf(renderfile_name, 1024,
                  "%s%s%s",
                  fdsprefix, view_suffix, renderfile_ext);
        printf("directory is: %s\n", renderfile_dir);
        printf("filename is: %s\n", renderfile_name);
    } else {
        snprintf(renderfile_name, 1024, "%s%s", basename, renderfile_ext);
    }
    return renderfile_name;
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

/* ------------------ SVimage2file ------------------------ */

int SVimage2file(char *directory, char *RENDERfilename, int rendertype, int woffset, int width, int hoffset, int height){
    printf("svimage2file: %s\n", RENDERfilename);
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

  if(test_smokesensors==1&&active_smokesensors==1&&show_smokesensors!=0){
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
