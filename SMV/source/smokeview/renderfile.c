// $Date$ 
// $Revision$
// $Author$

// svn revision character string
char renderfile_revision[]="$Revision$";

#include "options.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "string_util.h"
#include "smokeviewvars.h"

#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "gd.h"

#define PNG 0
#ifdef pp_JPEG
#define JPEG 1
#endif
#ifdef pp_GDGIF
#define GIF 2
#endif

/* ------------------ Render ------------------------ */

void Render(int view_mode){
  if(current_script_command!=NULL&&current_script_command->command==SCRIPT_VOLSMOKERENDERALL){
    if( (render_frame[itimes]>0&&showstereo==0)||(render_frame[itimes]>1&&showstereo!=0) ){
      if(itimes==0){
        current_script_command->remove_frame=itimes;
        current_script_command->exit=1;
        stept=0;
        return;
      }
    }
    render_frame[itimes]++;
    if( (render_frame[itimes]>0&&showstereo==0)||(render_frame[itimes]>1&&showstereo!=0) ){
      current_script_command->remove_frame=itimes;
    }
  }
  if(RenderOnceNow==0&&RenderGif !=0
    &&render_double==0
    ){
    if(plotstate==DYNAMIC_PLOTS && nglobal_times>0){
     if(itimes>=0&&itimes<nglobal_times&&
       ((render_frame[itimes] == 0&&showstereo==0)||(render_frame[itimes]<2&&showstereo!=0))
       ){
       render_frame[itimes]++;
       RenderFrame(view_mode);
     }
     else{
       ASSERT(RenderSkip>0);
       RenderState(0);
       RenderSkip=1;
     }
    }
    if(touring == 1 && nglobal_times == 0){
      if(rendertourcount % RenderSkip == 0)RenderFrame(view_mode);
      rendertourcount++;
      if(nglobal_times>0)tourangle_global += (float)(2.0*PI/((float)nglobal_times/(float)RenderSkip));
      if(nglobal_times==0)tourangle_global += (float)(2.0*PI/((float)maxtourframes/(float)RenderSkip));
      if(tourangle_global>2.0*PI){
        RenderState(0);
        RenderSkip=1;
        tourangle_global=0.0;
      }
    }
  }

  if(render_double==0){
    SNIFF_ERRORS("after render");
  }

  if(RenderOnceNow==1||RenderOnceNowL==1||RenderOnceNowR==1){
    if(render_double==0)RenderFrame(view_mode); 
    RenderOnceNow=0;
    if(view_mode==VIEW_LEFT)RenderOnceNowL=0;
    if(view_mode==VIEW_RIGHT)RenderOnceNowR=0;
    if(RenderOnceNowR==0&&RenderOnceNowL==0){
      RenderState(0);
      RenderSkip=1;
    }
  }
}

  /* ------------------ RenderFrame ------------------------ */

void RenderFrame(int view_mode){
  char renderfile_name[1024], renderfile_dir[1024], renderfile_full[1024], renderfile_suffix[1024], *renderfile_ext;
  int use_scriptfile;

  if(view_mode==VIEW_LEFT&&(showstereo==2||showstereo==3))return;

// construct filename for image to be rendered

  strcpy(renderfile_dir,"");
  strcpy(renderfile_suffix,"");
  use_scriptfile=0;

  // filename base

  if(current_script_command==NULL){
    strcpy(renderfile_name,fdsprefix);
  }
  else{
    if(
      (current_script_command->command==SCRIPT_RENDERONCE||
       current_script_command->command==SCRIPT_RENDERALL||
       current_script_command->command==SCRIPT_VOLSMOKERENDERALL
       )&&
       current_script_command->cval!=NULL
       ){
        strcpy(renderfile_name,current_script_command->cval);
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
      }
    }
  }
  
  // directory

  if(can_write_to_dir(renderfile_dir)==0){
    if(can_write_to_dir(smokeviewtempdir)==1){
      strcpy(renderfile_dir,smokeviewtempdir);
    }
    else{
      printf("unable to output render file\n");
      return;
    }
  }

  // filename suffix

  if(use_scriptfile==0||
    (current_script_command!=NULL&&
    (current_script_command->command==SCRIPT_RENDERALL||
     current_script_command->command==SCRIPT_VOLSMOKERENDERALL))){
    int image_num;
    char suffix[20];

    strcpy(renderfile_suffix,"_");
    if(RenderTime==0){
      image_num=seqnum;
    }
    else{
      image_num=itimes/RenderSkip;
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
    switch (view_mode){
    case VIEW_CENTER:
      if(RenderTime==0)seqnum++;
      break;
    case VIEW_LEFT:
      strcat(suffix,"_L");
      break;
    case VIEW_RIGHT:
      if(showstereo==0||showstereo==1){
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

  // filename extension

  switch (renderfiletype){
  case 0:
    renderfile_ext=ext_png;
    break;
#ifdef pp_JPEG
  case 1:
    renderfile_ext=ext_jpg;
    break;
#endif
#ifdef pp_GDGIF
  case 2:
    renderfile_ext=ext_gif;
    break;
#endif
  default:
    renderfiletype=2;
    renderfile_ext=ext_png;
    break;
  }

  // form full filename from parts

  strcpy(renderfile_full,"");
  if(strlen(renderfile_dir)>0){
    strcat(renderfile_full,renderfile_dir);
    if(renderfile_dir[strlen(renderfile_dir)-1]!=dirseparator[0]){
      strcat(renderfile_full,dirseparator);
    }
  }
  strcat(renderfile_full,renderfile_name);
  if(strlen(renderfile_suffix)>0)strcat(renderfile_full,renderfile_suffix);
  strcat(renderfile_full,renderfile_ext);

  // render image

  SVimage2file(renderfile_full,renderfiletype,screenWidth,screenHeight);
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

int mergescreenbuffers(GLubyte *screenbuffers[4]){

  char renderfile[1024];
  char renderfile2[1024];
  char *ext;
  FILE *RENDERfile;
  GLubyte *p;
  gdImagePtr RENDERimage;
  unsigned int r, g, b;
  int i,j,rgb_local;

  switch (renderfiletype){
  case PNG:
    ext=ext_png;
    break;
#ifdef pp_JPEG
  case JPEG:
    ext=ext_jpg;
    break;
#endif
#ifdef pp_GDGIF
  case GIF:
    ext=ext_gif;
    break;
#endif
  default:
    renderfiletype=0;
    ext=ext_png;
    break;
  }

  if(scriptoutstream!=NULL&&current_script_command!=NULL&&current_script_command->cval!=NULL){
    strcpy(renderfile2,current_script_command->cval);
  }
  else{
    strcpy(renderfile2,fdsprefix);
  }
  if(RenderTime==1){
      sprintf(renderfile,"%s_%04i",renderfile2,itimes/RenderSkip);
  }
  if(RenderTime==0){
      sprintf(renderfile,"%s_s%04i",renderfile2,seqnum);
      seqnum++;
  }
  strcat(renderfile,ext);

  // if there is a tempdir see if we need to use it

  if(smokeviewtempdir!=NULL){
    FILE *stream;

    stream=fopen(renderfile,"wb");
    if(stream==NULL){
      strcpy(renderfile2,smokeviewtempdir);
      strcat(renderfile2,renderfile);
      stream=fopen(renderfile2,"wb");
      if(stream!=NULL)strcpy(renderfile,renderfile2);
    }
    if(stream!=NULL)fclose(stream);
  }
  printf("Rendering to: %s .",renderfile);
  RENDERfile = fopen(renderfile, "wb");
  if (RENDERfile == NULL) {
    {
      char message[256];

      sprintf(message,"unable to write to %s",renderfile);
      warning_message(message);
    }
    return 1;
  }
  RENDERimage = gdImageCreateTrueColor(2*screenWidth,2*screenHeight);

  p=screenbuffers[0];
  for (i = 2*screenHeight-1 ; i>=screenHeight; i--) {
    for(j=0;j<screenWidth;j++){
      r=*p++; g=*p++; b=*p++;
      rgb_local = (r<<16)|(g<<8)|b;
      gdImageSetPixel(RENDERimage,j,i,rgb_local);

    }
  }
  p=screenbuffers[1];
  for (i = 2*screenHeight-1 ; i>=screenHeight; i--) {
    for(j=screenWidth;j<2*screenWidth;j++){
      r=*p++; g=*p++; b=*p++;
      rgb_local = (r<<16)|(g<<8)|b;
      gdImageSetPixel(RENDERimage,j,i,rgb_local);

    }
  }
  p=screenbuffers[2];
  for (i = screenHeight-1 ; i>=0; i--) {
    for(j=0;j<screenWidth;j++){
      r=*p++; g=*p++; b=*p++;
      rgb_local = (r<<16)|(g<<8)|b;
      gdImageSetPixel(RENDERimage,j,i,rgb_local);

    }
  }
  p=screenbuffers[3];
  for (i = screenHeight-1 ; i>=0; i--) {
    for(j=screenWidth;j<2*screenWidth;j++){
      r=*p++; g=*p++; b=*p++;
      rgb_local = (r<<16)|(g<<8)|b;
      gdImageSetPixel(RENDERimage,j,i,rgb_local);

    }
  }

  /* output the image */

  switch (renderfiletype){
  case PNG:
    gdImagePng(RENDERimage,RENDERfile);
    break;
#ifdef pp_JPEG
  case JPEG:
    gdImageJpeg(RENDERimage,RENDERfile,-1);
    break;
#endif
#ifdef pp_GDGIF
  case GIF:
    gdImageGif(RENDERimage,RENDERfile);
    break;
#endif
  default:
    ASSERT(0);
    break;
  }
  fclose(RENDERfile);

  /* free up memory used by both OpenGL and GIF images */

  gdImageDestroy(RENDERimage);
  printf(" Completed\n");
  return 0;
}



/* ------------------ SVimage2file ------------------------ */

int SVimage2file(char *RENDERfilename, int rendertype, int width, int height){

  FILE *RENDERfile;
  GLubyte *OpenGLimage, *p;
  gdImagePtr RENDERimage;
  unsigned int r, g, b;
  int i,j,rgb_local;
  int x=0, y=0;

  RENDERfile = fopen(RENDERfilename, "wb");
  if (RENDERfile == NULL) {
    char message[1024];

    strcpy(message,_("unable to write to "));
    strcat(message,RENDERfilename);
    warning_message(message);
    return 1;
  }
  NewMemory((void **)&OpenGLimage,width * height * sizeof(GLubyte) * 3);
  if(OpenGLimage == NULL){
    char message[1024];
    
    strcpy(message,_("error allocating render image: "));
    strcat(message,RENDERfilename);
    error_message(message);
    pauseSV();
    exit(1);
  }
  printf("Rendering to: %s .",RENDERfilename);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);

  /* get the image from the OpenGL frame buffer */

  glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE, OpenGLimage);

  /* copy the image from OpenGL memory to GIF memory */

  p = OpenGLimage;

  RENDERimage = gdImageCreateTrueColor(width,height);

  for (i = height-1 ; i>=0; i--) {
    for(j=0;j<width;j++){
      r=*p++; g=*p++; b=*p++;
      rgb_local = (r<<16)|(g<<8)|b;
      gdImageSetPixel(RENDERimage,j,i,rgb_local);
    }
  }

  if(test_smokesensors==1&&active_smokesensors==1&&show_smokesensors!=0){
    int idev;

    for(idev=0;idev<ndeviceinfo;idev++){
      device *devicei;
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

  switch (rendertype){
  case PNG:
    gdImagePng(RENDERimage,RENDERfile);
    break;
#ifdef pp_JPEG
  case JPEG:
    gdImageJpeg(RENDERimage,RENDERfile,-1);
    break;
#endif
#ifdef pp_GDGIF
  case GIF:
    gdImageGif(RENDERimage,RENDERfile);
    break;
#endif
  default:
    ASSERT(0);
    break;
  }

  /* free up memory used by both OpenGL and GIF images */

  fclose(RENDERfile);

  gdImageDestroy(RENDERimage);
  FREEMEMORY(OpenGLimage);
  printf(" Completed.\n");
  return 0;
}

/* ------------------ readpicture ------------------------ */

unsigned char *readpicture(char *filename, int *width, int *height){
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
      printf("Texture file:%s unavailable\n",filename);
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
        FREEMEMORY(filebuffer);
        printf("Texture file:%s unavailable\n",filebuffer);
        return NULL;
      }
      else{
        fclose(stream);
      }
    }
  }

  
  printf("Loading texture:%s ",filebuffer);
  ext = filebuffer + strlen(filebuffer) - 4;
#ifdef pp_JPEG
  if(strncmp(ext,".jpg",4)==0||strncmp(ext,".JPG",4)==0){
    returncode = readjpeg(filebuffer,width,height,pixel_skip);
  }
  else if(strncmp(ext,".png",4)==0||strncmp(ext,".PNG",4)==0){
    returncode = readpng(filebuffer,width,height);
  }
#else
  if(strncmp(ext,".png",4)==0||strncmp(ext,".PNG",4)==0){
    returncode = readpng(filebuffer,width,height);
  }
#endif
  else{
    if(allocated==1){
      FREEMEMORY(filebuffer);
    }
    return NULL;
  }
  if(allocated==1){
    FREEMEMORY(filebuffer);
  }
  if(returncode!=NULL){
    printf(" completed\n");
  }
  else{
    printf(" failed\n");
  }
  return returncode;

}

#ifdef pp_JPEG
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
  if (file == NULL)return NULL;
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
#endif

/* ------------------ readpng ------------------------ */

unsigned char *readpng(const char *filename,int *width, int *height){

  FILE *file;
  gdImagePtr image;
  unsigned char *dataptr,*dptr;
  int i,j;
  unsigned int intrgb;

  file = fopen(filename, "rb");
  if (file == NULL)return NULL;
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
