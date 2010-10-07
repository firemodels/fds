// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "smokeviewvars.h"

#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "MALLOC.h"
#include "gd.h"

#define PNG 0
#ifdef pp_JPEG
#define JPEG 1
#endif
#ifdef pp_GDGIF
#define GIF 2
#endif

// svn revision character string
char renderfile_revision[]="$Revision$";


/* ------------------ Render ------------------------ */

void Render(int view_mode){
  if(RenderOnceNow==0&&RenderGif !=0
    &&render_double==0
    ){
    if(plotstate==DYNAMIC_PLOTS && ntimes > 0){
     if(itime>=0&&itime<ntimes&&
       ((render_frame[itime] == 0&&showstereo==0)||(render_frame[itime]<2&&showstereo!=0))
       ){
       render_frame[itime]++;
       RenderFrame(view_mode);
     }
     else{
       ASSERT(RenderSkip>0);
       RenderState(0);
       RenderSkip=1;
     }
    }
    if(touring == 1 && ntimes == 0){
      if(rendertourcount % RenderSkip == 0)RenderFrame(view_mode);
      rendertourcount++;
      if(ntimes>0)tourangle += (float)(2.0*PI/((float)ntimes/(float)RenderSkip));
      if(ntimes==0)tourangle += (float)(2.0*PI/((float)maxtourframes/(float)RenderSkip));
      if(tourangle>2.0*PI){
        RenderState(0);
        RenderSkip=1;
        tourangle=0.0;
      }
    }
  }

  if(render_double==0)sniffErrors("after render");

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

//void pauseSV(void);


  /* ------------------ can_write_to_dir ------------------------ */

int can_write_to_dir(char *dir){
  char full_name[1024];
  char temp_name[1024];
  FILE *stream;

  if(dir==NULL)return 0;

  strcpy(temp_name,fdsprefix);
  strcat(temp_name,".write_test");

  if(strcmp(dir,".")==0||strlen(dir)==0){
    strcpy(full_name,temp_name);
  }
  else{
    strcpy(full_name,dir);
    strcat(full_name,dirseparator);
    strcat(full_name,temp_name);
  }
  stream=fopen(full_name,"wb");
  if(stream==NULL)return 0;
  fclose(stream);
  remove(full_name);

  return 1;
}

  /* ------------------ RenderFrame ------------------------ */

void RenderFrame(int view_mode){
  char renderfile[1024],renderfile2[1024];
  FILE *stream;
  char *ext;
  char *renderfile_prefix;
  int use_script_filename=0;
  char renderfile_name[1024], renderfile_dir[1024], renderfile_full[1024], renderfile_suffix[1024], *renderfile_ext;
  char *temp_name;
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
      (current_script_command->command==SCRIPT_RENDERONCE||current_script_command->command==SCRIPT_RENDERALL)&&
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

  if(use_scriptfile==0||(current_script_command!=NULL&&current_script_command->command==SCRIPT_RENDERALL)){
    int image_num;
    char suffix[20];

    if(RenderTime==0){
      image_num=seqnum;
      strcpy(renderfile_suffix,"_s");
    }
    else{
      image_num=itime/RenderSkip;
      strcpy(renderfile_suffix,"_");
    }
    switch (view_mode){
    case VIEW_CENTER:
      sprintf(suffix,"%04i",image_num);
      if(RenderTime==0)seqnum++;
      break;
    case VIEW_LEFT:
      sprintf(suffix,"%04i_L",image_num);
      break;
    case VIEW_RIGHT:
      if(showstereo!=0&&showstereo!=1){
        sprintf(suffix,"%04i",image_num);
      }
      else{
        sprintf(suffix,"%04i_R",image_num);
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

  GLubyte *OpenGLimage;

  int x=0, y=0;

  OpenGLimage = (GLubyte *) malloc(screenWidth * screenHeight * sizeof(GLubyte) * 3);

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
  int i,j,rgb;

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
      sprintf(renderfile,"%s_%04i",renderfile2,itime/RenderSkip);
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
#ifdef pp_MESSAGE
    {
      char message[256];

      sprintf(message,"*** warning: unable to write to %s\n",renderfile);
      warning_message(message);
    }
#else
    printf("*** warning: unable to write to %s\n",renderfile);
#endif
    return 1;
  }
  RENDERimage = gdImageCreateTrueColor(2*screenWidth,2*screenHeight);

  p=screenbuffers[0];
  for (i = 2*screenHeight-1 ; i>=screenHeight; i--) {
    for(j=0;j<screenWidth;j++){
      r=*p++; g=*p++; b=*p++;
      rgb = (r<<16)|(g<<8)|b;
      gdImageSetPixel(RENDERimage,j,i,rgb);

    }
  }
  p=screenbuffers[1];
  for (i = 2*screenHeight-1 ; i>=screenHeight; i--) {
    for(j=screenWidth;j<2*screenWidth;j++){
      r=*p++; g=*p++; b=*p++;
      rgb = (r<<16)|(g<<8)|b;
      gdImageSetPixel(RENDERimage,j,i,rgb);

    }
  }
  p=screenbuffers[2];
  for (i = screenHeight-1 ; i>=0; i--) {
    for(j=0;j<screenWidth;j++){
      r=*p++; g=*p++; b=*p++;
      rgb = (r<<16)|(g<<8)|b;
      gdImageSetPixel(RENDERimage,j,i,rgb);

    }
  }
  p=screenbuffers[3];
  for (i = screenHeight-1 ; i>=0; i--) {
    for(j=screenWidth;j<2*screenWidth;j++){
      r=*p++; g=*p++; b=*p++;
      rgb = (r<<16)|(g<<8)|b;
      gdImageSetPixel(RENDERimage,j,i,rgb);

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
  int i,j,rgb;
  int x=0, y=0;

  RENDERfile = fopen(RENDERfilename, "wb");
  if (RENDERfile == NULL) {
#ifdef pp_MESSAGE
    {
      char message[256];

      sprintf(message,"*** warning: unable to write to %s\n",RENDERfilename);
      warning_message(message);
    }
#else
    printf("*** warning:  unable to write to %s\n",RENDERfilename);
#endif
    return 1;
  }
  OpenGLimage = (GLubyte *) malloc(width * height * sizeof(GLubyte) * 3);
  if(OpenGLimage == NULL){
    printf("error allocating render image:%s\n",RENDERfilename);                                                                
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
      rgb = (r<<16)|(g<<8)|b;
      gdImageSetPixel(RENDERimage,j,i,rgb);
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
  free(OpenGLimage);
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
  else if(strncmp(ext,".rgb",4)==0||strncmp(ext,".RGB",4)==0){
    returncode = readrgb(filebuffer,width,height);
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

unsigned char *readjpeg(const char *filename,int *width, int *height, int skip){

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
  if(skip<0)skip=0;
  jump = skip + 1;
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

/* ------------------ bwtorgba ------------------------ */

void bwtorgba(const unsigned char *b,unsigned char *l,int n) {
    while(n--) {
        l[0] = *b;
        l[1] = *b;
        l[2] = *b;
        l[3] = 0xff;
        l += 4; b++;
    }
}

/* ------------------ rgbtorgba ------------------------ */

void rgbtorgba(const unsigned char *r,
               const unsigned char *g,
               const unsigned char *b,
               unsigned char *l,int n) {
    while(n--) {
        l[0] = r[0];
        l[1] = g[0];
        l[2] = b[0];
        l[3] = 0xff;
        l += 4; r++; g++; b++;
    }
}

/* ------------------ rgbatorgba ------------------------ */

void rgbatorgba(const unsigned char *r,
                const unsigned char *g,
                const unsigned char *b,
                const unsigned char *a,
                unsigned char *l,int n) {
    while(n--) {
        l[0] = r[0];
        l[1] = g[0];
        l[2] = b[0];
        l[3] = a[0];
        l += 4; r++; g++; b++; a++;
    }
}

typedef struct _ImageRec {
    unsigned short imagic;
    unsigned short type;
    unsigned short dim;
    unsigned short xsize, ysize, zsize;
    unsigned int min, max;
    unsigned int wasteBytes;
    char name[80];
    unsigned long colorMap;
    FILE *file;
    unsigned char *tmp, *tmpR, *tmpG, *tmpB;
    unsigned long rleEnd;
    unsigned int *rowStart;
    int *rowSize;
} ImageRec;

/* ------------------ ConvertShort ------------------------ */

static void ConvertShort(unsigned short *array, unsigned int length) {
    unsigned short b1, b2;
    unsigned char *ptr;

    ptr = (unsigned char *)array;
    while (length--) {
        b1 = *ptr++;
        b2 = *ptr++;
        *array++ = (b1 << 8) | (b2);
    }
}

/* ------------------ ConvertUint ------------------------ */

static void ConvertUint(unsigned *array, unsigned int length) {
    unsigned int b1, b2, b3, b4;
    unsigned char *ptr;

    ptr = (unsigned char *)array;
    while (length--) {
        b1 = *ptr++;
        b2 = *ptr++;
        b3 = *ptr++;
        b4 = *ptr++;
        *array++ = (b1 << 24) | (b2 << 16) | (b3 << 8) | (b4);
    }
}

/* ------------------ ImageOpen ------------------------ */

static ImageRec *ImageOpen(const char *fileName)
{
    union {
        int testWord;
        char testByte[4];
    } endianTest;
    ImageRec *image;
    int swapFlag;
    int x;

    endianTest.testWord = 1;
    if (endianTest.testByte[0] == 1) {
        swapFlag = 1;
    } else {
        swapFlag = 0;
    }

    image = (ImageRec *)malloc(sizeof(ImageRec));
    if (image == NULL) {
        fprintf(stderr, "Out of memory!\n");
        pauseSV();
        exit(1);
    }
    if ((image->file = fopen(fileName, "rb")) == NULL) {
        perror(fileName);
        pauseSV();
        exit(1);
    }

    fread(image, 1, 12, image->file);

    if (swapFlag) {
        ConvertShort(&image->imagic, 6);
    }

    image->tmp = (unsigned char *)malloc(image->xsize*256);
    image->tmpR = (unsigned char *)malloc(image->xsize*256);
    image->tmpG = (unsigned char *)malloc(image->xsize*256);
    image->tmpB = (unsigned char *)malloc(image->xsize*256);
    if (image->tmp == NULL || image->tmpR == NULL || image->tmpG == NULL ||
        image->tmpB == NULL) {
        fprintf(stderr, "Out of memory!\n");
        pauseSV();
        exit(1);
    }

    if ((image->type & 0xFF00) == 0x0100) {
        x = image->ysize * image->zsize * (int) sizeof(unsigned);
        image->rowStart = (unsigned *)malloc((unsigned int)x);
        image->rowSize = (int *)malloc((unsigned int)x);
        if (image->rowStart == NULL || image->rowSize == NULL) {
            fprintf(stderr, "Out of memory!\n");
            pauseSV();
            exit(1);
        }
        image->rleEnd = (unsigned long)(512 + (2 * x));
        fseek(image->file, 512, SEEK_SET);
        fread(image->rowStart, 1, (unsigned int)x, image->file);
        fread(image->rowSize, 1, (unsigned int)x, image->file);
        if (swapFlag) {
            ConvertUint(image->rowStart, (unsigned int)(x/(int) sizeof(unsigned)));
            ConvertUint((unsigned *)image->rowSize, (unsigned int)(x/(int) sizeof(int)));
        }
    }
    return image;
}

/* ------------------ ImageClose ------------------------ */

static void ImageClose(ImageRec *image) {
    fclose(image->file);
    free(image->tmp);
    free(image->tmpR);
    free(image->tmpG);
    free(image->tmpB);
    free(image);
}

/* ------------------ ImageGetRow ------------------------ */

static void ImageGetRow(const ImageRec *image, unsigned char *buf, int y, int z) {
    unsigned char *iPtr, *oPtr, pixel;
    int count;

    if ((image->type & 0xFF00) == 0x0100) {
        fseek(image->file, (long) image->rowStart[y+z*image->ysize], SEEK_SET);
        fread(image->tmp, 1, (unsigned int)image->rowSize[y+z*image->ysize],
              image->file);

        iPtr = image->tmp;
        oPtr = buf;
        for (;;) {
            pixel = *iPtr++;
            count = (int)(pixel & 0x7F);
            if (!count) {
                return;
            }
            if (pixel & 0x80) {
                while (count--) {
                    *oPtr++ = *iPtr++;
                }
            } else {
                pixel = *iPtr++;
                while (count--) {
                    *oPtr++ = pixel;
                }
            }
        }
    } else {
        fseek(image->file, 512+(y*image->xsize)+(z*image->xsize*image->ysize),
              SEEK_SET);
        fread(buf, 1, image->xsize, image->file);
    }
}

/* ------------------ imagergb ------------------------ */

unsigned char *readrgb(const char *name, int *width, int *height) {
    unsigned char *base, *lptr;
    unsigned char *rbuf, *gbuf, *bbuf, *abuf;
    int comp;
    int *components;
    ImageRec *image;
    int y;

    components = &comp;
    image = ImageOpen(name);
    
    if(!image)
        return NULL;
    (*width)=image->xsize;
    (*height)=image->ysize;
    (*components)=image->zsize;
    base=NULL;
    if(NewMemory((void **)&base,4*image->xsize*image->ysize*sizeof(unsigned char))==0){
      return NULL;
    }
    rbuf = (unsigned char *)malloc(image->xsize*sizeof(unsigned char));
    if(rbuf==NULL)return NULL;
    gbuf = (unsigned char *)malloc(image->xsize*sizeof(unsigned char));
    if(gbuf==NULL){
      free(rbuf);
      return NULL;
    }
    bbuf = (unsigned char *)malloc(image->xsize*sizeof(unsigned char));
    if(bbuf==NULL){
      free(rbuf);free(gbuf);
      return NULL;
    }
    abuf = (unsigned char *)malloc(image->xsize*sizeof(unsigned char));
    if(abuf==NULL){
      free(rbuf);free(gbuf);free(bbuf);
      return NULL;
    }
    lptr = base;
    for(y=0; y<image->ysize; y++) {
        if(image->zsize>=4) {
            ImageGetRow(image,rbuf,y,0);
            ImageGetRow(image,gbuf,y,1);
            ImageGetRow(image,bbuf,y,2);
            ImageGetRow(image,abuf,y,3);
            rgbatorgba(rbuf,gbuf,bbuf,abuf,(unsigned char *)lptr,image->xsize);
            lptr += 4*image->xsize;
        } else if(image->zsize==3) {
            ImageGetRow(image,rbuf,y,0);
            ImageGetRow(image,gbuf,y,1);
            ImageGetRow(image,bbuf,y,2);
            rgbtorgba(rbuf,gbuf,bbuf,(unsigned char *)lptr,image->xsize);
            lptr += 4*image->xsize;
        } else {
            ImageGetRow(image,rbuf,y,0);
            bwtorgba(rbuf,(unsigned char *)lptr,image->xsize);
            lptr += image->xsize;
        }
    }
    free(abuf);
    ImageClose(image);
    free(rbuf);
    free(gbuf);
    free(bbuf);
    return (unsigned char *) base;
}
