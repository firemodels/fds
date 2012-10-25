// $Date$ 
// $Revision$
// $Author$

// svn revision character string
char  viewports_revision[]="$Revision$";

#include "options.h"
#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "string_util.h"
#include "smokeviewvars.h"
#include "IOvolsmoke.h"

#define CONV(p,pl,pr,pxl,pxr) ( (pxl) + ((pxr)-(pxl))*((p)-(pl))/((pr)-(pl)) )

 /* ------------------------ SUB_portortho ------------------------- */
 
int SUB_portortho(int quad, 
                   GLint port_left, GLint port_down, GLsizei port_width, GLsizei port_height,
                   GLdouble portx_left, GLdouble portx_right, GLdouble portx_down, GLdouble portx_top,
                   GLint screen_left, GLint screen_down, GLsizei screen_width, GLsizei screen_height
                   ){
  
  GLint subport_left, subport_right, subport_down, subport_top;
  GLdouble subportx_left, subportx_right, subportx_down, subportx_top;
  GLsizei subport_width, subport_height;
  GLint subwindow_left, subwindow_right, subwindow_down, subwindow_top;
  GLint port_right, port_top;

  int irow, icol;

  switch (quad){
  case 0:            
    port_pixel_width = port_width;
    port_pixel_height = port_height;
    port_unit_width = portx_right - portx_left;
    port_unit_height = portx_top - portx_down;
    glViewport(port_left,port_down,port_width,port_height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(portx_left,portx_right,portx_down,portx_top);
    return 1;
  case 1:
    icol = screen_left/screen_width;
    irow = screen_down/screen_height;

    subwindow_left = icol*screen_width;
    subwindow_right = subwindow_left + screen_width;
    subwindow_down = irow*screen_height;
    subwindow_top = subwindow_down + screen_height;

    port_right = port_left + port_width;
    port_top = port_down + port_height;

    subport_left =  MAX( nrender_rows*port_left,subwindow_left);
    subport_right = MIN(nrender_rows*port_right,subwindow_right);
    subport_down =  MAX( nrender_rows*port_down,subwindow_down);
    subport_top =   MIN(  nrender_rows*port_top,subwindow_top);
    if(subport_left>=subport_right||subport_down>=subport_top)return 0;

#define CONV(p,pl,pr,pxl,pxr) ( (pxl) + ((pxr)-(pxl))*((p)-(pl))/((pr)-(pl)) )

    subportx_left = CONV(subport_left,nrender_rows*port_left,nrender_rows*port_right,portx_left,portx_right);
    subportx_right = CONV(subport_right,nrender_rows*port_left,nrender_rows*port_right,portx_left,portx_right);
    subportx_down = CONV(subport_down,nrender_rows*port_down,nrender_rows*port_top,portx_down,portx_top);
    subportx_top = CONV(subport_top,nrender_rows*port_down,nrender_rows*port_top,portx_down,portx_top);

    subport_left -= icol*screen_width;
    subport_right -= icol*screen_width;
    subport_down -= irow*screen_height;
    subport_top -= irow*screen_height;
    subport_width = subport_right - subport_left;
    subport_height = subport_top - subport_down;

    port_pixel_width = subport_width;
    port_pixel_height = subport_height;
    port_unit_width = subportx_right - subportx_left;
    port_unit_height = subportx_top - subportx_down;

    glViewport(subport_left,subport_down,subport_width,subport_height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(subportx_left,subportx_right,subportx_down,subportx_top);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  return 1;
}

/* ------------------------ SUB_portfrustum ------------------------- */
 
int SUB_portfrustum(int quad, 
                   GLint port_left, GLint port_down, GLsizei port_width, GLsizei port_height,
                   GLdouble portx_left, GLdouble portx_right, 
                   GLdouble portx_down, GLdouble portx_top,
                   GLdouble portx_near, GLdouble portx_far,
                   GLint screen_left, GLint screen_down, GLsizei screen_width, GLsizei screen_height
                   ){
  GLint subport_left, subport_right, subport_down, subport_top;
  GLdouble subportx_left, subportx_right, subportx_down, subportx_top;
  GLsizei subport_width, subport_height;
  GLint subwindow_left, subwindow_right, subwindow_down, subwindow_top;
  GLint port_right, port_top;

  int irow, icol;

  switch (quad){
  case 0:
    port_pixel_width = port_width;
    port_pixel_height = port_height;
    port_unit_width = portx_right - portx_left;
    port_unit_height = portx_top - portx_down;
    glViewport(port_left,port_down,port_width,port_height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if(camera_current->projection_type==0){
      glFrustum(
        (double)portx_left,(double)portx_right,
        (double)portx_down,(double)portx_top,
        (double)portx_near,(double)portx_far);
    }
    else{
      glOrtho(
        (double)portx_left,(double)portx_right,
        (double)portx_down,(double)portx_top,
        (double)portx_near,(double)portx_far);
    }
    return 1;
  case 1:
    icol = screen_left/screen_width;
    irow = screen_down/screen_height;

    subwindow_left = icol*screen_width;
    subwindow_right = subwindow_left + screen_width;
    subwindow_down = irow*screen_height;
    subwindow_top = subwindow_down + screen_height;

    port_right = port_left + port_width;
    port_top = port_down + port_height;

    subport_left =  MAX( nrender_rows*port_left,subwindow_left);
    subport_right = MIN(nrender_rows*port_right,subwindow_right);
    subport_down =  MAX( nrender_rows*port_down,subwindow_down);
    subport_top =   MIN(  nrender_rows*port_top,subwindow_top);
    if(subport_left>=subport_right||subport_down>=subport_top)return 0;

    subportx_left = CONV(subport_left,nrender_rows*port_left,nrender_rows*port_right,portx_left,portx_right);
    subportx_right = CONV(subport_right,nrender_rows*port_left,nrender_rows*port_right,portx_left,portx_right);
    subportx_down = CONV(subport_down,nrender_rows*port_down,nrender_rows*port_top,portx_down,portx_top);
    subportx_top = CONV(subport_top,nrender_rows*port_down,nrender_rows*port_top,portx_down,portx_top);

    subport_left -= icol*screen_width;
    subport_right -= icol*screen_width;
    subport_down -= irow*screen_height;
    subport_top -= irow*screen_height;
    subport_width = subport_right - subport_left;
    subport_height = subport_top - subport_down;

    port_pixel_width = subport_width;
    port_pixel_height = subport_height;
    port_unit_width = subportx_right - subportx_left;
    port_unit_height = subportx_top - subportx_down;

    glViewport(subport_left,subport_down,subport_width,subport_height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if(camera_current->projection_type==0){
      glFrustum(
        (double)subportx_left,(double)subportx_right,
        (double)subportx_down,(double)subportx_top,
        (double)portx_near,(double)portx_far);
    }
    else{
      glOrtho(
        (double)subportx_left,(double)subportx_right,
        (double)subportx_down,(double)subportx_top,
        (double)portx_near,(double)portx_far);
    }
    return 1;
  default:
    ASSERT(FFALSE);
    break;
  }
  return 1;
}


 /* ------------------------ CLIP_viewport ------------------------- */

void CLIP_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height){
  GLdouble x_left, x_right, x_down, x_top;
  float c_left, c_right, c_top, c_bottom;

  x_left=0.0;
  x_right=screenWidth;
  x_down=0.0;
  x_top=screenHeight;

  if(SUB_portortho(quad,0,0,screenWidth,screenHeight,
        x_left, x_right, x_down, x_top,
        s_left, s_down, s_width, s_height)==0){
          return;
   }

   c_left = render_clip_left-3;
   c_right = screenWidth + 3 - render_clip_right;
   c_bottom = render_clip_bottom -3;
   c_top = screenHeight + 3 - render_clip_top;

   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   glLineWidth(3.0);
   glBegin(GL_LINES);

   if(c_left>0){
     glVertex2f(c_left,c_bottom);
     glVertex2f(c_left,c_top);
   }

   if(c_right<screenWidth){
     glVertex2f(c_right,c_bottom);
     glVertex2f(c_right,c_top);
   }

   if(c_top<screenHeight){
     glVertex2f(c_left,c_top);
     glVertex2f(c_right,c_top);
   }

   if(c_bottom>0){
     glVertex2f(c_left,c_bottom);
     glVertex2f(c_right,c_bottom);
   }
   glEnd();
     

}

 /* ------------------------ BLOCK viewport ------------------------- */

void BLOCK_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height){
  float mesh_left;
  char slicelabel[255];
  float mesh_bot;
  int portview_left;
  float val_right,val_top;
  mesh *mesh_xyz=NULL;
  float xyz[3];

  mesh_left=0.9;
  if(fontindex==LARGE_FONT)mesh_left=0.7;
 
  portview_left=screenWidth-dwinW-fontWoffset-titlesafe_offset;
  if(screenWidth<screenHeight){
    val_right=1.0;
    val_top=ratio;
    if(SUB_portortho(quad,
      portview_left,
      titlesafe_offset,
      dwinW, 
      dwinH-fontHoffset,
      0.,1.0,0.,(double)ratio,
      s_left, s_down, s_width, s_height)==0)return;
  }
  else{
    val_right=ratio;
    val_top=1.0;
    if(SUB_portortho(quad,
      portview_left,
      titlesafe_offset,
      dwinW, 
      dwinH-fontHoffset,
      0.,(double)ratio,0.0,1.0,
      s_left, s_down, s_width, s_height)==0)return;
  }

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  if((showplot3d==1||visGrid!=noGridnoProbe)&&(visx_all==1||visy_all||visz_all)||visGrid==GridProbe||visGrid==noGridProbe){
    xyz[0]=DENORMALIZE_X(plotx_all[iplotx_all]);
    xyz[1]=DENORMALIZE_Y(ploty_all[iploty_all]);
    xyz[2]=DENORMALIZE_Z(plotz_all[iplotz_all]);
    mesh_xyz=getmesh(xyz);
  }
  if(visBlocklabel==1&&nmeshes>1){
    int labellength;
    char meshlabel[255];

    if(mesh_xyz==NULL){
      sprintf(meshlabel,"mesh: %i",highlight_mesh+1);
      mesh_xyz = meshinfo + highlight_mesh;
    }
    else{
      int imesh;

      imesh = mesh_xyz-meshinfo+1;
      sprintf(meshlabel,"mesh: %i",imesh);
    }
    labellength=glutBitmapLength(large_font, (const unsigned char *)meshlabel);
    mesh_left=val_right-val_right*labellength/(float)dwinW;
    mesh_bot=val_top-val_top*large_font_height/(float)(dwinH-fontHoffset);
    outputText(mesh_left,mesh_bot, meshlabel);
  }
  if(((showplot3d==1||visGrid!=noGridnoProbe)&&visx_all==1)||visGrid==noGridProbe||visGrid==GridProbe){
    float plotval;
    int iplotval;
    char buff_label[128];


    iplotval=mesh_xyz->iplotx_all[iplotx_all];
    plotval=xyz[0];
    if(plotval>0.0){
      plotval=(int)(plotval*100+0.5);
    }
    else{
      plotval=(int)(plotval*100-0.5);
    }
    plotval/=100;
          
    sprintf(buff_label,"%f",plotval);
    trimzeros(buff_label);
    strcat(buff_label," m");
    if(cursorPlot3D==1){
      sprintf(slicelabel,"*x: %i, ",iplotval);
    }
    else{
      sprintf(slicelabel,"x: %i, ",iplotval);
    }
    strcat(slicelabel,buff_label);
    if(visgridloc==1)outputText(mesh_left-0.5,0.6f, slicelabel);
  }
  if(((showplot3d==1||visGrid!=noGridnoProbe)&&visy_all==1)||visGrid==GridProbe||visGrid==noGridProbe){
    float plotval;
    int iplotval;
    char buff_label[128];

    iplotval=mesh_xyz->iploty_all[iploty_all];
    plotval=xyz[1];
    if(plotval>0.0){
      plotval=(int)(plotval*100+0.5);
    }
    else{
      plotval=(int)(plotval*100-0.5);
    }
    plotval/=100;
          
    sprintf(buff_label,"%f",plotval);
    trimzeros(buff_label);
    strcat(buff_label," m");
    if(cursorPlot3D==1){
      sprintf(slicelabel,"*y: %i, ",iplotval);
    }
    else{
      sprintf(slicelabel,"y: %i, ",iplotval);
    }
    strcat(slicelabel,buff_label);
    if(visgridloc==1)outputText(mesh_left-0.5,0.35f, slicelabel);
  }
  if(((showplot3d==1||visGrid!=noGridnoProbe)&&visz_all==1)||visGrid==GridProbe||visGrid==noGridProbe){
    float plotval;
    int iplotval;
    char buff_label[128];

    iplotval=mesh_xyz->iplotz_all[iplotz_all];
    plotval=xyz[2];
    if(plotval>0.0){
      plotval=(int)(plotval*100+0.5);
    }
    else{
      plotval=(int)(plotval*100-0.5);
    }
    plotval/=100;
          
    sprintf(buff_label,"%f",plotval);
    trimzeros(buff_label);
    strcat(buff_label," m");
    if(cursorPlot3D==1){
      sprintf(slicelabel,"*z: %i, ",iplotval);
    }
    else{
      sprintf(slicelabel,"z: %i, ",iplotval);
    }
    strcat(slicelabel,buff_label);
    if(visgridloc==1)outputText(mesh_left-.5,0.1f, slicelabel);
  }
}

#ifdef pp_MEMDEBUG
/* ------------------ getMemusage ------------------------ */

void getMemusage(MMsize totalmemory,char *MEMlabel){
  int size;
  float rsize;

  if(totalmemory<1000000000){
    size = totalmemory/1000000;
    sprintf(MEMlabel,"%i MB",size);
  }
  else{
    rsize = totalmemory/1000000000.0;
    sprintf(MEMlabel,"%4.2f GB",rsize);
  }

}
#endif

/* ------------------------ TIME BAR Viewport ------------------------- */

void TIMEBAR_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height){
  int timebarheight;
#ifdef pp_memstatus
  unsigned int availmemory;
  char percen[]="%";
#endif


  if(
    (visTimeLabels==1&&showtime==1)||
    (showtime==1&&
      (visFramerate==1||benchmark==1||
       (vis_slice_average==1&&show_slice_average&&slice_average_flag==1)
      )
      ||(hrrpuv_loaded==1&&show_hrrcutoff==1&&current_mesh!=NULL)
    )
#ifdef pp_memstatus
    ||visAvailmemory==1
#endif
    )
  {
    timebarheight=(int)(0.75*dwinH);
    if(screenWidth<screenHeight){
      if(SUB_portortho(quad,
        fontWoffset+titlesafe_offset, fontHoffset+titlesafe_offset, screenWidth-dwinWW-2*fontWoffset-2*titlesafe_offset, timebarheight,
        0.0,1.0,0.0,0.75*(double)ratio,
        s_left, s_down, s_width, s_height)==0){
          return;
      }

      xtemp=1.0;
    }
    else{
      if(SUB_portortho(quad,
        fontWoffset+titlesafe_offset, fontHoffset+titlesafe_offset, screenWidth-dwinWW-2*fontWoffset-2*titlesafe_offset, timebarheight,
        0.0,ratio,0.,0.75,
        s_left, s_down, s_width, s_height)==0){
          return;
      }
      xtemp=ratio;
    }

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

     xtimeleft=85.0f*xtemp/(screenWidth-dwinWW);
     xtimeright=xtimeleft+0.6*(xtemp-xtimeleft);

     if( visTimeLabels==1&&showtime==1){
      if(visTimelabel==1)outputText(0.0f,0.1f, timelabel);
      if(visFramelabel==1&&visHRRlabel==0)outputText(0.0f,0.4f, framelabel);
      if(visHRRlabel==1&&hrrinfo!=NULL)outputText(0.0f,0.4f, hrrinfo->hrrlabel);
      drawTimeBar();
     }

    if((benchmark==1||visFramerate==1)&&showtime==1
      ){
      if(benchmark==0){
        sprintf(frameratelabel," Frame rate:%4.1f",framerate);
      }
      else{
        if(framerate<0.0){
          strcpy(frameratelabel," *Frame rate:");
        }
        else{
          sprintf(frameratelabel," *Frame rate:%4.1f",framerate);
        }
      }
      outputText((float)(xtimeright+0.025),0.08f, frameratelabel);
    }
    if(show_slice_average==1&&vis_slice_average==1&&slice_average_flag==1){
      sprintf(frameratelabel," AVG: %4.1f",slice_average_interval);
      outputText((float)(xtimeright+0.025),0.56, frameratelabel); // test print
    }
    if(hrrpuv_loaded==1&&show_hrrcutoff==1&&current_mesh!=NULL){
      char hrrcut_label[256];
      int ihrrcut;
      float xxl, xxr, yyl, yyu, ddx=0.03, ddy=0.2;

      ihrrcut = (int)global_hrrpuv_cutoff;

      sprintf(hrrcut_label,">%i (kW/m3)",ihrrcut);
      outputText((float)(xtimeright+0.06),0.56f, hrrcut_label);
      xxl = xtimeright+0.025;
      xxr = xxl+ddx;
      yyl = 0.56;
      yyu = yyl + ddy;

      glBegin(GL_QUADS);
      glColor3f(fire_red/255.0,fire_green/255.0,fire_blue/255.0);
      glVertex3f(xxl,yyl,0.0);
      glVertex3f(xxr,yyl,0.0);
      glVertex3f(xxr,yyu,0.0);
      glVertex3f(xxl,yyu,0.0);
      glEnd();
    }
#ifdef pp_memstatus
    if(visAvailmemory==1
    ){
    
      MEMSTATUS(0,&availmemory,NULL,NULL);
      sprintf(frameratelabel," Mem Load:%u%s",availmemory,percen);
      if((benchmark==1||visFramerate==1)&&showtime==1){
        outputText((float)(xtimeright+0.025),0.32f, frameratelabel);
      }
      else{
        outputText((float)(xtimeright+0.025),0.08f, frameratelabel);
      }
    }
#endif
#ifdef pp_MEMDEBUG
    if(visUsagememory==1
#ifdef pp_memstatus
       &&visAvailmemory==0
#endif
      ){
        char MEMlabel[128];

        getMemusage(MMtotalmemory,MEMlabel);
        if((benchmark==1||visFramerate==1)&&showtime==1){
          outputText((float)(xtimeright+0.025),0.32f, MEMlabel);
        }
        else{
          outputText((float)(xtimeright+0.025),0.08f, MEMlabel);
        }
    }
#endif
  }
}

/* --------------------- COLOR BAR Viewport2 ------------------------- */

void COLORBAR_viewport2(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height){

  // visColorbarLabels
  // numColorbars
  // dwinH
  // screenHeight
  // screenWidth
  // titlesafe_offset

  if(visColorbarLabels==1&&numColorbars!=0){
    if(screenWidth<screenHeight){
  //    if(SUB_portortho(quad,)==0){
  //        return;
  //    }
    }
    else{
     // if(SUB_portortho(quad,)==0){
     //     return;
     // }
    }

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if( showtime==1 || showplot3d==1){
      drawColorBars();
    }
  }
}

/* --------------------- COLOR BAR Viewport ------------------------- */

void COLORBAR_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height){
  GLint temp;
  float xnum;

  if(visColorbarLabels==1&&numColorbars!=0){
    temp = (int)(1.2f*dwinH);
    xnum=numColorbars;
    if(fontindex==LARGE_FONT)xnum*=1.5;
    if(screenWidth<screenHeight){
      barright=xnum/3.0+0.1f;
      if(SUB_portortho(quad,
        screenWidth-2-dwinWW-fontWoffset-titlesafe_offset,
        temp+titlesafe_offset,
        dwinWW,
        screenHeight-temp-fontHoffset-2*titlesafe_offset,
        0.,(double)barright,-1.5,(double)(ratio*(nrgb+1)),
        s_left, s_down, s_width, s_height)==0){
          return;
      }
    }
    else{
      barright=(xnum/3.0+0.1f);
      if(SUB_portortho(quad,
        screenWidth-2-dwinWW-fontWoffset-titlesafe_offset,
        temp+titlesafe_offset,
        dwinWW,
        screenHeight-temp-fontHoffset-2*titlesafe_offset,
        0.,(double)barright,-1.5,(double)(nrgb+1),
        s_left, s_down, s_width, s_height)==0){
        return;
      }
    }

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if( showtime==1 || showplot3d==1){
      drawColorBars();
    }
  }
}

    /* -------------------------- TITLE Viewport -------------------------- */

void TITLE_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height){
  int left;
  float textdown;
  int box_width_pixels, box_height_pixels;
  float box_width, box_height;

  if(visTitle!=1)return;
  box_height_pixels=(int)(ntitles*dwinH/4);
  box_width_pixels=screenWidth-dwinWW-fontWoffset-2*titlesafe_offset;
  box_width=0.5;
  box_height=(ntitles*ratio);
  if(screenWidth<screenHeight){
    if(SUB_portortho(quad,
      fontWoffset+titlesafe_offset,
      (int)(screenHeight-1.1f*ntitles*dwinH/4.f-fontHoffset)-titlesafe_offset,
      box_width_pixels,box_height_pixels,
      0.,1.,0.,(double)(ntitles*ratio),
      s_left, s_down, s_width, s_height)==0)return;
    left=(int)((float)75/(float)(screenWidth-dwinWW));
  }
  else{
    if(SUB_portortho(quad,
      fontWoffset+titlesafe_offset,
      (int)(screenHeight-1.1f*ntitles*dwinH/4.f-fontHoffset)-titlesafe_offset,
      box_width_pixels,box_height_pixels,
      0.,1.,0.,(double)(ntitles*ratio),
      s_left, s_down, s_width, s_height)==0)return;
    left=(int)((float)75/(float)(screenWidth-dwinWW)*ratio);
  }
  textdown=ratio/5.0;

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  ititle=0;
  if(showtitle2==1){
    ititle++;
    outputText(left,(ititle-1)*ratio+textdown, TITLE2);
  }
  if(showtitle1==1){
    ititle++;
    outputText(left,(ititle-1)*ratio+textdown, TITLE1);
  }
  if(visTitle0==1){
    ititle++;
    if(visFullTitle==1&&showplot3d==1){
      outputText(left,(ititle-1)*ratio+textdown, FULLTITLE);
    }
    else{
      outputText(left,(ititle-1)*ratio+textdown, TITLE);
    }
  }
}

/* ----------------------- 3D scene Viewport ----------------------------- */

void Scene_viewport(int quad, int view_mode, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height){

  float up, down, left, right;
  float dh;
  float fleft, fright, fup, fdown;
  float StereoCameraOffset,FrustumAsymmetry;
  float aperture_temp;
  float widthdiv2;
  float eyexINI, eyeyINI, eyezINI;

  if(showstereo==2){
    down=0;
    left=s_left;
    up=down+s_height;
    right=left+s_width;
  }
  else{
    right=screenWidth;
    left=0.;
    up=screenHeight;
    down=0;
    if(visColorbarLabels==1||(visBlocklabel==1&&nmeshes>1))right=screenWidth-dwinWW;
    if(visTitle==1)up=screenHeight-1.1*ntitles*dwinH/4.0;
  }

  if((visTimeLabels==1&&showtime==1)||(showtime==1&&(visFramerate==1||benchmark==1))||(visGrid!=noGridnoProbe&&visgridloc==1)
#ifdef pp_memstatus
      ||visAvailmemory==1
#endif
      ||(hrrpuv_loaded==1&&show_hrrcutoff==1&&current_mesh!=NULL)
    ){
    down=0.75*dwinH;
  }
  down += fontHoffset;
  dh = (up-down)*2*fontWoffset/(right-left);
  
  aspect=(float)(up-down)/(float)(right-left);

  /* set view position for virtual tour */

  {
    tourdata *touri;
    pathdata *pj;
    if(plotstate==DYNAMIC_PLOTS&&selected_tour!=NULL&&selected_tour->timeslist!=NULL){
      if((viewtourfrompath==1&&selectedtour_index>=0)||keyframe_snap==1){
        touri = tourinfo + selectedtour_index;
        iframe = touri->timeslist[itimes];
        if(keyframe_snap==1&&selected_frame!=NULL){
          pj=&selected_frame->nodeval;
        }
        else{
          pj = touri->pathnodes + iframe;
        }

        camera_current->eye[0]=pj->eye[0];
        camera_current->eye[1]=pj->eye[1];
        camera_current->eye[2]=pj->eye[2];
        camera_current->az_elev[1]=0.0;
        camera_current->az_elev[0]=0.0;

      }
    }
  }

  if(plotstate==DYNAMIC_PLOTS&&select_avatar==1&&selected_avatar_tag>0&&view_from_selected_avatar==1){
    camera_current->eye[0]=selected_avatar_pos[0];
    camera_current->eye[1]=selected_avatar_pos[1];
    camera_current->eye[2]=selected_avatar_pos[2];
    camera_current->azimuth=selected_avatar_angle;
    camera_current->view_angle=0.0;
    update_camera(camera_current);
    //camera_current->az_elev[1]=0.0;
    //camera_current->az_elev[0]=0.0;
  }

  eyexINI = camera_current->eye[0];
  eyeyINI = camera_current->eye[1];
  eyezINI = camera_current->eye[2];

  fnear =  - eyeyINI-1.0;
  if(fnear<nearclip)fnear=nearclip;
  ffar = fnear + farclip;
  
  StereoCameraOffset=0.0;
  aperture_temp=aperture;
//  zoom = tan(PI*aperture/360.0)/sqrt(3.0)
//  aperture_temp = 360.0*atan(1.0/(zoom*sqrt(3.0)))/PI;
//  aperture_temp = (360.0/PI)*atan(0.5/zoom);
//  aperture_temp = (360.0/PI)*atan((8.5/20.0)/zoom);
  aperture_temp = zoom2aperture(zoom);

  {
    tourdata *touri;
    pathdata *pj;

    if(plotstate==DYNAMIC_PLOTS&&selected_tour!=NULL&&selected_tour->timeslist!=NULL){
      if((viewtourfrompath==1&&selectedtour_index>=0)||keyframe_snap==1){
        touri = tourinfo + selectedtour_index;
        iframe = touri->timeslist[itimes];
        if(keyframe_snap==1&&selected_frame!=NULL){
          pj=&selected_frame->nodeval;
        }
        else{
          pj = touri->pathnodes + iframe;
        }

        aperture_temp=zoom2aperture(pj->zoom);
      }
    }
  }


  widthdiv2 = fnear*tan(0.5*aperture_temp*DEG2RAD);
  fleft = -widthdiv2;
  fright = widthdiv2;
  fup = aspect*widthdiv2;
  fdown = -aspect*widthdiv2;
  
  if(showstereo==0||view_mode==VIEW_CENTER){
    StereoCameraOffset=0.0;
    FrustumAsymmetry=0.0;
  }
  else if(showstereo!=0&&(view_mode==VIEW_LEFT||view_mode==VIEW_RIGHT)){
    StereoCameraOffset = (fzero/xyzmaxdiff)/30.0;
    if(view_mode==VIEW_LEFT)StereoCameraOffset = -StereoCameraOffset;
    FrustumAsymmetry= -0.5*StereoCameraOffset*fnear/(fzero/xyzmaxdiff);
  }

  if(SUB_portfrustum(quad,
    (int)(left+fontWoffset+titlesafe_offset),
    (int)(down+titlesafe_offset),
    (int)(right-left-2*fontWoffset-2*titlesafe_offset),
    (int)(up-down-dh-2*titlesafe_offset),
    (double)(fleft+FrustumAsymmetry),(double)(fright+FrustumAsymmetry),
    (double)fdown,(double)fup,
    (double)fnear,(double)ffar,
      s_left, s_down, s_width, s_height)==0)return;

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  {
    float sin_dv_sum, cos_dv_sum;
    float sin_azimuth, cos_azimuth;
    float sn_view_angle, cs_view_angle;
    float *uup;
    float cos_elevation, sin_elevation;
    float xcen, ycen, zcen;
    float posx, posy, posz;
    float azimuth, elevation;

    sn_view_angle=sin(DEG2RAD*camera_current->view_angle);
    cs_view_angle=cos(DEG2RAD*camera_current->view_angle);

    sin_azimuth=sin(DEG2RAD*camera_current->azimuth);
    cos_azimuth=cos(DEG2RAD*camera_current->azimuth);

    xcen = camera_current->xcen;
    ycen = camera_current->ycen;
    zcen = camera_current->zcen;

    cos_elevation=cos(DEG2RAD*camera_current->elevation);
    sin_elevation=sin(DEG2RAD*camera_current->elevation);

    sin_dv_sum = sin_azimuth*cs_view_angle + cos_azimuth*sn_view_angle;
    cos_dv_sum = cos_azimuth*cs_view_angle - sin_azimuth*sn_view_angle;


    posx = eyexINI+StereoCameraOffset*cos_dv_sum;
    posy = eyeyINI-StereoCameraOffset*sin_dv_sum;
    posz = eyezINI;

    viewx = posx + sin_dv_sum*cos_elevation;
    viewy = posy + cos_dv_sum*cos_elevation;
    viewz = posz + sin_elevation;

    elevation = camera_current->az_elev[1];
    azimuth = camera_current->az_elev[0];
    
    /* set view direction for virtual tour */
    {
      tourdata *touri;
      pathdata *pj;

      if(plotstate==DYNAMIC_PLOTS&&selected_tour!=NULL&&selected_tour->timeslist!=NULL){
        if((viewtourfrompath==1&&selectedtour_index>=0)||keyframe_snap==1){
          touri = tourinfo + selectedtour_index;
          iframe = touri->timeslist[itimes];
          if(keyframe_snap==1&&selected_frame!=NULL){
            pj=&selected_frame->nodeval;
          }
          else{
            pj = touri->pathnodes + iframe;
          }

          viewx = pj->oview[0]+StereoCameraOffset*cos_dv_sum;
          viewy = pj->oview[1]-StereoCameraOffset*sin_dv_sum;
          viewz = pj->oview[2];
          elevation=0.0;
          azimuth=0.0;
        }
      }
    }

    uup = camera_current->up;
    gluLookAt(
      (double)(posx), (double)(posy), (double)posz,
      (double)viewx,  (double)viewy,  (double)viewz,
      uup[0], uup[1], uup[2]);

    glGetFloatv(GL_MODELVIEW_MATRIX,modelview_setup);
    getinverse(modelview_setup,inverse_modelview_setup);

    glMultMatrixf(modelview_identity);
    
    glTranslatef(xcen,ycen,zcen);

    // rotate scene
    
    if(rotation_type==ROTATION_3AXIS){
      glMultMatrixf(quat_rotation);
    }
    else{
      if(rotation_type==ROTATION_2AXIS){
        glRotatef(elevation,1.0,0.0,0.0);  /* rotate about x axis */
      }
      glRotatef(azimuth,0.0,0.0,1.0);      /* rotate about z axis */
    }
    
    glTranslatef(-xcen,-ycen,-zcen);

    glGetFloatv(GL_MODELVIEW_MATRIX,modelview_scratch);
    matmatmult(inverse_modelview_setup,modelview_scratch,modelview_current);

    get_world_eyepos(modelview_scratch, world_eyepos,scaled_eyepos);

    if(show_gslice_triangles==1||vis_gslice_data==1)update_gslice_planes();
    if(nrooms>0){
      getzonesmokedir(modelview_scratch);
    }
    if(nvolrenderinfo>0&&showvolrender==1){
      if(usevolrender==1){
        getvolsmokedir(modelview_scratch);
        SNIFF_ERRORS("after getvolsmokedir");
#ifdef pp_GPU
        if(usegpu==0){
          compute_all_smokecolors();
        }
#else
        compute_all_smokecolors();
#endif
      }
    }
    if(nsmoke3dinfo>0&&show3dsmoke==1){
      sort_smoke3dinfo();
      getsmokedir(modelview_scratch);
      SNIFF_ERRORS("after getsmokedir");
#ifdef pp_CULL
      if(showstereo==0){
        if(cullsmoke==1){
          getPixelCount();
          SNIFF_ERRORS("after getPixelCount");
        }
        if(cullactive==1&&update_initcullplane==1){
          initcullplane(cullsmoke);
        }
        SNIFF_ERRORS("after initcullplane");
      }
#endif
    }
    if(nface_transparent>0&&sort_transparent_faces==1)Sort_Transparent_Faces(modelview_scratch);
    if(showiso==1)Update_Isotris(0);
    FREEMEMORY(geominfoptrs);
    ngeominfoptrs=0;
    GetGeomInfoPtrs(&geominfoptrs,&ngeominfoptrs);
    if(ngeominfoptrs>0)Sort_Embedded_Geometry(modelview_scratch);
    if(showiso==1&&sort_iso_triangles==1&&niso_trans>0)Sort_Iso_Triangles(modelview_scratch);

    glScalef(mscale[0],mscale[1],mscale[2]);
    ExtractFrustum();
    set_cull_vis();
  }
}
