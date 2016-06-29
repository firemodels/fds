#include "options.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include GLUT_H

#include "smokeviewvars.h"
#include "IOvolsmoke.h"

#define CONV(p,pl,pr,pxl,pxr) ( (pxl) + ((pxr)-(pxl))*((p)-(pl))/((pr)-(pl)) )

/* ------------------------ getStringWidth ------------------------- */

int getStringWidth(char *string){
  char *c;
  int length=0;

  if(string==NULL)return 0;
  switch(fontindex){
    case SMALL_FONT:
      length = strlen(string);
      length *= (288.0/235.0)*glutBitmapWidth(GLUT_BITMAP_HELVETICA_10, 'a');
      break;
    case LARGE_FONT:
      length = strlen(string);
      length *= (416.0/423.0)*glutBitmapWidth(GLUT_BITMAP_HELVETICA_18, 'a');
      break;
    case SCALED_FONT:
      for(c=string;*c!='\0';c++){
        length += glutStrokeWidth(GLUT_STROKE_ROMAN, *c);
      }
      length *= (283.0/402.0)*scale_2d_x;
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
  return length;
}

/* ------------------------ get_viewport_info ------------------------- */

void get_viewport_info(void){
  int doit;
  float text_height;
  float text_width;
  int ninfo_lines=0;
  int info_width;

  info_width = getStringWidth("y: 115, 11.5 m");
  colorbar_label_width = getStringWidth("*10^-02");

  v_space = 2;
  text_height=18;
  text_width=18;
  if(fontindex==SCALED_FONT){
    scale_2d_x = (scaled_font2d_height2width*(float)scaled_font2d_height/(float)104.76);
    scale_2d_y = ((float)scaled_font2d_height/(float)152.38);

    text_height = MAX(18,(int)( (12.0/18.0)*(25.0/18.0)*(float)scaled_font2d_height));
    text_width =  MAX(18, (25.0/36.0)*(scaled_font2d_height2width*(float)scaled_font2d_height));
  }

  // full screen viewport dimensions

  VP_fullscreen.left = 0;
  VP_fullscreen.down = 0;
  VP_fullscreen.width = screenWidth;
  VP_fullscreen.height = screenHeight;
  VP_fullscreen.right = VP_fullscreen.left + VP_fullscreen.width;
  VP_fullscreen.top = VP_fullscreen.down + VP_fullscreen.height;
  VP_info.doit = 1;

  // INFO viewport dimensions

  doit=0;
  if(visMeshlabel==1){
    ninfo_lines++;
    doit=1;
  }
  if(((showplot3d==1||visGrid!=noGridnoProbe)&&visx_all==1)||visGrid==noGridProbe||visGrid==GridProbe){
    if(visgridloc==1){
      ninfo_lines++;
      doit=1;
    }
  }
  if(((showplot3d==1||visGrid!=noGridnoProbe)&&visy_all==1)||visGrid==GridProbe||visGrid==noGridProbe){
    if(visgridloc==1){
      ninfo_lines++;
      doit=1;
    }
  }
  if(((showplot3d==1||visGrid!=noGridnoProbe)&&visz_all==1)||visGrid==GridProbe||visGrid==noGridProbe){
    if(visgridloc==1){
      ninfo_lines++;
      doit=1;
    }
  }

  VP_info.left = screenWidth-info_width-titlesafe_offset;
  VP_info.down = titlesafe_offset;
  VP_info.doit = doit;
  VP_info.text_height = text_height;
  VP_info.text_width = text_width;
  if(doit==1){
    VP_info.width = info_width;
    VP_info.height = ninfo_lines*(text_height+v_space);
  }
  else{
    VP_info.width = 0;
    VP_info.height = 0;
  }
  VP_info.right = VP_info.left + VP_fullscreen.width;
  VP_fullscreen.top = VP_fullscreen.down + VP_info.height;

  // timebar viewport dimensions

  doit=0;
  if(
    (visTimebar==1&&showtime==1)||
    (showtime==1&&(visFramerate==1||(vis_slice_average==1&&show_slice_average&&slice_average_flag==1))||
    (hrrpuv_loaded==1&&show_hrrcutoff==1&&current_mesh!=NULL)
    )
#ifdef pp_memstatus
    ||visAvailmemory==1
#endif
    )doit=1;

  VP_timebar.left = titlesafe_offset;
  VP_timebar.down = titlesafe_offset;
  VP_timebar.doit=doit;
  VP_timebar.text_height=text_height;
  VP_timebar.text_width = text_width;
  if(doit==1){
    VP_timebar.width = screenWidth-VP_info.width-2*titlesafe_offset;
    VP_timebar.height=2*(text_height+v_space);
    if(hrrpuv_loaded==1&&show_hrrcutoff==1&&current_mesh!=NULL)VP_timebar.height=3*(text_height+v_space);
  }
  else{
    VP_timebar.width = 0;
    VP_timebar.height = 0;
  }
  VP_timebar.right = VP_timebar.left + VP_timebar.width;
  VP_timebar.top = VP_timebar.down + VP_timebar.height;

  // colorbar viewport dimensions

  doit=1;
  if(visColorbar==0||numColorbars==0||(showtime==0&&showplot3d==0))doit=0;
  VP_colorbar.left = screenWidth-colorbar_delta - numColorbars*(colorbar_label_width+2*h_space)-titlesafe_offset;
  VP_colorbar.down = MAX(VP_timebar.height,VP_info.height)+titlesafe_offset;
  VP_colorbar.doit = doit;
  VP_colorbar.text_height=text_height;
  VP_colorbar.text_width = text_width;
  if(doit==1){
    VP_colorbar.width = colorbar_delta + h_space+numColorbars*(colorbar_label_width+h_space);
    VP_colorbar.height = screenHeight-MAX(VP_timebar.height,VP_info.height)-2*titlesafe_offset;

  }
  else{
    VP_colorbar.width = 0;
    VP_colorbar.height = 0;
  }
  VP_colorbar.right = VP_colorbar.left+VP_colorbar.width;
  VP_colorbar.top = VP_colorbar.down+VP_colorbar.height;

    // title viewport dimensions

  if(visTitle==1){
    VP_title.width = screenWidth-VP_colorbar.width-2*titlesafe_offset;
    if(gversion==1){
      VP_title.height=3*text_height+2*v_space;
    }
    else{
      VP_title.height=text_height+v_space;
    }
    VP_title.doit = 1;
  }
  else{
    VP_title.width = 0;
    VP_title.height = 0;
    VP_title.doit = 0;
  }
  VP_title.text_height=text_height;
  VP_title.text_width = text_width;
  VP_title.left = titlesafe_offset;
  VP_title.down = (int)screenHeight-VP_title.height-titlesafe_offset;
  VP_title.right = VP_title.left + VP_title.width;
  VP_title.top = VP_title.down + VP_title.height;

  // scene viewport dimensions

  VP_scene.text_height = text_height;
  VP_scene.text_width = text_width;
  VP_scene.left=titlesafe_offset;
  VP_scene.down=titlesafe_offset+MAX(VP_timebar.height,VP_info.height);
  VP_scene.width=screenWidth-2*titlesafe_offset-VP_colorbar.width;
  VP_scene.height=screenHeight-MAX(VP_timebar.height,VP_info.height)-VP_title.height - 2*titlesafe_offset;
  VP_scene.right = VP_scene.left + VP_scene.width;
  VP_scene.top = VP_scene.down + VP_scene.height;

  scene_aspect_ratio = (float)VP_scene.height/(float)VP_scene.width;

  colorbar_right_pos = VP_colorbar.right-h_space;
  colorbar_left_pos = colorbar_right_pos - colorbar_delta;
  colorbar_top_pos = VP_colorbar.top - 4*(v_space + VP_colorbar.text_height) - colorbar_delta;
  colorbar_down_pos = VP_colorbar.down + colorbar_delta;

}

 /* ------------------------ SUB_portortho ------------------------- */

int SUB_portortho(int quad,
                  portdata *p,
                   GLdouble portx_left, GLdouble portx_right, GLdouble portx_down, GLdouble portx_top,
                   GLint screen_left, GLint screen_down
                   ){

  GLint subport_left, subport_right, subport_down, subport_top;
  GLdouble subportx_left, subportx_right, subportx_down, subportx_top;
  GLsizei subport_width, subport_height;
  GLint subwindow_left, subwindow_right, subwindow_down, subwindow_top;
  GLint port_right, port_top;

  int irow, icol;

  switch(quad){
  case 0:
    port_pixel_width = p->width;
    port_pixel_height = p->height;
    port_unit_width = portx_right - portx_left;
    port_unit_height = portx_top - portx_down;
    glViewport(p->left,p->down,p->width,p->height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(portx_left,portx_right,portx_down,portx_top);
    return 1;
  case 1:
    icol = screen_left/screenWidth;
    irow = screen_down/screenHeight;

    subwindow_left = icol*screenWidth;
    subwindow_right = subwindow_left + screenWidth;
    subwindow_down = irow*screenHeight;
    subwindow_top = subwindow_down + screenHeight;

    port_right = p->left + p->width;
    port_top = p->down + p->height;

    subport_left =  MAX( nrender_rows*p->left,subwindow_left);
    subport_right = MIN(nrender_rows*port_right,subwindow_right);
    subport_down =  MAX( nrender_rows*p->down,subwindow_down);
    subport_top =   MIN(  nrender_rows*port_top,subwindow_top);
    if(subport_left>=subport_right||subport_down>=subport_top)return 0;

    subportx_left = CONV(subport_left,nrender_rows*p->left,nrender_rows*port_right,portx_left,portx_right);
    subportx_right = CONV(subport_right,nrender_rows*p->left,nrender_rows*port_right,portx_left,portx_right);
    subportx_down = CONV(subport_down,nrender_rows*p->down,nrender_rows*port_top,portx_down,portx_top);
    subportx_top = CONV(subport_top,nrender_rows*p->down,nrender_rows*port_top,portx_down,portx_top);

    subport_left -= icol*screenWidth;
    subport_right -= icol*screenWidth;
    subport_down -= irow*screenHeight;
    subport_top -= irow*screenHeight;
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


/* ------------------------ SUB_portortho2 ------------------------- */

int SUB_portortho2(int quad,
                  portdata *p,
                  GLint screen_left, GLint screen_down
                  ){

  GLint subport_left, subport_right, subport_down, subport_top;
  GLdouble subportx_left, subportx_right, subportx_down, subportx_top;
  GLsizei subport_width, subport_height;
  GLint subwindow_left, subwindow_right, subwindow_down, subwindow_top;
  GLint port_right, port_top;

  int irow, icol;
  GLdouble portx_left, portx_right, portx_down, portx_top;

  portx_left = p->left;
  portx_right = p->left + p->width;
  portx_down = p->down;
  portx_top = p->down + p->height;
  switch(quad){
  case 0:
    port_pixel_width = p->width;
    port_pixel_height = p->height;
    port_unit_width = portx_right - portx_left;
    port_unit_height = portx_top - portx_down;
    glViewport(p->left,p->down,p->width,p->height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(portx_left,portx_right,portx_down,portx_top);
    return 1;
  case 1:
    icol = screen_left/screenWidth;
    irow = screen_down/screenHeight;

    subwindow_left = icol*screenWidth;
    subwindow_right = subwindow_left + screenWidth;
    subwindow_down = irow*screenHeight;
    subwindow_top = subwindow_down + screenHeight;

    port_right = p->left + p->width;
    port_top = p->down + p->height;

    subport_left =  MAX( nrender_rows*p->left,subwindow_left);
    subport_right = MIN(nrender_rows*port_right,subwindow_right);
    subport_down =  MAX( nrender_rows*p->down,subwindow_down);
    subport_top =   MIN(  nrender_rows*port_top,subwindow_top);
    if(subport_left>=subport_right||subport_down>=subport_top)return 0;

    subportx_left = CONV(subport_left,nrender_rows*p->left,nrender_rows*port_right,portx_left,portx_right);
    subportx_right = CONV(subport_right,nrender_rows*p->left,nrender_rows*port_right,portx_left,portx_right);
    subportx_down = CONV(subport_down,nrender_rows*p->down,nrender_rows*port_top,portx_down,portx_top);
    subportx_top = CONV(subport_top,nrender_rows*p->down,nrender_rows*port_top,portx_down,portx_top);

    subport_left -= icol*screenWidth;
    subport_right -= icol*screenWidth;
    subport_down -= irow*screenHeight;
    subport_top -= irow*screenHeight;
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
                   portdata *p,
                   GLdouble portx_left, GLdouble portx_right,
                   GLdouble portx_down, GLdouble portx_top,
                   GLdouble portx_near, GLdouble portx_far,
                   GLint screen_left, GLint screen_down
                   ){
  GLint subport_left, subport_right, subport_down, subport_top;
  GLdouble subportx_left, subportx_right, subportx_down, subportx_top;
  GLsizei subport_width, subport_height;
  GLint subwindow_left, subwindow_right, subwindow_down, subwindow_top;
  GLint port_right, port_top;

  int irow, icol;

  switch(quad){
  case 0:
    port_pixel_width = p->width;
    port_pixel_height = p->height;
    port_unit_width = portx_right - portx_left;
    port_unit_height = portx_top - portx_down;
    glViewport(p->left,p->down,p->width,p->height);
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
    icol = screen_left/screenWidth;
    irow = screen_down/screenHeight;

    subwindow_left = icol*screenWidth;
    subwindow_right = subwindow_left + screenWidth;
    subwindow_down = irow*screenHeight;
    subwindow_top = subwindow_down + screenHeight;

    port_right = p->left + p->width;
    port_top = p->down + p->height;

    subport_left =  MAX( nrender_rows*p->left,subwindow_left);
    subport_right = MIN(nrender_rows*port_right,subwindow_right);
    subport_down =  MAX( nrender_rows*p->down,subwindow_down);
    subport_top =   MIN(  nrender_rows*port_top,subwindow_top);
    if(subport_left>=subport_right||subport_down>=subport_top)return 0;

    subportx_left = CONV(subport_left,nrender_rows*p->left,nrender_rows*port_right,portx_left,portx_right);
    subportx_right = CONV(subport_right,nrender_rows*p->left,nrender_rows*port_right,portx_left,portx_right);
    subportx_down = CONV(subport_down,nrender_rows*p->down,nrender_rows*port_top,portx_down,portx_top);
    subportx_top = CONV(subport_top,nrender_rows*p->down,nrender_rows*port_top,portx_down,portx_top);

    subport_left -= icol*screenWidth;
    subport_right -= icol*screenWidth;
    subport_down -= irow*screenHeight;
    subport_top -= irow*screenHeight;
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

void CLIP_viewport(int quad, GLint screen_left, GLint screen_down){
  GLdouble x_left, x_right, x_down, x_top;
  float c_left, c_right, c_top, c_bottom;

  x_left=0.0;
  x_right=screenWidth;
  x_down=0.0;
  x_top=screenHeight;

  if(SUB_portortho(quad,&VP_fullscreen,x_left, x_right, x_down, x_top,screen_left, screen_down)==0)return;

   c_left = render_clip_left-3;
   c_right = screenWidth + 3 - render_clip_right;
   c_bottom = render_clip_bottom -3;
   c_top = screenHeight + 3 - render_clip_top;

   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   glLineWidth(3.0);
   glColor3fv(foregroundcolor);
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

 /* ------------------------ INFO_viewport ------------------------- */

void INFO_viewport(int quad, GLint screen_left, GLint screen_down){
  char slicelabel[255];
  meshdata *mesh_xyz=NULL;
  float xyz[3];
  int info_lines=0;

  if(SUB_portortho2(quad,&VP_info,screen_left, screen_down)==0)return;

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  if((showplot3d==1||visGrid!=noGridnoProbe)&&(visx_all==1||visy_all||visz_all)||visGrid==GridProbe||visGrid==noGridProbe){
    xyz[0]=DENORMALIZE_X(plotx_all[iplotx_all]);
    xyz[1]=DENORMALIZE_Y(ploty_all[iploty_all]);
    xyz[2]=DENORMALIZE_Z(plotz_all[iplotz_all]);
    mesh_xyz=getmesh_nofail(xyz);
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
    if(visgridloc==1){
      outputText(VP_info.left+h_space,VP_info.down+v_space, slicelabel);
      info_lines++;
    }
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
    if(visgridloc==1){
      outputText(VP_info.left+h_space,VP_info.down+v_space+info_lines*(v_space+VP_info.text_height), slicelabel);
      info_lines++;
    }
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
    if(visgridloc==1){
      outputText(VP_info.left+h_space,VP_info.down+v_space+info_lines*(v_space+VP_info.text_height), slicelabel);
      info_lines++;
    }
  }
  if(visMeshlabel==1){
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
    outputText(VP_info.left+h_space,VP_info.down+v_space+info_lines*(v_space+VP_info.text_height), meshlabel);
  }
}

/* ------------------------ TIMEBAR_viewport ------------------------- */

void TIMEBAR_viewport(int quad, GLint screen_left, GLint screen_down){
#ifdef pp_memstatus
  unsigned int availmemory;
  char percen[]="%";
#endif
  int right_label_pos,timebar_right_pos;
  int timebar_left_pos;

  if(SUB_portortho2(quad,&VP_timebar,screen_left,screen_down)==0)return;

  timebar_left_width = getStringWidth("Time: 1234.11");
  timebar_right_width = getStringWidth("Frame rate: 99.99");

  timebar_left_pos = VP_timebar.left+timebar_left_width;
  timebar_right_pos= VP_timebar.right-timebar_right_width-h_space;
  right_label_pos  = timebar_right_pos+h_space;

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  if( visTimebar==1&&showtime==1){
    if(visTimelabel==1){
      outputText(VP_timebar.left,v_space, timelabel);
    }
    if(visFramelabel==1&&(visHRRlabel==0||hrrinfo==NULL)){
      outputText(VP_timebar.left,v_space+VP_timebar.text_height+v_space, framelabel);
    }
    if(visHRRlabel==1&&hrrinfo!=NULL){
      outputText(VP_timebar.left,v_space+VP_timebar.text_height+v_space, hrrinfo->hrrlabel);
    }
    drawTimeBar(timebar_left_pos,timebar_right_pos,v_space+VP_timebar.down,v_space+(VP_timebar.down+20));
  }

  if(visFramerate==1&&showtime==1){
    sprintf(frameratelabel," Frame rate:%4.1f",framerate);
    outputText(right_label_pos,v_space,frameratelabel);
  }
  if(show_slice_average==1&&vis_slice_average==1&&slice_average_flag==1){
    sprintf(frameratelabel," AVG: %4.1f",slice_average_interval);
    outputText(right_label_pos,3*v_space+2*VP_timebar.text_height, frameratelabel); // test print
  }
  if(hrrpuv_loaded==1&&show_hrrcutoff==1&&current_mesh!=NULL){
    char hrrcut_label[256];
    int ihrrcut;

    ihrrcut = (int)global_hrrpuv_cutoff;

    sprintf(hrrcut_label,">%i (kW/m3)",ihrrcut);
    outputText(right_label_pos+5+h_space,3*v_space+2*VP_timebar.text_height,hrrcut_label);

    glBegin(GL_QUADS);
    glColor3f(fire_red/255.0,fire_green/255.0,fire_blue/255.0);
    glVertex3f(right_label_pos+h_space-20,5+2*VP_timebar.text_height   ,0.0);
    glVertex3f(right_label_pos+h_space   ,5+2*VP_timebar.text_height   ,0.0);
    glVertex3f(right_label_pos+h_space   ,5+2*VP_timebar.text_height+20,0.0);
    glVertex3f(right_label_pos+h_space-20,5+2*VP_timebar.text_height+20,0.0);
    glEnd();
  }
#ifdef pp_memstatus
  if(visAvailmemory==1){
    MEMSTATUS(0,&availmemory,NULL,NULL);
    sprintf(frameratelabel," Mem Load:%u%s",availmemory,percen);
    if(visFramerate==1&&showtime==1){
      outputText(right_label_pos,2*v_space+VP_timebar.text_height,frameratelabel);
    }
    else{
      outputText(right_label_pos,v_space,frameratelabel);
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
      if(visFramerate==1&&showtime==1){
        outputText(right_label_pos,2*v_space+VP_timebar.text_height,MEMlabel);
      }
      else{
        outputText(right_label_pos,v_space,MEMlabel);
      }
  }
#endif
}

/* --------------------- COLORBAR_viewport ------------------------- */

void COLORBAR_viewport(int quad, GLint screen_left, GLint screen_down){
  if(SUB_portortho2(quad,&VP_colorbar,screen_left, screen_down)==0)return;

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  drawColorBars();
}

    /* -------------------------- TITLE_viewport -------------------------- */

void TITLE_viewport(int quad, GLint screen_left, GLint screen_down){
  float left, textdown;

  if(SUB_portortho2(quad,&VP_title,screen_left,screen_down)==0)return;

  left=0;
  textdown=VP_title.down;

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  if(gversion==0){
    if(visFullTitle==1&&showplot3d==1){
      outputText(left,textdown, plot3d_title);
    }
    else{
      outputText(left,textdown, release_title);
    }
  }
  else{
    char label[256];
    int smv_top, smv_top2, fds_top;

    if(fds_githash!=NULL){
      fds_top=textdown;
      smv_top=fds_top+VP_title.text_height+v_space;
      smv_top2=smv_top+VP_title.text_height+v_space;
    }
    else{
      smv_top=textdown;
      smv_top2=smv_top+VP_title.text_height+v_space;
    }
    outputText(left,smv_top2,release_title);
    sprintf(label,"Smokeview (64 bit) build: %s",smv_githash);
    outputText(left,smv_top,label);
    if(fds_githash!=NULL){
      sprintf(label,"FDS build:%s",fds_githash);
      outputText(left,fds_top,label);
    }
  }
}

/* ----------------------- compare_meshes ----------------------------- */

int compare_meshes(const void *arg1, const void *arg2){
  smoke3ddata *smoke3di, *smoke3dj;
  meshdata *meshi, *meshj;
  float *xyzmini, *xyzmaxi;
  float *xyzminj, *xyzmaxj;
  int dir = 0;
  int returnval = 0;

  smoke3di = *(smoke3ddata **)arg1;
  smoke3dj = *(smoke3ddata **)arg2;
  meshi = meshinfo + smoke3di->blocknumber;
  meshj = meshinfo + smoke3dj->blocknumber;
  if(meshi == meshj)return 0;
  xyzmini = meshi->boxmin;
  xyzmaxi = meshi->boxmax;
  xyzminj = meshj->boxmin;
  xyzmaxj = meshj->boxmax;
  if(dir == 0){
    if(xyzmaxi[0] <= xyzminj[0])dir = 1;
    if(xyzmaxj[0] <= xyzmini[0])dir = -1;
  }
  if(dir == 0){
    if(xyzmaxi[1] <= xyzminj[1])dir = 2;
    if(xyzmaxj[1] <= xyzmini[1])dir = -2;
  }
  if(dir == 0){
    if(xyzmaxi[2] <= xyzminj[2])dir = 3;
    if(xyzmaxj[2] <= xyzmini[2])dir = -3;
  }
  switch(dir){
  case 0:
    returnval = 0;
    break;
  case XDIR:
    if(world_eyepos[0] < xyzmaxi[0]){
      returnval = 1;
    }
    else{
      returnval = -1;
    }
    break;
  case XDIRNEG:
    if(world_eyepos[0] < xyzmaxj[0]){
      returnval = -1;
    }
    else{
      returnval = 1;
    }
    break;
  case YDIR:
    if(world_eyepos[1] < xyzmaxi[1]){
      returnval = 1;
    }
    else{
      returnval = -1;
    }
    break;
  case YDIRNEG:
    if(world_eyepos[1] < xyzmaxj[1]){
      returnval = -1;
    }
    else{
      returnval = 1;
    }
    break;
  case ZDIR:
    if(world_eyepos[2] < xyzmaxi[2]){
      returnval = 1;
    }
    else{
      returnval = -1;
    }
    break;
  case ZDIRNEG:
    if(world_eyepos[2] < xyzmaxj[2]){
      returnval = -1;
    }
    else{
      returnval = 1;
    }
    break;
  default:
	  ASSERT(FFALSE);
	  break;
  }
  return returnval;
}

/* ------------------ sort_smoke3dinfo ------------------------ */

void sort_smoke3dinfo(void){
  if(nsmoke3dinfo > 1){
    qsort((meshdata **)smoke3dinfo_sorted, (size_t)nsmoke3dinfo, sizeof(smoke3ddata *), compare_meshes);
  }
}

/* ----------------------- Scene_viewport ----------------------------- */

void Scene_viewport(int quad, int view_mode, GLint screen_left, GLint screen_down, float *view_dir){

  float fleft, fright, fup, fdown;
  float StereoCameraOffset,FrustumAsymmetry;
  float aperture_temp;
  float widthdiv2;
  float eyexINI, eyeyINI, eyezINI;

  if(showstereo==STEREO_LR){
    VP_scene.left=screen_left;
    VP_scene.width=screenWidth;
  }
  if(plotstate==DYNAMIC_PLOTS&&selected_tour!=NULL&&selected_tour->timeslist!=NULL){
    if((viewtourfrompath==1&&selectedtour_index>=0)||keyframe_snap==1){
      tourdata *touri;
      pathdata *pj;

      touri = tourinfo + selectedtour_index;
      frame_index = touri->timeslist[itimes];
      if(keyframe_snap==1&&selected_frame!=NULL){
        pj=&selected_frame->nodeval;
      }
      else{
        pj = touri->pathnodes + frame_index;
      }

      camera_current->eye[0]=pj->eye[0];
      camera_current->eye[1]=pj->eye[1];
      camera_current->eye[2]=pj->eye[2];
      camera_current->az_elev[1]=0.0;
      camera_current->az_elev[0]=0.0;
     }
  }

  if(plotstate==DYNAMIC_PLOTS&&select_avatar==1&&selected_avatar_tag>0&&view_from_selected_avatar==1){
    camera_current->eye[0]=selected_avatar_pos[0];
    camera_current->eye[1]=selected_avatar_pos[1];
    camera_current->eye[2]=selected_avatar_pos[2];
    camera_current->azimuth=selected_avatar_angle;
    camera_current->view_angle=0.0;
    update_camera(camera_current);
  }

  eyexINI = camera_current->eye[0];
  eyeyINI = camera_current->eye[1];
  eyezINI = camera_current->eye[2];

  fnear =  - eyeyINI-1.0;
  if(fnear<nearclip)fnear=nearclip;
  ffar = fnear + farclip;

  FrustumAsymmetry=0.0;
  StereoCameraOffset=0.0;
  aperture_temp=aperture;
  aperture_temp = zoom2aperture(zoom);

  if(plotstate==DYNAMIC_PLOTS&&selected_tour!=NULL&&selected_tour->timeslist!=NULL){
    if((viewtourfrompath==1&&selectedtour_index>=0)||keyframe_snap==1){
      tourdata *touri;
      pathdata *pj;

      touri = tourinfo + selectedtour_index;
      frame_index = touri->timeslist[itimes];
      if(keyframe_snap==1&&selected_frame!=NULL){
        pj=&selected_frame->nodeval;
      }
      else{
        pj = touri->pathnodes + frame_index;
      }

      aperture_temp=zoom2aperture(pj->zoom);
    }
  }

  widthdiv2 = fnear*tan(0.5*aperture_temp*DEG2RAD);
  fleft = -widthdiv2;
  fright = widthdiv2;
  fup = scene_aspect_ratio*widthdiv2;
  fdown = -scene_aspect_ratio*widthdiv2;

  if(showstereo==STEREO_NONE||view_mode==VIEW_CENTER){
    StereoCameraOffset=0.0;
    FrustumAsymmetry=0.0;
  }
  else if(showstereo!=STEREO_NONE&&(view_mode==VIEW_LEFT||view_mode==VIEW_RIGHT)){
    StereoCameraOffset = SCALE2SMV(fzero)/30.0;
    if(view_mode==VIEW_LEFT)StereoCameraOffset = -StereoCameraOffset;
    FrustumAsymmetry= -0.5*StereoCameraOffset*fnear/SCALE2SMV(fzero);
  }

  if(SUB_portfrustum(quad,&VP_scene,
    (double)(fleft+FrustumAsymmetry),(double)(fright+FrustumAsymmetry),(double)fdown,(double)fup,(double)fnear,(double)ffar,
    screen_left, screen_down)==0)return;

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
    float az0=0.0, elev0=0.0;

    if (view_dir != NULL) {
      az0 = view_dir[0];
      elev0 = view_dir[1];
    }

    sn_view_angle=sin(DEG2RAD*camera_current->view_angle);
    cs_view_angle=cos(DEG2RAD*camera_current->view_angle);

    sin_azimuth=sin(DEG2RAD*(camera_current->azimuth + az0));
    cos_azimuth=cos(DEG2RAD*(camera_current->azimuth + az0));

    xcen = camera_current->xcen;
    ycen = camera_current->ycen;
    zcen = camera_current->zcen;

    cos_elevation=cos(DEG2RAD*(camera_current->elevation + elev0));
    sin_elevation=sin(DEG2RAD*(camera_current->elevation + elev0));

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
          frame_index = touri->timeslist[itimes];
          if(keyframe_snap==1&&selected_frame!=NULL){
            pj=&selected_frame->nodeval;
          }
          else{
            pj = touri->pathnodes + frame_index;
          }

          viewx = pj->tour_view[0]+StereoCameraOffset*cos_dv_sum;
          viewy = pj->tour_view[1]-StereoCameraOffset*sin_dv_sum;
          viewz = pj->tour_view[2];
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
    {
      float u[3], axis[3], angle;

      u[0] = 0.0;
      u[1] = 0.0;
      u[2] = 1.0;
      rotateu2v(user_zaxis, u, axis, &angle);
      glRotatef(RAD2DEG*angle, axis[0], axis[1], axis[2]);
      glRotatef(zaxis_angles[2], u[0], u[1], u[2]);
    }

    glTranslatef(-xcen,-ycen,-zcen);

    glGetFloatv(GL_MODELVIEW_MATRIX,modelview_scratch);
    matmatmult(inverse_modelview_setup,modelview_scratch,modelview_current);

    get_world_eyepos(modelview_scratch, world_eyepos,scaled_eyepos);

    if(show_gslice_triangles==1||SHOW_gslice_data==1){
      update_gslice_planes();
    }
    if(nrooms>0){
      getzonesmokedir(modelview_scratch);
    }
    if(nvolrenderinfo>0&&showvolrender==1&&usevolrender==1){
      getvolsmokedir(modelview_scratch);
      SNIFF_ERRORS("after getvolsmokedir");
#ifdef pp_GPU
      if(usegpu==0)compute_all_smokecolors();
#else
      compute_all_smokecolors();
#endif
    }
    if(nsmoke3dinfo>0&&show3dsmoke==1){
      sort_smoke3dinfo();
      getsmokedir(modelview_scratch);
      SNIFF_ERRORS("after getsmokedir");
#ifdef pp_CULL
      if(showstereo==STEREO_NONE){
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
    if(ngeominfoptrs>0)ShowHideSortGeometry(modelview_scratch);
    if(showiso==1&&sort_iso_triangles==1&&niso_trans>0)Sort_Iso_Triangles(modelview_scratch);

    glScalef(mscale[0],mscale[1],mscale[2]);
    ExtractFrustum();
    set_cull_vis();
  }
}
