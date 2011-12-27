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

#include "smokeviewvars.h"

 /* ------------------------ SUB_portortho ------------------------- */
 
int SUB_portortho(int quad, 
                   GLint i_left, GLint i_down, GLsizei i_width, GLsizei i_height,
                   GLdouble x_left, GLdouble x_right, GLdouble x_down, GLdouble x_top,
                   GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height
                   ){
  
  GLint n_left, n_down;
  GLsizei n_width, n_height;
  GLdouble nx_left, nx_right, nx_down, nx_top;

  switch (quad){
  case 0:            
    glViewport(i_left,i_down,i_width,i_height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(x_left,x_right,x_down,x_top);
    return 1;
  case 1:
    n_left = 2*i_left - s_left;
    n_down = 2*i_down - s_down;
    n_width = 2*i_width;
    n_height = 2*i_height;
    nx_left = x_left;
    nx_right = x_right;
    nx_down = x_down;
    nx_top = x_top;
    if(n_left+n_width<0||n_down+n_height<0)return 0;
    if(n_left>screenWidth||n_down>screenHeight)return 0;
    if(n_left<0){
      nx_left = x_left - n_left*(x_right-x_left)/n_width;
      n_width = n_left + n_width;
      n_left = 0;
    }
    if(n_left+n_width>screenWidth){
      nx_right = x_left + (screenWidth-n_left)*(x_right-x_left)/n_width;
      n_width = screenWidth-n_left;
    }
    if(n_down<0){
      nx_down = x_down - n_down*(x_top-x_down)/n_height;
      n_height = n_down + n_height;
      n_down = 0;
    }
    if(n_down+n_height>screenHeight){
      nx_top = x_down + (screenHeight-n_down)*(x_top-x_down)/n_height;
      n_height = screenHeight-n_down;
    }
    glViewport(n_left,n_down,n_width,n_height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(nx_left,nx_right,nx_down,nx_top);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  return 1;
}

/* ------------------------ SUB_portfrustum ------------------------- */
 
int SUB_portfrustum(int quad, 
                   GLint i_left, GLint i_down, GLsizei i_width, GLsizei i_height,
                   GLdouble x_left, GLdouble x_right, 
                   GLdouble x_down, GLdouble x_top,
                   GLdouble x_near, GLdouble x_far,
                   GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height
                   ){
  GLint n_left, n_down;
  GLsizei n_width, n_height;
  GLdouble nx_left, nx_right, nx_down, nx_top;

  switch (quad){
  case 0:
    glViewport(i_left,i_down,i_width,i_height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if(camera_current->projection_type==0){
      glFrustum(
        (double)x_left,(double)x_right,
        (double)x_down,(double)x_top,
        (double)x_near,(double)x_far);
    }
    else{
      glOrtho(
        (double)x_left,(double)x_right,
        (double)x_down,(double)x_top,
        (double)x_near,(double)x_far);
    }
    return 1;
  case 1:
    n_left = 2*i_left - s_left;
    n_down = 2*i_down - s_down;
    n_width = 2*i_width;
    n_height = 2*i_height;
    nx_left = x_left;
    nx_right = x_right;
    nx_down = x_down;
    nx_top = x_top;
    if(n_left+n_width<0||n_down+n_height<0)return 0;
    if(n_left>screenWidth||n_down>screenHeight)return 0;
    if(n_left<0){
      nx_left = x_left - n_left*(x_right-x_left)/n_width;
      n_width = n_left + n_width;
      n_left = 0;
    }
    if(n_left+n_width>screenWidth){
      nx_right = x_left + (screenWidth-n_left)*(x_right-x_left)/n_width;
      n_width = screenWidth-n_left;
    }
    if(n_down<0){
      nx_down = x_down - n_down*(x_top-x_down)/n_height;
      n_height = n_down + n_height;
      n_down = 0;
    }
    if(n_down+n_height>screenHeight){
      nx_top = x_down + (screenHeight-n_down)*(x_top-x_down)/n_height;
      n_height = screenHeight-n_down;
    }
    glViewport(n_left,n_down,n_width,n_height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    if(camera_current->projection_type==0){
      glFrustum(
        (double)nx_left,(double)nx_right,
        (double)nx_down,(double)nx_top,
        (double)x_near,(double)x_far);
    }
    else{
      glOrtho(
        (double)nx_left,(double)nx_right,
        (double)nx_down,(double)nx_top,
        (double)x_near,(double)x_far);
    }
    return 1;
  default:
    ASSERT(FFALSE);
    break;
  }
  return 1;
}

 /* ------------------------ BLOCK viewport ------------------------- */

void BLOCK_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height){
  float mesh_left;
  char slicelabel[255];
  float mesh_bot;
  int portview_left;
  float val_right,val_top;

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
  
    
  if(visBlocklabel==1&&nmeshes>1){
    int labellength;
    char meshlabel[255];

    sprintf(meshlabel,"mesh: %i",highlight_mesh+1);
    labellength=glutBitmapLength(large_font, (const unsigned char *)meshlabel);
    mesh_left=val_right-val_right*labellength/(float)dwinW;
    mesh_bot=val_top-val_top*large_font_height/(float)(dwinH-fontHoffset);
    outputText(mesh_left,mesh_bot, meshlabel);
  }
  if((showplot3d==1||visGrid==1)&&current_mesh->visx==1){
    {
      float plotval;
      int iplotval;
      char buff_label[128];

      iplotval=current_mesh->plotx;
      plotval=current_mesh->xplt_orig[iplotval];
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
    }
    if(visgridloc==1){
      outputText(mesh_left-0.5,0.6f, slicelabel);
    }
  }
  if((showplot3d==1||visGrid==1)&&current_mesh->visy==1){
    {
      float plotval;
      int iplotval;
      char buff_label[128];

      iplotval=current_mesh->ploty;
      plotval=current_mesh->yplt_orig[iplotval];
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
    }
    if(visgridloc==1)outputText(mesh_left-0.5,0.35f, slicelabel);
  }
  if((showplot3d==1||visGrid==1)&&current_mesh->visz==1){
    {
      float plotval;
      int iplotval;
      char buff_label[128];

      iplotval=current_mesh->plotz;
      plotval=current_mesh->zplt_orig[iplotval];
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
    }
    if(visgridloc==1)outputText(mesh_left-.5,0.1f, slicelabel);
  }
}

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
       (vis_slice_average==1&&show_slice_average&&(slice_average_flag==1||slice_turbprop_flag==1))
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
    if(show_slice_average==1&&vis_slice_average==1&&(slice_average_flag==1||slice_turbprop_flag==1)){
      sprintf(frameratelabel," AVG: %4.1f",slice_average_interval);
      outputText((float)(xtimeright+0.025),0.56, frameratelabel); // test print
    }
    if(hrrpuv_loaded==1&&show_hrrcutoff==1&&current_mesh!=NULL){
      char hrrcut_label[256];
      int ihrrcut;
      float xxl, xxr, yyl, yyu, ddx=0.03, ddy=0.2;

      ihrrcut = current_mesh->hrrpuv_cutoff;

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
    if(visAvailmemory==1){
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
  }
}

/* --------------------- COLOR BAR Viewport ------------------------- */

void COLORBAR_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height){
  GLint temp;
  float xnum;


  if(visColorLabels==1&&numColorbars!=0){
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

/* -------------------------- LOGO Viewport -------------------------- */

void LOGO_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height){
  if(use_nistlogo==1){   
      if(SUB_portortho(quad,
        titlesafe_offset,screenHeight-15-titlesafe_offset,75,100,
        0.0,1.4,0.0,1.75,
        s_left, s_down, s_width, s_height)==0)return;

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
//    outputLargeText(0.,0.,"NIST");
    {
      int n_loops=25,ipoint,ipoint2,iloop;
      float NIST[]={
      0.0,	0.0,                       
      0.0,	0.213947991,
      0.002354788,	0.226162333,
      0.010989011,	0.243498818,
      0.022762951,	0.253349094,
      0.030612245,	0.257289204,
      0.046310832,	0.261229314,
      0.062009419,	0.262017336,
      0.089481947,	0.256895193,
      0.105180534,	0.246256895,
      0.277472527,	0.064617809,
      0.285714286,	0.060283688,
      0.290816327,	0.064223798,
      0.293563579,	0.071315997,
      0.293563579,	0.262017336,
      0.348508634,	0.262017336,
      0.348508634,	0.052403467,
      0.340266876,	0.025216706,
      0.324175824,	0.008668243,
      0.30533752 ,	0.000394011,
      0.283751962,	0.0,
      0.266091052,	0.001182033,
      0.24686,	    0.010638,
      0.073783359,	0.201339638,
      0.065934066,	0.204491726,
      0.060047096,	0.201339638,
      0.058084772,	0.190307329,
      0.058084772,	0.0,
      0.390894819,	0.262017336,
      0.442700157,	0.262017336,
      0.442700157,	0.072892041,
      0.446624804,	0.06107171,
      0.458791209,	0.0571316,
      0.717817896,	0.0571316,
      0.733124019,	0.064223798,
      0.740580848,	0.079196217,
      0.737048666,	0.093380615,
      0.725274725,	0.102836879,
      0.710361068,	0.103624901,
      0.549058085,	0.103624901,
      0.520408163,	0.1107171,
      0.497645212,	0.124113475,
      0.478021978,	0.160756501,
      0.477237049,	0.198975571,
      0.491365777,	0.237194641,
      0.521193093,	0.260441292,
      0.54277865,	0.262017336,
      1.0,	0.262017336,
      1.0,	0.20606777,
      0.888147567,	0.20606777,
      0.888147567,	0.0,
      0.83477237,	0.0,
      0.83477237,	0.20606777,
      0.550627943,	0.20606777,
      0.537676609,	0.201339638,
      0.532574568,	0.189913318,
      0.533751962,	0.169818755,
      0.543956044,	0.159574468,
      0.553375196,	0.158786446,
      0.735086342,	0.158786446,
      0.78021978,	0.138691883,
      0.795918367,	0.100472813,
      0.791601256,	0.044523247,
      0.776295133,	0.017730496,
      0.760989011,	0.005516154,
      0.737048666,	0.0,
      0.454081633,	0.0,
      0.423469388,	0.005910165,
      0.407378336,	0.017336485,
      0.39599686,	0.038613081,
      0.391679749,	0.063829787,
      0.390894819,	0.095744681
      };
      int loops[25][9]={
      {8,0,1,2,3,4,25,26,27},
      {7,4,5,6,7,8,24,25},
      {4,8,9,23,24},
      {4,9,10,22,23},
      {7,10,11,18,19,20,21,22},
      {4,11,12,17,18},
      {4,12,13,16,17},
      {4,13,14,15,16},
      {4,28,29,31,69},
      {5,30,31,67,68,69},
      {4,31,32,66,67},
      {4,32,33,65,66},
      {5,33,34,63,64,65},
      {4,34,35,62,63},
      {4,35,36,61,62},
      {4,36,37,60,61},
      {4,37,38,59,60},
      {4,38,39,58,59},
      {5,39,40,41,57,58},
      {4,41,42,56,57},
      {4,42,43,55,56},
      {4,43,44,54,55},
      {4,44,45,53,54},
      {7,45,46,47,48,49,52,53},
      {4,52,49,50,51}
      };
      glColor3fv(foregroundcolor);
      for(iloop=0;iloop<n_loops;iloop++){
        glBegin(GL_POLYGON);
        for(ipoint2=1;ipoint2<=loops[iloop][0];ipoint2++){
          ipoint = loops[iloop][0]+1-ipoint2;
          glVertex2fv(NIST+2*loops[iloop][ipoint]);
        }
        glEnd();
      }
    }
  }
}

    /* -------------------------- TITLE Viewport -------------------------- */

void TITLE_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height){
  int left;
  float textdown;

  if(visTitle!=1)return;
  
  if(screenWidth<screenHeight){
    if(SUB_portortho(quad,
      fontWoffset+titlesafe_offset,(int)(screenHeight-1.1f*ntitles*dwinH/4.f-fontHoffset)-titlesafe_offset,screenWidth-dwinWW-fontWoffset-2*titlesafe_offset,(int)(ntitles*dwinH/4),
      0.,1.,0.,(double)(ntitles*ratio),
      s_left, s_down, s_width, s_height)==0)return;
    left=(float)75/(float)(screenWidth-dwinWW);
  }
  else{
    if(SUB_portortho(quad,
      fontWoffset+titlesafe_offset,(int)(screenHeight-1.1f*ntitles*dwinH/4.f-fontHoffset)-titlesafe_offset,screenWidth-dwinWW-fontWoffset-2*titlesafe_offset,(int)(ntitles*dwinH/4),
      0.,1.,0.,(double)(ntitles*ratio),
      s_left, s_down, s_width, s_height)==0)return;
    left=(float)75/(float)(screenWidth-dwinWW)*ratio;
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
    if(visColorLabels==1||(visBlocklabel==1&&nmeshes>1))right=screenWidth-dwinWW;
    if(visTitle==1)up=screenHeight-1.1*ntitles*dwinH/4.0;
  }

  if((visTimeLabels==1&&showtime==1)||(showtime==1&&(visFramerate==1||benchmark==1))||(visGrid==1&&visgridloc==1)
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
    if(plotstate==DYNAMIC_PLOTS&&selected_tour!=NULL&&selected_tour->path_timeslist!=NULL){
      if((viewtourfrompath==1&&selectedtour_index>=0)||keyframe_snap==1){
        touri = tourinfo + selectedtour_index;
        iframe = touri->path_timeslist[itimes];
        if(keyframe_snap==1&&selected_frame!=NULL){
          pj=&selected_frame->nodeval;
        }
        else{
          pj = touri->pathnodes + iframe;
        }

        camera_current->eye[0]=pj->eye[0];
        camera_current->eye[1]=pj->eye[1];
        camera_current->eye[2]=pj->eye[2];
        camera_current->angle_zx[1]=0.0;
        camera_current->angle_zx[0]=0.0;

      }
    }
  }

  if(plotstate==DYNAMIC_PLOTS&&select_avatar==1&&selected_avatar_tag>0&&view_from_selected_avatar==1){
    camera_current->eye[0]=selected_avatar_pos[0];
    camera_current->eye[1]=selected_avatar_pos[1];
    camera_current->eye[2]=selected_avatar_pos[2];
    camera_current->direction_angle=selected_avatar_angle;
    camera_current->view_angle=0.0;
    update_camera(camera_current);
    //camera_current->angle_zx[1]=0.0;
    //camera_current->angle_zx[0]=0.0;
  }

  eyexINI = camera_current->eye[0];
  eyeyINI = camera_current->eye[1];
  eyezINI = camera_current->eye[2];

  angleyzINI = camera_current->angle_zx[1];
  anglexyINI = camera_current->angle_zx[0];
  
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

    if(plotstate==DYNAMIC_PLOTS&&selected_tour!=NULL&&selected_tour->path_timeslist!=NULL){
      if((viewtourfrompath==1&&selectedtour_index>=0)||keyframe_snap==1){
        touri = tourinfo + selectedtour_index;
        iframe = touri->path_timeslist[itimes];
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


  widthdiv2 = fnear*tan(0.5*PI*aperture_temp/180.);
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
    float sn_direction_angle, cs_direction_angle;
    float sn_view_angle, cs_view_angle;
    float *uup;
    float *modelview_rotate;
    float cos_elevation_angle, sin_elevation_angle;
    float xcen, ycen, zcen;
    float posx, posy, posz;

    sn_view_angle=camera_current->sin_view_angle;
    cs_view_angle=camera_current->cos_view_angle;

    sn_direction_angle=camera_current->sin_direction_angle;
    cs_direction_angle=camera_current->cos_direction_angle;

    direction_angleINI = camera_current->direction_angle;

    modelview_rotate = camera_current->modelview;

    xcen = camera_current->xcen;
    ycen = camera_current->ycen;
    zcen = camera_current->zcen;

    cos_elevation_angle=camera_current->cos_elevation_angle;
    sin_elevation_angle=camera_current->sin_elevation_angle;

    sin_dv_sum = sn_direction_angle*cs_view_angle + cs_direction_angle*sn_view_angle;
    cos_dv_sum = cs_direction_angle*cs_view_angle - sn_direction_angle*sn_view_angle;


    posx = eyexINI+StereoCameraOffset*cos_dv_sum;
    posy = eyeyINI-StereoCameraOffset*sin_dv_sum;
    posz = eyezINI;

    viewx = posx + sin_dv_sum*cos_elevation_angle;
    viewy = posy + cos_dv_sum*cos_elevation_angle;
    viewz = posz + sin_elevation_angle;

    /* set view direction for virtual tour */
    {
      tourdata *touri;
      pathdata *pj;

      if(plotstate==DYNAMIC_PLOTS&&selected_tour!=NULL&&selected_tour->path_timeslist!=NULL){
        if((viewtourfrompath==1&&selectedtour_index>=0)||keyframe_snap==1){
          touri = tourinfo + selectedtour_index;
          iframe = touri->path_timeslist[itimes];
          if(keyframe_snap==1&&selected_frame!=NULL){
            pj=&selected_frame->nodeval;
          }
          else{
            pj = touri->pathnodes + iframe;
          }

          viewx = pj->oview[0]+StereoCameraOffset*cos_dv_sum;
          viewy = pj->oview[1]-StereoCameraOffset*sin_dv_sum;
          viewz = pj->oview[2];
          angleyzINI=0.0;
          anglexyINI=0.0;
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

    vecyz[0]=   cs_direction_angle; 
    vecyz[1] = -sn_direction_angle; 
    vecyz[2] = 0.0; 
    vecyz[3] = 1.0;

    glMultMatrixf(modelview_rotate);
    
    glTranslatef(xcen,ycen,zcen);
    if(eyeview==WORLD_CENTERED){
      glRotatef(angleyzINI,vecyz[0],vecyz[1],vecyz[2]);  /* rotate about the transformed x axis */
    }
    glRotatef(anglexyINI,0.0,0.0,1.0);                   /* rotate about z axis */
    glTranslatef(-xcen,-ycen,-zcen);

    glGetFloatv(GL_MODELVIEW_MATRIX,modelview_scratch);
    matmatmult(inverse_modelview_setup,modelview_scratch,modelview_current);

    get_world_eyepos(modelview_scratch, world_eyepos,scaled_eyepos);

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
    if(ngeominfo>0)Sort_Embedded_Geometry(modelview_scratch);
    if(showiso==1&&sort_iso_triangles==1&&niso_trans>0)Sort_Iso_Triangles(modelview_scratch);

    glScalef(mscale[0],mscale[1],mscale[2]);
    ExtractFrustum();
    set_cull_vis();
  }
}
