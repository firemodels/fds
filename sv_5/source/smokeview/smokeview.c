#include "options.h"
#ifdef pp_GPU
#include <GL/glew.h>
#endif
#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "contourdefs.h"
#include "isodefs.h"
#include "flowfiles.h"
#include "smokeviewapi.h"
#include "MALLOC.h"
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"
/* ksksksk */

#ifdef WIN32
#include <Commdlg.h>
#include <direct.h>
#endif

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
    break;
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
 /* ------------------------ SUB_portortho ------------------------- */
 
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
    if(projection_type==0){
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
    break;
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
    if(projection_type==0){
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
    break;
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

  mesh_left=0.9;
  if(fontindex==LARGE_FONT)mesh_left=0.7;
  if(visColorLabels==1){
    int portview_left;
    float val_right,val_top;

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

      sprintf(slicelabel,"mesh: %i",highlight_mesh+1);
      labellength=glutBitmapLength(large_font, slicelabel);
      mesh_left=val_right-val_right*labellength/(float)dwinW;
      mesh_bot=val_top-val_top*large_font_height/(float)(dwinH-fontHoffset);
      outputText(mesh_left,mesh_bot, slicelabel);
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
      if(visgridloc==1)outputText(mesh_left-0.5,0.6f, slicelabel);
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
}

/* ------------------------ TIME BAR Viewport ------------------------- */

void TIMEBAR_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height){
  int timebarheight;
#ifdef pp_memstatus
  unsigned int availmemory;
  char percen[]="%";
#endif


#ifdef pp_memstatus
  MEMSTATUS(0,&availmemory,NULL,NULL);
#endif

  if(
    (visTimeLabels==1&&showtime==1)||(showtime==1&&(visFramerate==1||benchmark==1||(vis_slice_average==1&&show_slice_average&&slice_average_flag==1)))
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
      s_left, s_down, s_width, s_height)==0)return;

      xtemp=1.0;
    }
     else{
    if(SUB_portortho(quad,
      fontWoffset+titlesafe_offset, fontHoffset+titlesafe_offset, screenWidth-dwinWW-2*fontWoffset-2*titlesafe_offset, timebarheight,
      0.0,ratio,0.,0.75,
      s_left, s_down, s_width, s_height)==0)return;

       xtemp=ratio;
     }

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

     xtimeleft=85.0f*xtemp/(screenWidth-dwinWW);
     xtimeright=xtimeleft+0.6*(xtemp-xtimeleft);

     if( visTimeLabels==1&&showtime==1){
      if(visTimelabel==1)outputText(0.0f,0.1f, timelabel);
#ifdef pp_HRR
      if(visFramelabel==1&&visHRRlabel==0)outputText(0.0f,0.4f, framelabel);
      if(visHRRlabel==1&&hrrinfo!=NULL)outputText(0.0f,0.4f, hrrinfo->hrrlabel);
#else
      if(visFramelabel==1)outputText(0.0f,0.4f, framelabel);
#endif
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
#ifdef pp_memstatus
    if(visAvailmemory==1){
      sprintf(frameratelabel," Mem Load:%i%s",availmemory,percen);
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
        screenWidth-dwinWW-fontWoffset-titlesafe_offset,
        temp+titlesafe_offset,
        dwinWW,
        screenHeight-temp-fontHoffset-2*titlesafe_offset,
        0.,(double)barright,-0.5,(double)(ratio*(nrgb+1)),
        s_left, s_down, s_width, s_height)==0)return;
    }
    else{
      barright=ratio*(xnum/3.0+0.1f);
      if(SUB_portortho(quad,
       screenWidth-dwinWW-fontWoffset-titlesafe_offset,temp+titlesafe_offset,dwinWW,screenHeight-temp-fontHoffset-2*titlesafe_offset,
       0.,(double)barright,-0.5,(double)(nrgb+1),
       s_left, s_down, s_width, s_height)==0)return;
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

  if(visTitle==1){
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
       if(trainer_mode==0){
         if(visFullTitle==1&&showplot3d==1){
           outputText(left,(ititle-1)*ratio+textdown, FULLTITLE);
         }
         else{
           outputText(left,(ititle-1)*ratio+textdown, TITLE);
         }
       }
       if(trainer_mode==1){
         outputText(left,(ititle-1)*ratio+textdown, TRAINERTITLE);
       }
     }
  }

} 

/* ----------------------- 3D scene Viewport ----------------------------- */
void Scene_viewport(int quad, int view_mode, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height){

  float up, down, left, right;
  float dh;
  float fnear, ffar, fleft, fright, fup, fdown, fzero;
  float StereoCameraOffset,FrustumAsymmetry;
  float aperture_temp;

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
    if(visColorLabels==1)right=screenWidth-dwinWW;
    if(visTitle==1)up=screenHeight-1.1*ntitles*dwinH/4.0;
  }

  if((visTimeLabels==1&&showtime==1)||(showtime==1&&(visFramerate==1||benchmark==1))
#ifdef pp_memstatus
      ||visAvailmemory==1
#endif
    ){
    down=0.75*dwinH;
  }
  down += fontHoffset;
  dh = (up-down)*2*fontWoffset/(right-left);
  
  aspect=(float)(right-left)/(float)(up-down);

  /* set view position for virtual tour */

  {
    tourdata *touri;
    pathdata *pj;
    if(plotstate==DYNAMIC_PLOTS&&selected_tour!=NULL&&selected_tour->path_timeslist!=NULL){
      if((viewtourfrompath==1&&selectedtour_index>=0)||keyframe_snap==1){
        touri = tourinfo + selectedtour_index;
        iframe = touri->path_timeslist[itime];
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

  eyexINI = camera_current->eye[0];
  eyeyINI = camera_current->eye[1];
  eyezINI = camera_current->eye[2];

  angleyzINI = camera_current->angle_zx[1];
  anglexyINI = camera_current->angle_zx[0];
  
  fnear =  - eyeyINI-1.0;
  if(fnear<nearclip)fnear=nearclip;
  ffar = fnear + farclip;
  fzero = camera_current->ycen-eyeyINI;
  if(fzero<fnear)fzero=fnear+0.5;
  
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
        iframe = touri->path_timeslist[itime];
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

  fup = fnear*tan(PI*aperture_temp/360.);
  if(showstereo==2)fup*=2.0;
  fdown = -fup;
  fright = fup*aspect;
  fleft = -fright;
  
  if(showstereo==0||view_mode==VIEW_CENTER){
    StereoCameraOffset=0.0;
    FrustumAsymmetry=0.0;
  }
  else if(showstereo!=0&&(view_mode==VIEW_LEFT||view_mode==VIEW_RIGHT)){
    StereoCameraOffset = (fright-fleft)*0.035*eoffset*fzero/fnear;
    if(view_mode==VIEW_LEFT)StereoCameraOffset = -StereoCameraOffset;
    FrustumAsymmetry=-StereoCameraOffset*pbalance*fnear/fzero;
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

    viewx = eyexINI + sin_dv_sum*cos_elevation_angle;
    viewy = eyeyINI + cos_dv_sum*cos_elevation_angle;
    viewz = eyezINI + sin_elevation_angle;

    /* set view direction for virtual tour */
    {
      tourdata *touri;
      pathdata *pj;

      if(plotstate==DYNAMIC_PLOTS&&selected_tour!=NULL&&selected_tour->path_timeslist!=NULL){
        if((viewtourfrompath==1&&selectedtour_index>=0)||keyframe_snap==1){
          touri = tourinfo + selectedtour_index;
          iframe = touri->path_timeslist[itime];
          if(keyframe_snap==1&&selected_frame!=NULL){
            pj=&selected_frame->nodeval;
          }
          else{
            pj = touri->pathnodes + iframe;
          }

          viewx = pj->oview[0];
          viewy = pj->oview[1];
          viewz = pj->oview[2];
          angleyzINI=0.0;
          anglexyINI=0.0;
        }
      }
    }

    uup = camera_current->up;
#ifdef pp_SPOTLIGHT
      if(eyeview==1&&spotlight==1){
        float spot_position[3], spot_direction[3];

        spot_position[0]=eyexINI;
        spot_position[1]=eyeyINI;
        spot_position[2]=eyezINI;
        spot_direction[0]=viewx;
        spot_direction[1]=viewy;
        spot_direction[2]=viewz;
        updateSpotLight(0,spot_position, spot_direction);
      }
#endif

    gluLookAt(
      (double)(eyexINI+StereoCameraOffset*cos_dv_sum),
      (double)(eyeyINI-StereoCameraOffset*sin_dv_sum),
      (double)eyezINI,
      (double)viewx,
      (double)viewy,
      (double)viewz,
      uup[0],
      uup[1],
      uup[2]);

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

    if(nsmoke3d>0)getsmokedir(modelview_scratch);
    if(nface_transparent>0)sort_transparent_faces(modelview_scratch);

    // calculate transparent distances and sort

    glScalef(mscale[0],mscale[1],mscale[2]);
    ExtractFrustum();
  }
}
  
/* ----------------------- setClipPlanes ----------------------------- */

void setClipPlanes(void){
  static GLdouble clipplane_x[4], clipplane_y[4], clipplane_z[4];
  static GLdouble clipplane_X[4], clipplane_Y[4], clipplane_Z[4];

  if(xyz_clipplane==1){
    if(clip_x==1){
      clipplane_x[0]=1.0;
      clipplane_x[1]=0.0;
      clipplane_x[2]=0.0;
      clipplane_x[3]=-(clip_x_val-xbar0)/xyzmaxdiff;
      glClipPlane(GL_CLIP_PLANE0,clipplane_x);
      glEnable(GL_CLIP_PLANE0);
    }
    else{
      glDisable(GL_CLIP_PLANE0);
    }

    if(clip_X==1){
      clipplane_X[0]=-1.0;
      clipplane_X[1]=0.0;
      clipplane_X[2]=0.0;
      clipplane_X[3]=(clip_X_val-xbar0)/xyzmaxdiff;
      glClipPlane(GL_CLIP_PLANE3,clipplane_X);
      glEnable(GL_CLIP_PLANE3);
    }
    else{
      glDisable(GL_CLIP_PLANE3);
    }

    if(clip_y==1){
      clipplane_y[0]=0.0;
      clipplane_y[1]=1.0;
      clipplane_y[2]=0.0;
      clipplane_y[3]=-(clip_y_val-ybar0)/xyzmaxdiff;
      glClipPlane(GL_CLIP_PLANE1,clipplane_y);
      glEnable(GL_CLIP_PLANE1);
    }
    else{
      glDisable(GL_CLIP_PLANE1);
    }

    if(clip_Y==1){
      clipplane_Y[0]=0.0;
      clipplane_Y[1]=-1.0;
      clipplane_Y[2]=0.0;
      clipplane_Y[3]=(clip_Y_val-ybar0)/xyzmaxdiff;
      glClipPlane(GL_CLIP_PLANE4,clipplane_Y);
      glEnable(GL_CLIP_PLANE4);
    }
    else{
      glDisable(GL_CLIP_PLANE4);
    }

    if(clip_z==1){
      clipplane_z[0]=0.0;
      clipplane_z[1]=0.0;
      clipplane_z[2]=1.0;
      clipplane_z[3]=-(clip_z_val-zbar0)/xyzmaxdiff;
      glClipPlane(GL_CLIP_PLANE2,clipplane_z);
      glEnable(GL_CLIP_PLANE2);
    }
    else{
      glDisable(GL_CLIP_PLANE2);
    }

    if(clip_Z==1){
      clipplane_Z[0]=0.0;
      clipplane_Z[1]=0.0;
      clipplane_Z[2]=-1.0;
      clipplane_Z[3]=(clip_Z_val-zbar0)/xyzmaxdiff;
      glClipPlane(GL_CLIP_PLANE5,clipplane_Z);
      glEnable(GL_CLIP_PLANE5);
    }
    else{
      glDisable(GL_CLIP_PLANE5);
    }
  }
  else{
    glDisable(GL_CLIP_PLANE0);
    glDisable(GL_CLIP_PLANE1);
    glDisable(GL_CLIP_PLANE2);
    glDisable(GL_CLIP_PLANE3);
    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);
  }
}


/* ------------------ ShowScene ------------------------ */

void ShowScene(int mode, int view_mode, int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height){
  int i;

  CheckMemory;

  if(xyz_clipplane==1){
    if(clip_x==1)glDisable(GL_CLIP_PLANE0);
    if(clip_y==1)glDisable(GL_CLIP_PLANE1);
    if(clip_z==1)glDisable(GL_CLIP_PLANE2);
    if(clip_X==1)glDisable(GL_CLIP_PLANE3);
    if(clip_Y==1)glDisable(GL_CLIP_PLANE4);
    if(clip_Z==1)glDisable(GL_CLIP_PLANE5);
  }

/* ++++++++++++++++++++++++ update variables as needed +++++++++++++++++++++++++ */

#ifdef pp_SVNET
    parse_svcom();
#endif
  if(loadfiles_at_startup&&update_load_startup==1){
    load_startup_smoke();
  }
  if(updategluiview==1){
    reset_glui_view(startup_view_ini);
    updategluiview=0;
  }
  if(unitclasses_ini!=NULL){
    nunitclasses=nunitclasses_ini;
    unitclasses=unitclasses_ini;
  }
  else{
    nunitclasses=nunitclasses_default;
    unitclasses=unitclasses_default;
  }
  //debug
  if(menusmooth==1&&smoothing_blocks==0&&updatesmoothblocks==1){
    smooth_blockages();
  }
  if(update_tourlist==1){
    Update_Tourlist();
  }
  if(camera_current->dirty==1){
    update_camera(camera_current);
  }
  if(update_selectedtour_index==1){
    update_tourindex();
  }
  if(trainer_mode==1&&fontindex!=LARGE_FONT)FontMenu(1);
  if(updateindexcolors==1){
    UpdateIndexColors();
  }
  if(force_isometric==1){
    force_isometric=0;
    projection_type=0;
    ZoomMenu(-2);
  }


  updateShow();
  if(times!=NULL&&updateUpdateFrameRateMenu==1)sv_FrameRateMenu(frameratevalue);
  if(updatefaces==1)update_faces();
  if(updatefacelists==1)update_facelists();
  if(showstereo!=2)ClearBuffers(mode);

/* ++++++++++++++++++++++++ draw viewports +++++++++++++++++++++++++ */
  if(mode==RENDER){
    /*
    BLOCK_viewport(0,0,0,screenWidth,screenHeight);
    TIMEBAR_viewport(0,0,0,screenWidth,screenHeight);
    COLORBAR_viewport(0,0,0,screenWidth,screenHeight);
    LOGO_viewport(0,0,0,screenWidth,screenHeight);
    TITLE_viewport(0,0,0,screenWidth,screenHeight);
    Scene_viewport(0,view_mode,0,0,screenWidth,screenHeight);
    */
       BLOCK_viewport(quad,          s_left,s_down,s_width,s_height);
     TIMEBAR_viewport(quad,          s_left,s_down,s_width,s_height);
    COLORBAR_viewport(quad,          s_left,s_down,s_width,s_height);
        LOGO_viewport(quad,          s_left,s_down,s_width,s_height);
       TITLE_viewport(quad,          s_left,s_down,s_width,s_height);
       Scene_viewport(quad,view_mode,s_left,s_down,s_width,s_height);
  }

  if(eyeview==1&&nskyboxinfo>0)draw_skybox();

  if(UpdateLIGHTS==1)updateLights(0);

  if(mode!=RENDER||viscolorbarpath!=1)setClipPlanes();
  if(mode==RENDER){
    if(viscolorbarpath==1){
      if(cb_hidesv==1){
        setColorbarClipPlanes(1);
      }
      else{
        setColorbarClipPlanes(0);
      }
    }
    glPointSize((float)1.0);

#ifdef pp_WUI

    /* ++++++++++++++++++++++++ draw trees +++++++++++++++++++++++++ */

    if(ntreeinfo>0){
      drawtrees();
    }
#endif

/* ++++++++++++++++++++++++ draw smoke +++++++++++++++++++++++++ */

    if(showsmoke==1){
      particle *parti;

      if(staticframe0==0||iframe!=0){
        for(i=0;i<npartinfo;i++){
          parti = partinfo + i;
          if(parti->loaded==0||parti->display==0)continue;
          if(parti->evac==1){
            drawEvac(parti);
          }
          else{
            drawPart(parti);
          }
        }
      }
      if(visStaticSmoke==1&&staticframe0==1){
        for(i=0;i<npartinfo;i++){
          parti = partinfo + i;
          if(parti->loaded==0||parti->display==0)continue;
          drawStaticPart(parti);
        }
      }
    }

/* ++++++++++++++++++++++++ draw evacuation +++++++++++++++++++++++++ */

    if(showevac==1){
        for(i=0;i<npartinfo;i++){
          particle *parti;

          parti = partinfo + i;
          if(parti->loaded==0||parti->display==0||parti->evac==0)continue;
          drawEvac(parti);
        }
    }

/* ++++++++++++++++++++++++ draw targets +++++++++++++++++++++++++ */

    if(showtarget==1){
      drawTargets();
    }

/* ++++++++++++++++++++++++ draw sensors/sprinklers/heat detectors +++++++++++++++++++++++++ */

    draw_devices();
    sniffErrors("after draw_devices");

    if(visaxislabels==1||showedit==1){
      outputAxisLabels();
    }


 /* ++++++++++++++++++++++++ draw ticks +++++++++++++++++++++++++ */

    if(visTicks==1&&nticks>0)drawticks();

    /* draw the box framing the simulation (corners at (0,0,0) (xbar,ybar,zbar) */


/* ++++++++++++++++++++++++ draw simulation frame (corners at (0,0,0) and (xbar,ybar,zbar) +++++++++++++++++++++++++ */

     if(isZoneFireModel==0&&visFrame==1&&highlight_flag==2){
       drawoutlines();
     }
     sniffErrors("after drawbox");


/* ++++++++++++++++++++++++ draw mesh +++++++++++++++++++++++++ */

     if(setPDIM==1){
       if(visGrid==1){
         int igrid;
         mesh *meshi;

         for(igrid=0;igrid<selected_case->nmeshes;igrid++){
           meshi=selected_case->meshinfo+igrid;
           drawgrid(meshi);
         }
       }
     }
   } /* end of if(mode==RENDER) code segment */


/* ++++++++++++++++++++++++ draw tours +++++++++++++++++++++++++ */

  if(mode==SELECT&&edittour==1&&ntours>0){
    drawselecttours();
    sniffErrors("after drawselecttours");
    return;
  }

/* ++++++++++++++++++++++++ draw "fancy" colorbar +++++++++++++++++++++++++ */

#ifdef pp_COLOR
  if(viscolorbarpath==1){
    if(cb_hidesv==1){
      setColorbarClipPlanes(0);
    }
    drawcolorbarpath();
    if(cb_hidesv==1){
      setColorbarClipPlanes(1);
    }
  }

#endif
/* ++++++++++++++++++++++++ draw blockages +++++++++++++++++++++++++ */

  drawBlockages(mode);

#ifdef pp_WUI
/* ++++++++++++++++++++++++ draw terrain +++++++++++++++++++++++++ */

  if(visTerrain==1){
    int i;

    for(i=0;i<nterraininfo;i++){
      terraindata *terri;

      terri = terraininfo + i;
      if(terri->loaded==1)drawterrain(terri);
    }
  }
#endif

/* ++++++++++++++++++++++++ draw boundary files +++++++++++++++++++++++++ */

  if(showpatch==1){
    patch *patchi;
    mesh *meshi;

    for(i=0;i<selected_case->nmeshes;i++){
      meshi=selected_case->meshinfo+i;
      if(meshi->npatches>0){
        int filenum;

        filenum=meshi->patchfilenum;
        if(filenum!=-1){
          patchi = patchinfo + filenum;
          if(patchi->loaded==0||patchi->display==0||patchi->type!=ipatchtype)continue;
          if(usetexturebar==1){
            if(vis_ignited==1&&do_ignited==1){
              drawpatch_texture_char(meshi);
            }
            else{
              drawpatch_texture(meshi);
            }
          }
          else{
            drawpatch(meshi);
          }
          if(vis_ignited==1&&vis_onlyignited==1&&do_ignited==1)drawonlyignited(meshi);
        }
      }
    }
  }

/* ++++++++++++++++++++++++ draw labels +++++++++++++++++++++++++ */

  if(visLabels==1){
    drawLabels();
  }

/* ++++++++++++++++++++++++ draw animated isosurfaces +++++++++++++++++++++++++ */

#ifdef pp_SPHERE
    if(isoinfo!=NULL)drawspherepoints(sphereinfo);
#endif
  if(showiso==1){
    iso *isoi;
    mesh *meshi;

    for(i=0;i<selected_case->nmeshes;i++){
      meshi=selected_case->meshinfo+i;
      if(meshi->isotimes==NULL||meshi->isofilenum<0)continue;
      isoi = isoinfo + meshi->isofilenum;
      if(isoi->loaded==0||isoi->display==0||isoi->type!=iisotype)continue;
      drawiso(meshi,SOLID);
    }

    /*  nothing transparent should be drawn before this portion of the code
    (ie draw all opaque opbjects first then 
        draw transparent objects
        */

    for(i=0;i<selected_case->nmeshes;i++){
      meshi=selected_case->meshinfo+i;
      if(meshi->isotimes==NULL||meshi->isofilenum<0)continue;
      isoi = isoinfo + meshi->isofilenum;
      if(isoi->loaded==0||isoi->display==0||isoi->type!=iisotype)continue;
      drawiso(meshi,SV_TRANSPARENT);
    }
  }

/* ++++++++++++++++++++++++ draw transparent faces +++++++++++++++++++++++++ */

    draw_transparent_faces();

/* ++++++++++++++++++++++++ draw 3D smoke +++++++++++++++++++++++++ */

  if(show3dsmoke==1){
    {
      smoke3d *smoke3di;

      CheckMemory;
#ifdef pp_GPU
      if(usegpu==1)useSmokeShaders();
#endif
      for(i=0;i<nsmoke3d;i++){
        smoke3di = smoke3dinfo + i;
        if(smoke3di->loaded==0||smoke3di->display==0)continue;
        if(smoke3di->d_display==0)continue;

#ifdef pp_GPU
        if(usegpu==1){
      //    getDepthTexture();
          drawsmoke3dGPU(smoke3di);
        }
        else{
          drawsmoke3d(smoke3di);
        }
#else
        drawsmoke3d(smoke3di);
#endif
      }
#ifdef pp_GPU
      if(usegpu==1)useOpenGLShaders();
#endif

    }
    sniffErrors("after drawsmoke");
  }

/* ++++++++++++++++++++++++ draw zone fire modeling info +++++++++++++++++++++++++ */

  if(nrooms>0){
    drawroomgeom();
    sniffErrors("after drawroomgeom");
  }
  if(showzone==1){
    drawroomdata();
    sniffErrors("after drawroomdata");
    if(ReadZoneFile==1&&nzvents>0){
      drawventdata();
      sniffErrors("after drawventdata");
    }
  }

/* ++++++++++++++++++++++++ draw slice files +++++++++++++++++++++++++ */
  
  if(showslice==1){
    int ii;

    for(ii=0;ii<nslice_loaded;ii++){
      slice *sd;

      i=slice_loaded_list[ii];
      sd = sliceinfo + i;
      if(sd->display==0||sd->type!=islicetype)continue;
      if(sd->slicetimes[0]>times[itime])continue;
      if(sd->compression_type==1||sd->compression_type==2){
        uncompress_slicedataframe(sd,sd->islice);
        sd->slicepoint=sd->slicecomplevel;
      }
      else{
        sd->slicepoint = sd->slicelevel + sd->islice*sd->nsliceii;
      }
      sd->slicedata=NULL;
      if(sd->compression_type==0){
        ASSERT(ValidPointer(sd->qslicedata,sizeof(float)*sd->nslicetotal));
      }

      if(sd->qslicedata!=NULL)sd->slicedata = sd->qslicedata + sd->islice*sd->nsliceii;
      if(usetexturebar==1){
        if(sd->volslice==1){
          drawvolslice_texture(sd);
       }
        else{
          drawslice_texture(sd);
        }
      }
      else{
        if(sd->volslice==1){
          drawvolslice(sd);
        }
        else{
          drawslice(sd);
        }
      }
    }
  } 
  sniffErrors("after drawslice");
//  draw_demo(20,20);
//  draw_demo2(1);

/* ++++++++++++++++++++++++ draw vector slice files +++++++++++++++++++++++++ */

  if(showvslice==1){
    for(i=0;i<nvslice;i++){
      vslice *vd;
      slice *u, *v, *w, *val;

      vd = vsliceinfo + i;
      if(vd->loaded==0||vd->display==0||sliceinfo[vd->ival].type!=islicetype)continue;
      val = vd->val;
      if(val==NULL)continue;
      u = vd->u;
      v = vd->v;
      w = vd->w;
      if(u==NULL&&v==NULL&&w==NULL)continue;
      if(sliceinfo[vd->ival].slicetimes[0]>times[itime])continue;
#define VAL val
      if(VAL->compression_type==1){
        uncompress_slicedataframe(VAL,VAL->islice);
        VAL->slicepoint=VAL->slicecomplevel;
      }
      else{
        if(VAL!=NULL)VAL->slicepoint = VAL->slicelevel + VAL->islice*VAL->nsliceii;
      }
#undef VAL
#define VAL u
      if(VAL!=NULL){
        if(VAL->compression_type==1){
          uncompress_slicedataframe(VAL,VAL->islice);
          VAL->slicepoint=VAL->slicecomplevel;
        }
        else{
          if(VAL!=NULL)VAL->slicepoint = VAL->slicelevel + VAL->islice*VAL->nsliceii;
        }
      }
#undef VAL
#define VAL v
      if(VAL!=NULL){
        if(VAL->compression_type==1){
          uncompress_slicedataframe(VAL,VAL->islice);
          VAL->slicepoint=VAL->slicecomplevel;
        }
        else{
          if(VAL!=NULL)VAL->slicepoint = VAL->slicelevel + VAL->islice*VAL->nsliceii;
        }
      }
#undef VAL
#define VAL w
      if(VAL!=NULL){
        if(VAL->compression_type==1){
          uncompress_slicedataframe(VAL,VAL->islice);
          VAL->slicepoint=VAL->slicecomplevel;
        }
        else{
          if(VAL!=NULL)VAL->slicepoint = VAL->slicelevel + VAL->islice*VAL->nsliceii;
        }
      }
      if(u!=NULL&&u->compression_type==0){
        u->qslice = u->qslicedata + u->islice*u->nsliceii;
      }
      if(v!=NULL&&v->compression_type==0){
        v->qslice = v->qslicedata + v->islice*v->nsliceii;
      }
      if(w!=NULL&&w->compression_type==0){
        w->qslice = w->qslicedata + w->islice*w->nsliceii;
      }
      if(vd->volslice==1){
        drawvvolslice(vd);
      }
      else{
        drawvslice(vd);
      }
    }
  }
  sniffErrors("after drawvslice");

/* ++++++++++++++++++++++++ draw tours +++++++++++++++++++++++++ */

  if(showtour==1){
    drawtours();
    sniffErrors("after drawTours");
  }


/* ++++++++++++++++++++++++ draw plot3d files +++++++++++++++++++++++++ */

  if(showplot3d==1){
    mesh *meshi;

    for(i=0;i<selected_case->nmeshes;i++){
      meshi=selected_case->meshinfo+i;
      if(meshi->plot3dfilenum==-1)continue;
      if(plot3dinfo[meshi->plot3dfilenum].display==0)continue;
      if(usetexturebar==1){
        drawplot3d_texture(meshi);
      }
      else{
        drawplot3d(meshi);
      }
    }
  }
  sniffErrors("after drawplot3d");

  /* ++++++++++++++++++++++++ draw cross hairs +++++++++++++++++++++++++ */

#ifdef _DEBUG
    if(eyeview==EYE_CENTERED)drawMovedir();
#endif

/* ++++++++++++++++++++++++ render scene +++++++++++++++++++++++++ */

#ifndef pp_nolibs
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
#endif
/* ++++++++++++++++++++++++ draw "fancy" colorbar +++++++++++++++++++++++++ */

#ifdef pp_COLOR
  if(viscolorbarpath==1){
    if(cb_hidesv==1){
      setColorbarClipPlanes(0);
    }
  }
#endif

  sniffErrors("end of loop");

}

/* ------------------ update_rotation_index ------------------------ */

void update_rotation_index(int val){
  mesh *meshi;
  int i;
  float *modelview_rotate;
  float *angle_zx;
  int *rotation_index;

  rotation_index = &camera_current->rotation_index;

  *rotation_index=val;
  if(*rotation_index==rotation_index_OLD)return;
  if(*rotation_index>=0&&*rotation_index<nmeshes){
    meshi = selected_case->meshinfo + *rotation_index;
    camera_current->xcen=meshi->xcen;
    camera_current->ycen=meshi->ycen;
    camera_current->zcen=meshi->zcen;
  }
  else{
    camera_current->xcen=xcenGLOBAL;
    camera_current->ycen=ycenGLOBAL;
    camera_current->zcen=zcenGLOBAL;
  }
  rotation_index_OLD=*rotation_index;
  modelview_rotate = camera_current->modelview;
  for(i=0;i<16;i++){
    modelview_rotate[i]=modelview_current[i];
  }

  angle_zx = camera_current->angle_zx;

  angle_zx[0]=0.; 
  angle_zx[1]=0.; 

  camera_current->direction_angle=0.0;
  camera_current->cos_direction_angle = 1.0;
  camera_current->sin_direction_angle = 0.0;

  camera_current->view_angle=0.0;
  camera_current->cos_view_angle = 1.0;
  camera_current->sin_view_angle = 0.0;

  update_meshlist1(val);

  glutPostRedisplay();

}

/* ------------------ matmatmult ------------------------ */

void matmatmult(float *m1, float *m2, float *m3){
  int i, j, k;
  int ij;

  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      ij = i+4*j;
      m3[ij]=0.0;
      for(k=0;k<4;k++){
        m3[ij]+=m1[i+4*k]*m2[k+4*j];
      }
    }
  }
}

/* ------------------ getinverse ------------------------ */

void getinverse(float *m, float *mi){
  int i,j;
  float *v,*vi;

  /*

  assume m is a 4x4 matrix parttioned as

  q00 q01 q02 v0
  q10 q11 q12 v1
  q20 q21 q22 v2
    0   0   0  a

  where v=(vi) and Q=(qij) is orthogonal ( Q*transpose(Q) = I )

  then inverse(m) =     transpose(Q)   -transpose(Q)*v/a
                            0                 1/a       

  note:  m_ij = m[i+4*j]
  */

  v=m+12;   /* fourth column of m */               
  vi=mi+12; /* fourth column of inverse(m) */
  for(i=0;i<3;i++){  /* compute transpose */
    for(j=0;j<3;j++){
      mi[i+4*j]=m[j+4*i];
    }
    mi[3+4*j]=0.0;
  }
  vi[3]=1.0/v[3];
  vi[0]=-(mi[0+4*0]*v[0]+mi[0+4*1]*v[1]+mi[0+4*2]*v[2])*vi[3];
  vi[1]=-(mi[1+4*0]*v[0]+mi[1+4*1]*v[1]+mi[1+4*2]*v[2])*vi[3];
  vi[2]=-(mi[2+4*0]*v[0]+mi[2+4*1]*v[1]+mi[2+4*2]*v[2])*vi[3];
}

/* ------------------ snifferrors ------------------------ */

void _sniffErrors(char *whereat){
  int error;
  while((error=glGetError())!=GL_NO_ERROR){
    fprintf(stderr,"GL Error:%s where:%s  %i\n",
      gluErrorString((unsigned int)error),whereat,snifferrornumber++);
  }
}

#ifdef pp_SPOTLIGHT
void updateSpotLight(int init, float *spot_position, float *spot_direction){
  GLfloat ambientlight2[4], diffuselight2[4];
  int i;

  return;
  if(init==1){
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);

    for(i=0;i<3;i++){
      ambientlight2[i]=ambientlight[i]/2.0;
      diffuselight2[i]=diffuselight[i]/2.0;
    }
    ambientlight2[3]=1.0;
    diffuselight2[3]=1.0;
    glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuselight2);
    glLightfv(GL_LIGHT0,GL_AMBIENT,ambientlight2);

    glLightfv(GL_LIGHT1,GL_DIFFUSE,diffuselight2);
    glLightfv(GL_LIGHT1,GL_AMBIENT,ambientlight2);
    glLightfv(GL_LIGHT1,GL_POSITION,light_position1);
  }
  glLightfv(GL_LIGHT0,GL_POSITION,spot_position);
  glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, spot_direction);
  glLightf(GL_LIGHT0, GL_SPOT_CUTOFF,10.0);  //xxx test

}
#endif

/* ------------------ updateLights ------------------------ */

void updateLights(int pos){
  GLfloat ambientlight2[4], diffuselight2[4];
  int i;

  if(visLIGHT0==1&&visLIGHT1==1){
    for(i=0;i<3;i++){
      ambientlight2[i]=ambientlight[i]/2.0;
      diffuselight2[i]=diffuselight[i]/2.0;
    }
    ambientlight2[3]=1.0;
    diffuselight2[3]=1.0;
    glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuselight2);
    glLightfv(GL_LIGHT0,GL_AMBIENT,ambientlight2);
    if(pos==1)glLightfv(GL_LIGHT0,GL_POSITION,light_position0);
    glEnable(GL_LIGHT0);

    glLightfv(GL_LIGHT1,GL_DIFFUSE,diffuselight2);
    glLightfv(GL_LIGHT1,GL_AMBIENT,ambientlight2);
    if(pos==1)glLightfv(GL_LIGHT1,GL_POSITION,light_position1);
    glEnable(GL_LIGHT1);

  }
  if(visLIGHT0==1&&visLIGHT1==0){
    diffuselight[3]=1.0;
    ambientlight[3]=1.0;
    glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuselight);
    glLightfv(GL_LIGHT0,GL_AMBIENT,ambientlight);
    if(pos==1)glLightfv(GL_LIGHT0,GL_POSITION,light_position0);
    glEnable(GL_LIGHT0);

    glDisable(GL_LIGHT1);

  }
  if(visLIGHT0==0&&visLIGHT1==1){
    diffuselight[3]=1.0;
    ambientlight[3]=1.0;
    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuselight);
    glLightfv(GL_LIGHT0,GL_AMBIENT,ambientlight);
    if(pos==1)glLightfv(GL_LIGHT0,GL_POSITION,light_position1);
    glDisable(GL_LIGHT1);
  }
  if(visLIGHT0==0&&visLIGHT1==0){
    glDisable(GL_LIGHT0);
    glDisable(GL_LIGHT1);
  }


  UpdateLIGHTS=0;

}

/* ------------------ updateShow ------------------------ */

void updateShow(void){
  int i,evacflag,sliceflag,vsliceflag,partflag,patchflag,isoflag,smoke3dflag;
  int ii;
  slice *sd;
  vslice *vd;
  iso *isoi;
  mesh *meshi;
  patch *patchi;
  particle *parti;
  showtime=0; showtime2=0; showplot3d=0; showpatch=0; 
  showslice=0; showvslice=0; showsmoke=0; showzone=0; showiso=0;
  showevac=0;
  showtarget=0;
  showtitle1=0; showtitle2=0;
  show3dsmoke=0;
  smoke3dflag=0;
  showtour=0;
#ifdef pp_WUI
  showterrain=0;
#endif
  ntitles=0;
  if(visTitle0==1)ntitles++;
  if(strlen(TITLE1)!=0&&visTitle1==1){ntitles++;showtitle1=1;}
  if(strlen(TITLE2)!=0&&visTitle2==1){ntitles++;showtitle2=1;}
  visTitle=0;
  if(visTitle0==1||showtitle1==1||showtitle2==1)visTitle=1;
  visTimeSmoke=1; visTimeSlice=1; visTimePatch=1; visTimeZone=1; visTimeIso=1;

  RenderTime=0;
  if(times!=NULL){
    if(settmin_p==1&&times[itime]<tmin_p)visTimeSmoke=0;
    if(settmax_p==1&&times[itime]>tmax_p)visTimeSmoke=0;

    if(settmin_s==1&&times[itime]<tmin_s)visTimeSlice=0;
    if(settmax_s==1&&times[itime]>tmax_s)visTimeSlice=0;

    if(settmin_i==1&&times[itime]<tmin_i)visTimeIso=0;
    if(settmax_i==1&&times[itime]>tmax_i)visTimeIso=0;

    if(settmin_b==1&&times[itime]<tmin_b)visTimePatch=0;
    if(settmax_b==1&&times[itime]>tmax_b)visTimePatch=0;

    if(settmin_z==1&&times[itime]<tmin_z)visTimeZone=0;
    if(settmax_z==1&&times[itime]>tmax_z)visTimeZone=0;

  }

  {
    tourdata *touri;

    if(ntours>0){
      for(i=0;i<ntours;i++){
        touri = tourinfo + i;
        if(touri->display==1){
          showtour=1;
          break;
        }
      }
    }
  }
#ifdef pp_WUI
  if(visTerrain==1){
    for(i=0;i<nterraininfo;i++){
      terraindata *terri;

      terri = terraininfo + i;
      if(terri->loaded==1){
        showterrain=1;
        break;
      }
    }
  }
#endif
  {
    smoke3d *smoke3di;

    for(ii=0;ii<nsmoke3d;ii++){
      smoke3di = smoke3dinfo + ii;
      if(smoke3di->loaded==1&&smoke3di->display==1){
        smoke3dflag=1;
        break;
      }
    }
  }
  sliceflag=0;
  if(visTimeSlice==1){
    for(ii=0;ii<nslice_loaded;ii++){
      i=slice_loaded_list[ii];
      sd = sliceinfo+i;
      if(sd->display==0||sd->type!=islicetype)continue;
      if(sd->nsteps>0){
        sliceflag=1;
        break;
      }
    }
  }
  isoflag=0;
  if(visTimeIso==1){
    for(i=0;i<niso;i++){
      isoi = isoinfo+i;
      if(isoi->loaded==0)continue;
 //xxx     if(isoi->display==0||isoi->type!=iisotype)continue;
 //     isoflag=1;
 //     break;
      if(isoi->display==1&&isoi->type==iisotype){
        isoflag=1;
        break;
      }
    }
  }
  vsliceflag=0;
  if(visTimeSlice==1){
    for(i=0;i<nvslice;i++){
      vd = vsliceinfo+i;
      if(vd->loaded==0||vd->display==0)continue;
      if(sliceinfo[vd->ival].type!=islicetype)continue;
      vsliceflag=1;
      break;
    }
  }
  patchflag=0;
  if(visTimePatch==1){
    for(ii=0;ii<npatch_loaded;ii++){
      i = patch_loaded_list[ii];
      patchi=patchinfo+i;
      if(patchi->display==0||patchi->type!=ipatchtype)continue;
      patchflag=1;
      break;
    }
  }
  partflag=0;
  if(visSmoke==1&&visTimeSmoke==1){
    for(i=0;i<npartinfo;i++){
      parti = partinfo + i;
      if(parti->evac==1)continue;
      if(parti->loaded==0||parti->display==0)continue;
      partflag=1;
      current_particle_type=parti->particle_type;
      if(current_particle_type!=last_particle_type)updatechopcolors();
      break;
    }
  }
  evacflag=0;
  if(visEvac==1&&visTimeEvac==1){
    for(i=0;i<npartinfo;i++){
      parti = partinfo + i;
      if(parti->evac==0)continue;
      if(parti->loaded==0||parti->display==0)continue;
      evacflag=1;
      break;
    }
  }
  if( plotstate==DYNAMIC_PLOTS && 
    ( sliceflag==1 || vsliceflag==1 || partflag==1 || patchflag==1 ||
    smoke3dflag==1|| showtour==1 || evacflag==1||
    (ReadZoneFile==1&&visZone==1&&visTimeZone==1)||
    (ReadTargFile==1&&visTarg==1)
#ifdef pp_WUI
    ||showterrain==1
#endif
    )
    )showtime=1;
    if(plotstate==DYNAMIC_PLOTS&&ReadIsoFile==1&&visAIso!=0&&isoflag==1)showtime2=1;
  if(plotstate==DYNAMIC_PLOTS){
    if(smoke3dflag==1)show3dsmoke=1;
    if(partflag==1)showsmoke=1;
    if(evacflag==1)showevac=1;
    if(patchflag==1)showpatch=1;
    for(i=0;i<selected_case->nmeshes;i++){
      meshi=selected_case->meshinfo+i;
      meshi->visInteriorPatches=0;
    }
    if(showpatch==1&&visPatchType[0]==1){
      for(i=0;i<selected_case->nmeshes;i++){
        meshi=selected_case->meshinfo+i;
        if(meshi->patchtimes==NULL)continue;
        patchi = patchinfo+meshi->patchfilenum;
        if(patchi->loaded==1&&patchi->display==1&&patchi->type==ipatchtype){
          meshi->visInteriorPatches=1;
        }
      }
    }
    if(sliceflag==1)showslice=1;
    if(vsliceflag==1)showvslice=1;
    if(ReadZoneFile==1&&visZone==1&&visTimeZone==1)showzone=1;
    if(ReadIsoFile==1&&visAIso!=0){
      showiso=1;
    }
    if(ReadTargFile==1&&visTarg==1)showtarget=1;
  }
  if(showsmoke==1||showpatch==1||showslice==1||showvslice==1||showzone==1||showiso==1||showevac==1)RenderTime=1;
  if(showtour==1||show3dsmoke==1)RenderTime=1;
  if(plotstate==STATIC_PLOTS&&ReadPlot3dFile==1&&plotn>0&&plotn<=numplot3dvars)showplot3d=1;

  numColorbars=0;
  if(ReadPartFile==1)numColorbars++;
  if(plotstate==DYNAMIC_PLOTS&&(sliceflag==1||vsliceflag==1))numColorbars++;
  if(plotstate==DYNAMIC_PLOTS&&patchflag==1)numColorbars++;
  if(plotstate==DYNAMIC_PLOTS&&ReadZoneFile==1)numColorbars++;
  if(ReadPlot3dFile==1&&numColorbars==0)numColorbars=1;
  /* note: animated iso-contours do not need a color bar,
           so we don't test for isosurface files */
  dwinWW = numColorbars*dwinW/3;
  if(fontindex==1)dwinWW=(int)(1.5*dwinWW);
  drawColorLabel=0;
  if((showtime==1||showplot3d==1)&&visColorLabels==1)drawColorLabel=1;
  if(drawColorLabel==1&&olddrawColorLabel==0)updatemenu=1;
  if(drawColorLabel==0&&olddrawColorLabel==1)updatemenu=1;
  olddrawColorLabel=drawColorLabel;
  if(showtime2==1)showtime=1;
  if(plotstate==DYNAMIC_PLOTS&&stept==1){
  glutIdleFunc(Idle);
  }
  else{
    glutIdleFunc(NULL);
  }

}

/* ------------------ antialias ------------------------ */

void antialias(int flag){
  if(antialiasflag==1){
    if(flag==1){
      glEnable(GL_LINE_SMOOTH);
      glEnable(GL_BLEND);
      glEnable(GL_POINT_SMOOTH);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glHint(GL_LINE_SMOOTH_HINT,GL_DONT_CARE);
    }
    if(flag==0){
      glDisable(GL_LINE_SMOOTH);
      glDisable(GL_POINT_SMOOTH);
      glDisable(GL_BLEND);
    }
  }
}

/* ------------------ transparenton ------------------------ */

void transparenton(void){
  glEnable(GL_BLEND);
  glDepthMask(GL_FALSE);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
//  glBlendFunc(GL_ZERO,GL_ONE_MINUS_SRC_ALPHA);
}

/* ------------------ transparentoff ------------------------ */

void transparentoff(void){
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
}

/* ------------------ updateclipbounds ------------------------ */

void updateclipbounds(int set_i0, int *i0, int set_i1, int *i1, int imax){ 

  if(set_i0==0&&set_i1==0)return;
  if(set_i0==1&&set_i1==1){
    if(*i0>imax-1){*i0=imax-1; *i1=imax;}
    if(*i1>imax)*i1=imax;
    if(*i1<1){*i1=1;*i0=0;}
    if(*i0<0)*i0=0;
    if(*i0>=*i1){*i0=*i1-1;}
  }
  if(set_i0==1&&set_i1==0){
    if(*i0<0)*i0=0;
    if(*i0>imax)*i0=imax;
  }
  if(set_i0==0&&set_i1==1){
    if(*i1<0)*i1=0;
    if(*i1>imax)*i1=imax;
  }
}

/* ------------------ updateclip ------------------------ */

void updateclip(int slicedir){
  stepclip_x=0; stepclip_y=0; stepclip_z=0; 
  stepclip_X=0; stepclip_Y=0; stepclip_Z=0;
  switch (slicedir){
  case 1:
    clip_x = 1 - clip_x;
    if(clip_x==1)printf("clip x on\n");
    if(clip_x==0)printf("clip x off\n");
    if(clip_x==1)stepclip_x=1;
    break;
  case 2:
    clip_y = 1 - clip_y;
    if(clip_y==1)printf("clip y on\n");
    if(clip_y==0)printf("clip y off\n");
    if(clip_y==1)stepclip_y=1;
    break;
  case 3:
    clip_z = 1 - clip_z;
    if(clip_z==1)printf("clip z on\n");
    if(clip_z==0)printf("clip z off\n");
    if(clip_z==1)stepclip_z=1;
    break;
  case -1:
    clip_X = 1 - clip_X;
    if(clip_X==1)printf("clip X on\n");
    if(clip_X==0)printf("clip X off\n");
    if(clip_X==1)stepclip_X=1;
    break;
  case -2:
    clip_Y = 1 - clip_Y;
    if(clip_Y==1)printf("clip Y on\n");
    if(clip_Y==0)printf("clip Y off\n");
    if(clip_Y==1)stepclip_Y=1;
    break;
  case -3:
    clip_Z = 1 - clip_Z;
    if(clip_Z==1)printf("clip Z on\n");
    if(clip_Z==0)printf("clip Z off\n");
    if(clip_Z==1)stepclip_Z=1;
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
}

/* ------------------ updatetimes ------------------------ */

void updatetimes(void){
  int n,n2,ntimes2;
  float *timescopy;
  int i,j,k;
  int listindex;
  slice *sd;
  iso *ib;
  mesh *meshi;
  blockagedata *bc;
  ventdata *vi;
  patch *patchi;
  particle *parti;
  tourdata *touri;
  int filenum;
#ifdef pp_HRR
  float dt_MIN=100000.0;
#endif

  updateShow();  
  ntimes = 0;
#ifdef pp_WUI
  if(visTerrain==1){
    for(i=0;i<nterraininfo;i++){
      terraindata *terri;

      terri = terraininfo + i;
      if(terri->loaded==1)ntimes+=terri->ntimes;
    }
  }
#endif
  for(i=0;i<ntours;i++){
    touri = tourinfo + i;
    if(touri->display==0)continue;
    ntimes += touri->npath;
  }
  for(i=0;i<npartinfo;i++){
    parti = partinfo + i;
    if(parti->loaded==0)continue;
    ntimes += parti->nframes;
  }
  for(i=0;i<nslice;i++){
    sd=sliceinfo+i;
    if(sd->loaded==1||sd->vloaded==1){
      ntimes+=sd->nsteps;
    }
  }
  if(ReadTargFile==1&&visTarg==1){
    ntimes+=ntargtimes;
  }
  for(i=0;i<selected_case->nmeshes;i++){
    meshi=selected_case->meshinfo+i;
    filenum =meshi->patchfilenum;
    if(filenum!=-1){
      patchi=patchinfo+filenum;
      if(patchi->loaded==1){
        ntimes+=meshi->npatch_frames;
      }
    }
  }
  if(ReadZoneFile==1&&visZone==1){
    ntimes+=nzonet;
  }
  if(ReadIsoFile==1&&visAIso!=0){
    for(i=0;i<selected_case->nmeshes;i++){
      meshi=selected_case->meshinfo+i;
      if(meshi->isofilenum<0)continue;
      ib = isoinfo + meshi->isofilenum;
      if(ib->loaded==0)continue;
      ntimes+=meshi->nisosteps;
    }
  }
  {
    smoke3d *smoke3di;

    if(Read3DSmoke3DFile==1&&vis3DSmoke3D==1){
      for(i=0;i<nsmoke3d;i++){
        smoke3di = smoke3dinfo + i;
        if(smoke3di->loaded==0)continue;
        ntimes += smoke3di->n_times;
      }
    }
  }

  if(ntimes==0)return;
  CheckMemory;
  FREEMEMORY(times);
  NewMemory((void **)&times,ntimes*sizeof(float));
  timescopy=times;
#ifdef pp_WUI
  if(visTerrain==1){
    for(i=0;i<nterraininfo;i++){
      terraindata *terri;

      terri = terraininfo + i;
      if(terri->loaded==0)continue;
      for(n=0;n<terri->ntimes;n++){
        *timescopy++=terri->times[n];
#ifdef pp_HRR
        if(n!=0&&timescopy[0]-timescopy[-1]<dt_MIN){
          dt_MIN=timescopy[0]-timescopy[-1];
        }
#endif
      }
    }
  }
#endif
  for(i=0;i<ntours;i++){
    touri = tourinfo + i;
    if(touri->display==0)continue;
    for(n=0;n<touri->npath;n++){
      *timescopy++=touri->path_times[n];
#ifdef pp_HRR
      if(n!=0&&timescopy[0]-timescopy[-1]<dt_MIN){
        dt_MIN=timescopy[0]-timescopy[-1];
      }
#endif
    }
  }
  for(i=0;i<npartinfo;i++){
    parti = partinfo + i;
    if(parti->loaded==0)continue;
    for(n=0;n<parti->nframes;n++){
      *timescopy++=parti->ptimes[n];
#ifdef pp_HRR
      if(n!=0&&timescopy[0]-timescopy[-1]<dt_MIN){
        dt_MIN=timescopy[0]-timescopy[-1];
      }
#endif
    }
  }
  for(i=0;i<nslice;i++){
    sd = sliceinfo + i;
    if(sd->loaded==1||sd->vloaded==1){
        for(n=0;n<sd->nsteps;n++){
          *timescopy++=sd->slicetimes[n];
#ifdef pp_HRR
          if(n!=0&&timescopy[0]-timescopy[-1]<dt_MIN){
            dt_MIN=timescopy[0]-timescopy[-1];
          }
#endif
        }
    }
  }
  if(ReadTargFile==1&&visTarg==1){
    for(n=0;n<ntargtimes;n++){
      *timescopy++=targtimes[n];
#ifdef pp_HRR
      if(n!=0&&timescopy[0]-timescopy[-1]<dt_MIN){
        dt_MIN=timescopy[0]-timescopy[-1];
      }
#endif
    }
  }
  for(i=0;i<selected_case->nmeshes;i++){
    meshi=selected_case->meshinfo + i;
    filenum=meshi->patchfilenum;
    if(filenum!=-1){
      patchi = patchinfo + filenum;
      if(patchi->loaded==1){
        for(n=0;n<meshi->npatch_frames;n++){
          *timescopy++=meshi->patchtimes[n];
#ifdef pp_HRR
          if(n!=0&&timescopy[0]-timescopy[-1]<dt_MIN){
            dt_MIN=timescopy[0]-timescopy[-1];
          }
#endif
        }
      }
    }
  }
  if(ReadZoneFile==1&&visZone==1){
    for(n=0;n<nzonet;n++){
      *timescopy++=zonet[n];
#ifdef pp_HRR
      if(n!=0&&timescopy[0]-timescopy[-1]<dt_MIN){
        dt_MIN=timescopy[0]-timescopy[-1];
      }
#endif
    }
  }
  if(ReadIsoFile==1&&visAIso!=0){
    for(i=0;i<niso;i++){
      ib = isoinfo+i;
      if(ib->loaded==0)continue;
      meshi=selected_case->meshinfo + ib->blocknumber;
      for(n=0;n<meshi->nisosteps;n++){
        *timescopy++=meshi->isotimes[n];
#ifdef pp_HRR
        if(n!=0&&timescopy[0]-timescopy[-1]<dt_MIN){
          dt_MIN=timescopy[0]-timescopy[-1];
        }
#endif
      }
    }
  }
  {
    smoke3d *smoke3di;

    if(Read3DSmoke3DFile==1&&vis3DSmoke3D==1){
      for(i=0;i<nsmoke3d;i++){
        smoke3di = smoke3dinfo + i;
        if(smoke3di->loaded==0)continue;
        for(n=0;n<smoke3di->n_times;n++){
          *timescopy++=smoke3di->times[n];
#ifdef pp_HRR
          if(n!=0&&timescopy[0]-timescopy[-1]<dt_MIN){
            dt_MIN=timescopy[0]-timescopy[-1];
          }
#endif
        }
      }
    }
  }

  qsort( (float *)times, (size_t)ntimes, sizeof( float ), compare );

  n2=1;ntimes2=ntimes;
  for(n=1;n<ntimes;n++){
#ifdef pp_HRR
    if(fabs(times[n]-times[n-1])>dt_MIN){
#else
    if(times[n]!=times[n-1]){
#endif
      times[n2]=times[n];
      n2++;
    }
    else{
      ntimes2--;
    }
  }
  ntimes=ntimes2;
  if(ntimes>0){
    FREEMEMORY(render_frame);
    NewMemory((void **)&render_frame,ntimes*sizeof(int));
    for(i=0;i<npartinfo;i++){
      parti=partinfo+i;
      FREEMEMORY(parti->ptimeslist);
      NewMemory((void **)&parti->ptimeslist,ntimes*sizeof(int));
    }
    for(i=0;i<ntours;i++){
      touri=tourinfo + i;
      if(touri->display==0)continue;
      FREEMEMORY(touri->path_timeslist);
      NewMemory((void **)&touri->path_timeslist,ntimes*sizeof(int));
    }
#ifdef pp_WUI
    if(visTerrain==1){
      for(i=0;i<nterraininfo;i++){
        terraindata *terri;

        terri = terraininfo + i;
        if(terri->loaded==0)continue;
        FREEMEMORY(terri->timeslist);
        NewMemory((void **)&terri->timeslist,ntimes*sizeof(int));
      }
    }
#endif
#ifdef pp_HRR
    if(hrrinfo!=NULL){
      if(hrrinfo->loaded==1&&hrrinfo->display==1&&ntimes>0){
        int jstart=0;

        FREEMEMORY(hrrinfo->timeslist);
        NewMemory((void **)&hrrinfo->timeslist,ntimes*sizeof(int));
        FREEMEMORY(hrrinfo->times);
        NewMemory((void **)&hrrinfo->times,ntimes*sizeof(float));
        FREEMEMORY(hrrinfo->hrrval);
        NewMemory((void **)&hrrinfo->hrrval,ntimes*sizeof(float));
        hrrinfo->ntimes=ntimes;
        for(i=0;i<ntimes;i++){
          int j, foundit;

          foundit=0;
          hrrinfo->times[i]=times[i];
          for(j=jstart;j<hrrinfo->ntimes_csv-1;j++){
            if(hrrinfo->times_csv[j]<=times[i]&&times[i]<hrrinfo->times_csv[j+1]){
              float f1, tbot;

              foundit=1;
              tbot = hrrinfo->times_csv[j+1]-hrrinfo->times_csv[j];
              if(tbot>0.0){
                f1 = (times[i]-hrrinfo->times_csv[j])/tbot;
              }
              else{
                f1=0.0;
              }
              hrrinfo->hrrval[i]=(1.0-f1)*hrrinfo->hrrval_csv[j]+f1*hrrinfo->hrrval_csv[j+1];
              jstart=j;
              break;
            }
          }
          if(foundit==0){
            hrrinfo->hrrval[i]=hrrinfo->hrrval_csv[hrrinfo->ntimes_csv-1];
          }
        }

      }
    }
#endif
    for(i=0;i<nslice;i++){
      sd = sliceinfo + i;
      FREEMEMORY(sd->slicetimeslist);
      NewMemory((void **)&sd->slicetimeslist,ntimes*sizeof(int));
    }
    {
      smoke3d *smoke3di;

      for(i=0;i<nsmoke3d;i++){
        smoke3di = smoke3dinfo + i;
        FREEMEMORY(smoke3di->timeslist);
        NewMemory((void **)&smoke3di->timeslist,ntimes*sizeof(int));
      }
    }
    for(i=0;i<selected_case->nmeshes;i++){
      meshi=selected_case->meshinfo+i;
      if(meshi->isotimes==NULL)continue;
      FREEMEMORY(meshi->isotimeslist);
      NewMemory((void **)&meshi->isotimeslist,  ntimes*sizeof(int));  
    }

    for(i=0;i<selected_case->nmeshes;i++){
      FREEMEMORY(selected_case->meshinfo[i].patchtimeslist); 
    }
    for(i=0;i<selected_case->nmeshes;i++){
      if(selected_case->meshinfo[i].patchtimes==NULL)continue;
      NewMemory((void **)&selected_case->meshinfo[i].patchtimeslist,ntimes*sizeof(int));
    }

    FREEMEMORY(zonetlist); 
    NewMemory((void **)&zonetlist,     ntimes*sizeof(int));

    FREEMEMORY(targtimeslist);
    NewMemory((void **)&targtimeslist,  ntimes*sizeof(int));

    if(ntotal_smooth_blockages>0){
      for(i=0;i<selected_case->nmeshes;i++){
        meshi=selected_case->meshinfo+i;
        FREEMEMORY(meshi->showsmoothtimelist);
        NewMemory((void **)&meshi->showsmoothtimelist,ntimes*sizeof(smoothblockage *));
      }
    }


  for(n=0;n<ntimes;n++){render_frame[n]=0;}
    ResizeMemory((void **)&times,ntimes*sizeof(float));
  }
  izone=0; itime=0;
  for(i=0;i<selected_case->nmeshes;i++){
    meshi=selected_case->meshinfo+i;
    meshi->ipatch=0;
  }
  for(i=0;i<nslice;i++){
    sd = sliceinfo + i;
    sd->islice=0; 
  }
  iframe=iframebeg; 
  for(i=0;i<selected_case->nmeshes;i++){
    meshi=selected_case->meshinfo+i;
    if(meshi->isotimes==NULL)continue;
    meshi->iiso=0;
  }

  /* determine visibility of each blockage at each time step */

  for(i=0;i<selected_case->nmeshes;i++){
    meshi=selected_case->meshinfo+i;
    for(j=0;j<meshi->nbptrs;j++){
      bc = meshi->blockageinfoptrs[j];
      if(bc->showtime==NULL)continue;
      FREEMEMORY(bc->showtimelist);
      if(ntimes>0){
        NewMemory((void **)&bc->showtimelist,ntimes*sizeof(int));
        for(k=0;k<ntimes;k++){
          bc->showtimelist[k]=1;
          listindex=getindex(times[k],bc->showtime,bc->nshowtime);
          bc->showtimelist[k]=bc->showhide[listindex];
        }
      }
    }
  }

  /* determine visibility of each vent at each time step */

  for(i=0;i<selected_case->nmeshes;i++){
    meshi=selected_case->meshinfo+i;
    if(meshi->ventinfo==NULL)continue;
    for(j=0;j<meshi->nvents;j++){
      vi = meshi->ventinfo+j;
      if(vi->showtime==NULL)continue;
      FREEMEMORY(vi->showtimelist);
      if(ntimes>0){
        NewMemory((void **)&vi->showtimelist,ntimes*sizeof(int));
        for(k=0;k<ntimes;k++){
          vi->showtimelist[k]=1;
          listindex=getindex(times[k],vi->showtime,vi->nshowtime);
          vi->showtimelist[k]=vi->showhide[listindex];
        }
      }
    }
  }

  /* determine visibility of each smooth blockage */

//  for(i=0;i<selected_case->nmeshes;i++){
//    meshi=selected_case->meshinfo+i;
//  }
  
  synctimes();
  updatefaces=1;
  UpdateTimeLabels();
}

/* ------------------ getindex ------------------------ */

int getindex(float key, const float *list, int nlist){
  int i;
  if(nlist==1)return 0;
  if(key<list[1])return 0;
  if(key>=list[nlist-1])return nlist-1;
  for(i=1;i<nlist-1;i++){
    if(list[i]<=key&&key<list[i+1])return i;
  }
  return 0;
}

/* ------------------ compare ------------------------ */

int compare( const void *arg1, const void *arg2 ){
  float x, y;
  x=*(float *)arg1;
  y=*(float *)arg2;
  if( x< y)return -1;
  if( x==y)return 0;
  return 1;
}

/* ------------------ drawTimeBar ------------------------ */

void drawTimeBar(void){
  float xleft=.175f, xright=1.0f, ybot=0.10f, ytop=.35f, xxright;

  xleft = xtimeleft;
  if(fontindex==LARGE_FONT)xleft=xtimeleft+0.11;
  xright = xtimeright;

  glLineWidth(linewidth);
  glBegin(GL_LINE_LOOP);
  glColor4fv(timebarcolor);
  glVertex2f(xleft,ybot);
  glVertex2f(xright,ybot);
  glVertex2f(xright,ytop);
  glVertex2f(xleft,ytop);
  glEnd();

  if(ntimes != 1){
    xxright = xleft + (float)itime*(xright-xleft)/(ntimes-1);
  }
  else{
    xxright=xright;
  }
  glBegin(GL_POLYGON);
  glColor4fv(timebarcolor);
  glVertex2f(xleft,ybot);
  glVertex2f(xxright,ybot);
  glVertex2f(xxright,ytop);
  glVertex2f(xleft,ytop);
  glEnd();

}

/* ------------------ getunitinfo ------------------------ */

void getunitinfo(const char *shortlabel, int *unitclass, int *unittype){
  int i,j;
  *unitclass=-1;
  *unittype=-1;
  for(i=0;i<nunitclasses;i++){
    for(j=0;j<unitclasses[i].nunittypes;j++){
      if(STRCMP(unitclasses[i].unittypes[j],shortlabel)==0){
        *unitclass=i;
        *unittype=unitclasses[i].active;
        return;
      }
    }
  }
}

/* ------------------ setsmokeviewvars ------------------------ */

void setsmokeviewvars(void){
  int i;

  for(i=0;i<20;i++){
    cputimes[i]=0.0;
  }
  eyexINI=0.0; 
  eyeyINI=0.0;
  eyezINI=0.0;
  for(i=0;i<16;i++){
    modelview_setup[i]=0.0;
  }
  for(i=0;i<4;i++){
    modelview_setup[i+4*i]=1.0;
  }
}

/* ------------------ Init ------------------------ */

void Init(void){
  int errorcode;
  mesh *meshi;

  int n,i;

  FREEMEMORY(plotiso);
  NewMemory((void **)&plotiso,mxplot3dvars*sizeof(int));

  for(n=0;n<mxplot3dvars;n++){
    plotiso[n]=nrgb/2;
  }

  for(i=0;i<16;i++){
    modelview_setup[i]=0.0;
  }
  for(i=0;i<4;i++){
    modelview_setup[i+4*i]=1.0;
  }

  for(i=0;i<selected_case->nmeshes;i++){
    meshi=selected_case->meshinfo+i;
    initcontour(&meshi->plot3dcontour1,rgbptr,nrgb+1);
    initcontour(&meshi->plot3dcontour2,rgbptr,nrgb+1);
    initcontour(&meshi->plot3dcontour3,rgbptr,nrgb+1);
  }

  if(set_no_part!=1&&nopart!=1&&npartinfo>0){
    readpart(partinfo[0].file,0,LOAD,&errorcode);
  }
  for(i=0;i<selected_case->nmeshes;i++){
    meshi=selected_case->meshinfo+i;
    meshi->currentsurf.defined=0;
    meshi->currentsurf2.defined=0;
  }

  /* initialize box sizes, lighting parameters */

  xyzbox = xbar;
  if(ybar>xyzbox){xyzbox=ybar;}
  if(zbar>xyzbox){xyzbox=zbar;}

  {
    char name_external[32];
    strcpy(name_external,"external");
    init_camera(camera_external,name_external);
    camera_external->view_id=1;
  }
  if(camera_ini->defined==1){
    copy_camera(camera_current,camera_ini);
  }
  else{
    camera_external->zoom=zoom;
    copy_camera(camera_current,camera_external);
  }
  strcpy(camera_label,camera_current->name);
  update_camera_label();
  {
    char name_internal[32];
    strcpy(name_internal,"internal");
    init_camera(camera_internal,name_internal);
  }
  camera_internal->eye[0]=0.5*xbar;
  camera_internal->eye[1]=0.5*ybar;
  camera_internal->eye[2]=0.5*zbar;
  camera_internal->view_id=0;
  copy_camera(camera_save,camera_current);
  copy_camera(camera_last,camera_current);

  init_camera_list();
  add_default_views();
  update_view_gluilist();

  //reset_glui_view(i_view_list);

  screenWidth2 = screenWidth - dwinW;
  screenHeight2 = screenHeight - dwinH;

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);
  if(cullfaces==1)glEnable(GL_CULL_FACE);

  glClearColor(backgroundcolor[0],backgroundcolor[1],backgroundcolor[2], 0.0f);
  glShadeModel(GL_SMOOTH); 
  glDisable(GL_DITHER);


  thistime=0;
  lasttime=0;

  /* define color bar */

  updatecolors(-1);

  block_ambient2[3] = 1.0;
  mat_ambient2[3] = 1.0;
  mat_specular2[3] = 1.0;

  reset_glui_view(startup_view_ini);
  updateShow();
}

/* ------------------ ResetView ------------------------ */

void ResetView(int option){
  int eyeview_save;

  switch (option){
  case RESTORE_EXTERIOR_VIEW_ZOOM:
  case RESTORE_EXTERIOR_VIEW:
    eyeview_save = camera_current->eyeview;
    copy_camera(camera_current,camera_external);
    camera_current->eyeview=eyeview_save;
    if(projection_type==1)camera_current->eye[1]=camera_current->isometric_y;
    break;
  case RESTORE_INTERIOR_VIEW:
    eyeview_save = camera_current->eyeview;
    copy_camera(camera_current,camera_internal);
    camera_current->eyeview=eyeview_save;
    break;
  case RESTORE_SAVED_VIEW:
    copy_camera(camera_current,camera_save);
    break;
  case 3:
  case 4:
    ASSERT(FFALSE);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  if(option==RESTORE_EXTERIOR_VIEW_ZOOM)camera_current->zoom=zooms[zoomindex];
  zoom=camera_current->zoom;
  update_glui_zoom();
}

/* ------------------ UpdateTimeLabels ------------------------ */

void UpdateTimeLabels(void){
  float time0;
  int hour, min, sec,sec10;
  time0 = timeoffset;


  if(times!=NULL)time0 = timeoffset + times[itime];
  if(vishmsTimelabel==1){
    hour = time0/3600;
    min = time0/60.0 - 60*hour;
    sec10 = 10*(time0 -  60*min - 3600*hour);
    sec = sec10/10;
    sec10 = sec10 - 10*sec;
    sprintf(timelabel,"  %i:%.2i:%.2i.%i",hour,min,sec,sec10);
  }
  else{
    float dt;

    dt=times[1]-times[0];
    if(dt<0.0)dt=-dt;
    if(dt<0.001){
      sprintf(timelabel,"Time: %4.4f",time0);
    }
    else if(dt>=0.001&&dt<0.01){
      sprintf(timelabel,"Time: %4.3f",time0);
    }
    else if(dt>=0.01&&dt<0.1){
      sprintf(timelabel,"Time: %4.2f",time0);
    }
    else{
      sprintf(timelabel,"Time: %4.1f",time0);
    }
  }
  sprintf(framelabel,"Frame: %i",itime);
#ifdef pp_HRR
  if(hrrinfo!=NULL&&hrrinfo->display==1&&hrrinfo->loaded==1){
    float hrr;

    hrr = hrrinfo->hrrval[hrrinfo->itime];
    if(hrr<1.0){
      sprintf(hrrinfo->hrrlabel,"HRR: %4.1f",hrr*1000.0);
    }
    else if(hrr>1000.0){
      sprintf(hrrinfo->hrrlabel,"HRR: %4.1f MW",hrr/1000.0);
    }
    else{
      sprintf(hrrinfo->hrrlabel,"HRR: %4.1f kW",hrr);
    }
  }
#endif
}

/* ------------------ synctimes ------------------------ */

void synctimes(void){
  int n,i,istart;
  int j,igrid,jj;
  slice *sd;
  particle *parti;
  mesh *meshi;


  /* synchronize smooth blockage times */

    //       meshi->nsmoothcolors;
    //       meshi->blockagesurfaces
  if(ntotal_smooth_blockages>0){
    for(igrid=0;igrid<selected_case->nmeshes;igrid++){
      meshi=selected_case->meshinfo+igrid;
      if(meshi->showsmoothtimelist==NULL)continue;
      for(n=0;n<ntimes;n++){
        smoothblockage *sb;

        sb = getsmoothblockage(meshi,times[n]);
        meshi->showsmoothtimelist[n] = sb;
      }
    }
  }

  for(n=0;n<ntimes;n++){

  /* synchronize tour times */

    {
      tourdata *tourj; 
      for(j=0;j<ntours;j++){
        tourj = tourinfo + j;
        if(tourj->display==0)continue;
        if(n==0){
          istart=0;
        }
        else{
          istart=tourj->path_timeslist[n-1];
        }
        i=istart;
        while(tourj->path_times[i]<times[n]&&i<tourj->npath){
          i++;
        }
        if(i>=tourj->npath){
          i--;
        }
        tourj->path_timeslist[n]=i;
      }
    }
#ifdef pp_WUI
  /* synchronize terrain times */

    for(j=0;j<nterraininfo;j++){
      terraindata *terri;

      terri = terraininfo + j;
      if(terri->loaded==0)continue;
      if(n==0){
        istart=0;
      }
      else{
        istart=terri->timeslist[n-1];
      }
      i=istart;
      while(terri->times[i]<times[n]&&i<terri->ntimes){
        i++;
      }
      if(i>=terri->ntimes){
        i--;
      }
      terri->timeslist[n]=i;
    }
#endif
#ifdef pp_HRR
    if(hrrinfo!=NULL&&hrrinfo->loaded==1&&hrrinfo->display==1){
      if(n==0){
        istart=0;
      }
      else{
        istart=hrrinfo->timeslist[n-1];
      }
      i=istart;
      while(hrrinfo->times[i]<times[n]&&i<hrrinfo->ntimes){
        i++;
      }
      if(i>=hrrinfo->ntimes){
        i--;
      }
      hrrinfo->timeslist[n]=i;
    }
#endif

  /* synchronize particle times */

    for(j=0;j<npartinfo;j++){
      parti=partinfo+j;
      if(parti->loaded==0)continue;
      if(n==0){
        istart=0;
      }
      else{
        istart=parti->ptimeslist[n-1];
      }
      i=istart;
      while(parti->ptimes[i]<times[n]&&i<parti->nframes){
        i++;
      }
      if(i>=parti->nframes){
        i--;
      }
      parti->ptimeslist[n]=i;
    }
       /* synchronize target times */

    if(ntarg_files>0){
      if(n==0){
        istart=0;
      }
      else{
        istart=targtimeslist[n-1];
      }
      i=istart;
      while(targtimes[i]<times[n]&&i<ntargtimes){
        i++;
      }
      if(i>=ntargtimes){
        i--;
      }
      targtimeslist[n]=i;
    }

  /* synchronize slice times */

    for(jj=0;jj<nslice_loaded;jj++){
      j = slice_loaded_list[jj];
      sd = sliceinfo + j;
      if(n==0){
        istart=0;
      }
      else{
        istart=sd->slicetimeslist[n-1];
      }
      i=istart;
      while(sd->slicetimes[i]<times[n]&&i<sd->nsteps){
        i++;
      }
      if(i>=sd->nsteps){
        i=sd->nsteps-1;
      }
      sd->slicetimeslist[n]=i;
    }

  /* synchronize smoke times */
    {
      smoke3d *smoke3di;

      for(jj=0;jj<nsmoke3d;jj++){
        smoke3di = smoke3dinfo + jj;
        if(smoke3di->loaded==0)continue;
         if(n==0){
          istart=0;
         }
        else{
          istart=smoke3di->timeslist[n-1];
        }
        i=istart;
        while(smoke3di->times[i]<times[n]&&i<smoke3di->n_times){
          i++;
        }
        if(i>=smoke3di->n_times){
          i=smoke3di->n_times-1;
        }
        smoke3di->timeslist[n]=i;
      }
    }

  /* synchronize patch times */

    for(j=0;j<selected_case->nmeshes;j++){
      meshi=selected_case->meshinfo+j;
      if(meshi->patchtimes==NULL)continue;
      if(n==0){
        istart=0;
      }
      else{
        istart=meshi->patchtimeslist[n-1];
      }
      i=istart;
      while(meshi->patchtimes[i]<times[n]&&i<meshi->npatch_frames){
        i++;
      }
      if(i>=meshi->npatch_frames){
        i=meshi->npatch_frames-1;
      }
      meshi->patchtimeslist[n]=i;
    }

  /* synchronize isosurface times */

    for(igrid=0;igrid<selected_case->nmeshes;igrid++){
      meshi=selected_case->meshinfo+igrid;
      if(meshi->isotimes==NULL)continue;
      if(n==0){
        istart=0;
      }
      else{
        istart=meshi->isotimeslist[n-1];
      }
      i=istart;
      while(meshi->isotimes[i]<times[n]&&i<meshi->nisosteps){
        i++;
      }
      if(i>=meshi->nisosteps){
        i=meshi->nisosteps-1;
      }
      meshi->isotimeslist[n]=i;
    }

    /* synchronize zone times */

    if(showzone==1){
      if(n==0){
        istart=0;
      }
      else{
        istart=zonetlist[n-1];
      }
      i=istart;
      while(zonet[i]<times[n]&&i<nzonet){
        i++;
      }
      if(i>=nzonet)i=nzonet-1;
      zonetlist[n]=i;
    }

  }
  reset_gltime();
}

/* ------------------ RenderFrame ------------------------ */

#ifndef pp_nolibs
void RenderFrame(int view_mode){
  char renderfile[1024],renderfile2[1024];
  FILE *stream;
  char *ext;

  if(view_mode==VIEW_LEFT&&showstereo==2)return;
  switch (renderfiletype){
  case 0:
    ext=ext_png;
    break;
  case 1:
    ext=ext_jpg;
    break;
#ifdef pp_GDGIF
  case 2:
    ext=ext_gif;
    break;
#endif
  default:
    renderfiletype=2;
    ext=ext_png;
    break;
  }

  glFlush();
  if(RenderTime==1){
    switch (view_mode){
    case VIEW_CENTER:
      sprintf(renderfile,"%s_%04i",fdsprefix,itime/RenderSkip);
      break;
    case VIEW_LEFT:
      sprintf(renderfile,"%s_%04i_L",fdsprefix,itime/RenderSkip);
      break;
    case VIEW_RIGHT:
      if(showstereo==2){
        sprintf(renderfile,"%s_%04i",fdsprefix,itime/RenderSkip);
      }
      else{
        sprintf(renderfile,"%s_%04i_R",fdsprefix,itime/RenderSkip);
      }
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
  }
  if(RenderTime==0){
    switch (view_mode){
    case VIEW_CENTER:
      sprintf(renderfile,"%s_s%04i",fdsprefix,seqnum);
      seqnum++;
      break;
    case VIEW_LEFT:
      sprintf(renderfile,"%s_s%04i_L",fdsprefix,seqnum);
      break;
    case VIEW_RIGHT:
      if(showstereo==2){
        sprintf(renderfile,"%s_s%04i",fdsprefix,seqnum);
      }
      else{
        sprintf(renderfile,"%s_s%04i_R",fdsprefix,seqnum);
      }
      seqnum++;
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
  }
  strcat(renderfile,ext);

  // if there is a tempdir see if we need to use it

  if(smokeviewtempdir!=NULL){
    stream=fopen(renderfile,"wb");
    if(stream==NULL){
      strcpy(renderfile2,smokeviewtempdir);
      strcat(renderfile2,renderfile);
      stream=fopen(renderfile2,"wb");
      if(stream!=NULL)strcpy(renderfile,renderfile2);
    }
    if(stream!=NULL)fclose(stream);
  }
  SVimage2file(renderfile,renderfiletype,screenWidth,screenHeight);
  if(RenderTime==1&&output_slicedata==1){
    output_Slicedata();
  }
}
#endif

/* ------------------ ClearBuffers ------------------------ */

void ClearBuffers(int mode){
  if(mode==RENDER){
    glClearColor(backgroundcolor[0],backgroundcolor[1],backgroundcolor[2], 0.0f);
  }
  else{
    glClearColor((float)0.0,(float)0.0,(float)0.0, (float)0.0);
  }
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

/* ------------------ Args ------------------------ */

void Args(int argc, char **argv){
  int i, len;
  char *temp;
  char buffer[255];
  int iarg;
  struct stat statbuffer;
  size_t len_memory;
  char *argi;
  char SMVFILENAME[1024];
  int smv_parse;

  STRCPY(TRAINERTITLE,TRAINERTITLEBASE);
  STRCAT(TRAINERTITLEBASE,__DATE__);
#ifdef _DEBUG
  STRCPY(TITLE,TITLEBASE);
  STRCAT(TITLE,__DATE__);
#else
  STRCPY(TITLE,TITLEBASE);
  STRCAT(TITLE,__DATE__);
#endif

  STRCPY(FULLTITLE,TITLE);

  STRCPY(TITLERELEASE,TITLE);

  CheckMemory;
  partscale=a_partscale;
  zonescale=a_zonescale;

  if(argc==1){
    usage(argv);
    exit(1);
  }
  if(strncmp(argv[1],"-ini",3)==0){
    InitOpenGL();
    updatecolors(-1);
    writeini(GLOBAL_INI);
    exit(0);
  }

  strcpy(SMVFILENAME,"");
  smv_parse=0;
  for(iarg=1;iarg<argc;iarg++){
    argi=argv[iarg];
    if(strncmp(argi,"-",1)==0){
      if(smv_parse==0)continue;
      if(smv_parse==1)break;
    }
    if(smv_parse==1)strcat(SMVFILENAME," ");
    smv_parse=1;
    strcat(SMVFILENAME,argi);
  }

  argi=SMVFILENAME;
#ifndef pp_OSX
  argi=lastname(argi);
#endif
  len = (int) strlen(argi);
  CheckMemory;
  FREEMEMORY(fdsprefix);
  len_memory=len+strlen(part_ext)+100;
  NewMemory((void **)&fdsprefix,(unsigned int)len_memory);
  STRCPY(fdsprefix,argi);
  FREEMEMORY(smvfilename);
  FREEMEMORY(trainer_filename);
  FREEMEMORY(shellfilename);

  strcpy(inputfilename_ext,"");

  if(len>4){
    char *c_ext;

    c_ext=strrchr(argi,'.');
    if(c_ext!=NULL){
      STRCPY(inputfilename_ext,c_ext);
      to_lower(inputfilename_ext);

      if(c_ext!=NULL&&(strcmp(inputfilename_ext,".smv")==0||strcmp(inputfilename_ext,".svt")==0)){
        c_ext[0]=0;
        STRCPY(fdsprefix,argi);
        FREEMEMORY(trainer_filename);
        NewMemory((void **)&trainer_filename,(unsigned int)(len+7));
        STRCPY(trainer_filename,argi);
        STRCAT(trainer_filename,".svt");
      }
    }
  }

  FREEMEMORY(casefilename);
  NewMemory((void **)&casefilename,len+strlen(ini_ext)+1);
  STRCPY(casefilename,fdsprefix);
  STRCAT(casefilename,ini_ext);

  if(smvfilename==NULL){
    NewMemory((void **)&smvfilename,(unsigned int)(len+6));
    FREEMEMORY(smvmenufile);
    NewMemory((void **)&smvmenufile,(unsigned int)(len+15));
    STRCPY(smvfilename,fdsprefix);
    STRCAT(smvfilename,".smv");
    STRCPY(smvmenufile,"Reload ");
    temp = strrchr(smvfilename,(int)(*dirseparator));
    if(temp!=NULL){
      STRCAT(smvmenufile,temp+1);
    }
    else{
      STRCAT(smvmenufile,smvfilename);
    }
  }
  if(smvfilename!=NULL){
    FREEMEMORY(fds_filein);
    NewMemory((void **)&fds_filein,strlen(fdsprefix)+6);
    STRCPY(fds_filein,fdsprefix);
    STRCAT(fds_filein,".data");
    if(stat(fds_filein,&statbuffer)!=0){
      FREEMEMORY(fds_filein);
    }
    if(fds_filein!=NULL){
      getnewfilename();
    }
  }

  if(trainer_filename==NULL){
    NewMemory((void **)&trainer_filename,(unsigned int)(len+6));
    STRCPY(trainer_filename,fdsprefix);
    STRCAT(trainer_filename,".svt");
  }
  if(shellfilename==NULL){
    NewMemory((void **)&shellfilename,(unsigned int)(len+6));
    STRCPY(shellfilename,fdsprefix);
    STRCAT(shellfilename,".shl");
  }


  set_no_part=0;
  for (i=1;i<argc;i++){
    if(strncmp(argv[i],"-",1)!=0)continue;
    if(strncmp(argv[i],"-ini",3)==0)writeini(GLOBAL_INI);
    else if(strncmp(argv[i],"-part",5)==0){
      set_no_part=0; 
    }
    else if(strncmp(argv[i],"-nopart",7)==0){
      set_no_part=1; 
    }
    else if(strncmp(argv[i],"-stereo",7)==0){
      if(benchmark==0){
        stereoactive=1;
        showstereo=1;
        printf("stereo option activated\n");
      }
    }
    else if(strncmp(argv[i],"-benchmark",10)==0){
      benchmark=1;
      buffertype=SINGLE_BUFFER;
      printf("benchmark option activated\n");
      if(stereoactive==1){
        stereoactive=0;
        showstereo=0;
        printf("stereo option deactivated\n");
      }
    }
    else if(strncmp(argv[i],"-points",7)==0){
      ++i;
      if(i<argc){
        mxpoints_comm = atol(argv[i]);mxpoints=mxpoints_comm;
      }
    }
    else if(strncmp(argv[i],"-frames",7)==0){
      ++i;
      if(i<argc){
        mxframes_comm = atol(argv[i]);mxframes=mxframes_comm;
      }
    }
    else if(strncmp(argv[i],"-isotest",8)==0){
      isotest=1;
    }
    else if(strncmp(argv[i],"-help",5)==0){
      usage(argv);
      exit(0);
    }
    else if(strncmp(argv[i],"-build",6)==0){
      showbuild=1;
      usage(argv);
      exit(0);
    }
    else {
      printf(" unknown option: %s\n",argv[i]);usage(argv);scanf("%s",buffer);exit(1);
    }
  }
}

/* ------------------ usage ------------------------ */

void usage(char **argv){
  printf("%s\n",TITLERELEASE);
  printf("Visualize fire/smoke flow simulations.  All parameters are optional.\n\n");
  printf("Usage:\n\n");
  printf("%s casename -points m -frames n -ini -part -nopart -stereo\n\n",argv[0]);
  printf("where \n\n");
  printf("  casename = project id (file names without the extension)\n");
  printf("         m = maximum number of particles.  Default=%i\n",MAXPOINTS);
  printf("         n = maximum number of particle frames.  Default=%i\n",MAXFRAMES);
  printf("      -ini = output default smokeview parameters to smokeview.ini\n");
  printf("     -part = load particle file if present \n");
  printf("   -nopart = do not load particle file \n");
  printf("   -stereo = activate stereo mode (if supported)\n");
  printf("    -build = show pre-preprocessing directives used to build smokeview\n");
  if(showbuild==1){
    printf("  \n");
    printf("  Smokeview was built with the following pre-processing directives set:\n");
#ifdef pp_release
    printf("   pp_release");
#endif
#ifdef _DEBUG
    printf("   _DEBUG");
#endif
#ifdef EGZ
    printf(", EGZ");
#endif
#ifdef ISO_DEBUG
    printf(", ISO_DEBUG");
#endif
#ifdef isotest
    printf(", isotest");
#endif
#ifdef pp_COLOR
    printf(", pp_COLOR");
#endif
#ifdef pp_COMPRESS
    printf(", pp_COMPRESS");
#endif
#ifdef pp_cvf
    printf(", pp_cvf");
#endif
#ifdef pp_DRAWISO
    printf(", pp_DRAWISO");
#endif
#ifdef pp_GPU
    printf(", pp_GPU");
#endif
#ifdef pp_LAHEY
    printf(", pp_LAHEY");
#endif
#ifdef pp_LIGHT
    printf(", pp_LIGHT");
#endif
#ifdef pp_memstatus
    printf(", pp_memstatus");
#endif
#ifdef pp_NISTREVIEW
    printf(", pp_NISTREVIEW");
#endif
#ifdef pp_noappend
    printf(", pp_noappend");
#endif
#ifdef pp_nofortran
    printf(", pp_nofortran");
#endif
#ifdef pp_NOSHOW
    printf(", pp_NOSHOW");
#endif
#ifdef pp_OSX
    printf(", pp_OSX");
#endif
#ifdef pp_SMOKETEST
    printf(", pp_SMOKETEST");
#endif
#ifdef pp_SPHERE
    printf(", pp_SPHERE");
#endif
#ifdef pp_SVNET
    printf(", pp_SVNET");
#endif
#ifdef pp_TEST
    printf(", pp_TEST");
#endif
#ifdef pp_TEXTUREPLOT
    printf(", pp_TEXTUREPLOT");
#endif
#ifdef pp_THREADS2
    printf(", pp_THREADS2");
#endif
#ifdef pp_WUI
    printf(", pp_WUI");
#endif
#ifdef USE_ZLIB
    printf(", USE_ZLIB");
#endif
#ifdef WIN32
    printf(", WIN32");
#endif
    printf("\n");
  }
}

/* ------------------ checktimebound ------------------------ */

void checktimebound(void){
  int i,j;
  slice *sd;
  mesh *meshi;
  blockagedata *bc;
  particle *parti;

  if(timedrag==0&&itime>ntimes-1||timedrag==1&&itime<0){
    izone=0;itime=0;iframe=iframebeg;
    for(i=0;i<nslice;i++){
      sd=sliceinfo+i;
      sd->islice=0;
    }
    for(i=0;i<selected_case->nmeshes;i++){
      meshi=selected_case->meshinfo+i;
      meshi->ipatch=0;
    }
    for(i=0;i<selected_case->nmeshes;i++){
      meshi=selected_case->meshinfo+i;
      if(meshi->isotimes==NULL)continue;
      meshi->iiso=0;
    }
  }
  if(timedrag==0&&itime<0||timedrag==1&&itime>ntimes-1){
    izone=nzonet-1;itime=ntimes-1;
    for(i=0;i<npartinfo;i++){
      parti=partinfo+i;
      parti->iframe=parti->nframes-1;
    }
    for(i=0;i<nslice;i++){
      sd=sliceinfo+i;
      sd->islice=sd->nsteps-1;
    }
    for(i=0;i<selected_case->nmeshes;i++){
      meshi=selected_case->meshinfo+i;
      meshi->ipatch=meshi->npatch_frames-1;
    }
    for(i=0;i<selected_case->nmeshes;i++){
      meshi=selected_case->meshinfo+i;
      if(meshi->isotimes==NULL)continue;
      meshi->iiso=meshi->nisosteps-1;
    }
  }
  /* set blockage visibility */

  for(i=0;i<selected_case->nmeshes;i++){
    meshi=selected_case->meshinfo+i;
    for(j=0;j<meshi->nbptrs;j++){
      bc=meshi->blockageinfoptrs[j];
      if(bc->showtimelist==NULL)continue;
      bc->show=bc->showtimelist[itime];
    }
  }
}

/* ------------------ getplotstate ------------------------ */

int getplotstate(int choice){
  int i;
  mesh *meshi;
  plot3d *ploti;
  slice *slicei;
  vslice *vslicei;
  patch *patchi;
  particle *parti;
  targ *targi;
  iso *isoi;
  zone *zonei;
  smoke3d *smoke3di;
  tourdata *touri;
  int ii;

  update_loaded_lists();
  switch (choice){
    case STATIC_PLOTS:
    case STATIC_PLOTS_NORECURSE:
      stept = 0;
      for(i=0;i<selected_case->nmeshes;i++){
        meshi=selected_case->meshinfo + i;
        if(meshi->plot3dfilenum==-1)continue;
        ploti = plot3dinfo + meshi->plot3dfilenum;
        if(ploti->loaded==0||ploti->display==0)continue;
        if(meshi->visx==0&&meshi->visy==0&&meshi->visz==0&&visiso==0)continue;
        return STATIC_PLOTS;
      }
      if(choice!=STATIC_PLOTS_NORECURSE){
        return getplotstate(DYNAMIC_PLOTS_NORECURSE);
      }
      break;
    case DYNAMIC_PLOTS:
    case DYNAMIC_PLOTS_NORECURSE:
      for(ii=0;ii<nslice_loaded;ii++){
        i = slice_loaded_list[ii];
        slicei = sliceinfo + i;
        if(slicei->display==0||slicei->type!=islicetype)continue;
        if(slicei->volslice==0&&visGrid==0)stept = 1; 
        return DYNAMIC_PLOTS;
      }
      if(visGrid==0)stept = 1;
      if(visGrid==1)stept = 0;
#ifdef pp_WUI
      if(visTerrain==1){
        for(i=0;i<nterraininfo;i++){
          terraindata *terri;

          terri = terraininfo + i;
          if(terri->loaded==1){
            return DYNAMIC_PLOTS;
          }
        }
      }
#endif
      for(i=0;i<nvslice;i++){
        vslicei = vsliceinfo + i;
        if(vslicei->display==0||vslicei->type!=islicetype)continue;
        return DYNAMIC_PLOTS;
      }
      for(ii=0;ii<npatch_loaded;ii++){
        i = patch_loaded_list[ii];
        patchi = patchinfo + i;
        if(patchi->display==0||patchi->type!=ipatchtype)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<npartinfo;i++){
        parti = partinfo + i;
        if(parti->loaded==0||parti->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<niso;i++){
        isoi = isoinfo + i;
        if(isoi->loaded==0)continue;
        if(isoi->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<nzone;i++){
        zonei = zoneinfo + i;
        if(zonei->loaded==0||zonei->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<ntarg_files;i++){
        targi = targinfo + i;
        if(targi->loaded==0||targi->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<ntours;i++){
        touri = tourinfo + i;
        if(touri->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<nsmoke3d;i++){
        smoke3di = smoke3dinfo + i;
        if(smoke3di->loaded==0||smoke3di->display==0)continue;
        return DYNAMIC_PLOTS;
      }

      if(choice!=DYNAMIC_PLOTS_NORECURSE)return getplotstate(STATIC_PLOTS_NORECURSE);
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
  stept = 0;
  return NO_PLOTS;
}

/* ------------------ ExtractFrustum ------------------------ */

void ExtractFrustum(void){

/* code from:  http://www.markmorley.com/opengl/frustumculling.html */
   float   proj[16];
   float   modl[16];
   float   clip[16];
   float   t;

   /* Get the current PROJECTION matrix from OpenGL */
   glGetFloatv( GL_PROJECTION_MATRIX, proj );

   /* Get the current MODELVIEW matrix from OpenGL */
   glGetFloatv( GL_MODELVIEW_MATRIX, modl );

   /* Combine the two matrices (multiply projection by modelview) */
   clip[ 0] = modl[ 0] * proj[ 0] + modl[ 1] * proj[ 4] + modl[ 2] * proj[ 8] + modl[ 3] * proj[12];
   clip[ 1] = modl[ 0] * proj[ 1] + modl[ 1] * proj[ 5] + modl[ 2] * proj[ 9] + modl[ 3] * proj[13];
   clip[ 2] = modl[ 0] * proj[ 2] + modl[ 1] * proj[ 6] + modl[ 2] * proj[10] + modl[ 3] * proj[14];
   clip[ 3] = modl[ 0] * proj[ 3] + modl[ 1] * proj[ 7] + modl[ 2] * proj[11] + modl[ 3] * proj[15];

   clip[ 4] = modl[ 4] * proj[ 0] + modl[ 5] * proj[ 4] + modl[ 6] * proj[ 8] + modl[ 7] * proj[12];
   clip[ 5] = modl[ 4] * proj[ 1] + modl[ 5] * proj[ 5] + modl[ 6] * proj[ 9] + modl[ 7] * proj[13];
   clip[ 6] = modl[ 4] * proj[ 2] + modl[ 5] * proj[ 6] + modl[ 6] * proj[10] + modl[ 7] * proj[14];
   clip[ 7] = modl[ 4] * proj[ 3] + modl[ 5] * proj[ 7] + modl[ 6] * proj[11] + modl[ 7] * proj[15];

   clip[ 8] = modl[ 8] * proj[ 0] + modl[ 9] * proj[ 4] + modl[10] * proj[ 8] + modl[11] * proj[12];
   clip[ 9] = modl[ 8] * proj[ 1] + modl[ 9] * proj[ 5] + modl[10] * proj[ 9] + modl[11] * proj[13];
   clip[10] = modl[ 8] * proj[ 2] + modl[ 9] * proj[ 6] + modl[10] * proj[10] + modl[11] * proj[14];
   clip[11] = modl[ 8] * proj[ 3] + modl[ 9] * proj[ 7] + modl[10] * proj[11] + modl[11] * proj[15];

   clip[12] = modl[12] * proj[ 0] + modl[13] * proj[ 4] + modl[14] * proj[ 8] + modl[15] * proj[12];
   clip[13] = modl[12] * proj[ 1] + modl[13] * proj[ 5] + modl[14] * proj[ 9] + modl[15] * proj[13];
   clip[14] = modl[12] * proj[ 2] + modl[13] * proj[ 6] + modl[14] * proj[10] + modl[15] * proj[14];
   clip[15] = modl[12] * proj[ 3] + modl[13] * proj[ 7] + modl[14] * proj[11] + modl[15] * proj[15];

   /* Extract the numbers for the RIGHT plane */
   frustum[0][0] = clip[ 3] - clip[ 0];
   frustum[0][1] = clip[ 7] - clip[ 4];
   frustum[0][2] = clip[11] - clip[ 8];
   frustum[0][3] = clip[15] - clip[12];

   /* Normalize the result */
   t = sqrt( frustum[0][0] * frustum[0][0] + frustum[0][1] * frustum[0][1] + frustum[0][2] * frustum[0][2] );
   frustum[0][0] /= t;
   frustum[0][1] /= t;
   frustum[0][2] /= t;
   frustum[0][3] /= t;

   /* Extract the numbers for the LEFT plane */
   frustum[1][0] = clip[ 3] + clip[ 0];
   frustum[1][1] = clip[ 7] + clip[ 4];
   frustum[1][2] = clip[11] + clip[ 8];
   frustum[1][3] = clip[15] + clip[12];

   /* Normalize the result */
   t = sqrt( frustum[1][0] * frustum[1][0] + frustum[1][1] * frustum[1][1] + frustum[1][2] * frustum[1][2] );
   frustum[1][0] /= t;
   frustum[1][1] /= t;
   frustum[1][2] /= t;
   frustum[1][3] /= t;

   /* Extract the BOTTOM plane */
   frustum[2][0] = clip[ 3] + clip[ 1];
   frustum[2][1] = clip[ 7] + clip[ 5];
   frustum[2][2] = clip[11] + clip[ 9];
   frustum[2][3] = clip[15] + clip[13];

   /* Normalize the result */
   t = sqrt( frustum[2][0] * frustum[2][0] + frustum[2][1] * frustum[2][1] + frustum[2][2] * frustum[2][2] );
   frustum[2][0] /= t;
   frustum[2][1] /= t;
   frustum[2][2] /= t;
   frustum[2][3] /= t;

   /* Extract the TOP plane */
   frustum[3][0] = clip[ 3] - clip[ 1];
   frustum[3][1] = clip[ 7] - clip[ 5];
   frustum[3][2] = clip[11] - clip[ 9];
   frustum[3][3] = clip[15] - clip[13];

   /* Normalize the result */
   t = sqrt( frustum[3][0] * frustum[3][0] + frustum[3][1] * frustum[3][1] + frustum[3][2] * frustum[3][2] );
   frustum[3][0] /= t;
   frustum[3][1] /= t;
   frustum[3][2] /= t;
   frustum[3][3] /= t;

   /* Extract the FAR plane */
   frustum[4][0] = clip[ 3] - clip[ 2];
   frustum[4][1] = clip[ 7] - clip[ 6];
   frustum[4][2] = clip[11] - clip[10];
   frustum[4][3] = clip[15] - clip[14];

   /* Normalize the result */
   t = sqrt( frustum[4][0] * frustum[4][0] + frustum[4][1] * frustum[4][1] + frustum[4][2] * frustum[4][2] );
   frustum[4][0] /= t;
   frustum[4][1] /= t;
   frustum[4][2] /= t;
   frustum[4][3] /= t;

   /* Extract the NEAR plane */
   frustum[5][0] = clip[ 3] + clip[ 2];
   frustum[5][1] = clip[ 7] + clip[ 6];
   frustum[5][2] = clip[11] + clip[10];
   frustum[5][3] = clip[15] + clip[14];

   /* Normalize the result */
   t = sqrt( frustum[5][0] * frustum[5][0] + frustum[5][1] * frustum[5][1] + frustum[5][2] * frustum[5][2] );
   frustum[5][0] /= t;
   frustum[5][1] /= t;
   frustum[5][2] /= t;
   frustum[5][3] /= t;
}

/* ------------------ PointInFrustum ------------------------ */

int PointInFrustum( float *xvec){
   int p;
   float x, y, z;
   x=xvec[0];
   y=xvec[1];
   z=xvec[2];

   for( p = 0; p < 6; p++ ){
     if( frustum[p][0]*x + frustum[p][1]*y + frustum[p][2]*z + frustum[p][3] <= 0 ){
         return 0;
     }
   }
   return 1;
}

/* ------------------ RectangleInFrustum ------------------------ */

int RectangleInFrustum( float *x11, float *x12, float *x22, float *x21){
   int p;

   for( p = 0; p < 6; p++ ){
      if( frustum[p][0]*x11[0] + frustum[p][1]*x11[1] + frustum[p][2]*x11[2] + frustum[p][3] > 0 )continue;
      if( frustum[p][0]*x12[0] + frustum[p][1]*x12[1] + frustum[p][2]*x12[2] + frustum[p][3] > 0 )continue;
      if( frustum[p][0]*x22[0] + frustum[p][1]*x22[1] + frustum[p][2]*x22[2] + frustum[p][3] > 0 )continue;
      if( frustum[p][0]*x21[0] + frustum[p][1]*x21[1] + frustum[p][2]*x21[2] + frustum[p][3] > 0 )continue;
      return 0;
   }
   return 1;
}

#include "menus.h"
