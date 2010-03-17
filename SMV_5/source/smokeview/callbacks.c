// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>  
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "contourdefs.h"
#include "flowfiles.h"
#include "smokeviewapi.h"
#include "MALLOC.h"
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

#define KEY_ALT 0
#define KEY_CTRL 1
#define KEY_SHIFT 3
#define KEY_NONE 2

#define TERRAIN_FIRE_LINE_UPDATE 39
void WUI_CB(int var);
char callbacks_revision[]="$Revision$";

#undef pp_GPU_CULL_STATE
#ifdef pp_GPU
#define pp_GPU_CULL_STATE
#endif
#ifdef pp_CULL
#ifndef pp_GPU_CULL_STATE
#define pp_GPU_CULL_STATE
#endif
#endif

#ifdef pp_GPU_CULL_STATE
void print_gpu_cull_state(void);
#endif
void ScriptMenu(int var);

void glui_script_enable(void);
void update_glui_viewlist(void);
void update_glui_cellcenter_interp(void);
float gmod(float x, float y);
void  OBJECT_CB(int flag);

/* ------------------ WindowStatus ------------------------ */

void WindowStatus(int state){
  printf("state=%i\n",state);
  switch (state){
  case GLUT_HIDDEN:
  case GLUT_FULLY_COVERED:
    break;
  case GLUT_FULLY_RETAINED:
  case GLUT_PARTIALLY_RETAINED:
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
}

/* ------------------ mouse_edit_tour ------------------------ */

void mouse_edit_tour(int button, int state, int x, int y){
  int val, val1;
  int mouse_x, mouse_y;
  GLubyte r, g, b;

  mouse_x=x; mouse_y=screenHeight-y;
  glDisable(GL_BLEND);
  glDisable(GL_DITHER);
  glDisable(GL_FOG);
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_TEXTURE_2D);
  glShadeModel(GL_FLAT);

  ShowScene(SELECT,VIEW_CENTER,0,0,0,screenWidth,screenHeight);
  glReadBuffer(GL_BACK);
  glReadPixels(mouse_x,mouse_y,1,1,GL_RED,   GL_UNSIGNED_BYTE, &r);
  glReadPixels(mouse_x,mouse_y,1,1,GL_GREEN, GL_UNSIGNED_BYTE, &g);
  glReadPixels(mouse_x,mouse_y,1,1,GL_BLUE,  GL_UNSIGNED_BYTE, &b);

  r = r>>nredshift;
  g = g>>ngreenshift;
  b = b>>nblueshift;

  val1 = (r << (nbluebits+ngreenbits)) | (g << nbluebits) | b;
  val = val1;
  if(val!=0&&itourknots>=0&&itourknots<ntourknots&&tourknotskeylist!=NULL){
    tourknotskeylist[itourknots]->selected=0;
  }
  if(val>0&&val<=ntourknots){
  
  /* need to start colors at 1 so that black (color 0,0,0) is not interpreted as a blockage */

    val--;
    itourknots=val;
    if(tourknotskeylist!=NULL){
      new_select(tourknotskeylist[itourknots]);
      selected_tour=tourknotstourlist[itourknots];
    }
    else{
      selected_tour=NULL;
      itourknots=-1;
    }
    set_glui_keyframe();
    update_tourcontrols();
  }
  glShadeModel(GL_SMOOTH);
  glEnable(GL_BLEND);
  glEnable(GL_LIGHTING);
}

/* ------------------ mouse_edit_blockage ------------------------ */

void mouse_edit_blockage(int button, int state, int x, int y){
  int val, val1;
  int mouse_x, mouse_y;
  GLubyte r, g, b;
  int i;
  mesh *meshi;
  selectdata *sd;

  mouse_x=x; mouse_y=screenHeight-y;
  glDisable(GL_BLEND);
  glDisable(GL_DITHER);
  glDisable(GL_FOG);
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_TEXTURE_2D);
  glShadeModel(GL_FLAT);

  ShowScene(SELECT,VIEW_CENTER,0,0,0,screenWidth,screenHeight);
  glReadBuffer(GL_BACK);
  glReadPixels(mouse_x,mouse_y,1,1,GL_RED,   GL_UNSIGNED_BYTE, &r);
  glReadPixels(mouse_x,mouse_y,1,1,GL_GREEN, GL_UNSIGNED_BYTE, &g);
  glReadPixels(mouse_x,mouse_y,1,1,GL_BLUE,  GL_UNSIGNED_BYTE, &b);

  r = r>>nredshift;
  g = g>>ngreenshift;
  b = b>>nblueshift;

  val1 = (r << (nbluebits+ngreenbits)) | (g << nbluebits) | b;
  val = val1;
  
  if(val>0&&val<=ntotalfaces){
      /* need to start colors at 1 so that black (color 0,0,0) is not
                interpreted as a blockage */
    val--;
    sd = selectfaceinfo + val;
    highlight_block=sd->blockage;
    highlight_mesh=sd->mesh;
    meshi = meshinfo + highlight_mesh;
    update_current_mesh(meshi);
    bchighlight_old=bchighlight;
    bchighlight = meshi->blockageinfoptrs[highlight_block];
    for(i=0;i<6;i++){
      surface_indices[i]=inv_sorted_surfidlist[bchighlight->surf_index[i]];
      surface_indices_bak[i]=inv_sorted_surfidlist[bchighlight->surf_index[i]];
    }

    glShadeModel(GL_SMOOTH);
    glEnable(GL_BLEND);
    glEnable(GL_LIGHTING);

    switch (sd->dir){
      case DOWN_X:
      case UP_X:
        xyz_dir=0;
        break;
      case DOWN_Y:
      case UP_Y:
        xyz_dir=1;
        break;
      case DOWN_Z:
      case UP_Z:
        xyz_dir=2;
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
    switch (sd->dir){
      case DOWN_X:
      case DOWN_Y:
      case DOWN_Z:
        which_face=0;
        break;
      case UP_X:
      case UP_Y:
      case UP_Z:
        which_face=1;
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
    update_blockvals(1);
  }
}

/* ------------------ mouse_select_device ------------------------ */

void mouse_select_device(int button, int state, int x, int y){
  int val;
  int mouse_x, mouse_y;
  GLubyte r, g, b;

  mouse_x=x; mouse_y=screenHeight-y;
  glDisable(GL_BLEND);
  glDisable(GL_DITHER);
  glDisable(GL_FOG);
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_TEXTURE_2D);
  glShadeModel(GL_FLAT);

  ShowScene(SELECT,VIEW_CENTER,0,0,0,screenWidth,screenHeight);
  glReadBuffer(GL_BACK);
  glReadPixels(mouse_x,mouse_y,1,1,GL_RED,   GL_UNSIGNED_BYTE, &r);
  glReadPixels(mouse_x,mouse_y,1,1,GL_GREEN, GL_UNSIGNED_BYTE, &g);
  glReadPixels(mouse_x,mouse_y,1,1,GL_BLUE,  GL_UNSIGNED_BYTE, &b);

  r = r>>nredshift;
  g = g>>ngreenshift;
  b = b>>nblueshift;

  val = (r << (nbluebits+ngreenbits)) | (g << nbluebits) | b;
  
  if(val>0){
    device *devicei;
    float *xyz;

    selected_device_tag=val;
    devicei = deviceinfo + val-1;
    xyz = devicei->xyz;

    if(devicei->labelptr!=NULL&&strcmp(devicei->labelptr,"null")!=0){
      printf("Selected Device: index=%i location:(%f,%f,%f) label:%s\n",val,xyz[0],xyz[1],xyz[2],devicei->labelptr);
    }
    else{
      printf("Selected Device: index=%i location:(%f,%f,%f)\n",val,xyz[0],xyz[1],xyz[2]);
    }
    glShadeModel(GL_SMOOTH);
    glEnable(GL_BLEND);
    glEnable(GL_LIGHTING);
  }
}

/* ------------------ select_avatar ------------------------ */

void mouse_select_avatar(int button, int state, int x, int y){
  int val;
  int mouse_x, mouse_y;
  GLubyte r, g, b;

  mouse_x=x; mouse_y=screenHeight-y;
  glDisable(GL_BLEND);
  glDisable(GL_DITHER);
  glDisable(GL_FOG);
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_TEXTURE_2D);
  glShadeModel(GL_FLAT);

  ShowScene(SELECT,VIEW_CENTER,0,0,0,screenWidth,screenHeight);
  glReadBuffer(GL_BACK);
  glReadPixels(mouse_x,mouse_y,1,1,GL_RED,   GL_UNSIGNED_BYTE, &r);
  glReadPixels(mouse_x,mouse_y,1,1,GL_GREEN, GL_UNSIGNED_BYTE, &g);
  glReadPixels(mouse_x,mouse_y,1,1,GL_BLUE,  GL_UNSIGNED_BYTE, &b);

  r = r>>nredshift;
  g = g>>ngreenshift;
  b = b>>nblueshift;

  val = (r << (nbluebits+ngreenbits)) | (g << nbluebits) | b;
  
  if(val>0){
    selected_avatar_tag=val;
    glShadeModel(GL_SMOOTH);
    glEnable(GL_BLEND);
    glEnable(GL_LIGHTING);
  }
}

/* ------------------ mouse ------------------------ */

void mouse(int button, int state, int x, int y){
  int mouse_x, mouse_y;
  GLubyte r, g, b;
  int val1, val;
  selectdata *sd;
  mesh *meshi;
  int i;
  float *eye_xyz;

  if(trainer_mode==1){
    update_glui_viewlist();
  }
  eye_xyz = camera_current->eye;
  if(selected_view!=-999){
    selected_view=-999;
    updatemenu=1;
  }
  glui_move_mode=-1;
  GLUTPOSTREDISPLAY
  if(state==GLUT_UP){
    eye_xyz0[0]=eye_xyz[0];
    eye_xyz0[1]=eye_xyz[1];
    eye_xyz0[2]=eye_xyz[2];
    update_translate();
    timedrag=0;
    colordrag=0;
    colorsplitdrag=0;
    glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
    update_trainer_moves();
    return;
  }
  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN){
    glutSetCursor(GLUT_CURSOR_INFO);

    /* edit blockages */

    if(blockageSelect==1){
      mouse_edit_blockage(button,state,x,y);
    }

    /* edit tours */

    if(edittour==1&&blockageSelect==0){
      mouse_edit_tour(button, state, x, y);
    }

    if(select_avatar==1){
      mouse_select_avatar(button, state, x, y);
    }

    if(select_device==1){
      mouse_select_device(button, state, x, y);
    }
    GLUTPOSTREDISPLAY
    if( showtime==1 || showplot3d==1){
      int temp;
      int yy;
      float factor;
      int ifactor;

      temp = (int)(1.2*dwinH);
      if(x>screenWidth-dwinWW){
        yy = screenHeight - y;
        factor=((float)(yy-temp)/(screenHeight-temp))*((nrgb+(float)1.0)/(nrgb-(float)0.5));
        if(screenHeight>screenWidth)factor *= (float)screenHeight/screenWidth;
        ifactor=(int)(255*factor);
        if(ifactor>=0&&ifactor<256){
          int valmax=255;
          int valmin=0;

          if(ifactor>valmax)ifactor=valmax;
          if(ifactor<valmin)ifactor=valmin;
        }
        else{
          ifactor=-1;
        }
        colorbar_select_index=ifactor;
        state=glutGetModifiers();
        if(state==GLUT_ACTIVE_CTRL&&(showtime==1 || showplot3d==1)&&current_colorbar!=NULL&&current_colorbar->nsplits==1){
          colorsplitdrag=1;
        }
        else{
          colordrag=1;
          updatecolors(ifactor);
        }
        return;
      }
    }
    if(screenHeight-y<50&&ntimes>0&&visTimeLabels==1&&showtime==1){
      float xleft;

      if(fontindex==LARGE_FONT){
        xleft=xtimeleft+0.11;
      }
      else{
        xleft=xtimeleft;
      }
      itime=(int)((xtemp*x/((screenWidth-dwinWW))-xleft)*(ntimes-1)/(xtimeright-xleft));
      checktimebound();
      timedrag=1;
      stept=0;
      IDLE();

      return;
    }
    copy_camera(camera_last,camera_current);
    if(canrestorelastview==0){
      updatemenu=1;
      canrestorelastview=1;
      enable_reset_saved_view();
    }
    state=glutGetModifiers();
    switch (state){
    case GLUT_ACTIVE_CTRL:
      key_state = KEY_CTRL;
      eye_xyz0[0]=eye_xyz[0];
      eye_xyz0[1]=eye_xyz[1];
      touring=0;
      break;
    case GLUT_ACTIVE_ALT:
      key_state = KEY_ALT;
      eye_xyz0[0]=eye_xyz[0];
      eye_xyz0[2]=eye_xyz[2];
      touring=0;
      break;
    case GLUT_ACTIVE_SHIFT:
    default:
      key_state = KEY_NONE;
      start_xyz0[0]=x;
      start_xyz0[1]=y;
      touring=0;
      break;
    }
    xm0=x; 
    ym0=y;
  }
  GLUTPOSTREDISPLAY
  if(blockageSelect==1){
    Display();
  }

}

/* ------------------ motion ------------------------ */

void motion(int xm, int ym){
  float xx, yy;
  float *eye_xyz, *angle_zx;
  float direction_angle;
  float elevation_angle;

  eye_xyz = camera_current->eye;
  angle_zx = camera_current->angle_zx;
  
  reset_glui_view(-1);

  GLUTPOSTREDISPLAY
  if( colordrag==1&&(showtime==1 || showplot3d==1)){
    int temp;
    int ifactor;
    float factor;
    int valmax=255;
    int valmin=0;

    temp = (int)(1.2*dwinH);
    if(xm>screenWidth-dwinWW){
      yy = screenHeight - ym;
      factor=(yy-temp)/(screenHeight-temp);
      factor *= (nrgb+1.0)/(nrgb-0.5);
      if(screenHeight>screenWidth)factor *= (float)screenHeight/screenWidth;
      ifactor=(int)(255*factor);
      if(ifactor<256||ifactor>=0){
        if(ifactor>valmax)ifactor=valmax;
        if(ifactor<valmin)ifactor=valmin;
      }
      else{
        ifactor=-1;
      }
      colorbar_select_index=ifactor;
      updatecolors(ifactor);
    }
    return;
  }
  if(colorsplitdrag==1&&(showtime==1 || showplot3d==1)&&current_colorbar!=NULL&&current_colorbar->nsplits==1){
    int temp;
    int ifactor;
    float factor;
    int valmax=255;
    int valmin=0;

    temp = (int)(1.2*dwinH);
    if(xm>screenWidth-dwinWW){
      int ii;
      unsigned char *cc;

      yy = screenHeight - ym;
      factor=(yy-temp)/(screenHeight-temp);
      factor *= (nrgb+1.0)/(nrgb-0.5);
      if(screenHeight>screenWidth)factor *= (float)screenHeight/screenWidth;
      ifactor=(int)(255*factor);

      cc=current_colorbar->index_node;

      if(ifactor>250)ifactor=250;
      if(ifactor<5)ifactor=5;
      ii=current_colorbar->splits[0];
      current_colorbar->index_node[ii]=ifactor;
      current_colorbar->index_node[ii-1]=ifactor;
      remapcolorbar(current_colorbar);
      updatecolors(-1);
      update_colorbar_splits(current_colorbar);
    }
    return;
  }
  if(timedrag==1){
	float xxleft;

	  xxleft = xtimeleft;
    if(fontindex==LARGE_FONT)xxleft=xtimeleft+0.11;
    if(screenHeight-ym<50&&ntimes>0&&visTimeLabels==1&&showtime==1){
      itime=(int)((xtemp*xm/((screenWidth-dwinWW))-xxleft)*(ntimes-1)/(xtimeright-xxleft));
      checktimebound();
      timedrag=1;
    }
  IDLE();

    return;
  }
  screenWidth2 = screenWidth - dwinWW;
  screenHeight2 = screenHeight - dwinH;

  switch (key_state){
    case KEY_NONE:
      switch (eyeview){
        case WORLD_CENTERED:
        case WORLD_CENTERED_LEVEL:
          angle_zx[0] += (xm - start_xyz0[0]);
          if(eyeview==WORLD_CENTERED){
            angle_zx[1] += (ym - start_xyz0[1]);
          }
          else{
            angle_zx[1]=0.0;
          }
          start_xyz0[0]=xm;
          start_xyz0[1]=ym;
          break;
        case EYE_CENTERED:
#define ANGLE_FACTOR 0.25
          camera_current->direction_angle += (xm - start_xyz0[0])*ANGLE_FACTOR;
          direction_angle=camera_current->direction_angle;
          camera_current->cos_direction_angle = cos(PI*direction_angle/180.0);
          camera_current->sin_direction_angle = sin(PI*direction_angle/180.0);
          start_xyz0[0]=xm;

          camera_current->elevation_angle -= (ym - start_xyz0[1])*ANGLE_FACTOR;
          elevation_angle=camera_current->elevation_angle;
          if(elevation_angle>80.0)elevation_angle=80.0;
          if(elevation_angle<-80.0)elevation_angle=-80.0;
          camera_current->elevation_angle=elevation_angle;
          camera_current->cos_elevation_angle = cos(PI*elevation_angle/180.0);
          camera_current->sin_elevation_angle = sin(PI*elevation_angle/180.0);
          start_xyz0[1]=ym;
          break;
        default:
          ASSERT(FFALSE);
          break;
      }
      break;

    case KEY_CTRL:
      {
        float dx, dy;

        xx = xm-xm0;
        xx = xx/(float)screenWidth2;
        yy = ym-ym0;
        yy = yy/(float)screenHeight2;
        if(eyeview==EYE_CENTERED){
          float xx2, yy2;

          xx2 = camera_current->cos_direction_angle*xx - camera_current->sin_direction_angle*yy;
          yy2 = camera_current->sin_direction_angle*xx + camera_current->cos_direction_angle*yy;
          xx = xx2;
          yy = yy2;
        }
        else{
          dx = (xyzbox+eye_xyz0[0])*xx;
          dy = -(xyzbox-eye_xyz0[1])*yy;
          eye_xyz[0] = eye_xyz0[0] + dx;
          eye_xyz[1] = eye_xyz0[1] + dy;
          eye_xyz0[0]=eye_xyz[0];
          eye_xyz0[1]=eye_xyz[1];
          xm0=xm;
          ym0=ym;
        }
      }
      break;

    case KEY_ALT:

      xx = xm-xm0;
      xx = xx/(float)screenWidth2;
      yy = ym-ym0;
      yy = yy/(float)screenHeight2;

      eye_xyz[0] = eye_xyz0[0]; /* disable horizontal motion */
      eye_xyz[2] = eye_xyz0[2] - 4*(xyzbox-eye_xyz0[2])*yy;
      viewx = eye_xyz[0] - delx;
      viewz = eye_xyz[2] - delz;
      break;
    case KEY_SHIFT:
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
}

/* ------------------ keyboard_up ------------------------ */

void keyboard_up(unsigned char key, int x, int y){
  resetclock=1;
}

/* ------------------ get_vecfactor ------------------------ */

float get_vecfactor(int *ivec){
  float vec;
  if(*ivec>NVECLENGTHS-1)*ivec=0;
  if(*ivec<0)*ivec=NVECLENGTHS-1;
  vec=VECFRACTION*veclengths[*ivec];
  return vec;
}
/* ------------------ keyboard ------------------------ */

void keyboard(unsigned char key, int x, int y){
  char key2;
  int skip2;
  char one='1';
  mesh *gbsave,*gbi;
  int i;

  GLUTPOSTREDISPLAY
  updatemenu=1;
  key2 = (char)key;
  if(key2!='H'&&key2!='N'&&key2!='R'&&key2!='P'&&key2!='T'&&key2!='G'&&key2!='S'&&key2!='M'
#ifdef pp_CULL
    &&key2!='C'
#endif
#ifdef pp_LIGHT
    &&key2!='L'
#endif
    &&isupper(key2))key2=tolower(key2); /* map upper case characters to lower */

#ifdef pp_LIGHT
  if(strncmp((const char *)&key2,"L",1)==0){
    show_smokelighting = 1 - show_smokelighting;
    update_showlight();
    updatemenu=1;
    return;
  }
#endif
#ifdef _DEBUG 
  if(nsmoke3d_files>0&&strncmp((const char *)&key2,"l",1)==0){
    smokecullflag=1-smokecullflag;
    if(smokecullflag==0){
      smokedrawtest=1-smokedrawtest;
    }
    printf("smokecullflag=%i\n smokedrawtest=%i\n",smokecullflag,smokedrawtest);
    update_smoke3dflags();
    return;
  }
  if(nsmoke3d_files>0&&strncmp((const char *)&key2,"n",1)==0){
    adjustalphaflag++;
    if(adjustalphaflag>3)adjustalphaflag=0;
    printf("adjustalphaflag=%i\n",adjustalphaflag);
    update_smoke3dflags();
    return;
  }
#endif
  if(strncmp((const char *)&key2,"u",1)==0){
    int state;

    state=glutGetModifiers();

    switch (state){
      case GLUT_ACTIVE_ALT:
        skip_slice_in_embedded_mesh = 1 - skip_slice_in_embedded_mesh;
        update_glui_cellcenter_interp();
        break;
      default:
        ReloadMenu(0);
        break;
    }
    return;
  }
  if(strncmp((const char *)&key2,"q",1)==0){
    blocklocation++;
    if(blocklocation>BLOCKlocation_cad||
       blocklocation>BLOCKlocation_exact&&ncadgeom==0){
       blocklocation=BLOCKlocation_grid;
    }
    if(showedit==1){
      if(blocklocation==BLOCKlocation_exact){
        blockage_as_input=1;
      }
      else{
        blockage_as_input=0;
      }
      OBJECT_CB(BLOCKAGE_AS_INPUT2);
    }
    return;
  }

  if(strncmp((const char *)&key,"#",1)==0){
    writeini(LOCAL_INI);
    updatemenu=1;
    return;
  }
  
  if(strncmp((const char *)&key2,"!",1)==0){
    snap_view_angles();
    return;
  }
  if(strncmp((const char *)&key,"$",1)==0){
    trainer_active=1-trainer_active;
    if(trainer_active==1){
      printf("Trainer mode active\n");
      trainer_mode=1;
      show_trainer();
    }
    if(trainer_active==0){
      printf("Trainer mode inactive\n");
      trainer_mode=0;
      hide_trainer();
    }
    updatemenu=1;
    return;
  }
  if(strncmp((const char *)&key2,"S",1)==0){
    showstereoOLD=showstereo;
    showstereo++;
    if(showstereo>5)showstereo=0;
    if(showstereo==1&&videoSTEREO!=1)showstereo=2;
    update_glui_stereo();
    return;
  }
  if(strncmp((const char *)&key2,"T",1)==0){
#ifdef pp_FRACTILE  
    usetexturebar++;
    if(usetexturebar==4)usetexturebar=0;
#else
    usetexturebar=1-usetexturebar;
#endif
    printf("usetexturebar=%i\n",usetexturebar);
    return;
  }
#ifdef pp_CULL
  if(strncmp((const char *)&key2,"C",1)==0){
    if(nsmoke3d_files>0&&cullactive==1){
      cullsmoke=1-cullsmoke;
      update_smoke3dflags();
      initcull(cullsmoke);
      print_gpu_cull_state();
    }
    return;    
  }
#endif
#ifdef pp_GPU
  if(nsmoke3d_files>0&&strncmp((const char *)&key2,"G",1)==0){
    if(gpuactive==1)usegpu=1-usegpu;
    update_smoke3dflags();
    print_gpu_cull_state();
    return;    
  }
#endif
  if(strncmp((const char *)&key2,"t",1)==0){
    int state;

    state=glutGetModifiers();

    switch (state){
    case GLUT_ACTIVE_ALT:
      DialogMenu(21); // display dialog
      break;
    case GLUT_ACTIVE_CTRL:
    case GLUT_ACTIVE_SHIFT:
    default:
      stept++;
      if(stept>1)stept=0;
      if(stept==1){
        plotstate=getplotstate(DYNAMIC_PLOTS);
        if(plotstate==DYNAMIC_PLOTS){
          reset_gltime();
        }
        else{
          stept=0;
        }
      }
    }
    return;
  }
  if(strncmp((const char *)&key2,"0",1)==0&&plotstate==DYNAMIC_PLOTS){
    updatetimes();
    reset_time_flag=1;
    return;
  }

  if(strncmp((const char *)&key2,"x",1)==0){
    int state;

    state=glutGetModifiers();

    switch (state){

    case GLUT_ACTIVE_ALT:
      DialogMenu(-2); // close all dialogs
      break;
    case GLUT_ACTIVE_CTRL:
      updateshowstep(1-current_mesh->visx,DIRX);
      break;
    case GLUT_ACTIVE_SHIFT:
      current_mesh->visx2 = 1-current_mesh->visx2;
      updateshowstep(1-current_mesh->visx,DIRX);
      break;
    default:
      for(i=0;i<nmeshes;i++){
        mesh *meshi;
        meshi = meshinfo + i;
        meshi->visx2 = 1 - current_mesh->visx;
      }
      updateshowstep(1-current_mesh->visx,DIRX);
    }
    return;
  }
  if(strncmp((const char *)&key2,"y",1)==0){
    int state;

    state=glutGetModifiers();

    switch (state){

    case GLUT_ACTIVE_ALT:
      cellcenter_interp = 1 - cellcenter_interp;
      update_glui_cellcenter_interp();
      break;
     default:
      updateshowstep(1-current_mesh->visy,DIRY);
    }
    return;
  }
  if(strncmp((const char *)&key2,"z",1)==0){
    int state;

    state=glutGetModifiers();

    switch (state){
    case GLUT_ACTIVE_ALT:
      DialogMenu(24); // compress dialog
      break;
    case GLUT_ACTIVE_CTRL:
    case GLUT_ACTIVE_SHIFT:
    default:
      updateshowstep(1-current_mesh->visz,DIRZ);
    }
    return;
  }
  if(strncmp((const char *)&key2,"k",1)==0){
    visTimeLabels = 1 - visTimeLabels;
    if(visTimeLabels==0)printf("Time bar hidden\n");
    if(visTimeLabels==1)printf("Time bar visible\n");
    updatemenu=1;
    return;
  }
  if(strncmp((const char *)&key2,"e",1)==0){
    int state;

    state=glutGetModifiers();

    switch (state){
    case GLUT_ACTIVE_ALT:
      DialogMenu(16); // edit geometry
      break;
    case GLUT_ACTIVE_CTRL:
    case GLUT_ACTIVE_SHIFT:
    default:
      eyeview++;
      if(eyeview>2)eyeview=0;
      handle_eyeview(0);
    }
    return;
  }
  if(strncmp((const char *)&key2,"h",1)==0){
    if(titlesafe_offset==0){
      titlesafe_offset=titlesafe_offsetBASE;
    }
    else{
      titlesafe_offset=0;
    }
    return;
  }
  if(strncmp((const char *)&key2,"f",1)==0){
    int state;

    state=glutGetModifiers();

    switch (state){
    case GLUT_ACTIVE_ALT:
      DialogMenu(14); // display dialog
      break;
    case GLUT_ACTIVE_CTRL:
    case GLUT_ACTIVE_SHIFT:
    default:
      pass_through=1-pass_through;
      update_blockpath();
    }
    return;
  }
  if(strncmp((const char *)&key2,"]",1)==0){
    int state;

    state=glutGetModifiers();
    if(state==GLUT_ACTIVE_ALT){
      printf("re-attaching menus to right mouse button\n");
      glutDetachMenu(GLUT_RIGHT_BUTTON);
      InitMenus(LOAD);
      glutAttachMenu(GLUT_RIGHT_BUTTON);
    }
    return;
  }

  if(strncmp((const char *)&key2,"j",1)==0&&eyeview==EYE_CENTERED){
    int state;

    state=glutGetModifiers();

    switch (state){
    case GLUT_ACTIVE_ALT:
    case GLUT_ACTIVE_CTRL:
    case GLUT_ACTIVE_SHIFT:
    default:
      eyeview_level = 1 - eyeview_level;
      update_glui_speed();
    }
    return;
  }
  if(strncmp((const char *)&key2,"g",1)==0){
    if(ntotal_blockages>0||isZoneFireModel==0){
      togglegridstate(1-visGrid);
    }
    return;
  }
  if(strncmp((const char *)&key2,"o",1)==0){
    highlight_flag++;
    if(highlight_flag>2&&noutlineinfo>0)highlight_flag=0;
    if(highlight_flag>1&&noutlineinfo==0)highlight_flag=0;
    printf("outline mode=%i\n",highlight_flag);
    return;
  }
  if(strncmp((const char *)&key2,"M",1)==0){
    updatemenu=1;
    return;
  }
  if(strncmp((const char *)&key2,"m",1)==0){
    int state;

    state=glutGetModifiers();

    switch (state){
    case GLUT_ACTIVE_ALT:
      DialogMenu(15); // display dialog
      break;
    case GLUT_ACTIVE_CTRL:
    case GLUT_ACTIVE_SHIFT:
    default:
      if(nmeshes>1){
        highlight_mesh++;
        if(highlight_mesh>nmeshes-1)highlight_mesh=0;
        update_current_mesh(meshinfo+highlight_mesh);
      }
    }
    return;
  }
  if(strncmp((const char *)&key2,"c",1)==0){
    int state;

    state=glutGetModifiers();
    switch (state){
    case GLUT_ACTIVE_ALT:
      DialogMenu(18); // clip dialog
      break;
    case GLUT_ACTIVE_CTRL:
    case GLUT_ACTIVE_SHIFT:
    default:
      p3cont2d++;
#ifdef pp_LINE
      if(p3cont2d>2)p3cont2d=0;
#else
      if(p3cont2d>1)p3cont2d=0;
#endif
      update_plot3d_display();
    }
    return;
  }
  if(strncmp((const char *)&key2,"w",1)==0){
    int state;

    state=glutGetModifiers();
    switch (state){
      case GLUT_ACTIVE_ALT:
        DialogMenu(26); // WUI dialog
        break;
      case GLUT_ACTIVE_CTRL:
      case GLUT_ACTIVE_SHIFT:
      default:
      if(eyeview==EYE_CENTERED){
        handle_move_keys(GLUT_KEY_UP);
      }
      else{
        xyz_clipplane++;
        if(xyz_clipplane>2)xyz_clipplane=0;
        update_clip_all();
      }
      break;
    }
    return;
  }
  if(strncmp((const char *)&key2,"a",1)==0&&(eyeview==EYE_CENTERED||(visVector==1&&ReadPlot3dFile==1)||showvslice==1)){
    if(eyeview==EYE_CENTERED){
      handle_move_keys(256+key2);
    }
    else{
      iveclengths += FlowDir;
      vecfactor = get_vecfactor(&iveclengths);
      printf("iveclengths=%i\n",iveclengths);
      if(visVector==1&&ReadPlot3dFile==1){
        gbsave=current_mesh;
        for(i=0;i<nmeshes;i++){
          gbi = meshinfo + i;
          if(gbi->plot3dfilenum==-1)continue;
          update_current_mesh(gbi);
          updateplotslice(1);
          updateplotslice(2);
          updateplotslice(3);
        }
        update_current_mesh(gbsave);
      }
    }
    return;
  }
  if(strncmp((const char *)&key2,"s",1)==0){
    int state;
    state=glutGetModifiers();

    switch (state){
    case GLUT_ACTIVE_ALT:
      DialogMenu(20); // display dialog
      break;
    case GLUT_ACTIVE_CTRL:
      snap_view_angles();
      break;
    case GLUT_ACTIVE_SHIFT:
    default:
      if(eyeview==EYE_CENTERED){
        handle_move_keys(GLUT_KEY_DOWN);
      }
      else{
        vectorskip++;
        if(vectorskip>4)vectorskip=1;
      }
    }
    return;
  }
  if(strncmp((const char *)&key2,"d",1)==0){
    int state;

    state=glutGetModifiers();

    switch (state){
    case GLUT_ACTIVE_ALT:
      DialogMenu(22); // display dialog
      break;
    case GLUT_ACTIVE_CTRL:
    case GLUT_ACTIVE_SHIFT:
    default:
      if(eyeview==EYE_CENTERED){
        handle_move_keys(256+key2);
      }
      else{
        demo_mode++;
        if(demo_mode>5)demo_mode=0;
      }
        break;
    }
    return;
  }
  if(strncmp((const char *)&key2,"&",1)==0){
    antialiasflag=1-antialiasflag;
    printf("antialiasflag=%i\n",antialiasflag);
    return;
  }
  if(strncmp((const char *)&key2,"i",1)==0&&unload_qdata==0){
    handleiso();
    return;
  }
  if(strncmp((const char *)&key2,"b",1)==0&&visiso==1&&unload_qdata==0){
    isooffset+=FlowDir;
    if(isooffset<1)isooffset=offsetmax;
    if(isooffset>offsetmax)isooffset=1;
    updatesurface();
    return;
  }
  if(strncmp((const char *)&key2,"v",1)==0){
    int state;
    
    state=glutGetModifiers();
    switch (state){
      case GLUT_ACTIVE_ALT:
        projection_type = 1 - projection_type;
        TRANSLATE_CB(PROJECTION);
        updatemenu=1;
        break;
      default:
        visVector=1-visVector;
        if(vectorspresent==0)visVector=0;
        updateglui();
        break;
    }
    return;
  }
  if(strncmp((const char *)&key2,"H",1)==0){
    int nslice_loaded=0, nvslice_loaded=0;

    for(i=0;i<nslice_files;i++){
      slice *sd;

      sd = sliceinfo + i;
      if(sd->loaded==1)nslice_loaded++;
    }
    for(i=0;i<nvslice;i++){
      vslice *vd;

      vd = vsliceinfo + i;
      if(vd->loaded==1)nvslice_loaded++;
    }
    stept=1;
    if(nvslice_loaded>0){
      if(show_all_slices==0){
        ShowVSliceMenu(SHOW_ALL);
        force_redisplay=1;
      }
      else{
        itime_save=itime;
        ShowVSliceMenu(HIDE_ALL);
      }
    }
    if(nvslice_loaded==0&&nslice_loaded>0){
      if(show_all_slices==0){
        ShowHideSliceMenu(SHOW_ALL);
        force_redisplay=1;
      }
      else{
        itime_save=itime;
        ShowHideSliceMenu(HIDE_ALL);
      }
    }
    return;
  }
  if(strncmp((const char *)&key2,"P",1)==0){
    cursorPlot3D=1-cursorPlot3D;
    update_cursor_checkbox();
    return;
  }
  if(strncmp((const char *)&key2,"p",1)==0){
    plotn += FlowDir;
    if(plotn<1){
      plotn=numplot3dvars;
    }
    if(plotn>numplot3dvars){
      plotn=1;
    }
    updateallplotslices();
    if(visiso==1&&unload_qdata==0)updatesurface();
    updateplot3dlistindex();
    return;
  }
  if(strncmp((const char *)&key2,"r",1)==0
    ||strncmp((const char *)&key2,"R",1)==0
    ){
#ifdef pp_RENDER
    int rflag=0;

    if(strncmp((const char *)&key2,"R",1)==0){
      render_double=2;
      rflag=1;
    }
    else{
      render_double=0;
    }
    if(scriptoutstream!=NULL){
      if(ntimes>0){
        float timeval;

        timeval=times[itime];
        fprintf(scriptoutstream,"SETTIMEVAL\n");
        fprintf(scriptoutstream," %f\n",timeval);
      }
      else{
        int show_plot3dkeywords=0;

        for(i=0;i<nmeshes;i++){
          mesh *meshi;
          plot3d *plot3di;
          float *xp, *yp, *zp;

          meshi = meshinfo  + i;
          if(meshi->plot3dfilenum==-1)continue;

          plot3di = plot3dinfo + meshi->plot3dfilenum;
          if(plot3di->display==0)continue;
          show_plot3dkeywords=1;
          xp = meshi->xplt_orig;
          yp = meshi->yplt_orig;
          zp = meshi->zplt_orig;
          fprintf(scriptoutstream,"SHOWPLOT3DDATA\n");
          fprintf(scriptoutstream," %i %i %i %i %f\n",i+1,1, plotn,meshi->visx,xp[meshi->plotx]);
          fprintf(scriptoutstream,"SHOWPLOT3DDATA\n");
          fprintf(scriptoutstream," %i %i %i %i %f\n",i+1,2, plotn,meshi->visy,yp[meshi->ploty]);
          fprintf(scriptoutstream,"SHOWPLOT3DDATA\n");
          fprintf(scriptoutstream," %i %i %i %i %f\n",i+1,3, plotn,meshi->visz,zp[meshi->plotz]);
          fprintf(scriptoutstream,"SHOWPLOT3DDATA\n");
          fprintf(scriptoutstream," %i %i %i %i %i\n",i+1,4, plotn,visiso,plotiso[plotn-1]);

        }
        if(show_plot3dkeywords==1){
          fprintf(scriptoutstream,"PLOT3DPROPS\n");
          fprintf(scriptoutstream," %i %i %i %i\n",plotn,visVector,iveclengths,p3cont2d);
        }
      }
      if(rflag==0){
        fprintf(scriptoutstream,"RENDERONCE\n");
      }
      else{
        fprintf(scriptoutstream,"RENDERDOUBLEONCE\n");
      }
      fprintf(scriptoutstream," %s\n",script_renderfile);
    }
    RenderOnceNow=1;
    if(showstereo!=0){
      RenderOnceNowL=1;
      RenderOnceNowR=1;
    }
    RenderState(1);
#else
    printf("*** rendering not supported in this version of Smokeview\n");
#endif
    return;
  }

  skip2=key2-one+1;
  if(skip2>0&&skip2<10)skip=skip2;

  /* if not a directional key then return */

  if(strncmp((const char *)&key2,"-",1)!=0&&strncmp((const char *)&key2," ",1)!=0&&
     strncmp((const char *)&key2,"+",1)!=0&&strncmp((const char *)&key2,"=",1)!=0&&
     strncmp((const char *)&key2,"<",1)!=0&&strncmp((const char *)&key2,">",1)!=0&&
     strncmp((const char *)&key2,",",1)!=0&&strncmp((const char *)&key2,".",1)!=0&&
     strncmp((const char *)&key2,"_",1)!=0)return;

  if(xyz_clipplane!=0&&(
    strncmp((const char *)&key2,"<",1)==0||strncmp((const char *)&key2,",",1)==0||
    strncmp((const char *)&key2,">",1)==0||strncmp((const char *)&key2,".",1)==0)){

    if(strncmp((const char *)&key2,"<",1)==0||strncmp((const char *)&key2,",",1)==0){ClipDir=-1;}
     else if(strncmp((const char *)&key2,">",1)==0||strncmp((const char *)&key2,".",1)==0){ClipDir=1;}

    if(stepclip_x==1  )clip_i += skip*ClipDir;
    if(stepclip_y==1  )clip_j += skip*ClipDir;
    if(stepclip_z==1  )clip_k += skip*ClipDir;
    if(stepclip_X==1  )clip_I += skip*ClipDir;
    if(stepclip_Y==1  )clip_J += skip*ClipDir;
    if(stepclip_Z==1  )clip_K += skip*ClipDir;

    updateclipbounds(clip_x,&clip_i,clip_X,&clip_I,current_mesh->ibar);
    updateclipbounds(clip_y,&clip_j,clip_Y,&clip_J,current_mesh->jbar);
    updateclipbounds(clip_z,&clip_k,clip_Z,&clip_K,current_mesh->kbar);
    return;
  }

  if(strncmp((const char *)&key2,"-",1)==0||strncmp((const char *)&key2,"_",1)==0){
    FlowDir=-1;
    updatemenu=0;
  }
   else if(strncmp((const char *)&key2," ",1)==0||
     strncmp((const char *)&key2,"=",1)==0||
     strncmp((const char *)&key2,"+",1)==0
     ){
     FlowDir=1;
     updatemenu=0;
   }

  if(plotstate==DYNAMIC_PLOTS){
    if(timedrag==0)itime += skip*FlowDir;
    checktimebound();
    IDLE();

    return;
  }
  if(current_mesh->visx != 0 && current_mesh->slicedir==1){
    current_mesh->plotx += skip*FlowDir;
    updateplotslice(1);
  }
  if(current_mesh->visy != 0 && current_mesh->slicedir==2){
    current_mesh->ploty += skip*FlowDir;
    updateplotslice(2);
  }
  if(current_mesh->visz != 0 && current_mesh->slicedir==3){
    current_mesh->plotz += skip*FlowDir;
    updateplotslice(3);
  }
  if(ReadPlot3dFile==1&&visiso !=0 && current_mesh->slicedir==4){
    plotiso[plotn-1] += FlowDir; 
    updatesurface(); 
  }
}
#ifdef pp_GPU_CULL_STATE
/* ------------------ print_gpu_cull_state ------------------------ */

void print_gpu_cull_state(void){
  char gpu_label[128];
#ifdef pp_CULL
  char cull_label[128];

  if(cullactive==1){
    if(cullsmoke==1&&usegpu==1){
      strcpy(cull_label,"Smoke culling in use.");
    }
    else if(cullsmoke==1&&usegpu==0){
      strcpy(cull_label,"Smoke culling not in use (available if GPU activates).");
    }
    else{
      strcpy(cull_label,"Smoke culling not in use.");
    }
  }
  else{
    strcpy(cull_label,"Smoke culling not available.");
  }
#endif
#ifdef pp_GPU
  if(gpuactive==1){
    if(usegpu==1){
      strcpy(gpu_label,"GPU in use.");
    }
    else{
      strcpy(gpu_label,"GPU not in use.");
    }
  }
  else{
    strcpy(gpu_label,"GPU not available.");
  }
  printf("%s ",gpu_label);
#endif
#ifdef pp_CULL
  printf("%s",cull_label);
#endif
  printf("\n");
}
#endif

/* ------------------ handle_eyeview ------------------------ */

void handle_eyeview(int flag){
  float *angle_zx;

  if(eyeview==eyeview_old)return;
  camera_current->eyeview=eyeview;
  angle_zx = camera_current->angle_zx;
  updatemenu=1;
  switch (eyeview){
  case WORLD_CENTERED:
      if(trainer_mode==0)printf("world centered\n");
      if(showtrainer==0&&flag==0&&eyeview_old==EYE_CENTERED){
        ResetView(RESTORE_EXTERIOR_VIEW);
      }
      break;
  case EYE_CENTERED:
       angle_zx[1]=0.0;
       if(showtrainer==0&&flag==0&&eyeview_old!=EYE_CENTERED){
         ResetView(RESTORE_EXTERIOR_VIEW);
       }
       update_glui_speed();
      if(trainer_mode==0)printf("eye centered\n");
      break;
  case WORLD_CENTERED_LEVEL:
    angle_zx[1]=0.0;
    if(trainer_mode==0)printf("world centered, level rotations\n");
    if(showtrainer==0&&flag==0&&eyeview_old==EYE_CENTERED){
      ResetView(RESTORE_EXTERIOR_VIEW);
    }
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  showhide_translate(eyeview);
  eyeview_old = eyeview;
  return;
}

/* ------------------ update_clipplanes ------------------------ */

void update_clipplanes(void){
  if(trainer_mode==0){
    if(xyz_clipplane!=xyz_clipplane_last){
      if(xyz_clipplane==0)printf("clipping off\n");
      if(xyz_clipplane==1)printf("clipping blockages + data\n");
      if(xyz_clipplane==2)printf("clipping blockages\n");
      xyz_clipplane_last=xyz_clipplane;
    }
  }
  if(xyz_clipplane==0){
    glDisable(GL_CLIP_PLANE0);
    glDisable(GL_CLIP_PLANE1);
    glDisable(GL_CLIP_PLANE2);
    glDisable(GL_CLIP_PLANE3);
    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);

  }
}

/* ------------------ handleiso ------------------------ */

void handleiso(void){
    if(ReadPlot3dFile==1){
      updateshowstep(1-visiso,ISO);
      if(visiso==1){
        updatesurface();
        plotstate=STATIC_PLOTS;
      }
    }
    updateglui();
    return;
}

/* ------------------ specialkeyoard ------------------------ */

void specialkeyboard_up(int key, int x, int y){
  resetclock=1;
}

/* ------------------ specialkeyoard ------------------------ */

void specialkeyboard(int key, int x, int y){

#define EYE_MODE 0
#define P3_MODE 1
  int keymode=EYE_MODE;

  GLUTPOSTREDISPLAY

  switch (cursorPlot3D){
    case 0:
      if(eyeview==EYE_CENTERED){
        keymode=EYE_MODE;
      }
      else{
        keymode=P3_MODE;
      }
      break;
    case 1:
      if(visGrid==1||plotstate==STATIC_PLOTS
  ||ReadVolSlice==1
        ){
        keymode=P3_MODE;
      }
      else{
        keymode=EYE_MODE;
      }
      break;
    default:
      ASSERT(FFALSE);
      break;
  }

  switch (keymode){
    case P3_MODE:
      handle_plot3d_keys(key);
      stept=0;
      break;
    case EYE_MODE:
      handle_move_keys(key);
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
} 

/* ------------------ handle_plot3d_keys ------------------------ */

void handle_plot3d_keys(int  key){
  int i;

  switch (key){
  case GLUT_KEY_LEFT:
  case GLUT_KEY_RIGHT:
    if(current_mesh->visx2==0){
      for(i=0;i<nmeshes;i++){
        mesh *meshi;
        meshi = meshinfo+i;
        if(meshi->visx2==1){
          update_current_mesh(meshi);
          break;
        }
      }
    }
    updateshowstep(1,DIRX);
    if(key==GLUT_KEY_LEFT){
      current_mesh->plotx -= skip;
    }
    if(key==GLUT_KEY_RIGHT){
      current_mesh->plotx += skip;
    }
    updateplotslice(1);
    break;
  case GLUT_KEY_DOWN:
  case GLUT_KEY_UP:
    updateshowstep(1,DIRY);
    if(key==GLUT_KEY_UP){
      current_mesh->ploty += skip;
    }
    if(key==GLUT_KEY_DOWN){
      current_mesh->ploty -= skip;
    }
    updateplotslice(2);
    break;
  case GLUT_KEY_PAGE_UP:
  case GLUT_KEY_PAGE_DOWN:
    updateshowstep(1,DIRZ);
    if(key==GLUT_KEY_PAGE_UP){
      current_mesh->plotz += skip;
    }
    if(key==GLUT_KEY_PAGE_DOWN){
      current_mesh->plotz -= skip;
    }
    updateplotslice(3);
    break;
  case GLUT_KEY_HOME:
    if(current_mesh->slicedir==1){
      current_mesh->plotx=0;
      updateplotslice(1);
    }
    if(current_mesh->slicedir==2){
      current_mesh->ploty=0;
      updateplotslice(2);
    }
    if(current_mesh->slicedir==3){
      current_mesh->plotz=0;
      updateplotslice(3);
    }
    break;
  case GLUT_KEY_END:
    if(current_mesh->slicedir==1){
      current_mesh->plotx=current_mesh->ibar;
      updateplotslice(1);
    }
    if(current_mesh->slicedir==2){
      current_mesh->ploty=current_mesh->jbar;
      updateplotslice(2);
    }
    if(current_mesh->slicedir==3){
      current_mesh->plotz=current_mesh->kbar;
      updateplotslice(3);
    }
    break;
  default:
    printf("warning key stroke=%i not handled\n",key);
    break;
  }
  plotstate=getplotstate(STATIC_PLOTS);
  //updatemenu=1;

} 

/* ------------------ training_move ------------------------ */
                                                 
void training_move(int  mode){
  float dx, dy;
  float *cos_direction_angle, *sin_direction_angle;
//  float *elevation_angle;
//  float *cos_elevation_angle, *sin_elevation_angle;

  float *eye_xyz;
  float *direction_angle;

  float INC_ANGLE=11.25;
  float INC_XY, INC_Z;

  INC_XY=meshinfo->cellsize/xyzmaxdiff;
  INC_Z=0.1/xyzmaxdiff;

  eye_xyz = camera_current->eye;
  direction_angle=&camera_current->direction_angle;

  cos_direction_angle=&camera_current->cos_direction_angle;
  sin_direction_angle=&camera_current->sin_direction_angle;

//  elevation_angle=&camera_current->elevation_angle;
//  cos_elevation_angle=&camera_current->cos_elevation_angle;
//  sin_elevation_angle=&camera_current->sin_elevation_angle;

#define ANGLE_LEFT -1
#define ANGLE_RIGHT 1
#define GO_FORWARD 0
#define GO_BACKWARD 2
#define MOVE_UP 3
#define MOVE_DOWN 4


  if(mode==ANGLE_LEFT||mode==ANGLE_RIGHT){
    switch (mode){
    case ANGLE_LEFT:
      *direction_angle-=INC_ANGLE;
      break;
    case ANGLE_RIGHT:
      *direction_angle+=INC_ANGLE;
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
    *cos_direction_angle = cos(PI*(*direction_angle)/180.0);
    *sin_direction_angle = sin(PI*(*direction_angle)/180.0);
  }

  dx = INC_XY*(*sin_direction_angle);
  dy = INC_XY*(*cos_direction_angle);

  if(mode==GO_FORWARD){
    getnewpos(eye_xyz,dx,dy,0.0,1.0);
  }
  else if(mode==GO_BACKWARD){
    getnewpos(eye_xyz,-dx,-dy,0.0,1.0);
  }
  else if(mode==MOVE_UP){
    eye_xyz[2] += INC_Z;
  }
  else if(mode==MOVE_DOWN){
    eye_xyz[2] -= INC_Z;
  }

  eye_xyz0[0]=eye_xyz[0];
  eye_xyz0[1]=eye_xyz[1];
  eye_xyz0[2]=eye_xyz[2];
  update_translate();
  GLUTPOSTREDISPLAY

} 

/* ------------------ handle_move_keys ------------------------ */
                                                 
void handle_move_keys(int  key){
  int state;
  float dx, dy;
  float *cos_direction_angle, *sin_direction_angle;
  float *elevation_angle, *cos_elevation_angle, *sin_elevation_angle;

  const float INC_ANGLE0=0.1;

  float INC_XY, INC_Z, INC_ANGLE;

  float *eye_xyz;
  float *direction_angle;
#define LOOKANGLE_CHANGE 11.25


  eye_xyz = camera_current->eye;
  direction_angle=&camera_current->direction_angle;

  cos_direction_angle=&camera_current->cos_direction_angle;
  sin_direction_angle=&camera_current->sin_direction_angle;

  elevation_angle=&camera_current->elevation_angle;
  cos_elevation_angle=&camera_current->cos_elevation_angle;
  sin_elevation_angle=&camera_current->sin_elevation_angle;


  glui_move_mode=-1;

  INC_XY=meshinfo->cellsize/xyzmaxdiff;
  INC_Z=INC_XY;
  INC_ANGLE = 5*INC_ANGLE0;

  state=glutGetModifiers();
  switch (state){
  case GLUT_ACTIVE_CTRL:
    key_state = KEY_CTRL;
    pass_through=1;
    update_blockpath();
    break;
  case GLUT_ACTIVE_ALT:
    key_state = KEY_ALT;
    break;
  case GLUT_ACTIVE_SHIFT:
    key_state = KEY_SHIFT;
    break;
  default:
    key_state = KEY_NONE;
    if(pass_through==1){
      pass_through=0;
      update_blockpath();
    }
    break;
  }
  switch (key){
  case GLUT_KEY_RIGHT:
    switch (key_state){
    case KEY_ALT:
      dx = INC_XY*(*cos_direction_angle);
      dy = INC_XY*(*sin_direction_angle);
      getnewpos(eye_xyz,dx,-dy,0.0,1.0);
//      eye_xyz[0] += dx;
//      eye_xyz[1] -= dy;
      break;
    case KEY_SHIFT:
    case KEY_CTRL:
    default:
      if(key_state==KEY_SHIFT){
        *direction_angle += 4.0*INC_ANGLE;
      }
      else{
        *direction_angle += INC_ANGLE;
      }
      *cos_direction_angle = cos(PI*(*direction_angle)/180.0);
      *sin_direction_angle = sin(PI*(*direction_angle)/180.0);
      break;
    }
    break;
  case 256+'d':
    dx = INC_XY*(*cos_direction_angle);
    dy = INC_XY*(*sin_direction_angle);
    {
      float local_speed_factor=1.0;

      switch (key_state){
      case KEY_SHIFT:
        local_speed_factor=4.0;
      default:
      getnewpos(eye_xyz,dx,-dy,0.0,local_speed_factor);
      break;
      }
    }
    break;
  case GLUT_KEY_LEFT:
    switch (key_state){
    case KEY_ALT:
      dx = INC_XY*(*cos_direction_angle);
      dy = INC_XY*(*sin_direction_angle);
      getnewpos(eye_xyz,-dx,dy,0.0,1.0);
      break;
    case KEY_SHIFT:
    case KEY_CTRL:
    default:
      if(key_state==KEY_SHIFT){
        *direction_angle -= 4.0*INC_ANGLE;
      }
      else{
        *direction_angle -= INC_ANGLE;
      }
      *cos_direction_angle = cos(PI*(*direction_angle)/180.0);
      *sin_direction_angle = sin(PI*(*direction_angle)/180.0);
      break;
    }
    break;
  case 256+'a':
    dx = INC_XY*(*cos_direction_angle);
    dy = INC_XY*(*sin_direction_angle);
    if(key_state==KEY_SHIFT){
      getnewpos(eye_xyz,-dx,dy,0.0,4.0);
    }
    else{
      getnewpos(eye_xyz,-dx,dy,0.0,1.0);
    }
    break;
  case GLUT_KEY_DOWN:
    {
      float local_speed_factor=1.0;

      switch (key_state){
      case KEY_ALT:
        eye_xyz[2] -= INC_Z;
        break;
      case KEY_SHIFT:
        speed_factor=4.0;
      case KEY_CTRL:
      default:
        dx = INC_XY*(*sin_direction_angle);
        dy = INC_XY*(*cos_direction_angle);
        getnewpos(eye_xyz,-dx,-dy,0.0,local_speed_factor);
        break;
      }
    }
    break;
  case GLUT_KEY_UP:
    {
      float local_speed_factor=1.0;

      switch (key_state){
      case KEY_ALT:
        eye_xyz[2] += INC_Z;
        break;
      case KEY_SHIFT:
        speed_factor=4.0;
      case KEY_CTRL:
      default:
        dx = INC_XY*(*sin_direction_angle);
        dy = INC_XY*(*cos_direction_angle);
        getnewpos(eye_xyz,dx,dy,0.0,local_speed_factor);
        break;
      }
    }
    break;
  case GLUT_KEY_PAGE_UP:
    *elevation_angle += LOOKANGLE_CHANGE;
    *cos_elevation_angle=cos(*elevation_angle*PI/180.0);
    *sin_elevation_angle=sin(*elevation_angle*PI/180.0);
    break;
  case GLUT_KEY_HOME:
    *elevation_angle=0.0;
    *cos_elevation_angle=1.0;
    *sin_elevation_angle=0.0;
    break;
  case GLUT_KEY_INSERT:
  case GLUT_KEY_PAGE_DOWN:
    *elevation_angle-=LOOKANGLE_CHANGE;
    *cos_elevation_angle=cos(*elevation_angle*PI/180.0);
    *sin_elevation_angle=sin(*elevation_angle*PI/180.0);
    break;
  case GLUT_KEY_END:
    ResetView(RESTORE_EXTERIOR_VIEW);
    break;
  case GLUT_KEY_F1:
    setspeed(speed_now/1.5);
    break;
  case GLUT_KEY_F2:
    setspeed(speed_walk);
    update_glui_speed();
    break;
  case GLUT_KEY_F3:
    setspeed(speed_now*1.5);
    break;
  case GLUT_KEY_F4:
    camera_current->view_angle-=LOOKANGLE_CHANGE;
    if(camera_current->view_angle<0.0)camera_current->view_angle+=360.0;
    camera_current->cos_view_angle=cos(camera_current->view_angle*PI/180.0);
    camera_current->sin_view_angle=sin(camera_current->view_angle*PI/180.0);
    break;
  case GLUT_KEY_F5:
    camera_current->view_angle=0.0;
    camera_current->cos_view_angle=1.0;
    camera_current->sin_view_angle=0.0;
    break;
  case GLUT_KEY_F6:
    camera_current->view_angle+=LOOKANGLE_CHANGE;
    if(camera_current->view_angle>360.0)camera_current->view_angle-=360.0;
    camera_current->cos_view_angle=cos(camera_current->view_angle*PI/180.0);
    camera_current->sin_view_angle=sin(camera_current->view_angle*PI/180.0);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  if(eyeview==EYE_CENTERED){
    eye_xyz0[0]=eye_xyz[0];
    eye_xyz0[1]=eye_xyz[1];
    eye_xyz0[2]=eye_xyz[2];
    update_translate();
  }
} 

/* ------------------ Display ------------------------ */

void Display(void){
  float *angle_zx;

  if(runscript==1&&default_script!=NULL){
    ScriptMenu(default_script->id);
    runscript=2;
  }
  script_render_flag=0;
  if(nscriptinfo>0&&current_script_command!=NULL){
    if(RenderGif==0){  // don't advance command if Smokeview is executing a RENDERALL command
      current_script_command++;
      script_render_flag=run_script();
      if(runscript==2&&noexit==0&&current_script_command==NULL){
        exit(0);
      }
      if(current_script_command==NULL){
        glui_script_enable();
      }
    }
    GLUTPOSTREDISPLAY
  }

  if(update_fire_line==1){
    WUI_CB(TERRAIN_FIRE_LINE_UPDATE);
    update_fire_line=0;
   }

  if(updatezoommenu==1||first_display>0){
     if(first_display>0)first_display--;
     updatezoommenu=0;
    ZoomMenu(zoomindex);
  }
  if(update_makeiblank_smoke3d==1){
    makeiblank_smoke3d();
  }
#ifdef pp_CULL
    if(update_initcull==1)initcull(cullsmoke);
#endif
   if(update_streaks==1&&ReadPartFile==1){
     void ParticleStreakShowMenu(int var);

     ParticleStreakShowMenu(streak_index);
     update_streaks=0;
   }
   if(update_screensize==1){
     update_screensize=0;
     update_windowsizelist();
     ResizeWindow(screenWidth,screenHeight);
   }
   if(update_colorbar_select_index==1&&colorbar_select_index>=0&&colorbar_select_index<=255){
     update_colorbar_select_index=0;
     updatecolors(colorbar_select_index);
   }
  if(updatemenu==1){
    if(menustatus==GLUT_MENU_NOT_IN_USE){
      glutDetachMenu(GLUT_RIGHT_BUTTON);
      InitMenus(LOAD);
      glutAttachMenu(GLUT_RIGHT_BUTTON);
    }
    else{
#ifdef _DEBUG
      printf("menus in use, will not be updated\n");
#endif
      /* 
      menus are being used used so keep re-displaying scene until
      user does something to cause menus to not be used
      */
    }
  }

  if(showtime==0&&ntotal_smooth_blockages>0){
    int i;

    for(i=0;i<nmeshes;i++){
      smoothblockage *sb;
      mesh *meshi;

      meshi = meshinfo+i;
      meshi->nsmoothblockagecolors=0;
      meshi->smoothblockagecolors=NULL;
      meshi->blockagesurfaces=NULL;

      if(meshi->smoothblockages_list!=NULL){
        sb=meshi->smoothblockages_list;
        if(sb!=NULL){
          meshi->nsmoothblockagecolors=sb->nsmoothblockagecolors;
          meshi->smoothblockagecolors=sb->smoothblockagecolors;
          meshi->blockagesurfaces=sb->smoothblockagesurfaces;
        }
      }
    }
  }
  glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
  if(showstereo==1&&videoSTEREO==1){  // temporal stereo (shuttered glasses)
    glDrawBuffer(GL_BACK_LEFT);
    if(showstereo_frame==0||showstereo_frame==2){
      ShowScene(RENDER,VIEW_LEFT,0,
        0,0,screenWidth,screenHeight);
    }
    glDrawBuffer(GL_BACK_RIGHT);
    if(showstereo_frame==1||showstereo_frame==2){
      ShowScene(RENDER,VIEW_RIGHT,0,
        0,0,screenWidth,screenHeight);
    }
    if(buffertype==DOUBLE_BUFFER&&benchmark_flag==0)glutSwapBuffers();
  }
  else if(showstereo==2){             // left/right stereo
    glDrawBuffer(GL_BACK);
    ClearBuffers(RENDER);
    if(showstereo_frame==0||showstereo_frame==2){
      ShowScene(RENDER,VIEW_LEFT,0,
        0,0,screenWidth/2,screenHeight);
    }
    if(showstereo_frame==1||showstereo_frame==2){
      ShowScene(RENDER,VIEW_RIGHT,0,
        screenWidth/2,0,screenWidth/2,screenHeight);
    }
    if(buffertype==DOUBLE_BUFFER&&benchmark_flag==0)glutSwapBuffers();
  }
  else if(showstereo==3){             // red/blue stereo
    glDrawBuffer(GL_BACK);
    glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
    glClearColor(1.0, 0.0, 0.0, 1.0); 
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if(showstereo_frame==0||showstereo_frame==2){
      glColorMask(GL_TRUE,GL_FALSE,GL_FALSE, GL_TRUE);
      ShowScene(RENDER,VIEW_LEFT,0,
        0,0,screenWidth,screenHeight);
      glFlush();
    }

    if(showstereo_frame==1||showstereo_frame==2){
      glDrawBuffer(GL_BACK);
      glColorMask(GL_FALSE,GL_FALSE,GL_TRUE,GL_TRUE);
      glClearColor(0.0, 0.0, 1.0, 1.0); 
      glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

      ShowScene(RENDER,VIEW_RIGHT,0,
        0,0,screenWidth,screenHeight);
      glFlush();
    }
    if(buffertype==DOUBLE_BUFFER&&benchmark_flag==0)glutSwapBuffers();
  }
  else if(showstereo==4){             // red/cyan stereo
    glDrawBuffer(GL_BACK);
    glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
    glClearColor(1.0, 0.0, 0.0, 1.0); 
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if(showstereo_frame==0||showstereo_frame==2){
      glColorMask(GL_TRUE,GL_FALSE,GL_FALSE, GL_TRUE);
      ShowScene(RENDER,VIEW_LEFT,0,
        0,0,screenWidth,screenHeight);
      glFlush();
    }

    if(showstereo_frame==1||showstereo_frame==2){
      glDrawBuffer(GL_BACK);
      glColorMask(GL_FALSE,GL_TRUE,GL_TRUE,GL_TRUE);
      glClearColor(0.0, 1.0, 1.0, 0.0); 
      glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

      ShowScene(RENDER,VIEW_RIGHT,0,
        0,0,screenWidth,screenHeight);
      glFlush();
    }
    if(buffertype==DOUBLE_BUFFER&&benchmark_flag==0)glutSwapBuffers();
  }
  else if(showstereo==5){             // custom red/blue stereo
    glDrawBuffer(GL_BACK);
    glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
    glClearColor(1.0, 1.0, 1.0, 1.0); 
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if(showstereo_frame==0||showstereo_frame==2){
      glColorMask(GL_TRUE,GL_FALSE,GL_FALSE, GL_TRUE);
      ShowScene(RENDER,VIEW_LEFT,0,0,0,screenWidth,screenHeight);
      glFlush();
    }
    if(showstereo_frame==1||showstereo_frame==2){
      glDrawBuffer(GL_BACK);
      glColorMask(GL_FALSE,GL_TRUE,GL_TRUE,GL_TRUE);
      glClearColor(0.0, 1.0, 1.0, 1.0); 
      glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

      ShowScene(RENDER,VIEW_RIGHT,0,0,0,screenWidth,screenHeight);

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();

      glEnable(GL_BLEND);
      glDisable(GL_LIGHTING);
      glDisable(GL_COLOR_MATERIAL);
      glDisable(GL_DITHER);
//      glShadeModel(GL_FLAT);

      glBlendFunc(GL_DST_COLOR,GL_ZERO);
      glBegin(GL_QUADS);
      glColor4f(0.0,right_green,right_blue,1.0);
      glVertex3f(-1.0,-1.0,0.1);
      glVertex3f(1.0,-1.0,0.1);
      glVertex3f(1.0,1.0,0.1);
      glVertex3f(-1.0,1.0,0.1);
      glEnd();

      glFlush();
    }
    if(buffertype==DOUBLE_BUFFER&&benchmark_flag==0)glutSwapBuffers();
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_DITHER);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
  }
  else{
    if(benchmark_flag==1){
      glDrawBuffer(GL_FRONT);
      ShowScene(RENDER,VIEW_CENTER,0,0,0,screenWidth,screenHeight);
    }
    else{
      if(render_double==0){
        glDrawBuffer(GL_BACK);
        ShowScene(RENDER,VIEW_CENTER,0,0,0,screenWidth,screenHeight);
        if(buffertype==DOUBLE_BUFFER)glutSwapBuffers();
      }
#ifdef pp_RENDER
      else{
        GLubyte *screenbuffers[4];
        int renderdoublenow=0;

        if(RenderOnceNow==1){
          renderdoublenow=1;
        }
    
        if(plotstate==DYNAMIC_PLOTS && ntimes > 0){
          if(itime>=0&&itime<ntimes&&
            ((render_frame[itime] == 0&&showstereo==0)||(render_frame[itime]<2&&showstereo!=0))
            ){
            render_frame[itime]++;
            renderdoublenow=1;
          }
        }

        if(renderdoublenow==1){
          glDrawBuffer(GL_BACK);

          ShowScene(RENDER,VIEW_CENTER,1,0,0,screenWidth,screenHeight);
          screenbuffers[0]=getscreenbuffer();
          if(buffertype==DOUBLE_BUFFER)glutSwapBuffers();

          ShowScene(RENDER,VIEW_CENTER,1,screenWidth,0,screenWidth,screenHeight);
          screenbuffers[1]=getscreenbuffer();
          if(buffertype==DOUBLE_BUFFER)glutSwapBuffers();

          ShowScene(RENDER,VIEW_CENTER,1,0,screenHeight,screenWidth,screenHeight);
          screenbuffers[2]=getscreenbuffer();
          if(buffertype==DOUBLE_BUFFER)glutSwapBuffers();

          ShowScene(RENDER,VIEW_CENTER,1,screenWidth,screenHeight,screenWidth,screenHeight);
          screenbuffers[3]=getscreenbuffer();
          if(buffertype==DOUBLE_BUFFER)glutSwapBuffers();

          mergescreenbuffers(screenbuffers);

          if(screenbuffers[0]!=NULL)free(screenbuffers[0]);
          if(screenbuffers[1]!=NULL)free(screenbuffers[1]);
          if(screenbuffers[2]!=NULL)free(screenbuffers[2]);
          if(screenbuffers[3]!=NULL)free(screenbuffers[3]);
        }
        if(renderdoublenow==0||RenderOnceNow==1){
          ASSERT(RenderSkip>0);
          RenderState(0);
          RenderSkip=1;
        }
      }
#endif
    }
  }
  if(touring == 1 ){
    if(RenderGif != 0){
      if(ntimes>0)angle += 2.0*PI/((float)ntimes/(float)RenderSkip);
      if(ntimes==0)angle += 2.0*PI/((float)maxtourframes/(float)RenderSkip);
    }
//    if(RenderGif == 0)angle += dang/lastcount;
    if(RenderGif == 0)angle += dang;
    if(angle>PI){angle -= -2.0f*PI;}
    if(eyeview==WORLD_CENTERED||eyeview==WORLD_CENTERED_LEVEL){
        angle_zx = camera_current->angle_zx;
        angle_zx[0] = anglexy0 + angle*180./PI;
    }
    else{          
      camera_current->direction_angle = direction_angle0 + angle*180./PI;
      camera_current->cos_direction_angle = cos(PI*camera_current->direction_angle/180.0);
      camera_current->sin_direction_angle = sin(PI*camera_current->direction_angle/180.0);
    }
    GLUTPOSTREDISPLAY
  }

}

/* ------------------ Idle ------------------------ */

void Idle(void){
  int changetime=0;
  float thisinterval;
  int oldcpuframe;
  float totalcpu;
  int redisplay=0;
  int ibenchrate;
  char buffer[256];
  float elapsed_time;

  glutSetWindow(mainwindow_id);
  updateShow();
  thistime = glutGet(GLUT_ELAPSED_TIME);
  thisinterval = thistime - lasttime;
//  printf("lasttime=%i thistime=%i thisinterval=%f\n",lasttime,thistime,thisinterval);
  count++;

  /* increment frame counter if the proper amount of time has passed
     or if we are rendering images or stepping by hand */

  if(showtime==1&&((stept==1&&(float)thisinterval>frameinterval)||RenderGif!=0||timedrag==1)){       /* ready for a new frame */

    cputimes[cpuframe]=thistime/1000.;
    
    oldcpuframe=cpuframe-10;
    if(oldcpuframe<0)oldcpuframe+=20;
    totalcpu=cputimes[cpuframe]-cputimes[oldcpuframe];
    if(benchmark==0){
      if(totalcpu==0.0){
   		  framerate=0.0;
      }
      else{
	      framerate=10.0/totalcpu;
      }
    }
    if(benchmark==1||benchmark_flag==1){
      if(itime==0)bench_starttime=thistime/1000.0;
      if(itime==ntimes-1){
        bench_stoptime=thistime/1000.0;
        ibenchrate=10*((float)ntimes/(bench_stoptime-bench_starttime))+0.5;
        framerate=(float)ibenchrate/10.0;
        sprintf(buffer,"%f",framerate);
        trim(buffer);
        trimzeros(buffer);
        printf("   frame rate=%s\n",buffer);
        if(benchmark_flag==1)bench_out(framerate);
      }
    }
    cpuframe++;
    if(cpuframe>=20)cpuframe=0;
   
    lastcount=count;
    count=1;
    lasttime = thistime;
    if(ntimes>0){
      changetime=1;
      if(stept ==1 && plotstate == DYNAMIC_PLOTS && timedrag==0 && RenderGif==0){
        /*  skip frames here if displaying in real time and frame rate is too slow*/
        if(times!=NULL&&realtime_flag!=0&&FlowDir>0){
          elapsed_time = (float)thistime/1000.0 - reset_time;
          elapsed_time *= (float)realtime_flag;
          elapsed_time += times[reset_frame];
          if(ntimes>1&&
            times[ntimes-1]>times[0]&&
            (elapsed_time>times[ntimes-1]||elapsed_time<0.0)
            ){
            elapsed_time = gmod(elapsed_time,times[ntimes-1]-times[0])+times[0];
          }
          itime = interval_search(times,ntimes,elapsed_time,itime);
        }
        else{
          if(script_render_flag==0){
            itime+=FlowDir;
          }
          else{
            itime=script_itime;
          }
        }
      }
      if(stept==1&&timedrag==0&&RenderGif!=0){
        itime+=RenderSkip*FlowDir;
      }

// if toggling time display with H then show the frame that was visible

      if(stept==0){
        itime_save=-1;
      }
      else{
        if(itime_save>=0){
          itime=itime_save;
        }
      }
#ifdef pp_SHOOTER
      if(shooter_firstframe==1&&visShooter!=0&&shooter_active==1){
        itime=0;
      }
#endif
      checktimebound();
      UpdateTimeLabels();
    }
    redisplay=1;
  }
  if(showtime==1&&stept==0&&itimeold!=itime){
    changetime=1;
    checktimebound();
    UpdateTimeLabels();
  }

  update_framenumber(changetime);
  if(redisplay==1){
    GLUTPOSTREDISPLAY
  }
}

/* ------------------ update_framenumber ------------------------ */

void update_framenumber(int changetime){
  int i,ii;
//  int redisplay;
  particle *parti;
  slice *sd;

  if(force_redisplay==1||(itimeold!=itime&&changetime==1)){
    force_redisplay=0;
    itimeold=itime;
//    redisplay=1;
    if(showsmoke==1||showevac==1){
      for(i=0;i<npart_files;i++){
        parti = partinfo+i;
        if(parti->loaded==1){
          if(parti->ptimeslist==NULL)continue;
          parti->iframe=parti->ptimeslist[itime];
        }
      }
    }
    if(hrrinfo!=NULL&&hrrinfo->loaded==1&&hrrinfo->display==1&&hrrinfo->timeslist!=NULL){
      hrrinfo->itime=hrrinfo->timeslist[itime];
    }
    if(showslice==1||showvslice==1){
      for(ii=0;ii<nslice_loaded;ii++){
        i = slice_loaded_list[ii];
        sd = sliceinfo+i;
        if(sd->slicetimeslist==NULL)continue;
        sd->islice=sd->slicetimeslist[itime];
      }
    }
    {
      smoke3d *smoke3di;

      if(show3dsmoke==1){
        for(i=0;i<nsmoke3d_files;i++){
          smoke3di = smoke3dinfo + i;
          if(smoke3di->loaded==0||smoke3di->display==0)continue;
          smoke3di->iframe=smoke3di->timeslist[itime];
          if(smoke3di->iframe!=smoke3di->lastiframe){
            smoke3di->lastiframe=smoke3di->iframe;
            updatesmoke3d(smoke3di);
          }
        }
        if(nsmoke3d_files>0)mergesmoke3dcolors();
      }
    }
    if(showpatch==1){
      for(i=0;i<nmeshes;i++){
        patch *patchi;
        mesh *meshi;

        meshi = meshinfo+i;
        patchi=patchinfo + meshi->patchfilenum;
        if(meshi->patchtimes==NULL)continue;
        if(meshi->patchtimeslist==NULL)continue;
        meshi->ipatch=meshi->patchtimeslist[itime];
        if(patchi->compression_type==0){
          meshi->ipqqi = meshi->ipqq + meshi->ipatch*meshi->npatchsize;
        }
        else{
#ifdef USE_ZLIB
          uncompress_patchdataframe(meshi,meshi->ipatch);
#endif
        }
      }
    }
    if(showiso==1){
      iso *isoi;
      mesh *meshi;

      CheckMemory;
      for(i=0;i<niso_files;i++){
        isoi = isoinfo + i;
        if(isoi->loaded==0)continue;
        meshi = meshinfo + isoi->blocknumber;

        if(meshi->isotimes==NULL)continue;
        if(meshi->isotimeslist==NULL)continue;
        meshi->iiso=meshi->isotimeslist[itime];

        if(isoi->compression_type==1){
          isosurface *asurface;
          asurface = meshi->animatedsurfaces + meshi->iiso*meshi->nisolevels-1;
          CheckMemory;
#ifdef USE_ZLIB
          uncompress_isodataframe(asurface,meshi->nisolevels);
#endif
          CheckMemory;
        }
      }
    }
    if(ntotal_smooth_blockages>0){
      for(i=0;i<nmeshes;i++){
        smoothblockage *sb;
        mesh *meshi;

        meshi = meshinfo+i;
        if(meshi->showsmoothtimelist!=NULL){
          sb=meshi->showsmoothtimelist[itime];
          if(sb==NULL)continue;
          meshi->nsmoothblockagecolors=sb->nsmoothblockagecolors;
          meshi->smoothblockagecolors=sb->smoothblockagecolors;
          meshi->blockagesurfaces=sb->smoothblockagesurfaces;
        }
      }
    }
    if(showzone==1){
      izone=zonetlist[itime];
    }
  }

}

/* ------------------ interval_search ------------------------ */

int interval_search(float *list, int nlist, float key, int guess){
  /* 
     find val such that list[val]<=key<list[val+1] 
     start with val=guess
  */

  int low, mid, high;

  if(nlist<=2||key<list[0])return 0;
  if(key>=list[nlist-2])return nlist-2;
  if(guess<0)guess=0;
  if(guess>nlist-2)guess=nlist-2;
  if(list[guess]<=key&&key<list[guess+1])return guess;

  low = 0;
  high = nlist - 1;
  while(high-low>1){
    mid = (low+high)/2;
    if(list[mid]>key){
      high=mid;
    }
    else{
      low=mid;
    }
  }
  if(list[high]==key)return high;
  return low;
}

/* ------------------ Reshape ------------------------ */
void update_camera_ypos(camera *camera_data);

void Reshape(int width, int height){
  updatemenu=1;
  ratio = (float)width/(float)height;
  aspect = ratio;
  if(ratio<1.0){ratio=1.0/ratio;}
  screenWidth = width;            
  screenHeight = height;
  screenWidth2 = screenWidth - dwinWW;   
  screenHeight2 = screenHeight - dwinH;
  windowresized=1;
  update_camera_ypos(camera_external);
  update_windowsizelist();
#ifdef pp_GPU
  createDepthTexture();
#endif
}

/* ------------------ togglegridstate ------------------------ */

void togglegridstate(int visg){
  int i;
  mesh *meshi;

  visGrid=visg;
  if(visGrid==1){
    if(current_mesh->plotx==-1){
      current_mesh->plotx=current_mesh->ibar/2; 
    }
    if(current_mesh->ploty==-1){
      current_mesh->ploty=current_mesh->jbar/2; 
    }
    if(current_mesh->plotz==-1){
      current_mesh->plotz=current_mesh->kbar/2;
    }
    for(i=0;i<nmeshes;i++){
      meshi=meshinfo + i;
      if(meshi->visx!=0||meshi->visy!=0||meshi->visz!=0)return;
    }
    updateshowstep(1-current_mesh->visy,DIRY);
  }
}

/* ------------------ cputime ------------------------ */

float cputime(void){
  return (float)clock()/(float)CLOCKS_PER_SEC;
}

/* ------------------ gmod ------------------------ */

float gmod(float x, float y){
  float returnval;

  if(y==0.0)return 0.0;
  returnval = x - (int)(x/y)*y;
  if(returnval<0.0)returnval+=y;
  return returnval;
}

/* ------------------ reset_gltime ------------------------ */

void reset_gltime(void){
  int inttime;

  if(showtime!=1)return;
  reset_frame=itime;
  inttime  = glutGet(GLUT_ELAPSED_TIME);
  reset_time = (float)inttime/1000.0;
  if(times!=NULL&&ntimes>0){
    start_frametime=times[0];
    stop_frametime=times[ntimes-1];
  }
}

/* ------------------ update_currentmesh ------------------------ */

void update_current_mesh(mesh *meshi){
  int i;
  iso *isoi;

  current_mesh=meshi;
  if(isoinfo!=NULL&&current_mesh->isofilenum==-1){
    loaded_isomesh=NULL;
    for(i=0;i<niso_files;i++){
      isoi = isoinfo + i;
      if(isoi->loaded==0)continue;
      loaded_isomesh = meshinfo+isoi->blocknumber;
      break;
    }
  }
  else{
    loaded_isomesh=NULL;
  }
  update_iso_showlevels();
  update_plot3dtitle();
}

