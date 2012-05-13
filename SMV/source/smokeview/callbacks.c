// $Date$ 
// $Revision$
// $Author$

// svn revision character string
char callbacks_revision[]="$Revision$";

#include "options.h"
#include <stdio.h>  
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "string_util.h"
#include "update.h"
#include "smokeviewvars.h"
#include "IOvolsmoke.h"

#define KEY_ALT 0
#define KEY_CTRL 1
#define KEY_SHIFT 3
#define KEY_NONE 2

#define TERRAIN_FIRE_LINE_UPDATE 39

#undef pp_GPU_CULL_STATE
#ifdef pp_GPU
#define pp_GPU_CULL_STATE
#endif
#ifdef pp_CULL
#ifndef pp_GPU_CULL_STATE
#define pp_GPU_CULL_STATE
#endif
#endif

/* ------------------ next_xindex ------------------------ */

void next_xindex(int inc,int flag){
  int i,j,first;

  first=1;
  if(flag==1)inc=-1;
  if(flag==-1)inc=1;
  for(j=0;j<nplotx_all;j++){
    if(first==1){
      first=0;
      if(flag==1)iplotx_all=nplotx_all-1;
      if(flag==-1)iplotx_all=0;
      if(flag==0)iplotx_all+=inc;
    }
    else{
      iplotx_all+=inc;
    }
    if(iplotx_all<0)iplotx_all=nplotx_all-1;
    if(iplotx_all>nplotx_all-1)iplotx_all=0;
    if(visGrid!=noGridnoProbe)return;
    if(plotstate==DYNAMIC_PLOTS){
      for(i=0;i<nsliceinfo;i++){
        slicedata *slicei;
        mesh *meshi;

        slicei = sliceinfo + i;
        if(slicei->loaded==0||slicei->display==0)continue;
        meshi = meshinfo + slicei->blocknumber;
        if(meshi->iplotx_all[iplotx_all]!=-1)return;
      }
      for(i=0;i<nvsliceinfo;i++){
        vslicedata *vslicei;
        mesh *meshi;

        vslicei = vsliceinfo + i;
        if(vslicei->loaded==0||vslicei->display==0)continue;
        meshi = meshinfo + vslicei->val->blocknumber;
        if(meshi->iploty_all[iploty_all]!=-1)return;
      }
    }
    else{
      for(i=0;i<nplot3dinfo;i++){
        plot3ddata *plot3di;
        mesh *meshi;

        plot3di = plot3dinfo + i;
        if(plot3di->loaded==0||plot3di->display==0)continue;
        meshi = meshinfo + plot3di->blocknumber;
        if(meshi->iplotx_all[iplotx_all]!=-1)return;
      }
    }
  }
}

/* ------------------ next_yindex ------------------------ */

void next_yindex(int inc,int flag){
  int i,j,first;

  first=1;
  if(flag==1)inc=-1;
  if(flag==-1)inc=1;
  for(j=0;j<nploty_all;j++){
    if(first==1){
      first=0;
      if(flag==1)iploty_all=nploty_all-1;
      if(flag==-1)iploty_all=0;
      if(flag==0)iploty_all+=inc;
    }
    else{
      iploty_all+=inc;
    }
    if(iploty_all<0)iploty_all=nploty_all-1;
    if(iploty_all>nploty_all-1)iploty_all=0;
    if(visGrid!=noGridnoProbe)return;
    if(plotstate==DYNAMIC_PLOTS){
      for(i=0;i<nsliceinfo;i++){
        slicedata *slicei;
        mesh *meshi;

        slicei = sliceinfo + i;
        if(slicei->loaded==0||slicei->display==0)continue;
        meshi = meshinfo + slicei->blocknumber;
        if(meshi->iploty_all[iploty_all]!=-1)return;
      }
      for(i=0;i<nvsliceinfo;i++){
        vslicedata *vslicei;
        mesh *meshi;

        vslicei = vsliceinfo + i;
        if(vslicei->loaded==0||vslicei->display==0)continue;
        meshi = meshinfo + vslicei->val->blocknumber;
        if(meshi->iploty_all[iploty_all]!=-1)return;
      }
    }
    else{
      for(i=0;i<nplot3dinfo;i++){
        plot3ddata *plot3di;
        mesh *meshi;

        plot3di = plot3dinfo + i;
        if(plot3di->loaded==0||plot3di->display==0)continue;
        meshi = meshinfo + plot3di->blocknumber;
        if(meshi->iploty_all[iploty_all]!=-1)return;
      }
    }
  }
}

/* ------------------ next_zindex ------------------------ */

void next_zindex(int inc,int flag){
  int i,j,first;

  first=1;
  if(flag==1)inc=-1;
  if(flag==-1)inc=1;
  for(j=0;j<nplotz_all;j++){
    if(first==1){
      first=0;
      if(flag==1)iplotz_all=nplotz_all-1;
      if(flag==-1)iplotz_all=0;
      if(flag==0)iplotz_all+=inc;
    }
    else{
      iplotz_all+=inc;
    }
    if(iplotz_all<0)iplotz_all=nplotz_all-1;
    if(iplotz_all>nplotz_all-1)iplotz_all=0;
    if(visGrid!=noGridnoProbe)return;
    if(plotstate==DYNAMIC_PLOTS){
      for(i=0;i<nsliceinfo;i++){
        slicedata *slicei;
        mesh *meshi;

        slicei = sliceinfo + i;
        if(slicei->loaded==0||slicei->display==0)continue;
        meshi = meshinfo + slicei->blocknumber;
        if(meshi->iplotz_all[iplotz_all]!=-1)return;
      }
      for(i=0;i<nvsliceinfo;i++){
        vslicedata *vslicei;
        mesh *meshi;

        vslicei = vsliceinfo + i;
        if(vslicei->loaded==0||vslicei->display==0)continue;
        meshi = meshinfo + vslicei->val->blocknumber;
        if(meshi->iploty_all[iploty_all]!=-1)return;
      }
    }
    else{
      for(i=0;i<nplot3dinfo;i++){
        plot3ddata *plot3di;
        mesh *meshi;

        plot3di = plot3dinfo + i;
        if(plot3di->loaded==0||plot3di->display==0)continue;
        meshi = meshinfo + plot3di->blocknumber;
        if(meshi->iplotz_all[iplotz_all]!=-1)return;
      }
    }
  }
}

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
    devicedata *devicei;
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

/* ------------------ checktimebound ------------------------ */

void checktimebound(void){
  int i,j;
  slicedata *sd;
  mesh *meshi;
  blockagedata *bc;
  partdata *parti;

  if(timedrag==0&&itimes>nglobal_times-1||timedrag==1&&itimes<0){
    izone=0;
    itimes=0;
    iframe=iframebeg;
    for(i=0;i<nsliceinfo;i++){
      sd=sliceinfo+i;
      sd->itime=0;
    }
    for(i=0;i<nmeshes;i++){
      meshi=meshinfo+i;
      meshi->patch_itime=0;
    }
    for(i=0;i<nmeshes;i++){
      meshi=meshinfo+i;
      if(meshi->iso_times==NULL)continue;
      meshi->iso_itime=0;
    }
  }
  if(timedrag==0&&itimes<0||timedrag==1&&itimes>nglobal_times-1){
    izone=nzone_times-1;
    itimes=nglobal_times-1;
    for(i=0;i<npartinfo;i++){
      parti=partinfo+i;
      parti->itime=parti->ntimes-1;
    }
    for(i=0;i<nsliceinfo;i++){
      sd=sliceinfo+i;
      sd->itime=sd->ntimes-1;
    }
    for(i=0;i<nmeshes;i++){
      meshi=meshinfo+i;
      meshi->patch_itime=meshi->npatch_times-1;
    }
    for(i=0;i<nmeshes;i++){
      meshi=meshinfo+i;
      if(meshi->iso_times==NULL)continue;
      meshi->iso_itime=meshi->niso_times-1;
    }
  }
  /* set blockage visibility */

  for(i=0;i<nmeshes;i++){
    meshi=meshinfo+i;
    for(j=0;j<meshi->nbptrs;j++){
      bc=meshi->blockageinfoptrs[j];
      if(bc->showtimelist==NULL)continue;
      bc->show=bc->showtimelist[itimes];
    }
  }
}

/* ------------------ setup_colorbar_drag ------------------------ */

int setup_colorbar_drag(int x, int y){
  int temp;
  int yy;
  float factor;
  int ifactor;
  int state;

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
    if(state==GLUT_ACTIVE_CTRL&&current_colorbar!=NULL&&current_colorbar->nsplits==1){
      colorsplitdrag=1;
    }
    else{
      colordrag=1;
      updatecolors(ifactor);
    }
    return 1;
  }
  return 0;
}

/* ------------------ setup_timebar_drag ------------------------ */

int setup_timebar_drag(int x, int y){
  if(screenHeight-y<50&&nglobal_times>0){
    float xleft;

    if(fontindex==LARGE_FONT){
      xleft=xtimeleft+0.11;
    }
    else{
      xleft=xtimeleft;
    }
    itimes=(int)((xtemp*x/((screenWidth-dwinWW))-xleft)*(nglobal_times-1)/(xtimeright-xleft));
    checktimebound();
    timedrag=1;
    stept=0;
    Idle_CB();
    return 1;
  }
  return 0;
}

/* ------------------ mouse ------------------------ */

void mouse_CB(int button, int state, int x, int y){
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
  move_gslice=0;
  glutPostRedisplay();

  if(state==GLUT_UP){
#ifdef pp_MOUSEDOWN
    mouse_down=0;
#endif
    show_gslice_normal_keyboard=0;
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

#ifdef pp_MOUSEDOWN
  mouse_down=1;
#endif
  
  // check for double click for translating/rotating 3D slice plane

  if(show_gslice_data==1||show_gslice_triangles==1||show_gslice_triangulation==1){
    this_mouse_time=glutGet(GLUT_ELAPSED_TIME)/1000.0;
    if(this_mouse_time-last_mouse_time<0.5){
      gslice_xyz0[0]=gslice_xyz[0];
      gslice_xyz0[1]=gslice_xyz[1];
      gslice_xyz0[2]=gslice_xyz[2];
      gslice_normal_azelev0[0]=gslice_normal_azelev[0];
      gslice_normal_azelev0[1]=gslice_normal_azelev[1];
      move_gslice=1;
      show_gslice_normal_keyboard=1;
    }
    last_mouse_time=this_mouse_time;
  }
  if(button==GLUT_LEFT_BUTTON||button==GLUT_MIDDLE_BUTTON||button==GLUT_RIGHT_BUTTON){
    glutSetCursor(GLUT_CURSOR_INFO);

    /* edit blockages */

    if(button==GLUT_LEFT_BUTTON){
      if(blockageSelect==1)mouse_edit_blockage(button,state,x,y);
      if(edittour==1&&blockageSelect==0)mouse_edit_tour(button,state,x,y);
      if(select_avatar==1)mouse_select_avatar(button,state,x,y);
      if(select_device==1)mouse_select_device(button,state,x,y);
    }
    glutPostRedisplay();
    if( showtime==1 || showplot3d==1){
      if(setup_colorbar_drag(x,y)==1)return;
    }
    if(visTimeLabels==1&&showtime==1){
      if(setup_timebar_drag(x,y)==1)return;
    }
    copy_camera(camera_last,camera_current);
    if(canrestorelastview==0){
      updatemenu=1;
      canrestorelastview=1;
      enable_reset_saved_view();
    }
    switch (button){
      case GLUT_MIDDLE_BUTTON:
        state=GLUT_ACTIVE_CTRL;
        break;
      case GLUT_RIGHT_BUTTON:
        state=GLUT_ACTIVE_ALT;
        break;
      default:
        state=glutGetModifiers();
        break;
    }
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
    mouse_down_xy0[0]=x; 
    mouse_down_xy0[1]=y;
  }
  glutPostRedisplay();
  if(blockageSelect==1){
    Display_CB();
  }
}

/* ------------------ drag_colorbar ------------------------ */

void drag_colorbar(int xm, int ym){
  int temp;
  int ifactor;
  float factor;
  int valmax=255;
  int valmin=0;

  temp = (int)(1.2*dwinH);
  if(xm>screenWidth-dwinWW){
    float yy;

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
}

/* ------------------ drag_colorsplit ------------------------ */

void drag_colorbarsplit(int xm, int ym){
  int temp;
  int ifactor;
  float factor;

  temp = (int)(1.2*dwinH);
  if(xm>screenWidth-dwinWW){
    int ii;
    float yy;

    yy = screenHeight - ym;
    factor=(yy-temp)/(screenHeight-temp);
    factor *= (nrgb+1.0)/(nrgb-0.5);
    if(screenHeight>screenWidth)factor *= (float)screenHeight/screenWidth;
    ifactor=(int)(255*factor);

    if(ifactor>250)ifactor=250;
    if(ifactor<5)ifactor=5;
    ii=current_colorbar->splits[0];
    current_colorbar->index_node[ii]=ifactor;
    current_colorbar->index_node[ii-1]=ifactor;
    remapcolorbar(current_colorbar);
    updatecolors(-1);
    update_colorbar_splits(current_colorbar);
  }
}

/* ------------------ drag_timebar ------------------------ */

void drag_timebar(int xm, int ym){
  float xxleft;

  xxleft = xtimeleft;
  if(fontindex==LARGE_FONT)xxleft=xtimeleft+0.11;
  if(screenHeight-ym<50&&nglobal_times>0&&visTimeLabels==1&&showtime==1){
    itimes=(int)((xtemp*xm/((screenWidth-dwinWW))-xxleft)*(nglobal_times-1)/(xtimeright-xxleft));
    checktimebound();
    timedrag=1;
  }
  Idle_CB();
}

/* ------------------ Move_Gen_Slice ------------------------ */

void Move_Gen_Slice(int xm, int ym){
  int dxm, dym;
  int screenWidth2, screenHeight2;

  screenWidth2 = screenWidth - dwinWW;
  screenHeight2 = screenHeight - dwinH;

  dxm = xm - start_xyz0[0];
  dym = ym - start_xyz0[1];
  switch (key_state){
    case KEY_NONE:
      {
        float daz, delev;

        daz = 360.0*dxm/(float)screenWidth2;
        delev = 360.0*dym/(float)screenHeight2;
        gslice_normal_azelev[0] += daz;
        gslice_normal_azelev[1] += delev;
        update_gslice_parms();
        start_xyz0[0]=xm;
        start_xyz0[1]=ym;
      }
      break;
    case KEY_CTRL:
      {
        float xx, yy;
        float dx, dy;

        xx = xm-mouse_down_xy0[0];
        xx = xx/(float)screenWidth2;
        yy = ym-mouse_down_xy0[1];
        yy = yy/(float)screenHeight2;
        dx = (xyzbox+gslice_xyz0[0])*xx;
        dy = -(xyzbox-gslice_xyz0[1])*yy;
        gslice_xyz[0] += dx;
        gslice_xyz[1] += dy;
        gslice_xyz0[0] = gslice_xyz[0];
        gslice_xyz0[1] = gslice_xyz[1];
        mouse_down_xy0[0]=xm;
        mouse_down_xy0[1]=ym;
        update_gslice_parms();
      }
      break;
    case KEY_ALT:
      {
        float yy;

        yy = ym-mouse_down_xy0[1];
        yy = yy/(float)screenHeight2;

        gslice_xyz[2] = gslice_xyz0[2] - DENORMALIZE_Z(4*(xyzbox-NORMALIZE_Z(gslice_xyz0[2]))*yy);
        update_gslice_parms();
      }
      break;
    case KEY_SHIFT:
    break;
    default:
      ASSERT(0);
      break;
  }
}

/* ------------------ Move_Scene ------------------------ */

void Move_Scene(int xm, int ym){
  float *eye_xyz, *angle_zx;
  int screenWidth2, screenHeight2;
  int dxm, dym;
  float direction_angle;
  float elevation_angle;
  float xx, yy;

  eye_xyz = camera_current->eye;
  angle_zx = camera_current->angle_zx;
  screenWidth2 = screenWidth - dwinWW;
  screenHeight2 = screenHeight - dwinH;

  dxm = xm - start_xyz0[0];
  dym = ym - start_xyz0[1];
  switch (key_state){
    case KEY_NONE:
      switch (eyeview){
        case WORLD_CENTERED:
        case WORLD_CENTERED_LEVEL:
          angle_zx[0] += dxm;
          if(eyeview==WORLD_CENTERED){
            angle_zx[1] += dym;
          }
          else{
            angle_zx[1]=0.0;
          }
          start_xyz0[0]=xm;
          start_xyz0[1]=ym;
          break;
        case EYE_CENTERED:
#define ANGLE_FACTOR 0.25
          camera_current->direction_angle += dxm*ANGLE_FACTOR;
          direction_angle=camera_current->direction_angle;
          camera_current->cos_direction_angle = cos(PI*direction_angle/180.0);
          camera_current->sin_direction_angle = sin(PI*direction_angle/180.0);
          start_xyz0[0]=xm;

          camera_current->elevation_angle -= dym*ANGLE_FACTOR;
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

        xx = xm-mouse_down_xy0[0];
        xx = xx/(float)screenWidth2;
        yy = ym-mouse_down_xy0[1];
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
          mouse_down_xy0[0]=xm;
          mouse_down_xy0[1]=ym;
        }
      }
      break;

    case KEY_ALT:
      xx = xm-mouse_down_xy0[0];
      xx = xx/(float)screenWidth2;
      yy = ym-mouse_down_xy0[1];
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

/* ------------------ throttle_gpu ------------------------ */

int throttle_gpu(void){
  float fps;

  thisMOTIONtime=glutGet(GLUT_ELAPSED_TIME)/1000.0;
  fps = MOTIONnframes/(thisMOTIONtime-lastMOTIONtime);
  if(fps>GPU_VOLframemax)return 1;
  MOTIONnframes++;
  if(thisMOTIONtime>lastMOTIONtime+0.25){
    printf("MOTION: %4.1f fps\n",fps);
    lastMOTIONtime=thisMOTIONtime;
    MOTIONnframes=0;
  }
  return 0;
}

/* ------------------ motion ------------------------ */

void motion_CB(int xm, int ym){

#ifdef pp_GPUTHROTTLE
  if(usegpu==1&&showvolrender==1&&show_volsmoke_moving==1){
    if(throttle_gpu()==1)return;
  }
#endif
  
  reset_glui_view(-1);

  glutPostRedisplay();

  if( colordrag==1&&(showtime==1 || showplot3d==1)){
    drag_colorbar(xm,ym);
    return;
  }
  if(
    colorsplitdrag==1&&
    (showtime==1 || showplot3d==1)&&current_colorbar!=NULL&&current_colorbar->nsplits==1){
    drag_colorbarsplit(xm,ym);
    return;
  }
  if(timedrag==1){
    drag_timebar(xm,ym);
    return;
  }
  if(move_gslice==1){
    Move_Gen_Slice(xm,ym);
    return;
  }

  Move_Scene(xm,ym);
}

/* ------------------ keyboard_up_CB ------------------------ */

void keyboard_up_CB(unsigned char key, int x, int y){
  resetclock=1;
}

/* ------------------ get_vecfactor ------------------------ */

float get_vecfactor(int *ivec){
  float vec;
  if(*ivec>NVECLENGTHS-1)*ivec=0;
  if(*ivec<0)*ivec=NVECLENGTHS-1;
  vec=veclengths[*ivec];
  return vec;
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

/* ------------------ keyboard ------------------------ */

void keyboard(unsigned char key, int flag){
  char key2;
  int skip2;
  mesh *gbsave,*gbi;
  int i;
  int keystate=0;

  if(flag==FROM_CALLBACK)keystate = glutGetModifiers();
  glutPostRedisplay();
  key2 = (char)key;

  switch (key2){
    case 'a':
      if((eyeview==EYE_CENTERED||(visVector==1&&ReadPlot3dFile==1)||showvslice==1||isZoneFireModel==1)){
        if(eyeview==EYE_CENTERED){
          handle_move_keys(256+key2);
        }
        else{
          iveclengths += FlowDir;
          vecfactor = get_vecfactor(&iveclengths);
          update_vector_widgets();

          if(isZoneFireModel==1){
            if(iveclengths==0){
              zone_ventfactor=1.0;
            }
            else{
              if(FlowDir>0){
                zone_ventfactor*=1.5;
              } 
              else{
                zone_ventfactor/=1.5;
              }
            }
          }
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
      break;
    case 'A':
      axislabels_smooth=1-axislabels_smooth;
      update_axislabels_smooth();
      break;
    case 'c':
      switch (keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(18); // clip dialog
        break;
      case GLUT_ACTIVE_CTRL:
      default:
        if(nrooms>0){
          zone_highlight_room++;
          if(zone_highlight_room>=nrooms)zone_highlight_room=0;
          printf("room %i\n",zone_highlight_room+1);
        }
        else{
          contour_type++;
          if(contour_type>2)contour_type=0;
          update_plot3d_display();
          updatecolors(-1);
        }
      }
      break;
    case 'C':
      if(nrooms>0){
        zone_highlight = 1 - zone_highlight;
        if(zone_highlight==1){
          printf("room %i\n",zone_highlight_room+1);
        }
      }
#ifdef pp_CULL
      else{
        if(nsmoke3dinfo>0&&cullactive==1&&gpuactive==1){
          cullsmoke=1-cullsmoke;
          update_smoke3dflags();
          initcull(cullsmoke);
          print_gpu_cull_state();
        }
        if(cullactive==0||gpuactive==0)cullsmoke=0;
      }
#endif
      break;
    case 'd':
    case 'D':
      switch (keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(22); // display dialog
        break;
      case GLUT_ACTIVE_CTRL:
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
      break;
    case 'e':
    case 'E':
      switch (keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(16); // edit geometry
        break;
      case GLUT_ACTIVE_CTRL:
      default:
        eyeview++;
        if(eyeview>2)eyeview=0;
        handle_eyeview(0);
      }
      break;
    case 'f':
      switch (keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(14); // file/bounds dialog
        break;
      case GLUT_ACTIVE_CTRL:
      default:
        break;
      }
      break;
    case 'F':
      hide_overlaps=1-hide_overlaps;
      updatehiddenfaces=1;
      UpdateHiddenFaces();
      update_showhidebuttons();
      glutPostRedisplay();
      break;
    case 'g':
      switch (keystate){
      case GLUT_ACTIVE_ALT:
      case GLUT_ACTIVE_CTRL:
      default:
        if(ntotal_blockages>0||isZoneFireModel==0){
          switch (visGrid){
            case noGridnoProbe:
              visGrid=GridnoProbe;
              break;
            case GridnoProbe:
              visGrid=GridProbe;
              break;
            case GridProbe:
              visGrid=noGridProbe;
              break;
            case noGridProbe:
              visGrid=noGridnoProbe;
              break;
            default:
              visGrid=noGridnoProbe;
              break;
          }
          if(visGrid==GridProbe||visGrid==noGridProbe)visgridloc=1;
        }
        break;
      }
      break;
#ifdef pp_GPU
    case 'G':
      if(gpuactive==1){
        usegpu=1-usegpu;
      }
      else{
        usegpu=0;
      }
      if(nsmoke3dinfo>0){
        update_smoke3dflags();
      }
      print_gpu_cull_state();
      return;    
#endif
    case 'h':
      if(titlesafe_offset==0){
        titlesafe_offset=titlesafe_offsetBASE;
      }
      else{
        titlesafe_offset=0;
      }
      break;
    case 'H':
      {
        int nslice_loaded_local=0, nvslice_loaded_local=0;

        for(i=0;i<nsliceinfo;i++){
          slicedata *sd;

          sd = sliceinfo + i;
          if(sd->loaded==1)nslice_loaded_local++;
        }
        for(i=0;i<nvsliceinfo;i++){
          vslicedata *vd;

          vd = vsliceinfo + i;
          if(vd->loaded==1)nvslice_loaded_local++;
        }
        stept=1;
        if(nvslice_loaded_local>0){
          if(show_all_slices==0){
            ShowVSliceMenu(SHOW_ALL);
            force_redisplay=1;
          }
          else{
            itime_save=itimes;
            ShowVSliceMenu(HIDE_ALL);
          }
        }
        if(nvslice_loaded_local==0&&nslice_loaded_local>0){
          if(show_all_slices==0){
            ShowHideSliceMenu(SHOW_ALL);
            force_redisplay=1;
          }
          else{
            itime_save=itimes;
            ShowHideSliceMenu(HIDE_ALL);
          }
        }
      }
      break;
    case 'i':
    case 'I':
      if(unload_qdata==0){
        handleiso();
        return;
      }
      break;
    case 'k':
    case 'K':
      visTimeLabels = 1 - visTimeLabels;
      if(visTimeLabels==0)printf("Time bar hidden\n");
      if(visTimeLabels==1)printf("Time bar visible\n");
      break;
#ifdef _DEBUG 
    case 'l':
      if(nsmoke3dinfo>0){
        smokecullflag=1-smokecullflag;
        if(smokecullflag==0){
          smokedrawtest=1-smokedrawtest;
        }
        printf("smokecullflag=%i\n smokedrawtest=%i\n",smokecullflag,smokedrawtest);
        update_smoke3dflags();
        return;
      }
      break;
#endif
    case 'L':
      UnloadSliceMenu(-2);
      break;
    case 'm':
    case 'M':
      switch (keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(15); // motion dialog
        break;
      case GLUT_ACTIVE_CTRL:
      default:
        if(nmeshes>1){
          highlight_mesh++;
          if(highlight_mesh>nmeshes-1)highlight_mesh=0;
          update_current_mesh(meshinfo+highlight_mesh);
        }
      }
      break;
#ifdef _DEBUG
    case 'n':
    case 'N':
      if(nsmoke3dinfo>0){
        adjustalphaflag++;
        if(adjustalphaflag>3)adjustalphaflag=0;
        printf("adjustalphaflag=%i\n",adjustalphaflag);
        update_smoke3dflags();
        return;
      }
      break;
#endif
    case 'o':
    case 'O':
      highlight_flag++;
      if(highlight_flag>2&&noutlineinfo>0)highlight_flag=0;
      if(highlight_flag>1&&noutlineinfo==0)highlight_flag=0;
      printf("outline mode=%i\n",highlight_flag);
      break;
    case 'p':
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
      break;
    case 'P':
      cursorPlot3D=1-cursorPlot3D;
      update_cursor_checkbox();
      break;
    case 'q':
    case 'Q':
      blocklocation++;
      if(blocklocation>BLOCKlocation_cad||
         blocklocation>BLOCKlocation_exact&&ncadgeom==0){
         blocklocation=BLOCKlocation_grid;
      }
      if(showedit_dialog==1){
        if(blocklocation==BLOCKlocation_exact){
          blockage_as_input=1;
        }
        else{
          blockage_as_input=0;
        }
        OBJECT_CB(BLOCKAGE_AS_INPUT2);
      }
      break;
    case 'r':
    case 'R':
      {
        int rflag=0;

        if(keystate==GLUT_ACTIVE_ALT){
          research_mode=1-research_mode;
          update_research_mode();
          return;
        }


        if(strncmp((const char *)&key2,"R",1)==0){
          render_double=2;
          rflag=1;
        }
        else{
          render_double=0;
          if(render_from_menu==0){
            renderW=0;
            renderH=0;
          }
        }
        if(scriptoutstream!=NULL){
          if(nglobal_times>0){
            float timeval;

            timeval=global_times[itimes];
            fprintf(scriptoutstream,"SETTIMEVAL\n");
            fprintf(scriptoutstream," %f\n",timeval);
            if(nvolrenderinfo>0&&load_at_rendertimes==1){
              for(i=0;i<nmeshes;i++){
                mesh *meshi;
                volrenderdata *vr;
                int j;
                int framenum;
                float timediffmin;

                meshi = meshinfo + i;
                vr = &meshi->volrenderinfo;
                if(vr->fire==NULL||vr->smoke==NULL)continue;
                if(vr->loaded==0||vr->display==0)continue;
                timediffmin = ABS(timeval-vr->times[0]);
                framenum=0;
                for(j=1;j<vr->ntimes;j++){
                  float timediff;

                  timediff = ABS(vr->times[j]-timeval);
                  if(timediff<timediffmin){
                    timediffmin=timediff;
                    framenum=j;
                  }
                }
                fprintf(scriptoutstream,"LOADVOLSMOKEFRAME\n");
                fprintf(scriptoutstream," %i %i\n",i,framenum);
              }
            }
          }
          else{
            int show_plot3dkeywords=0;

            for(i=0;i<nmeshes;i++){
              mesh *meshi;
              plot3ddata *plot3di;
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
              fprintf(scriptoutstream," %i %i %i %i %f\n",i+1,1, plotn,visx_all,xp[meshi->plotx]);
              fprintf(scriptoutstream,"SHOWPLOT3DDATA\n");
              fprintf(scriptoutstream," %i %i %i %i %f\n",i+1,2, plotn,visy_all,yp[meshi->ploty]);
              fprintf(scriptoutstream,"SHOWPLOT3DDATA\n");
              fprintf(scriptoutstream," %i %i %i %i %f\n",i+1,3, plotn,visz_all,zp[meshi->plotz]);
              fprintf(scriptoutstream,"SHOWPLOT3DDATA\n");
              fprintf(scriptoutstream," %i %i %i %i %i\n",i+1,4, plotn,visiso,plotiso[plotn-1]);
            }
            if(show_plot3dkeywords==1){
              fprintf(scriptoutstream,"PLOT3DPROPS\n");
              fprintf(scriptoutstream," %i %i %i %i\n",plotn,visVector,iveclengths,contour_type);
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
        render_from_menu=0;
      }
      break;
    case 's':
      switch (keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(20); // 3d smoke dialog
        break;
      case GLUT_ACTIVE_CTRL:
        snap_view_angles();
        break;
      default:
        if(eyeview==EYE_CENTERED){
          handle_move_keys(GLUT_KEY_DOWN);
        }
        else{
          vectorskip++;
          if(vectorskip>4)vectorskip=1;
        }
      }
      break;
    case 'S':
      showstereoOLD=showstereo;
      showstereo++;
      if(showstereo>5)showstereo=0;
      if(showstereo==1&&videoSTEREO!=1)showstereo=2;
      Update_Glui_Stereo();
      break;
    case 't':
      switch (keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(21); // tour dialog
        break;
      case GLUT_ACTIVE_CTRL:
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
      break;
    case 'T':
      usetexturebar=1-usetexturebar;
      printf("usetexturebar=%i\n",usetexturebar);
      break;
    case 'u':
    case 'U':
      switch (keystate){
        case GLUT_ACTIVE_ALT:
          skip_slice_in_embedded_mesh = 1 - skip_slice_in_embedded_mesh;
          break;
        default:
          ReloadMenu(0);
          break;
      }
      break;
    case 'v':
      switch (keystate){
        case GLUT_ACTIVE_ALT:
          projection_type = 1 - projection_type;
          TRANSLATE_CB(PROJECTION);
          break;
        default:
          visVector=1-visVector;
          if(vectorspresent==0)visVector=0;
          updateglui();
          break;
      }
      break;
    case 'V':
      if(nvolrenderinfo>0){
        usevolrender=1-usevolrender;
        update_smoke3dflags();
#ifdef pp_GPU
        print_gpu_cull_state();
#endif
        return;    
      }
      break;
    case 'w':
      switch (keystate){
        case GLUT_ACTIVE_ALT:
          DialogMenu(26); // WUI dialog
          break;
        case GLUT_ACTIVE_CTRL:
        default:
          if(eyeview==EYE_CENTERED){
            handle_move_keys(GLUT_KEY_UP);
          }
          else if(SHOW_gslice_data==1){
            show_gslice_data = 1 - show_gslice_data;
            update_gslice_parms();
          }
          break;
      }
      break;
    case 'W':
      xyz_clipplane++;
      if(xyz_clipplane>2)xyz_clipplane=0;
      update_clip_all();
      break;
    case 'x':
    case 'X':
      if(keystate==GLUT_ACTIVE_ALT){
        DialogMenu(-2); // close all dialogs
      }
      else{
        visx_all=1-visx_all;
      }
      break;
    case 'y':
    case 'Y':
      visy_all = 1-visy_all;
      break;
    case 'z':
    case 'Z':
      if(keystate==GLUT_ACTIVE_ALT){
        DialogMenu(24); // compress dialog
      }
      else{
        visz_all = 1 - visz_all;
      }
      break;
    case '0':
      if(plotstate==DYNAMIC_PLOTS){
        updatetimes();
        reset_time_flag=1;
        return;
      }
      break;
    case '@':
      cell_center_text = 1 - cell_center_text;
      break;
    case '!':
      snap_view_angles();
      break;
    case '#':
      writeini(LOCAL_INI);
      break;
    case '$':
      trainer_active=1-trainer_active;
      if(trainer_active==1){
        printf("Trainer mode active\n");
        trainer_mode=1;
        show_glui_trainer();
      }
      if(trainer_active==0){
        printf("Trainer mode inactive\n");
        trainer_mode=0;
        hide_glui_trainer();
      }
      break;
    case '&':
      antialiasflag=1-antialiasflag;
      printf("antialiasflag=%i\n",antialiasflag);
      break;
  }

  skip2=key2-'1'+1;
  if(skip2>0&&skip2<10)skip_global=skip2;

  /* if not a directional key then return */

  if(strncmp((const char *)&key2,"-",1)!=0&&strncmp((const char *)&key2," ",1)!=0&&
     strncmp((const char *)&key2,"+",1)!=0&&strncmp((const char *)&key2,"=",1)!=0&&
     strncmp((const char *)&key2,"<",1)!=0&&strncmp((const char *)&key2,">",1)!=0&&
     strncmp((const char *)&key2,",",1)!=0&&strncmp((const char *)&key2,".",1)!=0&&
     strncmp((const char *)&key2,"_",1)!=0&&(skip2<=0||skip2>=10))return;

  if(xyz_clipplane!=0&&(
    strncmp((const char *)&key2,"<",1)==0||strncmp((const char *)&key2,",",1)==0||
    strncmp((const char *)&key2,">",1)==0||strncmp((const char *)&key2,".",1)==0)){

    if(strncmp((const char *)&key2,"<",1)==0||strncmp((const char *)&key2,",",1)==0){ClipDir=-1;}
     else if(strncmp((const char *)&key2,">",1)==0||strncmp((const char *)&key2,".",1)==0){ClipDir=1;}

    if(stepclip_x==1  )clip_i += skip_global*ClipDir;
    if(stepclip_y==1  )clip_j += skip_global*ClipDir;
    if(stepclip_z==1  )clip_k += skip_global*ClipDir;
    if(stepclip_X==1  )clip_I += skip_global*ClipDir;
    if(stepclip_Y==1  )clip_J += skip_global*ClipDir;
    if(stepclip_Z==1  )clip_K += skip_global*ClipDir;

    updateclipbounds(clip_x,&clip_i,clip_X,&clip_I,current_mesh->ibar);
    updateclipbounds(clip_y,&clip_j,clip_Y,&clip_J,current_mesh->jbar);
    updateclipbounds(clip_z,&clip_k,clip_Z,&clip_K,current_mesh->kbar);
    return;
  }

  if(strncmp((const char *)&key2,"-",1)==0||strncmp((const char *)&key2,"_",1)==0){
    FlowDir=-1;
  }
  else if(strncmp((const char *)&key2," ",1)==0||
     strncmp((const char *)&key2,"=",1)==0||
     strncmp((const char *)&key2,"+",1)==0
     ){
     FlowDir=1;
  }

  if(plotstate==DYNAMIC_PLOTS){
    if(timedrag==0)itimes += skip_global*FlowDir;
    checktimebound();
    Idle_CB();

    return;
  }
  switch (iplot_state){
    case 1:
      next_xindex(skip_global*FlowDir,0);
      break;
    case 0:
    case 2:
      next_yindex(skip_global*FlowDir,0);
      break;
    case 3:
      next_zindex(skip_global*FlowDir,0);
      break;
    default:
      ASSERT(0);
      break;
  }
  if(ReadPlot3dFile==1&&visiso !=0 && current_mesh->slicedir==4){
    plotiso[plotn-1] += FlowDir; 
    updatesurface(); 
  }
  if(iplot_state!=0)updateplotslice(iplot_state);
}

/* ------------------ keyboard ------------------------ */

void keyboard_CB(unsigned char key, int x, int y){
  keyboard(key,FROM_CALLBACK);
  glutPostRedisplay();
  updatemenu=1;
}

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
      if(showtrainer_dialog==0&&flag==0&&eyeview_old==EYE_CENTERED){
        ResetView(RESTORE_EXTERIOR_VIEW);
      }
      break;
  case EYE_CENTERED:
       angle_zx[1]=0.0;
       if(showtrainer_dialog==0&&flag==0&&eyeview_old!=EYE_CENTERED){
         ResetView(RESTORE_EXTERIOR_VIEW);
       }
      if(trainer_mode==0)printf("eye centered\n");
      break;
  case WORLD_CENTERED_LEVEL:
    angle_zx[1]=0.0;
    if(trainer_mode==0)printf("world centered, level rotations\n");
    if(showtrainer_dialog==0&&flag==0&&eyeview_old==EYE_CENTERED){
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

void specialkeyboard_up_CB(int key, int x, int y){
  resetclock=1;
}

/* ------------------ specialkeyoard ------------------------ */

void specialkeyboard_CB(int key, int x, int y){

#define EYE_MODE 0
#define P3_MODE 1
  int keymode=EYE_MODE;

  glutPostRedisplay();

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
      if(visGrid!=noGridnoProbe||plotstate==STATIC_PLOTS||ReadVolSlice==1){
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
  switch (key){
  case GLUT_KEY_LEFT:
    visx_all=1;
    next_xindex(-1,0);
    iplot_state=1;
    break;
  case GLUT_KEY_RIGHT:
    visx_all=1;
    next_xindex(1,0);
    iplot_state=1;
    break;
  case GLUT_KEY_DOWN:
    visy_all=1;
    next_yindex(-1,0);
    iplot_state=2;
    break;
  case GLUT_KEY_UP:
    visy_all=1;
    next_yindex(1,0);
    iplot_state=2;
    break;
  case GLUT_KEY_PAGE_DOWN:
    visz_all=1;
    next_zindex(-1,0);
    iplot_state=3;
    break;
  case GLUT_KEY_PAGE_UP:
    visz_all=1;
    next_zindex(1,0);
    iplot_state=3;
    break;
  case GLUT_KEY_HOME:
    switch (iplot_state){
      case 0:
      case 1:
        next_xindex(0,-1);
        break;
      case 2:
        next_yindex(0,-1);
        break;
      case 3:
        next_zindex(0,-1);
        break;
      default:
        ASSERT(0);
        break;
    }
    break;
  case GLUT_KEY_END:
    switch (iplot_state){
      case 0:
      case 1:
        next_xindex(0,1);
        break;
      case 2:
        next_yindex(0,1);
        break;
      case 3:
        next_zindex(0,1);
        break;
      default:
        ASSERT(0);
        break;
    }
    break;
  default:
    ASSERT(0);
    break;
  }
  if(iplot_state!=0)updateplotslice(iplot_state);
  return;

//  plotstate=getplotstate(STATIC_PLOTS);

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
    break;
  case GLUT_ACTIVE_ALT:
    key_state = KEY_ALT;
    break;
  case GLUT_ACTIVE_SHIFT:
    key_state = KEY_SHIFT;
    break;
  default:
    key_state = KEY_NONE;
    break;
  }
  switch (key){
  case GLUT_KEY_RIGHT:
    switch (key_state){
    case KEY_ALT:
      dx = INC_XY*(*cos_direction_angle);
      dy = INC_XY*(*sin_direction_angle);
      getnewpos(eye_xyz,dx,-dy,0.0,1.0);
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

      if(key_state==KEY_SHIFT)local_speed_factor=4.0;
      getnewpos(eye_xyz,dx,-dy,0.0,local_speed_factor);
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
    if(key_state==KEY_ALT){
      eye_xyz[2] -= INC_Z;
    }
    else{
      float local_speed_factor=1.0;

      if(key_state==KEY_SHIFT)local_speed_factor=4.0;
      dx = INC_XY*(*sin_direction_angle);
      dy = INC_XY*(*cos_direction_angle);
      getnewpos(eye_xyz,-dx,-dy,0.0,local_speed_factor);
    }
    break;
  case GLUT_KEY_UP:
    if(key_state==KEY_ALT){  
      eye_xyz[2] += INC_Z;
    }
    else{
      float local_speed_factor=1.0;
  
      if(key_state==KEY_SHIFT)local_speed_factor=4.0;
      dx = INC_XY*(*sin_direction_angle);
      dy = INC_XY*(*cos_direction_angle);
      getnewpos(eye_xyz,dx,dy,0.0,local_speed_factor);
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

/* ------------------ gmod ------------------------ */

float gmod(float x, float y){
  float returnval;

  if(y==0.0)return 0.0;
  returnval = x - (int)(x/y)*y;
  if(returnval<0.0)returnval+=y;
  return returnval;
}

/* ------------------ UpdateFrame ------------------------ */

void UpdateFrame(float thisinterval, int *changetime, int *redisplay){
  int oldcpuframe;
  float totalcpu;
  int ibenchrate;
  char buffer[256];
  float elapsed_time;

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
      if(itimes==0)bench_starttime=thistime/1000.0;
      if(itimes==nglobal_times-1){
        bench_stoptime=thistime/1000.0;
        ibenchrate=10*((float)nglobal_times/(bench_stoptime-bench_starttime))+0.5;
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
   
    last_frame_count=frame_count;
    frame_count=1;
    lasttime = thistime;
    if(nglobal_times>0){
      *changetime=1;
      if(stept ==1 && plotstate == DYNAMIC_PLOTS && timedrag==0 && RenderGif==0){
        /*  skip frames here if displaying in real time and frame rate is too slow*/
        if(global_times!=NULL&&realtime_flag!=0&&FlowDir>0){
          elapsed_time = (float)thistime/1000.0 - reset_time;
          elapsed_time *= (float)realtime_flag;
          elapsed_time += global_times[reset_frame];
          if(nglobal_times>1&&
            global_times[nglobal_times-1]>global_times[0]&&
            (elapsed_time>global_times[nglobal_times-1]||elapsed_time<0.0)
            ){
            elapsed_time = gmod(elapsed_time,global_times[nglobal_times-1]-global_times[0])+global_times[0];
          }
          itimes = isearch(global_times,nglobal_times,elapsed_time,itimes);
        }
        else{
          if(script_render_flag==0){
            itimes+=FlowDir;
          }
          else{
            itimes=script_itime;
          }
        }
      }
      if(stept==1&&timedrag==0&&RenderGif!=0){
        itimes+=RenderSkip*FlowDir;
      }

// if toggling time display with H then show the frame that was visible

      if(stept==0){
        itime_save=-1;
      }
      else{
        if(itime_save>=0){
          itimes=itime_save;
        }
      }
#ifdef pp_SHOOTER
      if(shooter_firstframe==1&&visShooter!=0&&shooter_active==1){
        reset_itimes0();
      }
#endif
      checktimebound();
      UpdateTimeLabels();
    }
    *redisplay=1;
  }
}

/* ------------------ Idle_CB ------------------------ */

void Idle_CB(void){
  int changetime=0;
  float thisinterval;
  int redisplay=0;

  CheckMemory;
  glutSetWindow(mainwindow_id);
  updateShow();
  thistime = glutGet(GLUT_ELAPSED_TIME);
  thisinterval = thistime - lasttime;
  frame_count++;

  /* increment frame counter if the proper amount of time has passed
     or if we are rendering images or stepping by hand */

  UpdateFrame(thisinterval,&changetime,&redisplay);
  if(showtime==1&&stept==0&&itimeold!=itimes){
    changetime=1;
    checktimebound();
    UpdateTimeLabels();
  }
  update_framenumber(changetime);
  if(redisplay==1){
    glutPostRedisplay();
  }
}

/* ------------------ update_camera_ypos ------------------------ */

void Reshape_CB(int width, int height){
  int screenWidth2, screenHeight2;

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
#ifdef pp_GPUDEPTH
  createDepthTexture();
#endif
#endif
}

/* ------------------ reset_gltime ------------------------ */

void reset_gltime(void){
  int inttime;

  if(showtime!=1)return;
  reset_frame=itimes;
  inttime  = glutGet(GLUT_ELAPSED_TIME);
  reset_time = (float)inttime/1000.0;
  if(global_times!=NULL&&nglobal_times>0){
    start_frametime=global_times[0];
    stop_frametime=global_times[nglobal_times-1];
  }
}

/* ------------------ update_currentmesh ------------------------ */

void update_current_mesh(mesh *meshi){
  current_mesh=meshi;
  loaded_isomesh=get_loaded_isomesh();
  update_iso_showlevels();
  update_plot3dtitle();
}

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

/* ------------------ DoStereo ------------------------ */

int DoStereo(void){
  int return_code=0;
  
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
    return_code=1;
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
    return_code=2;
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
    return_code=3;
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
    return_code=4;
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
    return_code=5;
  }
  return return_code;
}

/* ------------------ DoScript ------------------------ */

void DoScript(void){
  if(runscript==1&&default_script!=NULL){
    ScriptMenu(default_script->id);
    runscript=2;
  }
  script_render_flag=0;
  if(nscriptinfo>0&&current_script_command!=NULL){
    if(current_script_command->command==SCRIPT_VOLSMOKERENDERALL){\
      if(current_script_command->exit==0){
        RenderState(1);
      }
      else{
        RenderState(0);
        current_script_command->first=1;
        current_script_command->exit=0;
      }
    }
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
    else{
      if(current_script_command->command==SCRIPT_VOLSMOKERENDERALL){
        int remove_frame;
  
        script_loadvolsmokeframe2();
        remove_frame=current_script_command->remove_frame;
        if(remove_frame>=0){
          unload_volsmoke_frame_allmeshes(remove_frame);
        }
        glutPostRedisplay();
      }
    }
    glutPostRedisplay();
  }
}

/* ------------------ update_smoothblockage_info ------------------------ */

void update_smoothblockage_info(void){
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

/* ------------------ Display ------------------------ */

void Display_CB(void){
  int dostereo;

  DoScript();
  if(update_glui_dialogs!=0){
    switch (update_glui_dialogs){
      case 1:
        update_glui_dialogs=2;
        break;
      case 2:
        Update_Glui_Dialogs();
        update_glui_dialogs=0;
        update_glui_dialogs=3;
        break;
      case 3:
        Show_Glui_Dialogs();
        update_glui_dialogs=0;
        break;
      default:
        ASSERT(0);
        break;
    }
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
  if(updatemenu==1&&usemenu==1&&menustatus==GLUT_MENU_NOT_IN_USE){
    glutDetachMenu(GLUT_RIGHT_BUTTON);
    InitMenus(LOAD);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    updatemenu=0;
  }
  if(update_fire_colorbar_index==1){
    SmokeColorBarMenu(fire_colorbar_index_ini);
    update_fire_colorbar_index=0;
  }
  if(showtime==0&&ntotal_smooth_blockages>0){
    update_smoothblockage_info();
  }
  glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
  if(showstereo==0){
    dostereo=0;
  }
  else{
    dostereo=DoStereo();
  }
  if(dostereo==0){
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
      else{
        GLubyte *screenbuffers[4];
        int renderdoublenow=0;

        if(RenderOnceNow==1){
          renderdoublenow=1;
        }
    
        if(plotstate==DYNAMIC_PLOTS && nglobal_times>0){
          if(itimes>=0&&itimes<nglobal_times&&
            ((render_frame[itimes] == 0&&showstereo==0)||(render_frame[itimes]<2&&showstereo!=0))
            ){
            render_frame[itimes]++;
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

          if(screenbuffers[0]!=NULL){
            FREEMEMORY(screenbuffers[0]);
          }
          if(screenbuffers[1]!=NULL){
            FREEMEMORY(screenbuffers[1]);
          }
          if(screenbuffers[2]!=NULL){
            FREEMEMORY(screenbuffers[2]);
          }
          if(screenbuffers[3]!=NULL){
            FREEMEMORY(screenbuffers[3]);
          }
        }
        if(renderdoublenow==0||RenderOnceNow==1){
          ASSERT(RenderSkip>0);
          RenderState(0);
          RenderSkip=1;
        }
      }
    }
  }
  if(touring == 1 ){
    if(RenderGif != 0){
      if(nglobal_times>0)angle_global += 2.0*PI/((float)nglobal_times/(float)RenderSkip);
      if(nglobal_times==0)angle_global += 2.0*PI/((float)maxtourframes/(float)RenderSkip);
    }
    if(RenderGif == 0)angle_global += dang_global;
    if(angle_global>PI){angle_global -= -2.0f*PI;}
    if(eyeview==WORLD_CENTERED||eyeview==WORLD_CENTERED_LEVEL){
      float *angle_zx;

      angle_zx = camera_current->angle_zx;
      angle_zx[0] = anglexy0 + angle_global*180./PI;
    }
    else{          
      camera_current->direction_angle = direction_angle0 + angle_global*180./PI;
      camera_current->cos_direction_angle = cos(PI*camera_current->direction_angle/180.0);
      camera_current->sin_direction_angle = sin(PI*camera_current->direction_angle/180.0);
    }
    glutPostRedisplay();
  }
}

/* ------------------ ResizeWindow ------------------------ */

void ResizeWindow(int width, int height){
  float wscaled, hscaled;

  if(render_double!=0)return;
  glutSetWindow(mainwindow_id);
  wscaled = (float)width/(float)max_screenWidth;
  hscaled = (float)height/(float)max_screenHeight;
  if(wscaled>1.0||hscaled>1.0){
    if(wscaled>hscaled){
      width/=wscaled;
      height/=wscaled;
    }
    else{
      width/=hscaled;
      height/=hscaled;
    }
  }
  glutReshapeWindow(width,height);
  glutPostRedisplay();
}


