#include "options.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include GLUT_H

#include "update.h"
#include "smokeviewvars.h"
#include "IOvolsmoke.h"

#ifdef pp_LUA
#include "lua_api.h"
#endif

#undef pp_GPU_CULL_STATE
#ifdef pp_GPU
#define pp_GPU_CULL_STATE
#endif
#ifdef pp_CULL
#ifndef pp_GPU_CULL_STATE
#define pp_GPU_CULL_STATE
#endif
#endif

/* ------------------ get_index ------------------------ */

int get_index(float x, int dir, float *plotxyz, int nplotxyz){
  int i;
  int return_index=0;
  float min_val,vali;

  switch(dir){
    case XDIR:
      x=NORMALIZE_X(x);
      break;
    case YDIR:
      x=NORMALIZE_Y(x);
      break;
    case ZDIR:
      x=NORMALIZE_X(x);
      break;
    default:
      ASSERT(FFALSE);
      break;
  }

  if(plotxyz!=NULL){
    min_val=ABS(x-plotxyz[0]);
    return_index=0;
    for(i=1;i<nplotxyz;i++){
      vali = ABS(x-plotxyz[i]);
      if(vali<min_val){
        return_index=i;
        min_val = vali;
      }
    }
  }
  return return_index;
}

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
        meshdata *meshi;

        slicei = sliceinfo + i;
        if(slicei->loaded==0||slicei->display==0)continue;
        meshi = meshinfo + slicei->blocknumber;
        if(meshi->iplotx_all[iplotx_all]!=-1)return;
      }
      for(i=0;i<nvsliceinfo;i++){
        vslicedata *vslicei;
        meshdata *meshi;

        vslicei = vsliceinfo + i;
        if(vslicei->loaded==0||vslicei->display==0)continue;
        meshi = meshinfo + vslicei->val->blocknumber;
        if(meshi->iploty_all[iploty_all]!=-1)return;
      }
    }
    else{
      for(i=0;i<nplot3dinfo;i++){
        plot3ddata *plot3di;
        meshdata *meshi;

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
        meshdata *meshi;

        slicei = sliceinfo + i;
        if(slicei->loaded==0||slicei->display==0)continue;
        meshi = meshinfo + slicei->blocknumber;
        if(meshi->iploty_all[iploty_all]!=-1)return;
      }
      for(i=0;i<nvsliceinfo;i++){
        vslicedata *vslicei;
        meshdata *meshi;

        vslicei = vsliceinfo + i;
        if(vslicei->loaded==0||vslicei->display==0)continue;
        meshi = meshinfo + vslicei->val->blocknumber;
        if(meshi->iploty_all[iploty_all]!=-1)return;
      }
    }
    else{
      for(i=0;i<nplot3dinfo;i++){
        plot3ddata *plot3di;
        meshdata *meshi;

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
        meshdata *meshi;

        slicei = sliceinfo + i;
        if(slicei->loaded==0||slicei->display==0)continue;
        meshi = meshinfo + slicei->blocknumber;
        if(meshi->iplotz_all[iplotz_all]!=-1)return;
      }
      for(i=0;i<nvsliceinfo;i++){
        vslicedata *vslicei;
        meshdata *meshi;

        vslicei = vsliceinfo + i;
        if(vslicei->loaded==0||vslicei->display==0)continue;
        meshi = meshinfo + vslicei->val->blocknumber;
        if(meshi->iploty_all[iploty_all]!=-1)return;
      }
    }
    else{
      for(i=0;i<nplot3dinfo;i++){
        plot3ddata *plot3di;
        meshdata *meshi;

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
  PRINTF("state=%i\n",state);
  switch(state){
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

  ShowScene(SELECTOBJECT,VIEW_CENTER,0,0,0,NULL);
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
  tour_drag=0;
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
    tour_drag=1;
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

  mouse_x=x; mouse_y=screenHeight-y;
  glDisable(GL_BLEND);
  glDisable(GL_DITHER);
  glDisable(GL_FOG);
  glDisable(GL_LIGHTING);
  glDisable(GL_TEXTURE_1D);
  glDisable(GL_TEXTURE_2D);
  glShadeModel(GL_FLAT);

  ShowScene(SELECTOBJECT,VIEW_CENTER,0,0,0,NULL);
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
    meshdata *meshi;
    selectdata *sd;
    int i;

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

    switch(sd->dir){
      case DOWN_X:
      case UP_X:
        xyz_dir=XDIR;
        break;
      case DOWN_Y:
      case UP_Y:
        xyz_dir=YDIR;
        break;
      case DOWN_Z:
      case UP_Z:
        xyz_dir=ZDIR;
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
    switch(sd->dir){
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
    Update_Blockvals(SELECT_BLOCKS);
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

  ShowScene(SELECTOBJECT,VIEW_CENTER,0,0,0,NULL);
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
      PRINTF("Selected Device: index=%i location:(%f,%f,%f) label:%s\n",val,xyz[0],xyz[1],xyz[2],devicei->labelptr);
    }
    else{
      PRINTF("Selected Device: index=%i location:(%f,%f,%f)\n",val,xyz[0],xyz[1],xyz[2]);
    }
    glShadeModel(GL_SMOOTH);
    glEnable(GL_BLEND);
    glEnable(GL_LIGHTING);
  }
}

/* ------------------ mouse_select_avatar ------------------------ */

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

  ShowScene(SELECTOBJECT,VIEW_CENTER,0,0,0,NULL);
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
  int i;

  if(timebar_drag==0&&itimes>nglobal_times-1||timebar_drag==1&&itimes<0){
    izone=0;
    itimes=first_frame_index;
    if(rendering_status==RENDER_ON){
      RenderMenu(RenderCancel);
      // following exits render command, do this a better way
      if(current_script_command!=NULL)current_script_command->exit=1;
    }
    frame_index=first_frame_index;
    for(i=0;i<nsliceinfo;i++){
      slicedata *sd;

      sd=sliceinfo+i;
      sd->itime=0;
    }
    for(i=0;i<nmeshes;i++){
      meshdata *meshi;

      meshi=meshinfo+i;

      meshi->patch_itime=0;
    }
    for(i=0;i<nmeshes;i++){
      meshdata *meshi;

      meshi=meshinfo+i;
      if(meshi->iso_times==NULL)continue;
      meshi->iso_itime=0;
    }
  }
  if(timebar_drag==0&&itimes<0||timebar_drag==1&&itimes>nglobal_times-1){
    izone=nzone_times-1;
    itimes=nglobal_times-1;
    for(i=0;i<npartinfo;i++){
      partdata *parti;

      parti=partinfo+i;
      parti->itime=parti->ntimes-1;
    }
    for(i=0;i<nsliceinfo;i++){
      slicedata *sd;

      sd=sliceinfo+i;
      sd->itime=sd->ntimes-1;
      if(sd->volslice==1)sd->itime--;
    }
    for(i=0;i<nmeshes;i++){
      meshdata *meshi;

      meshi=meshinfo+i;
      meshi->patch_itime=meshi->npatch_times-1;
    }
    for(i=0;i<nmeshes;i++){
      meshdata *meshi;

      meshi=meshinfo+i;
      if(meshi->iso_times==NULL)continue;
      meshi->iso_itime=meshi->niso_times-1;
    }
  }
  /* set blockage visibility */

  for(i=0;i<nmeshes;i++){
    meshdata *meshi;
    int j;

    meshi=meshinfo+i;
    for(j=0;j<meshi->nbptrs;j++){
      blockagedata *bc;

      bc=meshi->blockageinfoptrs[j];
      if(bc->showtimelist==NULL)continue;
      bc->show=bc->showtimelist[itimes];
    }
  }
}

#define CB_SELECT_CONTINUE -2
#define CB_SELECT_STOP -1
/* ------------------ get_colorbar_index ------------------------ */

int get_colorbar_index(int flag, int x, int y){

  if(flag==0||(colorbar_left_pos<=x&&x<=colorbar_right_pos)){
      y = screenHeight - y;
      if(colorbar_down_pos<=y&&y<=colorbar_top_pos){
        int index;
        index = CLAMP(255*(float)(y-colorbar_down_pos)/(float)(colorbar_top_pos - colorbar_down_pos),0,255);
        return index;
      }
      else{
        return CB_SELECT_STOP;
      }
  }
  return CB_SELECT_CONTINUE;
}

/* ------------------ colorbar_click ------------------------ */

int colorbar_click(int x, int y){
  int colorbar_index;

  colorbar_index = get_colorbar_index(1,x,y);
  if(colorbar_index>=0){
    int state;

    colorbar_select_index=colorbar_index;
    state=glutGetModifiers();
    if(state==GLUT_ACTIVE_CTRL&&current_colorbar!=NULL&&current_colorbar->nsplits==1){
      colorbar_splitdrag=1;
    }
    else{
      colorbar_drag=1;
      UpdateRGBColors(colorbar_index);
    }
    return 1;
  }
  else if(colorbar_index==CB_SELECT_CONTINUE){
    return 0;
  }
  else if(colorbar_index==CB_SELECT_STOP){
    colorbar_drag=0;
    UpdateRGBColors(COLORBAR_INDEX_NONE);
  }
  else{
    ASSERT(FFALSE);
  }
  return 0;
}

/* ------------------ timebar_click ------------------------ */

int timebar_click(int x, int y){
  if(screenHeight-y<VP_timebar.height&&nglobal_times>0){
    int timebar_right_pos;
    int timebar_left_pos;

    timebar_left_pos = VP_timebar.left+timebar_left_width;
    timebar_right_pos=VP_timebar.right-timebar_right_width;

    if(timebar_right_pos>timebar_left_pos){
      itimes = (float)nglobal_times*(float)(x-timebar_left_pos)/(float)(timebar_right_pos-timebar_left_pos);
    }
    else{
      itimes=0;
    }
//    PRINTF("ngt=%i xl=%i x=%i xr=%i\n",nglobal_times,timebar_left_pos,x,timebar_right_pos);
    checktimebound();
    timebar_drag=1;
    stept=0;
    Idle_CB();
    return 1;
  }
  return 0;
}

/* ------------------ update_mouseinfo ------------------------ */

#define MOUSE_DOWN 0
#define MOUSE_UP 1
#define MOUSE_MOTION 2

void update_mouseinfo(int flag, int xm, int ym){
  float thistime_local;
  float dx, dy;
  int maxWH;
  float quat_temp[4];
  float delta_angle,delta_angle2,delta_distance;
  float xymax;
  mousedata *mi;

  mi = &mouseinfo;

  ym = screenHeight - ym; // flip y

  maxWH = MAX(screenWidth,screenHeight);
  xymax = 0.3*(float)MIN(screenWidth,screenHeight)/(float)maxWH;
  switch(flag){
    case MOUSE_DOWN:
      mi->current[0]=xm;
      mi->current[1]=ym;
      mi->last[0]=xm;
      mi->last[1]=ym;
      mi->direction[0]=0;
      mi->direction[1]=0;
      mi->angle = 0.0;
      mi->lastangle=0.0;
      mi->lasttime = glutGet(GLUT_ELAPSED_TIME)/1000.0;
      mi->xcurrent[0] = (float)(mi->current[0]-screenWidth/2)/(float)maxWH;
      mi->xcurrent[1] = (float)(mi->current[1]-screenHeight/2)/(float)maxWH;
      mi->angle=atan2(mi->xcurrent[1],mi->xcurrent[0]);
      delta_angle2 = 0.0;
      mi->lastangle=mi->angle;
      delta_distance=0.0;
      if((ABS(mi->xcurrent[0])<xymax&&ABS(mi->xcurrent[1])<xymax)||key_state==KEY_SHIFT){
        mi->region=0;
      }
      else{
        mi->region=1;
      }
      break;
    case MOUSE_UP:
      delta_angle2=0.0;
      delta_distance=0.0;
      break;
    case MOUSE_MOTION:
      thistime_local=glutGet(GLUT_ELAPSED_TIME)/1000.0;
      mi->current[0]=xm;
      mi->current[1]=ym;
      if(thistime_local-mi->lasttime>0.2){
        mi->direction[0]=mi->current[0]-mi->last[0];
        mi->direction[1]=mi->current[1]-mi->last[1];
        mi->last[0]=mi->current[0];
        mi->last[1]=mi->current[1];
        mi->lasttime=thistime_local;
      }
      mi->xcurrent[0] = (float)(mi->current[0]-screenWidth/2)/(float)maxWH;
      mi->xcurrent[1] = (float)(mi->current[1]-screenHeight/2)/(float)maxWH;
      dx = (float)(mi->current[0]-mi->last[0])/(float)maxWH;
      dy = (float)(mi->current[1]-mi->last[1])/(float)maxWH;
      delta_distance=sqrt(dx*dx+dy*dy);
      mi->angle=atan2(mi->xcurrent[1],mi->xcurrent[0]);
      delta_angle2 = mi->angle - mi->lastangle;
      mi->lastangle = mi->angle;
      break;
    default:
      ASSERT(FFALSE);
      break;
  }

  if(mi->direction[0]==0&&mi->direction[1]==0){
    mi->xdirection[0]=0.0;
    mi->xdirection[1]=0.0;
  }
  else{
    float sum;

    sum = mi->direction[0]*mi->direction[0] + mi->direction[1]*mi->direction[1];
    sum = sqrt(sum);
    mi->xdirection[0]=mi->direction[0]/sum;
    mi->xdirection[1]=mi->direction[1]/sum;
  }
  if(mi->region==0){
    delta_angle = (delta_distance);
    if(delta_angle!=0.0&&ABS(mi->xdirection[0])+ABS(mi->xdirection[1])>0.0){
      float axis[3];

      if(key_state==KEY_SHIFT){
        float w, x, y, z;
        int inside=0;

        w = quat_general[0];
        x = quat_general[1];
        y = quat_general[2];
        z = quat_general[3];
        if(z*z<x*x+y*y)inside=1;
        if(inside==1&&dy<0.0||inside==0&&dx<0.0)delta_angle=-delta_angle;
        axis[0]=2.0*x*z + 2.0*w*y;
        axis[1]=2.0*y*z-2.0*w*x;
        axis[2]=1.0-2*x*x-2*y*y;
      }
      else{
        axis[0]=-mi->xdirection[1]; // direction orthogonal to mouse motion
        axis[1]=0.0;
        axis[2]=mi->xdirection[0];
      }
      angleaxis2quat(delta_angle,axis,quat_temp);
      mult_quat(quat_temp,quat_general,quat_general);
     // level_scene(0,1,quat_general);
      quat2rot(quat_general,quat_rotation);
    }
  }
  else{
    float axis[3];

    axis[0]=0.0;
    axis[1]=1.0;
    axis[2]=0.0;
    angleaxis2quat(-delta_angle2,axis,quat_temp);
    mult_quat(quat_temp,quat_general,quat_general);
   // level_scene(0,1,quat_general);
    quat2rot(quat_general,quat_rotation);
  }
#ifdef _DEBUG
  {
    float denom;

    denom = sqrt(1.0-quat_general[0]*quat_general[0]);
    if(denom>0.0){
      PRINTF("angle= %f axis=(%f,%f,%f)\n",2.0*RAD2DEG*acos(quat_general[0]),quat_general[1]/denom,quat_general[2]/denom,quat_general[3]/denom);
    }
    else{
      PRINTF("angle= %f axis=(*,*,*)\n",2.0*RAD2DEG*acos(quat_general[0]));
    }
  }
#endif
  camera_current->quaternion[0]=quat_general[0];
  camera_current->quaternion[1]=quat_general[1];
  camera_current->quaternion[2]=quat_general[2];
  camera_current->quaternion[3]=quat_general[3];
  camera_current->quat_defined=1;
}

/* ------------------ mouse_CB ------------------------ */

void mouse_CB(int button, int state, int xm, int ym){
  float *eye_xyz;

  if(rotation_type==ROTATION_3AXIS){
    if(state==GLUT_DOWN){
      update_mouseinfo(MOUSE_DOWN,xm,ym);
    }
    else if(state==GLUT_UP){
      update_mouseinfo(MOUSE_UP,xm,ym);
    }
  }

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
    tour_drag=0;
    mouse_down=0;
    show_gslice_normal_keyboard=0;
    eye_xyz0[0]=eye_xyz[0];
    eye_xyz0[1]=eye_xyz[1];
    eye_xyz0[2]=eye_xyz[2];
    update_translate();
    timebar_drag=0;
    colorbar_drag=0;
    colorbar_splitdrag=0;
    glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
    update_trainer_moves();
    return;
  }

  mouse_down=1;

  // check for double click for translating/rotating 3D slice plane

  if(vis_gslice_data==1||show_gslice_triangles==1||show_gslice_triangulation==1){
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
      if(blockageSelect == 1){
        get_geom_dialog_state();
        if(structured_isopen == 1 && unstructured_isopen == 0)mouse_edit_blockage(button, state, xm, ym);
      }
      if(edittour==1&&blockageSelect==0)mouse_edit_tour(button,state,xm,ym);
      if(select_avatar==1)mouse_select_avatar(button,state,xm,ym);
      if(select_device==1)mouse_select_device(button,state,xm,ym);
    }
    glutPostRedisplay();
    if( showtime==1 || showplot3d==1){
      if(colorbar_click(xm,ym)==1)return;
    }
    if(visTimebar==1&&showtime==1){
      if(timebar_click(xm,ym)==1)return;
    }
    CopyCamera(camera_last,camera_current);
    if(canrestorelastview==0){
      updatemenu=1;
      canrestorelastview=1;
      enable_reset_saved_view();
    }
    switch(button){
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
    switch(state){
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
        key_state = KEY_SHIFT;
        start_xyz0[0]=xm;
        start_xyz0[1]=ym;
        touring=0;
        break;
      default:
        key_state = KEY_NONE;
        start_xyz0[0]=xm;
        start_xyz0[1]=ym;
        touring=0;
        break;
    }
    mouse_down_xy0[0]=xm;
    mouse_down_xy0[1]=ym;
  }
  glutPostRedisplay();
  if(blockageSelect == 1){
    get_geom_dialog_state();
    if(structured_isopen == 1 && unstructured_isopen == 0)Display_CB();
  }
}

/* ------------------ Colorbar_Drag ------------------------ */

void Colorbar_Drag(int xm, int ym){
  int colorbar_index;

  colorbar_index = get_colorbar_index(0,xm,ym);
  if(colorbar_index>=0){
    colorbar_select_index=colorbar_index;
    UpdateRGBColors(colorbar_index);
  }
}

/* ------------------ Colorbar_SplitDrag ------------------------ */

void Colorbar_SplitDrag(int xm, int ym){
  int colorbar_index;

  colorbar_index = get_colorbar_index(0,xm,ym);
  if(colorbar_index>=0){
    int ii;

    if(colorbar_index>250)colorbar_index=250;
    if(colorbar_index<5)colorbar_index=5;
    ii=current_colorbar->splits[0];
    current_colorbar->index_node[ii]=colorbar_index;
    current_colorbar->index_node[ii-1]=colorbar_index;
    RemapColorbar(current_colorbar);
    UpdateRGBColors(COLORBAR_INDEX_NONE);
    UpdateColorbarSplits(current_colorbar);
  }
}

/* ------------------ Timebar_Drag ------------------------ */

void Timebar_Drag(int xm, int ym){
  if(nglobal_times>0){
    int timebar_right_pos;
    int timebar_left_pos;

    timebar_left_pos = VP_timebar.left+timebar_left_width;
    timebar_right_pos=VP_timebar.right-timebar_right_width;

    itimes=0;
    if(timebar_right_pos>timebar_left_pos){
      itimes = (float)nglobal_times*(float)(xm-timebar_left_pos)/(float)(timebar_right_pos-timebar_left_pos);
    }
    checktimebound();
    timebar_drag=1;
  }
  Idle_CB();
}

/* ------------------ Drag_Tour_Node ------------------------ */

void Drag_Tour_Node(int xm, int ym){
  float screen_perm[9];

  if(showtour_dialog==1&&edittour==1&&selected_frame!=NULL){
    get_screen_mapping(selected_frame->nodeval.eye,screen_perm);
  }
  else{
    return;
  }
  switch(key_state){
    case KEY_NONE:
    case KEY_ALT:
    case KEY_SHIFT:
    case KEY_CTRL:
      {
        float xy[2]={0.0,0.0},xyz[3];
        float dx, dy, dz;

// scale mouse coordinates

        if(key_state!=KEY_ALT)xy[0] = (float)(xm-mouse_down_xy0[0])/(float)screenWidth;
        if(key_state!=KEY_CTRL)xy[1] = -(float)(ym-mouse_down_xy0[1])/(float)screenHeight;

// permute coordinates according to how screen is rotated

        xyz[0]=DOT2(xy,screen_perm);
        xyz[1]=DOT2(xy,screen_perm+3);
        xyz[2]=DOT2(xy,screen_perm+6);

// compute mouse movement in physical coordinates

        dx =  (xbarORIG-xbar0ORIG)*xyz[0];
        dy =  (ybarORIG-ybar0ORIG)*xyz[1];
        dz =  (zbarORIG-zbar0ORIG)*xyz[2];

// offset tour node location

        tour_xyz[0] += dx;
        tour_xyz[1] += dy;
        tour_xyz[2] += dz;

        mouse_down_xy0[0]=xm;
        mouse_down_xy0[1]=ym;

// update tour data structures with new tour node location

        update_tour_parms();
        update_glui_keyframe();
      }
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
}

/* ------------------ Move_Gen_Slice ------------------------ */

void Move_Gen_Slice(int xm, int ym){
  int dxm, dym;

  dxm = xm - start_xyz0[0];
  dym = ym - start_xyz0[1];
  switch(key_state){
    case KEY_NONE:
      {
        float daz, delev;

        daz = 360.0*dxm/(float)screenWidth;
        delev = 360.0*dym/(float)screenHeight;
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
        xx = xx/(float)screenWidth;
        yy = ym-mouse_down_xy0[1];
        yy = yy/(float)screenHeight;
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
        yy = yy/(float)screenHeight;

        gslice_xyz[2] = gslice_xyz0[2] - DENORMALIZE_Z(4*(xyzbox-NORMALIZE_Z(gslice_xyz0[2]))*yy);
        update_gslice_parms();
      }
      break;
    case KEY_SHIFT:
    break;
    default:
      ASSERT(FFALSE);
      break;
  }
}

/* ------------------ Move_Scene ------------------------ */

void Move_Scene(int xm, int ym){
  float *eye_xyz, *az_elev;
  int dxm, dym;
  float elevation;
  float xx, yy;

  eye_xyz = camera_current->eye;
  az_elev = camera_current->az_elev;

  dxm = xm - start_xyz0[0];
  dym = ym - start_xyz0[1];
  switch(key_state){
    case KEY_NONE:
      switch(rotation_type){
        case ROTATION_3AXIS:
        break;
        case ROTATION_2AXIS:
        case ROTATION_1AXIS:
          az_elev[0] += dxm;
          if(rotation_type==ROTATION_2AXIS){
            az_elev[1] += dym;
          }
          else{
            az_elev[1]=0.0;
          }
          start_xyz0[0]=xm;
          start_xyz0[1]=ym;
          break;
        case EYE_CENTERED:
#define ANGLE_FACTOR 0.25
          camera_current->azimuth += dxm*ANGLE_FACTOR;
          start_xyz0[0]=xm;

          camera_current->elevation -= dym*ANGLE_FACTOR;
          elevation=camera_current->elevation;
          if(elevation>80.0)elevation=80.0;
          if(elevation<-80.0)elevation=-80.0;
          camera_current->elevation=elevation;
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
        xx = xx/(float)screenWidth;
        yy = ym-mouse_down_xy0[1];
        yy = yy/(float)screenHeight;
        if(rotation_type==EYE_CENTERED){
          float xx2, yy2;
          float cs_az, sn_az;

          cs_az = cos(DEG2RAD*camera_current->azimuth);
          sn_az = sin(DEG2RAD*camera_current->azimuth);

          xx2 = cs_az*xx - sn_az*yy;
          yy2 = sn_az*xx + cs_az*yy;
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
      xx = xx/(float)screenWidth;
      yy = ym-mouse_down_xy0[1];
      yy = yy/(float)screenHeight;

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
    PRINTF("MOTION: %4.1f fps\n",fps);
    lastMOTIONtime=thisMOTIONtime;
    MOTIONnframes=0;
  }
  return 0;
}

/* ------------------ motion_CB ------------------------ */

void motion_CB(int xm, int ym){

  in_external=0;
#ifdef pp_GPUTHROTTLE
  if(usegpu==1&&showvolrender==1&&show_volsmoke_moving==1){
    if(throttle_gpu()==1)return;
  }
#endif

  glutPostRedisplay();

  if( colorbar_drag==1&&(showtime==1 || showplot3d==1)){
    Colorbar_Drag(xm,ym);
    return;
  }
  if(
    colorbar_splitdrag==1&&
    (showtime==1 || showplot3d==1)&&current_colorbar!=NULL&&current_colorbar->nsplits==1){
    Colorbar_SplitDrag(xm,ym);
    return;
  }
  if(timebar_drag==1){
    Timebar_Drag(xm,ym);
    return;
  }
  if(move_gslice==1){
    Move_Gen_Slice(xm,ym);
    return;
  }
  if(tour_drag==1){
    Drag_Tour_Node(xm,ym);
    return;
  }

  if(rotation_type==ROTATION_3AXIS&&(key_state == KEY_NONE||key_state == KEY_SHIFT)){
    update_mouseinfo(MOUSE_MOTION,xm,ym);
  }
  if((rotation_type==ROTATION_2AXIS||rotation_type==ROTATION_3AXIS)&&gvec_down==1&&key_state == KEY_NONE){
    update_gvec_down(0);
  }
  Move_Scene(xm,ym);
}

/* ------------------ keyboard_up_CB ------------------------ */

void keyboard_up_CB(unsigned char key, int x, int y){
  resetclock=1;
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
  PRINTF("%s ",gpu_label);
#endif
#ifdef pp_CULL
  PRINTF("%s",cull_label);
#endif
  PRINTF("\n");
}
#endif

/* ------------------ keyboard ------------------------ */

void keyboard(unsigned char key, int flag){
  char key2;
  int skip2;
  meshdata *gbsave,*gbi;
  int i;
  int keystate=0;

  if(flag==FROM_CALLBACK){
    keystate = 6&glutGetModifiers();
    if(scriptoutstream!=NULL&&key!='t'&&key!='r'&&key!='R'&&key!=' '&&key!='-'){
      fprintf(scriptoutstream,"KEYBOARD\n");
      switch(keystate){
        case GLUT_ACTIVE_ALT:
          fprintf(scriptoutstream," ALT %c\n",key);
          break;
        case GLUT_ACTIVE_CTRL:
          fprintf(scriptoutstream," CTRL %c\n",key);
          break;
        default:
          fprintf(scriptoutstream," %c\n",key);
          break;
      }
    }
  }
  else if(flag==FROM_SCRIPT){
    keystate=script_keystate;
  }
  glutPostRedisplay();
  key2 = (char)key;

  switch(key2){
    case 'a':
      if(showtour_dialog==1&&edittour==1){
        add_delete_keyframe(ADD_KEYFRAME);
        break;
      }
      if(rotation_type==EYE_CENTERED){
        handle_move_keys(256+key2);
        break;
      }
      if((visVector==1&&ReadPlot3dFile==1)||showvslice==1||isZoneFireModel==1){
      }
      else{
        break;
      }
      if(isZoneFireModel==1){
        if(keystate==GLUT_ACTIVE_ALT){
          zone_ventfactor /= 1.5;
        }
        else{
          zone_ventfactor *= 1.5;
        }
        PRINTF("zone vent factor: %f\n",zone_ventfactor);
      }
      else{
        if(keystate==GLUT_ACTIVE_ALT){
          vecfactor/=1.5;
        }
        else{
          vecfactor*=1.5;
        }
        PRINTF("vector length factor: %f\n",vecfactor);
        update_glui_vecfactor();
      }
      if(visVector==1&&ReadPlot3dFile==1){
        gbsave=current_mesh;
        for(i=0;i<nmeshes;i++){
          gbi = meshinfo + i;
          if(gbi->plot3dfilenum==-1)continue;
          update_current_mesh(gbi);
          updateplotslice(XDIR);
          updateplotslice(YDIR);
          updateplotslice(ZDIR);
        }
        update_current_mesh(gbsave);
      }
      break;
    case 'A':
      axislabels_smooth=1-axislabels_smooth;
      update_axislabels_smooth();
      break;
    case 'b':
    case 'B':
      switch(keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(DIALOG_BOUNDS); // clip dialog
        break;
      case GLUT_ACTIVE_CTRL:
      default:
        showall_boundary = 1-showall_boundary;
        if(showall_boundary==1){
          ShowPatchMenu(SHOWALL_BOUNDARY);
          if(key2=='B'){
            VentMenu(HIDE_ALL_VENTS);
            BlockageMenu(visBLOCKHide);
          }
          }
        else{
          ShowPatchMenu(HIDEALL_BOUNDARY);
          if(key2=='B'){
            VentMenu(SHOW_ALL_VENTS);
            BlockageMenu(visBLOCKAsInput);
          }
        }
        break;
      }
      break;
    case 'c':
      switch(keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(DIALOG_CLIP); // clip dialog
        break;
      case GLUT_ACTIVE_CTRL:
      default:
        if(nrooms>0){
          zone_highlight_room++;
          if(zone_highlight_room>=nrooms)zone_highlight_room=0;
          PRINTF("room %i\n",zone_highlight_room+1);
        }
        else{
          contour_type++;
          if(contour_type>2)contour_type=0;
          update_plot3d_display();
          UpdateRGBColors(COLORBAR_INDEX_NONE);
        }
      }
      break;
    case 'C':
      switch(keystate){
        case GLUT_ACTIVE_ALT:
          DialogMenu(DIALOG_COLORBAR); // colorbar dialog
          break;
        case GLUT_ACTIVE_CTRL:
        default:
          if(nrooms>0){
            zone_highlight = 1 - zone_highlight;
            if(zone_highlight==1){
              PRINTF("room %i\n",zone_highlight_room+1);
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
      }
      break;
    case 'd':
    case 'D':
      if(key2=='d'&&showtour_dialog==1&&edittour==1){
        add_delete_keyframe(DELETE_KEYFRAME);
        break;
      }
      switch(keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(DIALOG_DISPLAY); // display dialog
        break;
      case GLUT_ACTIVE_CTRL:
      default:
        if(rotation_type==EYE_CENTERED){
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
      switch(keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(DIALOG_GEOMETRY); // edit geometry
        break;
      case GLUT_ACTIVE_CTRL:
      default:
        rotation_type++;
        if(rotation_type>3)rotation_type=0;
        rotation_type_CB(rotation_type);
        update_rotation_type(rotation_type);
        handle_rotation_type(ROTATION_2AXIS);
      }
      break;
    case 'f':
      switch(keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(DIALOG_BOUNDS); // file/bounds dialog
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
      switch(keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(DIALOG_VIEW);
        break;
      case GLUT_ACTIVE_CTRL:
      default:
        if(ntotal_blockages>0||isZoneFireModel==0||(isZoneFireModel==1&&ntrnx>0)){
          switch(visGrid){
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
      if(cache_qdata==1){
        handleiso();
        return;
      }
      break;
    case 'j':
    case 'J':
      if(keystate==GLUT_ACTIVE_ALT){
        sensorrelsize /= 1.5;
      }
      else{
        sensorrelsize *= 1.5;
      }
      update_device_size();
      break;
    case 'k':
    case 'K':
      visTimebar = 1 - visTimebar;
      if(visTimebar==0)PRINTF("Time bar hidden\n");
      if(visTimebar==1)PRINTF("Time bar visible\n");
      break;
#ifdef _DEBUG
    case 'l':
      if(nsmoke3dinfo>0){
        smokecullflag=1-smokecullflag;
        if(smokecullflag==0){
          smokedrawtest=1-smokedrawtest;
        }
        PRINTF("smokecullflag=%i\n smokedrawtest=%i\n",smokecullflag,smokedrawtest);
        update_smoke3dflags();
        return;
      }
      break;
#endif
    case 'L':
      UnloadSliceMenu(UNLOAD_LAST);
      break;
    case 'm':
    case 'M':
      switch(keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(DIALOG_MOTION); // motion dialog
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
        PRINTF("adjustalphaflag=%i\n",adjustalphaflag);
        update_smoke3dflags();
        return;
      }
      break;
#endif
    case 'O':
      switch(visBlocks){
        case visBLOCKAsInput:
        case visBLOCKAsInputOutline:
          BlockageMenu(visBLOCKHide);
          BlockageMenu(visBLOCKOnlyOutline);
          break;
        default:
          BlockageMenu(visBLOCKHide);
          BlockageMenu(visBLOCKAsInput);
          break;
      }
      break;
    case 'o':
      if(keystate==GLUT_ACTIVE_ALT){
        switch(visBlocks){
          case visBLOCKAsInput:
            BlockageMenu(visBLOCKHide);
            BlockageMenu(visBLOCKAsInput);
            BlockageMenu(visBLOCKAddOutline);
            break;
          case visBLOCKAsInputOutline:
            BlockageMenu(visBLOCKHide);
            BlockageMenu(visBLOCKNormal);
            break;
          case visBLOCKNormal:
            BlockageMenu(visBLOCKHide);
            BlockageMenu(visBLOCKOnlyOutline);
            break;
          case visBLOCKOutline:
          case visBLOCKOnlyOutline:
            BlockageMenu(visBLOCKHide);
            break;
          case visBLOCKHide:
            BlockageMenu(visBLOCKHide);
            BlockageMenu(visBLOCKAsInput);
            break;
          default:
            BlockageMenu(visBLOCKHide);
            BlockageMenu(visBLOCKAsInput);
            break;
        }
      }
      else{
        highlight_flag++;
        if(highlight_flag!=0&&visFrame==0){
          visFrame = 1;
          updatefacelists = 1;
          updatemenu = 1;
          glutPostRedisplay();
        }
        if(highlight_flag>2&&noutlineinfo>0)highlight_flag=0;
        if(highlight_flag>1&&noutlineinfo==0)highlight_flag=0;
        PRINTF("outline mode=%i\n",highlight_flag);
      }
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
      if(visiso==1&&cache_qdata==1)updatesurface();
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
      if(showedit_dialog==1&&geomtest_option==NO_TEST){
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

        if(render_360 == 1)render_mode = RENDER_360;
        if (render_mode!=RENDER_360) {
          if (strncmp((const char *)&key2, "r", 1) == 0) {
            render_mode = RENDER_XYSINGLE;
          }
          else {
            render_mode = RENDER_XYMULTI;
          }
        }
        render_times = RENDER_SINGLETIME;

        if(strncmp((const char *)&key2,"R",1)==0|| render_mode == RENDER_360){
          if(nrender_rows==1)nrender_rows=2;
          rflag=1;
        }
        else{
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
                meshdata *meshi;
                volrenderdata *vr;
                int j;
                int framenum;
                float timediffmin;

                meshi = meshinfo + i;
                vr = &meshi->volrenderinfo;
                if(vr->fireslice==NULL||vr->smokeslice==NULL)continue;
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
              meshdata *meshi;
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
              fprintf(scriptoutstream," %i %i 4 %i\n",plotn,visVector,contour_type);
            }
          }
          if(clip_rendered_scene!=0){
            fprintf(scriptoutstream,"RENDERCLIP\n");
            fprintf(scriptoutstream," %i %i %i %i %i\n",
              clip_rendered_scene,render_clip_left,render_clip_right,render_clip_bottom,render_clip_top);
          }
          if(vis_gslice_data==1||show_gslice_triangles==1||show_gslice_triangulation==1||show_gslice_normal){
            fprintf(scriptoutstream,"GSLICEVIEW\n");
            fprintf(scriptoutstream," %i %i %i %i\n",vis_gslice_data,show_gslice_triangles,show_gslice_triangulation,show_gslice_normal);
            fprintf(scriptoutstream,"GSLICEPOS\n");
            fprintf(scriptoutstream," %f %f %f\n",gslice_xyz[0],gslice_xyz[1],gslice_xyz[2]);
            fprintf(scriptoutstream,"GSLICEORIEN\n");
            fprintf(scriptoutstream," %f %f\n",gslice_normal_azelev[0],gslice_normal_azelev[1]);
          }
          if(rflag==0){
            fprintf(scriptoutstream,"RENDERONCE\n");
          }
          else{
            fprintf(scriptoutstream,"RENDERDOUBLEONCE\n");
          }
          fprintf(scriptoutstream," %s\n",script_renderfile);
        }
        RenderState(RENDER_ON);
        render_from_menu=0;
      }
      break;
    case 's':
      switch(keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(DIALOG_3DSMOKE); // 3d smoke dialog
        break;
      case GLUT_ACTIVE_CTRL:
        snap_scene();
        break;
      default:
        if(rotation_type==EYE_CENTERED){
          handle_move_keys(GLUT_KEY_DOWN);
        }
        else{
          vectorskip++;
          if(vectorskip>4)vectorskip=1;
        }
      }
      break;
    case 'S':
      stereotypeOLD=stereotype;
      stereotype++;
      if(stereotype>5)stereotype=0;
      if(stereotype==STEREO_TIME&&videoSTEREO!=1)stereotype=STEREO_LR;
      Update_Glui_Stereo();
      break;
    case 't':
      switch(keystate){
      case GLUT_ACTIVE_ALT:
        DialogMenu(DIALOG_TOUR); // tour dialog
        break;
      case GLUT_ACTIVE_CTRL:
      default:
        stept=(stept+1)%2;
        if(stept==1){
          plotstate=GetPlotState(DYNAMIC_PLOTS);
          if(plotstate==DYNAMIC_PLOTS){
            reset_gltime();
          }
          else{
            stept=0;
          }
        }
        if(stept==0)itime_save = -1;
      }
      break;
    case 'T':
      usetexturebar=1-usetexturebar;
      PRINTF("usetexturebar=%i\n",usetexturebar);
      break;
    case 'u':
    case 'U':
      switch(keystate){
        case GLUT_ACTIVE_ALT:
          skip_slice_in_embedded_mesh = 1 - skip_slice_in_embedded_mesh;
          break;
        default:
          ReloadMenu(RELOAD_NOW);
          break;
      }
      break;
    case 'v':
      switch(keystate){
        case GLUT_ACTIVE_ALT:
          projection_type = 1 - projection_type;
          Motion_CB(PROJECTION);
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
      switch(keystate){
        case GLUT_ACTIVE_ALT:
          DialogMenu(DIALOG_WUI); // WUI dialog
          break;
        case GLUT_ACTIVE_CTRL:
        default:
          if(rotation_type==EYE_CENTERED){
            handle_move_keys(GLUT_KEY_UP);
          }
          else{
            vis_gslice_data = 1 - vis_gslice_data;
            update_gslice_parms();
          }
          break;
      }
      break;
    case 'W':
      clip_mode++;
      if(clip_mode>3)clip_mode=0;
      update_clip_all();
      break;
    case 'x':
    case 'X':
      if(keystate==GLUT_ACTIVE_ALT){
        DialogMenu(DIALOG_HIDEALL); // close all dialogs
      }
      else{
        visx_all=1-visx_all;
        plotstate = GetPlotState(STATIC_PLOTS);
        updatemenu = 1;
      }
      break;
    case 'y':
    case 'Y':
      visy_all = 1-visy_all;
      plotstate = GetPlotState(STATIC_PLOTS);
      updatemenu = 1;
      break;
    case 'z':
    case 'Z':
      if(keystate==GLUT_ACTIVE_ALT){
        DialogMenu(DIALOG_SMOKEZIP); // compress dialog
      }
      else{
        visz_all = 1 - visz_all;
        plotstate = GetPlotState(STATIC_PLOTS);
        updatemenu = 1;
      }
      break;
    case '0':
      if(plotstate==DYNAMIC_PLOTS){
        UpdateTimes();
        reset_time_flag=1;
        return;
      }
      break;
    case '~':
      level_scene(1,1,quat_general);
      quat2rot(quat_general,quat_rotation);
      break;
    case '!':
      snap_scene();
      break;
    case '@':
      cell_center_text = 1 - cell_center_text;
      break;
    case '#':
      writeini(LOCAL_INI,NULL);
      break;
    case '$':
      trainer_active=1-trainer_active;
      if(trainer_active==1){
        PRINTF("Trainer mode active\n");
        trainer_mode=1;
        show_glui_trainer();
      }
      if(trainer_active==0){
        PRINTF("Trainer mode inactive\n");
        trainer_mode=0;
        hide_glui_trainer();
      }
      break;
    case '%':
      script_step=1-script_step;
      break;
    case '^':
      script_step_now=1;
      break;
    case '&':
      antialiasflag=1-antialiasflag;
      PRINTF("antialiasflag=%i\n",antialiasflag);
      break;
    case '*':
      visx_all=0;
      visy_all=0;
      visz_all=0;
      vis_gslice_data=0;
      show_gslice_triangles=0;
      show_gslice_triangulation=0;
      show_gslice_normal=0;
      break;
    case '(':
      clip_rendered_scene=1-clip_rendered_scene;
      break;
    case '[':
      edittour=1;
      update_edit_tour();
      break;
    case ']':
      edittour=0;
      update_edit_tour();
      break;
    case ';':
      ColorbarMenu(COLORBAR_FLIP);
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

  if(clip_mode!=CLIP_OFF&&(
    strncmp((const char *)&key2,"<",1)==0||strncmp((const char *)&key2,",",1)==0||
    strncmp((const char *)&key2,">",1)==0||strncmp((const char *)&key2,".",1)==0)){

    if(strncmp((const char *)&key2,"<",1)==0||strncmp((const char *)&key2,",",1)==0){ClipDir=-1;}
     else if(strncmp((const char *)&key2,">",1)==0||strncmp((const char *)&key2,".",1)==0){ClipDir=1;}

    if(stepclip_xmin==1  )clip_i += skip_global*ClipDir;
    if(stepclip_ymin==1  )clip_j += skip_global*ClipDir;
    if(stepclip_zmin==1  )clip_k += skip_global*ClipDir;
    if(stepclip_xmax==1  )clip_I += skip_global*ClipDir;
    if(stepclip_ymax==1  )clip_J += skip_global*ClipDir;
    if(stepclip_zmax==1  )clip_K += skip_global*ClipDir;

    UpdateClipbounds(clipinfo.clip_xmin,&clip_i,clipinfo.clip_xmax,&clip_I,current_mesh->ibar);
    UpdateClipbounds(clipinfo.clip_ymin,&clip_j,clipinfo.clip_ymax,&clip_J,current_mesh->jbar);
    UpdateClipbounds(clipinfo.clip_zmin,&clip_k,clipinfo.clip_zmax,&clip_K,current_mesh->kbar);
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
    if(timebar_drag==0)itimes += skip_global*FlowDir;
    checktimebound();
    Idle_CB();

    return;
  }
  switch(iplot_state){
    case XDIR:
      next_xindex(skip_global*FlowDir,0);
      break;
    case 0:
    case YDIR:
      next_yindex(skip_global*FlowDir,0);
      break;
    case ZDIR:
      next_zindex(skip_global*FlowDir,0);
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
  if(ReadPlot3dFile==1){
    plotstate = GetPlotState(STATIC_PLOTS);
    if(visiso!=0&&current_mesh->slicedir==ISO){
      plotiso[plotn-1] += FlowDir;
      updatesurface();
    }
    glutPostRedisplay();
  }
  if(iplot_state!=0)updateplotslice(iplot_state);
}

/* ------------------ keyboard_CB ------------------------ */

void keyboard_CB(unsigned char key, int x, int y){
  keyboard(key,FROM_CALLBACK);
  glutPostRedisplay();
  updatemenu=1;
}

/* ------------------ handle_rotation_type ------------------------ */

void handle_rotation_type(int flag){
  float *az_elev;

  if(rotation_type==rotation_type_old)return;
  camera_current->rotation_type=rotation_type;
  az_elev = camera_current->az_elev;
  updatemenu=1;
  switch(rotation_type){
  case ROTATION_3AXIS:
      if(trainer_mode==0)PRINTF("Scene centered (3 axis rotation)\n");
      if(showtrainer_dialog==0&&flag==ROTATION_2AXIS&&rotation_type_old==EYE_CENTERED){
        SetViewPoint(RESTORE_EXTERIOR_VIEW);
      }
      break;
  case ROTATION_2AXIS:
      if(trainer_mode==0)PRINTF("Scene centered (2 axis rotation)\n");
      if(showtrainer_dialog==0&&flag==ROTATION_2AXIS&&rotation_type_old==EYE_CENTERED){
        SetViewPoint(RESTORE_EXTERIOR_VIEW);
      }
      break;
  case EYE_CENTERED:
       az_elev[1]=0.0;
       if(showtrainer_dialog==0&&flag==ROTATION_2AXIS&&rotation_type_old!=EYE_CENTERED){
         SetViewPoint(RESTORE_EXTERIOR_VIEW);
       }
      if(trainer_mode==0)PRINTF("eye centered\n");
      break;
  case ROTATION_1AXIS:
    az_elev[1]=0.0;
    if(trainer_mode==0)PRINTF("Scene centered (level rotation)\n");
    if(showtrainer_dialog==0&&flag==ROTATION_2AXIS&&rotation_type_old==EYE_CENTERED){
      SetViewPoint(RESTORE_EXTERIOR_VIEW);
    }
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  showhide_translate(rotation_type);
  rotation_type_old = rotation_type;
  return;
}

/* ------------------ update_clipplanes ------------------------ */

void update_clipplanes(void){
  if(trainer_mode==0){
    if(clip_mode!=clip_mode_last){
      if(clip_mode==CLIP_OFF)PRINTF("clipping off\n");
      if(clip_mode==CLIP_BLOCKAGES_DATA)PRINTF("clipping blockages + data\n");
      if(clip_mode==CLIP_BLOCKAGES)PRINTF("clipping blockages\n");
      if(clip_mode==CLIP_DATA)PRINTF("clipping data\n");
      clip_mode_last=clip_mode;
    }
  }
  if(clip_mode==CLIP_OFF){
    setClipPlanes(NULL,CLIP_OFF);
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

/* ------------------ specialkeyboard_up_CB ------------------------ */

void specialkeyboard_up_CB(int key, int x, int y){
  resetclock=1;
}

/* ------------------ specialkeyboard_CB ------------------------ */

void specialkeyboard_CB(int key, int x, int y){

#define EYE_MODE 0
#define P3_MODE 1
  int keymode=EYE_MODE;

  glutPostRedisplay();

  switch(cursorPlot3D){
    case 0:
      if(rotation_type==EYE_CENTERED){
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

  switch(keymode){
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
  switch(key){
  case GLUT_KEY_LEFT:
    visx_all=1;
    next_xindex(-1,0);
    iplot_state=XDIR;
    break;
  case GLUT_KEY_RIGHT:
    visx_all=1;
    next_xindex(1,0);
    iplot_state=XDIR;
    break;
  case GLUT_KEY_DOWN:
    visy_all=1;
    next_yindex(-1,0);
    iplot_state=YDIR;
    break;
  case GLUT_KEY_UP:
    visy_all=1;
    next_yindex(1,0);
    iplot_state=YDIR;
    break;
  case GLUT_KEY_PAGE_DOWN:
    visz_all=1;
    next_zindex(-1,0);
    iplot_state=ZDIR;
    break;
  case GLUT_KEY_PAGE_UP:
    visz_all=1;
    next_zindex(1,0);
    iplot_state=ZDIR;
    break;
  case GLUT_KEY_HOME:
    switch(iplot_state){
      case 0:
      case XDIR:
        next_xindex(0,-1);
        break;
      case YDIR:
        next_yindex(0,-1);
        break;
      case ZDIR:
        next_zindex(0,-1);
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
    break;
  case GLUT_KEY_END:
    switch(iplot_state){
      case 0:
      case XDIR:
        next_xindex(0,1);
        break;
      case YDIR:
        next_yindex(0,1);
        break;
      case ZDIR:
        next_zindex(0,1);
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  if(iplot_state!=0)updateplotslice(iplot_state);
  return;

//  plotstate=GetPlotState(STATIC_PLOTS);

}

/* ------------------ handle_move_keys ------------------------ */

void handle_move_keys(int  key){
  int state;
  float dx, dy;
  float *azimuth, cos_azimuth, sin_azimuth;
  float *elevation;

  const float INC_ANGLE0=0.1;

  float INC_XY, INC_Z, INC_ANGLE;

  float *eye_xyz;
#define LOOKANGLE_CHANGE 11.25


  eye_xyz = camera_current->eye;

  azimuth=&camera_current->azimuth;
  cos_azimuth=cos(DEG2RAD*(*azimuth));
  sin_azimuth=sin(DEG2RAD*(*azimuth));

  elevation=&camera_current->elevation;

  glui_move_mode=-1;

  INC_XY=SCALE2SMV(meshinfo->cellsize);
  INC_Z=INC_XY;
  INC_ANGLE = 5*INC_ANGLE0;

  state=glutGetModifiers();
  switch(state){
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
  switch(key){
    case GLUT_KEY_RIGHT:
      switch(key_state){
        case KEY_ALT:
          dx = INC_XY*(cos_azimuth);
          dy = INC_XY*(sin_azimuth);
          getnewpos(eye_xyz,dx,-dy,0.0,1.0);
          break;
        case KEY_SHIFT:
        case KEY_CTRL:
        default:
          if(key_state==KEY_SHIFT){
            *azimuth += 4.0*INC_ANGLE;
          }
          else{
            *azimuth += INC_ANGLE;
          }
          break;
      }
      break;
    case 256+'d':
      dx = INC_XY*(cos_azimuth);
      dy = INC_XY*(sin_azimuth);
      {
        float local_speed_factor=1.0;

        if(key_state==KEY_SHIFT)local_speed_factor=4.0;
        getnewpos(eye_xyz,dx,-dy,0.0,local_speed_factor);
      }
      break;
    case GLUT_KEY_LEFT:
      switch(key_state){
        case KEY_ALT:
          dx = INC_XY*(cos_azimuth);
          dy = INC_XY*(sin_azimuth);
          getnewpos(eye_xyz,-dx,dy,0.0,1.0);
          break;
        case KEY_SHIFT:
        case KEY_CTRL:
        default:
          if(key_state==KEY_SHIFT){
            *azimuth -= 4.0*INC_ANGLE;
          }
          else{
            *azimuth -= INC_ANGLE;
          }
          break;
      }
      break;
    case 256+'a':
      dx = INC_XY*(cos_azimuth);
      dy = INC_XY*(sin_azimuth);
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
        dx = INC_XY*(sin_azimuth);
        dy = INC_XY*(cos_azimuth);
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
        dx = INC_XY*(sin_azimuth);
        dy = INC_XY*(cos_azimuth);
        getnewpos(eye_xyz,dx,dy,0.0,local_speed_factor);
      }
      break;
    case GLUT_KEY_PAGE_UP:
      *elevation += LOOKANGLE_CHANGE;
      break;
    case GLUT_KEY_HOME:
      *elevation=0.0;
      break;
    case GLUT_KEY_INSERT:
    case GLUT_KEY_PAGE_DOWN:
      *elevation-=LOOKANGLE_CHANGE;
      break;
    case GLUT_KEY_END:
      SetViewPoint(RESTORE_EXTERIOR_VIEW);
      break;
    case GLUT_KEY_F4:
      camera_current->view_angle-=LOOKANGLE_CHANGE;
      if(camera_current->view_angle<0.0)camera_current->view_angle+=360.0;
      break;
    case GLUT_KEY_F5:
      camera_current->view_angle=0.0;
      break;
    case GLUT_KEY_F6:
      camera_current->view_angle+=LOOKANGLE_CHANGE;
      if(camera_current->view_angle>360.0)camera_current->view_angle-=360.0;
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
  if(rotation_type==EYE_CENTERED){
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
  float elapsed_time;

  if(showtime==1&&((stept==1&&(float)thisinterval>frameinterval)||rendering_status==RENDER_ON||timebar_drag==1)){       /* ready for a new frame */
    cputimes[cpuframe]=thistime/1000.;

    oldcpuframe=cpuframe-10;
    if(oldcpuframe<0)oldcpuframe+=20;
    totalcpu=cputimes[cpuframe]-cputimes[oldcpuframe];
    if(totalcpu==0.0){
      framerate=0.0;
    }
    else{
     framerate=10.0/totalcpu;
    }
    cpuframe++;
    if(cpuframe>=20)cpuframe=0;

    last_frame_count=frame_count;
    frame_count=1;
    lasttime = thistime;
    if(nglobal_times>0){
      *changetime=1;
      if(stept ==1 && plotstate == DYNAMIC_PLOTS && timebar_drag==0 && rendering_status==RENDER_OFF){
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
          itimes = ISearch(global_times,nglobal_times,elapsed_time,itimes);
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
      if(stept==1&&timebar_drag==0&&rendering_status==RENDER_ON){
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
      if(shooter_firstframe==1&&visShooter!=0&&shooter_active==1){
        ResetItimes0();
      }
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
  UpdateShow();
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
  UpdateFrameNumber(changetime);
  if(redisplay==1){
    glutPostRedisplay();
  }
}

/* ------------------ setScreenSize ------------------------ */

void setScreenSize(int *width, int *height){
  if(width!=NULL){
    screenWidth=*width;
    if(screenWidth%2==1)screenWidth++;
  }
  if(height!=NULL){
    screenHeight=*height;
  }
}

/* ------------------ Reshape_CB ------------------------ */

void Reshape_CB(int width, int height){
  if(disable_reshape==1)return;
  updatemenu=1;
  window_aspect_ratio = (float)width/(float)height;
  if(window_aspect_ratio<1.0)window_aspect_ratio=1.0/window_aspect_ratio;
  setScreenSize(&width,&height);
  windowresized=1;
  UpdateCameraYpos(camera_external);
  if(strcmp(camera_current->name,"external")==0&&in_external==1){
    SetViewPoint(RESTORE_EXTERIOR_VIEW);
  }
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

/* ------------------ update_current_mesh ------------------------ */

void update_current_mesh(meshdata *meshi){
  current_mesh=meshi;
  loaded_isomesh=get_loaded_isomesh();
  update_iso_showlevels();
  update_plot3dtitle();
}

/* ------------------ ClearBuffers ------------------------ */

void ClearBuffers(int mode){
  if(mode==DRAWSCENE){
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

  if(stereotype==STEREO_TIME&&videoSTEREO==1){  // temporal stereo (shuttered glasses)
    glDrawBuffer(GL_BACK_LEFT);
    if(stereotype_frame==LEFT_EYE||stereotype_frame==BOTH_EYES){
      ShowScene(DRAWSCENE,VIEW_LEFT,0,0,0,NULL);
    }
    Render(VIEW_LEFT);
    glDrawBuffer(GL_BACK_RIGHT);
    if(stereotype_frame==RIGHT_EYE||stereotype_frame==BOTH_EYES){
      ShowScene(DRAWSCENE,VIEW_RIGHT,0,0,0,NULL);
    }
    Render(VIEW_RIGHT);
    if(buffertype==DOUBLE_BUFFER)glutSwapBuffers();
    return_code=1;
  }
  else if(stereotype==STEREO_LR){             // left/right stereo
    int i;
    int nscreens;

    glDrawBuffer(GL_BACK);
    ClearBuffers(DRAWSCENE);

    nscreens = 1;
    if(render_mode == RENDER_360&&rendering_status==RENDER_ON){
      nscreens = nscreeninfo;
      if (screeninfo == NULL || update_screeninfo == 1)setup_screeninfo();
    }

    for(i = 0; i < nscreens; i++){
      screendata *screeni;

      screeni = NULL;
      if(render_mode == RENDER_360 && rendering_status == RENDER_ON)screeni = screeninfo + i;
      if(stereotype_frame==LEFT_EYE||stereotype_frame==BOTH_EYES){
        int screenWidth_save;

        screenWidth_save=screenWidth;
        screenWidth/=2;
        ShowScene(DRAWSCENE,VIEW_LEFT,0,0,0,screeni);
        screenWidth=screenWidth_save;
      }
      if(stereotype_frame==RIGHT_EYE||stereotype_frame==BOTH_EYES){
        int screenWidth_save;

        screenWidth_save=screenWidth;
        screenWidth/=2;
        ShowScene(DRAWSCENE,VIEW_RIGHT,0,screenWidth,0,screeni);
        screenWidth=screenWidth_save;
      }
      if(render_mode == RENDER_360 && rendering_status == RENDER_ON)screeni->screenbuffer = getscreenbuffer();
      if (buffertype == DOUBLE_BUFFER)glutSwapBuffers();
    }
    if(rendering_status == RENDER_ON){
      if(render_mode == RENDER_360){
        MergeRenderScreenBuffers360();
        for(i = 0; i < nscreeninfo; i++){
          screendata *screeni;

          screeni = screeninfo + i;
          FREEMEMORY(screeni->screenbuffer);
        }
      }
      else{
        Render(VIEW_CENTER);
      }
    }
    return_code=2;
  }
  else if(stereotype==STEREO_RB){             // red/blue stereo
    glDrawBuffer(GL_BACK);
    glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
    glClearColor(1.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if(stereotype_frame==LEFT_EYE||stereotype_frame==BOTH_EYES){
      glColorMask(GL_TRUE,GL_FALSE,GL_FALSE, GL_TRUE);
      ShowScene(DRAWSCENE,VIEW_LEFT,0,0,0,NULL);
      glFlush();
    }

    if(stereotype_frame==RIGHT_EYE||stereotype_frame==BOTH_EYES){
      glDrawBuffer(GL_BACK);
      glColorMask(GL_FALSE,GL_FALSE,GL_TRUE,GL_TRUE);
      glClearColor(0.0, 0.0, 1.0, 1.0);
      glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

      ShowScene(DRAWSCENE,VIEW_RIGHT,0,0,0,NULL);
      glFlush();
    }
    Render(VIEW_CENTER);
    if(buffertype==DOUBLE_BUFFER)glutSwapBuffers();
    return_code=3;
  }
  else if(stereotype==STEREO_RC){             // red/cyan stereo
    glDrawBuffer(GL_BACK);
    glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
    glClearColor(1.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if(stereotype_frame==LEFT_EYE||stereotype_frame==BOTH_EYES){
      glColorMask(GL_TRUE,GL_FALSE,GL_FALSE, GL_TRUE);
      ShowScene(DRAWSCENE,VIEW_LEFT,0,0,0,NULL);
      glFlush();
    }

    if(stereotype_frame==RIGHT_EYE||stereotype_frame==BOTH_EYES){
      glDrawBuffer(GL_BACK);
      glColorMask(GL_FALSE,GL_TRUE,GL_TRUE,GL_TRUE);
      glClearColor(0.0, 1.0, 1.0, 0.0);
      glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

      ShowScene(DRAWSCENE,VIEW_RIGHT,0,0,0,NULL);
      glFlush();
    }
    Render(VIEW_CENTER);
    if(buffertype==DOUBLE_BUFFER)glutSwapBuffers();
    return_code=4;
  }
  else if(stereotype==STEREO_CUSTOM){             // custom red/blue stereo
    glDrawBuffer(GL_BACK);
    glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    if(stereotype_frame==LEFT_EYE||stereotype_frame==BOTH_EYES){
      glColorMask(GL_TRUE,GL_FALSE,GL_FALSE, GL_TRUE);
      ShowScene(DRAWSCENE,VIEW_LEFT,0,0,0,NULL);
      glFlush();
    }
    if(stereotype_frame==RIGHT_EYE||stereotype_frame==BOTH_EYES){
      glDrawBuffer(GL_BACK);
      glColorMask(GL_FALSE,GL_TRUE,GL_TRUE,GL_TRUE);
      glClearColor(0.0, 1.0, 1.0, 1.0);
      glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

      ShowScene(DRAWSCENE,VIEW_RIGHT,0,0,0,NULL);

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
    Render(VIEW_CENTER);
    if(buffertype==DOUBLE_BUFFER)glutSwapBuffers();
    glEnable(GL_LIGHTING);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_DITHER);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    return_code=5;
  }
  return return_code;
}

#ifdef pp_LUA
/* ------------------ DoScriptLua ------------------------ */
void DoScriptLua(void) {
  int script_return_code;
  if(runluascript == 1) {
    if(!luascript_loaded && strlen(luascript_filename)>0)
      load_script(luascript_filename);
    runluascript = 0;
    PRINTF("running lua script section\n");
    fflush(stdout);
    script_return_code = runLuaScript();
    if(script_return_code != LUA_OK && script_return_code != LUA_YIELD && exit_on_script_crash) {
        exit(1);
    }
  }
}
#endif

/* ------------------ DoScript ------------------------ */
#ifdef pp_LUA_SSF
void DoScript(void){
  int script_return_code;
  if(runscript == 1) {
      runscript = 0;
      PRINTF("running ssf script instruction\n");
      fflush(stdout);
      script_return_code = runSSFScript();
      if(script_return_code != LUA_OK && script_return_code != LUA_YIELD && exit_on_script_crash) {
          exit(1);
      }
    }
}
#else
void DoScript(void){
  if(runscript==1&&default_script!=NULL){
    ScriptMenu(default_script->id);
    runscript=2;
  }
  script_render_flag=0;
  if(nscriptinfo>0&&current_script_command!=NULL&&(script_step==0||(script_step==1&&script_step_now==1))){
    script_step_now=0;
#ifndef WIN32
    if(file_exists(stop_filename)){
      fprintf(stderr,"*** Warning: stop file found.  Remove before running smokeview script\n");
      exit(0);
    }
#endif
    if(current_script_command->command==SCRIPT_VOLSMOKERENDERALL){\
      if(current_script_command->exit==0){
        RenderState(RENDER_ON);
      }
      else{
        RenderState(RENDER_OFF);
        current_script_command->first=1;
        current_script_command->exit=0;
      }
    }
    if(current_script_command->command==SCRIPT_ISORENDERALL){\
      if(current_script_command->exit==0){
        RenderState(RENDER_ON);
      }
      else{
        RenderState(RENDER_OFF);
        current_script_command->first=1;
        current_script_command->exit=0;
      }
    }
    if(rendering_status==RENDER_OFF){  // don't advance command if Smokeview is executing a RENDERALL command
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
      }
      if(current_script_command->command==SCRIPT_ISORENDERALL){
        int remove_frame;

        script_loadisoframe2(current_script_command);
        remove_frame = current_script_command->remove_frame;
        if(remove_frame>=0){
          //unload_volsmoke_frame_allmeshes(remove_frame);
        }
      }
    }
    glutPostRedisplay();
  }
  else{
    first_frame_index=0;
    skip_render_frames=0;
    script_startframe=-1;
    script_skipframe=-1;
  }
}
#endif

/* ------------------ Display_CB ------------------------ */

void Display_CB(void){
  int dostereo;

  DoScript();
#ifdef pp_LUA
  DoScriptLua();
#endif
  UpdateDisplay();
  glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
  if(stereotype==STEREO_NONE){
    dostereo=0;
  }
  else{
    dostereo=DoStereo();
  }
  if(dostereo==0){
    if(render_mode == RENDER_XYSINGLE||rendering_status==RENDER_OFF){
      glDrawBuffer(GL_BACK);
      ShowScene(DRAWSCENE,VIEW_CENTER,0,0,0,NULL);
      if(render_mode != RENDER_360)Render(VIEW_CENTER);
      if(buffertype==DOUBLE_BUFFER)glutSwapBuffers();
    }
    else{
      int stop_rendering;

      stop_rendering = 1;
      if(plotstate==DYNAMIC_PLOTS && nglobal_times>0){
        if(itimes>=0&&itimes<nglobal_times&&
          ((render_frame[itimes] == 0&&stereotype==STEREO_NONE)||(render_frame[itimes]<2&&stereotype!=STEREO_NONE))
          ){
          render_frame[itimes]++;
          stop_rendering = 0;
        }
      }
      if(rendering_status == RENDER_ON&&render_mode==RENDER_XYMULTI){
        int nrender_cols;
        int i,ibuffer=0;
        GLubyte **screenbuffers;

        NewMemory((void **)&screenbuffers,nrender_rows*nrender_rows*sizeof(GLubyte *));

        glDrawBuffer(GL_BACK);

        nrender_cols=nrender_rows;
        for(i=0;i<nrender_rows;i++){
          int j;

          for(j=0;j<nrender_cols;j++){
            ShowScene(DRAWSCENE,VIEW_CENTER,1,j*screenWidth,i*screenHeight,NULL);
            screenbuffers[ibuffer++]=getscreenbuffer();
            if(buffertype==DOUBLE_BUFFER)glutSwapBuffers();
          }
        }

        MergeRenderScreenBuffers(nrender_rows,screenbuffers);

        for(i=0;i<nrender_rows*nrender_cols;i++){
          FREEMEMORY(screenbuffers[i]);
        }
        FREEMEMORY(screenbuffers);
      }
      if (render_mode == RENDER_360){
        int i;

        glDrawBuffer(GL_BACK);

        if (screeninfo == NULL||update_screeninfo==1)setup_screeninfo();

        for(i = 0; i < nscreeninfo; i++){
          screendata *screeni;

          screeni = screeninfo + i;
          ShowScene(DRAWSCENE, VIEW_CENTER, 0, 0, 0, screeni);
          screeni->screenbuffer = getscreenbuffer();
          if (buffertype == DOUBLE_BUFFER)glutSwapBuffers();
        }
        MergeRenderScreenBuffers360();

        for(i = 0; i < nscreeninfo; i++){
          screendata *screeni;

          screeni = screeninfo + i;
          FREEMEMORY(screeni->screenbuffer);
        }
      }
      if(stop_rendering==1){
        ASSERT(RenderSkip>0);
        RenderState(RENDER_OFF);
        RenderSkip=1;
      }
    }
  }
}

/* ------------------ ResizeWindow ------------------------ */

void ResizeWindow(int width, int height){
  float wscaled, hscaled;

  if(render_mode != RENDER_XYSINGLE)return;
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
