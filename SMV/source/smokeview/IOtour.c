// $Date$ 
// $Revision$
// $Author$

// svn revision character string
char IOtour_revision[]="$Revision$";

#include "options.h"
#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "update.h"
#include "smokeviewvars.h"

void drawcir(float *center, float rad, float *color);
float hermiteeye(float f1, int i, keyframe *kf1, keyframe *kf2, float *slope);
float hermiteview(float t, int i, keyframe *kf1, keyframe *kf2, float *slope);
void draw_SVOBJECT(sv_object *object, int frame_index_local,propdata *prop,int recurse_level);

/* ------------------ freetour ------------------------ */

void freetour(tourdata *touri){
  int i;
  keyframe *framei;

  for(i=0;i<touri->nkeyframes;i++){
    framei = touri->keyframe_list[i];
    FREEMEMORY(framei);
  }
  FREEMEMORY(touri->pathnodes);
  FREEMEMORY(touri->keyframe_list);
  FREEMEMORY(touri->keyframe_times);
  FREEMEMORY(touri->timeslist);
  FREEMEMORY(touri->path_times);
}

/* ------------------ inittour ------------------------ */

void inittour(tourdata *touri){
  touri->glui_avatar_index=0;
  touri->display2=0;
  touri->global_tension=0.0;
  touri->global_tension_flag=1;
  touri->display=0;
  touri->periodic=0;
  touri->first_frame.prev=NULL;
  touri->first_frame.next=&touri->last_frame;
  touri->first_frame.noncon_time=-1000000000.0;
  touri->first_frame.disp_time=touri->first_frame.noncon_time;

  touri->last_frame.prev=&touri->first_frame;
  touri->last_frame.next=NULL;
  touri->last_frame.noncon_time=1000000000.0;
  touri->last_frame.disp_time=touri->last_frame.noncon_time;

  touri->nkeyframes=0;
  touri->ntimes=view_ntimes;
  touri->pathnodes=NULL;
  touri->keyframe_times=NULL;
  touri->keyframe_list=NULL;
  touri->timeslist=NULL;
  touri->path_times=NULL;

  touri->global_dist=0.0;
  touri->local_dist=0.0;
  touri->startup=0;
  touri->isDefault=0;
}

/* ------------------ updateviewtour ------------------------ */

void updateviewtour(void){
  int i;
  tourdata *touri;

  viewalltours=1;
  viewanytours=0;
  for(i=0;i<ntours;i++){
    touri = tourinfo + i;
    if(touri->display==0){
      viewalltours=0;
    }
    else{
      viewanytours++;
    }
  }
}

/* ------------------ updatetourmenulabels ------------------------ */

void updatetourmenulabels(void){
  int i;
  tourdata *touri;

  if(ntours>0){
    for(i=0;i<ntours;i++){
      touri = tourinfo + i;
      STRCPY(touri->menulabel,touri->label);
    } 
  }
  updatemenu=1;
}

/* ------------------ drawtourpaths ------------------------ */

void drawtours(void){
  int i,j;
  tourdata *touri;
  pathdata *pj;
  keyframe *framej;
  int iframe_local;
  float *eye;

  float *tmp_tourcol_text,*tmp_tourcol_pathline,*tmp_tourcol_pathknots;

  tmp_tourcol_text=tourcol_text;
  if(tourcol_text[0]<0.0)tmp_tourcol_text=foregroundcolor;
 
  tmp_tourcol_pathline=tourcol_pathline;
  if(tourcol_pathline[0]<0.0)tmp_tourcol_pathline=foregroundcolor;

  tmp_tourcol_pathknots=tourcol_pathknots;
  if(tourcol_pathknots[0]<0.0)tmp_tourcol_pathknots=foregroundcolor;



  if(selected_frame!=NULL){
    selectedtour_index = selected_tour-tourinfo;
  }


  if(edittour==1){

    /* path line (non-selected)*/

    if(showtours_whenediting==1){
      antialias(1);
      glColor3fv(tmp_tourcol_pathline);
      glBegin(GL_LINES);
      for(i=0;i<ntours;i++){
        touri = tourinfo + i;
        if(touri->display==0||touri->nkeyframes<=1||selectedtour_index==i)continue;

        for(j=0;j<view_ntimes-1;j++){
          pj = touri->pathnodes + j;
          eye = pj->eye;
          glVertex3fv(eye);

          pj++;
          eye = pj->eye;
          glVertex3fv(eye);
        }
      } 
     glEnd();
     antialias(0);
    }

    /* path line (selected)*/

    touri = tourinfo + selectedtour_index;
    if(selectedtour_index!=-1&&touri->display==1&&touri->nkeyframes>1){
      glColor3fv(tourcol_selectedpathline);
      antialias(1);
      glBegin(GL_LINES);
      touri = tourinfo + selectedtour_index;

      for(j=0;j<view_ntimes-1;j++){
        pj = touri->pathnodes + j;
        eye = pj->eye;
        glVertex3fv(eye);

        pj++;
        eye = pj->eye;
        glVertex3fv(eye);
      }
      glEnd();
      antialias(0);

    }

    /* path knots */

    if(show_path_knots==1&&touri->nkeyframes>1){
      glColor3f(0.0f,0.0f,1.0f);
      glPointSize(5.0f);
      glBegin(GL_POINTS);
      for(i=0;i<ntours;i++){
        touri = tourinfo + i;
        if(touri->display==0||selectedtour_index!=i)continue;
   
        for(j=0;j<view_ntimes;j++){
          pj = touri->pathnodes + j;
          eye = pj->eye;
          glVertex3fv(eye);
        }
      }
      glEnd();

    }

    /* selected path - non-selected keyframe knots */

    touri = tourinfo + selectedtour_index;
    if(selectedtour_index!=-1&&selected_tour->display==1){
      glColor3fv(tourcol_selectedpathlineknots);
      glPointSize(10.0f);
      glBegin(GL_POINTS);
      for(j=0;j<touri->nkeyframes;j++){
        framej = touri->keyframe_list[j];
        if(framej->selected==1)continue;
        eye = framej->nodeval.eye;
        glVertex3fv(eye);
      }
      glEnd();
    }

    /* non-selcted path, keyframe knots */
    if(showtours_whenediting==1){
      glColor3fv(tmp_tourcol_pathknots);
      glPointSize(10.0f);
      glBegin(GL_POINTS);
      for(i=0;i<ntours;i++){
        touri = tourinfo + i;
        if(touri->display==0||i==selectedtour_index)continue;

        for(j=0;j<touri->nkeyframes;j++){
          framej = touri->keyframe_list[j];
          if(framej->selected==1)continue;
          eye = framej->nodeval.eye;
          glVertex3fv(eye);
        }
      }
      glEnd();
    }

    /* selected path, selected keyframe */

    if(selected_frame!=NULL&&selected_tour->display==1){
      glBegin(GL_POINTS);
      glColor3fv(tourcol_selectedknot);
      eye = selected_frame->nodeval.eye;
      glVertex3fv(eye);
      if(selected_frame->viewtype==1){
        glColor3fv(tourcol_selectedview);
        glVertex3fv(selected_frame->nodeval.xyz_view);
      }
      glEnd();

    }

    /* keyframe times */
    SNIFF_ERRORS("after select path, selected keyframe");
    CheckMemory;

    if(fontindex==SCALED_FONT)scale_3dfont();
    for(i=0;i<ntours;i++){
      touri = tourinfo + i;
      if(touri->display==0)continue;
      if(showtours_whenediting==0&&selectedtour_index!=i)continue;

      for(j=0;j<touri->nkeyframes;j++){
        framej = touri->keyframe_list[j];
        eye = framej->nodeval.eye;
        {
          char label[128];
          sprintf(label,"%8.2f",framej->disp_time+0.005);
          trimzeros(label);
          output3Text(tmp_tourcol_text,eye[0]+0.02f,eye[1]+0.015f,eye[2]+0.015f,label);
        }
      }
    }
  }

    /* keyframe avatar */

  if(show_tourlocus==1){
    switch (tourlocus_type){
      case 0:
        glColor3fv(tourcol_avatar);
        antialias(1);
        glBegin(GL_LINES);
        for(i=0;i<ntours;i++){
          float *tour_view;

          touri = tourinfo + i;
          if(touri->display==0||touri->nkeyframes<=1)continue;
          if(touri->timeslist==NULL)continue;

          iframe_local = touri->timeslist[itimes];
          pj = touri->pathnodes + iframe_local;
          if(keyframe_snap==1)pj = pj->keysnap;

          eye = pj->eye;
          tour_view = pj->tour_view;

          glVertex3fv(eye);
          glVertex3f(eye[0],eye[1],eye[2]+0.1f);
          glVertex3fv(eye);
          glVertex3f(tour_view[0],tour_view[1],tour_view[2]);
        }
        glEnd();
        antialias(0);
        break;
      case 1:
        for(i=0;i<ntours;i++){
          touri = tourinfo + i;
          if(touri->display==0||touri->nkeyframes<=1)continue;
          if(touri->timeslist==NULL)continue;

          iframe_local = touri->timeslist[itimes];
          pj = touri->pathnodes + iframe_local;
          if(keyframe_snap==1)pj = pj->keysnap;

          eye = pj->eye;

          drawcir(eye,tourrad_avatar,tourcol_avatar);

        }
        break;
      case 2:
        for(i=0;i<ntours;i++){
          float *tour_view, az_angle, dxy[2];

          touri = tourinfo + i;
          if(touri->display==0||touri->nkeyframes<=1)continue;
          if(touri->timeslist==NULL)continue;

          iframe_local = touri->timeslist[itimes];
          pj = touri->pathnodes + iframe_local;
          if(keyframe_snap==1)pj = pj->keysnap;
          eye = pj->eye;
          tour_view = pj->tour_view;
          dxy[0]=tour_view[0]-eye[0];
          dxy[1]=tour_view[1]-eye[1];
          if(dxy[0]!=0.0||dxy[1]!=0.0){
            az_angle=atan2(dxy[1],dxy[0])*RAD2DEG;
          }
          else{
            az_angle=0.0;
          }

          glPushMatrix();
          glTranslatef(eye[0],eye[1],eye[2]);
          glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));

          glRotatef(az_angle,0.0,0.0,1.0);

          draw_SVOBJECT(avatar_types[touri->glui_avatar_index],0,NULL,0);
          glPopMatrix();
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
  }

}

void drawcir(float *center, float rad, float *color){
  glColor3fv(color);
  glBegin(GL_QUADS);
  glVertex3f(center[0]-rad/2.0,center[1]-rad/2.0,center[2]);
  glVertex3f(center[0]+rad/2.0,center[1]-rad/2.0,center[2]);
  glVertex3f(center[0]+rad/2.0,center[1]+rad/2.0,center[2]);
  glVertex3f(center[0]-rad/2.0,center[1]+rad/2.0,center[2]);
  glEnd();
}


/* ------------------ drawselect_tours ------------------------ */

void drawselect_tours(void){
  int i,j;
  tourdata *touri;
  keyframe *framej;
  int color_index=0;
  unsigned char r, g, b;
  float *eye;

  glPointSize(10.0f);
  glBegin(GL_POINTS);
  for(i=0;i<ntours;i++){
    touri = tourinfo + i;
    for(j=0;j<touri->nkeyframes;j++){
  
      if(showtours_whenediting==1||selectedtour_index==i){

        if(touri->display==1){
          framej = touri->keyframe_list[j];
          eye = framej->nodeval.eye;

          getrgb(color_index+1,&r,&g,&b);
          glColor3ub(r,g,b);
          glVertex3fv(eye);
        }
      }
      color_index++;
    }
  }
  glEnd();
}


/* ------------------ createtourpaths ------------------------ */

void createtourpaths(void){
  int i,j,jj;
  int ntotal,ntotal2;
  keyframe *kf1, *kf2;
  float *lasteye, *thiseye, *nexteye;
  float *lastview, *thisview, *nextview;
  float vtime;
  int iframe_local;
  float f1, f2, dt;
  float dx, dy, denom, dx2, dy2;
  float a, b, c;
  float s1, s2, d1, d2;
  float az;
  float dist_total;
  float factor,total_time;
  float *tour_dist3a;
  float avgsum;
  int iframe_old, iframe_new;


  float vdist,vtime2,vdt;
  float tour_tstart, tour_tstop;
  float viewx_local, viewy_local, dummy;

  float dz, dist;

  keyframe **tourknotskeylist_copy;
  tourdata **tourknotstourlist_copy;
  int nframes;

  // construct keyframe list for selecting keyframes

  if(ntours==0)return;

  for(i=0;i<ntours;i++){
    tourdata *touri;

    touri = tourinfo + i;
    adjusttourtimes(touri);
  }

  ntourknots=0;
  for(i=0;i<ntours;i++){
    tourdata *touri;
    keyframe *keyj;

    touri = tourinfo + i;
    nframes=0;
    for(keyj=(touri->first_frame).next;keyj->next!=NULL;keyj=keyj->next){
      keyj->nodeval.time=keyj->noncon_time;
      keyj->disp_time=keyj->noncon_time;
      ntourknots++;
      nframes++;
    }
    touri->nkeyframes=nframes;
  }

  FREEMEMORY(tourknotskeylist);
  FREEMEMORY(tourknotstourlist);
  NewMemory((void **)&(tourknotskeylist),ntourknots*sizeof(keyframe*));
  NewMemory((void **)&(tourknotstourlist),ntourknots*sizeof(keyframe*));

  // construct spline slopes for each keyframe interval and each quantity being interpolated

  tourknotskeylist_copy=tourknotskeylist;
  tourknotstourlist_copy=tourknotstourlist;
  for(i=0;i<ntours;i++){
    tourdata *touri;
    keyframe *keyj;

    touri = tourinfo + i;
    if(viewalltours==1)touri->display=1;
    FREEMEMORY(touri->keyframe_list);
    NewMemory((void **)&(touri->keyframe_list),touri->nkeyframes*sizeof(keyframe*));
    FREEMEMORY(touri->keyframe_times);
    NewMemory((void **)&(touri->keyframe_times),touri->nkeyframes*sizeof(float));
    for(keyj=(touri->first_frame).next,j=0;keyj->next!=NULL;keyj=keyj->next,j++){
      touri->keyframe_list[j]=keyj;
      touri->keyframe_times[j]=keyj->nodeval.time;
    }
    if(touri->nkeyframes<=1)continue;
    for(keyj=(touri->first_frame).next,j=0;keyj->next!=NULL;keyj=keyj->next,j++){
      keyframe *lastkey, *thiskey, *nextkey;

      *tourknotskeylist_copy++ = keyj;
      *tourknotstourlist_copy++ = touri;
      
      keyj->selected=0;
      keyj->dist=0.0;
      keyj->npoints=0;

      lastkey=keyj->prev;
      thiskey=keyj;
      nextkey=keyj->next;

      if(lastkey==NULL)lastkey=thiskey;
      if(nextkey==NULL)nextkey=thiskey;
       

      if(j==0&&touri->periodic==1){
        lastkey=touri->keyframe_list[touri->nkeyframes-2];
      }
      if(j==touri->nkeyframes-1&&touri->periodic==1){
        nextkey=touri->keyframe_list[1];
      }
      if(touri->global_tension_flag==1){
        a=touri->global_tension;
      }
      else{
        a=keyj->tension;
      }
      b=keyj->bias;
      c=keyj->continuity;

      s1=(1.0-a)*(1.0+b)*(1-c)/2.0;
      s2=(1.0-a)*(1.0-b)*(1+c)/2.0;
      d1=(1.0-a)*(1.0+b)*(1+c)/2.0;
      d2=(1.0-a)*(1.0-b)*(1-c)/2.0;

      lasteye = lastkey->nodeval.eye;
      thiseye = thiskey->nodeval.eye;
      nexteye = nextkey->nodeval.eye;

      lastview = lastkey->nodeval.xyz_view;
      thisview = thiskey->nodeval.xyz_view;
      nextview = nextkey->nodeval.xyz_view;

      if(touri->periodic==0&&j==0){
        keyj->s_eye[0]=0.0;
        keyj->s_eye[1]=0.0;
        keyj->s_eye[2]=0.0;

        keyj->s_eye[3]=0.0;
        keyj->s_eye[4]=0.0;

        keyj->s_eye[5]=0.0;

        keyj->d_eye[0]=nexteye[0] - thiseye[0];
        keyj->d_eye[1]=nexteye[1] - thiseye[1];
        keyj->d_eye[2]=nexteye[2] - thiseye[2];

        keyj->d_eye[3]=nextkey->az_path           - thiskey->az_path;
        keyj->d_eye[4]=nextkey->nodeval.zoom      - thiskey->nodeval.zoom;
        keyj->d_eye[5]=nextkey->nodeval.elev_path - thiskey->nodeval.elev_path;

        keyj->s_xyz_view[0]=0.0;
        keyj->s_xyz_view[1]=0.0;
        keyj->s_xyz_view[2]=0.0;

        keyj->d_xyz_view[0]=nextview[0] - thisview[0];
        keyj->d_xyz_view[1]=nextview[1] - thisview[1];
        keyj->d_xyz_view[2]=nextview[2] - thisview[2];
      }
      else if(touri->periodic==0&&j==touri->nkeyframes-1){
        keyj->s_eye[0]=thiseye[0] - lasteye[0];
        keyj->s_eye[1]=thiseye[1] - lasteye[1];
        keyj->s_eye[2]=thiseye[2] - lasteye[2];

        keyj->s_eye[3]=thiskey->az_path           - lastkey->az_path;
        keyj->s_eye[4]=thiskey->nodeval.zoom      - lastkey->nodeval.zoom;
        keyj->s_eye[5]=thiskey->nodeval.elev_path - lastkey->nodeval.elev_path;

        keyj->d_eye[0]=0.0;
        keyj->d_eye[1]=0.0;
        keyj->d_eye[2]=0.0;

        keyj->d_eye[3]=0.0;
        keyj->d_eye[4]=0.0;
        keyj->d_eye[5]=0.0;

        keyj->s_xyz_view[0]=thisview[0] - lastview[0];
        keyj->s_xyz_view[1]=thisview[1] - lastview[1];
        keyj->s_xyz_view[2]=thisview[2] - lastview[2];

        keyj->d_xyz_view[0]=0.0;
        keyj->d_xyz_view[1]=0.0;
        keyj->d_xyz_view[2]=0.0;
      }
      else{

#define HERM1(s1,s2,vec0,vec1,vec2,vecout)\
  vecout[0]=(s1*(vec1[0] - vec0[0]) + s2*(vec2[0]-vec1[0]));\
  vecout[1]=(s1*(vec1[1] - vec0[1]) + s2*(vec2[1]-vec1[1]));\
  vecout[2]=(s1*(vec1[2] - vec0[2]) + s2*(vec2[2]-vec1[2]))
#define HERM2(s1,s2,vecout)\
  vecout[3]=(s1*(thiskey->az_path -           lastkey->az_path) +           s2*(nextkey->az_path-          thiskey->az_path));\
  vecout[4]=(s1*(thiskey->nodeval.zoom -      lastkey->nodeval.zoom) +      s2*(nextkey->nodeval.zoom-     thiskey->nodeval.zoom));\
  vecout[5]=(s1*(thiskey->nodeval.elev_path - lastkey->nodeval.elev_path) + s2*(nextkey->nodeval.elev_path-thiskey->nodeval.elev_path))

        HERM1(s1,s2,lasteye,thiseye,nexteye,keyj->s_eye);
        HERM2(s1,s2,keyj->s_eye);
        HERM1(s1,s2,lastview,thisview,nextview,keyj->s_xyz_view);

        HERM1(d1,d2,lasteye,thiseye,nexteye,keyj->d_eye);
        HERM2(d1,d2,keyj->d_eye);
        HERM1(d1,d2,lastview,thisview,nextview,keyj->d_xyz_view);
      }
    }

    for(keyj=(touri->first_frame).next;keyj->next!=NULL;keyj=keyj->next){
      keyj->keyview_x=keyj->d_eye[0];
      keyj->keyview_y=keyj->d_eye[1];
      if(keyj->viewtype==1)adjustviewangle(keyj,&dummy,&dummy);
    }

    // evaluate quantities along path - determine distances
    // define tour_t and tour_dist (tour_t is uniform )

    iframe_local=0;
    tour_dist[0]=0.0;
    for(j=0;j<view_ntimes;j++){
      float *xyz_view, *eye, *tour_view;
      pathdata *pj,*pjm1,*pjp1;

      pj = touri->pathnodes + j;
      if(j+1<view_ntimes){
        pjp1 = touri->pathnodes + j + 1;
      }
      else{
        pjp1 = pj;
      }
      f1 = (view_ntimes-1-j)/(float)(view_ntimes-1);
      f2 = 1-f1;
      vtime = view_tstart*f1 + view_tstop*f2;
      iframe_local = isearch(touri->keyframe_times,touri->nkeyframes,vtime,iframe_local);
      kf1 = touri->keyframe_list[iframe_local];
      kf2 = touri->keyframe_list[iframe_local+1];

      pj->keysnap=&kf1->nodeval;
      dt = kf2->nodeval.time - kf1->nodeval.time;
      f1 = CLAMP((vtime - kf1->nodeval.time)/dt,0.0,1.0);
      f2 = 1 - f1;
      pj->time=vtime;
      touri->path_times[j]=vtime;

      eye=pj->eye;
      xyz_view=pj->xyz_view;
      tour_view=pj->tour_view;

      if(kf1->nodeval.eye[0]==kf2->nodeval.eye[0]&&
         kf1->nodeval.eye[1]==kf2->nodeval.eye[1]&&
         kf1->nodeval.eye[2]==kf2->nodeval.eye[2]){
        eye[0] = hermiteeye(1.0,0,kf1->prev,kf1,&viewx_local);
        eye[1] = hermiteeye(1.0,1,kf1->prev,kf1,&viewy_local);
      }
      else{
        eye[0] = hermiteeye(f1,0,kf1,kf2,&viewx_local);
        eye[1] = hermiteeye(f1,1,kf1,kf2,&viewy_local);
      }
      eye[2] = hermiteeye(f1,2,kf1,kf2,NULL);
      eye[3] = hermiteeye(f1,3,kf1,kf2,NULL);
      //slerp(pj->xyz_view,pjp1->xyz_view,f1,xyz_view);
      xyz_view[0] = hermiteview(f1,0,kf1,kf2,NULL);
      xyz_view[1] = hermiteview(f1,1,kf1,kf2,NULL);
      xyz_view[2] = hermiteview(f1,2,kf1,kf2,NULL);

      pj->elev_path = hermiteeye(f1,5,kf1,kf2,NULL);
      pj->zoom   = hermiteeye(f1,4,kf1,kf2,NULL);

      tour_view[0]=viewx_local;
      tour_view[1]=viewy_local;
      tour_view[2]=0.0;
      tour_t[j]=vtime;
      if(j!=0){
        pjm1 = pj - 1;
        dx = eye[0]-pjm1->eye[0];
        dy = eye[1]-pjm1->eye[1];
        dz = eye[2]-pjm1->eye[2];
        dist = sqrt(dx*dx+dy*dy+dz*dz);
        tour_dist[j]=tour_dist[j-1] + dist;
        kf1->dist += dist;
      }
    }

    // construct running "cumdist" info

    kf1 = touri->first_frame.next;
    kf1->cumdist=0.0;
    dist_total=0.0;
    for(keyj=kf1->next;keyj->next!=NULL;keyj=keyj->next){
      keyj->cumdist = keyj->prev->cumdist + keyj->prev->dist;
      dist_total += keyj->prev->dist;
    }

    // decide fraction of global path vs. local path

    touri->global_dist=0.0;
    touri->local_dist=0.0;
    factor=0.0;
    total_time=0.0;
    for(keyj=(touri->first_frame).next;keyj->next!=NULL;keyj=keyj->next){
      if(keyj->next->next!=NULL)total_time += keyj->next->noncon_time - keyj->noncon_time;
      if(tour_constant_vel==1){
        touri->global_dist+=keyj->dist;
      }
      else{
        if(keyj->next->next!=NULL)touri->local_dist+=keyj->next->noncon_time-keyj->noncon_time;
      }
    }
    if(total_time!=0.0){
      factor = touri->local_dist/total_time;
    }

    // find number of points for each interval

    ntotal=0;
    for(keyj=(touri->first_frame).next;keyj->next!=NULL;keyj=keyj->next){
      if(tour_constant_vel==1){
        keyj->npoints=view_ntimes*(1.0-factor)*keyj->dist/touri->global_dist;
      }
      else{
        keyj->npoints=0;
        if(keyj->next->next!=NULL)keyj->npoints=view_ntimes*factor*(float)(keyj->next->noncon_time-keyj->noncon_time)/(float)touri->local_dist;
      }
      ntotal += keyj->npoints;
    }
    touri->first_frame.next->npoints += view_ntimes - ntotal;

    if(tour_constant_vel==1){
      ntotal2=0;
      touri->first_frame.next->disp_time=view_tstart;
      {
        float vtime_temp;

        for(keyj=(touri->first_frame).next;keyj->next->next!=NULL;keyj=keyj->next){
          ntotal2+=keyj->npoints;
          vtime_temp = view_tstart + (float)ntotal2/(float)view_ntimes*(view_tstop-view_tstart);
          keyj->next->disp_time=vtime_temp;
        }
      }
    }

   // construct distance array based on number of points in each interval

    jj = 0;
    tour_dist3a = tour_dist3 + 5;
    for(keyj=(touri->first_frame).next;keyj->next!=NULL;keyj=keyj->next){
      for(j=0;j<keyj->npoints;j++){
        tour_dist3a[jj] = keyj->cumdist + keyj->dist*(float)j/(float)keyj->npoints;
        jj++;
      }
    }

    // average tour_dist3 array and copy into tour_dist2

    for(j=0;j<5;j++){
      tour_dist3a[-1-j]=tour_dist3a[-j]-(tour_dist3a[1]-tour_dist3a[0]);
      tour_dist3a[view_ntimes+j]=tour_dist3a[view_ntimes-1+j]+(tour_dist3a[view_ntimes-1]-tour_dist3a[view_ntimes-2]);
    }

    for(j=0;j<view_ntimes;j++){
      avgsum=0.0;
      for(jj=-5;jj<6;jj++){
        avgsum += tour_dist3a[j+jj];
      }
      tour_dist2[j]=avgsum/11.0;
    }

    iframe_local = 0;
    tour_t2[0]=0.0;
    tour_dist2[0]=0.0;
    tour_tstart = touri->keyframe_list[0]->nodeval.time;
    tour_tstop = touri->keyframe_list[touri->nkeyframes-1]->nodeval.time;
    vdt = (tour_tstop - tour_tstart)/(float)(view_ntimes-1);
    for(j=1;j<view_ntimes;j++){
      vdist = tour_dist2[j];
      iframe_local = isearch(tour_dist,view_ntimes,vdist,iframe_local);
      f1 = (vdist-tour_dist[iframe_local])/(tour_dist[iframe_local+1]-tour_dist[iframe_local]);
      f2 = 1 - f1;
      tour_t2[j] = f2*tour_t[iframe_local] + f1*tour_t[iframe_local+1] ;
    }
    iframe_old=-1;
    for(j=0;j<view_ntimes;j++){
      float *xyz_view, *eye, *tour_view;
      pathdata *pj;

      pj = touri->pathnodes + j;
      vtime = tour_t2[j];
      vtime2 = touri->keyframe_list[0]->nodeval.time + j*vdt;
      iframe_new = isearch(touri->keyframe_times,touri->nkeyframes,vtime,iframe_old);
      kf1 = touri->keyframe_list[iframe_new];
      kf2 = touri->keyframe_list[iframe_new+1];
      dt = kf2->nodeval.time - kf1->nodeval.time;
      f1 = CLAMP((vtime - kf1->nodeval.time)/dt,0.0,1.0);
      f2 = 1 - f1;
      pj->time=vtime2;
      touri->path_times[j]=vtime2;

      eye=pj->eye;
      xyz_view=pj->xyz_view;
      tour_view=pj->tour_view;

      eye[0] = hermiteeye(f1,0,kf1,kf2,&viewx_local);
      eye[1] = hermiteeye(f1,1,kf1,kf2,&viewy_local);
      eye[2] = hermiteeye(f1,2,kf1,kf2,NULL);
      eye[3] = hermiteeye(f1,3,kf1,kf2,NULL);
      pj->zoom   = hermiteeye(f1,4,kf1,kf2,NULL);
      pj->elev_path = hermiteeye(f1,5,kf1,kf2,NULL);

      xyz_view[0] = hermiteview(f1,0,kf1,kf2,NULL);
      xyz_view[1] = hermiteview(f1,1,kf1,kf2,NULL);
      xyz_view[2] = hermiteview(f1,2,kf1,kf2,NULL);

      if(kf1->viewtype==0||kf2->viewtype==0){
        dx = viewx_local;
        dy = viewy_local;
        denom = 10.0*sqrt(dx*dx+dy*dy);
        if(denom==0.0)denom=1.0;
        dx /= denom;
        dy /= denom;
        az = pj->eye[3];
        az = az*DEG2RAD;
        dx2 = dx*cos(az) - dy*sin(az);
        dy2 = dx*sin(az) + dy*cos(az);
        dz = tan(pj->elev_path*DEG2RAD)/10.0;
        tour_view[0]=eye[0]+dx2;
        tour_view[1]=eye[1]+dy2;
        tour_view[2]=eye[2]+dz;
      }
      else{
        dx = xyz_view[0]-eye[0];
        dy = xyz_view[1]-eye[1];
        dz = xyz_view[2]-eye[2];
        denom = 10.0*sqrt(dx*dx+dy*dy+dz*dz);
        dx /= denom;
        dy /= denom;
        dz /= denom;
        tour_view[0]=eye[0]+dx;
        tour_view[1]=eye[1]+dy;
        tour_view[2]=eye[2]+dz;
      }
      if(iframe_old!=iframe_new){
        iframe_old=iframe_new;
        pj->keysnap->tour_view[0]=tour_view[0];
        pj->keysnap->tour_view[1]=tour_view[1];
        pj->keysnap->tour_view[2]=tour_view[2];
      }
    }
    for(keyj=kf1->next;keyj->next!=NULL;keyj=keyj->next){
      keyj->nodeval.time = tour_tstart + (tour_tstop-tour_tstart)*keyj->cumdist/dist_total;
    }
    if(selected_frame!=NULL)selected_frame->selected=1;
  }
}

/* ------------------ hermiteye ------------------------ */

float hermiteeye(float t, int i, keyframe *kf1, keyframe *kf2, float *slope){
  float p0, p1, m0, m1, val;
  float t3, t2;

  switch (i) {
  case 5:
    p0=kf1->nodeval.elev_path;
    p1=kf2->nodeval.elev_path;
    break;
  case 3:
    p0 = kf1->az_path;
    p1 = kf2->az_path;
    break;
  case 4:
    p0 = kf1->nodeval.zoom;
    p1 = kf2->nodeval.zoom;
    break;

  default:
    p0 = kf1->nodeval.eye[i];
    p1 = kf2->nodeval.eye[i];
  }
  m0 = kf1->d_eye[i];
  m1 = kf2->s_eye[i];
  t2 = t*t;
  t3 = t2*t;
  val = (2*t3-3*t2+1.0)*p0 + (t3-2*t2+t)*m0 + (t3-t2)*m1 + (-2*t3+3*t2)*p1;
  if(slope!=NULL)*slope = (6*t2-6*t)*p0 + (3*t2-4*t+1.0)*m0 + (3*t2-2.0*t)*m1 + (-6*t2+6*t)*p1;
  return val;

}

/* ------------------ hermiteye ------------------------ */

float hermiteview(float t, int i, keyframe *kf1, keyframe *kf2, float *slope){
  float p0, p1, m0, m1, val;
  float t3, t2;

  p0 = kf1->nodeval.xyz_view[i];
  p1 = kf2->nodeval.xyz_view[i];
  m0 = kf1->d_xyz_view[i];
  m1 = kf2->s_xyz_view[i];
  t2 = t*t;
  t3 = t2*t;
  val = (2*t3-3*t2+1.0)*p0 + (t3-2*t2+t)*m0 + (t3-t2)*m1 + (-2*t3+3*t2)*p1;
  if(slope!=NULL)*slope = (6*t2-6*t)*p0 + (3*t2-4*t+1.0)*m0 + (3*t2-2.0*t)*m1 + (-6*t2+6*t)*p1;
  return val;

}

/* ------------------ defaulttour ------------------------ */

void defaulttour(void){
  float *eye_xyz,*az_elev;

  touring=1;
  angle_global=0;
  eye_xyz = camera_current->eye;
  az_elev = camera_current->az_elev;

  anglexy0 = az_elev[0];
  azimuth0 = camera_current->azimuth;
  eyex0 = eye_xyz[0];
  eyey0 = eye_xyz[1];
  eyez0 = eye_xyz[2];
  selected_tour=NULL;
  selected_frame=NULL;
  selectedtour_index=-1;
  selectedtour_index_old=-1;

}

/* ------------------ addframe ------------------------ */
keyframe *add_frame(keyframe *framei, float time_local, float *eye, float key_az_path, float elev_path, float bank, float params[3],
                    int viewtype,float zoom_local,float view[3]){
  keyframe *frame,*framen;
  float *feye, *fxyz_view;

  NewMemory((void **)&frame,sizeof(keyframe));
  feye = frame->nodeval.eye;
  fxyz_view = frame->nodeval.xyz_view;
  if(viewtype!=0)viewtype=1;

  framen=framei->next;
  if(framen==NULL){
    return NULL;
  }

  framei->next=frame;
  frame->next=framen;

  framen->prev=frame;
  frame->prev=framei;

  frame->az_path=key_az_path;
  frame->nodeval.elev_path=elev_path;
  frame->bank=bank;
  NORMALIZE_XYZ(feye,eye);
  feye[3]=key_az_path;
  NORMALIZE_XYZ(fxyz_view,view);

  frame->noncon_time=time_local;
  frame->disp_time=time_local;

  frame->bias=params[1];
  frame->continuity=params[2];
  frame->bias=0.0;               // no longer using bias
  frame->continuity=0.0;         // no longer using continuity
  frame->tension=params[0];
  frame->viewtype=viewtype;
  frame->nodeval.zoom=zoom_local;
  frame->keyview_x=0.0;
  frame->keyview_y=0.0;

  CheckMemory;
  return frame;
}

/* ------------------ deleteframe ------------------------ */

keyframe *delete_frame(keyframe *frame){
  keyframe *prev, *next;
  tourdata *thistour;

  thistour=selected_tour;

  prev=frame->prev;
  next=frame->next;
  prev->next=next;
  next->prev=prev;
  FREEMEMORY(frame);
  if(thistour!=NULL){
    thistour->nkeyframes--;
    if(prev->prev!=NULL)return prev;
    if(next->next!=NULL)return next;
  }
  return NULL;
}

/* ------------------ new_select ------------------------ */

void new_select(keyframe *newselect){
  if(newselect!=selected_frame&&selected_frame!=NULL)selected_frame->selected=0;
  selected_frame=newselect;
  if(newselect!=NULL)selected_frame->selected=1;
}

/* ------------------ init_circulartour ------------------------ */

void init_circulartour(void){
  int nkeyframes,j;
  float key_az_path, elev_path, key_bank, params[3],key_view[3], key_xyz[3], zoom_local;
  int viewtype=0;
  tourdata *touri;
  float key_time;
  float angle_local;
  float f1;
  float rad, cosangle, sinangle;
  float max_xyz, dx, dy, dz;
  keyframe *thisframe,*addedframe;

  touri = tourinfo;
  inittour(touri);
  touri->isDefault=1;
  touri->startup=1;
  touri->periodic=1;
  strcpy(touri->label,"Circular");
  nkeyframes=16;
  NewMemory((void **)&touri->keyframe_times, nkeyframes*sizeof(float));
  NewMemory((void **)&touri->pathnodes,view_ntimes*sizeof(pathdata));
  NewMemory((void **)&touri->path_times,view_ntimes*sizeof(float));
  key_view[0]=(xbar0+xbarORIG)/2.0;
  key_view[1]=(ybar0+ybarORIG)/2.0;
  key_view[2]=(zbar0+zbarORIG)/2.0;
  dx = ABS(xbarORIG - xbar0)/2.0;
  dy = ABS(ybarORIG - ybar0)/2.0;
  dz = ABS(zbarORIG-zbar0)/2.0;
  max_xyz=dx;
  if(dy>max_xyz)max_xyz=dy;
  if(dz>max_xyz)max_xyz=dz;

  rad = max_xyz+max_xyz/tan(20.0*DEG2RAD);
  elev_path=0.0;

  thisframe=&touri->first_frame;
  for(j=0;j<nkeyframes;j++){
    key_az_path = 0.0;
    key_bank=0.0;
    params[0] = 0.0;
    params[1] = 0.0;
    params[2] = 0.0;
    angle_local = 2.0*PI*(float)j/(float)(nkeyframes-1);
    cosangle = cos(angle_local);
    sinangle = sin(angle_local);

    key_xyz[0] = key_view[0] + rad*cosangle;
    key_xyz[1] = key_view[1] + rad*sinangle;
    key_xyz[2] = key_view[2];
    f1 = (float)j/(float)(nkeyframes-1);
    key_time = view_tstart*(1.0-f1) + view_tstop*f1;

    viewtype=1;
    zoom_local=1.0;
    addedframe=add_frame(thisframe, key_time, key_xyz, 
      key_az_path, elev_path, key_bank, params, viewtype,zoom_local,key_view);
    thisframe=addedframe;
    touri->keyframe_times[j]=key_time;
  }

}

/* ------------------ add_tour  ------------------------ */

tourdata *add_tour(char *label){
  tourdata *tourtemp=NULL,*touri;
  int nkeyframes;
  float key_az_path, elev_path, key_bank, params[3],key_view[3], key_xyz[3], zoom_local;
  int viewtype=0;
  float key_time;
  int i;
  keyframe *thisframe, *addedframe;

  delete_tourlist();
  ntours++;
  NewMemory( (void **)&tourtemp, ntours*sizeof(tourdata));
  if(ntours>1)memcpy(tourtemp,tourinfo,(ntours-1)*sizeof(tourdata));
  FREEMEMORY(tourinfo);
  tourinfo=tourtemp;
  touri = tourinfo + ntours - 1;

  inittour(touri);
  if(label==NULL){
    sprintf(touri->label,"Tour %i",ntours);
  }
  else{
    strcpy(touri->label,label);
  }
  nkeyframes=2;
  NewMemory((void **)&touri->keyframe_times, nkeyframes*sizeof(float));
  NewMemory((void **)&touri->pathnodes,view_ntimes*sizeof(pathdata));
  NewMemory((void **)&touri->path_times,view_ntimes*sizeof(float));

  key_view[0]=0.0;
  key_view[1]=0.0;
  key_view[2]=0.0;

  key_az_path = 0.0;
  key_bank = 0.0;
  elev_path=0.0;
  params[0] = 0.0;
  params[1] = 0.0;
  params[2] = 0.0;
  viewtype=0;
  zoom_local=1.0;

  key_xyz[0] = xbar0 - 1.0;
  key_xyz[1] = ybar0 - 1.0;
  key_xyz[2] = (zbar0 + zbarORIG)/2.0;
  key_time = view_tstart;
  thisframe=&touri->first_frame;
  addedframe=add_frame(thisframe,key_time, key_xyz, key_az_path, elev_path, key_bank, 
    params, viewtype,zoom_local,key_view);
  touri->keyframe_times[0]=key_time;

  key_xyz[0] = xbarORIG + 1.0;
  key_xyz[1] = ybarORIG + 1.0;
  key_xyz[2] = (zbar0 + zbarORIG)/2.0;
  key_time = view_tstop;
  thisframe=addedframe;
  add_frame(thisframe,key_time, key_xyz, key_az_path, elev_path, key_bank,
    params, viewtype,zoom_local,key_view);
  touri->keyframe_times[1]=key_time;
  touri->display=1;

  for(i=0;i<ntours;i++){
    touri=tourinfo+i;
    touri->first_frame.next->prev=&touri->first_frame;
    touri->last_frame.prev->next=&touri->last_frame;
  }
  viewalltours=1;
  for(i=0;i<ntours;i++){
    touri = tourinfo + i;
    if(touri->display==0)viewalltours=0;
  }
  updatemenu=1;
  
  updatetourmenulabels();
  createtourpaths();
  Update_Times();
  create_tourlist();
  return tourinfo + ntours-1;
}


/* ------------------ delete_tour  ------------------------ */

void delete_tour(int tour_index){
  tourdata *touri,*tourtemp;
  int i;

  delete_tourlist();
  touri = tourinfo + tour_index;
  freetour(touri);
  ntours--;
  if(ntours>0){
    NewMemory( (void **)&tourtemp, ntours*sizeof(tourdata));
    if(tour_index>0)memcpy(tourtemp,tourinfo,tour_index*sizeof(tourdata));
    if(tour_index<ntours)memcpy(tourtemp+tour_index,tourinfo+tour_index+1,(ntours-tour_index)*sizeof(tourdata));
    FREEMEMORY(tourinfo);
    tourinfo=tourtemp;
    for(i=0;i<ntours;i++){
      touri=tourinfo+i;
      touri->first_frame.next->prev=&touri->first_frame;
      touri->last_frame.prev->next=&touri->last_frame;
    }
  }
  viewalltours=1;
  for(i=0;i<ntours;i++){
    touri = tourinfo + i;
    if(touri->display==0)viewalltours=0;
  }
  if(ntours==0){
    viewalltours=0;
  }
  updatemenu=1;
  createtourpaths();
  selectedtour_index=tour_index-1;
  selectedtour_index_old=selectedtour_index;
  if(ntours>0){
    if(selectedtour_index<0)selectedtour_index=0;
    selected_tour=tourinfo+selectedtour_index;
    selected_frame=selected_tour->first_frame.next;
  }
  else{
    selected_tour=NULL;
    selected_frame=NULL;
  }
  set_glui_keyframe();
  updatetourmenulabels();
  Update_Times();
  create_tourlist();

}

/* ------------------ ReallocTourMemory  ------------------------ */


void ReallocTourMemory(void){
  int i;
  tourdata *touri;

  if(view_ntimes>0){
    for(i=0;i<ntours;i++){
      touri = tourinfo + i;
      FREEMEMORY(touri->pathnodes);
      FREEMEMORY(touri->path_times);
      NewMemory((void **)&touri->pathnodes,view_ntimes*sizeof(pathdata));
      NewMemory((void **)&touri->path_times,view_ntimes*sizeof(float));
      touri->ntimes=view_ntimes;
      }
    FREEMEMORY(tour_t);
    FREEMEMORY(tour_t2);
    FREEMEMORY(tour_dist);
    FREEMEMORY(tour_dist2);
    FREEMEMORY(tour_dist3);
    NewMemory((void **)&tour_t,view_ntimes*sizeof(float));
    NewMemory((void **)&tour_t2,view_ntimes*sizeof(float));
    NewMemory((void **)&tour_dist,view_ntimes*sizeof(float));
    NewMemory((void **)&tour_dist2,view_ntimes*sizeof(float));
    NewMemory((void **)&tour_dist3,(view_ntimes+10)*sizeof(float));
  }
}

/* ------------------ setup_tour  ------------------------ */

void setup_tour(void){

  if(ntours==0){
    ReallocTourMemory();
    ntours=1;
    NewMemory( (void **)&tourinfo, ntours*sizeof(tourdata));
    init_circulartour();
    updatetourmenulabels();
    createtourpaths();
    Update_Times();
    plotstate=getplotstate(DYNAMIC_PLOTS);
    selectedtour_index=-1;
    selectedtour_index=-1;
    selected_frame=NULL;
    selected_tour=NULL;
    if(viewalltours==1)TourMenu(-3);
  }
}

/* ------------------ adjustviewangle ------------------------ */
void adjustviewangle(keyframe *kf, float *az_path, float *elev_path){
  float dx, dy, dz;
  float dx2, dy2;
  float distxy, distxy2;
  float angle_temp, angle_temp2;
  float az, elev;

  float *eye, *xyz_view;

  dx2 = kf->keyview_x;
  dy2 = kf->keyview_y;
  distxy2 = sqrt(dx2*dx2+dy2*dy2);
  if(distxy2<=0.0)return;
  
  eye = kf->nodeval.eye;
  xyz_view = kf->nodeval.xyz_view;

  dx = xyz_view[0] - eye[0];
  dy = xyz_view[1] - eye[1];
  dz = xyz_view[2] - eye[2];

  distxy = sqrt(dx*dx+dy*dy);
  if(distxy<=0.0)return;

  dx2 /= distxy2;
  dy2 /= distxy2;

  dx /= distxy;
  dy /= distxy;
  dz /= distxy;

  angle_temp = RAD2DEG*atan2(dy,dx);
  angle_temp2 = RAD2DEG*atan2(dy2,dx2);
  az = angle_temp - angle_temp2;
  if(az>180.0)az = 360.0 - az;
  if(az<-180.0)az = az + 360.0;
  elev=atan(dz)*RAD2DEG;

  kf->nodeval.elev_path=elev;
  kf->az_path=az;
  kf->az_path=az;
  *az_path=az;
  *elev_path = elev;
}

/* ------------------ adjusttourtimes ------------------------ */

void adjusttourtimes(tourdata *touri){
  if(touri->nkeyframes>1){
    keyframe *keyj;
    float tstart, tstop,dtmin;
    int small_flag;

    tstart = touri->first_frame.next->noncon_time;
    tstop = touri->last_frame.prev->noncon_time;
    dtmin = 4.0*(tstop-tstart)/(float)view_ntimes;
    small_flag=0;
    for(keyj=(touri->first_frame).next->next;keyj->next!=NULL;keyj=keyj->next){
      float dt;

      dt = keyj->noncon_time - keyj->prev->noncon_time;
      if(dt<dtmin){
        small_flag=1;
        keyj->noncon_time = keyj->prev->noncon_time + dtmin;
      }
    }
    if(small_flag==1&&tstop>view_tstop&&tstop>0.0){
      for(keyj=(touri->first_frame).next;keyj->next!=NULL;keyj=keyj->next){
        keyj->noncon_time = keyj->noncon_time*view_tstop/tstop;
      }
    }
  }
}

