// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "flowfiles.h"
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "MALLOC.h"
#include "smokeviewvars.h"

// svn revision character string
char IOtour_revision[]="$Revision$";

void drawcir(float *center, float rad, float *color);
//void TourMenu(int val);
float hermiteeye(float f1, int i, keyframe *kf1, keyframe *kf2, float *slope);
float hermiteview(float t, int i, keyframe *kf1, keyframe *kf2, float *slope);
void draw_SVOBJECT(sv_object *object, int iframe,propdata *prop);

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
  FREEMEMORY(touri->path_timeslist);
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
  touri->npath=view_ntimes;
  touri->pathnodes=NULL;
  touri->keyframe_times=NULL;
  touri->keyframe_list=NULL;
  touri->path_timeslist=NULL;
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
  int iframe;
  float *eye;

  float *tmp_tourcol_text,*tmp_tourcol_pathline,*tmp_tourcol_pathknots;

  tmp_tourcol_text=tourcol_text;
  if(tourcol_text[0]<0.0)tmp_tourcol_text=foregroundcolor;
 
  tmp_tourcol_pathline=tourcol_pathline;
  if(tourcol_pathline[0]<0.0)tmp_tourcol_pathline=foregroundcolor;

  tmp_tourcol_pathknots=tourcol_pathknots;
  if(tourcol_pathknots[0]<0.0)tmp_tourcol_pathknots=foregroundcolor;



  if(selected_frame!=NULL){
//    selected_tour= selected_frame->key_tour;
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
        glVertex3fv(selected_frame->nodeval.aview);
      }
      glEnd();

    }

    /* keyframe times */
    sniffErrors("after select path, selected keyframe");
    CheckMemory;

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

  //show_tourlocus=1;
  //tourlocus_type=2;
  if(show_tourlocus==1){
    switch (tourlocus_type){
      case 0:
        glColor3fv(tourcol_avatar);
        antialias(1);
        glBegin(GL_LINES);
        for(i=0;i<ntours;i++){
          float *oview;

          touri = tourinfo + i;
          if(touri->display==0||touri->nkeyframes<=1)continue;
          if(touri->path_timeslist==NULL)continue;

          iframe = touri->path_timeslist[itime];
          pj = touri->pathnodes + iframe;
          if(keyframe_snap==1)pj = pj->keysnap;

          eye = pj->eye;
          oview = pj->oview;

          glVertex3fv(eye);
          glVertex3f(eye[0],eye[1],eye[2]+0.1f);
          glVertex3fv(eye);
          glVertex3f(oview[0],oview[1],oview[2]);
        }
        glEnd();
        antialias(0);
        break;
      case 1:
        for(i=0;i<ntours;i++){
          touri = tourinfo + i;
          if(touri->display==0||touri->nkeyframes<=1)continue;
          if(touri->path_timeslist==NULL)continue;

          iframe = touri->path_timeslist[itime];
          pj = touri->pathnodes + iframe;
          if(keyframe_snap==1)pj = pj->keysnap;

          eye = pj->eye;

          drawcir(eye,tourrad_avatar,tourcol_avatar);

        }
        break;
      case 2:
        for(i=0;i<ntours;i++){
          float *oview, az_angle, dxy[2];

          touri = tourinfo + i;
          if(touri->display==0||touri->nkeyframes<=1)continue;
          if(touri->path_timeslist==NULL)continue;

          iframe = touri->path_timeslist[itime];
          pj = touri->pathnodes + iframe;
          if(keyframe_snap==1)pj = pj->keysnap;
          eye = pj->eye;
          oview = pj->oview;
          dxy[0]=oview[0]-eye[0];
          dxy[1]=oview[1]-eye[1];
          if(dxy[0]!=0.0||dxy[1]!=0.0){
            az_angle=atan2(dxy[1],dxy[0])*180.0/3.14159;
          }
          else{
            az_angle=0.0;
          }

          glPushMatrix();
          glTranslatef(eye[0],eye[1],eye[2]);
          glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);

          glRotatef(az_angle,0.0,0.0,1.0);

          draw_SVOBJECT(avatar_types[touri->glui_avatar_index],0,NULL);
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
  tourdata *touri;
  keyframe *keyj,*kf1,*kf2;
  keyframe *key0, *key1, *key2;
  pathdata *pj,*pjm1;
  float *eye0, *eye1, *eye2;
  float *aview0, *aview1, *aview2;
  float *eye, *aview, *oview;
  float vtime;
  int iframe;
  float f1, f2, dt;
  float dx, dy, denom, dx2, dy2;
  float a, b, c;
  float s1, s2, d1, d2;
  float del1, del2;
  float sfactor, dfactor, az;
  float dist_total;
  float factor,total_time;
  float *tour_dist3a;
  float avgsum;
  int iframe_old, iframe_new;


  float vdist,vtime2,vdt;
  float tour_tstart, tour_tstop;
  float viewx, viewy, dummy;

  float dz, dist, pi;
//  keyframe *framejm1;

  keyframe **tourknotskeylist_copy;
  tourdata **tourknotstourlist_copy;
  int nframes;

  // construct keyframe list for selecting keyframes

  if(ntours==0)return;

  for(i=0;i<ntours;i++){
    touri = tourinfo + i;
    adjusttourtimes(touri);
  }

  ntourknots=0;
  pi=4.0*atan(1.0);
  for(i=0;i<ntours;i++){
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
      *tourknotskeylist_copy++ = keyj;
      *tourknotstourlist_copy++ = touri;
      keyj->selected=0;
      keyj->dist=0.0;
      keyj->npoints=0;

      key0=keyj->prev;
      key1=keyj;
      key2=keyj->next;
      if(j==0&&touri->periodic==1){
        key0=touri->keyframe_list[touri->nkeyframes-2];
      }
      if(j==touri->nkeyframes-1&&touri->periodic==1){
        key2=touri->keyframe_list[1];
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

      eye0 = key0->nodeval.eye;
      eye1 = key1->nodeval.eye;
      eye2 = key2->nodeval.eye;
      aview0 = key0->nodeval.aview;
      aview1 = key1->nodeval.aview;
      aview2 = key2->nodeval.aview;


      if(touri->periodic==0&&j==0){
        keyj->s_eye[0]=0.0;
        keyj->s_eye[1]=0.0;
        keyj->s_eye[2]=0.0;

        keyj->s_eye[3]=0.0;
        keyj->s_eye[4]=0.0;

        keyj->s_eye[5]=0.0;

        keyj->d_eye[0]=eye2[0] - eye1[0];
        keyj->d_eye[1]=eye2[1] - eye1[1];
        keyj->d_eye[2]=eye2[2] - eye1[2];

        keyj->d_eye[3]=key2->az_path - key1->az_path;
        keyj->d_eye[4]=key2->nodeval.zoom - key1->nodeval.zoom;

        keyj->d_eye[5]=key2->nodeval.elev_path - key1->nodeval.elev_path;

        keyj->s_aview[0]=0.0;
        keyj->s_aview[1]=0.0;
        keyj->s_aview[2]=0.0;

        keyj->d_aview[0]=aview2[0] - aview1[0];
        keyj->d_aview[1]=aview2[1] - aview1[1];
        keyj->d_aview[2]=aview2[2] - aview1[2];
      }
      else if(touri->periodic==0&&j==touri->nkeyframes-1){
        keyj->s_eye[0]=eye1[0] - eye0[0];
        keyj->s_eye[1]=eye1[1] - eye0[1];
        keyj->s_eye[2]=eye1[2] - eye0[2];

        keyj->s_eye[3]=key1->az_path - key0->az_path;
        keyj->s_eye[4]=key1->nodeval.zoom - key0->nodeval.zoom;

        keyj->s_eye[5]=key1->nodeval.elev_path - key0->nodeval.elev_path;

        keyj->d_eye[0]=0.0;
        keyj->d_eye[1]=0.0;
        keyj->d_eye[2]=0.0;

        keyj->d_eye[3]=0.0;
        keyj->d_eye[4]=0.0;

        keyj->d_eye[5]=0.0;

        keyj->s_aview[0]=aview1[0] - aview0[0];
        keyj->s_aview[1]=aview1[1] - aview0[1];
        keyj->s_aview[2]=aview1[2] - aview0[2];

        keyj->d_aview[0]=0.0;
        keyj->d_aview[1]=0.0;
        keyj->d_aview[2]=0.0;
      }
      else{
        del1 = key1->nodeval.time - key0->nodeval.time;
        del2 = key2->nodeval.time - key1->nodeval.time;
        sfactor = 2*del2/(del1 + del2);
        dfactor = 2*del1/(del1 + del2);
        sfactor = 1.0;
        dfactor = 1.0;
        keyj->s_eye[0]=sfactor*(s1*(eye1[0] - eye0[0]) + s2*(eye2[0]-eye1[0]));
        keyj->s_eye[1]=sfactor*(s1*(eye1[1] - eye0[1]) + s2*(eye2[1]-eye1[1]));
        keyj->s_eye[2]=sfactor*(s1*(eye1[2] - eye0[2]) + s2*(eye2[2]-eye1[2]));

        keyj->s_eye[3]=sfactor*(s1*(key1->az_path -           key0->az_path) +           s2*(key2->az_path-          key1->az_path));
        keyj->s_eye[4]=sfactor*(s1*(key1->nodeval.zoom -      key0->nodeval.zoom) +      s2*(key2->nodeval.zoom-     key1->nodeval.zoom));
        keyj->s_eye[5]=sfactor*(s1*(key1->nodeval.elev_path - key0->nodeval.elev_path) + s2*(key2->nodeval.elev_path-key1->nodeval.elev_path));
        keyj->s_aview[0]=sfactor*(s1*(aview1[0] - aview0[0]) + s2*(aview2[0]-aview1[0]));
        keyj->s_aview[1]=sfactor*(s1*(aview1[1] - aview0[1]) + s2*(aview2[1]-aview1[1]));
        keyj->s_aview[2]=sfactor*(s1*(aview1[2] - aview0[2]) + s2*(aview2[2]-aview1[2]));

        keyj->d_eye[0]=dfactor*(d1*(eye1[0] - eye0[0]) + d2*(eye2[0]-eye1[0]));
        keyj->d_eye[1]=dfactor*(d1*(eye1[1] - eye0[1]) + d2*(eye2[1]-eye1[1]));
        keyj->d_eye[2]=dfactor*(d1*(eye1[2] - eye0[2]) + d2*(eye2[2]-eye1[2]));

        keyj->d_eye[3]=dfactor*(d1*(key1->az_path -           key0->az_path) +           d2*(key2->az_path-          key1->az_path));
        keyj->d_eye[4]=dfactor*(d1*(key1->nodeval.zoom -      key0->nodeval.zoom) +      d2*(key2->nodeval.zoom-     key1->nodeval.zoom));
        keyj->d_eye[5]=dfactor*(d1*(key1->nodeval.elev_path - key0->nodeval.elev_path) + d2*(key2->nodeval.elev_path-key1->nodeval.elev_path));

        keyj->d_aview[0]=dfactor*(d1*(aview1[0] - aview0[0]) + d2*(aview2[0]-aview1[0]));
        keyj->d_aview[1]=dfactor*(d1*(aview1[1] - aview0[1]) + d2*(aview2[1]-aview1[1]));
        keyj->d_aview[2]=dfactor*(d1*(aview1[2] - aview0[2]) + d2*(aview2[2]-aview1[2]));

      }
    }

    for(keyj=(touri->first_frame).next;keyj->next!=NULL;keyj=keyj->next){
      keyj->keyview_x=keyj->d_eye[0];
      keyj->keyview_y=keyj->d_eye[1];
      if(keyj->viewtype==1)adjustviewangle(keyj,&dummy,&dummy);
    }

    // evaluate quantities along path - determine distances
    // define tour_t and tour_dist (tour_t is uniform )

    iframe=0;
    tour_dist[0]=0.0;
    for(j=0;j<view_ntimes;j++){
      pj = touri->pathnodes + j;
      f1 = (view_ntimes-1-j)/(float)(view_ntimes-1);
      f2 = 1-f1;
      vtime = view_tstart*f1 + view_tstop*f2;
      iframe = interval_search(touri->keyframe_times,touri->nkeyframes,vtime,iframe);
      kf1 = touri->keyframe_list[iframe];
      kf2 = touri->keyframe_list[iframe+1];
      pj->keysnap=&kf1->nodeval;
      dt = kf2->nodeval.time - kf1->nodeval.time;
      f1 = (vtime - kf1->nodeval.time)/dt;
      if(f1<0.0)f1=0.0;
      if(f1>1.0)f1=1.0;
      f2 = 1 - f1;
      pj->time=vtime;
      touri->path_times[j]=vtime;

      eye=pj->eye;
      aview=pj->aview;
      oview=pj->oview;

      if(kf1->nodeval.eye[0]==kf2->nodeval.eye[0]&&
         kf1->nodeval.eye[1]==kf2->nodeval.eye[1]&&
         kf1->nodeval.eye[2]==kf2->nodeval.eye[2]){
        eye[0] = hermiteeye(1.0,0,kf1->prev,kf1,&viewx);
        eye[1] = hermiteeye(1.0,1,kf1->prev,kf1,&viewy);
      }
      else{
        eye[0] = hermiteeye(f1,0,kf1,kf2,&viewx);
        eye[1] = hermiteeye(f1,1,kf1,kf2,&viewy);
      }
      eye[2] = hermiteeye(f1,2,kf1,kf2,&dummy);
      eye[3] = hermiteeye(f1,3,kf1,kf2,&dummy);
      aview[0] = hermiteview(f1,0,kf1,kf2,&dummy);
      aview[1] = hermiteview(f1,1,kf1,kf2,&dummy);
      aview[2] = hermiteview(f1,2,kf1,kf2,&dummy);

      pj->elev_path = hermiteeye(f1,5,kf1,kf2,&dummy);
      pj->zoom   = hermiteeye(f1,4,kf1,kf2,&dummy);

      oview[0]=viewx;
      oview[1]=viewy;
      oview[2]=0.0;
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

    // constuct running "cumdist" info

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

    iframe = 0;
    tour_t2[0]=0.0;
    tour_dist2[0]=0.0;
    tour_tstart = touri->keyframe_list[0]->nodeval.time;
    tour_tstop = touri->keyframe_list[touri->nkeyframes-1]->nodeval.time;
    vdt = (tour_tstop - tour_tstart)/(float)(view_ntimes-1);
    for(j=1;j<view_ntimes;j++){
      vdist = tour_dist2[j];
      iframe = interval_search(tour_dist,view_ntimes,vdist,iframe);
      f1 = (vdist-tour_dist[iframe])/(tour_dist[iframe+1]-tour_dist[iframe]);
      f2 = 1 - f1;
      tour_t2[j] = f2*tour_t[iframe] + f1*tour_t[iframe+1] ;
    }
    iframe_old=-1;
    for(j=0;j<view_ntimes;j++){
      pj = touri->pathnodes + j;
      vtime = tour_t2[j];
      vtime2 = touri->keyframe_list[0]->nodeval.time + j*vdt;
      iframe_new = interval_search(touri->keyframe_times,touri->nkeyframes,vtime,iframe_old);
      kf1 = touri->keyframe_list[iframe_new];
      kf2 = touri->keyframe_list[iframe_new+1];
      dt = kf2->nodeval.time - kf1->nodeval.time;
      f1 = (vtime - kf1->nodeval.time)/dt;
      if(f1<0.0)f1=0.0;
      if(f1>1.0)f1=1.0;
      f2 = 1 - f1;
      pj->time=vtime2;
      touri->path_times[j]=vtime2;

      eye=pj->eye;
      aview=pj->aview;
      oview=pj->oview;

      eye[0] = hermiteeye(f1,0,kf1,kf2,&viewx);
      eye[1] = hermiteeye(f1,1,kf1,kf2,&viewy);
      eye[2] = hermiteeye(f1,2,kf1,kf2,&dummy);
      eye[3] = hermiteeye(f1,3,kf1,kf2,&dummy);
      pj->zoom   = hermiteeye(f1,4,kf1,kf2,&dummy);
      pj->elev_path = hermiteeye(f1,5,kf1,kf2,&dummy);

      aview[0] = hermiteview(f1,0,kf1,kf2,&dummy);
      aview[1] = hermiteview(f1,1,kf1,kf2,&dummy);
      aview[2] = hermiteview(f1,2,kf1,kf2,&dummy);

      if(kf1->viewtype==0||kf2->viewtype==0){
        dx = viewx;
        dy = viewy;
        denom = 10.0*sqrt(dx*dx+dy*dy);
        if(denom==0.0)denom=1.0;
        dx /= denom;
        dy /= denom;
        az = pj->eye[3];
        az = az*pi/180.0;
        dx2 = dx*cos(az) - dy*sin(az);
        dy2 = dx*sin(az) + dy*cos(az);
        dz = tan(pj->elev_path*pi/180.0)/10.0;
        oview[0]=eye[0]+dx2;
        oview[1]=eye[1]+dy2;
        oview[2]=eye[2]+dz;
      }
      else{
        dx = aview[0]-eye[0];
        dy = aview[1]-eye[1];
        dz = aview[2]-eye[2];
        denom = 10.0*sqrt(dx*dx+dy*dy+dz*dz);
        dx /= denom;
        dy /= denom;
        dz /= denom;
        oview[0]=eye[0]+dx;
        oview[1]=eye[1]+dy;
        oview[2]=eye[2]+dz;
      }
      if(iframe_old!=iframe_new){
        iframe_old=iframe_new;
        pj->keysnap->oview[0]=oview[0];
        pj->keysnap->oview[1]=oview[1];
        pj->keysnap->oview[2]=oview[2];
       // pj->keysnap->eye[0]=eye[0];
       // pj->keysnap->eye[1]=eye[1];
       // pj->keysnap->eye[2]=eye[2];
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
  *slope = (6*t2-6*t)*p0 + (3*t2-4*t+1.0)*m0 + (3*t2-2.0*t)*m1 + (-6*t2+6*t)*p1;
  return val;

}

/* ------------------ hermiteye ------------------------ */

float hermiteview(float t, int i, keyframe *kf1, keyframe *kf2, float *slope){
  float p0, p1, m0, m1, val;
  float t3, t2;

  p0 = kf1->nodeval.aview[i];
  p1 = kf2->nodeval.aview[i];
  m0 = kf1->d_aview[i];
  m1 = kf2->s_aview[i];
  t2 = t*t;
  t3 = t2*t;
  val = (2*t3-3*t2+1.0)*p0 + (t3-2*t2+t)*m0 + (t3-t2)*m1 + (-2*t3+3*t2)*p1;
  *slope = (6*t2-6*t)*p0 + (3*t2-4*t+1.0)*m0 + (3*t2-2.0*t)*m1 + (-6*t2+6*t)*p1;
  return val;

}

/* ------------------ defaulttour ------------------------ */

void defaulttour(void){
  float *eye_xyz,*angle_zx;

  touring=1;
  angle=0;
  eye_xyz = camera_current->eye;
  angle_zx = camera_current->angle_zx;

  anglexy0 = angle_zx[0];
  direction_angle0 = camera_current->direction_angle;
  eyex0 = eye_xyz[0];
  eyey0 = eye_xyz[1];
  eyez0 = eye_xyz[2];
  selected_tour=NULL;
  selected_frame=NULL;
  selectedtour_index=-1;
  selectedtour_index_old=-1;

}

/* ------------------ addframe ------------------------ */
keyframe *add_frame(keyframe *framei, float time, float *eye, float key_az_path, float elev_path, float bank, float params[3],
                    int viewtype,float zoom,float view[3]){
  keyframe *frame,*framen;
  float *feye, *faview;

  NewMemory((void **)&frame,sizeof(keyframe));
  feye = frame->nodeval.eye;
  faview = frame->nodeval.aview;
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
  feye[0]=(eye[0]-xbar0)/xyzmaxdiff;
  feye[1]=(eye[1]-ybar0)/xyzmaxdiff;
  feye[2]=(eye[2]-zbar0)/xyzmaxdiff;
  feye[3]=key_az_path;
  faview[0]=(view[0]-xbar0)/xyzmaxdiff;
  faview[1]=(view[1]-ybar0)/xyzmaxdiff;
  faview[2]=(view[2]-zbar0)/xyzmaxdiff;

  frame->noncon_time=time;
  frame->disp_time=time;

  frame->bias=params[1];
  frame->continuity=params[2];
  frame->bias=0.0;               // no longer using bias
  frame->continuity=0.0;         // no longer using continuity
  frame->tension=params[0];
  frame->viewtype=viewtype;
  frame->nodeval.zoom=zoom;
  frame->keyview_x=0.0;
  frame->keyview_y=0.0;

  CheckMemory;
  return frame;
}

/* ------------------ deleteframe ------------------------ */

keyframe *delete_frame(keyframe *frame){
  keyframe *prev, *next;
  tourdata *thistour;

  //thistour=frame->key_tour;
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
  float key_az_path, elev_path, key_bank, params[3],key_view[3], key_xyz[3], zoom;
  int viewtype=0;
  tourdata *touri;
  float key_time;
  float angle_local;
  float f1;
  float A, B, rad, cosangle, sinangle;
  float max_xyz, dx, dy, dz;
  float pi;
  keyframe *thisframe,*addedframe;

  pi=4.0*atan(1.0);

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
  dx = fabs(xbarORIG - xbar0)/2.0;
  dy = fabs(ybarORIG - ybar0)/2.0;
  dz = fabs(zbarORIG-zbar0)/2.0;
  max_xyz=dx;
  if(dy>max_xyz)max_xyz=dy;
  if(dz>max_xyz)max_xyz=dz;

  rad = max_xyz+max_xyz/tan(20.0/180.0*3.14159);
  elev_path=0.0;

  thisframe=&touri->first_frame;
  for(j=0;j<nkeyframes;j++){
    key_az_path = 0.0;
    key_bank=0.0;
    params[0] = 0.0;
    params[1] = 0.0;
    params[2] = 0.0;
    angle_local = 2.0*pi*(float)j/(float)(nkeyframes-1);
    cosangle = cos(angle_local);
    sinangle = sin(angle_local);

    key_xyz[0] = key_view[0] + rad*cosangle;
    key_xyz[1] = key_view[1] + rad*sinangle;
    key_xyz[2] = key_view[2];
    f1 = (float)j/(float)(nkeyframes-1);
    key_time = view_tstart*(1.0-f1) + view_tstop*f1;

    viewtype=1;
    zoom=1.0;
    addedframe=add_frame(thisframe, key_time, key_xyz, 
      key_az_path, elev_path, key_bank, params, viewtype,zoom,key_view);
    thisframe=addedframe;
    touri->keyframe_times[j]=key_time;
  }

}

/* ------------------ add_tour  ------------------------ */

tourdata *add_tour(char *label){
  tourdata *tourtemp=NULL,*touri;
  int nkeyframes;
  float key_az_path, elev_path, key_bank, params[3],key_view[3], key_xyz[3], zoom;
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
  zoom=1.0;

  key_xyz[0] = xbar0 - 1.0;
  key_xyz[1] = ybar0 - 1.0;
  key_xyz[2] = (zbar0 + zbarORIG)/2.0;
  key_time = view_tstart;
  thisframe=&touri->first_frame;
  addedframe=add_frame(thisframe,key_time, key_xyz, key_az_path, elev_path, key_bank, 
    params, viewtype,zoom,key_view);
  touri->keyframe_times[0]=key_time;

  key_xyz[0] = xbarORIG + 1.0;
  key_xyz[1] = ybarORIG + 1.0;
  key_xyz[2] = (zbar0 + zbarORIG)/2.0;
  key_time = view_tstop;
  thisframe=addedframe;
  add_frame(thisframe,key_time, key_xyz, key_az_path, elev_path, key_bank,
    params, viewtype,zoom,key_view);
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
  updatetimes();
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
  updatetimes();
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
      touri->npath=view_ntimes;
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
    updatetimes();
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

  float *eye, *aview;

  dx2 = kf->keyview_x;
  dy2 = kf->keyview_y;
  distxy2 = sqrt(dx2*dx2+dy2*dy2);
  if(distxy2<=0.0)return;
  
  eye = kf->nodeval.eye;
  aview = kf->nodeval.aview;

  dx = aview[0] - eye[0];
  dy = aview[1] - eye[1];
  dz = aview[2] - eye[2];

  distxy = sqrt(dx*dx+dy*dy);
  if(distxy<=0.0)return;

  dx2 /= distxy2;
  dy2 /= distxy2;

  dx /= distxy;
  dy /= distxy;
  dz /= distxy;

  angle_temp = 180.0*atan2(dy,dx)/PI;
  angle_temp2 = 180.0*atan2(dy2,dx2)/PI;
  az = angle_temp - angle_temp2;
  if(az>180.0)az = 360.0 - az;
  if(az<-180.0)az = az + 360.0;
  elev=atan(dz)*180.0/PI;

  kf->nodeval.elev_path=elev;
  kf->az_path=az;
  kf->az_path=az;
  *az_path=az;
  *elev_path = elev;
}

/* ------------------ adjusttourtimes ------------------------ */

void adjusttourtimes(tourdata *touri){
  // view_ntimes
  // view_tstart
  // view_tstop
  float dt;
  keyframe *keyj;
  float tstart, tstop,dtmin;
  int small_flag;


  if(touri->nkeyframes>1){
    tstart = touri->first_frame.next->noncon_time;
    tstop = touri->last_frame.prev->noncon_time;
    dtmin = (float)4*(tstop-tstart)/(float)view_ntimes;
    small_flag=0;
    for(keyj=(touri->first_frame).next->next;keyj->next!=NULL;keyj=keyj->next){
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

