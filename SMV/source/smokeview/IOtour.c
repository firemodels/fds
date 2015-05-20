#include "options.h"
#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include GLUT_H

#include "update.h"
#include "smokeviewvars.h"

void drawcir(float *center, float rad, float *color);
void hermiteeye(float f1, keyframe *kf1, keyframe *kf2, float *eye, float *slope);
void hermiteother(float f1, keyframe *kf1, keyframe *kf2, pathdata *pj);
void hermiteview(float t, keyframe *kf1, keyframe *kf2, float *view);
void draw_SVOBJECT(sv_object *object, int frame_index_local,propdata *prop,int recurse_level,float *rgbval,int vis_override);


/* ------------------ freetours ------------------------ */

void freetours(void){
  int i;

  if(ntours>0){
    for(i=0;i<ntours;i++){
      tourdata *touri;

      touri = tourinfo + i;
      freetour(touri);
    }
    FREEMEMORY(tourinfo);
  }
  ntours=0;
}

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
  int i;
  float *tmp_tourcol_text;
  float *tmp_tourcol_pathline;
  float *tmp_tourcol_pathknots;

  tmp_tourcol_text=tourcol_text;
  if(tourcol_text[0]<0.0)tmp_tourcol_text=foregroundcolor;
 
  tmp_tourcol_pathline=tourcol_pathline;
  if(tourcol_pathline[0]<0.0)tmp_tourcol_pathline=foregroundcolor;

  tmp_tourcol_pathknots=tourcol_pathknots;
  if(tourcol_pathknots[0]<0.0)tmp_tourcol_pathknots=foregroundcolor;

  if(selected_frame!=NULL)selectedtour_index = selected_tour-tourinfo;

  if(edittour==1){
    tourdata *tour_sel;

    /* path line (non-selected)*/

    if(showtours_whenediting==1){
      antialias(ON);
      glColor3fv(tmp_tourcol_pathline);
      glBegin(GL_LINES);
      for(i=0;i<ntours;i++){
        tourdata *touri;
        int j;

        touri = tourinfo + i;
        if(touri->display==0||touri->nkeyframes<=1||selectedtour_index==i)continue;

        for(j=0;j<view_ntimes-1;j++){
          pathdata *pj;

          pj = touri->pathnodes + j;
          glVertex3fv(pj->eye);

          pj++;
          glVertex3fv(pj->eye);
        }
      } 
     glEnd();
     antialias(OFF);
    }

    /* path line (selected)*/

    tour_sel = tourinfo + selectedtour_index;
    if(selectedtour_index!=-1&&tour_sel->display==1&&tour_sel->nkeyframes>1){
      int j;
      tourdata *touri;

      glColor3fv(tourcol_selectedpathline);
      if(tour_antialias==1)antialias(ON);
      glBegin(GL_LINES);
      touri = tourinfo + selectedtour_index;

      for(j=0;j<view_ntimes-1;j++){
        pathdata *pj;

        pj = touri->pathnodes + j;
        glVertex3fv(pj->eye);

        pj++;
        glVertex3fv(pj->eye);
      }
      glEnd();
      if(tour_antialias==1)antialias(OFF);

    }

    /* path knots */

    if(show_path_knots==1&&tour_sel->nkeyframes>1){
      glColor3f(0.0f,0.0f,1.0f);
      glPointSize(5.0f);
      glBegin(GL_POINTS);
      for(i=0;i<ntours;i++){
        int j;
        tourdata *touri;

        touri = tourinfo + i;
        if(touri->display==0||selectedtour_index!=i)continue;
   
        for(j=0;j<view_ntimes;j++){
          pathdata *pj;

          pj = touri->pathnodes + j;
          glVertex3fv(pj->eye);
        }
      }
      glEnd();

    }

    /* selected path - non-selected keyframe knots */

    if(selectedtour_index!=-1&&selected_tour->display==1){
      int j;

      glColor3fv(tourcol_selectedpathlineknots);
      glPointSize(10.0f);
      glBegin(GL_POINTS);
      for(j=0;j<tour_sel->nkeyframes;j++){
        keyframe *framej;

        framej = tour_sel->keyframe_list[j];
        if(framej->selected==1)continue;
        glVertex3fv(framej->nodeval.eye);
      }
      glEnd();
    }

    /* non-selected path, keyframe knots */

    if(showtours_whenediting==1){
      glColor3fv(tmp_tourcol_pathknots);
      glPointSize(10.0f);
      glBegin(GL_POINTS);
      for(i=0;i<ntours;i++){
        int j;
        tourdata *touri;

        touri = tourinfo + i;
        if(touri->display==0||i==selectedtour_index)continue;

        for(j=0;j<touri->nkeyframes;j++){
          keyframe *framej;

          framej = touri->keyframe_list[j];
          if(framej->selected==1)continue;
          glVertex3fv(framej->nodeval.eye);
        }
      }
      glEnd();
    }

    /* selected path, selected keyframe */

    if(selected_frame!=NULL&&selected_tour->display==1){
      glBegin(GL_POINTS);
      glColor3fv(tourcol_selectedknot);
      glVertex3fv(selected_frame->nodeval.eye);
      if(selected_frame->viewtype==ABS_VIEW){
        glColor3fv(tourcol_selectedview);
        glVertex3fv(selected_frame->nodeval.xyz_view_abs);
      }
      glEnd();

    }

    /* keyframe times */
    SNIFF_ERRORS("after select path, selected keyframe");
    CheckMemory;

    if(fontindex==SCALED_FONT)scale_3dfont();
    for(i=0;i<ntours;i++){
      int j;
      tourdata *touri;

      touri = tourinfo + i;
      if(touri->display==0)continue;
      if(showtours_whenediting==0&&selectedtour_index!=i)continue;

      for(j=0;j<touri->nkeyframes;j++){
        keyframe *framej;
        float *eye;

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
    switch(tourlocus_type){
      case 0:
        antialias(ON);
        glBegin(GL_LINES);
        glColor3fv(tourcol_avatar);
        for(i=0;i<ntours;i++){
          tourdata *touri;
          pathdata *pj;
          int iframe_local;

          touri = tourinfo + i;
          if(touri->display==0||touri->nkeyframes<=1)continue;
          if(touri->timeslist==NULL)continue;

          iframe_local = touri->timeslist[itimes];
          pj = touri->pathnodes + iframe_local;
          if(keyframe_snap==1)pj = pj->keysnap;

          glVertex3fv(pj->eye);
          glVertex3f(pj->eye[0],pj->eye[1],pj->eye[2]+0.1);

          glVertex3fv(pj->eye);
          glVertex3fv(pj->tour_view);
        }
        glEnd();
        antialias(OFF);
        break;
      case 1:
        for(i=0;i<ntours;i++){
          tourdata *touri;
          pathdata *pj;
          int iframe_local;

          touri = tourinfo + i;
          if(touri->display==0||touri->nkeyframes<=1)continue;
          if(touri->timeslist==NULL)continue;

          iframe_local = touri->timeslist[itimes];
          pj = touri->pathnodes + iframe_local;
          if(keyframe_snap==1)pj = pj->keysnap;


          drawcir(pj->eye,tourrad_avatar,tourcol_avatar);

        }
        break;
      case 2:
        for(i=0;i<ntours;i++){
          tourdata *touri;
          float *tour_view, az_angle, dxy[2];
          pathdata *pj;
          float *eye;
          int iframe_local;

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

          draw_SVOBJECT(avatar_types[touri->glui_avatar_index],0,NULL,0,NULL,0);
          glPopMatrix();
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
  }

}

/* ------------------ drawcir ------------------------ */

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
  float *eye;

  glPointSize(20.0f);
  glBegin(GL_POINTS);
  for(i=0;i<ntours;i++){
    touri = tourinfo + i;
    for(j=0;j<touri->nkeyframes;j++){
  
      if(showtours_whenediting==1||selectedtour_index==i){

        if(touri->display==1){
          unsigned char r, g, b;

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
  int i;
  
  //  keyframe *framejm1;

  keyframe **tourknotskeylist_copy;
  tourdata **tourknotstourlist_copy;

  // construct keyframe list for selecting keyframes

  ntourknots=0;
  if(ntours==0)return;

  for(i=0;i<ntours;i++){
    tourdata *touri;

    touri = tourinfo + i;
    adjusttourtimes(touri);
  }

  for(i=0;i<ntours;i++){
    tourdata *touri;
    keyframe *keyj;
    int nframes;

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
    keyframe *kf1, *kf2;
    float *tour_dist3a;
    float factor,total_time;
    float total_distance;
    float vdist,vtime2,vdt;
    float tour_tstart, tour_tstop;
    int j,jj;
    int iframe_local;
    int ntotal,ntotal2;
    int iframe_old, iframe_new;

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
      float *lasteye, *thiseye, *nexteye;
      float *xyz_view0, *xyz_view1, *xyz_view2;
      float s1, s2, d1, d2;
      float a, b, c;

      *tourknotskeylist_copy++ = keyj;
      *tourknotstourlist_copy++ = touri;
      keyj->selected=0;
      keyj->distance=0.0;
      keyj->npoints=0;

      lastkey=keyj->prev;
      thiskey=keyj;
      nextkey=keyj->next;
      if(touri->periodic==1){
        if(j==0){
          lastkey=touri->keyframe_list[touri->nkeyframes-2];
        }
        if(j==touri->nkeyframes-1){
          nextkey=touri->keyframe_list[1];
        }
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
      xyz_view0 = lastkey->nodeval.xyz_view_abs;
      xyz_view1 = thiskey->nodeval.xyz_view_abs;
      xyz_view2 = nextkey->nodeval.xyz_view_abs;

      if(touri->periodic==0&&j==0){
        keyj->s_eye[0]=0.0;
        keyj->s_eye[1]=0.0;
        keyj->s_eye[2]=0.0;

        keyj->s_az=0.0;
        keyj->s_elev=0.0;
        keyj->s_zoom=0.0;

        VECDIFF3(keyj->d_eye,nexteye,thiseye);

        keyj->d_az=nextkey->az_path - thiskey->az_path;
        keyj->d_zoom=nextkey->nodeval.zoom - thiskey->nodeval.zoom;
        keyj->d_elev=nextkey->nodeval.elev_path - thiskey->nodeval.elev_path;

        keyj->s_xyz_view[0]=0.0;
        keyj->s_xyz_view[1]=0.0;
        keyj->s_xyz_view[2]=0.0;

        VECDIFF3(keyj->d_xyz_view,xyz_view2,xyz_view1);
      }
      else if(touri->periodic==0&&j==touri->nkeyframes-1){
        VECDIFF3(keyj->s_eye,thiseye,lasteye);

        keyj->s_az  =thiskey->az_path           - lastkey->az_path;
        keyj->s_zoom=thiskey->nodeval.zoom      - lastkey->nodeval.zoom;
        keyj->s_elev=thiskey->nodeval.elev_path - lastkey->nodeval.elev_path;

        keyj->d_eye[0]=0.0;
        keyj->d_eye[1]=0.0;
        keyj->d_eye[2]=0.0;

        keyj->d_az=0.0;
        keyj->d_zoom=0.0;
        keyj->d_elev=0.0;

        VECDIFF3(keyj->s_xyz_view,xyz_view1,xyz_view0);

        keyj->d_xyz_view[0]=0.0;
        keyj->d_xyz_view[1]=0.0;
        keyj->d_xyz_view[2]=0.0;
      }
      else{
#ifdef xxx
        float del1, del2;
#endif
        float sfactor, dfactor;

#ifdef xxx
        del1 = thiskey->nodeval.time - lastkey->nodeval.time;
        del2 = nextkey->nodeval.time - thiskey->nodeval.time;
        sfactor = 2*del2/(del1 + del2);
        dfactor = 2*del1/(del1 + del2);
#endif
        sfactor = 1.0;
        dfactor = 1.0;

#define HERM1(sfactor,s1,s2,lastval,thisval,nextval,val)\
        val[0]=sfactor*(s1*(thisval[0] - lastval[0]) + s2*(nextval[0]-thisval[0]));\
        val[1]=sfactor*(s1*(thisval[1] - lastval[1]) + s2*(nextval[1]-thisval[1]));\
        val[2]=sfactor*(s1*(thisval[2] - lastval[2]) + s2*(nextval[2]-thisval[2]))

        HERM1(sfactor,s1,s2,lasteye,thiseye,nexteye,keyj->s_eye);
        keyj->s_az  =sfactor*(s1*(thiskey->az_path -           lastkey->az_path) +           s2*(nextkey->az_path-          thiskey->az_path));
        keyj->s_zoom=sfactor*(s1*(thiskey->nodeval.zoom -      lastkey->nodeval.zoom) +      s2*(nextkey->nodeval.zoom-     thiskey->nodeval.zoom));
        keyj->s_elev=sfactor*(s1*(thiskey->nodeval.elev_path - lastkey->nodeval.elev_path) + s2*(nextkey->nodeval.elev_path-thiskey->nodeval.elev_path));
        HERM1(sfactor,s1,s2,xyz_view0,xyz_view1,xyz_view2,keyj->s_xyz_view);

        HERM1(dfactor,d1,d2,lasteye,thiseye,nexteye,keyj->d_eye);
        keyj->d_az  =dfactor*(d1*(thiskey->az_path -           lastkey->az_path) +           d2*(nextkey->az_path-          thiskey->az_path));
        keyj->d_zoom=dfactor*(d1*(thiskey->nodeval.zoom -      lastkey->nodeval.zoom) +      d2*(nextkey->nodeval.zoom-     thiskey->nodeval.zoom));
        keyj->d_elev=dfactor*(d1*(thiskey->nodeval.elev_path - lastkey->nodeval.elev_path) + d2*(nextkey->nodeval.elev_path-thiskey->nodeval.elev_path));
        HERM1(dfactor,d1,d2,xyz_view0,xyz_view1,xyz_view2,keyj->d_xyz_view);
      }
    }

    for(keyj=(touri->first_frame).next;keyj->next!=NULL;keyj=keyj->next){
      float denom;

      keyj->keyview_xyz[0]=keyj->d_eye[0];
      keyj->keyview_xyz[1]=keyj->d_eye[1];
      keyj->keyview_xyz[2]=0.0;
      if(keyj->viewtype==ABS_VIEW)xyzview2azelev(keyj,NULL,NULL);

      ROTATE(keyj->keyview_xyz2,keyj->keyview_xyz,keyj->az_path*DEG2RAD);
      keyj->keyview_xyz2[2]=0.0;
      denom=NORM2(keyj->keyview_xyz2);
      if(denom==0.0)continue;
      VEC2MA(keyj->keyview_xyz2,10000.0/denom);
    }

    // evaluate quantities along path - determine distances
    // define tour_t and tour_dist (tour_t is uniform )

    iframe_local=0;
    tour_dist[0]=0.0;
    for(j=0;j<view_ntimes;j++){
      pathdata *pj,*pjm1;
      float *eye, *xyz_view, *tour_view;
      float f1, f2, dt;
      float view_local[3];
      float vtime;

      pj = touri->pathnodes + j;
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
      xyz_view=pj->xyz_view_abs;
      tour_view=pj->tour_view;

      if(kf1->nodeval.eye[0]==kf2->nodeval.eye[0]&&
         kf1->nodeval.eye[1]==kf2->nodeval.eye[1]&&
         kf1->nodeval.eye[2]==kf2->nodeval.eye[2]){
        hermiteeye(1.0,kf1->prev,kf1,eye,view_local);
        slerp(kf1->prev->keyview_xyz2,kf1->keyview_xyz2,1.0,view_local);
      }
      else{
        hermiteeye(f1,kf1,kf2,eye,view_local);
        slerp(kf1->keyview_xyz2,kf2->keyview_xyz2,f1,view_local);
      }

      hermiteview(f1,kf1,kf2, xyz_view);
      hermiteother(f1,kf1,kf2,pj);

      tour_view[0]=view_local[0];
      tour_view[1]=view_local[1];
      tour_view[2]=0.0;
      tour_t[j]=vtime;
      if(j!=0){
        float dx, dy;
        float dz, distance;

        pjm1 = pj - 1;
        dx = eye[0]-pjm1->eye[0];
        dy = eye[1]-pjm1->eye[1];
        dz = eye[2]-pjm1->eye[2];
        distance = sqrt(dx*dx+dy*dy+dz*dz);
        tour_dist[j]=tour_dist[j-1] + distance;
        kf1->distance += distance;
      }
    }

    // construct running "total_distance" info

    kf1 = touri->first_frame.next;
    kf1->total_distance=0.0;
    total_distance=0.0;
    for(keyj=kf1->next;keyj->next!=NULL;keyj=keyj->next){
      keyj->total_distance = keyj->prev->total_distance + keyj->prev->distance;
      total_distance += keyj->prev->distance;
    }

    // decide fraction of global path vs. local path

    touri->global_dist=0.0;
    touri->local_dist=0.0;
    factor=0.0;
    total_time=0.0;
    for(keyj=(touri->first_frame).next;keyj->next!=NULL;keyj=keyj->next){
      if(keyj->next->next!=NULL)total_time += keyj->next->noncon_time - keyj->noncon_time;
      if(tour_constant_vel==1){
        touri->global_dist+=keyj->distance;
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
        keyj->npoints=view_ntimes*(1.0-factor)*keyj->distance/touri->global_dist;
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
        tour_dist3a[jj] = keyj->total_distance + keyj->distance*(float)j/(float)keyj->npoints;
        jj++;
      }
    }

    // average tour_dist3 array and copy into tour_dist2

    for(j=0;j<5;j++){
      tour_dist3a[-1-j]=tour_dist3a[-j]-(tour_dist3a[1]-tour_dist3a[0]);
      tour_dist3a[view_ntimes+j]=tour_dist3a[view_ntimes-1+j]+(tour_dist3a[view_ntimes-1]-tour_dist3a[view_ntimes-2]);
    }
    for(j=0;j<view_ntimes;j++){
      float avgsum;

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
      float f1, f2;

      vdist = tour_dist2[j];
      iframe_local = isearch(tour_dist,view_ntimes,vdist,iframe_local);
      f1 = (vdist-tour_dist[iframe_local])/(tour_dist[iframe_local+1]-tour_dist[iframe_local]);
      f2 = 1 - f1;
      tour_t2[j] = f2*tour_t[iframe_local] + f1*tour_t[iframe_local+1] ;
    }
    iframe_old=-1;
    for(j=0;j<view_ntimes;j++){
      pathdata *pj;
      float *eye, *xyz_view, *tour_view;
      float f1, dt;
      float view_local[3];
      float vtime;

      pj = touri->pathnodes + j;
      vtime = tour_t2[j];
      vtime2 = touri->keyframe_list[0]->nodeval.time + j*vdt;
      iframe_new = isearch(touri->keyframe_times,touri->nkeyframes,vtime,iframe_old);
      kf1 = touri->keyframe_list[iframe_new];
      kf2 = touri->keyframe_list[iframe_new+1];
      dt = kf2->nodeval.time - kf1->nodeval.time;
      f1 = (vtime - kf1->nodeval.time)/dt;
      if(f1<0.0)f1=0.0;
      if(f1>1.0)f1=1.0;
      pj->time=vtime2;
      touri->path_times[j]=vtime2;

      eye=pj->eye;
      xyz_view=pj->xyz_view_abs;
      tour_view=pj->tour_view;

      hermiteeye(f1,kf1,kf2,eye,view_local);
      hermiteother(f1,kf1,kf2,pj);
      hermiteview(f1,kf1,kf2,xyz_view);

      if(kf1->viewtype==REL_VIEW||kf2->viewtype==REL_VIEW){
        float az;
        float dxyz[3], denom, dxyz2[3];

        dxyz[0] = view_local[0];
        dxyz[1] = view_local[1];
        denom = 10.0*NORM2(dxyz);
        if(denom==0.0)denom=1.0;
        dxyz[0] /= denom;
        dxyz[1] /= denom;
        az = pj->az_path*DEG2RAD;
        ROTATE(dxyz2,dxyz,az);
        dxyz2[2] = tan(pj->elev_path*DEG2RAD)/10.0;
        VECADD3(tour_view,eye,dxyz2);
      } 
      else{
        float dxyz[3], denom;

        VECDIFF3(dxyz,xyz_view,eye);
        denom = 10.0*NORM3(dxyz);
        dxyz[0] /= denom;
        dxyz[1] /= denom;
        dxyz[2] /= denom;
        VECADD3(tour_view,eye,dxyz);
      }
      if(iframe_old!=iframe_new){
        iframe_old=iframe_new;
        pj->keysnap->tour_view[0]=tour_view[0];
        pj->keysnap->tour_view[1]=tour_view[1];
        pj->keysnap->tour_view[2]=tour_view[2];
      }
    }
    for(keyj=kf1->next;keyj->next!=NULL;keyj=keyj->next){
      keyj->nodeval.time = tour_tstart + (tour_tstop-tour_tstart)*keyj->total_distance/total_distance;
    }
    if(selected_frame!=NULL)selected_frame->selected=1;
  }
}

#define HERMVAL() ((2.0*t3-3.0*t2+1.0)*p0 + (t3-2.0*t2+t)*m0 + (t3-t2)*m1 + (-2.0*t3+3.0*t2)*p1)
#define HERMDERIV() ((6.0*t2-6.0*t)*p0 + (3.0*t2-4.0*t+1.0)*m0 + (3.0*t2-2.0*t)*m1 + (-6.0*t2+6.0*t)*p1)

/* ------------------ hermiteye ------------------------ */

void hermiteeye(float t, keyframe *kf1, keyframe *kf2, float *eye, float *slope){
  int i;
  float t3, t2;

  t2 = t*t;
  t3 = t2*t;

  for(i=0;i<3;i++){
    float p0, p1, m0, m1;

    p0 = kf1->nodeval.eye[i];
    p1 = kf2->nodeval.eye[i];
    m0 = kf1->d_eye[i];
    m1 = kf2->s_eye[i];

    eye[i] = HERMVAL();
    if(i!=2)slope[i] = HERMDERIV();
  }
}


/* ------------------ hermiteother ------------------------ */

void hermiteother(float t, keyframe *kf1, keyframe *kf2, pathdata *pj){
  float p0, p1, m0, m1;
  float t3, t2;

  t2 = t*t;
  t3 = t2*t;

  p0 = kf1->az_path;
  p1 = kf2->az_path;
  m0 = kf1->d_az;
  m1 = kf2->s_az;
  pj->az_path = HERMVAL();

  p0 = kf1->nodeval.zoom;
  p1 = kf2->nodeval.zoom;
  m0 = kf1->d_zoom;
  m1 = kf2->s_zoom;
  pj->zoom = HERMVAL();

  p0=kf1->nodeval.elev_path;
  p1=kf2->nodeval.elev_path;
  m0 = kf1->d_elev;
  m1 = kf2->s_elev;
  pj->elev_path = HERMVAL();
}

/* ------------------ hermiteview ------------------------ */

void hermiteview(float t, keyframe *kf1, keyframe *kf2, float *view){
  int i;

  for(i=0;i<3;i++){
    float p0, p1, m0, m1;
    float t3, t2;

    p0 = kf1->nodeval.xyz_view_abs[i];
    p1 = kf2->nodeval.xyz_view_abs[i];
    m0 = kf1->d_xyz_view[i];
    m1 = kf2->s_xyz_view[i];
    t2 = t*t;
    t3 = t2*t;
    view[i] = HERMVAL();
  }
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
  fxyz_view = frame->nodeval.xyz_view_abs;
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
  frame->keyview_xyz[0]=0.0;
  frame->keyview_xyz[1]=0.0;
  frame->keyview_xyz[2]=0.0;

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
    strcpy(touri->label,trim_front(label));
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
    if(viewalltours==1)TourMenu(MENU_TOUR_SHOWALL);
  }
}

/* ------------------ xyzview2azelev ------------------------ */

void xyzview2azelev(keyframe *kf, float *az_path, float *elev_path){
  float dxyz[3];
  float dxy2[2];
  float distxy, distxy2;
  float angle_temp, angle_temp2;
  float az, elev;

  float *eye, *xyz_view;

  dxy2[0] = kf->keyview_xyz[0];
  dxy2[1] = kf->keyview_xyz[1];
  distxy2 = NORM2(dxy2);
  if(distxy2<=0.0)return;
  
  eye = kf->nodeval.eye;
  xyz_view = kf->nodeval.xyz_view_abs;

  VECDIFF3(dxyz,xyz_view,eye);

  distxy = NORM2(dxyz);
  if(distxy<=0.0)return;

  VEC2DA(dxy2,distxy2);
  VEC3DA(dxyz,distxy);

  angle_temp = RAD2DEG*atan2(dxyz[1],dxyz[0]);
  angle_temp2 = RAD2DEG*atan2(dxy2[1],dxy2[0]);
  az = angle_temp - angle_temp2;
  if(az>180.0)az = 360.0 - az;
  if(az<-180.0)az = az + 360.0;
  elev=atan(dxyz[2])*RAD2DEG;

  kf->az_path=az;
  kf->nodeval.elev_path=elev;
  if(az_path!=NULL)*az_path=az;
  if(elev_path!=NULL)*elev_path = elev;
}

/* ------------------ adjusttourtimes ------------------------ */

void adjusttourtimes(tourdata *touri){
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

