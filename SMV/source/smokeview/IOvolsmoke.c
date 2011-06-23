// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#ifdef pp_GPU
#include <GL/glew.h>
#endif
#include <stdio.h>  
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "flowfiles.h"
#include <sys/types.h>
#include <sys/stat.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <float.h>
#include "egz_stdio.h"
#include "MALLOC.h"
#include "smokeheaders.h"
#include "smokeviewvars.h"
#include "update.h"
#include "interp.h"

// svn revision character string
char IOvolsmoke_revision[]="$Revision$";


/* ------------------ compute_volvals ------------------------ */

void compute_volvals(void){
  int ii;

  for(ii=0;ii<nmeshes;ii++){
    mesh *meshi;
    volrenderdata *vr;
    int iwall;
    float *alpha;
    float dstep;
    float dx, dy, dz;
    float *x, *y, *z; 
    int i, j, k;
    int count=0;
    int ibar, jbar, kbar;

    meshi = meshinfo + ii;
    vr = &(meshi->volrenderinfo);
    if(vr->loaded==0||vr->loaded==0)continue;

    x = meshi->xplt;
    y = meshi->yplt;
    z = meshi->zplt;
    ibar = meshi->ibar;
    jbar = meshi->jbar;
    kbar = meshi->kbar;
    dx = x[1] - x[0];
    dy = y[1] - y[0];
    dz = z[1] - z[0];
    dstep = sqrt(dx*dx+dy*dy+dz*dz);
    
    if(vr->smoke==NULL)continue;
    for(iwall=-3;iwall<=3;iwall++){
      float *xyz,xyzarray[3];
      int i, j;

      xyz = xyzarray;
      if(iwall==0||meshi->drawsides[iwall+3]==0)continue;
      switch (iwall){
        case 1:
        case -1:
          if(iwall<0){
            alpha=vr->alpha_yz0;
            xyz[0] = meshi->x0;
          }
          else{
            alpha=vr->alpha_yz1;
            xyz[0] = meshi->x1;
          }
          for(i=0;i<=jbar;i++){
            xyz[1] = y[i];
            for(j=0;j<=kbar;j++){
              xyz[2] = z[j];
              *alpha=optical_depth(xyz,dstep,meshi,iwall);
              alpha++;
            }
          }
          break;
        case 2:
        case -2:
          if(iwall<0){
            alpha=vr->alpha_xz0;
            xyz[1] = meshi->y0;
          }
          else{
            alpha=vr->alpha_xz1;
            xyz[1] = meshi->y1;
          }
          for(i=0;i<=ibar;i++){
            xyz[0] = x[i];
            for(j=0;j<=kbar;j++){
              xyz[2] = z[j];
              *alpha=optical_depth(xyz,dstep,meshi,iwall);
              alpha++;
            }
          }
          break;
        case 3:
        case -3:
          if(iwall<0){
            alpha=vr->alpha_xy0;
            xyz[2]=meshi->z0;
          }
          else{
            alpha=vr->alpha_xy1;
            xyz[2]=meshi->z1;
          }
          for(i=0;i<=ibar;i++){
            xyz[0] = x[i];
            for(j=0;j<=jbar;j++){
              xyz[1] = y[j];
              *alpha=optical_depth(xyz,dstep,meshi,iwall);
              alpha++;
            }
          }
          break;
      }
    }
  }
}

/* ------------------ drawsmoke3dVOL ------------------------ */

void drawsmoke3dVOLdebug(void){
  int ii;

  for(ii=0;ii<nvolfacelistinfo;ii++){
    volfacelistdata *vi;
    mesh *meshi;
    volrenderdata *vr;
    int i,j;
    float x[2], y[2], z[2];
    float *xplt, *yplt, *zplt;
    int ibar, jbar, kbar;
    char label[256];
    int iwall;

    sprintf(label,"*** %i ***",ii);

    vi = volfacelistinfoptrs[ii];
    iwall=vi->iwall;
    meshi = vi->facemesh;
    xplt = meshi->xplt;
    yplt = meshi->yplt;
    zplt = meshi->zplt;
    ibar = meshi->ibar;
    jbar = meshi->jbar;
    kbar = meshi->kbar;
    vr = &(meshi->volrenderinfo);

    if(iwall==0||meshi->drawsides[iwall+3]==0)continue;
    switch (iwall){
      case 1:
      case -1:
        if(iwall<0){
          x[0] = meshi->x0;
          x[1] = x[0];
        }
        else{
          x[0]=meshi->x1;
          x[1]=x[0];
        }
        y[0] = yplt[0];
        y[1] = yplt[jbar];
        z[0] = zplt[0];
        z[1] = zplt[kbar];
        output3Text(foregroundcolor, (x[0]+x[1])/2.0,(y[0]+y[1])/2.0,(z[0]+z[1])/2.0, label);
        break;
      case 2:
      case -2:
        if(iwall<0){
          y[0] = meshi->y0;
          y[1] = y[0];
        }
        else{
          y[0] = meshi->y1;
          y[1] = y[0];
        }
        x[0] = xplt[0];
        x[1] = xplt[ibar];
        z[0] = zplt[0];
        z[1] = zplt[kbar];
        output3Text(foregroundcolor, (x[0]+x[1])/2.0,(y[0]+y[1])/2.0,(z[0]+z[1])/2.0, label);
        break;
      case 3:
      case -3:
        if(iwall<0){
          z[0] = meshi->z0;
          z[1] = z[0];
        }
        else{
          z[0] = meshi->z1;
          z[1] = z[0];
        }
        x[0] = xplt[0];
        x[1] = xplt[ibar];
        y[0] = yplt[0];
        y[1] = yplt[jbar];
        output3Text(foregroundcolor, (x[0]+x[1])/2.0,(y[0]+y[1])/2.0,(z[0]+z[1])/2.0, label);
        break;
    }
  }
  glBegin(GL_LINES);
  for(ii=0;ii<nvolfacelistinfo;ii++){
    volfacelistdata *vi;
    mesh *meshi;
    volrenderdata *vr;
    int i,j;
    float x[2], y[2], z[2];
    float *xplt, *yplt, *zplt;
    int ibar, jbar, kbar;
    char label[256];
    int iwall;

    sprintf(label,"*** %i %%%",ii);

    vi = volfacelistinfoptrs[ii];
    iwall=vi->iwall;
    meshi = vi->facemesh;
    xplt = meshi->xplt;
    yplt = meshi->yplt;
    zplt = meshi->zplt;
    ibar = meshi->ibar;
    jbar = meshi->jbar;
    kbar = meshi->kbar;
    vr = &(meshi->volrenderinfo);

    if(iwall==0||meshi->drawsides[iwall+3]==0)continue;
    switch (iwall){
      case 1:
      case -1:
        if(iwall<0){
          x[0] = meshi->x0;
          x[1] = x[0];
          glColor3f(1.0,0.0,0.0);
        }
        else{
          x[0]=meshi->x1;
          x[1]=x[0];
          glColor3f(0.0,0.0,1.0);
        }
        y[0] = yplt[0];
        y[1] = yplt[jbar];
        z[0] = zplt[0];
        z[1] = zplt[kbar];
      //  output3Text(foregroundcolor, (x[0]+x[1])/2.0,(y[0]+y[1])/2.0,(z[0]+z[1])/2.0, label);
        glVertex3f(x[0],y[0],z[0]);
        glVertex3f(x[0],y[1],z[1]);
        glVertex3f(x[0],y[1],z[0]);
        glVertex3f(x[0],y[0],z[1]);
        break;
      case 2:
      case -2:
        if(iwall<0){
          y[0] = meshi->y0;
          y[1] = y[0];
          glColor3f(1.0,0.0,0.0);
        }
        else{
          y[0] = meshi->y1;
          y[1] = y[0];
          glColor3f(0.0,0.0,1.0);
        }
        x[0] = xplt[0];
        x[1] = xplt[ibar];
        z[0] = zplt[0];
        z[1] = zplt[kbar];
      //  output3Text(foregroundcolor, (x[0]+x[1])/2.0,(y[0]+y[1])/2.0,(z[0]+z[1])/2.0, label);
        glVertex3f(x[0],y[0],z[0]);
        glVertex3f(x[1],y[0],z[1]);
        glVertex3f(x[0],y[0],z[1]);
        glVertex3f(x[1],y[0],z[0]);
        break;
      case 3:
      case -3:
        if(iwall<0){
          z[0] = meshi->z0;
          z[1] = z[0];
          glColor3f(1.0,0.0,0.0);
        }
        else{
          z[0] = meshi->z1;
          z[1] = z[0];
          glColor3f(0.0,0.0,1.0);
        }
        x[0] = xplt[0];
        x[1] = xplt[ibar];
        y[0] = yplt[0];
        y[1] = yplt[jbar];
      //  output3Text(foregroundcolor, (x[0]+x[1])/2.0,(y[0]+y[1])/2.0,(z[0]+z[1])/2.0, label);
        glVertex3f(x[0],y[0],z[0]);
        glVertex3f(x[1],y[1],z[0]);
        glVertex3f(x[0],y[1],z[0]);
        glVertex3f(x[1],y[0],z[0]);
        break;
        break;
    }
  }
  glEnd();
}

/* ------------------ drawsmoke3dVOL ------------------------ */

void drawsmoke3dVOL(void){
  int iwall;
  float xyz[3];
  float dx, dy, dz;
  int ii;


  if(smoke3dVoldebug==1){
    drawsmoke3dVOLdebug();
  }

  if(use_transparency_data==1)transparenton();
  for(ii=0;ii<nvolfacelistinfo;ii++){
    volfacelistdata *vi;
    mesh *meshi;
    volrenderdata *vr;
    int i,j;
    float xx, yy, zz;
    float x[2], y[2], z[2];
    float *alpha;
    int n00, n01, n10, n11;
    float *xplt, *yplt, *zplt;
    int ibar, jbar, kbar;

    vi = volfacelistinfoptrs[ii];
    iwall=vi->iwall;
    meshi = vi->facemesh;
    xplt = meshi->xplt;
    yplt = meshi->yplt;
    zplt = meshi->zplt;
    ibar = meshi->ibar;
    jbar = meshi->jbar;
    kbar = meshi->kbar;
    vr = &(meshi->volrenderinfo);

    if(iwall==0||meshi->drawsides[iwall+3]==0)continue;

    glBegin(GL_TRIANGLES);
    switch (iwall){
      case 1:
      case -1:
        if(iwall<0){
          xx = meshi->x0;
          alpha = vr->alpha_yz0;
        }
        else{
          xx=meshi->x1;
          alpha = vr->alpha_yz1;
        }
        n00 = 0;
        n01 = 1;
        n10 = kbar+1;
        n11 = 1 + kbar+1;
        for(i=0;i<jbar;i++){
          y[0] = yplt[i];
          y[1] = yplt[i+1];
          for(j=0;j<kbar;j++){
            z[0] = zplt[j];
            z[1] = zplt[j+1];

            if(meshi->inside==0&&iwall>0||meshi->inside!=0&&iwall<0){
              glColor4f(0.5,0.5,0.5,alpha[n00]);
              glVertex3f(xx,y[0],z[0]);
              glColor4f(0.5,0.5,0.5,alpha[n10]);
              glVertex3f(xx,y[1],z[0]);
              glColor4f(0.5,0.5,0.5,alpha[n11]);
              glVertex3f(xx,y[1],z[1]);

              glColor4f(0.5,0.5,0.5,alpha[n00]);
              glVertex3f(xx,y[0],z[0]);
              glColor4f(0.5,0.5,0.5,alpha[n11]);
              glVertex3f(xx,y[1],z[1]);
              glColor4f(0.5,0.5,0.5,alpha[n01]);
              glVertex3f(xx,y[0],z[1]);
            }
            else{
              glColor4f(0.5,0.5,0.5,alpha[n00]);
              glVertex3f(xx,y[0],z[0]);
              glColor4f(0.5,0.5,0.5,alpha[n11]);
              glVertex3f(xx,y[1],z[1]);
              glColor4f(0.5,0.5,0.5,alpha[n10]);
              glVertex3f(xx,y[1],z[0]);

              glColor4f(0.5,0.5,0.5,alpha[n00]);
              glVertex3f(xx,y[0],z[0]);
              glColor4f(0.5,0.5,0.5,alpha[n01]);
              glVertex3f(xx,y[0],z[1]);
              glColor4f(0.5,0.5,0.5,alpha[n11]);
              glVertex3f(xx,y[1],z[1]);
            }
            alpha++;
          }
          alpha++;
        }
        break;
      case 2:
      case -2:
        n00 = 0;
        n01 = 1;
        n10 = kbar+1;
        n11 = 1 + kbar+1;
        if(iwall<0){
          alpha = vr->alpha_xz0;
          yy=meshi->y0;
        }
        else{
          alpha = vr->alpha_xz1;
          yy=meshi->y1;
        }
        for(i=0;i<ibar;i++){
          x[0] = xplt[i];
          x[1] = xplt[i+1];
          for(j=0;j<kbar;j++){
            z[0] = zplt[j];
            z[1] = zplt[j+1];
            if(meshi->inside==0&&iwall>0||meshi->inside!=0&&iwall<0){
              glColor4f(0.5,0.5,0.5,alpha[n00]);
              glVertex3f(x[0],yy,z[0]);
              glColor4f(0.5,0.5,0.5,alpha[n11]);
              glVertex3f(x[1],yy,z[1]);
              glColor4f(0.5,0.5,0.5,alpha[n10]);
              glVertex3f(x[1],yy,z[0]);

              glColor4f(0.5,0.5,0.5,alpha[n00]);
              glVertex3f(x[0],yy,z[0]);
              glColor4f(0.5,0.5,0.5,alpha[n01]);
              glVertex3f(x[0],yy,z[1]);
              glColor4f(0.5,0.5,0.5,alpha[n11]);
              glVertex3f(x[1],yy,z[1]);
            }
            else{
              glColor4f(0.5,0.5,0.5,alpha[n00]);
              glVertex3f(x[0],yy,z[0]);
              glColor4f(0.5,0.5,0.5,alpha[n10]);
              glVertex3f(x[1],yy,z[0]);
              glColor4f(0.5,0.5,0.5,alpha[n11]);
              glVertex3f(x[1],yy,z[1]);

              glColor4f(0.5,0.5,0.5,alpha[n00]);
              glVertex3f(x[0],yy,z[0]);
              glColor4f(0.5,0.5,0.5,alpha[n11]);
              glVertex3f(x[1],yy,z[1]);
              glColor4f(0.5,0.5,0.5,alpha[n01]);
              glVertex3f(x[0],yy,z[1]);
            }
            alpha++;
          }
          alpha++;
        }
        break;
      case 3:
      case -3:
        n00 = 0;
        n01 = 1;
        n10 = jbar+1;
        n11 = 1 + jbar+1;
       if(iwall<0){
          alpha = vr->alpha_xy0;
          zz=meshi->z0;
        }
        else{
          alpha = vr->alpha_xy1;
          zz=meshi->z1;
        }
        for(i=0;i<ibar;i++){
          x[0] = xplt[i];
          x[1] = xplt[i+1];
          for(j=0;j<jbar;j++){
            y[0] = yplt[j];
            y[1] = yplt[j+1];
            if(meshi->inside==0&&iwall>0||meshi->inside!=0&&iwall<0){
              glColor4f(0.5,0.5,0.5,alpha[n00]);
              glVertex3f(x[0],y[0],zz);
              glColor4f(0.5,0.5,0.5,alpha[n10]);
              glVertex3f(x[1],y[0],zz);
              glColor4f(0.5,0.5,0.5,alpha[n11]);
              glVertex3f(x[1],y[1],zz);

              glColor4f(0.5,0.5,0.5,alpha[n00]);
              glVertex3f(x[0],y[0],zz);
              glColor4f(0.5,0.5,0.5,alpha[n11]);
              glVertex3f(x[1],y[1],zz);
              glColor4f(0.5,0.5,0.5,alpha[n01]);
              glVertex3f(x[0],y[1],zz);
            }
            else{
              glColor4f(0.5,0.5,0.5,alpha[n00]);
              glVertex3f(x[0],y[0],zz);
              glColor4f(0.5,0.5,0.5,alpha[n11]);
              glVertex3f(x[1],y[1],zz);
              glColor4f(0.5,0.5,0.5,alpha[n10]);
              glVertex3f(x[1],y[0],zz);

              glColor4f(0.5,0.5,0.5,alpha[n00]);
              glVertex3f(x[0],y[0],zz);
              glColor4f(0.5,0.5,0.5,alpha[n01]);
              glVertex3f(x[0],y[1],zz);
              glColor4f(0.5,0.5,0.5,alpha[n11]);
              glVertex3f(x[1],y[1],zz);
            }
            alpha++;
          }
          alpha++;
        }
        break;
    }
    glEnd();
  }
  if(use_transparency_data==1)transparentoff();
}

#ifdef pp_GPU
/* ------------------ drawsmoke3dGPUVOL ------------------------ */

void drawsmoke3dGPUVOL(volrenderdata *vr){

#define NROWS_GPU 2
#define NCOLS_GPU 2
  int iwall;
  float xyz[3];
  float dx, dy, dz;
  mesh *meshi;

  meshi = vr->rendermesh;
  
  glUniform3f(GPUvol_eyepos,xyzeyeorig[0],xyzeyeorig[1],xyzeyeorig[2]);
  glUniform1i(GPUvol_inside,meshi->inside);
  glUniform1f(GPUvol_xyzmaxdiff,xyzmaxdiff);
  glUniform3f(GPUvol_boxmin,meshi->x0,meshi->y0,meshi->z0);
  glUniform3f(GPUvol_boxmax,meshi->x1,meshi->y1,meshi->z1);

  for(iwall=-3;iwall<=3;iwall++){
    int i,j;
    float x1, x2, y1, y2, z1, z2;

    if(iwall==0||meshi->drawsides[iwall+3]==0)continue;

    glUniform1i(GPUvol_dir,iwall);
    glBegin(GL_TRIANGLES);

    switch (iwall){
      case 1:
      case -1:
        dy = meshi->dy01/(NCOLS_GPU-1);
        dz = meshi->dz01/(NROWS_GPU-1);
        if(iwall<0){
          x1 = meshi->x0;
        }
        else{
          x1=meshi->x1;
        }
        for(i=0;i<NCOLS_GPU-1;i++){
          y1 = meshi->y0 + i*dy;
          y2 = y1 + dy;
          for(j=0;j<NROWS_GPU-1;j++){
            z1 = meshi->z0 + j*dz;
            z2 = z1 + dz;
            
            if(meshi->inside==0){
              glVertex3f(x1,y1,z1);
              glVertex3f(x1,y2,z1);
              glVertex3f(x1,y2,z2);

              glVertex3f(x1,y1,z1);
              glVertex3f(x1,y2,z2);
              glVertex3f(x1,y1,z2);
            }
            else{
              glVertex3f(x1,y1,z1);
              glVertex3f(x1,y2,z2);
              glVertex3f(x1,y2,z1);

              glVertex3f(x1,y1,z1);
              glVertex3f(x1,y1,z2);
              glVertex3f(x1,y2,z2);
            }
          }
        }
        break;
      case 2:
      case -2:
        dx = meshi->dx01/(NCOLS_GPU-1);
        dz = meshi->dz01/(NROWS_GPU-1);
        if(iwall<0){
          y1=meshi->y0;
        }
        else{
          y1=meshi->y1;
        }
        for(i=0;i<NCOLS_GPU-1;i++){
          x1 = meshi->x0 + i*dx;
          x2 = x1 + dx;
          for(j=0;j<NROWS_GPU-1;j++){
            z1 = meshi->z0 + j*dz;
            z2 = z1 + dz;

            if(meshi->inside==0){
              glVertex3f(x1,y1,z1);
              glVertex3f(x2,y1,z1);
              glVertex3f(x2,y1,z2);

              glVertex3f(x1,y1,z1);
              glVertex3f(x2,y1,z2);
              glVertex3f(x1,y1,z2);
            }
            else{
              glVertex3f(x1,y1,z1);
              glVertex3f(x2,y1,z2);
              glVertex3f(x2,y1,z1);

              glVertex3f(x1,y1,z1);
              glVertex3f(x1,y1,z2);
              glVertex3f(x2,y1,z2);
            }
          }
        }
        break;
      case 3:
      case -3:
        dx = meshi->dx01/(NCOLS_GPU-1);
        dy = meshi->dy01/(NROWS_GPU-1);
        if(iwall<0){
          z1=meshi->z0;
        }
        else{
          z1=meshi->z1;
        }
        for(i=0;i<NCOLS_GPU-1;i++){
          x1 = meshi->x0 + i*dx;
          x2 = x1 + dx;
          for(j=0;j<NROWS_GPU-1;j++){
            y1 = meshi->y0 + j*dy;
            y2 = y1 + dy;

            if(meshi->inside==0){
              glVertex3f(x1,y1,z1);
              glVertex3f(x2,y1,z1);
              glVertex3f(x2,y2,z1);

              glVertex3f(x1,y1,z1);
              glVertex3f(x2,y2,z1);
              glVertex3f(x1,y2,z1);
            }
            else{
              glVertex3f(x1,y1,z1);
              glVertex3f(x2,y2,z1);
              glVertex3f(x2,y1,z1);

              glVertex3f(x1,y1,z1);
              glVertex3f(x1,y2,z1);
              glVertex3f(x2,y2,z1);
            }
          }
        }
        break;
    }
    glEnd();
  }

}
#endif

//  glUniform3f(GPUvol_eyepos,xyzeyeorig[0],xyzeyeorig[1],xyzeyeorig[2]);
//  glUniform1i(GPUvol_inside,meshi->inside);
//  glUniform1f(GPUvol_xyzmaxdiff,xyzmaxdiff);
//  glUniform3f(GPUvol_boxmin,meshi->x0,meshi->y0,meshi->z0);
//  glUniform3f(GPUvol_boxmax,meshi->x1,meshi->y1,meshi->z1);

/* ------------------ optical_depth ------------------------ */

float optical_depth(float *xyzvert, float dstep, mesh *meshi, int iwall){
  float t_intersect, t_intersect_min=FLT_MAX, *boxmin, *boxmax;
  float tmin, tmax;
  int i;
  int nsteps;
  float dist, dx, dy, dz;
  float distseg, dxseg, dyseg, dzseg;
  float xyz[3];
  float sootdensum;
  float opacity;
  float kfactor=8700.0;
  float *vert_beg, *vert_end;
  int iwall_min=0;
  float xyzvals[3];
  int isteps;
  char *blank;

  boxmin = meshi->boxmin_scaled;
  boxmax = meshi->boxmax_scaled;

  // xyz(t) = xyzvert + t*(xyzvert - xyzeyeorig )
  // integrate from t=0 to t=t_intersect_min  (if outside mesh)
  //     ie from vertex to nearest wall along a line from the eye position
  //        intersecting the vertex position
  // integrate from t=-1 to t=0 (if inside mesh)
  //     ie from the eye position to the vertex position

  if(meshi->inside==1){
    vert_beg=xyzeyeorig;
    vert_end=xyzvert;
  }
  else{
    vert_beg=xyzvert;
    vert_end=xyzvals;

    dx = xyzvert[0] - xyzeyeorig[0];
    dy = xyzvert[1] - xyzeyeorig[1];
    dz = xyzvert[2] - xyzeyeorig[2];
    for(i=1;i<4;i++){
      int ii;
      float diffmin,diffmax,denom;

      ii=i-1;
      diffmin = boxmin[ii]-xyzvert[ii];
      diffmax = boxmax[ii]-xyzvert[ii];
      denom = xyzvert[ii]-xyzeyeorig[ii];
      if(iwall!=-i&&denom<0.0){
        t_intersect = diffmin/denom;
        if(t_intersect<t_intersect_min){
          t_intersect_min=t_intersect;
          iwall_min=-i;
        }
      }
      if(iwall!=i&&denom>0.0){
        t_intersect = diffmax/denom;
        if(t_intersect<t_intersect_min){
          t_intersect_min=t_intersect;
          iwall_min=i;
        }
      }
    }
    switch (iwall_min){
      case -1:
        vert_end[0] = boxmin[0];
        vert_end[1] = CLAMP(xyzvert[1] + t_intersect_min*dy,boxmin[1],boxmax[1]);
        vert_end[2] = CLAMP(xyzvert[2] + t_intersect_min*dz,boxmin[2],boxmax[2]);
        break;
      case 1:
        vert_end[0] = boxmax[0];
        vert_end[1] = CLAMP(xyzvert[1] + t_intersect_min*dy,boxmin[1],boxmax[1]);
        vert_end[2] = CLAMP(xyzvert[2] + t_intersect_min*dz,boxmin[2],boxmax[2]);
        break;
      case -2:
        vert_end[0] = CLAMP(xyzvert[0] + t_intersect_min*dx,boxmin[0],boxmax[0]);
        vert_end[1] = boxmin[1];
        vert_end[2] = CLAMP(xyzvert[2] + t_intersect_min*dz,boxmin[2],boxmax[2]);
        break;
      case 2:
        vert_end[0] = CLAMP(xyzvert[0] + t_intersect_min*dx,boxmin[0],boxmax[0]);
        vert_end[1] = boxmax[1];
        vert_end[2] = CLAMP(xyzvert[2] + t_intersect_min*dz,boxmin[2],boxmax[2]);
        break;
      case -3:
        vert_end[0] = CLAMP(xyzvert[0] + t_intersect_min*dx,boxmin[0],boxmax[0]);
        vert_end[1] = CLAMP(xyzvert[1] + t_intersect_min*dy,boxmin[1],boxmax[1]);
        vert_end[2] = boxmin[2];
        break;
      case 3:
        vert_end[0] = CLAMP(xyzvert[0] + t_intersect_min*dx,boxmin[0],boxmax[0]);
        vert_end[1] = CLAMP(xyzvert[1] + t_intersect_min*dy,boxmin[1],boxmax[1]);
        vert_end[2] = boxmax[2];
        break;
    }
  }

  dxseg = vert_end[0] - vert_beg[0];
  dyseg = vert_end[1] - vert_beg[1];
  dzseg = vert_end[2] - vert_beg[2];
  distseg = sqrt(dxseg*dxseg+dyseg*dyseg+dzseg*dzseg);
  if(distseg<0.001)return 0.0;

  nsteps = distseg/dstep;
  if(nsteps<1){
    nsteps=1;
  }
  sootdensum=0.0;
  isteps=0;
  if(block_volsmoke==1){
    blank=meshi->c_iblank_cell;
  }
  else{
    blank=NULL;
  }
  for(i=0;i<nsteps;i++){
    float sootden, factor;
    int icell, jcell, kcell;
    int inobst;

    factor = (0.5 + (float)i)/(float)nsteps;

    xyz[0] = (1.0-factor)*vert_beg[0] + factor*vert_end[0];
    xyz[1] = (1.0-factor)*vert_beg[1] + factor*vert_end[1];
    xyz[2] = (1.0-factor)*vert_beg[2] + factor*vert_end[2];

    sootden = interp3d(xyz, meshi->volrenderinfo.smokedata, meshi, &inobst, blank);
    if(blank!=NULL&&inobst==1)break;
    isteps++;
    sootdensum += sootden;
  }
  if(isteps!=nsteps)distseg*=(float)isteps/(float)nsteps;
  sootdensum*=xyzmaxdiff*distseg/(float)nsteps;
  //opacity = 1.0 - exp(-sootdensum);
  opacity = 1.0 - exp(-kfactor*sootdensum);
  return opacity;
}


/* ------------------ init_volrender ------------------------ */

void init_volrender(void){
  int i;
  int ijkbarmax;

  nvolrenderinfo=0;
  for(i=0;i<nmeshes;i++){
    mesh *meshi;
    volrenderdata *vr;

    meshi = meshinfo + i;
    vr = &(meshi->volrenderinfo);
    vr->rendermesh=meshi;
    vr->fire=NULL;
    vr->smoke=NULL;
    vr->loaded=0;
    vr->show=0;
    vr->timeslist=NULL;
  }
#ifndef pp_VOLRENDER
  return;
#endif
  for(i=0;i<nsliceinfo;i++){
    slice *slicei;
    char *shortlabel;
    int blocknumber;
    mesh *meshi;
    volrenderdata *vr;

    slicei = sliceinfo + i;
    blocknumber = slicei->blocknumber;
    if(blocknumber<0||blocknumber>=nmeshes)continue;
    meshi = meshinfo + blocknumber;
    if(slicei->nslicei!=meshi->ibar+1||slicei->nslicej!=meshi->jbar+1||slicei->nslicek!=meshi->kbar+1)continue;
    vr = &(meshi->volrenderinfo);
    shortlabel = slicei->label.shortlabel;
   //*** turn off loading of 3d fire slice for now
   // if(STRCMP(shortlabel,"temp")==0){  
   //   vr->fire=slicei;
   //   continue;
   // }
    if(STRCMP(shortlabel,"rho_Soot")==0){
      vr->smoke=slicei;
      continue;
    }
  }
  ijkbarmax=0;
  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi = meshinfo + i;
    if(meshi->ibar>ijkbarmax)ijkbarmax=meshi->ibar;
    if(meshi->jbar>ijkbarmax)ijkbarmax=meshi->jbar;
    if(meshi->kbar>ijkbarmax)ijkbarmax=meshi->kbar;
  }
  for(i=0;i<nmeshes;i++){
    mesh *meshi;
    volrenderdata *vr;

    meshi = meshinfo + i;
    vr = &(meshi->volrenderinfo);
    if(vr->fire!=NULL||vr->smoke!=NULL){
      int nx, ny, nz;

      nvolrenderinfo++;
      nx = ijkbarmax+1;
      ny = ijkbarmax+1;
      nz = ijkbarmax+1;
      NewMemory((void **)&vr->alpha_yz0,ny*nz*sizeof(float));
      NewMemory((void **)&vr->alpha_yz1,ny*nz*sizeof(float));
      NewMemory((void **)&vr->alpha_xz0,nx*nz*sizeof(float));
      NewMemory((void **)&vr->alpha_xz1,nx*nz*sizeof(float));
      NewMemory((void **)&vr->alpha_xy0,nx*ny*sizeof(float));
      NewMemory((void **)&vr->alpha_xy1,nx*ny*sizeof(float));
    }
    else{
      vr->alpha_yz0=NULL;
      vr->alpha_yz1=NULL;
      vr->alpha_xz0=NULL;
      vr->alpha_xz1=NULL;
      vr->alpha_xy0=NULL;
      vr->alpha_xy1=NULL;
    }
  }
  if(nvolrenderinfo>0){
    NewMemory((void **)&volfacelistinfo,6*nmeshes*sizeof(volfacelistdata));
    NewMemory((void **)&volfacelistinfoptrs,6*nmeshes*sizeof(volfacelistdata *));
  }
}
