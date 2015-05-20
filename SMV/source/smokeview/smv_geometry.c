// $Date: 2015-05-01 19:51:23 -0400 (Fri, 01 May 2015) $ 
// $Revision: 22581 $
// $Author: gforney $

// svn revision character string
char smv_geometry_revision[]="$Revision: 22581 $";

#include "options.h"
#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include GLUT_H

#include "smokeviewvars.h"

/* ------------------ slerp ------------------------ */

void slerp(float *p0, float *p1, float t, float *pout){
  float cosangle,sinangle,denom,angle,factor1,factor2;

  denom = NORM3(p0)*NORM3(p1);
  if(denom==0.0){
    pout[0]=p0[0];
    pout[1]=p0[1];
    pout[2]=p0[2];
    return;
  }
  cosangle = DOT3(p0,p1)/denom;
  angle = acos(cosangle);
  sinangle = sin(angle);
  factor1 = sin((1.0-t)*angle)/sinangle;
  factor2 = sin(t*angle)/sinangle;
  pout[0]=factor1*p0[0]+factor2*p1[0];
  pout[1]=factor1*p0[1]+factor2*p1[1];
  pout[2]=factor1*p0[2]+factor2*p1[2];
}

/* ----------------------- drawtetra_outline ----------------------------- */

void drawtetra_outline(float *v1, float *v2, float *v3, float *v4, unsigned char *rgbcolor){
  glBegin(GL_LINES);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);
  glVertex3fv(v1);
  glVertex3fv(v2);
  glVertex3fv(v2);
  glVertex3fv(v3);
  glVertex3fv(v3);
  glVertex3fv(v1);
  glVertex3fv(v1);
  glVertex3fv(v4);
  glVertex3fv(v2);
  glVertex3fv(v4);
  glVertex3fv(v3);
  glVertex3fv(v4);
  glEnd();
}

/* ----------------------- drawfilledtetra ----------------------------- */

void drawfilledtetra(float *v1, float *v2, float *v3, float *v4, unsigned char *rgbcolor){
  float diff1[3],diff2[3],cross[3];

  glBegin(GL_TRIANGLES);
  if(rgbcolor!=NULL)glColor3ubv(rgbcolor);

  VECDIFF3(diff1,v1,v2);
  VECDIFF3(diff2,v4,v2);
  CROSS(cross,diff1,diff2);
  glNormal3f(cross[0],cross[1],cross[2]);
  glVertex3fv(v1);
  glVertex3fv(v2);
  glVertex3fv(v4);

  VECDIFF3(diff1,v1,v4);
  VECDIFF3(diff2,v2,v4);
  CROSS(cross,diff1,diff2);
  glNormal3f(cross[0],cross[1],cross[2]);
  glVertex3fv(v1);
  glVertex3fv(v4);
  glVertex3fv(v2);

  VECDIFF3(diff1,v2,v3);
  VECDIFF3(diff2,v4,v3);
  CROSS(cross,diff1,diff2);
  glNormal3f(cross[0],cross[1],cross[2]);
  glVertex3fv(v2);
  glVertex3fv(v3);
  glVertex3fv(v4);

  VECDIFF3(diff1,v2,v4);
  VECDIFF3(diff2,v3,v4);
  CROSS(cross,diff1,diff2);
  glNormal3f(cross[0],cross[1],cross[2]);
  glVertex3fv(v2);
  glVertex3fv(v4);
  glVertex3fv(v3);

  VECDIFF3(diff1,v4,v1);
  VECDIFF3(diff2,v3,v4);
  CROSS(cross,diff1,diff2);
  glNormal3f(cross[0],cross[1],cross[2]);
  glVertex3fv(v1);
  glVertex3fv(v4);
  glVertex3fv(v3);

  VECDIFF3(diff1,v1,v3);
  VECDIFF3(diff2,v4,v3);
  CROSS(cross,diff1,diff2);
  glNormal3f(cross[0],cross[1],cross[2]);
  glVertex3fv(v1);
  glVertex3fv(v3);
  glVertex3fv(v4);

  VECDIFF3(diff1,v1,v3);
  VECDIFF3(diff2,v2,v3);
  CROSS(cross,diff1,diff2);
  glNormal3f(cross[0],cross[1],cross[2]);
  glVertex3fv(v1);
  glVertex3fv(v3);
  glVertex3fv(v2);

  VECDIFF3(diff1,v1,v2);
  VECDIFF3(diff2,v3,v2);
  CROSS(cross,diff1,diff2);
  glNormal3f(cross[0],cross[1],cross[2]);
  glVertex3fv(v1);
  glVertex3fv(v2);
  glVertex3fv(v3);
  glEnd();
}

/* ----------------------- drawfilled2tetra ----------------------------- */

void drawfilled2tetra(float *v1, float *v2, float *v3, float *v4, 
                     unsigned char *rgb0color,
                     unsigned char *rgb1color,
                     unsigned char *rgb2color,
                     unsigned char *rgb3color,
                     int *vis_plane
                     ){
  float diff1[3],diff2[3],cross[3];

  glBegin(GL_TRIANGLES);
  if(vis_plane[0]==1){
     if(rgb0color!=NULL)glColor3ubv(rgb0color);
    VECDIFF3(diff1,v1,v2);
    VECDIFF3(diff2,v4,v2);
    CROSS(cross,diff2,diff1);
    glNormal3f(cross[0],cross[1],cross[2]);
    glVertex3fv(v1);
    glVertex3fv(v2);
    glVertex3fv(v4);

    VECDIFF3(diff1,v1,v4);
    VECDIFF3(diff2,v2,v4);
    CROSS(cross,diff2,diff1);
    glNormal3f(cross[0],cross[1],cross[2]);
    glVertex3fv(v1);
    glVertex3fv(v4);
    glVertex3fv(v2);
  }

  if(vis_plane[1]==1){
    if(rgb1color!=NULL)glColor3ubv(rgb1color);
    VECDIFF3(diff1,v2,v3);
    VECDIFF3(diff2,v4,v3);
    CROSS(cross,diff2,diff1);
    glNormal3f(cross[0],cross[1],cross[2]);
    glVertex3fv(v2);
    glVertex3fv(v3);
    glVertex3fv(v4);

    VECDIFF3(diff1,v2,v4);
    VECDIFF3(diff2,v3,v4);
    CROSS(cross,diff2,diff1);
    glNormal3f(cross[0],cross[1],cross[2]);
    glVertex3fv(v2);
    glVertex3fv(v4);
    glVertex3fv(v3);
  }

  if(vis_plane[2]==1){
    if(rgb2color!=NULL)glColor3ubv(rgb2color);
    VECDIFF3(diff1,v4,v1);
    VECDIFF3(diff2,v3,v4);
    CROSS(cross,diff2,diff1);
    glNormal3f(cross[0],cross[1],cross[2]);
    glVertex3fv(v1);
    glVertex3fv(v4);
    glVertex3fv(v3);

    VECDIFF3(diff1,v1,v3);
    VECDIFF3(diff2,v4,v3);
    CROSS(cross,diff2,diff1);
    glNormal3f(cross[0],cross[1],cross[2]);
    glVertex3fv(v1);
    glVertex3fv(v3);
    glVertex3fv(v4);
  }

  if(vis_plane[3]==1){
    if(rgb3color!=NULL)glColor3ubv(rgb3color);
    VECDIFF3(diff1,v1,v3);
    VECDIFF3(diff2,v2,v3);
    CROSS(cross,diff2,diff1);
    glNormal3f(cross[0],cross[1],cross[2]);
    glVertex3fv(v1);
    glVertex3fv(v3);
    glVertex3fv(v2);

    VECDIFF3(diff1,v1,v2);
    VECDIFF3(diff2,v3,v2);
    CROSS(cross,diff2,diff1);
    glNormal3f(cross[0],cross[1],cross[2]);
    glVertex3fv(v1);
    glVertex3fv(v2);
    glVertex3fv(v3);
  }
  glEnd();
}

/* ------------------ getmesh_zcell ------------------------ */

float getmesh_zcell(mesh *meshi, float xval, float yval, int *valid){
  float *xplt, *yplt,*zcell;
  float dx, dy;
  int ibar, jbar;
  int ival, jval;
  float zval;
  int nxcell;

  xplt = meshi->xplt_orig;
  yplt = meshi->yplt_orig;
  ibar = meshi->ibar;
  jbar = meshi->jbar;
  nxcell=ibar;
  *valid=0;
  if(xval<xplt[0]||xval>xplt[ibar])return 0.0;
  if(yval<yplt[0]||yval>yplt[jbar])return 0.0;

  dx = xplt[1]-xplt[0];
  dy = yplt[1]-yplt[0];
  ival = (xval-xplt[0])/dx;
  if(ival>=ibar)ival=ibar-1;
  jval = (yval-yplt[0])/dy;
  if(jval>=jbar)jval=jbar-1;
  zcell = meshi->zcell;
  zval = zcell[IJCELL2(ival,jval)];
  *valid=1;
  return zval;
}

/* ------------------ compare_float ------------------------ */

int compare_floats( const void *arg1, const void *arg2 ){
  float x, y;
  x=*(float *)arg1;
  y=*(float *)arg2;
  if( x< y)return -1;
  if( x> y)return 1;
  return 0;
}

/* ------------------ removedupfloats ------------------------ */

void removedupfloats(float **valsptr, int *nvals,int *ivals, float dval_min){
  int nv;
  int i,ii;
  float *vals,valmid;
  
  *ivals=0;
  if(*nvals==0)return;
  nv = *nvals;
  vals = *valsptr;
  qsort( (float *)vals, (size_t)nv, sizeof( float ), compare_floats );
  ii=1;
  for(i=1;i<nv;i++){
    if(ABS(vals[i]-vals[i-1])<=dval_min)continue;
    vals[ii]=vals[i];
    ii++;
  }
  valmid=(vals[0]+vals[*nvals-1])/2.0;
  if(*nvals!=ii){
    *nvals=ii;
    ResizeMemory((void **)&vals,*nvals*sizeof(float));
    *valsptr=vals;
  }
  for(i=1;i<*nvals;i++){
    if(vals[i-1]<=valmid&&valmid<=vals[i]){
      *ivals=i;
      break;
    }
  }
}

/* ------------------ closest_nodeindex ------------------------ */

int closest_nodeindex(float val,float *vals,int nvals, float eps){
  int j;

  if(val<vals[0])return -1;
  if(val>vals[nvals-1])return -1;
  for(j=0;j<nvals-1;j++){
    if(vals[j]<=val&&val<vals[j+1])return j;
  }
  return nvals-1;
}

/* ------------------ update_plot_alls ------------------------ */

void update_plotxyz_all(void){
  int i;
  float *xp, *yp, *zp;
  float dxyz_min=100000.0;

  FREEMEMORY(plotx_all);
  FREEMEMORY(ploty_all);
  FREEMEMORY(plotz_all);
  nplotx_all=0;
  nploty_all=0;
  nplotz_all=0;
  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi = meshinfo + i;
    nplotx_all+=(meshi->ibar+1);
    nploty_all+=(meshi->jbar+1);
    nplotz_all+=(meshi->kbar+1);
  }
  NewMemory((void **)&plotx_all,nplotx_all*sizeof(float));
  NewMemory((void **)&ploty_all,nploty_all*sizeof(float));
  NewMemory((void **)&plotz_all,nplotz_all*sizeof(float));
  xp = plotx_all;
  yp = ploty_all;
  zp = plotz_all;
  for(i=0;i<nmeshes;i++){
    int j;
    mesh *meshi;

    meshi = meshinfo + i;
    for(j=0;j<meshi->ibar+1;j++){
      *xp++ = meshi->xplt[j];
    }
    for(j=0;j<meshi->jbar+1;j++){
      *yp++ = meshi->yplt[j];
    }
    for(j=0;j<meshi->kbar+1;j++){
      *zp++ = meshi->zplt[j];
    }
    for(j=1;j<meshi->ibar+1;j++){
      float dxyz;

      dxyz = meshi->xplt[j]-meshi->xplt[j-1];
      dxyz_min = MIN(dxyz_min,dxyz);
    }
    for(j=1;j<meshi->jbar+1;j++){
      float dxyz;

      dxyz = meshi->yplt[j]-meshi->yplt[j-1];
      dxyz_min = MIN(dxyz_min,dxyz);
    }
    for(j=1;j<meshi->kbar+1;j++){
      float dxyz;

      dxyz = meshi->zplt[j]-meshi->zplt[j-1];
      dxyz_min = MIN(dxyz_min,dxyz);
    }
  }
  dxyz_min/=10.0;
  removedupfloats(&plotx_all,&nplotx_all,&iplotx_all,dxyz_min);
  removedupfloats(&ploty_all,&nploty_all,&iploty_all,dxyz_min);
  removedupfloats(&plotz_all,&nplotz_all,&iplotz_all,dxyz_min);
  for(i=0;i<nmeshes;i++){
    mesh *meshi;
    int j;

    meshi = meshinfo + i;
    NewMemory((void **)&meshi->iplotx_all,nplotx_all*sizeof(int));
    NewMemory((void **)&meshi->iploty_all,nploty_all*sizeof(int));
    NewMemory((void **)&meshi->iplotz_all,nplotz_all*sizeof(int));

    for(j=0;j<nplotx_all;j++){
      float val;
      int ival;

      meshi->iplotx_all[j]=-1;
      val = plotx_all[j];
      ival = closest_nodeindex(val,meshi->xplt,meshi->ibar+1,dxyz_min);
      if(ival<0)continue;
      meshi->iplotx_all[j]=ival;
    }
    for(j=0;j<nploty_all;j++){
      float val;
      int ival;

      meshi->iploty_all[j]=-1;
      val = ploty_all[j];
      ival = closest_nodeindex(val,meshi->yplt,meshi->jbar+1,dxyz_min);
      if(ival<0)continue;
      meshi->iploty_all[j]=ival;
    }
    for(j=0;j<nplotz_all;j++){
      float val;
      int ival;

      meshi->iplotz_all[j]=-1;
      val = plotz_all[j];
      ival = closest_nodeindex(val,meshi->zplt,meshi->kbar+1,dxyz_min);
      if(ival<0)continue;
      meshi->iplotz_all[j]=ival;
    }
  }
}

#define MESHEPS 0.001

/* ------------------ getmesh ------------------------ */

mesh *getmesh(float *xyz){
  int i;

  for(i=0;i<nmeshes;i++){
    mesh *meshi;
    int ibar, jbar, kbar;
    float *xplt, *yplt, *zplt;

    meshi = meshinfo+i;

    ibar = meshi->ibar;
    jbar = meshi->jbar;
    kbar = meshi->kbar;

    xplt = meshi->xplt_orig;
    yplt = meshi->yplt_orig;
    zplt = meshi->zplt_orig;

    if(
      xplt[0]<=xyz[0]&&xyz[0]<xplt[ibar]&&
      yplt[0]<=xyz[1]&&xyz[1]<yplt[jbar]&&
      zplt[0]<=xyz[2]&&xyz[2]<zplt[kbar]){
        return meshi;
    }
  }
  return NULL;
}

/* ------------------ on_mesh_boundary ------------------------ */

int on_mesh_boundary(float *xyz){
  int i;

  for(i = 0; i<nmeshes; i++){
    mesh *meshi;
    int ibar, jbar, kbar;
    float *xplt, *yplt, *zplt;

    meshi = meshinfo+i;

    ibar = meshi->ibar;
    jbar = meshi->jbar;
    kbar = meshi->kbar;

    xplt = meshi->xplt_orig;
    yplt = meshi->yplt_orig;
    zplt = meshi->zplt_orig;

    if(xyz[0]<xplt[0]-MESHEPS||xyz[0]>xplt[ibar]+MESHEPS)continue;
    if(xyz[1]<yplt[0]-MESHEPS||xyz[1]>yplt[jbar]+MESHEPS)continue;
    if(xyz[2]<zplt[0]-MESHEPS||xyz[2]>zplt[kbar]+MESHEPS)continue;

    // pt on xmin face
    
    if(ABS(xplt[0]-xyz[0])<=MESHEPS&&
      yplt[0]-MESHEPS<=xyz[1]&&xyz[1]<=yplt[jbar]+MESHEPS&&
      zplt[0]-MESHEPS<=xyz[2]&&xyz[2]<=zplt[kbar]+MESHEPS)return 1;

    // pt on xmax face
    
    if(ABS(xplt[ibar]-xyz[0])<=MESHEPS&&
      yplt[0]-MESHEPS<=xyz[1]&&xyz[1]<=yplt[jbar]+MESHEPS&&
      zplt[0]-MESHEPS<=xyz[2]&&xyz[2]<=zplt[kbar]+MESHEPS)return 1;

    // pt on ymin face
    
    if(ABS(yplt[0]-xyz[1])<=MESHEPS&&
      xplt[0]-MESHEPS<=xyz[0]&&xyz[0]<=xplt[ibar]+MESHEPS&&
      zplt[0]-MESHEPS<=xyz[2]&&xyz[2]<=zplt[kbar]+MESHEPS)return 1;

    // pt on ymax face
    
    if(ABS(yplt[jbar]-xyz[1])<=MESHEPS&&
      xplt[0]-MESHEPS<=xyz[0]&&xyz[0]<=xplt[ibar]+MESHEPS&&
      zplt[0]-MESHEPS<=xyz[2]&&xyz[2]<=zplt[kbar]+MESHEPS)return 1;

    // pt on zmin face
    
    if(ABS(zplt[0]-xyz[2])<=MESHEPS&&
      xplt[0]-MESHEPS<=xyz[0]&&xyz[0]<=xplt[ibar]+MESHEPS&&
      yplt[0]-MESHEPS<=xyz[1]&&xyz[1]<=yplt[jbar]+MESHEPS)return 1;

    // pt on zmax face
    
    if(ABS(zplt[kbar]-xyz[2])<=MESHEPS&&
      xplt[0]-MESHEPS<=xyz[0]&&xyz[0]<=xplt[ibar]+MESHEPS&&
      yplt[0]-MESHEPS<=xyz[1]&&xyz[1]<=yplt[jbar]+MESHEPS)return 1;
  }
  return 0;
}

/* ------------------ getmesh_nofail ------------------------ */

mesh *getmesh_nofail(float *xyz){
  int i;

  for(i=0;i<nmeshes;i++){
    mesh *meshi;
    int ibar, jbar, kbar;
    float *xplt, *yplt, *zplt;

    meshi = meshinfo+i;

    ibar = meshi->ibar;
    jbar = meshi->jbar;
    kbar = meshi->kbar;

    xplt = meshi->xplt_orig;
    yplt = meshi->yplt_orig;
    zplt = meshi->zplt_orig;

    if(
      xplt[0]<=xyz[0]&&xyz[0]<xplt[ibar]&&
      yplt[0]<=xyz[1]&&xyz[1]<yplt[jbar]&&
      zplt[0]<=xyz[2]&&xyz[2]<zplt[kbar]){
      return meshi;
    }
  }
  for(i=0;i<nmeshes;i++){
    mesh *meshi;
    int ibar, jbar, kbar;
    float *xplt, *yplt, *zplt;

    meshi = meshinfo+i;

    ibar = meshi->ibar;
    jbar = meshi->jbar;
    kbar = meshi->kbar;

    xplt = meshi->xplt_orig;
    yplt = meshi->yplt_orig;
    zplt = meshi->zplt_orig;

    if(
      xplt[0]<=xyz[0]+MESHEPS&&xyz[0]-MESHEPS<=xplt[ibar]&&
      yplt[0]<=xyz[1]+MESHEPS&&xyz[1]-MESHEPS<=yplt[jbar]&&
      zplt[0]<=xyz[2]+MESHEPS&&xyz[2]-MESHEPS<=zplt[kbar]){
      return meshi;
    }
  }
  return meshinfo;
}

/* ------------------ ExtractFrustum ------------------------ */

void ExtractFrustum(void){

/* code from:  http://www.crownandcutlass.com/features/technicaldetails/frustum.html */
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

int PointInFrustum( float x, float y, float z){
   if( frustum[0][0]*x + frustum[0][1]*y + frustum[0][2]*z + frustum[0][3] <= 0 )return 0;
   if( frustum[1][0]*x + frustum[1][1]*y + frustum[1][2]*z + frustum[1][3] <= 0 )return 0;
   if( frustum[2][0]*x + frustum[2][1]*y + frustum[2][2]*z + frustum[2][3] <= 0 )return 0;
   if( frustum[3][0]*x + frustum[3][1]*y + frustum[3][2]*z + frustum[3][3] <= 0 )return 0;
   if( frustum[4][0]*x + frustum[4][1]*y + frustum[4][2]*z + frustum[4][3] <= 0 )return 0;
   if( frustum[5][0]*x + frustum[5][1]*y + frustum[5][2]*z + frustum[5][3] <= 0 )return 0;
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
  vi[0]=-(mi[0]*v[0]+mi[4]*v[1]+ mi[8]*v[2])*vi[3];
  vi[1]=-(mi[1]*v[0]+mi[5]*v[1]+ mi[9]*v[2])*vi[3];
  vi[2]=-(mi[2]*v[0]+mi[6]*v[1]+mi[10]*v[2])*vi[3];
}

/* ------------------ compareisonodes ------------------------ */

int compare_volfacelistdata( const void *arg1, const void *arg2 ){
  volfacelistdata *vi, *vj;

  vi = *(volfacelistdata **)arg1;
  vj = *(volfacelistdata **)arg2;

  if(vi->dist2<vj->dist2)return 1;
  if(vi->dist2>vj->dist2)return -1;
  return 0;
}

/* ------------------ get_screen_mapping ------------------------ */

void get_screen_mapping(float *xyz0, float *screen_perm){
  GLdouble xyz[3];
  int viewport[4];
  GLdouble screen0[3];
  GLdouble screen[3];
  GLdouble modelview[16];
  GLdouble projection[16];
  int set;
  GLdouble maxvals[3];
  int min_index;

#define SETSCREEN(i1,i2,i3,dscreen)\
  if(i1==0)set=0;\
  if(set==0&&ABS((dscreen)[i1])>MAX( ABS((dscreen)[i2]) , ABS((dscreen)[i3]) ) ){\
    (dscreen)[i1]=SIGN((dscreen)[i1]);\
    (dscreen)[i2]=0.0;\
    (dscreen)[i3]=0.0;\
    set=1;\
  }

#define MAXABS3(x) (MAX(ABS((x)[0]),MAX(ABS((x)[1]),ABS((x)[2]))))

  glGetIntegerv(GL_VIEWPORT, viewport);
  glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
  glGetDoublev(GL_PROJECTION_MATRIX, projection);

  VECEQ3(xyz,xyz0);
  gluProject(xyz[0],xyz[1],xyz[2],modelview,projection,viewport,screen0,screen0+1,screen0+2);

  VECEQ3(xyz,xyz0);
  xyz[0]+=0.1;
  gluProject(xyz[0],xyz[1],xyz[2],modelview,projection,viewport,screen,screen+1,screen+2);
  VECDIFF3(screen_perm,screen,screen0);
  maxvals[0] = MAXABS3(screen_perm);

  VECEQ3(xyz,xyz0);
  xyz[1]+=0.1;
  gluProject(xyz[0],xyz[1],xyz[2],modelview,projection,viewport,screen,screen+1,screen+2);
  VECDIFF3(screen_perm+3,screen,screen0);
  maxvals[1] = MAXABS3(screen_perm+3);

  VECEQ3(xyz,xyz0);
  xyz[2]+=0.1;
  gluProject(xyz[0],xyz[1],xyz[2],modelview,projection,viewport,screen,screen+1,screen+2);
  VECDIFF3(screen_perm+6,screen,screen0);
  maxvals[2] = MAXABS3(screen_perm+6);

#ifdef _DEBUG
  PRINTF("%f %f %f\n",screen_perm[0],screen_perm[1],screen_perm[2]);
  PRINTF("%f %f %f\n",screen_perm[3],screen_perm[4],screen_perm[5]);
  PRINTF("%f %f %f\n",screen_perm[6],screen_perm[7],screen_perm[8]);
  PRINTF("\n");
#endif
 
  if(maxvals[0]<MIN(maxvals[1],maxvals[2])){
    min_index=0;
  }
  else if(maxvals[1]<MIN(maxvals[0],maxvals[2])){
    min_index=1;
  }
  else{
    min_index=2;
  }

  if(min_index==0){
    screen_perm[0]=0.0;
    screen_perm[1]=0.0;
    screen_perm[2]=0.0;
  }
  else{
    SETSCREEN(0,1,2,screen_perm);
    SETSCREEN(1,0,2,screen_perm);
    SETSCREEN(2,0,1,screen_perm);
  }

  if(min_index==1){
    screen_perm[3]=0.0;
    screen_perm[4]=0.0;
    screen_perm[5]=0.0;
  }
  else{
    SETSCREEN(0,1,2,screen_perm+3);
    SETSCREEN(1,0,2,screen_perm+3);
    SETSCREEN(2,0,1,screen_perm+3);
  }

  if(min_index==2){
    screen_perm[6]=0.0;
    screen_perm[7]=0.0;
    screen_perm[8]=0.0;
  }
  else{
    SETSCREEN(0,1,2,screen_perm+6);
    SETSCREEN(1,0,2,screen_perm+6);
    SETSCREEN(2,0,1,screen_perm+6);
  }
#ifdef _DEBUG
  PRINTF("%f %f %f\n",screen_perm[0],screen_perm[1],screen_perm[2]);
  PRINTF("%f %f %f\n",screen_perm[3],screen_perm[4],screen_perm[5]);
  PRINTF("%f %f %f\n",screen_perm[6],screen_perm[7],screen_perm[8]);
  PRINTF("\n");
#endif
}

  /* ------------------ getvolsmokedir ------------------------ */

void getvolsmokedir(float *mm){
    /*
      ( m0 m4 m8  m12 ) (x)    (0)
      ( m1 m5 m9  m13 ) (y)    (0)
      ( m2 m6 m10 m14 ) (z)  = (0)
      ( m3 m7 m11 m15 ) (1)    (1)

       ( m0 m4  m8 )      (m12)
   Q=  ( m1 m5  m9 )  u = (m13)
       ( m2 m6 m10 )      (m14)
      
       ( m0 m1  m2 )
 Q^T=  ( m4 m5  m6 )
       ( m8 m9 m10 )

           ( M_x  0    0  )
       M = ( 0   M_y   0  )
           ( 0    0   M_z )

      (Q   u) (M) (x)     (0)      
      (v^T 1) (1) (y)   = (1)
       
      m3=m7=m11=0, v^T=0, y=1   QMx+u=0 => x=-inv(M)Q^Tu

            ( m0 m1  m2 ) (m12)   ( m0*m12 + m1*m13 +  m2*m14 )/M_x
       x = -( m4 m5  m6 ) (m13) = ( m4*m12 + m5*m13 +  m6*m14 )/M_y
            ( m8 m9 m10 ) (m14)   ( m8*m12 + m9*m13 + m10*m14 )/M_z

    */
  int i,ii,j;
  float norm[3];
  float eyedir[3];
  float cosdir;
  float angles[7];

  volfacelistdata *vi;

  if(freeze_volsmoke==1)return;

  xyzeyeorig[0] = -DOT3(mm+0,mm+12)/mscale[0];
  xyzeyeorig[1] = -DOT3(mm+4,mm+12)/mscale[1];
  xyzeyeorig[2] = -DOT3(mm+8,mm+12)/mscale[2];
  
  for(j=0;j<nmeshes;j++){
    mesh *meshj;
    int *inside;
    int *drawsides;
    float x0, x1, yy0, yy1, z0, z1;
    float xcen, ycen, zcen;
    
    meshj = meshinfo + j;

    inside = &meshj->inside;
    drawsides = meshj->drawsides;

    x0 = meshj->x0;
    x1 = meshj->x1;
    yy0 = meshj->y0;
    yy1 = meshj->y1;
    z0 = meshj->z0;
    z1 = meshj->z1;
    xcen = meshj->xcen;
    ycen = meshj->ycen;
    zcen = meshj->zcen;

    *inside=0;
    if(
      xyzeyeorig[0]>x0&&xyzeyeorig[0]<x1&&
      xyzeyeorig[1]>yy0&&xyzeyeorig[1]<yy1&&
      xyzeyeorig[2]>z0&&xyzeyeorig[2]<z1
      ){
      for(i=-3;i<=3;i++){
        if(i==0)continue;
        drawsides[i+3]=1;
      }
      *inside=1;
      continue;
    }

    for(i=-3;i<=3;i++){
      if(i==0)continue;
      ii = i;
      if(i<0)ii=-i;
      norm[0]=0.0;
      norm[1]=0.0;
      norm[2]=0.0;
      switch(ii){
      case 1:
        if(i<0){
          norm[0]=-1.0;
          eyedir[0]=x0;
        }
        else{
          norm[0]=1.0;
          eyedir[0]=x1;
        }
        eyedir[1]=ycen;
        eyedir[2]=zcen;
        break;
      case 2:
        eyedir[0]=xcen;
        if(i<0){
          norm[1]=-1.0;
          eyedir[1]=yy0;
        }
        else{
          norm[1]=1.0;
          eyedir[1]=yy1;
        }
        eyedir[2]=zcen;
        break;
      case 3:
        eyedir[0]=xcen;
        eyedir[1]=ycen;
        if(i<0){
          norm[2]=-1.0;
          eyedir[2]=z0;
        }
        else{
          norm[2]=1.0;
          eyedir[2]=z1;
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
      }
      VECDIFF3(eyedir,xyzeyeorig,eyedir);
      normalize(eyedir,3);
      cosdir = CLAMP(DOT3(eyedir,norm),-1.0,1.0);
      cosdir=acos(cosdir)*RAD2DEG;
      if(cosdir<0.0)cosdir=-cosdir;
      angles[3+i]=cosdir;
    }
    for(i=-3;i<=3;i++){
      if(i==0)continue;
      if(angles[i+3]<90.0){
        drawsides[i+3]=1;
      }
      else{
        drawsides[i+3]=0;
      }
    }
  }

  // turn off drawing for mesh sides that are on the inside of a supermesh
  if(combine_meshes==1){
    for(i=0;i<nmeshes;i++){
      mesh *meshi;
      int *drawsides,*extsides;
      int jj;

      meshi = meshinfo + i;
      drawsides = meshi->drawsides;
      extsides = meshi->extsides;
      for(jj=0;jj<7;jj++){
        if(extsides[jj]==0){
          drawsides[jj]=0;
        }
      }
    }
    for(i=0;i<nsupermeshinfo;i++){
      supermesh *smesh;

      smesh = supermeshinfo + i;
      for(j=0;j<7;j++){
        smesh->drawsides[j]=0;
      }
      for(j=0;j<smesh->nmeshes;j++){
        mesh *meshj;
        int k;

        meshj = smesh->meshes[j];
        for(k=0;k<7;k++){
          if(meshj->extsides[k]==1&&meshj->drawsides[k]==1)smesh->drawsides[k]=1;
        }
      }
    }
  }

  vi = volfacelistinfo;
  nvolfacelistinfo=0;
  for(i=0;i<nmeshes;i++){
    mesh *meshi;
    int facemap[7]={12,6,0,0,3,9,15};
    volrenderdata *vr;
    int *drawsides;

    meshi = meshinfo + i;

    drawsides = meshi->drawsides;

    vr = &(meshi->volrenderinfo);
    if(vr->firedataptr==NULL&&vr->smokedataptr==NULL)continue;
    if(vr->loaded==0||vr->display==0)continue;
    for(j=-3;j<=3;j++){
      float dx, dy, dz;
      float *xyz;

      if(j==0)continue;
      if(drawsides[j+3]==0)continue;
      vi->facemesh=meshi;
      vi->iwall=j;
      xyz=meshi->face_centers+facemap[j+3];

      dx = xyz[0]-xyzeyeorig[0];
      dy = xyz[1]-xyzeyeorig[1];
      dz = xyz[2]-xyzeyeorig[2];
      vi->dist2=dx*dx+dy*dy+dz*dz;
      vi->xyz=xyz;
      vi++;
      nvolfacelistinfo++;
    }
  }
  if(nvolfacelistinfo>0){
    for(i=0;i<nvolfacelistinfo;i++){
      volfacelistinfoptrs[i]=volfacelistinfo + i;
    }
    qsort((volfacelistdata *)volfacelistinfoptrs,nvolfacelistinfo,sizeof(volfacelistdata *),compare_volfacelistdata);
  }
}

/* ------------------ getsmokedir ------------------------ */

void getsmokedir(float *mm){
    /*
      ( m0 m4 m8  m12 ) (x)    (0)
      ( m1 m5 m9  m13 ) (y)    (0)
      ( m2 m6 m10 m14 ) (z)  = (0)
      ( m3 m7 m11 m15 ) (1)    (1)

       ( m0 m4  m8 )      (m12)
   Q=  ( m1 m5  m9 )  u = (m13)
       ( m2 m6 m10 )      (m14)
      
      (Q   u) (x)     (0)      
      (v^T 1) (y)   = (1)
       
      m3=m7=m11=0, v^T=0, y=1   Qx+u=0 => x=-Q^Tu
    */
  int i,ii,j;
  mesh *meshj;
  float norm[3],scalednorm[3];
  float normdir[3];
  float absangle,cosangle,minangle;
  int iminangle;
  float dx, dy, dz;
  float factor;

  xyzeyeorig[0] = -DOT3(mm+0,mm+12)/mscale[0];
  xyzeyeorig[1] = -DOT3(mm+4,mm+12)/mscale[1];
  xyzeyeorig[2] = -DOT3(mm+8,mm+12)/mscale[2];
  
  for(j=0;j<nmeshes;j++){
    meshj = meshinfo + j;

    minangle=1000.0;
    iminangle=-10;
    meshj->dx=meshj->xplt_orig[1]-meshj->xplt_orig[0];
    meshj->dy=meshj->yplt_orig[1]-meshj->yplt_orig[0];
    meshj->dz=meshj->zplt_orig[1]-meshj->zplt_orig[0];
    meshj->dxy=meshj->dx*meshj->dx+meshj->dy*meshj->dy;
    meshj->dxy=sqrt(meshj->dxy)/2.0;
    meshj->dxz=meshj->dx*meshj->dx+meshj->dz*meshj->dz;
    meshj->dxz=sqrt(meshj->dxz)/2.0;
    meshj->dyz=meshj->dy*meshj->dy+meshj->dz*meshj->dz;
    meshj->dyz=sqrt(meshj->dyz)/2.0;

    meshj->dy/=meshj->dx;
    meshj->dz/=meshj->dx;
    meshj->dxy/=meshj->dx;
    meshj->dxz/=meshj->dx;
    meshj->dyz/=meshj->dx;
    meshj->dx=1.0;

    if(smokedrawtest2==1){
      meshj->norm[0]=1.0;
       meshj->norm[1]=0.0;
       meshj->norm[2]=0.0;
       meshj->smokedir=1;
       continue;
    }

    for(i=-9;i<=9;i++){
      if(i==0)continue;
      ii = i;
      if(i<0)ii=-i;
      norm[0]=0.0;
      norm[1]=0.0;
      norm[2]=0.0;
      switch(ii){
      case 1:
        if(i<0)norm[0]=-1.0;
        if(i>0)norm[0]=1.0;
        break;
      case 2:
        if(i<0)norm[1]=-1.0;
        if(i>0)norm[1]=1.0;
        break;
      case 3:
        if(i<0)norm[2]=-1.0;
        if(i>0)norm[2]=1.0;
        break;
      case 4:
        dx = meshj->xplt_orig[1]-meshj->xplt_orig[0];
        dy = meshj->yplt_orig[1]-meshj->yplt_orig[0];
        factor= dx*dx+dy*dy;
        if(factor==0.0){
          factor=1.0;
        }
        else{
          factor=1.0/sqrt(factor);
        }
        if(i<0){
          norm[0]=-dy*factor;
          norm[1]=-dx*factor;
        }
        else{
          norm[0]=dy*factor;
          norm[1]=dx*factor;
        }
        break;
      case 5:
        dx = meshj->xplt_orig[1]-meshj->xplt_orig[0];
        dy = meshj->yplt_orig[1]-meshj->yplt_orig[0];
        factor= dx*dx+dy*dy;
        if(factor==0.0){
          factor=1.0;
        }
        else{
          factor=1.0/sqrt(factor);
        }
        if(i<0){
          norm[0]= dy*factor;
          norm[1]=-dx*factor;
        }
        else{
          norm[0]=-dy*factor;
          norm[1]= dx*factor;
        }
        break;
      case 6:
        dy = meshj->yplt_orig[1]-meshj->yplt_orig[0];
        dz = meshj->zplt_orig[1]-meshj->zplt_orig[0];
        factor= dz*dz+dy*dy;
        if(factor==0.0){
          factor=1.0;
        }
        else{
          factor=1.0/sqrt(factor);
        }
        if(i<0){
          norm[1]=-dz*factor;
          norm[2]=-dy*factor;
        }
        else{
          norm[1]=dz*factor;
          norm[2]=dy*factor;
        }      
        break;
      case 7:
        dy = meshj->yplt_orig[1]-meshj->yplt_orig[0];
        dz = meshj->zplt_orig[1]-meshj->zplt_orig[0];
        factor= dz*dz+dy*dy;
        if(factor==0.0){
          factor=1.0;
        }
        else{
          factor=1.0/sqrt(factor);
        }
        if(i<0){
          norm[1]= dz*factor;
          norm[2]=-dy*factor;
        }
        else{
          norm[1]=-dz*factor;
          norm[2]= dy*factor;
        }
        break;
      case 8:
        dx = meshj->xplt_orig[1]-meshj->xplt_orig[0];
        dz = meshj->zplt_orig[1]-meshj->zplt_orig[0];
        factor= dz*dz+dx*dx;
        if(factor==0.0){
          factor=1.0;
        }
        else{
          factor=1.0/sqrt(factor);
        }
        if(i<0){
          norm[0]=-dz*factor;
          norm[2]=-dx*factor;
        }
        else{
          norm[0]=dz*factor;
          norm[2]=dx*factor;
        }      
        break;
      case 9:
        dx = meshj->xplt_orig[1]-meshj->xplt_orig[0];
        dz = meshj->zplt_orig[1]-meshj->zplt_orig[0];
        factor= dx*dx+dz*dz;
        if(factor==0.0){
          factor=1.0;
        }
        else{
          factor=1.0/sqrt(factor);
        }
        if(i<0){
          norm[0]= dz*factor;
          norm[2]=-dx*factor;
        }
        else{
          norm[0]=-dz*factor;
          norm[2]= dx*factor;
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
      }
      scalednorm[0]=norm[0]*mscale[0];
      scalednorm[1]=norm[1]*mscale[1];
      scalednorm[2]=norm[2]*mscale[2];

      normdir[0] = DOT3SKIP(mm,4,scalednorm,1);
      normdir[1] = DOT3SKIP(mm+1,4,scalednorm,1);
      normdir[2] = DOT3SKIP(mm+2,4,scalednorm,1);

      cosangle = CLAMP(normdir[2]/NORM3(normdir),-1.0,1.0);
      absangle=acos(cosangle)*RAD2DEG;
      if(absangle<0.0)absangle=-absangle;
      if(absangle<minangle){
        iminangle=i;
        minangle=absangle;
        meshj->norm[0]=norm[0];
        meshj->norm[1]=norm[1];
        meshj->norm[2]=norm[2];
      }
    }
    meshj->smokedir=iminangle;
#ifdef pp_CULL
    if(meshj->smokedir!=meshj->smokedir_old){
      meshj->smokedir_old=meshj->smokedir;
      update_initcullplane=1;
#ifdef _DEBUG
      PRINTF("mesh dir has changed\n");
#endif
    }
#endif
    if(demo_mode!=0){
      meshj->smokedir=1;
    }
  }
}

/* ------------------ get_interval ------------------------ */

int get_interval(float val, float *array, int n){
  int low, mid, high;

  if(val<array[0])return -1;
  if(val>array[n-1])return -1;

  low=0;
  high=n-1;
  while(high-low>1){
    mid=(low+high)/2;
    if(val>array[mid]){
      low=mid;
    }
    else{
      high=mid;
    }
  }
  ASSERT(low<n)
  return low;
}

/* ------------------ getnewpos ------------------------ */

void getnewpos(float *oldpos, float dx, float dy, float dz,float local_speed_factor){
  oldpos[0] += dx;
  oldpos[1] += dy;
  oldpos[2] += dz;
  from_glui_trainer=0;
}

/* ------------------ getblockage_distance ------------------------ */

float getblockage_distance(float x, float y, float z){
  int i;
  mesh *meshi;
  float *xplt, *yplt, *zplt;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  int ibar, jbar, kbar, nx, nxy;
  int ii, jj, kk;
  int ijknode,ijkcell;
  float view_height;
  char *iblank_cell;

  for(i=0;i<nmeshes;i++){
    meshi = meshinfo+i;

    iblank_cell = meshi->c_iblank_cell;

    ibar = meshi->ibar;
    jbar = meshi->jbar;
    kbar = meshi->kbar;
    nx = ibar+1;
    nxy = (ibar+1)*(jbar+1);

    xplt = meshi->xplt_orig;
    yplt = meshi->yplt_orig;
    zplt = meshi->zplt_orig;

    xmin = xplt[0];
    xmax = xplt[ibar];
    if(x<xmin||x>xmax)continue;

    ymin = yplt[0];
    ymax = yplt[jbar];
    if(y<ymin||y>ymax)continue;

    zmin = zplt[0];
    zmax = zplt[kbar];
    if(z<zmin||z>zmax)continue;

    ii = get_interval(x,xplt,ibar+1);
    jj = get_interval(y,yplt,jbar+1);
    kk = get_interval(z,zplt,kbar+1);
    if(ii!=-1&&jj!=-1&&kk!=-1){
      ijkcell=IJKCELL(ii,jj,kk);
      if(iblank_cell[ijkcell]==SOLID)return 0.0;
      ijknode=IJKNODE(ii,jj,kk);
      view_height = meshi->block_zdist[ijknode];
      if(view_height==0.0)return 0.0;
      view_height += (z-zplt[kk]);
      return view_height;
    }
  }
  return -1.0;
}

/* ------------------ init_blockdist  ------------------------ */

void init_blockage_distance(void){
  int ig,jg;
  float *b_zdist;
  mesh *meshi,*meshj;
  int minindex;
  int ibar,jbar,kbar, nx, nxy;
  int nnodes;
  int i,j,k;
  float zbottommin;
  float *xplt, *yplt, *zplt;
  float xx, yy, zz;
  int ijkm1cell, ijknode, ijkm1node;
  char *iblank_cell;

    for(ig=0;ig<nmeshes;ig++){
      meshi = meshinfo+ig;
      ibar = meshi->ibar;
      jbar = meshi->jbar;
      kbar = meshi->kbar;
      b_zdist=NULL;
      nnodes=(ibar+1)*(jbar+1)*(kbar+1);

      if(NewMemory((void **)&b_zdist,nnodes*sizeof(float))==0)return;
      for(i=0;i<nnodes;i++){
        b_zdist[i]=-1.0;
      }
      meshi->block_zdist=b_zdist;
      meshi->zdist_flag=0;
    }
    for(ig=0;ig<nmeshes;ig++){
      float dz;
      zbottommin=1000000000;
      for(jg=0;jg<nmeshes;jg++){
        meshj = meshinfo+jg;
        if(meshj->zdist_flag==1)continue;
        if(meshj->zplt[0]<zbottommin){
          zbottommin=meshj->zplt[0];
          minindex=jg;
        }
      }
      meshi=meshinfo+minindex;
      meshi->zdist_flag=1;

      xplt = meshi->xplt_orig;
      yplt = meshi->yplt_orig;
      zplt = meshi->zplt_orig;
      iblank_cell = meshi->c_iblank_cell;

      dz = zplt[1]-zplt[0];

      ibar = meshi->ibar;
      jbar = meshi->jbar;
      kbar = meshi->kbar;
      nx = ibar + 1;
      nxy = (ibar+1)*(jbar+1);

      b_zdist=meshi->block_zdist;
      nnodes=(ibar+1)*(jbar+1)*(kbar+1);


      // define first layer of b_zdist array
      //  if there is a mesh below first layer then add distance 
      //    to blockage in lower mesh

      k=0;
      zz = zplt[k];
      dz = zplt[1]-zplt[0];
      for(j=0;j<jbar;j++){
        yy = yplt[j];
        for(i=0;i<ibar;i++){
          float zdist;

          ijknode=IJKNODE(i,j,k);
          xx = xplt[i];
          zdist=getblockage_distance(xx,yy,zz-dz/2.0);
          if(zdist>0.0){
            b_zdist[ijknode]=zdist+dz/2.0;
          }
          else{
            b_zdist[ijknode]=0.0;
          }
        }
      }
      for(k=1;k<=kbar;k++){
        zz = zplt[k];
        dz = zplt[k]-zplt[k-1];
        for(j=0;j<jbar;j++){
          for(i=0;i<ibar;i++){
            ijknode=IJKNODE(i,j,k);
            ijkm1cell=IJKCELL(i,j,k-1);
            ijkm1node=IJKNODE(i,j,k-1);
            if(iblank_cell[ijkm1cell]==SOLID){
              b_zdist[ijknode]=0.0;
            }
            else{
              b_zdist[ijknode]=b_zdist[ijkm1node] + dz;
            }
          }
        }
      }
    }
}

/* ------------------ makeiblank_carve ------------------------ */

int makeiblank_carve(void){
  int i, j;
  int ibar, jbar, kbar;
  int ijksize;
  int nx, ny, nz, nxy;
  char *ib_embed;

  if(arg_iblank==0){
    if(autoterrain==1){
      use_iblank=0;
    }
    else{
      use_iblank=1;
    }
  }
  n_embedded_meshes=0;
  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi = meshinfo+i;
    meshi->c_iblank_embed=NULL;
  }
  if(nmeshes==1)return 0;


  for(i=0;i<nmeshes;i++){
    mesh *meshi;
    int n_embedded;

    meshi = meshinfo+i;
    ibar = meshi->ibar;
    jbar = meshi->jbar;
    kbar = meshi->kbar;
    nx = ibar + 1;
    ny = jbar + 1;
    nz = kbar + 1;
    nxy = nx*ny;
    ijksize=(ibar+1)*(jbar+1)*(kbar+1);

    meshi->c_iblank_embed=NULL;

    // check to see if there are any embedded meshes

    n_embedded=0;
    for(j=0;j<nmeshes;j++){
      mesh *meshj;

      if(i==j)continue;
      meshj = meshinfo + j;
      if(
        meshi->boxmin[0]<=meshj->boxmin[0]&&meshj->boxmax[0]<=meshi->boxmax[0]&&
        meshi->boxmin[1]<=meshj->boxmin[1]&&meshj->boxmax[1]<=meshi->boxmax[1]&&
        meshi->boxmin[2]<=meshj->boxmin[2]&&meshj->boxmax[2]<=meshi->boxmax[2]
      ){
        n_embedded++;
        n_embedded_meshes++;
      }
    }
    if(n_embedded==0)continue;

    ib_embed=NULL;
    if(use_iblank==1){
      if(NewMemory((void **)&ib_embed,ijksize*sizeof(char))==0)return 1;
    }
    meshi->c_iblank_embed=ib_embed;
    if(ib_embed==NULL)continue;
    for(j=0;j<ijksize;j++){
      ib_embed[j]=EMBED_NO;
    }
    for(j=0;j<nmeshes;j++){
      mesh *meshj;
      int i1, i2, jj1, j2, k1, k2;
      int ii, jj, kk;
      float *xplt, *yplt, *zplt;

      if(i==j)continue;
      meshj = meshinfo + j;
      // meshj is embedded inside meshi
      if(
        meshi->boxmin[0]>meshj->boxmin[0]||meshj->boxmax[0]>meshi->boxmax[0]||
        meshi->boxmin[1]>meshj->boxmin[1]||meshj->boxmax[1]>meshi->boxmax[1]||
        meshi->boxmin[2]>meshj->boxmin[2]||meshj->boxmax[2]>meshi->boxmax[2]
      )continue;

      xplt = meshi->xplt_orig;
      yplt = meshi->yplt_orig;
      zplt = meshi->zplt_orig;
      for(ii=0;ii<nx;ii++){
        if(xplt[ii]<=meshj->boxmin[0]&&meshj->boxmin[0]<xplt[ii+1]){
          i1=ii;
          break;
        }
      }
      for(ii=0;ii<nx;ii++){
        if(xplt[ii]<meshj->boxmax[0]&&meshj->boxmax[0]<=xplt[ii+1]){
          i2=ii;
          break;
        }
      }
      for(jj=0;jj<ny;jj++){
        if(yplt[jj]<=meshj->boxmin[1]&&meshj->boxmin[1]<yplt[jj+1]){
          jj1=jj;
          break;
        }
      }
      for(jj=0;jj<ny;jj++){
        if(yplt[jj]<meshj->boxmax[1]&&meshj->boxmax[1]<=yplt[jj+1]){
          j2=jj;
          break;
        }
      }
      for(kk=0;kk<nz;kk++){
        if(zplt[kk]<=meshj->boxmin[2]&&meshj->boxmin[2]<zplt[kk+1]){
          k1=kk;
          break;
        }
      }
      for(kk=0;kk<nz;kk++){
        if(zplt[kk]<meshj->boxmax[2]&&meshj->boxmax[2]<=zplt[kk+1]){
          k2=kk;
          break;
        }
      }

      for(kk=k1;kk<=k2;kk++){
        for(jj=jj1;jj<=j2;jj++){
          for(ii=i1;ii<=i2;ii++){
            ib_embed[IJKNODE(ii,jj,kk)]=EMBED_YES;
          }
        }
      }
    }
  }
  return 0;
}

/* ------------------ makeiblank ------------------------ */

int makeiblank(void){
  int ig;

  PRINTF("  initializing blanking array\n");
  if(use_iblank==0)return 0;
  for(ig=0;ig<nmeshes;ig++){
    mesh *meshi;
    int nx, ny, nxy;
    int ibar,jbar,kbar;
    float *fblank_cell=NULL;
    char *iblank_node=NULL,*iblank_cell=NULL,*c_iblank_x=NULL,*c_iblank_y=NULL,*c_iblank_z=NULL;
    int ii,ijksize;
    int i,j,k;

    meshi = meshinfo+ig;

    ibar = meshi->ibar;
    jbar = meshi->jbar;
    kbar = meshi->kbar;
    ijksize=(ibar+1)*(jbar+1)*(kbar+1);

    if(NewMemory((void **)&iblank_node,ijksize*sizeof(char))==0)return 1;
    if(NewMemory((void **)&iblank_cell,ibar*jbar*kbar*sizeof(char))==0)return 1;
    if(NewMemory((void **)&fblank_cell,ibar*jbar*kbar*sizeof(float))==0)return 1;
    if(NewMemory((void **)&c_iblank_x,ijksize*sizeof(char))==0)return 1;
    if(NewMemory((void **)&c_iblank_y,ijksize*sizeof(char))==0)return 1;
    if(NewMemory((void **)&c_iblank_z,ijksize*sizeof(char))==0)return 1;

    meshi->c_iblank_node=iblank_node;
    meshi->c_iblank_cell=iblank_cell;
    meshi->f_iblank_cell=fblank_cell;
    meshi->c_iblank_x=c_iblank_x;
    meshi->c_iblank_y=c_iblank_y;
    meshi->c_iblank_z=c_iblank_z;

    for(i=0;i<ibar*jbar*kbar;i++){
      iblank_cell[i]=GAS;
    }
    for(i=0;i<ijksize;i++){
      iblank_node[i]=GAS;
      c_iblank_x[i]=GAS;
      c_iblank_y[i]=GAS;
      c_iblank_z[i]=GAS;
    }

    nx = ibar+1;
    ny = jbar+1;
    nxy = nx*ny;

    for(ii=0;ii<meshi->nbptrs;ii++){
      blockagedata *bc;

      bc=meshi->blockageinfoptrs[ii];
      for(i=bc->ijk[IMIN];i<bc->ijk[IMAX];i++){
      for(j=bc->ijk[JMIN];j<bc->ijk[JMAX];j++){
      for(k=bc->ijk[KMIN];k<bc->ijk[KMAX];k++){
        iblank_cell[IJKCELL(i,j,k)]=SOLID;
      }
      }
      }
    }
    if(fblank_cell!=NULL){
      for(ii=0;ii<ibar*jbar*kbar;ii++){
        fblank_cell[ii]=iblank_cell[ii];
      }
    }
    for(i=0;i<ibar+1;i++){
      for(j=0;j<jbar+1;j++){
        for(k=0;k<kbar+1;k++){
          int test;

          test=0;
          if(i!=0&&j!=0&&k!=0)         test+=iblank_cell[IJKCELL(i-1,j-1,k-1)];
          if(i!=ibar&&j!=0&&k!=0)      test+=iblank_cell[IJKCELL(  i,j-1,k-1)];
          if(i!=0&&j!=jbar&&k!=0)      test+=iblank_cell[IJKCELL(i-1,  j,k-1)];
          if(i!=ibar&&j!=jbar&&k!=0)   test+=iblank_cell[IJKCELL(  i,  j,k-1)];
          if(i!=0&&j!=0&&k!=kbar)      test+=iblank_cell[IJKCELL(i-1,j-1,  k)];
          if(i!=ibar&&j!=0&&k!=kbar)   test+=iblank_cell[IJKCELL(  i,j-1,  k)];
          if(i!=0&&j!=jbar&&k!=kbar)   test+=iblank_cell[IJKCELL(i-1,  j,  k)];
          if(i!=ibar&&j!=jbar&&k!=kbar)test+=iblank_cell[IJKCELL(  i,  j,  k)];
          if(test==0)iblank_node[IJKNODE(i,j,k)]=0;
        }
      }
    }

    for(j=0;j<jbar;j++){
      for(k=0;k<kbar;k++){
        c_iblank_x[IJKNODE(0,j,k)]   =2*iblank_cell[IJKCELL(0,j,k)];
        for(i=1;i<ibar;i++){
          c_iblank_x[IJKNODE(i,j,k)]=iblank_cell[IJKCELL(i-1,j,k)]+iblank_cell[IJKCELL(i,j,k)];
        }
        c_iblank_x[IJKNODE(ibar,j,k)]=2*iblank_cell[IJKCELL(ibar-1,j,k)];
      }
    }
    for(i=0;i<ibar;i++){
      for(k=0;k<kbar;k++){
        c_iblank_y[IJKNODE(i,0,k)]=2*iblank_cell[IJKCELL(i,0,k)];
        for(j=1;j<jbar;j++){
          c_iblank_y[IJKNODE(i,j,k)]=iblank_cell[IJKCELL(i,j-1,k)]+iblank_cell[IJKCELL(i,j,k)];
        }
        c_iblank_y[IJKNODE(i,jbar,k)]=2*iblank_cell[IJKCELL(i,jbar-1,k)];
      }
    }

    for(i=0;i<ibar;i++){
      for(j=0;j<jbar;j++){
        c_iblank_z[IJKNODE(i,j,0)]=2*iblank_cell[IJKCELL(i,j,0)];
        for(k=1;k<kbar;k++){
          c_iblank_z[IJKNODE(i,j,k)]=iblank_cell[IJKCELL(i,j,k-1)]+iblank_cell[IJKCELL(i,j,k)];
        }
        c_iblank_z[IJKNODE(i,j,kbar)]=2*iblank_cell[IJKCELL(i,j,kbar-1)];
      }
    }
  }
  init_blockage_distance();
  PRINTF("  blanking array initialization completed\n");
  return 0;
}

/* ------------------ getmesh_in_smesh ------------------------ */

mesh *getmesh_in_smesh(mesh *mesh_guess, supermesh *smesh, float *xyz){
  int i;
  float *smin, *smax;

  smin = smesh->boxmin_scaled;
  smax = smesh->boxmax_scaled;

  if(xyz[0]<smin[0]||xyz[1]<smin[1]||xyz[2]<smin[2])return NULL;
  if(xyz[0]>smax[0]||xyz[1]>smax[1]||xyz[2]>smax[2])return NULL;
  for(i=-1;i<smesh->nmeshes;i++){
    mesh *meshi;
    float *bmin, *bmax;

    if(i==-1){
      if(mesh_guess==NULL)continue;
      meshi=mesh_guess;
    }
    else{
      meshi = smesh->meshes[i];
      if(meshi==mesh_guess)continue;
    }

    bmin = meshi->boxmin_scaled;
    bmax = meshi->boxmax_scaled;

    if(
      bmin[0]<=xyz[0]&&xyz[0]<=bmax[0]&&
      bmin[1]<=xyz[1]&&xyz[1]<=bmax[1]&&
      bmin[2]<=xyz[2]&&xyz[2]<=bmax[2]){
        return meshi;
    }
  }
  ASSERT(FFALSE);
  return NULL;
}


/* ------------------ init_clip ------------------------ */

void init_clip(void){
  clipdata *ci;

  clip_mode_last=-1;

  ci = &clipinfo;
  ci->clip_xmin=0;
  ci->clip_ymin=0;
  ci->clip_zmin=0;
  ci->clip_xmax=0;
  ci->clip_ymax=0;
  ci->clip_zmax=0;
  ci->xmin=0.0;
  ci->ymin=0.0;
  ci->zmin=0.0;
  ci->xmax=0.0;
  ci->ymax=0.0;
  ci->zmax=0.0;

  ci = &colorbar_clipinfo;
  ci->clip_xmin=1;
  ci->clip_ymin=1;
  ci->clip_zmin=1;
  ci->clip_xmax=1;
  ci->clip_ymax=1;
  ci->clip_zmax=1;
  ci->xmin=DENORMALIZE_X(2.0);
  ci->ymin=DENORMALIZE_X(2.0);
  ci->zmin=DENORMALIZE_Y(2.0);
  ci->xmax=DENORMALIZE_Y(2.0);
  ci->ymax=DENORMALIZE_Z(2.0);
  ci->zmax=DENORMALIZE_Z(2.0);

  clip_i=0;
  clip_j=0;
  clip_k=0;
  clip_I=0;
  clip_J=0;
  clip_K=0;

  stepclip_xmin=0,stepclip_ymin=0,stepclip_zmin=0;
  stepclip_xmax=0,stepclip_ymax=0,stepclip_zmax=0;
}


/* ------------------ volume_tetrahedron ------------------------ */

float volume_tetrahedron(float *v1, float *v2, float *v3, float *v4){
  float v2d[3], v3d[3], v4d[3], vcross[3];

  VECDIFF3(v2d,v2,v1);
  VECDIFF3(v3d,v3,v1);
  VECDIFF3(v4d,v4,v1);
  CROSS(vcross,v2d,v3d);
  return DOT3(v4d,vcross)/6.0;
}

/* ----------------------- initTetraClipInfo ----------------------------- */

void initTetraClipInfo(clipdata *ci,float *v1, float *v2, float *v3, float *v4){
  float v1d[3], v2d[3];
  GLdouble *clipvals;
  float vol;

  //     v4
  //     .  .
  //     .   v2
  //     .  /  \         v4           v4          v4             v3
  //     ./     \       /  \         /  \        /  \          /   \
  //    v1-------v3    v1---v3      v3---v2     v2---v1       v1---v2

  vol = volume_tetrahedron(v1,v2,v3,v4);
  
  clipvals = ci->clipvals;
  ci->option=TETRA_CLIPPLANES;

  VECDIFF3(v1d,v1,v3);
  VECDIFF3(v2d,v4,v3);
  CROSS(clipvals,v1d,v2d);
  if(vol>0.0){
    VEC3MA(clipvals,-1.0);
  }
  NORMALIZE3(clipvals);
  clipvals[3]=-DOT3(clipvals,v3);
  clipvals+=4;

  VECDIFF3(v1d,v3,v2);
  VECDIFF3(v2d,v4,v2);
  CROSS(clipvals,v1d,v2d);
  if(vol>0.0){
    VEC3MA(clipvals,-1.0);
  }
  NORMALIZE3(clipvals);
  clipvals[3]=-DOT3(clipvals,v2);
  clipvals+=4;

  VECDIFF3(v1d,v2,v1);
  VECDIFF3(v2d,v4,v1);
  CROSS(clipvals,v1d,v2d);
  if(vol>0.0){
    VEC3MA(clipvals,-1.0);
  }
  NORMALIZE3(clipvals);
  clipvals[3]=-DOT3(clipvals,v1);
  clipvals+=4;

  VECDIFF3(v1d,v3,v1);
  VECDIFF3(v2d,v2,v1);
  CROSS(clipvals,v1d,v2d);
  if(vol>0.0){
    VEC3MA(clipvals,-1.0);
  }
  NORMALIZE3(clipvals);
  clipvals[3]=-DOT3(clipvals,v1);
}

/* ----------------------- initBoxClipInfo ----------------------------- */

void initBoxClipInfo(clipdata *ci,float xmin, float xmax, float ymin, float ymax, float zmin, float zmax){
  ci->option=BOX_CLIPPLANES;
  ci->clip_xmin=1;
  ci->clip_xmax=1;
  ci->clip_ymin=1;
  ci->clip_ymax=1;
  ci->clip_zmin=1;
  ci->clip_zmax=1;
  ci->xmin=xmin;
  ci->xmax=xmax;
  ci->ymin=ymin;
  ci->ymax=ymax;
  ci->zmin=zmin;
  ci->zmax=zmax;
}


/* ----------------------- merge_max ----------------------------- */

float merge_max(int opti, float vali, int optj, float valj){
  if(opti==0&&optj==0)return vali;
  if(opti==1&&optj==0)return vali;
  if(opti==0&&optj==1)return valj;
  return MAX(vali,valj);
}

/* ----------------------- merge_min ----------------------------- */

float merge_min(int opti, float vali, int optj, float valj){
  if(opti==0&&optj==0)return vali;
  if(opti==1&&optj==0)return vali;
  if(opti==0&&optj==1)return valj;
  return MIN(vali,valj);
}

/* ----------------------- MergeClipPlanes ----------------------------- */

void MergeClipPlanes(clipdata *ci, clipdata *cj){
  ci->xmin = merge_max(ci->clip_xmin,ci->xmin,cj->clip_xmin,cj->xmin);
  ci->ymin = merge_max(ci->clip_ymin,ci->ymin,cj->clip_ymin,cj->ymin);
  ci->zmin = merge_max(ci->clip_zmin,ci->zmin,cj->clip_zmin,cj->zmin);
  ci->xmax = merge_min(ci->clip_xmax,ci->xmax,cj->clip_xmax,cj->xmax);
  ci->ymax = merge_min(ci->clip_ymax,ci->ymax,cj->clip_ymax,cj->ymax);
  ci->zmax = merge_min(ci->clip_zmax,ci->zmax,cj->clip_zmax,cj->zmax);

  ci->clip_xmin |= cj->clip_xmin;
  ci->clip_xmax |= cj->clip_xmax;
  ci->clip_ymin |= cj->clip_ymin;
  ci->clip_ymax |= cj->clip_ymax;
  ci->clip_zmin |= cj->clip_zmin;
  ci->clip_zmax |= cj->clip_zmax;
}

/* ----------------------- setClipPlanes ----------------------------- */

void setClipPlanes(clipdata *ci, int option){

  if(ci==NULL||option==CLIP_OFF){
    glDisable(GL_CLIP_PLANE0);
    glDisable(GL_CLIP_PLANE1);
    glDisable(GL_CLIP_PLANE2);
    glDisable(GL_CLIP_PLANE3);
    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);
    clipon=0;
    return;
  }
  clipon=1;

  if(ci->option==TETRA_CLIPPLANES){
    glClipPlane(GL_CLIP_PLANE0,ci->clipvals);
    glEnable(GL_CLIP_PLANE0);

    glClipPlane(GL_CLIP_PLANE1,ci->clipvals+4);
    glEnable(GL_CLIP_PLANE1);

    glClipPlane(GL_CLIP_PLANE2,ci->clipvals+8);
    glEnable(GL_CLIP_PLANE2);

    glClipPlane(GL_CLIP_PLANE3,ci->clipvals+12);
    glEnable(GL_CLIP_PLANE3);

    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);
    return;
  }

  if(ci->clip_xmin==1){
    GLdouble clipplane[4];

    clipplane[0]=1.0;
    clipplane[1]=0.0;
    clipplane[2]=0.0;
    if(option==CLIP_ON_DENORMAL)clipplane[3]=-ci->xmin;
    if(option==CLIP_ON)clipplane[3]=-NORMALIZE_X(ci->xmin);
    glClipPlane(GL_CLIP_PLANE0,clipplane);
    glEnable(GL_CLIP_PLANE0);
  }
  else{
    glDisable(GL_CLIP_PLANE0);
  }

  if(ci->clip_xmax==1){
    GLdouble clipplane[4];

    clipplane[0]=-1.0;
    clipplane[1]=0.0;
    clipplane[2]=0.0;
    if(option==CLIP_ON_DENORMAL)clipplane[3]=ci->xmax;
    if(option==CLIP_ON)clipplane[3]=NORMALIZE_X(ci->xmax);
    glClipPlane(GL_CLIP_PLANE3,clipplane);
    glEnable(GL_CLIP_PLANE3);
  }
  else{
    glDisable(GL_CLIP_PLANE3);
  }

  if(ci->clip_ymin==1){
    GLdouble clipplane[4];

    clipplane[0]=0.0;
    clipplane[1]=1.0;
    clipplane[2]=0.0;
    if(option==CLIP_ON_DENORMAL)clipplane[3]=-ci->ymin;
    if(option==CLIP_ON)clipplane[3]=-NORMALIZE_Y(ci->ymin);
    glClipPlane(GL_CLIP_PLANE1,clipplane);
    glEnable(GL_CLIP_PLANE1);
  }
  else{
    glDisable(GL_CLIP_PLANE1);
  }

  if(ci->clip_ymax==1){
    GLdouble clipplane[4];

    clipplane[0]=0.0;
    clipplane[1]=-1.0;
    clipplane[2]=0.0;
    if(option==CLIP_ON_DENORMAL)clipplane[3]=ci->ymax;
    if(option==CLIP_ON)clipplane[3]=NORMALIZE_Y(ci->ymax);
    glClipPlane(GL_CLIP_PLANE4,clipplane);
    glEnable(GL_CLIP_PLANE4);
  }
  else{
    glDisable(GL_CLIP_PLANE4);
  }

  if(ci->clip_zmin==1){
    GLdouble clipplane[4];

    clipplane[0]=0.0;
    clipplane[1]=0.0;
    clipplane[2]=1.0;
    if(option==CLIP_ON_DENORMAL)clipplane[3]=-ci->zmin;
    if(option==CLIP_ON)clipplane[3]=-NORMALIZE_Z(ci->zmin);
    glClipPlane(GL_CLIP_PLANE2,clipplane);
    glEnable(GL_CLIP_PLANE2);
  }
  else{
    glDisable(GL_CLIP_PLANE2);
  }

  if(ci->clip_zmax==1){
    GLdouble clipplane[4];

    clipplane[0]=0.0;
    clipplane[1]=0.0;
    clipplane[2]=-1.0;
    if(option==CLIP_ON_DENORMAL)clipplane[3]=ci->zmax;
    if(option==CLIP_ON)clipplane[3]=NORMALIZE_Z(ci->zmax);
    glClipPlane(GL_CLIP_PLANE5,clipplane);
    glEnable(GL_CLIP_PLANE5);
  }
  else{
    glDisable(GL_CLIP_PLANE5);
  }
}

