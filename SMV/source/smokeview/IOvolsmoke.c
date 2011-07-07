// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#ifdef pp_GPU
#include <GL/glew.h>
#endif
#include <stdio.h>  
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "flowfiles.h"
#include "MALLOC.h"
#include "smokeviewvars.h"
#include "interp.h"
#include "smv_endian.h"
#include "update.h"
#include "IOvolsmoke.h"

// svn revision character string
char IOvolsmoke_revision[]="$Revision$";

/* ----------------------- interp3d ----------------------------- */

#define INTERP1D(f0,f1,dx) (float)((f0) + ((f1)-(f0))*(dx))
void get_pt_smokecolor(float *smoke_tran, float **smoke_color, float dstep, float xyz[3], mesh *meshi, int *inobst, char *blank){
  int i, j, k;
  int ijk;
  float val000,val100,val010,val110;
  float val001,val101,val011,val111;
  float val00,val01,val10,val11;
  float val0, val1, val;
  int nx, ny, nz, nxy;
  float dx, dy, dz;
  float dxbar, dybar, dzbar;
  float *vv;
  int ijkcell;
  float *xplt, *yplt, *zplt;
  int ibar, jbar, kbar;
  float *smokedata, *firedata;
  float kfactor=8700.0;
  float soot_density, temperature;
  int index;

  smokedata = meshi->volrenderinfo.smokedata;
  firedata = meshi->volrenderinfo.firedata;

  xplt = meshi->xplt_cen;
  yplt = meshi->yplt_cen;
  zplt = meshi->zplt_cen;
  ibar = meshi->ibar;
  jbar = meshi->jbar;
  kbar = meshi->kbar;

  dxbar = xplt[1]-xplt[0];
  dybar = yplt[1]-yplt[0];
  dzbar = zplt[1]-zplt[0];

  nx = ibar + 1;
  ny = jbar + 1;
  nz = kbar + 1;
  nxy = nx*ny;

  GETINDEX(i,xyz[0],xplt[0],dxbar,ibar);
  GETINDEX(j,xyz[1],yplt[0],dybar,jbar);
  GETINDEX(k,xyz[2],zplt[0],dzbar,kbar);

  if(blank!=NULL){
    ijkcell=IJKCELL(i,j,k);
    if(blank[ijkcell]==0){
      *inobst=1;
      return;
    }
    else{
      *inobst=0;
    }
  }

  ijk = i + j*nx + k*nxy;

  dx = (xyz[0] - xplt[i])/dxbar;
  dx = CLAMP(dx,0.0,1.0);
  dy = (xyz[1] - yplt[j])/dybar;
  dy = CLAMP(dy,0.0,1.0);
  dz = (xyz[2] - zplt[k])/dzbar;
  dz = CLAMP(dz,0.0,1.0);

  if(firedata!=NULL){
    float dtemp;

    vv = firedata + ijk;
    val000 = vv[0]; // i,j,k
    val100 = vv[1]; // i+1,j,k

    vv += nx;
    val010 = vv[0]; // i,j+1,k
    val110 = vv[1]; // i+1,j+1,k

    vv += (nxy-nx);
    val001 = vv[0]; // i,j,k+1
    val101 = vv[1]; // i+1,j,k+1

    vv += nx;
    val011 = vv[0]; // i,j+1,k+1
    val111 = vv[1]; // i+1,j+1,k+1

    val00 = INTERP1D(val000,val100,dx);
    val10 = INTERP1D(val010,val110,dx);
    val01 = INTERP1D(val001,val101,dx);
    val11 = INTERP1D(val011,val111,dx);
     val0 = INTERP1D( val00, val10,dy);
     val1 = INTERP1D( val01, val11,dy);
    temperature = INTERP1D(  val0,  val1,dz);
    dtemp=(1200.0-20.0)/256;
    GETINDEX(index,temperature,20.0,dtemp,256);
    *smoke_color=rgb_smokecolormap+4*index;
  }
  if(smokedata!=NULL){
    vv = smokedata + ijk;
    val000 = vv[0]; // i,j,k
    val100 = vv[1]; // i+1,j,k

    vv += nx;
    val010 = vv[0]; // i,j+1,k
    val110 = vv[1]; // i+1,j+1,k

    vv += (nxy-nx);
    val001 = vv[0]; // i,j,k+1
    val101 = vv[1]; // i+1,j,k+1

    vv += nx;
    val011 = vv[0]; // i,j+1,k+1
    val111 = vv[1]; // i+1,j+1,k+1

    val00 = INTERP1D(val000,val100,dx);
    val10 = INTERP1D(val010,val110,dx);
    val01 = INTERP1D(val001,val101,dx);
    val11 = INTERP1D(val011,val111,dx);
     val0 = INTERP1D( val00, val10,dy);
     val1 = INTERP1D( val01, val11,dy);
     soot_density = INTERP1D(  val0,  val1,dz);
     if(firedata!=NULL&&index>128)soot_density*=5.0;
    *smoke_tran = exp(-kfactor*soot_density*dstep);
  }
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
    vr->display=0;
    vr->timeslist=NULL;
    vr->times_defined=0;
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

    if(STRCMP(shortlabel,"temp")==0){  
      vr->fire=slicei;
     continue;
    }
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
    vr->nframes=0;
    if(vr->smoke!=NULL){
      int nx, ny, nz, j;

      nvolrenderinfo++;
      vr->nframes=get_volsmoke_nframes(vr);
      nx = ijkbarmax+1;
      ny = ijkbarmax+1;
      nz = ijkbarmax+1;
      NewMemory((void **)&vr->smokecolor_yz0,4*ny*nz*sizeof(float));
      NewMemory((void **)&vr->smokecolor_yz1,4*ny*nz*sizeof(float));
      NewMemory((void **)&vr->smokecolor_xz0,4*nx*nz*sizeof(float));
      NewMemory((void **)&vr->smokecolor_xz1,4*nx*nz*sizeof(float));
      NewMemory((void **)&vr->smokecolor_xy0,4*nx*ny*sizeof(float));
      NewMemory((void **)&vr->smokecolor_xy1,4*nx*ny*sizeof(float));
      if(vr->nframes>0){
        NewMemory((void **)&vr->times,vr->nframes*sizeof(float));
        NewMemory((void **)&vr->firedataptrs,vr->nframes*sizeof(float *));
        NewMemory((void **)&vr->smokedataptrs,vr->nframes*sizeof(float *));
        for(j=0;j<vr->nframes;j++){
          vr->firedataptrs[j]=NULL;
          vr->smokedataptrs[j]=NULL;
        }
      }
    }
    else{
      vr->smokecolor_yz0=NULL;
      vr->smokecolor_yz1=NULL;
      vr->smokecolor_xz0=NULL;
      vr->smokecolor_xz1=NULL;
      vr->smokecolor_xy0=NULL;
      vr->smokecolor_xy1=NULL;
    }
  }
  if(nvolrenderinfo>0){
    NewMemory((void **)&volfacelistinfo,6*nmeshes*sizeof(volfacelistdata));
    NewMemory((void **)&volfacelistinfoptrs,6*nmeshes*sizeof(volfacelistdata *));
  }
}


//  glUniform3f(GPUvol_eyepos,xyzeyeorig[0],xyzeyeorig[1],xyzeyeorig[2]);
//  glUniform1i(GPUvol_inside,meshi->inside);
//  glUniform1f(GPUvol_xyzmaxdiff,xyzmaxdiff);
//  glUniform3f(GPUvol_boxmin,meshi->x0,meshi->y0,meshi->z0);
//  glUniform3f(GPUvol_boxmax,meshi->x1,meshi->y1,meshi->z1);

/* ------------------ get_cum_smokecolor ------------------------ */

void get_cum_smokecolor(float *cum_smokecolor, float *xyzvert, float dstep, mesh *meshi, int iwall){
  float t_intersect, t_intersect_min=FLT_MAX, *boxmin, *boxmax;
  float tmin, tmax;
  int i;
  int nsteps;
  float dist, dx, dy, dz;
  float distseg, dxseg, dyseg, dzseg;
  float xyz[3];
  float sootdensum;
  float opacity;
  float *vert_beg, *vert_end;
  int iwall_min=0;
  float xyzvals[3];
  int isteps;
  char *blank;
  float pt_smoketran, *pt_smokecolor;
  float cum_tran,tauhat,alphahat;

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
  if(distseg<0.001){
    cum_smokecolor[0]=0.0;
    cum_smokecolor[1]=0.0;
    cum_smokecolor[2]=0.0;
    cum_smokecolor[3]=0.0;
    return;
  }

  nsteps = distseg/dstep;
  if(nsteps<1){
    nsteps=1;
  }
  dstep=distseg/(float)nsteps;
  dstep*=xyzmaxdiff;
  sootdensum=0.0;
  isteps=0;
  if(block_volsmoke==1){
    blank=meshi->c_iblank_cell;
  }
  else{
    blank=NULL;
  }
  cum_tran = 1.0;
  cum_smokecolor[0]=0.0;
  cum_smokecolor[1]=0.0;
  cum_smokecolor[2]=0.0;
  cum_smokecolor[3]=0.0;
  tauhat=1.0;
  alphahat=0.0;
  for(i=0;i<nsteps;i++){
    float sootden, factor;
    int icell, jcell, kcell;
    int inobst;
    float alphai;

    factor = (0.5 + (float)i)/(float)nsteps;

    xyz[0] = (1.0-factor)*vert_beg[0] + factor*vert_end[0];
    xyz[1] = (1.0-factor)*vert_beg[1] + factor*vert_end[1];
    xyz[2] = (1.0-factor)*vert_beg[2] + factor*vert_end[2];

    get_pt_smokecolor(&pt_smoketran,&pt_smokecolor, dstep,xyz, meshi, &inobst, blank);
    if(blank!=NULL&&inobst==1)break;

    alphai = 1.0 - pt_smoketran;
    alphahat +=  alphai*tauhat;

    cum_smokecolor[0] += alphai*pt_smokecolor[0]*tauhat;
    cum_smokecolor[1] += alphai*pt_smokecolor[1]*tauhat;
    cum_smokecolor[2] += alphai*pt_smokecolor[2]*tauhat;
    tauhat *= pt_smoketran;
  }
  if(alphahat>0.0){
    cum_smokecolor[0]/=alphahat;
    cum_smokecolor[1]/=alphahat;
    cum_smokecolor[2]/=alphahat;
    cum_smokecolor[3]=alphahat;
  }
  else{
    cum_smokecolor[0]=0.0;
    cum_smokecolor[1]=0.0;
    cum_smokecolor[2]=0.0;
    cum_smokecolor[3]=0.0;
  }
}

/* ------------------ compute_all_smokecolors ------------------------ */

void compute_all_smokecolors(void){
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
    float *smokecolor;

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
            smokecolor=vr->smokecolor_yz0;
            xyz[0] = meshi->x0;
          }
          else{
            smokecolor=vr->smokecolor_yz1;
            xyz[0] = meshi->x1;
          }
          if(vr->firedata==NULL||vr->smokedata==NULL){
            for(i=0;i<=jbar;i++){
              for(j=0;j<=kbar;j++){
                smokecolor[0]=0.0;
                smokecolor[1]=0.0;
                smokecolor[2]=0.0;
                smokecolor[3]=0.0;
                smokecolor+=4;
              }
            }
          }
          else{
            for(i=0;i<=jbar;i++){
              xyz[1] = y[i];
              for(j=0;j<=kbar;j++){
                xyz[2] = z[j];
                get_cum_smokecolor(smokecolor,xyz,dstep,meshi,iwall);
                smokecolor+=4;
              }
            }
          }
          break;
        case 2:
        case -2:
          if(iwall<0){
            smokecolor=vr->smokecolor_xz0;
            xyz[1] = meshi->y0;
          }
          else{
            smokecolor=vr->smokecolor_xz1;
            xyz[1] = meshi->y1;
          }
          if(vr->firedata==NULL||vr->smokedata==NULL){
            for(i=0;i<=ibar;i++){
              for(j=0;j<=kbar;j++){
                smokecolor[0]=0.0;
                smokecolor[1]=0.0;
                smokecolor[2]=0.0;
                smokecolor[3]=0.0;
                smokecolor+=4;
              }
            }
          }
          else{
            for(i=0;i<=ibar;i++){
              xyz[0] = x[i];
              for(j=0;j<=kbar;j++){
                xyz[2] = z[j];
                get_cum_smokecolor(smokecolor,xyz,dstep,meshi,iwall);
                smokecolor+=4;
              }
            }
          }
          break;
        case 3:
        case -3:
          if(iwall<0){
            smokecolor=vr->smokecolor_xy0;
            xyz[2]=meshi->z0;
          }
          else{
            smokecolor=vr->smokecolor_xy1;
            xyz[2]=meshi->z1;
          }
          if(vr->firedata==NULL||vr->smokedata==NULL){
            for(i=0;i<=ibar;i++){
              for(j=0;j<=jbar;j++){
                smokecolor[0]=0.0;
                smokecolor[1]=0.0;
                smokecolor[2]=0.0;
                smokecolor[3]=0.0;
                smokecolor+=4;
              }
            }
          }
          else{
            for(i=0;i<=ibar;i++){
              xyz[0] = x[i];
              for(j=0;j<=jbar;j++){
                xyz[1] = y[j];
                get_cum_smokecolor(smokecolor,xyz,dstep,meshi,iwall);
                smokecolor+=4;
              }
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
    float *smokecolor;

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
          smokecolor = vr->smokecolor_yz0;
        }
        else{
          xx=meshi->x1;
          smokecolor = vr->smokecolor_yz1;
        }
        n00 = 0;
        n01 = 4;
        n10 = 4*(kbar+1);
        n11 = 4*(1 + kbar+1);
        for(i=0;i<jbar;i++){
          y[0] = yplt[i];
          y[1] = yplt[i+1];
          for(j=0;j<kbar;j++){
            z[0] = zplt[j];
            z[1] = zplt[j+1];

            if(meshi->inside==0&&iwall>0||meshi->inside!=0&&iwall<0){
              glColor4fv(smokecolor+n00);
              glVertex3f(xx,y[0],z[0]);
              glColor4fv(smokecolor+n10);
              glVertex3f(xx,y[1],z[0]);
              glColor4fv(smokecolor+n11);
              glVertex3f(xx,y[1],z[1]);

              glColor4fv(smokecolor+n00);
              glVertex3f(xx,y[0],z[0]);
              glColor4fv(smokecolor+n11);
              glVertex3f(xx,y[1],z[1]);
              glColor4fv(smokecolor+n01);
              glVertex3f(xx,y[0],z[1]);
            }
            else{
              glColor4fv(smokecolor+n00);
              glVertex3f(xx,y[0],z[0]);
              glColor4fv(smokecolor+n11);
              glVertex3f(xx,y[1],z[1]);
              glColor4fv(smokecolor+n10);
              glVertex3f(xx,y[1],z[0]);

              glColor4fv(smokecolor+n00);
              glVertex3f(xx,y[0],z[0]);
              glColor4fv(smokecolor+n01);
              glVertex3f(xx,y[0],z[1]);
              glColor4fv(smokecolor+n11);
              glVertex3f(xx,y[1],z[1]);
            }
            smokecolor+=4;
          }
         smokecolor+=4;
        }
        break;
      case 2:
      case -2:
        n00 = 0;
        n01 = 4;
        n10 = 4*(kbar+1);
        n11 = 4*(1 + kbar+1);
        if(iwall<0){
          smokecolor = vr->smokecolor_xz0;
          yy=meshi->y0;
        }
        else{
          smokecolor = vr->smokecolor_xz1;
          yy=meshi->y1;
        }
        for(i=0;i<ibar;i++){
          x[0] = xplt[i];
          x[1] = xplt[i+1];
          for(j=0;j<kbar;j++){
            z[0] = zplt[j];
            z[1] = zplt[j+1];
            if(meshi->inside==0&&iwall>0||meshi->inside!=0&&iwall<0){
              glColor4fv(smokecolor+n00);
              glVertex3f(x[0],yy,z[0]);
              glColor4fv(smokecolor+n11);
              glVertex3f(x[1],yy,z[1]);
              glColor4fv(smokecolor+n10);
              glVertex3f(x[1],yy,z[0]);

              glColor4fv(smokecolor+n00);
              glVertex3f(x[0],yy,z[0]);
              glColor4fv(smokecolor+n01);
              glVertex3f(x[0],yy,z[1]);
              glColor4fv(smokecolor+n11);
              glVertex3f(x[1],yy,z[1]);
            }
            else{
              glColor4fv(smokecolor+n00);
              glVertex3f(x[0],yy,z[0]);
              glColor4fv(smokecolor+n10);
              glVertex3f(x[1],yy,z[0]);
              glColor4fv(smokecolor+n11);
              glVertex3f(x[1],yy,z[1]);

              glColor4fv(smokecolor+n00);
              glVertex3f(x[0],yy,z[0]);
              glColor4fv(smokecolor+n11);
              glVertex3f(x[1],yy,z[1]);
              glColor4fv(smokecolor+n01);
              glVertex3f(x[0],yy,z[1]);
            }
            smokecolor+=4;
          }
          smokecolor+=4;
        }
        break;
      case 3:
      case -3:
        n00 = 0;
        n01 = 4;
        n10 = 4*(jbar+1);
        n11 = 4*(1 + jbar+1);
       if(iwall<0){
          smokecolor = vr->smokecolor_xy0;
          zz=meshi->z0;
        }
        else{
          smokecolor = vr->smokecolor_xy1;
          zz=meshi->z1;
        }
        for(i=0;i<ibar;i++){
          x[0] = xplt[i];
          x[1] = xplt[i+1];
          for(j=0;j<jbar;j++){
            y[0] = yplt[j];
            y[1] = yplt[j+1];
            if(meshi->inside==0&&iwall>0||meshi->inside!=0&&iwall<0){
              glColor4fv(smokecolor+n00);
              glVertex3f(x[0],y[0],zz);
              glColor4fv(smokecolor+n10);
              glVertex3f(x[1],y[0],zz);
              glColor4fv(smokecolor+n11);
              glVertex3f(x[1],y[1],zz);

              glColor4fv(smokecolor+n00);
              glVertex3f(x[0],y[0],zz);
              glColor4fv(smokecolor+n11);
              glVertex3f(x[1],y[1],zz);
              glColor4fv(smokecolor+n01);
              glVertex3f(x[0],y[1],zz);
            }
            else{
              glColor4fv(smokecolor+n00);
              glVertex3f(x[0],y[0],zz);
              glColor4fv(smokecolor+n11);
              glVertex3f(x[1],y[1],zz);
              glColor4fv(smokecolor+n10);
              glVertex3f(x[1],y[0],zz);

              glColor4fv(smokecolor+n00);
              glVertex3f(x[0],y[0],zz);
              glColor4fv(smokecolor+n01);
              glVertex3f(x[0],y[1],zz);
              glColor4fv(smokecolor+n11);
              glVertex3f(x[1],y[1],zz);
            }
            smokecolor+=4;
          }
          smokecolor+=4;
        }
        break;
    }
    glEnd();
  }
  if(use_transparency_data==1)transparentoff();
}

/* ------------------ drawsmoke3dGPUVOL ------------------------ */

#ifdef pp_GPU
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

#define HEADER_SIZE 4
#define TRAILER_SIZE 4
#define FORTSLICEREAD(var,size) fseek(SLICEFILE,HEADER_SIZE,SEEK_CUR);\
                           returncode=fread(var,4,size,SLICEFILE);\
                           if(endianswitch==1)endian_switch(var,size);\
                           fseek(SLICEFILE,TRAILER_SIZE,SEEK_CUR)

/* ------------------ get_volsmoke_sizes ------------------------ */

int get_volsmoke_nframes(volrenderdata *vr){
	slice *fireslice, *smokeslice;
  FILE *SLICEFILE;
  int framesize,skip,returncode;
  float time, *sliceframe_data;
  int endianswitch=0;
  int nf;
  int nframes,filesize;

  smokeslice=vr->smoke;
  fireslice=vr->fire;
  framesize = smokeslice->nslicei*smokeslice->nslicej*smokeslice->nslicek;
  framesize *= 4; // convert to bytes
  framesize += HEADER_SIZE + TRAILER_SIZE;

  skip =           (HEADER_SIZE+30        +TRAILER_SIZE); // long label
  skip +=          (HEADER_SIZE+30        +TRAILER_SIZE); // short label
  skip +=          (HEADER_SIZE+30        +TRAILER_SIZE); // unit label
  skip +=          (HEADER_SIZE+24        +TRAILER_SIZE); // is1, is2, js1, js2, ks1, ks2

  // nframes = (totalsize - skip)/(12 + framesize);

  nframes=0;
  filesize=getfilesize(smokeslice->file);
  if(filesize>0){
    nframes = (filesize-skip)/(12 + framesize);
  }
  return nframes;
}

/* ------------------ get_volsmoke_frame_time ------------------------ */

float get_volsmoke_frame_time(volrenderdata *vr, int framenum){
	slice *fireslice, *smokeslice;
  FILE *SLICEFILE;
  int framesize,skip,returncode;
  float time=0.0, *sliceframe_data;
  int endianswitch=0;
  char *meshlabel;

  if(framenum<0||framenum>=vr->nframes)return time;
  smokeslice=vr->smoke;
  framesize = smokeslice->nslicei*smokeslice->nslicej*smokeslice->nslicek;

  skip =           (HEADER_SIZE+30        +TRAILER_SIZE); // long label
  skip +=          (HEADER_SIZE+30        +TRAILER_SIZE); // short label
  skip +=          (HEADER_SIZE+30        +TRAILER_SIZE); // unit label
  skip +=          (HEADER_SIZE+24        +TRAILER_SIZE); // is1, is2, js1, js2, ks1, ks2
  skip += framenum*(HEADER_SIZE +4        +TRAILER_SIZE); // framenum time's
  skip += framenum*(HEADER_SIZE +4*framesize+TRAILER_SIZE); // framenum slice data's

  SLICEFILE=fopen(smokeslice->file,"rb");
  if(SLICEFILE==NULL)return time;

  returncode=fseek(SLICEFILE,skip,SEEK_SET); // skip from beginning of file

  FORTSLICEREAD(&time,1);
  fclose(SLICEFILE);
  return time;
}

/* ------------------ get_volsmoke_frame_time ------------------------ */

void get_volsmoke_all_times(volrenderdata *vr){
  int i;

  for(i=0;i<vr->nframes;i++){
    vr->times[i]=get_volsmoke_frame_time(vr,i);
  }
}

/* ------------------ read_volsmoke_frame ------------------------ */

void read_volsmoke_frame(volrenderdata *vr, int framenum, int *first){
	slice *fireslice, *smokeslice;
  FILE *SLICEFILE;
  int framesize,skip,returncode;
  float time, *sliceframe_data;
  int endianswitch=0;
  char *meshlabel;

  if(framenum<0||framenum>=vr->nframes)return;
  meshlabel=vr->rendermesh->label;
  smokeslice=vr->smoke;
  fireslice=vr->fire;
  framesize = smokeslice->nslicei*smokeslice->nslicej*smokeslice->nslicek;

  skip =           (HEADER_SIZE+30        +TRAILER_SIZE); // long label
  skip +=          (HEADER_SIZE+30        +TRAILER_SIZE); // short label
  skip +=          (HEADER_SIZE+30        +TRAILER_SIZE); // unit label
  skip +=          (HEADER_SIZE+24        +TRAILER_SIZE); // is1, is2, js1, js2, ks1, ks2
  skip += framenum*(HEADER_SIZE +4        +TRAILER_SIZE); // framenum time's
  skip += framenum*(HEADER_SIZE +4*framesize+TRAILER_SIZE); // framenum slice data's

  SLICEFILE=fopen(smokeslice->file,"rb");
  if(SLICEFILE==NULL)return;

  returncode=fseek(SLICEFILE,skip,SEEK_SET); // skip from beginning of file

  FORTSLICEREAD(&time,1);
  if(*first==1){
    *first=0;
    printf("time=%.2f %s: ",time,meshlabel);
  }
  else{
    if(time>=10.0)printf(" ");
    if(time>=100.0)printf(" ");
    if(time>=1000.0)printf(" ");
    printf("          %s: ",meshlabel);
  }

  vr->times[framenum]=time;
  NewMemory((void **)&sliceframe_data,framesize*sizeof(float));
  vr->smokedataptrs[framenum]=sliceframe_data;
  FORTSLICEREAD(sliceframe_data,framesize);
  printf("smoke");
  fclose(SLICEFILE);

  if(fireslice!=NULL){
    SLICEFILE=fopen(fireslice->file,"rb");
    if(SLICEFILE!=NULL){
      returncode=fseek(SLICEFILE,skip,SEEK_SET); // skip from beginning of file

      FORTSLICEREAD(&time,1);
      vr->times[framenum]=time;
      NewMemory((void **)&sliceframe_data,framesize*sizeof(float));
      vr->firedataptrs[framenum]=sliceframe_data;
      FORTSLICEREAD(sliceframe_data,framesize);
      printf(", fire");
      fclose(SLICEFILE);
    }
  }
  printf("\n");
}

/* ------------------ unload_volsmoke_allframes ------------------------ */

void unload_volsmoke_allframes(volrenderdata *vr){
  int i;

  printf("Unloading smoke %s - ",vr->rendermesh->label);
  for(i=0;i<vr->nframes;i++){
    FREEMEMORY(vr->firedataptrs[i]);
    FREEMEMORY(vr->smokedataptrs[i]);
  }
  vr->loaded=0;
  vr->display=0;
  plotstate = getplotstate(DYNAMIC_PLOTS);
  updatetimes();
  printf("completed\n");
}

/* ------------------ read_volsmoke_allframes ------------------------ */

void read_volsmoke_allframes(volrenderdata *vr){
  int nframes;
  int i;
  int first=1;

  nframes = vr->nframes;
  for(i=0;i<nframes;i++){
    read_volsmoke_frame(vr, i, &first);
  }
  vr->smokedata = vr->smokedataptrs[0];  //*** hack
  vr->firedata = vr->firedataptrs[0];
  vr->loaded=1;
  vr->display=1;
  plotstate=getplotstate(DYNAMIC_PLOTS);
  stept=1;
  updatetimes();
}

/* ------------------ read_volsmoke_frame_allmeshes ------------------------ */

void read_volsmoke_frame_allmeshes(int framenum){
  int i;
  float time_old;
  int ntimes_old;
  int first=1;

  for(i=0;i<nmeshes;i++){
    mesh *meshi;
    volrenderdata *vr;

    meshi = meshinfo + i;
    vr = &meshi->volrenderinfo;
    if(vr->fire==NULL||vr->smoke==NULL)continue;
    read_volsmoke_frame(vr,framenum,&first);
  }
  for(i=0;i<nmeshes;i++){
    mesh *meshi;
    volrenderdata *vr;

    meshi = meshinfo + i;
    vr = &meshi->volrenderinfo;
    if(vr->fire==NULL||vr->smoke==NULL)continue;
    if(framenum==0){
      vr->smokedata = vr->smokedataptrs[0];  //*** hack
      vr->firedata = vr->firedataptrs[0];
    }
    vr->loaded=1;
    vr->display=1;
  }
  plotstate=getplotstate(DYNAMIC_PLOTS);
  stept=1;
  updatetimes();
}

/* ------------------ read_volsmoke_allframes_allmeshes ------------------------ */

void read_volsmoke_allframes_allmeshes(void){
    int nframes=0;
    int i;

    for(i=0;i<nmeshes;i++){
      mesh *meshi;
      volrenderdata *vr;

      meshi = meshinfo + i;
      vr = &meshi->volrenderinfo;
      if(vr->fire==NULL||vr->smoke==NULL)continue;
      if(vr->nframes>0){
        nframes=vr->nframes;
        break;
      }
    }
    for(i=0;i<nframes;i++){
      read_volsmoke_frame_allmeshes(i);
    }
}
