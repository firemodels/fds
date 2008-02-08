// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#ifdef pp_LIGHT
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "zlib.h"
#include "egz_stdio.h"
#include "svzip.h"
#include "MALLOC.h"

// svn revision character string
char lightsmoke_revision[]="$Revision$";

float getlight_dist(float light_dist,float *xyz, float x, float y, float z);
float getalpha(mesh *smoke_mesh, float *xyz2);

#define IJKNODE(i,j,k) ((i)+(j)*nx+(k)*nxy)

/* ------------------ init_lightfield ------------------------ */

void init_lightfield(void){
  NewMemory((void **)&light_q_polar,NRAD*NTHETA*NPSI);
}

/* ------------------ update_lightfield ------------------------ */

void update_lightfield(smoke3d *smoke3di, unsigned char *lightingbuffer){
  int ilight;
  int nlight_q_rect;
  int i, j, k;

  if(smoke3di->smoke_mesh==NULL)return;
  nlight_q_rect = smoke3di->nx*smoke3di->ny*smoke3di->nz;
  if(smoke3di->light_q_rect==NULL){
    if(nlight_q_rect>0){
      NewMemory((void **)&smoke3di->light_q_rect,nlight_q_rect*sizeof(float));
    }
    else{
      return;
    }
  }
  for(i=0;i<nlight_q_rect;i++){
    smoke3di->light_q_rect[i]=0.0;
  }

  // accumulate hrr for each light

  for(ilight=0;ilight<nlightinfo;ilight++){
    lightdata *lighti;

    lighti = lightinfo + ilight;
    switch (lighti->type){
      float *xyz1, *xyz2;
      float dx,dy,dz,length;
      int npoint, nx, ny, nz;
      float dxx, dyy, dzz;
      float xyz[3];

      case 0:      // point
        set_lightfield(smoke3di,lighti->xyz1,lighti->q);
        break;
      case 1:      // line
        xyz1 = lighti->xyz1;
        xyz2 = lighti->xyz2;
        dx = xyz1[0]-xyz2[0];        
        dy = xyz1[1]-xyz2[1];        
        dz = xyz1[2]-xyz2[2];        
        length = sqrt(dx*dx+dy*dy+dz*dz);
        npoint = length/light_delta+1.5;
        if(npoint<2)npoint=2;
        for(i=0;i<npoint;i++){
          xyz[0] = ((float)(npoint-1-i)*xyz1[0] + (float)i*xyz2[0])/(float)(npoint-1);
          xyz[1] = ((float)(npoint-1-i)*xyz1[1] + (float)i*xyz2[1])/(float)(npoint-1);
          xyz[2] = ((float)(npoint-1-i)*xyz1[2] + (float)i*xyz2[2])/(float)(npoint-1);
          set_lightfield(smoke3di,xyz,lighti->q/(float)npoint);
        }
        break;
      case 2:      // region
        xyz1 = lighti->xyz1;
        xyz2 = lighti->xyz2;
        dx = abs(xyz1[0]-xyz2[0]);
        dy = abs(xyz1[1]-xyz2[1]);
        dz = abs(xyz1[2]-xyz2[2]);
        nx = dx/light_delta+1.5;
        ny = dy/light_delta+1.5;
        nz = dz/light_delta+1.5;
        dxx = 0.0;
        dyy = 0.0;
        dzz = 0.0;
        if(nx>1){
          dxx = (xyz2[0]-xyz1[0])/(nx-1);
        }
        if(ny>1){
          dyy = (xyz2[1]-xyz1[1])/(ny-1);
        }
        if(nz>1){
          dzz = (xyz2[2]-xyz1[2])/(nz-1);
        }
        for(k=0;k<nz;k++){
          xyz[2] = xyz1[2] + k*dzz;
          for(j=0;j<ny;j++){
            xyz[1] = xyz1[1] + j*dyy;
            for(i=0;i<nx;i++){
              xyz[0] = xyz1[0] + i*dxx;
              set_lightfield(smoke3di,xyz,lighti->q/(float)(nx*ny*nz));
            }
          }
        }
        break;
    }
  }

  // convert hrr field to colors
}

/* ------------------ set_lightfield ------------------------ */

void set_lightfield(smoke3d *smoke3di,float xyz[3], float light_q_source){
  int i,j,k;
  mesh *smoke_mesh;
  float rads[NRAD+1], area[NRAD+1];
  float cos_theta[NTHETA+1], sin_theta[NTHETA+1];
  float cos_psi[NPSI+1], sin_psi[NPSI+1];
  float light_dist;
  float PI, theta, psi;
  int ipsi, irad, itheta;
  int nrad, nradtheta;
  int inode;
  float xbar0, xbar;
  float ybar0, ybar;
  float zbar0, zbar;
  int ibar, jbar, kbar;
  int nx, ny, nxy;

  smoke_mesh = smoke3di->smoke_mesh;
  
  nrad=NRAD;
  nradtheta=NRAD*NTHETA;

  xbar0 = smoke_mesh->xbar0;
  xbar =  smoke_mesh->xbar;
  ybar0 = smoke_mesh->ybar0;
  ybar =  smoke_mesh->ybar;
  zbar0 = smoke_mesh->zbar0;
  zbar =  smoke_mesh->zbar;

  ibar = smoke_mesh->ibar;
  jbar = smoke_mesh->jbar;
  kbar = smoke_mesh->kbar;

  nx = smoke3di->nx;
  ny = smoke3di->ny;
  nxy = nx*ny;

  light_dist=0.0;
  light_dist=getlight_dist(light_dist,xyz,xbar0,ybar0,zbar0);
  light_dist=getlight_dist(light_dist,xyz, xbar,ybar0,zbar0);
  light_dist=getlight_dist(light_dist,xyz,xbar0, ybar,zbar0);
  light_dist=getlight_dist(light_dist,xyz, xbar, ybar,zbar0);
  light_dist=getlight_dist(light_dist,xyz,xbar0,ybar0, zbar);
  light_dist=getlight_dist(light_dist,xyz, xbar,ybar0, zbar);
  light_dist=getlight_dist(light_dist,xyz,xbar0, ybar, zbar);
  light_dist=getlight_dist(light_dist,xyz, xbar, ybar, zbar);

  PI=4.0*atan(1.0);
  for(i=0;i<NRAD;i++){
    float rad;
    rad=(float)(i+1)*light_dist/(float)NRAD;
    rads[i]=rad;
    area[i]=4.0*PI*rad*rad;
  }
  for(i=0;i<NTHETA+1;i++){
    theta = (float)i*2.0*PI/(float)NTHETA;
    cos_theta[i]=cos(theta);
    sin_theta[i]=sin(theta);
  }
  for(i=0;i<NPSI+1;i++){
    psi = (float)i*2.0*PI/(float)NPSI;
    cos_psi[i]=cos(psi);
    sin_psi[i]=sin(psi);
  }

  // set polar hrr field to zero

  for(i=0;i<NRAD*NTHETA*NPSI;i++){
    light_q_polar[i]=0.0;
  }

  // set center of field to hrr/area then
  //   set each successive shell using new_shell_hrrpua = old_shell_hrrpua*(r/(r+dr))^2 (1-alpha)

#define GETPOLARNODE(ipsi,itheta,irad) ((irad)+(itheta)*nrad+(ipsi)*nradtheta)
  inode=0;
  for(ipsi=0;ipsi<NPSI;ipsi++){
    for(itheta=0;itheta<NTHETA;itheta++){
      inode=GETPOLARNODE(ipsi,itheta,0);
      light_q_polar[inode]=light_q_source/area[0];
      for(irad=1;irad<NRAD;irad++){
        float xyz2[3];
        float alpha;

        if(light_q_polar[inode-1]==0.0){
          light_q_polar[inode]=0.0;
        }
        else{
          xyz2[0] = xyz[0] + rads[irad]*cos_theta[itheta]*cos_psi[ipsi];
          xyz2[1] = xyz[1] + rads[irad]*sin_theta[itheta]*cos_psi[ipsi];
          xyz2[2] = xyz[2] + rads[irad]*sin_psi[ipsi];
          alpha = getalpha(smoke_mesh,xyz2);
          light_q_polar[inode]=light_q_polar[inode-1]*area[irad-1]/area[irad]*(1.0-alpha);
        }
        inode++;
      }
    }
  }

  // add polar field to rectangular field
  
  for(k=0;k<kbar;k++){
    float dx, dy, dz;
    float r, theta, psi;
    float xy_length;
    int irad, ipsi, itheta;
    int polarnode, ijknode;

    dz = (zbar0*(float)(kbar-1-k) + (float)k*zbar)/(float)(kbar-1)-xyz[2];
    for(j=0;j<smoke_mesh->jbar;j++){
      dy = (ybar0*(float)(jbar-1-j) + (float)j*ybar)/(float)(jbar-1)-xyz[1];
      for(i=0;i<smoke_mesh->ibar;i++){
        dx = (xbar0*(float)(ibar-1-i) + (float)i*xbar)/(float)(ibar-1)-xyz[0];
        r = sqrt(dx*dx+dy*dy+dz*dz);
        irad=(NRAD-1)*(r/light_dist);
        if(irad<0)irad=0;
        if(irad>NRAD-1)irad=NRAD-1;

        xy_length = sqrt(dx*dx+dy*dy);
        psi = atan3(dz,xy_length)+PI/2.0;
        ipsi = (NPSI-1)*psi/PI;
        if(ipsi<0)ipsi=0;
        if(ipsi>NPSI-1)ipsi=NPSI-1;

        theta = 0.0;
        if(dx!=0.0||dy!=0.0){
          theta = atan2(dy,dx);
        }
        theta+=PI;
        itheta = theta*(NTHETA-1)*theta/(2.0*PI);
        if(itheta<0)itheta=0;
        if(itheta>NTHETA-1)itheta=NTHETA-1;

        polarnode = GETPOLARNODE(ipsi,itheta,irad);
        ijknode = IJKNODE(i,j,k);
        smoke3di->light_q_rect[ijknode]+=light_q_polar[polarnode];
      }
    }
  }
}

/* ------------------ getlight_dist ------------------------ */

float getlight_dist(float light_dist,float *xyz, float x, float y, float z){
  float dx, dy, dz;
  float dist;

  dx = x-xyz[0];
  dy = y-xyz[1];
  dz = z-xyz[2];
  dist = sqrt(dx*dx+dy*dy+dz*dz);
  if(light_dist>dist)dist=light_dist;
  return dist;
}

/* ------------------ getalpha ------------------------ */

float getalpha(mesh *smoke_mesh, float *xyz2){
  int i, j, k;
  int nx, ny, nxy;
  int ialpha;

  if(xyz2[0]<smoke_mesh->xbar0||xyz2[0]>smoke_mesh->xbar)return 1.0;
  if(xyz2[1]<smoke_mesh->ybar0||xyz2[1]>smoke_mesh->ybar)return 1.0;
  if(xyz2[2]<smoke_mesh->zbar0||xyz2[2]>smoke_mesh->zbar)return 1.0;

  nx = smoke_mesh->ibar;
  ny = smoke_mesh->jbar;
  nxy = nx*ny;

  i = (xyz2[0]-smoke_mesh->xbar0)/(smoke_mesh->xbar-smoke_mesh->xbar0)*smoke_mesh->ibar;
  if(i<0)i=0;
  if(i>smoke_mesh->ibar)i=smoke_mesh->ibar;

  j = (xyz2[1]-smoke_mesh->ybar0)/(smoke_mesh->ybar-smoke_mesh->ybar0)*smoke_mesh->jbar;
  if(j<0)j=0;
  if(j>smoke_mesh->jbar)j=smoke_mesh->jbar;

  k = (xyz2[2]-smoke_mesh->zbar0)/(smoke_mesh->zbar-smoke_mesh->zbar0)*smoke_mesh->kbar;
  if(k<0)k=0;
  if(k>smoke_mesh->kbar)k=smoke_mesh->kbar;

  ialpha = full_alphabuffer[IJKNODE(i,j,k)];
  return (float)ialpha/255.0;
}
#endif


