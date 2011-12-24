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
#include "datadefs.h"

char blockage_test_revision[]="$Revision$";

#define HIT 1
#define MISS 0

float llasttime=0.0;

void adjust_new_position(float oldpos[3], float newpos[3]);

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

/* ------------------ drawMovedir ------------------------ */

void drawMovedir(void){
  float dir_color[]={1.0,0.0,0.0,1.0};
  float dx, dy, dz;

#define DZ 0.005
#define DXZ (DZ/16.0)

  float xyz00[3], xyz01[3], xyz10[3], xyz11[3], *xyz0, xyz1[3];

  dx=-DXZ*movedir[1];
  dy= DXZ*movedir[0];
  dz= sqrt(dx*dx+dy*dy);
#ifdef _DEBUG
//  printf("dx=%f dy=%f dz=%f\n",dx,dy,dz);
#endif

  xyz0=eye_xyz0;
  xyz1[0]=xyz0[0]+DZ*movedir[0];
  xyz1[1]=xyz0[1]+DZ*movedir[1];
  xyz1[2]=xyz0[2]+DZ*movedir[2];

  xyz00[0]=xyz1[0]-dx;
  xyz00[1]=xyz1[1]-dy;
  xyz00[2]=xyz1[2]-dz;

  xyz01[0]=xyz1[0]-dx;
  xyz01[1]=xyz1[1]-dy;
  xyz01[2]=xyz1[2]+dz;

  xyz11[0]=xyz1[0]+dx;
  xyz11[1]=xyz1[1]+dy;
  xyz11[2]=xyz1[2]+dz;

  xyz10[0]=xyz1[0]+dx;
  xyz10[1]=xyz1[1]+dy;
  xyz10[2]=xyz1[2]-dz;

  glDepthMask(GL_FALSE);
  glPointSize(10.0);
  glBegin(GL_POINTS);
  glColor3fv(dir_color);
  glVertex3fv(xyz00);
  glVertex3fv(xyz01);
  glVertex3fv(xyz10);
  glVertex3fv(xyz11);
  glEnd();

  glBegin(GL_LINES);
  glColor3fv(dir_color);
  glVertex3fv(xyz00);
  glVertex3fv(xyz11);

  glVertex3fv(xyz01);
  glVertex3fv(xyz10);

  glVertex3fv(xyz00);
  glVertex3fv(xyz01);

  glVertex3fv(xyz01);
  glVertex3fv(xyz11);
  
  glVertex3fv(xyz11);
  glVertex3fv(xyz10);

  glVertex3fv(xyz10);
  glVertex3fv(xyz00);
  glEnd();
  glDepthMask(GL_TRUE);
}

/* ------------------ setspeed ------------------------ */

  void setspeed(float tospeed){
    speed_I=0;
    speed_desired=tospeed;
  }

/* ------------------ getnewpos ------------------------ */

void getnewpos(float *oldpos, float dx, float dy, float dz,float local_speed_factor){
  oldpos[0] += dx;
  oldpos[1] += dy;
  oldpos[2] += dz;
  from_glui_trainer=0;
}

/* ------------------ getmesh ------------------------ */

mesh *getmesh(float *xyz){
  mesh *meshi;
  int i;
  int ibar, jbar, kbar;
  float xmin, xmax;
  float ymin, ymax;
  float zmin, zmax;
  float *xplt, *yplt, *zplt;

  for(i=0;i<nmeshes;i++){
    meshi = meshinfo+i;

    ibar = meshi->ibar;
    jbar = meshi->jbar;
    kbar = meshi->kbar;

    xplt = meshi->xplt_orig;
    yplt = meshi->yplt_orig;
    zplt = meshi->zplt_orig;

    xmin = xplt[0];
    xmax = xplt[ibar];
    if(xyz[0]<xmin||xyz[0]>xmax)continue;

    ymin = yplt[0];
    ymax = yplt[jbar];
    if(xyz[1]<ymin||xyz[1]>ymax)continue;

    zmin = zplt[0];
    zmax = zplt[kbar];
    if(xyz[2]<zmin||xyz[2]>zmax)continue;

    return meshi;
  }
  return NULL;
}

/* ------------------ getblockdist ------------------------ */

float getblockdist(float x, float y, float z){
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
      if(iblank_cell[ijkcell]==0)return 0.0;
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

void init_blockdist(void){
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
          zdist=getblockdist(xx,yy,zz-dz/2.0);
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
            if(iblank_cell[ijkm1cell]==0){
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
  int n_embed;
  char *ib_embed;
  char *iblank_embed;

#define MESHIJ(i,j) (i)*nmeshes + (j)

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
      ib_embed[j]=1;
    }
    for(j=0;j<nmeshes;j++){
      mesh *meshj;
      int i1, i2, j1, j2, k1, k2;
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
          j1=jj;
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
        for(jj=j1;jj<=j2;jj++){
          for(ii=i1;ii<=i2;ii++){
            ib_embed[IJKNODE(ii,jj,kk)]=0;
          }
        }
      }
    }
  }
  return 0;
}

/* ------------------ makeiblank ------------------------ */

int makeiblank(void){
  blockagedata *bc;
  int ijksize,i,j,k;
  int ii,ig;
  int test;
  mesh *meshi;
  int ibar,jbar,kbar;
  int nx, ny, nxy;
  float *fblank_cell=NULL;
  char *iblank=NULL,*iblank_cell=NULL,*iblank_x=NULL,*iblank_y=NULL,*iblank_z=NULL;

  printf("initializing blanking arrays\n");
  if(use_iblank==0)return 0;
  for(ig=0;ig<nmeshes;ig++){
    meshi = meshinfo+ig;
    ibar = meshi->ibar;
    jbar = meshi->jbar;
    kbar = meshi->kbar;
    ijksize=(ibar+1)*(jbar+1)*(kbar+1);

    if(NewMemory((void **)&iblank,ijksize*sizeof(char))==0)return 1;
    if(NewMemory((void **)&iblank_cell,ibar*jbar*kbar*sizeof(char))==0)return 1;
    if(NewMemory((void **)&fblank_cell,ibar*jbar*kbar*sizeof(float))==0)return 1;
    if(NewMemory((void **)&iblank_x,ijksize*sizeof(char))==0)return 1;
    if(NewMemory((void **)&iblank_y,ijksize*sizeof(char))==0)return 1;
    if(NewMemory((void **)&iblank_z,ijksize*sizeof(char))==0)return 1;

    meshi->c_iblank=iblank;
    meshi->c_iblank_cell=iblank_cell;
    meshi->f_iblank_cell=fblank_cell;
    meshi->c_iblank_x=iblank_x;
    meshi->c_iblank_y=iblank_y;
    meshi->c_iblank_z=iblank_z;

    for(i=0;i<ibar*jbar*kbar;i++){
      iblank_cell[i]=1;
    }
    for(i=0;i<ijksize;i++){
      iblank[i]=1;
      iblank_x[i]=1;
      iblank_y[i]=1;
      iblank_z[i]=1;
    }

    nx = ibar+1;
    ny = jbar+1;
    nxy = nx*ny;

    for(ii=0;ii<meshi->nbptrs;ii++){
      bc=meshi->blockageinfoptrs[ii];
      for(i=bc->ijk[IMIN];i<bc->ijk[IMAX];i++){
      for(j=bc->ijk[JMIN];j<bc->ijk[JMAX];j++){
      for(k=bc->ijk[KMIN];k<bc->ijk[KMAX];k++){
        iblank_cell[IJKCELL(i,j,k)]=0;
      }
      }
      }
    }
    if(fblank_cell!=NULL){
      for(ii=0;ii<ibar*jbar*kbar;ii++){
        fblank_cell[ii]=iblank_cell[ii];
      }
    }
    i=0;
    for(j=1;j<jbar;j++){
    for(k=1;k<kbar;k++){
      test=0;
      test+=iblank_cell[IJKCELL(  i,j-1,k-1)];
      test+=iblank_cell[IJKCELL(  i,  j,k-1)];
      test+=iblank_cell[IJKCELL(  i,j-1,  k)];
      test+=iblank_cell[IJKCELL(  i,  j,  k)];
      if(test==0)iblank[IJKNODE(i,j,k)]=0;
    }
    }
    j=0;
    for(i=1;i<ibar;i++){
    for(k=1;k<kbar;k++){
      test=0;
      test+=iblank_cell[IJKCELL(i-1,  j,k-1)];
      test+=iblank_cell[IJKCELL(  i,  j,k-1)];
      test+=iblank_cell[IJKCELL(i-1,  j,  k)];
      test+=iblank_cell[IJKCELL(  i,  j,  k)];
      if(test==0)iblank[IJKNODE(i,j,k)]=0;
    }
    }
    k=0;
    for(i=1;i<ibar;i++){
    for(j=1;j<jbar;j++){
      test=0;
      test+=iblank_cell[IJKCELL(i-1,j-1,  k)];
      test+=iblank_cell[IJKCELL(  i,j-1,  k)];
      test+=iblank_cell[IJKCELL(i-1,  j,  k)];
      test+=iblank_cell[IJKCELL(  i,  j,  k)];
      if(test==0)iblank[IJKNODE(i,j,k)]=0;
    }
    }
    i=ibar;
    for(j=1;j<jbar;j++){
    for(k=1;k<kbar;k++){
      test=0;
      test+=iblank_cell[IJKCELL(i-1,j-1,k-1)];
      test+=iblank_cell[IJKCELL(i-1,  j,k-1)];
      test+=iblank_cell[IJKCELL(i-1,j-1,  k)];
      test+=iblank_cell[IJKCELL(i-1,  j,  k)];
      if(test==0)iblank[IJKNODE(i,j,k)]=0;
    }
    }
    j=jbar;
    for(i=1;i<ibar;i++){
    for(k=1;k<kbar;k++){
      test=0;
      test+=iblank_cell[IJKCELL(i-1,j-1,k-1)];
      test+=iblank_cell[IJKCELL(  i,j-1,k-1)];
      test+=iblank_cell[IJKCELL(i-1,j-1,  k)];
      test+=iblank_cell[IJKCELL(  i,j-1,  k)];
      if(test==0)iblank[IJKNODE(i,j,k)]=0;
    }
    }
    k=kbar;
    for(i=1;i<ibar;i++){
    for(j=1;j<jbar;j++){
      test=0;
      test+=iblank_cell[IJKCELL(i-1,j-1,k-1)];
      test+=iblank_cell[IJKCELL(  i,j-1,k-1)];
      test+=iblank_cell[IJKCELL(i-1,  j,k-1)];
      test+=iblank_cell[IJKCELL(  i,  j,k-1)];
      if(test==0)iblank[IJKNODE(i,j,k)]=0;
    }
    }

    for(i=1;i<ibar;i++){
    for(j=1;j<jbar;j++){
    for(k=1;k<kbar;k++){
      test=0;
      test+=iblank_cell[IJKCELL(i-1,j-1,k-1)];
      test+=iblank_cell[IJKCELL(  i,j-1,k-1)];
      test+=iblank_cell[IJKCELL(i-1,  j,k-1)];
      test+=iblank_cell[IJKCELL(  i,  j,k-1)];
      test+=iblank_cell[IJKCELL(i-1,j-1,  k)];
      test+=iblank_cell[IJKCELL(  i,j-1,  k)];
      test+=iblank_cell[IJKCELL(i-1,  j,  k)];
      test+=iblank_cell[IJKCELL(  i,  j,  k)];
      if(test==0)iblank[IJKNODE(i,j,k)]=0;
    }
    }
    }

    for(j=0;j<jbar;j++){
    for(k=0;k<kbar;k++){
      iblank_x[IJKNODE(0,j,k)]   =2*iblank_cell[IJKCELL(0,j,k)];
      for(i=1;i<ibar;i++){
        iblank_x[IJKNODE(i,j,k)]=iblank_cell[IJKCELL(i-1,j,k)]+iblank_cell[IJKCELL(i,j,k)];
      }
      iblank_x[IJKNODE(ibar,j,k)]=2*iblank_cell[IJKCELL(ibar-1,j,k)];
    }
    }
    for(i=0;i<ibar;i++){
    for(k=0;k<kbar;k++){
      iblank_y[IJKNODE(i,0,k)]=2*iblank_cell[IJKCELL(i,0,k)];
      for(j=1;j<jbar;j++){
        iblank_y[IJKNODE(i,j,k)]=iblank_cell[IJKCELL(i,j-1,k)]+iblank_cell[IJKCELL(i,j,k)];
      }
      iblank_y[IJKNODE(i,jbar,k)]=2*iblank_cell[IJKCELL(i,jbar-1,k)];
    }
    }

    for(i=0;i<ibar;i++){
    for(j=0;j<jbar;j++){
      iblank_z[IJKNODE(i,j,0)]=2*iblank_cell[IJKCELL(i,j,0)];
      for(k=1;k<kbar;k++){
        iblank_z[IJKNODE(i,j,k)]=iblank_cell[IJKCELL(i,j,k-1)]+iblank_cell[IJKCELL(i,j,k)];
      }
      iblank_z[IJKNODE(i,j,kbar)]=2*iblank_cell[IJKCELL(i,j,kbar-1)];
    }
    }
  }
   init_blockdist();
  return 0;
}
