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

float getalpha(mesh *smoke_mesh, float *xyz2);
int get_interval(float val, float *array, int n);
float interp_char(float xyz[3], mesh *smoke_mesh, unsigned char *full_alphabuffer);

#define IJKNODE(i,j,k) ((i)+(j)*nx+(k)*nxy)

/* ------------------ update_lightfield ------------------------ */

void update_lightfield(smoke3d *smoke3di, unsigned char *lightingbuffer){
}

/* ------------------ getalpha ------------------------ */

float getalpha(mesh *smoke_mesh, float *xyz){
  float val;

  if(xyz[0]<smoke_mesh->xbar0||xyz[0]>smoke_mesh->xbar)return 1.0;
  if(xyz[1]<smoke_mesh->ybar0||xyz[1]>smoke_mesh->ybar)return 1.0;
  if(xyz[2]<smoke_mesh->zbar0||xyz[2]>smoke_mesh->zbar)return 1.0;

  val = interp_char(xyz, smoke_mesh, full_alphabuffer);
  val /= 255.0;
  return val;
}

/* ------------------ interp_char ------------------------ */

float interp_char(float xyz[3], mesh *smoke_mesh, unsigned char *full_alphabuffer){
  int ialpha111, ialpha112, ialpha121, ialpha122;
  int ialpha211, ialpha212, ialpha221, ialpha222;
  float x1, x2, y1, y2, z1, z2;
  float f1=0.0, f2=0.0, g1=0.0, g2=0.0, h1=0.0, h2=0.0;
  float val;
  int i1, j1, k1;
  int i2, j2, k2;
  int nx, ny, nxy;
  float dx, dy, dz;

  nx = smoke_mesh->ibar;
  ny = smoke_mesh->jbar;
  nxy = nx*ny;

  i1 = get_interval(xyz[0],smoke_mesh->xplt,smoke_mesh->ibar+1);
  j1 = get_interval(xyz[1],smoke_mesh->yplt,smoke_mesh->jbar+1);
  k1 = get_interval(xyz[2],smoke_mesh->zplt,smoke_mesh->kbar+1);

  i2=i1+1;
  if(i2>smoke_mesh->ibar)i2=smoke_mesh->ibar;
  j2=j1+1;
  if(j2>smoke_mesh->jbar)j2=smoke_mesh->jbar;
  k2=k1+1;
  if(k2>smoke_mesh->kbar)k2=smoke_mesh->kbar;
  

  ialpha111 = (int)full_alphabuffer[IJKNODE(i1,j1,k1)];
  ialpha112 = (int)full_alphabuffer[IJKNODE(i1,j1,k2)];
  ialpha121 = (int)full_alphabuffer[IJKNODE(i1,j2,k1)];
  ialpha122 = (int)full_alphabuffer[IJKNODE(i1,j2,k2)];
  ialpha211 = (int)full_alphabuffer[IJKNODE(i2,j1,k1)];
  ialpha212 = (int)full_alphabuffer[IJKNODE(i2,j1,k2)];
  ialpha221 = (int)full_alphabuffer[IJKNODE(i2,j2,k1)];
  ialpha222 = (int)full_alphabuffer[IJKNODE(i2,j2,k2)];

  x1 = smoke_mesh->xplt[i1];
  y1 = smoke_mesh->yplt[j1];
  z1 = smoke_mesh->zplt[k1];
  x2 = smoke_mesh->xplt[i2];
  y2 = smoke_mesh->yplt[j2];
  z2 = smoke_mesh->zplt[k2];

  f1 = xyz[0]-x1;
  f2 = x2 - xyz[0];
  dx = x2 - x1;
  if(dx!=0.0){
    f1/=dx;
    f2/=dx;
  }
  g1 = xyz[1]-y1;
  g2 = y2 - xyz[1];
  dy = y2 - y1;
  if(dy!=0.0){
    g1/=dy;
    g2/=dy;
  }
  h1 = xyz[2]-z1;
  h2 = z2 - xyz[2];
  dz = z2 - z1;
  if(dz!=0.0){
    h1/=dz;
    h2/=dz;
  }

  val  = f2*g2*h2*ialpha111;
  val += f2*g2*h1*ialpha112;
  val += f2*g1*h2*ialpha121;
  val += f2*g1*h1*ialpha122;
  val += f1*g2*h2*ialpha211;
  val += f1*g2*h1*ialpha212;
  val += f1*g1*h2*ialpha221;
  val += f1*g1*h1*ialpha222;
  return val;
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


#endif


