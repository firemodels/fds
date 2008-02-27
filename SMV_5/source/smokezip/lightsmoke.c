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

float get_photon_step(smoke3d *smoke3di, float xyzpos[3], float xyzdir[3]);
int in_mesh(mesh *smoke_mesh,float xyzpos[3]);
float getalpha(mesh *smoke_mesh, float *xyz2);
int get_interval(float val, float *array, int n);
float interp_char(float xyz[3], mesh *smoke_mesh, unsigned char *full_alphabuffer);

#define IJKNODE(i,j,k) ((i)+(j)*nx+(k)*nxy)

/* ------------------ get_random_light ------------------------ */

lightdata *get_random_light(void){
  lightdata *light;
  int ilight;
  float val;

  if(nlightinfo==1)return lightinfo;
  val = rand_1d(0.0,1.0);
  ilight = get_interval(val, light_cdf, nlightinfo+1);
  light = lightinfo + ilight;
  return light;
}

/* ------------------ update_lightfield ------------------------ */

void update_lightfield(smoke3d *smoke3di, unsigned char *lightingbuffer){
 // accumulate hrr for each light

    lightdata *lighti;
    float xyzpos[3], xyzdir[3];
    float factor, *xyz1, *xyz2;
    int i;
    int *photon_bin;
    float photon_step;
    mesh *smoke_mesh;
    int nx, ny, nz, nxy;
    int binmax;

#define NPHOTONS 100000

    // pick a random light weighted by HRR  
    //  (ie a high wattage light will be picked proportionately more often than
    //          a low wattage one)

    smoke_mesh=smoke3di->smoke_mesh;
    nx = smoke_mesh->ibar;
    ny = smoke_mesh->jbar;
    nz = smoke_mesh->kbar;
    nxy = nx*ny;

// zero out bin used to collect photons

    photon_bin = smoke3di->smoke_mesh->photon_bin;
    for(i=0;i<nx*ny*nz;i++){
      photon_bin[i]=0;
    }
    CheckMemory;

    for(i=0;i<NPHOTONS;i++){
      lighti = get_random_light();

    //  get a random position within the light and a random direction

      switch (lighti->type){

        case 0:      // point
          xyzpos[0]=lighti->xyz1[0];
          xyzpos[1]=lighti->xyz1[1];
          xyzpos[2]=lighti->xyz1[2];
          rand_dir(xyzdir);
          break;

        case 1:      // line
          factor = rand_1d(0.0,1.0);
          xyz1 = lighti->xyz1;
          xyz2 = lighti->xyz2;
          xyzpos[0]= xyz1[0]*(1-factor) + xyz2[0]*factor;
          xyzpos[1]= xyz1[1]*(1-factor) + xyz2[1]*factor;
          xyzpos[2]= xyz1[2]*(1-factor) + xyz2[2]*factor;
          rand_dir(xyzdir);
          break;

        case 2:      // region
          xyz1 = lighti->xyz1;
          xyz2 = lighti->xyz2;
          factor = rand_1d(0.0,1.0);
          xyzpos[0]= xyz1[0]*(1-factor) + xyz2[0]*factor;
          factor = rand_1d(0.0,1.0);
          xyzpos[1]= xyz1[1]*(1-factor) + xyz2[1]*factor;
          factor = rand_1d(0.0,1.0);
          xyzpos[2]= xyz1[2]*(1-factor) + xyz2[2]*factor;
          rand_absdir(xyzdir,lighti->dir);
          break;
      }
      for(;;){
        int i1, j1, k1, ijk;

        photon_step=get_photon_step(smoke3di, xyzpos, xyzdir);
        xyzpos[0]+= photon_step*xyzdir[0];
        xyzpos[1]+= photon_step*xyzdir[1];
        xyzpos[2]+= photon_step*xyzdir[2];
        if(in_mesh(smoke_mesh,xyzpos)==0)break;

        i1 = get_interval(xyzpos[0],smoke_mesh->xplt,smoke_mesh->ibar+1);
        j1 = get_interval(xyzpos[1],smoke_mesh->yplt,smoke_mesh->jbar+1);
        k1 = get_interval(xyzpos[2],smoke_mesh->zplt,smoke_mesh->kbar+1);

        ijk = IJKNODE(i1,j1,k1);
        photon_bin[ijk]++;                 // record location of photon

        if(rand_1d(0.0,1.0)>albedo)break;   // if this if is true then the photon is absorbed
        rand_dir(xyzdir);
      }
      CheckMemory;
    }
    binmax = 0;
    for(i=0;i<nx*ny*nz;i++){
      if(photon_bin[i]>binmax)binmax=photon_bin[i];
    }
    CheckMemory;
    for(i=0;i<nx*ny*nz;i++){
      lightingbuffer[i]=(unsigned char)((float)photon_bin[i]/(float)binmax*254.0);
    }
    CheckMemory;
}

/* ------------------ in_mesh ------------------------ */

int in_mesh(mesh *smoke_mesh,float xyzpos[3]){
  if(xyzpos[0]<smoke_mesh->xbar0||xyzpos[0]>smoke_mesh->xbar)return 0;
  if(xyzpos[1]<smoke_mesh->ybar0||xyzpos[1]>smoke_mesh->ybar)return 0;
  if(xyzpos[2]<smoke_mesh->zbar0||xyzpos[2]>smoke_mesh->zbar)return 0;
  return 1;

}

/* ------------------ get_photon_step ------------------------ */

float get_photon_step(smoke3d *smoke3di, float xyzpos[3], float xyzdir[3]){
#define NPHOTON_PDF 50
  float photon_pdf[NPHOTON_PDF];
  int i;
  float xyz[3];
  float alpha;
  mesh *smoke_mesh;
  float var, dlength;
  int ivar, ibreak;
  float step;

  smoke_mesh = smoke3di->smoke_mesh;


  dlength = smoke_mesh->dxyzmax/NPHOTON_PDF;
  ibreak = 0;

  photon_pdf[0]=0.0;
  for(i=1;i<NPHOTON_PDF;i++){
    if(ibreak==0){
      xyz[0] = xyzpos[0] + i*dlength*xyzdir[0];
      xyz[1] = xyzpos[1] + i*dlength*xyzdir[1];
      xyz[2] = xyzpos[2] + i*dlength*xyzdir[2];
      if(in_mesh(smoke_mesh,xyz)==0){
        ibreak=1;
        break;
      }
      alpha = getalpha(smoke_mesh,xyz);
      if(alpha>0.999)alpha=0.999;
      photon_pdf[i]=photon_pdf[i-1]+log(1.0-alpha)/smoke_mesh->dx;
    }
    else{
      photon_pdf[i]=photon_pdf[i-1]-10.0;
    }
  }
  for(i=0;i<NPHOTON_PDF;i++){
    if(photon_pdf[i]<-50.0){
      photon_pdf[i]=1.0;
    }
    else{
      photon_pdf[i]=1.0-exp(photon_pdf[i]);
    }
  }
  var = rand_1d(0.0,1.0);
  ivar = get_interval(var,photon_pdf,NPHOTON_PDF);
  step = ivar*dlength;
  return step;
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
  float f1=0.0, f2=1.0, g1=0.0, g2=1.0, h1=0.0, h2=1.0;
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
  dx = x2 - x1;
  if(dx!=0.0)f1/=dx;
  f2 = 1.0 - f1;

  g1 = xyz[1]-y1;
  dy = y2 - y1;
  if(dy!=0.0)g1/=dy;
  g2 = 1.0 - g1;

  h1 = xyz[2]-z1;
  dz = z2 - z1;
  if(dz!=0.0)h1/=dz;
  h2 = 1.0 - h1;

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

  if(val<array[0])return 0;
  if(val>array[n-1])return n-1;

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


