// $Date: 2011-03-30 12:04:57 -0400 (Wed, 30 Mar 2011) $ 
// $Revision: 8021 $
// $Author: gforney $

#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "interpdata.h"

// svn revision character string
char interpdata_revision[]="$Revision: 8021 $";

/* ----------------------- create_lightmap ----------------------------- */

#define IJK(i,j,k) (i) + ny*(j) + nxy*(k)

void create_lightmap(lightdata *lightinfo){
  int i, j, k, nx, ny, nz, nxy;
  float *flightmap;
  unsigned char *lightmap, *opacity;

  nx = lightinfo->ibar+1;
  ny = lightinfo->jbar+1;
  nz = lightinfo->kbar+1;
  nxy = nx*ny;
  lightmap = lightinfo->lightmap;
  flightmap = lightinfo->flightmap;
  opacity = lightinfo->opacity;

  i=0;
  for(k=0;k<nz;k++){
  for(j=0;j<ny;j++){
    flightmap[IJK(i,j,k)]=1.0;
  }
  }

  for(k=0;k<nz;k++){
  for(j=0;j<ny;j++){
  for(i=1;i<nx;i++){
    int ijk,im1jk;

    ijk=IJK(i,j,k);
    im1jk=ijk-1;
    flightmap[ijk]=(float)(flightmap[im1jk]*(float)(opacity[im1jk])/255.0);
  }
  }
  }

  for(i=0;i<nx*ny*nz;i++){
    lightmap[i]=(unsigned char)(flightmap[i]*255.0);
  }

}

/* ----------------------- interp3d ----------------------------- */

unsigned char interp3d(lightdata *lightinfo, float xyz[3]){
  int i, j, k;
  int ijk;

  i = (xyz[0]-lightinfo->xbar0)/lightinfo->dx;
  j = (xyz[1]-lightinfo->ybar0)/lightinfo->dy;
  k = (xyz[2]-lightinfo->zbar0)/lightinfo->dz;

  ijk = i + j*lightinfo->ibar + k*lightinfo->ibar*lightinfo->jbar;

  return lightinfo->opacity[ijk];
}
