// $Date: 2011-03-30 12:04:57 -0400 (Wed, 30 Mar 2011) $ 
// $Revision: 8021 $
// $Author: gforney $

#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "interpdata.h"
#include "MALLOC.h"

// svn revision character string
char interpdata_revision[]="$Revision: 8021 $";

/* ----------------------- setup_radiancemap ----------------------------- */

void setup_radiancemap(radiancedata *radianceinfo, int ijkbar[3], float xyzbar0[3], float xyzbar[3], float dxyz[3], unsigned char *radiance, unsigned char *opacity){
  radianceinfo->ijkbar=ijkbar;
  radianceinfo->xyzbar0=xyzbar0;
  radianceinfo->xyzbar =xyzbar;
  radianceinfo->radiance=radiance;
  radianceinfo->opacity=opacity;
  radianceinfo->dxyz=dxyz;
}

/* ----------------------- build_radiancemap ----------------------------- */

#define IJKRAD(i,j,k) (i) + ny*(j) + nxy*(k)

void build_radiancemap(radiancedata *radianceinfo){
  int i, j, k, nx, ny, nz, nxy;
  float *fradiance;
  unsigned char *radiance, *opacity;

  nx = radianceinfo->ijkbar[0]+1;
  ny = radianceinfo->ijkbar[1]+1;
  nz = radianceinfo->ijkbar[2]+1;
  nxy = nx*ny;
  NewMemory((void **)&fradiance,nx*ny*nz*sizeof(float));

  radiance = radianceinfo->radiance;
  opacity = radianceinfo->opacity;

  i=0;
  for(k=0;k<nz;k++){
  for(j=0;j<ny;j++){
    fradiance[IJKRAD(i,j,k)]=1.0;
  }
  }

  for(k=0;k<nz;k++){
  for(j=0;j<ny;j++){
  for(i=1;i<nx;i++){
    int ijk,im1jk;

    ijk=IJKRAD(i,j,k);
    im1jk=ijk-1;
    fradiance[ijk]=(float)(fradiance[im1jk]*(float)(opacity[im1jk])/255.0);
  }
  }
  }

  for(i=0;i<nx*ny*nz;i++){
    radiance[i]=(unsigned char)(fradiance[i]*255.0);
  }
  FREEMEMORY(fradiance);
}

/* ----------------------- interp3d ----------------------------- */

unsigned char get_opacity(radiancedata *radianceinfo, float xyz[3]){
  int i, j, k;
  int ijk;

  i = (xyz[0]-radianceinfo->xyzbar0[0])/radianceinfo->dxyz[0];
  j = (xyz[1]-radianceinfo->xyzbar0[1])/radianceinfo->dxyz[1];
  k = (xyz[2]-radianceinfo->xyzbar0[2])/radianceinfo->dxyz[2];

  ijk = i + j*radianceinfo->ijkbar[0] + k*radianceinfo->ijkbar[0]*radianceinfo->ijkbar[1];

  return radianceinfo->opacity[ijk];
}
