#include "options.h"
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "svdiff.h"
#include "MALLOC.h"
#include "datadefs.h"

/* ------------------ mesh_match ------------------------ */

int mesh_match(mesh *mesh1, mesh *mesh2){
  int ibar, jbar, kbar;

  if(mesh1->ibar!=mesh2->ibar)return 0;
  if(mesh1->jbar!=mesh2->jbar)return 0;
  if(mesh1->kbar!=mesh2->kbar)return 0;
  ibar=mesh1->ibar;
  jbar=mesh1->jbar;
  kbar=mesh1->kbar;
  if(ABS(mesh1->xplt[0]-mesh2->xplt[0])>mesh1->dx)return 0;
  if(ABS(mesh1->yplt[0]-mesh2->yplt[0])>mesh1->dy)return 0;
  if(ABS(mesh1->zplt[0]-mesh2->zplt[0])>mesh1->dz)return 0;
  if(ABS(mesh1->xplt[ibar]-mesh2->xplt[ibar])>mesh1->dx)return 0;
  if(ABS(mesh1->yplt[jbar]-mesh2->yplt[jbar])>mesh1->dy)return 0;
  if(ABS(mesh1->zplt[kbar]-mesh2->zplt[kbar])>mesh1->dz)return 0;
  return 1;
}

/* ------------------ version ------------------------ */

void version(void){
    char smv_version[100];
    char revision[100];

    getPROGversion(smv_version);  // get Smokeview version (ie 5.x.z)
    getRevision(revision);    // get svn revision number
    PRINTF("\n");
    PRINTF("Smokediff\n\n");
    PRINTF("Version: %s\n",smv_version);
    PRINTF("Revision: %s\n",revision);
    PRINTF("Compile Date: %s\n",__DATE__);
#ifdef X64
    PRINTF("Platform: WIN64\n");
#endif
#ifdef WIN32
#ifndef X64
    PRINTF("Platform: WIN32\n");
#endif
#endif
#ifndef pp_OSX64
#ifdef pp_OSX
    PRINTF("Platform: OSX\n");
#endif
#endif
#ifdef pp_OSX64
    PRINTF("Platform: OSX64\n");
#endif
#ifndef pp_LINUX64
#ifdef pp_LINUX
    PRINTF("Platform: LINUX\n");
#endif
#endif
#ifdef pp_LINUX64
    PRINTF("Platform: LINUX64\n");
#endif


}

/* ------------------ similar_grid ------------------------ */

int similar_grid(mesh *mesh1, mesh *mesh2, int *factor){

  factor[0]=1;
  factor[1]=1;
  factor[2]=1;
  
  if(ABS( mesh1->xbar0-mesh2->xbar0)>mesh1->dx/2.0)return 0;
  if(ABS( mesh1->xbar- mesh2->xbar )>mesh1->dx/2.0)return 0;
  if(ABS( mesh1->ybar0-mesh2->ybar0)>mesh1->dy/2.0)return 0;
  if(ABS( mesh1->ybar- mesh2->ybar )>mesh1->dy/2.0)return 0;
  if(ABS( mesh1->zbar0-mesh2->zbar0)>mesh1->dz/2.0)return 0;
  if(ABS( mesh1->zbar- mesh2->zbar )>mesh1->dz/2.0)return 0;

  factor[0] = mesh2->ibar/mesh1->ibar;
  if(mesh1->ibar*factor[0]!=mesh2->ibar)return 0;

  factor[1] = mesh2->jbar/mesh1->jbar;
  if(mesh1->jbar*factor[1]!=mesh2->jbar)return 0;

  factor[2] = mesh2->kbar/mesh1->kbar;
  if(mesh1->kbar*factor[2]!=mesh2->kbar)return 0;

  return 1;
}

/* ------------------ exact_grid ------------------------ */

int exact_grid(mesh *mesh1, mesh *mesh2, int *factor){
  float eps;

  factor[0]=1;
  factor[1]=1;
  factor[2]=1;
  eps = mesh1->dx/1000.0;
  if(ABS(mesh1->dx-mesh2->dx)>eps)return 0;
  eps = mesh1->dy/1000.0;
  if(ABS(mesh1->dy-mesh2->dy)>eps)return 0;
  eps = mesh1->dz/1000.0;
  if(ABS(mesh1->dz-mesh2->dz)>eps)return 0;
  return 1;
}
