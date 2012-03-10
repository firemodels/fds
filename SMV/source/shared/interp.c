// $Date$ 
// $Revision$
// $Author$

// svn revision character string
char interp_revision[]="$Revision$";

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MALLOC.h"
#include "interp.h"

/* ------------------ get_z_interp_factors ------------------------ */

void get_z_interp_factors(float *zplt, int nz, float z, int *k1, int *k2, float *f1, float *f2){
  float dz;
  int ileft, iright;

  dz = zplt[1] - zplt[0];

  ileft = (z-zplt[0])/dz;
  if(ileft<0)ileft=0;
  if(ileft>nz-1)ileft=nz-1;
  iright = ileft + 1;

  *k1 = ileft;
  *k2 = iright;
  *f1 = (z-zplt[ileft])/dz;
  *f2 = (zplt[iright]-z)/dz;
  return;
}

/* ------------------ interp3dsliceindex ------------------------ */

int interp3dsliceindex(unsigned char *data, float *zplt, int nz, int n0, float z){
  int k1, k2;
  float dz;
  float val1, val2;
  float z1, z2;
  int ival;

  dz = zplt[1] - zplt[0];

  k1 = (z-zplt[0])/dz;
  if(k1<0)k1=0;
  if(k1>nz-1)k1=nz-1;
  k2 = k1 + 1;

  val1 = data[n0+k1];
  val2 = data[n0+k2];
  z1 = zplt[k1];
  z2 = zplt[k2];
  ival = ((z-z1)*val2 + (z2-z)*val1)/dz;
  if(ival<0)ival=0;
  if(ival>255)ival=255;
  return ival;
}


/* --------------------------  vertdata ------------------------------------ */

typedef struct _vertdata {
  float xyz[3];
  int nvertring,ntris;
  struct _vertdata **vertring;
  struct _tridata **tris;
  float best_xyz[3], best_norm[3];
} vertdata;

/* --------------------------  tridata ------------------------------------ */

typedef struct _tridata {
  vertdata *verts[3];
  int examined;
  float area,center[3],norm[3];
} tridata;

/* --------------------------  tridata ------------------------------------ */

typedef struct _geomdata {
  vertdata *verts;
  tridata *tris;
  int nverts, ntris;
} geomdata;

 
/* --------------------------  get_other_vertex ------------------------------------ */

vertdata *get_other_vertex(tridata *tri, vertdata *vert1, vertdata *vert2){
  vertdata **verts;
  
  verts = tri->verts;

  if(vert1==verts[0]&&vert2==verts[1]||vert1==verts[1]&&vert2==verts[0])return verts[2];
  if(vert1==verts[0]&&vert2==verts[2]||vert1==verts[2]&&vert2==verts[0])return verts[1];
  if(vert1==verts[1]&&vert2==verts[2]||vert1==verts[2]&&vert2==verts[1])return verts[0];
  return NULL;
}
 
/* --------------------------  MakeVertRing ------------------------------------ */

void MakeVertRing(vertdata *vert){
  int i;
  vertdata **vertring;
  tridata *tri;
  int nvertring;

  if(vert->ntris==0)return;
  NewMemory((void **)&vert->vertring,vert->ntris*sizeof(vertdata *));

  for(i=0;i<vert->ntris;i++){
    tridata *trii;

    trii = vert->tris[i];
    trii->examined=0;
  }
  vertring=vert->vertring;
  nvertring=0;
  tri=vert->tris[0];
  tri->examined=1;
  if(tri->verts[0]==vert){
    vertring[nvertring++]=tri->verts[1];
    vertring[nvertring++]=tri->verts[2];
  }
  else if(tri->verts[1]==vert){
    vertring[nvertring++]=tri->verts[2];
    vertring[nvertring++]=tri->verts[0];
  }
  else{
    vertring[nvertring++]=tri->verts[0];
    vertring[nvertring++]=tri->verts[1];
  }
  for(i=1;i<vert->ntris;i++){
    int j;

    for(j=1;j<vert->ntris;j++){
      tridata *trij;
      vertdata *other;

      trij = vert->tris[j];
      if(trij->examined==1)continue;
      other=get_other_vertex(trij,vert,vertring[nvertring-1]);
      if(other==NULL)continue;
      trij->examined=1;
      vertring[nvertring++]=other;
    }
  }
  vert->nvertring=nvertring;
}

/* ------------------ getVertInfo ------------------------ */

void getVertInfo(geomdata *geom){
  int i;

  for(i=0;i<geom->nverts;i++){
    vertdata *verti;

    verti = geom->verts+i;
    verti->ntris=0;
    verti->tris=NULL;
  }
  for(i=0;i<geom->ntris;i++){
    tridata *trii;
    int j;

    trii = geom->tris+i;
    for(j=0;j<3;j++){
      vertdata *verti;

      verti = trii->verts[j];
      verti->ntris++;
    }
  }
  for(i=0;i<geom->nverts;i++){
    vertdata *verti;

    verti = geom->verts+i;
    if(verti->ntris>0)NewMemory((void **)&verti->tris,verti->ntris*sizeof(tridata *));
    verti->ntris=0;
  }
  for(i=0;i<geom->ntris;i++){
    tridata *trii;
    int j;

    trii = geom->tris+i;
    for(j=0;j<3;j++){
      vertdata *verti;

      verti = trii->verts[j];
      verti->tris[verti->ntris]=trii;
      verti->ntris++;
    }
  }
  for(i=0;i<geom->nverts;i++){
    vertdata *verti;

    verti = geom->verts+i;
    MakeVertRing(verti);  
  }
}
/* ------------------ calcNormal3 ------------------------ */

void getNormal(float *v1, float *v2, float *v3, float *area, float *normal){
  float u[3], v[3], tri_area;
  int i;

  for(i=0;i<3;i++){
    u[i]=v2[i]-v1[i];
    v[i]=v3[i]-v1[i];
  }

  normal[0] = u[1]*v[2] - u[2]*v[1];
  normal[1] = u[2]*v[0] - u[0]*v[2];
  normal[2] = u[0]*v[1] - u[1]*v[0];

  tri_area = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
  normal[0]/=tri_area;
  normal[1]/=tri_area;
  normal[2]/=tri_area;
  *area = tri_area;
}

/* --------------------------  tridata ------------------------------------ */

void GetTriInfo(tridata *tri){
  float *xyz1, *xyz2, *xyz3, *center;

  xyz1 = tri->verts[0]->xyz;
  xyz2 = tri->verts[1]->xyz;
  xyz3 = tri->verts[2]->xyz;
  center = tri->center;
  center[0] = (xyz1[0]+xyz2[0]+xyz3[0])/3.0;
  center[1] = (xyz1[1]+xyz2[1]+xyz3[1])/3.0;
  center[2] = (xyz1[2]+xyz2[2]+xyz3[2])/3.0;
  getNormal(xyz1,xyz2,xyz3,&tri->area,tri->norm);
}

/* --------------------------  get_bestplane ------------------------------------ */

void GetBestPlane(vertdata *vert){
  int i;
  float *best_norm, *best_xyz;

  best_norm = vert->best_norm;
  best_xyz = vert->best_xyz;

  for(i=0;i<3;i++){
    best_norm[i]=0.0;
    best_xyz[i]=0.0;
  }
  for(i=0;i<vert->ntris;i++){
    int j;
    tridata *tri;

    tri = vert->tris[i];
    for(j=0;j<3;j++){
      best_norm[j] += tri->area*tri->norm[j];
      best_xyz[j] += tri->area*tri->center[j];
    }
  }
  for(i=0;i<3;i++){
    best_norm[i]/=vert->ntris;
    best_xyz[i]/=vert->ntris;
  }
}
