// $Date$ 
// $Revision$
// $Author$

// svn revision character string
char geometry_revision[]="$Revision$";

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MALLOC.h"
#include "geometry.h"
 
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

  if(vert->ntrifan==0)return;
  NewMemory((void **)&vert->vertring,vert->ntrifan*sizeof(vertdata *));

  for(i=0;i<vert->ntrifan;i++){
    tridata *trii;

    trii = vert->trifan[i];
    trii->examined=0;
  }
  vertring=vert->vertring;
  nvertring=0;
  tri=vert->trifan[0];
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
  for(i=1;i<vert->ntrifan;i++){
    int j;

    for(j=1;j<vert->ntrifan;j++){
      tridata *trij;
      vertdata *other;

      trij = vert->trifan[j];
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

void Insert_Edge(edgedata *edges, int *nedges, vertdata *vert1, vertdata *vert2){
  vertdata *vertmin, *vertmax;
  int i;
  int haveit=0;

  vertmin=vert1;
  vertmax=vert2;
  if(vert2<vert1)vertmin=vert2;
  if(vert2>vert1)vertmax=vert2;
  for(i=0;i<*nedges;i++){
    edgedata *edgei;

    edgei = edges+i;
    if(edgei->verts[0]==vertmin&&edgei->verts[1]==vertmax){
      haveit=1;
      break;
    }
  }
  if(haveit==0){
    edgedata *edgei;

    edgei = edges + *nedges;
    edgei->verts[0]=vertmin;
    edgei->verts[1]=vertmax;
    (*nedges)++;
  }
}

/* ------------------ getVertInfo ------------------------ */

void getVertInfo(geomdata *geom){
  int i;

  for(i=0;i<geom->nverts;i++){
    vertdata *verti;

    verti = geom->verts+i;
    verti->ntrifan=0;
    verti->trifan=NULL;
  }
  
  // count triangles connected to each vertex

  for(i=0;i<geom->ntris;i++){
    tridata *trii;
    int j;

    trii = geom->tris+i;
    trii->verts[0]->ntrifan++;
    trii->verts[1]->ntrifan++;
    trii->verts[2]->ntrifan++;
  }
  for(i=0;i<geom->nverts;i++){
    vertdata *verti;

    verti = geom->verts+i;
    if(verti->ntrifan>0)NewMemory((void **)&verti->trifan,verti->ntrifan*sizeof(tridata *));
    verti->ntrifan=0;
  }

  //  assign triangle connected to each vertex

  for(i=0;i<geom->ntris;i++){
    tridata *trii;
    vertdata *verti;
    int j;

    trii = geom->tris+i;
    verti = trii->verts[0];
    verti->trifan[verti->ntrifan++]=trii;
    verti = trii->verts[1];
    verti->trifan[verti->ntrifan++]=trii;
    verti = trii->verts[2];
    verti->trifan[verti->ntrifan++]=trii;
  }
  for(i=0;i<geom->nverts;i++){
    vertdata *verti;

    verti = geom->verts+i;
    MakeVertRing(verti);  
  }
  if(geom->ntris>0)NewMemory((void **)&geom->edges,3*geom->ntris*sizeof(edgedata *));
  geom->nedges=0;
  for(i=0;i<geom->ntris;i++){
    tridata *tri;

    tri = geom->tris+i;
    Insert_Edge(geom->edges,&geom->nedges,tri->verts[0],tri->verts[1]);
    Insert_Edge(geom->edges,&geom->nedges,tri->verts[1],tri->verts[2]);
    Insert_Edge(geom->edges,&geom->nedges,tri->verts[2],tri->verts[0]);
  }
  if(geom->nedges>0)ResizeMemory((void **)&geom->edges,geom->nedges*sizeof(edgedata *));
  for(i=0;i<geom->nedges;i++){
    edgedata *edgei;

    edgei = geom->edges+i;
    edgei->ntris=0;
  }
  for(i=0;i<geom->ntris;i++){
    tridata *trii;

  }
}
/* ------------------ calcNormal3 ------------------------ */

void getNormal(float *v1, float *v2, float *v3, float *area, float *normal){
  float u[3], v[3], tri_area;
  int i;

  for(i=0;i<3;i++){
    u[i]=v3[i]-v2[i];
    v[i]=v1[i]-v2[i];
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
  for(i=0;i<vert->ntrifan;i++){
    int j;
    tridata *tri;

    tri = vert->trifan[i];
    for(j=0;j<3;j++){
      best_norm[j] += tri->area*tri->norm[j];
      best_xyz[j] += tri->area*tri->center[j];
    }
  }
  for(i=0;i<3;i++){
    best_norm[i]/=vert->ntrifan;
    best_xyz[i]/=vert->ntrifan;
  }
}
