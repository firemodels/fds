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

/* ------------------ Insert_Edge ------------------------ */

void Insert_Edge(edgedata *edges, int *nedges, vertdata *vert1, vertdata *vert2){
  vertdata *vertmin, *vertmax;
  int i;
  int haveit=0;

  vertmin=vert1;
  vertmax=vert1;
  if(vert2<vert1)vertmin=vert2;
  if(vert2>vert1)vertmax=vert2;
  for(i=0;i<*nedges;i++){
    edgedata *edgei;

    edgei = edges+i;
    if(edgei->verts[0]==vertmin&&edgei->verts[1]==vertmax){
      edgei->ntris++;
      vertmin->nedges++;
      vertmax->nedges++;
      haveit=1;
      break;
    }
  }
  if(haveit==0){
    edgedata *edgei;

    edgei = edges + *nedges;
    edgei->verts[0]=vertmin;
    edgei->verts[1]=vertmax;
    vertmin->nedges=1;
    vertmax->nedges=1;
    edgei->ntris=1;
    (*nedges)++;
  }
}

/* ------------------ InsertVert ------------------------ */

void InsertVert(vertdata **verts, int *nverts,vertdata *vert){
  int i;
  int have_vert=0;

  for(i=0;i<*nverts;i++){
    if(verts[i]==vert){
      have_vert=1;
      break;
    }
  }
  if(have_vert==1){
    verts[i]=vert;
    (*nverts)++;
  }
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
  
  // count triangles connected to each vertex

  for(i=0;i<geom->ntris;i++){
    tridata *trii;
    int j;

    trii = geom->tris+i;
    trii->verts[0]->ntris++;
    trii->verts[1]->ntris++;
    trii->verts[2]->ntris++;
  }

  // allocate memory for triangle pointers connected to vertices

  for(i=0;i<geom->nverts;i++){
    vertdata *verti;

    verti = geom->verts+i;
    if(verti->ntris>0){
      NewMemory((void **)&verti->tris,verti->ntris*sizeof(tridata *));
      NewMemory((void **)&verti->verts,2*verti->ntris*sizeof(vertdata *));
      NewMemory((void **)&verti->verts_temp,2*verti->ntris*sizeof(vertdata *));
    }
    verti->ntris=0;
    verti->type=GEOM_INTERIOR;
  }

  //  associate triangles connected to each vertex

  for(i=0;i<geom->ntris;i++){
    tridata *trii;
    vertdata *verti;
    int j;

    trii = geom->tris+i;
    verti = trii->verts[0];
    verti->tris[verti->ntris++]=trii;
    verti = trii->verts[1];
    verti->tris[verti->ntris++]=trii;
    verti = trii->verts[2];
    verti->tris[verti->ntris++]=trii;
  }

  // allocate memory for edges (at most 3 * number of triangles)

  if(geom->ntris>0)NewMemory((void **)&geom->edges,3*geom->ntris*sizeof(edgedata *));

  // insert triangles edge pointers into edges data structure

  geom->nedges=0;
  for(i=0;i<geom->ntris;i++){
    tridata *tri;

    tri = geom->tris+i;
    Insert_Edge(geom->edges,&geom->nedges,tri->verts[0],tri->verts[1]);
    Insert_Edge(geom->edges,&geom->nedges,tri->verts[1],tri->verts[2]);
    Insert_Edge(geom->edges,&geom->nedges,tri->verts[2],tri->verts[0]);
  }
  if(geom->nedges>0)ResizeMemory((void **)&geom->edges,geom->nedges*sizeof(edgedata *));

  // determine vertex type

  //   interior: all connected edges are connected to exactly two triangles
  //   exterior: at least one connected edge is connected to only one triangle and
  //                    all other connected edges are connected to two triangles
  //    complex: at least one edge is connected to more than two triangles

  for(i=0;i<geom->nedges;i++){
    edgedata *edgei;
    vertdata *vert1,*vert2;

    edgei = geom->edges + i;
    vert1 = edgei->verts[0];
    vert2 = edgei->verts[1];
    if(edgei->ntris==1){
      if(vert1->type==GEOM_INTERIOR)vert1->type=GEOM_EXTERIOR;
      if(vert2->type==GEOM_INTERIOR)vert2->type=GEOM_EXTERIOR;
    }
    else if(edgei->ntris>2){
      vert1->type=GEOM_COMPLEX;
      vert2->type=GEOM_COMPLEX;
    }
  }

  for(i=0;i<geom->nverts;i++){
    int j;
    vertdata *verti;

    verti = geom->verts + i;
    for(j=0;j<verti->ntris;j++){
      tridata *trij;

      trij = verti->tris[j];
      if(trij->verts[0]!=verti)InsertVert(verti->verts_temp,&(verti->nverts),trij->verts[0]);
      if(trij->verts[1]!=verti)InsertVert(verti->verts_temp,&(verti->nverts),trij->verts[1]);
      if(trij->verts[2]!=verti)InsertVert(verti->verts_temp,&(verti->nverts),trij->verts[2]);
    }

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
