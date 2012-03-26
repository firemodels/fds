// $Date$ 
// $Revision$
// $Author$

// svn revision character string
char geometry_revision[]="$Revision$";

#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MALLOC.h"
#include "geometry.h"
#include "datadefs.h"

/* ------------------ normalize_vec3 ------------------------ */

void normalize_vec3(float *xyz){
  float norm2;
  int i;

  norm2 = 0.0;

  for(i=0;i<3;i++){
    norm2 += xyz[i]*xyz[i];
  }
  norm2=sqrt(norm2);
  if(norm2==0.0)norm2=1.0;
  for(i=0;i<3;i++){
    xyz[i]/=norm2;
  }
}

/* ------------------ InsertEdge ------------------------ */

void InsertEdge(edgedata *edges, int *nedges, vertdata *vert1, vertdata *vert2){
  vertdata *vertmin, *vertmax;
  int i;
  edgedata *edgei;

  vertmin=MIN(vert1,vert2);
  vertmax=MAX(vert1,vert2);
  for(i=0;i<*nedges;i++){
    edgei = edges+i;
    if(edgei->verts[0]==vertmin&&edgei->verts[1]==vertmax){
      edgei->ntris++;
      vertmin->nedges++;
      vertmax->nedges++;
      return;
    }
  }
  edgei = edges + *nedges;
  edgei->verts[0]=vertmin;
  edgei->verts[1]=vertmax;
  vertmin->nedges=1;
  vertmax->nedges=1;
  edgei->ntris=1;
  (*nedges)++;
}

/* ------------------ InsertVert ------------------------ */

void InsertVert(vertdata **verts, int *nverts,vertdata *vert){
  int i;
  int have_vert=0;

  for(i=0;i<*nverts;i++){
    if(verts[i]==vert)return;
  }
  verts[*nverts]=vert;
  (*nverts)++;
}

/* ------------------ getnextvert ------------------------ */

vertdata *getnexttrivert(tridata *tri, vertdata *vert1, vertdata *vert2){
  int this_index=-1,next_index=-1;
  int i;

  for(i=0;i<3;i++){
    if(tri->verts[i]==vert1){
      this_index=i;
      next_index=i+1;
      if(next_index==3)next_index=0;
      break;
    }
  }
  if(this_index==-1)return NULL;
  if(vert2==NULL)return tri->verts[next_index];
  if(vert2!=tri->verts[next_index])return NULL;
  next_index++;
  if(next_index==3)next_index=0;
  return tri->verts[next_index];
}

/* ------------------ InsertNextVert ------------------------ */

int InsertNextVert(vertdata *verti){
  int i;
  tridata *tri;
  vertdata *lastvert;

  if(verti->ntris==0||verti->type!=GEOM_INTERIOR)return 0;
  if(verti->nverts==0){
    for(i=0;i<verti->ntris;i++){
      tri = verti->tris[i];
      tri->examined=0;
    }
    tri = verti->tris[0];
    verti->verts[verti->nverts]=getnexttrivert(tri,verti,NULL);
    verti->nverts++;
    tri->examined=1;
    return 1;
  }
  lastvert=verti->verts[verti->nverts-1];
  for(i=0;i<verti->ntris;i++){
    vertdata *nextvert;

    tri = verti->tris[i];
    if(tri->examined==1)continue;
    nextvert=getnexttrivert(tri,verti,lastvert);
    if(nextvert==NULL)continue;
    verti->verts[verti->nverts]=nextvert;
    verti->nverts++;
    tri->examined=1;
    return 1;
  }
  return 0;
}

/* ------------------ getVertInfo ------------------------ */

void getVertInfo(geomdata *geom){
  int i;

  for(i=0;i<geom->nverts;i++){
    vertdata *verti;

    verti = geom->verts+i;
    verti->ntris=0;
    verti->tris=NULL;
    verti->nedges=0;
    verti->edges=NULL;
    verti->nverts=0;
    verti->verts=NULL;
    verti->type=GEOM_INTERIOR;
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
    }
    verti->ntris=0;
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
    InsertEdge(geom->edges,&geom->nedges,tri->verts[0],tri->verts[1]);
    InsertEdge(geom->edges,&geom->nedges,tri->verts[1],tri->verts[2]);
    InsertEdge(geom->edges,&geom->nedges,tri->verts[2],tri->verts[0]);
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

  // generate ordered list of connected vertices for each vertex

  for(i=0;i<geom->nverts;i++){
    int j;
    vertdata *verti;

    verti = geom->verts + i;
    for(;;){
      if(InsertNextVert(verti)==0)break;
    }
    if(verti->nverts>0)ResizeMemory((void **)&verti->verts,verti->nverts*sizeof(vertdata *));
  }
}
/* ------------------ calcNormal3 ------------------------ */

void getNormal(float *v1, float *v2, float *v3, float *area, float *normal){
  float u[3], v[3], norm;
  int i;

  for(i=0;i<3;i++){
    u[i]=v3[i]-v2[i];
    v[i]=v1[i]-v2[i];
  }

  normal[0] = u[1]*v[2] - u[2]*v[1];
  normal[1] = u[2]*v[0] - u[0]*v[2];
  normal[2] = u[0]*v[1] - u[1]*v[0];

  norm = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
  normal[0]/=norm;
  normal[1]/=norm;
  normal[2]/=norm;
  *area = norm/2.0;
}

/* --------------------------  GetTriInfo ------------------------------------ */

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

/* --------------------------  VertOffset ------------------------------------ */

float VertOffset(vertdata *vert){
  // find distance between vertex and "best" plane that passes through the vertices
  //      connected to this vertex
  int i;
  float best_norm[3], best_xyz[3], delta[3];
  float sum;

  // find best plane

  for(i=0;i<3;i++){
    best_norm[i]=0.0;
    best_xyz[i]=0.0;
  }
  sum=0.0;
  for(i=0;i<vert->ntris;i++){
    int j;
    tridata *tri;

    tri = vert->tris[i];
    sum+=tri->area;
    for(j=0;j<3;j++){
      best_norm[j] += tri->area*tri->norm[j];
      best_xyz[j] += tri->area*tri->center[j];
    }
  }
  for(i=0;i<3;i++){
    best_xyz[i]/=sum;
  }
  normalize_vec3(best_norm);

  // best plane: xyz pts satisfying (xyz-best_xyz).dot.best_norm = 0
  // find distance between plane and vertex coordinate

  sum=0.0;
  for(i=0;i<3;i++){
    delta[i] = best_xyz[i] - vert->xyz[i];
    sum += (delta[i]*best_norm[i]);
  }
  sum = ABS(sum);
  return sum;
}

/* --------------------------  InitGeom ------------------------------------ */

void InitGeom(geomdata *geomi, float *xyz, int nxyz, int *tris, int ntris){
  int i;
  vertdata *verti;
  tridata *trii;

  geomi->verts=NULL;
  geomi->tris=NULL;
  geomi->nverts=nxyz;
  geomi->ntris=nxyz;
  if(geomi->nverts>0)NewMemory((void **)&geomi->verts,(geomi->nverts+2)*sizeof(vertdata));
  if(geomi->ntris>0)NewMemory((void **)&geomi->tris,(geomi->ntris+2)*sizeof(tridata));

  geomi->verts++;
  geomi->first_vert=geomi->verts-1;
  geomi->last_vert=geomi->verts+geomi->nverts;
  for(i=0;i<geomi->nverts;i++){
    verti = geomi->verts + i;
    verti->next=verti+1;
    verti->prev=verti-1;
    verti->xyz[0]=xyz[3*i];
    verti->xyz[1]=xyz[3*i+1];
    verti->xyz[2]=xyz[3*i+2];
    verti->active=1;
  }
  verti=geomi->verts-1;
  verti->prev=NULL;
  verti->next=verti+1;
  verti=geomi->verts+geomi->nverts;
  verti->prev=verti-1;
  verti->next=NULL;

  geomi->tris++;
  geomi->first_tri=geomi->tris-1;
  geomi->last_tri=geomi->tris+geomi->ntris+1;
  for(i=0;i<geomi->ntris;i++){
    trii = geomi->tris+i;
    trii->verts[0]=geomi->verts+tris[3*i]-1;
    trii->verts[1]=geomi->verts+tris[3*i+1]-1;
    trii->verts[2]=geomi->verts+tris[3*i+2]-1;
    trii->prev=trii-1;
    trii->next=trii+1;
    trii->active=1;
  }
  trii=geomi->tris-1;
  trii->prev=NULL;
  trii->next=trii+1;
  trii=geomi->tris+geomi->nverts;
  trii->prev=trii-1;
  trii->next=NULL;
}

/* --------------------------  DeleteTri ------------------------------------ */

void RemoveTri(tridata *trii){
  tridata *prev, *next;

  prev=trii->prev;
  next=trii->next;
  prev->next=next;
  next->prev=prev;
  trii->active=0;
}

/* --------------------------  DeleteVert ------------------------------------ */

void RemoveVert(vertdata *verti){
  int i;
  vertdata *prev,*next;

  for(i=0;i<verti->ntris;i++){
    tridata *trii;

    trii = verti->tris[i];
    RemoveTri(trii);
  }
  prev=verti->prev;
  next=verti->next;
  prev->next=next;
  next->prev=prev;
  verti->active=0;
}

void Retriangulate(vertdata **verts, int nverts, tridata **tris){
}

/* --------------------------  DecimateMesh ------------------------------------ */

void DecimateMesh(geomdata *geomi, float delta){
  int i;
  float vert_offset_max;

  vert_offset_max=0.0;
  while(vert_offset_max>delta){
    vert_offset_max=0.0;
    for(i=0;i<geomi->nverts;i++){
      float vert_offset;
      vertdata *verti;

      verti = geomi->verts + i;
      if(verti->tris!=GEOM_INTERIOR)continue;
      vert_offset = VertOffset(verti);

      vert_offset_max=MAX(vert_offset_max,vert_offset);
      if(vert_offset>delta){
        RemoveVert(verti);
        Retriangulate(verti->verts,verti->nverts,verti->tris);

      }
    }
  }
}
