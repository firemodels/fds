// $Date$ 
// $Revision$
// $Author$
#define INKDTREE

#include "options.h"
#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "kdtree.h"
#include "MALLOC.h"

  // svn revision character string
char kdtree_revision[]="$Revision$";
int ntotal=0;
int nsort=1;

#define NKDPOINTS 1000000

/* ----------------------- test_kd ----------------------------- */

void test_kd(void){
  int i;
  kdpoint *points,*points2,**pointers;
  kd_data *kdtree;
  float xyztest[3];
  kd_data **bests;
  int nbests, nwanted;
  int nkdpoints=NKDPOINTS;

  NewMemory((void **)&points,NKDPOINTS*sizeof(kdpoint));
  NewMemory((void **)&points2,NKDPOINTS*sizeof(kdpoint));
  NewMemory((void **)&pointers,NKDPOINTS*sizeof(kdpoint *));
  for(i=0;i<NKDPOINTS;i++){
    points[i].xyz[0] = (float)rand()/RAND_MAX;
    points[i].xyz[1] = (float)rand()/RAND_MAX;
    points[i].xyz[2] = (float)rand()/RAND_MAX;
    points2[i].xyz[0] = points[i].xyz[0];
    points2[i].xyz[1] = points[i].xyz[1];
    points2[i].xyz[2] = points[i].xyz[2];
    pointers[i] = points2 + i;
  }
  kdtree = setup_kdtree(points,NKDPOINTS, NULL);
  printf("ntotal=%i nsort=%i\n",ntotal,nsort);
  xyztest[0] = 0.25;
  xyztest[1] = 0.25;
  xyztest[2] = 0.25;
  
  nwanted=10;
  nbests=0;
  NewMemory((void **)&bests,nwanted*sizeof(kd_data *));
  get_closest_nodes(kdtree, xyztest, bests, &nbests, nwanted);
  sort_closest_nodes(bests, nbests, xyztest);
  get_closest_points(pointers, nkdpoints, xyztest);
  for(i=0;i<nwanted;i++){
    kd_data *best;
    float *xyz,*xyzp;

    best = bests[i];
    xyz = best->median->xyz;
    xyzp = pointers[i]->xyz;
    printf("x1=%f y1=%f z1=%f dist=%f\n",xyz[0],xyz[1],xyz[2],best->median->dist2);
    printf("x2=%f y2=%f z2=%f dist=%f\n",xyzp[0],xyzp[1],xyzp[2],pointers[i]->dist2);
    printf("\n");
  }
}

/* ----------------------- compare_bests ----------------------------- */

int compare_bests( const void *arg1, const void *arg2 ){
  kd_data *besti, *bestj;

  besti = *(kd_data **)arg1;
  bestj = *(kd_data **)arg2;

  if(besti->median->dist2<bestj->median->dist2)return -1;
  if(besti->median->dist2>bestj->median->dist2)return 1;
  return 0;
}

/* ----------------------- compare_points ----------------------------- */

int compare_points( const void *arg1, const void *arg2 ){
  kdpoint *pointi, *pointj;

  pointi = *(kdpoint **)arg1;
  pointj = *(kdpoint **)arg2;

  if(pointi->dist2<pointj->dist2)return -1;
  if(pointi->dist2>pointj->dist2)return 1;
  return 0;
}

/* ----------------------- compare_pointx ----------------------------- */

int compare_pointx( const void *arg1, const void *arg2 ){
  kdpoint *pointi, *pointj;

  pointi = (kdpoint *)arg1;
  pointj = (kdpoint *)arg2;

  if(pointi->xyz[0]<pointj->xyz[0])return -1;
  if(pointi->xyz[0]>pointj->xyz[0])return 1;
  return 0;
}

/* ----------------------- compare_pointy ----------------------------- */

int compare_pointy( const void *arg1, const void *arg2 ){
  kdpoint *pointi, *pointj;

  pointi = (kdpoint *)arg1;
  pointj = (kdpoint *)arg2;

  if(pointi->xyz[1]<pointj->xyz[1])return -1;
  if(pointi->xyz[1]>pointj->xyz[1])return 1;
  return 0;
}

/* ----------------------- compare_pointz ----------------------------- */

int compare_pointz( const void *arg1, const void *arg2 ){
  kdpoint *pointi, *pointj;

  pointi = (kdpoint *)arg1;
  pointj = (kdpoint *)arg2;

  if(pointi->xyz[2]<pointj->xyz[2])return -1;
  if(pointi->xyz[2]>pointj->xyz[2])return 1;
  return 0;
}

/* ----------------------- distance_axis ----------------------------- */

float distance_axis(kd_data *node, float *xyz){
  float daxis, *xyzm;
  int axis;

  axis = node->axis;
  xyzm = node->median->xyz;

  daxis = xyzm[axis]-xyz[axis];
  return daxis*daxis;
}

/* ----------------------- distance2 ----------------------------- */

float distance2(kd_data *node, float *xyz){
  float dx, dy, dz, *xyzm;

  xyzm = node->median->xyz;
  dx = xyzm[0]-xyz[0];
  dy = xyzm[1]-xyz[1];
  dz = xyzm[2]-xyz[2];
  return dx*dx + dy*dy + dz*dz;
}

/* ----------------------- child_near ----------------------------- */

kd_data *child_near(kd_data *here, float *point){
  int axis;

  axis = here->axis;
  if(point[axis]<here->median->xyz[axis]){
    return here->left;
  }
  else{
    return here->right;
  }
}

/* ----------------------- child_far ----------------------------- */

kd_data *child_far(kd_data *here, float *point){
  int axis;

  axis = here->axis;
  if(point[axis]<here->median->xyz[axis]){
    return here->right;
  }
  else{
    return here->left;
  }
}

/* ----------------------- maxdistance ----------------------------- */

float maxdistance2(kd_data **bests,int nbests, float *xyz,int *nmax){
  int i;
  float maxdist2=-1.0;

  for(i=0;i<nbests;i++){
    float dist2;

    dist2 = distance2(bests[i],xyz);
    if(dist2>maxdist2){
      maxdist2=dist2;
      *nmax=i;
    }
  }
  return maxdist2;
}

/* ----------------------- get_closest_points ----------------------------- */

void get_closest_points(kdpoint **pointers, int npoints, float *point){
  int i;

  for(i=0;i<npoints;i++){
    float dx, dy, dz;
    float *xyz;

    xyz = pointers[i]->xyz;

    dx = xyz[0]-point[0];
    dy = xyz[1]-point[1];
    dz = xyz[2]-point[2];
    pointers[i]->dist2=dx*dx + dy*dy + dz*dz;
  }
  qsort((kdpoint **)pointers,npoints,sizeof(kdpoint *),compare_points);
}

/* ----------------------- sort_closest_nodes ----------------------------- */

void sort_closest_nodes(kd_data **bests, int nbests, float *point){
  int i;

  for(i=0;i<nbests;i++){
    float dx, dy, dz, *xyz;
    kd_data *best;

    best = bests[i];
    xyz = best->median->xyz;
    dx = xyz[0]-point[0];
    dy = xyz[1]-point[1];
    dz = xyz[2]-point[2];
    best->median->dist2=dx*dx + dy*dy + dz*dz;
  }
  qsort(bests,nbests,sizeof(kd_data *),compare_bests);
}

/* ----------------------- get_closest_nodes ----------------------------- */

void get_closest_nodes(kd_data *here, float *point, kd_data **bests, int *nbests, int nwanted){
  kd_data *child;
  int nmax;

  if(here==NULL)return;
  if(*nbests<nwanted){
    bests[*nbests]=here;
    (*nbests)++;
  }
  else if(distance2(here,point)<maxdistance2(bests,*nbests,point,&nmax)){
    bests[nmax]=here;
  }
  child = child_near(here,point);
  get_closest_nodes(child,point,bests,nbests,nwanted);
  if(distance_axis(here,point) < maxdistance2(bests,*nbests,point,&nmax)){
    child = child_far(here,point);
    get_closest_nodes(child,point,bests,nbests,nwanted);
  }
}

/* ----------------------- closest_node ----------------------------- */

kd_data *closest_node(kd_data *here, float *point, kd_data *best){
  kd_data *child;

  if(here==NULL){
    return best;
  }
  if(best==NULL||distance2(here,point)<distance2(best,point)){
    best=here;
  }
  child = child_near(here,point);
  best = closest_node(child,point,best);
  if(distance_axis(here,point) < distance2(best,point)){
    child = child_far(here,point);
    best = closest_node(child,point,best);
  }
  return best;
}

/* ----------------------- setup_kdtree ----------------------------- */

kd_data *setup_kdtree(kdpoint *points, int npoints, kd_data *parent){
  int axis;
  kd_data *kdptr;
  int median_index,nleft,nright;
  float xyzmin[3], xyzmax[3];
  float dx, dy, dz;
  int i;

  if(npoints<=0)return NULL;
  ntotal++;
  NewMemory((void **)&kdptr,sizeof(kd_data));

  xyzmin[0]=points[0].xyz[0];
  xyzmin[1]=points[0].xyz[1];
  xyzmin[2]=points[0].xyz[2];

  xyzmax[0]=xyzmin[0];
  xyzmax[1]=xyzmin[1];
  xyzmax[2]=xyzmin[2];

  for(i=1;i<npoints;i++){
    if(points[i].xyz[0]<xyzmin[0])xyzmin[0]=points[i].xyz[0];
    if(points[i].xyz[0]>xyzmax[0])xyzmax[0]=points[i].xyz[0];
    if(points[i].xyz[1]<xyzmin[1])xyzmin[1]=points[i].xyz[1];
    if(points[i].xyz[1]>xyzmax[1])xyzmax[1]=points[i].xyz[1];
    if(points[i].xyz[2]<xyzmin[2])xyzmin[2]=points[i].xyz[2];
    if(points[i].xyz[2]>xyzmax[2])xyzmax[2]=points[i].xyz[2];
  }
  dx = xyzmax[0] - xyzmin[0];
  dy = xyzmax[1] - xyzmin[1];
  dz = xyzmax[2] - xyzmin[2];

  if(dx>=dy&&dx>=dy)axis = 0;
  if(dy>=dx&&dy>=dz)axis = 1;
  if(dz>=dy&&dz>=dx)axis = 2;

  if(parent!=NULL&&axis!=parent->axis)nsort++;
  switch (axis) {
    case 0:
      qsort((kdpoint *)points,npoints,sizeof(kdpoint),compare_pointx);
      break;
    case 1:
      qsort((kdpoint *)points,npoints,sizeof(kdpoint),compare_pointy);
      break;
    case 2:
      qsort((kdpoint *)points,npoints,sizeof(kdpoint),compare_pointz);
      break;
  }
  median_index = npoints/2;
  kdptr->axis=axis;
  kdptr->median=points + median_index;
  kdptr->parent=parent;
  kdptr->left=NULL;
  kdptr->right=NULL;
  nleft = median_index;
  if(nleft>0){
    kdptr->left=setup_kdtree(points,nleft,kdptr);
  }
  nright = npoints - median_index - 1;
  if(nright>0){
    kdptr->right=setup_kdtree(points+median_index+1,nright,kdptr);
  }
  return kdptr;
}

/* ----------------------- free_kdtree ----------------------------- */

void free_kdtree(kd_data *kdtree){
  if(kdtree->left!=NULL)free_kdtree(kdtree->left);
  if(kdtree->right!=NULL)free_kdtree(kdtree->right);
  FREEMEMORY(kdtree);
}
