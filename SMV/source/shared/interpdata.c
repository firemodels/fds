// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "interpdata.h"
#include "MALLOC.h"

  // svn revision character string
char interpdata_revision[]="$Revision$";

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

#define IJKRAD(i,j,k) (i) + nx*(j) + nxy*(k)
#define IJKRAD2(i,j,k) ((i)+1) + (nx+2)*((j)+1) + (nx+2)*(ny+2)*((k)+1)

void build_radiancemap2(radiancedata *radianceinfo){
  int i, j, k, nx, ny, nz, nxy;
  float *fradiance;
  unsigned char *radiance, *opacity;

  nx = radianceinfo->ijkbar[0];
  ny = radianceinfo->ijkbar[1];
  nz = radianceinfo->ijkbar[2];

  radiance = radianceinfo->radiance;
  opacity = radianceinfo->opacity;

  nxy = nx*ny;
  NewMemory((void **)&fradiance,nx*ny*nz*sizeof(float));

  i=0;
  for(j=0;j<ny;j++){
  for(k=0;k<nz;k++){
    fradiance[IJKRAD(i,j,k)]=1.0;
  }
  }

  for(i=1;i<nx;i++){
  for(j=0;j<ny;j++){
  for(k=0;k<nz;k++){
    int ijk,im1jk;

    ijk=IJKRAD(i,j,k);
    im1jk=ijk-1;
    fradiance[ijk]=(float)(fradiance[im1jk]*(float)(255-opacity[im1jk])/255.0);
  }
  }
  }

  for(i=0;i<nx*ny*nz;i++){
    radiance[i]=(unsigned char)(fradiance[i]*255.0);
  }
  FREEMEMORY(fradiance);
}


/* ----------------------- build_radiancemap ----------------------------- */

void build_radiancemap(radiancedata *radianceinfo){
  int i, j, k, nx, ny, nz, nxy;
  float *fradiance,*total_fradiance;
  unsigned char *radiance, *opacity;

  nx = radianceinfo->ijkbar[0];
  ny = radianceinfo->ijkbar[1];
  nz = radianceinfo->ijkbar[2];

  radiance = radianceinfo->radiance;
  opacity = radianceinfo->opacity;

  nxy = nx*ny;
  NewMemory((void **)&fradiance,(nx+2)*(ny+2)*(nz+2)*sizeof(float));
  NewMemory((void **)&total_fradiance,(nx+2)*(ny+2)*(nz+2)*sizeof(float));

  for(i=0;i<(nx+2)*(ny+2)*(nz+2);i++){
    total_fradiance[i]=0.0;
  }
  
  for(j=-1;j<ny+1;j++){
  for(k=-1;k<nz+1;k++){
    fradiance[IJKRAD2(-1,j,k)]=0.2;
    fradiance[IJKRAD2(nx,j,k)]=0.2;
  }
  }
  for(i=-1;i<nx+1;i++){
  for(k=-1;k<nz+1;k++){
    fradiance[IJKRAD2(i,-1,k)]=0.2;
    fradiance[IJKRAD2(i,ny,k)]=0.2;
  }
  }
  for(i=-1;i<nx+1;i++){
  for(j=-1;j<ny+1;j++){
    fradiance[IJKRAD2(i,j,nz)]=0.2;
  }
  }
  CheckMemory;

  // increasing i
  
  for(i=0;i<nx;i++){
  for(j=0;j<ny;j++){
  for(k=0;k<nz;k++){
    int ijk,fijk,fim1jk;

    ijk=IJKRAD(i,j,k);
    fijk=IJKRAD2(i,j,k);
    fim1jk=fijk-1;
    fradiance[fijk]=(float)(fradiance[fim1jk]*(float)(255-opacity[ijk])/255.0);
  }
  }
  }
  for(i=0;i<(nx+2)*(ny+2)*(nz+2);i++){
    total_fradiance[i]+=fradiance[i];
  }
  CheckMemory;

  // decreasing i
  
  for(i=nx-1;i>=0;i--){
  for(j=0;j<ny;j++){
  for(k=0;k<nz;k++){
    int ijk,fijk,fip1jk;

    ijk=IJKRAD(i,j,k);
    fijk=IJKRAD2(i,j,k);
    fip1jk=fijk+1;
    fradiance[fijk]=(float)(fradiance[fip1jk]*(float)(255-opacity[ijk])/255.0);
  }
  }
  }
  for(i=0;i<(nx+2)*(ny+2)*(nz+2);i++){
    total_fradiance[i]+=fradiance[i];
  }
  CheckMemory;

  // increasing j
  
  for(j=0;j<ny;j++){
  for(i=0;i<nx;i++){
  for(k=0;k<nz;k++){
    int ijk,fijk,fijm1k;

    ijk=IJKRAD(i,j,k);
    fijk=IJKRAD2(i,j,k);
    fijm1k=ijk-(nx+2);
    fradiance[fijk]=(float)(fradiance[fijm1k]*(float)(255-opacity[ijk])/255.0);
  }
  }
  }
  for(i=0;i<(nx+2)*(ny+2)*(nz+2);i++){
    total_fradiance[i]+=fradiance[i];
  }
  CheckMemory;

  // decreasing j
  
  for(j=ny-1;j>=0;j--){
  for(i=0;i<nx;i++){
  for(k=0;k<nz;k++){
    int ijk,fijk,fijp1k;

    ijk=IJKRAD(i,j,k);
    fijk=IJKRAD2(i,j,k);
    fijp1k=ijk+(nx+2);
    fradiance[fijk]=(float)(fradiance[fijp1k]*(float)(255-opacity[ijk])/255.0);
  }
  }
  }
  for(i=0;i<(nx+2)*(ny+2)*(nz+2);i++){
    total_fradiance[i]+=fradiance[i];
  }
  CheckMemory;
  
  // decreasing k
  
  for(k=nz-1;k>=0;k--){
  for(j=0;j<ny;j++){
  for(i=0;i<nx;i++){
    int ijk,fijk,fijkp1;

    ijk=IJKRAD(i,j,k);
    fijk=IJKRAD2(i,j,k);
    fijkp1=ijk+(nx+2)*(ny+2);
    fradiance[fijk]=(float)(fradiance[fijkp1]*(float)(255-opacity[ijk])/255.0);
  }
  }
  }
  for(i=0;i<(nx+2)*(ny+2)*(nz+2);i++){
    total_fradiance[i]+=fradiance[i];
  }

  CheckMemory;

  for(k=0;k<nz;k++){
  for(j=0;j<ny;j++){
  for(i=0;i<nx;i++){
    radiance[IJKRAD(i,j,k)]=(unsigned char)(total_fradiance[IJKRAD2(i,j,k)]*255.0);
    CheckMemory;
  }
  }
  }
  CheckMemory;
  FREEMEMORY(fradiance);
  FREEMEMORY(total_fradiance);
}

/* ----------------------- interp3d ----------------------------- */

#define INTERP1D(f0,f1,dx) (float)((f0) + ((f1)-(f0))*(dx))
float interp3d(float *xplt, float *yplt, float *zplt, int ibar, int jbar, int kbar, float *vals, float xyz[3]){
  int i, j, k;
  int ijk;
  float val000,val100,val010,val110;
  float val001,val101,val011,val111;
  float val00,val01,val10,val11;
  float val0, val1, val;
  int nx, ny, nxy;
  float dx, dy, dz;
  float dxbar, dybar, dzbar;
  float *vv;

  dxbar = xplt[1]-xplt[0];
  dybar = yplt[1]-yplt[0];
  dzbar = zplt[1]-zplt[0];

  i = (xyz[0]-xplt[0])/dxbar;
  j = (xyz[1]-yplt[0])/dybar;
  k = (xyz[2]-zplt[0])/dzbar;

  dx = (xyz[0] - i*dxbar)/dxbar;
  dy = (xyz[1] - j*dybar)/dybar;
  dz = (xyz[2] - k*dzbar)/dzbar;

  nx = ibar;
  ny = jbar;
  nxy = nx*ny;

  ijk = i + j*nx + k*nxy;

  vv = vals + ijk;
  val000 = vv[0];
  val100 = vv[1];
  val010 = vv[nx];
  val110 = vv[1+nx];
  val001 = vv[nxy];
  val101 = vv[1+nxy];
  val011 = vv[nx+nxy];
  val111 = vv[1+nx+nxy];
  val00 = INTERP1D(val000,val100,dx);
  val10 = INTERP1D(val010,val110,dx);
  val01 = INTERP1D(val001,val101,dx);
  val11 = INTERP1D(val011,val111,dx);
   val0 = INTERP1D( val00, val10,dy);
   val1 = INTERP1D( val01, val11,dy);
    val = INTERP1D(  val0,  val1,dz);

  return val;
}

/*
typedef struct {
  float xyz[3];
} point;

typedef struct _kd_data {
  struct _kd_data *left, *right, *parent;
  int axis;
  point *median;
} kd_data;
*/

#define NKDPOINTS 1000000
void test_kd(void){
  int i;
  kdpoint *points,*points2,**pointers;
  kd_data *kdtree;
  float xyztest[3];
  kd_data **bests;
  int nbests, nwanted;
  int nkdpoints=NKDPOINTS;

  points = malloc(NKDPOINTS*sizeof(kdpoint));
  points2 = malloc(NKDPOINTS*sizeof(kdpoint));
  pointers = malloc(NKDPOINTS*sizeof(kdpoint *));
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
  xyztest[0] = 0.25;
  xyztest[1] = 0.25;
  xyztest[2] = 0.25;
  
  nwanted=10;
  nbests=0;
  bests = malloc(nwanted*sizeof(kd_data *));
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

/* ----------------------- child_away ----------------------------- */

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

/* ----------------------- vdistance2 ----------------------------- */

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

/* ----------------------- closest_nodes ----------------------------- */
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
  nleft = median_index;
  nright = npoints - median_index - 1;
  kdptr->parent=parent;
  kdptr->left=setup_kdtree(points,nleft,kdptr);
  kdptr->right=setup_kdtree(points+median_index+1,nright,kdptr);
  return kdptr;
}

/* ----------------------- free_kdtree ----------------------------- */

void free_kdtree(kd_data *kdtree){
  if(kdtree->left!=NULL)free_kdtree(kdtree->left);
  if(kdtree->right!=NULL)free_kdtree(kdtree->right);
  FREEMEMORY(kdtree);
}
