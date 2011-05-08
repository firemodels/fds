// $Date$ 
// $Revision$
// $Author$

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
  float x, y, z;
} point;

typedef struct _kd_data {
  struct _kd_data *left, *right;
  point *median;
} kd_data;
*/

/* ----------------------- compare_pointx ----------------------------- */

int compare_pointx( const void *arg1, const void *arg2 ){
  point *pointi, *pointj;

  pointi = *(point **)arg1;
  pointj = *(point **)arg2;

  if(pointi->xyz[0]<pointj->xyz[0])return -1;
  if(pointi->xyz[0]>pointj->xyz[0])return 1;
  return 0;
}

/* ----------------------- compare_pointy ----------------------------- */

int compare_pointy( const void *arg1, const void *arg2 ){
  point *pointi, *pointj;

  pointi = *(point **)arg1;
  pointj = *(point **)arg2;

  if(pointi->xyz[1]<pointj->xyz[1])return -1;
  if(pointi->xyz[1]>pointj->xyz[1])return 1;
  return 0;
}

/* ----------------------- compare_pointz ----------------------------- */

int compare_pointz( const void *arg1, const void *arg2 ){
  point *pointi, *pointj;

  pointi = *(point **)arg1;
  pointj = *(point **)arg2;

  if(pointi->xyz[2]<pointj->xyz[2])return -1;
  if(pointi->xyz[2]>pointj->xyz[2])return 1;
  return 0;
}

#define DIST2(result,xyz1,xyz2) \
  dx = xyz1[0]-xyz2[0];\
  dy = xyz1[1]-xyz2[1];\
  dz = xyz1[2]-xyz2[2];\
  result = dx*dx+dy*dy+dz*dz

/* ----------------------- closest_node_candidate ----------------------------- */

kd_data *closest_node_candidate(kd_data *kdtree, float *xyz, float *dist2){
  kd_data *cn;
  float dx, dy, dz;

  cn = kdtree;
  while(cn->left!=NULL&&cn->right!=NULL){
    int axis;

    axis = cn->axis;
    if(cn->right==NULL||xyz[axis]<cn->median->xyz[axis]){
      cn=cn->left;
    }
    else{
      cn=cn->right;
    }
  }
  DIST2(*dist2,cn->median->xyz,xyz);
  return cn;
}

/*
function kdsearchnn( here, point, best )
 
    if here == nil then
        return best
    end
 
    if best == nil then
        best = here
    end
 
    -- consider the current node --
    if distance(here,point) < distance(best,point) then
        best = here
    end
 
    -- search the near branch --
    child = child_near(here,point)
    best = kdsearchnn( child, point, best )
 
    -- search the away branch - maybe --
    -- (note that the following test should be <= if you let points equal to the median lie on either side --
    -- when the tree was being built) --
    if distance_axis(here,point) < distance(best,point) then
        child = child_away(here,point)
        best = kdsearchnn( child, point, best )
    end
 
    return best
 
end*/
/* ----------------------- closest_node ----------------------------- */

kd_data *closest_node(kd_data *kdtree, float *xyz){
  kd_data *node,*current_best,*parent, *other,*other_best;
  float current_distance2, other_distance2, split_dist2;

  current_best = closest_node_candidate(kdtree, xyz, &current_distance2);

  node = current_best;
  while(node->parent!=NULL){
    float node_dist2;
    float dx,dy,dz;

    DIST2(node_dist2,node->median->xyz,xyz);
    if(node_dist2<current_distance2){
      current_distance2 = node_dist2;
      current_best = node;
    }
    //difference between the splitting coordinate of the search point and current node 
    // is less than 
    // the distance (overall coordinates) from the search point to the current best.
    split_dist2 = node->median[node->axis]-xyz[node->axis];
    split_dist2 *= split_dist2;


    node = node->parent;
  }


  other_best = closest_node_candidate(other,xyz,&other_distance2);
  return current_best;
}

/* ----------------------- setup_kdtree ----------------------------- */

kd_data *setup_kdtree(point *points, int npoints, kd_data *parent){
  int axis;
  kd_data *kdptr;
  int median_index,nleft,nright;
  float xyzmin[3], xyzmax[3];
  float dx, dy, dz;
  int i;

  if(npoints<=0)return NULL;
  NewMemory((void **)&kdptr,sizeof(kd_data));
  xyzmin[0]=points[0].xyz[0];
  xyzmax[0]=xyzmin[0];

  xyzmin[1]=points[0].xyz[1];
  xyzmax[1]=xyzmin[1];

  xyzmin[2]=points[2].xyz[2];
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
      qsort((point *)points,npoints,sizeof(kd_data),compare_pointx);
      break;
    case 1:
      qsort((point *)points,npoints,sizeof(kd_data),compare_pointy);
      break;
    case 2:
      qsort((point *)points,npoints,sizeof(kd_data),compare_pointz);
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
