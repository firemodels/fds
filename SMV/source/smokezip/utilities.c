#include "options.h"
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "zlib.h"
#include "svzip.h"
#include "MALLOC.h"
#include "datadefs.h"

int iseed=0;

  /* ------------------ SmoothLabel ------------------------ */

void SmoothLabel(float *a, float *b, int n){
  float delta, factor, logdelta;
  int ndigits;

  delta = (*b-*a)/(n-2);
  if(delta==0.0)return;
  logdelta = log10((double)delta);
  ndigits=logdelta-1;
  if(logdelta<=1)ndigits--;
  factor = 5*pow(10,ndigits);
  delta = (int)(delta/factor + 0.5f)*factor;

  *a = factor*(int)(*a/factor+0.5f);
  *b = *a + (n-2)*delta;

}

/* ------------------ calcNormal2 ------------------------ */

void Normal(unsigned short *v1, unsigned short *v2, unsigned short *v3, float *normal, float *area){
  float u[3], v[3];

  float norm2;

  u[0]=v2[0]-v1[0];
  u[1]=v2[1]-v1[1];
  u[2]=v2[2]-v1[2];

  v[0]=v3[0]-v1[0];
  v[1]=v3[1]-v1[1];
  v[2]=v3[2]-v1[2];

  normal[0] = u[1]*v[2] - u[2]*v[1];
  normal[1] = u[2]*v[0] - u[0]*v[2];
  normal[2] = u[0]*v[1] - u[1]*v[0];

  norm2 = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
  *area = norm2/2.0;

  normal[0]/=norm2;
  normal[1]/=norm2;
  normal[2]/=norm2;
}

/* ------------------ atan3 ------------------------ */

float atan3(float dy,float dx){
  if(dx!=0.0)return atan(dy/dx);

  // dx is zero so atan(dy/dx) is PI/2 or -PI/2 depending on sign dy

  if(dy>0.0)return 2.0*atan(1.0);
  if(dy<0.0)return -2.0*atan(1.0);
  return 0.0;
}


/* ------------------ RandABsdir ------------------------ */

void RandABsdir(float xyz[3], int dir){
  float x=1.0, y=1.0, z=1.0;
  float sum;

  sum=x*x+y*y+z*z;
  while(sum>1.0||sum==0.0){
    x = rand_1d(0.0,1.0);
    y = rand_1d(0.0,1.0);
    z = rand_1d(0.0,1.0);
    sum=x*x+y*y+z*z;
  }
  xyz[0]=x/sqrt(sum);
  xyz[1]=y/sqrt(sum);
  xyz[2]=z/sqrt(sum);
  if(abs(dir)>=1&&abs(dir)<=3){
    if(dir>0){
      xyz[dir]=abs(xyz[dir]);
    }
    else{
      xyz[-dir]=-abs(xyz[-dir]);
    }
  }
}

/* ------------------ rand_dir ------------------------ */

void rand_cone_dir(float xyz[3], float conedir[3], float mincosangle){
  float cosangle=2.0;

  while(cosangle<mincosangle){
    rand_sphere_dir(xyz);
    cosangle = xyz[0]*conedir[0]+xyz[1]*conedir[1]+xyz[2]*conedir[2];
  }

  return;
}
/* ------------------ rand_dir ------------------------ */

void rand_sphere_dir(float xyz[3]){
  float x=1.0, y=1.0, z=1.0;
  float sum,sqsum;

  sum=x*x+y*y+z*z;
  while(sum>1.0||sum==0.0){
    x = rand_1d(-1.0,1.0);
    y = rand_1d(-1.0,1.0);
    z = rand_1d(-1.0,1.0);
    sum=x*x+y*y+z*z;
  }
  sqsum=sqrt(sum);
  xyz[0]=x/sqsum;
  xyz[1]=y/sqsum;
  xyz[2]=z/sqsum;
}

/* ------------------ rand_1d ------------------------ */

float rand_1d(float xmin, float xmax){
  float val;

  if(iseed==0){
    iseed=1;
    srand(iseed);
  }

  val = xmin + (xmax-xmin)*(float)rand()/(float)RAND_MAX;
  return val;
}

/* ------------------ rand_2d ------------------------ */

void rand_2d(float xy[2], float xmin, float xmax, float ymin, float ymax){
  xy[0]=rand_1d(xmin,xmax);
  xy[1]=rand_1d(ymin,ymax);
}

/* ------------------ rand_3d ------------------------ */

void rand_3d(float xyz[3], float xmin, float xmax, float ymin, float ymax, float zmin, float zmax){
  xyz[0]=rand_1d(xmin,xmax);
  xyz[1]=rand_1d(ymin,ymax);
  xyz[2]=rand_1d(zmin,zmax);
}
