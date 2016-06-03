#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MALLOC.h"
#include "datadefs.h"

/* ------------------ rotateu2v ------------------------ */

void rotateu2v(float *u, float *v, float *axis, float *angle){
  float sum,cosangle,normu,normv;

  /*
  i  j  k
  ux uy uz
  vx vy vz
  */

  CROSS(axis,u,v);
  sum = NORM3(axis);
  normu = NORM3(u);
  normv = NORM3(v);
  if(sum>0.0&&normu>0.0&&normv>0.0){
    axis[0]/=sum;
    axis[1]/=sum;
    axis[2]/=sum;
    cosangle = CLAMP(DOT3(u,v)/(normu*normv),-1.0,1.0);
    *angle=acos(cosangle);
  }
  else{
    axis[0]=0.0;
    axis[1]=0.0;
    axis[2]=1.0;
    *angle=0.0;
  }
}

/* ------------------ angleaxis2quat ------------------------ */

void angleaxis2quat(float angle, float *axis, float *quat){
  float sum;
  float cosang, sinang;

  // angle is in radians
  // axis is a vector

  sum = sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]);

  if(sum>0.0){
    cosang = cos(angle/2.0);
    sinang = sin(angle/2.0);

    quat[0] = cosang;
    quat[1] = axis[0]*sinang/sum;
    quat[2] = axis[1]*sinang/sum;
    quat[3] = axis[2]*sinang/sum;
  }
  else{
    quat[0] = 1.0;
    quat[1] = 0.0;
    quat[2] = 0.0;
    quat[3] = 0.0;
  }
}

/* ------------------ quat2rot------------------ */

void quat2rot(float quat[4],float rot[16]){
  float w, x, y, z,sum;

  sum = sqrt(quat[0]*quat[0]+quat[1]*quat[1]+quat[2]*quat[2]+quat[3]*quat[3]);
  w = quat[0]/sum;
  x = quat[1]/sum;
  y = quat[2]/sum;
  z = quat[3]/sum;

  rot[0] = 1.0 - 2.0*y*y - 2.0*z*z;
  rot[1] = 2.0*x*y + 2.0*w*z;
  rot[2] = 2.0*x*z - 2.0*w*y;
  rot[3] = 0.0;

  rot[4] = 2.0*x*y - 2.0*w*z;
  rot[5] = 1.0 - 2.0*x*x - 2.0*z*z;
  rot[6] = 2.0*y*z + 2.0*w*x;
  rot[7] = 0.0;

  rot[8] = 2.0*x*z + 2.0*w*y;
  rot[9] = 2.0*y*z - 2.0*w*x;
  rot[10] = 1.0 - 2.0*x*x - 2.0*y*y;
  rot[11] = 0.0;

  rot[12] = 0.0;
  rot[13] = 0.0;
  rot[14] = 0.0;
  rot[15] = 1.0;
}

/* ------------------ mult_quat ------------------------ */

void mult_quat(float x[4], float y[4], float z[4]){
  float z2[4];

  z2[0] = x[0]*y[0] - x[1]*y[1] - x[2]*y[2] - x[3]*y[3];
  z2[1] = x[0]*y[1] + x[1]*y[0] + x[2]*y[3] - x[3]*y[2];
  z2[2] = x[0]*y[2] - x[1]*y[3] + x[2]*y[0] + x[3]*y[1];
  z2[3] = x[0]*y[3] + x[1]*y[2] - x[2]*y[1] + x[3]*y[0];
  z[0]=z2[0];
  z[1]=z2[1];
  z[2]=z2[2];
  z[3]=z2[3];
}

/* ------------------ normalize_quat ------------------------ */

void normalize_quat(float x[4]){
  float sum;

  sum = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
  if(sum>0.0){
    x[0]/=sum;
    x[1]/=sum;
    x[2]/=sum;
    x[3]/=sum;
  }
}

/* ------------------ xyz2azelev ------------------------ */

void xyz2azelev(float *xyz,float *azimuth, float *elevation){
  float norm3;

  // x = ||xyz||cos(az)*cos(elev)
  // y = ||xyz||sin(az)*cos(elev)
  // z = ||xyz||sin(elev)
  // elev=asin(z/||xyz||)
  // az=atan(y/x)
  norm3=NORM3(xyz);
  if(norm3>0.00001&&ABS(xyz[2]/norm3)<=1.0){
    *elevation=RAD2DEG*asin(xyz[2]/norm3);
  }
  else{
    *elevation=0.0;
  }
  *azimuth=RAD2DEG*atan2(xyz[1],xyz[0]);
}
