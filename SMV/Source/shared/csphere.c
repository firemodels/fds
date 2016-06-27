#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef pp_DRAWISO
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#endif
#include "MALLOC.h"
#include <math.h>
#include "csphere.h"

/* ------------------ initspherepoints ------------------------ */

void initspherepoints(spherepoints *sphereinfo, int n){
  int i,j;
  float pi;
  float *normals;
  float phik, thetaj;
  float dphi, dtheta;
  float x, y, z;
  float cosphik;
  float maxerr=0.0;
  float *xyzfrom, *xyzto;
  float mindist;
  short *snormals;

  pi=4.0*atan(1.0);
  sphereinfo->n=n;
  sphereinfo->pi=pi;
  dphi=pi/n;
  sphereinfo->dphi=dphi;

  sphereinfo->dtheta=NULL;
  sphereinfo->nlong=NULL;
  sphereinfo->vallist=NULL;
  sphereinfo->normals=NULL;
  sphereinfo->snormals=NULL;

  NewMemory((void **)&sphereinfo->dtheta,(n+1)*sizeof(float));
  NewMemory((void **)&sphereinfo->nlong,(n+1)*sizeof(int));

  // allocate and define vallist

  NewMemory((void **)&sphereinfo->vallist,(n+1)*sizeof(int));
  sphereinfo->vallist[0]=1;
  sphereinfo->nlong[0]=1;
  for(i=1;i<n;i++){
    sphereinfo->nlong[i]=(int)(2.0*n*sin(i*sphereinfo->dphi)+0.5);
    sphereinfo->dtheta[i]=2.0*pi/sphereinfo->nlong[i];
    sphereinfo->vallist[i]=sphereinfo->vallist[i-1]+sphereinfo->nlong[i];
  }
  sphereinfo->nlong[n]=1;
  sphereinfo->vallist[n]=sphereinfo->vallist[n-1]+1;
  sphereinfo->npoints=sphereinfo->vallist[n];

  // allocate and define normal array

  NewMemory((void **)&sphereinfo->normals,3*(sphereinfo->npoints)*sizeof(float));
  NewMemory((void **)&sphereinfo->snormals,3*(sphereinfo->npoints)*sizeof(short));
  normals=sphereinfo->normals;
  snormals=sphereinfo->snormals;
  *normals++=0.0;
  *normals++=0.0;
  *normals++=1.0;
  *snormals++=0;
  *snormals++=0;
  *snormals++=32767;
  for(i=1;i<n;i++){
    phik = pi/2.0-i*dphi;
    dtheta=2.0*pi/sphereinfo->nlong[i];
    z=sin(phik);
    cosphik = cos(phik);
    for(j=0;j<sphereinfo->nlong[i];j++){
      thetaj=j*dtheta;
      x=cos(thetaj)*cosphik;
      y=sin(thetaj)*cosphik;
      *normals++=x;
      *normals++=y;
      *normals++=z;
      *snormals++=(short)(x*32767);
      *snormals++=(short)(y*32767);
      *snormals++=(short)(z*32767);
    }
  }
  *normals++=0.0;
  *normals++=0.0;
  *normals++=-1.0;
  *snormals++=0;
  *snormals++=0;
  *snormals++=-32767;

  mindist = 100000000.0;
  xyzfrom = sphereinfo->normals;
  for(i=0;i<n+1;i++){
    int ibeg, iend;
    int ii, jj;
    float xfrom, yfrom, zfrom;
    float xto, yto, zto;
    float dist;

    ibeg=i-1;
    if(ibeg<0)ibeg=0;
    iend=i+1;
    if(iend>n)iend=n;
    for(j=0;j<sphereinfo->nlong[i];j++){
      xfrom = *xyzfrom++;
      yfrom = *xyzfrom++;
      zfrom = *xyzfrom++;

      xyzto = sphereinfo->normals + 3*(sphereinfo->vallist[ibeg]-sphereinfo->nlong[ibeg]);
      for(ii=ibeg;ii<=iend;ii++){
        for(jj=0;jj<sphereinfo->nlong[ii];jj++){
          if(xyzto==xyzfrom-3){
            xyzto+=3;
            continue;
          }
          xto = *xyzto++;
          yto = *xyzto++;
          zto = *xyzto++;

          dist  = (xfrom-xto)*(xfrom-xto);
          dist += (yfrom-yto)*(yfrom-yto);
          dist += (zfrom-zto)*(zfrom-zto);
          dist = sqrt(dist);
          if(dist<mindist)mindist=dist;
        }
      }
      if(mindist>maxerr)maxerr=mindist;
    }
  }
  sphereinfo->maxerr_deg=asin(maxerr/2.0)*180.0/pi;
}

/* ------------------ freespherepoints ------------------------ */

void freespherepoints(spherepoints *sphereinfo){
  FREEMEMORY(sphereinfo->dtheta);
  FREEMEMORY(sphereinfo->nlong);
  FREEMEMORY(sphereinfo->vallist);
  FREEMEMORY(sphereinfo->normals);
}


/* ------------------ getnormalvectorptr ------------------------ */

float *getnormalvectorptr(spherepoints *sphereinfo, unsigned int index){
  float *normptr;

  if(index>sphereinfo->npoints)index=sphereinfo->npoints;
  normptr=sphereinfo->normals+3*index;
  return normptr;
}

/* ------------------ getnormalvector ------------------------ */

void getnormalvector(spherepoints *sphereinfo, unsigned int index, float *normal){
  float *normptr;

  if(index>sphereinfo->npoints)index=sphereinfo->npoints;
  normptr=sphereinfo->normals+3*index;
  normal[0]=normptr[0];
  normal[1]=normptr[1];
  normal[2]=normptr[2];
}

/* ------------------ getnormalindex2 ------------------------ */

unsigned int getnormalindex2(spherepoints *sphereinfo, float *normal){
  float norm;
  float x, y, z;
  float theta, phi;
  float pi;
  int j, k;
  int n;
  unsigned int returnval;

  pi=sphereinfo->pi;
  n=sphereinfo->n;
  x = normal[0];
  y = normal[1];
  z = normal[2];
  norm = sqrt(x*x + y*y + z*z);
  if(norm==0.0)norm=1.0;
  x/=norm;
  y/=norm;
  z/=norm;
  theta=atan2(y,x);
  if(theta<0.0)theta+=2.0*pi;
  phi=asin(z);

  k=(pi/2-phi)/sphereinfo->dphi+0.5;
  if(k<0)k=0;
  if(k>n)k=n;
  if(k==0){
    returnval=0;
  }
  else if(k==n){
    returnval=sphereinfo->npoints;
  }
  else{
    j=theta/sphereinfo->dtheta[k]+0.5;
    returnval=sphereinfo->vallist[k-1]+j;
    if(returnval>sphereinfo->npoints)returnval=sphereinfo->npoints;
  }
  return returnval;
}

/* ------------------ getnormalindex ------------------------ */

unsigned int getnormalindex(spherepoints *sphereinfo, float *normal){
  float norm;
  float x, y, z;
  unsigned int i;
  unsigned int returnval;
  float mindist2;
  float dist2;
  float *xyznorm;

  x = normal[0];
  y = normal[1];
  z = normal[2];
  norm = sqrt(x*x + y*y + z*z);
  if(norm==0.0)norm=1.0;
  x/=norm;
  y/=norm;
  z/=norm;

  returnval=0;
  mindist2 = 100000000.0;
  for(i=0;i<sphereinfo->npoints;i++){
    float xx, yy, zz;

    xyznorm = sphereinfo->normals + 3*i;
    xx = x - xyznorm[0];
    yy = y - xyznorm[1];
    zz = z - xyznorm[2];
    dist2 =  xx*xx + yy*yy + zz*zz;
    if(dist2<mindist2){
      mindist2 = dist2;
      returnval=i;
    }
  }
  return returnval;
}
#ifdef pp_DRAWISO

/* ------------------ drawspherepoints ------------------------ */

void drawspherepoints(spherepoints *spherei){
  float u[3], v[3], *w;
#define NPOINTS 64
  float cosang[NPOINTS], sinang[NPOINTS];
  float pi;
  int i,j,k;
  int index=-1;
  float rad;

  rad=spherei->dphi/2.0;

  pi=4.0*atan(1.0);
  for(i=0;i<NPOINTS;i++){
    float angle;

    angle = i*2.0*pi/NPOINTS;
    cosang[i] = cos(angle);
    sinang[i] = sin(angle);

  }

  for(i=0;i<=spherei->n;i++){
    for(j=0;j<spherei->nlong[i];j++){
      index++;
      w=spherei->normals+3*index;

      if(w[0]==0.0&&w[1]==0.0&&w[2]>0.0){
        u[0]=1.0;u[1]=0.0;u[2]=0.0;
        v[0]=0.0;v[1]=1.0;v[2]=0.0;
      }
      else if(w[0]==0.0&&w[1]==0.0&&w[2]<0.0){
        u[0]=0.0;u[1]=1.0;u[2]=0.0;
        v[0]=1.0;v[1]=0.0;v[2]=0.0;
      }
      else{
        float norm;

        u[0]=-w[1];
        u[1]=w[0];
        u[2]=0.0;
        norm = u[0]*u[0]+u[1]*u[1];
        norm=sqrt(norm);
        u[0]/=norm;
        u[1]/=norm;

        v[0]=-w[0]*w[2];
        v[1]=-w[1]*w[2];
        v[2]=w[0]*w[0]+w[1]*w[1];
        norm = v[0]*v[0]+v[1]*v[1]+v[2]*v[2];;
        norm=sqrt(norm);
        v[0]/=norm;
        v[1]/=norm;
        v[2]/=norm;
      }
      glBegin(GL_POLYGON);
      glColor3f(0.0,0.0,1.0);
      for(k=0;k<NPOINTS;k++){
        float xp, yp, zp;

        xp = 0.25 + (w[0] + rad*(cosang[k]*u[0] + sinang[k]*v[0]))/2.0;
        yp = 0.25 + (w[1] + rad*(cosang[k]*u[1] + sinang[k]*v[1]))/2.0;
        zp = 0.25 + (w[2] + rad*(cosang[k]*u[2] + sinang[k]*v[2]))/2.0;
        glVertex3f(xp,yp,zp);
      }
      glEnd();
    }
  }
}
#endif
