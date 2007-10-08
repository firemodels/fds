// $Date: 2007-10-07 22:08:47 -0400 (Sun, 07 Oct 2007) $ 
// $Revision: 800 $
// $Author: gforney $

#include "options.h"
#include <stdio.h>  
#include <math.h>

// svn revision character string
char segtri_revision[]="$Revision: 800 $";

#define HIT 1
#define MISS 0
#define USE_REJECT_TESTS

int is_seg_in_tri(const float p1[3], const float p2[3], 
                  const float t1[3], const float t2[3], const float t3[3],
                  const float fnorm[3], const int *inorm);
int seg_in_rect(float *p1, float *p2,  
                float xmin, float xmax, 
                float ymin, float ymax, 
                float zmin, float zmax,
                int checkbounds);

/* ------------------ seg_in_box ------------------------ */

int seg_in_box(float *p1, float *p2,  
               float xmin, float xmax, 
               float ymin, float ymax, 
               float zmin, float zmax){

  int checkbounds=0;

  if(p1[0]>xmax&&p2[0]>xmax)return MISS;
  if(p1[0]<xmin&&p2[0]<xmin)return MISS;

  if(p1[1]>ymax&&p2[1]>ymax)return MISS;
  if(p1[1]<ymin&&p2[1]<ymin)return MISS;

  if(p1[2]>zmax&&p2[2]>zmax)return MISS;
  if(p1[2]<zmin&&p2[2]<zmin)return MISS;

  if(seg_in_rect(p1,p2,xmin,xmin,ymin,ymax,zmin,zmax,checkbounds)==HIT)return HIT;
  if(seg_in_rect(p1,p2,xmin,xmax,ymin,ymin,zmin,zmax,checkbounds)==HIT)return HIT;
  if(seg_in_rect(p1,p2,xmin,xmax,ymin,ymax,zmin,zmin,checkbounds)==HIT)return HIT;

  if(seg_in_rect(p1,p2,xmax,xmax,ymin,ymax,zmin,zmax,checkbounds)==HIT)return HIT;
  if(seg_in_rect(p1,p2,xmin,xmax,ymax,ymax,zmin,zmax,checkbounds)==HIT)return HIT;
  if(seg_in_rect(p1,p2,xmin,xmax,ymin,ymax,zmax,zmax,checkbounds)==HIT)return HIT;

  return MISS;
}

/* ------------------ seg_in_rect ------------------------ */

int seg_in_rect(float *p1, float *p2,  
                float xmin, float xmax, 
                float ymin, float ymax, 
                float zmin, float zmax,
                int checkbounds){
  float r11[3], r22[3], r12[3], r21[3];

  if(checkbounds==1){
    if(p1[0]>xmax&&p2[0]>xmax)return MISS;
    if(p1[0]<xmin&&p2[0]<xmin)return MISS;
    if(p1[1]>ymax&&p2[1]>ymax)return MISS;
    if(p1[1]<ymin&&p2[1]<ymin)return MISS;
    if(p1[2]>zmax&&p2[2]>zmax)return MISS;
    if(p1[2]<zmin&&p2[2]<zmin)return MISS;
  }

  r11[0]=xmin;
  r11[1]=ymin;
  r11[2]=zmin;
  r22[0]=xmax;
  r22[1]=ymax;
  r22[2]=zmax;
  
  if(xmin==xmax){

    r12[0]=xmin;
    r12[1]=ymin;
    r12[2]=zmax;


    r21[0]=xmin;
    r21[1]=ymax;
    r21[2]=zmin;
  }
  else if(ymin==ymax){
    r12[0]=xmin;
    r12[1]=ymin;
    r12[2]=zmax;


    r21[0]=xmax;
    r21[1]=ymin;
    r21[2]=zmin;
  }
  else if(zmin==zmax){
    r12[0]=xmin;
    r12[1]=ymax;
    r12[2]=zmin;


    r21[0]=xmax;
    r21[1]=ymin;
    r21[2]=zmin;
  }
  else{
    return MISS;
  }
  if(is_seg_in_tri(p1,p2,r11,r12,r22,NULL,NULL)==HIT)return HIT;
  if(is_seg_in_tri(p1,p2,r11,r21,r22,NULL,NULL)==HIT)return HIT;
  return MISS;
}

/* ------------------ checksolve ------------------------ */

void checksolve(void){
  int i;
  float aa[9],a[9],b[3],x[3];
  float *a0, *a1, *a2;

  a0 = aa;
  a1 = aa+3;
  a2 = aa+6;
  for(i=0;i<9;i++){
    a[i]=1.0;
  }
  a[0]=2.0;
  a[4]=2.0;
  a[8]=2.0;
  for(i=0;i<9;i++){
    a[i]=1.1+sin(0.31+i*3.14/7.0);
    aa[i]=a[i];
  }
  b[0]=1.0;
  b[1]=2.0;
  b[2]=3.0;
  x[0] = a0[0]*b[0]+a0[1]*b[1]+a0[2]*b[2];
  x[1] = a1[0]*b[0]+a1[1]*b[1]+a1[2]*b[2];
  x[2] = a2[0]*b[0]+a2[1]*b[1]+a2[2]*b[2];
  printf("before %f %f %f\n after %f %f %f\n",b[0],b[1],b[2],x[0],x[1],x[2]);
}

/* ------------------ sge3x3slv ------------------------ */

int is_seg_in_tri(const float p1[3], const float p2[3], 
                  const float t1[3], const float t2[3], const float t3[3],
                  const float fnorm[3], const int *inorm){
  float *a0, *a1, *a2, *adummy;
  float dummy;
  float s1, s2;
  float factor;
//  float scale[3];
  float a[9],b[3];
  float sum;

  a0=a;
  a1=a+3;
  a2=a+6;

  /*
  scale[0]=fabs(a0[0]);
  if(fabs(a0[1])>scale[0])scale[0]=fabs(a0[1]);
  if(fabs(a0[2])>scale[0])scale[0]=fabs(a0[2]);
  if(scale[0]=0.0)return MISS;
  scale[1]=fabs(a1[0]);
  if(fabs(a1[1])>scale[1])scale[1]=fabs(a1[1]);
  if(fabs(a1[2])>scale[1])scale[1]=fabs(a1[2]);
  if(scale[1]=0.0)return MISS;
  scale[2]=fabs(a2[2]);
  if(fabs(a2[1])>scale[2])scale[2]=fabs(a2[1]);
  if(fabs(a2[2])>scale[2])scale[2]=fabs(a2[2]);
  if(scale[2]=0.0)return MISS;
  */

  a0[0]=t1[0]-t3[0];
  a1[0]=t1[1]-t3[1];
  a2[0]=t1[2]-t3[2];

  a0[1]=t2[0]-t3[0];
  a1[1]=t2[1]-t3[1];
  a2[1]=t2[2]-t3[2];
  a0[2]=p1[0]-p2[0];
  a1[2]=p1[1]-p2[1];
  a2[2]=p1[2]-p2[2];

  b[0]=p1[0]-t3[0];
  b[1]=p1[1]-t3[1];
  b[2]=p1[2]-t3[2];

#ifdef USE_REJECT_TESTS
  if(fnorm!=NULL){
    float b2[3];

    b2[0]=p2[0]-t3[0];
    b2[1]=p2[1]-t3[1];
    b2[2]=p2[2]-t3[2];

    s1 = fnorm[0]*b[0]+fnorm[1]*b[1]+fnorm[2]*b[2];
    s2 = fnorm[0]*b2[0]+fnorm[1]*b2[1]+fnorm[2]*b2[2];
    if(s1>0.0&&s2>0.0||s1<0.0&&s2<0.0)return MISS;
  }
  else if(inorm!=NULL){
    if(*inorm>=0&&*inorm<=2){
      s1=b[*inorm];
      s2=p2[*inorm]-t3[*inorm];
      if(s1>0.0&&s2>0.0||s1<0.0&&s2<0.0)return MISS;
    }
  }
  if(b[0]>=0.0){
    sum=0.0;
    if(a0[0]>0.0)sum+=a0[0];
    if(a0[1]>0.0)sum+=a0[1];
    if(a0[2]>0.0)sum+=a0[2];
    if(sum<b[0])return MISS;
  }
  else{
    sum=0.0;
    if(a0[0]<0.0)sum+=a0[0];
    if(a0[1]<0.0)sum+=a0[1];
    if(a0[2]<0.0)sum+=a0[2];
    if(sum>b[0])return MISS;
  }
  if(b[1]>=0.0){
    sum=0.0;
    if(a1[0]>0.0)sum+=a1[0];
    if(a1[1]>0.0)sum+=a1[1];
    if(a1[2]>0.0)sum+=a1[2];
    if(sum<b[1])return MISS;
  }
  else{
    sum=0.0;
    if(a1[0]<0.0)sum+=a1[0];
    if(a1[1]<0.0)sum+=a1[1];
    if(a1[2]<0.0)sum+=a1[2];
    if(sum>b[1])return MISS;
  }
  if(b[2]>=0.0){
    sum=0.0;
    if(a2[0]>0.0)sum+=a2[0];
    if(a2[1]>0.0)sum+=a2[1];
    if(a2[2]>0.0)sum+=a2[2];
    if(sum<b[2])return MISS;
  }
  else{
    sum=0.0;
    if(a2[0]<0.0)sum+=a2[0];
    if(a2[1]<0.0)sum+=a2[1];
    if(a2[2]<0.0)sum+=a2[2];
    if(sum>b[2])return MISS;
  }
#endif

  /* swap rows 0 and 1 (or 0 and 2) if needed to maximize pivot element 
     assume rows are reasonably scaled */
  if(fabs(a1[0])>fabs(a0[0])&&fabs(a1[0])>fabs(a2[0])){
    adummy=a1;
    a1=a0;
    a0=adummy;
    dummy=b[1];
    b[1]=b[0];
    b[0]=dummy;
  }
  else if(fabs(a2[0])>fabs(a0[0])&&fabs(a2[0])>fabs(a1[0])){
    adummy=a2;
    a2=a0;
    a0=adummy;
    dummy = b[2];
    b[2]=b[0];
    b[0]=dummy;
  }

  /* pivot element is 0 so matrix is singular */
  if(a0[0]==0.0)return MISS;
  factor=a1[0]/a0[0];
  a1[1]-=factor*a0[1];
  a1[2]-=factor*a0[2];
  b[1]-=factor*b[0];

  factor=a2[0]/a0[0];
  a2[1]-=factor*a0[1];
  a2[2]-=factor*a0[2];
  b[2]-=factor*b[0];

  /* swap rows 1 and 2 if needed to maximize pivot element 
     assume rows are reasonably scaled */
  if(fabs(a2[1])>fabs(a1[1])){
    adummy=a1;
    a1=a2;
    a2=adummy;
    dummy=b[1];
    b[1]=b[2];
    b[2]=dummy;
  }

  /* pivot element is 0 so matrix is singular */
  if(a1[1]==0.0)return MISS;
  factor=a2[1]/a1[1];
  a2[2]-=factor*a1[2];
  b[2]-=factor*b[1];

  /* pivot element is 0 so matrix is singular */
  if(a2[2]==0.0)return MISS;

#ifdef USE_REJECT_TESTS
  // 0 <= b[2]/a2[2] <= 1
  if(a2[2]>0.0&&(b[2]<0.0||b[2]>a2[2]))return MISS;
  if(a2[2]<0.0&&(b[2]>0.0||b[2]<a2[2]))return MISS;
  b[2]/=a2[2];

  b[1]-=a1[2]*b[2];
  // 0 <= b[1]/a1[1] <= 1
  if(a1[1]>0.0&&(b[1]<0.0||b[1]>a1[1]))return MISS;
  if(a1[1]<0.0&&(b[1]>0.0||b[1]<a1[1]))return MISS;
  b[1]/=a1[1];

  b[0]-=(a0[1]*b[1]+a0[2]*b[2]);
  // 0 <= b[0]/a0[0] <= 1
  if(a0[0]>0.0&&(b[0]<0.0||b[0]>a0[0]))return MISS;
  if(a0[0]<0.0&&(b[0]>0.0||b[0]<a0[0]))return MISS;
  // no need to compute b[0]/=a0[0]; since if we get here
  // it is within bounds
#else
  b[2]/=a2[2];
  b[1]-=a1[2]*b[2];
  b[1]/=a1[1];
  b[0]-=(a0[1]*b[1]+a0[2]*b[2]);
  b[0]/=a0[0];

  if(b[2]<0.0||b[2]>1.0)return MISS;
  if(b[1]<0.0||b[1]>1.0)return MISS;
  if(b[0]<0.0||b[0]>1.0)return MISS;
#endif

  return HIT;
}

