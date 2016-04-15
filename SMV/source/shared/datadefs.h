#ifndef DATADEFS_H_DEFINED
#define DATADEFS_H_DEFINED

#ifndef ADDPROCINFO
#define ADDPROCINFO(procinfo,nprocinfo,proc,proc_id) \
  procinfo[nprocinfo].rollout = proc; \
  procinfo[nprocinfo].rollout_id = proc_id; \
  nprocinfo++
#endif

#define ONEORZERO(val) if(val!=0)val=1

#define K2C(T) ((T)-273.15)
#define C2K(T) ((T)+273.15)

#define SCALE2FDS(x) ((x)*xyzmaxdiff)
#define SCALE2SMV(x) ((x)/xyzmaxdiff)

#define SCALE2FDSL(x) ((x)*xyzmaxdiff_local)

#define YES 1
#define NO 0

#define NORMALIZE_X(x) (((x)-xbar0)/xyzmaxdiff)
#define NORMALIZE_Y(y) (((y)-ybar0)/xyzmaxdiff)
#define NORMALIZE_Z(z) (((z)-zbar0)/xyzmaxdiff)
#define NORMALIZE_ZZ(z) (((z)-terrain_zmin)/(terrain_zmax-terrain_zmin))

#define DENORMALIZE_X(x) (xbar0+(x)*xyzmaxdiff)
#define DENORMALIZE_Y(y) (ybar0+(y)*xyzmaxdiff)
#define DENORMALIZE_Z(z) (zbar0+(z)*xyzmaxdiff)

#define DENORMALIZE_XX(x) (xbar0+(x)*(xbarORIG-xbar0))
#define DENORMALIZE_YY(y) (ybar0+(y)*(ybarORIG-ybar0))
#define DENORMALIZE_ZZ(z) (zbar0+(z)*(zbarORIG-zbar0))

#define DENORMALIZE_XYZ(XYZ_OUT,XYZ_IN)\
(XYZ_OUT)[0] = DENORMALIZE_X((XYZ_IN)[0]);\
(XYZ_OUT)[1] = DENORMALIZE_Y((XYZ_IN)[1]);\
(XYZ_OUT)[2] = DENORMALIZE_Z((XYZ_IN)[2])

#define NORMALIZE_XYZ(XYZ_OUT,XYZ_IN)\
(XYZ_OUT)[0] = NORMALIZE_X((XYZ_IN)[0]);\
(XYZ_OUT)[1] = NORMALIZE_Y((XYZ_IN)[1]);\
(XYZ_OUT)[2] = NORMALIZE_Z((XYZ_IN)[2])

#define INCIRCLE(x,y,z,incirc) \
{\
  float ddx, ddy, ddz;\
  ddx = DENORMALIZE_X(x)-cvi->origin[0];\
  ddy = DENORMALIZE_Y(y)-cvi->origin[1];\
  ddz = DENORMALIZE_Z(z)-cvi->origin[2];\
  incirc=( ddx*ddx + ddy*ddy + ddz*ddz <= cvi->radius*cvi->radius ? 1 : 0 );\
}

#define VEC2MA(vec,a)\
      (vec)[0] *= (a);\
      (vec)[1] *= (a)

#define VEC3MA(vec,a)\
  (vec)[0] *= (a);\
  (vec)[1] *= (a);\
  (vec)[2] *= (a)

#define VEC2DA(vec,a)\
  (vec)[0] /= (a);\
  (vec)[1] /= (a)

#define VEC3DA(vec,a)\
  (vec)[0] /= (a);\
  (vec)[1] /= (a);\
  (vec)[2] /= (a)

#ifndef DOT2
#define DOT2(x,y) ((x)[0]*(y)[0]+(x)[1]*(y)[1])
#endif

#ifndef ROTATE
#define ROTATE(xto,xfrom,az)\
  (xto)[0] = (xfrom)[0]*cos((az)) - (xfrom)[1]*sin((az));\
  (xto)[1] = (xfrom)[0]*sin((az)) + (xfrom)[1]*cos((az))
#endif

//   i    j    k
// x[0] x[1] x[2]
// y[0] y[1] y[2]

#ifndef CROSS
#define CROSS(xy,x,y) \
  (xy)[0] = (x)[1]*(y)[2] - (y)[1]*(x)[2];\
  (xy)[1] = (x)[2]*(y)[0] - (y)[2]*(x)[0];\
  (xy)[2] = (x)[0]*(y)[1] - (y)[0]*(x)[1]
#endif

#ifndef VECEQ3
#define VECEQ3(y,x)\
  (y)[0]=(x)[0];\
  (y)[1]=(x)[1];\
  (y)[2]=(x)[2]
#endif

#ifndef VECEQCONS
#define VECEQCONS(y,x)\
  (y)[0]=(x);\
  (y)[1]=(x);\
  (y)[2]=(x)
#endif

#ifndef VECADD3
#define VECADD3(ypx,x,y)\
  (ypx)[0]=(y)[0]+(x)[0];\
  (ypx)[1]=(y)[1]+(x)[1];\
  (ypx)[2]=(y)[2]+(x)[2]
#endif

#ifndef VECDIFF3
#define VECDIFF3(ymx,y,x)\
  (ymx)[0]=(y)[0]-(x)[0];\
  (ymx)[1]=(y)[1]-(x)[1];\
  (ymx)[2]=(y)[2]-(x)[2]
#endif

#ifndef DOT3
#define DOT3(x,y) ((x)[0]*(y)[0]+(x)[1]*(y)[1]+(x)[2]*(y)[2])
#endif

#ifndef DOT3SKIP
#define DOT3SKIP(x,ix,y,iy) ((x)[0]*(y)[0]+(x)[ix]*(y)[iy]+(x)[2*ix]*(y)[2*iy])
#endif

#ifndef DOT4SKIP
#define DOT4SKIP(x,ix,y,iy) ((x)[0]*(y)[0]+(x)[ix]*(y)[iy]+(x)[2*ix]*(y)[2*iy]+(x)[3*ix]*(y)[3*iy])
#endif

#ifndef MAXDIFF2
#define MAXDIFF2(x,y) MAX(ABS(x[0]-y[0]),ABS(x[1]-y[1]))
#endif

#ifndef MAXDIFF3
#define MAXDIFF3(x,y) MAX(  MAXDIFF2(x,y),ABS(x[2]-y[2])  )
#endif

#ifndef NORM3
#define NORM3(x) sqrt((x)[0]*(x)[0]+(x)[1]*(x)[1]+(x)[2]*(x)[2])
#endif

#ifndef NORMALIZE3
#define NORMALIZE3(x)\
  {\
  float denom;\
  denom=NORM3(x);\
  if(denom!=0.0){\
  (x)[0]/=denom;\
  (x)[1]/=denom;\
  (x)[2]/=denom;\
  }\
  }
#endif

#ifndef NORM2
#define NORM2(x) sqrt((x)[0]*(x)[0]+(x)[1]*(x)[1])
#endif

#ifndef DEG2RAD
#define DEG2RAD (3.14159265359f/180.0)
#endif

#ifndef RAD2DEG
#define RAD2DEG (180.0/3.14159265359f)
#endif

#ifndef PERCENT
#define PERCENT(num,denom)  ((int)(100.0*(float)(num)/(float)(denom)+0.5))
#endif

#ifndef MAX
#define MAX(a,b)  ((a)>(b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b)  ((a)<(b) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(a) ((a)>=0 ? (a) : (-(a)))
#endif

#ifndef SIGN
#define SIGN(x) ((x)>=0 ?  1  :  -1)
#endif

#ifndef MIX
#define MIX(f,a,b) ( (f)*(a) + (1.0-(f))*(b))
#endif

#ifndef MIX2
#define MIX2(i,n,a,b) ( ((float)(i)/(float)(n))*(a) + ( 1.0-((float)(i)/(float)(n)) )*(b))
#endif

#ifndef CLAMP
#define CLAMP(x,lo,hi)  MIN(MAX((x),(lo)),(hi))
#endif

#ifndef GETINDEX
#define GETINDEX(ival,xval,xmin,dx,nx) ival = ((xval)-(xmin))/(dx); ival = CLAMP(ival,0,(nx)-1)
#endif

#ifndef IJCIRC
#define IJCIRC(i,j) ((i)+1+((j)+1)*nx)
#endif

#ifndef IJKNODE
#define IJKNODE(i,j,k) ((i)+(j)*nx+(k)*nxy)
#endif

#ifndef IJKN
#define IJKN(i,j,k,n) (IJKNODE(i,j,k)+(n)*nxyz)
#endif

#ifndef GET_QVAL
#define GET_QVAL(i,j,k,n) \
  if(cache_qdata==1){\
    qval=qdata[IJKN(i,j,k,n)];\
  }\
  else{\
    float *qvals;\
    qvals=p3levels256[n];\
    qval=qvals[iqdata[IJKN(i,j,k,n)]];\
  }
#endif

#ifndef IJKCELL
#define IJKCELL(i,j,k) ((i)+ (j)*ibar+(k)*ibar*jbar)
#endif

#ifndef IJKCELL2
#define IJCELL2(i,j) (nxcell*(j) + (i))
#endif

#ifndef IJ
#define IJ(i,j) ((i)+(j)*nx)
#endif

#ifndef IJ2
#define IJ2(i,j) ((nycell+1)*(i) + (j))
#endif

#endif
