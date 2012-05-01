// $Date$ 
// $Revision$
// $Author$

#ifndef DATADEFS_H_DEFINED
#define DATADEFS_H_DEFINED

#ifndef PIFACTOR
#define PIFACTOR (3.14159/180.0)
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

#ifndef IJKNODE
#define IJKNODE(i,j,k) ((i)+(j)*nx+(k)*nxy)
#endif

#ifndef IJKN
#define IJKN(i,j,k,n) (IJKNODE(i,j,k)+(n)*nxyz)
#endif

#ifndef GET_QVAL
#define GET_QVAL(i,j,k,n) \
  if(unload_qdata==0){\
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

#ifndef MAXREV
#define MAXREV(cval) rev=getrevision(cval);max_revision=MAX(rev,max_revision)
#endif

#endif
