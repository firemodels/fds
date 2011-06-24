// $Date$ 
// $Revision$
// $Author$

#ifndef MAX
#define MAX(a,b)  ((a)>(b) ? (a) : (b))
#define MIN(a,b)  ((a)<(b) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(a) ((a)>=0 ? (a) : (-(a)))
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

#ifndef IJKCELL
#define IJKCELL(i,j,k) ((i)+ (j)*ibar+(k)*ibar*jbar)
#endif

