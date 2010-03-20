// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdlib.h>
#ifdef _DEBUG
#include <stdio.h>
#endif
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "contourdefs.h"
#include "MALLOC.h"

// svn revision character string
char scontour2d_revision[]="$Revision$";

/*  
  contouring algorithm
  --------------------

  task:  produce shaded color contours for values v_{i,j} at
         points (x_i,y_j): 0<=i<=ibar and 0<=j<=jbar.

  sub-task: shade the portion of each rectangle (formed by points 
           (x_i,y_j) (x_i+1,y_j+1) ) that lie between the levels
           contourlow and contourhigh

  there are 81 cases to consider since each of the four corners
  of the above rectangle can be 
      1) below contourlow, 
      2) between contourlow and contourhigh
      3) above contourhigh

  The 2-d array contourfill_list defined below enumerates each of the
  81 cases.  The 2-d array contourline_list describes the 81 cases for
  drawing contour lines (rather than shaded contours).



*/
int contourfill_list[81][9]={
  {0},{3,-5,3,-7},{4,-5,-6,-8,-7},
  {3,-3,2,-5},{4,-3,2,3,-7},{5,-3,2,-6,-8,-7},
  {4,-3,-4,-6,-5},{5,-3,-4,-6,3,-7},{4,-3,-4,-8,-7},
  {3,-1,1,-3},{6,-1,1,-3,-5,3,-7},{7,-1,1,-3,-5,-6,-8,-7},
  {4,-1,1,2,-5},{5,-1,1,2,3,-7},{6,-1,1,2,-6,-8,-7},
  {5,-1,1,-4,-6,-5},{6,-1,1,-4,-6,3,-7},{5,-1,1,-4,-8,-7},

  {4,-1,-2,-4,-3},{7,-1,-2,-4,-3,-5,3,-7},{8,-1,-2,-4,-3,-5,-6,-8,-7},
  {5,-1,-2,-4,2,-5},{6,-1,-2,-4,2,3,-7},{7,-1,-2,-4,2,-6,-8,-7},
  {4,-1,-2,-6,-5},{5,-1,-2,-6,3,-7},{4,-1,-2,-8,-7},
  {3,0,-1,-7},{4,0,-1,-5,3},{5,0,-1,-5,-6,-8},
  {6,0,-1,-3,2,-5,-7},{5,0,-1,-3,2,3},{6,0,-1,-3,2,-6,-8},
  {7,0,-1,-3,-4,-6,-5,-7},{6,0,-1,-3,-4,-6,3},{5,0,-1,-3,-4,-8},

  {4,0,1,-3,-7},{5,0,1,-3,-5,3},{6,0,1,-3,-5,-6,-8},
  {5,0,1,2,-5,-7},{4,0,1,2,3},{5,0,1,2,-6,-8},
  {6,0,1,-4,-6,-5,-7},{5,0,1,-4,-6,3},{4,0,1,-4,-8},
  {5,0,-2,-4,-3,-7},   {6,0,-2,-4,-3,-5,3},{7,0,-2,-4,-3,-5,-6,-8},
  {6,0,-2,-4,2,-5,-7},{5,0,-2,-4,2,3},{6,0,-2,-4,2,-6,-8},
  {5,0,-2,-6,-5,-7},{4,0,-2,-6,3},{3,0,-2,-8},

  {4,-2,-1,-7,-8},{5,-2,-1,-5,3,-8},{4,-2,-1,-5,-6},
  {7,-2,-1,-3,2,-5,-7,-8},{6,-2,-1,-3,2,3,-8},{5,-2,-1,-3,2,-6},
  {8,-2,-1,-3,-4,-6,-5,-7,-8},{7,-2,-1,-3,-4,-6,3,-8},{4,-2,-1,-3,-4},
  {5,-2,1,-3,-7,-8},{6,-2,1,-3,-5,3,-8},{5,-2,1,-3,-5,-6},
  {6,-2,1,2,-5,-7,-8},{5,-2,1,2,3,-8},{4,-2,1,2,-6},
  {7,-2,1,-4,-6,-5,-7,-8},{6,-2,1,-4,-6,3,-8},{3,-2,1,-4},

  {4,-4,-3,-7,-8},{5,-4,-3,-5,3,-8},{4,-4,-3,-5,-6},
  {5,-4,2,-5,-7,-8},{4,-4,2,3,-8},{3,-4,2,-6},
  {4,-6,-5,-7,-8},{3,-6,3,-8},{0}
};

int contourline_list[81][5]={
  {0},{2,5,7},{2,5,7},
  {2,3,5},{2,3,7},{2,3,7},
  {2,3,5},{2,3,7},{2,4,8},
  {2,1,3},{4,1,7,3,5},{4,3,5,7,1},
  {2,5,1},{2,1,7},{2,1,7},
  {2,1,5},{2,1,7},{2,1,7},

  {2,1,3},{4,1,7,3,5},{4,1,7,3,5},
  {2,1,5},{2,1,7},{2,1,7},
  {2,1,5},{2,1,7},{2,1,7},
  {2,1,7},{2,1,5},{2,1,5},
  {4,1,3,5,7},{2,1,3},{2,1,3},
  {4,1,3,5,7},{2,1,3},{2,1,3},

  {2,3,7},{2,3,5},{2,3,5},
  {2,5,7},{0},{0},
  {2,5,7},{0},{0},
  {2,3,7},{2,3,5},{2,3,5},
  {2,5,7},{0},{0},
  {2,5,7},{0},{0},

  {2,1,7},{2,1,5},{2,1,5},
  {4,1,3,5,7},{2,1,3},{2,1,3},
  {4,1,3,5,7},{2,1,3},{2,1,3},
  {2,3,7},{2,3,5},{2,3,5},
  {2,5,7},{0},{0},
  {2,5,7},{0},{0},

  {2,3,7},{2,3,5},{2,3,5},
  {2,5,7},{0},{0},
  {2,5,7},{0},{0}
};

/*  ------------------ initcontour ------------------------ */

void initcontour(contour *ci, float **rgbptr, int nlevels){

  int n;

  ci->nlevels=nlevels;
  ci->rgbptr=rgbptr;
  NewMemory((void **)&ci->levels,nlevels*sizeof(float));
  NewMemory((void **)&ci->nnodes,nlevels*sizeof(int));
  NewMemory((void **)&ci->npolys,nlevels*sizeof(int));
  NewMemory((void **)&ci->nlines,nlevels*sizeof(int));

  NewMemory((void **)&ci->polysize,nlevels*sizeof(int *));
  NewMemory((void **)&ci->xnode,nlevels*sizeof(float *));
  NewMemory((void **)&ci->ynode,nlevels*sizeof(float *));
  NewMemory((void **)&ci->xlines,nlevels*sizeof(float *));
  NewMemory((void **)&ci->ylines,nlevels*sizeof(float *));
  for(n=0;n<nlevels;n++){
    ci->polysize[n]=NULL;
    ci->xnode[n]=NULL;
    ci->ynode[n]=NULL;
    ci->xlines[n]=NULL;
    ci->ylines[n]=NULL;
  }
}

/*  ------------------ freecontour ------------------------ */

void freecontour(contour *ci){

  int n;

  FREEMEMORY(ci->levels);
  FREEMEMORY(ci->nnodes);
  FREEMEMORY(ci->npolys);
  FREEMEMORY(ci->nlines);

  for(n=0;n<ci->nlevels;n++){
    FREEMEMORY(ci->polysize[n]);
    FREEMEMORY(ci->xnode[n]);
    FREEMEMORY(ci->ynode[n]);
    FREEMEMORY(ci->xlines[n]);
    FREEMEMORY(ci->ylines[n]);
  }
  FREEMEMORY(ci->polysize);
  FREEMEMORY(ci->xnode);
  FREEMEMORY(ci->ynode);
  FREEMEMORY(ci->xlines);
  FREEMEMORY(ci->ylines);
  ci->nlevels=0;
}

/*  ------------------ setcontourslice ------------------------ */

void setcontourslice(contour *ci,int idir,float xyz){
  ci->idir=idir;
  ci->xyzval=xyz;

}

/*  ------------------ getcontours ------------------------ */

void getcontours(const  float *xgrid, const float *ygrid, int nx, int ny,  
                 const float *vals, const char *iblank, const float *levels,
                 const contour *ci){
  int n,i,j;
  double x[4],y[4],v[4];
  double contlow, conthigh;
  float *xnode=NULL, *ynode=NULL;
  float *xnodecopy, *ynodecopy;
  float *xline=NULL, *yline=NULL;
  float *xlinecopy, *ylinecopy;
  int *polysize;
  int nnodes, npolys;
  int mxnodes, mxpolys,mxlines;
  int nnode,nnode2,casen,doit;
  int ij0,ij2,i2j,i2j2;
  int lastcasenum;
  int nlinepts;
  int nlevels;
  int blankit=0;
  int minfill, maxfill;

#define ij(i,j) ((i)*ny+(j))
#define ijblank(i,j) ((i)*(ny-1)+(j))


  nlevels=ci->nlevels;
  for(n=0;n<nlevels-2;n++){
    ci->levels[n]=levels[n];
  }
  for(n=0;n<nlevels;n++){
    minfill=0;
    maxfill=0;
    if(n==nlevels-2){
      contlow=(double)levels[0]-(levels[1]-levels[0]);
      conthigh=(double)levels[0];
      minfill=1;

    }
    else if(n==nlevels-1){
      contlow=levels[nlevels-2];
      conthigh=levels[nlevels-2]+(levels[nlevels-3]-levels[nlevels-4]);
      maxfill=1;
    }
    else{
      contlow=(double)levels[n];
      conthigh=(double)levels[n+1];
    }

    mxpolys=(nx+1)*(ny+1)+1;
    mxnodes=8*mxpolys;
    mxlines=4*mxpolys;
    NewMemory((void **)&xnode,mxnodes*sizeof(float));
    xnodecopy=xnode;
    NewMemory((void **)&xline,mxlines*sizeof(float));
    xlinecopy=xline;
    NewMemory((void **)&ynode,mxnodes*sizeof(float));
    ynodecopy=ynode;
    NewMemory((void **)&yline,mxlines*sizeof(float));
    ylinecopy=yline;

    NewMemory((void **)&polysize,mxpolys*sizeof(int));

    npolys = 0;
    nnodes=0;
    nlinepts=0;
    for(j=0;j<ny-1;j++){
      y[0]=(double)ygrid[j];
      y[1]=(double)ygrid[j+1];
      y[2]=(double)ygrid[j+1];
      y[3]=(double)ygrid[j];
      lastcasenum=-1;
      for(i=0;i<nx-1;i++){
        x[0]=(double)xgrid[i];
        x[1]=(double)xgrid[i];
        x[2]=(double)xgrid[i+1];
        x[3]=(double)xgrid[i+1];
        doit=1;
        ij0=ij(i,j);
        ij2=ij(i,j+1);
        i2j2=ij(i+1,j+1);
        i2j=ij(i+1,j);
        if(iblank!=NULL&&iblank[ijblank(i,j)]!=2)doit=0;
        if(doit==1){
          v[0]=(double)vals[ij0];
          v[1]=(double)vals[ij2];
          v[2]=(double)vals[i2j2];
          v[3]=(double)vals[i2j];
          getcontournodes(n,nlevels,x,y,v,contlow,minfill,conthigh,maxfill,
            &nnode, xnodecopy, ynodecopy,
            &nnode2,xlinecopy,ylinecopy,
            &casen,blankit);
          if(nnode!=0){
          if(casen!=40||lastcasenum!=40){
              nnodes += nnode;
              xnodecopy += nnode;
              ynodecopy += nnode;
              polysize[npolys]=nnode;
              npolys++;
           }
            else{
              xnodecopy[-2]=xnodecopy[2];
              xnodecopy[-1]=xnodecopy[3];
              ynodecopy[-2]=ynodecopy[2];
              ynodecopy[-1]=ynodecopy[3];
           }
           lastcasenum=casen;
          }
          if(nnode2!=0){
            xlinecopy += nnode2;
            ylinecopy += nnode2;
            nlinepts += nnode2;
          }
        }
        else{lastcasenum=-1;}
      }
    }
    if(nnodes>0){
      ResizeMemory((void **)&xnode,nnodes*sizeof(float));
      ResizeMemory((void **)&ynode,nnodes*sizeof(float));
    }
    if(nlinepts>0){
      ResizeMemory((void **)&xline,nlinepts*sizeof(float));
      ResizeMemory((void **)&yline,nlinepts*sizeof(float));
    }
    if(npolys>0){
      ResizeMemory((void **)&polysize,npolys*sizeof(int));
    }
    ci->nnodes[n]=nnodes;
    ci->npolys[n]=npolys;
    ci->xnode[n]=xnode;
    ci->ynode[n]=ynode;
    ci->xlines[n]=xline;
    ci->ylines[n]=yline;
    ci->nlines[n]=nlinepts;
    ci->polysize[n]=polysize;
  }
  /* printf("polygon count=%i\n",npolystotal); */

}

/*  ------------------ getlinecontours ------------------------ */

void getlinecontours(const  float *xgrid, const float *ygrid, int nx, int ny,  
                 const float *vals, const char *iblank, const float *levels,
                 const contour *ci){
  int n,i,j;
  double x[4],y[4],v[4];
  double linelevel;
  float *xline=NULL, *yline=NULL;
  float *xlinecopy, *ylinecopy;
  int mxpolys,mxlines;
  int nnode2,casen,doit;
  int ij0,ij2,i2j,i2j2;
  int nlinepts;
  int nlevels;
  int blankit=0;
  int minfill, maxfill;

  nlevels=ci->nlevels;
  for(n=0;n<nlevels-2;n++){
    ci->levels[n]=levels[n];
  }
  for(n=0;n<nlevels;n++){
    linelevel=(double)levels[n];

    mxpolys=(nx+1)*(ny+1)+1;
    mxlines=4*mxpolys;
    NewMemory((void **)&xline,mxlines*sizeof(float));
    xlinecopy=xline;
    NewMemory((void **)&yline,mxlines*sizeof(float));
    ylinecopy=yline;

    nlinepts=0;
    for(j=0;j<ny-1;j++){
      y[0]=(double)ygrid[j];
      y[1]=(double)ygrid[j+1];
      y[2]=(double)ygrid[j+1];
      y[3]=(double)ygrid[j];
      for(i=0;i<nx-1;i++){
        x[0]=(double)xgrid[i];
        x[1]=(double)xgrid[i];
        x[2]=(double)xgrid[i+1];
        x[3]=(double)xgrid[i+1];
        doit=1;
        ij0=ij(i,j);
        ij2=ij(i,j+1);
        i2j2=ij(i+1,j+1);
        i2j=ij(i+1,j);
        if(iblank!=NULL&&iblank[ijblank(i,j)]!=2)doit=0;
        if(doit==1){
          v[0]=(double)vals[ij0];
          v[1]=(double)vals[ij2];
          v[2]=(double)vals[i2j2];
          v[3]=(double)vals[i2j];
          getlinecontournodes(linelevel,x,y,v,
            &nnode2,xlinecopy,ylinecopy,
            blankit);
          if(nnode2!=0){
            xlinecopy += nnode2;
            ylinecopy += nnode2;
            nlinepts += nnode2;
          }
        }
      }
    }
    if(nlinepts>0){
      ResizeMemory((void **)&xline,nlinepts*sizeof(float));
      ResizeMemory((void **)&yline,nlinepts*sizeof(float));
    }
    ci->xlines[n]=xline;
    ci->ylines[n]=yline;
    ci->nlines[n]=nlinepts;
  }
}

/*  ------------------ getcontournodes ------------------------ */

void getcontournodes(int ilev, int nlevels, const double x[4], const double y[4], const double val[4],
                     double  contlow, int minfill, double conthigh,int maxfill,
                     int *nnode, float *xnode, float *ynode,
                     int *nnode2,float *xline, float *yline,
                     int *casen,int blankit){
  int state[4]={1,1,1,1}, n,casenum;
  float xzero[9]={2.,2.,2.,2.,2.,2.,2.,2.,2.};
  float yzero[9]={2.,2.,2.,2.,2.,2.,2.,2.,2.};
  double xcopy[5],ycopy[5],valcopy[5],vallownet[5],valhighnet[5];
  int nn,edgenum;
  double factor;

  *nnode=0;
  *nnode2=0;
  if(conthigh<=contlow)return;
  if(blankit==0){
    for(n=0;n<4;n++){
      state[n]=1;
      if(minfill==0&&val[n]<=contlow){
        state[n]=0;
      }
      if(maxfill==0&&val[n]>=conthigh){
        state[n]=2;
      }
    }
    casenum=state[3]+3*state[2]+9*state[1]+27*state[0];
    if(maxfill==1&&casenum==80){
      casenum=40;
    }
    if(minfill==1&&casenum==0){
      casenum=40;
    }
  }
  else{
    casenum=40;
  }
  *casen=casenum;
  *nnode=contourfill_list[casenum][0];
  *nnode2=contourline_list[casenum][0];
  if(casenum==0)return;

  for(n=0;n<4;n++){
    xcopy[n]=x[n];
    ycopy[n]=y[n];
    valcopy[n]=val[n];
    vallownet[n]=val[n]-contlow;
    valhighnet[n]=val[n]-conthigh;
  }
  xcopy[4]=x[0]; 
  ycopy[4]=y[0];
  vallownet[4]=vallownet[0];
  valhighnet[4]=valhighnet[0];
  valcopy[4]=val[0];

  if(casenum!=40){
    for(n=0;n<4;n++){
      if(vallownet[n]*vallownet[n+1]<=0.0f){
        nn=2*n+1;
        factor = vallownet[n+1]/(valcopy[n+1]-valcopy[n]);
        xzero[nn] = xcopy[n]*factor + xcopy[n+1]*(1.0f-factor);
        yzero[nn] = ycopy[n]*factor + ycopy[n+1]*(1.0f-factor);
      }
      if(valhighnet[n]*valhighnet[n+1]<=0.0f){
        nn = 2*n+2;
        factor = valhighnet[n+1]/(valcopy[n+1]-valcopy[n]);
        xzero[nn] = xcopy[n]*factor + xcopy[n+1]*(1.0f-factor);
        yzero[nn] = ycopy[n]*factor + ycopy[n+1]*(1.0f-factor);
      }
    }
  }
  for(n=0;n<*nnode;n++){
    edgenum=contourfill_list[casenum][n+1];
    if(edgenum>=0){
      xnode[n]=x[edgenum]; 
      ynode[n]=y[edgenum];
    }
    else{
      xnode[n]=xzero[-edgenum];
      ynode[n]=yzero[-edgenum];
    }
  }
  for(n=0;n<*nnode2;n++){
    edgenum=contourline_list[casenum][n+1];
    xline[n]=xzero[edgenum];
    yline[n]=yzero[edgenum];
  }
}

/*  ------------------ getlinecontournodes ------------------------ */

void getlinecontournodes(double linelevel, const double x[4], const double y[4], const double val[4],
                     int *nline_nodes,float *xline, float *yline,
                     int blankit){
  int state[4]={1,1,1,1}, n,casenum;
  float xzero[9]={2.,2.,2.,2.,2.,2.,2.,2.,2.};
  float yzero[9]={2.,2.,2.,2.,2.,2.,2.,2.,2.};
  double xcopy[5],ycopy[5],valcopy[5],valnet[5];
  int nn,edgenum;
  double factor;

  *nline_nodes=0;
  if(blankit==0){
    for(n=0;n<4;n++){
      state[n]=1;
      if(val[n]<linelevel){
        state[n]=0;
      }
      if(val[n]>=linelevel){
        state[n]=2;
      }
    }
    casenum=state[3]+3*state[2]+9*state[1]+27*state[0];
  }
  else{
    casenum=40;
  }
  *nline_nodes=contourline_list[casenum][0];
  if(casenum==0)return;

  for(n=0;n<4;n++){
    xcopy[n]=x[n];
    ycopy[n]=y[n];
    valcopy[n]=val[n];
    valnet[n]=val[n]-linelevel;
  }
  xcopy[4]=x[0]; 
  ycopy[4]=y[0];
  valnet[4]=valnet[0];
  valcopy[4]=val[0];

  if(casenum!=40){
    for(n=0;n<4;n++){
      if(valnet[n]*valnet[n+1]<=0.0f){
        nn=2*n+1;
        factor = valnet[n+1]/(valcopy[n+1]-valcopy[n]);
        xzero[nn] = xcopy[n]*factor + xcopy[n+1]*(1.0f-factor);
        yzero[nn] = ycopy[n]*factor + ycopy[n+1]*(1.0f-factor);
      }
    }
  }
  for(n=0;n<*nline_nodes;n++){
    edgenum=contourline_list[casenum][n+1];
    xline[n]=xzero[edgenum];
    yline[n]=yzero[edgenum];
  }
}

/*  ------------------ drawcontours ------------------------ */

void DrawContours(const contour *ci){
  int nlevels, n, npolys, *polysize, ipoly,j,nnodes;
  float *xnode, *ynode,xyzval;
  float **rgb;

  int *npolysv;

  rgb=ci->rgbptr;
  nlevels=ci->nlevels;
  xyzval=ci->xyzval;
  for(n=0;n<nlevels;n++){
    xnode=ci->xnode[n];
    ynode=ci->ynode[n];
    polysize=ci->polysize[n];
    npolysv=ci->npolys;
    npolys=npolysv[n];
    if(ci->idir==1){
      glColor4fv(rgb[n]);
      for(ipoly=0;ipoly<npolys;ipoly++){
        glBegin(GL_POLYGON);
        for(j=0;j<polysize[ipoly];j++){
          glVertex3f(xyzval,*xnode++,*ynode++);
        }
        glEnd();
      }
    }
    else if(ci->idir==2){
      glColor4fv(rgb[n]);
      for(ipoly=0;ipoly<npolys;ipoly++){
        nnodes=polysize[ipoly];
        glBegin(GL_POLYGON);
        for(j=0;j<nnodes;j++){
          glVertex3f(*xnode++,xyzval,*ynode++);
        }
        glEnd();
      }
    }
    else if(ci->idir==3){
      glColor4fv(rgb[n]);
      for(ipoly=0;ipoly<npolys;ipoly++){
        glBegin(GL_POLYGON);
        for(j=0;j<polysize[ipoly];j++){
          glVertex3f(*xnode++,*ynode++,xyzval);
        }
        glEnd();
      }
    }
  }
}

/*  ------------------ DrawLineContours ------------------------ */

void DrawLineContours(const contour *ci, float linewidth){

  int nlevels, n, npolys, *polysize, ipoly,j,nnodes;
  float xyzval;
  float *xline, *yline;
  float **rgb;

  int nlinepts,iline;

  rgb=ci->rgbptr;
  nlevels=ci->nlevels;
  xyzval=ci->xyzval;
  for(n=0;n<nlevels-1;n++){
    xline=ci->xlines[n];
    yline=ci->ylines[n];
    nlinepts=ci->nlines[n];
    if(ci->idir==1){
      glColor4fv(rgb[n]);
      glLineWidth(linewidth);
      glBegin(GL_LINES);
      for(iline=0;iline<nlinepts;iline++){
        glVertex3f(xyzval,*xline++,*yline++);
      }
      glEnd();
    }
    else if(ci->idir==2){
      glColor4fv(rgb[n]);
      glLineWidth(linewidth);
      glBegin(GL_LINES);
      for(iline=0;iline<nlinepts;iline++){
        glVertex3f(*xline++,xyzval,*yline++);
      }
      glEnd();
    }
    else if(ci->idir==3){
      glColor4fv(rgb[n]);
      glLineWidth(linewidth);
      glBegin(GL_LINES);
      for(iline=0;iline<nlinepts;iline++){
        glVertex3f(*xline++,*yline++,xyzval);
      }
      glEnd();
    }
  }
}

