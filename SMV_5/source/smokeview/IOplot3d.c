// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>  
#include <stdlib.h>
#include <math.h>
#include "flowfiles.h"
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "contourdefs.h"
#include "MALLOC.h"
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

// svn revision character string
char IOplot3d_revision[]="$Revision$";

/* ------------------ plot3dcompare  ------------------------ */

int plot3dcompare( const void *arg1, const void *arg2 ){
  plot3d *plot3di, *plot3dj;

  plot3di = plot3dinfo + *(int *)arg1;
  plot3dj = plot3dinfo + *(int *)arg2;

  if(strcmp(plot3di->longlabel,plot3dj->longlabel)<0)return -1;
  if(strcmp(plot3di->longlabel,plot3dj->longlabel)>0)return 1;
  if(plot3di->time<plot3dj->time)return -1;
  if(plot3di->time>plot3dj->time)return 1;
  if(plot3di->blocknumber<plot3dj->blocknumber)return -1;
  if(plot3di->blocknumber>plot3dj->blocknumber)return 1;
  return 0;
}

#define ijknode(i,j,k) ((i)+(j)*nx+(k)*nxy)
#define ijkn(i,j,k,n) ((i)+(j)*nx+(k)*nxy+(n)*nxyz)

/* ------------------ readplot  ------------------------ */

void readplot(char *file, int ifile, int flag, int *errorcode){
  int n, nn, ntotal, i, nnn;
  float *udata, *vdata, *wdata, *sdata;
  float sum;
  int error;
  int ibar,jbar,kbar;
  mesh *meshi,*gbb,*gbi;
  plot3d *p;
  int nloaded=0;
  int nx, ny, nz;
  int pn;
  size_t plot3dfilelen;

  CheckMemory;
  STRCPY(FULLTITLE,TITLEBASE);
  *errorcode=0;

  p=plot3dinfo+ifile;
  if(flag==UNLOAD&&p->loaded==0)return;

  highlight_mesh=p->blocknumber;
  update_highlight_mesh();
  meshi=meshinfo+highlight_mesh;
  update_current_mesh(meshi);


  if(flag==UNLOAD){
    meshi->plot3dfilenum=-1;
  }
  else{
    pn = meshi->plot3dfilenum;
    if(pn!=-1){
      plot3dinfo[pn].loaded=0;
      plot3dinfo[pn].display=0;
    }
    meshi->plot3dfilenum=ifile;
  }

  FREEMEMORY(meshi->iqdata); FREEMEMORY(meshi->qdata);
  FREEMEMORY(meshi->yzcolorbase); FREEMEMORY(meshi->xzcolorbase); FREEMEMORY(meshi->xycolorbase); 
  FREEMEMORY(meshi->yzcolorfbase); FREEMEMORY(meshi->xzcolorfbase); FREEMEMORY(meshi->xycolorfbase);
  FREEMEMORY(meshi->yzcolortbase); FREEMEMORY(meshi->xzcolortbase); FREEMEMORY(meshi->xycolortbase);
  FREEMEMORY(meshi->dx_xy); FREEMEMORY(meshi->dy_xy); FREEMEMORY(meshi->dz_xy);
  FREEMEMORY(meshi->dx_xz); FREEMEMORY(meshi->dy_xz); FREEMEMORY(meshi->dz_xz);
  FREEMEMORY(meshi->dx_yz); FREEMEMORY(meshi->dy_yz); FREEMEMORY(meshi->dz_yz);

  freesurface(&meshi->currentsurf);
  freesurface(&meshi->currentsurf2);
  freecontour(&meshi->plot3dcontour1);
  freecontour(&meshi->plot3dcontour2);
  freecontour(&meshi->plot3dcontour3);
  initcontour(&meshi->plot3dcontour1,rgbptr,nrgb);
  initcontour(&meshi->plot3dcontour2,rgbptr,nrgb);
  initcontour(&meshi->plot3dcontour3,rgbptr,nrgb);


  for(i=0;i<nmeshes;i++){
    gbb=meshinfo+i;
    if(gbb->plot3dfilenum!=-1)nloaded++;
  }
  if(nloaded==0){
    if(colorlabelp3 != NULL){
      for(nn=0;nn<mxplot3dvars;nn++){
        for(nnn=0;nnn<MAXRGB;nnn++){
          FREEMEMORY((*(colorlabelp3+nn))[nnn]);
        }
        FREEMEMORY(colorlabelp3[nn]);
      }
      FREEMEMORY(colorlabelp3);
    }
    if(scalep3 != NULL){
      for(nn=0;nn<mxplot3dvars;nn++){FREEMEMORY(scalep3[nn]);}
      FREEMEMORY(scalep3);
    }
    if(p3levels!=NULL){
      for(nn=0;nn<mxplot3dvars;nn++){FREEMEMORY(p3levels[nn]);}
      FREEMEMORY(p3levels);
    }
    if(p3levels256!=NULL){
      for(nn=0;nn<mxplot3dvars;nn++){FREEMEMORY(p3levels256[nn]);}
      FREEMEMORY(p3levels256);
    }
  }

  if(flag==UNLOAD){
    p->loaded=0;
    p->display=0;
    plotstate=getplotstate(STATIC_PLOTS);
    meshi=meshinfo+p->blocknumber;
    meshi->plot3dfilenum=-1;
    if(nloaded==0){
      ReadPlot3dFile=0;
      numplot3dvars=0;
    }
    updatemenu=1;
#ifdef _DEBUG
    printf("After plot3d unload: ");
    PrintMemoryInfo;
#endif
    updatetimes();
    return;
  }
  if(ReadPlot3dFile==0){
    plotstate=getplotstate(STATIC_PLOTS);
    return;
  }

  ibar=meshi->ibar;
  jbar=meshi->jbar;
  kbar=meshi->kbar;
  nx = ibar+1; ny = jbar+1; nz = kbar+1;
  ntotal = nx*ny*nz;
  numplot3dvars=5;

  uindex = plot3dinfo[ifile].u;
  vindex = plot3dinfo[ifile].v;
  windex = plot3dinfo[ifile].w;
  if(uindex!=-1||vindex!=-1||windex!=-1)numplot3dvars=plot3dinfo[ifile].nvars;

  if(NewMemory((void **)&meshi->qdata,numplot3dvars*ntotal*sizeof(float))==0){
    *errorcode=1;
    readplot("",ifile,UNLOAD,&error);
    return;
  }

  if(NewMemory((void **)&meshi->yzcolorbase ,ny*nz*sizeof(int))==0||
     NewMemory((void **)&meshi->xzcolorbase ,nx*nz*sizeof(int))==0||
     NewMemory((void **)&meshi->xycolorbase ,nx*ny*sizeof(int))==0||
     NewMemory((void **)&meshi->yzcolorfbase,ny*nz*sizeof(int))==0||
     NewMemory((void **)&meshi->xzcolorfbase,nx*nz*sizeof(int))==0||
     NewMemory((void **)&meshi->xycolorfbase,nx*ny*sizeof(int))==0||
     NewMemory((void **)&meshi->yzcolortbase,ny*nz*sizeof(float))==0||
     NewMemory((void **)&meshi->xzcolortbase,nx*nz*sizeof(float))==0||
     NewMemory((void **)&meshi->xycolortbase,nx*ny*sizeof(float))==0||
     NewMemory((void **)&meshi->dx_xy       ,nx*ny*sizeof(int))==0||
     NewMemory((void **)&meshi->dy_xy       ,nx*ny*sizeof(int))==0||
     NewMemory((void **)&meshi->dz_xy       ,nx*ny*sizeof(int))==0||
     NewMemory((void **)&meshi->dx_xz       ,nx*nz*sizeof(int))==0||
     NewMemory((void **)&meshi->dy_xz       ,nx*nz*sizeof(int))==0||
     NewMemory((void **)&meshi->dz_xz       ,nx*nz*sizeof(int))==0||
     NewMemory((void **)&meshi->dx_yz       ,ny*nz*sizeof(int))==0||
     NewMemory((void **)&meshi->dy_yz       ,ny*nz*sizeof(int))==0||
     NewMemory((void **)&meshi->dz_yz       ,ny*nz*sizeof(int))==0){
     *errorcode=1;
     readplot("",ifile,UNLOAD,&error);
     return;
  }

  plot3dfilelen = strlen(file);
  printf("Loading plot3d data: %s\n",file);
    FORTgetplot3dq(file,&nx,&ny,&nz,meshi->qdata,&error,&endian,&isotest,plot3dfilelen);
  if(NewMemory((void **)&meshi->iqdata,numplot3dvars*ntotal*sizeof(int))==0){
    *errorcode=1;
    readplot("",ifile,UNLOAD,&error);
    return;
  }
  p->loaded=1;
  p->display=1;
  speedmax = -1.;
  if(uindex!=-1||vindex!=-1||windex!=-1){
    vectorspresent=1;
    p->nvars=mxplot3dvars;
    udata = meshi->qdata + ntotal*uindex;
    vdata = meshi->qdata + ntotal*vindex;
    wdata = meshi->qdata + ntotal*windex;
    sdata = meshi->qdata + ntotal*5;
    for(i=0;i<ntotal;i++){
      sum=0.0f;
      if(uindex!=-1){sum += *udata*(*udata);udata++;}
      if(vindex!=-1){sum += *vdata*(*vdata);vdata++;}
      if(windex!=-1){sum += *wdata*(*wdata);wdata++;}
      *sdata=sqrt((double)sum);
      sdata++; 
    }
  }                        


  if(NewMemory((void **)&colorlabelp3,mxplot3dvars*sizeof(char **))==0||
     NewMemory((void **)&scalep3     ,mxplot3dvars*sizeof(char *))==0){
    *errorcode=1;
    readplot("",ifile,UNLOAD,&error);
    return;
  }
  for(nn=0;nn<mxplot3dvars;nn++){
    colorlabelp3[nn]=NULL;
    scalep3[nn]=NULL;
  }
  for(nn=0;nn<mxplot3dvars;nn++){
    if(NewMemory((void **)&colorlabelp3[nn],MAXRGB*sizeof(char *))==0||
       NewMemory((void **)&scalep3[nn],30*sizeof(char))==0){
      *errorcode=1;
      readplot("",ifile,UNLOAD,&error);
      return;
    }
  }


  scalep3copy = scalep3;
  if(NewMemory((void **)&p3levels,mxplot3dvars*sizeof(float *))==0||
     NewMemory((void **)&p3levels256,mxplot3dvars*sizeof(float *))==0){
    *errorcode=1;
    readplot("",ifile,UNLOAD,&error);
    return;
  }
  for(nn=0;nn<mxplot3dvars;nn++){
    p3levels[nn]=NULL;
    p3levels256[nn]=NULL;
  }
  if(nloaded>1){
    for(nn=0;nn<numplot3dvars;nn++){
      if(setp3min[nn]!=SET_MIN&&setp3min[nn]!=CHOP_MIN){
        setp3min[nn]=SET_MIN;
      }
      if(setp3max[nn]!=SET_MAX&&setp3max[nn]!=CHOP_MAX){
        setp3max[nn]=SET_MAX;
      }
    }
  updateglui();
  }
  for (nn=0;nn<numplot3dvars;nn++){
    if(nplot3d>0){
      shortp3label[nn] = plot3dinfo[ifile].label[nn].shortlabel;
      unitp3label[nn] = plot3dinfo[ifile].label[nn].unit;
    }
    else{
      char numstring[4];
      
      sprintf(numstring,"%i",nn);
      strcpy(shortp3label[nn],numstring);
      unitp3label[nn] = blank;
    }

    for(n=0;n<MAXRGB;n++){(*(colorlabelp3+nn))[n]=NULL;}

    if(NewMemory((void **)&p3levels[nn],(nrgb+1)*sizeof(float))==0||
       NewMemory((void **)&p3levels256[nn],256*sizeof(float))==0){
      readplot("",ifile,UNLOAD,&error);
      *errorcode=1;
      return;
    }
    for(n=0;n<nrgb;n++){
      if(NewMemory((void **)&(*(colorlabelp3+nn))[n],11)==0){
        *errorcode=1;
        readplot("",ifile,UNLOAD,&error);
        return;
      }
    }
    getPlot3DColors(nn,
                  setp3min[nn],p3min+nn, setp3max[nn],p3max+nn, 
                  nrgb_full, nrgb, *(colorlabelp3+nn),scalep3copy,p3levels[nn],p3levels256[nn]);
    scalep3copy++;
  }
  if(meshi->plotx==-1)meshi->plotx=ibar/2; 
  if(meshi->ploty==-1)meshi->ploty=jbar/2; 
  if(meshi->plotz==-1)meshi->plotz=kbar/2;
  meshi->plot3d_speedmax=0.0f;
  if(uindex!=-1||vindex!=-1||windex!=-1||numplot3dvars>5)meshi->plot3d_speedmax=p3max[5];
  speedmax=-1000000.;
  for(i=0;i<nmeshes;i++){
    gbi=meshinfo+i;
    if(gbi->plot3dfilenum==-1)continue;
    if(speedmax<gbi->plot3d_speedmax)speedmax=gbi->plot3d_speedmax;
  }
  updateplotslice(1);
  updateplotslice(2);
  updateplotslice(3);
  visGrid=0;
  meshi->visInteriorPatches=0;
  if(meshi->visx==1){
    updateshowstep(1,DIRX);
  }
  if(meshi->visy==1){
    updateshowstep(1,DIRY);
  }
  if(meshi->visz==1){
    updateshowstep(1,DIRZ);
  }
  if(visiso==1){
    updatesurface();
  }
  STRCAT(FULLTITLE,", ");
  STRCAT(FULLTITLE,file);
  updateplot3dlistindex();
#ifdef _DEBUG
  printf("After plot3d load: ");
  PrintMemoryInfo;
#endif
  updatetimes();
  IDLE();

  GLUTPOSTREDISPLAY
}

/* ------------------ update_plot3dtitle ------------------------ */

void update_plot3dtitle(void){
  int filenum;
  plot3d *plot3di;
  mesh *meshi;

  STRCPY(FULLTITLE,TITLEBASE);
  meshi=current_mesh;
  if(meshi==NULL)meshi=meshinfo;
  filenum=meshi->plot3dfilenum;
  if(filenum!=-1){
    plot3di = plot3dinfo+meshi->plot3dfilenum;
    STRCAT(FULLTITLE,", ");
    STRCAT(FULLTITLE,plot3di->file);
  }
}
/* ------------------ drawplot3d_texture ------------------------ */

void drawplot3d_texture(mesh *meshi){
  int i,j,k;
  int colorindex;
  float *color1t, *color2t;
  float dx, dy, dz;
  int plotx, ploty, plotz;
  int visx, visy, visz;
  float *xplt, *yplt, *zplt;
  int ibar, jbar, kbar;
  isosurface *currentsurfptr,*currentsurf2ptr;
  contour *plot3dcontour1ptr, *plot3dcontour2ptr, *plot3dcontour3ptr;
  int *yzcolorbase, *xzcolorbase, *xycolorbase; 
  float *yzcolortbase, *xzcolortbase, *xycolortbase; 

  float *dx_xy, *dy_xy, *dz_xy;
  float *dx_xz, *dy_xz, *dz_xz;
  float *dx_yz, *dy_yz, *dz_yz;
  float *dx_yzcopy, *dy_yzcopy, *dz_yzcopy;
  float *dx_xzcopy, *dy_xzcopy, *dz_xzcopy;
  float *dx_xycopy, *dy_xycopy, *dz_xycopy;
  int nx, ny, nz,nxy;
  char *iblank_x, *iblank_y, *iblank_z, *iblank;
  float *vector_color;

  plotx = meshi->plotx;
  ploty = meshi->ploty;
  plotz = meshi->plotz;

  visx = meshi->visx;
  visy = meshi->visy;
  visz = meshi->visz;
  ibar = meshi->ibar;
  jbar = meshi->jbar;
  kbar = meshi->kbar;
  xplt = meshi->xplt;
  yplt = meshi->yplt;
  zplt = meshi->zplt;
  iblank_x = meshi->c_iblank_x;
  iblank_y = meshi->c_iblank_y;
  iblank_z = meshi->c_iblank_z;
  iblank = meshi->c_iblank;


  nx = ibar+1;
  ny = jbar+1;
  nz = kbar+1;
  nxy = nx*ny;

  currentsurfptr=&meshi->currentsurf;
  currentsurf2ptr=&meshi->currentsurf2;
  plot3dcontour1ptr=&meshi->plot3dcontour1;
  plot3dcontour2ptr=&meshi->plot3dcontour2;
  plot3dcontour3ptr=&meshi->plot3dcontour3;
  yzcolorbase=meshi->yzcolorbase;
  xzcolorbase=meshi->xzcolorbase;
  xycolorbase=meshi->xycolorbase;
  yzcolortbase=meshi->yzcolortbase;
  xzcolortbase=meshi->xzcolortbase;
  xycolortbase=meshi->xycolortbase;
  dx_xy=meshi->dx_xy;
  dx_xz=meshi->dx_xz;
  dx_yz=meshi->dx_yz;

  dy_xy=meshi->dy_xy;
  dy_xz=meshi->dy_xz;
  dy_yz=meshi->dy_yz;

  dz_xy=meshi->dz_xy;
  dz_xz=meshi->dz_xz;
  dz_yz=meshi->dz_yz;
  
  if(cullfaces==1)glDisable(GL_CULL_FACE);
  if(visiso==1){
    drawstaticiso(currentsurfptr,p3dsurfacetype,p3dsurfacesmooth,2,0);
    if(surfincrement!=0)drawstaticiso(currentsurf2ptr,p3dsurfacetype,p3dsurfacesmooth,2,0);
    if(visGrid!=0){
      if(transparentflag==1)transparentoff();
      if(cullfaces==1)glEnable(GL_CULL_FACE);
      return;
    }
  }

  if(transparentflag==1)transparenton();
  if(visVector==0&&p3cont2d==SHADED_CONTOURS){
    glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
    glEnable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D,texture_plot3d_colorbar_id);
  }

  /* +++++++++++++++++++++++++++   draw yz contours +++++++++++++++++++++++++++++++++++++ */

  if(visx!=0){
    if(visVector==0)DrawContours(plot3dcontour1ptr,p3cont2d);
    if(visVector==0&&p3cont2d==SHADED_CONTOURS){
      if(plotx<0){
        plotx=ibar;
        updateplotslice(1);
      }
      if(plotx>ibar){plotx=0;updateplotslice(1);}
      glBegin(GL_TRIANGLES);
      for(j=0; j<jbar; j++){
        color1t=yzcolortbase + j*nz;
        color2t=color1t+nz;
        for(k=0; k<kbar; k++){
          if(iblank_x[ijknode(plotx,j,k)]==2){
            if(abs(color1t[k]-color2t[k+1])<abs(color1t[k+1]-color2t[k])){
              glTexCoord1f(color1t[k]);  glVertex3f(xplt[plotx],yplt[j],zplt[k]);
              glTexCoord1f(color2t[k]);  glVertex3f(xplt[plotx],yplt[j+1],zplt[k]);
              glTexCoord1f(color2t[k+1]);glVertex3f(xplt[plotx],yplt[j+1],zplt[k+1]);

              glTexCoord1f(color1t[k]);  glVertex3f(xplt[plotx],yplt[j],zplt[k]);
              glTexCoord1f(color2t[k+1]);glVertex3f(xplt[plotx],yplt[j+1],zplt[k+1]);
              glTexCoord1f(color1t[k+1]);glVertex3f(xplt[plotx],yplt[j],zplt[k+1]);
            }
            else{
              glTexCoord1f(color1t[k]);  glVertex3f(xplt[plotx],yplt[j],zplt[k]);
              glTexCoord1f(color2t[k]);  glVertex3f(xplt[plotx],yplt[j+1],zplt[k]);
              glTexCoord1f(color1t[k+1]);glVertex3f(xplt[plotx],yplt[j],zplt[k+1]);

              glTexCoord1f(color2t[k]);  glVertex3f(xplt[plotx],yplt[j+1],zplt[k]);
              glTexCoord1f(color2t[k+1]);glVertex3f(xplt[plotx],yplt[j+1],zplt[k+1]);
              glTexCoord1f(color1t[k+1]);glVertex3f(xplt[plotx],yplt[j],zplt[k+1]);
            }
          }
        }
      }
      glEnd();
    } 
    /* draw yz vectors */

    if(visVector==1){
      int *yzcolor;

      yzcolor=yzcolorbase;
      dx_yzcopy=dx_yz; dy_yzcopy=dy_yz; dz_yzcopy=dz_yz;
      antialias(1);
      glLineWidth(vectorlinewidth);
      glBegin(GL_LINES);
      for(j=0; j<=jbar; j+=vectorskip){
        dx_yzcopy = dx_yz+j*(kbar+1);
        dy_yzcopy = dy_yz+j*(kbar+1);
        dz_yzcopy = dz_yz+j*(kbar+1);
        yzcolor = yzcolorbase + j*(kbar+1);
        for(k=0; k<=kbar; k+=vectorskip){
        colorindex=*yzcolor;
        vector_color = rgb_plot3d + 4*colorindex;
        if(iblank[ijknode(plotx,j,k)]==1&&vector_color[3]>0.5){
          glColor4fv(vector_color);
          dx=*dx_yzcopy/2.0;
          dy=*dy_yzcopy/2.0;
          dz=*dz_yzcopy/2.0;
          glVertex3f(xplt[plotx]-dx,yplt[j]-dy,zplt[k]-dz);
          glVertex3f(xplt[plotx]+dx,yplt[j]+dy,zplt[k]+dz);
        }
        dx_yzcopy+=vectorskip;
        dy_yzcopy+=vectorskip;
        dz_yzcopy+=vectorskip;
        yzcolor+=vectorskip;
      }} 
      glEnd();
      antialias(0);

      /* draw points for yz vectors */

      yzcolor=yzcolorbase;
      dx_yzcopy=dx_yz;
      dy_yzcopy=dy_yz;
      dz_yzcopy=dz_yz;
      glPointSize(vectorpointsize);
      glBegin(GL_POINTS);
      for(j=0; j<=jbar; j+=vectorskip){
        dx_yzcopy = dx_yz+j*(kbar+1);
        dy_yzcopy = dy_yz+j*(kbar+1);
        dz_yzcopy = dz_yz+j*(kbar+1);
        yzcolor = yzcolorbase + j*(kbar+1);
        for(k=0; k<=kbar; k+=vectorskip){
        colorindex=*yzcolor;
        vector_color = rgb_plot3d + 4*colorindex;
        if(iblank[ijknode(plotx,j,k)]==1&&vector_color[3]>0.5){
          glColor4fv(vector_color);
          glVertex3f(
            xplt[plotx]+*dx_yzcopy/(float)2.0,
            yplt[j]+*dy_yzcopy/(float)2.0,
            zplt[k]+*dz_yzcopy/(float)2.0
            );
        }
        dx_yzcopy+=vectorskip;
        dy_yzcopy+=vectorskip;
        dz_yzcopy+=vectorskip;
        yzcolor+=vectorskip;
      }} 
      glEnd();
    }
  }

  /* +++++++++++++++++++++++++++++++++  draw xz contours  ++++++++++++++++++++++++++++++++++++++++ */

  if(visy!=0){
    if(visVector==0){
      DrawContours(plot3dcontour2ptr,p3cont2d);
    }
    if(visVector==0&&p3cont2d==SHADED_CONTOURS){
      /*
      if(ploty<0){
        ploty=jbar;
        updateplotslice(2);
      }
      if(ploty>jbar){
        ploty=0;
        updateplotslice(2);
      }
      */
        glBegin(GL_TRIANGLES);
        for(i=0; i<ibar; i++){
          color1t=xzcolortbase + i*nz;
          color2t=color1t+nz;
          for(k=0; k<kbar; k++){
            if(iblank_y[ijknode(i,ploty,k)]==2){
              if(abs(color1t[k]-color2t[k+1])<abs(color1t[k+1]-color2t[k])){
                glTexCoord1f(color1t[k]);  glVertex3f(xplt[i],yplt[ploty],zplt[k]);
                glTexCoord1f(color2t[k]);  glVertex3f(xplt[i+1],yplt[ploty],zplt[k]);
                glTexCoord1f(color2t[k+1]);glVertex3f(xplt[i+1],yplt[ploty],zplt[k+1]);
  
                glTexCoord1f(color1t[k]);  glVertex3f(xplt[i],yplt[ploty],zplt[k]);
                glTexCoord1f(color2t[k+1]);glVertex3f(xplt[i+1],yplt[ploty],zplt[k+1]);
                glTexCoord1f(color1t[k+1]);glVertex3f(xplt[i],yplt[ploty],zplt[k+1]);
              }
              else{
                glTexCoord1f(color1t[k]);  glVertex3f(xplt[i],yplt[ploty],zplt[k]);
                glTexCoord1f(color2t[k]);  glVertex3f(xplt[i+1],yplt[ploty],zplt[k]);
                glTexCoord1f(color1t[k+1]);glVertex3f(xplt[i],yplt[ploty],zplt[k+1]);
              
                glTexCoord1f(color2t[k]);  glVertex3f(xplt[i+1],yplt[ploty],zplt[k]);
                glTexCoord1f(color2t[k+1]);glVertex3f(xplt[i+1],yplt[ploty],zplt[k+1]);
                glTexCoord1f(color1t[k+1]);glVertex3f(xplt[i],yplt[ploty],zplt[k+1]);
              }
            }   
          }      
        }         
        glEnd();
    } 

    /* draw xz vectors */

    if(visVector==1){
      int *xzcolor;

      xzcolor=xzcolorbase;
      antialias(1);
      glLineWidth(vectorlinewidth);
      glBegin(GL_LINES);
      for(i=0; i<=ibar; i+=vectorskip){
        dx_xzcopy=dx_xz+i*(kbar+1);
        dy_xzcopy=dy_xz+i*(kbar+1);
        dz_xzcopy=dz_xz+i*(kbar+1);
        xzcolor = xzcolorbase + i*(kbar+1);
        for(k=0; k<=kbar; k+=vectorskip){
        colorindex=*xzcolor;
        vector_color = rgb_plot3d + 4*colorindex;
        if(iblank[ijknode(i,ploty,k)]==1&&vector_color[3]>0.5){
          glColor4fv(vector_color);
          dx=*dx_xzcopy/2.0;
          dy=*dy_xzcopy/2.0;
          dz=*dz_xzcopy/2.0;
          glVertex3f(xplt[i]-dx,yplt[ploty]-dy,zplt[k]-dz);
          glVertex3f(xplt[i]+dx,yplt[ploty]+dy,zplt[k]+dz);
        }
        dx_xzcopy+=vectorskip;
        dy_xzcopy+=vectorskip;
        dz_xzcopy+=vectorskip;
        xzcolor+=vectorskip;
      }} 
      glEnd();
      antialias(0);

      /* draw points for xz vectors */

      xzcolor=xzcolorbase;
      glPointSize(vectorpointsize);
      glBegin(GL_POINTS);
      for(i=0; i<=ibar; i+=vectorskip){
        dx_xzcopy=dx_xz+i*(kbar+1);
        dy_xzcopy=dy_xz+i*(kbar+1);
        dz_xzcopy=dz_xz+i*(kbar+1);
        xzcolor = xzcolorbase + i*(kbar+1);
        for(k=0; k<=kbar; k+=vectorskip){
          colorindex=*xzcolor;
          vector_color = rgb_plot3d + 4*colorindex;
          if(iblank[ijknode(i,ploty,k)]==1&&vector_color[3]>0.5){
            glColor4fv(vector_color);
            dx=*dx_xzcopy/2.0;
            dy=*dy_xzcopy/2.0;
            dz=*dz_xzcopy/2.0;
            glVertex3f(xplt[i]+dx,yplt[ploty]+dy,zplt[k]+dz);
          }
          dx_xzcopy+=vectorskip;
          dy_xzcopy+=vectorskip;
          dz_xzcopy+=vectorskip;
          xzcolor+=vectorskip;
        }
      } 
      glEnd();
    }
  }

  /* ++++++++++++++++++++++++++++ draw xy contours ++++++++++++++++++++++++++++++++ */

  if(visz!=0){
    if(visVector==0){
      DrawContours(plot3dcontour3ptr,p3cont2d);
    }
    if(visVector==0&&p3cont2d==SHADED_CONTOURS){
      if(plotz<0){plotz=kbar;updateplotslice(3);}
      if(plotz>kbar){plotz=0;updateplotslice(3);}
      glBegin(GL_TRIANGLES);
      for(i=0; i<ibar; i++){
        color1t=xycolortbase + i*ny;
        color2t=color1t+ny;
        for(j=0; j<jbar; j++){
          if(iblank_z[ijknode(i,j,plotz)]==2){
            if(abs(color1t[j]-color2t[j+1])<abs(color1t[j+1]-color2t[j])){
              glTexCoord1f(  color1t[j]);glVertex3f(  xplt[i],  yplt[j],zplt[plotz]);
              glTexCoord1f(  color2t[j]);glVertex3f(xplt[i+1],  yplt[j],zplt[plotz]);
              glTexCoord1f(color2t[j+1]);glVertex3f(xplt[i+1],yplt[j+1],zplt[plotz]);

              glTexCoord1f(  color1t[j]);glVertex3f(  xplt[i],  yplt[j],zplt[plotz]);
              glTexCoord1f(color2t[j+1]);glVertex3f(xplt[i+1],yplt[j+1],zplt[plotz]);
              glTexCoord1f(color1t[j+1]);glVertex3f(  xplt[i],yplt[j+1],zplt[plotz]);
            }
            else{
              glTexCoord1f(  color1t[j]);glVertex3f(  xplt[i],  yplt[j],zplt[plotz]);
              glTexCoord1f(  color2t[j]);glVertex3f(xplt[i+1],  yplt[j],zplt[plotz]);
              glTexCoord1f(color1t[j+1]);glVertex3f(  xplt[i],yplt[j+1],zplt[plotz]);

              glTexCoord1f(  color2t[j]);glVertex3f(xplt[i+1],  yplt[j],zplt[plotz]);
              glTexCoord1f(color2t[j+1]);glVertex3f(xplt[i+1],yplt[j+1],zplt[plotz]);
              glTexCoord1f(color1t[j+1]);glVertex3f(  xplt[i],yplt[j+1],zplt[plotz]);
            }
          }
        }
      }
      glEnd();
    } 

    /* draw xy vectors */

    if(visVector==1){
      int *xycolor;

      xycolor=xycolorbase;
      antialias(1);
      glLineWidth(vectorlinewidth);
      glBegin(GL_LINES);
      for(i=0; i<=ibar; i+=vectorskip){
        dx_xycopy=dx_xy+i*(jbar+1);
        dy_xycopy=dy_xy+i*(jbar+1);
        dz_xycopy=dz_xy+i*(jbar+1);
        xycolor = xycolorbase + i*(jbar+1);
        for(j=0; j<=jbar; j+=vectorskip){
        colorindex=*xycolor;
        vector_color = rgb_plot3d + 4*colorindex;
        if(iblank[ijknode(i,j,plotz)]==1&&vector_color[3]>0.5){
          glColor4fv(vector_color);
          dx=*dx_xycopy/2.0;
          dy=*dy_xycopy/2.0;
          dz=*dz_xycopy/2.0;
          glVertex3f(xplt[i]-dx,yplt[j]-dy,zplt[plotz]-dz);
          glVertex3f(xplt[i]+dx,yplt[j]+dy,zplt[plotz]+dz);
        }
        dx_xycopy+=vectorskip;
        dy_xycopy+=vectorskip;
        dz_xycopy+=vectorskip;
        xycolor+=vectorskip;
      }} 
      glEnd();
      antialias(0);

      /* draw points for xy vectors */

      xycolor=xycolorbase;
      glPointSize(vectorpointsize);
      glBegin(GL_POINTS);
      for(i=0; i<=ibar; i+=vectorskip){
        dx_xycopy=dx_xy+i*(jbar+1);
        dy_xycopy=dy_xy+i*(jbar+1);
        dz_xycopy=dz_xy+i*(jbar+1);
        xycolor = xycolorbase + i*(jbar+1);
        for(j=0; j<=jbar; j+=vectorskip){
        colorindex=*xycolor;
        vector_color = rgb_plot3d + 4*colorindex;
        if(iblank[ijknode(i,j,plotz)]==1&&vector_color[3]>0.5){
          glColor4fv(vector_color);
          glVertex3f(
            xplt[i]+*dx_xycopy/(float)2.0,
            yplt[j]+*dy_xycopy/(float)2.0,
            zplt[plotz]+*dz_xycopy/(float)2.0
            );
        }
        dx_xycopy+=vectorskip;
        dy_xycopy+=vectorskip;
        dz_xycopy+=vectorskip;
        xycolor+=vectorskip;
      }} 
      glEnd();
    }
  }
  if(visVector==0&&p3cont2d==SHADED_CONTOURS){
    glDisable(GL_TEXTURE_1D);
  }
  if(transparentflag==1)transparentoff();
  if(cullfaces==1)glEnable(GL_CULL_FACE);

}

/* ------------------ drawplot3d ------------------------ */

void drawplot3d(mesh *meshi){
  int i,j,k;
  int colorindex;
  int *color1, *color2;
  float dx, dy, dz;
  int plotx, ploty, plotz;
  int visx, visy, visz;
  float *xplt, *yplt, *zplt;
  int ibar, jbar, kbar;
  isosurface *currentsurfptr,*currentsurf2ptr;
  contour *plot3dcontour1ptr, *plot3dcontour2ptr, *plot3dcontour3ptr;
  int *yzcolorbase, *xzcolorbase, *xycolorbase; 
  float *dx_xy, *dy_xy, *dz_xy;
  float *dx_xz, *dy_xz, *dz_xz;
  float *dx_yz, *dy_yz, *dz_yz;
  float *dx_yzcopy, *dy_yzcopy, *dz_yzcopy;
  float *dx_xzcopy, *dy_xzcopy, *dz_xzcopy;
  float *dx_xycopy, *dy_xycopy, *dz_xycopy;
  int nx, ny, nz,nxy;
  char *iblank_x, *iblank_y, *iblank_z, *iblank;
  float *vector_color;

  plotx = meshi->plotx;
  ploty = meshi->ploty;
  plotz = meshi->plotz;

  visx = meshi->visx;
  visy = meshi->visy;
  visz = meshi->visz;
  ibar = meshi->ibar;
  jbar = meshi->jbar;
  kbar = meshi->kbar;
  xplt = meshi->xplt;
  yplt = meshi->yplt;
  zplt = meshi->zplt;
  iblank_x = meshi->c_iblank_x;
  iblank_y = meshi->c_iblank_y;
  iblank_z = meshi->c_iblank_z;
  iblank = meshi->c_iblank;


  nx = ibar+1;
  ny = jbar+1;
  nz = kbar+1;
  nxy = nx*ny;

  currentsurfptr=&meshi->currentsurf;
  currentsurf2ptr=&meshi->currentsurf2;
  plot3dcontour1ptr=&meshi->plot3dcontour1;
  plot3dcontour2ptr=&meshi->plot3dcontour2;
  plot3dcontour3ptr=&meshi->plot3dcontour3;
  yzcolorbase=meshi->yzcolorbase;
  xzcolorbase=meshi->xzcolorbase;
  xycolorbase=meshi->xycolorbase;
  dx_xy=meshi->dx_xy;
  dx_xz=meshi->dx_xz;
  dx_yz=meshi->dx_yz;

  dy_xy=meshi->dy_xy;
  dy_xz=meshi->dy_xz;
  dy_yz=meshi->dy_yz;

  dz_xy=meshi->dz_xy;
  dz_xz=meshi->dz_xz;
  dz_yz=meshi->dz_yz;
  
  if(cullfaces==1)glDisable(GL_CULL_FACE);
  if(visiso==1){
    drawstaticiso(currentsurfptr,p3dsurfacetype,p3dsurfacesmooth,2,0);
    if(surfincrement!=0)drawstaticiso(currentsurf2ptr,p3dsurfacetype,p3dsurfacesmooth,2,0);
    if(visGrid!=0){
      if(transparentflag==1)transparentoff();
      if(cullfaces==1)glEnable(GL_CULL_FACE);
      return;
    }
  }

  if(transparentflag==1)transparenton();

  /* +++++++++++++++++++++++++++   draw yz contours +++++++++++++++++++++++++++++++++++++ */

  if(visx!=0){
    if(visVector==0)DrawContours(plot3dcontour1ptr,p3cont2d);
    if(visVector==0&&p3cont2d==SHADED_CONTOURS){
      if(plotx<0){
        plotx=ibar;
        updateplotslice(1);
      }
      if(plotx>ibar){plotx=0;updateplotslice(1);}
      glBegin(GL_TRIANGLES);
      for(j=0; j<jbar; j++){
        color1=yzcolorbase + j*nz;
        color2=color1+nz;
        for(k=0; k<kbar; k++){
          if(iblank_x[ijknode(plotx,j,k)]==2){
            if(abs(color1[k]-color2[k+1])<abs(color1[k+1]-color2[k])){
              glColor4fv(rgb_plot3d+4*color1[k]);  glVertex3f(xplt[plotx],yplt[j],zplt[k]);
              glColor4fv(rgb_plot3d+4*color2[k]);  glVertex3f(xplt[plotx],yplt[j+1],zplt[k]);
              glColor4fv(rgb_plot3d+4*color2[k+1]);glVertex3f(xplt[plotx],yplt[j+1],zplt[k+1]);

              glColor4fv(rgb_plot3d+4*color1[k]);  glVertex3f(xplt[plotx],yplt[j],zplt[k]);
              glColor4fv(rgb_plot3d+4*color2[k+1]);glVertex3f(xplt[plotx],yplt[j+1],zplt[k+1]);
              glColor4fv(rgb_plot3d+4*color1[k+1]);glVertex3f(xplt[plotx],yplt[j],zplt[k+1]);
            }
            else{
              glColor4fv(rgb_plot3d+4*color1[k]);  glVertex3f(xplt[plotx],yplt[j],zplt[k]);
              glColor4fv(rgb_plot3d+4*color2[k]);  glVertex3f(xplt[plotx],yplt[j+1],zplt[k]);
              glColor4fv(rgb_plot3d+4*color1[k+1]);glVertex3f(xplt[plotx],yplt[j],zplt[k+1]);

              glColor4fv(rgb_plot3d+4*color2[k]);  glVertex3f(xplt[plotx],yplt[j+1],zplt[k]);
              glColor4fv(rgb_plot3d+4*color2[k+1]);glVertex3f(xplt[plotx],yplt[j+1],zplt[k+1]);
              glColor4fv(rgb_plot3d+4*color1[k+1]);glVertex3f(xplt[plotx],yplt[j],zplt[k+1]);
            }
          }
        }
      }
      glEnd();
    } 
    /* draw yz vectors */

    if(visVector==1){
      int *yzcolor;

      yzcolor=yzcolorbase;
      dx_yzcopy=dx_yz; dy_yzcopy=dy_yz; dz_yzcopy=dz_yz;
      antialias(1);
      glLineWidth(vectorlinewidth);
      glBegin(GL_LINES);
      for(j=0; j<=jbar; j+=vectorskip){
        dx_yzcopy = dx_yz+j*(kbar+1);
        dy_yzcopy = dy_yz+j*(kbar+1);
        dz_yzcopy = dz_yz+j*(kbar+1);
        yzcolor = yzcolorbase + j*(kbar+1);
        for(k=0; k<=kbar; k+=vectorskip){
        colorindex=*yzcolor;
        vector_color = rgb_plot3d + 4*colorindex;
        if(iblank[ijknode(plotx,j,k)]==1&&vector_color[3]>0.5){
          glColor4fv(vector_color);
          dx=*dx_yzcopy/2.0;
          dy=*dy_yzcopy/2.0;
          dz=*dz_yzcopy/2.0;
          glVertex3f(xplt[plotx]-dx,yplt[j]-dy,zplt[k]-dz);
          glVertex3f(xplt[plotx]+dx,yplt[j]+dy,zplt[k]+dz);
        }
        dx_yzcopy+=vectorskip;
        dy_yzcopy+=vectorskip;
        dz_yzcopy+=vectorskip;
        yzcolor+=vectorskip;
      }} 
      glEnd();
      antialias(0);

      /* draw points for yz vectors */

      yzcolor=yzcolorbase;
      dx_yzcopy=dx_yz;
      dy_yzcopy=dy_yz;
      dz_yzcopy=dz_yz;
      glPointSize(vectorpointsize);
      glBegin(GL_POINTS);
      for(j=0; j<=jbar; j+=vectorskip){
        dx_yzcopy = dx_yz+j*(kbar+1);
        dy_yzcopy = dy_yz+j*(kbar+1);
        dz_yzcopy = dz_yz+j*(kbar+1);
        yzcolor = yzcolorbase + j*(kbar+1);
        for(k=0; k<=kbar; k+=vectorskip){
        colorindex=*yzcolor;
        vector_color = rgb_plot3d + 4*colorindex;
        if(iblank[ijknode(plotx,j,k)]==1&&vector_color[3]>0.5){
          glColor4fv(vector_color);
          glVertex3f(
            xplt[plotx]+*dx_yzcopy/(float)2.0,
            yplt[j]+*dy_yzcopy/(float)2.0,
            zplt[k]+*dz_yzcopy/(float)2.0
            );
        }
        dx_yzcopy+=vectorskip;
        dy_yzcopy+=vectorskip;
        dz_yzcopy+=vectorskip;
        yzcolor+=vectorskip;
      }} 
      glEnd();
    }
  }

  /* +++++++++++++++++++++++++++++++++  draw xz contours  ++++++++++++++++++++++++++++++++++++++++ */

  if(visy!=0){
    if(visVector==0){
      DrawContours(plot3dcontour2ptr,p3cont2d);
    }
    if(visVector==0&&p3cont2d==SHADED_CONTOURS){
      /*
      if(ploty<0){
        ploty=jbar;
        updateplotslice(2);
      }
      if(ploty>jbar){
        ploty=0;
        updateplotslice(2);
      }
      */
        glBegin(GL_TRIANGLES);
        for(i=0; i<ibar; i++){
          color1=xzcolorbase + i*nz;
          color2=color1+nz;
          for(k=0; k<kbar; k++){
            if(iblank_y[ijknode(i,ploty,k)]==2){
              if(abs(color1[k]-color2[k+1])<abs(color1[k+1]-color2[k])){
                glColor4fv(rgb_plot3d+4*color1[k]);  glVertex3f(xplt[i],yplt[ploty],zplt[k]);
                glColor4fv(rgb_plot3d+4*color2[k]);  glVertex3f(xplt[i+1],yplt[ploty],zplt[k]);
                glColor4fv(rgb_plot3d+4*color2[k+1]);glVertex3f(xplt[i+1],yplt[ploty],zplt[k+1]);
  
                glColor4fv(rgb_plot3d+4*color1[k]);  glVertex3f(xplt[i],yplt[ploty],zplt[k]);
                glColor4fv(rgb_plot3d+4*color2[k+1]);glVertex3f(xplt[i+1],yplt[ploty],zplt[k+1]);
                glColor4fv(rgb_plot3d+4*color1[k+1]);glVertex3f(xplt[i],yplt[ploty],zplt[k+1]);
              }
              else{
                glColor4fv(rgb_plot3d+4*color1[k]);  glVertex3f(xplt[i],yplt[ploty],zplt[k]);
                glColor4fv(rgb_plot3d+4*color2[k]);  glVertex3f(xplt[i+1],yplt[ploty],zplt[k]);
                glColor4fv(rgb_plot3d+4*color1[k+1]);glVertex3f(xplt[i],yplt[ploty],zplt[k+1]);
              
                glColor4fv(rgb_plot3d+4*color2[k]);  glVertex3f(xplt[i+1],yplt[ploty],zplt[k]);
                glColor4fv(rgb_plot3d+4*color2[k+1]);glVertex3f(xplt[i+1],yplt[ploty],zplt[k+1]);
                glColor4fv(rgb_plot3d+4*color1[k+1]);glVertex3f(xplt[i],yplt[ploty],zplt[k+1]);
              }
            }   
          }      
        }         
        glEnd();
    } 

    /* draw xz vectors */

    if(visVector==1){
      int *xzcolor;

      xzcolor=xzcolorbase;
      antialias(1);
      glLineWidth(vectorlinewidth);
      glBegin(GL_LINES);
      for(i=0; i<=ibar; i+=vectorskip){
        dx_xzcopy=dx_xz+i*(kbar+1);
        dy_xzcopy=dy_xz+i*(kbar+1);
        dz_xzcopy=dz_xz+i*(kbar+1);
        xzcolor = xzcolorbase + i*(kbar+1);
        for(k=0; k<=kbar; k+=vectorskip){
        colorindex=*xzcolor;
        vector_color = rgb_plot3d + 4*colorindex;
        if(iblank[ijknode(i,ploty,k)]==1&&vector_color[3]>0.5){
          glColor4fv(vector_color);
          dx=*dx_xzcopy/2.0;
          dy=*dy_xzcopy/2.0;
          dz=*dz_xzcopy/2.0;
          glVertex3f(xplt[i]-dx,yplt[ploty]-dy,zplt[k]-dz);
          glVertex3f(xplt[i]+dx,yplt[ploty]+dy,zplt[k]+dz);
        }
        dx_xzcopy+=vectorskip;
        dy_xzcopy+=vectorskip;
        dz_xzcopy+=vectorskip;
        xzcolor+=vectorskip;
      }} 
      glEnd();
      antialias(0);

      /* draw points for xz vectors */

      xzcolor=xzcolorbase;
      glPointSize(vectorpointsize);
      glBegin(GL_POINTS);
      for(i=0; i<=ibar; i+=vectorskip){
        dx_xzcopy=dx_xz+i*(kbar+1);
        dy_xzcopy=dy_xz+i*(kbar+1);
        dz_xzcopy=dz_xz+i*(kbar+1);
        xzcolor = xzcolorbase + i*(kbar+1);
        for(k=0; k<=kbar; k+=vectorskip){
          colorindex=*xzcolor;
          vector_color = rgb_plot3d + 4*colorindex;
          if(iblank[ijknode(i,ploty,k)]==1&&vector_color[3]>0.5){
            glColor4fv(vector_color);
            dx=*dx_xzcopy/2.0;
            dy=*dy_xzcopy/2.0;
            dz=*dz_xzcopy/2.0;
            glVertex3f(xplt[i]+dx,yplt[ploty]+dy,zplt[k]+dz);
          }
          dx_xzcopy+=vectorskip;
          dy_xzcopy+=vectorskip;
          dz_xzcopy+=vectorskip;
          xzcolor+=vectorskip;
        }
      } 
      glEnd();
    }
  }

  /* ++++++++++++++++++++++++++++ draw xy contours ++++++++++++++++++++++++++++++++ */

  if(visz!=0){
    if(visVector==0){
      DrawContours(plot3dcontour3ptr,p3cont2d);
    }
    if(visVector==0&&p3cont2d==SHADED_CONTOURS){
      if(plotz<0){plotz=kbar;updateplotslice(3);}
      if(plotz>kbar){plotz=0;updateplotslice(3);}
      glBegin(GL_TRIANGLES);
      for(i=0; i<ibar; i++){
        color1=xycolorbase + i*ny;
        color2=color1+ny;
        for(j=0; j<jbar; j++){
          if(iblank_z[ijknode(i,j,plotz)]==2){
            if(abs(color1[j]-color2[j+1])<abs(color1[j+1]-color2[j])){
              glColor4fv(rgb_plot3d+4*color1[j]);glVertex3f(xplt[i],yplt[j],zplt[plotz]);
              glColor4fv(rgb_plot3d+4*color2[j]);glVertex3f(xplt[i+1],yplt[j],zplt[plotz]);
              glColor4fv(rgb_plot3d+4*color2[j+1]);glVertex3f(xplt[i+1],yplt[j+1],zplt[plotz]);

              glColor4fv(rgb_plot3d+4*color1[j]);glVertex3f(xplt[i],yplt[j],zplt[plotz]);
              glColor4fv(rgb_plot3d+4*color2[j+1]);glVertex3f(xplt[i+1],yplt[j+1],zplt[plotz]);
              glColor4fv(rgb_plot3d+4*color1[j+1]);glVertex3f(xplt[i],yplt[j+1],zplt[plotz]);
            }
            else{
              glColor4fv(rgb_plot3d+4*color1[j]);glVertex3f(xplt[i],yplt[j],zplt[plotz]);
              glColor4fv(rgb_plot3d+4*color2[j]);glVertex3f(xplt[i+1],yplt[j],zplt[plotz]);
              glColor4fv(rgb_plot3d+4*color1[j+1]);glVertex3f(xplt[i],yplt[j+1],zplt[plotz]);

              glColor4fv(rgb_plot3d+4*color2[j]);glVertex3f(xplt[i+1],yplt[j],zplt[plotz]);
              glColor4fv(rgb_plot3d+4*color2[j+1]);glVertex3f(xplt[i+1],yplt[j+1],zplt[plotz]);
              glColor4fv(rgb_plot3d+4*color1[j+1]);glVertex3f(xplt[i],yplt[j+1],zplt[plotz]);
            }
          }
        }
      }
      glEnd();
    } 

    /* draw xy vectors */

    if(visVector==1){
      int *xycolor;

      xycolor=xycolorbase;
      antialias(1);
      glLineWidth(vectorlinewidth);
      glBegin(GL_LINES);
      for(i=0; i<=ibar; i+=vectorskip){
        dx_xycopy=dx_xy+i*(jbar+1);
        dy_xycopy=dy_xy+i*(jbar+1);
        dz_xycopy=dz_xy+i*(jbar+1);
        xycolor = xycolorbase + i*(jbar+1);
        for(j=0; j<=jbar; j+=vectorskip){
        colorindex=*xycolor;
        vector_color = rgb_plot3d + 4*colorindex;
        if(iblank[ijknode(i,j,plotz)]==1&&vector_color[3]>0.5){
          glColor4fv(vector_color);
          dx=*dx_xycopy/2.0;
          dy=*dy_xycopy/2.0;
          dz=*dz_xycopy/2.0;
          glVertex3f(xplt[i]-dx,yplt[j]-dy,zplt[plotz]-dz);
          glVertex3f(xplt[i]+dx,yplt[j]+dy,zplt[plotz]+dz);
        }
        dx_xycopy+=vectorskip;
        dy_xycopy+=vectorskip;
        dz_xycopy+=vectorskip;
        xycolor+=vectorskip;
      }} 
      glEnd();
      antialias(0);

      /* draw points for xy vectors */

      xycolor=xycolorbase;
      glPointSize(vectorpointsize);
      glBegin(GL_POINTS);
      for(i=0; i<=ibar; i+=vectorskip){
        dx_xycopy=dx_xy+i*(jbar+1);
        dy_xycopy=dy_xy+i*(jbar+1);
        dz_xycopy=dz_xy+i*(jbar+1);
        xycolor = xycolorbase + i*(jbar+1);
        for(j=0; j<=jbar; j+=vectorskip){
        colorindex=*xycolor;
        vector_color = rgb_plot3d + 4*colorindex;
        if(iblank[ijknode(i,j,plotz)]==1&&vector_color[3]>0.5){
          glColor4fv(vector_color);
          glVertex3f(
            xplt[i]+*dx_xycopy/(float)2.0,
            yplt[j]+*dy_xycopy/(float)2.0,
            zplt[plotz]+*dz_xycopy/(float)2.0
            );
        }
        dx_xycopy+=vectorskip;
        dy_xycopy+=vectorskip;
        dz_xycopy+=vectorskip;
        xycolor+=vectorskip;
      }} 
      glEnd();
    }
  }
  if(transparentflag==1)transparentoff();
  if(cullfaces==1)glEnable(GL_CULL_FACE);

}


/* ------------------ updatesurface ------------------------ */

void updatesurface(void){
  int colorindex,colorindex2;
  float level,level2;
  int ibar, jbar, kbar;
  float *xplt, *yplt, *zplt;
  isosurface *currentsurfptr,*currentsurf2ptr;
  float *qdata;
  char *iblank_cell;
  mesh *meshi;
  int plot3dsize;
  int i;

  for(i=0;i<nmeshes;i++){
    meshi = meshinfo+i;
    if(meshi->plot3dfilenum==-1)continue;

    ibar=meshi->ibar;
    jbar=meshi->jbar;
    kbar=meshi->kbar;
    plot3dsize=(ibar+1)*(jbar+1)*(kbar+1);
    xplt=meshi->xplt;
    yplt=meshi->yplt;
    zplt=meshi->zplt;
    iblank_cell=meshi->c_iblank_cell;

    currentsurfptr = &meshi->currentsurf;
    currentsurf2ptr = &meshi->currentsurf2;
    qdata=meshi->qdata;

    if(ReadPlot3dFile!=1)return;
    if(plotiso[plotn-1]<0){plotiso[plotn-1]=nrgb-2;}
    if(plotiso[plotn-1]>nrgb-2){plotiso[plotn-1]=0;}
    colorindex=plotiso[plotn-1];
    level = p3min[plotn-1] + colorindex*(p3max[plotn-1]-p3min[plotn-1])/((float)nrgb-2.0f);
    isolevelindex=colorindex;
    isolevelindex2=colorindex;
    freesurface(currentsurfptr);
    InitIsosurface(currentsurfptr, level, rgb[colorindex],-999);
    GetIsosurface(currentsurfptr,qdata+(plotn-1)*plot3dsize,NULL,iblank_cell,level,
      xplt,ibar+1,yplt,jbar+1,zplt,kbar+1,isooffset);
    GetNormalSurface(currentsurfptr);
    CompressIsosurface(currentsurfptr,1,
        xplt[0],xplt[ibar],
        yplt[0],yplt[jbar],
        zplt[0],zplt[kbar]);
    SmoothIsoSurface(currentsurfptr);

    if(surfincrement!=0){
      colorindex2=colorindex+surfincrement;
      if(colorindex2<0)colorindex2=nrgb-2;
      if(colorindex2>nrgb-2)colorindex2=0;
      level2 = p3min[plotn-1] + colorindex2*(p3max[plotn-1]-p3min[plotn-1])/((float)nrgb-2.0f);
      freesurface(currentsurf2ptr);
      InitIsosurface(currentsurf2ptr, level2, rgb[colorindex2],-999);
      GetIsosurface(currentsurf2ptr,qdata+(plotn-1)*plot3dsize,NULL,iblank_cell,level2,
        xplt,ibar+1,yplt,jbar+1,zplt,kbar+1,isooffset);
      GetNormalSurface(currentsurf2ptr);
      CompressIsosurface(currentsurf2ptr,1,
        xplt[0],xplt[ibar],
        yplt[0],yplt[jbar],
        zplt[0],zplt[kbar]);
      SmoothIsoSurface(currentsurf2ptr);

      isolevelindex2=colorindex2;
    }
  }
}

/* ------------------ updateallplotslices ------------------------ */

void updateallplotslices(void){
  int i;
  plot3d *p;
  mesh *gbsave;
  gbsave=current_mesh;
  for(i=0;i<nplot3d;i++){
    p=plot3dinfo+i;
    if(p->loaded==0)continue;
    update_current_mesh(meshinfo+p->blocknumber);
    if(current_mesh->visx==1)updateplotslice(1);
    if(current_mesh->visy==1)updateplotslice(2);
    if(current_mesh->visz==1)updateplotslice(3);
  }
  update_current_mesh(gbsave);
    
}

/* ------------------ update_plot_xyz ------------------------ */

void update_plot_xyz(mesh *current_mesh){
  int i;
  float xval, yval, zval;

  xval = current_mesh->xplt[current_mesh->plotx];
  yval = current_mesh->yplt[current_mesh->ploty];
  zval = current_mesh->zplt[current_mesh->plotz];

  for(i=0;i<nmeshes;i++){
    mesh *meshi;
    float xmin, xmax;
    float ymin, ymax;
    float zmin, zmax;

    meshi=meshinfo+i;
    if(meshi==current_mesh)continue;

    xmin = meshi->xplt[0];
    xmax = meshi->xplt[meshi->ibar];
    if(xmin<=xval&&xval<=xmax){
      int j;
      int iimin;
      float xxmin,xxdif;

      iimin=0;
      xxdif = xval - meshi->xplt[0];
      if(xxdif<0)xxdif=-xxdif;
      xxmin = xxdif;

      for(j=1;j<=meshi->ibar;j++){
        xxdif = xval - meshi->xplt[j];
        if(xxdif<0)xxdif=-xxdif;
        if(xxdif<xxmin){
          xxmin=xxdif;
          iimin=j;
        }
      }
      meshi->plotx=iimin;
    }

    ymin = meshi->yplt[0];
    ymax = meshi->yplt[meshi->jbar];
    if(ymin<=yval&&yval<=ymax){
      int j;
      int iimin;
      float yymin,yydif;

      iimin=0;
      yydif = yval - meshi->yplt[0];
      if(yydif<0)yydif=-yydif;
      yymin = yydif;

      for(j=1;j<=meshi->jbar;j++){
        yydif = yval - meshi->yplt[j];
        if(yydif<0)yydif=-yydif;
        if(yydif<yymin){
          yymin=yydif;
          iimin=j;
        }
      }
      meshi->ploty=iimin;
    }

    zmin = meshi->zplt[0];
    zmax = meshi->zplt[meshi->kbar];
    if(zmin<=zval&&zval<=zmax){
      int j;
      int iimin;
      float zzmin,zzdif;

      iimin=0;
      zzdif = zval - meshi->zplt[0];
      if(zzdif<0)zzdif=-zzdif;
      zzmin = zzdif;

      for(j=1;j<=meshi->kbar;j++){
        zzdif = zval - meshi->zplt[j];
        if(zzdif<0)zzdif=-zzdif;
        if(zzdif<zzmin){
          zzmin=zzdif;
          iimin=j;
        }
      }
      meshi->plotz=iimin;
    }



  }
}

/* ------------------ updateplotslice ------------------------ */

void updateplotslice(int slicedir){
  int i;
  if(current_mesh->plotx<0)current_mesh->plotx=current_mesh->ibar;
  if(current_mesh->plotx>current_mesh->ibar)current_mesh->plotx=0;

  if(current_mesh->ploty<0)                 current_mesh->ploty=current_mesh->jbar;
  if(current_mesh->ploty>current_mesh->jbar)current_mesh->ploty=0;

  if(current_mesh->plotz<0)                 current_mesh->plotz=current_mesh->kbar;
  if(current_mesh->plotz>current_mesh->kbar)current_mesh->plotz=0;
  update_plot_xyz(current_mesh);
  for(i=0;i<nmeshes;i++){
    mesh *meshjj;

    meshjj = meshinfo + i;
    if(meshjj->plot3dfilenum==-1)continue;
    updateplotslice_mesh(meshjj,slicedir);
  }
}

/* ------------------ updateplotslice_mesh ------------------------ */

void updateplotslice_mesh(mesh *mesh_in, int slicedir){
  int i, j, k;
  int plotx, ploty, plotz;
  int ibar, jbar, kbar;
  float *xplt, *yplt, *zplt;
  mesh *meshi;
  contour *plot3dcontour1ptr, *plot3dcontour2ptr, *plot3dcontour3ptr;
  int *yzcolorbase, *xzcolorbase, *xycolorbase; 
  float *yzcolorfbase, *xzcolorfbase, *xycolorfbase; 
  float *yzcolortbase, *xzcolortbase, *xycolortbase; 
  float *dx_xy, *dy_xy, *dz_xy;
  float *dx_xz, *dy_xz, *dz_xz;
  float *dx_yz, *dy_yz, *dz_yz;
  float *dx_yzcopy, *dy_yzcopy, *dz_yzcopy;
  float *dx_xzcopy, *dy_xzcopy, *dz_xzcopy;
  float *dx_xycopy, *dy_xycopy, *dz_xycopy;
  float *qdata;
  int *iqdata;
  char *iblank_xy=NULL, *iblank_xz=NULL, *iblank_yz=NULL;
  int nx, ny, nz, nxy, nxyz;
  char *iblank_x, *iblank_y, *iblank_z;

  meshi = mesh_in;
  plotx = meshi->plotx;
  ploty = meshi->ploty;
  plotz = meshi->plotz;
  ibar = meshi->ibar;
  jbar = meshi->jbar;
  kbar = meshi->kbar;
  xplt = meshi->xplt;
  yplt = meshi->yplt;
  zplt = meshi->zplt;
  iblank_x = meshi->c_iblank_x;
  iblank_y = meshi->c_iblank_y;
  iblank_z = meshi->c_iblank_z;

  yzcolorbase=meshi->yzcolorbase;
  xzcolorbase=meshi->xzcolorbase;
  xycolorbase=meshi->xycolorbase;

  yzcolorfbase=meshi->yzcolorfbase;
  xzcolorfbase=meshi->xzcolorfbase;
  xycolorfbase=meshi->xycolorfbase;

  yzcolortbase=meshi->yzcolortbase;
  xzcolortbase=meshi->xzcolortbase;
  xycolortbase=meshi->xycolortbase;

  dx_xy=meshi->dx_xy;
  dx_xz=meshi->dx_xz;
  dx_yz=meshi->dx_yz;

  dy_xy=meshi->dy_xy;
  dy_xz=meshi->dy_xz;
  dy_yz=meshi->dy_yz;

  dz_xy=meshi->dz_xy;
  dz_xz=meshi->dz_xz;
  dz_yz=meshi->dz_yz;

  plot3dcontour1ptr = &meshi->plot3dcontour1;
  plot3dcontour2ptr = &meshi->plot3dcontour2;
  plot3dcontour3ptr = &meshi->plot3dcontour3;
  nx = ibar+1;
  ny = jbar + 1;
  nz = kbar + 1;
  nxy = nx*ny;
  nxyz = nx*ny*nz;


  iqdata = meshi->iqdata;
  qdata = meshi->qdata;

  if(plotx<0){
    plotx=ibar;
  }
  if(plotx>ibar){
    plotx=0;
  }
  if(ploty<0){
    ploty=jbar;
  }
  if(ploty>jbar){
    ploty=0;
  }
  if(plotz<0){
    plotz=kbar;
  }
  if(plotz>kbar){
    plotz=0;
  }
  meshi->plotx=plotx;
  meshi->ploty=ploty;
  meshi->plotz=plotz;

  minfill=1;
  maxfill=1;
  if(setp3min[plotn-1]==CHOP_MIN)minfill=0;
  if(setp3max[plotn-1]==CHOP_MAX)maxfill=0;

  if(ReadPlot3dFile!=1)return;
  if(slicedir==1){
    float *yzcolorf;
    float *yzcolort;
    int *yzcolor;

    yzcolor=yzcolorbase;
    yzcolorf=yzcolorfbase;
    yzcolort=yzcolortbase;
    dx_yzcopy=dx_yz;
    dy_yzcopy=dy_yz;
    dz_yzcopy=dz_yz;
    NewMemory((void **)&iblank_yz,jbar*kbar*sizeof(char));
    for(j=0;j<jbar;j++){for(k=0;k<kbar;k++){
      iblank_yz[k+j*kbar]=iblank_x[ijknode(plotx,j,k)];
    }}
    for(j=0;j<=jbar;j++){for(k=0;k<=kbar;k++){
      *yzcolor++=iqdata[ijkn(plotx,j,k,plotn-1)];
      *yzcolort++=(float)iqdata[ijkn(plotx,j,k,plotn-1)]/255.0;
      *yzcolorf++=qdata[ijkn(plotx,j,k,plotn-1)];
      *dx_yzcopy = 0.0;
      *dy_yzcopy = 0.0;
      *dz_yzcopy = 0.0;
      if(uindex!=-1)*dx_yzcopy++=vecfactor*veclength*qdata[ijkn(plotx,j,k,uindex)]/speedmax;
      if(vindex!=-1)*dy_yzcopy++=vecfactor*veclength*qdata[ijkn(plotx,j,k,vindex)]/speedmax;
      if(windex!=-1)*dz_yzcopy++=vecfactor*veclength*qdata[ijkn(plotx,j,k,windex)]/speedmax;
    }}
    freecontour(plot3dcontour1ptr);
    initcontour(plot3dcontour1ptr,rgbptr,nrgb+1);
    setcontourslice(plot3dcontour1ptr,1,xplt[plotx]);
    getcontours(yplt,zplt,jbar+1,kbar+1, 
      yzcolorfbase, iblank_yz,p3levels[plotn-1],
      minfill, maxfill,
      plot3dcontour1ptr);
    FREEMEMORY(iblank_yz);
  }
  else if(slicedir==2){
    float *xzcolorf;
    float *xzcolort;
    int *xzcolor;

    xzcolor=xzcolorbase;
    xzcolorf=xzcolorfbase;
    xzcolort=xzcolortbase;
    dx_xzcopy=dx_xz;
    dy_xzcopy=dy_xz;
    dz_xzcopy=dz_xz;
    NewMemory((void **)&iblank_xz,ibar*kbar*sizeof(char));
    for(i=0;i<ibar;i++){for(k=0;k<kbar;k++){
      iblank_xz[k+i*kbar]=iblank_y[ijknode(i,ploty,k)];
    }}
    for(i=0;i<=ibar;i++){for(k=0;k<=kbar;k++){
      *xzcolor++=iqdata[ijkn(i,ploty,k,plotn-1)];
      *xzcolort++=(float)iqdata[ijkn(i,ploty,k,plotn-1)]/255.0;
      *xzcolorf++=qdata[ijkn(i,ploty,k,plotn-1)];
      *dx_xzcopy = 0.0;
      *dy_xzcopy = 0.0;
      *dz_xzcopy = 0.0;
      if(uindex!=-1)*dx_xzcopy++=vecfactor*veclength*qdata[ijkn(i,ploty,k,uindex)]/speedmax;
      if(vindex!=-1)*dy_xzcopy++=vecfactor*veclength*qdata[ijkn(i,ploty,k,vindex)]/speedmax;
      if(windex!=-1)*dz_xzcopy++=vecfactor*veclength*qdata[ijkn(i,ploty,k,windex)]/speedmax;
    }}
    freecontour(plot3dcontour2ptr);
    initcontour(plot3dcontour2ptr,rgbptr,nrgb+1);
    setcontourslice(plot3dcontour2ptr,2,yplt[ploty]);
    getcontours(xplt,zplt,ibar+1,kbar+1, 
      xzcolorfbase, iblank_xz,p3levels[plotn-1], 
      minfill, maxfill,
      plot3dcontour2ptr);
    FREEMEMORY(iblank_xz);
  }
  else if(slicedir==3){
    float *xycolorf;
    float *xycolort;
    int *xycolor;

    xycolor=xycolorbase;
    xycolorf=xycolorfbase;
    xycolort=xycolortbase;
    dx_xycopy=dx_xy;
    dy_xycopy=dy_xy;
    dz_xycopy=dz_xy;
    NewMemory((void **)&iblank_xy,ibar*jbar*sizeof(char));
    for(i=0;i<ibar;i++){for(j=0;j<jbar;j++){
      iblank_xy[j+i*jbar]=iblank_z[ijknode(i,j,plotz)];
    }}
    for(i=0;i<=ibar;i++){for(j=0;j<=jbar;j++){
      *xycolor++=iqdata[ijkn(i,j,plotz,plotn-1)];
      *xycolort++=(float)iqdata[ijkn(i,j,plotz,plotn-1)]/255.0;
      *xycolorf++=qdata[ijkn(i,j,plotz,plotn-1)];
      *dx_xycopy = 0.0;
      *dy_xycopy = 0.0;
      *dz_xycopy = 0.0;
      if(uindex!=-1)*dx_xycopy++=vecfactor*veclength*qdata[ijkn(i,j,plotz,uindex)]/speedmax;
      if(vindex!=-1)*dy_xycopy++=vecfactor*veclength*qdata[ijkn(i,j,plotz,vindex)]/speedmax;
      if(windex!=-1)*dz_xycopy++=vecfactor*veclength*qdata[ijkn(i,j,plotz,windex)]/speedmax;
    }}
    freecontour(plot3dcontour3ptr);
    initcontour(plot3dcontour3ptr,rgbptr,nrgb+1);
    setcontourslice(plot3dcontour3ptr,3,zplt[plotz]);
    getcontours(xplt,yplt,ibar+1,jbar+1, 
      xycolorfbase, iblank_xy,p3levels[plotn-1], 
      minfill, maxfill,
      plot3dcontour3ptr);
    FREEMEMORY(iblank_xy);
  }
}

/* ------------------ updateshowstep ------------------------ */

void updateshowstep(int val, int slicedir){

  if(ReadPlot3dFile==0&&visGrid==0&&ReadVolSlice==0)return;
  current_mesh->slicedir=slicedir;
  switch (slicedir){
  case 1:
    if(current_mesh->visx!=val)updatemenu=1;
    current_mesh->visx = val;
    break;
  case 2:
    if(current_mesh->visy!=val)updatemenu=1;
    current_mesh->visy = val;
    break;
  case 3:
    if(current_mesh->visz!=val)updatemenu=1;
    current_mesh->visz = val;
    break;
  case 4:
    if(ReadPlot3dFile==1){
      if(visiso!=val)updatemenu=1;
      visiso = val;
    }
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  plotstate=getplotstate(STATIC_PLOTS);
  stept=0;
  if(ReadVolSlice==0&&plotstate==DYNAMIC_PLOTS&&visGrid==0)updatetimes();
  {
    int i;
    float xmin, xmax;
    float ymin, ymax;
    float zmin, zmax;
#define MESHEPS 0.0001

    xmin = current_mesh->xplt[0];
    xmax = current_mesh->xplt[current_mesh->ibar];
    ymin = current_mesh->yplt[0];
    ymax = current_mesh->yplt[current_mesh->jbar];
    zmin = current_mesh->zplt[0];
    zmax = current_mesh->zplt[current_mesh->kbar];
    for(i=0;i<nmeshes;i++){
      mesh *meshi;
      float xmin2, xmax2;
      float ymin2, ymax2;
      float zmin2, zmax2;

      meshi = meshinfo+i;
      if(meshi==current_mesh)continue;

      xmin2 = meshi->xplt[0];
      xmax2 = meshi->xplt[meshi->ibar];
      ymin2 = meshi->yplt[0];
      ymax2 = meshi->yplt[meshi->jbar];
      zmin2 = meshi->zplt[0];
      zmax2 = meshi->zplt[meshi->kbar];


      if(slicedir==1&&(xmax-MESHEPS<xmin2||xmax2-MESHEPS<xmin))continue;
      if(slicedir==2&&(ymax-MESHEPS<ymin2||ymax2-MESHEPS<ymin))continue;
      if(slicedir==3&&(zmax-MESHEPS<zmin2||zmax2-MESHEPS<zmin))continue;

      meshi->visx=current_mesh->visx;
      meshi->visx=meshi->visx2;
      meshi->visy=current_mesh->visy;
      meshi->visz=current_mesh->visz;
    }
    {
      int nplanes=0;

      if(current_mesh->visx==1)nplanes++;
      if(current_mesh->visy==1)nplanes++;
      if(current_mesh->visz==1)nplanes++;
      if(nplanes<=1){
        transparentflagVOL=1;
      }
      else{
        transparentflagVOL=0;
      }
    }
  }
}

/* ------------------ drawgrid ------------------------ */

void drawgrid(const mesh *meshi){
  int i, j, k;
  float *xplt, *yplt, *zplt;
  int ibar, jbar, kbar;
  int plotx, ploty, plotz;

  xplt = meshi->xplt;
  yplt = meshi->yplt;
  zplt = meshi->zplt;
  ibar = meshi->ibar;
  jbar = meshi->jbar;
  kbar = meshi->kbar;
  plotx = meshi->plotx;
  ploty = meshi->ploty;
  plotz = meshi->plotz;


  if(meshi->visx==1||meshi->visy==1||meshi->visz==1){
    antialias(1);
    glLineWidth(vectorlinewidth);
    if(meshi->meshrgb_ptr!=NULL){
      const float *rgb;

      rgb=meshi->meshrgb;
      glColor3f(rgb[0],rgb[1],rgb[2]);
    }
    else{
      glColor4fv(foregroundcolor);
    }
    glBegin(GL_LINES);

    if(meshi->visx==1){
      for(j=0;j<jbar+1;j++){
        glVertex3f(xplt[plotx],yplt[j],zplt[0]);
        glVertex3f(xplt[plotx],yplt[j],zplt[kbar]);
      }
      for(k=0;k<kbar+1;k++){
        glVertex3f(xplt[plotx],yplt[0],zplt[k]);
        glVertex3f(xplt[plotx],yplt[jbar],zplt[k]);
      }
      /*
      if(plotx>0){
        for(j=0;j<jbar+1;j++){
          glVertex3f(xplt[plotx-1],yplt[j],zplt[0]);
          glVertex3f(xplt[plotx-1],yplt[j],zplt[kbar]);
        }
        for(k=0;k<kbar+1;k++){
          glVertex3f(xplt[plotx-1],yplt[0],zplt[k]);
          glVertex3f(xplt[plotx-1],yplt[jbar],zplt[k]);
        }
        for(j=0;j<jbar+1;j++){
          for(k=0;k<kbar+1;k++){
            glVertex3f(xplt[plotx-1],yplt[j],zplt[k]);
            glVertex3f(xplt[plotx],yplt[j],zplt[k]);
          }
        }
      }
      */
    }

    if(meshi->visy==1){
      for(i=0;i<ibar+1;i++){
        glVertex3f(xplt[i],yplt[ploty],zplt[0]);
        glVertex3f(xplt[i],yplt[ploty],zplt[kbar]);
      }
      for(k=0;k<kbar+1;k++){
        glVertex3f(   xplt[0],yplt[ploty],zplt[k]);
        glVertex3f(xplt[ibar],yplt[ploty],zplt[k]);
      }
    }
    if(meshi->visz==1){
      for(i=0;i<ibar+1;i++){
        glVertex3f(xplt[i],   yplt[0],zplt[plotz]);
        glVertex3f(xplt[i],yplt[jbar],zplt[plotz]);
      }
      for(j=0;j<jbar+1;j++){
        glVertex3f(xplt[0],yplt[j],zplt[plotz]);
        glVertex3f(xplt[ibar],yplt[j],zplt[plotz]);
      }
    }

    glEnd();
    antialias(0);
  }
}


/* ------------------ updateplot3dmenulabel ------------------------ */

void updateplot3dmenulabels(void){
  int i;
  plot3d *plot3di;
  char label[128];

  if(nplot3d>0){
    FREEMEMORY(plot3dorderindex);
    NewMemory((void **)&plot3dorderindex,sizeof(int)*nplot3d);
    for(i=0;i<nplot3d;i++){
      plot3dorderindex[i]=i;
    }
    qsort( (int *)plot3dorderindex, (size_t)nplot3d, sizeof(int), plot3dcompare );

    for(i=0;i<nplot3d;i++){
      plot3di = plot3dinfo + i;
      STRCPY(plot3di->menulabel,plot3di->longlabel);
      STRCPY(plot3di->menulabel,"");
      if(plot3di->time>=0.0){
        sprintf(label,"%f",plot3di->time);
        trimzeros(label);
        STRCAT(label," s");
        STRCAT(plot3di->menulabel,label);
      }
      if(nmeshes>1){
        if(plot3di->time>=0.0)STRCAT(plot3di->menulabel,", ");
        sprintf(label,"Mesh %i",1+plot3di->blocknumber);
        STRCAT(plot3di->menulabel,label);
      }
      if(showfiles==1||plot3di->time<0.0){
        if(plot3di->time>=0.0||nmeshes>1)STRCAT(plot3di->menulabel,", ");
        STRCAT(plot3di->menulabel,plot3di->file);
      }
    } 
  }


}

/* ------------------ plot3dlistcompare  ------------------------ */

int plot3dlistcompare( const void *arg1, const void *arg2 ){
  float val1, val2;

  val1 = *(float *)arg1;
  val2 = *(float *)arg2;

  if(val1<val2)return -1;
  if(val1>val2)return 1;
  return 0;
}
/* ------------------ init_plot3dtimelist  ------------------------ */

void init_plot3dtimelist(void){
  int i;
  plot3d *plot3di;
  float lasttime,val;

  nplot3dtimelist=0;
  if(plot3dtimelist==NULL&&nplot3d>0){
    NewMemory((void **)&plot3dtimelist,nplot3d*sizeof(float));
  }
  if(plot3dtimelist==NULL)return;

  for(i=0;i<nplot3d;i++){
    plot3di = plot3dinfo + i;
    plot3dtimelist[i]=plot3di->time;
  }
  qsort( (float *)plot3dtimelist, (size_t)nplot3d, sizeof(int), plot3dlistcompare );
  lasttime=-999999.0;
  for(i=0;i<nplot3d;i++){
    val=plot3dtimelist[i];
    if(fabs((double)(val-lasttime))>0.1){
      nplot3dtimelist++;
      plot3dtimelist[nplot3dtimelist-1]=val;
      lasttime=val;
    }
  }
}


