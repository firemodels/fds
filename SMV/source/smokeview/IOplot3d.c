#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include GLUT_H

#include "update.h"
#include "smokeviewvars.h"

void update_glui_plot3d(void);

/* ------------------ plot3dcompare  ------------------------ */

int plot3dcompare( const void *arg1, const void *arg2 ){
  plot3ddata *plot3di, *plot3dj;

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

/* ------------------ readplot3d  ------------------------ */

void readplot3d(char *file, int ifile, int flag, int *errorcode){
  int n, nn, ntotal, i, nnn;
  float *udata, *vdata, *wdata, *sdata;
  float sum;
  int error;
  int ibar,jbar,kbar;
  meshdata *meshi,*gbb,*gbi;
  plot3ddata *p;
  int nloaded=0;
  int nx, ny, nz;
  int pn;
  FILE_SIZE plot3dfilelen;
  int local_starttime=0, local_stoptime=0;
  FILE_SIZE file_size=0;
  int local_starttime0=0, local_stoptime0=0;
  float delta_time, delta_time0;

  CheckMemory;
  local_starttime0 = glutGet(GLUT_ELAPSED_TIME);
  *errorcode=0;

  ASSERT(ifile>=0&&ifile<nplot3dinfo);
  p=plot3dinfo+ifile;
  if(flag==UNLOAD&&p->loaded==0)return;

  highlight_mesh=p->blocknumber;
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
  initcontour(&meshi->plot3dcontour1,rgb_plot3d_contour,nrgb);
  initcontour(&meshi->plot3dcontour2,rgb_plot3d_contour,nrgb);
  initcontour(&meshi->plot3dcontour3,rgb_plot3d_contour,nrgb);


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
    if(colorlabeliso != NULL){
      for(nn=0;nn<mxplot3dvars;nn++){
        for(nnn=0;nnn<MAXRGB;nnn++){
          FREEMEMORY((*(colorlabeliso+nn))[nnn]);
        }
        FREEMEMORY(colorlabeliso[nn]);
      }
      FREEMEMORY(colorlabeliso);
    }
    if(scalep3 != NULL){
      for(nn=0;nn<mxplot3dvars;nn++){
        FREEMEMORY(scalep3[nn]);
      }
      FREEMEMORY(scalep3);
    }
    if(fscalep3!=NULL){
      FREEMEMORY(fscalep3);
    }
    if(p3levels!=NULL){
      for(nn=0;nn<mxplot3dvars;nn++){
        FREEMEMORY(p3levels[nn]);
      }
      FREEMEMORY(p3levels);
    }
    if(p3levels256!=NULL){
      for(nn=0;nn<mxplot3dvars;nn++){
        FREEMEMORY(p3levels256[nn]);
      }
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
#ifdef pp_MEMPRINT
    PRINTF("After plot3d unload: \n");
    PrintMemoryInfo;
#endif
    Update_Times();
    update_unit_defs();
    update_glui_plot3d();
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

  if(p->compression_type==UNCOMPRESSED){
    if(NewMemory((void **)&meshi->qdata,numplot3dvars*ntotal*sizeof(float))==0){
      *errorcode=1;
      readplot3d("",ifile,UNLOAD,&error);
      return;
    }
  }

  if(NewMemory((void **)&meshi->yzcolorbase ,ny*nz*sizeof(unsigned char))==0||
     NewMemory((void **)&meshi->xzcolorbase ,nx*nz*sizeof(unsigned char))==0||
     NewMemory((void **)&meshi->xycolorbase ,nx*ny*sizeof(unsigned char))==0||
     NewMemory((void **)&meshi->yzcolorfbase,ny*nz*sizeof(float))==0||
     NewMemory((void **)&meshi->xzcolorfbase,nx*nz*sizeof(float))==0||
     NewMemory((void **)&meshi->xycolorfbase,nx*ny*sizeof(float))==0||
     NewMemory((void **)&meshi->yzcolortbase,ny*nz*sizeof(float))==0||
     NewMemory((void **)&meshi->xzcolortbase,nx*nz*sizeof(float))==0||
     NewMemory((void **)&meshi->xycolortbase,nx*ny*sizeof(float))==0||
     NewMemory((void **)&meshi->dx_xy       ,nx*ny*sizeof(float))==0||
     NewMemory((void **)&meshi->dy_xy       ,nx*ny*sizeof(float))==0||
     NewMemory((void **)&meshi->dz_xy       ,nx*ny*sizeof(float))==0||
     NewMemory((void **)&meshi->dx_xz       ,nx*nz*sizeof(float))==0||
     NewMemory((void **)&meshi->dy_xz       ,nx*nz*sizeof(float))==0||
     NewMemory((void **)&meshi->dz_xz       ,nx*nz*sizeof(float))==0||
     NewMemory((void **)&meshi->dx_yz       ,ny*nz*sizeof(float))==0||
     NewMemory((void **)&meshi->dy_yz       ,ny*nz*sizeof(float))==0||
     NewMemory((void **)&meshi->dz_yz       ,ny*nz*sizeof(float))==0){
     *errorcode=1;
     readplot3d("",ifile,UNLOAD,&error);
     return;
  }

  file_size=get_filesize(file);
  plot3dfilelen = strlen(file);
  PRINTF("Loading plot3d data: %s\n",file);
  local_starttime = glutGet(GLUT_ELAPSED_TIME);
  if(p->compression_type==UNCOMPRESSED){
    FORTgetplot3dq(file,&nx,&ny,&nz,meshi->qdata,&error,&isotest,plot3dfilelen);
  }
  if(NewMemory((void **)&meshi->iqdata,numplot3dvars*ntotal*sizeof(unsigned char))==0){
    *errorcode=1;
    readplot3d("",ifile,UNLOAD,&error);
    return;
  }
  if(p->compression_type==COMPRESSED_ZLIB){
  }
  local_stoptime = glutGet(GLUT_ELAPSED_TIME);
  delta_time = (local_stoptime-local_starttime)/1000.0;
  p->loaded=1;
  p->display=1;
  speedmax = -1.;
  meshi->udata=NULL;
  meshi->vdata=NULL;
  meshi->wdata=NULL;
  if(p->compression_type==UNCOMPRESSED){
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
    if(uindex!=-1)meshi->udata=meshi->qdata + ntotal*uindex;
    if(vindex!=-1)meshi->vdata=meshi->qdata + ntotal*vindex;
    if(windex!=-1)meshi->wdata=meshi->qdata + ntotal*windex;
  }

  if(NewMemory((void **)&colorlabelp3,mxplot3dvars*sizeof(char **))==0||
     NewMemory((void **)&colorlabeliso,mxplot3dvars*sizeof(char **))==0||
     NewMemory((void **)&fscalep3     ,mxplot3dvars*sizeof(float))==0||
     NewMemory((void **)&scalep3     ,mxplot3dvars*sizeof(char *))==0){
    *errorcode=1;
    readplot3d("",ifile,UNLOAD,&error);
    return;
  }
  for(nn=0;nn<mxplot3dvars;nn++){
    colorlabelp3[nn]=NULL;
    colorlabeliso[nn]=NULL;
    scalep3[nn]=NULL;
    fscalep3[nn]=1.0;
  }
  for(nn=0;nn<mxplot3dvars;nn++){
    if(NewMemory((void **)&colorlabelp3[nn],MAXRGB*sizeof(char *))==0||
       NewMemory((void **)&colorlabeliso[nn],MAXRGB*sizeof(char *))==0||
       NewMemory((void **)&scalep3[nn],30*sizeof(char))==0){
      *errorcode=1;
      readplot3d("",ifile,UNLOAD,&error);
      return;
    }
  }

  scalep3copy = scalep3;
  if(NewMemory((void **)&p3levels,mxplot3dvars*sizeof(float *))==0||
     NewMemory((void **)&p3levels256,mxplot3dvars*sizeof(float *))==0){
    *errorcode=1;
    readplot3d("",ifile,UNLOAD,&error);
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
    if(nplot3dinfo>0){
      shortp3label[nn] = plot3dinfo[ifile].label[nn].shortlabel;
      unitp3label[nn] = plot3dinfo[ifile].label[nn].unit;
    }
    else{
      char numstring[4];

      sprintf(numstring,"%i",nn);
      strcpy(shortp3label[nn],numstring);
      unitp3label[nn] = blank_global;
    }

    for(n=0;n<MAXRGB;n++){
      (*(colorlabelp3+nn))[n]=NULL;
      (*(colorlabeliso+nn))[n]=NULL;
    }

    if(NewMemory((void **)&p3levels[nn],(nrgb+1)*sizeof(float))==0||
       NewMemory((void **)&p3levels256[nn],256*sizeof(float))==0){
      readplot3d("",ifile,UNLOAD,&error);
      *errorcode=1;
      return;
    }
    for(n=0;n<nrgb;n++){
      if(NewMemory((void **)&(*(colorlabelp3+nn))[n],11)==0){
        *errorcode=1;
        readplot3d("",ifile,UNLOAD,&error);
        return;
      }
      if(NewMemory((void **)&(*(colorlabeliso+nn))[n],11)==0){
        *errorcode=1;
        readplot3d("",ifile,UNLOAD,&error);
        return;
      }
    }
    getPlot3DColors(nn,
                  setp3min[nn],p3min+nn, setp3max[nn],p3max+nn,
                  nrgb_full, nrgb-1, *(colorlabelp3+nn),*(colorlabeliso+nn),scalep3copy,fscalep3+nn,p3levels[nn],p3levels256[nn],
                  plot3dinfo[ifile].extreme_min+nn,plot3dinfo[ifile].extreme_max+nn);
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
  updateplotslice(XDIR);
  updateplotslice(YDIR);
  updateplotslice(ZDIR);
  visGrid=0;
  meshi->visInteriorPatches=0;
  if(visx_all==1){
    updateshowstep(1,XDIR);
  }
  if(visy_all==1){
    updateshowstep(1,YDIR);
  }
  if(visz_all==1){
    updateshowstep(1,ZDIR);
  }
  if(visiso==1){
    updatesurface();
  }

  updateplot3dlistindex();
#ifdef pp_MEMPRINT
  PRINTF("After plot3d load: \n");
  PrintMemoryInfo;
#endif
  Update_Times();
  update_unit_defs();
  Idle_CB();
  local_stoptime0 = glutGet(GLUT_ELAPSED_TIME);
  delta_time0=(local_stoptime0-local_starttime0)/1000.0;

  if(file_size!=0&&delta_time>0.0){
    float loadrate;

    loadrate = ((float)file_size*8.0/1000000.0)/delta_time;
    PRINTF(" %.1f MB loaded in %.2f s - rate: %.1f Mb/s (overhead: %.2f s)\n",
    (float)file_size/1000000.,delta_time,loadrate,delta_time0-delta_time);
  }
  else{
    PRINTF(" %.1f MB downloaded in %.2f s (overhead: %.2f s)",
    (float)file_size/1000000.,delta_time,delta_time0-delta_time);
  }
  if(p->compression_type==COMPRESSED_ZLIB||cache_qdata==0){
    cache_qdata=0;
    FREEMEMORY(meshi->qdata);
  }
  update_glui_plot3d();
  update_plot3dtitle();
  glutPostRedisplay();
}

/* ------------------ update_plot3dtitle ------------------------ */

void update_plot3dtitle(void){
  int filenum;
  plot3ddata *plot3di;
  meshdata *meshi;
  char title_base[1024];

  getBaseTitle("Smokeview ", title_base);
  STRCPY(plot3d_title,title_base);
  meshi=current_mesh;
  if(meshi==NULL)meshi=meshinfo;
  filenum=meshi->plot3dfilenum;
  if(filenum!=-1){
    plot3di = plot3dinfo+meshi->plot3dfilenum;
    STRCAT(plot3d_title,", ");
    STRCAT(plot3d_title,plot3di->file);
  }
}
/* ------------------ drawplot3d_texture ------------------------ */

void drawplot3d_texture(meshdata *meshi){
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
  unsigned char *yzcolorbase, *xzcolorbase, *xycolorbase;
  float *yzcolortbase, *xzcolortbase, *xycolortbase;

  float *dx_xy, *dy_xy, *dz_xy;
  float *dx_xz, *dy_xz, *dz_xz;
  float *dx_yz, *dy_yz, *dz_yz;
  float *dx_yzcopy, *dy_yzcopy, *dz_yzcopy;
  float *dx_xzcopy, *dy_xzcopy, *dz_xzcopy;
  float *dx_xycopy, *dy_xycopy, *dz_xycopy;
  int nx, ny, nz,nxy;
  char *c_iblank_x, *c_iblank_y, *c_iblank_z, *iblank;
  float *vector_color;

  plotx = meshi->iplotx_all[iplotx_all];
  ploty = meshi->iploty_all[iploty_all];
  plotz = meshi->iplotz_all[iplotz_all];

  visx = visx_all;
  visy = visy_all;
  visz = visz_all;
  ibar = meshi->ibar;
  jbar = meshi->jbar;
  kbar = meshi->kbar;
  xplt = meshi->xplt;
  yplt = meshi->yplt;
  zplt = meshi->zplt;
  c_iblank_x = meshi->c_iblank_x;
  c_iblank_y = meshi->c_iblank_y;
  c_iblank_z = meshi->c_iblank_z;
  iblank = meshi->c_iblank_node;


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
    drawstaticiso(currentsurfptr,p3dsurfacetype,p3dsurfacesmooth,2,0,plot3dlinewidth);
    if(surfincrement!=0)drawstaticiso(currentsurf2ptr,p3dsurfacetype,p3dsurfacesmooth,2,0,plot3dlinewidth);
    if(visGrid!=noGridnoProbe){
      if(use_transparency_data==1)transparentoff();
      if(cullfaces==1)glEnable(GL_CULL_FACE);
      return;
    }
  }

  if(use_transparency_data==1){
    transparenton();
  }
  if(visVector==0&&contour_type!=STEPPED_CONTOURS){
    glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
    glEnable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D,texture_plot3d_colorbar_id);
  }

  /* +++++++++++++++++++++++++++   draw yz contours +++++++++++++++++++++++++++++++++++++ */

  if(plotx>=0&&visx!=0){
    if(visVector==0&&contour_type==STEPPED_CONTOURS){
      DrawContours(plot3dcontour1ptr);
    }
    if(visVector==0&&contour_type!=STEPPED_CONTOURS){
      if(plotx<0){
        plotx=ibar;
        updateplotslice(XDIR);
      }
      if(plotx>ibar){
        plotx=0;
        updateplotslice(XDIR);
      }
      glBegin(GL_TRIANGLES);
      for(j=0; j<jbar; j++){
        color1t=yzcolortbase + j*nz;
        color2t=color1t+nz;
        for(k=0; k<kbar; k++){
          if(c_iblank_x==NULL||c_iblank_x[IJKNODE(plotx,j,k)]==GASGAS){
            if(ABS(color1t[k]-color2t[k+1])<ABS(color1t[k+1]-color2t[k])){
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
      unsigned char *yzcolor;

      yzcolor=yzcolorbase;
      dx_yzcopy=dx_yz; dy_yzcopy=dy_yz; dz_yzcopy=dz_yz;
      antialias(ON);
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
          if((iblank==NULL||iblank[IJKNODE(plotx,j,k)]==GAS)&&vector_color[3]>0.5){
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
        }
      }
      glEnd();
      antialias(OFF);

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
          if((iblank==NULL||iblank[IJKNODE(plotx,j,k)]==GAS)&&vector_color[3]>0.5){
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
        }
      }
      glEnd();
    }
  }

  /* +++++++++++++++++++++++++++++++++  draw xz contours  ++++++++++++++++++++++++++++++++++++++++ */

  if(ploty>=0&&visy!=0){
    if(visVector==0&&contour_type==STEPPED_CONTOURS){
      DrawContours(plot3dcontour2ptr);
    }
    if(visVector==0&&contour_type!=STEPPED_CONTOURS){
      glBegin(GL_TRIANGLES);
      for(i=0; i<ibar; i++){
        color1t=xzcolortbase + i*nz;
        color2t=color1t+nz;
        for(k=0; k<kbar; k++){
          if(c_iblank_y==NULL||c_iblank_y[IJKNODE(i,ploty,k)]==GASGAS){
            if(ABS(color1t[k]-color2t[k+1])<ABS(color1t[k+1]-color2t[k])){
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
      unsigned char*xzcolor;

      xzcolor=xzcolorbase;
      antialias(ON);
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
          if((iblank==NULL||iblank[IJKNODE(i,ploty,k)]==GAS)&&vector_color[3]>0.5){
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
        }
      }
      glEnd();
      antialias(OFF);

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
          if((iblank==NULL||iblank[IJKNODE(i,ploty,k)]==GAS)&&vector_color[3]>0.5){
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

  if(plotz>=0&&visz!=0){
    if(visVector==0&&contour_type==STEPPED_CONTOURS){
      DrawContours(plot3dcontour3ptr);
    }
    if(visVector==0&&contour_type!=STEPPED_CONTOURS){
      if(plotz<0){
        plotz=kbar;
        updateplotslice(ZDIR);
      }
      if(plotz>kbar){
        plotz=0;
        updateplotslice(ZDIR);
      }
      glBegin(GL_TRIANGLES);
      for(i=0; i<ibar; i++){
        color1t=xycolortbase + i*ny;
        color2t=color1t+ny;
        for(j=0; j<jbar; j++){
          if(c_iblank_z==NULL||c_iblank_z[IJKNODE(i,j,plotz)]==GASGAS){
            if(ABS(color1t[j]-color2t[j+1])<ABS(color1t[j+1]-color2t[j])){
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
      unsigned char *xycolor;

      xycolor=xycolorbase;
      antialias(ON);
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
          if((iblank==NULL||iblank[IJKNODE(i,j,plotz)]==GAS)&&vector_color[3]>0.5){
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
        }
      }
      glEnd();
      antialias(OFF);

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
          if((iblank==NULL||iblank[IJKNODE(i,j,plotz)]==GAS)&&vector_color[3]>0.5){
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
        }
      }
      glEnd();
    }
  }
  if(visVector==0&&contour_type!=STEPPED_CONTOURS){
    glDisable(GL_TEXTURE_1D);
  }
  if(use_transparency_data==1){
    transparentoff();
  }
  if(cullfaces==1)glEnable(GL_CULL_FACE);

}

/* ------------------ draw_plot3dframe ------------------------ */

void draw_plot3dframe(void){
  int i;

  for(i=0;i<nmeshes;i++){
    meshdata *meshi;

    meshi=meshinfo+i;
    if(meshi->plot3dfilenum==-1)continue;
    if(plot3dinfo[meshi->plot3dfilenum].display==0)continue;
    if(usetexturebar!=0){
      drawplot3d_texture(meshi);
    }
    else{
      drawplot3d(meshi);
    }
  }
}

/* ------------------ drawplot3d ------------------------ */

void drawplot3d(meshdata *meshi){
  int i,j,k;
  int colorindex;
  unsigned char *color1, *color2;
  float dx, dy, dz;
  int plotx, ploty, plotz;
  int visx, visy, visz;
  float *xplt, *yplt, *zplt;
  int ibar, jbar, kbar;
  isosurface *currentsurfptr,*currentsurf2ptr;
  contour *plot3dcontour1ptr, *plot3dcontour2ptr, *plot3dcontour3ptr;
  unsigned char *yzcolorbase, *xzcolorbase, *xycolorbase;
  float *dx_xy, *dy_xy, *dz_xy;
  float *dx_xz, *dy_xz, *dz_xz;
  float *dx_yz, *dy_yz, *dz_yz;
  float *dx_yzcopy, *dy_yzcopy, *dz_yzcopy;
  float *dx_xzcopy, *dy_xzcopy, *dz_xzcopy;
  float *dx_xycopy, *dy_xycopy, *dz_xycopy;
  int nx, ny, nz,nxy;
  char *c_iblank_x, *c_iblank_y, *c_iblank_z, *iblank;
  float *vector_color;

  plotx = meshi->iplotx_all[iplotx_all];
  ploty = meshi->iploty_all[iploty_all];
  plotz = meshi->iplotz_all[iplotz_all];


  visx = visx_all;
  visy = visy_all;
  visz = visz_all;
  ibar = meshi->ibar;
  jbar = meshi->jbar;
  kbar = meshi->kbar;
  xplt = meshi->xplt;
  yplt = meshi->yplt;
  zplt = meshi->zplt;
  c_iblank_x = meshi->c_iblank_x;
  c_iblank_y = meshi->c_iblank_y;
  c_iblank_z = meshi->c_iblank_z;
  iblank = meshi->c_iblank_node;


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
    drawstaticiso(currentsurfptr,p3dsurfacetype,p3dsurfacesmooth,2,0,plot3dlinewidth);
    if(surfincrement!=0)drawstaticiso(currentsurf2ptr,p3dsurfacetype,p3dsurfacesmooth,2,0,plot3dlinewidth);
    if(visGrid!=noGridnoProbe){
      if(use_transparency_data==1)transparentoff();
      if(cullfaces==1)glEnable(GL_CULL_FACE);
      return;
    }
  }

  if(use_transparency_data==1)transparenton();

  /* +++++++++++++++++++++++++++   draw yz contours +++++++++++++++++++++++++++++++++++++ */

  if(plotx>=0&&visx!=0){
    if(visVector==0&&contour_type==STEPPED_CONTOURS){
      DrawContours(plot3dcontour1ptr);
    }
    if(visVector==0&&contour_type!=STEPPED_CONTOURS){
      if(plotx<0){
        plotx=ibar;
        updateplotslice(XDIR);
      }
      if(plotx>ibar){
        plotx=0;
        updateplotslice(XDIR);
      }
      glBegin(GL_TRIANGLES);
      for(j=0; j<jbar; j++){
        color1=yzcolorbase + j*nz;
        color2=color1+nz;
        for(k=0; k<kbar; k++){
          if(c_iblank_x==NULL||c_iblank_x[IJKNODE(plotx,j,k)]==GASGAS){
            if(ABS(color1[k]-color2[k+1])<ABS(color1[k+1]-color2[k])){
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
      unsigned char *yzcolor;

      yzcolor=yzcolorbase;
      dx_yzcopy=dx_yz; dy_yzcopy=dy_yz; dz_yzcopy=dz_yz;
      antialias(ON);
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
          if((iblank==NULL||iblank[IJKNODE(plotx,j,k)]==GAS)&&vector_color[3]>0.5){
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
        }
      }
      glEnd();
      antialias(OFF);

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
          if((iblank==NULL||iblank[IJKNODE(plotx,j,k)]==GAS)&&vector_color[3]>0.5){
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
        }
      }
      glEnd();
    }
  }

  /* +++++++++++++++++++++++++++++++++  draw xz contours  ++++++++++++++++++++++++++++++++++++++++ */

  if(ploty>=0&&visy!=0){
    if(visVector==0&&contour_type==STEPPED_CONTOURS){
      DrawContours(plot3dcontour2ptr);
    }
    if(visVector==0&&contour_type!=STEPPED_CONTOURS){
      glBegin(GL_TRIANGLES);
      for(i=0; i<ibar; i++){
        color1=xzcolorbase + i*nz;
        color2=color1+nz;
        for(k=0; k<kbar; k++){
          if(c_iblank_y==NULL||c_iblank_y[IJKNODE(i,ploty,k)]==GASGAS){
            if(ABS(color1[k]-color2[k+1])<ABS(color1[k+1]-color2[k])){
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
      unsigned char *xzcolor;

      xzcolor=xzcolorbase;
      antialias(ON);
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
          if((iblank==NULL||iblank[IJKNODE(i,ploty,k)]==GAS)&&vector_color[3]>0.5){
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
        }
      }
      glEnd();
      antialias(OFF);

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
          if((iblank==NULL||iblank[IJKNODE(i,ploty,k)]==GAS)&&vector_color[3]>0.5){
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

  if(plotz>=0&&visz!=0){
    if(visVector==0&&contour_type==STEPPED_CONTOURS){
      DrawContours(plot3dcontour3ptr);
    }
    if(visVector==0&&contour_type!=STEPPED_CONTOURS){
      if(plotz<0){
        plotz=kbar;
        updateplotslice(ZDIR);
      }
      if(plotz>kbar){
        plotz=0;
        updateplotslice(ZDIR);
      }
      glBegin(GL_TRIANGLES);
      for(i=0; i<ibar; i++){
        color1=xycolorbase + i*ny;
        color2=color1+ny;
        for(j=0; j<jbar; j++){
          if(c_iblank_z==NULL||c_iblank_z[IJKNODE(i,j,plotz)]==GASGAS){
            if(ABS(color1[j]-color2[j+1])<ABS(color1[j+1]-color2[j])){
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
      unsigned char *xycolor;

      xycolor=xycolorbase;
      antialias(ON);
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
          if(iblank[IJKNODE(i,j,plotz)]==GAS&&vector_color[3]>0.5){
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
        }
      }
      glEnd();
      antialias(OFF);

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
          if((iblank==NULL||iblank[IJKNODE(i,j,plotz)]==GAS)&&vector_color[3]>0.5){
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
        }
      }
      glEnd();
    }
  }
  if(use_transparency_data==1)transparentoff();
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
  int plot3dsize;
  int i;

  if(cache_qdata==0)return;
  for(i=0;i<nmeshes;i++){
    float dlevel=-1.0;
    meshdata *meshi;

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
    if(plotiso[plotn-1]<0){
      plotiso[plotn-1]=nrgb-3;
    }
    if(plotiso[plotn-1]>nrgb-3){
      plotiso[plotn-1]=0;
    }
    colorindex=plotiso[plotn-1];
    level = p3min[plotn-1] + (colorindex+0.5)*(p3max[plotn-1]-p3min[plotn-1])/((float)nrgb-2.0f);
    isolevelindex=colorindex;
    isolevelindex2=colorindex;
    freesurface(currentsurfptr);
    InitIsosurface(currentsurfptr, level, rgb_plot3d_contour[colorindex],-999);
    GetIsosurface(currentsurfptr,qdata+(plotn-1)*plot3dsize,NULL,iblank_cell,level,dlevel,
      xplt,ibar+1,yplt,jbar+1,zplt,kbar+1);
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
      InitIsosurface(currentsurf2ptr, level2, rgb_plot3d_contour[colorindex2],-999);
      GetIsosurface(currentsurf2ptr,qdata+(plotn-1)*plot3dsize,NULL,iblank_cell,level2,dlevel,
        xplt,ibar+1,yplt,jbar+1,zplt,kbar+1);
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
  if(visx_all==1)updateplotslice(XDIR);
  if(visy_all==1)updateplotslice(YDIR);
  if(visz_all==1)updateplotslice(ZDIR);
}

/* ------------------ get_plot3d_index ------------------------ */

int get_plot3d_index(meshdata *meshi, int dir, float val){
  float valmin;
  int i, ivalmin, nvals;
  float *xyz;

  switch(dir){
    case XDIR:
      xyz = meshi->xplt_orig;
      nvals = meshi->ibar;
      break;
    case YDIR:
      xyz = meshi->yplt_orig;
      nvals = meshi->jbar;
      break;
    case ZDIR:
      xyz = meshi->zplt_orig;
      nvals = meshi->kbar;
      break;
    default:
      ASSERT(FFALSE);
      break;
  }

  ivalmin=0;
  valmin = ABS(xyz[0]- val);
  for(i=1;i<=nvals;i++){
    if(ABS(xyz[i]-val)<valmin){
      valmin = ABS(xyz[i]-val);
      ivalmin = i;
    }
  }
  return ivalmin;
}

/* ------------------ update_plot_xyz ------------------------ */

void update_plot_xyz(meshdata *current_mesh_local){
  int i;
  float xval, yval, zval;

  xval = current_mesh_local->xplt[current_mesh_local->plotx];
  yval = current_mesh_local->yplt[current_mesh_local->ploty];
  zval = current_mesh_local->zplt[current_mesh_local->plotz];

  for(i=0;i<nmeshes;i++){
    meshdata *meshi;
    float xmin, xmax;
    float ymin, ymax;
    float zmin, zmax;

    meshi=meshinfo+i;
    if(meshi==current_mesh_local)continue;

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

  update_plot_xyz(current_mesh);
  for(i=0;i<nmeshes;i++){
    meshdata *meshjj;

    meshjj = meshinfo + i;
    if(meshjj->plot3dfilenum==-1)continue;
    updateplotslice_mesh(meshjj,slicedir);
  }
}

/* ------------------ updateplotslice_mesh ------------------------ */

void updateplotslice_mesh(meshdata *mesh_in, int slicedir){
  int i, j, k;
  int plotx, ploty, plotz;
  int ibar, jbar, kbar;
  float *xplt, *yplt, *zplt;
  meshdata *meshi;
  contour *plot3dcontour1ptr, *plot3dcontour2ptr, *plot3dcontour3ptr;
  unsigned char *yzcolorbase, *xzcolorbase, *xycolorbase;
  float *yzcolorfbase, *xzcolorfbase, *xycolorfbase;
  float *yzcolortbase, *xzcolortbase, *xycolortbase;
  float *dx_xy, *dy_xy, *dz_xy;
  float *dx_xz, *dy_xz, *dz_xz;
  float *dx_yz, *dy_yz, *dz_yz;
  float *dx_yzcopy, *dy_yzcopy, *dz_yzcopy;
  float *dx_xzcopy, *dy_xzcopy, *dz_xzcopy;
  float *dx_xycopy, *dy_xycopy, *dz_xycopy;
  float *qdata;
  unsigned char *iqdata;
  char *iblank_xy=NULL, *iblank_xz=NULL, *iblank_yz=NULL;
  int nx, ny, nz, nxy, nxyz;
  char *c_iblank_x, *c_iblank_y, *c_iblank_z;
  float qval;

  meshi = mesh_in;
  plotx = meshi->iplotx_all[iplotx_all];
  ploty = meshi->iploty_all[iploty_all];
  plotz = meshi->iplotz_all[iplotz_all];

  ibar = meshi->ibar;
  jbar = meshi->jbar;
  kbar = meshi->kbar;
  xplt = meshi->xplt;
  yplt = meshi->yplt;
  zplt = meshi->zplt;
  c_iblank_x = meshi->c_iblank_x;
  c_iblank_y = meshi->c_iblank_y;
  c_iblank_z = meshi->c_iblank_z;

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
  nx = ibar + 1;
  ny = jbar + 1;
  nz = kbar + 1;
  nxy = nx*ny;
  nxyz = nx*ny*nz;


  iqdata = meshi->iqdata;
  qdata = meshi->qdata;

  minfill=1;
  maxfill=1;
  if(setp3min[plotn-1]==CHOP_MIN)minfill=0;
  if(setp3max[plotn-1]==CHOP_MAX)maxfill=0;

  if(ReadPlot3dFile!=1)return;
  if(plotx>=0&&slicedir==XDIR){
    float *yzcolorf;
    float *yzcolort;
    unsigned char *yzcolor;

    yzcolor=yzcolorbase;
    yzcolorf=yzcolorfbase;
    yzcolort=yzcolortbase;
    dx_yzcopy=dx_yz;
    dy_yzcopy=dy_yz;
    dz_yzcopy=dz_yz;
    iblank_yz=NULL;
    if(use_iblank==1&&c_iblank_x!=NULL){
      NewMemory((void **)&iblank_yz,(jbar+1)*(kbar+1)*sizeof(char));
      for(j=0;j<jbar;j++){
        for(k=0;k<kbar;k++){
          iblank_yz[k+j*kbar]=c_iblank_x[IJKNODE(plotx,j,k)];
        }
      }
    }
    for(j=0;j<=jbar;j++){
      for(k=0;k<=kbar;k++){
        *yzcolor++=iqdata[IJKN(plotx,j,k,plotn-1)];
        *yzcolort++=(float)iqdata[IJKN(plotx,j,k,plotn-1)]/255.0;
        GET_QVAL(plotx,j,k,plotn-1)
        *yzcolorf++=qval;
        *dx_yzcopy = 0.0;
        *dy_yzcopy = 0.0;
        *dz_yzcopy = 0.0;
        if(uindex!=-1){
          GET_QVAL(plotx,j,k,uindex)
          *dx_yzcopy++=vecfactor*veclength*qval/speedmax;
        }
        if(vindex!=-1){
          GET_QVAL(plotx,j,k,vindex)
          *dy_yzcopy++=vecfactor*veclength*qval/speedmax;
        }
        if(windex!=-1){
          GET_QVAL(plotx,j,k,windex)
          *dz_yzcopy++=vecfactor*veclength*qval/speedmax;
        }
      }
    }
    freecontour(plot3dcontour1ptr);
    initcontour(plot3dcontour1ptr,rgb_plot3d_contour,nrgb);
    setcontourslice(plot3dcontour1ptr,1,xplt[plotx]);
    getcontours(yplt,zplt,jbar+1,kbar+1,yzcolorfbase, iblank_yz,p3levels[plotn-1],DONT_GET_AREAS,DATA_FORTRAN,plot3dcontour1ptr);
    FREEMEMORY(iblank_yz);
  }
  else if(ploty>=0&&slicedir==YDIR){
    float *xzcolorf;
    float *xzcolort;
    unsigned char *xzcolor;

    xzcolor=xzcolorbase;
    xzcolorf=xzcolorfbase;
    xzcolort=xzcolortbase;
    dx_xzcopy=dx_xz;
    dy_xzcopy=dy_xz;
    dz_xzcopy=dz_xz;
    iblank_xz=NULL;
    if(use_iblank==1&&c_iblank_y!=NULL){
      NewMemory((void **)&iblank_xz,(ibar+1)*(kbar+1)*sizeof(char));
      for(i=0;i<ibar;i++){
        for(k=0;k<kbar;k++){
          iblank_xz[k+i*kbar]=c_iblank_y[IJKNODE(i,ploty,k)];
        }
      }
    }
    for(i=0;i<=ibar;i++){for(k=0;k<=kbar;k++){
      *xzcolor++=iqdata[IJKN(i,ploty,k,plotn-1)];
      *xzcolort++=(float)iqdata[IJKN(i,ploty,k,plotn-1)]/255.0;
      GET_QVAL(i,ploty,k,plotn-1)
      *xzcolorf++=qval;
      *dx_xzcopy = 0.0;
      *dy_xzcopy = 0.0;
      *dz_xzcopy = 0.0;
      if(uindex!=-1){
        GET_QVAL(i,ploty,k,uindex)
        *dx_xzcopy++=vecfactor*veclength*qval/speedmax;
      }
      if(vindex!=-1){
        GET_QVAL(i,ploty,k,vindex)
        *dy_xzcopy++=vecfactor*veclength*qval/speedmax;
      }
      if(windex!=-1){
        GET_QVAL(i,ploty,k,windex)
        *dz_xzcopy++=vecfactor*veclength*qval/speedmax;
      }
    }}
    freecontour(plot3dcontour2ptr);
    initcontour(plot3dcontour2ptr,rgb_plot3d_contour,nrgb);
    setcontourslice(plot3dcontour2ptr,2,yplt[ploty]);
    getcontours(xplt,zplt,ibar+1,kbar+1,xzcolorfbase, iblank_xz,p3levels[plotn-1],DONT_GET_AREAS,DATA_FORTRAN,plot3dcontour2ptr);
    FREEMEMORY(iblank_xz);
  }
  else if(plotz>=0&&slicedir==ZDIR){
    float *xycolorf;
    float *xycolort;
    unsigned char *xycolor;

    xycolor=xycolorbase;
    xycolorf=xycolorfbase;
    xycolort=xycolortbase;
    dx_xycopy=dx_xy;
    dy_xycopy=dy_xy;
    dz_xycopy=dz_xy;
    iblank_xy=NULL;
    if(use_iblank==1&&c_iblank_z!=NULL){
      NewMemory((void **)&iblank_xy,(ibar+1)*(jbar+1)*sizeof(char));
      for(i=0;i<ibar;i++){
        for(j=0;j<jbar;j++){
          iblank_xy[j+i*jbar]=c_iblank_z[IJKNODE(i,j,plotz)];
        }
      }
    }
    for(i=0;i<=ibar;i++){for(j=0;j<=jbar;j++){
      *xycolor++=iqdata[IJKN(i,j,plotz,plotn-1)];
      *xycolort++=(float)iqdata[IJKN(i,j,plotz,plotn-1)]/255.0;
      GET_QVAL(i,j,plotz,plotn-1)
      *xycolorf++=qval;
      *dx_xycopy = 0.0;
      *dy_xycopy = 0.0;
      *dz_xycopy = 0.0;
      if(uindex!=-1){
        GET_QVAL(i,j,plotz,uindex)
        *dx_xycopy++=vecfactor*veclength*qval/speedmax;
      }
      if(vindex!=-1){
        GET_QVAL(i,j,plotz,vindex)
        *dy_xycopy++=vecfactor*veclength*qval/speedmax;
      }
      if(windex!=-1){
        GET_QVAL(i,j,plotz,windex)
        *dz_xycopy++=vecfactor*veclength*qval/speedmax;
      }
    }}
    freecontour(plot3dcontour3ptr);
    initcontour(plot3dcontour3ptr,rgb_plot3d_contour,nrgb);
    setcontourslice(plot3dcontour3ptr,3,zplt[plotz]);
    getcontours(xplt,yplt,ibar+1,jbar+1,xycolorfbase, iblank_xy,p3levels[plotn-1],DONT_GET_AREAS,DATA_FORTRAN,plot3dcontour3ptr);
    FREEMEMORY(iblank_xy);
  }
}

/* ------------------ updateshowstep ------------------------ */

void updateshowstep(int val, int slicedir){

  if(ReadPlot3dFile==0&&visGrid==0&&ReadVolSlice==0)return;
  current_mesh->slicedir=slicedir;
  switch(slicedir){
  case XDIR:
    if(visx_all!=val)updatemenu=1;
    break;
  case YDIR:
    if(visy_all!=val)updatemenu=1;
    break;
  case ZDIR:
    if(visz_all!=val)updatemenu=1;
    break;
  case ISO:
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
  if(ReadVolSlice==0&&plotstate==DYNAMIC_PLOTS&&visGrid==0)Update_Times();
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
      meshdata *meshi;
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


      if(slicedir==XDIR&&(xmax-MESHEPS<xmin2||xmax2-MESHEPS<xmin))continue;
      if(slicedir==YDIR&&(ymax-MESHEPS<ymin2||ymax2-MESHEPS<ymin))continue;
      if(slicedir==ZDIR&&(zmax-MESHEPS<zmin2||zmax2-MESHEPS<zmin))continue;
    }
  }
}

/* ------------------ drawgrid ------------------------ */

void drawsphere(float diameter, unsigned char *rgbcolor);
void drawgrid(const meshdata *meshi){
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
  plotx = meshi->iplotx_all[iplotx_all];
  ploty = meshi->iploty_all[iploty_all];
  plotz = meshi->iplotz_all[iplotz_all];

  // visGrid
  // 0 no grid, no pointer
  // 1 grid, no pointer
  // 2 grid, pointer
  // 3 no grid, pointer

  if(visGrid==noGridProbe||visGrid==GridProbe){
    if(plotx>=0&&ploty>=0&&plotz>=0){
      unsigned char pcolor[4];

      glPushMatrix();
      glTranslatef(plotx_all[iplotx_all],ploty_all[iploty_all],plotz_all[iplotz_all]);
      pcolor[0]=255*foregroundcolor[0];
      pcolor[1]=255*foregroundcolor[1];
      pcolor[2]=255*foregroundcolor[2];
      drawsphere(0.03,pcolor);
      glPopMatrix();
    }
  }
  if(visx_all==0&&visy_all==0&&visz_all==0)return;
  if(visGrid==GridProbe||visGrid==GridnoProbe){
    antialias(ON);
    glLineWidth(gridlinewidth);
    if(meshi->meshrgb_ptr!=NULL){
      glColor3fv(meshi->meshrgb);
    }
    else{
      glColor4fv(foregroundcolor);
    }
    glBegin(GL_LINES);

    if(visx_all==1&&plotx>=0){
      for(j=0;j<jbar+1;j++){
        glVertex3f(xplt[plotx],yplt[j],zplt[0]);
        glVertex3f(xplt[plotx],yplt[j],zplt[kbar]);
      }
      for(k=0;k<kbar+1;k++){
        glVertex3f(xplt[plotx],yplt[0],zplt[k]);
        glVertex3f(xplt[plotx],yplt[jbar],zplt[k]);
      }
    }

    if(visy_all==1&&ploty>=0){
      for(i=0;i<ibar+1;i++){
        glVertex3f(xplt[i],yplt[ploty],zplt[0]);
        glVertex3f(xplt[i],yplt[ploty],zplt[kbar]);
      }
      for(k=0;k<kbar+1;k++){
        glVertex3f(   xplt[0],yplt[ploty],zplt[k]);
        glVertex3f(xplt[ibar],yplt[ploty],zplt[k]);
      }
    }

    if(visz_all==1&&plotz>=0){
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
    antialias(OFF);
  }
}

/* ------------------ update_plot3d_menulabels ------------------------ */

void update_plot3d_menulabels(void){
  int i;
  plot3ddata *plot3di;
  char label[128];

  if(nplot3dinfo>0){
    FREEMEMORY(plot3dorderindex);
    NewMemory((void **)&plot3dorderindex,sizeof(int)*nplot3dinfo);
    for(i=0;i<nplot3dinfo;i++){
      plot3dorderindex[i]=i;
    }
    qsort( (int *)plot3dorderindex, (size_t)nplot3dinfo, sizeof(int), plot3dcompare );

    for(i=0;i<nplot3dinfo;i++){
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
        meshdata *plot3dmesh;

        plot3dmesh = meshinfo + plot3di->blocknumber;
        sprintf(label,"%s",plot3dmesh->label);
        if(plot3di->time>=0.0)STRCAT(plot3di->menulabel,", ");
        sprintf(label,"%s",plot3dmesh->label);
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
  plot3ddata *plot3di;
  float lasttime_local,val;

  FREEMEMORY(plot3dtimelist);
  nplot3dtimelist=0;
  if(nplot3dinfo>0){
    NewMemory((void **)&plot3dtimelist,nplot3dinfo*sizeof(float));
  }
  if(plot3dtimelist==NULL)return;

  for(i=0;i<nplot3dinfo;i++){
    plot3di = plot3dinfo + i;
    plot3dtimelist[i]=plot3di->time;
  }
  qsort( (float *)plot3dtimelist, (size_t)nplot3dinfo, sizeof(int), plot3dlistcompare );
  lasttime_local=-999999.0;
  for(i=0;i<nplot3dinfo;i++){
    val=plot3dtimelist[i];
    if(ABS((double)(val-lasttime_local))>0.1){
      nplot3dtimelist++;
      plot3dtimelist[nplot3dtimelist-1]=val;
      lasttime_local=val;
    }
  }
}

/* ------------------ get_plot3d_uvw  ------------------------ */

void get_plot3d_uvw(float xyz[3], float uvw[3]){
  int i;
  float *xplt, *yplt, *zplt;

  uvw[0]=0.0;
  uvw[1]=0.0;
  uvw[2]=0.0;
  for(i=0;i<nmeshes;i++){
    meshdata *meshi;
    int ibar, jbar, kbar;
    float *udata,  *vdata, *wdata, *qdata;
    int nx, ny, nxy;
    int ix, iy, iz;
    int ijk;

    meshi = meshinfo + i;

    xplt = meshi->xplt_orig;
    yplt = meshi->yplt_orig;
    zplt = meshi->zplt_orig;
    if(xyz[0]<xplt[0]||xyz[1]<yplt[0]||xyz[2]<zplt[0])continue;

    ibar = meshi->ibar;
    jbar = meshi->jbar;
    kbar = meshi->kbar;
    if(xyz[0]>xplt[ibar-1]||xyz[1]>yplt[jbar-1]||xyz[2]>zplt[kbar-1])continue;

    nx = ibar+1;
    ny = jbar+1;
    nxy = nx*ny;

    qdata = meshi->qdata;
    udata = meshi->udata;
    vdata = meshi->vdata;
    wdata = meshi->wdata;

    if(qdata==NULL)continue;
    if(udata==NULL&&vdata==NULL&&wdata==NULL)continue;

    ix = ibar*(xyz[0]-xplt[0])/(xplt[ibar]-xplt[0]);
    if(ix<0)ix=0;
    if(ix>ibar)ix=ibar;

    iy = jbar*(xyz[1]-yplt[0])/(yplt[jbar]-yplt[0]);
    if(iy<0)iy=0;
    if(iy>jbar)iy=jbar;

    iz = kbar*(xyz[2]-zplt[0])/(zplt[kbar]-zplt[0]);
    if(iz<0)iz=0;
    if(iz>kbar)iz=kbar;

    ijk=IJKNODE(ix,iy,iz);

    if(udata!=NULL){
      uvw[0]=udata[ijk];
    }

    if(vdata!=NULL){
      uvw[1]=vdata[ijk];
    }

    if(wdata!=NULL){
      uvw[2]=wdata[ijk];
    }
  }
}

  /* ------------------ getplot3dtime ------------------------ */

int getplot3dtime(float *time_local){
  int i;

  for(i=0;i<nplot3dinfo;i++){
    plot3ddata *plot3di;

    plot3di = plot3dinfo + i;
    if(plot3di->loaded==0||plot3di->display==0)continue;
    *time_local = plot3di->time;
    return 1;
  }
  return 0;
}
