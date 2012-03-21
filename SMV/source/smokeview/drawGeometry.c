// $Date$ 
// $Revision$
// $Author$

// svn revision character string
char drawGeometry_revision[]="$Revision$";

#include "options.h"
#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "string_util.h"
#include "update.h"
#include "smokeviewvars.h"

cadgeom *current_cadgeom;

/* ------------------ UpdateIndexolors ------------------------ */

void UpdateIndexColors(void){
  int i;
  int colorindex;
  blockagedata *bc;
  mesh *meshi;
  int j;
  ventdata *vi;
  float s_color[4];
  surfdata *surfi;

  updateindexcolors=0;

  if(strcmp(surfacedefault->surfacelabel,"INERT")==0){
    surfacedefault->color=block_ambient2;
  }
  for(i=0;i<nsurfinfo;i++){
    surfi = surfinfo + i;
    if(strcmp(surfi->surfacelabel,"INERT")==0){
      surfi->color=block_ambient2;
    }
    if(strcmp(surfi->surfacelabel,"OPEN")==0){
      surfi->color=ventcolor;
    }
  }

  for(i=0;i<nmeshes;i++){
    meshi = meshinfo + i;
    for(j=0;j<meshi->nbptrs;j++){
      bc = meshi->blockageinfoptrs[j];
      if(bc->usecolorindex==1){
        colorindex=bc->colorindex;
        if(colorindex>=0){
          bc->color = getcolorptr(rgb[nrgb+colorindex]);
        }
      }
    }
    for(j=0;j<meshi->nvents;j++){
      vi = meshi->ventinfo + j;
      if(vi->usecolorindex==1){
        colorindex=vi->colorindex;

        s_color[0]=rgb[nrgb+colorindex][0];
        s_color[1]=rgb[nrgb+colorindex][1];
        s_color[2]=rgb[nrgb+colorindex][2];
        s_color[3]=1.0;
        vi->color = getcolorptr(s_color);
      }
    }
  }
  updatefaces=1;
}

/* ------------------ drawoutlines ------------------------ */

void drawoutlines(void){
  int i,j;
  outline *outlinei;
  float *xx1, *yy1, *zz1;
  float *xx2, *yy2, *zz2;

  if(noutlineinfo<=0)return;
  antialias(1);
  glLineWidth(linewidth);
  glBegin(GL_LINES);
  glColor3fv(foregroundcolor);
  for(i=0;i<noutlineinfo;i++){
    outlinei = outlineinfo + i;
    xx1 = outlinei->x1;
    xx2 = outlinei->x2;
    yy1 = outlinei->y1;
    yy2 = outlinei->y2;
    zz1 = outlinei->z1;
    zz2 = outlinei->z2;
    for(j=0;j<outlinei->nlines;j++){
      glVertex3f(*xx1++,*yy1++,*zz1++);
      glVertex3f(*xx2++,*yy2++,*zz2++);
    }
  }
  glEnd();
  antialias(0);
}
/* ------------------ drawcbox ------------------------ */

void drawcbox(float x, float y, float z, float size){
  float xx[8], yy[8], zz[8];
  float xbound[2],ybound[2],zbound[2];
  int i;
  int ix[8]={0,1,1,0,0,1,1,0};
  int iy[8]={0,0,0,0,1,1,1,1};
  int iz[8]={0,0,1,1,0,0,1,1};
  float dsize=size/xyzmaxdiff;

  xbound[0]=x-dsize/2.0;
  ybound[0]=y-dsize/2.0;
  zbound[0]=z-dsize/2.0;
  xbound[1]=x+dsize/2.0;
  ybound[1]=y+dsize/2.0;
  zbound[1]=z+dsize/2.0;
  for(i=0;i<8;i++){
    xx[i]=xbound[ix[i]];
    yy[i]=ybound[iy[i]];
    zz[i]=zbound[iz[i]];
  }

  glVertex3f(xx[0],yy[0],zz[0]);
  glVertex3f(xx[1],yy[1],zz[1]);
  glVertex3f(xx[5],yy[5],zz[5]);
  glVertex3f(xx[4],yy[4],zz[4]);

  glVertex3f(xx[1],yy[1],zz[1]);
  glVertex3f(xx[2],yy[2],zz[2]);
  glVertex3f(xx[6],yy[6],zz[6]);
  glVertex3f(xx[5],yy[5],zz[5]);

  glVertex3f(xx[2],yy[2],zz[2]);
  glVertex3f(xx[3],yy[3],zz[3]);
  glVertex3f(xx[7],yy[7],zz[7]);
  glVertex3f(xx[6],yy[6],zz[6]);

  glVertex3f(xx[3],yy[3],zz[3]);
  glVertex3f(xx[0],yy[0],zz[0]);
  glVertex3f(xx[4],yy[4],zz[4]);
  glVertex3f(xx[7],yy[7],zz[7]);

  glVertex3f(xx[0],yy[0],zz[0]);
  glVertex3f(xx[3],yy[3],zz[3]);
  glVertex3f(xx[2],yy[2],zz[2]);
  glVertex3f(xx[1],yy[1],zz[1]);

  glVertex3f(xx[4],yy[4],zz[4]);
  glVertex3f(xx[5],yy[5],zz[5]);
  glVertex3f(xx[6],yy[6],zz[6]);
  glVertex3f(xx[7],yy[7],zz[7]);
}

/* ------------------ get_blockvals ------------------------ */

void get_blockvals(  float *xmin, float *xmax, 
                     float *ymin, float *ymax, 
                     float *zmin, float *zmax,
                     int *imin, int *jmin, int *kmin){
  blockagedata *bc;

  bc=bchighlight;
  if(bc==NULL){
    *xmin = 0.0;
    *xmax = 0.0;
    *ymin = 0.0;
    *ymax = 0.0;
    *zmin = 0.0;
    *zmax = 0.0;
    *imin = 0;
    *jmin = 0;
    *kmin = 0;
    return;
  }
  if(blockage_as_input==1){
    float *xyz;

    xyz = bc->xyzEXACT;
    *xmin = xyz[0];
    *xmax = xyz[1];
    *ymin = xyz[2];
    *ymax = xyz[3];
    *zmin = xyz[4];
    *zmax = xyz[5];
  }
  else{
    *xmin = current_mesh->xplt_orig[bc->ijk[IMIN]]-current_mesh->offset[0];
    *xmax = current_mesh->xplt_orig[bc->ijk[IMAX]]-current_mesh->offset[0];
    *ymin = current_mesh->yplt_orig[bc->ijk[JMIN]]-current_mesh->offset[1];
    *ymax = current_mesh->yplt_orig[bc->ijk[JMAX]]-current_mesh->offset[1];
    *zmin = current_mesh->zplt_orig[bc->ijk[KMIN]]-current_mesh->offset[2];
    *zmax = current_mesh->zplt_orig[bc->ijk[KMAX]]-current_mesh->offset[2];
  }
  *imin = bc->ijk[IMIN];
  *jmin = bc->ijk[JMIN];
  *kmin = bc->ijk[KMIN];

}
#define ijkcell(i,j,k) ((i)+(j)*nx+(k)*nxy)

/*  
      if(iv1==iv2){
          vi->dir2=1;
          ventdir=DOWN_X;
          offset=ventoffset_factor*(xplttemp[1]-xplttemp[0]);
          voffset=-offset;
          if(inblockage(meshi,xmid-offset,ymid,zmid)==1){
            voffset=offset;
            ventdir=UP_X;
          }
          if(inblockage(meshi,xmid+offset,ymid,zmid)==1){
            voffset=-offset;
            ventdir=DOWN_X;
          }
          if(iv1==0){
            ventdir=UP_X;
            voffset=offset;
          }
          if(iv1==ibartemp){
            ventdir=DOWN_X;
            voffset=-offset;
          }
          if(nn<nvents)vi->dir=ventdir;
          if(vi->dummy==0){
            vi->xvent1 += voffset;
            vi->xvent2 += voffset;
          }
        }
        if(jv1==jv2){
          vi->dir2=2;
          ventdir=DOWN_Y;
          offset=ventoffset_factor*(yplttemp[1]-yplttemp[0]);
          voffset=0.0;
          if(inblockage(meshi,xmid,ymid-offset,zmid)==1){
            ventdir = UP_Y;
            voffset=offset;
          }
          if(inblockage(meshi,xmid,ymid+offset,zmid)==1){
            ventdir=DOWN_Y;
            voffset=-offset;
          }
          if(jv1==0){
            ventdir=UP_Y;
            voffset=offset;
          }
          if(jv1==jbartemp){
            ventdir=DOWN_Y;
            voffset=-offset;
          }
          if(vi->dummy==0){
            vi->yvent1 += voffset;
            vi->yvent2 += voffset;
          }
          if(nn<nvents)vi->dir=ventdir;
        }
        if(kv1==kv2){
          vi->dir2=3;
          offset=ventoffset_factor*(zplttemp[1]-zplttemp[0]);
          ventdir = DOWN_Z;
          voffset=0.0;
          if(inblockage(meshi,xmid,ymid,zmid-offset)==1){
            ventdir = UP_Z;
            voffset=offset;
          }
          if(inblockage(meshi,xmid,ymid,zmid+offset)==1){
            ventdir = DOWN_Z;
            voffset=-offset;
          }
          if(kv1==0){
            ventdir = UP_Z;
            voffset=offset;
          }
          if(kv1==kbartemp){
            ventdir = DOWN_Z;
            voffset=-offset;
          }
          if(vi->dummy==0){
            vi->zvent1 += voffset;
            vi->zvent2 += voffset;
          }
          if(nn<nvents)vi->dir=ventdir;
        }
*/
void setventdirs(void){
  mesh *meshi;
  ventdata *vi;
  int orien;
  int ii;
  int iv;
  int dir;
  int i, j, k;
  int index1,index2, index3;
  int nx, ny, nxy;
  char *iblank_x, *iblank_y, *iblank_z;
  int state1, state2, state3;
  int breakloop;
  int ventdir;
  float voffset, offset;
  float *xplttemp;
  float *yplttemp;
  float *zplttemp;

  for(ii=0;ii<nmeshes;ii++){
    meshi=meshinfo+ii;

    nx = meshi->ibar+1;
    ny = meshi->jbar+1;
    nxy = nx*ny;
    iblank_x = meshi->c_iblank_x;
    iblank_y = meshi->c_iblank_y;
    iblank_z = meshi->c_iblank_z;
    xplttemp=meshi->xplt;
    yplttemp=meshi->yplt;
    zplttemp=meshi->zplt;


    for(iv=0;iv<meshi->nvents+12;iv++){
      vi=meshi->ventinfo+iv;

      dir=0;
      if(vi->imin==vi->imax)dir=1;
      if(vi->jmin==vi->jmax)dir=2;
      if(vi->kmin==vi->kmax)dir=3;
      orien=0;

      switch (dir){
      case 1:
        vi->dir2=1;
        offset=ventoffset_factor*(xplttemp[1]-xplttemp[0]);
        if(vi->imin==0){
          orien=1;
        }
        else if(vi->imin==meshi->ibar){
          orien=-1;
        }
        else{
          orien=1;
          i=vi->imin;
          breakloop=0;
          for(j=vi->jmin;j<=vi->jmax;j++){
            for(k=vi->kmin;k<=vi->kmax;k++){
              index1=ijkcell(i-1,j,k);
              index2=ijkcell(i,j,k);
              index3=ijkcell(i+1,j,k);
              if(use_iblank==1){
                state1=iblank_x[index1];
                state2=iblank_x[index2];
                state3=iblank_x[index3];
              }
              else{
                state1=2;
                state2=2;
                state3=2;
              }
              if(state1==2&&state3==2)continue; // air on both sides
              if(state1==0&&state3==0)continue; // solid on both sides
              if(state2==1){
                if(state1==2){
                  orien=-1;
                }
                if(state3==2){
                  orien=1;
                }
              }
              breakloop=1;
              break;
            }
            if(breakloop==1)break;
          }
        }
        if(orien==1){
          ventdir=UP_X;
          voffset=offset;
        }
        else{
          ventdir=DOWN_X;
          voffset=-offset;
        }
        if(iv<meshi->nvents)vi->dir=ventdir;
        if(vi->dummy==0){
          vi->xvent1 += voffset;
          vi->xvent2 += voffset;
        }
        break;
      case 2:
        vi->dir2=2;
        offset=ventoffset_factor*(yplttemp[1]-yplttemp[0]);
        if(vi->jmin==0){
          orien=1;
        }
        else if(vi->jmin==meshi->jbar){
          orien=-1;
        }
        else{
          orien=1;
          j=vi->jmin;
          breakloop=0;
          for(i=vi->imin;i<=vi->imax;i++){
            for(k=vi->kmin;k<=vi->kmax;k++){
              index1=ijkcell(i,j-1,k);
              index2=ijkcell(i,j,k);
              index3=ijkcell(i,j+1,k);
              if(use_iblank==1){
                state1=iblank_y[index1];
                state2=iblank_y[index2];
                state3=iblank_y[index3];
              }
              else{
                state1=2;
                state2=2;
                state3=2;
              }
              if(state1==2&&state3==2)continue; // air on both sides
              if(state1==0&&state3==0)continue; // solid on both sides
              if(state2==1){
                if(state1==2){
                  orien=-1;
                }
                if(state3==2){
                  orien=1;
                }
              }
              breakloop=1;
              break;
            }
            if(breakloop==1)break;
          }
        }
        if(orien==1){
          ventdir=UP_Y;
          voffset=offset;
        }
        else{
          ventdir=DOWN_Y;
          voffset=-offset;
        }
        if(iv<meshi->nvents)vi->dir=ventdir;
        if(vi->dummy==0){
          vi->yvent1 += voffset;
          vi->yvent2 += voffset;
        }
        break;
      case 3:
        vi->dir2=3;
        offset=ventoffset_factor*(zplttemp[1]-zplttemp[0]);
        if(vi->kmin==0){
          orien=1;
        }
        else if(vi->kmin==meshi->kbar){
          orien=-1;
        }
        else{
          orien=1;
          k=vi->kmin;
          breakloop=0;
          for(i=vi->imin;i<=vi->imax;i++){
            for(j=vi->jmin;j<=vi->jmax;j++){
              index1=ijkcell(i,j,k-1);
              index2=ijkcell(i,j,k);
              index3=ijkcell(i,j,k+1);
              if(use_iblank==1){
                state1=iblank_z[index1];
                state2=iblank_z[index2];
                state3=iblank_z[index3];
              }
              else{
                state1=2;
                state2=2;
                state3=2;
              }
              if(state1==2&&state3==2)continue; // air on both sides
              if(state1==0&&state3==0)continue; // solid on both sides
              if(state2==1){
                if(state1==2){
                  orien=-1;
                }
                if(state3==2){
                  orien=1;
                }
              }
              breakloop=1;
              break;
            }
            if(breakloop==1)break;
          }
        }
        if(orien==1){
          ventdir=UP_Z;
          voffset=offset;
        }
        else{
          ventdir=DOWN_Z;
          voffset=-offset;
        }
        if(iv<meshi->nvents)vi->dir=ventdir;
        if(vi->dummy==0){
          vi->zvent1 += voffset;
          vi->zvent2 += voffset;
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
      }
    }
  }
}

/* ------------------ inblockage ------------------------ */

int inblockage(const mesh *meshi,float x, float y, float z){
  int i;
  float xmin, xmax, ymin, ymax, zmin, zmax;
  blockagedata *bc;
  float *xplt, *yplt, *zplt;

  xplt=meshi->xplt;
  yplt=meshi->yplt;
  zplt=meshi->zplt;

  for(i=0;i<meshi->nbptrs;i++){
    bc=meshi->blockageinfoptrs[i];
    xmin = xplt[bc->ijk[IMIN]]; xmax = xplt[bc->ijk[IMAX]];
    ymin = yplt[bc->ijk[JMIN]]; ymax = yplt[bc->ijk[JMAX]];
    zmin = zplt[bc->ijk[KMIN]]; zmax = zplt[bc->ijk[KMAX]];
    if(xmin<x && x<xmax && ymin<y && y<ymax && zmin<z && z<zmax){return(1);}
  }
  return(0);
}

/* ------------------ freecadinfo ------------------------ */

void freecadinfo(void){
  int i;
  cadgeom *cd;

  if(cadgeominfo!=NULL)return;
  if(ncadgeom>0){
    for(i=0;i<ncadgeom;i++){
      cd = cadgeominfo + i;
      FREEMEMORY(cd->quad);
      FREEMEMORY(cd->order);
    }
    FREEMEMORY(cadgeominfo);
    ncadgeom=0;
  }
}


/* ------------------ calcQuadNormal ------------------------ */

void calcQuadNormal(float *xyz, float *out){
  float u[3],v[3];
  static const int x = 0;
  static const int y = 1;
  static const int z = 2;
  float *p1, *p2, *p3;
  float *pp1, *pp2, *pp3, *pp4;

  pp1=xyz;
  pp2=xyz+3;
  pp3=xyz+6;
  pp4 = xyz+9;

  p1=pp1;
  p2=pp2;
  p3=pp3;

  if(pp1[0]==pp2[0]&&pp1[1]==pp2[1]&&pp1[2]==pp2[2]){
    p1=pp2;
    p2=pp3;
    p3=pp4;
  }
  if(pp2[0]==pp3[0]&&pp2[1]==pp3[1]&&pp2[2]==pp3[2]){
    p1=pp1;
    p2=pp3;
    p3=pp4;
  }

  u[x] = p2[x] - p1[x];
  u[y] = p2[y] - p1[y];
  u[z] = p2[z] - p1[z];

  v[x] = p3[x] - p1[x];
  v[y] = p3[y] - p1[y];
  v[z] = p3[z] - p1[z];

  out[x] = u[y]*v[z] - u[z]*v[y];
  out[y] = u[z]*v[x] - u[x]*v[z];
  out[z] = u[x]*v[y] - u[y]*v[x];

  ReduceToUnit(out);

}


/* ------------------ readcadgeom ------------------------ */

void readcadgeom(cadgeom *cd){
  char buffer[255];
  char obstlabel[255];
  float lastcolor[3];
  FILE *stream;
  int nquads=0;
  float *xyzpoints;
  int colorindex;
  float *normal;
  char *colors;
  float rgbtemp[4]={(float)-1.,(float)-1.,(float)-1.,(float)1.};
  cadquad *quadi;

  if( (stream=fopen(cd->file,"r"))==NULL){
    return;
  }
  if(fgets(buffer,255,stream)==NULL)return;
  trim(buffer);
  if(strncmp(buffer,"[APPEARANCE]",12)==0){
    cd->version=2;
    fclose(stream);
    readcad2geom(cd);
    return;
  }
  else{
    cd->version=1;
    rewind(stream);
  }
  while(!feof(stream)){
    if(fgets(buffer,255,stream)==NULL)break;
    if(fgets(buffer,255,stream)==NULL)break;
    nquads++;
  }
  cd->nquads=nquads;
  rewind(stream);
  if(NewMemory((void **)&cd->quad,nquads*sizeof(cadquad))==0){
    freecadinfo();
    return;
  }
  nquads=0;
  colorindex=0;
  lastcolor[0]=-1.0;
  lastcolor[1]=-1.0;
  lastcolor[2]=-1.0;
  while(!feof(stream)){
    if(fgets(buffer,255,stream)==NULL)break;
    quadi = cd->quad + nquads;
    xyzpoints = quadi->xyzpoints;
    normal = quadi->normals;
    sscanf(buffer,"%f %f %f %f %f %f %f %f %f %f %f %f ",
      xyzpoints, xyzpoints+1, xyzpoints+2,
      xyzpoints+3, xyzpoints+4, xyzpoints+5,
      xyzpoints+6, xyzpoints+7, xyzpoints+8,
      xyzpoints+9,xyzpoints+10,xyzpoints+11
      );
    calcQuadNormal(xyzpoints, normal);

    if(fgets(buffer,255,stream)==NULL)break;
    colors=strstr(buffer," ");
    strcpy(obstlabel,buffer);
    rgbtemp[0]=(float)-1.0;
    rgbtemp[1]=(float)-1.0;
    rgbtemp[2]=(float)-1.0;
    rgbtemp[3]=(float)1.0;
    if(colors!=NULL){
      colors[0]='\0';
      sscanf(colors+1,"%f %f %f",rgbtemp,rgbtemp+1,rgbtemp+2);
    }
    if(lastcolor[0]!=rgbtemp[0]||lastcolor[1]!=rgbtemp[1]||lastcolor[2]!=rgbtemp[2]){
      quadi->colorindex=colorindex;
      colorindex++;
      lastcolor[0]=rgbtemp[0];
      lastcolor[1]=rgbtemp[1];
      lastcolor[2]=rgbtemp[2];
    }
    else{
      quadi->colorindex=colorindex;
    }
    if(colors!=NULL&&rgbtemp[0]>=0.0){
      quadi->colorindex=-1;
    }
    quadi->colors[0]=rgbtemp[0];
    quadi->colors[1]=rgbtemp[1];
    quadi->colors[2]=rgbtemp[2];
    quadi->colors[3]=rgbtemp[3];
    nquads++;
  }
  fclose(stream);

}

/* ------------------ quadcompare ------------------------ */

int quadcompare( const void *arg1, const void *arg2 ){
  int i1, i2;
  cadgeom *cd;
  cadquad *quadi, *quadj;
  cadlook *cli, *clj;

  cd=current_cadgeom;

  i1 = *(int *)arg1;
  i2 = *(int *)arg2;

  quadi = cd->quad + i1;
  cli = quadi->cadlookq;

  quadj = cd->quad + i2;
  clj = quadj->cadlookq;

  if(cli<clj)return 1;
  if(cli>clj)return -1;
  return 0;
}

/* ------------------ readcad2geom ------------------------ */

void readcad2geom(cadgeom *cd){
  char buffer[255];
  FILE *stream;
  int nquads=0;
  float *normal;
  int i,look_index;
  cadquad *quadi;
  cadlook *cl;
  float *rrgb;
  int iquad;
  float *xyzpoints;
  int have_textures=0;

  if( (stream=fopen(cd->file,"r"))==NULL){
    return;
  }

  /* read in [APPEARANCE] info */

  if(fgets(buffer,255,stream)==NULL)return;
  if(fgets(buffer,255,stream)==NULL)return;
  sscanf(buffer,"%i",&cd->ncadlookinfo);
  if(cd->ncadlookinfo<=0){
    cd->ncadlookinfo=0;
    return;
  }

  cd->cadlookinfo=NULL;
  NewMemory((void **)&cd->cadlookinfo,cd->ncadlookinfo*sizeof(cadlook));

  for(i=0;i<cd->ncadlookinfo;i++){
    cadlook *cdi;
    texture *texti;
    int errorcode;
    int ii;
    size_t lenbuffer,len;
    float *shininess;
    float *t_origin;

    cdi = cd->cadlookinfo + i;

    if(fgets(buffer,255,stream)==NULL)return; // material description (not used)

    if(fgets(buffer,255,stream)==NULL)return;
    rrgb=cdi->rgb;
    shininess=&cdi->shininess;
    cdi->texture_height=-1.0;
    cdi->texture_width=-1.0;
    t_origin = cdi->texture_origin;
    t_origin[0]=0.0;
    t_origin[1]=0.0;
    t_origin[2]=0.0;
    rrgb[0]=-255.0;
    rrgb[1]=-255.0;
    rrgb[2]=-255.0;
    rrgb[3]=1.0;
    *shininess=block_shininess;
    lenbuffer=strlen(buffer);
    for(ii=0;ii<lenbuffer;ii++){
      if(buffer[ii]==',')buffer[ii]=' ';
    }
    sscanf(buffer,"%i %f %f %f %f %f %f %f %f %f %f",
      &cdi->index,rrgb,rrgb+1,rrgb+2,
      &cdi->texture_width,&cdi->texture_height,
      rrgb+3,shininess,
      t_origin,t_origin+1,t_origin+2
      );
      
    rrgb[0]/=255.0;
    rrgb[1]/=255.0;
    rrgb[2]/=255.0;

    if(fgets(buffer,255,stream)==NULL)return;
    trim(buffer);
    len=strlen(buffer);

    texti = &cdi->textureinfo;
    texti->file=NULL;
    if(len>0){
      NewMemory((void **)&texti->file,len+1);
      strcpy(texti->file,buffer);
    }
    texti->display=0;
    texti->loaded=0;
    texti->used=0;
    texti->name=0;

    if(texti->file!=NULL){
      int texwid, texht;
      unsigned char *floortex;

      if(have_textures==0){
        printf("     Loading CAD textures\n");
        have_textures=1;
      }
      printf("       Loading texture: %s",texti->file);
      glGenTextures(1,&texti->name);
      glBindTexture(GL_TEXTURE_2D,texti->name);
      floortex=readpicture(texti->file,&texwid,&texht,0);
      if(floortex==NULL){
        printf(" - failed\n");
        continue;
      }
      errorcode=gluBuild2DMipmaps(GL_TEXTURE_2D,4, texwid, texht, GL_RGBA, GL_UNSIGNED_BYTE, floortex);
      if(errorcode!=0){
        FREEMEMORY(floortex);
        printf(" - failed\n");
        continue;
      }
      FREEMEMORY(floortex);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
      texti->loaded=1;
      printf(" - completed\n");
    }
  }

  /* read in [FACES] info */

  if(fgets(buffer,255,stream)==NULL)return;
  if(fgets(buffer,255,stream)==NULL)return;
  sscanf(buffer,"%i",&nquads);
  if(nquads<=0){
    cd->nquads=0;
    return;
  }
  cd->nquads=nquads;
  cd->order=NULL;
  if(NewMemory((void **)&cd->quad,nquads*sizeof(cadquad))==0||
     NewMemory((void **)&cd->order,nquads*sizeof(int))==0
    ){
    freecadinfo();
    return;
  }

  iquad=0;
  for(i=0;i<nquads;i++){
    if(fgets(buffer,255,stream)==NULL)break;
    iquad++;
    quadi = cd->quad + i;
    xyzpoints = quadi->xyzpoints;
    normal = quadi->normals;
    sscanf(buffer,"%f %f %f %f %f %f %f %f %f %f %f %f %i",
      xyzpoints,xyzpoints+1,xyzpoints+2,
      xyzpoints+3,xyzpoints+4,xyzpoints+5,
      xyzpoints+6,xyzpoints+7,xyzpoints+8,
      xyzpoints+9,xyzpoints+10,xyzpoints+11,&look_index
      );
    if(look_index<0||look_index>cd->ncadlookinfo-1)look_index=0;
    quadi->cadlookq=cd->cadlookinfo+look_index;
    cl = quadi->cadlookq;
    quadi->colors[0]=cl->rgb[0];
    quadi->colors[1]=cl->rgb[1];
    quadi->colors[2]=cl->rgb[2];
    quadi->colors[3]=1.0;
    calcQuadNormal(xyzpoints, normal);
  }
  if(iquad<nquads){
    printf("  *** warning: number of faces expected=%i number of faces found=%i\n",cd->nquads,iquad);
    cd->nquads=iquad;
  }
  for(i=0;i<cd->nquads;i++){
    cd->order[i]=i;
  }
  current_cadgeom=cd;
  qsort(cd->order,(size_t)cd->nquads,sizeof(int),quadcompare);
  fclose(stream);
  if(have_textures==1){
    printf("     CAD textures loading completed\n");
  }
}

/* ------------------ updaate_cadtextcoords ------------------------ */

void update_cadtextcoords(cadquad *quadi){
  float twidth, theight;
  float nx, ny, nz;
  float l1, l2;
  int i;
  float qx, qy, qz;
  float *xyz;
  float *txy;
  float *t_origin;

  xyz = quadi->xyzpoints;
  txy = quadi->txypoints;

  twidth = quadi->cadlookq->texture_width;
  t_origin = quadi->cadlookq->texture_origin;
  if(twidth==0.0)twidth=1.0;
  theight = quadi->cadlookq->texture_height;
  if(theight==0.0)theight=1.0;
  nx=quadi->normals[0];
  ny=quadi->normals[1];
  nz=quadi->normals[2];
  l1 = sqrt(nx*nx+ny*ny);
  l2 = l1*sqrt(nx*nx+ny*ny+nz*nz);


  for(i=0;i<4;i++){
    qx=xbar0+xyz[3*i]*xyzmaxdiff - t_origin[0];
    qy=ybar0+xyz[3*i+1]*xyzmaxdiff - t_origin[1];
    qz=zbar0+xyz[3*i+2]*xyzmaxdiff - t_origin[2];

    if(l1!=0.0){
      txy[2*i]=(ny*qx-nx*qy)/l1;
      txy[2*i+1]=(nx*nz*qx+ny*nz*qy-(nx*nx+ny*ny)*qz)/l2;
    }
    else{
      txy[2*i]=-qx;
      txy[2*i+1]=-qy;
    }
    txy[2*i]/=twidth;
    txy[2*i+1]/=theight;
    txy[2*i]=1.0-txy[2*i];
    txy[2*i+1]=1.0-txy[2*i+1];
  }




}

/* ------------------ newtextptr ------------------------ */

char *newtextptr(char ***list,int *nlist,char *this,char *last){
  int i,n;
  char *thisname;
  unsigned int n_thisname;
  char **from_list, **to_list;

  if(last!=NULL&&strcmp(this,last)==0)return last;
  if(this==NULL)return NULL;

  from_list=*list;
  n=*nlist;
  for(i=0;i<n;i++){
    if(strcmp(from_list[i],this)==0)return from_list[i];
  }
  
  n++;
  NewMemory((void **)&to_list,n*sizeof(char **));
  for(i=0;i<n-1;i++){
    to_list[i]=from_list[i];
  }
  n_thisname=strlen(this);
  NewMemory((void *)&thisname,n_thisname+1);
  strcpy(thisname,this);
  to_list[n-1]=thisname;
  FREEMEMORY(from_list);

  *nlist=n;
  *list=to_list;
  return thisname;
}

/* ------------------ drawcadgeom ------------------------ */

void drawcadgeom(const cadgeom *cd){
  int i;
  float *xyzpoint;
  float *normal;
  int last_colorindex=-999;
  int colorindex;
  int colorindex2;
  float *thiscolor,*lastcolor; 
  float rgbtemp[4]={(float)-1.0,(float)-1.0,(float)-1.0,(float)-1.0};
  cadquad *quadi;

  lastcolor=rgbtemp;
  if(cullfaces==1)glDisable(GL_CULL_FACE);

  glEnable(GL_LIGHTING); 
  glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
  glEnable(GL_COLOR_MATERIAL);
  glBegin(GL_QUADS);
  for(i=0;i<cd->nquads;i++){
    quadi = cd->quad+i;
    colorindex=quadi->colorindex;
    thiscolor=quadi->colors;
    if(colorindex!=last_colorindex||
      thiscolor[0]!=lastcolor[0]||
      thiscolor[1]!=lastcolor[1]||
      thiscolor[2]!=lastcolor[2]
      ){
      if(colorindex==-1){
        thiscolor=quadi->colors;
      }
      else{
        colorindex2 = 15 + (15*colorindex % 230);
        thiscolor=rgb_cad[colorindex2];
      }
      glColor4fv(thiscolor);
    }
    last_colorindex=colorindex;
    lastcolor=thiscolor;
    xyzpoint = quadi->xyzpoints;
    normal = quadi->normals;

    glNormal3fv(normal);
    glVertex3fv(xyzpoint);
    glVertex3fv(xyzpoint+3);
    glVertex3fv(xyzpoint+6);
    glVertex3fv(xyzpoint+9);
  }
  glEnd();
  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);
  SNIFF_ERRORS("drawcadgeom");
  if(cullfaces==1)glEnable(GL_CULL_FACE);

}

/* ------------------ drawcadgeom ------------------------ */

void drawcad2geom(const cadgeom *cd, int trans_flag){
  int ii,i;
  float *xyzpoint;
  float *txypoint;
  float *normal;
  float *thiscolor,*lastcolor; 
  cadquad *quadi;
  int colorindex,colorindex2;
  texture *lasttexture;
  float last_block_shininess;

  lastcolor=NULL;
  last_block_shininess=-1.0;
  if(cullfaces==1)glDisable(GL_CULL_FACE);

  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_LIGHTING); 
  if(trans_flag==DRAW_TRANSPARENT)transparenton();
  glBegin(GL_QUADS);
  colorindex=0;
  for(ii=0;ii<cd->nquads;ii++){
    texture *texti;
    float this_block_shininess;

    i=cd->order[ii];
    ASSERT(i>=0&&i<cd->nquads);
    quadi = cd->quad+i;
    xyzpoint = quadi->xyzpoints;
    texti = &quadi->cadlookq->textureinfo;

    if(visCadTextures==1&&texti->loaded==1)continue;

    thiscolor=quadi->cadlookq->rgb;
    if(thiscolor!=lastcolor){
      if(thiscolor[0]<0.0){
        GLfloat *colorptr;

        colorindex2 = 15 + (15*colorindex % 230);
        colorptr = &rgb_cad[colorindex2][0];
        glColor4fv(colorptr);
        colorindex++;
      }
      else{
        if((thiscolor[3]<1.0&&trans_flag!=DRAW_TRANSPARENT)||
           (thiscolor[3]>=1.0&&trans_flag==DRAW_TRANSPARENT)){
           continue;
        }
        glColor4fv(thiscolor);
      }
      lastcolor=thiscolor;
    }

    if(RectangleInFrustum(xyzpoint,xyzpoint+3,xyzpoint+6,xyzpoint+9)==0)continue;

    this_block_shininess = quadi->cadlookq->shininess;
    if(last_block_shininess!=this_block_shininess){
      last_block_shininess=this_block_shininess;
      glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&this_block_shininess);
    }
    
    normal = quadi->normals;

    glNormal3fv(normal);
    glVertex3fv(xyzpoint);
    glVertex3fv(xyzpoint+3);
    glVertex3fv(xyzpoint+6);
    glVertex3fv(xyzpoint+9);
  }
  glEnd();
  if(visCadTextures==1){
    glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
    glEnable(GL_TEXTURE_2D);

    lasttexture=NULL;
    glBegin(GL_QUADS);
    for(ii=0;ii<cd->nquads;ii++){
      texture *texti;
    
      i=cd->order[ii];
      ASSERT(i>=0&&i<cd->nquads);
      quadi = cd->quad+i;
      xyzpoint = quadi->xyzpoints;

      texti=&quadi->cadlookq->textureinfo;
      if(texti->loaded==0)continue;

      txypoint = quadi->txypoints;
      normal = quadi->normals;

      if(lasttexture!=texti){
        glEnd();
        glBindTexture(GL_TEXTURE_2D,texti->name);
        glBegin(GL_QUADS);
        lasttexture=texti;
      }

      if(RectangleInFrustum(xyzpoint,xyzpoint+3,xyzpoint+6,xyzpoint+9)==0)continue;

      glNormal3fv(normal);
      glTexCoord2fv(txypoint);
      glVertex3fv(xyzpoint);

      glTexCoord2fv(txypoint+2);
      glVertex3fv(xyzpoint+3);

      glTexCoord2fv(txypoint+4);
      glVertex3fv(xyzpoint+6);

      glTexCoord2fv(txypoint+6);
      glVertex3fv(xyzpoint+9);
    }
    glEnd();
    glDisable(GL_TEXTURE_2D);
  }

  glDisable(GL_LIGHTING);
  glDisable(GL_COLOR_MATERIAL);
  if(trans_flag==DRAW_TRANSPARENT)transparentoff();
  if(cullfaces==1)glEnable(GL_CULL_FACE);

}

/* ------------------ drawcad2ageom ------------------------ */

void drawcad2ageom(const cadgeom *cd){
  int i;
  float *xyzpoint;
  float *thiscolor,*lastcolor; 
  cadquad *quadi;
  int colorindex,colorindex2;
  float *color2;

  lastcolor=NULL;

  glBegin(GL_LINES);
  colorindex=0;
  for(i=0;i<cd->nquads;i++){

    quadi = cd->quad+i;

    thiscolor=quadi->cadlookq->rgb;
    if(thiscolor!=lastcolor){
      if(thiscolor[0]<0.0){
        colorindex2 = 15 + (15*colorindex % 230);
        color2=rgb_cad[colorindex2];
        glColor3f(color2[0],color2[1],color2[2]);
        colorindex++;
      }
      else{
        glColor3f(thiscolor[0],thiscolor[1],thiscolor[2]);
      }
      lastcolor=thiscolor;
    }
    xyzpoint = quadi->xyzpoints;

    glVertex3fv(xyzpoint);
    glVertex3fv(xyzpoint+3);
    glVertex3fv(xyzpoint+3);
    glVertex3fv(xyzpoint+6);
    glVertex3fv(xyzpoint+6);
    glVertex3fv(xyzpoint+9);
    glVertex3fv(xyzpoint+9);
    glVertex3fv(xyzpoint);
    glVertex3fv(xyzpoint);
  }
  glEnd();
}

/* ------------------ updatefaces ------------------------ */

void update_faces(void){
  int i,j;
  mesh *meshi;
  blockagedata *bc;
  facedata *faceptr;
  ventdata *vi;
  allocate_faces();
  updatefaces=0;
  have_vents_int=0;
  for(i=0;i<nmeshes;i++){
    meshi = meshinfo + i;
    faceptr = meshi->faceinfo;
    for(j=0;j<meshi->nbptrs;j++){
      bc = meshi->blockageinfoptrs[j];
      if(visTerrainType!=TERRAIN_HIDDEN&&bc->is_wuiblock==1)continue;
      obst_or_vent2faces(meshi,bc,NULL,faceptr,BLOCK_face);
      faceptr += 6;
    }
    for(j=0;j<meshi->nvents;j++){
      vi = meshi->ventinfo+j;
      obst_or_vent2faces(meshi,NULL,vi,faceptr,VENT_face);
      faceptr++;
    }
    for(j=meshi->nvents;j<meshi->nvents+6;j++){
      vi = meshi->ventinfo+j;
      obst_or_vent2faces(meshi,NULL,vi,faceptr,OUTLINE_FRAME_face);
      ASSERT(faceptr->color!=NULL);
      faceptr++;
    }
    for(j=meshi->nvents+6;j<meshi->nvents+12;j++){
      vi = meshi->ventinfo+j;
      obst_or_vent2faces(meshi,NULL,vi,faceptr,SHADED_FRAME_face);
      ASSERT(faceptr->color!=NULL);
      faceptr++;
    }
    meshi->nfaces=faceptr-meshi->faceinfo;
  }
  UpdateHiddenFaces();
  update_facelists();
  update_selectfaces();
  update_selectblocks();
}

/* ------------------ obst_or_vent2faces ------------------------ */

void obst_or_vent2faces(const mesh *meshi,blockagedata *bc, 
                        ventdata *vi, facedata *faceptr, int facetype){
  /*
         
       7---------6
     /         /
   /         /
  4--------5
                 
       3 ------  2
      /         /
    /         /
  0 ------ 1

  */
  int n,j,k;
  int jend=1,jjj;
  float xminmax[2]={0.0,1.0}, xminmax2[2]={0.0,1.0};
  float yminmax[2]={0.0,1.0}, yminmax2[2]={0.0,1.0};
  float zminmax[2]={0.0,1.0}, zminmax2[2]={0.0,1.0};
  float xx[8], yy[8], zz[8];
  float xx2[8], yy2[8], zz2[8];
  float *xplt, *yplt, *zplt;
  int ii[8] = {0,1,1,0,0,1,1,0};
  int jj[8] = {0,0,1,1,0,0,1,1};
  int kk[8] = {0,0,0,0,1,1,1,1};
  float offset[3];

  int blockfaceindex[]={
               0,1,5,4,  /* down y */
               1,2,6,5,  /*   up x */
               2,3,7,6,  /*   up y */
               3,0,4,7,  /* down x */
               3,2,1,0,  /* down z */
               4,5,6,7}; /*   up z */
  int *bfi;
  float *xtex, *ytex;
  float *xtex2, *ytex2;
  float t_width, t_height;
  float *xstart, *ystart;

  surfdata *face_surf;

  ASSERT((bc==NULL&&vi!=NULL)||(bc!=NULL&&vi==NULL));

  xplt = meshi->xplt;
  yplt = meshi->yplt;
  zplt = meshi->zplt;
  ASSERT(bc!=NULL&&vi==NULL||bc==NULL&&vi!=NULL);
  if(bc!=NULL){
    jend=6;
    xminmax[0] = xplt[bc->ijk[IMIN]]; 
    xminmax[1] = xplt[bc->ijk[IMAX]];
    yminmax[0] = yplt[bc->ijk[JMIN]]; 
    yminmax[1] = yplt[bc->ijk[JMAX]];
    zminmax[0] = zplt[bc->ijk[KMIN]]; 
    zminmax[1] = zplt[bc->ijk[KMAX]];

    xminmax2[0] = bc->xmin; 
    xminmax2[1] = bc->xmax;
    yminmax2[0] = bc->ymin; 
    yminmax2[1] = bc->ymax;
    zminmax2[0] = bc->zmin; 
    zminmax2[1] = bc->zmax;
  }
  if(vi!=NULL){
    jend=1;
    xminmax[0] = xplt[vi->imin]; 
    xminmax[1] = xplt[vi->imax];
    yminmax[0] = yplt[vi->jmin]; 
    yminmax[1] = yplt[vi->jmax];
    zminmax[0] = zplt[vi->kmin]; 
    zminmax[1] = zplt[vi->kmax];

    xminmax2[0] = vi->xmin; 
    xminmax2[1] = vi->xmax;
    yminmax2[0] = vi->ymin; 
    yminmax2[1] = vi->ymax;
    zminmax2[0] = vi->zmin; 
    zminmax2[1] = vi->zmax;
  }


  for(n=0;n<8;n++){
    xx[n]=xminmax[ii[n]];
    yy[n]=yminmax[jj[n]];
    zz[n]=zminmax[kk[n]];
    xx2[n]=xminmax2[ii[n]];
    yy2[n]=yminmax2[jj[n]];
    zz2[n]=zminmax2[kk[n]];
  }

  for(j=0;j<jend;j++){
    faceptr->meshindex=meshi-meshinfo;
    faceptr->type2=facetype;
    faceptr->thinface=0;
    faceptr->is_interior=0;
    faceptr->show_bothsides=0;
    faceptr->bc=NULL;
    
    if(bc!=NULL){
      faceptr->bc=bc;
      faceptr->hidden=0;
      faceptr->patchpresent=0;
      faceptr->blockageindex=-2;
      if(visBlocks==visBLOCKSolidOutline){
        faceptr->linewidth=&solidlinewidth;
      }
      else{
        faceptr->linewidth=&linewidth;
      }
      faceptr->showtimelist_handle=&bc->showtimelist;
      faceptr->del=bc->del;
      faceptr->invisible=bc->invisible;
      faceptr->surfinfo=bc->surf[j];
      faceptr->texture_origin=bc->texture_origin;
      faceptr->transparent=bc->transparent;
    }
    if(vi!=NULL){
      faceptr->blockageindex=-2;
      faceptr->hidden=0;
      faceptr->patchpresent=0;
      faceptr->del=0;
      faceptr->invisible=0;
      faceptr->texture_origin=vi->texture_origin;
      faceptr->transparent=vi->transparent;
      if(faceptr->type2==OUTLINE_FRAME_face){
        faceptr->linewidth=&linewidth;
      }
      else{
        faceptr->linewidth=&ventlinewidth;
      }
      faceptr->showtimelist_handle=&vi->showtimelist;
      faceptr->surfinfo=vi->surf[j];
    }
    face_surf = faceptr->surfinfo;
    if(face_surf!=NULL){
      t_width =face_surf->t_width;
      t_height=face_surf->t_height;
      if(bc!=NULL){
        if(bc->type==-1){
          faceptr->type=face_surf->type;
        }
        else{
          faceptr->type=bc->type;
        }
        if(face_surf->textureinfo!=NULL){
          faceptr->type=face_surf->type;
        }
      }
      if(vi!=NULL){
        if(vi->type==99||vi->type==-99){
          faceptr->type=face_surf->type;
        }
        else{
          faceptr->type=vi->type;
        }
      }
    }
    else{
      if(bc!=NULL)faceptr->type=bc->type;
      if(vi!=NULL)faceptr->type=vi->type;
      t_width=1.0;
      t_height=1.0;
    }
    if(t_width==0.0)t_width=1.0;
    if(t_height==0.0)t_height=1.0;
    if(bc!=NULL){
      switch (bc->useblockcolor){
      case 1:
        faceptr->color=bc->color;
        faceptr->transparent=bc->transparent;
        break;
      case 0:
        if(bc->surf[j]==surfacedefault){
         // faceptr->color=block_ambient2;
          faceptr->color=surfacedefault->color;  /* fix ?? */
          faceptr->invisible=surfacedefault->invisible;
          faceptr->transparent=surfacedefault->transparent;
        }
        else{
          faceptr->color=bc->surf[j]->color;
          faceptr->invisible=bc->surf[j]->invisible;
          faceptr->transparent=bc->surf[j]->transparent;
        }
        break;
      default:
        ASSERT(FFALSE);
        break;
      }
      if(visBlocks==visBLOCKSolidOutline){
        faceptr->linecolor=foregroundcolor;
      }
      else{
        faceptr->linecolor=faceptr->color;
      }
    }
    if(vi!=NULL){
      if(vi->useventcolor==1){
        faceptr->color=vi->color;
        if(facetype!=OUTLINE_FRAME_face)faceptr->transparent=vi->transparent;
      }
      else if(facetype==OUTLINE_FRAME_face){
        if(meshi->meshrgb_ptr!=NULL){
          faceptr->color=meshi->meshrgb_ptr;
        }
        else{
          faceptr->color=vi->color;
        }
      }
      else{
        faceptr->color=vi->surf[j]->color;
        faceptr->transparent=vi->surf[j]->transparent;
      }
      faceptr->linecolor=faceptr->color;
    }


    if(bc!=NULL){
      faceptr->textureinfo=bc->surf[j]->textureinfo;
    }
    if(vi!=NULL){
      faceptr->textureinfo=vi->surf[j]->textureinfo;
    }
    if(bc!=NULL){
      faceptr->dir=j;
    }
    if(vi!=NULL){
      faceptr->dir=vi->dir;
    }
    faceptr->normal[0]=0.0;
    faceptr->normal[1]=0.0;
    faceptr->normal[2]=0.0;
    if(bc!=NULL){
      faceptr->imin=bc->ijk[IMIN];
      faceptr->imax=bc->ijk[IMAX];
      faceptr->jmin=bc->ijk[JMIN];
      faceptr->jmax=bc->ijk[JMAX];
      faceptr->kmin=bc->ijk[KMIN];
      faceptr->kmax=bc->ijk[KMAX];
    }
    if(vi!=NULL){
      faceptr->imin=vi->imin;
      faceptr->imax=vi->imax;
      faceptr->jmin=vi->jmin;
      faceptr->jmax=vi->jmax;
      faceptr->kmin=vi->kmin;
      faceptr->kmax=vi->kmax;
      if(faceptr->imin==faceptr->imax){
        if(faceptr->imin>0&&faceptr->imin<meshi->ibar)faceptr->is_interior=1;
      }
      if(faceptr->jmin==faceptr->jmax){
        if(faceptr->jmin>0&&faceptr->jmin<meshi->jbar)faceptr->is_interior=1;
      }
      if(faceptr->kmin==faceptr->kmax){
        if(faceptr->kmin>0&&faceptr->kmin<meshi->kbar)faceptr->is_interior=1;
      }
      if(faceptr->is_interior==1)have_vents_int=1;
      if(faceptr->is_interior==1&&show_bothsides_int==1)faceptr->show_bothsides=1;
      if(faceptr->is_interior==0&&show_bothsides_ext==1)faceptr->show_bothsides=1;
    }
    offset[0]=(float)0.0;
    offset[1]=(float)0.0;
    offset[2]=(float)0.0;
    switch (faceptr->dir) {
     case DOWN_Y: 
       faceptr->normal[1]=(float)-1.0;
       if(facetype==VENT_face&&vi!=NULL&&vi->dummy==0)offset[1] = -meshi->vent_offset[1];
       faceptr->jmax=faceptr->jmin;
       if(xminmax2[0]==xminmax2[1]||zminmax2[0]==zminmax2[1])faceptr->thinface=1;
       xtex = xx;
       ytex = zz;
       xtex2 = xx2;
       ytex2 = zz2;
       xstart = &xbar0;
       ystart = &zbar0;
       break;
     case UP_X:    
       faceptr->normal[0]=(float)1.0;
       if(facetype==VENT_face&&vi!=NULL&&vi->dummy==0)offset[0] = meshi->vent_offset[0];
       faceptr->imin=faceptr->imax;
       if(yminmax2[0]==yminmax2[1]||zminmax2[0]==zminmax2[1])faceptr->thinface=1;
       xtex = yy;
       ytex = zz;
       xtex2 = yy2;
       ytex2 = zz2;
       xstart = &ybar0;
       ystart = &zbar0;
       break;
     case UP_Y:   
       faceptr->normal[1]=(float)1.0;
       if(facetype==VENT_face&&vi!=NULL&&vi->dummy==0)offset[1] = meshi->vent_offset[1];
       faceptr->jmin=faceptr->jmax;
       if(xminmax2[0]==xminmax2[1]||zminmax2[0]==zminmax2[1])faceptr->thinface=1;
       xtex = xx;
       ytex = zz;
       xtex2 = xx2;
       ytex2 = zz2;
       xstart = &xbar0;
       ystart = &zbar0;
       break;
     case DOWN_X:  
       if(facetype==VENT_face&&vi!=NULL&&vi->dummy==0)offset[0] = -meshi->vent_offset[0];
       xtex = yy;
       ytex = zz;
       xtex2 = yy2;
       ytex2 = zz2;
       faceptr->normal[0]=(float)-1.0;
       faceptr->imax=faceptr->imin;
       if(yminmax2[0]==yminmax2[1]||zminmax2[0]==zminmax2[1])faceptr->thinface=1;
       xstart = &ybar0;
       ystart = &zbar0;
       break;
     case DOWN_Z: 
       if(facetype==VENT_face&&vi!=NULL&&vi->dummy==0)offset[2] = -meshi->vent_offset[2];
       xtex = xx;
       ytex = yy;
       xtex2 = xx2;
       ytex2 = yy2;
       faceptr->normal[2]=(float)-1.0;
       faceptr->kmax=faceptr->kmin;
       if(xminmax2[0]==xminmax2[1]||yminmax2[0]==yminmax2[1])faceptr->thinface=1;
       xstart = &xbar0;
       ystart = &ybar0;
       break;
     case UP_Z:   
       if(facetype==VENT_face&&vi!=NULL&&vi->dummy==0)offset[2] = meshi->vent_offset[2];
       xtex = xx;
       ytex = yy;
       xtex2 = xx2;
       ytex2 = yy2;
       faceptr->normal[2]=(float)1.0;
       faceptr->kmin=faceptr->kmax;
       if(xminmax2[0]==xminmax2[1]||yminmax2[0]==yminmax2[1])faceptr->thinface=1;
       xstart = &xbar0;
       ystart = &ybar0;
       break;
     default:
       ASSERT(FFALSE);
       break;
    }
    if(facetype==SHADED_FRAME_face){
      faceptr->normal[0]*=-1;
      faceptr->normal[1]*=-1;
      faceptr->normal[2]*=-1;
    }
    bfi=blockfaceindex + 4*faceptr->dir;
    faceptr->approx_center_coord[0]=0.0;
    faceptr->approx_center_coord[1]=0.0;
    faceptr->approx_center_coord[2]=0.0;
    faceptr->dist2eye=0.0;
    faceptr->xmin=xx2[bfi[0]]+offset[0];
    faceptr->ymin=yy2[bfi[0]]+offset[1];
    faceptr->zmin=zz2[bfi[0]]+offset[2];
    faceptr->xmax=xx2[bfi[0]]+offset[0];
    faceptr->ymax=yy2[bfi[0]]+offset[1];
    faceptr->zmax=zz2[bfi[0]]+offset[2];
    for(k=0;k<4;k++){
      float xvert, yvert, zvert;

      jjj = bfi[k];
      
      xvert=xx[jjj]+offset[0];
      yvert=yy[jjj]+offset[1];
      zvert=zz[jjj]+offset[2];
      faceptr->approx_vertex_coords[3*k]=xvert;
      faceptr->approx_vertex_coords[3*k+1]=yvert;
      faceptr->approx_vertex_coords[3*k+2]=zvert;

      faceptr->approx_center_coord[0]+=xvert;
      faceptr->approx_center_coord[1]+=yvert;
      faceptr->approx_center_coord[2]+=zvert;


      faceptr->exact_vertex_coords[3*k]=xx2[jjj]+offset[0];
      faceptr->exact_vertex_coords[3*k+1]=yy2[jjj]+offset[1];
      faceptr->exact_vertex_coords[3*k+2]=zz2[jjj]+offset[2];
      if(faceptr->exact_vertex_coords[3*k]<faceptr->xmin)faceptr->xmin=faceptr->exact_vertex_coords[3*k];
      if(faceptr->exact_vertex_coords[3*k]>faceptr->xmax)faceptr->xmax=faceptr->exact_vertex_coords[3*k];
      if(faceptr->exact_vertex_coords[3*k+1]<faceptr->ymin)faceptr->ymin=faceptr->exact_vertex_coords[3*k+1];
      if(faceptr->exact_vertex_coords[3*k+1]>faceptr->ymax)faceptr->ymax=faceptr->exact_vertex_coords[3*k+1];
      if(faceptr->exact_vertex_coords[3*k+2]<faceptr->zmin)faceptr->zmin=faceptr->exact_vertex_coords[3*k+2];
      if(faceptr->exact_vertex_coords[3*k+2]>faceptr->zmax)faceptr->zmax=faceptr->exact_vertex_coords[3*k+2];
    }
    faceptr->approx_center_coord[0]/=4.0;
    faceptr->approx_center_coord[1]/=4.0;
    faceptr->approx_center_coord[2]/=4.0;


    {
      float xa_texture[4], ya_texture[4];
      float xe_texture[4], ye_texture[4];
      float dx_e, dy_e, dx_a, dy_a;

      xa_texture[0]=(*xstart+xyzmaxdiff*xtex[bfi[0]]);
      ya_texture[0]=(*ystart+xyzmaxdiff*ytex[bfi[0]]);

      xe_texture[0]=(*xstart+xyzmaxdiff*xtex2[bfi[0]]);
      ye_texture[0]=(*ystart+xyzmaxdiff*ytex2[bfi[0]]);

      switch (faceptr->dir){
        case DOWN_X:
        case UP_X:
          xe_texture[0] -= faceptr->texture_origin[1];
          ye_texture[0] -= faceptr->texture_origin[2];
          xa_texture[0] -= faceptr->texture_origin[1];
          ya_texture[0] -= faceptr->texture_origin[2];
          break;
        case DOWN_Y:
        case UP_Y:
          xe_texture[0] -= faceptr->texture_origin[0];
          ye_texture[0] -= faceptr->texture_origin[2];
          xa_texture[0] -= faceptr->texture_origin[0];
          ya_texture[0] -= faceptr->texture_origin[2];
          break;
        case DOWN_Z:
        case UP_Z:
          xe_texture[0] -= faceptr->texture_origin[0];
          ye_texture[0] -= faceptr->texture_origin[1];
          xa_texture[0] -= faceptr->texture_origin[0];
          ya_texture[0] -= faceptr->texture_origin[1];
          break;
        default:
          ASSERT(FFALSE);
          break;
      }

      dx_a = xyzmaxdiff*(xtex[bfi[0]]-xtex[bfi[1]]);
      dy_a = xyzmaxdiff*(ytex[bfi[0]]-ytex[bfi[2]]);
      dx_e = xyzmaxdiff*(xtex2[bfi[0]]-xtex2[bfi[1]]);
      dy_e = xyzmaxdiff*(ytex2[bfi[0]]-ytex2[bfi[2]]);
      if(dx_a<0.0)dx_a=-dx_a;
      if(dy_a<0.0)dy_a=-dy_a;
      if(dx_e<0.0)dx_e=-dx_e;
      if(dy_e<0.0)dy_e=-dy_e;

      xe_texture[1] = xe_texture[0] + dx_e;
      xe_texture[2] = xe_texture[1];
      xe_texture[3] = xe_texture[0];
      ye_texture[1] = ye_texture[0];
      ye_texture[2] = ye_texture[1] + dy_e;
      ye_texture[3] = ye_texture[2];

      xa_texture[1] = xa_texture[0] + dx_a;
      xa_texture[2] = xa_texture[1];
      xa_texture[3] = xa_texture[0];
      ya_texture[1] = ya_texture[0];
      ya_texture[2] = ya_texture[1] + dy_a;
      ya_texture[3] = ya_texture[2];

      for(k=0;k<4;k++){
        faceptr->approx_texture_coords[2*k]=xa_texture[k]/t_width; 
        faceptr->approx_texture_coords[2*k+1]=ya_texture[k]/t_height;
        faceptr->exact_texture_coords[2*k]=xe_texture[k]/t_width; 
        faceptr->exact_texture_coords[2*k+1]=ye_texture[k]/t_height;
      }
    }
    faceptr++;
  }
}

/* ------------------ clip_face ------------------------ */

int clip_face(facedata *facei){
  if(xyz_clipplane==0)return 0;
  if(clip_x==1&&xbar0+xyzmaxdiff*facei->xmax<clip_x_val)return 1;
  if(clip_X==1&&xbar0+xyzmaxdiff*facei->xmin>clip_X_val)return 1;
  if(clip_y==1&&ybar0+xyzmaxdiff*facei->ymax<clip_y_val)return 1;
  if(clip_Y==1&&ybar0+xyzmaxdiff*facei->ymin>clip_Y_val)return 1;
  if(clip_z==1&&zbar0+xyzmaxdiff*facei->zmax<clip_z_val)return 1;
  if(clip_Z==1&&zbar0+xyzmaxdiff*facei->zmin>clip_Z_val)return 1;
  return 0;
}

/* ------------------ set_cull_vis ------------------------ */

void set_cull_vis(void){
  int imesh;
#ifdef pp_GEOMPRINT
  int nports=0;
  int ntotal=0;
#endif

  if(update_initcullgeom==1){
    initcullgeom(cullgeom);
    update_facelists();
  }
  for(imesh=0;imesh<nmeshes;imesh++){
    int iport;
    mesh *meshi;

    meshi = meshinfo + imesh;
#ifdef pp_GEOMPRINT
    ntotal+=meshi->ncullgeominfo;
#endif    
    for(iport=0;iport<meshi->ncullgeominfo;iport++){
      culldata *culli;
      float xx[2], yy[2], zz[2];

      culli = meshi->cullgeominfo+iport;
      culli->vis=0;

      xx[0] = (culli->xbeg-xbar0)/xyzmaxdiff;
      xx[1] = (culli->xend-xbar0)/xyzmaxdiff;
      yy[0] = (culli->ybeg-ybar0)/xyzmaxdiff;
      yy[1] = (culli->yend-ybar0)/xyzmaxdiff;
      zz[0] = (culli->zbeg-zbar0)/xyzmaxdiff;
      zz[1] = (culli->zend-zbar0)/xyzmaxdiff;
      
      if(PointInFrustum(xx[0],yy[0],zz[0])==1){
        culli->vis=1;
#ifdef pp_GEOMPRINT
        nports++;
#endif        
        continue;
      }
      if(PointInFrustum(xx[1],yy[0],zz[0])==1){
        culli->vis=1;
#ifdef pp_GEOMPRINT
        nports++;
#endif        
        continue;
      }
      if(PointInFrustum(xx[0],yy[1],zz[0])==1){
        culli->vis=1;
#ifdef pp_GEOMPRINT
        nports++;
#endif        
        continue;
      }
      if(PointInFrustum(xx[1],yy[1],zz[0])==1){
        culli->vis=1;
#ifdef pp_GEOMPRINT
        nports++;
#endif        
        continue;
      }
      if(PointInFrustum(xx[0],yy[0],zz[1])==1){
        culli->vis=1;
#ifdef pp_GEOMPRINT
        nports++;
#endif        
        continue;
      }
      if(PointInFrustum(xx[1],yy[0],zz[1])==1){
        culli->vis=1;
#ifdef pp_GEOMPRINT
        nports++;
#endif        
        continue;
      }
      if(PointInFrustum(xx[0],yy[1],zz[1])==1){
        culli->vis=1;
#ifdef pp_GEOMPRINT
        nports++;
#endif        
        continue;
      }
      if(PointInFrustum(xx[1],yy[1],zz[1])==1){
        culli->vis=1;
#ifdef pp_GEOMPRINT
        nports++;
#endif        
        continue;
      }
    }
  }
#ifdef pp_GEOMPRINT
  printf("ports=%i ports visible=%i\n",ntotal,nports); 
#endif
}


/* ------------------ comparesinglefaces ------------------------ */

int comparesinglefaces0( const void *arg1, const void *arg2 ){
  facedata *facei, *facej;
  int dirs[6];
  facei = *(facedata **)arg1;
  facej = *(facedata **)arg2;

  dirs[DOWN_X]=1;
  dirs[UP_X]=1;
  dirs[DOWN_Y]=2;
  dirs[UP_Y]=2;
  dirs[DOWN_Z]=3;
  dirs[UP_Z]=3;
  if(dirs[facei->dir]<dirs[facej->dir])return -1;
  if(dirs[facei->dir]>dirs[facej->dir])return 1;
  switch(facei->dir){
    case DOWN_X:
    case UP_X:
      if(facei->imin<facej->imin)return -1;
      if(facei->imin>facej->imin)return 1;
      if(facei->jmin<facej->jmin)return -1;
      if(facei->jmin>facej->jmin)return 1;
      if(facei->kmin<facej->kmin)return -1;
      if(facei->kmin>facej->kmin)return 1;
      if(facei->imax<facej->imax)return -1;
      if(facei->imax>facej->imax)return 1;
      if(facei->jmax<facej->jmax)return -1;
      if(facei->jmax>facej->jmax)return 1;
      if(facei->kmax<facej->kmax)return -1;
      if(facei->kmax>facej->kmax)return 1;
      break;
    case DOWN_Y:
    case UP_Y:
      if(facei->jmin<facej->jmin)return -1;
      if(facei->jmin>facej->jmin)return 1;
      if(facei->imin<facej->imin)return -1;
      if(facei->imin>facej->imin)return 1;
      if(facei->kmin<facej->kmin)return -1;
      if(facei->kmin>facej->kmin)return 1;
      if(facei->jmax<facej->jmax)return -1;
      if(facei->jmax>facej->jmax)return 1;
      if(facei->imax<facej->imax)return -1;
      if(facei->imax>facej->imax)return 1;
      if(facei->kmax<facej->kmax)return -1;
      if(facei->kmax>facej->kmax)return 1;
      break;
    case DOWN_Z:
    case UP_Z:
      if(facei->kmin<facej->kmin)return -1;
      if(facei->kmin>facej->kmin)return 1;
      if(facei->imin<facej->imin)return -1;
      if(facei->imin>facej->imin)return 1;
      if(facei->jmin<facej->jmin)return -1;
      if(facei->jmin>facej->jmin)return 1;
      if(facei->kmax<facej->kmax)return -1;
      if(facei->kmax>facej->kmax)return 1;
      if(facei->imax<facej->imax)return -1;
      if(facei->imax>facej->imax)return 1;
      if(facei->jmax<facej->jmax)return -1;
      if(facei->jmax>facej->jmax)return 1;
      break;
    default:
      ASSERT(0);
      break;
  }
  if(facei->color<facej->color)return  1;
  if(facei->color>facej->color)return -1;
  return 0;
}

/* ------------------ comparesinglefaces ------------------------ */

int comparesinglefaces( const void *arg1, const void *arg2 ){
  facedata *facei, *facej;

  facei = *(facedata **)arg1;
  facej = *(facedata **)arg2;

  if(facei->dir<facej->dir)return -1;
  if(facei->dir>facej->dir)return 1;
  switch(facei->dir){
    case DOWN_X:   // sort DOWN data from big to small
      if(facei->imin<facej->imin)return 1;
      if(facei->imin>facej->imin)return -1;
      break;
    case UP_X:     // sort UP data from small to big
      if(facei->imin<facej->imin)return -1;
      if(facei->imin>facej->imin)return 1;
      break;
    case DOWN_Y:
      if(facei->jmin<facej->jmin)return 1;
      if(facei->jmin>facej->jmin)return -1;
      break;
    case UP_Y:
      if(facei->jmin<facej->jmin)return -1;
      if(facei->jmin>facej->jmin)return 1;
      break;
    case DOWN_Z:
      if(facei->kmin<facej->kmin)return 1;
      if(facei->kmin>facej->kmin)return -1;
      break;
    case UP_Z:
      if(facei->kmin<facej->kmin)return -1;
      if(facei->kmin>facej->kmin)return 1;
      break;
    default:
      ASSERT(0);
      break;
  }
  if(facei->color<facej->color)return  1;
  if(facei->color>facej->color)return -1;
  return 0;
}

/* ------------------ comparecolorfaces ------------------------ */

int comparecolorfaces( const void *arg1, const void *arg2 ){
  facedata *facei, *facej;

  facei = *(facedata **)arg1;
  facej = *(facedata **)arg2;

  if(facei->color<facej->color)return  1;
  if(facei->color>facej->color)return -1;
  return 0;
}

/* ------------------ update_facelists ------------------------ */

void update_facelists(void){
  int n_textures, n_outlines;
  int n_normals_single, n_normals_double, n_transparent_double;
  int i,j,k;
  mesh *meshi;
  facedata *facej;
  int local_showpatch, loadpatch;
  patch *patchi;
  int patchfilenum;
  int patch_dir[6]={2,1,3,0,4,5};
  int vent_offset, outline_offset, exteriorsurface_offset;
  ventdata *vi;
  int drawing_smooth;
  int drawing_transparent, drawing_blockage_transparent, drawing_vent_transparent;

//     visBlocks (visBLOCKNormal, visBLOCKAsInput )
//     visSmoothAsNormal
//     visTransparentBlockage

  get_drawing_parms(&drawing_smooth, &drawing_transparent, &drawing_blockage_transparent, &drawing_vent_transparent);

  if(updatehiddenfaces==1)UpdateHiddenFaces();
  updatefacelists=0;
  nface_normals_single=0;
  nface_normals_double=0;
  nface_transparent_double=0;
  nface_textures=0;
  nface_outlines=0;
  nface_transparent=0;
  if(opengldefined==1){
    glutPostRedisplay();
  }
  for(i=0;i<nmeshes;i++){
    meshi = meshinfo + i;

    for(j=0;j<meshi->nfaces;j++){
      facej = meshi->faceinfo + j;
      facej->patchpresent=0;
      facej->cullport=NULL;
    }

    local_showpatch=0;
    loadpatch=0;
    patchfilenum=meshi->patchfilenum;
    if(hidepatchsurface==1&&patchfilenum>=0&&patchfilenum<npatchinfo){
      patchi = patchinfo + patchfilenum;
      if(patchi->loaded==1)loadpatch=1;
      if(patchi->display==1)local_showpatch=1;
    }
    else{
      patchi=NULL;
    }

    if(
      local_showpatch==1&&loadpatch==1){
      for(j=0;j<meshi->nbptrs;j++){
        blockagedata *bc;

        bc=meshi->blockageinfoptrs[j];
        if(bc->prop!=NULL&&bc->prop->blockvis==0)continue;
        facej = meshi->faceinfo + 6*j;
        for(k=0;k<6;k++){
          facej->patchpresent=1-bc->patchvis[patch_dir[k]];
          facej++;
        }
      }
    }

    n_normals_single=0;
    n_normals_double=0;
    n_transparent_double=0;
    n_textures=0;
    n_outlines=0;

    for(j=0;j<meshi->nfaces;j++){
      facej = meshi->faceinfo + j;

      if(showonly_hiddenfaces==0&&facej->hidden==1)continue;
      if(showonly_hiddenfaces==1&&facej->hidden==0)continue;
      if(facej->thinface==1)continue;
      if(facej->bc!=NULL&&facej->bc->prop!=NULL&&facej->bc->prop->blockvis==0)continue;
      if(clip_face(facej)==1)continue;

      vent_offset = 6*meshi->nbptrs;
      outline_offset = vent_offset + meshi->nvents;
      exteriorsurface_offset = outline_offset + 6;
      if(showedit_dialog==1&&j<vent_offset){
        if(facej->show_bothsides==0)meshi->face_normals_single[n_normals_single++]=facej;
        if(facej->show_bothsides==1)meshi->face_normals_double[n_normals_double++]=facej;
        continue;
      }
      if(j<vent_offset){
        if(visBlocks==visBLOCKHide)continue;
      }
      if(j>=outline_offset&&j<outline_offset+6&&visFrame==0){
        continue;
      }
      if(j>=vent_offset&&j<vent_offset+meshi->nvents){
        vi = meshi->ventinfo+j-vent_offset;
        if(visOpenVents==0&&vi->isOpenvent==1)continue;
        if(visDummyVents==0&&vi->dummy==1)continue;
        if(visOtherVents==0&&vi->isOpenvent==0&&vi->dummy==0)continue;
        if(patchi!=NULL&&patchi->loaded==1&&patchi->display==1&&
          (vis_threshold==0||vis_onlythreshold==0||do_threshold==0)&&
          (vi->dummy==1||vi->hideboundary==0)){
          continue;
        }
        if(facej->transparent==1&&drawing_vent_transparent==1){
          if(facej->show_bothsides==1){
            meshi->face_transparent_double[n_transparent_double++]=facej;
          }
          else{
            face_transparent[nface_transparent++]=facej;
          }
          continue;
        }
        if(facej->textureinfo!=NULL&&facej->textureinfo->display==1){
          meshi->face_textures[n_textures++]=facej;
          continue;
        }
      }
      if(j>=exteriorsurface_offset){
        switch (j-exteriorsurface_offset){
         case DOWN_Z:
           if(visFloor==0){
             continue;
           }
          break;
         case UP_Z:
           if(visCeiling==0)continue;
          break;
         case UP_X:
         case DOWN_X:
         case UP_Y:
         case DOWN_Y:
           if(visWalls==0)continue;
          break;
         default:
           ASSERT(FFALSE);
           break;
        }
      }
      if((
         (visBlocks==visBLOCKOutline||visBlocks==visBLOCKSolidOutline)&&j<vent_offset)||
         (facej->patchpresent==1&&(vis_threshold==0||vis_onlythreshold==0||do_threshold==0))||
         (facej->type==BLOCK_outline&&visBlocks==visBLOCKAsInput)||
         ((j>=vent_offset&&j<vent_offset+meshi->nvents)&&vi->isOpenvent==1&&visOpenVentsAsOutline==1)
        ){
        meshi->face_outlines[n_outlines++]=facej;
        if(visBlocks!=visBLOCKSolidOutline)continue;
      }
      if(j<vent_offset){
        int drawing_texture=0;

        if(facej->type==BLOCK_texture&&facej->textureinfo!=NULL&&facej->textureinfo->display==1){
          drawing_texture=1;
        }

        if(drawing_smooth==1&&facej->type==BLOCK_smooth)continue;
        if(facej->transparent==0||drawing_blockage_transparent==0){
          if(drawing_texture==0){
            if(facej->show_bothsides==0){
              meshi->face_normals_single[n_normals_single++]=facej;
            }
            if(facej->show_bothsides==1){
              meshi->face_normals_double[n_normals_double++]=facej;
            }
            continue;
          }
        }
        if(facej->transparent==1&&drawing_blockage_transparent==1){
          face_transparent[nface_transparent++]=facej;
          continue;
        }
      }

      switch (facej->type){
       case BLOCK_regular:
        if(facej->show_bothsides==0)meshi->face_normals_single[n_normals_single++]=facej;
        if(facej->show_bothsides==1)meshi->face_normals_double[n_normals_double++]=facej;
        break;
       case BLOCK_texture:
        if(facej->textureinfo!=NULL){
          if(facej->textureinfo->display==1){
            meshi->face_textures[n_textures++]=facej;
          }
          else{
            if(facej->type2==BLOCK_face){
              if(facej->show_bothsides==0)meshi->face_normals_single[n_normals_single++]=facej;
              if(facej->show_bothsides==1)meshi->face_normals_double[n_normals_double++]=facej;
            }
            if(facej->type2==VENT_face)meshi->face_outlines[n_outlines++]=facej;
          }
          continue;
        }
        break;
       case BLOCK_outline:
        meshi->face_outlines[n_outlines++]=facej;
         break;
       case BLOCK_smooth:
         if(updatesmoothblocks!=0||visSmoothAsNormal==1){
           if(facej->show_bothsides==0)meshi->face_normals_single[n_normals_single++]=facej;
           if(facej->show_bothsides==1)meshi->face_normals_double[n_normals_double++]=facej;
         }
         break;
       case BLOCK_hidden:
         break;
       default:
         printf("facej->type=%i\n",facej->type);
         ASSERT(FFALSE);
         break;
      }
    }

    meshi->nface_textures = n_textures;
    meshi->nface_normals_single  = n_normals_single;
    meshi->nface_normals_double  = n_normals_double;
    meshi->nface_transparent_double  = n_transparent_double;
    meshi->nface_outlines = n_outlines;
    nface_textures += n_textures;
    nface_normals_single += n_normals_single;
    nface_normals_double += n_normals_double;
    nface_transparent_double += n_transparent_double;
    nface_outlines += n_outlines;

    if(use_new_drawface==0)continue;

    meshi->nface_normals_single_DOWN_X=0;
    meshi->nface_normals_single_UP_X=0;
    meshi->nface_normals_single_DOWN_Y=0;
    meshi->nface_normals_single_UP_Y=0;
    meshi->nface_normals_single_DOWN_Z=0;
    meshi->nface_normals_single_UP_Z=0;
    if(n_normals_single>1){
      int iface;
      int istartD=-1,istartU=-1;
      int jstartD=-1,jstartU=-1;
      int kstartD=-1,kstartU=-1;
      int nhidden;

      nhidden=0;
      qsort((facedata **)meshi->face_normals_single,(size_t)n_normals_single,sizeof(facedata *),comparesinglefaces0);
      meshi->face_normals_single[0]->dup=0;
      for(iface=1;iface<meshi->nface_normals_single;iface++){
        facedata *facei;
        facedata *faceim1;

        facei=meshi->face_normals_single[iface];
        facei->dup=0;
        faceim1 = meshi->face_normals_single[iface-1];
        if(
          facei->imax-facei->imin<=1&&facei->jmax-facei->jmin<=1&&facei->kmax-facei->kmin<=1&& // only hide duplicate one cell faces
          faceim1->imin==facei->imin&&faceim1->imax==facei->imax&&
          faceim1->jmin==facei->jmin&&faceim1->jmax==facei->jmax&&
          faceim1->kmin==facei->kmin&&faceim1->kmax==facei->kmax&&faceim1->dir!=facei->dir&&facei->thinface==0){
          if(*(faceim1->showtimelist_handle)==NULL)faceim1->dup=1;
          if(*(facei->showtimelist_handle)==NULL){
            facei->dup=1;
            nhidden++;
          }
        }
      }
      if(nhidden>0){
        n_normals_single=0;
        for(iface=0;iface<meshi->nface_normals_single;iface++){
          facedata *facei  ;

          facei=meshi->face_normals_single[iface];
          if(facei->dup==0)meshi->face_normals_single[n_normals_single++]=facei;
        }
        meshi->nface_normals_single=n_normals_single;
      }
#ifdef pp_GEOMPRINT
      printf("faces removed=%i\n",nhidden);
#endif      

      qsort((facedata **)meshi->face_normals_single,(size_t)n_normals_single,sizeof(facedata *),comparesinglefaces);
      for(iface=0;iface<meshi->nface_normals_single;iface++){
        facedata *facei;

        facei=meshi->face_normals_single[iface];
        facei->cullport=get_face_port(meshi,facei);
        switch(facei->dir){
          case DOWN_X:
            if(istartD==-1){
              meshi->face_normals_single_DOWN_X=meshi->face_normals_single+iface;
              istartD=0;
            }
            meshi->nface_normals_single_DOWN_X++;
            break;
          case UP_X:
            if(istartU==-1){
              meshi->face_normals_single_UP_X=meshi->face_normals_single+iface;
              istartU=0;
            }
            meshi->nface_normals_single_UP_X++;
            break;
          case DOWN_Y:
            if(jstartD==-1){
              meshi->face_normals_single_DOWN_Y=meshi->face_normals_single+iface;
              jstartD=0;
            }
            meshi->nface_normals_single_DOWN_Y++;
            break;
          case UP_Y:
            if(jstartU==-1){
              meshi->face_normals_single_UP_Y=meshi->face_normals_single+iface;
              jstartU=0;
            }
            meshi->nface_normals_single_UP_Y++;
            break;
          case DOWN_Z:
            if(kstartD==-1){
              meshi->face_normals_single_DOWN_Z=meshi->face_normals_single+iface;
              kstartD=0;
            }
            meshi->nface_normals_single_DOWN_Z++;
            break;
          case UP_Z:
            if(kstartU==-1){
              meshi->face_normals_single_UP_Z=meshi->face_normals_single+iface;
              kstartU=0;
            }
            meshi->nface_normals_single_UP_Z++;
            break;
          default:
            ASSERT(0);
            break;
        }
      }
    }
    if(n_normals_double>1){
      qsort((facedata **)meshi->face_normals_double,(size_t)n_normals_double,sizeof(facedata *),comparecolorfaces);
    }
    if(n_outlines>1){
      qsort((facedata **)meshi->face_outlines,(size_t)n_outlines,sizeof(facedata *),comparecolorfaces);
    }
  }
}

/* ------------------ drawselect_faces ------------------------ */

void drawselect_faces(){
  int i,j,k;
  mesh *meshi;
  float *vertices;
  int color_index=0;
  unsigned char r, g, b;
  int sides[]={DOWN_Y,UP_X,UP_Y,DOWN_X,DOWN_Z,UP_Z};
  facedata *facek;

  glDisable(GL_LIGHTING);
  glBegin(GL_QUADS);
  for(i=0;i<nmeshes;i++){
    meshi=meshinfo + i;
    for(j=0;j<meshi->nbptrs;j++){
      for(k=0;k<6;k++){
        facek = meshi->faceinfo + 6*j + sides[k];

        vertices = facek->approx_vertex_coords;
        getrgb(color_index+1,&r,&g,&b);
        glColor3ub(r,g,b);
        color_index++;
        glVertex3fv(vertices);
        glVertex3fv(vertices+3);
        glVertex3fv(vertices+6);
        glVertex3fv(vertices+9);
      }
    }
  }
  glEnd();
  return;
}

#ifdef pp_GEOMPRINT
#define COLOR_SWAPS color_swaps++;
#define FACES_DRAWN faces_drawn++;
#else
#define COLOR_SWAPS
#define FACES_DRAWN
#endif
#define DRAWFACE(DEFfacetest,DEFeditcolor)    \
        float *facepos;\
        culldata *cullport;\
        facei = face_START[i];\
        cullport=facei->cullport;\
        if(cullport!=NULL&&cullport->vis==0)continue;\
        if(blocklocation==BLOCKlocation_grid){\
          vertices = facei->approx_vertex_coords;\
        }\
        else{\
          vertices = facei->exact_vertex_coords;\
        }\
        facepos=vertices;\
        if(DEFfacetest)break;\
 \
        showtimelist_handle = facei->showtimelist_handle;\
        showtimelist = *showtimelist_handle;\
        if(showtimelist!=NULL&&showtimelist[itimes]==0)continue;\
        if(showedit_dialog==0){\
          new_color=facei->color;\
        }\
        else{\
          if(visNormalEditColors==0)new_color=block_ambient2;\
          if(visNormalEditColors==1)new_color=facei->color;\
          if(highlight_block==facei->blockageindex&&highlight_mesh==facei->meshindex){\
            new_color=DEFeditcolor;\
          }\
        }\
        if(new_color!=old_color){\
          old_color=new_color;\
          glColor4fv(old_color);\
          COLOR_SWAPS\
        }\
        glVertex3fv(vertices);\
        glVertex3fv(vertices+3);\
        glVertex3fv(vertices+6);\
        glVertex3fv(vertices);\
        glVertex3fv(vertices+6);\
        glVertex3fv(vertices+9);\
        FACES_DRAWN

/* ------------------ drawfaces ------------------------ */

void draw_faces(){
  int i,j;
  mesh *meshi;
  facedata *facei;
  float *vertices,*tvertices;
  texture *texti;
  float *new_color,*old_color=NULL;
  int **showtimelist_handle, *showtimelist;
  float up_color[4]={0.9,0.9,0.9,1.0};
  float down_color[4]={0.1,0.1,0.1,1.0};
  float highlight_color[4]={1.0,0.0,0.0,1.0};
#ifdef pp_GEOMPRINT
  int color_swaps=0;
  int faces_drawn=0;
#endif

  if(nface_normals_single>0){
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glEnable(GL_COLOR_MATERIAL);
    glBegin(GL_TRIANGLES);
    for(j=0;j<nmeshes;j++){
      facedata **face_START;
      meshi=meshinfo + j;
      if(meshi->blockvis==0)continue;

      // DOWN_X faces

      glNormal3f(-1.0,0.0,0.0);
      face_START=meshi->face_normals_single_DOWN_X;
      for(i=0;i<meshi->nface_normals_single_DOWN_X;i++){
        DRAWFACE(scaled_eyepos[0]>facepos[0],down_color)
      }

      // UP_X faces

      glNormal3f(1.0,0.0,0.0);
      face_START=meshi->face_normals_single_UP_X;
      for(i=0;i<meshi->nface_normals_single_UP_X;i++){
        DRAWFACE(scaled_eyepos[0]<facepos[0],up_color)
      }

      // DOWN_Y faces

      glNormal3f(0.0,-1.0,0.0);
      face_START=meshi->face_normals_single_DOWN_Y;
      for(i=0;i<meshi->nface_normals_single_DOWN_Y;i++){
        DRAWFACE(scaled_eyepos[1]>facepos[1],down_color)
      }

      // UP_Y faces

      glNormal3f(0.0,1.0,0.0);
      face_START=meshi->face_normals_single_UP_Y;
      for(i=0;i<meshi->nface_normals_single_UP_Y;i++){
        DRAWFACE(scaled_eyepos[1]<facepos[1],up_color)
      }

      // DOWN_Z faces

      glNormal3f(0.0,0.0,-1.0);
      face_START=meshi->face_normals_single_DOWN_Z;
      for(i=0;i<meshi->nface_normals_single_DOWN_Z;i++){
        DRAWFACE(scaled_eyepos[2]>facepos[2],down_color)
      }

      // UP_Z faces

      glNormal3f(0.0,0.0,1.0);
      face_START=meshi->face_normals_single_UP_Z;
      for(i=0;i<meshi->nface_normals_single_UP_Z;i++){
        DRAWFACE(scaled_eyepos[2]<facepos[2],up_color)
      }
    }
    glEnd();
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
  }
#ifdef pp_GEOMPRINT
  printf("faces=%i, faces drawn=%i, color switches=%i\n",nface_normals_single,faces_drawn,color_swaps);
#endif  
  if(nface_normals_double>0){
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glEnable(GL_COLOR_MATERIAL);
    if(cullfaces==1)glDisable(GL_CULL_FACE);
    glBegin(GL_QUADS);
    for(j=0;j<nmeshes;j++){
      meshi=meshinfo + j;
      for(i=0;i<meshi->nface_normals_double;i++){
        facei = meshi->face_normals_double[i];
        if(blocklocation==BLOCKlocation_grid){
          vertices = facei->approx_vertex_coords;
        }
        else{
          vertices = facei->exact_vertex_coords;
        }
        showtimelist_handle = facei->showtimelist_handle;
        showtimelist = *showtimelist_handle;
        if(showtimelist!=NULL&&showtimelist[itimes]==0)continue;
        if(showedit_dialog==0){
          new_color=facei->color;
        }
        else{
          if(visNormalEditColors==0)new_color=block_ambient2;
          if(visNormalEditColors==1)new_color=facei->color;
          if(highlight_block==facei->blockageindex&&highlight_mesh==facei->meshindex){
            new_color=highlight_color;
            switch (xyz_dir){
             case XDIR:
              if(facei->dir==UP_X)new_color=up_color;
              if(facei->dir==DOWN_X)new_color=down_color;
              break;
             case YDIR:
              if(facei->dir==UP_Y)new_color=up_color;
              if(facei->dir==DOWN_Y)new_color=down_color;
              break;
             case ZDIR:
              if(facei->dir==UP_Z)new_color=up_color;
              if(facei->dir==DOWN_Z)new_color=down_color;
              break;
             default:
              ASSERT(FFALSE);
              break;
            }
          }
        }
        if(new_color!=old_color){
          old_color=new_color;
          glColor4fv(old_color);
        }
        glNormal3fv(facei->normal);
        glVertex3fv(vertices);
        glVertex3fv(vertices+3);
        glVertex3fv(vertices+6);
        glVertex3fv(vertices+9);
      }
    }
    glEnd();
    if(cullfaces==1)glEnable(GL_CULL_FACE);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
  }
  if(nface_outlines>0){
    glDisable(GL_LIGHTING);
    antialias(1);
    glLineWidth(linewidth);
    glBegin(GL_LINES);
    for(j=0;j<nmeshes;j++){
      meshi = meshinfo + j;
      if(meshi->blockvis==0)continue;
      for(i=0;i<meshi->nface_outlines;i++){
        facei = meshi->face_outlines[i];
        showtimelist_handle = facei->showtimelist_handle;
        showtimelist = *showtimelist_handle;
        if(showtimelist!=NULL&&showtimelist[itimes]==0&&facei->type2==BLOCK_face)continue;
        if(blocklocation==BLOCKlocation_grid){
          vertices = facei->approx_vertex_coords;
        }
        else{
          vertices = facei->exact_vertex_coords;
        }
        if(facei->type2!=OUTLINE_FRAME_face||highlight_flag==1){
          glEnd();
          if(nmeshes>1&&facei->type2==OUTLINE_FRAME_face&&
            highlight_mesh==facei->meshindex&&highlight_flag==1){
            glLineWidth(highlight_linewidth);
          }
          else{
            glLineWidth(*facei->linewidth);
          }
          glBegin(GL_LINES);
          glColor3fv(facei->linecolor);
          glVertex3fv(vertices);
          glVertex3fv(vertices+3);
          glVertex3fv(vertices+3);
          glVertex3fv(vertices+6);
          glVertex3fv(vertices+6);
          glVertex3fv(vertices+9);
          glVertex3fv(vertices+9);
          glVertex3fv(vertices);
          if(showtimelist!=NULL&&showtimelist[itimes]==0){
            glVertex3fv(vertices);
            glVertex3fv(vertices+6);
            glVertex3fv(vertices+3);
            glVertex3fv(vertices+9);
          }
        }
      }
    }
    glEnd();
    antialias(0);
  }
  if(nface_textures>0){
    glEnable(GL_LIGHTING);
    glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
    glEnable(GL_TEXTURE_2D);
    for(j=0;j<nmeshes;j++){
      meshi = meshinfo + j;
      if(meshi->blockvis==0)continue;
      for(i=0;i<meshi->nface_textures;i++){
        facei=meshi->face_textures[i];
        showtimelist_handle = facei->showtimelist_handle;
        showtimelist = *showtimelist_handle;
        if(showtimelist!=NULL&&showtimelist[itimes]==0)continue;
        texti=facei->textureinfo;
        if(blocklocation==BLOCKlocation_grid){
           vertices = facei->approx_vertex_coords;
          tvertices = facei->approx_texture_coords;
        }
        else{
           vertices = facei->exact_vertex_coords;
          tvertices = facei->exact_texture_coords;
        }

        if(facei->type2==BLOCK_face&&cullfaces==0)glDisable(GL_CULL_FACE);


        glBindTexture(GL_TEXTURE_2D,texti->name);
        glBegin(GL_QUADS);

        glNormal3fv(facei->normal);
        glTexCoord2fv(tvertices);
        glVertex3fv(vertices);

        glTexCoord2fv(tvertices+2);
        glVertex3fv(vertices+3);

        glTexCoord2fv(tvertices+4);
        glVertex3fv(vertices+6);

        glTexCoord2fv(tvertices+6);
        glVertex3fv(vertices+9);
        glEnd();
      }
      if(cullfaces==1)glEnable(GL_CULL_FACE);


    }
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);
  }
}

/* ------------------ compareisonodes ------------------------ */

int comparetransparentfaces( const void *arg1, const void *arg2 ){
  facedata *facei, *facej;

  facei = *(facedata **)arg1;
  facej = *(facedata **)arg2;

  if(facei->dist2eye<facej->dist2eye)return  1;
  if(facei->dist2eye>facej->dist2eye)return -1;
  return 0;
}

/* ------------------ Sort_Transparent_Faces ------------------------ */

void Sort_Transparent_Faces(float *mm){
  int i;
  float *xyzface;
  float xyzeye[3];

  for(i=0;i<nface_transparent;i++){
    facedata *facei;

    facei = face_transparent[i];
    xyzface = facei->approx_center_coord;
    xyzeye[0] = mm[0]*xyzface[0] + mm[4]*xyzface[1] +  mm[8]*xyzface[2] + mm[12];
    xyzeye[1] = mm[1]*xyzface[0] + mm[5]*xyzface[1] +  mm[9]*xyzface[2] + mm[13];
    xyzeye[2] = mm[2]*xyzface[0] + mm[6]*xyzface[1] + mm[10]*xyzface[2] + mm[14];
    xyzeye[0]/=mscale[0];
    xyzeye[1]/=mscale[1];
    xyzeye[2]/=mscale[2];
    facei->dist2eye=xyzeye[0]*xyzeye[0]+xyzeye[1]*xyzeye[1]+xyzeye[2]*xyzeye[2];
  }
  qsort((facedata **)face_transparent,(size_t)nface_transparent,sizeof(facedata *),comparetransparentfaces);


}

/* ------------------ draw_transparent_faces ------------------------ */

void draw_transparent_faces(){
  int i,j;
  facedata *facei;
  mesh *meshi;
  float *vertices;
  float old_color[4]={(float)-1.0,(float)-1.0,(float)-1.0,(float)-1.0};
  float *new_color;
  int **showtimelist_handle, *showtimelist;
  float up_color[4]={0.9,0.9,0.9,1.0};
  float down_color[4]={0.1,0.1,0.1,1.0};
  float highlight_color[4]={1.0,0.0,0.0,1.0};
  int drawing_smooth, drawing_transparent, drawing_blockage_transparent, drawing_vent_transparent;

  if(blocklocation==BLOCKlocation_cad||(ncadgeom!=0&&show_cad_and_grid==1))return;

  get_drawing_parms(&drawing_smooth, &drawing_transparent, &drawing_blockage_transparent, &drawing_vent_transparent);

  if(nface_transparent<=0&&nface_transparent_double<=0)return;

  if(drawing_transparent==1)transparenton();

  if(nface_transparent>0){
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glEnable(GL_COLOR_MATERIAL);
    glBegin(GL_QUADS);
    for(i=0;i<nface_transparent;i++){
      facei = face_transparent[i];
      if(blocklocation==BLOCKlocation_grid){
        vertices = facei->approx_vertex_coords;
      }
      else{
        vertices = facei->exact_vertex_coords;
      }
      showtimelist_handle = facei->showtimelist_handle;
      showtimelist = *showtimelist_handle;
      if(showtimelist!=NULL&&showtimelist[itimes]==0)continue;
      if(showedit_dialog==0){
        new_color=facei->color;
      }
      else{
        if(visNormalEditColors==0)new_color=block_ambient2;
        if(visNormalEditColors==1)new_color=facei->color;
        if(highlight_block==facei->blockageindex&&highlight_mesh==facei->meshindex){
          new_color=highlight_color;
          switch (xyz_dir){
           case XDIR:
            if(facei->dir==UP_X)new_color=up_color;
            if(facei->dir==DOWN_X)new_color=down_color;
            break;
           case YDIR:
            if(facei->dir==UP_Y)new_color=up_color;
            if(facei->dir==DOWN_Y)new_color=down_color;
            break;
           case ZDIR:
            if(facei->dir==UP_Z)new_color=up_color;
            if(facei->dir==DOWN_Z)new_color=down_color;
            break;
           default:
            ASSERT(FFALSE);
            break;
          }
        }
      }
      if(
         ABS(new_color[0]-old_color[0])>0.0001||
         ABS(new_color[1]-old_color[1])>0.0001||
         ABS(new_color[2]-old_color[2])>0.0001||
         ABS(new_color[3]-old_color[3])>0.0001||
         use_transparency_geom==1
         ){
        old_color[0]=new_color[0];
        old_color[1]=new_color[1];
        old_color[2]=new_color[2];
        old_color[3]=new_color[3];
        if(use_transparency_geom==1)old_color[3]=transparency_geom;
        glColor4fv(old_color);
      }
      glNormal3fv(facei->normal);
      glVertex3fv(vertices);
      glVertex3fv(vertices+3);
      glVertex3fv(vertices+6);
      glVertex3fv(vertices+9);
    }
    glEnd();
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
  }

  if(nface_transparent_double>0){
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glEnable(GL_COLOR_MATERIAL);
    if(cullfaces==1)glDisable(GL_CULL_FACE);
    glBegin(GL_QUADS);
    for(j=0;j<nmeshes;j++){
      meshi=meshinfo + j;
      for(i=0;i<meshi->nface_transparent_double;i++){
        facei = meshi->face_transparent_double[i];
        if(blocklocation==BLOCKlocation_grid){
          vertices = facei->approx_vertex_coords;
        }
        else{
          vertices = facei->exact_vertex_coords;
        }
        showtimelist_handle = facei->showtimelist_handle;
        showtimelist = *showtimelist_handle;
        if(showtimelist!=NULL&&showtimelist[itimes]==0)continue;
        new_color=facei->color;
        if(
         ABS(new_color[0]-old_color[0])>0.0001||
         ABS(new_color[1]-old_color[1])>0.0001||
         ABS(new_color[2]-old_color[2])>0.0001||
         ABS(new_color[3]-old_color[3])>0.0001
         ){
          old_color[0]=new_color[0];
          old_color[1]=new_color[1];
          old_color[2]=new_color[2];
          old_color[3]=new_color[3];
          glColor4fv(old_color);
        }
        glNormal3fv(facei->normal);
        glVertex3fv(vertices);
        glVertex3fv(vertices+3);
        glVertex3fv(vertices+6);
        glVertex3fv(vertices+9);
      }
    }
    glEnd();
    if(cullfaces==1)glEnable(GL_CULL_FACE);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
  }

  if(drawing_transparent==1)transparentoff();
}

/* ------------------ UpdateHiddenFaces ------------------------ */

void UpdateHiddenFaces(){
  int i;


  updatehiddenfaces=0;
  if(hide_overlaps!=0)printf("  identifying hidden faces -");
  for(i=0;i<nmeshes;i++){
    int j;
    mesh *meshi;

    meshi=meshinfo + i;

    for(j=0;j<6*meshi->nbptrs;j++){
      int k;
      facedata *facej;

      facej = meshi->faceinfo+j;
      facej->hidden=0;

    }
    if(hide_overlaps==0)continue;
    for(j=0;j<6*meshi->nbptrs;j++){
      int k;
      facedata *facej;

      facej = meshi->faceinfo+j;

      for(k=0;k<6*meshi->nbptrs;k++){
        facedata *facek;

        if(j==k)continue;
        facek = meshi->faceinfo+k;
        if(facek->hidden==1)continue;
        if(facej->xmin<facek->xmin||facej->xmax>facek->xmax)continue;
        if(facej->ymin<facek->ymin||facej->ymax>facek->ymax)continue;
        if(facej->zmin<facek->zmin||facej->zmax>facek->zmax)continue;
        facej->hidden=1;
        break;
      }
    }
  }
  if(hide_overlaps!=0)printf(" complete\n");
}

/* ------------------ allocate_faces ------------------------ */

void allocate_faces(){
  int i;
  mesh *meshi;
  int ntotal2=0;
  int abortflag=0;

  FREEMEMORY(face_transparent);
  for(i=0;i<nmeshes;i++){
    int ntotal;

    meshi = meshinfo + i;
    ntotal = 6*meshi->nbptrs + meshi->nvents+12;
    ntotal2 += ntotal;

    FREEMEMORY(meshi->faceinfo);
    FREEMEMORY(meshi->face_normals_single);
    FREEMEMORY(meshi->face_normals_double);
    FREEMEMORY(meshi->face_transparent_double);
    FREEMEMORY(meshi->face_textures);
    FREEMEMORY(meshi->face_outlines);
    if(ntotal>0){
      if(abortflag==0&&NewMemory((void **)&meshi->faceinfo,sizeof(facedata)*ntotal)==0){
        abortflag=1;
      }
      if(abortflag==0&&NewMemory((void **)&meshi->face_normals_single,sizeof(facedata *)*ntotal)==0){
        abortflag=1;
      }
      if(abortflag==0&&NewMemory((void **)&meshi->face_normals_double,sizeof(facedata *)*ntotal)==0){
        abortflag=1;
      }
      if(abortflag==0&&NewMemory((void **)&meshi->face_transparent_double,sizeof(facedata *)*ntotal)==0){
        abortflag=1;
      }
      if(abortflag==0&&NewMemory((void **)&meshi->face_textures,sizeof(facedata *)*ntotal)==0){
        abortflag=1;
      }
      if(abortflag==0&&NewMemory((void **)&meshi->face_outlines,sizeof(facedata *)*ntotal)==0){
        abortflag=1;
      }
    }
  }
  if(ntotal2>0){
    if(abortflag==0&&NewMemory((void **)&face_transparent,sizeof(facedata *)*ntotal2)==0){
      abortflag=1;
    }
  }
  if(abortflag==1){
    int mem_sum;
    float rmem;
    int nfaces_temp;

    mem_sum=0;
    nfaces_temp=0;
    ntotal2=0;
    for(i=0;i<nmeshes;i++){
      int ntotal;

      meshi = meshinfo + i;
 
      ntotal = 6*meshi->nbptrs + meshi->nvents+12;
      nfaces_temp+=(6*meshi->nbptrs);
      mem_sum+= ntotal*(sizeof(facedata)+5*sizeof(facedata *));
      ntotal2 += ntotal;
    }
    mem_sum+= ntotal2*sizeof(facedata *);
    printf("*** Fatal error.  Unable to allocate ");
    if(mem_sum>=1000000000){
      rmem=(float)mem_sum/1000000000.0;
      printf("%4.2f GB of memory\n",rmem);
      printf("                  for %i blockage faces.\n",nfaces_temp);
    }
    else if(mem_sum>1000000&&mem_sum<1000000000){
      rmem=(float)mem_sum/1000000.0;
      printf("%4.2f MB of memory\n",rmem);
      printf("                  for %i blockage faces.\n",nfaces_temp);
    }
    else{
      printf("%i bytes of memory\n",mem_sum);
      printf("                  for %i blockage faces.\n",nfaces_temp);
    }
    abortSV(_("memory allocation error."));
  }
  printf("\n");
}

/* ------------------ blockcompare ------------------------ */

int blockcompare( const void *arg1, const void *arg2 ){
  blockagedata *bc1,*bc2;
  int i1, i2;

  i1 = *(int *)arg1;
  i2 = *(int *)arg2;
  bc1 = selectblockinfo[i1];
  bc2 = selectblockinfo[i2];
  if(bc1->id!=-1&&bc2->id!=-1){
    if(bc1->id<bc2->id)return -1;
    if(bc1->id>bc2->id)return 1;
    return 0;
  }
  if(bc1->id==-1&&bc2->id==-1)return 0;
  if(bc1->id==-1)return 1;
  if(bc2->id==-1)return -1;
  return 0;
}

/* ------------------ update_selectblocks ------------------------ */

void update_selectblocks(void){
  mesh *meshi;
  blockagedata *bc;
  int i,j;
  int ntotal=0;
  int local_count=0;

  for(i=0;i<nmeshes;i++){
    meshi = meshinfo + i;
    ntotal += meshi->nbptrs;
  }
  if(ntotal==0)return;
  FREEMEMORY(selectblockinfo);
  FREEMEMORY(sortedblocklist);
//  FREEMEMORY(changed_idlist);

  NewMemory((void **)&selectblockinfo,sizeof(blockagedata *)*ntotal);
  NewMemory((void **)&sortedblocklist,sizeof(int)*ntotal);
 // NewMemory((void **)&changed_idlist,sizeof(int)*(ntotal+1));

  for(i=0;i<ntotal;i++){
    sortedblocklist[i]=i;
 //   changed_idlist[i]=0;
  }
  for(i=0;i<nmeshes;i++){
    meshi = meshinfo + i;
    for(j=0;j<meshi->nbptrs;j++){
      bc = meshi->blockageinfoptrs[j];
      selectblockinfo[local_count++]=bc;
    }
  }
  qsort(sortedblocklist,(size_t)ntotal,sizeof(int),blockcompare);
  nselectblocks=ntotal;
 // nchanged_idlist=ntotal;

}


/* ------------------ update_selectfaces ------------------------ */

void update_selectfaces(void){

  /* store info about faces that could be selected */

  int i,j,k;
  mesh *meshi;
  facedata *facek;
  selectdata *sd;
  int sides[]={DOWN_Y,UP_X,UP_Y,DOWN_X,DOWN_Z,UP_Z};
  
  FREEMEMORY(selectfaceinfo);

  ntotalfaces=0;
  for(i=0;i<nmeshes;i++){
    meshi=meshinfo + i;
    ntotalfaces += 6*meshi->nbptrs;
  }
  if(ntotalfaces==0)return;

  NewMemory((void **)&selectfaceinfo,ntotalfaces*sizeof(selectdata));

/* down y 
     up x 
     up y 
   down x 
   down z 
     up z */
  ntotalfaces=0;
  sd = selectfaceinfo;
  for(i=0;i<nmeshes;i++){
    meshi=meshinfo + i;
    for(j=0;j<meshi->nbptrs;j++){
      for(k=0;k<6;k++){
        facek = meshi->faceinfo + 6*j + sides[k];
        facek->blockageindex=j;
        facek->meshindex=i;
        ntotalfaces++;
        sd->mesh=i;
        sd->dir=facek->dir;
        sd->side=sides[k];
        sd->blockage=j;
        sd->facei=facek;
        sd->type=BLOCK_face;
        sd++;
      }
    }
  }
}

/* ------------------ update_smooth_blockages ------------------------ */

void update_smooth_blockages(void){
  int i, blocktotal;
  mesh *meshi;
  smoothblockage *sb;
  int j;

  blocktotal=0;

  if(menusmooth==0){
    STREAM_SB=fopen(filename_sb,"rb");
  }
  if(STREAM_SB!=NULL){
    time_t sb_modtime;

    sb_modtime=file_modtime(filename_sb);
    if(sb_modtime!=0&&smv_modtime!=0&&smv_modtime>sb_modtime){
      fclose(STREAM_SB);
      STREAM_SB=NULL;
    }
  }
  if(STREAM_SB!=NULL){
    int version;
    read_smoothobst=1;
    if(fread(&version,4,1,STREAM_SB)==1){
      rewind(STREAM_SB);
    }
    else{
      read_smoothobst=0;
    }
  }
  if(STREAM_SB==NULL){
    read_smoothobst=0;
    STREAM_SB=fopen(filename_sb,"wb");
  }

  for(i=0;i<nmeshes;i++){
    meshi = meshinfo + i;
    blocktotal += meshi->nbptrs;
  }
  if(blocktotal>0){
    printf("Initializing smooth blockage data - ");
    for(i=0;i<nmeshes;i++){
      meshi=meshinfo+i;

      for(j=0;j<meshi->nsmoothblockages_list;j++){
        if(read_smoothobst==1){
          printf("Reading smooth blockages %i of %i in mesh %i\n",j+1,meshi->nsmoothblockages_list,i+1);
        }
        else{
          printf("Smoothing blockages %i of %i in mesh %i\n",j+1,meshi->nsmoothblockages_list,i+1);
        }
        sb=meshi->smoothblockages_list+j;

        getsmoothblockparms(meshi,sb);
        MakeIsoBlockages(meshi,sb);
      }
    }
    printf(" - completed \n");
  }
  fclose(STREAM_SB);
  STREAM_SB=NULL;
  blocksneedsmoothing=0;
  updatesmoothblocks=0;
  smoothing_blocks=0;
}

/* ------------------ getsmoothblockage ------------------------ */

smoothblockage *getsmoothblockage(mesh *meshi,float tt){
  int j;
  smoothblockage *sb,*sb2;


  sb=meshi->smoothblockages_list;
  if(sb==NULL)return NULL;
  if(tt<0.0)return sb;

  for(j=1;j<meshi->nsmoothblockages_list-1;j++){
    sb=meshi->smoothblockages_list+j;
    sb2=sb+1;
    if(sb->time<=tt&&tt<sb2->time)return sb;
  }
  sb=meshi->smoothblockages_list+meshi->nsmoothblockages_list-1;
  return sb;
}

/* ------------------ isblockagevisible ------------------------ */

int isblockagevisible(blockagedata *bc, float local_time){
  int listindex,val;

  if(bc->showhide==NULL||local_time<0.0)return 1;
  listindex=getindex(local_time,bc->showtime,bc->nshowtime);
  val = bc->showhide[listindex];
  return val;
}


/* ------------------ getsmoothblockparams ------------------------ */

void getsmoothblockparms(mesh *meshi, smoothblockage *sb){
  int i,j;
  blockagedata *bc,*bc2;
  int nsmoothcolors=0;
  int fail;
  
  /* number of unique smooth block colors */

  for(i=0;i<meshi->nbptrs;i++){
    bc = meshi->blockageinfoptrs[i];
    if(bc->type!=BLOCK_smooth||bc->del==1)continue;
    if(isblockagevisible(bc,sb->time)!=1){
      continue;
    }
    fail=0;
    for(j=0;j<i;j++){
      bc2=meshi->blockageinfoptrs[j];
      if(bc2->type!=BLOCK_smooth)continue;
      if(bc2->del==1)continue;
      if(isblockagevisible(bc2,sb->time)!=1)continue;
      if(ABS(bc->color[0]-bc2->color[0])>0.0001)continue;
      if(ABS(bc->color[1]-bc2->color[1])>0.0001)continue;
      if(ABS(bc->color[2]-bc2->color[2])>0.0001)continue;
      if(ABS(bc->color[3]-bc2->color[3])>0.0001)continue;
      fail=1;
      break;
    }
    if(fail==0){
      nsmoothcolors++;
    }
  }

  meshi->nsmoothblockagecolors=nsmoothcolors;

  /* free and allocate memory */

  FREEMEMORY(sb->smoothblockagecolors);
  FREEMEMORY(sb->smoothblockagesurfaces);
//  FREEMEMORY(meshi->smoothblockagecolors);
//  FREEMEMORY(meshi->blockagesurfaces);
  if(nsmoothcolors>0){
    NewMemory((void **)&meshi->smoothblockagecolors,4*nsmoothcolors*sizeof(float));
    NewMemory((void **)&meshi->blockagesurfaces,nsmoothcolors*sizeof(isosurface *));
  }
  
  sb->nsmoothblockagecolors=meshi->nsmoothblockagecolors;
  sb->smoothblockagecolors=meshi->smoothblockagecolors;
  sb->smoothblockagesurfaces=meshi->blockagesurfaces;
  
  for(i=0;i<nsmoothcolors;i++){
    meshi->blockagesurfaces[i]=NULL;
  }
  nsmoothcolors=0;
  
  /* smooth block colors */

  for(i=0;i<meshi->nbptrs;i++){
    bc = meshi->blockageinfoptrs[i];
    if(bc->type!=BLOCK_smooth||bc->del==1)continue;
    if(isblockagevisible(bc,sb->time)!=1)continue;
    fail=0;
    for(j=0;j<i;j++){
      bc2=meshi->blockageinfoptrs[j];
      if(bc2->type!=BLOCK_smooth)continue;
      if(bc2->del==1)continue;
      if(isblockagevisible(bc2,sb->time)!=1)continue;
      if(ABS(bc->color[0]-bc2->color[0])>0.0001)continue;
      if(ABS(bc->color[1]-bc2->color[1])>0.0001)continue;
      if(ABS(bc->color[2]-bc2->color[2])>0.0001)continue;
      if(ABS(bc->color[3]-bc2->color[3])>0.0001)continue;
      fail=1;
      break;
    }
    if(fail==0){
      meshi->smoothblockagecolors[4*nsmoothcolors]=bc->color[0];
      meshi->smoothblockagecolors[4*nsmoothcolors+1]=bc->color[1];
      meshi->smoothblockagecolors[4*nsmoothcolors+2]=bc->color[2];
      meshi->smoothblockagecolors[4*nsmoothcolors+3]=bc->color[3];
      nsmoothcolors++;
    }
  }
}
/* ------------------ ReadSmoothIsoSurface ------------------------ */

int ReadSmoothIsoSurface(isosurface *asurface){
  // use STREAM_SB
  int one;
  float color[4];

  if(STREAM_SB==NULL)return 1;

  fread(&one,4,1,STREAM_SB);
  if(feof(STREAM_SB)!=0)return 1;

  fread(&asurface->nvertices,4,1,STREAM_SB);
  if(feof(STREAM_SB)!=0)return 1;
  fread(&asurface->ntriangles,4,1,STREAM_SB);
  if(feof(STREAM_SB)!=0)return 1;

  fread(&asurface->xmin,4,1,STREAM_SB);
  if(feof(STREAM_SB)!=0)return 1;
  fread(&asurface->ymin,4,1,STREAM_SB);
  if(feof(STREAM_SB)!=0)return 1;
  fread(&asurface->zmin,4,1,STREAM_SB);
  if(feof(STREAM_SB)!=0)return 1;
  fread(&asurface->xyzmaxdiff,4,1,STREAM_SB);
  if(feof(STREAM_SB)!=0)return 1;

  fread(color,4,4,STREAM_SB);
  if(feof(STREAM_SB)!=0)return 1;
  asurface->color=getcolorptr(color);

  asurface->vertices=NULL;
  asurface->vertexnorm=NULL;
  if(asurface->nvertices>0){
    NewMemory((void **)&asurface->vertices,3*asurface->nvertices*sizeof(short));
    NewMemory((void **)&asurface->vertexnorm,3*asurface->nvertices*sizeof(short));
    fread(asurface->vertices,2,3*asurface->nvertices,STREAM_SB); // vertices scaled between 0 and 2**16-1
    if(feof(STREAM_SB)!=0){
      FREEMEMORY(asurface->vertices);
      FREEMEMORY(asurface->vertexnorm);
      return 1;
    }
    fread(asurface->vertexnorm,2,3*asurface->nvertices,STREAM_SB); // norms scaled between 0 and 2**16-1
    if(feof(STREAM_SB)!=0){
      FREEMEMORY(asurface->vertices);
      FREEMEMORY(asurface->vertexnorm);
      return 1;
    }
  }
  asurface->triangles=NULL;
  if(asurface->ntriangles>0){
    NewMemory((void **)&asurface->triangles,asurface->ntriangles*sizeof(int));
    if(fread(asurface->triangles,4,asurface->ntriangles,STREAM_SB)<asurface->ntriangles){
      FREEMEMORY(asurface->triangles);
      return 1;
    }
  }
  return 0;
}

/* ------------------ WriteSmoothIsoSurface ------------------------ */

void WriteSmoothIsoSurface(isosurface *asurface){
  // use STREAM_SB
  int one=1;

  if(STREAM_SB==NULL)return;

  fwrite(&one,4,1,STREAM_SB);

  fwrite(&asurface->nvertices,4,1,STREAM_SB);
  fwrite(&asurface->ntriangles,4,1,STREAM_SB);

  fwrite(&asurface->xmin,4,1,STREAM_SB);
  fwrite(&asurface->ymin,4,1,STREAM_SB);
  fwrite(&asurface->zmin,4,1,STREAM_SB);
  fwrite(&asurface->xyzmaxdiff,4,1,STREAM_SB);

  fwrite(asurface->color,4,4,STREAM_SB);

  if(asurface->nvertices>0){
    fwrite(asurface->vertices,2,3*asurface->nvertices,STREAM_SB); // vertices scaled between 0 and 2**16-1
    fwrite(asurface->vertexnorm,2,3*asurface->nvertices,STREAM_SB); // norms scaled between 0 and 2**16-1
  }
  if(asurface->ntriangles>0){
    fwrite(asurface->triangles,4,asurface->ntriangles,STREAM_SB); // triangle indices
  }

}

/* ------------------ MakeIsoBlockages ------------------------ */

void MakeIsoBlockages(mesh *meshi, smoothblockage *sb){
  blockagedata *bc;
  float *cellcopy,*cell=NULL,*node=NULL,*nodecopy;
  int ib,i,j,k,iblockcolor;
  int imin, imax, jmin, jmax, kmin, kmax;
  float val;
  isosurface *asurface;
  float level;
  float vals[8];
  float *xplt2, *yplt2, *zplt2;
  float *xplt,*yplt,*zplt;
  float *rgbtemp,*rgbtemp2;
  int ibar,jbar,kbar;
  
  int ii, jj, kk;
  int im1, jm1, km1;
  int read_error=0;

#undef cellindex
#define cellindex(i,j,k) ((i)+(j)*(ibar+2)+(k)*(ibar+2)*(jbar+2))
#define nodeindex(i,j,k) ((i)+(j)*(ibar+3)+(k)*(ibar+3)*(jbar+3))

  xplt=meshi->xplt;
  yplt=meshi->yplt;
  zplt=meshi->zplt;
  ibar=meshi->ibar;
  jbar=meshi->jbar;
  kbar=meshi->kbar;

  NewMemory((void **)&cell,(ibar+2)*(jbar+2)*(kbar+2)*sizeof(float));
  NewMemory((void **)&node,(ibar+3)*(jbar+3)*(kbar+3)*sizeof(float));

  NewMemory((void **)&xplt2,(ibar+3)*sizeof(float));
  NewMemory((void **)&yplt2,(jbar+3)*sizeof(float));
  NewMemory((void **)&zplt2,(kbar+3)*sizeof(float));
  for(i=0;i<ibar+1;i++){
    xplt2[i+1]=xplt[i];
  }
  xplt2[0]=xplt[0]-(xplt[1]-xplt[0])/10.0;
  xplt2[ibar+2]=xplt[ibar]+(xplt[ibar]-xplt[ibar-1])/10.0;
  for(j=0;j<jbar+1;j++){
    yplt2[j+1]=yplt[j];
  }
  yplt2[0]=yplt[0]-(yplt[1]-yplt[0])/10.0;
  yplt2[jbar+2]=yplt[jbar]+(yplt[jbar]-yplt[jbar-1])/10.0;
  for(k=0;k<kbar+1;k++){
    zplt2[k+1]=zplt[k];
  }
  zplt2[0]=zplt[0]-(zplt[1]-zplt[0])/10.0;
  zplt2[kbar+2]=zplt[kbar]+(zplt[kbar]-zplt[kbar-1])/10.0;

  meshi->nsmoothblockagecolors=sb->nsmoothblockagecolors;
  meshi->smoothblockagecolors=sb->smoothblockagecolors;
  meshi->blockagesurfaces=sb->smoothblockagesurfaces;

  for(iblockcolor=0;iblockcolor<meshi->nsmoothblockagecolors;iblockcolor++){

    rgbtemp=meshi->smoothblockagecolors + 4*iblockcolor;
    if(read_smoothobst==0){
      cellcopy=cell;
      nodecopy=node;
      for(i=0;i<(ibar+2)*(jbar+2)*(kbar+2);i++){
        *cellcopy++=0.0;
      }
      for(ib=0;ib<meshi->nbptrs;ib++){
        bc=meshi->blockageinfoptrs[ib];
        if(bc->type!=BLOCK_smooth||bc->del==1)continue;
        if(isblockagevisible(bc,sb->time)!=1)continue;
        rgbtemp2=bc->color;
        if(ABS(rgbtemp[0]-rgbtemp2[0])<0.0001&&
           ABS(rgbtemp[1]-rgbtemp2[1])<0.0001&&
           ABS(rgbtemp[2]-rgbtemp2[2])<0.0001&&
           ABS(rgbtemp[3]-rgbtemp2[3])<0.0001
           ){
          imin = bc->ijk[IMIN];
          imax = bc->ijk[IMAX];
          jmin = bc->ijk[JMIN];
          jmax = bc->ijk[JMAX];
          kmin = bc->ijk[KMIN];
          kmax = bc->ijk[KMAX];
          for(k=kmin;k<kmax;k++){
            for(j=jmin;j<jmax;j++){
              for(i=imin;i<imax;i++){
                cell[cellindex(i+1,j+1,k+1)]=1.0;
              }
            }
          }
        }
      }
      for(i=0;i<(ibar+3)*(jbar+3)*(kbar+3);i++){
        *nodecopy++=0.0;
      }
      for(kk=1;kk<kbar+2;kk++){
        if(kk==1){
          km1=kk;
          k=kk;
        }
        else if(kk==kbar+1){
          km1=kk-1;
          k=kk-1;
        }
        else{
          k=kk; 
          km1=kk-1;
        }
        for(jj=1;jj<jbar+2;jj++){
          if(jj==1){
            jm1=jj;
            j=jj;
          }
          else if(jj==jbar+1){
            jm1=jj-1;
            j=jj-1;
          }
          else{
            j=jj; 
            jm1=jj-1;
          }
          for(ii=1;ii<ibar+2;ii++){
            if(ii==1){
              im1=ii;
              i=ii;
            }
            else if(ii==ibar+1){
              im1=ii-1;
              i=ii-1;
            }
            else{
              i=ii; 
              im1=ii-1;
            }
            vals[0]=cell[cellindex(im1,jm1,km1)];
            vals[1]=cell[cellindex(im1,jm1,k)];
            vals[2]=cell[cellindex(im1,j  ,km1)];
            vals[3]=cell[cellindex(im1,j  ,k)];
            vals[4]=cell[cellindex(i  ,jm1,km1)];
            vals[5]=cell[cellindex(i  ,jm1,k)];
            vals[6]=cell[cellindex(i  ,j  ,km1)];
            vals[7]=cell[cellindex(i  ,j  ,k)];
  
            val = (vals[0]+vals[1]+vals[2]+vals[3]+vals[4]+vals[5]+vals[6]+vals[7]+0.01)/8.0;
            node[nodeindex(ii,jj,kk)]=val;
          }
        }
      }
    }
    asurface=NULL;
    NewMemory((void **)&asurface,sizeof(isosurface));
    level=0.250;
    InitIsosurface(asurface, level, rgbtemp,0);
    if(read_smoothobst==1){
      // read in smoothed iso info here
      read_error=ReadSmoothIsoSurface(asurface);
      if(read_error!=0){
        read_smoothobst=0;
        printf(" *** warning: unexpected end of file encountered while\n");
        printf("              reading the smooth blockage file.\n");
      }
    }
    if(read_smoothobst==0){
      float dlevel=-1.0;
      
      GetIsosurface(asurface, node, NULL, NULL, level, dlevel,
                     xplt2, ibar+3, yplt2, jbar+3, zplt2, kbar+3);
      GetNormalSurface(asurface);
      CompressIsosurface(asurface,1,
          xplt2[0],xplt2[ibar+2],
          yplt2[0],yplt2[jbar+2],
          zplt2[0],zplt2[kbar+2]);
      SmoothIsoSurface(asurface);
      // write out smoothed iso info here
      if(read_error==0)WriteSmoothIsoSurface(asurface);
    }

    if(meshi->blockagesurfaces!=NULL)meshi->blockagesurfaces[iblockcolor]=asurface;
    meshi->blockagesurface=asurface;
    if(sb->smoothblockagesurfaces!=NULL)sb->smoothblockagesurfaces[iblockcolor]=asurface;
  }
  FREEMEMORY(node); FREEMEMORY(cell);
  FREEMEMORY(xplt2);FREEMEMORY(yplt2);FREEMEMORY(zplt2);
  return;
}
/* ------------------ MakeIsoBlockages2 ------------------------ */

void MakeIsoBlockages2(mesh *meshi, smoothblockage *sb){
//xxx experimental smooth blockage generation routine
  blockagedata *bc;
  float *cell=NULL;
  int ib,i,j,k,iblockcolor;
  int imin, imax, jmin, jmax, kmin, kmax;
  isosurface *asurface;
  float level;
  float *xplt2, *yplt2, *zplt2;
  float *XPLT2, *YPLT2, *ZPLT2;
  float *xplt,*yplt,*zplt;
  float *rgbtemp,*rgbtemp2;
  int ibar,jbar,kbar;
  int nx2, ny2, nz2;
  
  int read_error=0;

#define cellindex2(i,j,k) ((i+1)+(j+1)*nx2+(k+1)*nx2*ny2)

  xplt=meshi->xplt;
  yplt=meshi->yplt;
  zplt=meshi->zplt;
  ibar=meshi->ibar;
  jbar=meshi->jbar;
  kbar=meshi->kbar;
  nx2=ibar+2;
  ny2=jbar+2;
  nz2=kbar+2;

  NewMemory((void **)&cell,nx2*ny2*nz2*sizeof(float));

  NewMemory((void **)&XPLT2,nx2*sizeof(float));
  NewMemory((void **)&YPLT2,ny2*sizeof(float));
  NewMemory((void **)&ZPLT2,nz2*sizeof(float));

  xplt2 = XPLT2+1;
  yplt2 = YPLT2+1;
  zplt2 = ZPLT2+1;

  for(i=0;i<ibar;i++){
    xplt2[i]=(xplt[i]+xplt[i+1])/2.0;
  }
  xplt2[-1]  =xplt2[0]-(xplt[1]-xplt[0]);
  xplt2[ibar]=xplt2[ibar-1]+(xplt[ibar]-xplt[ibar-1]);

  for(i=0;i<jbar;i++){
    yplt2[i]=(yplt[i]+yplt[i+1])/2.0;
  }
  yplt2[-1]=yplt2[0]-(yplt[1]-yplt[0]);
  yplt2[jbar]=yplt2[jbar-1]+(yplt[jbar]-yplt[jbar-1]);

  for(i=0;i<kbar;i++){
    zplt2[i]=(zplt[i]+zplt[i+1])/2.0;
  }
  zplt2[-1]=zplt2[0]-(zplt[1]-zplt[0]);
  zplt2[kbar]=zplt2[kbar-1]+(zplt[kbar]-zplt[kbar-1]);

  meshi->nsmoothblockagecolors=sb->nsmoothblockagecolors;
  meshi->smoothblockagecolors=sb->smoothblockagecolors;
  meshi->blockagesurfaces=sb->smoothblockagesurfaces;

  for(iblockcolor=0;iblockcolor<meshi->nsmoothblockagecolors;iblockcolor++){

    rgbtemp=meshi->smoothblockagecolors + 4*iblockcolor;
    if(read_smoothobst==0){
      for(i=0;i<nx2*ny2*nz2;i++){
        cell[i]=0.0;
      }
      for(ib=0;ib<meshi->nbptrs;ib++){
        bc=meshi->blockageinfoptrs[ib];
        if(bc->type!=BLOCK_smooth||bc->del==1)continue;
        if(isblockagevisible(bc,sb->time)!=1)continue;
        rgbtemp2=bc->color;
        if(ABS(rgbtemp[0]-rgbtemp2[0])<0.0001&&
           ABS(rgbtemp[1]-rgbtemp2[1])<0.0001&&
           ABS(rgbtemp[2]-rgbtemp2[2])<0.0001&&
           ABS(rgbtemp[3]-rgbtemp2[3])<0.0001
           ){
          imin = bc->ijk[IMIN];
          imax = bc->ijk[IMAX];
          jmin = bc->ijk[JMIN];
          jmax = bc->ijk[JMAX];
          kmin = bc->ijk[KMIN];
          kmax = bc->ijk[KMAX];
          for(k=kmin;k<kmax;k++){
            for(j=jmin;j<jmax;j++){
              for(i=imin;i<imax;i++){
                cell[cellindex2(i,j,k)]=1.0;
              }
            }
          }
        }
      }
    }
    asurface=NULL;
    NewMemory((void **)&asurface,sizeof(isosurface));
    level=0.50;
    InitIsosurface(asurface, level, rgbtemp,0);
    if(read_smoothobst==1){
      // read in smoothed iso info here
      read_error=ReadSmoothIsoSurface(asurface);
      if(read_error!=0){
        read_smoothobst=0;
        printf(" *** warning: unexpected end of file encountered while\n");
        printf("              reading the smooth blockage file.\n");
      }
    }
    else{
      float dlevel=-1.0;
      
      GetIsosurface(asurface, cell, NULL, NULL, level, dlevel,
                     xplt2-1, nx2, yplt2-1, ny2, zplt2-1, nz2);
      GetNormalSurface(asurface);
      CompressIsosurface(asurface,1,
          xplt2[-1],xplt2[ibar],
          yplt2[-1],yplt2[jbar],
          zplt2[-1],zplt2[kbar]);
      SmoothIsoSurface(asurface);
      if(read_error==0)WriteSmoothIsoSurface(asurface);
    }

    if(meshi->blockagesurfaces!=NULL)meshi->blockagesurfaces[iblockcolor]=asurface;
    meshi->blockagesurface=asurface;
    if(sb->smoothblockagesurfaces!=NULL)sb->smoothblockagesurfaces[iblockcolor]=asurface;
  }
  FREEMEMORY(cell); 
  FREEMEMORY(XPLT2);
  FREEMEMORY(YPLT2);
  FREEMEMORY(ZPLT2);
  return;
}

/* ------------------ update_demo ------------------------ */

void init_demo(float rad, int nlat, int nlong){
  int i,j;
  float phi, psi;
  extern float *sphere_xyz;
  extern int update_demo;
  float *s_xyz;

  if(nlat<=0||nlong<=0)return;
  update_demo=0;
  FREEMEMORY(sphere_xyz);
  NewMemory((void **)&sphere_xyz,3*nlat*(nlong+1)*sizeof(float));
  s_xyz=sphere_xyz;
  for(j=0;j<nlong+1;j++){
    phi=-PI + 2.0*PI*j/nlong;
    for(i=0;i<nlat;i++){
      psi = -PI/2.0 + i*PI/(nlat-1);
      *s_xyz++ = 0.2143 + rad*cos(psi)*cos(phi);
      *s_xyz++ = 0.2143 + rad*cos(psi)*sin(phi);
      *s_xyz++ = 0.5 + rad*sin(psi);
    }
  }
}

/* ------------------ calcNormal3 ------------------------ */

void calcNormal3(const float *v1, 
                 const float *v2, 
                 const float *v3, 
                 float *out){
  float u[3], v[3];
  int i;


  for(i=0;i<3;i++){
    u[i]=v2[i]-v1[i];
    v[i]=v3[i]-v1[i];
  }


  out[0] = u[1]*v[2] - u[2]*v[1];
  out[1] = u[2]*v[0] - u[0]*v[2];
  out[2] = u[0]*v[1] - u[1]*v[0];

  ReduceToUnit(out);

}

/* ------------------ calcNormal4 ------------------------ */

void calcNormal4(const float *v1, 
                 float *out){
  out[0]=v1[0]-0.2143;
  out[1]=v1[1]-0.2143;
  out[2]=v1[2]-0.5;
  


  ReduceToUnit(out);

}

/* ------------------ draw_demo ------------------------ */

void draw_demo2(int option){
      demo_mode++;
      glBegin(GL_QUADS);
      if(demo_mode%2==0){
        glColor3f(1.0,0.0,0.0);
      }
      else{
        glColor3f(0.0,0.0,1.0);
      }
      glVertex3f(0.0,0.3,0.0);
      glVertex3f(0.0,0.6,0.0);
      glVertex3f(1.0,0.6,0.0);
      glVertex3f(1.0,0.3,0.0);
      glEnd();
}
void draw_demo(int nlat, int nlong){
  int i, j;
  extern float *sphere_xyz;
  extern int update_demo,demo_mode;
  float red, green;
//  float blue;
  float *xyz;
  float *xyz00,*xyz01,*xyz10,*xyz11;
  float norm[3];
//  float norm1[3],norm2[3],norm3[3];
//  float denom;
  float specular[4]={0.8,0.8,0.8,1.0};
#define sphere_index(ilat,ilong) (3*((ilong)*nlat + (ilat)))

  if(nlat<=0||nlong<=0)return;
  if(update_demo==1)init_demo(0.4,nlat,nlong);
  switch (demo_mode){
    case 0:
      glPointSize(6.0);
      glColor3f(0.0,0.0,1.0);
      glBegin(GL_POINTS);
      for(j=0;j<nlong;j++){
        for(i=0;i<nlat;i++){
          xyz = sphere_xyz + sphere_index(i,j);
          glVertex3fv(xyz);
        }
      }
      glEnd();
      break;
    case 1:
      glLineWidth(2.0);
      glBegin(GL_LINES);
      glColor3f(0.0,0.0,1.0);
      for(j=0;j<nlong;j++){
        for(i=0;i<nlat-1;i++){
          xyz00 = sphere_xyz + sphere_index(i,j);
          xyz10 = sphere_xyz + sphere_index(i,j+1);
          xyz01 = sphere_xyz + sphere_index(i+1,j);
          xyz11 = sphere_xyz + sphere_index(i+1,j+1);
          glVertex3fv(xyz00);
          glVertex3fv(xyz01);
          glVertex3fv(xyz00);
          glVertex3fv(xyz10);
        }
      }
      glEnd();
      break;
    case 2:
      glBegin(GL_TRIANGLES);
      glColor3f(0.0,0.0,1.0);
      for(j=0;j<nlong;j++){
        for(i=0;i<nlat-1;i++){
          xyz00 = sphere_xyz + sphere_index(i,j);
          xyz10 = sphere_xyz + sphere_index(i,j+1);
          xyz01 = sphere_xyz + sphere_index(i+1,j);
          xyz11 = sphere_xyz + sphere_index(i+1,j+1);
          glVertex3fv(xyz00);
          glVertex3fv(xyz01);
          glVertex3fv(xyz11);
          glVertex3fv(xyz00);
          glVertex3fv(xyz11);
          glVertex3fv(xyz10);
        }
      }
      glEnd();
      break;
    case 3:
    case 4:
    case 5:
//#define COLOR(x) (1.0+((x)-0.2143)/0.3)/2.0 
#define COLOR(x) 0.0
      glEnable(GL_LIGHTING);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,specular);
      glEnable(GL_COLOR_MATERIAL);
      for(j=0;j<nlong;j++){
        for(i=0;i<nlat-1;i++){
          xyz00 = sphere_xyz + sphere_index(i,j);
          xyz10 = sphere_xyz + sphere_index(i,j+1);
          xyz01 = sphere_xyz + sphere_index(i+1,j);
          xyz11 = sphere_xyz + sphere_index(i+1,j+1);

          if(demo_mode==3)calcNormal3(xyz00,xyz11,xyz01,norm);
          glBegin(GL_TRIANGLES);
          if(demo_mode!=3)calcNormal4(xyz00,norm);
          glNormal3fv(norm);
          red = COLOR(xyz00[0]);
          green = COLOR(xyz00[1]);
          glColor3f(red,green,1.0);
          glVertex3fv(xyz00);

          red = COLOR(xyz11[0]);
          green = COLOR(xyz11[1]);
          if(demo_mode!=3)calcNormal4(xyz11,norm);
          glNormal3fv(norm);
          glColor3f(red,green,1.0);
          glVertex3fv(xyz11);

          red = COLOR(xyz01[0]);
          green = COLOR(xyz01[1]);
          glColor3f(red,green,1.0);
          if(demo_mode!=3)calcNormal4(xyz01,norm);
          glNormal3fv(norm);
          glVertex3fv(xyz01);

          glEnd();
          if(demo_mode==5){
            glLineWidth(2.0);
            glBegin(GL_LINES);
            glColor3f(0.0,0.0,0.0);
            glVertex3fv(xyz00);
            glVertex3f(xyz00[0]+norm[0]/10.0,xyz00[1]+norm[1]/10.0,xyz00[2]+norm[2]/10.0);
            glEnd();
          }

          if(demo_mode==3)calcNormal3(xyz00,xyz11,xyz01,norm);
          glBegin(GL_TRIANGLES);
          if(demo_mode==3)calcNormal3(xyz00,xyz11,xyz01,norm);
          red = COLOR(xyz00[0]);
          green = COLOR(xyz00[1]);
          glColor3f(red,green,1.0);
          if(demo_mode!=3)calcNormal4(xyz00,norm);
          glNormal3fv(norm);
          glVertex3fv(xyz00);

          red = COLOR(xyz10[0]);
          green = COLOR(xyz10[1]);
          glColor3f(red,green,1.0);
          if(demo_mode!=3)calcNormal4(xyz10,norm);
          glNormal3fv(norm);
          glVertex3fv(xyz10);

          red = COLOR(xyz11[0]);
          green = COLOR(xyz11[1]);
          glColor3f(red,green,1.0);
          if(demo_mode!=3)calcNormal4(xyz11,norm);
          glNormal3fv(norm);
          glVertex3fv(xyz11);

          glEnd();
        }
      }
      glDisable(GL_LIGHTING);
      glDisable(GL_COLOR_MATERIAL);
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
}

/* ------------------ innitticks ------------------------ */

void init_user_ticks(void){
  int i;

  user_tick_min[0]=1000000000.0;
  user_tick_min[1]=1000000000.0;
  user_tick_min[2]=1000000000.0;
  user_tick_max[0]=-1000000000.0;
  user_tick_max[1]=-1000000000.0;
  user_tick_max[2]=-1000000000.0;

  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi = meshinfo + i;
    if(meshi->boxmin[0]<user_tick_min[0])user_tick_min[0]=meshi->boxmin[0];
    if(meshi->boxmin[1]<user_tick_min[1])user_tick_min[1]=meshi->boxmin[1];
    if(meshi->boxmin[2]<user_tick_min[2])user_tick_min[2]=meshi->boxmin[2];

    if(meshi->boxmax[0]>user_tick_max[0])user_tick_max[0]=meshi->boxmax[0];
    if(meshi->boxmax[1]>user_tick_max[1])user_tick_max[1]=meshi->boxmax[1];
    if(meshi->boxmax[2]>user_tick_max[2])user_tick_max[2]=meshi->boxmax[2];
  }

  user_tick_origin[0]=user_tick_min[0];
  user_tick_origin[1]=user_tick_min[1];
  user_tick_origin[2]=user_tick_min[2];

  user_tick_step[0]=1.0;
  user_tick_step[1]=1.0;
  user_tick_step[2]=1.0;

  user_tick_sub=5;

  user_tick_length=0.1;
  user_tick_width=2.0;

}

/* ------------------ drawticks ------------------------ */

void draw_user_ticks(void){
  int i;
  float xyz[3],xyz2[3];
  float tick_origin[3], step[3];
  int show_tick_x, show_tick_y, show_tick_z;

#define MIN_DTICK 0.0
#define TEXT_FACTOR 1.5

  user_tick_option=get_tick_dir(modelview_scratch);

  if(auto_user_tick_placement==0){
    tick_origin[0]=user_tick_origin[0];
    tick_origin[1]=user_tick_origin[1];
    tick_origin[2]=user_tick_origin[2];
    step[0]=user_tick_step[0];
    step[1]=user_tick_step[1];
    step[2]=user_tick_step[2];
    show_tick_x = user_tick_show_x;
    show_tick_y = user_tick_show_y;
    show_tick_z = user_tick_show_z;
  }
  if(auto_user_tick_placement==1){
    step[0]=ABS(user_tick_step[0]);
    step[1]=ABS(user_tick_step[1]);
    step[2]=ABS(user_tick_step[2]);
    show_tick_x=0;
    show_tick_y=0;
    show_tick_z=0;
    switch (user_tick_option){
      case -1:
        tick_origin[0] = user_tick_origin[0];
        tick_origin[1] = user_tick_origin[1];
        tick_origin[2] = user_tick_origin[2];
        show_tick_x = 1;
        show_tick_z = 1;
        break;
      case 1:
        tick_origin[0] = user_tick_origin[0];
        tick_origin[1] = user_tick_max[1];
        step[1] = -step[1];
        tick_origin[2] = user_tick_origin[2];
        show_tick_x = 1;
        show_tick_z = 1;
        break;
      case -2:
        tick_origin[0] = user_tick_origin[0];
        tick_origin[1] = user_tick_origin[1];
        tick_origin[2] = user_tick_origin[2];
        show_tick_y = 1;
        show_tick_z = 1;
        break;
      case 2:
        tick_origin[0] = user_tick_max[0];
        step[0] = -step[0];
        tick_origin[1] = user_tick_origin[1];
        tick_origin[2] = user_tick_origin[2];
        show_tick_y = 1;
        show_tick_z = 1;
        break;
      case -3:
        tick_origin[0] = user_tick_origin[0];
        tick_origin[1] = user_tick_origin[1];
        tick_origin[2] = user_tick_origin[2];
        show_tick_x = 1;
        show_tick_y = 1;
        break;
      case 3:
        tick_origin[0] = user_tick_origin[0];
        tick_origin[1] = user_tick_origin[1];
        tick_origin[2] = user_tick_max[2];
        step[2] = -step[2];
        show_tick_x = 1;
        show_tick_y = 1;
        break;
      default:
        ASSERT(0);
        break;
    }
  }
  if(step[0]>MIN_DTICK){
    user_tick_nxyz[0]=ABS((user_tick_max[0]+1.0-tick_origin[0])/step[0]);
  }
  else if(step[0]<-MIN_DTICK){
    user_tick_nxyz[0]=ABS((tick_origin[0]+1.0-user_tick_min[0])/step[0]);
  }
  else{
    user_tick_nxyz[0]=0;
  }
  if(step[1]>MIN_DTICK){
    user_tick_nxyz[1]=ABS((user_tick_max[1]+1.0-tick_origin[1])/step[1]);
  }
  else if(step[1]<-MIN_DTICK){
    user_tick_nxyz[1]=ABS((tick_origin[1]+1.0-user_tick_min[1])/step[1]);
  }
  else{
    user_tick_nxyz[1]=0;
  }
  if(step[2]>MIN_DTICK){
    user_tick_nxyz[2]=ABS((user_tick_max[2]+1.0-tick_origin[2])/step[2]);
  }
  else if(step[2]<-MIN_DTICK){
    user_tick_nxyz[2]=ABS((tick_origin[2]+1.0-user_tick_min[2])/step[2]);
  }
  else{
    user_tick_nxyz[2]=0;
  }
  if(user_tick_option<0){
    user_tick_option=-user_tick_option;
  }
  glPushMatrix();
  glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);
  glTranslatef(-xbar0,-ybar0,-zbar0);
  glLineWidth(user_tick_width);

  //glPointSize(20.0);
  //glBegin(GL_POINTS);
  //glVertex3fv(tick_origin);
  //glEnd();
    
 //*** x axis tick/lables

 // major ticks
  if(show_tick_x==1){
    glBegin(GL_LINES);
    glColor3fv(foregroundcolor);
    for(i=0;i<user_tick_nxyz[0];i++){
      xyz[0]=tick_origin[0] + i*step[0];
      if(
        (step[0]>0.0&&xyz[0]>user_tick_max[0])||
        (step[0]<0.0&&xyz[0]<user_tick_min[0])
        )continue;
      xyz[1]=tick_origin[1];
      xyz[2]=tick_origin[2];
      if(user_tick_option==3){
        xyz2[0]=xyz[0];
        xyz2[1]=xyz[1]-user_tick_length*xyzmaxdiff;
        xyz2[2]=xyz[2];
      }
      else{
        xyz2[0]=xyz[0];
        xyz2[1]=xyz[1];
        xyz2[2]=xyz[2]-user_tick_length*xyzmaxdiff;
      }
      if(i==0){
        glVertex3fv(xyz);
        if(step[0]>0.0){
          glVertex3f(user_tick_max[0],xyz[1],xyz[2]);
        }
        else{
          glVertex3f(user_tick_min[0],xyz[1],xyz[2]);
        }
      }
      glVertex3fv(xyz);
      glVertex3fv(xyz2);
    }

// minor ticks

    if(user_tick_sub>1){
      for(i=1;i<user_tick_nxyz[0]*user_tick_sub;i++){
        if(i%user_tick_sub==0)continue;
        xyz[0]=tick_origin[0] + i*step[0]/(float)user_tick_sub;
        if(
          (step[0]>0.0&&xyz[0]>user_tick_max[0])||
          (step[0]<0.0&&xyz[0]<user_tick_min[0])
          )continue;
        xyz[1]=tick_origin[1];
        xyz[2]=tick_origin[2];
        if(user_tick_option==3){
          xyz2[0]=xyz[0];
          xyz2[1]=xyz[1]-user_tick_length*xyzmaxdiff/2.0;
          xyz2[2]=xyz[2];
        }
        else{
          xyz2[0]=xyz[0];
          xyz2[1]=xyz[1];
          xyz2[2]=xyz[2]-user_tick_length*xyzmaxdiff/2.0;
        }
        glVertex3fv(xyz);
        glVertex3fv(xyz2);
      }
    }
    glEnd();
    for(i=0;i<user_tick_nxyz[0];i++){
      char label[128];

      xyz[0]=tick_origin[0] + i*step[0];
      if((step[0]>0.0&&xyz[0]>user_tick_max[0])||(step[0]<0.0&&xyz[0]<user_tick_min[0]))continue;
      xyz[1]=tick_origin[1];
      xyz[2]=tick_origin[2];
      if(user_tick_option==3){
        xyz2[0]=xyz[0];
        xyz2[1]=xyz[1]-TEXT_FACTOR*user_tick_length*xyzmaxdiff;
        xyz2[2]=xyz[2];
      }
      else{
        xyz2[0]=xyz[0];
        xyz2[1]=xyz[1];
        xyz2[2]=xyz[2]-TEXT_FACTOR*user_tick_length*xyzmaxdiff;
      }
      sprintf(label,"%f",xyz[0]);
      trimzeros(label);
      output3Text(foregroundcolor,xyz2[0],xyz2[1],xyz2[2],label);
    }
  }

 //*** y axis tick/lables

 // major ticks

  if(show_tick_y){
    glLineWidth(user_tick_width);
    glBegin(GL_LINES);
    glColor3fv(foregroundcolor);
    for(i=0;i<user_tick_nxyz[1];i++){
      xyz[0]=tick_origin[0];
      xyz[1]=tick_origin[1] + i*step[1];
      if(
        (step[1]>0.0&&xyz[1]>user_tick_max[1])||
        (step[1]<0.0&&xyz[1]<user_tick_min[1])
        )continue;
      xyz[2]=tick_origin[2];
      if(user_tick_option==3){
        xyz2[0]=xyz[0]-user_tick_length*xyzmaxdiff;
        xyz2[1]=xyz[1];
        xyz2[2]=xyz[2];
      }
      else{
        xyz2[0]=xyz[0];
        xyz2[1]=xyz[1];
        xyz2[2]=xyz[2]-user_tick_length*xyzmaxdiff;
      }
      if(i==0){
        glVertex3fv(xyz);
        if(step[1]>0.0){
          glVertex3f(xyz[0],user_tick_max[1],xyz[2]);
        }
        else{
          glVertex3f(xyz[0],user_tick_min[1],xyz[2]);
        }
      }
      glVertex3fv(xyz);
      glVertex3fv(xyz2);
    }

// minor ticks

    if(user_tick_sub>1){
      for(i=1;i<user_tick_nxyz[1]*user_tick_sub;i++){
        if(i%user_tick_sub==0)continue;
        xyz[0]=tick_origin[0];
        xyz[1]=tick_origin[1] + i*step[1]/(float)user_tick_sub;
        if(
          (step[1]>0.0&&xyz[1]>user_tick_max[1])||
          (step[1]<0.0&&xyz[1]<user_tick_min[1])
          )continue;
        xyz[2]=tick_origin[2];
        if(user_tick_option==3){
          xyz2[0]=xyz[0]-user_tick_length*xyzmaxdiff/2.0;
          xyz2[1]=xyz[1];
          xyz2[2]=xyz[2];
        }
        else{
          xyz2[0]=xyz[0];
          xyz2[1]=xyz[1];
          xyz2[2]=xyz[2]-user_tick_length*xyzmaxdiff/2.0;
        }
        glVertex3fv(xyz);
        glVertex3fv(xyz2);
      }
    }
    glEnd();
    for(i=0;i<user_tick_nxyz[1];i++){
      char label[128];

      xyz[0]=tick_origin[0];
      xyz[1]=tick_origin[1] + i*step[1];
      if(
        (step[1]>0.0&&xyz[1]>user_tick_max[1])||
        (step[1]<0.0&&xyz[1]<user_tick_min[1])
        )continue;
      xyz[2]=tick_origin[2];
      xyz2[0]=xyz[0];
      if(user_tick_option==3){
        xyz2[0]=xyz[0]-TEXT_FACTOR*user_tick_length*xyzmaxdiff;
        xyz2[1]=xyz[1];
        xyz2[2]=xyz[2];
      }
      else{
        xyz2[0]=xyz[0];
        xyz2[1]=xyz[1];
        xyz2[2]=xyz[2]-TEXT_FACTOR*user_tick_length*xyzmaxdiff;
      }
      sprintf(label,"%f",xyz[1]);
      trimzeros(label);
      output3Text(foregroundcolor,xyz2[0],xyz2[1],xyz2[2],label);
    }
  }

 //*** z axis tick/lables

 // major ticks
  if(show_tick_z){
    glLineWidth(user_tick_width);
    glBegin(GL_LINES);
    glColor3fv(foregroundcolor);
    for(i=0;i<user_tick_nxyz[2];i++){
      xyz[0]=tick_origin[0];
      xyz[1]=tick_origin[1];
      xyz[2]=tick_origin[2] + i*step[2];
      if(
        (step[2]>0.0&&xyz[2]>user_tick_max[2])||
        (step[2]<0.0&&xyz[2]<user_tick_min[2])
        )continue;
      if(user_tick_option==2){
        xyz2[0]=xyz[0];
        xyz2[1]=xyz[1]-user_tick_length*xyzmaxdiff;
      }
      else{
        xyz2[0]=xyz[0]-user_tick_length*xyzmaxdiff;
        xyz2[1]=xyz[1];
      }
      xyz2[2]=xyz[2];
      if(i==0){
        glVertex3fv(xyz);
        if(step[2]>0.0){
          glVertex3f(xyz[0],xyz[1],user_tick_max[2]);
        }
        else{
          glVertex3f(xyz[0],xyz[1],user_tick_min[2]);
        }
      }
      glVertex3fv(xyz);
      glVertex3fv(xyz2);
    }

// minor ticks

    if(user_tick_sub>1){
      for(i=1;i<user_tick_nxyz[2]*user_tick_sub;i++){
        if(i%user_tick_sub==0)continue;
        xyz[0]=tick_origin[0];
        xyz[1]=tick_origin[1];
        xyz[2]=tick_origin[2] + i*step[2]/(float)user_tick_sub;
        if(
          (step[2]>0.0&&xyz[2]>user_tick_max[2])||
          (step[2]<0.0&&xyz[2]<user_tick_min[2])
          )continue;
        if(user_tick_option==2){
          xyz2[0]=xyz[0];
          xyz2[1]=xyz[1]-user_tick_length*xyzmaxdiff/2.0;
        }
        else{
          xyz2[0]=xyz[0]-user_tick_length*xyzmaxdiff/2.0;
          xyz2[1]=xyz[1];
        }
        xyz2[2]=xyz[2];
        glVertex3fv(xyz);
        glVertex3fv(xyz2);
      }
    }
    glEnd();
    for(i=0;i<user_tick_nxyz[2];i++){
      char label[128];

      xyz[0]=tick_origin[0];
      xyz[1]=tick_origin[1];
      xyz[2]=tick_origin[2] + i*step[2];
      if(
        (step[2]>0.0&&xyz[2]>user_tick_max[2])||
        (step[2]<0.0&&xyz[2]<user_tick_min[2])
        )continue;
      if(user_tick_option==2){
        xyz2[0]=xyz[0];
        xyz2[1]=xyz[1]-TEXT_FACTOR*user_tick_length*xyzmaxdiff;
      }
      else{
        xyz2[0]=xyz[0]-TEXT_FACTOR*user_tick_length*xyzmaxdiff;
        xyz2[1]=xyz[1];
      }
      xyz2[2]=xyz[2];
      sprintf(label,"%f",xyz[2]);
      trimzeros(label);
      output3Text(foregroundcolor,xyz2[0],xyz2[1],xyz2[2],label);
    }
  }

  glPopMatrix();
}


/* ------------------ getdir ------------------------ */

int get_tick_dir(float *mm){
    /*
      ( m0 m4 m8  m12 ) (x)    (0)
      ( m1 m5 m9  m13 ) (y)    (0)
      ( m2 m6 m10 m14 ) (z)  = (0)
      ( m3 m7 m11 m15 ) (1)    (1)

       ( m0 m4  m8 )      (m12)
   Q=  ( m1 m5  m9 )  u = (m13)
       ( m2 m6 m10 )      (m14)
      
      (Q   u) (x)     (0)      
      (v^T 1) (y)   = (1)
       
      m3=m7=m11=0, v^T=0, y=1   Qx+u=0 => x=-Q^Tu
    */
  int i,ii;
  float norm[3],scalednorm[3];
  float normdir[3];
  float absangle,cosangle,minangle;
  int iminangle;
  float pi;

  pi=4.0*atan(1.0);

  xyzeyeorig[0] = -(mm[0]*mm[12]+mm[1]*mm[13]+ mm[2]*mm[14])/mscale[0];
  xyzeyeorig[1] = -(mm[4]*mm[12]+mm[5]*mm[13]+ mm[6]*mm[14])/mscale[1];
  xyzeyeorig[2] = -(mm[8]*mm[12]+mm[9]*mm[13]+mm[10]*mm[14])/mscale[2];
  
  minangle=1000000.0;

  for(i=-3;i<=3;i++){
    if(i==0)continue;
    ii = i;
    if(i<0)ii=-i;
    norm[0]=0.0;
    norm[1]=0.0;
    norm[2]=0.0;
    switch (ii){
    case 1:
      if(i<0)norm[1]=-1.0;
      if(i>0)norm[1]=1.0;
      break;
    case 2:
      if(i<0)norm[0]=-1.0;
      if(i>0)norm[0]=1.0;
      break;
    case 3:
      if(i<0)norm[2]=-1.0;
      if(i>0)norm[2]=1.0;
      break;
    default:
      ASSERT(0);
      break;
    }
    scalednorm[0]=norm[0]*mscale[0];
    scalednorm[1]=norm[1]*mscale[1];
    scalednorm[2]=norm[2]*mscale[2];

    normdir[0] = mm[0]*scalednorm[0] + mm[4]*scalednorm[1] + mm[8]*scalednorm[2];
    normdir[1] = mm[1]*scalednorm[0] + mm[5]*scalednorm[1] + mm[9]*scalednorm[2];
    normdir[2] = mm[2]*scalednorm[0] + mm[6]*scalednorm[1] + mm[10]*scalednorm[2];

    cosangle = normdir[2]/sqrt(normdir[0]*normdir[0]+normdir[1]*normdir[1]+normdir[2]*normdir[2]);
    if(cosangle>1.0)cosangle=1.0;
    if(cosangle<-1.0)cosangle=-1.0;
    absangle=acos(cosangle)*180.0/pi;
    if(absangle<0.0)absangle=-absangle;
    if(absangle<minangle){
      iminangle=i;
      minangle=absangle;
    }
  }
  return iminangle;
}

/* ------------------ drawticks ------------------------ */

void drawticks(void){
  int i,j;
  tickdata *ticki;
  float *dxyz,xyz[3],xyz2[3],*begt,*endt,dbar[3];

  for(i=0;i<nticks;i++){
    ticki = tickinfo + i;
    begt = ticki->begin;
    endt = ticki->end;

    glLineWidth(ticki->width);
    glBegin(GL_LINES);
    if(ticki->useforegroundcolor==1){
      glColor3fv(foregroundcolor);
    }
    else{
      glColor3fv(ticki->rgb);
    }

    dxyz=ticki->dxyz;
    if(ticki->nbars>1){
      dbar[0]=(endt[0]-begt[0])/(float)(ticki->nbars-1);
      dbar[1]=(endt[1]-begt[1])/(float)(ticki->nbars-1);
      dbar[2]=(endt[2]-begt[2])/(float)(ticki->nbars-1);
    }
    else{
      dbar[0] = 0.0;
      dbar[1] = 0.0;
      dbar[2] = 0.0;
    }

    for(j=0;j<ticki->nbars;j++){
      xyz[0]=begt[0] + j*dbar[0];
      xyz[1]=begt[1] + j*dbar[1];
      xyz[2]=begt[2] + j*dbar[2];
      xyz2[0]=xyz[0]+dxyz[0];
      xyz2[1]=xyz[1]+dxyz[1];
      xyz2[2]=xyz[2]+dxyz[2];
      xyz[0] = (xyz[0] - xbar0)/xyzmaxdiff;
      xyz[1] = (xyz[1] - ybar0)/xyzmaxdiff;
      xyz[2] = (xyz[2] - zbar0)/xyzmaxdiff;
      xyz2[0] = (xyz2[0] - xbar0)/xyzmaxdiff;
      xyz2[1] = (xyz2[1] - ybar0)/xyzmaxdiff;
      xyz2[2] = (xyz2[2] - zbar0)/xyzmaxdiff;
      glVertex3fv(xyz);
      glVertex3fv(xyz2);
    }
    glEnd();

  }
}

/* ------------------ drawBlockages ------------------------ */

void drawBlockages(int mode, int trans_flag){

  mesh *meshi;
  int smoothnorms;
  int i,j;
  cadgeom *cd;
  int drawing_smooth, drawing_transparent, drawing_blockage_transparent, drawing_vent_transparent;

  get_drawing_parms(&drawing_smooth, &drawing_transparent, &drawing_blockage_transparent, &drawing_vent_transparent);

  if(drawing_smooth==1&&showedit_dialog==0){
    if(xyz_clipplane!=0)glDisable(GL_CULL_FACE);
    for(i=0;i<nmeshes;i++){
      meshi = meshinfo + i;

      for(j=0;j<meshi->nsmoothblockagecolors;j++){
        isosurface *bsurface;
        smoothnorms=1;
        if(meshi->blockagesurface!=NULL){
          bsurface=meshi->blockagesurfaces[j];
          drawstaticiso(bsurface,1,smoothnorms,trans_flag,1,plot3dlinewidth);
        }
      }
    }
    SNIFF_ERRORS("after drawblocks");
    if(xyz_clipplane!=0)glEnable(GL_CULL_FACE);
  }
  if(trans_flag!=DRAW_TRANSPARENT&&blocklocation!=BLOCKlocation_cad){
    if(mode==SELECT){
      if(blockageSelect==1){
        drawselect_faces();
        return;
      }
    }
    else{
      if(use_new_drawface==1){
        draw_faces();
      }
      else{
        draw_facesOLD();
      }
    }
  }

  if(blocklocation==BLOCKlocation_cad||(ncadgeom!=0&&show_cad_and_grid==1)){
    for(i=0;i<ncadgeom;i++){
      cd=cadgeominfo+i;
      if(cd->version==1){
        if(trans_flag==DRAW_TRANSPARENT)continue;
        if(xyz_clipplane==2){
          setClipPlanes(1);
        }
        drawcadgeom(cd);
        if(xyz_clipplane==2){
          unsetClipPlanes();
        }
      }
      else if(cd->version==2){
        if(xyz_clipplane==2){
          setClipPlanes(1);
        }
        drawcad2geom(cd,trans_flag);
        if(xyz_clipplane==2){
          unsetClipPlanes();
        }
      }
    }
  }
}

/* ------------------ snap_view_angles ------------------------ */

void snap_view_angles(void){
  float *az, *elev;
  int iaz, ielev;

#define DELTA 45.0

  az = camera_current->angle_zx;
  elev = camera_current->angle_zx+1;

  if(*az>0.0){
    iaz = (*az+DELTA/2.0)/DELTA;
  }
  else{
    iaz = (*az-DELTA/2.0)/DELTA;
  }
  *az = (int)(DELTA*iaz);

  if(*elev>0.0){
    ielev = (*elev+DELTA/2.0)/DELTA;
  }
  else{
    ielev = (*elev-DELTA/2.0)/DELTA;
  }
  *elev = (int)(DELTA*ielev);
  update_trainer_moves();
  camera_current->dirty=1;

}

/* ------------------ get_drawing_parms ------------------------ */

void get_drawing_parms(int *drawing_smooth, int *drawing_transparent, int *drawing_blockage_transparent, int *drawing_vent_transparent){
  *drawing_smooth=0;
  *drawing_transparent=0;
  *drawing_blockage_transparent=0;
  *drawing_vent_transparent=0;
  if(ntotal_smooth_blockages>0&&updatesmoothblocks==0){
    if(visSmoothAsNormal==0||visBlocks==visBLOCKAsInput){
      if(visBlocks!=visBLOCKOutline&&visBlocks!=visBLOCKSolidOutline&&visBlocks!=visBLOCKHide)*drawing_smooth=1;
    }
  }
  if(ntransparentblocks>0){
    if(visTransparentBlockage==1||visBlocks==visBLOCKAsInput){
      if(visBlocks!=visBLOCKOutline&&visBlocks!=visBLOCKHide){
        *drawing_transparent=1;
        *drawing_blockage_transparent=1;
      }
    }
  }
  if(ntransparentvents>0&&show_transparent_vents==1){
    *drawing_transparent=1;
    *drawing_vent_transparent=1;
  }
}

/* ------------------ drawfacesBAK ------------------------ */

void draw_facesOLD(){
  int i,j;
  mesh *meshi;
  facedata *facei;
  float *vertices,*tvertices;
  texture *texti;
  float *new_color,*old_color=NULL;
  int **showtimelist_handle, *showtimelist;
  float up_color[4]={0.9,0.9,0.9,1.0};
  float down_color[4]={0.1,0.1,0.1,1.0};
  float highlight_color[4]={1.0,0.0,0.0,1.0};
#ifdef pp_GEOMPRINT
  int color_swaps=0;
  int faces_drawn=0;
#endif
  if(nface_normals_single>0){
    glEnable(GL_CULL_FACE);
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glEnable(GL_COLOR_MATERIAL);
    glBegin(GL_TRIANGLES);
    for(j=0;j<nmeshes;j++){
      meshi=meshinfo + j;
      if(meshi->blockvis==0)continue;
      for(i=0;i<meshi->nface_normals_single;i++){
        facei = meshi->face_normals_single[i];
        if(blocklocation==BLOCKlocation_grid){
          vertices = facei->approx_vertex_coords;
        }
        else{
          vertices = facei->exact_vertex_coords;
        }
        showtimelist_handle = facei->showtimelist_handle;
        showtimelist = *showtimelist_handle;
        if(showtimelist!=NULL&&showtimelist[itimes]==0)continue;
        if(showedit_dialog==0){
          new_color=facei->color;
        }
        else{
          if(visNormalEditColors==0)new_color=block_ambient2;
          if(visNormalEditColors==1)new_color=facei->color;
          if(highlight_block==facei->blockageindex&&highlight_mesh==facei->meshindex){
            new_color=highlight_color;
            switch (xyz_dir){
             case XDIR:
              if(facei->dir==UP_X)new_color=up_color;
              if(facei->dir==DOWN_X)new_color=down_color;
              break;
             case YDIR:
              if(facei->dir==UP_Y)new_color=up_color;
              if(facei->dir==DOWN_Y)new_color=down_color;
              break;
             case ZDIR:
              if(facei->dir==UP_Z)new_color=up_color;
              if(facei->dir==DOWN_Z)new_color=down_color;
              break;
             default:
              ASSERT(FFALSE);
              break;
            }
          }
        }
        if(new_color!=old_color){
          old_color=new_color;
          glColor4fv(old_color);
#ifdef pp_GEOMPRINT
          color_swaps++;
#endif
        }
        glNormal3fv(facei->normal);
        glVertex3fv(vertices);
        glVertex3fv(vertices+3);
        glVertex3fv(vertices+6);
        glVertex3fv(vertices);
        glVertex3fv(vertices+6);
        glVertex3fv(vertices+9);
#ifdef pp_GEOMPRINT
        faces_drawn++;
#endif
      }
    }
    glEnd();
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
  }
#ifdef pp_GEOMPRINT
  printf("faces=%i, faces drawn=%i, color switches=%i\n",nface_normals_single,faces_drawn,color_swaps);
#endif
  if(nface_normals_double>0){
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glEnable(GL_COLOR_MATERIAL);
    if(cullfaces==1)glDisable(GL_CULL_FACE);
    glBegin(GL_QUADS);
    for(j=0;j<nmeshes;j++){
      meshi=meshinfo + j;
      for(i=0;i<meshi->nface_normals_double;i++){
        facei = meshi->face_normals_double[i];
        if(blocklocation==BLOCKlocation_grid){
          vertices = facei->approx_vertex_coords;
        }
        else{
          vertices = facei->exact_vertex_coords;
        }
        showtimelist_handle = facei->showtimelist_handle;
        showtimelist = *showtimelist_handle;
        if(showtimelist!=NULL&&showtimelist[itimes]==0)continue;
        if(showedit_dialog==0){
          new_color=facei->color;
        }
        else{
          if(visNormalEditColors==0)new_color=block_ambient2;
          if(visNormalEditColors==1)new_color=facei->color;
          if(highlight_block==facei->blockageindex&&highlight_mesh==facei->meshindex){
            new_color=highlight_color;
            switch (xyz_dir){
             case XDIR:
              if(facei->dir==UP_X)new_color=up_color;
              if(facei->dir==DOWN_X)new_color=down_color;
              break;
             case YDIR:
              if(facei->dir==UP_Y)new_color=up_color;
              if(facei->dir==DOWN_Y)new_color=down_color;
              break;
             case ZDIR:
              if(facei->dir==UP_Z)new_color=up_color;
              if(facei->dir==DOWN_Z)new_color=down_color;
              break;
             default:
              ASSERT(FFALSE);
              break;
            }
          }
        }
        if(new_color!=old_color){
          old_color=new_color;
          glColor4fv(old_color);
        }
        glNormal3fv(facei->normal);
        glVertex3fv(vertices);
        glVertex3fv(vertices+3);
        glVertex3fv(vertices+6);
        glVertex3fv(vertices+9);
      }
    }
    glEnd();
    if(cullfaces==1)glEnable(GL_CULL_FACE);
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
  }
  if(nface_outlines>0){
    glDisable(GL_LIGHTING);
    antialias(1);
    glLineWidth(linewidth);
    glBegin(GL_LINES);
    for(j=0;j<nmeshes;j++){
      meshi = meshinfo + j;
      if(meshi->blockvis==0)continue;
      for(i=0;i<meshi->nface_outlines;i++){
        facei = meshi->face_outlines[i];
        showtimelist_handle = facei->showtimelist_handle;
        showtimelist = *showtimelist_handle;
        if(showtimelist!=NULL&&showtimelist[itimes]==0&&facei->type2==BLOCK_face)continue;
        if(blocklocation==BLOCKlocation_grid){
          vertices = facei->approx_vertex_coords;
        }
        else{
          vertices = facei->exact_vertex_coords;
        }
        if(facei->type2!=OUTLINE_FRAME_face||highlight_flag==1){
          glEnd();
          if(nmeshes>1&&facei->type2==OUTLINE_FRAME_face&&
            highlight_mesh==facei->meshindex&&highlight_flag==1){
            glLineWidth(highlight_linewidth);
          }
          else{
            glLineWidth(*facei->linewidth);
          }
          glBegin(GL_LINES);
          glColor3fv(facei->linecolor);
          glVertex3fv(vertices);
          glVertex3fv(vertices+3);
          glVertex3fv(vertices+3);
          glVertex3fv(vertices+6);
          glVertex3fv(vertices+6);
          glVertex3fv(vertices+9);
          glVertex3fv(vertices+9);
          glVertex3fv(vertices);
          if(showtimelist!=NULL&&showtimelist[itimes]==0){
            glVertex3fv(vertices);
            glVertex3fv(vertices+6);
            glVertex3fv(vertices+3);
            glVertex3fv(vertices+9);
          }
        }
      }
    }
    glEnd();
    antialias(0);
  }
  if(nface_textures>0){
    glEnable(GL_LIGHTING);
    glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
    glEnable(GL_TEXTURE_2D);
    for(j=0;j<nmeshes;j++){
      meshi = meshinfo + j;
      if(meshi->blockvis==0)continue;
      for(i=0;i<meshi->nface_textures;i++){
        facei=meshi->face_textures[i];
        showtimelist_handle = facei->showtimelist_handle;
        showtimelist = *showtimelist_handle;
        if(showtimelist!=NULL&&showtimelist[itimes]==0)continue;
        texti=facei->textureinfo;
        if(blocklocation==BLOCKlocation_grid){
           vertices = facei->approx_vertex_coords;
          tvertices = facei->approx_texture_coords;
        }
        else{
           vertices = facei->exact_vertex_coords;
          tvertices = facei->exact_texture_coords;
        }

        if(facei->type2==BLOCK_face&&cullfaces==0)glDisable(GL_CULL_FACE);


        glBindTexture(GL_TEXTURE_2D,texti->name);
        glBegin(GL_QUADS);

        glNormal3fv(facei->normal);
        glTexCoord2fv(tvertices);
        glVertex3fv(vertices);

        glTexCoord2fv(tvertices+2);
        glVertex3fv(vertices+3);

        glTexCoord2fv(tvertices+4);
        glVertex3fv(vertices+6);

        glTexCoord2fv(tvertices+6);
        glVertex3fv(vertices+9);
        glEnd();
      }
      if(cullfaces==1)glEnable(GL_CULL_FACE);


    }
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_LIGHTING);
  }
}


/* ------------------ initcullgeom ------------------------ */

void initcullgeom(int cullgeomflag){
  culldata *culli;
  int imesh;

  update_initcullgeom=0;
  updatefacelists=1;
  for(imesh=0;imesh<nmeshes;imesh++){
    mesh *meshi;
    int iskip, jskip, kskip;
    int ibeg, iend, jbeg, jend, kbeg, kend;
    float xbeg, xend, ybeg, yend, zbeg, zend;
    int i, j, k;
    int nx, ny, nz;
    int *nxyzgeomcull, *nxyzskipgeomcull;

    meshi=meshinfo+imesh;

    get_cullskips(meshi,cullgeomflag,cullgeom_portsize,&iskip,&jskip,&kskip);
    nx = (meshi->ibar-1)/iskip + 1;
    ny = (meshi->jbar-1)/jskip + 1;
    nz = (meshi->kbar-1)/kskip + 1;
    meshi->ncullgeominfo = nx*ny*nz;

    nxyzgeomcull=meshi->nxyzgeomcull;
    nxyzskipgeomcull=meshi->nxyzskipgeomcull;

    nxyzgeomcull[0]=nx;
    nxyzgeomcull[1]=ny;
    nxyzgeomcull[2]=nz;

    nxyzskipgeomcull[0]=iskip;
    nxyzskipgeomcull[1]=jskip;
    nxyzskipgeomcull[2]=kskip;


    FREEMEMORY(meshi->cullgeominfo);
    NewMemory( (void **)&meshi->cullgeominfo,nx*ny*nz*sizeof(culldata));
    culli=meshi->cullgeominfo;

    for(k=0;k<nz;k++){
      kbeg = k*kskip;
      kend = kbeg + kskip;
      if(kend>meshi->kbar)kend=meshi->kbar;
      zbeg = meshi->zplt[kbeg];
      zend = meshi->zplt[kend];
      for(j=0;j<ny;j++){
        jbeg = j*jskip;
        jend = jbeg + jskip;
        if(jend>meshi->jbar)jend=meshi->jbar;
        ybeg = meshi->yplt[jbeg];
        yend = meshi->yplt[jend];
        for(i=0;i<nx;i++){
          ibeg = i*iskip;
          iend = ibeg + iskip;
          if(iend>meshi->ibar)iend=meshi->ibar;
          xbeg = meshi->xplt[ibeg];
          xend = meshi->xplt[iend];

          culli->ibeg=ibeg;
          culli->iend=iend;

          culli->jbeg=jbeg;
          culli->jend=jend;

          culli->kbeg=kbeg;
          culli->kend=kend;

          culli->xbeg=xbeg;
          culli->xend=xend;

          culli->ybeg=ybeg;
          culli->yend=yend;

          culli->zbeg=zbeg;
          culli->zend=zend;

          culli->iskip=iskip;
          culli->jskip=jskip;
          culli->kskip=kskip;

          culli->npixels=0;
          culli->npixels_old=-1;

          culli++;
        }
      }
    }
  }
}

/* ------------------ get_cullskips ------------------------ */

void get_cullskips(mesh *meshi, int cullflag, int cull_portsize_local, int *iiskip, int *jjskip, int *kkskip){
  int iskip, jskip, kskip;

  if(cullflag==1){
    iskip = cull_portsize_local;
    if(iskip<3)iskip=3;
    if(iskip>meshi->ibar+1)iskip=meshi->ibar+1;

    jskip = cull_portsize_local;
    if(jskip<3)jskip=3;
    if(jskip>meshi->jbar+1)jskip=meshi->jbar+1;

    kskip = cull_portsize_local;
    if(kskip<3)kskip=3;
    if(kskip>meshi->kbar+1)kskip=meshi->kbar+1;
  }
  else{
    iskip = meshi->ibar+1;
    jskip = meshi->jbar+1;
    kskip = meshi->kbar+1;
  }
  *iiskip=iskip;
  *jjskip=jskip;
  *kkskip=kskip;
}

/* ------------------ get_face_port ------------------------ */

culldata *get_face_port(mesh *meshi, facedata *facei){
  int ii1, jj1, kk1;
  int ii2, jj2, kk2;
  int nx, ny, nz;
  int ixyz;
  culldata *return_cull;
  int *skip2,*nxyz;

  skip2=meshi->nxyzskipgeomcull;
  nxyz=meshi->nxyzgeomcull;
  nx=nxyz[0];
  ny=nxyz[1];
  nz=nxyz[2];

  ii1=facei->imin/skip2[0];
  if(facei->imax!=facei->imin){
    ii2=(facei->imax-1)/skip2[0];
    if(ii1!=ii2)return NULL;
  }
  if(ii1<0||ii1>nx)return NULL;

  jj1=facei->jmin/skip2[1];
  if(facei->jmin!=facei->jmax){
    jj2=(facei->jmax-1)/skip2[1];
    if(jj1!=jj2)return NULL;
  }
  if(jj1<0||jj1>ny)return NULL;

  kk1=facei->kmin/skip2[2];
  if(facei->kmin!=facei->kmax){
    kk2=(facei->kmax-1)/skip2[2];
    if(kk1!=kk2)return NULL;
  }
  if(kk1<0||kk1>nz)return NULL;


  ixyz = ii1 + jj1*nx + kk1*nx*ny;
  return_cull = meshi->cullgeominfo + ixyz;
  
  return return_cull;
}

/* ------------------ blockagecompare ------------------------ */

int blockagecompare( const void *arg1, const void *arg2 ){
  blockagedata *bc1, *bc2;
  int *ijk1, *ijk2;

  bc1 = *(blockagedata **)arg1;
  bc2 = *(blockagedata **)arg2;

  ijk1 = bc1->ijk;
  ijk2 = bc2->ijk;

  if(ijk1[0]<ijk2[0])return -1;
  if(ijk1[0]>ijk2[0])return 1;
  if(ijk1[1]<ijk2[1])return -1;
  if(ijk1[1]>ijk2[1])return 1;
  if(ijk1[2]<ijk2[2])return -1;
  if(ijk1[2]>ijk2[2])return 1;
  if(ijk1[3]<ijk2[3])return -1;
  if(ijk1[3]>ijk2[3])return 1;
  if(ijk1[4]<ijk2[4])return -1;
  if(ijk1[4]>ijk2[4])return 1;
  if(ijk1[5]<ijk2[5])return -1;
  if(ijk1[5]>ijk2[5])return 1;
  return 0;
}

/* ------------------ remove_dup_blockages ------------------------ */

void remove_dup_blockages(void){
  int i;

  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi = meshinfo + i;

    if(meshi->nbptrs>1){
      blockagedata **bclist;
      int jj,j;

      bclist=meshi->blockageinfoptrs;
      qsort(bclist,(size_t)meshi->nbptrs,sizeof(blockagedata *),blockagecompare);
      for(j=1;j<meshi->nbptrs;j++){
        blockagedata *bc, *bcm1;
        int *ijk1, *ijk2;
      
        bc = bclist[j];
        bcm1 = bclist[j-1];
        if(bc->nshowtime>0)continue;
        if(bcm1->nshowtime>0)continue;
        ijk1=bcm1->ijk;
        ijk2=bc->ijk;
        if(ijk1[1]-ijk1[0]>1)continue; // only consider removing one cell blockages
        if(ijk1[3]-ijk1[2]>1)continue;
        if(ijk1[5]-ijk1[4]>1)continue;
        if(ijk1[0]!=ijk2[0]||ijk1[1]!=ijk2[1])continue;
        if(ijk1[2]!=ijk2[2]||ijk1[3]!=ijk2[3])continue;
        if(ijk1[4]!=ijk2[4]||ijk1[5]!=ijk2[5])continue;
        bcm1->dup=1;
        bc->dup=2;
      }
      jj=0;
      for(j=0;j<meshi->nbptrs;j++){
        blockagedata *bc;
      
        bc=bclist[j];
        if(bc->dup==1)continue;
        bclist[jj++]=bc;
      }
      meshi->nbptrs=jj;
    }
  }
  updatefacelists=1;
}

/* ------------------ getlabels ------------------------ */

void getobstlabels(const char *filein){

  FILE *stream_in;
  char buffer[1000];
  int fdsobstcount=0;
  int i,j;
  char *obstlabel;
  mesh *meshi;
  blockagedata *bc;
  int id;
  size_t lenlabel;
  char **obstlabels=NULL;
  int nobstlabels=0;

  if(filein==NULL)return;
  stream_in = fopen(filein,"r");
  if(stream_in==NULL)return;

  while(!feof(stream_in)){
    if(fgets(buffer,1000,stream_in)==NULL)break;

    if(STRSTR(buffer,"&OBST")==NULL)continue;
    fdsobstcount++;
  }
  nobstlabels=fdsobstcount;
  if(nobstlabels>0){
    NewMemory((void **)&obstlabels,nobstlabels*sizeof(char *));
    for(i=0;i<nobstlabels;i++){
      obstlabels[i]=NULL;
    }
  }
  rewind(stream_in);
  fdsobstcount=0;
  while(!feof(stream_in)){
    if(fgets(buffer,1000,stream_in)==NULL)break;

    if(STRSTR(buffer,"&OBST")==NULL)continue;
    fdsobstcount++;
    while((obstlabel=strstr(buffer,"/"))==NULL){
      fgets(buffer,1000,stream_in);
    }
    obstlabel++;
    lenlabel=strlen(obstlabel);
    obstlabel=trim_front(obstlabel);
    trim(obstlabel);
    lenlabel=strlen(obstlabel);
    if(lenlabel>0){
      NewMemory((void **)&obstlabels[fdsobstcount-1],(unsigned int)(lenlabel+1));
      strcpy(obstlabels[fdsobstcount-1],obstlabel);
    }
  }
  fclose(stream_in);

  for(i=0;i<nmeshes;i++){
    meshi = meshinfo + i;
    for(j=0;j<meshi->nbptrs;j++){
      bc = meshi->blockageinfoptrs[j];
      id = bc->id-1;
      if(id>=0&&id<nobstlabels){
        if(obstlabels[id]!=NULL){
          lenlabel=strlen(obstlabels[id]);
          ResizeMemory((void **)&bc->label,(unsigned int)(lenlabel+1));
          strcpy(bc->label,obstlabels[id]);
        }
      }
    }
  }
  for(i=0;i<nobstlabels;i++){
    FREEMEMORY(obstlabels[i]);
  }
  FREEMEMORY(obstlabels);
}
