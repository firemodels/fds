#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include GLUT_H

#include "update.h"
#include "smv_endian.h"
#include "smokeviewvars.h"

float get_z_terrain(float x, float y);
void init_tnode(terraindata *terri);
void init_tnorm(terraindata *terri);
void init_terraincell(terraindata *terri);
void free_terraincell(terraindata *terri);

#define ijnode2(i,j) ((nxcell+1)*(j) + (i))
#define ijnode3(i,j) ((nycell+1)*(i) + (j))
#define FORTWUIREAD(var,size) FSEEK(WUIFILE,4,SEEK_CUR);\
                           returncode=fread(var,4,size,WUIFILE);\
                           if(endianswitch==1)endian_switch(var,size);\
                           FSEEK(WUIFILE,4,SEEK_CUR)

float *get_terraincolor(terraincell *ti);
int getterrain_data(char *file,terraindata *terri);
int getterrain_size(char *file,float *xmin, float *xmax, int *nx, float *ymin, float *ymax, int *ny, int *ntimes);
void drawcone(float d1, float height, float *rgbcolor);
void drawdisk(float diameter, float height, float *rgbcolor);


/* ------------------ drawnorth ------------------------ */

void drawnorth(void){
  glPushMatrix();
  glTranslatef(northangle_position[0], northangle_position[1], northangle_position[2]);
  glRotatef(-northangle, 0.0, 0.0, 1.0);
  glBegin(GL_LINES);
  glColor3fv(foregroundcolor);
  glVertex3f(0.0, 0.0, 0.0);
  glVertex3f(0.0, 0.1, 0.0);
  glEnd();
  glPopMatrix();
}

  /*
typedef struct {
  float xyz[3];
  float trunk_diam;
  float tree_height;
  float base_diam;
  float base_height;
*/

/* ------------------ drawtrees ------------------------ */

void drawtrees(void){
  int i;

 glEnable(GL_LIGHTING);

 glEnable(GL_COLOR_MATERIAL);

  glPushMatrix();
  glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
  glTranslatef(-xbar0,-ybar0,-zbar0);
  for(i=0;i<ntreeinfo;i++){
    treedata *treei;
    float crown_height;
    int state;

    treei = treeinfo + i;

    state=0;
    if(showtime==1&&global_times!=NULL){
      ASSERT(itimes>=0)
      if(treei->time_char>0.0&&global_times[itimes]>treei->time_char)state=1;
      if(treei->time_complete>0.0&&global_times[itimes]>treei->time_complete)state=2;
    }

    glPushMatrix();
    glTranslatef(treei->xyz[0],treei->xyz[1],treei->xyz[2]);

    switch(state){
      case 0:
        glColor4fv(trunccolor);
        drawdisk(treei->trunk_diam,treei->base_height,trunccolor);

        crown_height=treei->tree_height-treei->base_height;
        glTranslatef(0.0,0.0,treei->base_height);
        glColor4fv(treecolor);
        drawcone(treei->base_diam,crown_height,treecolor);
        break;
      case 1:
        glColor4fv(treecharcolor);
        drawdisk(treei->trunk_diam,treei->base_height,trunccolor);

        crown_height=treei->tree_height-treei->base_height;
        glTranslatef(0.0,0.0,treei->base_height);
        drawcone(treei->base_diam,crown_height,treecolor);
        break;
      case 2:
        glColor4fv(treecharcolor);
        drawdisk(treei->trunk_diam,treei->base_height,trunccolor);
        crown_height=treei->tree_height-treei->base_height;
        glTranslatef(0.0,0.0,treei->base_height);
        drawcone(treei->trunk_diam,crown_height,trunccolor);
        break;
      default:
        ASSERT(FFALSE);
        break;
    }
    glPopMatrix();


  }
  glPopMatrix();


  glDisable(GL_COLOR_MATERIAL);

}

/* ------------------ get_zcell_val ------------------------ */

float get_zcell_val(meshdata *meshi,float xval, float yval, float *zval_offset, int *loc){
  int imesh;
  int meshstart=-1;

  if(meshi==NULL)meshstart=0;
  if(zval_offset!=NULL)*zval_offset=0.0;
  for(imesh=meshstart;imesh<nmeshes;imesh++){
    meshdata *meshj;
    float *xplt, *yplt;
    int ibar, jbar;

    if(imesh==-1){
      meshj=meshi;
    }
    else{
      meshj=meshinfo+imesh;
      if(meshi==meshj)continue;
    }
    xplt = meshj->xplt_orig;
    yplt = meshj->yplt_orig;
    ibar = meshj->ibar;
    jbar = meshj->jbar;
    if(xplt[0]<=xval&&xval<=xplt[ibar]&&yplt[0]<=yval&&yval<=yplt[jbar]){
      float dx, dy;
      terraindata *terri;
      int ival, jval;
      int nxcell;
      float zval;
      float *zcell;
      float zval_return;

      dx = xplt[1]-xplt[0];
      dy = yplt[1]-yplt[0];
      ival = (xval-xplt[0])/dx;
      if(ival>=ibar)ival=ibar-1;
      jval = (yval-yplt[0])/dy;
      if(jval>=jbar)jval=jbar-1;
      terri=meshj->terrain;
      nxcell = terri->nx;
      zcell = terri->zcell;
      zval = zcell[IJCELL2(ival,jval)];
      *loc=1;
      zval_return = zterrain_min+vertical_factor*(zval-zterrain_min);
      if(zval_offset!=NULL){
        *zval_offset=(vertical_factor-1.0)*(zval-zterrain_min);
      }
      return zval_return;
    }
  }
  *loc=0;
  return 0.0;
}

/* ------------------ get_zcell_val_offset ------------------------ */

float get_zcell_val_offset(meshdata *meshi,float xval, float yval, int *loc){
  int imesh;
  int meshstart=-1;

  if(meshi==NULL)meshstart=0;

  for(imesh=meshstart;imesh<nmeshes;imesh++){
    meshdata *meshj;
    float *xplt, *yplt;
    int ibar, jbar;

    if(imesh==-1){
      meshj=meshi;
    }
    else{
      meshj=meshinfo+imesh;
      if(meshi==meshj)continue;
    }

// convert xval and yval to "user" units
    xval = DENORMALIZE_X(xval);
    yval = DENORMALIZE_Y(yval);

    xplt = meshj->xplt_orig;
    yplt = meshj->yplt_orig;
    ibar = meshj->ibar;
    jbar = meshj->jbar;
    if(xplt[0]<=xval&&xval<=xplt[ibar]&&yplt[0]<=yval&&yval<=yplt[jbar]){
      terraindata *terri;
      int ival, jval;
      float dx, dy;
      int nxcell;
      float *zcell,zval;
      float zvaloffset;

      dx = xplt[1]-xplt[0];
      dy = yplt[1]-yplt[0];
      ival = (xval-xplt[0])/dx;
      if(ival>=ibar)ival=ibar-1;
      jval = (yval-yplt[0])/dy;
      if(jval>=jbar)jval=jbar-1;
      terri=meshj->terrain;
      nxcell = terri->nx;
      zcell = terri->zcell;
      zval = zcell[IJCELL2(ival,jval)];
      *loc=1;
      zvaloffset = zterrain_min+vertical_factor*(zval-zterrain_min)-zval;

    // convert zoffset back to smokeview/scaled units

      zvaloffset = SCALE2SMV(zvaloffset);
      return zvaloffset;
    }
  }
  *loc=0;
  return 0.0;
}

/* ------------------ update_terrain_colors ------------------------ */

void update_terrain_colors(void){
  int i;

  for(i=0;i<MAXRGB;i++){
    float f1;

    f1 = (float)i/(float)(MAXRGB-1);
    rgbterrain[4*i  ]=((1.0-f1)*terrain_rgba_zmin[0] + f1*terrain_rgba_zmax[0])/255.0;
    rgbterrain[4*i+1]=((1.0-f1)*terrain_rgba_zmin[1] + f1*terrain_rgba_zmax[1])/255.0;
    rgbterrain[4*i+2]=((1.0-f1)*terrain_rgba_zmin[2] + f1*terrain_rgba_zmax[2])/255.0;
    rgbterrain[4*i+3]=1.0;
  }
}

/* ------------------ terrain2geom ------------------------ */

void terrain2geom(float xmin, float xmax, float ymin, float ymax, int nx, int ny, float (*comp_func)(float, float)){
  int i, j;
  float dx, dy;
  float x, y, z;
  int nverts, nfaces;
  float *verts;
  int *faces;
  int ivert, iface;

  if(nx<2||ny<2)return;
  dx = (xmax-xmin)/(nx-1);
  dy = (ymax-ymin)/(ny-1);

  nverts = nx*ny;
  nfaces = 2*(nx-1)*(ny-1);

  NewMemory((void **)&verts,3*nverts*sizeof(float));
  NewMemory((void **)&faces,3*nfaces*sizeof(int));

  ivert=0;
  for(j=0;j<ny;j++){
    y = ymin + j*dy;
    for(i=0;i<nx;i++){
      x = xmin + i*dx;
      z = comp_func(x,y);
      verts[ivert++]=x;
      verts[ivert++]=y;
      verts[ivert++]=z;
    }
  }
#define IJgeom(i,j) (1+(i) + (j)*nx)

  iface=0;
  for(j=0;j<ny-1;j++){
    for(i=0;i<nx-1;i++){
      int i11, i12, i22, i21;

      i11 = IJgeom(i,j);
      i12 = IJgeom(i,j+1);
      i22 = IJgeom(i+1,j+1);
      i21 = IJgeom(i+1,j);
      faces[iface++]=i11;
      faces[iface++]=i21;
      faces[iface++]=i22;
      faces[iface++]=i11;
      faces[iface++]=i22;
      faces[iface++]=i12;
    }
  }
  FORTgeomout(verts,&nverts,faces,&nfaces);
  FREEMEMORY(verts);
  FREEMEMORY(faces);
}

/* ------------------ initterrain_all ------------------------ */

void initterrain_all(void){
  int imesh;
  float zmin, zmax, dz;

  for(imesh=0;imesh<nmeshes;imesh++){
    meshdata *meshi;
    terraindata *terri;
    int j;
    float dx, dy;
    float *znode, *znode_offset;
    int nycell;
    unsigned char *uc_znormal;

    meshi = meshinfo + imesh;

    terri = meshi->terrain;

    dx = terri->x[1]-terri->x[0];
    dy = terri->y[1]-terri->y[0];

    znode = terri->znode;
    znode_offset = terri->znode_offset;
    nycell = terri->ny;

    uc_znormal = terri->uc_znormal;
    for(j=0;j<=terri->ny;j++){
      int i;
      float ynode;

      ynode = terri->y[j];

      for(i=0;i<=terri->nx;i++){
        float xnode;
        int count, loc1, loc2, loc3, loc4;
        float val1, val2, val3, val4;
        float val1_offset, val2_offset, val3_offset, val4_offset;
        float valx1a, valx2a, valx3a, valx4a;
        float valx1b, valx2b, valx3b, valx4b;
        float valy1a, valy2a, valy3a, valy4a;
        float valy1b, valy2b, valy3b, valy4b;
        float zval, zval_offset;
        float zvalxa, zvalxb;
        float zvalya, zvalyb;
        float dxa, dxb, dya, dyb;
        float dzdx, dzdy;
        float sum;
        float znormal3[3];
        float denom;

        xnode = terri->x[i];

        val1 =  get_zcell_val(meshi,xnode-dx/2.0,ynode-dy/2.0,&val1_offset,&loc1);
        val2 =  get_zcell_val(meshi,xnode+dx/2.0,ynode-dy/2.0,&val2_offset,&loc2);
        val3 =  get_zcell_val(meshi,xnode+dx/2.0,ynode+dy/2.0,&val3_offset,&loc3);
        val4 =  get_zcell_val(meshi,xnode-dx/2.0,ynode+dy/2.0,&val4_offset,&loc4);
        count = loc1 + loc2 + loc3 + loc4;

        zval = (val1*loc1 + val2*loc2 + val3*loc3 + val4*loc4)/(float)MAX(1,count);

        znode[ijnode3(i,j)]=zval;
        zval_offset = SCALE2SMV((val1_offset*loc1 + val2_offset*loc2 + val3_offset*loc3 + val4_offset*loc4))/(float)MAX(1,count);

        *znode_offset++=zval_offset;

 // compute (f(x+dx,y) - f(x-dx,y))/(2*dx)

        valx1a =  get_zcell_val(meshi,xnode-dx-dx/2.0,ynode-dy/2.0,NULL,&loc1);
        valx2a =  get_zcell_val(meshi,xnode-dx+dx/2.0,ynode-dy/2.0,NULL,&loc2);
        valx3a =  get_zcell_val(meshi,xnode-dx+dx/2.0,ynode+dy/2.0,NULL,&loc3);
        valx4a =  get_zcell_val(meshi,xnode-dx-dx/2.0,ynode+dy/2.0,NULL,&loc4);
        count = loc1 + loc2 + loc3 + loc4;
        zvalxa = valx1a*loc1 + valx2a*loc2 + valx3a*loc3 + valx4a*loc4;
        if(count==0){
          zvalxa = zval;
          dxa = 0.0;
        }
        else{
          zvalxa /= (float)count;
          dxa = dx;
        }
        valx1b =  get_zcell_val(meshi,xnode+dx-dx/2.0,ynode-dy/2.0,NULL,&loc1);
        valx2b =  get_zcell_val(meshi,xnode+dx+dx/2.0,ynode-dy/2.0,NULL,&loc2);
        valx3b =  get_zcell_val(meshi,xnode+dx+dx/2.0,ynode+dy/2.0,NULL,&loc3);
        valx4b =  get_zcell_val(meshi,xnode+dx-dx/2.0,ynode+dy/2.0,NULL,&loc4);
        count = loc1 + loc2 + loc3 + loc4;
        zvalxb = valx1b*loc1 + valx2b*loc2 + valx3b*loc3 + valx4b*loc4;
        if(count==0){
          zvalxb = zval;
          dxb = 0.0;
        }
        else{
          zvalxb /= (float)count;
          dxb = dx;
        }
        denom = dxa+dxb;
        if(denom==0.0){
          dzdx=1.0;
        }
        else{
          dzdx = (zvalxb - zvalxa)/denom;
        }

 // compute (f(x,y+dy) - f(x,y-dy))/(2*dy)

        valy1a =  get_zcell_val(meshi,xnode-dx/2.0,ynode-dy-dy/2.0,NULL,&loc1);
        valy2a =  get_zcell_val(meshi,xnode+dx/2.0,ynode-dy-dy/2.0,NULL,&loc2);
        valy3a =  get_zcell_val(meshi,xnode+dx/2.0,ynode-dy+dy/2.0,NULL,&loc3);
        valy4a =  get_zcell_val(meshi,xnode-dx/2.0,ynode-dy+dy/2.0,NULL,&loc4);
        count = loc1 + loc2 + loc3 + loc4;
        zvalya = valy1a*loc1 + valy2a*loc2 + valy3a*loc3 + valy4a*loc4;
        if(count==0){
          zvalya = zval;
          dya = 0.0;
        }
        else{
          zvalya /= (float)count;
          dya = dy;
        }
        valy1b =  get_zcell_val(meshi,xnode-dx/2.0,ynode+dy-dy/2.0,NULL,&loc1);
        valy2b =  get_zcell_val(meshi,xnode+dx/2.0,ynode+dy-dy/2.0,NULL,&loc2);
        valy3b =  get_zcell_val(meshi,xnode+dx/2.0,ynode+dy+dy/2.0,NULL,&loc3);
        valy4b =  get_zcell_val(meshi,xnode-dx/2.0,ynode+dy+dy/2.0,NULL,&loc4);
        count = loc1 + loc2 + loc3 + loc4;
        zvalyb = valy1b*loc1 + valy2b*loc2 + valy3b*loc3 + valy4b*loc4;
        if(count==0){
          zvalyb = zval;
          dyb = 0.0;
        }
        else{
          zvalyb /= (float)count;
          dyb = dy;
        }
        denom = dya + dyb;
        if(denom==0.0){
          dzdy=1.0;
        }
        else{
          dzdy = (zvalyb - zvalya)/denom;
        }

     //     i  j  k
     //     1  0 dzdx
     //     0  1 dzdy

     //     -dzdx -dzdy 1

        //znormal = terri->znormal + 3*ijnode2(i,j);
        uc_znormal = terri->uc_znormal + ijnode3(i,j);
        znormal3[0] = -dzdx;
        znormal3[1] = -dzdy;
        znormal3[2] = 1.0;

        sum  = znormal3[0]*znormal3[0];
        sum += znormal3[1]*znormal3[1];
        sum += znormal3[2]*znormal3[2];
        sum = sqrt(sum);
        znormal3[0]/=sum;
        znormal3[1]/=sum;
        znormal3[2]/=sum;
        *uc_znormal = getnormalindex(wui_sphereinfo, znormal3);
      }
    }
  }

  zmin = meshinfo->terrain->znode[0];
  zmax = zmin;
  for(imesh=0;imesh<nmeshes;imesh++){
    meshdata *meshi;
    terraindata *terri;
    int i;

    meshi = meshinfo + imesh;
    terri = meshi->terrain;

    for(i=0;i<(terri->nx+1)*(terri->ny+1);i++){
      float *znode;

      znode = terri->znode+i;
      zmin = MIN(zmin,*znode);
      zmax = MAX(zmax,*znode);
    }
  }
  dz = (zmax - zmin)/12.0;
  for(imesh=0;imesh<nmeshes;imesh++){
    meshdata *meshi;
    terraindata *terri;
    int i;

    meshi = meshinfo + imesh;
    terri = meshi->terrain;
    for(i=0;i<13;i++){
      terri->levels[i]=zmin + i*dz;
    }
    terri->levels[12]=zmax;

    freecontour(&meshi->terrain_contour);
    initcontour(&meshi->terrain_contour,rgbptr,nrgb);

    meshi->terrain_contour.idir=3;
    meshi->terrain_contour.xyzval=zmin;

    getcontours(meshi->xplt_orig,meshi->yplt_orig,terri->nx+1,terri->ny+1,
      terri->znode, NULL, terri->levels,DONT_GET_AREAS,DATA_FORTRAN,
      &meshi->terrain_contour);

  }
  {
    int nx, ny;

    nx = (xbarORIG-xbar0ORIG)/(meshinfo[0].xplt_orig[1]-meshinfo[0].xplt_orig[0])+1;
    ny = (ybarORIG-ybar0ORIG)/(meshinfo[0].yplt_orig[1]-meshinfo[0].yplt_orig[0])+1;

    terrain2geom(xbar0ORIG, xbarORIG, ybar0ORIG, ybarORIG, nx, ny,get_z_terrain);
  }
}

/* ------------------ initterrain_znode ------------------------ */

void initterrain_znode(meshdata *meshi, terraindata *terri, float xmin, float xmax, int nx, float ymin, float ymax, int ny,
                       int allocate_memory){
  float dx, dy;
  float *x, *y, *z;
  int nxcell;
  int i,j;
  int ij;

  if(meshi!=NULL){
    meshi->terrain=terri;
  }

  if(allocate_memory==1){
    terri->display=0;
    terri->loaded=0;
    terri->autoload=0;
    terri->ntimes=0;
    terri->x=NULL;
    terri->y=NULL;
    terri->times=NULL;
    terri->zcell=NULL;
    terri->znode=NULL;
    terri->uc_znormal=NULL;
    terri->tcell=NULL;
    terri->ter_texture=NULL;
    terri->state=NULL;
    terri->timeslist=NULL;
  }

  terri->xmin=xmin;
  terri->xmax=xmax;
  terri->ymin=ymin;
  terri->ymax=ymax;
  terri->ny=ny;
  if(nx<0){
    nx=-nx;
  }
  terri->nx=nx;

  if(allocate_memory==1){
    NewMemory((void **)&terri->x,(nx+1)*sizeof(float));
    NewMemory((void **)&terri->y,(ny+1)*sizeof(float));
    NewMemory((void **)&terri->zcell,nx*ny*sizeof(float));
    NewMemory((void **)&terri->state,nx*ny);
    NewMemory((void **)&terri->znode,(nx+1)*(ny+1)*sizeof(float));
    NewMemory((void **)&terri->znode_offset,(nx+1)*(ny+1)*sizeof(float));
    NewMemory((void **)&terri->znode_scaled,(nx+1)*(ny+1)*sizeof(float));
    NewMemory((void **)&terri->uc_znormal,(nx+1)*(ny+1)*sizeof(unsigned char));
  }

  x = terri->x;
  y = terri->y;
  dx = (xmax-xmin)/nx;
  dy = (ymax-ymin)/ny;
  for(i=0;i<nx;i++){
    x[i] = xmin + dx*i;
  }
  x[nx] = xmax;

  for(i=0;i<ny;i++){
    y[i] = ymin + dy*i;
  }
  y[ny] = ymax;

  z=terri->zcell;

  nxcell = nx;
  for(j=0;j<ny;j++){
    for(i=0;i<nx;i++){
      float zval;

      ij = IJCELL2(i,j);
      zval=meshi->zcell[ij];
      zterrain_min = MIN(zval,zterrain_min);
      zterrain_max = MAX(zval,zterrain_max);
      z[ij]=zval;
    }
  }
}

/* ------------------ drawterrain ------------------------ */

void drawterrain(terraindata *terri, int only_geom){
  float *znode, *zn;
  unsigned char *uc_znormal;
  int nycell;
  int i, j;
  float *x, *y;
  terraincell *ti;
  float terrain_color[4];
  float terrain_shininess=100.0;
  float terrain_specular[4]={0.8,0.8,0.8,1.0};
  //float zt_min, zt_max;

#define ZOFFSET 0.001

//xxx is_bottom is not defined correctly
//    comment out following line to show all meshes
  if(terri->terrain_mesh->is_bottom==0)return;

  terrain_color[0]=0.47843;
  terrain_color[1]=0.45882;
  terrain_color[2]=0.18824;
  terrain_color[3]=1.0;

  glPushMatrix();
  glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
  glTranslatef(-xbar0,-ybar0,-zbar0);

  glEnable(GL_LIGHTING);
  glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&terrain_shininess);
  glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,rgbterrain);
  glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,terrain_specular);
  glEnable(GL_COLOR_MATERIAL);

  glBegin(GL_QUADS);
  uc_znormal = terri->uc_znormal;
  znode = terri->znode;
  nycell = terri->ny;
  x = terri->x;
  y = terri->y;
  ti = terri->tcell;
  glColor4fv(terrain_color);
  for(j=0;j<terri->ny;j++){
    int jp1;

    jp1 = j + 1;

    for(i=0;i<terri->nx;i++){
      unsigned char *uc_zn;
      int ip1;
      float *ter_rgbptr;
      float zval;

      ip1 = i + 1;
      if(only_geom==0){
        ter_rgbptr = get_terraincolor(ti);
        glColor4fv(ter_rgbptr);
      }
      uc_zn = uc_znormal+ijnode3(i,j);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));

      glNormal3fv(zn);
      zval = znode[ijnode3(i,j)]+ZOFFSET;
      glVertex3f(x[i],y[j],zval);

      uc_zn = uc_znormal+ijnode3(ip1,j);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));
      glNormal3fv(zn);
      zval = znode[ijnode3(ip1,j)]+ZOFFSET;
      glVertex3f(x[i+1],y[j],zval);

      uc_zn = uc_znormal+ijnode3(ip1,jp1);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));
      glNormal3fv(zn);
      zval = znode[ijnode3(ip1,jp1)]+ZOFFSET;
      glVertex3f(x[i+1],y[j+1],zval);

      uc_zn = uc_znormal+ijnode3(i,jp1);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));
      glNormal3fv(zn);
      zval = znode[ijnode3(i,jp1)]+ZOFFSET;
      glVertex3f(x[i],y[j+1],zval);

      ti++;
    }
  }
  glEnd();

  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);

  glPopMatrix();

}

/* ------------------ drawterrain_texture ------------------------ */

void drawterrain_texture(terraindata *terri, int only_geom){
  float *znode;
  unsigned char *uc_znormal, *uc_zn;
  int nxcell,nycell;
  int i, j;
  float *x, *y;
  float terrain_color[4];

  if(terri->terrain_mesh->is_bottom==0)return;

  terrain_color[0]=1.0;
  terrain_color[1]=1.0;
  terrain_color[2]=1.0;
  terrain_color[3]=1.0;

  glPushMatrix();
  glScalef(SCALE2SMV(mscale[0]),SCALE2SMV(mscale[1]),SCALE2SMV(mscale[2]));
  glTranslatef(-xbar0,-ybar0,-zbar0);

  glEnable(GL_LIGHTING);
  glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D,terrain_texture->name);

  glEnable(GL_COLOR_MATERIAL);
  glColor4fv(terrain_color);
  glBegin(GL_QUADS);
  uc_znormal = terri->uc_znormal;
  znode = terri->znode;
  nxcell = terri->nx;
  nycell = terri->ny;
  x = terri->x;
  y = terri->y;
  for(j=0;j<terri->ny;j++){
    int jp1;
    float ty,typ1;

    jp1 = j + 1;
    ty = (y[j]-ybar0ORIG)/(ybarORIG-ybar0ORIG);
    typ1 = (y[j+1]-ybar0ORIG)/(ybarORIG-ybar0ORIG);

    for(i=0;i<terri->nx;i++){
      float *zn;
      int ip1;
      float tx,txp1;

      ip1 = i + 1;
      tx = (x[i]-xbar0ORIG)/(xbarORIG-xbar0ORIG);
      txp1 = (x[i+1]-xbar0ORIG)/(xbarORIG-xbar0ORIG);

      uc_zn = uc_znormal+ijnode2(i,j);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));
      glNormal3fv(zn);
      glTexCoord2f(tx,ty);
      glVertex3f(x[i],y[j],znode[ijnode3(i,j)]);

      uc_zn = uc_znormal+ijnode2(ip1,j);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));
      glNormal3fv(zn);
      glTexCoord2f(txp1,ty);
      glVertex3f(x[i+1],y[j],znode[ijnode3(ip1,j)]);

      uc_zn = uc_znormal+ijnode2(ip1,jp1);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));
      glNormal3fv(zn);
      glTexCoord2f(txp1,typ1);
      glVertex3f(x[i+1],y[j+1],znode[ijnode3(ip1,jp1)]);

      uc_zn = uc_znormal+ijnode2(i,jp1);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));
      glNormal3fv(zn);
      glTexCoord2f(tx,typ1);
      glVertex3f(x[i],y[j+1],znode[ijnode3(i,jp1)]);
    }
  }
  glEnd();

  glDisable(GL_TEXTURE_2D);

  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);

  glPopMatrix();

}
/* ------------------ get_terraincolor ------------------------ */

float *get_terraincolor(terraincell *ti){
  int index;
  int i, ileft;
  float sv_time;
  float *ter_time;
  float wuicolor[4]={1.0,0.0,0.0,1.0};

  if(ti==NULL)return getcolorptr(wuicolor);

  if(global_times==NULL||ti->time==NULL){
    index = ti->state[0]%10;
    return rgb_terrain[index];
  }

  sv_time = global_times[itimes];
  ter_time = ti->time;
  ileft = ti->interval;

  if(ter_time[ileft]<=sv_time&&sv_time<ter_time[ileft+1]){
    return rgb_terrain[ileft%10];
  }

  for(i=ileft+1;i<ti->nstates-1;i++){
    if(ter_time[i]<=sv_time&&sv_time<ter_time[i+1]){
      ti->interval=i;
      return rgb_terrain[i%10];
    }
  }
  if(sv_time>=ter_time[ti->nstates-1]){
    ileft = ti->nstates-1;
    ti->interval=ileft;
    return rgb_terrain[ileft%10];
  }
  ileft = 0;
  ti->interval=ileft;
  return rgb_terrain[ileft%10];
}

/* ------------------ readterrain ------------------------ */

void readterrain(char *file, int ifile, int flag, int *errorcode){
  terraindata *terri=NULL;
  float xmin, xmax;
  int nx;
  float ymin, ymax;
  int ny;
  float *x, *y;
  float dx, dy;
  int i;

  if(ifile<0||ifile>=nterraininfo)return;
  terri = terraininfo + ifile;

  if(flag==UNLOAD){
    FREEMEMORY(terri->x);
    FREEMEMORY(terri->y);
    FREEMEMORY(terri->zcell);
    FREEMEMORY(terri->znode);
    FREEMEMORY(terri->uc_znormal);
    FREEMEMORY(terri->times);
    free_terraincell(terri);
    terri->loaded=0;
    terri->display=0;
    Update_Times();
    return;
  }

  if(getterrain_size(file,&xmin, &xmax, &nx, &ymin, &ymax, &ny, &nglobal_times)!=0)return;

  terri->xmin = xmin;
  terri->xmax = xmax;
  terri->nx = nx;
  terri->ymin = ymin;
  terri->ymax = ymax;
  terri->ny = ny;
  terri->ntimes=nglobal_times;

  NewMemory((void **)&terri->times,nglobal_times*sizeof(float));
  NewMemory((void **)&terri->x,(nx+1)*sizeof(float));
  NewMemory((void **)&terri->y,(ny+1)*sizeof(float));
  NewMemory((void **)&terri->zcell,nx*ny*sizeof(float));
  NewMemory((void **)&terri->state,nx*ny);
  NewMemory((void **)&terri->znode,(nx+1)*(ny+1)*sizeof(float));
//  NewMemory((void **)&terri->znormal,3*(nx+1)*(ny+1)*sizeof(float));
  init_terraincell(terri);

  x = terri->x;
  y = terri->y;
  dx = (xmax-xmin)/nx;
  dy = (ymax-ymin)/ny;
  for(i=0;i<nx;i++){
    x[i] = xmin + dx*i;
  }
  x[nx] = xmax;
  for(i=0;i<ny;i++){
    y[i] = ymin + dy*i;
  }
  y[ny] = ymax;

  if(getterrain_data(file,terri)!=0){
    readterrain("",ifile,UNLOAD,errorcode);
    return;
  }
  terri->loaded=1;
  visTerrainType=TERRAIN_3D;
  plotstate=getplotstate(DYNAMIC_PLOTS);
  Update_Times();
#ifdef pp_MEMPRINT
  PRINTF("After terrain file load: \n");
  PrintMemoryInfo;
#endif
  Idle_CB();
  glutPostRedisplay();
}

/* ------------------ getterrain_size ------------------------ */

int getterrain_size(char *file,float *xmin, float *xmax, int *nx, float *ymin, float *ymax, int *ny, int *times_local){
  FILE *WUIFILE;
  int one;
  float xyminmax[4];
  int nxy[2];
  int endianswitch=0;
  size_t returncode;
  int version;
  float time_local;
  int nchanges;
  int nt=0;

  WUIFILE = fopen(file,"rb");
  if(WUIFILE==NULL)return 1;

  FSEEK(WUIFILE,4,SEEK_CUR);fread(&one,4,1,WUIFILE);FSEEK(WUIFILE,4,SEEK_CUR);
  if(one!=1)endianswitch=1;

  FORTWUIREAD(&version,1);
  FORTWUIREAD(xyminmax,4);
  *xmin=xyminmax[0];
  *xmax=xyminmax[1];
  *ymin=xyminmax[2];
  *ymax=xyminmax[3];

  FORTWUIREAD(nxy,2);
  *nx=nxy[0];
  *ny=nxy[1];

  FSEEK(WUIFILE,16+5*(*nx)*(*ny),SEEK_CUR); // skip over zelev and state

  for(;;){

    FORTWUIREAD(&time_local,1);
    if(returncode==0)break;

    FORTWUIREAD(&nchanges,1);
    if(returncode==0)break;

    if(nchanges>0)FSEEK(WUIFILE,16+5*nchanges,SEEK_CUR);

    nt++;

  }
  *times_local=nt;

  fclose(WUIFILE);

  return 0;
}

/* ------------------ getterrain_data ------------------------ */

int getterrain_data(char *file,terraindata *terri){
  FILE *WUIFILE;
  int one;
  int endianswitch=0;
  size_t returncode;
  float time_local;
  int nchanges;
  int nt;
  int nx, ny;
  int *cellindex_buffer;
  unsigned char *cellstate_buffer;
  float *times_local;
  int ntotal;

  WUIFILE = fopen(file,"rb");
  if(WUIFILE==NULL)return 1;

  FSEEK(WUIFILE,4,SEEK_CUR);fread(&one,4,1,WUIFILE);FSEEK(WUIFILE,4,SEEK_CUR);
  if(one!=1)endianswitch=1;

  FSEEK(WUIFILE,12,SEEK_CUR);    // skip over version
  FSEEK(WUIFILE,8+4*4,SEEK_CUR); // skip over xmin,xmax,ymin,ymax
  FSEEK(WUIFILE,8+2*4,SEEK_CUR); // skip over nx, ny

  nx = terri->nx;
  ny = terri->ny;
  ntotal = nx*ny;
  times_local=terri->times;

  NewMemory((void **)&cellindex_buffer,nx*ny*sizeof(int));
  NewMemory((void **)&cellstate_buffer,nx*ny);


  FORTWUIREAD(terri->zcell,ntotal);
  FSEEK(WUIFILE,4,SEEK_CUR);fread(terri->state,1,ntotal,WUIFILE);FSEEK(WUIFILE,4,SEEK_CUR);
  init_tnode(terri);
  init_tnorm(terri);

  for(nt=0;nt<terri->ntimes;nt++){

    FORTWUIREAD(&time_local,1);
    PRINTF("terrain time=%f\n",time_local);
    if(returncode==0)break;
    *times_local++ = time_local;

    FORTWUIREAD(&nchanges,1);
    if(returncode==0)break;

    if(nchanges>0){
      int i;

      FORTWUIREAD(cellindex_buffer,nchanges);
      if(returncode==0)break;
      FSEEK(WUIFILE,4,SEEK_CUR);returncode=fread(cellstate_buffer,1,nchanges,WUIFILE);FSEEK(WUIFILE,4,SEEK_CUR);
      if(returncode==0)break;
      for(i=0;i<nchanges;i++){
        terraincell *ti;
        int ii;

        ti = terri->tcell + cellindex_buffer[i];
        if(ti->nstates+1>ti->nallocated){
          ti->nallocated=ti->nstates+5;
          ResizeMemory((void **)&ti->state,ti->nallocated);
          ResizeMemory((void **)&ti->time,ti->nallocated*sizeof(float));
        }
        ii = ti->nstates;
        ti->state[ii] = cellstate_buffer[i];
        ti->time[ii] = time_local;
        ti->nstates++;
      }
    }

  }

  fclose(WUIFILE);
  FREEMEMORY(cellindex_buffer);
  FREEMEMORY(cellstate_buffer);
  return 0;
}

/* ------------------ init_terraincell ------------------------ */

void init_terraincell(terraindata *terri){
  int i;
  int nx, ny;

  nx = terri->nx;
  ny = terri->ny;

  NewMemory((void **)&terri->tcell,nx*ny*sizeof(terraincell));
  for(i=0;i<terri->nx*terri->ny;i++){
    terraincell *ti;
    int nalloc;

    nalloc=5;
    ti = terri->tcell + i;
    ti->nallocated=nalloc;
    ti->nstates=0;
    ti->state=NULL;
    ti->interval=0;
    ti->time=NULL;
    NewMemory((void **)&ti->state,nalloc);
    NewMemory((void **)&ti->time,nalloc*sizeof(float));
  }

}

/* ------------------ free_terraincell ------------------------ */

void free_terraincell(terraindata *terri){
  int i;
  if(terri->tcell!=NULL){
    for(i=0;i<terri->nx*terri->ny;i++){
      terraincell *ti;

      ti = terri->tcell+i;

      FREEMEMORY(ti->time);
      FREEMEMORY(ti->state);
    }
    FREEMEMORY(terri->tcell);
  }
}

/* ------------------ init_tnode ------------------------ */

void init_tnode(terraindata *terri){
  float *znode, *zcell;
  int i, j;
  int nxcell;

  znode = terri->znode;
  zcell = terri->zcell;
  nxcell = terri->nx;
  for(j=0;j<=terri->ny;j++){
    int jm1, im1, ii, jj;
    float zz;

    jm1 = j - 1;
    if(jm1<0)jm1=0;
    jj = j;
    if(jj==terri->ny)jj--;

    for(i=0;i<=terri->nx;i++){
      im1 = i - 1;
      if(im1<0)im1 = 0;
      ii = i;
      if(ii==terri->nx)ii--;

      zz =  zcell[IJCELL2(im1,jm1)];
      zz += zcell[IJCELL2(im1,jj)];
      zz += zcell[IJCELL2(ii,jm1)];
      zz += zcell[IJCELL2(ii,jj)];
      zz *= 0.25;
      *znode++=zz;
    }
  }
}

/* ------------------ init_tnorm ------------------------ */

void init_tnorm(terraindata *terri){
  float *znode;
  unsigned char *uc_znormal;
  float znormal3[3];
  int i, j;
  int nycell;
  float dx, dy;

  //znormal = terri->znormal;
  uc_znormal = terri->uc_znormal;
  znode = terri->znode;
  nycell = terri->ny;
  dx = (terri->xmax-terri->xmin)/terri->nx;
  dy = (terri->ymax-terri->ymin)/terri->ny;

  for(j=0;j<=terri->ny;j++){
    int jp1, ip1;
    float dzdx, dzdy;
    float sum;

    jp1 = j + 1;
    if(jp1>terri->ny)jp1=terri->ny;

    for(i=0;i<=terri->nx;i++){
      ip1 = i + 1;
      if(ip1>terri->nx)ip1=terri->nx;
      dzdx = (znode[ijnode3(ip1,j)] - znode[ijnode3(i,j)])/dx;
      dzdy = (znode[ijnode3(i,jp1)] - znode[ijnode3(i,j)])/dy;

 //     i  j  k
 //     1  0 dzdx           uu
 //     0  1 dzdy           vv

 //     -dzdx -dzdy 1       uu x vv


//      znormal = terri->znormal + 3*ijnode2(i,j);
      uc_znormal = terri->uc_znormal + ijnode3(i,j);
      znormal3[0] = -dzdx;
      znormal3[1] = -dzdy;
      znormal3[2] = 1.0;

      sum  = znormal3[0]*znormal3[0];
      sum += znormal3[1]*znormal3[1];
      sum += znormal3[2]*znormal3[2];
      sum = sqrt(sum);
      znormal3[0]/=sum;
      znormal3[1]/=sum;
      znormal3[2]/=sum;
      *uc_znormal = getnormalindex(wui_sphereinfo, znormal3);
    }
  }
}

/* ------------------ update_terrain ------------------------ */

void update_terrain(int allocate_memory, float vertical_factor_local){
  if(autoterrain==1){
    int i;

    nterraininfo = nmeshes;
    if(allocate_memory==1){
      NewMemory((void **)&terraininfo,nterraininfo*sizeof(terraindata));
    }

    for(i=0;i<nmeshes;i++){
      meshdata *meshi;
      terraindata *terri;
      float xmin, xmax, ymin, ymax;
      int nx, ny;

      meshi=meshinfo + i;
      terri = terraininfo + i;

      nx = meshi->ibar;
      ny = meshi->jbar;
      xmin = meshi->xplt_orig[0];
      xmax = meshi->xplt_orig[nx];
      ymin = meshi->yplt_orig[0];
      ymax = meshi->yplt_orig[ny];

      initterrain_znode(meshi, terri, xmin, xmax, nx, ymin, ymax, ny, allocate_memory);
      initcontour(&meshi->terrain_contour,rgbptr,nrgb);
    }
    initterrain_all();
  }
  if(nterraininfo>0){
    int imesh;

    for(imesh=0;imesh<nmeshes;imesh++){
      meshdata *meshi;
      terraindata *terri;
      float *znode, *znode_scaled;
      int i, j;

      meshi=meshinfo + imesh;
      terri = meshi->terrain;
      if(terri==NULL)continue;
      terri->terrain_mesh = meshi;
      znode = terri->znode;
      znode_scaled = terri->znode_scaled;

      for(j=0;j<=terri->ny;j++){
        for(i=0;i<=terri->nx;i++){
          *znode_scaled = NORMALIZE_Z(*znode);
          znode++;
          znode_scaled++;
        }
      }
    }
  }
}

/* ------------------ have_terrain_slice ------------------------ */

int have_terrain_slice(void){
  int i;

  for(i=0;i<nsliceinfo;i++){
    slicedata *slicei;

    slicei = sliceinfo + i;

    if(slicei->loaded==1&&slicei->slicetype==SLICE_TERRAIN)return 1;

  }
  return 0;
}

/* ------------------ if_slice_terrain ------------------------ */

void update_terrain_options(void){
  if(nterraininfo>0||autoterrain==1){
    visOpenVents=0;
    visDummyVents=0;
    visFrame=0;
    updatemenu=1;
  }
}

/* ------------------ getmesh_zcell ------------------------ */

float getmesh_zcell(meshdata *meshi, float xval, float yval, int *valid){
  float *xplt, *yplt, *zcell;
  float dx, dy;
  int ibar, jbar;
  int ival, jval;
  float zval;
  int nxcell;

  xplt = meshi->xplt_orig;
  yplt = meshi->yplt_orig;
  ibar = meshi->ibar;
  jbar = meshi->jbar;
  nxcell = ibar;
  *valid = 0;
  if(xval<xplt[0]||xval>xplt[ibar])return 0.0;
  if(yval<yplt[0]||yval>yplt[jbar])return 0.0;

  dx = xplt[1]-xplt[0];
  dy = yplt[1]-yplt[0];
  ival = (xval-xplt[0])/dx;
  if(ival>=ibar)ival = ibar-1;
  jval = (yval-yplt[0])/dy;
  if(jval>=jbar)jval = jbar-1;
  zcell = meshi->zcell;
  zval = zcell[IJCELL2(ival, jval)];
  *valid = 1;
  return zval;
}

/* ------------------ update_mesh_terrain ------------------------ */

void update_mesh_terrain(void){
  int i;

  if(nterraininfo<=0)return;
  for(i=0;i<nmeshes;i++){
    meshdata *meshi;
    meshdata *meshj;
    int ii, jj;
    float *x, *y, *z;
    float xyz[3];

    meshi = meshinfo + i;

    meshi->is_bottom=1;

    x = meshi->xplt_orig;
    y = meshi->yplt_orig;
    z = meshi->zplt_orig;

    xyz[2] = z[0] - (z[1]-z[0])/2.0;

    for(ii=0;ii<meshi->ibar;ii++){
      xyz[0]=x[ii];
      for(jj=0;jj<meshi->jbar;jj++){
        xyz[1]=y[jj];
        meshj = getmesh(xyz);
        if(meshj==NULL||meshj==meshi)continue;
        meshi->is_bottom=0;
        break;
      }
      if(meshi->is_bottom==0)break;
    }
  }

  // compute z level above bottom mesh

  for(i=0;i<nmeshes;i++){
    meshdata *meshi;
    int ii, jj;
    float xyz[3], *x, *y;
    float *zcell;
    int nxcell;

    meshi = meshinfo + i;
    if(meshi->is_bottom==0)continue;
    x = meshi->xplt_orig;
    y = meshi->yplt_orig;
    nxcell = meshi->ibar;
    zcell = meshi->zcell;
    for(ii=0;ii<meshi->ibar;ii++){
      xyz[0]=x[ii];
      for(jj=0;jj<meshi->jbar;jj++){
        int j;

        xyz[1]=y[jj];
        for(j=0;j<nmeshes;j++){
          meshdata *meshj,*mesh_above;

          meshj = meshinfo + j;
          if(meshi==meshj)continue;
          xyz[2]=meshj->zplt_orig[1];
          mesh_above=getmesh(xyz);
          if(mesh_above!=NULL){
            float zz;
            int valid;
            int ij;

            ij = IJCELL2(ii,jj);
            zz=getmesh_zcell(mesh_above, xyz[0],xyz[1], &valid);
            ij=IJCELL2(ii,jj);
            if(valid==1&&zz>zcell[ij]){
              zcell[ij]=zz;
            }
          }
        }
      }
    }
  }
}

/* ------------------ get_z_terrain ------------------------ */

float get_z_terrain(float x, float y){
  int loc;
  float zterrain;

  zterrain = get_zcell_val(NULL,x,y,NULL,&loc);
  return zterrain;
}

