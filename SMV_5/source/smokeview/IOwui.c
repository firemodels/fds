// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>  
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "flowfiles.h"
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "MALLOC.h"
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

// svn revision character string
char IOwui_revision[]="$Revision$";

void init_tnode(terraindata *terri);
void init_tnorm(terraindata *terri);
void init_terraincell(terraindata *terri);
void free_terraincell(terraindata *terri);
void endian_switch(void *val, int nval);

#define ijcell2(i,j) nxcell*(j) + (i)
#define ijnode2(i,j) ((nxcell+1)*(j) + (i))
#define IJKCELL(i,j,k) ((i)+ (j)*ibar+(k)*ibar*jbar)
#define FORTWUIREAD(var,size) fseek(WUIFILE,4,SEEK_CUR);\
                           returncode=fread(var,4,size,WUIFILE);\
                           if(endianswitch==1)endian_switch(var,size);\
                           fseek(WUIFILE,4,SEEK_CUR)

float *get_terraincolor(terraincell *ti);
int getterrain_data(char *file,terraindata *terri);
int getterrain_size(char *file,float *xmin, float *xmax, int *nx, float *ymin, float *ymax, int *ny, int *ntimes);
void drawcone(float d1, float height, float *rgbcolor);
void drawtrunccone(float d1, float d2, float height, float *rgbcolor);
void drawdisk(float diameter, float height, float *rgbcolor);
static float specular[4]={0.4,0.4,0.4,1.0};

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
  glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);
  glTranslatef(-xbar0,-ybar0,-zbar0);
  for(i=0;i<ntreeinfo;i++){
    treedata *treei;
    float crown_height;
    int state;

    treei = treeinfo + i;

    state=0;
    if(showtime==1&&times!=NULL){
      ASSERT(itime>=0)
      if(treei->time_char>0.0&&times[itime]>treei->time_char)state=1;
      if(treei->time_complete>0.0&&times[itime]>treei->time_complete)state=2;
    }

    glPushMatrix();
    glTranslatef(treei->xyz[0],treei->xyz[1],treei->xyz[2]);
    
    switch (state){
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
    }
    glPopMatrix();

 
  }
  glPopMatrix();


  glDisable(GL_COLOR_MATERIAL);

}

/* ------------------ get_zcell ------------------------ */

float get_zcell_val(mesh *meshi,float xval, float yval, int *loc){
  terraindata *terri;
  mesh *meshj;
  int ival, jval;
  float *xplt, *yplt;
  int ibar, jbar;
  float dx, dy;
  float *znode;
  int nxcell;
  float *zcell,zval;
  int imesh;

  for(imesh=-1;imesh<nmeshes;imesh++){
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
      dx = xplt[1]-xplt[0];
      dy = yplt[1]-yplt[0];
      ival = (xval-xplt[0])/dx;
      if(ival>=ibar)ival=ibar-1;
      jval = (yval-yplt[0])/dy;
      if(jval>=jbar)jval=jbar-1;
      terri=meshj->terrain;
      nxcell = terri->nx;
      zcell = terri->zcell;
      zval = zcell[ijcell2(ival,jval)];
      *loc=1;
      zval = zterrain_min+vertical_factor*(zval-zterrain_min);
      return zval;
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
    rgbterrain[4*i+0]=(1.0-f1)*terrain_rgba_zmin[0] + f1*terrain_rgba_zmax[0];
    rgbterrain[4*i+1]=(1.0-f1)*terrain_rgba_zmin[1] + f1*terrain_rgba_zmax[1];
    rgbterrain[4*i+2]=(1.0-f1)*terrain_rgba_zmin[2] + f1*terrain_rgba_zmax[2];
    rgbterrain[4*i+3]=1.0;
  }
}

/* ------------------ initterrain_all ------------------------ */

void initterrain_all(void){
  char buffer[1024];
  float dx, dy;
  float *x, *y;
  float *znode;
  int nxcell;
  //float *znormal;
  float znormal3[3];
  int i,j,k;
  int nz;
  int ibar, jbar;
  mesh *meshi;
  int imesh;
  terraindata *terri;
  float *xplt, *yplt;
  float denom;
  unsigned char *uc_znormal;

  for(imesh=0;imesh<nmeshes;imesh++){
    
    meshi = meshinfo + imesh;

    terri = meshi->terrain;

    dx = terri->x[1]-terri->x[0];
    dy = terri->y[1]-terri->y[0];

    znode = terri->znode;
    nxcell = terri->nx;

    uc_znormal = terri->uc_znormal;
    znode = terri->znode;
    for(j=0;j<=terri->ny;j++){
      int jm1, im1, ii, jj;
      float zz;
      float ynode;

      ynode = terri->y[j];

      for(i=0;i<=terri->nx;i++){
        int ii, jj;
        float xnode;
        int count, loc1, loc2, loc3, loc4;
        float val1, val2, val3, val4;
        float valx1, valx2, valx3, valx4;
        float valx1a, valx2a, valx3a, valx4a;
        float valx1b, valx2b, valx3b, valx4b;
        float valy1a, valy2a, valy3a, valy4a;
        float valy1b, valy2b, valy3b, valy4b;
        float zval;
        float zvalxa, zvalxb;
        float zvalya, zvalyb;
        float dxa, dxb, dya, dyb;
        float dzdx, dzdy;
        float sum;

        xnode = terri->x[i];

        val1 =  get_zcell_val(meshi,xnode-dx/2.0,ynode-dy/2.0,&loc1);
        val2 =  get_zcell_val(meshi,xnode+dx/2.0,ynode-dy/2.0,&loc2);
        val3 =  get_zcell_val(meshi,xnode+dx/2.0,ynode+dy/2.0,&loc3);
        val4 =  get_zcell_val(meshi,xnode-dx/2.0,ynode+dy/2.0,&loc4);
        count = loc1 + loc2 + loc3 + loc4;
        zval = val1*loc1 + val2*loc2 + val3*loc3 + val4*loc4;
        if(count==0)count=1;
        zval /= (float)count;

        *znode++=zval;

 // compute (f(x+dx,y) - f(x-dx,y))/(2*dx)

        valx1a =  get_zcell_val(meshi,xnode-dx-dx/2.0,ynode-dy/2.0,&loc1);
        valx2a =  get_zcell_val(meshi,xnode-dx+dx/2.0,ynode-dy/2.0,&loc2);
        valx3a =  get_zcell_val(meshi,xnode-dx+dx/2.0,ynode+dy/2.0,&loc3);
        valx4a =  get_zcell_val(meshi,xnode-dx-dx/2.0,ynode+dy/2.0,&loc4);
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
        valx1b =  get_zcell_val(meshi,xnode+dx-dx/2.0,ynode-dy/2.0,&loc1);
        valx2b =  get_zcell_val(meshi,xnode+dx+dx/2.0,ynode-dy/2.0,&loc2);
        valx3b =  get_zcell_val(meshi,xnode+dx+dx/2.0,ynode+dy/2.0,&loc3);
        valx4b =  get_zcell_val(meshi,xnode+dx-dx/2.0,ynode+dy/2.0,&loc4);
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

        valy1a =  get_zcell_val(meshi,xnode-dx/2.0,ynode-dy-dy/2.0,&loc1);
        valy2a =  get_zcell_val(meshi,xnode+dx/2.0,ynode-dy-dy/2.0,&loc2);
        valy3a =  get_zcell_val(meshi,xnode+dx/2.0,ynode-dy+dy/2.0,&loc3);
        valy4a =  get_zcell_val(meshi,xnode-dx/2.0,ynode-dy+dy/2.0,&loc4);
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
        valy1b =  get_zcell_val(meshi,xnode-dx/2.0,ynode+dy-dy/2.0,&loc1);
        valy2b =  get_zcell_val(meshi,xnode+dx/2.0,ynode+dy-dy/2.0,&loc2);
        valy3b =  get_zcell_val(meshi,xnode+dx/2.0,ynode+dy+dy/2.0,&loc3);
        valy4b =  get_zcell_val(meshi,xnode-dx/2.0,ynode+dy+dy/2.0,&loc4);
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
        uc_znormal = terri->uc_znormal + ijnode2(i,j);
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
}

/* ------------------ initterrain_znode ------------------------ */

void initterrain_znode(mesh *meshi, terraindata *terri, float xmin, float xmax, int nx, float ymin, float ymax, int ny, 
                       int allocate_memory){
  char buffer[1024];
  float dx, dy;
  float *x, *y, *z;
  float *znode, *zcell;
  int nxcell;
  float *znormal;
  int i,j,k;
  int ijkcell;
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

      ij = ijcell2(i,j);
      zval=meshi->zcell[ij];
      if(zval<zterrain_min)zterrain_min=zval;
      if(zval>zterrain_max)zterrain_max=zval;
      z[ij]=zval;
    }
  }
}

/* ------------------ drawterrain ------------------------ */

void drawterrain(terraindata *terri, int only_geom){
  float *znode, *zn;
  unsigned char *uc_znormal;
  int nxcell;
  int i, j;
  float *x, *y;
  terraincell *ti;
  float terrain_color[4];
  float terrain_shininess=100.0;
  float terrain_specular[4]={0.8,0.8,0.8,1.0};
  float zt_min, zt_max;

  zt_min = zterrain_min;
  zt_max = zterrain_min + vertical_factor*(zterrain_max-zterrain_min);

  terrain_color[0]=0.47843;
  terrain_color[1]=0.45882;
  terrain_color[2]=0.18824; 
  terrain_color[3]=1.0;

  glPushMatrix();
  glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);
  glTranslatef(-xbar0,-ybar0,-zbar0);

  glColor4f(1.0,0.0,0.0,1.0);
  glEnable(GL_LIGHTING);
  glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&terrain_shininess);
  glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,rgbterrain);
  glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,terrain_specular);
  glEnable(GL_COLOR_MATERIAL);

  glBegin(GL_QUADS);
//  znormal = terri->znormal;
  uc_znormal = terri->uc_znormal;
  znode = terri->znode;
  nxcell = terri->nx;
  x = terri->x;
  y = terri->y;
  ti = terri->tcell;
//  glColor4fv(block_ambient2);
  glColor4fv(terrain_color);
  for(j=0;j<terri->ny;j++){
    int jp1;

    jp1 = j + 1;

    for(i=0;i<terri->nx;i++){
      //float *zn;
      unsigned char *uc_zn;
      int ip1;
      float *ter_rgbptr;
      float zval;
      unsigned char izval;

      ip1 = i + 1;

      if(only_geom==0){
        ter_rgbptr = get_terraincolor(ti);
        glColor4fv(ter_rgbptr);
      }
      //zn = znormal+3*ijnode2(i,j);
      uc_zn = uc_znormal+ijnode2(i,j);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));

      glNormal3fv(zn);
      zval = znode[ijnode2(i,j)];
      izval = (MAXRGB-1)*(zval-zt_min)/(zt_max-zt_min);
      glColor4fv(rgbterrain+4*izval);
      glVertex3f(x[i],y[j],zval);

//      zn = znormal+3*ijnode2(ip1,j);
      uc_zn = uc_znormal+ijnode2(ip1,j);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));
      glNormal3fv(zn);
      zval = znode[ijnode2(ip1,j)];
      izval = (MAXRGB-1)*(zval-zt_min)/(zt_max-zt_min);
      glColor4fv(rgbterrain+4*izval);
      glVertex3f(x[i+1],y[j],zval);

//      zn = znormal+3*ijnode2(ip1,jp1);
      uc_zn = uc_znormal+ijnode2(ip1,jp1);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));
      glNormal3fv(zn);
      zval = znode[ijnode2(ip1,jp1)];
      izval = (MAXRGB-1)*(zval-zt_min)/(zt_max-zt_min);
      glColor4fv(rgbterrain+4*izval);
      glVertex3f(x[i+1],y[j+1],zval);

      //zn = znormal+3*ijnode2(i,jp1);
      uc_zn = uc_znormal+ijnode2(i,jp1);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));
      glNormal3fv(zn);
      zval = znode[ijnode2(i,jp1)];
      izval = (MAXRGB-1)*(zval-zt_min)/(zt_max-zt_min);
      glColor4fv(rgbterrain+4*izval);
      glVertex3f(x[i],y[j+1],zval);

      ti++;
    }
  }
  glEnd();
    
  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_LIGHTING);

  glPopMatrix();

}

/* ------------------ drawterrain ------------------------ */

void drawterrain_texture(terraindata *terri, int only_geom){
  float *znode, *znormal, *zn;
  unsigned char *uc_znormal, *uc_zn;
  int nxcell;
  int i, j;
  float *x, *y;
  terraincell *ti;
  float terrain_color[4];

  terrain_color[0]=1.0;
  terrain_color[1]=1.0;
  terrain_color[2]=1.0; 
  terrain_color[3]=1.0;

  glPushMatrix();
  glScalef(mscale[0]/xyzmaxdiff,mscale[1]/xyzmaxdiff,mscale[2]/xyzmaxdiff);
  glTranslatef(-xbar0,-ybar0,-zbar0);

  glEnable(GL_LIGHTING);
  glEnable(GL_NORMALIZE);
  glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D,terrain_texture->name);

  glEnable(GL_COLOR_MATERIAL);
  glColor4fv(terrain_color);
  glBegin(GL_QUADS);
  //znormal = terri->znormal;
  uc_znormal = terri->uc_znormal;
  znode = terri->znode;
  nxcell = terri->nx;
  x = terri->x;
  y = terri->y;
  ti = terri->tcell;
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

//      zn = znormal+3*ijnode2(i,j);
      uc_zn = uc_znormal+ijnode2(i,j);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));
      glNormal3fv(zn);
      glTexCoord2f(tx,ty);
      glVertex3f(x[i],y[j],znode[ijnode2(i,j)]);

//      zn = znormal+3*ijnode2(ip1,j);
      uc_zn = uc_znormal+ijnode2(ip1,j);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));
      glNormal3fv(zn);
      glTexCoord2f(txp1,ty);
      glVertex3f(x[i+1],y[j],znode[ijnode2(ip1,j)]);

//      zn = znormal+3*ijnode2(ip1,jp1);
      uc_zn = uc_znormal+ijnode2(ip1,jp1);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));
      glNormal3fv(zn);
      glTexCoord2f(txp1,typ1);
      glVertex3f(x[i+1],y[j+1],znode[ijnode2(ip1,jp1)]);

//      zn = znormal+3*ijnode2(i,jp1);
      uc_zn = uc_znormal+ijnode2(i,jp1);
      zn = getnormalvectorptr(wui_sphereinfo, (unsigned int)(*uc_zn));
      glNormal3fv(zn);
      glTexCoord2f(tx,typ1);
      glVertex3f(x[i],y[j+1],znode[ijnode2(i,jp1)]);

      ti++;
    }
  }
  glEnd();

  glDisable(GL_TEXTURE_2D);
   
  glDisable(GL_COLOR_MATERIAL);
  glDisable(GL_NORMALIZE);
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

  if(times==NULL||ti->time==NULL){
    index = ti->state[0]%10;
    return rgb_terrain[index];
  }
  
  sv_time = times[itime];
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
    int ileft;

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

  if(ifile>=0&&ifile<nterraininfo)terri = terraininfo + ifile;

  if(flag==UNLOAD){
    if(terri!=NULL){
      FREEMEMORY(terri->x);
      FREEMEMORY(terri->y);
      FREEMEMORY(terri->zcell);
      FREEMEMORY(terri->znode);
      //FREEMEMORY(terri->znormal);
      FREEMEMORY(terri->uc_znormal);
      FREEMEMORY(terri->times);
      free_terraincell(terri);
      terri->loaded=0;
      terri->display=0;
      updatetimes();
    }
    return;
  }
  if(terri==NULL)return;

  if(getterrain_size(file,&xmin, &xmax, &nx, &ymin, &ymax, &ny, &ntimes)!=0)return;

  terri->xmin = xmin;
  terri->xmax = xmax;
  terri->nx = nx;
  terri->ymin = ymin;
  terri->ymax = ymax;
  terri->ny = ny;
  terri->ntimes=ntimes;

  NewMemory((void **)&terri->times,ntimes*sizeof(float));
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
  visTerrain=1;
  plotstate=getplotstate(DYNAMIC_PLOTS);
  updatetimes();
#ifdef _DEBUG
  printf("After terrain file load: ");
  PrintMemoryInfo;
#endif
  IDLE();
  GLUTPOSTREDISPLAY
}

/* ------------------ getterrain_size ------------------------ */

int getterrain_size(char *file,float *xmin, float *xmax, int *nx, float *ymin, float *ymax, int *ny, int *ntimes){
  FILE *WUIFILE;
  int one;
  float xyminmax[4];
  int nxy[2];
  int endianswitch=0;
  size_t returncode;
  int version;
  float time;
  int nchanges;
  int nt=0;

  WUIFILE = fopen(file,"rb");
  if(WUIFILE==NULL)return 1;

  fseek(WUIFILE,4,SEEK_CUR);fread(&one,4,1,WUIFILE);fseek(WUIFILE,4,SEEK_CUR);
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
  
  fseek(WUIFILE,16+5*(*nx)*(*ny),SEEK_CUR); // skip over zelev and state

  for(;;){
    
    FORTWUIREAD(&time,1);
    if(returncode==0)break;

    FORTWUIREAD(&nchanges,1);
    if(returncode==0)break;

    if(nchanges>0)fseek(WUIFILE,16+5*nchanges,SEEK_CUR);

    nt++;

  }
  *ntimes=nt;

  fclose(WUIFILE);

  return 0;


}

/* ------------------ getterrain_data ------------------------ */

int getterrain_data(char *file,terraindata *terri){
  FILE *WUIFILE;
  int one;
  int endianswitch=0;
  size_t returncode;
  float time;
  int nchanges;
  int nt;
  int nx, ny;
  int *cellindex_buffer;
  unsigned char *cellstate_buffer;
  float *times;
  int ntotal;

  WUIFILE = fopen(file,"rb");
  if(WUIFILE==NULL)return 1;

  fseek(WUIFILE,4,SEEK_CUR);fread(&one,4,1,WUIFILE);fseek(WUIFILE,4,SEEK_CUR);
  if(one!=1)endianswitch=1;

  fseek(WUIFILE,12,SEEK_CUR);    // skip over version
  fseek(WUIFILE,8+4*4,SEEK_CUR); // skip over xmin,xmax,ymin,ymax
  fseek(WUIFILE,8+2*4,SEEK_CUR); // skip over nx, ny

  nx = terri->nx;
  ny = terri->ny;
  ntotal = nx*ny;
  times=terri->times;

  NewMemory((void **)&cellindex_buffer,nx*ny*sizeof(int));
  NewMemory((void **)&cellstate_buffer,nx*ny);


  FORTWUIREAD(terri->zcell,ntotal); 
  fseek(WUIFILE,4,SEEK_CUR);fread(terri->state,1,ntotal,WUIFILE);fseek(WUIFILE,4,SEEK_CUR);
  init_tnode(terri);
  init_tnorm(terri);
  
  for(nt=0;nt<terri->ntimes;nt++){
    
    FORTWUIREAD(&time,1);
    printf("terrain time=%f\n",time);
    if(returncode==0)break;
    *times++ = time;

    FORTWUIREAD(&nchanges,1);
    if(returncode==0)break;

    if(nchanges>0){
      int i;

      FORTWUIREAD(cellindex_buffer,nchanges);
      if(returncode==0)break;
      fseek(WUIFILE,4,SEEK_CUR);returncode=fread(cellstate_buffer,1,nchanges,WUIFILE);fseek(WUIFILE,4,SEEK_CUR);
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
        ti->time[ii] = time;
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

      zz =  zcell[ijcell2(im1,jm1)];
      zz += zcell[ijcell2(im1,jj)];
      zz += zcell[ijcell2(ii,jm1)];
      zz += zcell[ijcell2(ii,jj)];
      zz *= 0.25;
      *znode++=zz;
    }
  }
}

/* ------------------ init_tnorm ------------------------ */

void init_tnorm(terraindata *terri){
  float *znormal, *znode;
  unsigned char *uc_znormal;
  float znormal3[3];
  int i, j;
  int nxcell;
  float dx, dy;

  //znormal = terri->znormal;
  uc_znormal = terri->uc_znormal;
  znode = terri->znode;
  nxcell = terri->nx;
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
      dzdx = (znode[ijnode2(ip1,j)] - znode[ijnode2(i,j)])/dx;
      dzdy = (znode[ijnode2(i,jp1)] - znode[ijnode2(i,j)])/dy;

 //     i  j  k
 //     1  0 dzdx           uu
 //     0  1 dzdy           vv

 //     -dzdx -dzdy 1       uu x vv

      
//      znormal = terri->znormal + 3*ijnode2(i,j);
      uc_znormal = terri->uc_znormal + ijnode2(i,j);
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

void update_terrain(int allocate_memory, float vertical_factor){
  int i, j;

  if(autoterrain==1){

    zterrain_min=1000000000.0;
    zterrain_max=-zterrain_min;

    nterraininfo = nmeshes;
    if(allocate_memory==1){
      NewMemory((void **)&terraininfo,nterraininfo*sizeof(terraindata));
    }

    for(i=0;i<nmeshes;i++){
      mesh *meshi;
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
    }
    initterrain_all();
  }
  if(nterraininfo>0){
    int imesh;

    for(imesh=0;imesh<nmeshes;imesh++){
      mesh *meshi;
      terraindata *terri;
      float *znode, *znode_scaled;
      int i, j;

      meshi=meshinfo + imesh;
      terri = meshi->terrain;
      if(terri==NULL)continue;
      znode = terri->znode;
      znode_scaled = terri->znode_scaled;

      for(j=0;j<=terri->ny;j++){
        for(i=0;i<=terri->nx;i++){
          *znode_scaled++ = (*znode++-zbar0)/xyzmaxdiff;
        }
      }

    }
  }
}

