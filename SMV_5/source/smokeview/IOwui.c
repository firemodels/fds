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

#define ijnode2(i,j) ((nxcell+1)*(j) + (i))

/* ------------------ getterrain_z ------------------------ */

float getterrain_z(float x, float y){
  int i;
  float returnval=0.0;
  int ii, jj;
  int iip1, jjp1;
  float dx, dy;
  float z11, z12, z21, z22;
  int nxcell;
  float *xnode, *ynode, *znode;

  for(i=0;i<nterraininfo;i++){
    terraindata *terri;

    terri = terraininfo + i;
    if(x<terri->xmin||x>terri->xmax)continue;
    if(y<terri->ymin||y>terri->ymax)continue;
    nxcell = terri->nx;

    dx = (terri->xmax-terri->xmin)/terri->nx;
    dy = (terri->ymax-terri->ymin)/terri->ny;
    ii = (x-terri->xmin)/dx;
    if(ii<0)ii=0;
    if(ii>terri->nx-1)ii=terri->nx-1;

    jj = (y-terri->ymin)/dy;
    if(jj<0)jj=0;
    if(jj>terri->ny-1)jj=terri->ny-1;
    iip1 = ii+1;
    jjp1 = jj+1;

    znode = terri->znode;
    xnode = terri->x;
    ynode = terri->y;

    z11 = znode[ijnode2(  ii,  jj)];
    z12 = znode[ijnode2(  ii,jjp1)];
    z21 = znode[ijnode2(iip1,  jj)];
    z22 = znode[ijnode2(iip1,jjp1)];

    returnval  =  (xnode[iip1]-x)*(ynode[jjp1]-y)*z11;
    returnval  += (xnode[iip1]-x)*(y-ynode[jj])*z12;
    returnval  += (x-xnode[ii])  *(ynode[jjp1]-y)*z21;
    returnval  += (x-xnode[ii])  *(y-ynode[jj])*z22;
    returnval /= (dx*dy);
  }
  return returnval;
}

/* ------------------ initterrain ------------------------ */
#define IJKCELL(i,j,k) ((i)+ (j)*ibar+(k)*ibar*jbar)

void initterrain(FILE *stream, mesh *meshi, terraindata *terri, float xmin, float xmax, int nx, float ymin, float ymax, int ny){
  char buffer[1024];
  float dx, dy;
  float *x, *y, *z;
  float *znode, *zcell;
  int nxcell;
  float *znormal;
  int i,j,k;
  int nz;
  int ibar, jbar;

  terri->x=NULL;
  terri->y=NULL;
  terri->display=0;
  terri->loaded=0;
  terri->autoload=0;
  terri->times=NULL;
  terri->ntimes=0;
  terri->times=NULL;
  terri->zcell=NULL;
  terri->znode=NULL;
  terri->znormal=NULL;
  terri->tcell=NULL;
  terri->ter_texture=NULL;
  terri->state=NULL;
  terri->timeslist=NULL;
  terri->zcell=NULL;
  terri->znode=NULL;
  terri->znormal=NULL;

  terri->xmin=xmin;
  terri->xmax=xmax;
  terri->ymin=ymin;
  terri->ymax=ymax;
  terri->ny=ny;
  if(nx<0){
    nx=-nx;
  }
  terri->nx=nx;

  NewMemory((void **)&terri->x,(nx+1)*sizeof(float));
  NewMemory((void **)&terri->y,(ny+1)*sizeof(float));
  NewMemory((void **)&terri->zcell,nx*ny*sizeof(float));
  NewMemory((void **)&terri->state,nx*ny);
  NewMemory((void **)&terri->znode,(nx+1)*(ny+1)*sizeof(float));
  NewMemory((void **)&terri->znormal,3*(nx+1)*(ny+1)*sizeof(float));

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

#define ijcell2(i,j) nxcell*(j) + (i)
  z=terri->zcell;
  if(stream==NULL){
    int *iblank_cell;
    int ijkcell;
    int ij;

    iblank_cell = meshi->iblank_cell;
    ibar = nx;
    jbar = ny;
    nxcell = nx;
    nz = meshi->kbar;
    for(j=0;j<ny;j++){
      for(i=0;i<nx;i++){
        ij = ijcell2(i,j);
        z[ij]=meshi->zplt_orig[0];
        for(k=nz-1;k>=0;k--){
          ijkcell = IJKCELL(i,j,k);
          if(iblank_cell==NULL||iblank_cell[ijkcell]==0){
            z[ij]=meshi->zplt_orig[k];
            break;
          }
        }
      }
    }
  }
  else{
    for(i=0;i<nx*ny;i++){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f",z);
      z++;
    }
  }
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
#define ijnode2(i,j) ((nxcell+1)*(j) + (i))
  znormal = terri->znormal;  ;
  znode = terri->znode;
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
     //     1  0 dzdx
     //     0  1 dzdy

     //     -dzdx -dzdy 1

          
      znormal = terri->znormal + 3*ijnode2(i,j);
      znormal[0] = -dzdx;
      znormal[1] = -dzdy;
      znormal[2] = 1.0;

      sum  = znormal[0]*znormal[0];
      sum += znormal[1]*znormal[1];
      sum += znormal[2]*znormal[2];
      sum = sqrt(sum);
      znormal[0]/=sum;
      znormal[1]/=sum;
      znormal[2]/=sum;
    }
  }
}

/* ------------------ drawterrain ------------------------ */

#define ijnode2(i,j) ((nxcell+1)*(j) + (i))
void drawterrain(terraindata *terri, int only_geom){
  float *znode, *znormal;
  int nxcell;
  int i, j;
  float *x, *y;
  terraincell *ti;
  float terrain_color[4];

  terrain_color[0]=0.47843;
  terrain_color[1]=0.45882;
  terrain_color[2]=0.18824; 
  terrain_color[3]=1.0;

  glPushMatrix();
  glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);
  glTranslatef(-xbar0,-ybar0,-zbar0);

  glColor4f(1.0,0.0,0.0,1.0);
  glEnable(GL_LIGHTING);
  glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
//  glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
  glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,terrain_color);
  glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,specular);
  glEnable(GL_COLOR_MATERIAL);

  glBegin(GL_QUADS);
  znormal = terri->znormal;  ;
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
      float *zn;
      int ip1;
      float *ter_rgbptr;

      ip1 = i + 1;

      if(only_geom==0){
        ter_rgbptr = get_terraincolor(ti);
        glColor4fv(ter_rgbptr);
      }
      zn = znormal+3*ijnode2(i,j);
      glNormal3fv(zn);
      glVertex3f(x[i],y[j],znode[ijnode2(i,j)]);

      zn = znormal+3*ijnode2(ip1,j);
      glNormal3fv(zn);
      glVertex3f(x[i+1],y[j],znode[ijnode2(ip1,j)]);

      zn = znormal+3*ijnode2(ip1,jp1);
      glNormal3fv(zn);
      glVertex3f(x[i+1],y[j+1],znode[ijnode2(ip1,jp1)]);

      zn = znormal+3*ijnode2(i,jp1);
      glNormal3fv(zn);
      glVertex3f(x[i],y[j+1],znode[ijnode2(i,jp1)]);

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
  float *znode, *znormal;
  int nxcell;
  int i, j;
  float *x, *y;
  terraincell *ti;

  glPushMatrix();
  glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);
  glTranslatef(-xbar0,-ybar0,-zbar0);

  glEnable(GL_LIGHTING);
  glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
  glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_MODULATE);
  glEnable(GL_TEXTURE_2D);
  glBindTexture(GL_TEXTURE_2D,terrain_texture->name);

  glBegin(GL_QUADS);
  znormal = terri->znormal;  ;
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
      float *ter_rgbptr;
      float tx,txp1;

      ip1 = i + 1;
      tx = (x[i]-xbar0ORIG)/(xbarORIG-xbar0ORIG);
      txp1 = (x[i+1]-xbar0ORIG)/(xbarORIG-xbar0ORIG);

      zn = znormal+3*ijnode2(i,j);
      glNormal3fv(zn);
      glTexCoord2f(tx,ty);
      glVertex3f(x[i],y[j],znode[ijnode2(i,j)]);

      zn = znormal+3*ijnode2(ip1,j);
      glNormal3fv(zn);
      glTexCoord2f(txp1,ty);
      glVertex3f(x[i+1],y[j],znode[ijnode2(ip1,j)]);

      zn = znormal+3*ijnode2(ip1,jp1);
      glNormal3fv(zn);
      glTexCoord2f(txp1,typ1);
      glVertex3f(x[i+1],y[j+1],znode[ijnode2(ip1,jp1)]);

      zn = znormal+3*ijnode2(i,jp1);
      glNormal3fv(zn);
      glTexCoord2f(tx,typ1);
      glVertex3f(x[i],y[j+1],znode[ijnode2(i,jp1)]);

      ti++;
    }
  }
  glEnd();

  glDisable(GL_TEXTURE_2D);
   
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

/*
typedef struct {
  int nallocated, nstates;
  float *time, *tcurrent;
  unsigned char *state;
} terraincell;

typedef struct {
  char *file;
  int loaded, display;
  int autoload;
  texture *ter_texture;
  int nx, ny;
  float xmin, xmax, ymin, ymax;
  float *x, *y;
  float *zcell, *znode, *znormal;
  float *times;
  terraincell *tcell;
  int ntimes;
} terraindata;
*/

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
      FREEMEMORY(terri->znormal);
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
  NewMemory((void **)&terri->znormal,3*(nx+1)*(ny+1)*sizeof(float));
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
  glutPostRedisplay();
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
  int i, j;
  int nxcell;
  float dx, dy;

  znormal = terri->znormal;
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

      
      znormal = terri->znormal + 3*ijnode2(i,j);
      znormal[0] = -dzdx;
      znormal[1] = -dzdy;
      znormal[2] = 1.0;

      sum  = znormal[0]*znormal[0];
      sum += znormal[1]*znormal[1];
      sum += znormal[2]*znormal[2];
      sum = sqrt(sum);
      znormal[0]/=sum;
      znormal[1]/=sum;
      znormal[2]/=sum;
    }
  }
}
