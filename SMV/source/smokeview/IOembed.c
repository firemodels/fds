// $Date$ 
// $Revision$
// $Author$

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

#include "smokeviewvars.h"
#include "smv_endian.h"

// svn revision character string
char IOembed_revision[]="$Revision$";

/* ------------------ CalcTriNormal ------------------------ */

void CalcTriNormal(float *v1, float *v2, float *v3, float *norm){
  float u[3], v[3];
  int i;

  for(i=0;i<3;i++){
    u[i]=v2[i]-v1[i];
    v[i]=v3[i]-v1[i];
  }
  /*
     i   j  k
     ux uy uz
     vx vy vz
  */
  norm[0]=u[1]*v[2]-u[2]*v[1];
  norm[1]=u[2]*v[0]-u[0]*v[2];
  norm[2]=u[0]*v[1]-u[1]*v[0];
  ReduceToUnit(norm);
}

/* ----------------------- compare_verts ----------------------------- */

int compare_verts( const void *arg1, const void *arg2 ){
  point *pointi, *pointj;
  float *xyzi, *xyzj;

  pointi = (point *)arg1;
  pointj = (point *)arg2;
  xyzi = pointi->xyz;
  xyzj = pointj->xyz;

  if(xyzi[0]<xyzj[0])return -1;
  if(xyzi[0]>xyzj[0])return 1;
  if(xyzi[1]<xyzj[1])return -1;
  if(xyzi[1]>xyzj[1])return 1;
  if(xyzi[2]<xyzj[2])return -1;
  if(xyzi[2]>xyzj[2])return 1;
  return 0;
}

/* ------------------ distxy ------------------------ */

float distxy(float *x, float *y){
  float r1, r2, r3;

  r1 = x[0]-y[0];
  r2 = x[1]-y[1];
  r3 = x[2]-y[2];
  return sqrt(r1*r1+r2*r2+r3*r3);
}

/* ------------------ get_angle ------------------------ */

float get_angle(float d1, float d2, float d3){
  float angle;
  float arg;

  arg = (d2*d2+d3*d3-d1*d1)/(2.0*d2*d3);
  if(arg<-1.0)arg=-1.0;
  if(arg>1.0)arg=1.0;
  angle = acos(arg)*180.0/(4.0*atan(1.0));
  return angle;
}

/* ------------------ get_minangle ------------------------ */

float get_minangle(triangle *trii){
  float minangle;
  float d1, d2, d3;
  float *xyz1, *xyz2, *xyz3;
  float angle1, angle2, angle3;

  xyz1 = trii->points[0]->xyz;
  xyz2 = trii->points[1]->xyz;
  xyz3 = trii->points[2]->xyz;
  d1 = distxy(xyz1,xyz2);
  d2 = distxy(xyz1,xyz3);
  d3 = distxy(xyz2,xyz3);
  angle1 = get_angle(d1,d2,d3);
  angle2 = get_angle(d2,d1,d3);
  angle3 = get_angle(d3,d1,d2);
  minangle = angle1;
  if(angle2<minangle)minangle=angle2;
  if(angle3<minangle)minangle=angle3;
  return minangle;
}

/* ------------------ draw_faceinfo ------------------------ */

void get_faceinfo(void){
  int i;

  for(i=0;i<ngeominfo;i++){
    geomdata *geomi;
    geomlistdata *geomlisti;
    point **points;
    int j;
    int ndups=0,nused=0,nskinny=0;

    geomi = geominfo + i;
    geomlisti = geomi->geomlistinfo;

    if(geomlisti->npoints>0){
      NewMemory((void **)&points,geomlisti->npoints*sizeof(point *));
      for(j=0;j<geomlisti->npoints;j++){
        points[j]=geomlisti->points+j;
        points[j]->nused=0;
      }
      qsort(points,geomlisti->npoints,sizeof(point *),compare_verts);
      for(j=1;j<geomlisti->npoints;j++){
        if(compare_verts(points[j-1],points[j])==0)ndups++;
      }
      for(j=0;j<geomlisti->ntriangles;j++){
        triangle *trii;
        float min_angle;

        trii = geomlisti->triangles + j;
        trii->points[0]->nused++;
        trii->points[1]->nused++;
        trii->points[2]->nused++;
        if(get_minangle(trii)<=10.0){
          trii->skinny=1;
          nskinny++;
        }
        else{
          trii->skinny=0;
        }

      }
      for(j=0;j<geomlisti->npoints;j++){
        if(points[j]->nused>0)nused++;
      }
      printf("Face/Vertex Summary\n");
      printf("      Faces: %i\n",geomlisti->ntriangles);
      printf(" slim faces: %i\n",nskinny);
      printf("   Vertices: %i\n",geomlisti->npoints);
      printf("     unused: %i\n",geomlisti->npoints-nused);
      printf(" duplicates: %i\n\n",ndups);
      FREEMEMORY(points);
    }
  }
}

/* ------------------ draw_geom ------------------------ */

void draw_geom(int flag){
  int i;
  float black[]={0.0,0.0,0.0,1.0};
  float blue[]={0.0,0.0,1.0,1.0};
  float skinny_color[]={1.0,0.0,0.0,1.0};
  float *last_color=NULL;


  if(patchembedded==0&&showtrisurface==1){
    int ntris;
    triangle **tris;
    float *color;

    if(flag==DRAW_OPAQUE){
      ntris=nopaque_triangles;
      tris=opaque_triangles;
    }
    if(flag==DRAW_TRANSPARENT){
      if(use_transparency_data==1)transparenton();
      ntris=ntransparent_triangles;
      tris=transparent_triangles;
    }
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glEnable(GL_COLOR_MATERIAL);
    
    glPushMatrix();
    glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);
    glTranslatef(-xbar0,-ybar0,-zbar0);
    glBegin(GL_TRIANGLES);
    for(i=0;i<ntris;i++){
      triangle *trianglei;
      float *xyzptr[3];
      float *xyznorm;

      trianglei = tris[i];

      if(hilight_skinny==1&&trianglei->skinny==1){
        color=skinny_color;
      }
      else{
        color = trianglei->surf->color;
      }
      if(color!=last_color){
        glColor4fv(color);
        last_color=color;
      }
      if(smoothtrinormal==0){

        xyznorm=trianglei->normal;
        glNormal3fv(xyznorm);

        xyzptr[0] = trianglei->points[0]->xyz;
        glVertex3fv(xyzptr[0]);

        xyzptr[1] = trianglei->points[1]->xyz;
        glVertex3fv(xyzptr[1]);

        xyzptr[2] = trianglei->points[2]->xyz;
        glVertex3fv(xyzptr[2]);
      }
      else{
        xyznorm = trianglei->points[0]->norm;
        glNormal3fv(xyznorm);
        xyzptr[0] = trianglei->points[0]->xyz;
        glVertex3fv(xyzptr[0]);

        xyznorm = trianglei->points[1]->norm;
        glNormal3fv(xyznorm);
        xyzptr[1] = trianglei->points[1]->xyz;
        glVertex3fv(xyzptr[1]);

        xyznorm = trianglei->points[2]->norm;
        glNormal3fv(xyznorm);
        xyzptr[2] = trianglei->points[2]->xyz;
        glVertex3fv(xyzptr[2]);
      }
    }
    glEnd();
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
    glPopMatrix();
  }
  if(flag==DRAW_TRANSPARENT){
    if(use_transparency_data==1)transparentoff();
    return;
  }

  for(i=0;i<ngeominfo;i++){
    geomdata *geomi;
    geomlistdata *geomlisti;
    int ntris,npoints;
    int j;
    float *color;

    geomi = geominfo + i;
    if(geomi->loaded==0||geomi->display==0)continue;
    geomlisti = geomi->geomlistinfo;
    ntris = geomlisti->ntriangles;
    npoints = geomlisti->npoints;
    if(showtrioutline==1){
      glPushMatrix();
      glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);
      glTranslatef(-xbar0,-ybar0,-zbar0);
      glBegin(GL_LINES);
      for(j=0;j<ntris;j++){
        float *xyzptr[3];
        float *xyznorm;
        triangle *trianglei;

        trianglei = geomlisti->triangles+j;
       
        xyznorm=trianglei->normal;
        glNormal3fv(xyznorm);

        xyzptr[0] = trianglei->points[0]->xyz;
        xyzptr[1] = trianglei->points[1]->xyz;
        xyzptr[2] = trianglei->points[2]->xyz;

        if(showtrisurface==1){
          color = black;
        }
        else{
          color = trianglei->surf->color;
        }
        glColor3fv(color);
        glVertex3fv(xyzptr[0]);
        glVertex3fv(xyzptr[1]);
        glVertex3fv(xyzptr[1]);
        glVertex3fv(xyzptr[2]);
        glVertex3fv(xyzptr[2]);
        glVertex3fv(xyzptr[0]);
      }
      glEnd();
      glPopMatrix();
    }
    if(showtrinormal==1){
      if(smoothtrinormal==0){
        glPushMatrix();
        glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);
        glTranslatef(-xbar0,-ybar0,-zbar0);
        glBegin(GL_LINES);
        for(j=0;j<ntris;j++){
          float *p1, *p2, *p3;
          float *xyznorm;
          triangle *trianglei;
          float xyz1[3], xyz2[3];

          trianglei = geomlisti->triangles+j;
       
          xyznorm=trianglei->normal;

          p1 = trianglei->points[0]->xyz;
          p2 = trianglei->points[1]->xyz;
          p3 = trianglei->points[2]->xyz;

          xyz1[0] = (p1[0] + p2[0] + p3[0])/3.0;
          xyz1[1] = (p1[1] + p2[1] + p3[1])/3.0;
          xyz1[2] = (p1[2] + p2[2] + p3[2])/3.0;
          xyz2[0] = xyz1[0] + 0.1*xyznorm[0];
          xyz2[1] = xyz1[1] + 0.1*xyznorm[1];
          xyz2[2] = xyz1[2] + 0.1*xyznorm[2];

          if(trianglei->fdsnorm==1){
            color = black;
          }
          else{
            color=blue;
          }
          glColor3fv(color);
          glVertex3fv(xyz1);
          glVertex3fv(xyz2);
        }
        glEnd();
        
        glPointSize(6.0);
        glBegin(GL_POINTS);
        for(j=0;j<ntris;j++){
          float *p1, *p2, *p3;
          float *xyznorm;
          triangle *trianglei;
          float xyz1[3], xyz2[3];

          trianglei = geomlisti->triangles+j;
       
          xyznorm=trianglei->normal;

          p1 = trianglei->points[0]->xyz;
          p2 = trianglei->points[1]->xyz;
          p3 = trianglei->points[2]->xyz;

          xyz1[0] = (p1[0] + p2[0] + p3[0])/3.0;
          xyz1[1] = (p1[1] + p2[1] + p3[1])/3.0;
          xyz1[2] = (p1[2] + p2[2] + p3[2])/3.0;
          xyz2[0] = xyz1[0] + 0.1*xyznorm[0];
          xyz2[1] = xyz1[1] + 0.1*xyznorm[1];
          xyz2[2] = xyz1[2] + 0.1*xyznorm[2];

          if(trianglei->fdsnorm==1){
            color = black;
          }  
          else{
            color=blue;
          }
          glColor3fv(color);
          glVertex3fv(xyz2);
        }
        glEnd();
        glPopMatrix();
      }
      if(smoothtrinormal==1){
        glPushMatrix();
        glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);
        glTranslatef(-xbar0,-ybar0,-zbar0);
        glBegin(GL_LINES);
        for(j=0;j<npoints;j++){
          float *xyznorm;
          point *pointi;
          float *xyz1, xyz2[3];

          pointi = geomlisti->points+j;
          xyznorm = pointi->norm;       
          xyz1 = pointi->xyz;

          xyz2[0] = xyz1[0] + 0.1*xyznorm[0];
          xyz2[1] = xyz1[1] + 0.1*xyznorm[1];
          xyz2[2] = xyz1[2] + 0.1*xyznorm[2];

          color = black;
          glColor3fv(color);
          glVertex3fv(xyz1);
          glVertex3fv(xyz2);
        }
        glEnd();

        glPointSize(6.0);
        glBegin(GL_POINTS);
        for(j=0;j<npoints;j++){
          float *xyznorm;
          point *pointi;
          float *xyz1, xyz2[3];

          pointi = geomlisti->points+j;
          xyznorm = pointi->norm;       
          xyz1 = pointi->xyz;

          xyz2[0] = xyz1[0] + 0.1*xyznorm[0];
          xyz2[1] = xyz1[1] + 0.1*xyznorm[1];
          xyz2[2] = xyz1[2] + 0.1*xyznorm[2];

          color = black;
          glColor3fv(color);
          glVertex3fv(xyz2);
        }
        glEnd();
        glPopMatrix();
      }
    }
  }
}

/* ------------------ update_triangles ------------------------ */

void update_triangles(void){
  int j;

  for(j=0;j<ngeominfo;j++){
    geomlistdata *geomlisti;
    geomdata *geomi;
    float *xyzptr[3];
    float *xyznorm;
    float *xyz;
    int i;

    geomi = geominfo + j;
    if(geomi->loaded==0||geomi->display==0)continue;
    geomlisti = geomi->geomlistinfo;
    
    for(i=0;i<geomlisti->npoints;i++){
      xyz = geomlisti->points[i].xyz;
      xyz[0] = (xyz[0] - xbar0)/xyzmaxdiff;
      xyz[1] = (xyz[1] - ybar0)/xyzmaxdiff;
      xyz[2] = (xyz[2] - zbar0)/xyzmaxdiff;
      xyz += 3;
    }
    
    for(i=0;i<geomlisti->ntriangles;i++){
      triangle *trianglei;

      trianglei = geomlisti->triangles+i;
      if(geomlisti->triangles[i].fdsnorm==0){
        xyzptr[0] = trianglei->points[0]->xyz;
        xyzptr[1] = trianglei->points[1]->xyz;
        xyzptr[2] = trianglei->points[2]->xyz;
        xyznorm = trianglei->normal;
        CalcTriNormal(xyzptr[0],xyzptr[1],xyzptr[2],xyznorm);
      }
    }

    for(i=0;i<geomlisti->npoints;i++){
      point *pointi;

      pointi = geomlisti->points + i;
      pointi->ntriangles=0;
      pointi->itriangle=0;
    }
    for(i=0;i<geomlisti->ntriangles;i++){
      triangle *trianglei;

      trianglei = geomlisti->triangles+i;
      trianglei->points[0]->ntriangles++;
      trianglei->points[1]->ntriangles++;
      trianglei->points[2]->ntriangles++;
    }
    for(i=0;i<geomlisti->npoints;i++){
      point *pointi;

      pointi = geomlisti->points + i;
      if(pointi->ntriangles>0){
        NewMemory((void **)&pointi->triangles,pointi->ntriangles*sizeof(triangle *));
      }
    }
    for(i=0;i<geomlisti->ntriangles;i++){
      triangle *trianglei;
      point *pointi;

      trianglei = geomlisti->triangles+i;
      pointi = trianglei->points[0];
      pointi->triangles[pointi->itriangle++]=trianglei;
      pointi = trianglei->points[1];
      pointi->triangles[pointi->itriangle++]=trianglei;
      pointi = trianglei->points[2];
      pointi->triangles[pointi->itriangle++]=trianglei;
    }
    for(i=0;i<geomlisti->npoints;i++){
      point *pointi;
      int k;
      float *norm;

      pointi = geomlisti->points + i;
      norm=pointi->norm;
      norm[0]=0.0;
      norm[1]=0.0;
      norm[2]=0.0;
      for(k=0;k<pointi->ntriangles;k++){
        float *norm2;
        triangle *trianglei;

        trianglei = pointi->triangles[k];
        norm2 = trianglei->normal;
        norm[0]+=norm2[0];
        norm[1]+=norm2[1];
        norm[2]+=norm2[2];
      }
      ReduceToUnit(norm);
    }
  }
}

#define FORTREAD(var,count,STREAM) fseek(STREAM,4,SEEK_CUR);\
                           returncode=fread(var,4,count,STREAM);\
                           if(returncode!=count)returncode=0;\
                           if(endianswitch==1&&returncode!=0)endian_switch(var,count);\
                           fseek(STREAM,4,SEEK_CUR)

#define FORTREADBR(var,count,STREAM) FORTREAD(var,count,STREAM);if(returncode==0)break;

/* ------------------ get_geom_header ------------------------ */

void get_geom_header(char *file, int *ntimes){
  FILE *stream;
  int one=0,endianswitch=0;
  int nvertfaces[4];
  float times[2];
  int first=1;
  int nt;
  int returncode;
  int version;
  int size;

  stream = fopen(file,"rb");
  if(stream==NULL){
    *ntimes=-1;
    return;
  }
  fseek(stream,4,SEEK_CUR);fread(&one,4,1,stream);fseek(stream,4,SEEK_CUR);
  if(one!=1)endianswitch=1;
  FORTREAD(&version,1,stream);
  nt=0;
  for(;;){
    if(first==1){
      FORTREADBR(times,1,stream);
      FORTREADBR(&nvertfaces,4,stream);
      if(nvertfaces[0]!=0)fseek(stream,4+3*nvertfaces[0]*4+4,SEEK_CUR);    
      if(nvertfaces[1]!=0)fseek(stream,4+3*nvertfaces[1]*4+4,SEEK_CUR);    
      if(nvertfaces[2]!=0)fseek(stream,4+3*nvertfaces[2]*4+4,SEEK_CUR);    
      if(nvertfaces[3]!=0)fseek(stream,4+3*nvertfaces[3]*4+4,SEEK_CUR);    
    }
    else{
      int *geom_type;

      FORTREADBR(times,2,stream);
      geom_type = (int *)(times+1);

      if(*geom_type==0){
        FORTREADBR(&nvertfaces,2,stream);
        if(nvertfaces[0]!=0)fseek(stream,4+3*nvertfaces[0]*4+4,SEEK_CUR);    
        if(nvertfaces[1]!=0)fseek(stream,4+3*nvertfaces[1]*4+4,SEEK_CUR);    
      }
      else{
        fseek(stream,4+8*4+4,SEEK_CUR);
      }
    }

    nt++;
  }
  *ntimes=nt;
  fclose(stream);
}

/* ------------------ get_geomdata_header ------------------------ */

void get_geomdata_header(char *file, int *ntimes, int *nvals){
  FILE *stream;
  int one=1,endianswitch=0;
  int nface_static,nface_dynamic;
  float time;
  int nt,nv;
  int returncode;

  stream = fopen(file,"r");
  if(stream==NULL){
    *ntimes=-1;
    return;
  }
  fseek(stream,4,SEEK_CUR);fread(&one,4,1,stream);fseek(stream,4,SEEK_CUR);
  if(one!=1)endianswitch=1;
  nt=-1;
  nv=0;
  for(;;){
    FORTREADBR(&time,1,stream);
    FORTREADBR(&nface_static,1,stream);
    if(nface_static!=0)fseek(stream,4+nface_static*4+4,SEEK_CUR);    
    FORTREADBR(&nface_dynamic,1,stream);
    if(nface_dynamic!=0)fseek(stream,4+nface_dynamic*4+4,SEEK_CUR);    
    nt++;
    nv+=(nface_static+nface_dynamic);
  }
  *ntimes=nt;
  *nvals=nv;
  fclose(stream);
}

/* ------------------ read_all_geom ------------------------ */

void read_all_geom(void){
  int i, errorcode;

  for(i=0;i<ngeominfo;i++){
    read_geom(i,LOAD,&errorcode);
  }
}

/* ------------------ read_geom ------------------------ */

void read_geom(int ifile, int flag, int *errorcode){
  geomdata *geomi;
  char *file;
  FILE *stream;
  int one=1, endianswitch=0;
  int returncode;
  int ntimes;
  float *xyz=NULL;
  int nxyz=0;
  int i;
  point *points;
  triangle *triangles;
  int first=1;
  int version;

  // 1
  // version
  // stime
  // nvert_s, nface_s, nvert_d, nface_d
  // x1 y1 z1 ... xnvert_s ynvert_s znvert_s
  // x1 y1 z1 ... xnvert_d ynvert_d znvert_d
  // i1 j1 k1 ... inface_s jnface_s knface_s
  // i1 j1 k1 ... inface_d jnface_d knface_d

  // time geom_type
  // if geom_type==0
  //   x1 y1 z1 ... xnvert_d ynvert_d znvert_d  
  //   i1 j1 k1 ... inface_d jnface_d knface_d
  // if geom_type==1
  //   xtran, ytran, ztran, xrot0, yrot0, zrot0, rot_az, rot_elev


  geomi = geominfo + ifile;

  if(geomi->ntimes>0){
    for(i=0;i<geomi->ntimes;i++){
      geomlistdata *geomlisti;

      geomlisti = geomi->geomlistinfo+i;
      FREEMEMORY(geomlisti->points);
      FREEMEMORY(geomlisti->triangles);
    }
    FREEMEMORY(geomi->times);
    FREEMEMORY(geomi->geomlistinfo);
  }
  if(flag==UNLOAD){
    geomi->loaded=0;
    geomi->display=0;
  }

  file = geomi->file;

  get_geom_header(file,&ntimes);
  if(ntimes<0)return;
  stream = fopen(file,"rb");

  fseek(stream,4,SEEK_CUR);fread(&one,4,1,stream);fseek(stream,4,SEEK_CUR);
  if(one!=1)endianswitch=1;
  FORTREAD(&version,1,stream);

  geomi->ntimes=ntimes;
  NewMemory((void **)&geomi->geomlistinfo,ntimes*sizeof(geomlistdata));

  for(i=0;i<ntimes;i++){
    float times[2];
    geomlistdata *geomlisti;
    int *geom_typeptr,geom_type=0;
    int nverts[4];
    int nvert_s, nvert_d, ntri_s, ntri_d;

    geomlisti = geomi->geomlistinfo+i;
    if(first==1){
      FORTREADBR(times,1,stream);
      FORTREADBR(nverts,4,stream);
      geom_typeptr=&geom_type;
      first=0;
    }
    else{
      nvert_s=0;
      ntri_s=0;
      FORTREADBR(times,2,stream);
      FORTREADBR(nverts+2,2,stream);
      geom_typeptr=(int *)(times+1);
    }
    nvert_s=nverts[0];
    ntri_s=nverts[1];
    nvert_d=nverts[2];
    ntri_d=nverts[3];
    geomlisti->points=NULL;
    if(*geom_typeptr==0&&nvert_s+nvert_d>0){
      int nvert;

      nvert = nvert_s + nvert_d;
      NewMemory((void **)&xyz,3*nvert*sizeof(float));
      NewMemory((void **)&points,nvert*sizeof(point));
      geomlisti->points=points;
      geomlisti->npoints=nvert;
      if(nvert_s>0){
        FORTREADBR(xyz,3*nvert_s,stream);
      }
      if(nvert_d>0){
        FORTREADBR(xyz+3*nvert_s,3*nvert_d,stream);
      }
      for(i=0;i<nvert;i++){
        points[i].xyz[0]=xyz[3*i+0];
        points[i].xyz[1]=xyz[3*i+1];
        points[i].xyz[2]=xyz[3*i+2];
      }
      FREEMEMORY(xyz);
    }
    if(*geom_typeptr==1){
      float tran_rot[8];

      FORTREADBR(tran_rot,8,stream);
      memcpy(geomlisti->translate,tran_rot,3*sizeof(float));
      memcpy(geomlisti->rot0,tran_rot+3,3*sizeof(float));
      memcpy(geomlisti->rot,tran_rot+6,2*sizeof(float));
    }
    if(*geom_typeptr==0&&ntri_s+ntri_d>0){
      int ntris;
      int *surf_ind=NULL,*ijk=NULL;

      ntris = ntri_s+ntri_d;
      NewMemory((void **)&triangles,ntris*sizeof(triangle));
      NewMemory((void **)&ijk,3*ntris*sizeof(int));
      NewMemory((void **)&surf_ind,ntris*sizeof(int));
      geomlisti->triangles=triangles;
      geomlisti->ntriangles=ntris;
      if(ntri_s>0){
        FORTREADBR(ijk,3*ntri_s,stream);
      }
      if(ntri_d>0){
        FORTREADBR(ijk+3*ntri_s,3*ntri_d,stream);
      }
     // if(ntri_s>0){
     //   FORTREADBR(surf_ind,ntri_s,stream);
     // }
     // if(ntri_d>0){
     //   FORTREADBR(ijk+ntri_s,ntri_d,stream);
     // }
      for(i=0;i<ntris;i++){
        triangles[i].points[0]=points+ijk[3*i+0]-1;
        triangles[i].points[1]=points+ijk[3*i+1]-1;
        triangles[i].points[2]=points+ijk[3*i+2]-1;
        //triangles[i].surf=surfaceinfo + surf_ind[i] - 1;
        triangles[i].surf=surfacedefault;
      }
      FREEMEMORY(ijk);
      FREEMEMORY(surf_ind);
    }
  }
  geomi->loaded=1;
  geomi->display=1;
}

/* ------------------ read_geomdata ------------------------ */

void read_geomdata(int ifile, int flag, int *errorcode){
  patch *patchi;
  char *file;
  FILE *stream;
  int one=1, endianswitch=0;
  int returncode;
  int nvert, ntris;
  int ntimes;
  float *xyz=NULL;
  int nxyz=0;
  int *ijk=NULL;
  int i;
  point *points;
  triangle *triangles;
  int first=1;
  float *val_buffer=NULL;
  int nval_buffer=0;
  int nvals;
  float patchmin_global, patchmax_global;
  int n;
  int error;
  FILE_SIZE lenfile;

  // 1
  // time
  // nstatic
  // vals_1, ...vals_nstatic
  // ndynamic
  // vals_1, ... vals_ndyamic

  patchi = patchinfo + ifile;
  if(patchi->filetype!=2)return;
  file = patchi->file;

  patchi->loaded=0;
  patchi->display=0;

  FREEMEMORY(patchi->geom_nstatics);
  FREEMEMORY(patchi->geom_ndynamics);
  FREEMEMORY(patchi->geom_ivals_static);
  FREEMEMORY(patchi->geom_ivals_dynamic);
  FREEMEMORY(patchi->geom_vals);
  FREEMEMORY(patchi->geom_ivals);
  FREEMEMORY(patchi->geom_times);
  if(flag==UNLOAD){
    plotstate=getplotstate(DYNAMIC_PLOTS);
    update_patchtype();
    update_unit_defs();
    updatetimes();
    return;
  }

  //get_geomdata_header(file,&ntimes,&nvals);
  endian = getendian();
  lenfile = strlen(file);

  FORTgetembeddatasize(file, &endian, &ntimes, &nvals, &error, lenfile);

  if(nvals>0){
    NewMemory((void **)&patchi->geom_nstatics,ntimes*sizeof(int));
    NewMemory((void **)&patchi->geom_ndynamics,ntimes*sizeof(int));
    NewMemory((void **)&patchi->geom_times,ntimes*sizeof(float));
    NewMemory((void **)&patchi->geom_ivals_static,ntimes*sizeof(int *));
    NewMemory((void **)&patchi->geom_ivals_dynamic,ntimes*sizeof(int *));
    NewMemory((void **)&patchi->geom_vals,nvals*sizeof(float));
    NewMemory((void **)&patchi->geom_ivals,nvals*sizeof(char));
  }
  FORTgetembeddata(file, &endian, &ntimes, &nvals, patchi->geom_times, 
    patchi->geom_nstatics, patchi->geom_ndynamics, patchi->geom_vals, &error, lenfile);

  init_histogram(patchi->histogram);

  update_histogram(patchi->geom_vals,nvals,patchi->histogram);

  complete_histogram(patchi->histogram);

  patchi->geom_ntimes=ntimes;
  patchi->geom_nvals=nvals;
  patchi->geom_ivals_static[0] = patchi->geom_ivals;
  patchi->geom_ivals_dynamic[0] = patchi->geom_ivals_static[0]+patchi->geom_nstatics[0];
  for(i=1;i<ntimes;i++){
    patchi->geom_ivals_static[i] = patchi->geom_ivals_dynamic[i-1]+patchi->geom_ndynamics[i-1];
    patchi->geom_ivals_dynamic[i] = patchi->geom_ivals_static[i] + patchi->geom_nstatics[i];
  }
  if(colorlabelpatch!=NULL){
    for(n=0;n<MAXRGB;n++){
      FREEMEMORY(colorlabelpatch[n]);
    }
    FREEMEMORY(colorlabelpatch);
  }
  if(NewMemory((void **)&colorlabelpatch,MAXRGB*sizeof(char *))==0){
    readpatch(ifile,UNLOAD,&error);
    return;
  }
  for(n=0;n<MAXRGB;n++){
    colorlabelpatch[n]=NULL;
  }
  for(n=0;n<nrgb;n++){
    if(NewMemory((void **)&colorlabelpatch[n],11)==0){
      readpatch(ifile,UNLOAD,&error);
      return;
    }
  }
  getBoundaryColors3(patchi,patchi->geom_vals, patchi->geom_nvals, patchi->geom_ivals,
    setpatchmin,&patchmin, setpatchmax,&patchmax, 
    &patchmin_global, &patchmax_global,
    nrgb, colorlabelpatch,patchi->scale,boundarylevels256,
    &patchi->extreme_min,&patchi->extreme_max);
  FREEMEMORY(patchi->geom_vals);
  patchi->loaded=1;
  patchi->display=1;
  ipatchtype=getpatchtype(patchinfo+ifile);
  plotstate=getplotstate(DYNAMIC_PLOTS);
  update_patchtype();
  update_unit_defs();
  updatetimes();
  update_framenumber(1);
}

/* ------------------ draw_geomdata ------------------------ */

void draw_geomdata(patch *patchi){
  int i;
  float black[]={0.0,0.0,0.0,1.0};
  float blue[]={0.0,0.0,1.0,1.0};
  float skinny_color[]={1.0,0.0,0.0,1.0};
  float *last_color=NULL;

  for(i=0;i<ngeominfo;i++){
    geomdata *geomi;
    geomlistdata *geomlisti;
    int ntris,npoints;
    int j;
    float *color;

    geomi = geominfo + i;
    if(geomi->display==0||geomi->loaded==0)continue;
    geomlisti = geomi->geomlistinfo;
    ntris = geomlisti->ntriangles;
    npoints = geomlisti->npoints;

    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glEnable(GL_COLOR_MATERIAL);

    glPushMatrix();
    glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);
    glTranslatef(-xbar0,-ybar0,-zbar0);
    glBegin(GL_TRIANGLES);
    if(smoothtrinormal==0){
      for(j=0;j<ntris;j++){
        float *xyzptr[3];
        float *xyznorm;
        triangle *trianglei;
        int color_index;

        trianglei = geomlisti->triangles+j;
       
        xyznorm=trianglei->normal;
        glNormal3fv(xyznorm);

        color_index = patchi->geom_ival_static[j];
        color=rgb_patch+4*color_index;
        glColor3fv(color);

        xyzptr[0] = trianglei->points[0]->xyz;
        glVertex3fv(xyzptr[0]);

        xyzptr[1] = trianglei->points[1]->xyz;
        glVertex3fv(xyzptr[1]);

        xyzptr[2] = trianglei->points[2]->xyz;
        glVertex3fv(xyzptr[2]);
      }
    }
    else{
      for(j=0;j<ntris;j++){
        float *xyzptr[3];
        float *xyznorm;
        triangle *trianglei;
        int color_index;

        trianglei = geomlisti->triangles+j;
       
        color_index = patchi->geom_ival_static[j];
        color=rgb_patch+4*color_index;
        glColor3fv(color);

        xyznorm = trianglei->points[0]->norm;
        glNormal3fv(xyznorm);
        xyzptr[0] = trianglei->points[0]->xyz;
        glVertex3fv(xyzptr[0]);

        xyznorm = trianglei->points[1]->norm;
        glNormal3fv(xyznorm);
        xyzptr[1] = trianglei->points[1]->xyz;
        glVertex3fv(xyzptr[1]);

        xyznorm = trianglei->points[2]->norm;
        glNormal3fv(xyznorm);
        xyzptr[2] = trianglei->points[2]->xyz;
        glVertex3fv(xyzptr[2]);
      }
    }
    glEnd();
    glPopMatrix();
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
  }

}

/* ------------------ compare_trdata ------------------------ */

int compare_transparent_triangles( const void *arg1, const void *arg2 ){
  triangle *tri, *trj;

  tri = (triangle *)arg1; 
  trj = (triangle *)arg2;

  if(tri->distance<trj->distance)return -1;
  if(tri->distance>trj->distance)return 1;
  return 0;
}

/* ------------------ sort_triangles ------------------------ */

void Sort_Embedded_Geometry(float *mm){
  int itri;
  int newflag;
  int i;
  int count_transparent,count_all;

  CheckMemory;
  count_transparent=0;
  count_all=0;
  ntransparent_triangles=count_transparent;
  nopaque_triangles=count_all-count_transparent;
  for(i=0;i<ngeominfo;i++){
    geomlistdata *geomlisti;
    int ntris,npoints;
    int j;
    float *color;
    geomdata *geomi;

    geomi = geominfo + i;
    if(geomi->loaded==0||geomi->display==0)continue;
    geomlisti = geomi->geomlistinfo;

    count_all+=geomlisti->ntriangles;
    if(use_transparency_data==0)continue;
    for(j=0;j<geomlisti->ntriangles;j++){
      triangle *tri;
      float xyz[3];
      float *xyz1, *xyz2, *xyz3;
      float xyzeye[3];

      tri = geomlisti->triangles + j;
      if(hilight_skinny==1&&tri->skinny==1)continue;
      if(tri->surf->color[3]>=1.0)continue;
      count_transparent++;
      if(sort_embedded_geometry==0)continue;
      xyz1 = tri->points[0]->xyz;
      xyz2 = tri->points[1]->xyz;
      xyz3 = tri->points[2]->xyz;
      xyz[0] = xyz1[0]+xyz2[0]+xyz3[0];
      xyz[1] = xyz1[1]+xyz2[1]+xyz3[1];
      xyz[2] = xyz1[2]+xyz2[2]+xyz3[2];

      xyzeye[0] = mm[0]*xyz[0] + mm[4]*xyz[1] +   mm[8]*xyz[2] + mm[12];
      xyzeye[1] = mm[1]*xyz[0] + mm[5]*xyz[1] +   mm[9]*xyz[2] + mm[13];
      xyzeye[2] = mm[2]*xyz[0] + mm[6]*xyz[1] +  mm[10]*xyz[2] + mm[14];
      xyzeye[0]/=mscale[0];
      xyzeye[1]/=mscale[1];
      xyzeye[2]/=mscale[2];
      tri->distance=xyzeye[0]*xyzeye[0]+xyzeye[1]*xyzeye[1]+xyzeye[2]*xyzeye[2];
      CheckMemory;
    }
  }
  CheckMemory;
  if(count_all==0)return;
  FREEMEMORY(alltriangles);
  NewMemory((void **)&alltriangles,count_all*sizeof(triangle **));
  transparent_triangles=alltriangles;
  opaque_triangles=alltriangles+count_transparent;
  ntransparent_triangles=count_transparent;
  nopaque_triangles=count_all-count_transparent;
  count_transparent=0;
  count_all=0;
  for(i=0;i<ngeominfo;i++){
    geomlistdata *geomlisti;
    int j;
    float *color;
    geomdata *geomi;

    geomi = geominfo + i;
    if(geomi->loaded==0||geomi->display==0)continue;
    geomlisti = geomi->geomlistinfo;
    for(j=0;j<geomlisti->ntriangles;j++){
      triangle *tri;

      tri = geomlisti->triangles + j;
      if(use_transparency_data==0||(hilight_skinny==1&&tri->skinny==1)||tri->surf->color[3]>=1.0){
        opaque_triangles[count_all++]=tri;
      }
      else{
        transparent_triangles[count_transparent++]=tri;
      }
    }
  }
  if(sort_embedded_geometry==1&&ntransparent_triangles>0){
    qsort((isotri **)transparent_triangles,(size_t)ntransparent_triangles,sizeof(triangle **),compare_transparent_triangles);
  }
}


