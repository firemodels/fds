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
#include "flowfiles.h"
#include "MALLOC.h"
#include "smokeviewvars.h"
#include "translate.h"
#include "update.h"
#include "string_util.h"
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

  for(i=0;i<ntrilistinfo;i++){
    trilistdata *trilisti;
    pointlistdata *pointlisti;
    point **points;
    int j;
    int ndups=0,nused=0,nskinny=0;

    trilisti = trilistinfo + i;
    pointlisti = pointlistinfo + i;

    if(pointlisti->npoints>0){
      NewMemory((void **)&points,pointlisti->npoints*sizeof(point *));
      for(j=0;j<pointlisti->npoints;j++){
        points[j]=pointlisti->points+j;
        points[j]->nused=0;
      }
      qsort(points,pointlisti->npoints,sizeof(point *),compare_verts);
      for(j=1;j<pointlisti->npoints;j++){
        if(compare_verts(points[j-1],points[j])==0)ndups++;
      }
      for(j=0;j<trilisti->ntriangles;j++){
        triangle *trii;
        float min_angle;

        trii = trilisti->triangles + j;
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
      for(j=0;j<pointlisti->npoints;j++){
        if(points[j]->nused>0)nused++;
      }
      printf("Face/Vertex Summary\n");
      printf("      Faces: %i\n",trilisti->ntriangles);
      printf(" slim faces: %i\n",nskinny);
      printf("   Vertices: %i\n",pointlisti->npoints);
      printf("     unused: %i\n",pointlisti->npoints-nused);
      printf(" duplicates: %i\n\n",ndups);
      FREEMEMORY(points);
    }
  }
}

/* ------------------ draw_geom ------------------------ */

void draw_geom(void){
  int i;
  float black[]={0.0,0.0,0.0,1.0};
  float blue[]={0.0,0.0,1.0,1.0};
  float skinny_color[]={1.0,0.0,0.0,1.0};
  float *last_color=NULL;

  for(i=0;i<ntrilistinfo;i++){
    trilistdata *trilisti;
    pointlistdata *pointlisti;
    int ntris,npoints;
    int j;
    float *color;

    trilisti = trilistinfo + i;
    pointlisti = pointlistinfo + i;
    ntris = trilisti->ntriangles;
    npoints = pointlisti->npoints;
    if(patchembedded==0&&showtrisurface==1){
      glEnable(GL_LIGHTING);
      glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
      glEnable(GL_COLOR_MATERIAL);
      glBegin(GL_TRIANGLES);
      if(smoothtrinormal==0){
        for(j=0;j<ntris;j++){
          float *xyzptr[3];
          float *xyznorm;
          triangle *trianglei;

          trianglei = trilisti->triangles+j;
       
          xyznorm=trianglei->normal;
          glNormal3fv(xyznorm);

          if(hilight_skinny==1&&trianglei->skinny==1){
            color=skinny_color;
          }
          else{
            color = trianglei->surf->color;
          }
          if(color!=last_color){
            glColor3fv(color);
            last_color=color;
          }

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

          trianglei = trilisti->triangles+j;
       
          if(hilight_skinny==1&&trianglei->skinny==1){
            color=skinny_color;
          }
          else{
            color = trianglei->surf->color;
          }
          if(color!=last_color){
            glColor3fv(color);
            last_color=color;
          }
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
      glDisable(GL_COLOR_MATERIAL);
      glDisable(GL_LIGHTING);
    }
    if(showtrioutline==1){
      glBegin(GL_LINES);
      for(j=0;j<ntris;j++){
        float *xyzptr[3];
        float *xyznorm;
        triangle *trianglei;

        trianglei = trilisti->triangles+j;
       
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
    }
    if(showtrinormal==1){
      if(smoothtrinormal==0){
        glBegin(GL_LINES);
        for(j=0;j<ntris;j++){
          float *p1, *p2, *p3;
          float *xyznorm;
          triangle *trianglei;
          float xyz1[3], xyz2[3];

          trianglei = trilisti->triangles+j;
       
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

          trianglei = trilisti->triangles+j;
       
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
      }
      if(smoothtrinormal==1){
        glBegin(GL_LINES);
        for(j=0;j<npoints;j++){
          float *xyznorm;
          point *pointi;
          float *xyz1, xyz2[3];

          pointi = pointlisti->points+j;
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

          pointi = pointlisti->points+j;
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
      }
    }
  }
}

/* ------------------ update_triangles ------------------------ */

void update_triangles(void){
  int j;

  for(j=0;j<ntrilistinfo;j++){
    trilistdata *trilisti;
    pointlistdata *pointlisti;
    float *xyzptr[3];
    float *xyznorm;
    float *xyz;
    int i;

    trilisti = trilistinfo + j;
    pointlisti = pointlistinfo + j;
    
    for(i=0;i<pointlisti->npoints;i++){
      xyz = pointlisti->points[i].xyz;
      xyz[0] = (xyz[0] - xbar0)/xyzmaxdiff;
      xyz[1] = (xyz[1] - ybar0)/xyzmaxdiff;
      xyz[2] = (xyz[2] - zbar0)/xyzmaxdiff;
      xyz += 3;
    }
    
    for(i=0;i<trilisti->ntriangles;i++){
      triangle *trianglei;

      trianglei = trilisti->triangles+i;
      if(trilisti->triangles[i].fdsnorm==0){
        xyzptr[0] = trianglei->points[0]->xyz;
        xyzptr[1] = trianglei->points[1]->xyz;
        xyzptr[2] = trianglei->points[2]->xyz;
        xyznorm = trianglei->normal;
        CalcTriNormal(xyzptr[0],xyzptr[1],xyzptr[2],xyznorm);
      }
    }

    for(i=0;i<pointlisti->npoints;i++){
      point *pointi;

      pointi = pointlisti->points + i;
      pointi->ntriangles=0;
      pointi->itriangle=0;
    }
    for(i=0;i<trilisti->ntriangles;i++){
      triangle *trianglei;

      trianglei = trilisti->triangles+i;
      trianglei->points[0]->ntriangles++;
      trianglei->points[1]->ntriangles++;
      trianglei->points[2]->ntriangles++;
    }
    for(i=0;i<pointlisti->npoints;i++){
      point *pointi;

      pointi = pointlisti->points + i;
      if(pointi->ntriangles>0){
        NewMemory((void **)&pointi->triangles,pointi->ntriangles*sizeof(triangle *));
      }
    }
    for(i=0;i<trilisti->ntriangles;i++){
      triangle *trianglei;
      point *pointi;

      trianglei = trilisti->triangles+i;
      pointi = trianglei->points[0];
      pointi->triangles[pointi->itriangle++]=trianglei;
      pointi = trianglei->points[1];
      pointi->triangles[pointi->itriangle++]=trianglei;
      pointi = trianglei->points[2];
      pointi->triangles[pointi->itriangle++]=trianglei;
    }
    for(i=0;i<pointlisti->npoints;i++){
      point *pointi;
      int k;
      float *norm;

      pointi = pointlisti->points + i;
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
/*
typedef struct {
  float xyz[3],norm[3];
  int itriangle,ntriangles,nused;
  struct _triangle **triangles;
} point;

typedef struct _triangle {
  struct _surface *surf;
  int interior;
  float *color;
  int fdsnorm,skinny;
  point *points[3];
  float normal[3];
} triangle;

typedef struct {
  int npoints;
  point *points;
} pointlistdata;

typedef struct {
  int ntriangles;
  triangle *triangles;
} trilistdata;

 typedef struct {
  char *file;
  pointlistdata *pointlistinfo;
  trilistdata *trilistinfo;
  float *times;
  int ntimes;
} geomdata;
*/

/* ------------------ get_geom_header ------------------------ */

void get_geom_header(char *file, int *ntimes){
  FILE *stream;
  int one=1,endianswitch=0;
  int ntris,nvert;
  float time;
  int first=1;
  int nt;
  int returncode;

  stream = fopen(file,"r");
  if(stream==NULL){
    *ntimes=-1;
    return;
  }
  fseek(stream,4,SEEK_CUR);fread(&one,4,1,stream);fseek(stream,4,SEEK_CUR);
  if(one!=1)endianswitch=1;
  nt=-1;
  for(;;){
    if(first!=1){
      FORTREADBR(&time,1,stream);
    }
    first=0;
    FORTREADBR(&nvert,1,stream);
    if(nvert!=0)fseek(stream,4+3*nvert*4+4,SEEK_CUR);    
    FORTREADBR(&ntris,1,stream);
    if(ntris!=0)fseek(stream,4+3*ntris*4+4,SEEK_CUR);    
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

/* ------------------ read_geom ------------------------ */

void read_geom(int ifile, int flag, int *errorcode){
  geomdata *geomi;
  char *file;
  FILE *stream;
  int one=1, endianswitch=0;
  int returncode;
  int nvert, ntris;
  int ntimes;
  float *xyz=NULL;
  int nxyz=0;
  int *ijk=NULL,nijk=0;
  int i;
  point *points;
  triangle *triangles;
  int first=1;

  // 1
  //   static geometry
  // nverts
  // x1 y1 z1 ... xnverts ynverts znverts
  // nfaces
  // i1 j1 k1 ... infaces jnfaces knfaces
  //  dynamic geometry (one entry for each time step)
  // time
  // nverts
  // x1 y1 z1 ... xnverts ynverts znverts
  // nfaces
  // i1 j1 k1 ... infaces jnfaces knfaces

  geomi = geominfo + ifile;
  file = geomi->file;

  get_geom_header(file,&ntimes);
  if(ntimes<0)return;
  stream = fopen(file,"r");

  fseek(stream,4,SEEK_CUR);fread(&one,4,1,stream);fseek(stream,4,SEEK_CUR);
  if(one!=1)endianswitch=1;

  geomi->ntimes=ntimes;
  NewMemory((void **)&geomi->pointlistinfo,(ntimes+1)*sizeof(pointlistdata));
  NewMemory((void **)&geomi->trilistinfo,(ntimes+1)*sizeof(trilistdata));
  NewMemory((void **)&geomi->times,(ntimes+1)*sizeof(float));

  for(i=0;i<=ntimes;i++){
    float time;
    pointlistdata *pointlisti;
    trilistdata *trilisti;

    pointlisti = geomi->pointlistinfo+i;
    trilisti = geomi->trilistinfo+i;
    if(first!=1)FORTREAD(&time,1,stream);
    first=0;
    FORTREAD(&nvert,1,stream);
    FORTREAD(&nvert,1,stream);
    if(nvert>0){
      FREEMEMORY(xyz);
      NewMemory((void **)&xyz,3*nvert*sizeof(float));
      NewMemory((void **)&points,nvert*sizeof(point));
      pointlisti->points=points;
      FORTREAD(xyz,3*nvert,stream);
      for(i=0;i<nvert;i++){
        points[i].xyz[0]=xyz[3*i+0];
        points[i].xyz[1]=xyz[3*i+1];
        points[i].xyz[2]=xyz[3*i+2];
      }
    }
    FORTREAD(&ntris,1,stream);
    if(ntris>0){
      FREEMEMORY(ijk);
      NewMemory((void **)&ijk,3*ntris*sizeof(int));
      NewMemory((void **)&triangles,ntris*sizeof(triangle));
      trilisti->triangles=triangles;
      FORTREAD(ijk,3*ntris,stream);
      for(i=0;i<nvert;i++){
        triangles[i].points[0]=points+ijk[3*i+0];
        triangles[i].points[1]=points+ijk[3*i+1];
        triangles[i].points[2]=points+ijk[3*i+2];
      }
    }
  }
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
  int *ijk=NULL,nijk=0;
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

  FORTgetembedsize(file, &endian, &ntimes, &nvals, &error, lenfile);

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

  for(i=0;i<ntrilistinfo;i++){
    trilistdata *trilisti;
    pointlistdata *pointlisti;
    int ntris,npoints;
    int j;
    float *color;

    trilisti = trilistinfo + i;
    pointlisti = pointlistinfo + i;
    ntris = trilisti->ntriangles;
    npoints = pointlisti->npoints;

    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glEnable(GL_COLOR_MATERIAL);
    glBegin(GL_TRIANGLES);
    if(smoothtrinormal==0){
      for(j=0;j<ntris;j++){
        float *xyzptr[3];
        float *xyznorm;
        triangle *trianglei;
        int color_index;

        trianglei = trilisti->triangles+j;
       
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

        trianglei = trilisti->triangles+j;
       
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
    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
  }

}
