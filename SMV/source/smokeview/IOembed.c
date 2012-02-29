// $Date$ 
// $Revision$
// $Author$

// svn revision character string
char IOembed_revision[]="$Revision$";

#include "options.h"
#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "string_util.h"
#include "smv_endian.h"
#include "update.h"
#include "smokeviewvars.h"

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
  float angle_local;
  float arg;

  arg = (d2*d2+d3*d3-d1*d1)/(2.0*d2*d3);
  if(arg<-1.0)arg=-1.0;
  if(arg>1.0)arg=1.0;
  angle_local = acos(arg)*180.0/(4.0*atan(1.0));
  return angle_local;
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

  for(i=0;i<ngeominfoptrs;i++){
    geomdata *geomi;
    geomlistdata *geomlisti;
    point **points;
    int j;
    int ndups=0,nused=0,nskinny=0;

    geomi = geominfoptrs[i];
    geomlisti = geomi->geomlistinfo;

    if(geomlisti->npoints>0){
      NewMemory((void **)&points,geomlisti->npoints*sizeof(point *));
      for(j=0;j<geomlisti->npoints;j++){
        points[j]=geomlisti->points+j;
        points[j]->nused=0;
      }
      for(j=0;j<geomlisti->ntriangles;j++){
        triangle *trii;

        trii = geomlisti->triangles + j;
        trii->points[0]->nused=0;
        trii->points[1]->nused=0;
        trii->points[2]->nused=0;
      }
      qsort(points,geomlisti->npoints,sizeof(point *),compare_verts);
      for(j=1;j<geomlisti->npoints;j++){
        if(compare_verts(points[j-1],points[j])==0)ndups++;
      }
      for(j=0;j<geomlisti->ntriangles;j++){
        triangle *trii;

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

void draw_geom(int flag, int frameflag){
  int i;
  float black[]={0.0,0.0,0.0,1.0};
  float blue[]={0.0,0.0,1.0,1.0};
  float skinny_color[]={1.0,0.0,0.0,1.0};
  float *last_color=NULL;
  float last_transparent_level=-1.0;

  if(patchembedded==0&&showtrisurface==1&&frameflag==0){
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
    
    if(cullfaces==1)glDisable(GL_CULL_FACE);
    glEnable(GL_NORMALIZE);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,iso_specular);
    glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,iso_shininess);
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
      float transparent_level_local;

      trianglei = tris[i];
      if(hilight_skinny==1&&trianglei->skinny==1){
        color=skinny_color;
        transparent_level_local=1.0;
      }
      else{
        color = trianglei->surf->color;
        transparent_level_local=trianglei->surf->transparent_level;
      }
      if(color!=last_color||ABS(last_transparent_level-transparent_level_local)>0.001){
        glColor4f(color[0],color[1],color[2],transparent_level_local);
        last_color=color;
        last_transparent_level=transparent_level_local;
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
    if(flag==DRAW_TRANSPARENT){
      if(use_transparency_data==1)transparentoff();
      return;
    }
    if(cullfaces==1)glEnable(GL_CULL_FACE);
  }

#define VECFACTOR 0.01

  for(i=0;i<ngeominfoptrs;i++){
    geomdata *geomi;
    geomlistdata *geomlisti;
    int ntris,npoints;
    int j;
    float *color;

    geomi = geominfoptrs[i];
    if(geomi->loaded==0||geomi->display==0)continue;
    if(frameflag==0){
      geomlisti = geomi->geomlistinfo-1;
    }
    else{
      geomlisti = geomi->geomlistinfo+geomi->iframe;
    }
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
    if(showtripoints==1){
      glPushMatrix();
      glScalef(1.0/xyzmaxdiff,1.0/xyzmaxdiff,1.0/xyzmaxdiff);
      glTranslatef(-xbar0,-ybar0,-zbar0);
      glPointSize(6.0);
      glBegin(GL_POINTS);
      for(j=0;j<geomlisti->npoints;j++){
        point *pointi;

        pointi = geomlisti->points+j;
        color = black;
        glColor3fv(color);
        glVertex3fv(pointi->xyz);
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
          xyz2[0] = xyz1[0] + VECFACTOR*xyzmaxdiff*xyznorm[0];
          xyz2[1] = xyz1[1] + VECFACTOR*xyzmaxdiff*xyznorm[1];
          xyz2[2] = xyz1[2] + VECFACTOR*xyzmaxdiff*xyznorm[2];

          glColor3fv(blue);
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
          xyz2[0] = xyz1[0] + VECFACTOR*xyzmaxdiff*xyznorm[0];
          xyz2[1] = xyz1[1] + VECFACTOR*xyzmaxdiff*xyznorm[1];
          xyz2[2] = xyz1[2] + VECFACTOR*xyzmaxdiff*xyznorm[2];

          glColor3fv(blue);
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

          xyz2[0] = xyz1[0] + VECFACTOR*xyzmaxdiff*xyznorm[0];
          xyz2[1] = xyz1[1] + VECFACTOR*xyzmaxdiff*xyznorm[1];
          xyz2[2] = xyz1[2] + VECFACTOR*xyzmaxdiff*xyznorm[2];

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

          xyz2[0] = xyz1[0] + VECFACTOR*xyzmaxdiff*xyznorm[0];
          xyz2[1] = xyz1[1] + VECFACTOR*xyzmaxdiff*xyznorm[1];
          xyz2[2] = xyz1[2] + VECFACTOR*xyzmaxdiff*xyznorm[2];

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
  int ii,j;
   
  for(j=0;j<ngeominfoptrs;j++){
    geomlistdata *geomlisti;
    geomdata *geomi;
    float *xyzptr[3];
    float *xyznorm;
    int i;

    geomi = geominfoptrs[j];
    if(geomi->loaded==0||geomi->display==0)continue;

    for(ii=-1;ii<geomi->ntimes;ii++){
      geomlisti = geomi->geomlistinfo+ii;
    
      for(i=0;i<geomlisti->ntriangles;i++){
        triangle *trianglei;

        trianglei = geomlisti->triangles+i;

        xyzptr[0] = trianglei->points[0]->xyz;
        xyzptr[1] = trianglei->points[1]->xyz;
        xyzptr[2] = trianglei->points[2]->xyz;
        xyznorm = trianglei->normal;
        CalcTriNormal(xyzptr[0],xyzptr[1],xyzptr[2],xyznorm);
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
}

#define FORTREAD(var,count,STREAM) fseek(STREAM,4,SEEK_CUR);\
                           returncode=fread(var,4,count,STREAM);\
                           if(returncode!=count)returncode=0;\
                           if(endianswitch==1&&returncode!=0)endian_switch(var,count);\
                           fseek(STREAM,4,SEEK_CUR)

#define FORTREADBR(var,count,STREAM) FORTREAD(var,(count),STREAM);if(returncode==0)break;

/* ------------------ read_geom_header ------------------------ */

void read_geom_header(geomdata *geomi, int *ntimes_local){
  FILE *stream;
  int one=0,endianswitch=0;
  int nvertfaces[2];
  float times_local[2];
  int nt;
  int returncode;
  int version;
  int nfloat_vals, nint_vals;
  int *int_vals;
  float *float_vals;
  int nverts, ntris;

  stream = fopen(geomi->file,"rb");
  if(stream==NULL){
    *ntimes_local=-1;
    return;
  }
  fseek(stream,4,SEEK_CUR);fread(&one,4,1,stream);fseek(stream,4,SEEK_CUR);
  if(one!=1)endianswitch=1;
  FORTREAD(&version,1,stream);

// floating point header

  FORTREAD(&nfloat_vals,1,stream);
  if(nfloat_vals>0){
    NewMemory((void **)&float_vals,nfloat_vals*sizeof(float));
    FORTREAD(float_vals,nfloat_vals,stream);
    geomi->float_vals=float_vals;
    geomi->nfloat_vals=nfloat_vals;
  }

// integer header

  FORTREAD(&nint_vals,1,stream);
  if(nint_vals>0){
    NewMemory((void **)&int_vals,nint_vals*sizeof(float));
    FORTREAD(int_vals,nint_vals,stream);
    geomi->int_vals=int_vals;
    geomi->nint_vals=nint_vals;
  }

// static verts

  FORTREAD(&nvertfaces,2,stream);
  nverts=nvertfaces[0];
  if(nverts>0){
    fseek(stream,4+3*nverts*4+4,SEEK_CUR);
  }

// static triangles

  ntris=nvertfaces[1];
  if(ntris>0){
    fseek(stream,4+3*ntris*4+4,SEEK_CUR);
    fseek(stream,4+ntris*4+4,SEEK_CUR);
  }

  nt=0;
  for(;;){
    int geom_type;

    FORTREADBR(times_local,2,stream);
    geom_type = *((int *)(times_local+1));

    if(geom_type==0){
      FORTREADBR(&nvertfaces,2,stream);

      nverts=nvertfaces[0];
      ntris=nvertfaces[1];

      if(nverts>0){
        fseek(stream,4+3*nverts*4+4,SEEK_CUR);
      }
      if(ntris>0){
        fseek(stream,4+3*ntris*4+4,SEEK_CUR);
        fseek(stream,4+ntris*4+4,SEEK_CUR);
      }
    }
    else{
      fseek(stream,4+8*4+4,SEEK_CUR);
    }

    nt++;
  }
  *ntimes_local=nt;
  fclose(stream);
}

/* ------------------ get_geomdata_header ------------------------ */

void get_geomdata_header(char *file, int *ntimes_local, int *nvals){
  FILE *stream;
  int one=1,endianswitch=0;
  int nface_static,nface_dynamic;
  float time_local;
  int nt,nv;
  int returncode;

  stream = fopen(file,"r");
  if(stream==NULL){
    *ntimes_local=-1;
    return;
  }
  fseek(stream,4,SEEK_CUR);fread(&one,4,1,stream);fseek(stream,4,SEEK_CUR);
  if(one!=1)endianswitch=1;
  nt=-1;
  nv=0;
  for(;;){
    FORTREADBR(&time_local,1,stream);
    FORTREADBR(&nface_static,1,stream);
    if(nface_static!=0)fseek(stream,4+nface_static*4+4,SEEK_CUR);    
    FORTREADBR(&nface_dynamic,1,stream);
    if(nface_dynamic!=0)fseek(stream,4+nface_dynamic*4+4,SEEK_CUR);    
    nt++;
    nv+=(nface_static+nface_dynamic);
  }
  *ntimes_local=nt;
  *nvals=nv;
  fclose(stream);
}

/* ------------------ read_all_geom ------------------------ */

void read_all_geom(void){
  int i, errorcode;

  for(i=0;i<ngeominfo;i++){
    geomdata *geomi;

    geomi = geominfo + i;
    read_geom(geomi,LOAD,&errorcode);
  }
}

/* ------------------ read_geom ------------------------ */

void read_geom(geomdata *geomi, int flag, int *errorcode){
  FILE *stream;
  int one=1, endianswitch=0;
  int returncode;
  int ntimes_local;
  float *xyz=NULL;
  int i;
  point *points;
  triangle *triangles;
  int version;
  int nvertfaces[2];
  int nfloat_vals, nint_vals;

  if(geomi->geomlistinfo!=NULL){
    for(i=-1;i<geomi->ntimes;i++){
      geomlistdata *geomlisti;

      geomlisti = geomi->geomlistinfo+i;
      FREEMEMORY(geomlisti->points);
      FREEMEMORY(geomlisti->triangles);
    }  
  }
  FREEMEMORY(geomi->times);
  FREEMEMORY(geomi->geomlistinfo_0);
  FREEMEMORY(geomi->float_vals);
  FREEMEMORY(geomi->int_vals);
  geomi->nfloat_vals=0;
  geomi->nint_vals=0;

  if(flag==UNLOAD){
    geomi->loaded=0;
    geomi->display=0;
    return;
  }

  read_geom_header(geomi,&ntimes_local);
  if(ntimes_local<0)return;
  stream = fopen(geomi->file,"rb");
  if(stream==NULL)return;

  fseek(stream,4,SEEK_CUR);fread(&one,4,1,stream);fseek(stream,4,SEEK_CUR);
  if(one!=1)endianswitch=1;
  FORTREAD(&version,1,stream);

  FORTREAD(&nfloat_vals,1,stream);
  if(nfloat_vals>0)fseek(stream,4+nfloat_vals*4+4,SEEK_CUR);
  FORTREAD(&nint_vals,1,stream);
  if(nint_vals>0)fseek(stream,4+nint_vals*4+4,SEEK_CUR);

  geomi->ntimes=ntimes_local;
  geomi->iframe=0;
  NewMemory((void **)&geomi->geomlistinfo_0,(ntimes_local+1)*sizeof(geomlistdata));
  geomi->geomlistinfo=geomi->geomlistinfo_0+1;
  NewMemory((void **)&geomi->times,ntimes_local*sizeof(float));

  for(i=-1;i<ntimes_local;i++){
    float times_local[2];
    geomlistdata *geomlisti;
    int geom_type;
    int nverts, ntris;

    geom_type=0;
    geomlisti = geomi->geomlistinfo+i;
    geomlisti->points=NULL;
    geomlisti->triangles=NULL;
    geomlisti->npoints=0;
    geomlisti->ntriangles=0;
    if(i>=0){
      FORTREADBR(times_local,2,stream);
      geom_type=*((int *)(times_local+1));
      geomi->times[i]=times_local[0];
    }
    if(geom_type==0){
      FORTREADBR(nvertfaces,2,stream);
      nverts=nvertfaces[0];
      ntris=nvertfaces[1];
      if(i>=0){
        printf("time=%.2f triangles: %i\n",times_local[0],ntris);
      }
      else{
      }
      if(nverts>0){
        int ii;

        if(i<0)printf("static geometry\n");
        NewMemory((void **)&xyz,3*nverts*sizeof(float));
        NewMemory((void **)&points,nverts*sizeof(point));
        geomlisti->points=points;
        geomlisti->npoints=nverts;
        if(nverts>0){
          FORTREADBR(xyz,3*nverts,stream);
        }
        for(ii=0;ii<nverts;ii++){
          points[ii].xyz[0]=xyz[3*ii];
          points[ii].xyz[1]=xyz[3*ii+1];
          points[ii].xyz[2]=xyz[3*ii+2];
        }
        FREEMEMORY(xyz);
      }
      if(ntris>0){
        int *surf_ind=NULL,*ijk=NULL;
        int ii;

        NewMemory((void **)&triangles,ntris*sizeof(triangle));
        NewMemory((void **)&ijk,3*ntris*sizeof(int));
        NewMemory((void **)&surf_ind,ntris*sizeof(int));
        geomlisti->triangles=triangles;
        geomlisti->ntriangles=ntris;
        if(ntris>0){
          FORTREADBR(ijk,3*ntris,stream);
        }
        if(ntris>0){
          FORTREADBR(surf_ind,ntris,stream);
        }
        for(ii=0;ii<ntris;ii++){
          triangles[ii].points[0]=points+ijk[3*ii]-1;
          triangles[ii].points[1]=points+ijk[3*ii+1]-1;
          triangles[ii].points[2]=points+ijk[3*ii+2]-1;
          triangles[ii].surf=surfinfo + surf_ind[ii]+nsurfinfo;
        }
        FREEMEMORY(ijk);
        FREEMEMORY(surf_ind);
      }
    }
    else{  // geom_type==1
      float tran_rot[8];

      FORTREADBR(tran_rot,8,stream);
      memcpy(geomlisti->translate,tran_rot,3*sizeof(float));
      memcpy(geomlisti->rot0,tran_rot+3,3*sizeof(float));
      memcpy(geomlisti->rot,tran_rot+6,2*sizeof(float));
      if(i>=0)printf("time=%.2f",times_local[0]);
    }
  }
  geomi->loaded=1;
  geomi->display=1;
}

/* ------------------ read_geomdata ------------------------ */

void read_geomdata(int ifile, int flag, int *errorcode){
  patch *patchi;
  char *file;
  int ntimes_local;
  int i;
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
  endian_smv = getendian();
  lenfile = strlen(file);

  FORTgetembeddatasize(file, &endian_smv, &ntimes_local, &nvals, &error, lenfile);

  if(nvals>0&&ntimes_local>0){
    NewMemory((void **)&patchi->geom_nstatics,ntimes_local*sizeof(int));
    NewMemory((void **)&patchi->geom_ndynamics,ntimes_local*sizeof(int));
    NewMemory((void **)&patchi->geom_times,ntimes_local*sizeof(float));
    NewMemory((void **)&patchi->geom_ivals_static,ntimes_local*sizeof(int *));
    NewMemory((void **)&patchi->geom_ivals_dynamic,ntimes_local*sizeof(int *));
    NewMemory((void **)&patchi->geom_vals,nvals*sizeof(float));
    NewMemory((void **)&patchi->geom_ivals,nvals*sizeof(char));
  }
  FORTgetembeddata(file, &endian_smv, &ntimes_local, &nvals, patchi->geom_times, 
    patchi->geom_nstatics, patchi->geom_ndynamics, patchi->geom_vals, &error, lenfile);

  init_histogram(patchi->histogram);

  update_histogram(patchi->geom_vals,nvals,patchi->histogram);

  complete_histogram(patchi->histogram);

  patchi->geom_ntimes=ntimes_local;
  patchi->geom_nvals=nvals;
  patchi->geom_ivals_static[0] = patchi->geom_ivals;
  patchi->geom_ivals_dynamic[0] = patchi->geom_ivals_static[0]+patchi->geom_nstatics[0];
  for(i=1;i<ntimes_local;i++){
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

  for(i=0;i<ngeominfoptrs;i++){
    geomdata *geomi;
    geomlistdata *geomlisti;
    int ntris;
    int j;
    float *color;

    geomi = geominfoptrs[i];
    if(geomi->display==0||geomi->loaded==0)continue;
    geomlisti = geomi->geomlistinfo;
    ntris = geomlisti->ntriangles;

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

/* ------------------ compare_transparent_triangles ------------------------ */

int compare_transparent_triangles( const void *arg1, const void *arg2 ){
  triangle *tri, *trj;

  tri = *(triangle **)arg1; 
  trj = *(triangle **)arg2;

  if(tri->distance<trj->distance)return 1;
  if(tri->distance>trj->distance)return -1;
  return 0;
}

/* ------------------ GetGeomInfoPtrs ------------------------ */

void GetGeomInfoPtrs(geomdata ***geominfoptrs_local,int *ngeominfoptrs_local){
  geomdata **gptr;
  int i,count=0;

  count=0;
  for(i=0;i<ngeominfo;i++){
    geomdata *geomi;

    geomi = geominfo + i;
    if(geomi->loaded==0||geomi->display==0)continue;
    count++;
  }
  for(i=0;i<nisoinfo;i++){
    iso *isoi;
    geomdata *geomi;
    
    isoi = isoinfo + i;
    if(isoi->loaded==0||isoi->display==0)continue;
    geomi = isoi->geominfo;
    if(geomi->loaded==0||geomi->display==0)continue;
    count++;
  }
  *ngeominfoptrs_local=count;
  if(count>0){
    NewMemory((void **)&gptr,count*sizeof(geomdata *));
  }
  else{
    *geominfoptrs_local=NULL;
    return;
  }
  *geominfoptrs_local=gptr;
  for(i=0;i<ngeominfo;i++){
    geomdata *geomi;

    geomi = geominfo + i;
    if(geomi->loaded==0||geomi->display==0)continue;
    *gptr++=geomi;
  }
  for(i=0;i<nisoinfo;i++){
    iso *isoi;
    geomdata *geomi;
    
    isoi = isoinfo + i;
    if(isoi->loaded==0||isoi->display==0)continue;
    geomi = isoi->geominfo;
    if(geomi->loaded==0||geomi->display==0)continue;
    *gptr++=geomi;
  }
}

/* ------------------ sort_triangles ------------------------ */

void Sort_Embedded_Geometry(float *mm){
  int i;
  int count_transparent,count_all;
  int itime;
  int *showlevels;

  CheckMemory;
  count_transparent=0;
  count_all=0;
  ntransparent_triangles=count_transparent;
  nopaque_triangles=count_all-count_transparent;
  showlevels=loaded_isomesh->showlevels;

  for(i=0;i<ngeominfoptrs;i++){
    geomlistdata *geomlisti;
    int j;
    geomdata *geomi;

    geomi = geominfoptrs[i];
    for(itime=0;itime<2;itime++){
      if(itime==0){
        geomlisti = geomi->geomlistinfo-1;
      }
      else{
        geomlisti = geomi->geomlistinfo+geomi->iframe;
      }

      count_all+=geomlisti->ntriangles;
      if(use_transparency_data==0)continue;
      for(j=0;j<geomlisti->ntriangles;j++){
        triangle *tri;
        float xyz[3];
        float *xyz1, *xyz2, *xyz3;
        float xyzeye[3];
        int isurf;

        tri = geomlisti->triangles + j;
        if(hilight_skinny==1&&tri->skinny==1)continue;
        if(tri->surf->transparent_level>=1.0)continue;
        isurf=tri->surf-surfinfo-nsurfinfo-1;
        if(showlevels[isurf]==0||tri->surf->transparent_level<=0.0){
          count_all--;
          continue;
        }
        count_transparent++;
        if(sort_embedded_geometry==0)continue;
        xyz1 = tri->points[0]->xyz;
        xyz2 = tri->points[1]->xyz;
        xyz3 = tri->points[2]->xyz;
        xyz[0] = ((xyz1[0]+xyz2[0]+xyz3[0])/3.0-xbar0)/xyzmaxdiff;
        xyz[1] = ((xyz1[1]+xyz2[1]+xyz3[1])/3.0-ybar0)/xyzmaxdiff;
        xyz[2] = ((xyz1[2]+xyz2[2]+xyz3[2])/3.0-zbar0)/xyzmaxdiff;

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
  for(i=0;i<ngeominfoptrs;i++){
    geomlistdata *geomlisti;
    int j;
    geomdata *geomi;

    geomi = geominfoptrs[i];
    for(itime=0;itime<2;itime++){
      if(itime==0){
        geomlisti = geomi->geomlistinfo-1;
      }
      else{
        geomlisti = geomi->geomlistinfo+geomi->iframe;
      }

      for(j=0;j<geomlisti->ntriangles;j++){
        triangle *tri;
        int isurf;

        tri = geomlisti->triangles + j;

        isurf=tri->surf-surfinfo-nsurfinfo-1;
        if(showlevels[isurf]==0){
          continue;
        }
        if(use_transparency_data==0||(hilight_skinny==1&&tri->skinny==1)||tri->surf->transparent_level>=1.0){
          opaque_triangles[count_all++]=tri;
        }
        else{
          transparent_triangles[count_transparent++]=tri;
        }
      }
    }
  }
  if(sort_embedded_geometry==1&&ntransparent_triangles>0){
    qsort((isotri **)transparent_triangles,(size_t)ntransparent_triangles,sizeof(triangle **),compare_transparent_triangles);
  }
}

/* ------------------ init_geom ------------------------ */

void init_geom(geomdata *geomi){
  geomi->display=0;
  geomi->loaded=0;
  geomi->geomlistinfo_0=NULL;
  geomi->geomlistinfo=NULL;
  geomi->times=NULL;
  geomi->ntimes=0;
  geomi->times=NULL;
  geomi->timeslist=NULL;
  geomi->float_vals=NULL;
  geomi->int_vals=NULL;
  geomi->nfloat_vals=0;
  geomi->nint_vals=0;
}

