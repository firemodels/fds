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
  angle_local = acos(arg)*RAD2DEG;
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
      PRINTF("Face/Vertex Summary\n");
      PRINTF("      Faces: %i\n",geomlisti->ntriangles);
      PRINTF(" slim faces: %i\n",nskinny);
      PRINTF("   Vertices: %i\n",geomlisti->npoints);
      PRINTF("     unused: %i\n",geomlisti->npoints-nused);
      PRINTF(" duplicates: %i\n\n",ndups);
      FREEMEMORY(points);
    }
  }
}

/* ------------------ draw_geom ------------------------ */

void draw_geom(int flag, int geomtype){
  int i;
  float black[]={0.0,0.0,0.0,1.0};
  float blue[]={0.0,0.0,1.0,1.0};
  float skinny_color[]={1.0,0.0,0.0,1.0};
  float *last_color=NULL;
  float last_transparent_level=-1.0;
  int ntris;
  triangle **tris;

  if(flag==DRAW_OPAQUE){
    ntris=nopaque_triangles;
    tris=opaque_triangles;
  }
  if(flag==DRAW_TRANSPARENT){
    ntris=ntransparent_triangles;
    tris=transparent_triangles;
  }

  if(ntris>0&&patchembedded==0&&showtrisurface==1&&geomtype==GEOM_STATIC){
    float *color;

    if(flag==DRAW_TRANSPARENT&&use_transparency_data==1)transparenton();

  // draw geometry surface

    if(flag==DRAW_TRANSPARENT&&use_transparency_data==1)transparenton();
    if(cullfaces==1)glDisable(GL_CULL_FACE);
    glEnable(GL_NORMALIZE);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,iso_specular);
    glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,iso_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glEnable(GL_COLOR_MATERIAL);
    
    glPushMatrix();
    glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
    glTranslatef(-xbar0,-ybar0,-zbar0);
    glBegin(GL_TRIANGLES);
    for(i=0;i<ntris;i++){
      triangle *trianglei;
      float transparent_level_local;
      texturedata *ti;
      int  j;

      trianglei = tris[i];
      ti = trianglei->textureinfo;
      if(visGeomTextures==1&&ti!=NULL&&ti->loaded==1)continue;
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
        glNormal3fv(trianglei->tri_norm);
        for(j=0;j<3;j++){
          point *pointj;

          pointj = trianglei->points[j];
          glVertex3fv(pointj->xyz);
        }
      }
      else{
        for(j=0;j<3;j++){
          point *pointj;

          pointj = trianglei->points[j];
          glNormal3fv(pointj->point_norm);
          glVertex3fv(pointj->xyz);
        }
      }
    }
    glEnd();

    if(visGeomTextures==1){
      texturedata *lasttexture;

      glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
      glEnable(GL_TEXTURE_2D);

      lasttexture=NULL;
      glBegin(GL_TRIANGLES);
      for(i=0;i<ntris;i++){
        triangle *trianglei;
        texturedata *texti;
        int j;

        trianglei = tris[i];
        texti = trianglei->textureinfo;
        if(texti==NULL||texti->loaded!=1)continue;
        if(lasttexture!=texti){
          glEnd();
          glBindTexture(GL_TEXTURE_2D,texti->name);
          glBegin(GL_TRIANGLES);
          lasttexture=texti;
        }
        for(j=0;j<3;j++){
          point *pointj;
          float *tpointj;

          pointj = trianglei->points[j];
          tpointj = trianglei->tpoints+2*j;
          glNormal3fv(pointj->point_norm);
          glTexCoord2fv(tpointj);
          glVertex3fv(pointj->xyz);
        }
      }
      glEnd();
      glDisable(GL_TEXTURE_2D);
    }

    glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
    glPopMatrix();
    if(flag==DRAW_TRANSPARENT){
      if(use_transparency_data==1)transparentoff();
      return;
    }
    if(cullfaces==1)glEnable(GL_CULL_FACE);
  }

#define VECFACTOR 0.03

  for(i=0;i<ngeominfoptrs;i++){
    geomdata *geomi;
    geomlistdata *geomlisti;
    int npoints;
    int nvolus;
    int j;
    float *color;

    geomi = geominfoptrs[i];
    if(geomi->loaded==0||geomi->display==0)continue;
    if(geomtype==GEOM_STATIC){
      geomlisti = geomi->geomlistinfo-1;
    }
    else{
      geomlisti = geomi->geomlistinfo+geomi->itime;
    }
    ntris = geomlisti->ntriangles;
    npoints = geomlisti->npoints;
    nvolus = geomlisti->nvolus;

  // draw volume outline
  
    if(nvolus>0){
      glPushMatrix();
      glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
      glTranslatef(-xbar0,-ybar0,-zbar0);
      glBegin(GL_LINES);
      for(j=0;j<nvolus;j++){
        tetrahedron *volumei;
        float *xyzptr[4];
        int *exterior,*duplicate;
        //             0
        //               \
        //           /   .3
        //             .  3
        //         / .  ./  
        //         1---2
        //
        int facelist[12]={0,1,2, 0,2,3, 0,3,1, 1,3,2};
        int k;

        volumei = geomlisti->volumes+j;
        exterior = volumei->exterior;
        duplicate = volumei->duplicate;
        xyzptr[0] = volumei->points[0]->xyz;
        xyzptr[1] = volumei->points[1]->xyz;
        xyzptr[2] = volumei->points[2]->xyz;
        xyzptr[3] = volumei->points[3]->xyz;

        color = volumei->surf->color;
        glColor3fv(color);

        for(k=0;k<4;k++){
          if(exterior[k]==1&&show_geometry_exterior==1||
             exterior[k]==0&&show_geometry_interior==1||
             duplicate[k]==1&&show_geometry_duplicates==1){
               glVertex3fv(xyzptr[facelist[3*k]]);
               glVertex3fv(xyzptr[facelist[3*k+1]]);
               glVertex3fv(xyzptr[facelist[3*k+1]]);
               glVertex3fv(xyzptr[facelist[3*k+2]]);
               glVertex3fv(xyzptr[facelist[3*k+2]]);
               glVertex3fv(xyzptr[facelist[3*k+0]]);
          }
        }
      }

      glEnd();
      glPopMatrix();
    }

    // draw geometry (faces) outline

    if(ntris>0&&showtrioutline==1){
      glPushMatrix();
      glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
      glTranslatef(-xbar0,-ybar0,-zbar0);
      glBegin(GL_LINES);
      for(j=0;j<ntris;j++){
        float *xyzptr[3];
        float *xyznorm;
        triangle *trianglei;

        trianglei = geomlisti->triangles+j;
       
        xyznorm=trianglei->tri_norm;
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

    // draw geometry points

    if(showtripoints==1&&geomlisti->npoints>0){
      glPushMatrix();
      glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
      glTranslatef(-xbar0,-ybar0,-zbar0);
      glPointSize(6.0);
      glBegin(GL_POINTS);
      for(j=0;j<geomlisti->npoints;j++){
        point *pointi;

        pointi = geomlisti->points+j;
        color = pointi->triangles[0]->surf->color;
        glColor3fv(color);
        glVertex3fv(pointi->xyz);
      }
      glEnd();
      glPopMatrix();
    }

    // draw geometry normal vectors

    if(showtrinormal==1){
      if(ntris>0&&smoothtrinormal==0){  // draw faceted normals
        glPushMatrix();
        glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
        glTranslatef(-xbar0,-ybar0,-zbar0);
        glBegin(GL_LINES);
        for(j=0;j<ntris;j++){
          float *p1, *p2, *p3;
          float *xyznorm;
          triangle *trianglei;
          float xyz1[3], xyz2[3];

          trianglei = geomlisti->triangles+j;
       
          xyznorm=trianglei->tri_norm;

          p1 = trianglei->points[0]->xyz;
          p2 = trianglei->points[1]->xyz;
          p3 = trianglei->points[2]->xyz;

          xyz1[0] = (p1[0] + p2[0] + p3[0])/3.0;
          xyz1[1] = (p1[1] + p2[1] + p3[1])/3.0;
          xyz1[2] = (p1[2] + p2[2] + p3[2])/3.0;
          xyz2[0] = xyz1[0] + SCALE2FDS(VECFACTOR)*xyznorm[0];
          xyz2[1] = xyz1[1] + SCALE2FDS(VECFACTOR)*xyznorm[1];
          xyz2[2] = xyz1[2] + SCALE2FDS(VECFACTOR)*xyznorm[2];

          glColor3fv(blue);
          glVertex3fv(xyz1);
          glVertex3fv(xyz2);
        }
        glEnd();
        
        glPointSize(6.0);  // draw points at end of vector
        glBegin(GL_POINTS);
        for(j=0;j<ntris;j++){
          float *p1, *p2, *p3;
          float *xyznorm;
          triangle *trianglei;
          float xyz1[3], xyz2[3];

          trianglei = geomlisti->triangles+j;
       
          xyznorm=trianglei->tri_norm;

          p1 = trianglei->points[0]->xyz;
          p2 = trianglei->points[1]->xyz;
          p3 = trianglei->points[2]->xyz;

          xyz1[0] = (p1[0] + p2[0] + p3[0])/3.0;
          xyz1[1] = (p1[1] + p2[1] + p3[1])/3.0;
          xyz1[2] = (p1[2] + p2[2] + p3[2])/3.0;
          xyz2[0] = xyz1[0] + SCALE2FDS(VECFACTOR)*xyznorm[0];
          xyz2[1] = xyz1[1] + SCALE2FDS(VECFACTOR)*xyznorm[1];
          xyz2[2] = xyz1[2] + SCALE2FDS(VECFACTOR)*xyznorm[2];

          glColor3fv(blue);
          glVertex3fv(xyz2);
        }
        glEnd();
        glPopMatrix();
      }
      if(npoints>0&&smoothtrinormal==1){ // draw smooth normals
        glPushMatrix();
        glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
        glTranslatef(-xbar0,-ybar0,-zbar0);
        glBegin(GL_LINES);
        for(j=0;j<npoints;j++){
          float *xyznorm;
          point *pointi;
          float *xyz1, xyz2[3];

          pointi = geomlisti->points+j;
          xyznorm = pointi->point_norm;       
          xyz1 = pointi->xyz;

          xyz2[0] = xyz1[0] + SCALE2FDS(VECFACTOR)*xyznorm[0];
          xyz2[1] = xyz1[1] + SCALE2FDS(VECFACTOR)*xyznorm[1];
          xyz2[2] = xyz1[2] + SCALE2FDS(VECFACTOR)*xyznorm[2];

          color = black;
          glColor3fv(color);
          glVertex3fv(xyz1);
          glVertex3fv(xyz2);
        }
        glEnd();

        glPointSize(6.0);  // draw points at end of vector
        glBegin(GL_POINTS);
        for(j=0;j<npoints;j++){
          float *xyznorm;
          point *pointi;
          float *xyz1, xyz2[3];

          pointi = geomlisti->points+j;
          xyznorm = pointi->point_norm;       
          xyz1 = pointi->xyz;

          xyz2[0] = xyz1[0] + SCALE2FDS(VECFACTOR)*xyznorm[0];
          xyz2[1] = xyz1[1] + SCALE2FDS(VECFACTOR)*xyznorm[1];
          xyz2[2] = xyz1[2] + SCALE2FDS(VECFACTOR)*xyznorm[2];

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
        xyznorm = trianglei->tri_norm;
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
        norm=pointi->point_norm;
        norm[0]=0.0;
        norm[1]=0.0;
        norm[2]=0.0;
        for(k=0;k<pointi->ntriangles;k++){
          float *norm2;
          triangle *trianglei;

          trianglei = pointi->triangles[k];
          norm2 = trianglei->tri_norm;
          norm[0]+=norm2[0];
          norm[1]+=norm2[1];
          norm[2]+=norm2[2];
        }
        ReduceToUnit(norm);
      }
    }  
  }
#ifdef pp_GEOMTEST
  box_bounds2[0]=DENORMALIZE_XX(0.25);
  box_bounds2[1]=DENORMALIZE_XX(0.75);
  box_bounds2[2]=DENORMALIZE_YY(0.25);
  box_bounds2[3]=DENORMALIZE_YY(0.75);
  box_bounds2[4]=DENORMALIZE_ZZ(0.25);
  box_bounds2[5]=DENORMALIZE_ZZ(0.75);
  box_translate[0]=0.0;
  box_translate[1]=0.0;
  box_translate[2]=0.0;

  tetra_vertices[0]=DENORMALIZE_XX(0.2);
  tetra_vertices[1]=DENORMALIZE_YY(0.2);
  tetra_vertices[2]=DENORMALIZE_ZZ(0.2);

  tetra_vertices[3]=DENORMALIZE_XX(0.8);
  tetra_vertices[4]=DENORMALIZE_YY(0.2);
  tetra_vertices[5]=DENORMALIZE_ZZ(0.2);

  tetra_vertices[6]=DENORMALIZE_XX(0.5);
  tetra_vertices[7]=DENORMALIZE_YY(0.8);
  tetra_vertices[8]=DENORMALIZE_ZZ(0.2);

  tetra_vertices[9]=DENORMALIZE_XX(0.5);
  tetra_vertices[10]=DENORMALIZE_YY(0.5);
  tetra_vertices[11]=DENORMALIZE_ZZ(0.8);

#endif

}

#define FORTREAD(var,count,STREAM) FSEEK(STREAM,4,SEEK_CUR);\
                           returncode=fread(var,4,count,STREAM);\
                           if(returncode!=count)returncode=0;\
                           if(endianswitch==1&&returncode!=0)endian_switch(var,count);\
                           FSEEK(STREAM,4,SEEK_CUR)

#define FORTREADBR(var,count,STREAM) FORTREAD(var,(count),STREAM);if(returncode==0)break;

/* ------------------ read_geom_header0 ------------------------ */

void read_geom_header0(geomdata *geomi, int *ntimes_local){
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
  int nverts=0, ntris=0;
  int first_frame_static;

  stream = fopen(geomi->file,"rb");
  if(stream==NULL){
    *ntimes_local=-1;
    return;
  }
  FSEEK(stream,4,SEEK_CUR);fread(&one,4,1,stream);FSEEK(stream,4,SEEK_CUR);
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

  FORTREAD(nvertfaces,2,stream);
  nverts=nvertfaces[0];
  ntris=nvertfaces[1];

  // static vertices

  if(nverts>0){
    FSEEK(stream,4+3*nverts*4+4,SEEK_CUR);
  }

  // static triangles

  if(ntris>0){
    FSEEK(stream,4+3*ntris*4+4,SEEK_CUR);
    FSEEK(stream,4+ntris*4+4,SEEK_CUR);
  }

  nt=0;
  for(;;){
    FORTREADBR(times_local,2,stream);
    FORTREADBR(nvertfaces,2,stream);
    nverts=nvertfaces[0];
    ntris=nvertfaces[1];

    // dynamic vertices

    if(nverts>0){
      FSEEK(stream,4+3*nverts*4+4,SEEK_CUR);
    }

    // dynamic faces

    if(ntris>0){
      FSEEK(stream,4+3*ntris*4+4,SEEK_CUR);
      FSEEK(stream,4+ntris*4+4,SEEK_CUR);
    }

    nt++;
  }
  *ntimes_local=nt;
  fclose(stream);
}

/* ------------------ read_geom_header2 ------------------------ */

void read_geom_header2(geomdata *geomi, int *ntimes_local){
  FILE *stream;
  int one=0,endianswitch=0;
  int nvertfacesvolus[3];
  float times_local[2];
  int nheaders[3];
  int nt;
  int returncode;
  int version;
  int nfloat_vals, nint_vals;
  int *int_vals;
  float *float_vals;
  int nverts=0, ntris=0, nvolus=0;
  int first_frame_static;
  int haves[3], have_surf=0, have_matl=0, have_texture=0;
  int header[3];
  float time_local;

  stream = fopen(geomi->file,"rb");
  if(stream==NULL){
    *ntimes_local=-1;
    return;
  }
  FSEEK(stream,4,SEEK_CUR);fread(&one,4,1,stream);FSEEK(stream,4,SEEK_CUR);
  if(one!=1)endianswitch=1;
  FORTREAD(&version,1,stream);

  FORTREAD(header,3,stream);
  first_frame_static=header[2];

  nt=0;
  if(first_frame_static==1)nt=-1;
  for(;;){
    FORTREADBR(&time_local,1,stream);
    FORTREADBR(nvertfacesvolus,3,stream);
    nverts=nvertfacesvolus[0];
    ntris=nvertfacesvolus[1];
    nvolus=nvertfacesvolus[2];

    FORTREADBR(haves,3,stream);
    have_surf=haves[0];
    have_matl=haves[1];
    have_texture=haves[2];

    // vertices

    if(nverts>0){
      FSEEK(stream,4+3*nverts*4+4,SEEK_CUR);
    }

    // faces

    if(ntris>0){
      FSEEK(stream,4+3*ntris*4+4,SEEK_CUR); 
      if(have_surf==1)FSEEK(stream,4+ntris*4+4,SEEK_CUR);
      if(have_texture==1)FSEEK(stream,4+6*ntris*4+4,SEEK_CUR);
    }

    // volumes

    if(nvolus>0){
      FSEEK(stream,4+4*nvolus*4+4,SEEK_CUR);
      if(have_matl==1)FSEEK(stream,4+nvolus*4+4,SEEK_CUR);
    }
    nt++;
  }
  *ntimes_local=nt;
  fclose(stream);
}

/* ------------------ read_geom_header ------------------------ */

void read_geom_header(geomdata *geomi, int *ntimes_local){
  FILE *stream;
  int version;
  int returncode;
  int one=0,endianswitch=0;

  stream = fopen(geomi->file,"rb");
  if(stream==NULL){
    *ntimes_local=-1;
    return;
  }
  FSEEK(stream,4,SEEK_CUR);fread(&one,4,1,stream);FSEEK(stream,4,SEEK_CUR);
  if(one!=1)endianswitch=1;
  FORTREAD(&version,1,stream);
  fclose(stream);

  if(version<=1){
    read_geom_header0(geomi, ntimes_local);
  }
  else{
    read_geom_header2(geomi, ntimes_local);
  }
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
  FSEEK(stream,4,SEEK_CUR);fread(&one,4,1,stream);FSEEK(stream,4,SEEK_CUR);
  if(one!=1)endianswitch=1;
  nt=-1;
  nv=0;
  for(;;){
    FORTREADBR(&time_local,1,stream);
    FORTREADBR(&nface_static,1,stream);
    if(nface_static!=0)FSEEK(stream,4+nface_static*4+4,SEEK_CUR);    
    FORTREADBR(&nface_dynamic,1,stream);
    if(nface_dynamic!=0)FSEEK(stream,4+nface_dynamic*4+4,SEEK_CUR);    
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
    read_geom(geomi,LOAD,GEOM_NORMAL,&errorcode);
  }
}

/* ------------------ read_geom0 ------------------------ */

void read_geom0(geomdata *geomi, int load_flag, int type, int *errorcode){
  FILE *stream;
  int one=1, endianswitch=0;
  int returncode;
  int ntimes_local;
  int i;
  point *points;
  triangle *triangles;
  int version;
  int nvertfacesvolus[3];
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
  geomi->geomlistinfo=NULL;
  FREEMEMORY(geomi->float_vals);
  FREEMEMORY(geomi->int_vals);
  geomi->nfloat_vals=0;
  geomi->nint_vals=0;

  if(load_flag==UNLOAD){
    geomi->loaded=0;
    geomi->display=0;
    return;
  }

  read_geom_header(geomi,&ntimes_local);
  if(ntimes_local<0)return;
  stream = fopen(geomi->file,"rb");
  if(stream==NULL)return;

  FSEEK(stream,4,SEEK_CUR);fread(&one,4,1,stream);FSEEK(stream,4,SEEK_CUR);
  if(one!=1)endianswitch=1;
  FORTREAD(&version,1,stream);

  FORTREAD(&nfloat_vals,1,stream);
  if(nfloat_vals>0)FSEEK(stream,4+nfloat_vals*4+4,SEEK_CUR);
  FORTREAD(&nint_vals,1,stream);
  if(nint_vals>0)FSEEK(stream,4+nint_vals*4+4,SEEK_CUR);

  geomi->ntimes=ntimes_local;
  geomi->itime=0;
  NewMemory((void **)&geomi->geomlistinfo_0,(ntimes_local+1)*sizeof(geomlistdata));
  geomi->geomlistinfo=geomi->geomlistinfo_0+1;
  NewMemory((void **)&geomi->times,ntimes_local*sizeof(float));

  for(i=-1;i<ntimes_local;i++){
    float times_local[2];
    geomlistdata *geomlisti;
    int nverts, ntris, nvolus;

    geomlisti = geomi->geomlistinfo+i;
    geomlisti->points=NULL;
    geomlisti->triangles=NULL;
    geomlisti->volumes=NULL;
    geomlisti->npoints=0;
    geomlisti->ntriangles=0;
    geomlisti->nvolus=0;
    if(i>=0){
      FORTREADBR(times_local,2,stream);
      geomi->times[i]=times_local[0];
    }
    FORTREADBR(nvertfacesvolus,2,stream);
    nvolus=0;
    nverts=nvertfacesvolus[0];
    ntris=nvertfacesvolus[1];
    if(i>=0){
      PRINTF("time=%.2f triangles: %i\n",times_local[0],ntris);
    }
    if(nverts>0){
      int ii;
      float *xyz=NULL;

      if(i<0)PRINTF("static geometry\n");
      NewMemory((void **)&xyz,3*nverts*sizeof(float));
      NewMemory((void **)&points,nverts*sizeof(point));
      geomlisti->points=points;
      geomlisti->npoints=nverts;
      FORTREADBR(xyz,3*nverts,stream);
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
      int offset=0;

      NewMemory((void **)&triangles,ntris*sizeof(triangle));
      NewMemory((void **)&ijk,3*ntris*sizeof(int));
      NewMemory((void **)&surf_ind,ntris*sizeof(int));
      geomlisti->triangles=triangles;
      geomlisti->ntriangles=ntris;
      FORTREADBR(ijk,3*ntris,stream);
      FORTREADBR(surf_ind,ntris,stream);
      if(type==GEOM_ISO)offset=nsurfinfo;
      for(ii=0;ii<ntris;ii++){
        surfdata *surfi;

        triangles[ii].points[0]=points+ijk[3*ii]-1;
        triangles[ii].points[1]=points+ijk[3*ii+1]-1;
        triangles[ii].points[2]=points+ijk[3*ii+2]-1;
        surfi=surfinfo + surf_ind[ii]+offset;
        triangles[ii].surf=surfi;
        triangles[ii].textureinfo=NULL;
      }
      FREEMEMORY(ijk);
      FREEMEMORY(surf_ind);
    }
  }
  geomi->loaded=1;
  geomi->display=1;
}

/* ------------------ read_geom2 ------------------------ */

void read_geom2(geomdata *geomi, int load_flag, int type, int *errorcode){
  FILE *stream;
  int one=1, endianswitch=0;
  int returncode;
  int ntimes_local;
  int i;
  point *points;
  triangle *triangles;
  tetrahedron *volumes;
  int version;
  int nvertfacesvolus[3];
  int nheaders[3], nfloat_vals, nint_vals, first_frame_static;
  int has_vals[3], has_surf, has_matl, has_texture;

  if(geomi->geomlistinfo!=NULL){
    for(i=-1;i<geomi->ntimes;i++){
      geomlistdata *geomlisti;

      geomlisti = geomi->geomlistinfo+i;
      FREEMEMORY(geomlisti->points);
      FREEMEMORY(geomlisti->triangles);
      FREEMEMORY(geomlisti->volumes);
    }  
  }
  FREEMEMORY(geomi->times);
  FREEMEMORY(geomi->geomlistinfo_0);
  geomi->geomlistinfo=NULL;
  FREEMEMORY(geomi->float_vals);
  FREEMEMORY(geomi->int_vals);
  geomi->nfloat_vals=0;
  geomi->nint_vals=0;

  if(load_flag==UNLOAD){
    geomi->loaded=0;
    geomi->display=0;
    return;
  }

  read_geom_header(geomi,&ntimes_local);
  if(ntimes_local<0)return;
  stream = fopen(geomi->file,"rb");
  if(stream==NULL)return;

  FSEEK(stream,4,SEEK_CUR);fread(&one,4,1,stream);FSEEK(stream,4,SEEK_CUR);
  if(one!=1)endianswitch=1;
  FORTREAD(&version,1,stream);

  FORTREAD(nheaders,3,stream);
  nfloat_vals=nheaders[0];
  nint_vals=nheaders[1];
  first_frame_static=nheaders[2];

  if(nfloat_vals>0)FSEEK(stream,4+nfloat_vals*4+4,SEEK_CUR);
  if(nint_vals>0)FSEEK(stream,4+nint_vals*4+4,SEEK_CUR);

  geomi->ntimes=ntimes_local;
  geomi->itime=0;
  NewMemory((void **)&geomi->geomlistinfo_0,(ntimes_local+1)*sizeof(geomlistdata));
  geomi->geomlistinfo=geomi->geomlistinfo_0+1;
  if(ntimes_local>0)NewMemory((void **)&geomi->times,ntimes_local*sizeof(float));

  for(i=-1;i<ntimes_local;i++){
    float time_local;
    geomlistdata *geomlisti;
    int nverts, ntris, nvolus;

    geomlisti = geomi->geomlistinfo+i;
    geomlisti->points=NULL;
    geomlisti->triangles=NULL;
    geomlisti->volumes=NULL;
    geomlisti->npoints=0;
    geomlisti->ntriangles=0;
    geomlisti->nvolus=0;

    FORTREADBR(&time_local,1,stream);
    if(i>=0)geomi->times[i]=time_local;

    FORTREADBR(nvertfacesvolus,3,stream);
    nverts=nvertfacesvolus[0];
    ntris=nvertfacesvolus[1];
    nvolus=nvertfacesvolus[2];

    FORTREADBR(has_vals,3,stream);
    has_surf=has_vals[0];
    has_matl=has_vals[1];
    has_texture=has_vals[2];

    if(i>=0){
      PRINTF("time=%.2f triangles: %i\n",time_local,ntris);
    }
    if(nverts>0){
      int ii;
      float *xyz=NULL;

      if(i<0)PRINTF("static geometry\n");
      NewMemory((void **)&xyz,3*nverts*sizeof(float));
      NewMemory((void **)&points,nverts*sizeof(point));
      geomlisti->points=points;
      geomlisti->npoints=nverts;
      FORTREADBR(xyz,3*nverts,stream);
      for(ii=0;ii<nverts;ii++){
        points[ii].xyz[0]=xyz[3*ii];
        points[ii].xyz[1]=xyz[3*ii+1];
        points[ii].xyz[2]=xyz[3*ii+2];
      }
      FREEMEMORY(xyz);
    }
    if(ntris>0){
      int *surf_ind=NULL,*ijk=NULL;
      float *texture_coords=NULL;
      int ii;
      int offset=0;

      NewMemory((void **)&triangles,ntris*sizeof(triangle));
      NewMemory((void **)&ijk,3*ntris*sizeof(int));
      NewMemory((void **)&surf_ind,ntris*sizeof(int));
      if(has_texture==1)NewMemory((void **)&texture_coords,6*ntris*sizeof(float));
      geomlisti->triangles=triangles;
      geomlisti->ntriangles=ntris;
      FORTREADBR(ijk,3*ntris,stream);
      if(has_surf==1){
        FORTREADBR(surf_ind,ntris,stream);
      }
      if(has_texture==1){
        FORTREADBR(texture_coords,6*ntris,stream);
      }
      if(type==GEOM_ISO)offset=nsurfinfo;
      for(ii=0;ii<ntris;ii++){
        surfdata *surfi;

        triangles[ii].points[0]=points+ijk[3*ii]-1;
        triangles[ii].points[1]=points+ijk[3*ii+1]-1;
        triangles[ii].points[2]=points+ijk[3*ii+2]-1;
        if(has_texture==1){
          triangles[ii].tpoints[0]=texture_coords[6*ii];
          triangles[ii].tpoints[1]=texture_coords[6*ii+1];
          triangles[ii].tpoints[2]=texture_coords[6*ii+2];
          triangles[ii].tpoints[3]=texture_coords[6*ii+3];
          triangles[ii].tpoints[4]=texture_coords[6*ii+4];
          triangles[ii].tpoints[5]=texture_coords[6*ii+5];
          CheckMemory;
        }
        if(has_surf==1){
          surfi=surfinfo + surf_ind[ii]+offset;
          triangles[ii].surf=surfi;
          triangles[ii].textureinfo=surfi->textureinfo;
        }
        else{
          triangles[ii].surf=surfacedefault;
          triangles[ii].textureinfo=NULL;
        }
      }

      FREEMEMORY(ijk);
      FREEMEMORY(surf_ind);
      if(has_texture==1)FREEMEMORY(texture_coords);
    }
    if(nvolus>0){ 
      int ii;
      int *ijk;
      int *surf_ind=NULL;
      int offset=0;

      NewMemory((void **)&volumes,nvolus*sizeof(tetrahedron));
      geomlisti->volumes=volumes;
      NewMemory((void **)&ijk,4*nvolus*sizeof(int));
      FORTREADBR(ijk,4*nvolus,stream);
      for(ii=0;ii<nvolus;ii++){
        volumes[ii].points[0]=points+ijk[4*ii]-1;
        volumes[ii].points[1]=points+ijk[4*ii+1]-1;
        volumes[ii].points[2]=points+ijk[4*ii+2]-1;
        volumes[ii].points[3]=points+ijk[4*ii+3]-1;
      }
      FREEMEMORY(ijk);
      NewMemory((void **)&surf_ind,nvolus*sizeof(int));
      if(has_matl==1){
        FORTREADBR(surf_ind,nvolus,stream);
      }
      if(type==GEOM_ISO)offset=nsurfinfo;
      for(ii=0;ii<nvolus;ii++){
        surfdata *surfi;

        surfi=surfinfo + offset;
        if(has_matl==1)surfi+=surf_ind[ii];
        volumes[ii].surf=surfi;
        volumes[ii].textureinfo=surfi->textureinfo;
      }
      FREEMEMORY(surf_ind);
      geomlisti->nvolus=nvolus;
    }
  }
  geomi->loaded=1;
  geomi->display=1;
}

/* ------------------ reorder_face ------------------------ */

void reorder_face(int *faces){
  int face_temp[5];

  face_temp[0]=faces[0];
  face_temp[1]=faces[1];
  face_temp[2]=faces[2];
  face_temp[3]=faces[0];
  face_temp[4]=faces[1];
  if(faces[0]<=MIN(faces[1],faces[2]))return;
  if(faces[1]<=MIN(faces[0],faces[2])){
    faces[0]=face_temp[1];
    faces[1]=face_temp[2];
    faces[2]=face_temp[3];
    return;
  }
  faces[0]=face_temp[2];
  faces[1]=face_temp[3];
  faces[2]=face_temp[4];
}

/* ------------------ find_common_face ------------------------ */

void find_common_face(tetrahedron *tetra1, tetrahedron *tetra2){
  int *verts1,*verts2, *face1, *face2;
  int i;
  int ncommon=0;

  verts1=tetra1->vert_index;
  verts2=tetra2->vert_index;
  face1=tetra1->faces;
  face2=tetra2->faces;
  for(i=0;i<4;i++){
    int j;

    for(j=0;j<4;j++){
      if(verts1[i]==verts2[j]){
        ncommon++;
        break;
      }
    }
  }
  if(ncommon<3)return; // two tetrahedrons need at least 3 common vertices to have a common face

  for(i=0;i<4;i++){
    int j;
    int *facei;

    facei = face1 + 3*i;

    for(j=0;j<4;j++){
      int *facej;

      facej = face2 + 3*j;
      if(facei[0]!=facej[0])continue;
      if(facei[1]==facej[1]&&facei[2]==facej[2]){ // duplicate face
        tetra1->duplicate[i]=1;
        tetra2->duplicate[j]=1;
      }
      if(facei[1]==facej[2]&&facei[2]==facej[1]){ // exterior face
        tetra1->exterior[i]=0;
        tetra2->exterior[j]=0;
      }

    }
  }

}

/* ------------------ classify_geom ------------------------ */

void classify_geom(geomdata *geomi){
  int ntimes,i;
  
  ntimes = geomi->ntimes;
  for(i=-1;i<ntimes;i++){
    float time_local;
    geomlistdata *geomlisti;
    int nverts, ntris, nvolus;
    int j;
    point *pointbase;

    geomlisti = geomi->geomlistinfo+i;
    nverts=geomlisti->npoints;
    nvolus=geomlisti->nvolus;
    if(nverts==0||nvolus==0||geomlisti->points==NULL)continue;
    pointbase = geomlisti->points;
    for(j=0;j<nvolus;j++){
      tetrahedron *tetrai;
      int *vert_index;
      point **points;
      int *faces;

      tetrai = geomlisti->volumes+j;
      vert_index = tetrai->vert_index;
      points = tetrai->points;
      faces = tetrai->faces;
      tetrai->exterior[0]=1;
      tetrai->exterior[1]=1;
      tetrai->exterior[2]=1;
      tetrai->exterior[3]=1;
      tetrai->duplicate[0]=0;
      tetrai->duplicate[1]=0;
      tetrai->duplicate[2]=0;
      tetrai->duplicate[3]=0;
      vert_index[0] = points[0]-pointbase;
      vert_index[1] = points[1]-pointbase;
      vert_index[2] = points[2]-pointbase;
      vert_index[3] = points[3]-pointbase;

      faces[0]=vert_index[0];
      faces[1]=vert_index[1];
      faces[2]=vert_index[2];
      reorder_face(faces);

      faces[3]=vert_index[0];
      faces[4]=vert_index[2];
      faces[5]=vert_index[3];
      reorder_face(faces+3);

      faces[6]=vert_index[0];
      faces[7]=vert_index[3];
      faces[8]=vert_index[1];
      reorder_face(faces+6);

      faces[9]=vert_index[1];
      faces[10]=vert_index[3];
      faces[11]=vert_index[2];
      reorder_face(faces+9);
    }
    for(j=0;j<nvolus;j++){
      int k;
      tetrahedron *tetraj;

      tetraj = geomlisti->volumes+j;

      for(k=0;k<nvolus;k++){
        tetrahedron *tetrak;

          if(j==k)continue;
          tetrak = geomlisti->volumes+k;
          find_common_face(tetraj,tetrak);
      }
    }
  }
}

/* ------------------ read_geom ------------------------ */

void read_geom(geomdata *geomi, int load_flag, int type, int *errorcode){
  FILE *stream;
  int version;
  int returncode;
  int one=0,endianswitch=0;

  stream = fopen(geomi->file,"rb");
  if(stream==NULL){
    return;
  }
  FSEEK(stream,4,SEEK_CUR);fread(&one,4,1,stream);FSEEK(stream,4,SEEK_CUR);
  if(one!=1)endianswitch=1;
  FORTREAD(&version,1,stream);
  fclose(stream);

  if(version<=1){
    read_geom0(geomi,load_flag,type,errorcode);
  }
  else{
    read_geom2(geomi,load_flag,type,errorcode);
  }
  if(load_flag==LOAD)classify_geom(geomi);
}

/* ------------------ read_geomdata ------------------------ */

void read_geomdata(int ifile, int load_flag, int *errorcode){
  patchdata *patchi;
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
  if(load_flag==UNLOAD){
    plotstate=getplotstate(DYNAMIC_PLOTS);
    update_patchtype();
    update_unit_defs();
    Update_Times();
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
    patchi->geom_nstatics, patchi->geom_ndynamics, patchi->geom_vals, &redirect, &error, lenfile);

  init_histogram(patchi->histogram);

  update_histogram(patchi->geom_vals,nvals,patchi->histogram);

  complete_histogram(patchi->histogram);

  patchi->ngeom_times=ntimes_local;
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
  Update_Times();
  Update_Framenumber(1);
}

/* ------------------ draw_geomdata ------------------------ */

void draw_geomdata(patchdata *patchi, int geomtype){
  int i;

  for(i=0;i<ngeominfoptrs;i++){
    geomdata *geomi;
    geomlistdata *geomlisti;
    int ntris;
    int j;
    float *color;

    geomi = geominfoptrs[i];
    if(geomi->display==0||geomi->loaded==0)continue;
    if(geomtype==GEOM_STATIC){
      geomlisti = geomi->geomlistinfo-1;
    }
    else{
      geomlisti = geomi->geomlistinfo+geomi->itime;
    }

    ntris = geomlisti->ntriangles;
    if(ntris==0)continue;

    glEnable(GL_NORMALIZE);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,iso_specular);
    glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,iso_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glEnable(GL_COLOR_MATERIAL);

    glPushMatrix();
    glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
    glTranslatef(-xbar0,-ybar0,-zbar0);
    glBegin(GL_TRIANGLES);
    if(smoothtrinormal==0){
      for(j=0;j<ntris;j++){
        float *xyzptr[3];
        float *xyznorm;
        triangle *trianglei;
        int color_index;

        trianglei = geomlisti->triangles+j;
       
        xyznorm=trianglei->tri_norm;
        glNormal3fv(xyznorm);

        //color_index = patchi->geom_ival_static[j];
        color_index = patchi->geom_ival_dynamic[j];
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
       
//        color_index = patchi->geom_ival_static[j];
        color_index = patchi->geom_ival_dynamic[j];
        color=rgb_patch+4*color_index;
        glColor3fv(color);

        xyznorm = trianglei->points[0]->point_norm;
        glNormal3fv(xyznorm);
        xyzptr[0] = trianglei->points[0]->xyz;
        glVertex3fv(xyzptr[0]);

        xyznorm = trianglei->points[1]->point_norm;
        glNormal3fv(xyznorm);
        xyzptr[1] = trianglei->points[1]->xyz;
        glVertex3fv(xyzptr[1]);

        xyznorm = trianglei->points[2]->point_norm;
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
    isodata *isoi;
    geomdata *geomi;
    
    isoi = isoinfo + i;
    if(isoi->loaded==0||isoi->display==0)continue;
    geomi = isoi->geominfo;
    if(geomi==NULL)continue;
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
    isodata *isoi;
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
  int *showlevels=NULL;

  CheckMemory;
  count_transparent=0;
  count_all=0;
  ntransparent_triangles=count_transparent;
  nopaque_triangles=count_all-count_transparent;
  if(loaded_isomesh!=NULL)showlevels=loaded_isomesh->showlevels;

  for(i=0;i<ngeominfoptrs;i++){
    geomlistdata *geomlisti;
    int j;
    geomdata *geomi;

    geomi = geominfoptrs[i];
    for(itime=0;itime<2;itime++){ //xxx was itime<2 check this
      if(itime==0){
        geomlisti = geomi->geomlistinfo-1;
      }
      else{
        geomlisti = geomi->geomlistinfo+geomi->itime;
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
        if((showlevels!=NULL&&showlevels[isurf]==0)||tri->surf->transparent_level<=0.0){
          count_all--;
          continue;
        }
        count_transparent++;
        if(sort_embedded_geometry==0)continue;
        xyz1 = tri->points[0]->xyz;
        xyz2 = tri->points[1]->xyz;
        xyz3 = tri->points[2]->xyz;
        xyz[0] = NORMALIZE_X((xyz1[0]+xyz2[0]+xyz3[0])/3.0);
        xyz[1] = NORMALIZE_Y((xyz1[1]+xyz2[1]+xyz3[1])/3.0);
        xyz[2] = NORMALIZE_Z((xyz1[2]+xyz2[2]+xyz3[2])/3.0);

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
        geomlisti = geomi->geomlistinfo+geomi->itime;
      }

      for(j=0;j<geomlisti->ntriangles;j++){
        triangle *tri;
        int isurf;

        tri = geomlisti->triangles + j;

        isurf=tri->surf-surfinfo-nsurfinfo-1;
        if(showlevels!=NULL&&showlevels[isurf]==0){
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
  geomi->surf=NULL;
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

