#include "options.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "smv_endian.h"
#include "update.h"
#include "smokeviewvars.h"

tetdata *volume_list;
tridata *triangle_list;
void Volume_CB(int var);

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
  vertdata *verti, *vertj;
  float *xyzi, *xyzj;

  verti = (vertdata *)arg1;
  vertj = (vertdata *)arg2;
  xyzi = verti->xyz;
  xyzj = vertj->xyz;

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
  float denom;


//         v1
//        /  \
//       d3   d2
//      /      \
//    v2---d1---v3

//       d2^2 + d3^2 - d1^2
// arg = ------------------
//           2*d2*d3

// d2==0.0 ==> d3==d1 ==> arg=0.0
// d3==0.0 ==> d2==d1 ==> arg=0.0


  denom = 2.0*d2*d3;
  if(ABS(denom) > 0.0){
    arg = CLAMP((d2*d2 + d3*d3 - d1*d1) / denom,-1.0,1.0);
    angle_local = acos(arg)*RAD2DEG;
  }
  else{
    angle_local = acos(0.0)*RAD2DEG;
  }
  return angle_local;
}

/* ------------------ get_minangle ------------------------ */

float get_minangle(tridata *trii){
  float minangle;
  float d1, d2, d3;
  float *xyz1, *xyz2, *xyz3;
  float angle1, angle2, angle3;

  xyz1 = trii->verts[0]->xyz;
  xyz2 = trii->verts[1]->xyz;
  xyz3 = trii->verts[2]->xyz;
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

/* ------------------ get_faceinfo ------------------------ */

void get_faceinfo(void){
  int i;

  for(i=0;i<ngeominfoptrs;i++){
    geomdata *geomi;
    geomlistdata *geomlisti;
    vertdata **verts;
    int j;
    int ndups=0,nused=0,nskinny=0;

    geomi = geominfoptrs[i];
    geomlisti = geomi->geomlistinfo;

    if(geomlisti->nverts>0){
      NewMemory((void **)&verts,geomlisti->nverts*sizeof(vertdata *));
      for(j=0;j<geomlisti->nverts;j++){
        verts[j]=geomlisti->verts+j;
        verts[j]->nused=0;
      }
      for(j=0;j<geomlisti->ntriangles;j++){
        tridata *trii;

        trii = geomlisti->triangles + j;
        trii->verts[0]->nused=0;
        trii->verts[1]->nused=0;
        trii->verts[2]->nused=0;
      }
      qsort(verts,geomlisti->nverts,sizeof(vertdata *),compare_verts);
      for(j=1;j<geomlisti->nverts;j++){
        if(compare_verts(verts[j-1],verts[j])==0)ndups++;
      }
      for(j=0;j<geomlisti->ntriangles;j++){
        tridata *trii;

        trii = geomlisti->triangles + j;
        trii->verts[0]->nused++;
        trii->verts[1]->nused++;
        trii->verts[2]->nused++;
        if(get_minangle(trii)<=10.0){
          trii->skinny=1;
          nskinny++;
        }
        else{
          trii->skinny=0;
        }
      }
      for(j=0;j<geomlisti->nverts;j++){
        if(verts[j]->nused>0)nused++;
      }
      PRINTF("Face/Vertex Summary\n");
      PRINTF("      Faces: %i\n",geomlisti->ntriangles);
      PRINTF(" slim faces: %i\n",nskinny);
      PRINTF("   Vertices: %i\n",geomlisti->nverts);
      PRINTF("     unused: %i\n",geomlisti->nverts-nused);
      PRINTF(" duplicates: %i\n\n",ndups);
      FREEMEMORY(verts);
    }
  }
}

/* ------------------ draw_geomdiag ------------------------ */

void draw_geomdiag(void){
  int i;

    glPushMatrix();
    glScalef(SCALE2SMV(1.0), SCALE2SMV(1.0), SCALE2SMV(1.0));
    glTranslatef(-xbar0, -ybar0, -zbar0);
    glBegin(GL_TRIANGLES);
    for(i = 0; i < ngeomdiaginfo; i++){
      geomdiagdata *geomdiagi;
      geomlistdata *geomframe;
      int ntriangles;
      int j;

      geomdiagi = geomdiaginfo + i;
      geomframe = geomdiagi->geom->geomlistinfo_0;
      ntriangles = geomframe->ntriangles;
      for(j = 0; j < ntriangles; j++){
        tridata *trianglej;

        trianglej = geomframe->triangles + j;
        glVertex3fv(trianglej->verts[0]->xyz);
        glVertex3fv(trianglej->verts[1]->xyz);
        glVertex3fv(trianglej->verts[2]->xyz);
      }
    }
    glEnd();
    glPopMatrix();
}

  /* ------------------ get_geom_zbounds ------------------------ */

void get_geom_zbounds(float *zmin, float *zmax){
  int j;
  int first = 1;
  
  for(j = 0; j < ngeominfoptrs; j++){
    geomdata *geomi;
    int iend, ii;

    geomi = geominfoptrs[j];
    if(geomi->loaded == 0 || geomi->display == 0)continue;
    if(geomi->geomtype != GEOM_GEOM&&geomi->geomtype != GEOM_ISO)continue;

    iend = geomi->ntimes;
    if(geomi->currentframe != NULL)iend = 1;

    for(ii = -1; ii < iend; ii++){
      geomlistdata *geomlisti;
      int k;

      if(ii == -1 || geomi->currentframe == NULL){
        geomlisti = geomi->geomlistinfo + ii;
      }
      else{
        geomlisti = geomi->currentframe;
      }
      for(k = 0; k < geomlisti->nverts; k++){
        float zval;
        vertdata *vertk;

        vertk = geomlisti->verts + k;
        zval = vertk->xyz[2];
        if(first == 1){
          *zmin = zval;
          *zmax = zval;
          first = 0;
        }
        else{
          *zmin = MIN(*zmin, zval);
          *zmax = MAX(*zmax, zval);
        }
      }
    }
  }
}

  /* ------------------ draw_geom ------------------------ */

void draw_geom(int flag, int timestate){
  int i;
  float black[]={0.0,0.0,0.0,1.0};
  float blue[]={0.0,0.0,1.0,1.0};
  float skinny_color[]={1.0,0.0,0.0,1.0};
  float *last_color=NULL;
  float last_transparent_level=-1.0;
  int ntris;
  tridata **tris;

  if(flag == DRAW_OPAQUE){
    ntris=nopaque_triangles;
    tris=opaque_triangles;
  }
  if(flag==DRAW_TRANSPARENT){
    ntris=ntransparent_triangles;
    tris=transparent_triangles;
  }

  if(ntris>0&&timestate==GEOM_STATIC){
    float *color;

  // draw geometry surface

    if(flag==DRAW_TRANSPARENT&&use_transparency_data==1)transparenton();
    if(cullfaces == 1)glDisable(GL_CULL_FACE);
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
      tridata *trianglei;
      float transparent_level_local;
      texturedata *ti;
      int  j;

      trianglei = tris[i];
      if(trianglei->exterior == 1 && show_faces_exterior == 0)continue;
      if(trianglei->exterior == 0 && show_faces_interior == 0)continue;
      if(trianglei->geomtype == GEOM_GEOM&&show_faces_solid == 0)continue;
      if(trianglei->geomtype == GEOM_ISO&&show_iso_solid == 0)continue;

      ti = trianglei->textureinfo;
      if(show_texture_1dimage==1)continue;
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
      if((trianglei->geomtype == GEOM_GEOM&&smooth_geom_normal == 0) ||
         (trianglei->geomtype == GEOM_ISO &&smooth_iso_normal == 0)){
        glNormal3fv(trianglei->tri_norm);
        for(j=0;j<3;j++){
          vertdata *vertj;

          vertj = trianglei->verts[j];
          glVertex3fv(vertj->xyz);
        }
      }
      else{
        for(j=0;j<3;j++){
          vertdata *vertj;

          vertj = trianglei->verts[j];
          glNormal3fv(trianglei->vert_norm+3*j);
          glVertex3fv(vertj->xyz);
        }
      }
    }
    glEnd();

    if(visGeomTextures == 1 || show_texture_1dimage == 1){
      texturedata *lasttexture;

      if(show_texture_1dimage == 1){
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
        glEnable(GL_TEXTURE_1D);
        glBindTexture(GL_TEXTURE_1D, terrain_colorbar_id);
      }
      else{
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
        glEnable(GL_TEXTURE_2D);
      }

      lasttexture=NULL;
      glBegin(GL_TRIANGLES);
      for(i=0;i<ntris;i++){
        tridata *trianglei;
        texturedata *texti;
        int j;

        trianglei = tris[i];
        if(trianglei->exterior == 1 && show_faces_exterior == 0)continue;
        if(trianglei->exterior == 0 && show_faces_interior == 0)continue;
        if(trianglei->geomtype == GEOM_ISO &&show_iso_outline == 0)continue;

        if(show_texture_1dimage == 1){
          for(j = 0; j < 3; j++){
            vertdata *vertj;
            float *xyz, texture_z;

            vertj = trianglei->verts[j];
            xyz = vertj->xyz;
            //       znew = terrain_zmin + geom_vert_exag*(zold-terrain_zmin);
            //    zmaxnew = terrain_zmin + geom_vert_exag*(terrain_zmax-terrain_zmin)
            //       zold = terrain_zmin + (znew-terrain_zmin)/geom_vert_exag

            texture_z = (xyz[2] - terrain_zmin)/(geom_vert_exag*(terrain_zmax-terrain_zmin));

            glNormal3fv(vertj->vert_norm);
            glTexCoord1f(texture_z);
            glVertex3fv(xyz);
          }
        }
        else{
          texti = trianglei->textureinfo;
          if(texti == NULL || texti->loaded != 1)continue;
          if(lasttexture != texti){
            glEnd();
            glBindTexture(GL_TEXTURE_2D, texti->name);
            glBegin(GL_TRIANGLES);
            lasttexture = texti;
          }
          for(j = 0; j < 3; j++){
            vertdata *vertj;
            float *tvertj;

            vertj = trianglei->verts[j];
            tvertj = trianglei->tverts + 2 * j;
            glNormal3fv(vertj->vert_norm);
            glTexCoord2fv(tvertj);
            glVertex3fv(vertj->xyz);
          }
        }
      }
      glEnd();
      if(show_texture_1dimage == 1){
        glDisable(GL_TEXTURE_1D);
      }
      else{
        glDisable(GL_TEXTURE_2D);
      }
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

  for(i=0;i<ngeominfoptrs;i++){
    geomdata *geomi;
    geomlistdata *geomlisti;
    int nvolumes;
    int j;
    float *color;

    geomi = geominfoptrs[i];
    if(geomi->loaded==0||geomi->display==0)continue;
    if(geomi->geomtype!=GEOM_GEOM&&geomi->geomtype!=GEOM_ISO)continue;
    if(timestate==GEOM_STATIC){
      geomlisti = geomi->geomlistinfo-1;
    }
    else{
      geomlisti = geomi->geomlistinfo+geomi->itime;
    }
    ntris = geomlisti->ntriangles;
    nvolumes = geomlisti->nvolumes;

    if(nvolumes > 0 && show_volumes_solid == 1){

      // draw volume solid

      last_color = NULL;
      glPushMatrix();
      glScalef(SCALE2SMV(1.0), SCALE2SMV(1.0), SCALE2SMV(1.0));
      glTranslatef(-xbar0, -ybar0, -zbar0);

      glEnable(GL_NORMALIZE);
      glShadeModel(GL_SMOOTH);
      glEnable(GL_LIGHTING);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, iso_specular);
      glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, iso_shininess);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, block_ambient2);
      glEnable(GL_COLOR_MATERIAL);

      last_color = NULL;
      glBegin(GL_TRIANGLES);
      for(j = 0; j < nvolumes; j++){
        tetdata *volumei;
        float *xyzptr[4];
        int *exterior;
        //
        //             0
        //            /  \
        //           /   .3
        //             .   \
        //         / .      \
        //         1--------2
        //
        int facelist[12] = {0, 1, 2, 0, 2, 3, 0, 3, 1, 1, 3, 2};
        int k;

        volumei = geomlisti->volumes + j;
        exterior = volumei->exterior;
        xyzptr[0] = volumei->verts[0]->xyz;
        xyzptr[1] = volumei->verts[1]->xyz;
        xyzptr[2] = volumei->verts[2]->xyz;
        xyzptr[3] = volumei->verts[3]->xyz;


        color = volumei->matl->color;
        if(last_color != color){
          glColor3fv(color);
          last_color = color;
        }

        for(k = 0; k < 4; k++){
          int kk;
          float *v0, *v1, *v2;
          float v1m0[3], v2m0[3], v2m1[3], vcross[3];
          float v0delta[3], v1delta[3], v2delta[3];

          if(show_volumes_exterior == 0 && exterior[k] == 1)continue;
          if(show_volumes_interior == 0 && exterior[k] == 0)continue;
          v0 = xyzptr[facelist[3 * k]];
          v1 = xyzptr[facelist[3 * k + 1]];
          v2 = xyzptr[facelist[3 * k + 2]];
          VECDIFF3(v1m0, v1, v0);
          VECDIFF3(v2m0, v2, v0);
          VECDIFF3(v2m1, v2, v1);
          CROSS(vcross, v1m0, v2m0);

          for(kk = 0; kk < 3; kk++){
            v0delta[kk] = v0[kk] + face_factor*v1m0[kk] + face_factor*v2m0[kk];
            v1delta[kk] = v1[kk] - face_factor*v1m0[kk] + face_factor*v2m1[kk];
            v2delta[kk] = v2[kk] - face_factor*v2m0[kk] - face_factor*v2m1[kk];
          }
          glNormal3fv(vcross);
          glVertex3fv(v0delta);
          glVertex3fv(v1delta);
          glVertex3fv(v2delta);
        }
      }
      glEnd();
      glDisable(GL_COLOR_MATERIAL);
	  glDisable(GL_LIGHTING);
	  glPopMatrix();
    }

      // draw volume outline

    if(nvolumes > 0 && show_volumes_outline == 1){
      last_color = NULL;
      glPushMatrix();
      glScalef(SCALE2SMV(1.0), SCALE2SMV(1.0), SCALE2SMV(1.0));
      glTranslatef(-xbar0, -ybar0, -zbar0);
      glDisable(GL_COLOR_MATERIAL);
      glDisable(GL_LIGHTING);
      glLineWidth(20.0);
      glBegin(GL_LINES);
      for(j=0;j<nvolumes;j++){
        tetdata *volumei;
        float *xyzptr[4];
        int *exterior;
        //
        //             0
        //            /  \
        //           /   .3
        //             .   \
        //         / .      \
        //         1--------2
        //
        int facelist[12]={0,1,2, 0,2,3, 0,3,1, 1,3,2};
        int k;

        volumei = geomlisti->volumes+j;
        exterior = volumei->exterior;
        xyzptr[0] = volumei->verts[0]->xyz;
        xyzptr[1] = volumei->verts[1]->xyz;
        xyzptr[2] = volumei->verts[2]->xyz;
        xyzptr[3] = volumei->verts[3]->xyz;

        for(k=0;k<4;k++){
          if(show_volumes_exterior == 0 && exterior[k] == 1)continue;
          if(show_volumes_interior == 0 && exterior[k] == 0)continue;
          if(show_volumes_solid==1){
             color=black;
          }
          else{
            color = volumei->matl->color;
          }
          if(last_color!=color){
            glColor3fv(color);
            last_color=color;
          }
          glVertex3fv(xyzptr[facelist[3*k]]);
          glVertex3fv(xyzptr[facelist[3*k+1]]);
          glVertex3fv(xyzptr[facelist[3*k+1]]);
          glVertex3fv(xyzptr[facelist[3*k+2]]);
          glVertex3fv(xyzptr[facelist[3*k+2]]);
          glVertex3fv(xyzptr[facelist[3*k+0]]);
        }
      }
      glEnd();

      glPopMatrix();
    }

    // draw geometry (faces) outline

    last_color=NULL;
    if(ntris>0){
      glPushMatrix();
      glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
      glTranslatef(-xbar0,-ybar0,-zbar0);
      glLineWidth(2.0);
      glBegin(GL_LINES);
      for(j=0;j<ntris;j++){
        float *xyzptr[3];
        float *xyznorm;
        tridata *trianglei;

        trianglei = geomlisti->triangles+j;
        if(trianglei->exterior == 1 && show_faces_exterior == 0)continue;
        if(trianglei->exterior == 0 && show_faces_interior == 0)continue;
        if(trianglei->geomtype == GEOM_GEOM&&show_faces_outline == 0)continue;
        if(trianglei->geomtype == GEOM_ISO&&show_iso_outline == 0)continue;

        xyznorm=trianglei->tri_norm;
        glNormal3fv(xyznorm);

        xyzptr[0] = trianglei->verts[0]->xyz;
        xyzptr[1] = trianglei->verts[1]->xyz;
        xyzptr[2] = trianglei->verts[2]->xyz;

        if(show_iso_solid==1){
          color = black;
        }
        else{
          color = trianglei->surf->color;
        }
        if(last_color!=color){
          glColor3fv(color);
          last_color=color;
        }
        {
          int ind[6] = {0, 1, 1, 2, 2, 0};
          int k;

          for(k = 0; k < 6; k++){
            float *xyzval, *pknorm, xyzval2[3];
            vertdata *pk;

            pk = trianglei->verts[ind[k]];
            pknorm = pk->vert_norm;
            xyzval = xyzptr[ind[k]];

            VECEQ3(xyzval2, pknorm);
            VEC3MA(xyzval2, geom_outline_offset);
            VECADD3(xyzval2, xyzval2, xyzval);
            glVertex3fv(xyzval2);
          }
        }
      }
      glEnd();
      glPopMatrix();
    }

    // draw geometry verts

    last_color=NULL;
    if(geomlisti->nverts>0){
      glPushMatrix();
      glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
      glTranslatef(-xbar0,-ybar0,-zbar0);
      glPointSize(6.0);
      glBegin(GL_POINTS);
      for(j=0;j<geomlisti->nverts;j++){
        vertdata *verti;

        verti = geomlisti->verts+j;
        if(verti->geomtype == GEOM_GEOM&&show_geom_verts == 0)continue;
        if(verti->geomtype == GEOM_ISO&&show_iso_verts == 0)continue;
        if(verti->ntriangles==0)continue;
        color = verti->triangles[0]->surf->color;
        if(last_color!=color){
          glColor3fv(color);
          last_color=color;
        }
        glVertex3fv(verti->xyz);
      }
      glEnd();
      glPopMatrix();
    }

    // draw geometry normal vectors

    if(ntris>0){  // draw faceted normals
      glPushMatrix();
      glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
      glTranslatef(-xbar0,-ybar0,-zbar0);
      glBegin(GL_LINES);
      glColor3fv(blue);
      for(j=0;j<ntris;j++){
        float *v1, *v2, *v3;
        float *xyznorm;
        tridata *trianglei;
        float xyz1[3], xyz2[3];

        trianglei = geomlisti->triangles+j;
        if(trianglei->exterior==1&&show_faces_exterior==0)continue;
        if(trianglei->exterior==0)continue;

        if(trianglei->geomtype==GEOM_GEOM&&(show_geom_normal==0||smooth_geom_normal==1))continue;
        if(trianglei->geomtype == GEOM_ISO &&(show_iso_normal == 0||smooth_iso_normal==1))continue;

        xyznorm=trianglei->tri_norm;

        v1 = trianglei->verts[0]->xyz;
        v2 = trianglei->verts[1]->xyz;
        v3 = trianglei->verts[2]->xyz;

        xyz1[0] = (v1[0] + v2[0] + v3[0])/3.0;
        xyz1[1] = (v1[1] + v2[1] + v3[1])/3.0;
        xyz1[2] = (v1[2] + v2[2] + v3[2])/3.0;
        xyz2[0] = xyz1[0] + SCALE2FDS(geom_vecfactor)*xyznorm[0];
        xyz2[1] = xyz1[1] + SCALE2FDS(geom_vecfactor)*xyznorm[1];
        xyz2[2] = xyz1[2] + SCALE2FDS(geom_vecfactor)*xyznorm[2];

        glVertex3fv(xyz1);
        glVertex3fv(xyz2);
      }
      glEnd();

      glPointSize(6.0);  // draw verts at end of vector
      glBegin(GL_POINTS);
      glColor3fv(black);
      for(j=0;j<ntris;j++){
        float *v1, *v2, *v3;
        float *xyznorm;
        tridata *trianglei;
        float xyz1[3], xyz2[3];

        trianglei = geomlisti->triangles+j;
        if(trianglei->exterior==1&&show_faces_exterior==0)continue;
        if(trianglei->exterior==0)continue;
        if(trianglei->geomtype == GEOM_GEOM && (show_geom_normal == 0 || smooth_geom_normal == 1))continue;
        if(trianglei->geomtype == GEOM_ISO&&(show_iso_normal == 0||smooth_iso_normal==1))continue;

        xyznorm=trianglei->tri_norm;

        v1 = trianglei->verts[0]->xyz;
        v2 = trianglei->verts[1]->xyz;
        v3 = trianglei->verts[2]->xyz;

        xyz1[0] = (v1[0] + v2[0] + v3[0])/3.0;
        xyz1[1] = (v1[1] + v2[1] + v3[1])/3.0;
        xyz1[2] = (v1[2] + v2[2] + v3[2])/3.0;
        xyz2[0] = xyz1[0] + SCALE2FDS(geom_vecfactor)*xyznorm[0];
        xyz2[1] = xyz1[1] + SCALE2FDS(geom_vecfactor)*xyznorm[1];
        xyz2[2] = xyz1[2] + SCALE2FDS(geom_vecfactor)*xyznorm[2];

        glVertex3fv(xyz2);
      }
      glEnd();
      glPopMatrix();
    }
    if(ntris > 0){  // draw smooth normals
      glPushMatrix();
      glScalef(SCALE2SMV(1.0), SCALE2SMV(1.0), SCALE2SMV(1.0));
      glTranslatef(-xbar0, -ybar0, -zbar0);
      glBegin(GL_LINES);
      glColor3fv(blue);
      for(j = 0; j < ntris; j++){
        tridata *trianglei;
        int k;

        trianglei = geomlisti->triangles + j;
        if(trianglei->exterior==1&&show_faces_exterior==0)continue;
        if(trianglei->exterior==0)continue;
        if(trianglei->geomtype == GEOM_GEOM && (show_geom_normal == 0 || smooth_geom_normal == 0))continue;
        if(trianglei->geomtype == GEOM_ISO && (show_iso_normal == 0 || smooth_iso_normal == 0))continue;

        for(k = 0; k < 3; k++){
          float *pk;
          float *pknorm;
          float xyz2[3];
          vertdata *vertk;

          vertk = trianglei->verts[k];
          pk = vertk->xyz;
          pknorm = trianglei->vert_norm+3*k;
          xyz2[0] = pk[0] + SCALE2FDS(geom_vecfactor)*pknorm[0];
          xyz2[1] = pk[1] + SCALE2FDS(geom_vecfactor)*pknorm[1];
          xyz2[2] = pk[2] + SCALE2FDS(geom_vecfactor)*pknorm[2];
          glVertex3fv(pk);
          glVertex3fv(xyz2);
        }
      }
      glEnd();
      glPointSize(6.0);  // draw verts at end of vector
      glBegin(GL_POINTS);
      glColor3fv(black);
      for(j = 0; j < ntris; j++){
        tridata *trianglei;
        int k;

        trianglei = geomlisti->triangles + j;
        if(trianglei->exterior==1&&show_faces_exterior==0)continue;
        if(trianglei->exterior==0)continue;
        if(trianglei->geomtype == GEOM_GEOM && (show_geom_normal == 0 || smooth_geom_normal == 0))continue;
        if(trianglei->geomtype == GEOM_ISO && (show_iso_normal == 0 || smooth_iso_normal == 0))continue;

        for(k = 0; k < 3; k++){
          float *pk;
          float *pknorm;
          float xyz2[3];
          vertdata *vertk;

          vertk = trianglei->verts[k];
          pk = vertk->xyz;
          pknorm = trianglei->vert_norm+3*k;
          xyz2[0] = pk[0] + SCALE2FDS(geom_vecfactor)*pknorm[0];
          xyz2[1] = pk[1] + SCALE2FDS(geom_vecfactor)*pknorm[1];
          xyz2[2] = pk[2] + SCALE2FDS(geom_vecfactor)*pknorm[2];
          glVertex3fv(xyz2);
        }
      }
      glEnd();
      glPopMatrix();
    }
  }
}

/* ------------------ smooth_geom_normals ------------------------ */

void smooth_geom_normals(geomlistdata *geomlisti, int geomtype){
  int i;
  float zmin, *zORIG;

      // compute average normals - method 1

  if(geomlisti->nverts > 0){
    zORIG = geomlisti->zORIG;
    zmin = zORIG[0];
    for(i = 1; i < geomlisti->nverts; i++){
      zmin = MIN(zmin, zORIG[i]);
    }
  }
  for(i = 0; i < geomlisti->nverts; i++){
    vertdata *verti;
    int k;
    float *norm, *xyz;

    verti = geomlisti->verts + i;
    norm = verti->vert_norm;
    norm[0] = 0.0;
    norm[1] = 0.0;
    norm[2] = 0.0;
    for(k = 0; k < verti->ntriangles; k++){
      float *normk;
      tridata *trianglei;

      trianglei = verti->triangles[k];
      normk = trianglei->tri_norm;
      norm[0] += normk[0];
      norm[1] += normk[1];
      norm[2] += normk[2];
    }
    ReduceToUnit(norm);

    xyz = verti->xyz;
    xyz[2] = zmin + geom_vert_exag*(zORIG[i] - zmin);
  }

      // compute average normals - method 2

  for(i = 0; i < geomlisti->ntriangles; i++){
    tridata *trianglei;
    int j;
    float *tri_normi;

    trianglei = geomlisti->triangles + i;
    tri_normi = trianglei->tri_norm;
    for(j = 0; j<3; j++){
      vertdata *vertj;
      int k;
      float *norm;

      vertj = trianglei->verts[j];
      norm = trianglei->vert_norm + 3 * j;
      if(vertj->ntriangles>0){
        norm[0] = 0.0;
        norm[1] = 0.0;
        norm[2] = 0.0;
      }
      else{
        norm[0] = 0.0;
        norm[1] = 0.0;
        norm[2] = 1.0;
      }
      for(k = 0; k<vertj->ntriangles; k++){
        tridata *trianglek;
        float *tri_normk, cosang;

        trianglek = vertj->triangles[k];
        if(trianglek->exterior == 0)continue;
        tri_normk = trianglek->tri_norm;
        cosang = DOT3(tri_normk, tri_normi)/(NORM3(tri_normk)*NORM3(tri_normi));
        if(use_max_angle==0||geomtype==GEOM_ISO||cosang>cos_geom_max_angle){ // smooth using all triangles if an isosurface
          norm[0] += tri_normk[0];
          norm[1] += tri_normk[1];
          norm[2] += tri_normk[2];
        }
      }
      ReduceToUnit(norm);
    }
  }
}

/* ------------------ update_geom_normals ------------------------ */

void update_geom_normals(void){
  int j, ii;

  for(j = 0; j < ngeominfoptrs; j++){
    geomdata *geomi;
    int iend;

    geomi = geominfoptrs[j];
    if(geomi->loaded == 0 || geomi->display == 0)continue;
    if(geomi->geomtype != GEOM_GEOM&&geomi->geomtype != GEOM_ISO)continue;

    iend = geomi->ntimes;
    if(geomi->currentframe != NULL)iend = 1;

    for(ii = -1; ii < iend; ii++){
      geomlistdata *geomlisti;

      if(ii == -1 || geomi->currentframe == NULL){
        geomlisti = geomi->geomlistinfo + ii;
      }
      else{
        geomlisti = geomi->currentframe;
      }
      smooth_geom_normals(geomlisti,geomi->geomtype);
    }
  }
}

/* ------------------ update_triangles ------------------------ */

void update_triangles(int flag,int update){
  int j, ii, ntimes;

  if(update==GEOM_UPDATE_NORMALS){
    update_geom_normals();
    return;
  }
  for(j=0;j<ngeominfoptrs;j++){
    geomdata *geomi;
    float *xyzptr[3];
    float *xyznorm;
    int i;
    int iend;

    geomi = geominfoptrs[j];
    if(geomi->loaded==0||geomi->display==0)continue;
    if(geomi->geomtype != GEOM_GEOM&&geomi->geomtype!=GEOM_ISO)continue;

    iend = geomi->ntimes;
    if(geomi->currentframe != NULL)iend = 1;

    for(ii=-1;ii<iend;ii++){
      geomlistdata *geomlisti;
      int ntriangles;
      tridata **triangles;

      if(ii==-1||geomi->currentframe==NULL){
        geomlisti = geomi->geomlistinfo + ii;
      }
      else{
        geomlisti = geomi->currentframe;
      }
      for(i=0;i<geomlisti->ntriangles;i++){
        tridata *trianglei;

        trianglei = geomlisti->triangles+i;

        xyzptr[0] = trianglei->verts[0]->xyz;
        xyzptr[1] = trianglei->verts[1]->xyz;
        xyzptr[2] = trianglei->verts[2]->xyz;
        xyznorm = trianglei->tri_norm;
        CalcTriNormal(xyzptr[0],xyzptr[1],xyzptr[2],xyznorm);
      }

      for(i=0;i<geomlisti->nverts;i++){
        vertdata *verti;

        verti = geomlisti->verts + i;
        verti->ntriangles=0;
        verti->itriangle=0;
        verti->on_mesh_boundary = 0;
      }
      for(i=0;i<geomlisti->ntriangles;i++){
        tridata *trianglei;

        trianglei = geomlisti->triangles+i;
        trianglei->verts[0]->ntriangles++;
        trianglei->verts[1]->ntriangles++;
        trianglei->verts[2]->ntriangles++;
      }

      // count number of triangles

      ntriangles = 0;
      for(i = 0; i<geomlisti->nverts; i++){
        vertdata *verti;

        verti = geomlisti->verts + i;
        ntriangles += verti->ntriangles;
      }

      // allocate triangle pointers

      FREEMEMORY(geomlisti->triangleptrs);
      if(ntriangles>0){
        NewMemoryMemID((void **)&triangles, ntriangles*sizeof(tridata *), geomi->memory_id);
        geomlisti->triangleptrs = triangles;
      }

      // assign triangle pointers to verts

      for(i = 0; i<geomlisti->nverts; i++){
        vertdata *verti;

        verti = geomlisti->verts + i;
        if(verti->ntriangles>0){
          verti->triangles = triangles;
          triangles += verti->ntriangles;
        }
      }
      for(i=0;i<geomlisti->ntriangles;i++){
        tridata *trianglei;
        vertdata *verti;

        trianglei = geomlisti->triangles+i;
        verti = trianglei->verts[0];
        verti->triangles[verti->itriangle++]=trianglei;
        verti = trianglei->verts[1];
        verti->triangles[verti->itriangle++]=trianglei;
        verti = trianglei->verts[2];
        verti->triangles[verti->itriangle++]=trianglei;
      }

      smooth_geom_normals(geomlisti,geomi->geomtype);

    }
  }

  // smooth normals at mesh boundaries

  if(ngeominfoptrs>0){
    vertdata **surface_verts = NULL;
    int *match_verts = NULL;

    if(flag == GEOM_STATIC){
      ntimes = 0;
    }
    else{
      ntimes = geominfoptrs[0]->ntimes;
    }
    for(ii = -1; ii<ntimes; ii++){
      int nsurface_verts;

  // identify and count verts on mesh surfaces

      nsurface_verts = 0;
      for(j = 0; j<ngeominfoptrs; j++){
        geomlistdata *geomlisti;
        int  i;
        geomdata *geomj;

        geomj = geominfoptrs[j];
        if(geomj->geomtype != GEOM_GEOM&&geomj->geomtype!=GEOM_ISO)continue;
        geomlisti = geomj->geomlistinfo+ii;
        for(i = 0; i<geomlisti->nverts; i++){
          vertdata *verti;

          verti = geomlisti->verts+i;
          verti->on_mesh_boundary = on_mesh_boundary(verti->xyz);
          if(verti->on_mesh_boundary==1)nsurface_verts++;
        }
      }

  // copy surface verts into an array

      if(nsurface_verts>0){
        int isurf,iii;

        isurf = 0;
        FREEMEMORY(surface_verts);
        FREEMEMORY(match_verts);
        NewMemory((void **)&surface_verts, nsurface_verts*sizeof(vertdata *));
        NewMemory((void **)&match_verts, nsurface_verts*sizeof(int));
        for(j = 0; j<ngeominfoptrs; j++){
          geomlistdata *geomlisti;
          int  i;
          geomdata *geomj;

          geomj = geominfoptrs[j];
          if(geomj->geomtype != GEOM_GEOM&&geomj->geomtype != GEOM_ISO)continue;
          geomlisti = geomj->geomlistinfo + ii;
          for(i = 0; i<geomlisti->nverts; i++){
            vertdata *verti;

            verti = geomlisti->verts+i;
            if(verti->on_mesh_boundary==1){
              if(isurf<nsurface_verts){
                surface_verts[isurf] = verti;
                match_verts[isurf] = -1;
                isurf++;
              }
            }
          }
        }

        // average normals

        for(iii = 0; iii<nsurface_verts; iii++){
          int jjj;
          vertdata *verti;
          float *xyzi, *normi;
          float avgnorm[3];

          if(match_verts[iii]>=0)continue;
          verti = surface_verts[iii];
          xyzi = verti->xyz;
          normi = verti->vert_norm;
          avgnorm[0] = normi[0];
          avgnorm[1] = normi[1];
          avgnorm[2] = normi[2];
          match_verts[iii] = iii;
          for(jjj = iii+1; jjj<nsurface_verts; jjj++){
            vertdata *vertj;
            float *xyzj, *normj;

            if(match_verts[jjj]>=0)continue;
            vertj = surface_verts[jjj];
            xyzj = vertj->xyz;
            normj = vertj->vert_norm;
#define POINTEPS 0.001
            if(ABS(xyzi[0]-xyzj[0])<POINTEPS&&ABS(xyzi[1]-xyzj[1])<POINTEPS&&ABS(xyzi[2]-xyzj[2])<POINTEPS){
              match_verts[jjj] = iii;
              avgnorm[0] += normj[0];
              avgnorm[1] += normj[1];
              avgnorm[2] += normj[2];
            }
          }
          ReduceToUnit(avgnorm);
          for(jjj = iii; jjj<nsurface_verts; jjj++){
            if(match_verts[jjj] == match_verts[iii]){
              vertdata *vertj;
              float *normj;

              vertj = surface_verts[jjj];
              normj = vertj->vert_norm;
              normj[0] = avgnorm[0];
              normj[1] = avgnorm[1];
              normj[2] = avgnorm[2];
            }
          }
        }
      }
    }
    FREEMEMORY(surface_verts);
    FREEMEMORY(match_verts);
  }

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
}

#define FORTREAD(var,count,STREAM) FSEEK(STREAM,4,SEEK_CUR);\
                           returncode=fread(var,4,count,STREAM);\
                           if(returncode!=count)returncode=0;\
                           if(endianswitch==1&&returncode!=0)endian_switch(var,count);\
                           FSEEK(STREAM,4,SEEK_CUR)

#define FORTREADBR(var,count,STREAM) FORTREAD(var,(count),STREAM);if(returncode==0)break;

/* ------------------ read_geom_header0 ------------------------ */

void read_geom_header0(geomdata *geomi, int *geom_frame_index, int *ntimes_local){
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
  int icount;

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
    NewMemoryMemID((void **)&float_vals,nfloat_vals*sizeof(float),geomi->memory_id);
    FORTREAD(float_vals,nfloat_vals,stream);
    geomi->float_vals=float_vals;
    geomi->nfloat_vals=nfloat_vals;
  }

  // integer header

  FORTREAD(&nint_vals,1,stream);
  if(nint_vals>0){
    NewMemoryMemID((void **)&int_vals,nint_vals*sizeof(float),geomi->memory_id);
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
  icount=-1;
  for(;;){
    FORTREADBR(times_local,2,stream);
    icount++;
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

    if(geom_frame_index==NULL){
      if(use_tload_begin == 1 && times_local[0] < tload_begin)continue;
      if(use_tload_skip == 1 && tload_skip>1 && icount%tload_skip!=0)continue;
      if(use_tload_end == 1 && times_local[0] > tload_end)break;
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
  int nvertfacesvolumes[3];
  int nt;
  int returncode;
  int version;
  int nverts=0, ntris=0, nvolumes=0;
  int first_frame_static;
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
    FORTREADBR(nvertfacesvolumes,3,stream);
    nverts=nvertfacesvolumes[0];
    ntris=nvertfacesvolumes[1];
    nvolumes=nvertfacesvolumes[2];

    // vertices

    if(nverts>0){
      FSEEK(stream,4+3*nverts*4+4,SEEK_CUR); // skip vertices
    }

    // faces

    if(ntris>0){
      FSEEK(stream,4+3*ntris*4+4,SEEK_CUR); // skip triangles
      FSEEK(stream,4+ntris*4+4,SEEK_CUR);   // skip surf
      FSEEK(stream,4+6*ntris*4+4,SEEK_CUR); // skip textures
    }

    // volumes

    if(nvolumes>0){
      FSEEK(stream,4+4*nvolumes*4+4,SEEK_CUR); // skip volumes
      FSEEK(stream,4+nvolumes*4+4,SEEK_CUR);   // skip matl
    }
    nt++;
  }
  *ntimes_local=nt;
  fclose(stream);
}

/* ------------------ read_geom_header ------------------------ */

void read_geom_header(geomdata *geomi, int *geom_frame_index, int *ntimes_local){
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
    read_geom_header0(geomi, geom_frame_index, ntimes_local);
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
    read_geom(geomi,LOAD,GEOM_GEOM,NULL,&errorcode);
  }
  for(i = 0; i < ngeomdiaginfo; i++){
    geomdiagdata *geomdiagi;

    geomdiagi = geomdiaginfo + i;
    read_geom(geomdiagi->geom, LOAD, GEOM_GEOM, NULL, &errorcode);
  }
}

/* ------------------ read_geom0 ------------------------ */

void read_geom0(geomdata *geomi, int load_flag, int type, int *geom_frame_index, int *errorcode){
  FILE *stream;
  int one=1, endianswitch=0;
  int returncode;
  int ntimes_local;
  int version;
  int nvertfacesvolumes[3];
  int nfloat_vals, nint_vals;
  int iframe, icount;

  FreeAllMemory(geomi->memory_id);
  geomi->geomlistinfo = NULL;
  geomi->currentframe = NULL;
  geomi->nfloat_vals=0;
  geomi->nint_vals=0;

  if(load_flag==UNLOAD){
    geomi->loaded=0;
    geomi->display=0;
    return;
  }

  read_geom_header(geomi,geom_frame_index,&ntimes_local);
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
  NewMemoryMemID((void **)&geomi->geomlistinfo_0,(ntimes_local+1)*sizeof(geomlistdata),geomi->memory_id);
  geomi->geomlistinfo=geomi->geomlistinfo_0+1;
  NewMemoryMemID((void **)&geomi->times,ntimes_local*sizeof(float),geomi->memory_id);

  icount=-1;
  for(iframe=-1;iframe<ntimes_local;){
    float times_local[2];
    geomlistdata *geomlisti;
    int nverts, ntris;
    int  skipframe;
    vertdata *verts;

    geomlisti = geomi->geomlistinfo+iframe;
    geomlisti->verts=NULL;
    geomlisti->triangles=NULL;
    geomlisti->triangleptrs = NULL;
    geomlisti->volumes=NULL;
    geomlisti->nverts=0;
    geomlisti->ntriangles=0;
    geomlisti->nvolumes=0;
    skipframe = 0;

    if(iframe>=0){
      FORTREADBR(times_local,2,stream);
      icount++;
      if(geom_frame_index == NULL){
        if(use_tload_begin == 1 && times_local[0] < tload_begin)skipframe = 1;
        if(use_tload_skip == 1 && tload_skip>1 && icount%tload_skip != 0)skipframe = 1;
        if(use_tload_end == 1 && times_local[0] > tload_end)skipframe = 1;
        if(skipframe == 0)geomi->times[iframe] = times_local[0];
      }
      else{
        if(iframe!=*geom_frame_index)skipframe = 1;
        geomi->times[iframe] = times_local[0];
        if(skipframe == 0)geomi->currentframe = geomlisti;
      }
    }
    FORTREADBR(nvertfacesvolumes,2,stream);
    nverts=nvertfacesvolumes[0];
    ntris=nvertfacesvolumes[1];
    if(skipframe==0&&iframe>=0){
      PRINTF("time=%.2f triangles: %i\n",times_local[0],ntris);
    }
    if(skipframe==1){
      int file_offset = 0;
      if(nverts>0)file_offset += 4+3*nverts*4+4;
      if(ntris>0)file_offset += (4+3*ntris*4+4)+(4+ntris*4+4);
      if(file_offset>0)FSEEK(stream, file_offset, SEEK_CUR);
    }
    if(skipframe==0&&nverts>0){
      int ii;
      float *xyz=NULL;
      float *zORIG;

      if(iframe<0)PRINTF("static geometry\n");
      NewMemory((void **)&xyz,3*nverts*sizeof(float));
      NewMemoryMemID((void **)&verts,nverts*sizeof(vertdata),geomi->memory_id);
      NewMemory((void **)&zORIG, nverts*sizeof(float));
      geomlisti->zORIG = zORIG;
      geomlisti->verts = verts;
      geomlisti->nverts=nverts;
      FORTREADBR(xyz,3*nverts,stream);
      for(ii=0;ii<nverts;ii++){
        verts[ii].xyz[0]=xyz[3*ii];
        verts[ii].xyz[1]=xyz[3*ii+1];
        verts[ii].xyz[2]=xyz[3*ii+2];
        zORIG[ii] = xyz[3 * ii+2];
      }
      FREEMEMORY(xyz);
    }
    if(skipframe==0&&ntris>0){
      int *surf_ind=NULL,*ijk=NULL;
      int ii;
      int offset=0;
      tridata *triangles;

      NewMemoryMemID((void **)&triangles,ntris*sizeof(tridata),geomi->memory_id);
      NewMemory((void **)&ijk,3*ntris*sizeof(int));
      NewMemory((void **)&surf_ind,ntris*sizeof(int));
      geomlisti->triangles=triangles;
      geomlisti->ntriangles=ntris;
      FORTREADBR(ijk,3*ntris,stream);
      FORTREADBR(surf_ind,ntris,stream);
      if(type==GEOM_ISO)offset=nsurfinfo;
      for(ii=0;ii<ntris;ii++){
        surfdata *surfi;

        triangles[ii].verts[0]=verts+ijk[3*ii]-1;
        triangles[ii].verts[1]=verts+ijk[3*ii+1]-1;
        triangles[ii].verts[2]=verts+ijk[3*ii+2]-1;
        surfi = surfinfo+CLAMP(surf_ind[ii]+offset, nsurfinfo+1, nsurfinfo+MAX_ISO_COLORS);
        triangles[ii].surf=surfi;
        triangles[ii].textureinfo=NULL;
      }
      FREEMEMORY(ijk);
      FREEMEMORY(surf_ind);
    }

    if(skipframe==0||geom_frame_index!=NULL){
      // add decimation code here
      iframe++;
    }
    if(geom_frame_index==NULL&&use_tload_end == 1 && times_local[0] > tload_end)break;
  }
  geomi->loaded = 1;
  geomi->display=1;
}

/* ------------------ read_geom2 ------------------------ */

void read_geom2(geomdata *geomi, int load_flag, int type, int *errorcode){
  FILE *stream;
  int one=1, endianswitch=0;
  int returncode;
  int ntimes_local;
  int i;
  vertdata *verts;
  tridata *triangles;
  tetdata *volumes;
  int version;
  int nvertfacesvolumes[3];
  int nheaders[3], nfloat_vals, nint_vals, first_frame_static;

  FreeAllMemory(geomi->memory_id);
  geomi->geomlistinfo=NULL;
  geomi->currentframe = NULL;
  geomi->nfloat_vals=0;
  geomi->nint_vals=0;

  if(load_flag==UNLOAD){
    geomi->loaded=0;
    geomi->display=0;
    return;
  }

  read_geom_header(geomi,NULL,&ntimes_local);
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
  NewMemoryMemID((void **)&geomi->geomlistinfo_0,(ntimes_local+1)*sizeof(geomlistdata),geomi->memory_id);
  geomi->geomlistinfo=geomi->geomlistinfo_0+1;
  if(ntimes_local>0)NewMemoryMemID((void **)&geomi->times,ntimes_local*sizeof(float),geomi->memory_id);

  for(i=-1;i<ntimes_local;i++){
    float time_local;
    geomlistdata *geomlisti;
    int nverts, ntris, nvolumes;

    geomlisti = geomi->geomlistinfo+i;
    geomlisti->verts=NULL;
    geomlisti->triangles=NULL;
    geomlisti->triangleptrs = NULL;
    geomlisti->volumes=NULL;
    geomlisti->nverts=0;
    geomlisti->ntriangles=0;
    geomlisti->nvolumes=0;

    if(first_frame_static==0&&i==-1)continue;

    FORTREADBR(&time_local,1,stream);
    if(i>=0)geomi->times[i]=time_local;

    FORTREADBR(nvertfacesvolumes,3,stream);
    nverts=nvertfacesvolumes[0];
    ntris=nvertfacesvolumes[1];
    nvolumes=nvertfacesvolumes[2];
    if(nvolumes>0)have_volume=1;

    if(i>=0){
      PRINTF("time=%.2f triangles: %i\n",time_local,ntris);
    }
    if(nverts>0){
      int ii;
      float *xyz=NULL;
      float *zORIG;

      if(i<0)PRINTF("static geometry\n");
      NewMemory((void **)&xyz,3*nverts*sizeof(float));
      NewMemory((void **)&zORIG,nverts*sizeof(float));
      NewMemoryMemID((void **)&verts,nverts*sizeof(vertdata),geomi->memory_id);
      geomlisti->verts=verts;
      geomlisti->zORIG=zORIG;
      geomlisti->nverts=nverts;
      FORTREADBR(xyz,3*nverts,stream);
      for(ii=0;ii<nverts;ii++){
        verts[ii].xyz[0]=xyz[3*ii];
        verts[ii].xyz[1]=xyz[3*ii+1];
        verts[ii].xyz[2]=xyz[3*ii+2];
        zORIG[ii] = xyz[3*ii+2];
      }
      FREEMEMORY(xyz);
    }
    if(ntris>0){
      int *surf_ind=NULL,*ijk=NULL;
      float *texture_coords=NULL;
      int ii;
      int offset=0;

      NewMemoryMemID((void **)&triangles,ntris*sizeof(tridata),geomi->memory_id);
      NewMemory((void **)&ijk,3*ntris*sizeof(int));
      NewMemory((void **)&surf_ind,ntris*sizeof(int));
      NewMemory((void **)&texture_coords,6*ntris*sizeof(float));
      geomlisti->triangles=triangles;
      geomlisti->ntriangles=ntris;
      FORTREADBR(ijk,3*ntris,stream);
      FORTREADBR(surf_ind,ntris,stream);
      FORTREADBR(texture_coords,6*ntris,stream);
      if(type==GEOM_ISO)offset=nsurfinfo;
      for(ii=0;ii<ntris;ii++){
        surfdata *surfi;
        int k;

        for(k=0;k<3;k++){
          triangles[ii].verts[k]=verts+ijk[3*ii+k]-1;
        }

        for(k=0;k<6;k++){
          triangles[ii].tverts[k]=texture_coords[6*ii+k];
        }

        surfi=surfinfo + surf_ind[ii]+offset;
        triangles[ii].surf=surfi;
        triangles[ii].insolid = surf_ind[ii];
        triangles[ii].textureinfo=surfi->textureinfo;
      }

      FREEMEMORY(ijk);
      FREEMEMORY(surf_ind);
      FREEMEMORY(texture_coords);
    }
    if(nvolumes>0){
      int ii;
      int *ijk;
      int *matl_ind=NULL;

      NewMemoryMemID((void **)&volumes,nvolumes*sizeof(tetdata),geomi->memory_id);
      geomlisti->volumes=volumes;
      NewMemory((void **)&ijk,4*nvolumes*sizeof(int));
      FORTREADBR(ijk,4*nvolumes,stream);
      for(ii=0;ii<nvolumes;ii++){
        int k;

        for(k=0;k<4;k++){
          volumes[ii].verts[k]=verts+ijk[4*ii+k]-1;
        }
      }
      FREEMEMORY(ijk);
      NewMemory((void **)&matl_ind,nvolumes*sizeof(int));
      FORTREADBR(matl_ind,nvolumes,stream);
      for(ii=0;ii<nvolumes;ii++){
        matldata *matli;
        int index;

        index = CLAMP(matl_ind[ii],0,nmatlinfo-1);
        matli=matlinfo + index;
        volumes[ii].matl=matli;
      }
      FREEMEMORY(matl_ind);
      geomlisti->nvolumes=nvolumes;
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

/* ------------------ compare_faces ------------------------ */

int compare_faces(const void *arg1, const void *arg2){
  tridata *face1, *face2;
  int *verts1, *verts2;
  int v1[3], v2[3];

  face1 = triangle_list + *(int *)arg1;
  face2 = triangle_list + *(int *)arg2;
  verts1 = face1->vert_index;
  verts2 = face2->vert_index;

  v1[0] = MIN(verts1[0], MIN(verts1[1], verts1[2]));
  v1[2] = MAX(verts1[0], MAX(verts1[1], verts1[2]));
  v1[1] = verts1[0]+verts1[1]+verts1[2]-v1[0]-v1[2];

  v2[0] = MIN(verts2[0], MIN(verts2[1], verts2[2]));
  v2[2] = MAX(verts2[0], MAX(verts2[1], verts2[2]));
  v2[1] = verts2[0]+verts2[1]+verts2[2]-v2[0]-v2[2];

  if(v1[0]<v2[0])return -1;
  if(v1[0]>v2[0])return 1;

  if(v1[1]<v2[1])return -1;
  if(v1[1]>v2[1])return 1;

  if(v1[2]<v2[2])return -1;
  if(v1[2]>v2[2])return 1;
  return 0;
}

/* ------------------ compare_volume_faces ------------------------ */

int compare_volume_faces(const void *arg1, const void *arg2){
  int face1, face2;
  tetdata *vol1, *vol2;
  int *verts1, *verts2;
  int v1[3], v2[3];

  face1=*(int *)arg1;
  face2=*(int *)arg2;
  vol1 = volume_list + face1/4;
  vol2 = volume_list + face2/4;
  face1 %= 4;
  face2 %= 4;
  verts1 = vol1->faces+3*face1;
  verts2 = vol2->faces+3*face2;

  v1[1]=MIN(verts1[1],verts1[2]);
  v1[2]=MAX(verts1[1],verts1[2]);

  v2[1]=MIN(verts2[1],verts2[2]);
  v2[2]=MAX(verts2[1],verts2[2]);

  if(verts1[0]<verts2[0])return -1;
  if(verts1[0]>verts2[0])return 1;

  if(v1[1]<v2[1])return -1;
  if(v1[1]>v2[1])return 1;

  if(v1[2]<v2[2])return -1;
  if(v1[2]>v2[2])return 1;
  return 0;
}

/* ------------------ classify_geom ------------------------ */

void classify_geom(geomdata *geomi,int *geom_frame_index){
  int i, iend;

  iend = geomi->ntimes;
  if(geom_frame_index!=NULL)iend=1;

  for(i = -1; i<iend; i++){
    geomlistdata *geomlisti;
    int nverts, nvolumes, ntriangles;
    int j;
    vertdata *vertbase;

    geomlisti = geomi->geomlistinfo+i;
    if(i!=-1&&geom_frame_index!=NULL)geomlisti = geomi->geomlistinfo+(*geom_frame_index);

    nverts=geomlisti->nverts;
    nvolumes=geomlisti->nvolumes;
    ntriangles = geomlisti->ntriangles;
    if(nverts==0||geomlisti->verts==NULL)continue;
    vertbase = geomlisti->verts;
    for(j=0;j<nvolumes;j++){
      tetdata *tetrai;
      int *vert_index;
      vertdata **verts;
      int *faces;

      tetrai = geomlisti->volumes+j;
      vert_index = tetrai->vert_index;
      verts = tetrai->verts;
      faces = tetrai->faces;
      tetrai->exterior[0]=1;
      tetrai->exterior[1]=1;
      tetrai->exterior[2]=1;
      tetrai->exterior[3]=1;
      vert_index[0] = verts[0]-vertbase;
      vert_index[1] = verts[1]-vertbase;
      vert_index[2] = verts[2]-vertbase;
      vert_index[3] = verts[3]-vertbase;

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
    if(nvolumes>0){
      int *facelist=NULL,nfacelist;

      nfacelist=4*nvolumes;
      volume_list=geomlisti->volumes;
      NewMemory((void **)&facelist,4*nfacelist*sizeof(int));
      for(j=0;j<nfacelist;j++){
        facelist[j]=j;
      }
      qsort(facelist,nfacelist,sizeof(int),compare_volume_faces);
      for(j=1;j<nfacelist;j++){
        int face1, face2;
        tetdata *vol1, *vol2;
        int *verts1, *verts2;

        face1=facelist[j-1];
        face2=facelist[j];
        vol1 = volume_list + face1/4;
        vol2 = volume_list + face2/4;
        face1 %= 4;
        face2 %= 4;
        verts1 = vol1->faces+3*face1;
        verts2 = vol2->faces+3*face2;
        if(verts1[0]!=verts2[0])continue;
        if(MIN(verts1[1],verts1[2])!=MIN(verts2[1],verts2[2]))continue;
        if(MAX(verts1[1],verts1[2])!=MAX(verts2[1],verts2[2]))continue;
        vol1->exterior[face1]=0;
        vol2->exterior[face2]=0;
      }

      FREEMEMORY(facelist);
    }
    if(ntriangles > 0){
      int *facelist = NULL, nfacelist;

      nfacelist = ntriangles;
      triangle_list = geomlisti->triangles;
      NewMemory((void **)&facelist, nfacelist*sizeof(int));
      for(j = 0; j < nfacelist; j++){
        tridata *trij;
        int *vert_index;

        trij = geomlisti->triangles + j;
        trij->exterior = 1;
        facelist[j] = j;
        vert_index = trij->vert_index;
        vert_index[0] = trij->verts[0] - vertbase;
        vert_index[1] = trij->verts[1] - vertbase;
        vert_index[2] = trij->verts[2] - vertbase;
      }
      qsort(facelist, nfacelist, sizeof(int), compare_faces);
      for(j = 1; j < nfacelist; j++){
        if(compare_faces(facelist + j, facelist + j - 1) == 0){
          tridata *trij, *trijm1;

          trij = geomlisti->triangles + facelist[j];
          trij->exterior = 0;

          trijm1 = geomlisti->triangles + facelist[j - 1];
          trijm1->exterior = 0;
         }
      }

      FREEMEMORY(facelist);
    }
  }
}

/* ------------------ read_geom ------------------------ */

void read_geom(geomdata *geomi, int load_flag, int type, int *geom_frame_index, int *errorcode){
  FILE *stream;
  int version;
  int returncode;
  int one=0,endianswitch=0;

  if(geomi->file==NULL)return;
  stream = fopen(geomi->file,"rb");
  if(stream==NULL)return;
  FSEEK(stream,4,SEEK_CUR);fread(&one,4,1,stream);FSEEK(stream,4,SEEK_CUR);
  if(one!=1)endianswitch=1;
  FORTREAD(&version,1,stream);
  fclose(stream);

  if(version<=1){
    read_geom0(geomi,load_flag,type,geom_frame_index,errorcode);
  }
  else{
    read_geom2(geomi,load_flag,type,errorcode);
  }
  if(load_flag==LOAD)classify_geom(geomi,geom_frame_index);
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
  if(patchi->filetype!=PATCH_GEOMETRY)return;
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

  FORTgetembeddatasize(file, &ntimes_local, &nvals, &error, lenfile);

  if(nvals>0&&ntimes_local>0){
    NewMemory((void **)&patchi->geom_nstatics,ntimes_local*sizeof(int));
    NewMemory((void **)&patchi->geom_ndynamics,ntimes_local*sizeof(int));
    NewMemory((void **)&patchi->geom_times,ntimes_local*sizeof(float));
    NewMemory((void **)&patchi->geom_ivals_static,ntimes_local*sizeof(int *));
    NewMemory((void **)&patchi->geom_ivals_dynamic,ntimes_local*sizeof(int *));
    NewMemory((void **)&patchi->geom_vals,nvals*sizeof(float));
    NewMemory((void **)&patchi->geom_ivals,nvals*sizeof(char));
  }
  FORTgetembeddata(file, &ntimes_local, &nvals, patchi->geom_times,
    patchi->geom_nstatics, patchi->geom_ndynamics, patchi->geom_vals, &redirect, &error, lenfile);

  reset_histogram(patchi->histogram);

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

/* ------------------ draw_test_clip ------------------------ */

void draw_test_clip(void){
  float *xmin, *xmax, *ymin, *ymax, *zmin, *zmax;
  unsigned char cube0color[4] ={255,  0,  0,255};
  unsigned char cube1color[4] ={128,  0,  0,255};
  unsigned char cube2color[4] ={  0,255,  0,255};
  unsigned char cube3color[4] ={  0,128,  0,255};
  unsigned char cube4color[4] ={  0,  0,255,255};
  unsigned char cube5color[4] ={  0,  0,128,255};
  unsigned char tetra0color[4]={  0,255,255,255};
  unsigned char tetra1color[4]={255,  0,255,255};
  unsigned char tetra2color[4]={255,255,  0,255};
  unsigned char tetra3color[4]={ 64, 64, 64,255};
  clipdata tetra_clipinfo, box_clipinfo;
  float *v1, *v2, *v3, *v4;
  int nverts;
  int faces[600], npolys, nfaces;
  int which_poly[200];
  float verts[600];

  box_state = b_state+1;
  v1 = tetra_vertices;
  v2 = v1 + 3;
  v3 = v2 + 3;
  v4 = v3 + 3;

  {
    float specular[4]={0.4,0.4,0.4,1.0};

    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SHININESS,&block_shininess);
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,block_ambient2);
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,specular);
    glEnable(GL_COLOR_MATERIAL);
  }

  initTetraClipInfo(&tetra_clipinfo,v1,v2,v3,v4);

  xmin = box_bounds;
  xmax = box_bounds+1;
  ymin = box_bounds+2;
  ymax = box_bounds+3;
  zmin = box_bounds+4;
  zmax = box_bounds+5;
  {
    int i;

    FORTgetverts(box_bounds, v1, v2, v3, v4, verts, &nverts, faces, face_id, which_poly, &nfaces, &npolys, b_state);
    if(update_volbox_controls==1){
      for(i=0;i<10;i++){
        face_vis[i]=0;
      }
      for(i=0;i<npolys;i++){
        face_vis[CLAMP(face_id[i],0,9)]=1;
      }
      for(i=0;i<npolys;i++){
        int face;
        int tet_face;

        face = face_id[i];
        if(face<6){
          tet_face=box_state[face];
          if(tet_face>-1){
            ASSERT(tet_face>=0&&tet_face<4);
            face_vis[tet_face+6]=1;
          }
        }
      }
      Volume_CB(2);
    }

    if(npolys>10){
      PRINTF("***error: nface=%i should not be bigger than 10\n",npolys);
    }
    initBoxClipInfo(&box_clipinfo,*xmin,*xmax,*ymin,*ymax,*zmin,*zmax);
    if(box_state[0]!=-1)box_clipinfo.clip_xmin=0;
    if(box_state[1]!=-1)box_clipinfo.clip_xmax=0;
    if(box_state[2]!=-1)box_clipinfo.clip_ymin=0;
    if(box_state[3]!=-1)box_clipinfo.clip_ymax=0;
    if(box_state[4]!=-1)box_clipinfo.clip_zmin=0;
    if(box_state[5]!=-1)box_clipinfo.clip_zmax=0;
#define  EPSBOX (-0.0001)
    box_clipinfo.xmin+=EPSBOX;
    box_clipinfo.xmax-=EPSBOX;
    box_clipinfo.ymin+=EPSBOX;
    box_clipinfo.ymax-=EPSBOX;
    box_clipinfo.zmin+=EPSBOX;
    box_clipinfo.zmax-=EPSBOX;
  }

  glPushMatrix();
  glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
  glTranslatef(-xbar0,-ybar0,-zbar0);

  // draw box

  setClipPlanes(&tetra_clipinfo,CLIP_ON_DENORMAL);
  glPushMatrix();
  glTranslatef(*xmin,*ymin,*zmin);
  glScalef(ABS(*xmax-*xmin),ABS(*ymax-*ymin),ABS(*zmax-*zmin));
  {
    glBegin(GL_QUADS);

    if(box_state[4]==-1&&tetrabox_vis[4]==1){
      glNormal3f( 0.0, 0.0,-1.0);
      glColor3ubv(cube4color);
      glVertex3f( 0.0,0.0,0.0);  // 1
      glVertex3f( 0.0,1.0,0.0);  // 4
      glVertex3f( 1.0,1.0,0.0);  // 3
      glVertex3f( 1.0,0.0,0.0);  // 2

      glVertex3f( 0.0,0.0,0.0);  // 1
      glVertex3f( 1.0,0.0,0.0);  // 2
      glVertex3f( 1.0,1.0,0.0);  // 3
      glVertex3f( 0.0,1.0,0.0);  // 4
    }

    if(box_state[5]==-1&&tetrabox_vis[5]==1){
      glNormal3f(0.0,0.0,1.0);
      glColor3ubv(cube5color);
      glVertex3f(0.0,0.0,1.0);  // 5
      glVertex3f(1.0,0.0,1.0);  // 6
      glVertex3f(1.0,1.0,1.0);  // 7
      glVertex3f(0.0,1.0,1.0);  // 8

      glVertex3f(0.0,0.0,1.0);  // 5
      glVertex3f(0.0,1.0,1.0);  // 8
      glVertex3f(1.0,1.0,1.0);  // 7
      glVertex3f(1.0,0.0,1.0);  // 6
    }

    if(box_state[2]==-1&&tetrabox_vis[2]==1){
      glNormal3f(0.0,-1.0,0.0);
      glColor3ubv(cube2color);
      glVertex3f(0.0,0.0,0.0);  // 1
      glVertex3f(1.0,0.0,0.0);  // 2
      glVertex3f(1.0,0.0,1.0);  // 6
      glVertex3f(0.0,0.0,1.0);  // 5

      glVertex3f(0.0,0.0,0.0);  // 1
      glVertex3f(0.0,0.0,1.0);  // 5
      glVertex3f(1.0,0.0,1.0);  // 6
      glVertex3f(1.0,0.0,0.0);  // 2
    }

    if(box_state[3]==-1&&tetrabox_vis[3]==1){
      glNormal3f(0.0,1.0,0.0);
      glColor3ubv(cube3color);
      glVertex3f(1.0,1.0,0.0);  // 3
      glVertex3f(0.0,1.0,0.0);  // 4
      glVertex3f(0.0,1.0,1.0);  // 8
      glVertex3f(1.0,1.0,1.0);  // 7

      glVertex3f(1.0,1.0,0.0);  // 3
      glVertex3f(1.0,1.0,1.0);  // 7
      glVertex3f(0.0,1.0,1.0);  // 8
      glVertex3f(0.0,1.0,0.0);  // 4
    }

    if(box_state[0]==-1&&tetrabox_vis[0]==1){
      glNormal3f(-1.0,0.0,0.0);
      glColor3ubv(cube0color);
      glVertex3f(0.0,0.0,0.0);  // 1
      glVertex3f(0.0,0.0,1.0);  // 5
      glVertex3f(0.0,1.0,1.0);  // 8
      glVertex3f(0.0,1.0,0.0);  // 4

      glVertex3f(0.0,0.0,0.0);  // 1
      glVertex3f(0.0,1.0,0.0);  // 4
      glVertex3f(0.0,1.0,1.0);  // 8
      glVertex3f(0.0,0.0,1.0);  // 5
    }

    if(box_state[1]==-1&&tetrabox_vis[1]==1){
      glNormal3f(1.0,0.0,0.0);
      glColor3ubv(cube1color);
      glVertex3f(1.0,0.0,0.0);  // 2
      glVertex3f(1.0,1.0,0.0);  // 3
      glVertex3f(1.0,1.0,1.0);  // 7
      glVertex3f(1.0,0.0,1.0);  // 6

      glVertex3f(1.0,0.0,0.0);  // 2
      glVertex3f(1.0,0.0,1.0);  // 6
      glVertex3f(1.0,1.0,1.0);  // 7
      glVertex3f(1.0,1.0,0.0);  // 3
    }
    glEnd();
  }
#define EPS 0.02
  glPopMatrix();
  glPopMatrix();

  // draw tetrahedron

  glPushMatrix();
  glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
  glTranslatef(-xbar0,-ybar0,-zbar0);
  setClipPlanes(&box_clipinfo,CLIP_ON_DENORMAL);
  drawfilled2tetra(v1,v3,v2,v4,tetra0color,tetra1color,tetra2color,tetra3color,tetrabox_vis+6);

  glPopMatrix();
  // tetrahedron

  setClipPlanes(NULL,CLIP_OFF);
  glDisable(GL_LIGHTING);
  glDisable(GL_COLOR_MATERIAL);
}

/* ------------------ draw_test_outline ------------------------ */

void draw_test_outline(void){
  float *xmin, *xmax, *ymin, *ymax, *zmin, *zmax;
  unsigned char cubecolor[4]={0,0,0,255};
  unsigned char tetracoloroutline[4]={0,0,0,255};
  float *v1, *v2, *v3, *v4;
  float areas[6],cent_solid[3];
  int nverts;
  int faces[600], npolys, nfaces;
  int which_poly[200];
  float verts[600];

  box_state = b_state+1;
  v1 = tetra_vertices;
  v2 = v1 + 3;
  v3 = v2 + 3;
  v4 = v3 + 3;

  xmin = box_bounds;
  xmax = box_bounds+1;
  ymin = box_bounds+2;
  ymax = box_bounds+3;
  zmin = box_bounds+4;
  zmax = box_bounds+5;

  // tetrahedron

  glPushMatrix();
  glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
  glTranslatef(-xbar0,-ybar0,-zbar0);

  if(show_tetratest_labels == 1){
    output3Text(foregroundcolor, v1[0] - EPS, v1[1] - EPS, v1[2] - EPS, "v1");
    output3Text(foregroundcolor, v2[0] + EPS, v2[1] - EPS, v2[2] - EPS, "v2");
    output3Text(foregroundcolor, v3[0], v3[1] + EPS, v3[2] - EPS, "v3");
    output3Text(foregroundcolor, v4[0], v4[1], v4[2] + EPS, "v4");
  }

  antialias(ON);
  glLineWidth(tetra_line_thickness);
  drawtetra_outline(v1,v2,v3,v4,tetracoloroutline);
  antialias(OFF);

  glPopMatrix();
  // tetrahedron

  {
    float vsolid;

    FORTgetverts(box_bounds, v1, v2, v3, v4, verts, &nverts, faces, face_id, which_poly, &nfaces, &npolys, b_state);
    FORTgettetravol(box_bounds,v1,v2,v3,v4,&vsolid,areas,cent_solid);
    PRINTF("volume=%f\n",vsolid);
    PRINTF("centroid (solid)=%f %f %f\n",cent_solid[0],cent_solid[1],cent_solid[2]);
    {
      float vcell,vgas,cent_cell[3],cent_gas[3],*cc,*cs,*cg;

      cs=cent_solid;
      cg=cent_gas;
      cc=cent_cell;
      cc[0]=(box_bounds[0]+box_bounds[1])/2.0;
      cc[1]=(box_bounds[2]+box_bounds[3])/2.0;
      cc[2]=(box_bounds[4]+box_bounds[5])/2.0;
      vcell = box_bounds[1]-box_bounds[0];
      vcell *= (box_bounds[3]-box_bounds[2]);
      vcell *= (box_bounds[5]-box_bounds[4]);
      vcell = ABS(vcell);
      vgas = vcell - vsolid;
      if(vgas>0.0){
        cg[0] = (cc[0]*vcell-cs[0]*vsolid)/vgas;
        cg[1] = (cc[1]*vcell-cs[1]*vsolid)/vgas;
        cg[2] = (cc[2]*vcell-cs[2]*vsolid)/vgas;
      }
      else{
        cg[0]=-1.0;
        cg[0]=-1.0;
        cg[0]=-1.0;
      }
      PRINTF("centroid (gas)  =%f %f %f\n",cg[0],cg[1],cg[2]);
    }
    if(npolys>10){
      PRINTF("***error: nface=%i should not be bigger than 10\n",npolys);
    }
    if(nverts>0){
      int j;

      glPushMatrix();
      glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
      glTranslatef(-xbar0,-ybar0,-zbar0);
      glPointSize(tetra_point_size);
      glBegin(GL_POINTS);
      if(show_test_in_tetra==1){
        float green[4]={0.0,1.0,0.0,1.0};
        int in_tetra, tetra_state[4];

        glColor3fv(green);
        glVertex3fv(tetra_xyz);
        FORTtest_in_tetra(tetra_xyz,&in_tetra,tetra_state);
        PRINTF("in tetra:%i tetra state: %i %i %i %i\n",in_tetra,tetra_state[0],tetra_state[1],tetra_state[2],tetra_state[3]);
      }
      glColor3fv(foregroundcolor);
      for(j=0;j<nfaces;j++){
        if(tetrabox_vis[which_poly[j]]==1){
          glVertex3fv(verts+3*faces[3*j]);
          glVertex3fv(verts+3*faces[3*j+1]);
          glVertex3fv(verts+3*faces[3*j+2]);
        }
      }
      glEnd();
      glPopMatrix();
    }
  }

  glPushMatrix();
  glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
  glTranslatef(-xbar0,-ybar0,-zbar0);

  glPushMatrix();
  glTranslatef(*xmin,*ymin,*zmin);
  glScalef(ABS(*xmax-*xmin),ABS(*ymax-*ymin),ABS(*zmax-*zmin));
  if(show_tetratest_labels==1){
    char label[30];

    sprintf(label,"xmin area=%f",areas[0]);
    trimzeros(label);
    output3Text(foregroundcolor, -EPS, 0.5, 0.5, label);

    sprintf(label,"xmax area=%f",areas[1]);
    trimzeros(label);
    output3Text(foregroundcolor, 1.0+EPS, 0.5, 0.5, label);

    sprintf(label,"ymin area=%f",areas[2]);
    trimzeros(label);
    output3Text(foregroundcolor, 0.5, -EPS, 0.5, label);

    sprintf(label,"ymax area=%f",areas[3]);
    trimzeros(label);
    output3Text(foregroundcolor, 0.5, 1.0+EPS, 0.5, label);

    sprintf(label,"zmin area=%f",areas[4]);
    trimzeros(label);
    output3Text(foregroundcolor, 0.5, 0.5, -EPS, label);

    sprintf(label,"zmax area=%f",areas[5]);
    trimzeros(label);
    output3Text(foregroundcolor, 0.5, 0.5, 1.0+EPS, label);
  }

  antialias(ON);
  glLineWidth(tetra_line_thickness);
  drawcubec_outline(1.0,cubecolor);
  antialias(OFF);

  glPopMatrix();
  glPopMatrix();

}

/* ------------------ draw_geom_cutcells ------------------------ */

void draw_geom_cutcells(void){
  int i;

  glPushMatrix();
  glScalef(SCALE2SMV(1.0),SCALE2SMV(1.0),SCALE2SMV(1.0));
  glTranslatef(-xbar0,-ybar0,-zbar0);
  for(i=0;i<nmeshes;i++){
    meshdata *meshi;
    int j;
    int nx, nxy;
    float *x, *y, *z;

    meshi = meshinfo + i;
    nx = meshi->ibar;
    nxy = meshi->ibar*meshi->jbar;
    x = meshi->xplt_orig;
    y = meshi->yplt_orig;
    z = meshi->zplt_orig;

    if(meshi->ncutcells==0)continue;
    for(j=0;j<meshi->ncutcells;j++){
      int ijk, ii, jj, kk;

      ijk = meshi->cutcells[j];
      kk = ijk/nxy;
      jj = (ijk-kk*nxy)/nx;
      ii = ijk%nx;
      drawbox_outline(x[ii],x[ii+1],y[jj],y[jj+1],z[kk],z[kk+1],foregroundcolor);
    }
  }
  glPopMatrix();
}

/* ------------------ draw_test_triangle ------------------------ */

void draw_test_triangle(void){
  unsigned char trianglecolor[4] = {0, 0, 255, 255};
  unsigned char incolor[4] = {0, 255, 0, 255};
  unsigned char outcolor[4] = {255, 0, 0, 255};
  float *v1, *v2, *v3, *v4;
  int flag,flag2;

  v1 = tetra_vertices;
  v2 = v1 + 3;
  v3 = v2 + 3;
  v4 = v3 + 3;

  glPushMatrix();
  glScalef(SCALE2SMV(1.0), SCALE2SMV(1.0), SCALE2SMV(1.0));
  glTranslatef(-xbar0, -ybar0, -zbar0);

  antialias(ON);

  FORTget_in_triangle(v4, v1, v2, v3, &flag);
  FORTget_is_angle_ge_180(v1, v2, v3, &flag2);

  glLineWidth(tetra_line_thickness);
  glBegin(GL_LINES);
  if(flag2==1){
    glColor3ubv(outcolor);
  }
  else{
    glColor3ubv(incolor);
  }
  glVertex3f(v1[0],v1[1],0.0);
  glVertex3f(v2[0],v2[1],0.0);

  glVertex3f(v2[0],v2[1],0.0);
  glVertex3f(v3[0],v3[1],0.0);

  glColor3ubv(trianglecolor);
  glVertex3f(v3[0],v3[1],0.0);
  glVertex3f(v1[0],v1[1],0.0);
  glEnd();

  glPointSize(tetra_point_size);
  glBegin(GL_POINTS);
  if(flag==1){
    glColor3ubv(incolor);
  }
  else{
    glColor3ubv(outcolor);
  }
  glVertex3f(v4[0],v4[1],0.0);
  glEnd();
  antialias(OFF);
  output3Text(foregroundcolor, v1[0], v1[1], 0.0, "1");
  output3Text(foregroundcolor, v2[0], v2[1], 0.0, "2");
  output3Text(foregroundcolor, v3[0], v3[1], 0.0, "3");
  output3Text(foregroundcolor, v4[0], v4[1], 0.0, "4");

  glPopMatrix();
}

/* ------------------ draw_test_polygon ------------------------ */

void draw_test_polygon(void){
  float *v1, *v2, *v3, *v4;
  float verts[8];
  int nverts,poly[4], npoly, tris[12], ntris;
  int i;

  v1 = tetra_vertices;
  v2 = v1 + 3;
  v3 = v2 + 3;
  v4 = v3 + 3;

  glPushMatrix();
  glScalef(SCALE2SMV(1.0), SCALE2SMV(1.0), SCALE2SMV(1.0));
  glTranslatef(-xbar0, -ybar0, -zbar0);

  antialias(ON);

  verts[0] = v1[0];
  verts[1] = v1[1];
  verts[2] = v2[0];
  verts[3] = v2[1];
  verts[4] = v3[0];
  verts[5] = v3[1];
  verts[6] = v4[0];
  verts[7] = v4[1];
  nverts = 4;
  poly[0] = 1;
  poly[1] = 2;
  poly[2] = 3;
  poly[3] = 4;
  npoly = 4;


  FORTfpoly2tri(verts, &nverts, poly, &npoly, tris, &ntris);
  printf("triangles:\n");
  for(i = 0; i < ntris; i++){
    printf("%i: %i %i %i\n",i+1, tris[3 * i], tris[3 * i + 1], tris[3 * i + 2]);
  }
  printf("\n");

  glLineWidth(tetra_line_thickness);
  glBegin(GL_LINES);
  glColor3fv(foregroundcolor);
  for(i = 0; i < npoly; i++){
    int ii,iip1;

    ii = poly[i]-1;
    iip1 = poly[(i + 1) % npoly] - 1;
    glVertex3f(verts[2 * ii], verts[2 * ii + 1], 0.0);
    glVertex3f(verts[2 * iip1], verts[2 * iip1 + 1], 0.0);
  }
  glEnd();
  antialias(OFF);
  output3Text(foregroundcolor, v1[0], v1[1], 0.0, "1");
  output3Text(foregroundcolor, v2[0], v2[1], 0.0, "2");
  output3Text(foregroundcolor, v3[0], v3[1], 0.0, "3");
  output3Text(foregroundcolor, v4[0], v4[1], 0.0, "4");

  glPopMatrix();
}

/* ------------------ draw_geomdata ------------------------ */

void draw_geomdata(int flag, patchdata *patchi, int geom_type){
  int i;
  unsigned char *ivals;

  if(geom_type==GEOM_STATIC){
    ivals = patchi->geom_ival_static;
  }
  else{
    ivals = patchi->geom_ival_dynamic;
  }
  if(show_patch_solid == 1){
    for(i = 0; i < 1; i++){
      geomdata *geomi;
      geomlistdata *geomlisti;
      int ntris;
      int j;
      float *color;

      geomi = patchi->geominfo;
      if(geomi == NULL || geomi->display == 0 || geomi->loaded == 0)continue;
      if(geom_type == GEOM_STATIC){
        geomlisti = geomi->geomlistinfo - 1;
      }
      else{
        geomlisti = geomi->geomlistinfo + geomi->itime;
      }

      ntris = geomlisti->ntriangles;
      if(ntris == 0)continue;

      if(flag == DRAW_TRANSPARENT&&use_transparency_data == 1 && patchi->slice == 1)transparenton();

      glEnable(GL_NORMALIZE);
      glShadeModel(GL_SMOOTH);
      if(patchi->slice == 0)glEnable(GL_LIGHTING);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, iso_specular);
      glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, iso_shininess);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, block_ambient2);
      glEnable(GL_COLOR_MATERIAL);

      glPushMatrix();
      glScalef(SCALE2SMV(1.0), SCALE2SMV(1.0), SCALE2SMV(1.0));
      glTranslatef(-xbar0, -ybar0, -zbar0);
      glBegin(GL_TRIANGLES);
      if(smooth_iso_normal == 0){
        for(j = 0; j < ntris; j++){
          float *xyzptr[3];
          float *xyznorm;
          tridata *trianglei;
          int color_index;

          trianglei = geomlisti->triangles + j;

          xyznorm = trianglei->tri_norm;
          glNormal3fv(xyznorm);

          color_index = ivals[j];
          color = rgb_patch + 4 * color_index;
          if(patchi->slice == 1){
            if(trianglei->insolid == IN_CUTCELL && show_patch_incutcell == 0)continue;
            if(trianglei->insolid == IN_SOLID && show_patch_insolid == 0)continue;
            if(trianglei->insolid==IN_GAS&&show_patch_ingas==0)continue;
            glColor4f(color[0], color[1], color[2], transparent_level);
          }
          else{
            glColor3fv(color);
          }

          xyzptr[0] = trianglei->verts[0]->xyz;
          xyzptr[1] = trianglei->verts[1]->xyz;
          xyzptr[2] = trianglei->verts[2]->xyz;

          glVertex3fv(xyzptr[0]);
          glVertex3fv(xyzptr[1]);
          glVertex3fv(xyzptr[2]);

          if(patchi->slice == 1){
            glVertex3fv(xyzptr[0]);
            glVertex3fv(xyzptr[2]);
            glVertex3fv(xyzptr[1]);
          }
        }
      }
      else{
        for(j = 0; j < ntris; j++){
          float *xyzptr[3];
          float *xyznorm[3];
          tridata *trianglei;
          int color_index;

          trianglei = geomlisti->triangles + j;

          color_index = ivals[j];
          color = rgb_patch + 4 * color_index;
          if(patchi->slice == 1){
            if(trianglei->insolid==IN_CUTCELL&&show_patch_incutcell==0)continue;
            if(trianglei->insolid == IN_SOLID && show_patch_insolid == 0)continue;
            if(trianglei->insolid == IN_GAS && show_patch_ingas == 0)continue;
            glColor4f(color[0], color[1], color[2], transparent_level);
          }
          else{
            glColor3fv(color);
          }

          xyzptr[0] = trianglei->verts[0]->xyz;
          xyzptr[1] = trianglei->verts[1]->xyz;
          xyzptr[2] = trianglei->verts[2]->xyz;

          xyznorm[0] = trianglei->verts[0]->vert_norm;
          xyznorm[1] = trianglei->verts[1]->vert_norm;
          xyznorm[2] = trianglei->verts[2]->vert_norm;

          glNormal3fv(xyznorm[0]);
          glVertex3fv(xyzptr[0]);

          glNormal3fv(xyznorm[1]);
          glVertex3fv(xyzptr[1]);

          glNormal3fv(xyznorm[2]);
          glVertex3fv(xyzptr[2]);

          if(patchi->slice == 1){
            glNormal3fv(xyznorm[0]);
            glVertex3fv(xyzptr[0]);

            glNormal3fv(xyznorm[1]);
            glVertex3fv(xyzptr[2]);

            glNormal3fv(xyznorm[2]);
            glVertex3fv(xyzptr[1]);
          }
        }
      }
      glEnd();
      glPopMatrix();
      glDisable(GL_COLOR_MATERIAL);
      if(patchi->slice == 0)glDisable(GL_LIGHTING);
      if(flag == DRAW_TRANSPARENT&&use_transparency_data == 1 && patchi->slice == 1)transparentoff();
    }
  }
  if(show_patch_outline == 1){
    for(i = 0; i < 1; i++){
      geomdata *geomi;
      geomlistdata *geomlisti;
      int ntris;
      int j;
      float *color;

      geomi = patchi->geominfo;
      if(geomi == NULL || geomi->display == 0 || geomi->loaded == 0)continue;
      if(geom_type == GEOM_STATIC){
        geomlisti = geomi->geomlistinfo - 1;
      }
      else{
        geomlisti = geomi->geomlistinfo + geomi->itime;
      }

      ntris = geomlisti->ntriangles;
      if(ntris == 0)continue;

      glPushMatrix();
      glScalef(SCALE2SMV(1.0), SCALE2SMV(1.0), SCALE2SMV(1.0));
      glTranslatef(-xbar0, -ybar0, -zbar0);
      glBegin(GL_LINES);
        for(j = 0; j < ntris; j++){
          float *xyzptr[3];
          tridata *trianglei;
          int color_index;

          trianglei = geomlisti->triangles + j;
          if(patchi->slice==1){
            if(trianglei->insolid == IN_CUTCELL && show_patch_incutcell == 0)continue;
            if(trianglei->insolid == IN_SOLID && show_patch_insolid == 0)continue;
            if(trianglei->insolid==IN_GAS&&show_patch_ingas==0)continue;
          }

          color_index = ivals[j];
          color = rgb_patch + 4 * color_index;
          if(show_patch_solid == 1){
            glColor4fv(foregroundcolor);
          }
          else{
            glColor3fv(color);
          }

          xyzptr[0] = trianglei->verts[0]->xyz;
          xyzptr[1] = trianglei->verts[1]->xyz;
          xyzptr[2] = trianglei->verts[2]->xyz;

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
  }
  if(show_patch_verts == 1){
    for(i = 0; i < 1; i++){
      geomdata *geomi;
      geomlistdata *geomlisti;
      int ntris;
      int j;
      float *color;

      geomi = patchi->geominfo;
      if(geomi == NULL || geomi->display == 0 || geomi->loaded == 0)continue;
      if(geom_type == GEOM_STATIC){
        geomlisti = geomi->geomlistinfo - 1;
      }
      else{
        geomlisti = geomi->geomlistinfo + geomi->itime;
      }

      ntris = geomlisti->ntriangles;
      if(ntris == 0)continue;

      glPushMatrix();
      glScalef(SCALE2SMV(1.0), SCALE2SMV(1.0), SCALE2SMV(1.0));
      glTranslatef(-xbar0, -ybar0, -zbar0);
      glBegin(GL_POINTS);
      for(j = 0; j < ntris; j++){
        float *xyzptr[3];
        tridata *trianglei;

        trianglei = geomlisti->triangles + j;

        if(patchi->slice==1){
          if(trianglei->insolid == IN_CUTCELL && show_patch_incutcell == 0)continue;
          if(trianglei->insolid == IN_SOLID && show_patch_insolid == 0)continue;
          if(trianglei->insolid==IN_GAS&&show_patch_ingas==0)continue;
        }
        if(show_patch_solid == 1||show_patch_outline==1){
          glColor4fv(foregroundcolor);
        }
        else{
          int color_index;

          color_index = ivals[j];
          color = rgb_patch + 4 * color_index;
          glColor3fv(color);
        }

        xyzptr[0] = trianglei->verts[0]->xyz;
        xyzptr[1] = trianglei->verts[1]->xyz;
        xyzptr[2] = trianglei->verts[2]->xyz;

        glVertex3fv(xyzptr[0]);
        glVertex3fv(xyzptr[1]);
        glVertex3fv(xyzptr[2]);
      }
      glEnd();
      glPopMatrix();
    }
  }

}

/* ------------------ CompareTransparentTriangles ------------------------ */

int CompareTransparentTriangles(const void *arg1, const void *arg2){
  tridata *tri, *trj;

  tri = *(tridata **)arg1;
  trj = *(tridata **)arg2;

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

/* ------------------ ShowHideSortGeometry ------------------------ */

void ShowHideSortGeometry(float *mm){
  int i;
  int count_transparent,count_opaque;
  int itime;
  int *showlevels=NULL;
  int iter;

  if(loaded_isomesh!=NULL)showlevels=loaded_isomesh->showlevels;

  for(iter = 0; iter < 2; iter++){
    CheckMemory;
    count_transparent = 0;
    count_opaque = 0;
    ntransparent_triangles = count_transparent;
    nopaque_triangles = count_opaque;
    for(i = 0; i < ngeominfoptrs; i++){
      geomlistdata *geomlisti;
      int j;
      geomdata *geomi;

      geomi = geominfoptrs[i];
      if( (geomi->fdsblock == NOT_FDSBLOCK && geomi->geomtype!=GEOM_ISO)|| geomi->patchactive == 1)continue;
      for(itime = 0; itime < 2; itime++){
        if(itime == 0){
          geomlisti = geomi->geomlistinfo - 1;
        }
        else{
          geomlisti = geomi->geomlistinfo + geomi->itime;
          if(geomi->currentframe != NULL)geomlisti = geomi->currentframe;
        }

        for(j = 0; j < geomlisti->ntriangles; j++){
          tridata *tri;
          float xyz[3];
          float *xyz1, *xyz2, *xyz3;
          float xyzeye[3];
          int isurf;
          int is_opaque;

          is_opaque = 0;
          tri = geomlisti->triangles + j;
          if(hilight_skinny == 1 && tri->skinny == 1)is_opaque = 1;
          if(tri->surf->transparent_level >= 1.0)is_opaque = 1;
          isurf = tri->surf - surfinfo - nsurfinfo - 1;
          if((geomi->geomtype==GEOM_ISO&&showlevels != NULL&&showlevels[isurf] == 0) || tri->surf->transparent_level <= 0.0){
            continue;
          }
          if(iter == 1){
             tri->geomtype = geomi->geomtype;
             tri->verts[0]->geomtype = geomi->geomtype;
             tri->verts[1]->geomtype = geomi->geomtype;
             tri->verts[2]->geomtype = geomi->geomtype;
          }
          if(is_opaque == 1){
            if(iter==1)opaque_triangles[count_opaque] = tri;
            count_opaque++;
            if(iter==0)continue;
          }
          else{
            if(iter==1)transparent_triangles[count_transparent] = tri;
            count_transparent++;
          }
          if(iter==0&&sort_geometry == 1){
            xyz1 = tri->verts[0]->xyz;
            xyz2 = tri->verts[1]->xyz;
            xyz3 = tri->verts[2]->xyz;
            xyz[0] = NORMALIZE_X((xyz1[0] + xyz2[0] + xyz3[0]) / 3.0);
            xyz[1] = NORMALIZE_Y((xyz1[1] + xyz2[1] + xyz3[1]) / 3.0);
            xyz[2] = NORMALIZE_Z((xyz1[2] + xyz2[2] + xyz3[2]) / 3.0);

            xyzeye[0] = mm[0] * xyz[0] + mm[4] * xyz[1] + mm[8] * xyz[2] + mm[12];
            xyzeye[1] = mm[1] * xyz[0] + mm[5] * xyz[1] + mm[9] * xyz[2] + mm[13];
            xyzeye[2] = mm[2] * xyz[0] + mm[6] * xyz[1] + mm[10] * xyz[2] + mm[14];
            xyzeye[0] /= mscale[0];
            xyzeye[1] /= mscale[1];
            xyzeye[2] /= mscale[2];
            tri->distance = xyzeye[0] * xyzeye[0] + xyzeye[1] * xyzeye[1] + xyzeye[2] * xyzeye[2];
            CheckMemory;
          }
        }
      }
    }
    if(iter == 0){
      CheckMemory;
      if(count_transparent == 0 && count_opaque == 0)return;
      FREEMEMORY(alltriangles);
      NewMemory((void **)&alltriangles, (count_opaque + count_transparent)*sizeof(tridata **));
      transparent_triangles = alltriangles;
      opaque_triangles = alltriangles + count_transparent;
    }
  }
  ntransparent_triangles = count_transparent;
  nopaque_triangles = count_opaque;
  if(sort_geometry==1&&ntransparent_triangles>0){
    qsort((isotri **)transparent_triangles, (size_t)ntransparent_triangles, sizeof(tridata **), CompareTransparentTriangles);
  }
}

/* ------------------ init_geom ------------------------ */

void init_geom(geomdata *geomi,int geomtype, int fdsblock){
  geomi->file=NULL;
  geomi->display=0;
  geomi->loaded=0;
  geomi->geomlistinfo_0=NULL;
  geomi->surf=NULL;
  geomi->geomlistinfo=NULL;
  geomi->currentframe = NULL;
  geomi->times=NULL;
  geomi->ntimes=0;
  geomi->times=NULL;
  geomi->timeslist=NULL;
  geomi->float_vals=NULL;
  geomi->int_vals=NULL;
  geomi->nfloat_vals=0;
  geomi->nint_vals=0;
  geomi->geomtype = geomtype;
  geomi->fdsblock = fdsblock;
}
/* ------------------ rotateu2v ------------------------ */

void rotateu2v(float *u, float *v, float *axis, float *angle) {
  float sum, cosangle, normu, normv;

  /*
  i  j  k
  ux uy uz
  vx vy vz
  */

  CROSS(axis, u, v);
  sum = NORM3(axis);
  normu = NORM3(u);
  normv = NORM3(v);
  if (sum>0.0&&normu>0.0&&normv>0.0) {
    axis[0] /= sum;
    axis[1] /= sum;
    axis[2] /= sum;
    cosangle = CLAMP(DOT3(u, v) / (normu*normv), -1.0, 1.0);
    *angle = acos(cosangle);
  }
  else {
    axis[0] = 0.0;
    axis[1] = 0.0;
    axis[2] = 1.0;
    *angle = 0.0;
  }
}

/* ------------------ angleaxis2quat ------------------------ */

void angleaxis2quat(float angle, float *axis, float *quat) {
  float sum;
  float cosang, sinang;

  // angle is in radians
  // axis is a vector

  sum = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);

  if (sum>0.0) {
    cosang = cos(angle / 2.0);
    sinang = sin(angle / 2.0);

    quat[0] = cosang;
    quat[1] = axis[0] * sinang / sum;
    quat[2] = axis[1] * sinang / sum;
    quat[3] = axis[2] * sinang / sum;
  }
  else {
    quat[0] = 1.0;
    quat[1] = 0.0;
    quat[2] = 0.0;
    quat[3] = 0.0;
  }
}

/* ------------------ quat2rot------------------ */

void quat2rot(float quat[4], float rot[16]) {
  float w, x, y, z, sum;

  sum = sqrt(quat[0] * quat[0] + quat[1] * quat[1] + quat[2] * quat[2] + quat[3] * quat[3]);
  w = quat[0] / sum;
  x = quat[1] / sum;
  y = quat[2] / sum;
  z = quat[3] / sum;

  rot[0] = 1.0 - 2.0*y*y - 2.0*z*z;
  rot[1] = 2.0*x*y + 2.0*w*z;
  rot[2] = 2.0*x*z - 2.0*w*y;
  rot[3] = 0.0;

  rot[4] = 2.0*x*y - 2.0*w*z;
  rot[5] = 1.0 - 2.0*x*x - 2.0*z*z;
  rot[6] = 2.0*y*z + 2.0*w*x;
  rot[7] = 0.0;

  rot[8] = 2.0*x*z + 2.0*w*y;
  rot[9] = 2.0*y*z - 2.0*w*x;
  rot[10] = 1.0 - 2.0*x*x - 2.0*y*y;
  rot[11] = 0.0;

  rot[12] = 0.0;
  rot[13] = 0.0;
  rot[14] = 0.0;
  rot[15] = 1.0;
}

/* ------------------ mult_quat ------------------------ */

void mult_quat(float x[4], float y[4], float z[4]) {
  float z2[4];

  z2[0] = x[0] * y[0] - x[1] * y[1] - x[2] * y[2] - x[3] * y[3];
  z2[1] = x[0] * y[1] + x[1] * y[0] + x[2] * y[3] - x[3] * y[2];
  z2[2] = x[0] * y[2] - x[1] * y[3] + x[2] * y[0] + x[3] * y[1];
  z2[3] = x[0] * y[3] + x[1] * y[2] - x[2] * y[1] + x[3] * y[0];
  z[0] = z2[0];
  z[1] = z2[1];
  z[2] = z2[2];
  z[3] = z2[3];
}

/* ------------------ normalize_quat ------------------------ */

void normalize_quat(float x[4]) {
  float sum;

  sum = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
  if (sum>0.0) {
    x[0] /= sum;
    x[1] /= sum;
    x[2] /= sum;
    x[3] /= sum;
  }
}

/* ------------------ xyz2azelev ------------------------ */

void xyz2azelev(float *xyz, float *azimuth, float *elevation) {
  float norm3;

  // x = ||xyz||cos(az)*cos(elev)
  // y = ||xyz||sin(az)*cos(elev)
  // z = ||xyz||sin(elev)
  // elev=asin(z/||xyz||)
  // az=atan(y/x)
  norm3 = NORM3(xyz);
  if (norm3>0.00001&&ABS(xyz[2] / norm3) <= 1.0) {
    *elevation = RAD2DEG*asin(xyz[2] / norm3);
  }
  else {
    *elevation = 0.0;
  }
  *azimuth = RAD2DEG*atan2(xyz[1], xyz[0]);
}


