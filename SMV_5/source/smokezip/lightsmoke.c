// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#ifdef pp_LIGHT
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "zlib.h"
#include "egz_stdio.h"
#include "svzip.h"
#include "MALLOC.h"

// svn revision character string
char lightsmoke_revision[]="$Revision$";

float average(int ijknode, float *celldata, mesh *smoke_mesh);
float get_photon_step(smoke3d *smoke3di, float xyzpos[3], float xyzdir[3]);
int in_mesh(mesh *smoke_mesh,float xyzpos[3]);
float getlog_1_m_alpha(mesh *smoke_mesh, float *xyz);
int get_interval(float val, float *array, int n);
float interp(float xyz[3], mesh *smoke_mesh, float *full_logalphabuffer);

#define IJKNODE(i,j,k) ((i)+(j)*nx+(k)*nxy)
#define IJKCELL(i,j,k) ((i)+(j)*nxcell+(k)*nxycell)

/* ------------------ get_random_light ------------------------ */

lightdata *get_random_light(void){
  lightdata *light;
  int ilight;
  float val;

  if(nlightinfo==1)return lightinfo;
  val = rand_1d(0.0,1.0);
  ilight = get_interval(val, light_cdf, nlightinfo+1);
  light = lightinfo + ilight;
  return light;
}

/* ------------------ integrate_alpha ------------------------ */

float integrate_alpha(mesh *smoke_mesh,float *xyz_light,float *xyz2){
  float dx, dy, dz, dist, dint, xyz[3];
  int i, n;
  float f1, sum;
  float returnval;

  dx = xyz2[0]-xyz_light[0];
  dy = xyz2[1]-xyz_light[1];
  dz = xyz2[2]-xyz_light[2];
  dist = dx*dx + dy*dy + dz*dz;
  dist = sqrt(dist);
  n = dist/smoke_mesh->dx + 2;
  dint = dist/n;

  sum = 0.0;
  for(i=0;i<n+1;i++){
    float val;

    f1 = (float)(n-i)/(float)n;
    xyz[0] = xyz_light[0]*f1 + xyz2[0]*(1.0-f1);
    xyz[1] = xyz_light[1]*f1 + xyz2[1]*(1.0-f1);
    xyz[2] = xyz_light[2]*f1 + xyz2[2]*(1.0-f1);
    val = getlog_1_m_alpha(smoke_mesh, xyz);
    if(i==0||i==n-1){
      sum += val/2.0;
    }
    else{
      sum += val;
    }
    if(sum<-5.0)return 0.0;
  }
  sum *= dint;
  returnval = exp(sum);
  if(returnval>1.0)returnval=1.0;
  if(returnval<0.0)returnval=0.0;
  return returnval;
}

/* ------------------ remove_trial ------------------------ */

void remove_trial(nodedata *node){
  nodedata *prev, *next;

  next=node->next;
  prev=node->prev;

  prev->next=next;
  next->prev=prev;
}

/* ------------------ insert_trial ------------------------ */

void insert_trial(nodedata *node,nodedata *before){
  nodedata *prev, *next;

  next=before;
  prev=before->prev;

  prev->next=node;
  node->prev=prev;

  node->next=next;
  next->prev=node;
}

/* ------------------ update_lightfield ------------------------ */

void update_lightfield(float time, smoke3d *smoke3di, unsigned char *lightingbuffer){
    int i, j, k, ijk, ii;
    mesh *smoke_mesh;
    int nxcell, nycell, nzcell;
    int nx, ny, nxy, nz;
    nodedata *nodeinfo;
    nodedata *first_trial, *last_trial, first_t, last_t;
    int naccepted;
    

    printf("time=%f\n",time);
    CheckMemory;
    smoke_mesh=smoke3di->smoke_mesh;

    nxcell = smoke_mesh->ibar;
    nycell = smoke_mesh->jbar;
    nzcell = smoke_mesh->kbar;
    nx=nxcell+1;
    ny=nycell+1;
    nz=nzcell+1;
    nxy=nx*ny;

    nodeinfo=smoke3di->nodeinfo;
    if(nodeinfo==NULL){
    NewMemory((void **)&nodeinfo, nx*ny*nz*sizeof(nodedata));
    smoke3di->nodeinfo=nodeinfo;

    // define node data structures

    for(k=0;k<nz;k++){
      for(j=0;j<ny;j++){
        for(i=0;i<nx;i++){
          nodedata *nodei;

          ijk = IJKNODE(i,j,k);
          nodei = nodeinfo + ijk;

          nodei->state=UNKNOWN;
          nodei->next=NULL;
          nodei->prev=NULL;

          nodei->nabors[0]=nodeinfo+ijk-1;
          nodei->nabors[1]=nodeinfo+ijk+1;
          nodei->nabors[2]=nodeinfo+ijk-nx;
          nodei->nabors[3]=nodeinfo+ijk+nx;
          nodei->nabors[4]=nodeinfo+ijk-nxy;
          nodei->nabors[5]=nodeinfo+ijk+nxy;

          for(ii=0;ii<6;ii++){
            if(nodei->nabors[ii]-nodeinfo<0||nodei->nabors[ii]-nodeinfo>nx*ny*nz-1){
              nodei->nabors[ii]=NULL;
            }
          }
        }
      }
    }
    }
    naccepted=0;


    for(i=0;i<nx*ny*nz;i++){
      nodedata *nodei;

      nodei = nodeinfo + i;
      nodei->totaldist=255.0;
      nodei->nodedist=(float)(255-full_alphabuffer[i])/255.0;
      nodei->state=UNKNOWN;
    }
    
    // assume boundaries light comes in through open boundaries

    for(ii=0;ii<smoke3di->smoke_mesh->nvents;ii++){
      ventdata *venti;
      int *ib;
      int jjj,kkk;

      venti = smoke3di->smoke_mesh->ventinfo + ii;
      if(venti->is_open==0)continue;

      ib = venti->ib;
      kkk=(ib[4]-1)*nxy;
      for(k=ib[4];k<=ib[5];k++){
        kkk+=nxy;
        jjj=(ib[2]-1)*nx;
      for(j=ib[2];j<=ib[3];j++){
        jjj+=nx;
      for(i=ib[0];i<=ib[1];i++){
        nodedata *nodei;

        ijk=kkk+jjj+i;
        nodei = nodeinfo + ijk;
        nodei->state=ACCEPTED;
        naccepted++;
      }
      }
      }
    }

    // define trial list

    first_trial=&first_t;
    last_trial=&last_t;

    first_trial->prev=NULL;
    first_trial->next=last_trial;

    last_trial->prev=first_trial;
    last_trial->next=NULL;
    for(i=0;i<nx*ny*nz;i++){
      nodedata *nodei;

      nodei = nodeinfo+i;
      if(nodei->state!=UNKNOWN)continue;

      for(j=0;j<6;j++){
        nodedata *nodej;

        nodej=nodei->nabors[j];
        if(nodej==NULL||nodej->state!=ACCEPTED)continue;
        nodei->state=TRIAL;
        insert_trial(nodei,last_trial);
        nodei->totaldist=nodej->totaldist*(nodej->nodedist+nodei->nodedist)/2.0;
        break;
      }
    }

    // define trial linked list
    // also define initial distance for trial nodes


    // now it begins
    // 1.  find the trial node with the minimum distance from adjacent accepted nodes
    // 2.     remove it from the trial list
    // 3.     add all unknown neighbors of trial_x, to the trial list
    // 4.     recompute distance to all neighbors of trial list
    //        repeat 1-4 until there are no more trial  nodes

    while(first_trial!=last_trial->prev){
      nodedata *maxti;
//      float maxdist;

      maxti=first_trial->next;
  //    maxdist=maxti->totaldist;
      /*
      if(maxdist!=255.0){
        for(ti=maxti->next;ti!=last_trial;ti=ti->next){
          if(ti->totaldist>maxdist){
            maxdist=ti->totaldist;
            maxti=ti;
          }
        }
      }
      */
      remove_trial(maxti);
      maxti->state=ACCEPTED;
      for(j=0;j<6;j++){
        nodedata *tj;

        tj=maxti->nabors[j];
        if(tj!=NULL&&tj->state!=ACCEPTED){
          if(tj->state==UNKNOWN){
            tj->state=TRIAL;
            insert_trial(tj,last_trial);
          }
          tj->totaldist=maxti->totaldist*(tj->nodedist+maxti->nodedist)/2.0;
        }
      }
    }
    for(i=0;i<nx*ny*nz;i++){
      nodedata *ti;

      ti = nodeinfo + i;
      lightingbuffer[i]=albedo*(unsigned char)ti->totaldist;
    }
}

/* ------------------ update_lightfield ------------------------ */

void update_lightfield_simple(float time, smoke3d *smoke3di, unsigned char *lightingbuffer){
 // accumulate hrr for each light

    int i, j, k;
    mesh *smoke_mesh;
    int nxcell, nycell, nzcell;
    float *tau, *f_light;
    int nx, ny, nxy, nz;

    tau = full_logalphabuffer;

    CheckMemory;
    smoke_mesh=smoke3di->smoke_mesh;
    f_light=smoke_mesh->light_cell_radiance;

    nxcell = smoke_mesh->ibar;
    nycell = smoke_mesh->jbar;
    nzcell = smoke_mesh->kbar;
    nx=nxcell+1;
    ny=nycell+1;
    nz=nzcell+1;
    nxy=nx*ny;

    for(i=0;i<nx*ny*nz;i++){
      f_light[i]=255.0;
      tau[i] = (float)(255-full_alphabuffer[i])/255.0;
    }

    for(k=0;k<nz;k++){
      for(j=0;j<ny;j++){
        for(i=1;i<nx;i++){
          int ijk, ijk2;

          ijk=IJKNODE(i,j,k);
          ijk2=ijk-1;
          f_light[ijk]=f_light[ijk2]*tau[ijk2];
        }
      }
    }
    for(k=0;k<nz;k++){
      for(j=0;j<ny;j++){
        for(i=nx-2;i>=0;i--){
          int ijk, ijk2;
          float val;

          ijk=IJKNODE(i,j,k);

          ijk2=ijk+1;
          val=f_light[ijk2]*tau[ijk2];
          if(val<f_light[ijk])break;
          f_light[ijk]=val;
        }
      }
    }
    for(k=0;k<nz;k++){
      for(i=0;i<nx;i++){
        for(j=1;j<ny;j++){
          int ijk, ijk2;
          float val;

          ijk=IJKNODE(i,j,k);
          ijk2=ijk-nx;
          val=f_light[ijk2]*tau[ijk2];
          if(val<f_light[ijk])break;
          f_light[ijk]=val;
        }
      }
    }
    for(k=0;k<nz;k++){
      for(i=0;i<nx;i++){
        for(j=ny-2;j>=0;j--){
          int ijk, ijk2;
          float val;

          ijk=IJKNODE(i,j,k);
          ijk2=ijk+nx;
          val=f_light[ijk2]*tau[ijk2];
          if(val<f_light[ijk])break;
          f_light[ijk]=val;
        }
      }
    }
    for(i=0;i<nx*ny*nz;i++){
      lightingbuffer[i]=albedo*(unsigned char)f_light[i];
    }
    CheckMemory;
}

/* ------------------ update_lightfield ------------------------ */

void update_lightfield_lights(float time, smoke3d *smoke3di, unsigned char *lightingbuffer){
 // accumulate hrr for each light

    lightdata *lighti;
    int ilight;
    int i, j, k;
    mesh *smoke_mesh;
    int nx, ny, nz;
    int nxcell, nycell, nzcell;
    int nxycell;
    float *xplt, *yplt, *zplt;
    float xyz_cell_pos[3];
    float logindex[256];
    float *light_cell_radiance;
    float four_pi;
    float radiance_min, radiance_max;
    float radiance_mean, radiance_dev;
    int ijk_node, ijk_cell;

    // pick a random light weighted by HRR  
    //  (ie a high wattage light will be picked proportionately more often than
    //          a low wattage one)

    four_pi = 16.0*atan(1.0);
    CheckMemory;
    smoke_mesh=smoke3di->smoke_mesh;
    xplt = smoke_mesh->xplt;
    yplt = smoke_mesh->yplt;
    zplt = smoke_mesh->zplt;
    light_cell_radiance = smoke_mesh->light_cell_radiance;

    nxcell = smoke_mesh->ibar;
    nycell = smoke_mesh->jbar;
    nxycell = nxcell*nycell;
    nzcell = smoke_mesh->kbar;
    nx = nxcell+1;
    ny = nycell+1;
    nz = nzcell+1;

    for(i=1;i<256;i++){
      logindex[i]=log((float)i/255.0)/smoke_mesh->dx;
    }
    logindex[0]=log(0.001)/smoke_mesh->dx;

    for(i=0;i<nxcell*nycell*nzcell;i++){
      full_logalphabuffer[i] = logindex[255-full_alphabuffer[i]];
      lightingbuffer[i]=0;
      light_cell_radiance[i]=0.0;
    }

    ijk_cell = 0;
    for(k=0;k<nzcell;k++){
      xyz_cell_pos[2]=(zplt[k]+zplt[k+1])/2.0;
      for(j=0;j<nycell;j++){
        xyz_cell_pos[1]=(yplt[j]+yplt[j+1])/2.0;
        for(i=0;i<nxcell;i++){
          float val2;

          xyz_cell_pos[0]=(xplt[i]+xplt[i+1])/2.0;
          val2 = getlog_1_m_alpha(smoke_mesh, xyz_cell_pos);
          if(val2<0.0){
            for(ilight=0;ilight<nlightinfo;ilight++){
              float radiance, dx, dy, dz;
              float *xyz_e1, *xyz_e2;
              float *xyz_light;
              int ilight2;
              float light_dist2;

              lighti = lightinfo + ilight;

              // for now, lights do not light smoke in other meshes
              if(in_mesh(smoke_mesh,lighti->xyz1)==0)continue;

              switch (lighti->type){
                case 0:
                xyz_light=lighti->xyz1;
                dx = xyz_light[0]-xyz_cell_pos[0];
                dy = xyz_light[1]-xyz_cell_pos[1];
                dz = xyz_light[2]-xyz_cell_pos[2];
                light_dist2 = dx*dx+dy*dy+dz*dz;
                if(light_dist2<0.1)light_dist2=0.1;
//
// qgas = qlight*optical_thickness*(cell_cross_section/(4*pi*lightdist))*albedo
// radiance = qgas/(cell_area*pi)
//

                radiance = lighti->q*integrate_alpha(smoke_mesh,xyz_light,xyz_cell_pos)*albedo/(four_pi*four_pi*light_dist2);
                light_cell_radiance[ijk_cell] += radiance;
                break;

                case 1:

                xyz_e1=lighti->xyz1;
                xyz_e2=lighti->xyz2;

                for(ilight2=0;ilight2<lighti->nstep;ilight2++){
                  float f1;

                  f1 = (float)(lighti->nstep-1-ilight2)/(float)(lighti->nstep-1);
                  xyz_light[0] = f1*xyz_e1[0] + (1.0-f1)*xyz_e2[0];
                  xyz_light[1] = f1*xyz_e1[1] + (1.0-f1)*xyz_e2[1];
                  xyz_light[2] = f1*xyz_e1[2] + (1.0-f1)*xyz_e2[2];

                  dx = xyz_light[0]-xyz_cell_pos[0];
                  dy = xyz_light[1]-xyz_cell_pos[1];
                  dz = xyz_light[2]-xyz_cell_pos[2];
                  light_dist2 = dx*dx+dy*dy+dz*dz;
                  if(light_dist2<0.1)light_dist2=0.1;

                  radiance = lighti->q*integrate_alpha(smoke_mesh,xyz_light,xyz_cell_pos)*albedo/(four_pi*four_pi*light_dist2);
                  light_cell_radiance[ijk_cell] += radiance;
                }
                break;
              }
            }
          }

          ijk_cell++;
        }
      }
    }
    radiance_min=light_cell_radiance[0];
    radiance_max=radiance_min;
    radiance_mean=0.0;
  //  factor = 1.0/log(light_max/light_min);
    for(i=0;i<nxcell*nycell*nzcell;i++){
      radiance_mean += light_cell_radiance[i];
      if(light_cell_radiance[i]<radiance_min)radiance_min=light_cell_radiance[i];
      if(light_cell_radiance[i]>radiance_max)radiance_max=light_cell_radiance[i];
    }
    radiance_mean /= (nxcell*nycell*nzcell);

    radiance_dev=0.0;
    for(i=0;i<nxcell*nycell*nzcell;i++){
      float dd;

      dd = light_cell_radiance[i] - radiance_mean;
      radiance_dev += dd*dd;
    }
    radiance_dev = sqrt(radiance_dev/(nxcell*nycell*nzcell));
    printf("radiance min=%f mean=%f dev=%f max=%f\n", radiance_min, radiance_mean, radiance_dev, radiance_max);

    ijk_node=0;
    for(k=0;k<nz;k++){
      int k1, k2;

      k1=k-1;
      if(k1<0)k1=0;
      k2=k;
      if(k2==nz-1)k2=k1;

      for(j=0;j<ny;j++){
        int j1, j2;

        j1=j-1;
        if(j1<0)j1=0;
        j2=j;
        if(j2==ny-1)j2=j1;

        for(i=0;i<nx;i++){
          int i1, i2;
          float val_111, val_112, val_121, val_122;
          float val_211, val_212, val_221, val_222;
          int cell_111, cell_112, cell_121, cell_122;
          int cell_211, cell_212, cell_221, cell_222;
          float radiance;

          i1=i-1;
          if(i1<0)i1=0;
          i2=i;
          if(i2==nx-1)i2=i1;

          //#define IJKCELL(i,j,k) ((i)+(j)*nxcell+(k)*nxycell)

          cell_111 = IJKCELL(i1,j1,k1);
          if(i2==i1){
            cell_211 = cell_111;
          }
          else{
            cell_211 = cell_111 + 1;
          }
          if(j2==j1){
            cell_121 = cell_111;
          }
          else{
            cell_121 = cell_111 + nxcell;
          }
          if(i2==i1){
            cell_221 = cell_121;
          }
          else{
            cell_221 = cell_121 + 1;
          }

          if(k2==k1){
            cell_112 = cell_111;
            cell_122 = cell_121;
            cell_212 = cell_211;
            cell_222 = cell_221;
          }
          else{
            cell_112 = cell_111 + nxycell;
            cell_122 = cell_121 + nxycell;
            cell_212 = cell_211 + nxycell;
            cell_222 = cell_221 + nxycell;
          }

          val_111 = light_cell_radiance[cell_111];
          val_112 = light_cell_radiance[cell_112];
          val_121 = light_cell_radiance[cell_121];
          val_122 = light_cell_radiance[cell_122];
          val_211 = light_cell_radiance[cell_211];
          val_212 = light_cell_radiance[cell_212];
          val_221 = light_cell_radiance[cell_221];
          val_222 = light_cell_radiance[cell_222];

          radiance  = (val_111+val_112+val_121+val_122);
          radiance += (val_211+val_212+val_221+val_222);
          radiance /=8.0;
          if(radiance<light_min)radiance=light_min;
          if(radiance>light_max)radiance=light_max;
//          lightingbuffer[ijk_node]=254*log(radiance/light_min)*factor;
          lightingbuffer[ijk_node]=254*(radiance-light_min)/(light_max-light_min);

          ijk_node++;

        }
      }
    }
  //  printf("fluxmin=%f fluxmax=%f\n",fluxmin,fluxmax);
    CheckMemory;
}

/* ------------------ update_lightfield ------------------------ */

void update_lightfield2(float time, smoke3d *smoke3di, unsigned char *lightingbuffer){
 // accumulate hrr for each light

    lightdata *lighti;
    float xyzpos[3], xyzdir[3];
    float factor, *xyz1, *xyz2;
    int i;
    float *photon_cell;
    float photon_step;
    mesh *smoke_mesh;
    int nx, ny, nz;
    int nxcell, nycell, nzcell, nxycell;
    float binmax,binsum;
    float logindex[256];

    // pick a random light weighted by HRR  
    //  (ie a high wattage light will be picked proportionately more often than
    //          a low wattage one)

    CheckMemory;
    smoke_mesh=smoke3di->smoke_mesh;
    nxcell = smoke_mesh->ibar;
    nycell = smoke_mesh->jbar;
    nzcell = smoke_mesh->kbar;
    nxycell = nxcell*nycell;
    nx = nxcell+1;
    ny = nycell+1;
    nz = nzcell+1;

// zero out bin used to collect photons

    photon_cell = smoke3di->smoke_mesh->photon_cell;
    for(i=1;i<256;i++){
      logindex[i]=log((float)i/255.0)/smoke_mesh->dx;
    }
    logindex[0]=log(0.001)/smoke_mesh->dx;

    for(i=0;i<nxcell*nycell*nzcell;i++){
      photon_cell[i]=0.0;
    }
    for(i=0;i<nx*ny*nz;i++){
      float val;
      int ival;

      ival = (int)full_alphabuffer[i];
      val = logindex[255-ival];
      full_logalphabuffer[i]=val;
    }
    CheckMemory;

    for(i=0;i<NPHOTONS;i++){
      lighti = get_random_light();

    //  get a random position within the light and a random direction

      switch (lighti->type){

        case 0:      // point
          if(lighti->move==0){
            xyzpos[0]=lighti->xyz1[0];
            xyzpos[1]=lighti->xyz1[1];
            xyzpos[2]=lighti->xyz1[2];
          }
          else{
            float f1, f2, dt;

            dt = lighti->t2 - lighti->t1;
            if(time<lighti->t1){
              f1 = 0.0;
            }
            else if(time>lighti->t2){
              f1 = 1.0;
            }
            else{
              if(dt>0.0){
                f1 = (time - lighti->t1)/dt;
              }
              else{
                f1 = 0.0;
              }
            }
            f2 = 1.0 - f1;
            xyzpos[0]=f2*lighti->xyz1[0]+f1*lighti->xyz2[0];
            xyzpos[1]=f2*lighti->xyz1[1]+f1*lighti->xyz2[1];
            xyzpos[2]=f2*lighti->xyz1[2]+f1*lighti->xyz2[2];
            if(i==0){
              printf("(%f, %f) (%f, %f, %f)\n",f1,f2,xyzpos[0],xyzpos[1],xyzpos[2]);
            }
          }
          rand_sphere_dir(xyzdir);
          break;

        case 1:      // line
          factor = rand_1d(0.0,1.0);
          xyz1 = lighti->xyz1;
          xyz2 = lighti->xyz2;
          xyzpos[0]= xyz1[0]*(1-factor) + xyz2[0]*factor;
          xyzpos[1]= xyz1[1]*(1-factor) + xyz2[1]*factor;
          xyzpos[2]= xyz1[2]*(1-factor) + xyz2[2]*factor;
          rand_sphere_dir(xyzdir);
          break;

        case 2:      // region
          xyz1 = lighti->xyz1;
          xyz2 = lighti->xyz2;
          factor = rand_1d(0.0,1.0);
          xyzpos[0]= xyz1[0]*(1-factor) + xyz2[0]*factor;
          factor = rand_1d(0.0,1.0);
          xyzpos[1]= xyz1[1]*(1-factor) + xyz2[1]*factor;
          factor = rand_1d(0.0,1.0);
          xyzpos[2]= xyz1[2]*(1-factor) + xyz2[2]*factor;
          rand_absdir(xyzdir,lighti->dir);
          break;
      }
      for(;;){
        int i1, j1, k1, ijk;

        photon_step=get_photon_step(smoke3di, xyzpos, xyzdir);
        CheckMemory;
        xyzpos[0]+= photon_step*xyzdir[0];
        xyzpos[1]+= photon_step*xyzdir[1];
        xyzpos[2]+= photon_step*xyzdir[2];
        if(in_mesh(smoke_mesh,xyzpos)==0)break;

//        i1 = get_interval(xyzpos[0],smoke_mesh->xplt,nxcell);
//        j1 = get_interval(xyzpos[1],smoke_mesh->yplt,nycell);
//        k1 = get_interval(xyzpos[2],smoke_mesh->zplt,nzcell);
        i1 = GET_INTERVAL(xyzpos[0],smoke_mesh->xplt[0],smoke_mesh->dx);
        i1 = BOUND(i1,0,nxcell-1);
        j1 = GET_INTERVAL(xyzpos[1],smoke_mesh->yplt[0],smoke_mesh->dy);
        j1 = BOUND(j1,0,nycell-1);
        k1 = GET_INTERVAL(xyzpos[2],smoke_mesh->zplt[0],smoke_mesh->dz);
        k1 = BOUND(k1,0,nzcell-1);

        ijk = IJKCELL(i1,j1,k1);
        photon_cell[ijk]++;                 // record location of photon

        if(rand_1d(0.0,1.0)>albedo)break;   // if this if is true then the photon is absorbed
        rand_sphere_dir(xyzdir);
      }
      CheckMemory;
    }
    binmax = 0.0;
    binsum=0.0;
    for(i=0;i<nxcell*nycell*nzcell;i++){
      if(photon_cell[i]>binmax)binmax=photon_cell[i];
      binsum+=photon_cell[i];
    }
    printf("binmax=%f binsum=%f\n",binmax,binsum);
    if(binmax<1.0)binmax=1.0;
    CheckMemory;
    for(i=0;i<nx*ny*nz;i++){
      float val;
      int ival;
    
      val = average(i, photon_cell, smoke_mesh);
      ival = 254.0*val/(float)binmax;
      lightingbuffer[i]=(unsigned char)ival;
    }
    CheckMemory;
}

/* ------------------ in_mesh ------------------------ */

int in_mesh(mesh *smoke_mesh,float xyzpos[3]){
  if(xyzpos[0]<smoke_mesh->xbar0||xyzpos[0]>smoke_mesh->xbar)return 0;
  if(xyzpos[1]<smoke_mesh->ybar0||xyzpos[1]>smoke_mesh->ybar)return 0;
  if(xyzpos[2]<smoke_mesh->zbar0||xyzpos[2]>smoke_mesh->zbar)return 0;
  return 1;

}

/* ------------------ get_photon_step ------------------------ */

float get_photon_step(smoke3d *smoke3di, float xyzpos[3], float xyzdir[3]){
#define NPHOTON_PDF 50
  float photon_pdf[NPHOTON_PDF], photon_cdf[NPHOTON_PDF];
  int i;
  float xyz[3];
  float logalpha;
  mesh *smoke_mesh;
  float var, dlength;
  int ibreak;
  float step;

  smoke_mesh = smoke3di->smoke_mesh;


  dlength = smoke_mesh->dxyzmax/NPHOTON_PDF;
  ibreak = 0;

  photon_pdf[0]=0.0;
  photon_cdf[0]=0.0;
  var = rand_1d(0.0,1.0);
  for(i=1;i<NPHOTON_PDF;i++){
    if(ibreak==0){
      xyz[0] = xyzpos[0] + i*dlength*xyzdir[0];
      xyz[1] = xyzpos[1] + i*dlength*xyzdir[1];
      xyz[2] = xyzpos[2] + i*dlength*xyzdir[2];
      if(in_mesh(smoke_mesh,xyz)==0){
        ibreak=1;
      }
      logalpha = getlog_1_m_alpha(smoke_mesh,xyz);
      photon_pdf[i]=photon_pdf[i-1]+dlength*logalpha;
    }
    else{
      photon_pdf[i]=photon_pdf[i-1]-10.0;
    }
    if(photon_pdf[i]<-50.0){
      photon_cdf[i]=1.0;
    }
    else{
      photon_cdf[i]=1.0-exp(photon_pdf[i]);
    }
    if(photon_cdf[i-1]<var&&var<photon_cdf[i]){
      step = (i-1)*dlength;
      return step;
    }
  }
  step = (NPHOTON_PDF-2)*dlength;
  return step;
}

/* ------------------ getalpha ------------------------ */

float getlog_1_m_alpha(mesh *smoke_mesh, float *xyz){
  float val;
  float logalpha111, logalpha112, logalpha121, logalpha122;
  float logalpha211, logalpha212, logalpha221, logalpha222;
  float x1, x2, y1, y2, z1, z2;
  float f1=0.0, f2=1.0, g1=0.0, g2=1.0, h1=0.0, h2=1.0;
  int i1, j1, k1;
  int i2, j2, k2;
  int i2mi1, j2mj1;
  int nx, ny, nxy;
  float dx, dy, dz;
  int ijk111, ijk112, ijk121, ijk122;
  int ijk211, ijk212, ijk221, ijk222;

  if(xyz[0]<smoke_mesh->xbar0||xyz[0]>smoke_mesh->xbar)return 0.0;
  if(xyz[1]<smoke_mesh->ybar0||xyz[1]>smoke_mesh->ybar)return 0.0;
  if(xyz[2]<smoke_mesh->zbar0||xyz[2]>smoke_mesh->zbar)return 0.0;

  nx = smoke_mesh->ibar+1;
  ny = smoke_mesh->jbar+1;
  nxy = nx*ny;

  i1 = GET_INTERVAL(xyz[0],smoke_mesh->xplt[0],smoke_mesh->dx);
  j1 = GET_INTERVAL(xyz[1],smoke_mesh->yplt[0],smoke_mesh->dy);
  k1 = GET_INTERVAL(xyz[2],smoke_mesh->zplt[0],smoke_mesh->dz);

  i2=i1+1;
  if(i2>smoke_mesh->ibar)i2=smoke_mesh->ibar;
  i2mi1=i2-i1;

  j2=j1+1;
  if(j2>smoke_mesh->jbar)j2=smoke_mesh->jbar;
  j2mj1 = j2 - j1;

  k2=k1+1;
  if(k2>smoke_mesh->kbar)k2=smoke_mesh->kbar;
  
//#define IJKNODE(i,j,k) ((i)+(j)*nx+(k)*nxy)

  ijk111 = IJKNODE(i1,j1,k1);
  ijk211 = ijk111 + i2mi1;

  ijk121 = ijk111 + j2mj1*nx;
  ijk221 = ijk121 + i2mi1;

  if(k1!=k2){
    ijk112 = ijk111 + nxy;
    ijk212 = ijk211 + nxy;
    ijk122 = ijk121 + nxy;
    ijk222 = ijk221 + nxy;
  }
  else{
    ijk112 = ijk111;
    ijk212 = ijk211;
    ijk122 = ijk121;
    ijk222 = ijk221;
  }

  logalpha111 = full_logalphabuffer[ijk111];
  logalpha211 = full_logalphabuffer[ijk211];

  logalpha121 = full_logalphabuffer[ijk121];
  logalpha221 = full_logalphabuffer[ijk221];

  logalpha112 = full_logalphabuffer[ijk112];
  logalpha212 = full_logalphabuffer[ijk212];

  logalpha122 = full_logalphabuffer[ijk122];
  logalpha222 = full_logalphabuffer[ijk222];

  x1 = smoke_mesh->xplt[i1];
  x2 = smoke_mesh->xplt[i2];
  y1 = smoke_mesh->yplt[j1];
  y2 = smoke_mesh->yplt[j2];
  z1 = smoke_mesh->zplt[k1];
  z2 = smoke_mesh->zplt[k2];

  dx = x2 - x1;
  dy = y2 - y1;
  dz = z2 - z1;

  if(dx!=0.0&&dy!=0.0&&dz!=0.0){
    f1 = xyz[0]-x1;
    f2 = x2 - xyz[0];

    g1 = xyz[1]-y1;
    g2 = y2 - xyz[1];

    h1 = xyz[2]-z1;
    h2 = z2 - xyz[2];

    val  = f2*(g2*(h2*logalpha111 + h1*logalpha112) + g1*(h2*logalpha121 + h1*logalpha122));
    val += f1*(g2*(h2*logalpha211 + h1*logalpha212) + g1*(h2*logalpha221 + h1*logalpha222));
    val /= (dx*dy*dz);
  }
  else{
    val = logalpha111;
  }
  return val;
}

/* ------------------ interp ------------------------ */

float interp(float xyz[3], mesh *smoke_mesh, float *full_logalphabuffer){
  float logalpha111, logalpha112, logalpha121, logalpha122;
  float logalpha211, logalpha212, logalpha221, logalpha222;
  float x1, x2, y1, y2, z1, z2;
  float f1=0.0, f2=1.0, g1=0.0, g2=1.0, h1=0.0, h2=1.0;
  float val;
  int i1, j1, k1;
  int i2, j2, k2;
  int i2mi1, j2mj1;
  int nx, ny, nxy;
  float dx, dy, dz;
  int ijk111, ijk112, ijk121, ijk122;
  int ijk211, ijk212, ijk221, ijk222;

  nx = smoke_mesh->ibar+1;
  ny = smoke_mesh->jbar+1;
  nxy = nx*ny;

  i1 = GET_INTERVAL(xyz[0],smoke_mesh->xplt[0],smoke_mesh->dx);
  i1= BOUND(i1,0,smoke_mesh->ibar);
  j1 = GET_INTERVAL(xyz[1],smoke_mesh->yplt[0],smoke_mesh->dy);
  j1= BOUND(j1,0,smoke_mesh->ibar);
  k1 = GET_INTERVAL(xyz[2],smoke_mesh->zplt[0],smoke_mesh->dz);
  k1= BOUND(k1,0,smoke_mesh->ibar);

  i2=i1+1;
  if(i2>smoke_mesh->ibar)i2=smoke_mesh->ibar;
  i2mi1=i2-i1;

  j2=j1+1;
  if(j2>smoke_mesh->jbar)j2=smoke_mesh->jbar;
  j2mj1 = j2 - j1;

  k2=k1+1;
  if(k2>smoke_mesh->kbar)k2=smoke_mesh->kbar;
  
//#define IJKNODE(i,j,k) ((i)+(j)*nx+(k)*nxy)

  ijk111 = IJKNODE(i1,j1,k1);
  ijk211 = ijk111 + i2mi1;

  ijk121 = ijk111 + j2mj1*nx;
  ijk221 = ijk121 + i2mi1;

  if(k1!=k2){
    ijk112 = ijk111 + nxy;
    ijk212 = ijk211 + nxy;
    ijk122 = ijk121 + nxy;
    ijk222 = ijk221 + nxy;
  }
  else{
    ijk112 = ijk111;
    ijk212 = ijk211;
    ijk122 = ijk121;
    ijk222 = ijk221;
  }

  logalpha111 = full_logalphabuffer[ijk111];
  logalpha211 = full_logalphabuffer[ijk211];

  logalpha121 = full_logalphabuffer[ijk121];
  logalpha221 = full_logalphabuffer[ijk221];

  logalpha112 = full_logalphabuffer[ijk112];
  logalpha212 = full_logalphabuffer[ijk212];

  logalpha122 = full_logalphabuffer[ijk122];
  logalpha222 = full_logalphabuffer[ijk222];

  x1 = smoke_mesh->xplt[i1];
  x2 = smoke_mesh->xplt[i2];
  y1 = smoke_mesh->yplt[j1];
  y2 = smoke_mesh->yplt[j2];
  z1 = smoke_mesh->zplt[k1];
  z2 = smoke_mesh->zplt[k2];

  dx = x2 - x1;
  dy = y2 - y1;
  dz = z2 - z1;

  if(dx!=0.0&&dy!=0.0&&dz!=0.0){
    f1 = xyz[0]-x1;
    f2 = x2 - xyz[0];

    g1 = xyz[1]-y1;
    g2 = y2 - xyz[1];

    h1 = xyz[2]-z1;
    h2 = z2 - xyz[2];

    val  = f2*(g2*(h2*logalpha111 + h1*logalpha112) + g1*(h2*logalpha121 + h1*logalpha122));
    val += f1*(g2*(h2*logalpha211 + h1*logalpha212) + g1*(h2*logalpha221 + h1*logalpha222));
    val /= (dx*dy*dz);
  }
  else{
    val = logalpha111;
  }
  return val;
}

/* ------------------ average ------------------------ */

float average(int ijknode, float *celldata, mesh *smoke_mesh){
  int nx, ny, nxy;
  int nxcell, nycell, nxycell, nzcell;
  int i1, j1, k1;
  int i2, j2, k2;
  int ijk;
  int i111, i112, i121, i122;
  int i211, i212, i221, i222;
  int v111, v112, v121, v122;
  int v211, v212, v221, v222;
  float val;

  nxcell = smoke_mesh->ibar;
  nycell = smoke_mesh->jbar;
  nzcell = smoke_mesh->kbar;
  nxycell = nxcell*nycell;
  nx = nxcell + 1;
  ny = nycell + 1;
  nxy = nx*ny;

//#define IJKNODE(i,j,k) ((i)+(j)*nx+(k)*nxy)
  ijk = ijknode;
  k2 = ijk/nxy;
  ijk-=k2*nxy;
  j2 = ijk/nx;
  i2 = ijk-j2*nx;

  if(i2>nxcell-1)i2=nxcell-1;
  if(j2>nycell-1)j2=nycell-1;
  if(k2>nzcell-1)k2=nzcell-1;

  k1=k2-1;
  if(k1<0)k1=0;

  j1=j2-1;
  if(j1<0)j1=0;

  i1=i2-1;
  if(i1<0)i1=0;
  
  i111 = IJKCELL(i1,j1,k1);
  i211 = i111 + (i2-i1);
  i121 = i111 + (j2-j1)*nxcell;
  i221 = i121 + (i2-i1);

  if(k1!=k2){
    i112 = i111 + nxycell;
    i212 = i211 + nxycell;
    i122 = i121 + nxycell;
    i222 = i221 + nxycell;
  }
  else{
    i112 = i111;
    i212 = i211;
    i122 = i121;
    i222 = i221;
  }

  v111 = celldata[i111];
  v112 = celldata[i112];
  v121 = celldata[i121];
  v122 = celldata[i122];
  v211 = celldata[i211];
  v212 = celldata[i212];
  v221 = celldata[i221];
  v222 = celldata[i222];


  val  = (v111 + v112 + v121 + v122 + v211 + v212 + v221 + v222)/8.0;
  return val;
}

/* ------------------ get_interval ------------------------ */

int get_interval(float val, float *array, int n){
  int low, mid, high;

  if(val<array[0])return 0;
  if(val>array[n-1])return n-1;

  low=0;
  high=n-1;
  while(high-low>1){
    mid=(low+high)/2;
    if(val>array[mid]){
      low=mid;
    }
    else{
      high=mid;
    }
  }
  ASSERT(low<n)
  return low;
}


#endif


