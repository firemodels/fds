#define IN_UPDATE
#include "options.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include GLUT_H

#include "update.h"
#include "smokeviewvars.h"
#include "compress.h"

void ParticleStreakShowMenu(int var);

/* ------------------ CompareFloat ------------------------ */

int CompareFloat( const void *arg1, const void *arg2 ){
  float x, y;

  x=*(float *)arg1;
  y=*(float *)arg2;
  if( x< y)return -1;
  if( x> y)return 1;
  return 0;
}

/* ------------------ UpdateHrrinfo ------------------------ */

void UpdateHrrinfo(int vis) {
  if(hrrinfo != NULL) {
    hrrinfo->display = vis;
	UpdateTimes();
  }
}

/* ------------------ UpdateFrameNumber ------------------------ */

void UpdateFrameNumber(int changetime){
  if(force_redisplay==1||(itimeold!=itimes&&changetime==1)){
    int i;

    force_redisplay=0;
    itimeold=itimes;
    if(showsmoke==1||showevac==1){
      for(i=0;i<npartinfo;i++){
        partdata *parti;

        parti = partinfo+i;
        if(parti->loaded==0||parti->timeslist==NULL)continue;
        parti->itime=parti->timeslist[itimes];
      }
    }
    if(hrrinfo!=NULL&&hrrinfo->loaded==1&&hrrinfo->display==1&&hrrinfo->timeslist!=NULL){
      hrrinfo->itime=hrrinfo->timeslist[itimes];
    }
    if(showvolrender==1){
      int imesh;

      for(imesh=0;imesh<nmeshes;imesh++){
        meshdata *meshi;
        volrenderdata *vr;
        slicedata *fireslice, *smokeslice;
        int j;

        meshi = meshinfo + imesh;
        vr = &(meshi->volrenderinfo);
        fireslice=vr->fireslice;
        smokeslice=vr->smokeslice;
        if(fireslice==NULL||smokeslice==NULL)continue;
        if(vr->loaded==0||vr->display==0)continue;
        vr->itime = vr->timeslist[itimes];
        for(j=vr->itime;j>=0;j--){
          if(vr->dataready[j]==1)break;
        }
        vr->itime=j;
        if(smokeslice!=NULL&&vr->itime>=0){
          if(vr->is_compressed==1||load_volcompressed==1){
            unsigned char *c_smokedata_compressed;
            uLongf framesize;
            float timeval;

            c_smokedata_compressed = vr->smokedataptrs[vr->itime];
            framesize = smokeslice->nslicei*smokeslice->nslicej*smokeslice->nslicek;
            uncompress_volsliceframe(c_smokedata_compressed,
                           vr->smokedata_view, framesize, &timeval,
                           vr->c_smokedata_view);

            vr->smokedataptr = vr->smokedata_view;
          }
          else{
            if(runscript==0)vr->smokedataptr = vr->smokedataptrs[vr->itime];
          }
          CheckMemory;
        }
        if(fireslice!=NULL&&vr->itime>=0){
          if(vr->is_compressed==1||load_volcompressed==1){
            unsigned char *c_firedata_compressed;
            uLongf framesize;
            float timeval;

            c_firedata_compressed = vr->firedataptrs[vr->itime];
            framesize = fireslice->nslicei*fireslice->nslicej*fireslice->nslicek;
            uncompress_volsliceframe(c_firedata_compressed,
                           vr->firedata_view, framesize, &timeval,
                           vr->c_firedata_view);

            vr->firedataptr = vr->firedata_view;
            CheckMemory;
          }
          else{
            if(runscript==0)vr->firedataptr = vr->firedataptrs[vr->itime];
          }
          CheckMemory;
        }
      }
    }
    for(i=0;i<ngeominfoptrs;i++){
      geomdata *geomi;

      geomi = geominfoptrs[i];
      if(geomi->loaded==0||geomi->timeslist==NULL)continue;
      geomi->itime=geomi->timeslist[itimes];
    }
    if(showslice==1||showvslice==1){
      int showfed=0;
      int ii;

      for(ii=0;ii<nslice_loaded;ii++){
        slicedata *sd;

        i = slice_loaded_list[ii];
        sd = sliceinfo+i;
        if(sd->timeslist==NULL)continue;
        sd->itime=sd->timeslist[itimes];
        if(sd->is_fed==1)showfed=1;
      }
      if(showfed==1){
        fed_areas=NULL;
        for(i=0;i<nmultisliceinfo;i++){
          multislicedata *mslicei;
          int itime;
          int j;

          mslicei = multisliceinfo + i;
          if(mslicei->contour_areas_percen==NULL)continue;
          if(fed_areas!=NULL){
            fed_areas=NULL;
            break;
          }
          for(j=0;j<mslicei->nslices;j++){
            slicedata *slicej;

            slicej = sliceinfo + mslicei->islices[j];
            if(slicej->loaded==0||slicej->display==0)continue;
            itime = slicej->itime;
            if(mslicei->nslices==1){
              fed_areas=slicej->contour_areas_percen+4*itime;
            }
            else{
              fed_areas=mslicei->contour_areas_percen+4*itime;
            }
            break;
          }
        }
      }
    }
    if(show3dsmoke==1){
      for(i=0;i<nsmoke3dinfo;i++){
        smoke3ddata *smoke3di;

        smoke3di = smoke3dinfo + i;
        if(smoke3di->loaded==0||smoke3di->display==0)continue;
        smoke3di->ismoke3d_time=smoke3di->timeslist[itimes];
        if(smoke3di->ismoke3d_time!=smoke3di->lastiframe){
          smoke3di->lastiframe=smoke3di->ismoke3d_time;
          updatesmoke3d(smoke3di);
        }
      }
      if(nsmoke3dinfo>0)mergesmoke3dcolors(NULL);
    }
    if(showpatch==1){
      for(i=0;i<npatchinfo;i++){
        patchdata *patchi;

        patchi = patchinfo + i;
        if(patchi->filetype!=PATCH_GEOMETRY||patchi->geom_times==NULL||patchi->geom_timeslist==NULL)continue;
        patchi->geom_itime=patchi->geom_timeslist[itimes];
        patchi->geom_ival_static = patchi->geom_ivals_static[patchi->geom_itime];
        patchi->geom_ival_dynamic = patchi->geom_ivals_dynamic[patchi->geom_itime];
        patchi->geom_nval_static = patchi->geom_nstatics[patchi->geom_itime];
        patchi->geom_nval_dynamic = patchi->geom_ndynamics[patchi->geom_itime];
      }
      for(i=0;i<nmeshes;i++){
        patchdata *patchi;
        meshdata *meshi;

        meshi = meshinfo+i;
        patchi=patchinfo + meshi->patchfilenum;
        if(patchi->filetype==PATCH_GEOMETRY||meshi->patch_times==NULL||meshi->patch_timeslist==NULL)continue;
        meshi->patch_itime=meshi->patch_timeslist[itimes];
        if(patchi->compression_type==UNCOMPRESSED){
          meshi->cpatchval_iframe = meshi->cpatchval + meshi->patch_itime*meshi->npatchsize;
        }
        else{
          uncompress_patchdataframe(meshi,meshi->patch_itime);
        }
      }
    }
    if(showiso==1){
      isodata *isoi;
      meshdata *meshi;

      CheckMemory;
      for(i=0;i<nisoinfo;i++){
        isoi = isoinfo + i;
        meshi = meshinfo + isoi->blocknumber;
        if(isoi->loaded==0||meshi->iso_times==NULL||meshi->iso_timeslist==NULL)continue;
        meshi->iso_itime=meshi->iso_timeslist[itimes];
      }
    }
    if(showzone==1){
      izone=zone_timeslist[itimes];
    }
  }
}

/* ------------------ UpdateShow ------------------------ */

void UpdateShow(void){
  int i,evacflag,sliceflag,vsliceflag,partflag,patchflag,isoflag,smoke3dflag,tisoflag;
  int slicecolorbarflag;
  int shooter_flag;

  showtime=0;
  showtime2=0;
  showplot3d=0;
  showpatch=0;
  showslice=0;
  showvslice=0;
  showsmoke=0;
  showzone=0;
  showiso=0;
  showvolrender=0;
  have_extreme_mindata=0;
  have_extreme_maxdata=0;
  showshooter=0;
  showevac=0;
  showevac_colorbar=0;
  show3dsmoke=0;
  smoke3dflag=0;
  showtours=0;
  showterrain=0;
  visTimeParticles=1; visTimeSlice=1; visTimePatch=1; visTimeZone=1; visTimeIso=1;

  RenderTime=0;
  if(global_times!=NULL){
    if(settmin_p==1&&global_times[itimes]<tmin_p)visTimeParticles=0;
    if(settmax_p==1&&global_times[itimes]>tmax_p)visTimeParticles=0;

    if(settmin_s==1&&global_times[itimes]<tmin_s)visTimeSlice=0;
    if(settmax_s==1&&global_times[itimes]>tmax_s)visTimeSlice=0;

    if(settmin_i==1&&global_times[itimes]<tmin_i)visTimeIso=0;
    if(settmax_i==1&&global_times[itimes]>tmax_i)visTimeIso=0;

    if(settmin_b==1&&global_times[itimes]<tmin_b)visTimePatch=0;
    if(settmax_b==1&&global_times[itimes]>tmax_b)visTimePatch=0;

    if(settmin_z==1&&global_times[itimes]<tmin_z)visTimeZone=0;
    if(settmax_z==1&&global_times[itimes]>tmax_z)visTimeZone=0;

  }

  {
    tourdata *touri;

    if(ntours>0){
      for(i=0;i<ntours;i++){
        touri = tourinfo + i;
        if(touri->display==1){
          showtours=1;
          break;
        }
      }
    }
  }
  if(visTerrainType!=TERRAIN_HIDDEN){
    for(i=0;i<nterraininfo;i++){
      terraindata *terri;

      terri = terraininfo + i;
      if(terri->loaded==1){
        showterrain=1;
        break;
      }
    }
  }
  {
    smoke3ddata *smoke3di;
    int ii;

    for(ii=0;ii<nsmoke3dinfo;ii++){
      smoke3di = smoke3dinfo + ii;
      if(smoke3di->loaded==1&&smoke3di->display==1){
        smoke3dflag=1;
        break;
      }
    }
  }
  if(nvolrenderinfo>0&&usevolrender==1){
    for(i=0;i<nmeshes;i++){
      meshdata *meshi;
      volrenderdata *vr;

      meshi = meshinfo + i;
      vr = &(meshi->volrenderinfo);
      if(vr->fireslice==NULL||vr->smokeslice==NULL)continue;
      if(vr->loaded==0||vr->display==0)continue;
      showvolrender=1;
      break;
    }
  }
  sliceflag=0;
  slicecolorbarflag=0;
  SHOW_gslice_data=0;
  if(visTimeSlice==1){
    int ii;

    for(ii=0;ii<nslice_loaded;ii++){
      slicedata *sd;

      i=slice_loaded_list[ii];
      sd = sliceinfo+i;
      if(sd->display==0||sd->type!=islicetype)continue;
      if(sd->volslice==1&&sd->slicetype==SLICE_NODE_CENTER&&vis_gslice_data==1)SHOW_gslice_data=1;
      if(sd->ntimes>0){
        sliceflag=1;
        break;
      }
    }
    if(SHOW_gslice_data!=SHOW_gslice_data_old){
      SHOW_gslice_data_old=SHOW_gslice_data;
      updatemenu=1;
    }
    for(ii=0;ii<nslice_loaded;ii++){
      meshdata *slicemesh;
      slicedata *sd;

      i=slice_loaded_list[ii];
      sd = sliceinfo+i;
      slicemesh = meshinfo + sd->blocknumber;
      if(sd->display==0||sd->type!=islicetype)continue;
      if(sd->constant_color==NULL&&show_evac_colorbar==0&&slicemesh->mesh_type!=0)continue;
      if(sd->constant_color!=NULL)continue;
      if(sd->ntimes>0){
        slicecolorbarflag=1;
        break;
      }
    }
    if(have_extreme_maxdata==0){
      for(ii=0;ii<nslice_loaded;ii++){
        slicedata *sd;

        i=slice_loaded_list[ii];
        sd = sliceinfo+i;
        if(sd->display==0||sd->type!=islicetype)continue;
        if(sd->extreme_max==1){
          have_extreme_maxdata=1;
          break;
        }
      }
    }
    if(have_extreme_mindata==0){
      for(ii=0;ii<nslice_loaded;ii++){
        slicedata *sd;

        i=slice_loaded_list[ii];
        sd = sliceinfo+i;
        if(sd->display==0||sd->type!=islicetype)continue;
        if(sd->extreme_min==1){
          have_extreme_mindata=1;
          break;
        }
      }
    }
  }
  isoflag=0;
  tisoflag=0;
  if(visTimeIso==1){
    for(i=0;i<nisoinfo;i++){
      isodata *isoi;

      isoi = isoinfo+i;
      if(isoi->loaded==0)continue;
      if(isoi->display==1&&isoi->type==iisotype){
        isoflag=1;
        if(isoi->dataflag==1){
          tisoflag=1;
          break;
        }
      }
    }
  }
  vsliceflag=0;
  vslicecolorbarflag=0;
  if(visTimeSlice==1){
    for(i=0;i<nvsliceinfo;i++){
      vslicedata *vd;
      slicedata *sd;

      vd = vsliceinfo+i;
      if(vd->loaded==0||vd->display==0)continue;
      sd = sliceinfo + vd->ival;

      if(sd->type!=islicetype)continue;
      if(sd->volslice==1&&sd->slicetype==SLICE_NODE_CENTER&&vis_gslice_data==1)SHOW_gslice_data=1;
      vsliceflag=1;
      break;
    }
    for(i=0;i<nvsliceinfo;i++){
      meshdata *slicemesh;
      slicedata *sd;
      vslicedata *vd;

      vd = vsliceinfo+i;
      sd = sliceinfo + vd->ival;
      slicemesh = meshinfo + sd->blocknumber;
      if(vd->loaded==0||vd->display==0)continue;
      if(sliceinfo[vd->ival].type!=islicetype)continue;
      if(sd->constant_color==NULL&&show_evac_colorbar==0&&slicemesh->mesh_type!=0)continue;
      if(sd->constant_color!=NULL)continue;
      vslicecolorbarflag=1;
      break;
    }
  }
  patchflag=0;
  if(visTimePatch==1){
    int ii;

    wc_flag=0;
    for(ii=0;ii<npatch_loaded;ii++){
      patchdata *patchi;

      i = patch_loaded_list[ii];
      patchi=patchinfo+i;
      if(patchi->display==0||patchi->type!=ipatchtype)continue;
      if(strcmp(patchi->label.shortlabel,"wc")==0)wc_flag=1;
      patchflag=1;
      break;
    }
    for(ii=0;ii<npatch_loaded;ii++){
      patchdata *patchi;

      i = patch_loaded_list[ii];
      patchi=patchinfo+i;
      if(patchi->display==0||patchi->type!=ipatchtype)continue;
      if(patchi->extreme_max==1){
        have_extreme_maxdata=1;
        break;
      }
    }
    for(ii=0;ii<npatch_loaded;ii++){
      patchdata *patchi;

      i = patch_loaded_list[ii];
      patchi=patchinfo+i;
      if(patchi->display==0||patchi->type!=ipatchtype)continue;
      if(patchi->extreme_min==1){
        have_extreme_mindata=1;
        break;
      }
    }
    for(i = 0; i < ngeominfo; i++){
      geomdata *geomi;

      geomi = geominfo + i;
      geomi->patchactive = 0;
    }
    for(ii=0;ii<npatch_loaded;ii++){
      patchdata *patchi;

      i = patch_loaded_list[ii];
      patchi=patchinfo+i;
      if(patchi->geominfo!=NULL&&patchi->display == 1 && patchi->type == ipatchtype){
        patchi->geominfo->patchactive = 1;
      }
    }
  }
  partflag=0;
  if(visParticles==1&&visTimeParticles==1){
    for(i=0;i<npartinfo;i++){
      partdata *parti;

      parti = partinfo + i;
      if(parti->evac==0&&parti->loaded==1&&parti->display==1){
        partflag=1;
        break;
      }
    }
    if(current_property!=NULL){
      if(current_property->extreme_max==1)have_extreme_maxdata=1;
      if(current_property->extreme_min==1)have_extreme_mindata=1;
    }
  }
  evacflag=0;
  if(visEvac==1&&visTimeEvac==1){
    for(i=0;i<npartinfo;i++){
      partdata *parti;

      parti = partinfo + i;
      if(parti->evac==1&&parti->loaded==1&&parti->display==1){
        evacflag=1;
        break;
      }
    }
  }
  shooter_flag=0;
  if(visShooter!=0&&shooter_active==1){
    shooter_flag=1;
  }

  if( plotstate==DYNAMIC_PLOTS &&
    ( sliceflag==1 || vsliceflag==1 || partflag==1 || patchflag==1 ||
    shooter_flag==1||
    smoke3dflag==1|| showtours==1 || evacflag==1||
    (ReadZoneFile==1&&visZone==1&&visTimeZone==1)||
    showterrain==1||showvolrender==1
    )
    )showtime=1;
  if(plotstate==DYNAMIC_PLOTS&&ReadIsoFile==1&&visAIso!=0&&isoflag==1)showtime2=1;
  if(plotstate==DYNAMIC_PLOTS){
    if(smoke3dflag==1)show3dsmoke=1;
    if(partflag==1)showsmoke=1;
    if(evacflag==1)showevac=1;
    if(showevac==1&&parttype>0){
      showevac_colorbar=1;
      if(current_property!=NULL&&strcmp(current_property->label->longlabel,"HUMAN_COLOR")==0){
        showevac_colorbar=0;
      }
    }
    if(patchflag==1)showpatch=1;
    for(i=0;i<nmeshes;i++){
      meshdata *meshi;

      meshi=meshinfo+i;
      meshi->visInteriorPatches=0;
    }
    if(showpatch==1&&visPatchType[0]==1){
      for(i=0;i<nmeshes;i++){
        patchdata *patchi;
        meshdata *meshi;

        meshi=meshinfo+i;
        if(meshi->patch_times==NULL)continue;
        patchi = patchinfo+meshi->patchfilenum;
        if(patchi->loaded==1&&patchi->display==1&&patchi->type==ipatchtype){
          meshi->visInteriorPatches=1;
        }
      }
    }
    if(sliceflag==1)showslice=1;
    if(vsliceflag==1)showvslice=1;
    if(ReadZoneFile==1&&visZone==1&&visTimeZone==1)showzone=1;
    if(ReadIsoFile==1&&visAIso!=0){
      showiso=1;
    }
    if(shooter_flag==1)showshooter=1;
  }
  if(showsmoke==1||showevac==1||showpatch==1||showslice==1||showvslice==1||showzone==1||showiso==1||showevac==1)RenderTime=1;
  if(showtours==1||show3dsmoke==1||touring==1||showvolrender==1)RenderTime=1;
  if(showshooter==1)RenderTime=1;
  if(plotstate==STATIC_PLOTS&&ReadPlot3dFile==1&&plotn>0&&plotn<=numplot3dvars)showplot3d=1;
  if(showplot3d==1){
    for(i=0;i<nmeshes;i++){
      meshdata *meshi;
      int ii;

      meshi=meshinfo+i;
      ii=meshi->plot3dfilenum;
      if(ii==-1)continue;
      if(plot3dinfo[ii].loaded==0)continue;
      if(plot3dinfo[ii].display==0)continue;
      if(plot3dinfo[ii].extreme_min[plotn-1]==1)have_extreme_mindata=1;
    }
    for(i=0;i<nmeshes;i++){
      meshdata *meshi;
      int ii;

      meshi=meshinfo+i;
      ii=meshi->plot3dfilenum;
      if(ii==-1)continue;
      if(plot3dinfo[ii].loaded==0)continue;
      if(plot3dinfo[ii].display==0)continue;
      if(plot3dinfo[ii].extreme_max[plotn-1]==1)have_extreme_maxdata=1;
    }
  }

  numColorbars=0;
  if(ReadEvacFile==1)numColorbars++;
  if(ReadPartFile==1)numColorbars++;
  if(plotstate==DYNAMIC_PLOTS&&(slicecolorbarflag==1||vslicecolorbarflag==1))numColorbars++;
  if(plotstate==DYNAMIC_PLOTS&&patchflag==1&&wc_flag==0)numColorbars++;
  if(plotstate==DYNAMIC_PLOTS&&ReadZoneFile==1)numColorbars++;
  if(plotstate==DYNAMIC_PLOTS&&tisoflag==1){
    showiso_colorbar=1;
    numColorbars++;
  }
  if(ReadPlot3dFile==1&&numColorbars==0)numColorbars=1;
  /* note: animated iso-contours do not need a color bar,
           so we don't test for isosurface files */
  drawColorLabel=0;
  if((showtime==1||showplot3d==1)&&visColorbar==1)drawColorLabel=1;
  if(drawColorLabel==1&&olddrawColorLabel==0)updatemenu=1;
  if(drawColorLabel==0&&olddrawColorLabel==1)updatemenu=1;
  olddrawColorLabel=drawColorLabel;
  if(showtime2==1)showtime=1;
  if(plotstate==DYNAMIC_PLOTS&&stept==1){
    glutIdleFunc(Idle_CB);
  }
  else{
    glutIdleFunc(NULL);
  }
}

/* ------------------ GetItime ------------------------ */

int GetItime(int n, int *timeslist, float *times, int ntimes){
  int istart=0;

  if(n>0)istart=timeslist[n-1];
  while(istart<ntimes-1&&times[istart+1]<=global_times[n]){
    istart++;
  }
  istart=CLAMP(istart,0,ntimes-1);
  return istart;
}

/* ------------------ SynchTimes ------------------------ */

void SynchTimes(void){
  int n,i,istart,igrid;

  /* synchronize smooth blockage times */

  for(n=0;n<nglobal_times;n++){
    int j,jj;

  /* synchronize tour times */

    for(j=0;j<ntours;j++){
      tourdata *tourj;

      tourj = tourinfo + j;
      if(tourj->display==0)continue;
      tourj->timeslist[n]=GetItime(n,tourj->timeslist,tourj->path_times,tourj->ntimes);
    }

    /* synchronize terrain times */

    for(j=0;j<nterraininfo;j++){
      terraindata *terri;

      terri = terraininfo + j;
      if(terri->loaded==0)continue;
      terri->timeslist[n]=GetItime(n,terri->timeslist,terri->times,terri->ntimes);
    }
    if(hrrinfo!=NULL&&hrrinfo->loaded==1&&hrrinfo->display==1){
      hrrinfo->timeslist[n]=GetItime(n,hrrinfo->timeslist,hrrinfo->times,hrrinfo->ntimes);
    }

  /* synchronize geometry times */

    for(j=0;j<ngeominfoptrs;j++){
      geomdata *geomi;

      geomi = geominfoptrs[j];
      if(geomi->loaded==0||geomi->display==0)continue;
      geomi->timeslist[n]=GetItime(n,geomi->timeslist,geomi->times,geomi->ntimes);
    }

  /* synchronize particle times */

    for(j=0;j<npartinfo;j++){
      partdata *parti;

      parti=partinfo+j;
      if(parti->loaded==0)continue;
      parti->timeslist[n]=GetItime(n,parti->timeslist,parti->times,parti->ntimes);
    }

  /* synchronize shooter times */
    if(visShooter!=0&&shooter_active==1){
      if(n==0){
        istart=0;
      }
      else{
        istart=shooter_timeslist[n-1];
      }
      i=istart;
      while(shoottimeinfo[i].time<global_times[n]&&i<nshooter_frames){
        i++;
      }
      if(i>=nshooter_frames){
        i=nshooter_frames-1;
      }
      shooter_timeslist[n]=i;
    }

  /* synchronize slice times */

    for(jj=0;jj<nslice_loaded;jj++){
      slicedata *sd;

      j = slice_loaded_list[jj];
      sd = sliceinfo + j;
      sd->timeslist[n]=GetItime(n,sd->timeslist,sd->times,sd->ntimes);
    }

  /* synchronize smoke times */
    {
      smoke3ddata *smoke3di;

      for(jj=0;jj<nsmoke3dinfo;jj++){
        smoke3di = smoke3dinfo + jj;
        if(smoke3di->loaded==0)continue;
        smoke3di->timeslist[n]=GetItime(n,smoke3di->timeslist,smoke3di->times,smoke3di->ntimes);
      }
    }

  /* synchronize patch times */

    for(j=0;j<npatchinfo;j++){
      patchdata *patchi;

      patchi = patchinfo + j;
      if(patchi->loaded==0)continue;
      if(patchi->filetype != PATCH_GEOMETRY)continue;
      patchi->geom_timeslist[n]=GetItime(n,patchi->geom_timeslist,patchi->geom_times,patchi->ngeom_times);
    }
    for(j=0;j<nmeshes;j++){
      patchdata *patchi;
      meshdata *meshi;

      meshi=meshinfo+j;
      if(meshi->patchfilenum<0||meshi->patch_times==NULL)continue;
      patchi=patchinfo+meshi->patchfilenum;
      if(patchi->filetype==PATCH_GEOMETRY)continue;
      meshi->patch_timeslist[n]=GetItime(n,meshi->patch_timeslist,meshi->patch_times,meshi->npatch_times);
    }

  /* synchronize isosurface times */

    for(igrid=0;igrid<nmeshes;igrid++){
      meshdata *meshi;

      meshi=meshinfo+igrid;
      if(meshi->iso_times==NULL)continue;
      meshi->iso_timeslist[n]=GetItime(n,meshi->iso_timeslist,meshi->iso_times,meshi->niso_times);
    }

  /* synchronize volume render times */

    if(nvolrenderinfo>0){
      for(igrid=0;igrid<nmeshes;igrid++){
        volrenderdata *vr;
        meshdata *meshi;

        meshi=meshinfo+igrid;
        vr = &meshi->volrenderinfo;
        if(vr->smokeslice==NULL)continue;
        if(vr->loaded==0||vr->display==0)continue;
        if(vr->times==NULL)continue;
        vr->timeslist[n]=GetItime(n,vr->timeslist,vr->times,vr->ntimes);
      }
    }
    /* synchronize zone times */

    if(showzone==1){
      zone_timeslist[n]=GetItime(n,zone_timeslist,zone_times,nzone_times);
    }

  }
  reset_gltime();
}

/* ------------------ GetLoadvfileinfo ------------------------ */

int GetLoadvfileinfo(FILE *stream, char *filename){
  int i;
  char *fileptr;

  trim_back(filename);
  fileptr = trim_front(filename);
  for(i = 0; i<nsliceinfo; i++){
    slicedata *slicei;

    slicei = sliceinfo+i;
    if(strcmp(fileptr, slicei->file)==0){
      fprintf(stream, "// LOADVFILE\n");
      fprintf(stream, "//  %s\n", slicei->file);
      fprintf(stream, "LOADVSLICEM\n");
      fprintf(stream, " %s\n", slicei->label.longlabel);
      if(slicei->volslice==1){
        fprintf(stream, " %i %f\n", 0, slicei->position_orig);
      }
      else{
        fprintf(stream, " %i %f\n", slicei->idir, slicei->position_orig);
      }
      fprintf(stream, " %i\n", slicei->blocknumber+1);
      return 1;
    }
  }
  return 0;
}

/* ------------------ GetLoadfileinfo ------------------------ */

int GetLoadfileinfo(FILE *stream, char *filename){
  int i;
  char *fileptr;

  trim_back(filename);
  fileptr = trim_front(filename);
  for(i = 0; i<nsliceinfo; i++){
    slicedata *slicei;

    slicei = sliceinfo+i;
    if(strcmp(fileptr, slicei->file)==0){
      fprintf(stream, "// LOADFILE\n");
      fprintf(stream, "//  %s\n", slicei->file);
      fprintf(stream, "LOADSLICEM\n");
      fprintf(stream, " %s\n", slicei->label.longlabel);
      if(slicei->volslice==1){
        fprintf(stream, " %i %f\n", 0, slicei->position_orig);
      }
      else{
        fprintf(stream, " %i %f\n", slicei->idir, slicei->position_orig);
      }
      fprintf(stream, " %i\n", slicei->blocknumber+1);
      return 1;
    }
  }
  for(i = 0; i < nisoinfo; i++){
    isodata *isoi;

    isoi = isoinfo + i;
    if(strcmp(fileptr, isoi->file) == 0){
      fprintf(stream, "// LOADFILE\n");
      fprintf(stream, "//  %s\n", isoi->file);
      fprintf(stream, "LOADISOM\n");
      fprintf(stream, " %s\n", isoi->surface_label.longlabel);
      fprintf(stream, " %i\n", isoi->blocknumber+1);
      return 1;
    }

  }
  for(i = 0; i < npatchinfo; i++){
    patchdata *patchi;

    patchi = patchinfo + i;
    if(strcmp(fileptr, patchi->file) == 0){
      fprintf(stream, "// LOADFILE\n");
      fprintf(stream, "//  %s\n", patchi->file);
      fprintf(stream, "LOADBOUNDARYM\n");
      fprintf(stream, " %s\n", patchi->label.longlabel);
      fprintf(stream, " %i\n", patchi->blocknumber+1);
      return 1;
    }

  }
  return 0;
}

  /* ------------------ ConvertSsf ------------------------ */

void ConvertSsf(void){
  FILE *stream_from, *stream_to;
  int outeqin = 0;
#define LENTEMP 20
  char tempfile[LENTEMP+1];
  char *template = "tempssf";

  if(ssf_from==NULL||ssf_to==NULL)return;
  stream_from = fopen(ssf_from, "r");
  if(stream_from==NULL)return;

  if(strcmp(ssf_from, ssf_to)==0){
    strcpy(tempfile, template);
    if(randstr(tempfile+strlen(template), LENTEMP-strlen(template))==NULL||strlen(tempfile)==0){
      fclose(stream_from);
      return;
    }
    stream_to = fopen(tempfile, "w");
    outeqin = 1;
  }
  else{
    stream_to = fopen(ssf_to, "w");
  }
  if(stream_to==NULL){
    fclose(stream_from);
    return;
  }

  while(!feof(stream_from)){
    char buffer[255], filename[255];

    CheckMemory;
    if(fgets(buffer, 255, stream_from)==NULL)break;
    trim_back(buffer);
    if(strlen(buffer)>=8 && strncmp(buffer, "LOADFILE", 8)==0){
      if(fgets(filename, 255, stream_from)==NULL)break;
      if(GetLoadfileinfo(stream_to,filename)==0){
        fprintf(stream_to, "%s\n", buffer);
        fprintf(stream_to, "%s\n", filename);
      }
    }
    else if(strlen(buffer)>=9&&strncmp(buffer, "LOADVFILE", 9)==0){
      if(fgets(filename, 255, stream_from)==NULL)break;
      if(GetLoadvfileinfo(stream_to, filename)==0){
        fprintf(stream_to, "%s\n", buffer);
        fprintf(stream_to, "%s\n", filename);
      }
    }
    else{
      fprintf(stream_to, "%s\n", buffer);
    }
  }
  fclose(stream_from);
  fclose(stream_to);
  if(outeqin == 1){
    char *tofile;

    tofile = ssf_from;
    copyfile(".",tempfile,tofile,REPLACE_FILE);
    UNLINK(tempfile);
  }
}

  /* ------------------ UpdateTimes ------------------------ */

void UpdateTimes(void){
  int ntimes2;
  float *timescopy;
  int i;
  float dt_MIN=100000.0;

  FREEMEMORY(geominfoptrs);
  ngeominfoptrs=0;
  GetGeomInfoPtrs(&geominfoptrs,&ngeominfoptrs);

  // pass 1 - determine ntimes

  UpdateShow();
  CheckMemory;
  nglobal_times = 0;

  for(i=0;i<ngeominfoptrs;i++){
    geomdata *geomi;

    geomi = geominfoptrs[i];
    if(geomi->loaded==0||geomi->display==0)continue;
    nglobal_times+=geomi->ntimes;
  }
  if(visTerrainType!=TERRAIN_HIDDEN){
    for(i=0;i<nterraininfo;i++){
      terraindata *terri;

      terri = terraininfo + i;
      if(terri->loaded==1)nglobal_times+=terri->ntimes;
    }
  }
  if(visShooter!=0&&shooter_active==1){
    nglobal_times+=nshooter_frames;
  }
  for(i=0;i<ntours;i++){
    tourdata *touri;

    touri = tourinfo + i;
    if(touri->display==0)continue;
    nglobal_times+= touri->ntimes;
  }
  for(i=0;i<npartinfo;i++){
    partdata *parti;

    parti = partinfo + i;
    if(parti->loaded==0)continue;
    nglobal_times+= parti->ntimes;
  }
  for(i=0;i<nsliceinfo;i++){
    slicedata *sd;

    sd=sliceinfo+i;
    if(sd->loaded==1||sd->vloaded==1){
      nglobal_times+=sd->ntimes;
    }
  }
  for(i=0;i<npatchinfo;i++){
    patchdata *patchi;

    patchi = patchinfo + i;
    if(patchi->loaded==1&&patchi->filetype==PATCH_GEOMETRY){
      nglobal_times+=patchi->ngeom_times;
    }
  }
  for(i=0;i<nmeshes;i++){
    patchdata *patchi;
    meshdata *meshi;
    int filenum;

    meshi=meshinfo+i;
    filenum =meshi->patchfilenum;
    if(filenum!=-1){
      patchi=patchinfo+filenum;
      if(patchi->loaded==1&&patchi->filetype!=PATCH_GEOMETRY){
        nglobal_times+=meshi->npatch_times;
      }
    }
  }
  if(ReadZoneFile==1&&visZone==1){
    nglobal_times+=nzone_times;
  }
  if(ReadIsoFile==1&&visAIso!=0){
    for(i=0;i<nmeshes;i++){
      meshdata *meshi;
      isodata *ib;

      meshi=meshinfo+i;
      if(meshi->isofilenum<0)continue;
      ib = isoinfo + meshi->isofilenum;
      if(ib->geomflag==1||ib->loaded==0)continue;
      nglobal_times+=meshi->niso_times;
    }
  }
  if(nvolrenderinfo>0){
    for(i=0;i<nmeshes;i++){
      volrenderdata *vr;
      meshdata *meshi;

      meshi=meshinfo+i;
      vr = &meshi->volrenderinfo;
      if(vr->fireslice==NULL||vr->smokeslice==NULL)continue;
      if(vr->loaded==0||vr->display==0)continue;
      nglobal_times+=vr->ntimes;
    }
  }
  {
    smoke3ddata *smoke3di;

    if(Read3DSmoke3DFile==1&&vis3DSmoke3D==1){
      for(i=0;i<nsmoke3dinfo;i++){
        smoke3di = smoke3dinfo + i;
        if(smoke3di->loaded==0)continue;
        nglobal_times+= smoke3di->ntimes;
      }
    }
  }

  // end pass 1

  CheckMemory;
  FREEMEMORY(global_times);
  if(nglobal_times>0)NewMemory((void **)&global_times,nglobal_times*sizeof(float));
  timescopy=global_times;

  // pass 2 - merge times arrays

  for(i=0;i<ngeominfoptrs;i++){
    geomdata *geomi;
    int n;

    geomi = geominfoptrs[i];
    if(geomi->loaded==0||geomi->display==0)continue;
    for(n=0;n<geomi->ntimes;n++){
      float t_diff;

      *timescopy++=geomi->times[n];
      t_diff = timescopy[-1]-timescopy[-2];
      if(n>0&&t_diff<dt_MIN&&t_diff>0.0){
        dt_MIN=t_diff;
      }
    }
  }
  if(visTerrainType!=TERRAIN_HIDDEN){
    for(i=0;i<nterraininfo;i++){
      terraindata *terri;
      int n;

      terri = terraininfo + i;
      if(terri->loaded==0)continue;
      for(n=0;n<terri->ntimes;n++){
        float t_diff;

        *timescopy++=terri->times[n];
        t_diff = timescopy[-1]-timescopy[-2];
        if(n>0&&t_diff<dt_MIN&&t_diff>0.0){
          dt_MIN=t_diff;
        }
      }
    }
  }
  if(visShooter!=0&&shooter_active==1){
    for(i=0;i<nshooter_frames;i++){
      float t_diff;

      *timescopy++=shoottimeinfo[i].time;

      t_diff = timescopy[-1]-timescopy[-2];
      if(i>0&&t_diff<dt_MIN&&t_diff>0.0){
        dt_MIN=t_diff;
      }
    }
    CheckMemory;
  }

  for(i=0;i<ntours;i++){
    tourdata *touri;
    int n;

    touri = tourinfo + i;
    if(touri->display==0)continue;
    for(n=0;n<touri->ntimes;n++){
      float t_diff;

      *timescopy++=touri->path_times[n];
      t_diff = timescopy[-1]-timescopy[-2];
      if(n>0&&t_diff<dt_MIN&&t_diff>0.0){
        dt_MIN=t_diff;
      }
    }
  }

  tmax_part=0.0;
  for(i=0;i<npartinfo;i++){
    partdata *parti;
    int n;

    parti = partinfo + i;
    if(parti->loaded==0)continue;
    for(n=0;n<parti->ntimes;n++){
      float t_diff;

      *timescopy++=parti->times[n];
      t_diff = timescopy[-1]-timescopy[-2];
      if(n>0&&t_diff<dt_MIN&&t_diff>0.0){
        dt_MIN=t_diff;
      }
    }
    tmax_part=MAX(parti->times[parti->ntimes-1],tmax_part);
  }

  for(i=0;i<nsliceinfo;i++){
    slicedata *sd;
    int n;

    sd = sliceinfo + i;
    if(sd->loaded==1||sd->vloaded==1){
      for(n=0;n<sd->ntimes;n++){
        float t_diff;

        *timescopy++=sd->times[n];
        t_diff = timescopy[-1]-timescopy[-2];
        if(n>0&&t_diff<dt_MIN&&t_diff>0.0){
          dt_MIN=t_diff;
        }
      }
    }
  }

  for(i=0;i<npatchinfo;i++){
    patchdata *patchi;
    int n;

    patchi = patchinfo + i;
    if(patchi->loaded==1&&patchi->filetype==PATCH_GEOMETRY){
      for(n=0;n<patchi->ngeom_times;n++){
        float t_diff;

        *timescopy++=patchi->geom_times[n];
        t_diff = timescopy[-1]-timescopy[-2];
        if(n>0&&t_diff<dt_MIN&&t_diff>0.0){
          dt_MIN=t_diff;
        }
      }
    }
  }
  for(i=0;i<nmeshes;i++){
    patchdata *patchi;
    meshdata *meshi;
    int filenum;

    meshi=meshinfo + i;
    filenum=meshi->patchfilenum;
    if(filenum!=-1){
      patchi = patchinfo + filenum;
      if(patchi->loaded==1&&patchi->filetype!=PATCH_GEOMETRY){
        int n;

        for(n=0;n<meshi->npatch_times;n++){
          float t_diff;

          *timescopy++=meshi->patch_times[n];
          t_diff = timescopy[-1]-timescopy[-2];
          if(n>0&&t_diff<dt_MIN&&t_diff>0.0){
            dt_MIN=t_diff;
          }
        }
      }
    }
  }
  if(nvolrenderinfo>0){
    for(i=0;i<nmeshes;i++){
      volrenderdata *vr;
      meshdata *meshi;
      int n;

      meshi=meshinfo + i;
      vr = &meshi->volrenderinfo;
      if(vr->smokeslice==NULL)continue;
      if(vr->loaded==0||vr->display==0)continue;
      for(n=0;n<vr->ntimes;n++){
        float t_diff;

        *timescopy++=vr->times[n];
        if(n>0){
          t_diff = timescopy[-1]-timescopy[-2];
          if(t_diff<dt_MIN&&t_diff>0.0){
            dt_MIN=t_diff;
          }
        }
      }
    }
  }
  if(ReadZoneFile==1&&visZone==1){
    int n;

    for(n=0;n<nzone_times;n++){
      float t_diff;

      *timescopy++=zone_times[n];
      t_diff = timescopy[-1]-timescopy[-2];
      if(n>0&&t_diff<dt_MIN&&t_diff>0.0){
        dt_MIN=t_diff;
      }
    }
  }
  if(ReadIsoFile==1&&visAIso!=0){
    for(i=0;i<nisoinfo;i++){
      meshdata *meshi;
      isodata *ib;
      int n;

      ib = isoinfo+i;
      if(ib->geomflag==1||ib->loaded==0)continue;
      meshi=meshinfo + ib->blocknumber;
      for(n=0;n<meshi->niso_times;n++){
        float t_diff;

        *timescopy++=meshi->iso_times[n];
        t_diff = timescopy[-1]-timescopy[-2];
        if(n>0&&t_diff<dt_MIN&&t_diff>0.0){
          dt_MIN=t_diff;
        }
      }
    }
  }
  {
    smoke3ddata *smoke3di;

    if(Read3DSmoke3DFile==1&&vis3DSmoke3D==1){
      for(i=0;i<nsmoke3dinfo;i++){
        int n;

        smoke3di = smoke3dinfo + i;
        if(smoke3di->loaded==0)continue;
        for(n=0;n<smoke3di->ntimes;n++){
          float t_diff;

          *timescopy++=smoke3di->times[n];
          t_diff = timescopy[-1]-timescopy[-2];
          if(n>0&&t_diff<dt_MIN&&t_diff>0.0){
            dt_MIN=t_diff;
          }
        }
      }
    }
  }

  // end pass 2

  // sort times array and remove duplicates

  if(nglobal_times>0)qsort( (float *)global_times, (size_t)nglobal_times, sizeof( float ), CompareFloat );

  {
    int n,n2;

    for(n2=1,ntimes2=nglobal_times,n=1;n<nglobal_times;n++){
      if(ABS(global_times[n]-global_times[n-1])>dt_MIN/10.0){
        global_times[n2]=global_times[n];
        n2++;
      }
      else{
        ntimes2--;
      }
    }
  }
  nglobal_times=ntimes2;

  // pass 3 - allocate memory for individual times array

  if(nglobal_times>ntimes_old){
    ntimes_old=nglobal_times;
    FREEMEMORY(render_frame);
    if(nglobal_times>0)NewMemory((void **)&render_frame,nglobal_times*sizeof(int));
  }

  for(i=0;i<ngeominfoptrs;i++){
    geomdata *geomi;

    geomi = geominfoptrs[i];
    if(geomi->loaded==0||geomi->display==0)continue;
    FREEMEMORY(geomi->timeslist);
    if(nglobal_times>0)NewMemory((void **)&geomi->timeslist,nglobal_times*sizeof(int));
  }
  for(i=0;i<npartinfo;i++){
    partdata *parti;

    parti=partinfo+i;
    FREEMEMORY(parti->timeslist);
    if(nglobal_times>0)NewMemory((void **)&parti->timeslist,nglobal_times*sizeof(int));
  }
  for(i=0;i<ntours;i++){
    tourdata *touri;

    touri=tourinfo + i;
    if(touri->display==0)continue;
    FREEMEMORY(touri->timeslist);
    if(nglobal_times>0)NewMemory((void **)&touri->timeslist,nglobal_times*sizeof(int));
  }
  if(visTerrainType!=TERRAIN_HIDDEN){
    for(i=0;i<nterraininfo;i++){
      terraindata *terri;

      terri = terraininfo + i;
      if(terri->loaded==0)continue;
      FREEMEMORY(terri->timeslist);
      if(nglobal_times>0)NewMemory((void **)&terri->timeslist,nglobal_times*sizeof(int));
    }
  }
  if(hrrinfo!=NULL){
    FREEMEMORY(hrrinfo->timeslist);
    FREEMEMORY(hrrinfo->times);
    FREEMEMORY(hrrinfo->hrrval);
    if(hrrinfo->loaded==1&&hrrinfo->display==1&&nglobal_times>0){
      int jstart=0;

      NewMemory((void **)&hrrinfo->timeslist,nglobal_times*sizeof(int));
      NewMemory((void **)&hrrinfo->times,nglobal_times*sizeof(float));
      NewMemory((void **)&hrrinfo->hrrval,nglobal_times*sizeof(float));
      hrrinfo->ntimes=nglobal_times;
      for(i=0;i<nglobal_times;i++){
        int j, foundit;

        foundit=0;
        hrrinfo->times[i]=global_times[i];
        for(j=jstart;j<hrrinfo->ntimes_csv-1;j++){
          if(hrrinfo->times_csv[j]<=global_times[i]&&global_times[i]<hrrinfo->times_csv[j+1]){
            float f1, tbot;

            foundit=1;
            tbot = hrrinfo->times_csv[j+1]-hrrinfo->times_csv[j];
            if(tbot>0.0){
              f1 = (global_times[i]-hrrinfo->times_csv[j])/tbot;
            }
            else{
              f1=0.0;
            }
            hrrinfo->hrrval[i]=(1.0-f1)*hrrinfo->hrrval_csv[j]+f1*hrrinfo->hrrval_csv[j+1];
            jstart=j;
            break;
          }
        }
        if(foundit==0){
          hrrinfo->hrrval[i]=hrrinfo->hrrval_csv[hrrinfo->ntimes_csv-1];
        }
      }
    }
  }
  FREEMEMORY(shooter_timeslist);
  if(visShooter!=0&&shooter_active==1){
    NewMemory((void **)&shooter_timeslist,nshooter_frames*sizeof(int));
  }

  for(i=0;i<nsliceinfo;i++){
    slicedata *sd;

    sd = sliceinfo + i;
    FREEMEMORY(sd->timeslist);
    if(nglobal_times>0)NewMemory((void **)&sd->timeslist,nglobal_times*sizeof(int));
  }
  if(nvolrenderinfo>0){
    for(i=0;i<nmeshes;i++){
      meshdata *meshi;
      volrenderdata *vr;

      meshi = meshinfo + i;
      vr = &(meshi->volrenderinfo);
      if(vr->fireslice==NULL||vr->smokeslice==NULL)continue;
      if(vr->loaded==0||vr->display==0)continue;
      FREEMEMORY(vr->timeslist);
      if(nglobal_times>0)NewMemory((void **)&vr->timeslist,nglobal_times*sizeof(int));
    }
  }
  {
    smoke3ddata *smoke3di;

    for(i=0;i<nsmoke3dinfo;i++){
      smoke3di = smoke3dinfo + i;
      FREEMEMORY(smoke3di->timeslist);
      if(nglobal_times>0)NewMemory((void **)&smoke3di->timeslist,nglobal_times*sizeof(int));
    }
  }
  for(i=0;i<nmeshes;i++){
    meshdata *meshi;

    meshi=meshinfo+i;
    if(meshi->iso_times==NULL)continue;
    FREEMEMORY(meshi->iso_timeslist);
    if(nglobal_times>0)NewMemory((void **)&meshi->iso_timeslist,  nglobal_times*sizeof(int));
  }

  for(i=0;i<npatchinfo;i++){
    patchdata *patchi;

    patchi = patchinfo + i;
    FREEMEMORY(patchi->geom_timeslist);
    if(patchi->filetype!=PATCH_GEOMETRY)continue;
    if(patchi->geom_times==NULL)continue;
    if(nglobal_times>0)NewMemory((void **)&patchi->geom_timeslist,nglobal_times*sizeof(int));
  }
  for(i=0;i<nmeshes;i++){
    FREEMEMORY(meshinfo[i].patch_timeslist);
  }
  for(i=0;i<nmeshes;i++){
    if(meshinfo[i].patch_times==NULL)continue;
    if(nglobal_times>0)NewMemory((void **)&meshinfo[i].patch_timeslist,nglobal_times*sizeof(int));
  }

  FREEMEMORY(zone_timeslist);
  if(nglobal_times>0)NewMemory((void **)&zone_timeslist,     nglobal_times*sizeof(int));

  FREEMEMORY(targtimeslist);
  if(nglobal_times>0)NewMemory((void **)&targtimeslist,  nglobal_times*sizeof(int));

  // end pass 3

  // reset render_frame array

  if(current_script_command!=NULL&&
    (current_script_command->command==SCRIPT_VOLSMOKERENDERALL||current_script_command->command==SCRIPT_ISORENDERALL)
    ){
    if(current_script_command->first==1){
      int n;

      for(n=0;n<nglobal_times;n++){
        render_frame[n]=0;
      }
      current_script_command->first=0;
    }
  }
  else{
    int n;

    for(n=0;n<nglobal_times;n++){
      render_frame[n]=0;
    }
  }

  // reallocate times array

  if(nglobal_times==0){
    FREEMEMORY(global_times);
  }
  if(nglobal_times>0)ResizeMemory((void **)&global_times,nglobal_times*sizeof(float));

  // pass 4 - initialize individual time pointers

  izone=0;
  ResetItimes0();
  for(i=0;i<ngeominfoptrs;i++){
    geomdata *geomi;

    geomi = geominfoptrs[i];
    if(geomi->loaded==0||geomi->display==0)continue;
    geomi->itime=0;
  }
  for(i=0;i<nmeshes;i++){
    meshdata *meshi;

    meshi=meshinfo+i;
    meshi->patch_itime=0;
  }
  for(i=0;i<nsliceinfo;i++){
    slicedata *sd;

    sd = sliceinfo + i;
    sd->itime=0;
  }
  frame_index=first_frame_index;
  for(i=0;i<nmeshes;i++){
    meshdata *meshi;

    meshi=meshinfo+i;
    if(meshi->iso_times==NULL)continue;
    meshi->iso_itime=0;
  }
  for(i=0;i<npartinfo;i++){
    partdata *parti;

    parti = partinfo + i;
    parti->itime=0;
  }

  /* determine visibility of each blockage at each time step */

  for(i=0;i<nmeshes;i++){
    int j;
    meshdata *meshi;

    meshi=meshinfo+i;
    for(j=0;j<meshi->nbptrs;j++){
      blockagedata *bc;

      bc = meshi->blockageinfoptrs[j];
      if(bc->showtime==NULL)continue;
      FREEMEMORY(bc->showtimelist);
      if(nglobal_times>0){
        int k;

        NewMemory((void **)&bc->showtimelist,nglobal_times*sizeof(int));
        for(k=0;k<nglobal_times;k++){
          int listindex;

          bc->showtimelist[k]=1;
          listindex=GetIndex(global_times[k],bc->showtime,bc->nshowtime);
          bc->showtimelist[k]=bc->showhide[listindex];
        }
      }
    }
  }

  /* determine state of each device at each time step */

  for(i=0;i<ndeviceinfo;i++){
    devicedata *devicei;

    devicei = deviceinfo + i;
    if(devicei->object->visible==0)continue;
    if(devicei->nstate_changes==0)continue;
    FREEMEMORY(devicei->showstatelist);
    if(nglobal_times>0){
      int k;

      NewMemory((void **)&devicei->showstatelist,nglobal_times*sizeof(int));
      for(k=0;k<nglobal_times;k++){
        int listindex;

        listindex=GetIndex(global_times[k],devicei->act_times,devicei->nstate_changes);
        devicei->showstatelist[k]=devicei->state_values[listindex];
      }
    }
  }

  /* determine visibility of each vent at each time step */

  for(i=0;i<nmeshes;i++){
    int j;
    meshdata *meshi;

    meshi=meshinfo+i;
    if(meshi->ventinfo==NULL)continue;
    for(j=0;j<meshi->nvents;j++){
      ventdata *vi;

      vi = meshi->ventinfo+j;
      if(vi->showtime==NULL)continue;
      FREEMEMORY(vi->showtimelist);
      if(nglobal_times>0){
        int k;

        NewMemory((void **)&vi->showtimelist,nglobal_times*sizeof(int));
        for(k=0;k<nglobal_times;k++){
          int listindex;

          vi->showtimelist[k]=1;
          listindex=GetIndex(global_times[k],vi->showtime,vi->nshowtime);
          vi->showtimelist[k]=vi->showhide[listindex];
        }
      }
    }
  }

  /* determine visibility of each circular vent at each time step */

  for (i = 0; i<nmeshes; i++) {
    int j;
    meshdata *meshi;

    meshi = meshinfo + i;
    if(meshi->cventinfo == NULL)continue;
    for (j = 0; j<meshi->ncvents; j++) {
      cventdata *cvi;

      cvi = meshi->cventinfo + j;
      if(cvi->showtime == NULL)continue;
      FREEMEMORY(cvi->showtimelist);
      if(nglobal_times>0) {
        int k;

        NewMemory((void **)&cvi->showtimelist, nglobal_times * sizeof(int));
        for (k = 0; k<nglobal_times; k++) {
          int listindex;

          cvi->showtimelist[k] = 1;
          listindex = GetIndex(global_times[k], cvi->showtime, cvi->nshowtime);
          cvi->showtimelist[k] = cvi->showhide[listindex];
        }
      }
    }
  }

  if(nglobal_times>0)SynchTimes();
  updatefaces=1;
  if(nglobal_times>0){
    UpdateTimeLabels();
    updateGluiTimeBounds(global_times[0],global_times[nglobal_times-1]);
  }
  show_slice_terrain=0;
  if(visTerrainType==TERRAIN_3D_MAP){
    for(i=0;i<nsliceinfo;i++){
      slicedata *sd;

      sd = sliceinfo + i;
      if(sd->loaded==0||sd->display==0||sd->slicetype!=SLICE_TERRAIN)continue;
      show_slice_terrain=1;
      break;
    }
  }
}

/* ------------------ GetPlotState ------------------------ */

int GetPlotState(int choice){
  int i;

  update_loaded_lists();
  switch(choice){
    case STATIC_PLOTS:
    case STATIC_PLOTS_NORECURSE:
      stept = 0;
      for(i=0;i<nmeshes;i++){
        plot3ddata *ploti;
        meshdata *meshi;

        meshi=meshinfo + i;
        if(meshi->plot3dfilenum==-1)continue;
        ploti = plot3dinfo + meshi->plot3dfilenum;
        if(ploti->loaded==0||ploti->display==0)continue;
        if(visx_all==0&&visy_all==0&&visz_all==0&&visiso==0)continue;
        return STATIC_PLOTS;
      }
      if(choice!=STATIC_PLOTS_NORECURSE){
        return GetPlotState(DYNAMIC_PLOTS_NORECURSE);
      }
      break;
    case DYNAMIC_PLOTS:
    case DYNAMIC_PLOTS_NORECURSE:
      for(i=0;i<nslice_loaded;i++){
        slicedata *slicei;

        slicei = sliceinfo + slice_loaded_list[i];
        if(slicei->display==0||slicei->type!=islicetype)continue;
        stept = 1;
        return DYNAMIC_PLOTS;
      }
      if(visGrid==0)stept = 1;
      if(visTerrainType!=TERRAIN_HIDDEN){
        for(i=0;i<nterraininfo;i++){
          terraindata *terri;

          terri = terraininfo + i;
          if(terri->loaded==1){
            return DYNAMIC_PLOTS;
          }
        }
      }
      for(i=0;i<nvsliceinfo;i++){
        vslicedata *vslicei;

        vslicei = vsliceinfo + i;
        if(vslicei->display==0||vslicei->type!=islicetype)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<npatch_loaded;i++){
        patchdata *patchi;

        patchi = patchinfo + patch_loaded_list[i];
        if(patchi->display==0||patchi->type!=ipatchtype)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<npartinfo;i++){
        partdata *parti;

        parti = partinfo + i;
        if(parti->loaded==0||parti->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<nisoinfo;i++){
        isodata *isoi;

        isoi = isoinfo + i;
        if(isoi->loaded==0)continue;
        if(isoi->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<nzoneinfo;i++){
        zonedata *zonei;

        zonei = zoneinfo + i;
        if(zonei->loaded==0||zonei->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<ntours;i++){
        tourdata *touri;

        touri = tourinfo + i;
        if(touri->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<nsmoke3dinfo;i++){
        smoke3ddata *smoke3di;

        smoke3di = smoke3dinfo + i;
        if(smoke3di->loaded==0||smoke3di->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      if(nvolrenderinfo>0){
        for(i=0;i<nmeshes;i++){
          meshdata *meshi;
          volrenderdata *vr;

          meshi = meshinfo + i;
          vr = &(meshi->volrenderinfo);
          if(vr->fireslice==NULL||vr->smokeslice==NULL)continue;
          if(vr->loaded==0||vr->display==0)continue;
          return DYNAMIC_PLOTS;
        }
      }
      if(visShooter!=0&&shooter_active==1){
        return DYNAMIC_PLOTS;
      }
      if(choice!=DYNAMIC_PLOTS_NORECURSE)return GetPlotState(STATIC_PLOTS_NORECURSE);
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
  stept = 0;
  return NO_PLOTS;
}

/* ------------------ GetIndex ------------------------ */

int GetIndex(float key, const float *list, int nlist){
  int i;

  if(nlist==1)return 0;
  if(key<list[1])return 0;
  if(key>=list[nlist-1])return nlist-1;
  for(i=1;i<nlist-1;i++){
    if(list[i]<=key&&key<list[i+1])return i;
  }
  return 0;
}

/* ------------------ ISearch ------------------------ */

int ISearch(float *list, int nlist, float key, int guess){
  /*
     find val such that list[val]<=key<list[val+1]
     start with val=guess
  */

  int low, mid, high;

  if(nlist<=2||key<list[0])return 0;
  if(key>=list[nlist-2])return nlist-2;
  if(guess<0)guess=0;
  if(guess>nlist-2)guess=nlist-2;
  if(list[guess]<=key&&key<list[guess+1])return guess;

  low = 0;
  high = nlist - 1;
  while(high-low>1){
    mid = (low+high)/2;
    if(list[mid]>key){
      high=mid;
    }
    else{
      low=mid;
    }
  }
  if(list[high]==key)return high;
  return low;
}

/* ------------------ ResetItimes0 ------------------------ */

void ResetItimes0(void){
  if(current_script_command==NULL||current_script_command->command!=SCRIPT_VOLSMOKERENDERALL||current_script_command->command!=SCRIPT_ISORENDERALL){
    itimes=first_frame_index;
  }
}

/* ------------------ UpdateClipbounds ------------------------ */

void UpdateClipbounds(int set_i0, int *i0, int set_i1, int *i1, int imax){

  if(set_i0==0&&set_i1==0)return;
  if(set_i0==1&&set_i1==1){
    if(*i0>imax-1){*i0=imax-1; *i1=imax;}
    if(*i1>imax)*i1=imax;
    if(*i1<1){*i1=1;*i0=0;}
    if(*i0<0)*i0=0;
    if(*i0>=*i1){*i0=*i1-1;}
  }
  if(set_i0==1&&set_i1==0){
    if(*i0<0)*i0=0;
    if(*i0>imax)*i0=imax;
  }
  if(set_i0==0&&set_i1==1){
    if(*i1<0)*i1=0;
    if(*i1>imax)*i1=imax;
  }
}

/* ------------------ UpdateColorTable ------------------------ */

void UpdateColorTable(colortabledata *ctableinfo, int nctableinfo){
  int i, ncolortableinfo_old;

  if(nctableinfo<=0)return;

  ncolortableinfo_old=ncolortableinfo;
  ResizeMemory((void **)&colortableinfo, (ncolortableinfo+nctableinfo)*sizeof(colortabledata));
  for(i = 0; i<nctableinfo; i++){
    colortabledata *newentryi, *fromi;

    fromi = ctableinfo+i;
    newentryi=get_colortable(fromi->label);
    if(newentryi==NULL){
      newentryi = colortableinfo + ncolortableinfo;
      ncolortableinfo++;
    }
    newentryi->color[0] = fromi->color[0];
    newentryi->color[1] = fromi->color[1];
    newentryi->color[2] = fromi->color[2];
    newentryi->color[3] = fromi->color[3];
    strcpy(newentryi->label, fromi->label);
  }
  ResizeMemory((void **)&colortableinfo, ncolortableinfo*sizeof(colortabledata));
  UpdateColorTableList(ncolortableinfo_old);
}

/* ------------------ UpdateShowScene ------------------------ */

void UpdateShowScene(void){
  if(update_playmovie==1){
    enable_disable_playmovie();
    update_playmovie = 0;
  }
  update_render_start_button();
  if(update_makemovie == 1)MakeMovie();
  if(compute_fed == 1)DefineAllFEDs();
  if(restart_time == 1){
    restart_time = 0;
    ResetItimes0();
  }
  if(loadfiles_at_startup==1&&update_load_Files == 1){
    load_Files();
  }
  if(update_startup_view == 1){
    cameradata *ca;

    ca = GetCamera(label_startup_view);
    if(ca != NULL){
      ResetMenu(ca->view_id);
    }
    update_rotation_center = 0;
    update_rotation_center_ini = 0;
    update_startup_view = 0;
  }
  if(update_tourlist == 1){
    Update_Tourlist();
  }
  if(update_gslice == 1){
    update_gslice_parms();
  }
#define MESH_LIST 4
  if(update_rotation_center == 1){
    camera_current->rotation_index = glui_rotation_index;
    Motion_CB(MESH_LIST);
    update_rotation_center = 0;
  }
  if(update_rotation_center_ini == 1){
    camera_current->rotation_index = glui_rotation_index_ini;
    Motion_CB(MESH_LIST);
    update_rotation_center_ini = 0;
  }
  if(camera_current->dirty == 1){
    UpdateCamera(camera_current);
  }
  if(updateclipvals == 1){
    Clip2Cam(camera_current);
    update_clip_all();
    updateclipvals = 0;
  }
  if(update_selectedtour_index == 1){
    update_tourindex();
  }
  if(trainer_mode == 1 && fontindex != LARGE_FONT)FontMenu(LARGE_FONT);
  if(updateindexcolors == 1){
    UpdateIndexColors();
  }
  if(force_isometric == 1){
    force_isometric = 0;
    projection_type = 1;
    camera_current->projection_type = projection_type;
    ZoomMenu(UPDATE_PROJECTION);
  }
  if(convert_ini == 1){
    writeini(SCRIPT_INI, ini_to);
    exit(0);
  }
  if(convert_ssf==1||update_ssf==1){
    ConvertSsf();
    exit(0);
  }
  UpdateShow();
  if(global_times!=NULL&&updateUpdateFrameRateMenu==1)FrameRateMenu(frameratevalue);
  if(updatefaces==1)UpdateFaces();
  if(updatefacelists==1)UpdateFacelists();
}

/* ------------------ UpdateDisplay ------------------------ */
#define TERRAIN_FIRE_LINE_UPDATE 39

void UpdateDisplay(void){

  LOCK_IBLANK
  if(update_setvents==1){
    SetVentDirs();
    update_setvents=0;
  }
  UNLOCK_IBLANK
  if(update_have_gvec == 1){
    update_have_gvec = 0;
    update_gvec_down(1);
  }
  if(update_smokecolorbar == 1){
    update_smokecolorbar = 0;
    SmokeColorbarMenu(fire_colorbar_index);
  }
  if(update_colorbartype == 1){
    colorbardata *cb;

    cb = GetColorbar(colorbarname);
    if(cb != NULL){
      colorbartype = cb - colorbarinfo;
      UpdateCurrentColorbar(cb);
      if(colorbartype != colorbartype_default){
        colorbartype_ini = colorbartype;
      }
    }
    update_colorbartype = 0;
  }
  if(update_fire_line == 1){
    WUI_CB(TERRAIN_FIRE_LINE_UPDATE);
    update_fire_line = 0;
  }
  if(updatezoommenu == 1 || first_display > 0){
    if(first_display > 0)first_display--;
    updatezoommenu = 0;
    ZoomMenu(zoomindex);
  }
  if(update_makeiblank_smoke3d == 1){
    makeiblank_smoke3d();
  }
#ifdef pp_CULL
  if(update_initcull == 1)initcull(cullsmoke);
#endif
  if(update_streaks == 1 && ReadPartFile == 1){
    ParticleStreakShowMenu(streak_index);
    update_streaks = 0;
  }
  if(update_screensize == 1){
    update_screensize = 0;
    update_windowsizelist();
    ResizeWindow(screenWidthINI, screenHeightINI);
  }
  if(updatemenu == 1 && usemenu == 1 && menustatus == GLUT_MENU_NOT_IN_USE){
    glutDetachMenu(GLUT_RIGHT_BUTTON);
    InitMenus(LOAD);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    updatemenu = 0;
  }
  if(update_fire_colorbar_index == 1){
    SmokeColorbarMenu(fire_colorbar_index_ini);
    update_fire_colorbar_index = 0;
  }
  if(update_colorbar_select_index == 1 && colorbar_select_index >= 0 && colorbar_select_index <= 255){
    update_colorbar_select_index = 0;
    UpdateRGBColors(colorbar_select_index);
  }
}
