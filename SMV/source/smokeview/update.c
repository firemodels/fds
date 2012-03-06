// $Date$ 
// $Revision$
// $Author$

// svn revision character string
char update_revision[]="$Revision$";

#define IN_UPDATE
#include "options.h"
#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include "string_util.h"
#include "update.h"
#include "smokeviewvars.h"
#include "compress.h"

/* ------------------ compare ------------------------ */

int compare_float( const void *arg1, const void *arg2 ){
  float x, y;
  x=*(float *)arg1;
  y=*(float *)arg2;
  if( x< y)return -1;
  if( x> y)return 1;
  return 0;
}

/* ------------------ update_framenumber ------------------------ */

void update_framenumber(int changetime){
  int i,ii;
//  int redisplay;

  if(force_redisplay==1||(itimeold!=itimes&&changetime==1)){
    force_redisplay=0;
    itimeold=itimes;
    if(showsmoke==1||showevac==1){
      for(i=0;i<npartinfo;i++){
        particle *parti;

        parti = partinfo+i;
        if(parti->loaded==1){
          if(parti->timeslist==NULL)continue;
          parti->iframe=parti->timeslist[itimes];
        }
      }
    }
    if(hrrinfo!=NULL&&hrrinfo->loaded==1&&hrrinfo->display==1&&hrrinfo->timeslist!=NULL){
      hrrinfo->itime=hrrinfo->timeslist[itimes];
    }
    if(showvolrender==1){
      int imesh;

      for(imesh=0;imesh<nmeshes;imesh++){
        mesh *meshi;
        volrenderdata *vr;
        slice *fire, *smoke;
        int j;

        meshi = meshinfo + imesh;
        vr = &(meshi->volrenderinfo);
        fire=vr->fire;
        smoke=vr->smoke;
        vr->smokedataptr=NULL;
        vr->firedataptr=NULL;
        if(fire==NULL||smoke==NULL)continue;
        if(vr->loaded==0||vr->display==0)continue;
        vr->iframe = vr->timeslist[itimes];
        for(j=vr->iframe;j>=0;j--){
          if(vr->dataready[j]==1)break;
        }
        vr->iframe=j;
        if(smoke!=NULL&&vr->iframe>=0){
          if(vr->is_compressed==1||load_volcompressed==1){
            unsigned char *c_smokedata_compressed;
            uLongf framesize;
            float timeval;

            c_smokedata_compressed = vr->smokedataptrs[vr->iframe];
            framesize = smoke->nslicei*smoke->nslicej*smoke->nslicek;
            uncompress_volsliceframe(c_smokedata_compressed,
                           vr->smokedata_view, framesize, &timeval,
                           vr->c_smokedata_view);

            vr->smokedataptr = vr->smokedata_view;
          }
          else{
            vr->smokedataptr = vr->smokedataptrs[vr->iframe];
          }
          CheckMemory;
        }
        if(fire!=NULL&&vr->iframe>=0){
          if(vr->is_compressed==1||load_volcompressed==1){
            unsigned char *c_firedata_compressed;
            uLongf framesize;
            float timeval;

            c_firedata_compressed = vr->firedataptrs[vr->iframe];
            framesize = fire->nslicei*fire->nslicej*fire->nslicek;
            uncompress_volsliceframe(c_firedata_compressed,
                           vr->firedata_view, framesize, &timeval,
                           vr->c_firedata_view);

            vr->firedataptr = vr->firedata_view;
            CheckMemory;
          }
          else{
            vr->firedataptr = vr->firedataptrs[vr->iframe];
          }
          CheckMemory;
        }
      }
    }
    for(i=0;i<ngeominfoptrs;i++){
      geomdata *geomi;

      geomi = geominfoptrs[i];
      if(geomi->loaded==0)continue;
      if(geomi->timeslist==NULL)continue;
      geomi->iframe=geomi->timeslist[itimes];
    }
    if(showslice==1||showvslice==1){
      for(ii=0;ii<nslice_loaded;ii++){
        slice *sd;

        i = slice_loaded_list[ii];
        sd = sliceinfo+i;
        if(sd->timeslist==NULL)continue;
        sd->islice=sd->timeslist[itimes];
      }
    }
    if(show3dsmoke==1){
      for(i=0;i<nsmoke3dinfo;i++){
        smoke3d *smoke3di;

        smoke3di = smoke3dinfo + i;
        if(smoke3di->loaded==0||smoke3di->display==0)continue;
        smoke3di->iframe=smoke3di->timeslist[itimes];
        if(smoke3di->iframe!=smoke3di->lastiframe){
          smoke3di->lastiframe=smoke3di->iframe;
          updatesmoke3d(smoke3di);
        }
      }
      if(nsmoke3dinfo>0)mergesmoke3dcolors(NULL);
    }
    if(showpatch==1){
      for(i=0;i<npatchinfo;i++){
        patch *patchi;

        patchi = patchinfo + i;
        if(patchi->filetype!=2)continue;
        if(patchi->geom_times==NULL)continue;
        if(patchi->geom_timeslist==NULL)continue;
        patchi->geom_ipatch=patchi->geom_timeslist[itimes];
        patchi->geom_ival_static = patchi->geom_ivals_static[patchi->geom_ipatch];
        patchi->geom_ival_dynamic = patchi->geom_ivals_dynamic[patchi->geom_ipatch];
        patchi->geom_nval_static = patchi->geom_nstatics[patchi->geom_ipatch];
        patchi->geom_nval_dynamic = patchi->geom_ndynamics[patchi->geom_ipatch];
      }
      for(i=0;i<nmeshes;i++){
        patch *patchi;
        mesh *meshi;

        meshi = meshinfo+i;
        patchi=patchinfo + meshi->patchfilenum;
        if(patchi->filetype==2)continue;
        if(meshi->patch_times==NULL)continue;
        if(meshi->patch_timeslist==NULL)continue;
        meshi->ipatch=meshi->patch_timeslist[itimes];
        if(patchi->compression_type==0){
          meshi->ipqqi = meshi->ipqq + meshi->ipatch*meshi->npatchsize;
        }
        else{
#ifdef USE_ZLIB
          uncompress_patchdataframe(meshi,meshi->ipatch);
#endif
        }
      }
    }
    if(showiso==1){
      iso *isoi;
      mesh *meshi;

      CheckMemory;
      for(i=0;i<nisoinfo;i++){
        isoi = isoinfo + i;
        if(isoi->loaded==0)continue;
        meshi = meshinfo + isoi->blocknumber;

        if(meshi->iso_times==NULL)continue;
        if(meshi->iso_timeslist==NULL)continue;
        meshi->iiso=meshi->iso_timeslist[itimes];
      }
    }
    if(ntotal_smooth_blockages>0){
      for(i=0;i<nmeshes;i++){
        smoothblockage *sb;
        mesh *meshi;

        meshi = meshinfo+i;
        if(meshi->showsmoothtimelist!=NULL){
          sb=meshi->showsmoothtimelist[itimes];
          if(sb==NULL)continue;
          meshi->nsmoothblockagecolors=sb->nsmoothblockagecolors;
          meshi->smoothblockagecolors=sb->smoothblockagecolors;
          meshi->blockagesurfaces=sb->smoothblockagesurfaces;
        }
      }
    }
    if(showzone==1){
      izone=zonetlist[itimes];
    }
  }
}

/* ------------------ updateShow ------------------------ */

void updateShow(void){
  int i,evacflag,sliceflag,vsliceflag,partflag,patchflag,isoflag,smoke3dflag,tisoflag;
  int slicecolorbarflag;

#ifdef pp_SHOOTER
  int shooter_flag;
#endif
  int ii;
  vslice *vd;
  iso *isoi;
  patch *patchi;
  particle *parti;
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
  show_extreme_below=0;
  show_extreme_above=0;
#ifdef pp_SHOOTER
  showshooter=0;
#endif
  showevac=0;
  showevac_colorbar=0;
  showtarget=0;
  showtitle1=0; showtitle2=0;
  show3dsmoke=0;
  smoke3dflag=0;
  showtour_dialog=0;
  showterrain=0;
  ntitles=0;
  if(visTitle0==1)ntitles++;
  if(strlen(TITLE1)!=0&&visTitle1==1){ntitles++;showtitle1=1;}
  if(strlen(TITLE2)!=0&&visTitle2==1){ntitles++;showtitle2=1;}
  visTitle=0;
  if(visTitle0==1||showtitle1==1||showtitle2==1)visTitle=1;
  visTimeSmoke=1; visTimeSlice=1; visTimePatch=1; visTimeZone=1; visTimeIso=1;

  RenderTime=0;
  if(global_times!=NULL){
    if(settmin_p==1&&global_times[itimes]<tmin_p)visTimeSmoke=0;
    if(settmax_p==1&&global_times[itimes]>tmax_p)visTimeSmoke=0;

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
          showtour_dialog=1;
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
    smoke3d *smoke3di;

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
      mesh *meshi;
      volrenderdata *vr;

      meshi = meshinfo + i;
      vr = &(meshi->volrenderinfo);
      if(vr->fire==NULL||vr->smoke==NULL)continue;
      if(vr->loaded==0||vr->display==0)continue;
      showvolrender=1;
      break;
    }
  }
  sliceflag=0;
  slicecolorbarflag=0;
  if(visTimeSlice==1){
    for(ii=0;ii<nslice_loaded;ii++){
      slice *sd;

      i=slice_loaded_list[ii];
      sd = sliceinfo+i;
      if(sd->display==0||sd->type!=islicetype)continue;
      if(sd->ntimes>0){
        sliceflag=1;
        break;
      }
    }
    for(ii=0;ii<nslice_loaded;ii++){
      mesh *slicemesh;
      slice *sd;

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
    if(show_extreme_above==0){
      for(ii=0;ii<nslice_loaded;ii++){
        slice *sd;

        i=slice_loaded_list[ii];
        sd = sliceinfo+i;
        if(sd->display==0||sd->type!=islicetype)continue;
        if(sd->extreme_max==1){
          show_extreme_above=1;
          break;
        }
      }
    }
    if(show_extreme_below==0){
      for(ii=0;ii<nslice_loaded;ii++){
        slice *sd;

        i=slice_loaded_list[ii];
        sd = sliceinfo+i;
        if(sd->display==0||sd->type!=islicetype)continue;
        if(sd->extreme_min==1){
          show_extreme_below=1;
          break;
        }
      }
    }
  }
  isoflag=0;
  tisoflag=0;
  if(visTimeIso==1){
    for(i=0;i<nisoinfo;i++){
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
    for(i=0;i<nvslice;i++){
      vd = vsliceinfo+i;
      if(vd->loaded==0||vd->display==0)continue;
      if(sliceinfo[vd->ival].type!=islicetype)continue;
      vsliceflag=1;
      break;
    }
    for(i=0;i<nvslice;i++){
      mesh *slicemesh;
      slice *sd;

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
    for(ii=0;ii<npatch_loaded;ii++){
      i = patch_loaded_list[ii];
      patchi=patchinfo+i;
      if(patchi->display==0||patchi->type!=ipatchtype)continue;
      patchflag=1;
      break;
    }
    for(ii=0;ii<npatch_loaded;ii++){
      i = patch_loaded_list[ii];
      patchi=patchinfo+i;
      if(patchi->display==0||patchi->type!=ipatchtype)continue;
      if(patchi->extreme_max==1){
        show_extreme_above=1;
        break;
      }
    }
    for(ii=0;ii<npatch_loaded;ii++){
      i = patch_loaded_list[ii];
      patchi=patchinfo+i;
      if(patchi->display==0||patchi->type!=ipatchtype)continue;
      if(patchi->extreme_min==1){
        show_extreme_below=1;
        break;
      }
    }
    patchembedded=0;
    for(ii=0;ii<npatch_loaded;ii++){
      i = patch_loaded_list[ii];
      patchi=patchinfo+i;
      if(patchi->filetype!=2||patchi->display==0||patchi->type!=ipatchtype)continue;
      patchembedded=1;
    }
  }
  partflag=0;
  if(visSmoke==1&&visTimeSmoke==1){
    for(i=0;i<npartinfo;i++){
      parti = partinfo + i;
      if(parti->evac==1)continue;
      if(parti->loaded==0||parti->display==0)continue;
      partflag=1;
      current_particle_type=parti->particle_type;
      if(current_particle_type!=last_particle_type)updatechopcolors();
      break;
    }
    if(current_property!=NULL){
      if(current_property->extreme_max==1)show_extreme_above=1;
      if(current_property->extreme_min==1)show_extreme_below=1;
    }
  }
  evacflag=0;
  if(visEvac==1&&visTimeEvac==1){
    for(i=0;i<npartinfo;i++){
      parti = partinfo + i;
      if(parti->evac==0)continue;
      if(parti->loaded==0||parti->display==0)continue;
      evacflag=1;
      break;
    }
  }
#ifdef pp_SHOOTER
  shooter_flag=0;
  if(visShooter!=0&&shooter_active==1){
    shooter_flag=1;
  }
#endif

  if( plotstate==DYNAMIC_PLOTS && 
    ( sliceflag==1 || vsliceflag==1 || partflag==1 || patchflag==1 ||
#ifdef pp_SHOOTER
    shooter_flag==1||
#endif
    smoke3dflag==1|| showtour_dialog==1 || evacflag==1||
    (ReadZoneFile==1&&visZone==1&&visTimeZone==1)||
    (ReadTargFile==1&&visTarg==1)
    ||showterrain==1||showvolrender==1
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
      mesh *meshi;

      meshi=meshinfo+i;
      meshi->visInteriorPatches=0;
    }
    if(showpatch==1&&visPatchType[0]==1){
      for(i=0;i<nmeshes;i++){
        mesh *meshi;

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
    if(ReadTargFile==1&&visTarg==1)showtarget=1;
#ifdef pp_SHOOTER
    if(shooter_flag==1)showshooter=1;
#endif    
  }
  if(showsmoke==1||showevac==1||showpatch==1||showslice==1||showvslice==1||showzone==1||showiso==1||showevac==1)RenderTime=1;
  if(showtour_dialog==1||show3dsmoke==1||touring==1||showvolrender==1)RenderTime=1;
#ifdef pp_SHOOTER
  if(showshooter==1)RenderTime=1;
#endif
  if(plotstate==STATIC_PLOTS&&ReadPlot3dFile==1&&plotn>0&&plotn<=numplot3dvars)showplot3d=1;
  if(showplot3d==1){
    for(i=0;i<nmeshes;i++){
      mesh *meshi;

      meshi=meshinfo+i;
      ii=meshi->plot3dfilenum;
      if(ii==-1)continue;
      if(plot3dinfo[ii].loaded==0)continue;
      if(plot3dinfo[ii].display==0)continue;
      if(plot3dinfo[ii].extreme_min[plotn-1]==1)show_extreme_below=1;
    }
    for(i=0;i<nmeshes;i++){
      mesh *meshi;

      meshi=meshinfo+i;
      ii=meshi->plot3dfilenum;
      if(ii==-1)continue;
      if(plot3dinfo[ii].loaded==0)continue;
      if(plot3dinfo[ii].display==0)continue;
      if(plot3dinfo[ii].extreme_max[plotn-1]==1)show_extreme_above=1;
    }
  }

  numColorbars=0;
  if(ReadEvacFile==1)numColorbars++;
  if(ReadPartFile==1)numColorbars++;
  if(plotstate==DYNAMIC_PLOTS&&(slicecolorbarflag==1||vslicecolorbarflag==1))numColorbars++;
  if(plotstate==DYNAMIC_PLOTS&&patchflag==1)numColorbars++;
  if(plotstate==DYNAMIC_PLOTS&&ReadZoneFile==1)numColorbars++;
  if(plotstate==DYNAMIC_PLOTS&&tisoflag==1){
    showiso_colorbar=1;
    numColorbars++;
  }
  if(ReadPlot3dFile==1&&numColorbars==0)numColorbars=1;
  /* note: animated iso-contours do not need a color bar,
           so we don't test for isosurface files */
  dwinWW = numColorbars*dwinW/3;
  if(fontindex==1)dwinWW=(int)(1.5*dwinWW);
  drawColorLabel=0;
  if((showtime==1||showplot3d==1)&&visColorLabels==1)drawColorLabel=1;
  if(drawColorLabel==1&&olddrawColorLabel==0)updatemenu=1;
  if(drawColorLabel==0&&olddrawColorLabel==1)updatemenu=1;
  olddrawColorLabel=drawColorLabel;
  if(showtime2==1)showtime=1;
  if(plotstate==DYNAMIC_PLOTS&&stept==1){
    glutIdleFunc(Idle);
  }
  else{
    glutIdleFunc(NULL);
  }

}

/* ------------------ synctimes ------------------------ */

void synctimes(void){
  int n,i,istart;
  int j,igrid,jj;


  /* synchronize smooth blockage times */

  if(ntotal_smooth_blockages>0){
    for(igrid=0;igrid<nmeshes;igrid++){
      mesh *meshi;

      meshi=meshinfo+igrid;
      if(meshi->showsmoothtimelist==NULL)continue;
      for(n=0;n<nglobal_times;n++){
        smoothblockage *sb;

        sb = getsmoothblockage(meshi,global_times[n]);
        meshi->showsmoothtimelist[n] = sb;
      }
    }
  }

  for(n=0;n<nglobal_times;n++){

  /* synchronize tour times */

    {
      tourdata *tourj; 
      for(j=0;j<ntours;j++){
        tourj = tourinfo + j;
        if(tourj->display==0)continue;
        if(n==0){
          istart=0;
        }
        else{
          istart=tourj->timeslist[n-1];
        }
        i=istart;
        while(tourj->path_times[i]<global_times[n]&&i<tourj->ntimes){
          i++;
        }
        if(i>=tourj->ntimes){
          i--;
        }
        tourj->timeslist[n]=i;
      }
    }

    /* synchronize terrain times */

    for(j=0;j<nterraininfo;j++){
      terraindata *terri;

      terri = terraininfo + j;
      if(terri->loaded==0)continue;
      if(n==0){
        istart=0;
      }
      else{
        istart=terri->timeslist[n-1];
      }
      i=istart;
      while(terri->times[i]<global_times[n]&&i<terri->ntimes){
        i++;
      }
      if(i>=terri->ntimes){
        i--;
      }
      terri->timeslist[n]=i;
    }
    if(hrrinfo!=NULL&&hrrinfo->loaded==1&&hrrinfo->display==1){
      if(n==0){
        istart=0;
      }
      else{
        istart=hrrinfo->timeslist[n-1];
      }
      i=istart;
      while(hrrinfo->times[i]<global_times[n]&&i<hrrinfo->ntimes){
        i++;
      }
      if(i>=hrrinfo->ntimes){
        i--;
      }
      hrrinfo->timeslist[n]=i;
    }

  /* synchronize geometry times */

    for(j=0;j<ngeominfoptrs;j++){
      geomdata *geomi;

      geomi = geominfoptrs[j];
      if(geomi->loaded==0||geomi->display==0)continue;
      if(n==0){
        istart=0;
      }
      else{
        istart=geomi->timeslist[n-1];
      }
      i=istart;
      while(geomi->times[i]<global_times[n]&&i<geomi->ntimes){
        i++;
      }
      if(i>=geomi->ntimes){
        i--;
      }
      geomi->timeslist[n]=i;
    }

  /* synchronize particle times */

    for(j=0;j<npartinfo;j++){
      particle *parti;

      parti=partinfo+j;
      if(parti->loaded==0)continue;
      if(n==0){
        istart=0;
      }
      else{
        istart=parti->timeslist[n-1];
      }
      i=istart;
      while(parti->ptimes[i]<global_times[n]&&i<parti->ntimes){
        i++;
      }
      if(i>=parti->ntimes){
        i--;
      }
      parti->timeslist[n]=i;
    }

    /* synchronize target times */

    if(ntarginfo>0){
      if(n==0){
        istart=0;
      }
      else{
        istart=targtimeslist[n-1];
      }
      i=istart;
      while(targtimes[i]<global_times[n]&&i<ntargtimes){
        i++;
      }
      if(i>=ntargtimes){
        i--;
      }
      targtimeslist[n]=i;
    }

  /* synchronize shooter times */
#ifdef pp_SHOOTER
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
#endif

  /* synchronize slice times */

    for(jj=0;jj<nslice_loaded;jj++){
      slice *sd;

      j = slice_loaded_list[jj];
      sd = sliceinfo + j;
      if(n==0){
        istart=0;
      }
      else{
        istart=sd->timeslist[n-1];
      }
      i=istart;
      while(sd->times[i]<global_times[n]&&i<sd->ntimes){
        i++;
      }
      if(i>=sd->ntimes){
        i=sd->ntimes-1;
      }
      sd->timeslist[n]=i;
    }

  /* synchronize smoke times */
    {
      smoke3d *smoke3di;

      for(jj=0;jj<nsmoke3dinfo;jj++){
        smoke3di = smoke3dinfo + jj;
        if(smoke3di->loaded==0)continue;
        if(n==0){
         istart=0;
        }
        else{
          istart=smoke3di->timeslist[n-1];
        }
        i=istart;
        while(smoke3di->times[i]<global_times[n]&&i<smoke3di->ntimes){
          i++;
        }
        if(i>=smoke3di->ntimes){
          i=smoke3di->ntimes-1;
        }
        smoke3di->timeslist[n]=i;
      }
    }

  /* synchronize patch times */

    for(j=0;j<npatchinfo;j++){
      patch *patchi;

      patchi = patchinfo + j;
      if(patchi->loaded==0)continue;
      if(patchi->filetype!=2)continue;
      if(n==0){
        istart=0;
      }
      else{
        istart=patchi->geom_timeslist[n-1];
      }
      i=istart;
      while(patchi->geom_times[i]<global_times[n]&&i<patchi->ngeom_times){
        i++;
      }
      if(i>=patchi->ngeom_times){
        i=patchi->ngeom_times-1;
      }
      patchi->geom_timeslist[n]=i;
    }
    for(j=0;j<nmeshes;j++){
      patch *patchi;
      mesh *meshi;

      meshi=meshinfo+j;
      if(meshi->patchfilenum<0||meshi->patch_times==NULL)continue;
      patchi=patchinfo+meshi->patchfilenum;
      if(patchi->filetype==2)continue;
      if(n==0){
        istart=0;
      }
      else{
        istart=meshi->patch_timeslist[n-1];
      }
      i=istart;
      while(meshi->patch_times[i]<global_times[n]&&i<meshi->npatch_times){
        i++;
      }
      if(i>=meshi->npatch_times){
        i=meshi->npatch_times-1;
      }
      meshi->patch_timeslist[n]=i;
    }

  /* synchronize isosurface times */

    for(igrid=0;igrid<nmeshes;igrid++){
      mesh *meshi;

      meshi=meshinfo+igrid;
      if(meshi->iso_times==NULL)continue;
      if(n==0){
        istart=0;
      }
      else{
        istart=meshi->iso_timeslist[n-1];
      }
      i=istart;
      while(meshi->iso_times[i]<global_times[n]&&i<meshi->niso_times){
        i++;
      }
      if(i>=meshi->niso_times){
        i=meshi->niso_times-1;
      }
      meshi->iso_timeslist[n]=i;
    }

  /* synchronize volume render times */

    if(nvolrenderinfo>0){
      for(igrid=0;igrid<nmeshes;igrid++){
        volrenderdata *vr;
        mesh *meshi;

        meshi=meshinfo+igrid;
        vr = &meshi->volrenderinfo;
        if(vr->smoke==NULL)continue;
        if(vr->loaded==0||vr->display==0)continue;
        if(vr->times==NULL)continue;
        if(n==0){
          istart=0;
        }
        else{
          istart=vr->timeslist[n-1];
        }
        i=istart;
        while(vr->times[i]<global_times[n]&&i<vr->ntimes){
          i++;
        }
        if(i>=vr->ntimes){
          i=vr->ntimes-1;
        }
        vr->timeslist[n]=i;
      }
    }
    /* synchronize zone times */

    if(showzone==1){
      if(n==0){
        istart=0;
      }
      else{
        istart=zonetlist[n-1];
      }
      i=istart;
      while(zonet[i]<global_times[n]&&i<nzonet){
        i++;
      }
      if(i>=nzonet)i=nzonet-1;
      zonetlist[n]=i;
    }

  }
  reset_gltime();
}

/* ------------------ updatetimes ------------------------ */

void updatetimes(void){
  int n,n2,ntimes2;
  float *timescopy;
  int i,k;
  slice *sd;
  iso *ib;
  blockagedata *bc;
  ventdata *vi;
  particle *parti;
  tourdata *touri;
  int filenum;
  float dt_MIN=100000.0;

  FREEMEMORY(geominfoptrs);
  ngeominfoptrs=0;
  GetGeomInfoPtrs(&geominfoptrs,&ngeominfoptrs);

  // pass 1 - determine ntimes

  updateShow();  
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
#ifdef pp_SHOOTER
  if(visShooter!=0&&shooter_active==1){
    nglobal_times+=nshooter_frames;
  }
#endif
  for(i=0;i<ntours;i++){
    touri = tourinfo + i;
    if(touri->display==0)continue;
    nglobal_times+= touri->ntimes;
  }
  for(i=0;i<npartinfo;i++){
    parti = partinfo + i;
    if(parti->loaded==0)continue;
    nglobal_times+= parti->ntimes;
  }
  for(i=0;i<nsliceinfo;i++){
    sd=sliceinfo+i;
    if(sd->loaded==1||sd->vloaded==1){
      nglobal_times+=sd->ntimes;
    }
  }
  if(ReadTargFile==1&&visTarg==1){
    nglobal_times+=ntargtimes;
  }
  for(i=0;i<npatchinfo;i++){
    patch *patchi;

    patchi = patchinfo + i;
    if(patchi->loaded==1&&patchi->filetype==2){
      nglobal_times+=patchi->ngeom_times;
    }
  }
  for(i=0;i<nmeshes;i++){
    patch *patchi;
    mesh *meshi;

    meshi=meshinfo+i;
    filenum =meshi->patchfilenum;
    if(filenum!=-1){
      patchi=patchinfo+filenum;
      if(patchi->loaded==1&&patchi->filetype!=2){
        nglobal_times+=meshi->npatch_times;
      }
    }
  }
  if(ReadZoneFile==1&&visZone==1){
    nglobal_times+=nzonet;
  }
  if(ReadIsoFile==1&&visAIso!=0){
    for(i=0;i<nmeshes;i++){
      mesh *meshi;

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
      mesh *meshi;

      meshi=meshinfo+i;
      vr = &meshi->volrenderinfo;
      if(vr->fire==NULL||vr->smoke==NULL)continue;
      if(vr->loaded==0||vr->display==0)continue;
      nglobal_times+=vr->ntimes;
    }
  }
  {
    smoke3d *smoke3di;

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
#ifdef pp_SHOOTER
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
#endif

  for(i=0;i<ntours;i++){
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

  for(i=0;i<npartinfo;i++){
    parti = partinfo + i;
    if(parti->loaded==0)continue;
    for(n=0;n<parti->ntimes;n++){
      float t_diff;

      *timescopy++=parti->ptimes[n];
      t_diff = timescopy[-1]-timescopy[-2];
      if(n>0&&t_diff<dt_MIN&&t_diff>0.0){
        dt_MIN=t_diff;
      }
    }
  }

  for(i=0;i<nsliceinfo;i++){
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

  if(ReadTargFile==1&&visTarg==1){
    for(n=0;n<ntargtimes;n++){
      float t_diff;

      *timescopy++=targtimes[n];
      t_diff = timescopy[-1]-timescopy[-2];
      if(n>0&&t_diff<dt_MIN&&t_diff>0.0){
        dt_MIN=t_diff;
      }
    }
  }
  for(i=0;i<npatchinfo;i++){
    patch *patchi;

    patchi = patchinfo + i;
    if(patchi->loaded==1&&patchi->filetype==2){
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
    patch *patchi;
    mesh *meshi;

    meshi=meshinfo + i;
    filenum=meshi->patchfilenum;
    if(filenum!=-1){
      patchi = patchinfo + filenum;
      if(patchi->loaded==1&&patchi->filetype!=2){
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
      mesh *meshi;


      meshi=meshinfo + i;
      vr = &meshi->volrenderinfo;
      if(vr->smoke==NULL)continue;
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
    for(n=0;n<nzonet;n++){
      float t_diff;

      *timescopy++=zonet[n];
      t_diff = timescopy[-1]-timescopy[-2];
      if(n>0&&t_diff<dt_MIN&&t_diff>0.0){
        dt_MIN=t_diff;
      }
    }
  }
  if(ReadIsoFile==1&&visAIso!=0){
    for(i=0;i<nisoinfo;i++){
      mesh *meshi;

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
    smoke3d *smoke3di;

    if(Read3DSmoke3DFile==1&&vis3DSmoke3D==1){
      for(i=0;i<nsmoke3dinfo;i++){
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

  if(nglobal_times>0)qsort( (float *)global_times, (size_t)nglobal_times, sizeof( float ), compare_float );
  n2=1;
  ntimes2=nglobal_times;
  for(n=1;n<nglobal_times;n++){
    if(ABS(global_times[n]-global_times[n-1])>dt_MIN/10.0){
      global_times[n2]=global_times[n];
      n2++;
    }
    else{
      ntimes2--;
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
    parti=partinfo+i;
    FREEMEMORY(parti->timeslist);
    if(nglobal_times>0)NewMemory((void **)&parti->timeslist,nglobal_times*sizeof(int));
  }
  for(i=0;i<ntours;i++){
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
#ifdef pp_SHOOTER
  FREEMEMORY(shooter_timeslist);
  if(visShooter!=0&&shooter_active==1){
    NewMemory((void **)&shooter_timeslist,nshooter_frames*sizeof(int));
  }
#endif

  for(i=0;i<nsliceinfo;i++){
    sd = sliceinfo + i;
    FREEMEMORY(sd->timeslist);
    if(nglobal_times>0)NewMemory((void **)&sd->timeslist,nglobal_times*sizeof(int));
  }
  if(nvolrenderinfo>0){
    for(i=0;i<nmeshes;i++){
      mesh *meshi;
      volrenderdata *vr;

      meshi = meshinfo + i;
      vr = &(meshi->volrenderinfo);
      if(vr->fire==NULL||vr->smoke==NULL)continue;
      if(vr->loaded==0||vr->display==0)continue;
      FREEMEMORY(vr->timeslist);
      if(nglobal_times>0)NewMemory((void **)&vr->timeslist,nglobal_times*sizeof(int));
    }
  }
  {
    smoke3d *smoke3di;

    for(i=0;i<nsmoke3dinfo;i++){
      smoke3di = smoke3dinfo + i;
      FREEMEMORY(smoke3di->timeslist);
      if(nglobal_times>0)NewMemory((void **)&smoke3di->timeslist,nglobal_times*sizeof(int));
    }
  }
  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi=meshinfo+i;
    if(meshi->iso_times==NULL)continue;
    FREEMEMORY(meshi->iso_timeslist);
    if(nglobal_times>0)NewMemory((void **)&meshi->iso_timeslist,  nglobal_times*sizeof(int));  
  }

  for(i=0;i<npatchinfo;i++){
    patch *patchi;

    patchi = patchinfo + i;
    FREEMEMORY(patchi->geom_timeslist);
    if(patchi->filetype!=2)continue;
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

  FREEMEMORY(zonetlist); 
  if(nglobal_times>0)NewMemory((void **)&zonetlist,     nglobal_times*sizeof(int));

  FREEMEMORY(targtimeslist);
  if(nglobal_times>0)NewMemory((void **)&targtimeslist,  nglobal_times*sizeof(int));

  if(ntotal_smooth_blockages>0){
    for(i=0;i<nmeshes;i++){
      mesh *meshi;

      meshi=meshinfo+i;
      FREEMEMORY(meshi->showsmoothtimelist);
      if(nglobal_times>0)NewMemory((void **)&meshi->showsmoothtimelist,nglobal_times*sizeof(smoothblockage *));
    }
  }

  // end pass 3

  // reset render_fram array

  if(current_script_command!=NULL&&current_script_command->command==SCRIPT_VOLSMOKERENDERALL){
    if(current_script_command->first==1){
      for(n=0;n<nglobal_times;n++){
        render_frame[n]=0;
      }
      current_script_command->first=0;
    }
  }
  else{
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
  reset_itimes0();
  for(i=0;i<ngeominfoptrs;i++){
    geomdata *geomi;

    geomi = geominfoptrs[i];
    if(geomi->loaded==0||geomi->display==0)continue;
    geomi->iframe=0;
  }
  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi=meshinfo+i;
    meshi->ipatch=0;
  }
  for(i=0;i<nsliceinfo;i++){
    sd = sliceinfo + i;
    sd->islice=0; 
  }
  iframe=iframebeg; 
  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi=meshinfo+i;
    if(meshi->iso_times==NULL)continue;
    meshi->iiso=0;
  }
  for(i=0;i<npartinfo;i++){
    parti = partinfo + i;
    parti->iframe=0;
  }

  /* determine visibility of each blockage at each time step */

  for(i=0;i<nmeshes;i++){
    int j;
    mesh *meshi;

    meshi=meshinfo+i;
    for(j=0;j<meshi->nbptrs;j++){
      bc = meshi->blockageinfoptrs[j];
      if(bc->showtime==NULL)continue;
      FREEMEMORY(bc->showtimelist);
      if(nglobal_times>0){
        NewMemory((void **)&bc->showtimelist,nglobal_times*sizeof(int));
        for(k=0;k<nglobal_times;k++){
          int listindex;

          bc->showtimelist[k]=1;
          listindex=getindex(global_times[k],bc->showtime,bc->nshowtime);
          bc->showtimelist[k]=bc->showhide[listindex];
        }
      }
    }
  }

  /* determine state of each device at each time step */

  for(i=0;i<ndeviceinfo;i++){
    device *devicei;

    devicei = deviceinfo + i;
    if(devicei->object->visible==0)continue;
    if(devicei->nstate_changes==0)continue;
    FREEMEMORY(devicei->showstatelist);
    if(nglobal_times>0){
      NewMemory((void **)&devicei->showstatelist,nglobal_times*sizeof(int));
      for(k=0;k<nglobal_times;k++){
        int listindex;

        listindex=getindex(global_times[k],devicei->act_times,devicei->nstate_changes);
        devicei->showstatelist[k]=devicei->state_values[listindex];
      }
    }
  }

  /* determine visibility of each vent at each time step */

  for(i=0;i<nmeshes;i++){
    int j;
    mesh *meshi;

    meshi=meshinfo+i;
    if(meshi->ventinfo==NULL)continue;
    for(j=0;j<meshi->nvents;j++){
      vi = meshi->ventinfo+j;
      if(vi->showtime==NULL)continue;
      FREEMEMORY(vi->showtimelist);
      if(nglobal_times>0){
        NewMemory((void **)&vi->showtimelist,nglobal_times*sizeof(int));
        for(k=0;k<nglobal_times;k++){
          int listindex;

          vi->showtimelist[k]=1;
          listindex=getindex(global_times[k],vi->showtime,vi->nshowtime);
          vi->showtimelist[k]=vi->showhide[listindex];
        }
      }
    }
  }

  if(nglobal_times>0)synctimes();
  updatefaces=1;
  if(nglobal_times>0){
    UpdateTimeLabels();
    updateGluiTimeBounds(global_times[0],global_times[nglobal_times-1]);
  }
  show_slice_terrain=0;
  if(visTerrainType==TERRAIN_3D_MAP){
    for(i=0;i<nsliceinfo;i++){
      sd = sliceinfo + i;
      if(sd->loaded==0||sd->display==0||sd->terrain==0)continue;
      show_slice_terrain=1;
      break;
    }
  }
}

/* ------------------ getplotstate ------------------------ */

int getplotstate(int choice){
  int i;
  plot3d *ploti;
  slice *slicei;
  vslice *vslicei;
  patch *patchi;
  particle *parti;
  targ *targi;
  iso *isoi;
  zonedata *zonei;
  smoke3d *smoke3di;
  tourdata *touri;
  int ii;

  update_loaded_lists();
  switch (choice){
    case STATIC_PLOTS:
    case STATIC_PLOTS_NORECURSE:
      stept = 0;
      for(i=0;i<nmeshes;i++){
        mesh *meshi;

        meshi=meshinfo + i;
        if(meshi->plot3dfilenum==-1)continue;
        ploti = plot3dinfo + meshi->plot3dfilenum;
        if(ploti->loaded==0||ploti->display==0)continue;
        if(meshi->visx==0&&meshi->visy==0&&meshi->visz==0&&visiso==0)continue;
        return STATIC_PLOTS;
      }
      if(choice!=STATIC_PLOTS_NORECURSE){
        return getplotstate(DYNAMIC_PLOTS_NORECURSE);
      }
      break;
    case DYNAMIC_PLOTS:
    case DYNAMIC_PLOTS_NORECURSE:
      for(ii=0;ii<nslice_loaded;ii++){
        i = slice_loaded_list[ii];
        slicei = sliceinfo + i;
        if(slicei->display==0||slicei->type!=islicetype)continue;
        if(slicei->volslice==0&&visGrid==0)stept = 1; 
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
      for(i=0;i<nvslice;i++){
        vslicei = vsliceinfo + i;
        if(vslicei->display==0||vslicei->type!=islicetype)continue;
        return DYNAMIC_PLOTS;
      }
      for(ii=0;ii<npatch_loaded;ii++){
        i = patch_loaded_list[ii];
        patchi = patchinfo + i;
        if(patchi->display==0||patchi->type!=ipatchtype)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<npartinfo;i++){
        parti = partinfo + i;
        if(parti->loaded==0||parti->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<nisoinfo;i++){
        isoi = isoinfo + i;
        if(isoi->loaded==0)continue;
        if(isoi->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<nzone;i++){
        zonei = zoneinfo + i;
        if(zonei->loaded==0||zonei->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<ntarginfo;i++){
        targi = targinfo + i;
        if(targi->loaded==0||targi->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<ntours;i++){
        touri = tourinfo + i;
        if(touri->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      for(i=0;i<nsmoke3dinfo;i++){
        smoke3di = smoke3dinfo + i;
        if(smoke3di->loaded==0||smoke3di->display==0)continue;
        return DYNAMIC_PLOTS;
      }
      if(nvolrenderinfo>0){
        for(i=0;i<nmeshes;i++){
          mesh *meshi;
          volrenderdata *vr;

          meshi = meshinfo + i;
          vr = &(meshi->volrenderinfo);
          if(vr->fire==NULL||vr->smoke==NULL)continue;
          if(vr->loaded==0||vr->display==0)continue;
          return DYNAMIC_PLOTS;
        }
      }
#ifdef pp_SHOOTER
      if(visShooter!=0&&shooter_active==1){
        return DYNAMIC_PLOTS;
      }
#endif
      if(choice!=DYNAMIC_PLOTS_NORECURSE)return getplotstate(STATIC_PLOTS_NORECURSE);
      break;
    default:
      ASSERT(FFALSE);
      break;
  }
  stept = 0;
  return NO_PLOTS;
}

/* ------------------ getindex ------------------------ */

int getindex(float key, const float *list, int nlist){
  int i;
  if(nlist==1)return 0;
  if(key<list[1])return 0;
  if(key>=list[nlist-1])return nlist-1;
  for(i=1;i<nlist-1;i++){
    if(list[i]<=key&&key<list[i+1])return i;
  }
  return 0;
}

/* ------------------ isearch ------------------------ */

int isearch(float *list, int nlist, float key, int guess){
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

/* ------------------ isearch ------------------------ */

void reset_itimes0(void){
  if(current_script_command==NULL||current_script_command->command!=SCRIPT_VOLSMOKERENDERALL){
    itimes=0;
  }
}


/* ------------------ updateclipbounds ------------------------ */

void updateclipbounds(int set_i0, int *i0, int set_i1, int *i1, int imax){ 

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

/* ------------------ updateclip ------------------------ */

void updateclip(int slicedir){
  stepclip_x=0; stepclip_y=0; stepclip_z=0; 
  stepclip_X=0; stepclip_Y=0; stepclip_Z=0;
  switch (slicedir){
  case 1:
    clip_x = 1 - clip_x;
    if(clip_x==1)printf("clip x on\n");
    if(clip_x==0)printf("clip x off\n");
    if(clip_x==1)stepclip_x=1;
    break;
  case 2:
    clip_y = 1 - clip_y;
    if(clip_y==1)printf("clip y on\n");
    if(clip_y==0)printf("clip y off\n");
    if(clip_y==1)stepclip_y=1;
    break;
  case 3:
    clip_z = 1 - clip_z;
    if(clip_z==1)printf("clip z on\n");
    if(clip_z==0)printf("clip z off\n");
    if(clip_z==1)stepclip_z=1;
    break;
  case -1:
    clip_X = 1 - clip_X;
    if(clip_X==1)printf("clip X on\n");
    if(clip_X==0)printf("clip X off\n");
    if(clip_X==1)stepclip_X=1;
    break;
  case -2:
    clip_Y = 1 - clip_Y;
    if(clip_Y==1)printf("clip Y on\n");
    if(clip_Y==0)printf("clip Y off\n");
    if(clip_Y==1)stepclip_Y=1;
    break;
  case -3:
    clip_Z = 1 - clip_Z;
    if(clip_Z==1)printf("clip Z on\n");
    if(clip_Z==0)printf("clip Z off\n");
    if(clip_Z==1)stepclip_Z=1;
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
}

/* ------------------ update_glui_names ------------------------ */

void update_glui_names(void){
  update_glui_dialogs=1;

  hide_glui_colorbar();
  update_glui_colorbar=1;

  hide_glui_bounds();
  hide_glui_3dsmoke();
  update_glui_bounds=1;

  hide_glui_edit();
  update_glui_edit=1;
  
  hide_glui_clip();
  update_glui_clip=1;

  hide_glui_device();
  update_glui_device=1;

  hide_glui_display();
  update_glui_labels=1;
  
  hide_glui_message();
  update_glui_message=1;
  
  hide_glui_motion();
  update_glui_motion=1;

#ifdef pp_SHOOTER
  hide_glui_shooter();
  update_glui_shooter=1;
#endif
  
  hide_glui_tour();
  update_glui_tour=1;
  
  hide_glui_trainer();
  update_glui_trainer=1;

  hide_glui_alert();
  update_glui_alert=1;
  
  hide_glui_stereo();
  update_glui_stereo=1;

  hide_glui_wui();
  update_glui_wui=1;
}

void glui_advancedtour_setup(int main_window);
void glui_tour_setup(int main_window);
void glui_colorbar_setup(int main_window);
void glui_bounds_setup(int main_window);
void glui_3dsmoke_setup(int main_window);
void glui_edit_setup(int main_window);
void glui_clip_setup(int main_window);
void glui_device_setup(int main_window);
void glui_labels_setup(int main_window);
void glui_message_setup(int main_window);
void glui_motion_setup(int main_window);
void glui_shooter_setup(int main_window);
void glui_trainer_setup(int main_window);
void glui_stereo_setup(int main_window);
void glui_wui_setup(int main_window);

/* ------------------ update_glui_dialogs ------------------------ */

void Update_Glui_Dialogs(void){
  if(update_glui_colorbar==1)glui_colorbar_setup(mainwindow_id);
  if(update_glui_bounds==1){
    glui_bounds_setup(mainwindow_id);
    glui_3dsmoke_setup(mainwindow_id);
  }
  if(update_glui_edit==1)glui_edit_setup(mainwindow_id);
  if(update_glui_clip==1)glui_clip_setup(mainwindow_id);
  if(update_glui_device==1)glui_device_setup(mainwindow_id);
  if(update_glui_labels==1)glui_labels_setup(mainwindow_id);
  if(update_glui_message==1)glui_message_setup(mainwindow_id);
  if(update_glui_motion==1)glui_motion_setup(mainwindow_id);
#ifdef pp_SHOOTER
  if(update_glui_shooter==1)glui_shooter_setup(mainwindow_id);
#endif
  if(update_glui_tour==1){
    glui_tour_setup(mainwindow_id);
    glui_advancedtour_setup(mainwindow_id);
  }
  if(update_glui_trainer==1)glui_trainer_setup(mainwindow_id);
  if(update_glui_alert==1)glui_alert_setup(mainwindow_id);
  if(update_glui_stereo==1)glui_stereo_setup(mainwindow_id);
  if(update_glui_wui==1)glui_wui_setup(mainwindow_id);
}

/* ------------------ Show_Glui_Dialogs ------------------------ */

void Show_Glui_Dialogs(void){
  if(showcolorbar_dialog_save==1)show_glui_colorbar();

  if(showbounds_dialog_save==1)show_glui_bounds();
  if(showmotion_dialog_save==1)show_glui_motion();
  if(showedit_dialog_save==1)show_glui_edit();
  if(showclip_dialog_save==1)show_glui_clip();
  if(showstereo_dialog_save==1)show_glui_stereo();
  if(showtour_dialog_save==1)show_glui_tour();
  if(showdisplay_dialog_save==1)show_glui_display();
  if(showwui_dialog_save==1)show_glui_wui();
  if(showdevice_dialog_save==1)show_glui_device();
  if(show3dsmoke_dialog_save==1)show_glui_3dsmoke();
  if(showtrainer_dialog_save==1)show_glui_trainer();
}






