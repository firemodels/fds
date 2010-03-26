// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "flowfiles.h"
#include "MALLOC.h"
#include <pthread.h>
#ifdef pp_SPHERE
#include "csphere.h"
#endif

#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

// svn revision character string
char readsmv_revision[]="$Revision$";

#define DEVICE_DEVICE 0
#define DEVICE_THCP 1
#define DEVICE_HEAT 2     
#define DEVICE_SPRK 3
#define DEVICE_SMOKE 4

#define ijcell2(i,j) nxcell*(j) + (i)
int GeometryMenu(int var);
propdata *get_prop_id(char *prop_id);


/* ------------------ update_inilist ------------------------ */

void update_inilist(void){
  inifiledata *inifile;
  int first_time=1;
  FILE *filein_inilist=NULL;
  char ini_listfile[1024];
  char buffer[255],buffer2[255];

  strcpy(ini_listfile,caseinifilename);
  strcat(ini_listfile,"list");
  filein_inilist=fopen(ini_listfile,"r");
  if(filein_inilist!=NULL){

    while(!feof(filein_inilist)){
      CheckMemory;
      if(fgets(buffer,255,filein_inilist)==NULL)break;
      if(match(buffer,"INIFILE",7)==1){
        if(fgets(buffer2,255,filein_inilist)==NULL)break;
        cleanbuffer(buffer,buffer2);
        insert_inifile(buffer);
        continue;
        }
    }
    fclose(filein_inilist);
  }
}

/* ------------------ propi ------------------------ */

void init_prop(propdata *propi, int nsmokeview_ids, char *label){
  int nlabel;

  nlabel = strlen(label);
  if(nlabel==0){
    NewMemory((void **)&propi->label,5);
    strcpy(propi->label,"null");
  }
  else{
    NewMemory((void **)&propi->label,nlabel+1);
    strcpy(propi->label,label);
  }
  
  NewMemory((void **)&propi->smokeview_ids,nsmokeview_ids*sizeof(char *));
  NewMemory((void **)&propi->smv_objects,nsmokeview_ids*sizeof(sv_object *));

  propi->nsmokeview_ids=nsmokeview_ids;
  propi->blockvis=1;
  propi->inblockage=0;
  propi->ntextures=0;
  propi->nvars_dep=0;
  propi->nvars_evac=0;
  propi->nvars_indep=0;
  propi->vars_indep=NULL;
  propi->svals=NULL;
  propi->texturefiles=NULL;
  propi->rotate_axis=NULL;
  propi->rotate_angle=0.0;
}

/* ------------------ readsmv ------------------------ */

int readsmv(char *file){

/* read the .smv file */

  device *devicecopy;
  int do_pass4=0;
  outline *outlinei;
  int roomdefined=0;
  int haveSHELL=0;
  float *x1, *x2, *yy1, *yy2, *z1, *z2;
  float temp_ignition, emis, t_width, t_height;
  float s_color[4];
  int *ijk;
  int s_type;
  int errorcode;
  int noGRIDpresent=1,startpass;

  int ipart=0, islice=0, ipatch=0, iplot3d=0, iroom=0,izone=0,ifire=0,iiso=0;
  int ismoke3d=0;
  int  itarg=0;
  int factor=256*128;
  int setGRID=0;
  int idummy;
  int nbtemp,nvents;
  int tempval;
  int iv1, iv2;
  int jv1, jv2;
  int kv1, kv2;
  float *xp=NULL, *yp=NULL, *zp=NULL;
  float *xp2=NULL, *yp2=NULL, *zp2=NULL;
  float *xpltcopy, *ypltcopy, *zpltcopy;
  float dxbar, dybar, dzbar;
  float dxsbar, dysbar, dzsbar;
  int showobst;
  float time;
  int blocktemp;
  ventdata *vinfo,*vi;
  int showvent;
  int colorindex, blocktype;
  int ventindex,venttype;
  int dataflag;
  int roomnumber;
  float width,ventoffset,bottom,top;
  blockagedata *bc;
  plot3d *p;
  int igrid;
  int ioffset;
  float *xplttemp,*yplttemp,*zplttemp;
  float *xplt_origtemp,*yplt_origtemp,*zplt_origtemp;
  particle *parti;
  int itrnx, itrny, itrnz, ipdim, iobst, ivent;
  int ibartemp=2, jbartemp=2, kbartemp=2;
  size_t len;
  int  n;
  FILE *stream;
  FILE *ENDIANfile;
  char buffer[255],buffer2[255],*bufptr;
  char *buffer3;
  int blocknumber;
  float *xsprcopy, *ysprcopy, *zsprcopy;
  float *xheatcopy, *yheatcopy, *zheatcopy;
  int nn;
  char *UVEL="U-VEL";
  char *VVEL="V-VEL";
  char *WVEL="W-VEL";
  int i, j, k;
  STRUCTSTAT statbuffer,statbuffer2;
  texture *texti,*textj;
  cadgeom *cd;
  surface *surfi;
  int dup_texture;
  int version;

  int nn_smoke3d=0;
  int nn_patch=0;
  int nn_iso=0;
  int nn_part=0;
  int nn_plot3d=0;
  int nn_slice=0;

  npropinfo=0;
  navatar_colors=0;
  FREEMEMORY(avatar_colors);

  FREEMEMORY(treeinfo);
  ntreeinfo=0;
  for(i=0;i<nterraininfo;i++){
    terraindata *terri;

    terri = terraininfo + i;
    FREEMEMORY(terri->x);
    FREEMEMORY(terri->y);
    FREEMEMORY(terri->zcell);
    FREEMEMORY(terri->znode);
//    FREEMEMORY(terri->znormal);
  }
  FREEMEMORY(terraininfo);
  nterraininfo=0;
  niso_compressed=0;
#ifdef pp_SPHERE
  if(sphereinfo==NULL){
    NewMemory((void **)&sphereinfo,sizeof(spherepoints));
    initspherepoints(sphereinfo,14);
  }
  if(wui_sphereinfo==NULL){
    NewMemory((void **)&wui_sphereinfo,sizeof(spherepoints));
    initspherepoints(wui_sphereinfo,14);
  }
#endif

  ntotal_blockages=0;
  ntotal_smooth_blockages=0;

  FREEMEMORY(tickinfo);
  nticks=0;
  ntickssmv=0;

  FREEMEMORY(labelinfo);
  nlabels=0;
  nlabelssmv=0;

  FREEMEMORY(camera_external);
  if(file!=NULL)NewMemory((void **)&camera_external,sizeof(camera));

  FREEMEMORY(camera_external_save);
  if(file!=NULL)NewMemory((void **)&camera_external_save,sizeof(camera));

  FREEMEMORY(camera_ini);
  if(file!=NULL){
    NewMemory((void **)&camera_ini,sizeof(camera));
    camera_ini->defined=0;
  }

  FREEMEMORY(camera_current);
  if(file!=NULL)NewMemory((void **)&camera_current,sizeof(camera));

  FREEMEMORY(camera_internal);
  if(file!=NULL)NewMemory((void **)&camera_internal,sizeof(camera));

  FREEMEMORY(camera_save);
  if(file!=NULL)NewMemory((void **)&camera_save,sizeof(camera));

  FREEMEMORY(camera_last);
  if(file!=NULL)NewMemory((void **)&camera_last,sizeof(camera));
 
  updatefaces=1;
  STRCPY(TITLE1,"");
  STRCPY(TITLE2,"");
  nfires=0;
  nrooms=0;

  if(file!=NULL){
    initsurface(&sdefault);
    NewMemory((void **)&sdefault.surfacelabel,(5+1));
    strcpy(sdefault.surfacelabel,"INERT");

    initventsurface(&v_surfacedefault);
    NewMemory((void **)&v_surfacedefault.surfacelabel,(4+1));
    strcpy(v_surfacedefault.surfacelabel,"VENT");
  
    initsurface(&e_surfacedefault);
    NewMemory((void **)&e_surfacedefault.surfacelabel,(8+1));
    strcpy(e_surfacedefault.surfacelabel,"EXTERIOR");
    e_surfacedefault.color=mat_ambient2;
  }

  // free memory for particle class

  if(partclassinfo!=NULL){
    for(i=0;i<npartclassinfo+1;i++){
      part5class *partclassi;

      partclassi = partclassinfo + i;
      FREEMEMORY(partclassi->name);
      if(partclassi->ntypes>0){
        for(j=0;j<partclassi->ntypes;j++){
          flowlabels *labelj;
          
          labelj = partclassi->labels+j;
          freelabels(labelj);
        }
        FREEMEMORY(partclassi->labels);
        partclassi->ntypes=0;
      }
    }
    FREEMEMORY(partclassinfo);
  }
  npartclassinfo=0;

  freeall_objects();
  if(ndeviceinfo>0){
    for(i=0;i<ndeviceinfo;i++){
    }
    FREEMEMORY(deviceinfo);
    ndeviceinfo=0;
  }

  if(noutlineinfo>0){
    for(i=0;i<noutlineinfo;i++){
      outlinei = outlineinfo + i;
      FREEMEMORY(outlinei->x1);
      FREEMEMORY(outlinei->y1);
      FREEMEMORY(outlinei->z1);
      FREEMEMORY(outlinei->x2);
      FREEMEMORY(outlinei->y2);
      FREEMEMORY(outlinei->z2);
    }
    FREEMEMORY(outlineinfo);
    noutlineinfo=0;
  }

  if(nzone>0){
    for(i=0;i<nzone;i++){
      for(n=0;n<4;n++)freelabels(&zoneinfo[i].label[n]);
      FREEMEMORY(zoneinfo[i].file);
    }
    FREEMEMORY(zoneinfo);
  }
  nzone=0;

  if(nplot3d_files>0){
    for(i=0;i<nplot3d_files;i++){
      for(n=0;n<6;n++)freelabels(&plot3dinfo[i].label[n]);
      FREEMEMORY(plot3dinfo[i].file);
    }
    FREEMEMORY(plot3dinfo);
  }
  nplot3d_files=0;
  if(nsmoke3d_files>0){
    {
      smoke3d *smoke3di;

      for(i=0;i<nsmoke3d_files;i++){
        smoke3di = smoke3dinfo + i;
        freesmoke3d(smoke3di);
#ifdef pp_LIGHT
        FREEMEMORY(smoke3di->light_file);
#endif
        FREEMEMORY(smoke3di->comp_file);
        FREEMEMORY(smoke3di->reg_file);
      }
      FREEMEMORY(smoke3dinfo);
      nsmoke3d_files=0;
    }
  }

  if(npart_files>0){
    for(i=0;i<npart_files;i++){
      freelabels(&partinfo[i].label);
      FREEMEMORY(partinfo[i].partclassptr);
      FREEMEMORY(partinfo[i].reg_file);
      FREEMEMORY(partinfo[i].comp_file);
      FREEMEMORY(partinfo[i].size_file);
    }
    FREEMEMORY(partinfo);
  }
  npart_files=0;

  ntarg_files=0;

  FREEMEMORY(surfaceinfo);


  //*** free slice data

  if(nslice_files>0){
    for(i=0;i<nslice_files;i++){
      slice *sd;
      sd = sliceinfo + i;
      freelabels(&sliceinfo[i].label);
      FREEMEMORY(sd->reg_file);
      FREEMEMORY(sd->comp_file);
      FREEMEMORY(sd->rle_file);
      FREEMEMORY(sd->size_file);
    }
    FREEMEMORY(sliceorderindex);
    for(i=0;i<nmultislices;i++){
      multislice *mslicei;

      mslicei = multisliceinfo + i;
      FREEMEMORY(mslicei->islices);
    }
    FREEMEMORY(multisliceinfo);
    nmultislices=0;
    FREEMEMORY(sliceinfo);
  }
  nslice_files=0;

  //*** free multi-vector slice data

  if(nvslice>0){
    FREEMEMORY(vsliceorderindex);
    for(i=0;i<nmultivslices;i++){
      multivslice *mvslicei;

      mvslicei = multivsliceinfo + i;
      FREEMEMORY(mvslicei->ivslices);
    }
    FREEMEMORY(multivsliceinfo);
    nmultivslices=0;
  }

  if(npatch_files>0){
    for(i=0;i<npatch_files;i++){
      freelabels(&patchinfo[i].label);
      FREEMEMORY(patchinfo[i].reg_file);
      FREEMEMORY(patchinfo[i].comp_file);
      FREEMEMORY(patchinfo[i].size_file);
    }
    FREEMEMORY(patchinfo);
  }
  npatch_files=0;
  
  if(niso_files>0){
    for(i=0;i<niso_files;i++){
      freelabels(&isoinfo[i].surface_label);
      FREEMEMORY(isoinfo[i].file);
    }
    FREEMEMORY(isoinfo);
  }
  niso_files=0;

  freecadinfo();

  updateindexcolors=0;
  ntrnx=0;ntrny=0;ntrnz=0;nmeshes=0;npdim=0;nvent=0;nobst=0,noffset=0;nsurfaces=0;
  nvent_transparent=0;

  nvents=0; setPDIM=0;
  endian = getendian();
  endian_native = getendian();
  endian_data=endian_native;
  FREEMEMORY(LESsystem);
  FREEMEMORY(LESendian);

  FREEMEMORY(databasefilename);

  FREEMEMORY(targinfo);

  FREEMEMORY(vsliceinfo);
  FREEMEMORY(sliceinfo);
  FREEMEMORY(slicetypes);
  FREEMEMORY(vslicetypes);

  FREEMEMORY(plot3dinfo);
  FREEMEMORY(patchinfo);
  FREEMEMORY(patchtypes);
  FREEMEMORY(isoinfo);
  FREEMEMORY(isotypes);
  FREEMEMORY(roominfo);
  FREEMEMORY(fireinfo);
  FREEMEMORY(zoneinfo);
  FREEMEMORY(zventinfo);

  FREEMEMORY(textureinfo);
  FREEMEMORY(surfaceinfo);
  FREEMEMORY(terrain_texture);

  if(cadgeominfo!=NULL)freecadinfo();

  if(file==NULL){
    initvars1();
    return -1;  // finished  unloading memory from previous case
  }

  if(NewMemory((void **)&LESsystem,4)==0)return 2;
  STRCPY(LESsystem,"");
  if(NewMemory((void **)&LESendian,4)==0)return 2;
  STRCPY(LESendian,"");

  if( (stream=fopen(file,"r"))==NULL)return 1;
  getfile_modtime(file, &smv_modtime);
  
  printf("\nReading: %s\n",file);

/* 
   ************************************************************************
   ************************ start of pass 1 ********************************* 
   ************************************************************************
 */

  nbtemp=0; nvents=0; igrid=0; ioffset=0;
  ntc_total=0, nspr_total=0, nheat_total=0;
  printf("reading input file\n");
  printf("   pass 1\n");
  while(!feof(stream)){
    if(fgets(buffer,255,stream)==NULL)break;
    if(strncmp(buffer," ",1)==0)continue;

    /* 
      The keywords TRNX, TRNY, TRNZ, GRID, PDIM, OBST and VENT are not required 
      BUT if any one these keywords are present then the number of each MUST be equal 
    */


    if(match(buffer,"PROP",4) == 1){
      npropinfo++;
      continue;
    }
    if(match(buffer,"SMOKEDIFF",9) == 1){
      smokediff=1;
      continue;
    }
    if(match(buffer,"AVATAR_COLOR",12) == 1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&navatar_colors);
      if(navatar_colors<0)navatar_colors=0;
      if(navatar_colors>0){
        float *acolor;

        NewMemory((void **)&avatar_colors,3*navatar_colors*sizeof(float));
        acolor=avatar_colors;
        for(i=0;i<navatar_colors;i++){
          int irgb[3];
          fgets(buffer,255,stream);
          irgb[0]=0;
          irgb[1]=0;
          irgb[2]=0;
          sscanf(buffer,"%i %i %i",irgb,irgb+1,irgb+2);
          acolor[0]=(float)irgb[0]/255.0;
          acolor[1]=(float)irgb[1]/255.0;
          acolor[2]=(float)irgb[2]/255.0;
          acolor+=3;
        }
      }
      continue;
    }

    if(match(buffer,"TERRAIN",7) == 1){
      nterraininfo++;
      continue;
    }
    if(match(buffer,"CLASS_OF_PARTICLES",18) == 1||
       match(buffer,"CLASS_OF_HUMANS",15) == 1){
      npartclassinfo++;
      continue;
    }
    if(match(buffer,"PL3D",4) == 1){
      nplot3d_files++;
      continue;
   
    }
    if(match(buffer,"AUTOTERRAIN",11) == 1){
      int len_buffer;
      texture *tt;
      char *buff2;

      NewMemory((void **)&terrain_texture,sizeof(texture));
      tt = terrain_texture;
      autoterrain=1;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visTerrainType);
    //  if(visTerrain!=1)visTerrain=0;
  
      fgets(buffer,255,stream);
      buff2 = trim_front(buffer);
      trim(buff2);
      len_buffer = strlen(buff2);

      NewMemory((void **)&tt->file,(len_buffer+1)*sizeof(char));
      strcpy(tt->file,buff2);

      continue;
    }
    if(
      (match(buffer,"DEVICE",6) == 1)&&
      (match(buffer,"DEVICE_ACT",10) != 1)
      ){
      fgets(buffer,255,stream);
      fgets(buffer,255,stream);
      ndeviceinfo++;
      continue;
    }
    if(
       matchonly(buffer,"SPRK",4) == 1||
       matchonly(buffer,"HEAT",4) == 1||
       matchonly(buffer,"SMOD",4) == 1||
       matchonly(buffer,"THCP",4) == 1
    ){
      int local_ntc;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&local_ntc);
      if(local_ntc<0)local_ntc=0;
      ndeviceinfo+=local_ntc;
      for(i=0;i<local_ntc;i++){
        fgets(buffer,255,stream);
      }
      continue;
    }
    if(match(buffer,"DATABASE",8)==1){
      if(fgets(buffer,255,stream)==NULL)break;
      buffer3=trim_front(buffer);
      trim(buffer3);
      len=strlen(buffer3);
      if(STAT(buffer3,&statbuffer)==0){
        FREEMEMORY(databasefilename);
        if(NewMemory((void **)&databasefilename,(unsigned int)(len+1))==0)break;
        strcpy(databasefilename,buffer3);
      }
      continue;
    }

    if(match(buffer,"REVISION",8)==1){
      revision_fds=-1;
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%i",&revision_fds);
      if(revision_fds<0)revision_fds=-1;
      continue;
    }
    if(match(buffer,"TOFFSET",7)==1){
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%f %f %f",texture_origin,texture_origin+1,texture_origin+2);
      continue;
    }

    if(match(buffer,"USETEXTURES",11) == 1){
      usetextures=1;
      continue;
    }

    if(match(buffer,"CADTEXTUREPATH",14) == 1||
       match(buffer,"TEXTUREDIR",10) == 1){
      if(fgets(buffer,255,stream)==NULL)break;
      trim(buffer);
      {
        size_t texturedirlen;

        texturedirlen=strlen(trim_front(buffer));
        if(texturedirlen>0){
          FREEMEMORY(texturedir);
          NewMemory( (void **)&texturedir,texturedirlen+1);
          strcpy(texturedir,trim_front(buffer));
        }
      }
      continue;
    }

    if(match(buffer,"VIEWTIMES",9) == 1){
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%f %f %i",&view_tstart,&view_tstop,&view_ntimes);
      if(view_ntimes<2)view_ntimes=2;
      ReallocTourMemory();
      continue;
    }
    if(match(buffer,"OUTLINE",7) == 1){
      noutlineinfo++;
      continue;
    }
    if(match(buffer,"TICKS",5) == 1){
      nticks++;
      ntickssmv++;
      continue;
    }
    if(match(buffer,"LABEL",5) == 1){
      nlabels++;
      nlabelssmv++;
      continue;
    }
    if(match(buffer,"TRNX",4) == 1){
      ntrnx++;
      continue;
    }
    if(match(buffer,"TRNY",4) == 1){
      ntrny++;
      continue;
    }
    if(match(buffer,"TRNZ",4) == 1){
      ntrnz++;
      continue;
    }
    if(match(buffer,"SURFACE",7) ==1&&match(buffer,"SURFACE DENSITY",15)!=1){
      nsurfaces++;
      continue;
    }
    if(match(buffer,"GRID",4) == 1){
//      int lenbuffer;

//      trim(buffer);
//      lenbuffer=strlen(buffer);
//      if(lenbuffer>4){
//        if(buffer[5]!=' ')continue;
//      }
      
      noGRIDpresent=0;
      nmeshes++;
      continue;
    }
    if(match(buffer,"OFFSET",6) == 1){
      noffset++;
      continue;
    }
    if(match(buffer,"PDIM",4) == 1){
      npdim++;
      setPDIM=1;
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f %f %f %f",&xbar0,&xbar,&ybar0,&ybar,&zbar0,&zbar);
      continue;
    }
    if(match(buffer,"OBST",4) == 1){
      nobst++;
      continue;
    }
    if(match(buffer,"CADGEOM",7) == 1){
      ncadgeom++;
      continue;
    }
    if(match(buffer,"VENT",4) == 1 && match(buffer,"VENTGEOM",8) != 1){
      nvent++;
      continue;
    }
    if(match(buffer,"PART",4) == 1||match(buffer,"EVAC",4)==1
      ||match(buffer,"PRT5",4)==1||match(buffer,"EVA5",4)==1
      ){
      npart_files++;
      continue;
    }
    if( (match(buffer,"SLCF",4) == 1)||
        (match(buffer,"SLCC",4) == 1)||
        (match(buffer,"SLFL",4) == 1)||
        (match(buffer,"SLCT",4) == 1)
      ){
      nslice_files++;
      continue;
    }
    if(match(buffer,"SMOKE3D",7) == 1){
      nsmoke3d_files++;
      continue;
    }
    if(
      match(buffer,"MINMAXBNDF",10) == 1||
      match(buffer,"MINMAXPL3D",10) == 1||
      match(buffer,"MINMAXSLCF",10) == 1
      ){
      do_pass4=1;
      continue;
    }
    if(match(buffer,"BNDF",4) == 1|| match(buffer,"BNDC",4) == 1){
      npatch_files++;
      continue;
    }
    if(match(buffer,"ISOF",4) == 1||match(buffer,"TISOF",5)==1){
      niso_files++;
      continue;
    }
    if(match(buffer,"ROOM",4) == 1){
      isZoneFireModel=1;
      nrooms++;
      continue;
    }
    if(match(buffer,"FIRE",4) == 1){
      nfires++;
      continue;
    }
    if(match(buffer,"ZONE",4) == 1){
      nzone++;
      continue;
    }
    if(
      match(buffer,"TARG",4) ==1||
      match(buffer,"FTARG",5)==1
      ){
      ntarg_files++;
      continue;
    }
    if(match(buffer,"VENTGEOM",8) == 1||match(buffer,"VFLOWGEOM",9)==1){
      nzvents++;
      continue;
    }

  }

/* 
   ************************************************************************
   ************************ end of pass 1 ********************************* 
   ************************************************************************
 */

 if(npropinfo>0){
   NewMemory((void **)&propinfo,npropinfo*sizeof(propdata));
   npropinfo=0;
 }
 if(nterraininfo>0){
   NewMemory((void **)&terraininfo,nterraininfo*sizeof(terraindata));
   nterraininfo=0;
 }
 if(npartclassinfo>=0){
   float rgb_class[4];
   part5class *partclassi;


   NewMemory((void **)&partclassinfo,(npartclassinfo+1)*sizeof(part5class));

   // define a dummy class

   partclassi = partclassinfo + npartclassinfo;
   strcpy(buffer,"Default");
   trim(buffer);
   len=strlen(buffer);
   partclassi->name=NULL;
   if(len>0){
     NewMemory((void **)&partclassi->name,len+1);
     STRCPY(partclassi->name,trim_front(buffer));
   }

   rgb_class[0]=1.0;
   rgb_class[1]=0.0;
   rgb_class[2]=0.0;
   rgb_class[3]=1.0;
   partclassi->rgb=getcolorptr(rgb_class);

   partclassi->ntypes=0;
   partclassi->xyz=NULL;
   partclassi->maxpoints=0;
   partclassi->labels=NULL;

   NewMemory((void **)&partclassi->labels,sizeof(flowlabels));
   createnulllabel(partclassi->labels);

   npartclassinfo=0;


 }

  ibartemp=2;
  jbartemp=2;
  kbartemp=2;

  /* --------- set up multi-block data structures ------------- */

  /* 
     The keywords TRNX, TRNY, TRNZ, GRID, PDIM, OBST and VENT are not required 
     BUT if any one is present then the number of each must be equal 
  */

  if(nmeshes==0&&ntrnx==0&&ntrny==0&&ntrnz==0&&npdim==0&&nobst==0&&nvent==0&&noffset==0){
    nmeshes=1;ntrnx=1;ntrny=1;ntrnz=1;npdim=1;nobst=1;noffset=1;
  }
  else{
  if(nmeshes>1){
      if(nmeshes!=ntrnx||nmeshes!=ntrny||nmeshes!=ntrnz||
         nmeshes!=npdim||nmeshes!=nobst||nmeshes!=nvent||
         nmeshes!=noffset){
        printf("*** fatal error:\n");
        if(nmeshes!=ntrnx)printf("  found %i TRNX keywords, was expecting %i\n",ntrnx,nmeshes);
        if(nmeshes!=ntrny)printf("  found %i TRNY keywords, was expecting %i\n",ntrny,nmeshes);
        if(nmeshes!=ntrnz)printf("  found %i TRNZ keywords, was expecting %i\n",ntrnz,nmeshes);
        if(nmeshes!=npdim)printf("  found %i PDIM keywords, was expecting %i\n",npdim,nmeshes);
        if(nmeshes!=nobst)printf("  found %i OBST keywords, was expecting %i\n",nobst,nmeshes);
        if(nmeshes!=nvent)printf("  found %i VENT keywords, was expecting %i\n",nvent,nmeshes);
        if(nmeshes!=noffset)printf("  found %i OFFSET keywords, was expecting %i\n",noffset,nmeshes);
        return 2;
      }
     }
  }
  FREEMEMORY(meshinfo);
  if(NewMemory((void **)&meshinfo,nmeshes*sizeof(mesh))==0)return 2;
  meshinfo->plot3dfilenum=-1;
  update_current_mesh(meshinfo);
  for(n=0;n<nmeshes;n++){
    mesh *meshi;

    meshi=meshinfo+n;
    meshi->ibar=0;
    meshi->jbar=0;
    meshi->kbar=0;
    meshi->nbptrs=0;
    meshi->nvents=0;
    meshi->plotn=1;
    meshi->itextureoffset=0;

    meshi->nsmoothblockages_list=0;
    meshi->smoothblockages_list=NULL;
    meshi->nsmoothblockages_list++;
  }
  if(setPDIM==0){
    mesh *meshi;

    if(roomdefined==0){
      xbar0 = 0.0;    xbar = 1.0;
      ybar0 = 0.0;    ybar = 1.0;
      zbar0 = 0.0;    zbar = 1.0; 
    }
    meshi=meshinfo;
    meshi->xbar0=xbar0;
    meshi->xbar =xbar;
    meshi->xcen=(xbar+xbar0)/2.0;
    meshi->ybar0=ybar0;
    meshi->ybar =ybar;
    meshi->ycen=(ybar+ybar0)/2.0;
    meshi->zbar0=zbar0;
    meshi->zbar =zbar;
    meshi->zcen=(zbar+zbar0)/2.0;
  }

  // define labels and memory for default colorbars

  FREEMEMORY(partinfo);
  if(npart_files!=0){
    if(NewMemory((void **)&partinfo,npart_files*sizeof(particle))==0)return 2;
  }

  FREEMEMORY(targinfo);
  if(ntarg_files!=0){
    if(NewMemory((void **)&targinfo,ntarg_files*sizeof(targ))==0)return 2;
  }

  FREEMEMORY(vsliceinfo);
  FREEMEMORY(sliceinfo);
  FREEMEMORY(slicetypes);
  FREEMEMORY(vslicetypes);
  if(nslice_files>0){
    if(NewMemory( (void **)&vsliceinfo, 3*nslice_files*sizeof(vslice) )==0||
       NewMemory( (void **)&sliceinfo,  nslice_files*sizeof(slice)    )==0||
       NewMemory( (void **)&slicetypes, nslice_files*sizeof(int)      )==0||
       NewMemory( (void **)&vslicetypes,3*nslice_files*sizeof(int)    )==0){return 2;}
  }
  if(nsmoke3d_files>0){
    if(NewMemory( (void **)&smoke3dinfo, nsmoke3d_files*sizeof(smoke3d))==0)return 2;
  }

#ifndef pp_OSX
//  updatecolors(-1);
#endif

  FREEMEMORY(plot3dinfo);
  if(nplot3d_files>0){
    if(NewMemory((void **)&plot3dinfo,nplot3d_files*sizeof(plot3d))==0)return 2;
  }
  FREEMEMORY(patchinfo);
  FREEMEMORY(patchtypes);
  if(npatch_files!=0){
    if(NewMemory((void **)&patchinfo,npatch_files*sizeof(patch))==0)return 2;
    for(i=0;i<npatch_files;i++){
      patch *patchi;

      patchi = patchinfo + i;
      patchi->reg_file=NULL;
      patchi->comp_file=NULL;
      patchi->file=NULL;
      patchi->size_file=NULL;
    }
    if(NewMemory((void **)&patchtypes,npatch_files*sizeof(int))==0)return 2;
  }
  FREEMEMORY(isoinfo);
  FREEMEMORY(isotypes);
  if(niso_files>0){
    if(NewMemory((void **)&isoinfo,niso_files*sizeof(iso))==0)return 2;
    if(NewMemory((void **)&isotypes,niso_files*sizeof(int))==0)return 2;
  }
  FREEMEMORY(roominfo);
  if(nrooms>0){
    if(NewMemory((void **)&roominfo,(nrooms+1)*sizeof(roomdata))==0)return 2;
  }
  FREEMEMORY(fireinfo);
  if(nfires>0){
    if(NewMemory((void **)&fireinfo,nfires*sizeof(firedata))==0)return 2;
  }
  FREEMEMORY(zoneinfo);
  if(nzone>0){
    if(NewMemory((void **)&zoneinfo,nzone*sizeof(zone))==0)return 2;
  }
  FREEMEMORY(zventinfo);
  if(nzvents>0){
    if(NewMemory((void **)&zventinfo,nzvents*sizeof(zvent))==0)return 2;
  }
  nzvents=0;

  FREEMEMORY(textureinfo);
  FREEMEMORY(surfaceinfo);
  if(nsurfaces>0){
    if(NewMemory((void **)&surfaceinfo,nsurfaces*sizeof(surface))==0)return 2;
  }

  if(cadgeominfo!=NULL)freecadinfo();
  if(ncadgeom>0){
    if(NewMemory((void **)&cadgeominfo,ncadgeom*sizeof(cadgeom))==0)return 2;
  }

  if(noutlineinfo>0){
    if(NewMemory((void **)&outlineinfo,noutlineinfo*sizeof(outline))==0)return 2;
    for(i=0;i<noutlineinfo;i++){
      outlinei = outlineinfo + i;
      outlinei->x1=NULL;
      outlinei->x2=NULL;
      outlinei->y1=NULL;
      outlinei->y2=NULL;
      outlinei->z1=NULL;
      outlinei->z2=NULL;
    }
  }
  if(nticks>0){
    if(NewMemory((void **)&tickinfo,nticks*sizeof(tickdata))==0)return 2;
    nticks=0;
    ntickssmv=0;
  }

  if(nlabels>0){
    if(NewMemory((void **)&labelinfo,nlabels*sizeof(labeldata))==0)return 2;
    nlabels=0;
    nlabelssmv=0;
  }
  if(ndeviceinfo>0){
    if(NewMemory((void **)&deviceinfo,ndeviceinfo*sizeof(device))==0)return 2;
    devicecopy=deviceinfo;
  }

  // read in device (.svo) definitions

    init_object_defs();


/* 
   ************************************************************************
   ************************ start of pass 2 ********************************* 
   ************************************************************************
 */

  startpass=1;
  ioffset=0;
  iobst=0;
  ncadgeom=0;
  nsurfaces=0;
  ndeviceinfo=0;
  noutlineinfo=0;
  if(noffset==0)ioffset=1;
  rewind(stream);
  printf("   pass 1 completed\n");
  printf("   pass 2\n");
  while(!feof(stream)){

    if(noGRIDpresent==1&&startpass==1){
      strcpy(buffer,"GRID");
      startpass=0;
    }
    else{
      if(fgets(buffer,255,stream)==NULL)break;
      if(strncmp(buffer," ",1)==0||buffer[0]==0)continue;
    }

    /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++ PROP ++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

    if(match(buffer,"PROP",4) == 1){
      propdata *propi;
      char *fbuffer;
      char proplabel[255];
      int lenbuf;
      int ntextures;
      int nsmokeview_ids;
      char *smokeview_id;

      propi = propinfo + npropinfo;

      if(fgets(proplabel,255,stream)==NULL)break; // prop label
      trim(proplabel);
      fbuffer=trim_front(proplabel);

      if(fgets(buffer,255,stream)==NULL)break; // number of smokeview_id's
      sscanf(buffer,"%i",&nsmokeview_ids);

      init_prop(propi,nsmokeview_ids,fbuffer);
      for(i=0;i<nsmokeview_ids;i++){
        if(fgets(buffer,255,stream)==NULL)break; // smokeview_id
        trim(buffer);
        fbuffer=trim_front(buffer);
        lenbuf=strlen(fbuffer);
        NewMemory((void **)&smokeview_id,lenbuf+1);
        strcpy(smokeview_id,fbuffer);
        propi->smokeview_ids[i]=smokeview_id;
        propi->smv_objects[i]=get_SVOBJECT_type(propi->smokeview_ids[i],missing_device);
      }
      propi->smv_object=propi->smv_objects[0];
      propi->smokeview_id=propi->smokeview_ids[0];

      if(fgets(buffer,255,stream)==NULL)break; // keyword_values
      sscanf(buffer,"%i",&propi->nvars_indep);
      propi->vars_indep=NULL;
      propi->svals=NULL;
      propi->texturefiles=NULL;
      ntextures=0;
      if(propi->nvars_indep>0){
        NewMemory((void **)&propi->vars_indep,propi->nvars_indep*sizeof(char *));
        NewMemory((void **)&propi->svals,propi->nvars_indep*sizeof(char *));
        NewMemory((void **)&propi->fvals,propi->nvars_indep*sizeof(float));
        NewMemory((void **)&propi->vars_indep_index,propi->nvars_indep*sizeof(int));
        NewMemory((void **)&propi->texturefiles,propi->nvars_indep*sizeof(char *));

        for(i=0;i<propi->nvars_indep;i++){
          char *equal;

          propi->svals[i]=NULL;
          propi->vars_indep[i]=NULL;
          propi->fvals[i]=0.0;
          fgets(buffer,255,stream);
          equal=strchr(buffer,'=');
          if(equal!=NULL){
            char *buf1, *buf2, *keyword, *val;
            int lenkey, lenval;
            char *texturefile;

            buf1=buffer;
            buf2=equal+1;
            *equal=0;

            trim(buf1);
            keyword=trim_front(buf1);
            lenkey=strlen(keyword);

            trim(buf2);
            val=trim_front(buf2);
            lenval=strlen(val);

            if(lenkey==0||lenval==0)continue;

            if(val[0]=='"'){
              val[0]=' ';
              if(val[lenval-1]=='"')val[lenval-1]=' ';
              trim(val);
              val=trim_front(val);
              NewMemory((void **)&propi->svals[i],lenval+1);
              strcpy(propi->svals[i],val);
              texturefile=strstr(val,"t%");
              if(texturefile!=NULL){
                int lentexture;

                texturefile+=2;
                texturefile=trim_front(texturefile);
                propi->texturefiles[ntextures]=propi->svals[i];
                strcpy(propi->svals[i],texturefile);

                ntextures++;
              }
            }

            NewMemory((void **)&propi->vars_indep[i],lenkey+1);
            strcpy(propi->vars_indep[i],keyword);

            sscanf(val,"%f",propi->fvals+i);
          }
        }
        get_indep_var_indices(propi->smv_object,
          propi->vars_indep,propi->nvars_indep,
          propi->vars_indep_index);

        get_evac_indices(propi->smv_object,
          propi->vars_evac_index,&propi->nvars_evac);

      }
      propi->ntextures=ntextures;
      npropinfo++;
      continue;
    }

    /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++ TERRAIN +++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

    if(match(buffer,"TERRAIN",7) == 1){
//      terraindata *terri;
      float xmin, xmax, ymin, ymax;
      int nx, ny;

      manual_terrain=1;

  //    terri = terraininfo + nterraininfo;

      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %i %f %f %i",&xmin, &xmax, &nx, &ymin, &ymax, &ny);
      // must implement new form for defining terrain surfaces
      //initterrain(stream, NULL, terri, xmin, xmax, nx, ymin, ymax, ny);

      nterraininfo++;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++ CLASS_OF_PARTICLES +++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"CLASS_OF_PARTICLES",18) == 1||
       match(buffer,"CLASS_OF_HUMANS",15) == 1){
      float rgb_class[4];
      part5class *partclassi;
      char *device_ptr;
      char *smokeview_id, *prop_id;

      partclassi = partclassinfo + npartclassinfo;
      partclassi->kind=PARTICLES;
      if(match(buffer,"CLASS_OF_HUMANS",15) == 1)partclassi->kind=HUMANS;
      fgets(buffer,255,stream);

      get_labels(buffer,&device_ptr,&prop_id);
      if(prop_id!=NULL){
        device_ptr=NULL;
      }
      partclassi->prop=NULL;

      partclassi->sphere=NULL;
      partclassi->smv_device=NULL;
      partclassi->device_name=NULL;
      if(device_ptr!=NULL){
        partclassi->sphere=get_SVOBJECT_type("SPHERE",missing_device);

        partclassi->smv_device=get_SVOBJECT_type(device_ptr,missing_device);
        if(partclassi->smv_device!=NULL){
          len = strlen(device_ptr);
          NewMemory((void **)&partclassi->device_name,len+1);
          STRCPY(partclassi->device_name,device_ptr);
        }
        else{
          char tube[10];

          strcpy(tube,"TUBE");
          len = strlen(tube);
          NewMemory((void **)&partclassi->device_name,len+1);
          STRCPY(partclassi->device_name,tube);
          partclassi->smv_device=get_SVOBJECT_type(tube,missing_device);
        }
      }

      trim(buffer);
      len=strlen(buffer);
      partclassi->name=NULL;
      if(len>0){
        NewMemory((void **)&partclassi->name,len+1);
        STRCPY(partclassi->name,trim_front(buffer));
      }

      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",rgb_class,rgb_class+1,rgb_class+2);
      rgb_class[3]=1.0;
      partclassi->rgb=getcolorptr(rgb_class);

      partclassi->ntypes=0;
      partclassi->xyz=NULL;
      partclassi->maxpoints=0;
      partclassi->labels=NULL;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&partclassi->ntypes);
      partclassi->ntypes+=2;
      partclassi->nvars_dep=partclassi->ntypes-2+3; // subtract off two "dummies" at beginning and add 3 at end for r,g,b
      if(partclassi->ntypes>0){
        flowlabels *labelj;
        char shortdefaultlabel[]="Uniform";
        char longdefaultlabel[]="Uniform color";

        NewMemory((void **)&partclassi->labels,partclassi->ntypes*sizeof(flowlabels));
       
        labelj = partclassi->labels+0; // placeholder for hidden

        labelj->longlabel=NULL;
        labelj->shortlabel=NULL;
        labelj->unit=NULL;

        labelj = partclassi->labels+1;  // placeholder for default

        labelj->longlabel=NULL;
        NewMemory((void **)&labelj->longlabel,strlen(longdefaultlabel)+1);
        strcpy(labelj->longlabel,longdefaultlabel);
        labelj->shortlabel=NULL;
        NewMemory((void **)&labelj->shortlabel,strlen(shortdefaultlabel)+1);
        strcpy(labelj->shortlabel,shortdefaultlabel);
        labelj->unit=NULL;

        partclassi->col_azimuth=-1;
        partclassi->col_diameter=-1;
        partclassi->col_elevation=-1;
        partclassi->col_length=-1;
        partclassi->col_u_vel=-1;
        partclassi->col_v_vel=-1;
        partclassi->col_w_vel=-1;
        partclassi->vis_type=PART_POINTS;
        for(j=2;j<partclassi->ntypes;j++){
          labelj = partclassi->labels+j;
          labelj->longlabel=NULL;
          labelj->shortlabel=NULL;
          labelj->unit=NULL;
          readlabels(labelj,stream);
          partclassi->vars_dep[j-2]=labelj->shortlabel;
          if(strcmp(labelj->shortlabel,"DIAMETER")==0){
            partclassi->col_diameter=j-2;
          }
          if(strcmp(labelj->shortlabel,"LENGTH")==0){
            partclassi->col_length=j-2;
          }
          if(strcmp(labelj->shortlabel,"AZIMUTH")==0){
            partclassi->col_azimuth=j-2;
          }
          if(strcmp(labelj->shortlabel,"ELEVATION")==0){
            partclassi->col_elevation=j-2;
          }
          if(strcmp(labelj->shortlabel,"U-VEL")==0){
            partclassi->col_u_vel=j-2;
          }
          if(strcmp(labelj->shortlabel,"V-VEL")==0){
            partclassi->col_v_vel=j-2;
          }
          if(strcmp(labelj->shortlabel,"W-VEL")==0){
            partclassi->col_w_vel=j-2;
          }
        }
      }
      partclassi->diameter=1.0;
      partclassi->length=1.0;
      partclassi->azimuth=0.0;
      partclassi->elevation=0.0;
      partclassi->dx=0.0;
      partclassi->dy=0.0;
      partclassi->dz=0.0;
      if(device_ptr!=NULL){
        float diameter, length, azimuth, elevation;

        fgets(buffer,255,stream);
        sscanf(buffer,"%f %f %f %f",&diameter,&length,&azimuth,&elevation);
        partclassi->diameter=diameter;
        partclassi->length=length;
        partclassi->azimuth=azimuth;
        partclassi->elevation=elevation;
      }
      npartclassinfo++;
      continue;
    }

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ LABEL ++++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

    if(match(buffer,"LABEL",5) == 1){
      nlabels++;
      nlabelssmv++;

      /*
      LABEL
      x y z r g b tstart tstop  
      label

      */
      {
        float *xyz, *rgbtemp, *tstart_stop;
        labeldata *labeli;

        labeli = labelinfo + nlabels-1;

        xyz = labeli->xyz;
        rgbtemp = labeli->rgb;
        tstart_stop = labeli->tstart_stop;

        fgets(buffer,255,stream);
        rgbtemp[0]=-1.0;
        rgbtemp[1]=-1.0;
        rgbtemp[2]=-1.0;
        rgbtemp[3]=1.0;
        tstart_stop[0]=-1.0;
        tstart_stop[1]=-1.0;
        sscanf(buffer,"%f %f %f %f %f %f %f %f",
          xyz,xyz+1,xyz+2,
          rgbtemp,rgbtemp+1,rgbtemp+2,
          tstart_stop,tstart_stop+1);
        if(rgbtemp[0]<0.0||rgbtemp[1]<0.0||rgbtemp[2]<0.0||rgbtemp[0]>1.0||rgbtemp[1]>1.0||rgbtemp[2]>1.0){
          labeli->useforegroundcolor=1;
        }
        else{
          labeli->useforegroundcolor=0;
        }
        fgets(buffer,255,stream);
        strcpy(labeli->label,buffer);
      }
      continue;
    }

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ TICKS ++++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

/*
typedef struct {
  float begin[3],end[3],length;
  float dxyz[3],dlength;
  int dir,nbars;
} tickdata;
*/

    if(match(buffer,"TICKS",5) == 1){
      nticks++;
      ntickssmv++;
      {
        tickdata *ticki;
        float *begt, *endt;
        int *nbarst;
        float term;
        float length=0.0;
        float *dxyz;
        float sum;

        ticki = tickinfo + nticks - 1;
        begt = ticki->begin;
        endt = ticki->end;
        nbarst=&ticki->nbars;
        dxyz = ticki->dxyz;
        

        /*
        TICKS
        b1 b2 b3 e1 e2 e3 nb
        ticklength tickdir tickcolor (r g b) tickwidth 
        */
        if(fgets(buffer,255,stream)==NULL)break;
        *nbarst=0;
        sscanf(buffer,"%f %f %f %f %f %f %i",begt,begt+1,begt+2,endt,endt+1,endt+2,nbarst);
        if(*nbarst<1)*nbarst=1;
        if(fgets(buffer,255,stream)==NULL)break;
        {
          float *rgbtemp;

          rgbtemp=ticki->rgb;
          rgbtemp[0]=-1.0;
          rgbtemp[1]=-1.0;
          rgbtemp[2]=-1.0;
          ticki->width=-1.0;
          sscanf(buffer,"%f %i %f %f %f %f",&ticki->dlength,&ticki->dir,rgbtemp,rgbtemp+1,rgbtemp+2,&ticki->width);
          if(rgbtemp[0]<0.0||rgbtemp[0]>1.0||
             rgbtemp[1]<0.0||rgbtemp[1]>1.0||
             rgbtemp[2]<0.0||rgbtemp[2]>1.0){
            ticki->useforegroundcolor=1;
          }
          else{
            ticki->useforegroundcolor=0;
          }
          if(ticki->width<0.0)ticki->width=1.0;
        }
        for(i=0;i<3;i++){
          term = endt[i]-begt[i];
          length += term*term;
        }
        if(length<=0.0){
          endt[0]=begt[0]+1.0;
          length = 1.0;
        }
        ticki->length=sqrt(length);
        dxyz[0] =  0.0;
        dxyz[1] =  0.0;
        dxyz[2] =  0.0;
        switch (ticki->dir){
        case 1:
        case -1:
          dxyz[0]=1.0;
          break;
        case 2:
        case -2:
          dxyz[1]=1.0;
          break;
        case 3:
        case -3:
          dxyz[2]=1.0;
          break;
        default:
          ASSERT(FFALSE);
          break;
        }
        if(ticki->dir<0){
          for(i=0;i<3;i++){
            dxyz[i]=-dxyz[i];
          }
        }
        sum = 0.0;
        sum = dxyz[0]*dxyz[0] + dxyz[1]*dxyz[1] + dxyz[2]*dxyz[2];
        if(sum>0.0){
          sum=sqrt(sum);
          dxyz[0] *= (ticki->dlength/sum);
          dxyz[1] *= (ticki->dlength/sum);
          dxyz[2] *= (ticki->dlength/sum);
        }
      }
      continue;
    }


  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ OUTLINE ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

    if(match(buffer,"OUTLINE",7) == 1){
      noutlineinfo++;
      outlinei = outlineinfo + noutlineinfo - 1;
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%i",&outlinei->nlines);
      if(outlinei->nlines>0){
        x1=NULL;
        x2=NULL;
        yy1=NULL;
        yy2=NULL;
        z1=NULL;
        z2=NULL;
        NewMemory((void **)&outlinei->x1,outlinei->nlines*sizeof(float));
        NewMemory((void **)&outlinei->y1,outlinei->nlines*sizeof(float));
        NewMemory((void **)&outlinei->z1,outlinei->nlines*sizeof(float));
        NewMemory((void **)&outlinei->x2,outlinei->nlines*sizeof(float));
        NewMemory((void **)&outlinei->y2,outlinei->nlines*sizeof(float));
        NewMemory((void **)&outlinei->z2,outlinei->nlines*sizeof(float));
        for(i=0;i<outlinei->nlines;i++){
          fgets(buffer,255,stream);
          sscanf(buffer,"%f %f %f %f %f %f",
            outlinei->x1+i,outlinei->y1+i, outlinei->z1+i,
            outlinei->x2+i,outlinei->y2+i, outlinei->z2+i
            );
        }
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ CADGEOM ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"CADGEOM",7) == 1){
      if(fgets(buffer,255,stream)==NULL)break;
      bufptr=trim_string(buffer);
      len=strlen(bufptr);
      cadgeominfo[ncadgeom].order=NULL;
      cadgeominfo[ncadgeom].quad=NULL;
      cadgeominfo[ncadgeom].file=NULL;
      if(STAT(buffer,&statbuffer)==0){
        if(NewMemory((void **)&cadgeominfo[ncadgeom].file,(unsigned int)(len+1))==0)return 2;
        STRCPY(cadgeominfo[ncadgeom].file,bufptr);
        printf("   reading cad file: %s\n",bufptr);
        readcadgeom(cadgeominfo+ncadgeom);
        printf("   completed\n");
        ncadgeom++;
      }
      else{
        printf("   CAD geometry file: %s could not be opened\n",bufptr);
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ OFFSET ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"OFFSET",6) == 1){
      ioffset++;
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SURFDEF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"SURFDEF",7) == 1){
      fgets(buffer,255,stream);
      bufptr=trim_string(buffer);
      strcpy(surfacedefaultlabel,trim_front(bufptr));
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SMOKE3D ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"SMOKE3D",7) == 1){
      size_t lenbuffer;

      nn_smoke3d++;
      trim(buffer);
      len=strlen(buffer);
      if(nmeshes>1){
        blocknumber=ioffset-1;
      }
      else{
        blocknumber=0;
      }
      if(len>8){
        buffer3=buffer+8;
        sscanf(buffer3,"%i",&blocknumber);
        blocknumber--;
      }
      if(fgets(buffer,255,stream)==NULL){
        nsmoke3d_files--;
        break;
      }
      bufptr=trim_string(buffer);
      len=strlen(buffer);
      lenbuffer=len;
      {
        smoke3d *smoke3di;
        
        smoke3di = smoke3dinfo + ismoke3d;

        if(NewMemory((void **)&smoke3di->reg_file,(unsigned int)(len+1))==0)return 2;
        STRCPY(smoke3di->reg_file,bufptr);

#ifdef pp_LIGHT
        smoke3di->use_lighting_file=0;
        smoke3di->light_file=NULL;
#endif
        smoke3di->seq_id=nn_smoke3d;
        smoke3di->autoload=0;
        smoke3di->version=-1;
        smoke3di->hrrpuv_color=NULL;
        smoke3di->water_color=NULL;
        smoke3di->soot_color=NULL;
        smoke3di->file=NULL;
        smoke3di->smokeframe_in=NULL;
        smoke3di->smokeframe_comp_list=NULL;
        smoke3di->smokeframe_out=NULL;
        smoke3di->timeslist=NULL;
        smoke3di->smoke_comp_all=NULL;
        smoke3di->smokeview_tmp=NULL;
        smoke3di->times=NULL;
        smoke3di->use_smokeframe=NULL;
        smoke3di->nchars_compressed_smoke=NULL;
        smoke3di->nchars_compressed_smoke_full=NULL,
#ifdef pp_LIGHT
        smoke3di->nchars_compressed_light=NULL;
        smoke3di->light_comp_all=NULL;
        smoke3di->lightframe_in=NULL;
        smoke3di->lightframe_out=NULL;
        smoke3di->lightframe_comp_list=NULL;
        smoke3di->lightview_tmp=NULL;
#endif

        smoke3di->display=0;
        smoke3di->loaded=0;
        smoke3di->d_display=0;
        smoke3di->blocknumber=blocknumber;
        smoke3di->lastiframe=-999;
        smoke3di->soot_index=-1;
        smoke3di->water_index=-1;
        smoke3di->hrrpuv_index=-1;

        STRCPY(buffer2,bufptr);
        STRCAT(buffer2,".svz");

        len=lenbuffer+4;
        if(NewMemory((void **)&smoke3di->comp_file,(unsigned int)(len+1))==0)return 2;
        STRCPY(smoke3di->comp_file,buffer2);

        if(STAT(smoke3di->comp_file,&statbuffer)==0){
          smoke3di->file=smoke3di->comp_file;
        }
        else{
          smoke3di->file=smoke3di->reg_file;
        }
        if(STAT(smoke3di->file,&statbuffer)==0){
          if(readlabels(&smoke3di->label,stream)==2)return 2;
          if(strcmp(smoke3di->label.longlabel,"HRRPUV")==0){
            show_hrrcutoff_active=1;
          }
          ismoke3d++;
        }
        else{
          if(readlabels(&smoke3di->label,stream)==2)return 2;
          nsmoke3d_files--;
        }
        smoke3di->version=getsmoke3dversion(smoke3di);
        if(strncmp(smoke3di->label.shortlabel,"soot",4)==0){
          smoke3di->type=1;
        }
        else if(strncmp(smoke3di->label.shortlabel,"hrrpuv",6)==0){
          smoke3di->type=2;
        }
        else if(strncmp(smoke3di->label.shortlabel,"water",5)==0){
          smoke3di->type=3;
        }
        else{
          smoke3di->type=1;
        }
#ifdef pp_LIGHT
        STRCPY(buffer2,buffer);
        STRCAT(buffer2,".lvz");

        len=lenbuffer+4;
        if(NewMemory((void **)&smoke3di->light_file,(unsigned int)(len+1))==0)return 2;
        STRCPY(smoke3di->light_file,buffer2);

        if(STAT(smoke3di->light_file,&statbuffer)==0){
          FILE *light_stream;

          light_stream=fopen(smoke3di->light_file,"rb");
          if(light_stream!=NULL){
            fclose(light_stream);
            smoke3di->use_lighting_file=1;
          }
        }

#endif

      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SURFACE ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"SURFACE",7) ==1&&match(buffer,"SURFACE DENSITY",15)!=1){
      surfi = surfaceinfo + nsurfaces;
      initsurface(surfi);
      fgets(buffer,255,stream);
      trim(buffer);
      len=strlen(buffer);
      NewMemory((void **)&surfi->surfacelabel,(len+1)*sizeof(char));
      strcpy(surfi->surfacelabel,trim_front(buffer));
      if(strstr(surfi->surfacelabel,"MIRROR")!=NULL||
         strstr(surfi->surfacelabel,"INTERPOLATED")!=NULL||
         strstr(surfi->surfacelabel,"OPEN")!=NULL){
        surfi->obst_surface=0;
      }
      if(strstr(surfi->surfacelabel,"INERT")!=NULL){
        surfaceinfo[0].obst_surface=1;
      }

      temp_ignition=TEMP_IGNITION_MAX;
      emis = 1.0;
      t_width=1.0;
      t_height=1.0;
      s_type=0;
      s_color[0]=surfi->color[0];
      s_color[1]=surfi->color[1];
      s_color[2]=surfi->color[2];
      //s_color[3]=1.0;
      s_color[3]=surfi->color[3];
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f",&temp_ignition,&emis);
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %f %f %f %f %f",
        &s_type,&t_width, &t_height,s_color,s_color+1,s_color+2,s_color+3);

      surfi->type=s_type;
      if(s_color[0]<=0.0&&s_color[1]<=0.0&&s_color[2]<=0.0){
        if(s_color[0]!=0.0||s_color[1]!=0.0||s_color[2]!=0.0){
          surfi->invisible=1;
          surfi->type=BLOCK_hidden;
          s_color[0]=-s_color[0];
          s_color[1]=-s_color[1];
          s_color[2]=-s_color[2];
        }
      }
      surfi->color = getcolorptr(s_color);
      if(s_color[3]<0.99){
        surfi->transparent=1;
      }
      surfi->temp_ignition=temp_ignition;
      surfi->emis=emis;
      surfi->t_height=t_height;
      surfi->t_width=t_width;
      surfi->texturefile=NULL;
      surfi->textureinfo=NULL;

      fgets(buffer,255,stream);
      trim(buffer);
      buffer3 = trim_front(buffer);
#ifdef pp_RENDER
      {
        int found_texture;
        char texturebuffer[1024];

        found_texture=0;
        if(smokeviewbindir!=NULL&&STAT(buffer3,&statbuffer)!=0){
          STRCPY(texturebuffer,smokeviewbindir);
          STRCAT(texturebuffer,buffer3);
          if(STAT(texturebuffer,&statbuffer)==0){
            if(NewMemory((void **)&surfi->texturefile,strlen(texturebuffer)+1)==0)return 2;
            STRCPY(surfi->texturefile,texturebuffer);
            found_texture=1;
          }
        }
        if(STAT(buffer3,&statbuffer)==0){
          len=strlen(buffer3);
          if(NewMemory((void **)&surfi->texturefile,(unsigned int)(len+1))==0)return 2;
          STRCPY(surfi->texturefile,buffer3);
          found_texture=1;
        }
        if(texturebuffer!=NULL&&buffer3!=NULL&&found_texture==0&&strncmp(buffer3,"null",4)!=0){
          printf("*** Warning: The texture file: %s was not found in either \n",buffer3);
          printf("             the current working directory or in %s\n",smokeviewbindir);
        }
      }
#endif
      nsurfaces++;
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ GRID ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"GRID",4) == 1){
      mesh *meshi;
      int mesh_type=0;

//      int lenbuffer;

//      trim(buffer);
//      lenbuffer=strlen(buffer);
//      if(lenbuffer>4){
//        if(buffer[5]!=' ')continue;
//      }

      igrid++;
      if(meshinfo!=NULL){
        meshi=meshinfo+igrid-1;
        initmesh(meshi);
        {
          size_t len_meshlabel;
          char *meshlabel;

          len_meshlabel=0;
          if(strlen(buffer)>5){
            meshlabel=trim_front(buffer+5);
            trim(meshlabel);
            len_meshlabel=strlen(meshlabel);
          }
          if(len_meshlabel>0){
            NewMemory((void **)&meshi->label,(len_meshlabel+1));
            strcpy(meshi->label,meshlabel);
          }
          else{
            sprintf(buffer,"%i",igrid);
            NewMemory((void **)&meshi->label,strlen(buffer)+1);
            strcpy(meshi->label,buffer);
          }
        }

      }
      setGRID=1;
      if(noGRIDpresent==1){
        ibartemp=2;
        jbartemp=2;
        kbartemp=2;
      }
      else{
        fgets(buffer,255,stream);
        sscanf(buffer,"%i %i %i %i",&ibartemp,&jbartemp,&kbartemp,&mesh_type);
      }
      if(ibartemp<1)ibartemp=1;
      if(jbartemp<1)jbartemp=1;
      if(kbartemp<1)kbartemp=1;
      xp=NULL; yp=NULL; zp=NULL;
      xp2=NULL; yp2=NULL; zp2=NULL;
      if(NewMemory((void **)&xp,sizeof(float)*(ibartemp+1))==0||
         NewMemory((void **)&yp,sizeof(float)*(jbartemp+1))==0||
         NewMemory((void **)&zp,sizeof(float)*(kbartemp+1))==0||
         NewMemory((void **)&xp2,sizeof(float)*(ibartemp+1))==0||
         NewMemory((void **)&yp2,sizeof(float)*(jbartemp+1))==0||
         NewMemory((void **)&zp2,sizeof(float)*(kbartemp+1))==0
         )return 2;
      if(meshinfo!=NULL){
        meshi->mesh_type=mesh_type;
        meshi->xplt=xp;
        meshi->yplt=yp;
        meshi->zplt=zp;
        meshi->xplt_orig=xp2;
        meshi->yplt_orig=yp2;
        meshi->zplt_orig=zp2;
        meshi->ibar=ibartemp;
        meshi->jbar=jbartemp;
        meshi->kbar=kbartemp;
        meshi->plotx=ibartemp/2;
        meshi->ploty=jbartemp/2;
        meshi->plotz=kbartemp/2;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ ZONE ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"ZONE",4) == 1){
      char *last_name,*full_name,*filename;
      zone *zonei;

      zonei = zoneinfo + izone;
      if(fgets(buffer,255,stream)==NULL){
        nzone--;
        break;
      }
      bufptr=trim_string(buffer);
      len=strlen(bufptr);
      zonei->loaded=0;
      zonei->display=0;

      full_name=bufptr;
      if(STAT(full_name,&statbuffer)!=0)full_name=NULL;

      last_name=lastname(buffer);
      if(STAT(last_name,&statbuffer)!=0)last_name=NULL;

      if(full_name!=NULL&&last_name!=NULL){
        if(strcmp(last_name,full_name)==0){
          last_name=NULL;
        }
      }



      if(NewMemory((void **)&zonei->file,(unsigned int)(len+1))==0)return 2;
      if(last_name!=NULL&&full_name!=NULL){
        filename=last_name;
      }
      else if(last_name==NULL&&full_name!=NULL){
        filename=full_name;
      }
      else if(last_name!=NULL&&full_name==NULL){
        filename=last_name;
      }
      else{
        filename=NULL;
      }
      if(filename!=NULL)STRCPY(zonei->file,filename);

      if(filename==NULL){
        for(n=0;n<4;n++){
          if(readlabels(&zonei->label[n],stream)==2){
            return 2;
          }
        }
        nzone--;
      }
      else{
        for(n=0;n<4;n++){
          if(readlabels(&zonei->label[n],stream)==2){
            return 2;
          }
        }
        izone++;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ AMBIENT ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"AMBIENT",7)==1){
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%f %f %f",&pref,&pamb,&tamb);
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ ROOM ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"ROOM",4) == 1){
      roomdata *roomi;

      isZoneFireModel=1;
      visFrame=0;
      roomdefined=1;
      iroom++;
      roomi = roominfo + iroom - 1;
      roomi->valid=0;
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%f %f %f",&roomi->dx,&roomi->dy,&roomi->dz);
      roomi->valid=1;
      if(fgets(buffer,255,stream)==NULL){
        roomi->x0=0.0;
        roomi->y0=0.0;
        roomi->z0=0.0;
      }
      else{
        sscanf(buffer,"%f %f %f",&roomi->x0,&roomi->y0,&roomi->z0);
      }
      roomi->x1=roomi->x0+roomi->dx;
      roomi->y1=roomi->y0+roomi->dy;
      roomi->z1=roomi->z0+roomi->dz;

      if(setPDIM==0){
        if(roomi->x0<xbar0)xbar0=roomi->x0;
        if(roomi->y0<ybar0)ybar0=roomi->y0;
        if(roomi->z0<zbar0)zbar0=roomi->z0;
        if(roomi->x1>xbar)xbar=roomi->x1;
        if(roomi->y1>ybar)ybar=roomi->y1;
        if(roomi->z1>zbar)zbar=roomi->z1;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ DEVICE +++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    DEVICE
    label
    x y z xn yn zn state nparams ntextures 
    p0 p1 ... p5
    p6 ...    p11
    texturefile1
    ...
    texturefilen

    */
    if(
      (match(buffer,"DEVICE",6) == 1)&&
      (match(buffer,"DEVICE_ACT",10) != 1)
      ){
      device *devicei;
      float xyz[3]={0.0,0.0,0.0}, xyzn[3]={0.0,0.0,0.0};
      int state0=0;
      int nparams=0, nparams_textures=0;
      char *labelptr, *prop_id;

      devicei = deviceinfo + ndeviceinfo;
      devicei->type=DEVICE_DEVICE;
      fgets(buffer,255,stream);
      trim(buffer);
      strcpy(devicei->label,trim_front(buffer));
      devicei->object = get_SVOBJECT_type(buffer,missing_device);
      devicei->params=NULL;
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f %f %f %f %i %i %i",
        xyz,xyz+1,xyz+2,xyzn,xyzn+1,xyzn+2,&state0,&nparams,&nparams_textures);
      get_labels(buffer,&prop_id,NULL);
      devicei->prop=get_prop_id(prop_id);
      if(prop_id!=NULL&&devicei->prop!=NULL&&devicei->prop->smv_object!=NULL){
        devicei->object=devicei->prop->smv_object;
      }
      else{
        NewMemory((void **)&devicei->prop,sizeof(propdata));
        init_prop(devicei->prop,1,devicei->label);
        devicei->prop->smv_object=devicei->object;
        devicei->prop->smv_objects[0]=devicei->prop->smv_object;
      }
      if(nparams_textures<0)nparams_textures=0;
      if(nparams_textures>1)nparams_textures=1;
      devicei->ntextures=nparams_textures;
      if(nparams_textures>0){
         NewMemory((void **)&devicei->textureinfo,sizeof(texture));
      }
      else{
        devicei->textureinfo=NULL;
        devicei->texturefile=NULL;
      }

      labelptr=strchr(buffer,'%');
      if(labelptr!=NULL){
        trim(labelptr);
        if(strlen(labelptr)>1){
          labelptr++;
          labelptr=trim_front(labelptr);
          if(strlen(labelptr)==0)labelptr=NULL;
        }
        else{
          labelptr=NULL;
        }
      }

      if(nparams<=0){
        init_device(devicei,xyz,xyzn,state0,0,NULL,labelptr);
      }
      else{
        float *params,*pc;
        int nsize;

        nsize = 6*((nparams-1)/6+1);
        NewMemory((void **)&params,(nsize+devicei->ntextures)*sizeof(float));
        pc=params;
        for(i=0;i<nsize/6;i++){
          fgets(buffer,255,stream);
          sscanf(buffer,"%f %f %f %f %f %f",pc,pc+1,pc+2,pc+3,pc+4,pc+5);
          pc+=6;
        }
        init_device(devicei,xyz,xyzn,state0,nparams,params,labelptr);
      }
      get_elevaz(devicei->xyznorm,&devicei->dtheta,devicei->rotate_axis);
      if(nparams_textures>0){
        fgets(buffer,255,stream);
        trim(buffer);
        buffer3=trim_front(buffer);
        NewMemory((void **)&devicei->texturefile,strlen(buffer3)+1);
        strcpy(devicei->texturefile,buffer3);
      }

      CheckMemory;
      ndeviceinfo++;
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ DEVICE_ACT ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"DEVICE_ACT",10) == 1){
      device *devicei;
      int idevice;
      float act_time;
      int act_state;

      do_pass4=1;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i",&idevice,&act_time,&act_state);
      idevice--;
      if(idevice>=0&&idevice<ndeviceinfo){
        devicei = deviceinfo + idevice;
        devicei->act_time=act_time;
        devicei->nstate_changes++;
      }

      continue;
    }
  }
/* 
   ************************************************************************
   ************************ end of pass 2 ********************************* 
   ************************************************************************
 */

  CheckMemory;
  parsedatabase(databasefilename);


  if(setGRID==0){
    mesh *meshi;

    xp=NULL; yp=NULL; zp=NULL;
    xp2=NULL; yp2=NULL; zp2=NULL;
    if(NewMemory((void **)&xp,sizeof(float)*(ibartemp+1))==0||
       NewMemory((void **)&yp,sizeof(float)*(jbartemp+1))==0||
       NewMemory((void **)&zp,sizeof(float)*(kbartemp+1))==0||
       NewMemory((void **)&xp2,sizeof(float)*(ibartemp+1))==0||
       NewMemory((void **)&yp2,sizeof(float)*(jbartemp+1))==0||
       NewMemory((void **)&zp2,sizeof(float)*(kbartemp+1))==0
       )return 2;
    for(nn=0;nn<=ibartemp;nn++){
      xp[nn]=xbar0+(float)nn*(xbar-xbar0)/(float)ibartemp;
    }
    for(nn=0;nn<=jbartemp;nn++){
      yp[nn]=ybar0+(float)nn*(ybar-ybar0)/(float)jbartemp;
    }
    for(nn=0;nn<=kbartemp;nn++){
      zp[nn]=zbar0+(float)nn*(zbar-zbar0)/(float)kbartemp;
    }
    meshi=meshinfo;
    initmesh(meshi);
    meshi->xplt=xp;
    meshi->yplt=yp;
    meshi->zplt=zp;
    meshi->xplt_orig=xp2;
    meshi->yplt_orig=yp2;
    meshi->zplt_orig=zp2;
    meshi->ibar=ibartemp;
    meshi->jbar=jbartemp;
    meshi->kbar=kbartemp;
    meshi->plotx=ibartemp/2;
    meshi->ploty=jbartemp/2;
    meshi->plotz=kbartemp/2;
  }
  if(setPDIM==0&&roomdefined==1){
    mesh *meshi;

    meshi=meshinfo;
    meshi->xbar0=xbar0;
    meshi->xbar =xbar;
    meshi->xcen=(xbar+xbar0)/2.0;
    meshi->ybar0=ybar0;
    meshi->ybar =ybar;
    meshi->ycen=(ybar+ybar0)/2.0;
    meshi->zbar0=zbar0;
    meshi->zbar =zbar;
    meshi->zcen=(zbar+zbar0)/2.0;
  }
  

  // define texture data structures by constructing a list of unique file names from surfaceinfo and devices   

  update_device_textures();
  if(nsurfaces>0||ndevice_texture_list>0){
    if(NewMemory((void **)&textureinfo,(nsurfaces+ndevice_texture_list)*sizeof(texture))==0)return 2;
  }

  // get texture filename from SURF and device info

  ntextures = 0;
  for(i=0;i<nsurfaces;i++){
    surfi = surfaceinfo + i;
    if(surfi->texturefile==NULL)continue;
    texti = textureinfo + ntextures;
    len = strlen(surfi->texturefile);
    NewMemory((void **)&texti->file,(len+1)*sizeof(char));
    strcpy(texti->file,surfi->texturefile);
    texti->loaded=0;
    texti->used=0;
    texti->display=0;
    ntextures++;
    surfi->textureinfo=textureinfo+ntextures-1;
  }

  for(i=0;i<ndevice_texture_list;i++){
    char *texturefile;

    texturefile = device_texture_list[i];
    texti = textureinfo + ntextures;
    len = strlen(texturefile);
    NewMemory((void **)&texti->file,(len+1)*sizeof(char));
    device_texture_list_index[i]=ntextures;
    strcpy(texti->file,texturefile);
    texti->loaded=0;
    texti->used=0;
    texti->display=0;
    ntextures++;
  }

  // check to see if texture files exist .
  // If so, then convert to OpenGL format 
#ifdef pp_RENDER
  for(i=0;i<ntextures;i++){
    unsigned char *floortex;
    int texwid, texht;

    texti = textureinfo + i;
    texti->loaded=0;
    if(texti->file!=NULL){
      printf("      Loading textures: ");
    }
    else{
      continue;
    }
    dup_texture=0;
    for(j=0;j<i;j++){
      textj = textureinfo + j;
      if(textj->loaded==0)continue;
      if(strcmp(texti->file,textj->file)==0){
        texti->name=textj->name;
        texti->loaded=1;
        dup_texture=1;
      }
    }
    if(dup_texture==1){
      printf("%s - duplicate\n",texti->file);
      continue;
    }
    CheckMemory;
    glGenTextures(1,&texti->name);
    glBindTexture(GL_TEXTURE_2D,texti->name);
    floortex=readpicture(texti->file,&texwid,&texht);
    if(floortex==NULL){
      printf(" - failed\n");
      continue;
    }
    errorcode=gluBuild2DMipmaps(GL_TEXTURE_2D,4, texwid, texht, GL_RGBA, GL_UNSIGNED_BYTE, floortex);
    if(errorcode!=0){
      FREEMEMORY(floortex);
      printf(" - failed\n");
      continue;
    }
    FREEMEMORY(floortex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    texti->loaded=1;
    printf(" - completed\n");
  }
  CheckMemory;
#endif
  if(ntextures==0)FREEMEMORY(textureinfo);

  // define colobar textures

  printf("      Loading colorbar texture: ");

  glGenTextures(1,&texture_colorbar_id);
  glBindTexture(GL_TEXTURE_1D,texture_colorbar_id);
  glTexImage1D(GL_TEXTURE_1D,0,4,256,0,GL_RGBA,GL_FLOAT,rgb_full);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);

  glGenTextures(1,&texture_slice_colorbar_id);
  glBindTexture(GL_TEXTURE_1D,texture_slice_colorbar_id);
  glTexImage1D(GL_TEXTURE_1D,0,4,256,0,GL_RGBA,GL_FLOAT,rgb_slice);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);

  glGenTextures(1,&texture_plot3d_colorbar_id);
  glBindTexture(GL_TEXTURE_1D,texture_plot3d_colorbar_id);
  glTexImage1D(GL_TEXTURE_1D,0,4,256,0,GL_RGBA,GL_FLOAT,rgb_plot3d);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);

  glGenTextures(1,&texture_iso_colorbar_id);
  glBindTexture(GL_TEXTURE_1D,texture_iso_colorbar_id);
  glTexImage1D(GL_TEXTURE_1D,0,4,256,0,GL_RGBA,GL_FLOAT,rgb_iso);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);

  CheckMemory;

  printf(" - completed\n");
#ifdef pp_GPU
  createDepthTexture();
#endif

  if(autoterrain==1){
    texture *tt;
    unsigned char *floortex;
    int texwid, texht;

#ifdef pp_RENDER
    printf("      Loading terrain texture: ");
#endif
    tt = terrain_texture;
    tt->loaded=0;
    tt->used=0;
    tt->display=0;

#ifdef pp_RENDER
    glGenTextures(1,&tt->name);
    glBindTexture(GL_TEXTURE_2D,tt->name);
    floortex=NULL;
    errorcode=1;
    if(tt->file!=NULL){
      floortex=readpicture(tt->file,&texwid,&texht);
    }
    if(floortex!=NULL){
      errorcode=gluBuild2DMipmaps(GL_TEXTURE_2D,4, texwid, texht, GL_RGBA, GL_UNSIGNED_BYTE, floortex);
    }
    if(errorcode!=0){
      printf(" - failed\n");
    }
    FREEMEMORY(floortex);
    if(errorcode==0){
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
      tt->loaded=1;
      printf(" - completed\n");
    }
#endif
  }

/* 
    Initialize blockage labels and blockage surface labels

    Define default surface for each block.
    Define default vent surface for each block.
  
  */

  surfacedefault=&sdefault;
  for(k=0;k<nsurfaces;k++){
    if(strcmp(surfacedefaultlabel,surfaceinfo[k].surfacelabel)==0){
      surfacedefault=surfaceinfo+k;
      break;
    }
  }
  vent_surfacedefault=&v_surfacedefault;
  for(k=0;k<nsurfaces;k++){
    if(strcmp(vent_surfacedefault->surfacelabel,surfaceinfo[k].surfacelabel)==0){
      vent_surfacedefault=surfaceinfo+k;
      break;
    }
  }

  exterior_surfacedefault=&e_surfacedefault;
  for(k=0;k<nsurfaces;k++){
    if(strcmp(exterior_surfacedefault->surfacelabel,surfaceinfo[k].surfacelabel)==0){
      exterior_surfacedefault=surfaceinfo+k;
      break;
    }
  }

  nbtemp=0; nvents=0;
  itrnx=0, itrny=0, itrnz=0, igrid=0, ipdim=0, iobst=0, ivent=0;
  ioffset=0;
  npartclassinfo=0;
  if(noffset==0)ioffset=1;

/* 
   ************************************************************************
   ************************ start of pass 3 ********************************* 
   ************************************************************************
 */

  rewind(stream);
  printf("   pass 2 completed\n");
  printf("   pass 3\n");
  startpass=1;

  while(!feof(stream)){

    if(nvent==0){
      nvent=0;
      strcpy(buffer,"VENT");
    }
    else{
      if(startpass==1&&noGRIDpresent==1){
        strcpy(buffer,"GRID");
        startpass=0;
      }
      else{
        if(fgets(buffer,255,stream)==NULL)break;
        if(strncmp(buffer," ",1)==0||buffer[0]==0)continue;
      }
    }
    CheckMemory;

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++ CLASS_OF_PARTICLES +++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"CLASS_OF_PARTICLES",18) == 1||
       match(buffer,"CLASS_OF_HUMANS",15) == 1){
      float rgb_class[4];
      part5class *partclassi;
      char *device_ptr;
      char buffer_copy[1024];
      char *smokeview_id, *prop_id;
      int nvar;

      partclassi = partclassinfo + npartclassinfo;
      partclassi->kind=PARTICLES;
      if(match(buffer,"CLASS_OF_HUMANS",15) == 1)partclassi->kind=HUMANS;
      fgets(buffer,255,stream);

      get_labels(buffer,&device_ptr,&prop_id);
      partclassi->prop=get_prop_id(prop_id);
      update_partclass_depend(partclassi);

      npartclassinfo++;
      continue;
    }
    
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ HRRPUVCUT ++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"HRRPUVCUT",9) == 1){
      int nhrrpuvcut;
      float hrrpuvcut;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&nhrrpuvcut);
      if(nhrrpuvcut>=1){
        fgets(buffer,255,stream);
        sscanf(buffer,"%f",&hrrpuvcut);
        for(i=0;i<nmeshes;i++){
          mesh *meshi;

          meshi=meshinfo+i;
          meshi->hrrpuv_cutoff=hrrpuvcut;
        }
        for(i=1;i<nhrrpuvcut;i++){
          mesh *meshi;

          fgets(buffer,255,stream);
          if(i>=nmeshes)continue;
          sscanf(buffer,"%f",&hrrpuvcut);
          meshi=meshinfo+i;
          meshi->hrrpuv_cutoff=hrrpuvcut;
        }
      }
      continue;
    }

    /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ PL3D ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"PL3D",4) == 1){
      nn_plot3d++;
      trim(buffer);
      len=strlen(buffer);
      blocknumber = 0;
      if(nmeshes>1){
        blocknumber=ioffset-1;
      }
      else{
        blocknumber=0;
      }
      if(strlen(buffer)>5){
        sscanf(buffer,"%s %f %i",buffer2,&time,&blocktemp);
        if(blocktemp>0&&blocktemp<=nmeshes)blocknumber = blocktemp-1;
      }
      else{
        time=-1.0;
      }
      if(fgets(buffer,255,stream)==NULL){
        nplot3d_files--;
        break;
      }
      bufptr=trim_string(buffer);
      len=strlen(bufptr);

      p=plot3dinfo+iplot3d;
      p->blocknumber=blocknumber;
      p->seq_id=nn_plot3d;
      p->autoload=0;
      p->time=time;
      p->loaded=0;
      p->display=0;

      NewMemory((void **)&p->reg_file,(unsigned int)(len+1));
      STRCPY(p->reg_file,bufptr);

      NewMemory((void **)&p->comp_file,(unsigned int)(len+4+1));
      STRCPY(p->comp_file,bufptr);
      STRCAT(p->comp_file,".svz");

      if(STAT(p->comp_file,&statbuffer)==0){
        p->compression_type=1;
        p->file=p->comp_file;
      }
      else{
        p->compression_type=0;
        p->file=p->reg_file;
      }
      //disable compression for now
      p->compression_type=0;
      p->file=p->reg_file;

      if(STAT(p->file,&statbuffer)!=0){
        for(n=0;n<5;n++){
          if(readlabels(&p->label[n],stream)==2)return 2;
        }
        nplot3d_files--;
      }
      else{
        p->u = -1;
        p->v = -1;
        p->w = -1;
        for(n=0;n<5;n++){
          if(readlabels(&p->label[n],stream)==2)return 2;
          if(match(p->label[n].shortlabel,UVEL,5) == 1){
            p->u = n;
          }
          if(match(p->label[n].shortlabel,VVEL,5) == 1){
            p->v = n;
          }
          if(match(p->label[n].shortlabel,WVEL,5) == 1){
            p->w = n;
          }
        }
        if(p->u>-1||p->v>-1||p->w>-1){
          p->nvars=mxplot3dvars;
        }
        else{
          p->nvars=5;
        }
        if(NewMemory((void **)&p->label[5].longlabel,6)==0)return 2;
        if(NewMemory((void **)&p->label[5].shortlabel,6)==0)return 2;
        if(NewMemory((void **)&p->label[5].unit,4)==0)return 2;

        STRCPY(p->label[5].longlabel,"Speed");
        STRCPY(p->label[5].shortlabel,"Speed");
        STRCPY(p->label[5].unit,"m/s");

        STRCPY(p->longlabel,"");
        for(n=0;n<5;n++){
          STRCAT(p->longlabel,p->label[n].shortlabel);
          if(n!=4)STRCAT(p->longlabel,", ");
        }

        iplot3d++;
      }
      continue;
    }

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ OFFSET ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"OFFSET",6) == 1){
      mesh *meshi;

      ioffset++;
      meshi=meshinfo+ioffset-1;
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",meshi->offset,meshi->offset+1,meshi->offset+2);
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ VENTGEOM ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"VENTGEOM",8)==1||match(buffer,"VFLOWGEOM",9)==1){
      int zonevent_orien=0;
      int hvent_type=0;
      zvent *zvi;
      float vent_area;
      int roomfrom, roomto, face;
      roomdata *roomi;
      float color[4];

      nzvents++;
      zvi = zventinfo + nzvents - 1;
      if(match(buffer,"VFLOWGEOM",9)==1)zonevent_orien=1;
      zvi->vent_orien=zonevent_orien;
      if(fgets(buffer,255,stream)==NULL)break;
      color[0]=1.0;
      color[1]=0.0;
      color[2]=1.0;
      color[3]=1.0;
      if(zonevent_orien==0){
        sscanf(buffer,"%i %i %i %f %f %f %f %f %f %f",
          &roomfrom,&roomto, &face,&width,&ventoffset,&bottom,&top,
          color,color+1,color+2
          );
        
        if(roomfrom<1||roomfrom>nrooms)roomfrom=nrooms+1;
        roomi = roominfo + roomfrom-1;
        zvi->room1 = roomi;

        if(roomto<1||roomto>nrooms)roomto=nrooms+1;
        zvi->room2=roominfo+roomto-1;

        zvi->dir=face;
        zvi->z1=roomi->z0+bottom;
        zvi->z2=roomi->z0+top;
        zvi->face=face;
        switch (face){
        case 1:
          zvi->yy=roomi->y0;
          zvi->x1=roomi->x0+ventoffset;
          zvi->x2=roomi->x0+ventoffset+width;
          break;
        case 2:
          zvi->yy=roomi->x1;
          zvi->x1=roomi->y0+ventoffset;
          zvi->x2=roomi->y0+ventoffset+width;
          break;
        case 3:
          zvi->yy=roomi->y1;
          zvi->x1=roomi->x0+ventoffset;
          zvi->x2=roomi->x0+ventoffset+width;
          break;
        case 4:
          zvi->yy=roomi->x0;
          zvi->x1=roomi->y0+ventoffset;
          zvi->x2=roomi->y0+ventoffset+width;
          break;
        default:
          ASSERT(FFALSE);
        }
      }
      else{
        float ventside;
        float xcen, ycen;
        int r_from, r_to;

        sscanf(buffer,"%i %i %i %f %i %f %f %f",
          &r_from,&r_to,&face,&vent_area,&hvent_type,
          color,color+1,color+2
          );
        zvi->dir=face;
        roomfrom=r_from;
        roomto=r_to;
        if(roomfrom<1||roomfrom>nrooms){
          roomfrom=r_to;
          roomto=r_from;
          if(roomfrom<1||roomfrom>nrooms){
            roomfrom=1;
          }
        }
        roomi = roominfo + roomfrom - 1;
        if(vent_area<0.0)vent_area=-vent_area;
        ventside=sqrt(vent_area);
        xcen = (roomi->x0+roomi->x1)/2.0;
        ycen = (roomi->y0+roomi->y1)/2.0;
        if(face==5){
          zvi->zz=roomi->z0;
        }
        else{
          zvi->zz=roomi->z1;
        }
        switch (hvent_type){
        case 1:
        case 2:
          zvi->x1=xcen-ventside/2.0;
          zvi->x2=xcen+ventside/2.0;
          zvi->y1=ycen-ventside/2.0;
          zvi->y2=ycen+ventside/2.0;
          break;
        default:
          ASSERT(FFALSE);
          break;
        }

        zvi->vent_type=hvent_type;
        zvi->area=vent_area;
        zvi->face=face;
      }
      zvi->color=getcolorptr(color);
      CheckMemory;
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ TITLE1 ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"TITLE1",6)==1){
      if(fgets(buffer,255,stream)==NULL)break;
      bufptr=trim_string(buffer);
      strcpy(TITLE1,bufptr);
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ TITLE2 ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"TITLE2",6)==1){
      if(fgets(buffer,255,stream)==NULL)break;
      bufptr=trim_string(buffer);
      strcpy(TITLE2,bufptr);
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ FIRE ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"FIRE",4)==1){
      ifire++;
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%i %f %f %f",&roomnumber,&fireinfo[ifire-1].x,
        &fireinfo[ifire-1].y,&fireinfo[ifire-1].z);
      if(roomnumber>=1&&roomnumber<=nrooms){
        fireinfo[ifire-1].valid=1;
        fireinfo[ifire-1].roomnumber=roomnumber;
        fireinfo[ifire-1].absx=roominfo[roomnumber-1].x0+fireinfo[ifire-1].x;
        fireinfo[ifire-1].absy=roominfo[roomnumber-1].y0+fireinfo[ifire-1].y;
        fireinfo[ifire-1].absz=roominfo[roomnumber-1].z0+fireinfo[ifire-1].z;
      }
      else{
        fireinfo[ifire-1].valid=0;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ PDIM ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"PDIM",4) == 1){
      mesh *meshi;
      float *meshrgb;

      ipdim++;
      meshi=meshinfo+ipdim-1;
      meshrgb = meshi->meshrgb;

      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f %f %f %f %f %f %f",&xbar0,&xbar,&ybar0,&ybar,&zbar0,&zbar,meshrgb,meshrgb+1,meshrgb+2);

      if(meshrgb[0]!=0.0||meshrgb[1]!=0.0||meshrgb[2]!=0.0){
        meshi->meshrgb_ptr=meshi->meshrgb;
      }
      else{
        meshi->meshrgb_ptr=NULL;
      }

      meshi->xbar0=xbar0;
      meshi->xbar =xbar;
      meshi->xcen=(xbar+xbar0)/2.0;
      meshi->ybar0=ybar0;
      meshi->ybar =ybar;
      meshi->ycen =(ybar+ybar0)/2.0;
      meshi->zbar0=zbar0;
      meshi->zbar =zbar;
      meshi->zcen =(zbar+zbar0)/2.0;
      if(ntrnx==0){
        for(nn=0;nn<=meshi->ibar;nn++){
          meshi->xplt[nn]=xbar0+(float)nn*(xbar-xbar0)/(float)meshi->ibar;
        }
      }
      if(ntrny==0){
        for(nn=0;nn<=meshi->jbar;nn++){
          meshi->yplt[nn]=ybar0+(float)nn*(ybar-ybar0)/(float)meshi->jbar;
        }
      }
      if(ntrnz==0){
        for(nn=0;nn<=meshi->kbar;nn++){
          meshi->zplt[nn]=zbar0+(float)nn*(zbar-zbar0)/(float)meshi->kbar;
        }
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ TRNX ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"TRNX",4)==1){
      itrnx++;
      xpltcopy=meshinfo[itrnx-1].xplt;
      ibartemp=meshinfo[itrnx-1].ibar;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&idummy);
      for(nn=0;nn<idummy;nn++){
        fgets(buffer,255,stream);
      }
      for(nn=0;nn<=ibartemp;nn++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%i %f",&idummy,xpltcopy);xpltcopy++;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ TRNY ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"TRNY",4) == 1){
      itrny++;
      ypltcopy=meshinfo[itrny-1].yplt;
      jbartemp=meshinfo[itrny-1].jbar;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&idummy);
      for(nn=0;nn<idummy;nn++){
        fgets(buffer,255,stream);
      }
      for(nn=0;nn<=jbartemp;nn++){
//        if(jbartemp==2&&nn==2)continue;
        fgets(buffer,255,stream);
        sscanf(buffer,"%i %f",&idummy,ypltcopy);
        ypltcopy++;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ TRNZ ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"TRNZ",4) == 1){
      itrnz++;
      zpltcopy=meshinfo[itrnz-1].zplt;
      kbartemp=meshinfo[itrnz-1].kbar;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&idummy);
      for(nn=0;nn<idummy;nn++){
        fgets(buffer,255,stream);
      }
      for(nn=0;nn<=kbartemp;nn++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%i %f",&idummy,zpltcopy);zpltcopy++;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ TREE ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    /*
typedef struct {
  float xyz[3];
  float trunk_diam;
  float tree_height;
  float base_diam;
  float base_height;
*/
    if(match(buffer,"TREE",4) == 1&&match(buffer,"TREESTATE",9) !=1){
      fgets(buffer,255,stream);
      if(ntreeinfo!=0)continue;
      sscanf(buffer,"%i",&ntreeinfo);
      if(ntreeinfo>0){
        NewMemory((void **)&treeinfo,sizeof(treedata)*ntreeinfo);
        for(i=0;i<ntreeinfo;i++){
          treedata *treei;
          float *xyz;

          treei = treeinfo + i;
          treei->time_char=-1.0;
          treei->time_complete=-1.0;
          xyz = treei->xyz;
          fgets(buffer,255,stream);
          sscanf(buffer,"%f %f %f %f %f %f %f",
            xyz,xyz+1,xyz+2,
            &treei->trunk_diam,&treei->tree_height,
            &treei->base_diam,&treei->base_height);
        }
      }
      continue;
    }
    if(match(buffer,"TREESTATE",9) ==1){
      int tree_index, tree_state;
      float tree_time;
      treedata *treei;

      if(ntreeinfo==0)continue;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %i %f",&tree_index,&tree_state,&tree_time);
      if(tree_index>=1&&tree_index<=ntreeinfo){
        treei = treeinfo + tree_index - 1;
        if(tree_state==1){
          treei->time_char = tree_time;
        }
        else if(tree_state==2){
          treei->time_complete = tree_time;
        }
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ OBST ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

    // note - OBST lines for autoterran==1 are processed in pass 4
    //        because Smokeview needs info from pass 3

    if(match(buffer,"OBST",4) == 1&&autoterrain==0){
      mesh *meshi;
      propdata *prop;
      char *proplabel;

      CheckMemoryOff;
      iobst++;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&tempval);

      meshi=meshinfo+iobst-1;
      if(tempval<=0)tempval=0;
      meshi->nbptrs=tempval;
      if(tempval==0)continue;
      nbtemp=tempval;
      meshi->blockageinfoptrs=NULL;
      if(
        NewMemory((void **)&meshi->blockageinfoptrs,sizeof(blockagedata *)*nbtemp)==0||
        NewMemory((void **)&meshi->deletelist,sizeof(blockagedata *)*nbtemp)==0
        )return 2;


      ntotal_blockages+=nbtemp;
      for(nn=0;nn<nbtemp;nn++){
        int s_num[6];

        meshi->blockageinfoptrs[nn]=NULL;
        meshi->deletelist[nn]=NULL;
        NewMemory((void **)&meshi->blockageinfoptrs[nn],sizeof(blockagedata));
        bc=meshi->blockageinfoptrs[nn];
        initobst(bc,surfacedefault,nn+1,iobst-1);
        fgets(buffer,255,stream);
        trim(buffer);
        for(i=0;i<6;i++){
          s_num[i]=-1;
        }
        proplabel=strchr(buffer,'%');
        prop=NULL;
        if(proplabel!=NULL){
          proplabel++;
          trim(proplabel);
          proplabel = trim_front(proplabel);
          for(i=0;i<npropinfo;i++){
            propdata *propi;

            propi = propinfo + i;
            if(STRCMP(proplabel,propi->label)==0){
              prop = propi;
              propi->inblockage=1;
              break;
            }
          }
        }
        bc->prop=prop;
        {
          float t_origin[3];
          float *xyzEXACT;

          t_origin[0]=texture_origin[0];
          t_origin[1]=texture_origin[1];
          t_origin[2]=texture_origin[2];
          xyzEXACT = bc->xyzEXACT;
          sscanf(buffer,"%f %f %f %f %f %f %i %i %i %i %i %i %i %f %f %f",
            xyzEXACT,xyzEXACT+1,xyzEXACT+2,xyzEXACT+3,xyzEXACT+4,xyzEXACT+5,
            &(bc->id),s_num+DOWN_X,s_num+UP_X,s_num+DOWN_Y,s_num+UP_Y,s_num+DOWN_Z,s_num+UP_Z,
            t_origin,t_origin+1,t_origin+2);
          bc->xmin=xyzEXACT[0];
          bc->xmax=xyzEXACT[1];
          bc->ymin=xyzEXACT[2];
          bc->ymax=xyzEXACT[3];
          bc->zmin=xyzEXACT[4];
          bc->zmax=xyzEXACT[5];
          bc->texture_origin[0]=t_origin[0];
          bc->texture_origin[1]=t_origin[1];
          bc->texture_origin[2]=t_origin[2];
          if(bc->id<0){
            bc->changed=1;
            bc->id=-bc->id;
          }
        }

        /* define block label */

        sprintf(buffer,"**blockage %i",bc->id);
        len=strlen(buffer);
        ResizeMemory((void **)&bc->label,(len+1)*sizeof(char));
        strcpy(bc->label,buffer);

        for(i=0;i<6;i++){
          if(surfaceinfo==NULL||s_num[i]<0||s_num[i]>=nsurfaces)continue;
          surfi=surfaceinfo+s_num[i];
          bc->surf[i]=surfi;
        }
        for(i=0;i<6;i++){
          bc->surf[i]->used_by_obst=1;
        }
        setsurfaceindex(bc);
      }

#define COLOR_INVISIBLE -2

#define BLOCK_OUTLINE 2

      for(nn=0;nn<nbtemp;nn++){
        bc=meshi->blockageinfoptrs[nn];
        colorindex=-1;
        blocktype=-1;

        /* 
        OBST format:
        i1 i2 j1 j2 k1 k2 colorindex blocktype r g b : ignore rgb if blocktype != -3
        
        int colorindex, blocktype;
        colorindex: -1 default color
                    -2 invisible
                    -3 use r g b color
                    >=0 color/color2/texture index
        blocktype: 0 regular block
                   2 outline
                   3 smoothed block
                   (note: if blocktype&8 == 1 then this is a "terrain" blockage
                         if so then subtract 8 and set bc->is_wuiblock=1)
        r g b           colors used if colorindex==-3
        */

        fgets(buffer,255,stream);
        ijk = bc->ijk;
        sscanf(buffer,"%i %i %i %i %i %i %i %i",
          ijk,ijk+1,ijk+2,ijk+3,ijk+4,ijk+5,
          &colorindex,&blocktype);
        if(blocktype>0){
          if((blocktype&8)==1){
            bc->is_wuiblock=1;
            blocktype -= 8;
          }
        }

        /* custom color */
        
        if(colorindex==0||colorindex==7)colorindex=-3;

        if(colorindex==-3){
          ijk = bc->ijk;

          s_color[0]=-1.0;
          s_color[1]=-1.0;
          s_color[2]=-1.0;
          s_color[3]=1.0;

          sscanf(buffer,"%i %i %i %i %i %i %i %i %f %f %f %f",
            ijk,ijk+1,ijk+2,ijk+3,ijk+4,ijk+5,
            &colorindex,&blocktype,s_color,s_color+1,s_color+2,s_color+3);
          if(blocktype>0){
            if((blocktype&8)==1){
              bc->is_wuiblock=1;
              blocktype -= 8;
            }
          }

          if(s_color[3]<0.999){
            bc->transparent=1;
          }
          if(colorindex==0||colorindex==7){
            switch (colorindex){
            case 0:
              s_color[0]=1.0;
              s_color[1]=1.0;
              s_color[2]=1.0;
              break;
            case 7:
              s_color[0]=0.0;
              s_color[1]=0.0;
              s_color[2]=0.0;
              break;
            default:
              ASSERT(FFALSE);
              break;
            }
            colorindex=-3;
          }
          if(s_color[0]>=0.0&&s_color[1]>=0.0&&s_color[2]>=0.0){
            bc->color=getcolorptr(s_color);
          }
          bc->nnodes=(ijk[1]+1-ijk[0])*(ijk[3]+1-ijk[2])*(ijk[5]+1-ijk[4]);
          bc->useblockcolor=1;
        }
        else{
          if(colorindex>=0){
            bc->color = getcolorptr(rgb[nrgb+colorindex]);
            bc->useblockcolor=1;
            bc->usecolorindex=1;
            bc->colorindex=colorindex;
            updateindexcolors=1;
          }
        }

        bc->colorindex=colorindex;
        bc->type=blocktype;

        if(colorindex==COLOR_INVISIBLE){
          bc->type=BLOCK_hidden;
//          bc->del=1;
          bc->invisible=1;
        }
        if(bc->useblockcolor==0){
          bc->color=bc->surf[0]->color;
        }
        else{
          if(colorindex==-1){
            updateindexcolors=1;
          }
        }
        if(bc->type==BLOCK_smooth){
          ntotal_smooth_blockages++;
          use_menusmooth=1;
        }

      }
      CheckMemoryOn;
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ VENT ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
        /* 
        new VENT format:
        i1 i2 j1 j2 k1 k2 ventindex venttype r g b a
        
        int ventindex, venttype;
        vent index: -99 or 99 : use default color
                     +n or -n : use n'th pallette color
                    < 0       : DO NOT draw boundary file over this vent
                    > 0       : DO draw boundary file over this vent
        vent type: 0 solid surface
                   2 outline
                  -2 hidden 

        r g b           colors - only specify if you wish to over-ride surface or default
        */

    if(match(buffer,"VENT",4) == 1&&match(buffer,"VENTGEOM",8) != 1){
      mesh *meshi;

      ivent++;
      meshi=meshinfo+ivent-1;
      xplttemp=meshi->xplt;
      yplttemp=meshi->yplt;
      zplttemp=meshi->zplt;
      if(nvent==0){
        strcpy(buffer,"0 0");
        nvent=1;
      }
      else{
        fgets(buffer,255,stream);
      }
      ndummyvents=0;
      sscanf(buffer,"%i %i",&nvents,&ndummyvents);
      if(ndummyvents!=0)dummyvents=1;
      meshi->nvents=nvents;
      meshi->ndummyvents=ndummyvents;
      vinfo=NULL;
      meshi->ventinfo=vinfo;
      if(NewMemory((void **)&vinfo,(nvents+12)*sizeof(ventdata))==0)return 2;
      meshi->ventinfo=vinfo;

      for(nn=0;nn<nvents+12;nn++){
        int s_num[6];

        vi=vinfo+nn;
        vi->transparent=0;
        vi->useventcolor=0;
        vi->usecolorindex=0;
        vi->nshowtime=0;
        vi->isOpenvent=0;
		    vi->hideboundary=0;
        vi->surf[0]=vent_surfacedefault;
        vi->textureinfo[0]=NULL;
        vi->texture_origin[0]=texture_origin[0];
        vi->texture_origin[1]=texture_origin[1];
        vi->texture_origin[2]=texture_origin[2];
        vi->colorindex=-1;
        if(nn>nvents-ndummyvents-1&&nn<nvents){
          vi->dummy=1;
        }
        else{
          vi->dummy=0;
        }
        s_num[0]=-1;
        if(nn<nvents){
          {
            float t_origin[3];
            t_origin[0]=texture_origin[0];
            t_origin[1]=texture_origin[1];
            t_origin[2]=texture_origin[2];
            fgets(buffer,255,stream);
            sscanf(buffer,"%f %f %f %f %f %f %i %i %f %f %f",
                   &vi->xmin,&vi->xmax,&vi->ymin,&vi->ymax,&vi->zmin,&vi->zmax,
                   &vi->id,s_num,t_origin,t_origin+1,t_origin+2);
            vi->texture_origin[0]=t_origin[0];
            vi->texture_origin[1]=t_origin[1];
            vi->texture_origin[2]=t_origin[2];
          }
        }
        else{
          vi->xmin=meshi->xbar0+meshi->offset[0];
          vi->xmax=meshi->xbar +meshi->offset[0];
          vi->ymin=meshi->ybar0+meshi->offset[1];
          vi->ymax=meshi->ybar +meshi->offset[1];
          vi->zmin=meshi->zbar0+meshi->offset[2];
          vi->zmax=meshi->zbar +meshi->offset[2];
          s_num[0]=-1;
          switch (nn-nvents){
          case DOWN_Y:
          case DOWN_Y+6:
            vi->ymax=vi->ymin;
            break;
          case UP_X:
          case UP_X+6:
            vi->xmin=vi->xmax;
            break;
          case UP_Y:
          case UP_Y+6:
            vi->ymin=vi->ymax;
            break;
          case DOWN_X:
          case DOWN_X+6:
            vi->xmax=vi->xmin;
            break;
          case DOWN_Z:
          case DOWN_Z+6:
            vi->zmax=vi->zmin;
            break;
          case UP_Z:
          case UP_Z+6:
            vi->zmin=vi->zmax;
            break;
          default:
            ASSERT(FFALSE);
            break;
          }
          if(nn>=nvents+6){
            vi->surf[0]=exterior_surfacedefault;
          }
        }
        if(surfaceinfo!=NULL&&s_num[0]>=0&&s_num[0]<nsurfaces){
          vi->surf[0]=surfaceinfo+s_num[0];
          if(vi->surf[0]!=NULL&&strncmp(vi->surf[0]->surfacelabel,"OPEN",4)==0){
            vi->isOpenvent=1;
          }
          vi->surf[0]->used_by_vent=1;
        }
      }
      for(nn=0;nn<nvents+12;nn++){
        vi = vinfo+nn;
        vi->type=vi->surf[0]->type;
        vi->color=vi->surf[0]->color;
        s_color[0]=vi->surf[0]->color[0];
        s_color[1]=vi->surf[0]->color[1];
        s_color[2]=vi->surf[0]->color[2];
        s_color[3]=vi->surf[0]->color[3];
//        s_color[3]=1.0;
        venttype=-99;
        if(nn<nvents){
          float s2_color[4];

          s2_color[0]=-1.0;
          s2_color[1]=-1.0;
          s2_color[2]=-1.0;
          s2_color[3]=1.0;


          fgets(buffer,255,stream);
          sscanf(buffer,"%i %i %i %i %i %i %i %i %f %f %f %f",
               &iv1,&iv2,&jv1,&jv2,&kv1,&kv2,
               &ventindex,&venttype,
               s2_color,s2_color+1,s2_color+2,s2_color+3);
          if(s2_color[0]>=0.0&&s2_color[1]>=0.0&&s2_color[2]>=0.0){
            s_color[0]=s2_color[0];
            s_color[1]=s2_color[1];
            s_color[2]=s2_color[2];
            if(s2_color[3]<0.99){
              vi->transparent=1;
              s_color[3]=s2_color[3];
            }
            vi->useventcolor=1;
          }
		      if(ventindex<0)vi->hideboundary=1;
          if(venttype!=-99){
            vi->type=venttype;
          }
          //vi->colorindex=ventindex;
          if(ventindex!=-99&&ventindex!=99){
            if(ventindex<0)ventindex=-ventindex;
            if(ventindex>nrgb2-1)ventindex=nrgb2-1;
            s_color[0]=rgb[nrgb+ventindex][0];
            s_color[1]=rgb[nrgb+ventindex][1];
            s_color[2]=rgb[nrgb+ventindex][2];
            s_color[3]=1.0;
            vi->colorindex=ventindex;
            vi->usecolorindex=1;
            vi->useventcolor=1;
            updateindexcolors=1;
          }
          else{
          }
          vi->color = getcolorptr(s_color);
        }
        else{
          iv1=0; iv2 = meshi->ibar;
          jv1=0; jv2 = meshi->jbar;
          kv1=0; kv2 = meshi->kbar;
          ventindex=-99;
          vi->dir=nn-nvents;
          if(vi->dir>5)vi->dir-=6;
          vi->dir2=0;
          switch (nn-nvents){
          case DOWN_Y:
          case DOWN_Y+6:
            jv2=jv1;
            if(nn>=nvents+6)vi->dir=UP_Y;
            vi->dir2=2;
            break;
          case UP_X:
          case UP_X+6:
            iv1=iv2;
            if(nn>=nvents+6)vi->dir=DOWN_X;
            vi->dir2=1;
            break;
          case UP_Y:
          case UP_Y+6:
            jv1=jv2;
            if(nn>=nvents+6)vi->dir=DOWN_Y;
            vi->dir2=2;
            break;
          case DOWN_X:
          case DOWN_X+6:
            iv2=iv1;
            if(nn>=nvents+6)vi->dir=UP_X;
            vi->dir2=1;
            break;
          case DOWN_Z:
          case DOWN_Z+6:
            kv2=kv1;
            if(nn>=nvents+6)vi->dir=UP_Z;
            vi->dir2=3;
            break;
          case UP_Z:
          case UP_Z+6:
            kv1=kv2;
            if(nn>=nvents+6)vi->dir=DOWN_Z;
            vi->dir2=3;
            break;
          default:
            ASSERT(FFALSE);
            break;
          }
        }
        if(vi->transparent==1)nvent_transparent++;
        vi->linewidth=&ventlinewidth;
        vi->showhide=NULL;
        vi->showtime=NULL;
        vi->showtimelist=NULL;
        vi->xvent1 = xplttemp[iv1];
        vi->xvent2 = xplttemp[iv2]; 
        vi->yvent1 = yplttemp[jv1]; 
        vi->yvent2 = yplttemp[jv2]; 
        vi->zvent1 = zplttemp[kv1]; 
        vi->zvent2 = zplttemp[kv2];
        vi->imin = iv1;
        vi->imax = iv2;
        vi->jmin = jv1;
        vi->jmax = jv2;
        vi->kmin = kv1;
        vi->kmax = kv2;
        if(nn>=nvents&&nn<nvents+6){
          vi->color=foregroundcolor;
        }
        ASSERT(vi->color!=NULL);
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ PART ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"PART",4) == 1 || match(buffer,"EVAC",4)==1
      ||match(buffer,"PRT5",4)==1||match(buffer,"EVA5",4)==1
      ){
      unsigned int lenkey;

      nn_part++;

      parti = partinfo + ipart;

      parti->version=0;
      if(match(buffer,"PRT5",4)==1||match(buffer,"EVA5",4)==1){
        parti->version=1;
      }
      lenkey=4;
      parti->evac=0;
      if(match(buffer,"EVAC",4)==1
        ||match(buffer,"EVA5",4)==1
        ){
        parti->evac=1;
        nevac++;
      }
      len=strlen(buffer);
      if(nmeshes>1){
        blocknumber=ioffset-1;
      }
      else{
        blocknumber=0;
      }
      if(len>lenkey+1){
        buffer3=buffer+lenkey;
        if(parti->evac==1){
          float zoffset=0.0;

          sscanf(buffer3,"%i %f",&blocknumber,&zoffset);
          parti->zoffset=zoffset;
        }
        else{
          sscanf(buffer3,"%i",&blocknumber);
        }
        blocknumber--;
      }

      parti->blocknumber=blocknumber;
      parti->seq_id=nn_part;
      parti->autoload=0;
      if(fgets(buffer,255,stream)==NULL){
        npart_files--;
        break;
      }

      bufptr=trim_string(buffer);
      len=strlen(bufptr);
      parti->reg_file=NULL;
      if(NewMemory((void **)&parti->reg_file,(unsigned int)(len+1))==0)return 2;
      STRCPY(parti->reg_file,bufptr);

      parti->size_file=NULL;
      if(NewMemory((void **)&parti->size_file,(unsigned int)(len+1+3))==0)return 2;
      STRCPY(parti->size_file,bufptr);
      STRCAT(parti->size_file,".sz");

      parti->comp_file=NULL;
      if(NewMemory((void **)&parti->comp_file,(unsigned int)(len+1+4))==0)return 2;
      STRCPY(parti->comp_file,bufptr);
      STRCAT(parti->comp_file,".svz");

      if(STAT(parti->comp_file,&statbuffer)==0){
        parti->compression_type=1;
        parti->file=parti->comp_file;
      }
      else{
        parti->compression_type=0;
        if(STAT(parti->reg_file,&statbuffer)==0){
          parti->file=parti->reg_file;
        }
        else{
          FREEMEMORY(parti->reg_file);
          FREEMEMORY(parti->comp_file);
          FREEMEMORY(parti->size_file);
          parti->file=NULL;
        }
      }
      parti->compression_type=0;
      parti->sort_tags_loaded=0;
      parti->loaded=0;
      parti->display=0;
      parti->ptimes=NULL; 
      parti->xpart=NULL;  
      parti->ypart=NULL;  
      parti->zpart=NULL;  
      parti->xpartb=NULL; 
      parti->ypartb=NULL; 
      parti->zpartb=NULL; 
      parti->xparts=NULL; 
      parti->yparts=NULL; 
      parti->zparts=NULL; 
      parti->tpart=NULL;  
      parti->itpart=NULL; 
      parti->isprink=NULL;
      parti->sframe=NULL; 
      parti->bframe=NULL;
      parti->sprframe=NULL;
      parti->ptimeslist=NULL;
      parti->particle_type=0;
      parti->droplet_type=0;

      parti->data5=NULL;
      parti->partclassptr=NULL;

      if(parti->version==1){
        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&parti->nclasses);
        if(parti->nclasses>0){
//          createnulllabel(&parti->label);
          if(parti->file!=NULL)NewMemory((void **)&parti->partclassptr,parti->nclasses*sizeof(part5class *));
          for(i=0;i<parti->nclasses;i++){
            int iclass;
            int ic,iii;

            fgets(buffer,255,stream);
            if(parti->file==NULL)continue;
            sscanf(buffer,"%i",&iclass);
            if(iclass<1)iclass=1;
            if(iclass>npartclassinfo)iclass=npartclassinfo;
            ic=0;
            for(iii=0;iii<npartclassinfo;iii++){
              part5class *pci;

              pci = partclassinfo + iii;
              if(parti->evac==1&&pci->kind!=HUMANS)continue;
              if(parti->evac==0&&pci->kind!=PARTICLES)continue;
              if(iclass-1==ic){
                parti->partclassptr[i]=pci;
                break;
              }
              ic++;
            }
          }
        }
        // if no classes were specifed for the prt5 entry then assign it the default class
        if(parti->file!=NULL&&parti->nclasses==0){
          NewMemory((void **)&parti->partclassptr,sizeof(part5class *));
            parti->partclassptr[i]=partclassinfo + parti->nclasses;
        }
        if(parti->file==NULL||STAT(parti->file,&statbuffer)!=0){
          npart_files--;
        }
        else{
          ipart++;
        }
      }
      else{
        if(STAT(buffer,&statbuffer)==0){
          if( readlabels(&parti->label,stream)==2 )return 2;
          ipart++;
        }
        else{
          if(readlabels(&parti->label,stream)==2)return 2;
          npart_files--;
        }
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ TARG ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(
      match(buffer,"TARG",4) == 1||
      match(buffer,"FTARG",5) == 1
      ){
      targinfo[itarg].type=1;
      if(match(buffer,"FTARG",5)==1)targinfo[itarg].type=2;

      if(fgets(buffer,255,stream)==NULL){
        ntarg_files--;
        break;
      }
      len=strlen(buffer);
      buffer[len-1]='\0';
      trim(buffer);
      len=strlen(buffer);
      targinfo[itarg].loaded=0;
      targinfo[itarg].display=0;
      if(NewMemory((void **)&targinfo[itarg].file,(unsigned int)(len+1))==0)return 2;
      STRCPY(targinfo[itarg].file,buffer);
      if(targfilename!=NULL){
        FREEMEMORY(targfilename);
        if(NewMemory((void **)&targfilename,(unsigned int)(len+1))==0)return 2;
        STRCPY(targfilename,buffer);
      }
      if(STAT(buffer,&statbuffer)==0){
        itarg++;
      }
      else{
        ntarg_files--;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SYST ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"SYST",4) == 1){
      if(fgets(buffer,255,stream)==NULL)break;
      len=strlen(buffer);
      buffer[len-1]='\0';
      STRCPY(LESsystem,buffer);
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ ENDIAN ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"ENDIAN",6) == 1){
      if(fgets(buffer,255,stream)==NULL)break;
      len=strlen(buffer);
      buffer[len-1]='\0';
      strncpy(LESendian,buffer,1);
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ ENDF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    /* ENDF superscedes ENDIAN */
    if(match(buffer,"ENDF",4) == 1){
      if(fgets(buffer,255,stream)==NULL)break;
      bufptr=trim_string(buffer);
      len=strlen(bufptr);
      NewMemory((void **)&endianfilename,(unsigned int)(len+1));
      strcpy(endianfilename,bufptr);
      ENDIANfile = fopen(endianfilename,"rb");
      if(ENDIANfile!=NULL){
        endian_native = getendian();
        fseek(ENDIANfile,4,SEEK_SET);
        fread(&endian_data,4,1,ENDIANfile);
        fclose(ENDIANfile);
        endian=endian_native;
        if(endian_data!=1)endian=1-endian_native;
        setendian=1;
      }
      continue;
    }

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ INPF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"INPF",4) == 1){
      if(fgets(buffer,255,stream)==NULL)break;
      bufptr=trim_string(buffer);

      len=strlen(bufptr);
      FREEMEMORY(fds_filein);
      if(NewMemory((void **)&fds_filein,(unsigned int)(len+1))==0)return 2;
      STRCPY(fds_filein,bufptr);
      if(STAT(fds_filein,&statbuffer)!=0){
        FreeMemory(fds_filein);
        fds_filein=NULL;
      }

      if(fds_filein!=NULL){
        {
          int returnval;
          returnval=getnewfilename();
          if(returnval!=0)return returnval;
        }
      }
      if(chidfilebase==NULL){
        char *chidptr=NULL;

        if(fds_filein!=NULL)chidptr=get_chid(fds_filein);
        if(chidptr!=NULL){
          NewMemory((void **)&chidfilebase,(unsigned int)(strlen(chidptr)+1));
          STRCPY(chidfilebase,chidptr);
        }
      }
      if(chidfilebase!=NULL){
        NewMemory((void **)&hrrfilename,(unsigned int)(strlen(chidfilebase)+8+1));
        STRCPY(hrrfilename,chidfilebase);
        STRCAT(hrrfilename,"_hrr.csv");
        if(STAT(hrrfilename,&statbuffer)!=0){
          FREEMEMORY(hrrfilename);
        }
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ CHID +++++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"CHID",4) == 1){
      if(fgets(buffer,255,stream)==NULL){
        break;
      }
      bufptr=trim_string(buffer);
      len=strlen(bufptr);
      FREEMEMORY(chidfilebase);
      NewMemory((void **)&chidfilebase,(unsigned int)(len+1));
      STRCPY(chidfilebase,bufptr);

      if(chidfilebase!=NULL){
        NewMemory((void **)&hrrfilename,(unsigned int)(strlen(chidfilebase)+8+1));
        STRCPY(hrrfilename,chidfilebase);
        STRCAT(hrrfilename,"_hrr.csv");
        if(STAT(hrrfilename,&statbuffer)!=0){
          FREEMEMORY(hrrfilename);
        }
      }
      continue;
    }

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SLCF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if( (match(buffer,"SLCF",4) == 1)||
        (match(buffer,"SLCC",4) == 1)||
        (match(buffer,"SLFL",4) == 1)||
        (match(buffer,"SLCT",4) == 1)
      ){
      int terrain=0, cellcenter=0, fire_line=0;
      float above_ground_level=0.0;
      slice *sd;

      nn_slice++;
      if(match(buffer,"SLCT",4) == 1){
        terrain=1;
      }
      if(match(buffer,"SLFL",4) == 1){
        terrain=1;
        fire_line=1;
      }
      if(match(buffer,"SLCC",4) == 1){
        cellcenter_slice_active=1;
        cellcenter=1;
      }
      trim(buffer);
      len=strlen(buffer);
      if(nmeshes>1){
        blocknumber=ioffset-1;
      }
      else{
        blocknumber=0;
      }
      if(len>5){
        buffer3=buffer+4;
        sscanf(buffer3,"%i %f",&blocknumber,&above_ground_level);
        blocknumber--;
      }
      if(fgets(buffer,255,stream)==NULL){
        nslice_files--;
        break;
      }

      bufptr=trim_string(buffer);
      len=strlen(bufptr);
      sd = sliceinfo+islice;
      NewMemory((void **)&sd->reg_file,(unsigned int)(len+1));
      STRCPY(sd->reg_file,bufptr);

      NewMemory((void **)&sd->comp_file,(unsigned int)(len+4+1));
      STRCPY(sd->comp_file,bufptr);
      STRCAT(sd->comp_file,".svz");

      NewMemory((void **)&sd->rle_file,(unsigned int)(len+4+1));
      STRCPY(sd->rle_file,bufptr);
      STRCAT(sd->rle_file,".rle");

      NewMemory((void **)&sd->size_file,(unsigned int)(len+3+1));
      STRCPY(sd->size_file,bufptr);
      STRCAT(sd->size_file,".sz");

      sd->compression_type=0;
      if(STAT(sd->rle_file,&statbuffer)==0){
        sd->compression_type=2;
        sd->file=sd->rle_file;
      }
      if(STAT(sd->comp_file,&statbuffer)==0){
        sd->compression_type=1;
        sd->file=sd->comp_file;
      }
      if(sd->compression_type==0){
        sd->file=sd->reg_file;
      }
      sd->fire_line=fire_line;
      sd->terrain=terrain;
      sd->cellcenter=cellcenter;
      sd->above_ground_level=above_ground_level;
      sd->seq_id=nn_slice;
      sd->autoload=0;
      sd->display=0;
      sd->loaded=0;
      sd->qslicedata=NULL;
      sd->compindex=NULL;
      sd->slicecomplevel=NULL;
      sd->qslicedata_compressed=NULL;
      sd->volslice=0;
      sd->slicetimes=NULL;
      sd->slicelevel=NULL;
      sd->slicepoint=NULL;
      sd->slicedata=NULL;
      sd->slicetimeslist=NULL;
      sd->c_iblank=NULL;
      sd->blocknumber=blocknumber;
      sd->vloaded=0;
      sd->reload=0;
      sd->nline_contours=0;
      sd->line_contours=NULL;
      {
        mesh *meshi;

        meshi = meshinfo + blocknumber;
        sd->mesh_type=meshi->mesh_type;
      }
      if(STAT(sd->file,&statbuffer)==0){
        if(sd->terrain==1){
          if(readlabels_terrain(&sd->label,stream)==2)return 2;
        }
        else if(sd->cellcenter==1){
          if(readlabels_cellcenter(&sd->label,stream)==2)return 2;
        }
        else{
          if(readlabels(&sd->label,stream)==2)return 2;
        }
        islice++;
      }
      else{
        if(readlabels(&sd->label,stream)==2)return 2;
        nslice_files--;
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ BNDF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"BNDF",4) == 1||match(buffer,"BNDC",4) == 1){
      patch *patchi;

      nn_patch++;

      trim(buffer);
      len=strlen(buffer);

      if(nmeshes>1){
        blocknumber=ioffset-1;
      }
      else{
        blocknumber=0;
      }
      version=0;
      if(len>5){
        buffer3=buffer+4;
        sscanf(buffer3,"%i %i",&blocknumber,&version);
        blocknumber--;
      }
      patchi = patchinfo + ipatch;

      patchi->version=version;
  
      if(match(buffer,"BNDC",4) == 1){
        patchi->cellcenter=1;
        cellcenter_bound_active=1;
      }
      else{
        patchi->cellcenter=0;
      }

      if(fgets(buffer,255,stream)==NULL){
        npatch_files--;
        break;
      }

      bufptr=trim_string(buffer);
      len=strlen(bufptr);
      NewMemory((void **)&patchi->reg_file,(unsigned int)(len+1));
      STRCPY(patchi->reg_file,bufptr);

      NewMemory((void **)&patchi->comp_file,(unsigned int)(len+4+1));
      STRCPY(patchi->comp_file,bufptr);
      STRCAT(patchi->comp_file,".svz");

      NewMemory((void **)&patchi->size_file,(unsigned int)(len+4+1));
      STRCPY(patchi->size_file,bufptr);
//      STRCAT(patchi->size_file,".szz"); when we actully use file check both .sz and .szz extensions

      if(STAT(patchi->comp_file,&statbuffer)==0){
        patchi->compression_type=1;
        patchi->file=patchi->comp_file;
      }
      else{
        patchi->compression_type=0;
        patchi->file=patchi->reg_file;
      }
      patchi->blocknumber=blocknumber;
      patchi->seq_id=nn_patch;
      patchi->autoload=0;
      patchi->loaded=0;
      patchi->display=0;
      meshinfo[blocknumber].patchfilenum=-1;
      if(STAT(patchi->file,&statbuffer)==0){
        if(patchi->cellcenter==1){
          if(readlabels_cellcenter(&patchi->label,stream)==2)return 2;
        }
        else{
          if(readlabels(&patchi->label,stream)==2)return 2;
        }
#ifdef pp_HIST
        NewMemory((void **)&patchi->histogram,sizeof(histogramdata));
#endif
        ipatch++;
      }
      else{
        if(trainer_mode==0)printf("*** Warning: the file, %s, does not exist.\n",buffer);
        if(readlabels(&patchi->label,stream)==2)return 2;
        npatch_files--;
      }
      continue;
    }

    /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ ISOF ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

    if(match(buffer,"ISOF",4) == 1||match(buffer,"TISOF",5)==1){
      iso *isoi;
      int get_isolevels;

      isoi = isoinfo + iiso;
      nn_iso++;
      dataflag=0;
      if(match(buffer,"TISOF",5)==1)dataflag=1;
      trim(buffer);
      len=strlen(buffer);

      if(nmeshes>1){
        blocknumber=ioffset-1;
      }
      else{
        blocknumber=0;
      }
      if(len>5&&dataflag==0){
        buffer3=buffer+4;
        sscanf(buffer3,"%i",&blocknumber);
        blocknumber--;
      }
      if(len>6&&dataflag==1){
        buffer3=buffer+5;
        sscanf(buffer3,"%i",&blocknumber);
        blocknumber--;
      }
      if(fgets(buffer,255,stream)==NULL){
        niso_files--;
        break;
      }

      isoi->seq_id=nn_iso;
      isoi->autoload=0;
      isoi->blocknumber=blocknumber;
      isoi->loaded=0;
      isoi->display=0;
      isoi->dataflag=dataflag;
      isoi->nlevels=0;
      isoi->levels=NULL;

      isoi->normaltable=NULL;
      isoi->comp_buffer=NULL;
      isoi->comp_bufferframe=NULL;
      isoi->full_bufferframe=NULL;
      isoi->color_label.longlabel=NULL;
      isoi->color_label.shortlabel=NULL;
      isoi->color_label.unit=NULL;

      bufptr=trim_string(buffer);

      len=strlen(bufptr);

      NewMemory((void **)&isoi->reg_file,(unsigned int)(len+1));
      STRCPY(isoi->reg_file,bufptr);

      NewMemory((void **)&isoi->comp_file,(unsigned int)(len+4+1));
      STRCPY(isoi->comp_file,bufptr);
      STRCAT(isoi->comp_file,".svz");

      NewMemory((void **)&isoi->size_file,(unsigned int)(len+3+1));
      STRCPY(isoi->size_file,bufptr);
      STRCAT(isoi->size_file,".sz");

      if(STAT(isoi->comp_file,&statbuffer)==0){
        get_isolevels=1;
        isoi->compression_type=1;
        niso_compressed++;
        isoi->file=isoi->comp_file;
        if(readlabels(&isoi->surface_label,stream)==2)return 2;
        getcisolevels(isoi->file,&isoi->levels,&isoi->nlevels);
        if(dataflag==1){
          if(readlabels(&isoi->color_label,stream)==2)return 2;
        }
        iiso++;
      }
      else if(STAT(isoi->reg_file,&statbuffer2)==0){
        get_isolevels=1;
        isoi->compression_type=0;
        isoi->file=isoi->reg_file;
        if(readlabels(&isoi->surface_label,stream)==2)return 2;
        getisolevels(isoi->file,dataflag,&isoi->levels,&isoi->nlevels);
        if(dataflag==1){
          if(readlabels(&isoi->color_label,stream)==2)return 2;
        }
        iiso++;
      }
      else{
        get_isolevels=0;
        if(trainer_mode==0)printf("*** Warning: the file, %s, does not exist.\n",buffer);
        if(readlabels(&isoi->surface_label,stream)==2)return 2;
        if(dataflag==1){
          if(readlabels(&isoi->color_label,stream)==2)return 2;
        }
        niso_files--;
      }
      if(get_isolevels==1){
        int len_clevels;
        char clevels[1024];

        array2string(isoi->levels,isoi->nlevels,clevels);
        len_clevels = strlen(clevels);
        if(len_clevels>0){
          int len_long;
          char *long_label, *unit_label;

          long_label = isoi->surface_label.longlabel;
          unit_label = isoi->surface_label.unit;
          len_long = strlen(long_label)+strlen(unit_label)+len_clevels+3+1;
          if(dataflag==1)len_long+=(strlen(isoi->color_label.longlabel)+15+1);
          ResizeMemory((void **)&long_label,(unsigned int)len_long);
          isoi->surface_label.longlabel=long_label;
          strcat(long_label,": ");
          strcat(long_label,clevels);
          strcat(long_label," ");
          strcat(long_label,unit_label);
          if(dataflag==1){
            strcat(long_label," (Colored by: ");
            strcat(long_label,isoi->color_label.longlabel);
            strcat(long_label,")");
          }
        }
      }
      continue;
    }

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ THCP ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"THCP",4) == 1){
      mesh *meshi;
      float normdenom;
      char *device_label;

      if(ioffset==0)ioffset=1;
      meshi=meshinfo + ioffset - 1;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&tempval);
      if(tempval<0)tempval=0;
      meshi->ntc=tempval;
      ntc_total += meshi->ntc;
      hasSensorNorm=0;
      if(meshi->ntc>0){
        for(nn=0;nn<meshi->ntc;nn++){
          float *xyz, *xyznorm;
          fgets(buffer,255,stream);

          xyz = devicecopy->xyz;
          xyznorm = devicecopy->xyznorm;
          xyz[0]=0.0;
          xyz[1]=0.0;
          xyz[2]=0.0;
          xyznorm[0]=0.0;
          xyznorm[1]=0.0;
          xyznorm[2]=-1.0;
          device_label=get_device_label(buffer);
          sscanf(buffer,"%f %f %f %f %f %f",xyz,xyz+1,xyz+2,xyznorm,xyznorm+1,xyznorm+2);
          normdenom=0.0;
          normdenom+=xyznorm[0]*xyznorm[0];
          normdenom+=xyznorm[1]*xyznorm[1];
          normdenom+=xyznorm[2]*xyznorm[2];
          if(normdenom>0.1){
            hasSensorNorm=1;
            normdenom=sqrt(normdenom);
            xyznorm[0]/=normdenom;
            xyznorm[1]/=normdenom;
            xyznorm[2]/=normdenom;
          }
          if(device_label==NULL){
            if(isZoneFireModel==1){
              devicecopy->object = get_SVOBJECT_type("target",thcp_object_backup);
            }
            else{
              devicecopy->object = get_SVOBJECT_type("thermoc4",thcp_object_backup);
            }
          }
          else{
            devicecopy->object = get_SVOBJECT_type(device_label,thcp_object_backup);
          }
          get_elevaz(xyznorm,&devicecopy->dtheta,devicecopy->rotate_axis);
    
          init_device(devicecopy,xyz,xyznorm,0,0,NULL,NULL);

          devicecopy++;
          ndeviceinfo++;
        }
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SPRK ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"SPRK",4) == 1 && match(buffer,"SPRK_ACT",8) != 1){
      mesh *meshi;
      char *device_label;
      meshi=meshinfo + ioffset - 1;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&tempval);
      if(tempval<0)tempval=0;
      meshi->nspr=tempval;
      nspr_total += meshi->nspr;
      if(meshi->nspr>0){
        FREEMEMORY(meshi->xspr); FREEMEMORY(meshi->yspr); FREEMEMORY(meshi->zspr); FREEMEMORY(meshi->tspr);
        if(NewMemory((void **)&meshi->xspr,meshi->nspr*sizeof(float))==0||
           NewMemory((void **)&meshi->yspr,meshi->nspr*sizeof(float))==0||
           NewMemory((void **)&meshi->zspr,meshi->nspr*sizeof(float))==0||
           NewMemory((void **)&meshi->tspr,meshi->nspr*sizeof(float))==0||
           NewMemory((void **)&meshi->xsprplot,meshi->nspr*sizeof(float))==0||
           NewMemory((void **)&meshi->ysprplot,meshi->nspr*sizeof(float))==0||
           NewMemory((void **)&meshi->zsprplot,meshi->nspr*sizeof(float))==0)return 2;
        for(nn=0;nn<meshi->nspr;nn++){meshi->tspr[nn]=99999.;}
        xsprcopy=meshi->xspr;
        ysprcopy=meshi->yspr;
        zsprcopy=meshi->zspr;
        for(nn=0;nn<meshi->nspr;nn++){
          float *xyznorm;
          float normdenom;

          fgets(buffer,255,stream);
          xyznorm = devicecopy->xyznorm;
          xyznorm[0]=0.0;
          xyznorm[1]=0.0;
          xyznorm[2]=-1.0;
          device_label=get_device_label(buffer);
          sscanf(buffer,"%f %f %f %f %f %f",xsprcopy,ysprcopy,zsprcopy,xyznorm,xyznorm+1,xyznorm+2);
          devicecopy->act_time=-1.0;
          devicecopy->type = DEVICE_SPRK;
          devicecopy->xyz[0]=*xsprcopy;
          devicecopy->xyz[1]=*ysprcopy;
          devicecopy->xyz[2]=*zsprcopy;
          normdenom=0.0;
          normdenom+=xyznorm[0]*xyznorm[0];
          normdenom+=xyznorm[1]*xyznorm[1];
          normdenom+=xyznorm[2]*xyznorm[2];
          normdenom=sqrt(normdenom);
          if(normdenom>0.001){
            xyznorm[0]/=normdenom;
            xyznorm[1]/=normdenom;
            xyznorm[2]/=normdenom;
          }
          if(device_label==NULL){
            devicecopy->object = get_SVOBJECT_type("sprinkler_upright",sprinkler_upright_object_backup);
          }
          else{
            devicecopy->object = get_SVOBJECT_type(device_label,sprinkler_upright_object_backup);
          }
          get_elevaz(xyznorm,&devicecopy->dtheta,devicecopy->rotate_axis);
    
          init_device(devicecopy,NULL,xyznorm,0,0,NULL,NULL);

          devicecopy++;
          ndeviceinfo++;

          xsprcopy++; ysprcopy++; zsprcopy++;
        }
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SPRK_ACT ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"SPRK_ACT",8) == 1){
      mesh *meshi;

      if(nmeshes>1){
        blocknumber=ioffset-1;
      }
      else{
        blocknumber=0;
      }
      if(strlen(buffer)>9){
        sscanf(buffer,"%s %i",buffer2,&blocktemp);
        if(blocktemp>0&&blocktemp<=nmeshes)blocknumber = blocktemp-1;
      }
      meshi=meshinfo + blocknumber;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f",&nn,&time);
      if(meshi->tspr!=NULL && nn <= meshi->nspr && nn > 0){
        meshi->tspr[nn-1]=time;
        {
          int idev;
          int count=0;

          for(idev=0;idev<ndeviceinfo;idev++){
            device *devicei;

            devicei = deviceinfo + idev;
            if(devicei->type==DEVICE_SPRK){
              count++;
              if(nn==count){
                devicei->act_time=time;
                break;
              }
            }
          }
        }
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ HEAT ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(matchonly(buffer,"HEAT",4) == 1){
      mesh *meshi;
      char *device_label;

      meshi=meshinfo + ioffset - 1;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&tempval);
      if(tempval<0)tempval=0;
      meshi->nheat=tempval;
      nheat_total += meshi->nheat;
      if(meshi->nheat>0){
        FREEMEMORY(meshi->xheat); FREEMEMORY(meshi->yheat); FREEMEMORY(meshi->zheat); FREEMEMORY(meshi->theat);
        FREEMEMORY(meshi->xheatplot); FREEMEMORY(meshi->yheatplot); FREEMEMORY(meshi->zheatplot);
        if(NewMemory((void **)&meshi->xheat,meshi->nheat*sizeof(float))==0||
           NewMemory((void **)&meshi->yheat,meshi->nheat*sizeof(float))==0||
           NewMemory((void **)&meshi->zheat,meshi->nheat*sizeof(float))==0||
           NewMemory((void **)&meshi->theat,meshi->nheat*sizeof(float))==0||
           NewMemory((void **)&meshi->xheatplot,meshi->nheat*sizeof(float))==0||
           NewMemory((void **)&meshi->yheatplot,meshi->nheat*sizeof(float))==0||
           NewMemory((void **)&meshi->zheatplot,meshi->nheat*sizeof(float))==0)return 2;
        for(nn=0;nn<meshi->nheat;nn++){
          meshi->theat[nn]=99999.;
        }
        xheatcopy=meshi->xheat;
        yheatcopy=meshi->yheat;
        zheatcopy=meshi->zheat;
        for(nn=0;nn<meshi->nheat;nn++){
          float *xyznorm;
          float normdenom;
          fgets(buffer,255,stream);
          xyznorm=devicecopy->xyznorm;
          xyznorm[0]=0.0;
          xyznorm[1]=0.0;
          xyznorm[2]=-1.0;
          device_label=get_device_label(buffer);
          sscanf(buffer,"%f %f %f %f %f %f",xheatcopy,yheatcopy,zheatcopy,xyznorm,xyznorm+1,xyznorm+2);
          devicecopy->type = DEVICE_HEAT;
          devicecopy->act_time=-1.0;
          devicecopy->xyz[0]=*xheatcopy;
          devicecopy->xyz[1]=*yheatcopy;
          devicecopy->xyz[2]=*zheatcopy;
          normdenom=0.0;
          normdenom+=xyznorm[0]*xyznorm[0];
          normdenom+=xyznorm[1]*xyznorm[1];
          normdenom+=xyznorm[2]*xyznorm[2];
          normdenom=sqrt(normdenom);
          if(normdenom>0.001){
            xyznorm[0]/=normdenom;
            xyznorm[1]/=normdenom;
            xyznorm[2]/=normdenom;
          }
          if(device_label==NULL){
            devicecopy->object = get_SVOBJECT_type("heat_detector",heat_detector_object_backup);
          }
          else{
            devicecopy->object = get_SVOBJECT_type(device_label,heat_detector_object_backup);
          }
          get_elevaz(xyznorm,&devicecopy->dtheta,devicecopy->rotate_axis);

          init_device(devicecopy,NULL,xyznorm,0,0,NULL,NULL);

          devicecopy++;
          ndeviceinfo++;
          xheatcopy++; yheatcopy++; zheatcopy++;

        }
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ HEAT_ACT ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"HEAT_ACT",8) == 1){
      mesh *meshi;

      if(nmeshes>1){
        blocknumber=ioffset-1;
      }
      else{
        blocknumber=0;
      }
      if(strlen(buffer)>9){
        sscanf(buffer,"%s %i",buffer2,&blocktemp);
        if(blocktemp>0&&blocktemp<=nmeshes)blocknumber = blocktemp-1;
      }
      meshi=meshinfo + blocknumber;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f",&nn,&time);
      if(meshi->theat!=NULL && nn <= meshi->nheat && nn > 0){
        meshi->theat[nn-1]=time;
        {
          int idev;
          int count=0;

          for(idev=0;idev<ndeviceinfo;idev++){
            device *devicei;

            devicei = deviceinfo + idev;
            if(devicei->type==DEVICE_HEAT){
              count++;
              if(nn==count){
                devicei->act_time=time;
                break;
              }
            }
          }
        }
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SMOD ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"SMOD",4) == 1 && match(buffer,"SMOD_ACT",8) != 1){
      float xyz[3];
      int sdnum;
      char *device_label;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&sdnum);
      if(sdnum<0)sdnum=0;
      for(nn=0;nn<sdnum;nn++){
        float *xyznorm;
        float normdenom;

        xyznorm = devicecopy->xyznorm;
        xyznorm[0]=0.0;
        xyznorm[1]=0.0;
        xyznorm[2]=-1.0;
        device_label=get_device_label(buffer);
        fgets(buffer,255,stream);
        sscanf(buffer,"%f %f %f %f %f %f",xyz,xyz+1,xyz+2,xyznorm,xyznorm+1,xyznorm+2);
        devicecopy->type = DEVICE_SMOKE;
        devicecopy->act_time=-1.0;
        devicecopy->xyz[0]=xyz[0];
        devicecopy->xyz[1]=xyz[1];
        devicecopy->xyz[2]=xyz[2];
        normdenom=0.0;
        normdenom+=xyznorm[0]*xyznorm[0];
        normdenom+=xyznorm[1]*xyznorm[1];
        normdenom+=xyznorm[2]*xyznorm[2];
        normdenom=sqrt(normdenom);
        if(normdenom>0.001){
          xyznorm[0]/=normdenom;
          xyznorm[1]/=normdenom;
          xyznorm[2]/=normdenom;
        }
        if(device_label==NULL){
          devicecopy->object = get_SVOBJECT_type("smoke_detector",smoke_detector_object_backup);
        }
        else{
          devicecopy->object = get_SVOBJECT_type(device_label,smoke_detector_object_backup);
        }
        get_elevaz(xyznorm,&devicecopy->dtheta,devicecopy->rotate_axis);
    
        init_device(devicecopy,xyz,xyznorm,0,0,NULL,NULL);

        devicecopy++;
        ndeviceinfo++;

      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SMOD_ACT ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"SMOD_ACT",8) == 1){
      int idev;
      int count=0;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f",&nn,&time);
      for(idev=0;idev<ndeviceinfo;idev++){
        device *devicei;

        devicei = deviceinfo + idev;
        if(devicei->type==DEVICE_SMOKE){
          count++;
          if(nn==count){
            devicei->act_time=time;
            break;
          }
        }
      }
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SHOW_OBST ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"SHOW_OBST",9) == 1||match(buffer,"HIDE_OBST",9)==1){
      mesh *meshi;

      do_pass4=1;
      if(nmeshes>1){
        blocknumber=ioffset-1;
      }
      else{
        blocknumber=0;
      }
      if(strlen(buffer)>10){
        sscanf(buffer,"%s %i",buffer2,&blocktemp);
        if(blocktemp>0&&blocktemp<=nmeshes)blocknumber = blocktemp-1;
      }
      showobst=0;
      if(match(buffer,"SHOW_OBST",9) == 1)showobst=1;
      meshi=meshinfo + blocknumber;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f",&tempval,&time);
      tempval--;
      if(tempval<0||tempval>=meshi->nbptrs)continue;
      nn=tempval;
      bc=meshi->blockageinfoptrs[nn];
      bc->nshowtime++;
      meshi->nsmoothblockages_list++;
      if(meshi->nsmoothblockages_list>0)use_menusmooth=1;
      continue;
    }
  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ OPEN_VENT ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"OPEN_VENT",9) == 1||match(buffer,"CLOSE_VENT",10)==1){
      mesh *meshi;

      do_pass4=1;
      showvent=1;
      if(match(buffer,"CLOSE_VENT",10) == 1)showvent=0;
      if(nmeshes>1){
        blocknumber=ioffset-1;
      }
      else{
        blocknumber=0;
      }
      trim(buffer);
      len=strlen(buffer);
      if(showvent==1){
        if(len>10){
          sscanf(buffer+10,"%i",&blocknumber);
          blocknumber--;
          if(blocknumber<0)blocknumber=0;
          if(blocknumber>nmeshes-1)blocknumber=nmeshes-1;
        }
      }
      else{
        if(len>11){
          sscanf(buffer+11,"%i",&blocknumber);
          blocknumber--;
          if(blocknumber<0)blocknumber=0;
          if(blocknumber>nmeshes-1)blocknumber=nmeshes-1;
        }
      }
      meshi=meshinfo + blocknumber;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f",&tempval,&time);
      tempval--;
      if(meshi->ventinfo==NULL)continue;
      if(tempval<0||tempval>=meshi->nvents)continue;
      nn=tempval;
      vi=meshi->ventinfo+nn;
      vi->nshowtime++;
      continue;
    }
  }

/* 
   ************************************************************************
   ************************ end of pass 3 ********************************* 
   ************************************************************************
 */

  if(autoterrain==1){
    nobst=0;
    iobst=0;
    for(i=0;i<nmeshes;i++){
      mesh *meshi;
      int ibar, jbar, kbar;

      meshi = meshinfo + i;
      ibar = meshi->ibar;
      jbar = meshi->jbar;
      if(ibar>0&&jbar>0){
        float *zcell;

        NewMemory((void **)&meshi->zcell,ibar*jbar*sizeof(float));
        zcell = meshi->zcell;
        for(j=0;j<ibar*jbar;j++){
          zcell[j]=meshi->zbar0; // assume initially that zcell (terrain height) is at base of the domain
        }
      }
    }
  }

  for(i=0;i<nmeshes;i++){
    mesh *meshi;
    int nlist;

    meshi=meshinfo+i;

    nlist=meshi->nsmoothblockages_list+2; // add an entry for t=0.0
    NewMemory((void **)&meshi->smoothblockages_list,nlist*sizeof(smoothblockage));

    meshi->nsmoothblockages_list++;

    for(j=0;j<nlist;j++){
      smoothblockage *sb;

      sb=meshi->smoothblockages_list+j;
      sb->smoothblockagecolors=NULL;
      sb->smoothblockagesurfaces=NULL;
      sb->nsmoothblockagecolors=0;
    }
    meshi->nsmoothblockages_list=1;
    meshi->smoothblockages_list[0].time=-1.0;
    meshi->nsmoothblockages_list++;
  }

  /* 
   ************************************************************************
   ************************ start of pass 4 ********************************* 
   ************************************************************************
 */

  rewind(stream);
  printf("   pass 3 completed\n");
  if(do_pass4==1||autoterrain==1)printf("   pass 4\n");

  while((autoterrain==1||do_pass4==1)&&!feof(stream)){
    if(fgets(buffer,255,stream)==NULL)break;
    if(strncmp(buffer," ",1)==0||buffer[0]==0)continue;

    if(match(buffer,"MINMAXBNDF",10) == 1){
      char *file_ptr, file2[1024];
      float valmin, valmax;
      float percentile_min, percentile_max;

      fgets(buffer,255,stream);
      strcpy(file2,buffer);
      file_ptr = file2;
      trim(file2);
      file_ptr = trim_front(file2);

      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f %f",&valmin,&valmax,&percentile_min,&percentile_max);

      for(i=0;i<npatch_files;i++){
        patch *patchi;

        patchi = patchinfo + i;
        if(strcmp(file_ptr,patchi->file)==0){
          patchi->diff_valmin=percentile_min;
          patchi->diff_valmax=percentile_max;
          break;
        }
      }
      continue;
    }
    if(match(buffer,"MINMAXPL3D",10) == 1){
      char *file_ptr, file2[1024];
      float valmin[5], valmax[5];
      float percentile_min[5], percentile_max[5];

      fgets(buffer,255,stream);
      strcpy(file2,buffer);
      file_ptr = file2;
      trim(file2);
      file_ptr = trim_front(file2);

      for(i=0;i<5;i++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%f %f %f %f",valmin +i,valmax+i, percentile_min+i,percentile_max+i);
      }

      for(i=0;i<nplot3d_files;i++){
        plot3d *plot3di;

        plot3di = plot3dinfo + i;
        if(strcmp(file_ptr,plot3di->file)==0){
          for(j=0;j<5;j++){
            plot3di->diff_valmin[j]=percentile_min[j];
            plot3di->diff_valmax[j]=percentile_max[j];
          }
          break;
        }
      }
      continue;
    }
    if(match(buffer,"MINMAXSLCF",10) == 1){
      char *file_ptr, file2[1024];
      float valmin, valmax;
      float percentile_min, percentile_max;

      fgets(buffer,255,stream);
      strcpy(file2,buffer);
      file_ptr = file2;
      trim(file2);
      file_ptr = trim_front(file2);

      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f %f",&valmin,&valmax,&percentile_min,&percentile_max);

      for(i=0;i<nslice_files;i++){
        slice *slicei;

        slicei = sliceinfo + i;
        if(strcmp(file_ptr,slicei->file)==0){
          slicei->diff_valmin=percentile_min;
          slicei->diff_valmax=percentile_max;
          break;
        }
      }
      continue;
    }
    if(match(buffer,"OBST",4) == 1&&autoterrain==1){
      mesh *meshi;
      int nxcell;

      nobst++;
      iobst++;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&tempval);


      if(tempval<=0)tempval=0;
      if(tempval==0)continue;
      nbtemp=tempval;

      meshi=meshinfo+iobst-1;

      for(nn=0;nn<nbtemp;nn++){
        fgets(buffer,255,stream);
      }
      nxcell = meshi->ibar;
      for(nn=0;nn<nbtemp;nn++){
        int ijk2[5],kmax;
        int ii, jj;
        float zval;

        fgets(buffer,255,stream);
        sscanf(buffer,"%i %i %i %i %i %i",ijk2,ijk2+1,ijk2+2,ijk2+3,ijk2+4,&kmax);
        zval = meshi->zplt[kmax];
        for(ii=ijk2[0];ii<ijk2[1];ii++){
          for(jj=ijk2[2];jj<ijk2[3];jj++){
            int ij;
            float zcell;

            ij = ijcell2(ii,jj);
            zcell = meshi->zcell[ij];
            if(zval>zcell){
              meshi->zcell[ij]=zval;
            }
          }
        }
      }
      continue;
    }

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ DEVICE_ACT ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

    if(match(buffer,"DEVICE_ACT",10) == 1){
      device *devicei;
      int idevice;
      float act_time;
      int act_state=1;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i",&idevice,&act_time,&act_state);
      idevice--;
      if(idevice>=0&&idevice<ndeviceinfo){
        int istate;

        devicei = deviceinfo + idevice;
        devicei->act_time=act_time;
        if(devicei->act_times==NULL){
          devicei->nstate_changes++;
          NewMemory((void **)&devicei->act_times,devicei->nstate_changes*sizeof(int));
          NewMemory((void **)&devicei->state_values,devicei->nstate_changes*sizeof(int));
          devicei->act_times[0]=0.0;
          devicei->state_values[0]=devicei->state0;
          devicei->istate_changes++;
        }
        istate = devicei->istate_changes++;
        devicei->act_times[istate]=act_time;
        devicei->state_values[istate]=act_state;
      }

      continue;
    }

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ SHOW_OBST ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

    if(match(buffer,"SHOW_OBST",9) == 1||match(buffer,"HIDE_OBST",9)==1){
      mesh *meshi;
      int nlist;

      if(nmeshes>1){
        blocknumber=ioffset-1;
      }
      else{
        blocknumber=0;
      }
      if(strlen(buffer)>10){
        sscanf(buffer,"%s %i",buffer2,&blocktemp);
        if(blocktemp>0&&blocktemp<=nmeshes)blocknumber = blocktemp-1;
      }
      showobst=0;
      if(match(buffer,"SHOW_OBST",9) == 1)showobst=1;
/*
    implement a better strategy for setting this variable

      if(match(buffer,"HIDE_OBST",9) == 1)show_slice_in_obst=1;
      */
      meshi=meshinfo + blocknumber;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f",&tempval,&time);
      tempval--;
      if(tempval<0||tempval>=meshi->nbptrs)continue;
      nn=tempval;
      bc=meshi->blockageinfoptrs[nn];

      meshi->nsmoothblockages_list++;
      nlist=meshi->nsmoothblockages_list;
      if(nlist==2&&time!=0.0){          //  insert time=0.0 into list if not there
        meshi->smoothblockages_list[1].time=0.0;
        meshi->nsmoothblockages_list++;
        nlist++;
      }
      meshi->smoothblockages_list[nlist-1].time=time;
      if(time==meshi->smoothblockages_list[nlist-2].time){
        nlist--;
        meshi->nsmoothblockages_list=nlist;
      }

      if(bc->showtime==NULL){
        if(time!=0.0)bc->nshowtime++;
        NewMemory((void **)&bc->showtime,bc->nshowtime*sizeof(float));
        NewMemory((void **)&bc->showhide,bc->nshowtime*sizeof(unsigned char));
        bc->nshowtime=0;
        if(time!=0.0){
          bc->nshowtime=1;
          bc->showtime[0]=0.0;
          if(showobst==1){
            bc->showhide[0]=0;
          }
          else{
            bc->showhide[0]=1;
          }
        }
      }
      bc->nshowtime++;
      if(showobst==1){
        bc->showhide[bc->nshowtime-1]=1;
      }
      else{
        bc->showhide[bc->nshowtime-1]=0;
      }
      bc->showtime[bc->nshowtime-1]=time;
      continue;
    }

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ OPEN_VENT ++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */
    if(match(buffer,"OPEN_VENT",9) == 1||match(buffer,"CLOSE_VENT",10)==1){
      mesh *meshi;

      showvent=1;
      if(match(buffer,"CLOSE_VENT",10) == 1)showvent=0;
      if(nmeshes>1){
        blocknumber=ioffset-1;
      }
      else{
        blocknumber=0;
      }
      trim(buffer);
      len=strlen(buffer);
      if(showvent==1){
        if(len>10){
          sscanf(buffer+10,"%i",&blocknumber);
          blocknumber--;
          if(blocknumber<0)blocknumber=0;
          if(blocknumber>nmeshes-1)blocknumber=nmeshes-1;
        }
      }
      else{
        if(len>11){
          sscanf(buffer+11,"%i",&blocknumber);
          blocknumber--;
          if(blocknumber<0)blocknumber=0;
          if(blocknumber>nmeshes-1)blocknumber=nmeshes-1;
        }
      }
      meshi=meshinfo + blocknumber;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f",&tempval,&time);
      tempval--;
      if(meshi->ventinfo==NULL)continue;
      if(tempval<0||tempval>=meshi->nvents)continue;
      nn=tempval;
      vi=meshi->ventinfo+nn;
      if(vi->showtime==NULL){
        NewMemory((void **)&vi->showtime,(vi->nshowtime+1)*sizeof(float));
        NewMemory((void **)&vi->showhide,(vi->nshowtime+1)*sizeof(unsigned char));
        vi->nshowtime=1;
        vi->showtime[0]=0.0;
        vi->showhide[0]=1;
      }
      if(showvent==1){
        vi->showhide[vi->nshowtime]=1;
      }
      else{
        vi->showhide[vi->nshowtime]=0;
      }
      vi->showtime[vi->nshowtime++]=time;
      continue;
    }

  }

  if(do_pass4==1)printf("   pass 4 completed\n");

  printf("reading input file completed\n");
  printf("beginning wrap up - \n");


/* 
   ************************************************************************
   ************************ wrap up *************************************** 
   ************************************************************************
 */

  CheckMemory;

  update_inilist();

  if(meshinfo!=NULL&&meshinfo->jbar==1){
    force_isometric=1;
  }
  if(hrrfilename!=NULL){
    readhrr(LOAD, &errorcode);
  }

#ifdef pp_THREAD
  if(mt_compress==1)pthread_mutex_init(&mutexCOMPRESS,NULL);
#endif

  init_part5prop();
  init_plot3dtimelist();

  if(noutlineinfo>0){
    highlight_flag=2;
  }
  else{
    highlight_flag=1;
  }
  if(dummyvents==1){
    visFloor=0;
    visCeiling=0;
    visWalls=0;
  }
  initcadcolors();

  /* read in FDS input file and copy OBST labels to obstlabels */
/*
  if(fds_filein!=NULL){
    CheckMemoryOff;
    getlabels(fds_filein);
    CheckMemoryOn;
  }
  */

  FREEMEMORY(slice_loaded_list);
  if(nslice_files>0){
    NewMemory((void **)&slice_loaded_list,nslice_files*sizeof(int));
  }
  FREEMEMORY(patch_loaded_list);
  if(npatch_files>0){
    NewMemory((void **)&patch_loaded_list,npatch_files*sizeof(int));
  }
  update_loaded_lists();

  if(setPDIM==0){
    for(nn=0;nn<=ibartemp;nn++){
      current_mesh->xplt[nn]=xbar0+(float)nn*(xbar-xbar0)/(float)ibartemp;
    }
    for(nn=0;nn<=jbartemp;nn++){
      current_mesh->yplt[nn]=ybar0+(float)nn*(ybar-ybar0)/(float)jbartemp;
    }
    for(nn=0;nn<=kbartemp;nn++){
      current_mesh->zplt[nn]=zbar0+(float)nn*(zbar-zbar0)/(float)kbartemp;
    }
  }

  /* define highlighted block */

  /* add in offsets */
  for(i=0;i<nmeshes;i++){
    mesh *meshi;
    int ii;

    meshi=meshinfo+i;
    meshi->xbar += meshi->offset[0];
    meshi->ybar += meshi->offset[1];
    meshi->zbar += meshi->offset[2];
    meshi->xbar0 += meshi->offset[0];
    meshi->ybar0 += meshi->offset[1];
    meshi->zbar0 += meshi->offset[2];
    {
      float dx, dy, dz;

      dx = meshi->xbar - meshi->xbar0;
      dx /= meshi->ibar;
      dy = meshi->ybar - meshi->ybar0;
      dy /= meshi->jbar;
      dz = meshi->zbar - meshi->zbar0;
      dz /= meshi->kbar;
      meshi->cellsize=sqrt(dx*dx+dy*dy+dz*dz);
    }
    for(ii=0;ii<meshi->ibar+1;ii++){
      meshi->xplt[ii] += meshi->offset[0];
    }
    for(ii=0;ii<meshi->jbar+1;ii++){
      meshi->yplt[ii] += meshi->offset[1];
    }
    for(ii=0;ii<meshi->kbar+1;ii++){
      meshi->zplt[ii] += meshi->offset[2];
    }
    meshi->xcen+=meshi->offset[0];
    meshi->ycen+=meshi->offset[1];
    meshi->zcen+=meshi->offset[2];
  }
  for(i=1;i<nmeshes;i++){
    mesh *meshi;

    meshi=meshinfo+i;
    if(meshi->zbar0!=meshinfo->zbar0){
      visFloor=0;
      updatefacelists=1;
      updatemenu=1;
      break;
    }
  }


  /*
    Associate a surface with each block.
  */
  updateusetextures();

  /* compute global bar's and box's */

  {
    mesh *meshi;

    meshi=meshinfo;

    xbar  =  meshi->xbar;   ybar=meshi->ybar;   zbar=meshi->zbar;
    xbar0 =  meshi->xbar0; ybar0=meshi->ybar0; zbar0=meshi->zbar0;
  }

  ijkbarmax=meshinfo->ibar;
  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi=meshinfo+i;

    if(meshi->ibar>ijkbarmax)ijkbarmax=meshi->ibar;
    if(meshi->jbar>ijkbarmax)ijkbarmax=meshi->jbar;
    if(meshi->kbar>ijkbarmax)ijkbarmax=meshi->kbar;
  }

  for(i=1;i<nmeshes;i++){
    mesh *meshi;

    meshi=meshinfo+i;

    if(xbar <meshi->xbar )xbar =meshi->xbar;
    if(ybar <meshi->ybar )ybar =meshi->ybar;
    if(zbar <meshi->zbar )zbar =meshi->zbar;
    if(xbar0>meshi->xbar0)xbar0=meshi->xbar0;
    if(ybar0>meshi->ybar0)ybar0=meshi->ybar0;
    if(zbar0>meshi->zbar0)zbar0=meshi->zbar0;
  }
  dxbar = (xbar-xbar0)/(float)256;
  dybar = (ybar-ybar0)/(float)256;
  dzbar = (zbar-zbar0)/(float)256;

  factor = 256*128;
  dxsbar = (xbar-xbar0)/factor;
  dysbar = (ybar-ybar0)/factor;
  dzsbar = (zbar-zbar0)/factor;


  for(nn=0;nn<=255;nn++)  {xpltb[nn]=xbar0+((float)nn+0.5)*dxbar;}
  for(nn=0;nn<=255;nn++)  {ypltb[nn]=ybar0+((float)nn+0.5)*dybar;}
  for(nn=0;nn<=255;nn++)  {zpltb[nn]=zbar0+((float)nn+0.5)*dzbar;}
  for(nn=0;nn<factor;nn++){xplts[nn]=xbar0+((float)nn+0.5)*dxsbar;}
  for(nn=0;nn<factor;nn++){yplts[nn]=ybar0+((float)nn+0.5)*dysbar;}
  for(nn=0;nn<factor;nn++){zplts[nn]=zbar0+((float)nn+0.5)*dzsbar;}


  /* compute scaling factors */

  {
    float dxclip, dyclip, dzclip;

    dxclip = (xbar-xbar0)/1000.0;
    dyclip = (ybar-ybar0)/1000.0;
    dzclip = (zbar-zbar0)/1000.0;
    xclip_min = xbar0-dxclip;
    yclip_min = ybar0-dyclip;
    zclip_min = zbar0-dzclip;
    xclip_max = xbar+dxclip;
    yclip_max = ybar+dyclip;
    zclip_max = zbar+dzclip;
  }

  xyzmaxdiff=xbar-xbar0;
  if(ybar-ybar0>xyzmaxdiff)xyzmaxdiff=ybar-ybar0;
  if(zbar-zbar0>xyzmaxdiff)xyzmaxdiff=zbar-zbar0;

  for(i=0;i<npartclassinfo;i++){
    part5class *partclassi;

    partclassi = partclassinfo + i;

    if(partclassi->device_name!=NULL){
        float diameter, length, azimuth, elevation;
       
        partclassi->diameter/=xyzmaxdiff;
        partclassi->length/=xyzmaxdiff;
        length=partclassi->length;
        azimuth = partclassi->azimuth*PI/180.0;
        elevation = partclassi->elevation*PI/180.0;
        partclassi->dx = cos(azimuth)*cos(elevation)*length/2.0;
        partclassi->dy = sin(azimuth)*cos(elevation)*length/2.0;
        partclassi->dz =              sin(elevation)*length/2.0;
    }
  }

#ifdef pp_SHOOTER
  shooter_xyz[0]=xbar/2.0;
  shooter_xyz[1] = 0.0;
  shooter_xyz[2] = zbar/2.0;
  shooter_dxyz[0]=xbar/4.0;
  shooter_dxyz[1]=0.0;
  shooter_dxyz[2]=0.0;
  shooter_nparts=100;
  shooter_velmag=1.0;
  shooter_veldir=0.0;
  shooter_fps=10;
  shooter_vel_type=1;
#endif

  /* rescale both global and local xbar, ybar and zbar */

  xbar0ORIG = xbar0;
  ybar0ORIG = ybar0;
  zbar0ORIG = zbar0;
  xbarORIG = xbar;
  ybarORIG = ybar;
  zbarORIG = zbar;
  xbar = (xbar-xbar0)/xyzmaxdiff;
  ybar = (ybar-ybar0)/xyzmaxdiff;
  zbar = (zbar-zbar0)/xyzmaxdiff;
  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi=meshinfo+i;
    /* compute a local scaling factor for each block */
    meshi->xyzmaxdiff=meshi->xbar-meshi->xbar0;
    if(meshi->ybar-meshi->ybar0>meshi->xyzmaxdiff)meshi->xyzmaxdiff=meshi->ybar-meshi->ybar0;
    if(meshi->zbar-meshi->zbar0>meshi->xyzmaxdiff)meshi->xyzmaxdiff=meshi->zbar-meshi->zbar0;

    meshi->xbar = (meshi->xbar-xbar0)/xyzmaxdiff;
    meshi->ybar = (meshi->ybar-ybar0)/xyzmaxdiff;
    meshi->zbar = (meshi->zbar-zbar0)/xyzmaxdiff;
    meshi->xcen = (meshi->xcen-xbar0)/xyzmaxdiff;
    meshi->ycen = (meshi->ycen-ybar0)/xyzmaxdiff;
    meshi->zcen = (meshi->zcen-zbar0)/xyzmaxdiff;
  }

  for(i=0;i<noutlineinfo;i++){
    outlinei = outlineinfo + i;
    x1 = outlinei->x1;
    x2 = outlinei->x2;
    yy1 = outlinei->y1;
    yy2 = outlinei->y2;
    z1 = outlinei->z1;
    z2 = outlinei->z2;
    for(j=0;j<outlinei->nlines;j++){
      x1[j]=(x1[j]-xbar0)/xyzmaxdiff;
      x2[j]=(x2[j]-xbar0)/xyzmaxdiff;
      yy1[j]=(yy1[j]-ybar0)/xyzmaxdiff;
      yy2[j]=(yy2[j]-ybar0)/xyzmaxdiff;
      z1[j]=(z1[j]-zbar0)/xyzmaxdiff;
      z2[j]=(z2[j]-zbar0)/xyzmaxdiff;
    }
  }

  {
    mesh *meshi;
    meshi=meshinfo;
    veclength = meshi->xplt[1]-meshi->xplt[0];
  }
  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi=meshinfo+i;
    if(veclength>meshi->xplt[1]-meshi->xplt[0])veclength=meshi->xplt[1]-meshi->xplt[0];
    if(veclength>meshi->yplt[1]-meshi->yplt[0])veclength=meshi->yplt[1]-meshi->yplt[0];
    if(veclength>meshi->zplt[1]-meshi->zplt[0])veclength=meshi->zplt[1]-meshi->zplt[0];
  }
  veclength = veclength/xyzmaxdiff;
  veclength = 0.01;

  for(igrid=0;igrid<nmeshes;igrid++){
    mesh *meshi;

    meshi=meshinfo+igrid;
    ibartemp=meshi->ibar; jbartemp=meshi->jbar; kbartemp=meshi->kbar;
    xplt_origtemp = meshi->xplt_orig;
    yplt_origtemp = meshi->yplt_orig;
    zplt_origtemp = meshi->zplt_orig;
    xplttemp = meshi->xplt;
    yplttemp = meshi->yplt;
    zplttemp = meshi->zplt;

    for(i=0;i<ibartemp+1;i++){
      xplt_origtemp[i]=xplttemp[i];
      xplttemp[i]=(xplttemp[i]-xbar0)/xyzmaxdiff;
    }
    for(j=0;j<jbartemp+1;j++){
      yplt_origtemp[j]=yplttemp[j];
      yplttemp[j]=(yplttemp[j]-ybar0)/xyzmaxdiff;
    }
    for(k=0;k<kbartemp+1;k++){
      zplt_origtemp[k]=zplttemp[k];
      zplttemp[k]=(zplttemp[k]-zbar0)/xyzmaxdiff;
    }
    meshi->boxoffset=-(zplttemp[1]-zplttemp[0])/10.0;
    meshi->boxmin[0]=xplt_origtemp[0];
    meshi->boxmin[1]=yplt_origtemp[0];
    meshi->boxmin[2]=zplt_origtemp[0];
    meshi->boxmax[0]=xplt_origtemp[ibartemp];
    meshi->boxmax[1]=yplt_origtemp[jbartemp];
    meshi->boxmax[2]=zplt_origtemp[kbartemp];

  }

  active_smokesensors=0;
  for(i=0;i<ndeviceinfo;i++){
    device *devicei;
    char *label;

    devicei = deviceinfo + i;
    devicei->device_mesh=get_mesh(devicei->xyz);
    label = devicei->object->label;
    if(strcmp(label,"smokesensor")==0){
      active_smokesensors=1;
    }
    if(devicei->plane_surface!=NULL){
      init_device_plane(devicei);
    }
  }

  nsmoothblocks=0;
  ntransparentblocks=0;
  ntransparentvents=0;
  nopenvents=0;
  nopenvents_nonoutline=0;
  ndummyvents=0;
  for(igrid=0;igrid<nmeshes;igrid++){
    mesh *meshi;

    meshi=meshinfo+igrid;
    for(i=0;i<meshi->nbptrs;i++){
      bc=meshi->blockageinfoptrs[i];
      if(bc->type==BLOCK_smooth)nsmoothblocks++;
      if(bc->color[3]<0.99)ntransparentblocks++;
    }
    for(i=0;i<meshi->nvents;i++){
      vi = meshi->ventinfo + i;
      if(vi->isOpenvent==1){
        nopenvents++;
        if(vi->type!=BLOCK_OUTLINE)nopenvents_nonoutline++;
      }
      if(vi->dummy==1)ndummyvents++;
      if(vi->color[3]<0.99)ntransparentvents++;
    }
  }

  for(igrid=0;igrid<nmeshes;igrid++){
    mesh *meshi;

    meshi=meshinfo+igrid;
    for(i=0;i<meshi->nbptrs;i++){
      bc=meshi->blockageinfoptrs[i];
      bc->xmin += meshi->offset[0];
      bc->xmax += meshi->offset[0];
      bc->ymin += meshi->offset[1];
      bc->ymax += meshi->offset[1];
      bc->zmin += meshi->offset[2];
      bc->zmax += meshi->offset[2];
      bc->xmin = (bc->xmin-xbar0)/xyzmaxdiff;
      bc->xmax = (bc->xmax-xbar0)/xyzmaxdiff;
      bc->ymin = (bc->ymin-ybar0)/xyzmaxdiff;
      bc->ymax = (bc->ymax-ybar0)/xyzmaxdiff;
      bc->zmin = (bc->zmin-zbar0)/xyzmaxdiff;
      bc->zmax = (bc->zmax-zbar0)/xyzmaxdiff;
      bc->xyzORIG[0]=bc->xmin;
      bc->xyzORIG[1]=bc->xmax;
      bc->xyzORIG[2]=bc->ymin;
      bc->xyzORIG[3]=bc->ymax;
      bc->xyzORIG[4]=bc->zmin;
      bc->xyzORIG[5]=bc->zmax;
      bc->ijkORIG[0]=bc->ijk[0];
      bc->ijkORIG[1]=bc->ijk[1];
      bc->ijkORIG[2]=bc->ijk[2];
      bc->ijkORIG[3]=bc->ijk[3];
      bc->ijkORIG[4]=bc->ijk[4];
      bc->ijkORIG[5]=bc->ijk[5];
    }
    for(i=0;i<meshi->nvents+12;i++){
      vi=meshi->ventinfo+i;
      vi->xmin = (vi->xmin-xbar0)/xyzmaxdiff;
      vi->xmax = (vi->xmax-xbar0)/xyzmaxdiff;
      vi->ymin = (vi->ymin-ybar0)/xyzmaxdiff;
      vi->ymax = (vi->ymax-ybar0)/xyzmaxdiff;
      vi->zmin = (vi->zmin-zbar0)/xyzmaxdiff;
      vi->zmax = (vi->zmax-zbar0)/xyzmaxdiff;
    }
  }
  for(igrid=0;igrid<nmeshes;igrid++){
    mesh *meshi;

    meshi=meshinfo+igrid;
    for(i=0;i<meshi->nbptrs;i++){
      bc=meshi->blockageinfoptrs[i];
      backup_blockage(bc);
    }
  }

  {
    cadquad *quadi;

     for(i=0;i<ncadgeom;i++){
      cd=cadgeominfo+i;
      for(j=0;j<cd->nquads;j++){
        quadi = cd->quad+j;
        for(k=0;k<4;k++){
          quadi->xyzpoints[3*k+0] = (quadi->xyzpoints[3*k+0]-xbar0)/xyzmaxdiff;
          quadi->xyzpoints[3*k+1] = (quadi->xyzpoints[3*k+1]-ybar0)/xyzmaxdiff;
          quadi->xyzpoints[3*k+2] = (quadi->xyzpoints[3*k+2]-zbar0)/xyzmaxdiff;
        }
        if(cd->version==2&&quadi->cadlookq->textureinfo.loaded==1){
          update_cadtextcoords(quadi);
        }
      }
    }
  }

  for(nn=0;nn<=255;nn++)  {xpltb[nn]=(xpltb[nn]-xbar0)/xyzmaxdiff;}
  for(nn=0;nn<=255;nn++)  {ypltb[nn]=(ypltb[nn]-ybar0)/xyzmaxdiff;}
  for(nn=0;nn<=255;nn++)  {zpltb[nn]=(zpltb[nn]-zbar0)/xyzmaxdiff;}

  for(nn=0;nn<factor;nn++){xplts[nn]=(xplts[nn]-xbar0)/xyzmaxdiff;}
  for(nn=0;nn<factor;nn++){yplts[nn]=(yplts[nn]-ybar0)/xyzmaxdiff;}
  for(nn=0;nn<factor;nn++){zplts[nn]=(zplts[nn]-zbar0)/xyzmaxdiff;}

  for(n=0;n<nrooms;n++){
    roomdata *roomi;

    roomi = roominfo + n;
    roomi->x0=(roomi->x0-xbar0)/xyzmaxdiff;
    roomi->y0=(roomi->y0-ybar0)/xyzmaxdiff;
    roomi->z0=(roomi->z0-zbar0)/xyzmaxdiff;
    roomi->x1=(roomi->x1-xbar0)/xyzmaxdiff;
    roomi->y1=(roomi->y1-ybar0)/xyzmaxdiff;
    roomi->z1=(roomi->z1-zbar0)/xyzmaxdiff;
    roomi->dx=roomi->dx/xyzmaxdiff;
    roomi->dy=roomi->dy/xyzmaxdiff;
    roomi->dz=roomi->dz/xyzmaxdiff;
  }
  for(n=0;n<nfires;n++){
    fireinfo[n].absx=(fireinfo[n].absx-xbar0)/xyzmaxdiff;
    fireinfo[n].absy=(fireinfo[n].absy-ybar0)/xyzmaxdiff;
    fireinfo[n].absz=(fireinfo[n].absz-zbar0)/xyzmaxdiff;
    fireinfo[n].dz=fireinfo[n].dz/xyzmaxdiff;
  }
  for(n=0;n<nzvents;n++){
    zvent *zvi;

    zvi = zventinfo + n;

    switch (zvi->dir){
    case 1:
    case 3:
      zvi->x1 = (zvi->x1-xbar0)/xyzmaxdiff;
      zvi->x2 = (zvi->x2-xbar0)/xyzmaxdiff;
      zvi->z1 = (zvi->z1-zbar0)/xyzmaxdiff;
      zvi->z2 = (zvi->z2-zbar0)/xyzmaxdiff;
      zvi->yy = (zvi->yy-ybar0)/xyzmaxdiff;
      break;
    case 2:
    case 4:
      zvi->x1 = (zvi->x1-ybar0)/xyzmaxdiff;
      zvi->x2 = (zvi->x2-ybar0)/xyzmaxdiff;
      zvi->z1 = (zvi->z1-zbar0)/xyzmaxdiff;
      zvi->z2 = (zvi->z2-zbar0)/xyzmaxdiff;
      zvi->yy = (zvi->yy-xbar0)/xyzmaxdiff;
      break;
    case 5:
    case 6:
      zvi->x1 = (zvi->x1-xbar0)/xyzmaxdiff;
      zvi->x2 = (zvi->x2-xbar0)/xyzmaxdiff;
      zvi->y1 = (zvi->y1-ybar0)/xyzmaxdiff;
      zvi->y2 = (zvi->y2-ybar0)/xyzmaxdiff;
      zvi->zz = (zvi->zz-zbar0)/xyzmaxdiff;
      break;
    default:
      ASSERT(FFALSE);
      break;
    }
  }

  /* if we have to create smooth block structures do it now,
     otherwise postpone job until we view all blockages as smooth */

//  the code below was moved to after readini in startup (so that sb_atstart will be defined)
//  if(sb_atstart==1){
//    smooth_blockages();
//  }

  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi=meshinfo + i;
    for(n=0;n<meshi->nspr;n++){
      meshi->xsprplot[n]=(meshi->offset[0]+meshi->xspr[n]-xbar0)/xyzmaxdiff;
      meshi->ysprplot[n]=(meshi->offset[1]+meshi->yspr[n]-ybar0)/xyzmaxdiff;
      meshi->zsprplot[n]=(meshi->offset[2]+meshi->zspr[n]-zbar0)/xyzmaxdiff;
    }
    for(n=0;n<meshi->nheat;n++){
      meshi->xheatplot[n]=(meshi->offset[0]+meshi->xheat[n]-xbar0)/xyzmaxdiff;
      meshi->yheatplot[n]=(meshi->offset[1]+meshi->yheat[n]-ybar0)/xyzmaxdiff;
      meshi->zheatplot[n]=(meshi->offset[2]+meshi->zheat[n]-zbar0)/xyzmaxdiff;
    }
    for(n=0;n<meshi->nvents+12;n++){
      vi = meshi->ventinfo+n;
      vi->xvent1plot=(meshi->offset[0]+vi->xvent1-xbar0)/xyzmaxdiff;
      vi->xvent2plot=(meshi->offset[0]+vi->xvent2-xbar0)/xyzmaxdiff;
      vi->yvent1plot=(meshi->offset[1]+vi->yvent1-ybar0)/xyzmaxdiff;
      vi->yvent2plot=(meshi->offset[1]+vi->yvent2-ybar0)/xyzmaxdiff;
      vi->zvent1plot=(meshi->offset[2]+vi->zvent1-zbar0)/xyzmaxdiff;
      vi->zvent2plot=(meshi->offset[2]+vi->zvent2-zbar0)/xyzmaxdiff;
    }
  }

  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi=meshinfo + i;
    meshi->vent_offset[0] = ventoffset_factor*(meshi->xplt[1]-meshi->xplt[0]);
    meshi->vent_offset[1] = ventoffset_factor*(meshi->yplt[1]-meshi->yplt[0]);
    meshi->vent_offset[2] = ventoffset_factor*(meshi->zplt[1]-meshi->zplt[0]);
  }
  makeiblank();
  makeiblank_carve();
  makeiblank_smoke3d();
  setventdirs();
  update_faces();

  xcenGLOBAL=xbar/2.0;  ycenGLOBAL=ybar/2.0; zcenGLOBAL=zbar/2.0;


  xxmax = xbar;
  if(ybar>xxmax)xxmax=ybar;
  if(zbar>xxmax)xxmax=zbar;

  if(setendian==0){
    if(match(LESsystem,"AIX",3)==1){endian_data=1;}
    if(match(LESsystem,"SGI",3)==1){endian_data=1;}
    if(match(LESsystem,"DVF",3)==1){endian_data=0;}
    if(match(LESendian,"l",1)==1||match(LESendian,"L",1)==1){endian_data=0;}
    if(match(LESendian,"b",1)==1||match(LESendian,"B",1)==1){endian_data=1;}
    endian = endian_data;
  }

  if(niso_files>0){
    FREEMEMORY(isoindex);
    FREEMEMORY(isobounds);
    if(NewMemory((void*)&isoindex,niso_files*sizeof(int))==0)return 2;
    if(NewMemory((void*)&isobounds,niso_files*sizeof(databounds))==0)return 2;
    niso_bounds=0;
    for(i=0;i<niso_files;i++){
      iso *isoi;

      isoi = isoinfo + i;
      if(isoi->dataflag==0)continue;
      isoi->firstshort=1;
      isoi->setvalmin=0;
      isoi->setvalmax=0;
      isoi->valmin=1.0;
      isoi->valmax=0.0;
      isoindex[niso_bounds]=i;
      isobounds[niso_bounds].datalabel=isoi->color_label.shortlabel;
      isobounds[niso_bounds].setvalmin=0;
      isobounds[niso_bounds].setvalmax=0;
      isobounds[niso_bounds].valmin=1.0;
      isobounds[niso_bounds].valmax=0.0;
      isobounds[niso_bounds].setchopmax=0;
      isobounds[niso_bounds].setchopmin=0;
      isobounds[niso_bounds].chopmax=0.0;
      isobounds[niso_bounds].chopmin=1.0;
      isobounds[niso_bounds].label=&isoi->color_label;
      niso_bounds++;
      for(n=0;n<i;n++){
        iso *ison;

        ison = isoinfo + n;
        if(ison->dataflag==0)continue;
        if(strcmp(isoi->color_label.shortlabel,ison->color_label.shortlabel)==0){
          isoi->firstshort=0;
          niso_bounds--;
          break;
        }
      }
    }
  }

  if(nslice_files>0){
    FREEMEMORY(sliceindex);
    FREEMEMORY(slicebounds);
    if(NewMemory((void*)&sliceindex,nslice_files*sizeof(int))==0)return 2;
    if(NewMemory((void*)&slicebounds,nslice_files*sizeof(databounds))==0)return 2;
    nslice2=0;
    for(i=0;i<nslice_files;i++){
      slice *slicei;

      slicei = sliceinfo + i;
      slicei->firstshort=1;
      slicei->valmin=1.0;
      slicei->valmax=0.0;
      slicei->setvalmin=0;
      slicei->setvalmax=0;
      sliceindex[nslice2]=i;
      slicebounds[nslice2].datalabel=slicei->label.shortlabel;
      slicebounds[nslice2].setvalmin=0;
      slicebounds[nslice2].setvalmax=0;
      slicebounds[nslice2].valmin=1.0;
      slicebounds[nslice2].valmax=0.0;
      slicebounds[nslice2].chopmax=0.0;
      slicebounds[nslice2].chopmin=1.0;
      slicebounds[nslice2].setchopmax=0;
      slicebounds[nslice2].setchopmin=0;
#ifdef pp_SLICECONTOURS
      slicebounds[nslice2].line_contour_min=0.0;
      slicebounds[nslice2].line_contour_max=1.0;
      slicebounds[nslice2].line_contour_num=1;
#endif
      nslice2++;
      for(n=0;n<i;n++){
        slice *slicen;

        slicen = sliceinfo + n;
        if(strcmp(slicei->label.shortlabel,slicen->label.shortlabel)==0){
          slicei->firstshort=0;
          nslice2--;
          break;
        }
      }
    }
  }
  canshow_threshold=0;
  if(npatch_files>0){
    npatch2=0;
    FREEMEMORY(patchlabellist);
    if(NewMemory((void **)&patchlabellist,npatch_files*sizeof(char *))==0)return 2;
    for(i=0;i<npatch_files;i++){
      patchinfo[i].firstshort=1;
      patchinfo[i].valmin=1.0;
      patchinfo[i].valmax=0.0;
      patchinfo[i].setvalmin=0;
      patchinfo[i].setvalmax=0;
      if(strncmp(patchinfo[i].label.shortlabel,"temp",4)==0||
         strncmp(patchinfo[i].label.shortlabel,"TEMP",4)==0){
        canshow_threshold=1;
      }
      patchlabellist[npatch2]=patchinfo[i].label.shortlabel;
      npatch2++;
      for(n=0;n<i;n++){
        if(strcmp(patchinfo[i].label.shortlabel,patchinfo[n].label.shortlabel)==0){
          patchinfo[i].firstshort=0;
          npatch2--;
          break;
        }
      }
    }
  }
  updatechar();

  /* make static blockage iso-surfaces */

#ifndef WIN32
  if(endian!=getendian()){
    printf("*** Warning: Smokeview is running on a ");
    if(getendian()==1){
      printf(" little endian computer\n");
    }
    else{
      printf(" big endian computer\n");
    }
    printf("    but the data being visualized was generated on a ");
    if(endian==1){
      printf(" little endian computer\n");
    }
    else{
      printf(" big endian computer\n");
    }
  }
#endif


  fclose(stream);
  stream=NULL;
  update_selectfaces();
  updateslicetypes();
  updatesliceboundlabels();
  updateisotypes();
  updatepatchtypes();
  if(autoterrain==1){
    for(i=0;i<nmeshes;i++){
      mesh *meshi;
      float *zcell;

      meshi = meshinfo + i;
      zcell = meshi->zcell;

      zterrain_min = zcell[0];
      zterrain_max = zterrain_min;
      for(j=1;j<meshi->ibar*meshi->jbar;j++){
        float zval;

        zval=*zcell++;
        if(zval<zterrain_min)zterrain_min=zval;
        if(zval>zterrain_max)zterrain_max=zval;
      }
    }
  }
  update_terrain(1,vertical_factor);
  update_terrain_colors();
  updatevslices();
  updatesmoke3dmenulabels();
  updatevslicetypes();
  updatepatchmenulabels();
  updateisomenulabels();
  updateplot3dmenulabels();
  updatepartmenulabels();
  updatetourmenulabels();
  init_user_ticks();
  clip_I=ibartemp; clip_J=jbartemp; clip_K=kbartemp;

  // define changed_idlist used for blockage editing

  {
    int ntotal=0;

    for(i=0;i<nmeshes;i++){
      mesh *meshi;

      meshi = meshinfo + i;
      ntotal += meshi->nbptrs;
    }
    FREEMEMORY(changed_idlist);

    NewMemory((void **)&changed_idlist,sizeof(int)*(ntotal+1));

    for(i=0;i<ntotal;i++){
      changed_idlist[i]=0;
    }
    nchanged_idlist=ntotal;
  }
#ifdef pp_CULL

  // define data structures used to speed up 3d smoke drawing (by culling non-visible smoke planes)

  if(cullactive==1)initcull(cullsmoke);
#endif


  printf("wrap up completed\n");

  
  return 0;
}


/* ------------------ parsedatabase ------------------------ */

void parsedatabase(char *file){
  FILE *stream;
  char buffer[1000],*buffer2=NULL,*buffer3,*slashptr;
  size_t lenbuffer,lenbuffer2;
  size_t sizebuffer2;
  char *surf_id=NULL,*start,*surf_id2;
  char *c;
  int i,j;
  surface *surfj;
  char *labeli, *labelj;
  int nexti;
  int nsurfids_shown;

  /* free memory called before */

  for(i=0;i<nsurfids;i++){
    surf_id=surfids[i].label;
    FREEMEMORY(surf_id);
  }
  FREEMEMORY(surfids);
  nsurfids=0;


  if( file==NULL||strlen(file)==0||(stream=fopen(file,"r"))==NULL){
    NewMemory((void **)&surfids,(nsurfids+1)*sizeof(surfid));
    surf_id=NULL;
    NewMemory((void **)&surf_id,6);
    strcpy(surf_id,"INERT");
    surfids[0].label=surf_id;
    surfids[0].location=0;
    surfids[0].show=1;
    nsurfids=1;
  }

  else{
    sizebuffer2=1001;
    NewMemory((void **)&buffer2,sizebuffer2);

  /* find out how many surfs are in database file so memory can be allocated */
  
    while(!feof(stream)){
      if(fgets(buffer,1000,stream)==NULL)break;
      if(STRSTR(buffer,"&SURF")==NULL)continue;


      slashptr=strstr(buffer,"/");
      if(slashptr!=NULL)strcpy(buffer2,buffer);
      buffer3=buffer;
      while(slashptr!=NULL){
        fgets(buffer,1000,stream);
        lenbuffer=strlen(buffer);
        lenbuffer2=strlen(buffer2);
        if(lenbuffer2+lenbuffer+2>sizebuffer2){
          sizebuffer2 = lenbuffer2+lenbuffer+2+1000;
          ResizeMemory((void **)&buffer2,(unsigned int)sizebuffer2);
        }
        strcat(buffer2,buffer);
        slashptr=strstr(buffer,"/");
        buffer3=buffer2;
      }
      start=STRSTR(buffer3,"ID");
      if(start!=NULL)nsurfids++;
    }

  /* allocate memory */

    NewMemory((void **)&surfids,(nsurfids+1)*sizeof(surfid));
    surf_id=NULL;
    NewMemory((void **)&surf_id,6);
    strcpy(surf_id,"INERT");
    surfids[0].label=surf_id;
    surfids[0].location=0;
    surfids[0].show=1;


  /* now look for IDs and copy them into an array */

    rewind(stream);
    nsurfids=1;
    while(!feof(stream)){
      if(fgets(buffer,1000,stream)==NULL)break;
      if(STRSTR(buffer,"&SURF")==NULL)continue;
   

      slashptr=strstr(buffer,"/");
      if(slashptr!=NULL)strcpy(buffer2,buffer);
      while(slashptr!=NULL){
        fgets(buffer,1000,stream);
        lenbuffer=strlen(buffer);
        lenbuffer2=strlen(buffer2);
        if(lenbuffer2+lenbuffer+2>sizebuffer2){
          sizebuffer2 = lenbuffer2+lenbuffer+2+1000;
          ResizeMemory((void **)&buffer2,(unsigned int)sizebuffer2);
        }
        strcat(buffer2,buffer);
        slashptr=strstr(buffer,"/");
        buffer3=buffer2;
      }
      start=STRSTR(buffer3+3,"ID");
      if(start!=NULL)nsurfids++;
      surf_id=NULL;
      surf_id2=NULL;
      for(c=start;*c!='\0';c++){
        if(surf_id==NULL&&*c=='\''){
          surf_id=c+1;
          continue;
        }
        if(surf_id!=NULL&&*c=='\''){
          *c='\0';
          NewMemory((void **)&surf_id2,strlen(surf_id)+1);
          strcpy(surf_id2,surf_id);
          surfids[nsurfids-1].label=surf_id2;
          surfids[nsurfids-1].location=1;
          surfids[nsurfids-1].show=1;
          break;
        }
      }

    }
  }

  /* identify duplicate surfaces */
  /*** debug: make sure ->show is defined for all cases ***/

  nsurfids_shown=0;
  for(i=0;i<nsurfids;i++){
    labeli = surfids[i].label;
    nexti = 0;
    for(j=0;j<nsurfaces;j++){
      surfj = surfaceinfo + j;
      labelj = surfj->surfacelabel;
      if(strcmp(labeli,labelj)==0){
        nexti = 1;
        break;
      }
    }
    if(nexti==1){
      surfids[i].show=0;
      continue;
    }
    for(j=0;j<i;j++){
      labelj = surfids[j].label;
      if(strcmp(labeli,labelj)==0){
        nexti = 1;
        break;
      }
    }
    if(nexti==1){
      surfids[i].show=0;
      continue;
    }
    nsurfids_shown++;

  }

  /* add surfaces found in database to those surfaces defined in previous SURF lines */

  if(nsurfids_shown>0){
    if(nsurfaces==0){
      FREEMEMORY(surfaceinfo);
      FREEMEMORY(textureinfo);
      NewMemory((void **)&surfaceinfo,nsurfids_shown*sizeof(surface));
      NewMemory((void **)&textureinfo,nsurfids_shown*sizeof(surface));
    }
    if(nsurfaces>0){
      if(surfaceinfo==NULL){
        NewMemory((void **)&surfaceinfo,(nsurfids_shown+nsurfaces)*sizeof(surface));
      }
      else{
        ResizeMemory((void **)&surfaceinfo,(nsurfids_shown+nsurfaces)*sizeof(surface));
      }
      if(textureinfo==NULL){
        NewMemory((void **)&textureinfo,(nsurfids_shown+nsurfaces)*sizeof(surface));
      }
      else{
        ResizeMemory((void **)&textureinfo,(nsurfids_shown+nsurfaces)*sizeof(surface));
      }
    }
    surfj = surfaceinfo + nsurfaces - 1;
    for(j=0;j<nsurfids;j++){
      if(surfids[j].show==0)continue;
      surfj++;
      initsurface(surfj);
      surfj->surfacelabel=surfids[j].label;
    }
    nsurfaces += nsurfids_shown;
  }
  update_sorted_surfidlist();
}

/* ------------------ surfid_compare ------------------------ */

int surfid_compare( const void *arg1, const void *arg2 ){
  surface *surfi, *surfj;

  int i, j;
  i=*(int *)arg1;
  j=*(int *)arg2;

  surfi = surfaceinfo + i;
  surfj = surfaceinfo + j;

  return(strcmp(surfi->surfacelabel,surfj->surfacelabel));
}

/* ------------------ updated_sorted_surfidlist ------------------------ */

void update_sorted_surfidlist(){
  int i;

  FREEMEMORY(sorted_surfidlist);
  FREEMEMORY(inv_sorted_surfidlist);
  NewMemory((void **)&sorted_surfidlist,nsurfaces*sizeof(int));
  NewMemory((void **)&inv_sorted_surfidlist,nsurfaces*sizeof(int));


  nsorted_surfidlist=nsurfaces;
  for(i=0;i<nsorted_surfidlist;i++){
    sorted_surfidlist[i]=i;
  }


  qsort( (int *)sorted_surfidlist, (size_t)nsurfaces, sizeof(int), surfid_compare );
  for(i=0;i<nsorted_surfidlist;i++){
    inv_sorted_surfidlist[sorted_surfidlist[i]]=i;
  }


}


        /* 
        new OBST format:
        i1 i2 j1 j2 k1 k2 colorindex blocktype       : if blocktype!=1&&colorindex!=-3
        i1 i2 j1 j2 k1 k2 colorindex blocktype r g b : if blocktype!=1&&colorindex==-3
        i1 i2 j1 j2 k1 k2 ceiling texture blocktype [wall texture floor texture] : if blocktype==1
        int colorindex, blocktype;
        colorindex: -1 default color
                    -2 invisible
                    -3 use r g b color
                    >=0 color/color2/texture index
        blocktype: 0 regular block non-textured
                   1 regular block textured
                   2 outline
                   3 smoothed block
        r g b           colors used if colorindex==-3
        */

/* ------------------ backup_blockage ------------------------ */

void backup_blockage(blockagedata *bc){
  int i;

  bc->xyzORIG[0]=bc->xmin;
  bc->xyzORIG[1]=bc->xmax;
  bc->xyzORIG[2]=bc->ymin;
  bc->xyzORIG[3]=bc->ymax;
  bc->xyzORIG[4]=bc->zmin;
  bc->xyzORIG[5]=bc->zmax;
  bc->ijkORIG[0]=bc->ijk[0];
  bc->ijkORIG[1]=bc->ijk[1];
  bc->ijkORIG[2]=bc->ijk[2];
  bc->ijkORIG[3]=bc->ijk[3];
  bc->ijkORIG[4]=bc->ijk[4];
  bc->ijkORIG[5]=bc->ijk[5];
  for(i=0;i<6;i++){
    bc->surfORIG[i]=bc->surf[i];
    bc->surf_indexORIG[i]=bc->surf_index[i];
  }
}

/* ------------------ ifsmoothblock ------------------------ */

int ifsmoothblock(void){
  int i,j;
  mesh *meshi;
  blockagedata *bc;

  for(i=0;i<nmeshes;i++){
    meshi = meshinfo + i;
    for(j=0;j<meshi->nbptrs;j++){
      bc = meshi->blockageinfoptrs[j];
      if(bc->type==BLOCK_smooth&&bc->del!=1)return 1;
    }
  }
  return 0;
}
/* ------------------ updateusetextures ------------------------ */

void updateusetextures(void){
  int i,j,k;
  mesh *meshi;
  blockagedata *bc;
  texture *texti;
  ventdata *vi;

  for(i=0;i<ntextures;i++){
    texti=textureinfo + i;
    texti->used=0;
  }
  for(i=0;i<nmeshes;i++){
    meshi=meshinfo + i;
    if(textureinfo!=NULL){
      for(j=0;j<meshi->nbptrs;j++){
        bc=meshi->blockageinfoptrs[j];
        for(k=0;k<6;k++){
          texti = bc->surf[k]->textureinfo;
          if(texti!=NULL&&texti->loaded==1){
            if(texti->loaded==1&&bc->type==BLOCK_smooth)bc->type=BLOCK_texture;
            if(usetextures==1)texti->display=1;
            texti->used=1;
          }
        }
      }
    }
    for(j=0;j<meshi->nvents+12;j++){
      vi = meshi->ventinfo + j;
      if(vi->surf[0]==NULL){
        if(vent_surfacedefault!=NULL){
          if(j>=meshi->nvents+6){
            vi->surf[0]=exterior_surfacedefault;
          }
          else{
            vi->surf[0]=vent_surfacedefault;
          }
        }
      }
      if(textureinfo!=NULL){
        if(vi->surf[0]!=NULL){
          texti = vi->surf[0]->textureinfo;
          if(texti!=NULL&&texti->loaded==1){
            if(usetextures==1)texti->display=1;
            texti->used=1;
          }
        }
      }
    }
  }
  for(i=0;i<ndevice_texture_list;i++){
    int texture_index;

    texture_index  = device_texture_list_index[i];
    texti=textureinfo + texture_index;
    if(texti!=NULL&&texti->loaded==1){
      if(usetextures==1)texti->display=1;
      texti->used=1;
    }
  }
  ntextures_loaded_used=0;
  for(i=0;i<ntextures;i++){
    texti = textureinfo + i;
    if(texti->loaded==0)continue;
    if(texti->used==0)continue;
    ntextures_loaded_used++;
  }
}

/* ------------------ initsurface ------------------------ */

void initsurface(surface *surf){
  surf->used_by_obst=0;
  surf->used_by_vent=0;
  surf->emis=1.0;
  surf->temp_ignition=TEMP_IGNITION_MAX;
  surf->surfacelabel=NULL;
  surf->texturefile=NULL;
  surf->textureinfo=NULL;
  surf->color=block_ambient2;
  surf->t_width=1.0;
  surf->t_height=1.0;
  surf->type=BLOCK_regular;
  surf->obst_surface=1;
  surf->location=0;
  surf->invisible=0;
  surf->transparent=0;
}

/* ------------------ initventsurface ------------------------ */

void initventsurface(surface *surf){
  surf->emis=1.0;
  surf->temp_ignition=TEMP_IGNITION_MAX;
  surf->surfacelabel=NULL;
  surf->texturefile=NULL;
  surf->textureinfo=NULL;
  surf->color=ventcolor;
  surf->t_width=1.0;
  surf->t_height=1.0;
  surf->type=BLOCK_outline;
  surf->obst_surface=0;

}

/* ------------------ setsurfaceindex ------------------------ */

void setsurfaceindex(blockagedata *bc){
  int i,j,jj;
  surface *surfj;
  char *surflabel, *bclabel;
  int *surf_index;
  int wall_flag;

  for(i=0;i<6;i++){
    bc->surf_index[i]=-1;
    bclabel=bc->surf[i]->surfacelabel;
    if(bc->surf[i]==NULL)continue;
    for(jj=1;jj<nsurfaces+1;jj++){
      j=jj;
      if(jj==nsurfaces)j=0;
      surfj = surfaceinfo + j;
      surflabel=surfj->surfacelabel;
      if(strcmp(bclabel,surflabel)!=0)continue;
      bc->surf_index[i]=j;
      break;
    }
  }
  surf_index = bc->surf_index;
  wall_flag = 1;
  for(i=1;i<6;i++){
    if(surf_index[i]!=surf_index[0]){
      wall_flag=0;
      break;
    }
  }
  if(wall_flag==1){
    bc->walltype=WALL_1;
    bc->walltypeORIG=WALL_1;
    return;
  }
  if(
      surf_index[UP_X]==surf_index[DOWN_X]&&
    surf_index[DOWN_Y]==surf_index[DOWN_X]&&
      surf_index[UP_Y]==surf_index[DOWN_X]
    ){
    bc->walltype=WALL_3;
    bc->walltypeORIG=WALL_3;
    return;
  }
  bc->walltype=WALL_6;
  bc->walltypeORIG=WALL_6;
  
}

/* ------------------ getnewfilename ------------------------ */

int getnewfilename(void){
  size_t len;
  int incnumber;
  char *c,*base,*index,*ext;
  char buffer[1025];
  char buffer2[1025];
  int numbers;

  if(fds_filein==NULL)return 2;
  len=strlen(fds_filein);
  if(len==0||len>1023)return 2;
  strcpy(buffer,fds_filein);
  base=buffer;

  ext=strrchr(fds_filein,'.');
  if(ext==NULL)return 2;
  base[ext-fds_filein]='\0';

  index=strrchr(base,'_');
  numbers=0;
  if(index!=NULL){
    numbers=1;
    for(c=index+1;*c!='\0';c++){
      if(isdigit(*c)==0){
        numbers=0;
        break;
      }
    }
  }
  if(numbers==1&&index!=NULL){
    base[index-base]='\0';
    index++;
    sscanf(index,"%i",&incnumber);
    sprintf(buffer2,"_%.3i",incnumber+1);
  }
  else{
    index=NULL;
    strcpy(buffer2,"_001");
  }
  len=strlen(base)+strlen(buffer2)+strlen(ext);
  FREEMEMORY(fds_fileout);
  FREEMEMORY(fds_fileout2);
  if(NewMemory((void **)&fds_fileout,(unsigned int)(len+1))==0||
    NewMemory((void **)&fds_fileout2,(unsigned int)(len+5))==0){
    FREEMEMORY(fds_fileout);
    FREEMEMORY(fds_fileout2);
    return 2;
  }
  STRCPY(fds_fileout,base);
  STRCAT(fds_fileout,buffer2);
  STRCAT(fds_fileout,ext);
  STRCPY(fds_fileout2,fds_fileout);
  STRCAT(fds_fileout2,"_chg");


/*
  len=strlen(fds_filein);
  FREEMEMORY(fds_fileout);
  FREEMEMORY(fds_fileout2);
  if(NewMemory((void **)&fds_fileout,(unsigned int)(len+5))==0||
    NewMemory((void **)&fds_fileout2,(unsigned int)(len+5))==0){
    FREEMEMORY(fds_fileout);
    FREEMEMORY(fds_fileout2);
    return 2;
  }
  STRCPY(fds_fileout,fds_filein);
  STRCAT(fds_fileout,"_new");
  STRCPY(fds_fileout2,fds_filein);
  STRCAT(fds_fileout2,"_chg");
  */
  return 0;
}

/* ------------------ initobst ------------------------ */

void initobst(blockagedata *bc, surface *surf,int index,int meshindex){
  int colorindex, blocktype;
  int i;
  char blocklabel[255];
  size_t len;

  bc->prop=NULL;
  bc->is_wuiblock=0;
  bc->transparent=0;
  bc->usecolorindex=0;
  bc->colorindex=-1;
  bc->nshowtime=0;
  bc->hole=0;
  bc->showtime=NULL;
  bc->showhide=NULL;
  bc->showtimelist=NULL;
  bc->show=1;
  bc->id=index;
  bc->meshindex=meshindex;
  bc->hidden=0;
  bc->invisible=0;
  bc->texture_origin[0]=texture_origin[0];
  bc->texture_origin[1]=texture_origin[1];
  bc->texture_origin[2]=texture_origin[2];

  /* 
  new OBST format:
  i1 i2 j1 j2 k1 k2 colorindex blocktype       : if blocktype!=1&&colorindex!=-3
  i1 i2 j1 j2 k1 k2 colorindex blocktype r g b : if blocktype!=1&&colorindex==-3
  i1 i2 j1 j2 k1 k2 ceiling texture blocktype [wall texture floor texture] : if blocktype==1
  int colorindex, blocktype;
  colorindex: -1 default color
              -2 invisible
              -3 use r g b color
              >=0 color/color2/texture index
  blocktype: 0 regular block non-textured
             1 regular block textured
             2 outline
             3 smoothed block
  r g b           colors used if colorindex==-3
  */
  colorindex=-1;
  blocktype=0;
  bc->ijk[IMIN]=0;
  bc->ijk[IMAX]=1;
  bc->ijk[JMIN]=0;
  bc->ijk[JMAX]=1;
  bc->ijk[KMIN]=0;
  bc->ijk[KMAX]=1;
  bc->colorindex=colorindex;
  bc->type=blocktype;

  bc->del=0;
  bc->changed=0;
  bc->changed_surface=0;
  bc->walltype=WALL_1;

  bc->color=surf->color;
  bc->useblockcolor=0;
  for(i=0;i<6;i++){
    bc->surf_index[i]=-1;
    bc->surf[i]=surf;
    bc->faceinfo[i]=NULL;
  }
  for(i=0;i<7;i++){
    bc->patchvis[i]=1;
  }
  sprintf(blocklabel,"**blockage %i",index);
  len=strlen(blocklabel);
  NewMemory((void **)&bc->label,((unsigned int)(len+1))*sizeof(char));
  strcpy(bc->label,blocklabel);
  
}

/* ------------------ initmesh ------------------------ */

void initmesh(mesh *meshi){

  meshi->mesh_offset_ptr=NULL;
#ifdef pp_CULL
  meshi->cullinfo=NULL;
  meshi->culldefined=0;
  meshi->cullQueryId=NULL;
  meshi->cull_smoke3d=NULL;
  meshi->smokedir_old=-100;
#endif
  meshi->blockvis=1;
  meshi->zcell=NULL;
  meshi->terrain=NULL;
  meshi->meshrgb[0]=0.0;
  meshi->meshrgb[1]=0.0;
  meshi->meshrgb[2]=0.0;
  meshi->meshrgb_ptr=NULL;
  meshi->showsmoothtimelist=NULL;
  meshi->cellsize=0.0;
  meshi->smokeloaded=0;
  meshi->hrrpuv_cutoff=600.0;
  meshi->smokedir=1;
  meshi->merge_alpha=NULL;
  meshi->merge_color=NULL;
  meshi->dx=1.0;
  meshi->dy=1.0;
  meshi->dz=1.0;
  meshi->dxy=1.0;
  meshi->dxz=1.0;
  meshi->dyz=1.0;
  meshi->label=NULL;
  meshi->mxpatch_frames=0;
  meshi->visx=0;
  meshi->visy=1;
  meshi->visz=0;
  meshi->visx2=0;
  meshi->visy2=1;
  meshi->visz2=0;
  meshi->slicedir=2;
  meshi->visInteriorPatches=0;
  meshi->plot3dfilenum=-1;
  meshi->patchfilenum=-1;
  meshi->obst_bysize=NULL;
  meshi->iqdata=NULL;      meshi->qdata=NULL;
  meshi->yzcolorbase=NULL; meshi->xzcolorbase=NULL; meshi->xycolorbase=NULL; 
  meshi->yzcolorfbase=NULL;meshi->xzcolorfbase=NULL;meshi->xycolorfbase=NULL;
  meshi->yzcolortbase=NULL;meshi->xzcolortbase=NULL;meshi->xycolortbase=NULL;
  meshi->dx_xy=NULL;       meshi->dy_xy=NULL;       meshi->dz_xy=NULL;
  meshi->dx_xz=NULL;       meshi->dy_xz=NULL;       meshi->dz_xz=NULL;
  meshi->dx_yz=NULL;       meshi->dy_yz=NULL;       meshi->dz_yz=NULL;
  meshi->c_iblank_xy=NULL; meshi->c_iblank_xz=NULL; meshi->c_iblank_yz=NULL;
  meshi->iblank_smoke3d=NULL;
  meshi->animatedsurfaces=NULL;
  meshi->blockagesurface=NULL;
  meshi->blockagesurfaces=NULL;
  meshi->smoothblockagecolors=NULL;
  meshi->showlevels=NULL;
  meshi->isolevels=NULL;
  meshi->nisolevels=0;
  meshi->isotimes=NULL;
  meshi->isotimeslist=NULL;
  meshi->isofilenum=-1;
  meshi->nvents=0;
  meshi->ndummyvents=0;
  meshi->npatches=0;
  meshi->patchfacevis2=0;
  meshi->patchtype=NULL;
  meshi->offset[0]=0.0;
  meshi->offset[1]=0.0;
  meshi->offset[2]=0.0;
  meshi->patchtype=NULL;
//  meshi->patchfacevis=NULL;
  meshi->patchdir=NULL;
  meshi->patch_surfindex=NULL;
  meshi->pi1=NULL; meshi->pi2=NULL;
  meshi->pj1=NULL; meshi->pj2=NULL; 
  meshi->pk1=NULL; meshi->pk2=NULL;
  meshi->meshonpatch=NULL;
  meshi->blockonpatch=NULL;
  meshi->ptype=NULL;
  meshi->patchrow=NULL, meshi->patchcol=NULL, meshi->blockstart=NULL;
  meshi->zipoffset=NULL, meshi->zipsize=NULL;
  meshi->visPatches=NULL;
  meshi->xyzpatch=NULL;
  meshi->xyzpatch_threshold=NULL;
  meshi->ipqq=NULL, meshi->ipqqi=NULL;
  meshi->ipqqi_zlib=NULL, meshi->ipqq_zlib=NULL;
  meshi->patchtimes=NULL, meshi->pqq=NULL, meshi->pqqi=NULL;
  meshi->thresholdtime=NULL;
  meshi->patchblank=NULL;
  meshi->patch_contours=NULL;
  meshi->patchtimeslist=NULL;
  meshi->ntc=0;
  meshi->nspr=0;
  meshi->xsprplot=NULL;
  meshi->ysprplot=NULL;
  meshi->zsprplot=NULL;
  meshi->xspr=NULL;
  meshi->yspr=NULL;
  meshi->zspr=NULL;
  meshi->tspr=NULL;
  meshi->nheat=0;
  meshi->xheatplot=NULL;
  meshi->yheatplot=NULL;
  meshi->zheatplot=NULL;
  meshi->xheat=NULL;
  meshi->yheat=NULL;
  meshi->zheat=NULL;
  meshi->theat=NULL;
  meshi->blockageinfoptrs=NULL;
  meshi->deletelist=NULL;
  meshi->carveblockptrs=NULL;
  meshi->ncarveblocks=0;
  meshi->ndeletelist=0;
  meshi->nsmoothblockagecolors=0;

  meshi->surface_tempmax=SURFACE_TEMPMAX;
  meshi->surface_tempmin=SURFACE_TEMPMIN;
 
  meshi->faceinfo=NULL;
  meshi->face_normals_single=NULL;
  meshi->face_normals_double=NULL;
  meshi->face_transparent_double=NULL;
  meshi->face_textures=NULL;
  meshi->face_outlines=NULL;

  meshi->xplt=NULL;
  meshi->yplt=NULL;
  meshi->zplt=NULL;
  meshi->xplt_orig=NULL;
  meshi->yplt_orig=NULL;
  meshi->zplt_orig=NULL;
  meshi->c_iblank_cell=NULL;
  meshi->c_iblank_x=NULL;
  meshi->c_iblank_y=NULL;
  meshi->c_iblank_z=NULL;

  meshi->xbar=1.0;
  meshi->xbar0=0.0;
  meshi->ybar=1.0;
  meshi->ybar0=0.0;
  meshi->zbar=1.0;
  meshi->zbar0=0.0;

  meshi->xyzmaxdiff=1.0;

  meshi->xcen=0.5;
  meshi->ycen=0.5;
  meshi->zcen=0.5;

  meshi->plotx=1;
  meshi->ploty=1;
  meshi->plotz=1;

  meshi->boxoffset=0.0;
  meshi->c_iblank=NULL;
  meshi->ventinfo=NULL;
  meshi->select_min=0;
  meshi->select_max=0;
}

/* ------------------ freelabels ------------------------ */

void freelabels(flowlabels *flowlabel){
  FREEMEMORY(flowlabel->longlabel);
  FREEMEMORY(flowlabel->shortlabel);
  FREEMEMORY(flowlabel->unit);
}

/* ------------------ createnulllaels ------------------------ */

int createnulllabel(flowlabels *flowlabel){
  char buffer[255];
  size_t len;

  len=strlen("Particles");
  strcpy(buffer,"Particles");
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void **)&flowlabel->longlabel,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->longlabel,buffer);


  len=strlen("Particles");
  strcpy(buffer,"Particles");
  len=strlen(buffer);
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void **)&flowlabel->shortlabel,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->shortlabel,buffer);

  len=strlen("Particles");
  strcpy(buffer,"Particles");
  len=strlen(buffer);
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void *)&flowlabel->unit,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->unit,buffer);
  return 0;
}

/* ------------------ readlabels ------------------------ */

int readlabels(flowlabels *flowlabel, FILE *stream){
  char buffer2[255], *buffer;
  size_t len;

  if(fgets(buffer2,255,stream)==NULL){
    strcpy(buffer2,"*");
  }

  len=strlen(buffer2);
  buffer2[len-1]='\0';
  buffer=trim_front(buffer2);
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void **)&flowlabel->longlabel,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->longlabel,buffer);


  if(fgets(buffer2,255,stream)==NULL){
    strcpy(buffer2,"**");
  }

  len=strlen(buffer2);
  buffer2[len-1]='\0';
  buffer=trim_front(buffer2);
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void **)&flowlabel->shortlabel,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->shortlabel,buffer);

  if(fgets(buffer2,255,stream)==NULL){
    strcpy(buffer2,"***");
  }

  len=strlen(buffer2);
  buffer2[len-1]='\0';
  buffer=trim_front(buffer2);
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void *)&flowlabel->unit,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->unit,buffer);
  return 0;
}

/* ------------------ readlabels_cellcenter ------------------------ */

int readlabels_cellcenter(flowlabels *flowlabel, FILE *stream){
  char buffer2[255], *buffer;
  size_t len;

  if(fgets(buffer2,255,stream)==NULL){
    strcpy(buffer2,"*");
  }

  len=strlen(buffer2);
  buffer2[len-1]='\0';
  buffer=trim_front(buffer2);
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void **)&flowlabel->longlabel,(unsigned int)(len+1+15))==0)return 2;
  STRCPY(flowlabel->longlabel,buffer);
  STRCAT(flowlabel->longlabel,"(cell centered)");

  if(fgets(buffer2,255,stream)==NULL){
    strcpy(buffer2,"**");
  }

  len=strlen(buffer2);
  buffer2[len-1]='\0';
  buffer=trim_front(buffer2);
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void **)&flowlabel->shortlabel,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->shortlabel,buffer);

  if(fgets(buffer2,255,stream)==NULL){
    strcpy(buffer2,"***");
  }

  len=strlen(buffer2);
  buffer2[len-1]='\0';
  buffer=trim_front(buffer2);
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void *)&flowlabel->unit,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->unit,buffer);
  return 0;
}

/* ------------------ readlabels_terrain ------------------------ */

int readlabels_terrain(flowlabels *flowlabel, FILE *stream){
  char buffer2[255],*buffer;
  size_t len;

  if(fgets(buffer2,255,stream)==NULL){
    strcpy(buffer2,"*");
  }

  len=strlen(buffer2);
  buffer2[len-1]='\0';
  buffer=trim_front(buffer2);
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void **)&flowlabel->longlabel,(unsigned int)(len+1+9))==0)return 2;
  STRCPY(flowlabel->longlabel,buffer);
  STRCAT(flowlabel->longlabel,"(terrain)");

  if(fgets(buffer2,255,stream)==NULL){
    strcpy(buffer2,"**");
  }

  len=strlen(buffer2);
  buffer2[len-1]='\0';
  buffer=trim_front(buffer2);
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void **)&flowlabel->shortlabel,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->shortlabel,buffer);

  if(fgets(buffer2,255,stream)==NULL){
    strcpy(buffer2,"***");
  }

  len=strlen(buffer2);
  buffer2[len-1]='\0';
  buffer=trim_front(buffer2);
  trim(buffer);
  len=strlen(buffer);
  if(NewMemory((void *)&flowlabel->unit,(unsigned int)(len+1))==0)return 2;
  STRCPY(flowlabel->unit,buffer);
  return 0;
}

/* ------------------ match_upper ------------------------ */

int match_upper(char *buffer, const char *key, unsigned int lenkey){
  size_t lenbuffer;
  int i;

  ASSERT(strlen(key)==lenkey);
  trim(buffer);
  lenbuffer=strlen(buffer);

  if(lenbuffer<lenkey)return 0;
  for(i=0;i<lenkey;i++){
    if(toupper(buffer[i])!=toupper(key[i]))return 0;
  }
  if(lenbuffer>lenkey&&buffer[lenkey]==':')return 2;
  if(lenbuffer>lenkey&&buffer[lenkey]!=' ')return 0;
  return 1;
}

/* ------------------ match ------------------------ */

int match(char *buffer, const char *key, unsigned int lenkey){
  size_t lenbuffer;

  ASSERT(strlen(key)==lenkey);
//  if(strncmp(buffer,key,lenkey) == 0)return(1);
//  return(0);
    if(strncmp(buffer,key,lenkey) != 0)return(0);
    trim(buffer);
    lenbuffer=strlen(buffer);
    if(lenbuffer==lenkey)return 1;
    if(lenbuffer>lenkey){
      if(buffer[lenkey]==' ')return 1;
    }
    return 0;
}

/* ------------------ matchonly ------------------------ */

int matchonly(char *buffer, const char *key, unsigned int lenkey){
  size_t lenbuffer;

  ASSERT(strlen(key)==lenkey);
  
  if(strncmp(buffer,key,lenkey) != 0)return 0;
  trim(buffer);
  lenbuffer=strlen(buffer);
  if(lenbuffer==lenkey){
    return 1;
  }
  else{
    return 0;
  }
}
/* ------------------ readini ------------------------ */

int readini(int scriptconfigfile){
  STRUCTSTAT statbuff1, statbuff2, statbuff3, statbuff4;
  int statfile1=-1, statfile2=-1, statfile3=-1, statfile4=-1;
  char smvprogini[1024];
  char *smvprogini_ptr=NULL;

  nlabels=nlabelssmv;
  nticks=ntickssmv;
  strcpy(smvprogini,"");
  if(smvprogdir!=NULL)strcat(smvprogini,smvprogdir);
  strcat(smvprogini,"smokeview.ini");
  if(smokeviewini!=NULL){
#ifdef WIN32
    if(STRCMP(smvprogini,smokeviewini)!=0)smvprogini_ptr=smvprogini;
#else
    if(strcmp(smvprogini,smokeviewini)!=0)smvprogini_ptr=smvprogini;
#endif
  }
  
  // check if configuration files exist

  if(smokeviewini!=NULL)statfile1=STAT(smokeviewini,&statbuff1);
  if(smvprogini_ptr!=NULL){
    statfile2=STAT(smvprogini_ptr,&statbuff2);
    if(statfile2!=0)smvprogini_ptr=NULL;
  }
  if(INIfile!=NULL)statfile3=STAT(INIfile,&statbuff3);
  if(caseinifilename!=NULL)statfile4=STAT(caseinifilename,&statbuff4);

  // check if config files read in earlier were modifed later

  if(statfile1==0&&statfile2==0&&statbuff1.st_mtime>statbuff2.st_mtime){
    printf("*** warning: The initialization file, %s,\n is newer than %s \n",smokeviewini,smvprogini_ptr);
  }
  if(statfile1==0&&statfile3==0&&statbuff1.st_mtime>statbuff3.st_mtime){
    printf("*** warning: The initialization file, %s, is newer than %s \n",smokeviewini,INIfile);
  }
  if(statfile1==0&&statfile4==0&&statbuff1.st_mtime>statbuff4.st_mtime){
    printf("*** warning: The initialization file, %s, is newer than %s \n",smokeviewini,caseinifilename);
  }

  if(statfile2==0&&statfile3==0&&statbuff2.st_mtime>statbuff3.st_mtime){
    printf("*** warning: The initialization file, %s, is newer than %s \n",smvprogini_ptr,INIfile);
  }
  if(statfile2==0&&statfile4==0&&statbuff2.st_mtime>statbuff4.st_mtime){
    printf("*** warning: The initialization file, %s, is newer than %s \n",smvprogini_ptr,caseinifilename);
  }

  if(statfile3==0&&statfile4==0&&statbuff3.st_mtime>statbuff4.st_mtime){
    printf("*** warning: The initialization file, %s, is newer than  %s \n",INIfile,caseinifilename);
  }

  // read in config files if they exist

  // smokeview.ini ini in install directory

  if(statfile1==0&&smokeviewini!=NULL){
    if(readini2(smokeviewini,0)==2)return 2;
  }

  // smokeview.ini in smokeview directory (could be different than above)

  if(statfile2==0&&smvprogini_ptr!=NULL){
    if(readini2(smvprogini_ptr,0)==2)return 2;
  }

  // smokeview.ini in case directory

  if(statfile3==0&&INIfile!=NULL){
    if(readini2(INIfile,0)==2)return 2;
  }

  // read in casename.ini

  if(statfile4==0&&caseinifilename!=NULL){
    if(readini2(caseinifilename,1)==2)return 2;
  }

  // read in ini file specified in script

  if(scriptinifilename2!=NULL&&scriptconfigfile==2){
    if(readini2(scriptinifilename2,1)==2)return 2;
    updatecolors(-1);
    updateshowtitles();
    scriptinifilename2=NULL;
  }
  updateglui();
  if(showall_textures==1)TextureShowMenu(-1);
  if(ncolorbars<=ndefaultcolorbars){
    initdefaultcolorbars();
  }
  updatezoommenu=1;
  return 0;
}

/* ------------------ readini2 ------------------------ */

int readini2(char *inifile, int localfile){
  char buffer[255],buffer2[255];

  int tours_flag;
  int mxframes_ini=0;
  int mxpoints_ini=0;
  int nn;
  int n3d;
  int iplot3d, isetmin, isetmax;
  float p3mintemp, p3maxtemp;
  float *scale;
  float valmin, valmax;
  int setvalmin, setvalmax;
  FILE *stream;
  int i;
  int j;
  int tempval;

  updatemenu=1;
  updatefacelists=1;

  if( (stream=fopen(inifile,"r"))==NULL)return 1;

  for(i=0;i<nunitclasses_ini;i++){
    f_units *uc;

    uc=unitclasses_ini+i;
    FREEMEMORY(uc->units);
  }
  FREEMEMORY(unitclasses_ini);
  nunitclasses_ini=0;

  if(localfile==1){
    update_inilist();
  }

  printf("reading: %s\n",inifile);
  if(localfile==1){
    update_selectedtour_index=0;
  }

  /* find number of each kind of file */

  while(!feof(stream)){
    CheckMemory;
    if(fgets(buffer,255,stream)==NULL)break;

    if(match(buffer,"MESHOFFSET",10)==1){
      int meshnum;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&meshnum);
      if(meshnum>=0&&meshnum<nmeshes){
        mesh *meshi;

        meshi = meshinfo + meshnum;
        meshi->mesh_offset_ptr=meshi->mesh_offset;
      }
    }
    if(match(buffer,"MESHVIS",7)==1){
      int nm;
      mesh *meshi;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&nm);
      for(i=0;i<nm;i++){
        int vis;

        if(i>nmeshes-1)break;
        meshi = meshinfo + i;
        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&meshi->blockvis);
        if(meshi->blockvis!=0)meshi->blockvis=1;
      }
      continue;
    }
    if(match(buffer,"SPHERESEGS",10)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&device_sphere_segments);
      if(device_sphere_segments<3)device_sphere_segments=3;
      initspheresegs(device_sphere_segments,2*device_sphere_segments);
      continue;
    }
    if(match(buffer,"OFFSETSLICE",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&offset_slice);
      if(offset_slice!=0)offset_slice=1;
      continue;
    }

    if(match(buffer,"VECLENGTH",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&iveclengths);
      vecfactor = get_vecfactor(&iveclengths);
      continue;
    }

    if(match(buffer,"ISOTRAN2",8)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&transparent_state);
      continue;
    }

    if(match(buffer,"SHOWSTREAK",10)==1){
      void ParticleStreakShowMenu(int var);

      fgets(buffer,255,stream);
      sscanf(buffer,"%i %i %i %i",&streak5show,&streak5step,&showstreakhead,&streak_index);
      if(streak5show!=1){
        streak5show=0;
        streak_index=-2;
      }
      if(showstreakhead!=1)showstreakhead=0;
      update_streaks=1;
      continue;
    }

    if(match(buffer,"SHOWTERRAIN",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visTerrainType);
      /*
      if(visTerrain!=0)visTerrain=1;
      if(visTerrain==1){
        visTerrain=1-visTerrain;
        GeometryMenu(17);
      }
      */
      continue;
    }
    if(match(buffer,"STEREO",6)==1){
      fgets(buffer,255,stream);
      showstereoOLD=showstereo;
      sscanf(buffer,"%i",&showstereo);
      if(showstereo<0)showstereo=0;
      if(showstereo>5)showstereo=5;
      if(showstereo==1&&videoSTEREO!=1)showstereo=0;
      update_glui_stereo();
      continue;
    }
    if(match(buffer,"TERRAINPARMS",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %i %i",terrain_rgba_zmin,terrain_rgba_zmin+1,terrain_rgba_zmin+2);

      fgets(buffer,255,stream);
      sscanf(buffer,"%i %i %i",terrain_rgba_zmax,terrain_rgba_zmax+1,terrain_rgba_zmax+2);

      fgets(buffer,255,stream);
      sscanf(buffer,"%f",&vertical_factor);

      for(i=0;i<3;i++){
        if(terrain_rgba_zmin[i]<0)terrain_rgba_zmin[i]=0;
        if(terrain_rgba_zmin[i]>255)terrain_rgba_zmin[i]=255;
        if(terrain_rgba_zmax[i]<0)terrain_rgba_zmax[i]=0;
        if(terrain_rgba_zmax[i]>255)terrain_rgba_zmax[i]=255;
      }
      if(vertical_factor<0.25)vertical_factor=0.25;
      if(vertical_factor>4.0)vertical_factor=4.0;
      update_terrain(0,vertical_factor);
      update_terrain_colors();
    
      continue;
    }
    if(match(buffer,"SMOKESENSORS",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %i",&show_smokesensors,&test_smokesensors);
      continue;
    }
    if(match(buffer,"SBATSTART",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&sb_atstart);
      if(sb_atstart!=0)sb_atstart=1;
      continue;
    }
#ifdef pp_GPU
    if(match(buffer,"USEGPU",6)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&usegpu);
      if(usegpu!=0)usegpu=1;
      continue;
    }
#endif
    if(match(buffer,"V_PLOT3D",8)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&tempval);
      if(tempval<0)tempval=0;
      n3d=tempval;
      if(n3d>mxplot3dvars)n3d=mxplot3dvars;
      for(i=0;i<n3d;i++){  
        fgets(buffer,255,stream);
        sscanf(buffer,"%i %i %f %i %f",&iplot3d,&isetmin,&p3mintemp,&isetmax,&p3maxtemp);
        iplot3d--;
        if(iplot3d>=0&&iplot3d<mxplot3dvars){
          setp3min[iplot3d]=isetmin;
          setp3max[iplot3d]=isetmax;
          p3min[iplot3d]=p3mintemp;
          p3max[iplot3d]=p3maxtemp;
        }
      }
      continue;
    }
    if(match(buffer,"UNLOAD_QDATA",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&unload_qdata);
      if(unload_qdata!=1)unload_qdata=0;
      continue;
    }
    if(match(buffer,"TREECOLORS",10)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",trunccolor,trunccolor+1,trunccolor+2);
      sscanf(buffer,"%f %f %f",treecolor,treecolor+1,treecolor+2);
      sscanf(buffer,"%f %f %f",treecharcolor,treecharcolor+1,treecharcolor+2);
      for(i=0;i<4;i++){
        if(treecolor[i]<0.0)treecolor[i]=0.0;
        if(treecharcolor[i]<0.0)treecharcolor[i]=0.0;
        if(trunccolor[i]<0.0)trunccolor[i]=0.0;
        if(treecolor[i]>1.0)treecolor[i]=1.0;
        if(treecharcolor[i]>1.0)treecharcolor[i]=1.0;
        if(trunccolor[i]>1.0)trunccolor[i]=1.0;
      }
      treecolor[3]=1.0;
      treecharcolor[3]=1.0;
      trunccolor[3]=1.0;
      continue;
    }
    if(match(buffer,"TRAINERVIEW",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&trainerview);
      if(trainerview!=2&&trainerview!=3)trainerview=1;
      continue;
    }
#ifdef pp_LIGHT
    if(match(buffer,"SMOKELIGHTING",13)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&show_smokelighting);
      if(show_smokelighting!=0)show_smokelighting=1;
      continue;
    }
#endif
    if(match(buffer,"SMOOTHBLOCKSOLID",16)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&smooth_block_solid);
      if(smooth_block_solid!=1)smooth_block_solid=0;
      continue;
    }
    if(match(buffer,"SHOWTRANSPARENTVENTS",20)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&show_transparent_vents);
      if(show_transparent_vents!=0)show_transparent_vents=1;
      continue;
    }
    if(match(buffer,"COLORBARTYPE",12)==1){
      colorbardata *cbi;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&colorbartype);
      continue;
    }
    if(match(buffer,"SHOWEXTREMEDATA",15)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&show_extremedata);
      if(show_extremedata!=1)show_extremedata=0;
      continue;
    }
    if(match(buffer,"EXTREMECOLORS",13) == 1){
      int mmin[3],mmax[3];

      fgets(buffer,255,stream);
      sscanf(buffer,"%i %i %i %i %i %i",
        mmin,mmin+1,mmin+2,
        mmax,mmax+1,mmax+2);
      for(i=0;i<3;i++){
        rgb_below_min[i]=mmin[i];
        rgb_above_max[i]=mmax[i];
      }
    }
    if(match(buffer,"SLICEAVERAGE",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i %i",&slice_average_flag,&slice_average_interval,&vis_slice_average,&slice_turbprop_flag);
      if(slice_average_flag!=1)slice_average_flag=0;
      if(slice_turbprop_flag!=1)slice_turbprop_flag=0;
      if(slice_average_flag==1)slice_turbprop_flag=0;
      if(slice_average_interval<0.0)slice_average_interval=0.0;
      continue;
    }
#ifdef pp_RENDER
    if(match(buffer,"SKYBOX",6)==1){
      skyboxdata *skyi;

      free_skybox();
      nskyboxinfo=1;
      NewMemory((void **)&skyboxinfo,nskyboxinfo*sizeof(skyboxdata));
      skyi = skyboxinfo;

      for(i=0;i<6;i++){
        fgets(buffer,255,stream);
        loadskytexture(buffer,skyi->face+i);
      }
    }
#endif
    if(match(buffer,"C_PLOT3D",8)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&tempval);
      if(tempval<0)tempval=0;
      n3d=tempval;
      if(n3d>mxplot3dvars)n3d=mxplot3dvars;
      for(i=0;i<n3d;i++){  
        fgets(buffer,255,stream);
        sscanf(buffer,"%i %i %f %i %f",&iplot3d,&isetmin,&p3mintemp,&isetmax,&p3maxtemp);
        iplot3d--;
        if(iplot3d>=0&&iplot3d<mxplot3dvars){
          setp3chopmin[iplot3d]=isetmin;
          setp3chopmax[iplot3d]=isetmax;
          p3chopmin[iplot3d]=p3mintemp;
          p3chopmax[iplot3d]=p3maxtemp;
        }
      }
      continue;
    }
    if(match(buffer,"DEVICENORMLENGTH",16)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f",&devicenorm_length);
      if(devicenorm_length<0.0||devicenorm_length>1.0)devicenorm_length=0.1;
      continue;
    }

    if(match(buffer,"SHOWHRRCUTOFF",13)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&show_hrrcutoff);
      if(show_hrrcutoff!=1)show_hrrcutoff=0;
      continue;
    }
    if(match(buffer,"TWOSIDEDVENTS",13)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %i",&show_bothsides_int,&show_bothsides_ext);
      if(show_bothsides_int!=1)show_bothsides_int=0;
      if(show_bothsides_ext!=1)show_bothsides_ext=0;
      continue;
    }
    if(match(buffer,"SHOWSLICEINOBST",15)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&show_slice_in_obst);
      if(show_slice_in_obst!=1)show_slice_in_obst=0;
      continue;
    }
    if(match(buffer,"SKIPEMBEDSLICE",14)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&skip_slice_in_embedded_mesh);
      if(skip_slice_in_embedded_mesh!=1)skip_slice_in_embedded_mesh=0;
      continue;
    }
    if(match(buffer,"CELLCENTERINTERP",16)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&cellcenter_interp);
      if(cellcenter_interp!=1)cellcenter_interp=0;
      continue;
    }
    if(match(buffer,"PERCENTILELEVEL",15)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f",&percentile_level);
      if(percentile_level<0.0)percentile_level=0.01;
      if(percentile_level>0.5)percentile_level=0.01;
      continue;
    }
    if(match(buffer,"TRAINERMODE",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&trainer_mode);
      continue;
    }
    if(match(buffer,"COMPRESSAUTO",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&compress_autoloaded);
      continue;
    }
    if(match(buffer,"PLOT3DAUTO",10)==1){
      int n3dsmokes=0;
      int seq_id;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&n3dsmokes);
      for(i=0;i<n3dsmokes;i++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&seq_id);
        get_startup_plot3d(seq_id);
      }
      update_load_startup=1;
      continue;
    }
    if(match(buffer,"VSLICEAUTO",10)==1){
      int n3dsmokes=0;
      int seq_id;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&n3dsmokes);
      for(i=0;i<n3dsmokes;i++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&seq_id);
        get_startup_vslice(seq_id);
      }
      update_load_startup=1;
      continue;
    }
    if(match(buffer,"SLICEAUTO",9)==1){
      int n3dsmokes=0;
      int seq_id;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&n3dsmokes);
      for(i=0;i<n3dsmokes;i++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&seq_id);
        get_startup_slice(seq_id);
      }
      update_load_startup=1;
      continue;
    }
    if(match(buffer,"PARTAUTO",8)==1){
      int n3dsmokes=0;
      int seq_id;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&n3dsmokes);
      for(i=0;i<n3dsmokes;i++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&seq_id);
        get_startup_part(seq_id);
      }
      update_load_startup=1;
      continue;
    }
    if(match(buffer,"ISOAUTO",7)==1){
      int n3dsmokes=0;
      int seq_id;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&n3dsmokes);
      for(i=0;i<n3dsmokes;i++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&seq_id);
        get_startup_iso(seq_id);
      }
      update_load_startup=1;
      continue;
    }
    if(match(buffer,"S3DAUTO",7)==1){
      int n3dsmokes=0;
      int seq_id;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&n3dsmokes);
      for(i=0;i<n3dsmokes;i++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&seq_id);
        get_startup_smoke(seq_id);
      }
      update_load_startup=1;
      continue;
    }
    if(match(buffer,"PATCHAUTO",9)==1){
      int n3dsmokes=0;
      int seq_id;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&n3dsmokes);
      for(i=0;i<n3dsmokes;i++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&seq_id);
        get_startup_patch(seq_id);
      }
      update_load_startup=1;
      continue;
    }
    if(match(buffer,"LOADFILESATSTARTUP",18)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&loadfiles_at_startup);
      continue;
    }
    if(match(buffer,"SHOWALLTEXTURES",15)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&showall_textures);
      continue;
    }
    if(match(buffer,"PIXELSKIP",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&pixel_skip);
      continue;
    }
    if(match(buffer,"PROJECTION",10)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&projection_type);
      update_projection_type();
      continue;
    }
    if(match(buffer,"T_PARTICLES",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i %f",&settmin_p,&tmin_p,&settmax_p,&tmax_p);
      continue;
    }
    if(match(buffer,"T_SLICE",7)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i %f",&settmin_s,&tmin_s,&settmax_s,&tmax_s);
      continue;
    }
    if(match(buffer,"T_ISO",5)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i %f",&settmin_i,&tmin_i,&settmax_i,&tmax_i);
      continue;
    }
    if(match(buffer,"T_BOUNDARY",10)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i %f",&settmin_b,&tmin_b,&settmax_b,&tmax_b);
      continue;
    }
    if(match(buffer,"V_PARTICLES",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i %f",&setpartmin,&partmin,&setpartmax,&partmax);
      continue;
    }
    if(match(buffer,"V5_PARTICLES",12)==1){
      int ivmin, ivmax;
      float vmin, vmax;
      char short_label[256],*s1;

      strcpy(short_label,"");
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i %f %s",&ivmin,&vmin,&ivmax,&vmax,short_label);
  
      if(npart5prop>0){
        int label_index;
        trim(short_label);
        s1=trim_front(short_label);
        label_index=get_part5prop_index_s(s1);
        if(label_index>=0&&label_index<npart5prop){
          part5prop *propi;

          propi = part5propinfo + label_index;
          propi->setvalmin=ivmin;
          propi->setvalmax=ivmax;
          propi->valmin=vmin;
          propi->valmax=vmax;
          switch (ivmin){
            case PERCENTILE_MIN:
              propi->percentile_min=vmin;
              break;
            case GLOBAL_MIN:
              propi->global_min=vmin;
              break;
            case SET_MIN:
              propi->user_min=vmin;
              break;
          }
          switch (ivmax){
            case PERCENTILE_MAX:
              propi->percentile_max=vmax;
              break;
            case GLOBAL_MAX:
              propi->global_max=vmax;
              break;
            case SET_MAX:
              propi->user_max=vmax;
              break;
          }
        }
      }
      continue;
    }
    if(match(buffer,"C_PARTICLES",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i %f",&setpartchopmin,&partchopmin,&setpartchopmax,&partchopmax);
      continue;
    }
    if(match(buffer,"V_SLICE",7)==1){
#ifdef pp_SLICECONTOURS
      char *level_val;
#endif

      fgets(buffer,255,stream);
      strcpy(buffer2,"");
      sscanf(buffer,"%i %f %i %f %s",&setvalmin,&valmin,&setvalmax,&valmax,buffer2);
#ifdef pp_SLICECONTOURS
      {
        char *colen;

        colen=strstr(buffer,":");
        level_val=NULL;
        if(colen!=NULL){
          level_val=colen+1;
          trim(level_val);
          *colen=0;
          if(strlen(level_val)>1){
            sscanf(level_val,"%f %f %i",&slice_line_contour_min,&slice_line_contour_max,&slice_line_contour_num);
          }
          {
            level_val=NULL;
          }
        }
      }
#endif
      if(strcmp(buffer2,"")!=0){
        for(i=0;i<nslice2;i++){
          if(strcmp(slicebounds[i].datalabel,buffer2)!=0)continue;
          slicebounds[i].setvalmin=setvalmin;
          slicebounds[i].setvalmax=setvalmax;
          slicebounds[i].valmin=valmin;
          slicebounds[i].valmax=valmax;
#ifdef pp_SLICECONTOURS
          if(level_val!=NULL){
            slicebounds[i].line_contour_min=slice_line_contour_min;  
            slicebounds[i].line_contour_max=slice_line_contour_max;  
            slicebounds[i].line_contour_num=slice_line_contour_num;  
          }
#endif
          break;
        }
      }
      else{
        for(i=0;i<nslice2;i++){
          slicebounds[i].setvalmin=setvalmin;
          slicebounds[i].setvalmax=setvalmax;
          slicebounds[i].valmin=valmin;
          slicebounds[i].valmax=valmax;
#ifdef pp_SLICECONTOURS
          slicebounds[i].line_contour_min=slice_line_contour_min;  
          slicebounds[i].line_contour_max=slice_line_contour_max;  
          slicebounds[i].line_contour_num=slice_line_contour_num;  
#endif
        }
      }
      continue;
    }
    if(match(buffer,"C_SLICE",7)==1){
      fgets(buffer,255,stream);
      strcpy(buffer2,"");
      sscanf(buffer,"%i %f %i %f %s",&setvalmin,&valmin,&setvalmax,&valmax,buffer2);
      if(strcmp(buffer,"")!=0){
        for(i=0;i<nslice2;i++){
          if(strcmp(slicebounds[i].datalabel,buffer2)!=0)continue;
          slicebounds[i].setchopmin=setvalmin;
          slicebounds[i].setchopmax=setvalmax;
          slicebounds[i].chopmin=valmin;
          slicebounds[i].chopmax=valmax;
          break;
        }
      }
      else{
        for(i=0;i<nslice2;i++){
          slicebounds[i].setchopmin=setvalmin;
          slicebounds[i].setchopmax=setvalmax;
          slicebounds[i].chopmin=valmin;
          slicebounds[i].chopmax=valmax;
        }
      }
      continue;
    }
    if(match(buffer,"V_ISO",5)==1){
      fgets(buffer,255,stream);
      strcpy(buffer2,"");
      sscanf(buffer,"%i %f %i %f %s",&setvalmin,&valmin,&setvalmax,&valmax,buffer2);
      if(strcmp(buffer2,"")!=0){
        for(i=0;i<niso_bounds;i++){
          if(strcmp(isobounds[i].datalabel,buffer2)!=0)continue;
          isobounds[i].setvalmin=setvalmin;
          isobounds[i].setvalmax=setvalmax;
          isobounds[i].valmin=valmin;
          isobounds[i].valmax=valmax;
          break;
        }
      }
      else{
        for(i=0;i<niso_bounds;i++){
          isobounds[i].setvalmin=setvalmin;
          isobounds[i].setvalmax=setvalmax;
          isobounds[i].valmin=valmin;
          isobounds[i].valmax=valmax;
        }
      }
      continue;
    }
    if(match(buffer,"C_ISO",5)==1){
      fgets(buffer,255,stream);
      strcpy(buffer2,"");
      sscanf(buffer,"%i %f %i %f %s",&setvalmin,&valmin,&setvalmax,&valmax,buffer2);
      if(strcmp(buffer,"")!=0){
        for(i=0;i<niso_bounds;i++){
          if(strcmp(isobounds[i].datalabel,buffer2)!=0)continue;
          isobounds[i].setchopmin=setvalmin;
          isobounds[i].setchopmax=setvalmax;
          isobounds[i].chopmin=valmin;
          isobounds[i].chopmax=valmax;
          break;
        }
      }
      else{
        for(i=0;i<niso_bounds;i++){
          isobounds[i].setchopmin=setvalmin;
          isobounds[i].setchopmax=setvalmax;
          isobounds[i].chopmin=valmin;
          isobounds[i].chopmax=valmax;
        }
      }
      continue;
    }
    if(match(buffer,"V_BOUNDARY",10)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i %f %s",&setpatchmin,&patchmin,&setpatchmax,&patchmax,buffer2);
      if(strcmp(buffer2,"")!=0)local2globalpatchbounds(buffer2);
      continue;
    }
    if(match(buffer,"V_ZONE",6)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i %f",&setzonemin,&zonemin,&setzonemax,&zonemax);
      continue;
    }
    if(match(buffer,"V_TARGET",8)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i %f %s",&settargetmin,&targetmin,&settargetmax,&targetmax,buffer2);
      continue;
    }
    if(match(buffer,"OUTLINEMODE",11)==1){
	    fgets(buffer,255,stream);
	    sscanf(buffer,"%i",&highlight_flag);
      if(nmeshes<2&&highlight_flag!=0)highlight_flag=1;
      continue;
    }
    if(match(buffer,"PARTFRAMESTEP",13)==1){
	    fgets(buffer,255,stream);
	    sscanf(buffer,"%i",&partframestep);
	    if(partframestep<1)partframestep=1;
      partframeskip=partframestep-1;
      continue;
    }
    if(match(buffer,"EVACFRAMESTEP",13)==1){
	    fgets(buffer,255,stream);
	    sscanf(buffer,"%i",&evacframestep);
	    if(evacframestep<1)evacframestep=1;
      evacframeskip=evacframestep-1;
      continue;
    }
    if(match(buffer,"USENISTLOGO",11)==1){
	    fgets(buffer,255,stream);
	    sscanf(buffer,"%i",&use_nistlogo);
      continue;
    }
    if(match(buffer,"SLICEDATAOUT",12)==1){
      {
        int sliceoutflag=0;
  	    fgets(buffer,255,stream);
	      sscanf(buffer,"%i",&sliceoutflag);
        if(sliceoutflag!=0){
          output_slicedata=1;
        }
        else{
          output_slicedata=0;
        }
      }
      continue;
    }
    if(match(buffer,"SLICEFRAMESTEP",14)==1){
	    fgets(buffer,255,stream);
	    sscanf(buffer,"%i",&sliceframestep);
	    if(sliceframestep<1)sliceframestep=1;
      sliceframeskip=sliceframestep-1;
      continue;
    }
    if(match(buffer,"SMOKE3DFRAMESTEP",16)==1){
	    fgets(buffer,255,stream);
	    sscanf(buffer,"%i",&smoke3dframestep);
	    if(smoke3dframestep<1)smoke3dframestep=1;
      smoke3dframeskip=smoke3dframestep-1;
      continue;
    }
    if(match(buffer,"SMOKE3DZIPSTEP",14)==1){
	    fgets(buffer,255,stream);
	    sscanf(buffer,"%i",&smoke3dzipstep);
	    if(smoke3dzipstep<1)smoke3dzipstep=1;
      smoke3dzipskip=smoke3dzipstep-1;
      continue;
    }
    if(match(buffer,"SLICEZIPSTEP",12)==1){
	    fgets(buffer,255,stream);
	    sscanf(buffer,"%i",&slicezipstep);
	    if(slicezipstep<1)slicezipstep=1;
      slicezipskip=slicezipstep-1;
      continue;
    }
    if(match(buffer,"ISOZIPSTEP",10)==1){
	    fgets(buffer,255,stream);
	    sscanf(buffer,"%i",&isozipstep);
	    if(isozipstep<1)isozipstep=1;
      isozipskip=isozipstep-1;
      continue;
    }
    if(match(buffer,"BOUNDFRAMESTEP",14)==1){
	    fgets(buffer,255,stream);
	    sscanf(buffer,"%i",&boundframestep);
	    if(boundframestep<1)boundframestep=1;
      boundframeskip=boundframestep-1;
      continue;
    }
    if(match(buffer,"BOUNDZIPSTEP",12)==1){
	    fgets(buffer,255,stream);
	    sscanf(buffer,"%i",&boundzipstep);
	    if(boundzipstep<1)boundzipstep=1;
      boundzipskip=boundzipstep-1;
      continue;
    }
    if(match(buffer,"PARTPOINTSTEP",13)==1){
		  fgets(buffer,255,stream);
		  sscanf(buffer,"%i",&partpointstep);
		  if(partpointstep<1)partpointstep=1;
      partpointskip=partpointstep-1;
      continue;
    }
    if(match(buffer,"MAXPOINTS",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&mxpoints_ini);
      continue;
    }
    if(match(buffer,"MAXFRAMES",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&mxframes);
      continue;
    }
    if(match(buffer,"MSCALE",6)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",mscale,mscale+1,mscale+2);
      continue;
    }
    if(match(buffer,"CLIP",4)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f",&nearclip,&farclip);
      continue;
    }

    if(match(buffer,"SHOWTRACERSALWAYS",17)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&show_tracers_always);
      if(show_tracers_always!=1)show_tracers_always=0;
      continue;
    }
    if(localfile==1&&match(buffer,"PROPINDEX",9)==1){
      int nvals;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&nvals);
      for(i=0;i<nvals;i++){
        propdata *propi;
        int ind, val;

        fgets(buffer,255,stream);
        sscanf(buffer,"%i %i",&ind,&val);
        if(ind<0||ind>npropinfo-1)continue;
        propi = propinfo + ind;
        if(val<0||val>propi->nsmokeview_ids-1)continue;
        propi->smokeview_id=propi->smokeview_ids[val];
        propi->smv_object=propi->smv_objects[val];
      }
      for(i=0;i<npartclassinfo;i++){
        part5class *partclassi;

        partclassi = partclassinfo + i;
        update_partclass_depend(partclassi);

      }
    }
    if(localfile==1&&match(buffer,"PART5CLASSVIS",13)==1){
      int ntemp;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&ntemp);

      for(j=0;j<ntemp;j++){
        part5class *partclassj;

        if(j>npartclassinfo)break;

        partclassj = partclassinfo + j;
        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&partclassj->vis_type);
      }
    }
    if(match(buffer,"PART5COLOR",10)==1){
      for(i=0;i<npart5prop;i++){
        part5prop *propi;

        propi = part5propinfo + i;
        propi->display=0;
      }
      part5colorindex=0;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&i);
      if(i>=0&&i<npart5prop){
        part5prop *propi;

        part5colorindex=i;
        propi = part5propinfo + i;
        propi->display=1;
      }
      continue;
    }

    if(match(buffer,"PART5PROPDISP",13)==1){
      char *token;

      for(i=0;i<npart5prop;i++){
        part5prop *propi;

        propi = part5propinfo + i;
        fgets(buffer,255,stream);
    
        trim(buffer);
        token=strtok(buffer," ");
        j=0;
        while(token!=NULL&&j<npartclassinfo){
          int visval;

          sscanf(token,"%i",&visval);
          propi->class_vis[j]=visval;
          token=strtok(NULL," ");
          j++;
        }
      }
      CheckMemory;
      continue;
    }

    if(match(buffer,"COLORBAR",8) == 1 && match(buffer,"COLORBARFLIP",12)!=1){
      float *rgb_ini_copy;
      
      CheckMemory;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %i %i",&nrgb_ini,&usetexturebar,&colorbar_select_index);
      FREEMEMORY(rgb_ini);
      if(NewMemory((void **)&rgb_ini,4*nrgb_ini*sizeof(float))==0)return 2;
      rgb_ini_copy=rgb_ini;
      for(nn=0;nn<nrgb_ini;nn++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%f %f %f ",rgb_ini_copy,rgb_ini_copy+1,rgb_ini_copy+2);
        rgb_ini_copy+=3;
      }
      initrgb();
      if(colorbar_select_index>=0&&colorbar_select_index<=255){
        update_colorbar_select_index=1;
      }
      continue;
    }
    if(match(buffer,"COLOR2BAR",9) == 1){
      float *rgb_ini_copy;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&nrgb2_ini);
      if(nrgb2_ini<8){
        printf("*** fatal error: must have at lease 8 colors in COLOR2BAR\n");
        exit(1);
      }
      FREEMEMORY(rgb2_ini);
      if(NewMemory((void **)&rgb2_ini,4*nrgb_ini*sizeof(float))==0)return 2;
      rgb_ini_copy=rgb2_ini;
      for(nn=0;nn<nrgb2_ini;nn++){
        fgets(buffer,255,stream);
        sscanf(buffer,"%f %f %f ",rgb_ini_copy,rgb_ini_copy+1,rgb_ini_copy+2);
        rgb_ini_copy+=3;
      }
      continue;
    }
    if(match(buffer,"P3DSURFACETYPE",14)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&p3dsurfacetype);
      continue;
    }
    if(match(buffer,"P3DSURFACESMOOTH",16)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&p3dsurfacesmooth);
      continue;
    }
    if(match(buffer,"CULLFACES",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&cullfaces);
      continue;
    }
    if(match(buffer,"PARTPOINTSIZE",13)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&partpointsize);
      continue;
    }
    if(match(buffer,"ISOPOINTSIZE",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&isopointsize);
      continue;
    }
    if(match(buffer,"ISOLINEWIDTH",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&isolinewidth);
      continue;
    }
    if(match(buffer,"PLOT3DPOINTSIZE",15)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&plot3dpointsize);
      continue;
    }
    if(match(buffer,"PLOT3DLINEWIDTH",15)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&plot3dlinewidth);
      continue;
    }
    if(match(buffer,"VECTORPOINTSIZE",15)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&vectorpointsize);
      continue;
    }
    if(match(buffer,"VECTORLINEWIDTH",15)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&vectorlinewidth);
      continue;
    }

    if(match(buffer,"STREAKLINEWIDTH",15)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&streaklinewidth);
      continue;
    }
    if(match(buffer,"LINEWIDTH",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&linewidth);
      continue;
    }
    if(match(buffer,"VENTLINEWIDTH",13)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&ventlinewidth);
      continue;
    }
    if(match(buffer,"SLICEOFFSET",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&sliceoffset_factor);
      continue;
    }
    if(match(buffer,"TITLESAFE",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&titlesafe_offset);
      continue;
    }
    if(match(buffer,"VENTOFFSET",10)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&ventoffset_factor);
      continue;
    }
    if(match(buffer,"AXISSMOOTH",10)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&axissmooth);
      continue;
    }
    if(match(buffer,"AXISNUM",7)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&axisnum);
      continue;
    }
    if(match(buffer,"SHOWBLOCKS",10)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&visBlocks);
      continue;
    }
    if(match(buffer,"SHOWSENSORS",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %i ",&visSensor,&visSensorNorm);
      if(visSensor!=0)visSensor=1;
      if(visSensorNorm!=0)visSensorNorm=1;
      continue;
    }
    if(match(buffer,"AVATAREVAC",10)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&iavatar_evac);
      continue;
    }
    if(match(buffer,"SHOWVENTS",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %i %i ",&visVents,&visVentLines,&visVentSolid);
      if(visVents==0){
        visVentLines=0;
        visVentSolid=0;
      }
      if(visVentSolid==1)visVentLines=0;
      if(visVentLines==1)visVentSolid=0;
      continue;
    }
    if(match(buffer,"SHOWTIMELABEL",13)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&visTimelabel);
      continue;
    }
    if(match(buffer,"SHOWHMSTIMELABEL",16)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&vishmsTimelabel);
      continue;
    }
    if(match(buffer,"SHOWFRAMELABEL",14)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&visFramelabel);
      continue;
    }
    if(match(buffer,"SHOWHRRLABEL",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&visHRRlabel);
      continue;
    }
    if(match(buffer,"RENDERFILETYPE",14)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&renderfiletype);
#ifndef pp_JPEG
      if(renderfiletype==1){
        printf("*** warning: JPEG not supported, render filetype changed to PNG\n");
        renderfiletype=0;
      }
#endif
      continue;
    }
    if(match(buffer,"SHOWGRIDLOC",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&visgridloc);
      continue;
    }
//    if(match(buffer,"SHOWCADANDGRID",14)==1){
//      fgets(buffer,255,stream);
//      sscanf(buffer,"%i ",&show_cad_and_grid);
//      continue;
//    }
    if(match(buffer,"SHOWFLOOR",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&visFloor);
      continue;
    }
    if(match(buffer,"SPEED",5)==1){
      fgets(buffer,255,stream);
    //  sscanf(buffer,"%f %f",&speed_crawl,&speed_walk);
      continue;
    }
    if(match(buffer,"FONTSIZE",8)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&fontindex);
      if(fontindex<0)fontindex=0;
      if(fontindex>1)fontindex=1;
      FontMenu(fontindex);
      continue;
    }
    if(match(buffer,"ZOOM",4)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f ",&zoomindex,&zoom);
      if(zoomindex!=-1){
        if(zoomindex<0)zoomindex=2;
        if(zoomindex>4)zoomindex=2;
        zoom=zooms[zoomindex];
      }
      else{
        if(zoom<zooms[0]){
          zoom=zooms[0];
          zoomindex=0;
        }
        if(zoom>zooms[4]){
          zoom=zooms[4];
          zoomindex=4;
        }
      }
      updatezoommenu=1;
      ZoomMenu(zoomindex);
      continue;
    }
    if(match(buffer,"APERATURE",9)==1||match(buffer,"APERTURE",8)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&apertureindex);
      if(apertureindex<0)apertureindex=0;
      if(apertureindex>4)apertureindex=4;
      ApertureMenu(apertureindex);
      continue;
    }
    if(match(buffer,"SHOWWALLS",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&visWalls);
      continue;
      }
    if(match(buffer,"SHOWCEILING",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&visCeiling);
      continue;
      }
    if(match(buffer,"SHOWTITLE",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&visTitle0);
      continue;
      }
    if(match(buffer,"SHOWNORMALWHENSMOOTH",20)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&visSmoothAsNormal);
      continue;
      }
    if(match(buffer,"SHOWTRANSPARENT",15)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&visTransparentBlockage);
      continue;
      }
    if(match(buffer,"VECTORPOINTSIZE",15)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&vectorpointsize);
      continue;
      }
    if(match(buffer,"VECTORSKIP",10)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&vectorskip);
      if(vectorskip<1)vectorskip=1;
      continue;
      }
    if(match(buffer,"SPRINKLERABSSIZE",16)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&sprinklerabssize);
      continue;
      }
    if(match(buffer,"SENSORABSSIZE",13)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&sensorabssize);
      continue;
      }
    if(match(buffer,"SENSORRELSIZE",13)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f ",&sensorrelsize);
      if(sensorrelsize<sensorrelsizeMIN)sensorrelsize=sensorrelsizeMIN;
      continue;
      }
    if(match(buffer,"SETBW",5)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&setbw);
      continue;
      }
    if(match(buffer,"FLIP",4)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&background_flip);
      continue;
      }
    if(match(buffer,"COLORBARFLIP",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&colorbarflip);
      continue;
      }
    if(match(buffer,"TRANSPARENT",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i ",&transparentflag);
//      transparentflagSAVE=transparentflag;
      continue;
      }
    if(match(buffer,"VENTCOLOR",9)==1){
      float ventcolor_temp[4];

      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",ventcolor_temp,ventcolor_temp+1,ventcolor_temp+2);
      ventcolor_temp[3]=1.0;
      ventcolor=getcolorptr(ventcolor_temp);
      updatefaces=1;
      updateindexcolors=1;
      continue;
      }
    if(match(buffer,"STATICPARTCOLOR",15)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",static_color,static_color+1,static_color+2);
      continue;
      }
    if(match(buffer,"HEATOFFCOLOR",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",heatoffcolor,heatoffcolor+1,heatoffcolor+2);
      continue;
      }
    if(match(buffer,"HEATONCOLOR",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",heatoncolor,heatoncolor+1,heatoncolor+2);
      continue;
      }
    if(match(buffer,"SENSORCOLOR",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",sensorcolor,sensorcolor+1,sensorcolor+2);
      continue;
      }
    if(match(buffer,"SENSORNORMCOLOR",15)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",sensornormcolor,sensornormcolor+1,sensornormcolor+2);
      continue;
      }
    if(match(buffer,"SPRINKONCOLOR",13)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",sprinkoncolor,sprinkoncolor+1,sprinkoncolor+2);
      continue;
      }
    if(match(buffer,"SPRINKOFFCOLOR",14)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",sprinkoffcolor,sprinkoffcolor+1,sprinkoffcolor+2);
      continue;
      }
    if(match(buffer,"BACKGROUNDCOLOR",15)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",backgroundbasecolor,backgroundbasecolor+1,backgroundbasecolor+2);
      continue;
      }
    if(match(buffer,"FOREGROUNDCOLOR",15)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",foregroundbasecolor,foregroundbasecolor+1,foregroundbasecolor+2);
      continue;
      }
    if(match(buffer,"BLOCKCOLOR",10)==1){
      float blockcolor_temp[4];

      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",blockcolor_temp,blockcolor_temp+1,blockcolor_temp+2);
      blockcolor_temp[3]=1.0;
      block_ambient2=getcolorptr(blockcolor_temp);
      updatefaces=1;
      updateindexcolors=1;
      continue;
      }
    if(match(buffer,"BLOCKLOCATION",13)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&blocklocation);
      continue;
      }
    if(match(buffer,"SHOWOPENVENTS",13)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %i",&visOpenVents,&visOpenVentsAsOutline);
      continue;
      }
    if(match(buffer,"SHOWDUMMYVENTS",14)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visDummyVents);
      continue;
      }
    if(match(buffer,"SHOWTICKS",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visTicks);
      continue;
      }
    if(match(buffer,"USERTICKS",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %i %i %i %i %i",&vis_user_ticks,&auto_user_tick_placement,&user_tick_sub,
        &user_tick_show_x,&user_tick_show_y,&user_tick_show_z);
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",user_tick_origin,user_tick_origin+1,user_tick_origin+2);
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",user_tick_min,user_tick_min+1,user_tick_min+2);
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",user_tick_max,user_tick_max+1,user_tick_max+2);
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",user_tick_step,user_tick_step+1,user_tick_step+2);
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %i %i",&user_tick_show_x,&user_tick_show_y,&user_tick_show_z);
      continue;
      }
    if(match(buffer,"SHOWLABELS",10)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visLabels);
      continue;
      }
    if(match(buffer,"BOUNDCOLOR",10)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",boundcolor,boundcolor+1,boundcolor+2);
      continue;
      }
    if(match(buffer,"TIMEBARCOLOR",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",timebarcolor,timebarcolor+1,timebarcolor+2);
      continue;
      }
    if(match(buffer,"P3CONT2D",8)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&p3cont2d);
#ifdef pp_LINE
      if(p3cont2d>2)p3cont2d=2;
#else
      if(p3cont2d>1)p3cont2d=1;
#endif
      continue;
      }
    if(match(buffer,"P3VIEW",6)==1){
      for(i=0;i<nmeshes;i++){
        mesh *meshi;

        meshi = meshinfo + i;
        fgets(buffer,255,stream);
        sscanf(buffer,"%i %i %i %i %i %i",
          &meshi->visx,&meshi->plotx,&meshi->visy,&meshi->ploty,&meshi->visz,&meshi->plotz);
        if(meshi->visx!=0)meshi->visx=1;
        if(meshi->visy!=0)meshi->visy=1;
        if(meshi->visy!=0)meshi->visz=1;
        if(meshi->plotx<0)meshi->plotx=0;
        if(meshi->plotx>meshi->ibar)meshi->plotx=meshi->ibar;
        if(meshi->ploty>meshi->jbar)meshi->ploty=meshi->jbar;
        if(meshi->plotz>meshi->kbar)meshi->plotz=meshi->kbar;
      }
      continue;
      }
    if(match(buffer,"P3CONT3DSMOOTH",14)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&p3cont3dsmooth);
      continue;
      }
    if(match(buffer,"VECTORLENGTH",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f",&VECFRACTION);
      continue;
      }
    if(match(buffer,"SURFINC",7)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&surfincrement);
      continue;
      }
    if(match(buffer,"FRAMERATE",9)==1&&match(buffer,"FRAMERATEVALUE",14)!=1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visFramerate);
      continue;
      }
    if(match(buffer,"SHOWFRAMERATE",13)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visFramerate);
      continue;
      }
    if(match(buffer,"SHOWFRAME",9)==1&&match(buffer,"SHOWFRAMERATE",13)!=1&&match(buffer,"SHOWFRAMELABEL",14)!=1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visFrame);
      if(isZoneFireModel==1)visFrame=0;
      continue;
    }
    if(match(buffer,"FRAMERATEVALUE",14)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&frameratevalue);
      FrameRateMenu(frameratevalue);
      continue;
      }
    if(match(buffer,"SHOWSPRINKPART",14)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visSprinkPart);
      continue;
      }
    if(match(buffer,"SHOWAXISLABELS",14)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visaxislabels);
      continue;
      }
#ifdef pp_memstatus
    if(match(buffer,"SHOWMEMLOAD",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visAvailmemory);
      continue;
      }
#endif
    if(match(buffer,"SHOWBLOCKLABEL",14)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visBlocklabel);
      continue;
      }
    if(match(buffer,"SHOWVZONE",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visVZone);
      continue;
      }
    if(match(buffer,"SHOWHZONE",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visHZone);
      continue;
      }
    if(match(buffer,"SHOWHAZARDCOLORS",16)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&sethazardcolor);
      continue;
      }
    if(match(buffer,"SHOWSMOKEPART",13)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visSmokePart);
      continue;
      }
    if(match(buffer,"RENDEROPTION",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&render_option);
      RenderMenu(render_option);
      continue;
      }
    if(match(buffer,"SHOWISONORMALS",14)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&showisonormals);
      if(showisonormals!=1)showisonormals=0;
      continue;
    }
    if(match(buffer,"SHOWISO",7)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visAIso);
      if(visAIso<0||visAIso>3){
        visAIso=1;
      }
      continue;
      }
    if(trainer_mode==0&&windowresized==0){
      if(match(buffer,"WINDOWWIDTH",11)==1){
        int scrWidth;

        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&scrWidth);
        if(scrWidth<=0){
          scrWidth = glutGet(GLUT_SCREEN_WIDTH);
        }
        if(scrWidth!=screenWidth){
          screenWidth=scrWidth;
          update_screensize=1;
        }
        continue;
        }
      if(match(buffer,"WINDOWHEIGHT",12)==1){
        int scrHeight;

        fgets(buffer,255,stream);
        sscanf(buffer,"%i",&scrHeight);
        if(scrHeight<=0){
          scrHeight = glutGet(GLUT_SCREEN_HEIGHT);
        }
        if(scrHeight!=screenHeight){
          screenHeight=scrHeight;
          update_screensize=1;
        }
        continue;
        }
    }
    if(match(buffer,"SHOWTIMEBAR",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visTimeLabels);
      continue;
      }
    if(match(buffer,"SHOWCOLORBARS",13)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visColorLabels);
      continue;
      }
    if(match(buffer,"EYEVIEW",7)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&eyeview);
      continue;
      }
    if(localfile==1&&match(buffer,"SCRIPTFILE",10)==1){
      if(fgets(buffer2,255,stream)==NULL)break;
      cleanbuffer(buffer,buffer2);
      insert_scriptfile(buffer);
      updatemenu=1;
      continue;
      }
    if(localfile==1&&match(buffer,"SHOWDEVICES",11)==1){
      sv_object *obj_typei;
      char *dev_label;
      int ndevices_ini;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&ndevices_ini);

      for(i=0;i<nobject_defs;i++){
        obj_typei = object_defs[i];
        obj_typei->visible=0;
      }
      for(i=0;i<ndevices_ini;i++){
        fgets(buffer,255,stream);
        trim(buffer);
        dev_label=trim_front(buffer);
        obj_typei=get_object(dev_label);
        if(obj_typei!=NULL){
          obj_typei->visible=1;
        }
      }
    }
    if(localfile==1&&match(buffer,"XYZCLIP",7)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&xyz_clipplane);
      if(xyz_clipplane<0)xyz_clipplane=0;
      if(xyz_clipplane>2)xyz_clipplane=2;
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i %f",&clip_x, &clip_x_val, &clip_X, &clip_X_val);
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i %f",&clip_y, &clip_y_val, &clip_Y, &clip_Y_val);
      fgets(buffer,255,stream);
      sscanf(buffer,"%i %f %i %f",&clip_z, &clip_z_val, &clip_Z, &clip_Z_val);
      updateclipvals=1;
      continue;
      }
    if(match(buffer,"NOPART",6)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&nopart);
      continue;
      }
    if(match(buffer,"WINDOWOFFSET",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&titlesafe_offsetBASE);
      continue;
      }
    if(match(buffer,"SHOWLIGHT0",10)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visLIGHT0);
      UpdateLIGHTS=1;
      continue;
    }
    if(match(buffer,"SHOWLIGHT1",10)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&visLIGHT1);
      UpdateLIGHTS=1;
      continue;
    }
    if(match(buffer,"SHOWLIGHTMENU",13)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&showlightmenu);
      updatemenu=1;
      continue;
    }
    if(match(buffer,"AMBIENTLIGHT",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",ambientlight,ambientlight+1,ambientlight+2);
      UpdateLIGHTS=1;
      continue;
    }
    if(match(buffer,"DIFFUSELIGHT",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",diffuselight,diffuselight+1,diffuselight+2);
      UpdateLIGHTS=1;
      continue;
    }
    if(match(buffer,"ISOFRAMESTEP",12)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&isoframestep);
      if(isoframestep<1)isoframestep=1;
      isoframeskip=isoframestep-1;
      continue;
    }
    if(match(buffer,"LABELSTARTUPVIEW",16)==1){
      char *front;
      camera *ca;

      fgets(buffer,255,stream);
      front=trim_front(buffer);
      trim(front);
      strcpy(label_startup_view,front);
      updategluiview=1;
      continue;
    }
    if(match(buffer,"INPUT_FILE",10) == 1){
      {
        size_t len;

        if(fgets(buffer,255,stream)==NULL)break;
        len=strlen(buffer);
        buffer[len-1]='\0';
        trim(buffer);
        len=strlen(buffer);
 
        FREEMEMORY(INI_fds_filein);
        if(NewMemory((void **)&INI_fds_filein,(unsigned int)(len+1))==0)return 2;
        STRCPY(INI_fds_filein,buffer);
        continue;
      }

    }


      {
        float *eye,*mat,*angle_zx;
        int is_viewpoint4=0;
        int is_viewpoint5=0;

    if(match(buffer,"VIEWPOINT3",10)==1
      ||match(buffer,"VIEWPOINT4",10)==1
      ||match(buffer,"VIEWPOINT5",10)==1
      ){
        int p_type;

        if(match(buffer,"VIEWPOINT4",10)==1)is_viewpoint4=1;
        if(match(buffer,"VIEWPOINT5",10)==1){
          is_viewpoint4=1;
          is_viewpoint5=1;
        }
        eye=camera_ini->eye;
        mat=camera_ini->modelview;
        angle_zx=camera_ini->angle_zx;

        {
          char name_ini[32];
          strcpy(name_ini,"ini");
          init_camera(camera_ini,name_ini);
        }

		    fgets(buffer,255,stream);
		    sscanf(buffer,"%i %i %i",&camera_ini->eyeview,&camera_ini->rotation_index,&camera_ini->view_id);

        {
          float zoom_in;
          int zoomindex_in;

          zoom_in=zoom;
          zoomindex_in=zoomindex;
          fgets(buffer,255,stream);
	  	    sscanf(buffer,"%f %f %f %f %i",eye,eye+1,eye+2,&zoom_in,&zoomindex_in);
          zoom=zoom_in;
          zoomindex=zoomindex_in;
          if(zoomindex!=-1){
           if(zoomindex<0)zoomindex=2;
           if(zoomindex>4)zoomindex=2;
           zoom=zooms[zoomindex];
          }
          else{
            if(zoom<zooms[0]){
              zoom=zooms[0];
              zoomindex=0;
            }
            if(zoom>zooms[4]){
              zoom=zooms[4];
              zoomindex=4;
            }
          }
          updatezoommenu=1;
        }

        p_type=0;
		    fgets(buffer,255,stream);
		    sscanf(buffer,"%f %f %f %i",
          &camera_ini->view_angle, 
          &camera_ini->direction_angle,
          &camera_ini->elevation_angle,
          &p_type);
        if(p_type!=1)p_type=0;
        camera_ini->projection_type=p_type;

		    fgets(buffer,255,stream);
		    sscanf(buffer,"%f %f %f",&camera_ini->xcen,&camera_ini->ycen,&camera_ini->zcen);

		    fgets(buffer,255,stream);
        sscanf(buffer,"%f %f",angle_zx,angle_zx+1);

		    fgets(buffer,255,stream);
        sscanf(buffer,"%f %f %f %f",mat+0,mat+1,mat+2,mat+3);

		    fgets(buffer,255,stream);
        sscanf(buffer,"%f %f %f %f",mat+4,mat+5,mat+6,mat+7);

		    fgets(buffer,255,stream);
        sscanf(buffer,"%f %f %f %f",mat+8,mat+9,mat+10,mat+11);

		    fgets(buffer,255,stream);
        sscanf(buffer,"%f %f %f %f",mat+12,mat+13,mat+14,mat+15);
        if(is_viewpoint5==1){
          camera *ci;

          ci = camera_ini;
  		    fgets(buffer,255,stream);
          sscanf(buffer,"%i %i %i %i %i %i %i",
            &ci->xyz_clipplane,
            &ci->clip_x,&ci->clip_y,&ci->clip_z,
            &ci->clip_X,&ci->clip_Y,&ci->clip_Z);
  		    fgets(buffer,255,stream);
          sscanf(buffer,"%f %f %f %f %f %f",
            &ci->clip_x_val,&ci->clip_y_val,&ci->clip_z_val,
            &ci->clip_X_val,&ci->clip_Y_val,&ci->clip_Z_val);

        }
        if(is_viewpoint4==1){
  		    fgets(buffer,255,stream);
          trim(buffer);
          strcpy(camera_ini->name,buffer);
          init_camera_list();
          {
            //*** following code shouldn't be here but leaving it here commented
            //    in case it is really necessary

            //  camera *cam;

            insert_camera(&camera_list_first,camera_ini,buffer);
           // if(cam!=NULL){
           //   cam->view_id=camera_ini->view_id;
           // }
          }
        }

        enable_reset_saved_view();
        camera_ini->dirty=1;
        camera_ini->defined=1;
        continue;
      }
    }

    if(match(buffer,"ISOCOLORS",9)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f",&iso_shininess,&iso_transparency);
      fgets(buffer,255,stream);
      sscanf(buffer,"%f %f %f",iso_specular,iso_specular+1,iso_specular+2);
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&n_iso_ambient_ini);
      if(n_iso_ambient_ini>0){
        FREEMEMORY(iso_ambient_ini);
        if(NewMemory((void**)&iso_ambient_ini,n_iso_ambient_ini*4*sizeof(float))==0)return 2;
        for(nn=0;nn<n_iso_ambient_ini;nn++){
          fgets(buffer,255,stream);
          sscanf(buffer,"%f %f %f",iso_ambient_ini+4*nn,iso_ambient_ini+4*nn+1,iso_ambient_ini+4*nn+2);
          iso_ambient_ini[4*nn+3]=iso_transparency;
        }
        iso_ambient_ini[3]=1.0;
      }
      continue;
    }
    if(match(buffer,"UNITCLASSES",11)==1){
      int nuc;

      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&nuc);
      for(i=0;i<nuc;i++){
        int active;

        fgets(buffer,255,stream);
        if(i>nunitclasses-1)continue;
        sscanf(buffer,"%i",&active);
        unitclasses[i].active=active;
      }
      continue;
    }
    if(match(buffer,"SMOOTHLINES",11)==1){
      fgets(buffer,255,stream);
      sscanf(buffer,"%i",&antialiasflag);
      continue;
      }
    if(localfile==1&&match(buffer,"VIEWTIMES",9) == 1){
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%f %f %i",&view_tstart,&view_tstop,&view_ntimes);
      if(view_ntimes<2)view_ntimes=2;
      ReallocTourMemory();
      continue;
    }
#ifdef pp_SHOOTER
    if(localfile==1&&match(buffer,"SHOOTER",7) == 1){
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%f %f %f",shooter_xyz,shooter_xyz+1,shooter_xyz+2);
      
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%f %f %f",shooter_dxyz,shooter_dxyz+1,shooter_dxyz+2);
      
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%f %f %f",shooter_uvw,shooter_uvw+1,shooter_uvw+2);

      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%f %f %f",&shooter_velmag,&shooter_veldir,&shooterpointsize);
      
      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%i %i %i %i %i",&shooter_fps,&shooter_vel_type,&shooter_nparts,&visShooter,&shooter_cont_update);

      if(fgets(buffer,255,stream)==NULL)break;
      sscanf(buffer,"%f %f",&shooter_duration,&shooter_v_inf);
      continue;
    }
#endif
    {
      int nkeyframes;
      float key_time, key_xyz[3], key_az_path, key_view[3], params[3], zzoom, key_elev_path;
      float t_globaltension, key_bank;
      int t_globaltension_flag;
      int viewtype,uselocalspeed;
      float *col;

      if(match(buffer,"ADJUSTALPHA",11)==1){
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%i",&adjustalphaflag);
        continue;
      }
      if(match(buffer,"SMOKECULL",9)==1){
        if(fgets(buffer,255,stream)==NULL)break;
#ifdef pp_CULL
        sscanf(buffer,"%i",&cullsmoke);
        if(cullsmoke!=0)cullsmoke=1;
#else
        sscanf(buffer,"%i",&smokecullflag);
#endif
        continue;
      }
      if(match(buffer,"SMOKESKIP",9)==1){
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%i",&smokeskipm1);
        continue;
      }
      if(match(buffer,"SMOKESHADE",10)==1){
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%i",&smoke_shade);
        continue;
      }
      if(match(buffer,"SMOKETHICK",10)==1){
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%i",&smoke3d_thick);
        continue;
      }
#ifdef pp_GPU
      if(match(buffer,"SMOKERTHICK",11)==1){
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%f",&smoke3d_rthick);
        if(smoke3d_rthick<1.0)smoke3d_rthick=1.0;
        if(smoke3d_rthick>255.0)smoke3d_rthick=255.0;
        smoke3d_thick=log_base2(smoke3d_rthick);
        continue;
      }
#endif
      if(match(buffer,"FIRECOLOR",9)==1){
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%i %i %i",&fire_red,&fire_green,&fire_blue);
        continue;
      }
      if(match(buffer,"FIREDEPTH",9)==1){
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%f",&fire_halfdepth);
        continue;
      }


      if(match(buffer,"VIEWTOURFROMPATH",16)==1){
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%i",&viewtourfrompath);
        continue;
      }
      if(match(buffer,"VIEWALLTOURS",12)==1){
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%i",&viewalltours);
        continue;
      }
      if(match(buffer,"SHOWTOURROUTE",13)==1){
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%i",&edittour);
        continue;
      }
      if(match(buffer,"TIMEOFFSET",10)==1){
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%f",&timeoffset);
        continue;
      }
      if(match(buffer,"SHOWPATHNODES",13)==1){
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%i",&show_path_knots);
        continue;
      }
      if(match(buffer,"SHOWIGNITION",12)==1){
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%i %i",&vis_threshold,&vis_onlythreshold);
        continue;
      }
      if(match(buffer,"SHOWTHRESHOLD",13)==1){
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%i %i %f",&vis_threshold,&vis_onlythreshold,&temp_threshold);
        continue;
      }
      if(match(buffer,"TOURCONSTANTVEL",15)==1){
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%i",&tour_constant_vel);
        continue;
      }
      if(match(buffer,"TOUR_AVATAR",11)==1){
        if(fgets(buffer,255,stream)==NULL)break;
//        sscanf(buffer,"%i %f %f %f %f",&tourlocus_type,tourcol_avatar,tourcol_avatar+1,tourcol_avatar+2,&tourrad_avatar);
//        if(tourlocus_type!=0)tourlocus_type=1;
        continue;
      }

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ GENCOLORBAR ++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

      if(match(buffer,"GCOLORBAR",9) == 1){
        colorbardata *cbi;
        int r1, g1, b1;
        int n;
        int ncolorbarini;

        fgets(buffer,255,stream);
        ncolorbarini=0;
        sscanf(buffer,"%i",&ncolorbarini);

        initdefaultcolorbars();

        ncolorbars=ndefaultcolorbars+ncolorbarini;
        if(ncolorbarini>0)ResizeMemory((void **)&colorbarinfo,ncolorbars*sizeof(colorbardata));

        for(n=ndefaultcolorbars;n<ncolorbars;n++){
          int extreme, rgbmin[3],rgbmax[3];

          cbi = colorbarinfo + n;

          fgets(buffer,255,stream);
          trim(buffer);
          strcpy(cbi->label,buffer);

          fgets(buffer,255,stream);
          sscanf(buffer,"%i %i",&cbi->nnodes,&cbi->nodehilight);
          if(cbi->nnodes<0)cbi->nnodes=0;
          if(cbi->nodehilight<0||cbi->nodehilight>=cbi->nnodes){
            cbi->nodehilight=0;
          }

          cbi->label_ptr=cbi->label;
          for(i=0;i<cbi->nnodes;i++){
            int icbar;

            fgets(buffer,255,stream);
            r1=-1; g1=-1; b1=-1; 
            sscanf(buffer,"%i %i %i %i",&icbar,&r1,&g1,&b1);
            cbi->index_node[i]=icbar;
            nn = 3*i;
            cbi->rgb_node[nn  ]=r1;
            cbi->rgb_node[nn+1]=g1;
            cbi->rgb_node[nn+2]=b1;
          }

          remapcolorbar(cbi);
        }
    }

      if(match(buffer,"TOURCOLORS",10)==1){
        col=tourcol_selectedpathline;
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%f %f %f",col,col+1,col+2);

        col=tourcol_selectedpathlineknots;
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%f %f %f",col,col+1,col+2);

        col=tourcol_selectedknot;
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%f %f %f",col,col+1,col+2);

        col=tourcol_pathline;
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%f %f %f",col,col+1,col+2);

        col=tourcol_pathknots;
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%f %f %f",col,col+1,col+2);

        col=tourcol_text;
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%f %f %f",col,col+1,col+2);

        col=tourcol_avatar;
        if(fgets(buffer,255,stream)==NULL)break;
        sscanf(buffer,"%f %f %f",col,col+1,col+2);

        continue;
      }

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ LABEL ++++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

    if(localfile==1&&match(buffer,"LABEL",5) == 1){
      nlabels++;
      ResizeMemory((void **)&labelinfo,(nlabels)*sizeof(labeldata));

      /*
      LABEL
      x y z r g b tstart tstop  
      label

      */
      {
        float *xyz, *rgbtemp, *tstart_stop;
        labeldata *labeli;

        labeli = labelinfo + nlabels-1;

        xyz = labeli->xyz;
        rgbtemp = labeli->rgb;
        tstart_stop = labeli->tstart_stop;

        fgets(buffer,255,stream);
        rgbtemp[0]=-1.0;
        rgbtemp[1]=-1.0;
        rgbtemp[2]=-1.0;
        rgbtemp[3]=1.0;
        tstart_stop[0]=-1.0;
        tstart_stop[1]=-1.0;
        sscanf(buffer,"%f %f %f %f %f %f %f %f",
          xyz,xyz+1,xyz+2,
          rgbtemp,rgbtemp+1,rgbtemp+2,
          tstart_stop,tstart_stop+1);
        if(rgbtemp[0]<0.0||rgbtemp[1]<0.0||rgbtemp[2]<0.0||rgbtemp[0]>1.0||rgbtemp[1]>1.0||rgbtemp[2]>1.0){
          labeli->useforegroundcolor=1;
        }
        else{
          labeli->useforegroundcolor=0;
        }
        fgets(buffer,255,stream);
        strcpy(labeli->label,buffer);
      }
      continue;
    }

  /*
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ++++++++++++++++++++++ TICKS ++++++++++++++++++++++++++++++++
    +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  */

/*
typedef struct {
  float begin[3],end[3],length;
  float dxyz[3],dlength;
  int dir,nbars;
} tickdata;
*/

    if(localfile==1&&match(buffer,"TICKS",5) == 1){
      nticks++;
      ResizeMemory((void **)&tickinfo,(nticks)*sizeof(tickdata));

      {
        tickdata *ticki;
        float *begt, *endt;
        int *nbarst;
        float term;
        float length=0.0;
        float *dxyz;
        float sum;

        ticki = tickinfo + nticks - 1;
        begt = ticki->begin;
        endt = ticki->end;
        nbarst=&ticki->nbars;
        dxyz = ticki->dxyz;
        

        /*
        TICKS
        b1 b2 b3 e1 e2 e3 nb
        ticklength tickdir tickcolor (r g b) tickwidth 
        */
        if(fgets(buffer,255,stream)==NULL)break;
        *nbarst=0;
        sscanf(buffer,"%f %f %f %f %f %f %i",begt,begt+1,begt+2,endt,endt+1,endt+2,nbarst);
        if(*nbarst<1)*nbarst=1;
        if(fgets(buffer,255,stream)==NULL)break;
        {
          float *rgbtemp;

          rgbtemp=ticki->rgb;
          rgbtemp[0]=-1.0;
          rgbtemp[1]=-1.0;
          rgbtemp[2]=-1.0;
          ticki->width=-1.0;
          sscanf(buffer,"%f %i %f %f %f %f",&ticki->dlength,&ticki->dir,rgbtemp,rgbtemp+1,rgbtemp+2,&ticki->width);
          if(rgbtemp[0]<0.0||rgbtemp[0]>1.0||
             rgbtemp[1]<0.0||rgbtemp[1]>1.0||
             rgbtemp[2]<0.0||rgbtemp[2]>1.0){
            ticki->useforegroundcolor=1;
          }
          else{
            ticki->useforegroundcolor=0;
          }
          if(ticki->width<0.0)ticki->width=1.0;
        }
        for(i=0;i<3;i++){
          term = endt[i]-begt[i];
          length += term*term;
        }
        if(length<=0.0){
          endt[0]=begt[0]+1.0;
          length = 1.0;
        }
        ticki->length=sqrt(length);
        dxyz[0] =  0.0;
        dxyz[1] =  0.0;
        dxyz[2] =  0.0;
        switch (ticki->dir){
        case 1:
        case -1:
          dxyz[0]=1.0;
          break;
        case 2:
        case -2:
          dxyz[1]=1.0;
          break;
        case 3:
        case -3:
          dxyz[2]=1.0;
          break;
        default:
          ASSERT(FFALSE);
          break;
        }
        if(ticki->dir<0){
          for(i=0;i<3;i++){
            dxyz[i]=-dxyz[i];
          }
        }
        sum = 0.0;
        sum = dxyz[0]*dxyz[0] + dxyz[1]*dxyz[1] + dxyz[2]*dxyz[2];
        if(sum>0.0){
          sum=sqrt(sum);
          dxyz[0] *= (ticki->dlength/sum);
          dxyz[1] *= (ticki->dlength/sum);
          dxyz[2] *= (ticki->dlength/sum);
        }
      }
      continue;
    }

      if(localfile==1){
        tours_flag=0;
        if(match(buffer,"TOURS",5)==1)tours_flag=1;
        if( tours_flag==1){
          if(ntours>0){
            for(i=0;i<ntours;i++){
              tourdata *touri;

              touri = tourinfo + i;
              freetour(touri);
            }
            FREEMEMORY(tourinfo);
          }
          ntours=0;

          fgets(buffer,255,stream);
          sscanf(buffer,"%i",&ntours);
          ntours++;
          if(ntours>0){
            if(NewMemory( (void **)&tourinfo, ntours*sizeof(tourdata))==0)return 2;
            for(i=0;i<ntours;i++){
              tourdata *touri;

              touri=tourinfo+i;
              touri->path_times=NULL;
              touri->pathnodes=NULL;
            }
          }
          ReallocTourMemory();
          init_circulartour();
          {
            keyframe *thisframe, *addedframe;
            tourdata *touri;
            int glui_avatar_index;

            for(i=1;i<ntours;i++){
              touri = tourinfo + i;
              inittour(touri);
              fgets(buffer,255,stream);
              trim(buffer);
              strcpy(touri->label,buffer);

              t_globaltension=touri->global_tension;
              t_globaltension_flag=touri->global_tension_flag;
              fgets(buffer,255,stream);
              glui_avatar_index=0;
              sscanf(buffer,"%i %i %f %i %i",
                &nkeyframes,&t_globaltension_flag,&t_globaltension,&glui_avatar_index,&touri->display2);
              if(glui_avatar_index<0)glui_avatar_index=0;
              if(glui_avatar_index>navatar_types-1)glui_avatar_index=navatar_types-1;
              touri->glui_avatar_index=glui_avatar_index;
              if(touri->display2!=1)touri->display2=0;
              touri->global_tension_flag=t_globaltension_flag;
              touri->global_tension=t_globaltension;
              touri->nkeyframes=nkeyframes;

              if(NewMemory( (void **)&touri->keyframe_times, nkeyframes*sizeof(float))==0)return 2;
              if(NewMemory((void **)&touri->pathnodes,view_ntimes*sizeof(pathdata))==0)return 2;
              if(NewMemory((void **)&touri->path_times,view_ntimes*sizeof(float))==0)return 2;

              thisframe=&touri->first_frame;
              for(j=0;j<nkeyframes;j++){
                params[0] = 0.0;
                params[1] = 0.0;
                params[2] = 0.0;
                key_view[0]=0.0;
                key_view[1]=0.0;
                key_view[2]=0.0;
                key_bank=0.0;
                key_az_path = 0.0;
                key_elev_path=0.0;
                viewtype=0;
                zzoom=1.0;
                uselocalspeed=0;
                fgets(buffer,255,stream);

                sscanf(buffer,"%f %f %f %f %i",
                  &key_time,
                  key_xyz,key_xyz+1,key_xyz+2,
                  &viewtype);

                if(viewtype==0){
                  sscanf(buffer,"%f %f %f %f %i %f %f %f %f %f %f %f %i",
                  &key_time,
                  key_xyz,key_xyz+1,key_xyz+2,
                  &viewtype, &key_az_path, &key_elev_path,&key_bank, 
                  params,params+1,params+2,
                  &zzoom,&uselocalspeed);
                }
                else{
                  sscanf(buffer,"%f %f %f %f %i %f %f %f %f %f %f %f %i",
                  &key_time,
                  key_xyz,key_xyz+1,key_xyz+2,
                  &viewtype, key_view, key_view+1, key_view+2,
                  params,params+1,params+2,
                  &zzoom,&uselocalspeed);
                }
                if(zzoom<0.25)zzoom=0.25;
                if(zzoom>4.00)zzoom=4.0;
                addedframe=add_frame(thisframe, key_time, key_xyz, key_az_path, key_elev_path, 
                  key_bank, params, viewtype,zzoom,key_view);
                thisframe=addedframe;
                touri->keyframe_times[j]=key_time;
              }
            }
          }
          if(tours_flag==1){
            for(i=0;i<ntours;i++){
              tourdata *touri;

              touri=tourinfo+i;
              touri->first_frame.next->prev=&touri->first_frame;
              touri->last_frame.prev->next=&touri->last_frame;
            }
            updatetourmenulabels();
            createtourpaths();
            updatetimes();
            plotstate=getplotstate(DYNAMIC_PLOTS);
            selectedtour_index=-1;
            selected_frame=NULL;
            selected_tour=NULL;
            if(viewalltours==1)TourMenu(-3);
          }
          else{
            ntours=0;
          }
          strcpy(buffer,"1.00000 1.00000 2.0000 0");
          trimmzeros(buffer);
          continue;
        }
        if(match(buffer,"TOURINDEX",9)){
          if(fgets(buffer,255,stream)==NULL)break;
          sscanf(buffer,"%i",&selectedtour_index_ini);
          if(selectedtour_index_ini<0)selectedtour_index_ini=-1;
          update_selectedtour_index=1;
        }
      }
    }

  } 
  fclose(stream);
  if(mxframes_ini!=0&&mxframes_comm==0)mxframes=mxframes_ini;
  if(mxpoints_ini!=0&&mxpoints_comm==0)mxpoints=mxpoints_ini;
  return 0;

}

/* ------------------ writeini ------------------------ */

void writeini(int flag){
  char buffer[1024];

  float *iso_colors_temp;
  int n_iso_colors_temp;
  int n3d;
  int nn;
  extern int numplot3dvars;
  FILE *fileout;
  int i;
  int j;

  fileout=NULL;
  switch (flag) {
  case GLOBAL_INI:
    fileout=fopen(INIfile,"w");
    break;
  case STDOUT_INI:
    fileout=stdout;
    break;
  case SCRIPT_INI:
    fileout=fopen(scriptinifilename,"w");
    break;
  case LOCAL_INI:
    fileout=fopen(caseinifilename,"w");
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  if(flag==SCRIPT_INI)flag=LOCAL_INI;
  if(fileout==NULL){
    printf("error: unable to open %s for writing\n",caseinifilename);
    return;
  }


  fprintf(fileout,"# NIST Smokeview configuration file, Release %s\n\n",__DATE__);
  fprintf(fileout,"COLORS\n");
  fprintf(fileout,"------\n\n");
  fprintf(fileout,"COLORBAR\n");
  fprintf(fileout," %i %i %i\n",nrgb,usetexturebar,colorbar_select_index);
  for(nn=0;nn<nrgb;nn++){
    fprintf(fileout," %f %f %f\n",rgb[nn][0],rgb[nn][1],rgb[nn][2]);
  }
  fprintf(fileout,"COLOR2BAR\n");
  fprintf(fileout," %i\n",8);
  fprintf(fileout," %f %f %f :white  \n",rgb2[0][0],rgb2[0][1],rgb2[0][2]);
  fprintf(fileout," %f %f %f :yellow \n",rgb2[1][0],rgb2[1][1],rgb2[1][2]);
  fprintf(fileout," %f %f %f :blue   \n",rgb2[2][0],rgb2[2][1],rgb2[2][2]);
  fprintf(fileout," %f %f %f :red    \n",rgb2[3][0],rgb2[3][1],rgb2[3][2]);
  fprintf(fileout," %f %f %f :green  \n",rgb2[4][0],rgb2[4][1],rgb2[4][2]);
  fprintf(fileout," %f %f %f :magenta\n",rgb2[5][0],rgb2[5][1],rgb2[5][2]);
  fprintf(fileout," %f %f %f :cyan   \n",rgb2[6][0],rgb2[6][1],rgb2[6][2]);
  fprintf(fileout," %f %f %f :black  \n",rgb2[7][0],rgb2[7][1],rgb2[7][2]);
  if(iso_ambient_ini!=NULL){
    n_iso_colors_temp=n_iso_ambient_ini;
    iso_colors_temp=iso_ambient_ini;
  }
  else{
    n_iso_colors_temp=n_iso_ambient;
    iso_colors_temp=iso_ambient;
  }
  fprintf(fileout,"ISOCOLORS\n");
	fprintf(fileout," %f %f : shininess, transparency\n",iso_shininess, iso_transparency);
	fprintf(fileout," %f %f %f : specular\n",iso_specular[0],iso_specular[1],iso_specular[2]);
  fprintf(fileout," %i\n",n_iso_colors_temp);
  for(nn=0;nn<n_iso_colors_temp;nn++){
    fprintf(fileout," %f %f %f\n",iso_colors_temp[4*nn],iso_colors_temp[4*nn+1],iso_colors_temp[4*nn+2]);
  }
  fprintf(fileout,"VENTCOLOR\n");
  fprintf(fileout," %f %f %f\n",ventcolor[0],ventcolor[1],ventcolor[2]);
  fprintf(fileout,"SENSORCOLOR\n");
  fprintf(fileout," %f %f %f\n",sensorcolor[0],sensorcolor[1],sensorcolor[2]);
  fprintf(fileout,"SENSORNORMCOLOR\n");
  fprintf(fileout," %f %f %f\n",sensornormcolor[0],sensornormcolor[1],sensornormcolor[2]);
  fprintf(fileout,"HEATONCOLOR\n");
  fprintf(fileout," %f %f %f\n",heatoncolor[0],heatoncolor[1],heatoncolor[2]);
  fprintf(fileout,"HEATOFFCOLOR\n");
  fprintf(fileout," %f %f %f\n",heatoffcolor[0],heatoffcolor[1],heatoffcolor[2]);
  fprintf(fileout,"SPRINKONCOLOR\n");
  fprintf(fileout," %f %f %f\n",sprinkoncolor[0],sprinkoncolor[1],sprinkoncolor[2]);
  fprintf(fileout,"SPRINKOFFCOLOR\n");
  fprintf(fileout," %f %f %f\n",sprinkoffcolor[0],sprinkoffcolor[1],sprinkoffcolor[2]);
  fprintf(fileout,"BLOCKCOLOR\n");
  fprintf(fileout," %f %f %f\n",block_ambient2[0],block_ambient2[1],block_ambient2[2]);
  fprintf(fileout,"BOUNDCOLOR\n");
  fprintf(fileout," %f %f %f\n",boundcolor[0],boundcolor[1],boundcolor[2]);
  fprintf(fileout,"STATICPARTCOLOR\n");
  fprintf(fileout," %f %f %f\n",static_color[0],static_color[1],static_color[2]);
  fprintf(fileout,"BACKGROUNDCOLOR\n");
  fprintf(fileout," %f %f %f\n",backgroundbasecolor[0],backgroundbasecolor[1],backgroundbasecolor[2]);
  fprintf(fileout,"FOREGROUNDCOLOR\n");
  fprintf(fileout," %f %f %f\n",foregroundbasecolor[0],foregroundbasecolor[1],foregroundbasecolor[2]);
/*  extern GLfloat iso_ambient[4], iso_specular[4], iso_shininess;*/

  fprintf(fileout,"FLIP\n");
  fprintf(fileout," %i\n",background_flip);
  fprintf(fileout,"TIMEBARCOLOR\n");
  fprintf(fileout," %f %f %f\n",timebarcolor[0],timebarcolor[1],timebarcolor[2]);
  fprintf(fileout,"SETBW\n");
  fprintf(fileout," %i\n",setbw);
  fprintf(fileout,"COLORBARFLIP\n");
  fprintf(fileout," %i\n",colorbarflip);

  fprintf(fileout,"\n LIGHTING\n");
  fprintf(fileout,"--------\n\n");
#ifdef pp_SHOWLIGHT
  fprintf(fileout,"SHOWLIGHT0\n");
  fprintf(fileout," %i\n",visLIGHT0);
  fprintf(fileout,"SHOWLIGHT1\n");
  fprintf(fileout," %i\n",visLIGHT1);
  fprintf(fileout,"SHOWLIGHTMENU\n");
  fprintf(fileout," %i\n",showlightmenu);
#endif
  fprintf(fileout,"AMBIENTLIGHT\n");
  fprintf(fileout," %f %f %f\n",ambientlight[0],ambientlight[1],ambientlight[2]);
  fprintf(fileout,"DIFFUSELIGHT\n");
  fprintf(fileout," %f %f %f\n",diffuselight[0],diffuselight[1],diffuselight[2]);

  fprintf(fileout,"\n SIZES\n");
  fprintf(fileout,"-----------\n\n");
  fprintf(fileout,"VECTORPOINTSIZE\n");
  fprintf(fileout," %f\n",vectorpointsize);
  fprintf(fileout,"VECTORLINEWIDTH\n");
  fprintf(fileout," %f\n",vectorlinewidth);
  fprintf(fileout,"PARTPOINTSIZE\n");
  fprintf(fileout," %f\n",partpointsize);
  fprintf(fileout,"STREAKLINEWIDTH\n");
  fprintf(fileout," %f\n",streaklinewidth);
  fprintf(fileout,"VECTORPOINTSIZE\n");
  fprintf(fileout," %f\n",vectorpointsize);
  fprintf(fileout,"ISOPOINTSIZE\n");
  fprintf(fileout," %f\n",isopointsize);
  fprintf(fileout,"ISOLINEWIDTH\n");
  fprintf(fileout," %f\n",isolinewidth);
  fprintf(fileout,"PLOT3DPOINTSIZE\n");
  fprintf(fileout," %f\n",plot3dpointsize);
  fprintf(fileout,"PLOT3DLINEWIDTH\n");
  fprintf(fileout," %f\n",plot3dlinewidth);
  fprintf(fileout,"VECTORLENGTH\n");
  fprintf(fileout," %f\n",VECFRACTION);
  fprintf(fileout,"SENSORABSSIZE\n");
  fprintf(fileout," %f\n",sensorabssize);
  fprintf(fileout,"SENSORRELSIZE\n");
  fprintf(fileout," %f\n",sensorrelsize);
  fprintf(fileout,"SPRINKLERABSSIZE\n");
  fprintf(fileout," %f\n",sprinklerabssize);
  fprintf(fileout,"SPHERESEGS\n");
  fprintf(fileout," %i\n",device_sphere_segments);
  if(no_graphics==0&&
     (screenWidth == glutGet(GLUT_SCREEN_WIDTH)||screenHeight == glutGet(GLUT_SCREEN_HEIGHT))
    ){
    fprintf(fileout,"WINDOWWIDTH\n");
    fprintf(fileout," %i\n",-1);
    fprintf(fileout,"WINDOWHEIGHT\n");
    fprintf(fileout," %i\n",-1);
  }
  else{
    fprintf(fileout,"WINDOWWIDTH\n");
    fprintf(fileout," %i\n",screenWidth);
    fprintf(fileout,"WINDOWHEIGHT\n");
    fprintf(fileout," %i\n",screenHeight);
  }
  fprintf(fileout,"WINDOWOFFSET\n");
  fprintf(fileout," %i\n",titlesafe_offsetBASE);
  fprintf(fileout,"WINDOWOFFSET\n");
  fprintf(fileout," %i\n",titlesafe_offsetBASE);
  fprintf(fileout,"RENDEROPTION\n");
  fprintf(fileout," %i\n",render_option);
  fprintf(fileout,"\n LINES\n");
  fprintf(fileout,"-----------\n\n");
  fprintf(fileout,"LINEWIDTH\n");
  fprintf(fileout," %f\n",linewidth);
  fprintf(fileout,"VENTLINEWIDTH\n");
  fprintf(fileout," %f\n",ventlinewidth);
  fprintf(fileout,"SMOOTHLINES\n");
  fprintf(fileout," %i\n",antialiasflag);

  fprintf(fileout,"\nOFFSETS\n");
  fprintf(fileout,"-------\n\n");
  fprintf(fileout,"VENTOFFSET\n");
  fprintf(fileout," %f\n",ventoffset_factor);
  fprintf(fileout,"SLICEOFFSET\n");
  fprintf(fileout," %f\n",sliceoffset_factor);

  fprintf(fileout,"\nTIME MIN/MAX\n");
  fprintf(fileout,"------------\n");
  fprintf(fileout,"(0/1 min 0/1 max (1=set, 0=unset)\n\n");
  fprintf(fileout,"T_PARTICLES\n");
  fprintf(fileout," %i %f %i %f\n",settmin_p,tmin_p,settmax_p,tmax_p);
  fprintf(fileout,"T_SLICE\n");
  fprintf(fileout," %i %f %i %f\n",settmin_s,tmin_s,settmax_s,tmax_s);
  fprintf(fileout,"T_ISO\n");
  fprintf(fileout," %i %f %i %f\n",settmin_i,tmin_i,settmax_i,tmax_i);
  fprintf(fileout,"T_BOUNDARY\n");
  fprintf(fileout," %i %f %i %f\n",settmin_b,tmin_b,settmax_b,tmax_b);

  fprintf(fileout,"\nVALUE MIN/MAX\n");
  fprintf(fileout,"-------------\n");
  fprintf(fileout,"(0/1 min 0/1 max (1=set, 0=unset)\n\n");
  if(npart5prop>0){
    for(i=0;i<npart5prop;i++){
      part5prop *propi;

      propi = part5propinfo + i;
      fprintf(fileout,"V5_PARTICLES\n");
      fprintf(fileout," %i %f %i %f %s\n",
        propi->setvalmin,propi->valmin,propi->setvalmax,propi->valmax,propi->label->shortlabel);
    }
  }
  fprintf(fileout,"V_PARTICLES\n");
  fprintf(fileout," %i %f %i %f\n",setpartmin,partmin,setpartmax,partmax);
  fprintf(fileout,"C_PARTICLES\n");
  fprintf(fileout," %i %f %i %f\n",setpartchopmin,partchopmin,setpartchopmax,partchopmax);
  if(nslice2>0){
    for(i=0;i<nslice2;i++){
      fprintf(fileout,"V_SLICE\n");
#ifdef pp_SLICECONTOURS
      fprintf(fileout," %i %f %i %f %s : %f %f %i\n",
#else
      fprintf(fileout," %i %f %i %f %s\n",
#endif
        slicebounds[i].setvalmin,slicebounds[i].valmin,
        slicebounds[i].setvalmax,slicebounds[i].valmax,
        slicebounds[i].label->shortlabel
#ifdef pp_SLICECONTOURS
        ,slicebounds[i].line_contour_min,slicebounds[i].line_contour_max,slicebounds[i].line_contour_num
#endif
        );
    }
    for(i=0;i<nslice2;i++){
      fprintf(fileout,"C_SLICE\n");
      fprintf(fileout," %i %f %i %f %s\n",
        slicebounds[i].setchopmin,slicebounds[i].chopmin,
        slicebounds[i].setchopmax,slicebounds[i].chopmax,
        slicebounds[i].label->shortlabel
        );
    }
    fprintf(fileout,"SLICEDATAOUT\n");
    fprintf(fileout," %i \n",output_slicedata);
  }
  if(niso_bounds>0){
    for(i=0;i<niso_bounds;i++){
      fprintf(fileout,"V_ISO\n");
      fprintf(fileout," %i %f %i %f %s\n",
        isobounds[i].setvalmin,isobounds[i].valmin,
        isobounds[i].setvalmax,isobounds[i].valmax,
        isobounds[i].label->shortlabel
        );
    }
    for(i=0;i<niso_bounds;i++){
      fprintf(fileout,"C_ISO\n");
      fprintf(fileout," %i %f %i %f %s\n",
        isobounds[i].setchopmin,isobounds[i].chopmin,
        isobounds[i].setchopmax,isobounds[i].chopmax,
        isobounds[i].label->shortlabel
        );
    }
  }
  for(i=0;i<npatch_files;i++){
    if(patchinfo[i].firstshort==1){
      fprintf(fileout,"V_BOUNDARY\n");
      fprintf(fileout," %i %f %i %f %s\n",
        patchinfo[i].setvalmin,patchinfo[i].valmin,
        patchinfo[i].setvalmax,patchinfo[i].valmax,
        patchinfo[i].label.shortlabel
        );
    }
  }
  fprintf(fileout,"V_ZONE\n");
  fprintf(fileout," %i %f %i %f\n",setzonemin,zonemin,setzonemax,zonemax);

  fprintf(fileout,"V_PLOT3D\n");
  n3d = 5;
  if(n3d<numplot3dvars)n3d=numplot3dvars;
  if(n3d>mxplot3dvars)n3d=mxplot3dvars;
  fprintf(fileout," %i\n",n3d);
  for(i=0;i<n3d;i++){
    fprintf(fileout," %i %i %f %i %f\n",i+1,setp3min[i],p3min[i],setp3max[i],p3max[i]);
  }
  fprintf(fileout,"C_PLOT3D\n");
  n3d = 5;
  if(n3d<numplot3dvars)n3d=numplot3dvars;
  if(n3d>mxplot3dvars)n3d=mxplot3dvars;
  fprintf(fileout," %i\n",n3d);
  for(i=0;i<n3d;i++){
    fprintf(fileout," %i %i %f %i %f\n",i+1,setp3chopmin[i],p3chopmin[i],setp3chopmax[i],p3chopmax[i]);
  }

  fprintf(fileout,"UNLOAD_QDATA\n");
  fprintf(fileout," %i\n",unload_qdata);

  fprintf(fileout,"V_TARGET\n");
  fprintf(fileout," %i %f %i %f\n",settargetmin,targetmin,settargetmax,targetmax);
  fprintf(fileout,"PERCENTILELEVEL\n");
  fprintf(fileout," %f\n",percentile_level);


  fprintf(fileout,"\nDATA LOADING\n");
  fprintf(fileout,"------------\n\n");
  fprintf(fileout,"MAXFRAMES\n");
  fprintf(fileout," %i\n",mxframes);
  fprintf(fileout,"MAXPOINTS\n");
  fprintf(fileout," %i\n",mxpoints);
  fprintf(fileout,"NOPART\n");
  fprintf(fileout," %i\n",nopart);
  fprintf(fileout,"PARTPOINTSTEP\n");
  fprintf(fileout," %i\n",partpointstep);
  fprintf(fileout,"PARTFRAMESTEP\n");
  fprintf(fileout," %i\n",partframestep);
  fprintf(fileout,"EVACFRAMESTEP\n");
  fprintf(fileout," %i\n",evacframestep);
  fprintf(fileout,"SLICEFRAMESTEP\n");
  fprintf(fileout," %i\n",sliceframestep);
  fprintf(fileout,"SMOKE3DFRAMESTEP\n");
  fprintf(fileout," %i\n",smoke3dframestep);
  fprintf(fileout,"SLICEAVERAGE\n");
  fprintf(fileout," %i %f %i %i\n",slice_average_flag,slice_average_interval,vis_slice_average,slice_turbprop_flag);
  fprintf(fileout,"SMOKE3DZIPSTEP\n");
  fprintf(fileout," %i\n",smoke3dzipstep);
  fprintf(fileout,"ISOZIPSTEP\n");
  fprintf(fileout," %i\n",isozipstep);
  fprintf(fileout,"SLICEZIPSTEP\n");
  fprintf(fileout," %i\n",slicezipstep);
  fprintf(fileout,"BOUNDFRAMESTEP\n");
  fprintf(fileout," %i\n",boundframestep);
  fprintf(fileout,"BOUNDZIPSTEP\n");
  fprintf(fileout," %i\n",boundzipstep);
  fprintf(fileout,"ISOFRAMESTEP\n");
  fprintf(fileout," %i\n",isoframestep);
  fprintf(fileout,"SHOWTRACERSALWAYS\n");
  fprintf(fileout," %i\n",show_tracers_always);
 
  if(flag==LOCAL_INI){
    fprintf(fileout,"AVATAREVAC\n");
    fprintf(fileout," %i\n",iavatar_evac);
  }
 
  if(flag==LOCAL_INI){
    {
      int ndevice_vis=0;
      sv_object *obj_typei;

      for(i=0;i<nobject_defs;i++){
        obj_typei = object_defs[i];
        if(obj_typei->used==1&&obj_typei->visible==1){
          ndevice_vis++;
        }
      }
      fprintf(fileout,"SHOWDEVICES\n");
      fprintf(fileout," %i\n",ndevice_vis);
      for(i=0;i<nobject_defs;i++){
        obj_typei = object_defs[i];
        if(obj_typei->used==1&&obj_typei->visible==1){
          fprintf(fileout," %s\n",obj_typei->label);
        }
      }
    }

    put_startup_smoke3d(fileout);
    fprintf(fileout,"LOADFILESATSTARTUP\n");
    fprintf(fileout," %i\n",loadfiles_at_startup);
    if(npart5prop>0){
      fprintf(fileout,"PART5PROPDISP\n");
      for(i=0;i<npart5prop;i++){
        part5prop *propi;

        propi = part5propinfo + i;
        fprintf(fileout," ");
        for(j=0;j<npartclassinfo;j++){
          fprintf(fileout,"%i ",propi->class_vis[j]);
        }
        fprintf(fileout,"\n");
      }
      fprintf(fileout,"PART5COLOR\n");
      for(i=0;i<npart5prop;i++){
        part5prop *propi;

        propi = part5propinfo + i;
        if(propi->display==1){
          fprintf(fileout," %i\n",i);
          break;
        }
      }

    }
  }
  if(flag==LOCAL_INI&&npartclassinfo>0){
    fprintf(fileout,"PART5CLASSVIS\n");
    fprintf(fileout," %i\n",npartclassinfo);
    for(j=0;j<npartclassinfo;j++){
      part5class *partclassj;

      partclassj = partclassinfo + j;
      fprintf(fileout," %i\n",partclassj->vis_type);
    }
  }
  if(flag==LOCAL_INI&&npropinfo>0){
    fprintf(fileout,"PROPINDEX\n");
    fprintf(fileout," %i\n",npropinfo);
    for(i=0;i<npropinfo;i++){
      propdata *propi;
      int offset;
      int jj;

      propi = propinfo + i;
      offset=-1;
      for(jj=0;jj<propi->nsmokeview_ids;jj++){
        if(strcmp(propi->smokeview_id,propi->smokeview_ids[jj])==0){
          offset=jj;
          break;
        }
      }
      fprintf(fileout," %i %i\n",i,offset);
    }
  }

  fprintf(fileout,"\nCONTOURS\n");
  fprintf(fileout,"--------\n\n");
  fprintf(fileout,"P3CONT2D\n");
  fprintf(fileout," %i\n",p3cont2d);
  fprintf(fileout,"P3VIEW\n");
  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi = meshinfo + i;
    fprintf(fileout," %i %i %i %i %i %i \n",meshi->visx,meshi->plotx,meshi->visy,meshi->ploty,meshi->visz,meshi->plotz);
  }
  for(i=0;i<nmeshes;i++){
    mesh *meshi;

    meshi = meshinfo + i;
    if(meshi->mesh_offset_ptr!=NULL){
      fprintf(fileout,"MESHOFFSET\n");
      fprintf(fileout," %i\n",i);
    }
  }
  fprintf(fileout,"TRANSPARENT\n");
  fprintf(fileout," %i\n",transparentflag);
  fprintf(fileout,"SURFINC\n");
  fprintf(fileout," %i\n",surfincrement);
  fprintf(fileout,"P3DSURFACETYPE\n");
  fprintf(fileout," %i\n",p3dsurfacetype);
  fprintf(fileout,"P3DSURFACESMOOTH\n");
  fprintf(fileout," %i\n",p3dsurfacesmooth);

  fprintf(fileout,"\nVISIBILITY\n");
  fprintf(fileout,"----------\n\n");
  if(nmeshes>1){
    fprintf(fileout,"MESHVIS\n");
    fprintf(fileout," %i\n",nmeshes);

    for(i=0;i<nmeshes;i++){
      mesh *meshi;

      meshi = meshinfo + i;
      fprintf(fileout," %i\n",meshi->blockvis);
    }
  }
  fprintf(fileout,"SHOWTITLE\n");
  fprintf(fileout," %i\n",visTitle0);
  fprintf(fileout,"SHOWCOLORBARS\n");
  fprintf(fileout," %i\n",visColorLabels);
  fprintf(fileout,"SHOWBLOCKS\n");
  fprintf(fileout," %i\n",visBlocks);
  fprintf(fileout,"SHOWNORMALWHENSMOOTH\n");
  fprintf(fileout," %i\n",visSmoothAsNormal);
  fprintf(fileout,"SMOOTHBLOCKSOLID\n");
  fprintf(fileout," %i\n",smooth_block_solid);
  fprintf(fileout,"SBATSTART\n");
  fprintf(fileout," %i\n",sb_atstart);
  fprintf(fileout,"SHOWTRANSPARENT\n");
  fprintf(fileout," %i\n",visTransparentBlockage);
  fprintf(fileout,"SHOWVENTS\n");
  fprintf(fileout," %i %i %i\n",visVents,visVentLines,visVentSolid);
  fprintf(fileout,"SHOWTRANSPARENTVENTS\n");
  fprintf(fileout," %i\n",show_transparent_vents);
  fprintf(fileout,"SHOWSENSORS\n");
  fprintf(fileout," %i %i\n",visSensor,visSensorNorm);
  fprintf(fileout,"SHOWTIMEBAR\n");
  fprintf(fileout," %i\n",visTimeLabels);
  fprintf(fileout,"SHOWTIMELABEL\n");
  fprintf(fileout," %i\n",visTimelabel);
  fprintf(fileout,"SHOWFRAMELABEL\n");
  fprintf(fileout," %i\n",visFramelabel);
  fprintf(fileout,"SHOWFRAMELABEL\n");
  fprintf(fileout," %i\n",visFramelabel);
  fprintf(fileout,"SHOWFLOOR\n");
  fprintf(fileout," %i\n",visFloor);
  fprintf(fileout,"SHOWWALLS\n");
  fprintf(fileout," %i\n",visWalls);
  fprintf(fileout,"SHOWCEILING\n");
  fprintf(fileout," %i\n",visCeiling);
  fprintf(fileout,"SHOWSMOKEPART\n");
  fprintf(fileout," %i\n",visSmokePart);
  fprintf(fileout,"SHOWSPRINKPART\n");
  fprintf(fileout," %i\n",visSprinkPart);
#ifdef pp_memstatus
  fprintf(fileout,"SHOWMEMLOAD\n");
  fprintf(fileout," %i\n",visAvailmemory);
#endif
  fprintf(fileout,"SHOWBLOCKLABEL\n");
  fprintf(fileout," %i\n",visBlocklabel);
  fprintf(fileout,"SHOWAXISLABELS\n");
  fprintf(fileout," %i\n",visaxislabels);
  fprintf(fileout,"SHOWFRAME\n");
  fprintf(fileout," %i\n",visFrame);
  fprintf(fileout,"SHOWALLTEXTURES\n");
  fprintf(fileout," %i\n",showall_textures);
  fprintf(fileout,"SHOWTHRESHOLD\n");
  fprintf(fileout," %i %i %f\n",vis_threshold,vis_onlythreshold,temp_threshold);
  fprintf(fileout,"SHOWHRRCUTOFF\n");
  fprintf(fileout," %i\n",show_hrrcutoff);
  fprintf(fileout,"TWOSIDEDVENTS\n");
  fprintf(fileout," %i %i\n",show_bothsides_int,show_bothsides_ext);
  fprintf(fileout,"TRAINERVIEW\n");
  fprintf(fileout," %i\n",trainerview);
  fprintf(fileout,"SHOWTERRAIN\n");
  fprintf(fileout," %i\n",visTerrainType);
  fprintf(fileout,"TERRAINPARMS\n");
  fprintf(fileout,"%i %i %i\n",terrain_rgba_zmin[0],terrain_rgba_zmin[1],terrain_rgba_zmin[2]);
  fprintf(fileout,"%i %i %i\n",terrain_rgba_zmax[0],terrain_rgba_zmax[1],terrain_rgba_zmax[2]);
  fprintf(fileout,"%f\n",vertical_factor);
  fprintf(fileout,"OFFSETSLICE\n");
  fprintf(fileout," %i\n",offset_slice);
  fprintf(fileout,"SHOWSTREAK\n");
  fprintf(fileout," %i %i %i %i\n",streak5show,streak5step,showstreakhead,streak_index);
  fprintf(fileout,"VECLENGTH\n");
  fprintf(fileout," %i\n",iveclengths);
  fprintf(fileout,"ISOTRAN2\n");
  fprintf(fileout," %i\n",transparent_state);
  fprintf(fileout,"SHOWISO\n");
  fprintf(fileout," %i\n",visAIso);
  fprintf(fileout,"SHOWISONORMALS\n");
  fprintf(fileout," %i\n",showisonormals);
  fprintf(fileout,"SMOKESENSORS\n");
  fprintf(fileout," %i %i\n",show_smokesensors,test_smokesensors);


  fprintf(fileout,"\nMISC\n");
  fprintf(fileout,"----\n\n");
  if(use_nistlogo==1){
  fprintf(fileout,"USENISTLOGO\n");
  fprintf(fileout," %i\n",use_nistlogo);
  }
  if(trainer_mode==1){
    fprintf(fileout,"TRAINERMODE\n");
    fprintf(fileout," %i\n",trainer_mode);
  }
  fprintf(fileout,"SHOWOPENVENTS\n");
  fprintf(fileout," %i %i\n",visOpenVents,visOpenVentsAsOutline);
  fprintf(fileout,"SHOWDUMMYVENTS\n");
  fprintf(fileout," %i\n",visDummyVents);
  fprintf(fileout,"SHOWSLICEINOBST\n");
  fprintf(fileout," %i\n",show_slice_in_obst);
  fprintf(fileout,"SKIPEMBEDSLICE\n");
  fprintf(fileout," %i\n",skip_slice_in_embedded_mesh);
  fprintf(fileout,"CELLCENTERINTERP\n");
  fprintf(fileout," %i\n",cellcenter_interp);
  fprintf(fileout,"SHOWTICKS\n");
  fprintf(fileout," %i\n",visTicks);
  if(flag==LOCAL_INI){
    fprintf(fileout,"USERTICKS\n");
    fprintf(fileout," %i %i %i %i %i %i\n",vis_user_ticks,auto_user_tick_placement,user_tick_sub,
      user_tick_show_x,user_tick_show_y,user_tick_show_z);
    fprintf(fileout," %f %f %f\n",user_tick_origin[0],user_tick_origin[1],user_tick_origin[2]);
    fprintf(fileout," %f %f %f\n",user_tick_min[0],user_tick_min[1],user_tick_min[2]);
    fprintf(fileout," %f %f %f\n",user_tick_max[0],user_tick_max[1],user_tick_max[2]);
    fprintf(fileout," %f %f %f\n",user_tick_step[0],user_tick_step[1],user_tick_step[2]);
    fprintf(fileout," %i %i %i\n",user_tick_show_x,user_tick_show_y,user_tick_show_z);
  }
#ifdef pp_SHOOTER
  if(flag==LOCAL_INI){
    fprintf(fileout,"SHOOTER\n");
    fprintf(fileout," %f %f %f\n",shooter_xyz[0],shooter_xyz[1],shooter_xyz[2]);
    fprintf(fileout," %f %f %f\n",shooter_dxyz[0],shooter_dxyz[1],shooter_dxyz[2]);
    fprintf(fileout," %f %f %f\n",shooter_uvw[0],shooter_uvw[1],shooter_uvw[2]);
    fprintf(fileout," %f %f %f\n",   shooter_velmag, shooter_veldir, shooterpointsize);
    fprintf(fileout," %i %i %i %i %i\n",shooter_fps,shooter_vel_type,shooter_nparts,visShooter,shooter_cont_update);
    fprintf(fileout," %f %f\n",shooter_duration,shooter_v_inf);
  }
#endif
  fprintf(fileout,"SHOWLABELS\n");
  fprintf(fileout," %i\n",visLabels);
  fprintf(fileout,"SHOWFRAMERATE\n");
  fprintf(fileout," %i\n",visFramerate);
  fprintf(fileout,"FRAMERATEVALUE\n");
  fprintf(fileout," %i\n",frameratevalue);
  fprintf(fileout,"VECTORSKIP\n");
  fprintf(fileout," %i\n",vectorskip);
  fprintf(fileout,"AXISSMOOTH\n");
  fprintf(fileout," %i\n",axissmooth);
  fprintf(fileout,"AXISNUM\n");
  fprintf(fileout," %i\n",axisnum);
  fprintf(fileout,"BLOCKLOCATION\n");
  fprintf(fileout," %i\n",blocklocation);
  fprintf(fileout,"SHOWCADANDGRID\n");
  fprintf(fileout," %i\n",show_cad_and_grid);
  fprintf(fileout,"OUTLINEMODE\n");
  fprintf(fileout," %i\n",highlight_flag);
  fprintf(fileout,"TITLESAFE\n");
  fprintf(fileout," %i\n",titlesafe_offset);
  fprintf(fileout,"FONTSIZE\n");
  fprintf(fileout," %i\n",fontindex);
  fprintf(fileout,"ZOOM\n");
  fprintf(fileout," %i %f\n",zoomindex,zoom);
  fprintf(fileout,"APERTURE\n");
  fprintf(fileout," %i\n",apertureindex);
  fprintf(fileout,"RENDERFILETYPE\n");
  fprintf(fileout," %i\n",renderfiletype);
  fprintf(fileout,"SHOWGRIDLOC\n");
  fprintf(fileout," %i\n",visgridloc);
  fprintf(fileout,"PIXELSKIP\n");
  fprintf(fileout," %i\n",pixel_skip);
  fprintf(fileout,"PROJECTION\n");
  fprintf(fileout," %i\n",projection_type);
  fprintf(fileout,"STEREO\n");
  fprintf(fileout," %i\n",showstereo);

  if(nskyboxinfo>0){
    int iskybox;
    skyboxdata *skyi;
    char *filei;
    char *nullfile="NULL";

    for(iskybox=0;iskybox<nskyboxinfo;iskybox++){
      skyi = skyboxinfo + iskybox;
      fprintf(fileout,"SKYBOX\n");
      for(i=0;i<6;i++){
        filei = skyi->face[i].file;
        if(filei==NULL)filei=nullfile;
        if(strcmp(filei,"NULL")==0){
          fprintf(fileout,"NULL\n");
        }
        else{
          fprintf(fileout,"%s\n",filei);
        }
      }
    }
  }

  fprintf(fileout,"UNITCLASSES\n");
  fprintf(fileout," %i\n",nunitclasses);
  for(i=0;i<nunitclasses;i++){
    fprintf(fileout, "%i\n",unitclasses[i].active);
  }
  if(flag==LOCAL_INI){
    fprintf(fileout,"MSCALE\n");
    fprintf(fileout," %f %f %f\n",mscale[0],mscale[1],mscale[2]);
  }
  fprintf(fileout,"CLIP\n");
  fprintf(fileout," %f %f\n",nearclip,farclip);

  if(flag==LOCAL_INI){
    fprintf(fileout,"XYZCLIP\n");
    fprintf(fileout," %i\n",xyz_clipplane);
    fprintf(fileout," %i %f %i %f\n",clip_x, clip_x_val, clip_X, clip_X_val);
    fprintf(fileout," %i %f %i %f\n",clip_y, clip_y_val, clip_Y, clip_Y_val);
    fprintf(fileout," %i %f %i %f\n",clip_z, clip_z_val, clip_Z, clip_Z_val);

    for(i=nlabelssmv;i<nlabels;i++){
      labeldata *labeli;
      float *xyz, *rgbtemp, *tstart_stop;

      labeli = labelinfo + i;
      xyz = labeli->xyz;
      rgbtemp = labeli->rgb;
      tstart_stop = labeli->tstart_stop;

      fprintf(fileout,"LABEL\n");
      fprintf(fileout," %f %f %f %f %f %f %f %f\n",xyz[0],xyz[1],xyz[2],rgbtemp[0],rgbtemp[1],rgbtemp[2],tstart_stop[0],tstart_stop[1]);
      fprintf(fileout,"%s\n",labeli->label);
    }

    
    for(i=ntickssmv;i<nticks;i++){
      float *begt;
      float *endt;
      float *rgbtemp;
      tickdata *ticki;

      ticki = tickinfo + i;
      begt = ticki->begin;
      endt = ticki->end;
      rgbtemp = ticki->rgb;

      fprintf(fileout,"TICKS\n");
      fprintf(fileout," %f %f %f %f %f %f %i\n",begt[0],begt[1],begt[2],endt[0],endt[1],endt[2],ticki->nbars);
      fprintf(fileout," %f %i %f %f %f %f\n",ticki->dlength,ticki->dir,rgbtemp[0],rgbtemp[1],rgbtemp[2],ticki->width);
    }
    
  }
  if(fds_filein!=NULL&&strlen(fds_filein)>0){
    fprintf(fileout,"INPUT_FILE\n");
    fprintf(fileout,"%s\n",fds_filein);
  }

  fprintf(fileout,"EYEX\n");
  fprintf(fileout," %f\n",eyexfactor);
  fprintf(fileout,"EYEY\n");
  fprintf(fileout," %f\n",eyeyfactor);
  fprintf(fileout,"EYEZ\n");
  fprintf(fileout," %f\n",eyezfactor);
  fprintf(fileout,"EYEVIEW\n");
  fprintf(fileout," %i\n",eyeview);
  {
    char *label;

    label = get_camera_label(startup_view_ini);
    if(label!=NULL){
      fprintf(fileout,"LABELSTARTUPVIEW\n");
      fprintf(fileout," %s\n",label);
    }
  }
  fprintf(fileout,"VIEWTIMES\n");
  fprintf(fileout," %f %f %i\n",view_tstart,view_tstop,view_ntimes);
  fprintf(fileout,"TIMEOFFSET\n");
  fprintf(fileout," %f\n",timeoffset);
  fprintf(fileout,"SHOWHMSTIMELABEL\n");
  fprintf(fileout," %i\n",vishmsTimelabel);
  fprintf(fileout,"SPEED\n");
  fprintf(fileout," %f %f\n",speed_crawl,speed_walk);

  fprintf(fileout,"CULLFACES\n");
  fprintf(fileout," %i\n",cullfaces);
  fprintf(fileout,"\nZone\n");
  fprintf(fileout,"----\n\n");
  fprintf(fileout,"SHOWHZONE\n");
  fprintf(fileout," %i\n",visHZone);
  fprintf(fileout,"SHOWVZONE\n");
  fprintf(fileout," %i\n",visVZone);
  fprintf(fileout,"SHOWHAZARDCOLORS\n");
  fprintf(fileout," %i\n",sethazardcolor);
  if(
    ((INI_fds_filein!=NULL&&fds_filein!=NULL&&strcmp(INI_fds_filein,fds_filein)==0)||
    flag==LOCAL_INI)){
    {
      float *eye, *angle_zx, *mat;
      camera *ca;

      for(ca=camera_list_first.next;ca->next!=NULL;ca=ca->next){
        if(strcmp(ca->name,"internal")==0)continue;
        if(strcmp(ca->name,"external")==0)continue;
        fprintf(fileout,"VIEWPOINT5\n");
        eye = ca->eye;
        angle_zx = ca->angle_zx;
        mat = ca->modelview;

		    fprintf(fileout," %i %i %i\n",
          ca->eyeview,
          ca->rotation_index,
          ca->view_id);
		    fprintf(fileout," %f %f %f %f %i\n",
          eye[0],eye[1],eye[2],
          zoom,zoomindex);
  		  fprintf(fileout," %f %f %f %i\n",
          ca->view_angle, 
          ca->direction_angle,
          ca->elevation_angle,
          ca->projection_type);
		    fprintf(fileout," %f %f %f\n",
          ca->xcen,
          ca->ycen,
          ca->zcen);

        fprintf(fileout," %f %f\n",angle_zx[0],angle_zx[1]);
        fprintf(fileout," %f %f %f %f\n",mat[0],mat[1],mat[2],mat[3]);
        fprintf(fileout," %f %f %f %f\n",mat[4],mat[5],mat[6],mat[7]);
        fprintf(fileout," %f %f %f %f\n",mat[8],mat[9],mat[10],mat[11]);
        fprintf(fileout," %f %f %f %f\n",mat[12],mat[13],mat[14],mat[15]);
        fprintf(fileout," %i %i %i %i %i %i %i\n",
            ca->xyz_clipplane,
            ca->clip_x,ca->clip_y,ca->clip_z,
            ca->clip_X,ca->clip_Y,ca->clip_Z);
        fprintf(fileout," %f %f %f %f %f %f\n",
            ca->clip_x_val,ca->clip_y_val,ca->clip_z_val,
            ca->clip_X_val,ca->clip_Y_val,ca->clip_Z_val);
        fprintf(fileout,"%s\n",ca->name);
      }
    }
  }

  fprintf(fileout,"\n3D SMOKE INFO\n");
  fprintf(fileout,"-------------\n\n");
  fprintf(fileout,"ADJUSTALPHA\n");
  fprintf(fileout," %i\n",adjustalphaflag);
#ifdef pp_GPU
  fprintf(fileout,"USEGPU\n");
  fprintf(fileout," %i\n",usegpu);
#endif
  fprintf(fileout,"SMOKECULL\n");
#ifdef pp_CULL
  fprintf(fileout," %i\n",cullsmoke);
#else
  fprintf(fileout," %i\n",smokecullflag);
#endif
  fprintf(fileout,"SMOKESKIP\n");
  fprintf(fileout," %i\n",smokeskipm1);
  fprintf(fileout,"SMOKESHADE\n");
  fprintf(fileout," %i\n",smoke_shade);
#ifdef pp_GPU
  fprintf(fileout,"SMOKERTHICK\n");
  fprintf(fileout," %f\n",smoke3d_rthick);
#else
  fprintf(fileout,"SMOKETHICK\n");
  fprintf(fileout," %i\n",smoke3d_thick);
#endif
  fprintf(fileout,"FIRECOLOR\n");
  fprintf(fileout," %i %i %i\n",fire_red,fire_green,fire_blue);
  fprintf(fileout,"FIREDEPTH\n");
  fprintf(fileout," %f\n",fire_halfdepth);
#ifdef pp_LIGHT
  fprintf(fileout,"SMOKELIGHTING\n");
  fprintf(fileout," %i\n",show_smokelighting);
#endif

  fprintf(fileout,"COLORBARTYPE\n");
  fprintf(fileout," %i\n",colorbartype);
  fprintf(fileout,"SHOWEXTREMEDATA\n");
  fprintf(fileout," %i\n",show_extremedata);
  {
    int mmin[3],mmax[3];
    for(i=0;i<3;i++){
      mmin[i]=rgb_below_min[i];
      mmax[i]=rgb_above_max[i];
    }
    fprintf(fileout,"EXTREMECOLORS\n");
    fprintf(fileout," %i %i %i %i %i %i\n",
     mmin[0],mmin[1],mmin[2],
     mmax[0],mmax[1],mmax[2]);
  }
  if(ncolorbars>ndefaultcolorbars){
	  colorbardata *cbi;
	  unsigned char *rrgb;
	  int n;

    fprintf(fileout,"GCOLORBAR\n");
    fprintf(fileout," %i\n",ncolorbars-ndefaultcolorbars);
    for(n=ndefaultcolorbars;n<ncolorbars;n++){
      cbi = colorbarinfo + n;
      fprintf(fileout,"%s\n",cbi->label);
      fprintf(fileout," %i %i\n",cbi->nnodes,cbi->nodehilight);
      for(i=0;i<cbi->nnodes;i++){
        rrgb = cbi->rgb_node+3*i;
        fprintf(fileout," %i %i %i %i\n",cbi->index_node[i],(int)rrgb[0],(int)rrgb[1],(int)rrgb[2]);
      }
    }
  }
  fprintf(fileout,"\nTOUR INFO\n");
  fprintf(fileout,"---------\n\n");
  fprintf(fileout,"VIEWTOURFROMPATH\n");
  fprintf(fileout," %i\n",viewtourfrompath);
  fprintf(fileout,"VIEWALLTOURS\n");
  fprintf(fileout," %i\n",viewalltours);
  fprintf(fileout,"SHOWTOURROUTE\n");
  fprintf(fileout," %i\n",edittour);
  fprintf(fileout,"SHOWPATHNODES\n");
  fprintf(fileout," %i\n",show_path_knots);
  fprintf(fileout,"TOURCONSTANTVEL\n");
  fprintf(fileout," %i\n",tour_constant_vel);
//  fprintf(fileout,"TOUR_AVATAR\n");
//  fprintf(fileout," %i %f %f %f %f\n",
//    tourlocus_type,
//    tourcol_avatar[0],tourcol_avatar[1],tourcol_avatar[2],
//    tourrad_avatar);
  {
    keyframe *framei;
    float *col;
    int startup_count=0;
    int ii,uselocalspeed=0;

    fprintf(fileout,"TOURCOLORS\n");
    col=tourcol_selectedpathline;
    fprintf(fileout," %f %f %f   :selected path line\n",col[0],col[1],col[2]);
    col=tourcol_selectedpathlineknots;
    fprintf(fileout," %f %f %f   :selected path line knots\n",col[0],col[1],col[2]);
    col=tourcol_selectedknot;
    fprintf(fileout," %f %f %f   :selected knot\n",col[0],col[1],col[2]);
    col=tourcol_pathline;
    fprintf(fileout," %f %f %f   :path line\n",col[0],col[1],col[2]);
    col=tourcol_pathknots;
    fprintf(fileout," %f %f %f   :path knots\n",col[0],col[1],col[2]);
    col=tourcol_text;
    fprintf(fileout," %f %f %f   :text\n",col[0],col[1],col[2]);
    col=tourcol_avatar;
    fprintf(fileout," %f %f %f   :avatar\n",col[0],col[1],col[2]);

    if(flag==LOCAL_INI){
      fprintf(fileout,"TOURINDEX\n");
      fprintf(fileout," %i\n",selectedtour_index);
      startup_count=0;
      for(i=0;i<ntours;i++){
        tourdata *touri;

        touri = tourinfo + i;
        if(touri->startup==1)startup_count++;
      }
      for(ii=0;ii<2;ii++){
        if(ii==0&&ntours-startup_count>0){
          fprintf(fileout,"TOURS\n");
          fprintf(fileout," %i\n",ntours-startup_count);
        }
        else if(ii==1&&startup_count>0){
          fprintf(fileout,"DEFAULTTOURS\n");
          fprintf(fileout," %i\n",startup_count);
        }
        else{
          continue;
        }
        for(i=0;i<ntours;i++){
          tourdata *touri;

          touri = tourinfo + i;
          if(ii==1&&touri->startup==0)continue;
          if(ii==0&&touri->startup==1)continue;

          trim(touri->label);
          fprintf(fileout,"%s\n",touri->label);
          fprintf(fileout," %i %i %f %i %i\n",
            touri->nkeyframes,touri->global_tension_flag,touri->global_tension,touri->glui_avatar_index,touri->display);
          for(framei=&touri->first_frame;framei!=&touri->last_frame;framei=framei->next){
            if(framei==&touri->first_frame)continue;
            sprintf(buffer,"%f %f %f %f ",
              framei->noncon_time,
              xbar0+framei->nodeval.eye[0]*xyzmaxdiff,
              ybar0+framei->nodeval.eye[1]*xyzmaxdiff,
              zbar0+framei->nodeval.eye[2]*xyzmaxdiff);
            trimmzeros(buffer);
            fprintf(fileout,"%s %i ",buffer,framei->viewtype);
            if(framei->viewtype==0){
              sprintf(buffer,"%f %f %f %f %f %f %f ",
                framei->az_path,framei->nodeval.elev_path,framei->bank,
                framei->tension, framei->bias, framei->continuity,
                framei->nodeval.zoom);
            }
            else{
              sprintf(buffer,"%f %f %f %f %f %f %f ",
                xbar0+framei->nodeval.aview[0]*xyzmaxdiff,
                ybar0+framei->nodeval.aview[1]*xyzmaxdiff,
                zbar0+framei->nodeval.aview[2]*xyzmaxdiff,
                framei->tension, framei->bias, framei->continuity,
                framei->nodeval.zoom);
            }
            trimmzeros(buffer);
            fprintf(fileout,"%s %i\n",buffer,uselocalspeed);
          }
        }
      }
    }
  }
  
  if(flag==LOCAL_INI){
    scriptfiledata *scriptfile;
    inifiledata *inifile;
    int first_time=1;
    FILE *fileout_inilist=NULL;

    for(scriptfile=first_scriptfile.next;scriptfile->next!=NULL;scriptfile=scriptfile->next){
      char *file;

      file=scriptfile->file;
      if(file!=NULL){
        fprintf(fileout,"SCRIPTFILE\n");
        fprintf(fileout," %s\n",file);
      }
    }
    for(inifile=first_inifile.next;inifile->next!=NULL;inifile=inifile->next){
      char ini_listfile[1024];

      if(inifile->file==NULL)continue;
      if(first_time==1){
        first_time=0;
        strcpy(ini_listfile,caseinifilename);
        strcat(ini_listfile,"list");
        fileout_inilist=fopen(ini_listfile,"w");
      }
      if(fileout_inilist!=NULL){
        fprintf(fileout_inilist,"INIFILE\n");
        fprintf(fileout_inilist," %s\n",inifile->file);
      }
    }
    if(fileout_inilist!=NULL)fclose(fileout_inilist);
  }
  {
    int svn_num;

    svn_num=getmaxrevision();    // get svn revision number
    fprintf(fileout,"\n\n# Development Environment\n");
    fprintf(fileout,"# -----------------------\n\n");
    fprintf(fileout,"# Smokeview Version: %s\n",SMVVERSION);
    fprintf(fileout,"# Smokeview Revision Number: %i\n",svn_num);
    fprintf(fileout,"# Smokeview Compile Date: %s\n",__DATE__);
    if(no_graphics==0){
      char version_label[256];
      char *glversion=NULL;

      glversion=(char *)glGetString(GL_VERSION);
      if(glversion!=NULL){
        strcpy(version_label,"OpenGL Version: "); 
        strcat(version_label,glversion);
        fprintf(fileout,"# %s\n",version_label);
      }
    }
    if(revision_fds>0){
      fprintf(fileout,"# FDS Revision Number: %i\n",revision_fds);
    }
#ifdef X64
    fprintf(fileout,"# Platform: WIN64\n");
#endif
#ifdef WIN32
#ifndef X64
    fprintf(fileout,"# Platform: WIN32\n");
#endif
#endif
#ifndef pp_OSX64
#ifdef pp_OSX
    printf("Platform: OSX\n");
#endif
#endif
#ifdef pp_OSX64
    printf("Platform: OSX64\n");
#endif
#ifndef pp_LINUX64
#ifdef pp_LINUX
    fprintf(fileout,"# Platform: LINUX\n");
#endif
#endif
#ifdef pp_LINUX64
    fprintf(fileout,"# Platform: LINUX64\n");
#endif

    if(no_graphics==0){
      GLint nred, ngreen, nblue, ndepth, nalpha;

      glGetIntegerv(GL_RED_BITS,&nred);    
      glGetIntegerv(GL_GREEN_BITS,&ngreen);
      glGetIntegerv(GL_BLUE_BITS,&nblue); 
      glGetIntegerv(GL_DEPTH_BITS,&ndepth);
      glGetIntegerv(GL_ALPHA_BITS,&nalpha);
      fprintf(fileout,"\n\n# Graphics Environment\n");
      fprintf(fileout,"# --------------------\n\n");
      fprintf(fileout,"#   Red bits:%i\n",nred);
      fprintf(fileout,"# Green bits:%i\n",ngreen);
      fprintf(fileout,"#  Blue bits:%i\n",nblue);
      fprintf(fileout,"# Alpha bits:%i\n",nalpha);
      fprintf(fileout,"# Depth bits:%i\n\n",ndepth);
    }
  }

  if(fileout!=stdout)fclose(fileout);
}

/* ------------------ update_loaded_lists ------------------------ */

void update_loaded_lists(void){
  int i;
  slice *slicei;
  patch *patchi;

  nslice_loaded=0;
  for(i=0;i<nslice_files;i++){
    slicei = sliceinfo + i;
    if(slicei->loaded==1){
      slice_loaded_list[nslice_loaded]=i;
      nslice_loaded++;
    }
  }

  npatch_loaded=0;
  for(i=0;i<npatch_files;i++){
    patchi = patchinfo + i;
    if(patchi->loaded==1){
      patch_loaded_list[npatch_loaded]=i;
      npatch_loaded++;
    }
  }

}

/* ------------------ get_elevaz ------------------------ */

void get_elevaz(float *xyznorm,float *dtheta,float *rotate_axis){
  float norm2;
  float pi;
  float az, elev;

  // cos(dtheta) = (xyznorm .dot. vec3(0,0,1))/||xyznorm||
  // rotate_axis = xyznorm .cross. vec3(0,0,1)

  pi=4.0*atan(1.0);
  normalize(xyznorm,3);
  *dtheta = 180.0*acos(xyznorm[2])/pi;
  rotate_axis[0]=-xyznorm[1];
  rotate_axis[1]= xyznorm[0];
  rotate_axis[2]=0.0;
  normalize(rotate_axis,2);
}

/* ------------------ getfile_modtime ------------------------ */

void getfile_modtime(char *filename, time_t *modtime){
  STRUCTSTAT statbuffer;
  int statfile;

  *modtime=0;
  if(filename==NULL)return;
  statfile=STAT(filename,&statbuffer);
  if(statfile!=0)return;
  *modtime = statbuffer.st_mtime;
  return;
}

/* ------------------ getfile_size ------------------------ */

void getfile_size(const char *filename, FILE_SIZE *filesize){
  STRUCTSTAT statbuffer;
  int statfile;

  *filesize=0;
  if(filename==NULL)return;
  statfile=STAT(filename,&statbuffer);
  if(statfile!=0)return;
  *filesize = statbuffer.st_size;
  return;
}

/* ------------------ file_exit ------------------------ */

int file_exist(char *file){
  STRUCTSTAT statbuffer;
  int statfile;

  if(file==NULL)return 0;
  statfile=STAT(file,&statbuffer);
  if(statfile!=0)return 0;
  return 1;
}

/* ------------------ get_labels ------------------------ */

void get_labels(char *buffer, char **label1, char **label2){
  char *tok0, *tok1, *tok2;

  tok0=NULL;
  tok1=NULL;
  tok2=NULL;
  tok0 = strtok(buffer,"%");
  if(tok0!=NULL)tok1=strtok(NULL,"%");
  if(tok1!=NULL)tok2=strtok(NULL,"%");
  if(tok1!=NULL){
    trim(tok1);
    tok1=trim_front(tok1);
    if(strlen(tok1)==0)tok1=NULL;
  }
  if(tok2!=NULL){
    trim(tok2);
    tok2=trim_front(tok2);
    if(strlen(tok2)==0)tok2=NULL;
  }
  if(label1!=NULL)*label1=tok1;
  if(label2!=NULL)*label2=tok2;
}

/* ------------------ get_prop_id ------------------------ */

propdata *get_prop_id(char *prop_id){
  int i;

  if(prop_id==NULL||strlen(prop_id)==0)return NULL;
  for(i=0;i<npropinfo;i++){
    propdata *propi;

    propi = propinfo + i;

    if(strcmp(propi->label,prop_id)==0)return propi;
  }
  return NULL;
}
