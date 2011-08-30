// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#ifdef pp_GPU
#include "glew.h"
#endif
#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "contourdefs.h"
#include "isodefs.h"

#include "flowfiles.h"
#include "smokeviewapi.h"
#include "MALLOC.h"
#include "smokeviewvars.h"
#include "translate.h"
#include "update.h"

/* dummy change to bump revision number to 5.1.5 */

#ifdef WIN32
#include <Commdlg.h>
#include <direct.h>
#endif

// svn revision character string
char smokeview_revision[]="$Revision$";
int can_write_to_dir(char *dir);

/* ----------------------- unsetClipPlanes ----------------------------- */

void unsetClipPlanes(void){
  if(xyz_clipplane==2){
    glDisable(GL_CLIP_PLANE0);
    glDisable(GL_CLIP_PLANE1);
    glDisable(GL_CLIP_PLANE2);
    glDisable(GL_CLIP_PLANE3);
    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);
  }
}

/* ----------------------- setClipPlanes ----------------------------- */

void setClipPlanes(int mode){
  static GLdouble clipplane_x[4], clipplane_y[4], clipplane_z[4];
  static GLdouble clipplane_X[4], clipplane_Y[4], clipplane_Z[4];

  if(mode==0&&xyz_clipplane==2)return;
  if(mode==1&&xyz_clipplane!=2)return;
  if(xyz_clipplane==0){
    glDisable(GL_CLIP_PLANE0);
    glDisable(GL_CLIP_PLANE1);
    glDisable(GL_CLIP_PLANE2);
    glDisable(GL_CLIP_PLANE3);
    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);
    return;
  }

  if(clip_x==1){
    clipplane_x[0]=1.0;
    clipplane_x[1]=0.0;
    clipplane_x[2]=0.0;
    clipplane_x[3]=-(clip_x_val-xbar0)/xyzmaxdiff;
    glClipPlane(GL_CLIP_PLANE0,clipplane_x);
    glEnable(GL_CLIP_PLANE0);
  }
  else{
    glDisable(GL_CLIP_PLANE0);
  }

  if(clip_X==1){
    clipplane_X[0]=-1.0;
    clipplane_X[1]=0.0;
    clipplane_X[2]=0.0;
    clipplane_X[3]=(clip_X_val-xbar0)/xyzmaxdiff;
    glClipPlane(GL_CLIP_PLANE3,clipplane_X);
    glEnable(GL_CLIP_PLANE3);
  }
  else{
    glDisable(GL_CLIP_PLANE3);
  }

  if(clip_y==1){
    clipplane_y[0]=0.0;
    clipplane_y[1]=1.0;
    clipplane_y[2]=0.0;
    clipplane_y[3]=-(clip_y_val-ybar0)/xyzmaxdiff;
    glClipPlane(GL_CLIP_PLANE1,clipplane_y);
    glEnable(GL_CLIP_PLANE1);
  }
  else{
    glDisable(GL_CLIP_PLANE1);
  }

  if(clip_Y==1){
    clipplane_Y[0]=0.0;
    clipplane_Y[1]=-1.0;
    clipplane_Y[2]=0.0;
    clipplane_Y[3]=(clip_Y_val-ybar0)/xyzmaxdiff;
    glClipPlane(GL_CLIP_PLANE4,clipplane_Y);
    glEnable(GL_CLIP_PLANE4);
  }
  else{
    glDisable(GL_CLIP_PLANE4);
  }

  if(clip_z==1){
    clipplane_z[0]=0.0;
    clipplane_z[1]=0.0;
    clipplane_z[2]=1.0;
    clipplane_z[3]=-(clip_z_val-zbar0)/xyzmaxdiff;
    glClipPlane(GL_CLIP_PLANE2,clipplane_z);
    glEnable(GL_CLIP_PLANE2);
  }
  else{
    glDisable(GL_CLIP_PLANE2);
  }

  if(clip_Z==1){
    clipplane_Z[0]=0.0;
    clipplane_Z[1]=0.0;
    clipplane_Z[2]=-1.0;
    clipplane_Z[3]=(clip_Z_val-zbar0)/xyzmaxdiff;
    glClipPlane(GL_CLIP_PLANE5,clipplane_Z);
    glEnable(GL_CLIP_PLANE5);
  }
  else{
    glDisable(GL_CLIP_PLANE5);
  }
}

/* ------------------ update_rotation_index ------------------------ */

void update_rotation_index(int val){
  mesh *meshi;
  int i;
  float *modelview_rotate;
  float *angle_zx;
  int *rotation_index;

  rotation_index = &camera_current->rotation_index;

  *rotation_index=val;
  if(*rotation_index==rotation_index_OLD)return;
  if(*rotation_index>=0&&*rotation_index<nmeshes){
    meshi = meshinfo + *rotation_index;
    camera_current->xcen=meshi->xcen;
    camera_current->ycen=meshi->ycen;
    camera_current->zcen=meshi->zcen;
  }
  else{
    camera_current->xcen=xcenGLOBAL;
    camera_current->ycen=ycenGLOBAL;
    camera_current->zcen=zcenGLOBAL;
  }
  rotation_index_OLD=*rotation_index;
  modelview_rotate = camera_current->modelview;
  for(i=0;i<16;i++){
    modelview_rotate[i]=modelview_current[i];
  }

  angle_zx = camera_current->angle_zx;

  angle_zx[0]=0.; 
  angle_zx[1]=0.; 

  camera_current->direction_angle=0.0;
  camera_current->cos_direction_angle = 1.0;
  camera_current->sin_direction_angle = 0.0;

  camera_current->view_angle=0.0;
  camera_current->cos_view_angle = 1.0;
  camera_current->sin_view_angle = 0.0;

  update_meshlist1(val);

  glutPostRedisplay();

}

/* ------------------ snifferrors ------------------------ */

void _sniffErrors(char *whereat){
  int error;
  char *glu_error;
  while((error=glGetError())!=GL_NO_ERROR){
    glu_error=(char *)gluErrorString((unsigned int)error);
    fprintf(stderr,"GL Error:%s, where:%s %i\n",
      glu_error,whereat,snifferrornumber);
      snifferrornumber++;
  }
}

/* ------------------ updateLights ------------------------ */

void updateLights(int pos){
  GLfloat ambientlight2[4], diffuselight2[4];
  int i;

  if(visLIGHT0==1&&visLIGHT1==1){
    for(i=0;i<3;i++){
      ambientlight2[i]=ambientlight[i]/2.0;
      diffuselight2[i]=diffuselight[i]/2.0;
    }
    ambientlight2[3]=1.0;
    diffuselight2[3]=1.0;
    glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuselight2);
    glLightfv(GL_LIGHT0,GL_AMBIENT,ambientlight2);
    if(pos==1)glLightfv(GL_LIGHT0,GL_POSITION,light_position0);
    glEnable(GL_LIGHT0);

    glLightfv(GL_LIGHT1,GL_DIFFUSE,diffuselight2);
    glLightfv(GL_LIGHT1,GL_AMBIENT,ambientlight2);
    if(pos==1)glLightfv(GL_LIGHT1,GL_POSITION,light_position1);
    glEnable(GL_LIGHT1);

  }
  if(visLIGHT0==1&&visLIGHT1==0){
    diffuselight[3]=1.0;
    ambientlight[3]=1.0;
    glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuselight);
    glLightfv(GL_LIGHT0,GL_AMBIENT,ambientlight);
    if(pos==1)glLightfv(GL_LIGHT0,GL_POSITION,light_position0);
    glEnable(GL_LIGHT0);

    glDisable(GL_LIGHT1);

  }
  if(visLIGHT0==0&&visLIGHT1==1){
    diffuselight[3]=1.0;
    ambientlight[3]=1.0;
    glEnable(GL_LIGHT0);
    glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuselight);
    glLightfv(GL_LIGHT0,GL_AMBIENT,ambientlight);
    if(pos==1)glLightfv(GL_LIGHT0,GL_POSITION,light_position1);
    glDisable(GL_LIGHT1);
  }
  if(visLIGHT0==0&&visLIGHT1==0){
    glDisable(GL_LIGHT0);
    glDisable(GL_LIGHT1);
  }


  UpdateLIGHTS=0;

}

/* ------------------ antialias ------------------------ */

void antialias(int flag){
  if(antialiasflag==1){
    if(flag==1){
      glEnable(GL_LINE_SMOOTH);
      glEnable(GL_BLEND);
      glEnable(GL_POINT_SMOOTH);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glHint(GL_LINE_SMOOTH_HINT,GL_DONT_CARE);
    }
    if(flag==0){
      glDisable(GL_LINE_SMOOTH);
      glDisable(GL_POINT_SMOOTH);
      glDisable(GL_BLEND);
    }
  }
}

/* ------------------ transparenton ------------------------ */

void transparenton(void){
  glEnable(GL_BLEND);
  glDepthMask(GL_FALSE);
  glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
}

/* ------------------ transparentoff ------------------------ */

void transparentoff(void){
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
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

/* ------------------ drawTimeBar ------------------------ */

void drawTimeBar(void){
  float xleft=.175f, xright=1.0f, ybot=0.10f, ytop=.35f, xxright;

  glDisable(GL_LIGHTING);
  xleft = xtimeleft;
  if(fontindex==LARGE_FONT)xleft=xtimeleft+0.11;
  xright = xtimeright;

  glLineWidth(linewidth);
  glBegin(GL_LINE_LOOP);
  glColor4fv(timebarcolor);
  glVertex2f(xleft,ybot);
  glVertex2f(xright,ybot);
  glVertex2f(xright,ytop);
  glVertex2f(xleft,ytop);
  glEnd();

  if(ntimes != 1){
    xxright = xleft + (float)itimes*(xright-xleft)/(ntimes-1);
  }
  else{
    xxright=xright;
  }
  glBegin(GL_POLYGON);
  glColor4fv(timebarcolor);
  glVertex2f(xleft,ybot);
  glVertex2f(xxright,ybot);
  glVertex2f(xxright,ytop);
  glVertex2f(xleft,ytop);
  glEnd();

}

/* ------------------ setsmokeviewvars ------------------------ */

void setsmokeviewvars(void){
  int i;

  for(i=0;i<20;i++){
    cputimes[i]=0.0;
  }
  eyexINI=0.0; 
  eyeyINI=0.0;
  eyezINI=0.0;
  for(i=0;i<16;i++){
    modelview_setup[i]=0.0;
  }
  for(i=0;i<4;i++){
    modelview_setup[i+4*i]=1.0;
  }
}

/* ------------------ Init ------------------------ */

void Init(void){
  int errorcode;
  mesh *meshi;

  int n,i;

  FREEMEMORY(plotiso);
  NewMemory((void **)&plotiso,mxplot3dvars*sizeof(int));

  for(n=0;n<mxplot3dvars;n++){
    plotiso[n]=nrgb/2;
  }

  for(i=0;i<16;i++){
    modelview_setup[i]=0.0;
  }
  for(i=0;i<4;i++){
    modelview_setup[i+4*i]=1.0;
  }

  for(i=0;i<nmeshes;i++){
    meshi=meshinfo+i;
    initcontour(&meshi->plot3dcontour1,rgb_plot3d_contour,nrgb);
    initcontour(&meshi->plot3dcontour2,rgb_plot3d_contour,nrgb);
    initcontour(&meshi->plot3dcontour3,rgb_plot3d_contour,nrgb);
  }

  for(i=0;i<nmeshes;i++){
    meshi=meshinfo+i;
    meshi->currentsurf.defined=0;
    meshi->currentsurf2.defined=0;
  }

  /* initialize box sizes, lighting parameters */

  xyzbox = xbar;
  if(ybar>xyzbox){xyzbox=ybar;}
  if(zbar>xyzbox){xyzbox=zbar;}

  {
    char name_external[32];
    strcpy(name_external,"external");
    init_camera(camera_external,name_external);
    camera_external->view_id=1;

    copy_camera(camera_external_save,camera_external);
  }
  if(camera_ini->defined==1){
    copy_camera(camera_current,camera_ini);
  }
  else{
    camera_external->zoom=zoom;
    copy_camera(camera_current,camera_external);
  }
  strcpy(camera_label,camera_current->name);
  update_camera_label();
  {
    char name_internal[32];
    strcpy(name_internal,"internal");
    init_camera(camera_internal,name_internal);
  }
  camera_internal->eye[0]=0.5*xbar;
  camera_internal->eye[1]=0.5*ybar;
  camera_internal->eye[2]=0.5*zbar;
  camera_internal->view_id=0;
  copy_camera(camera_save,camera_current);
  copy_camera(camera_last,camera_current);

  init_camera_list();
  add_default_views();
  update_view_gluilist();

  //reset_glui_view(i_view_list);

  screenWidth2 = screenWidth - dwinW;
  screenHeight2 = screenHeight - dwinH;

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_NORMALIZE);
  if(cullfaces==1)glEnable(GL_CULL_FACE);

  glClearColor(backgroundcolor[0],backgroundcolor[1],backgroundcolor[2], 0.0f);
  glShadeModel(GL_SMOOTH); 
  glDisable(GL_DITHER);


  thistime=0;
  lasttime=0;

  /* define color bar */

  updatecolors(-1);

  block_ambient2[3] = 1.0;
  mat_ambient2[3] = 1.0;
  mat_specular2[3] = 1.0;

  reset_glui_view(startup_view_ini);
  updateShow();
}

/* ------------------ ResetView ------------------------ */

void ResetView(int option){
  int eyeview_save;

  switch (option){
  case RESTORE_EXTERIOR_VIEW_ZOOM:
  case RESTORE_EXTERIOR_VIEW:
    eyeview_save = camera_current->eyeview;
    copy_camera(camera_current,camera_external);
    camera_current->eyeview=eyeview_save;
    if(camera_current->projection_type==1){
      camera_current->eye[1]=camera_current->isometric_y;
    }
    break;
  case RESTORE_INTERIOR_VIEW:
    eyeview_save = camera_current->eyeview;
    copy_camera(camera_current,camera_internal);
    camera_current->eyeview=eyeview_save;
    break;
  case RESTORE_SAVED_VIEW:
    copy_camera(camera_current,camera_save);
    break;
  case 3:
  case 4:
    ASSERT(FFALSE);
    break;
  default:
    ASSERT(FFALSE);
    break;
  }
  if(option==RESTORE_EXTERIOR_VIEW_ZOOM)camera_current->zoom=zooms[zoomindex];
  zoom=camera_current->zoom;
  update_glui_zoom();
}

/* ------------------ UpdateTimeLabels ------------------------ */

void UpdateTimeLabels(void){
  float time0;
  int hour, min, sec,sec10;
 
  time0 = timeoffset;
  if(times!=NULL)time0 = timeoffset + times[itimes];
  if(vishmsTimelabel==1){
    hour = time0/3600;
    min = time0/60.0 - 60*hour;
    sec10 = 10*(time0 -  60*min - 3600*hour);
    sec = sec10/10;
    sec10 = sec10 - 10*sec;
    sprintf(timelabel,"  %i:%.2i:%.2i.%i",hour,min,sec,sec10);
  }
  else{
    float dt;
    char timeval[30], *timevalptr;

    if(ntimes>1){
      dt=times[1]-times[0];
    }
    else{
      dt=0.0;
    }
    if(dt<0.0)dt=-dt;
    timevalptr=time2timelabel(time0,dt,timeval);
    strcpy(timelabel,"Time: ");
    strcat(timelabel,timevalptr);
  }
  sprintf(framelabel,"Frame: %i",itimes);
  if(hrrinfo!=NULL&&hrrinfo->display==1&&hrrinfo->loaded==1){
    float hrr;

    hrr = hrrinfo->hrrval[hrrinfo->itime];
    if(hrr<1.0){
      sprintf(hrrinfo->hrrlabel,"HRR: %4.1f",hrr*1000.0);
    }
    else if(hrr>1000.0){
      sprintf(hrrinfo->hrrlabel,"HRR: %4.1f MW",hrr/1000.0);
    }
    else{
      sprintf(hrrinfo->hrrlabel,"HRR: %4.1f kW",hrr);
    }
  }
}

/* ------------------ ClearBuffers ------------------------ */

void ClearBuffers(int mode){
  if(mode==RENDER){
    glClearColor(backgroundcolor[0],backgroundcolor[1],backgroundcolor[2], 0.0f);
  }
  else{
    glClearColor((float)0.0,(float)0.0,(float)0.0, (float)0.0);
  }
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

/* ------------------ Args ------------------------ */

void Args(int argc, char **argv){
  int i, len;
  char *temp;
  char buffer[255];
  int iarg;
  size_t len_memory;
  char *argi;
  char SMVFILENAME[1024];
  int smv_parse;

  CheckMemory;
  partscale=a_partscale;
  zonescale=a_zonescale;

  if(argc==1){
  //  usage(argv);
    exit(1);
  }
  if(strncmp(argv[1],"-ini",3)==0){
    init_camera_list();
    InitOpenGL();
    updatecolors(-1);
    writeini(GLOBAL_INI);
    exit(0);
  }

  if(strncmp(argv[1],"-ng_ini",6)==0){
    init_camera_list();
    use_graphics=0;
    updatecolors(-1);
    writeini(GLOBAL_INI);
    exit(0);
  }
  strcpy(SMVFILENAME,"");
  smv_parse=0;
  for(iarg=1;iarg<argc;iarg++){
    argi=argv[iarg];
    if(strncmp(argi,"-",1)==0){
      if(
        strncmp(argv[1],"-points",7)==0||
        strncmp(argv[1],"-frames",7)==0||
#ifdef pp_LANG        
        strncmp(argv[1],"-lang",5)==0||
#endif        
        strncmp(argv[1],"-script",7)==0
        ){
        iarg++;
      }

      if(smv_parse==0)continue;
      if(smv_parse==1)break;
    }
    if(smv_parse==1)strcat(SMVFILENAME," ");
    smv_parse=1;
    strcat(SMVFILENAME,argi);
  }

  argi=SMVFILENAME;
#ifndef pp_OSX
  argi=lastname(argi);
#endif
  len = (int) strlen(argi);
  CheckMemory;
  FREEMEMORY(fdsprefix);
  len_memory=len+strlen(part_ext)+100;
  NewMemory((void **)&fdsprefix,(unsigned int)len_memory);
  STRCPY(fdsprefix,argi);
  FREEMEMORY(trainer_filename);
  FREEMEMORY(test_filename);
  FREEMEMORY(filename_sb);

  strcpy(inputfilename_ext,"");

  if(len>4){
    char *c_ext;

    c_ext=strrchr(argi,'.');
    if(c_ext!=NULL){
      STRCPY(inputfilename_ext,c_ext);
      to_lower(inputfilename_ext);

      if(c_ext!=NULL&&(strcmp(inputfilename_ext,".smv")==0||strcmp(inputfilename_ext,".svd")==0||strcmp(inputfilename_ext,".smt")==0)){
        c_ext[0]=0;
        STRCPY(fdsprefix,argi);
        FREEMEMORY(trainer_filename);
        NewMemory((void **)&trainer_filename,(unsigned int)(len+7));
        STRCPY(trainer_filename,argi);
        STRCAT(trainer_filename,".svd");
        FREEMEMORY(test_filename);
        NewMemory((void **)&test_filename,(unsigned int)(len+7));
        STRCPY(test_filename,argi);
        STRCAT(test_filename,".smt");
      }
    }
  }

  FREEMEMORY(logfilename);
  NewMemory((void **)&logfilename,len+4+1);
  STRCPY(logfilename,fdsprefix);
  STRCAT(logfilename,".log");

  FREEMEMORY(caseinifilename);
  NewMemory((void **)&caseinifilename,len+strlen(ini_ext)+1);
  STRCPY(caseinifilename,fdsprefix);
  STRCAT(caseinifilename,ini_ext);

  FREEMEMORY(boundinifilename);
  NewMemory((void **)&boundinifilename,len+5+1);
  STRCPY(boundinifilename,fdsprefix);
  STRCAT(boundinifilename,".bini");

  if(smvfilename==NULL){
    STRUCTSTAT statbuffer;

    NewMemory((void **)&smvfilename,(unsigned int)(len+6));
    FREEMEMORY(smvmenufile);
    NewMemory((void **)&smvmenufile,(unsigned int)(len+15));
    STRCPY(smvfilename,fdsprefix);
    STRCAT(smvfilename,".smv");
    {
      char scriptbuffer[1024];

      STRCPY(scriptbuffer,fdsprefix);
      STRCAT(scriptbuffer,".ssf");
      if(default_script==NULL&&STAT(scriptbuffer,&statbuffer)==0){
        default_script = insert_scriptfile(scriptbuffer);
      }
    }
    STRCPY(smvmenufile,"Reload ");
    temp = strrchr(smvfilename,(int)(*dirseparator));
    if(temp!=NULL){
      STRCAT(smvmenufile,temp+1);
    }
    else{
      STRCAT(smvmenufile,smvfilename);
    }
  }
  if(smvfilename!=NULL){
    STRUCTSTAT statbuffer;

    FREEMEMORY(fds_filein);
    NewMemory((void **)&fds_filein,strlen(fdsprefix)+6);
    STRCPY(fds_filein,fdsprefix);
    STRCAT(fds_filein,".fds");
    if(STAT(fds_filein,&statbuffer)!=0){
      FREEMEMORY(fds_filein);
    }
  }

  // if smokezip created part2iso files then concatenate .smv entries found in the .isosmv file 
  // to the end of the .smv file creating a new .smv file.  Then read in that .smv file.

  {
    FILE *stream_iso=NULL;

    NewMemory((void **)&smvisofilename,len+7+1);
    STRCPY(smvisofilename,fdsprefix);
    STRCAT(smvisofilename,".isosmv");
    stream_iso=fopen(smvisofilename,"r");
    if(stream_iso!=NULL){
      fclose(stream_iso);
    }
    else{
      FREEMEMORY(smvisofilename);
    }
  }

  if(trainer_filename==NULL){
    NewMemory((void **)&trainer_filename,(unsigned int)(len+6));
    STRCPY(trainer_filename,fdsprefix);
    STRCAT(trainer_filename,".svd");
  }
  if(test_filename==NULL){
    NewMemory((void **)&test_filename,(unsigned int)(len+6));
    STRCPY(test_filename,fdsprefix);
    STRCAT(test_filename,".svd");
  }
  if(filename_sb==NULL){
    NewMemory((void **)&filename_sb,(unsigned int)(len+6));
    STRCPY(filename_sb,fdsprefix);
    STRCAT(filename_sb,".sb");
  }

  for (i=1;i<argc;i++){
    if(strncmp(argv[i],"-",1)!=0)continue;
    if(strncmp(argv[i],"-ini",3)==0){
      writeini(GLOBAL_INI);
    }
    else if(strncmp(argv[i],"-update_bounds",14)==0){
      use_graphics=0;
      update_bounds=1;
    }
    else if(strncmp(argv[i],"-demo",5)==0){
       demo_option=1;
    }
    else if(strncmp(argv[i],"-stereo",7)==0){
      if(benchmark==0){
        stereoactive=1;
        showstereo=1;
        printf("stereo option activated\n");
      }
    }
    else if(strncmp(argv[i],"-benchmark",10)==0){
      benchmark=1;
      buffertype=SINGLE_BUFFER;
      printf("benchmark option activated\n");
      if(stereoactive==1){
        stereoactive=0;
        showstereo=0;
        printf("stereo option deactivated\n");
      }
    }
    else if(strncmp(argv[i],"-points",7)==0){
      ++i;
      if(i<argc){
        mxpoints_comm = atol(argv[i]);
        mxpoints=mxpoints_comm;
      }
    }
    else if(strncmp(argv[i],"-frames",7)==0){
      ++i;
      if(i<argc){
        mxframes_comm = atol(argv[i]);
        mxframes=mxframes_comm;
      }
    }
#ifdef pp_LANG
    else if(strncmp(argv[i],"-lang",5)==0){
      ++i;
      if(i<argc){
        int langlen;
        char *lang;

        FREEMEMORY(tr_name);
        lang=argv[i];
        langlen=strlen(lang);
        NewMemory((void **)&tr_name,langlen+48+1);
        strcpy(tr_name,lang);
      }
    }
#endif    
    else if(strncmp(argv[i],"-isotest",8)==0){
      isotest=1;
    }
    else if(strncmp(argv[i],"-h",2)==0){
      usage(argv);
      exit(0);
    }
    else if(strncmp(argv[i],"-noblank",8)==0){
      arg_iblank=1;
      use_iblank=0;
    }
    else if(strncmp(argv[i],"-blank",6)==0){
      arg_iblank=1;
      use_iblank=1;
    }
    else if(
      strncmp(argv[i],"-version",8)==0||
      strncmp(argv[i],"-v",2)==0
      ){
      version();
      exit(0);
    }
    else if(strncmp(argv[i],"-runscript",10)==0){
      runscript=1;
    }
    else if(strncmp(argv[i],"-script",7)==0){
      ++i;
      if(i<argc){
        char scriptbuffer[256];
        scriptfiledata *sfd;

        strcpy(scriptbuffer,argv[i]);
        sfd = insert_scriptfile(scriptbuffer);
        if(sfd!=NULL)default_script=sfd;
        runscript=1;
      }
    }
    else if(strncmp(argv[i],"-noexit",6)==0){
      noexit=1;
    }
    else if(strncmp(argv[i],"-build",6)==0){
      showbuild=1;
      usage(argv);
      exit(0);
    }
    else {
      printf(" unknown option: %s\n",argv[i]);
      usage(argv);
      exit(1);
    }
  }
}

/* ------------------ version ------------------------ */

void version(void){
    int svn_num;

    svn_num=getmaxrevision();    // get svn revision number
    printf("\n");
    printf("%s\n\n",TITLERELEASE);
    printf("Version: %s\n",SMVVERSION);
#ifdef BIT64
    printf("Smokeview (64 bit) Revision Number: %i\n",svn_num);
#else
    printf("Smokeview (32 bit) Revision Number: %i\n",svn_num);
#endif
    if(revision_fds>0){
      printf("FDS Revision Number: %i\n",revision_fds);
    }
    printf("Compile Date: %s\n",__DATE__);
#ifdef WIN32
#ifdef X64
    printf("Platform: WIN64 ");
#else
    printf("Platform: WIN32 ");
#endif
#ifdef pp_WIN_INTEL
    printf(" (Intel C/C++)\n");
#else
    printf(" (MSVS C/C++)\n");
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
    printf("Platform: LINUX\n");
#endif
#endif
#ifdef pp_LINUX64
    printf("Platform: LINUX64\n");
#endif

}

/* ------------------ usage ------------------------ */

void usage(char **argv){
  printf("%s\n",TITLERELEASE);
  printf(_("Visualize fire/smoke flow simulations."));
  printf("\n\n");
  printf(_("Usage: %s [options] casename"),argv[0]);
  printf("\n\n");
  printf(_("where "));
  printf("\n\n");
  printf(_(" casename       - project id (file names without the extension)"));
  printf("\n");
  printf(_(" -build         - show directives used in this build of Smokeview"));
  printf("\n");
  printf(_(" -demo          - use demonstrator mode of Smokeview"));
  printf("\n");
#ifdef pp_DEPRECATED
  printf(_(" -frame nframes - maximum number of particle frames.  Default=%i",MAXFRAMES));
  printf("\n");
#endif
  printf(_(" -help          - display this message"));
  printf("\n");
  printf(_(" -ini           - output default smokeview parameters to smokeview.ini"));
  printf("\n");
  printf(_(" -ng_ini        - No graphics version of -ini."));
  printf("\n");
#ifdef pp_DEPRECATED
  printf(_(" -points npoints - maximum number of particles.  Default=%i"),MAXPOINTS);
  printf("\n");
#endif
  printf(_(" -runscript     - run the script file casename.ssf"));
  printf("\n");
  printf(_(" -script scriptfile - run the script file scriptfile"));
  printf("\n");
  printf(_(" -stereo        - activate stereo mode"));
  printf("\n");
  printf(_(" -update_bounds - calculate boundary file bounds and save to casename.bini"));
  printf("\n");
  printf(_(" -version       - display version information"));
  printf("\n");

  if(showbuild==1){
    char label[1024],*labelptr;

    labelptr=label+2;
    strcpy(label,"");
#ifdef pp_BENCHMARK
    strcat(label,", pp_BENCHMARK");
#endif
#ifdef pp_COMPRESS
    strcat(label,", pp_COMPRESS");
#endif
#ifdef pp_CULL
    strcat(label,", pp_CULL");
#endif
#ifdef _DEBUG
    strcat(label," _DEBUG");
#endif
#ifdef pp_DRAWISO
    strcat(label,", pp_DRAWISO");
#endif
#ifdef EGZ
    strcat(label,", EGZ");
#endif
#ifdef pp_GPU
    strcat(label,", pp_GPU");
#endif
#ifdef pp_GPUDEPTH
    strcat(label,", pp_GPUDEPTH");
#endif
#ifdef ISO_DEBUG
    strcat(label,", ISO_DEBUG");
#endif
#ifdef pp_JPEG
    strcat(label,", pp_JPEG");
#endif
#ifdef pp_LIGHT
    strcat(label,", pp_LIGHT");
#endif
#ifdef pp_LINUX64
    strcat(label,", pp_LINUX64");
#endif
#ifdef pp_memstatus
    strcat(label,", pp_memstatus");
#endif
#ifdef pp_MESSAGE
    strcat(label,", pp_MESSAGE");
#endif
#ifdef pp_MOUSEDOWN
    strcat(label,", pp_MOUSEDOWN");
#endif
#ifdef pp_OPEN
    strcat(label,", pp_OPEN");
#endif
#ifdef pp_OSX64
    strcat(label,", pp_OSX64");
#endif
#ifdef pp_noappend
    strcat(label,", pp_noappend");
#endif
#ifdef pp_OSX
    strcat(label,", pp_OSX");
#endif
#ifdef pp_release
    strcat(label," pp_release");
#endif
#ifdef pp_SHOOTER
    strcat(label,", pp_SHOOTER");
#endif
#ifdef pp_SHOWLIGHT
    strcat(label,", pp_SHOWLIGHT");
#endif
#ifdef pp_SLICECONTOURS
    strcat(label,", pp_SLICECONTOURS");
#endif
#ifdef pp_SMOKETEST
    strcat(label,", pp_SMOKETEST");
#endif
#ifdef pp_THREAD
    strcat(label,", pp_THREAD");
#endif
#ifdef X64
    strcat(label,", X64");
#endif
#ifdef USE_ZLIB
    strcat(label,", USE_ZLIB");
#endif
#ifdef WIN32
    strcat(label,", WIN32");
#endif
    printf("  \n");
    printf(_("  Smokeview was built with the following pre-processing directives set:"));
    printf("\n\n");
    printf("%s \n",labelptr);
  }
}

/* ------------------ checktimebound ------------------------ */

void checktimebound(void){
  int i,j;
  slice *sd;
  mesh *meshi;
  blockagedata *bc;
  particle *parti;

  if(timedrag==0&&itimes>ntimes-1||timedrag==1&&itimes<0){
    izone=0;
    itimes=0;
    iframe=iframebeg;
    for(i=0;i<nsliceinfo;i++){
      sd=sliceinfo+i;
      sd->islice=0;
    }
    for(i=0;i<nmeshes;i++){
      meshi=meshinfo+i;
      meshi->ipatch=0;
    }
    for(i=0;i<nmeshes;i++){
      meshi=meshinfo+i;
      if(meshi->isotimes==NULL)continue;
      meshi->iiso=0;
    }
  }
  if(timedrag==0&&itimes<0||timedrag==1&&itimes>ntimes-1){
    izone=nzonet-1;
    itimes=ntimes-1;
    for(i=0;i<npartinfo;i++){
      parti=partinfo+i;
      parti->iframe=parti->nframes-1;
    }
    for(i=0;i<nsliceinfo;i++){
      sd=sliceinfo+i;
      sd->islice=sd->nsteps-1;
    }
    for(i=0;i<nmeshes;i++){
      meshi=meshinfo+i;
      meshi->ipatch=meshi->npatch_frames-1;
    }
    for(i=0;i<nmeshes;i++){
      meshi=meshinfo+i;
      if(meshi->isotimes==NULL)continue;
      meshi->iiso=meshi->nisosteps-1;
    }
  }
  /* set blockage visibility */

  for(i=0;i<nmeshes;i++){
    meshi=meshinfo+i;
    for(j=0;j<meshi->nbptrs;j++){
      bc=meshi->blockageinfoptrs[j];
      if(bc->showtimelist==NULL)continue;
      bc->show=bc->showtimelist[itimes];
    }
  }
}

/* ------------------ filecat ------------------------ */

int filecat(char *file_in1, char *file_in2, char *file_out){
#define SIZEBUFFER 10000
  char buffer[SIZEBUFFER];
  FILE *stream_in1, *stream_in2, *stream_out;
  int chars_in;

  if(file_in1==NULL||file_in2==NULL)return -1;
  if(file_out==NULL)return -2;

  stream_in1=fopen(file_in1,"r");
  if(stream_in1==NULL)return -1;

  stream_in2=fopen(file_in2,"r");
  if(stream_in2==NULL){
    fclose(stream_in1);
    return -1;
  }

  stream_out=fopen(file_out,"w");
  if(stream_out==NULL){
    fclose(stream_in1);
    fclose(stream_in2);
    return -2;
  }

  for(;;){
    int eof;
       
    eof=0;
    chars_in=fread(buffer,1,SIZEBUFFER,stream_in1);
    if(chars_in!=SIZEBUFFER)eof=1;
    if(chars_in>0)fwrite(buffer,chars_in,1,stream_out);
    if(eof==1)break;
  }
  fclose(stream_in1);

  for(;;){
    int eof;
       
    eof=0;
    chars_in=fread(buffer,1,SIZEBUFFER,stream_in2);
    if(chars_in!=SIZEBUFFER)eof=1;
    if(chars_in>0)fwrite(buffer,chars_in,1,stream_out);
    if(eof==1)break;
  }
  fclose(stream_in2);
  fclose(stream_out);
  return 0;

}
