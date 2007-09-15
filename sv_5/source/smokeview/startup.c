#include "options.h"
#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "flowfiles.h"
#include "MALLOC.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"


//void get_smokezippath(char *progdir, char **zippath);
void glui_colorbar_setup(int main_window);
void glui_motion_setup(int main_window);
void glui_bounds_setup(int main_window);
void glui_labels_setup(int main_window);
void glui_edit_setup(int main_window);
void glui_clip_setup(int main_window);
void glui_tour_setup(int main_window);
void glui_advancedtour_setup(int main_window);
void glui_stereo_setup(int main_window);
void glui_trainer_setup(int main_window);
void glui_3dsmoke_setup2(int main_window);
void RenderMenu(int value);
//void InitOpenGL(void);

/* ------------------ initcase_c ------------------------ */

void to_lower(char *string){
   size_t lenstring;
   char c;
   size_t i;

   if(string==NULL)return;
   lenstring=strlen(string);
   for(i=0;i<lenstring;i++){
     c=string[i];
     if(c>='A'&&c<='Z')string[i]=c+'a'-'A';
   }
}

/* ------------------ initcase_c ------------------------ */

int initcase_c(int argc, char **argv){
  int return_code;
  int input_type=0;

  /* 
  warning: the following line was commented out!! (perhaps it broke something)
     this line is necessary in order to define smvfilename and trainer_filename
  */
  Args(argc, argv); 
  return_code=-1;
  if(strcmp(inputfilename_ext,".svt")==0){
    input_type=1;
    return_code=readsmv(trainer_filename);
    if(return_code==0){
      trainer_mode=1;
      trainer_active=1;
      if(trainer_mode==1){
        show_trainer();
        show_load_alert();
      }
    }
  }
  else{
    input_type=2;
    return_code=readsmv(smvfilename);
  }
  switch (return_code){
  case 1:
    if(input_type==1){
      printf("Training input files %s not found.\n",trainer_filename);
    }
    else{
      printf("Smokeview input files %s not found.\n",smvfilename);
    }
    pauseSV();
    return 1;
  case 2:
      printf("*** Fatal error: unable to allocate necessary memory\n");
      pauseSV();
      return 2;
  case 0:
    break;
  default:
    ASSERT(FFALSE);
  }
  readini(0);
  if(ntours==0)setup_tour();
  glui_colorbar_setup(mainwindow_id);
  glui_motion_setup(mainwindow_id);
  glui_bounds_setup(mainwindow_id);
  glui_edit_setup(mainwindow_id);
  glui_clip_setup(mainwindow_id);
  glui_labels_setup(mainwindow_id);
  glui_tour_setup(mainwindow_id);
  glui_alert_setup(mainwindow_id);
  glui_advancedtour_setup(mainwindow_id);
  glui_stereo_setup(mainwindow_id);
  glui_3dsmoke_setup2(mainwindow_id);


  if(UpdateLIGHTS==1)updateLights(0);

  glutReshapeWindow(screenWidth,screenHeight);

  glutSetWindow(mainwindow_id);
  glutShowWindow();
  glutSetWindowTitle(fdsprefix);
  Init();
  glui_trainer_setup(mainwindow_id);
  glutDetachMenu(GLUT_RIGHT_BUTTON);
  InitMenus(LOAD);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
  if(trainer_mode==1){
    show_trainer();
    show_load_alert();
  }
  return 0;
}


/* ------------------ sv_startup_c ------------------------ */

void sv_startup_c(int argc, char **argv){
  int i;
  f_unit *units;
  f_units *ut;
  char *smoketempdir;
  size_t lensmoketempdir,lensmokebindir;
#ifdef pp_OSX
  char workingdir[1000];
#endif

// get smokeview bin directory from argv[0] which contains the full path of the smokeview binary

  lensmokebindir = strlen(argv[0]);
  NewMemory((void **)&smokeviewbindir, (unsigned int)(lensmokebindir+1));

  strcpy(smokeviewbindir,argv[0]);
  for(i=lensmokebindir-1;i>=0;i--){
    char c;

    c = smokeviewbindir[i];
    if(strncmp(&c,dirseparator,1)==0){
      smokeviewbindir[i+1]=0;
      break;
    }
  }

  NewMemory((void **)&smokeviewini,    (unsigned int)(lensmokebindir+14));
  STRCPY(smokeviewini,smokeviewbindir);
  STRCAT(smokeviewini,"smokeview.ini");
  
  startup_pass=2;

  smoketempdir=getenv("SVTEMPDIR");
  if(smoketempdir==NULL)smoketempdir=getenv("svtempdir");
  if(smoketempdir==NULL)smoketempdir=getenv("TEMP");
  if(smoketempdir==NULL)smoketempdir=getenv("temp");

  if(smoketempdir != NULL){
    lensmoketempdir = strlen(smoketempdir);
    if(NewMemory((void **)&smokeviewtempdir,(unsigned int)(lensmoketempdir+2))!=0){
      STRCPY(smokeviewtempdir,smoketempdir);
      if(strncmp(smokeviewtempdir+lensmoketempdir-1,dirseparator,1)!=0){
        STRCAT(smokeviewtempdir,dirseparator);
      }
      printf("Scratch directory: %s\n",smokeviewtempdir);
    }
  }

  if(texturedir==NULL){
    char *texture_buffer;
    size_t texture_len;

    texture_buffer=getenv("texturedir");
    if(texture_buffer!=NULL){
      texture_len=strlen(texture_buffer);
      NewMemory((void **)&texturedir,texture_len+1);
      strcpy(texturedir,texture_buffer);
    }
  }

  printf("*********************************************************************\n");
  printf("*********************************************************************\n");
#ifdef pp_NISTREVIEW
  printf("This version of Smokeview is intended for review and testing ONLY.\n");
  printf("\n");
#endif
  printf("The US Department of Commerce makes no warranty, expressed or\n");
  printf("implied, to users of Smokeview, and accepts no responsibility\n");
  printf("for its use. Users of Smokeview assume sole responsibility under\n");
  printf("Federal law for determining the appropriateness of its use in any\n");
  printf("particular application; for any conclusions drawn from the results\n"); 
  printf("of its use; and for any actions taken or not taken as a result of\n"); 
  printf("analysis performed using this tools.\n");
  printf("\n");
  printf("Smokeview and the companion program FDS is intended for use only\n");
  printf("by those competent in the fields of fluid dynamics, thermodynamics,\n");
  printf("combustion, and heat transfer, and is intended only to supplement\n");
  printf("the informed judgment of the qualified user. These software packages\n");
  printf("may or may not have predictive capability when applied to a specific\n");
  printf("set of factual circumstances.  Lack of accurate predictions could lead\n");
  printf("to erroneous conclusions with regard to fire safety.  All results\n");
  printf("should be evaluated by an informed user.\n\n");
  printf("*********************************************************************\n");
  printf("*********************************************************************\n");

#ifdef pp_OSX
  getcwd(workingdir,1000);
#endif
  glutInit(&argc, argv);
#ifdef pp_OSX
  chdir(workingdir);
#endif

  mxpoints=mxpoints_orig;
  mxframes=mxframes_orig;
  glutInitWindowSize(screenWidth, screenHeight);
  max_screenWidth = glutGet(GLUT_SCREEN_WIDTH);
  max_screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
  if(trainer_mode==1){
    int TRAINER_WIDTH;
    TRAINER_WIDTH=300;
    screenWidth = glutGet(GLUT_SCREEN_WIDTH)-TRAINER_WIDTH;
    screenHeight = glutGet(GLUT_SCREEN_HEIGHT)-50; 
    max_screenWidth = screenWidth;
    max_screenHeight = screenHeight;
  }
  InitOpenGL();

#ifdef pp_SVNET
  init_svcom();
#endif
  /* initialize units */

#ifdef COMMENT
typedef struct {
  char unit[10];   /* m/s, mph etc - appears in color bar */
  float scale[2];  /* newval=scale[0]*oldval+scale[1] */
} f_unit;

typedef struct {
  int nunits;
  char unitclass[30]; /* ie: velocity, temperature */
  char unittype[80]; /* list of equivalent units separated by ';' ie: U-VEL;V-VEL;W-VEL */
  f_unit *units;
} f_units;
#endif


  nunitclasses_default=2;
  NewMemory((void **)&unitclasses_default,nunitclasses_default*sizeof(f_units));

  
  for(i=0;i<nunitclasses_default;i++){
    unitclasses_default[i].submenuid=0;
  }
  ut=unitclasses_default;

  ut->nunits=3;
  ut->active=0;
  strcpy(ut->unitclass,"temperature");
  ut->nunittypes=2;
  strcpy(ut->unittypes[0],"temp");
  strcpy(ut->unittypes[1],"d_temp");

  NewMemory((void **)&(ut->units),ut->nunits*sizeof(f_unit));
  units=ut->units;
  strcpy(units[0].unit,"C");
  units[0].scale[0]=1.0;
  units[0].scale[1]=0.0;
  strcpy(units[1].unit,"F");
  units[1].scale[0]=1.8;
  units[1].scale[1]=32.0;
  strcpy(units[2].unit,"K");
  units[2].scale[0]=1.0;
  units[2].scale[1]=273.15;


  ut = unitclasses_default+1;
  ut->active=0;
  ut->nunits=3;
  strcpy(ut->unitclass,"velocity");
  ut->nunittypes=4;
  strcpy(ut->unittypes[0],"U-VEL");
  strcpy(ut->unittypes[1],"V-VEL");
  strcpy(ut->unittypes[2],"W-VEL");
  strcpy(ut->unittypes[3],"Speed");

  NewMemory((void **)&(ut->units),ut->nunits*sizeof(f_unit));
  units=ut->units;
  strcpy(units[0].unit,"m/s");
  units[0].scale[0]=1.0;
  units[0].scale[1]=0.0;
  strcpy(units[1].unit,"mph");
  units[1].scale[0]=2.236931818;
  units[1].scale[1]=0.0;
  strcpy(units[2].unit,"f/s");
  units[2].scale[0]=3.2808333333;
  units[2].scale[1]=0.0;

  NewMemory((void **)&rgbptr,MAXRGB*sizeof(float *));
  for(i=0;i<MAXRGB;i++){
    rgbptr[i]=&rgb[i][0];
  }
}

/* ------------------ InitOpenGL ------------------------ */

void InitOpenGL(void){
  int type;
  
  type = GLUT_RGB|GLUT_DEPTH;
  if(buffertype==GLUT_DOUBLE){
    type |= GLUT_DOUBLE;
  }
  else{
    type |= GLUT_SINGLE;
  }

//  glutInitDisplayMode(GLUT_STEREO);
  if(stereoactive==1){
    if(glutGet(GLUT_DISPLAY_MODE_POSSIBLE)==1){
      videoSTEREO=1;
      type |= GLUT_STEREO;
    }
    else{
      videoSTEREO=0;
      printf("***warning: video hardware does not support stereo\n");
    }
  }
  if(videoSTEREO==1){
    stereo_frame=1;
    stereo_leftright=0;
    stereo_off=0;
  }
  else{
    stereo_frame=0;
    stereo_leftright=0;
    stereo_off=1;
  }

  glutInitDisplayMode(type);

  mainwindow_id=glutCreateWindow("");

  glutSpecialUpFunc(specialkeyboard_up);
  glutKeyboardUpFunc(keyboard_up);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutSpecialFunc(specialkeyboard);
  glutMotionFunc(motion);
  glutReshapeFunc(Reshape);
  glutDisplayFunc(Display);
  glutVisibilityFunc(NULL);
  glutMenuStatusFunc(MenuStatus);
//  glutWindowStatusFunc(WindowStatus);

#ifdef pp_GPU
  init_shaders();
#endif

  light_position0[0]=1.0f;
  light_position0[1]=1.0f;
  light_position0[2]=4.0f; 
  light_position0[3]=0.f;

  light_position1[0]=-1.0f;
  light_position1[1]=1.0f;
  light_position1[2]=4.0f;
  light_position1[3]=0.f;

  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,GL_TRUE);
  updateLights(1);

  {
    glGetIntegerv(GL_RED_BITS,&nredbits);    
    glGetIntegerv(GL_GREEN_BITS,&ngreenbits);
    glGetIntegerv(GL_BLUE_BITS,&nbluebits);

    nredshift = 8 - nredbits;
    if(nredshift<0)nredshift=0;
    ngreenshift = 8 - ngreenbits;
    if(ngreenshift<0)ngreenshift=0;
    nblueshift=8-nbluebits;
    if(nblueshift<0)nblueshift=0;
  }
  opengldefined=1;
}

/* ------------------ set_3dsmoke_startup ------------------------ */

 void set_3dsmoke_startup(void){
   int i;

    for(i=0;i<nvslice;i++){
      vslice *vslicei;

      vslicei = vsliceinfo + i;

      if(vslicei->loaded==1){
        vslicei->autoload=1;
      }
      else{
        vslicei->autoload=0;
      }
    }
    for(i=0;i<npartinfo;i++){
      particle *parti;

      parti = partinfo + i;

      if(parti->loaded==1){
        parti->autoload=1;
      }
      else{
        parti->autoload=0;
      }
    }
    for(i=0;i<nplot3d;i++){
      plot3d *plot3di;

      plot3di = plot3dinfo + i;

      if(plot3di->loaded==1){
        plot3di->autoload=1;
      }
      else{
        plot3di->autoload=0;
      }
    }
    for(i=0;i<nsmoke3d;i++){
      smoke3d *smoke3di;

      smoke3di = smoke3dinfo + i;

      if(smoke3di->loaded==1){
        smoke3di->autoload=1;
      }
      else{
        smoke3di->autoload=0;
      }
    }
    for(i=0;i<npatch_files;i++){
      patch *patchi;

      patchi = patchinfo + i;

      if(patchi->loaded==1){
        patchi->autoload=1;
      }
      else{
        patchi->autoload=0;
      }
    }
    for(i=0;i<niso;i++){
      iso *isoi;

      isoi = isoinfo + i;

      if(isoi->loaded==1){
        isoi->autoload=1;
      }
      else{
        isoi->autoload=0;
      }
    }
    for(i=0;i<nslice;i++){
      slice *slicei;

      slicei = sliceinfo + i;

      if(slicei->loaded==1){
        slicei->autoload=1;
      }
      else{
        slicei->autoload=0;
      }
    }
 }

 /* ------------------ clear_3dsmoke_startup ------------------------ */

 void clear_3dsmoke_startup(void){
   int i;

    for(i=0;i<nvslice;i++){
      vslice *vslicei;

      vslicei = vsliceinfo + i;

      vslicei->autoload=0;
    }

    for(i=0;i<npartinfo;i++){
      particle *parti;

      parti = partinfo + i;

      parti->autoload=0;
    }

    for(i=0;i<nplot3d;i++){
      plot3d *plot3di;

      plot3di = plot3dinfo + i;

      plot3di->autoload=0;
    }

    for(i=0;i<nsmoke3d;i++){
      smoke3d *smoke3di;

      smoke3di = smoke3dinfo + i;

      smoke3di->autoload=0;
    }

    for(i=0;i<npatch_files;i++){
      patch *patchi;

      patchi = patchinfo + i;

      patchi->autoload=0;
    }

    for(i=0;i<niso;i++){
      iso *isoi;

      isoi = isoinfo + i;

      isoi->autoload=0;
    }
    for(i=0;i<nslice;i++){
      slice *slicei;

      slicei = sliceinfo + i;

      slicei->autoload=0;
    }
 }

 /* ------------------ put_startup_smoke3d ------------------------ */

  void put_startup_smoke3d(FILE *fileout){
   int i;
   int nstartup;

   if(fileout==NULL)return;


   // startup particle

   nstartup=0;
   for(i=0;i<npartinfo;i++){
      particle *parti;

      parti = partinfo + i;

      if(parti->loaded==1)nstartup++;
   }
   if(nstartup!=0){
     fprintf(fileout,"PARTAUTO\n");
     fprintf(fileout," %i \n",nstartup);
     for(i=0;i<npartinfo;i++){
        particle *parti;

        parti = partinfo + i;

        if(parti->loaded==1)fprintf(fileout," %i\n",parti->seq_id);
     }
   }

   // startup plot3d

   nstartup=0;
   for(i=0;i<nplot3d;i++){
      plot3d *plot3di;

      plot3di = plot3dinfo + i;

      if(plot3di->loaded==1)nstartup++;
   }
   if(nstartup!=0){
     fprintf(fileout,"PLOT3DAUTO\n");
     fprintf(fileout," %i \n",nstartup);
     for(i=0;i<nplot3d;i++){
        plot3d *plot3di;

        plot3di = plot3dinfo + i;

        if(plot3di->loaded==1)fprintf(fileout," %i\n",plot3di->seq_id);
     }
   }

   // startup iso

   nstartup=0;
   for(i=0;i<niso;i++){
      iso *isoi;

      isoi = isoinfo + i;

      if(isoi->loaded==1)nstartup++;
   }
   if(nstartup!=0){
     fprintf(fileout,"ISOAUTO\n");
     fprintf(fileout," %i \n",nstartup);
     for(i=0;i<niso;i++){
        iso *isoi;

        isoi = isoinfo + i;

        if(isoi->loaded==1)fprintf(fileout," %i\n",isoi->seq_id);
     }
   }

   // startup vslice

   nstartup=0;
   for(i=0;i<nvslice;i++){
      vslice *vslicei;

      vslicei = vsliceinfo + i;

      if(vslicei->loaded==1)nstartup++;
   }
   if(nstartup!=0){
     fprintf(fileout,"VSLICEAUTO\n");
     fprintf(fileout," %i \n",nstartup);
     for(i=0;i<nvslice;i++){
        vslice *vslicei;

        vslicei = vsliceinfo + i;

        if(vslicei->loaded==1)fprintf(fileout," %i\n",vslicei->seq_id);
     }
   }

   // startup slice

   nstartup=0;
   for(i=0;i<nslice;i++){
      slice *slicei;

      slicei = sliceinfo + i;

      if(slicei->loaded==1)nstartup++;
   }
   if(nstartup!=0){
     fprintf(fileout,"SLICEAUTO\n");
     fprintf(fileout," %i \n",nstartup);
     for(i=0;i<nslice;i++){
        slice *slicei;

        slicei = sliceinfo + i;
        if(slicei->loaded==1)fprintf(fileout," %i\n",slicei->seq_id);
     }
   }

   // startup smoke

   nstartup=0;
   for(i=0;i<nsmoke3d;i++){
      smoke3d *smoke3di;

      smoke3di = smoke3dinfo + i;

      if(smoke3di->loaded==1)nstartup++;
   }
   if(nstartup!=0){
     fprintf(fileout,"S3DAUTO\n");
     fprintf(fileout," %i \n",nstartup);
     for(i=0;i<nsmoke3d;i++){
        smoke3d *smoke3di;

        smoke3di = smoke3dinfo + i;

        if(smoke3di->loaded==1)fprintf(fileout," %i\n",smoke3di->seq_id);
     }
   }

   // startup patch

   nstartup=0;
   for(i=0;i<npatch_files;i++){
      patch *patchi;

      patchi = patchinfo + i;

      if(patchi->loaded==1)nstartup++;
   }
   if(nstartup!=0){
     fprintf(fileout,"PATCHAUTO\n");
     fprintf(fileout," %i \n",nstartup);
     for(i=0;i<npatch_files;i++){
        patch *patchi;

        patchi = patchinfo + i;

        if(patchi->loaded==1)fprintf(fileout," %i\n",patchi->seq_id);
     }
   }

 }

 /* ------------------ get_startup_part ------------------------ */

  void get_startup_part(int seq_id){
    int i;
    for(i=0;i<npartinfo;i++){
      particle *parti;

      parti = partinfo + i;
      if(parti->seq_id==seq_id){
        parti->autoload=1;
        return;
      }
    }
  }

 /* ------------------ get_startup_plot3d ------------------------ */

  void get_startup_plot3d(int seq_id){
    int i;
    for(i=0;i<nplot3d;i++){
      plot3d *plot3di;

      plot3di = plot3dinfo + i;
      if(plot3di->seq_id==seq_id){
        plot3di->autoload=1;
        return;
      }
    }
  }

 /* ------------------ get_startup_patch ------------------------ */

  void get_startup_patch(int seq_id){
    int i;
    for(i=0;i<npatch_files;i++){
      patch *patchi;

      patchi = patchinfo + i;
      if(patchi->seq_id==seq_id){
        patchi->autoload=1;
        return;
      }
    }
  }

 /* ------------------ get_startup_smoke3d ------------------------ */

  void get_startup_smoke(int seq_id){
    int i;
    for(i=0;i<nsmoke3d;i++){
      smoke3d *smoke3di;

      smoke3di = smoke3dinfo + i;

      if(smoke3di->seq_id==seq_id){
        smoke3di->autoload=1;
        return;
      }
    }
  }


 /* ------------------ get_startup_iso ------------------------ */

  void get_startup_iso(int seq_id){
    int i;
    for(i=0;i<niso;i++){
      iso *isoi;

      isoi = isoinfo + i;

      if(isoi->seq_id==seq_id){
        isoi->autoload=1;
        return;
      }
    }
  }

 /* ------------------ get_startup_slice ------------------------ */

  void get_startup_slice(int seq_id){
    int i;
    for(i=0;i<nslice;i++){
      slice *slicei;

      slicei = sliceinfo + i;

      if(slicei->seq_id==seq_id){
        slicei->autoload=1;
        return;
      }
    }
  }


 /* ------------------ get_startup_vslice ------------------------ */

  void get_startup_vslice(int seq_id){
    int i;
    for(i=0;i<nvslice;i++){
      vslice *vslicei;

      vslicei = vsliceinfo + i;

      if(vslicei->seq_id==seq_id){
        vslicei->autoload=1;
        return;
      }
    }
  }

 /* ------------------ load_startup_smoke3d ------------------------ */

  void load_startup_smoke(void){
    int i;
    int errorcode;

//    show_load_alert();
    for(i=0;i<nplot3d;i++){
      plot3d *plot3di;

      plot3di = plot3dinfo + i;
      if(plot3di->autoload==0&&plot3di->loaded==1){
        readplot(plot3di->file,i,UNLOAD,&errorcode);
      }
      if(plot3di->autoload==1){
        ReadPlot3dFile=1;
        readplot(plot3di->file,i,LOAD,&errorcode);
      }
    }
    for(i=0;i<npartinfo;i++){
      particle *parti;

      parti = partinfo + i;
      if(parti->autoload==0&&parti->loaded==1)readpart(parti->file,i,UNLOAD,&errorcode);
      if(parti->autoload==1)readpart(parti->file,i,LOAD,&errorcode);
    }
    for(i=0;i<niso;i++){
      iso *isoi;

      isoi = isoinfo + i;
      if(isoi->autoload==0&&isoi->autoload==1)readiso(isoi->file,i,UNLOAD,&errorcode);
      if(isoi->autoload==1)readiso(isoi->file,i,LOAD,&errorcode);
    }
    for(i=0;i<nvslice;i++){
      vslice *vslicei;

      vslicei = vsliceinfo + i;
      if(vslicei->autoload==0&&vslicei->loaded==1)readvslice(i,UNLOAD,&errorcode);
      if(vslicei->autoload==1){
        readvslice(i,LOAD,&errorcode);
      }
    }
    // note:  only slices that are NOT a part of a vector slice will be loaded here
    for(i=0;i<nslice;i++){
      slice *slicei;

      slicei = sliceinfo + i;
      if(slicei->autoload==0&&slicei->loaded==1)readslice(slicei->file,i,UNLOAD,&errorcode);
      if(slicei->autoload==1&&slicei->loaded==0){
        readslice(slicei->file,i,LOAD,&errorcode);
      }
    }
#ifdef pp_WUI
    for(i=0;i<nterraininfo;i++){
      terraindata *terri;

      terri = terraininfo + i;
      if(terri->autoload==0&&terri->loaded==1)readterrain(terri->file,i,UNLOAD,&errorcode);
      if(terri->autoload==1&&terri->loaded==0)readslice(terri->file,i,LOAD,&errorcode);
    }
#endif
    for(i=0;i<nsmoke3d;i++){
      smoke3d *smoke3di;

      smoke3di = smoke3dinfo + i;
      if(smoke3di->autoload==0&&smoke3di->loaded==1)readsmoke3d(i,UNLOAD,&errorcode);
      if(smoke3di->autoload==1)readsmoke3d(i,LOAD,&errorcode);
    }
    for(i=0;i<npatch_files;i++){
      patch *patchi;

      patchi = patchinfo + i;
      if(patchi->autoload==0&&patchi->loaded==1)readpatch(i,UNLOAD,&errorcode);
      if(patchi->autoload==1)readpatch(i,LOAD,&errorcode);
    }
    force_redisplay=1;
    update_framenumber(0);
    updatemenu=1;
    update_load_startup=0;
    hide_load_alert();
    TrainerViewMenu(trainerview);
  }


/* ------------------ initvars1 ------------------------ */

void initvars1(void){
  force_isometric=0;
#ifdef pp_SPOTLIGHT
  spotlight=0;
#endif
  cb_valmin=0.0;
  cb_valmax=100.0;
  cb_val=50.0;
  cb_colorindex=128;

#ifdef pp_WUI
  rgb_terrain[0][0]=1.0;
  rgb_terrain[0][1]=0.0;
  rgb_terrain[0][2]=0.0;
  rgb_terrain[0][3]=1.0;

  rgb_terrain[1][0]=0.5;
  rgb_terrain[1][1]=0.5;
  rgb_terrain[1][2]=0.0;
  rgb_terrain[1][3]=1.0;

  rgb_terrain[2][0]=0.0;
  rgb_terrain[2][1]=1.0;
  rgb_terrain[2][2]=0.0;
  rgb_terrain[2][3]=1.0;

  rgb_terrain[3][0]=0.0;
  rgb_terrain[3][1]=0.5;
  rgb_terrain[3][2]=0.0;
  rgb_terrain[3][3]=1.0;

  rgb_terrain[4][0]=0.0;
  rgb_terrain[4][1]=0.5;
  rgb_terrain[4][2]=0.5;
  rgb_terrain[4][3]=1.0;

  rgb_terrain[5][0]=0.0;
  rgb_terrain[5][1]=0.0;
  rgb_terrain[5][2]=1.0;
  rgb_terrain[5][3]=1.0;

  rgb_terrain[6][0]=0.5;
  rgb_terrain[6][1]=0.0;
  rgb_terrain[6][2]=0.5;
  rgb_terrain[6][3]=1.0;

  rgb_terrain[7][0]=1.0;
  rgb_terrain[7][1]=0.5;
  rgb_terrain[7][2]=0.0;
  rgb_terrain[7][3]=1.0;

  rgb_terrain[8][0]=1.0;
  rgb_terrain[8][1]=0.5;
  rgb_terrain[8][2]=0.5;
  rgb_terrain[8][3]=1.0;

  rgb_terrain[9][0]=1.0;
  rgb_terrain[9][1]=0.25;
  rgb_terrain[9][2]=0.5;
  rgb_terrain[9][3]=1.0;
#endif

  render_double=0;
  render_double_state=0;
  usetexturebar=1;
  percentile_level=0.01;

  trainerview=1;
  show_bothsides_int=1;
  show_bothsides_ext=0;
  show_slice_in_obst=0;
  updategluiview=1;
  trainer_pause=0;
  trainee_location=0;
  trainer_inside=0;
  from_glui_trainer=0;
  trainer_path_old=-3;
  trainer_outline=1;
  trainer_viewpoints=-1;
  trainer_realtime=1;
  trainer_path=0;
  trainer_xzy[0]=0.0;
  trainer_xzy[1]=0.0;
  trainer_xzy[2]=0.0;
  trainer_ab[0]=0.0;
  trainer_ab[1]=0.0;
  motion_ab[0]=0.0;
  motion_ab[1]=0.0;
  motion_dir[0]=0.0;
  motion_dir[1]=0.0;
  trainerload=0;
  fontsize_save=0;
  showtrainer=0;
//#ifdef pp_TRAINER
//  trainer_mode=1;
//  trainer_active=1;
//#else
  trainer_mode=0;
  trainer_active=0;
//#endif
  slice_average_flag=0;
  show_slice_average=0;
  vis_slice_average=1;
  slice_average_interval=10.0;

  angle=0.0, dang=0.025f, tourangle=0.0;
  maxtourframes=500;
  blockageSelect=0;
  ntourknots=0;
  itourknots=-1;
  stretch_var_black=0; 
  stretch_var_white=0; 
  move_var=0;

  showhide_option=0;
  snifferrornumber=0;
  xyz_dir=0;
  xyz_blockage_dir=0;
  which_face=2;
  showfontmenu=1;
  showlightmenu=0;

  VECFRACTION=1.0;
  vecfactor=0.25f;
  veclength=0.0f;
  iveclengths=0;

  glui_active=0;

  drawColorLabel=0;
  olddrawColorLabel=0;
  staticframe0=0;
  visStaticSmoke=1;
  vis3DSmoke3D=1;
  smokeskip=1;
  smokeskipm1=0;
  nrooms=0;
  nzone=0;
  nzvents=0;
  nfires=0;
  ratio=1.0;
  aspect=1.;
  visLIGHT0=1;
  visLIGHT1=1;
  visLIGHTMENU=1;
  UpdateLIGHTS=1;

  screenWidth = 640, screenHeight = 480;
  renderW = 640, renderH=480;
  glui_screenWidth=640, glui_screenHeight=480;
  windowsize_pointer=0;
  sethazardcolor=0;
  mxpoints=MAXPOINTS,mxframes=MAXFRAMES,mxframepoints;
  mxpoints_orig=MAXPOINTS,mxframes_orig=MAXFRAMES;
  mxpoints_comm=0, mxframes_comm=0;
  timedrag=0,colordrag=0;
  isonormtype=1,showisonormals=0;
  global_changecolorindex=-1;
  fontindex=0,fontWoffset=0,fontHoffset=0;

  xcenGLOBAL=0.5, ycenGLOBAL=0.5, zcenGLOBAL=0.5;
  xbar=1.0, ybar=1.0, zbar=1.0;
  xbar0=0.0, ybar0=0.0, zbar0=0.0;
  xbarORIG=1.0, ybarORIG=1.0, zbarORIG=1.0;
  xbar0ORIG=0.0, ybar0ORIG=0.0, zbar0ORIG=0.0;
  ReadPlot3dFile=0, ReadIsoFile=0;

  ReadVolSlice=0;
  Read3DSmoke3DFile=0;
  ReadZoneFile=0, ReadPartFile=0, ReadEvacFile=0;;

  editwindow_status=-1;
  startup_pass=1;
  ntargtimes=500;
  showtitle1=0, showtitle2=0;

  slicefilenumber=0;
  exportdata=0;
  count=1, lastcount=1;
  nspr=0;
  RenderGif=0, RenderSkip=1;
  isoframestep=1;
  isoframeskip=0;
  smoke3dframestep=1;
  smoke3dframeskip=0;
  vectorskip=1;
  iframe=0, iframebeg=0, izone=0;
  eyeview=0,eyeview_level=1;
  eyeview_old=0,eyeview_SAVE=0,eyeview_last=0;
  frameratevalue=1000;
  setpartmin=PERCENTILE_MIN, setpartmax=PERCENTILE_MAX, setslicemin=PERCENTILE_MIN, setslicemax=PERCENTILE_MAX, endian=0;
  setpartmin_old=setpartmin;
  setpartmax_old=setpartmax;
  setpatchmin=GLOBAL_MIN, setpatchmax=GLOBAL_MAX, setzonemin=0, setzonemax=0;
  loadpatchbysteps=0;
  settargetmin=0, settargetmax=0;
  setpartchopmin=0, setpartchopmax=0;
  setslicechopmin=0, setslicechopmax=0;
  partchopmin=1.0,  partchopmax=0.;
  slicechopmin=0, slicechopmax=0;
  setpatchchopmin=0, setpatchchopmax=0;
  patchchopmin=0,  patchchopmax=0;

  vis_onlyignited=0, vis_ignited=0, canshow_ignited=1, activate_ignited=1;
  settmin_p=0, settmin_b=0, settmin_s=0, settmin_z=0, settmin_i=0;
  settmax_p=0, settmax_b=0, settmax_s=0, settmax_z=0, settmax_i=0;
  set_no_part=0;
  tmin_p=1., tmin_b=1., tmin_s=1., tmin_z=1., tmin_i=1.;
  tmax_p=0., tmax_b=0., tmax_s=0., tmax_z=0., tmax_i=0.;
  patchmin=1., patchmax=0.;
  targetmin=1.0, targetmax=0.0;
  partmin=1., partmax=0., slicemin=1., slicemax=0.;
  zonemin=1.0, zonemax=0.0;
  speedmax=0.0;
  axissmooth=1,axisnum=1;
  FlowDir=1,ClipDir=1;
  plotn=1;
  stept=0;
  plotstate=NO_PLOTS;
  visVector=0;
  visSmokePart=2, visSprinkPart=1, havesprinkpart=0;
  visaxislabels=0;
  numplot3dvars=0;
  skip=1;
  p3dsurfacesmooth=1;
  p3dsurfacetype=1;
  parttype=0;
  allexterior=1,showexterior=1;
  allinterior=1;
  showbounds=0,showmotion=0,showedit=0, showclip=0, showgluistereo=0, showtour=0, showlabels=0, showcolorbar=0;
#ifdef pp_WUI
  showterrain=0;
#endif
  showgluitrainer=0;
  colorbarcycle=0;
  colorbartype=1;
  colorbartype_save=1;
  colorbarpoint=0;
  vectorspresent=0;

  visTarg = 0, ReadTargFile;
  showtarget=0;
  visAIso=0;
  surfincrement=0,visiso=0;
  isotest=0;
  isooffset=1,offsetmax=5;
  isolevelindex=0, isolevelindex2=0;
  pref=101325.,pamb=0.,tamb=288.;
  ntc_total=0, nspr_total=0, nheat_total=0;
  n_devices=0;

  npartinfo=0, nslice=0, nvslice=0, nslice2=0, npatch2=0, nplot3d=0, npatch_files=0;
  nevac=0;
  current_particle_type=-1,last_particle_type=-2;
  nsmoke3d=0;
  niso=0;
  ntrnx=0, ntrny=0, ntrnz=0,npdim=0,nmeshes=0,clip_mesh=0;
  nobst=0,nvent=0,noffset=0;
  nlabels=0,visLabels=0,nlabelssmv=0;
  ntarg_files=0;
  showallslicevectors=0;
  framerate=-1.0;
  ntimes=0, itime=0, itimeold=-999, seqnum=0,RenderTime=0;
  npqq=0, nopart=1;
  uindex=-1, vindex=-1, windex=-1;

  p3cont2d=STEPPED_CONTOURS, p3cont3dsmooth=0;
  cullfaces=1;
  showonly_hiddenfaces=0;


  windowresized=0;

  updatemenu=0;
  updatezoommenu=0;
  updatemenu_count=0;

  updatefaces=0,updatefacelists=0;
  updateOpenSMVFile=0;

  periodic_reloads=0;
  periodic_value=-2;

  slicefilenum=-1;
  partfilenum=-1,zonefilenum=-1;
  targfilenum=-1;

  setPDIM=0;
  menustatus=GLUT_MENU_NOT_IN_USE;
  visTimeZone=1, visTimeSmoke=1, visTimeSlice=1, visTimePatch=1, visTimeIso=1, visTimeEvac=1;
  vishmsTimelabel=0, visTimeLabels=1, visColorLabels=1;
  visTitle=1, visFullTitle=1, visFramerate=0, visFramelabel=1, visTimelabel=1;
#ifdef pp_memstatus
  visAvailmemory=0;
#endif
#ifdef pp_HRR
  visHRRlabel=0;
#endif
  visBlocklabel=1;
  visOpenVents=1,visDummyVents=1;
  visOpenVentsAsOutline=0;
  visTitle0=1, visTitle1=1, visTitle2=1;
  ntitles=0,ititle=0;
  visSmoke=1, visZone=0;
  visEvac=1;
  visBlocks=visBLOCKAsInput;
  visSmoothAsNormal=0;
  visBlocksSave=visBLOCKAsInput;
  blocklocation=BLOCKlocation_grid;
  ncadgeom=0;
  visFloor=0, visFrame=1;
  visNormalEditColors=1;
  visWalls=0, visGrid=0, visCeiling=0, cursorPlot3D=0;
  visVZone=1, visHZone=0;
  visSensor=1, visSensorNorm=1, hasSensorNorm=0, visSprink=1, visHeat=1;
  visVents=1;
  partframestep=1, sliceframestep=1, boundframestep=1;
  partframeskip=0, sliceframeskip=0, boundframeskip=0;
  boundzipstep=1, boundzipskip=0;
  smoke3dzipstep=1, smoke3dzipskip=0;
  evacframeskip=0, evacframestep=1;
  partpointstep=1;
  partpointstep_old=0;
  partpointskip=0;
  render_option=RenderWindow;
  RenderMenu(render_option);
  viewoption=0;
  xyz_clipplane=0;
  clip_x=0,clip_y=0,clip_z=0,clip_i=0,clip_j=0,clip_k=0;
  clip_X=0,clip_Y=0,clip_Z=0,clip_I=0,clip_J=0,clip_K=0;
  clip_x_val=0.0, clip_y_val=0.0, clip_z_val=0.0;
  clip_X_val=0.0, clip_Y_val=0.0, clip_Z_val=0.0;
  stepclip_x=0,stepclip_y=0,stepclip_z=0;
  stepclip_X=0,stepclip_Y=0,stepclip_Z=0;
  partpointsize=2.0,vectorpointsize=3.0,streaklinewidth=1.0;
  vectorpointsize=2.0, vectorlinewidth=1.0;
  sprinklerabssize=0.076f, sensorabssize=0.038f, heatabssize=0.076f;

  linewidth=2.0, ventlinewidth=2.0, highlight_linewidth=4.0;
  sliceoffset_factor=0.1f, ventoffset_factor=0.1f;

  nrgb=NRGB;
  nrgb_ini=0;
  nrgb2_ini=0;
  rgb_white=NRGB, rgb_yellow=NRGB+1, rgb_blue=NRGB+2, rgb_red=NRGB+3;
  rgb_green=NRGB+4, rgb_magenta=NRGB+5, rgb_cyan=NRGB+6, rgb_black=NRGB+7;
  numColorbars=0;
  setbw=0,colorbarflip=0;
  flip=1;
  transparentlevel=0.8f;
  transparentflag=1;
  transparentflagVOL=1;
  transparentflagSAVE=1;
  antialiasflag=1;
  nrgb_full=256;
  nrgb_cad=256;
  eyexfactor=0.5f, eyeyfactor=-0.9f, eyezfactor=0.5f;
  transparent_state=3;

  frameinterval=1.0;

  defaulttour_loaded=0;
  blockages_dirty=0;
  usetextures=0;
  canrestorelastview=0;
  ntargets=0;
  endian_data=0, endian_native=0, setendian=0;

  mainwindow_id=0,dwinWW=dwinW;
  rendertourcount=0;

  static_color[0]=0.0;
  static_color[1]=1.0;
  static_color[2]=0.0;
  static_color[3]=1.0;

  sensorcolor[0]=1.0;
  sensorcolor[1]=1.0;
  sensorcolor[2]=0.0;
  sensorcolor[3]=1.0;
  
  
  sensornormcolor[0]=1.0;
  sensornormcolor[1]=1.0;
  sensornormcolor[2]=0.0;
  sensornormcolor[3]=1.0;


  sprinkoncolor[0]=0.0;
  sprinkoncolor[1]=1.0;
  sprinkoncolor[2]=0.0;
  sprinkoncolor[3]=1.0;

  sprinkoffcolor[0]=1.0; //xxxx check
  sprinkoffcolor[1]=0.0;
  sprinkoffcolor[2]=0.0;
  sprinkoffcolor[3]=1.0;


  heatoncolor[0]=1.0; //xxx check
  heatoncolor[1]=0.0;
  heatoncolor[2]=0.0;
  heatoncolor[3]=1.0;

  heatoffcolor[0]=1.0;
  heatoffcolor[1]=0.0;
  heatoffcolor[2]=0.0;
  heatoffcolor[3]=1.0;

  backgroundbasecolor[0]=0.0;
  backgroundbasecolor[1]=0.0;
  backgroundbasecolor[2]=0.0;
  backgroundbasecolor[3]=1.0;

  backgroundcolor[0]=0.0;
  backgroundcolor[1]=0.0;
  backgroundcolor[2]=0.0;
  backgroundcolor[3]=1.0;

  foregroundbasecolor[0]=1.0;
  foregroundbasecolor[1]=1.0;
  foregroundbasecolor[2]=1.0;
  foregroundbasecolor[3]=1.0;

  foregroundcolor[0]=1.0;
  foregroundcolor[1]=1.0;
  foregroundcolor[2]=1.0;
  foregroundcolor[3]=1.0;

  boundcolor[0]=0.5;
  boundcolor[1]=0.5;
  boundcolor[2]=0.5;
  boundcolor[3]=1.0;

  timebarcolor[0]=0.6;
  timebarcolor[1]=0.6;
  timebarcolor[2]=0.6;
  timebarcolor[3]=1.0;
 
  redcolor[0]=1.0;
  redcolor[1]=0.0;
  redcolor[2]=0.0;
  redcolor[3]=1.0;

  loadfiles_at_startup=0;

  nmenus=0;
  showbuild=0;

  strcpy(TITLE1,"");
  strcpy(TITLE2,"");
  strcpy(TITLERELEASE,"");
  strcpy(TITLE,"");
  strcpy(TRAINERTITLE,"");
  strcpy(FULLTITLE,"");
  strcpy(emptylabel,"");
  large_font=GLUT_BITMAP_HELVETICA_12;
  small_font=GLUT_BITMAP_HELVETICA_10;

  veclengths[0]=0.25;
  veclengths[1]=0.40;
  veclengths[2]=0.63;
  veclengths[3]=1.00;
  veclengths[4]=1.59;
  veclengths[5]=2.52;
  veclengths[6]=4.00;

  texture_origin[0]=0.0;
  texture_origin[1]=0.0;
  texture_origin[2]=0.0;

  lock_allsmoke=0;

  visVentLines=0, visVentSolid=1;
  isZoneFireModel=0;
  output_slicedata=0,init_slicedata=1;
  nunitclasses=0,nunitclasses_default=0,nunitclasses_ini=0;
  callfrom_tourglui=0;
  showtours_whenediting=0;

  pbalanceindex=9, eoffsetindex=9;
  showstereo=0;
  showglui3dsmoke=0;
  showgluitour=0;
  showalert=0;
  stereoactive=0;
  apertureindex=2;
  zoomindex=2;
  projection_type=0;
  apertures[0]=30.;
  apertures[1]=45.;
  apertures[2]=60.;
  apertures[3]=75.;
  apertures[3]=90.;

  aperture=60.,aperture_glui,aperture_default;
  zooms[0]=0.25;
  zooms[1]=0.5;
  zooms[2]=1.0;
  zooms[3]=2.0;
  zooms[4]=4.0;
  zoom=1.0;

  {
    int ii;
    rgbmask[0]=1;
    for(ii=1;ii<16;ii++){
      rgbmask[ii]=2*rgbmask[ii-1]+1;
    }
  }

  sv_age=0;
#ifdef pp_NISTREVIEW
  nistreview=1;
  atNIST=0;
#else
  nistreview=0;
  atNIST=1;
#endif
  titlesafe_offset=0;
  titlesafe_offsetBASE=45;
  reset_frame=0;
  reset_time=0.0,start_frametime=0.0,stop_frametime=0.0;
  reset_time_flag=0;
  velocity_range=0;
  RenderOnceNow=0, RenderOnceNowR=0, RenderOnceNowL=0;

  pass_through=0;
  nsorted_surfidlist=0;

  overwrite_all=0,erase_all=0;
  strcpy(ext_png,".png");
  strcpy(ext_jpg,".jpg");
#ifdef pp_GDGIF
  strcpy(ext_gif,".gif");
#endif
  renderfiletype=0;
  strcpy(part_ext,".part");
  strcpy(ini_ext,".ini");

  updatehiddenfaces=0;
  isShell=0;

  start_xyz0[0]=0.0;
  start_xyz0[1]=0.0;
  start_xyz0[2]=0.0;
  glui_move_mode=-1;

  timeoffset=0.0;
  motion_factor=1.0,speed_factor=1.0;
  speed_desired=1;
  speed_I=-1;
  speed_crawl=0.22,speed_walk=0.45,speed_now=0.447;
  status_now=0,old_status_now=-1;
  motion_flag=0;
  update_tourlist=0;
  desired_view_height=1.5;
  resetclock=1,initialtime=0;
  realtime_flag=0;
  islicetype=-1,islicetype_save=-1,ipatchtype=-1;
  iisotype=-1;


  cpuframe=0;

  direction_angleINI=0.0;

  direction_angleINI0=0.0;

  adjustalphaflag=2;

  highlight_block=-1, highlight_mesh=0, highlight_flag=2;
  updatesmoothblocks=1;
  menusmooth=0;
  use_menusmooth=0;
  smoothing_blocks=0;
  blocksneedsmoothing=0;
  updategetlabels=1;

  pixel_skip=0;
  smoke_extinct=7.600,smoke_dens=.50,smoke_pathlength=1.0;
#ifdef pp_SMOKETEST
  smoketest=1,show_smoketest=1;
#else
  smoketest=0,show_smoketest=0;
#endif
  showall_textures=0;

  ncolorbars=0;
  ndefaultcolorbars=3;
  update_load_startup=0;
  do_ignited=0;
#ifdef pp_THREADS2
  mt_compress=1;
#else
  mt_compress=0;
#endif
  updateindexcolors=0;
  show_path_knots=0;
  keyframe_snap=0;
  tourviewtype=0;
  show_tourlocus=1;
  tourlocus_type=0;
  tourrad_avatar=0.1;
  dirtycircletour=0;
  view_tstart=0.0, view_tstop=100.0;
  tour_constant_vel=0;
  tour_bias=0.0,tour_continuity=0.0;
  view_ntimes=1000;
  ntours=0,selectedtour_index=-1,selectedtour_index_old=-1,selectedtour_index_ini=-1;
  update_selectedtour_index=0;
  viewtourfrompath=0,viewalltours=0,viewanytours=0,edittour=0;
  tour_usecurrent=0;
  rotation_index_OLD=-1;
  nticks=0,ntickssmv=0;
  visTicks=0;
  visCadTextures=1;
  cb_hidesv=0;
#ifdef pp_COLOR
  viscolorbarpath=0;
  showzerosplit=1;
#else
  viscolorbarpath=0;
  showzerosplit=0;
#endif
  nselectblocks=0;
  surface_indices[0]=0.0;
  surface_indices[1]=0.0;
  surface_indices[2]=0.0;
  wall_case=0;
  strcpy(surfacedefaultlabel,"");
  mscale[0]=1.0;
  mscale[1]=1.0;
  mscale[2]=1.0;
  nearclip=0.001,farclip=3.0;
  updateclipvals=0;
  updateUpdateFrameRateMenu=0;
  ntextures=0;
  nskyboxinfo=0;
  part5show=1;
  streak5show=0;
//SVEXTERN int nstreak_value; // 5
//SVEXTERN char *streak_values[5]; // "1","2","4","8","16"
//SVEXTERN float streak_rvalue[5]; // 1.0, 2.0 4.0, 8.0, 16.0 
//SVEXTERN int streak_index;       // 0
//SVEXTERN float float_streak5value;// 1.0

  streak_rvalue[0]=0.25;
  streak_rvalue[1]=0.5;
  streak_rvalue[2]=1.0;
  streak_rvalue[3]=2.0;
  streak_rvalue[4]=4.0;
  streak_rvalue[5]=8.0;
  streak_rvalue[6]=16.0;
  nstreak_value=7;
  streak_index=-1;
  float_streak5value=0.0;
  if(streak_index>=0)float_streak5value=streak_rvalue[streak_index];

  streak5step=0;
  showstreakhead=1;
  npartclassinfo=0;
  prop_index=1;
  dummyvents=0;
  noutlineinfo=0;
  nmultislices=0;
  nmultivslices=0;

  svofile_exists=0;
  devicenorm_length = 0.1;
  ndeviceinfo=0;
  ndevice_defs=0;
  strcpy(device_def_first.label,"first");
  device_def_first.next=&device_def_last;
  device_def_first.prev=NULL;

  strcpy(device_def_last.label,"last");
  device_def_last.next=NULL;
  device_def_last.prev=&device_def_first;
  device_defs=NULL;
 
  showfiles=0;
  smoke3d_external=0;
 
  ncases=0,case_number=0;
  smoke_shade=0, fire_red=255, fire_green=128, fire_blue=0;
  fire_halfdepth=2.0;
  hrrpuv_cutoff=600.0;

  smokecullflag=1;
  smokedrawtest=0,smokedrawtest2=0;
  visMAINmenus=0;
  smoke3d_thick=0;
  smokedrawtest_nummin=1;
  smokedrawtest_nummax=1;
  ijkbarmax=5;
  blockage_as_input=0;
  show_cad_and_grid=0;
  use_nistlogo=0;
  benchmark_flag=0;
  nplot3dtimelist=0;

  buffertype=DOUBLE_BUFFER;
  benchmark=0;
  opengldefined=0;
  transparency_level=0.0;
  transparency_override=0;

  dwinHbase=60;
  dwinH=60;

#ifndef pp_cvf
  strcpy(dirseparator,"/");
#else
  strcpy(dirseparator,"\\");
#endif

  {
    char svn[1024];
    char *svnnum;
#ifdef pp_TEST
    char sv_version[100]="5.1.x";
#else
    char sv_version[100]="5.0.0";
#endif
  
    strcpy(svn,"$Revision$");

    svnnum=strchr(svn,':');
    if(svnnum!=NULL&&strlen(svnnum)>4){
      svnnum++;
      svnnum=trim_front(svnnum);
      svnnum[strlen(svnnum)-1]=0;
      trim(svnnum);
      strcat(sv_version,"_");
      strcat(sv_version,svnnum);
    }

    strcpy(TITLEBASE,"Smokeview ");

#define pp_BETA
#ifdef pp_BETA
    strcat(TITLEBASE," Beta ");
#endif
#ifdef pp_TEST
    strcat(TITLEBASE," Test");
#endif
    strcat(TITLEBASE,sv_version);
    strcat(TITLEBASE," - ");

    strcpy(TRAINERTITLEBASE,"Smokeview Demonstrator");

#ifdef pp_BETA
    strcat(TRAINERTITLEBASE," Beta ");
#endif
#ifdef pp_TEST
    strcat(TRAINERTITLEBASE," Test");
#endif
    strcat(TRAINERTITLEBASE,sv_version);
    strcat(TRAINERTITLEBASE," - ");

  }

  strcpy(INIfile,"smokeview.ini");
  strcpy(WRITEINIfile,"Write smokeview.ini");

  tourcol_selectedpathline[0]=1.0;
  tourcol_selectedpathline[1]=0.0;
  tourcol_selectedpathline[2]=0.0;


  tourcol_selectedpathlineknots[0]=1.0;
  tourcol_selectedpathlineknots[1]=0.0;
  tourcol_selectedpathlineknots[2]=0.0;


  tourcol_selectedknot[0]=0.0;
  tourcol_selectedknot[1]=1.0;
  tourcol_selectedknot[2]=0.0;


  tourcol_selectedview[0]=1.0;
  tourcol_selectedview[1]=1.0;
  tourcol_selectedview[2]=0.0;


  tourcol_pathline[0]=-1.0;
  tourcol_pathline[1]=-1.0;
  tourcol_pathline[2]=-1.0;

  tourcol_pathknots[0]=-1.0;
  tourcol_pathknots[1]=-1.0;
  tourcol_pathknots[2]=-1.0;

  tourcol_text[0]=-1.0;
  tourcol_text[1]=-1.0;
  tourcol_text[2]=-1.0;


  tourcol_avatar[0]=1.0;
  tourcol_avatar[1]=0.0;
  tourcol_avatar[2]=0.0;

  iso_specular[0] = 0.7;
  iso_specular[1] = 0.7;
  iso_specular[2] = 0.7;
  iso_specular[3] = 1.7;
  iso_shininess = 10.0f;

  block_shininess = 100.;

  light_position0[0] = 1.0f;
  light_position0[1] = 1.0f;
  light_position0[2] = 1.0f;
  light_position0[3] = 0.0f;

  light_position1[0] = -1.0f;
  light_position1[1] = -1.0f;
  light_position1[2] =  1.0f;
  light_position1[3] =  0.0f;


  ambientlight[0] = 0.6f;
  ambientlight[1] = 0.6f;
  ambientlight[2] = 0.6f;
  ambientlight[3] = 1.0f;


  diffuselight[0] = 0.50f;
  diffuselight[1] = 0.50f;
  diffuselight[2] = 0.50f;
  diffuselight[3] = 1.00f;


  list_p3_index_old=0, list_slice_index_old=0, list_patch_index_old=0;

  videoSTEREO=0;
  stereo_frame=0;
  stereo_leftright=0;
  stereo_off=1;

  pbalance=-0.5, eoffset=0.7;
  pbalanceORIG=-0.5, eoffsetORIG=0.7;
  
   pbalances[0]=0.000;
   pbalances[1]=0.500;
   pbalances[2]=0.545;
   pbalances[3]=0.595;
   pbalances[4]=0.648;
   pbalances[5]=0.707;
   pbalances[6]=0.771;
   pbalances[7]=0.841;
   pbalances[8]=0.917;
   pbalances[9]=1.000;
  pbalances[10]=1.090;
  pbalances[11]=1.190;
  pbalances[12]=1.300;
  pbalances[13]=1.410;
  pbalances[14]=1.540;
  pbalances[15]=1.680;
  pbalances[16]=1.830;
  pbalances[17]=2.000;
  
   eoffsets[0]=0.000;
   eoffsets[1]=0.500;
   eoffsets[2]=0.545;
   eoffsets[3]=0.595;
   eoffsets[4]=0.648;
   eoffsets[5]=0.707;
   eoffsets[6]=0.771;
   eoffsets[7]=0.841;
   eoffsets[8]=0.917;
   eoffsets[9]=1.000;
  eoffsets[10]=1.090;
  eoffsets[11]=1.190;
  eoffsets[12]=1.300;
  eoffsets[13]=1.410;
  eoffsets[14]=1.540;
  eoffsets[15]=1.680;
  eoffsets[16]=1.830;
  eoffsets[17]=2.000;

  strcpy(blank,"");

  demo_mode=0;
  update_demo=1;
  menu_view_number=0;
  mxplot3dvars=MAXPLOT3DVARS;

  valindex=0;

   iso_ambient[0] = 0.96;
   iso_ambient[1] = 0.28;
   iso_ambient[2] = 0.00;
   iso_ambient[3] = 1.00;
   iso_ambient[4] = 0.75;
   iso_ambient[5] = 0.80;
   iso_ambient[6] = 0.80;
   iso_ambient[7] = 3,00;
   iso_ambient[8] = 0.00;
   iso_ambient[9] = 0.96;
  iso_ambient[10] = 0.28;
  iso_ambient[11] = 0.80;


  iso_transparency=0.8;
  n_iso_ambient=3;
  n_iso_ambient_ini=0;
  mouse_deltax=0.0, mouse_deltay=0.0;

  char_color[0]=0.0;
  char_color[1]=0.0;
  char_color[2]=0.0;
  char_color[3]=0.0;

  movedir[0]=0.0;
  movedir[1]=1.0;
  movedir[2]=0.0;

  memcpy(rgb_base,rgb_baseBASE,MAXRGB*4*sizeof(float));
  memcpy(bw_base,bw_baseBASE,MAXRGB*4*sizeof(float));
  memcpy(rgb2,rgb2BASE,MAXRGB*3*sizeof(float));
  memcpy(bw_base,bw_baseBASE,MAXRGB*4*sizeof(float));
  
  nrgb2=8;

  ncamera_list=0,i_view_list=1,init_camera_list_flag=1;
  camera_max_id=2;
  startup=0,startup_view_ini=1,selected_view=-999;
  

  {
    int iii;

    for(iii=0;iii<7;iii++){
      visPatchType[iii]=0;
    }
    visPatchType[0]=1;
    for(iii=0;iii<MAXPLOT3DVARS;iii++){
      setp3min[iii]=PERCENTILE_MIN;
      p3min[iii]=1.0f;
      p3chopmin[iii]=1.0f;
      setp3max[iii]=PERCENTILE_MAX;
      p3max[iii]=1.0f;
      p3chopmax[iii]=-.0f;
    }
  }

}

/* ------------------ initvars0 ------------------------ */

void initvars0(void){

#ifdef pp_GPU
  GPU_depthtexture=0;
#endif

  camera_current=NULL, camera_save=NULL, camera_last=NULL;
  camera_external=NULL, camera_internal=NULL, camera_ini=NULL;
  camera_list=NULL;

  camera_label=NULL;

  iso_ambient_ini=NULL;
  rgb_ini=NULL;

  sphere_xyz=NULL;
  LESsystem=NULL,LESendian=NULL;

  zonet=NULL, zoneylay=NULL, zonetl=NULL, zonetu=NULL, zonepr=NULL, zoneqfire=NULL;
  hazardcolor=NULL;
  izonetu=NULL;
  tspr=NULL;

#ifdef pp_SPHERE
  sphereinfo=NULL;
#endif

  tour_t=NULL, tour_t2=NULL, tour_dist=NULL, tour_dist2=NULL, tour_dist3=NULL;
  colorbarinfo=NULL;

  selectfaceinfo=NULL;
  selectblockinfo=NULL;
  tickinfo=NULL;
  firstcolor=NULL;
  textureinfo=NULL;

  sortedblocklist=NULL,changed_idlist=NULL;
  surfaceinfo=NULL;
  skyboxinfo=NULL;
  fireinfo=NULL;
  roominfo=NULL;
  zventinfo=NULL;
  zoneinfo=NULL;
  activezone=NULL;
  partinfo=NULL;
   current_property=NULL;
   partclassinfo=NULL;
   part5propinfo=NULL;
   npart5prop=0;
  targinfo=NULL;
  sliceinfo=NULL;
  caminfo=NULL;
  multisliceinfo=NULL;
  multivsliceinfo=NULL;
  outlineinfo=NULL;
  sliceorderindex=NULL,vsliceorderindex=NULL,partorderindex=NULL;
  patchorderindex=NULL,isoorderindex=NULL,plot3dorderindex=NULL;
  slicebounds=NULL;
  vsliceinfo=NULL;
  smoke3dinfo=NULL;
  caseinfo=NULL,selected_case=NULL;
  labelinfo=NULL;
  slicetypes=NULL, isotypes=NULL, vslicetypes=NULL, patchtypes=NULL;
  plot3dinfo=NULL;
  plot3dtimelist=NULL;
  patchinfo=NULL;
  isoinfo=NULL;
  target_positions=NULL;

  bchighlight=NULL,bchighlight_old=NULL;
  cadgeominfo=NULL;

  ventcolor=NULL;
  meshinfo=NULL,current_mesh=NULL, mesh_save=NULL, mesh_last=NULL, loaded_isomesh=NULL;
  tourinfo=NULL,default_tour;
#ifdef pp_WUI
  treeinfo=NULL;
  ntreeinfo=0;
  terraininfo=NULL;
  nterraininfo=0;
  visTerrain=0;
  treecolor[0]=0.09;
  treecolor[1]=0.5;
  treecolor[2]=0.015;
  treecolor[3]=1.0;
  treecharcolor[0]=0.3;
  treecharcolor[1]=0.3;
  treecharcolor[2]=0.3;
  treecharcolor[3]=1.0;
  trunccolor[0]=0.6;
  trunccolor[1]=0.2;
  trunccolor[2]=0.0;
  trunccolor[3]=1.0;
#endif
  tourknotskeylist=NULL;
  tourknotstourlist=NULL;
  selected_frame=NULL;
  selected_tour=NULL;
  unitclasses=NULL,unitclasses_default=NULL,unitclasses_ini=NULL;
  smokeviewbindir=NULL;
  smokeviewtempdir=NULL;
  partshortlabel=NULL,partunitlabel=NULL;
  texturedir=NULL;
  targtimes=NULL;
  targtimeslist=NULL;
  zonetlist=NULL;
  slice_loaded_list=NULL;
  render_frame=NULL;
  fdsprefix=NULL, fdsprefix2=NULL;
  endianfilename=NULL;
  targfilename=NULL;
  sorted_surfidlist=NULL,inv_sorted_surfidlist=NULL;
  trainer_filename=NULL;
  smvfilename=NULL, smvmenufile=NULL,databasefilename=NULL,smvprogdir=NULL;
  flushfile=NULL, chidfilebase=NULL;
#ifdef pp_HRR
  hrrfilename=NULL;
  hrrinfo=NULL;
#endif
  smokezippath=NULL;
  shellfilename=NULL;
  INI_fds_filein=NULL, fds_filein=NULL, fds_fileout=NULL,fds_fileout2=NULL;
  casefilename=NULL;
  zonelonglabels=NULL, zoneshortlabels=NULL, zoneunits=NULL;
  smokeviewini=NULL;
  surfids=NULL;
  colorlabelpart=NULL;
  colorlabelpatch=NULL;
  colorlabelzone=NULL;
  p3levels=NULL, zonelevels=NULL;
  p3levels256=NULL;
  colorlabelp3=NULL;
  partscale=NULL;
  zonescale=NULL;
  plotiso=NULL;
  times=NULL;
  patchlabellist=NULL;
  sliceindex=NULL;
  face_transparent=NULL;
  deviceinfo=NULL;
}

void init_default_devices(void){


  /*
  devices
  -------
  1 - sensor
  2 - sprinkler
  3 - heat detector
  4 - smoke detector
  */



}
