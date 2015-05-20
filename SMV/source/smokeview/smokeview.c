#include "options.h"
#include <stdio.h>  
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "glew.h"
#include GLUT_H

#include "smokeviewvars.h"

#ifdef WIN32
#include <direct.h>
#endif

/* ------------------ snifferrors ------------------------ */

void _Sniff_Errors(char *whereat){
  int error;

  while((error=glGetError())!=GL_NO_ERROR){
    char *glu_error;

    glu_error=(char *)gluErrorString((unsigned int)error);
    fprintf(stderr,"*** Error: OpenGL error:%s, where:%s %i\n",
      glu_error,whereat,snifferrornumber);
      snifferrornumber++;
  }
}

/* ------------------ updateLights ------------------------ */

void updateLights(float *pos1, float *pos2){
  int i;
  GLfloat ambientlight2[4], diffuselight2[4];  
  int lightCount;
  float div;
        
  glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, lightmodel_localviewer == 0? GL_FALSE : GL_TRUE);
  glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, lightmodel_separatespecularcolor == 0? GL_SINGLE_COLOR : GL_SEPARATE_SPECULAR_COLOR);
  
  lightCount = 0;
  if(light_enabled0){
    ++lightCount;
  }
  if(light_enabled1){
    ++lightCount;
  }

  div = lightCount > 0? 1.0f/(float)lightCount : 1.0f;  
  for(i=0;i<3;i++){
    ambientlight2[i]=ambientlight[i]*div;
    diffuselight2[i]=diffuselight[i]*div;
  }
  ambientlight2[3]=1.0;
  diffuselight2[3]=1.0;
  if(light_enabled0){  
    glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuselight2);
    glLightfv(GL_LIGHT0,GL_AMBIENT,ambientlight2);
    if(pos1!=NULL)glLightfv(GL_LIGHT0,GL_POSITION,pos1);
    glEnable(GL_LIGHT0);
  }
  else{
    glDisable(GL_LIGHT0);
  }

  if(light_enabled1){  
    glLightfv(GL_LIGHT1,GL_DIFFUSE,diffuselight2);
    glLightfv(GL_LIGHT1,GL_AMBIENT,ambientlight2);
    if(pos2!=NULL)glLightfv(GL_LIGHT1,GL_POSITION,pos2);
    glEnable(GL_LIGHT1);
  }
  else{
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

/* ------------------ smv2quat ------------------------ */

void camera2quat(camera *ca, float *quat, float *rotation){
  if(ca->quat_defined==1){
    quat[0]=ca->quaternion[0];
    quat[1]=ca->quaternion[1];
    quat[2]=ca->quaternion[2];
    quat[3]=ca->quaternion[3];
  }
  else{
    float quat_temp[4];
    float azimuth, elevation,axis[3];

    azimuth = ca->az_elev[0]*DEG2RAD;
    elevation = (ca->az_elev[1])*DEG2RAD;

    axis[0]=1.0;
    axis[1]=0.0;
    axis[2]=0.0;

    angleaxis2quat(elevation,axis,quat_temp);

    axis[0]=0.0;
    axis[1]=0.0;
    axis[2]=1.0;

    angleaxis2quat(azimuth,axis,quat);

    mult_quat(quat_temp,quat,quat);
  }

  if(rotation!=NULL)quat2rot(quat,rotation);
}


/* ------------------ ResetView ------------------------ */

void ResetView(int option){
  in_external=0;
  switch(option){
    int rotation_type_save;
    int projection_type_save;

  case RESTORE_EXTERIOR_VIEW_ZOOM:
    break;
  case RESTORE_EXTERIOR_VIEW:
    in_external=1;
    rotation_type_save = camera_current->rotation_type;
    projection_type_save = camera_current->projection_type;
    copy_camera(camera_current,camera_external);
    camera_current->rotation_type=rotation_type_save;
    camera_current->projection_type=projection_type_save;
    if(camera_current->projection_type==1){
      camera_current->eye[1]=camera_current->isometric_y;
    }
    break;
  case RESTORE_INTERIOR_VIEW:
    rotation_type_save = camera_current->rotation_type;
    projection_type_save = camera_current->projection_type;
    copy_camera(camera_current,camera_internal);
    camera_current->rotation_type=rotation_type_save;
    camera_current->projection_type=projection_type_save;
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
  if(rotation_type==ROTATION_3AXIS){
    float azimuth, elevation,axis[3];
    float quat_temp[4];
    float x, y, z;
   
    azimuth = camera_current->az_elev[0]*DEG2RAD;
    elevation = camera_current->az_elev[1]*DEG2RAD;

    x = cos(azimuth);
    y = sin(azimuth);
    z = cos(elevation);

    axis[0]=0.0;
    axis[1]=0.0;
    axis[2]=1.0;

    angleaxis2quat(azimuth,axis,quat_temp);

    axis[0]=x;
    axis[1]=y;
    axis[2]=0.0;

    angleaxis2quat(acos(z),axis,quat_general);

    mult_quat(quat_temp,quat_general,quat_general);

    quat2rot(quat_general,quat_rotation);
  }
  if(option==RESTORE_EXTERIOR_VIEW_ZOOM)camera_current->zoom=zooms[zoomindex];
  zoom=camera_current->zoom;
  update_glui_zoom();
}

/* ------------------ init_volrender_script ------------------------ */

void init_volrender_script(char *prefix, char *tour_label, int startframe, int skipframe){
  scriptfiledata *sfd;
  FILE *script_stream;

  if(volrender_scriptname==NULL){
    int len;

    len = strlen(fdsprefix)+strlen("_volrender.ssf")+1;
    NewMemory((void **)&volrender_scriptname,(unsigned int)(len));
    STRCPY(volrender_scriptname,fdsprefix);
    STRCAT(volrender_scriptname,"_volrender.ssf");
  }

  sfd = insert_scriptfile(volrender_scriptname);
  if(sfd!=NULL)default_script=sfd;
  script_stream=fopen(volrender_scriptname,"w");
  if(script_stream!=NULL){
    fprintf(script_stream,"RENDERDIR\n");
    fprintf(script_stream," .\n");
    if(tour_label!=NULL&&strcmp(tour_label,"Manual")!=0){
      fprintf(script_stream,"LOADTOUR\n");
      fprintf(script_stream," %s\n",tour_label);
    }
    fprintf(script_stream,"VOLSMOKERENDERALL\n");
    fprintf(script_stream," %i %i\n",skipframe,startframe);
    fprintf(script_stream," %s\n",prefix);
    runscript=1;
    fclose(script_stream);
  }
} 

/* ------------------ parse_commandline ------------------------ */

void parse_commandline(int argc, char **argv){
  int i, len_casename;
  int iarg;
  size_t len_memory;
  char *argi;
  char SMVFILENAME[1024];
  int smv_parse;

  CheckMemory;
  partscale=a_partscale;
  zonescale=a_zonescale;

  if(argc==1){
    exit(1);
  }
  if(strncmp(argv[1],"-ini",3)==0){
    init_camera_list();
    InitOpenGL();
    UpdateRGBColors(COLORBAR_INDEX_NONE);
    writeini(GLOBAL_INI,NULL);
    exit(0);
  }

  if(strncmp(argv[1],"-ng_ini",6)==0){
    init_camera_list();
    use_graphics=0;
    UpdateRGBColors(COLORBAR_INDEX_NONE);
    writeini(GLOBAL_INI,NULL);
    exit(0);
  }
  strcpy(SMVFILENAME,"");
  smv_parse=0;
  for(iarg=1;iarg<argc;iarg++){
    argi=argv[iarg];
    if(strncmp(argi,"-",1)==0){
      if(
        strncmp(argi,"-points",7)==0||
        strncmp(argi,"-frames",7)==0||
#ifdef pp_LANG        
        strncmp(argi,"-lang",5)==0||
#endif        
        strncmp(argi,"-script",7)==0||
        strncmp(argi,"-startframe",11)==0||
        strncmp(argi,"-skipframe",10)==0||
        strncmp(argi,"-bindir",7)==0||
        strncmp(argi,"-update_ini",11)==0
        ){
        iarg++;
      }
      if(strncmp(argi,"-convert_ini",12)==0)iarg+=2;

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
  len_casename = (int) strlen(argi);
  CheckMemory;
  FREEMEMORY(fdsprefix);
  len_memory=len_casename+strlen(part_ext)+100;
  NewMemory((void **)&fdsprefix,(unsigned int)len_memory);
  STRCPY(fdsprefix,argi);
  strcpy(movie_name, fdsprefix);
  strcpy(render_file_base, fdsprefix);
  FREEMEMORY(trainer_filename);
  FREEMEMORY(test_filename);
  FREEMEMORY(smoothblockage_filename);

  strcpy(input_filename_ext,"");

  if(len_casename>4){
    char *c_ext;

    c_ext=strrchr(argi,'.');
    if(c_ext!=NULL){
      STRCPY(input_filename_ext,c_ext);
      to_lower(input_filename_ext);

      if(c_ext!=NULL&&
        (strcmp(input_filename_ext,".smv")==0||
         strcmp(input_filename_ext,".svd")==0||
         strcmp(input_filename_ext,".smt")==0)
         ){
        c_ext[0]=0;
        STRCPY(fdsprefix,argi);
        strcpy(movie_name, fdsprefix);
        strcpy(render_file_base, fdsprefix);
        FREEMEMORY(trainer_filename);
        NewMemory((void **)&trainer_filename,(unsigned int)(len_casename+7));
        STRCPY(trainer_filename,argi);
        STRCAT(trainer_filename,".svd");
        FREEMEMORY(test_filename);
        NewMemory((void **)&test_filename,(unsigned int)(len_casename+7));
        STRCPY(test_filename,argi);
        STRCAT(test_filename,".smt");
      }
    }
  }

  FREEMEMORY(log_filename);
  NewMemory((void **)&log_filename,len_casename+7+1);
  STRCPY(log_filename,fdsprefix);
  STRCAT(log_filename,".smvlog");

  FREEMEMORY(caseini_filename);
  NewMemory((void **)&caseini_filename,len_casename+strlen(ini_ext)+1);
  STRCPY(caseini_filename,fdsprefix);
  STRCAT(caseini_filename,ini_ext);

  FREEMEMORY(boundini_filename);
  NewMemory((void **)&boundini_filename,len_casename+5+1);
  STRCPY(boundini_filename,fdsprefix);
  STRCAT(boundini_filename,".bini");

  if(smv_filename==NULL){
    STRUCTSTAT statbuffer;

    NewMemory((void **)&smv_filename,(unsigned int)(len_casename+6));
    STRCPY(smv_filename,fdsprefix);
    STRCAT(smv_filename,".smv");
    {
      char scriptbuffer[1024];

      STRCPY(scriptbuffer,fdsprefix);
      STRCAT(scriptbuffer,".ssf");
      if(default_script==NULL&&STAT(scriptbuffer,&statbuffer)==0){
        default_script = insert_scriptfile(scriptbuffer);
      }
    }
  }
  if(smv_filename!=NULL){
    STRUCTSTAT statbuffer;

    FREEMEMORY(fds_filein);
    NewMemory((void **)&fds_filein,strlen(fdsprefix)+6);
    STRCPY(fds_filein,fdsprefix);
    STRCAT(fds_filein,".fds");
    if(STAT(fds_filein,&statbuffer)!=0){
      FREEMEMORY(fds_filein);
    }
  }
  if(fed_filename==NULL){
    STRCPY(fed_filename_base,fdsprefix);
    STRCAT(fed_filename_base,".fed_smv");
    fed_filename=get_filename(smokeviewtempdir,fed_filename_base,tempdir_flag);
  }
  if(stop_filename==NULL){
    NewMemory((void **)&stop_filename,(unsigned int)(len_casename+6));
    STRCPY(stop_filename,fdsprefix);
    STRCAT(stop_filename,".stop");
  }
  if(sliceinfo_filename==NULL){
    NewMemory((void **)&sliceinfo_filename,strlen(fdsprefix)+11+1);
    STRCPY(sliceinfo_filename,fdsprefix);
    STRCAT(sliceinfo_filename,"_slice.info");
  }

  // if smokezip created part2iso files then concatenate .smv entries found in the .isosmv file 
  // to the end of the .smv file creating a new .smv file.  Then read in that .smv file.

  {
    FILE *stream_iso=NULL;

    NewMemory((void **)&iso_filename,len_casename+7+1);
    STRCPY(iso_filename,fdsprefix);
    STRCAT(iso_filename,".isosmv");
    stream_iso=fopen(iso_filename,"r");
    if(stream_iso!=NULL){
      fclose(stream_iso);
    }
    else{
      FREEMEMORY(iso_filename);
    }
  }

  if(trainer_filename==NULL){
    NewMemory((void **)&trainer_filename,(unsigned int)(len_casename+6));
    STRCPY(trainer_filename,fdsprefix);
    STRCAT(trainer_filename,".svd");
  }
  if(test_filename==NULL){
    NewMemory((void **)&test_filename,(unsigned int)(len_casename+6));
    STRCPY(test_filename,fdsprefix);
    STRCAT(test_filename,".svd");
  }
  if(smoothblockage_filename==NULL){
    char filename_base[1024];

    STRCPY(filename_base,fdsprefix);
    STRCAT(filename_base,".sb");
    smoothblockage_filename=get_filename(smokeviewtempdir,filename_base,tempdir_flag);
  }

  for (i=1;i<argc;i++){
    if(strncmp(argv[i],"-",1)!=0)continue;
    if(strncmp(argv[i],"-ini",3)==0){
      writeini(GLOBAL_INI,NULL);
    }
    else if(strncmp(argv[i],"-update_bounds",14)==0){
      use_graphics=0;
      update_bounds=1;
    }
    else if(strncmp(argv[i],"-nogpu",6)==0){
      disable_gpu=1;
    }
    else if(strncmp(argv[i],"-demo",5)==0){
       demo_option=1;
    }
    else if(strncmp(argv[i],"-stereo",7)==0){
      stereoactive=1;
      showstereo=STEREO_TIME;
      PRINTF("stereo option activated\n");
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
        show_lang_menu=1;
      }
    }
#endif    
    else if(strncmp(argv[i],"-convert_ini",12)==0){
      char *local_ini_from=NULL, *local_ini_to=NULL;

      if(++i<argc)local_ini_from = argv[i];
      if(++i<argc)local_ini_to = argv[i];
      if(local_ini_from!=NULL&&local_ini_to!=NULL){
        NewMemory((void **)&ini_from,strlen(local_ini_from)+1);
        strcpy(ini_from,local_ini_from);
        
        NewMemory((void **)&ini_to,strlen(local_ini_to)+1);
        strcpy(ini_to,local_ini_to);
        convert_ini=1;
      }
    }
    else if(strncmp(argv[i],"-update_ini",11)==0){
      char *local_ini_from=NULL, *local_ini_to=NULL;

      if(++i<argc)local_ini_from = argv[i];
      local_ini_to = local_ini_from;
      if(local_ini_from!=NULL){
        NewMemory((void **)&ini_from,strlen(local_ini_from)+1);
        strcpy(ini_from,local_ini_from);
        
        NewMemory((void **)&ini_to,strlen(local_ini_to)+1);
        strcpy(ini_to,local_ini_to);
        convert_ini=1;
      }
    }
    else if(strncmp(argv[i],"-isotest",8)==0){
      isotest=1;
    }
#ifdef _DEBUG
    else if(strncmp(argv[i],"-tempdir",8)==0){
      tempdir_flag=1;
    }
#endif
    else if(strncmp(argv[i],"-time",5)==0){
      time_flag=1;
    }
    else if(strncmp(argv[i],"-h",2)==0){
      usage(argv);
      exit(0);
    }
    else if(strncmp(argv[i],"-noblank",8)==0){
      arg_iblank=1;
      use_iblank=0;
    }
    else if(strncmp(argv[i],"-fed",4)==0){
      compute_fed=1;
    }
    else if(strncmp(argv[i],"-blank",6)==0){
      arg_iblank=1;
      use_iblank=1;
    }
    else if(strncmp(argv[i],"-gversion",9)==0){
      gversion=1;
    }
    else if(
      strncmp(argv[i],"-volrender",10)!=0&&(strncmp(argv[i],"-version",8)==0||strncmp(argv[i],"-v",2)==0)
      ){
      display_version_info();
      exit(0);
    }
    else if(
      strncmp(argv[i],"-redirect",9)==0
      ){
        LOG_FILENAME=fopen(log_filename,"w");
        if(LOG_FILENAME!=NULL){
          redirect=1;
          set_stdout(LOG_FILENAME);
        }
    }
    else if(strncmp(argv[i],"-runscript",10)==0){
      from_commandline=1;
      runscript=1;
    }
    else if(strncmp(argv[i],"-skipframe",10)==0){
      from_commandline=1;
      ++i;
      if(i<argc){
        sscanf(argv[i],"%i",&skipframe0);
      }
    }
    else if(strncmp(argv[i],"-startframe",11)==0){
      from_commandline=1;
      ++i;
      if(i<argc){
        sscanf(argv[i],"%i",&startframe0);
      }
    }
    else if(strncmp(argv[i],"-volrender",10)==0){
      from_commandline=1;
      make_volrender_script=1;
    }
    else if(strncmp(argv[i],"-script",7)==0){
      from_commandline=1;
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
    else if(strncmp(argv[i],"-noexit",7)==0){
      noexit=1;
    }
    else if(strncmp(argv[i],"-bindir",7)==0){
      ++i;
      if(i<argc){
        int len2;

        len2 = strlen(argv[i]);
        NewMemory((void **)&smokeview_bindir,len2+2);
        strcpy(smokeview_bindir,argv[i]);
        if(smokeview_bindir[len2-1]!=dirseparator[0])strcat(smokeview_bindir,dirseparator);
      }
    }
    else if(strncmp(argv[i],"-build",6)==0){
      showbuild=1;
      usage(argv);
      exit(0);
    }
    else {
      fprintf(stderr,"*** Error: unknown option: %s\n",argv[i]);
      usage(argv);
      exit(1);
    }
  }
  if(make_volrender_script==1){

    NewMemory((void **)&volrender_scriptname,(unsigned int)(len_casename+14+1));
    STRCPY(volrender_scriptname,fdsprefix);
    STRCAT(volrender_scriptname,"_volrender.ssf");

    init_volrender_script(fdsprefix, NULL, startframe0, skipframe0);
  }
#ifndef pp_BETA
  if(time_flag==1){
    STRCAT(TITLE," - ");
    STRCAT(TITLE,__TIME__);
  }
#endif
}

/* ------------------ version ------------------------ */

void display_version_info(void){
    char version[256];
    char revision[256];

    getPROGversion(version);
    getRevision(revision);    // get svn revision number
    PRINTF("\n");
    PRINTF("%s\n\n",TITLERELEASE);
    PRINTF("Version: %s\n",version);
#ifdef BIT64
    PRINTF("Smokeview (64 bit) Revision Number: %s\n",revision);
#else
    PRINTF("Smokeview (32 bit) Revision Number: %s\n",revision);
#endif
#ifdef WIN32
#ifdef X64
    PRINTF("Platform: WIN64 ");
#else
    PRINTF("Platform: WIN32 ");
#endif
#ifdef pp_INTEL
    PRINTF(" (Intel C/C++)\n");
#else
#ifdef WIN32
    PRINTF(" (MSVS C/C++)\n");
#endif
#endif
#endif
#ifndef pp_OSX64
#ifdef pp_OSX
    PRINTF("Platform: OSX\n");
#endif
#endif
#ifdef pp_OSX64
    PRINTF("Platform: OSX64\n");
#endif
#ifndef pp_LINUX64
#ifdef pp_LINUX
    PRINTF("Platform: LINUX\n");
#endif
#endif
#ifdef pp_LINUX64
    PRINTF("Platform: LINUX64\n");
#endif
    PRINTF("Build Date: %s\n",__DATE__);
    if(revision_fds>0){
      PRINTF("FDS Revision Number: %i\n",revision_fds);
    }
    if(smokeviewpath!=NULL){
      PRINTF("Smokeview path: %s\n",smokeviewpath);
    }
    if(smokezippath!=NULL){
      PRINTF("Smokezip path: %s\n",smokezippath);
    }
    if(texturedir!=NULL){
      PRINTF("Texture directory path: %s\n",texturedir);
    }
}

/* ------------------ usage ------------------------ */

void usage(char **argv){
  char buffer[1000];

  PRINTF("%s\n",TITLERELEASE);
  PRINTF("%s\n\n",_("Visualize fire/smoke flow simulations."));
  PRINTF("Usage: %s [options] casename",get_basefilename(buffer,argv[0]));
  PRINTF("%s\n\n",_("where "));
  PRINTF("%s\n",_(" casename       - project id (file names without the extension)"));
  PRINTF("%s\n",_(" -bindir dir    - specify location of smokeview bin directory"));
  PRINTF("%s\n",_(" -build         - show directives used in this build of Smokeview"));
  PRINTF("%s\n",_(" -convert_ini case1.ini case2.ini - update case1.ini to the current format"));
  PRINTF("%s\n",_("                  and save results into case2.ini"));
  PRINTF("%s\n",_(" -demo          - use demonstrator mode of Smokeview"));
  PRINTF("%s\n",_(" -fed            - pre-calculate all FED slice files"));
  PRINTF("%s\n",_(" -help          - display this message"));
  PRINTF("%s\n",_(" -ini           - output default smokeview parameters to smokeview.ini"));
  PRINTF("%s\n",_(" -ng_ini        - No graphics version of -ini."));
  PRINTF("%s\n",_(" -runscript     - run the script file casename.ssf"));
  PRINTF("%s\n",_(" -script scriptfile - run the script file scriptfile"));
  PRINTF("%s\n",_(" -skipframe n   - render every n frames"));
  PRINTF("%s\n",_(" -startframe n  - start rendering at frame n"));
  PRINTF("%s\n",_(" -stereo        - activate stereo mode"));
  PRINTF("%s\n",_(" -tempdir       - forces output files to be written to the temporary directory"));
  PRINTF("%s\n",_(" -update_bounds - calculate boundary file bounds and save to casename.bini"));
  PRINTF("%s\n",_(" -update_ini case.ini - update case.ini to the current format"));
  PRINTF("%s\n",_(" -version       - display version information"));
  PRINTF("%s\n",_(" -volrender     - generate images of volume rendered smoke and fire"));

  if(showbuild==1){
    char label[1024],*labelptr;

    labelptr=label+2;
    strcpy(label,"");
#ifdef _DEBUG
    strcat(label,", _DEBUG");
#endif
#ifdef pp_BETA
    strcat(label,", pp_BETA");
#endif
#ifdef pp_COMPRESS
    strcat(label,", pp_COMPRESS");
#endif
#ifdef pp_CULL
    strcat(label,", pp_CULL");
#endif
#ifdef pp_DEG
    strcat(label,", pp_DEG");
#endif
#ifdef pp_GEOMTEST
    strcat(label,", pp_GEOMTEST");
#endif
#ifdef pp_GPU
    strcat(label,", pp_GPU");
#endif
#ifdef pp_GPUDEPTH
    strcat(label,", pp_GPUDEPTH");
#endif
#ifdef pp_GPUTHROTTLE
    strcat(label,", pp_GPUTHROTTLE");
#endif
#ifdef pp_INTEL
    strcat(label,", pp_INTEL");
#endif
#ifdef pp_LANG
    strcat(label,", pp_LANG");
#endif
#ifdef pp_LINUX
    strcat(label,", pp_LINUX");
#endif
#ifdef pp_LINUX64
    strcat(label,", pp_LINUX64");
#endif
#ifdef pp_MEMDEBUG
    strcat(label,", pp_MEMDEBUG");
#endif
#ifdef pp_memstatus
    strcat(label,", pp_memstatus");
#endif
#ifdef pp_noappend
    strcat(label,", pp_noappend");
#endif
#ifdef pp_OFFICIAL_RELEASE
    strcat(label,", pp_OFFICIAL_RELEASE");
#endif
#ifdef pp_OSX
    strcat(label,", pp_OSX");
#endif
#ifdef pp_OSX64
    strcat(label,", pp_OSX64");
#endif
#ifdef pp_PILOT
    strcat(label,", pp_PILOT");
#endif
#ifdef pp_release
    strcat(label,", pp_release");
#endif
#ifdef pp_THREAD
    strcat(label,", pp_THREAD");
#endif
#ifdef WIN32
    strcat(label,", WIN32");
#endif
#ifdef X64
    strcat(label,", X64");
#endif
    PRINTF("  \n");
    PRINTF("%s\n\n",_("  Smokeview was built using the following pre-processing directives:"));
    PRINTF("%s \n",labelptr);
  }
}
