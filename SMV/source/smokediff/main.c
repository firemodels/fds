#include "options.h"
#define INMAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "svdiff.h"
#include "MALLOC.h"

//dummy change to bump version number to 1.0.10
//dummy change to force githash change

/* ------------------ main ------------------------ */

int main(int argc, char **argv){

  char *smv1=NULL, *smv2=NULL, *arg;
  char smv1_out[1024];
  char svdlogfile[1024];
  char *smoke1, *smoke2, smv_out[1024];
  char smoke1a[1024], smoke2a[1024];
  char smoke1b[1024], smoke2b[1024];
  char fed_smoke1[1024], fed_smoke2[1024];

  FILE *stream_out, *stream_in1, *stream_in2;
  int no_plot3d=0, no_slice=0, no_boundary=0;
  int i;
  int open_smokeview=0;
  int redirect=0;

  display_warnings=1;
  set_stdout(stdout);
  initMALLOC();
#ifdef WIN32
  strcpy(dirseparator,"\\");
#else
  strcpy(dirseparator,"/");
#endif
  strcpy(pp,"%");

  NewMemory((void **)&caseinfo,2*sizeof(casedata));


 // check_histogram();
  test_mode=0;
  sourcedir1=NULL;
  sourcedir2=NULL;
  destdir=NULL;
  strcpy(type_label,"");

  if(argc==1){
    PRINTversion("Smokediff ");
    return 1;
  }

/* -e{850} loop index i is modified within loop */
  for(i=1;i<argc;i++){
    arg=argv[i];
    if(arg[0]=='-'&&strlen(arg)>1){
      char *key;

      key = arg+1;
      switch(arg[1]){
      case 't':
        if(strcmp(key,"type")==0){
          char *label;

          i++;
          strcpy(type_label,"");
          if(i<argc){
            label=argv[i];
            if(label!=NULL&&strlen(label)>0){
              strcpy(type_label,label);
            }
          }
        }
        else{
          test_mode=1;
        }
        break;
      case 'h':
        Usage();
        return 1;
      case 'n':
        if(arg[2]=='p'){
          no_plot3d=1;
        }
        else if(arg[2]=='s'){
          no_slice=1;
        }
        else if(arg[2]=='b'){
          no_boundary=1;
        }
        else{
          Usage();
          return 1;
        }
        break;
      case 'r':
        redirect=1;
        break;
      case 's':
        if(arg[2]=='m'&&arg[3]=='v'){
          open_smokeview=1;
          break;
        }
        if(i+1>=argc)break;
        if(arg[2]=='1'){
          sourcedir1=setdir(argv[i+1]);
          if(sourcedir1==NULL)return 1;
          i++;
        }
        if(arg[2]=='2'){
          sourcedir2=setdir(argv[i+1]);
          if(sourcedir2==NULL)return 1;
          i++;
        }
        break;
      case 'd':
        if(i+1>=argc)break;
        destdir=setdir(argv[i+1]);
        if(destdir==NULL)return 1;
        i++;
        break;
      case 'v':
        PRINTversion("Smokediff ");
        return 1;
      case 'w':
        display_warnings=0;
        break;
      default:
        Usage();
        return 1;
      }
    }
    else{
      if(smv1==NULL){
        smv1=argv[i];
      }
      else{
        smv2=argv[i];
      }
    }
  }

  strcpy(smv1_out,"");
  if(smv1!=NULL){
    strcat(smv1_out,smv1);
    strcat(smv1_out,".smv");
    fullfile(smoke1a,sourcedir1,smv1);

    strcpy(fed_smoke1,smoke1a);
    strcat(fed_smoke1,".fed_smv");

    strcpy(smoke1b,smoke1a);
    strcat(smoke1b,".smvtmp");

    strcat(smoke1a,".smv");
    smoke1 = smoke1a;

    if(file_exists(fed_smoke1)==1){
      copyfile(".",smoke1a, smoke1b, REPLACE_FILE);
      copyfile(".",fed_smoke1, smoke1b, APPEND_FILE);
      smoke1=smoke1b;
    }
  }
  if(smv2!=NULL){
    fullfile(smoke2a,sourcedir2,smv2);

    strcpy(fed_smoke2,smoke2a);
    strcat(fed_smoke2,".fed_smv");

    strcpy(smoke2b,smoke2a);
    strcat(smoke2b,".smvtmp");
    strcat(smoke2a,".smv");
    smoke2 = smoke2a;

    if(file_exists(fed_smoke2)==1){
      copyfile(".",smoke2a, smoke2b, REPLACE_FILE);
      copyfile(".",fed_smoke2, smoke2b, APPEND_FILE);
      smoke2=smoke2b;
    }
  }
  // make sure smv file names exists

  if(redirect==1){
    strcpy(svdlogfile,"");
    if(destdir!=NULL)strcat(svdlogfile,destdir);
    strcat(svdlogfile,smv1);
    strcat(svdlogfile,"_diff.svdlog");
    LOG_FILENAME=fopen(svdlogfile,"w");
    if(LOG_FILENAME!=NULL){
      set_stdout(LOG_FILENAME);
    }
  }
  if(getfileinfo(smoke1,NULL,NULL)!=0||getfileinfo(smoke2,NULL,NULL)!=0){
    if(getfileinfo(smoke1,NULL,NULL)!=0){
      fprintf(stderr,"*** Error The .smv file, %s, does not exist\n",smoke1);
    }
    if(getfileinfo(smoke2,NULL,NULL)!=0){
      fprintf(stderr,"*** Error The .smv file, %s, does not exist\n",smoke2);
    }
    return 1;
  }
  make_outfile(smv_out,destdir,smv1_out,".smv");

  stream_out=fopen(smv_out,"w");
  if(stream_out==NULL){
    fprintf(stderr,"*** Error The .smv file, %s, could not be opened for output.\n",smv_out);
  }
  stream_in1=fopen(smoke1,"r");
  if(stream_in1==NULL){
    fprintf(stderr,"*** Error The .smv file, %s, could not be opened for input\n",smoke1);
  }
  stream_in2=fopen(smoke2,"r");
  if(stream_in2==NULL){
    fprintf(stderr,"*** Error The .smv file, %s, could not be opened for input.\n",smoke2);
  }
  if(stream_out==NULL||stream_in1==NULL||stream_in2==NULL)return 1;

  caseinfo[0].dir=sourcedir1;
  caseinfo[0].endian=0;
  caseinfo[1].dir=sourcedir2;
  caseinfo[1].endian=0;

  PRINTF("reading %s\n",smoke1);
  FFLUSH();
  ReadSMV(stream_in1, stream_out, caseinfo);
  fclose(stream_in1);

  PRINTF("reading %s\n",smoke2);
  FFLUSH();
  ReadSMV(stream_in2, NULL, caseinfo+1);
  fclose(stream_in2);

  if(no_plot3d==0){
    setup_plot3d(stream_out);
    diff_plot3ds(stream_out);
  }
  if(no_slice==0){
    setup_slice(stream_out);
    diff_slices(stream_out);
  }
  if(no_boundary==0){
    setup_boundary(stream_out);
    diff_boundaryes(stream_out);
  }

  fclose(stream_out);
  if(open_smokeview==1){
    char command[1024];

    strcpy(command,"smokeview ");
    strcat(command,smv_out);
    system(command);
  }

  return 0;
}

/* ------------------ usage ------------------------ */

void Usage(void){
  char smv_version[100];
  char githash[100];
  char gitdate[100];

  getPROGversion(smv_version);  // get Smokeview version (ie 5.x.z)
  getGitInfo(githash,gitdate);    // get githash

  PRINTF("\n");
  PRINTF("  smokediff [options] smv_case1 smv_case2\n");
  PRINTF("    version: %s (githash %s) - %s\n\n",smv_version,githash,__DATE__);

  PRINTF("  smokediff compares two FDS cases by subtracting data referenced in smv_case2 from\n");
  PRINTF("  corresponding data referenced in smv_case1 (smv_case1 - smv_case2).  Slice, PLOT3d\n");
  PRINTF("  and boundary files are supported.  Differenced results may be viewed by opening\n");
  PRINTF("  smv_case1_diff.smv in Smokeview or by using the -smv option when running smokediff.\n\n");

  PRINTF("  Mesh bounds must be identical for corresponding meshes.  Mesh resolutions must be\n");
  PRINTF("  identical when differencing boundary and PLOT3D files.  The x, y, and z mesh\n");
  PRINTF("  resolutions in smv_case2 must be integer multiples of the corresponding x, y, z mesh\n");
  PRINTF("  resolutions in smv_case1 when differencing slice files.\n\n");

  PRINTF("  -h  - display this message\n");
  PRINTF("  -v  - display version information\n");
  PRINTF("  -s1 dir1 - directory containing case smv_case1.smv\n");
  PRINTF("  -s2 dir2 - directory containing case smv_case2.smv\n");
  PRINTF("  -d  dir  - directory containing created differenced files\n");
  PRINTF("  -nb      - do not difference boundary files\n");
  PRINTF("  -np      - do not difference Plot3d files\n");
  PRINTF("  -ns      - do not difference slice files\n");
  PRINTF("  -smv     - view case in smokeview when differencing is complete\n");
  PRINTF("  -type label - difference only data of type label (in boundary and slice files)\n");
  PRINTF("  smv_case1,smv_case2 - Two smokeview cases to compare.\n");
}


