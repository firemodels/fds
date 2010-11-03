// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#define INMAIN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "svdiff.h"
#include "MALLOC.h"

//dummy change to bump version number to 1.0.3

// svn revision character string
char main_revision[]="$Revision$";

/* ------------------ main ------------------------ */

int main(int argc, char **argv){

  char *smv1=NULL, *smv2=NULL, *arg;
  char smv1_out[1024];
  char smoke1[1024], smoke2[1024], smv_out[1024];
  FILE *stream_out, *stream_in1, *stream_in2;
  int no_plot3d=0, no_slice=0, no_boundary=0;
  int i;
  int open_smokeview=0;

  initMM();
#ifdef WIN32
  strcpy(dirseparator,"\\");
#else
  strcpy(dirseparator,"/");
#endif
  strcpy(pp,"%");

 // check_histogram();  
  test_mode=0;
  sourcedir1=NULL;
  sourcedir2=NULL;
  destdir=NULL;
  strcpy(type_label,"");

  if(argc==1){
    version();
    return 1;
  }

  for(i=1;i<argc;i++){
    int lenarg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
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
        usage();
        return 1;
        break;
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
          usage();
          return 1;
        }
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
        version();
        return 1;
        break;
      default:
        usage();
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
    fullfile(smoke1,sourcedir1,smv1);
    strcat(smoke1,".smv");
  }
  if(smv2!=NULL){
    fullfile(smoke2,sourcedir2,smv2);
    strcat(smoke2,".smv");
  }
  // make sure smv file names exists

  if(getfileinfo(smoke1,NULL,NULL)!=0||getfileinfo(smoke2,NULL,NULL)!=0){
    if(getfileinfo(smoke1,NULL,NULL)!=0){
      printf("***error The .smv file, %s, does not exist\n",smoke1);
    }
    if(getfileinfo(smoke2,NULL,NULL)!=0){
      printf("***error The .smv file, %s, does not exist\n",smoke2);
    }
    return 1;
  }
  make_outfile(smv_out,destdir,smv1_out,".smv");

  stream_out=fopen(smv_out,"w");
  if(stream_out==NULL){
    printf("***error The .smv file, %s, could not be opened for output.\n",smv_out);
  }
  stream_in1=fopen(smoke1,"r");
  if(stream_in1==NULL){
    printf("***error The .smv file, %s, could not be opened for input\n",smoke1);
  }
  stream_in2=fopen(smoke2,"r");
  if(stream_in2==NULL){
    printf("***error The .smv file, %s, could not be opened for input.\n",smoke2);
  }
  if(stream_out==NULL||stream_in1==NULL||stream_in2==NULL)return 1;

  caseinfo[0].dir=sourcedir1;
  caseinfo[0].endian=0;
  caseinfo[1].dir=sourcedir2;
  caseinfo[1].endian=0;

  printf("reading %s\n",smoke1);
  fflush(stdout);
  readsmv(stream_in1, stream_out, caseinfo);
  fclose(stream_in1);

  printf("reading %s\n",smoke2);
  fflush(stdout);
  readsmv(stream_in2, NULL, caseinfo+1);
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
       
/* ------------------ setdir ------------------------ */

char *setdir(char *argdir){
  int lendir;
  char *dir;

  lendir=strlen(argdir);
  NewMemory((void **)&dir,lendir+2);
  strcpy(dir,argdir);
  if(dir[lendir-1]!=dirseparator[0]){
    strcat(dir,dirseparator);
  }
  return dir;
}

/* ------------------ usage ------------------------ */

void usage(void){
  char pp[2];
  char smv_version[100];
  int svn_num;

  getSMDiffversion(smv_version);  // get Smokeview version (ie 5.x.z)
  svn_num=getmaxrevision();    // get svn revision number

  strcpy(pp,"%");
  printf("\n");
  printf("  smokediff [options] smv_case1 smv_case2\n");
  printf("    version: %s (revision %i) - %s\n\n",smv_version,svn_num,__DATE__);

  printf("  smokediff compares two FDS cases by subtracting data referenced in smv_case2 from\n");
  printf("  corresponding data referenced in smv_case1 (smv_case1 - smv_case2).  Slice, PLOT3d\n");
  printf("  and boundary files are supported.  Differenced results may be viewed by opening\n"); 
  printf("  smv_case1_diff.smv in Smokeview or by using the -smv option when running smokediff.\n\n");

  printf("  Mesh bounds must be identical for corresponsing meshes.  Mesh resolutions must be\n");
  printf("  identical when differencing boundary or PLOT3D files.  The x, y, and/or z mesh\n");
  printf("  resolution in smv_case1 must be an integer multiple of the corresponding x, y, z mesh\n");
  printf("  resolution in smv_case2 when differencing slice files.\n\n");

  printf("  -h  - display this message\n");
  printf("  -v  - display version information\n");
  printf("  -s1 dir1 - directory containing case smv_case1.smv\n");
  printf("  -s2 dir2 - directory containing case smv_case2.smv\n");
  printf("  -d  dir  - directory containing created differenced files\n");
  printf("  -nb      - do not difference boundary files\n");
  printf("  -np      - do not difference Plot3d files\n");
  printf("  -ns      - do not difference slice files\n");
  printf("  -smv     - view case in smokeview when differencing is complete\n");
  printf("  -type label - difference only data of type label (in boundary and slice files)\n");
  printf("  smv_case1,smv_case2 - Two smokeview cases to compare.\n");
}


