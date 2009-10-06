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

//dummy change to bump version number to 1.0.0

// svn revision character string
char main_revision[]="$Revision$";


/* ------------------ main ------------------------ */

int main(int argc, char **argv){

  char *smv1=NULL, *smv2=NULL, *prog=NULL, *arg;
  char smoke1[1024], smoke2[1024], smv_out[1024];
  FILE *stream_out, *stream_in1, *stream_in2;
  int i;

#ifdef WIN32
  strcpy(dirseparator,"\\");
#else
  strcpy(dirseparator,"/");
#endif

  sourcedir1=NULL;
  sourcedir2=NULL;
  destdir=NULL;

  prog=argv[0];
  if(argc==1){
    version();
    return 1;
  }

  for(i=1;i<argc;i++){
    int lenarg;

    arg=argv[i];
    lenarg=strlen(arg);
    if(arg[0]=='-'&&lenarg>1){
      switch(arg[1]){
      case 'h':
        usage(prog);
        return 1;
        break;
      case 's':
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
      default:
        usage(prog);
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

  if(smv1!=NULL){
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
  make_outfile(smv_out,destdir,smoke1,".smv");

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

  readsmv(stream_in1, stream_out, caseinfo);
  fclose(stream_in1);
  readsmv(stream_in2, NULL, caseinfo+1);
  fclose(stream_in2);
  setup_slice(stream_out);
  diff_slices();
  fclose(stream_out);
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

void usage(char *prog){
  char pp[2];
  char smv_version[100];
  int svn_num;

  getSMDiffversion(smv_version);  // get Smokeview version (ie 5.x.z)
  svn_num=getmaxrevision();    // get svn revision number

  strcpy(pp,"%");
  printf("\n");
  printf("  smokediff %s (%i) - %s\n\n",smv_version,svn_num,__DATE__);
  printf("  smokediff compares two FDS cases by subtracting corresponding slice\n");
  printf("  files referenced in two .smv files specified on the command line.\n");
  printf("  PLOT3D files are also differenced.  A new .smv file is generated.\n");
  printf("  Mesh dimensions in the two cases must be identical, though the number\n");
  printf("  of grid cells need not be.  Usage:\n\n");
  printf("  %s [-h] smv_case1 smv_case2\n",prog);
  printf("  -h  - display this message\n\n");
  printf("  smv_case1,smv_case2 - Two smokeview cases to compare.  smokediff\n");
  printf("        generates an .smv file named smv_case1_diff.smv\n\n");
}
