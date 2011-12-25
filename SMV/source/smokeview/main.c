// $Date$ 
// $Revision$
// $Author$

#define INMAIN
#include "options.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "MALLOC.h"
#include "smokeviewvars.h"

// svn revision character string
char main_revision[]="$Revision$";
//  dummy change to update version to  5.6.3
//  dummy change to force revision update

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char **argv_sv;
  int startup_flag;
  size_t len;
  char *progname;

  //listdir(".");
  initMALLOC();
  initcolors();
  initvars();
  if(argc==1){
    display_version_info();
  }
#ifdef WIN32
  copy_args(&argc, argv, &argv_sv);
  if(argc==1)return 0;
#else
  argv_sv=argv;
#endif

  if(argc==0)return 0;

  progname=argv_sv[0];
  smvprogdir=getprogdir(progname);
  init_texturedir();
  smokezippath=get_smokezippath(smvprogdir);

  CheckMemory;
  parse_commandline(argc, argv_sv);
  display_version_info();
  if(smokezippath!=NULL)printf("Smokezip file: %s found\n",smokezippath);
  setup_glut(argc,argv_sv);
  startup_flag=setup_case(argc,argv_sv);
  if(startup_flag==0&&update_bounds==1)startup_flag=Update_Bounds();
  if(startup_flag!=0)return 1;

  glutMainLoop();
  return 0;
}	 
