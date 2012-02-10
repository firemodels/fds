// $Date$ 
// $Revision$
// $Author$

// svn revision character string
char main_revision[]="$Revision$";

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

#include "string_util.h"
#include "smokeviewvars.h"

//  dummy change to update version to  5.6.3
//  dummy change to force revision update

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char **argv_sv;
  int return_code;
  char *progname;

  initMALLOC();
  initvars();
  if(argc==1)display_version_info();
  copy_args(&argc, argv, &argv_sv);
  if(argc==0||argc==1)return 0;

  progname=argv_sv[0];
  smokeview_bindir=getprogdir(progname,&smokeviewpath);
  init_texturedir();
  smokezippath=get_smokezippath(smokeview_bindir);
  parse_commandline(argc, argv_sv);
  display_version_info();
  setup_glut(argc,argv_sv);
  return_code=setup_case(argc,argv_sv);
  if(return_code==0&&update_bounds==1)return_code=Update_Bounds();
  if(return_code!=0)return 1;

  glutMainLoop();
  return 0;
}	 
