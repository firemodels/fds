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
#include GLUT_H

#include "string_util.h"
#include "smokeviewvars.h"

//  dummy change to update version to 6.1.14
// dummy change  to force revision update

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  char **argv_sv;
  int return_code;
  char *progname;

  set_stdout(stdout);
  initMALLOC();
  init_rand_ab(1000000);
  initvars();
  if(argc==1)display_version_info();
  copy_args(&argc, argv, &argv_sv);
  if(argc==0||argc==1)return 0;

  progname=argv_sv[0];
  parse_commandline(argc, argv_sv);
  if(smokeview_bindir==NULL){
    smokeview_bindir=getprogdir(progname,&smokeviewpath);
  }
  init_texturedir();
  smokezippath=get_smokezippath(smokeview_bindir);
  display_version_info();
  setup_glut(argc,argv_sv);
  return_code=setup_case(argc,argv_sv);
  if(return_code==0&&update_bounds==1)return_code=Update_Bounds();
  if(return_code!=0)return 1;
  if(convert_ini==1){
    readini(ini_from);
  }

  glutMainLoop();
  return 0;
}	 
