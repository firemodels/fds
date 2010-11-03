// $Date$ 
// $Revision$
// $Author$

#define CPP
#include "options.h"
#ifdef pp_MESSAGE
#include <string.h>
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <math.h>
#include "glui.h"
#include "flowfiles.h"
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

// svn revision character string
extern "C" char glui_message_revision[]="$Revision$";

#define warning_OK 1
#define error_OK 2
#define abort_OK 3

void Message_CB(int var);
GLUI *glui_warning=NULL,*glui_error=NULL,*glui_abort=NULL;
GLUI_StaticText *warning_text=NULL;
GLUI_StaticText *warning_text2=NULL;
GLUI_StaticText *error_text=NULL;
GLUI_StaticText *error_text2=NULL;
GLUI_StaticText *abort_text=NULL;
GLUI_StaticText *abort_text2=NULL;


/* ------------------ glui_message_setup ------------------------ */

extern "C" void glui_message_setup(int main_window){

  if(glui_warning!=NULL)glui_warning->close();
  glui_warning = GLUI_Master.create_glui("Warning:",0,0,0);
  glui_warning->add_statictext("");
  warning_text=glui_warning->add_statictext("");
  warning_text2=glui_warning->add_statictext("");
  glui_warning->add_separator();
  glui_warning->add_button("OK",warning_OK,Message_CB);
  glui_warning->set_main_gfx_window( main_window );
  glui_warning->hide();

  if(glui_error!=NULL)glui_error->close();
  glui_error = GLUI_Master.create_glui("Error",0,0,0);
  glui_error->add_statictext("");
  error_text=glui_error->add_statictext("");
  error_text2=glui_error->add_statictext("");
  glui_error->add_separator();
  glui_error->add_button("OK",error_OK,Message_CB);
  glui_error->set_main_gfx_window( main_window );
  glui_error->hide();

  if(glui_abort!=NULL)glui_abort->close();
  glui_abort = GLUI_Master.create_glui("Abort",0,0,0);
  glui_abort->add_statictext("");
  abort_text=glui_abort->add_statictext("");
  abort_text2=glui_abort->add_statictext("");
  glui_abort->add_separator();
  glui_abort->add_button("OK",abort_OK,Message_CB);
  glui_abort->set_main_gfx_window( main_window );
  glui_abort->hide();

}

/* ------------------ warning_message ------------------------ */

extern "C" void warning_message(char *message){
  char message2[sizeof(GLUI_String)];

  strcpy(message2,"*** Warning: ");
  strcat(message2,message);
  printf("%s\n",message2);
  if(glui_warning==NULL)return;

  if(show_glui_warning==1){
    warning_text2->set_name("*** Additional warnings have occurred.  See command shell for details.");
  }
  else{
    warning_text->set_name(message2);
    warning_text2->set_name("");
  }
  show_glui_warning=1;
  glui_warning->show();
}


/* ------------------ error_message ------------------------ */

extern "C" void error_message(char *message){
  char message2[sizeof(GLUI_String)];

  strcpy(message2,"*** Error: ");
  strcat(message2,message);
  printf("%s\n",message2);
  if(glui_error==NULL)return;

  if(show_glui_error==1){
    error_text2->set_name("*** Additional errors have occurred.  See command shell for details.");
  }
  else{
    error_text->set_name(message2);
    error_text2->set_name("");
  }
  show_glui_error=1;
  glui_error->show();
}


/* ------------------ abort_message ------------------------ */

extern "C" void abort_message(char *message){
  char message2[sizeof(GLUI_String)];

  strcpy(message2,"*** Fatal error: ");
  strcat(message2,message);
  printf("%s\n",message2);
  if(glui_abort==NULL)return;

  if(show_glui_abort==0){
    abort_text->set_name(message2);
    abort_text2->set_name("Press OK to terminate Smokeview");
  }
  show_glui_abort=1;
  glui_abort->show();
}

/* ------------------ MESSAGE_CB ------------------------ */

void Message_CB(int var){
  int i;

  switch (var){
    case warning_OK:
      glui_warning->hide();
      break;
    case error_OK:
      glui_error->hide();
      show_glui_error=0;
      break;
    case abort_OK:
      glui_abort->hide();
      exit(1);
      break;
  }
}


#endif
