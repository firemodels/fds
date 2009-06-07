// $Date$ 
// $Revision$
// $Author$

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
#define CPP
#include "ASSERT.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

// svn revision character string
extern "C" char glui_message_revision[]="$Revision$";

#define message_OK 1
#define warning_OK 2
#define error_OK 3
#define abort_OK 4

void Message_CB(int var);
GLUI *glui_message=NULL,*glui_warning=NULL,*glui_error=NULL,*glui_abort=NULL;
GLUI_StaticText *message_text=NULL;
GLUI_StaticText *warning_text=NULL;
GLUI_StaticText *error_text=NULL;
GLUI_StaticText *abort_text=NULL;


/* ------------------ glui_message_setup ------------------------ */

extern "C" void glui_message_setup(int main_window){

  if(glui_message!=NULL)glui_message->close();
  glui_message = GLUI_Master.create_glui("Message:",0,0,0);
  glui_message->add_statictext("");
  message_text=glui_message->add_statictext("");
  glui_message->add_statictext("");
  glui_message->add_separator();
  glui_message->add_button("OK",message_OK,Message_CB);
  glui_message->set_main_gfx_window( main_window );
  glui_message->hide();

  if(glui_warning!=NULL)glui_warning->close();
  glui_warning = GLUI_Master.create_glui("Warning:",0,0,0);
  glui_warning->add_statictext("");
  warning_text=glui_warning->add_statictext("");
  glui_warning->add_statictext("");
  glui_warning->add_separator();
  glui_warning->add_button("OK",warning_OK,Message_CB);
  glui_warning->set_main_gfx_window( main_window );
  glui_warning->hide();

  if(glui_error!=NULL)glui_error->close();
  glui_error = GLUI_Master.create_glui("Error",0,0,0);
  glui_error->add_statictext("");
  error_text=glui_error->add_statictext("");
  glui_error->add_statictext("");
  glui_error->add_separator();
  glui_error->add_button("OK",error_OK,Message_CB);
  glui_error->set_main_gfx_window( main_window );
  glui_error->hide();

  if(glui_abort!=NULL)glui_abort->close();
  glui_abort = GLUI_Master.create_glui("Abort",0,0,0);
  glui_abort->add_statictext("");
  abort_text=glui_abort->add_statictext("");
  glui_abort->add_statictext("");
  glui_abort->add_separator();
  glui_abort->add_button("OK",abort_OK,Message_CB);
  glui_abort->set_main_gfx_window( main_window );
  glui_abort->hide();

}

/* ------------------ message_message ------------------------ */

extern "C" void message_message(char *message){
  printf("%s\n",message);
  if(message_text!=NULL)message_text->set_name(message);
  if(glui_message!=NULL)glui_message->show();
}


/* ------------------ warning_message ------------------------ */

extern "C" void warning_message(char *message){
  printf("***Warning: %s\n",message);
  if(warning_text!=NULL)warning_text->set_name(message);
  if(glui_warning!=NULL)glui_warning->show();
}


/* ------------------ message_message ------------------------ */

extern "C" void error_message(char *message){
  printf("***Error: %s\n",message);
  if(error_text!=NULL)error_text->set_name(message);
  if(glui_error!=NULL)glui_error->show();
}


/* ------------------ message_message ------------------------ */

extern "C" void abort_message(char *message){
  printf("***Fatal Error: %s\n",message);
  if(abort_text!=NULL)abort_text->set_name(message);
  if(glui_abort!=NULL)glui_abort->show();
}

/* ------------------ MESSAGE_CB ------------------------ */

void Message_CB(int var){
  int i;

  switch (var){
    case message_OK:
      glui_message->hide();
      break;
    case warning_OK:
      glui_warning->hide();
      break;
    case error_OK:
      glui_error->hide();
      break;
    case abort_OK:
      glui_abort->hide();
      break;
  }
}


#endif
