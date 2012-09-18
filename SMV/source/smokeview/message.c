// $Date$ 
// $Revision$
// $Author$

#include "options.h"

// svn revision character string
char message_revision[]="$Revision$";


#include <stdio.h>
#include <string.h>
/*
#ifdef pp_OSX
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <math.h>
*/
#include "string_util.h"
#include "smokeviewvars.h"


#define LINE_LENGTH 60

/* ------------------ wrap_message ------------------------ */

void wrap_message(char *message){
  char buffer[LINE_LENGTH+1];

  for(;;){
    if(strlen(message)<LINE_LENGTH){
      strcpy(buffer,message);
      fprintf(stderr,"%s\n",buffer);
      break;
    }
    else{
      strncpy(buffer,message,LINE_LENGTH);
      buffer[LINE_LENGTH]=0;
      fprintf(stderr,"%s",buffer);
      message+=LINE_LENGTH;
    }
  }
}

/* ------------------ warning_message ------------------------ */

void warning_message(char *message){
  if(message==NULL||strlen(message)==0)return;
  fprintf(stderr,_("*** Warning: "));
  wrap_message(message);
}


/* ------------------ error_message ------------------------ */

void error_message(char *message){
  if(message==NULL||strlen(message)==0)return;
  fprintf(stderr,_("*** Error: "));
  wrap_message(message);
}

/* ------------------ abort_message ------------------------ */

void abort_message(char *message){
  if(message==NULL||strlen(message)==0)return;
  fprintf(stderr,_("*** Error (fatal): "));
  wrap_message(message);
}


