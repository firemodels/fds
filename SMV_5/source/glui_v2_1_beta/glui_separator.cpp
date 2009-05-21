/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_separator.cpp - GLUI_Separator control class


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/

#include "glui.h"
#include "stdinc.h"

/****************************** GLUI_Separator::draw() **********/

void    GLUI_Separator::draw( int x, int y )
{
  int width, indent, orig;
  int           cont_x, cont_y, cont_w, cont_h, cont_x_off, cont_y_off;

  if ( NOT can_draw() )
    return;

  orig = set_to_glut_window();

  if ( parent() != NULL ) {
    get_this_column_dims(&cont_x, &cont_y, &cont_w, &cont_h, 
			 &cont_x_off, &cont_y_off);

    width = cont_w - cont_x_off*2;
  }
  else {
    width = this->w;
  }

  indent = width * .05;

  glLineWidth( 1.0 );
  glBegin( GL_LINES );
  glColor3f( .5, .5, .5 );
  glVertex2i( indent,       GLUI_SEPARATOR_HEIGHT/2-1 );    
  glVertex2i( width-indent, GLUI_SEPARATOR_HEIGHT/2-1 );    

  glColor3f( 1., 1., 1. );
  glVertex2i( indent,       GLUI_SEPARATOR_HEIGHT/2 );    
  glVertex2i( width-indent, GLUI_SEPARATOR_HEIGHT/2 );    
  glEnd();

  restore_window(orig);
}


