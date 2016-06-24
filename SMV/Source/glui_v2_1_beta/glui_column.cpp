/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_column.cpp - GLUI_Column control class


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/

#include "glui.h"
#include "stdinc.h"

/**************************************** GLUI_Column::draw() ************/

void  GLUI_Column::draw( int x, int y )
{
  int   orig;
  int   panel_x, panel_y, panel_w, panel_h, panel_x_off, panel_y_off;
  int   y_diff;
  
  if ( NOT can_draw() )
    return;

  if ( int_val == 1 ) {  /* Draw a vertical bar */
    orig = set_to_glut_window();

    if ( parent() != NULL ) {
      get_this_column_dims(&panel_x, &panel_y, &panel_w, &panel_h, 
			   &panel_x_off, &panel_y_off);

      y_diff = y_abs - panel_y;

      if ( 0 ) {
	glLineWidth(1.0);
	glBegin( GL_LINES );
	glColor3f( .5, .5, .5 );
	glVertex2i( -GLUI_XOFF+1, -y_diff + GLUI_SEPARATOR_HEIGHT/2 );
	glVertex2i( -GLUI_XOFF+1, -y_diff + panel_h - GLUI_SEPARATOR_HEIGHT/2);

	glColor3f( 1.0, 1.0, 1.0 );
	glVertex2i( -GLUI_XOFF+2, -y_diff + GLUI_SEPARATOR_HEIGHT/2 );
	glVertex2i( -GLUI_XOFF+2, -y_diff + panel_h - GLUI_SEPARATOR_HEIGHT/2);
	glEnd();
      }
      else {
	glLineWidth(1.0);
	glBegin( GL_LINES );
	glColor3f( .5, .5, .5 );
	glVertex2i( -2, 0 );
	glVertex2i( -2, h );
	/*glVertex2i( 0, -y_diff + GLUI_SEPARATOR_HEIGHT/2 );              */
	/*glVertex2i( 0, -y_diff + panel_h - GLUI_SEPARATOR_HEIGHT/2);              */

	glColor3f( 1.0, 1.0, 1.0 );
	glVertex2i( -1, 0 );
	glVertex2i( -1, h );
	/*glVertex2i( 1, -y_diff + GLUI_SEPARATOR_HEIGHT/2 );              */
	/*glVertex2i( 1, -y_diff + panel_h - GLUI_SEPARATOR_HEIGHT/2);              */
	glEnd();
      }
			
    }    

    restore_window(orig);
  }
}

