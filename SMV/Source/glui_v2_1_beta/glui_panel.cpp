/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_panel.cpp - GLUI_Panel control class


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/

#include "glui.h"
#include "stdinc.h"

/****************************** GLUI_Panel::draw() **********/

void    GLUI_Panel::draw( int x, int y )
{
  int top, orig;

  if ( NOT can_draw() )
    return;

  orig = set_to_glut_window();

  if ( int_val == GLUI_PANEL_RAISED ) {
    top = 0;
    glLineWidth( 1.0 );
    glColor3f( 1.0, 1.0, 1.0 );
    glBegin( GL_LINE_LOOP );
    glVertex2i( 0, top );    glVertex2i( w, top );
    glVertex2i( 0, top );    glVertex2i( 0, h );
    glEnd();
    
    glColor3f( .5, .5, .5 );
    glBegin( GL_LINE_LOOP );
    glVertex2i( w, top );
    glVertex2i( w, h );
    glVertex2i( 0, h );
    glVertex2i( w, h );
    glEnd();

    /** ORIGINAL RAISED PANEL METHOD - A LITTLE TOO HIGH **
    glLineWidth(1.0);
    glBegin( GL_LINES );
    glColor3f( 1.0, 1.0, 1.0 );
    glVertex2i( 1, 1 );    glVertex2i( w-2, 1 );
    glVertex2i( 1, 1 );    glVertex2i( 1, h-2 );
    
    glColor3f( .5, .5, .5 );
    glVertex2i( w-1, 1 );    glVertex2i( w-1, h-1 );
    glVertex2i( 1, h-1 );    glVertex2i( w-1, h-1 );
    
    glColor3f( 0.0, 0.0, 0.0 );
    glVertex2i( 0, h );    glVertex2i( w, h );
    glVertex2i( w, 0 );    glVertex2i( w, h );
    glEnd();
    
    -- Touch up the lines a bit (needed in some opengl implementations   
    glBegin( GL_POINTS );
    glColor3f( .5, .5, .5 );
    glVertex2i( w-1, h-1 );
    glColor3f( 0.0, 0.0, 0.0 );
    glVertex2i( w, h );
    glEnd();
    **/    
      }
  else if ( int_val == GLUI_PANEL_EMBOSSED ) {
    if ( name[0] == '\0' ) {
      top = 0;
    }
    else {
      top = GLUI_PANEL_EMBOSS_TOP;
    }

    glLineWidth( 1.0 );
    glColor3f( 1.0, 1.0, 1.0 );
    glBegin( GL_LINE_LOOP );
    glVertex2i( 0, top );    glVertex2i( w, top );
    glVertex2i( w, h );    glVertex2i( 0, h );

    glVertex2i( 1, top+1 );    glVertex2i( w-1, top+1 );
    glVertex2i( w-1, h-1 );    glVertex2i( 1, h-1 );
    glEnd();
    
    glColor3f( .5, .5, .5 );
    glBegin( GL_LINE_LOOP );
    glVertex2i( 0, top );
    glVertex2i( w-1, top );
    glVertex2i( w-1, h-1 );
    glVertex2i( 0, h-1 );
    glEnd();

    /**** Only display text in embossed panel ****/
    if ( name[0] != '\0' ) { /* Only  draw non-null strings */
      int left = 7, height=GLUI_PANEL_NAME_DROP+1;
      int str_width;

      str_width = string_width(name);

      if ( glui )
	glColor3ub(glui->bkgd_color.r,glui->bkgd_color.g,glui->bkgd_color.b);
      glDisable( GL_CULL_FACE );
      glBegin( GL_QUADS );
      glVertex2i( left-3, 0 );               glVertex2i( left+str_width+3, 0 );
      glVertex2i( left+str_width+3, height );  glVertex2i( left-3, height );
      glEnd();

      draw_name( left, GLUI_PANEL_NAME_DROP );
    }
  }

  glLineWidth( 1.0 );

  restore_window(orig);
}


/****************************** GLUI_Panel::set_name() **********/

void    GLUI_Panel::set_name( char *new_name )
{
  strncpy(name,new_name,sizeof(GLUI_String));

  update_size();

  if ( glui )
    glui->refresh();
}


/****************************** GLUI_Panel::set_type() **********/

void    GLUI_Panel::set_type( int new_type )
{
  int old_window;

  if ( new_type != int_val ) {
    int_val = new_type;

    /*    translate_and_draw_front();             */
    update_size();

    old_window = set_to_glut_window();
    glutPostRedisplay( );
    restore_window( old_window );
  }
}


/************************************** GLUI_Panel::update_size() **********/

void   GLUI_Panel::update_size( void )
{
  int text_size;

  if ( NOT glui )
    return;

  text_size = string_width(name);

  if ( w < text_size + 16 )
    w = text_size + 16 ;

  if ( name[0] != '\0' AND int_val == GLUI_PANEL_EMBOSSED ) {
    this->y_off_top = GLUI_YOFF + 8;
  }
  else {
    this->y_off_top = GLUI_YOFF;
  }
}
