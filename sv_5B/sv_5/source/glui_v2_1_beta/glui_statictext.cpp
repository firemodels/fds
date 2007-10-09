/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_statictext.cpp - GLUI_StaticText Control


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/

#include "glui.h"
#include "stdinc.h"

/****************************** GLUI_StaticText::draw() **********/

void    GLUI_StaticText::draw( int x, int y )
{
  int orig;

  if ( NOT can_draw() )
    return;

  orig = set_to_glut_window();

  draw_text();

  restore_window( orig );
}


/****************************** GLUI_StaticText::set_text() **********/

void    GLUI_StaticText::set_text( char *text )
{
  int orig;

  /**** Erase old text first *****/
  glMatrixMode( GL_MODELVIEW );
  glPushMatrix();
  translate_to_origin();
  erase_text();
  glPopMatrix();

  set_name( text );

  if ( NOT can_draw() )
    return;

  orig = set_to_glut_window();
  /**** Redraw the text in the window ****/
  glMatrixMode( GL_MODELVIEW );
  glPushMatrix();
  translate_to_origin();
  draw_text();
  glPopMatrix();

  restore_window( orig );
}


/************************************ GLUI_StaticText::update_size() **********/

void   GLUI_StaticText::update_size( void )
{
  int text_size;

  if ( NOT glui )
    return;

  text_size = string_width( name );

  if ( w < text_size )
    w = text_size;    
}


/****************************** GLUI_StaticText::draw_text() **********/

void    GLUI_StaticText::draw_text( void )
{
  if ( NOT can_draw() )
    return;

  erase_text();
  draw_name( 0, 9 );
}


/****************************** GLUI_StaticText::erase_text() **********/

void    GLUI_StaticText::erase_text( void )
{
  if ( NOT can_draw() )
    return;

  set_to_bkgd_color();
  glDisable( GL_CULL_FACE );
  glBegin( GL_TRIANGLES );
  glVertex2i( 0,0 );   glVertex2i( w, 0 );  glVertex2i( w, h );  
  glVertex2i( 0, 0 );  glVertex2i( w, h );  glVertex2i( 0, h );   
  glEnd();
}



