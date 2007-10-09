/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_button.cpp - GLUI_Button control class


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/

#include "glui.h"
#include "stdinc.h"

/****************************** GLUI_Button::mouse_down_handler() **********/

int    GLUI_Button::mouse_down_handler( int local_x, int local_y )
{
  int_val = 1;   /** A button always in unpressed before here, so
		   now we invariably set it to 'depressed' **/

  currently_inside = true;

  draw_pressed();

  return false;
}


/****************************** GLUI_Button::mouse_up_handler() **********/

int    GLUI_Button::mouse_up_handler( int local_x, int local_y, int inside )
{
  set_int_val( 0 );   /** A button always turns off after you press it **/

  draw_unpressed();

  if ( NOT inside ) {    
  }
  else {
    /*** Invoke the callback ***/
    execute_callback();

    /** Tell the main gfx window to update itself **/
    if( glui )
      glui->post_update_main_gfx();
  }

  return false;
}


/****************************** GLUI_Button::mouse_held_down_handler() ******/

int    GLUI_Button::mouse_held_down_handler( int local_x, int local_y,
					     int new_inside)
{
  if ( NOT new_inside AND currently_inside == true ) {
    draw_unpressed();
  } 
  else if ( new_inside AND currently_inside == false ) {
    draw_pressed();
  }

  currently_inside = new_inside;
  
  return false;
}


/****************************** GLUI_Button::key_handler() **********/

int    GLUI_Button::key_handler( unsigned char key,int modifiers )
{
  return false;
}


/************************************** GLUI_Button::draw_pressed() ******/

void   GLUI_Button::draw_pressed( void )
{
  int state, orig;

  if ( NOT can_draw() )
    return;

  orig  = set_to_glut_window();
  state = glui->set_front_draw_buffer();
  
  glColor3f( 0.0, 0.0, 0.0 );
  glPushMatrix();
  translate_to_origin();

  draw_text( 1 );

  glBegin( GL_LINE_LOOP );
  glVertex2i( 0, 0 );         glVertex2i( w, 0 );
  glVertex2i( w, h );         glVertex2i( 0, h );
  glEnd();

  glBegin( GL_LINE_LOOP );
  glVertex2i( 1, 1 );         glVertex2i( w-1, 1 );
  glVertex2i( w-1, h-1 );     glVertex2i( 1, h-1 );
  glEnd();

  glPopMatrix();

  glui->restore_draw_buffer(state);
  restore_window(orig);
}


/************************************** GLUI_Button::draw_unpressed() ******/

void   GLUI_Button::draw_unpressed( void )
{
  if ( NOT can_draw() )
    return;

  translate_and_draw_front();
}




/********************************************** GLUI_Button::draw() **********/

void    GLUI_Button::draw( int x, int y )
{
  if ( NOT can_draw() )
    return;

  if ( glui )
    glui->draw_raised_box( 0, 0, w, h );

  draw_text( 0 );
}


/**************************************** GLUI_Button::draw_text() **********/

void     GLUI_Button::draw_text( int sunken )
{
  int string_width, orig;

  if ( NOT can_draw() )
    return;

  orig = set_to_glut_window();

  glColor3ub( glui->bkgd_color.r, glui->bkgd_color.g, glui->bkgd_color.b );
  glDisable( GL_CULL_FACE );
  glBegin( GL_QUADS );
  glVertex2i( 2, 2 );         glVertex2i( w-2, 2 );
  glVertex2i( w-2, h-2 );     glVertex2i( 2, h-2 );
  glEnd();

  glColor3ub( 0,0,0 );
  
  string_width = _glutBitmapWidthString( glui->font,
					 this->name );
  if ( NOT sunken ) {
    draw_name( MAX((w-string_width),0)/2, 13);
  }
  else {
    draw_name( MAX((w-string_width),0)/2 + 1, 13 + 1);
  }

  if ( active ) {
    glEnable( GL_LINE_STIPPLE );
    glLineStipple( 1, 0x5555 );
    
    glColor3f( 0., 0., 0. );
    
    glBegin( GL_LINE_LOOP );
    glVertex2i( 3, 3 );         glVertex2i( w-3, 3 );
    glVertex2i( w-3, h-3 );     glVertex2i( 3, h-3 );
    glEnd();
    
    glDisable( GL_LINE_STIPPLE );
  }

  restore_window(orig);
}


/************************************** GLUI_Button::update_size() **********/

void   GLUI_Button::update_size( void )
{
  int text_size;

  if ( NOT glui )
    return;

  text_size = string_width( name );

  if ( w < text_size + 16 )
    w = text_size + 16 ;
}
