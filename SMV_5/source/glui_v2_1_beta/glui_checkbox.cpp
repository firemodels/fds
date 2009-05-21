/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_checkbox - GLUI_Checkbox control class


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/

#include "glui.h"
#include "stdinc.h"

/****************************** GLUI_Checkbox::mouse_down_handler() **********/

int    GLUI_Checkbox::mouse_down_handler( int local_x, int local_y )
{
  orig_value = int_val;
  
  TOGGLE_BOOL( int_val );

  currently_inside = true;

  if ( int_val )
    draw_checked();
  else
    draw_unchecked();

  return false;
}


/****************************** GLUI_Checkbox::mouse_up_handler() **********/

int    GLUI_Checkbox::mouse_up_handler( int local_x, int local_y, int inside )
{
  if ( NOT inside ) {
    int_val = orig_value;    
  }
  else {
    set_int_val( int_val );

    /*** Invoke the callback ***/
    execute_callback();
  }

  if ( int_val )
    draw_checked();
  else
    draw_unchecked();

  return false;
}


/****************************** GLUI_Checkbox::mouse_held_down_handler() ******/

int    GLUI_Checkbox::mouse_held_down_handler( int local_x, int local_y,
					       int inside)
{
  /********** Toggle checked and unchecked bitmap if we're entering or
    leaving the checkbox area **********/

  /** oops, this was wrong!
    if ( NOT inside AND currently_inside == true )
    draw_unchecked();
  else if ( inside AND currently_inside == false ) 
    draw_checked();*/

  if ( inside != currently_inside )
    TOGGLE_BOOL( int_val );

  currently_inside = inside;

  if ( int_val )
    draw_checked();
  else
    draw_unchecked();
  
  return false;
}


/****************************** GLUI_Checkbox::key_handler() **********/

int    GLUI_Checkbox::key_handler( unsigned char key,int modifiers )
{
  return false;
}


/****************************** GLUI_Checkbox::draw() **********/

void    GLUI_Checkbox::draw( int x, int y )
{
  int orig;

  if ( NOT can_draw() )
    return;

  orig = set_to_glut_window();

  if ( int_val != 0 ) {
    if ( enabled ) 
      glui->std_bitmaps.draw( GLUI_STDBITMAP_CHECKBOX_ON, 0, 0 );
    else
      glui->std_bitmaps.draw( GLUI_STDBITMAP_CHECKBOX_ON_DIS, 0, 0 );
  }
  else {
    if ( enabled )
      glui->std_bitmaps.draw( GLUI_STDBITMAP_CHECKBOX_OFF, 0, 0 );
    else
      glui->std_bitmaps.draw( GLUI_STDBITMAP_CHECKBOX_OFF_DIS, 0, 0 );      
  }

  draw_active_area();

  draw_name( text_x_offset, 10);

  restore_window(orig);
}


/************************************** GLUI_Checkbox::draw_checked() ******/

void   GLUI_Checkbox::draw_checked( void )
{
  int state, orig;

  if ( NOT can_draw() )
    return;

  orig = set_to_glut_window();
  state = glui->set_front_draw_buffer();
  
  glColor3f( 0.0, 0.0, 0.0 );
  glPushMatrix();
  translate_to_origin();
  if ( enabled )
    glui->std_bitmaps.draw( GLUI_STDBITMAP_CHECKBOX_ON, 0, 0 );
  else
    glui->std_bitmaps.draw( GLUI_STDBITMAP_CHECKBOX_ON_DIS, 0, 0 );
  draw_active_area();
  glPopMatrix();

  glui->restore_draw_buffer(state);
  restore_window(orig);
}


/************************************** GLUI_Checkbox::draw_unchecked() ******/

void   GLUI_Checkbox::draw_unchecked( void )
{
  int state, orig;

  if ( NOT can_draw() )
    return;

  orig = set_to_glut_window();
  state = glui->set_front_draw_buffer();

  glColor3f( 1.0, 1.0, 1.0 );
  glPushMatrix();
  translate_to_origin();
  if ( enabled )
    glui->std_bitmaps.draw( GLUI_STDBITMAP_CHECKBOX_OFF, 0, 0 );
  else
    glui->std_bitmaps.draw( GLUI_STDBITMAP_CHECKBOX_OFF_DIS, 0, 0 );
  draw_active_area();
  glPopMatrix();

  glui->restore_draw_buffer(state);
  restore_window(orig);  
}


/**************************************** GLUI_Checkbox::draw_X() ************/

void   GLUI_Checkbox::draw_X( void )
{
  int orig;

  orig = set_to_glut_window();

  /*  glPointSize( 1.0 );
      glBegin( GL_POINTS );
      for( i=2; i<=GLUI_CHECKBOX_SIZE-2; i++ ) {
      glVertex2i( i,i );
      glVertex2i( i,GLUI_CHECKBOX_SIZE-i );
      }
      glEnd();
      */

  /*  glColor3f( 0.0, 0.0, 0.0 );              */

  /*
    glBegin( GL_LINES );
    glVertex2i( 0,                    0 );
    glVertex2i( GLUI_CHECKBOX_SIZE-0, GLUI_CHECKBOX_SIZE-0 );
    glVertex2i( GLUI_CHECKBOX_SIZE-0, 0  );
    glVertex2i( 0,                    GLUI_CHECKBOX_SIZE-0 );
    glEnd();*/

  restore_window(orig);
}


/********************************** GLUI_Checkbox::draw_empty_box() **********/

void     GLUI_Checkbox::draw_empty_box( void )
{
  int orig;

  if ( NOT can_draw())
    return;
  
  orig = set_to_glut_window();

  glColor3f( 1.0, 1.0, 1.0 );
  glBegin( GL_QUADS );
  glVertex2i( 0, 0 );
  glVertex2i( GLUI_CHECKBOX_SIZE, 0 );
  glVertex2i( GLUI_CHECKBOX_SIZE, GLUI_CHECKBOX_SIZE );
  glVertex2i( 0,                  GLUI_CHECKBOX_SIZE );
  glEnd();

  glColor3f( 0.0, 0.0, 0.0 );
  glBegin( GL_LINE_LOOP );
  glVertex2i( 0, 0 );
  glVertex2i( GLUI_CHECKBOX_SIZE, 0 );
  glVertex2i( GLUI_CHECKBOX_SIZE, GLUI_CHECKBOX_SIZE );
  glVertex2i( 0,                  GLUI_CHECKBOX_SIZE );
  glEnd();

  restore_window(orig);
}


/************************************ GLUI_Checkbox::update_size() **********/

void   GLUI_Checkbox::update_size( void )
{
  int text_size;

  if ( NOT glui )
    return;

  text_size = _glutBitmapWidthString( glui->font, name );

  /*  if ( w < text_x_offset + text_size + 6 )              */
  w = text_x_offset + text_size + 6 ;
}


/**************************** GLUI_Checkbox::draw_active_area() **************/

void    GLUI_Checkbox::draw_active_area( void )
{
  int text_width, left, right, orig;

  if ( NOT can_draw())
    return;

  orig = set_to_glut_window();

  text_width = _glutBitmapWidthString( glui->font, name );
  left       = text_x_offset-3;
  right      = left + 7 + text_width;

  if ( active ) {
    glEnable( GL_LINE_STIPPLE );
    glLineStipple( 1, 0x5555 );
    glColor3f( 0., 0., 0. );
  } else {
    glColor3ub( glui->bkgd_color.r, glui->bkgd_color.g, glui->bkgd_color.b );
  }

  glBegin( GL_LINE_LOOP );
  glVertex2i(left,0);     glVertex2i( right,0);
  glVertex2i(right,h+1);   glVertex2i( left,h+1);
  glEnd();
  
  glDisable( GL_LINE_STIPPLE );

  restore_window(orig);
}


/********************************* GLUI_Checkbox::set_int_val() **************/

void    GLUI_Checkbox::set_int_val( int new_val )
{
  int_val = new_val;

  /*** Update the variable we're (possibly) pointing to ***/
  output_live(true);

  if ( can_draw() ) {
    if ( int_val )
      draw_checked();
    else
      draw_unchecked();
  }
}
