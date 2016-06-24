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

/****************************** GLUI_Rollout::open() **********/

void    GLUI_Rollout::open( void )
{
  int orig;

  if ( NOT glui )
    return;

  if ( is_open )
    return;

  is_open = true;

  orig = set_to_glut_window();

  child_head = collapsed_node.child_head;
  child_tail = collapsed_node.child_tail;

  collapsed_node.child_head = NULL;
  collapsed_node.child_tail = NULL;

  if ( child_head != NULL ) {
    ((GLUI_Control*) child_head)->unhide_internal( true );
  }

  glui->refresh();

#ifndef pp_GLUI_ORIG  
  execute_callback();
#endif

  restore_window(orig);
}


/****************************** GLUI_Rollout::close() **********/

void    GLUI_Rollout::close( void )
{
  int orig;

  if ( NOT glui )
    return;

  if ( NOT is_open )
    return;

  orig = set_to_glut_window();

  if ( child_head != NULL ) {
    ((GLUI_Control*) child_head)->hide_internal( true );
  }

  collapsed_node.child_head = first_child();
  collapsed_node.child_tail = last_child();

  child_head = NULL;
  child_tail = NULL;

  restore_window(orig);

  this->h = GLUI_DEFAULT_CONTROL_HEIGHT + 7;

  is_open = false;

  glui->refresh();
}


/**************************** GLUI_Rollout::mouse_down_handler() **********/


int   GLUI_Rollout::mouse_down_handler( int local_x, int local_y )
{
  if ( local_y - y_abs > 18 ) {
    initially_inside = currently_inside = false;
    return false;
  }

  currently_inside = true;
  initially_inside = true;

  draw_pressed();

  return false;
}


/**************************** GLUI_Rollout::mouse_down_handler() **********/

int   GLUI_Rollout::mouse_up_handler( int local_x, int local_y, int inside )
{
  draw_unpressed();

  if ( currently_inside ) {    
    if ( is_open )
      close();
    else
      open();
  }

  initially_inside = false;

  return false;
}


/********************************* GLUI_Rollout::draw() ***********/

void   GLUI_Rollout::draw( int x, int y )
{
  int orig, left, right, top, bottom;

  if ( NOT can_draw() )
    return;

  orig = set_to_glut_window();
	
  left   = 5;
  right  = w-left;
  top    = 3;
  bottom = 3+16;

  if ( is_open ) 
    draw_emboss_box( 0, w, top+3, h );
  else
    draw_emboss_box( 0, w, top+3, h-7 );

  glui->draw_raised_box( left, top, w-left*2, 16 );

  if ( glui )
    glColor3ub(glui->bkgd_color.r,glui->bkgd_color.g,glui->bkgd_color.b);
  glDisable( GL_CULL_FACE );
  glBegin( GL_QUADS );
  glVertex2i( left+1, top+1 );      glVertex2i( right-1, top+1 );
  glVertex2i( right-1, bottom-1 );  glVertex2i( left+1, bottom-1 );
  glEnd();

  draw_name( left+8, top+11 );

  if ( active ) 
    /*draw_active_box( left+4, left+string_width( name.string )+12,       */
    draw_active_box( left+4, right-17, 
		     top+2, bottom-2 );


  /**   Draw '+' or '-'  **/

  glBegin( GL_LINES );
  if ( is_open ) {
    if ( enabled )		glColor3f( 0.0, 0.0, 0.0 );
    else			glColor3f( 0.5, 0.5, 0.5 );
    glVertex2i(right-14,(top+bottom)/2);  glVertex2i(right-5,(top+bottom)/2);

    glColor3f( 1.0, 1.0, 1.0 );
    glVertex2i(right-14,1+(top+bottom)/2);glVertex2i(right-5,1+(top+bottom)/2);
  }
  else
  {
    glColor3f( 1.0, 1.0, 1.0 );
    glVertex2i(right-9,top+3);							glVertex2i(right-9,bottom-4);
    glVertex2i(right-14,(top+bottom)/2);		glVertex2i(right-5,(top+bottom)/2);

    if ( enabled )		glColor3f( 0.0, 0.0, 0.0 );
    else			glColor3f( 0.5, 0.5, 0.5 );
    glVertex2i(right-14,-1+(top+bottom)/2);
    glVertex2i(right-5,-1+(top+bottom)/2);
    glVertex2i(right-10,top+3);
    glVertex2i(right-10,bottom-4);
  }
  glEnd();

  glLineWidth( 1.0 );

  restore_window(orig);
}


/***************************** GLUI_Rollout::update_size() **********/

void   GLUI_Rollout::update_size( void )
{
  int text_size;

  if ( NOT glui )
    return;

  text_size = string_width(name);

  if ( w < text_size + 36 )
    w = text_size + 36;
}


/**************************** GLUI_Rollout::draw_pressed() ***********/

void   GLUI_Rollout::draw_pressed( void )
{
  int state, orig;
  int left, right, top, bottom;

  left   = 5;
  right  = w-left;
  top    = 3;
  bottom = 3+16;

  if ( NOT can_draw() )
    return;

  orig  = set_to_glut_window();
  state = glui->set_front_draw_buffer();
  
  glColor3f( 0.0, 0.0, 0.0 );
  glPushMatrix();
  translate_to_origin();

  glBegin( GL_LINE_LOOP );
  glVertex2i( left, top );         glVertex2i( right, top );
  glVertex2i( right, bottom );     glVertex2i( left,bottom );
  glEnd();

  glBegin( GL_LINE_LOOP );
  glVertex2i( left+1, top+1 );         glVertex2i( right-1, top+1 );
  glVertex2i( right-1, bottom-1 );     glVertex2i( left+1,bottom-1 );
  glEnd();

  glPopMatrix();

  glui->restore_draw_buffer(state);
  restore_window(orig);

}


/**************************** GLUI_Rollout::draw_unpressed() ***********/

void   GLUI_Rollout::draw_unpressed( void )
{
  if ( NOT can_draw() )
    return;

  translate_and_draw_front();
}


/**************************** GLUI_Rollout::mouse_held_down_handler() ****/

int  GLUI_Rollout::mouse_held_down_handler( 
					   int local_x, int local_y, 
					   int new_inside )
{
  if ( NOT initially_inside )
    return false;

  if ( local_y - y_abs> 18 )
    new_inside = false;

  if ( NOT new_inside AND currently_inside == true ) {
    draw_unpressed();
  } 
  else if ( new_inside AND currently_inside == false ) {
    draw_pressed();
  }

  currently_inside = new_inside;
  
  return false;
}
