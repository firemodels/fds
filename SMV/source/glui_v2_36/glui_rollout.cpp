/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_panel.cpp - GLUI_Panel control class


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  WWW:    http://sourceforge.net/projects/glui/
  Forums: http://sourceforge.net/forum/?group_id=92496

  This software is provided 'as-is', without any express or implied 
  warranty. In no event will the authors be held liable for any damages 
  arising from the use of this software. 

  Permission is granted to anyone to use this software for any purpose, 
  including commercial applications, and to alter it and redistribute it 
  freely, subject to the following restrictions: 

  1. The origin of this software must not be misrepresented; you must not 
  claim that you wrote the original software. If you use this software 
  in a product, an acknowledgment in the product documentation would be 
  appreciated but is not required. 
  2. Altered source versions must be plainly marked as such, and must not be 
  misrepresented as being the original software. 
  3. This notice may not be removed or altered from any source distribution. 

*****************************************************************************/

#include "glui_internal_control.h"

enum {rollout_height_pixels=GLUI_DEFAULT_CONTROL_HEIGHT + 7};

/****************************** GLUI_Rollout::GLUI_Rollout() **********/

GLUI_Rollout::GLUI_Rollout( GLUI_Node *parent, const char *name, 
                            int open, int type )
{
  common_init();
  set_name( name );
  user_id    = -1;
  int_val    = type;
		
  if ( NOT open ) {
    is_open = false;
    h = rollout_height_pixels;
  }

  parent->add_control( this );
}

/****************************** GLUI_Rollout::open() **********/

void    GLUI_Rollout::open( void )
{
  if ( NOT glui )
    return;

  if ( is_open )
    return;
  is_open = true;

  GLUI_DRAWINGSENTINAL_IDIOM

  /* Copy hidden children into our private list "collapsed_node" */
  child_head = collapsed_node.child_head;
  child_tail = collapsed_node.child_tail;
  collapsed_node.child_head = NULL;
  collapsed_node.child_tail = NULL;

  if ( child_head != NULL ) {
    ((GLUI_Control*) child_head)->unhide_internal( true );
  }

  glui->refresh();
}


/****************************** GLUI_Rollout::close() **********/

void    GLUI_Rollout::close( void )
{
  if ( NOT glui )
    return;

  if ( NOT is_open )
    return;
  is_open = false;
  
  GLUI_DRAWINGSENTINAL_IDIOM

  if ( child_head != NULL ) {
    ((GLUI_Control*) child_head)->hide_internal( true );
  }

  /* Move all children into a private list of hidden children */
  collapsed_node.child_head = first_child();
  collapsed_node.child_tail = last_child();
  child_head = NULL;
  child_tail = NULL;

  this->h = rollout_height_pixels;

  glui->refresh();
}


/**************************** GLUI_Rollout::mouse_down_handler() **********/


int   GLUI_Rollout::mouse_down_handler( int local_x, int local_y )
{
  if ( local_y - y_abs > rollout_height_pixels ) {
    initially_inside = currently_inside = false;
    return false;
  }

  currently_inside = true;
  initially_inside = true;
  redraw();

  return false;
}


/**************************** GLUI_Rollout::mouse_held_down_handler() ****/

int  GLUI_Rollout::mouse_held_down_handler( 
					   int local_x, int local_y, 
					   bool new_inside )
{
  if ( NOT initially_inside )
    return false;

  if ( local_y - y_abs> rollout_height_pixels )
    new_inside = false;
  
  if (new_inside != currently_inside) {
     currently_inside = new_inside;
     redraw();
  }
  
  return false;
}


/**************************** GLUI_Rollout::mouse_down_handler() **********/

int   GLUI_Rollout::mouse_up_handler( int local_x, int local_y, bool inside )
{
  if ( currently_inside ) {    
    if ( is_open )
      close();
    else
      open();
  }

  currently_inside = false;
  initially_inside = false;
  redraw();

  return false;
}


/********************************* GLUI_Rollout::draw() ***********/

void   GLUI_Rollout::draw( int x, int y )
{
  GLUI_DRAWINGSENTINAL_IDIOM
  
  int left, right, top, bottom;

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
    glColor3ubv(glui->bkgd_color);
  glDisable( GL_CULL_FACE );
  glBegin( GL_QUADS );
  glVertex2i( left+1, top+1 );      glVertex2i( right-1, top+1 );
  glVertex2i( right-1, bottom-1 );  glVertex2i( left+1, bottom-1 );
  glEnd();

  draw_name( left+8, top+11 );

  if ( active ) 
    /*draw_active_box( left+4, left+string_width( name.c_str() )+12,       */
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
  
  if (currently_inside) {draw_pressed(); /* heavy black outline when pressed */ }
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
  int left, right, top, bottom;

  left   = 5;
  right  = w-left;
  top    = 3;
  bottom = 3+16;

  
  glColor3f( 0.0, 0.0, 0.0 );

  glBegin( GL_LINE_LOOP );
  glVertex2i( left, top );         glVertex2i( right, top );
  glVertex2i( right, bottom );     glVertex2i( left,bottom );
  glEnd();

  glBegin( GL_LINE_LOOP );
  glVertex2i( left+1, top+1 );         glVertex2i( right-1, top+1 );
  glVertex2i( right-1, bottom-1 );     glVertex2i( left+1,bottom-1 );
  glEnd();
}
