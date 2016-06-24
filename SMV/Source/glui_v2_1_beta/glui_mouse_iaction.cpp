/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_mouse_iaction - GLUI Mouse Interaction control class


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/

#include "glui.h"
#include "stdinc.h"

/********************** GLUI_Mouse_Interaction::mouse_down_handler() ******/

int    GLUI_Mouse_Interaction::mouse_down_handler( int local_x, int local_y )
{
  int win_h = glutGet( GLUT_WINDOW_HEIGHT );

  /*	iaction_mouse_down_handler( local_x, local_y );              */
  iaction_mouse_down_handler( local_x-x_abs, local_y-y_abs );
  /*local_x-x_abs, ((glui->h-local_y)-y_abs) );              */
  
  return false;
}


/**************************** GLUI_Mouse_Interaction::mouse_up_handler() */

int    GLUI_Mouse_Interaction::mouse_up_handler( int local_x, int local_y, int inside )
{
  iaction_mouse_up_handler( local_x-x_abs, local_y-y_abs, inside );
  return false;
}


/****************************** GLUI_Mouse_Interaction::mouse_held_down_handler() ******/

int    GLUI_Mouse_Interaction::mouse_held_down_handler( int local_x, int local_y,
							int inside)
{  
  iaction_mouse_held_down_handler( local_x-x_abs, local_y-y_abs , inside );

  draw_active_area();

  /** Tell the main graphics window to update iteself **/
  if( glui )
    glui->post_update_main_gfx();

  execute_callback();

  return false;
}



/****************************** GLUI_Mouse_Interaction::draw() **********/

void    GLUI_Mouse_Interaction::draw( int x, int y )
{
  int orig;
  int text_width	= string_width( this->name );
  int x_left			= this->w/2 - text_width/2;
  
  if ( NOT glui )
    return;

  if ( NOT draw_active_area_only ) {
    orig = set_to_glut_window();
    draw_name( x_left, h-4 );
    restore_window(orig);

    draw_active_box( x_left-4, x_left+string_width( name.string )+4, 
		     h, h-14 );
  }

  draw_active_area();
}


/************************************ GLUI_Mouse_Interaction::update_size() **********/

void   GLUI_Mouse_Interaction::update_size( void )
{
  if ( NOT glui )
    return;

  int text_width = string_width( this->name );

  if ( w < text_width+6 )
    w = text_width+6;

  if ( h - 18 > w )
    w = h - 18;

  iaction_init();
}


/****************************** GLUI_Mouse_Interaction::special_handler() **********/

int    GLUI_Mouse_Interaction::special_handler( int key,int modifiers )
{
  int center_x, center_y;
  int drag_x, drag_y;

  center_x = w/2;
  center_y = (h-18)/2;
  drag_x   = 0;
  drag_y   = 0;

	
  if ( key == GLUT_KEY_LEFT )
    drag_x = -6;
  else if ( key == GLUT_KEY_RIGHT )
    drag_x = 6;
  else if ( key == GLUT_KEY_UP )
    drag_y = -6;
  else if ( key == GLUT_KEY_DOWN )
    drag_y = 6;

  if ( drag_x != 0 OR drag_y != 0 ) {
    mouse_down_handler( center_x, center_y );
    mouse_held_down_handler( center_x + drag_x, center_y + drag_y,true );
    mouse_up_handler( center_x + drag_x, center_y + drag_y, true );
  }

  return false;
}


/****************************** GLUI_Mouse_Interaction::draw_active_area() **********/

void    GLUI_Mouse_Interaction::draw_active_area( void )
{
  int orig;
  int win_h = glutGet( GLUT_WINDOW_HEIGHT ), win_w = glutGet(GLUT_WINDOW_WIDTH);
  
  if ( NOT glui )
    return;

  /*putchar( 'X' ); flushout;              */

  orig = set_to_glut_window();
	
  int text_height = 18; /* what a kludge              */

  int viewport_size = h-text_height;  /*MIN(w,h);              */

  glMatrixMode( GL_MODELVIEW );
  glPushMatrix();
  glLoadIdentity();
  glTranslatef( (float) win_w/2.0, (float) win_h/2.0, 0.0 );
  glRotatef( 180.0, 0.0, 1.0, 0.0 );
  glRotatef( 180.0, 0.0, 0.0, 1.0 );
  glTranslatef( (float) -win_w/2.0, (float) -win_h/2.0, 0.0 );

  glTranslatef( (float) this->x_abs + .5, (float) this->y_abs + .5, 0.0 );

  glTranslatef( (float)this->w/2.0, (float)viewport_size/2.0 + 2.0 , 0.0  );

  /***   Draw the interaction control's orthographic elements   ***/
  iaction_draw_active_area_ortho();	 
 
  /***   Setup and draw the interaction control's perspective elements   ***/

  /***  Set the viewport to just the square of the drawing area  ***/
  /* glViewport( this->x_abs , glui->main_panel->h - this->y_abs - this->h,*/
  /*glViewport( this->x_abs+1+(this->w/2-viewport_size/2),
    this->h-this->y_abs-viewport_size-1, 
    viewport_size, viewport_size );*/
	
  viewport_size -= 4;
  int offset = 0;
  if ( ((this->w-viewport_size) % 2) == 1 ) 
    offset = 1;

  glViewport( this->x_abs + (this->w-viewport_size)/2 + offset, 
	      win_h - this->y_abs - this->h + text_height, 
	      viewport_size, viewport_size );

  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glFrustum( -1.0*.08, 1.0*.08, -.08, .08, .1, 8.0 );
  glMatrixMode( GL_MODELVIEW );
  glPushMatrix();
  glLoadIdentity();
  glTranslatef( 0.0, 0.0, -3.2f );
  
  /*	glutSolidTeapot( 1.0 );              */
  iaction_draw_active_area_persp();

  glMatrixMode( GL_MODELVIEW );
  glPopMatrix();

  glui->set_viewport();
  glui->set_ortho_projection();

  glMatrixMode( GL_MODELVIEW );
  glPopMatrix();

  restore_window(orig);
}

