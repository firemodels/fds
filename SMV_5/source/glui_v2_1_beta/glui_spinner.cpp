/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_spinner.cpp - GLUI_Spinner class


  notes: 
     spinner does not explicitly keep track of the current value - this is all
        handled by the underlying edittext control
        -> thus, spinner->sync_live() has no meaning, nor spinner->output_live
	-> BUT, edittext will alter this spinner's float_val and int_val,
	   so that spinner->get/set will work


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/

  
#include "glui.h"
#include "stdinc.h"
#include <math.h>

/*static int __debug=0;              */

#define  GLUI_SPINNER_GROWTH_STEPS         800
#define  GLUI_SPINNER_MIN_GROWTH_STEPS     100
#define  GLUI_SPINNER_CALLBACK_INTERVAL    1

 
/****************************** spinner_edittext_callback() ******************/
/*   This function is not used anymore.  It has been replaced by directly    */
/*   Including an optional pointer to a spinner from an edittext box         */

void  spinner_edittext_callback( int id )
{
  GLUI_Spinner *spinner;

  putchar( '.' ); flushout;
  
  spinner = (GLUI_Spinner*) id;

  if ( NOT spinner )
    return;

  spinner->do_callbacks();
}


/****************************** GLUI_Spinner::mouse_down_handler() **********/

int    GLUI_Spinner::mouse_down_handler( int local_x, int local_y )
{
  this->state = find_arrow( local_x, local_y );

  /*  printf( "spinner: mouse down  : %d/%d   arrow:%d\n", local_x, local_y,
      find_arrow( local_x, local_y ));
      */

  if ( state != GLUI_SPINNER_STATE_UP AND state != GLUI_SPINNER_STATE_DOWN )
    return true;

  reset_growth();
  if ( can_draw() )
    draw_arrows();  

  /*** ints and floats behave a bit differently.  When you click on
    an int spinner, you expect the value to immediately go up by 1, whereas
    for a float it'll go up only by a fractional amount.  Therefore, we
    go ahead and increment by one for int spinners ***/
  if ( data_type == GLUI_SPINNER_INT ) {
    if ( state == GLUI_SPINNER_STATE_UP )
      edittext->set_float_val( edittext->float_val + 1.0 );
    else if ( state == GLUI_SPINNER_STATE_DOWN )
      edittext->set_float_val( edittext->float_val - .9 );
  }
  
  do_click();  
  
  return false;
}


/******************************** GLUI_Spinner::mouse_up_handler() **********/

int    GLUI_Spinner::mouse_up_handler( int local_x, int local_y, int inside )
{
  state = GLUI_SPINNER_STATE_NONE;

  /*  printf("spinner: mouse up  : %d/%d    inside: %d\n",local_x,local_y,inside);              */

  /*glutSetCursor( GLUT_CURSOR_INHERIT );              */
  glutSetCursor( GLUT_CURSOR_LEFT_ARROW );

  draw_arrows();

  /*  do_callbacks(); --- stub               */
  /*  if ( callback )               */
  /*  callback( this->user_id );              */
  
  return false;
}


/***************************** GLUI_Spinner::mouse_held_down_handler() ******/

int    GLUI_Spinner::mouse_held_down_handler( int local_x, int local_y,
					      int new_inside)
{
  int new_state;

  if ( state == GLUI_SPINNER_STATE_NONE )
    return false;

  /*  printf("spinner: mouse held: %d/%d    inside: %d\n",local_x,local_y,
      new_inside);
      */

  if ( state == GLUI_SPINNER_STATE_BOTH ) {   /* dragging? */
    do_drag( local_x, local_y );
  }
  else {                                      /* not dragging */
    new_state = find_arrow( local_x, local_y );
    
    if ( new_state == state ) {
      /** Still in same arrow **/
      do_click();
    }
    else {
      if ( new_inside OR 1) {
	/** The state changed, but we're still inside - that
	  means we moved off the arrow: begin dragging **/
	state = GLUI_SPINNER_STATE_BOTH;
      }
      else {
	/*** Here check y of mouse position to determine whether to 
	  drag ***/

	/* ... */
      }
    }

    /*** We switched to up/down dragging ***/
    if ( state == GLUI_SPINNER_STATE_BOTH ) {
      glutSetCursor( GLUT_CURSOR_UP_DOWN );
      last_x = local_x;
      last_y = local_y;

      /** If the spinner has limits, we reset the growth value, since
	reset_growth() will compute a new growth value for dragging
	vs. clicking.  If the spinner has no limits, then we just let the
	growth remain at whatever the user has incremented it up to **/
      if ( edittext->has_limits != GLUI_LIMIT_NONE )
	reset_growth();
    }

    if ( can_draw() )
      draw_arrows();
  }

  return false;
}


/****************************** GLUI_Spinner::key_handler() **********/

int    GLUI_Spinner::key_handler( unsigned char key,int modifiers )
{
  

  return true;
}


/****************************** GLUI_Spinner::draw() **********/

void    GLUI_Spinner::draw( int x, int y )
{
  int orig;

  if ( NOT can_draw() )
    return;

  orig = set_to_glut_window();

  if ( enabled ) {
    /*** Draw the up arrow either pressed or unrpessed ***/
    if ( state == GLUI_SPINNER_STATE_UP OR state == GLUI_SPINNER_STATE_BOTH )
      glui->std_bitmaps.draw( GLUI_STDBITMAP_SPINNER_UP_ON, 
			      w-GLUI_SPINNER_ARROW_WIDTH-1, 
			      GLUI_SPINNER_ARROW_Y);
    else
      glui->std_bitmaps.draw( GLUI_STDBITMAP_SPINNER_UP_OFF, 
			      w-GLUI_SPINNER_ARROW_WIDTH-1,
			      GLUI_SPINNER_ARROW_Y);

    /*** Draw the down arrow either pressed or unrpessed ***/
    if (state == GLUI_SPINNER_STATE_DOWN OR state == GLUI_SPINNER_STATE_BOTH)
      glui->std_bitmaps.draw( GLUI_STDBITMAP_SPINNER_DOWN_ON, 
			      w-GLUI_SPINNER_ARROW_WIDTH-1, 
			      GLUI_SPINNER_ARROW_HEIGHT+GLUI_SPINNER_ARROW_Y);
    else
      glui->std_bitmaps.draw( GLUI_STDBITMAP_SPINNER_DOWN_OFF, 
			      w-GLUI_SPINNER_ARROW_WIDTH-1,
			      GLUI_SPINNER_ARROW_HEIGHT+GLUI_SPINNER_ARROW_Y);
  }
  else {  /**** The spinner is disabled ****/
    glui->std_bitmaps.draw( GLUI_STDBITMAP_SPINNER_UP_DIS, 
			    w-GLUI_SPINNER_ARROW_WIDTH-1, 
			    GLUI_SPINNER_ARROW_Y);
    glui->std_bitmaps.draw( GLUI_STDBITMAP_SPINNER_DOWN_DIS, 
			    w-GLUI_SPINNER_ARROW_WIDTH-1, 
			    GLUI_SPINNER_ARROW_HEIGHT+GLUI_SPINNER_ARROW_Y);
  }

  if ( active ) {
    glColor3ub( 0, 0, 0 );
    glEnable( GL_LINE_STIPPLE );
    glLineStipple( 1, 0x5555 );
  }
  else {
    glColor3ub( glui->bkgd_color.r,glui->bkgd_color.g,glui->bkgd_color.b );
  } 

  glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
  glDisable( GL_CULL_FACE );
  glBegin( GL_QUADS );
  glVertex2i( w-GLUI_SPINNER_ARROW_WIDTH-2, 0 );
  glVertex2i( w, 0 );
  glVertex2i( w, h );
  glVertex2i( w-GLUI_SPINNER_ARROW_WIDTH-2, h );
  glEnd();
  glDisable( GL_LINE_STIPPLE );  
  glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

  restore_window( orig );
}





/********************************* GLUI_Spinner::special_handler() **********/

int    GLUI_Spinner::special_handler( int key,int modifiers )
{
  if ( key == GLUT_KEY_UP ) {    /** Simulate a click in the up arrow **/
    mouse_down_handler( x_abs + w - GLUI_SPINNER_ARROW_WIDTH + 1,
			y_abs + GLUI_SPINNER_ARROW_Y+1 );
    mouse_up_handler( x_abs + w - GLUI_SPINNER_ARROW_WIDTH + 1,
		      y_abs + GLUI_SPINNER_ARROW_Y+1, true );
  }
  else if ( key == GLUT_KEY_DOWN ) {  /** Simulate a click in the up arrow **/
    mouse_down_handler(x_abs + w - GLUI_SPINNER_ARROW_WIDTH + 1,
		       y_abs+GLUI_SPINNER_ARROW_Y+1+GLUI_SPINNER_ARROW_HEIGHT);
    mouse_up_handler( x_abs + w - GLUI_SPINNER_ARROW_WIDTH + 1,
		      y_abs+GLUI_SPINNER_ARROW_Y+1 +GLUI_SPINNER_ARROW_HEIGHT,
		      true );
  }
  else if ( key == GLUT_KEY_HOME ) {  /** Set value to limit top - 
					or increment by 10 **/
  }
  else if ( key == GLUT_KEY_END ) {  
  }

  return true;
}


/******************************* GLUI_Spinner::set_float_val() ************/

void   GLUI_Spinner::set_float_val( float new_val )
{ 
  if ( NOT edittext )
    return;

  edittext->set_float_val( new_val );
}


/********************************** GLUI_Spinner::set_int_val() ************/

void   GLUI_Spinner::set_int_val( int new_val )
{
  if ( NOT edittext )
    return;

  edittext->set_int_val( new_val );
}


/************************************ GLUI_Spinner::update_size() **********/

void   GLUI_Spinner::update_size( void )
{
  /*edittext->w = this->w - GLUI_SPINNER_ARROW_WIDTH-3;              */
  this->w = edittext->w + GLUI_SPINNER_ARROW_WIDTH + 3;
}
 

/************************************ GLUI_Spinner::find_arrow() ************/

int    GLUI_Spinner::find_arrow( int local_x, int local_y )
{
  local_x -= x_abs; 
  local_y -= y_abs;
 
  if ( local_x >= (w - GLUI_SPINNER_ARROW_WIDTH) AND
       local_x <= w ) {

    if ( local_y >= GLUI_SPINNER_ARROW_Y AND 
	 local_y <= (GLUI_SPINNER_ARROW_Y+GLUI_SPINNER_ARROW_HEIGHT) )
      return GLUI_SPINNER_STATE_UP;

    if ( local_y >= GLUI_SPINNER_ARROW_Y+GLUI_SPINNER_ARROW_HEIGHT AND 
	 local_y <= (GLUI_SPINNER_ARROW_Y+GLUI_SPINNER_ARROW_HEIGHT*2) )
      return GLUI_SPINNER_STATE_DOWN;

  }

  return GLUI_SPINNER_STATE_NONE;
}


/****************************** GLUI_Spinner::draw_arrows() **********/

void    GLUI_Spinner::draw_arrows( void )
{
  if ( NOT can_draw() )
    return;

  translate_and_draw_front();
}


/***************************************** GLUI_Spinner::do_click() **********/

void    GLUI_Spinner::do_click( void )
{
  int    direction;
  float  incr;
  float  modifier_factor;

  if ( state == GLUI_SPINNER_STATE_UP )
    direction = +1;
  else if ( state == GLUI_SPINNER_STATE_DOWN )
    direction = -1;

  increase_growth();

  modifier_factor = 1.0;
  if ( glui ) {
    if ( glui->curr_modifiers & GLUT_ACTIVE_SHIFT ) 
      modifier_factor = 100.0f;
    else if ( glui->curr_modifiers & GLUT_ACTIVE_CTRL ) 
      modifier_factor = .01f;
  }

  if ( this->data_type == GLUI_SPINNER_FLOAT OR 1) {
    incr = growth * direction * modifier_factor * user_speed;
    edittext->set_float_val( edittext->float_val + incr );
    /** Remember, edittext mirrors the float and int values ***/
  }

  /*** Now update live variable and do callback.  We don't want
    to do the callback on each iteration of this function, just on every 
    i^th iteration, where i is given by GLUI_SPINNER_CALLBACK_INTERVAL ****/
  callback_count++;
  if ( (callback_count % GLUI_SPINNER_CALLBACK_INTERVAL ) == 0 )
    do_callbacks();
}


/***************************************** GLUI_Spinner::do_drag() **********/

void    GLUI_Spinner::do_drag( int x, int y )
{
  int   delta_y;
  float incr, modifier_factor;
  /* int delta_x;              */

  modifier_factor = 1.0f;
  if ( glui ) {
    if ( glui->curr_modifiers & GLUT_ACTIVE_SHIFT ) 
      modifier_factor = 100.0f;
    else if ( glui->curr_modifiers & GLUT_ACTIVE_CTRL ) 
      modifier_factor = .01f;
  }

  /*  delta_x = x - last_x;              */
  delta_y = -(y - last_y);
 
  if ( this->data_type == GLUI_SPINNER_FLOAT OR 1 ) {
    incr = growth * delta_y * modifier_factor * user_speed;
    edittext->set_float_val( edittext->float_val + incr );
    /** Remember, edittext mirrors the float and int values ***/
  }

  last_x = x;
  last_y = y;

  /*** Now update live variable and do callback.  We don't want
    to do the callback on each iteration of this function, just on every 
    i^th iteration, where i is given by GLUI_SPINNER_CALLBACK_INTERVAL ****/

  callback_count++;
  if ( (callback_count % GLUI_SPINNER_CALLBACK_INTERVAL ) == 0 )
    do_callbacks();
}


/***************************************** GLUI_Spinner::needs_idle() ******/

int     GLUI_Spinner::needs_idle( void )
{
  if  (state == GLUI_SPINNER_STATE_UP OR state == GLUI_SPINNER_STATE_DOWN ) {
    return true;
  }
  else {
    return false;
  }
}

/***************************************** GLUI_Spinner::idle() **********/

void    GLUI_Spinner::idle( void )
{
  if ( NOT needs_idle() )
    return;
  else
    do_click();
}


/************************************ GLUI_Spinner::do_callbacks() **********/

void    GLUI_Spinner::do_callbacks( void )
{
  /*** This is not necessary, b/c edittext automatically updates us ***/
  if ( NOT edittext )
    return;
  this->float_val = edittext->float_val;
  this->int_val   = edittext->int_val;
  /*    *******************************************/

  if ( NOT first_callback ) {
    if ( data_type == GLUI_SPINNER_INT AND int_val == last_int_val ) {
      return;
    }
    
    if ( data_type == GLUI_SPINNER_FLOAT AND float_val == last_float_val ) {
      return;
    }
  }
  
  this->execute_callback();

  last_int_val   = int_val;
  last_float_val = float_val;
  first_callback = false;
}


/********************************* GLUI_Spinner::set_float_limits() *********/

void GLUI_Spinner::set_float_limits( float low, float high, int limit_type )
{
  if ( NOT edittext ) 
    return;

  edittext->set_float_limits( low, high, limit_type );
}


/*********************************** GLUI_Spinner::set_int_limits() *********/

void   GLUI_Spinner::set_int_limits( int low, int high, int limit_type )
{
  if ( NOT edittext ) 
    return;
  
  edittext->set_int_limits( low, high, limit_type );
}
 

/*********************************** GLUI_Spinner:reset_growth() *************/

void    GLUI_Spinner::reset_growth( void )
{
  float lo, hi;

  if ( edittext->has_limits == GLUI_LIMIT_NONE ) {
    if ( data_type == GLUI_SPINNER_FLOAT )
      growth = sqrt(ABS(edittext->float_val)) * .05f;
    else if ( data_type == GLUI_SPINNER_INT )
      growth = .4f; 
  }
  else {
    if ( data_type == GLUI_SPINNER_FLOAT ) {
      lo = edittext->float_low;
      hi = edittext->float_high;
      growth = (hi-lo) / GLUI_SPINNER_GROWTH_STEPS;
    }
    else if ( data_type == GLUI_SPINNER_INT ) {
      lo = (float) edittext->int_low;
      hi = (float) edittext->int_high;
      
      growth = (hi-lo) / GLUI_SPINNER_GROWTH_STEPS;
    }
  }

  if ( growth == 0.0f )
    growth = .001f;
}


/******************************* GLUI_Spinner:increase_growth() *************/

void    GLUI_Spinner::increase_growth( void )
{
  float hi,lo;

  if ( data_type == GLUI_SPINNER_FLOAT ) {
    lo = edittext->float_low;
    hi = edittext->float_high;
  }
  else if ( data_type == GLUI_SPINNER_INT ) {
    lo = (float) edittext->int_low;
    hi = (float) edittext->int_high;
  }
 
  if ( growth < (hi-lo) / GLUI_SPINNER_MIN_GROWTH_STEPS )
    growth *= growth_exp;

  /*  printf( "growth: %f\n", growth );              */
}


/*************************************** GLUI_Spinner:get_text() *************/

char    *GLUI_Spinner::get_text( void )
{ 
  if (edittext) 
    return edittext->text; 
  else 
    return ""; 
}


/********************************** GLUI_Spinner:get_float_val() *************/

float    GLUI_Spinner::get_float_val( void )
{
  if (edittext) 
    return edittext->float_val; 
  else 
    return 0.0f; 
}


/********************************** GLUI_Spinner:get_int_val() *************/

int    GLUI_Spinner::get_int_val( void )
{
  if (edittext) 
    return edittext->int_val; 
  else 
    return 0; 
}


