/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_control.cpp - top-level GLUI_Control class


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/

#include "glui.h"
#include "stdinc.h"

int _glui_draw_border_only = 0;


/**************************************** GLUI_Control::align() **************/

void    GLUI_Control::align( void )
{
  int  col_x, col_y, col_w, col_h, col_x_off, col_y_off;
  int  orig_x_abs;

  orig_x_abs = x_abs;

  /* Fix alignment bug relating to columns              */
  /*return;              */

  if ( NOT parent() )
    return;  /* Clearly this shouldn't happen, though */

  get_this_column_dims(&col_x, &col_y, &col_w, &col_h, 
		       &col_x_off, &col_y_off);

  if ( type == GLUI_CONTROL_COLUMN ) {
    /*		if ( this->prev() != NULL ) {
		((GLUI_Control*)prev())->get_this_column_dims(&col_x, &col_y, &col_w, &col_h, 
		&col_x_off, &col_y_off);
		
		x_abs = col_x + col_w;
		}
		else {
		x_abs = ((GLUI_Control*)parent())->x_abs;
		}
		*/
    return;
  }

  if ( alignment == GLUI_ALIGN_LEFT ) {
    x_abs = col_x + col_x_off;
  }
  else if ( alignment == GLUI_ALIGN_RIGHT ) {
    x_abs = col_x + col_w - col_x_off - this->w;
  }
  else if ( alignment == GLUI_ALIGN_CENTER ) {
    x_abs = col_x + (col_w - this->w) / 2;
  }

  if ( this->is_container ) {
    /***   Shift all child columns   ***/
    int delta = x_abs - orig_x_abs;

    GLUI_Control *node;
		
    node = (GLUI_Control*) this->first_child();
    while( node != NULL ) {
      if ( node->type == GLUI_CONTROL_COLUMN ) { 
	node->x_abs += delta;
      } 

      node = (GLUI_Control*) node->next();
    }
  }

}


/**************************************** GLUI_Control::pack() ************/
/* Recalculate positions and offsets  */

void    GLUI_Control::pack_old( int x, int y )
{
  GLUI_Control  *node;
  int            max_w, curr_y, curr_x, max_y;
  int            x_in = x, y_in =y;
  int            x_margin, y_margin_top, y_margin_bot;
  int            y_top_column;
  int            column_x;
  GLUI_Column   *curr_column = NULL;
  this->update_size();
  x_margin     = this->x_off;
  y_margin_top = this->y_off_top;
  y_margin_bot = this->y_off_bot;
  this->x_abs = x_in;
  this->y_abs = y_in;
  max_w  = -1;
  max_y  = -1;
  curr_x = this->x_abs + x_margin;
  curr_y = this->y_abs + y_margin_top;
  /*** Record start of this set of columns ***/
  y_top_column = curr_y;
  column_x     = 0;
  if ( this == glui->main_panel ) {
    x=x;
  }
  /*** Iterate over children, packing them first ***/
  node = (GLUI_Control*) this->first_child();
  while( node != NULL ) {
    if ( node->type == GLUI_CONTROL_PANEL ) { /* Pad some space above panels */
      curr_y += GLUI_ITEMSPACING;
    } 
    else if ( node->type == GLUI_CONTROL_COLUMN ) {
      curr_column = (GLUI_Column*) node;
      if ( 1 ) {
	column_x += max_w + 2 * x_margin;
	curr_x  += max_w + 2 * x_margin;
      }
      else {
	column_x += max_w + 0 * x_margin;
	curr_x  += max_w + 0 * x_margin;
      }
      /*node->pack( curr_x, curr_y );              */
      node->x_abs = curr_x;
      node->y_abs = y_top_column;
      node->w     = 2;
      node->h     = curr_y - y_top_column;
      curr_x  += x_margin * 3 + 40;
      curr_y  = y_top_column;
      max_w = 0;
      node = (GLUI_Control*) node->next();
      continue;
    }
    node->pack( curr_x, curr_y );
    if ( node->type == GLUI_CONTROL_PANEL )  /* Pad some space below panels */
      curr_y += GLUI_ITEMSPACING;
    curr_y  += node->h;
    if ( node->w > max_w ) {
      max_w = node->w;
      if ( curr_column != NULL )
	curr_column->w = max_w;
    }
    node = (GLUI_Control*) node->next();
    if ( node ) {
      curr_y += GLUI_ITEMSPACING;
    }
    if ( curr_y > max_y )
      max_y = curr_y;
  }
  if ( this->is_container ) {
    max_y += y_margin_bot;  /*** Add bottom border inside box */
    if ( this->first_child() ) {
      if ( this->type == GLUI_CONTROL_ROLLOUT ) {	
	/**  We don't want the rollout to shrink in width when it's
	  closed **/
	this->w = MAX(this->w, column_x + max_w + 2 * x_margin );
      }
      else {
	this->w        = column_x + max_w + 2 * x_margin;
      }
      this->h        = (max_y - y_in);
    }
    else  {            /* An empty container, so just assign default w & h */
      this->w        = GLUI_DEFAULT_CONTROL_WIDTH;
      this->h        = GLUI_DEFAULT_CONTROL_HEIGHT;
    }
    /** Expand panel if necessary (e.g., to include all the text in 
      a panel label) **/
    this->update_size();   
  }
}


/********************************* GLUT_Control::draw_recursive() **********/

void   GLUI_Control::draw_recursive( int x, int y )
{
  GLUI_Control *node;

  /*  printf( "%s %d\n", this->name.string, this->hidden );*/
  if ( NOT can_draw() )
    return;

  /*if ( 1 ) {  --  Debugging to check control width  
    glColor3f( 1.0, 0.0, 0.0 );
    glBegin( GL_LINES );
    glVertex2i( x_abs, y_abs );00
    glVertex2i( x_abs+w, y_abs );

    glEnd();
    }*/

  glMatrixMode( GL_MODELVIEW );
  glPushMatrix();

  glTranslatef( (float) this->x_abs + .5, 
		(float) this->y_abs + .5, 
		0.0 );

  if ( NOT _glui_draw_border_only ) {
    if ( NOT strcmp( name.string, "Rollout" ) ) {
    }

    this->draw( this->x_off, this->y_off_top );
  } 
  else 
  {
    if ( this->type == GLUI_CONTROL_COLUMN ) {
      /*   printf( "%s w/h:   %d/%d\n", (char*) name, w, h );              */
      /*w = 2;              */
    }

    /* The following draws the area of each control              */
    glColor3f( 1.0, 0.0, 0.0 );
    glBegin( GL_LINE_LOOP );
    glVertex2i( 0, 0 ); glVertex2i( w, 0 );
    glVertex2i( w, h ); glVertex2i( 0, h );
    glEnd();
  }
  glPopMatrix();
  
  node = (GLUI_Control*) first_child();
  while( node ) {
    node->draw_recursive( node->x_abs, node->y_abs );
    node = (GLUI_Control*) node->next();
  }
}


/******************************** GLUI_Control::set_to_glut_window() *********/
/*  Sets the current window to the glut window associated with this control  */

int   GLUI_Control::set_to_glut_window( void )
{
  int orig_window;

  if ( NOT glui) 
    return 1;

  orig_window = glutGetWindow();

  glutSetWindow( glui->get_glut_window_id());
  glDrawBuffer( GL_FRONT );

  return orig_window;
}


/************************************ GLUI_Control::restore_window() *********/

void   GLUI_Control::restore_window( int orig )
{
  if ( orig > 0 )
    glutSetWindow( orig );
}


/************************************* GLUI_Control::enable() ****************/

void   GLUI_Control::enable( void )
{
  GLUI_Control *node;

  enabled = true;
   
  if ( NOT glui )
    return;

  translate_and_draw_front();

  /*** Now recursively enable all buttons below it ***/
  node = (GLUI_Control*) first_child();
  while(node) {
    node->enable();
    node = (GLUI_Control*) node->next();
  }
}


/************************************ GLUI_Control::disable() ****************/

void   GLUI_Control::disable( void )
{
  GLUI_Control *node;

  enabled = false;
  
  if ( NOT glui )
    return;

  if ( glui->active_control == this )
    glui->disactivate_current_control();

  translate_and_draw_front();

  /*** Now recursively disable all buttons below it ***/
  node = (GLUI_Control*) first_child();
  while(node) {
    node->disable();
    node = (GLUI_Control*) node->next();
  }
}


/***************************************** GLUI_Control::set_font() **********/

void    GLUI_Control::set_font( void *new_font )
{
  int orig, state;

  font = new_font;

  /** Now redraw **/
  if ( NOT glui )
    return;

  orig = set_to_glut_window();
  state = glui->set_front_draw_buffer();
  glMatrixMode( GL_MODELVIEW );
  glPushMatrix();
  translate_to_origin();
  glDisable( GL_CULL_FACE );
  set_to_bkgd_color();
  glBegin( GL_QUADS );
  glVertex2i( 0,0 );   glVertex2i( w,0 );
  glVertex2i( w,h );   glVertex2i( 0,h );
  glEnd();
  glPopMatrix();
  glui->restore_draw_buffer(state);  
  restore_window(orig);

  translate_and_draw_front();
}


/************************************ GLUI_Control::draw_string() ************/

void   GLUI_Control::draw_string( char *text )
{
  _glutBitmapString( get_font(), text );
}


/****************************************** GLUI_Control::draw_char() ********/

void         GLUI_Control::draw_char( char c )
{
  glutBitmapCharacter( get_font(), (unsigned char)c );
}


/************************************* GLUI_Control::string_width() **********/

int          GLUI_Control::string_width( char *text )
{
  return _glutBitmapWidthString( get_font(), text );
}


/*************************************** GLUI_Control::char_width() **********/

int          GLUI_Control::char_width( char c )
{
  return glutBitmapWidth( get_font(), c );
}


/***************************************** GLUI_Control::get_font() **********/

void    *GLUI_Control::get_font( void )
{
  /*** Does this control have its own font? ***/
  if ( this->font != NULL )
    return this->font;
  
  /*** Does the parent glui have a font? ***/
  if ( glui )
    return glui->font;

  /*** Return the default font ***/
  return GLUT_BITMAP_HELVETICA_12;
}


/*************************************** GLUI_Control::draw_name() ***********/
/* This draws the name of the control as either black (if enabled), or       */
/* embossed if disabled.                                                     */

void         GLUI_Control::draw_name( int x, int y )
{
  if ( NOT can_draw() )
    return;

  if ( enabled ) {
    set_to_bkgd_color();
    glRasterPos2i(x+1, y+1);
    draw_string(name);
    glColor3b( 0, 0, 0 );
    glRasterPos2i(x, y);
    draw_string(name);
  }
  else {   /* Control is disabled - emboss the string */
    glColor3f( 1.0f, 1.0f, 1.0f );
    glRasterPos2i(x+1, y+1);
    draw_string(name);
    glColor3f( .4f, .4f, .4f );
    glRasterPos2i(x, y);
    draw_string(name);
  }
}


/*************************** GLUI_Control::translate_and_draw_front() ********/

void    GLUI_Control::translate_and_draw_front( void )
{
  int orig,state;

  if ( NOT can_draw() )
    return;

  orig = set_to_glut_window();
  state = glui->set_front_draw_buffer();
  glMatrixMode( GL_MODELVIEW );
  glPushMatrix();
  translate_to_origin();
  draw(0,0);
  glPopMatrix();
  glui->restore_draw_buffer(state);  
  restore_window(orig);
}


/**************************************** GLUI_Control::set_w() **************/

void      GLUI_Control::set_w( int new_w )
{
  w = new_w;
  update_size();  /* Make sure control is big enough to fit text */

  if ( NOT glui )
    return;
  
  glui->pack_controls();

  if ( glui->get_glut_window_id() != -1 ) {
    int orig = set_to_glut_window();
    glutReshapeWindow( glui->main_panel->w, glui->main_panel->h );
    glutPostRedisplay();
    /*    printf( "Requesting a reshape to window %d: %d %d\n", 
	  glutGetWindow(),
	  glui->main_panel->w, glui->main_panel->h );*/
    restore_window(orig);
  }
}


/**************************************** GLUI_Control::set_h() **************/

void      GLUI_Control::set_h( int new_h )
{
  h = new_h;
  update_size();  /* Make sure control is big enough to fit text */

  if ( NOT glui )
    return;
  
  glui->pack_controls();

  if ( glui->get_glut_window_id() != -1 ) {
    int orig = set_to_glut_window();
    glutReshapeWindow( glui->main_panel->w, glui->main_panel->h );
    glutPostRedisplay();
    restore_window(orig);
  }
}


/**************************************** GLUI_Control::set_alignment() ******/

void      GLUI_Control::set_alignment( int new_align )
{
  alignment = new_align;

  if ( glui ) {
    glui->align_controls( this );  

    if ( glui->get_glut_window_id() != -1 ) {
      int orig = set_to_glut_window();
      glutPostRedisplay();
      restore_window(orig);
    }
  }
}


/************************************* GLUI_Control::sync_live() ************/
/* Reads live variable and sets control to its current value                */
/* This function is recursive, and operates on control's children           */

void      GLUI_Control::sync_live( int recurse, int draw_it )
{
  GLUI_Node *node;
  int        sync_it=true;
  int        i;
  float      *fp;
  int        changed = false;

  /*** If this is currently active control, and mouse button is down,
    don't sync ***/
  if ( glui ) {
    if ( this == glui->active_control AND glui->mouse_button_down )
      sync_it = false;

    /*** Actually, just disable syncing if button is down ***/
    /*** Nope, go ahead and sync if mouse is down - this allows syncing in
      callbacks ***/
    if ( 0 ) { /* THIS CODE BELOW SHOULD NOT BE EXECUTED */
      if ( glui->mouse_button_down ) {
	/* printf( "Can't sync\n" );              */
	return;
      }
    }
  }

  /***  If this control has a live variable, we check its current value
    against the stored value in the control  ***/

  if ( ptr_val != NULL ) {
    if ( live_type == GLUI_LIVE_NONE OR NOT sync_it ) {
    }
    else if ( live_type == GLUI_LIVE_INT ) {
      if ( *((int*)ptr_val) != last_live_int ) {
	set_int_val( *((int*)ptr_val) );
	last_live_int = *((int*)ptr_val);
	changed = true;
      }
    }   
    else if ( live_type == GLUI_LIVE_FLOAT ) {
      if ( *((float*)ptr_val) != last_live_float ) {
	set_float_val( *((float*)ptr_val) );
	last_live_float = *((float*)ptr_val);
	changed = true;
      }
    } 
    else if ( live_type == GLUI_LIVE_TEXT ) {
      if ( strncmp((char*)ptr_val,last_live_text,sizeof(GLUI_String)) != 0 ) {
	set_text( (char*) ptr_val );
	strncpy( last_live_text, (char*) ptr_val, sizeof(GLUI_String));
	changed = true;
      }
    } 
    else if ( live_type == GLUI_LIVE_FLOAT_ARRAY ) { 
      /***  Step through the arrays, and see if they're the same  ***/
      
      fp = (float*) ptr_val;
      for ( i=0; i<float_array_size; i++ ) {
	if ( *fp != last_live_float_array[i] ) {
	  changed = true;
	  break;
	}
		
	fp++;
      }
      
      if ( changed == true) {
	fp = (float*) ptr_val;
	set_float_array_val( fp );
	for ( i=0; i<float_array_size; i++ ) {
	  last_live_float_array[i] = *fp;
	  fp++;
	}
      }
    }
    else if ( live_type == GLUI_LIVE_DOUBLE ) {
    }
  }

  /***  If this control is changed and we're supposed to be drawing, then
    draw it now    ***/
  if ( changed == true AND draw_it ) {
    translate_and_draw_front();
  }

  if ( recurse ) {
    /*** Now recursively output live vars for all children ***/
    
    node = this->first_child();
    while( node ) {
      ((GLUI_Control*) node)->sync_live(true, true);
      node = node->next();
    }

    if ( collapsible == true AND is_open == false ) {
      /** Here we have a collapsed control (e.g., a rollout that is closed **/
      /** We need to go in and sync all the collapsed controls inside      **/
       
      node = this->collapsed_node.first_child();
      while( node ) {
	((GLUI_Control*) node)->sync_live(true, false);
	node = node->next();
      }
    }
  }
}


/************************************ GLUI_Control::output_live() ************/
/* Writes current value of control to live variable.                         */

void      GLUI_Control::output_live( int update_main_gfx )
{
  int    i;
  float *fp;

  if ( ptr_val == NULL )
    return;

  if ( NOT live_inited ) 
    return;
  
  if ( live_type == GLUI_LIVE_NONE ) {
  }
  else if ( live_type == GLUI_LIVE_INT ) {
    *((int*)ptr_val) = int_val;
    last_live_int    = int_val;
  } 
  else if ( live_type == GLUI_LIVE_FLOAT ) {
    *((float*)ptr_val) = float_val;
    last_live_float    = float_val;
  } 
  else if ( live_type == GLUI_LIVE_TEXT ) {
    strncpy( (char*) ptr_val, text, sizeof(GLUI_String));
    strncpy( last_live_text,  text, sizeof(GLUI_String));
  } 
  else if ( live_type == GLUI_LIVE_FLOAT_ARRAY ) {
    fp = (float*) ptr_val;

    for( i=0; i<float_array_size; i++ ) {
      *fp                      = float_array_val[i];
      last_live_float_array[i] = float_array_val[i];
      fp++;
    }
  }
  else if ( live_type == GLUI_LIVE_DOUBLE ) {
  }

  /** Update the main gfx window? **/
  if ( update_main_gfx AND this->glui != NULL ) {
    this->glui->post_update_main_gfx();
  }
}


/******************************** GLUI_Control::execute_callback() **********/

void     GLUI_Control::execute_callback( void )
{
  int old_window;
  
  old_window = glutGetWindow();

  if ( glui AND glui->main_gfx_window_id != -1 ) 
    glutSetWindow( glui->main_gfx_window_id );

  if ( this->callback ) 
    this->callback( this->user_id );

  glutSetWindow( old_window );
}


/***************************** GLUI_Control::get_this_column_dims() **********/
/* Gets the x,y,w,h,and x/y offsets of the column to which a control belongs */

void    GLUI_Control::get_this_column_dims( int *col_x, int *col_y, 
					    int *col_w, int *col_h, 
					    int *col_x_off, int *col_y_off )
{
  GLUI_Control *node, *parent_ptr;
  int           parent_h, parent_y_abs;

  parent_ptr = (GLUI_Control*) parent();

  if ( parent_ptr==NULL )
    return;

  parent_h     = parent_ptr->h;
  parent_y_abs = parent_ptr->y_abs;
  
  if ( parent_ptr->type == GLUI_CONTROL_PANEL AND
       parent_ptr->int_val == GLUI_PANEL_EMBOSSED AND
       parent_ptr->name[0] != '\0' ) {
    parent_h -= GLUI_PANEL_EMBOSS_TOP;
    parent_y_abs += GLUI_PANEL_EMBOSS_TOP;
  }

  if ( 0 ) {
    GLUI_Node *first, *last, *curr;

    /**   Look for first control in this column   **/
    first = this;
    while (first->prev() != NULL AND 
	   ((GLUI_Control*)first->prev())->type != GLUI_CONTROL_COLUMN )
      first = first->prev();

    /**   Look for last control in this column    **/
    last = this;
    while ( last->next() != NULL AND 
	    ((GLUI_Control*)first->next())->type != GLUI_CONTROL_COLUMN )
      last = last->next();

    curr = first;
    int max_w = -1;
    do {
      if ( ((GLUI_Control*)curr)->w > max_w )
	max_w = ((GLUI_Control*)curr)->w;

      if ( curr == last )
	break;

      curr = curr->next();
    } while( curr != NULL );

    *col_x     = ((GLUI_Control*)first)->x_abs;
    *col_y     = ((GLUI_Control*)first)->y_abs;
    *col_w     = max_w;
    if ( parent() ) {
      *col_h     = ((GLUI_Control*)parent())->h;
      *col_x_off = ((GLUI_Control*)parent())->x_off;
    }
    else {
      *col_h = 10;
      *col_x_off = 0;
    }
    *col_y_off = 0;

    return;
  }

  if ( 1 ) {    /* IS THIS WRONG? */
    /*** Look for preceding column ***/
    node = (GLUI_Control*) this->prev();
    while( node ) {
      if ( node->type == GLUI_CONTROL_COLUMN ) {
	*col_x     = node->x_abs;
	*col_y     = parent_y_abs;
	*col_w     = node->w;
	*col_h     = parent_h;
	*col_x_off = node->x_off;
	*col_y_off = 0;

	return;
      }

      node = (GLUI_Control*) node->prev();
    }

    /*** Nope, Look for next column ***/
    node = (GLUI_Control*) this->next();
    while( node ) {
      if ( node->type == GLUI_CONTROL_COLUMN ) {
	*col_x     = parent_ptr->x_abs;
	*col_y     = parent_y_abs;
	*col_w     = node->x_abs - parent_ptr->x_abs;
	*col_h     = parent_h;
	*col_x_off = node->x_off;
	*col_y_off = 0;

	return;
      }

      node = (GLUI_Control*) node->next();
    }

    /***   This is single-column panel, so return panel dims   ***/
    *col_x     = parent_ptr->x_abs;
    *col_y     = parent_y_abs;
    *col_w     = parent_ptr->w;
    *col_h     = parent_h;
    *col_x_off = parent_ptr->x_off;
    *col_y_off = 0;
  }
}


/**************************************** GLUI_Control::init_live() **********/
/* Reads in  value of a live variable.  Called once, when ctrl is created   */

void     GLUI_Control::init_live( void )
{
  int    i;
  float *fp;

  if ( ptr_val == NULL )
    return;

  if ( live_type == GLUI_LIVE_NONE ) {
  }
  else if ( live_type == GLUI_LIVE_INT ) {
    set_int_val( *((int*)ptr_val) );
    last_live_int = *((int*)ptr_val);
  } 
  else if ( live_type == GLUI_LIVE_FLOAT ) {
    set_float_val( *((float*)ptr_val) );
    last_live_float = *((float*)ptr_val);
  } 
  else if ( live_type == GLUI_LIVE_TEXT ) {
    set_text( (char*) ptr_val );
    strncpy( last_live_text, (char*) ptr_val, sizeof(GLUI_String));
  }
  else if ( live_type == GLUI_LIVE_FLOAT_ARRAY ) {
    set_float_array_val( (float*) ptr_val );

    fp = (float*) ptr_val;

    for( i=0; i<float_array_size; i++ ) {
      last_live_float_array[i] = *fp;
      fp++;
    }

  }
  else if ( live_type == GLUI_LIVE_DOUBLE ) {
  }

  live_inited = true;
}


/**************************************** GLUI_Control::needs_idle() *********/
/* This method gets overloaded by specific classes, e.g. Spinner.            */
/* It returns whether or not a control needs to receive an idle event or not */
/* For example, a spinner only needs idle events when the user is holding    */
/* the mouse down in one of the arrows.  Otherwise, don't waste cycles       */
/* and OpenGL context switching by calling its idle.                         */

int  GLUI_Control::needs_idle( void )
{ 
  return false; 
};


/*********************************** GLUI_Control::~GLUI_Control() **********/

GLUI_Control::~GLUI_Control()
{
  GLUI_Control *node, *node_tmp;
  
  /*  printf( "destroying %s\n", this->name );              */

  node = (GLUI_Control*) this->first_child();
  while( node != NULL ) {
    /*    printf( "recursively destroying: '%s'\n", node->name );              */

    node_tmp = node;
    node = (GLUI_Control*) node->next();
    delete node_tmp;
  }
}


/************************** GLUI_Control::draw_box_inwards_outline() ********/

void     GLUI_Control::draw_box_inwards_outline( int x_min, int x_max, 
						 int y_min, int y_max )
{
  glBegin( GL_LINES );
  glColor3f( .5, .5, .5 );
  glVertex2i( x_min, y_min );     glVertex2i( x_max, y_min );
  glVertex2i( x_min, y_min );     glVertex2i( x_min, y_max );     

  glColor3f( 1., 1., 1. );
  glVertex2i( x_min, y_max );     glVertex2i( x_max, y_max );
  glVertex2i( x_max, y_max );     glVertex2i( x_max, y_min );

  if ( enabled )
    glColor3f( 0., 0., 0. );
  else
    glColor3f( .25, .25, .25 );
  glVertex2i( x_min+1, y_min+1 );     glVertex2i( x_max-1, y_min+1 );
  glVertex2i( x_min+1, y_min+1 );     glVertex2i( x_min+1, y_max-1 );

  glColor3f( .75, .75, .75 );
  glVertex2i( x_min+1, y_max-1 );     glVertex2i( x_max-1, y_max-1 );
  glVertex2i( x_max-1, y_max-1 );     glVertex2i( x_max-1, y_min+1 );
  glEnd();  
}


/************************************* GLUI_Control::draw_box() **********/

void     GLUI_Control::draw_box( int x_min, int x_max, int y_min, int y_max,
				 float r, float g, float b)
{
  if ( r == 1.0 AND g == 1.0 AND b == 1.0 AND NOT enabled AND glui ) {
    draw_bkgd_box( x_min, x_max, y_min, y_max );
    return;
  }

  glColor3f( r, g, b );
  glBegin( GL_QUADS );
  glVertex2i( x_min, y_min );       glVertex2i( x_max, y_min );
  glVertex2i( x_max, y_max );       glVertex2i( x_min, y_max );
  glEnd();
}


/************************************* GLUI_Control::draw_bkgd_box() **********/

void     GLUI_Control::draw_bkgd_box( int x_min, int x_max, int y_min, int y_max )
{
  set_to_bkgd_color();

  glBegin( GL_QUADS );
  glVertex2i( x_min, y_min );       glVertex2i( x_max, y_min );
  glVertex2i( x_max, y_max );       glVertex2i( x_min, y_max );
  glEnd();
}


/***************************** GLUI_Control::draw_active_area() ********/

void    GLUI_Control::draw_active_box( int x_min, int x_max, 
				       int y_min, int y_max )
{
  int orig;

  if ( NOT glui )
    return;

  orig = set_to_glut_window();

  if ( active ) {
    glEnable( GL_LINE_STIPPLE );
    glLineStipple( 1, 0x5555 );
    glColor3f( 0., 0., 0. );
  } else {
    set_to_bkgd_color();
  }

  glBegin( GL_LINE_LOOP );
  glVertex2i(x_min, y_min);      glVertex2i( x_max, y_min );
  glVertex2i(x_max, y_max);      glVertex2i( x_min, y_max );
  glEnd();
  
  glDisable( GL_LINE_STIPPLE );

  restore_window(orig);
}


/************************************** GLUI_Control::set_to_bkgd_color() ********/

void GLUI_Control::set_to_bkgd_color( void )
{
  if ( NOT glui )
    return;

  glColor3ub( glui->bkgd_color.r, glui->bkgd_color.g, glui->bkgd_color.b );
}


/******************************* GLUI_Control::set_float_array_val() ********/

void  GLUI_Control::set_float_array_val( float *array_ptr )
{
  int i;

  if ( array_ptr == NULL )
    return;

  for( i=0; i<float_array_size; i++ ) {
    float_array_val[i] = array_ptr[i];
  }

  /*** Output the live var, without updating the main gfx window ***/
  output_live(false);
}


/******************************* GLUI_Control::get_float_array_val() ********/

void  GLUI_Control::get_float_array_val( float *array_ptr )
{
  int i;

  if ( array_ptr == NULL )
    return;

  for( i=0; i<float_array_size; i++ ) {
    array_ptr[i] = float_array_val[i];
  }
}



/****************************** GLUI_Control::set_name() ********************/

void   GLUI_Control::set_name(char *string )
{
  strncpy((char*)name,string,sizeof(GLUI_String)); 

  if ( glui) 
    glui->refresh(); 
};





void    GLUI_Control::pack( int x, int y )
{
  GLUI_Control  *node;
  int            max_w, curr_y, curr_x, max_y;
  int            x_in = x, y_in =y;
  int            x_margin, y_margin_top, y_margin_bot;
  int            y_top_column;
  int            column_x;
  GLUI_Column   *curr_column = NULL;

  this->update_size();

  x_margin     = this->x_off;
  y_margin_top = this->y_off_top;
  y_margin_bot = this->y_off_bot;
  
  this->x_abs = x_in;
  this->y_abs = y_in;

  max_w  = 0;
  max_y  = 0;
  curr_x = this->x_abs + x_margin;
  curr_y = this->y_abs + y_margin_top;

  /*** Record start of this set of columns ***/

  y_top_column = curr_y;
  column_x     = curr_x;

  /*** Iterate over children, packing them first ***/

  node = (GLUI_Control*) this->first_child();
  while( node != NULL ) {
    if ( node->type == GLUI_CONTROL_PANEL ) { /* Pad some space above panels */
      curr_y += GLUI_ITEMSPACING;
    } 
    else if ( node->type == GLUI_CONTROL_COLUMN ) {
      curr_column = (GLUI_Column*) node;
      curr_x   += max_w + 1 * x_margin;
      column_x  = curr_x;

      node->x_abs = curr_x;
      node->y_abs = y_top_column;
      node->w     = 2;
      node->h     = curr_y - y_top_column;

      curr_x  += x_margin * 1;
      curr_y  = y_top_column;
      max_w = 0;

      node = (GLUI_Control*) node->next();
      continue;
    }
		
    node->pack( curr_x, curr_y );

    if ( node->type == GLUI_CONTROL_PANEL )  /* Pad some space below panels */
      curr_y += GLUI_ITEMSPACING;
    
    curr_y  += node->h;

    if ( node->w > max_w ) {
      max_w = node->w;
      if ( curr_column != NULL )
	curr_column->w = max_w + x_margin;
    }

    if ( curr_y > max_y ) {
      max_y = curr_y;
      if ( curr_column != NULL )
	curr_column->h = max_y - y_top_column;
    }

    node = (GLUI_Control*) node->next();
    
    if ( node ) {
      curr_y += GLUI_ITEMSPACING;
    }

  }

  if ( this->is_container ) {
    max_y += y_margin_bot;  /*** Add bottom border inside box */

    if ( this->first_child() ) {
      this->w        = column_x + max_w + 2 * x_margin - x_in;
      this->h        = (max_y - y_in);
    }
    else  {            /* An empty container, so just assign default w & h */
      if ( this->type != GLUI_CONTROL_ROLLOUT ) {
	this->w        = GLUI_DEFAULT_CONTROL_WIDTH;
	this->h        = GLUI_DEFAULT_CONTROL_HEIGHT;
      }
    }

    /** Expand panel if necessary (e.g., to include all the text in 
      a panel label) **/
    this->update_size();   


    /*** Now we step through the GLUI_Columns, setting the 'h'  ***/
    node = (GLUI_Control*) this->first_child();
    while( node != NULL ) {
      if ( node->type == GLUI_CONTROL_COLUMN ) {
	node->h = this->h - y_margin_bot - y_margin_top;
      }

      node = (GLUI_Control*) node->next();
    }
  }
}


/****************************** GLUI_Control::draw_emboss_box() ********/

void   GLUI_Control::draw_emboss_box(int x_min,int x_max,int y_min,int y_max)
{
  glLineWidth( 1.0 );
  glColor3f( 1.0, 1.0, 1.0 );

  glBegin( GL_LINE_LOOP );
  glVertex2i( x_min, y_min );    glVertex2i( x_max, y_min );
  glVertex2i( x_max, y_max );    glVertex2i( x_min, y_max );
  glEnd();
 
  glBegin( GL_LINE_LOOP );
  glVertex2i( x_min+1, y_min+1 );    glVertex2i( x_max-1, y_min+1 );
  glVertex2i( x_max-1, y_max-1 );    glVertex2i( x_min+1, y_max-1 );
  glEnd();
  
  glColor3f( .5, .5, .5 );
  glBegin( GL_LINE_LOOP );
  glVertex2i( x_min, y_min );
  glVertex2i( x_max-1, y_min );
  glVertex2i( x_max-1, y_max-1 );
  glVertex2i( x_min, y_max-1 );
  glEnd();
}


/*********************************** GLUI_Control::hide_internal() ********/
/** Sets hidden==true for this  control and all its siblings.             */
/**  If recurse is true, we go to children as well                       */

void         GLUI_Control::hide_internal( int recurse )
{
  GLUI_Node *node;

  node = (GLUI_Node *) this;
  while( node != NULL ) {
    ((GLUI_Control*)node)->hidden = true;

    if ( recurse AND node->first_child() != NULL )  
      ((GLUI_Control*) node->first_child())->hide_internal(true);
      
    node = node->next();
  }

  node = this->collapsed_node.first_child();
  while( node != NULL ) {
    ((GLUI_Control*)node)->hidden = true;

    if ( recurse AND node->first_child() != NULL )  
      ((GLUI_Control*) node->first_child())->hide_internal(true);
      
    node = node->next();
  }
}


/*********************************** GLUI_Control::unhide_internal() ********/
/** Sets hidden==false for this  control and all its siblings.             */
/**  If recurse is true, we go to children as well                       */

void         GLUI_Control::unhide_internal( int recurse )
{
  GLUI_Node *node;

  node = (GLUI_Node *) this;
  while( node != NULL ) {
    /*    printf( "unhide: %s [%d]\n", ((GLUI_Control*)node)->name.string, 
	    ((GLUI_Control*)node)->hidden );*/
    ((GLUI_Control*)node)->hidden = false;

    if ( recurse AND node->first_child() != NULL )  
      ((GLUI_Control*) node->first_child())->unhide_internal(true);
      
    node = node->next();
  }

  node = this->collapsed_node.first_child();
  while( node != NULL ) {
    ((GLUI_Control*)node)->hidden = false;

    if ( recurse AND node->first_child() != NULL )  
      ((GLUI_Control*) node->first_child())->unhide_internal(true);
      
    node = node->next();
  }
}
