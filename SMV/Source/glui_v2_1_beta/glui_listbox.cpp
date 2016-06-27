/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_listbox - GLUI_ListBox control class


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/

#include "glui.h"
#include "stdinc.h"

/****************************** GLUI_Listbox::mouse_down_handler() **********/

int    GLUI_Listbox::mouse_down_handler( int local_x, int local_y )
{
  return false;
}


/****************************** GLUI_Listbox::mouse_up_handler() **********/

int    GLUI_Listbox::mouse_up_handler( int local_x, int local_y, int inside )
{

  return false;
}


/****************************** GLUI_Listbox::mouse_held_down_handler() ******/

int    GLUI_Listbox::mouse_held_down_handler( int local_x, int local_y,
					      int inside)
{
  
  return false;
}


/****************************** GLUI_Listbox::key_handler() **********/

int    GLUI_Listbox::key_handler( unsigned char key,int modifiers )
{
  return false;
}


/****************************** GLUI_Listbox::draw() **********/

void    GLUI_Listbox::draw( int x, int y )
{
  int orig, name_x;

  if ( NOT can_draw() )
    return;

  orig = set_to_glut_window();

  /*  draw_active_area();              */

  name_x = MAX(text_x_offset - string_width(this->name) - 3,0);
  draw_name( name_x , 13);
  draw_box_inwards_outline( text_x_offset, w,
			    0, h );


  if ( NOT active ) {
    draw_box( text_x_offset+3, w-2, 2, h-2, 1.0, 1.0, 1.0 );
    if ( NOT enabled )
      glColor3b( 32, 32, 32 );
    else
      glColor3f( 0.0, 0.0, 0.0 );
    glRasterPos2i( text_x_offset+5, 13 );
    draw_string( curr_text.string );
  }
  else {
    draw_box( text_x_offset+3, w-2, 2, h-2, .0, .0, .6 );
    glColor3f( 1.0, 1.0, 1.0 );
    glRasterPos2i( text_x_offset+5, 13 );
    draw_string( curr_text.string );
  }


  if ( enabled ) {
    glui->std_bitmaps.
      draw(GLUI_STDBITMAP_LISTBOX_UP,
	   w-glui->std_bitmaps.bitmaps[GLUI_STDBITMAP_LISTBOX_UP].w-1,
	   2 );
  }
  else {
    glui->std_bitmaps.
      draw(GLUI_STDBITMAP_LISTBOX_UP_DIS,
	   w-glui->std_bitmaps.bitmaps[GLUI_STDBITMAP_LISTBOX_UP].w-1,
	   2 );
  }

  restore_window(orig);
}


/************************************ GLUI_Listbox::update_si() **********/

void   GLUI_Listbox::update_size( void )
{
  int text_size, delta;
  int item_text_size;

  if ( NOT glui )
    return;

  text_size = string_width( name );

  /*** Find the longest item string ***/
  item_text_size = 0;
  delta = 0;

  if ( text_x_offset < text_size +2 )
    delta = text_size+2-text_x_offset;

  text_x_offset += delta;
  if ( w < text_x_offset+MAX(GLUI_EDITTEXT_MIN_TEXT_WIDTH,item_text_size)+20)
    w = text_x_offset + MAX( GLUI_EDITTEXT_MIN_TEXT_WIDTH,item_text_size)+20;
}



/********************************* GLUI_Listbox::set_int_val() **************/

void    GLUI_Listbox::set_int_val( int new_val )
{
  /*  int_val = new_val;              */

  do_selection( new_val );

  /*** Update the variable we're (possibly) pointing to, and update the main gfx ***/
  output_live(true);
}


/**************************** GLUI_Listbox::draw_active_area() **************/

void    GLUI_Listbox::draw_active_area( void )
{
  int orig;

  if ( NOT can_draw() )
    return;

  orig = set_to_glut_window();

  restore_window(orig);
}


/**************************************** GLUI_Listbox::add_item() **********/

int  GLUI_Listbox::add_item( int id, char *new_text )
{
  GLUI_Listbox_Item *new_node = new GLUI_Listbox_Item;
  GLUI_Listbox_Item *head;

  strcpy( new_node->text, new_text );
  new_node->id = id;

  head = (GLUI_Listbox_Item*) items_list.first_child();
  new_node->link_this_to_parent_last( &items_list );

  if ( head == NULL ) {
    /***   This is first item added   ***/

    int_val       = id+1;  /** Different than id **/
    do_selection( id );
    last_live_int = id;

    if( glui )
      glui->post_update_main_gfx();
  }

  /*** Check if we need to increase control size ***/
  if ( w < text_x_offset + MAX( GLUI_EDITTEXT_MIN_TEXT_WIDTH, string_width( new_text ) ) + 20 ) {
    w = text_x_offset + MAX( GLUI_EDITTEXT_MIN_TEXT_WIDTH, string_width( new_text ) ) + 20;

    if ( glui )
      glui->refresh();

    /*		printf( "%s\n", new_text );              */
  }


  return true;
}


/************************************** GLUI_Listbox::delete_item() **********/

int  GLUI_Listbox::delete_item( char *text )
{
  GLUI_Listbox_Item *node = get_item_ptr( text );

  if ( node ) {
    node->unlink();
    delete node;
    return true;
  }
  else {
    return false;
  }
}


/************************************** GLUI_Listbox::delete_item() **********/

int  GLUI_Listbox::delete_item( int id )
{
  GLUI_Listbox_Item *node = get_item_ptr( id );

  if ( node ) {
    node->unlink();
    delete node;
    return true;
  }
  else {
    return false;
  }
}


/************************************** GLUI_Listbox::sort_items() **********/

int  GLUI_Listbox::sort_items( void )
{
  return false;
}


/********************************************* GLUI_Listbox::dump() **********/

void     GLUI_Listbox::dump( FILE *output )
{
  GLUI_Listbox_Item *item;

  /*  printf( "%p\n", (char*) name );              */

  fprintf( output, "Listbox: %s\n", (char*) name );

  item = (GLUI_Listbox_Item *) items_list.first_child();
  while( item ) {
    fprintf( output, "         %3d : %s\n", item->id, (char*) item->text );
    
    item = (GLUI_Listbox_Item *) item->next();
  }
}


/************************************ GLUI_Listbox::get_item_ptr() **********/

GLUI_Listbox_Item *GLUI_Listbox::get_item_ptr( char *text )
{
  GLUI_Listbox_Item *item;

  item = (GLUI_Listbox_Item *) items_list.first_child();
  while( item ) {
    if ( NOT strcmp( item->text, text ))
      return item;
    
    item = (GLUI_Listbox_Item *) item->next();
  }

  return NULL;
}


/************************************ GLUI_Listbox::get_item_ptr() **********/

GLUI_Listbox_Item *GLUI_Listbox::get_item_ptr( int id )
{
  GLUI_Listbox_Item *item;

  item = (GLUI_Listbox_Item *) items_list.first_child();
  while( item ) {
    if ( item->id == id )
      return item;
    
    item = (GLUI_Listbox_Item *) item->next();
  }

  return NULL;
}


/************************************ GLUI_Listbox::mouse_over() **********/

static void listbox_callback( int i )
{
  int old_val;

  if ( NOT GLUI_Master.curr_left_button_glut_menu OR 
       GLUI_Master.curr_left_button_glut_menu->type != GLUI_CONTROL_LISTBOX ) 
    return;

  old_val = ((GLUI_Listbox*)GLUI_Master.curr_left_button_glut_menu)->int_val;
  ((GLUI_Listbox*)GLUI_Master.curr_left_button_glut_menu)->set_int_val(i);

  /****   If value changed, execute callback   ****/
//  if ( old_val != 
//       ((GLUI_Listbox*)GLUI_Master.curr_left_button_glut_menu)->int_val ) {
//    ((GLUI_Listbox*)GLUI_Master.curr_left_button_glut_menu)->execute_callback();
//  }
  /****   always call execute callback   ****/
    ((GLUI_Listbox*)GLUI_Master.curr_left_button_glut_menu)->execute_callback();
}


/*************************************** GLUI_Listbox::mouse_over() **********/

int     GLUI_Listbox::mouse_over( int state, int x, int y )
{
  GLUI_Listbox_Item *item;

  /*  printf( "x/y:   %d/%d\n", x, y );              */

  if ( state AND enabled AND x > x_abs + text_x_offset) {
    /****  Build a GLUT menu for this listbox   ***/
    
    /*	printf( "%d %d\n", x, y );              */

    glut_menu_id = glutCreateMenu(listbox_callback);

    item = (GLUI_Listbox_Item *) items_list.first_child();
    while( item ) {
      glutAddMenuEntry( item->text, item->id );
      item = (GLUI_Listbox_Item *) item->next();
    }

    glutAttachMenu( GLUT_LEFT_BUTTON);
    
    GLUI_Master.set_left_button_glut_menu_control( this );
  }
  else if ( glut_menu_id != -1 ) {
    /*    printf( "OUT\n" );              */
    glutDetachMenu( GLUT_LEFT_BUTTON );
    glutDestroyMenu( glut_menu_id );
    glut_menu_id = -1;
  }

  return true;
}


/************************************ GLUI_Listbox::do_selection() **********/

int    GLUI_Listbox::do_selection( int item_num )
{
  GLUI_Listbox_Item *item, *sel_item;

  /***  Is this item already selected?  ***/
  if ( item_num == int_val )
    return false;

  sel_item = NULL;
  item     = (GLUI_Listbox_Item *) items_list.first_child();
  while( item ) {
    if ( item->id == item_num ) {
      sel_item = item;
      break;
    }
    
    item = (GLUI_Listbox_Item *) item->next();
  }

  if ( NOT sel_item )
    return false;

  /*  printf( "-> %s\n", (char*) sel_item->text );              */

  int_val = item_num;
  strcpy( curr_text.string, sel_item->text.string );

  translate_and_draw_front();

  return true;
}


/*********************************** GLUI_Listbox::~GLUI_Listbox() **********/

GLUI_Listbox::~GLUI_Listbox( )
{
  GLUI_Listbox_Item *item, *tmp_item;

  item = (GLUI_Listbox_Item *) items_list.first_child();
  while( item ) {
    tmp_item = item;

    delete item;
    
    item = (GLUI_Listbox_Item *) tmp_item->next();
  }
}


/****************************** GLUI_Listbox::special_handler() **********/

int    GLUI_Listbox::special_handler( int key,int modifiers )
{
  GLUI_Listbox_Item *node, *new_node;

  node     = get_item_ptr( int_val );
  new_node = NULL;

  if ( key == GLUT_KEY_DOWN ) {
    new_node = (GLUI_Listbox_Item*) node->next();
  }
  else if ( key == GLUT_KEY_UP ) {
    new_node = (GLUI_Listbox_Item*) node->prev();
  }
  else if ( key == GLUT_KEY_HOME ) {
    new_node = (GLUI_Listbox_Item*) items_list.first_child();
  }
  else if ( key == GLUT_KEY_END ) {
    new_node = (GLUI_Listbox_Item*) items_list.last_child();
  }

  if ( new_node != NULL AND new_node != node ) {
    node = new_node;
    set_int_val( node->id );
    execute_callback();
    return true;
  }
  else {
    return false;
  }
}


/************************* GLUI_Listbox::increase_width( void ) ***********/

void    GLUI_Listbox::increase_width( void )
{
}
