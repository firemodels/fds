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

GLUI_Panel::GLUI_Panel( GLUI_Node *parent, const char *name, int type )
{
  common_init();
  set_name( name );
  user_id    = -1;
  int_val    = type;

  parent->add_control( this );
}

/****************************** GLUI_Panel::draw() **********/

void    GLUI_Panel::draw( int x, int y )
{
  int top;
  GLUI_DRAWINGSENTINAL_IDIOM

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
    if ( parent_node == NULL || name == "" ) {
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
    if ( parent_node != NULL && name != "" ) { /* Only  draw non-null strings */
      int left = 7, height=GLUI_PANEL_NAME_DROP+1;
      int str_width;

      str_width = string_width(name);

      if ( glui )
	glColor3ubv(glui->bkgd_color);
      glDisable( GL_CULL_FACE );
      glBegin( GL_QUADS );
      glVertex2i( left-3, 0 );               glVertex2i( left+str_width+3, 0 );
      glVertex2i( left+str_width+3, height );  glVertex2i( left-3, height );
      glEnd();

      draw_name( left, GLUI_PANEL_NAME_DROP );
    }
  }

  glLineWidth( 1.0 );
}

/****************************** GLUI_Panel::set_name() **********/

void    GLUI_Panel::set_name( const char *new_name )
{
  name = new_name ? new_name : "";

  update_size();

  if ( glui )
    glui->refresh();
}


/****************************** GLUI_Panel::set_type() **********/

void    GLUI_Panel::set_type( int new_type )
{
  if ( new_type != int_val ) {
    int_val = new_type;
    update_size();
    redraw();
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

  if ( name != "" AND int_val == GLUI_PANEL_EMBOSSED ) {
    this->y_off_top = GLUI_YOFF + 8;
  }
  else {
    this->y_off_top = GLUI_YOFF;
  }
}
