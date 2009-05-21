/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_rotation - GLUI_Rotation control class


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/

#include "glui.h"
#include "arcball.h"
#include "algebra3.h"

/*************************** GLUI_Rotation::iaction_mouse_down_handler() ***/

int    GLUI_Rotation::iaction_mouse_down_handler( int local_x, int local_y )
{
  copy_float_array_to_ball();

  init_ball();

  local_y = 2.0 * ball->center[1] - local_y;

  ball->mouse_down( local_x, local_y );

  /*	printf( "%d %d - %f %f\n", local_x, local_y, ball->center[0], ball->center[1] );              */

  copy_ball_to_float_array();

  spinning = false;

  return false;
}


/*********************** GLUI_Rotation::iaction_mouse_up_handler() **********/

int    GLUI_Rotation::iaction_mouse_up_handler( int local_x, int local_y, 
						int inside )
{
  copy_float_array_to_ball();

  ball->mouse_up();

  return false;
}


/******************* GLUI_Rotation::iaction_mouse_held_down_handler() ******/

int    GLUI_Rotation::iaction_mouse_held_down_handler( int local_x, int local_y,
						       int inside)
{  
  if ( NOT glui )
    return 0;

  copy_float_array_to_ball();

  local_y = 2.0 * ball->center[1] - local_y;

  /*	printf( "%d %d\n", local_x, local_y );              */

  ball->mouse_motion( local_x, local_y, 0, 
		     (glui->curr_modifiers & GLUT_ACTIVE_ALT) != 0, 
		     (glui->curr_modifiers & GLUT_ACTIVE_CTRL) != 0 );
 
  copy_ball_to_float_array();

  if ( can_spin )
    spinning = true;

  return false;
}


/******************** GLUI_Rotation::iaction_draw_active_area_persp() **************/

void    GLUI_Rotation::iaction_draw_active_area_persp( void )
{
  if ( NOT can_draw() )
    return;

  copy_float_array_to_ball();

  setup_texture();
  setup_lights();
	
  glEnable(GL_CULL_FACE );

  glMatrixMode( GL_MODELVIEW );
  glPushMatrix();

  mat4 tmp_rot = *ball->rot_ptr;
  glMultMatrixf( (float*) &tmp_rot[0][0] );

  /*** Draw the checkered box ***/
  /*glDisable( GL_TEXTURE_2D );              */
  draw_ball( 1.96 );

  glPopMatrix();

  glDisable( GL_TEXTURE_2D );
  glDisable( GL_LIGHTING );
  glDisable( GL_CULL_FACE );
}


/******************** GLUI_Rotation::iaction_draw_active_area_ortho() **********/

void    GLUI_Rotation::iaction_draw_active_area_ortho( void )
{
  if ( NOT can_draw() )
    return;

  /********* Draw emboss circles around arcball control *********/
  int k;
  float radius;
  radius = (float)(h-22)/2.0;  /*MIN((float)w/2.0, (float)h/2.0);              */
  glLineWidth( 1.0 );
  glBegin( GL_LINE_LOOP);
  for( k=0; k<60; k++ ) {
    float phi = 2*M_PI*(float)k/60.0;
    vec2 p( cos(phi) * (2.0 + radius), sin(phi) * (2.0 + radius));
    if ( p[1] < -p[0] ) 			glColor3ub( 128,128,128 );
    else					glColor3ub( 255,255,255 );
    glVertex2fv((float*)&p[0]);
  }
  glEnd();

  glBegin( GL_LINE_LOOP);
  for( k=0; k<60; k++ ) {
    float phi = 2*M_PI*(float)k/60.0;
    vec2 p( cos(phi) * (1.0 + radius), sin(phi) * (1.0 + radius));
    if ( enabled ) {
      if ( p[1] < -p[0] ) 			glColor3ub( 0,0,0);
      else					glColor3ub( 192,192,192);
    }
    else
    {
      if ( p[1] < -p[0] ) 			glColor3ub( 180,180,180);
      else					glColor3ub( 192,192,192);
    }
    glVertex2fv((float*)&p[0]);
  }
  glEnd();
}


/******************************** GLUI_Rotation::iaction_dump() **********/

void     GLUI_Rotation::iaction_dump( FILE *output )
{
}


/******************** GLUI_Rotation::iaction_special_handler() **********/

int    GLUI_Rotation::iaction_special_handler( int key,int modifiers )
{

  return false;
}

/********************************** GLUI_Rotation::init_ball() **********/

void  GLUI_Rotation::init_ball( void )
{
  /*printf( "%f %f %f", float( MIN(w/2,h/2)), (float) w/2, (float) h/2 );              */

  ball->set_params( vec2( (float)(w/2), (float)((h-18)/2)), 
		   (float) 2.0*(h-18) );
  /*ball->set_damping( .05 );              */
  /*float( MIN(w/2,h/2))*2.0  );              */
  /*	ball->reset_mouse();              */
}


/****************************** GLUI_Rotation::setup_texture() *********/

void    GLUI_Rotation::setup_texture( void )
{
  int i, j;
  int dark, light;   /*** Dark and light colors for ball checkerboard  ***/

#define CHECKBOARD_SIZE 64
  unsigned char texture_image[CHECKBOARD_SIZE] [CHECKBOARD_SIZE] [3];
  unsigned char c;
  for( i=0; i<CHECKBOARD_SIZE; i++ ) {
    for( j=0; j<CHECKBOARD_SIZE; j++ ) {
      if ( enabled ) {
	dark = 110;
	light = 220;
      }
      else {
	dark = glui->bkgd_color.r - 30;
	light = glui->bkgd_color.r;
      }

      c = ((((i&0x8)==0) ^ ((j&0x8))==0)) * light;
      if ( c == 0 )
	c = dark;
      texture_image[i][j][0] = c;
      texture_image[i][j][1] = c;
      texture_image[i][j][2] = c;
    }    
  }
  
  glColor3f( 1.0, 1.0, 1.0 );
  glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
  glEnable( GL_TEXTURE_2D);
  glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
  glTexImage2D( GL_TEXTURE_2D, 0, GL_RGB, CHECKBOARD_SIZE, 
		CHECKBOARD_SIZE, 0, GL_RGB, GL_UNSIGNED_BYTE,
		texture_image );

}


/****************************** GLUI_Rotation::setup_lights() ***********/

void    GLUI_Rotation::setup_lights( void )
{
  glEnable( GL_LIGHTING );
  /*  if ( enabled ) 
      glEnable( GL_LIGHTING );
      else
      glDisable( GL_LIGHTING );*/
  glEnable(GL_LIGHT0);
  glColorMaterial(GL_AMBIENT_AND_DIFFUSE, GL_FRONT_AND_BACK );
  glEnable(GL_COLOR_MATERIAL);
  GLfloat light0_ambient[] =  {0.2f, 0.2f, 0.2f, 1.0f};
  GLfloat light0_diffuse[] =  {1.f, 1.f, 1.0f, 1.0f};
  GLfloat light0_position[] = {-1.f, 1.f, 1.0f, 0.0f};
  glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
}


/****************************** GLUI_Rotation::draw_ball() **************/

void    GLUI_Rotation::draw_ball( float radius )
{
  if ( NOT can_draw() )
    return;

  if (quadObj == NULL)	quadObj = gluNewQuadric();
  if (quadObj) {
    gluQuadricDrawStyle(quadObj, GLU_FILL);
    gluQuadricNormals(quadObj, GLU_SMOOTH);
    gluQuadricTexture(quadObj, true );
    gluSphere(quadObj, radius, 16, 16);
  }
}


/****************************** GLUI_Rotation::reset() **********/

void  GLUI_Rotation::reset( void )
{
  ball->init(); /** reset quaternion, etc. **/
  ball->set_params( vec2( (float)(w/2), (float)((h-18)/2)), 
		   (float) 2.0*(h-18) );

  set_spin( this->damping );	

  copy_ball_to_float_array();

  translate_and_draw_front();

  output_live(true); /*** Output live and draw main grx window ***/
}


/****************************** GLUI_Rotation::needs_idle() *********/

int       GLUI_Rotation::needs_idle( void )
{
  return can_spin;
}


/****************************** GLUI_Rotation::idle() ***************/

void        GLUI_Rotation::idle( void )
{
  spinning = ball->is_spinning;

  if ( can_spin == true AND spinning == true ) {
    copy_float_array_to_ball();
    ball->idle();

    *ball->rot_ptr = *ball->rot_ptr * ball->rot_increment;

    mat4 tmp_rot;
    tmp_rot = *ball->rot_ptr;

    copy_ball_to_float_array();

    draw_active_area_only = true;
    translate_and_draw_front();
    draw_active_area_only = false;

    output_live(true); /** output live and update gfx **/
  }
  else { 
  }
}


/********************** GLUI_Rotation::copy_float_array_to_ball() *********/

void     GLUI_Rotation::copy_float_array_to_ball( void )
{
  int i;
  float *fp_src, *fp_dst;

  fp_src = &float_array_val[0];
  fp_dst = &((*ball->rot_ptr)[0][0]);

  for( i=0; i<16; i++ ) {
    *fp_dst = *fp_src;

    fp_src++;
    fp_dst++;
  }
}


/********************** GLUI_Rotation::copy_ball_to_float_array() *********/

void     GLUI_Rotation::copy_ball_to_float_array( void )
{
  mat4 tmp_rot;
  tmp_rot = *ball->rot_ptr;

  set_float_array_val( (float*) &tmp_rot[0][0] );
}


/************************ GLUI_Rotation::set_spin() **********************/

void   GLUI_Rotation::set_spin( float damp_factor )
{
  if ( damp_factor == 0.0 ) 
    can_spin = false;
  else
    can_spin = true;

  ball->set_damping( 1.0 - damp_factor );

  this->damping = damp_factor;
}


/************** GLUI_Rotation::GLUI_Rotation() ********************/


GLUI_Rotation::GLUI_Rotation( void ) 
{
  sprintf( name, "Rotation: %p", this );
  type                = GLUI_CONTROL_ROTATION;
  w                   = GLUI_ROTATION_WIDTH;
  h                   = GLUI_ROTATION_HEIGHT;
  can_activate        = true;
  live_type           = GLUI_LIVE_FLOAT_ARRAY;
  float_array_size    = 16;
  quadObj             = NULL;
  alignment           = GLUI_ALIGN_CENTER;
  can_spin            = false;
  spinning            = false;
  damping             = 0.0;
  ball                = new Arcball;

  reset();
}
