/**********************************************************************

  arcball.cpp


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

**********************************************************************/


#include "arcball.h"
#include <stdio.h>


/**************************************** Arcball::Arcball() ****/
/* Default (void) constructor for Arcball                         */

Arcball::Arcball( void ) 
{
  rot_ptr = &rot;

  init();
}


/**************************************** Arcball::Arcball() ****/
/* Takes as argument a mat4 to use instead of the internal rot  */

Arcball::Arcball( mat4 *mtx ) 
{
  rot_ptr = mtx;
}


/**************************************** Arcball::Arcball() ****/
/* A constructor that accepts the screen center and arcball radius*/

Arcball::Arcball( vec2 _center, float _radius )
{
  rot_ptr = &rot;

  init();
  set_params( _center, _radius );
}


/************************************** Arcball::set_params() ****/

void Arcball::set_params( vec2 _center, float _radius )
{
  center      = _center;
  radius      = _radius;
}


/*************************************** Arcball::init() **********/

void Arcball::init( void )
{
  center.set( 0.0, 0.0 );
  radius        = 1.0;
  q_now         = quat_identity();
  *rot_ptr      = identity3D();
  q_increment   = quat_identity();
  rot_increment = identity3D();
  is_mouse_down = false;
  is_spinning   = false;
  damp_factor   = 0.0;
  zero_increment = true;
}


/*********************************** Arcball::mouse_to_sphere() ****/

vec3 Arcball::mouse_to_sphere( vec2 p )
{
  float mag;
  vec2  v2 = (p - center) / radius;
  vec3  v3( v2[0], v2[1], 0.0 );
  vec3  axis;
  
  mag = v2*v2;
  
  if ( mag > 1.0 ) {
    v3.normalize();
  }
  else {
    v3[VZ] = sqrt( 1.0 - mag );
  }

  /* Now we add constraints - X takes precedence over Y */
  if ( constraint_x ) {
    v3 = constrain_vector( v3, vec3( 1.0, 0.0, 0.0 ));
  } else if ( constraint_y ) {
    v3 = constrain_vector( v3, vec3( 0.0, 1.0, 0.0 ));
  }
  
  return v3;
}


/************************************ Arcball::constrain_vector() ****/

vec3 Arcball::constrain_vector( vec3 vector, vec3 axis )
{
  return (vector-(vector*axis)*axis).normalize();
}

/************************************ Arcball::mouse_down() **********/

void Arcball::mouse_down( int x, int y )
{
  down_pt.set( (float)x, (float) y );
  is_mouse_down = true;

  q_increment   = quat_identity();
  rot_increment = identity3D();
  zero_increment = true;
}


/************************************ Arcball::mouse_up() **********/

void Arcball::mouse_up( void )
{
  q_now = q_drag * q_now;
  is_mouse_down = false;
}


/********************************** Arcball::mouse_motion() **********/

void Arcball::mouse_motion( int x, int y, int shift, int ctrl, int alt )
{
  /* Set the X constraint if CONTROL key is pressed, Y if ALT key */
  set_constraints( ctrl != 0, alt != 0 );

  vec2 new_pt( (float)x, (float) y );
  vec3 v0 = mouse_to_sphere( down_pt );
  vec3 v1 = mouse_to_sphere( new_pt );

  vec3 cross = v0^v1;
  
  q_drag.set( cross, v0 * v1 );
  
  //    *rot_ptr = (q_drag * q_now).to_mat4();
  mat4 temp = q_drag.to_mat4();
  *rot_ptr = *rot_ptr * temp;
    
  down_pt = new_pt;

	/* We keep a copy of the current incremental rotation (= q_drag) */
  q_increment   = q_drag;
  rot_increment = q_increment.to_mat4();

  set_constraints( false, false );

  if ( q_increment.s < .999999 ) {
    is_spinning = true;

    zero_increment = false;
  }
  else {
    is_spinning = false;
    zero_increment = true;
  }
}


/********************************** Arcball::mouse_motion() **********/

void Arcball::mouse_motion( int x, int y )
{
  mouse_motion( x, y, 0, 0, 0 );
}


/***************************** Arcball::set_constraints() **********/

void Arcball::set_constraints( Bool _constraint_x, Bool _constraint_y )
{
  constraint_x = _constraint_x;
  constraint_y = _constraint_y;
}

/***************************** Arcball::idle() *********************/

void  Arcball::idle( void )
{
  if ( is_mouse_down ) {
    is_spinning = false;
    zero_increment = true;
  }

  if ( damp_factor < 1.0 ) {
    q_increment.scale_angle( 1.0 - damp_factor );
  }

  rot_increment = q_increment.to_mat4();

  if ( q_increment.s >= .999999 ) {
    is_spinning = false;
    zero_increment = true;
  }
}


/************************ Arcball::set_damping() *********************/

void  Arcball::set_damping( float d )
{
  damp_factor = d;
}





