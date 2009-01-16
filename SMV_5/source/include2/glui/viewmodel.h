/*********************************************************************
  
  ViewModel.h

  The ViewModel class implements a camera model, minus the perspective
  transformation.  It maintains the coordinate system of the camera
  as three vectors (forward, up, and side), which are themselves
  derived from two points (eye and lookat).

  Apart from simplifying camera translation, this class provides
  three rotation methods: yaw (rotation about the up axis),
  roll (about the forward axis), and pitch (about the side axis).
  Also, these rotations can take place about the eye (in which
  case the eye location is fixed but the lookat point changes -
  like rotating your head), or about the lookat point (in which
  case your eye rotates around the point of interest).

  There is also a routine for rotating the eye about an arbitrary
  axis, which is quite useful in conjuction with the world up axis.
  
  This class is heavily-dependent on the vec3 class in
  the file algebra3.h.

  The function update() sets the side and forward vectors based
  on the lookat, eye, and up vectors.  Remember to call this
  function if you change the side or forward directly.

  Sample use:
     ViewModel vm;
     vm.eye.set( 0.0, 0.0, 20.0 );
     vm.lookat.set( 0.0, 0.0, 0.0 );
     vm.update();  // updates the three vectors 

     vm.move( 0.0, 2.0, 0.0 );  // Moves the eye and lookat 
				   by 2 units in the world 
				   Y direction 
     vm.roll( 5.0 );            // rolls 5 degrees about the forward
				   axis 
     vm.lookat_yaw( -25.0 );    // Yaws about the eye (lookat is
				   fixed) by -25 degrees 
     vm.load_to_openGL();       // Sets OpenGL modelview matrix 
     
     .. render ...


  ---------- List of member functions -----------------------------
  set_distance()   - Sets the distance from the eye to the lookat
  set_up()         - Sets the current up vector
  set_eye()        - Sets the current eye point
  set_lookat()     - Sets the current lookat point
  roll()           - Rolls about forward axis
  eye_yaw()        - Rotates eye about up vector
  eye_yaw_abs()    - Rotates eye about arbitrary axis
  eye_pitch()      - Rotates eye about side vector
  lookat_yaw()     - Rotates lookat about up vector
  lookat_pitch()   - Rotates lookat about side vector
  reset_up_level() - Resets the up vector (after unwanted rolls), and 
	                   sets the eye level with the lookat point
  move()           - Moves eye and lookat by some amount
  move_by_eye()    - Moves eye to new position, lookat follows
  move_by_lookat() - Moves lookat to new position, eye follows
  move_abs()       - Moves eye and lookat in world coordinates
  rot_about_eye()  - Rotates about the eye point by a given 4x4 matrix
  rot_about_lookat() - Rotates about the lookat point by a given 4x4 matrix  
  make_mtx()       - Constructs 4x4 matrix, used by load_to_openGL()
  load_to_openGL() - Loads current camera description in openGL
  load_to_openGL_noident() - Loads camera into OpenGL without first
                     resetting the OpenGL matrix to identity
  reset()          - Resets class values
  ViewModel()      - constructor 
  update()         - Recalculates side and forward based on eye,
                     lookat, and up
  dump()           - Prints class contents to a file


  - 1996 Paul Rademacher (rademach@cs.unc.edu)

*********************************************************************/

#ifndef _VIEW_MODEL_H_
#define _VIEW_MODEL_H_

#include "algebra3.h" /* Include algebra3.h first to avoid Windows 
			 problems */
#include <GL/glut.h>
#include <stdio.h>
#include <stdlib.h>

class ViewModel {
public:
  vec3    eye, lookat;
  vec3    up, side, forward;
  mat4    mtx;
  float   distance;

  /******************************* set_distance() ***********/
  /* This readjusts the distance from the eye to the lookat */
  /* (changing the eye point in the process)                */
  /* The lookat point is unaffected                         */
  void set_distance( float new_distance ) {
    if ( new_distance <= 0.0 )  /* Distance has to be positive */
      return;

    /* We find the current forward vector */
    forward = lookat - eye;
    forward.normalize();
    
    /* Set distance */
    distance = new_distance;

    /* Find new eye point */
    eye = lookat - forward * distance;
  }

  /******************************* set_up() ***************/
  void set_up( vec3 new_up ) {
    up = new_up;
    update();
  }

  void set_up( float x, float y, float z ) {
    set_up(vec3(x,y,z));
  }

  /******************************* set_eye() ***************/
  void set_eye( vec3 new_eye ) {
    eye = new_eye;
    update();
  }

  void set_eye( float x, float y, float z ) {
    set_eye(vec3(x,y,z));
  }

  /******************************* set_lookat() ***************/
  void set_lookat( vec3 new_lookat ) {
    lookat = new_lookat;
    update();
  }

  void set_lookat( float x, float y, float z ) {
    set_lookat(vec3(x,y,z));
  }


  /******************************* roll() *****************/
  /* Rotates about the forward vector                     */
  /* eye and lookat remain unchanged                      */
  void roll( float angle ) {
    mat4 rot = rotation3D( forward, angle );

    up = rot * up;

    update();
  }

  /******************************* eye_yaw() *********************/
  /* Rotates the eye about the lookat point, using the up vector */
  /* Lookat is unaffected                                 */
  void eye_yaw( float angle ) {
    vec3 eye_pt = eye - lookat; /* eye w/lookat at center */
    mat4 rot    = rotation3D( up, angle );
    
    eye_pt = rot * eye_pt;
    eye    = lookat + eye_pt;
     
    update();
  }

  /******************************* eye_yaw_abs() ******************/
  /* Rotates the eye about the lookat point, with a specific axis */
  /* Lookat is unaffected                                 */
  void eye_yaw_abs( float angle, vec3 axis ) {
    vec3 eye_pt   = eye  - lookat;   /* eye w/lookat at center */
    mat4 rot      = rotation3D( axis, angle );
    
    eye_pt  = rot * eye_pt;
    eye     = lookat + eye_pt;

    up      = rot * up;

    update();
  }


  /******************************* eye_pitch() ************/
  /* Rotates the eye about the side vector                */
  /* Lookat is unaffected                                 */
  void eye_pitch( float angle ) {
    vec3 eye_pt = eye - lookat; /* eye w/lookat at center */
    mat4 rot    = rotation3D( side, angle );
    
    eye_pt = rot * eye_pt;
    eye    = lookat + eye_pt;

    up = rot * up;
     
    update();
  }


  /******************************* lookat_yaw()************/
  /* This assumes the up vector is correct.               */
  /* Rotates the lookat about the side vector             */
  /* Eye point is unaffected                              */
  void lookat_yaw( float angle ) {
    vec3 lookat_pt = lookat - eye;   /* lookat w/eye at 
					center */
    mat4 rot = rotation3D( up, -angle );
    
    lookat_pt = rot * lookat_pt;
    lookat = eye + lookat_pt;

    update();
  }

  

  /******************************* lookat_pitch() *********/
  /* Rotates the lookat point about the side vector       */
  /* This assumes the side vector is correct.             */
  /* Eye point is unaffected                              */
  void lookat_pitch( float angle ) {
    vec3 lookat_pt = lookat - eye;   /* lookat w/eye at 
					center */
    mat4 rot = rotation3D( side, -angle );
    
    lookat_pt = rot * lookat_pt;
    lookat = eye + lookat_pt;

    up = rot * up;

    update();
  }


  /******************************* reset_up() ******************/
  /* Resets the up vector to a specified axis (0=X, 1=Y, 2=Z)  */
  /* Also sets the eye point level with the lookat point,      */
  /* along the specified axis                                  */
  void reset_up( int axis_num ) {
    float eye_distance = (lookat - eye).length();
    eye[axis_num] = lookat[axis_num]; 
    /* Bring eye to same level as lookat */

    vec3 vector = eye - lookat;
    vector.normalize();
    vector *= eye_distance;

    eye = lookat + vector;
    up.set( 0.0, 0.0, 0.0 );
    up[axis_num] = 1.0;

    update();
  }

  void reset_up( void ) {
    reset_up( VY ); /* Resets to the Y axis */
  }

  /******************************* move() *****************/
  /* Moves a specified distance in the forward, side, and up */
  /* directions.  This function does NOT move by world       */
  /* coordinates.  To move by world coords, use the move_abs */
  /* function.                                               */
  void move( float side_move, float up_move, float forw_move ) {
    eye += forward * forw_move;
    eye += side    * side_move;
    eye += up      * up_move;
    lookat += forward * forw_move;
    lookat += side    * side_move;
    lookat += up      * up_move;
    update();
  }

  void move( vec3 v ) { /* A vector version of the above command */
    move( v[VX], v[VY], v[VZ] );
  }

  /******************************* move_by_eye() ***********/
  /* Sets the eye point, AND moves the lookat point by the */
  /* same amount as the eye is moved.                      */
  void move_by_eye( vec3 new_eye ) {
    vec3 diff = new_eye - eye;

    eye    += diff;
    lookat += diff;

    update();
  }

  /******************************* move_by_lookat() *********/
  /* Sets the lookat point, AND moves the eye point by the  */
  /* same amount as the lookat is moved.                    */
  void move_by_lookat( vec3 new_lookat ) {
    vec3 diff = new_lookat - lookat;

    lookat += diff;
    eye    += diff;

    update();
  }

  /******************************* move_abs() *****************/
  /* Move the eye and lookat in world coordinates             */
  void move_abs( vec3 v ) {
    eye    += v;
    lookat += v;

    update();
  }

  /****************************** rot_about_eye() ************/
  /* Rotates the lookat point about the eye, based on a 4x4  */
  /* (pure) rotation matrix                                  */
  void rot_about_eye( mat4 rot ) {
    vec3  view = lookat - eye;

    view = rot * view;
    up   = rot * up;

    lookat = eye + view;

    update();
  }

  /****************************** rot_about_lookat() ************/
  /* Rotates the lookat point about the lookat, based on a 4x4  */
  /* (pure) rotation matrix                                  */
  void rot_about_lookat( mat4 rot ) {
    // NOT QUITE RIGHT YET


	 vec3 view = eye - lookat;

    view = rot * view;
    up   = rot * up;

    eye = lookat + view;

    update();
  }

  /******************************* make_mtx() *************/
  /* Constructs a 4x4 matrix - used by load_to_openGL()   */
  void make_mtx( void ) {
    update();

    mtx[0][0]=side[VX]; mtx[0][1]=up[VX]; mtx[0][2]=forward[VX]; mtx[0][3]=0.0;
    mtx[1][0]=side[VY]; mtx[1][1]=up[VY]; mtx[1][2]=forward[VY]; mtx[1][3]=0.0;
    mtx[2][0]=side[VZ]; mtx[2][1]=up[VZ]; mtx[2][2]=forward[VZ]; mtx[2][3]=0.0;
    mtx[3][0]=0.0;     mtx[3][1]=0.0;   mtx[3][2]= 0.0;        mtx[3][3]=1.0;
  }

  /******************************* load_to_openGL() *******/
  /* Sets the OpenGL modelview matrix based on the current */
  /* camera coordinates                                    */
  void load_to_openGL( void ) {
    int   i;
    mat4  m;

    make_mtx();

    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    glMultMatrixf( (float*) &mtx[0][0]);
    glTranslatef( -eye[VX], -eye[VY], -eye[VZ] );
  }

  /******************************* load_to_openGL_noident() ******/
  /* Multiplies the current camera matrix by the existing openGL */
  /* modelview matrix.  This is same as above function, but      */
  /* does not set the OpenGL matrix to identity first            */
  void load_to_openGL_noident( void ) {
    int   i;
    mat4  m;

    make_mtx();

    glMatrixMode( GL_MODELVIEW );
    glMultMatrixf( (float*) &mtx[0][0]);
    glTranslatef( -eye[VX], -eye[VY], -eye[VZ] );
  }

  /******************************* reset() ****************/
  /* Resets the parameters of this class                  */
  void reset( void ) {
    up.set( 0.0, 1.0, 0.0 );

    eye.set(0.0,0.0,10.0);
    lookat.set(0.0,0.0,0.0);

    mtx = identity3D();

    update();
  }

  /******************************* ViewModel() ************/
  /* Constructor                                          */
  ViewModel( void ) {
    reset();
  }

  /******************************* update() ****************/
  /* updates the view params.  Call this after making      */
  /* direct changes to the vectors or points of this class */
  void update( void ) {
    /* get proper side and forward vectors, and distance  */
    forward  = -(lookat - eye);
    distance = forward.length();
    forward /= distance;

    side     = up ^ forward;
    up       = forward ^ side;

    forward.normalize();
    up.normalize();
    side.normalize();
  }


  /******************************* dump() *******************/
  /* Prints the contents of this class to a file, typically */
  /* stdin or stderr                                        */
  void dump( FILE *output ) {
    fprintf( output, "Viewmodel: \n" );
    eye.print(    output, "  eye"    );
    lookat.print( output, "  lookat" );
    up.print(     output, "  up"     );
    side.print(   output, "  side"   );
    forward.print(output, "  forward");
    mtx.print(    output, "  mtx"    );
  }
};


#endif
