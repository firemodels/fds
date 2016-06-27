/****************************************************************************

  example5.cpp

  A GLUI program demonstrating subwindows, rotation controls,
	translation controls, and listboxes

  -----------------------------------------------------------------------
	   
  7/10/98 Paul Rademacher (rademach@cs.unc.edu)

****************************************************************************/

#include <string.h>
#include <GL/glut.h>
#include "glui.h"

float xy_aspect;
int   last_x, last_y;
float rotationX = 0.0, rotationY = 0.0;

/** These are the live variables passed into GLUI ***/
int   wireframe = 0;
int   obj_type = 1;
int   segments = 8;
int   segments2 = 8;
int   light0_enabled = 1;
int   light1_enabled = 1;
float light0_intensity = 1.0;
float light1_intensity = .4;
int   main_window;
float scale = 1.0;
int   show_sphere=1;
int   show_torus=1;
int   show_axes = 1;
int   show_text = 1;
float sphere_rotate[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
float torus_rotate[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
float view_rotate[16] = { 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };
float obj_pos[] = { 0.0, 0.0, 0.0 };
char *string_list[] = { "Hello World!", "Foo", "Testing...", "Bounding box: on" };
int   curr_string = 0;

/** Pointers to the windows and some of the controls we'll create **/
GLUI *glui, *glui2;
GLUI_Checkbox   *checkbox;
GLUI_Spinner    *spinner, *light0_spinner, *light1_spinner;
GLUI_RadioGroup *radio;
GLUI_Panel      *obj_panel;

/********** User IDs for callbacks ********/
#define LIGHT0_ENABLED_ID    200
#define LIGHT1_ENABLED_ID    201
#define LIGHT0_INTENSITY_ID  250
#define LIGHT1_INTENSITY_ID  251
#define ENABLE_ID            300
#define DISABLE_ID           301
#define SHOW_ID              302
#define HIDE_ID              303


/********** Miscellaneous global variables **********/

GLfloat light0_ambient[] =  {0.1f, 0.1f, 0.3f, 1.0f};
GLfloat light0_diffuse[] =  {.6f, .6f, 1.0f, 1.0f};
GLfloat light0_position[] = {.5f, .5f, 1.0f, 0.0f};

GLfloat light1_ambient[] =  {0.1f, 0.1f, 0.3f, 1.0f};
GLfloat light1_diffuse[] =  {.9f, .6f, 0.0f, 1.0f};
GLfloat light1_position[] = {-1.0f, -1.0f, 1.0f, 0.0f};

GLfloat lights_rotation[16] = {1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1 };

/**************************************** control_cb() *******************/
/* GLUI control callback                                                 */

void control_cb( int control )
{
  if ( control == LIGHT0_ENABLED_ID ) {
    if ( light0_enabled ) {
      glEnable( GL_LIGHT0 );
      light0_spinner->enable();
    }
    else {
      glDisable( GL_LIGHT0 ); 
      light0_spinner->disable();
    }
  }
  else if ( control == LIGHT1_ENABLED_ID ) {
    if ( light1_enabled ) {
      glEnable( GL_LIGHT1 );
      light1_spinner->enable();
    }
    else {
      glDisable( GL_LIGHT1 ); 
      light1_spinner->disable();
    }
  }
  else if ( control == LIGHT0_INTENSITY_ID ) {
    float v[] = { light0_diffuse[0],  light0_diffuse[1],
		  light0_diffuse[2],  light0_diffuse[3] };
    
    v[0] *= light0_intensity;
    v[1] *= light0_intensity;
    v[2] *= light0_intensity;

    glLightfv(GL_LIGHT0, GL_DIFFUSE, v );
  }
  else if ( control == LIGHT1_INTENSITY_ID ) {
    float v[] = { light1_diffuse[0],  light1_diffuse[1],
		  light1_diffuse[2],  light1_diffuse[3] };
    
    v[0] *= light1_intensity;
    v[1] *= light1_intensity;
    v[2] *= light1_intensity;

    glLightfv(GL_LIGHT1, GL_DIFFUSE, v );
  }
  else if ( control == ENABLE_ID )
  {
    glui2->enable();
  }
  else if ( control == DISABLE_ID )
  {
    glui2->disable();
  }
  else if ( control == SHOW_ID )
  {
    glui2->show();
  }
  else if ( control == HIDE_ID )
  {
    glui2->hide();
  }
}

/**************************************** myGlutKeyboard() **********/

void myGlutKeyboard(unsigned char Key, int x, int y)
{
  switch(Key)
  {
  case 27: 
  case 'q':
    exit(0);
    break;
  };
  
  glutPostRedisplay();
}


/***************************************** myGlutMenu() ***********/

void myGlutMenu( int value )
{
  myGlutKeyboard( value, 0, 0 );
}


/***************************************** myGlutIdle() ***********/

void myGlutIdle( void )
{
  /* According to the GLUT specification, the current window is 
     undefined during an idle callback.  So we need to explicitly change
     it if necessary */
  if ( glutGetWindow() != main_window ) 
    glutSetWindow(main_window);  

  /*  GLUI_Master.sync_live_all();  -- not needed - nothing to sync in this
                                       application  */

  glutPostRedisplay();
}

/***************************************** myGlutMouse() **********/

void myGlutMouse(int button, int button_state, int x, int y )
{
}


/***************************************** myGlutMotion() **********/

void myGlutMotion(int x, int y )
{
  glutPostRedisplay(); 
}

/**************************************** myGlutReshape() *************/

void myGlutReshape( int x, int y )
{
  int tx, ty, tw, th;
  GLUI_Master.get_viewport_area( &tx, &ty, &tw, &th );
  glViewport( tx, ty, tw, th );

  xy_aspect = (float)tw / (float)th;

  glutPostRedisplay();
}


/************************************************** draw_axes() **********/
/* Disables lighting, then draws RGB axes                                */

void draw_axes( float scale )
{
  glDisable( GL_LIGHTING );

  glPushMatrix();
  glScalef( scale, scale, scale );

  glBegin( GL_LINES );
 
  glColor3f( 1.0, 0.0, 0.0 );
  glVertex3f( .8f, 0.05f, 0.0 );  glVertex3f( 1.0, 0.25f, 0.0 ); /* Letter X */
  glVertex3f( 0.8f, .25f, 0.0 );  glVertex3f( 1.0, 0.05f, 0.0 );
  glVertex3f( 0.0, 0.0, 0.0 );  glVertex3f( 1.0, 0.0, 0.0 ); /* X axis      */

  glColor3f( 0.0, 1.0, 0.0 );
  glVertex3f( 0.0, 0.0, 0.0 );  glVertex3f( 0.0, 1.0, 0.0 ); /* Y axis      */

  glColor3f( 0.0, 0.0, 1.0 );
  glVertex3f( 0.0, 0.0, 0.0 );  glVertex3f( 0.0, 0.0, 1.0 ); /* Z axis    */
  glEnd();

  glPopMatrix();

  glEnable( GL_LIGHTING );
}


/***************************************** myGlutDisplay() *****************/

void myGlutDisplay( void )
{
  glClearColor( .9f, .9f, .9f, 1.0f );
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glFrustum( -xy_aspect*.04, xy_aspect*.04, -.04, .04, .1, 15.0 );

  glMatrixMode( GL_MODELVIEW );

  glLoadIdentity();
  glMultMatrixf( lights_rotation );
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
  
  glLoadIdentity();
  glTranslatef( 0.0, 0.0, -2.6f );
  glTranslatef( obj_pos[0], obj_pos[1], -obj_pos[2] ); 
  glMultMatrixf( view_rotate );

  glScalef( scale, scale, scale );

  /*** Now we render object, using the variables 'obj_type', 'segments', and
    'wireframe'.  These are _live_ variables, which are transparently 
    updated by GLUI ***/

  glPushMatrix();
  glTranslatef( -.5, 0.0, 0.0 );
  glMultMatrixf( sphere_rotate );
  if ( wireframe && show_sphere)
    glutWireSphere( .4, segments, segments );
  else if ( show_sphere )
    glutSolidSphere( .4, segments, segments );
  if ( show_axes )
    draw_axes(.52f);
  glPopMatrix();

  glPushMatrix();
  glTranslatef( .5, 0.0, 0.0 );
  glMultMatrixf( torus_rotate );
  if ( wireframe && show_torus )
    glutWireTorus( .15,.3,16,segments );
  else if ( show_torus )
    glutSolidTorus( .15,.3,16,segments );
  if ( show_axes )
    draw_axes(.52f);
  glPopMatrix();

  if ( show_text ) 
  {
    glDisable( GL_LIGHTING );  /* Disable lighting while we render text */
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluOrtho2D( 0.0, 100.0, 0.0, 100.0  );
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();
    glColor3ub( 0, 0, 0 );
    glRasterPos2i( 10, 10 );

    /*  printf( "text: %s\n", text );              */

    /*** Render the live character array 'text' ***/
    int i;
    for( i=0; i<(int)strlen( string_list[curr_string] ); i++ )
      glutBitmapCharacter( GLUT_BITMAP_HELVETICA_18, string_list[curr_string][i] );
  }

  glEnable( GL_LIGHTING );


  glutSwapBuffers(); 
}


/**************************************** main() ********************/

void main(int argc, char* argv[])
{
  /****************************************/
  /*   Initialize GLUT and create window  */
  /****************************************/

  glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
  glutInitWindowPosition( 50, 50 );
  glutInitWindowSize( 800, 600 );
 
  main_window = glutCreateWindow( "GLUI Example 5" );
  glutDisplayFunc( myGlutDisplay );
  GLUI_Master.set_glutReshapeFunc( myGlutReshape );  
  GLUI_Master.set_glutKeyboardFunc( myGlutKeyboard );
  GLUI_Master.set_glutSpecialFunc( NULL );
  GLUI_Master.set_glutMouseFunc( myGlutMouse );
  glutMotionFunc( myGlutMotion );

  /****************************************/
  /*       Set up OpenGL lights           */
  /****************************************/

  glEnable(GL_LIGHTING);
  glEnable( GL_NORMALIZE );

  glEnable(GL_LIGHT0);
  glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);

  glEnable(GL_LIGHT1);
  glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
  glLightfv(GL_LIGHT1, GL_POSITION, light1_position);

  /****************************************/
  /*          Enable z-buferring          */
  /****************************************/

  glEnable(GL_DEPTH_TEST);

  /****************************************/
  /*         Here's the GLUI code         */
  /****************************************/

  printf( "GLUI version: %3.2f\n", GLUI_Master.get_version() );

  /*** Create the side subwindow ***/
  glui = GLUI_Master.create_glui_subwindow( main_window, 
					    GLUI_SUBWINDOW_RIGHT );

  obj_panel = glui->add_rollout( "Properties", false );

  /***** Control for object params *****/

  checkbox = 
    glui->add_checkbox_to_panel( obj_panel, "Wireframe", &wireframe, 1, 
				 control_cb );
  spinner  = glui->add_spinner_to_panel( obj_panel, "Segments:",
					 GLUI_SPINNER_INT, &segments);
  spinner->set_int_limits( 3, 60 );
  spinner->set_alignment( GLUI_ALIGN_RIGHT );

  GLUI_Spinner *scale_spinner = 
    glui->add_spinner_to_panel( obj_panel, "Scale:",
				GLUI_SPINNER_FLOAT, &scale);
  scale_spinner->set_float_limits( .2f, 4.0 );
  scale_spinner->set_alignment( GLUI_ALIGN_RIGHT );


  /******** Add some controls for lights ********/

  GLUI_Rollout *roll_lights = glui->add_rollout( "Lights", false );

  GLUI_Panel *light0 = glui->add_panel_to_panel( roll_lights, "Light 1" );
  GLUI_Panel *light1 = glui->add_panel_to_panel( roll_lights, "Light 2" );

  glui->add_checkbox_to_panel( light0, "Enabled", &light0_enabled,
			       LIGHT0_ENABLED_ID, control_cb );
  light0_spinner = 
    glui->add_spinner_to_panel( light0, "Intensity:", GLUI_SPINNER_FLOAT,
				&light0_intensity, LIGHT0_INTENSITY_ID,
				control_cb );
  light0_spinner->set_float_limits( 0.0, 1.0 );

  glui->add_checkbox_to_panel( light1, "Enabled", &light1_enabled,
			       LIGHT1_ENABLED_ID, control_cb );
  light1_spinner = 
    glui->add_spinner_to_panel( light1, "Intensity:", GLUI_SPINNER_FLOAT,
				&light1_intensity, LIGHT1_INTENSITY_ID,
				control_cb );
  light1_spinner->set_float_limits( 0.0, 1.0 );


  /*** Add another rollout ***/
  GLUI_Rollout *options = glui->add_rollout( "Options", true );
  glui->add_checkbox_to_panel( options, "Draw sphere", &show_sphere );
  glui->add_checkbox_to_panel( options, "Draw torus", &show_torus );
  glui->add_checkbox_to_panel( options, "Draw axes", &show_axes );
  glui->add_checkbox_to_panel( options, "Draw text", &show_text );


  /**** Add listbox ****/
  glui->add_statictext( "" );
  GLUI_Listbox *list = glui->add_listbox( "Text:", &curr_string );
  int i;
  for( i=0; i<4; i++ )
    list->add_item( i, string_list[i] );

  glui->add_statictext( "" );


  /*** Disable/Enable buttons ***/
  glui->add_button( "Disable movement", DISABLE_ID, control_cb );
  glui->add_button( "Enable movement", ENABLE_ID, control_cb );
  glui->add_button( "Hide", HIDE_ID, control_cb );
  glui->add_button( "Show", SHOW_ID, control_cb );

  glui->add_statictext( "" );

  /****** A 'quit' button *****/
  glui->add_button( "Quit", 0,(GLUI_Update_CB)exit );


  /**** Link windows to GLUI, and register idle callback ******/
  
  glui->set_main_gfx_window( main_window );


  /*** Create the bottom subwindow ***/
  glui2 = GLUI_Master.create_glui_subwindow( main_window, 
					     GLUI_SUBWINDOW_BOTTOM );
  glui2->set_main_gfx_window( main_window );

  GLUI_Rotation *view_rot = glui2->add_rotation( "Objects", view_rotate );
  view_rot->set_spin( 1.0 );
  glui2->add_column( false );
  GLUI_Rotation *sph_rot = glui2->add_rotation( "Sphere", sphere_rotate );
  sph_rot->set_spin( .98 );
  glui2->add_column( false );
  GLUI_Rotation *tor_rot = glui2->add_rotation( "Torus", torus_rotate );
  tor_rot->set_spin( .98 );
  glui2->add_column( false );
  GLUI_Rotation *lights_rot = glui2->add_rotation( "Blue Light", lights_rotation );
  lights_rot->set_spin( .82 );
  glui2->add_column( false );
  GLUI_Translation *trans_xy = 
    glui2->add_translation( "Objects XY", GLUI_TRANSLATION_XY, obj_pos );
  trans_xy->set_speed( .005 );
  glui2->add_column( false );
  GLUI_Translation *trans_x = 
    glui2->add_translation( "Objects X", GLUI_TRANSLATION_X, obj_pos );
  trans_x->set_speed( .005 );
  glui2->add_column( false );
  GLUI_Translation *trans_y = 
    glui2->add_translation( "Objects Y", GLUI_TRANSLATION_Y, &obj_pos[1] );
  trans_y->set_speed( .005 );
  glui2->add_column( false );
  GLUI_Translation *trans_z = 
    glui2->add_translation( "Objects Z", GLUI_TRANSLATION_Z, &obj_pos[2] );
  trans_z->set_speed( .005 );

  /**** We register the idle callback with GLUI, *not* with GLUT ****/
  GLUI_Master.set_glutIdleFunc( myGlutIdle );

  /**** Regular GLUT main loop ****/
  
  glutMainLoop();
}

