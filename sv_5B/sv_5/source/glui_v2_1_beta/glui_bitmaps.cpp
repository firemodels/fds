/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui_bitmaps.cpp


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher

  This program is freely distributable without licensing fees and is
  provided without guarantee or warrantee expressed or implied. This
  program is -not- in the public domain.

*****************************************************************************/

#include "glui.h"
#include "stdinc.h"
#include <GL/glut.h>

int *bitmap_arrays[] = {
  glui_img_checkbox_0,
  glui_img_checkbox_1,
  glui_img_radiobutton_0,
  glui_img_radiobutton_1,
  glui_img_uparrow,
  glui_img_downarrow,
  glui_img_leftarrow,
  glui_img_rightarrow,
  glui_img_spinup_0,
  glui_img_spinup_1,
  glui_img_spindown_0,
  glui_img_spindown_1,
  glui_img_checkbox_0_dis,
  glui_img_checkbox_1_dis,
  glui_img_radiobutton_0_dis,
  glui_img_radiobutton_1_dis,
  glui_img_spinup_dis,
  glui_img_spindown_dis,
  glui_img_listbox_up,
  glui_img_listbox_down,
  glui_img_listbox_up_dis,
};


/************************************ GLUI_Bitmap::load_from_array() ********/

void GLUI_Bitmap::load_from_array( int *array )
{
  int i;

  w = 0;
  h = 0;

  if ( array == NULL )
    return;

  w = array[0];
  h = array[1];
  
  pixels = (unsigned char*) malloc( sizeof(unsigned char) * w * h * 3);

  for( i = 0; i<w*h*3; i++ )
    pixels[i] = (unsigned char) array[i+2];
}


/*********************************** GLUI_StdBitmaps::draw() *****************/

void GLUI_StdBitmaps::draw( int bitmap_num, int x, int y )
{
  //  int i;

  if ( bitmaps[bitmap_num].pixels != NULL ) {
    glRasterPos2f( (float)x+.5, (float)y + (float)bitmaps[bitmap_num].h + .5);
    glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );

    glDrawPixels( bitmaps[bitmap_num].w, bitmaps[bitmap_num].h,
		  GL_RGB, GL_UNSIGNED_BYTE, bitmaps[bitmap_num].pixels ); 

    /*    printf( "%d %d %d %d %d %d %d\n", 
	  bitmaps[bitmap_num].w, bitmaps[bitmap_num].h, 0, 0, 0, 0,
	  bitmaps[bitmap_num].pixels );

	  for( i=0; i<bitmaps[bitmap_num].w*bitmaps[bitmap_num].h*3; i++ ) {
	  printf( "%d ", bitmaps[bitmap_num].pixels[i] );
	  }
	  printf( "\n" );
	  */
    /*    for( i=0; i<bitmaps[bitmap_num].w; i++ ) {
	  for( j=0; j<bitmaps[bitmap_num].h; j++ ) {
	  if ( bitmaps[bitmap_num].pixels[(i+j*bitmaps[bitmap_num].w)*3]
	  > 128 )
	  putchar( '#' );
	  else 
	  putchar( ' ' );
	  }
	  putchar( '\n' );
	  }
	  fflush( stdout );
	  */
  }
}
