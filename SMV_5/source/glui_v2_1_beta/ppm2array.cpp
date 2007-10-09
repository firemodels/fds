#include "imagetype.h"
#include <stdio.h>

void main( int argc, char *argv[] )
{
  int        i;
  ImagePPM   img;
  FILE      *output;
  char       basename[200];

  if ( argc != 3 ) {
    fprintf( stderr, "USAGE: %s input.ppm output.cpp\n", argv[0] );
    return;
  }

  img.read( argv[1] );

  if ( img.valid ) {
    strcpy( basename, argv[1] );
    basename[ strlen(basename)-4 ] = '\0';
    
    output = fopen( argv[2], "w" );
    if ( NOT output ) {
      fprintf( stderr, "ERROR: File '%s' could not be opened for writing\n",
	       argv[2] );
      return;
    }

    img.flip_v(); /* Opengl bitmaps are specified bottom-to-top */

    fprintf( output, "\n\n");
    fprintf( output, "int %s[] = {", basename );
    fprintf( output, "    %d, %d,   /* width, height */\n", 
	     img.x_size, img.y_size );

    fprintf( output, "    " );  
    for( i=0; i<img.x_size * img.y_size; i++ ) {
      fprintf( output, "%3d,%3d,%3d,  ", 
	       img.pixels[i*3+0],
	       img.pixels[i*3+1],
	       img.pixels[i*3+2] );
      if ( (i%5) == 4 ) {
	fprintf( output, "\n    " );
      }
    }

    fprintf( output, "\n};\n" );
    fclose( output );
  }
  else {
    fprintf( stderr, "ERROR: Image '%s' invalid\n", argv[1] );
  }
  
}


