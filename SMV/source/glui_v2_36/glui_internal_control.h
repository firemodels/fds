/*
 Header file for use by GLUI controls.  
 Everything you need is right here.

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

*/
#ifndef __GLUI_INTERNAL_CONTROL_H
#define __GLUI_INTERNAL_CONTROL_H

/* This is the main GLUI external header */
#include "GL/glui.h"

/* Here's some utility routines */
#include "glui_internal.h"


/**
  A GLUI_Control-drawing sentinal object.
  On creation, saves the current draw buffer and window.
  On destruction, restores draw buffer and window.
  This is way nicer than calling save/restore manually.
*/
class GLUI_DrawingSentinal {
	int orig_buf, orig_win;
	GLUI_Control *c;
public:
	/** The constructor sets up the drawing system */
	GLUI_DrawingSentinal(GLUI_Control *c_);
	/** The destructor cleans up drawing back how it was */
	~GLUI_DrawingSentinal();
	
	// Do-nothing routine to avoid compiler warning about unused variable
	inline void avoid_warning(void) {}
};
/** Just drop a GLUI_DRAWINGSENTINAL_IDIOM at the start of your draw methods,
and they'll return if we can't be drawn, and 
automatically save and restore all needed state. 
*/
#define GLUI_DRAWINGSENTINAL_IDIOM  if (NOT can_draw()) return; GLUI_DrawingSentinal drawSentinal(this); drawSentinal.avoid_warning();


/** Return the time, in seconds. */
inline double GLUI_Time(void) {return 0.001*glutGet(GLUT_ELAPSED_TIME);}

#endif
