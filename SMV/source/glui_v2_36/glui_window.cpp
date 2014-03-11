/*

  glui_window.cpp - GLUI_Button control class

  GLUI User Interface Toolkit 
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

*/

#include "GL/glui.h"
#include "glui_internal.h"

GLUI_Glut_Window::GLUI_Glut_Window()
:   GLUI_Node(),

	glut_window_id(0),
	glut_keyboard_CB(NULL),
	glut_special_CB(NULL),
	glut_reshape_CB(NULL),
	glut_passive_motion_CB(NULL),
	glut_mouse_CB(NULL),
	glut_visibility_CB(NULL),
	glut_motion_CB(NULL),
	glut_display_CB(NULL),
	glut_entry_CB(NULL)
{
}
