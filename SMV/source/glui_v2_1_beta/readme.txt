Welcome to the GLUI User Interface Library, v2.0 beta!
-------------------------------------------------

This distribution contains the full GLUI sources, as well as 5 example
programs.  You'll find the full manual under "glui_manual.pdf".  The
GLUI web page is at 

	http://www.cs.unc.edu/~rademach/glui


		    ---------- Windows ----------

The directory 'msvc' contains a Visual C++ workspace entitled
'glui.dsw'.  To recompile the library and examples, open this
workspace and run the menu command "Build:Batch Build:Build".  The 3
executables will be in the 'bin' directory, and the library in the
'lib' directory.

To create a new Windows executable using GLUI, create a "Win32 Console
Application" in VC++, add the GLUI library (in 'msvc/lib/glui32.lib'),
and add the OpenGL libs:

	glui32.lib glut32.lib glu32.lib opengl32.lib   (Microsoft OpenGL)

Include the file "glui.h" in any file that uses the GLUI library.


		      ---------- Unix ----------

An SGI/HP makefile is found in the file 'makefile' (certain lines may need 
to be commented/uncommented).

To include GLUI in your own apps, add the glui library to your
makefile (before the glut library 'libglut.a'), and include "glui.h"
in your sources.



----------------------------------------------------------------------

Please let me know what you think, what you'd like to change or add,
and especially what bugs you encounter.  Also, please send me your
e-mail so I can add you to a mailing list for updates.

Good luck, and thanks for trying this out!

Paul Rademacher
rademach@cs.unc.edu