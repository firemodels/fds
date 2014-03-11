/****************************************************************************
  
  GLUI User Interface Toolkit
  ---------------------------

     glui.cpp


          --------------------------------------------------

  Copyright (c) 1998 Paul Rademacher (this file, Bill Baxter 2005)

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

#include "GL/glui.h"
#include <stdarg.h>

#ifdef _MSC_VER
#define vsnprintf _vsnprintf
#endif

GLUI_String& glui_format_str(GLUI_String& str, const char* fmt, ...)
{
  const size_t ISIZE = 128;
  char stackbuf[ISIZE];
  size_t bufsz = ISIZE;
  char *buf = stackbuf;
  str = "";
  va_list arg;
  while (1) {
    va_start(arg, fmt);
    int ret = vsnprintf(buf,bufsz-1,fmt,arg);
    va_end(arg);
    if (ret>=0) {
      break;
    }
    // else make a bigger buf, try again
    bufsz <<= 1;
    if (buf==stackbuf) buf = (char*)malloc(sizeof(char)*bufsz);
    else buf = (char*)realloc(buf, sizeof(char)*bufsz);
  }
  if (buf!=stackbuf) free(buf);
  str=buf;
  return str;
}
