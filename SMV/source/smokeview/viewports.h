// $Date: 2012-10-25 14:45:44 -0400 (Thu, 25 Oct 2012) $ 
// $Revision: 13520 $
// $Author: gforney $

#ifndef VIEWPORTS_H_DEFINED
#define VIEWPORTS_H_DEFINED

void CLIP_viewport(int quad, GLint s_left, GLint s_down);
void INFO_viewport(int quad, GLint s_left, GLint s_down);
void TIMEBAR_viewport(int quad, GLint s_left, GLint s_down);
void COLORBAR_viewport(int quad, GLint s_left, GLint s_down);
void TITLE_viewport(int quad, GLint s_left, GLint s_down);
void Scene_viewport(int quad, int view_mode, GLint s_left, GLint s_down);
#endif
