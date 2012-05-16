// $Date$ 
// $Revision$
// $Author$

#ifndef VIEWPORTS_H_DEFINED
#define VIEWPORTS_H_DEFINED

void CLIP_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height);
void BLOCK_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height);
void TIMEBAR_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height);
void COLORBAR_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height);
void LOGO_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height);
void TITLE_viewport(int quad, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height);
void Scene_viewport(int quad, int view_mode, GLint s_left, GLint s_down, GLsizei s_width, GLsizei s_height);
#endif
