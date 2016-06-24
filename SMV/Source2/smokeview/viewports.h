#ifndef VIEWPORTS_H_DEFINED
#define VIEWPORTS_H_DEFINED

void CLIP_viewport(int quad, GLint s_left, GLint s_down);
void INFO_viewport(int quad, GLint s_left, GLint s_down);
void TIMEBAR_viewport(int quad, GLint s_left, GLint s_down);
void COLORBAR_viewport(int quad, GLint s_left, GLint s_down);
void TITLE_viewport(int quad, GLint s_left, GLint s_down);
void Scene_viewport(int quad, int view_mode, GLint s_left, GLint s_down);
#endif
