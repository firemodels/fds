// $Date: 2007-10-07 22:08:47 -0400 (Sun, 07 Oct 2007) $ 
// $Revision: 800 $
// $Author: gforney $

#include "options.h"
#ifdef pp_GPU
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <GL/glew.h>
#include "MALLOC.h"
#include "flowfiles.h"
#include "smokeviewdefs.h"
#include "smokeviewvars.h"
#include "smokeheaders.h"

GLhandleARB v,f,f2,p_smoke, p_slice ;

#define printOpenGLError() printOglError(__FILE__, __LINE__)

/* ------------------ printfOglError ------------------------ */

int printOglError(char *file, int line)
{
    //
    // Returns 1 if an OpenGL error occurred, 0 otherwise.
    //
    GLenum glErr;
    int    retCode = 0;

    glErr = glGetError();
    while (glErr != GL_NO_ERROR)
    {
        printf("glError in file %s @ line %d: %s\n", file, line, gluErrorString(glErr));
        retCode = 1;
        glErr = glGetError();
    }
    return retCode;
}

/* ------------------ printfInfoLog ------------------------ */

void printInfoLog(GLhandleARB obj)
{
    int infologLength = 0;
    int charsWritten  = 0;
    char *infoLog;

	glGetObjectParameterivARB(obj, GL_OBJECT_INFO_LOG_LENGTH_ARB,
                                         &infologLength);

    if (infologLength > 0)
    {
        infoLog = (char *)malloc(infologLength);
        glGetInfoLogARB(obj, infologLength, &charsWritten, infoLog);
		printf("%s\n",infoLog);
        free(infoLog);
    }
}

/* ------------------ setSmokeShaders ------------------------ */

int setSmokeShaders() {

  char *vs = NULL,*fs = NULL;
  const char * vv;
  const char * ff;
  GLint error_code;
  char fragbuffer[1024];

  v = glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);
  f = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);
  f2 = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);


//  NewMemory((void **)&vs,NSHADERLINES*80*sizeof(char));
//  NewMemory((void **)&fs,NSHADERLINES*80*sizeof(char));
  
  strcpy(fragbuffer,smokeviewbindir);
  strcat(fragbuffer,"smoke.vert");
  vs = textFileRead(fragbuffer);

  strcpy(fragbuffer,smokeviewbindir);
  strcat(fragbuffer,"smoke.frag");
  fs = textFileRead(fragbuffer);
  
  vv = vs;
  ff = fs;

  // vertex shader

  // pass norm to GPU as uniform vec3
  // pass xyzeyeorig to GPU as uniform vec3
  // pass aspectratio to GPU as uniform float
  // pass skip to GPU as uniform int


  if(vv!=NULL){
    glShaderSourceARB(v, 1, &vv,NULL);
    FREEMEMORY(vs);
    glCompileShaderARB(v);
  }
  if(ff!=NULL){
    glShaderSourceARB(f, 1, &ff,NULL);
    FREEMEMORY(fs);
    glCompileShaderARB(f);
  }

#ifdef _DEBUG
  printInfoLog(v);
  printInfoLog(f);
  printInfoLog(f2);
#endif

  p_smoke = glCreateProgramObjectARB();
  glAttachObjectARB(p_smoke,v);
  glAttachObjectARB(p_smoke,f);

  glLinkProgramARB(p_smoke);
  glGetObjectParameterivARB(p_smoke,GL_OBJECT_LINK_STATUS_ARB,&error_code);
#ifdef _DEBUG
  printf("  Smoke shader completion code:");
  switch (error_code){
  case GL_INVALID_VALUE:
    printf(" INVALID VALUE\n");
    break;
  case GL_INVALID_OPERATION:
    printf(" INVALID OPERATION\n");
    break;
  case GL_INVALID_ENUM:
    printf(" INVALID ENUM\n");
    break;
  case 0:
    printf(" Link failed\n");
    break;
  case 1:
    printf(" Link succeeded\n");
    break;
  default:
    printf(" unknown error\n");
    break;
  }
  printInfoLog(p_smoke);
#endif
  GPU_hrrcutoff = glGetUniformLocationARB(p_smoke,"hrrcutoff");
  GPU_hrr = glGetAttribLocationARB(p_smoke,"hrr");
  GPU_smokealpha = glGetAttribLocationARB(p_smoke,"smoke_alpha");
  GPU_skip = glGetUniformLocationARB(p_smoke,"skip");
  GPU_smoke3d_thick = glGetUniformLocationARB(p_smoke,"smoke3d_thick");
  GPU_smokeshade = glGetUniformLocationARB(p_smoke,"smoke_shade");
  GPU_firealpha = glGetUniformLocationARB(p_smoke,"fire_alpha");
  GPU_firered = glGetUniformLocationARB(p_smoke,"fire_red");
  GPU_firegreen = glGetUniformLocationARB(p_smoke,"fire_green");
  GPU_fireblue = glGetUniformLocationARB(p_smoke,"fire_blue");
  GPU_aspectratio = glGetUniformLocationARB(p_smoke,"aspectratio");
  GPU_normx = glGetUniformLocationARB(p_smoke,"normx");
  GPU_normy = glGetUniformLocationARB(p_smoke,"normy");
  GPU_normz = glGetUniformLocationARB(p_smoke,"normz");
  GPU_eyex = glGetUniformLocationARB(p_smoke,"eyex");
  GPU_eyey = glGetUniformLocationARB(p_smoke,"eyey");
  GPU_eyez = glGetUniformLocationARB(p_smoke,"eyez");
  return error_code;

}

/* ------------------ useSmokeShaders ------------------------ */

void useSmokeShaders(void){
  glUseProgramObjectARB(p_smoke);
}

/* ------------------ useOpenGLShaders ------------------------ */

void useOpenGLShaders(void){
  glUseProgramObjectARB(0);
}

/* ------------------ init_slice_shader ------------------------ */

void init_shaders(void){
	glewInit();
  gpuactive=0;
  if (GLEW_ARB_vertex_shader && GLEW_ARB_fragment_shader){
    if(setSmokeShaders()==1){
  		printf("***GPU smoke shader successfully compiled, linked and loaded.\n");
      gpuactive=1;
    }
    else{
  		printf("***GPU smoke shader failed to load.\n");
    }
  }
	else {
		printf("GPU smoke shader not supported.\n");
    usegpu=0;
	}
}

/* ------------------ textFileRead ------------------------ */

char *textFileRead(char *fn) {


	FILE *fp;
	char *content = NULL;

	int count=0;

	if (fn != NULL) {
		fp = fopen(fn,"rt");

		if (fp != NULL) {
      
      fseek(fp, 0, SEEK_END);
      count = ftell(fp);
      rewind(fp);

			if (count > 0) {
//				content = (char *)malloc(sizeof(char) * (count+1));
        NewMemory((void **)&content,count+1);
				count = fread(content,sizeof(char),count,fp);
				content[count] = '\0';
			}
			fclose(fp);
		}
	}
	return content;
}

/* ------------------ createDepthTexture ------------------------ */

void createDepthTexture( void ){
  if ( GPU_depthtexture!=0 ){
		glDeleteTextures( 1, &GPU_depthtexture );
		GPU_depthtexture = 0;
	}
	
	glGenTextures(1, &GPU_depthtexture);
	sniffErrors("after createDepthTextures 1");
	glBindTexture(GL_TEXTURE_RECTANGLE_EXT, GPU_depthtexture);
	sniffErrors("after createDepthTextures 2");
	glTexParameteri(GL_TEXTURE_RECTANGLE_EXT, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	sniffErrors("after createDepthTextures 3");
	glTexParameteri(GL_TEXTURE_RECTANGLE_EXT, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_RECTANGLE_EXT, 0, GL_DEPTH_COMPONENT, screenWidth, screenHeight,
				  0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);	
					  	  
	sniffErrors("after createDepthTextures 4");
}

/* ------------------ getDepthTexture ------------------------ */

void getDepthTexture( void ){
	if ( GPU_depthtexture==0 ) createDepthTexture();
	glBindTexture(GL_TEXTURE_RECTANGLE_EXT, GPU_depthtexture);
	glCopyTexSubImage2D(GL_TEXTURE_RECTANGLE_EXT, 0, 0,0,
						 0, 0, screenWidth, screenHeight);
	glBindTexture( GL_TEXTURE_RECTANGLE_EXT, 0);	
	sniffErrors("after getDepthTexture");
}

#endif

