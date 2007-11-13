// $Date$ 
// $Revision$
// $Author$

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

GLhandleARB v,f,p_smoke, p_slice ;
char *textFileRead(char *fn);
void printInfoLog(GLhandleARB obj);

/* ------------------ setSmokeShaders ------------------------ */

int setSmokeShaders() {
  GLint error_code;

  const GLchar *FragmentShaderSource[]={
    "varying vec4 newcolor;"
    "uniform sampler2D depthmap;"
    "uniform float nearPlane;"
    "uniform float farPlane;"

//float convertZ( in float near, in float far, in float depthBufferValue )
//{
//	float clipZ = ( depthBufferValue - 0.5 ) * 2.0;
//	return -(2.0 * far * near) / ( clipZ * ( far - near ) - ( far + near ));
//}

    "void main(){"

//   float sceneDepth = texture2DRect( depthmap, gl_FragCoord.st ).z;
//   float sceneDepth = shadow2D( depthmap, vec3(gl_FragCoord) ).z;
//   float thickness = convertZ( nearPlane, farPlane, sceneDepth ) - gl_FragCoord.z;

// if(gl_FragCoord.z<=0.75)discard;
     "gl_FragColor = newcolor;"
    "}"
  };

  const GLchar *VertexShaderSource[]={
    "uniform float aspectratio,normx,normy,normz;"
    "uniform float eyex,eyey,eyez;"
    "uniform float fire_red, fire_green, fire_blue, fire_alpha;"
    "uniform int skip;"
    "varying vec4 newcolor;"
    "uniform float hrrcutoff;"
    "uniform float smoke_shade;"
    "uniform float smoke3d_rthick;"
    "attribute float hrr, smoke_alpha;"
    "void main()"
    "{"
    "  float bottom,top,alpha,r;"
    "  float term1, term2, term3, term4;"
    "  vec3 rel_pos,eyexyz;"
    "  if(hrrcutoff>0.0&&hrr>hrrcutoff){"
    "    newcolor = vec4(fire_red,fire_green,fire_blue,fire_alpha);"
    "  }"
    "  else{"
    "    eyexyz = vec3(eyex,eyey,eyez);"
    "    rel_pos=vec3(gl_Vertex)-eyexyz;"
    "    bottom = abs(rel_pos.x*normx+rel_pos.y*normy+rel_pos.z*normz);"
    "    top=length(rel_pos);"
    "    r=aspectratio*top/bottom;"
    "    alpha=smoke_alpha/256.0;"
    "    term1 = alpha*r;"
    "    term2 = -term1*alpha*(r-1.0)/2.0;"
    "    term3 = -term2*alpha*(r-2.0)/3.0;"
    "    term4 = -term3*alpha*(r-3.0)/4.0;"
    "    alpha = term1+term2+term3+term4;"
    // newcolor.a *= (1.0 - pow(1.0-gl_Color.a,aspectratio*top/bottom));
    "    if(skip==2){"
    "      alpha = 2.0*alpha*(1.0-alpha);"
    "    }"
    "    else if(skip==3){"
    "      alpha = 3.0*alpha*(1.0-alpha-alpha*alpha/3.0);"
    "    }"
    "    alpha /= smoke3d_rthick;"
    "    newcolor = vec4(smoke_shade,smoke_shade,smoke_shade,alpha);"
    "  }"
    "  gl_Position = ftransform();"
    "}"
  };

  v = glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);
  f = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);

  // vertex shader

  // pass norm to GPU as uniform vec3
  // pass xyzeyeorig to GPU as uniform vec3
  // pass aspectratio to GPU as uniform float
  // pass skip to GPU as uniform int


    glShaderSource(v,1, VertexShaderSource,NULL);
    glCompileShaderARB(v);

    glShaderSource(f, 1, FragmentShaderSource,NULL);
    glCompileShaderARB(f);

#ifdef _DEBUG
  printInfoLog(v);
  printInfoLog(f);
#endif

  p_smoke = glCreateProgramObjectARB();
  glAttachObjectARB(p_smoke,v);
  glAttachObjectARB(p_smoke,f);

  glLinkProgram(p_smoke);
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
  GPU_hrrcutoff = glGetUniformLocation(p_smoke,"hrrcutoff");
  GPU_hrr = glGetAttribLocation(p_smoke,"hrr");
  GPU_smokealpha = glGetAttribLocation(p_smoke,"smoke_alpha");
  GPU_skip = glGetUniformLocation(p_smoke,"skip");
  GPU_smoke3d_rthick = glGetUniformLocation(p_smoke,"smoke3d_rthick");
  GPU_smokeshade = glGetUniformLocation(p_smoke,"smoke_shade");
  GPU_firealpha = glGetUniformLocation(p_smoke,"fire_alpha");
  GPU_firered = glGetUniformLocation(p_smoke,"fire_red");
  GPU_firegreen = glGetUniformLocation(p_smoke,"fire_green");
  GPU_fireblue = glGetUniformLocation(p_smoke,"fire_blue");
  GPU_aspectratio = glGetUniformLocation(p_smoke,"aspectratio");
  GPU_normx = glGetUniformLocation(p_smoke,"normx");
  GPU_normy = glGetUniformLocation(p_smoke,"normy");
  GPU_normz = glGetUniformLocation(p_smoke,"normz");
  GPU_eyex = glGetUniformLocation(p_smoke,"eyex");
  GPU_eyey = glGetUniformLocation(p_smoke,"eyey");
  GPU_eyez = glGetUniformLocation(p_smoke,"eyez");
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
#endif

