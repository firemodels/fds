// $Date$ 
// $Revision$
// $Author$

#include "options.h"
char shaders_revision[]="$Revision$";

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
#ifdef pp_GPU_BLANK
    "varying float fblank;"
#endif
    "void main(){"
#ifdef pp_GPU_BLANK
     "float a,base;"
     "if(fblank>0.95)discard;"
     "if(fblank>0.001){"
     "  a=1.0-fblank;"
     "  base=1.0-newcolor.a/a;"
     "  if(base>0.0){"
     "    newcolor.a=(1.0-pow(base,a));"
     "  }"
     "  else{"
     "    newcolor.a=0.0;"
     "  }"
     "}"
#endif
      "gl_FragColor = newcolor;"
    "}"
  };

  const GLchar *VertexShaderSource[]={
    "uniform vec4 firecolor;"
    "uniform vec3 norm, eye;"
    "uniform float aspectratio;"
    "uniform float hrrcutoff;"
    "uniform float smoke_shade;"
#ifdef pp_LIGHT
    "attribute float smoke_color;"
#endif
    "uniform float smoke3d_rthick;"
    "uniform int adjustalphaflag;"
    "varying vec4 newcolor;"
#ifdef pp_GPU_BLANK
    "varying float fblank;"
    "attribute float blank;"
#endif
    "attribute float hrr, smoke_alpha;"
    "void main()"
    "{"
    "  float bottom,top,alpha,r;"
    "  float term1, term2, term3, term4;"
#ifdef pp_LIGHT
    "  float s_color;"
#else
    "  float s_shade;"
#endif
    "  vec3 rel_pos;"
    "  if(hrrcutoff>0.0&&hrr>hrrcutoff){"
    "    newcolor = firecolor;"
    "  }"
    "  else{"
    "    alpha=smoke_alpha/256.0;"
    "    if(adjustalphaflag==1||adjustalphaflag==2){"
    "      rel_pos=vec3(gl_Vertex)-eye;"
    "      bottom = abs(dot(rel_pos,norm));"
    "      top=length(rel_pos);"
    "      r=aspectratio*top/bottom;"
    "      term1 = alpha*r;"
    "      term2 = -term1*alpha*(r-1.0)/2.0;"
    "      term3 = -term2*alpha*(r-2.0)/3.0;"
    "      term4 = -term3*alpha*(r-3.0)/4.0;"
    "      alpha = term1+term2+term3+term4;"
    "    }"
    // newcolor.a *= (1.0 - pow(1.0-gl_Color.a,aspectratio*top/bottom));
    "    alpha /= smoke3d_rthick;"
#ifdef pp_LIGHT
    "    s_color = smoke_color/255.0;"
    "    newcolor = vec4(s_color,s_color,s_color,alpha);"
#else
    "    s_shade = smoke_shade/255.0;"
    "    newcolor = vec4(s_shade,s_shade,s_shade,alpha);"
#endif
    "  }"
#ifdef pp_GPU_BLANK
    "  fblank=blank;"
#endif
    "  gl_Position = ftransform();"
    "}"
  };

  v = glCreateShaderObjectARB(GL_VERTEX_SHADER_ARB);
  f = glCreateShaderObjectARB(GL_FRAGMENT_SHADER_ARB);

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
  if(error_code!=1)return error_code;
  GPU_hrrcutoff = glGetUniformLocation(p_smoke,"hrrcutoff");
  GPU_hrr = glGetAttribLocation(p_smoke,"hrr");
#ifdef pp_LIGHT
  GPU_smokecolor = glGetAttribLocation(p_smoke,"smoke_color");
#endif
  GPU_smokealpha = glGetAttribLocation(p_smoke,"smoke_alpha");
  GPU_blank = glGetAttribLocation(p_smoke,"blank");
  GPU_skip = glGetUniformLocation(p_smoke,"skip");
  GPU_smoke3d_rthick = glGetUniformLocation(p_smoke,"smoke3d_rthick");
  GPU_smokeshade = glGetUniformLocation(p_smoke,"smoke_shade");
  GPU_firecolor = glGetUniformLocation(p_smoke,"firecolor");
  GPU_aspectratio = glGetUniformLocation(p_smoke,"aspectratio");
  GPU_norm = glGetUniformLocation(p_smoke,"norm");
  GPU_eye = glGetUniformLocation(p_smoke,"eye");
  GPU_adjustalphaflag = glGetUniformLocation(p_smoke,"adjustalphaflag");
  return error_code;

}

/* ------------------ LoadSmokeShaders ------------------------ */

void LoadSmokeShaders(void){
  glUseProgramObjectARB(p_smoke);
}

/* ------------------ UnloadShadeShaders ------------------------ */

void UnloadSmokeShaders(void){
  glUseProgramObjectARB(0);
}

/* ------------------ init_shaders ------------------------ */

void init_shaders(void){
  char version_label[256];
  int i, major, minor;

	glewInit();
  gpuactive=0;
  strcpy(version_label,(char *)glGetString(GL_VERSION));
  for(i=0;i<strlen(version_label);i++){
    if(version_label[i]=='.')version_label[i]=' ';
  }
  sscanf(version_label,"%i %i",&major,&minor);
  if(major<2){
    trim(version_label);
    printf("Smokeview is running on a system using OpenGL %s\n",version_label);
    printf("OpenGL 2.0 or later is required to use the GPU.\n");
		printf("GPU smoke shader not supported.\n");
    usegpu=0;
    return;
  }

  if (GLEW_ARB_vertex_shader && GLEW_ARB_fragment_shader){
    if(setSmokeShaders()==1){
  		printf("GPU smoke shader successfully compiled, linked and loaded.\n");
      gpuactive=1;
    }
    else{
  		printf("*** GPU smoke shader failed to load.\n");
      usegpu=0;
    }
  }
	else {
		printf("*** GPU smoke shader not supported.\n");
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
        printf("glError in file %s @ line %d: %s\n", file, line, (char *)gluErrorString(glErr));
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

