// $Date$ 
// $Revision$
// $Author$

#include "options.h"
char shaders_revision[]="$Revision$";

#ifdef pp_GPU
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "glew.h"
#include "MALLOC.h"
#include "flowfiles.h"
#include "smokeviewvars.h"
#include "MALLOC.h"

GLhandleARB v,f,p_smoke, p_zonesmoke, p_volsmoke;
char *textFileRead(char *fn);
void printInfoLog(GLhandleARB obj);

/* ------------------ setZoneSmokeShaders ------------------------ */

int setZoneSmokeShaders() {
  GLint error_code;
  const GLchar *FragmentShaderSource[]={
    "uniform int zonedir,zoneinside;"
    "uniform float xyzmaxdiff,zlay,odl,odu;"
    "uniform vec3 eyepos,boxmin, boxmax;"
    "varying vec3 fragpos;"
    
    "void main(){"
    "  float L,opacity,alpha,alpha_min,alpha_zlay;"
    "  float factor_U, factor_L;"
    
    "  vec3 dalphamin,dalphamax;"
    "  L=distance(fragpos,eyepos)*xyzmaxdiff;"
    "  alpha_min=1000000.0;"
    "  dalphamin=-(boxmin-fragpos)/(eyepos-fragpos);"
    "  dalphamax=-(boxmax-fragpos)/(eyepos-fragpos);"
    "  alpha_zlay = -(zlay-fragpos.z)/(eyepos.z-fragpos.z);"
    "  if(zoneinside==0){"
    "    if(zonedir!=-1&&dalphamin.x>0.0&&dalphamin.x<alpha_min)alpha_min=dalphamin.x;"
    "    if(zonedir!=1 &&dalphamax.x>0.0&&dalphamax.x<alpha_min)alpha_min=dalphamax.x;"
    "    if(zonedir!=-2&&dalphamin.y>0.0&&dalphamin.y<alpha_min)alpha_min=dalphamin.y;"
    "    if(zonedir!=2 &&dalphamax.y>0.0&&dalphamax.y<alpha_min)alpha_min=dalphamax.y;"
    "    if(zonedir!=-3&&dalphamin.z>0.0&&dalphamin.z<alpha_min)alpha_min=dalphamin.z;"
    "    if(zonedir!=3 &&dalphamax.z>0.0&&dalphamax.z<alpha_min)alpha_min=dalphamax.z;"
    "    if(eyepos.z>zlay&&fragpos.z>zlay){"
    "      if(alpha_zlay>0.0&&alpha_zlay<alpha_min){"
    "        factor_U=alpha_zlay/odu;"
    "        factor_L=(alpha_min-alpha_zlay)/odl;"
    "      }"
    "      else{"
    "        factor_U=alpha_min/odu;"
    "        factor_L=0.0;"
    "      }"
    "    }"
    "    else if(eyepos.z>zlay&&fragpos.z<=zlay){"
    "      factor_U=0.0;"
    "      factor_L=alpha_min/odl;"
    "    }"
    "    else if(eyepos.z<=zlay&&fragpos.z>zlay){"
    "      factor_U=alpha_min/odu;"
    "      factor_L=0.0;"
    "    }"
    "    else if(eyepos.z<=zlay&&fragpos.z<=zlay){"
    "      if(alpha_zlay>0.0&&alpha_zlay<alpha_min){"
    "        factor_U=(alpha_min-alpha_zlay)/odu;"
    "        factor_L=alpha_zlay/odl;"
    "      }"
    "      else{"
    "        factor_U=0.0;"
    "        factor_L=alpha_min/odl;"
    "      }"
    "    }"
    "  }" // end inside=0
    "  if(zoneinside==1){"
    "    if(eyepos.z>zlay&&fragpos.z>zlay){"
    "      factor_U=1.0/odu;"
    "      factor_L=0.0;"
    "    }"
    "    else if(eyepos.z>zlay&&fragpos.z<=zlay){"
    "      factor_U=(1.0+alpha_zlay)/odu;"
    "      factor_L=-alpha_zlay/odl;"
    "    }"
    "    else if(eyepos.z<=zlay&&fragpos.z>zlay){"
    "      factor_U=-alpha_zlay/odu;"
    "      factor_L=(1.0+alpha_zlay)/odl;"
    "    }"
    "    else if(eyepos.z<=zlay&&fragpos.z<=zlay){"
    "      factor_U=0.0;"
    "      factor_L=1.0/odl;"
    "    }"
    "  }" // end inside=1
    "  opacity = 1.0-exp(-(factor_L+factor_U)*L);"
    "  gl_FragColor = vec4(0.3,0.3,0.3,opacity);"
    "}" // end of main
  };

  const GLchar *VertexShaderSource[]={
    "varying vec3 fragpos;"
    "void main(){"
    "  fragpos=gl_Vertex;"
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

  p_zonesmoke = glCreateProgramObjectARB();
  glAttachObjectARB(p_zonesmoke,v);
  glAttachObjectARB(p_zonesmoke,f);

  glLinkProgram(p_zonesmoke);
  glGetObjectParameterivARB(p_zonesmoke,GL_OBJECT_LINK_STATUS_ARB,&error_code);
#ifdef _DEBUG
  printf("  Zone Smoke shader completion code:");
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
  printInfoLog(p_zonesmoke);
#endif
  GPUzone_zoneinside = glGetUniformLocation(p_zonesmoke,"zoneinside");
  GPUzone_zonedir = glGetUniformLocation(p_zonesmoke,"zonedir");
  GPUzone_eyepos = glGetUniformLocation(p_zonesmoke,"eyepos");
  GPUzone_xyzmaxdiff = glGetUniformLocation(p_zonesmoke,"xyzmaxdiff");
  GPUzone_boxmin = glGetUniformLocation(p_zonesmoke,"boxmin");
  GPUzone_boxmax = glGetUniformLocation(p_zonesmoke,"boxmax");
  GPUzone_zlay = glGetUniformLocation(p_zonesmoke,"zlay");
  GPUzone_odl = glGetUniformLocation(p_zonesmoke,"odl");
  GPUzone_odu = glGetUniformLocation(p_zonesmoke,"odu");

  if(error_code!=1)return error_code;
  return error_code;

}

/* ------------------ setVolSmokeShaders ------------------------ */

int setVolSmokeShaders() {
  GLint error_code;
  const GLchar *FragmentShaderSource[]={
    "uniform sampler1D smokecolormap;"
    "uniform sampler3D soot_density_texture,fire_texture;"
    "uniform int dir,inside;"
    "uniform float xyzmaxdiff,dcell;"
    "uniform vec3 eyepos,boxmin, boxmax;"
    "varying vec3 fragpos;"
    
    "void main(){"
    "  vec3 dalphamin,dalphamax,fragmaxpos,posi;"
    "  vec3 pt_soot, pt_color,cum_color;"
    "  float opacity,alpha_min,factor,pathdist,k;"
    "  float colorindex;"
    "  float tauhat, alphahat, taui, tauterm;"
    "  float dstep;"
    "  int i,n_iter;"

    "  alpha_min=1000000.0;"
    "  dalphamin=-(boxmin-fragpos)/(eyepos-fragpos);"
    "  dalphamax=-(boxmax-fragpos)/(eyepos-fragpos);"
    "  if(inside==0){"
    "    if(dir!=-1&&dalphamin.x>0.0&&dalphamin.x<alpha_min)alpha_min=dalphamin.x;"
    "    if(dir!=1 &&dalphamax.x>0.0&&dalphamax.x<alpha_min)alpha_min=dalphamax.x;"
    "    if(dir!=-2&&dalphamin.y>0.0&&dalphamin.y<alpha_min)alpha_min=dalphamin.y;"
    "    if(dir!=2 &&dalphamax.y>0.0&&dalphamax.y<alpha_min)alpha_min=dalphamax.y;"
    "    if(dir!=-3&&dalphamin.z>0.0&&dalphamin.z<alpha_min)alpha_min=dalphamin.z;"
    "    if(dir!=3 &&dalphamax.z>0.0&&dalphamax.z<alpha_min)alpha_min=dalphamax.z;"
    "  }" // end inside=0
    "  if(inside==1){"
    "  }" // end inside=1
    "  fragmaxpos = mix(fragpos,eyepos,-alpha_min);"
    "  pathdist = distance(fragpos,fragmaxpos);"
    "  n_iter = max(1,int(pathdist/dcell+0.5));"
    "  dstep = pathdist*xyzmaxdiff/n_iter;"
    "  tauhat=1.0;"
    "  alphahat=0.0;"
    "  cum_color=vec3(0.0,0.0,0.0);"
    "  k=8700;"
    "  for(i=0;i<n_iter;i++){"
    "    factor = (0.5+i)/n_iter;"
    "    posi = (mix(fragpos,fragmaxpos,factor)-boxmin)/(boxmax-boxmin);"
    "    colorindex = texture3D(fire_texture,posi)/1200.0;" 
    "    pt_color = texture1D(smokecolormap,colorindex).rgb;"
    "    pt_soot = texture3D(soot_density_texture,posi);"
    "    if(colorindex>0.5){"
    "      pt_soot *= 3.0;"
    "    };"
    "    taui = exp(-k*pt_soot*dstep);"
    "    tauterm = (1.0-taui)*tauhat;"
    "    alphahat  += tauterm;"
    "    cum_color += tauterm*pt_color;"
    "    tauhat *= taui;"
    "  }"
    "  cum_color /= alphahat;"
    "  gl_FragColor = vec4(cum_color,alphahat);"
    "}" // end of main
  };

  const GLchar *VertexShaderSource[]={
    "varying vec3 fragpos;"
    "void main(){"
    "  fragpos=gl_Vertex;"
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

  p_volsmoke = glCreateProgramObjectARB();
  glAttachObjectARB(p_volsmoke,v);
  glAttachObjectARB(p_volsmoke,f);

  glLinkProgram(p_volsmoke);
  glGetObjectParameterivARB(p_volsmoke,GL_OBJECT_LINK_STATUS_ARB,&error_code);
#ifdef _DEBUG
  printf("  Zone Smoke shader completion code:");
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
  printInfoLog(p_volsmoke);
#endif
  GPUvol_inside = glGetUniformLocation(p_volsmoke,"inside");
  GPUvol_dir    = glGetUniformLocation(p_volsmoke,"dir");
  GPUvol_eyepos = glGetUniformLocation(p_volsmoke,"eyepos");
  GPUvol_dcell = glGetUniformLocation(p_volsmoke,"dcell");
  GPUvol_xyzmaxdiff = glGetUniformLocation(p_volsmoke,"xyzmaxdiff");
  GPUvol_boxmin = glGetUniformLocation(p_volsmoke,"boxmin");
  GPUvol_boxmax = glGetUniformLocation(p_volsmoke,"boxmax");
  GPUvol_soot_density = glGetUniformLocation(p_volsmoke,"soot_density_texture");
  GPUvol_fire = glGetUniformLocation(p_volsmoke,"fire_texture");
  GPUvol_smokecolormap = glGetUniformLocation(p_volsmoke,"smokecolormap");

  if(error_code!=1)return error_code;
  return error_code;

}

/* ------------------ setSmokeShaders ------------------------ */

int setSmokeShaders() {
  GLint error_code;

  const GLchar *FragmentShaderSource[]={
    "varying vec4 newcolor;"
    "void main(){"
    "  gl_FragColor = newcolor;"
    "}"
  };

  const GLchar *VertexShaderSource[]={
    "uniform sampler1D smokecolormap;"
    "uniform float hrrpuv_max_smv, hrrpuv_cutoff;"
    "uniform float aspectratio, smoke3d_rthick, fire_alpha;"
    "uniform int is_smoke, adjustalphaflag;"

    "attribute float hrr, smoke_alpha;"

    "varying vec4 newcolor;"
    "void main()"
    "{"
    "  float alpha,r;"
    "  float term1, term2, term3, term4;"
    "  vec4 hrrcolor;"
    "  float colorindex;"
    "  float hrrlocal;"

    "  alpha=0.0;"
    "  if(is_smoke==1){"
    "    alpha=smoke_alpha/255.0;"
//    f(alpha) = 1 - (1-alpha)^r = f(0) + f'(0)alpha + f''(0)alpha^2/2 + f'''(0)alpha^3/6 + ...
//             f(0) = 0, f'(0) = r, f''(0) = -r(r-1), f'''(0) = r(r-1)(r-2), f''''(0)=-r(r-1)(r-2)(r-3)
    "    if(adjustalphaflag==1||adjustalphaflag==3){"
    "      r=aspectratio;"
    "      term1 = alpha*r;"
    "      term2 = -term1*alpha*(r-1.0)/2.0;"
    "      term3 = -term2*alpha*(r-2.0)/3.0;"
    "      term4 = -term3*alpha*(r-3.0)/4.0;"
    "      alpha = term1+term2+term3+term4;"
    "    }"
    "    alpha /= smoke3d_rthick;"
    "  }"
    "  hrrlocal=(hrr/254.0)*hrrpuv_max_smv;"
    "  if(hrrlocal>hrrpuv_cutoff){"
    "    colorindex=0.51+(hrrlocal-hrrpuv_cutoff)/(hrrpuv_max_smv-hrrpuv_cutoff);"
    "    colorindex=clamp(colorindex,0.5,1.0);"
    "    alpha*=3.0;"
    "  }"
    "  else{"
    "    colorindex=hrrlocal/hrrpuv_cutoff;"
    "    colorindex=clamp(colorindex,0.0,0.49);"
    "  }"
    "  hrrcolor = texture1D(smokecolormap,colorindex);"
    "  newcolor=vec4(vec3(hrrcolor),alpha);"
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
  GPU_hrrpuv_max_smv = glGetUniformLocation(p_smoke,"hrrpuv_max_smv");
  GPU_hrrpuv_cutoff =  glGetUniformLocation(p_smoke,"hrrpuv_cutoff");
  GPU_smokecolormap =  glGetUniformLocation(p_smoke,"smokecolormap");
  GPU_smoke3d_rthick = glGetUniformLocation(p_smoke,"smoke3d_rthick");
  GPU_is_smoke =       glGetUniformLocation(p_smoke,"is_smoke");
  GPU_aspectratio =    glGetUniformLocation(p_smoke,"aspectratio");
  GPU_adjustalphaflag =glGetUniformLocation(p_smoke,"adjustalphaflag");
  GPU_hrr =            glGetAttribLocation(p_smoke,"hrr");
  GPU_smokealpha =     glGetAttribLocation(p_smoke,"smoke_alpha");
  return error_code;

}

/* ------------------ LoadZoneSmokeShaders ------------------------ */

void LoadZoneSmokeShaders(void){
  glUseProgramObjectARB(p_zonesmoke);
}

/* ------------------ LoadSmokeShaders ------------------------ */

void LoadSmokeShaders(void){
  glUseProgramObjectARB(p_smoke);
}

/* ------------------ LoadVolSmokeShaders ------------------------ */

void LoadVolSmokeShaders(void){
  glUseProgramObjectARB(p_volsmoke);
}

/* ------------------ UnloadShaders ------------------------ */

void UnLoadShaders(void){
  glUseProgramObjectARB(0);
}

/* ------------------ init_shaders ------------------------ */

int init_shaders(void){
  char version_label[256];
  char version_label2[256];
  int i, major, minor;
  GLenum err;
  const GLubyte *version_string;

  gpuactive=0;
  usegpu=0;
  err=glewInit();
  if(err!=GLEW_OK){
    printf("   *** Warning: Initialization of GLEW extenstion library failed\n");
    return 1;
  }
  version_string=glGetString(GL_VERSION);
  if(version_string==NULL){
    printf("   *** warning: GL_VERSION string is NULL in init_shaders()\n");
    err = 1;
    return err;
  }
  strcpy(version_label,(char *)version_string);
  strcpy(version_label2,version_label);
  for(i=0;i<strlen(version_label);i++){
    if(version_label[i]=='.')version_label[i]=' ';
  }
  sscanf(version_label,"%i %i",&major,&minor);
  if(major<2){
    trim(version_label);
    printf("   Smokeview is running on a system using OpenGL %s\n",version_label2);
    printf("   OpenGL 2.0 or later is required to use the GPU.\n");
    printf("   GPU smoke shader not supported.\n");
    return 1;
  }

  if (GLEW_ARB_vertex_shader && GLEW_ARB_fragment_shader){
    if(setSmokeShaders()==1){
#ifdef _DEBUG
  		printf("   GPU smoke shader successfully compiled, linked and loaded.\n");
#endif
      gpuactive=1;
      err=0;
    }
    else{
      printf("   *** GPU smoke shader failed to load.\n");
      usegpu=0;
      err=1;
    }
    if(setVolSmokeShaders()==1){
#ifdef _DEBUG
  		printf("   GPU smoke Volume shader successfully compiled, linked and loaded.\n");
#endif
      gpuactive=1;
      err=0;
    }
    else{
      printf("   *** GPU smoke volume shader failed to load.\n");
      usegpu=0;
      err=1;
    }
    if(err==0){
      if(setZoneSmokeShaders()==1){
#ifdef _DEBUG
  		printf("   GPU zone smoke shader successfully compiled, linked and loaded.\n");
#endif
        gpuactive=1;
        err=0;
      }
      else{
        printf("   *** GPU zone smoke shader failed to load.\n");
        usegpu=0;
        err=1;
      }
    }
  }
	else {
    printf("   *** GPU smoke shader not supported.\n");
    usegpu=0;
    err=1;
	}
  return err;
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

void printInfoLog(GLhandleARB obj){
  int infologLength = 0;
  int charsWritten  = 0;
  char *infoLog;

	glGetObjectParameterivARB(obj, GL_OBJECT_INFO_LOG_LENGTH_ARB,&infologLength);
  if (infologLength > 0){
    NewMemory((void **)&infoLog,infologLength);
    glGetInfoLogARB(obj, infologLength, &charsWritten, infoLog);
    printf("%s\n",infoLog);
    FREEMEMORY(infoLog);
  }
}
#endif

