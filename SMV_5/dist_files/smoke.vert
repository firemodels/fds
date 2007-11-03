// $Date$ 
// $Revision$
// $Author$
uniform float aspectratio,normx,normy,normz;
uniform float eyex,eyey,eyez;
uniform float fire_red, fire_green, fire_blue, fire_alpha;
uniform int skip;
varying vec4 newcolor;
uniform float hrrcutoff;
uniform float smoke_shade;
uniform int smoke3d_thick;
attribute float hrr, smoke_alpha;

void main(){
  float bottom,top,alpha,r;
  float term1, term2, term3, term4;
  vec3 rel_pos,eyexyz;
  int drawfire=0;

  if(hrrcutoff>0.0&&hrr>hrrcutoff){
    newcolor = vec4(fire_red,fire_green,fire_blue,fire_alpha);
  }
  else{
    eyexyz = vec3(eyex,eyey,eyez);
    rel_pos=vec3(gl_Vertex)-eyexyz;
    bottom = abs(rel_pos.x*normx+rel_pos.y*normy+rel_pos.z*normz);
    top=length(rel_pos);
    r=aspectratio*top/bottom;
    alpha=smoke_alpha/256.0;
    term1 = alpha*r;
    term2 = -term1*alpha*(r-1.0)/2.0;
    term3 = -term2*alpha*(r-2.0)/3.0;
    term4 = -term3*alpha*(r-3.0)/4.0;
    alpha = term1+term2+term3+term4;
// newcolor.a *= (1.0 - pow(1.0-gl_Color.a,aspectratio*top/bottom));
    if(skip==2){
      alpha = 2.0*alpha*(1.0-alpha);
    }
    else if(skip==3){
      alpha = 3.0*alpha*(1.0-alpha-alpha*alpha/3.0);
    }
    if(smoke3d_thick==1){
      alpha=alpha/2;
      }
    else if (smoke3d_thick==2){
      alpha=alpha/4;
    }
    else if (smoke3d_thick==3){
      alpha=alpha/8;
    }
    else if (smoke3d_thick==4){
      alpha=alpha/16;
    }
    else if (smoke3d_thick==5){
      alpha=alpha/32;
    }
    else if (smoke3d_thick==6){
      alpha=alpha/64;
    }
    else if (smoke3d_thick==7){
      alpha=alpha/128;
    }
    newcolor = vec4(smoke_shade,smoke_shade,smoke_shade,alpha);
  }
  gl_Position = ftransform();
}

