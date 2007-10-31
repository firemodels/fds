varying vec4 newcolor;
uniform sampler2D depthmap;
uniform float nearPlane;
uniform float farPlane;

float convertZ( in float near, in float far, in float depthBufferValue )
{
	float clipZ = ( depthBufferValue - 0.5 ) * 2.0;
	return -(2.0 * far * near) / ( clipZ * ( far - near ) - ( far + near ));
}

void main(){

//   float sceneDepth = texture2DRect( depthmap, gl_FragCoord.st ).z;
//   float sceneDepth = shadow2D( depthmap, vec3(gl_FragCoord) ).z;
//   float thickness = convertZ( nearPlane, farPlane, sceneDepth ) - gl_FragCoord.z;

 // if(gl_FragCoord.z<=0.75)discard;
  gl_FragColor = newcolor;
}
