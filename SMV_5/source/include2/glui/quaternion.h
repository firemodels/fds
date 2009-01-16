/**************************************************************************
  
  quaternion.h

  A quaternion class

  ---------------------------------------------------------------------

  Feb 1998, Paul Rademacher (rademach@cs.unc.edu)  
                
**************************************************************************/

#ifndef _QUATERNION_H_
#define _QUATERNION_H_

#include "algebra3.h"
#include <stdio.h>
#include <stdlib.h>


/* this line defines a new type: pointer to a function which returns a */
/* float and takes as argument a float */
typedef float (*V_FCT_PTR)(float);


/****************************************************************
*                    Quaternion                                 *
****************************************************************/

class quat
{
  /*protected: */
public:

  vec3  v;  /* vector component */
  float s;  /* scalar component */

  /*public: */
  
  /* Constructors */

  quat(void);
  quat(const float x, const float y, const float z, const float w);
  quat( vec3 v, float s ); 
  quat( float s, vec3 v );
  quat(const float *d );      /* copy from four-element float array */
  quat(const double *f );     /* copy from four-element double array */
  quat(const quat &q );        /* copy from other quat */

  /* Assignment operators */

  quat &operator  = ( const quat &v );      /* assignment of a quat */
  quat &operator += ( const quat &v );      /* incrementation by a quat */
  quat &operator -= ( const quat &v );      /* decrementation by a quat */
  quat &operator *= ( const float d );      /* multiplication by a constant */
  quat &operator /= ( const float d );      /* division by a constant */
  float &operator [] ( int i);              /* indexing */
  
  /* special functions */
  
  float length(void);                    /* length of a quat */
  float length2(void);                   /* squared length of a quat */
  quat &normalize(void);                 /* normalize a quat */
  quat &apply(V_FCT_PTR fct);            /* apply a func. to each component */
  void  set( float x, float y, float z );   /* set quat */
  void  set( vec3 v, float s );             /* set quat */
  void  print( FILE *file, char *name );    /* print quat to a file */
  vec3  xform( const vec3 &v );             /* q*v*q-1 */
  mat4  to_mat4( void );
  void  set_angle( float f );               /* set rot angle (degrees) */
  void  scale_angle( float f );             /* scale rot angle (degrees) */
  float get_angle( void );                  /* set rot angle (degrees) */
  vec3  get_axis( void );                   /* get axis */

  /* friends */

  friend quat operator - (const quat &v);                   /* -q1 */
  friend quat operator + (const quat &a, const quat &b);    /* q1 + q2 */
  friend quat operator - (const quat &a, const quat &b);    /* q1 - q2 */
  friend quat operator * (const quat &a, const float d);    /* q1 * 3.0 */
  friend quat operator * (const float d, const quat &a);    /* 3.0 * q1 */
  friend quat operator * (const quat &a, const quat &b);    /* q1 * q2 */
  friend quat operator / (const quat &a, const float d);    /* q1 / 3.0 */
  friend int operator == (const quat &a, const quat &b);    /* q1 == q2 ? */
  friend int operator != (const quat &a, const quat &b);    /* q1 != q2 ? */
  friend void swap(quat &a, quat &b);                       /* swap q1  &q2 */
  /*friend quat min(const quat &a, const quat &b);          -- min(q1, q2) */
  /*friend quat max(const quat &a, const quat &b);          -- max(q1, q2) */
  friend quat prod(const quat &a, const quat &b);      /* term by term mult */
							    }; 

/* Utility functions */

quat quat_identity( void );        /* Returns quaternion identity element */
quat quat_slerp( quat from, quat to, float t );

#endif
