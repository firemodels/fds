/****************************************************************************

  quaternion.h - A quaternion class

  GLUI User Interface Toolkit 
  Copyright (c) 1998 Paul Rademacher

  ---------------------------------------------------------------------

  WWW:    http://sourceforge.net/projects/glui/
  Forums: http://sourceforge.net/forum/?group_id=92496

  This software is provided 'as-is', without any express or implied 
  warranty. In no event will the authors be held liable for any damages 
  arising from the use of this software. 

  Permission is granted to anyone to use this software for any purpose, 
  including commercial applications, and to alter it and redistribute it 
  freely, subject to the following restrictions: 

  1. The origin of this software must not be misrepresented; you must not 
  claim that you wrote the original software. If you use this software 
  in a product, an acknowledgment in the product documentation would be 
  appreciated but is not required. 
  2. Altered source versions must be plainly marked as such, and must not be 
  misrepresented as being the original software. 
  3. This notice may not be removed or altered from any source distribution. 

*****************************************************************************/

#ifndef GLUI_QUATERNION_H
#define GLUI_QUATERNION_H

#include "algebra3.h"
#include <cstdio>

/* this line defines a new type: pointer to a function which returns a */
/* float and takes as argument a float */
typedef float (*V_FCT_PTR)(float);

/****************************************************************
 *                    Quaternion                                *
 ****************************************************************/

class quat
{
  /*protected: */
public:

  vec3  v;  /* vector component */
  float s;  /* scalar component */

  /*public: */
  
  /* Constructors */

  quat();
  quat(float x, float y, float z, float w);
  quat(const vec3 &v, float s); 
  quat(float   s, const vec3 &v);
  quat(const float  *d);     /* copy from four-element float array  */
  quat(const double *f);     /* copy from four-element double array */
  quat(const quat   &q);     /* copy from other quat                */

  /* Assignment operators */

  quat  &operator  = (const quat &v);      /* assignment of a quat            */
  quat  &operator += (const quat &v);      /* incrementation by a quat        */
  quat  &operator -= (const quat &v);      /* decrementation by a quat        */
  quat  &operator *= (float d);      /* multiplication by a constant    */
  quat  &operator /= (float d);      /* division by a constant          */
  
  /* special functions */
  
  float  length() const;                   /* length of a quat                */
  float  length2() const;                  /* squared length of a quat        */
  quat  &normalize();                      /* normalize a quat                */
  quat  &apply(V_FCT_PTR fct);             /* apply a func. to each component */
  vec3   xform(const vec3 &v );            /* q*v*q-1                         */
  mat4   to_mat4() const;
  void   set_angle(float f);               /* set rot angle (degrees)         */
  void   scale_angle(float f);             /* scale rot angle (degrees)       */
  float  get_angle() const;                /* set rot angle (degrees)         */
  vec3   get_axis()  const;                /* get axis                        */

  void   print( FILE *file, const char *name ) const;  /* print to a file     */

        float &operator [] (int i);        /* indexing                        */
  const float &operator [] (int i) const;  /* indexing                        */

  void   set(float x, float y, float z);   /* set quat                        */
  void   set(const vec3 &v, float s);      /* set quat                        */

  /* friends */

  friend quat operator - (const quat &v);                   /* -q1            */
  friend quat operator + (const quat &a, const quat &b);    /* q1 + q2        */
  friend quat operator - (const quat &a, const quat &b);    /* q1 - q2        */
  friend quat operator * (const quat &a, float d);          /* q1 * 3.0       */
  friend quat operator * (float d, const quat &a);          /* 3.0 * q1       */
  friend quat operator * (const quat &a, const quat &b);    /* q1 * q2        */
  friend quat operator / (const quat &a, float d);          /* q1 / 3.0       */
  friend int operator == (const quat &a, const quat &b);    /* q1 == q2 ?     */
  friend int operator != (const quat &a, const quat &b);    /* q1 != q2 ?     */
  friend void swap(quat &a, quat &b);                       /* swap q1  &q2   */
  /*friend quat min(const quat &a, const quat &b);          -- min(q1, q2)    */
  /*friend quat max(const quat &a, const quat &b);          -- max(q1, q2)    */
  friend quat prod(const quat &a, const quat &b);          /* term by term mult*/
}; 

/* Utility functions */

quat quat_identity();        /* Returns quaternion identity element */
quat quat_slerp(const quat &from, const quat &to, float t);

#endif
