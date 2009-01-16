#ifndef _STDINC_H_
#define _STDINC_H_

#include <stdio.h>
#include <math.h>

#ifndef AND
#define AND &&
#define OR  ||
#define NOT !
#endif

#ifndef true
#define true 1
#define false 0
#endif

#ifndef Bool
typedef char Bool;
#endif

#ifndef MAX
#define MAX(a,b)  ((a)>(b) ? (a) : (b))
#define MIN(a,b)  ((a)<(b) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(a) ((a)>=0 ? (a) : (-(a)))
#endif

/********* TOGGLE_BOOL(boolean) : toggles values of 'boolean' ******/
#ifndef TOGGLE_BOOL
#define TOGGLE_BOOL(a)  ((a)=1-(a))
#endif


/********************  bit comparisons and operations ***************/
#ifndef TEST_BIT
#define TEST_BIT( x, b )   (((x) & (1<<(b))) != 0 )
#define SET_BIT( x, b )    ((x) |= (1 << (b)))
#define CLEAR_BIT( x, b )  ((x) &= ~(1 << (b)))
#define TOGGLE_BIT( x, b ) ((TEST_BIT(x,b)) ?(CLEAR_BIT(x,b)):(SET_BIT(x,b)))
#endif

#ifndef TEST_AND
#define TEST_AND( a, b ) ((a&b)==b)
#endif


#ifndef M_PI
#define M_PI 3.141592654
#endif

typedef  char String[81];


/*********** flush the stdout and stderr output streams *************/
#ifndef flushout
#define flushout fflush(stdout)
#define flusherr fflush(stderr)
#endif


/********** Debugging functions *************************************/
#ifndef error_return
#define error_return( c ); {fprintf(stderr,c);return;}
#endif


/************************* floating-point random ********************/
#ifndef randf
#define randf()  ((float) rand() / (float)RAND_MAX )
#endif

#ifndef SIGN
#define SIGN(x) ((x)>=0 ?  1  :  -1)
#endif


/****************** conversion between degrees and radians **********/
#ifndef DEG2RAD
#define DEG2RAD(x) ((x)/180.0*M_PI)
#define RAD2DEG(x) ((x)/M_PI*180.0)
#endif


/***************** clamp a value to some fixed interval *************/
#ifndef CLAMP
#define CLAMP(x,lo,hi)  {if ((x) < (lo)) {(x)=(lo);} else if((x) > (hi)) {(x)=(hi);}}
#endif



/************ check if a value lies within a closed interval *********/
#ifndef IN_BOUNDS
#define IN_BOUNDS( x, lo, hi ) ( (x) >= (lo) AND (x) <= (hi) )
#endif


/************ check if a 2D point lies within a 2D box ***************/
#ifndef PT_IN_BOX
#define PT_IN_BOX( x, y, lo_x, hi_x, lo_y, hi_y ) \
( IN_BOUNDS(x,lo_x,hi_x) AND IN_BOUNDS(y,lo_y,hi_y) )
#endif

     /****** check if value lies on proper side of another value     *****/
     /*** if side is positive => proper side is positive, else negative **/
#ifndef CHECK_PROPER_SIDE
#define CHECK_PROPER_SIDE(x,val,side) ((side) > 0 ? (x) > (val) : (x) < (val))
#endif


     /***** Small value when we want to do a comparison to 'close to zero' *****/
#ifndef FUDGE
#define FUDGE .00001
#endif


     /******************* swap two values, using a temp variable *********/
#ifndef SWAP2
#define SWAP2(a,b,t) {t=a;a=b;b=t;}     
#endif

#define VEC3_TO_ARRAY(v,a)  a[0]=v[0], a[1]=v[1], a[2]=v[2]

#endif /* _STDINC_H_ */




     
