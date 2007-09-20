// #define pp_LINT // uncomment if using LINT
#define _CRT_SECURE_NO_DEPRECATE
#ifdef WIN32
#define pp_cvf
#define pp_memstatus
#define pp_COMPRESS
#define pp_noappend

#pragma warning (disable:4100)		/* disable bogus conversion warnings */
#pragma warning (disable:4115)		/* disable bogus conversion warnings */
#pragma warning (disable:4127)		/* disable bogus conversion warnings */
#pragma warning (disable:4201)		/* disable bogus conversion warnings */
#pragma warning (disable:4244)		/* disable bogus conversion warnings */
#pragma warning (disable:4305)		/* disable bogus conversion warnings */
#pragma warning (disable:4505)		/* disable bogus conversion warnings */
#pragma warning (disable:4701)		/* disable bogus conversion warnings */
#pragma warning (disable:4018)		/* signied/unsigned match */
#pragma warning (disable:4206)		/* translation unit empty */
#pragma warning (disable:4267)		/* size_t to int possible loss of data */
#pragma warning (disable:4389)		/* signed/unsigned mis-match */


#endif

#ifndef pp_OSX
//#define pp_GDGIF
#endif

// #define pp_FINALxxxRELEASE
#ifndef _DEBUG
#define pp_RELEASE
#endif

#define pp_glui
#define pp_DRAWISO
#define EGZ
#define USE_ZLIB
#define pp_THREADS2
#define pp_HRR

#define pp_WUI
//#define pp_BETA  // uncomment when building beta versions of SV

#ifdef pp_TEST2
#define pp_TEST
#define pp_SVNET
#define pp_SPOTLIGHT
#endif

#ifdef pp_TEST
#define pp_TOUR
#define pp_COLOR
#define pp_LIGHT
#define pp_GPU
#define pp_MEM2
#endif

#undef pp_OPEN
#ifdef WIN32
#ifdef _DEBUG
#define pp_OPEN
#endif
#ifdef pp_MEM2
#define pp_OPEN
#endif
#endif


#ifdef pp_LINT
/*lint -e534 */
/*lint -e553 */
/*lint -e506 */
/*lint -e578 */  //hide symbol
/*lint -e537 */
/*lint -e732 loss of sign */
/*lint -e774 boolean within if always evaluates to True */
/*lint -e785 */
/*lint -e736 */
/*lint -e818 */
/*lint -e790 */
/*lint -e747 */
/*lint -e737 */
/*lint -e524 */
/*lint -e834 */
/*lint -e725 */
/*lint -e539 */
/*lint -e613 */
/*lint -e525 */
/*lint -e19 Expecting ; */
/*lint -e10 Useless declaration */
/*lint -e746 */
/*lint -e644 */
/*lint -e529 */
/*lint -e795 */
/*lint -e526 */
/*lint -e628 */
/*lint -e777 */
/*lint -e1065 */
/*lint -e1776 */

/*lint -e701 */
/*lint -e713 */
/*lint -e734 */
/*lint -e765 */
/*lint -e773 */

#endif








