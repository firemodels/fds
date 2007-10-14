// $Date$ 
// $Revision$
// $Author$

#ifdef pp_BETA
#define SMVVERSION "Test"
#else
#ifdef pp_ALPHA
#define SMVVERSION "Experimental"
#else
#define SMVVERSION "5.0.2"
#endif
#endif

// #define pp_LINT // uncomment if using LINT
#define _CRT_SECURE_NO_DEPRECATE
#ifdef WIN32
#define pp_cvf
#define pp_memstatus
#define pp_COMPRESS
#define pp_noappend

#ifdef pp_WIN_INTEL
#pragma warning (disable:810)		/* conversion from xx to yy may lose sig bits */
#pragma warning (disable:869)		/* parameer was never referenced */
#pragma warning (disable:4018)		/* signied/unsigned match */
#pragma warning (disable:4206)		/* translation unit empty */
#pragma warning (disable:4267)		/* size_t to int possible loss of data */
#pragma warning (disable:4389)		/* signed/unsigned mis-match */
#pragma warning (disable:1418)		/* external function definition with no prior declaration */
#pragma warning (disable:1599)		/* declaration hides variable */
#pragma warning (disable:981)		/* operands are evaluated in unspecified order */
#pragma warning (disable:1419)		/* external declaration in primary source file */
#pragma warning (disable:1572)		/* floating-point equality and inequality comparisons are unreliable */
#pragma warning (disable:494)		/* omission of "class"is nonstandard */
#pragma warning (disable:444)		/* destructor for base class */
#else
#pragma warning (disable:4305)		/* truncation from 'double' to 'GLfloat' */
#pragma warning (disable:4244)		/* truncation from '__w64' to 'int' */
#pragma warning (disable:4267)		/* conversion from size_t to int */
#pragma warning (disable:4018)		/* signed/unsigned mismatch */
#pragma warning (disable:4100)		/* unreferenced formal parameter */
#pragma warning (disable:4505)		/* unreferenced local function */
#pragma warning (disable:4701)		/* potentially unitialized local variable */
#pragma warning (disable:4389)		/* signed/unsigned mismatch */
#pragma warning (disable:4189)		/* local variable set but not referenced */
#pragma warning (disable:4206)		/* translation unit empty */

#endif

#endif

#ifndef _DEBUG
#define pp_RELEASE
#endif

#define pp_DRAWISO
#define EGZ
#define USE_ZLIB
#define pp_THREADS2

#ifdef pp_ALPHA
#define pp_COLOR
#define pp_TOUR
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








