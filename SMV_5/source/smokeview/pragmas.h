// $Date: 2007-11-06 16:39:27 -0500 (Tue, 06 Nov 2007) $ 
// $Revision: 938 $
// $Author: gforney $

#ifdef pp_WIN_INTEL
#pragma warning (disable:810)		/* conversion from xx to yy may lose sig bits */
#pragma warning (disable:869)		/* parameter was never referenced */
//#pragma warning (disable:4018)		/* signied/unsigned match */
#pragma warning (disable:4206)		/* translation unit empty */
//#pragma warning (disable:4267)		/* size_t to int possible loss of data */
//#pragma warning (disable:4389)		/* signed/unsigned mis-match */
#pragma warning (disable:1418)		/* external function definition with no prior declaration */
// #pragma warning (disable:1599)		/* declaration hides variable */
#pragma warning (disable:981)		/* operands are evaluated in unspecified order */
#pragma warning (disable:1419)		/* external declaration in primary source file */
#pragma warning (disable:1572)		/* floating-point equality and inequality comparisons are unreliable */
#pragma warning (disable:494)		/* omission of "class"is nonstandard */
#pragma warning (disable:444)		/* destructor for base class */
#pragma warning (disable:2259)
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
