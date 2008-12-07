// $Date: 2008-10-30 14:40:35 -0400 (Thu, 30 Oct 2008) $ 
// $Revision: 2576 $
// $Author: gforney $


 /* ------------------ options ------------------------ */
#ifdef WIN32
#ifdef _DEBUG
#pragma float_control( precise, on)
#pragma float_control( except, on )
#endif
#pragma warning (disable:4996)		
#pragma warning (disable:4701)		
#pragma warning (disable:4310)		
#pragma warning (disable:4127)		
#pragma warning (disable:4244)		

#pragma warning (disable:810)		/* conversion from xx to yy may lose sig bits */
#pragma warning (disable:869)		/* parameer was never referenced */
#pragma warning (disable:4018)		/* signied/unsigned match */
#pragma warning (disable:4206)		/* translation unit empty */
#pragma warning (disable:4267)		/* size_t to int possible loss of data */
#pragma warning (disable:4389)		/* signed/unsigned mis-match */
#pragma warning (disable:981)		/* operands are evaluated in unspecified order */
#pragma warning (disable:494)		/* omission of "class"is nonstandard */
#pragma warning (disable:444)		/* destructor for base class */

#define _CRT_SECURE_NO_WARNINGS
#endif

#undef pp_release
#define pp_release

#ifdef pp_ALPHA
#undef pp_release
#define SPRVERSION "experimental"
#endif

#ifdef pp_BETA
#undef pp_release
#define SPRVERSION "test"
#endif

#ifdef pp_release
#define SPRVERSION "1.0.0"
#endif
