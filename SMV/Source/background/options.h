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
#pragma warning (disable:1478)
#pragma warning (disable:1786)

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
#pragma warning (disable:2259)		/* non-pointer conversion from "double" to "float" ma lose significant bits */

#define _CRT_SECURE_NO_WARNINGS
#endif

#undef pp_release
#define pp_release

#ifdef pp_release
#define PROGVERSION "1.1"
#endif

#ifdef CPP
#define CCC "C"
#else
#define CCC
#endif

#ifdef INMAIN
#define SVEXTERN
#define SVDECL(var,val)  var=val
#else
#define SVEXTERN extern CCC
#define SVDECL(var,val)  var
#endif

#ifdef CPP
#define EXTERNCPP extern "C"
#else
#define EXTERNCPP
#endif

#define FILE_SIZE unsigned int

#ifdef X64
#define STRUCTSTAT struct __stat64
#define STAT _stat64
#else
#define STRUCTSTAT struct stat
#define STAT stat
#endif


