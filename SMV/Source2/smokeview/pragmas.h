#ifndef PRAGMAS_H_DEFINED
#define PRAGMAS_H_DEFINED
#ifdef pp_INTEL
#pragma warning (disable:777) // Testing floats for equality
#pragma warning (disable:776) // Possible truncation of addition
#pragma warning (disable:712) // Loss of precision
#pragma warning (disable:1786)
#pragma warning (disable:177)
#pragma warning (disable:2557)		/* comparison between signed and unsigned operands */
#pragma warning (disable:695)		/* calling convenetion if ignored for this type */
#pragma warning (disable:810)		/* conversion from xx to yy may lose sig bits */
#pragma warning (disable:869)		/* parameter was never referenced */
#pragma warning (disable:4206)		/* translation unit empty */
#pragma warning (disable:1418)		/* external function definition with no prior declaration */
#pragma warning (disable:981)		/* operands are evaluated in unspecified order */
#pragma warning (disable:1419)		/* external declaration in primary source file */
#pragma warning (disable:1572)		/* floating-point equality and inequality comparisons are unreliable */
#pragma warning (disable:494)		/* omission of "class"is nonstandard */
#pragma warning (disable:444)		/* destructor for base class */
#pragma warning (disable:2259)
#pragma warning (disable:1678)		/* cannot enable speculation unless fenv_access and exception_semantics are disabled */
#else
#pragma warning (disable:47)
#pragma warning (disable:4305)		/* truncation from 'double' to 'GLfloat' */
#pragma warning (disable:4244)		/* truncation from '__w64' to 'int' */
#pragma warning (disable:4267)		/* conversion from size_t to int */
#pragma warning (disable:4018)		/* signed/unsigned mismatch */
#pragma warning (disable:4100)		/* unreferenced formal parameter */
#pragma warning (disable:4127)		/* conditional expression is constant */
#pragma warning (disable:4505)		/* unreferenced local function */
#pragma warning (disable:4701)		/* potentially unitialized local variable */
#pragma warning (disable:4389)		/* signed/unsigned mismatch */
#pragma warning (disable:4189)		/* local variable set but not referenced */
#pragma warning (disable:4206)		/* translation unit empty */
#pragma warning (disable:4310)
#pragma warning (disable:4005)
#pragma warning (disable:4245)
#endif
#endif
