// $Date$ 
// $Revision$
// $Author$

#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNING


 /* ------------------ options ------------------------ */


#pragma warning (disable:1419)		/* external declaration in primary source file */

#define REG_GET 0
#define REG_SET 1
#define REG_USER_PATH 2
#define REG_SYSTEM_PATH 3

int reg_path(int setget, int pathtype, char *path);
