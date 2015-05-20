// $Date: 2012-08-20 19:39:01 -0400 (Mon, 20 Aug 2012) $ 
// $Revision: 12156 $
// $Author: koverholt $
#ifndef SMV_ENDIAN_H_DEFINED
#define SMV_ENDIAN_H_DEFINED
int getendian(void);
float float_switch(float val);
int int_switch(int val);
void endian_switch(void *val, int nval);
#endif

