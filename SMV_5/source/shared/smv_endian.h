// $Date$ 
// $Revision$
// $Author$
#ifndef SMV_ENDIAN
#define SMV_ENDIAN
int getendian(void);
float float_switch(float val);
int int_switch(int val);
void endian_switch(void *val, int nval);
#endif

