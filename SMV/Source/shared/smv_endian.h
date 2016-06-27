#ifndef SMV_ENDIAN_H_DEFINED
#define SMV_ENDIAN_H_DEFINED
int getendian(void);
float float_switch(float val);
int int_switch(int val);
void endian_switch(void *val, int nval);
#endif

