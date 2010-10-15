// $Date$ 
// $Revision$
// $Author$

#define BIG_ENDIAN 1
#define LITTLE_ENDIAN 0

int getendian(void);
float float_switch(float val);
int int_switch(int val);
void endian_switch(void *val, int nval);

