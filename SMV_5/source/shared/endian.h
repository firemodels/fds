// $Date: 2010-06-21 14:35:03 -0400 (Mon, 21 Jun 2010) $ 
// $Revision: 6361 $
// $Author: gforney $

#define BIG_ENDIAN 1
#define LITTLE_ENDIAN 0

int getendian(void);
float float_switch(float val);
int int_switch(int val);
void endian_switch(void *val, int nval);

