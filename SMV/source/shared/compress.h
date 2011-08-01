// $Date$ 
// $Revision$
// $Author$

#include "zlib.h"
unsigned int rle(unsigned char *buffer_in, int nchars_in, unsigned char *buffer_out);
unsigned int irle(unsigned char *buffer_in, int nchars_in, unsigned char *buffer_out);
void compress_volsliceframe(float *data, int n_data, 
                float timeval, float *vmin_in, float *vmax_in,
                unsigned char **compressed_data, uLongf *ncompressed_data);
int uncompress_volsliceframe(unsigned char *compressed_data,
                           float *data, int n_data, float *timeval,
                           unsigned char *fullbuffer);

