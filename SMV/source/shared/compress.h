// $Date: 2014-02-14 22:53:03 -0500 (Fri, 14 Feb 2014) $ 
// $Revision: 18397 $
// $Author: gforney $

#ifndef COMPRESS_H_DEFINED
#define COMPRESS_H_DEFINED
#ifndef DEF_ZLIB
#define DEF_ZLIB
#include <zlib.h>
#endif
unsigned int rle(unsigned char *buffer_in, int nchars_in, unsigned char *buffer_out);
unsigned int irle(unsigned char *buffer_in, int nchars_in, unsigned char *buffer_out);
void compress_volsliceframe(float *data_in, int n_data_in, 
                float timeval_in, float *vmin_in, float *vmax_in,
                unsigned char **compressed_data_out, uLongf *ncompressed_data_out);
int uncompress_volsliceframe(unsigned char *compressed_data_in,
                           float *data_out, int n_data_in, float *timeval_out,
                           unsigned char *fullbuffer);
#endif

