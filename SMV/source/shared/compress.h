#ifndef COMPRESS_H_DEFINED
#define COMPRESS_H_DEFINED
#ifndef DEF_ZLIB
#define DEF_ZLIB
#include <zlib.h>
#endif
int compress_zlib(unsigned char *dest, uLongf *destLen, unsigned char *source, int sourceLen);
int uncompress_zlib(unsigned char *dest, uLongf *destLen, unsigned char *source, int sourceLen);
unsigned int compress_rle(unsigned char *buffer_in, int nchars_in, unsigned char *buffer_out);
unsigned int uncompress_rle(unsigned char *buffer_in, int nchars_in, unsigned char *buffer_out);
void compress_volsliceframe(float *data_in, int n_data_in,
                float timeval_in, float *vmin_in, float *vmax_in,
                unsigned char **compressed_data_out, uLongf *ncompressed_data_out);
int uncompress_volsliceframe(unsigned char *compressed_data_in,
                           float *data_out, int n_data_in, float *timeval_out,
                           unsigned char *fullbuffer);
#endif

