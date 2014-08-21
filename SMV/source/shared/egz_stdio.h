// $Date$ 
// $Revision$
// $Author$

#ifndef EGZ_H_DEFINED
#define EGZ_H_DEFINED

#ifdef WIN32

#define _WINDOWS

/* #define ZLIB_DLL */
#endif   /* end WIN32 */

#ifndef DEF_ZLIB
#define DEF_ZLIB
#include <zlib.h>
#endif

typedef struct {
  int endianswitch;
  int compression;
  gzFile *stream;
  FILE *stream2;
} EGZ_FILE;
int EGZ_FCLOSE(EGZ_FILE *egz_stream);
EGZ_FILE *EGZ_FOPEN(const char *file, const char *mode, int compress, int endian);
size_t EGZ_FREAD( void *buffer, size_t size, size_t count, const EGZ_FILE *stream );
size_t EGZ_FWRITE( void *buffer, size_t size, size_t count, const EGZ_FILE *egz_stream );
int EGZ_FEOF( const EGZ_FILE *egz_stream );
int EGZ_FSEEK( const EGZ_FILE *stream, long offset, int origin );
long EGZ_FTELL( const EGZ_FILE *stream );
void EGZ_REWIND(const EGZ_FILE *egz_stream);
#endif
