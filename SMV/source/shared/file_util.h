#ifndef FILE_UTIL_H_DEFINED
#define FILE_UTIL_H_DEFINED

#include <time.h>
#ifdef __MINGW32__
#include <stdio.h>
#include "options.h"
#endif

#ifdef WIN32
#define UNLINK _unlink
#else
#define UNLINK unlink
#endif

typedef struct {
  char *file;
  int type;
} filelistdata;

#ifdef X64
#define FSEEK(a,b,c) _fseeki64(a,b,c)
#define FTELL(a) _ftelli64(a)
#else
#define FSEEK(a,b,c) fseeko(a,b,c)
#define FTELL(a) ftello(a)
#endif

#define REPLACE_FILE 0
#define APPEND_FILE 1

EXTERNCPP int FFLUSH(void);
EXTERNCPP int PRINTF(const char * format, ...);
EXTERNCPP void set_stdout(FILE *stream);
EXTERNCPP void getfilesizelabel(int size, char *sizelabel);
EXTERNCPP void copyfile(char *destdir, char *filein, char *fileout, int mode);
EXTERNCPP char *get_smokezippath(char *progdir);
EXTERNCPP int have_prog(char *prog);
EXTERNCPP int filecat(char *file_in1, char *file_in2, char *file_out);
EXTERNCPP void make_outfile(char *outfile, char *destdir, char *file1, char *ext);
EXTERNCPP void fullfile(char *fileout, char *dir, char *file);
EXTERNCPP char *get_filename(char *temp_dir, char *file, int flag);
EXTERNCPP char *get_basefilename(char *buffer, char *file);

EXTERNCPP char *setdir(char *argdir);
EXTERNCPP int getfileinfo(char *filename, char *sourcedir, FILE_SIZE *filesize);
EXTERNCPP char *get_zonefilename(char *buffer);
EXTERNCPP int can_write_to_dir(char *dir);
EXTERNCPP int file_exists(char *filename);

EXTERNCPP void free_filelist(filelistdata *filelist, int *nfilelist);
EXTERNCPP int get_nfilelist(const char *path, char *key) ;
EXTERNCPP int get_filelist(const char *path, char *key, int maxfiles, filelistdata **filelist);
EXTERNCPP char *which(char *progname);
EXTERNCPP FILE_SIZE get_filesize(const char *filename);
EXTERNCPP time_t file_modtime(char *filename);
EXTERNCPP int is_file_newer(char *file1, char *file2);
EXTERNCPP char *getprogdir(char *progname, char **svpath);
#ifdef pp_LUA
EXTERNCPP char *getprogdirabs(char *progname, char **svpath);
#endif

EXTERNCPP char *lastname(char *argi);

#ifndef STREXTERN
#ifdef WIN32
STREXTERN char STRDECL(dirseparator[],"\\");
#else
STREXTERN char STRDECL(dirseparator[],"/");
#endif
#endif

#define DPRINTF(_fmt, ...)  fprintf(stderr, "[file %s, line %d]: " _fmt, __FILE__, __LINE__, __VA_ARGS__)

#endif
