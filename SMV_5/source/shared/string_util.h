// $Date$ 
// $Revision$
// $Author$

EXTERNCPP int can_write_to_dir(char *dir);
EXTERNCPP int file_exists(char *filename);
EXTERNCPP char *which(char *progname);
EXTERNCPP FILE_SIZE get_filesize(const char *filename);
EXTERNCPP time_t file_modtime(char *filename);
EXTERNCPP int is_file_newer(char *file1, char *file2);
EXTERNCPP char *getdir(char *progname);

EXTERNCPP char *lastname(char *argi);
EXTERNCPP void trim(char *line);
EXTERNCPP char *trim_front(char *line);
EXTERNCPP void trimzeros(char *line);
EXTERNCPP void trimmzeros(char *line);
EXTERNCPP char *STRSTR(char *c, const char *key);
EXTERNCPP void scalestring(const char *stringfrom, char *stringto, const float *scale, float range);
EXTERNCPP void scalefloat2string(float floatfrom, char *stringto, const float *scale, float range);
EXTERNCPP void num2string(char *string, float tval,float range);
EXTERNCPP char *trim_string(char *buffer);
EXTERNCPP int STRCMP(const char *s1, const char *s2);
EXTERNCPP char *get_chid(char *file);
#ifdef pp_GPU
EXTERNCPP int log_base2(float xx);
#endif
EXTERNCPP void array2string(float *vals, int nvals, char *string);
EXTERNCPP void parse_string(char *string, char **tokens, int *ntokens);
EXTERNCPP float MIN(float x,float y);
EXTERNCPP float MAX(float x,float y);
EXTERNCPP float frexp10(float x, int *exp10);


SVEXTERN char dirseparator[3];


