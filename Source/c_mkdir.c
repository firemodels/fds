#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void c_mkdir(const char*directory){
  struct stat st = {0};
  if (stat(directory, &st) == -1) {
      mkdir(directory, 0700);
  }
}

