#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* ------------------ main ------------------------ */

int main(int argc, char **argv){
  struct tm *tm;
  time_t t;
  char str_time[100];
  char str_date[100];

  t = time(NULL);
  tm = localtime(&t);

  strftime(str_time, sizeof(str_time), "%H:%M:%S", tm);
  strftime(str_date, sizeof(str_date), "%b %d, %Y", tm);
  printf("%s %s\n",str_date,str_time);
}
