#include <stdio.h>
#include <time.h>

int main ()
{
  time_t rawtime;

  time ( &rawtime );
  printf ( "%i", rawtime);

  return 0;
}