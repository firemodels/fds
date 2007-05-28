#include "options.h"
#include <stdio.h>
#include <math.h>

void smoothlabel(float *a, float *b, int n){
  float delta, factor, logdelta;
  int ndigits;

  delta = (*b-*a)/(n-2);
  if(delta==0.0)return;
  logdelta = log10((double)delta);
  ndigits=logdelta-1;
  if(logdelta<=1)ndigits--;
  factor = 5*pow(10,ndigits);
  delta = (int)(delta/factor + 0.5f)*factor;

  *a = factor*(int)(*a/factor+0.5f);
  *b = *a + (n-2)*delta;

}

