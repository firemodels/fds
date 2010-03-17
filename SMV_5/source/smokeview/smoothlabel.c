// $Date$ 
// $Revision$
// $Author$

#include "options.h"
#include <stdio.h>
#include <math.h>

// svn revision character string
char smoothlabel_revision[]="$Revision$";

void smoothlabel(float *a, float *b, int n){
  double delta, factor, logdelta;
  int ndigits;

  delta = ((double)*b-(double)*a)/(double)(n-2);
  if(delta==0.0)return;
  logdelta = log10((double)delta);
  ndigits=logdelta-1;
  if(logdelta<=1)ndigits--;
  factor = 5*pow(10,ndigits);
  delta = (int)(delta/factor + 0.5d)*factor;

  *a = factor*(int)(*a/factor+0.5d);
  *b = *a + (n-2)*delta;

}

