#include "survPresmooth.h"

void nadarayawatson(double *x, int *nx, double *t, int *delta, int *nt, double *bw, int *nkernel, double *phat){
  int i, j;
  double temp, *num, *den;
	
  num = calloc(*nx, sizeof(double));
  den = calloc(*nx, sizeof(double));
  
  for (i = 0; i < *nx; i++){
    for (j = 0; j < *nt; j++)
      if (fabs(x[i] - t[j]) < *bw){
	temp = kernel((x[i] - t[j])/(*bw), *nkernel);
	if (delta[j] == 1) num[i] += temp;
	den[i] += temp;
      }
    if (den[i] < 0.00000000001) phat[i] = 0;
    else phat[i] = num[i]/den[i];	
  }
  free(num);
  free(den);
}

