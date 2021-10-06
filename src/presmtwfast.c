#include "survPresmooth.h"

void presmtwfast(double *x, int *nx, double *t, int *nt, double *bw, int *nkernel, int *dup, double *phat, double *ptw) {
  int i, j, counter = 0;

  for (i = 0; i < *nx; i++) {
    ptw[i] = 0;
    for (j = 0; j < *nt; j++) {
      if (dup[j] == 1)
	counter++;
      else 
	counter = 0;
      if (fabs(x[i] - t[j]) < (*bw))
	ptw[i] += kernel((x[i] - t[j]) / (*bw), *nkernel) * phat[j] / (*nt - j + counter);
    }
    ptw[i] = ptw[i] / (*bw);
  }
}

