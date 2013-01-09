#include "survPresmooth.h"

void densuncens(double *x, int *nx, double *t, int *nt, double *bw, int *nkernel, int *deriv, double *duncens){
  int i, j;
	         
  for (i = 0; i < *nx; i++){
    duncens[i] = 0;
    for (j = 0; j < *nt; j++)
      if (fabs(x[i] - t[j]) < (*bw))
	duncens[i] += kernelder((x[i] - t[j])/(*bw), *nkernel, *deriv);
    duncens[i] = duncens[i]/pow(*bw, *deriv + 1)/(*nt);
  }
}

