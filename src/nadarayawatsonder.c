#include "survPresmooth.h"

void nadarayawatsonder(double *x, int *nx, double *t, int *delta, int *nt, double *bw, int *nkernel, double *phat, double *p1hat, double *p2hat){
  int i, j;
  double *psi, *psi1, *psi2, *h, *h1, *h2;
	
  psi = calloc(*nx, sizeof(double));
  psi1 = calloc(*nx, sizeof(double));
  psi2 = calloc(*nx, sizeof(double));
  h = calloc(*nx, sizeof(double));
  h1 = calloc(*nx, sizeof(double));
  h2 = calloc(*nx, sizeof(double));

  for (i = 0; i < *nx; i++){
    for (j = 0; j < *nt; j++)
      if (fabs(x[i] - t[j]) < *bw){
	h[i] += kernel((x[i] - t[j])/(*bw), *nkernel)/(*bw)/(*nt);
	h1[i] += kernelder((x[i] - t[j])/(*bw), *nkernel, 1)/(*bw)/(*bw)/(*nt);
	h2[i] += kernelder((x[i] - t[j])/(*bw), *nkernel, 2)/(*bw)/(*bw)/(*bw)/(*nt);
	if (delta[j] == 1) {
	  psi[i] += kernel((x[i] - t[j])/(*bw), *nkernel)/(*bw)/(*nt);
	  psi1[i] += kernelder((x[i] - t[j])/(*bw), *nkernel, 1)/(*bw)/(*bw)/(*nt);
	  psi2[i] += kernelder((x[i] - t[j])/(*bw), *nkernel, 2)/(*bw)/(*bw)/(*bw)/(*nt);
	}
      }  
    if (h[i] < 0.00000000001){
      phat[i] = 0;
      p1hat[i] = 0;
      p2hat[i] = 0;
    }
    else{
      phat[i] = psi[i]/h[i];	
      p1hat[i] = (psi1[i]*h[i] - psi[i]*h1[i])/h[i]/h[i];
      p2hat[i] = (psi2[i]*h[i]*h[i] - psi[i]*h2[i]*h[i] - 2*psi1[i]*h1[i]*h[i] + 2*psi[i]*h1[i]*h1[i])/h[i]/h[i]/h[i];	
    }
  }
  free(psi);
  free(psi1);
  free(psi2);
  free(h);
  free(h1);
  free(h2);
}

