#include "survPresmooth.h"

void lscv(double *t, int *delta, int *n, double *gridbw, int *legridbw, int *nkernel, double *cv){
  int i, j, k, l, *temp1, *temp2, *temp3;
  double *temp4, *phat;

  phat = malloc(sizeof(double));
  temp1 = malloc((*n-1) * sizeof(int));
  temp2 = malloc(sizeof(int));
  temp3 = malloc(sizeof(int));
  temp4 = malloc((*n-1) * sizeof(double));
  *temp2 = 1;
  *temp3 = *n-1;
  for(i = 0; i < *legridbw; i++)
    for(j = 0; j < *n; j++){
      for(k = 0; k < (*n-1); k++){
	if (k == j) l = k+1;
	else l = k;
	temp1[k] = delta[l];
	temp4[k] = t[l];
      }
      nadarayawatson(&(t[j]), temp2, temp4, temp1, temp3, &(gridbw[i]), nkernel, phat);
      cv[i] += (delta[j] - *phat) * (delta[j] - *phat);	
    }
  free(temp1);
  free(temp2);
  free(temp3);
  free(temp4);
  free(phat);
}

