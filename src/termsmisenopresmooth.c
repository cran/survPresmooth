#include "survPresmooth.h"

void termsmisenopresmooth(double *t, int *delta, int *n, double *esf, double *grid, int *legrid, double *step, double *bw, int *nkernel, int *nestimand, double *integral1, double *integral2){

  int i, *temp, *pnull;
  double *pnull2, *p, *p1, *p2, *pt, *S, *f2, *h, *h1, *h2, *integrand1, *integrand2, *deltadbl;
  
  temp = calloc(1, sizeof(int));
  p = malloc(*legrid * sizeof(double));
  p1 = malloc(*legrid * sizeof(double));
  p2 = malloc(*legrid * sizeof(double));
  h = malloc(*legrid * sizeof(double));
  integrand1 = malloc(*legrid * sizeof(double));
  integrand2 = malloc(*legrid * sizeof(double));

  nadarayawatsonder(grid, legrid, t, delta, n, &(bw[0]), nkernel, p, p1, p2);
  densuncens(grid, legrid, t, n, &(bw[1]), nkernel, temp, h);
  *temp = 1;
  if(*nestimand == 3){
//f
    pnull = calloc(1, sizeof(int));
    pnull2 = calloc(1, sizeof(double));
    pt = malloc(*n * sizeof(double));
    S = malloc(*legrid * sizeof(double));
    f2 = malloc(*legrid * sizeof(double));
    deltadbl = malloc(*n * sizeof(double));
    nadarayawatson(t, n, t, delta, n, &(bw[0]), nkernel, pt);
// without presmoothing, deltas instead of pt
    for (i = 0; i < *n; i++) // type conversion from integer to double
      deltadbl[i] = (double)delta[i];
    presmestim(grid, legrid, t, n, pnull2, pnull, pnull, deltadbl, pnull, temp, S);
    presmdens2der(grid, legrid, t, n, &(bw[2]), nkernel, pt, f2);
    for (i = 0; i < *legrid; i++){
      integrand1[i] = pow(f2[i], 2);
      integrand2[i] = h[i] * p[i] * pow(S[i] / esf[i], 2);
    }
    free(pnull);
    free(pnull2);
    free(pt);
    free(S);
    free(f2);
    free(deltadbl);
  }
  else 
    if (*nestimand == 4){
// h
      h1 = malloc(*legrid * sizeof(double));
      h2 = malloc(*legrid * sizeof(double)); 
      densuncens(grid, legrid, t, n, &(bw[1]), nkernel, temp, h1);
      *temp = 2;
      densuncens(grid, legrid, t, n, &(bw[1]), nkernel, temp, h2);
      for (i = 0; i < *legrid; i++){
	integrand1[i] = pow(((h2[i] + 3*h[i]*h1[i]/esf[i] + 2*pow(h[i],3)/pow(esf[i],2))*p[i] + 2*(h1[i] + pow(h[i],2)/esf[i])*p1[i] + h[i]*p2[i])/esf[i], 2);
	integrand2[i] = h[i] * p[i] / pow(esf[i],2);
      }
      free(h1);
      free(h2);
    }
  simpson(integrand1, legrid, integral1);
  simpson(integrand2, legrid, integral2);
  *integral1 *= *step;
  *integral2 *= *step;
  free(temp);
  free(p);
  free(p1);
  free(p2);
  free(h);
  free(integrand1);
  free(integrand2);
}


