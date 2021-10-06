#include "survPresmooth.h"

void c1integrand2(double *grid, int *legrid, double *t, int *n, double *coef, double *par, double *c1int2) {
  int i, *deriv;
  double *d0, *d1, *p1, *p2, *ecdf;

  deriv = calloc(1, sizeof(int));
  d0 = malloc(*legrid * sizeof(double));
  d1 = malloc(*legrid * sizeof(double));
  p1 = malloc(*legrid * sizeof(double));
  p2 = malloc(*legrid * sizeof(double));
  ecdf = malloc(*legrid * sizeof(double));

  dweibullder(grid, legrid, par, deriv, d0);
  *deriv = 1;
  plogistder(grid, legrid, coef, deriv, p1);
  dweibullder(grid, legrid, par, deriv, d1);
  *deriv = 2;
  plogistder(grid, legrid, coef, deriv, p2);
  ecdfuncens(grid, legrid, t, n, ecdf);
  for (i = 0; i < *legrid; i++)
    c1int2[i] = (d1[i] * p1[i] + d0[i] * p2[i]/2) / (1 - ecdf[i] + 1 / (*n));
  free(deriv);
  free(d0);
  free(d1);
  free(p1);
  free(p2);
  free(ecdf);
}

