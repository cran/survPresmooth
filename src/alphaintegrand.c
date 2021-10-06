#include "survPresmooth.h"

void alphaintegrand(double *t, int *delta, int *n, double *grid, int *legrid, double *bw1, double *bw2, int *nkernel, double *alphaint) {
  int i, *temp;
  double *p, *p1, *p2, *d, *d1, *ecdf;

  temp = calloc(1, sizeof(int));
  p = malloc(*legrid * sizeof(double));
  p1 = malloc(*legrid * sizeof(double));
  p2 = malloc(*legrid * sizeof(double));
  d = malloc(*legrid * sizeof(double));
  d1 = malloc(*legrid * sizeof(double));
  ecdf = malloc(*legrid * sizeof(double));
  nadarayawatsonder(grid, legrid, t, delta, n, bw1, nkernel, p, p1, p2);
  densuncens(grid, legrid, t, n, bw2, nkernel, temp, d);
  ecdfuncens(grid, legrid, t, n, ecdf);
  *temp = 1;
  densuncens(grid, legrid, t, n, bw2, nkernel, temp, d1);
  for (i = 0; i < *legrid; i++)
    alphaint[i] = (d1[i] * p1[i] + d[i] * p2[i]/2) / (1 - ecdf[i] + 1 / (*n));
  free(temp);
  free(p);
  free(p1);
  free(p2);
  free(d);
  free(d1);
  free(ecdf);
}
