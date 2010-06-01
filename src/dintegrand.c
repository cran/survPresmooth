#include "survPresmooth.h"

void dintegrand(double *grid, int *legrid, double *t, int *n, double *coef, double *par, double *p0, double *d1int, double *d2int){
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
	for (i = 0; i < *legrid; i++){
			d1int[i] = (1 - 2*p0[i])*(d1[i] * p1[i] + d0[i] * p2[i]/2) / pow(1 - ecdf[i] + 1/(*n), 2);
			d2int[i] = p0[i] * (1 - p0[i]) / pow(1 - ecdf[i] + 1/(*n), 2);
	}
	free(deriv);
	free(d0);
	free(d1);
	free(p1);
	free(p2);
	free(ecdf);
}

