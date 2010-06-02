#include "survPresmooth.h"

void c1integrand1(double *grid, int *legrid, double *t, int *n, double *coef, double *par, double *c1int1){
	int i, *deriv;
	double *d0, *d1, *d2, *d3, *p1, *p2, *p3, *p4, *ecdf;

	deriv = calloc(1, sizeof(int));
	d0 = malloc(*legrid * sizeof(double));
	d1 = malloc(*legrid * sizeof(double));
	d2 = malloc(*legrid * sizeof(double));
	d3 = malloc(*legrid * sizeof(double));
	p1 = malloc(*legrid * sizeof(double));
	p2 = malloc(*legrid * sizeof(double));
	p3 = malloc(*legrid * sizeof(double));
	p4 = malloc(*legrid * sizeof(double));
	ecdf = malloc(*legrid * sizeof(double));

	dweibullder(grid, legrid, par, deriv, d0);
	*deriv = 1;
	plogistder(grid, legrid, coef, deriv, p1);
	dweibullder(grid, legrid, par, deriv, d1);
	*deriv = 2;
	plogistder(grid, legrid, coef, deriv, p2);
	dweibullder(grid, legrid, par, deriv, d2);
	*deriv = 3;
	plogistder(grid, legrid, coef, deriv, p3);
	dweibullder(grid, legrid, par, deriv, d3);
	*deriv = 4;
	plogistder(grid, legrid, coef, deriv, p4);
	ecdfuncens(grid, legrid, t, n, ecdf);
	for (i = 0; i < *legrid; i++)
			c1int1[i] = (p4[i]*d0[i] + 4*p3[i]*d1[i] + 5*p2[i]*d2[i] + 4*p1[i]*d3[i] - 2*p1[i]*d1[i]*d2[i]/d0[i]) / (1 - ecdf[i] + 1/(*n));
	free(deriv);
	free(d0);
	free(d1);
	free(d2);
	free(d3);
	free(p1);
	free(p2);
	free(p3);
	free(p4);
	free(ecdf);
}

