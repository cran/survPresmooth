#include "survPresmooth.h"

void termsmise(double *t, int *delta, int *n, double *esf, double *grid, int *legrid, double *step, double *bw, int *nkernel, int *nestimand, double *alpha, double *integral1, double *integral2, double *integral3, double *integral4, double *integral5){

	int i, *temp, *pnull;
	double *pnull2, *p, *p1, *p2, *pt, *S, *f2, *h, *h1, *h2, *integrand1, *integrand2, *integrand3, *integrand4, *integrand5;

	temp = calloc(1, sizeof(int));
	pnull = calloc(1, sizeof(int));
	pnull2 = calloc(1, sizeof(double));
	p = malloc(*legrid * sizeof(double));
	p1 = malloc(*legrid * sizeof(double));
	p2 = malloc(*legrid * sizeof(double));
	pt = malloc(*n * sizeof(double));
	S = malloc(*legrid * sizeof(double));
	f2 = malloc(*legrid * sizeof(double));
	h = malloc(*legrid * sizeof(double));
	h1 = malloc(*legrid * sizeof(double));
	h2 = malloc(*legrid * sizeof(double)); 
	integrand1 = malloc(*legrid * sizeof(double));
	integrand2 = malloc(*legrid * sizeof(double));
	integrand3 = malloc(*legrid * sizeof(double));
	integrand4 = malloc(*legrid * sizeof(double));
	integrand5 = malloc(*legrid * sizeof(double));
	nadarayawatsonder(grid, legrid, t, delta, n, &(bw[0]), nkernel, p, p1, p2);
	densuncens(grid, legrid, t, n, &(bw[1]), nkernel, temp, h);
	*temp = 1;
	densuncens(grid, legrid, t, n, &(bw[1]), nkernel, temp, h1);
    if(*nestimand == 3){
// f
		nadarayawatson(t, n, t, delta, n, &(bw[0]), nkernel, pt);
		presmestim(grid, legrid, t, n, pnull2, pnull, pnull, pt, pnull, temp, S);
		presmdens2der(grid, legrid, t, n, &(bw[2]), nkernel, pt, f2);
		for (i = 0; i < *legrid; i++){
			integrand1[i] = pow(f2[i], 2);
			integrand2[i] = f2[i] * 2 * S[i] / esf[i] * (p1[i] * h1[i] + p2[i] * h[i] / 2  - p[i] * h[i] * alpha[i]);
			integrand3[i] = pow(2 * S[i] / esf[i] * (p1[i] * h1[i] + p2[i] * h[i] / 2 - p[i] * h[i] * alpha[i]), 2);
			integrand4[i] = h[i] * pow(p[i] * S[i] / esf[i], 2);
			integrand5[i] = h[i] * p[i] * (1-p[i]) * pow(S[i]/esf[i], 2);
		}
	}
    else if (*nestimand == 4){
// h
		*temp = 2;
		densuncens(grid, legrid, t, n, &(bw[1]), nkernel, temp, h2);
		for (i = 0; i < *legrid; i++){
			integrand1[i] = pow(((h2[i] + 3*h[i]*h1[i]/esf[i] + 2*pow(h[i],3)/pow(esf[i],2))*p[i] + 2*(h1[i] + pow(h[i],2)/esf[i])*p1[i] + h[i]*p2[i])/esf[i], 2);
			integrand2[i] = ((h2[i] + 3*h[i]*h1[i]/esf[i] + 2*pow(h[i],3)/pow(esf[i],2))*p[i] + 2*(h1[i] + pow(h[i],2)/esf[i])*p1[i] + h[i]*p2[i]) * (h[i]*p2[i] + 2*h1[i]*p1[i]/pow(esf[i],2));
			integrand3[i] = pow((h[i]*p2[i] + 2*h1[i]*p1[i])/esf[i],2);
			integrand4[i] = h[i] * pow(p[i],2) / pow(esf[i],2);
			integrand5[i] = h[i]*p[i]*(1 - p[i]) / pow(esf[i],2);
		}
	}
	simpson(integrand1, legrid, integral1);
	simpson(integrand2, legrid, integral2);
	simpson(integrand3, legrid, integral3);
	simpson(integrand4, legrid, integral4);
	simpson(integrand5, legrid, integral5);
    *integral1 *= *step;
	*integral2 *= *step;
	*integral3 *= *step;
	*integral4 *= *step;
	*integral5 *= *step;
	free(temp);
	free(pnull);
	free(pnull2);
	free(p);
	free(p1);
	free(p2);
	free(pt);
	free(S);
	free(f2);
	free(h);
	free(h1);
    free(h2);
	free(integrand1);
	free(integrand2);
	free(integrand3);
	free(integrand4);
	free(integrand5);
}





