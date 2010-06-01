#ifndef PRESMOOTSURV_H
#define PRESMOOTSURV_H

#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

void Ak(double *x, int *nkernel, double *value);

void alphaintegrand(double *t, int *delta, int *n, double *grid, int *legrid, double *bw1, double *bw2, int *nkernel, double *alphaint);

void densuncens(double *x, int *nx, double *t, int *nt, double *bw, int *nkernel, int *deriv, double *duncens);

void ecdfuncens(double *x, int *nx, double *t, int *nt, double *ecdf);

void dweibullder(double *x, int *nx, double *par, int *deriv, double *dder);

double kernelboundary(double x, double q, int nkernel);

double kernel(double x, int nkernel);

double kernelder(double x, int nkernel, int deriv);

void nadarayawatson(double *x, int *nx, double *t, int *delta, int *nt, double *bw, int *nkernel, double *phat);

void nadarayawatsonder(double *x, int *nx, double *t, int *delta, int *nt, double *bw, int *nkernel, double *phat, double *p1hat, double *p2hat);

void plogistder(double *x, int *nx, double *coef, int *deriv, double *pder);

void presmestim(double *x, int *nx, double *t, int *nt, double *bw, int *nkernel, int *nbound, double *phat, int *nestimand, double *pest);

void presmtwfast(double *x, int *nx, double *t, int *nt, double *bw, int *nkernel, double *phat, double *ptw);

void presmdensfast(double *x, int *nx, double *t, int *nt, double *bw, int *nkernel, double *phat, double *pd);

void nadarayawatson(double *x, int *nx, double *t, int *delta, int *nt, double *bw, int *nkernel, double *phat);

void simpson(double *integrand, int *lintegrand, double *step, double *integral);

void weightspresmkm(double *t, int *nt, double *phat, double *w);

#endif


