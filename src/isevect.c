#include "survPresmooth.h"

void isevect(double *t, int *delta, int *n, int *nboot, double *gridise, int *legridise, double *gridbw1, int *legridbw1, double *gridbw2, int *legridbw2, int *nkernel, int * dup, int *nestimand, double *phat, double *estim, int* presmoothing, double *isev){
  int i, j, k, boot, *indices, *deltaboot, *pnull;
  double *pnull2, *ptemp, *estimboot, *tboot, *integrand, *isecomp, *deltabootdbl;

  indices = malloc(*n * sizeof(int));
  ptemp = malloc(*n * sizeof(double));
  estimboot = malloc(*legridise * sizeof(double));
  tboot = malloc(*n * sizeof(double));
  integrand = malloc(*legridise * sizeof(double));
  isecomp = malloc(sizeof(double));

  GetRNGstate();
  if(*presmoothing == 1){ // with presmoothing
    deltaboot = malloc(*n * sizeof(int));
    switch(*nestimand){
// S
    case 1:		
      pnull = calloc(1, sizeof(int));
      pnull2 = calloc(1, sizeof(double));
      for (boot = 0; boot < *nboot; boot++){
	R_FlushConsole();
	R_ProcessEvents();
	for (i = 0; i < *n; i++)
	  indices[i] = (int)ftrunc(runif(0, 1) * (*n));
	R_isort(indices, *n);
	for (i = 0; i < *n; i++){
	  tboot[i] = t[indices[i]];
	  deltaboot[i] = (int)rbinom(1, phat[indices[i]]);
	}
	for (i = 0; i < *legridbw1; i++){
	  nadarayawatson(tboot, n, tboot, deltaboot, n, &(gridbw1[i]), nkernel, ptemp);
	  presmestim(gridise, legridise, tboot, n, pnull2, pnull, pnull, ptemp, pnull, nestimand, estimboot);
	  for (j = 0; j < *legridise; j++)
	    integrand[j] = (estimboot[j] - estim[j])*(estimboot[j] - estim[j]);
	  simpson(integrand, legridise, isecomp);
	  isev[i] += *isecomp;
	}
      }
      free(pnull);
      free(pnull2);
      break;
// H
    case 2:
      pnull = calloc(1, sizeof(int));
      for (boot = 0; boot < *nboot; boot++){
	R_FlushConsole();
	R_ProcessEvents();
	for (i = 0; i < *n; i++)
	  indices[i] = (int)ftrunc(runif(0, 1) * (*n));
	R_isort(indices, *n);
	for (i = 0; i < *n; i++){
	  tboot[i] = t[indices[i]];
	  deltaboot[i] = (int)rbinom(1, phat[indices[i]]);
	}
	for (i = 0; i < *legridbw1; i++){
	  nadarayawatson(tboot, n, tboot, deltaboot, n, &(gridbw1[i]), nkernel, ptemp);
	  presmestim(gridise, legridise, tboot, n, gridbw2, nkernel, pnull, ptemp, dup, nestimand, estimboot);
	  for (j = 0; j < *legridise; j++)
	    integrand[j] = (estimboot[j] - estim[j])*(estimboot[j] - estim[j]);
	  simpson(integrand, legridise, isecomp);
	  isev[i] += *isecomp;
	}
      }
      free(pnull);
      break;
// f
    case 3:
      for (boot = 0; boot < *nboot; boot++){
	R_FlushConsole();
	R_ProcessEvents();
	for (i = 0; i < *n; i++)
	  indices[i] = (int)ftrunc(runif(0, 1) * (*n));
	R_isort(indices, *n);
	for (i = 0; i < *n; i++){
	  tboot[i] = t[indices[i]];
	  deltaboot[i] = (int)rbinom(1, phat[indices[i]]);
	}
	for (j = 0; j < *legridbw2; j++)
	  for (i = 0; i < *legridbw1; i++){
	    nadarayawatson(tboot, n, tboot, deltaboot, n, &(gridbw1[i]), nkernel, ptemp);
	    presmdensfast(gridise, legridise, tboot, n, &(gridbw2[j]), nkernel, ptemp, estimboot);
	    for (k = 0; k < *legridise; k++)
	      integrand[k] = (estimboot[k] - estim[k])*(estimboot[k] - estim[k]);
	    simpson(integrand, legridise, isecomp);
	    isev[j * (*legridbw1) + i] += *isecomp;
	  }
      }
      break;
// h
    case 4:
      for (boot = 0; boot < *nboot; boot++){
	R_FlushConsole();
	R_ProcessEvents();
	for (i = 0; i < *n; i++)
	  indices[i] = (int)ftrunc(runif(0, 1) * (*n));
	R_isort(indices, *n);
	for (i = 0; i < *n; i++){
	  tboot[i] = t[indices[i]];
	  deltaboot[i] = (int)rbinom(1, phat[indices[i]]);
	}
	for (j = 0; j < *legridbw2; j++)
	  for (i = 0; i < *legridbw1; i++){
	    nadarayawatson(tboot, n, tboot, deltaboot, n, &(gridbw1[i]), nkernel, ptemp);
	    presmtwfast(gridise, legridise, tboot, n, &(gridbw2[j]), nkernel, dup, ptemp, estimboot);
	    for (k = 0; k < *legridise; k++)
	      integrand[k] = (estimboot[k] - estim[k])*(estimboot[k] - estim[k]);
	    simpson(integrand, legridise, isecomp);
	    isev[j * (*legridbw1) + i] += *isecomp;
	  }
      }
      break;
    default:
      break;
    }
    free(deltaboot);
  }
  else{ // without presmoothing
    deltabootdbl = malloc(*n * sizeof(double));
    if(*nestimand == 3){
// f
      for (boot = 0; boot < *nboot; boot++){
	R_FlushConsole();
	R_ProcessEvents();
	for (i = 0; i < *n; i++)
	  indices[i] = (int)ftrunc(runif(0, 1) * (*n));
	R_isort(indices, *n);
	for (i = 0; i < *n; i++){
	  tboot[i] = t[indices[i]];
	  deltabootdbl[i] = (double)delta[indices[i]];
	}
	for (i = 0; i < *legridbw2; i++){
	  presmdensfast(gridise, legridise, tboot, n, &(gridbw2[i]), nkernel, deltabootdbl, estimboot);
	  for (j = 0; j < *legridise; j++)
	    integrand[j] = (estimboot[j] - estim[j])*(estimboot[j] - estim[j]);
	  simpson(integrand, legridise, isecomp);
	  isev[i] += *isecomp;
	}
      }
    }
    else{
// h
      for (boot = 0; boot < *nboot; boot++){
	R_FlushConsole();
	R_ProcessEvents();
	for (i = 0; i < *n; i++)
	  indices[i] = (int)ftrunc(runif(0, 1) * (*n));
	R_isort(indices, *n);
	for (i = 0; i < *n; i++){
	  tboot[i] = t[indices[i]];
	  deltabootdbl[i] = (double)delta[indices[i]];
	}
	for (i = 0; i < *legridbw2; i++){
	  presmtwfast(gridise, legridise, tboot, n, &(gridbw2[i]), nkernel, dup, deltabootdbl, estimboot);
	  for (j = 0; j < *legridise; j++)
	    integrand[j] = (estimboot[j] - estim[j])*(estimboot[j] - estim[j]);
	  simpson(integrand, legridise, isecomp);
	  isev[i] += *isecomp;
	}
      }
    }
    free(deltabootdbl);
  }
  PutRNGstate();
  free(indices);
  free(ptemp);
  free(estimboot);
  free(tboot);
  free(integrand);
  free(isecomp);
}
