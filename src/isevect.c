#include "survPresmooth.h"
void isevect(double *t, int *delta, int *n, int *nboot, double *gridise, int *legridise, double *gridbw1, int *legridbw1, double *gridbw2, int *legridbw2,
int *nkernel, int * dup, int *nestimand, double *phat, double *estim, double *isev){
	int i, j, k, boot, *indices, *deltaboot, *pnull;
	double *pnull2, *ptemp, *estimboot, *tboot;
	
	pnull = calloc(1, sizeof(int));
	pnull2 = calloc(1, sizeof(double));
	indices = malloc(*n * sizeof(int));
    ptemp = malloc(*n * sizeof(double));
	estimboot = malloc(*legridise * sizeof(double));
	tboot = malloc(*n * sizeof(double));
	deltaboot = malloc(*n * sizeof(int));

    GetRNGstate();
	switch(*nestimand) {
// S
	case 1:		
		for (boot = 0; boot < *nboot; boot++){
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
				isev[boot * (*legridbw1) + i] = (estimboot[0] - estim[0])*(estimboot[0] - estim[0])/2;
				for (k = 1; k < (*legridise - 1); k++)
					isev[boot * (*legridbw1) + i] += (estimboot[k] - estim[k])*(estimboot[k] - estim[k]);
				isev[boot * (*legridbw1) + i] += (estimboot[k] - estim[k])*(estimboot[k] - estim[k])/2;
			}
		}
		break;
// H
	case 2:
		for (boot = 0; boot < *nboot; boot++){
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
				isev[boot * (*legridbw1) + i] = (estimboot[0] - estim[0])*(estimboot[0] - estim[0])/2;
				for (k = 1; k < (*legridise - 1); k++)
					isev[boot * (*legridbw1) + i] += (estimboot[k] - estim[k])*(estimboot[k] - estim[k]);
				isev[boot * (*legridbw1) + i] += (estimboot[k] - estim[k])*(estimboot[k] - estim[k])/2;
			}
		}
		break;
// f

	case 3:
		for (boot = 0; boot < *nboot; boot++){
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
					isev[boot * (*legridbw1) * (*legridbw2) + j * (*legridbw1) + i] = (estimboot[0] - estim[0])*(estimboot[0] - estim[0])/2;
					for (k = 1; k < (*legridise - 1); k++)
						isev[boot * (*legridbw1) * (*legridbw2) + j * (*legridbw1) + i] += (estimboot[k] - estim[k])*(estimboot[k] - estim[k]);
					isev[boot * (*legridbw1) * (*legridbw2) + j * (*legridbw1) + i] += (estimboot[k] - estim[k])*(estimboot[k] - estim[k])/2;
				}
		}
		break;
// h
	case 4:
		for (boot = 0; boot < *nboot; boot++){
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
					isev[boot * (*legridbw1) * (*legridbw2) + j * (*legridbw1) + i] = (estimboot[0] - estim[0])*(estimboot[0] - estim[0])/2;
					for (k = 1; k < (*legridise - 1); k++)
						isev[boot * (*legridbw1) * (*legridbw2) + j * (*legridbw1) + i] += (estimboot[k] - estim[k])*(estimboot[k] - estim[k]);
					isev[boot * (*legridbw1) * (*legridbw2) + j * (*legridbw1) + i] += (estimboot[k] - estim[k])*(estimboot[k] - estim[k])/2;
				}
		}
		break;
	default:
		break;
	}
	PutRNGstate();
	free(indices);
	free(pnull);
	free(pnull2);
	free(ptemp);
	free(estimboot);
	free(tboot);
	free(deltaboot);
}

