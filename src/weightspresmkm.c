#include "survPresmooth.h"

void weightspresmkm(double *t, int *nt, double *phat, double *w){
	int i;
	double *cum;

    cum = malloc(*nt * sizeof(double));
    w[0] = phat[0] / *nt;
	cum[0] = 1;
	for (i = 1; i < *nt; i++){
        w[i] = phat[i] / (*nt - i) * (cum[i-1]);
		cum[i] = cum[i-1] * (1 - phat[i] / (*nt - i));
	}
	free(cum);
}

