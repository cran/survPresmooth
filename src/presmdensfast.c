#include "survPresmooth.h"

void presmdensfast(double *x, int *nx, double *t, int *nt, double *bw, int *nkernel, double *phat, double *pd){
	int i, j;
	double *w;

    w = malloc(*nt * sizeof(double));
	weightspresmkm(t, nt, phat, w);
	for (i = 0; i < *nx; i++){
		pd[i] = 0;
		for (j = 0; j < *nt; j++)
			if (fabs(x[i] - t[j]) < *bw)
				pd[i] += kernel((x[i] - t[j])/(*bw), *nkernel) * w[j];
		pd[i] = pd[i]/ (*bw);
	}
	free(w);
}
