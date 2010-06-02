#include "survPresmooth.h"

void presmtwfast(double *x, int *nx, double *t, int *nt, double *bw, int *nkernel, double *phat, double *ptw){
	int i, j;

	for (i = 0; i < *nx; i++){
		ptw[i] = 0;
		for (j = 0; j < *nt; j++)
			if (fabs(x[i] - t[j]) < (*bw))
			   ptw[i] += kernel((x[i] - t[j])/(*bw), *nkernel) * phat[j] / (*nt - j);
		ptw[i] = ptw[i] / (*bw);
	}
}

