#include "survPresmooth.h"

void ecdfuncens(double *x, int *nx, double *t, int *nt, double *ecdf){
	int i, j;
	         
	for (i = 0; i < *nx; i++){
		ecdf[i] = 0;
		for (j = 0; j < *nt; j++)
			if (t[j] <= x[i])
				ecdf[i] += 1.0 / (*nt);
	}
}

