#include "survPresmooth.h"

void pilot2forhintegrand(double *grid, int *legrid, double *par, int *lepar, int *modif, double *pilot2hint){
	int i, *deriv;
	double *d31, *d32, *d33;

	deriv = calloc(1, sizeof(int));
	d31 = malloc(*legrid * sizeof(double));
	d32 = malloc(*legrid * sizeof(double));
	d33 = malloc(*legrid * sizeof(double));
//f
    if(*modif == 1){
		*deriv = 2;
		dweibullder(grid, legrid, &(par[0]), deriv, d31);
		switch(*lepar){
			case 2:{
				for (i = 0; i < *legrid; i++)
					pilot2hint[i] = pow(d31[i], 2);
				}
				break;
			case 5:{
		 		dweibullder(grid, legrid, &(par[2]), deriv, d32);
				for (i = 0; i < *legrid; i++)
					pilot2hint[i] = pow(par[4] * d31[i] + (1 - par[4]) * d32[i], 2);
				}
				break;
			case 8:{
			 	dweibullder(grid, legrid, &(par[2]), deriv, d32);
				dweibullder(grid, legrid, &(par[4]), deriv, d33);
				for (i = 0; i < *legrid; i++)
					pilot2hint[i] = pow(par[6] * d31[i] + par[7] * d32[i] + (1 - par[6] - par[7]) * d33[i], 2);
				}
				break;
			default:
				break;
		}
	}
//h
	else{
		*deriv = 3;
		dweibullder(grid, legrid, &(par[0]), deriv, d31);
		switch(*lepar){
			case 2:{
				for (i = 0; i < *legrid; i++)
					pilot2hint[i] = pow(d31[i], 2);
				}
				break;
			case 5:{
		 		dweibullder(grid, legrid, &(par[2]), deriv, d32);
				for (i = 0; i < *legrid; i++)
					pilot2hint[i] = pow(par[4] * d31[i] + (1 - par[4]) * d32[i], 2);
				}
				break;
			case 8:{
		 		dweibullder(grid, legrid, &(par[2]), deriv, d32);
				dweibullder(grid, legrid, &(par[4]), deriv, d33);
				for (i = 0; i < *legrid; i++)
					pilot2hint[i] = pow(par[6] * d31[i] + par[7] * d32[i] + (1 - par[6] - par[7]) * d33[i], 2);
				}
				break;
			default:
				break;
		}
	}
	free(deriv);
	free(d31);
	free(d32);
	free(d33);
}

