#include "survPresmooth.h"

void dweibullder(double *x, int *nx, double *par, int *deriv, double *dder){

	int i;
	double y[*nx];

	for(i = 0; i < *nx; i++)
		y[i] = x[i]/par[1];
	switch(*deriv + 1){
		case 1:
			for (i = 0; i < *nx; i++)
				dder[i] = par[0]/par[1]*pow(y[i], par[0] - 1)*exp(-pow(y[i], par[0]));
			break;
		case 2:
			for (i = 0; i < *nx; i++)
				dder[i] = -par[0]/pow(x[i], 2)*pow(y[i], par[0])*exp(-pow(y[i], par[0]))*(par[0]*pow(y[i], par[0]) - par[0] + 1);
			break;
		case 3:
			for (i = 0; i < *nx; i++){
				dder[i] = par[0]/pow(x[i], 3)*pow(y[i], par[0])*exp(-pow(y[i], par[0]))*(pow(par[0], 2)*pow(y[i], 2*par[0]) - 3*par[0] + 3*par[0]*pow(y[i], par[0]) - 3*pow(par[0], 2)*pow(y[i], par[0]) + pow(par[0], 2) + 2);
			}
			break;
		case 4:
			for (i = 0; i < *nx; i++)
				dder[i] = -par[0]/pow(x[i], 4)*pow(y[i], par[0])*exp(-pow(y[i], par[0]))*(6*pow(par[0], 2)*pow(y[i], 2*par[0]) - 11*par[0] - 6*pow(par[0], 3)*pow(y[i], 2*par[0]) + pow(par[0], 3)*pow(y[i], 3*par[0]) + 11*par[0]*pow(y[i], par[0]) - 18*pow(par[0], 2)*pow(y[i], par[0]) + 7*pow(par[0], 3)*pow(y[i], par[0]) + 6*pow(par[0], 2) - pow(par[0], 3) + 6);
			break;
		default:
			break;
	}
}
