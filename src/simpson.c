#include "survPresmooth.h"

void simpson(double *integrand, int *lintegrand, double *step, double *integral){
	int i;
	div_t resdiv;

	*integral = integrand[0];
	for (i = 1; i < (*lintegrand - 1); i++){
		resdiv = div(i, 2);
    	if(resdiv.rem == 0)
			*integral += 2 * integrand[i];
		else 
			*integral += 4 * integrand[i];
	}
	*integral += integrand[i];
	*integral *= *step / 3;
}

