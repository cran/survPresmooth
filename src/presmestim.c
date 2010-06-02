#include "survPresmooth.h"

void presmestim(double *x, int *nx, double *t, int *nt, double *bw, int *nkernel, int *nbound, double *phat, int *nestimand, double *psest){
	int i, j;
	double q, R = t[*nt -1]/(*bw), *w;

	w = malloc(*nt * sizeof(double));
	switch(*nestimand){
// S
		case 1:{
			for (i = 0; i < *nx; i++){
				psest[i] = 1;
				for (j = 0; (t[j] <= x[i]) && (j < *nt); j++)
					psest[i] *= 1 - phat[j] / (*nt - j);
			}
			break;
		}
// H  
		case 2:{
			for (i = 0; i < *nx; i++){
				psest[i] = 0;
				for (j = 0; (t[j] <= x[i]) && (j < *nt); j++)
					psest[i] += phat[j] / (*nt - j);
			}
			break;
		}
// f
	    case 3:{
			weightspresmkm(t, nt, phat, w);
			switch(*nbound){
				case 1:
					for (i = 0; i < *nx; i++){
						psest[i] = 0;
						for (j = 0; j < *nt; j++)
							if (fabs(x[i] - t[j]) < *bw)
								psest[i] += kernel((x[i] - t[j])/(*bw), *nkernel) * w[j];
						psest[i] = psest[i]/ (*bw);
					}
					break;
				case 2:
					for (i = 0; i < *nx; i++){
						psest[i] = 0;
						q = x[i]/(*bw);
						if (q < 1){
							for (j = 0; j < *nt; j++)
								if (fabs(x[i] - t[j]) < *bw)
									psest[i] += kernelboundary((x[i] - t[j])/(*bw), q, *nkernel) * w[j];
							if (psest[i] < 0) psest[i] = 0;
						}
						else{
							for (j = 0; j < *nt; j++)
								if (fabs(x[i] - t[j]) < *bw)
									psest[i] += kernel((x[i] - t[j])/(*bw), *nkernel) * w[j];
						}
						psest[i] = psest[i]/ (*bw);
					}
					break;
				case 3:
					for (i = 0; i < *nx; i++){
						psest[i] = 0;
						q = x[i]/(*bw);
						if (q < (R -1)){
							for (j = 0; j < *nt; j++)
								if (fabs(x[i] - t[j]) < *bw)
									psest[i] += kernel((x[i] - t[j])/(*bw), *nkernel) * w[j];
						}
						else{
							for (j = 0; j < *nt; j++)
								if (fabs(x[i] - t[j]) < *bw)
									psest[i] += kernelboundary(-(x[i] - t[j])/(*bw), R - q, *nkernel) * w[j];
							if (psest[i] < 0) psest[i] = 0;
						}
						psest[i] = psest[i]/ (*bw);
					}
					break;
				case 4:
					for (i = 0; i < *nx; i++){
						psest[i] = 0;
						q = x[i]/(*bw);
						if (q < 1){
							for (j = 0; j < *nt; j++)
								if (fabs(x[i] - t[j]) < *bw)
									psest[i] += kernelboundary((x[i] - t[j])/(*bw), q, *nkernel) * w[j];
							if (psest[i] < 0) psest[i] = 0;
						}
						else{ 
							if (q < (R -1)){
								for (j = 0; j < *nt; j++)
									if (fabs(x[i] - t[j]) < *bw)
										psest[i] += kernel((x[i] - t[j])/(*bw), *nkernel) * w[j];
							}
							else{
								for (j = 0; j < *nt; j++)
									if (fabs(x[i] - t[j]) < *bw)
										psest[i] += kernelboundary(-(x[i] - t[j])/(*bw), R - q, *nkernel) * w[j];
								if (psest[i] < 0) psest[i] = 0;
							}
						}
						psest[i] = psest[i]/ (*bw);
					}
					break;
				default:
					break;
			}
			break;
		}
// h
	    case 4:{
			switch(*nbound){
				case 1:
					for (i = 0; i < *nx; i++){
						psest[i] = 0;		
						for (j = 0; j < *nt; j++)
							if (fabs(x[i] - t[j]) < *bw)
								psest[i] += kernel((x[i] - t[j])/(*bw), *nkernel) * phat[j] / (*nt - j);
						psest[i] = psest[i] / (*bw);
					}
					break;
				case 2:
					for (i = 0; i < *nx; i++){
						psest[i] = 0;
						q = x[i]/(*bw);
						if (q < 1){
							for (j = 0; j < *nt; j++)
								if (fabs(x[i] - t[j]) < *bw)
									psest[i] += kernelboundary((x[i] - t[j])/(*bw), q, *nkernel) * phat[j] / (*nt - j);
							if (psest[i] < 0) psest[i] = 0;
						}
						else{
							for (j = 0; j < *nt; j++)
								if (fabs(x[i] - t[j]) < *bw)
									psest[i] += kernel((x[i] - t[j])/(*bw), *nkernel) * phat[j] / (*nt - j);
						}
						psest[i] = psest[i] / (*bw);
					}
					break;
				case 3:
					for (i = 0; i < *nx; i++){
						psest[i] = 0;
						q = x[i]/(*bw);
						if (q < (R -1)){
							for (j = 0; j < *nt; j++)
								if (fabs(x[i] - t[j]) < *bw)
									psest[i] += kernel((x[i] - t[j])/(*bw), *nkernel) * phat[j] / (*nt - j);
						}
						else{
							for (j = 0; j < *nt; j++)
								if (fabs(x[i] - t[j]) < *bw)
									psest[i] += kernelboundary(-(x[i] - t[j])/(*bw), R - q, *nkernel) * phat[j] / (*nt - j);
							if (psest[i] < 0) psest[i] = 0;
						}
						psest[i] = psest[i] / (*bw);
					}
					break;
				case 4:
					for (i = 0; i < *nx; i++){
						psest[i] = 0;
						q = x[i]/(*bw);
						if (q < 1){
							for (j = 0; j < *nt; j++)
								if (fabs(x[i] - t[j]) < *bw)
									psest[i] += kernelboundary((x[i] - t[j])/(*bw), q, *nkernel) * phat[j] / (*nt - j);
								if (psest[i] < 0) psest[i] = 0;
						}
						else{ 
							if (q < (R -1)){
								for (j = 0; j < *nt; j++)
									if (fabs(x[i] - t[j]) < *bw)
										psest[i] += kernel((x[i] - t[j])/(*bw), *nkernel) * phat[j] / (*nt - j);
							}
							else{
								for (j = 0; j < *nt; j++)
									if (fabs(x[i] - t[j]) < *bw)
										psest[i] += kernelboundary(-(x[i] - t[j])/(*bw), R - q, *nkernel) * phat[j] / (*nt - j);
								if (psest[i] < 0) psest[i] = 0;
							}
						}
						psest[i] = psest[i] / (*bw);
					}
					break;
				default:
					break;
			}
			break;
		default:
			break;
		}
	}
	free(w);
}

