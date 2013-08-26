#include "survPresmooth.h"

double auxfunplugin(double x, int nkernel);

void funplugin(double *x1, double *x2, int *n, int *nkernel, double *const1, double *const2, int *nestimand, double *int1, double *int2, double *int3, double *int4, double *int5, double *result){
  double temp;
  
  if(*nestimand == 3)
    temp = auxfunplugin(*x1 / *x2, *nkernel) / *x2;
  else
    temp = auxfunplugin(*x2 / *x1, *nkernel) / *x1;
  *result = pow(*const1, 2) / 4.0 * (*int1 * pow(*x2, 4) + 2 * *int2 * pow(*x1, 2) * pow(*x2, 2) + *int3 * pow(*x1, 4)) + (*const2 * *int4 / *x2 + temp * *int5) / *n;
}

double auxfunplugin(double x, int nkernel){
  switch(nkernel){
  case 1:
    if(x < 1) return - 10.0/15827.0*pow(x, 9) + 15.0/1547.0*pow(x, 7) - 15.0/143.0*pow(x, 5) + 10.0/49.0*pow(x, 4) - 15.0/49.0*pow(x, 2) + 5.0/7.0;
    else return 5.0/2263261.0 * (323323.0/x - 138567.0/pow(x, 3) + 92378.0/pow(x, 5) - 47481.0/pow(x, 6) + 4389.0/pow(x, 8) - 286.0/pow(x, 10));
    break;
  case 2:
    if(x < 1) return - 35.0/799227.0*pow(x, 13) + 35.0/43263.0*pow(x, 11) - 175.0/22287.0*pow(x, 9) + 70.0/969.0*pow(x, 7) - 175.0/1287.0*pow(x, 6) + 175.0/891.0*pow(x, 4) - 35.0/99.0*pow(x, 2) + 350.0/429.0;
    else return 350.0/429.0/x - 35.0/99.0/pow(x, 3) + 175.0/891.0/pow(x, 5) - 175.0/1287.0/pow(x, 7) + 70.0/969.0/pow(x, 8) - 175.0/22287.0/pow(x, 10) + 35.0/43263.0/pow(x, 12) - 35.0/799227.0/pow(x, 14);
    break;
  default:
    return 0;
    break;
  }
}


