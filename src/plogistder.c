#include "survPresmooth.h"

void plogistder(double *x, int *nx, double *coef, int *deriv, double *pder) {
  int i;

  switch (*deriv) {
  case 1:
    for (i = 0; i < *nx; i++)
      pder[i] = coef[1] * exp(coef[0] + coef[1] * x[i]) / pow(1 + exp(coef[0] + coef[1] * x[i]), 2);
    break;
  case 2:
    for (i = 0; i < *nx; i++)
      pder[i] = pow(coef[1], 2) * exp(coef[0] + coef[1] * x[i]) * (1 - exp(coef[0] + coef[1] * x[i])) / pow(1 + exp(coef[0] + coef[1] * x[i]), 3);
    break;
  case 3:
    for (i = 0; i < *nx; i++)
      pder[i] = pow(coef[1], 3) * exp(coef[0] + coef[1] * x[i]) * (1 - 4 * exp(coef[0] + coef[1] * x[i]) + pow(exp(coef[0] + coef[1] * x[i]), 2)) / pow(1 + exp(coef[0] + coef[1] * x[i]), 4);
    break;
  case 4:
    for (i = 0; i < *nx; i++)
      pder[i] = pow(coef[1], 4) * exp(coef[0] + coef[1] * x[i]) * (1 - 11 * exp(coef[0] + coef[1] * x[i]) + 11 * pow(exp(coef[0] + coef[1] * x[i]), 2) - pow(exp(coef[0] + coef[1] * x[i]), 3))/pow(1 + exp(coef[0] + coef[1] * x[i]), 5);
    break;
  default:
    break;
  }
}

