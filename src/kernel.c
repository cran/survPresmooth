double kernel(double x, int nkernel) {
  switch(nkernel){
  case 1:
    return 15.0/16.0 * (1.0 - x * x) * (1.0 - x * x);
    break;
  case 2:
    return 35.0/32.0 * (1.0 - x * x) * (1.0 - x * x) * (1.0 - x * x);
    break;
  default:
    return 15.0/16.0 * (1.0 - x * x) * (1.0 - x * x);
    break;
  }
}

