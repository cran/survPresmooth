double kernelder(double x, int nkernel, int deriv) {
	switch(nkernel){
	   	case 1:
			switch(deriv){
				case 0: 
					return 15.0/16.0 * (1.0 - x * x) * (1.0 - x * x);
					break;
				case 1: 
					return -15.0/4.0 * x * (1.0 - x * x);
					break;
				case 2: 
					return 15.0/4.0 * (3.0 * x * x - 1.0);
					break;
				default:
					break;
				}
		      	break;
	   	case 2:
			switch(deriv){
				case 0: 
					return 35.0/32.0 * (1.0 - x * x) * (1.0 - x * x) * (1.0 - x * x);
					break;
				case 1: 
					return -105.0/16.0 * x * (1.0 - x * x) * (1.0 - x * x) * (1.0 - x * x);
					break;
				case 2: 
					return 105.0/16.0 * (1.0 - x * x) * (1.0 - x * x) * (7.0 * x * x - 1.0);
					break;
				default:
					break;
			}
		      	break;
	   	default:
			switch(deriv){
				case 0: 
					return 15.0/16.0 * (1.0 - x * x) * (1.0 - x * x);
					break;
				case 1: 
					return -15.0/4.0 * x * (1.0 - x * x);
					break;
				case 2: 
					return 15.0/4.0 * (3.0 * x * x - 1.0);
					break;
				default:
					break;
			}
		      	break;
	}
}

