dcomplex B1ScattTheta(double kprop, dcomplex deltanm1, dcomplex gammanm1, dcomplex h2nkaHankel_m,dcomplex Dh2nkaHankel_m, double XnmsTheta0, double XnmcPhi0)
{
	dcomplex B1in1 = kprop*kprop*h2nkaHankel_m*deltanm1;
	dcomplex B1in2 = Dh2nkaHankel_m*gammanm1;
	return B1in1*XnmsTheta0 - B1in2*XnmcPhi0;
}

dcomplex B1ScattPhi(double kprop, dcomplex deltanm1, dcomplex gammanm1, dcomplex h2nkaHankel_m,dcomplex Dh2nkaHankel_m, double XnmsPhi0, double XnmcTheta0)
{
	dcomplex B1in1 = kprop*kprop*h2nkaHankel_m*deltanm1;
	dcomplex B1in2 = Dh2nkaHankel_m*gammanm1;
	return B1in1*XnmsPhi0 + B1in2*XnmcTheta0;
}

dcomplex B1ScattRadial(int n, dcomplex gammanm1, dcomplex h2nkaHankel_m, double Ynmc0, double R_m)
{
	dcomplex B1in3 = h2nkaHankel_m*gammanm1;
	return n*(n+1.0)/R_m*B1in3*Ynmc0;
}
