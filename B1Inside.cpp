dcomplex B1InsideTheta(double kprop,double N,dcomplex C1, dcomplex cnm1,dcomplex dnm1,  dcomplex gnm1,dcomplex pnm1, dcomplex qnm1,dcomplex Dfnm1, double jnNkaBessel_m, double DjnNkaBessel_m, double XnmsTheta0, double XnmcPhi0)
{
	dcomplex B1in1 = kprop*kprop*N*N*(gnm1 + jnNkaBessel_m*qnm1 + (1.0*j)*cnm1*C1/(N*kprop)*jnNkaBessel_m);
	dcomplex B1in2 = Dfnm1 + DjnNkaBessel_m*pnm1 + C1*dnm1*DjnNkaBessel_m;
	return B1in1*XnmsTheta0 - B1in2*XnmcPhi0;
}

dcomplex B1InsidePhi(double kprop,double N,dcomplex C1, dcomplex cnm1,dcomplex dnm1,  dcomplex gnm1,dcomplex pnm1, dcomplex qnm1,dcomplex Dfnm1, double jnNkaBessel_m, double DjnNkaBessel_m, double XnmsPhi0, double XnmcTheta0)
{
	dcomplex B1in1 = kprop*kprop*N*N*(gnm1 + jnNkaBessel_m*qnm1 + (1.0*j)*cnm1*C1/(N*kprop)*jnNkaBessel_m);
	dcomplex B1in2 = Dfnm1 + DjnNkaBessel_m*pnm1 + C1*dnm1*DjnNkaBessel_m;
	return B1in1*XnmsPhi0 + B1in2*XnmcTheta0;
}

dcomplex B1InsideRadial(int n, dcomplex C1, dcomplex dnm1, dcomplex fnm1, dcomplex pnm1, double jnNkaBessel_m, double Ynmc0, double R_m)
{
	dcomplex B1in3 = fnm1 + jnNkaBessel_m*pnm1 + C1*dnm1*jnNkaBessel_m;
	return n*(n+1.0)/R_m*B1in3*Ynmc0;
}
