dcomplex D1InsideTheta(dcomplex pnm1,dcomplex qnm1, dcomplex fnm1, dcomplex Dgnm1, double jnNkaBessel_m, double DjnNkaBessel_m, double  XnmcTheta0, double XnmsPhi0)
{
	dcomplex D1in1 = fnm1 + pnm1*jnNkaBessel_m;
	dcomplex D1in2 = Dgnm1 + qnm1*DjnNkaBessel_m;
	return D1in1*XnmcTheta0 - D1in2*XnmsPhi0;
}

dcomplex D1InsidePhi(dcomplex pnm1,dcomplex qnm1, dcomplex fnm1, dcomplex Dgnm1, double jnNkaBessel_m, double DjnNkaBessel_m, double  XnmcPhi0, double XnmsTheta0)
{
	dcomplex D1in1 = fnm1 + pnm1*jnNkaBessel_m;
	dcomplex D1in2 = Dgnm1 + qnm1*DjnNkaBessel_m;
	return D1in1*XnmcPhi0 + D1in2*XnmsTheta0;
}

dcomplex D1InsideRadial(int n, dcomplex gnm1, dcomplex qnm1, double jnNkaBessel_m, double Ynms0, double R_m)
{
	dcomplex D1in3 = gnm1 + qnm1*jnNkaBessel_m;
	return n*1.0*(n+1.0)/R_m*D1in3*Ynms0;
}
