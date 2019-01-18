dcomplex E0InsideTheta(double N,double kprop, double jnNkaBessel_m, dcomplex cnm0,dcomplex dnm0,dcomplex alphanm0, dcomplex betanm0, double XnmsTheta0, double curl01Theta0)
{
	dcomplex E0in1 = dnm0*betanm0*jnNkaBessel_m;
	dcomplex E0in2 = j/(kprop*N)*cnm0*alphanm0;
	
	return E0in1*XnmsTheta0 + E0in2*curl01Theta0;
}

dcomplex E0InsidePhi(double N,double kprop, double jnNkaBessel_m, dcomplex cnm0,dcomplex dnm0,dcomplex alphanm0, dcomplex betanm0, double XnmsPhi0, double curl01Phi0)
{
	dcomplex E0in1 = dnm0*betanm0*jnNkaBessel_m;
	dcomplex E0in2 = j/(kprop*N)*cnm0*alphanm0;

	return E0in1*XnmsPhi0+ E0in2*curl01Phi0;
}

dcomplex E0InsideRadial(double N,double kprop, dcomplex cnm0,dcomplex alphanm0, double curl01Radial0)
{
	dcomplex E0in2 = j/(kprop*N)*cnm0*alphanm0;

	return E0in2*curl01Radial0;
}
