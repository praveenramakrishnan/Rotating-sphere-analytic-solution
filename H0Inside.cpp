dcomplex H0InsideTheta(double N,double kprop, double jnNkaBessel_m, dcomplex cnm0,dcomplex dnm0,dcomplex alphanm0, dcomplex betanm0, double XnmcTheta0, double curl02Theta0)
{
	dcomplex H0in1 = -cnm0*alphanm0*jnNkaBessel_m;
	dcomplex H0in2 = (1.0*j)/(kprop*N)*dnm0*betanm0;
	
	return H0in1*XnmcTheta0 + H0in2*curl02Theta0;
}

dcomplex H0InsidePhi(double N,double kprop, double jnNkaBessel_m, dcomplex cnm0,dcomplex dnm0,dcomplex alphanm0, dcomplex betanm0, double XnmcPhi0, double curl02Phi0)
{
	dcomplex H0in1 = -cnm0*alphanm0*jnNkaBessel_m;
	dcomplex H0in2 = (1.0*j)/(kprop*N)*dnm0*betanm0;

	return H0in1*XnmcPhi0 + H0in2*curl02Phi0;
}

dcomplex H0InsideRadial(double N,double kprop, dcomplex dnm0,dcomplex betanm0, double curl02Radial0)
{
	dcomplex H0in2 = (1.0*j)/(kprop*N)*dnm0*betanm0;

	return H0in2*curl02Radial0;
}
