dcomplex H0ScattTheta(double kprop, dcomplex h2nkaHankel_m, dcomplex anm0,dcomplex bnm0,dcomplex alphanm0, dcomplex betanm0, double XnmcTheta0, dcomplex curl04Theta0)
{
	dcomplex H0sc1 = anm0*alphanm0*h2nkaHankel_m;
	dcomplex H0sc2 = -(1.0*j)/kprop*bnm0*betanm0;
	
	return H0sc1*XnmcTheta0 + H0sc2*curl04Theta0;
}

dcomplex H0ScattPhi(double kprop, dcomplex h2nkaHankel_m, dcomplex anm0,dcomplex bnm0,dcomplex alphanm0, dcomplex betanm0, double XnmcPhi0, dcomplex curl04Phi0)
{
	dcomplex H0sc1 = anm0*alphanm0*h2nkaHankel_m;
	dcomplex H0sc2 = -(1.0*j)/kprop*bnm0*betanm0;

	return H0sc1*XnmcPhi0 + H0sc2*curl04Phi0;
}

dcomplex H0ScattRadial(double kprop, dcomplex bnm0,dcomplex betanm0, dcomplex curl04Radial0)
{
	dcomplex H0sc2 = -(1.0*j)/kprop*bnm0*betanm0;

	return H0sc2*curl04Radial0;
}
