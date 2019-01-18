dcomplex E0ScattTheta(double kprop, dcomplex h2nkaHankel_m, dcomplex anm0, dcomplex bnm0, dcomplex alphanm0, dcomplex betanm0, double XnmsTheta0, dcomplex curl03Theta0 )
{
	dcomplex E0sc1 = -bnm0*betanm0*h2nkaHankel_m;
	dcomplex E0sc2 = -j/kprop*anm0*alphanm0;
	
	return E0sc1*XnmsTheta0 + E0sc2*curl03Theta0;
}

dcomplex E0ScattPhi(double kprop, dcomplex h2nkaHankel_m, dcomplex anm0, dcomplex bnm0, dcomplex alphanm0, dcomplex betanm0, double XnmsPhi0, dcomplex curl03Phi0)
{
	dcomplex E0sc1 = -bnm0*betanm0*h2nkaHankel_m;
	dcomplex E0sc2 = -j/kprop*anm0*alphanm0;

	return E0sc1*XnmsPhi0 + E0sc2*curl03Phi0;
}

dcomplex E0ScattRadial(double kprop, dcomplex anm0,dcomplex alphanm0, dcomplex curl03Radial0)
{
	dcomplex E0sc2 = -j/kprop*anm0*alphanm0;

	return E0sc2*curl03Radial0;
}
