dcomplex D1ScattTheta(dcomplex gammanm1,dcomplex deltanm1,dcomplex h2nkaHankel_m,dcomplex Dh2nkaHankel_m,double XnmcTheta0,double XnmsPhi0)
{
	dcomplex D1in1 = h2nkaHankel_m*gammanm1;
	dcomplex D1in2 = Dh2nkaHankel_m*deltanm1;
	return D1in1*XnmcTheta0 - D1in2*XnmsPhi0;
}

dcomplex D1ScattPhi(dcomplex gammanm1,dcomplex deltanm1,dcomplex h2nkaHankel_m,dcomplex Dh2nkaHankel_m,double XnmcPhi0,double XnmsTheta0)
{
	dcomplex D1in1 = h2nkaHankel_m*gammanm1;
	dcomplex D1in2 = Dh2nkaHankel_m*deltanm1;
	return D1in1*XnmcPhi0 + D1in2*XnmsTheta0;
}

dcomplex D1ScattRadial(int n,dcomplex deltanm1, dcomplex h2nkaHankel_m, double Ynms0, double R_m)
{
	dcomplex D1in2 = h2nkaHankel_m*deltanm1;
	return n*(n+1.0)/R_m*D1in2*Ynms0;
}
