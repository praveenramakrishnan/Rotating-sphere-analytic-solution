dcomplex alphanm(int l, int m, double theta_i, dcomplex knm0)
{

	return knm0*PlmPrime(l,m,theta_i);
	
}

dcomplex betanm(int l, int m, double theta_i, dcomplex knm0)
{
	return m*gsl_sf_legendre_Plm(l,m,cos(theta_i))*knm0/sin(theta_i);
	
}


long double fact(long double n)
{
	assert(n>=0);
	if(n == 0)
	{
		return 1.0;
	}
	else
	{
		return n*fact(n-1);
	}
}
dcomplex Knm(int n,int m)
{
	assert(n>0 && m>=0);
	double epsilon=0;
	if(m==0)
	{
		epsilon = 1.0;
	}
	else
	{
		epsilon = 2.0;
	}
	double mu = (2.0*n+1.0)*fact(n-m)*1.0/fact(n+m);
	return -pow(-j,n)*epsilon*mu*1.0/(n*(n+1.0));	

}
