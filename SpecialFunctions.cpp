dcomplex h2n(int n,double x)
{
	return gsl_sf_bessel_jl(n,x) - j*gsl_sf_bessel_yl(n,x);
}

double jnPrime(int n, double x)
{
	return n/x*gsl_sf_bessel_jl(n,x) - gsl_sf_bessel_jl(n+1,x);
}

double ynPrime(int n, double x)
{
	return n/x*gsl_sf_bessel_yl(n,x) - gsl_sf_bessel_yl(n+1,x);
}

dcomplex h2nPrime(int n, double x)
{
	return n/x*h2n(n,x) - h2n(n+1,x);
}
double PlmPrime(int l,int m,double theta)
{
	double result_array[l-m+1];
	double result_deriv_array[l-m+1];
	gsl_sf_legendre_Plm_deriv_array(l,m,cos(theta),result_array,result_deriv_array);
	return result_deriv_array[l-m]*(-sin(theta)); 
}
double PlmDoublePrime(int l,int m,double theta)
{
	long double delta = 0.0;
	if(theta == 0.0)
	{
		long double temp = 1.0/1.0e6;
		delta = temp - theta;
	}
	else
	{
		long double temp = theta+theta/1.0e6;
		delta = temp - theta;
	}
	return (PlmPrime(l, m, theta+delta)-PlmPrime(l, m, theta-delta))/(2.0*delta); 
}
double Djn(int n,double x)
{	
	return jnPrime(n,x) + gsl_sf_bessel_jl(n,x)/x;
}

dcomplex Dh2n(int n,double x)
{	
	return h2nPrime(n,x) + h2n(n,x)/x;
}
