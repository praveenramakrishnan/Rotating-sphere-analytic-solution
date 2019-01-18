dcomplex F0ScattTheta(int n, double m,dcomplex anm0,dcomplex bnm0,dcomplex anmPrime,dcomplex bnmPrime,double epsilon, double mu, double kprop, double theta_m, double phi_m)
{
	
	dcomplex T0Theta = cos(m*phi_m)*(anmPrime*PlmPrime(n,m,theta_m) + bnmPrime*gsl_sf_legendre_Plm(n , m , cos(theta_m))*1.0*m/sin(theta_m));
	dcomplex F0Theta = j/kprop*mu*epsilon*T0Theta;
	return F0Theta;
	
}

dcomplex F0ScattPhi(int n, int m,dcomplex anm0,dcomplex bnm0,dcomplex anmPrime,dcomplex bnmPrime,double epsilon, double mu, double kprop, double theta_m, double phi_m)
{
	dcomplex T0Phi = -sin(m*phi_m)*(m*1.0*anmPrime*gsl_sf_legendre_Plm(n , m , cos(theta_m))/sin(theta_m) + bnmPrime*PlmPrime(n,m,theta_m));
	dcomplex F0Phi = j/kprop*epsilon*mu*1.0*T0Phi;
	return F0Phi;
}

void G0Scatt(int nMax,double N, double mu_r, double epsilon_r,double kprop, double Rsphere, double theta_i,double theta_m,double phi_m)
{
	
	dcomplex F0Theta = 0.0;
	dcomplex F0Phi = 0.0;
	
	for(int n=1;n<nMax;n++)
	{
		//Special functions
		double jnNkaBessel = gsl_sf_bessel_jl(n, N*kprop*Rsphere);
		double jnkaBessel = gsl_sf_bessel_jl(n, kprop*Rsphere);
		double DjnNkaBessel = Djn(n, N*kprop*Rsphere);
		double DjnkaBessel = Djn(n, kprop*Rsphere);
		
		dcomplex h2nkaHankel =  h2n(n, kprop*Rsphere);
		dcomplex Dh2nkaHankel =  Dh2n(n, kprop*Rsphere);
		
		for(int m=0;m<n+1;m++)
		{
		//Coefficients
			
			dcomplex anm0 = anm(N,mu_r,epsilon_r, h2nkaHankel,jnNkaBessel,jnkaBessel,Dh2nkaHankel,DjnNkaBessel,DjnkaBessel);
			dcomplex bnm0 = bnm(N,mu_r,epsilon_r, h2nkaHankel,jnNkaBessel,jnkaBessel,Dh2nkaHankel,DjnNkaBessel,DjnkaBessel);
		
			dcomplex anmPrime = anm0*PlmPrime(n,m,theta_i);
			dcomplex bnmPrime = m*1.0*bnm0*gsl_sf_legendre_Plm(n , m , cos(theta_i))/sin(theta_i);
			double epsilon=0.0;
			if(m==0)
			{
				epsilon = 1.0;
			}
			else
			{
				epsilon = 2.0;
			}
			double mu = (2.0*n+1.0)*fact(n-m)*1.0/(fact(n+m)*n*(n+1.0));
			
			F0Theta += F0ScattTheta(n,m,anm0,bnm0,anmPrime,bnmPrime,epsilon,mu,kprop,theta_m,phi_m);
			F0Phi += F0ScattPhi(n,m,anm0,bnm0,anmPrime,bnmPrime,epsilon,mu,kprop,theta_m,phi_m);
		}
	}
	dcomplex G0 = (abs(F0Theta)*abs(F0Theta) + abs(F0Phi)*abs(F0Phi))*4.0/(Rsphere*Rsphere);
	cout << G0 << endl;
}



