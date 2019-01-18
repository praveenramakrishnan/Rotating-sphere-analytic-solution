dcomplex cnmF(int m, dcomplex alphanm0,dcomplex cnm0 )
{
		return m*1.0*alphanm0*cnm0;	
}
dcomplex dnmF(int m, dcomplex betanm0,dcomplex dnm0 )
{
		return m*1.0*betanm0*dnm0;	
}

dcomplex tau(int n,double N,double kprop,dcomplex C1,double R_m)
{
	return N*kprop*C1*R_m*0.5*(gsl_sf_bessel_jl(n+1, N*kprop*R_m)-gsl_sf_bessel_jl(n-1, N*kprop*R_m));
}

dcomplex Dtau(int n,double N,double kprop,dcomplex C1,double R_m)
{
	dcomplex term1 = N*kprop*C1*0.5*(gsl_sf_bessel_jl(n+1, N*kprop*R_m)-gsl_sf_bessel_jl(n-1, N*kprop*R_m));
	dcomplex term2 = N*kprop*R_m*C1*0.5*(Djn(n+1, N*kprop*R_m)-Djn(n-1, N*kprop*R_m))*N*kprop;
	return term1 + term2;
}

dcomplex fnm(dcomplex dnm1, dcomplex tau1)
{
	return dnm1*tau1;
}

dcomplex gnm(dcomplex cnm1, dcomplex tau1, double N, double kprop)
{
	return -(1.0*j)*cnm1*tau1/(N*kprop);
}

dcomplex Dfnm(dcomplex dnm1, dcomplex Dtau1)
{
	return dnm1*Dtau1;
}

dcomplex Dgnm(dcomplex cnm1, dcomplex Dtau1, double N, double kprop)
{
	return -(1.0*j)*cnm1*Dtau1/(N*kprop);
}

double K1(int n, int m, double N, double kprop, double R_m)
{
	return (n-1.0)*(n-m)/((2.0*n-1.0))*gsl_sf_bessel_jl(n-1, N*kprop*R_m);//avoided division by m
}

double K2(int n, int m, double N, double kprop, double R_m)
{
	return (n+2.0)*(n+m+1.0)/((2.0*n+3.0))*gsl_sf_bessel_jl(n+1, N*kprop*R_m);//avoided division by m
}


dcomplex Anm(int n, int m, double N,double epsilon_r, double mu_r, double kprop, double Rsphere, double theta_i)
{
	double jnNkaBessel02 = gsl_sf_bessel_jl(n+1, N*kprop*Rsphere);
	double jnkaBessel02 = gsl_sf_bessel_jl(n+1, kprop*Rsphere);
	double DjnNkaBessel02 = Djn(n+1, N*kprop*Rsphere)*N*kprop;
	double DjnkaBessel02 = Djn(n+1, kprop*Rsphere)*kprop;	
	dcomplex h2nkaHankel02 =  h2n(n+1, kprop*Rsphere);
	dcomplex Dh2nkaHankel02 =  Dh2n(n+1, kprop*Rsphere)*kprop;
	dcomplex knm02 = Knm(n+1,m);
	dcomplex dnm02 = dnm(N, mu_r, epsilon_r, h2nkaHankel02, jnNkaBessel02, jnkaBessel02, Dh2nkaHankel02, DjnNkaBessel02, DjnkaBessel02);
	dcomplex betanm02 = betanm(n+1,m,theta_i,knm02);
	
	if(n==1 || n==m)
	{
		return -K2(n, m, N, kprop, Rsphere)*betanm02*dnm02;//since K1=0 for n=1 and n=m
	}
	else
	{
		double jnNkaBessel01 = gsl_sf_bessel_jl(n-1, N*kprop*Rsphere);
		double jnkaBessel01 = gsl_sf_bessel_jl(n-1, kprop*Rsphere);
		double DjnNkaBessel01 = Djn(n-1, N*kprop*Rsphere)*N*kprop;
		double DjnkaBessel01 = Djn(n-1, kprop*Rsphere)*kprop;	
		dcomplex h2nkaHankel01 =  h2n(n-1, kprop*Rsphere);
		dcomplex Dh2nkaHankel01 =  Dh2n(n-1, kprop*Rsphere)*kprop;
		dcomplex knm01 = Knm(n-1,m);
		dcomplex dnm01 = dnm(N, mu_r, epsilon_r, h2nkaHankel01, jnNkaBessel01, jnkaBessel01, Dh2nkaHankel01, DjnNkaBessel01, DjnkaBessel01);
		dcomplex betanm01 = betanm(n-1,m,theta_i,knm01);
		return K1(n, m, N, kprop, Rsphere)*betanm01*dnm01 - K2(n, m, N, kprop, Rsphere)*betanm02*dnm02;
	}	
}

dcomplex Cnm(int n, int m, double N,double epsilon_r, double mu_r, double kprop, double Rsphere, double theta_i)
{
	double jnNkaBessel02 = gsl_sf_bessel_jl(n+1, N*kprop*Rsphere);
	double jnkaBessel02 = gsl_sf_bessel_jl(n+1, kprop*Rsphere);
	double DjnNkaBessel02 = Djn(n+1, N*kprop*Rsphere)*N*kprop;
	double DjnkaBessel02 = Djn(n+1, kprop*Rsphere)*kprop;	
	dcomplex h2nkaHankel02 =  h2n(n+1, kprop*Rsphere);
	dcomplex Dh2nkaHankel02 =  Dh2n(n+1, kprop*Rsphere)*kprop;
	dcomplex knm02 = Knm(n+1,m);
	dcomplex cnm02 = cnm(N, mu_r, epsilon_r, h2nkaHankel02, jnNkaBessel02, jnkaBessel02, Dh2nkaHankel02, DjnNkaBessel02, DjnkaBessel02);
	dcomplex alphanm02 = alphanm(n+1,m,theta_i,knm02);

	if(n==1 || n==m)
	{
		return -K2(n, m, N, kprop, Rsphere)*alphanm02*cnm02;//since K1=0 for n=1 and n=m
	}
	else
	{
		double jnNkaBessel01 = gsl_sf_bessel_jl(n-1, N*kprop*Rsphere);
		double jnkaBessel01 = gsl_sf_bessel_jl(n-1, kprop*Rsphere);
		double DjnNkaBessel01 = Djn(n-1, N*kprop*Rsphere)*N*kprop;
		double DjnkaBessel01 = Djn(n-1, kprop*Rsphere)*kprop;	
		dcomplex h2nkaHankel01 =  h2n(n-1, kprop*Rsphere);
		dcomplex Dh2nkaHankel01 =  Dh2n(n-1, kprop*Rsphere)*kprop;
		dcomplex knm01 = Knm(n-1,m);
		dcomplex cnm01 = cnm(N, mu_r, epsilon_r, h2nkaHankel01, jnNkaBessel01, jnkaBessel01, Dh2nkaHankel01, DjnNkaBessel01, DjnkaBessel01);
		dcomplex alphanm01 = alphanm(n-1,m,theta_i,knm01);
		return K1(n, m, N, kprop, Rsphere)*alphanm01*cnm01 - K2(n, m, N, kprop, Rsphere)*alphanm02*cnm02;
	}
}
