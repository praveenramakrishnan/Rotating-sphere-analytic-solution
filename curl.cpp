double curl01Theta(int n, int m,double N,double kprop, double R_m, double theta_m, double phi_m)
{
	return -XnmcPhi(n,m,theta_m,phi_m)*(jnPrime(n,kprop*N*R_m)*kprop*N+gsl_sf_bessel_jl(n,kprop*N*R_m)/R_m);
}
double curl01Phi(int n, int m,double N,double kprop, double R_m, double theta_m, double phi_m)
{
	return XnmcTheta(n,m,theta_m,phi_m)*(jnPrime(n,kprop*N*R_m)*kprop*N+gsl_sf_bessel_jl(n,kprop*N*R_m)/R_m);
}
double curl01Radial(int n, int m,double N,double kprop, double R_m, double theta_m, double phi_m)
{
	return (sin(theta_m)*XnmcPhiPrime(n,m,theta_m,phi_m) + cos(theta_m)*XnmcPhi(n,m,theta_m,phi_m) - XnmcThetaPrime(n,m,theta_m,phi_m) )/(R_m*sin(theta_m))*gsl_sf_bessel_jl(n,kprop*N*R_m);
}


double curl02Theta(int n, int m,double N,double kprop, double R_m, double theta_m, double phi_m)
{
	return -XnmsPhi(n,m,theta_m,phi_m)*(jnPrime(n,kprop*N*R_m)*kprop*N+gsl_sf_bessel_jl(n,kprop*N*R_m)/R_m);
}
double curl02Phi(int n, int m,double N,double kprop, double R_m, double theta_m, double phi_m)
{
	return XnmsTheta(n,m,theta_m,phi_m)*(jnPrime(n,kprop*N*R_m)*kprop*N+gsl_sf_bessel_jl(n,kprop*N*R_m)/R_m);
}
double curl02Radial(int n, int m,double N,double kprop, double R_m, double theta_m, double phi_m)
{
	return (sin(theta_m)*XnmsPhiPrime(n,m,theta_m,phi_m) + cos(theta_m)*XnmsPhi(n,m,theta_m,phi_m)  - XnmsThetaPrime(n,m,theta_m,phi_m))/(R_m*sin(theta_m))*gsl_sf_bessel_jl(n,kprop*N*R_m);
}



dcomplex curl03Theta(int n, int m,double kprop, double R_m, double theta_m, double phi_m)
{
	return -XnmcPhi(n,m,theta_m,phi_m)*(h2nPrime(n,kprop*R_m)*kprop + h2n(n,kprop*R_m)/R_m);
}
dcomplex curl03Phi(int n, int m,double kprop, double R_m, double theta_m, double phi_m)
{
	return XnmcTheta(n,m,theta_m,phi_m)*(h2nPrime(n,kprop*R_m)*kprop + h2n(n,kprop*R_m)/R_m);
}
dcomplex curl03Radial(int n, int m,double kprop, double R_m, double theta_m, double phi_m)
{
	return (sin(theta_m)*XnmcPhiPrime(n,m,theta_m,phi_m) + cos(theta_m)*XnmcPhi(n,m,theta_m,phi_m) - XnmcThetaPrime(n,m,theta_m,phi_m) )/(R_m*sin(theta_m))*h2n(n,kprop*R_m);
}



dcomplex curl04Theta(int n, int m,double kprop, double R_m, double theta_m, double phi_m)
{
	return -XnmsPhi(n,m,theta_m,phi_m)*(h2nPrime(n,kprop*R_m)*kprop + h2n(n,kprop*R_m)/R_m);
}
dcomplex curl04Phi(int n, int m,double kprop, double R_m, double theta_m, double phi_m)
{
	return XnmsTheta(n,m,theta_m,phi_m)*(h2nPrime(n,kprop*R_m)*kprop + h2n(n,kprop*R_m)/R_m);
}
dcomplex curl04Radial(int n, int m,double kprop, double R_m, double theta_m, double phi_m)
{
	return (sin(theta_m)*XnmsPhiPrime(n,m,theta_m,phi_m) + cos(theta_m)*XnmsPhi(n,m,theta_m,phi_m) - XnmsThetaPrime(n,m,theta_m,phi_m) )/(R_m*sin(theta_m))*h2n(n,kprop*R_m);
}
