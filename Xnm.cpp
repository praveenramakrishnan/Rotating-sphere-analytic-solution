double XnmsTheta(int n, int m, double theta_m, double phi_m)
{
	return m*cos(m*phi_m)/sin(theta_m)*gsl_sf_legendre_Plm(n , m , cos(theta_m));
}
double XnmsPhi(int n, int m, double theta_m, double phi_m)
{
	return -sin(m*phi_m)*PlmPrime(n, m, theta_m);
}

double XnmcTheta(int n, int m, double theta_m, double phi_m)
{
	return -m*sin(m*phi_m)/sin(theta_m)*gsl_sf_legendre_Plm(n , m , cos(theta_m));
}

double XnmcPhi(int n, int m, double theta_m, double phi_m)
{
	return -cos(m*phi_m)*PlmPrime(n , m, theta_m);
}

double XnmsPhiPrime(int n, int m, double theta_m, double phi_m)
{
	return -sin(m*phi_m)*PlmDoublePrime(n , m, theta_m);
}
double XnmsThetaPrime(int n, int m, double theta_m, double phi_m)
{
	return -m*m*sin(m*phi_m)/sin(theta_m)*gsl_sf_legendre_Plm(n , m , cos(theta_m));
}

double XnmcPhiPrime(int n, int m, double theta_m, double phi_m)
{
	return -cos(m*phi_m)*PlmDoublePrime(n , m, theta_m);
}
double XnmcThetaPrime(int n, int m, double theta_m, double phi_m)
{
	return -m*m*cos(m*phi_m)/sin(theta_m)*gsl_sf_legendre_Plm(n , m , cos(theta_m));
}

