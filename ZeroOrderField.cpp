void ZeroOrderField(int nMax,double N, double mu_r, double epsilon_r,double kprop, double Rsphere, double theta_i, double R_m, double theta_m, double phi_m, string out_filename_Efield, string out_filename_Hfield)
{
/*
DATE:  15/03/2016
AUTHOR: PRAVEEN K R
This function prints the E and H field for the stationary dielectric sphere in the output files in the specified format.
First the E and H field are calculated in spherical co-ordinate and then converted into rectangular co-ordinates.
*/
	//E and H in spherical co-ordinates
	dcomplex E0Radial = 0.0;
	dcomplex E0Theta = 0.0;
	dcomplex E0Phi = 0.0;
	dcomplex H0Radial = 0.0;
	dcomplex H0Theta = 0.0;
	dcomplex H0Phi = 0.0;

	//E and H in rectangular co-ordinates
	dcomplex Ex, Ey,  Ez, Hx, Hy, Hz;


	for(int n=1;n<nMax;n++)
	{
		//Special functions at the sphere surface
		double jnNkaBessel_s = gsl_sf_bessel_jl(n, N*kprop*Rsphere);
		double jnkaBessel_s = gsl_sf_bessel_jl(n, kprop*Rsphere);
		double DjnNkaBessel_s = Djn(n, N*kprop*Rsphere)*N*kprop;
		double DjnkaBessel_s = Djn(n, kprop*Rsphere)*kprop;
		dcomplex h2nkaHankel_s =  h2n(n, kprop*Rsphere);
		dcomplex Dh2nkaHankel_s =  Dh2n(n, kprop*Rsphere)*kprop;
		
		//Special functions at the measurement point
		double jnNkaBessel_m = gsl_sf_bessel_jl(n, N*kprop*R_m);
		double jnkaBessel_m = gsl_sf_bessel_jl(n, kprop*R_m);
		double DjnNkaBessel_m = Djn(n, N*kprop*R_m)*N*kprop;
		double DjnkaBessel_m = Djn(n, kprop*R_m)*kprop;
		dcomplex h2nkaHankel_m =  h2n(n, kprop*R_m);
		dcomplex Dh2nkaHankel_m =  Dh2n(n, kprop*R_m)*kprop;

		//abcd coefficients
		dcomplex anm0 = anm(N,mu_r,epsilon_r, h2nkaHankel_s,jnNkaBessel_s,jnkaBessel_s,Dh2nkaHankel_s,DjnNkaBessel_s,DjnkaBessel_s);
		dcomplex bnm0 = bnm(N,mu_r,epsilon_r, h2nkaHankel_s,jnNkaBessel_s,jnkaBessel_s,Dh2nkaHankel_s,DjnNkaBessel_s,DjnkaBessel_s);
		dcomplex cnm0 = cnm(N,mu_r,epsilon_r, h2nkaHankel_s,jnNkaBessel_s,jnkaBessel_s,Dh2nkaHankel_s,DjnNkaBessel_s,DjnkaBessel_s);
		dcomplex dnm0 = dnm(N,mu_r,epsilon_r, h2nkaHankel_s,jnNkaBessel_s,jnkaBessel_s,Dh2nkaHankel_s,DjnNkaBessel_s,DjnkaBessel_s);

		for(int m=0;m<n+1;m++)
		{
		//Coefficients
			dcomplex knm0 = Knm(n,m);
			dcomplex alphanm0 = alphanm(n,m,theta_i,knm0);
			dcomplex betanm0 = betanm(n,m,theta_i,knm0);
			
			double XnmcTheta0 = XnmcTheta(n, m, theta_m, phi_m);
			double XnmcPhi0 = XnmcPhi(n, m, theta_m, phi_m);
			double XnmsTheta0 = XnmsTheta(n, m, theta_m, phi_m);
			double XnmsPhi0 = XnmsPhi(n, m, theta_m, phi_m);

			double XnmcPhiPrime0 = XnmcPhiPrime(n, m, theta_m, phi_m);
			double XnmcThetaPrime0 = XnmcThetaPrime(n, m, theta_m, phi_m);
			double XnmsPhiPrime0 = XnmsPhiPrime(n, m, theta_m, phi_m);
			double XnmsThetaPrime0 = XnmsThetaPrime(n, m, theta_m, phi_m);

		//Curl
			double curl01Theta0 = curl01Theta(n,m,N,kprop,R_m,theta_m,phi_m);
			double curl01Phi0 = curl01Phi(n, m, N, kprop, R_m, theta_m, phi_m);
			double curl01Radial0 = curl01Radial(n, m, N, kprop, R_m, theta_m, phi_m);
			double curl02Theta0 = curl02Theta(n, m, N, kprop, R_m, theta_m, phi_m);
			double curl02Phi0 = curl02Phi(n, m, N, kprop, R_m, theta_m, phi_m);
			double curl02Radial0 = curl02Radial(n, m, N, kprop, R_m, theta_m, phi_m);
			dcomplex curl03Theta0 = curl03Theta(n, m, kprop, R_m, theta_m, phi_m);
			dcomplex curl03Phi0 = curl03Phi(n, m, kprop, R_m, theta_m, phi_m);
			dcomplex curl03Radial0 = curl03Radial(n, m, kprop,R_m, theta_m, phi_m);
			dcomplex curl04Theta0 = curl04Theta(n, m, kprop, R_m, theta_m, phi_m);
			dcomplex curl04Phi0 = curl04Phi(n, m, kprop, R_m, theta_m, phi_m);
			dcomplex curl04Radial0 = curl04Radial(n,m, kprop, R_m, theta_m, phi_m);

			if(R_m <= Rsphere)
			{
			//Electric field insdie the sphere
			
			   //Zero Order Electric field
				//Radial component
				E0Radial += E0InsideRadial(N, kprop, cnm0, alphanm0, curl01Radial0);	
				//Theta component (polar)
				E0Theta += E0InsideTheta(N, kprop, jnNkaBessel_m, cnm0, dnm0, alphanm0, betanm0,  XnmsTheta0, curl01Theta0);		
				//Phi component (azimuthal)
				E0Phi += E0InsidePhi(N, kprop, jnNkaBessel_m, cnm0, dnm0, alphanm0, betanm0,  XnmsPhi0, curl01Phi0);
				
			//Magnetic field insdie the sphere

			   //Zero Order Magnetic field
				//Radial component
				H0Radial += H0InsideRadial(N, kprop, dnm0, betanm0, curl02Radial0);
				//Theta component (polar)
				H0Theta += H0InsideTheta(N, kprop, jnNkaBessel_m, cnm0, dnm0, alphanm0, betanm0, XnmcTheta0, curl02Theta0);
				//Phi component (azimuthal)
				H0Phi += H0InsidePhi(N, kprop, jnNkaBessel_m, cnm0, dnm0, alphanm0, betanm0, XnmcPhi0, curl02Phi0);
			}
			else
			{
			//Electric field scattered
			
				//Radial component
				E0Radial += E0ScattRadial(kprop, anm0, alphanm0, curl03Radial0);
				//Theta component (polar)
				E0Theta += E0ScattTheta(kprop, h2nkaHankel_m, anm0, bnm0, alphanm0, betanm0, XnmsTheta0, curl03Theta0);
				//Phi component (azimuthal)
				E0Phi += E0ScattPhi(kprop, h2nkaHankel_m, anm0, bnm0, alphanm0, betanm0, XnmsPhi0, curl03Phi0);

			//Magnetic field scattered
			
				//Radial component
				H0Radial += H0ScattRadial(kprop, bnm0, betanm0, curl04Radial0);
				//Theta component (polar)
				H0Theta += H0ScattTheta(kprop, h2nkaHankel_m, anm0, bnm0, alphanm0, betanm0, XnmcTheta0, curl04Theta0);
				//Phi component (azimuthal)
				H0Phi += H0ScattPhi(kprop, h2nkaHankel_m, anm0, bnm0, alphanm0, betanm0, XnmcPhi0, curl04Phi0);

			}
		}
	}
	
	//Adjust the magnitudes of magnetic fields with impedances
	double eps_in = epsilon_r*epsilon_0;
	double mu_in = mu_0*mu_r;
	double Rc0 = sqrt(mu_in/eps_in);
	double Z0 = sqrt(mu_0/epsilon_0);
	if(R_m <= Rsphere)
	{
		H0Radial /= Rc0;
		H0Theta /= Rc0;
		H0Phi /= Rc0;
	}
	else
	{
		H0Radial /= Z0;
		H0Theta /= Z0;
		H0Phi /= Z0;
	}
	//Add incident field to get total field
	if(R_m > Rsphere)
	{
		double Psi = kprop*R_m*(cos(theta_m)*cos(theta_i) + sin(theta_m)*cos(phi_m)*sin(theta_i));
		dcomplex Ei = exp(-j*Psi);
		dcomplex Exi = -Ei*cos(theta_i);
		dcomplex Ezi = Ei*sin(theta_i);

		E0Radial += Ezi*cos(theta_m) + Exi*cos(phi_m)*sin(theta_m);
		E0Theta  += -Ezi*sin(theta_m) + Exi*cos(phi_m)*cos(theta_m);
		E0Phi += -Exi*sin(phi_m);
		
		dcomplex Hyi = -Ei/Z0;
		H0Radial += Hyi*sin(phi_m)*sin(theta_m);
		H0Theta += Hyi*sin(phi_m)*cos(theta_m);
		H0Phi += Hyi*cos(phi_m);
	}

	//Conversion into rectangular co-ordinates
  
 	double xm = R_m*sin(theta_m)*cos(phi_m);
	double ym = R_m*sin(theta_m)*sin(phi_m);
	double zm = R_m*cos(theta_m);
	Ex = E0Radial*sin(theta_m)*cos(phi_m) + E0Theta*cos(theta_m)*cos(phi_m) - E0Phi*sin(phi_m);
	Ey = E0Radial*sin(theta_m)*sin(phi_m) + E0Theta*cos(theta_m)*sin(phi_m) + E0Phi*cos(phi_m);
	Ez = E0Radial*cos(theta_m) - E0Theta*sin(theta_m);
	Hx = H0Radial*sin(theta_m)*cos(phi_m) + H0Theta*cos(theta_m)*cos(phi_m) - H0Phi*sin(phi_m);
	Hy = H0Radial*sin(theta_m)*sin(phi_m) + H0Theta*cos(theta_m)*sin(phi_m) + H0Phi*cos(phi_m);
	Hz = H0Radial*cos(theta_m) - H0Theta*sin(theta_m);

	//Print to output files
        ofstream file1;
	file1.open(out_filename_Efield.c_str(),ios::app);
	file1 << scientific << setprecision(12) << xm << "\t" << ym << "\t" << zm << "\t" << real(Ex) << "\t" << imag(Ex) << "\t" << abs(Ex) << "\t" << atan2(imag(Ex),real(Ex)) << "\t" << real(Ey) << "\t" << imag(Ey) << "\t" << abs(Ey) << "\t" << atan2(imag(Ey),real(Ey)) << "\t" << real(Ez) << "\t" << imag(Ez) << "\t" << abs(Ez) << "\t" << atan2(imag(Ez),real(Ez)) <<endl;
	
        ofstream file2;
	file2.open(out_filename_Hfield.c_str(),ios::app);
	file2 << scientific << setprecision(12) << xm << "\t" << ym << "\t" << zm << "\t" << real(Hx) << "\t" << imag(Hx) << "\t" << abs(Hx) << "\t" << atan2(imag(Hx),real(Hx)) << "\t" << real(Hy) << "\t" << imag(Hy) << "\t" << abs(Hy) << "\t" << atan2(imag(Hy),real(Hy)) << "\t" << real(Hz) << "\t" << imag(Hz) << "\t" << abs(Hz) << "\t" << atan2(imag(Hz),real(Hz)) <<endl;
}
