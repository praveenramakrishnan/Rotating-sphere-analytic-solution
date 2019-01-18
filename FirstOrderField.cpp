void TotalField(int nMax,double N, double mu_r, double epsilon_r,double kprop, double Rsphere, double theta_i, double R_m, double theta_m, double phi_m, dcomplex C1, double beta_a, double K01, string out_filename_Efield, string out_filename_Hfield)
{
/*
DATE:  15/03/2016
AUTHOR: PRAVEEN K R
This function prints the E and H field for the rotating dielectric sphere in the output files in the specified format.
First the E and H field are calculated in spherical co-ordinate and then converted into rectangular co-ordinates.
The zero-order and first-order fields are calculated separately and then combined.
In De Zutter's paper, in the first order case, the formulas are for D and B fields.
This is converted to E and H fields using constitutive relations.
Also, for points ouside sphere, the incident field is added to this to get the total field.
*/

	//Zero order E and H fields
	dcomplex E0Radial = 0.0;
	dcomplex E0Theta = 0.0;
	dcomplex E0Phi = 0.0;
	dcomplex H0Radial = 0.0;
	dcomplex H0Theta = 0.0;
	dcomplex H0Phi = 0.0;

	//First order D and B fields
	dcomplex D1Radial = 0.0;
	dcomplex D1Theta = 0.0;
	dcomplex D1Phi = 0.0;
	dcomplex B1Radial = 0.0;
	dcomplex B1Theta = 0.0;
	dcomplex B1Phi = 0.0;

	//First order E and H fields
	dcomplex E1Radial = 0.0;
	dcomplex E1Theta = 0.0;
	dcomplex E1Phi = 0.0;
	dcomplex H1Radial = 0.0;
	dcomplex H1Theta = 0.0;
	dcomplex H1Phi = 0.0;

	//E and H fields by combining zero and first order fields
	dcomplex ETotalRadial = 0.0;
	dcomplex ETotalTheta = 0.0;
	dcomplex ETotalPhi = 0.0;
	dcomplex HTotalRadial = 0.0;
	dcomplex HTotalTheta = 0.0;
	dcomplex HTotalPhi = 0.0;

	//E nad H fields converted to rectangular co-ordinates
	dcomplex ETx, ETy, ETz, HTx, HTy, HTz;



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

		dcomplex tau1_s = tau(n, N, kprop, C1, Rsphere);
		dcomplex Dtau1_s = Dtau( n, N, kprop, C1, Rsphere);
		dcomplex tau1_m = tau(n, N, kprop, C1, R_m);
		dcomplex Dtau1_m = Dtau( n, N, kprop, C1, R_m);

		for(int m=0;m<n+1;m++)
		{
		//Coefficients
			dcomplex knm0 = Knm(n,m);
			dcomplex alphanm0 = alphanm(n,m,theta_i,knm0);
			dcomplex betanm0 = betanm(n,m,theta_i,knm0);

			dcomplex cnm1 = cnmF(m, alphanm0, cnm0 );
			dcomplex dnm1 = dnmF(m, betanm0, dnm0 );

			dcomplex fnm1_s = fnm(dnm1, tau1_s);
			dcomplex gnm1_s = gnm(cnm1, tau1_s, N, kprop);
			dcomplex Dfnm1_s = Dfnm(dnm1, Dtau1_s);
			dcomplex Dgnm1_s = Dgnm(cnm1, Dtau1_s, N, kprop);
			dcomplex Anm1_s = Anm(n, m, N, epsilon_r, mu_r, kprop, Rsphere, theta_i);
			dcomplex Cnm1_s = Cnm(n, m, N,epsilon_r, mu_r, kprop, Rsphere, theta_i);

			//p,q,gamma,delta coefficients
			dcomplex pnm1 = pnm( N, epsilon_r, kprop, C1, Cnm1_s, h2nkaHankel_s, jnNkaBessel_s, dnm1, fnm1_s, Dh2nkaHankel_s, DjnNkaBessel_s, Dfnm1_s);
			dcomplex qnm1 = qnm(N, epsilon_r, C1, Anm1_s, h2nkaHankel_s, jnNkaBessel_s, gnm1_s, Dh2nkaHankel_s, DjnNkaBessel_s,  Dgnm1_s);
			dcomplex gammanm1 = gammanm(N, epsilon_r, kprop, C1, Cnm1_s, h2nkaHankel_s, jnNkaBessel_s, dnm1, fnm1_s, Dh2nkaHankel_s, DjnNkaBessel_s, Dfnm1_s);
			dcomplex deltanm1 = deltanm( N, epsilon_r, C1, Anm1_s, h2nkaHankel_s, jnNkaBessel_s, gnm1_s, Dh2nkaHankel_s, DjnNkaBessel_s, Dgnm1_s);

			dcomplex fnm1_m = fnm( dnm1, tau1_m);
			dcomplex gnm1_m = gnm(cnm1, tau1_m, N, kprop);
			dcomplex Dfnm1_m = Dfnm(dnm1, Dtau1_m);
			dcomplex Dgnm1_m = Dgnm(cnm1, Dtau1_m, N, kprop);

		//Ynm and Xnm	
			double Ynmc0 = cos(m*phi_m)*gsl_sf_legendre_Plm(n , m , cos(theta_m));
			double Ynms0 = sin(m*phi_m)*gsl_sf_legendre_Plm(n , m , cos(theta_m));
			
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

			   //First Order Electric flux density
				//Radial component
				D1Radial += D1InsideRadial(n, gnm1_m, qnm1, jnNkaBessel_m, Ynms0, R_m);
				//Theta component (polar)
				D1Theta += D1InsideTheta(pnm1, qnm1, fnm1_m, Dgnm1_m, jnNkaBessel_m, DjnNkaBessel_m, XnmcTheta0, XnmsPhi0);
				//Phi component (azimuthal)
				D1Phi += D1InsidePhi(pnm1, qnm1, fnm1_m, Dgnm1_m, jnNkaBessel_m, DjnNkaBessel_m, XnmcPhi0, XnmsTheta0);
				

			//Magnetic field insdie the sphere

			   //Zero Order Magnetic field
				//Radial component
				H0Radial += H0InsideRadial(N, kprop, dnm0, betanm0, curl02Radial0);
				//Theta component (polar)
				H0Theta += H0InsideTheta(N, kprop, jnNkaBessel_m, cnm0, dnm0, alphanm0, betanm0, XnmcTheta0, curl02Theta0);
				//Phi component (azimuthal)
				H0Phi += H0InsidePhi(N, kprop, jnNkaBessel_m, cnm0, dnm0, alphanm0, betanm0, XnmcPhi0, curl02Phi0);
			   //First Order Magnetic flux density
				//Radial component
				B1Radial += B1InsideRadial(n, C1, dnm1, fnm1_m, pnm1, jnNkaBessel_m, Ynmc0, R_m);
				//Theta component (polar)
				B1Theta += B1InsideTheta(kprop, N, C1, cnm1, dnm1, gnm1_m, pnm1, qnm1, Dfnm1_m, jnNkaBessel_m, DjnNkaBessel_m, XnmsTheta0,XnmcPhi0);
				//Phi component (azimuthal)
				B1Phi += B1InsidePhi(kprop,N,C1, cnm1, dnm1, gnm1_m, pnm1, qnm1, Dfnm1_m, jnNkaBessel_m, DjnNkaBessel_m, XnmsPhi0, XnmcTheta0);
			}
			else
			{
			//Electric field scattered
			
			  //Zero Order Electric field
				//Radial component 
				E0Radial += E0ScattRadial(kprop, anm0, alphanm0, curl03Radial0);
				//Theta component (polar)
				E0Theta += E0ScattTheta(kprop, h2nkaHankel_m, anm0, bnm0, alphanm0, betanm0, XnmsTheta0, curl03Theta0);
				//Phi component (azimuthal)
				E0Phi += E0ScattPhi(kprop, h2nkaHankel_m, anm0, bnm0, alphanm0, betanm0, XnmsPhi0, curl03Phi0);
			  //First Order Electric flux density
				//Radial component
				D1Radial += D1ScattRadial(n, deltanm1, h2nkaHankel_m, Ynms0, R_m);
				//Theta component (polar)
				D1Theta += D1ScattTheta(gammanm1,deltanm1, h2nkaHankel_m, Dh2nkaHankel_m, XnmcTheta0, XnmsPhi0);
				//Phi component (azimuthal)
				D1Phi += D1ScattPhi(gammanm1,deltanm1, h2nkaHankel_m, Dh2nkaHankel_m, XnmcPhi0, XnmsTheta0);

			//Magnetic field scattered
			
			   //Zero Order Magnetic field
				//Radial component
				H0Radial += H0ScattRadial(kprop, bnm0, betanm0, curl04Radial0);
				//Theta component (polar)
				H0Theta += H0ScattTheta(kprop, h2nkaHankel_m, anm0, bnm0, alphanm0, betanm0, XnmcTheta0, curl04Theta0);
				//Phi component (azimuthal)
				H0Phi += H0ScattPhi(kprop, h2nkaHankel_m, anm0, bnm0, alphanm0, betanm0, XnmcPhi0, curl04Phi0);
			   //First Order Magnetic flux density
				//Radial component
				B1Radial += B1ScattRadial(n, gammanm1, h2nkaHankel_m, Ynmc0, R_m);
				//Theta component (polar)
				B1Theta += B1ScattTheta(kprop, deltanm1, gammanm1, h2nkaHankel_m, Dh2nkaHankel_m, XnmsTheta0, XnmcPhi0);
				//Phi component (azimuthal)
				B1Phi += B1ScattPhi(kprop, deltanm1, gammanm1, h2nkaHankel_m, Dh2nkaHankel_m, XnmsPhi0, XnmcTheta0);
			}
		}
	}
	
	//Adjust the magnitudes of magnetic fields with impedances
	double eps_in = epsilon_r*epsilon_0;
	double mu_in = mu_0*mu_r;
	double Rc0 = sqrt(mu_in/eps_in);
	dcomplex Rc1 = -(1.0*j)*kprop*C0*eps_in;
	double Z0 = sqrt(mu_0/epsilon_0);
	dcomplex Z1 = -(1.0*j)*kprop*C0*epsilon_0;
	if(R_m <= Rsphere)
	{
		H0Radial /= Rc0;
		H0Theta /= Rc0;
		H0Phi /= Rc0;
		B1Radial /= Rc1;
		B1Theta /= Rc1;
		B1Phi /= Rc1;
	}
	else
	{
		H0Radial /= Z0;
		H0Theta /= Z0;
		H0Phi /= Z0;
		B1Radial /= Z1;
		B1Theta /= Z1;
		B1Phi /= Z1;
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
	//Calculate E1,H1 from D1,B1
	double Const = K01*beta_a*C0;
	if(R_m <= Rsphere)
	{
		E1Radial = (D1Radial - Const*(-H0Theta*sin(theta_m)))/eps_in;
		E1Theta = (D1Theta - Const*H0Radial*sin(theta_m))/eps_in;
		E1Phi = (D1Phi)/eps_in;
		H1Radial = (B1Radial + Const*(-E0Theta*sin(theta_m)))/mu_in;
		H1Theta = (B1Theta + Const*E0Radial*sin(theta_m))/mu_in;
		H1Phi = (B1Phi)/mu_in;
	}
	else
	{
		E1Radial = D1Radial/epsilon_0;
		E1Theta = D1Theta/epsilon_0;
		E1Phi = D1Phi/epsilon_0;
		H1Radial = B1Radial/mu_0;
		H1Theta = B1Theta/mu_0; 
		H1Phi = B1Phi/mu_0;
	}
	
	//Sum of 0 and 1st order fields
	ETotalRadial = E0Radial + beta_a*E1Radial;
	ETotalTheta = E0Theta + beta_a*E1Theta;
	ETotalPhi = E0Phi + beta_a*E1Phi;
	HTotalRadial = H0Radial + beta_a*H1Radial;
	HTotalTheta = H0Theta + beta_a*H1Theta;
	HTotalPhi = H0Phi + beta_a*H1Phi;

	//Conversion to rectangular co-ordinates
 	double xm = R_m*sin(theta_m)*cos(phi_m);
	double ym = R_m*sin(theta_m)*sin(phi_m);
	double zm = R_m*cos(theta_m);

	//The fields along rectangular co-ordinates
	  //Total
		ETx = ETotalRadial*sin(theta_m)*cos(phi_m) + ETotalTheta*cos(theta_m)*cos(phi_m) - ETotalPhi*sin(phi_m);
		ETy = ETotalRadial*sin(theta_m)*sin(phi_m) + ETotalTheta*cos(theta_m)*sin(phi_m) + ETotalPhi*cos(phi_m);
		ETz = ETotalRadial*cos(theta_m) - ETotalTheta*sin(theta_m);
		HTx = HTotalRadial*sin(theta_m)*cos(phi_m) + HTotalTheta*cos(theta_m)*cos(phi_m) - HTotalPhi*sin(phi_m);
		HTy = HTotalRadial*sin(theta_m)*sin(phi_m) + HTotalTheta*cos(theta_m)*sin(phi_m) + HTotalPhi*cos(phi_m);
		HTz = HTotalRadial*cos(theta_m) - HTotalTheta*sin(theta_m);

//Print to output files
        ofstream file1;
	file1.open(out_filename_Efield.c_str(),ios::app);
	file1 << scientific << setprecision(12) << xm << "\t" << ym << "\t" << zm << "\t" << real(ETx) << "\t" << imag(ETx) << "\t" << abs(ETx) << "\t" << atan2(imag(ETx),real(ETx)) << "\t" << real(ETy) << "\t" << imag(ETy) << "\t" << abs(ETy) << "\t" << atan2(imag(ETy),real(ETy)) << "\t" << real(ETz) << "\t" << imag(ETz) << "\t" << abs(ETz) << "\t" << atan2(imag(ETz),real(ETz)) << endl;

	ofstream file2;
	file2.open(out_filename_Hfield.c_str(),ios::app);
	file2 << scientific << setprecision(12) << xm << "\t" << ym << "\t" << zm << "\t" << real(HTx) << "\t" << imag(HTx) << "\t" << abs(HTx) << "\t" << atan2(imag(HTx),real(HTx)) << "\t" << real(HTy) << "\t" << imag(HTy) << "\t" << abs(HTy) << "\t" << atan2(imag(HTy),real(HTy)) << "\t" << real(HTz) << "\t" << imag(HTz) << "\t" << abs(HTz) << "\t" << atan2(imag(HTz),real(HTz))  << endl;
}
