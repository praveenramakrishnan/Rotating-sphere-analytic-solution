/*
DATE:  15/03/2016
AUTHOR: PRAVEEN K R

This program is for the calculation of electric and magnetic fields E and H when a rotating sphere is illuminated by a plane wave.
This is an implementation of the paper titled "Scattering by a rotating dielectric sphere" by Daniel De Zutter,  publised in IEEE transactions in antennas and peopagation, vol.AP-28,no.5, September 1980.
Refer to Fig 1 of the paper for the configuration used. 
Equations 15 and 16 are used for calculating the zero order fields.
Equations 22 and 23 are used for calculating first order field.
In the end they are combined to calculate the total field.

The program takes as input the following parameters: mu_r, epsilon_r, Rsphere, frequency, beta_a, theta_i, in_filename_points, out_filename_Efield,  out_filename_Hfield

Format of input and output files:
1)in_filename_points contain the list of points at which the field is calculated
The file should contain the number of points in its first line followed by the lines containing the x, y, z co-ordinates of each point 
2)out_filename_Efield
The output of the calculated E field is written in this file
Each line of this file contains 15 entries which are x,y,z,real(Ex),Imag(Ex),mag(Ex),angle(Ex),real(Ey),Imag(Ey),mag(Ey),angle(Ey), real(Ez),Imag(Ez),mag(Ez),angle(Ez)
3)out_filename_Hfield
This file has same format as above with values realted to H field.

*Note that an error was detected in equation 26 of the paper and equation 5.67 of De Zutter's PhD thesis was used instead.*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <cassert>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>
#include <cstring>

using namespace std;
typedef complex<double> dcomplex;
#define j dcomplex(0.0,1.0)
const double pi = 3.1415926535897932;
const double C0 = 299792458.0;//speed of light
const double mu_0 = 4.0*pi*1.0e-7; 
const double epsilon_0 = 1.0/(mu_0*C0*C0);

#include "SpecialFunctions.cpp"
#include "coefficients.cpp"
#include "abcdConstants.cpp"
#include "coefficientsFirst.cpp"
#include "pqgadeConstants.cpp"
#include "Xnm.cpp"
#include "curl.cpp"
#include "E0Inside.cpp"
#include "H0Inside.cpp"
#include "E0Scatt.cpp"
#include "H0Scatt.cpp"
#include "D1Inside.cpp"
#include "D1Scatt.cpp"
#include "B1Inside.cpp"
#include "B1Scatt.cpp"
#include "ZeroOrderField.cpp"
#include "FirstOrderField.cpp"

// g++ -o ZeroOrderFields.o ZeroOrderFields.cpp -lgsl -lgslcblas -lm
// ./ZeroOrderFields.o

int main()
{

	double mu_r, epsilon_r, Rsphere, frequency, beta_a, theta_i, R_m, theta_m, phi_m, Npoints ;
	string in_filename_points, out_filename_Efield, out_filename_Hfield;
	cin >> mu_r >> epsilon_r >> Rsphere >> frequency >> beta_a >> theta_i >> in_filename_points >> out_filename_Efield >> out_filename_Hfield;

	double N = sqrt(mu_r*epsilon_r);//refractive index

	double kprop = 2.0*pi*frequency/C0;//propagation constant
		
	double K01 = 0.0;
	dcomplex C1 = 0.0;
	if(beta_a != 0.0)
	{
		K01 = (N*N -1.0)/(beta_a*C0*C0);
		double Omega = beta_a*C0/Rsphere;
		C1 = j*Omega*K01/(N*kprop*sqrt(mu_r*mu_0/(epsilon_r*epsilon_0)));
	}
	
	int nMax = 40;//number of terms used in the series expansion

	//Read the file containing the list of points
	double xm, ym, zm;
	ifstream file_points;
	file_points.open(in_filename_points.c_str());
	file_points >> Npoints;

	//This loop calculate and prints the fields for each point
	for(int i=0; i<Npoints; i++ )
	{
		// Read the co-ordinates of point and convert from rectangular to spherical
		file_points >> xm >> ym >> zm;		
		R_m = sqrt(xm*xm+ym*ym+zm*zm);
		theta_m = atan2(sqrt(xm*xm+ym*ym),zm);
		phi_m = atan2(ym,xm);
		//convert the range of phi_m from 0 to 2*pi
		if(phi_m < 0)
		{
			phi_m += 2*pi;
		}
		//skip calculations for points closer to z-axis than 0.001
		if(sqrt(xm*xm+ym*ym) < 0.001) 
		{
			continue;
		}
		//call the programs for calculating field
		if(beta_a != 0.0)
		{
			TotalField(nMax,N,mu_r,epsilon_r,kprop,Rsphere,theta_i,R_m,theta_m,phi_m, C1, beta_a, K01, out_filename_Efield, out_filename_Hfield);
		}
		else
		{
			ZeroOrderField(nMax,N,mu_r,epsilon_r,kprop,Rsphere,theta_i,R_m,theta_m,phi_m,out_filename_Efield, out_filename_Hfield);
		}
	}
	return 0;
}
