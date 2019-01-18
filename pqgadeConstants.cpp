dcomplex pnm(double N, double epsilon_r, double kprop, dcomplex C1, dcomplex Cnm1, dcomplex h2nkaHankel, double jnNkaBessel,dcomplex dnm1, dcomplex fnm1,  dcomplex Dh2nkaHankel,  double DjnNkaBessel,  dcomplex Dfnm1)
{
	dcomplex A11 = jnNkaBessel;
	dcomplex A12 = -epsilon_r*h2nkaHankel;
	dcomplex A21 = DjnNkaBessel;
	dcomplex A22 = -N*N*Dh2nkaHankel;
	dcomplex B1 = -C1*jnNkaBessel*dnm1 - fnm1;
	dcomplex B2 = -C1*DjnNkaBessel*dnm1-Dfnm1-j*kprop*N*C1*Cnm1;

	return (B1*A22-B2*A12)/(A11*A22-A21*A12);
}

dcomplex qnm(double N, double epsilon_r, dcomplex C1, dcomplex Anm1, dcomplex h2nkaHankel, double jnNkaBessel,dcomplex gnm1,  dcomplex Dh2nkaHankel,  double DjnNkaBessel,  dcomplex Dgnm1)
{
	dcomplex A11 = jnNkaBessel;
	dcomplex A12 = -h2nkaHankel;
	dcomplex A21 = DjnNkaBessel;
	dcomplex A22 = -epsilon_r*Dh2nkaHankel;
	dcomplex B1 = -gnm1;
	dcomplex B2 = -C1*Anm1-Dgnm1;

	return (B1*A22-B2*A12)/(A11*A22-A21*A12);
}

dcomplex gammanm(double N, double epsilon_r, double kprop, dcomplex C1, dcomplex Cnm1, dcomplex h2nkaHankel, double jnNkaBessel,dcomplex dnm1, dcomplex fnm1,  dcomplex Dh2nkaHankel,  double DjnNkaBessel,  dcomplex Dfnm1)
{
	dcomplex A11 = jnNkaBessel;
	dcomplex A12 = -epsilon_r*h2nkaHankel;
	dcomplex A21 = DjnNkaBessel;
	dcomplex A22 = -N*N*Dh2nkaHankel;
	dcomplex B1 = -C1*jnNkaBessel*dnm1 - fnm1;
	dcomplex B2 = -C1*DjnNkaBessel*dnm1-Dfnm1-j*kprop*N*C1*Cnm1;

	return (A11*B2-A21*B1)/(A11*A22-A21*A12);
}

dcomplex deltanm(double N, double epsilon_r, dcomplex C1, dcomplex Anm1, dcomplex h2nkaHankel, double jnNkaBessel,dcomplex gnm1,  dcomplex Dh2nkaHankel,  double DjnNkaBessel,  dcomplex Dgnm1)
{
	dcomplex A11 = jnNkaBessel;
	dcomplex A12 = -h2nkaHankel;
	dcomplex A21 = DjnNkaBessel;
	dcomplex A22 = -epsilon_r*Dh2nkaHankel;
	dcomplex B1 = -gnm1;
	dcomplex B2 = -C1*Anm1-Dgnm1;

	return (A11*B2-A21*B1)/(A11*A22-A21*A12);
}
