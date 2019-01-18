dcomplex anm(double N, double mu_r, double epsilon_r, dcomplex h2nkaHankel, double jnNkaBessel, double jnkaBessel, dcomplex Dh2nkaHankel,  double DjnNkaBessel,  double DjnkaBessel)
{
	dcomplex A11 = h2nkaHankel*sqrt(mu_r/epsilon_r);
	dcomplex A12 = jnNkaBessel;
	dcomplex A21 = Dh2nkaHankel*N;
	dcomplex A22 = DjnNkaBessel;
	dcomplex B1 = jnkaBessel*sqrt(mu_r/epsilon_r);
	dcomplex B2 = DjnkaBessel*N;
	return (B1*A22-B2*A12)/(A11*A22-A21*A12);
}

dcomplex bnm(double N, double mu_r, double epsilon_r, dcomplex h2nkaHankel, double jnNkaBessel, double jnkaBessel, dcomplex Dh2nkaHankel,  double DjnNkaBessel,  double DjnkaBessel)
{
	dcomplex A11 = h2nkaHankel;
	dcomplex A12 = jnNkaBessel;
	dcomplex A21 = mu_r*Dh2nkaHankel;
	dcomplex A22 = DjnNkaBessel;
	dcomplex B1 = jnkaBessel;
	dcomplex B2 = mu_r*DjnkaBessel;
	return (B1*A22-B2*A12)/(A11*A22-A21*A12);
}

dcomplex cnm(double N, double mu_r, double epsilon_r, dcomplex h2nkaHankel, double jnNkaBessel, double jnkaBessel, dcomplex Dh2nkaHankel,  double DjnNkaBessel,  double DjnkaBessel)
{
	dcomplex A11 = h2nkaHankel*sqrt(mu_r/epsilon_r);
	dcomplex A12 = jnNkaBessel;
	dcomplex A21 = Dh2nkaHankel*N;
	dcomplex A22 = DjnNkaBessel;
	dcomplex B1 = jnkaBessel*sqrt(mu_r/epsilon_r);
	dcomplex B2 = DjnkaBessel*N;
	return (A11*B2-A21*B1)/(A11*A22-A21*A12);
}

dcomplex dnm(double N, double mu_r, double epsilon_r, dcomplex h2nkaHankel, double jnNkaBessel, double jnkaBessel, dcomplex Dh2nkaHankel,  double DjnNkaBessel,  double DjnkaBessel)
{
	dcomplex A11 = h2nkaHankel;
	dcomplex A12 = jnNkaBessel;
	dcomplex A21 = mu_r*Dh2nkaHankel;
	dcomplex A22 = DjnNkaBessel;
	dcomplex B1 = jnkaBessel;
	dcomplex B2 = mu_r*DjnkaBessel;
	return (A11*B2-A21*B1)/(A11*A22-A21*A12);
}
