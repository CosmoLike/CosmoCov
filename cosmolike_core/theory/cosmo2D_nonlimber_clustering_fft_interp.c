/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by CosmoLike developers
******************************************************************************/

void f_chi_for_Psi_cl(double* chi_ar, int Nchi, double* f_chi_ar, int ni);
void f_chi_for_Psi_cl_RSD(double* chi_ar, int Nchi, double* f_chi_RSD_ar, int ni);
void f_chi_for_Psi_cl_Mag(double* chi_ar, int Nchi, double* f_chi_Mag_ar, int ni);
void C_cl_mixed(int L, int LMAX, int ni, int nj, double *Cl, double dev, double tolerance);

double f_growth(double z){
	double aa = 1./(1+z);
	double gamma = 0.55;
	return pow(cosmology.Omega_m /(cosmology.Omega_m +omv_vareos(aa) *aa*aa*aa),gamma);
}


double int_for_C_cl_lin(double a, void *params)
{
	double res,ell, fK, k;
	double *ar = (double *) params;
	ell       = ar[2]+0.5;
	fK     = f_K(chi(a));
	k      = ell/fK;
	
	res=W_gal(a,ar[0])*W_gal(a,ar[1])*dchi_da(a)/fK/fK;
	res= res*p_lin(k,a);
	return res;
}


double C_cl_lin_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{
	double array[3] = {1.0*ni,1.0*nj,l};
	// return int_gsl_integrate_medium_precision(int_for_C_cl_lin,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),fmin(amax_lens(ni),amax_lens(nj)),NULL,1000);
	return int_gsl_integrate_medium_precision(int_for_C_cl_lin,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),0.99999,NULL,1000);
}

/////// Integrand for galaxy density
void f_chi_for_Psi_cl(double* chi_ar, int Nchi, double* f_chi_ar, int ni){
	double g0 =1./growfac(1.);
	double a, z;
	int i;
	double real_coverH0 = cosmology.coverH0 / cosmology.h0; // unit Mpc
	double pf;
	for(i=0;i<Nchi;i++) {
		a = a_chi(chi_ar[i] / real_coverH0) ; // first convert unit of chi from Mpc to c/H0
		z = 1./a - 1.;
		if( (z<tomo.clustering_zmin[ni]) || (z>tomo.clustering_zmax[ni]) )
		{
			f_chi_ar[i] = 0.;
		}
		else
		{
			pf = (pf_photoz(z,ni)<0.)? 0:pf_photoz(z,ni); // get rid of unphysical negatives
			f_chi_ar[i] = chi_ar[i] * pf*growfac(a)*g0*gbias.b1_function(z,ni)*hoverh0(a)/real_coverH0;
		}
	}
}

// Integrand for galaxy density RSD
void f_chi_for_Psi_cl_RSD(double* chi_ar, int Nchi, double* f_chi_RSD_ar, int ni){
	double g0 =1./growfac(1.);
	double a, z;
	int i;
	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double pf;
	for(i=0;i<Nchi;i++) {
		a = a_chi(chi_ar[i] / real_coverH0) ; // first convert unit of chi from Mpc to c/H0
		z = 1./a - 1.;
		if( (z<tomo.clustering_zmin[ni]) || (z>tomo.clustering_zmax[ni]) )
		{
			f_chi_RSD_ar[i] = 0.;
		}
		else
		{
			pf = (pf_photoz(z,ni)<0.)? 0:pf_photoz(z,ni); // get rid of unphysical negatives
			f_chi_RSD_ar[i] = -chi_ar[i] * pf*growfac(a)*g0*f_growth(z)*hoverh0(a)/real_coverH0;
		}
	}
}

// Integrand for lensing magnification of galaxy density
void f_chi_for_Psi_cl_Mag(double* chi_ar, int Nchi, double* f_chi_Mag_ar, int ni){
	double g0 =1./growfac(1.);
	double a, z, fK;
	int i;
	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double window_M, wmag;
	for(i=0;i<Nchi;i++) {
		a = a_chi(chi_ar[i] / real_coverH0) ; // first convert unit of chi from Mpc to c/H0
		z = 1./a - 1.;
		fK = f_K(chi_ar[i]/real_coverH0);
		if( (z>tomo.clustering_zmax[ni]) )
		{
			f_chi_Mag_ar[i] = 0.;
		}
		else
		{
			// printf("Here! a, fK, ni: %lg,%lg,%d\n", a, fK, ni);
			wmag = W_mag(a, fK, (double)ni);
			window_M = wmag/ fK / (real_coverH0*real_coverH0);
			// printf("bmag, wkappa, f_K, real_coverH0, %lg %lg %lg %lg\n", gbias.b_mag[ni], wkappa, fK,real_coverH0);
			// pf = (pf_photoz(z,ni)<0.)? 0:pf_photoz(z,ni); // get rid of unphysical negatives
			// f_chi_Mag_ar[i] = chi_ar[i]/a * window_M*growfac(a)*g0;
			f_chi_Mag_ar[i] = window_M*growfac(a)*g0; // unit [Mpc^-2]
		}
		// printf("%lg\n", f_chi_Mag_ar[i]);
	}
}

// Mixture of non-Limber and Limber of C_cl (galaxy clustering)
void C_cl_interp(int LMAX, int ni, int nj, double *Cl) {

	int i,j,i_block;
	long l;
	// run 100 ells at a time, and see if switching to Limber is needed.
	// Save runtime for Limber, and save re-creation time of fftw_plan.
	int Nell_block = 100, Nchi = 1000;
	int ell_ar[Nell_block];
	double **k1_ar, **k2_ar, **Fk1_ar, **Fk2_ar;
	double **Fk1_Mag_ar, **Fk2_Mag_ar;

	k1_ar = malloc(Nell_block * sizeof(double *));
	k2_ar = malloc(Nell_block * sizeof(double *));
	Fk1_ar = malloc(Nell_block * sizeof(double *));
	Fk2_ar = malloc(Nell_block * sizeof(double *));

	Fk1_Mag_ar = malloc(Nell_block * sizeof(double *));
	Fk2_Mag_ar = malloc(Nell_block * sizeof(double *));
	for(i=0;i<Nell_block;i++) {
		k1_ar[i] = malloc(Nchi * sizeof(double));
		k2_ar[i] = malloc(Nchi * sizeof(double));
		Fk1_ar[i] = malloc(Nchi * sizeof(double));
		Fk2_ar[i] = malloc(Nchi * sizeof(double));
		Fk1_Mag_ar[i] = malloc(Nchi * sizeof(double));
		Fk2_Mag_ar[i] = malloc(Nchi * sizeof(double));
		for(j=0;j<Nchi;j++) {
			Fk1_ar[i][j] = 0.;
			Fk2_ar[i][j] = 0.;
			Fk1_Mag_ar[i][j] = 0.;
			Fk2_Mag_ar[i][j] = 0.;
		}
	}

	double chi_ar[Nchi], f1_chi_ar[Nchi], f2_chi_ar[Nchi];
	double f1_chi_RSD_ar[Nchi], f2_chi_RSD_ar[Nchi];
	double f1_chi_Mag_ar[Nchi], f2_chi_Mag_ar[Nchi];

	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double chi_min = chi(1./(1.+0.002))*real_coverH0, chi_max = chi(1./(1.+4.))*real_coverH0;
	double dlnchi = log(chi_max/chi_min) / (Nchi - 1.);
	double dlnk = dlnchi;

	for(i=0; i<Nchi; i++) {
		chi_ar[i] = chi_min * exp(dlnchi*i);
	}

	f_chi_for_Psi_cl(chi_ar, Nchi, f1_chi_ar, ni);
	if(ni != nj) {f_chi_for_Psi_cl(chi_ar, Nchi, f2_chi_ar, nj);}

	f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f1_chi_RSD_ar, ni);
	if(ni != nj) {f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f2_chi_RSD_ar, nj);}

	f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f1_chi_Mag_ar, ni);
	if(ni != nj) {f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f2_chi_Mag_ar, nj);}

	i_block = 0;
	double cl_temp;

	config my_config, my_config_RSD, my_config_Mag;
	my_config.nu = 1.;
	my_config.c_window_width = 0.25;
	my_config.derivative = 0;
	my_config.N_pad = 200;
	my_config.N_extrap_low = 0;
	my_config.N_extrap_high = 0;

	my_config_RSD.nu = 1.01;
	my_config_RSD.c_window_width = 0.25;
	my_config_RSD.derivative = 2;
	my_config_RSD.N_pad = 500;
	my_config_RSD.N_extrap_low = 0;
	my_config_RSD.N_extrap_high = 0;

	my_config_Mag.nu = 1.;
	my_config_Mag.c_window_width = 0.25;
	my_config_Mag.derivative = 0;
	my_config_Mag.N_pad = 500;
	my_config_Mag.N_extrap_low = 0;
	my_config_Mag.N_extrap_high = 0;

	double ell_prefactor;

	double k1_cH0;

	// compute L:0-99
	for(i=0;i<Nell_block;i++) {ell_ar[i]=i+i_block*Nell_block;}

	cfftlog_ells(chi_ar, f1_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k1_ar, Fk1_ar);
	if(ni != nj) {cfftlog_ells(chi_ar, f2_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k2_ar, Fk2_ar);}

	cfftlog_ells_increment(chi_ar, f1_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, k1_ar, Fk1_ar);
	if(ni != nj) {cfftlog_ells_increment(chi_ar, f2_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, k2_ar, Fk2_ar);}

	// Add in lensing magnification contribution
	cfftlog_ells(chi_ar, f1_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar, Nell_block, k1_ar, Fk1_Mag_ar);
	if(ni != nj) {cfftlog_ells(chi_ar, f2_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar, Nell_block, k2_ar, Fk2_Mag_ar);}
	for(i=0;i<Nell_block;i++) {
		ell_prefactor = ell_ar[i]*(ell_ar[i]+1.);
		for(j=0;j<Nchi;j++) {
			Fk1_ar[i][j]+= (ell_prefactor / (k1_ar[i][j]*k1_ar[i][j]) * (gbias.b_mag[ni]) *  Fk1_Mag_ar[i][j]);
			if(ni != nj) {Fk2_ar[i][j]+= (ell_prefactor / (k2_ar[i][j]*k2_ar[i][j])* (gbias.b_mag[nj]) *  Fk2_Mag_ar[i][j]);}
		}
	}

	for(i=0;i<Nell_block;i++) {
		cl_temp = 0.;
		for(j=0;j<Nchi;j++) {
			// printf("k,Fk: %d,%d, %lf,%lf\n", i,j, k1_ar[i][j], Fk1_ar[i][j]);
			k1_cH0 = k1_ar[i][j] * real_coverH0;
			if(ni == nj) {
				cl_temp += (Fk1_ar[i][j]) * (Fk1_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0,1.0);
				// cl_temp += (Fk1_ar[i][j]) * (Fk1_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *Pdelta(k1_cH0,0.9999);
				// printf("plin,%lg, %lg\n", k1_ar[i][j],p_lin(k1_cH0,1.0));
			}
			else {
				cl_temp += (Fk1_ar[i][j])*(Fk2_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0,1.0);
			}
		}
		Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI + C_cl_tomo_nointerp(1.*ell_ar[i],ni,nj) - C_cl_lin_nointerp(1.*ell_ar[i],ni,nj);
		// printf("ell_ar[i],%d, %lg\n", ell_ar[i],Cl[ell_ar[i]]);
	}

	////// L: 100 - 100000 log-spaced
	for(i=0;i<Nell_block;i++) {ell_ar[i]=round(100.*exp(i*log((LMAX-1)/100.)/(Nell_block-1)));}

	cfftlog_ells(chi_ar, f1_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k1_ar, Fk1_ar);
	if(ni != nj) {cfftlog_ells(chi_ar, f2_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k2_ar, Fk2_ar);}

	cfftlog_ells_increment(chi_ar, f1_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, k1_ar, Fk1_ar);
	if(ni != nj) {cfftlog_ells_increment(chi_ar, f2_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, k2_ar, Fk2_ar);}

	// Add in lensing magnification contribution
	cfftlog_ells(chi_ar, f1_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar, Nell_block, k1_ar, Fk1_Mag_ar);
	if(ni != nj) {cfftlog_ells(chi_ar, f2_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar, Nell_block, k2_ar, Fk2_Mag_ar);}
	for(i=0;i<Nell_block;i++) {
		ell_prefactor = ell_ar[i]*(ell_ar[i]+1.);
		for(j=0;j<Nchi;j++) {
			Fk1_ar[i][j]+= (ell_prefactor / (k1_ar[i][j]*k1_ar[i][j]) * (gbias.b_mag[ni]) *  Fk1_Mag_ar[i][j]);
			if(ni != nj) {Fk2_ar[i][j]+= (ell_prefactor / (k2_ar[i][j]*k2_ar[i][j])* (gbias.b_mag[nj]) *  Fk2_Mag_ar[i][j]);}
		}
	}
	
	for(i=0;i<Nell_block;i++) {
		cl_temp = 0.;
		for(j=0;j<Nchi;j++) {
			k1_cH0 = k1_ar[i][j] * real_coverH0;
			if(ni == nj) {
				cl_temp += (Fk1_ar[i][j]) * (Fk1_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0,1.0);
			}
			else {
				cl_temp += (Fk1_ar[i][j])*(Fk2_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0,1.0);
			}
		}
		Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI + C_cl_tomo_nointerp(1.*ell_ar[i],ni,nj) - C_cl_lin_nointerp(1.*ell_ar[i],ni,nj);
		// printf("ell_ar[i],%d, %lg\n", ell_ar[i],Cl[ell_ar[i]]);
	}

	int ell_interp, ell0, ell1;
	double Cl0, Cl1;
	for(i=0;i<Nell_block-1;i++) {
		ell0 = ell_ar[i]; ell1 = ell_ar[i+1];
		Cl0 = Cl[ell0]; Cl1 = Cl[ell1];
		for(ell_interp=ell0+1; ell_interp<ell1; ell_interp++){
			Cl[ell_interp] = log(ell_interp/(float)ell0) * (Cl1 - Cl0)/log(ell1/(float)ell0) + Cl0;
			// printf("ell_ar[i],%d, %lg, %lg,%lg, %lg\n", ell_interp,log(ell_interp/ell0), Cl0, (Cl1 - Cl0)/log(ell1/ell0), (ell1/ell0));
		}
	}

	// printf("finished bin %d\n", ni);
	for(i=0;i<Nell_block;i++) {
		free(k1_ar[i]);free(k2_ar[i]);
		free(Fk1_ar[i]);free(Fk2_ar[i]);
		free(Fk1_Mag_ar[i]);free(Fk2_Mag_ar[i]);
	}
	free(k1_ar);free(k2_ar);
	free(Fk1_ar);free(Fk2_ar);
	free(Fk1_Mag_ar);free(Fk2_Mag_ar);

}

double C_cl_tomo_nonlimber_interp(double l, int ni, int nj){
	static int LMAX = 100000;
	static double *Cl =0;
	static double *Cl_nz =0;

	int nz;
	if(Cl==0){
		int i;
		int N_cltomo = tomo.clustering_Nbin*(tomo.clustering_Nbin+1)/2;
		Cl = create_double_vector(0,LMAX-1);
		Cl_nz = create_double_vector(0,N_cltomo*LMAX-1);
		for (nz = 0; nz<N_cltomo; nz++){
			C_cl_interp(LMAX, Zcl1(nz),Zcl2(nz), Cl);
			for(i=0;i<LMAX;i++){
				Cl_nz[nz*LMAX+i] = Cl[i];
			}
		}
	}

	nz = N_clustering_tomo(ni,nj);
	if((int)l == l){
		return Cl_nz[nz*LMAX+(int)l];
	}else{
		int l_floor = floor(l);
		return (Cl_nz[nz*LMAX+l_floor+1]-Cl_nz[nz*LMAX+l_floor])/log((l_floor+1)/l_floor) * log(l/l_floor) + Cl_nz[nz*LMAX+l_floor];
	}
}
