/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by CosmoLike developers
******************************************************************************/

/*relations for converting mass -> radius, using \Delta = 200 \rho_m for consistency with mass & bias function */
double m_Delta(double r, double a);
double r_Delta(double m, double a);

/* fitting functions for mass and bias */
double massfunc(double m, double a);
double B1 (double m, double a);
/* halo density profiles */
double conc(double m, double a);//mass-concentration relation
double u_nfw_c(double c,double k, double m, double a); //analytic expression for Fourier transform of NFW density profile
double u_nfw(double c,double k, double m, double a);   // Fourier transform of NFW density profile, truncated at r_Delta, through direct integration
// use this routine in inner I[0-1]j and modify int_rho_nfw if to work with different halo profiles, e.g. to include adiabatic contraction, AGN feedback, etc.
/*halo model building blocks */
double I0j (int j, double k1, double k2, double k3, double k4, double a);
double I_11 (double k,double a);
double I1j (int j, double k1, double k2, double k3, double a);
/*halo model matter power spectrum, bispectrum, trispectrum*/
double p_1h(double k, double a);
double p_2h(double k, double a);
double Pdelta_halo(double k, double a);
/*look up table for 1-h halo sample variance term */
double I12_SSC (double k,double a);

/*==============================================================*/
/************** begin halo properties *****************/

/********* routines to convert scales/radii to halo masses, and vice versa *********/

double delta_c(double a) /*set to 1.686 for consistency with Tinker mass & bias definition */
{
  return 1.686;
}
double delta_Delta(double a)
{
  return 200.0;//using Tinker et al mass & bias functions with \Delta = 200 rho_{mean}
}

double rho_Delta(double a) //virial density in solar masses/h/(H0/c)^3
{
  //using Tinker et al mass & bias functions with \Delta = 200 rho_{mean}
  return (delta_Delta(a) *cosmology.rho_crit*cosmology.Omega_m);
}

double m_Delta(double r, double a){
  return 4.*M_PI/3.0*pow(r,3.0)*rho_Delta(a);
}
double r_Delta(double m, double a) //calculate r_Delta in c/H0 given m in (solar masses/h)
{
  return pow(3./(4.*M_PI)*(m/rho_Delta(a)),1./3.);
}

double r_s(double m, double a)
{
  return r_Delta(m,a)/conc(m,a);
}

double radius(double m)
{
  return pow(3./4.*m/(M_PI*cosmology.rho_crit*cosmology.Omega_m),1./3.);
}

/*++++++++++++++++++++++++++++++++++++++++*
*  Variance of the density field         *
*++++++++++++++++++++++++++++++++++++++++*/

double sigma2_integrand(double x, void * params)   // inner integral
{
  double *array = (double*)params;
  double k= x/array[0];
  //refactored FT of spherical top-hat to avoid numerica divergence of 1/x
  return p_lin(k,1.0)*pow(3.*gsl_sf_bessel_j1(x)/array[0],2.)/(array[0]*2.*M_PI*M_PI);
}
double sigma2(double m)
{
  static cosmopara C;

  static double *table_S2;
  static double dm = .0, logmmin =1., logmmax = 1.;
  double mlog,result,array[1],x1;
  int i,j;
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);
    if (table_S2==0){
      table_S2  = create_double_vector(0, Ntable.N_S2-1);
      logmmin = log(limits.M_min/2.0);
      logmmax = log(limits.M_max*2.0);
      dm = (logmmax - logmmin)/(Ntable.N_S2);
    }
    mlog = logmmin;

    for (i=0; i<Ntable.N_S2; i++, mlog += dm) {
      array[0] = radius(exp(mlog));
      result  = int_gsl_integrate_medium_precision(sigma2_integrand,(void*)array,0.,14.1,NULL,1000);
      table_S2[i]=log(result);
    }
  }
  return exp(interpol(table_S2, Ntable.N_S2, logmmin, logmmax, dm,log(m), 1.0,1.0 ));
}

double nu(double m, double a){
  static int init = 1;
  if ((cosmology.Omega_nu > 0 || cosmology.M_nu >0) && init){
    fprintf(stderr,"halo.c does not support cosmologies with massive neutrinos\n");
//    exit(1);
    init =0;
  }
	return delta_c(a)/(sqrt(sigma2(m))*growfac(a)/growfac(1.));
}

double nulogm_a1(double lgm,void * params)
{
  (void)(params);
  return log(nu(exp(lgm),1.));
}


double dlognudlogm(double m)
{
  static cosmopara C;

  static double *table_DS;
  static double dm = .0, logmmin = 1.0,logmmax = 1.0;
  double mlog;
  int i;
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);

    if (table_DS==0){
      table_DS  = create_double_vector(0, Ntable.N_DS-1);
      logmmin = log(limits.M_min/2.0);
      logmmax = log(limits.M_max*2.0);
      dm = (logmmax - logmmin)/(Ntable.N_DS);
    }
    gsl_function F;
    double result, abserr;

    F.function = &nulogm_a1;
    F.params = 0;
    mlog = logmmin;
    for (i=0; i<Ntable.N_DS; i++, mlog += dm) {
      gsl_deriv_central (&F,mlog, 0.1*mlog, &result, &abserr);
      table_DS[i]=result;
    }
  }
  return interpol(table_DS, Ntable.N_DS, logmmin, logmmax, dm,log(m), 1.0,1.0);
}

/*************** mass function & halo bias (Tinker et al.) **********/
double f_tinker(double n, double a_in)
{  //Eqs. (8-12) + Table 4 from Tinker et al. 2010
  //aa = alpha, b = beta, c = gamma, f = phi, e = eta
  double a = fmax(0.25,a_in); //limit fit range of mass function evolution to z <= 3, as discussed after Eq. 12 of http://arxiv.org/pdf/1001.3162v2.pdf
	double aa,b,c,f,e;
  aa = 0.368;
  b = 0.589*pow(a,-0.2); c = 0.864*pow(a,0.01); f = -0.729*pow(a,.08); e = -0.243*pow(a,-0.27);
	return aa*(1.0+pow(b*n,-2.0*f))*pow(n,2.0*e)*exp(-c*n*n/2.0);
}
double fnu_tinker(double n, double a)
{
  return f_tinker(n,a)*n;
}

double B1_nu (double n,double a){
  // Eq (6) + Table 2 from Tinker et al. 2010
	double A,aa,B,b,C,c,y;
	y= log10(200.0);
	A = 1.0+0.24*y*exp(-pow(4.0/y,4.0));
	aa = 0.44*y-0.88; B = 0.183; b = 1.5;
	C = 0.019+0.107*y+0.19*exp(-pow(4.0/y,4.0)); c =2.4;
	return 1.0-A*pow(n,aa)/(pow(n,aa)+pow(delta_c(a),aa)) + B*pow(n,b)+C*pow(n,c);
}



//correct bias for halo mass cuts (so that large-scale 2-h term matches PT results at all redshifts)

double bias_norm_integrand (double n,void * params){
  double *array = (double*)params;
  double a = array[0];
  return B1_nu(n,a)*f_tinker(n,a);
}

double bias_norm(double a)
{
	static cosmopara C;
	static double *table_BN;
	static double da = .0, amin = 0., amax =0.;
	double aa,result,array[1];
	int i;
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);
  	if (table_BN==0){
      table_BN  = create_double_vector(0, Ntable.N_a-1);
      amin = limits.a_min; //1./(1+redshift.shear_zdistrpar_zmax);
      amax = 1.;
      da = (amax - amin)/(Ntable.N_a-1.);
    }
		aa= amin;
		for (i=0; i<Ntable.N_a-1; i++, aa += da) {
			array[0] = aa;
			result = int_gsl_integrate_medium_precision(bias_norm_integrand,(void*)array,nu(limits.M_min,aa),nu(limits.M_max,aa),NULL,5000);
			table_BN[i]=result;
		}
    table_BN[Ntable.N_a-1] = 1.;
	}
	return interpol(table_BN, Ntable.N_a, amin,amax, da,fmin(a,amax-da), 1.0,1.0 );
}

double massfunc(double m, double a){
	return fnu_tinker(nu(m,a),a)*cosmology.rho_crit*cosmology.Omega_m/m/m*dlognudlogm(m);
}


double B1 (double m,double a){ //b(m,a) based on redshift evolution fits in Tinker et al. paper, no additional normalization
	return B1_nu(nu(m,a),a);
}
double B1_normalized (double m,double a){ //divide by bias norm only in matter spectra, not in HOD modeling/cluster analyses
	return B1_nu(nu(m,a),a)/bias_norm(a);
}
/***************** halo profiles ***************/
/******** mass-concentration relation **********/
double conc(double m, double a)
{
	return 9.*pow(nu(m,a),-.29)*pow(growfac(a)/growfac(1.),1.15);// Bhattacharya et al. 2013, Delta = 200 rho_{mean} (Table 2)
  //return 10.14*pow(m/2.e+12,-0.081)*pow(a,1.01); //Duffy et al. 2008 (Delta = 200 mean)
}

/***********  FT of NFW Profile **************/

double int_rho_nfw(double r, void *params){//Fourier kernel * NFW profile, integrand for u_nfw
  double *array = (double*)params;
  double k = array[0];
  double c = array[3];
  double rv = array[4];
  double rs =rv/c;
  return r*sinl(k*r)/k*pow(rs,-3.)*1./(log(1.+c)-c/(1.+c))/(r/rs)*pow(1+r/rs,-2.);
}

double u_nfw(double c, double k, double m, double a){
  // FFT density of NFW profile, truncated at r_Delta, through direct integration
  // use this routine in inner I[0-1]j and modify int_rho_nfw to work with different halo profiles

  double array[5] ={k,m,a,c,r_Delta(m,a)};
  return int_gsl_integrate_medium_precision(int_rho_nfw, (void*)array, 0,array[4],NULL, 5000);
}


double u_nfw_c(double c,double k, double m, double aa){// analytic FT of NFW profile, from Cooray & Sheth 01
  double x, xu;
  x = k * r_Delta(m,aa)/c;
  xu = (1.+c)*x;
  return (sin(x)*(gsl_sf_Si(xu)-gsl_sf_Si(x))- sinl(c*x)/xu +cos(x)*(gsl_sf_Ci(xu)-gsl_sf_Ci(x)))*1./(log(1.+c)-c/(1.+c));
}


double inner_I0j (double logm, void *para){
  double *array = (double *) para;
  double m = exp(logm);
  long double u = 1.0;
  double a= array[6];
  double c = conc(m,a);
  int l;
  int j = (int)(array[5]);
  for (l = 0; l< j; l++){
    u = u*u_nfw_c(c,array[l],m,a);
  }
  return massfunc(m,a)*m*pow(m/(cosmology.rho_crit*cosmology.Omega_m),(double)j)*u;
}

double I0j (int j, double k1, double k2, double k3, double k4,double a){
  double array[7] = {k1,k2,k3,k4,0.,(double)j,a};
  return int_gsl_integrate_medium_precision(inner_I0j,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000);
}

double inner_I1j (double logm, void *para){
  double *array = (double *) para;
  double m = exp(logm);
  long double u = 1.0;
  double a= array[6];
  double c = conc(m,a);
  int l;
  int j = (int)(array[5]);
  for (l = 0; l< j; l++){
    u = u*u_nfw_c(c,array[l],m,a);
  }
  return massfunc(m,a)*m*pow(m/(cosmology.rho_crit*cosmology.Omega_m),(double)j)*u*B1_normalized(m,a);
}

double I_11 (double k,double a){//look-up table for I11 integral
  static cosmopara C;
  static double amin = 0, amax = 0,logkmin = 0., logkmax = 0., dk = 0., da = 0.;

  static double **table_I1=0;

  double aa,klog;
  int i,j;

  if (recompute_cosmo3D(C)){ //extend this by halo model parameters if these become part of the model
    update_cosmopara(&C);
    if (table_I1==0) table_I1 = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);
    double array[7];
    array[5]=1.0;

    amin = limits.a_min;
    amax = 1.;
    da = (amax - amin)/(Ntable.N_a);
    aa = amin;
    logkmin = log(limits.k_min_cH0);
    logkmax = log(limits.k_max_cH0);
    dk = (logkmax - logkmin)/(Ntable.N_k_nlin);
    for (i=0; i<Ntable.N_a; i++, aa +=da) {
      array[6] = fmin(aa,0.999);
      klog  = logkmin;
      for (j=0; j<Ntable.N_k_nlin; j++, klog += dk) {
        array[0]= exp(klog);
        table_I1[i][j] = log(int_gsl_integrate_medium_precision(inner_I1j,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000));
      }
    }
  }
  aa = fmin(a,amax-1.1*da);//to avoid interpolation errors near z=0
  return exp(interpol2d(table_I1, Ntable.N_a, amin, amax, da, aa, Ntable.N_k_nlin, logkmin, logkmax, dk, log(k), 1.0, 1.0));
}

double I1j (int j, double k1, double k2, double k3,double a){
  if (j ==1) {return I_11(k1,a);}
  double array[7] = {k1,k2,k3,0.,0.,(double)j,a};
  return int_gsl_integrate_medium_precision(inner_I1j,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000);
}

/*++++++++++++++++++++++++++++++++++++++++*
* 1Halo Terms                             *
*++++++++++++++++++++++++++++++++++++++++*/

double p_1h(double k, double a)
{
  return I0j(2,k,k,0.,0.,a);
}

/*++++++++++++++++++++++++++++++++++++++++*
* 2Halo Terms                         	   *
*++++++++++++++++++++++++++++++++++++++++*/
double p_2h(double k, double a)
{
  return pow(I_11(k,a),2.0)*p_lin(k,a);
}

/*********** Look-up table for halomodel matter power spectrum ************/
double Pdelta_halo(double k,double a)
{
  static cosmopara C;
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;

  static double **table_P_NL=0;

  double klog,val,aa,kk;
  int i,j;

  if (recompute_cosmo3D(C)){ //extend this by halo model parameters if these become part of the model
    update_cosmopara(&C);
    if (table_P_NL!=0) free_double_matrix(table_P_NL,0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);
    table_P_NL = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);

    da = (1. - limits.a_min)/(Ntable.N_a*1.0);
    logkmin = log(limits.k_min_cH0);
    logkmax = log(limits.k_max_cH0);
    dk = (logkmax - logkmin)/(Ntable.N_k_nlin*1.0);
    aa= limits.a_min;
    for (i=0; i<Ntable.N_a; i++, aa +=da) {
      if(aa>0.999) aa=.999;
      klog  = logkmin;
      for (j=0; j<Ntable.N_k_nlin; j++, klog += dk) {
        kk = exp(klog);
        table_P_NL[i][j] = log(p_1h(kk,aa) + p_2h(kk,aa));
      }
    }
  }
  klog = log(k);
  val = interpol2d(table_P_NL, Ntable.N_a, limits.a_min, 1., da, a, Ntable.N_k_nlin, logkmin, logkmax, dk, klog, 0.0, 0.0);
  return exp(val);
}

/****************** Lookup table for 1-halo term super-sample covariance/halo sample variance term**********/
double I12_SSC (double k,double a){//one-halo term contribution to super sample covariance
  static cosmopara C;
  static double amin = 0, amax = 0,logkmin = 0., logkmax = 0., dk = 0., da = 0.;

  static double **table_I1=0;

  double aa,klog;
  int i,j;

  if (recompute_cosmo3D(C)){ //extend this by halo model parameters if these become part of the model
    update_cosmopara(&C);
    if (table_I1==0) table_I1 = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);
    double array[7];
    array[5]=2.0;

    amin = limits.a_min;
    amax = 1.-1.e-5;
    da = (amax - amin)/(Ntable.N_a);
    aa = amin;
    logkmin = log(limits.k_min_cH0);
    logkmax = log(limits.k_max_cH0);
    dk = (logkmax - logkmin)/(Ntable.N_k_nlin);
    for (i=0; i<Ntable.N_a; i++, aa +=da) {
      array[6] = fmin(aa,0.999);
      klog  = logkmin;
      for (j=0; j<Ntable.N_k_nlin; j++, klog += dk) {
        array[0]= exp(klog);array[1]= exp(klog);
        table_I1[i][j] = log(int_gsl_integrate_medium_precision(inner_I1j,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000));
      }
    }
  }
  if (log(k) <logkmin){return 1.0;};
  if (log(k) >logkmax){return 0.0;};
  aa = fmin(a,amax-1.1*da);
  return exp(interpol2d(table_I1, Ntable.N_a, amin, amax, da, aa, Ntable.N_k_nlin, logkmin, logkmax, dk, log(k), 1.0, 1.0));
}
