/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by CosmoLike developers
******************************************************************************/

double W_kappa(double a, double fK, double nz);//complete lens efficiency weight
double W_source(double a, double nz); //source redshift distribution (radial weight for IA,source clustering)
double W_gal(double a, double nz); //complete weight for galaxy statistics

double C_cl_tomo(double l, int ni, int nj);  //galaxy clustering power spectrum of galaxies in bins ni, nj
double C_cl_tomo_nointerp(double l, int ni, int nj);

double C_gl_tomo(double l, int ni, int nj);  //G-G lensing power spectrum, lens bin ni, source bin nj
double C_gl_tomo_nointerp(double l, int ni, int nj);

double C_shear_tomo(double l, int ni, int nj); //shear tomography power spectra
double C_shear_tomo_nointerp(double l, int ni, int nj);

/*** fast routines for total shear-shear + ggl signal for like.IA = 1****/
double C_shear_shear_IA(double s, int ni, int nj);
double C_ggl_IA(double s, int ni, int nj);
/*** tabulated versions of these for angular correlation functions ****/
double C_shear_shear_IA_tab(double l, int ni, int nj);
double C_ggl_IA_tab(double l, int ni, int nj);
/**********************************************************/

double MG_Sigma(double a)
{
  double omegav=omv_vareos(a);
  double hub=hoverh0(a);
  hub = hub*hub;
 
  return cosmology.MGSigma*omegav/hub/cosmology.Omega_v;
}

double dchi_da(double a){
  return 1./(a*a*hoverh0(a));
}
double W_kappa(double a, double fK, double nz){
  double wkappa = 1.5*cosmology.Omega_m*fK/a*g_tomo(a,(int)nz);
  if(cosmology.MGSigma != 0.){
    wkappa *= (1.+MG_Sigma(a));
  }
  return wkappa;
}

double W_mag(double a, double fK, double nz){
  double wmag = 1.5*cosmology.Omega_m*fK/a*g_lens(a,(int)nz);
  if(cosmology.MGSigma != 0.){
    wmag *= (1.+MG_Sigma(a));
  }
  return wmag;
}

double W_gal(double a, double nz){
  double wgal = gbias.b1_function(1./a-1.,(int)nz)*pf_photoz(1./a-1.,(int)nz)*hoverh0(a);
  double wmag = gbias.b_mag[(int)nz]*1.5*cosmology.Omega_m*f_K(chi(a))/a*g_lens(a,(int)nz);
  if(cosmology.MGSigma != 0.){
    wmag *= (1.+MG_Sigma(a));
  }
  return wgal + wmag;
  // return wgal;
}

double W_source(double a, double nz){
  return zdistr_photoz(1./a-1.,(int)nz)*hoverh0(a);
}

double f_rsd (double aa){
  double gamma = 0.55;
  return pow(cosmology.Omega_m /(cosmology.Omega_m +omv_vareos(aa) *aa*aa*aa),gamma);
}
double W_RSD(double l, double a0, double a1, double nz){
  double w = (1+8.*l)/((2*l+1)*(2*l+1))*pf_photoz(1./a0-1.,(int)nz)*hoverh0(a0)*f_rsd(a0);
  w -= 4./(2*l+3)*sqrt((2*l+1.)/(2*l+3.))*pf_photoz(1./a1-1.,(int)nz)*hoverh0(a1)*f_rsd(a1);
  return w;
}

/*********** Limber approximation integrands for angular power spectra *************/
// use W_i W_j/fK^2 dchi_da(a) as Limber weight for power spectra

double int_for_C_cl_tomo(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_cl_tomo");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  res=W_gal(a,ar[0])*W_gal(a,ar[1])*dchi_da(a)/fK/fK;
  res= res*Pdelta(k,a);
  return res;
}

double int_for_C_cl_tomo_RSD(double k, void *params)
{
  double res,ell, chi_0, chi_1, a_0, a_1;
//  k = exp(lgk);
  double *ar = (double *) params;
  ell       = ar[2]+0.5;
  chi_0 = f_K(ell/k);
  chi_1 = f_K((ell+1.)/k);
  if (chi_1 > chi(limits.a_min)){
    return 0;}
  a_0 = a_chi(chi_0);
  a_1 = a_chi(chi_1);
//  res=(W_gal(a_0,ar[0]))*(W_gal(a_0,ar[1]));
  res=(W_gal(a_0,ar[0])+W_RSD(ell,a_0,a_1,ar[0]))*(W_gal(a_0,ar[1])+W_RSD(ell,a_0,a_1,ar[1]));
  res= res*Pdelta(k,a_0);
  return res;
}

double int_for_C_gl_tomo(double a, void *params) // Add RSD
{
  double *ar = (double *) params;
  double res,ell, fK, k;
  if (a >= 1.0) error("a>=1 in int_for_C_gl_tomo");

  double ell_prefactor1 = (ar[2])*(ar[2]+1.);
  double ell_prefactor2 = (ar[2]-1.)*ell_prefactor1*(ar[2]+2.);
  if(ell_prefactor2<=0.) 
    ell_prefactor2=0.;
  else
    ell_prefactor2=sqrt(ell_prefactor2);

  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  double chi_0,chi_1,a_0,a_1;
  chi_0 = f_K(ell/k);
  chi_1 = f_K((ell+1.)/k);
  if (chi_1 > chi(limits.a_min)){
    return 0;}
  a_0 = a_chi(chi_0);
  a_1 = a_chi(chi_1);

  double wgal = W_gal(a,ar[0]);
  wgal += W_mag(a,fK,ar[0])*(ell_prefactor1/ell/ell -1.) ;
  // res= (wgal+W_RSD(ell, a_0, a_1, ar[0]))*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK  * ell_prefactor2/ell/ell;
  res= (wgal)*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK  * ell_prefactor2/ell/ell;
  // res= (W_gal(a,ar[0]) + W_mag(a,fK,ar[0])*gbias.b_mag[(int)ar[0]] *(ell_prefactor1/ell/ell-1.) )*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK /fK/fK * ell_prefactor2/k/k;
  res= res*Pdelta(k,a);
  return res;
}

double int_for_C_shear_tomo(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell, fK, k;
  if (a >= 1.0) error("a>=1 in int_for_C_shear_tomo");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  res= W_kappa(a,fK,ar[0])*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK;
  res= res*Pdelta(k,a); 
  return res;
}

/*********** angular power spectra - without look-up tables ******************/
double C_cl_RSD_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{ 
  double array[3] = {1.0*ni,1.0*nj,l};
  double chi_max = chi(fmax(amin_lens(ni),amin_lens(nj)));
  double chi_min = chi(fmin(amax_lens(ni),amax_lens(nj)));
  double k_min = (l+0.5)/f_K(chi_max);
  double k_max = (l+0.5)/f_K(chi_min);
  return int_gsl_integrate_medium_precision(int_for_C_cl_tomo_RSD,(void*)array,k_min,k_max,NULL,1000)*2./(2*l+1.);
}

double C_cl_tomo_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{ static int init =-1;
  double array[3] = {1.0*ni,1.0*nj,l};
  if (ni == nj){
    return int_gsl_integrate_medium_precision(int_for_C_cl_tomo,(void*)array,amin_lens(ni),0.999999,NULL,1000);
  }
  // return int_gsl_integrate_medium_precision(int_for_C_cl_tomo,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),fmin(amax_lens(ni),amax_lens(nj)),NULL,1000);
  return int_gsl_integrate_medium_precision(int_for_C_cl_tomo,(void*)array,amin_lens(nj),0.99999,NULL,1000); // zi<=zj
}


double C_gl_tomo_nointerp(double l, int ni, int nj)  //G-G lensing power spectrum, lens bin ni, source bin nj
{
  // l=1.0*2;
  if(l==1.) {return 0.;}
  double array[3] = {(double)ni,(double)nj,l};
  double res = int_gsl_integrate_medium_precision(int_for_C_gl_tomo,(void*)array,amin_lens(ni),0.99999,NULL,1000);
  return res;
}

double C_shear_tomo_nointerp(double l, int ni, int nj) //shear tomography power spectra of source galaxy bins ni, nj
{
  double res;
  double array[3] = {(double)ni, (double)nj,l};
  // printf("ni,nj:%d,%d\n", ni,nj );
  int j,k;
  if (ni <= nj){j =nj; k = ni;}
  else{j = ni; k = nj;}
  res=int_gsl_integrate_medium_precision(int_for_C_shear_tomo,(void*)array,amin_source(j),amax_source(k),NULL,1000);
  // printf("res:%lg\n", res);
  return res;
  
}

/*********** angular power spectra - with look-up tables ******************/

double C_cl_tomo(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxies in bins ni, nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.clustering_Nbin ||nj < 0 || nj >= tomo.clustering_Nbin){
    printf("C_cl_tomo(l,%d,%d) outside tomo.clustering_Nbin range\nEXIT\n",ni,nj); exit(1);
  }

  if (recompute_clustering(C,G,N,ni,nj))
  {
    if (table==0) {
      table   = create_double_matrix(0, tomo.clustering_Nbin*tomo.clustering_Nbin-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    for (int i = 0; i < tomo.clustering_Nbin*tomo.clustering_Nbin; i++){table[i][0] = 123456789.0;}
    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
  }
  int j = ni*tomo.clustering_Nbin+nj;
  if(table[j][0] > 123456780.0){ //still need to recompute this tomography bin combination
    //printf("recompute C_cl_tomo, %d, %d\n", ni, nj) ;
    double llog = logsmin;
    double result;
    for (int i=0; i<Ntable.N_ell; i++, llog+=ds) {

//      table[j][i]= log(C_cl_RSD_nointerp(exp(llog),ni,nj));
      result = C_cl_tomo_nointerp(exp(llog),ni,nj);
      if(result<=0) table[j][i] = -100; 
      else table[j][i] = log(result);
      table[nj*tomo.clustering_Nbin+ni][i]=table[j][i];
    }
  }
  double f1 = exp(interpol_fitslope(table[j], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.));
  if (isnan(f1)){f1 = 0.;}
  return f1;
}

double C_gl_tomo(double l, int ni, int nj)  //G-G lensing power spectrum, lens bin ni, source bin nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1 = 0.;
  
  if (ni < 0 || ni >= tomo.clustering_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_gl_tomo(l,%d,%d) outside tomo.X_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_ggl(C,G,N,ni)){
    if (table==0){
      table   = create_double_matrix(0, tomo.ggl_Npowerspectra-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell-1.);
    }
    int i,k;
    double llog;
    
    for (k=0; k<tomo.ggl_Npowerspectra; k++) {
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
          table[k][i]= log(C_gl_tomo_nointerp(exp(llog),ZL(k),ZS(k)));
      }
    }
    
    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
    
  }
  if(test_zoverlap(ni,nj)){f1 = exp(interpol_fitslope(table[N_ggl(ni,nj)], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.));}
  if (isnan(f1)){f1 = 0;}
  return f1;
}

double C_shear_tomo(double l, int ni, int nj)  //shear power spectrum of source galaxies in bins ni, nj
{
  static cosmopara C;
  static nuisancepara N;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.shear_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_shear_tomo(l,%d,%d) outside tomo.shear_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_shear(C,N)){
    if (table==0) {
      table   = create_double_matrix(0, tomo.shear_Npowerspectra-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell-1.);
    }
    
    double llog;
    int i,k;
    
    for (k=0; k<tomo.shear_Npowerspectra; k++) {
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
        table[k][i]= log(C_shear_tomo_nointerp(exp(llog),Z1(k),Z2(k)));
       // printf("table[k][i]: %lg:\n",table[k][i]);
      }
    }
    update_cosmopara(&C); update_nuisance(&N); 
  }
  double f1 = exp(interpol_fitslope(table[N_shear(ni,nj)], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.));
  if (isnan(f1)){f1 = 0.;}
  // printf("%le %d %d %le\n", l, ni, nj, f1);
  return f1;
}


/******** faster shear/ggl + IA routines *********/

double int_for_C_ggl_IA(double a, void *params)
{
  double res, ell, fK, k,norm;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");

  double ell_prefactor1 = (ar[2])*(ar[2]+1.);
  double ell_prefactor2 = (ar[2]-1.)*ell_prefactor1*(ar[2]+2.);
  if(ell_prefactor2<=0.) 
    ell_prefactor2=0.;
  else
    ell_prefactor2=sqrt(ell_prefactor2);

  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;

  double chi_0,chi_1,a_0,a_1;
  chi_0 = f_K(ell/k);
  chi_1 = f_K((ell+1.)/k);
  if (chi_1 > chi(limits.a_min)){
    return 0;}
  a_0 = a_chi(chi_0);
  a_1 = a_chi(chi_1);

  norm = cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(1.)/growfac(a)*nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
  res= (W_gal(a,ar[0]) +W_RSD(ell, a_0, a_1, ar[0]) +W_mag(a,fK,ar[0])*(ell_prefactor1/ell/ell -1.))*(W_kappa(a,fK,ar[1])-W_source(a,ar[1])*norm);
  return res*Pdelta(k,a)*dchi_da(a)/fK/fK * ell_prefactor2/ell/ell;
}

double int_for_C_shear_shear_IA(double a, void *params)
{
  double res, ell, fK, k,ws1,ws2,wk1,wk2, norm;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  ws1 = W_source(a,ar[0]);
  ws2 = W_source(a,ar[1]);
  wk1 = W_kappa(a,fK,ar[0]);
  wk2 = W_kappa(a,fK,ar[1]);
  norm = cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(1.)/growfac(a)*nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
  res= ws1*ws2*norm*norm - (ws1*wk2+ws2*wk1)*norm+wk1*wk2;
  
  return res*Pdelta(k,a)*dchi_da(a)/fK/fK;
}

double C_shear_shear_IA(double s, int ni, int nj)
{
 double array[3] = {(double) ni, (double) nj,s};
  int j=nj,k=ni;
  if (ni > nj){error("cosmo2D_fourier.c: C_shear_shear_IA: ni>nj not supported!");}
  return int_gsl_integrate_medium_precision(int_for_C_shear_shear_IA,(void*)array,amin_source(j),amax_source(k),NULL,1000);
}

double C_ggl_IA(double s, int nl, int ns)
{
  double array[3] = {(double) nl, (double) ns,s};
  return int_gsl_integrate_medium_precision(int_for_C_ggl_IA,(void*)array,amin_lens(nl),0.999999,NULL,1000);
}
/********** tabulated versions of C_ggl_IA and C_shear_shear_IA******/
/************for calculation of angular corelation functions ********/
double C_ggl_IA_tab(double l, int ni, int nj)  //G-G lensing power spectrum, lens bin ni, source bin nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  
  if (ni < 0 || ni >= tomo.clustering_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_gl_tomo(l,%d,%d) outside tomo.X_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_ggl(C,G,N,ni)){
      //printf("calculating C_ggl_IA_tab  %e %e %e %e %e\n", nuisance.A_z[0], nuisance.A_z[1], nuisance.A_z[2], nuisance.A_z[3], nuisance.A_z[4]);
    if (table==0){
      table   = create_double_matrix(0, tomo.ggl_Npowerspectra-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    int i,k;
    double llog;
    
    for (k=0; k<tomo.ggl_Npowerspectra; k++) {
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
        table[k][i] = C_ggl_IA(exp(llog),ZL(k),ZS(k));
      }
    }

    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
    
  }
  int k = N_ggl(ni,nj);
  double f1 = 0.;
  if(test_zoverlap(ni,nj)){f1 = interpol_fitslope(table[k], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.);}
  if (isnan(f1)){f1 = 0;}
  return f1;
}
double C_shear_shear_IA_tab(double l, int ni, int nj)  //shear power spectrum of source galaxies in bins ni, nj
{
  static cosmopara C;
  static nuisancepara N;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.shear_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_shear_tomo(l,%d,%d) outside tomo.shear_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_shear(C,N)){
    //printf("calculating C_shear_shear_IA_tab  %e %e %e %e %e\n", nuisance.A_z[0], nuisance.A_z[1], nuisance.A_z[2], nuisance.A_z[3], nuisance.A_z[4]);
    if (table==0) {
      table   = create_double_matrix(0, tomo.shear_Npowerspectra-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    
    double llog;
    int i,k;
    
    for (k=0; k<tomo.shear_Npowerspectra; k++) {
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
        table[k][i]= (C_shear_shear_IA(exp(llog),Z1(k),Z2(k)));
      }
    }
    update_cosmopara(&C); update_nuisance(&N); 
  }
  double f1 = (interpol_fitslope(table[N_shear(ni,nj)], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.));
  if (isnan(f1)){f1 = 0.;}
  // printf("%le %d %d %le\n", l, ni, nj, f1);
  return f1;
}
