/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by CosmoLike developers
******************************************************************************/

/*bispectrum & trispectrum routines for arbitrary configurations*/
double fs_2(double k1x, double k1y, double k2x, double k2y);
double b_lin (double k1x, double k1y, double k2x, double k2y, double a);
double t_lin(double k1x, double k1y, double k2x, double k2y, double k3x, double k3y, double a);
double b_lin_cov (double k1,double k2, double a);
double b_lin_cov_mu(double k1, double k2, double a);
double bi_1h(double k1,double k2, double k3, double a);
double tri_1h(double k1, double k2, double k3, double k4,double a);
double tri_2h (double k1x,double k1y,double k2x,double k2y,double k3x,double k3y,double a);

/*matter trispectrum for averaged all parallelogram configurations contributing to covariance T(k1,k2)*/
double tri_1h_cov(double k1, double k2, double a);
double tri_2h_cov(double k1, double k2, double a);
double tri_3h_cov(double k1, double k2, double a);
double tri_4h_cov(double k1, double k2, double a);
double tri_matter_cov(double k1, double k2, double a); //4h+3h+2h+1h trispectrum terms in covariance configuration
double tri_multih_cov(double k1, double k2, double a); //4h+3h+2h terms in covariance configuration
/***** routines for survey variance ********/
double survey_variance (double a, double fsky);
double delP_SSC(double k, double a);
double w_mask(double theta_min);

/* mode coupling functions */

double alpha_c(double k1x,double k1y,double k2x,double k2y)    {

  double k1 = sqrt(k1x*k1x+k1y*k1y);
  double res  = 1.+(k1x*k2x+k1y*k2y)/(k1*k1);
  if(isnan(res)) {res=0.0;}
  return res;
}

double beta_c(double k1x,double k1y,double k2x,double k2y)
{
  double k1 = sqrt(k1x*k1x+k1y*k1y);
  double k2 = sqrt(k2x*k2x+k2y*k2y);
  double k1k2 = (k1x*k2x+k1y*k2y);
  double beta = k1k2*(k1*k1+k2*k2+2.*k1k2)/(2.*k1*k1*k2*k2);

  if (isnan(beta)) {beta=0.0;}
  return beta;
}

/* F & G Kernels */


double fs_2(double k1x, double k1y, double k2x, double k2y)
{
  double k1,k2,res;
  k1 = sqrt(k1x*k1x+k1y*k1y);
  k2 = sqrt(k2x*k2x+k2y*k2y);
  res = 5./7.+2./7.*pow((k1x*k2x+k1y*k2y)/(k1*k2),2.)+.5*((k1x*k2x+k1y*k2y)/(k1*k2))*(k1/k2+k2/k1);
  if(isnan(res)) {res=0.0;}
  //  printf("fs_2 %le\n",res);
  return res;
}

double gs_2(double k1x, double k1y,double k2x, double k2y)
{
  double k1,k2,res;
  k1 = sqrt(k1x*k1x+k1y*k1y);
  k2 = sqrt(k2x*k2x+k2y*k2y);
  res = 3./7.+4./7.*pow((k1x*k2x+k1y*k2y)/(k1*k2),2.)+.5*((k1x*k2x+k1y*k2y)/(k1*k2))*(k1/k2+k2/k1);
  if(isnan(res)) {res=0.0;}
  return res;
}

double fs_3(double k1x, double k1y, double k2x, double k2y, double k3x, double k3y)
{
  double term1,term2,term3;
  term1=7./54.*(alpha_c(k1x,k1y,k2x+k3x,k2y+k3y)*fs_2(k2x,k2y,k3x,k3y)
                + alpha_c(k2x,k2y,k1x+k3x,k1y+k3y)*fs_2(k1x,k1y,k3x,k3y)
                + alpha_c(k3x,k3y,k1x+k2x,k1y+k2y)*fs_2(k1x,k1y,k2x,k2y));

  term2=7./54.*(alpha_c(k1x+k2x,k1y+k2y,k3x,k3y)*gs_2(k1x,k1y,k2x,k2y)
                +     alpha_c(k1x+k3x,k1y+k3y,k2x,k2y)*gs_2(k1x,k1y,k3x,k3y)
                +     alpha_c(k2x+k3x,k2y+k3y,k1x,k1y)*gs_2(k2x,k2y,k3x,k3y));

  term3=4./54.*(beta_c(k1x,k1y,k2x+k3x,k2y+k3y)*gs_2(k2x,k2y,k3x,k3y)
                + beta_c(k2x,k2y,k1x+k3x,k1y+k3y)*gs_2(k1x,k1y,k3x,k3y)
                + beta_c(k3x,k3y,k1x+k2x,k1y+k2y)*gs_2(k1x,k1y,k2x,k2y));
  return term1+term2+term3;
}

double b_lin (double k1x, double k1y, double k2x, double k2y, double a)
{
  double k1,k2,k3,k3x,k3y;
  k3x = -k1x-k2x;
  k3y = -k1y-k2y;
  k1 = sqrt(k1x*k1x+k1y*k1y);
  k2 = sqrt(k2x*k2x+k2y*k2y);
  k3 = sqrt(k3x*k3x+k3y*k3y);
  return 2.*(fs_2(k1x,k1y,k2x,k2y)*p_lin(k1,a)*p_lin(k2,a)+fs_2(k2x,k2y,k3x,k3y)*p_lin(k2,a)*p_lin(k3,a)+fs_2(k3x,k3y,k1x,k1y)*p_lin(k3,a)*p_lin(k1,a));
}
double t_lin(double k1x, double k1y, double k2x, double k2y, double k3x, double k3y, double a)
{

  double k1,k2,k3,k4,k4x,k4y,p1,p2,p3,p4,p12,p13,p14;
  double k12,k13,k14;
  k4x = -k1x-k2x-k3x;
  k4y = -k1y-k2y-k3y;
  k1 = sqrt(k1x*k1x+k1y*k1y);
  k2 = sqrt(k2x*k2x+k2y*k2y);
  k3 = sqrt(k3x*k3x+k3y*k3y);
  k4 = sqrt(k4x*k4x+k4y*k4y);

  p1 = p_lin(k1,a);
  p2 = p_lin(k2,a);
  p3 = p_lin(k3,a);
  p4 = p_lin(k4,a);

  //  printf("t_lin %le %le %le %le %le %le %le %le %le\n",1./a-1.,k1,k2,k3,k4,p1,p2,p3,p4);
  k12 =sqrt(k1*k1+k2*k2+2.*k1x*k2x+2.*k1y*k2y);
  k13 =sqrt(k1*k1+k3*k3+2.*k1x*k3x+2.*k1y*k3y);
  k14 =sqrt(k1*k1+k4*k4+2.*k1x*k4x+2.*k1y*k4y);


  if(isnan(k12)) {k12=0.0; p12 = 0.0;}
  else {p12 = p_lin(k12,a);}
  if(isnan(k13)) {k13=0.0; p13 = 0.0;}
  else {p13 = p_lin(k13,a);}
  if(isnan(k14)) {k14=0.0; p14 = 0.0;}
  else {p14 = p_lin(k14,a);}

  // new test routine since p12 and p13 became nans
  if(isnan(p12) || isnan(p13) || isnan(p14)){
    printf("covariances_3D.c: NAN in tri_lin %le %le %le %le %le %le\n", k12,k13,k14,p12,p13,p14);
    exit(1);
  }
  //I introduced a return 0.0 in p_lin if k is below/above a certain range
  //--------------------------------------------------//
  double fs3_term,fs2_term;
  fs3_term=  fs_3(k1x,k1y,k2x,k2y,k3x,k3y)*p1*p2*p3 +fs_3(k1x,k1y,k2x,k2y,k4x,k4y)*p1*p2*p4 +fs_3(k1x,k1y,k3x,k3y,k4x,k4y)*p1*p3*p4 +fs_3(k2x,k2y,k3x,k3y,k4x,k4y)*p2*p3*p4;

  fs2_term=  fs_2(-k1x,-k1y,(k1x+k2x),(k1y+k2y))*fs_2(k3x,k3y,(k1x+k2x),(k1y+k2y))*p1*p12*p3
  +fs_2(-k1x,-k1y,(k1x+k2x),(k1y+k2y))*fs_2(k4x,k4y,(k1x+k2x),(k1y+k2y))*p1*p12*p4

  +fs_2(-k2x,-k2y,(k1x+k2x),(k1y+k2y))*fs_2(k3x,k3y,(k1x+k2x),(k1y+k2y))*p2*p12*p3
  +fs_2(-k2x,-k2y,(k1x+k2x),(k1y+k2y))*fs_2(k4x,k4y,(k1x+k2x),(k1y+k2y))*p2*p12*p4

  +fs_2(-k1x,-k1y,(k1x+k3x),(k1y+k3y))*fs_2(k2x,k2y,(k1x+k3x),(k1y+k3y))*p1*p13*p2
  +fs_2(-k1x,-k1y,(k1x+k3x),(k1y+k3y))*fs_2(k4x,k4y,(k1x+k3x),(k1y+k3y))*p1*p13*p4

  +fs_2(-k3x,-k3y,(k1x+k3x),(k1y+k3y))*fs_2(k2x,k2y,(k1x+k3x),(k1y+k3y))*p3*p13*p2
  +fs_2(-k3x,-k3y,(k1x+k3x),(k1y+k3y))*fs_2(k4x,k4y,(k1x+k3x),(k1y+k3y))*p1*p13*p4

  +fs_2(-k1x,-k1y,(k1x+k4x),(k1y+k4y))*fs_2(k2x,k2y,(k1x+k4x),(k1y+k4y))*p1*p14*p2
  +fs_2(-k1x,-k1y,(k1x+k4x),(k1y+k4y))*fs_2(k3x,k3y,(k1x+k4x),(k1y+k4y))*p1*p14*p3

  +fs_2(-k4x,-k4y,(k1x+k4x),(k1y+k4y))*fs_2(k2x,k2y,(k1x+k4x),(k1y+k4y))*p4*p14*p2
  +fs_2(-k4x,-k4y,(k1x+k4x),(k1y+k4y))*fs_2(k3x,k3y,(k1x+k4x),(k1y+k4y))*p4*p14*p3;

  //printf("fsterm %le %le %le %le %le %le %le\n",p1,p2,p3,p4,p12,p13,p14);
  return 6.*fs3_term+4.*fs2_term;
}
/*++++++++++++++++++++++++++++++++++++++++*
 * 1Halo Terms                            *
 *++++++++++++++++++++++++++++++++++++++++*/
double bi_1h(double k1,double k2, double k3, double a)
{
  return I0j(3,k1,k2,k3,0.,a);
}
double tri_1h(double k1, double k2, double k3, double k4,double a)
{
  return I0j(4,k1,k2,k3,k4,a);
}

/*++++++++++++++++++++++++++++++++++++++++*
 * 2Halo Terms                            					   *
 *++++++++++++++++++++++++++++++++++++++++*/
double tri_2h_13 (double k1x,double k1y,double k2x,double k2y,double k3x,double k3y,double a){
  double k1,k2,k3,k4,k4x,k4y;
  k4x = -k1x-k2x-k3x;
  k4y = -k1y-k2y-k3y;
  k1 = sqrt(k1x*k1x+k1y*k1y);
  k2 = sqrt(k2x*k2x+k2y*k2y);
  k3 = sqrt(k3x*k3x+k3y*k3y);
  k4 = sqrt(k4x*k4x+k4y*k4y);
  return (I1j(3,k1,k2,k3,a)*I1j(1,k4,0,0,a)*p_lin(k4,a) + I1j(3,k1,k2,k4,a)*I1j(1,k3,0,0,a)*p_lin(k3,a) + I1j(3,k1,k3,k4,a)*I1j(1,k2,0,0,a)*p_lin(k2,a) + I1j(3,k4,k2,k3,a)*I1j(1,k1,0,0,a)*p_lin(k1,a));
}

double tri_2h_22 (double k1x,double k1y,double k2x,double k2y,double k3x,double k3y,double a){
  double k1,k2,k3,k4,k4x,k4y,k12,k13,k14;
  k4x = -k1x-k2x-k3x;
  k4y = -k1y-k2y-k3y;
  k1 = sqrt(k1x*k1x+k1y*k1y);
  k2 = sqrt(k2x*k2x+k2y*k2y);
  k3 = sqrt(k3x*k3x+k3y*k3y);
  k4 = sqrt(k4x*k4x+k4y*k4y);
  k12 = sqrt(k1*k1 + 2.*k1x*k2x+2.*k1y*k2y+k2*k2);
  k13 = sqrt(k1*k1 + 2.*k1x*k3x+2.*k1y*k3y+k3*k3);
  k14 = sqrt(k1*k1 + 2.*k1x*k4x+2.*k1y*k4y+k4*k4);

  if(isnan(k12)) k12=0.0;
  if(isnan(k13)) k13=0.0;
  if(isnan(k14)) k14=0.0;

  return (I1j(2,k1,k2,0,a)*I1j(2,k3,k4,0,a)*p_lin(k12,a) + I1j(2,k1,k3,0,a)*I1j(2,k2,k4,0,a)*p_lin(k13,a) + I1j(2,k1,k4,0,a)*I1j(2,k3,k2,0,a)*p_lin(k14,a));
}

double tri_2h (double k1x,double k1y,double k2x,double k2y,double k3x,double k3y,double a){
  return tri_2h_13(k1x, k1y, k2x, k2y, k3x, k3y, a) + tri_2h_22 (k1x, k1y, k2x, k2y, k3x, k3y, a);
}

/*++++++++++++++++++++++++++++++++++++++++*
 * 3Halo Terms                            					   *
 *++++++++++++++++++++++++++++++++++++++++*/
double b_lin_mu(double k1,double k2, double m, double a){
  return b_lin(k1,0.,k2*m,k2*sqrt(1.-m*m),a);
}

double int_b_lin_mu(double m, void *params){
  double *ar = (double *) params;
  return b_lin(ar[0],0.,ar[1]*m,ar[1]*sqrt(1.-m*m),ar[2])/sqrt(1.-m*m);
}
double int_b_lin(double theta, void *params){
  double *ar = (double *) params;
  return b_lin(ar[0],0, ar[1]*cos(theta), ar[1]*sin(theta), ar[2]);
}
double b_lin_cov(double k1, double k2, double a){ //linear bispectrum, averaged over angle between k1 and k2
  double array[3] = {k1,k2,a};
  return int_gsl_integrate_medium_precision(int_b_lin,(void*)array,0.*M_PI,(2.*M_PI),NULL,1000)/(2.*M_PI);
}
double b_lin_cov_mu(double k1, double k2, double a){ //linear bispectrum, averaged over angle between k1 and k2
  double array[3] = {k1,k2,a};
  return int_gsl_integrate_medium_precision(int_b_lin_mu,(void*)array,-.999999999,0.999999999,NULL,1000)/2.;
}
double tri_3h_cov(double k1, double k2, double a){
  return 4.0*I1j(2,k1,k2,0,a)*I1j(1,k1,0,0,a)*I1j(1,k2,0,0,a)*b_lin_cov(k1,k2,a);
}
/********** survey variance ***************/
double C_survey_window(int l){
  static double *Cl = 0;
  if (Cl ==0){
    FILE *F1;
    F1 = fopen(covparams.C_FOOTPRINT_FILE,"r");
    if (F1 != NULL) { //covparams.C_FOOTPRINT_FILE exists, use healpix C_mask(l) to compute mask correlation function
      fclose(F1);
      int lbins = line_count(covparams.C_FOOTPRINT_FILE);
      Cl = create_double_vector(0,lbins-1);
      F1=fopen(covparams.C_FOOTPRINT_FILE,"r");
      for (int i = 0; i < lbins; i++){
        int tmp;
        double tmp2;
        fscanf(F1,"%d %le\n",&tmp, &tmp2);
        Cl[i] = tmp2;
        if(i>0) {Cl[i] /= Cl[0];}
      }
      Cl[0] = 1.;
      fclose(F1);
    }
    else{
      printf("covariances_3D.c:C_survey_window: covparams.C_FOOTPRINT_FILE =%s not found\nEXIT\n",covparams.C_FOOTPRINT_FILE);
      exit(1);
    }
  }
  return Cl[l];
}
double int_for_variance (double logk, void *params){
  double *ar = (double *) params;
  double k = exp(logk);
  double x = pow(4.0*ar[1],0.5)*k*chi(ar[0]); //theta_s*k*chi(a)
  return k*k/constants.twopi*p_lin(k,ar[0])*pow(2.*gsl_sf_bessel_J1(x)/x,2.0);
}

//curved sky, healpix window function
double sum_variance_healpix(double a){
  double res = 0.;
  double r = f_K(chi(a));
  for (int l = 0; l < 1000; l++){
    res+= (2.*l+1.)/(r*r)*C_survey_window(l)*p_lin((l+0.5)/r,a)/(4.*M_PI);
  }
  return res;
}

double survey_variance (double a, double fsky){
  static cosmopara C;
  static double FSKY = -42.;

  static double *table_SV;
  static double da =0, amin = 1./(1+4.01), amax =0.999;
  double aa,result,array[2];
  int i;
  if (recompute_cosmo3D(C) || FSKY != fsky){
    update_cosmopara(&C);
    FSKY = fsky;
    if (table_SV==0){
      table_SV  = create_double_vector(0, Ntable.N_a-1);
      amin  = limits.a_min;
      da = (amax - amin)/(Ntable.N_a-1.);
    }
    aa= amin;
    array[1] = fsky;
    FILE *F1;
    F1 = fopen(covparams.C_FOOTPRINT_FILE,"r");
    if (F1 != NULL) { //covparams.C_FOOTPRINT_FILE exists, use healpix C_mask(l) to compute mask correlation function
      fclose(F1);
      for (i=0; i<Ntable.N_a; i++, aa += da) {
         table_SV[i]=sum_variance_healpix(aa);
      }
    }
    else{ //covparams.C_FOOTPRINT_FILE doesn't exist, use analytic window
      for (i=0; i<Ntable.N_a; i++, aa += da) {
        array[0] = aa;
        result = int_gsl_integrate_high_precision(int_for_variance,(void*)array,log(1.e-6),log(1.e+6),NULL,2000);
        table_SV[i]=result;
      }
    }
  }
  return interpol(table_SV, Ntable.N_a, amin, amax, da,a, 1.0,1.0 );
}

/******* trispectrum terms for covariance (parallelogram) configurations *******/
double inner_P_bin(double theta, void *params){
  double *ar = (double *) params;
  double k = sqrt(ar[0]*ar[0]+2.0*cos(theta)*ar[0]*ar[1]+ar[1]*ar[1]);
  if(isnan(k)) k=0.;
  return p_lin(k,ar[2]);
}
double tri_2h_13_cov (double k1,double k2,double a){
   return (I1j(3,k1,k2,k2,a)*I1j(1,k1,0,0,a)*p_lin(k1,a)+I1j(3,k1,k1,k2,a)*I1j(1,k2,0,0,a)*p_lin(k2,a));
}

double tri_2h_22_cov (double k1, double k2, double a){
  double array[3];
  array[0] = k1; array[1] = k2;
  array[2] = a;
  return 2.0/M_PI*pow(I1j(2,k1,k2,0,a),2.0)*int_gsl_integrate_low_precision(inner_P_bin,(void*)array,0,M_PI,NULL,1000);
}
double tri_1h_cov(double k1, double k2, double a){
  return I0j(4,k1,k1,k2,k2,a);
}
double tri_2h_cov (double k1,double k2,double a){
  return tri_2h_22_cov(k1,k2, a)+ tri_2h_13_cov(k1,k2, a);
}
double inner_tri_lin_cov(double theta, void *params){
  double *ar = (double *) params;
  return t_lin(ar[0],0, ar[1]*cos(theta), ar[1]*sin(theta), -ar[0],0,ar[2]);
}
double tri_4h_cov(double k1, double k2, double a){
  double array[3] = {k1,k2,a};
  return pow(I1j(1,k1,0,0,a)*I1j(1,k2,0,0,a),2.0)*int_gsl_integrate_low_precision(inner_tri_lin_cov,(void*)array,0,M_PI,NULL,1000)/M_PI;
}

double tri_matter_cov(double k1, double k2, double a){
  return tri_1h_cov(k1,k2,a)+tri_2h_cov(k1,k2,a)+tri_3h_cov(k1,k2,a)+tri_4h_cov(k1,k2,a);
}
double tri_multih_cov(double k1, double k2, double a){
  return tri_2h_cov(k1,k2,a)+tri_3h_cov(k1,k2,a)+tri_4h_cov(k1,k2,a);
}
/*********************** super-sample covariance routines ***********************/

// MANUWARNING: shouldn't it be plin and not p2h?
double Delta_LD(double logk,void * params)
{
  (void)(params);
//  return log(p_2h(exp(logk),1.0));
  return log(Pdelta(exp(logk),0.999));
}

//linear dilation factor, cf. Eq. 27 in http://arxiv.org/pdf/1401.0385v2.pdf
// MANUWARNING: do we want to use the equation of Chiang+14, which uses Plin and not P2h?
// Eq 4.18, 4.32 in Chiang+14 http://arxiv.org/abs/1403.3411v2
double LD_term(double k)
{
  static cosmopara C;

  static double *table_LD;
  static double dlogk = .0, logkmin = 1.0,logkmax = 1.0;
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);
    double klog, abserr,result;
    int i;

    gsl_function F;
    F.function = &Delta_LD;
    F.params = 0;

    if (table_LD==0){
      table_LD  = create_double_vector(0, Ntable.N_k_nlin-1);
      logkmin = log(limits.k_min_cH0);
      logkmax = log(limits.k_max_cH0);
      dlogk = (logkmax - logkmin)/(Ntable.N_k_nlin);
    }
    klog = logkmin;
    for (i=0; i<Ntable.N_k_nlin; i++, klog += dlogk) {
      gsl_deriv_central (&F,klog, 0.1*klog, &result, &abserr);
      table_LD[i]=(result/3.+1.);
    }
  }
  return -interpol(table_LD, Ntable.N_k_nlin, logkmin, logkmax, dlogk,log(k), 1.0,1.0);
}

// the HSV on the 2h term is not included (only beat coupling), as in Li+14 http://arxiv.org/pdf/1401.0385v2.pdf
// as discussed in Chiang+14 (after Eq 4.32) http://arxiv.org/abs/1403.3411v2
// the inclusion of the HSV on P2h is ambiguous, requires b2, and only leads to percent differences.
double delP_SSC(double k, double a){
  return (68./21.+LD_term(k))*Pdelta(k,a);
//  return (68./21.+LD_term(k))*p_2h(k,a)+ I12_SSC(k,a);
}


/*********************** super-sample covariance: linear term only ***********************/


double DeltaLin_LD(double logk, void * params)
{
   (void)(params);
   return log(p_lin(exp(logk),1.0));
}

// Eq 4.18, 4.32 in Chiang+14 http://arxiv.org/abs/1403.3411v2
double LDlin_term(double k)
{
   static cosmopara C;

   static double *table_LD;
   static double dlogk = .0, logkmin = 1.0,logkmax = 1.0;
   if (recompute_cosmo3D(C)){
      update_cosmopara(&C);
      double klog, abserr,result;
      int i;

      gsl_function F;
      F.function = &DeltaLin_LD;
      F.params = 0;

      if (table_LD==0){
         table_LD  = create_double_vector(0, Ntable.N_k_nlin-1);
         logkmin = log(limits.k_min_cH0);
         logkmax = log(limits.k_max_cH0);
         dlogk = (logkmax - logkmin)/(Ntable.N_k_nlin);
      }
      klog = logkmin;
      for (i=0; i<Ntable.N_k_nlin; i++, klog += dlogk) {
         gsl_deriv_central (&F,klog, 0.1*klog, &result, &abserr);
         table_LD[i]=(result/3.+1.);
      }
   }
   return - interpol(table_LD, Ntable.N_k_nlin, logkmin, logkmax, dlogk,log(k), 1.0,1.0);
}

double delPlin_SSC(double k, double a){
   double res = 68./21.  + LDlin_term(k);
   res *= p_lin(k, a);
   return res;
}

double w_mask(double theta_min){
  static int NTHETA = 0;
  static double *w_vec =0;
  int i,l;
  if (like.theta_min == 0 || like.Ntheta < 1){
    printf("covariances_real_binned.c:w_mask: like.theta_min or like.Ntheta not initialized\nEXIT\n");
    exit(1);
  }
  if (w_vec ==0){
    w_vec = create_double_vector(0,like.Ntheta-1);
    NTHETA = like.Ntheta;
    FILE *F1;
    F1 = fopen(covparams.C_FOOTPRINT_FILE,"r");
    if (F1 != NULL) { //covparams.C_FOOTPRINT_FILE exists, use healpix C_mask(l) to compute mask correlation function
      fclose(F1);
      int lbins = line_count(covparams.C_FOOTPRINT_FILE);
      if (lbins < 2.*M_PI/like.vtmin){
        printf("\n\n\ncovariances_3D.c:w_mask:WARNING\nCompute covariance without accounting to survey boundary effect on shape/shot noise:\n");
        printf("Survey mask power spectrum in %s tabulated only to l_max = %d\n",covparams.C_FOOTPRINT_FILE,lbins-1);
        printf("This is insufficient to calculate masking effect at theta_min = %.2f arcmin\n\n\n",like.vtmin/constants.arcmin);
        for (i = 0; i<NTHETA; i ++){
          w_vec[i] = 1.0;
        }
      }
      else{
        double *Cl;
        Cl = create_double_vector(0,lbins-1);
        F1=fopen(covparams.C_FOOTPRINT_FILE,"r");
        for (int i = 0; i < lbins; i++){
          int tmp;
          double tmp2;
          fscanf(F1,"%d %le\n",&tmp, &tmp2);
          Cl[i] = tmp2;
        }
        fclose(F1);

        printf("\nTabulating w_mask(theta) from mask power spectrum %s\n",covparams.C_FOOTPRINT_FILE);
        double **Pl, *xmin, *xmax, *Pmin, *Pmax;
        Pl =create_double_matrix(0, like.Ntheta-1, 0, lbins);
        xmin= create_double_vector(0, like.Ntheta-1);
        xmax= create_double_vector(0, like.Ntheta-1);
        Pmin= create_double_vector(0, lbins);
        Pmax= create_double_vector(0, lbins);
        for(i=0; i<like.Ntheta ; i++){
          xmin[i]=cos(like.theta_min[i]);
          xmax[i]=cos(like.theta_min[i+1]);
        }
        for (i = 0; i < like.Ntheta; i++){
          gsl_sf_legendre_Pl_array(lbins, xmin[i],Pmin);
          gsl_sf_legendre_Pl_array(lbins, xmax[i],Pmax);
          Pl[i][0] = 1.0/(4.*M_PI);
          for (int l = 1; l < lbins; l ++){
            Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/(xmin[i]-xmax[i]);
          }
        }
        free_double_vector(xmin,0,like.Ntheta-1);
        free_double_vector(xmax,0,like.Ntheta-1);
        free_double_vector(Pmin,0,lbins);
        free_double_vector(Pmax,0,lbins);
        for (i = 0; i < NTHETA; i++){
          w_vec[i] =0.;
          for (l = 0; l < lbins; l++){
            w_vec[i]+=Cl[l]*Pl[i][l];
          }
          printf("w_mask[%d] = %e\n",i, w_vec[i]);
        }
        free_double_vector(Cl,0,lbins-1);
        free_double_matrix(Pl, 0, like.Ntheta-1, 0, lbins);
      }
    }
    else{ //covparams.C_FOOTPRINT_FILE does not exit, ignore boundary effects
      printf("covparams.C_FOOTPRINT_FILE = %s not found\nNo footprint boundary effect correction applied to shape/shot noise\n",covparams.C_FOOTPRINT_FILE);
      for (i = 0; i<NTHETA; i ++){
        w_vec[i] = 1.0;
      }
    }
  }
  i = 0;
  while(like.theta_min[i]< theta_min){
    i ++;
  }
  return w_vec[i];
}
