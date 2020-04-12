/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by CosmoLike developers
******************************************************************************/

/************************* covariance routines for angular power spectra *********************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
//     No T(l1,l2) look-up tables; most useful for small number of l-bins
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/************************* routines for angular trispectrum terms ************************/
/********************           AUTO-COVARIANCE BLOCKS              **********************/
/* Note: Gaussian covs are defined in each covariances_real routine */
/*       Only NG covs are here */

double cov_NG_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4);
double cov_NG_gl_gl_tomo(double l1,double l2, int z1l, int z1s, int z2l, int z2s);
double cov_NG_cl_cl_tomo(double l1,double l2, int z1, int z2, int z3, int z4);

/********************          CROSS-COVARIANCE BLOCKS              **********************/
double cov_NG_gl_shear_tomo(double l1,double l2, int zl, int zs, int z1, int z2);//zl,zs g-g lensing bins; z3,z4 shear bins
double cov_NG_cl_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4); //z1,z2 clustering bins; z3,z4 shear bins
double cov_NG_cl_gl_tomo(double l1,double l2, int z1, int z2, int zl, int zs); //z1,z2 clustering bins; zl,zs g-g lensing bins

int cNG = 1;
/************** shear x shear routines ***************/
double inner_project_tri_cov_shear_shear_tomo(double a,void *params)
{
  cNG = covparams.cng;
  double k1,k2,fK,weights,res = 0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_kappa(a,fK,ar[2])*W_kappa(a,fK,ar[3])*W_kappa(a,fK,ar[4])*W_kappa(a,fK,ar[5])*dchi_da(a);
  if (weights >0.){
    if (cNG) {res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);}
    res += delP_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}

double cov_NG_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  double a1,a2,array[7];
  int zmin;
  zmin = (z2 < z1 ? z2 : z1);
  zmin = (z3 < zmin ? z3 : zmin);
  zmin = (z4 < zmin ? z4 : zmin);
  a1 = amin_source(zmin);
  a2 = amax_source(zmin);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) z1;
  array[3] = (double) z2;
  array[4] = (double) z3;
  array[5] = (double) z4;
  array[6] = survey.area/41253.0;
  return int_gsl_integrate_low_precision(inner_project_tri_cov_shear_shear_tomo,(void*)array,a1,a2,NULL,1000);
}

/***** gg-lensing x gg-lensing routines ****/
double bgal_a(double a, double nz){
  return gbias.b1_function(1./a-1.,(int)nz);
}
double inner_project_tri_cov_gl_gl_tomo(double a,void *params)
{
  cNG = covparams.cng;
  double k1,k2,fK,weights,res = 0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_gal(a,ar[2])*W_kappa(a,fK,ar[3])*W_gal(a,ar[4])*W_kappa(a,fK,ar[5])*dchi_da(a);
  if (weights >0.){
    if (cNG) {res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);}
    res +=(delP_SSC(k1,a)-bgal_a(a,ar[2])*Pdelta(k1,a))*(delP_SSC(k2,a)-bgal_a(a,ar[4])*Pdelta(k2,a))*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}

double cov_NG_gl_gl_tomo(double l1,double l2, int z1l, int z1s, int z2l, int z2s){
  double a1,a2,array[7];
  int zmax, zmin;
  zmax = (z2l > z1l ? z2l : z1l);
  a1 = amin_lens(zmax);
  zmin = (z2l < z1l ? z2l : z1l);
  a2 = amax_lens(zmin);
//  printf("%le %le\n",l1,l2);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) z1l;
  array[3] = (double) z1s;
  array[4] = (double) z2l;
  array[5] = (double) z2s;
  array[6] = survey.area/41253.0;

  return int_gsl_integrate_low_precision(inner_project_tri_cov_gl_gl_tomo,(void*)array,a1,a2,NULL,1000);
}

/***** clustering x clustering routines ****/

double inner_project_tri_cov_cl_cl_tomo(double a,void *params)
{
  cNG = covparams.cng;
  double k1,k2,fK,weights,res = 0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_gal(a,ar[2])*W_gal(a,ar[3])*W_gal(a,ar[4])*W_gal(a, ar[5])*dchi_da(a);
  if (weights >0.){
    if (cNG) {res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);}
    res += (delP_SSC(k1,a)-2.*bgal_a(a,ar[2])*Pdelta(k1,a))*(delP_SSC(k2,a)-2.*bgal_a(a,ar[4])*Pdelta(k2,a))*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}


double cov_NG_cl_cl_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  double a1,a2,array[7];
  int zmin, zmax;
  zmax = (z2 > z1 ? z2 : z1);
  zmax = (z3 > zmax ? z3 : zmax);
  zmax = (z4 > zmax ? z4 : zmax);
  a1 = amin_lens(zmax);
  zmin = (z2 < z1 ? z2 : z1);
  zmin = (z3 < zmin ? z3 : zmin);
  zmin = (z4 < zmin ? z4 : zmin);
  a2 = amax_lens(zmin);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) z1;
  array[3] = (double) z2;
  array[4] = (double) z3;
  array[5] = (double) z4;
  array[6] = survey.area/41253.0;

  return int_gsl_integrate_low_precision(inner_project_tri_cov_cl_cl_tomo,(void*)array,a1,a2,NULL,1000);
}

/************** g-g lensing x shear routines ***************/

double inner_project_tri_cov_gl_shear_tomo(double a,void *params)
{
  cNG = covparams.cng;
  double k1,k2,fK,weights,res = 0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_gal(a, ar[2])*W_kappa(a,fK,ar[3])*W_kappa(a,fK,ar[4])*W_kappa(a,fK,ar[5])*dchi_da(a);
  if (weights >0.){
    if (cNG) {res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);}
    res += (delP_SSC(k1,a)-bgal_a(a,ar[2])*Pdelta(k1,a))*(delP_SSC(k2,a))*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}

double cov_NG_gl_shear_tomo(double l1,double l2, int zl, int zs, int z3, int z4){ //zl,zs g-g lensing bins; z3,z4 shear bins
  double a1,a2,array[7];
  a1 = amin_lens(zl);
  a2 = amax_lens(zl);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) zl;
  array[3] = (double) zs;
  array[4] = (double) z3;
  array[5] = (double) z4;
  array[6] = survey.area/41253.0;
  return int_gsl_integrate_low_precision(inner_project_tri_cov_gl_shear_tomo,(void*)array,a1,a2,NULL,1000);
}

/************** clustering x shear routines ***************/

double inner_project_tri_cov_cl_shear_tomo(double a,void *params)
{
  cNG = covparams.cng;
  double k1,k2,fK,weights,res = 0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_gal(a,ar[2])*W_gal(a, ar[3])*W_kappa(a,fK,ar[4])*W_kappa(a,fK,ar[5])*dchi_da(a);
  if (weights >0.){
    if (cNG) {res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);}
    // printf("res cNG:%lg , ", res);
    res += (delP_SSC(k1,a)-2.*bgal_a(a,ar[2])*Pdelta(k1,a))*(delP_SSC(k2,a))*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  // printf("res +SSC, weights:%lg, %lg\n", res, weights);
  res *= weights;
  return res;
}

double cov_NG_cl_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){ //z1,z2 clustering bins; z3,z4 shear bins
  double a1,a2,array[7];
  int zmin, zmax;
  zmax = (z1 > z2 ? z1 : z2);
  zmin = (z1 < z2 ? z1 : z2);
  a1 = amin_lens(zmax);
  a2 = amax_lens(zmin);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) z1;
  array[3] = (double) z2;
  array[4] = (double) z3;
  array[5] = (double) z4;
  array[6] = survey.area/41253.0;
  return int_gsl_integrate_low_precision(inner_project_tri_cov_cl_shear_tomo,(void*)array,a1,a2,NULL,1000);
}

/************** clustering x shear routines ***************/

double inner_project_tri_cov_cl_gl_tomo(double a,void *params)
{
  cNG = covparams.cng;
  double k1,k2,fK,weights,res = 0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_gal(a,ar[2])*W_gal(a,ar[3])*W_gal(a,ar[4])*W_kappa(a,fK,ar[5])*dchi_da(a);
  if (weights >0.){
    if (cNG) {res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);}
    res += (delP_SSC(k1,a)-2.*bgal_a(a,ar[2])*Pdelta(k1,a))*(delP_SSC(k2,a)-bgal_a(a,ar[4])*Pdelta(k2,a))*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}

double cov_NG_cl_gl_tomo(double l1,double l2, int z1, int z2, int zl, int zs){ //z1,z2 clustering bins; zl,zs g-g lensing bins
  double a1,a2,array[7];
  a1 = amin_lens(z1);
  a2 = amax_lens(z2);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) z1;
  array[3] = (double) z2;
  array[4] = (double) zl;
  array[5] = (double) zs;
  array[6] = survey.area/41253.0;
  return int_gsl_integrate_low_precision(inner_project_tri_cov_cl_gl_tomo,(void*)array,a1,a2,NULL,1000);
}

// Binned NG funcs


double filter_cov_fourier(double l1, double l2, double lmax, double lpivot) {
  long i,j;
  double W;
  double l_interval = lmax - lpivot;
  if(l_interval<=0) {return 1.;}

  if(l1<=lpivot){
    W = 1.;
  }
  else if(l1>=lmax){
    W = 0.;
  }
  else{
    W = (lmax - l1) / l_interval - 1./(2.*M_PI) * sin(2.*(lmax - l1)*M_PI/l_interval);
  }

  if(l2<=lpivot){
    W *= 1.;
  }
  else if(l2>=lmax){
    W *= 0.;
  }
  else{
    W *= (lmax - l2) / l_interval - 1./(2.*M_PI) * sin(2.*(lmax - l2)*M_PI/l_interval);
  }
  return W;
}




/**************** look-up tables for covariance  *********************/
double bin_cov_NG_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 20;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(5.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=log(cov_NG_shear_shear_tomo(ll1,ll2,z1,z2,z3,z4));
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = exp(interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0));}
    res *= filter_cov_fourier(llog1, llog2, logsmax, log(4.e4));
  return res;
}
double bin_cov_NG_gl_gl_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(5.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]= cov_NG_gl_gl_tomo(ll1,ll2,z1,z2,z3,z4);
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}
double bin_cov_NG_cl_cl_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(5.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=cov_NG_cl_cl_tomo(ll1,ll2,z1,z2,z3,z4);
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
  }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}
double bin_cov_NG_cl_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(5.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=cov_NG_cl_shear_tomo(ll1,ll2,z1,z2,z3,z4);
        // printf("cov_NG_cl_shear_tomo(%lg,%lg),%lg\n",ll1,ll2, table[i][j]);
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}
double bin_cov_NG_cl_gl_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(5.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=cov_NG_cl_gl_tomo(ll1,ll2,z1,z2,z3,z4);
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}
double bin_cov_NG_gl_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(5.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=cov_NG_gl_shear_tomo(ll1,ll2,z1,z2,z3,z4);
        // printf("cov_NG_gl_shear_tomo(%lg,%lg),%lg\n",ll1,ll2, cov_NG_gl_shear_tomo(ll1,ll2,z1,z2,z3,z4));
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}
