/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by CosmoLike developers
******************************************************************************/

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <time.h>
#include <stdio.h>

void cov_G_cl_cl_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int z3, int z4, double **cov_g_interp);
void cov_NG_cl_cl_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int z3, int z4, double **cov_ng_interp);
void cov_G_cl_gl_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int zl, int zs, double **cov_g_interp);
void cov_NG_cl_gl_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int zl, int zs, double **cov_ng_interp);
void cov_G_cl_shear_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int z3, int z4, int pm, double **cov_g_interp);
void cov_NG_cl_shear_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int z3, int z4, int pm, double **cov_ng_interp);
void cov_G_gl_gl_real_fft_binned(double *theta, int Ntheta, double *dtheta, int zl1, int zs1, int zl2, int zs2, double **cov_g_interp);
void cov_NG_gl_gl_real_fft_binned(double *theta, int Ntheta, double *dtheta, int zl1, int zs1, int zl2, int zs2, double **cov_ng_interp);
void cov_G_gl_shear_real_fft_binned(double *theta, int Ntheta, double *dtheta, int zl, int zs, int z3, int z4, int pm, double **cov_g_interp);
void cov_NG_gl_shear_real_fft_binned(double *theta, int Ntheta, double *dtheta, int zl, int zs, int z3, int z4, int pm, double **cov_ng_interp);
void cov_G_shear_shear_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int z3, int z4, int pm1, int pm2, double **cov_g_interp);
void cov_NG_shear_shear_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int z3, int z4, int pm1, int pm2, double **cov_ng_interp);

void cov_G_real_fft_bin_template(double *theta, int Ntheta, double **cov_g_interp, config my_config, double *N, int *array, double *task_type);
void cov_NG_real_fft_bin_template(double *theta, int Ntheta, double **cov_ng_interp, config my_config, int *array, double *task_type);

double func_for_cov_G_shear(double l, int *ar);
double func_for_cov_NG_shear(double l1, double l2, int *ar);
double func_for_cov_G_cl(double l, int *ar);
double func_for_cov_NG_cl(double l1, double l2, int *ar);
double func_for_cov_G_cl_gl(double l, int *ar);
double func_for_cov_NG_cl_gl(double l1, double l2, int *ar);
double func_for_cov_G_cl_shear(double l, int *ar);
double func_for_cov_NG_cl_shear(double l1, double l2, int *ar);
double func_for_cov_G_gl(double l, int *ar);
double func_for_cov_NG_gl(double l1, double l2, int *ar);
double func_for_cov_G_gl_shear(double l, int *ar);
double func_for_cov_NG_gl_shear(double l1, double l2, int *ar);

void cov_G_shear_shear_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int z3, int z4, int pm1, int pm2, double **cov_g_interp) {

  // pure noise term
  double N[Ntheta],res = 0.;
  int array[5];
  array[1] = z1; array[2] = z2; array[3] = z3; array[4] = z4;
  
  int i,j;
  for(i=0;i<Ntheta;i++) {N[i] = 0.;}

  if (z1 ==z3 && z2 ==z4 && pm1 == pm2){ //&& pm1 == pm2 required as C+- doesn't have the diagonal shot noise term
    for(i=0;i<Ntheta;i++) {
      N[i] += pow(survey.sigma_e,4.0)/(M_PI*(theta[i+1]*theta[i+1]-theta[i]*theta[i])*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    } //(number of galaxy pairs in the survey contributing to annulus of width Dtheta centered at theta1)^-1
  }
  if (z1 ==z4 && z2 ==z3 && pm1 == pm2){ //&& pm1 == pm2 required as C+- doesn't have the diagonal shot noise term
    for(i=0;i<Ntheta;i++) {
      N[i] += pow(survey.sigma_e,4.0)/(M_PI*(theta[i+1]*theta[i+1]-theta[i]*theta[i])*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    }
  }
  for(i=0;i<Ntheta;i++) {
    if(N[i]) {N[i]/= w_mask(theta[i]); printf("N[i]: %lg\n", N[i]);}
  }

  // get big cov_G matrix
  config my_config;
  if (pm1 == 1){my_config.l1 = -0.5;}
  else{my_config.l1 = 3.5;}
  
  if (pm2 == 1){my_config.l2 = -0.5;}
  else{my_config.l2 = 3.5;}
  
  // my_config.nu1 = 1.01-my_config.l1/2.;
  // my_config.nu2 = 1.01-my_config.l2/2.;

  my_config.nu1 = 1.01;
  my_config.nu2 = 1.01;


  double task_type[] = {2,2};
  cov_G_real_fft_bin_template(theta, Ntheta, cov_g_interp, my_config, N, array, task_type);
}

void cov_NG_shear_shear_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int z3, int z4, int pm1, int pm2, double **cov_ng_interp) {
  int array[5];
  array[1] = z1; array[2] = z2; array[3] = z3; array[4] = z4;
  int i,j;

  // get big cov_NG matrix
  config my_config;
  if (pm1 == 1){my_config.l1 = -0.5;}
  else{my_config.l1 = 3.5;}
  
  if (pm2 == 1){my_config.l2 = -0.5;}
  else{my_config.l2 = 3.5;}
  
  my_config.nu1 = 1.01;
  my_config.nu2 = 1.01;

  double task_type[] = {2,2};
  cov_NG_real_fft_bin_template(theta, Ntheta, cov_ng_interp, my_config, array, task_type);

}

void cov_G_cl_cl_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int z3, int z4, double **cov_g_interp) {

  // pure noise term
  double N[Ntheta],res = 0.;
  int array[5];
  array[1] = z1; array[2] = z2; array[3] = z3; array[4] = z4;

  int i,j;
  for(i=0;i<Ntheta;i++) {N[i] = 0.;}


  if (z1 ==z3 && z2 ==z4){
    for(i=0;i<Ntheta;i++) {
      N[i] += 1./(M_PI*(theta[i+1]*theta[i+1]-theta[i]*theta[i])*nlens(z1)*nlens(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor); //(number of galaxy pairs in the survey contributing to annulus of width Dtheta centered at theta1)^-1
    }
  }
  if (z1 ==z4 && z2 ==z3){
    for(i=0;i<Ntheta;i++) {
      N[i] += 1./(M_PI*(theta[i+1]*theta[i+1]-theta[i]*theta[i])*nlens(z1)*nlens(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    }
  } 
  for(i=0;i<Ntheta;i++) {
    if(N[i]) N[i]/= w_mask(theta[i]);
  }

  // get big cov_G matrix
  config my_config;
  my_config.l1 = -0.5;
  my_config.l2 = -0.5;

  my_config.nu1 = 1.01;
  my_config.nu2 = 1.01;

  double task_type[] = {0,0};
  cov_G_real_fft_bin_template(theta, Ntheta, cov_g_interp, my_config, N, array, task_type);
}

void cov_NG_cl_cl_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int z3, int z4, double **cov_ng_interp) {
  int array[5];
  array[1] = z1; array[2] = z2; array[3] = z3; array[4] = z4;
  int i,j;

  // get big cov_NG matrix
  config my_config;
  my_config.l1 = -0.5;
  my_config.l2 = -0.5;
  my_config.nu1 = 1.01;
  my_config.nu2 = 1.01;
  double task_type[] = {0,0};
  cov_NG_real_fft_bin_template(theta, Ntheta, cov_ng_interp, my_config, array, task_type);
}

void cov_G_cl_gl_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int zl, int zs, double **cov_g_interp) {

  double N[Ntheta],res = 0.;
  int array[5];
  array[1] = z1; array[2] = z2; array[3] = zl; array[4] = zs;
  
  int i,j;

  // get big cov_G matrix
  config my_config;
  my_config.l1 = -0.5;
  my_config.l2 = 1.5;
  my_config.nu1 = 1.01;
  my_config.nu2 = 1.01;

  double task_type[] = {0,1};
  cov_G_real_fft_bin_template(theta, Ntheta, cov_g_interp, my_config, N, array, task_type);
}

void cov_NG_cl_gl_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int zl, int zs, double **cov_ng_interp) {
  int array[5];
  array[1] = z1; array[2] = z2; array[3] = zl; array[4] = zs;
  int i,j;

  // get big cov_NG matrix
  config my_config;
  my_config.l1 = -0.5;
  my_config.l2 = 1.5;
  my_config.nu1 = 1.01;
  my_config.nu2 = 1.01;

  double task_type[] = {0,1};
  cov_NG_real_fft_bin_template(theta, Ntheta, cov_ng_interp, my_config, array, task_type);
}

void cov_G_cl_shear_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int z3, int z4, int pm, double **cov_g_interp) {

  double N[Ntheta],res = 0.;
  int array[5];
  array[1] = z1; array[2] = z2; array[3] = z3; array[4] = z4;
  
  int i,j;

  // get big cov_G matrix
  config my_config;
  my_config.l1 = -0.5;
  if (pm == 1){my_config.l2 = -0.5;}
  else{my_config.l2 = 3.5;}
  my_config.nu1 = 1.01;
  my_config.nu2 = 1.01;

  double task_type[] = {0,2};
  cov_G_real_fft_bin_template(theta, Ntheta, cov_g_interp, my_config, N, array, task_type);
}

void cov_NG_cl_shear_real_fft_binned(double *theta, int Ntheta, double *dtheta, int z1, int z2, int z3, int z4, int pm, double **cov_ng_interp) {
  int array[5];
  array[1] = z1; array[2] = z2; array[3] = z3; array[4] = z4;
  int i,j;

  // get big cov_NG matrix
  config my_config;
  my_config.l1 = -0.5;
  if (pm == 1){my_config.l2 = -0.5;}
  else{my_config.l2 = 3.5;}
  my_config.nu1 = 1.01;
  my_config.nu2 = 1.01;

  double task_type[] = {0,2};
  cov_NG_real_fft_bin_template(theta, Ntheta, cov_ng_interp, my_config, array, task_type);
}

void cov_G_gl_gl_real_fft_binned(double *theta, int Ntheta, double *dtheta, int zl1, int zs1, int zl2, int zs2, double **cov_g_interp) {

  // pure noise term
  double N[Ntheta],res = 0.;
  int array[5];
  array[1] = zl1; array[2] = zs1; array[3] = zl2; array[4] = zs2;
  
  int i,j;
  for(i=0;i<Ntheta;i++) {N[i] = 0.;}

  if (zl1 ==zl2 && zs1 ==zs2){
    for(i=0;i<Ntheta;i++) {
      N[i] += pow(survey.sigma_e,2.0)/(2.0*M_PI*(theta[i+1]*theta[i+1]-theta[i]*theta[i])*nlens(zl1)*nsource(zs2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
      if(N[i]) N[i]/= w_mask(theta[i]);
    }
  }

  // get big cov_G matrix
  config my_config;
  my_config.l1 = 1.5;
  my_config.l2 = 1.5;
  my_config.nu1 = 1.01;
  my_config.nu2 = 1.01;

  double task_type[] = {1,1};
  cov_G_real_fft_bin_template(theta, Ntheta, cov_g_interp, my_config, N, array, task_type);
}

void cov_NG_gl_gl_real_fft_binned(double *theta, int Ntheta, double *dtheta, int zl1, int zs1, int zl2, int zs2, double **cov_ng_interp) {
  int array[5];
  array[1] = zl1; array[2] = zs1; array[3] = zl2; array[4] = zs2;
  int i,j;

  // get big cov_NG matrix
  config my_config;
  my_config.l1 = 1.5;
  my_config.l2 = 1.5;
  my_config.nu1 = 1.01;
  my_config.nu2 = 1.01;

  double task_type[] = {1,1};
  cov_NG_real_fft_bin_template(theta, Ntheta, cov_ng_interp, my_config, array, task_type);
}

void cov_G_gl_shear_real_fft_binned(double *theta, int Ntheta, double *dtheta, int zl, int zs, int z3, int z4, int pm, double **cov_g_interp) {

  double N[Ntheta],res = 0.;
  int array[5];
  array[1] = zl; array[2] = zs; array[3] = z3; array[4] = z4;
  
  int i,j;

  // get big cov_G matrix
  config my_config;
  my_config.l1 = 1.5;
  
  if (pm == 1){my_config.l2 = -0.5;}
  else{my_config.l2 = 3.5;}
  my_config.nu1 = 1.01;
  my_config.nu2 = 1.01;

  double task_type[] = {1,2};
  cov_G_real_fft_bin_template(theta, Ntheta, cov_g_interp, my_config, N, array, task_type);
}

void cov_NG_gl_shear_real_fft_binned(double *theta, int Ntheta, double *dtheta, int zl, int zs, int z3, int z4, int pm, double **cov_ng_interp) {
  int array[5];
  array[1] = zl; array[2] = zs; array[3] = z3; array[4] = z4;
  int i,j;

  // get big cov_NG matrix
  config my_config;
  my_config.l1 = 1.5;
  
  if (pm == 1){my_config.l2 = -0.5;}
  else{my_config.l2 = 3.5;}
  my_config.nu1 = 1.01;
  my_config.nu2 = 1.01;

  double task_type[] = {1,2};
  cov_NG_real_fft_bin_template(theta, Ntheta, cov_ng_interp, my_config, array, task_type);
}




/*********** Templates *************************************/

void cov_G_real_fft_bin_template(double *theta, int Ntheta, double **cov_g_interp, config my_config, double *N, int *array, double *task_type) {

  my_config.sym_Flag = 0;
  my_config.c_window_width = 0.25;

  double res = 0.;
  int i,j;
  // get input func log-sampled in ell
  const int Nsub = 20; // split each (theta[i],theta[i+1]) range into Nsub sub-ranges, i.e. fill in (Nsub-1) sampling points within the range

  // clock_t start_t, end_t;
  // start_t = clock();

  double ellmax = 50000.;
  double ell_first = 1./theta[Ntheta]; // include the last bin edge
  double ell_last = 1./theta[0];
  double dlnell = log(theta[1]/theta[0])/Nsub;
  int ind_first = floor(log(ell_first)/dlnell); // index of ell_first in log-sampled array (from range [1,ellmax]) containing our orginal desired array (1/theta)
  int Nell = floor(log(ellmax/ell_first)/dlnell) + ind_first + 1; // total sampling points in range [1,ellmax]
  double ell[Nell], func[Nell];
  // printf("ell_first:%lg\n", ell_first);
  // printf("Nell:%d\n", Nell);

  // end_t = clock();
  // printf("Total time taken by CPU: %f\n", (double)(end_t - start_t) / CLOCKS_PER_SEC  );

  double noise_factor = 0.;
  double (*func_for_cov_G)(double, int*);
  if(task_type[0]==0) {
    if(task_type[1]==0) {func_for_cov_G = &func_for_cov_G_cl; noise_factor=1.;}
    else if(task_type[1]==1) {func_for_cov_G = &func_for_cov_G_cl_gl;}
    else if(task_type[1]==2) {func_for_cov_G = &func_for_cov_G_cl_shear;}
  }
  else if(task_type[0]==1) {
    if(task_type[1]==1) {func_for_cov_G = &func_for_cov_G_gl; noise_factor=1.;}
    else if(task_type[1]==2) {func_for_cov_G = &func_for_cov_G_gl_shear;}
  }
  else {
    if(task_type[1]==2) {func_for_cov_G = &func_for_cov_G_shear; noise_factor=2.;}
  }


  for(i=0; i<Nell; i++) {
    ell[i] = exp(dlnell * (i-ind_first)) * ell_first;
    func[i] = func_for_cov_G(ell[i], array);
  }

  // end_t = clock();
  // printf("Total time taken by CPU: %f\n", (double)(end_t - start_t) / CLOCKS_PER_SEC  );


  // extrapolate
  // const int N_extra = 800-246;
  const int N_extra = 800;
  int Ntot = N_extra*2+Nell;
  double large_ell[Ntot], large_fl[Ntot];
  extrap_log_linear_cfastcov(ell, Nell, N_extra, large_ell);
  extrap_log_linear_cfastcov(func, Nell, N_extra, large_fl);

  double **fk1k2;
  fk1k2 = malloc(Ntot*sizeof(double *));
  for(i=0;i<Ntot;i++) {
    fk1k2[i] = malloc(Ntot*sizeof(double));
  }
  mk_diag_g_to_ng(large_fl, Ntot, dlnell, fk1k2);

  // end_t = clock();
  // printf("Total time taken by CPU: %f\n", (double)(end_t - start_t) / CLOCKS_PER_SEC  );

  // do integration
  double theta_fine[Ntot];
  double **result;
  result = malloc(Ntot*sizeof(double *));
  for(i=0;i<Ntot;i++) {
    result[i] = malloc(Ntot*sizeof(double));
  }
  // twobessel(large_ell, large_ell, fk1k2, Ntot, Ntot, &my_config, theta_fine, theta_fine, result);
  double smooth_dlnr = dlnell*Nsub;
  printf("Ntot,%d\n",Ntot);
  twobessel_binave(large_ell, large_ell, fk1k2, Ntot, Ntot, &my_config, smooth_dlnr, 2.5, theta_fine, theta_fine, result);
  // printf("Ntot,%d\n",Ntot);
  // end_t = clock();
  // printf("Total time taken by CPU: %f\n", (double)(end_t - start_t) / CLOCKS_PER_SEC  );

  int ind_theta_first = Ntot - 1 - N_extra-ind_first - Nsub*(Ntheta); // different from nobin case because ind_first is set different!
  for(i=0; i<Ntheta; i++) {
    for(j=0; j<Ntheta; j++) {
      cov_g_interp[i][j] = result[ind_theta_first+i*Nsub][ind_theta_first+j*Nsub];
      cov_g_interp[i][j] /= ( (theta[i+1]*theta[i+1]-theta[i]*theta[i])*(theta[j+1]*theta[j+1]-theta[j]*theta[j]) *M_PI*M_PI*M_PI*survey.area/41253.);
    } 
    // printf("theta[%d]=%lg, cov_g_interp[%d][%d]: %lg\n",i, theta[i], i,i, cov_g_interp[i][i]);
    cov_g_interp[i][i] += noise_factor * N[i];
  }


  // end_t = clock();
  // printf("Total time taken by CPU: %f\n", (double)(end_t - start_t) / CLOCKS_PER_SEC  );

  free(result);free(fk1k2);

}

void cov_NG_real_fft_bin_template(double *theta, int Ntheta, double **cov_ng_interp, config my_config, int *array, double *task_type) {
  int i,j;
  my_config.sym_Flag = 0;
  my_config.c_window_width = 0.25;

  // get input func log-sampled in ell
  const int Nsub = 20; // split each (theta[i],theta[i+1]) range into Nsub sub-ranges, i.e. fill in (Nsub-1) sampling points within the range

  // clock_t start_t, end_t;
  // start_t = clock();

  double ellmax = 50000.;
  double ell_first = 1./theta[Ntheta]; // include the last bin edge
  double ell_last = 1./theta[0];
  double dlnell = log(theta[1]/theta[0])/Nsub;
  int ind_first = floor(log(ell_first)/dlnell); // index of ell_first in log-sampled array (from range [1,ellmax]) containing our orginal desired array (1/theta)
  int Nell = floor(log(ellmax/ell_first)/dlnell) + ind_first + 1; // total sampling points in range [1,ellmax]
  double ell[Nell];
  // printf("ell_first:%lg\n", ell_first);
  // printf("Nell:%d\n", Nell);
  for(i=0; i<Nell; i++) {
    ell[i] = exp(dlnell * (i-ind_first)) * ell_first;
  }

  double (*func_for_cov_NG)(double, double, int*);
  if(task_type[0]==0) {
    if(task_type[1]==0) {func_for_cov_NG = &func_for_cov_NG_cl;}
    else if(task_type[1]==1) {func_for_cov_NG = &func_for_cov_NG_cl_gl;}
    else if(task_type[1]==2) {func_for_cov_NG = &func_for_cov_NG_cl_shear;}
  }
  else if(task_type[0]==1) {
    if(task_type[1]==1) {func_for_cov_NG = &func_for_cov_NG_gl;}
    else if(task_type[1]==2) {func_for_cov_NG = &func_for_cov_NG_gl_shear;}
  }
  else {
    if(task_type[1]==2) {func_for_cov_NG = &func_for_cov_NG_shear;}
  }


  double **func;
  func = malloc(Nell*sizeof(double *));
  for(i=0; i<Nell; i++) {
    func[i] = malloc(Nell * sizeof(double));
    for(j=0; j<Nell; j++) {
      func[i][j] = func_for_cov_NG(ell[i], ell[j], array);
    }
  }

  // FILE * output = fopen("func_ng.txt", "w");
  // for(i=0;i<Nell;i++){
  //   for(j=0;j<Nell;j++){
  //   fprintf(output, "%lg ", func[i][j]);
  // }
  // fprintf(output, "\n");
  // }
  // fclose(output);
  // exit(0);

  // extrapolate
  const int N_extra = 800;
  int Ntot = N_extra*2+Nell;
  double large_ell[Ntot];
  double **fk1k2;
  fk1k2 = malloc(Ntot*sizeof(double *));
  for(i=0;i<Ntot;i++) {
    fk1k2[i] = malloc(Ntot*sizeof(double));
  }  

  extrap_log_linear_cfastcov(ell, Nell, N_extra, large_ell);
  extrap_2dzeros(func, Nell, N_extra, fk1k2);

  // do integration
  double theta_fine[Ntot];
  double **result;
  result = malloc(Ntot*sizeof(double *));
  for(i=0;i<Ntot;i++) {
    result[i] = malloc(Ntot*sizeof(double));
  }
  // twobessel(large_ell, large_ell, fk1k2, Ntot, Ntot, &my_config, theta_fine, theta_fine, result);
  double smooth_dlnr = dlnell*Nsub;
  twobessel_binave(large_ell, large_ell, fk1k2, Ntot, Ntot, &my_config, smooth_dlnr, 2.5, theta_fine, theta_fine, result);
  printf("Ntot,%d\n",Ntot);

  int ind_theta_first = Ntot - 1 - N_extra-ind_first - Nsub*(Ntheta); // different from nobin case because ind_first is set different!
  for(i=0; i<Ntheta; i++) {
    for(j=0; j<Ntheta; j++) {
      cov_ng_interp[i][j] = result[ind_theta_first+i*Nsub][ind_theta_first+j*Nsub];
      cov_ng_interp[i][j] /= ( (theta[i+1]*theta[i+1]-theta[i]*theta[i])*(theta[j+1]*theta[j+1]-theta[j]*theta[j]) /2. *M_PI*M_PI*M_PI);
    }
  }


  free(result);free(fk1k2);free(func);
}




/********** Functions for differnt covariances ***************************************/
double func_for_cov_G_shear(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[1]; n2 = ar[2]; n3 = ar[3]; n4 = ar[4];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_shear_shear_IA_tab(l,n1,n3);
  C24 = C_shear_shear_IA_tab(l,n2,n4);
  C14 = C_shear_shear_IA_tab(l,n1,n4);
  C23 = C_shear_shear_IA_tab(l,n2,n3);
  
  if (n1 == n3){N13= pow(survey.sigma_e,2.0)/(2.*nsource(n1)*survey.n_gal_conversion_factor);}
  if (n1 == n4){N14= pow(survey.sigma_e,2.0)/(2.*nsource(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= pow(survey.sigma_e,2.0)/(2.*nsource(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.*nsource(n2)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*l*l*l;
}

double func_for_cov_NG_shear(double l1, double l2, int *ar){
  double tri = 0.,res =0;
  int n1,n2,n3,n4;
  n1 = ar[1]; n2 = ar[2]; n3 = ar[3]; n4 = ar[4];
  tri= bin_cov_NG_shear_shear_tomo(l1,l2,n1,n2,n3,n4);
  // if(isnan(tri)) {tri=0.;}
  return tri * pow(l1, 2.5) * pow(l2, 2.5);

}


double func_for_cov_G_cl(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[1]; n2 = ar[2]; n3 = ar[3]; n4 = ar[4];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo(l,n1,n3);
  C24 = C_cl_tomo(l,n2,n4);
  C14 = C_cl_tomo(l,n1,n4);
  C23 = C_cl_tomo(l,n2,n3);

  N13= tomo.n_lens_ij[n1][n3]/(nlens(n1)*nlens(n3)*survey.n_gal_conversion_factor);
  N14= tomo.n_lens_ij[n1][n4]/(nlens(n1)*nlens(n4)*survey.n_gal_conversion_factor);
  N23= tomo.n_lens_ij[n2][n3]/(nlens(n2)*nlens(n3)*survey.n_gal_conversion_factor);
  N24= tomo.n_lens_ij[n2][n4]/(nlens(n2)*nlens(n4)*survey.n_gal_conversion_factor);

  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*l*l*l;
}

double func_for_cov_NG_cl(double l1, double l2, int *ar){
  double tri = 0.,res =0;
  int n1,n2,n3,n4;
  n1 = ar[1]; n2 = ar[2]; n3 = ar[3]; n4 = ar[4];
  tri= bin_cov_NG_cl_cl_tomo(l1,l2,n1,n2,n3,n4);
  // if(isnan(tri)) {tri=0.;}
  return tri * pow(l1, 2.5) * pow(l2, 2.5);
}


double func_for_cov_G_cl_gl(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[1]; n2 = ar[2]; n3 = ar[3]; n4 = ar[4];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo(l,n1,n3);
  C24 = C_ggl_IA_tab(l,n2,n4);
  C14 = C_ggl_IA_tab(l,n1,n4);
  C23 = C_cl_tomo(l,n2,n3);
  
  N13= tomo.n_lens_ij[n1][n3]/(nlens(n1)*nlens(n3)*survey.n_gal_conversion_factor);
  N23= tomo.n_lens_ij[n2][n3]/(nlens(n2)*nlens(n3)*survey.n_gal_conversion_factor);
  //printf("%lg, %lg, %lg, %lg, %lg, %lg, %lg\n", l, C13, C24, C14, C23, N13, N23);
  return (C13*C24+C13*N24+C24*N13+C14*C23+C14*N23+C23*N14)*l*l*l;
}

double func_for_cov_NG_cl_gl(double l1, double l2, int *ar){
  double tri = 0.,res =0;
  int n1,n2,n3,n4;
  n1 = ar[1]; n2 = ar[2]; n3 = ar[3]; n4 = ar[4];
  tri= bin_cov_NG_cl_gl_tomo(l1,l2,n1,n2,n3,n4);
  // if(isnan(tri)) {tri=0.;}
  return tri * pow(l1, 2.5) * pow(l2, 2.5);
}


double func_for_cov_G_cl_shear(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[1]; n2 = ar[2]; n3 = ar[3]; n4 = ar[4];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_ggl_IA_tab(l,n1,n3);
  C24 = C_ggl_IA_tab(l,n2,n4);
  C14 = C_ggl_IA_tab(l,n1,n4);
  C23 = C_ggl_IA_tab(l,n2,n3);

  return (C13*C24+ C14*C23)*l*l*l;
}

double func_for_cov_NG_cl_shear(double l1, double l2, int *ar){
  double tri = 0.,res =0;
  int n1,n2,n3,n4;
  n1 = ar[1]; n2 = ar[2]; n3 = ar[3]; n4 = ar[4];
  tri= bin_cov_NG_cl_shear_tomo(l1,l2,n1,n2,n3,n4);
  // if(isnan(tri)) {tri=0.;}
  return tri * pow(l1, 2.5) * pow(l2, 2.5);
}


double func_for_cov_G_gl(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[1]; n2 = ar[2]; n3 = ar[3]; n4 = ar[4];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo(l,n1,n3);
  C24 = C_shear_shear_IA_tab(l,n2,n4);
  C14 = C_ggl_IA_tab(l,n1,n4);
  C23 = C_ggl_IA_tab(l,n3,n2);
  
  N13= tomo.n_lens_ij[n1][n3]/(nlens(n1)*nlens(n3)*survey.n_gal_conversion_factor);
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.0*nsource(n2)*survey.n_gal_conversion_factor);}
  // printf("%lg, %lg, %lg, %lg, %lg, %lg, %lg\n", l, C13, C24, C14, C23, N13, N24);
  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*l*l*l;
}

double func_for_cov_NG_gl(double l1, double l2, int *ar){
  double tri = 0.,res =0;
  int n1,n2,n3,n4;
  n1 = ar[1]; n2 = ar[2]; n3 = ar[3]; n4 = ar[4];
  tri= bin_cov_NG_gl_gl_tomo(l1,l2,n1,n2,n3,n4);
  // if(isnan(tri)) {tri=0.;}
  return tri * pow(l1, 2.5) * pow(l2, 2.5);
}


double func_for_cov_G_gl_shear(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[1]; n2 = ar[2]; n3 = ar[3]; n4 = ar[4];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_ggl_IA_tab(l,n1,n3);
  C24 = C_shear_shear_IA_tab(l,n2,n4);
  C14 = C_ggl_IA_tab(l,n1,n4);
  C23 = C_shear_shear_IA_tab(l,n2,n3);
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.0*nsource(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= pow(survey.sigma_e,2.0)/(2.0*nsource(n2)*survey.n_gal_conversion_factor);}
 
  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*l*l*l;
}

double func_for_cov_NG_gl_shear(double l1, double l2, int *ar){
  double tri = 0.,res =0;
  int n1,n2,n3,n4;
  n1 = ar[1]; n2 = ar[2]; n3 = ar[3]; n4 = ar[4];
  tri= bin_cov_NG_gl_shear_tomo(l1,l2,n1,n2,n3,n4);
  // if(isnan(tri)) {tri=0.;}
  return tri * pow(l1, 2.5) * pow(l2, 2.5);
}
