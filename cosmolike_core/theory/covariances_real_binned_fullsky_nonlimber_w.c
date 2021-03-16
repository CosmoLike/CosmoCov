/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by CosmoLike developers
******************************************************************************/

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
// Note the look-up tables for power spectrum covariances are recomputed if one the redshift bins changes
//      Loop over theta1, theta2 first, before computing the next combination of redshift bins
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
// The vector called "theta" contains theta_min, the lower bin boundaries
// This is used in the bin-average and shot/shape noise computation
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
double cov_G_shear_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2,int z3,int z4,int pm1,int pm2); //Version of Gaussian cov calculation for wide bins
double cov_NG_shear_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2,int z3,int z4,int pm1,int pm2);

double cov_NG_cl_cl_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2, int z3, int z4);
double cov_G_cl_cl_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2, int z3, int z4);

double cov_NG_gl_gl_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2, int z3, int z4);
double cov_G_gl_gl_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1l,int z1s, int z2l, int z2s);

double cov_NG_cl_gl_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2, int zl, int zs);//z1,z2 clustering bins; zl,zs g-g lensing bins
double cov_G_cl_gl_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2, int zl, int zs);

double cov_NG_gl_shear_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int zl,int zs, int z3, int z4, int pm);
double cov_G_gl_shear_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max,int zl,int zs, int z3, int z4,int pm);

double cov_NG_cl_shear_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2, int z3, int z4, int pm); //z1,z2 clustering bins; z3,z4 shear bins
double cov_G_cl_shear_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2, int z3, int z4, int pm);

double func_for_cov_G_shear(double l, int *ar);
double func_for_cov_G_cl(double l, int *ar);
double func_for_cov_G_gl(double l, int *ar);
double func_for_cov_G_cl_shear(double l, int *ar);
double func_for_cov_G_cl_gl(double l, int *ar);
double func_for_cov_G_gl_shear(double l, int *ar);

/////// full sky covs

void cov_shear_shear_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4,int pm1,int pm2, int FLAG_NG, double *theta, double *dtheta){

  int i,j;
  static int LMAX = 50000;
  static double **Glplus =0;
  static double **Glminus =0;
  if (Glplus ==0){
    Glplus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    Glminus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    double *xmin, *xmax, *Pmin, *Pmax, *dPmin, *dPmax;
    xmin= create_double_vector(0, like.Ntheta-1);
    xmax= create_double_vector(0, like.Ntheta-1);
    for(i=0; i<like.Ntheta ; i++){
      xmin[i]=cos(theta[i]);
      xmax[i]=cos(theta[i+1]);
    }
    Pmin= create_double_vector(0, LMAX+1);
    Pmax= create_double_vector(0, LMAX+1);
    dPmin= create_double_vector(0, LMAX+1);
    dPmax= create_double_vector(0, LMAX+1);
    for (i = 0; i<like.Ntheta; i ++){
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmin[i],Pmin,dPmin);
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmax[i],Pmax,dPmax);
      // printf("xmax[%d]:%lg\n", i,xmax[i]);
      for (int l = 2; l < LMAX; l ++){
        /*double plm = gsl_sf_legendre_Plm(l,2,x);
        double plm_1 = gsl_sf_legendre_Plm(l-1,2,x);
        Glplus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
        *(plm*((4-l+2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
        +plm_1*(l-1,2,x)*(l+2)*(x-2)/(1-x*x));


        Glminus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
        *(plm*(l,2,x)*((4-l-2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
        +plm_1*(l-1,2,x)*(l+2)*(x+2)/(1-x*x));*/

        Glplus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1.)/(2.*l+1)         * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        +2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        -2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);

        Glminus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1.)/(2.*l+1)         * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        -2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        +2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);

        // printf("Pmin[%d], Pmax[%d]:%lg, %lg\n", l,l, Pmin[l], Pmax[l]);
        // printf("Glplus[%d][%d],%lg\n", i,l, Glplus[i][l]);
      }
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }

  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ntheta ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  // printf("i,j:%d,%d\n", i,j);
  double triGl1;
  double covGl1, cov_g_l;
  int ar[4];
  ar[0] = z1; ar[1] = z2; ar[2] = z3; ar[3] = z4;

  double N[like.Ntheta];
  for(i=0;i<like.Ntheta;i++) {N[i] = 0.;}
  if (z1 ==z3 && z2 ==z4 && pm1 == pm2){ //&& pm1 == pm2 required as C+- doesn't have the diagonal shot noise term
    for(i=0;i<like.Ntheta;i++) {
      N[i] += pow(survey.sigma_e,4.0)/(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    } //(number of galaxy pairs in the survey contributing to annulus of width Dtheta centered at theta1)^-1
  }
  if (z1 ==z4 && z2 ==z3 && pm1 == pm2){ //&& pm1 == pm2 required as C+- doesn't have the diagonal shot noise term
    for(i=0;i<like.Ntheta;i++) {
      N[i] += pow(survey.sigma_e,4.0)/(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    }
  }
  for(i=0;i<like.Ntheta;i++) {if(N[i]) N[i] /= w_mask(theta[i]); printf("N[i]: %lg\n", N[i]);}

  if(pm1>0 && pm2>0){
    for (l1 = 2; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Glplus[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glplus[j][l1];
        }
      }
      // printf("l1,%d\n", l1);
      if(FLAG_NG){
        for (l2 = 2; l2 < LMAX; l2++){
          tri = bin_cov_NG_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          for(i=0; i<like.Ntheta ; i++){
            triGl1 = tri * Glplus[i][l1];
            // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triGl1 * Glplus[j][l2];
            }
          }
        }
      }
    }
  }
  else if(pm1>0 && pm2==0){
    for (l1 = 2; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Glplus[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glminus[j][l1];
        }
      }
      if(FLAG_NG){
        for (l2 = 2; l2 < LMAX; l2++){
          tri = bin_cov_NG_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          for(i=0; i<like.Ntheta ; i++){
            triGl1 = tri * Glplus[i][l1];
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triGl1 * Glminus[j][l2];
            }
          }
        }
      }
    }
  }
  else if(pm1==0 && pm2>0){
    for (l1 = 2; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Glminus[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glplus[j][l1];
        }
      }
      if(FLAG_NG){
        for (l2 = 2; l2 < LMAX; l2++){
          tri = bin_cov_NG_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          for(i=0; i<like.Ntheta ; i++){
            triGl1 = tri * Glminus[i][l1];
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triGl1 * Glplus[j][l2];
            }
          }
        }
      }
    }
  }
  else{
    for (l1 = 2; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Glminus[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glminus[j][l1];
        }
      }
      if(FLAG_NG){
        for (l2 = 2; l2 < LMAX; l2++){
          tri = bin_cov_NG_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          for(i=0; i<like.Ntheta ; i++){
            triGl1 = tri * Glminus[i][l1];
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triGl1 * Glminus[j][l2];
            }
          }
        }
      }
    }
  }
  for(i=0; i<like.Ntheta; i++){cov[i][i] += 2.* N[i];}
}

void cov_gl_shear_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta){

  int i,j;
  static int LMAX = 50000;
  static double **Glplus =0;
  static double **Glminus =0;
  static double **Pl =0;
  if (Glplus ==0){
    Glplus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    Glminus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);

    double *xmin, *xmax, *Pmin, *Pmax, *dPmin, *dPmax;
    xmin= create_double_vector(0, like.Ntheta-1);
    xmax= create_double_vector(0, like.Ntheta-1);
    for(i=0; i<like.Ntheta ; i++){
      xmin[i]=cos(theta[i]);
      xmax[i]=cos(theta[i+1]);
    }
    Pmin= create_double_vector(0, LMAX+1);
    Pmax= create_double_vector(0, LMAX+1);
    dPmin= create_double_vector(0, LMAX+1);
    dPmax= create_double_vector(0, LMAX+1);
    for (i = 0; i<like.Ntheta; i ++){
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmin[i],Pmin,dPmin);
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmax[i],Pmax,dPmax);
      // printf("xmax[%d]:%lg\n", i,xmax[i]);
      for (int l = 1; l < LMAX; l ++){
        /*double plm = gsl_sf_legendre_Plm(l,2,x);
        double plm_1 = gsl_sf_legendre_Plm(l-1,2,x);
        Glplus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
        *(plm*((4-l+2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
        +plm_1*(l-1,2,x)*(l+2)*(x-2)/(1-x*x));


        Glminus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
        *(plm*(l,2,x)*((4-l-2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
        +plm_1*(l-1,2,x)*(l+2)*(x+2)/(1-x*x));*/

        Glplus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1.)/(2.*l+1)         * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        +2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        -2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);

        Glminus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1.)/(2.*l+1)         * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        -2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        +2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);

        // printf("Pmin[%d], Pmax[%d]:%lg, %lg\n", l,l, Pmin[l], Pmax[l]);
        // printf("Glplus[%d][%d],%lg\n", i,l, Glplus[i][l]);
        Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
        *((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
        +(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        -2./(2*l+1.)*(Pmin[l+1]-Pmax[l+1]));
      }
      Glplus[i][1] = 0.; Glminus[i][1] = 0.;
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }

  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ntheta ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  // printf("i,j:%d,%d\n", i,j);
  double triP;
  double covGl1, cov_g_l;
  int ar[4];
  ar[0] = z1; ar[1] = z2; ar[2] = z3; ar[3] = z4;

  if(pm>0){
    for (l1 = 1; l1 < LMAX; l1++){
      l1_double = (double)l1;
      // printf("l1,%d\n", l1);
      cov_g_l = func_for_cov_G_gl_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Pl[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glplus[j][l1];
        }
      }

      if(FLAG_NG){
        for (l2 = 2; l2 < LMAX; l2++){
          tri = bin_cov_NG_gl_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          for(i=0; i<like.Ntheta ; i++){
            triP = tri * Pl[i][l1];
            // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triP * Glplus[j][l2];
            }
          }
        }
      }
    }
  }
  else{
    for (l1 = 1; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_gl_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Pl[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glminus[j][l1];
        }
      }

      if(FLAG_NG){
        for (l2 = 2; l2 < LMAX; l2++){
          tri = bin_cov_NG_gl_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          for(i=0; i<like.Ntheta ; i++){
            triP = tri * Pl[i][l1];
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triP * Glminus[j][l2];
            }
          }
        }
      }
    }
  }
}

void cov_cl_shear_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta){

  int i,j;
  static int LMAX = 10000;
  static double **Glplus =0;
  static double **Glminus =0;
  static double **Pl =0;
  if (Glplus ==0){
    Glplus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    Glminus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);

    double *xmin, *xmax, *Pmin, *Pmax, *dPmin, *dPmax;
    xmin= create_double_vector(0, like.Ntheta-1);
    xmax= create_double_vector(0, like.Ntheta-1);
    for(i=0; i<like.Ntheta ; i++){
      xmin[i]=cos(theta[i]);
      xmax[i]=cos(theta[i+1]);
    }
    Pmin= create_double_vector(0, LMAX+1);
    Pmax= create_double_vector(0, LMAX+1);
    dPmin= create_double_vector(0, LMAX+1);
    dPmax= create_double_vector(0, LMAX+1);
    for (i = 0; i<like.Ntheta; i ++){
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmin[i],Pmin,dPmin);
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmax[i],Pmax,dPmax);
      // printf("xmax[%d]:%lg\n", i,xmax[i]);
      for (int l = 1; l < LMAX; l ++){
        /*double plm = gsl_sf_legendre_Plm(l,2,x);
        double plm_1 = gsl_sf_legendre_Plm(l-1,2,x);
        Glplus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
        *(plm*((4-l+2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
        +plm_1*(l-1,2,x)*(l+2)*(x-2)/(1-x*x));


        Glminus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
        *(plm*(l,2,x)*((4-l-2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
        +plm_1*(l-1,2,x)*(l+2)*(x+2)/(1-x*x));*/

        Glplus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1.)/(2.*l+1)         * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        +2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        -2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);

        Glminus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1.)/(2.*l+1)         * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        -2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        +2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);

        Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/(xmin[i]-xmax[i]);
        // printf("Pmin[%d], Pmax[%d]:%lg, %lg\n", l,l, Pmin[l], Pmax[l]);
        // printf("Glplus[%d][%d],%lg\n", i,l, Glplus[i][l]);
      }
      Pl[i][0] = 1./(4.*M_PI);
      Glplus[i][1] = 0.; Glplus[i][0] = 0.;
      Glminus[i][1] = 0.; Glminus[i][0] = 0.;
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }

  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ntheta ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  // printf("i,j:%d,%d\n", i,j);
  double triP;
  double covGl1, cov_g_l;
  int ar[4];
  ar[0] = z1; ar[1] = z2; ar[2] = z3; ar[3] = z4;

  if(pm>0){
    for (l1 = 1; l1 < LMAX; l1++){
      l1_double = (double)l1;
      // printf("l1,%d\n", l1);
      cov_g_l = func_for_cov_G_cl_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Pl[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glplus[j][l1];
        }
      }

      if(FLAG_NG){
        for (l2 = 2; l2 < LMAX; l2++){
          tri = bin_cov_NG_cl_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          for(i=0; i<like.Ntheta ; i++){
            triP = tri * Pl[i][l1];
            // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triP * Glplus[j][l2];
            }
          }
        }
      }
    }
  }
  else{
    for (l1 = 1; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_cl_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Pl[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glminus[j][l1];
        }
      }

      if(FLAG_NG){
        for (l2 = 2; l2 < LMAX; l2++){
          tri = bin_cov_NG_cl_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          // printf("tri:%lg\n", tri);
          for(i=0; i<like.Ntheta ; i++){
            triP = tri * Pl[i][l1];
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triP * Glminus[j][l2];
            }
          }
        }
      }
    }
  }
}

void cov_cl_gl_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *theta, double *dtheta){

  int i,j;
  static int LMAX = 50000;
  static double **Pl =0;
  static double **Pl2=0;
  if (Pl ==0){
    Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    Pl2=create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);

    double *xmin, *xmax, *Pmin, *Pmax, *dPmin, *dPmax;
    xmin= create_double_vector(0, like.Ntheta-1);
    xmax= create_double_vector(0, like.Ntheta-1);
    for(i=0; i<like.Ntheta ; i++){
      xmin[i]=cos(theta[i]);
      xmax[i]=cos(theta[i+1]);
    }
    Pmin= create_double_vector(0, LMAX+1);
    Pmax= create_double_vector(0, LMAX+1);
    dPmin= create_double_vector(0, LMAX+1);
    dPmax= create_double_vector(0, LMAX+1);
    for (i = 0; i<like.Ntheta; i ++){
      gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
      Pl[i][0] = 1./(4.*M_PI);
      for (int l = 1; l < LMAX; l ++){
        Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/(xmin[i]-xmax[i]);
        Pl2[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
        *((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
        +(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        -2./(2*l+1.)*(Pmin[l+1]-Pmax[l+1]));
      }
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }

  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ntheta ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  double triP;
  double covGl1, cov_g_l;
  int ar[4];
  ar[0] = z1; ar[1] = z2; ar[2] = z3; ar[3] = z4;

  for (l1 = 1; l1 < LMAX; l1++){
    l1_double = (double)l1;
    cov_g_l = func_for_cov_G_cl_gl(l1_double, ar);
    for(i=0; i<like.Ntheta ; i++){
      covGl1 = cov_g_l * Pl[i][l1];
      // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
      for(j=0; j<like.Ntheta ; j++){
        cov[i][j] += covGl1 * Pl2[j][l1];
      }
    }

    if(FLAG_NG){
      for (l2 = 1; l2 < LMAX; l2++){
        tri = bin_cov_NG_cl_gl_tomo(l1_double,(double)l2,z1,z2,z3,z4);
        for(i=0; i<like.Ntheta ; i++){
          triP = tri * Pl[i][l1];
          for(j=0; j<like.Ntheta ; j++){
            covNG[i][j] += triP * Pl2[j][l2];
          }
        }
      }
    }
  }
}

void cov_gl_gl_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *theta, double *dtheta){

  int i,j;
  static int LMAX = 50000;
  static double **Pl =0;
  if (Pl ==0){
    Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);

    double *xmin, *xmax, *Pmin, *Pmax, *dPmin, *dPmax;
    xmin= create_double_vector(0, like.Ntheta-1);
    xmax= create_double_vector(0, like.Ntheta-1);
    for(i=0; i<like.Ntheta ; i++){
      xmin[i]=cos(theta[i]);
      xmax[i]=cos(theta[i+1]);
    }
    Pmin= create_double_vector(0, LMAX+1);
    Pmax= create_double_vector(0, LMAX+1);
    dPmin= create_double_vector(0, LMAX+1);
    dPmax= create_double_vector(0, LMAX+1);
    for (i = 0; i<like.Ntheta; i ++){
      gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
      for (int l = 1; l < LMAX; l ++){
        Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
        *((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
        +(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        -2./(2*l+1.)*(Pmin[l+1]-Pmax[l+1]));
      }
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }

  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ntheta ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  double N[like.Ntheta];
  for(i=0;i<like.Ntheta;i++) {N[i] = 0.;}

  int zl1, zs1, zl2, zs2;
  zl1 = z1; zs1 = z2; zl2 = z3; zs2 = z4;
  if (zl1 ==zl2 && zs1 ==zs2){
    for(i=0;i<like.Ntheta;i++) {
      N[i] += pow(survey.sigma_e,2.0)/(2.0*M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*nlens(zl1)*nsource(zs2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    }
  }
  for(i=0;i<like.Ntheta;i++) {if(N[i]) N[i] /= w_mask(theta[i]);}

  double triP;
  double covGl1, cov_g_l;
  int ar[4];
  ar[0] = z1; ar[1] = z2; ar[2] = z3; ar[3] = z4;

  for (l1 = 1; l1 < LMAX; l1++){
    l1_double = (double)l1;
    cov_g_l = func_for_cov_G_gl(l1_double, ar);
    for(i=0; i<like.Ntheta ; i++){
      covGl1 = cov_g_l * Pl[i][l1];
      // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
      for(j=0; j<like.Ntheta ; j++){
        cov[i][j] += covGl1 * Pl[j][l1];
      }
    }

    if(FLAG_NG){
      for (l2 = 1; l2 < LMAX; l2++){
        tri = bin_cov_NG_gl_gl_tomo(l1_double,(double)l2,z1,z2,z3,z4);
        for(i=0; i<like.Ntheta ; i++){
          triP = tri * Pl[i][l1];
          for(j=0; j<like.Ntheta ; j++){
            covNG[i][j] += triP * Pl[j][l2];
          }
        }
      }
    }
  }
  for(i=0; i<like.Ntheta; i++){cov[i][i] += N[i];}
}

void cov_cl_cl_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *theta, double *dtheta){

  int i,j;
  static int LMAX = 50000;
  static double **Pl =0;
  if (Pl ==0){
    Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);

    double *xmin, *xmax, *Pmin, *Pmax, *dPmin, *dPmax;
    xmin= create_double_vector(0, like.Ntheta-1);
    xmax= create_double_vector(0, like.Ntheta-1);
    for(i=0; i<like.Ntheta ; i++){
      xmin[i]=cos(theta[i]);
      xmax[i]=cos(theta[i+1]);
    }
    Pmin= create_double_vector(0, LMAX+1);
    Pmax= create_double_vector(0, LMAX+1);
    dPmin= create_double_vector(0, LMAX+1);
    dPmax= create_double_vector(0, LMAX+1);
    for (i = 0; i<like.Ntheta; i ++){
      gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
      Pl[i][0] = 1./(4.*M_PI);
      for (int l = 1; l < LMAX; l ++){
        Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/(xmin[i]-xmax[i]);
      }
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }

  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ntheta ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  double N[like.Ntheta];
  for(i=0;i<like.Ntheta;i++) {N[i] = 0.;}

  if (z1 ==z3 && z2 ==z4){
    for(i=0;i<like.Ntheta;i++) {
      N[i] += 1./(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*nlens(z1)*nlens(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor); //(number of galaxy pairs in the survey contributing to annulus of width Dtheta centered at theta1)^-1
    }
  }
  if (z1 ==z4 && z2 ==z3){
    for(i=0;i<like.Ntheta;i++) {
      N[i] += 1./(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*nlens(z1)*nlens(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    }
  }
  for(i=0;i<like.Ntheta;i++) {if(N[i]) N[i] /= w_mask(theta[i]);}

  double triP;
  double covGl1, cov_g_l;
  int ar[4];
  ar[0] = z1; ar[1] = z2; ar[2] = z3; ar[3] = z4;

  for (l1 = 0; l1 < LMAX; l1++){
    l1_double = (double)l1;
    cov_g_l = func_for_cov_G_cl(l1_double, ar);
    for(i=0; i<like.Ntheta ; i++){
      covGl1 = cov_g_l * Pl[i][l1];
      // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
      for(j=0; j<like.Ntheta ; j++){
        cov[i][j] += covGl1 * Pl[j][l1];
      }
    }

    if(FLAG_NG && (l1>0)){
      for (l2 = 1; l2 < LMAX; l2++){
        tri = bin_cov_NG_cl_cl_tomo(l1_double,(double)l2,z1,z2,z3,z4);
        for(i=0; i<like.Ntheta ; i++){
          triP = tri * Pl[i][l1];
          for(j=0; j<like.Ntheta ; j++){
            covNG[i][j] += triP * Pl[j][l2];
          }
        }
      }
    }
  }
  for(i=0; i<like.Ntheta; i++){cov[i][i] += N[i];}
}



/********** Functions for differnt covariances ***************************************/
double func_for_cov_G_shear(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  // printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_shear_shear_IA_tab(l,n1,n3);
  C24 = C_shear_shear_IA_tab(l,n2,n4);
  C14 = C_shear_shear_IA_tab(l,n1,n4);
  C23 = C_shear_shear_IA_tab(l,n2,n3);

  if (n1 == n3){N13= pow(survey.sigma_e,2.0)/(2.*nsource(n1)*survey.n_gal_conversion_factor);}
  if (n1 == n4){N14= pow(survey.sigma_e,2.0)/(2.*nsource(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= pow(survey.sigma_e,2.0)/(2.*nsource(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.*nsource(n2)*survey.n_gal_conversion_factor);}

  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}

double func_for_cov_G_cl(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo_nonlimber_interp(l,n1,n3);
  C24 = C_cl_tomo_nonlimber_interp(l,n2,n4);
  C14 = C_cl_tomo_nonlimber_interp(l,n1,n4);
  C23 = C_cl_tomo_nonlimber_interp(l,n2,n3);

  N13= tomo.n_lens_ij[n1][n3]/(nlens(n1)*nlens(n3)*survey.n_gal_conversion_factor);
  N14= tomo.n_lens_ij[n1][n4]/(nlens(n1)*nlens(n4)*survey.n_gal_conversion_factor);
  N23= tomo.n_lens_ij[n2][n3]/(nlens(n2)*nlens(n3)*survey.n_gal_conversion_factor);
  N24= tomo.n_lens_ij[n2][n4]/(nlens(n2)*nlens(n4)*survey.n_gal_conversion_factor);

  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}


double func_for_cov_G_cl_gl(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo_nonlimber_interp(l,n1,n3);
  C24 = C_ggl_IA_tab(l,n2,n4);
  C14 = C_ggl_IA_tab(l,n1,n4);
  C23 = C_cl_tomo_nonlimber_interp(l,n2,n3);

  N13= tomo.n_lens_ij[n1][n3]/(nlens(n1)*nlens(n3)*survey.n_gal_conversion_factor);
  N23= tomo.n_lens_ij[n2][n3]/(nlens(n2)*nlens(n3)*survey.n_gal_conversion_factor);

  return (C13*C24+C13*N24+C24*N13+C14*C23+C14*N23+C23*N14)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}


double func_for_cov_G_cl_shear(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_ggl_IA_tab(l,n1,n3);
  C24 = C_ggl_IA_tab(l,n2,n4);
  C14 = C_ggl_IA_tab(l,n1,n4);
  C23 = C_ggl_IA_tab(l,n2,n3);

  return (C13*C24+ C14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}

double func_for_cov_G_gl(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo_nonlimber_interp(l,n1,n3);
  C24 = C_shear_shear_IA_tab(l,n2,n4);
  C14 = C_ggl_IA_tab(l,n1,n4);
  C23 = C_ggl_IA_tab(l,n3,n2);

  N13= tomo.n_lens_ij[n1][n3]/(nlens(n1)*nlens(n3)*survey.n_gal_conversion_factor);
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.0*nsource(n2)*survey.n_gal_conversion_factor);}

  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}

double func_for_cov_G_gl_shear(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_ggl_IA_tab(l,n1,n3);
  C24 = C_shear_shear_IA_tab(l,n2,n4);
  C14 = C_ggl_IA_tab(l,n1,n4);
  C23 = C_shear_shear_IA_tab(l,n2,n3);
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.0*nsource(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= pow(survey.sigma_e,2.0)/(2.0*nsource(n2)*survey.n_gal_conversion_factor);}

  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}
