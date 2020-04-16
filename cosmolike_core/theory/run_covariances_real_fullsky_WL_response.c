/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by CosmoLike developers
******************************************************************************/
void cov_shear_shear_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4,int pm1,int pm2, int FLAG_NG, double *theta, double *dtheta);
double cov_G_shear_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2,int z3,int z4,int pm1,int pm2); //Version of Gaussian cov calculation for wide bins
double cov_NG_shear_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2,int z3,int z4,int pm1,int pm2);

double func_for_cov_G_shear(double l, int *ar);
void run_cov_shear_shear_real_bin(char *OUTFILE, char *PATH, double *theta, double *dtheta, int Ntheta, int n1, int n2, int pm1, int pm2, int start);


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
      for (int l = 3; l < LMAX; l ++){

        Glplus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1)/(2*l+1)           * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        +2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        -2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);

        Glminus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1)/(2*l+1)           * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        -2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        +2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);

      }
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }
  printf("Pl tabulated\n");
  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ntheta ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

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
  for(i=0;i<like.Ntheta;i++) {if(N[i]) N[i] /= w_mask(theta[i]);}

  if(pm1>0 && pm2>0){
    for (l1 = 3; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Glplus[i][l1];
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glplus[j][l1];
        }
      }
      if(FLAG_NG){
        for (l2 = 3; l2 < LMAX; l2++){
          tri = bin_cov_SSC_response_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          if(covparams.cng) tri += bin_cov_NG_response_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          for(i=0; i<like.Ntheta ; i++){
            triGl1 = tri * Glplus[i][l1];
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triGl1 * Glplus[j][l2];
            }
          }
        }
      }
    }
  }
  else if(pm1>0 && pm2==0){
    for (l1 = 3; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Glplus[i][l1];
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glminus[j][l1];
        }
      }
      if(FLAG_NG){
        for (l2 = 3; l2 < LMAX; l2++){
          tri = bin_cov_SSC_response_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          if(covparams.cng) tri += bin_cov_NG_response_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
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
    for (l1 = 3; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Glminus[i][l1];
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glplus[j][l1];
        }
      }
      if(FLAG_NG){
        for (l2 = 3; l2 < LMAX; l2++){
          tri = bin_cov_SSC_response_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4)\
                +bin_cov_NG_response_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
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
    for (l1 = 3; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Glminus[i][l1];
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glminus[j][l1];
        }
      }
      if(FLAG_NG){
        for (l2 = 3; l2 < LMAX; l2++){
          tri = bin_cov_SSC_response_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          if(covparams.cng) tri += bin_cov_NG_response_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
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


/********** Functions for differnt covariances ***************************************/
double func_for_cov_G_shear(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
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


void print_citations(FILE *F1){
  fprintf(F1, "# Please cite the following papers in research using this covariance:\n");
  fprintf(F1, "# arXiv: 1503.03487, https://arxiv.org/abs/1503.03487\n");
  fprintf(F1, "# arXiv: 1601.05779, https://arxiv.org/abs/1601.05779\n");
  fprintf(F1, "# arXiv: 1703.09212, https://arxiv.org/abs/1703.09212\n");
  fprintf(F1, "# arXiv: 1705.01092, https://arxiv.org/abs/1705.01092\n");
  fprintf(F1, "# arXiv: 1711.07467, https://arxiv.org/abs/1711.07467\n");
  fprintf(F1, "# arXiv: 1803.03274, https://arxiv.org/abs/1803.03274\n");
  fprintf(F1, "# arXiv: 2004.04833, https://arxiv.org/abs/2004.04833\n");
  if (w_mask(like.vtmin) < 1.0){
    fprintf(F1, "# arXiv: 1804.10663, https://arxiv.org/abs/1804.10663\n");
  }
  fprintf(F1, "########################################################\n");
}

void print_ordering(FILE *F2, int i, double t, char *probe, char *sample_type, int z1, char *sample_type2, int z2){
    fprintf(F2, "%d %e %s %s%d %s%d\n", i, t, probe, sample_type, z1, sample_type2, z2);
}

void run_cov_shear_shear_real_bin(char *OUTFILE, char *PATH,double *theta, double *dtheta, int Ntheta, int n1, int n2, int pm1, int pm2, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  double c_ng, c_g,sn;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  z1 = Z1(n1); z2 = Z2(n1);
  printf("N_shear1 = %d (%d,%d)\n", n1,z1,z2);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear2 = %d (%d, %d)\n",n2,z3,z4);
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");


  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_shear_shear_real_binned_fullsky(cov_fullsky_G,cov_fullsky_NG,z1,z2,z3,z4,pm1,pm2, covparams.ng, theta, dtheta);

  print_citations(F1);

  FILE *F2;
  char filename2[300];
  int if_write = 0, i1;
  if(pm1==0) {i1 = Ntheta*(tomo.shear_Npowerspectra+n1);}
  else {i1 = Ntheta*n1;}
  sprintf(filename2,"%sorder_%s_i_%d-%d",PATH,covparams.filename,i1,i1+Ntheta-1);
  if (fopen(filename2, "r") == NULL){
    if_write=1;
    F2 =fopen(filename2,"w");
    fprintf(F2, "# i, bin-averaged theta_i, 2pt func, {s: source}tomo bin index 1, {s, l}tomo bin index 2\n");
  }

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
  	double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
	  double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;sn = 0.;
      if(covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      if(pm1==1 && pm2==1) fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*n1+nl1,Ntheta*(n2)+nl2, t1,t2,z1,z2,z3,z4,c_g,c_ng);
      if(pm1==0 && pm2==0) fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(tomo.shear_Npowerspectra+n2)+nl2,t1,t2,z1,z2,z3,z4,c_g,c_ng);
      if(pm1==1 && pm2==0) fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*n1+nl1,Ntheta*(tomo.shear_Npowerspectra+n2)+nl2,t1,t2,z1,z2,z3,z4,c_g,c_ng);
    }
    if(if_write) {
      if(pm1==0){
        fprintf(F2, "%d %e xi- s%d s%d\n", i1+nl1, t1, z1, z2);
      }
      else{
        fprintf(F2, "%d %e xi+ s%d s%d\n", i1+nl1, t1, z1, z2);
      }
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ntheta-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ntheta-1, 0, like.Ntheta-1);
  fclose(F1);
  if(if_write) {fclose(F2);}
}
