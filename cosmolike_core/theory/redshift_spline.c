/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by CosmoLike developers
******************************************************************************/

//#define Z_SPLINE_TYPE gsl_interp_akima
#define Z_SPLINE_TYPE gsl_interp_cspline
// lens efficiencies
double g_tomo(double a, int zbin); // lens efficiency of source galaxies in tomography bin zbin
double g_lens(double a, int zbin); // lens efficiency of *lens* galaxies in tomography bin zbin - used for magnification calculations

//shear routines
double zdistr_photoz(double zz,int j); //returns n(ztrue | j), works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j

//clustering routines
double pf_photoz(double zz,int j); //returns n(ztrue | j), works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j

//lens + source galaxy counts
double int_nsource(double z, void* param);
double int_nlens(double z, void* param);
double nsource(int j); //returns n_gal for shear tomography bin j, works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
double nlens(int j); //returns n_gal for clustering tomography bin j, works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
double zmean(int j);//mean true redshift of (clustering/lens) galaxies within redshift bin j
double zmean_source(int j); //mean true redshift of source galaxies in tomography bin j


/***** redshift integration boundaries **********/
double amin_source(int i);
double amax_source(int i);
double amax_source_IA(int i);
double amin_lens(int i);
double amax_lens(int i);

//routines for association of a pair redshift bin numbers and power spectrum tomography bin
int test_kmax(double l, int zl); //valid l + lens/clustering redshift combination?
int test_zoverlap (int zl, int zs); //valid lens + source redshift bin combination?
int N_ggl(int zl, int zs);//(z_l,z_s) -> N_tomo_ggl
int ZL(int Nbin);//N_tomo_ggl -> z_l
int ZS(int Nbin);//N_tomo_ggl -> z_s
int N_shear (int z1, int z2);//(z_1,z_2)-> N_tomo_shear, assuming z1<=z2
int Z1(int Nbin);//N_tomo_shear -> z1, assuming z1<=z2
int Z2(int Nbin);//N_tomo_shear -> z2, assuming z1<=z2
int Zcl1(int Nbin);// find zcl1 of tomography combination (zcl1,zcl2) constituting galaxy clustering tomography bin Nbin
int Zcl2(int Nbin);// find zcl2 of tomography combination (zcl1,zcl2) constituting galaxy clustering tomography bin Nbin


/********** integration boundary routines *************/
double amin_source(int i){
  return 1./(1+tomo.shear_zmax[i]+2.*fabs(nuisance.bias_zphot_shear[i]));
}
double amax_source(int i){
return 1./(1.+fmax(redshift.shear_zdistrpar_zmin,0.001));
}

double amax_source_IA(int i){
  return 1./(1.+fmax(redshift.shear_zdistrpar_zmin,0.001));
}
double amin_lens(int i){
  return 1./(1+tomo.clustering_zmax[i]+2.*fabs(nuisance.bias_zphot_clustering[i]));
}
double amax_lens(int i){
  if (gbias.b_mag[i] != 0){return 1./(1.+fmax(redshift.shear_zdistrpar_zmin,0.001));}
  return 1./(1+fmax(tomo.clustering_zmin[i]-2.*fabs(nuisance.bias_zphot_clustering[i]),0.001));
}

/************ redshift overlap tests, allowed tomography combinations **********/
int test_kmax(double l, int zl){ //test whether the (l,zl) bin is in the linear clustering regime - return 1 if true, 0 otherwise
  static double chiref[20] = {-1.};
  if (chiref[0] < 0){
    int i;
    for (i = 0; i < tomo.clustering_Nbin; i++){
      chiref[i] = chi(1./(1.+0.5*(tomo.clustering_zmin[i]+tomo.clustering_zmax[i])));
    }
  }
  double R_min = like.Rmin_bias; //set minimum scale to which we trust our bias model, in Mpc/h
  double kmax = constants.twopi/R_min*cosmology.coverH0;
  if((l+0.5)/chiref[zl] < kmax){ return 1;}
  return 0;
}

int test_zoverlap(int zl, int zs){ //test whether source bin zs is behind lens bin zl
  if (tomo.shear_zmax[zs] >= tomo.clustering_zmin[zl]){return 1;}
  return 0;
}

int N_ggl(int zl, int zs){
  static int N[20][20] = {-42};
  if (N[0][0] < 0){
    int i, j,n;
    n = 0;
    for (i = 0; i < tomo.clustering_Nbin; i ++){
      for (j = 0; j < tomo.shear_Nbin; j++){
        if (test_zoverlap(i,j)){N[i][j] = n; n++;}
        else{N[i][j] = -1;}
      }
    }
  }
  return N[zl][zs];
}

int ZL(int Nbin){
  static int N[400] = {-42};
  if (N[0] < -1){
    int i,j,n = 0;
    for (i = 0; i < tomo.clustering_Nbin; i ++){
      for (j = 0; j < tomo.shear_Nbin; j++){
        if (test_zoverlap(i,j)){N[n] = i; n++;}
      }
    }
  }
  return N[Nbin];
}
int ZS(int Nbin){
  static int N[400] = {-42};
  if (N[0] < -1){
    int i,j,n = 0;
    for (i = 0; i < tomo.clustering_Nbin; i ++){
      for (j = 0; j < tomo.shear_Nbin; j++){
        if (test_zoverlap(i,j)){N[n] = j; n++;}
      }
    }
  }
  return N[Nbin];
}

int N_shear (int z1, int z2){ //find shear tomography bin number N_shear of tomography combination (z1,z2)
  static int N[20][20] = {-42};
  if (N[0][0] < -1){
    int i, j,n = 0;
    for (i = 0; i < tomo.shear_Nbin; i ++){
      for (j = i; j < tomo.shear_Nbin; j++){
        N[i][j] = n;
        N[j][i] = n;
        n++;
      }
    }
  }
  return N[z1][z2];
}
int Z1(int Nbin){// find z1 of tomography combination (z1,z2) constituting shear tomography bin Nbin
  static int N[210] = {-42};
  if (N[0] < -1){
    int i, j,n = 0;
    for (i = 0; i < tomo.shear_Nbin; i ++){
      for (j = i; j < tomo.shear_Nbin; j++){
        N[n] = i;
        n++;
      }
    }
  }
  return N[Nbin];
}

int Z2(int Nbin){ // find z2 of tomography combination (z1,z2) constituting shear tomography bin Nbin
  static int N[210]={-42};
  if (N[0] < -1){
    int i, j,n = 0;
    for (i = 0; i < tomo.shear_Nbin; i ++){
      for (j = i; j < tomo.shear_Nbin; j++){
        N[n] = j;
        n++;
      }
    }
  }
  return N[Nbin];
}


int Zcl1(int Nbin){// find zcl1 of tomography combination (zcl1,zcl2) constituting galaxy clustering tomography bin Nbin
  static int N[210] = {-42};
  if (N[0] < -1){
    int i, j,n = 0;
    for (i = 0; i < tomo.clustering_Nbin; i ++){
      for (j = i; j < tomo.clustering_Nbin; j++){
        N[n] = i;
        n++;
      }
    }
  }
  return N[Nbin];
}

int Zcl2(int Nbin){ // find zcl2 of tomography combination (zcl1,zcl2) constituting galaxy clustering tomography bin Nbin
  static int N[210]={-42};
  if (N[0] < -1){
    int i, j,n = 0;
    for (i = 0; i < tomo.clustering_Nbin; i ++){
      for (j = i; j < tomo.clustering_Nbin; j++){
        N[n] = j;
        n++;
      }
    }
  }
  return N[Nbin];
}

int N_clustering_tomo (int z1, int z2){ //find shear tomography bin number N_shear of tomography combination (z1,z2)
  static int N[20][20] = {-42};
  if (N[0][0] < -1){
    int i, j,n = 0;
    for (i = 0; i < tomo.clustering_Nbin; i ++){
      for (j = i; j < tomo.clustering_Nbin; j++){
        N[i][j] = n;
        N[j][i] = n;
        n++;
      }
    }
  }
  return N[z1][z2];
}

/******************** routines for redshift distributions, including photo-zs (optional) ********************/

double zdistr_histo_n(double z,  void *params) // return nz(z,j) based on redshift file with structure z[i] nz[0][i] .. nz[tomo.shear_Nbin-1][i]
{
  double *array = (double*)params;
  static double **tab;
  FILE *ein;
  static double zhisto_max,zhisto_min,dz;

  if (tab==0){
    double *z_v;
    int i,k,zbins;
    zbins = line_count(redshift.shear_REDSHIFT_FILE);
    tab=create_double_matrix(0,tomo.shear_Nbin-1,0, zbins-1);
    z_v=create_double_vector(0, zbins-1);
    ein=fopen(redshift.shear_REDSHIFT_FILE,"r");
    for (i=0;i<zbins;i++){
      fscanf(ein, "%le", &z_v[i]);
      if (i > 0 && z_v[i] < z_v[i-1]){
        break;
      }
      for (k = 0; k < tomo.shear_Nbin; k++){
        fscanf(ein,"%le",&tab[k][i]);
      }
    }
    fclose(ein);
    dz = (z_v[i-1]-z_v[0])/(1.*i -1.);
    zhisto_max=z_v[i-1]+dz;
    zhisto_min=z_v[0];
    zbins = i;
    // now, set tomography bin boundaries
    for (k = 0; k < tomo.shear_Nbin; k++){
      double max = tab[k][0];
      for (i = 1; i < zbins; i++){
        if (tab[k][i]> max){max = tab[k][i];}
      }
      i = 0;
      while (tab[k][i] <1.e-8*max && i < zbins-2){i++;}
      tomo.shear_zmin[k] = z_v[i];
      i = zbins-1;
      while (tab[k][i] <1.e-8*max && i > 0){i--;}
      tomo.shear_zmax[k] = z_v[i];
      printf("tomo.shear_zmin[%d] = %.3f,tomo.shear_zmax[%d] = %.3f\n",k,tomo.shear_zmin[k],k,tomo.shear_zmax[k]);
    }
    free_double_vector(z_v,0,zbins-1);
    if (zhisto_max < tomo.shear_zmax[tomo.shear_Nbin-1] || zhisto_min > tomo.shear_zmin[0]){
      printf("zhisto_min = %e,zhisto_max = %e\n",zhisto_min,zhisto_max);
      printf("tomo.shear_zmin[0] = %e, tomo.shear_zmax[N-1] = %e\n", tomo.shear_zmin[0],tomo.shear_zmax[tomo.shear_Nbin-1]);
      printf("Error in redshift.c:zdistr_histo_n: %s parameters incompatible with tomo.shear bin choice\nEXIT!\n",redshift.shear_REDSHIFT_FILE);
      exit(1);
    }
  }

  if ((z>=zhisto_min) &&(z<zhisto_max)){
    return tab[(int)array[0]][(int)floor((z-zhisto_min)/dz)];
  }
  return 0.;
}

double zdistr_photoz(double zz,int j) //returns n(ztrue | j), works only with binned distributions; j>= 0 -> tomography bin j
{
  static double **table = 0;
  static double *z_v = 0;
  static double da = 0.0;
  static double zhisto_max,zhisto_min;
  static nuisancepara N;
  static int zbins = -1;
  static gsl_spline * photoz_splines[21];
  static gsl_interp_accel * photoz_accel[21];

  if (table==0){
    update_nuisance(&N);
    int zbins1 = line_count(redshift.shear_REDSHIFT_FILE);
    if(redshift.shear_photoz !=4){zbins = zbins1*20;}
    else {zbins = zbins1;}
    table   = create_double_matrix(0, tomo.shear_Nbin, 0, zbins-1);
    z_v=create_double_vector(0, zbins-1);
    for (int i = 0; i < tomo.shear_Nbin+1; i++){
      photoz_splines[i] = gsl_spline_alloc(Z_SPLINE_TYPE, zbins);
      photoz_accel[i] = gsl_interp_accel_alloc();
    }

    FILE *ein;
    double space;
    int i,k;
    ein=fopen(redshift.shear_REDSHIFT_FILE,"r");
    for (i=0;i<zbins1;i++){
      fscanf(ein, "%le", &z_v[i]);
      if (i > 0 && z_v[i] < z_v[i-1]){break;}
      for (k = 0; k < tomo.shear_Nbin; k++){fscanf(ein,"%le",&space);}
    }
    fclose(ein);
    redshift.shear_zdistrpar_zmin = fmax(z_v[0],1.e-5);
    redshift.shear_zdistrpar_zmax = z_v[i-1] +(z_v[i-1]-z_v[0])/(zbins1-1.);
    printf("redshift_spline.c %d %e %e\n", zbins,redshift.shear_zdistrpar_zmin,redshift.shear_zdistrpar_zmax);


    zhisto_max =redshift.shear_zdistrpar_zmax;
    zhisto_min = redshift.shear_zdistrpar_zmin;
    da = (zhisto_max-zhisto_min)/(1.0*zbins);
    for (int i = 0;i < zbins; i++){
      z_v[i] = zhisto_min+(i+0.5)*da;
    }
    double array[4], NORM[21],norm,x1,x2,eta,outfrac;
    //the outlier fraction (outfrac) should be specified externally. This is a temporary hack.
    for (i = 0; i< tomo.shear_Nbin; i++){
      array[0] =tomo.shear_zmin[i];
      array[1] =tomo.shear_zmax[i];
      if(redshift.shear_photoz == 4){// histogram file contains n(z) estimates for each bin
        array[0] = 1.0*i;
        //call zdistr_histo_n once to update bin boundaries
        zdistr_histo_n(0.,(void*) array);
   /*     norm = int_gsl_integrate_medium_precision(zdistr_histo_n, (void*)array, tomo.shear_zmin[i],tomo.shear_zmax[i],NULL, 1024);// *(1.-nuisance.bias_zphot_shear[i]);
        if (norm == 0){
          printf("redshift.c:zdistr_photoz:norm(nz=%d)=0\nEXIT\n",i);
          exit(1);
        }*/
        norm = 0.0;
        for (k = 0;k<zbins; k++){
          table[i+1][k] = zdistr_histo_n(z_v[k],(void*)array);//norm;
          norm += table[i+1][k]*da;
        }
        for (k = 0;k<zbins; k++){table[i+1][k]/= norm;}
      }
      else{
          printf("redshift.shear_photoz = %d not supported in this cosmolike version\n",redshift.shear_photoz);
          exit(1);
      }
      NORM[i] = norm;
    }
    // calculate normalized overall redshift distribution (without bins), store in table[0][:]
    norm = 0;
    for (i = 0; i <tomo.shear_Nbin; i++){norm += NORM[i];}

    for (k = 0;k<zbins; k++){
      table[0][k] = 0;
      for (i = 0; i <tomo.shear_Nbin; i++){table[0][k]+= table[i+1][k]*NORM[i]/norm;}
    }
    for (i = -1; i < tomo.shear_Nbin; i++){
      gsl_spline_init(photoz_splines[i+1], z_v, table[i+1], zbins);
//      printf("spline init shear %d, %e\n",i,gsl_spline_eval(photoz_splines[i+1],1.0,NULL));
    }
  }
  if (j >= tomo.shear_Nbin){
    printf("redshift.c: zdistr_photoz(z,%d) outside tomo.shear_Nbin range\n", j);
    exit(1);
  }
  zz = zz -nuisance.bias_zphot_shear[j];
  if (zz <= z_v[0] || zz >= z_v[zbins-1]) return 0.0;
  return gsl_spline_eval(photoz_splines[j+1],zz,photoz_accel[j+1]);
}



////////////// End shear routines //////////////////////

double pf_histo(double z, void *params) //return pf(z) based on redshift file with one redshift distribution
{
  static double *tab = 0;
  FILE *ein;

  static double zhisto_max,zhisto_min,dz;

  if (tab==0){
    double *z_v,space1,space2;
    int i,zbins;
    zbins = line_count(redshift.clustering_REDSHIFT_FILE);
    tab=create_double_vector(0, zbins-1);
    z_v=create_double_vector(0, zbins-1);
    ein=fopen(redshift.clustering_REDSHIFT_FILE,"r");
    for (i=0;i<zbins;i++){
      fscanf(ein,"%le %le %le %le\n",&z_v[i],&space1,&space2,&tab[i]);
      if (i > 0 && z_v[i] < z_v[i-1]){break;}
    }
    fclose(ein);
    dz = (z_v[i-1]-z_v[0])/(1.*i-1.);
    zhisto_max=z_v[i-1]+dz;
    zhisto_min=z_v[0];
    // redshift.clustering_zdistrpar_zmin = zhisto_min;
    // redshift.clustering_zdistrpar_zmax = zhisto_max;
    free_double_vector(z_v,0,zbins-1);
    if (zhisto_max < tomo.clustering_zmax[tomo.clustering_Nbin-1] || zhisto_min > tomo.clustering_zmin[0]){
      printf("Error in redshift_spline.c:pf_histo.c: %s parameters incompatible with tomo.clustering bin choice\nEXIT!\n",redshift.clustering_REDSHIFT_FILE);
      exit(1);
    }
  }

  if ((z>=zhisto_min) &&(z<zhisto_max)){
    return tab[(int)floor((z-zhisto_min)/dz)];
  }
  return 0.;
}


double pf_histo_n(double z,  void *params) //return pf(z,j) based on redshift file with structure z[i] nz[0][i] .. nz[tomo.clustering_Nbin-1][i]
{
  double *array = (double*)params;
  static double **tab;
  FILE *ein;
  static double zhisto_max,zhisto_min,dz;

  if (tab==0){
    double *z_v;
    int i,k,zbins;
    zbins = line_count(redshift.clustering_REDSHIFT_FILE);
    tab=create_double_matrix(0,tomo.clustering_Nbin-1,0, zbins-1);
    z_v=create_double_vector(0, zbins-1);
    ein=fopen(redshift.clustering_REDSHIFT_FILE,"r");
    for (i=0;i<zbins;i++){
      fscanf(ein, "%le", &z_v[i]);
      if (i > 0 && z_v[i] < z_v[i-1]){break;}
      for (k = 0; k < tomo.clustering_Nbin; k++){
        fscanf(ein," %le",&tab[k][i]);
      }
    }
    fclose(ein);
    dz =(z_v[i-1]-z_v[0])/(1.*i-1.);
    zhisto_max=z_v[i-1]+dz;
    zhisto_min=z_v[0];
    // now, set tomography bin boundaries
    for (k = 0; k < tomo.clustering_Nbin; k++){
      double max = tab[k][0];
      for (i = 1; i < zbins; i++){
        if (tab[k][i]> max){max = tab[k][i];}
      }
      i = 0;
      while (tab[k][i] <1.e-8*max && i < zbins-2){i++;}
      tomo.clustering_zmin[k] = z_v[i];
      i = zbins-1;
      while (tab[k][i] <1.e-8*max && i > 0){i--;}
      tomo.clustering_zmax[k] = z_v[i];
      printf("tomo.clustering_zmin[%d] = %.3f,tomo.clustering_zmax[%d] = %.3f\n",k,tomo.clustering_zmin[k],k,tomo.clustering_zmax[k]);
    }

    free_double_vector(z_v,0,zbins-1);
    if (zhisto_max < tomo.clustering_zmax[tomo.clustering_Nbin-1] || zhisto_min > tomo.clustering_zmin[0]){
      printf("%e %e   %e %e\n",zhisto_min,tomo.clustering_zmin[0],zhisto_max,tomo.clustering_zmax[tomo.clustering_Nbin-1]);
      printf("Error in redshift.c:pf_histo_n.c: %s parameters incompatible with tomo.clustering bin choice\nEXIT!\n",redshift.clustering_REDSHIFT_FILE);
      exit(1);
    }
  }

  if ((z>=zhisto_min) &&(z<zhisto_max)){
    return tab[(int)array[0]][(int)floor((z-zhisto_min)/dz)];
  }
  return 0.;
}
double pf_photoz(double zz,int j) //returns n(ztrue, j), works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
{
  static double **table = 0;
  static double *z_v = 0;
  static double da = 0.0;
  static double zhisto_max,zhisto_min;
  static nuisancepara N;
  static int zbins = -1;
  static gsl_spline * photoz_splines[21];
  static gsl_interp_accel * photoz_accel[21];

  if (table==0){
    update_nuisance(&N);
    zbins = line_count(redshift.clustering_REDSHIFT_FILE);
    table   = create_double_matrix(0, tomo.clustering_Nbin, 0, zbins-1);
    z_v=create_double_vector(0, zbins-1);
    for (int i = 0; i < tomo.clustering_Nbin+1; i++){
      photoz_splines[i] = gsl_spline_alloc(Z_SPLINE_TYPE, zbins);
      photoz_accel[i] = gsl_interp_accel_alloc();
    }

    FILE *ein;
    double space;
    int i,k;
    ein=fopen(redshift.clustering_REDSHIFT_FILE,"r");
    for (i=0;i<zbins;i++){
      fscanf(ein, "%le", &z_v[i]);
      if (i > 0 && z_v[i] < z_v[i-1]){break;}
      for (k = 0; k < tomo.clustering_Nbin; k++){fscanf(ein,"%le",&space);}
    }
    fclose(ein);
    redshift.clustering_zdistrpar_zmin = fmax(z_v[0],1.e-5);
    redshift.clustering_zdistrpar_zmax = z_v[i-1] +(z_v[i-1]-z_v[0])/(zbins-1.);

    zhisto_max = redshift.clustering_zdistrpar_zmax;
    zhisto_min = redshift.clustering_zdistrpar_zmin;
    da = (zhisto_max-zhisto_min)/(1.0*zbins);
    for (int i = 0;i < zbins; i++){
      z_v[i] = zhisto_min+(i+0.5)*da;
    }

    double array[4], NORM[21],norm,x1,x2,eta,outfrac,zi;
    //the outlier fraction (outfrac) should be specified externally. This is a temporary hack.

    for (i = 0; i< tomo.clustering_Nbin; i++){
      array[0] =tomo.clustering_zmin[i];
      array[1] =tomo.clustering_zmax[i];
      if(redshift.clustering_photoz == 4){// histogram file contains n(z) estimates for each bin
        array[0] = 1.0*i;
        pf_histo_n(0.,(void*) array);
        norm = int_gsl_integrate_medium_precision(pf_histo_n, (void*)array, tomo.clustering_zmin[i],tomo.clustering_zmax[i],NULL, 1024);
        if (norm == 0){
          printf("redshift.c:pf_photoz:norm(nz=%d)=0\nEXIT\n",i);
          exit(1);
        }
        for (k = 0;k<zbins; k++){ table[i+1][k] = pf_histo_n(z_v[k],(void*)array)/norm;}
      }
      else{
        printf("redshift.clustering_photoz = %d not supported in this cosmolike version\n",redshift.clustering_photoz);
        exit(1);
      }
      NORM[i] = norm;
    }
    // calculate normalized overall redshift distribution (without bins), store in table[0][:]
    norm = 0;
    for (i = 0; i <tomo.clustering_Nbin; i++){norm += NORM[i];}
    for (k = 0,zi = zhisto_min;k<zbins; k++,zi+=da){
      table[0][k] = 0;
      for (i = 0; i <tomo.clustering_Nbin; i++){table[0][k]+= table[i+1][k]*NORM[i]/norm;}
    }
      for (i = -1; i < tomo.clustering_Nbin; i++){
      //for(k = 0; k<zbins; k++){printf("%d %d %e %e\n",i,k,z_v[k],table[i+1][k]);}
      gsl_spline_init(photoz_splines[i+1], z_v, table[i+1], zbins);
    }

  }
  if (j >= tomo.clustering_Nbin){
    printf("redshift.c: pf_photoz(z,%d) outside tomo.clustering_Nbin range\n", j);
    exit(1);
  }
  if (redshift.clustering_photoz == 4){ zz = zz -nuisance.bias_zphot_clustering[j];}
  if (zz <= z_v[0] || zz >= z_v[zbins-1]) return 0.0;
  return gsl_spline_eval(photoz_splines[j+1],zz,photoz_accel[j+1]);
}

/*********** routines calculating the number of source and lens galaxies per bin ****************/

double nsource(int j) //returns n_gal for shear tomography bin j, works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
{
  static double **table = 0;

  double array[3];
  int i;
  if (table ==0|| table[0][0] != survey.n_gal){
    array[0] = zdistr_photoz(0.,0);
    table   = create_double_matrix(0, tomo.shear_Nbin, 0, 1);
    printf("redshift.shear_photoz = 4, using tabulated tomo.n_source =");
    for (i = 0; i< tomo.shear_Nbin; i++){
      printf(" %e,", tomo.n_source[i]);
      table[i+1][0] = tomo.n_source[i];
      if (tomo.n_source[i] < 0.01 || tomo.n_source[i] > 100.){
        printf("\n\n!!!!!!!!!\n\nredshift.shear_photoz = 4, tomo.n_source[%d] = %e, EXIT\n\n!!!!!!!!!\n",i,tomo.n_source[i]);
        exit(0);
      }
    }
    printf("\n");
    table[0][0] = survey.n_gal;
  }
  return table[j+1][0];
}


double nlens(int j) //returns n_gal for clustering tomography bin j, works only with binned distributions; j =-1 -> no tomography; j>= 0 -> tomography bin j
{
  static double **table =0;
  double array[3];
  int i;
  if (table ==0 ||table[0][0] != survey.n_lens){
    array[0] = pf_photoz(0.,0);
    table   = create_double_matrix(0, tomo.clustering_Nbin, 0, 1);
    printf("redshift.clustering_photoz = 4, using tabulated tomo.n_lens_ij[i][i] =");
    for (i = 0; i< tomo.clustering_Nbin; i++){
      printf(" %e,", tomo.n_lens_ij[i][i]);
      table[i+1][0] = tomo.n_lens_ij[i][i];
    }
    printf("\n");
    table[0][0] = survey.n_lens;
  }
  return table[j+1][0];
}


double int_for_zmean(double z, void *params){
  double *array = (double*)params;
  return z*pf_photoz(z,(int)array[0]);
}

double norm_for_zmean(double z, void *params){
  double *array = (double*)params;
  return pf_photoz(z,(int)array[0]);
}
double zmean(int j){ //mean true redshift of galaxies in tomography bin j
  static double **table = 0;
  if (table ==0){
    double array[1];
    array[0] = pf_photoz(0.,0);
    table   = create_double_matrix(0, tomo.clustering_Nbin, 0, 1);
    for (int i = 0; i< tomo.clustering_Nbin; i++){
     array[0]  = 1.0*i;
     table[i][0] = int_gsl_integrate_low_precision(int_for_zmean, (void*)array, tomo.clustering_zmin[i],tomo.clustering_zmax[i],NULL, 1024)/int_gsl_integrate_low_precision(norm_for_zmean, (void*)array, tomo.clustering_zmin[i],tomo.clustering_zmax[i],NULL, 1024);
    }
  }
  return table[j][0];
}

double int_for_zmean_source(double z, void *params){

  double *array = (double*)params;
  return z*zdistr_photoz(z,(int)array[0]);
}
double zmean_source(int j){ //mean true redshift of source galaxies in tomography bin j
  static double **table = 0;
  if (table ==0){

    double array[1];
    array[0] = zdistr_photoz(0.,0);
    int i;
    table   = create_double_matrix(0, tomo.shear_Nbin, 0, 1);
    for (i = 0; i< tomo.shear_Nbin; i++){
     array[0]  = 1.0*i;
     table[i][0] = int_gsl_integrate_medium_precision(int_for_zmean_source, (void*)array, tomo.shear_zmin[i],tomo.shear_zmax[i],NULL, 1024);
    }
  }
  return table[j][0];
}

/************ integrands for bin-averaged lens efficiencies **************/

double int_for_g_tomo(double aprime,void *params)
{
  double chi1, chi_prime,val;
  double *ar = (double *) params;
  int zbin= (int) ar[0];
  chi1 = chi(ar[1]);
  chi_prime = chi(aprime);

  val=zdistr_photoz(1./aprime-1.,zbin)*f_K(chi_prime-chi1)/f_K(chi_prime)/(aprime*aprime);
  return val;
}

double g_tomo(double a, int zbin) // for tomography bin zbin
{
  static nuisancepara N;
  static cosmopara C;

  static double **table = 0;
  static double da = 0.0;
  double aa;
  int i,j;
  double array[2];
  if (table ==0 || recompute_zphot_shear(N) || recompute_expansion(C)){
    if (table==0) table   = create_double_matrix(0, tomo.shear_Nbin, 0, Ntable.N_a-1);
    da = (0.999999-1./(redshift.shear_zdistrpar_zmax+1.))/(Ntable.N_a-1);
    for (j=-1;j<tomo.shear_Nbin;j++) {
      array[0]=(double) j; //if j=-1, no tomography is being done
      aa = 1./(redshift.shear_zdistrpar_zmax+1.);
      for (i=0;i<Ntable.N_a;i++,aa+=da) {
        array[1] = aa;
        table[j+1][i] = int_gsl_integrate_medium_precision(int_for_g_tomo,(void*)array,1./(redshift.shear_zdistrpar_zmax+1.),aa,NULL,4000);
      }
    }
    update_nuisance(&N);
    update_cosmopara(&C);
  }
  if (a<=1./(redshift.shear_zdistrpar_zmax+1.) || a>1.0-da) return 0.0;
  return interpol(table[zbin+1], Ntable.N_a, 1./(redshift.shear_zdistrpar_zmax+1.), 0.999999, da, a, 1.0, 1.0); //zbin =-1 is non-tomography
}


double int_for_g_lens(double aprime,void *params)
{
  double chi1, chi_prime,val;
  double *ar = (double *) params;
  int zbin= (int) ar[0];
  chi1 = chi(ar[1]);
  chi_prime = chi(aprime);

  val=pf_photoz(1./aprime-1.,zbin)*f_K(chi_prime-chi1)/f_K(chi_prime)/(aprime*aprime);
  return val;
}

double g_lens(double a, int zbin) // for *lens* tomography bin zbin
{
  static nuisancepara N;
  static cosmopara C;

  static double **table = 0;
  static double da = 0.0;
  double aa;
  int i,j;
  double array[2];
  if (table ==0 || recompute_zphot_clustering(N) || recompute_expansion(C)){
    if (table==0) table   = create_double_matrix(0, tomo.clustering_Nbin, 0, Ntable.N_a-1);
    da = (0.999999-1./(redshift.clustering_zdistrpar_zmax+1.))/(Ntable.N_a-1);
    for (j=-1;j<tomo.clustering_Nbin;j++) {
      array[0]=(double) j; //if j=-1, no tomography is being done
      aa = 1./(redshift.clustering_zdistrpar_zmax+1.);
      for (i=0;i<Ntable.N_a;i++,aa+=da) {
        array[1] = aa;
        table[j+1][i] = int_gsl_integrate_medium_precision(int_for_g_lens,(void*)array,1./(redshift.shear_zdistrpar_zmax+1.),aa,NULL,4000);
      }
    }
    update_nuisance(&N);
    update_cosmopara(&C);
  }
  if (a<=1./(redshift.clustering_zdistrpar_zmax+1.) || a>1.0-da) return 0.0;
  return interpol(table[zbin+1], Ntable.N_a, 1./(redshift.clustering_zdistrpar_zmax+1.), 0.999999, da, a, 1.0, 1.0); //zbin =-1 is non-tomography
}
