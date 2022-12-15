/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by CosmoLike developers
******************************************************************************/

typedef struct {
  int Ncl;
  int Ntheta;
  int Ncos;
  int Ndata;
  double lmin;
  double lmax;
  double vtmax;
  double vtmin;
  double *theta_min;
  double cosmax;
  double Rmin_bias;
  double Rmin_shear;
  double lmax_shear;
  int baryons;
  int IA;
  int bias;
  int wlphotoz;
  int clphotoz;
  int shearcalib;
  char DATA_FILE[500];
  char INV_FILE[500];
  char COV_FILE[500];
  char BARY_FILE[500];
  char MASK_FILE[500];
  int shear_shear;
  int shear_pos;
  int pos_pos;
  char probes[500];
  char ext_data[500];
  int theta_s;
}likepara;
likepara like ={.theta_min = 0, .baryons = 0, .IA = 0., .bias = 0, .wlphotoz = 0, .clphotoz = 0, .shearcalib = 0,.theta_s =0};

typedef struct {
     double Omega_m;  /* matter density parameter                       */
     double Omega_v;  /* cosmogical constant parameter                  */
     double sigma_8;  /* power spectrum normalization                   */
     double A_s;
     double n_spec;   /* spectral index of initial power spectrum       */
     double alpha_s;   /* running of spectral index of initial power spectrum       */
     double w0; //time dependent Dark energy parametrization zero order
     double wa; //time dependent Dark energy parametrization first order
     double omb; //Omega baryon
     double h0; //Hubble constant
     double M_nu;
     double Omega_nu; //density parameter of massive neutrinos; Omega_m = Omega_cdm+ Omega_nu + omb
     double coverH0; //units for comoving distances - speeds up code
     double rho_crit;      /* = 3 H_0^2/(8 pi G), critical comoving density */
     double f_NL;
     double MGSigma;
     double MGmu;
     double theta_s;
}cosmopara;
cosmopara cosmology = {.A_s = 0., .sigma_8=0., .alpha_s =0.0, .M_nu =0., .Omega_nu =0.,.coverH0= 2997.92458, .rho_crit = 7.4775e+21,.MGSigma=0.0,.MGmu=0.0,.theta_s =0.0};

typedef struct {
  int shear_Nbin; // number of tomography bins
  int shear_Npowerspectra;// number of tomography power spectra+2+3+...+Nbin
  double shear_zmax[20]; // code needs modification if more than 10 zbins
  double shear_zmin[20];
  double n_source[20];
  int clustering_Nbin; // number of tomography bins
  int clustering_Npowerspectra;// number of tomography power spectra+2+3+...+Nbin
  double clustering_zmax[20];
  double clustering_zmin[20];
  // double n_lens[20];
  double n_lens_ij[20][20];
  int ggl_Npowerspectra;// number of ggl tomography combinations
}tomopara;
tomopara tomo = {.n_source = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}};

typedef struct {
  int shear_photoz;
  double shear_zdistrpar_zmin;
  double shear_zdistrpar_zmax;
  int shear_histogram_zbins;
  char shear_REDSHIFT_FILE[200];

  int clustering_photoz;
  double clustering_zdistrpar_zmin;
  double clustering_zdistrpar_zmax;
  int clustering_histogram_zbins;
  char clustering_REDSHIFT_FILE[200];
}redshiftpara;
redshiftpara redshift;

typedef struct {
     double area;/* survey_area in deg^2. */
     double n_gal;/* galaxy density per arcmin^2 */
     double sigma_e;/* rms inrinsic ellipticity noise*/
     double area_conversion_factor; /*factor from deg^2 to radian^2: 60*60*constants.arcmin*constants.arcmin */
     double n_gal_conversion_factor; /*factor from n_gal/arcmin^2 to n_gal/radian^2: 1.0/constants.arcmin/constants.arcmin */
     double n_lens;/* lens galaxy density per arcmin^2 */
     double m_lim;
     char name[500];
     int surveystage;
     char sourcephotoz[256];
     char lensphotoz[256];
     char galsample[256];
}sur;

sur survey = {.area_conversion_factor = 60.0*60.0*2.90888208665721580e-4*2.90888208665721580e-4, .n_gal_conversion_factor = 1.0/2.90888208665721580e-4/2.90888208665721580e-4};

double bgal_z(double z, int nz);
double b1_per_bin(double z, int nz);

typedef  double (*B1_model)(double z, int nz);
typedef struct{
  double b[20]; /* linear galaxy bias paramter in clustering bin i*/
  double b2[20]; /* quadratic bias parameter for redshift bin i */
  double bs2[20]; /* leading order tidal bias for redshift bin i */
  double b_mag[20]; /*amplitude of magnification bias, b_mag[i] = 5*s[i]+beta[i] -2 */
  B1_model b1_function;
}galpara;
galpara gbias ={.b2 ={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                .bs2 ={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
                .b1_function = &b1_per_bin, 
                .b_mag ={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0}}; //default: point to old bgal_z routin
double b1_per_bin(double z, int ni){
  return gbias.b[ni];
}

typedef struct {
  char runmode[300];
  char baryons[300];
  double DIFF_n; //difference fucntion describing the constant uncertainty in Pdelta for k>0.01
  double DIFF_A; //difference fucntion describing the scale dependent uncertainty in Pdelta for k>0.01
}pdeltapara;
pdeltapara pdeltaparams = {.runmode = "Halofit", .DIFF_n = 0., .DIFF_A = 0.};

typedef struct {
  double A_ia; //A IA see Joachimi2012
  double eta_ia; //eta_other IA see Joachimi2012
  double oneplusz0_ia; //oneplusz0-ia MegaZ
  double c1rhocrit_ia;
  double fred[20];
  double shear_calibration_m[20];
  double sigma_zphot_shear[20];
  double bias_zphot_shear[20];
  double sigma_zphot_clustering[20];
  double bias_zphot_clustering[20];
  double sigma_zphot_magnification[20];
  double bias_zphot_magnification[20];
}
nuisancepara;
nuisancepara nuisance ={.c1rhocrit_ia = 0.013873073650776856,
  .shear_calibration_m = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
  .sigma_zphot_shear = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
  .bias_zphot_shear = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
  .sigma_zphot_clustering = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
  .bias_zphot_clustering = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}};

typedef struct {
    double tmin; /* Theta min (arcmin) */
    double tmax; /* Theta max (arcmin) */
    int ntheta;  /* number of theta bins */
    int lin_bins;/* switch between log-binning (lin_bins = 0, default) and linear binning (lin_bins =1)*/
    int full_tomo;/* switch between auto-clustering-correlation (full_tomo = 0, default) and cross-clustering-correlation (full_tomo=1) */
    double lmin; /* ell min  */
    double lmax; /* ell max  */
    int ncl;/* number of ell bins */
    int ng;/* ng covariance? */
    int cng;/* cng covariance */
    char outdir[200]; /* output directory */
    char filename[200]; /* output file name prefix */
    char C_FOOTPRINT_FILE[200]; /*angular power spectrum of survey footprint, in healpix format */
    char ss[8]; /* Calculate shear-shear components */
    char ls[8]; /* Calculate shear-position components */
    char ll[8]; /* Calculate position-position components */
} covpar;
covpar covparams = {.lin_bins = 0, .full_tomo = 0, .ss ="false", .ls ="false", .ll ="false"};
