/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by Alex Barreira and Elisabeth Krause
******************************************************************************/

#include <dirent.h>
#include <python2.7/Python.h>
/*******look-up tables for Fourier-space WL covariances using responses********/
/******************************************************************************/
/******************************************************************************/
//WL super-sample covariance, using response model + non-Limber calculation of 1711.07467
//for computational efficiency, switch to Limber calculation above L = L_max (set below)
double bin_cov_SSC_response_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4);
//WL cNG covariance, using response model of 1703.09212 and 1705.01092
double bin_cov_NG_response_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4);

// L_max for non-Limber calculation; L_max = 0 for full Limber
int L_MAX = 5;

// These values are the maximum values for which there are response measurements (uncommenting these will use these zmax values, even for z values beyond the simulation measurements)
double z_max_covNG_response = 3.0;
double z_max_R_1_coef       = 3.0;
double z_max_R_K_coef       = 2.629;
double z_max_plpm           = 3.; //This is set to three just because it was the default value for runs using the response values in this set.
/******************** input for python response calculation *******************/
/******************************************************************************/
/******************************************************************************/
void write_pl_multiz(char *filename,char *filenamez,int Na,int Nk_per_dec){
  //create filename from filebas, processid, Omega_m and A_ia value to ensure unique file names
  FILE *F;
  F = fopen(filename,"w");
  double a[100];
  fprintf(F,"#k[h/Mpc] P_L(k,z)\n#z={");
  for (int i = 0; i < Na; i++){
    a[i] = 1. -(z_max_plpm/(z_max_plpm+1))*i/(Na-1.);
    fprintf(F, "%e",1./a[i]-1.);
    if (i < Na-1){fprintf(F, ", ");}
  }
   fprintf(F,"}\n");

  double dlgk = log(10.)/(double) Nk_per_dec;
  for (double lgk = log(1.e-5), i = 0; lgk < log(1.e+3); i++,lgk+=dlgk){
    double k = exp(lgk);
    fprintf(F,"%e",k);
    for (int i = 0; i < Na; i++){
      fprintf(F," %e", p_lin(k*cosmology.coverH0,a[i])*pow(cosmology.coverH0,3.0));
    }
   fprintf(F,"\n");
  }
  fclose(F);
  F = fopen(filenamez,"w");
  for (int i = 0; i < Na; i++){
    fprintf(F, "%e\n",1./a[i]-1.);
  }
  fclose(F);
}

void write_pm_multiz(char *filename,char *filenamez, int Na,int Nk_per_dec){
  //create filename from filebas, processid, Omega_m and A_ia value to ensure unique file names
  FILE *F;
  F = fopen(filename,"w");
  double a[100];
  fprintf(F,"#k[h/Mpc] P_m(k,z)\n#z={");
  for (int i = 0; i < Na; i++){
    a[i] = 1. -(z_max_plpm/(z_max_plpm+1))*i/(Na-1.);
    fprintf(F, "%e",1./a[i]-1.);
    if (i < Na-1){fprintf(F, ", ");}
  }
   fprintf(F,"}\n");

  double dlgk = log(10.)/(double) Nk_per_dec;
  for (double lgk = log(1.e-5), i = 0; lgk < log(1.e+3); i++,lgk+=dlgk){
    double k = exp(lgk);
    fprintf(F,"%e",k);
    for (int i = 0; i < Na; i++){
      fprintf(F," %e", Pdelta(k*cosmology.coverH0,a[i])*pow(cosmology.coverH0,3.0));
    }
    fprintf(F,"\n");
  }
  fclose(F);
  F = fopen(filenamez,"w");
  for (int i = 0; i < Na; i++){
    fprintf(F, "%e\n",1./a[i]-1.);
  }
  fclose(F);
}

/********************* interfaces to response calculation *********************/
/******************************************************************************/
/******************************************************************************/
//interface to response covNG python module, based on arXiv:1703.09212 and arXiv:1705.01092
double covNG_response(double k1, double k2, double a_ext){
  static int init = 1;
  static PyObject * covNG_response_py = NULL;
  static cosmopara C;
  double a = fmax(a_ext,1./(1.+z_max_covNG_response));
  if (init){
    DIR* dir = opendir("../cosmolike_core/theory/lookup_tables");
    if (dir){closedir(dir);}
    else{mkdir("../cosmolike_core/theory/lookup_tables",0777);}
    //write power spectra for current cosmology
    write_pl_multiz("../cosmolike_core/theory/lookup_tables/P_L_class.dat","../cosmolike_core/theory/lookup_tables/P_L_zvalues.dat",15,20);
    write_pm_multiz("../cosmolike_core/theory/lookup_tables/P_m_runmode.dat","../cosmolike_core/theory/lookup_tables/P_m_zvalues.dat",15,20);
    int ret = system("python ../cosmolike_core/theory/write_Pderivs.py");
    if (ret){
      printf("cov_response_WL.c:covNG_response: call to write_Pderivs.py failed\nEXIT\n");
      exit(1);
    }
    // Start python
    Py_Initialize();
    PyRun_SimpleString("import sys; sys.path.append('../cosmolike_core/theory')");
    init = 0;
    // Name of module and function to import and the number of args it has
    char * module_name_c = "interface_cov_response";
    char * function_name_c = "call_cov_response";
    // Convert C string to python string.
    // PyRun_SimpleString("import sys; sys.path.append('.')");
    PyObject * module_name_py = PyString_FromString(module_name_c);

    // Import the module
    PyObject * module = PyImport_Import(module_name_py);

    // Free module_name_py
    Py_DECREF(module_name_py);
    if (module == NULL) {
      fprintf(stderr,"cov_response_WL.c:covNG_response: Error: Module %s not found\n",module_name_c);
      fprintf(stderr,"add CosmoCov/cosmolike_core/theory to PYTHONPATH\n");
      PyErr_Print();
      exit(1);
    }
    // Get the function out of the module
    covNG_response_py = PyObject_GetAttrString(module, function_name_c);

    // Check that the function was found.
    if (covNG_response_py == NULL) {
      fprintf(stderr,"cov_response_WL.c:covNG_response: Error: Function %s not found in module %s\n",
            function_name_c,module_name_c);
      exit(1);
    }
    update_cosmopara(&C);
  }
  //cosmology changed
  if (recompute_cosmo3D(C)){
    //write power spectra for currently cosmology
    write_pl_multiz("../cosmolike_core/theory/lookup_tables/P_L_class.dat","../cosmolike_core/theory/lookup_tables/P_L_zvalues.dat",15,20);
    write_pm_multiz("../cosmolike_core/theory/lookup_tables/P_m_runmode.dat","../cosmolike_core/theory/lookup_tables/P_m_zvalues.dat",15,20);
    int ret = system("python ../cosmolike_core/theory/write_Pderivs.py");
    if (ret){
      printf("cov_response_WL.c:covNG_response: call to write_Pderivs.py failed\nEXIT\n");
      exit(1);
    }
    printf("cov_response_WL.c:covNG_response: cosmology changed, still need to implement update of look-up tables in covNG_resp.py\nEXIT\n");
    exit(1);
    update_cosmopara(&C);
  }
  // Create the arguments to be passed into the function
  int number_arguments = 3;
  PyObject * arguments = PyTuple_New(number_arguments);
  PyTuple_SetItem(arguments, 0, PyFloat_FromDouble(k1/cosmology.coverH0));
  PyTuple_SetItem(arguments, 1, PyFloat_FromDouble(k2/cosmology.coverH0));
  PyTuple_SetItem(arguments, 2, PyFloat_FromDouble(1./a-1.));

  PyObject * result = PyObject_CallObject(covNG_response_py, arguments);
  // Free the arguments
  Py_DECREF(arguments);

  // Check that the function worked and print the python
  // traceback if it fails
  if (result == NULL) {
    printf("cov_response_WL.c:covNG_response: Error: call to python function failed\n");
    PyErr_Print();
    exit(1);
  }
  return PyFloat_AsDouble(result)*pow(cosmology.coverH0,-9.0);
}
//interface to R_1 python module, which G_1 results from https://arxiv.org/abs/1503.03487
double R_1_coef(double k1, double a_ext){
  static int init = 1;
  static PyObject * R_int_py = NULL;
  static cosmopara C;
  double a = fmax(a_ext,1./(1.+z_max_R_1_coef));
  if (init){
    DIR* dir = opendir("../cosmolike_core/theory/lookup_tables");
    if (dir){closedir(dir);}
    else{mkdir("../cosmolike_core/theory/lookup_tables",0777);}
    //write power spectra for current cosmology
    write_pl_multiz("../cosmolike_core/theory/lookup_tables/P_L_class.dat","../cosmolike_core/theory/lookup_tables/P_L_zvalues.dat",15,20);
    write_pm_multiz("../cosmolike_core/theory/lookup_tables/P_m_runmode.dat","../cosmolike_core/theory/lookup_tables/P_m_zvalues.dat",15,20);
    // Call python script that writes the derivative files (check if failed)
    int ret = system("python ../cosmolike_core/theory/write_Pderivs.py");
    if (ret){
      printf("cov_response_WL.c:R_1_coef: call to write_Pderivs.py failed\nEXIT\n");
      exit(1);
    }
    // Start python
    Py_Initialize();
    PyRun_SimpleString("import sys; sys.path.append('../cosmolike_core/theory')");
    init = 0;
    // Name of module and function to import
    char * module_name_c = "interface_R1RK";
    char * function_name_c = "call_R_1_int";

    // Convert C string to python string
    PyObject * module_name_py = PyString_FromString(module_name_c);
    // Import the python module
    PyObject * module = PyImport_Import(module_name_py);
    // Free module_name_py from memory
    Py_DECREF(module_name_py);
    // Check if module was loaded
    if (module == NULL){
      fprintf(stderr,"cov_response_WL.c:R_1_coef: Error: Module %s not found\n",module_name_c);
      fprintf(stderr,"add CosmoCov/cosmolike_core/theory to PYTHONPATH\n");
      PyErr_Print();
      exit(1);
    }
    // Get the functions out of the module
    R_int_py = PyObject_GetAttrString(module, function_name_c);
    // Check if the functions were found
    if (R_int_py == NULL){
            fprintf(stderr,"cov_response_WL.c:R_1_coef: Error: Function %s not found in module %s\n", function_name_c,module_name_c);
      exit(1);
    }
    update_cosmopara(&C);
  }
  // If cosmology changes
  if (recompute_cosmo3D(C)){
    //write power spectra for currently cosmology
    write_pl_multiz("../cosmolike_core/theory/lookup_tables/P_L_class.dat","../cosmolike_core/theory/lookup_tables/P_L_zvalues.dat",15,20);
    write_pm_multiz("../cosmolike_core/theory/lookup_tables/P_m_runmode.dat","../cosmolike_core/theory/lookup_tables/P_m_zvalues.dat",15,20);
    // Call python script that writes the derivative files (check if failed)
    int ret = system("python ../cosmolike_core/theory/write_Pderivs.py");
    if (ret){
      printf("cov_response_WL.c:R_1_coef: call to write_Pderivs.py failed\nEXIT\n");
      exit(1);
      }
      printf("cov_response_WL.c:R_1_coef: cosmology changed, still need to implement update of look-up tables in R1RK_SSC.py\nEXIT\n");
      exit(1);
      update_cosmopara(&C);
    }
    // Create the argument to be passed into the function
    int number_arguments = 2;
    PyObject * arguments = PyTuple_New(number_arguments);
    PyTuple_SetItem(arguments, 0, PyFloat_FromDouble(1./a - 1.));
    PyTuple_SetItem(arguments, 1, PyFloat_FromDouble(k1/cosmology.coverH0));
    // Get the result
    PyObject * result = PyObject_CallObject(R_int_py, arguments);
    //printf("python call returned %e for k1 = %e k2 = %e z = %e\n",PyFloat_AsDouble(result), k1/cosmology.coverH0, k2/cosmology.coverH0, 1./a-1.);
    // Free the arguments from memory
    Py_DECREF(arguments);
    // Check that the function worked; print the python traceback if it didn't
    if (result == NULL){
    printf("cov_response_WL.c:R_1_coef: Error: call to R_1_int python function failed\n");
    PyErr_Print();
    exit(1);
    }
    return  PyFloat_AsDouble(result);
}

//interface to R_K python module, , which uses G_K results from https://arxiv.org/abs/1803.03274
double R_K_coef(double k1, double a_ext){
  static int init = 1;
  static PyObject * R_int_py = NULL;
  static cosmopara C;
  double a = fmax(a_ext,1./(1.+z_max_R_K_coef));
  if (init){
    DIR* dir = opendir("../cosmolike_core/theory/lookup_tables");
    if (dir){closedir(dir);}
    else{mkdir("../cosmolike_core/theory/lookup_tables",0777);}
    //write power spectra for current cosmology
    write_pl_multiz("../cosmolike_core/theory/lookup_tables/P_L_class.dat","../cosmolike_core/theory/lookup_tables/P_L_zvalues.dat",15,20);
    write_pm_multiz("../cosmolike_core/theory/lookup_tables/P_m_runmode.dat","../cosmolike_core/theory/lookup_tables/P_m_zvalues.dat",15,20);
    // Call python script that writes the derivative files (check if failed)
    int ret = system("python ../cosmolike_core/theory/write_Pderivs.py");
    if (ret){
      printf("cov_response_WL.c:R_K_coef: call to write_Pderivs.py failed\nEXIT\n");
      exit(1);
    }
    // Start python
    Py_Initialize();
    PyRun_SimpleString("import sys; sys.path.append('../cosmolike_core/theory')");
    init = 0;
    // Name of module and function to import
    char * module_name_c = "interface_R1RK";
    char * function_name_c = "call_R_K_int";

    // Convert C string to python string
    PyObject * module_name_py = PyString_FromString(module_name_c);
    // Import the python module
    PyObject * module = PyImport_Import(module_name_py);
    // Free module_name_py from memory
    Py_DECREF(module_name_py);
    // Check if module was loaded
    if (module == NULL){
      fprintf(stderr,"cov_response_WL.c:R_K_coef: Error: Module %s not found\n",module_name_c);
      fprintf(stderr,"add CosmoCov/cosmolike_core/theory to PYTHONPATH\n");
      PyErr_Print();
      exit(1);
    }
    // Get the functions out of the module
    R_int_py = PyObject_GetAttrString(module, function_name_c);
    // Check if the functions were found
    if (R_int_py == NULL){
            fprintf(stderr,"cov_response_WL.c:R_K_coef: Error: Function %s not found in module %s\n", function_name_c,module_name_c);
      exit(1);
    }
    update_cosmopara(&C);
  }
  // If cosmology changes
  if (recompute_cosmo3D(C)){
    //write power spectra for currently cosmology
    write_pl_multiz("../cosmolike_core/theory/lookup_tables/P_L_class.dat","../cosmolike_core/theory/lookup_tables/P_L_zvalues.dat",15,20);
    write_pm_multiz("../cosmolike_core/theory/lookup_tables/P_m_runmode.dat","../cosmolike_core/theory/lookup_tables/P_m_zvalues.dat",15,20);
    // Call python script that writes the derivative files (check if failed)
    int ret = system("python ../cosmolike_core/theory/write_Pderivs.py");
    if (ret){
      printf("cov_response_WL.c:R_K_coef: call to write_Pderivs.py failed\nEXIT\n");
      exit(1);
      }
      printf("cov_response_WL.c:R_K_coef: cosmology changed, still need to implement update of look-up tables in R1RK_SSC.py\nEXIT\n");
      exit(1);
      update_cosmopara(&C);
    }
    // Create the argument to be passed into the function
    int number_arguments = 2;
    PyObject * arguments = PyTuple_New(number_arguments);
    PyTuple_SetItem(arguments, 0, PyFloat_FromDouble(1./a - 1.));
    PyTuple_SetItem(arguments, 1, PyFloat_FromDouble(k1/cosmology.coverH0));
    // Get the result
    PyObject * result = PyObject_CallObject(R_int_py, arguments);
    //printf("python call returned %e for k1 = %e k2 = %e z = %e\n",PyFloat_AsDouble(result), k1/cosmology.coverH0, k2/cosmology.coverH0, 1./a-1.);
    // Free the arguments from memory
    Py_DECREF(arguments);
    // Check that the function worked; print the python traceback if it didn't
    if (result == NULL){
    printf("cov_response_WL.c:R_K_coef: Error: call to R_K_int python function failed\n");
    PyErr_Print();
    exit(1);
    }
    return PyFloat_AsDouble(result);
}



/******************************************************************************/
/********************* variance calculation************************************/
/******************************************************************************/
//curved sky, healpix window function
// variance calculation, restricted to modes with L > L_MAX
// modes with L<= L_MAX are calculated separately in non-Limber calculation
double sum_variance_healpix_LMIN(double a,int L_MIN){
  double res = 0.;
  double r = f_K(chi(a));
  for (int l = L_MIN; l < 500; l++){
    res+= (2.*l+1.)/(r*r)*C_survey_window(l)*p_lin((l+0.5)/r,a)/(4.*M_PI);
  }
  return res;
}

/******************************************************************************/
/********************* full sky SSC calculation           *********************/
/******************************************************************************/
//  (R_1+R_K/6) part of Eq. 2.14 in in 1807.04266
double flLp_chi_response(double chi1, double l1, int z1, int z2){
  double k1 = (l1+0.5)/chi1;
  double aa = a_chi(chi1);
  double res;
  res = 1./6.*R_K_coef(k1,aa);
  res += R_1_coef(k1,aa);
  res *= growfac(aa)/growfac(1.0)*Pdelta(k1,aa)*W_kappa(aa,f_K(chi1),z1)*W_kappa(aa,f_K(chi1),z2)/(chi1*chi1);
  return res;
}
//  R_K/2\partial^2_x part of Eq. 2.14 in in 1807.04266
double flLp_chi_response_deriv(double chi1, double l1, int z1, int z2){
  double k1 = (l1+0.5)/chi1;
  double aa = a_chi(chi1);
  double res;
  res = 1./2.*R_K_coef(k1,aa);
  res *= growfac(aa)/growfac(1.0)*Pdelta(k1,aa)*W_kappa(aa,f_K(chi1),z1)*W_kappa(aa,f_K(chi1),z2)/(chi1*chi1);
  return res;
}
double flLp_chi_response_spline(double chi1, double l1, int z1, int z2){
  static double L1 = -1;
  static int Z1 = -1, Z2 =-1;
  static gsl_spline * f_spline = NULL;
  static gsl_interp_accel * f_accel = NULL;
  static cosmopara C;
  static double chi_max = -1;
  if (!f_spline){
    f_spline = gsl_spline_alloc(gsl_interp_cspline, Ntable.N_a);
    f_accel = gsl_interp_accel_alloc();
  }
  if (L1 != l1 || Z1 != z1 || Z2 != z2 || recompute_cosmo3D(C)){
    L1 = l1; Z1 = z1; Z2 = z2; update_cosmopara(&C);
    double aa, k1,*table_f,*table_chi;
    table_f  = create_double_vector(0, Ntable.N_a-1);
    table_chi  = create_double_vector(0, Ntable.N_a-1);
    chi_max = chi(amin_source(z1));
    //skip [0] elements, as lens efficiency is zero at chi = 0
    table_chi[0] = 0.0; table_f[0] = 0.0;
    for (int i = 1; i < Ntable.N_a; i++){
      table_chi[i] = (double)i/(Ntable.N_a-1.)*chi_max;
      table_f[i] =flLp_chi_response(table_chi[i],L1,Z1,Z2);
    }
    gsl_interp_accel_reset(f_accel);
    gsl_spline_init(f_spline, table_chi, table_f, Ntable.N_a);
    free_double_vector(table_f,0, Ntable.N_a-1);
    free_double_vector(table_chi,0, Ntable.N_a-1);
  }
  if (chi1 <=0.0 || chi1 >= chi_max){
    return 0.0;
  }
  return gsl_spline_eval(f_spline,chi1,f_accel);
}

double flLp_chi_response_deriv_spline(double chi1, double l1, int z1, int z2){
  static double L1 = -1;
  static int Z1 = -1, Z2 =-1;
  static gsl_spline * f_spline = NULL;
  static gsl_interp_accel * f_accel = NULL;
  static cosmopara C;
  static double chi_max = -1;
  if (!f_spline){
    f_spline = gsl_spline_alloc(gsl_interp_cspline, Ntable.N_a);
    f_accel = gsl_interp_accel_alloc();
  }
  if (L1 != l1 || Z1 != z1 || Z2 != z2 || recompute_cosmo3D(C)){
    L1 = l1; Z1 = z1; Z2 = z2; update_cosmopara(&C);
    double aa, k1,*table_f,*table_chi;
    table_f  = create_double_vector(0, Ntable.N_a-1);
    table_chi  = create_double_vector(0, Ntable.N_a-1);
    chi_max = chi(amin_source(z1));
    //skip [0] elements, as lens efficiency is zero at chi = 0
    table_chi[0] = 0.0; table_f[0] = 0.0;
    for (int i = 1; i < Ntable.N_a; i++){
      table_chi[i] = (double)i/(Ntable.N_a-1.)*chi_max;
      table_f[i] =flLp_chi_response_deriv(table_chi[i],L1,Z1,Z2);
    }
    gsl_interp_accel_reset(f_accel);
    gsl_spline_init(f_spline, table_chi, table_f, Ntable.N_a);
    free_double_vector(table_f,0, Ntable.N_a-1);
    free_double_vector(table_chi,0, Ntable.N_a-1);
  }
  if (chi1 <=0.0 || chi1 >= chi_max){
    return 0.0;
  }

  return gsl_spline_eval(f_spline,chi1,f_accel);
}
//integrand of Eq. 2.14 in 1807.04266
double integrand_flLp_chi(double chi1, void *params){
  double *ar = (double *) params;
  double l1 = ar[0];
  int L = (int)ar[1];
  int z1 = (int)ar[2];
  int z2 = (int)ar[3];
  double p = ar[4];
  double x = p*chi1;
  double res = 0.0;
  res = flLp_chi_response_deriv_spline(chi1, l1,z1,z2)*((L*L-L-x*x)*gsl_sf_bessel_jl(L,x)+2.*x*gsl_sf_bessel_jl(L+1,x))/(x*x);
  res += flLp_chi_response_spline(chi1, l1,z1,z2)*gsl_sf_bessel_jl(L,x);
  return res;
}

double flL_p(double p,double l1, int L, int z1, int z2){
  int zmin;
  zmin = (z2 < z1 ? z2 : z1);
  double chi_max = chi(amin_source(zmin));
  double ar[5] = {l1, (double)L, (double)z1, (double)z2,p};
  double chi1,f,delta_f;
  int N_zero = 2;
  chi1 = 0.;
  f = int_gsl_integrate_medium_precision(integrand_flLp_chi,(void*)ar,chi1,chi_max,NULL,1000);
  return f;
}
//build lookup table of flL_p, interpolate in p
//first index of table combines (l1,L,z1,z2) into unique integer
/*only compute table elements when (z1,z2) is first called
  to avoid unnecessary overhead for parallel compuations */
double flL_p_tabulated(double p, int n_l1, int L, int z1, int z2,double *ell){
  static int N1 = -1;
  static double **table_f = 0;
  static double logpmin = 0., logpmax = 0., dp = 0.;
  static int table_Np = 0;
  static cosmopara C;
  int i;
  double logp;
  if (N1 < 0){
    N1 = like.Ncl*L_MAX*tomo.shear_Npowerspectra;
    table_Np = 100;
    table_f = create_double_matrix(0, N1, 0, table_Np-1);
    for (int I = 0; I< N1; I++){
      table_f[I][0] = -123.0;
    }
    logpmin = log(limits.k_min_cH0*10.); logpmax = log(limits.k_max_cH0/100.);  //byalex: note scalings of this extrema to improve accuracy of interpolation
    dp = (logpmax - logpmin)/(table_Np-1.);
  }
  i =N_shear(z1,z2)*like.Ncl*L_MAX + n_l1*L_MAX + L;
  if (table_f[i][0]< -100.){
    for (int il1 =0; il1 < like.Ncl; il1 ++){
      for (int iL =0; iL < L_MAX; iL++){
        int I =N_shear(z1,z2)*like.Ncl*L_MAX + il1*L_MAX + iL;
        logp = logpmin;
        for (int j = 0; j < table_Np; j ++, logp += dp){
          table_f[I][j] = flL_p(exp(logp),ell[il1], iL, z1,z2);
        }
      }
    }
  }
  logp = log(p);
  if (n_l1 >= like.Ncl){printf("flL_p_tabulated: argument n_l1 = %d does not match bining like.Ncl = %d\nEXIT\n", n_l1,like.Ncl); exit(1);}
  if (L >=L_MAX){printf("flL_p_tabulated: argument L = %d does not match table bining L_MAX = %d\nEXIT\n",L, L_MAX); exit(1);}
  if (logp <= logpmin || logp >= logpmax){return 0.;}
  return interpol(table_f[i], table_Np, logpmin, logpmax, dp,logp, 0.0,0.0);
}

double integrand_sigma_L_p(double logp, void *params){
  static double *ell = 0;
  if (!ell){
    double logdl = (log(like.lmax)-log(like.lmin))/(like.Ncl-1);
    ell=create_double_vector(0,like.Ncl-1);
    for(int i=0;i<like.Ncl;i++){ ell[i]=exp(log(like.lmin)+i*logdl);}
  }
  int *ar = (int *) params;
  int L = ar[0];
  int l1 = ar[1]; int l2 = ar[2];
  int z1 = ar[3]; int z2 = ar[4]; int z3 = ar[5]; int z4 = ar[6];
  double p = exp(logp);
  return p*p*p*p_lin(p,1.0)*flL_p_tabulated(p,l1,L,z1,z2,ell)*flL_p_tabulated(p,l2,L,z3,z4,ell);
}
double integrand_sigma_L_Limber(double chi1, void *params){
  static double *ell = 0;
  if (!ell){
    double logdl = (log(like.lmax)-log(like.lmin))/(like.Ncl-1.0);
    ell=create_double_vector(0,like.Ncl-1);
    for(int i=0;i<like.Ncl;i++){ ell[i]=exp(log(like.lmin)+i*logdl);}
  }
  int *ar = (int *) params;
  int L = ar[0];
  int l1 = ar[1]; int l2 = ar[2];
  int z1 = ar[3]; int z2 = ar[4]; int z3 = ar[5]; int z4 = ar[6];
  double aa = a_chi(chi1);
  return flLp_chi_response(chi1, ell[l1], z1,z2)*flLp_chi_response(chi1, ell[l2], z3,z4)/(chi1*chi1)*p_lin((L+0.5)/chi1,1.0);
}

/********************* pieces for Limber SSC calculation  *********************/
/******************************************************************************/
/******************************************************************************/
double weights_cov_shear_shear(double a, int z1, int z2,int z3, int z4){
  double fK=f_K(chi(a));
  return W_kappa(a,fK,z1)*W_kappa(a,fK,z2)*W_kappa(a,fK,z3)*W_kappa(a,fK,z4)*dchi_da(a);

}

//  integrand for Fourier space convergence SSC covariance at (l1,l2,z1,z2,z3,z4)
//  (index(l1),index(l2),z1,z2,z3,z4) passed in params
double integrand_cov_SSC_Limber_LMAX(double a, void *params){
  static double *ell = 0;
  if (!ell){
    double logdl = (log(like.lmax)-log(like.lmin))/(like.Ncl-1.0);
    ell=create_double_vector(0,like.Ncl-1);
    for(int i=0;i<like.Ncl;i++){ ell[i]=exp(log(like.lmin)+i*logdl);}
  }
  int *ar = (int *) params;
  int l1 = ar[0]; int l2 = ar[1];
  int z1 = ar[2]; int z2 = ar[3]; int z3 = ar[4]; int z4 = ar[5];
  double fK,k1,k2,z,weights,R_perp1,R_perp2,var;
  fK = f_K(chi(a));
  k1 = (ell[l1]+0.5)/fK;
  k2 = (ell[l2]+0.5)/fK;
  z     = 1./a - 1.;
  weights = weights_cov_shear_shear(a,z1,z2,z3,z4);
  if (weights >0){
    R_perp1 = R_1_coef(k1,a)+R_K_coef(k1,a)/6.;
    R_perp2 = R_1_coef(k2,a)+R_K_coef(k2,a)/6.;
    var = sum_variance_healpix_LMIN(a, L_MAX);
    return R_perp1*R_perp2*Pdelta(k1,a)*Pdelta(k2,a)*weights*var*pow(fK,-4.);
  }
  return 0.;
}
//Limber integral for SSC covariance, restricting the variance calculation to p > (L_MAX+1)/chi
double cov_SSC_Limber_LMAX(int l1, int l2, int z1, int z2, int z3, int z4){
  int ar[6] ={l1,l2,z1,z2,z3,z4};
  return int_gsl_integrate_low_precision(integrand_cov_SSC_Limber_LMAX,(void*)ar,amin_source(z1),amax_source(z1),NULL,1000);
}

double sigma_Ll1l2_Limber(int L, int l1, int l2, int z1, int z2, int z3, int z4){
  int ar[7] ={L,l1,l2,z1,z2,z3,z4};
  return int_gsl_integrate_low_precision(integrand_sigma_L_Limber,(void*)ar,0.0,chi(amin_source(z1)),NULL,1000);
}
double sigma_Ll1l2(int L, int l1, int l2, int z1, int z2, int z3, int z4){
  int ar[7] ={L,l1,l2,z1,z2,z3,z4};
  if (L >= L_MAX){return sigma_Ll1l2_Limber(L,l1,l2,z1,z2,z3,z4);}
  return 2./M_PI*int_gsl_integrate_low_precision(integrand_sigma_L_p,(void*)ar,log(limits.k_min_cH0),log(3.e+3),NULL,1000);
}


/********************* Limber integrals for NG covariance  *********************/
/******************************************************************************/
/******************************************************************************/
//  lookup table for cNG_response(k1=(l1+0.5)/chi(a),k2=(l2+0.5)/chi(a),a)
//  retabulated when l1 or l2 change
double cNG_reponse_bin(double a, double l1, double l2, int nz)
{
  static double *table_cNG;
  static double L1=-1., L2=-1.;
  static int NZ=-1;
  static double da = 0., amin = 0., amax =0.;
  if (L1 != l1 || L2 !=l2 || NZ != nz){
    L1 =l1; L2 = l2; NZ =nz;
    if (table_cNG==0){
      table_cNG  = create_double_vector(0, Ntable.N_a_halo-1);
       amin = 1./(1. + tomo.shear_zmax[NZ]);
      amax = .999;
      da = (amax - amin)/(Ntable.N_a_halo-1.);
    }
    double aa= amin;
    for (int i=0; i<Ntable.N_a_halo; i++, aa += da) {
      table_cNG[i]=log(covNG_response((l1+0.5)/chi(aa),(l2+0.5)/chi(aa),aa));
      }
    }
  return exp(interpol(table_cNG, Ntable.N_a_halo, amin, amax, da,a, 1.0,1.0 ));
}
//  integrand for Fourier space convergence cNG covariance at (l1,l2,z1,z2,z3,z4)
//  (l1,l2,z1,z2,z3,z4) passed in params
double inner_project_tri_cov_shear_shear_tomo(double a,void *params)
{
  double fK,weights,res = 0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  weights = weights_cov_shear_shear(a,ar[2],ar[3],ar[4],ar[5]);

  if (weights >0.){
    res = cNG_reponse_bin(a,ar[0],ar[1],(int)ar[2])*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
  }
  res *= weights;
  return res;
}

/********************* final covariance expressions ***************************/
/******************************************************************************/
/******************************************************************************/
//  Fourier space convergence cNG covariance at (l1,l2,z1,z2,z3,z4)
double cov_NG_response_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
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
  return int_gsl_integrate_low_precision(inner_project_tri_cov_shear_shear_tomo,(void*)array,a1,a2,NULL,1000);
}

//  Fourier space convergence SSC covariance at (l1,l2,z1,z2,z3,z4)
//  pieces together the non-Limber summation for L < L_MAX
//  with the restructured Limber integral cov_SSC_Limber_LMAX
double cov_SSC_response_shear_shear(int l1,int l2, int z1, int z2, int z3, int z4){
  double cov_SSC = 0.;
  for (int L = 0; L < L_MAX; L ++){
    cov_SSC += C_survey_window(L)*(2*L+1)*sigma_Ll1l2(L,l1,l2,z1,z2,z3,z4)/(4.*M_PI);
  }
  cov_SSC +=cov_SSC_Limber_LMAX(l1,l2,z1,z2,z3,z4);
  return cov_SSC;
}


/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
// Note the look-up tables for power spectrum covariances are recomputed if one the redshift bins changes
//      Loop over theta1, theta2 first, before computing the next combination of redshift bins
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/**************** look-up tables for covariance  *********************/
//look-up table in (l1,l2) of Fourier space convergence cNG covariance
//recomputes when one of the tomography indices changes
double bin_cov_NG_response_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 25;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(2.);
    logsmax = log(5.1e+4);
    ds = (logsmax - logsmin)/(Ntab -1.0);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=log(cov_NG_response_shear_shear_tomo(ll1,ll2,z1,z2,z3,z4));
      }
    }
    printf("cNG tabulated\n");
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 >= logsmin && llog2 >= logsmin && llog1 <= logsmax && llog2 <= logsmax){
    res = exp(interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,1.0,1.0));}
  return res;
}

//look-up table in (l1,l2) of Fourier space convergence SSC covariance
//recomputes when one of the tomography indices changes
double bin_cov_SSC_response_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 25;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    //these are required to tabulate l-binning in subroutines of cov_SSC_shear_shear
    like.Ncl = Ntab;
    like.lmin = 2.;
    like.lmax = 5.1e+4;
    logsmin = log(like.lmin);
    logsmax = log(like.lmax);
    ds = (logsmax - logsmin)/(Ntab-1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=log(cov_SSC_response_shear_shear(i,j,z1,z2,z3,z4));
      }
    }
    printf("SSC tabulated\n");
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 >= logsmin && llog2 >= logsmin && llog1 <= logsmax && llog2 <= logsmax){
    res = exp(interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,1.0,1.0));}
  return res;
}
