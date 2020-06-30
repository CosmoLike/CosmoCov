/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by CosmoLike developers
******************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
// #include "../class/include/class.h"

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

//void omega_a(double aa,double *om_m,double *om_v);
double omv_vareos(double a);
static inline double hoverh0(double a);
double growfac(double a);
//int func_for_growfac(double a,const double y[],double f[],void *params);
double Tsqr_EH_wiggle(double khoverMPC);
//double int_for_sigma_r_sqr(double k, void * args);
double sigma_r_sqr();
//double Delta_L_wiggle(double k);
//double Delta_lin_wiggle(double k,double a);
double p_lin(double k,double a);
//double int_sig_R_knl(double logk, void *args);
//double int_neff(double lnk, void *args);
//double int_cur(double lnk, void *args);
//void nonlin_scale(double amp, double *R_NL, double *neff, double *Curv);
//double Halofit(double k, double amp, double omm, double omv,double w_z, double R_NL, double neff,double Curv, double P_delta_Lin);
//void Delta_halofit(double **table_P_NL,double logkmin, double logkmax, double dk, double da);
//double Delta_NL_Halofit(double k_NL, double a); //k in h/Mpc
double Pdelta(double k_NL,double a); //k in coverH0 units

//double int_for_chi(double a,void * args);
double f_K(double chi);
double chi(double a);
double a_chi(double chi1);
// extern void emu(double *xstar, double *ystar, double *kstar);

//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//variable Omega_v
//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double omv_vareos(double a)
{
  return(cosmology.Omega_v*exp(-3.*((cosmology.w0+cosmology.wa+1.)*log(a)+cosmology.wa*(1.-a))));
}

//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//c evolution of omega matter and omega lamda with expansion factor

void omega_a(double aa,double *om_m,double *om_v)
{
  double a2,omega_curv;
  a2=aa*aa;
  omega_curv=1.0-cosmology.Omega_m- cosmology.Omega_v;
  *om_m=cosmology.Omega_m /(cosmology.Omega_m +aa*(omv_vareos(aa) *a2 +omega_curv));
  *om_v=omv_vareos(aa)*a2*aa/(cosmology.Omega_m+aa*(a2*omv_vareos(aa) +omega_curv));
}
//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//growth factor including Dark energy parameters w0, wa
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//function for growfac (DGL)
int func_for_growfac(double a,const double y[],double f[],void *params)
{
  //double *p=(double *)params;
  if (a == 0) {
    printf("a=0 in function 'func_for_growfac'!\n");
    exit(1);
  }
  double aa=a*a;
  double omegam=cosmology.Omega_m/(aa*a);
  double omegav=omv_vareos(a);
  double hub = hoverh0(a);
  double one_plus_mg_mu = 1.;
  hub = hub*hub;
  f[0]=y[1];
  if(cosmology.MGmu != 0){
    one_plus_mg_mu += cosmology.MGmu*omegav/hub/cosmology.Omega_v;
  }
  f[1]=y[0]*3.*cosmology.Omega_m/(2.*hub*aa*aa*a)*one_plus_mg_mu-y[1]/a*(2.-(omegam+(3.*(cosmology.w0+cosmology.wa*(1.-a))+1)*omegav)/(2.*hub));
  return GSL_SUCCESS;
}

static inline double hoverh0(double a){
  return sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
}

double growfac(double a)
{
  const double MINA=1.e-8;
  static cosmopara C;
  static double *ai;
  static double *table;
  double res;
  
  gsl_interp *intf=gsl_interp_alloc(gsl_interp_linear,Ntable.N_a);
  gsl_interp_accel *acc=gsl_interp_accel_alloc();
  
  if (recompute_expansion(C))
  {

    if(table!=0) free_double_vector(table,0, Ntable.N_a-1);
    if(ai!=0) free_double_vector(ai,0, Ntable.N_a-1);    
    ai=create_double_vector(0, Ntable.N_a-1);
    table=create_double_vector(0, Ntable.N_a-1);
    
    int i;

    const gsl_odeiv_step_type *T=gsl_odeiv_step_rkf45;
    gsl_odeiv_step *s=gsl_odeiv_step_alloc(T,2);
    gsl_odeiv_control *c=gsl_odeiv_control_y_new(1.e-6,0.0);
    gsl_odeiv_evolve *e=gsl_odeiv_evolve_alloc(2);
    
    double t=MINA;            //start a
    double t1=1.1;                //final a
    double h=1.e-6;              //initial step size
    double y[2]={MINA,MINA};   //initial conditions
    double norm;
    double par[0]={};
    gsl_odeiv_system sys={func_for_growfac,NULL,2,&par};
    
    for (i=1;i<=Ntable.N_a;i++) {
      ai[i-1]=i*t1/(1.*Ntable.N_a);
      while(t<ai[i-1]) 
        gsl_odeiv_evolve_apply(e,c,s,&sys,&t,ai[i-1],&h,y);
      if (i==1) norm=y[0]/ai[i-1];
      table[i-1]=y[0]/norm;
    }
    
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
    update_cosmopara(&C);
	 
  }
  gsl_interp_init(intf,ai,table,Ntable.N_a);
  res=gsl_interp_eval(intf,ai,table,a,acc);
  gsl_interp_accel_free(acc);
  gsl_interp_free(intf);
  return(res);
}



// ---------------------------- Transfer Function from EH98 ---------------------- 
//Input: k -- Wavenumber at which to calculate transfer function, in Mpc^-1. Output: Returns the value of the full transfer function fitting formula. This is the form given in Section 3 of Eisenstein & Hu (1997). Notes: Units are Mpc, not h^-1 Mpc. 

double Tsqr_EH_wiggle(double khoverMPC)
{
  static double omhh=-123.;
  static double obhh=-123.;
  static double OMEGA_V = -123.;
  static double f_baryon;
  
  static double k_equality;
  static double sound_horizon;
  static double beta_c;
  static double alpha_c;
  static double beta_node;
  static double alpha_b;
  static double beta_b;
  static double k_silk;
  
  //if (omhh != cosmology.Omega_m*cosmology.h0*cosmology.h0 || obhh != cosmology.omb*cosmology.h0*cosmology.h0|| OMEGA_V != cosmology.Omega_v){  
  double theta_cmb,z_equality,z_drag,R_drag,R_equality;
  double z_drag_b1, z_drag_b2;
  double alpha_c_a1, alpha_c_a2, beta_c_b1, beta_c_b2, alpha_b_G, y;

  omhh = cosmology.Omega_m*cosmology.h0*cosmology.h0;  
  obhh = cosmology.omb*cosmology.h0*cosmology.h0;
  OMEGA_V = cosmology.Omega_v;
  f_baryon = obhh/omhh; 
    //printf("%le\n",f_baryon);

    theta_cmb=2.728/2.7;// Tcmb in units of 2.7 K 
    z_equality= 2.50e4*omhh/POW4(theta_cmb);//Redshift of matter-radiation equality, really 1+z 
    
    k_equality=0.0746*omhh/SQR(theta_cmb);//Scale of equality, in Mpc^-1 
    z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
    z_drag_b2 = 0.238*pow(omhh,0.223);
    z_drag=1291*pow(omhh,0.251)/(1+0.659*pow(omhh,0.828))*(1+z_drag_b1*pow(obhh,z_drag_b2));//Redshift of drag epoch 
    R_drag=31.5*obhh/POW4(theta_cmb)*(1000/(1+z_drag));//Photon-baryon ratio at drag epoch 
    
    R_equality= 31.5*obhh/POW4(theta_cmb)*(1000/z_equality);//Photon-baryon ratio at equality epoch 
    sound_horizon=2./3./k_equality*sqrt(6./R_equality)*log((sqrt(1+R_drag)+sqrt(R_drag+R_equality))/(1+sqrt(R_equality)));//Sound horizon at drag epoch, in Mpc 
    k_silk= 1.6*pow(obhh,0.52)*pow(omhh,0.73)*(1+pow(10.4*omhh,-0.95));//Silk damping scale, in Mpc^-1 
    alpha_c_a1 = pow(46.9*omhh,0.670)*(1+pow(32.1*omhh,-0.532));
    alpha_c_a2 = pow(12.0*omhh,0.424)*(1+pow(45.0*omhh,-0.582));
    alpha_c=pow(alpha_c_a1,-f_baryon)*pow(alpha_c_a2,-CUBE(f_baryon)); //CDM suppression 
    
    beta_c_b1 = 0.944/(1+pow(458*omhh,-0.708));
    beta_c_b2 = pow(0.395*omhh, -0.0266);
    beta_c=1.0/(1+beta_c_b1*(pow(1-f_baryon, beta_c_b2)-1));//CDM log shift 
    
    y = z_equality/(1+z_drag);
    alpha_b_G = y*(-6.*sqrt(1+y)+(2.+3.*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)));
    alpha_b=2.07*k_equality*sound_horizon*pow(1+R_drag,-0.75)*alpha_b_G;//Baryon suppression 
    beta_node = 8.41*pow(omhh, 0.435);//Sound horizon shift 
    beta_b= 0.5+f_baryon+(3.-2.*f_baryon)*sqrt(pow(17.2*omhh,2.0)+1);//Baryon envelope shift 
  //}  
  // Start of TFfit_onek from the original tf_fit.c at http://background.uchicago.edu/~whu/transfer/transferpage.html
    double T_c_ln_beta,T_c_ln_nobeta , T_c_C_alpha, T_c_C_noalpha;
    double q,qsqr, xx, xx_tilde;
    double T_c_f, T_c, s_tilde, T_b_T0, T_b;

  double k=khoverMPC*cosmology.h0; //internally this routine uses Mpc^-1 not h/Mpc
  q = k/13.41/k_equality; 
  qsqr =SQR(q);
  xx = k*sound_horizon;
  
  T_c_ln_beta = log(2.718282+1.8*beta_c*q);
  T_c_ln_nobeta = log(2.718282+1.8*q);
  
  T_c_ln_beta = log(2.718282+1.8*beta_c*q);
  T_c_C_alpha = 14.2/alpha_c + 386.0/(1+69.9*pow(q,1.08));
  T_c_C_noalpha = 14.2 + 386.0/(1+69.9*pow(q,1.08));
  
  T_c_f = 1.0/(1.0+POW4(xx/5.4));
  T_c = T_c_f*T_c_ln_beta/(T_c_ln_beta+T_c_C_noalpha*qsqr) +(1-T_c_f)*T_c_ln_beta/(T_c_ln_beta+T_c_C_alpha*qsqr);
  
  s_tilde = sound_horizon*pow(1+CUBE(beta_node/xx),-1./3.);
  xx_tilde = k*s_tilde;
  
  T_b_T0 = T_c_ln_nobeta/(T_c_ln_nobeta+T_c_C_noalpha*qsqr);
  
  T_b = sin(xx_tilde)/(xx_tilde)*(T_b_T0/(1+SQR(xx/5.2))+alpha_b/(1+CUBE(beta_b/xx))*exp(-pow(k/k_silk,1.4)));
  
  return SQR(f_baryon*T_b + (1-f_baryon)*T_c); 
}


//Calculate Normalization see Cosmology Notes 8.105  
double int_for_sigma_r_sqr(double k, void * args)
{
  double kR, res, x;
  kR = k*8.; // r=8 Mpc/h
  x = (sin(kR) - kR*cos(kR))/(kR*kR*kR);
  res = pow(k,2.+cosmology.n_spec+ 0.5*cosmology.alpha_s*log(k/0.05))*Tsqr_EH_wiggle(k)*x*x;
  return res;
}

double sigma_r_sqr()   
{
  static double res = -123.;
  static cosmopara C;
  double integral,array[1];
  
  if (recompute_Delta(C)) //strictly speaking, this is recomputed unnecessarily if only sigma_8 changes
  {
    integral = int_gsl_integrate_medium_precision(int_for_sigma_r_sqr,(void*)array,1e-4,1e6,NULL,512);
    res = 9.0*integral;   //see Peackock97, eq. 29
    update_cosmopara(&C);

  }
  if (!(res>0.0)){
    fprintf(stderr,"failed with sigma_r_sqr = %le\n", res);
  }
  assert(res>0.0);
  return res;
}



double Delta_L_wiggle(double k)
{
  static cosmopara C;
  
  static double *table_P;
  static double dk = .0, logkmin = .0, logkmax = .0;
  
  double klog,f1,norm;
  int i;     
  
  if (k < limits.k_min_mpc || k > limits.k_max_mpc){
    norm=cosmology.sigma_8*cosmology.sigma_8/sigma_r_sqr(); 
    
    return norm*pow(k,cosmology.n_spec+ 0.5*cosmology.alpha_s*log(k/0.05)+3.0)*Tsqr_EH_wiggle(k);
    //printf("outside Delta_L_tab\n");   
  }
  else{  
    if (recompute_Delta(C))
    {
      if (cosmology.M_nu > 0){
        printf("Implementation of EH transfer function does not support massive neutrinos\n EXIT\n");
      }
      update_cosmopara(&C);
      norm=cosmology.sigma_8*cosmology.sigma_8/sigma_r_sqr();
      
      if(table_P!=0) free_double_vector(table_P,0, Ntable.N_k_lin-1);
      table_P=create_double_vector(0, Ntable.N_k_lin-1);
      
      logkmin = log(limits.k_min_mpc);
      logkmax = log(limits.k_max_mpc);
      dk = (logkmax - logkmin)/(Ntable.N_k_lin-1.);
      klog = logkmin;
      
      for (i=0; i<Ntable.N_k_lin; i++, klog += dk) {  
       table_P[i]=log(norm*pow(exp(klog),cosmology.n_spec+ 0.5*cosmology.alpha_s*log(k/0.05)+3.0)*Tsqr_EH_wiggle(exp(klog)));
     }
      //printf("finished Delta_L_wiggle\n");   
   }
 }
 klog=log(k);
 f1=interpol(table_P, Ntable.N_k_lin, logkmin, logkmax, dk,klog, 1.0,1.0 );
 return exp(f1);  
}

// linear power spectrum routine with k in units H_0/c; used in covariances.c for beat coupling and in halo.c
double p_lin(double k,double a)
{
  static cosmopara C;
  static double **table_P_Lz = 0;
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  int status;
  // if (strcmp(pdeltaparams.runmode,"CLASS")==0 || strcmp(pdeltaparams.runmode,"class")==0) return p_class(k,a,0, &status);
  
  double amp,ampsqr,grow0,aa,klog,val;
  
  int i,j;
  if (a >= 0.99999){a =0.99999;}
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);
    if (table_P_Lz!=0) free_double_matrix(table_P_Lz,0, Ntable.N_a-1, 0, Ntable.N_k_lin-1);
    table_P_Lz = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_lin-1);
    grow0=growfac(1.);
    da = (1. - limits.a_min)/(Ntable.N_a-1.);
    aa = limits.a_min;
    for (i=0; i<Ntable.N_a; i++, aa +=da) {
      if(aa>1.0) aa=1.0;
      amp=growfac(aa)/grow0;
      ampsqr=amp*amp;
      
      logkmin = log(limits.k_min_mpc);
      logkmax = log(limits.k_max_mpc);
      dk = (logkmax - logkmin)/(Ntable.N_k_lin-1.);
      klog = logkmin;
      for (j=0; j<Ntable.N_k_lin; j++, klog += dk) {
        table_P_Lz[i][j] = log(ampsqr*Delta_L_wiggle(exp(klog)));
        //printf("%le %le",exp(klog),Delta_L_wiggle(exp(klog));
      }
    }
  }
  if (k/cosmology.coverH0 > exp(logkmax)) return 0.0;
  klog = log(k/cosmology.coverH0);
  val = interpol2d(table_P_Lz, Ntable.N_a, limits.a_min, 1., da, a, Ntable.N_k_lin, logkmin, logkmax, dk, klog, 3.0+cosmology.n_spec, 0.0);
  if(isnan(val) || (k==0)) return 0.0;
  return 2.0*constants.pi_sqr*exp(val)/k/k/k;     
}


double int_sig_R_knl(double lnk, void *args) // tak12 A4
{
  double krsqr;
  
  double *params= (double *) args;
  double Rscale=params[0];
  //printf("Rscale %le k %le\n",Rscale,exp(lnk));
  krsqr= SQR(exp(lnk)*Rscale);
  return Delta_L_wiggle(exp(lnk))*exp(-krsqr);
}


double int_neff(double lnk, void *args) //tak12 A5
{
  double krsqr;
  double *params= (double *) args;
  double Rscale=params[0];
  krsqr= SQR(exp(lnk)*Rscale);
  return Delta_L_wiggle(exp(lnk))*2.0*krsqr*exp(-krsqr); //see S03 eq. 59
}


double int_cur(double lnk, void *args) //tak12 A5
{
  double krsqr;
  double *params= (double *) args;
  double Rscale=params[0];
  krsqr= SQR(exp(lnk)*Rscale);
  return Delta_L_wiggle(exp(lnk))*4.0*krsqr*(1.0-krsqr)*exp(-krsqr); // S03 eq.60
}


//iterative calculation of the nonlinear scale as defined in tak12 A4
void nonlin_scale(double amp, double *R_NL, double *neff, double *Curv)
{
  double sig_R,kmax,logkmax,sig_R_noamp,neffplus3;
  int iterstep; 
  const int itermax  = 40;
  int converged=0;
  double array[1];
  double logRmin = -3.0;
  double logRmax =  4.0; 

  iterstep=0;  
  while(converged==0)
  {        
    array[0]=pow(10.,(logRmin+logRmax)/2.0);
    
    //flexible upper limit of integration depending on drop-off of filter function
    kmax  = sqrt(5.*log(10.))/array[0];
    if (kmax<8000.0) logkmax = log(8000.0);
    sig_R_noamp=sqrt(int_gsl_integrate_medium_precision(int_sig_R_knl,(void*)array,-4.5,logkmax,NULL,512)); //integral goes over ln k exponent correspond to k_min~0.011
    
    sig_R=amp*sig_R_noamp; 
    if (sig_R>1.0)  logRmin=log10(array[0]);
    if (sig_R<1.0)  logRmax=log10(array[0]);
    iterstep=iterstep+1;
    if(fabs(sig_R-1.0) < 0.0001 || iterstep>itermax) converged=1;    
  }  
  *R_NL=array[0]; //R where sig_R==1
  neffplus3=int_gsl_integrate_medium_precision(int_neff,(void*)array,-4.5,logkmax,NULL,512)/sig_R_noamp/sig_R_noamp;
  *neff= neffplus3 - 3.0;
  *Curv= int_gsl_integrate_medium_precision(int_cur,(void*)array,-4.5,logkmax,NULL,512)/sig_R_noamp/sig_R_noamp + SQR(neffplus3);
  
  //printf("%d %le\n",iterstep,amp);  
}                  


double Halofit(double k, double amp, double omm, double omv,double w_z, double R_NL, double neff,double Curv, double P_delta_Lin)
{
  double y_scale,n2eff,n3eff,n4eff;
  double a_n,b_n,c_n,gamma_n,alpha_n,beta_n,nu_n,f1,f2,f3;
  
  double Delta_H,Delta_H_Prime,Delta_Q;
  //determine nonlinear scale, neff and curvature, see tak12 A4, A5
  y_scale=k*R_NL;
  
  n2eff=neff*neff;
  n3eff=n2eff*neff;
  n4eff=n2eff*n2eff;
  
  //calculate coefficients 
  a_n = pow(10.,1.5222+2.8553*neff + 2.3706*n2eff+0.9903*n3eff+0.2250*n4eff-0.6038*Curv+0.1749*omv*(1.0+w_z));
  b_n = pow(10., -0.5642+0.5864*neff + 0.5716*n2eff-1.5474*Curv +0.2279*omv*(1.0+w_z));
  c_n = pow(10., 0.3698+ 2.0404*neff + 0.8161*n2eff+0.5869*Curv);
  gamma_n = 0.1971-0.0843*neff + 0.8460*Curv;
  alpha_n = fabs(6.0835 + 1.3373*neff - 0.1959*n2eff - 5.5274*Curv);
  beta_n = 2.0379 - 0.7354*neff + 0.3157*n2eff + 1.2490*n3eff + 0.3980*n4eff - 0.1682*Curv;
  nu_n = pow(10,5.2105+3.6902*neff);
  
  f1 = pow(omm,(-0.0307));
  f2 = pow(omm,(-0.0585));
  f3 = pow(omm,(0.0743));  
  
  //TwoHaloTerm
  Delta_Q=P_delta_Lin*(pow((1.0+P_delta_Lin),beta_n)/(1.0+alpha_n*P_delta_Lin))*exp(-(y_scale/4.0+y_scale*y_scale/8.0));
  //OneHaloterm
  Delta_H_Prime=(a_n*pow(y_scale,3.0*f1))/(1.0+b_n*pow(y_scale,f2)+pow(c_n*f3*y_scale,3.0-gamma_n));
  Delta_H=Delta_H_Prime/(1.0+nu_n*pow(y_scale,-2.0)); // using mu=0.0 Tak A12
  //printf("Delta_Q %le Delta_H %le\n",Delta_Q,Delta_H);
  return Delta_H+Delta_Q;
}




void Delta_halofit(double **table_P_NL,double logkmin, double logkmax, double dk, double da)
{       
  double rk,omm,omv,w_z,amp,grow0,aa,klog;
  double R_NL,Curv,neff,P_delta,P_delta_Lin;
  int i,j;
  
  grow0=growfac(1.);
  aa = limits.a_min;
  //binning in k and a must be the same as in emu
  for (i=0; i<Ntable.N_a; i++, aa +=da) { 
    if(aa>1.0) aa=1.0;
    omega_a(aa,&omm,&omv);
    w_z=cosmology.w0+cosmology.wa*(1.-aa);
    amp=growfac(aa)/grow0;
    nonlin_scale(amp, &R_NL, &neff, &Curv);
    //printf("%le %le %le %le\n",aa,R_NL,neff,Curv);
    klog = logkmin;
    for (j=0; j<Ntable.N_k_nlin; j++, klog += dk) { 
      rk=exp(klog);
      P_delta_Lin=amp*amp*Delta_L_wiggle(rk);
      P_delta=Halofit(rk, amp, omm, omv, w_z, R_NL, neff, Curv, P_delta_Lin);
      table_P_NL[i][j]=log(P_delta);
    }
  }
}


double Delta_NL_Halofit(double k_NL, double a)
{     
  static cosmopara C;
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  
  static double **table_P_NL=0;
  double klog,val; 

  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);
    if (table_P_NL!=0) free_double_matrix(table_P_NL,0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);
    table_P_NL = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);     
    
    da = (1. - limits.a_min)/(Ntable.N_a-1.);
    logkmin = log(limits.k_min_mpc);
    logkmax = log(limits.k_max_mpc);
    dk = (logkmax - logkmin)/(Ntable.N_k_nlin-1.);
    
    Delta_halofit(table_P_NL,logkmin, logkmax, dk, da);
  }
  klog = log(k_NL);
  val = interpol2d(table_P_NL, Ntable.N_a, limits.a_min, 1., da, a, Ntable.N_k_nlin, logkmin, logkmax, dk, klog, cosmology.n_spec, 0.0);
  return exp(val);
  // returns the dimensionless power spectrum as a function of scale factor a and k  
}


double nonlinear_scale_computation(double a)
{       
  static cosmopara C;
  static double da = 0.;
  static double *table=0;
  
  double omm,omv,amp,grow0,aa,res;
  double R_NL,Curv,neff;
  int i;
  
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);
    grow0=growfac(1.);
    da = (1.-limits.a_min)/(Ntable.N_a-1.);
    aa = limits.a_min;
    if (table!=0) free_double_vector(table, 0, Ntable.N_a-1);
    table   = create_double_vector(0, Ntable.N_a-1);
    for (i=0; i<Ntable.N_a; i++, aa+=da) {
      if(aa>1.0) aa=1.0;
      omega_a(aa,&omm,&omv);
      amp=growfac(aa)/grow0;
      nonlin_scale(amp, &R_NL, &neff, &Curv);
      table[i] = 1./R_NL;
    }
  }
  res = interpol(table, Ntable.N_a, limits.a_min, 1., da, a, 0.0, 0.0); 
  return res;
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
//Pdelta is called with k in units H0/c since the comoving distance chi is in units c/H0. Upstream Pdelta all routines are in h/mpc 
double Pdelta(double k_NL,double a)
{ 
  static int P_type = -1;
  if (P_type == -1){
    if (strcmp(pdeltaparams.runmode,"Halofit")==0) P_type = 0;
    if (strcmp(pdeltaparams.runmode,"halofit")==0) P_type = 0;
    // if (strcmp(pdeltaparams.runmode,"emu")==0) P_type = 1;
    // if (strcmp(pdeltaparams.runmode,"emu_only")==0) P_type = 2;
    if (strcmp(pdeltaparams.runmode,"linear")==0) P_type = 3;
    // if (strcmp(pdeltaparams.runmode,"CLASS")==0) P_type = 4;
    // if (strcmp(pdeltaparams.runmode,"class")==0) P_type = 4;
    if (strcmp(pdeltaparams.runmode,"cosmo_sim_test") ==0) P_type = 5;

  }

  //printf("%s set\n",pdeltaparams.runmode);
  double pdelta = 0.,kintern=k_NL/cosmology.coverH0,error,k_nonlin,res;
  int status;
  switch (P_type){
    case 0: pdelta=2.0*constants.pi_sqr*Delta_NL_Halofit(kintern,a)/k_NL/k_NL/k_NL; break;
    // case 1: pdelta=2.0*constants.pi_sqr*Delta_NL_emu(kintern,a)/k_NL/k_NL/k_NL; break;
    // case 2: pdelta=2.0*constants.pi_sqr*Delta_NL_emu_only(kintern,a)/k_NL/k_NL/k_NL; break;
    case 3: pdelta=p_lin(k_NL,a); break;
    // case 4: pdelta=p_class(k_NL,a,1, &status); break;
    case 5: k_nonlin=nonlinear_scale_computation(a);
    if (kintern<0.01) pdelta=2.0*constants.pi_sqr*Delta_NL_Halofit(kintern,a)/k_NL/k_NL/k_NL;
    else{ 
      error=0.01*pow((pdeltaparams.DIFF_A*kintern/k_nonlin),pdeltaparams.DIFF_n);
      pdelta=2.0*constants.pi_sqr*Delta_NL_Halofit(kintern,a)*(1.0+error)/k_NL/k_NL/k_NL;
    }
    break;
    default: 
    printf("cosmo3D:Pdelta: %s Pdelta runmode not defined\n",pdeltaparams.runmode);
    printf("using Halofit (standard)\n");
    pdelta=2.0*constants.pi_sqr*Delta_NL_Halofit(kintern,a)/k_NL/k_NL/k_NL;
    break;
  }

//      double z = 1./a -1 ; 
//	  if (z > 4.) { 
//		 printf("	z:%lf    a:%lf   	k_Mpc:%lf    PkRatio:%lf \n",z,a,kintern,PkRatio_baryons(kintern, a));
//	  }
  
  // if (bary.isPkbary==1) pdelta = pdelta*PkRatio_baryons(kintern, a);
  return pdelta;  
}    

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*============================================================
 *see BS 2.41 bzw Logbook for detailed calculation of chi from a.*/

double int_for_chi(double a, void * args){
  //double res,asqr;
  //asqr=a*a;
  //res= 1./sqrt(a*cosmology.Omega_m + asqr*(1.-cosmology.Omega_m -cosmology.Omega_v ) + asqr*asqr*omv_vareos(a));
  //return res;
  return 1./(a*a*hoverh0(a)); //changed to call of hoverh0 to be ready for other parametrizations
}

/*for the calculation of chi we have to integrate from a(z2)=a up to a(z1)=1, which means todays expansion factor*/	
double chi(double a)
{
  static cosmopara C;
  static double *table;
  static double da = 0.;
  double aa,res;
  int i;
  double array[1];
  
  if (recompute_expansion(C)){
    update_cosmopara(&C);
    da = (1.-limits.a_min)/(Ntable.N_a-1.);
    aa = limits.a_min;
    if (table!=0) free_double_vector(table, 0, Ntable.N_a-1);
    table   = create_double_vector(0, Ntable.N_a-1);
    for (i=0; i<Ntable.N_a-1; i++, aa+=da) {
      table[i] = int_gsl_integrate_medium_precision(int_for_chi,(void*)array, aa, 1.,NULL,1000);
    }
    table[Ntable.N_a-1] =.0;
  }
  res = interpol(table, Ntable.N_a, limits.a_min, 1., da, a, 0.0, 0.0); // comoving distance in c/H_0
  if (res < 0){printf ("interpolation error in chi(%e)\n",a); res=0.01;}
  return res;
}
//auxilary function to look up a(chi)
double a_chi(double chi1){
  static gsl_spline * a_spline = NULL;
  static gsl_interp_accel * a_accel = NULL;
  static cosmopara C;
  static double chi_max =-1.;
  if (!a_spline){
    a_spline = gsl_spline_alloc(gsl_interp_cspline, Ntable.N_a);
    a_accel = gsl_interp_accel_alloc();
  }
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);
    double *table_a,*table_chi;
    table_a  = create_double_vector(0, Ntable.N_a-1);
    table_chi  = create_double_vector(0, Ntable.N_a-1);
    for (int i = 0; i < Ntable.N_a; i++){
      table_a[i] = 1.0 - (1.0 - 0.99*limits.a_min)/(Ntable.N_a-1.)*(double)i;
      table_chi[i] = int_gsl_integrate_medium_precision(int_for_chi,NULL, table_a[i], 1.,NULL,1000);
   //   printf("%d %e %e\n",i,table_a[i],table_chi[i]);
    }
    chi_max = int_gsl_integrate_medium_precision(int_for_chi,NULL, limits.a_min, 1.,NULL,1000);
    gsl_spline_init(a_spline, table_chi, table_a, Ntable.N_a);
    free_double_vector(table_a,0, Ntable.N_a-1);
    free_double_vector(table_chi,0, Ntable.N_a-1);
  }
  if (chi1 <=0.0){return 1.0;}
  if (chi1 > chi_max){printf("called a_chi(chi) with chi > chi(limits.a_min\nEXIT\n");exit(1);}
  return gsl_spline_eval(a_spline,chi1,a_accel);
}

/*===============================calculating the angular diameter distance f_K BS01 2.4, 2.30: f_K is a radial function that, depending on the curvature of the Universe, is a trigonometric, linear, or hyperbolic function of chi  */
double f_K(double chi)
{  
  double K, K_h, f;
  K = (cosmology.Omega_m   + cosmology.Omega_v  - 1.);
  if (K > precision.medium) {           /* open */
    K_h = sqrt(K); // K in units H0/c see BS eq. 2.30
    f = 1./K_h*sin(K_h*chi);
    //printf("open\n");
  } else if (K < -precision.medium) {   /* closed */
    K_h = sqrt(-K); 
    f = 1./K_h*sinh(K_h*chi);
    //printf("closed K=%le %le %le\n",K,cosmology.Omega_m,cosmology.Omega_v);
  } else {                     /* flat */
    f = chi;
    //printf("flatK=%le %le %le\n",K,cosmology.Omega_m,cosmology.Omega_v);
  }
  return f;
}




