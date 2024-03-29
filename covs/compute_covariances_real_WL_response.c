/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by CosmoLike developers
******************************************************************************/

#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <dirent.h>
#include <sys/types.h>
#include <unistd.h>
#include <python2.7/Python.h>


#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <fftw3.h>

#include "../cosmolike_core/theory/basics.c"
#include "../cosmolike_core/theory/structs.c"
#include "../cosmolike_core/theory/recompute.c"
#include "../cosmolike_core/theory/cosmo3D.c"
#include "../cosmolike_core/theory/redshift_spline.c"
#include "../cosmolike_core/theory/halo.c"
#include "../cosmolike_core/theory/cosmo2D_fourier.c"

#include "../cosmolike_core/theory/covariances_3D.c"

#include "init.c"
#include "../cosmolike_core/theory/cov_response_WL.c"
#include "../cosmolike_core/theory/run_covariances_real_fullsky_WL_response.c"

// used in checking existence of output directory
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
//

int main(int argc, char** argv)
{
  int hit=atoi(argv[1]);
  char *inifile;
  inifile = argv[2];

  FILE *F1,*F2;
  int i,l,m,n,o,s,p,output;
  double ktmp;
  char OUTFILE[400],filename[400];

  // set this to one to output details about inputs for diagnostics
  output = 0;
  FILE *F;

  set_cosmological_parameters(inifile,1);
  set_survey_parameters(inifile,1);
  set_cov_parameters(inifile,1);

  init_source_sample(redshift.shear_REDSHIFT_FILE,tomo.shear_Nbin);
  //here: setting values internally

  int NG, cNG;
  if (covparams.ng==1){
    NG = 1;
    if (covparams.cng==1){
      cNG = 1;
    }
    else{
      cNG = 0;
    }
  }
  else {
    NG = 0;
    cNG = 0;
  }

  printf("running multi_covariance_real with NG = %d, cNG = %d\n",NG, cNG);

  like.Ntheta=covparams.ntheta;
  like.vtmin = covparams.tmin;
  like.vtmax = covparams.tmax;

  int Ntheta = like.Ntheta;
  double *thetamin,*dtheta,*theta;
  thetamin=create_double_vector(0,Ntheta);
  dtheta=create_double_vector(0,Ntheta-1);
  theta=create_double_vector(0,Ntheta-1);
  set_angular_binning(thetamin,dtheta);

  printf("numbers of power spectra:shear %d, clustering %d, ggl %d \n",tomo.shear_Npowerspectra, tomo.clustering_Npowerspectra,tomo.ggl_Npowerspectra);

  // if outdir doesn't exist, create one
  struct stat st = {0};
  if (stat(covparams.outdir, &st) == -1) {
      mkdir(covparams.outdir, 0700);
  }

  int k=1;
  if (strcmp(covparams.ss,"true")==0)
  {
    sprintf(OUTFILE,"%s_ssss_++_cov_Ntheta%d_Ntomo%d",covparams.filename,like.Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Npowerspectra; l++){
      for (m=l;m<tomo.shear_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_shear_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,like.Ntheta,l,m,1,1,k);
        }
        k=k+1;
      }
    }
    sprintf(OUTFILE,"%s_ssss_--_cov_Ntheta%d_Ntomo%d",covparams.filename,like.Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Npowerspectra; l++){
      for (m=l;m<tomo.shear_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_shear_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,like.Ntheta,l,m,0,0,k);
        }
        k=k+1;
      }
    }
    sprintf(OUTFILE,"%s_ssss_+-_cov_Ntheta%d_Ntomo%d",covparams.filename,like.Ntheta,tomo.shear_Nbin);
    for (l=0;l<tomo.shear_Npowerspectra; l++){
      for (m=0;m<tomo.shear_Npowerspectra; m++){
        if(k==hit){
          sprintf(filename,"%s%s_%d",covparams.outdir,OUTFILE,k);
          // if (fopen(filename, "r") != NULL){exit(1);}
          run_cov_shear_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,like.Ntheta,l,m,1,0,k);
        }
        k=k+1;
      }
    }
  }
  if (hit==0)
  {
    sprintf(OUTFILE,"%s",covparams.filename);
    // write_gglensing_zbins(OUTFILE);

    sprintf(OUTFILE,"%s%s.blocks",covparams.outdir,covparams.filename);
    F1 = fopen(OUTFILE,"w");
    fprintf(F1,"%d\n",k-1);
    fclose(F1);
  }

  printf("number of cov blocks for parallelization: %d\n",k-1);
  printf("-----------------\n");
  printf("PROGRAM EXECUTED\n");
  printf("-----------------\n");
  return 0;
}
