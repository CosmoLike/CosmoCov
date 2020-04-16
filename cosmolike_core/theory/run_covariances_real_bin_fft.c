/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by CosmoLike developers
******************************************************************************/

#include <time.h>

void run_cov_shear_shear_real_binned(char *OUTFILE, char *PATH, double *t, double *theta, double *dtheta, int Ntheta, int n1, int n2, int pm1, int pm2, int start);
void run_cov_clustering_real_binned(char *OUTFILE, char *PATH, double *t, double *theta, double *dtheta,int Ntheta, int n1, int n2, int start);
void run_cov_ggl_real_binned(char *OUTFILE, char *PATH, double *t, double *theta, double *dtheta,int Ntheta, int n1, int n2, int start);
void run_cov_ggl_shear_real_binned(char *OUTFILE, char *PATH, double *t, double *theta, double *dtheta, int Ntheta, int n1, int n2, int pm, int start);
void run_cov_clustering_shear_real_binned(char *OUTFILE, char *PATH, double *t, double *theta, double *dtheta,int Ntheta, int n1, int n2, int pm, int start);
void run_cov_clustering_ggl_real_binned(char *OUTFILE, char *PATH, double *t, double *theta, double *dtheta,int Ntheta, int n1, int n2, int start);

void print_citations(FILE *F1){
  fprintf(F1, "# Please cite the following papers in research using this covariance:\n");
  fprintf(F1, "# arXiv: 1601.05779, https://arxiv.org/abs/1601.05779\n");
  fprintf(F1, "# arXiv: 2004.04833, https://arxiv.org/abs/2004.04833\n");
  if (w_mask(like.vtmin) < 1.0){
    fprintf(F1, "# arXiv: 1804.10663, https://arxiv.org/abs/1804.10663\n");
  }
  fprintf(F1, "########################################################\n");
}

void run_cov_shear_shear_real_binned(char *OUTFILE, char *PATH, double *t, double *theta, double *dtheta, int Ntheta, int n1, int n2, int pm1, int pm2, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  z1 = Z1(n1); z2 = Z2(n1);
  printf("N_shear1 = %d (%d,%d)\n", n1,z1,z2);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear2 = %d (%d, %d)\n",n2,z3,z4);
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  printf("%s\n",filename);
  F1 =fopen(filename,"w");

  double **cov_g, **cov_ng;
  int i;
  cov_g = malloc(Ntheta*sizeof(double *));
  cov_ng = malloc(Ntheta*sizeof(double *));
  for(i=0;i<Ntheta;i++) {
    cov_g[i] = malloc(Ntheta*sizeof(double));
    cov_ng[i] = malloc(Ntheta*sizeof(double));
  }

  cov_G_shear_shear_real_fft_binned(theta, Ntheta, dtheta, z1,z2,z3,z4,pm1,pm2, cov_g);
  if(covparams.ng){cov_NG_shear_shear_real_fft_binned(theta, Ntheta, dtheta, z1,z2,z3,z4,pm1,pm2, cov_ng);}

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
    fprintf(F2, "# i, bin-averaged theta_i, 2pt func, {s: source, l:lens}tomo bin index 1, {s, l}tomo bin index 2\n");
  }

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      c_ng = 0.; c_g = 0.;
      c_g =  cov_g[nl1][nl2]; // modified here
      if(covparams.ng){ c_ng = cov_ng[nl1][nl2];} // modified here
      if(pm1==1 && pm2==1) fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*n1+nl1,Ntheta*(n2)+nl2, t[nl1],t[nl2],z1,z2,z3,z4,c_g,c_ng);
      if(pm1==0 && pm2==0) fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(tomo.shear_Npowerspectra+n2)+nl2,t[nl1],t[nl2],z1,z2,z3,z4,c_g,c_ng);
      if(pm1==1 && pm2==0) fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*n1+nl1,Ntheta*(tomo.shear_Npowerspectra+n2)+nl2,t[nl1],t[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
    if(if_write) {
      if(pm1==0){
        fprintf(F2, "%d %e xi- s%d s%d\n", i1+nl1, t[nl1], z1, z2);
      }
      else{
        fprintf(F2, "%d %e xi+ s%d s%d\n", i1+nl1, t[nl1], z1, z2);
      }
    }
  }
  fclose(F1);
  if(if_write) {fclose(F2);}
  free(cov_g);free(cov_ng);
}


void run_cov_clustering_real_binned(char *OUTFILE, char *PATH, double *t, double *theta, double *dtheta,int Ntheta, int n1, int n2, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  //  printf("\nN_cl_1 = %d \n", n1);
  z3 = n2; z4 = n2;
  printf("N_cl_1 = %d, N_cl_2 = %d\n", n1,n2);

  double **cov_g, **cov_ng;
  int i;
  cov_g = malloc(Ntheta*sizeof(double *));
  cov_ng = malloc(Ntheta*sizeof(double *));
  for(i=0;i<Ntheta;i++) {
    cov_g[i] = malloc(Ntheta*sizeof(double));
    cov_ng[i] = malloc(Ntheta*sizeof(double));
  }

  cov_G_cl_cl_real_fft_binned(theta, Ntheta, dtheta, z1,z2,z3,z4, cov_g);
  if (z1 == z3){
      if (covparams.ng){cov_NG_cl_cl_real_fft_binned(theta, Ntheta, dtheta, z1,z2,z3,z4, cov_ng);}
  }

  print_citations(F1);

  FILE *F2;
  char filename2[300];

  int if_write = 0, i1;
  i1 = Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1);
  sprintf(filename2,"%sorder_%s_i_%d-%d",PATH,covparams.filename,i1,i1+Ntheta-1);
  if (fopen(filename2, "r") == NULL){
    if_write=1;
    F2 =fopen(filename2,"w");
    fprintf(F2, "# i, bin-averaged theta_i, 2pt func, {s: source, l:lens}tomo bin index 1, {s, l}tomo bin index 2\n");
  }

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      c_ng = 0.;
      c_g = cov_g[nl1][nl2];
      if (z1 == z3){
          if (covparams.ng){c_ng = cov_ng[nl1][nl2];}
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra + n2)+nl2,t[nl1],t[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
    if(if_write) {
      fprintf(F2, "%d %e w l%d l%d\n", i1+nl1, t[nl1], z1, z2);
    }
  }
  fclose(F1);
  if(if_write) {fclose(F2);}
  free(cov_g);free(cov_ng);
}



void run_cov_ggl_real_binned(char *OUTFILE, char *PATH, double *t, double *theta, double *dtheta,int Ntheta, int n1, int n2, int start)
{
  int zl1,zl2,zs1,zs2,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");

  zl1 = ZL(n1); zs1 = ZS(n1);
  printf("\nN_tomo_1 = %d (%d, %d)\n", n1,zl1,zs1);
  zl2 = ZL(n2); zs2 = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl2,zs2);

  double **cov_g, **cov_ng;
  int i;
  cov_g = malloc(Ntheta*sizeof(double *));
  cov_ng = malloc(Ntheta*sizeof(double *));
  for(i=0;i<Ntheta;i++) {
    cov_g[i] = malloc(Ntheta*sizeof(double));
    cov_ng[i] = malloc(Ntheta*sizeof(double));
  }

  cov_G_gl_gl_real_fft_binned(theta, Ntheta, dtheta, zl1,zs1,zl2,zs2, cov_g);
  int run_ng = 0;
  if(covparams.ng && zl1 == zl2 && test_zoverlap(zl1,zs1)*test_zoverlap(zl2,zs2)){
    cov_NG_gl_gl_real_fft_binned(theta, Ntheta, dtheta, zl1,zs1,zl2,zs2, cov_ng);
    run_ng = 1;
  }

  print_citations(F1);

  FILE *F2;
  char filename2[300];
  int if_write = 0, i1;
  i1 = Ntheta*(2*tomo.shear_Npowerspectra+n1);
  sprintf(filename2,"%sorder_%s_i_%d-%d",PATH,covparams.filename,i1,i1+Ntheta-1);
  if (fopen(filename2, "r") == NULL){
    if_write=1;
    F2 =fopen(filename2,"w");
    fprintf(F2, "# i, bin-averaged theta_i, 2pt func, {s: source, l:lens}tomo bin index 1, {s, l}tomo bin index 2\n");
  }

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      c_ng = 0.;
      if(run_ng){c_ng = cov_ng[nl1][nl2];}
      c_g = cov_g[nl1][nl2];
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(2*tomo.shear_Npowerspectra+n2)+nl2,t[nl1],t[nl2],zl1,zs1,zl2,zs2,c_g,c_ng);
      //printf("%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(2*tomo.shear_Npowerspectra+n2)+nl2,theta[nl1],theta[nl2],zl1,zs1,zl2,zs2,c_g,c_ng);
    }
    if(if_write) {
      fprintf(F2, "%d %e gammat l%d s%d\n", i1+nl1, t[nl1], zl1, zs1);
    }
  }
  fclose(F1);
  if(if_write) {fclose(F2);}
  free(cov_g);free(cov_ng);
}


void run_cov_ggl_shear_real_binned(char *OUTFILE, char *PATH, double *t, double *theta, double *dtheta, int Ntheta, int n1, int n2, int pm, int start)
{
  int zl,zs,z3,z4,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  zl = ZL(n1); zs = ZS(n1);
  printf("\nN_ggl = %d (%d, %d)\n", n1,zl,zs);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);

  double **cov_g, **cov_ng;
  int i;
  cov_g = malloc(Ntheta*sizeof(double *));
  cov_ng = malloc(Ntheta*sizeof(double *));
  for(i=0;i<Ntheta;i++) {
    cov_g[i] = malloc(Ntheta*sizeof(double));
    cov_ng[i] = malloc(Ntheta*sizeof(double));
  }


  cov_G_gl_shear_real_fft_binned(theta, Ntheta, dtheta, zl,zs,z3,z4,pm, cov_g);
  int run_ng = 0;
  if (test_zoverlap(zl,zs)*test_zoverlap(zl,z3)*test_zoverlap(zl,z4) && covparams.ng){
    cov_NG_gl_shear_real_fft_binned(theta, Ntheta, dtheta, zl,zs,z3,z4,pm, cov_ng);
    run_ng=1;
  }

  print_citations(F1);

  FILE *F2;
  char filename2[300];
  int if_write = 0, i1;
  i1 = Ntheta*(2*tomo.shear_Npowerspectra+n1);
  sprintf(filename2,"%sorder_%s_i_%d-%d",PATH,covparams.filename,i1,i1+Ntheta-1);
  if (fopen(filename2, "r") == NULL){
    if_write=1;
    F2 =fopen(filename2,"w");
    fprintf(F2, "# i, bin-averaged theta_i, 2pt func, {s: source, l:lens}tomo bin index 1, {s, l}tomo bin index 2\n");
  }

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (run_ng){c_ng = cov_ng[nl1][nl2];}
      c_g = cov_g[nl1][nl2];
      if(pm==1) fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(n2)+nl2,t[nl1],t[nl2],zl,zs,z3,z4,c_g,c_ng);
      if(pm==0) fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+n1)+nl1,Ntheta*(tomo.shear_Npowerspectra+n2)+nl2,t[nl1],t[nl2],zl,zs,z3,z4,c_g,c_ng);
    }
    if(if_write) {
      fprintf(F2, "%d %e gammat l%d s%d\n", i1+nl1, t[nl1], zl, zs);
    }
  }
  fclose(F1);
  if(if_write) {fclose(F2);}
  free(cov_g);free(cov_ng);
}

void run_cov_clustering_shear_real_binned(char *OUTFILE, char *PATH, double *t, double *theta, double *dtheta,int Ntheta, int n1, int n2, int pm, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl = %d \n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);

  double **cov_g, **cov_ng;
  int i;
  cov_g = malloc(Ntheta*sizeof(double *));
  cov_ng = malloc(Ntheta*sizeof(double *));
  for(i=0;i<Ntheta;i++) {
    cov_g[i] = malloc(Ntheta*sizeof(double));
    cov_ng[i] = malloc(Ntheta*sizeof(double));
  }

  cov_G_cl_shear_real_fft_binned(theta, Ntheta, dtheta, z1,z2,z3,z4,pm, cov_g);
  int run_ng = 0;
  if (test_zoverlap(z1,z3)*test_zoverlap(z1,z4) && covparams.ng){
    cov_NG_cl_shear_real_fft_binned(theta, Ntheta, dtheta, z1,z2,z3,z4,pm, cov_ng);
    run_ng=1;
  }

  print_citations(F1);

  FILE *F2;
  char filename2[300];
  int if_write = 0, i1;
  i1 = Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1);
  sprintf(filename2,"%sorder_%s_i_%d-%d",PATH,covparams.filename,i1,i1+Ntheta-1);
  if (fopen(filename2, "r") == NULL){
    if_write=1;
    F2 =fopen(filename2,"w");
    fprintf(F2, "# i, bin-averaged theta_i, 2pt func, {s: source, l:lens}tomo bin index 1, {s, l}tomo bin index 2\n");
  }

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (run_ng){c_ng = cov_ng[nl1][nl2];}
      c_g = cov_g[nl1][nl2];
      if(pm==1)fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ntheta*(n2)+nl2,t[nl1],t[nl2],z1,z2,z3,z4,c_g,c_ng);
      if(pm==0)fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ntheta*(tomo.shear_Npowerspectra+n2)+nl2, t[nl1],t[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
    if(if_write) {
      fprintf(F2, "%d %e w l%d l%d\n", i1+nl1, t[nl1], z1, z2);
    }
  }
  fclose(F1);
  if(if_write) {fclose(F2);}
  free(cov_g);free(cov_ng);
}



void run_cov_clustering_ggl_real_binned(char *OUTFILE, char *PATH, double *t, double *theta, double *dtheta,int Ntheta, int n1, int n2, int start)
{
  int z1,z2,zl,zs,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl_1 = %d \n", n1);
  zl = ZL(n2); zs = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl,zs);

  double **cov_g, **cov_ng;
  int i;
  cov_g = malloc(Ntheta*sizeof(double *));
  cov_ng = malloc(Ntheta*sizeof(double *));
  for(i=0;i<Ntheta;i++) {
    cov_g[i] = malloc(Ntheta*sizeof(double));
    cov_ng[i] = malloc(Ntheta*sizeof(double));
  }

  cov_G_cl_gl_real_fft_binned(theta, Ntheta, dtheta, z1,z2,zl,zs, cov_g);
  int run_ng = 0;
  if (z1 == zl && covparams.ng && test_zoverlap(z1,zs)*test_zoverlap(zl,zs)){
    cov_NG_cl_gl_real_fft_binned(theta, Ntheta, dtheta, z1,z2,zl,zs, cov_ng);
    run_ng=1;
  }

  print_citations(F1);

  FILE *F2;
  char filename2[300];
  int if_write = 0, i1;
  i1 = Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1);
  sprintf(filename2,"%sorder_%s_i_%d-%d",PATH,covparams.filename,i1,i1+Ntheta-1);
  if (fopen(filename2, "r") == NULL){
    if_write=1;
    F2 =fopen(filename2,"w");
    fprintf(F2, "# i, bin-averaged theta_i, 2pt func, {s: source, l:lens}tomo bin index 1, {s, l}tomo bin index 2\n");
  }

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      c_ng = 0.;
      if (run_ng){c_ng = cov_ng[nl1][nl2];}
      c_g = cov_g[nl1][nl2];
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ntheta*(2*tomo.shear_Npowerspectra+n2)+nl2,t[nl1],t[nl2],z1,z2,zl,zs,c_g,c_ng);
    }
    if(if_write) {
      fprintf(F2, "%d %e w l%d l%d\n", i1+nl1, t[nl1], z1, z2);
    }
  }
  fclose(F1);
  if(if_write) {fclose(F2);}
  free(cov_g);free(cov_ng);
}

///////////
