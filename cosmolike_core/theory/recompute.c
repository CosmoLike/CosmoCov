/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by CosmoLike developers
******************************************************************************/

void update_cosmopara (cosmopara *C);
void update_gal (galpara *G);
void update_nuisance (nuisancepara *N);
//conditions for recomputing look-up tables
int recompute_expansion(cosmopara C);
int recompute_Delta(cosmopara C);
int recompute_cosmo3D(cosmopara C);
int recompute_cosmo3D_CLASS(cosmopara C);
int recompute_zphot_shear(nuisancepara N);
int recompute_zphot_clustering(nuisancepara N);
int recompute_shear(cosmopara C, nuisancepara N); //for shear 2-pt statics
int recompute_ii(cosmopara C, nuisancepara N); //for shear 2-pt statics
int recompute_ggl(cosmopara C, galpara G, nuisancepara N,int i);//for gg-lensing statistics
int recompute_clustering(cosmopara C, galpara G, nuisancepara N, int i, int j);//clustering

void update_cosmopara (cosmopara *C){
  C->Omega_m = cosmology.Omega_m;
  C->Omega_v = cosmology.Omega_v;
  C->Omega_nu = cosmology.Omega_nu;
  C->M_nu = cosmology.M_nu;
  C->sigma_8 = cosmology.sigma_8;
  C->A_s = cosmology.A_s;
  C->n_spec = cosmology.n_spec;
  C->alpha_s = cosmology.alpha_s;
  C->w0 = cosmology.w0;
  C->wa = cosmology.wa;
  C->omb = cosmology.omb;
  C->h0 = cosmology.h0;
  C->f_NL = cosmology.f_NL;
  C->MGSigma = cosmology.MGSigma;
  C->MGmu = cosmology.MGmu;
  C->M_nu = cosmology.M_nu;
  C->theta_s = cosmology.theta_s;
}

void update_galpara (galpara *G){
  int i,j;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    if (gbias.b[i]> 0.2 && gbias.b[i] < 20){
      G->b[i] = gbias.b[i];
      G->b2[i] = gbias.b2[i];
      G->bs2[i] = gbias.bs2[i];
    }
    else{ printf("lens bin %d: neither HOD nor linear bias set, exit\n",i); exit(EXIT_FAILURE);}
  }
}

void update_nuisance (nuisancepara *N){
  int i;
  N->A_ia = nuisance.A_ia;
  N->eta_ia = nuisance.eta_ia;
  for(i = 0; i < tomo.clustering_Nbin; i++){
    N-> fred[i] = nuisance.fred[i];
    N->sigma_zphot_clustering[i] = nuisance.sigma_zphot_clustering[i];
    N->bias_zphot_clustering[i] = nuisance.bias_zphot_clustering[i];
  }
  for(i = 0; i < tomo.shear_Nbin; i++){
    N->sigma_zphot_shear[i] = nuisance.sigma_zphot_shear[i];
    N->bias_zphot_shear[i] = nuisance.bias_zphot_shear[i];
  }
}
int recompute_expansion(cosmopara C){ //rules for recomputing growth factor & comoving distance
  if (C.Omega_m != cosmology.Omega_m || C.Omega_v != cosmology.Omega_v || C.w0 != cosmology.w0 || C.wa != cosmology.wa || C.MGmu != cosmology.MGmu || C.M_nu != cosmology.M_nu){return 1;}
  if (cosmology.theta_s > 0 && C.theta_s != cosmology.theta_s){return 1;}
  else{return 0;}
}

int recompute_Delta(cosmopara C){ //rules for recomputing early time power spectrum Delta_L
  if (C.Omega_m != cosmology.Omega_m || C.Omega_v != cosmology.Omega_v || C.Omega_nu != cosmology.Omega_nu || C.M_nu != cosmology.M_nu || C.h0 != cosmology.h0 || C.omb != cosmology.omb || C.n_spec != cosmology.n_spec|| C.alpha_s != cosmology.alpha_s){return 1;}
  if (cosmology.A_s){
    if(C.A_s != cosmology.A_s){return 1;}
  }
  else{
    if (C.sigma_8 != cosmology.sigma_8){return 1;}
  }
  return 0;
}

int recompute_cosmo3D(cosmopara C){
  if (C.Omega_m != cosmology.Omega_m || C.Omega_v != cosmology.Omega_v || C.Omega_nu != cosmology.Omega_nu || C.M_nu != cosmology.M_nu || C.h0 != cosmology.h0 || C.omb != cosmology.omb || C.n_spec != cosmology.n_spec|| C.alpha_s != cosmology.alpha_s ||  C.w0 != cosmology.w0 || C.wa != cosmology.wa || C.MGSigma != cosmology.MGSigma || C.MGmu != cosmology.MGmu || C.M_nu != cosmology.M_nu){return 1;}
  if (cosmology.A_s){
     if(C.A_s != cosmology.A_s){return 1;}
  }
  else{
     if (C.sigma_8 != cosmology.sigma_8){return 1;}
  }
  if (cosmology.theta_s > 0 && C.theta_s != cosmology.theta_s){return 1;}
  return 0;
}
int recompute_cosmo3D_CLASS(cosmopara C){
  if (C.Omega_m != cosmology.Omega_m || C.Omega_v != cosmology.Omega_v || C.Omega_nu != cosmology.Omega_nu || C.M_nu != cosmology.M_nu || C.h0 != cosmology.h0 || C.omb != cosmology.omb || C.n_spec != cosmology.n_spec|| C.alpha_s != cosmology.alpha_s ||  C.w0 != cosmology.w0 || C.wa != cosmology.wa || C.MGSigma != cosmology.MGSigma || C.MGmu != cosmology.MGmu || C.M_nu != cosmology.M_nu){return 1;}
  if (cosmology.A_s > 0){
     if(C.A_s != cosmology.A_s){return 1;}
  }
  else{
     if (C.sigma_8 != cosmology.sigma_8){return 1;}
  }
  if (cosmology.theta_s > 0 && C.theta_s != cosmology.theta_s){return 1;}
  return 0;
}

int recompute_zphot_shear(nuisancepara N){
  static int photoz = -1;
  if (photoz != redshift.shear_photoz){photoz = redshift.shear_photoz; return 1;}
  if (redshift.shear_photoz != 3 && redshift.shear_photoz != 4){return 0;}
  int i, res = 0;
  for(i = 0; i < tomo.shear_Nbin; i++){
    if (N.sigma_zphot_shear[i]!= nuisance.sigma_zphot_shear[i] || N.bias_zphot_shear[i]!= nuisance.bias_zphot_shear[i]){ res = 1;}
  }
  return res;
}
int recompute_zphot_clustering(nuisancepara N){
  static int photoz = -1;
  if (photoz != redshift.clustering_photoz){photoz = redshift.clustering_photoz; return 1;}
  if (redshift.clustering_photoz != 3 && redshift.clustering_photoz != 4){return 0;}
  int i, res = 0;
  for(i = 0; i < tomo.clustering_Nbin; i++){
    if (N.sigma_zphot_clustering[i]!= nuisance.sigma_zphot_clustering[i] || N.bias_zphot_clustering[i]!= nuisance.bias_zphot_clustering[i]){ res = 1;}
  }
  return res;
}

int recompute_IA(nuisancepara N){
  if (N.A_ia != nuisance.A_ia || N.eta_ia != nuisance.eta_ia) return 1;
  return 0;
}

int recompute_shear(cosmopara C, nuisancepara N){
  if (recompute_cosmo3D(C) || recompute_zphot_shear(N)||recompute_IA(N)){return 1;}
  else{return 0;}
}

int recompute_ii(cosmopara C, nuisancepara N){
  if (recompute_cosmo3D(C) || recompute_zphot_clustering(N)|| recompute_zphot_shear(N)|| N.A_ia != nuisance.A_ia || N.eta_ia != nuisance.eta_ia){return 1;}
  else{return 0;} 
}

int recompute_galaxies(galpara G, int i){
 int j;
  if (i == -1){return 0;}
  if(G.b[i] != gbias.b[i] || G.b2[i] != gbias.b2[i] || G.bs2[i] != gbias.bs2[i]){return 1;}
  return 0;
}

int recompute_ggl(cosmopara C, galpara G, nuisancepara N, int i){
  if (recompute_cosmo3D(C) || recompute_zphot_clustering(N) || recompute_zphot_shear(N) || recompute_galaxies(G,i) ||recompute_IA(N) ){return 1;}
  else{return 0;}
}

int recompute_clustering(cosmopara C, galpara G, nuisancepara N, int i, int j){
  if (recompute_cosmo3D(C) || recompute_zphot_clustering(N) || recompute_galaxies(G,i)|| recompute_galaxies(G,j)){return 1;}
  else{return 0;}
 
}
