/******************************************************************************
CosmoLike Configuration Space Covariances for Projected Galaxy 2-Point Statistics
https://github.com/CosmoLike/CosmoCov
by CosmoLike developers
******************************************************************************/

// routines internal to this file
void set_cosmological_parameters();
void set_cov_parameters();
void set_survey_parameters();

int count_rows(char* filename,const char delimiter);
void init_source_sample(char *multihisto_file, int Ntomo);
void init_ggl_tomo();
void init_lens_sample(char *multihisto_file, int Ntomo);
void set_angular_binning(double *thetamin, double *dtheta);

void set_angular_binning(double *thetamin, double *dtheta){
  double *thetamax;
  thetamax=create_double_vector(0,like.Ntheta-1);
  if(covparams.lin_bins){ /* linear binning */
    double dt = (like.vtmax-like.vtmin)/like.Ntheta;
    for(int i=0; i<like.Ntheta; i++){
      thetamin[i]=like.vtmin+i*dt;
      thetamax[i]=like.vtmin+(i+1.0)*dt;
      dtheta[i]=thetamax[i]-thetamin[i];
    }
  }
  else{  /* logarithmic binning */
    double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
    for(int i=0; i<like.Ntheta; i++){
      thetamin[i]=exp(log(like.vtmin)+(i+0.0)*logdt);
      thetamax[i]=exp(log(like.vtmin)+(i+1.0)*logdt);
      dtheta[i]=thetamax[i]-thetamin[i];
    }
  }
  thetamin[like.Ntheta] = thetamax[like.Ntheta-1];
  free_double_vector(thetamax,0,like.Ntheta-1);
  like.theta_min = thetamin;
}

void init_source_sample(char *multihisto_file, int Ntomo)
{
  sprintf(redshift.shear_REDSHIFT_FILE,"%s",multihisto_file);
  redshift.shear_photoz=4;
  tomo.shear_Nbin = Ntomo;
  tomo.shear_Npowerspectra = tomo.shear_Nbin*(tomo.shear_Nbin+1)/2;
  printf("Source redshifts: multi-histo file %s, containing %d tomography bins\n",multihisto_file,tomo.shear_Nbin);
  for (int i=0;i<tomo.shear_Nbin; i++)
  {
    printf("bin %d: <z_s>=%f\n",i,zmean_source(i));
    //tomo.n_source[i]= n_source[i];
    nuisance.bias_zphot_shear[i]=0.0;
  }
  if(like.IA==0){printf("No IA \n");}
  else if(like.IA==1){
    printf("like.IA: %d, using parameters A_ia, eta_ia\n",like.IA);
    printf("A_ia = %lf, eta_ia = %lf \n", nuisance.A_ia, nuisance.eta_ia);
  }
  printf("init_source_sample complete\n");
}

void init_ggl_tomo(){
  if (tomo.clustering_Nbin ==0){
    printf("WARNING! init.c: init_ggl_tomo called while tomo.clustering_Nbin =0\n");
  }
  if (tomo.shear_Nbin ==0){
    printf("WARNING! init.c: init_ggl_tomo called while tomo.shear_Nbin =0\n");
  }
  int n = 0;
  for (int i = 0; i < tomo.clustering_Nbin; i++){
    for(int j = 0; j<tomo.shear_Nbin;j++){
      n += test_zoverlap(i,j);
      //printf("GGL combinations zl=%d zs=%d accept=%d; <z_l> = %.3f, <z_s> = %.3f\n",i,j,test_zoverlap(i,j), zmean(i),zmean_source(j));
    }
  }
  tomo.ggl_Npowerspectra = n;
  printf("%d GGL Powerspectra\n",tomo.ggl_Npowerspectra);
}

void init_lens_sample(char *multihisto_file, int Ntomo)
{
  sprintf(redshift.clustering_REDSHIFT_FILE,"%s",multihisto_file);
  redshift.clustering_photoz=4;
  tomo.clustering_Nbin = Ntomo;
  if(covparams.full_tomo){
    tomo.clustering_Npowerspectra = tomo.clustering_Nbin*(tomo.clustering_Nbin+1)/2;
  } else{
    tomo.clustering_Npowerspectra = tomo.clustering_Nbin;
  }
  printf("Lens redshifts: multi-histo file %s, containing %d tomography bins\n",multihisto_file,tomo.clustering_Nbin);
  printf("%d Clustering Powerspectra\n", tomo.clustering_Npowerspectra);
  pf_photoz(0.1,0);
  // for (int i=0;i<tomo.clustering_Nbin; i++)
  // {
  //   gbias.b1_function = & b1_per_bin;
  //   gbias.b[i] = b1[i];
  //   gbias.b_mag[i] = b_mag[i];
  //   nuisance.bias_zphot_clustering[i]=0.0;
  // }
  init_ggl_tomo();
  printf("init_lens_sample complete\n");
}



int count_rows(char* filename,const char delimiter){
  FILE *file = fopen (filename, "r" );
  char line[1000];
  if (file != NULL) {
    fgets(line,sizeof line,file);
    fclose(file);
  }
  else{
    printf("count_rows: file %s not found.\nEXIT\n",filename);
    exit(1);
  }
  int count = 1;
  char *p;

  p = line;
  while (*p != '\0')
  {
    if (*p == delimiter){
        while (*p == delimiter){p++;}
        count++;
    }
      p++;
    }
   return count;
}

void set_cov_parameters(char *covparamfile, int output)
{

  char line[256];
  int iline=0;

  FILE* input = fopen(covparamfile, "r");
  while(fgets(line, 256, input) != NULL)
  {
    char name[128],val[128];

    iline++;
    if(line[0] == '#') continue;

    sscanf(line, "%128s : %128s", name, val);
    if(strcmp(name, "tmin")==0)
    {
      sscanf(val, "%lf", &covparams.tmin);
      covparams.tmin*=constants.arcmin;
      if(output==1)
      {
        printf("tmin %f \n",covparams.tmin);
      }
      continue;
    }
    else if(strcmp(name, "tmax")==0)
    {
      sscanf(val, "%lf", &covparams.tmax);
      covparams.tmax*=constants.arcmin;
      if(output==1)
      {
        printf("tmax %f \n",covparams.tmax);
      }
      continue;
    }
    else if(strcmp(name, "ntheta")==0)
    {
      sscanf(val, "%d", &covparams.ntheta);
      if(output==1)
      {
        printf("ntheta %d \n",covparams.ntheta);
      }
      continue;
    }
    else if(strcmp(name, "linear_binning")==0)
    {
      sscanf(val, "%d", &covparams.lin_bins);
      if(output==1)
      {
        printf("linear binning %d \n",covparams.lin_bins);
      }
      continue;
    }
    else if(strcmp(name, "full_tomo")==0)
    {
      sscanf(val, "%d", &covparams.full_tomo);
      if(output==1)
      {
        printf("full tomo %d \n",covparams.full_tomo);
      }
      continue;
    }
    if(strcmp(name, "lmin")==0)
    {
      sscanf(val, "%lf", &covparams.lmin);
      if(output==1)
      {
        printf("lmin %f \n",covparams.lmin);
      }
      continue;
    }
    else if(strcmp(name, "lmax")==0)
    {
      sscanf(val, "%lf", &covparams.lmax);
      if(output==1)
      {
        printf("lmax %f \n",covparams.lmax);
      }
      continue;
    }
    else if(strcmp(name, "ncl")==0)
    {
      sscanf(val, "%d", &covparams.ncl);
      if(output==1)
      {
        printf("Ncl %d \n",covparams.ncl);
      }
      continue;
    }
    else if(strcmp(name, "ng")==0)
    {
      sscanf(val, "%d", &covparams.ng);
      if(output==1)
      {
        printf("ng %d \n",covparams.ng);
      }
      continue;
    }
    else if(strcmp(name, "cng")==0)
    {
      sscanf(val, "%d", &covparams.cng);
      if(output==1)
      {
        printf("cng %d \n",covparams.cng);
      }
      continue;
    }
    else if(strcmp(name, "outdir")==0)
    {
      sprintf(covparams.outdir,"%s",val);
      if(output==1)
      {
        printf("outdir %s \n",covparams.outdir);
      }
      continue;
    }
    else if(strcmp(name, "C_FOOTPRINT_FILE")==0)
    {
      sprintf(covparams.C_FOOTPRINT_FILE,"%s",val);
      if(output==1)
      {
        printf("C_FOOTPRINT_FILE %s \n",covparams.C_FOOTPRINT_FILE);
      }
      continue;
    }
    else if(strcmp(name, "c_footprint_file")==0)
    {
      sprintf(covparams.C_FOOTPRINT_FILE,"%s",val);
      if(output==1)
      {
        printf("C_FOOTPRINT_FILE %s \n",covparams.C_FOOTPRINT_FILE);
      }
      continue;
    }
    else if(strcmp(name, "filename")==0)
    {
      sprintf(covparams.filename,"%s",val);
      if(output==1)
      {
        printf("filename %s \n",covparams.filename);
      }
      continue;
    }
    else if(strcmp(name, "ss")==0)
    {
      sprintf(covparams.ss,"%s",val);
      if(output==1)
      {
        printf("ss %s \n",covparams.ss);
      }
      continue;
    }
    else if(strcmp(name, "ls")==0)
    {
      sprintf(covparams.ls,"%s",val);
      if(output==1)
      {
        printf("ls %s \n",covparams.ls);
      }
      continue;
    }
    else if(strcmp(name, "ll")==0)
    {
      sprintf(covparams.ll,"%s",val);
      if(output==1)
      {
        printf("ll %s \n",covparams.ll);
      }
      continue;
    }
  }
}

//various cosmological models

void set_cosmological_parameters(char *cosmofile, int output)
{
  char line[256];
  int iline=0;

  FILE* input = fopen(cosmofile, "r");
  while(fgets(line, 256, input) != NULL)
  {
    char name[128],val[128];

    iline++;
    if(line[0] == '#') continue;

    sscanf(line, "%128s : %128s", name, val);
    if(strcmp(name, "Omega_m")==0)
    {
      sscanf(val, "%lf", &cosmology.Omega_m);
      if(output==1)
      {
        printf("Omega_m %f \n",cosmology.Omega_m);
      }
      continue;
    }
    else if(strcmp(name, "Omega_v")==0)
    {
      sscanf(val, "%lf", &cosmology.Omega_v);
      if(output==1)
      {
        printf("Omega_v %f \n",cosmology.Omega_v);
      }
      continue;
    }
    else if(strcmp(name, "sigma_8")==0)
    {
      sscanf(val, "%lf", &cosmology.sigma_8);
      if(output==1)
      {
        printf("sigma_8 %f \n",cosmology.sigma_8);
      }
      continue;
    }
    else if(strcmp(name, "n_spec")==0)
    {
      sscanf(val, "%lf", &cosmology.n_spec);
      if(output==1)
      {
        printf("n_spec %f \n",cosmology.n_spec);
      }
      continue;
    }
    else if(strcmp(name, "w0")==0)
    {
      sscanf(val, "%lf", &cosmology.w0);
      if(output==1)
      {
        printf("w0 %f \n",cosmology.w0);
      }
      continue;
    }
    else if(strcmp(name, "wa")==0)
    {
      sscanf(val, "%lf", &cosmology.wa);
      if(output==1)
      {
        printf("wa %f \n",cosmology.wa);
      }
      continue;
    }
    else if(strcmp(name, "omb")==0)
    {
      sscanf(val, "%lf", &cosmology.omb);
      if(output==1)
      {
        printf("omb %f \n",cosmology.omb);
      }
      continue;
    }
    else if(strcmp(name, "h0")==0)
    {
      sscanf(val, "%lf", &cosmology.h0);
      if(output==1)
      {
        printf("h0 %f \n",cosmology.h0);
      }
      continue;
    }
    else if(strcmp(name, "coverH0")==0)
    {
      sscanf(val, "%lf", &cosmology.coverH0);
      if(output==1)
      {
        printf("coverH0 %f \n",cosmology.coverH0);
      }
      continue;
    }
    else if(strcmp(name, "rho_crit")==0)
    {
      sscanf(val, "%lf", &cosmology.rho_crit);
      if(output==1)
      {
        printf("rho_crit %f \n",cosmology.rho_crit);
      }
      continue;
    }
    else if(strcmp(name, "f_NL")==0)
    {
      sscanf(val, "%lf", &cosmology.f_NL);
      if(output==1)
      {
        printf("f_NL %f \n",cosmology.f_NL);
      }
      continue;
    }
    else if(strcmp(name, "pdelta_runmode")==0)
    {
      sprintf(pdeltaparams.runmode,"%s",val);
      if(output==1)
      {
        printf("runmode %s \n",pdeltaparams.runmode);
      }
      continue;
    }
    else if(strcmp(name, "A_s")==0)
    {
      sscanf(val, "%lf", &cosmology.A_s);
      if(output==1)
      {
        printf("f_NL %f \n",cosmology.A_s);
      }
      continue;
    }
  }
}


/////////// Survey parameters //////////////

void set_survey_parameters(char *surveyfile, int output)
{

  char line[256];
  int iline=0,i,j;
  double dummy_var;

  for(int i=0; i<20; i++){
    for(int j=0; j<20; j++){
      tomo.n_lens_ij[i][j]=0.;
    }
  }
  
  FILE* input = fopen(surveyfile, "r");
  while(fgets(line, 256, input) != NULL)
  {
    char name[128],val[256];

    iline++;

    if(line[0] == '#') continue;

    sscanf(line, "%128s : %128s", name, val);
    if(strcmp(name, "area")==0)
    {
      sscanf(val, "%lf", &survey.area);
      if(output==1)
      {
        printf("area %f \n",survey.area);
      }
      continue;
    }
    else if(strcmp(name, "m_lim")==0)
    {
      sscanf(val, "%lf", &survey.m_lim);
      if(output==1)
      {
        printf("mlim %f \n",survey.m_lim);
      }
      continue;
    }
    else if(strcmp(name, "name")==0)
    {
      sprintf(survey.name,"%s",val);
      if(output==1)
      {
        printf("sheartomo %s \n",survey.name);
      }
      continue;
    }
    else if(strcmp(name, "shear_REDSHIFT_FILE")==0)
    {
      sprintf(redshift.shear_REDSHIFT_FILE,"%s",val);
      if(output==1)
      {
        printf("shear_REDSHIFT_FILE %s \n",redshift.shear_REDSHIFT_FILE);
      }
      continue;
    }
    else if(strcmp(name, "clustering_REDSHIFT_FILE")==0)
    {
      sprintf(redshift.clustering_REDSHIFT_FILE,"%s",val);
      if(output==1)
      {
        printf("clustering_REDSHIFT_FILE %s \n",redshift.clustering_REDSHIFT_FILE);
      }
      continue;
    }
    else if(strcmp(name, "sourcephotoz")==0)
    {
      sprintf(survey.sourcephotoz,"%s",val);
      if(output==1)
      {
        printf("sourcephotoz %s \n",survey.sourcephotoz);
      }
      continue;
    }
    else if(strcmp(name, "lensphotoz")==0)
    {
      sprintf(survey.lensphotoz,"%s",val);
      if(output==1)
      {
        printf("lensphotoz %s \n",survey.lensphotoz);
      }
      continue;
    }
    else if(strcmp(name, "galsample")==0)
    {
      sprintf(survey.galsample,"%s",val);
      if(output==1)
      {
        printf("galsample %s \n",survey.galsample);
      }
      continue;
    }
    else if(strcmp(name, "source_tomobins")==0)
    {
      sscanf(val, "%d", &tomo.shear_Nbin);
      if(output==1)
      {
        printf("number of source tomo bins %d \n",tomo.shear_Nbin);
      }
      continue;
    }
    else if(strcmp(name, "lens_tomobins")==0)
    {
      sscanf(val, "%d", &tomo.clustering_Nbin);
      tomo.clustering_Npowerspectra=tomo.clustering_Nbin;
      if(output==1)
      {
        printf("number of lens tomo bins %d \n",tomo.clustering_Nbin);
      }
      continue;
    }
    else if(strcmp(name, "sigma_e")==0)
    {
      sscanf(val, "%lf", &survey.sigma_e);
      if(output==1)
      {
        printf("sigmae %f \n",survey.sigma_e);
      }
      continue;
    }
    else if(strcmp(name, "source_n_gal")==0)
    {
      if(strcmp(survey.sourcephotoz, "multihisto")==0)
      {
        i=0;
        for (char *p = strtok(val,","); p != NULL; p = strtok(NULL, ","))
        {
          double var =0;
          if (i<tomo.shear_Nbin)
          {
            sscanf(p, "%lf", &var);
            tomo.n_source[i]=var;
            if(output==1)
            {
              printf("tomo.n_source[%d]=%f \n",i,tomo.n_source[i]);
            }
          }
/*          if (i>0)
          {
            sscanf(p, "%lf", &var);
            tomo.n_source[i-1]=var;
            printf("tomo.n_source[%d]=%f \n",i-1,tomo.n_source[i-1]);
          }*/
          survey.n_gal+=var;
          i++;
        }
      }
      else
      {
        sscanf(val, "%lf", &survey.n_gal);
      }
      if(output==1)
      {
        printf("ngal %f \n",survey.n_gal);
      }
      continue;
    }
    else if(strcmp(name, "lens_n_gal")==0)
    {
      if(strcmp(survey.lensphotoz, "multihisto")==0)
      {
        i=0;
        for (char *p = strtok(val,","); p != NULL; p = strtok(NULL, ","))
        {
          double var =0.;
          if (i<tomo.clustering_Nbin)
          {
            sscanf(p, "%lf", &var);
            tomo.n_lens_ij[i][i]=var;
            if(output==1)
            {
              printf("tomo.n_lens_ij[%d][%d]=%f \n",i,i,tomo.n_lens_ij[i][i]);
            }
          }
/*          if (i>0)
          {
            sscanf(p, "%lf", &var);
            tomo.n_lens[i-1]=var;
            printf("tomo.n_lens[%d]=%f \n",i-1,tomo.n_lens[i-1]);
          }*/
          survey.n_lens+=var;
          i++;
        }
      }
      else
      {
        sscanf(val, "%lf", &survey.n_lens);
      }
      if(output==1)
      {
        printf("nlens %f \n",survey.n_lens);
      }
      continue;
    }
    else if(strcmp(name, "lens_n_ij")==0)
    {
      if(strcmp(survey.lensphotoz, "multihisto")==0)
      {
        sscanf(val, "%d[^,],%d[^,],%lf", &i, &j, &dummy_var);
        tomo.n_lens_ij[i][j] = dummy_var;
        tomo.n_lens_ij[j][i] = dummy_var;
        if(output==1)
        {
          printf("tomo.n_lens_ij[%d][%d]=%f \n",i,j, tomo.n_lens_ij[i][j]);
        }
      }
      else
      {
        printf("Error: lens_n_ij requires multihisto!\n"); exit(1);
      }
      continue;
    }
    else if(strcmp(name, "lens_tomogbias")==0)
    {
      i=0;
      for (char *p = strtok(val,","); p != NULL; p = strtok(NULL, ","))
      {
        double var;
        sscanf(p, "%lf", &var);
        gbias.b[i]=var;
        if(output==1)
        {
          printf("gbias[%d]=%f \n",i,gbias.b[i]);
        }
        i++;
      }
    }
    else if(strcmp(name, "lens_tomo_bmag")==0)
    {
      i=0;
      for (char *p = strtok(val,","); p != NULL; p = strtok(NULL, ","))
      {
        double var;
        sscanf(p, "%lf", &var);
        gbias.b_mag[i]=var;
        if(output==1)
        {
          printf("b_mag[%d]=%f \n",i,gbias.b_mag[i]);
        }
        i++;
      }
    }
    else if(strcmp(name, "lens_zphot_sigma")==0)
    {
      i=0;
      for (char *p = strtok(val,","); p != NULL; p = strtok(NULL, ","))
      {
        double var;
        sscanf(p, "%lf", &var);
        nuisance.sigma_zphot_clustering[i]=var;
        if(output==1)
        {
          printf("lens sigma[%d]=%f \n",i,nuisance.sigma_zphot_clustering[i]);
        }
        i++;
      }
    }
    else if(strcmp(name, "source_zphot_sigma")==0)
    {

      i=0;
      for (char *p = strtok(val,","); p != NULL; p = strtok(NULL, ","))
      {
        double var;
        sscanf(p, "%lf", &var);
        nuisance.sigma_zphot_shear[i]=var;
        if(output==1)
        {
          printf("source sigma[%d]=%f \n",i,nuisance.sigma_zphot_shear[i]);
        }
        i++;
      }
    }
    else if(strcmp(name, "IA")==0)
    {
      sscanf(val, "%d", &like.IA);
      if(like.IA!=0 && like.IA!=1){
        printf("like.IA = %d not supported in des_mpp\nEXIT\n", like.IA);
        exit(1);
      }
      if(output==1)
      {
        printf("IA = %d \n",like.IA);
      }
      continue;
    }
    else if(strcmp(name, "A_ia")==0 && like.IA==1)
    {
      sscanf(val, "%lf", &nuisance.A_ia);
      if(output==1)
      {
        printf("A_ia = %lf \n",nuisance.A_ia);
      }
      continue;
    }
    else if(strcmp(name, "eta_ia")==0 && like.IA==1)
    {
      sscanf(val, "%lf", &nuisance.eta_ia);
      if(output==1)
      {
        printf("eta_ia = %lf \n",nuisance.eta_ia);
      }
      continue;
    }
  }
  survey.area_conversion_factor = 60.0*60.0*constants.arcmin*constants.arcmin;
  survey.n_gal_conversion_factor=1.0/constants.arcmin/constants.arcmin;
}
