/**
 * Computes the nonlinear correction on the linear power spectrum via
 * the method presented in Mead et al. 2009.01858
 *
 * @param ppr Input: pointer to precision structure
 * @param pba Input: pointer to background structure
 * @param ppt Input: pointer to perturbation structure
 * @param ppm Input: pointer to primordial structure
 * @param pnl Input: pointer to nonlinear structure
 * @param index_pk   Input: index of the pk type, either index_m or index_cb
 * @param index_tau  Input: index of tau, at which to compute the nl correction
 * @param tau        Input: tau, at which to compute the nl correction
 * @param pk_nl      Output:nonlinear power spectrum
 * @param lnpk_l     Input: logarithm of the linear power spectrum for both index_m and index_cb
 * @param ddlnpk_l   Input: spline of the logarithm of the linear power spectrum for both index_m and index_cb
 * @param nl_corr_not_computable_at_this_k Ouput: was the computation doable?
 * @param k_nl       Output: nonlinear scale for index_m and index_cb
 * @param pnw        Input/Output: pointer to nonlinear workspace
 * @return the error status
 */

int nonlinear_hmcode2020(
                     struct precision *ppr,
                     struct background *pba,
                     struct perturbs *ppt,
                     struct primordial *ppm,
                     struct nonlinear *pnl,
                     int index_pk,
                     int index_tau,
                     double tau,
                     double *pk_nl,
                     double **lnpk_l,
                     double **ddlnpk_l,
                     double *k_nl,
                     short * nl_corr_not_computable_at_this_k,
                     struct nonlinear_workspace * pnw
                     ) {
      /* integers */
      int index_mass, i, ng, nsig;
      int index_k, index_ncol;
      int last_index=0;
      int index_pk_cb;
      int counter, index_nl;

      int index_nu, index_cut;
      int index_y;
      int index_ddy;

      /* Background parameters */
      double Omega_m,fnu,Omega0_m Omega0_b,Omega0_c;
      double z_at_tau;
      double rho_crit_today_in_msun_mpc3;
      double growth;
      double intgrowth;
      double anorm;

      /* temporary numbers */
      double m, r, nu, sig, sigf;
      double diff, r1, r2;

      /* HMcode parameters */
      double mmin, mmax, nu_min;
      double dark_energy_correction;

      double sigma_disp, sigma_disp100, sigma8, sigma8cc;
      double delta_c, Delta_v;
      double fraction;

      double sigma_nl, nu_nl, r_nl;
      double sigma_prime;
      double dlnsigdlnR;
      double n_eff;
      double alpha;

      double z_form, g_form;

      double eta;
      double gst, window_nfw;
      double nu_cut;
      double fac, k_star, fdamp, nd, f, kd;
      double pk_lin, pk_2h, pk_1h, pk_wiggle;

      /* data fields */
      double * pvecback;
      double * conc;
      double * mass;
      double * sigma_r;
      double * sigmaf_r;
      double * r_virial;
      double * r_real;
      double * nu_arr;

      double * p1h_integrand;


      /** include precision parameters that control the number of entries in the growth and sigma tables */
      ng = ppr->n_hmcode_tables;
      nsig = ppr->n_hmcode_tables;

      /** Compute background quantitites today */

      Omega0_m = pba->Omega0_m;
      Omega0_b = pba->Omega0_b; 
      Omega0_c = pba->Omega0_cdm;
      fnu      = pba->Omega0_ncdm_tot/Omega0_m;

      /** If index_pk_cb, choose Omega0_cb as the matter density parameter.
       * If index_pk_m, choose Omega0_cbn as the matter density parameter. */
      if (index_pk==pnl->index_pk_cb){
        Omega0_m = Omega0_m - pba->Omega0_ncdm_tot;
      }

      anorm    = 1./(2*pow(_PI_,2));

      /** Call all the relevant background parameters at this tau */
      class_alloc(pvecback,pba->bg_size*sizeof(double),pnl->error_message);

      class_call(background_at_tau(pba,tau,pba->long_info,pba->inter_normal,&last_index,pvecback),
                 pba->error_message,
                 pnl->error_message);

      Omega_m = pvecback[pba->index_bg_Omega_m];//TBC (i.e. check if for P_cb here we should use Omega_cb) here the total time varying Omega_m is used for delta_c and for Delta_v according to the Mead fit of the Massara simulations.

      growth = pvecback[pba->index_bg_D];

      z_at_tau = 1./pvecback[pba->index_bg_a]-1.;
      class_call(array_interpolate_two_arrays_one_column(pnw->ztable,
                                                           pnw->intgrowtable,
                                                           1,
                                                           0,
                                                           z_at_tau,
                                                           &intgrowth,
                                                           pnl->error_message),
                   pnl->error_message, pnl->error_message);



      if (pnl->feedback==hmcode2020tagn){
        double theta_tagn = pnl->hmcode2020log10tagn-7.8;
        double B_0 = 3.44-0.496 * theta_tagn;
        double B_z = -0.0671-0.0371 * theta_tagn;
        pnl->c_min = B_0*pow(10, z_at_tau*B_z); 
        double fstar = (0.01*(2.01-0.030*theta_tagn))*pow(10, z_at_tau*(0.409+0.0224*theta_tagn));
        if(ftar>Omega0_b/Omega0_m) ftar = Omega0_b/Omega0_m;
        double Mb = pow(10, 13.87+1.81*theta_tagn+z_at_tau*(-0.108+0.195*theta_tagn))/pba->h; // unit Msun
        
      }
      else{
        pnl->c_min = 5.196;
      }




      /* The number below is the critical density today, rho_c = 3 * H0^2 / 8*pi*G, in units of M_sun over Mpc^3 */
      rho_crit_today_in_msun_mpc3 = 3.*pow(1.e5*pba->h, 2)/8./_PI_/_G_*_Mpc_over_m_/_M_SUN_;

      free(pvecback);

      /** Test whether pk_cb has to be taken into account (only if we have massive neutrinos)*/
      if (pba->has_ncdm==_TRUE_){
        index_pk_cb = pnl->index_pk_cb;
      }
      else {
        index_pk_cb = index_pk;
      }


      /** Get sigma(R=8 Mpc/h), sigma_disp(R=0), sigma_disp(R=100 Mpc/h) and write them into pnl structure */

      class_call(nonlinear_sigmas(pnl,
                                  8./pba->h,
                                  lnpk_l[index_pk],ddlnpk_l[index_pk],
                                  pnl->k_size_extra,
                                  ppr->sigma_k_per_decade,
                                  out_sigma,
                                  &sigma8),
                 pnl->error_message,
                 pnl->error_message);

        class_call(nonlinear_sigmas(pnl,
                                    8./pba->h,
                                    lnpk_l[index_pk_cb],ddlnpk_l[index_pk_cb],
                                    pnl->k_size_extra,
                                    ppr->sigma_k_per_decade,
                                    out_sigma,
                                    &sigma8cc),
                   pnl->error_message, pnl->error_message);


      class_call(nonlinear_sigmas(pnl,
                                  0.,
                                  lnpk_l[index_pk],ddlnpk_l[index_pk],
                                  pnl->k_size_extra,
                                  ppr->sigma_k_per_decade,
                                  out_sigma_disp,
                                  &sigma_disp),
                 pnl->error_message,
                 pnl->error_message);

      class_call(nonlinear_sigmas(pnl,
                                  100./pba->h,
                                  lnpk_l[index_pk],ddlnpk_l[index_pk],
                                  pnl->k_size_extra,
                                  ppr->sigma_k_per_decade,
                                  out_sigma_disp,
                                  &sigma_disp100),
                 pnl->error_message,
                 pnl->error_message);

      pnw->sigma_8[index_pk][index_tau] = sigma8;
      pnw->sigma_disp[index_pk][index_tau] = sigma_disp;
      pnw->sigma_disp_100[index_pk][index_tau] = sigma_disp100;

      /** Initialisation steps for the 1-Halo Power Integral */
      mmin=ppr->mmin_for_p1h_integral/pba->h; //Minimum mass for integration; (unit conversion from  m[Msun/h] to m[Msun]  )
      mmax=ppr->mmax_for_p1h_integral/pba->h; //Maximum mass for integration;

      class_alloc(mass,ppr->nsteps_for_p1h_integral*sizeof(double),pnl->error_message);
      class_alloc(r_real,ppr->nsteps_for_p1h_integral*sizeof(double),pnl->error_message);
      class_alloc(r_virial,ppr->nsteps_for_p1h_integral*sizeof(double),pnl->error_message);
      class_alloc(sigma_r,ppr->nsteps_for_p1h_integral*sizeof(double),pnl->error_message);
      class_alloc(sigmaf_r,ppr->nsteps_for_p1h_integral*sizeof(double),pnl->error_message);
      class_alloc(nu_arr,ppr->nsteps_for_p1h_integral*sizeof(double),pnl->error_message);

      // Linear theory density perturbation threshold for spherical collapse
      //delta_c = 1.59+0.0314*log(sigma8); //Mead et al. (2015; arXiv 1505.07833)
      //delta_c = delta_c*(1.+0.0123*log10(Omega_m)); //Nakamura & Suto (1997) fitting formula for LCDM models (as in Mead 2016)
      //delta_c = delta_c*(1.+0.262*fnu); //Mead et al. (2016; arXiv 1602.02154) neutrino addition
      delta_c = dc_Mead(pvecback[pba->index_bg_a], Omega_m, fnu, growth, intgrowth); // Mead et al. (2017; arXiv 1606.05345)

      // virialized overdensity
      //Delta_v=418.*pow(Omega_m, -0.352); //Mead et al. (2015; arXiv 1505.07833)
      //Delta_v=Delta_v*(1.+0.916*fnu); //Mead et al. (2016; arXiv 1602.02154) neutrino addition
      Delta_v= Dv_Mead(pvecback[pba->index_bg_a], Omega_m, fnu, growth, intgrowt); 

      // mass or radius fraction respectively
      fraction = pow(0.01, 1./3.);

      /* Fill the arrays needed for the P1H Integral: mass, r_real, r_virial, nu_arr, sigma_r, sigmaf_r
       * The P1H Integral is an integral over nu=delta_c/sigma(M), where M is connected to R via R=(3M)/(4*pi*rho_m).
       * The Integrand is M*Window^2{nu(M)*k, Rv(M), c(M)}*f(nu) with the window being the fouriertransformed
       * NFW profile, Rv = R/Delta_v^(1/3) and Sheth-Thormen halo mass function f.
       * The halo concentration-mass-relation c(M) will be found later.  */

      for (index_mass=0;index_mass<ppr->nsteps_for_p1h_integral;index_mass++){

        m = exp(log(mmin)+log(mmax/mmin)*(index_mass)/(ppr->nsteps_for_p1h_integral-1));
        r = pow((3.*m/(4.*_PI_*rho_crit_today_in_msun_mpc3*Omega0_m)), (1./3.));
        mass[index_mass] = m;
        r_real[index_mass] = r;
        r_virial[index_mass] = r_real[index_mass]/pow(Delta_v, 1./3.);

        class_call(array_interpolate_spline(pnw->rtab,
                                            nsig,
                                            pnw->stab,
                                            pnw->ddstab,
                                            1,
                                            r,
                                            &last_index,
                                            &sig,
                                            1,
                                            pnl->error_message),
                   pnl->error_message, pnl->error_message);

        class_call(array_interpolate_spline(pnw->rtab,
                                            nsig,
                                            pnw->stab,
                                            pnw->ddstab,
                                            1,
                                            r*fraction,
                                            &last_index,
                                            &sigf,
                                            1,
                                            pnl->error_message),
                   pnl->error_message, pnl->error_message);

        nu=delta_c/sig;
        sigma_r[index_mass] = sig;
        sigmaf_r[index_mass] = sigf;
        nu_arr[index_mass] = nu;
      }

      /** find nonlinear scales k_nl and r_nl and the effective spectral index n_eff */
      nu_nl = 1.;
      nu_min = nu_arr[0];

      /* stop calculating the nonlinear correction if the nonlinear scale is not reached in the table: */
      if (nu_min > nu_nl) {
        if (pnl->nonlinear_verbose>0) fprintf(stdout, " -> [WARNING:] the minimum mass in the mass-table is too large to find the nonlinear scale at this redshift.\n   Decrease mmin_for_p1h_integral\n");
        * nl_corr_not_computable_at_this_k = _TRUE_;
        free(mass);
        free(r_real);
        free(r_virial);
        free(sigma_r);
        free(sigmaf_r);
        free(nu_arr);
        return _SUCCESS_;
      }

      /* make a first guess for the nonlinear scale */
      class_call(array_interpolate_two_arrays_one_column(
                                                         nu_arr,
                                                         r_real,
                                                         1,
                                                         0,
                                                         ppr->nsteps_for_p1h_integral,
                                                         nu_nl,
                                                         &r_nl,
                                                         pnl->error_message),
                 pnl->error_message, pnl->error_message);

      class_call(array_search_bisect(ppr->nsteps_for_p1h_integral,nu_arr,nu_nl,&index_nl,pnl->error_message), pnl->error_message, pnl->error_message);

      r1 = r_real[index_nl-1];
      r2 = r_real[index_nl+2];

      /* // for debugging: (if it happens that r_nl is not between r1 and r2, which should never be the case)
         fprintf(stdout, "%e %e %e %e\n", r1, nu_arr[index_nl-1], r2, nu_arr[index_nl+2]);
      */

      /* do bisectional iteration between r1 and r2 to find the precise value of r_nl */
      counter = 0;
      do {
        r_nl = (r1+r2)/2.;
        counter ++;

        class_call(nonlinear_sigmas(pnl,
                                    r_nl,
                                    lnpk_l[index_pk_cb],ddlnpk_l[index_pk_cb],
                                    pnl->k_size_extra,
                                    ppr->sigma_k_per_decade,
                                    out_sigma,
                                    &sigma_nl),
                   pnl->error_message, pnl->error_message);

        diff = sigma_nl - delta_c;

        if (diff > ppr->hmcode_tol_sigma){
          r1=r_nl;
        }
        else if (diff < -ppr->hmcode_tol_sigma) {
          r2 = r_nl;
        }

        class_test(counter > _MAX_IT_,
                   pnl->error_message,
                   "could not converge within maximum allowed number of iterations");

      } while (fabs(diff) > ppr->hmcode_tol_sigma);

      if (pnl->nonlinear_verbose>5){
        fprintf(stdout, "number of iterations for r_nl at z = %e: %d\n", z_at_tau, counter);
      }
      *k_nl = 1./r_nl;

      if (*k_nl > pnl->k[pnl->k_size-1]) {
        * nl_corr_not_computable_at_this_k = _TRUE_;
        free(mass);
        free(r_real);
        free(r_virial);
        free(sigma_r);
        free(sigmaf_r);
        free(nu_arr);
        return _SUCCESS_;
      }
      else {
        * nl_corr_not_computable_at_this_k = _FALSE_;
      }

      /* call sigma_prime function at r_nl to find the effective spectral index n_eff */

      class_call(nonlinear_sigmas(pnl,
                                  r_nl,
                                  lnpk_l[index_pk_cb],ddlnpk_l[index_pk_cb],
                                  pnl->k_size_extra,
                                  ppr->sigma_k_per_decade,
                                  out_sigma_prime,
                                  &sigma_prime),
                 pnl->error_message,
                 pnl->error_message);

      dlnsigdlnR = r_nl*pow(sigma_nl, -2)*sigma_prime;
      n_eff = -3.- dlnsigdlnR;
      alpha = 1.875* pow(1.603, n_eff);
      pnw->sigma_prime[index_pk][index_tau] = sigma_prime;


      /** Compute the nonlinear correction */
      eta = 0.1281*pow(sigma8cc, -0.3644);; // halo bloating parameter table 2 of 2009.01858 (Mead)
      k_star = 0.05618*pow(sigma8cc, -1.013); //table 2 of 2009.01858 (Mead)
      nd = 2.853;
      kd = 0.05699*pow(sigma8cc, -1.089);
      f = 0.2696*pow(sigma8cc, 0.9403);



      /** Calculate halo concentration-mass relation conc(mass) (Bullock et al. 2001) */
      class_alloc(conc,ppr->nsteps_for_p1h_integral*sizeof(double),pnl->error_message);

      for (index_mass=0;index_mass<ppr->nsteps_for_p1h_integral;index_mass++){
        //find growth rate at formation
        g_form = delta_c*growth/sigmaf_r[index_mass];
        if (g_form > 1.) g_form = 1.;

        //
        class_call(array_interpolate_two_arrays_one_column(
                                                           pnw->growtable,
                                                           pnw->ztable,
                                                           1,
                                                           0,
                                                           ng,
                                                           g_form,
                                                           &z_form,
                                                           pnl->error_message),
                   pnl->error_message, pnl->error_message);
        if (z_form < z_at_tau){
          conc[index_mass] = pnl->c_min;
        } else {

             class_call(array_interpolate_two_arrays_one_column(
                                                           pnw->ztable,
                                                           pnw->dark_energy_correction_hmcode2020_table,
                                                           1,
                                                           0,
                                                           ng,
                                                           z_at_tau,
                                                           &dark_energy_correction,
                                                           pnl->error_message),






          conc[index_mass] = pnl->c_min*(1.+z_form)/(1.+z_at_tau)*dark_energy_correction;
        }
      }



    // growth factor

      /* the 1h integral contains the halo mass function proportional to exp(-nu^2).
       * To save time, the integration loop cuts, when nu exceeds a large value,
       * where the integrand is 0 anyhow. This cut index is found here. */
      nu_cut = 10.;
      if (nu_cut < nu_arr[ppr->nsteps_for_p1h_integral-1]){
        class_call(array_search_bisect(ppr->nsteps_for_p1h_integral,nu_arr,nu_cut,&index_cut,pnl->error_message), pnl->error_message, pnl->error_message);
      }
      else {
        index_cut = ppr->nsteps_for_p1h_integral;
      }

      i=0;
      index_nu=i;
      i++;
      index_y=i;
      i++;
      index_ddy=i;
      i++;
      index_ncol=i;

      for (index_k = 0; index_k < pnl->k_size; index_k++){

        class_alloc(p1h_integrand,index_cut*index_ncol*sizeof(double),pnl->error_message);

        pk_lin = exp(lnpk_l[index_pk][index_k])*pow(pnl->k[index_k],3)*anorm; //convert P_k to Delta_k^2

        class_call(array_interpolate_two_arrays_one_column(pnw->lnk_wiggle_tab,
                                                           pnw->Pk_wiggle_tab,
                                                           1,
                                                           0,
                                                           log(pnl->k[index_k]),
                                                           &pk_wiggle,
                                                           pnl->error_message),
                   pnl->error_message, pnl->error_message);




        pk_wiggle *= growth*growth; // growth factor
        pk_wiggle *= pow(pnl->k[index_k],3)*anorm; // P_k to Delta_k^2

        for (index_mass=0; index_mass<index_cut; index_mass++){ //Calculates the integrand for the ph1 integral at all nu values
          //get the nu^eta-value of the window
          class_call(nonlinear_hmcode_window_nfw(
                                                 pnl,
                                                 pow(nu_arr[index_mass], eta)*pnl->k[index_k],
                                                 r_virial[index_mass],
                                                 conc[index_mass],
                                                 &window_nfw),
                     pnl->error_message, pnl->error_message);

          if (pnl->feedback==hmcode2020tagn){
            fg = (Omega0_b/Omega0_m-fstar)/(1+pow(mass[index_mass]/Mb, -2)); //eq 24 in Mead 2009.01858
            window_nfw = window_nfw*(Omega0_c/Omega0_m+fg) + fstar;
          } 




          //get the value of the halo mass function
          class_call(nonlinear_hmcode_halomassfunction(
                                                       nu_arr[index_mass],
                                                       &gst),
                     pnl->error_message, pnl->error_message);

          p1h_integrand[index_mass*index_ncol+index_nu] = nu_arr[index_mass];

          p1h_integrand[index_mass*index_ncol+index_y] = mass[index_mass]*gst*pow(window_nfw, 2.);
          //if ((tau==pba->conformal_age) && (index_k == 0)) {
          //fprintf(stdout, "%d %e %e\n", index_cut, p1h_integrand[index_mass*index_ncol+index_nu], p1h_integrand[index_mass*index_ncol+index_y]);
          //}
        }
        class_call(array_spline(p1h_integrand,
                                index_ncol,
                                index_cut,
                                index_nu,
                                index_y,
                                index_ddy,
                                _SPLINE_EST_DERIV_,
                                pnl->error_message),
                   pnl->error_message,
                   pnl->error_message);

        class_call(array_integrate_all_trapzd_or_spline(
                                                        p1h_integrand,
                                                        index_ncol,
                                                        index_cut,
                                                        index_cut-1, //0 or n-1
                                                        index_nu,
                                                        index_y,
                                                        index_ddy,
                                                        &pk_1h,
                                                        pnl->error_message),
                   pnl->error_message,
                   pnl->error_message);


        fac = 1/(1+pow(pnl->k[index_k]/k_star,-4));

        pk_1h = pk_1h*anorm*pow(pnl->k[index_k],3)*(fac)/(rho_crit_today_in_msun_mpc3*Omega0_m);  // dimensionless power


        fdamp =  exp(-1*pnl->k[index_k]*pnl->k[index_k]*sigma_disp*sigma_disp)-1;
        if (fdamp<1.e-3) fdamp=1.e-3;
        if (fdamp>0.99)  fdamp=0.99;

        if (fdamp==0){
          pk_2h=pk_lin;
        }else{
          pk_2h= pk_lin + fdamp*pk_wiggle;   
        }
        if (pk_2h<0.) pk_2h=0.;

        pk_2h *= 1-f/(1+pow(pnl->k[index_k]/kd, -1*nd));


        pk_nl[index_k] = pow((pow(pk_1h, alpha) + pow(pk_2h, alpha)), (1./alpha))/pow(pnl->k[index_k],3)/anorm; //converted back to P_k

        free(p1h_integrand);
      }

      // print parameter values
      if ((pnl->nonlinear_verbose > 1 && tau==pba->conformal_age) || pnl->nonlinear_verbose > 3){
        fprintf(stdout, " -> Parameters at redshift z = %e:\n", z_at_tau);
        fprintf(stdout, "    fnu:		%e\n", fnu);
        fprintf(stdout, "    sigd [Mpc/h]:	%e\n", sigma_disp*pba->h);
        fprintf(stdout, "    sigd100 [Mpc/h]:    %e\n", sigma_disp100*pba->h);
        fprintf(stdout, "    sigma8:		%e\n", sigma8);
        fprintf(stdout, "    nu min:		%e\n", nu_arr[0]);
        fprintf(stdout, "    nu max:		%e\n", nu_arr[ppr->nsteps_for_p1h_integral-1]);
        fprintf(stdout, "    r_v min [Mpc/h]:    %e\n", r_virial[0]*pba->h);
        fprintf(stdout, "    r_v max [Mpc/h]:    %e\n", r_virial[ppr->nsteps_for_p1h_integral-1]*pba->h);
        fprintf(stdout, "    r_nl [Mpc/h]:	%e\n", r_nl*pba->h);
        fprintf(stdout, "    k_nl [h/Mpc]:	%e\n", *k_nl/pba->h);
        fprintf(stdout, "    sigma_nl:		%e\n", sigma_nl/delta_c);
        fprintf(stdout, "    neff:		%e\n", n_eff);
        fprintf(stdout, "    c min:		%e\n", conc[ppr->nsteps_for_p1h_integral-1]);
        fprintf(stdout, "    c max:		%e\n", conc[0]);
        fprintf(stdout, "    Dv:			%e\n", Delta_v);
        fprintf(stdout, "    dc:			%e\n", delta_c);
        fprintf(stdout, "    eta:		%e\n", eta);
        fprintf(stdout, "    k*:			%e\n", k_star/pba->h);
        fprintf(stdout, "    Abary:		%e\n", pnl->c_min);
        fprintf(stdout, "    fdamp:		%e\n", fdamp);
        fprintf(stdout, "    alpha:		%e\n", alpha);
        fprintf(stdout, "    ksize, kmin, kmax:   %d, %e, %e\n", pnl->k_size, pnl->k[0]/pba->h, pnl->k[pnl->k_size-1]/pba->h);
    }
      free(conc);
      free(mass);
      free(r_real);
      free(r_virial);
      free(sigma_r);
      free(sigmaf_r);
      free(nu_arr);

      return _SUCCESS_;
}

double _f_Mead(double x, double y, double p10, double p11, double p12, double p13) {
    return p10 + p11 * x + p12 * x * x + p13 * y;
}

double dc_Mead(double a, double Om_m, double f_nu, double g, double G) {
    double p10 = -0.0069, p11 = -0.0208, p12 = 0.0312, p13 = 0.0021;
    double p20 = 0.0001, p21 = -0.0647, p22 = -0.0417, p23 = 0.0646;
    double a1 = 1;
    double dc0 = 1;
    double dc_Mead = dc0 * (1. + log10(Om_m) * a1 * _f_Mead(g/a, G/a, p10, p11, p12, p13) + _f_Mead(g/a, G/a, p20, p21, p22, p23)) * (1. - 0.041 * f_nu);
    return dc_Mead;
}

double Dv_Mead(double a, double Om_m, double f_nu, double g, double G) {
    /*
    Delta_v fitting function from Mead (2017; 1606.05345)
    All input parameters should be evaluated as functions of a/z
    */
    // See Appendix A of Mead (2017) for naming convention
    double p30 = -0.79, p31 = -10.17, p32 = 2.51, p33 = 6.51;
    double p40 = -1.89, p41 = 0.38, p42 = 18.8, p43 = -15.87;
    double a3 = 1, a4 = 2;
    double Dv0 = 18.0 * pow(Om_m, 0.45);

    // Halo virial overdensity
    double Dv = 1.0;
    Dv = Dv + _f_Mead(g / a, G / a, p30, p31, p32, p33) * log10(Om_m) * a3;
    Dv = Dv + _f_Mead(g / a, G / a, p40, p41, p42, p43) * log10(Om_m) * a4;
    Dv = Dv * Dv0 * (1.0 + 0.763 * f_nu);
    return Dv;
}







