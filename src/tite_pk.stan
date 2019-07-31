/*
  A time-to-event pharmacokinetic modeel for
  phase I dose-escalation trials
  
  TITE-PK One Parameter model
  Authors: Burak Kuersad Guenhan and Sebastian Weber
 
 */

functions {
#include "../lib/utils.stan"
#include "../lib/model_lib.stan"

  // we fit a 1-cmt oral dosing situation
  matrix pk_system(vector lref, vector Dt, vector theta, real[] x_r, int[] x_i) {
    // as we are using an effect cmt, the in and out rate into the
    // effct cmt must be equal, hence k12=k20
    //return(pk_1cmt_metabolite_depot(lref[1:3], Dt, theta[1], theta[2], theta[2], 0, 0));
    return(pk_1cmt_eff_auc(lref[1:3], Dt, theta[1], theta[2], 0, 0));

  }

  // note that n are the additional doses to be added such that in total
  // n+1 are added
  matrix pk_system_addl(vector lref, vector Dt, int cmt, real lamt, real tau, int n, vector theta, real[] x_r, int[] x_i) {
    matrix[num_elements(Dt), num_elements(lref)] lstate;
    matrix[num_elements(Dt), num_elements(lref)] lstate_mdose;
    vector[num_elements(lref)] lref_mdose;
    int S;

    // evolve reference state freely...
    lstate = pk_system(lref, Dt, theta, x_r, x_i);

    // ... and add the extra doses correctly time-shifted
    S = num_elements(lref);
    lref_mdose = rep_vector(-35, S);
    lref_mdose[cmt] = lamt;
    // debugging code is disabled
    //if(prod(Dt - tau * n) < 0) reject("All requested times must be past the last dosing addl event.");
    /*
    for(i in 1:num_elements(Dt))
      if((Dt[i] - tau * n) < -1E-3) {
        print("tau = ", tau);
        print("n = ", n);
        print("Dt[i] - tau * n = ", Dt[i] - tau * n);
        reject("All requested times must be past the last dosing addl event.");
      }
    */
//    lstate_mdose = pk_1cmt_metabolite_depot(lref_mdose[1:3], Dt - tau * n, theta[1], theta[2], theta[2], tau, n+1);
    lstate_mdose = pk_1cmt_eff_auc(lref_mdose[1:3], Dt - tau * n, theta[1], theta[2], tau, n+1);
    
    
    for(s in 1:S)
      for(t in 1:num_elements(Dt))
        lstate[t,s] = log_sum_exp(lstate_mdose[t,s], lstate[t,s]);
    return(lstate);
  }


  /*
   * used to calculate the standardized cmt values
   */
  vector calc_ref_lscale(real ref_dose, real ref_tau, vector ref_theta) {
    vector[3] ref_lscale;
    vector[3] metrics;

    metrics = pk_1cmt_metabolite_metrics(ref_tau, ref_theta[1], ref_theta[2], ref_theta[2]);
    ref_lscale[1] = log(ref_dose) + metrics[1];
    ref_lscale[2] = log(ref_dose) + metrics[2];
    // scale the depot cmt by it's input rate such that it is the AUC
    // of the effect cmt (in addition to the re-scaling of the effect conc)
    // use as time unit 28 days (one month) instead of h

    ref_lscale[3] = ref_lscale[2] + ref_theta[2] + log(24 * 7 * 4);

    return(ref_lscale);
  }

  void pretty_print(real[,] x) {
    int d[2];
    d = dims(x);
    for (m in 1:d[1]) {
      row_vector[d[2]] rv;
      for (n in 1:d[2])
	rv[n] = round(1000*x[m,n])/1000.;
      print("row ", m, " = ", rv);
    }
  }

}


data {
  #include "../lib/nm_data.stan"
  int<lower=0, upper=1> dv[N]; // observations, 1 = DLT, 0 = censoring for mdv == 0 records

  // length of cyce 1 (in months) 
  real<lower=0> tref_month;

  // log of PK parameters (First element is half elimination rate (Te), 2 is for k_e)
  vector[2] theta;
  
  // Specifying refernce regimen (reference dose and frequency (tau))
  real<lower=0> ref_dose;
  real<lower=0> ref_tau;
  // Provisional regimens (NONMEM-style arguments)
  int Nregimens;                   // number of provisional regimens    
  vector[Nregimens] regimens_lamt; // log(amt)
  vector[Nregimens] regimens_tau;  // frequencies
  int<lower=0> addls[Nregimens];
  // parameter for the prior distribution (logbeta)
  vector[2] params_prior;
}

transformed data {
  #include "../lib/nm_defs.stan"
  vector[J] K;
  int Nobs = count_elem(cmt, 10);
  int Ncens = count_elem(cmt, 11);
  // pk concentrations/AUC at the events
  vector[Nobs] lpk;
  vector[Nobs] lpk_auc;
  vector[Nobs] time_obs;
  vector[Nobs] ltime_obs;
  int cid_obs[Nobs];
  int ind_obs[Nobs];
  // auc's needed for censoring
  vector[Ncens] lauc;
  vector[Ncens] time_cens;
  vector[Ncens] ltime_cens;
  int ind_cens[Ncens];
  // cumulative number of events time-point per subject
  vector[Nobs] Nev;
  vector[J] Ntot;
  real time_unit;
  real ltime_unit;
  real ltref;
  // ref scalings which are calculated using the reference
  // regimen. The scaling is chosen such that at the reference
  // time-point the AUC is equal to 1.
  // Calculating ref_lscale (Rescaling factor)
  vector[1] dose_lamt_ref = rep_vector(log(ref_dose), 1);
  vector[1] dose_tau_ref = rep_vector(ref_tau, 1);
  int dose_addl_ref[1] = rep_array(730, 1);
  vector[3] lscale = rep_vector(0, 3);
  vector[3] ref_lscale;

  // For calculating End-of-cycle-1 log(AUC) of each dose
  vector[1] dose_lamt_tref[Nregimens];
  int dose_cmt_tref[1] = rep_array(1, 1);
  vector[1] dose_time_tref = rep_vector(0, 1);
  vector[1] dose_tau_tref[Nregimens];
  vector[3] init_lstate_tref = rep_vector(-35, 3);
  real init_time_tref = 0;
  vector[1] obs_time_tref;
  matrix[1,3] lauc_tref[Nregimens];
  vector[Nregimens] laucs_tref;
  
    // we use month as time unit; time is given as h
  time_unit = 24.*7.*4.;
  obs_time_tref = rep_vector(time_unit * tref_month, 1);
  ltime_unit = log(time_unit);


{
  // Calculating ref_lscale
    ref_lscale = to_vector(pk_model(dose_lamt_ref, dose_cmt_tref, dose_time_tref, dose_tau_ref, dose_addl_ref,
                                    init_lstate_tref, init_time_tref, obs_time_tref, 
                                    theta, 
                                    lscale,
                                    x_r, 
                                    x_i));

  }
  //print("ref_lscale    = ", ref_lscale);

  #include "../lib/nm_checks.stan"

  if(Ncens != J)
    reject("Censoring times must be given for all J patients!");

  K = rep_vector(0, J);

  ind_obs  = which_elem(cmt, 10);
  ind_cens = which_elem(cmt, 11);

  // subject id for each observation; taken from the continuously
  // labelled ID set
  cid_obs = cid[ind_obs];
  
  // added for the case when no event observed
  if(Nobs > 0)
  {
    int cur_j;
    cur_j = cid_obs[1];
    Nev = rep_vector(1.0, Nobs);
    Ntot = rep_vector(0.0, J);
    Ntot[cur_j] = 1;
    for(o in 2:Nobs) {
      if(cid_obs[o] == cur_j) {
	Nev[o] = Nev[o-1] + 1;
	Ntot[cur_j] = Nev[o];
      } else {
	cur_j = cid_obs[o];
	Ntot[cur_j] = 1;
      }
    }
  }
  //print("Nev    = ", Nev);
  //print("id_obs = ", id[ind_obs]);
  //print("cid_obs = ", cid_obs);
  //print("Ntot   = ", Ntot);
  print("Total events = ", sum(Ntot));
  //print("Nobs    = ", Nobs);
  //print("Patient number =", J);



  // censoring times (log), converted to month
  time_cens  = time[ind_cens] / time_unit;
  ltime_cens = log(time[ind_cens]) - ltime_unit;
  time_obs = time[ind_obs] / time_unit;
  ltime_obs = log(time[ind_obs]) - ltime_unit;

  ltref = log(tref_month);

  //print("theta[1:2] = ", theta[1:2]);

  // calculate PK metrics (the time-varying exposure metric covariate)
  {
    matrix[O,3] lpk_all;
    matrix[J,3] Init_lstate = rep_matrix(-25, J, 3);
    vector[J]   init_time   = rep_vector(0, J);
    int obs_ind_ev[Nobs]    = which_elem(cmt[obs_ind], 10);
    int obs_ind_cens[Ncens] = which_elem(cmt[obs_ind], 11);

    lpk_all = evaluate_model_fast(dose_M, dose_lamt, dose_cmt, dose_time, dose_tau, dose_addl, dose_next_obs,
                                   Init_lstate, init_time,
                                   obs_M, obs_time, obs_time_rank, obs_dose_given,
                                   rep_matrix(to_row_vector(theta), J),
                                   rep_matrix(to_row_vector(ref_lscale), J),
                                   x_r,
                                   x_i);

    // split concentration records from those only used to get the
    // AUC_inf and assign appropiate columns to output arrays
    lpk = col(subset_matrix(lpk_all, obs_ind_ev), 2);
    lpk_auc = col(subset_matrix(lpk_all, obs_ind_ev), 3);
    lauc  = col(subset_matrix(lpk_all, obs_ind_cens), 3);
  }
  
  // Calculate End-of-cycle-1 log(AUC) to calculate DLT rate of each dose
   for(i in 1:Nregimens) {  
     dose_lamt_tref[i] = rep_vector(regimens_lamt[i], 1);
     dose_tau_tref[i] = rep_vector(regimens_tau[i], 1);
   }
  //print("dose_lamt_tref = ", dose_lamt_tref);
  //print("dose_tau_tref = ", dose_tau_tref);

  {
   for(i in 1:Nregimens) {
      int dose_addl_tref[1] = rep_array(addls[i], 1);

      lauc_tref[i] = pk_model(dose_lamt_tref[i], dose_cmt_tref, dose_time_tref, dose_tau_tref[i], dose_addl_tref,
                         init_lstate_tref, init_time_tref, obs_time_tref, 
                         theta, 
                         ref_lscale,
                         x_r, x_i);
    }
  }
  
  // Log AUC values for calculating dose specific DLT rates
  for(i in 1:Nregimens)
    laucs_tref[i] = lauc_tref[i, 1, 3];
  //print("laucs_tref = ", laucs_tref);


}
parameters {
  real logbeta;
}
transformed parameters {
  real<lower=0, upper=1> P_ref;      // End-of-cycle-1 DLT rate of reference regimen
  real lH_ref;                       // End-of-cycle-1 Log cumulative hazard of reference regimen
  vector[Nobs] lh;
  vector[J] lHc;

  // note: 1-F(T>t_ref) = 1-exp(-Href) = Pref
  // => log(Href) = log(-log(1-Pref)) = cloglog(Pref)
  lH_ref = logbeta;
  P_ref = inv_cloglog(lH_ref);

  // calculate log hazard per observation
  for(o in 1:Nobs) lh[o] = logbeta + lpk[o];

  // calculate log(H(t)) at censoring times per subject
  for(j in 1:J) lHc[j] = logbeta + lauc[j];

}

model {

  // Putting prior on logbeta
  logbeta  ~ normal(params_prior[1], params_prior[2]);

  // add in terms of the log-lik
  // log-hazard per event
  target += lh;
  // per subject censoring, i.e. F(t_ev > t_c); F=exp(-H(t))
  target += -exp(lHc);
}

generated quantities {
  vector[Nregimens] P_dose;
  matrix[Nregimens,3] P_cat;
 
  // The log of cumulative hazard of dose 30 at the end of first cycle
  for(i in 1:Nregimens) 
    P_dose[i] = inv_cloglog(logbeta + laucs_tref[i]); 

  // calculating Posterior interval probabilities
  for (i in 1:Nregimens) {
    P_cat[i,1] = step(0.16 - P_dose[i]); 
    P_cat[i,2] = step(0.33 - P_dose[i]) - step(0.16 - P_dose[i]); 
    P_cat[i,3] = step(1 - P_dose[i]) - step(0.33 - P_dose[i]);
  }

}
