/*
 * Using PK Stan library
 */

functions {
#include "../lib/utils.stan"
#include "../lib/model_lib.stan"

  // we fit a 1-cmt oral dosing situation
  matrix pk_system(vector lref, vector Dt, vector theta, real[] x_r, int[] x_i) {
    // as we are using an effect cmt, the in and out rate into the
    // effct cmt must be equal, hence k12=k20
    return(pk_1cmt_metabolite_depot(lref[1:3], Dt, theta[1], theta[2], theta[2], 0, 0));
    //return(pk_1cmt_eff_auc(lref[1:3], Dt, theta[1], theta[2], 0, 0));

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
    lstate_mdose = pk_1cmt_metabolite_depot(lref_mdose[1:3], Dt - tau * n, theta[1], theta[2], theta[2], tau, n+1);
//    lstate_mdose = pk_1cmt_eff_auc(lref_mdose[1:3], Dt - tau * n, theta[1], theta[2], tau, n+1);
    
    
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
  // reference time point (in months) 
  real<lower=0> tref_month;

  // log of pk parameters (1 is T_e half elimination rate, 2 is for k_e)
  vector[2] theta;

  real<lower=0> ref_dose;
  real<lower=0> ref_tau;
  // Number of provisional regimens
  int Nregimens;
  // Provisional regimens
  vector[Nregimens] regimens_lamt; // log(amt)
  vector[Nregimens] regimens_tau;  // frequencies
  int<lower=0> addls[Nregimens];
  // parameter for the prior distribution (logbeta)
  vector[2] params_prior;
}
transformed data {
  real x_r[0];
  int x_i[0];
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
  real time_unit;

  
  // we use month as time unit; time is given as h
  time_unit = 24.*7.*4.;
  obs_time_tref = rep_vector(time_unit * tref_month, 1);

//  {
  // Calculating ref_lscale
//  ref_lscale = to_vector(pk_model(dose_lamt_ref, dose_cmt_tref, dose_time_tref, dose_tau_ref, dose_addl_ref,
//                         init_lstate_tref, init_time_tref, obs_time_tref, 
//                         theta, 
//                         lscale,
//                         x_r, 
//                         x_i));
//  }
//  print("ref_lscale    = ", ref_lscale);
// Code from Sebastian
{
  // Calculating ref_lscale
    ref_lscale = to_vector(pk_model(dose_lamt_ref, dose_cmt_tref, dose_time_tref, dose_tau_ref, dose_addl_ref,
                                    init_lstate_tref, init_time_tref, obs_time_tref, 
                                    theta, 
                                    lscale,
                                    x_r, 
                                    x_i));

    //ref_lscale[2] = ref_lscale[3] - theta[2];
    //ref_lscale[1] = ref_lscale[3] - theta[2];
    //ref_lscale[3] = ref_lscale[3];
  }
  print("ref_lscale    = ", ref_lscale);


  // Calculate End-of-cycle-1 log(AUC) to calculate DLT rate of each dose
   for(i in 1:Nregimens) {  
     dose_lamt_tref[i] = rep_vector(regimens_lamt[i], 1);
     dose_tau_tref[i] = rep_vector(regimens_tau[i], 1);
   }
  print("dose_lamt_tref = ", dose_lamt_tref);
  print("dose_tau_tref = ", dose_tau_tref);

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
  print("laucs_tref = ", laucs_tref);

}

model {
}

generated quantities {
  real logbeta_prior;
  real lH_ref;
  real<lower=0, upper=1> P_ref;
  vector[Nregimens] P_dose;
  matrix[Nregimens,3] P_cat;
    
  logbeta_prior = normal_rng(params_prior[1], params_prior[2]);

  // note: 1-F(T>t_ref) = 1-exp(-Href) = Pref
  // => log(Href) = log(-log(1-Pref)) = cloglog(Pref)
  //lH_dose_ref = log(-log1m(P_dose_ref));
  lH_ref = logbeta_prior;
  P_ref = inv_cloglog(lH_ref);
  
  // The log of cumulative hazard of dose 30 at the end of first cycle
  for(i in 1:Nregimens) 
    P_dose[i] = inv_cloglog(logbeta_prior + laucs_tref[i]); 

  // calculating Posterior interval probabilities
  for (i in 1:Nregimens) {
    P_cat[i,1] = step(0.16 - P_dose[i]); 
    P_cat[i,2] = step(0.33 - P_dose[i]) - step(0.16 - P_dose[i]); 
    P_cat[i,3] = step(1 - P_dose[i]) - step(0.33 - P_dose[i]);
  }
}  
