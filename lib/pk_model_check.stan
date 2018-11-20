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

