/*
 * Using PK Stan library
 */

functions {

/*
 * Author: S. Weber
 *
 * Various utlities
 */

// helper function for rle_int
int rle_elem_count(int[] set) {
  int U;
    U = 1;
    for(i in 2:num_elements(set)) {
      if(set[i-1] != set[i])
	U = U + 1;
    }
    return(U);
}

// repeated length encoding, see rle in R
int[] rle_int(int[] set) {
  int res[rle_elem_count(set)];
  int c;
  res[1] = 1;
  c = 1;
  for(i in 2:num_elements(set)) {
    if(set[i-1] == set[i]) {
      res[c] = res[c] + 1;
    } else {
      c = c + 1;
      res[c] = 1;
    }
  }
  return(res);
}

// create an integer sequence
int[] seq_int(int start, int end) {
  int N = end - start + 1;
  int seq[N];
  for(i in 1:N) seq[i] = i + start - 1;
  return(seq);
}

// return an int vector with each element in the input repeated as
// many times as given in the each argument
int[] rep_each(int[] set, int[] each) {
  int N = sum(each);
  int replicated[N];
  int p = 1;

  for(i in 1:size(set)) {
    replicated[p:p+each[i]-1] = rep_array(set[i], each[i]);
    p = p + each[i];
  }

  return(replicated);
}

/* calculate the absolute value of a - b in log-space with log(a)
   and log(b) given. Does so by squaring and taking the root, i.e.
   
   la = log(a)
   lb = log(b)
   
   sqrt( (a - b)^2 ) = sqrt( a^2 - 2 * a * b + b^2 )
   
   <=> 0.5 * log_diff_exp(log_sum_exp(2*la, 2*lb), log(2) + la + lb)
*/
real log_diff_exp_abs(real la, real lb) {
  return(0.5 * log_diff_exp(log_sum_exp(2*la, 2*lb), log(2) + la + lb));
}


/* find_interval, see findInterval from R
 * i.e. find the ranks of x in sorted; sorted is assumed to be weakly-sorted.
 */
int[] find_interval_slow(vector x, vector sorted) {
  int res[num_elements(x)];
  // very brutal and ineffcient way of doing this, but well, it's
  // C++ speed...
  for(i in 1:num_elements(x)) {
    res[i] = rank(append_row(rep_vector(x[i], 1), sorted), 1);
  }
  return(res);
}

/* faster version which uses bisectioning search
 */
int find_interval_elem(real x, vector sorted, int start_ind) {
  int res;
  int N;
  int max_iter;
  real left;
  real right;
  int left_ind;
  int right_ind;
  int iter;
    
  N = num_elements(sorted);
  
  if(N == 0) return(0);
  
  left_ind  = start_ind;
  right_ind = N;
  
  max_iter = 100 * N;
  left  = sorted[left_ind ] - x;
  right = sorted[right_ind] - x;
  
  if(0 <= left)  return(left_ind-1);
  if(0 == right) return(N-1);
  if(0 >  right) return(N);
  
  iter = 1;
  while((right_ind - left_ind) > 1  && iter != max_iter) {
    int mid_ind;
    real mid;
    // is there a controlled way without being yelled at with a
    // warning?
    mid_ind = (left_ind + right_ind) / 2;
    mid = sorted[mid_ind] - x;
    if (mid == 0) return(mid_ind-1);
    if (left  * mid < 0) { right = mid; right_ind = mid_ind; }
    if (right * mid < 0) { left  = mid; left_ind  = mid_ind; }
    iter = iter + 1;
  }
  if(iter == max_iter)
    print("Maximum number of iterations reached.");
  return(left_ind);
}

int[] find_interval(vector x, vector sorted) {
  int res[num_elements(x)];
  for(i in 1:num_elements(x)) {
    res[i] = find_interval_elem(x[i], sorted, 1);
  }
  return(res);
}

// takes as input x an ascending sorted vector x which allows to
// move the left starting index to be moved
int[] find_interval_asc(vector x, vector sorted) {
  int res[num_elements(x)];
  int last;
  last = 1;
  for(i in 1:num_elements(x)) {
    res[i] = find_interval_elem(x[i], sorted, last);
    if(res[i] > 0) last = res[i];
  }
  return(res);
}

int[] find_interval_blocked(int[] vals_M, vector vals, int[] sorted_M, vector sorted) {
  int res[num_elements(vals)];
  int M;
  int v;
  int s;
  M = num_elements(vals_M);
  v = 1;
  s = 1;
  for(m in 1:M) {
    int temp[vals_M[m]];
    temp = find_interval(segment(vals, v, vals_M[m]), segment(sorted, s, sorted_M[m]));
    for(n in 1:vals_M[m])
      res[v + n - 1] = temp[n];
    v = v + vals_M[m];
    s = s + sorted_M[m];
  }
  return(res);
}

// count number times elem appears in test set
int count_elem(int[] test, int elem) {
  int count;
  count = 0;
  for(i in 1:num_elements(test))
    if(test[i] == elem)
      count = count + 1;
  return(count);
}

// count number times elems appears in test set
int[] count_elems(int[] test, int[] elems) {
  int counts[num_elements(elems)];
  for(i in 1:num_elements(elems))
    counts[i] = count_elem(test, elems[i]);
  return(counts);
}

// find elements in test which are equal to elem
int[] which_elem(int[] test, int elem) {
  int res[count_elem(test, elem)];
  int ci;
  ci = 1;
  for(i in 1:num_elements(test))
    if(test[i] == elem) {
      res[ci] = i;
      ci = ci + 1;
    }
  return(res);
}

// divide fac by div and return the rounded down integer
int floor_div_int(real fac, real div) {
  int count;
  if(fac < 0)
    reject("floor_div_int only works for positive values.");
  count = 1;
  while(count * div <= fac) { count = count + 1; }
  count = count - 1;
  return count;
}

int[] count_obs_event_free(int[] obs_timeRank, int ndose) {
  int dose_next_obs[ndose];
  int o;
  int O;
  dose_next_obs = rep_array(0, ndose);
  o = 0;
  O = size(obs_timeRank);
  while (o < O && obs_timeRank[o+1] == 0) { o = o + 1; }
  for (i in 1:ndose) {
    int count;
    count = 0;
    while(o < O && obs_timeRank[o+1] == i) {
      o = o + 1;
      count = count + 1;
    }
    dose_next_obs[i] = count;
  }
  return(dose_next_obs);
}

int[] count_obs_event_free_blocked(int[] M, int[] obs_timeRank, int[] ndose) {
  int dose_next_obs[sum(ndose)];
  int l;
  int ld;
  dose_next_obs = rep_array(0, sum(ndose));
  l = 1;
  ld = 1;
  for (i in 1:size(M)) {
    int u;
    int ud;
    u = l + M[i] - 1;
    ud = ld + ndose[i] - 1;
    dose_next_obs[ld:ud] = count_obs_event_free(obs_timeRank[l:u], ndose[i]);
    l = u + 1;
    ld = ud + 1;
  }
  return(dose_next_obs);
}

/*
 * subsets the input data structure to the indices given as second
 * argument
 */
int[] subset_int(int[] cand, int[] ind_set) {
  int out[size(ind_set)];
  for(i in 1:size(ind_set))
    out[i] = cand[ind_set[i]];
  return out;
}
  
vector subset_vec(vector cand, int[] ind_set) {
  vector[size(ind_set)] out;
  for(i in 1:size(ind_set))
    out[i] = cand[ind_set[i]];
  return out;
}
  
matrix subset_matrix(matrix cand, int[] ind_set) {
  matrix[size(ind_set),cols(cand)] out;
  for(i in 1:size(ind_set))
    out[i] = cand[ind_set[i]];
  return out;
}
  
// check that assumption of 1...J labeling of rows hold and warn if
// not
void check_ids(int[] id) {
  int cid = 0;
  int warned = 0;
  cid = 0;
  for(n in 1:num_elements(id)) {
    if(id[n] != cid) {
      if(id[n] != cid + 1) {
        if(!warned)
          print("WARNING: id vector not correctly sorted, i.e. not in range 1..J. Consider using the cid vector internally.");
        warned = 1;
      } else {
        cid = cid + 1;
      }
    }
  }
  if(max(id) != cid)
    print("WARNING: Last patient's id not equal to max(id).");
}

// check that addl dose coding is correct. That is, we can only have a
// single active addl dose such that addl dosing must be off whenever
// the next dose occurs.
void check_addl_dosing(vector dose_time, vector dose_tau, int[] dose_addl) {
  int D = num_elements(dose_time);

  for(d in 2:D) {
    if(dose_time[d] < (dose_time[d-1] + dose_tau[d-1] * dose_addl[d-1]))
      reject("Forbidden overlapping dosing records found.");
  }
}

void check_addl_dosing_blocked(int[] dose_M, vector dose_time, vector dose_tau, int[] dose_addl) {
  int M = num_elements(dose_M);
  int b = 1;
  for(m in 1:M) {
    check_addl_dosing(segment(dose_time, b, dose_M[m]),
                      segment(dose_tau, b, dose_M[m]),
                      segment(dose_addl, b, dose_M[m]));
    b = b + dose_M[m];
  }
}

  // turn a slicing variable for a ragged array
  // S = {5, 6, 11}
  // into
  // Si = {0, 5, 5 + 6, 5+6+11} + 1
  // such that we can index the ragged array A as
  // A[Si[u] : Si[u+1]-1]
  // for the uth unit
  int[] make_slice_index(int[] S) {
    int Si[size(S)+1];
    int cv = 1;
    Si[1] = cv;
    for(i in 1:size(S)) {
      cv = cv + S[i];
      Si[i+1] = cv;
    }
    return(Si);
  }

int[] count_dose_given(vector time, vector dose_time, vector dose_tau, int[] dose_addl) {
  int dose_count[num_elements(time)];
  int time_rank[num_elements(time)];
  int o;
  int O;
  o = 1;
  O = num_elements(time);
  time_rank = find_interval(time, dose_time);
  dose_count = rep_array(0, O);
  //print("time_rank = ", time_rank);
  // first skip all dosings before the first time
  while(o < O && time_rank[o] == 0) { o = o + 1; }
  //print("o = ", o);
  for(i in o:O) {
    int d;
    d = time_rank[i];
    if(dose_tau[d] > 0)
      dose_count[i] = min(floor_div_int(time[i] - dose_time[d], dose_tau[d]), dose_addl[d]);
  }
  return dose_count;
}

int[] count_dose_given_blocked(int[] M, vector time, int[] M_dose, vector dose_time, vector dose_tau, int[] dose_addl) {
  int dose_count[num_elements(time)];
  int B;
  int tl;
  int dl;
  B = num_elements(M);
  tl = 1;
  dl = 1;
  for(b in 1:B) {
    int tu;
    int du;
    tu = tl + M[b] - 1;
    du = dl + M_dose[b] - 1;
    dose_count[tl:tu] = count_dose_given(time[tl:tu], dose_time[dl:du], dose_tau[dl:du], dose_addl[dl:du]);
    tl = tu + 1;
    dl = du + 1;
  }
  return dose_count;
}

/*
 * Need to define an operator which develops a state forward in time
 * by an amount Dt. Input is the state at t, the log-coefficients of
 * the exponentials and the exponential exponents and Dt. Output is
 * the new state after Dt time units for one of the cmts (depends on
 * given coefs).
 *
 * lstate_ref is the initial
 *
 * coefsN is a 2-row matrix with the first row being the
 * log-coefficients and the second row the exponents.
 */
matrix evolve_lsystem(int S, vector Dt, matrix coefs, int[] coefs_map) {
  matrix[num_elements(Dt),S] lsystem;
  int T;
  T = num_elements(Dt);

  // initialize to zero
  lsystem = rep_matrix(-500, T, S);

  for(o in 1:cols(coefs)) {
    int s;
    vector[T] term;
    s    = coefs_map[o];
    term = coefs[1,o] + Dt * coefs[2,o];
    for(t in 1:T)
      lsystem[t,s] = log_sum_exp(lsystem[t,s], term[t]);
  }
  return(lsystem);
}

real lgeometric_series(real la, int n) {
  return( log1m_exp(la * n) - log1m_exp(la) );
}

/*
 * models 2 cmts 1 and 2 which each have an elimination rate k1 / k20
 * and we have a flow from 1 to 2, k12. This models a 1cmt oral dosing
 * (k1=k12) and/or a building block of a metabolite. Finally all mass
 * exiting cmt 2 is recorded in cmt 3.
 *
 * cmt1 -- k12 --> cmt2
 *  |               |
 *  k1              k20
 *  |               |
 * None            cmt3 / Depot
 */
matrix pk_1cmt_metabolite_depot(vector lref, vector Dt, real lk1, real lk12, real lk20, real tau, int n) {
  matrix[2,8] coefs1;
  int coefs1_map[8];
  int coefs1_zero[8];
  matrix[2,6] coefs2;
  int coefs2_map[6];
  int coefs2_zero[6];
  // positive terms
  matrix[num_elements(Dt),3] lsystem1;
  // negative terms (due to Bateman function)
  matrix[num_elements(Dt),2] lsystem2;
  real ldeltaSq;
  real nk1;
  real nk20;

  coefs1_zero = rep_array(0, 8);
  coefs2_zero = rep_array(0, 6);

  ldeltaSq = 2*log_diff_exp_abs(lk1, lk20);

  nk1  = -exp(lk1);
  nk20 = -exp(lk20);

  // setup coefficient matrix and coef index vectors
  coefs1_map[1] = 1;
  coefs1[1,1] = lref[1];
  coefs1[2,1] = nk1;

  coefs1_map[2] = 2;
  coefs1[1,2] = lref[1] + lk12 - ldeltaSq + lk1;
  coefs1[2,2] = nk20;
  coefs1_map[3] = 2;
  coefs1[1,3] = lref[1] + lk12 - ldeltaSq + lk20;
  coefs1[2,3] = nk1;

  coefs1_map[4] = 2;
  coefs1[1,4] = lref[2];
  coefs1[2,4] = nk20;

  // whatever is in the depot cmt doesn´t go away
  coefs1_map[5] = 3;
  coefs1[1,5] = log_sum_exp(lref[3], lref[2]);
  coefs1[2,5] = 0;
  coefs1_zero[5] = 1;

  coefs1_map[6] = 3;
  coefs1[1,6] = lref[1] + lk12 + lk20 - ldeltaSq + log_sum_exp(lk1 - lk20, lk20 - lk1);
  coefs1[2,6] = 0;
  coefs1_zero[6] = 1;

  coefs1_map[7] = 3;
  coefs1[1,7] = lref[1] + lk12 + lk20 - ldeltaSq;
  coefs1[2,7] = nk1;

  coefs1_map[8] = 3;
  coefs1[1,8] = lref[1] + lk12 + lk20 - ldeltaSq;
  coefs1[2,8] = nk20;

  // for the negative terms we only use a two cmts; hence 2 is
  // relabeled to 1, and 3 to 2
  coefs2_map[1] = 1;
  coefs2[1,1] = lref[1] + lk12 - ldeltaSq + lk1;
  coefs2[2,1] = nk1;
  coefs2_map[2] = 1;
  coefs2[1,2] = lref[1] + lk12 - ldeltaSq + lk20;
  coefs2[2,2] = nk20;

  coefs2_map[3] = 2;
  coefs2[1,3] = lref[2];
  coefs2[2,3] = nk20;

  coefs2_map[4] = 2;
  coefs2[1,4] = lref[1] + lk12 - ldeltaSq + lk20 + log(2);
  coefs2[2,4] = 0;
  coefs2_zero[4] = 1;

  coefs2_map[5] = 2;
  coefs2[1,5] = lref[1] + lk12 - ldeltaSq + lk20 + lk1 - lk20;
  coefs2[2,5] = nk20;

  coefs2_map[6] = 2;
  coefs2[1,6] = lref[1] + lk12 - ldeltaSq + lk20 + lk20 - lk1;
  coefs2[2,6] = nk1;

  // in case the initial state is dosed in a regular pattern, we can
  // take advantage of the geometric series here by modifing the
  // coefficients
  if(n>1) {
    real logn;
    logn = log(n);
    for(i in 1:8) {
      if(coefs1_zero[i]) {
        coefs1[1,i] = coefs1[1,i] + logn;
      } else {
        coefs1[1,i] = coefs1[1,i] + lgeometric_series(coefs1[2,i] * tau, n);
      }
    }
    for(i in 1:6) {
      if(coefs2_zero[i]) {
        coefs2[1,i] = coefs2[1,i] + logn;
      } else {
        coefs2[1,i] = coefs2[1,i] + lgeometric_series(coefs2[2,i] * tau, n);
      }
    }
  }

  //print("AFTER: coefs1 = ", coefs1);
  //print("AFTER: coefs2 = ", coefs2);

  lsystem1 = evolve_lsystem(3, Dt, coefs1, coefs1_map);
  lsystem2 = evolve_lsystem(2, Dt, coefs2, coefs2_map);

  //print("lsystem1 = ", lsystem1);
  //print("lsystem2 = ", lsystem2);

  // final system is the difference of the two solutions
  for(t in 1:num_elements(Dt)) {
    lsystem1[t,2] = log_diff_exp(lsystem1[t,2], lsystem2[t,1]);
    lsystem1[t,3] = log_diff_exp(lsystem1[t,3], lsystem2[t,2]);
  }

  return(lsystem1);
}

// same as above, but no depot cmt
matrix pk_1cmt_metabolite(vector lref, vector Dt, real lk1, real lk12, real lk20, real tau, int n) {
  matrix[2,4] coefs1;
  int coefs1_map[4];
  matrix[2,2] coefs2;
  int coefs2_map[2];
  // positive terms
  matrix[num_elements(Dt),2] lsystem1;
  // negative terms (due to Bateman function)
  matrix[num_elements(Dt),1] lsystem2;
  real ldeltaSq;
  real nk1;
  real nk20;

  ldeltaSq = 2*log_diff_exp_abs(lk1, lk20);

  nk1  = -exp(lk1);
  nk20 = -exp(lk20);

  // setup coefficient matrix and coef index vectors
  coefs1_map[1] = 1;
  coefs1[1,1] = lref[1];
  coefs1[2,1] = nk1;

  coefs1_map[2] = 2;
  coefs1[1,2] = lref[1] + lk12 - ldeltaSq + lk1;
  coefs1[2,2] = nk20;
  coefs1_map[3] = 2;
  coefs1[1,3] = lref[1] + lk12 - ldeltaSq + lk20;
  coefs1[2,3] = nk1;

  coefs1_map[4] = 2;
  coefs1[1,4] = lref[2];
  coefs1[2,4] = nk20;

  // for the negative terms we only use a two cmts; hence 2 is
  // relabeled to 1, and 3 to 2
  coefs2_map[1] = 1;
  coefs2[1,1] = lref[1] + lk12 - ldeltaSq + lk1;
  coefs2[2,1] = nk1;
  coefs2_map[2] = 1;
  coefs2[1,2] = lref[1] + lk12 - ldeltaSq + lk20;
  coefs2[2,2] = nk20;

  // in case the initial state is dosed in a regular pattern, we can
  // take advantage of the geometric series here by modifing the
  // coefficients
  if(n>1) {
    for(i in 1:4) {
      coefs1[1,i] = coefs1[1,i] + lgeometric_series(coefs1[2,i] * tau, n);
    }
    for(i in 1:2) {
      coefs2[1,i] = coefs2[1,i] + lgeometric_series(coefs2[2,i] * tau, n);
    }
  }

  //print("AFTER: coefs1 = ", coefs1);
  //print("AFTER: coefs2 = ", coefs2);

  lsystem1 = evolve_lsystem(2, Dt, coefs1, coefs1_map);
  lsystem2 = evolve_lsystem(1, Dt, coefs2, coefs2_map);

  //print("lsystem1 = ", lsystem1);
  //print("lsystem2 = ", lsystem2);

  // final system is the difference of the two solutions
  for(t in 1:num_elements(Dt)) {
    lsystem1[t,2] = log_diff_exp(lsystem1[t,2], lsystem2[t,1]);
  }

  return(lsystem1);
}

/*
 * models 1 cmts, an effect cmt and the AUC over the effect cmt. The
 * cmt3 is NOT the mass in there, but divided by ke such that it is
 * equal to the AUC in cmt E
 *
 * cmt1 -- ke --> cmtE
 *  |               |
 *  k1              ke
 *  |               |
 * None            AUC_cmtE
 */
// BELOW FUNCTION NOT CORRECT!!!
//matrix pk_1cmt_eff_auc(vector lref, vector Dt, real lk1, real lke, real tau, int n) {
//  int T = num_elements(Dt);
//  matrix[T,3] lsystem = pk_1cmt_metabolite_depot(lref, Dt, lk1, lke, lke, tau, n);
//  
//  for(i in 1:T) {
//    lsystem[i,3] = lsystem[i,3] - lke;
//  }
//  return(lsystem);
//}

// calculates PK metrics, for now only SS concentration (per unit of dose administered)
vector pk_1cmt_metabolite_metrics(real tau, real lk1, real lk12, real lk20) {
  vector[3] metrics;
  real k1;
  real k20;

  k1  = exp(lk1);
  k20 = exp(lk20);

  // SS in main cmt
  metrics[1] = - log1m_exp(-k1*tau) - k1*tau;

  // SS in metabolite cmt
  metrics[2] = lk12
    - log_diff_exp_abs(lk1, lk20)
    + log_diff_exp_abs(-k20*tau  - log1m_exp(-k20*tau), -k1 * tau - log1m_exp(-k1*tau));

  // Rac in main cmt
  metrics[3] = - log1m_exp(-k1*tau);

  return(metrics);
}

// Oral 2cmt analytical solutions: Needs alignment with other
// functions!

/** ANALYTICAL SOLUTION:
 * Calculates the 2-cmt model for one patient given as input nonmem
 * type data and as input parameters the logarith of the micro rate
 * and micro constants
 *
 * returns the log-concentration of the central compartement
 *
 * No ADDL dosing supported right now; needs update
 **/
matrix pk_oral_2cmt(vector state0, vector Dt,
                    real lka, real lalphaR, real lbetaR, real lA, real lB) {
  real lstateRefOral; // ref state for the 2-cmt with oral cmt (only the oral cmt)
  real lstateRef[2];  // ref state for the 2-cmt without oral cmt
  int N;
  real alphaR;
  real betaR;
  real ka;
  real lAt;
  real lBt;
  matrix[num_elements(Dt),3] lstate;
  real lk12;
  real lk21;
  real lD;
  real lad2;
  real lbd2;
  real ltemp;

  N = num_elements(Dt);

  ka = exp(lka);
  alphaR = exp(lalphaR);
  betaR  = exp(lbetaR);

  // Bateman coefficients
  lAt = lA + lka - log_diff_exp_abs(lka, lalphaR);
  lBt = lB + lka - log_diff_exp_abs(lka, lbetaR );

  // needed constant for the unobserved peripheral cmt C
  lD = log_diff_exp(lalphaR, lbetaR);   // discriminant which is always positive
  ltemp = log_sum_exp(lB + lD, lbetaR);
  lk12 = log_diff_exp(log_sum_exp(lalphaR, lbetaR), log_sum_exp(2*ltemp, lalphaR + lbetaR) - ltemp );
  lk21 = log_diff_exp(lalphaR, lA + lD);

  lad2 = 2 * log_diff_exp_abs(lalphaR, lka);
  lbd2 = 2 * log_diff_exp_abs(lbetaR , lka);

  // by convention time starts just at the first observation
  lstateRefOral = state0[1];
  lstateRef[1]  = state0[2];
  lstateRef[2]  = state0[3];
  for(i in 1:N) {
    lstate[i,1] = lstateRefOral - ka * Dt[i];
    // solution for the concentration which is in the central and
    // peripheral cmt
    lstate[i,2] = lstateRef[1] + log_sum_exp(lA - alphaR * Dt[i], lB - betaR * Dt[i]);
    lstate[i,3] = lstateRef[2] + log_sum_exp(lB - alphaR * Dt[i], lA - betaR * Dt[i]);

    // other changes in the state can only meaningful be calculated
    // if Dt[i] is large enough to allow for diffusion to occur
    if(Dt[i] > 0.) {
      lstate[i,2] = log_sum_exp(lstate[i,2], lstateRef[2] + lD - lk12 + lA + lB + log_diff_exp(- betaR * Dt[i], - alphaR * Dt[i]) );
      lstate[i,3] = log_sum_exp(lstate[i,3], lstateRef[1] + lk12 - lD + log_diff_exp(-betaR * Dt[i], -alphaR * Dt[i]));


      // add in the part which stems from oral cmt which results in
      // the superposition of Bateman functions in the main cmt
      lstate[i,2] = log_sum_exp(lstate[i,2], lstateRefOral + log_sum_exp(lAt + log_diff_exp_abs( -alphaR * Dt[i], -ka * Dt[i]),
                                                                         lBt + log_diff_exp_abs( -betaR  * Dt[i], -ka * Dt[i]))
                                );
      // last, we add into the peripheral cmt the effect of the oral
      // cmt dosing
      //lstate[i,3] = log_sum_exp(lstate[i,3], lstateRefOral + lk12 + lka - lD - ka * Dt[i] + log( D*A2i*B2i + A2i * exp(- (alphaR-ka) * Dt[i]) - B2i * exp(-(betaR-ka)*Dt[i])   ) );
      // the huge expression below is a sign-sorted version of (brute force)
      // k12 * ka /[ (alpha-ka) * (beta-ka) ] * [ exp(-ka * t) - (ka-beta)/(alpha-beta) * exp(-alpha*t) + (ka-alpha)/(alpha-beta) * exp(-beta*t) ]
      lstate[i,3] = log_sum_exp(lstate[i,3], lstateRefOral + lk12 + lka - lD - lad2 - lbd2 +
                                log_diff_exp(log_sum_exp(log_sum_exp(lD + log_sum_exp(lalphaR + lbetaR, 2*lka) - ka * Dt[i], lalphaR + lbd2 - alphaR * Dt[i]), lka    + lad2 - betaR * Dt[i] ),
                                             log_sum_exp(log_sum_exp(lD + lka + log_sum_exp(lalphaR,   lbetaR) - ka * Dt[i], lka     + lbd2 - alphaR * Dt[i]), lbetaR + lad2 - betaR * Dt[i] )
                                             )
                                );
    }
  }

  return lstate;
}

// forward declare pk system functions
matrix pk_system(vector lref, vector Dt, vector theta, real[] x_r, int[] x_i);


matrix pk_system_addl(vector lref, vector Dt, int cmt, real lamt, real tau, int n, vector theta, real[] x_r, int[] x_i);

// model evaluation function taking dosing (and addl dosing) into
// account for a single patient. The function calculates always with
// respect to a reference state. The initial reference state is the
// initial time and initial_lstate. Upon the first dosing event past
// the initial time, the reference time and state is set to the
// respective dosing event.
matrix pk_model_fast(vector dose_lamt, int[] dose_cmt, vector dose_time, vector dose_tau, int[] dose_addl, int[] dose_next_obs,
                     vector init_lstate, real init_time,
                     vector obs_time, int[] obs_timeRank, int[] obs_dose_given,
                     vector theta,
                     vector lscale,
                     real[] x_r, int[] x_i)
{
  int D = num_elements(dose_lamt);
  int O = num_elements(obs_time);
  int d = 0;
  int o = 1;
  int active_addl = 0; // 0 = FALSE
  int init_ref = 1;
  vector[num_elements(init_lstate)] lstate_ref = init_lstate;
  matrix[num_elements(obs_time),num_elements(init_lstate)] lstate;

  // skip all dosing records prior to init_time
  while(d < D && dose_time[d+1] < init_time) { d = d + 1; }
  // next, process all elements which are past active dosing
  while(o <= O) {
    // first update reference state to be just after the last
    // dosing
    while(d != obs_timeRank[o]) {
      int nd;
      // the dose of the reference state is not yet the very
      // last dose given before this observation, add it
      nd = d + 1;
      if(init_ref) {
        lstate_ref = to_vector(pk_system(lstate_ref, rep_vector(dose_time[nd] - init_time, 1), theta,
                                          x_r, x_i)[1]);
      } else if(active_addl) {
        // in case of an active addl record, we have to super-impose
        // the developed reference state with the emerged dosing which
        // is handled by the _addl function
        lstate_ref = to_vector(pk_system_addl(lstate_ref, rep_vector(dose_time[nd] - dose_time[d], 1), dose_cmt[d], dose_lamt[d], dose_tau[d], dose_addl[d], theta,
                                               x_r, x_i)[1]);
      } else {
        lstate_ref = to_vector(pk_system(lstate_ref, rep_vector(dose_time[nd] - dose_time[d], 1), theta,
                                          x_r, x_i)[1]);
      }
      // the new reference time point is the dosing event; time is now
      // measure as time-after-dose

      // add in the dosing, but only if we have a simple dosing
      // event, i.e. no additional dosings
      active_addl = dose_addl[nd] > 0;
      if(!active_addl) {
        lstate_ref[dose_cmt[nd]] = log_sum_exp(lstate_ref[dose_cmt[nd]], dose_lamt[nd]);
      }
      d = nd;
      // at this point the initial time is not any more the reference
      // state
      init_ref = 0;
    }
    // ok, evolve from reference (last dose or initial) to current
    // observation...
    if(init_ref) {
      lstate[o] = pk_system(lstate_ref, segment(obs_time, o, 1) - init_time, theta,
                             x_r, x_i)[1];
      o = o + 1;
    } else if(active_addl) {
      int ndose;
      int event_free;
      // ...in case of addl dosing, the effect of the multiple
      // dosing has not yet been added
      // number of dosings given from the active addl records
      ndose = obs_dose_given[o];
      // advance as far as we can by counting the number of
      // observations which have the same number of doses given
      event_free = 0;
      while((o + event_free) < O && obs_dose_given[o + event_free + 1] == ndose)
        event_free = event_free + 1;

      lstate[o:(o+event_free)] = pk_system_addl(lstate_ref, segment(obs_time, o, event_free + 1) - dose_time[d], dose_cmt[d], dose_lamt[d], dose_tau[d], ndose, theta,
                                                 x_r, x_i);
      o = o + event_free + 1;
    } else {
      // ... which is simple for non-addl dosing as dose is
      // already merged, evolve as long as no other dosing occurs
      int event_free;
      event_free = dose_next_obs[d];
      lstate[o:(o+event_free-1)] = pk_system(lstate_ref, segment(obs_time, o, event_free) - dose_time[d], theta,
                                              x_r, x_i);
      o = o + event_free;
    }
  }
  return(lstate - rep_matrix(to_row_vector(lscale), O));
}

matrix pk_model(vector dose_lamt, int[] dose_cmt, vector dose_time, vector dose_tau, int[] dose_addl,
                vector init_lstate, real init_time,
                vector obs_time,
                vector theta,
                vector lscale,
                real[] x_r, int[] x_i) {
  int obs_timeRank[num_elements(obs_time)] = find_interval(obs_time, dose_time);
  check_addl_dosing(dose_time, dose_tau, dose_addl);
  if (init_time > dose_time[1])
    reject("Initial time must be at or before first dose!");
  return(pk_model_fast(dose_lamt, dose_cmt, dose_time, dose_tau, dose_addl, count_obs_event_free(obs_timeRank, size(dose_cmt)),
                       init_lstate, init_time,
                       obs_time, obs_timeRank,
                       count_dose_given(obs_time, dose_time, dose_tau, dose_addl),
                       theta,
                       lscale,
                       x_r, x_i));
}

matrix evaluate_model_fast(int[] dose_M, vector dose_lamt, int[] dose_cmt, vector dose_time, vector dose_tau, int[] dose_addl, int[] dose_next_obs,
                           matrix init_lstate, vector init_time,
                           int[] obs_M, vector obs_time, int[] obs_timeRank, int[] obs_dose_given,
                           matrix theta,
                           matrix lscale,
                           real[] x_r, int[] x_i) {
  matrix[num_elements(obs_time), cols(init_lstate)] lstate;
  int J = num_elements(dose_M);
  int d = 1;
  int o = 1;
  int S = cols(lscale);

  for(j in 1:J) {
    matrix[obs_M[j],S] lstate_j;
    int d_m = dose_M[j];
    int o_m = obs_M[j];
    //print("Processing patient ", j);
    lstate_j = pk_model_fast(segment(dose_lamt, d, d_m), segment(dose_cmt, d, d_m), segment(dose_time, d, d_m), segment(dose_tau, d, d_m), segment(dose_addl, d, d_m), segment(dose_next_obs, d, d_m)
                              ,to_vector(init_lstate[j]), init_time[j]
                              ,segment(obs_time, o, o_m), segment(obs_timeRank, o, o_m), segment(obs_dose_given, o, o_m)
                              ,to_vector(theta[j])
                              ,to_vector(lscale[j])
                              ,x_r, x_i);

    for(i in 1:o_m)
      lstate[i + o - 1] = lstate_j[i];

    d = d + d_m;
    o = o + o_m;
  }
  return(lstate);
}

// returns only at the observed cmts
vector evaluate_model_fast_cmt(int[] dose_M, vector dose_lamt, int[] dose_cmt, vector dose_time, vector dose_tau, int[] dose_addl, int[] dose_next_obs,
                               matrix init_lstate, vector init_time,
                               int[] obs_M, vector obs_time, int[] obs_timeRank, int[] obs_dose_given, int[] obs_cmt,
                               matrix theta,
                               matrix lscale,
                               real[] x_r, int[] x_i) {
  vector[num_elements(obs_time)] lstate;
  int J = num_elements(dose_M);
  int d = 1;
  int o = 1;
  int S = cols(lscale);

  for(j in 1:J) {
    matrix[obs_M[j],S] lstate_j;
    int d_m = dose_M[j];
    int o_m = obs_M[j];
    //print("Processing patient ", j);
    lstate_j = pk_model_fast(segment(dose_lamt, d, d_m), segment(dose_cmt, d, d_m), segment(dose_time, d, d_m), segment(dose_tau, d, d_m), segment(dose_addl, d, d_m), segment(dose_next_obs, d, d_m)
                              ,to_vector(init_lstate[j]), init_time[j]
                              ,segment(obs_time, o, o_m), segment(obs_timeRank, o, o_m), segment(obs_dose_given, o, o_m)
                              ,to_vector(theta[j])
                              ,to_vector(lscale[j])
                              ,x_r, x_i);

    for(i in 1:o_m)
      lstate[i + o - 1] = lstate_j[i, obs_cmt[i + o - 1] ];

    d = d + d_m;
    o = o + o_m;
  }
  return(lstate);
}

matrix evaluate_model(int[] dose_M, vector dose_lamt, int[] dose_cmt, vector dose_time, vector dose_tau, int[] dose_addl,
                      matrix init_lstate, vector init_time,
                      int[] obs_M, vector obs_time,
                      matrix theta,
                      matrix lscale,
                      real[] x_r, int[] x_i) {
  int obs_timeRank[num_elements(obs_time)] = find_interval_blocked(obs_M, obs_time, dose_M, dose_time);
  check_addl_dosing_blocked(dose_M, dose_time, dose_tau, dose_addl);
  return(evaluate_model_fast(dose_M, dose_lamt, dose_cmt, dose_time, dose_tau, dose_addl, count_obs_event_free_blocked(obs_M, obs_timeRank, dose_M),
                             init_lstate, init_time,
                             obs_M, obs_time, obs_timeRank,
                             count_dose_given_blocked(obs_M, obs_time, dose_M, dose_time, dose_tau, dose_addl),
                             theta,
                             lscale,
                             x_r, x_i));
}

matrix evaluate_model_nm(int[] id, vector time, int[] cmt, int[] evid, vector amt, vector tau, int[] addl, int[] mdv,
                         matrix init_lstate, vector init_time,
                         matrix theta, matrix lscale,
                         real[] x_r, int[] x_i) {
  int dose_ind[count_elem(evid, 1)] = which_elem(evid, 1);

  return(evaluate_model(rle_int(id[dose_ind]), log(amt[dose_ind]), cmt[dose_ind], time[dose_ind], tau[dose_ind], addl[dose_ind],
                        init_lstate, init_time,
                        rle_int(id), time,
                        theta,
                        lscale,
                        x_r, x_i));
}

// returns only the states for which the respective row was specified
// for
vector evaluate_model_nm_cmt(int[] id, vector time, int[] cmt, int[] evid, vector amt, vector tau, int[] addl, int[] mdv,
                             matrix init_lstate, vector init_time,
                             matrix theta, matrix lscale,
                             real[] x_r, int[] x_i) {
  int dose_ind[count_elem(evid, 1)] = which_elem(evid, 1);
  int N = num_elements(time);
  matrix[N,cols(init_lstate)] states;
  vector[N] obs_states;

  states = evaluate_model(rle_int(id[dose_ind]), log(amt[dose_ind]), cmt[dose_ind], time[dose_ind], tau[dose_ind], addl[dose_ind],
                          init_lstate, init_time,
                          rle_int(id), time,
                          theta,
                          lscale,
                          x_r, x_i);
  for(i in 1:N)
    obs_states[i] = states[i,cmt[i]];
  
  return(obs_states);
}

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
