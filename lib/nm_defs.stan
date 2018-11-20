// standard NM dosing declarations
int D = count_elem(evid, 1);
int O = count_elem(mdv,  0);
int J = rle_elem_count(id);
int cid[N] = rep_each(seq_int(1, J), rle_int(id));

int dose_ind[D] = which_elem(evid, 1);
int dose_M[J] = rle_int(id[dose_ind]);
vector[D] dose_time = time[dose_ind];
vector[D] dose_tau = tau[dose_ind];
int dose_addl[D] = addl[dose_ind];
vector[D] dose_lamt = log(amt[dose_ind]);
int dose_cmt[D] = cmt[dose_ind];

// standard NM defs for observations
int obs_ind[O]  = which_elem(mdv , 0);
int obs_M[J] = rle_int(id[obs_ind]);
vector[O] obs_time = time[obs_ind];
int obs_cmt[O] = cmt[obs_ind];

// define an index vector pointing to the first entry of each subject,
// used for baseline covariates
int first_ind[J] = make_slice_index(rle_int(id))[1:J];

int obs_time_rank[O] = find_interval_blocked(obs_M, obs_time, dose_M, dose_time);
int obs_dose_given[O] = count_dose_given_blocked(obs_M, obs_time, dose_M, dose_time, dose_tau, dose_addl);

int dose_next_obs[D] = count_obs_event_free_blocked(obs_M, obs_time_rank, dose_M);
real x_r[0];
int x_i[0];
