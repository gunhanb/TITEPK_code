// Initialize data structures from NM input data

// check that id vector is labelled from 1...J subsequently, but only
// warn the user about it as he then needs to use the cid vector
check_ids(id);

// these two checks are obsolete as the assignment in the nm_defs file
// will simply fail
if(J != rle_elem_count(id[dose_ind]))
  reject("Some patient(s) have no dosing event at all!");

if(J != rle_elem_count(id[obs_ind]))
  reject("Some patient(s) have no observation at all!");

// check that we have no overlapping doses due to addl coding
check_addl_dosing_blocked(dose_M, dose_time, dose_tau, dose_addl);
