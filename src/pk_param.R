# These parameters should be predefined
# T_e
# k_e
# ref_dose
# ref_tau
# tref_month (length of first cycle)
# tau_drug (drug frequency in hours, eg. daily = 24)

# Setting the reference regimen (dose and frequency)
ref_regimen <- list(ref_dose = ref_dose, ref_tau = ref_tau)
# Fixing the PK model parameters (T_e and k_e)
pk_theta <- log(c(T_e, k_e))
params   <- list(theta = pk_theta, ref_dose = ref_regimen$ref_dose, ref_tau = ref_regimen$ref_tau)
## calculations for the scalings
## NOTE: PK constants are set in units of [h] while TTE-PK model
## uses month as time-unit which is accounted for by appropiate scalings
time_unit  <- 24 * 7 * 4
tref_h     <- tref_month * time_unit
dose_addl  <- addl
# Creating refernce regimen to be used by pk_model function
dosing_ref <- data.frame(dose_lamt=log(ref_dose), dose_cmt=1, dose_time=0, dose_tau=ref_tau, 
                         dose_addl=dose_addl)

ref_lscale <- do.call(pk_model, c(dosing_ref, list(init_lstate = rep(-35, 3), init_time = 0, 
                                                   obs_time = tref_h, theta = params$theta, 
                                                   lscale=rep(0,3), 
                                                   x_r = 1.0, x_i = 1L )))

ref_lscale <- rep(ref_lscale[3], 3)
#ref_lscale[1] <-  ref_lscale[3] - params$theta[2]
#ref_lscale[2] <-  ref_lscale[3] - params$theta[2]
