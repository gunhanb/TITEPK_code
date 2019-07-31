#####################################################################
### 31.07.2019
### Authors: Burak Kürsad Günhan and Sebastian Weber
#####################################################################
### This is the accompanying code to paper 
### "A Bayesian time-to-event pharmacokinetic model for phase I
### dose-escalation trials with multiple schedules"
### See the paper for details
## Simulations for TITE-PK model
## Scenarios inspired by Wages et al 2018 
#####################################################################


##############################################################################
### LICENCE
### 
### THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
### IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
### FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
### AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
### LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
### OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
### SOFTWARE.
################################################################################

##### Libraries
library(assertthat)
library(rstan)
rstan_options(auto_write = TRUE)
options(rcpp.cache.dir="cache")

library(functional)
library(plyr)

# resproducibility
set.seed(356567)
#######################################################################
# Load oc functions from CECB 
#######################################################################
source("simulation/TitePK_Sim.R")
source("simulation/TitePKSimDecisionRule.R")
source("simulation/TitePKSimSummary.R")
#######################################################################
# Analysis
#######################################################################

# ------------------------------------------------
# Reference dose and schedule
# ------------------------------------------------
ref_dose = 24
ref_tau  = 24 * 4 #(hours) reference frequency of dosing
TRef    = 24 * 7 * 4 # reference time point
# ------------------------------------------------
# Priors for TITE PK
# ------------------------------------------------
cloglog <- binomial(link="cloglog")$linkfun
Prior_titepk = c(cloglog(0.30), 1.75)

# ------------------------------------------------
# Doses - 
# ------------------------------------------------

Doses.schedule1 = Doses.schedule2 = Doses.schedule3 = Doses.schedule4 = c(8, 16, 24)

# ------------------------------------------------
# Settings specific to TITE-PK model
# ------------------------------------------------
## compile Stan model functions via stanc_builder which avoids model
## name obfuscation and then admits caching in case nothing changes
stan_model_code_titepk <- stanc_builder("src/tite_pk.stan")
cat(stan_model_code_titepk$model_code, file="src/tite_pk_run.stan")
cat("\n", file="src/tite_pk_run.stan", append = TRUE)

cat("Exposing Stan functions...")
rstan::expose_stan_functions(stan_model_code_titepk)
cat("done.\n")
# Stan code for estimation
model_titepk <- stan_model("src/tite_pk_run.stan")

source("src/utils_tite_pk.R")
source("src/utils.R")
source("lib/tools.R")

T_e <- log(2) / 4
k_e <- 1.32

# Since cycle length is 4 weeks
tref_month = 1
# Assume the weekly dose frequency
tau_schedule1 = 24 * 8
addl = 91
# Calculate refernce scale
source("src/pk_param.R")
pk_calc_schedule1 <- Curry(pk_model, theta=params$theta, lscale=ref_lscale)
tau_drug <- tau_schedule1
# Calculating aucs
source("src/regimen_ref.R")
dosings <- regimen_ref_packed
dosings <- dosings[rep(seq_len(nrow(dosings)), each=length(Doses.schedule1)),]
dosings$amt  <- Doses.schedule1
dosings$lamt <- log(Doses.schedule1)
schedules1 = dosings
pksys <- ddply(dosings, .(amt), pk_subject, time = tref_h, pk_fun = pk_calc_schedule1)
# Only taking log(AUC_E(tref)) values for calculating dose specific DLT rates
laucs.schedule1 <- pksys[, c(4)]
#########################
### Calculating schedules2
#########################
tau_schedule2 = 24 * 4
addl = 182
# Calculate refernce scale
source("src/pk_param.R")
pk_calc_schedule2 <- Curry(pk_model, theta=params$theta, lscale=ref_lscale)
tau_drug <- tau_schedule2
# Calculating aucs
source("src/regimen_ref.R")
dosings <- regimen_ref_packed
dosings <- dosings[rep(seq_len(nrow(dosings)), each=length(Doses.schedule2)),]
dosings$amt  <- Doses.schedule2
dosings$lamt <- log(Doses.schedule2)
schedules2 = dosings
pksys <- ddply(dosings, .(amt), pk_subject, time = 24 * 7 * 4, pk_fun = pk_calc_schedule2)
laucs.schedule2 <- pksys[, c(4)]

### schedule 3
tau_schedule3 = 24 * 2
addl = 365
# Calculate refernce scale
source("src/pk_param.R")
pk_calc_schedule3 <- Curry(pk_model, theta=params$theta, lscale=ref_lscale)
tau_drug <- tau_schedule3
# Calculating aucs
source("src/regimen_ref.R")
dosings <- regimen_ref_packed
dosings <- dosings[rep(seq_len(nrow(dosings)), each=length(Doses.schedule3)),]
dosings$amt  <- Doses.schedule3
dosings$lamt <- log(Doses.schedule3)
schedules3 = dosings
pksys <- ddply(dosings, .(amt), pk_subject, time = 24 * 7 * 4, pk_fun = pk_calc_schedule3)
laucs.schedule3 <- pksys[, c(4)]


tau_schedule4 = 24 * 1
addl = 730
# Calculate refernce scale
source("src/pk_param.R")
pk_calc_schedule4 <- Curry(pk_model, theta=params$theta, lscale=ref_lscale)
tau_drug <- tau_schedule4
# Calculating aucs
source("src/regimen_ref.R")
dosings <- regimen_ref_packed
dosings <- dosings[rep(seq_len(nrow(dosings)), each=length(Doses.schedule4)),]
dosings$amt  <- Doses.schedule4
dosings$lamt <- log(Doses.schedule4)
schedules4 = dosings
pksys <- ddply(dosings, .(amt), pk_subject, time = 24 * 7 * 4, pk_fun = pk_calc_schedule4)
laucs.schedule4 <- pksys[, c(4)]

## Combining two schedules
schedules <- rbind(schedules1, schedules2, schedules3, schedules4)
laucs    <- c(laucs.schedule1, laucs.schedule2, laucs.schedule3, laucs.schedule4)
aucs = matrix(NA, ncol = 3 * 4, nrow = 3)
aucs[1,] = c(Doses.schedule1, Doses.schedule2, Doses.schedule3, Doses.schedule4)
aucs[2,] = c(rep(tau_schedule1, 3), rep(tau_schedule2, 3), rep(tau_schedule3, 3), rep(tau_schedule4, 3))
aucs[3,] = exp(laucs)
aucs = data.frame(aucs)
row.names(aucs) = c("Dose", "Freq", "AUC")
colnames(aucs) = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11", "L12")


# ------------------------------------------------
# Fixed - 
# ------------------------------------------------
Fixed = list(DoseRef = ref_dose,
             schedules = schedules,
             theta = log(c(T_e, k_e)),
             ref_tau = ref_tau,
             Prior   = Prior_titepk,
             AUC = aucs)

Pcutoffs <- c(0.20, 0.40)  # from Wages et al 2013
Crit = c(1, 1, 0.25)
# ---------------------------------------------------
# Considering Scenario 1 from the paper
# ---------------------------------------------------
True.Probs <- vector("list", 4)
True.Probs[[1]] <- c(0.05, 0.07, 0.11)
True.Probs[[2]] <- c(0.09, 0.12, 0.18)
True.Probs[[3]] <- c(0.16, 0.18, 0.23)
True.Probs[[4]] <- c(0.22, 0.26, 0.30)

True.Probs.all <- c(True.Probs[[1]], True.Probs[[2]], True.Probs[[3]], True.Probs[[4]])


True.logbetas <- vector("list", 4)

for(i in 1:length(Doses.schedule1)) {
  True.logbetas[[1]][i] <- cloglog(True.Probs[[1]][i]) - laucs.schedule1[i]
  True.logbetas[[2]][i] <- cloglog(True.Probs[[2]][i]) - laucs.schedule2[i]
  True.logbetas[[3]][i] <- cloglog(True.Probs[[3]][i]) - laucs.schedule3[i]
  True.logbetas[[4]][i] <- cloglog(True.Probs[[4]][i]) - laucs.schedule4[i]
}

True.logbetas <- round(c(True.logbetas[[1]], True.logbetas[[2]] , True.logbetas[[3]], True.logbetas[[4]]), 2)

## number of simulations
nTrials = 3


## Running the simulations
run.1 <- TitePKSim(Fixed = Fixed,
                   Decision = NULL,
                   nTrials = nTrials,
                   TrueProbs = True.Probs.all,
                   Truelogbeta = True.logbetas,
                   StartDose = 8,
                   StartFreq = tau_schedule1,
                   StartDoselevel = "L1",
                   CohortRnd = 1,
                   MaxN = 60,
                   mnMTD = 9,
                   mnTOT = 21,
                   mnTAR = 0.5,
                   Pcutoffs = Pcutoffs,
                   niter = 2000,
                   SeedParam = 3222,
                   ### Pay attention to EWOC!!!!!
                   Crit = Crit,
                   StanSeed = 111,
                   saveRes = "result/TITEPK_Scenario_1.Rdata",
                   Escalation = "maximal") 

###
###  The first realization 
###
print(run.1$out[[1]])

### Calculating summary statistics
out <- TitePKSimSummary(Simulation = run.1$out,
                        Fixed = Fixed, 
                        TrueProbs = True.Probs.all,
                        td.pos = 2,
                        Pcutoffs = Pcutoffs,
                        Outfile= "result/TITEPK_Scenario_1_summary.Rdata",
                        saveRes = "result/TITEPK_Scenario_1_summary.Rdata",
                        Ntrials = nTrials)

# ------------------------------------------------
# End of the program
# ------------------------------------------------
sessionInfo()



