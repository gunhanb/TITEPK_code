## Operating characteristics for TITE-PK

TitePKSim <- function(Fixed = NULL,
                      Decision = NULL,
                      Escalation = "maximal",
                      nTrials = NULL,
                      TrueProbs = NULL,
                      Truelogbeta = NULL,
                      StartDose = NULL,
                      StartFreq = NULL,
                      StartDoselevel = NULL,
                      CohortRnd = list(
                        c(3,     4,    5,    6),
                        c(0.8, 0.1, 0.05, 0.05)),
                      MaxN = NULL,
                      mnMTD = NULL, 
                      mnTOT = NULL, 
                      mnTAR = NUll,
                      Pcutoffs = NULL,
                      saveRes = NULL,
                      SeedParam = 1234,
                      Crit = NULL,
                      niter = niter,
                      StanSeed = 11){
  
  
  if(!is.null(Pcutoffs)){
    
    print("For Pcutoffs, the highest category will always be interpreted as overdose category.\n 
          The second highest category will always be interpreted as target category.")
    
    if(is.null(Crit)){
      stop("If Pcutoffs is not NULL, then Crit cannot be NULL.")
    }
    
    
    od.pos <- length(Pcutoffs) + 1
    td.pos <- od.pos - 1
    
  } else{
    
    od.pos <- 3
    td.pos <- 2
  }
  
  
  # ----------------------------------------------------
  #  Auxiliary steps
  #
  #  1) change all names of Fixed input to upper case
  #  2) initate object "simulation", a list object
  #     where results for each trial are stored
  #  3) set up seed for each possible cohort for the
  #     simulation (object all.seeds), length is 
  #     maximal number of trials * maximal number of
  #     cohorts * maximal number of attempts per cohort
  #  4) counter for each cohort, each attempt
  # ----------------------------------------------------
  
  names(Fixed) <- toupper(names(Fixed))
  Simulation <- list()
  simsummary <- list()
  NMATRX <- list()
  NNMATRX <- list()
  DLTMATRX <- list()
  MTDsim <- list()
  
  set.seed(SeedParam)
  #all.seeds <- sample(10*nTrials*ceiling(MaxN/min(CohortRnd[[1]])))
  all.seeds <- sample(10*nTrials*ceiling(MaxN/min(c(3, 4, 5))))
  
  counter <- 1
  
  #  sink(paste(Outfile, ".log", sep = ""))
  cat(paste("Simulation for ", nTrials, " trials started on ", Sys.time(), sep = ""))
  cat("\n\n")
  #  sink()
  
  # ----------------------------------------------------
  # Loop 1: over all trials (i = 1 to nTrials)
  # ----------------------------------------------------
  
  for(i in 1:nTrials){
    
    
    cat(paste("Trial: ", i, "\n", sep = ""))
    
    # ----------------------------------------------------
    #  Auxiliary steps
    #  1) initiate Ntot counting total number of patients in 
    #     the DDS for current trial
    #  2) indicator if MTD has been reached or not
    #  3) indicator if all doses are considered
    #     too toxic
    #  4) matrix for the i-th trial, containing output for
    #     each cohort of the simulated trial
    #  5) set current dose to start dose
    #  6) set cohort number to 1
    #  7) set crashed to FALSE
    # ----------------------------------------------------
    
    Ntot <- 0
    is.MTD <- FALSE
    too.toxic <- FALSE 
    #Simulation[[i]] <- matrix(nrow = ceiling(MaxN/min(CohortRnd[[1]])), ncol = 29)
    Simulation[[i]] <- matrix(nrow = 60, ncol = 29)
    
    
    colnames(Simulation[[i]]) = c("Trial", "Cohort", "Ntox", "Npat","DoseAdm", "FreqAdm", "Doselevel",
                                  "under", "target", "over", "true.prob", 
                                  "nextDoseAdm",  "nextFreqAdm", 
                                  "nextunder", "nexttarget", "nextover", "nextprob", 
                                  "MTD", "too.toxic", "Rhat.max", 
                                  "crash", "Rseed", "WBseed", 
                                  "T1", "T2", "T3", "T4", "T5", "T6") 
    
    Simulation[[i]] <- as.data.frame.matrix(Simulation[[i]])
    currentDoseAdm = StartDose
    currentDoselevel = StartDoselevel
    currentFreqAdm = StartFreq
    cohort <- 1
    crashed <- FALSE
    
    # ----------------------------------------------------
    # Loop 2: for i-th trial, loop over all cohorts
    # this loop stops if either of the following is reached
    # A) all dose combinations too toxic
    # B) MTD reached
    # C) Number of patients >= maximum number of patients
    #    (note: if MaxN / size of cohort not an integer,
    #    then one can end up with actual Ntot > MaxN)
    # D) The analysis did not crash (fully crash)
    # ----------------------------------------------------
    
    previous_sim_data <- list()
    while(!too.toxic & !is.MTD & (Ntot < MaxN) & !crashed){
      # ----------------------------------------------------
      #  Auxiliary steps
      #  1) store trial number
      #  2) store cohort number
      #  3) set seed
      #  4) store current seed 
      #  5) increase counter by 1
      #  6) sample cohort size or set to fixed value if fixed
      #  7) store cohort size
      #  8) update total N
      #  9) store current administered dose 
      # ---------------------------------------------------
      
      Simulation[[i]][cohort,"Trial"] <- i
      Simulation[[i]][cohort,"Cohort"] <- cohort
      set.seed(all.seeds[counter])
      Simulation[[i]][cohort,"Rseed"] <- all.seeds[counter]
      counter <- counter + 1
      
      #      if(length(CohortRnd[[1]]) == 1){
      #        CohortSize <- CohortRnd[[1]] 
      #      } else {
      #        CohortSize <- sample(CohortRnd[[1]],1,prob = CohortRnd[[2]])
      #      }
      CohortSize = CohortRnd
      
      Simulation[[i]][cohort,"Npat"] <- CohortSize
      Ntot <- Ntot + CohortSize     
      Simulation[[i]][cohort,"DoseAdm"] <- currentDoseAdm
      Simulation[[i]][cohort,"FreqAdm"] <- currentFreqAdm
      Simulation[[i]][cohort,"Doselevel"] <- currentDoselevel
      
      # ----------------------------------------------------
      #  Data generating step
      #  using TITE-PK model
      #  with true underlying probability
      # ----------------------------------------------------
      amount <- Simulation[[i]][cohort,"DoseAdm"]
      J      <- Simulation[[i]][cohort,"Npat"]
      tau    <- Simulation[[i]][cohort,"FreqAdm"]
      if(tau == 24 * 8) {addl = 91}
      if(tau == 24 * 4) {addl = 182}
      if(tau == 24 * 2) {addl = 365}
      if(tau == 24 * 1) {addl = 730}
      
      dummy <- rbind(data.frame(cmt = 1, time = 1, dv = 0, mdv = 1, evid = 1, 
                                amt = amount, lamt = log(amount), tau = tau, addl = addl), 
                     data.frame(cmt = 10, time = time_unit * tref_month, dv = 1, 
                                mdv = 0, evid = 0, amt = -15, lamt = 0, tau = 0, addl = 0))
      dummy <- do.call("rbind", replicate(J, dummy, simplify = FALSE))
      dummy$id <- rep(1:J, each = 2) + CohortSize * (cohort - 1)
      
      
      sim_data <- plyr::ddply(dummy, .(id), simulate_patient, pk_fun = pk_model, 
                              lbeta_coef = Truelogbeta[which(colnames(Fixed$AUC)==currentDoselevel)],  
                              leta = rep(0, J), cache = FALSE,
                              theta = params$theta, lscale = ref_lscale)
      
      
      
      Simulation[[i]][cohort,"Ntox"] <- nrow(subset(sim_data, dv == 1))
      # Maximum number of patients 6 
      # Saving the time-to-first-DLTs
      if(!is.na(subset(sim_data, dv == 0)$time[1])) {Simulation[[i]][cohort,"T1"] <- subset(sim_data, dv == 1)$time[1]}
      if(!is.na(subset(sim_data, dv == 0)$time[2])) {Simulation[[i]][cohort,"T2"] <- subset(sim_data, dv == 1)$time[2]}
      if(!is.na(subset(sim_data, dv == 0)$time[3])) {Simulation[[i]][cohort,"T3"] <- subset(sim_data, dv == 1)$time[3]}
      if(!is.na(subset(sim_data, dv == 0)$time[4])) {Simulation[[i]][cohort,"T4"] <- subset(sim_data, dv == 1)$time[4]}
      if(!is.na(subset(sim_data, dv == 0)$time[5])) {Simulation[[i]][cohort,"T5"] <- subset(sim_data, dv == 1)$time[5]}
      if(!is.na(subset(sim_data, dv == 0)$time[6])) {Simulation[[i]][cohort,"T6"] <- subset(sim_data, dv == 1)$time[6]}
      
      # Adding PK arguments to the simulated data frame
      sim_data <- transform(sim_data, amt = 0, cmt = ifelse(dv == 1, 10, 11), 
                            week = floor(time / (24 * 7)), mdv = 0, evid = 0, 
                            lamt = -15, tau = 0, addl = 0)[, names(dummy)]
      sim_data <- arrange(rbind(subset(dummy, cmt %in% c(1)), sim_data), id, time, cmt)
      sim_data <- rbind(previous_sim_data, sim_data)      
      
      
      
      # ----------------------------------------------------
      #  Auxiliary steps
      #  1) store true probability of current dose
      #  2) create list object "curent data" which
      #  contains the current data
      # ----------------------------------------------------
      Simulation[[i]][cohort,"true.prob"] <- TrueProbs[which(colnames(Fixed$AUC)==currentDoselevel)]
      
      
      current.data <- list(
        Npat = Simulation[[i]][,"Npat"][1:cohort],
        Ntox = Simulation[[i]][,"Ntox"][1:cohort],
        DosesAdm = Simulation[[i]][,"DoseAdm"][1:cohort],
        FreqAdm = Simulation[[i]][,"FreqAdm"][1:cohort]
      )
      
      # current data to be used in Stan
      current.data.stan  <- c(sim_data, 
                              list(N = nrow(sim_data),
                                   tref_month = tref_month,
                                   Nregimens = nrow(Fixed$SCHEDULES),
                                   regimens_lamt = Fixed$SCHEDULES[, 4],
                                   regimens_tau = Fixed$SCHEDULES[, 8],
                                   addls = Fixed$SCHEDULES[, 7],
                                   params_prior =  Prior_titepk,                              
                                   theta =  params$theta,
                                   ref_dose = Fixed$DOSEREF,
                                   ref_tau = Fixed$REF_TAU))
      
      
      currentRun <- sampling(model_titepk,
                             current.data.stan,
                             cores = 1,
                             init = 1,
                             chains = 4,
                             seed = StanSeed,
                             refresh = 0,
                             iter = niter,
                             control = list(adapt_delta=0.975),
                             open_progress = FALSE)
      
      
      # ----------------------------------------------------
      #  1) get next dose and if MTD reached from
      #  decision rule
      #  2) store different values of interest
      # ----------------------------------------------------
      current_prefs_titepk <- matrix(summary(currentRun)$summary[c("P_dose[1]", "P_dose[2]", "P_dose[3]", "P_dose[4]",
                                                                   "P_dose[5]", "P_dose[6]", "P_dose[7]", "P_dose[8]",
                                                                   "P_dose[9]", "P_dose[10]", "P_dose[11]", "P_dose[12]"), 
                                                                 c(4, 6, 8)], ncol = 3)
      row.names(current_prefs_titepk) <- c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11", "L12")
      current_pcats_titepk <- matrix(summary(currentRun)$summary[c("P_cat[1,1]","P_cat[1,2]","P_cat[1,3]","P_cat[2,1]","P_cat[2,2]","P_cat[2,3]",
                                                                   "P_cat[3,1]","P_cat[3,2]","P_cat[3,3]","P_cat[4,1]","P_cat[4,2]","P_cat[4,3]",
                                                                   "P_cat[5,1]","P_cat[5,2]","P_cat[5,3]","P_cat[6,1]","P_cat[6,2]","P_cat[6,3]",
                                                                   "P_cat[7,1]","P_cat[7,2]","P_cat[7,3]","P_cat[8,1]","P_cat[8,2]","P_cat[8,3]",
                                                                   "P_cat[9,1]","P_cat[9,2]","P_cat[9,3]","P_cat[10,1]","P_cat[10,2]","P_cat[10,3]",
                                                                   "P_cat[11,1]","P_cat[11,2]","P_cat[11,3]","P_cat[12,1]","P_cat[12,2]","P_cat[12,3]"),
                                                                 1], 
                                     ncol = 3, nrow = ncol(Fixed$AUC), byrow = TRUE)
      row.names(current_pcats_titepk) <- c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11", "L12")
      # estimated logbeta
      
      Simulation[[i]][cohort,"Rhat.max"] <- max(summary(currentRun)$summary[, "Rhat"], na.rm = TRUE)
      
      Simulation[[i]][cohort,"under"]  <- current_pcats_titepk[, 1][which(colnames(Fixed$AUC)==currentDoselevel)]
      Simulation[[i]][cohort,"target"] <- current_pcats_titepk[, 2][which(colnames(Fixed$AUC)==currentDoselevel)]
      Simulation[[i]][cohort,"over"]   <- current_pcats_titepk[, 3][which(colnames(Fixed$AUC)==currentDoselevel)]
      
      # ----------------------------------------------------
      # create datsummary object
      # suppress output via sink
      # ----------------------------------------------------
      decision <- TitePKSimDecisionRule(Fixed = Fixed, 
                                        datsummary = data.frame(current.data), 
                                        Ntot = Ntot, 
                                        mnMTD = mnMTD, 
                                        mnTOT = mnTOT, 
                                        mnTAR = mnTAR,
                                        Pall = cbind(current_prefs_titepk, current_pcats_titepk),
                                        Decision = Decision, 
                                        Escalation = Escalation,
                                        DoseAdm = currentDoseAdm,
                                        DoselevelAdm = currentDoselevel,
                                        FreqAdm = currentFreqAdm,
                                        Pcutoffs = Pcutoffs,
                                        Crit = Crit)
      
      too.toxic <- decision$too.toxic
      is.MTD <- decision$is.MTD
      nextDoseAdm <- decision$DoseAdm
      nextFreqAdm <- decision$FreqAdm
      nextDoselevel <- decision$DoselevelAdm
      
      if(!too.toxic){
        
        
        
        
        Simulation[[i]][cohort,"nextunder"]  <- as.numeric(current_pcats_titepk[, 1][which(colnames(Fixed$AUC)==nextDoselevel)])
        Simulation[[i]][cohort,"nexttarget"] <- as.numeric(current_pcats_titepk[, 2][which(colnames(Fixed$AUC)==nextDoselevel)])
        Simulation[[i]][cohort,"nextover"]   <- as.numeric(current_pcats_titepk[, 3][which(colnames(Fixed$AUC)==nextDoselevel)])
        
        Simulation[[i]][cohort,"nextprob"] <-  TrueProbs[which(colnames(Fixed$AUC)==nextDoselevel)]
        
        Simulation[[i]][cohort,"MTD"] <- is.MTD
        Simulation[[i]][cohort,"too.toxic"] <- too.toxic
        
        Simulation[[i]][cohort,"nextDoseAdm"] <- nextDoseAdm
        Simulation[[i]][cohort,"nextFreqAdm"] <- nextFreqAdm
        
        
      } else {
        
        Simulation[[i]][cohort,"nextunder"] <- NA
        Simulation[[i]][cohort,"nexttarget"] <- NA
        Simulation[[i]][cohort,"nextover"] <- NA
        Simulation[[i]][cohort,"MTD"] <- is.MTD
        Simulation[[i]][cohort,"too.toxic"] <- too.toxic
        Simulation[[i]][cohort,"nextprob"] <- NA
        
        Simulation[[i]][cohort,"nextDoseAdm"] <- nextDoseAdm
        Simulation[[i]][cohort,"nextFreqAdm"] <- nextFreqAdm
        
      }
      
      
      currentDoseAdm <- nextDoseAdm
      currentFreqAdm <- nextFreqAdm
      currentDoselevel <- nextDoselevel
      
      previous_sim_data <- sim_data
      
      
      cat(paste("Cohort: ", cohort,"\n",sep=""))
      
      cohort <- cohort + 1 
      
    }
    
  }
  
  cat("\n All simulations successfully done! ")
  cat(paste(Sys.time()))
  
  outlist <- list("Fixed" = Fixed,
                  "nTrials" = nTrials,
                  "TrueProbs" = TrueProbs,
                  "Pcutoffs" = Pcutoffs,
                  "td.pos" = td.pos,
                  "out" = Simulation)
  
  save(outlist, file = saveRes)
  
  return(outlist)
} 

