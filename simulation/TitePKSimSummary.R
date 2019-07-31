#########################################################################
# 
# Code to summarize OC's for Log2BayesSim
# Developer: Simon Wandel
# Adapted by: Burak Gunhan
#########################################################################
TitePKSimSummary <- function(Simulation = NULL, 
                             Fixed = NULL, 
                             Pcutoffs = NULL, 
                             TrueProbs = NULL,
                             td.pos = NULL, 
                             Outfile = NULL,
                             saveRes = NULL,
                             Ntrials = NULL){
  
  
  # ----------------------------------------------------
  # create matrices 
  # ----------------------------------------------------
  mtds        <- matrix(NA, nrow = 3, ncol = length(colnames(Fixed$AUC)) + 2)
  mtds.region <- matrix(NA, nrow = 1, ncol = 5)  
  avg.dlt     <- matrix(NA, nrow = 1, ncol = 3)
  avg.pat     <- matrix(NA, nrow = 1, ncol = 3) 
  stdev.pat   <- matrix(NA, nrow = 1, ncol = 1)
  # -----------------------------------------------
  # Additional summary
  #
  # Summary table of the number of patients enrolled by dose level
  # Summary table of the number of DLTs by dose level
  # -----------------------------------------------
  avg.dlt.dose <- matrix(NA, nrow = 1, ncol =  length(colnames(Fixed$AUC)))
  avg.pat.dose <- matrix(NA, nrow = 1, ncol =  length(colnames(Fixed$AUC))) 
  
  # ----------------------------------------------------
  # create the under, target, over doses
  # ----------------------------------------------------
  under.dose.regions  <- TrueProbs <  Pcutoffs[td.pos-1]
  target.dose.regions <- TrueProbs >= Pcutoffs[td.pos-1] & TrueProbs <= Pcutoffs[td.pos]
  over.dose.regions   <- TrueProbs > Pcutoffs[td.pos]
  
  # -----------------------------------------------
  # derive the usual statistics
  # Average proportion of patients, DLT in target dose
  # Average proportion of patients, DLT in over dose
  # Average proportion of patients, DLT in under dose
  # -----------------------------------------------
  
  
  avg.pat[,1] <- mean(sapply(Simulation, FUN = function(x){
    sum(x[x[,"Doselevel"] %in% colnames(Fixed$AUC)[under.dose.regions],"Npat"])
  }
  ))
  
  avg.pat[,2] <- mean(sapply(Simulation, FUN = function(x){
    sum(x[x[,"Doselevel"] %in% colnames(Fixed$AUC)[target.dose.regions],"Npat"])}
  ))
  
  avg.pat[,3] <- mean(sapply(Simulation, FUN = function(x){
    sum(x[x[,"Doselevel"] %in% colnames(Fixed$AUC)[over.dose.regions],"Npat"])}
  ))

  ## Variability patients overall
  stdev.pat <-  round(sqrt(var(sapply(Simulation, FUN = function(x){
    sum(x[x[,"Doselevel"] %in% colnames(Fixed$AUC),"Npat"])
  }
  ))), 1)
  
  
  
  avg.dlt[,1] <- mean(sapply(Simulation, FUN = function(x){
    sum(x[x[,"Doselevel"] %in% colnames(Fixed$AUC)[under.dose.regions],"Ntox"])}
  ))
  
  avg.dlt[,2] <- mean(sapply(Simulation, FUN = function(x){
    sum(x[x[,"Doselevel"] %in% colnames(Fixed$AUC)[target.dose.regions],"Ntox"])}
  ))
  
  avg.dlt[,3] <- mean(sapply(Simulation, FUN = function(x){
    sum(x[x[,"Doselevel"] %in% colnames(Fixed$AUC)[over.dose.regions],"Ntox"])}
  ))
  
  mtds[1,1:length(colnames(Fixed$AUC))] <- TrueProbs
  
  mtd.doses <- sapply(Simulation, FUN = function(x){x[which(x[,"MTD"] == 1),"Doselevel"]})
  
  for(m in 1:length(colnames(Fixed$AUC))){
    
    mtds[2, m] <- sum(mtd.doses == colnames(Fixed$AUC)[m], na.rm = TRUE)
    mtds[3, m] <- 100*mtds[2, m] / length(Simulation)
    
  }
  
  mtds[2, m + 1] <- sum(sapply(Simulation, FUN = function(x){sum(x[!is.na(x[,"too.toxic"]),"too.toxic"] == 1)})) 
  mtds[2, m + 2] <- length(Simulation) - sum(mtds[2, 1:(m + 1)])
  
  mtds[3, m + 1] <- 100*mtds[2, m + 1] / length(Simulation)
  mtds[3, m + 2] <- 100*mtds[2, m + 2] / length(Simulation)
  mtds.region[, 1] <- sum(mtds[2, which(under.dose.regions)])
  mtds.region[, 2] <- sum(mtds[2, which(target.dose.regions)])
  mtds.region[, 3] <- sum(mtds[2, which(over.dose.regions)])
  mtds.region[, 4] <- mtds[2 ,m + 1]   
  mtds.region[, 5] <- length(Simulation) - sum(mtds.region[, 1:4]) 
  
  colnames(mtds) <- c(Fixed$AUC["Dose",], "stopped (tox)", "stopped (max N), or never started")
  rownames(mtds) <- c("TrueProb", "Ntrials", "%") 
  
  colnames(mtds.region) <- c("underdose", "targetdose", "overdose", "stopped (tox)", "stopped (max N), or never started")
  
  mtds.region.perc <- 100*mtds.region/length(Simulation)
  
  colnames(avg.pat) <- c("underdose", "targetdose", "overdose")
  colnames(avg.dlt) <- c("underdose", "targetdose", "overdose")
  
  
  # -----------------------------------------------
  # Additional summary
  #
  # Summary table of the number of patients enrolled by dose level
  # Summary table of the variability of patients enroolled
  # Summary table of the number of DLTs by dose level
  # -----------------------------------------------
  for (j in 1:length(colnames(Fixed$AUC))){
    avg.pat.dose[1,j] <- mean(sapply(Simulation, FUN = function(x){
      sum(x[x[,"Doselevel"] %in% colnames(Fixed$AUC)[j], "Npat"])}
    ))
    
    avg.dlt.dose[1,j] <- mean(sapply(Simulation, FUN = function(x){
      sum(x[x[,"Doselevel"] %in% colnames(Fixed$AUC)[j], "Ntox"])}
    ))
    
  }
  colnames(avg.pat.dose) <- colnames(Fixed$AUC)
  colnames(avg.dlt.dose) <- colnames(Fixed$AUC)
  
  outlist <- list("avg.pat" = avg.pat, 
                  "stdev.pat" = stdev.pat,
                  "avg.dlt" = avg.dlt, 
                  "mtds.collapsed" = mtds.region, 
                  "mtds.collapsed.perc" = mtds.region.perc,
                  "avg.pat.dose" = avg.pat.dose,
                  "avg.dlt.dose" = avg.dlt.dose,
                  "mtds.dose" = mtds)
  
  save(outlist, file = saveRes)
  
  print(outlist)
  
}


