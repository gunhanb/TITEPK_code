TitePKSimDecisionRule <- function(Fixed = NULL,
                                  datsummary = NULL,
                                  Ntot = NULL,
                                  mnMTD = NULL,
                                  mnTOT = NULL,
                                  mnTAR= NULL,
                                  Pall = NULL,
                                  Decision = 1, 
                                  Escalation = "maximal",
                                  DoseAdm = NULL,
                                  FreqAdm = NULL,
                                  DoselevelAdm = NULL,
                                  Pcutoffs = NULL,
                                  Crit = NULL){
  
  
  # ---------------------------------------------
  # Escalation rule: either maximal or
  # target
  # ---------------------------------------------
  escalation.in <- toupper(Escalation)
  if(!escalation.in %in% c("MAXIMAL", "TARGET")){
    stop("Argument Escalation can only take values maximal or target.")
  }
  
  
  # ---------------------------------------------
  # CDoseAdm: current administered dose
  # if Pcutoffs / Crit NULL: set to default
  # ---------------------------------------------
  if(is.null(Pcutoffs)){
    Pcutoffs <- c(0.16, 0.33)
    Crit <- c(1, 1, 0.25)
  }
  
  od.pos <- length(Pcutoffs) + 1
  td.pos <- od.pos - 1
  
  
  # ---------------------------------------------
  # Preparation steps:
  # 1) collecting number of patients in current dose level
  #    datsummary is not collapsed by dose level
  # 2) define eligible doses (those respecting EWOC)
  # 3) define possible doses (those with <= 100%
  #    escalation)
  # ---------------------------------------------
  totNPAT <- sum(datsummary[,"Npat"][datsummary[,"DosesAdm"] == DoseAdm & datsummary[,"FreqAdm"] == FreqAdm])
  
  odprob <-  Pall[, 6]
  tgprob <-  Pall[, 5]
  names(odprob) = names(tgprob) = c("L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "L11", "L12")
  
  dose.eligible <- odprob < Crit[length(Crit)]
  #possible.doses <- (Fixed$DOSES/DoseAdm) <= 3.5
  #possible.doses <- (Fixed$DOSES/DoseAdm) <= 2
  ## Using AUCs instead of dose levels
  #possible.doses = aucs["AUC",]/aucs["AUC",aucs["Dose",] == DoseAdm & aucs["Freq",] == FreqAdm]  <= 2
  possible.doses = aucs["AUC",]
  
  # ---------------------------------------------
  # if all too toxic: stop
  # ---------------------------------------------
  if(sum(dose.eligible) == 0){
    too.toxic     <- TRUE
    is.MTD        <- FALSE
    nextDose      <- NA
    nextFreq      <- NA
    nextDoselevel <- NA
  }else{
    too.toxic <- FALSE
    
    #-----------------------------------------------------------------        
    # if escalation = maximal: next = highest possible dose
    # if escalation = target: next = dose maximizing target interval
    #----------------------------------------------------------------        
    if(escalation.in == "MAXIMAL") {
      if(length(as.numeric(Fixed$AUC[1,dose.eligible & possible.doses])) == 1) {
        if(Fixed$AUC[2,dose.eligible & possible.doses] == 24 * 1) {
          nextDoselevel       = "L10"
          nextDose            = 8
          nextFreq            = 24 * 1
        }
        if(Fixed$AUC[2,dose.eligible & possible.doses] == 24 * 2) {
          nextDoselevel       = "L7"
          nextDose            = 8
          nextFreq            = 24 * 2
        }
        
        if(Fixed$AUC[2,dose.eligible & possible.doses] == 24 * 4) {
          nextDoselevel       = "L4"
          nextDose            = 8
          nextFreq            = 24 * 4
        }
        if(Fixed$AUC[2,dose.eligible & possible.doses] == 24 * 8) {
          nextDoselevel       = "L1"
          nextDose            = 8
          nextFreq            = 24 * 8
        }
        
      }
      
      else {
        selected.doseLevels = colnames(Fixed$AUC[1,dose.eligible & possible.doses])
        selected.aucs       = aucs["AUC", colnames(aucs) %in% selected.doseLevels]
        nextDoselevel       = colnames(selected.aucs)[apply(selected.aucs,1,which.max)]
        nextDose            = aucs["Dose", nextDoselevel]
        nextFreq            = aucs["Freq", nextDoselevel]
      }
    }
    next.equal <- nextDoselevel == DoselevelAdm
    
    # ---------------------------------------------
    # check if current dose level = MTD
    # ---------------------------------------------
    
    mnMTD.check <- (totNPAT >= mnMTD)
    MTD.rule    <- (Ntot >= mnTOT || tgprob[which(colnames(aucs) == DoselevelAdm)] >= mnTAR)
    #MTD.rule    <- (Ntot >= mnTOT)
    
    is.MTD      <- (mnMTD.check & MTD.rule & next.equal)
  }
  
  # ---------------------------------------------
  # Output
  # ---------------------------------------------
  
  return(list("too.toxic"= too.toxic,
              "is.MTD"   = is.MTD,
              "DoseAdm" = nextDose,
              "FreqAdm" = nextFreq,
              "DoselevelAdm" = nextDoselevel))
  
}

