## defining the dosing history for the refernce regimen
regimen_ref <- data.frame(cmt = 1, time = seq(0, 365*2*24, by = tau_drug), amt = ref_dose, lamt = log(ref_dose), 
                          mdv = 1, evid = 1)
regimen_ref_packed <- arrange(ddply(regimen_ref, .(label_dose_periods(cmt, amt)), pack_dosing), time, cmt)
regimen_ref_packed[,1] <- NULL
