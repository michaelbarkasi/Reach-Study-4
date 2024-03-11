
# Preprocessing & Data Cleaning:
# Functions to process data into dataframes we can use for the analysis

sample_slide_limit <- -100 # For fine-tune adjustments of synchronization between OT and SU data

# Dependency function: 
#   ParticipantChopOT is used in chop_data. It partitions ("chops") raw OptiTrack data into discrete reaches using Son Unit data;
#   function cannot be run on its own; must be run inside chop_data
ParticipantChopOT <- function(index_) {
  
  model <- pTable$model[index_]
  NOWtime <- pTable$SUStartTime[index_]
  ppp <- pTable$number[index_]
  
  float_cutoff <- float__distance
  
  noisy_print <- TRUE

  # Begin chopping
  if (noisy_print) cat("P: ", ppp) 
  psv <- c( ppp ) # add participant number to psv vector
  
  # SU = Sonification system, OT = OptiTrack data
  
  SUreaches <- motionSU$reach[which(motionSU$participant==ppp)]
  SUtimes <- motionSU$time[which(motionSU$participant==ppp)] - NOWtime
  
  OTtimes <- motionOT$time[which(motionOT$participant==ppp)]
  tx <- motionOT$x[which(motionOT$participant==ppp)]
  ty <- motionOT$y[which(motionOT$participant==ppp)]
  tz <- motionOT$z[which(motionOT$participant==ppp)]
  
  # Code should check all steps for unexpected values
  if ( any(is.na(SUtimes)) ) {
    cat("\n WARNING! Unexpected NA in SUtimes, participant:", ppp,"\n")
  }
  
  # Go through SU data (motionSU) and find the rows where a new reach starts
  SUstarts <- 1
  for ( i in 1:(length(SUreaches)-1) ) {
    if (is.na(SUreaches[i]) || is.na(SUreaches[i+1]) ) {
      cat("\n WARNING! Unexpected NA in SUreaches, participant:", ppp, "index:", i,"\n")
    }
    if ( SUreaches[i] != SUreaches[i+1] ) {
      SUstarts <- c( SUstarts, i+1 )
    }
  }
  SUnumofreaches <- length(SUstarts)
  if ( length(unique(SUreaches)) != SUnumofreaches || 
       SUnumofreaches != max(SUreaches,na.rm=TRUE) ) {
    cat("\n WARNING! Something went wrong counting reaches for particippant:", ppp,"\n")
  }
  if (noisy_print) cat(", R: ", SUnumofreaches) 
  psv <- c( psv, SUnumofreaches )
  
  # Go through SU data (motionSU) and find the rows where each reach ends
  SUends <- rep(NA,SUnumofreaches)
  for ( i in 1:(length(SUends)-1) ) {
    SUends[i] <- SUstarts[i+1] - 1
  }
  SUends[length(SUends)] <- length(SUreaches)
  if ( any(is.na(SUends)) ) {
    cat("\n WARNING! Unexpected NA in SUends, participant:", ppp,"\n")
  }
  
  # Find velocities in SU data
  SUqax <- motionSU$qax[which(motionSU$participant==ppp)]
  SUqay <- motionSU$qay[which(motionSU$participant==ppp)]
  SUqaz <- motionSU$qaz[which(motionSU$participant==ppp)]
  SUqar <- motionSU$qar[which(motionSU$participant==ppp)]
  SUqbx <- motionSU$qbx[which(motionSU$participant==ppp)]
  SUqby <- motionSU$qby[which(motionSU$participant==ppp)]
  SUqbz <- motionSU$qbz[which(motionSU$participant==ppp)]
  SUqbr <- motionSU$qbr[which(motionSU$participant==ppp)]
  SUvelocity <- rep(NA,length(SUreaches))
  for ( i in 1:SUnumofreaches ) {
    for ( j in (SUstarts[i]+1):SUends[i] ) {
      pa2 <- c( SUqax[j], SUqay[j], SUqaz[j], SUqar[j] )
      pb2 <- c( SUqbx[j], SUqby[j], SUqbz[j], SUqbr[j] )
      pa1 <- c( SUqax[j-1], SUqay[j-1], SUqaz[j-1], SUqar[j-1] )
      pb1 <- c( SUqbx[j-1], SUqby[j-1], SUqbz[j-1], SUqbr[j-1] )
      SUvelocity[j] <- (vdist(pa1,pa2) + vdist(pb1,pb2)) * 1000 # put into qu/s
      SUvelocity[j] <- round(SUvelocity[j],significant_digitsSU)
    }
    firstvel <- SUvelocity[SUstarts[i]+1] - abs(SUvelocity[SUstarts[i]+2] - SUvelocity[SUstarts[i]+1])
    if ( is.na(firstvel) ) {
      SUvelocity[SUstarts[i]] <- NA
    } else if ( firstvel > 0 ) {
      SUvelocity[SUstarts[i]] <- round(firstvel,significant_digitsSU)
    } else {
      SUvelocity[SUstarts[i]] <- 0
    }
  }
  
  # check for bad samples in SU data.
  bad_SU_samples <- NA
  if ( any(is.na(SUqax)) || any(!is.numeric(SUqax)) ) {
    problem <- which(is.na(SUqax) | !is.numeric(SUqax) )
    bad_SU_samples <- c( bad_SU_samples, problem )
  }
  if ( any(is.na(SUqay)) || any(!is.numeric(SUqay)) ) {
    problem <- which(is.na(SUqay) | !is.numeric(SUqay) )
    bad_SU_samples <- c( bad_SU_samples, problem )
  }
  if ( any(is.na(SUqaz)) || any(!is.numeric(SUqaz)) ) {
    problem <- which(is.na(SUqaz) | !is.numeric(SUqaz) )
    bad_SU_samples <- c( bad_SU_samples, problem )
  }
  if ( any(is.na(SUqar)) || any(!is.numeric(SUqar)) ) {
    problem <- which(is.na(SUqar) | !is.numeric(SUqar) )
    bad_SU_samples <- c( bad_SU_samples, problem )
  }
  if ( any(is.na(SUqbx)) || any(!is.numeric(SUqbx)) ) {
    problem <- which(is.na(SUqbx) | !is.numeric(SUqbx) )
    bad_SU_samples <- c( bad_SU_samples, problem )
  }
  if ( any(is.na(SUqby)) || any(!is.numeric(SUqby)) ) {
    problem <- which(is.na(SUqby) | !is.numeric(SUqby) )
    bad_SU_samples <- c( bad_SU_samples, problem )
  }
  if ( any(is.na(SUqbz)) || any(!is.numeric(SUqbz)) ) {
    problem <- which(is.na(SUqbz) | !is.numeric(SUqbz) )
    bad_SU_samples <- c( bad_SU_samples, problem )
  }
  if ( any(is.na(SUqbr)) || any(!is.numeric(SUqbr)) ) {
    problem <- which(is.na(SUqbr) | !is.numeric(SUqbr) )
    bad_SU_samples <- c( bad_SU_samples, problem )
  }
  if ( length(bad_SU_samples) > 1 ) {
    bad_SU_samples <- bad_SU_samples[2:length(bad_SU_samples)]
    bad_SU_samples <- unique(bad_SU_samples)
    if (noisy_print) cat(", BSUSs: ", length(bad_SU_samples))
    psv <- c( psv, length(bad_SU_samples) )
  } else {
    if (noisy_print) cat(", BSUSs: ", 0)
    psv <- c( psv, 0 )
  }
  
  # This vector gives rows where we expect bad SU velocity samples
  bad_SU_samplesV <- c(bad_SU_samples, bad_SU_samples+1)
  # Now check to make sure there are no other bad SU velocity samples
  if ( any(is.na(SUvelocity)) ) {
    problem <- which(is.na(SUvelocity)) 
    if ( !all(problem %in% bad_SU_samplesV) ) {
      cat("\n WARNING! Unexpected NA in SUvelocity, participant:", ppp,"samples:",problem,"known bad SU samples:",bad_SU_samples,"\n")
    }
  }
  if ( any(!is.numeric(SUvelocity)) ) {
    problem <- which(!(is.numeric(SUvelocity)))
    if ( !all(problem %in% bad_SU_samplesV) ) {
      cat("\n WARNING! At least one nonnumeric entry in SUvelocity, participant:", ppp,"sample:",problem,"known bad SU samples:",bad_SU_samples,"\n")
    }
  }
  
  # First, smooth SU velocity data w/ 50ms rolling average
  SUvelocity <- rolling_mean(SUvelocity,50,significant_digitsSU)
  # Save velocity to dataframe; This puts SU velocity in unit quats / ms
  motionSU$RotVelocity[which(motionSU$participant==ppp)] <<- SUvelocity
  
  # Compute lengths of reaches in SU data (rows)
  SUlengths <- rep(NA,SUnumofreaches)
  for ( i in 1:length(SUends) ) {
    SUlengths[i] <- (SUends[i] - SUstarts[i]) + 1
  }
  if ( any(is.na(SUlengths)) ) {
    cat("\n WARNING! Unexpected NA in SUlengths, participant:", ppp,"\n")
  }
  
  # Find and fix overruns in SU data.
  overrun_count <- 0
  overrun_slides <- rep(0, length(SUlengths))
  if ( length(SUlengths) != length(SUends) ) {
    cat("\n WARNING! about to check for overruns, but lengths of SUlengths and SUends doesn't match, participant:", ppp,"\n")
  }
  model_overran <- FALSE
  for ( i in 1:length(SUlengths) ) {
    if ( SUlengths[i] >= overrun_length ) {
      # User most likely never paused, so, we're looking for the inflection point where they turned around.
      # Even if user really did overrun just from going very slow or far, min velocity should still be near the end.
      if ( i != model ) {
        oldend <- SUends[i]
        SUends[i] <- (SUstarts[i]+499) + which.min(SUvelocity[(SUstarts[i]+500):(SUends[i]-5)]) # Note: which.min is insensitive to NAs.
        if ( !is.numeric(SUends[i]) || is.na(SUends[i]) ) {
          cat("\n WARNING! attempting to fix overrun, but new SUends is not numeric or is NA, participant:", ppp, "reach:", i,"\n")
        }
        SUlengths[i] <- (SUends[i] - SUstarts[i]) + 1
        if ( !is.numeric(SUlengths[i]) || is.na(SUlengths[i]) ) {
          cat("\n WARNING! attempting to fix overrun, but new SUlengths is not numeric or is NA, participant:", ppp, "reach:", i,"\n")
        }
        overrun_count <- overrun_count+1
        if ( oldend > SUends[i] ) {
          overrun_slides[i] <- length((SUends[i]+1):oldend)
          SUreaches[(SUends[i]+1):oldend] <- NA
        }
      } else {
        model_overran <- TRUE
      }
    }
  }
  if (noisy_print && !model_overran) cat(", ORs: ", overrun_count) 
  if (noisy_print && model_overran) cat(", ORs:* ", overrun_count) 
  psv <- c( psv, overrun_count )
  
  # Estimate OTsampleslide
  OTsampleslide <- -20 # default guess ( = 200ms ) 
  R1fts <- SUtimes[SUstarts[1]] # Not the true reach 1 SU start time
  R1l <- SUlengths[1] / 1000 
  max_initial_ty <- max(ty[1:IPrun],na.rm = TRUE)
  checkindex <- 1
  found_OT_R1_guess <- FALSE
  while ( OTtimes[checkindex] < (R1fts+5) && 
          !found_OT_R1_guess && 
          checkindex < 30000 ) {
    if ( ty[checkindex] > (max_initial_ty+1) ) {
      found_OT_R1_guess <- TRUE
    }
    checkindex <- checkindex + 1
  }
  OT_R1_guess_index <- checkindex - 1
  OT_R1_guess_time <- OTtimes[OT_R1_guess_index]
  slide_time_guess <- ( R1fts - R1l ) - OT_R1_guess_time
  if ( !(is.na(slide_time_guess) || 
         !is.numeric(slide_time_guess) || 
         slide_time_guess < 0) ) {
    OTsampleslide <- as.integer(-slide_time_guess * 100)
  }
  if ( OTsampleslide < sample_slide_limit ) {
    OTsampleslide <- -20 # something went wrong, set back to default
  }
  if (noisy_print) cat(", Sl: ", OTsampleslide) 
  psv <- c( psv, OTsampleslide )
  
  # Find and remove short reaches in SU data.
  short_reaches <- NA
  if ( length(short_reaches) != 1 ) {
    cat("\n WARNING! expected short_reach length of 1 at this point, but found otherwise, participant:", ppp,"\n")
  }
  for ( i in 1:length(SUlengths) ) {
    if ( SUlengths[i] < reach_min_length ) {
      short_reaches <- c( short_reaches, i )
      SUreaches[SUstarts[i]:SUends[i]] <- NA
      SUstarts[i] <- NA
      SUends[i] <- NA
      SUlengths[i] <- NA
    }
  }
  num_of_short_reaches <- 1
  if ( length(short_reaches) == 1 ) {
    num_of_short_reaches <- 0 
  } else {
    short_reaches <- short_reaches[2:length(short_reaches)] 
    num_of_short_reaches <- length(short_reaches)
  }
  if (noisy_print) cat(", Ss: ", num_of_short_reaches) 
  psv <- c( psv, num_of_short_reaches )
  if ( model %in% short_reaches ) {
    cat("\n WARNING! model is in short_reaches, participant:", ppp,"\n")
  }
  
  # Find reach endpoints in OT data (initial guess).
  if ( any(is.na(OTtimes)) ) {
    cat("\n WARNING! Unexpected NA in OTtimes, participant:", ppp,"\n")
  }
  OTends <- rep(NA,SUnumofreaches)
  if ( length(OTends) != length(overrun_slides) ) {
    cat("\n WARNING! about to find OTends, but lengths of OTends and overrun_slides doesn't match, participant:", ppp,"\n")
  }
  for ( i in 1:last_no_feedback_trial_num ) {
    if ( !(i %in% short_reaches) ) {
      tSU <- SUtimes[SUstarts[i]] # the reach ended at the time the print out started
      if ( !is.numeric(tSU) || is.na(tSU) ) {
        cat("\n WARNING! attempting to find OTends, but tSU is not numeric or is NA, participant:", ppp, "reach:", i,"\n")
        cat("\n delayfactor: ", delayfactor, "SUtimes: ", SUtimes[SUstarts[i]], "SUstart: ", SUstarts[i], "SUlengths: ", SUlengths[i], "\n")
      }
      SEARCHING <- TRUE
      for ( j in 1:length(OTtimes) ) {
        if ( SEARCHING == TRUE && OTtimes[j] > tSU ) {
          OTends[i] <- ((j-1) + OTsampleslide) - (overrun_slides[i]/10)
          if ( !is.numeric(OTends[i]) || is.na(OTends[i]) ) {
            cat("\n WARNING! attempting to set OTends, but OTends is not numeric or is NA, participant:", ppp, "reach:", i,"\n")
          }
          SEARCHING <- FALSE
        }
      }
    }
  }
  for ( i in (last_no_feedback_trial_num+1):(pb_start-1) ) {
    if ( !(i %in% short_reaches) ) {
      delayfactor <- 0 
      if ( FBcondition2 == "terminal" & pTable$Condition[which(pTable$number==ppp)] == FBcondition2 ) {
        delayfactor <- ( 0.5 + (SUlengths[model]/1000) ) + (SUlengths[i]/1000)
      } else {
        delayfactor <- 0.5 + (SUlengths[model]/1000) # the reach ended at the time the print out started, minus 0.5 seconds and the time of the model
      }
      tSU <- ( SUtimes[SUstarts[i]] - delayfactor )
      if ( !is.numeric(tSU) || is.na(tSU) ) {
        cat("\n WARNING! attempting to find OTends, but tSU is not numeric or is NA, participant:", ppp, "reach:", i,"\n")
        cat("\n delayfactor: ", delayfactor, "SUtimes: ", SUtimes[SUstarts[i]], "SUstart: ", SUstarts[i], "SUlengths: ", SUlengths[i], "\n")
      }
      SEARCHING <- TRUE
      for ( j in 1:length(OTtimes) ) {
        if ( SEARCHING == TRUE && OTtimes[j] > tSU ) {
          OTends[i] <- ((j-1) + OTsampleslide) - (overrun_slides[i]/10)
          if ( !is.numeric(OTends[i]) || is.na(OTends[i]) ) {
            cat("\n WARNING! attempting to set OTends, but OTends is not numeric or is NA, participant:", ppp, "reach:", i,"\n")
          }
          SEARCHING <- FALSE
        }
      }
    }
  }
  for ( i in pb_start:SUnumofreaches ) {
    if ( !(i %in% short_reaches) ) {
      delayfactor <- 0 
      if ( FBcondition2 == "terminal" & pTable$Condition[which(pTable$number==ppp)] == FBcondition2 ) {
        delayfactor <- 0.5 + (SUlengths[model]/1000) # the reach ended at the time the print out started, minus 0.5 seconds and the time of the model
      } else if ( FBcondition2 == "terminal" & pTable$Condition[which(pTable$number==ppp)] == FBcondition1 ) {
        delayfactor <- ( 0.5 + (SUlengths[model]/1000) ) + (SUlengths[i]/1000)
      } else {
        delayfactor <- 0.5 + (SUlengths[model]/1000)
      }
      tSU <- ( SUtimes[SUstarts[i]] - delayfactor )
      if ( !is.numeric(tSU) || is.na(tSU) ) {
        cat("\n WARNING! attempting to find OTends, but tSU is not numeric or is NA, participant:", ppp, "reach:", i,"\n")
        cat("\n delayfactor: ", delayfactor, "SUtimes: ", SUtimes[SUstarts[i]], "SUstart: ", SUstarts[i], "SUlengths: ", SUlengths[i], "\n")
      }
      SEARCHING <- TRUE
      for ( j in 1:length(OTtimes) ) {
        if ( SEARCHING == TRUE && OTtimes[j] > tSU ) {
          OTends[i] <- ((j-1) + OTsampleslide) - (overrun_slides[i]/10)
          if ( !is.numeric(OTends[i]) || is.na(OTends[i]) ) {
            cat("\n WARNING! attempting to set OTends, but OTends is not numeric or is NA, participant:", ppp, "reach:", i,"\n")
          }
          SEARCHING <- FALSE
        }
      }
    }
  }
  
  # Find reach start points in OT data (initial guess).
  OTstarts <- rep(NA,SUnumofreaches)
  for ( i in 1:SUnumofreaches ) {
    if ( !(i %in% short_reaches) ) {
      starttime <- OTtimes[OTends[i]] - ( SUlengths[i] / 1000 )
      if ( !is.numeric(starttime) || is.na(starttime) ) {
        cat("\n WARNING! attempting to find OTstarts, but starttime is not numeric or is NA, participant:", ppp, "reach:", i,"\n")
        cat("\n OTends: ", OTends[i], "SUlengths: ", SUlengths[i], "\n")
      }
      if (length(starttime) != 1 ) {
        cat("\n WARNING! attempting to find OTstarts, but starttime is not of length 1, participant:", ppp, "reach:", i,"\n")
        cat("\n OTends: ", OTends[i], "SUlengths: ", SUlengths[i], "\n")
      }
      SEARCHING <- TRUE
      for ( j in 1:length(OTtimes) ) {
        if ( SEARCHING == TRUE && OTtimes[j] > starttime ) {
          OTstarts[i] <- j-1
          if ( !is.numeric(OTstarts[i]) || is.na(OTstarts[i]) ) {
            cat("\n WARNING! attempting to set OTstarts, but OTstarts is not numeric or is NA, participant:", ppp, "reach:", i,"\n")
          }
          SEARCHING <- FALSE
        }
      }
    }
  }
  
  # Refine start points (and so endpoints too, length is fixed) in OT data 
  #  by moving them back until doing so more no longer gets us closer to the true initial position.
  avg_tx <- mean(tx[1:IPrun],na.rm = TRUE) # na.rm b/c OT sometimes blinks out and leaves blank rows
  avg_ty <- mean(ty[1:IPrun],na.rm = TRUE)
  avg_tz <- mean(tz[1:IPrun],na.rm = TRUE) 
  IP <- c( avg_tx, avg_ty, avg_tz )
  #  The random reaches are handled separately, as they don't have the model playback affecting SU printout time
  k_vec <- 1
  for ( i in 1:last_no_feedback_trial_num ) {
    if ( !(i %in% short_reaches) ) {
      no_count <- 0
      k <- OTstarts[i]
      SEARCHINGforK <- TRUE
      for ( j in k:(k-100) ) {
        if ( SEARCHINGforK ) {
          if ( no_count < 4 ) {
            rIP <- c( tx[j], ty[j], tz[j] ) 
            rIPback <- c( tx[j-1], ty[j-1], tz[j-1] ) 
            no_closer <- vdist(IP,rIPback) >= vdist(IP,rIP) * 0.995
            if ( isTRUE(no_closer) ) {
              no_count <- no_count+1
            }
          } else {
            k <- j+2
            SEARCHINGforK <- FALSE
          }
        }
      }
      k_vec <- c( k_vec, OTstarts[i] - k )
    }
  }
  k_vec <- k_vec[2:length(k_vec)]
  sync_tune <- as.integer(median(k_vec))
  # There might be NA in OTstarts and OTends, from removed shorts, but that's okay,
  #     since NA - sync_tune is still NA. 
  OTstarts[1:last_no_feedback_trial_num] <- OTstarts[1:last_no_feedback_trial_num] - sync_tune
  OTends[1:last_no_feedback_trial_num] <- OTends[1:last_no_feedback_trial_num] - sync_tune
  if ( is.na(sync_tune) ) {
    cat("\n WARNING! sync_tune for Random is NA, participant:", ppp,"\n")
  }
  if (noisy_print) cat(", S-R: ", sync_tune) 
  psv <- c( psv, sync_tune )
  # Now do reaches w/ model playback (feedback): 
  k_vec <- 1
  for ( i in (last_no_feedback_trial_num+1):(pb_start-1) ) {
    if ( !(i %in% short_reaches) ) {
      no_count <- 0
      k <- OTstarts[i]
      SEARCHINGforK <- TRUE
      for ( j in k:(k-100) ) {
        if ( SEARCHINGforK ) {
          if ( no_count < 4 ) {
            rIP <- c( tx[j], ty[j], tz[j] ) 
            rIPback <- c( tx[j-1], ty[j-1], tz[j-1] ) 
            no_closer <- vdist(IP,rIPback) >= vdist(IP,rIP) * 0.995
            if ( isTRUE(no_closer) ) {
              no_count <- no_count+1
            }
          } else {
            k <- j+2
            SEARCHINGforK <- FALSE
          }
        }
      }
      k_vec <- c( k_vec, OTstarts[i] - k )
    }
  }
  k_vec <- k_vec[2:length(k_vec)]
  sync_tune <- as.integer(median(k_vec))
  # There might be NA in OTstarts and OTends, from removed shorts, but that's okay,
  #     since NA - sync_tune is still NA. 
  OTstarts[(last_no_feedback_trial_num+1):(pb_start-1)] <- OTstarts[(last_no_feedback_trial_num+1):(pb_start-1)] - sync_tune
  OTends[(last_no_feedback_trial_num+1):(pb_start-1)] <- OTends[(last_no_feedback_trial_num+1):(pb_start-1)] - sync_tune
  if ( is.na(sync_tune) ) {
    cat("\n WARNING! sync_tune for Feedback is NA, participant:", ppp,"\n")
  }
  if (noisy_print) cat(", S-F: ", sync_tune) 
  psv <- c( psv, sync_tune )
  # Now do reaches w/ model playback (post-feedback block): 
  k_vec <- 1
  for ( i in pb_start:SUnumofreaches ) {
    if ( !(i %in% short_reaches) ) {
      no_count <- 0
      k <- OTstarts[i]
      SEARCHINGforK <- TRUE
      for ( j in k:(k-100) ) {
        if ( SEARCHINGforK ) {
          if ( no_count < 4 ) {
            rIP <- c( tx[j], ty[j], tz[j] ) 
            rIPback <- c( tx[j-1], ty[j-1], tz[j-1] ) 
            no_closer <- vdist(IP,rIPback) >= vdist(IP,rIP) * 0.995
            if ( isTRUE(no_closer) ) {
              no_count <- no_count+1
            }
          } else {
            k <- j+2
            SEARCHINGforK <- FALSE
          }
        }
      }
      k_vec <- c( k_vec, OTstarts[i] - k )
    }
  }
  k_vec <- k_vec[2:length(k_vec)]
  sync_tune <- as.integer(median(k_vec))
  # There might be NA in OTstarts and OTends, from removed shorts, but that's okay,
  #     since NA - sync_tune is still NA. 
  OTstarts[pb_start:SUnumofreaches] <- OTstarts[pb_start:SUnumofreaches] - sync_tune
  OTends[pb_start:SUnumofreaches] <- OTends[pb_start:SUnumofreaches] - sync_tune
  if ( is.na(sync_tune) ) {
    cat("\n WARNING! sync_tune for Feedback is NA, participant:", ppp,"\n")
  }
  if (noisy_print) cat(", S-PB: ", sync_tune) 
  psv <- c( psv, sync_tune )
  
  # Find and save velocities in OT data.
  OTvelocity <- rep(NA,length(OTtimes))
  for ( i in 1:SUnumofreaches ) {
    if ( !(i %in% short_reaches) ) {
      OTvelocity[OTstarts[i]] <- 0
      for ( j in (OTstarts[i]+1):OTends[i] ) {
        p2 <- c( tx[j], ty[j], tz[j] )
        p1 <- c( tx[j-1], ty[j-1], tz[j-1] )
        if ( !any(is.na(p1)) && !any(is.na(p2)) ) {
          OTvelocity[j] <- vdist(p1,p2) # this is naturally in cm/cs, which is the same as m/s
          OTvelocity[j] <- round(OTvelocity[j],significant_digitsOT)
          if ( !is.numeric(OTvelocity[j]) || is.na(OTvelocity[j]) || length(OTvelocity[j]) != 1 ) {
            cat("\n WARNING! OTvelocity is not numeric, is NA, or isn't length 1, participant:", ppp, "reach:", i,"sample:", j,"\n")
          }
        }
      }
      firstvel <- OTvelocity[OTstarts[i]+1] - abs(OTvelocity[OTstarts[i]+2] - OTvelocity[OTstarts[i]+1])
      if ( is.na(firstvel) ) {
        OTvelocity[OTstarts[i]] <- NA
      } else if ( firstvel > 0 ) {
        OTvelocity[OTstarts[i]] <- round(firstvel,significant_digitsOT)
      } else {
        OTvelocity[OTstarts[i]] <- 0
      }
    }
  }
  
  # First, smooth OT velocity data w/ 50ms rolling average
  OTvelocity <- rolling_mean(OTvelocity,5,significant_digitsOT)
  # OT velocity is in m/s
  motionOT$SpatialVelocity[which(motionOT$participant==ppp)] <<- OTvelocity
  
  # One more sync check using velocity and adjusting the end position
  end_point_slide <- 0.25 # don't slide below this percentage from the end sample when trying to estimate true reach end point.
  # Start w/ first 25 reaches
  k_vec <- 1
  for ( i in 1:last_no_feedback_trial_num ) {
    if ( !(i %in% short_reaches) ) {
      rlength <- OTends[i] - OTstarts[i]
      end_velocities <- OTvelocity[(OTends[i]-(end_point_slide*rlength)):OTends[i]]
      k <- OTends[i]
      if ( any(!is.na(end_velocities)) ) {
        min_Evel <- min(OTvelocity[(OTends[i]-(end_point_slide*rlength)):OTends[i]],na.rm = TRUE)
        Evel <- OTvelocity[OTends[i]]
        ii <- OTends[i]
        SEARCHINGforII <- TRUE
        while ( is.na(Evel) && SEARCHINGforII ) {
          ii <- ii-1
          Evel <- OTvelocity[ii]
          if ( ii < (OTends[i] - (end_point_slide*rlength)) ) SEARCHINGforII <- FALSE
        }
        if ( !is.na(Evel) && (Evel > (min_Evel*1.02)) ) {
          k <- ii
          SEARCHINGforK <- TRUE
          while(SEARCHINGforK) {
            k <- k-1
            if ( !is.na(OTvelocity[k]) && (OTvelocity[k] <= (min_Evel*1.02)) ) SEARCHINGforK <- FALSE
            if ( k < (OTends[i] - (end_point_slide*rlength)) ) SEARCHINGforK <- FALSE
          }
        } else {
          k <- OTends[i]
        }
      } 
      k_vec <- c( k_vec, OTends[i] - k )
    }
  }
  k_vec <- k_vec[2:length(k_vec)]
  sync_tune <- as.integer(median(k_vec))
  # There might be NA in OTstarts and OTends, from removed shorts, but that's okay,
  #     since NA - sync_tune is still NA. 
  OTstarts[1:last_no_feedback_trial_num] <- OTstarts[1:last_no_feedback_trial_num] - sync_tune
  OTends[1:last_no_feedback_trial_num] <- OTends[1:last_no_feedback_trial_num] - sync_tune
  if ( is.na(sync_tune) ) {
    cat("\n WARNING! end sync_tune for Random is NA, participant:", ppp,"\n")
  }
  if (noisy_print) cat(", ES-R: ", sync_tune) 
  psv <- c( psv, sync_tune )
  # Now the reaches w/ model playback (feedback)
  k_vec <- 1
  for ( i in (last_no_feedback_trial_num+1):(pb_start-1) ) {
    if ( !(i %in% short_reaches) ) {
      rlength <- OTends[i] - OTstarts[i]
      end_velocities <- OTvelocity[(OTends[i]-(end_point_slide*rlength)):OTends[i]]
      k <- OTends[i]
      if ( any(!is.na(end_velocities)) ) {
        min_Evel <- min(OTvelocity[(OTends[i]-(end_point_slide*rlength)):OTends[i]],na.rm = TRUE)
        Evel <- OTvelocity[OTends[i]]
        ii <- OTends[i]
        SEARCHINGforII <- TRUE
        while ( is.na(Evel) && SEARCHINGforII ) {
          ii <- ii-1
          Evel <- OTvelocity[ii]
          if ( ii < (OTends[i] - (end_point_slide*rlength)) ) SEARCHINGforII <- FALSE
        }
        if ( !is.na(Evel) && (Evel > (min_Evel*1.02)) ) {
          k <- ii
          SEARCHINGforK <- TRUE
          while(SEARCHINGforK) {
            k <- k-1
            if ( !is.na(OTvelocity[k]) && (OTvelocity[k] <= (min_Evel*1.02)) ) SEARCHINGforK <- FALSE
            if ( k < (OTends[i] - (end_point_slide*rlength)) ) SEARCHINGforK <- FALSE
          }
        } else {
          k <- OTends[i]
        }
      } 
      k_vec <- c( k_vec, OTends[i] - k )
    }
  }
  k_vec <- k_vec[2:length(k_vec)]
  sync_tune <- as.integer(median(k_vec))
  # There might be NA in OTstarts and OTends, from removed shorts, but that's okay,
  #     since NA - sync_tune is still NA. 
  OTstarts[(last_no_feedback_trial_num+1):(pb_start-1)] <- OTstarts[(last_no_feedback_trial_num+1):(pb_start-1)] - sync_tune
  OTends[(last_no_feedback_trial_num+1):(pb_start-1)] <- OTends[(last_no_feedback_trial_num+1):(pb_start-1)] - sync_tune
  if ( is.na(sync_tune) ) {
    cat("\n WARNING! end sync_tune for Feedback is NA, participant:", ppp,"\n")
  }
  if (noisy_print) cat(", ES-F: ", sync_tune) 
  psv <- c( psv, sync_tune )
  # Now the reaches w/ model playback (post-feedback block)
  k_vec <- 1
  for ( i in pb_start:SUnumofreaches ) {
    if ( !(i %in% short_reaches) ) {
      rlength <- OTends[i] - OTstarts[i]
      end_velocities <- OTvelocity[(OTends[i]-(end_point_slide*rlength)):OTends[i]]
      k <- OTends[i]
      if ( any(!is.na(end_velocities)) ) {
        min_Evel <- min(OTvelocity[(OTends[i]-(end_point_slide*rlength)):OTends[i]],na.rm = TRUE)
        Evel <- OTvelocity[OTends[i]]
        ii <- OTends[i]
        SEARCHINGforII <- TRUE
        while ( is.na(Evel) && SEARCHINGforII ) {
          ii <- ii-1
          Evel <- OTvelocity[ii]
          if ( ii < (OTends[i] - (end_point_slide*rlength)) ) SEARCHINGforII <- FALSE
        }
        if ( !is.na(Evel) && (Evel > (min_Evel*1.02)) ) {
          k <- ii
          SEARCHINGforK <- TRUE
          while(SEARCHINGforK) {
            k <- k-1
            if ( !is.na(OTvelocity[k]) && (OTvelocity[k] <= (min_Evel*1.02)) ) SEARCHINGforK <- FALSE
            if ( k < (OTends[i] - (end_point_slide*rlength)) ) SEARCHINGforK <- FALSE
          }
        } else {
          k <- OTends[i]
        }
      } 
      k_vec <- c( k_vec, OTends[i] - k )
    }
  }
  k_vec <- k_vec[2:length(k_vec)]
  sync_tune <- as.integer(median(k_vec))
  # There might be NA in OTstarts and OTends, from removed shorts, but that's okay,
  #     since NA - sync_tune is still NA. 
  OTstarts[pb_start:SUnumofreaches] <- OTstarts[pb_start:SUnumofreaches] - sync_tune
  OTends[pb_start:SUnumofreaches] <- OTends[pb_start:SUnumofreaches] - sync_tune
  if ( is.na(sync_tune) ) {
    cat("\n WARNING! end sync_tune for Feedback is NA, participant:", ppp,"\n")
  }
  if (noisy_print) cat(", ES-PB: ", sync_tune) 
  psv <- c( psv, sync_tune )
  
  # If there are problems, or you suspect the initial guess for OTsampleslide is way off, these plots useful
  #   plot(ty,pch=19)
  #   abline(v=OTstarts)
  
  # Build the OT reach col
  OTreachcol <- rep(NA,length(OTtimes))
  OTtrialcol <- rep(NA,length(OTtimes))
  for ( i in 1:SUnumofreaches ) {
    if ( !(i %in% short_reaches) ) {
      OTreachcol[OTstarts[i]:OTends[i]] <- i
    }
  }
  OTtrialcol <- OTreachcol
  
  # Check for floating reaches
  num_of_floats <- 0
  floating_reaches <- NA
  for ( i in 1:SUnumofreaches ) {
    if ( !(i %in% short_reaches) ) {
      rStartPoint <- c( tx[OTstarts[i]], ty[OTstarts[i]], tz[OTstarts[i]] )
      iii <- OTstarts[i]
      while ( any(is.na(rStartPoint)) ) {
        iii <- iii+1
        rStartPoint <- c( tx[iii], ty[iii], tz[iii] )
      }
      if ( !is.numeric(vdist(IP,rStartPoint)) || is.na(vdist(IP,rStartPoint)) || length(vdist(IP,rStartPoint)) != 1 ) {
        cat("\n WARNING! float check vdist is not numeric, is NA, or isn't length 1, participant:", ppp, "reach:", i,"\n")
      }
      if ( vdist(IP,rStartPoint) > float_cutoff ) {
        SUreaches[SUstarts[i]:SUends[i]] <- NA
        SUstarts[i] <- NA
        SUends[i] <- NA
        SUlengths[i] <- NA
        OTreachcol[OTstarts[i]:OTends[i]] <- NA
        OTstarts[i] <- NA
        OTends[i] <- NA
        num_of_floats <- num_of_floats + 1
        floating_reaches <- c( floating_reaches, i ) 
      }
    }
  }
  floating_reaches <- floating_reaches[2:length(floating_reaches)]
  if (noisy_print) cat(", Fs: ", num_of_floats) 
  psv <- c( psv, num_of_floats )
  
  #if (noisy_print) plot(ty,type="p",pch=19)
  #for ( i in 1:SUnumofreaches ) {
  #  if (noisy_print) abline(v=OTends[i])
  #  if (noisy_print) abline(v=OTstarts[i])
  #}
  
  # Count how many reaches we have left
  OTreachcolNONA <- OTreachcol[!is.na(OTreachcol)]
  SUreachesNONA <- SUreaches[!is.na(SUreaches)]
  OTreachcolNONA <- unique(OTreachcolNONA)
  SUreachesNONA <- unique(SUreachesNONA)
  if ( length(OTreachcolNONA) != length(SUreachesNONA) ) {
    cat("\n WARNING! OT and SU reach columns don't contain the same reaches, participant:",ppp,"\n")
  }
  if (noisy_print) cat(", R (final): ", length(OTreachcolNONA),"\n") 
  psv <- c( psv, length(OTreachcolNONA) )
  psv <- c( psv, length(SUreachesNONA[which(SUreachesNONA<=last_no_feedback_trial_num)]) ) 
  psv <- c( psv, length(SUreachesNONA[which(SUreachesNONA>last_no_feedback_trial_num & SUreachesNONA<pb_start)]) )
  psv <- c( psv, length(SUreachesNONA[which(SUreachesNONA>=pb_start)]) )
  
  # Now find the average length of reaches, without and with feedback
  tally <- NA
  for ( i in SUreachesNONA ) {
    if ( i <= last_no_feedback_trial_num ) {
      tally <- c( tally, length(SUreaches[which(SUreaches==i)]) )
    }
  }
  tally <- tally[2:length(tally)]
  psv <- c( psv, mean(tally) )
  tally <- NA
  for ( i in SUreachesNONA ) {
    if ( i > last_no_feedback_trial_num && i < pb_start ) {
      tally <- c( tally, length(SUreaches[which(SUreaches==i)]) )
    }
  }
  tally <- tally[2:length(tally)]
  psv <- c( psv, mean(tally) )
  tally <- NA
  for ( i in SUreachesNONA ) {
    if ( i >= pb_start ) {
      tally <- c( tally, length(SUreaches[which(SUreaches==i)]) )
    }
  }
  tally <- tally[2:length(tally)]
  psv <- c( psv, mean(tally) )
  
  # Save reach columns
  motionOT$reach[which(motionOT$participant==ppp)] <<- OTreachcol
  motionOT$trials[which(motionOT$participant==ppp)] <<- OTtrialcol
  motionSU$reach[which(motionSU$participant==ppp)] <<- SUreaches
  
  # clean velocity columns: 
  OTindex <- which(motionOT$participant==ppp & 
                     is.na(motionOT$reach))
  motionOT$SpatialVelocity[OTindex] <<- NA
  SUindex <- which(motionSU$participant==ppp & 
                     is.na(motionSU$reach))
  motionSU$RotVelocity[SUindex] <<- NA
  
  part_chop_summary[index_,] <<- psv
  
}

# Begin by checking sonification unit data for corrupted sample lines 
SU_data_check <- function() {
  
  # Count here will include rows already found by the check in the setup script. 
  # This count should be (reasonably) accurate. 
  
  cat("\nChecking Sonification System Data: \n") 
  
  if ( length(SU_corrupted_rows) > 0 ) {
    for ( i in 1:length(SU_corrupted_rows) ) {
      SU_corrupted_rows[i] <- SU_corrupted_rows[i] - ( i - 1 ) 
    }
  }
  
  reach_col <- motionSU$reach
  sample_col <- motionSU$sample
  part_col <- motionSU$participant
  new_reach_row <- c(1) 
  corrupted_rows <- c()
  count <- 0 
  
  percentages <- seq(length(reach_col) * 0.1, length(reach_col) * 1, by = length(reach_col) * 0.1)
  percentages <- as.integer(percentages)
  
  cat("\nFinding reach start rows: ") 
  for ( i in 2:length(reach_col) ) {
    if ( is.na(i) ) {
      cat("\nWARNING! Unexppected i = NA\n")
    }
    if ( i %in% percentages ) {
      cat("*")
    }
    if ( reach_col[i] != reach_col[i-1] ) {
      if ( reach_col[i] == (reach_col[i-1]+1) ||  part_col[i] != part_col[i-1] ) {
        new_reach_row <- c( new_reach_row, i )
      } else {
        corrupted_rows <- c( corrupted_rows, i )
      }
    }
  }
  
  if ( any(is.na(new_reach_row)) ) {
    cat("\nWARNING! NA found in new_reach_row at index: ", which(is.na(new_reach_row)),"\n")
  }
  
  percentages <- seq(length(new_reach_row) * 0.1, length(new_reach_row) * 1, by = length(new_reach_row) * 0.1)
  percentages <- as.integer(percentages)
  
  cat("\nSearching reaches for bad rows: ") 
  for ( i in 1:length(new_reach_row) ) {
    if ( i %in% percentages ) {
      cat("*")
      if ( i == length(new_reach_row) ) {
        cat("\n")
      }
    }
    rows <- 1
    if ( i == 1 ) {
      rows <- 1:(new_reach_row[2]-1)
    } else if ( i == length(new_reach_row) ) {
      rows <- new_reach_row[i]:length(reach_col)
    } else {
      rows <- new_reach_row[i]:(new_reach_row[i+1]-1)
    }
    if ( length(rows) < 2 ) {
      cat("\nWARNING! Length of rows < 2")
    }
    for ( j in rows[2:length(rows)] ) {
      if ( is.na(j) ) {
        cat("\nWARNING! Unexppected j = NA\n")
      }

      if ( sample_col[j] != (sample_col[j-1]+1) ) {
        
        if ( length(SU_corrupted_rows) > 0 && !(j %in% SU_corrupted_rows) ) {
          corrupted_rows <- c( corrupted_rows, j )
        } else if ( length(SU_corrupted_rows) == 0 ) {
          corrupted_rows <- c( corrupted_rows, j )
        }
        
      } 
      
    }
    count <- count + ((sample_col[rows[length(rows)]]+1)-length(rows)) # sample count from SU starts at zero
  }
  
  #cat("\nCorrupted Rows:",corrupted_rows,"\n")
  if ( length(corrupted_rows) > 0 ) motionSU <<- motionSU[-corrupted_rows,]
  
  # Corrupted rows are literal corrupted rows in the SU dataframe.
  # Lost rows are rows of data that should have been recorded, but weren't.
  # Can be multiple lost rows per corrupted row. 
  
  num_of_corrupted_rows_motionSU <<- length(corrupted_rows)
  percentage_of_corrupted_rows_motionSU <<- (length(corrupted_rows)/length(reach_col)) * 100
  
  cat("\nNumber of corrupted rows: ", length(corrupted_rows), "\n")
  cat("Ratio of corrupted rows: ", length(corrupted_rows), "/", length(reach_col), "\n")
  cat("Percentage of corrupted rows: ", (length(corrupted_rows)/length(reach_col)) * 100, "%\n")
  cat("\n")
  cat("Corrupted Rows:", corrupted_rows, "\n") 
  
  num_of_lost_rows_motionSU <<- count
  percentage_of_lost_rows_motionSU <<- (count/length(reach_col)) * 100
  
  cat("\nNumber of lost rows: ", count, "\n")
  cat("Ratio of lost rows: ", count, "/", length(reach_col), "\n")
  cat("Percentage of lost rows: ", (count/length(reach_col)) * 100, "%\n")
  
}
SU_data_check() 

# Run this to chop all OptiTrack data and check all reaches for all participants as a batch
chop_data <- function() {
  
  # Raw OT data is in meters; put in cm
  motionOT$x <<- motionOT$x * 100
  motionOT$y <<- motionOT$y * 100
  motionOT$z <<- motionOT$z * 100
  
  # Copy SU reach column into a trial column; won't be cleaned out with NAs once bad reaches removed
  motionSU$trials <<- motionSU$reach
  motionOT$trials <<- rep(NA,length(motionOT$reach))
  
  cat("\nChopping and checking reaches\n")
  
  participant_vector <- pTable$number
  
  cat("\n")
  
  cat("Number of participants: ", length(participant_vector), "\n")
  cat("Participants: ", participant_vector, "\n") 
  
  cat("\n")
  
  part_chop_summary <- matrix( NA, nrow=length(participant_vector), ncol=20 )
  colnames(part_chop_summary) <- c( "Participant Number", 
                                    "RR", 
                                    "BSUS", 
                                    "ORs", 
                                    "slide",
                                    "Ss", 
                                    "S-R", 
                                    "S-F", 
                                    "S-PB",
                                    "ES-R", 
                                    "ES-F",
                                    "ES-PB",
                                    "Fs", 
                                    "R(f)", 
                                    "R(r)", 
                                    "R(fb)", 
                                    "R(L25)", 
                                    "MRLr(ms)", 
                                    "MRLfb(ms)",
                                    "MRLpb(ms)")
  part_chop_summary <<- data.frame(part_chop_summary)
  
  # This is the part of the function that chops the data.
  # Runs ParticipantChopOT for all participants.      
  l_ply(1:length(participant_vector),ParticipantChopOT)

  if ( model_check(FALSE) ) cat("\n WARNING! At least one model lost in chopping.\n")
  
  cat("\npart_chop_summary column names:\n
        RR = Raw Reaches, slide = OTsampleslide, BSUS = Bad SU Samples, ORs = Overruns, Ss = Shorts, 
        S-R = Random Sync, S-F = Feedback Sync, S-SB = Post-feedback block Sync, ES-R = Random Sync (end), 
        ES-F = Feedback Sync (end), ES-SB = Post-feedback block Sync (end), Fs = Floats, R(f) = Final Reaches, 
        R(r) = Random Reaches, R(fb) = Reaches w/ Feedback (#26-75), R(L25) = Final 25 Reache (post-feedback block)s,
        MRLr(ms) = Mean Reach Length, Random (ms), MRLfb(ms) = Mean Reach Length, Feedback (ms),
        MRLpb(ms) = Mean Reach Length, Post-feedback Block (ms)\n","\n" )
  
} 
chop_data()

# Remove bad rows in OT data found by hand on previous runs
if ( remove_bad_OT_rows ) {
  
  # Note: These bad points were found by hand on previous runs of processing the data, using 3D plots.
  
  cat("\nRemoving bad OT rows found on previous runs of data.\n")
  bad_60 <- which(motionOT$participant==60 & motionOT$x< -10 & motionOT$y<20 & !is.na(motionOT$trials))
  bad_64 <- which(motionOT$participant==64 & motionOT$x< -20 & !is.na(motionOT$trials))
  bad_58 <- which(motionOT$participant==58 & motionOT$x< -20 & !is.na(motionOT$trials))
  bad_56 <- which(motionOT$participant==56 & motionOT$x< -18 & !is.na(motionOT$trials))
  bad_rows <- c( bad_60, bad_64, bad_58, bad_56 )
  motionOT$x[bad_rows] <- NA
  motionOT$y[bad_rows] <- NA
  motionOT$z[bad_rows] <- NA
  
  num_OT_bad_rows_removed <<- length(bad_rows)
  percentage_OT_bad_rows_removed <<- (length(bad_rows)/length(which(!is.na(motionOT$trials))))*100
  
  cat("Number of rows removed:", length(bad_rows),"\n")
  cat("Total rows of OT data containing reaches:", length(which(!is.na(motionOT$trials))),"\n")
  cat("Rows removed as percentage of total rows containing reaches:", (length(bad_rows)/length(which(!is.na(motionOT$trials))))*100)
  
  removed_bad_rows <- TRUE # To ensure setup is always run before preprocessing.
  
}
