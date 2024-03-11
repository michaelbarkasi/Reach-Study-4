
# Initialize a dataframe to save error information for printing later
part_error_summary <- data.frame(
  Num = pTable$number,
  mpp.OT = rep(NA, length(pTable$number)),
  mpp.SU = rep(NA, length(pTable$number)),
  mrp.OT = rep(NA, length(pTable$number)),
  rep.OT = rep(NA, length(pTable$number)),
  rep.SUa = rep(NA, length(pTable$number)),
  rep.SUb = rep(NA, length(pTable$number))
)

# Major dependency functions
cat("Loading helper functions.\n")
cat("\n")
find_spatial_path_error <- function(printplots,ppp,useDTW) {
  
  # For all reaches R of participant ppp, this function finds the distance between each point of R
  #   and ppp's model, defined as the distance between that point and a matching point of R, where
  #   "matching" is done either by simple rescaling, or by DTW of velocity curves. At the end, finds
  #   the mean of these distances for each reach, called "path error". Note that, for comparing OT and SU data,
  #   this function also normalizes these values by dividing them by the length of ppp's model. 
  #   This normalization was used in previous studies, but not in RS$ (online vs terminal feedback). 
  
  # normalized error for the OT/SU correlation analysis (RS3) doesn't use DTW.
  
  model <- pTable$model[which(pTable$number==ppp)]
  MRD <- rep(NA,max(motionOT$reach,na.rm=TRUE)) # mean reach distance
  
  # Grab data
  tx <- motionOT$x[which(motionOT$participant==ppp)]
  ty <- motionOT$y[which(motionOT$participant==ppp)]
  tz <- motionOT$z[which(motionOT$participant==ppp)]
  reach <- motionOT$reach[which(motionOT$participant==ppp)]
  mx <- tx[which(reach==model)] 
  my <- ty[which(reach==model)] 
  mz <- tz[which(reach==model)] 
  
  avg_distances <- rep(NA,max(reach,na.rm = TRUE))
  model_length <- length(reach[which(reach==model)])
  checkifnumber(model_length)
  
  # Find length of the model for normalization
  model_path_length <- 0 
  for ( k in 2:model_length ) {
    model_path_length <- c( model_path_length, 
                            vdist(
                              c( mx[k],my[k],mz[k] ),
                              c( mx[k-1],my[k-1],mz[k-1] )
                            )
    )
  }
  model_path_length <- sum(model_path_length,na.rm=TRUE)
  checkifnumber(model_path_length)
  
  # Now find distances
  missed_points <- 0 
  found_points <- 0
  nr <- max(reach,na.rm = TRUE) # Number of reaches
  for ( i in 1:nr ) {
    if ( i %in% reach ) {
      r_start <- min(which(motionOT$reach==i & motionOT$participant==ppp),na.rm = TRUE)
      reach_length <- length(reach[which(reach==i)])
      checkifnumber(reach_length)
      distances <- rep(NA,reach_length)
      rtx <- tx[which(reach==i)]
      rty <- ty[which(reach==i)]
      rtz <- tz[which(reach==i)]
      for ( j in 1:reach_length ) {
        if ( j < 2 ) {
          model_point <- j
        } else {
          model_point <- as.integer( model_length * ( j / reach_length ) )
          checkifnumber(model_point)
        }
        if ( model_point == 0 ) model_point <- 1
        p <- c(rtx[j],rty[j],rtz[j])
        m <- c(mx[model_point],my[model_point],mz[model_point])
        if ( !any(is.na(p)) && !any(is.na(m)) ) {
          checkifnumericvector(p,3)
          checkifnumericvector(m,3)
          distances[j] <- round(vdist(p,m),significant_digitsOT)
          checkifnumber(distances[j])
          found_points <- found_points+1
        } else {
          missed_points <- missed_points+1
        }
      }
      avg_distances[i] <- round(mean(distances,na.rm=TRUE),significant_digitsOT)
      normalized_distances <- distances / model_path_length
      motionOT$Rerror[r_start:((r_start+reach_length)-1)] <<- distances
      motionOT$Nerror[r_start:((r_start+reach_length)-1)] <<- normalized_distances
      if (useDTW) {
        d <- DTW_model_indexesOT_Vel(ppp,i)
        mi <- d$index1
        ri <- d$index2
        distances <- rep(NA,length(ri))
        if ( length(mi) != length(ri) ) {
          cat("Something went wrong w/ DTW")
          return()
        }
        for ( j in 1:length(ri) ) {
          pj <- ri[j]
          mj <- mi[j]
          p <- c(rtx[pj],rty[pj],rtz[pj])
          m <- c(mx[mj],my[mj],mz[mj])
          if ( !any(is.na(p)) && !any(is.na(m)) ) {
            checkifnumericvector(p,3)
            checkifnumericvector(m,3)
            distances[j] <- round(vdist(p,m),significant_digitsOT)
            checkifnumber(distances[j])
          } 
        }
        avg_distances[i] <- round(mean(distances,na.rm=TRUE),significant_digitsOT)
      }
      checkifnumber(avg_distances[i])
      MRD[i] <- avg_distances[i]
    } 
  }
  
  # Print plots, if desired
  if ( printplots ) {
    #for coloring plots
    cv <- 1
    for ( i in 1:(model-1) ) {
      cv <- c( cv, R_color )
    }
    cv <- c( cv, M_color ) 
    for ( i in (model+1):last_no_feedback_trial_num ) {
      cv <- c( cv, R_color )
    }
    if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition1 ) {
      for ( i in (last_no_feedback_trial_num+1):(pb_start-1) ) {
        cv <- c( cv, C1_color ) 
      }
      for ( i in pb_start:max(reach,na.rm = TRUE) ) {
        cv <- c( cv, C1_color_dark ) 
      }
    } else if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition2 ) {
      for ( i in (last_no_feedback_trial_num+1):(pb_start-1) ) {
        cv <- c( cv, C2_color ) 
      }
      for ( i in pb_start:max(reach,na.rm = TRUE) ) {
        cv <- c( cv, C2_color_dark ) 
      }
    }
    cv <- cv[2:length(cv)]
    #plots
    plot(1:max(reach,na.rm = TRUE),avg_distances,pch=16,col=cv,
         main="Mean Distance from Model by Reach (XZY)", xlab="reach number",ylab="distance (meters)")
  }
  
  cat("Missed point percentage (OT):",(missed_points/(found_points+missed_points))*100,"\n")
  part_error_summary$mpp.OT[ppp] <<- (missed_points/(found_points+missed_points))*100
  
  return(MRD)
  
}
find_spatial_target_error <- function(printplots,ppp) {
  
  # For all reaches R of participant ppp, this function finds the distance between the last point of R
  #   and the last point of ppp's model, or between the nearest valid points to these points. 
  
  model <- pTable$model[which(pTable$number==ppp)]
  MRD <- rep(NA,max(motionOT$reach,na.rm=TRUE))
  
  # load data
  tx <- motionOT$x[which(motionOT$participant==ppp)]
  ty <- motionOT$y[which(motionOT$participant==ppp)]
  tz <- motionOT$z[which(motionOT$participant==ppp)]
  reach <- motionOT$reach[which(motionOT$participant==ppp)]
  mx <- tx[which(reach==model)] 
  my <- ty[which(reach==model)] 
  mz <- tz[which(reach==model)] 
  
  # Find valid model point nearest to its end
  model_end <- c(mx[length(mx)],my[length(my)],mz[length(mz)])
  if ( any(is.na(model_end)) || length(model_end) != 3 ) {
    NEEDMODELEND <- TRUE
    check <- length(mx)
    while (NEEDMODELEND) {
      check <- check-1
      model_end <- c(mx[check],my[check],mz[check])
      if ( !any(is.na(model_end)) && length(model_end) == 3 ) {
        NEEDMODELEND <- FALSE
      }
    }
    if ( check < length(mx) * 0.5 ) {
      cat("WARNING! model end (OT) found more 50% from last model sample\n")
    }
    cat("Model end (OT) found at %:",(check/length(mx))*100,"\n")
  }
  checkifnumericvector(model_end,3)
  
  # Find target distances
  missed_reach <- 0 
  found_reach <- 0
  missed_reach_list <- NA
  distances <- rep(NA,max(reach,na.rm = TRUE))
  check_avg <- NA
  for ( i in 1:max(reach,na.rm = TRUE) ) {
    if ( i %in% reach ) {
      rx <- tx[which(reach==i)] 
      ry <- ty[which(reach==i)] 
      rz <- tz[which(reach==i)] 
      r_end <- c(rx[length(rx)],ry[length(ry)],rz[length(rz)])
      SEARCHINGTOOLONG <- FALSE
      if ( any(is.na(r_end)) || length(r_end) != 3 ) {
        NEEDREND <- TRUE
        check <- length(rx)
        checkifnumber(check)
        while (NEEDREND) {
          check <- check-1
          r_end <- c(rx[check],ry[check],rz[check])
          if ( !any(is.na(r_end)) && length(r_end) == 3 ) {
            NEEDREND <- FALSE
          }
          if ( check < length(rx) * 0.35 ) SEARCHINGTOOLONG <- TRUE
          if ( SEARCHINGTOOLONG ) NEEDEND <- FALSE
        }
        check_avg <- c( check_avg, (check/length(rx))*100 )
      }
      check_avg <- c( check_avg, 100 )
      if ( SEARCHINGTOOLONG ) {
        missed_reach_list <- c( missed_reach_list, i ) 
        missed_reach <- missed_reach+1
      } else {
        checkifnumericvector(r_end,3)
        distances[i] <- round(vdist(model_end,r_end),significant_digitsOT)
        checkifnumber(distances[i])
        MRD[i] <- distances[i]
        found_reach <- found_reach+1
      }
    } 
  }
  if ( any(!is.na(check_avg)) ) {
    check_avg <- mean(check_avg,na.rm=TRUE)
    cat("Reach ends found on avg at %:",check_avg,"\n")
  }
  part_error_summary$rep.OT[ppp] <<- check_avg
  
  # Print plot of results, if desired
  if ( printplots ) {
    #for coloring plots
    cv <- 1
    for ( i in 1:(model-1) ) {
      cv <- c( cv, R_color )
    }
    cv <- c( cv, M_color ) 
    for ( i in (model+1):last_no_feedback_trial_num ) {
      cv <- c( cv, R_color )
    }
    if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition1 ) {
      for ( i in (last_no_feedback_trial_num+1):(pb_start-1) ) {
        cv <- c( cv, C1_color ) 
      }
      for ( i in pb_start:max(reach,na.rm = TRUE) ) {
        cv <- c( cv, C1_color_dark ) 
      }
    } else if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition2 ) {
      for ( i in (last_no_feedback_trial_num+1):(pb_start-1) ) {
        cv <- c( cv, C2_color ) 
      }
      for ( i in pb_start:max(reach,na.rm = TRUE) ) {
        cv <- c( cv, C2_color_dark ) 
      }
    }
    cv <- cv[2:length(cv)]
    plot(1:max(reach,na.rm = TRUE),distances,pch=16,col=cv,
         main="Distance to Model End by Reach (XZY space)", xlab="reach number",ylab="distance (meters)")
  }
  
  if ( length(missed_reach_list) > 1 ) {
    missed_reach_list <- missed_reach_list[2:length(missed_reach_list)]
    cat("Missed reach percentage (OT):",(missed_reach/(found_reach+missed_reach))*100,"Reaches:",missed_reach_list,"\n")
    part_error_summary$mrp.OT[ppp] <<- (missed_reach/(found_reach+missed_reach))*100
  } else {
    cat("Missed reach percentage (OT):",(missed_reach/(found_reach+missed_reach))*100,"\n")
    part_error_summary$mrp.OT[ppp] <<- (missed_reach/(found_reach+missed_reach))*100
  }
  
  return(MRD)
  
}
find_rotation_path_error <- function(printplots,ppp,useDTW) {
  
  # For all reaches R of participant ppp, this function finds the distance between each point of R
  #   and ppp's model, defined as the distance between that point and a matching point of R, where
  #   "matching" is done either by simple rescaling, or by DTW of velocity curves. At the end, finds
  #   the mean of these distances for each reach, called "path error". Note that, for comparing OT and SU data,
  #   this function also normalizes these values by dividing them by the length of ppp's model. 
  #   This normalization was used in previous studies, but not in RS$ (online vs terminal feedback). 
  
  # normalized error for the OT/SU correlation analysis (RS3) doesn't use DTW.
  
  model <- pTable$model[which(pTable$number==ppp)]
  MRD <- rep(NA,max(motionOT$reach,na.rm=TRUE))
  
  # Grab data
  atx <- motionSU$qax[which(motionSU$participant==ppp)]
  aty <- motionSU$qay[which(motionSU$participant==ppp)]
  atz <- motionSU$qaz[which(motionSU$participant==ppp)]
  atr <- motionSU$qar[which(motionSU$participant==ppp)]
  btx <- motionSU$qbx[which(motionSU$participant==ppp)]
  bty <- motionSU$qby[which(motionSU$participant==ppp)]
  btz <- motionSU$qbz[which(motionSU$participant==ppp)]
  btr <- motionSU$qbr[which(motionSU$participant==ppp)]
  reach <- motionSU$reach[which(motionSU$participant==ppp)]
  amx <- atx[which(reach==model)] 
  amy <- aty[which(reach==model)] 
  amz <- atz[which(reach==model)] 
  amr <- atr[which(reach==model)] 
  bmx <- btx[which(reach==model)] 
  bmy <- bty[which(reach==model)] 
  bmz <- btz[which(reach==model)] 
  bmr <- btr[which(reach==model)] 
  
  avg_distances <- rep(NA,max(reach,na.rm = TRUE))
  model_length <- length(reach[which(reach==model)])
  checkifnumber(model_length)
  
  # Find length of the model for normalization
  model_path_length <- 0 
  for ( k in 2:model_length ) {
    model_path_length <- c( model_path_length, 
                            vdist(
                              c( amx[k],amy[k],amz[k],amr[k] ),
                              c( amx[k-1],amy[k-1],amz[k-1],amr[k-1] )
                            ) + 
                              vdist(
                                c( bmx[k],bmy[k],bmz[k],bmr[k] ),
                                c( bmx[k-1],bmy[k-1],bmz[k-1],bmr[k-1] )
                              )
    )
  }
  model_path_length <- sum(model_path_length,na.rm=TRUE)
  checkifnumber(model_path_length)
  
  # Now find distances
  missed_points <- 0 
  found_points <- 
    nr <- max(reach,na.rm = TRUE) # number of reaches
  for ( i in 1:nr ) {
    if ( i %in% reach ) {
      r_start <- min(which(motionSU$reach==i & motionSU$participant==ppp),na.rm = TRUE)
      reach_length <- length(reach[which(reach==i)])
      checkifnumber(reach_length)
      distances <- rep(NA,reach_length)
      artx <- atx[which(reach==i)]
      arty <- aty[which(reach==i)]
      artz <- atz[which(reach==i)]
      artr <- atr[which(reach==i)]
      brtx <- btx[which(reach==i)]
      brty <- bty[which(reach==i)]
      brtz <- btz[which(reach==i)]
      brtr <- btr[which(reach==i)]
      
      for ( j in 1:reach_length ) {
        if ( j < 2 ) {
          model_point <- j
        } else {
          model_point <- as.integer( model_length * ( j / reach_length ) )
          checkifnumber(model_point)
        }
        if ( model_point == 0 ) model_point <- 1
        p1 <- c(artx[j],arty[j],artz[j],artr[j])
        p2 <- c(brtx[j],brty[j],brtz[j],brtr[j])
        m1 <- c(amx[model_point],amy[model_point],amz[model_point],amr[model_point])
        m2 <- c(bmx[model_point],bmy[model_point],bmz[model_point],bmr[model_point])
        if ( !any(is.na(p1)) && !any(is.na(p2)) && !any(is.na(m1)) && !any(is.na(m2)) ) {
          checkifnumericvector(p1,4)
          checkifnumericvector(p2,4)
          checkifnumericvector(m1,4)
          checkifnumericvector(m2,4)
          distances[j] <- round(vdist(p1,m1) + vdist(p2,m2),significant_digitsSU)
          checkifnumber(distances[j])
          found_points <- found_points+1
        } else {
          missed_points <- missed_points+1
        }
      }
      avg_distances[i] <- round(mean(distances,na.rm=TRUE),significant_digitsSU)
      normalized_distances <- distances / model_path_length
      motionSU$Rerror[r_start:((r_start+reach_length)-1)] <<- distances
      motionSU$Nerror[r_start:((r_start+reach_length)-1)] <<- normalized_distances
      if (useDTW) {
        d <- DTW_model_indexesSU_Vel(ppp,i)
        mi <- d$index1
        ri <- d$index2
        distances <- rep(NA,length(ri))
        if ( length(mi) != length(ri) ) {
          cat("Something went wrong w/ DTW")
          return()
        }
        for ( j in 1:length(ri) ) {
          pj <- ri[j]
          mj <- mi[j]
          p1 <- c(artx[pj],arty[pj],artz[pj],artr[pj])
          p2 <- c(brtx[pj],brty[pj],brtz[pj],brtr[pj])
          m1 <- c(amx[mj],amy[mj],amz[mj],amr[mj])
          m2 <- c(bmx[mj],bmy[mj],bmz[mj],bmr[mj])
          if ( !any(is.na(p1)) && !any(is.na(p2)) && !any(is.na(m1)) && !any(is.na(m2)) ) {
            checkifnumericvector(p1,4)
            checkifnumericvector(p2,4)
            checkifnumericvector(m1,4)
            checkifnumericvector(m2,4)
            distances[j] <- round(vdist(p1,m1) + vdist(p2,m2),significant_digitsSU)
            checkifnumber(distances[j])
          } 
        }
        avg_distances[i] <- round(mean(distances,na.rm=TRUE),significant_digitsSU)
      }
      checkifnumber(avg_distances[i])
      MRD[i] <- avg_distances[i]
      
    } 
    
  }
  
  # Print plots, if desired
  if ( printplots ) {
    #for coloring plots
    cv <- 1
    for ( i in 1:(model-1) ) {
      cv <- c( cv, R_color )
    }
    cv <- c( cv, M_color ) 
    for ( i in (model+1):last_no_feedback_trial_num ) {
      cv <- c( cv, R_color )
    }
    if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition1 ) {
      for ( i in (last_no_feedback_trial_num+1):(pb_start-1) ) {
        cv <- c( cv, C1_color ) 
      }
      for ( i in pb_start:max(reach,na.rm = TRUE) ) {
        cv <- c( cv, C1_color_dark ) 
      }
    } else if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition2 ) {
      for ( i in (last_no_feedback_trial_num+1):(pb_start-1) ) {
        cv <- c( cv, C2_color ) 
      }
      for ( i in pb_start:max(reach,na.rm = TRUE) ) {
        cv <- c( cv, C2_color_dark ) 
      }
    }
    cv <- cv[2:length(cv)]
    #plots
    plot(1:max(reach,na.rm = TRUE),avg_distances,pch=16,col=cv,
         main="Distance to Model End by Reach (QUAT space)", 
         xlab="reach number",ylab="distance (quaternion units)")
    
  }
  
  cat("Missed point percentage (SU):",(missed_points/(found_points+missed_points))*100,"\n")
  part_error_summary$mpp.SU[ppp] <<- (missed_points/(found_points+missed_points))*100
  
  return(MRD)
  
}
find_rotation_target_error <- function(printplots,ppp) {
  
  # For all reaches R of participant ppp, this function finds the distance between the last point of R
  #   and the last point of ppp's model, or between the nearest valid points to these points. 
  
  model <- pTable$model[which(pTable$number==ppp)]
  MRD <- rep(NA,max(motionOT$reach,na.rm=TRUE))
  
  # load data
  atx <- motionSU$qax[which(motionSU$participant==ppp)]
  aty <- motionSU$qay[which(motionSU$participant==ppp)]
  atz <- motionSU$qaz[which(motionSU$participant==ppp)]
  atr <- motionSU$qar[which(motionSU$participant==ppp)]
  btx <- motionSU$qbx[which(motionSU$participant==ppp)]
  bty <- motionSU$qby[which(motionSU$participant==ppp)]
  btz <- motionSU$qbz[which(motionSU$participant==ppp)]
  btr <- motionSU$qbr[which(motionSU$participant==ppp)]
  reach <- motionSU$reach[which(motionSU$participant==ppp)]
  amx <- atx[which(reach==model)] 
  amy <- aty[which(reach==model)] 
  amz <- atz[which(reach==model)] 
  amr <- atr[which(reach==model)] 
  bmx <- btx[which(reach==model)] 
  bmy <- bty[which(reach==model)] 
  bmz <- btz[which(reach==model)] 
  bmr <- btr[which(reach==model)] 
  
  # Find valid model point nearest to its end
  model_end1 <- c(amx[length(amx)],amy[length(amy)],amz[length(amz)],amr[length(amr)])
  model_end2 <- c(bmx[length(bmx)],bmy[length(bmy)],bmz[length(bmz)],bmr[length(bmr)])
  if ( any(is.na(model_end1)) || length(model_end1) != 4 ) {
    NEEDMODELEND <- TRUE
    check <- length(amx)
    while (NEEDMODELEND) {
      check <- check-1
      model_end1 <- c(amx[check],amy[check],amz[check],amr[check])
      if ( !any(is.na(model_end1)) && length(model_end1) == 4 ) {
        NEEDMODELEND <- FALSE
      }
    }
    if ( check < length(amx) * 0.5 ) {
      cat("WARNING! model end (SUa) found more 50% from last model sample\n")
    }
    cat("Model end (SUa) found at %:",(check/length(mx))*100,"\n")
  }
  if ( any(is.na(model_end2)) || length(model_end2) != 4 ) {
    NEEDMODELEND <- TRUE
    check <- length(bmx)
    while (NEEDMODELEND) {
      check <- check-1
      model_end2 <- c(bmx[check],bmy[check],bmz[check],bmr[check])
      if ( !any(is.na(model_end2)) && length(model_end2) == 4 ) {
        NEEDMODELEND <- FALSE
      }
    }
    if ( check < length(bmx) * 0.5 ) {
      cat("WARNING! model end (SUb) found more 50% from last model sample\n")
    }
    cat("Model end (SUb) found at %:",(check/length(mx))*100,"\n")
  }
  checkifnumericvector(model_end1,4)
  checkifnumericvector(model_end2,4)
  
  # Find target distances
  missed_reach <- 0 
  found_reach <- 0
  missed_reach_list <- NA
  distances <- rep(NA,max(reach,na.rm = TRUE))
  check_avgA <- NA
  check_avgB <- NA
  for ( i in 1:max(reach,na.rm = TRUE) ) {
    if ( i %in% reach ) {
      arx <- atx[which(reach==i)] 
      ary <- aty[which(reach==i)] 
      arz <- atz[which(reach==i)] 
      arr <- atr[which(reach==i)]
      brx <- btx[which(reach==i)] 
      bry <- bty[which(reach==i)] 
      brz <- btz[which(reach==i)] 
      brr <- btr[which(reach==i)]
      r_end1 <- c(arx[length(arx)],ary[length(ary)],arz[length(arz)],arr[length(arr)])
      r_end2 <- c(brx[length(brx)],bry[length(bry)],brz[length(brz)],brr[length(brr)])
      SEARCHINGTOOLONG <- FALSE
      if ( any(is.na(r_end1)) || length(r_end1) != 4 ) {
        NEEDREND <- TRUE
        check <- length(arx)
        checkifnumber(check)
        while (NEEDREND) {
          check <- check-1
          r_end1 <- c(arx[check],ary[check],arz[check],arr[check])
          if ( !any(is.na(r_end1)) && length(r_end1) == 4 ) {
            NEEDREND <- FALSE
          }
          if ( check < length(arx) * 0.35 ) SEARCHINGTOOLONG <- TRUE
          if ( SEARCHINGTOOLONG ) NEEDEND <- FALSE
        }
        check_avgA <- c( check_avgA, ((check/length(arx))*100) ) 
      }
      check_avgA <- c( check_avgA, 100 )
      if ( any(is.na(r_end2)) || length(r_end2) != 4 ) {
        NEEDREND <- TRUE
        check <- length(brx)
        checkifnumber(check)
        if ( is.na(length(brx)) || !is.numeric(length(brx)) || !(length(brx)>0) ) {
          cat("WARNING! length of brx isn't a number:",length(brx),"\n")
        }
        while (NEEDREND) {
          check <- check-1
          r_end2 <- c(brx[check],bry[check],brz[check],brr[check])
          if ( !any(is.na(r_end2)) && length(r_end2) == 4 ) {
            NEEDREND <- FALSE
          }
          if ( check < length(brx) * 0.35 ) SEARCHINGTOOLONG <- TRUE
          if ( SEARCHINGTOOLONG ) NEEDEND <- FALSE
        }
        check_avgB <- c( check_avgB, ((check/length(brx))*100) )
      }
      check_avgB <- c( check_avgB, 100 )
      if ( SEARCHINGTOOLONG ) {
        missed_reach_list <- c( missed_reach_list, i ) 
        missed_reach <- missed_reach+1
      } else {
        checkifnumericvector(r_end1,4)
        checkifnumericvector(r_end2,4)
        distances[i] <- round(vdist(model_end1,r_end1) + vdist(model_end2,r_end2),significant_digitsSU)
        checkifnumber(distances[i])
        MRD[i] <- distances[i]
        found_reach <- found_reach+1
      }
    } 
  }
  if ( any(!is.na(check_avgA)) ) {
    check_avgA <- mean(check_avgA,na.rm=TRUE)
    cat("Reach ends (SUa) found on avg at %:",check_avgA,"\n")
    part_error_summary$rep.SUa[ppp] <<- check_avgA
  }
  if ( any(!is.na(check_avgB)) ) {
    check_avgB <- mean(check_avgB,na.rm=TRUE)
    cat("Reach ends (SUb) found on avg at %:",check_avgB,"\n")
    part_error_summary$rep.SUb[ppp] <<- check_avgB
  }
  
  # Print plot of results, if desired
  if ( printplots ) {
    #for coloring plots
    cv <- 1
    for ( i in 1:(model-1) ) {
      cv <- c( cv, R_color )
    }
    cv <- c( cv, M_color ) 
    for ( i in (model+1):last_no_feedback_trial_num ) {
      cv <- c( cv, R_color )
    }
    if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition1 ) {
      for ( i in (last_no_feedback_trial_num+1):(pb_start-1) ) {
        cv <- c( cv, C1_color ) 
      }
      for ( i in pb_start:max(reach,na.rm = TRUE) ) {
        cv <- c( cv, C1_color_dark ) 
      }
    } else if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition2 ) {
      for ( i in (last_no_feedback_trial_num+1):(pb_start-1) ) {
        cv <- c( cv, C2_color ) 
      }
      for ( i in pb_start:max(reach,na.rm = TRUE) ) {
        cv <- c( cv, C2_color_dark ) 
      }
    }
    cv <- cv[2:length(cv)]
    #plots
    plot(1:max(reach,na.rm = TRUE),distances,pch=16,col=cv,
         main="Distance to Model End by Reach (QUAT space)", 
         xlab="reach number",ylab="distance (quaternion units)")
    
  }
  
  if ( length(missed_reach_list) > 1 ) {
    missed_reach_list <- missed_reach_list[2:length(missed_reach_list)]
    cat("Missed reach percentage (SU):",(missed_reach/(found_reach+missed_reach))*100,"Reaches:",missed_reach_list,"\n")
  } else {
    cat("Missed reach percentage (SU):",(missed_reach/(found_reach+missed_reach))*100,"\n")
  }
  
  return(MRD)
  
}
error_curve <- function( participant_, expected_length = 100 ) {
  
  # output: row 1 = SPE, row 2 = STE, row 3 = RPE, row 4 = RTE
  
  # grab raw error measurements (error for each of the 100 reaches)
  SPE <- rTable$SpatialPathError[which( rTable$Participants == participant_ )]
  STE <- rTable$SpatialTargetError[which( rTable$Participants == participant_ )]
  RPE <- rTable$RotationPathError[which( rTable$Participants == participant_ )]
  RTE <- rTable$RotationTargetError[which( rTable$Participants == participant_ )]
  
  # sanity check to make sure these are the right length
  if ( length(SPE) != expected_length || length(STE) != expected_length || length(RPE) != expected_length || length(RTE) != expected_length) {
    cat("\nWarning! error vectors not the expected length")
  }
  
  # Find mean error in random/no FB reaches (remove model)
  index__ <- 1:last_no_feedback_trial_num
  if ( remove_model ) {
    m_num <- pTable$model[which(pTable$number==participant_)]
    index__ <- index__[-m_num]
  }
  mean_nFB_SPE <- mean( SPE[index__], na.rm = TRUE )
  mean_nFB_STE <- mean( STE[index__], na.rm = TRUE )
  mean_nFB_RPE <- mean( RPE[index__], na.rm = TRUE )
  mean_nFB_RTE <- mean( RTE[index__], na.rm = TRUE )
  
  # Now normalize all reaches as a percentage change from mean jerk during random reaching
  SPE <- ( SPE - mean_nFB_SPE ) / mean_nFB_SPE
  STE <- ( STE - mean_nFB_STE ) / mean_nFB_STE
  RPE <- ( RPE - mean_nFB_RPE ) / mean_nFB_RPE
  RTE <- ( RTE - mean_nFB_RTE ) / mean_nFB_RTE
  
  output <- array(NA, dim = c(4,expected_length) )
  
  output[1,] <- SPE
  output[2,] <- STE
  output[3,] <- RPE
  output[4,] <- RTE
  
  return(output)
  
}
fill_rTable_NormError <- function() {
  
  for ( p in unique(rTable$Participants) ) {
    normalized_error <- error_curve(p)
    rTable$SpatialPathErrorNorm[which(rTable$Participants==p)] <<- normalized_error[1,]
    rTable$SpatialTargetErrorNorm[which(rTable$Participants==p)] <<- normalized_error[2,]
    rTable$RotationPathErrorNorm[which(rTable$Participants==p)] <<- normalized_error[3,]
    rTable$RotationTargetErrorNorm[which(rTable$Participants==p)] <<- normalized_error[4,]
  }
  
}

# Compute error and put into rTable
fill_rTable <-function(useDTW) {
  
  participant_vector <- pTable$number
  
  Participants <- c()
  ReachNum <- c()
  SpatialPathError <- c() # = Mean Reach Spatial Error, i.e mean "path" error in physical space
  SpatialTargetError <- c() # = Mean Final Spatial Error, i.e. mean "target" error in physical space
  RotationPathError <- c() # = Mean Reach Rotation Error, i.e. mean "path" error in rotation space
  RotationTargetError <- c() # = Mean Final Rotation Error, i.e. mean "target" error in rotation space
  
  mr <- max(motionOT$reach,na.rm=TRUE) # Max reach number
  
  # For computing and saving "real" error and normalized error, 
  #   at each point in a reach, defined as distance to its scaled matching model point. 
  motionOT$Rerror <<- rep(NA,length(motionOT$reach))
  motionOT$Nerror <<- rep(NA,length(motionOT$reach))
  motionSU$Rerror <<- rep(NA,length(motionSU$reach))
  motionSU$Nerror <<- rep(NA,length(motionSU$reach))
  
  # Compute spatial and rotation, path and target errors
  for ( i in participant_vector ) {
    
    temp_rv <- motionOT$reach[which(motionOT$participant==i)]
    
    Participants <- c( Participants, rep(i,mr) )
    
    for ( j in 1:mr ) {
      if ( j %in% temp_rv ) {
        ReachNum <- c( ReachNum, j )
      } else {
        ReachNum <- c( ReachNum, NA ) 
      }
    }
    
    cat("Computing reach errors for participant: ", i,"\n")
    SpatialPathError <- c( SpatialPathError, find_spatial_path_error(FALSE,i,useDTW) )
    SpatialTargetError <- c( SpatialTargetError, find_spatial_target_error(FALSE,i) )
    RotationPathError <- c( RotationPathError, find_rotation_path_error(FALSE,i,useDTW) )
    RotationTargetError <- c( RotationTargetError, find_rotation_target_error(FALSE,i) )
    
  }
  
  # Build rTable (contains data on reaches)
  Condition <- rep(NA,length(Participants))
  Condition2 <- rep(NA,length(Participants))
  Feedback <- rep(NA,length(Participants))
  for ( i in 1:length(Participants) ) {
    if ( !is.na(ReachNum[i]) && 
         ReachNum[i] <= last_no_feedback_trial_num && 
         pTable$model[which(pTable$number==Participants[i])] != ReachNum[i] 
    ) {
      Condition2[i] <- NF_condition
      Feedback[i] <- "No"
    } else if ( !is.na(ReachNum[i]) && ReachNum[i] > last_no_feedback_trial_num && ReachNum[i] < pb_start ) {
      Condition2[i] <- pTable$Condition[which(pTable$number==Participants[i])] 
      Feedback[i] <- "Yes"
    } else if ( !is.na(ReachNum[i]) && ReachNum[i] > pb_start ) {
      Condition2[i] <- pb_condition 
      Feedback[i] <- "PB"
    }
  }
  for ( i in 1:length(Participants) ) {
    Condition[i] <- pTable$Condition[which(pTable$number==Participants[i])] 
  }
  rTable <<- data.frame(Participants,
                        Condition,
                        Condition2,
                        Feedback,
                        ReachNum,
                        SpatialPathError,
                        SpatialTargetError,
                        RotationPathError,
                        RotationTargetError)
  rTable$trials <<- rep(NA,length(rTable$Participants))
  for ( p in pTable$number ) {
    rTable$trials[which(rTable$Participants==p)] <<- 1:100
  }
  
} 
fill_rTable(TRUE) 

# Create normalized error columns in rTable
cat("\nCreating normalized error columns in rTable.\n")
create_norm_jerk_rTable_columns <- function() {
  
  rTable$SpatialPathErrorNorm <<- rep(NA,length(rTable$Participants))
  rTable$SpatialTargetErrorNorm <<- rep(NA,length(rTable$Participants))
  rTable$RotationPathErrorNorm <<- rep(NA,length(rTable$Participants))
  rTable$RotationTargetErrorNorm <<- rep(NA,length(rTable$Participants))
  
}
create_norm_jerk_rTable_columns()

# Filling rTable with normalized error
cat("\nFilling rTable with normalized error values.\n")
fill_rTable_NormError()

# Need to run this to convert Condition and Condition2 columns to factor for plotting
convert_Condition <- function() {
  pTable$Condition <<- as.factor(pTable$Condition)
  rTable$Condition <<- as.factor(rTable$Condition)
  rTable$Condition2 <<- as.factor(rTable$Condition2)
  rTable$Feedback <<- as.factor(rTable$Feedback)
}
convert_Condition()

# keep only trial rows in motionOT
motionOT <<- motionOT[!is.na(motionOT$trials),]

#count reaches left after shorts and floats are removed: 
cat("\nReaches left after shorts and floats are removed:", length(which(!is.na(rTable$ReachNum))))
reach_after_SF_removed <- length(which(!is.na(rTable$ReachNum)))