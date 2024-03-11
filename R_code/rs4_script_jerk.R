# Jerk Analysis RS4
cat("\nJerk analysis for study data.\n")

# Dependency functions 
cat("\nLoading helper functions.\n")
first_derivative <- function(v) {
  output <- diff(v)
  return(c(NA, output))
}
Nyquist_normalized_from_Hz <- function(hz_,sample_rate_) {
  nn <- (hz_*2) / sample_rate_
  return(nn)
}
jerk_diff <- function() {
  running_diff <- mean(abs(rTable$SUjerkNorm - rTable$OTjerkNorm),na.rm=TRUE)*100
  return(running_diff) 
}
fill_rTable_jerk <- function(FO_OT, # filter order, OT
                             FC_OT, # filter cutoff, OT
                             FO_SU, # filter order, SU
                             FC_SU,  # filter cutoff, SU
                             print_=TRUE,
                             sample_num_=0
) {
  
  #cat("\nFilling rTable with mean directional jerk magnitude numbers\n")
  
  OTjerk <- rep(NA,length(rTable$ReachNum))
  SUjerk <- rep(NA,length(rTable$ReachNum))
  
  participant_vector = pTable$number
  
  if (print_) cat("\nParticipant: ")
  
  for ( p in participant_vector ) {
    
    if (print_) {
      if ( p == participant_vector[1]) {
        cat(p)
      } else {
        cat(",",p)
      }
    }
    
    reaches <- rTable$ReachNum[which(rTable$Participants==p)]
    
    reaches_before <- (which(pTable$number==p)-1)*100 # 100 is number of reaches done by each participant
    counter <- 0
    
    for ( i in 1:length(reaches) ) {
      
      counter <- counter+1
      r_row_num <- reaches_before + counter
      
      r <- reaches[i]
      
      if ( !is.na(r) ) {
        
        OTjerk[r_row_num] <- mean_reach_jerk_OT(p,r,PLOT=FALSE,FILTER=TRUE,forder=FO_OT,fcutoff=FC_OT)
        SUjerk[r_row_num] <- mean_reach_jerk_SU(p,r,PLOT=FALSE,FILTER=TRUE,forder=FO_SU,fcutoff=FC_SU)
        
      }
      
    }
    
  }
  
  rTable$OTjerk <<- OTjerk
  rTable$SUjerk <<- SUjerk
  
}
fill_rTable_NormJerk <- function(print_=TRUE) {
  
  if (print_) cat("\nParticipant: ")
  
  participant_vector = pTable$number
  
  for ( p in participant_vector ) {
    
    if (print_) {
      if ( p == participant_vector[1]) {
        cat(p)
      } else {
        cat(",",p)
      }
    }
    
    # rTable
    normalized_jerk <- jerk_curve(p)
    rTable$OTjerkNorm[which(rTable$Participants==p)] <<- normalized_jerk[1,]
    rTable$SUjerkNorm[which(rTable$Participants==p)] <<- normalized_jerk[2,]
    
  }
  
}
jerk_curve <- function( participant_, expected_length = 100 ) {
  
  # output: row 1 = OT, row 2 = SU
  
  # grab raw jerk measurements (mean directional jerk magnetude for each of the 100 reaches)
  OTjerk <- rTable$OTjerk[which( rTable$Participants == participant_ )]
  SUjerk <- rTable$SUjerk[which( rTable$Participants == participant_ )]
  
  # sanity check to make sure these are the right length
  if ( length(OTjerk) != expected_length || length(SUjerk) != expected_length ) {
    cat("\nWarning! jerk vectors not the expected length")
  }
  
  # Find mean jerk in random/no FB reaches 
  mean_nFBjerk_OT <- mean( OTjerk[1:last_no_feedback_trial_num], na.rm = TRUE )
  mean_nFBjerk_SU <- mean( SUjerk[1:last_no_feedback_trial_num], na.rm = TRUE )
  
  # Now normalize all reaches as a percentage change from mean jerk during random reaching
  OTjerk <- ( OTjerk - mean_nFBjerk_OT ) / mean_nFBjerk_OT
  SUjerk <- ( SUjerk - mean_nFBjerk_SU ) / mean_nFBjerk_SU
  
  output <- array(NA, dim = c(2,expected_length) )
  
  output[1,] <- OTjerk
  output[2,] <- SUjerk
  
  return(output)
  
}

# Create normalized jerk columns in rTable
cat("\nCreating normalized jerk columns in rTable.\n")
create_norm_jerk_rTable_columns <- function() {
  
  rTable$SUjerkNorm <<- rep(NA,length(rTable$Participants))
  rTable$OTjerkNorm <<- rep(NA,length(rTable$Participants))

}
create_norm_jerk_rTable_columns()

# Functions for finding mean directional jerk magnitude 
cat("\nLoading functions for finding mean directional jerk magnitude (for RS4).\n")
mean_reach_jerk_OT <- function(p,r,PLOT,FILTER=TRUE,forder,fcutoff,cutoff=as.integer(cutoffSU/10),return_output=TRUE,PLOTmag=FALSE) {
  
  index <- which( motionOT$participant==p & motionOT$reach==r )
  
  x <- motionOT$x[index]
  y <- motionOT$y[index]
  z <- motionOT$z[index]
  
  x_na <- which(is.na(x))
  y_na <- which(is.na(y))
  z_na <- which(is.na(z))
  
  na_index <- unique(c(x_na, y_na, z_na))
  na_index_low <- 0
  na_index_high <- 0
  
  # Filter
  if (FILTER) {

    if (length(na_index) >= 1) {
      if (any(na_index < cutoff)) na_index_low <- length(na_index < cutoff)
      if (any(na_index > (length(index)-cutoff))) na_index_high <- length(na_index > (length(index)-cutoff))
      x <- x[-na_index]
      y <- y[-na_index]
      z <- z[-na_index]
    }
    
    x <- filter(butter(forder, fcutoff, type = "low"), x)
    y <- filter(butter(forder, fcutoff, type = "low"), y)
    z <- filter(butter(forder, fcutoff, type = "low"), z)
    
  }

  # cutoff first and last 200ms of recording
  if ( length(x) != length(y) ||
       length(x) != length(z) ||
       length(z) != length(y) ) {
    cat("\n\nWarning! length of x, y, and z not all the same!\n\n")
  }
  cut_index <- ((1+cutoff)-na_index_low):((length(x)-cutoff)+na_index_high)
  
  x <- x[cut_index]
  y <- y[cut_index]
  z <- z[cut_index]
  
  vx <- first_derivative(x) * OT_sample_rate # cm/s
  vy <- first_derivative(y) * OT_sample_rate
  vz <- first_derivative(z) * OT_sample_rate
  
  ax <- first_derivative(vx) * OT_sample_rate # cm/s^2
  ay <- first_derivative(vy) * OT_sample_rate
  az <- first_derivative(vz) * OT_sample_rate
  
  jx <- first_derivative(ax) * OT_sample_rate # cm/s^3
  jy <- first_derivative(ay) * OT_sample_rate
  jz <- first_derivative(az) * OT_sample_rate
  
  if (PLOT) {
    fig_num <- fig_num + 1
    plot_title <- "Spatial dimension x, Unfiltered"
    if (FILTER) {
      plot_title <- "Spatial dimension x, Filtered"
    }
    if (print_all_R_output) plot_title <- paste("Fig ", sec_num, ".", fig_num, ": ", plot_title, sep = "")
    xplot <- ggplot(data = data.frame(cut_index = cut_index, x = x)) +
      geom_point(aes(x = cut_index, y = x), color = "blue", size = 1.5, shape = 16) +
      labs(x = "sample", y = expression(cm)) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      theme_minimal() +
      theme(text=element_text(size=7))
    vplot <- ggplot(data = data.frame(cut_index = cut_index, vx = vx)) +
      geom_point(aes(x = cut_index, y = vx), color = "gray", size = 1.5, shape = 16) +
      labs(x = "sample", y = expression(cm/s)) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      theme_minimal() +
      theme(text=element_text(size=7))
    aplot <- ggplot(data = data.frame(cut_index = cut_index, ax = ax)) +
      geom_point(aes(x = cut_index, y = ax), color = "black", size = 1.5, shape = 16) +
      labs(x = "sample", y = expression(cm/s^2)) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      theme_minimal() +
      theme(text=element_text(size=7))
    jplot <- ggplot(data = data.frame(cut_index = cut_index, jx = jx)) +
      geom_point(aes(x = cut_index, y = jx), color = "red", size = 1.5, shape = 16) +
      labs(x = "sample", y = expression(cm/s^3)) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      theme_minimal() +
      theme(text=element_text(size=7))
    return_plot <- list(list(xplot), list(vplot), list(aplot), list(jplot))
  
  }
  
  jx <- abs(jx)
  jy <- abs(jy)
  jz <- abs(jz) 
  
  matrix_jerk <- cbind(jx, jy, jz)
  sample_magnitudes <- sqrt(rowSums(matrix_jerk^2))
  
  if (PLOTmag) {
    fig_num <- fig_num + 1
    plot_title <- "Dir. Spat. Jerk Mag, Unfiltered"
    if (FILTER) {
      plot_title <- "Dir. Spat. Jerk Mag, Filtered"
    }
    if (print_all_R_output) plot_title <- paste("Fig ", sec_num, ".", fig_num, ": ", plot_title, sep = "")
    return_plot <- ggplot(data = data.frame(cut_index = cut_index, sample_magnitudes = sample_magnitudes)) +
      geom_point(aes(x = cut_index, y = sample_magnitudes), color = "red4", size = 2, shape = 16) +
      labs(x = "sample", y = expression(cm/s^3)) +
      theme_minimal() +
      theme(text=element_text(size=9), aspect.ratio = 1) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black")
  }
  
  output <- mean( sample_magnitudes, na.rm=TRUE ) 
  
  if ( return_output && (!PLOT || !PLOTmag) ) {
    return(output)
  } else if (PLOT || PLOTmag) {
    return(return_plot)
  } else if (print_all_R_output) {
    cat(c("Participant: ", p, "Reach: ", r))
  }
  
}
mean_reach_jerk_SU <- function(p,r,PLOT,FILTER=TRUE,forder,fcutoff,cutoff=cutoffSU,return_output=TRUE,PLOTmag=FALSE) {
  
  index <- which( motionSU$participant==p & motionSU$reach==r )
  
  qax <- motionSU$qax[index]
  qay <- motionSU$qay[index]
  qaz <- motionSU$qaz[index]
  qar <- motionSU$qar[index]
  
  qbx <- motionSU$qbx[index]
  qby <- motionSU$qby[index]
  qbz <- motionSU$qbz[index]
  qbr <- motionSU$qbr[index]
  
  qax_na <- which(is.na(qax))
  qay_na <- which(is.na(qay))
  qaz_na <- which(is.na(qaz))
  qar_na <- which(is.na(qar))
  
  qbx_na <- which(is.na(qbx))
  qby_na <- which(is.na(qby))
  qbz_na <- which(is.na(qbz))
  qbr_na <- which(is.na(qbr))
  
  na_index <- unique(c(qax_na,qay_na,qaz_na,qar_na,qbx_na,qby_na,qbz_na,qbr_na))
  na_index_low <- 0
  na_index_high <- 0
  
  # Filter
  if (FILTER) {
    
    if (length(na_index) >= 1) {
      if (any(na_index < cutoff)) na_index_low <- length(na_index < cutoff)
      if (any(na_index > (length(index)-cutoff))) na_index_high <- length(na_index > (length(index)-cutoff))
      qax <- qax[-na_index]
      qay <- qay[-na_index]
      qaz <- qaz[-na_index]
      qar <- qar[-na_index]
      qbx <- qbx[-na_index]
      qby <- qby[-na_index]
      qbz <- qbz[-na_index]
      qbr <- qbr[-na_index]
    }
    
    qax <- filter(butter(forder, fcutoff, type = "low"), qax)
    qay <- filter(butter(forder, fcutoff, type = "low"), qay)
    qaz <- filter(butter(forder, fcutoff, type = "low"), qaz)
    qar <- filter(butter(forder, fcutoff, type = "low"), qar)
    
    qbx <- filter(butter(forder, fcutoff, type = "low"), qbx)
    qby <- filter(butter(forder, fcutoff, type = "low"), qby)
    qbz <- filter(butter(forder, fcutoff, type = "low"), qbz)
    qbr <- filter(butter(forder, fcutoff, type = "low"), qbr)
    
  }
  
  # cutoff first and last 200ms of recording
  cut_index <- ((1+cutoff)-na_index_low):((length(qax)-cutoff)+na_index_high)
  
  qax <- qax[cut_index]
  qay <- qay[cut_index]
  qaz <- qaz[cut_index]
  qar <- qar[cut_index]
  
  qbx <- qbx[cut_index]
  qby <- qby[cut_index]
  qbz <- qbz[cut_index]
  qbr <- qbr[cut_index]
  
  vqax <- first_derivative(qax) * SU_sample_rate # uq/s
  vqay <- first_derivative(qay) * SU_sample_rate
  vqaz <- first_derivative(qaz) * SU_sample_rate
  vqar <- first_derivative(qar) * SU_sample_rate
  
  vqbx <- first_derivative(qbx) * SU_sample_rate 
  vqby <- first_derivative(qby) * SU_sample_rate
  vqbz <- first_derivative(qbz) * SU_sample_rate
  vqbr <- first_derivative(qbr) * SU_sample_rate
  
  aqax <- first_derivative(vqax) * SU_sample_rate # uq/s^2
  aqay <- first_derivative(vqay) * SU_sample_rate
  aqaz <- first_derivative(vqaz) * SU_sample_rate
  aqar <- first_derivative(vqar) * SU_sample_rate
  
  aqbx <- first_derivative(vqbx) * SU_sample_rate 
  aqby <- first_derivative(vqby) * SU_sample_rate
  aqbz <- first_derivative(vqbz) * SU_sample_rate
  aqbr <- first_derivative(vqbr) * SU_sample_rate
  
  jqax <- first_derivative(aqax) * SU_sample_rate # uq/s^3
  jqay <- first_derivative(aqay) * SU_sample_rate
  jqaz <- first_derivative(aqaz) * SU_sample_rate
  jqar <- first_derivative(aqar) * SU_sample_rate
  
  jqbx <- first_derivative(aqbx) * SU_sample_rate
  jqby <- first_derivative(aqby) * SU_sample_rate
  jqbz <- first_derivative(aqbz) * SU_sample_rate
  jqbr <- first_derivative(aqbr) * SU_sample_rate
  
  if (PLOT) {
    fig_num <<- fig_num + 1
    plot_title <- "Quat dimension qax, Unfiltered"
    if (FILTER) {
      plot_title <- "Quat dimension qax, Filtered"
    }
    if (print_all_R_output) plot_title <- paste("Fig ", sec_num, ".", fig_num, ": ", plot_title, sep = "")
    xplot <- ggplot(data = data.frame(cut_index = cut_index, qax = qax)) +
      geom_point(aes(x = cut_index, y = qax), color = "blue", size = 1.5, shape = 16) +
      labs(x = "sample", y = expression(uq)) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      theme_minimal() +
      theme(text=element_text(size=7))
    vplot <- ggplot(data = data.frame(cut_index = cut_index, vqax = vqax)) +
      geom_point(aes(x = cut_index, y = vqax), color = "gray", size = 1.5, shape = 16) +
      labs(x = "sample", y = expression(uq/s)) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      theme_minimal() +
      theme(text=element_text(size=7))
    aplot <- ggplot(data = data.frame(cut_index = cut_index, aqax = aqax)) +
      geom_point(aes(x = cut_index, y = aqax), color = "black", size = 1.5, shape = 16) +
      labs(x = "sample", y = expression(uq/s^2)) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      theme_minimal() +
      theme(text=element_text(size=7))
    jplot <- ggplot(data = data.frame(cut_index = cut_index, jqax = jqax)) +
      geom_point(aes(x = cut_index, y = jqax), color = "red", size = 1.5, shape = 16) +
      labs(x = "sample", y = expression(uq/s^3)) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      theme_minimal() +
      theme(text=element_text(size=7))
    return_plot <- list(list(xplot), list(vplot), list(aplot), list(jplot))
  }
  
  jqax <- abs(jqax)
  jqay <- abs(jqay)
  jqaz <- abs(jqaz) 
  jqar <- abs(jqar)
  
  jqbx <- abs(jqbx)
  jqby <- abs(jqby)
  jqbz <- abs(jqbz) 
  jqbr <- abs(jqbr)
  
  matrix_jerk <- cbind(jqax, jqay, jqaz, jqar, jqbx, jqby, jqbz, jqbr)
  sample_magnitudes <- sqrt(rowSums(matrix_jerk^2))
  
  if (PLOTmag) {
    fig_num <- fig_num + 1
    plot_title <- "Dir. Rot. Jerk Mag, Unfiltered"
    if (FILTER) {
      plot_title <- "Dir. Rot. Jerk Mag, Filtered"
    }
    if (print_all_R_output) plot_title <- paste("Fig ", sec_num, ".", fig_num, ": ", plot_title, sep = "")
    return_plot <- ggplot(data = data.frame(cut_index = cut_index, sample_magnitudes = sample_magnitudes)) +
      geom_point(aes(x = cut_index, y = sample_magnitudes), color = "red4", size = 2, shape = 16) +
      labs(x = "sample", y = expression(uq/s^3)) +
      theme_minimal() +
      theme(text=element_text(size=9), aspect.ratio = 1) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black")
  }

  output <- mean( sample_magnitudes, na.rm=TRUE ) 
  
  if ( return_output && (!PLOT || !PLOTmag) ) {
    return(output)
  } else if (PLOT || PLOTmag) {
    return(return_plot)
  } else if (print_all_R_output) {
    cat(c("Participant: ", p, "Reach: ", r))
  }
  
}

cat("\n")
cat("\nFilter Settings:\n") 
cat("\nSet cutoff for OT (Hz): ", filter_cutoff_OT)
cat("\nSet cutoff for SU (Hz): ", filter_cutoff_SU, "\n")
cat("\nSet cutoff for OT (NN): ", Nyquist_normalized_from_Hz(filter_cutoff_OT,OT_sample_rate))
cat("\nSet cutoff for SU (NN): ", Nyquist_normalized_from_Hz(filter_cutoff_SU,SU_sample_rate), "\n")
cat("\nSet order for OT: ", filter_order_OT)
cat("\nSet order for SU: ", filter_order_SU, "\n")
  
# Filling rTable with mean directional jerk magnitude numbers.
cat("\nFilling rTable with mean directional jerk magnitudes with set filter settings.\n")
fill_rTable_jerk(filter_order_OT,
                 Nyquist_normalized_from_Hz(filter_cutoff_OT,OT_sample_rate),
                 filter_order_SU,
                 Nyquist_normalized_from_Hz(filter_cutoff_SU,SU_sample_rate))

# Filling rTable with normalized mean directional jerk magnitude numbers
cat("\n")
cat("\nFilling rTable with normalized mean directional jerk magnitudes.\n")
fill_rTable_NormJerk()

cat("\n")
cat("\nMean normalized jerk difference (percentage points) between systems: ",jerk_diff(),"\n")

# Sanity check: Make sure (as expected) that, 
#   (1) there's a difference between taking the magnitude of position change, then the third derivative, vs 
#       taking the third derivative of each position axis, then taking jerk magnitude, and 
#   (2) the latter procedure results in higher values. 
compare_jerk_mag_first_and_end <- function(PLOT=FALSE,
                                           FILTER=TRUE,
                                           forder=filter_order_OT,
                                           fcutoff=Nyquist_normalized_from_Hz(filter_cutoff_OT,OT_sample_rate),
                                           cutoff=as.integer(cutoffSU/10)) {
  
  mag_first <- 0
  mag_end <- 0
  equal <- 0
  diff <- c()
  
  mf_all <- c()
  me_all <- c()
  
  cat("\nComparing taking magnetude of position first then differentiating three times,
      vs differentiating three times, then taking jerk magnetude (directional jerk) for OT data...\n")
  
  cat("\nParticipant: 1")
  
  for ( p in pTable$number ) {
    
    if (p!=1) cat(",", p)
    
    reaches <- rTable$ReachNum[which(rTable$Participants==p)]
    reaches <- reaches[which(!is.na(reaches))]
    
    for ( r in reaches ) {
      
      index <- which( motionOT$participant==p & motionOT$reach==r )
      
      x <- motionOT$x[index]
      y <- motionOT$y[index]
      z <- motionOT$z[index]
      
      sm <- sqrt(rowSums(cbind(x,y,z)^2))
      
      x_na <- which(is.na(x))
      y_na <- which(is.na(y))
      z_na <- which(is.na(z))
      
      na_index <- unique(c(x_na, y_na, z_na))
      na_index_low <- 0
      na_index_high <- 0
      
      # Filter
      if (FILTER) {
        
        if (length(na_index) >= 1) {
          if (any(na_index < cutoff)) na_index_low <- length(na_index < cutoff)
          if (any(na_index > (length(index)-cutoff))) na_index_high <- length(na_index > (length(index)-cutoff))
          x <- x[-na_index]
          y <- y[-na_index]
          z <- z[-na_index]
          sm <- sm[-na_index]
        }
        
        x <- filter(butter(forder, fcutoff, type = "low"), x)
        y <- filter(butter(forder, fcutoff, type = "low"), y)
        z <- filter(butter(forder, fcutoff, type = "low"), z)
        sm <- filter(butter(forder, fcutoff, type = "low"), sm)
        
      }
      
      # cutoff first and last 200ms of recording
      if ( length(x) != length(y) ||
           length(x) != length(z) ||
           length(z) != length(y) ) {
        cat("\n\nWarning! length of x, y, and z not all the same!\n\n")
      }
      
      cut_index <- ((1+cutoff)-na_index_low):((length(x)-cutoff)+na_index_high)

      x <- x[cut_index]
      y <- y[cut_index]
      z <- z[cut_index]
      sm <- sm[cut_index]
      
      vx <- first_derivative(x) * OT_sample_rate # cm/s
      vy <- first_derivative(y) * OT_sample_rate
      vz <- first_derivative(z) * OT_sample_rate
      vsm <- first_derivative(sm) * OT_sample_rate
      
      ax <- first_derivative(vx) * OT_sample_rate # cm/s^2
      ay <- first_derivative(vy) * OT_sample_rate
      az <- first_derivative(vz) * OT_sample_rate
      asm <- first_derivative(vsm) * OT_sample_rate
      
      jx <- first_derivative(ax) * OT_sample_rate # cm/s^3
      jy <- first_derivative(ay) * OT_sample_rate
      jz <- first_derivative(az) * OT_sample_rate
      jsm <- first_derivative(asm) * OT_sample_rate
      
      jx <- abs(jx)
      jy <- abs(jy)
      jz <- abs(jz) 
      
      matrix_jerk <- cbind(jx, jy, jz)
      sample_magnitudes <- sqrt(rowSums(matrix_jerk^2))
      
      if (PLOT) cat("\nNAs, mag first:", which(is.na(jsm)))
      if (PLOT) cat("\nNAs, mag end:", which(is.na(sample_magnitudes)))
      
      if (PLOT) plot(abs(jsm), col = "red", main="mag first")
      if (PLOT) plot(sample_magnitudes, col = "blue", main="mag end")
      
      # Plot the first set of points
      if (PLOT) plot(abs(jsm), col = "red", main = "Overlay of mag first and mag end")
      
      # Add the second set of points to the same plot
      if (PLOT) points(sample_magnitudes, col = "blue")
      
      # Add legend
      if (PLOT) legend("topright", legend = c("mag first", "mag end"), col = c("red", "blue"), pch = 1)
      
      if (PLOT) cat("\n\nMean, mag first:", mean(abs(jsm), na.rm=TRUE))
      if (PLOT) cat("\nMean, mag end:", mean(sample_magnitudes, na.rm=TRUE),"\n")
      
      mf <- mean(abs(jsm), na.rm=TRUE)
      me <- mean(sample_magnitudes, na.rm=TRUE)
      mf_all <- c( mf_all, mf )
      me_all <- c( me_all, me )
      
      diff <- c( diff, abs(mf - me))
      
      if ( mf > me ) {
        if (PLOT) cat("\n\nMag first larger.\n")
        mag_first <- mag_first + 1
      } else if ( mf == me ) {
        if (PLOT) cat("\n\nTwo equal.\n")
        equal <- equal + 1
      } else {
        if (PLOT) cat("\n\nMag end larger.\n")
        mag_end <- mag_end + 1
      }
      
      if (PLOT) cat("\nParticipant: ", p, "Reach: ", r)
      
    }
  }
  
  cat("\n\nNumber of reaches with larger mean jerk:\n") 
  cat("\nMag first:", mag_first)
  cat("\nMag end:", mag_end)
  cat("\nEqual:", equal)
  
  cat("\n\nMean mag first: ", mean(mf_all, na.rm = TRUE))
  cat("\nMean mag end: ", mean(me_all, na.rm = TRUE))
  
  cat("\n\n... as percent of mean mag first:")
  
  cat("\n\nMean difference: ", (mean(diff, na.rm = TRUE)/mean(mf_all, na.rm = TRUE))*100 )
  cat("\nsd difference: ", (sd(diff, na.rm = TRUE)/mean(mf_all, na.rm = TRUE))*100 )
  cat("\nMin difference: ", (min(diff, na.rm = TRUE)/mean(mf_all, na.rm = TRUE))*100 )
  cat("\nMax difference: ", (max(diff, na.rm = TRUE)/mean(mf_all, na.rm = TRUE))*100 )

}
compare_jerk_mag_first_and_end()
# Why is the mean jerk magnitude of a reach usually higher when magnitude of a sample is taken after
#   jerk is calculated along each dimension (instead of taking magnitude of each position sample and then
#   taking the jerk of that changing value)? 
# Answer: Because there can be changes in jerk along an axis that are offset by proportional changes in jerk
#   along another axis. 




