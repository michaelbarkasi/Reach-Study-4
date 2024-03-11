
# Setup: 

# Set the seed for as much reproducibility as possible
set.seed(9482, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

# Load packages
cat("Loading packages.\n")
load_packages <- function() {
  
  if(!require(hms)) {
    install.packages("hms")
    library(hms)
  }
  if(!require(readr)) {
    install.packages("readr")
    library(readr)
  }
  if(!require(dtw)) {
    install.packages("dtw")
    library(dtw)
  }
  if(!require(optimx)) {
    install.packages("optimx")
    library(optimx)
  }
  if(!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
  }
  if(!require(ggpubr)) {
    install.packages("ggpubr")
    library(ggpubr)
  }
  if (!require(rstatix)) {
    install.packages("rstatix")
    library(rstatix)
  }
  if(!require(plyr)) {
    install.packages("plyr")
    library(plyr)
  }
  if(!require(signal)) {
    install.packages("signal")
    library(signal)
  }

}
load_packages()

# Basic parameters of data set:
cat("\nLoading constants and other basic data set parameters.\n")
load_dataset_parameters <- function() {
  
  OT_sample_rate <<- 100 # OptiTrack samples at 100Hz
  SU_sample_rate <<- 1000 # Sonification system samples at 1000Hz
  last_no_feedback_trial_num <<- 25 # first how many reaches are random motion with no auditory feedback?
  avg_last_x_no_saturation <<- 10 # if participant did not saturate (learning), average last how many reaches to find Sat block
  avg_first_x_low_saturation <<- 10 # if participant saturates really fast ( < this number ), average first (this number) together for PreSat block
  pb_start <<- 76 # reach number which starts the post-block feedback (i.e., flipping online to terminal and vice-versa)
  significant_digitsOT <<- 1 # raw data is in meters and will be imported as cm, so this is down to the mm, mostly used in preprocessing
  significant_digitsSU <<- 5 # SU returns floats, 5 ensures all digits are significant, mostly used in preprocessing
  float__distance <<- 15 # if reach starts further than this (cm) from IP in OT data, discard as float
  IPrun <<- 150 # how many OT samples at beginning of OT run should be averaged to estimate IP?
  reach_min_length <<- 450 # SU recordings below this number of samples are likely twitches, not real reaches
  overrun_length <<- 1750 # SU records up to 2500 samples, but if more than this number is recorded, 
                          # participant likely failed to pause at end of reach, and reach
                          # end needs to be estimated.
  
  SU_corrupted_rows <<- c() # for collected corrupted rows in SU data
  
  FBcondition1 <<- "online" # in RS3, was "S" / "Sonification"
  FBcondition2 <<- "terminal" # in RS3, was "C" / "Faux Sonification" 
  pb_condition <<- "alt feedback" # in RS4, the post-feedback block is the opposite feedback, e.g., online switches to terminal
  NF_condition <<- "random"
  
  # Colors for plotting
  M_color <<- "orange" # model color
  R_color <<- "gray" # random reach color (first 25)
  R_color_dark <<- "gray"
  C1_color <<- "blue3" # feedback color (reaches 26-75) for condition 1
  C1_color_light <<- "azure"
  C1_color_dark <<- "blue4"
  C2_color <<- "red3" # feedback color (reaches 26-75) for condition 2
  C2_color_light <<- "lightpink"
  C2_color_dark <<- "red4"
  pb_color <<- "green4" # post-fb block color
  
}
load_dataset_parameters()

# Import data
cat("Importing Data.\n")
suppressWarnings(
motionSU <- read_csv("reach-study-4-SUdata.csv", 
                     col_types = cols(.default = col_double()))
)
motionOT <- read_csv("reach-study-4-OTdata.csv", 
                     col_types = cols(.default = col_double()))
pTable <- read_csv("reach-study-4-participant-table.csv", 
                   col_types = cols(Condition = col_character(), 
                                    sex = col_character(),
                                    .default = col_double()))

if (!public_data) {
  
  pTableD <- read_csv("pTable-demographic-data/reach-study-4-participant-table-demographics.csv", 
                     col_types = cols(Condition = col_character(), 
                                      sex = col_character(),
                                      .default = col_double()))
  
  # Participant demographics
  cat("\nDemographic Information:\n")
  
  count_M <- length(which(pTableD$sex == "M"))
  count_F <- length(which(pTableD$sex == "F"))
  mean_age <- mean(pTableD$age)
  min_age <- min(pTableD$age)
  max_age <- max(pTableD$age)
  std_dev_age <- sd(pTableD$age)
  
  # Print the counts for overall
  cat("\nMean age:", mean_age, "\n")
  cat("Minimum age:", min_age, "\n")
  cat("Maximum age:", max_age, "\n")
  cat("Standard deviation:", std_dev_age, "\n")
  cat("Number of females:", count_F, "\n")
  cat("Number of males:", count_M, "\n")
  
  count_M_online <- length(which(pTableD$sex[which(pTableD$Condition == "online")] == "M"))
  count_F_online <- length(which(pTableD$sex[which(pTableD$Condition == "online")] == "F"))
  mean_age_online <- mean(pTableD$age[which(pTableD$Condition == "online")])
  min_age_online <- min(pTableD$age[which(pTableD$Condition == "online")])
  max_age_online <- max(pTableD$age[which(pTableD$Condition == "online")])
  std_dev_age_online <- sd(pTableD$age[which(pTableD$Condition == "online")])
  
  # Print the counts for online 
  cat("\nMean age for online condition:", mean_age_online, "\n")
  cat("Minimum age for online condition:", min_age_online, "\n")
  cat("Maximum age for online condition:", max_age_online, "\n")
  cat("Standard deviation for age for online condition:", std_dev_age_online, "\n")
  cat("Number of females in online condition:", count_F_online, "\n")
  cat("Number of males in online condition:", count_M_online, "\n")
  
  count_M_terminal <- length(which(pTableD$sex[which(pTableD$Condition == "terminal")] == "M"))
  count_F_terminal <- length(which(pTableD$sex[which(pTableD$Condition == "terminal")] == "F"))
  mean_age_terminal <- mean(pTableD$age[which(pTableD$Condition == "terminal")])
  min_age_terminal <- min(pTableD$age[which(pTableD$Condition == "terminal")])
  max_age_terminal <- max(pTableD$age[which(pTableD$Condition == "terminal")])
  std_dev_age_terminal <- sd(pTableD$age[which(pTableD$Condition == "terminal")])
  
  # Print the counts for terminal
  cat("\nMean age for terminal condition:", mean_age_terminal, "\n")
  cat("Minimum age for terminal condition:", min_age_terminal, "\n")
  cat("Maximum age for terminal condition:", max_age_terminal, "\n")
  cat("Standard deviation for age for terminal condition:", std_dev_age_terminal, "\n")
  cat("Number of females in terminal condition:", count_F_terminal, "\n")
  cat("Number of males in terminal condition:", count_M_terminal, "\n")
  
}

# Clean sonification unit data by removing NA rows
remove_na_motionSU <- function() {
  
  columns_to_ignore <- c("high.gyros", "total.samples", "total.time", "rate") # These have lots of NAs in raw data, but not corrupted. 
  motionSU_temp <- motionSU[complete.cases(motionSU[, !names(motionSU) %in% columns_to_ignore]), ]
  num_of_NA_rows <- nrow(motionSU) - nrow(motionSU_temp)
  
  cat("\nRemoving rows with NA from motionSU: ", num_of_NA_rows, "\n")
  
  num_of_NA_rows_motionSU <<- num_of_NA_rows
  percentage_of_NA_rows_motionSU <<- (num_of_NA_rows / nrow(motionSU)) * 100
  
  motionSU <<- motionSU_temp
  
  # Not related to corrupted rows, but add the velocity column now.
  motionSU$RotVelocity <<- rep(NA,nrow(motionSU))
  
}
remove_na_motionSU() 

while ( length( problems(motionSU)$row ) > 0 ) {

  cat("\nCorrupted Rows found in Sonification System Data: ", length(problems(motionSU)$row), "/",
                                                              length(motionSU$reach), "=",
                                                              ( length(problems(motionSU)$row) / length(motionSU$reach) ) * 100, "%\n")

  cat("Corrupted Rows: ", problems(motionSU)$row, "\n")

  SU_corrupted_rows <<- unlist(problems(motionSU)$row)

  motionSU <<- motionSU[-SU_corrupted_rows,]

}

cat("\n")
cat("Loading functions for plotting data.\n")
# Plotting functions to check reaches after chopping
# Later number trials are darker
plotOTchopped <- function(ppp) {
  
  if ( !is.null(.check3d()) ) dev.off(which = .check3d())
  
  x_ <- c()
  y_ <- c()
  z_ <- c()
  
  model <- pTable$model[which(pTable$number==ppp)]
  
  firstcol <- R_color
  
  if ( pTable$model[which(pTable$number==ppp)]==1 ) {
    firstcol <- M_color
  }
  plot3d(
    motionOT$x[which(motionOT$reach==1 & motionOT$participant==ppp)],
    motionOT$y[which(motionOT$reach==1 & motionOT$participant==ppp)],
    motionOT$z[which(motionOT$reach==1 & motionOT$participant==ppp)],
    col = firstcol,xlab="x (forward)",ylab="y (vertical)",zlab="z (left-right)"
  )
  if ( pTable$model[which(pTable$number==ppp)]!=1 ) {
    x_ <- motionOT$x[motionOT$reach < model & motionOT$participant == ppp]
    y_ <- motionOT$y[motionOT$reach < model & motionOT$participant == ppp]
    z_ <- motionOT$z[motionOT$reach < model & motionOT$participant == ppp]
    plot3d(
      x_,
      y_,
      z_,
      col = R_color,add=TRUE
    )
  }
  plot3d(
    motionOT$x[which(motionOT$reach==model & motionOT$participant==ppp)],
    motionOT$y[which(motionOT$reach==model & motionOT$participant==ppp)],
    motionOT$z[which(motionOT$reach==model & motionOT$participant==ppp)],
    col = M_color,add=TRUE
  )
  x_ <- motionOT$x[motionOT$reach > model & motionOT$reach <= last_no_feedback_trial_num & motionOT$participant == ppp]
  y_ <- motionOT$y[motionOT$reach > model & motionOT$reach <= last_no_feedback_trial_num & motionOT$participant == ppp]
  z_ <- motionOT$z[motionOT$reach > model & motionOT$reach <= last_no_feedback_trial_num & motionOT$participant == ppp]
  plot3d(
    x_,
    y_,
    z_,
    col = R_color,add=TRUE
  )
  
  if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition1 ) {
    colorgradFunc <- colorRampPalette(c(C1_color_light,C1_color_dark))
    print_colors <- C1_color
  } else if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition2 ) {
    colorgradFunc <- colorRampPalette(c(C2_color_light,C2_color_dark))
    print_colors <- C2_color
  }
  
  x_ <- motionOT$x[motionOT$reach > last_no_feedback_trial_num & motionOT$reach < pb_start & motionOT$participant == ppp]
  y_ <- motionOT$y[motionOT$reach > last_no_feedback_trial_num & motionOT$reach < pb_start & motionOT$participant == ppp]
  z_ <- motionOT$z[motionOT$reach > last_no_feedback_trial_num & motionOT$reach < pb_start & motionOT$participant == ppp]
  colorgrad <- colorgradFunc(length(x_))
  plot3d(
    x_,
    y_,
    z_,
    col = colorgrad,add=TRUE
  )
  
  x_ <- motionOT$x[motionOT$reach >= pb_start & motionOT$participant == ppp]
  y_ <- motionOT$y[motionOT$reach >= pb_start & motionOT$participant == ppp]
  z_ <- motionOT$z[motionOT$reach >= pb_start & motionOT$participant == ppp]
  plot3d(
    x_,
    y_,
    z_,
    col = pb_color,add=TRUE
  )
  
  aspect3d(1, 1, 1)
  
}
plotOTchoppedReach <- function(ppp,rnum,with_model=FALSE,with_error_jerk=FALSE) {
  
  #if ( !is.null(.check3d()) ) dev.off(which = .check3d())
  model <- pTable$model[which(pTable$number==ppp)]
  
  if ( rnum == model ) {
    c <- M_color
  } else if ( rnum > last_no_feedback_trial_num ) {
    if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition1 ) {
      c <- C1_color
    } else if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition2 ) {
      c <- C2_color
    }
  } else {
    c <- R_color
  }
  
  c_ <- c
  x_ <- motionOT$x[which(motionOT$reach==rnum & motionOT$participant==ppp)]
  y_ <- motionOT$y[which(motionOT$reach==rnum & motionOT$participant==ppp)]
  z_ <- motionOT$z[which(motionOT$reach==rnum & motionOT$participant==ppp)]
  
  if (with_model) {
    mc <- c() 
    for ( i in which((motionOT$reach==rnum | motionOT$reach==model) & motionOT$participant==ppp) ) {
      if ( motionOT$reach[i] == rnum ) {
        mc <- c(mc,c)
      } else {
        mc <- c(mc,M_color)
      }
    }
    c_ <- mc
    x_ <- motionOT$x[which((motionOT$reach==rnum | motionOT$reach==model) & motionOT$participant==ppp)]
    y_ <- motionOT$y[which((motionOT$reach==rnum | motionOT$reach==model) & motionOT$participant==ppp)]
    z_ <- motionOT$z[which((motionOT$reach==rnum | motionOT$reach==model) & motionOT$participant==ppp)]
  } 
  
  plot3d(
    x_,y_,z_,
    col = c_, size = 7,
    xlab="x (forward)",ylab="y (vertical)",zlab="z (left-right)"
  )
  
  if (with_error_jerk) {
    x_ <- motionOT$x[which(motionOT$reach==rnum & motionOT$participant==ppp)]
    y_ <- motionOT$y[which(motionOT$reach==rnum & motionOT$participant==ppp)]
    z_ <- motionOT$z[which(motionOT$reach==rnum & motionOT$participant==ppp)]
    jerk_matrix <- compute_jerk_matrix(x_,y_,z_)
    s1 <- as.integer(nrow(jerk_matrix)*0.25)
    arrow_length <- 4
    arrow_components <- c(0,0,0)
    # Hard to see if/when there are big differences in jerk component magnitude
    arrow_components_x <- (jerk_matrix[s1,1]/max(c(abs(jerk_matrix[s1,1]),abs(jerk_matrix[s1,2]),abs(jerk_matrix[s1,3])), na.rm = TRUE))*arrow_length
    arrow_components_y <- (jerk_matrix[s1,2]/max(c(abs(jerk_matrix[s1,1]),abs(jerk_matrix[s1,2]),abs(jerk_matrix[s1,3])), na.rm = TRUE))*arrow_length
    arrow_components_z <- (jerk_matrix[s1,3]/max(c(abs(jerk_matrix[s1,1]),abs(jerk_matrix[s1,2]),abs(jerk_matrix[s1,3])), na.rm = TRUE))*arrow_length
    # better for demonstration purposes: 
    sorted <- order(jerk_matrix[s1,1:3])
    arrow_components[sorted[1]] <- arrow_length*0.67 # smallest component
    arrow_components[sorted[2]] <- arrow_length*0.85 # middle component
    arrow_components[sorted[3]] <- arrow_length # largest component
    thickness_ <- 2.5
    style_ <- "lines"
    b_ <- 100
    # Jerk magnitude
    arrow3d(p0 = c(x_[s1],y_[s1],z_[s1]),
            p1 = c(x_[s1]+arrow_components[1],
                   y_[s1]+arrow_components[2],
                   z_[s1]+arrow_components[3]),
            type = style_, s = 1/6, thickness = thickness_, n = b_, theta = pi/12, col = "cyan4"
    )
    # x, y, and z directional jerk
    arrow3d(p0 = c(x_[s1],y_[s1],z_[s1]),
            p1 = c(x_[s1]+ arrow_components[1],y_[s1],z_[s1]),
            type = style_, s = 1/6, thickness = thickness_, n = b_, theta = pi/12
    )
    arrow3d(p0 = c(x_[s1],y_[s1],z_[s1]),
            p1 = c(x_[s1],y_[s1]+ arrow_components[2],z_[s1]),
            type = style_, s = 1/6, thickness = thickness_, n = b_, theta = pi/12
    )
    arrow3d(p0 = c(x_[s1],y_[s1],z_[s1]),
            p1 = c(x_[s1],y_[s1],z_[s1]+ arrow_components[3]),
            type = style_, s = 1/6, thickness = thickness_, n = b_, theta = pi/12
    )
    d <- DTW_model_indexesOT_Vel(ppp,rnum)
    mi <- d$index1
    ri <- d$index2
    indices <- seq(from = 1, to = length(ri), length.out = 15)
    midpoint <- length(indices)/2
    for ( i in 1:length(indices) ) {
      if (indices[i] < midpoint ) {
        indices[i] <- as.integer(indices[i] + (midpoint-indices[i])/3)
      } else {
        indices[i] <- as.integer(indices[i] - (indices[i]-midpoint)/3)
      }
    }
    mx_ <- motionOT$x[which(motionOT$reach==model& motionOT$participant==ppp)]
    my_ <- motionOT$y[which(motionOT$reach==model & motionOT$participant==ppp)]
    mz_ <- motionOT$z[which(motionOT$reach==model & motionOT$participant==ppp)]
    ri_indices <- c(ri[indices], max(ri))
    mi_indices <- c(mi[indices], max(mi))
    xseg <- rbind(x_[ri_indices], mx_[mi_indices])
    yseg <- rbind(y_[ri_indices], my_[mi_indices])
    zseg <- rbind(z_[ri_indices], mz_[mi_indices])
    segments3d(x=xseg[,1:(ncol(xseg)-1)],y=yseg[,1:(ncol(yseg)-1)],z=zseg[,1:(ncol(zseg)-1)],col="gray",lwd=2)
    segments3d(x=xseg[,ncol(xseg)],y=yseg[,ncol(yseg)],z=zseg[,ncol(zseg)],col="black",lwd=3)
  }
  
}
plotSUReach <- function(ppp,rnum) {
  
  if ( !is.null(.check3d()) ) dev.off(which = .check3d())
  
  model <- pTable$model[which(pTable$number==ppp)]
  if ( rnum == model ) {
    c <- M_color
    c2 <- M_color
  } else if ( rnum > 25 ) {
    if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition1 ) {
      c <- C1_color
      c2 <- C1_color_dark
    } else if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition2 ) {
      c <- C2_color
      c2 <- C2_color_dark
    }
  } else {
    c <- R_color
    c2 <- R_color_dark
  }
  cl <- length(motionSU$qax[which(motionSU$reach==rnum & motionSU$participant==ppp)])
  cl2 <- length(motionSU$qbx[which(motionSU$reach==rnum & motionSU$participant==ppp)])
  plot3d(
    c(motionSU$qax[which(motionSU$reach==rnum & motionSU$participant==ppp)],
      motionSU$qbx[which(motionSU$reach==rnum & motionSU$participant==ppp)]),
    c(motionSU$qay[which(motionSU$reach==rnum & motionSU$participant==ppp)],
      motionSU$qby[which(motionSU$reach==rnum & motionSU$participant==ppp)]),
    c(motionSU$qaz[which(motionSU$reach==rnum & motionSU$participant==ppp)],
      motionSU$qbz[which(motionSU$reach==rnum & motionSU$participant==ppp)]),
    col = c( rep(c,cl), rep(c2,cl2) ),xlab="x",ylab="y",zlab="z"
  )
  
} # for plotting the "x/y/z" coordinates of quaternions (specific reach)
plotSU <- function(ppp) {
  
  if ( !is.null(.check3d()) ) dev.off(which = .check3d())
  
  model <- pTable$model[which(pTable$number==ppp)]
  
  plot3d(
    motionSU$qax[which(motionSU$reach==1 & motionSU$participant==ppp)],
    motionSU$qay[which(motionSU$reach==1 & motionSU$participant==ppp)],
    motionSU$qaz[which(motionSU$reach==1 & motionSU$participant==ppp)],
    col = R_color,xlab="x",ylab="y",zlab="z"
  )
  for ( i in 1:(model-1) ) {
    plot3d(
      motionSU$qax[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qay[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qaz[which(motionSU$reach==i & motionSU$participant==ppp)],
      col = R_color,add=TRUE
    )
  }
  plot3d(
    motionSU$qax[which(motionSU$reach==model & motionSU$participant==ppp)],
    motionSU$qay[which(motionSU$reach==model & motionSU$participant==ppp)],
    motionSU$qaz[which(motionSU$reach==model & motionSU$participant==ppp)],
    col = M_color,add=TRUE
  )
  for ( i in (model-1):last_no_feedback_trial_num ) {
    plot3d(
      motionSU$qax[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qay[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qaz[which(motionSU$reach==i & motionSU$participant==ppp)],
      col = R_color,add=TRUE
    )
  }
  
  if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition1 ) {
    colorgradFunc <- colorRampPalette(c(C1_color_light,C1_color_dark))
  } else if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition2 ) {
    colorgradFunc <- colorRampPalette(c(C2_color_light,C2_color_dark))
  }
  
  colorgrad <- colorgradFunc(length((last_no_feedback_trial_num+1):(pb_start-1)))
  
  for ( i in (last_no_feedback_trial_num+1):(pb_start-1) ) {
    plot3d(
      motionSU$qax[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qay[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qaz[which(motionSU$reach==i & motionSU$participant==ppp)],
      col = colorgrad[i-last_no_feedback_trial_num],add=TRUE
    )
  }
  
  for ( i in pb_start:max(motionSU$reach,na.rm=TRUE) ) {
    plot3d(
      motionSU$qax[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qay[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qaz[which(motionSU$reach==i & motionSU$participant==ppp)],
      col = pb_color,add=TRUE
    )
  }
  
  for ( i in 1:(model-1) ) {
    plot3d(
      motionSU$qbx[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qby[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qbz[which(motionSU$reach==i & motionSU$participant==ppp)],
      col = R_color,add=TRUE
    )
  }
  plot3d(
    motionSU$qbx[which(motionSU$reach==model & motionSU$participant==ppp)],
    motionSU$qby[which(motionSU$reach==model & motionSU$participant==ppp)],
    motionSU$qbz[which(motionSU$reach==model & motionSU$participant==ppp)],
    col = M_color,add=TRUE
  )
  for ( i in (model-1):last_no_feedback_trial_num ) {
    plot3d(
      motionSU$qbx[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qby[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qbz[which(motionSU$reach==i & motionSU$participant==ppp)],
      col = R_color,add=TRUE
    )
  }
  
  if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition1 ) {
    colorgradFunc <- colorRampPalette(c(C1_color_light,C1_color_dark))
  } else if ( pTable$Condition[which(pTable$number==ppp)] == FBcondition2 ) {
    colorgradFunc <- colorRampPalette(c(C2_color_light,C2_color_dark))
  }
  
  colorgrad <- colorgradFunc(length((last_no_feedback_trial_num+1):(pb_start-1)))
  
  for ( i in (last_no_feedback_trial_num+1):(pb_start-1) ) {
    plot3d(
      motionSU$qbx[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qby[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qbz[which(motionSU$reach==i & motionSU$participant==ppp)],
      col = colorgrad[i-last_no_feedback_trial_num],add=TRUE
    )
  }
  
  for ( i in pb_start:max(motionSU$reach,na.rm=TRUE) ) {
    plot3d(
      motionSU$qbx[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qby[which(motionSU$reach==i & motionSU$participant==ppp)],
      motionSU$qbz[which(motionSU$reach==i & motionSU$participant==ppp)],
      col = pb_color,add=TRUE
    )
  }
  
  aspect3d(1, 1, 1)
  
} # for plotting the "x/y/z" coordinates of quaternions
compute_jerk_matrix <- function(x,y,z,FILTER=TRUE) {
  
  filter_order_OT <- 9 
  filter_cutoff_OT <- 10 # in Hz
  filter_cutoff_SU <- 10
  
  forder=filter_order_OT
  fcutoff=Nyquist_normalized_from_Hz(filter_cutoff_OT,OT_sample_rate)
  cutoff=as.integer(cutoffSU/10)
  
  x_na <- which(is.na(x))
  y_na <- which(is.na(y))
  z_na <- which(is.na(z))
  
  na_index <- unique(c(x_na, y_na, z_na))

  # Filter
  if (FILTER) {
    
    if (length(na_index) >= 1) {
      x <- x[-na_index]
      y <- y[-na_index]
      z <- z[-na_index]
      sm <- sm[-na_index]
    }
    
    x <- filter(butter(forder, fcutoff, type = "low"), x)
    y <- filter(butter(forder, fcutoff, type = "low"), y)
    z <- filter(butter(forder, fcutoff, type = "low"), z)

  }
  
  vx <- first_derivative(x) * OT_sample_rate # cm/s
  vy <- first_derivative(y) * OT_sample_rate
  vz <- first_derivative(z) * OT_sample_rate
  
  ax <- first_derivative(vx) * OT_sample_rate # cm/s^2
  ay <- first_derivative(vy) * OT_sample_rate
  az <- first_derivative(vz) * OT_sample_rate
  
  jx <- first_derivative(ax) * OT_sample_rate # cm/s^3
  jy <- first_derivative(ay) * OT_sample_rate
  jz <- first_derivative(az) * OT_sample_rate
  
  jx <- abs(jx)
  jy <- abs(jy)
  jz <- abs(jz) 
  
  matrix_jerk <- cbind(jx, jy, jz)
  sample_magnitudes <- sqrt(rowSums(matrix_jerk^2))
  matrix_jerk <- cbind(matrix_jerk, sample_magnitudes)
  
  return(matrix_jerk)
  
} # Used in plotOTchoppedReach

# Dependency Functions
cat("Loading other dependency functions.\n")

# Misc Functions
vdist <- function(a,b) sqrt(sum((a - b)^2)) 
rolling_mean <- function(v, n, roundto) {
  
  if (n >= length(v)) {
    stop(cat("\n\nWarning! n must be less than the length of v, length v =", length(v), ", n = ", n))
  }
  
  output <- rep(NA, length(v))
  
  for ( i in 1:n ) {
    output[i] <- mean(v[1:i],na.rm=TRUE)
    if ( is.nan(output[i]) ) output[i] <- NA
  }
  
  for (i in (n+1):length(v)) {
    output[i] <- mean(v[((i-n)+1):i],na.rm=TRUE)
    if ( is.nan(output[i]) ) output[i] <- NA
  }
  
  output <- round(output,roundto)
  
  return(output)
  
}
model_check <- function(noisyprint) {
  if (noisyprint) cat("Model check started\n") 
  found <- FALSE
  for ( p in pTable$number ) {
    m <- pTable$model[which(pTable$number==p)]
    if (!(m %in% motionSU$reach[which(motionSU$participant==p)])) {
      if (noisyprint) cat("WARNING! model missing from participant:",p,"\n")
      found <- TRUE
    }
  }
  if (noisyprint) cat("Model check complete")
  if (found) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
checkifnumber <- function(x) {
  if ( is.na(x) || !is.numeric(x) || length(x) != 1 ) {
    cat("\n WARNING! expected number, but got NA, nonnumeric, or length", length(x),"\n")
  }
}
checkifnumericvector <- function(x,l) {
  if ( any(is.na(x)) || any(!is.numeric(x)) || length(x) != l ) {
    cat("\n WARNING! expected numeric vector of length", l,"but got NAs, nonnumeric items, or length", length(x),"\n")
  }
}

# For DTW
DTW_model_indexesOT_Vel <- function(p,r) {
  
  m_num <- pTable$model[which(pTable$number==p)]
  r1i <- which(motionOT$participant==p & motionOT$reach==r)
  rv <- motionOT$SpatialVelocity[r1i]
  mi <- which(motionOT$participant==p & motionOT$reach==m_num)
  mv <- motionOT$SpatialVelocity[mi]
  
  rv <- na.omit(rv)
  mv <- na.omit(mv)
  
  d <- dtw(mv,rv)
  
  # the reach velocity series is y (reference), the model velocity series is x (query).
  
  return(d)
  
} # For OptiTrack data: takes a participant and reach num 
                                               # and outputs the DTW structure of that reach (reference)
                                               # warped by velocity to the model (query).
DTW_model_indexesSU_Vel <- function(p,r) {
  
  m_num <- pTable$model[which(pTable$number==p)]
  r1i <- which(motionSU$participant==p & motionSU$reach==r)
  rv <- motionSU$RotVelocity[r1i]
  mi <- which(motionSU$participant==p & motionSU$reach==m_num)
  mv <- motionSU$RotVelocity[mi]
  
  rv <- na.omit(rv)
  mv <- na.omit(mv)
  
  d <- dtw(mv,rv)
  
  # the reach velocity series is y (reference), the model velocity series is x (query).
  
  return(d)
  
} # ibid, for Son System data


