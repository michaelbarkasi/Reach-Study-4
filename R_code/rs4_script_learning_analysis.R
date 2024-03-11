
sample_learning_plots <- list() # Variable to collect sample learning plots for printing

# Learning Analysis: 
# Functions to fit exponential decay model to participant reach sequences, 
#   to calculate learning rate (LR), learning amount / asymptote (Asym), and saturation trial (SatT). 
cat("\nRunning learning analyis on study data.\n")

cat("\nSaturation factor: ", sat_factor)

# Initiate data frame to store learning variables 
part_LV_summary <- data.frame(
  Num = pTable$number,
  STE.L = rep(NA, length(pTable$number)),
  STE.A = rep(NA, length(pTable$number)),
  STE.S = rep(NA, length(pTable$number)),
  STE.MSE = rep(NA, length(pTable$number)),
  SPE.L = rep(NA, length(pTable$number)),
  SPE.A = rep(NA, length(pTable$number)),
  SPE.S = rep(NA, length(pTable$number)),
  SPE.MSE = rep(NA, length(pTable$number)),
  RTE.L = rep(NA, length(pTable$number)),
  RTE.A = rep(NA, length(pTable$number)),
  RTE.S = rep(NA, length(pTable$number)),
  RTE.MSE = rep(NA, length(pTable$number)),
  RPE.L = rep(NA, length(pTable$number)),
  RPE.A = rep(NA, length(pTable$number)),
  RPE.S = rep(NA, length(pTable$number)),
  RPE.MSE = rep(NA, length(pTable$number))
)

# Dependency functions: 
cat("\nLoading helper functions.\n")
compute_learning_curve <- function(trial_error,participant__) {
  
  # Remove model? Set when computing error
  index__ <- 1:last_no_feedback_trial_num
  if ( remove_model ) {
    m_num <- pTable$model[which(pTable$number==participant__)]
    index__ <- index__[-m_num]
  }
  
  # Compute initial error 
  IE <- mean(trial_error[index__],na.rm=TRUE)
  
  # Compute "learning" = improvement as a percentage of IE
  learning <- (IE - trial_error) / IE
  
  # Apply floor
  learning[which( learning < 0 )] <- 0
  
  return(learning)
  
}
# These functions adapted from the R code provided by Ruttle et al. 2021 (GNU GENERAL PUBLIC LICENSE v3), by Jennifer E. Ruttle, Bernard Marius 't Hart, and Denise Henriques
# EDF = Exponential Decay Function
EDF <- function(L,A,l) {
  
  # This function models an exponential decay learning process P with asymptote Asym, 
  #   where the value P_next of P on the next step is: 
  #   P_next = P_prev + L * ( Asym - P_prev )
  
  # the process and error states are initialized at 0:
  Pt <- 0
  Et <- 0
  
  # the total output is stored here:
  output <- c()
  
  for (t in c(1:l)) {
    
    # compute the process state at this step
    Pt <- Pt + (L * Et)
    # compute error for the next trial:
    Et <- A - Pt
    # save the process state in our vector:
    output <- c(output, Pt)
    
  }
  
  return(output)
  
}
MSE <- function(v,w) {
  
  mse <- mean((v-w)^2, na.rm=TRUE)
  
  return( mse )
  
}
MSE_EDF <- function(input,s) {
  
  l <- input[1]
  a <- input[2]
  output <- MSE(EDF(l,a,length(s)),s)
  
  return(output)
  
}
Fit_EDF <- function(learning,gridpoints=11,gridfits=10) {
  
  # set the search grid:
  gridbase <- seq(1/gridpoints/2,1-(1/gridpoints/2),1/gridpoints)
  
  #pos_learning <- c(learning[which(learning>0)],0)
  maxAsymptote <- 2*max(abs(learning), na.rm=TRUE)
  
  # define the search grid:
  searchgrid <- expand.grid('L' = gridbase, 
                            'A' = gridbase * maxAsymptote)
  
  MSEs <- c()
  
  # First pass, compute MSE for the coarse search grid just created
  for ( i in 1:nrow(searchgrid) ) {
    
    L <- searchgrid[i,'L']
    A <- searchgrid[i,'A']
    
    MSEi <- MSE_EDF(c(L,A),learning)
    
    MSEs <- c( MSEs, MSEi )
    
  }
  
  # Grab the gridfits best fits from the coarse search grid
  startpoints <- order(MSEs, decreasing = FALSE)[1:gridfits]
  searchgrid <- searchgrid[startpoints,]
  
  LA <- array(NA,dim=c(length(startpoints),3))
  
  for ( i in 1:length(startpoints) ) {
    
    L <- searchgrid[i,1]
    A <- searchgrid[i,2]
    
    suppressWarnings( # optimx can't compute gradient and throws warnings, but returns okay results
      results <- optimx(par = c(L,A), 
                        fn = MSE_EDF, 
                        s = learning, 
                        lower = c(0,0),
                        upper = c(1,maxAsymptote),
                        method = 'L-BFGS-B')
    )
    
    if ( results$convcode==0 ) { # if optimx is successful at converging, use its results
      LA[i,] <- c( results$p1[1], results$p2[1], results$value[1] )
    } else {
      LA[i,] <- c(L,A,MSE_EDF(c(L,A),learning))
    }
    
  }
  
  # Find row of L, A with smallest MSE
  min <- order(LA[,3], decreasing = FALSE)
  min <- min[1]
  
  best_L <- LA[min,1]
  best_A <- LA[min,2]
  
  return(c(best_L,best_A))
  
}

# Compute learning variables 
compute_learning_variables <- function() {
  
  cat("\nComputing learning variables for all participants, all types of error.\n")
  cat("\n")
  if ( !print_learning_results ) {
    cat("\n... Printing three sample participants: ", print_participants, "\n")
    cat("\n")
  }
  
  # Create new rows for pTable
  
  pTable$LR_STE <<- rep(NA,length(pTable$number))
  pTable$Asym_STE <<- rep(NA,length(pTable$number))
  pTable$SatT_STE <<- rep(NA,length(pTable$number))
  pTable$MSE_STE <<- rep(NA,length(pTable$number))
  
  pTable$LR_SPE <<- rep(NA,length(pTable$number))
  pTable$Asym_SPE <<- rep(NA,length(pTable$number))
  pTable$SatT_SPE <<- rep(NA,length(pTable$number))
  pTable$MSE_SPE <<- rep(NA,length(pTable$number))
  
  pTable$LR_RTE <<- rep(NA,length(pTable$number))
  pTable$Asym_RTE <<- rep(NA,length(pTable$number))
  pTable$SatT_RTE <<- rep(NA,length(pTable$number))
  pTable$MSE_RTE <<- rep(NA,length(pTable$number))
  
  pTable$LR_RPE <<- rep(NA,length(pTable$number))
  pTable$Asym_RPE <<- rep(NA,length(pTable$number))
  pTable$SatT_RPE <<- rep(NA,length(pTable$number))
  pTable$MSE_RPE <<- rep(NA,length(pTable$number))
  
  saved_learning <- array(NA,dim=c(length(pTable$number),length((last_no_feedback_trial_num+1):(pb_start-1))))
  saved_models <- array(NA,dim=c(length(pTable$number),length((last_no_feedback_trial_num+1):(pb_start-1))))
  saved_asymptote <- rep(NA,length(pTable$number))
  saved_sat_trial <- rep(NA,length(pTable$number))
  
  for ( i in 1:length(pTable$number) ) {
    
    p <- pTable$number[i]
    
    if ( print_learning_results || p %in% print_participants ) {
      cat("Participant: ", p, "\n")
    } 
    
    # Grab the error series for the participant's trials
    index <- which( rTable$Participants==p )
    STE <- rTable$SpatialTargetError[index]
    SPE <- rTable$SpatialPathError[index]
    RTE <- rTable$RotationTargetError[index]
    RPE <- rTable$RotationPathError[index]
    
    # Compute learning curve for that full error series
    STElearning <- compute_learning_curve(STE,p)
    SPElearning <- compute_learning_curve(SPE,p)
    RTElearning <- compute_learning_curve(RTE,p)
    RPElearning <- compute_learning_curve(RPE,p)
    
    # Cut learning curve down to just the feedback trials
    index <- (last_no_feedback_trial_num+1):(pb_start-1)
    STElearning <- STElearning[index]
    SPElearning <- SPElearning[index]
    RTElearning <- RTElearning[index]
    RPElearning <- RPElearning[index]
    
    # Find best fit EDF for each cut-down learning curve: 
    
    LA_STE <- Fit_EDF(STElearning)
    LA_SPE <- Fit_EDF(SPElearning)
    LA_RTE <- Fit_EDF(RTElearning)
    LA_RPE <- Fit_EDF(RPElearning)
    
    pTable$LR_STE[i] <<- LA_STE[1]
    pTable$LR_SPE[i] <<- LA_SPE[1]
    pTable$LR_RTE[i] <<- LA_RTE[1]
    pTable$LR_RPE[i] <<- LA_RPE[1]
    
    pTable$Asym_STE[i] <<- LA_STE[2]
    pTable$Asym_SPE[i] <<- LA_SPE[2]
    pTable$Asym_RTE[i] <<- LA_RTE[2]
    pTable$Asym_RPE[i] <<- LA_RPE[2]
    
    part_LV_summary$STE.L[i] <<- LA_STE[1]
    part_LV_summary$SPE.L[i] <<- LA_SPE[1]
    part_LV_summary$RTE.L[i] <<- LA_RTE[1]
    part_LV_summary$RPE.L[i] <<- LA_RPE[1]
    
    part_LV_summary$STE.A[i] <<- LA_STE[2]
    part_LV_summary$SPE.A[i] <<- LA_SPE[2]
    part_LV_summary$RTE.A[i] <<- LA_RTE[2]
    part_LV_summary$RPE.A[i] <<- LA_RPE[2]
    
    # Find saturation trial and MSE:
    
    # STE
    model <- EDF(LA_STE[1],LA_STE[2],length(STElearning))
    sat_STE <- Inf
    if ( length(which(model >= LA_STE[2] * sat_factor)) > 0 ) {
      sat_STE <- min(which(model >= LA_STE[2] * sat_factor),na.rm = TRUE)
    }
    mse_STE <- MSE(model,STElearning)
    saved_learning[i,] <- STElearning
    saved_models[i,] <- model
    saved_asymptote[i] <- LA_STE[2]
    saved_sat_trial[i] <- sat_STE
    
    # SPE
    model <- EDF(LA_SPE[1],LA_SPE[2],length(SPElearning))
    sat_SPE <- Inf
    if ( length(which(model >= LA_SPE[2] * sat_factor)) > 0 ) {
      sat_SPE <- min(which(model >= LA_SPE[2] * sat_factor),na.rm = TRUE)
    }
    mse_SPE <- MSE(model,SPElearning)
    
    # RTE
    model <- EDF(LA_RTE[1],LA_RTE[2],length(RTElearning))
    sat_RTE <- Inf
    if ( length(which(model >= LA_RTE[2] * sat_factor)) > 0 ) {
      sat_RTE <- min(which(model >= LA_RTE[2] * sat_factor),na.rm = TRUE)
    }
    mse_RTE <- MSE(model,RTElearning)
    
    # RPE
    model <- EDF(LA_RPE[1],LA_RPE[2],length(RPElearning))
    sat_RPE <- Inf
    if ( length(which(model >= LA_RPE[2] * sat_factor)) > 0 ) {
      sat_RPE <- min(which(model >= LA_RPE[2] * sat_factor),na.rm = TRUE)
    }
    mse_RPE <- MSE(model,RPElearning)
    
    # Save
    pTable$SatT_STE[i] <<- sat_STE
    pTable$SatT_SPE[i] <<- sat_SPE
    pTable$SatT_RTE[i] <<- sat_RTE
    pTable$SatT_RPE[i] <<- sat_RPE
    
    pTable$MSE_STE[i] <<- mse_STE
    pTable$MSE_SPE[i] <<- mse_SPE
    pTable$MSE_RTE[i] <<- mse_RTE
    pTable$MSE_RPE[i] <<- mse_RPE
    
    part_LV_summary$STE.S[i] <<- sat_STE
    part_LV_summary$SPE.S[i] <<- sat_SPE
    part_LV_summary$RTE.S[i] <<- sat_RTE
    part_LV_summary$RPE.S[i] <<- sat_RPE
    
    part_LV_summary$STE.MSE[i] <<- mse_STE
    part_LV_summary$SPE.MSE[i] <<- mse_SPE
    part_LV_summary$RTE.MSE[i] <<- mse_RTE
    part_LV_summary$RPE.MSE[i] <<- mse_RPE
    
    if ( print_learning_results || p %in% print_participants ) {
      cat("STE, L:", LA_STE[1] , "A:", LA_STE[2] , "S:", sat_STE , "MSE:", mse_STE , "\n")
      cat("SPE, L:", LA_SPE[1] , "A:", LA_SPE[2] , "S:", sat_SPE , "MSE:", mse_SPE , "\n")
      cat("RTE, L:", LA_RTE[1] , "A:", LA_RTE[2] , "S:", sat_RTE , "MSE:", mse_RTE , "\n")
      cat("RPE, L:", LA_RPE[1] , "A:", LA_RPE[2] , "S:", sat_RPE , "MSE:", mse_RPE , "\n")
    }
    
  }
  
  # Print plots
  sample_plot_num <- 0
  cat("\nPrinting plots for Spatial Target Error Learning:\n")
  for ( i in 1:length(pTable$number) ) {
    p <- pTable$number[i]
    if ( print_learning_plots || p %in% print_participants ) {
      
      fig_num <- fig_num + 1
      plot_title <- paste("Fig ", sec_num, ".", fig_num, ": Spatial Target Learning Curve, Participant: ", p, sep = "")
      plot(saved_learning[i,],
           pch=19,
           col="black",
           ylim=c(0,1),
           main=plot_title,
           xlab="feedback trial",
           ylab="learning") # plot results of STE so we can see
      points(saved_models[i,],col="blue",pch=19)
      abline(v = saved_sat_trial[i], col = "red", lty = 2)
      abline(h = saved_asymptote[i], col = "green", lty = 2)
      
      # Save sample nice plots
      if ( p %in% print_participants ) {
        
        sample_plot_num <- sample_plot_num + 1
        
        df__ <- data.frame(
          feedback_trial = 1:length(saved_learning[i, ])+25,
          learning = saved_learning[i, ],
          saved_models = saved_models[i, ]
        )
        
        i_asym <- saved_asymptote[i]
        if ( i_asym > 1 ) {
          i_asym <- 1.0
        }
        
        title__ <- paste("ST Learning Curve, Ex", sample_plot_num)
        
        plot_ <- ggplot(df__, aes(x = feedback_trial, y = learning)) +
          geom_point(shape = 19, color = "black", size = 0.75) +
          geom_point(aes(x = feedback_trial, y = saved_models), shape = 19, color = "blue", size = 0.75) +
          geom_vline(xintercept = saved_sat_trial[i]+last_no_feedback_trial_num, linetype = "dashed", color = "red") +
          geom_hline(yintercept = i_asym, linetype = "dashed", color = "green") +
          ylim(0, 1) +
          xlim(min(df__$feedback_trial)-1, max(df__$feedback_trial)+1) +
          labs(
            title = title__,
            x = "feedback trial",
            y = "learning"
          ) +
          theme_minimal() +
          theme(text=element_text(size=9))
        
        sample_learning_plots <<- c( sample_learning_plots, list(plot_) )
        
      }
    } 
    
  }
  
}
compute_learning_variables()
