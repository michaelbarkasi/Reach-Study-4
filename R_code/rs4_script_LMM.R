
# Jerk Stats
cat("\nSetting up stat analysis on normalized jerk and normalized error.\n")

cat("\nLoading packages.\n")
load_packages_LLM <- function() {
  
  cat("\nLoading tidyverse, supressing messages to allow PDF knitting.\n")
  suppressMessages(
    if(!require(tidyverse)) {
      install.packages("tidyverse")
      library(tidyverse)
    })
  if(!require(lme4)) {
    install.packages("lme4")
    library(lme4)
  }
  if(!require(boot)) {
    install.packages("boot")
    library(boot)
  }
  # For checking if t-distribution:
  if(!require(nortest)) {
    install.packages("nortest")
    library(nortest)
  }
  # For model diagnostics 
  if(!require(JWileymisc)) {
    install.packages("JWileymisc")
    library(JWileymisc)
  }
  if(!require(multilevelTools)) {
    install.packages("multilevelTools")
    library(multilevelTools)
  }
  # For plots
  if(!require(patchwork)) {
    install.packages("patchwork")
    library(patchwork)
  }
  
}
load_packages_LLM()

# Bonferroni Correction
if (Bonferroni_Correction) {
  
  # Check correlation between spatial and target variables.
  corJ <<- cor(na.omit(rTable$SUjerkNorm),na.omit(rTable$OTjerkNorm))
  corT <<- cor(na.omit(rTable$RotationTargetErrorNorm),na.omit(rTable$SpatialTargetErrorNorm))
  corP <<- cor(na.omit(rTable$RotationPathErrorNorm),na.omit(rTable$SpatialPathErrorNorm))
  cat("\nChecking correlation between spatial and target variables.")
  cat("\nNormalized Jerk:", corJ)
  cat("\nNormalized Target Error:", corT)
  cat("\nNormalized Path Error:", corP, "\n")
  
  # Adjust CI for multiple comparisons
  cat("\nAdjusting confidence inteveral to account for multiple comparisons.")
  cat("\nDesired level:", confidence_interval)
  confidence_interval <- 1-((1-confidence_interval)/(regression_coefficients*dependent_variables_tested))
  cat("\nAdjusted level:", confidence_interval, "\n")
  
}

# Helper functions
cat("\nLoading helper functions.\n")
`%!in%` <- function(x, table) !(x %in% table)
infer_F1T_effect_in_F2T <- function(F1T_effect_in_F2R, interaction_estimate) {
  
  # F1 = block (factor 1)
  # R = Sat (reference level)
  # T = Sw (treatment level)
  # F2 = learningfeedback (factor 2)
  # R = terminal (reference level)
  # T = online (treatment level)
  
  # Then...
  
  # F1T_effect_in_F2R = effect of Sw block for terminal learningfeedback
  # F2T_effect_in_F1R = effect of online learningfeedback in Sat block
  # and what we want is: effect of Sw block (F1T_effect) for online learningfeedback (_in_F2T)
  
  # Have: 
  F1T_effect_in_F2T <- F1T_effect_in_F2R + interaction_estimate
  # So, we can immediately infer what we want from the model output.
  
  return(F1T_effect_in_F2T)
  
}
plot_implied_means <- function(ef_df,dvN_,dv_) {
  
  plot_data <- data.frame(
    block = rep(c("Sat", "Sw"), each = 2),
    LearningFeedback = rep(c("Terminal", "Online"), 2),
    value = NA,
    real_value = NA,
    CI_lower = NA,
    CI_upper = NA
  )
  
  # Order of effects in ef_df: Intercept, Sw (Term), Online (Sat), Interaction, Sw (Online)
  plot_data$value[plot_data$block == "Sat" & 
                    plot_data$LearningFeedback == "Terminal"] <- ef_df$Estimate[which(ef_df$Effect == "Intercept")]
  plot_data$value[plot_data$block == "Sw" & 
                    plot_data$LearningFeedback == "Terminal"] <- ef_df$Estimate[which(ef_df$Effect == "Intercept")] + 
                                                                  ef_df$Estimate[which(ef_df$Effect == "Sw (Term)")]
  plot_data$value[plot_data$block == "Sat" & 
                    plot_data$LearningFeedback == "Online"] <- ef_df$Estimate[which(ef_df$Effect == "Intercept")] + 
                                                                ef_df$Estimate[which(ef_df$Effect == "Online (Sat)")]
  plot_data$value[plot_data$block == "Sw" & 
                    plot_data$LearningFeedback == "Online"] <- ef_df$Estimate[which(ef_df$Effect == "Intercept")] + 
                                                                ef_df$Estimate[which(ef_df$Effect == "Online (Sat)")] + 
                                                                ef_df$Estimate[which(ef_df$Effect == "Sw (Online)")]
  
  plot_data$real_value[plot_data$block == "Sat" & plot_data$LearningFeedback == "Terminal"] <- summary_stats$mean[which(summary_stats$Block == "Sat" & 
                                                                                                                        summary_stats$LearningFeedback == "terminal" &
                                                                                                                        summary_stats$variable == dv_)]
  plot_data$real_value[plot_data$block == "Sw" & plot_data$LearningFeedback == "Terminal"] <- summary_stats$mean[which(summary_stats$Block == "Sw" & 
                                                                                                                          summary_stats$LearningFeedback == "terminal" &
                                                                                                                          summary_stats$variable == dv_)]
  plot_data$real_value[plot_data$block == "Sat" & plot_data$LearningFeedback == "Online"] <- summary_stats$mean[which(summary_stats$Block == "Sat" & 
                                                                                                                          summary_stats$LearningFeedback == "online" &
                                                                                                                          summary_stats$variable == dv_)]
  plot_data$real_value[plot_data$block == "Sw" & plot_data$LearningFeedback == "Online"] <- summary_stats$mean[which(summary_stats$Block == "Sw" & 
                                                                                                                         summary_stats$LearningFeedback == "online" &
                                                                                                                         summary_stats$variable == dv_)]
  
  # Grab and adjust CI values
  plot_data$CI_lower[plot_data$block == "Sat" & 
                       plot_data$LearningFeedback == "Terminal"] <- ef_df$Lower_CI[which(ef_df$Effect == "Intercept")]
  plot_data$CI_upper[plot_data$block == "Sat" & 
                       plot_data$LearningFeedback == "Terminal"] <- ef_df$Upper_CI[which(ef_df$Effect == "Intercept")]
  
  plot_data$CI_lower[plot_data$block == "Sw" &
                       plot_data$LearningFeedback == "Terminal"] <- ef_df$Lower_CI[which(ef_df$Effect == "Intercept")] + 
    ef_df$Lower_CI[which(ef_df$Effect == "Sw (Term)")]
  plot_data$CI_upper[plot_data$block == "Sw" &
                       plot_data$LearningFeedback == "Terminal"] <- ef_df$Upper_CI[which(ef_df$Effect == "Intercept")] + 
    ef_df$Upper_CI[which(ef_df$Effect == "Sw (Term)")]
  
  plot_data$CI_lower[plot_data$block == "Sat" &
                       plot_data$LearningFeedback == "Online"] <- ef_df$Lower_CI[which(ef_df$Effect == "Intercept")] + 
    ef_df$Lower_CI[which(ef_df$Effect == "Online (Sat)")]
  plot_data$CI_upper[plot_data$block == "Sat" &
                       plot_data$LearningFeedback == "Online"] <- ef_df$Upper_CI[which(ef_df$Effect == "Intercept")] + 
    ef_df$Upper_CI[which(ef_df$Effect == "Online (Sat)")]
  
  plot_data$CI_lower[plot_data$block == "Sw" &
                       plot_data$LearningFeedback == "Online"] <- ef_df$Lower_CI[which(ef_df$Effect == "Intercept")] + 
    ef_df$Lower_CI[which(ef_df$Effect == "Online (Sat)")] + 
    ef_df$Lower_CI[which(ef_df$Effect == "Sw (Online)")]
  plot_data$CI_upper[plot_data$block == "Sw" &
                       plot_data$LearningFeedback == "Online"] <- ef_df$Upper_CI[which(ef_df$Effect == "Intercept")] + 
    ef_df$Upper_CI[which(ef_df$Effect == "Online (Sat)")] + 
    ef_df$Upper_CI[which(ef_df$Effect == "Sw (Online)")]
  
  # Plot
  fig_num <<- fig_num + 1
  final_plot_title <- paste("Fig ", sec_num, ".", subsec_num, ".", fig_num, ": ",
                         "Implied Means and Confidence Intervals", sep="")
  final_plot <- ggplot(plot_data, aes(x = block, y = value, color = LearningFeedback, group = LearningFeedback)) +
    geom_point(position = position_dodge(width = 0.0), size = 4) +
    geom_line(position = position_dodge(width = 0.0), linewidth = 1.5) +
    geom_errorbar(
      aes(ymin = CI_lower, ymax = CI_upper), 
      width = 0.2, 
      position = position_dodge(width = 0.0)) +
    labs(title = final_plot_title,
         x = "Block",
         y = dvN_) +
    scale_color_manual(values = c("Terminal" = C2_color, "Online" = C1_color)) +
    theme_minimal()
  
  print(final_plot)
  
  final_plot_title <- paste("Implied Means and CIs", sep="")
  final_plot <- ggplot(plot_data, aes(x = block, y = value, color = LearningFeedback, group = LearningFeedback)) +
    geom_point(position = position_dodge(width = 0.0), size = 4) +
    geom_line(position = position_dodge(width = 0.0), linewidth = 1.5) +
    geom_errorbar(
      aes(ymin = CI_lower, ymax = CI_upper), 
      width = 0.2, 
      position = position_dodge(width = 0.0)) +
    labs(title = final_plot_title,
         x = "Block",
         y = dvN_,
         color = "LFB") +
    scale_color_manual(values = c("Terminal" = C2_color, "Online" = C1_color)) +
    theme_minimal() + 
    theme(aspect.ratio = 1, legend.position = "bottom")
  
  BS_IM_plot_saved <<- final_plot
  
  # Save values 
  # colnames(plot_data) <- c( "Block", "Learning.Feedback", "Implied.Mean", "Real.Mean", "Implied.CI.Low", "Implied.CI.High" )
  saved_implied_means <<- plot_data
  
  cat("\nImplied means:\n")
  cat("\nSat, Term: ", plot_data$value[plot_data$block == "Sat" & plot_data$LearningFeedback == "Terminal"])
  cat("\nSw, Term: ", plot_data$value[plot_data$block == "Sw" & plot_data$LearningFeedback == "Terminal"])
  cat("\nSat, Online: ", plot_data$value[plot_data$block == "Sat" & plot_data$LearningFeedback == "Online"])
  cat("\nSw, Online: ", plot_data$value[plot_data$block == "Sw" & plot_data$LearningFeedback == "Online"])
  
  cat("\n\nImplied CIs for means:\n")
  cat("\nSat, Term: ", plot_data$CI_lower[plot_data$block == "Sat" & plot_data$LearningFeedback == "Terminal"], 
      ", ", plot_data$CI_upper[plot_data$block == "Sat" & plot_data$LearningFeedback == "Terminal"])
  cat("\nSw, Term: ", plot_data$CI_lower[plot_data$block == "Sw" & plot_data$LearningFeedback == "Terminal"], 
      ", ", plot_data$CI_upper[plot_data$block == "Sw" & plot_data$LearningFeedback == "Terminal"])
  cat("\nSat, Online: ", plot_data$CI_lower[plot_data$block == "Sat" & plot_data$LearningFeedback == "Online"],
      ", ", plot_data$CI_upper[plot_data$block == "Sat" & plot_data$LearningFeedback == "Online"])
  cat("\nSw, Online: ", plot_data$CI_lower[plot_data$block == "Sw" & plot_data$LearningFeedback == "Online"],
      ", ", plot_data$CI_upper[plot_data$block == "Sw" & plot_data$LearningFeedback == "Online"])
  
}

# Used for defining the four trial blocks
Sat_trial_replacement <- pb_start - avg_last_x_no_saturation
Sat_trial_replacement_low <- last_no_feedback_trial_num + avg_first_x_low_saturation
cat("\nDefining Sat trial replacement level:", Sat_trial_replacement)
cat("\nIf sat trial is higher than this or infinite, use this value as the sat trial.\n")
cat("\nDefining Sat trial replacement level (low):", Sat_trial_replacement_low)
cat("\nIf sat trial is lower than this, use this value as the sat trial.\n")

# Build reduced rTable for LMM
build_reduced_rTable <- function(sat_trial_vec,
                                 sat_type,
                                 BlockRef = "Sat", 
                                 LFeedbackRef = "online",
                                 FeedbackRef = "online",
                                 rescale_dv = FALSE,
                                 use_sat = TRUE,
                                 noisy_print = FALSE) {
  
  # Find rTable indexes
  index_Sat <- c()
  index_Sw <- c()
  
  for ( p in pTable$number ) {
    
    sat_trial <- sat_trial_vec[which(pTable$number==p)] + last_no_feedback_trial_num
    if ( is.infinite(sat_trial) || sat_trial > Sat_trial_replacement ) {
      sat_trial <- Sat_trial_replacement
    } else if ( sat_trial < Sat_trial_replacement_low ) {
      sat_trial <- Sat_trial_replacement_low
    }
    
    if (!use_sat) {
      sat_trial <- 51
    }
    
    index_Sat <- c( index_Sat, 
                    which(rTable$trials>=sat_trial & 
                            rTable$trials<pb_start & 
                            rTable$Participants==p & 
                            !is.na(rTable$SUjerkNorm) & 
                            !is.na(rTable$OTjerkNorm) & 
                            !is.na(rTable$RotationTargetErrorNorm) & 
                            !is.na(rTable$RotationPathErrorNorm) & 
                            !is.na(rTable$SpatialTargetErrorNorm) & 
                            !is.na(rTable$SpatialPathErrorNorm)) )
    index_Sw <- c( index_Sw, 
                   which(rTable$trials>=pb_start & 
                           rTable$Participants==p & 
                           !is.na(rTable$SUjerkNorm) & 
                           !is.na(rTable$OTjerkNorm) & 
                           !is.na(rTable$RotationTargetErrorNorm) & 
                           !is.na(rTable$RotationPathErrorNorm) & 
                           !is.na(rTable$SpatialTargetErrorNorm) & 
                           !is.na(rTable$SpatialPathErrorNorm)) )
    
  }
  
  # Print what we're doing and number of reaches found
  if (rescale_dv) {
    cat("\n\nBuilding reshaped rTable for lmer modelling, with rescaling of dv.\n")
  } else {
    cat("\n\nBuilding reshaped rTable for lmer modelling.\n")
  }
  cat("\nSaturation type used to define Sat block:", sat_type, "\n")
  cat("Total Sat reaches found: ", length(index_Sat), "\n")
  cat("Total Sw reaches found: ", length(index_Sw), "\n")
  cat("Average Sat reaches found per participant:", length(index_Sat)/length(pTable$number), "\n")
  cat("Average Sw reaches found per participant:", length(index_Sw)/length(pTable$number), "\n")
  
  num_Sat_reaches_found <<- length(index_Sat)
  num_Sw_reaches_found <<- length(index_Sw)
  avg_Sat_reaches_found_pp <<- length(index_Sat)/length(pTable$number)
  avg_Sw_reaches_found_pp <<- length(index_Sw)/length(pTable$number)
  
  # Reshape columns for new data frame
  NormalizedRotationalJerk <- c(rTable$SUjerkNorm[index_Sat], rTable$SUjerkNorm[index_Sw])
  NormalizedSpatialJerk <- c(rTable$OTjerkNorm[index_Sat], rTable$OTjerkNorm[index_Sw])
  NormalizedRotationalTargetError <- c(rTable$RotationTargetErrorNorm[index_Sat], rTable$RotationTargetErrorNorm[index_Sw])
  NormalizedRotationalPathError <- c(rTable$RotationPathErrorNorm[index_Sat], rTable$RotationPathErrorNorm[index_Sw])
  NormalizedSpatialTargetError <- c(rTable$SpatialTargetErrorNorm[index_Sat], rTable$SpatialTargetErrorNorm[index_Sw])
  NormalizedSpatialPathError <- c(rTable$SpatialPathErrorNorm[index_Sat], rTable$SpatialPathErrorNorm[index_Sw])
  Block <- c( rep("Sat", length(rTable$SUjerkNorm[index_Sat])), rep("Sw", length(rTable$SUjerkNorm[index_Sw])))
  LearningFeedback <- c(rTable$Condition[index_Sat], rTable$Condition[index_Sw])
  Participant <- c(rTable$Participants[index_Sat], rTable$Participants[index_Sw])
  
  Feedback <- c() 
  
  for ( i in 1:length(LearningFeedback) ) {
    
    if (LearningFeedback[i]=="online") {
      if (Block[i]=="Sat") {
        Feedback <- c( Feedback, "online" )
      } else {
        Feedback <- c( Feedback, "terminal" )
      }
    } else {
      if (Block[i]=="Sat") {
        Feedback <- c( Feedback, "terminal" )
      } else {
        Feedback <- c( Feedback, "online" )
      }
    }
  }
  
  # Rescale 
  if (rescale_dv) {
    
    NormalizedRotationalJerk <- scale(NormalizedRotationalJerk)
    NormalizedSpatialJerk <- scale(NormalizedSpatialJerk)
    NormalizedRotationalTargetError <- scale(NormalizedRotationalTargetError)
    NormalizedRotationalPathError <- scale(NormalizedRotationalPathError)
    NormalizedSpatialTargetError <- scale(NormalizedSpatialTargetError)
    NormalizedSpatialPathError <- scale(NormalizedSpatialPathError)
    
  }
  
  stratum <- rep(NA, length(NormalizedRotationalJerk))
  
  # Build new reshaped data frame
  lmm_df <- data.frame(NormalizedRotationalJerk,
                       NormalizedSpatialJerk,
                       NormalizedRotationalTargetError,
                       NormalizedRotationalPathError,
                       NormalizedSpatialTargetError,
                       NormalizedSpatialPathError,
                       Block,
                       LearningFeedback, 
                       Feedback, 
                       Participant,
                       stratum)
  
  # Ensure numeric
  lmm_df$NormalizedRotationalJerk <- as.numeric(lmm_df$NormalizedRotationalJerk)
  lmm_df$NormalizedSpatialJerk <- as.numeric(lmm_df$NormalizedSpatialJerk)
  lmm_df$NormalizedRotationalTargetError <- as.numeric(lmm_df$NormalizedRotationalTargetError)
  lmm_df$NormalizedRotationalPathError <- as.numeric(lmm_df$NormalizedRotationalPathError)
  lmm_df$NormalizedSpatialTargetError <- as.numeric(lmm_df$NormalizedSpatialTargetError)
  lmm_df$NormalizedSpatialPathError <- as.numeric(lmm_df$NormalizedSpatialPathError)
  
  # Convert to factors 
  lmm_df$Block <- as.factor(lmm_df$Block)
  lmm_df$LearningFeedback <- as.factor(lmm_df$LearningFeedback)
  lmm_df$Participant <- as.factor(lmm_df$Participant)
  lmm_df$Feedback <- as.factor(lmm_df$Feedback)
  
  lmm_df$Block <- relevel(lmm_df$Block, ref = BlockRef)
  lmm_df$LearningFeedback <- relevel(lmm_df$LearningFeedback, ref = LFeedbackRef)
  lmm_df$Feedback <- relevel(lmm_df$Feedback, ref = FeedbackRef)
  
  lmm_df$stratum <- interaction(lmm_df$Block, lmm_df$LearningFeedback)
  
  # Check
  o_Sat <- c()
  o_Sw <- c()
  t_Sat <- c()
  t_Sw <- c() 
  for ( p in pTable$number ) {
    
    n_Sat <- length(which(lmm_df$Block == "Sat" & lmm_df$Participant == p))
    n_Sw <- length(which(lmm_df$Block == "Sw" & lmm_df$Participant == p))
    
    if ( pTable$Condition[which(pTable$number==p)]=="online") {
      o_Sat <- c( o_Sat, n_Sat )
      o_Sw <- c( o_Sw, n_Sw )
    } else {
      t_Sat <- c( t_Sat, n_Sat )
      t_Sw <- c( t_Sw, n_Sw )
    }
    
    if (noisy_print) {
      cat("\nParticipant:", p, "\n")
      cat("Number of Sat reaches:", n_Sat, "\n")
      cat("Number of Sw reaches:", n_Sw, "\n")
    }
    
  }
  cat("\nMean num of Sat reaches in online:", mean(o_Sat))
  cat("\nMean num of Sw reaches in online:", mean(o_Sw))
  cat("\nMean num of Sat reaches in terminal:", mean(t_Sat))
  cat("\nMean num of Sw reaches in terminal:", mean(t_Sw))
  
  avg_num_of_Sat_online <<- mean(o_Sat)
  avg_num_of_Sw_online <<- mean(o_Sw)
  avg_num_of_Sat_terminal <<- mean(t_Sat)
  avg_num_of_Sw_terminal <<- mean(t_Sw)
  
  # Return reshaped dataframe
  return(lmm_df)
  
}

cat("\nTaking mean of Spatial Target Error and Rotation Target Error sat trials to define sat trial.\n")
pTable$SatT <- apply(data.frame(pTable$SatT_STE, pTable$SatT_RTE), 1, max, na.rm = TRUE)    

RrTable <- build_reduced_rTable(sat_trial_vec = pTable$SatT,
                                sat_type = "Max Spatial/Rotational Target Error",
                                BlockRef = "Sat", 
                                LFeedbackRef = "terminal",
                                FeedbackRef = "terminal",
                                rescale_dv = FALSE)

RrTable_rescaled <- build_reduced_rTable(sat_trial_vec = pTable$SatT,
                                         sat_type = "Max Spatial/Rotational Target Error",
                                         BlockRef = "Sat", 
                                         LFeedbackRef = "terminal",
                                         FeedbackRef = "terminal",
                                         rescale_dv = TRUE)

# Check summary stats of new table
cat("\n\nCompute summary stats (mean and sd) of all six dependent variables.\n")
cat("\n")
summary_stats<-RrTable %>%
  group_by(LearningFeedback,Block) %>%
  get_summary_stats(c(NormalizedRotationalJerk,
                      NormalizedSpatialJerk,
                      NormalizedRotationalTargetError,
                      NormalizedRotationalPathError,
                      NormalizedSpatialTargetError,
                      NormalizedSpatialPathError), type = "mean_sd")
print(summary_stats, n = 24)

# Load function to run LMM
cat("\nLoading function to run linear mixed-effects models.\n")
cat("\nNote: Optimizing with ML criterion, as we're primarily interested in fixed-effects.\n")
LMM_analysis <- function(dv,dvN,
                         iv_formula = IV_FORMULA,
                         num_bootstraps = num_of_bs,
                         REML_val = REML_VAL, # optimize with REML or ML criterion? Default is TRUE.
                         conf_level_ = confidence_interval, # confidence interval to compute for LMM fixed-effect estimates.
                         reference_block = "Sat", # "Sat" or "Sw"
                         rescale_dv_ = RESCALE_DV
                         ) {
  
  cat("\n*************************************************************************")
  cat("\n*************************************************************************")
  cat("\n*************************************************************************\n")
  
  # Set reference block 
  RrTable$Block <- relevel(RrTable$Block, ref = reference_block)
  RrTable_rescaled$Block <- relevel(RrTable_rescaled$Block, ref = reference_block)
  
  # Say what we're doing
  cat("\nRunning LMM analysis for variable:", dv, "\nReference block:", reference_block, "\nReference LearningFeedback: terminal\n")
  if (rescale_dv_) {
    cat("Rescaling of dv: TRUE\n")
  } else {
    cat("Rescaling of dv: FALSE\n")
  }
  
  # Set and save for use later if rescaling
  dv_sd <- sd(RrTable[,dv], na.rm = TRUE)
  dv_mean <- mean(RrTable[,dv], na.rm = TRUE)
  
  cat("\n*************************************************************************")
  cat("\n*************************************************************************\n")
  
  cat("\nPrinting prelim plot of data\n")
  color_mapping <- c("online" = C1_color, "terminal" = C2_color)
  
  # Prelim plot of data: 
  fig_num <- fig_num + 1
  prelim_plot <- ggplot(RrTable, aes(x=Block,y=.data[[dv]],color=LearningFeedback)) +
    geom_boxplot() +
    scale_color_manual(values = color_mapping) + 
    labs(
      title = paste("Fig ", sec_num, ".", subsec_num, ".", fig_num, ": ","Reach ", dvN, sep = ""),
      x = "Block",
      y = dvN
    ) 
  print(prelim_plot)
  prelim_plot_reaches <<- ggplot(RrTable, aes(x=Block,y=.data[[dv]],color=LearningFeedback)) +
    geom_boxplot() +
    scale_color_manual(values = color_mapping) + 
    labs(
      title = paste("Reach ", dvN, sep = ""),
      x = "Block",
      y = dvN,
      color = "LFB"
    ) + 
    theme_minimal()  +
    theme(legend.position = "bottom")
  
  cat("\n*************************************************************************")
  cat("\n*************************************************************************\n")
  
  # Construct model
  cat("\nModelling with both Block and LearningFeedback (Block*LFB).\n")
  model_data <- RrTable
  if (rescale_dv_) {
    cat("\nRescaling data.\n")
    model_data <- RrTable_rescaled
  }
  formula_block <- paste(dv, iv_formula)
  cat("\nModelling using formula:", formula_block, "\n")
  formula_block <- as.formula(formula_block)
  cat("\n")
  lmm_block <- lmer(formula_block, data = model_data, REML = REML_val,
                    control = lmerControl(optimizer = "nlminbwrap",
                                          check.conv.grad = .makeCC("warning", tol = fit_tolerance, relTol = NULL)))
  print(summary(lmm_block))
  
  # Save summary information in data frames
  # AIC tab
  AICtab_saved <<- as.data.frame(summary(lmm_block)$AICtab)
  AICtab_saved <<- t(AICtab_saved)
  AICtab_saved <<- as.data.frame(AICtab_saved)
  row.names(AICtab_saved) <<- NULL
  # Fixed Effects
  FE_saved <<- as.data.frame(summary(lmm_block)$coefficients)
  FE_saved <<- rownames_to_column(FE_saved, var = "Fixed Effect")
  FE_saved[1,1] <<- "$\\textrm{(Intercept)}$"
  row.names(FE_saved) <<- NULL
  # Scaled residuals:
  residual_summary <- summary(summary(lmm_block)$residuals)
  residual_summary_saved <<- data.frame(
    Min = as.numeric(residual_summary[1]),
    Q1 = as.numeric(residual_summary[2]),
    Median = as.numeric(residual_summary[3]),
    Q3 = as.numeric(residual_summary[5]),
    Max = as.numeric(residual_summary[6]),
    Mean = mean(summary(lmm_block)$residuals, na.rm = TRUE),
    SD = sd(summary(lmm_block)$residuals, na.rm = TRUE)
  )
  # Random Effects
  RanEff_saved <<- as.data.frame(VarCorr(lmm_block))
  RanEff_saved <<- RanEff_saved[-3]
  colnames(RanEff_saved) <<- c("Groups", "Name", "Variance", "Std.Dev.")
  # Correlation of Fixed Effects
  Cor_FE_saved <<- as.data.frame.matrix(summary(lmm_block)$vcov@factors$correlation)
  Cor_FE_saved <<- rownames_to_column(Cor_FE_saved, var = "Cor of FE")
  Cor_FE_saved[1,1] <<- "$\\textrm{(Intercept)}$"
  row.names(Cor_FE_saved) <<- NULL
  
  # Unscale residuals
  unscaled_res <<- residuals(lmm_block, scaled = FALSE)
  cat("\nUnscaled residuals, std.dv: ", sd(unscaled_res, na.rm = TRUE), "\n")
  
  # Summarize the fixed effects of interest: 
  if (rescale_dv_) {
    cat("\nFixed effects of interest (undoing rescaling):")
    intercept__ <- fixef(lmm_block)[1] * dv_sd + dv_mean
    F1T_effect_in_F2R_ <- fixef(lmm_block)[2] * dv_sd 
    F2T_effect_in_F1R_ <- fixef(lmm_block)[3] * dv_sd 
    interaction_estimate_ <- fixef(lmm_block)[4] * dv_sd 
    F1T_effect_in_F2T_ <- infer_F1T_effect_in_F2T(F1T_effect_in_F2R_, interaction_estimate_)
  } else {
    cat("\nFixed effects of interest:")
    intercept__ <- fixef(lmm_block)[1]
    F1T_effect_in_F2R_ <- fixef(lmm_block)[2]
    F2T_effect_in_F1R_ <- fixef(lmm_block)[3]
    interaction_estimate_ <- fixef(lmm_block)[4]
    F1T_effect_in_F2T_ <- infer_F1T_effect_in_F2T(F1T_effect_in_F2R_, interaction_estimate_)
  }
  cat("\nIntercept = mean terminal value in Sat block:", intercept__)
  cat("\nF2T_effect_in_F1R = effect of online learningfeedback in Sat block:", F2T_effect_in_F1R_)
  cat("\nF1T_effect_in_F2R = effect of Sw block for terminal learningfeedback:", F1T_effect_in_F2R_)
  cat("\nF1T_effect_in_F2T = effect of Sw block for online learningfeedback:", F1T_effect_in_F2T_)
  cat("\n")
  Fixed.Effect <- c("Int.", "Sw (Tm)", "On (Sat)", "Inter.X", "Sw (On)")
  FE.Estimate <- c(intercept__, F1T_effect_in_F2R_, F2T_effect_in_F1R_, interaction_estimate_, F1T_effect_in_F2T_)
  FEestimates_saved <<- data.frame(Fixed.Effect, FE.Estimate)
  
  # Extract predicted values
  cat("\nExtracting predicted values\n")
  model_data$actual <- model_data[,dv]
  if (rescale_dv_) {
    model_data$actual <- (model_data$actual * dv_sd) + dv_mean
  }
  model_data$PredictedValues_block <- predict(lmm_block, newdata = model_data)
  if (rescale_dv_) {
    model_data$PredictedValues_block <- (model_data$PredictedValues_block * dv_sd) + dv_mean
  }

  y_limits <- range(c(model_data$actual, 
                      model_data$PredictedValues_block))
  
  # Plot actual and LLM-predicted values 
  fig_num <- fig_num + 1
  color_mapping <- c("online" = C1_color, "terminal" = C2_color, "predicted" = "black")
  cat("\nPlotting actual and LLM-predicted normalized reach jerk...\n")
  combined_plot <- ggplot(model_data, aes(x = LearningFeedback, color = LearningFeedback)) +
    geom_point(aes(y = actual), position = position_jitter(width = 0.4), size = 2, alpha = 0.2) +
    geom_point(aes(y = PredictedValues_block, color = "predicted"), position = position_jitter(width = 0.15), alpha = 0.2) +
    facet_wrap(~Block) +
    scale_color_manual(values = color_mapping) + 
    labs(title = paste("Fig ", sec_num, ".", subsec_num, ".", fig_num, ": ",
                       "Actual and Predicted ", dvN, sep = ""),
         x = "Block",
         y = dvN) +
    theme_minimal() +
    theme(legend.position = "bottom") 
    coord_cartesian(ylim = y_limits)
  print(combined_plot)
  combined_plot <- ggplot(model_data, aes(x = LearningFeedback, color = LearningFeedback)) +
    geom_point(aes(y = actual), position = position_jitter(width = 0.4), size = 2, alpha = 0.2) +
    geom_point(aes(y = PredictedValues_block, color = "predicted"), position = position_jitter(width = 0.15), alpha = 0.2) +
    facet_wrap(~Block) +
    scale_color_manual(values = color_mapping) + 
    labs(title = "Actual and Predicted Values",
         x = "Block",
         y = dvN) +
    theme_minimal() +
    theme(legend.position = "bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
    labs(color = "LFB")  +
    coord_cartesian(ylim = y_limits)
  Predicted_values_plot_reaches_saved <<- combined_plot
  
  # Analyze participant means by block and LearningFeedback.
  
  cat("\n*************************************************************************")
  cat("\n*************************************************************************\n")
  
  cat("\nPlotting participant means by block and learning feedback...\n")
  
  cat("\nExtracting and printing participant mean values for each block.\n")
  
  # Extract participant mean values for each block
  mean_data <- aggregate(actual ~ Participant + Block + LearningFeedback,
                         data = model_data,
                         FUN = mean,
                         na.rm = TRUE)
  mean_data_predicted_block <- aggregate(PredictedValues_block ~ Participant + Block + LearningFeedback,
                         data =  model_data,
                         FUN = mean,
                         na.rm = TRUE)

  y_limits2 <- range(c(mean_data$actual, 
                       mean_data_predicted_block$PredictedValues_block))
  
  # Find MSE of model
  mse_model <- mean((mean_data_predicted_block$PredictedValues_block - mean_data$actual)^2, na.rm = TRUE)
  
  # Plot participant mean values connected by lines, colored by LearningFeedback 
  fig_num <- fig_num + 1
  color_mapping_LFB <- c("online" = C1_color_dark, "terminal" = C2_color_dark)
  mean_plot <- ggplot(mean_data, aes(x = Block, y = actual, group = Participant, color = LearningFeedback)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = color_mapping_LFB) + 
    labs(title = "Actual",
         x = "Block",
         y = paste("Mean", dvN)) +
    theme_minimal() +
    theme(legend.position='none') +
    coord_cartesian(ylim = y_limits2)
  mean_plot_p <- ggplot(mean_data_predicted_block, aes(x = Block, y = PredictedValues_block, group = Participant, color = LearningFeedback)) +
    geom_line(alpha = 0.7) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = color_mapping_LFB) + 
    labs(title = "Predicted",
         x = "Block",
         y = paste("Mean", dvN)) +
    theme_minimal() +
    theme(legend.position='none') +
    coord_cartesian(ylim = y_limits2)
  combined_mean_plot <- mean_plot + mean_plot_p +
    plot_layout(ncol = 2) +
    plot_annotation(title = paste("Fig ", sec_num, ".", subsec_num, ".", fig_num, ": ",
                                  "Participant Means by Block", sep = ""),
                    theme = theme_minimal(base_size = 15)) +
    theme(legend.position = "bottom") 
  print(combined_mean_plot)
  combined_mean_plot <- mean_plot +
    theme(legend.position = "bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
    labs(color = "LFB") + 
    mean_plot_p +
    plot_layout(ncol = 2) +
    plot_annotation(title = paste("Participant Means by Block", sep = ""),
                    theme = theme_minimal(base_size = 15)) 
  
  Predicted_values_plot_pmeans_saved <<- combined_mean_plot
  
  cat("\n*************************************************************************")
  cat("\n*************************************************************************\n")
  
  # Investigate normality
  cat("\nChecking normality of residuals.\n")
  fig_num <- fig_num + 1
  residuals <- resid(lmm_block)
  hist(residuals, main = paste("Fig ", sec_num, ".", subsec_num, ".", fig_num, ": ","LMM Residuals,", dvN), breaks = 20)
  # Shapiro-Wilk test for normality
  print(shapiro.test(residuals))
  # Investigate t-distribution
  cat("\nChecking t-distribution of residuals.\n")
  print(ad.test(residuals))
  
  # Run model diagnostics 
  cat("\nPlotting model diagnostics.\n")
  md <- modelDiagnostics(lmm_block, ev.perc = .001)
  md_plots <- plot(md, plot = FALSE, ask = FALSE)
  
  res_plots <- md_plots[1]
  ran_plots <- md_plots[2]
  
  ResFittedPlot_reshaped <- res_plots$Residuals$ResFittedPlot + 
    theme(aspect.ratio = 0.15)
  
  ResPlot_reshaped <- res_plots$Residuals$ResPlot + 
    theme(aspect.ratio = 0.75)
  
  fig_num <- fig_num + 1
  rp_combined <- res_plots$Residuals$ResPlot + res_plots$Residuals$ResFittedPlot +
    plot_layout(ncol = 2) +
    plot_annotation(title = paste("Fig ", sec_num, ".", subsec_num, ".", fig_num, ": ",
                                  "Residual Diagnostics", sep = ""),
                    theme = theme_minimal(base_size = 15))
  print(rp_combined)
  rp_combined <- ResPlot_reshaped + ResFittedPlot_reshaped +
    plot_layout(ncol = 2) +
    theme(aspect.ratio = 1) +
    plot_annotation(title = paste("Residual Diagnostics", sep = ""),
                    theme = theme_minimal(base_size = 12))
  
  residuals_diagnostics_plot_saved <<- rp_combined
  
  fig_num <- fig_num + 1
  ran_combined <- ran_plots$RandomEffects[[1]] + ran_plots$RandomEffects[[2]] + ran_plots$RandomEffects[[3]] +
    plot_layout(ncol = 2, guides = 'collect') +
    plot_annotation(title = paste("Fig ", sec_num, ".", subsec_num, ".", fig_num, ": ",
                                  "Random Variable Diagnostics", sep = ""),
                    theme = theme_minimal(base_size = 15))
  print(ran_combined)
  ran_combined <- ran_plots$RandomEffects[[1]] + ran_plots$RandomEffects[[2]] + ran_plots$RandomEffects[[3]] +
    plot_layout(ncol = 3, widths = c(2,2,3), guides = 'collect') +
    theme(aspect.ratio = 1) +
    plot_annotation(title = paste("Random Variable Diagnostics", sep = ""),
                    theme = theme_minimal(base_size = 12))
  
  random_variable_diagnostics_plot_saved <<- ran_combined
  
  cat("\nRunning model performance check.\n")
  print(modelPerformance(lmm_block))
  
  cat("\nThe data/residuals is not normal, nor t-distributed. Hence, we will use semiparametric and 
      nonparametric bootstrapping to estimate the fixed effects.\n")
  
  cat("\n*************************************************************************")
  cat("\n*************************************************************************\n")
  
  cat("\nComputing confidence intervals for the LMM estimates...\n")
  
  cat("\nConfidence Intervals (semiparametric bootstrap,", num_bootstraps, "resamples):\n")
  # See: https://www.rdocumentation.org/packages/lme4/versions/1.1-35.1/topics/bootMer
  lmm_block_CI <- confint(lmm_block,
                          level = conf_level_,
                          method = "boot",
                          nsim = num_bootstraps,
                          .progress = "txt",
                          use.u = TRUE, # uses the random effects from the model, doesn't randomize them
                          type = "semiparametric",
                          verbose = FALSE)
  if (rescale_dv_) {
    cat("\n\nNote: These are scaled values.\n\n")
  } else {
    cat("\n\n")
  }
  print(lmm_block_CI)
  
  # Use nonparametric bootstrapping to recalculate p-values of Fixed effects
  cat("\nConfidence Intervals (nonparametric bootstrap,", num_bootstraps, "resamples):\n")
  saved_mse <- c() # initiate to save mse of bootstrapped models
  
  # Function to fit the mixed-effects model to a bootstrap sample
  fit_model <- function(data, indices) {
    
    # Subset data to resamples
    bootstrap_data <- data[indices, ]
    
    if (rescale_dv_) {
      bs_dv_sd <- sd(bootstrap_data[,dv], na.rm = TRUE)
      bs_dv_mean <- mean(bootstrap_data[,dv], na.rm = TRUE)
      bootstrap_data[,dv] <- ( bootstrap_data[,dv] - bs_dv_mean ) / bs_dv_sd
    }
    
    # Fit model
    model <- lmer(formula_block, data = bootstrap_data, REML = REML_val,
                  control = lmerControl(optimizer = "nlminbwrap",
                                        check.conv.grad = .makeCC("warning", tol = fit_tolerance, relTol = NULL),
                                        calc.derivs = FALSE))
    
    # Extract fixed effects
    fixef_expanded <- fixef(model)
    if (rescale_dv_) {
      fixef_expanded <- fixef_expanded * bs_dv_sd
      fixef_expanded[1] <- fixef_expanded[1] + bs_dv_mean
    }
    F1T_effect_in_F2T__ <- infer_F1T_effect_in_F2T(fixef_expanded[2], fixef_expanded[4])
    fixef_expanded <- c( fixef_expanded, F1T_effect_in_F2T__ )
    
    # Find bootstrapped model mse for diagnostics 
    bootstrap_data$PredictedValues_block_bs <- predict(model, newdata = bootstrap_data)
    if (rescale_dv_) {
      bootstrap_data$PredictedValues_block_bs <- (bootstrap_data$PredictedValues_block_bs * bs_dv_sd) + bs_dv_mean
      bootstrap_data[,dv] <- (bootstrap_data[,dv] * bs_dv_sd) + bs_dv_mean
    }
    mean_data_bs <- aggregate(formula(paste(dv, "~ Participant + Block + LearningFeedback")),
                           data = bootstrap_data,
                           FUN = mean,
                           na.rm = TRUE)
    mean_data_predicted_block_bs <- aggregate(PredictedValues_block_bs ~ Participant + Block + LearningFeedback,
                                              data = bootstrap_data,
                                              FUN = mean,
                                              na.rm = TRUE)
    mse_model_bs <- mean((mean_data_predicted_block_bs$PredictedValues_block_bs - mean_data_bs[,dv])^2, na.rm = TRUE)
    saved_mse <<- c( saved_mse, mse_model_bs ) # Need to treat as global variable.
    
    return(fixef_expanded)  # Extract fixed effects estimates
    
  }

  # Apply bootstrapping
  # Note: lmm_df$stratum <- interaction(lmm_df$Block, lmm_df$LearningFeedback)
  bootstrap_results <- boot(data = RrTable, 
                            statistic = fit_model, 
                            R = num_bootstraps,
                            sim = "ordinary",
                            strata = RrTable$stratum)
  
  # Display raw bootstrap results
  print(bootstrap_results)
  
  # Save
  effect_names_raw <- c( "Intercept", "SwTm", "OnSat", "InterX", "SwOn" )
  saved_bs_results_raw <<- as.data.frame( bootstrap_results$t)
  colnames(saved_bs_results_raw) <<- effect_names_raw
  
  # Fill out bootstrap results with effect names and t-values
  treatment_block_name <- "Sat"
  if (reference_block=="Sat") treatment_block_name <- "Sw"
  effect_names <- c( paste("Inter.", sep = ""),
                     paste(treatment_block_name, " (Tm)", sep = ""),
                     paste("On (", reference_block, ")", sep = ""),
                     paste("InterX.", sep =""),
                     paste(treatment_block_name, " (On)", sep = "") )
  mean_bs_effects <- c() 
  bs_t_value <- c()
  fixed_effects_estimates <- fixef(lmm_block)
  if (rescale_dv_) {
    fixed_effects_estimates <- fixed_effects_estimates * dv_sd 
    fixed_effects_estimates[1] <- fixed_effects_estimates[1] + dv_mean
  }
  F1T_effect_in_F2T__ <- infer_F1T_effect_in_F2T(fixed_effects_estimates[2], fixed_effects_estimates[4])
  fixed_effects_estimates <- c( fixed_effects_estimates, F1T_effect_in_F2T__ )
  bs_original <- fixed_effects_estimates
  bs_bias <- c()
  bs_std_error <- c()
  p_observed <- c()
  p_estimated <- c()
  for (i in 1:5) {
    bs_r_mean <- mean(bootstrap_results$t[,i], na.rm = TRUE)
    #bs_r_sd <- sd(bootstrap_results$t[,i], na.rm = TRUE)
    bs_std_error_i <- sd(bootstrap_results$t[,i], na.rm = TRUE)
    bs_std_error <- c( bs_std_error, bs_std_error_i )
    bs_bias <- c( bs_bias, bs_r_mean - bs_original[i] )
    mean_bs_effects <- c( mean_bs_effects, bs_r_mean )
    bs_t_value <- c( bs_t_value, bs_r_mean / bs_std_error_i )
    p_observed <- c( p_observed, 1.0 - mean( abs(bs_original[i]) >= abs(bootstrap_results$t[,i]-bs_original[i]) ) )
    p_estimated <- c( p_estimated, 1.0 - mean( abs(bs_r_mean) >= abs(bootstrap_results$t[,i]-bs_r_mean) ) )
  }
  
  # Save these as a sanity check
  p_observed_uncorrected <<- p_observed
  p_estimated_uncorrected <<- p_estimated
  
  cat("\n")

  # Calculate confidence intervals based on bootstrap results
  lows_ <- c()
  highs_ <- c()
  for (i in 1:5) {
    param_ci <- boot.ci(bootstrap_results, conf = conf_level_, type = "perc", index = i)
    para_ <- paste("intercept (", reference_block, " BK, terminal LF)", sep = "")
    if ( i == 2 ) para_ <- paste("Bk ", treatment_block_name, " for terminal LFB", sep = "")
    if ( i == 3 ) para_ <- paste("LFB online in ", reference_block, " Block", sep = "")
    if ( i == 4 ) para_ <- paste("LFB online:Bk ", treatment_block_name, sep ="")
    if ( i == 5 ) para_ <- paste("Bk ", treatment_block_name, " for online LFB", sep = "")
    cat("Confidence intervals for", para_, ": ")
    low_ <- param_ci$percent[1,4]
    high_ <- param_ci$percent[1,5]
    lows_ <- c( lows_ , low_ )
    highs_ <- c( highs_, high_ )
    cat(low_, ", ", high_, "\n")
  }
  
  p_correction <- 1.0
  if (Bonferroni_Correction) {
    p_correction <- regression_coefficients*dependent_variables_tested
  }
  p_observed <- p_observed * p_correction
  p_estimated <- p_estimated * p_correction
  bootstrap_results_final <- data.frame( effect_names, 
                                         mean_bs_effects, 
                                         bs_original, 
                                         bs_bias, 
                                         bs_std_error, 
                                         bs_t_value, 
                                         lows_, 
                                         highs_,
                                         p_observed,
                                         p_estimated ) 

  new_names <- c( "Effect",
                  "BS.Estimate",
                  "Estimate",
                  "Bias",
                  "Std.Error",
                  "t.value",
                  "CI.Low",
                  "CI.High",
                  "p.Observed",
                  "p.Estimated" )
  
  # Rename bootstrap results
  colnames(bootstrap_results_final) <- new_names
  
  # Display the bootstrap results
  cat("\nPrinting final nonparametric bootstrapping results:\n")
  print(bootstrap_results_final)
  
  saved_bs_results <<- bootstrap_results_final

  cat("\nPlotting BS estimated fixed effects, with nonparametric bootstrapped confidence intervals.\n")
 
  effects_df <- data.frame(
    Effect = c( "Intercept", "Sw (Term)", "Online (Sat)", "Interaction", "Sw (Online)" ),
    Estimate = bootstrap_results_final$BS.Estimate,
    Lower_CI = lows_,
    Upper_CI = highs_
  )
  desired_order <- c("Intercept", "Sw (Term)", "Sw (Online)", "Online (Sat)", "Interaction")
  effects_df$Effect <- factor(effects_df$Effect, levels = desired_order)
  
  fig_num <- fig_num + 1
  CI_plot_title <- paste("Fig ", sec_num, ".", subsec_num, ".", fig_num, ": ",
                         "Estimated Fixed Effects, ", as.character(round(conf_level_,3)), "% CI, ", dvN, sep="")
  CI_plot <- ggplot(effects_df, aes(x = Effect, y = Estimate)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, position = position_dodge(width = 0.2)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0, hjust=0)) +
    labs(x = "Fixed Effects", y = "Estimate", title = CI_plot_title)

  print(CI_plot)
  
  CI_plot_title <- paste("Est. FE, ", as.character(round(conf_level_,3)), "% CI",  sep="")
  CI_plot <- ggplot(effects_df, aes(x = Effect, y = Estimate)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, position = position_dodge(width = 0.2)) +
    theme_minimal() +
    theme(aspect.ratio = 1) +
    theme(axis.text.x = element_text(angle = -90, vjust = 0, hjust=0)) +
    labs(x = "Fixed Effects", y = "Estimate", title = CI_plot_title) 
  
  BS_CI_plot_saved <<- CI_plot
  
  cat("\nFinding and plotting implied means and confidence intervals.\n")
  plot_implied_means(effects_df,dvN,dv)
  
  # Run bootstrap model diagnostics 
  cat("\n\nRunning bootstrap model diagnostics.\n")
  fig_num <- fig_num + 1
  BS_diag_title <- paste("Fig ", sec_num, ".", subsec_num, ".", fig_num, ": ",
                         "Bootstraped models, ", dvN, sep="")
  mean_mse_bs <- mean(saved_mse, na.rm = TRUE)
  sd_mse_bs <- sd(saved_mse, na.rm = TRUE)
  num_of_na <- length(which(is.na(saved_mse)))
  num_of_saved_mse <- length(saved_mse)
  cat("\nMSE of LMM of full data: ", mse_model)
  cat("\nNumber of saved mse for boostraped models: ", num_of_saved_mse)
  cat("\nNumber of NA in saved mse: ", num_of_na)
  cat("\nPercent NA in saved mse: ", (num_of_na/num_of_saved_mse)*100)
  cat("\nMean of mse for bootstrapped models: ", mean_mse_bs)
  cat("\nSD of mse for bootstrapped models: ", sd_mse_bs)
  hist(saved_mse, breaks = 100, main = BS_diag_title, xlab = "MSE", ylab = "Number of Models")
  abline(v = mse_model, col = "red", lty = 2)
  
  saved_mse_df <- data.frame(saved_mse)
  BS_diag_title <- paste("Bootstraped models, ", dvN, sep="")
  BS_MSE_diagnostics_plot_saved <<- ggplot(data = saved_mse_df, aes(x = saved_mse)) +
    geom_histogram(bins = 100, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = mse_model, color = "red", linetype = "dashed") +
    labs(
      title = BS_diag_title,
      x = "MSE",
      y = "Number of Models"
    ) +
    theme_minimal() + 
    annotate(
      "text",
      x = Inf,
      y = max(hist(saved_mse, breaks = 100, plot = FALSE)$counts)*0.975,
      label = paste("MSE full model =", format(mse_model, digits = 2)),
      hjust = 1,
      size = 5
    ) +
    annotate(
      "text",
      x = Inf,
      y = max(hist(saved_mse, breaks = 100, plot = FALSE)$counts)*0.90,
      label = paste("Mean MSE, BS models =", format(mean_mse_bs, digits = 2)),
      hjust = 1,
      size = 5
    ) +
    annotate(
      "text",
      x = Inf,
      y = max(hist(saved_mse, breaks = 100, plot = FALSE)$counts)*0.825,
      label = paste("SD MSE, BS models =", format(sd_mse_bs, digits = 2)),
      hjust = 1,
      size = 5
    )
  
}

# For testing and inspection
# DONT DELETE! SHOWS HOW TO RECONSTRUCT LME4 MODEL PREDICTIONS
Inspect_LMM_analysis <- function(dv,dvN,
                                 iv_formula = "~ Block * LearningFeedback + (Block | Participant)",
                                 REML_val = FALSE, # optimize with REML or ML criterion? Default is TRUE.
                                 reference_block = "Sat", # "Sat" or "Sw"
                                 rescale_dv_ = FALSE
) {
  
  # Set reference block 
  RrTable$Block <- relevel(RrTable$Block, ref = reference_block)
  RrTable_rescaled$Block <- relevel(RrTable_rescaled$Block, ref = reference_block)
  
  cat("\nRunning LMM analysis for variable:", dv, "\nReference block:", reference_block, "\nReference LearningFeedback: terminal\n")

  # Construct model
  cat("\nModelling with both Block and LearningFeedback (Block*LFB).\n")
  model_data <- RrTable
  formula_block <- paste(dv, iv_formula)
  cat("\nModelling using formula:", formula_block, "\n")
  formula_block <- as.formula(formula_block)
  cat("\n")
  lmm_block <- lmer(formula_block, data = model_data, REML = REML_val,
                    control = lmerControl(optimizer = "nlminbwrap",
                                          check.conv.grad = .makeCC("warning", tol = fit_tolerance, relTol = NULL)))
  print(summary(lmm_block))
  
  # Print fixed effects design matrix
  inspect_fixed_matrix <<- as.data.frame(as.matrix(model.matrix(lmm_block, type = "fixed")))
  inspect_fixed_matrix <<- getME(lmm_block, "X")
  
  #inspect_random_matrix <<- as.data.frame(as.matrix(model.matrix(lmm_block, type = "randomListRaw")))
  inspect_random_matrix <<- model.matrix(lmm_block, type = "randomListRaw")
  inspect_random_matrix <<- getME(lmm_block, "Z")
  
  inspect_random_effects <<- ranef(lmm_block)
  
  inspect_fixed_effects <<- fixef(lmm_block)
  
  inspect_residuals <<- resid(lmm_block)
  
  inspect_predictions <<- predict(lmm_block, newdata = model_data) + inspect_residuals
  
}
my_prediction <- function(r,print_values=FALSE) {
  
  p <- as.integer(RrTable$Participant[r])
  
  res <- inspect_residuals[r]
  
  beta_I <- inspect_fixed_effects[1]
  beta_ST <- inspect_fixed_effects[2]
  beta_OS <- inspect_fixed_effects[3]
  beta_X <- inspect_fixed_effects[4]
  
  BK <- inspect_fixed_matrix[r,2]
  LFB <- inspect_fixed_matrix[r,3]
  X <- inspect_fixed_matrix[r,4]
  
  #u_rI <- inspect_random_matrix$`Block | Participant`[r,1]
  #u_rST <- inspect_random_matrix$`Block | Participant`[r,2]
  
  u_rI <- as.matrix(inspect_random_matrix)[r,(2*p-1)]
  u_rST <- as.matrix(inspect_random_matrix)[r,(2*p)]
  
  rI <- inspect_random_effects$Participant[p,1]
  rBK <- inspect_random_effects$Participant[p,2]
  
  if (print_values) {
    cat("Values:\n")
    cat("p: ", p, "\n")
    cat("res: ", res, "\n")
    cat("beta_I: ", beta_I, "\n")
    cat("beta_ST: ", beta_ST, "\n")
    cat("beta_OS: ", beta_OS, "\n")
    cat("beta_X: ", beta_X, "\n")
    cat("BK: ", BK, "\n")
    cat("LFB: ", LFB, "\n")
    cat("X: ", X, "\n")
    cat("u_rI: ", u_rI, "\n")
    cat("u_rST: ", u_rST, "\n")
    cat("rI: ", rI, "\n")
    cat("rBK: ", rBK, "\n")
  }
  
  my_prediction <- ( beta_I + beta_ST * BK + beta_OS * LFB + beta_X * X ) + ( u_rI * rI + u_rST * rBK ) + res
  
  return(my_prediction)
  
}
check_all <- function() {
  
  for ( r_ in 1:nrow(inspect_fixed_matrix) ) {
    if ( round(inspect_predictions[r_], digits = 5) != round(my_prediction(r_), digits = 5) ) {
      cat("Error at row: ", r_, "\n")
    }
  }
  cat("\nCheck complete.\n")
  
}
#Inspect_LMM_analysis(dv = "NormalizedRotationalJerk", dvN = "Normalized Rotational Jerk")
#check_all()


