
##################################################################################
##################################################################################
# Code to simulate RS4 data using Linear Mixed-Effects Modelling (LMM)
cat("\nSimulating RS4 data using Linear Mixed-Effects Modelling (LMM).\n")

# libraries
cat("\nLoading necessary libaries.\n")
running_alone <- TRUE # Set to TRUE if running this script alone, FALSE if running from main RS4 markdown file
if ( exists("rTable") ) running_alone <- FALSE
if (running_alone) {
  
  if(!require(lme4)) {
    install.packages("lme4")
    library(lme4)
  }
  if(!require(boot)) {
    install.packages("boot")
    library(boot)
  }
  if(!require(patchwork)) {
    install.packages("patchwork")
    library(patchwork)
  }
  if(!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
  }
  if(!require(JWileymisc)) {
    install.packages("JWileymisc")
    library(JWileymisc)
  }
  if(!require(multilevelTools)) {
    install.packages("multilevelTools")
    library(multilevelTools)
  }
  # For checking if t-distribution:
  if(!require(nortest)) {
    install.packages("nortest")
    library(nortest)
  }
  
  set.seed(19293)
  fig_num <- 0
  sec_num <- 0
  subsec_num <- 0
  
}
if(!require(detectnorm)) {
  install.packages("detectnorm")
  library(detectnorm)
}

##################################################################################
##################################################################################
# Part 1: Simulated data using LMM

# Simulation parameters
num_trials <- 50 # number of trials to simulate (half will be F1R and half F1T)
normality_failure <- TRUE
p_num <- 70 # Number of participants to simulate
cat("\nSimulation parameters:\n")
cat("Number of trials: ", num_trials, "\n")
cat("Normality failure: ", normality_failure, "\n")
cat("Number of participants: ", p_num, "\n")

# Factor levels
F1R <- "Sat" # Factor 1, reference level
F1T <- "Sw" # Factor 1, treatment level
F2R <- "terminal" # Factor 2, reference level
F2T <- "online" # Factor 2, treatment level

# Factors 
F1 <- c( F1R, F1T )
F2 <- c( F2R, F2T )

# Set fixed effects (stipulate values)
F1RF2R <- -0.49 # Intercept
F1TeF2R <- -0.0375 # effect of F1T on F2R
F2TeF1R <- -0.0541 # effect of F2T on F1R
IntX <- -0.0091 # Interaction
cat("\nStipulated fixed effects:\n")
cat("F1RF2R: ", F1RF2R, "\n")
cat("F1TeF2R: ", F1TeF2R, "\n")
cat("F2TeF1R: ", F2TeF1R, "\n")
cat("IntX: ", IntX, "\n")

# Set random effects 
p_variation_intercept <- 0.267 # Standard deviation of random effects, taken from real study data for Spatial Target Error, for intercept (F1R)
p_variation_slope <- 0.1372 # Standard deviation of random effects, taken from real study data for Spatial Target Error, for slope
p_intercept_skew <- 10
p_intercept_kurt <- 0
p_slope_skew <- p_intercept_skew*2
p_slope_kurt <- p_intercept_kurt*2
R_intercept <- rnorm(p_num, 0, p_variation_intercept) 
R_slope <- rnorm(p_num, 0, p_variation_slope) 
if ( normality_failure ) {
  
  R_intercept_skew <- rnonnorm(p_num, mean = 0, sd = p_variation_intercept, skew = p_intercept_skew, kurt = p_intercept_kurt)$dat
  R_intercept_T <- rt(p_num, 2.15)/10 # This value was hand-tuned to have a std.dv matching p_variation_skew
  
  R_slope_skew <- rnonnorm(p_num, mean = 0, sd = p_variation_slope, skew = p_slope_skew, kurt = p_slope_kurt)$dat
  R_slope_T <- rt(p_num, 5)/10 # This value was hand-tuned to have a std.dv matching p_variation_slope
  
  R_intercept <- R_intercept_skew * runif(p_num, 3, 5) + R_intercept_T * 0.5
  R_slope <- R_slope_skew + R_slope_T
  
  R_intercept <- (R_intercept / sd(R_intercept, na.rm = TRUE)) * p_variation_intercept
  R_slope <- (R_slope / sd(R_slope, na.rm = TRUE)) * p_variation_slope
  
  # plot(density(R_intercept)) # Check distribution
  
} 
cat("\nStipulated random effect distribution:\n")
cat("p_variation_intercept (std.dv): ", p_variation_intercept, "\n")
cat("p_variation_slope (std.dv): ", p_variation_slope, "\n")
if ( normality_failure ) {
  cat("p_intercept_skew: ", p_intercept_skew, "\n")
  cat("p_intercept_kurt: ", p_intercept_kurt, "\n")
  cat("p_slope_skew: ", p_slope_skew, "\n")
  cat("p_slope_kurt: ", p_slope_kurt, "\n")
}

# Generate residuals
res_sd <- 0.221 # standard deviation of residuals, taken from real study data for Spatial Target Error
res_skew <- 0
res_kurt <- 100
res_cap <- 7
residuals <- rnorm(num_trials*p_num, 0, res_sd)
if ( normality_failure ) {
  residuals <- rt(num_trials*p_num, 5) # This value was hand-tuned to have a std.dv matching res_sd
  residuals[which(residuals > 0)] <- residuals[which(residuals > 0)] + (residuals[which(residuals > 0)]*0.5)^2 # fatten tails
  residuals[which(residuals < 0)] <- residuals[which(residuals < 0)] - (residuals[which(residuals < 0)]*0.5)^2 # fatten tails
  residuals[which(residuals > res_cap)] <- res_cap # Cap 
  residuals[which(residuals < -res_cap)] <- -res_cap # Cap 
  residuals <- (residuals / sd(residuals, na.rm = TRUE)) * res_sd # normalize to desired std.dv
} 
cat("\nStipulated residual distribution:\n")
cat("res_sd (std.dv): ", res_sd, "\n")
if ( normality_failure ) {
  cat("res_skew: ", res_skew, "\n")
  cat("res_kurt: ", res_kurt, "\n")
}
# Shapiro-Wilk test for normality
# cat("\nChecking normality of residuals.\n")
print(shapiro.test(residuals))
# Investigate t-distribution
cat("\nChecking t-distribution of residuals.\n")
print(ad.test(residuals))

# Assign even number participants to F2R and odd number participants to F2T
F2_participants <- rep( c( F2T, F2R ), p_num/2 )

# Build fixed-effects design matrix
X_design <- c() 
for ( p in 1:p_num ) {
  for ( t in 1:(num_trials/2) ) {
    F1 <- 0
    F2 <- 0 
    if ( F2_participants[p] == F2T ) F2 <- 1
    X_design <- c( X_design, 1, F1, F2, F1*F2 )
  }
  for ( t in (num_trials/2+1):num_trials ) {
    F1 <- 1
    F2 <- 0
    if ( F2_participants[p] == F2T ) F2 <- 1
    X_design <- c( X_design, 1, F1, F2, F1*F2 )
  }
}
X_design <- matrix( data = X_design,
                    nrow = num_trials*p_num,
                    ncol = 4,
                    byrow = TRUE )

# Build random-effects design matrix
Z_design <- c()
tn <- 0 
for ( p in rep( 1:p_num, each = num_trials ) ) {
  tn <- tn+1
  row <- rep( 0, 2*p_num )
  row[(2*p)-1] <- 1
  if ( tn > 25 ) row[2*p] <- 1
  Z_design <- c( Z_design, row )
  if ( tn == num_trials ) tn <- 0
}
Z_design <- matrix( data = Z_design,
                    nrow = num_trials*p_num,
                    ncol = 2*p_num,
                    byrow = TRUE )

# Set effects vectors
Fixed_effects <- c( F1RF2R, F1TeF2R, F2TeF1R, IntX )
Random_effects <- c(rbind( R_intercept, R_slope ))

# Generate data
cat("\nGenerating simulated data using LMM.\n")
simulate_missed_reaches <- TRUE
Participant <- rep( 1:p_num, each = num_trials )
LearningFeedback <- c() 
Block <- c()
for ( i in 1:(p_num*num_trials) ) {
  if ( X_design[i,3] == 0 ) LearningFeedback <- c( LearningFeedback, F2R )
  else LearningFeedback <- c( LearningFeedback, F2T )
  if ( X_design[i,2] == 0 ) Block <- c( Block, F1R )
  else Block <- c( Block, F1T )
}
LearningFeedback <- as.factor(LearningFeedback)
Block <- as.factor(Block)
NormalizedSpatialTargetError <- X_design %*% Fixed_effects + Z_design %*% Random_effects + residuals
noise <- runif(length(which(NormalizedSpatialTargetError < -1)), min = 0.01, max = 0.25)
NormalizedSpatialTargetError[which(NormalizedSpatialTargetError < -1)] <- noise - 1
if ( simulate_missed_reaches ) {
  p_missed <- rnorm(p_num, mean = 0.1, sd = 0.05)
  p_missed[which(p_missed < 0)] <- 0
  missed_index <- c() 
  for ( p_ in 1:p_num ) {
    missed_index <- c( missed_index, sample( which(Participant == p_), as.integer(p_missed[p_]*num_trials), replace = FALSE ) )
  }
  NormalizedSpatialTargetError <- NormalizedSpatialTargetError[-missed_index]
  Participant <- Participant[-missed_index]
  LearningFeedback <- LearningFeedback[-missed_index]
  Block <- Block[-missed_index]
}
NormalizedSpatialTargetError_sd <- sd(NormalizedSpatialTargetError, na.rm = TRUE) # In case scaling is used, to unscale
NormalizedSpatialTargetError_mean <- mean(NormalizedSpatialTargetError, na.rm = TRUE) # In case scaling is used, to unscale

# Format data in data frame for LMM fit
stratum <- rep(NA, length(NormalizedSpatialTargetError))
ndata <- data.frame(NormalizedSpatialTargetError, Participant, Block, LearningFeedback, stratum)
ndata$stratum <- interaction(ndata$Block, ndata$LearningFeedback)
ndata$Block <- relevel(ndata$Block, ref = F1R)
ndata$LearningFeedback <- relevel(ndata$LearningFeedback, ref = F2R)

##################################################################################
##################################################################################
# Part 2: Model data and evaluate model

# LMM modelling parameters
form <- "NormalizedSpatialTargetError ~ Block * LearningFeedback + (Block | Participant)"
reml <-FALSE
rescale <- FALSE
if ( rescale ) ndata$NormalizedSpatialTargetError <- scale(ndata$NormalizedSpatialTargetError)
cat("\nLMM modelling parameters:\n")
cat("LMM formula: ", form, "\n")
cat("REML: ", reml, "\n")

# Model data
cat("\nFitting LMM to simulated data.\n")
model <- lmer(form, data = ndata, REML = reml,
              control = lmerControl(optimizer = "nlminbwrap",
                                    check.conv.grad = .makeCC("warning", tol = 4e-3, relTol = NULL),
                                    calc.derivs = FALSE))
# Print model summary
summary(model)

# Run model diagnostics 
cat("\nPlotting model diagnostics.\n")
md <- modelDiagnostics(model, ev.perc = .001)
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

residuals_diagnostics_plot_saved <- rp_combined

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

random_variable_diagnostics_plot_saved <- ran_combined

# Decompose the modal for inspection
cat("\nDecomposing model for inspection.\n")
inspect_fixed_matrix <- as.data.frame(as.matrix(model.matrix(model, type = "fixed")))
inspect_fixed_matrix <- getME(model, "X")
inspect_random_matrix <- model.matrix(model, type = "randomListRaw")
inspect_random_matrix <- getME(model, "Z")
inspect_random_effects <- ranef(model)
inspect_fixed_effects <- fixef(model)
inspect_residuals <- resid(model)
inspect_predictions <- predict(model, newdata = ndata) + inspect_residuals

# Make sure these components come together as expected
cat("\nChecking to make sure model predictions from package methods match expected predictions from model components.\n")
my_prediction <- function(r,print_values=FALSE) {
  
  Participant <- as.integer(ndata$Participant[r])
  
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
  
  u_rI <- as.matrix(inspect_random_matrix)[r,(2*Participant-1)]
  u_rST <- as.matrix(inspect_random_matrix)[r,(2*Participant)]
  
  rI <- inspect_random_effects$Participant[Participant,1]
  rBK <- inspect_random_effects$Participant[Participant,2]
  
  if (print_values) {
    cat("Values:\n")
    cat("Participant: ", Participant, "\n")
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
check_all()

##################################################################################
##################################################################################
# Part 3: Plot model results

# Extract fitted values from the model
cat("\nExtracting fitted values from the model.\n")
ndata$predicted <- fitted(model)

# Plot raw data as box plots
cat("\nPlotting raw data as box plots.\n")
color_mapping <- c( "online" = "blue", "terminal" = "red")
prelim_plot_reaches <- ggplot(ndata, aes(x=Block,y=.data[["NormalizedSpatialTargetError"]],color=LearningFeedback)) +
  geom_boxplot() +
  scale_color_manual(values = color_mapping) +
  labs(
    title = "Reach NSTE (Simulated Data)",
    x = "Block",
    y = "Normalized Spatial Target Error",
    color = "LFB"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
print(prelim_plot_reaches)

# Plot raw and fitted data as scatter plots
cat("\nPlotting raw and predicted data as scatter plots.\n")
color_mapping <- c( "online" = "blue", "terminal" = "red", "predicted" = "black" )
Predicted_values_plot_reaches_saved <- ggplot(data = ndata, aes(x = LearningFeedback, color = LearningFeedback)) +
  geom_point(aes(y=NormalizedSpatialTargetError), position = position_jitter(width = 0.3), size = 3, alpha = 0.2) +
  geom_point(aes(y=predicted, color = "predicted"), position = position_jitter(width = 0.15), size = 3, alpha = 0.2) +
  facet_wrap(~Block) +
  scale_color_manual(values = color_mapping) + 
  labs(title = "Actual and Predicted Values (Simu.)",
       x = "Block",
       y = "Normalized Spatial Target Error") +
  theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
  labs(color = "LFB")
print(Predicted_values_plot_reaches_saved)

# Find actual and predicted mean values per participant
cat("\nFinding actual and predicted mean values per participant.\n")
mean_data <- aggregate(NormalizedSpatialTargetError ~ Participant + Block + LearningFeedback,
                       data = ndata,
                       FUN = mean,
                       na.rm = TRUE)
mean_data_predicted_block <- aggregate(predicted ~ Participant + Block + LearningFeedback,
                                       data =  ndata,
                                       FUN = mean,
                                       na.rm = TRUE)

# Plot actual and predicted mean values per participant
cat("\nPlotting actual and predicted mean values per participant.\n")
color_mapping_LFB <- c( "online" = "blue", "terminal" = "red")
mean_plot <- ggplot(mean_data, aes(x = Block, y = NormalizedSpatialTargetError, group = Participant, color = LearningFeedback)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = color_mapping_LFB) + 
  labs(title = "Actual",
       x = "Block",
       y = "Mean Normalized Spatial Target Error") +
  theme_minimal() +
  theme(legend.position='none') 
mean_plot_p <- ggplot(mean_data_predicted_block, aes(x = Block, y = predicted, group = Participant, color = LearningFeedback)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = color_mapping_LFB) + 
  labs(title = "Predicted",
       x = "Block",
       y = "Mean Normalized Spatial Target Error") +
  theme_minimal() +
  theme(legend.position='none') 
combined_mean_plot <- mean_plot +
  theme(legend.position = "bottom", legend.title = element_text(size = 8), legend.text = element_text(size = 8)) +
  labs(color = "LFB") + 
  mean_plot_p +
  plot_layout(ncol = 2) +
  plot_annotation(title = "Participant Means by Block (Simulated Data)",
                  theme = theme_minimal(base_size = 15)) 
Predicted_values_plot_pmeans_saved <- combined_mean_plot
print(combined_mean_plot)

##################################################################################
##################################################################################
# Part 4: Perform bootstrapping
cat("\nPerforming bootstrapping.\n")

# Bootstrapping parameters
resamples_ <- 1000
cat("\nBootstrapping parameters:\n")
cat("Number of resamples: ", resamples_, "\n")

# Define stat function for bootstrapping
cat("\nDefining stat function for bootstrapping.\n")
fit <- function(data,indices) {
  bootstrap_data <- data[indices, ]
  model <- lmer(form, data = bootstrap_data, REML = reml,
                control = lmerControl(optimizer = "nlminbwrap",
                                      check.conv.grad = .makeCC("warning", tol = 4e-3, relTol = NULL),
                                      calc.derivs = FALSE))
  return(fixef(model))
}

# Perform bootstrapping (with stratified data)
cat("\nPerforming bootstrapping with stratafied data.\n")
bootstrap_results <- boot(data = ndata,
                          statistic = fit,
                          R = resamples_,
                          sim = "ordinary",
                          strata = ndata$stratum)

# Extract statistics from bootstrapping
mean_bs_effects1 <- c()
p_estimated1 <- c()
for (i in 1:4) {
  bs_r_mean <- mean(bootstrap_results$t[,i], na.rm = TRUE)
  mean_bs_effects1 <- c( mean_bs_effects1, bs_r_mean )
  p_estimated1 <- c( p_estimated1, 1.0 - mean( abs(bs_r_mean) >= abs(bootstrap_results$t[,i]-bs_r_mean) ) )
}

# Perform bootstrapping (without stratified data)
cat("\nPerforming bootstrapping without stratafied data.\n")
bootstrap_results <- boot(data = ndata,
                          statistic = fit,
                          R = resamples_,
                          sim = "ordinary")

# Extract statistics from bootstrapping
mean_bs_effects2 <- c()
p_estimated2 <- c()
for (i in 1:4) {
  bs_r_mean <- mean(bootstrap_results$t[,i], na.rm = TRUE)
  mean_bs_effects2 <- c( mean_bs_effects2, bs_r_mean )
  p_estimated2 <- c( p_estimated2, 1.0 - mean( abs(bs_r_mean) >= abs(bootstrap_results$t[,i]-bs_r_mean) ) )
}

# Print results
cat("Results with stratified bootstrapping:\n")
cat(mean_bs_effects1, "\n")
cat(p_estimated1, "\n")
cat("Results without stratified bootstrapping:\n")
cat(mean_bs_effects2, "\n")
cat(p_estimated2, "\n")

# Save
saved_fee <- round( mean_bs_effects1, digits = 3 )
saved_p <- round( p_estimated1, digits = 3 )
stipulated_fe <- round( c( F1RF2R, F1TeF2R, F2TeF1R, IntX ), digits = 3 )
sim_results <- rbind( stipulated_fe, saved_fee, saved_p )
rownames(sim_results) <- c( "Stipulated Fixed Effects", "Bootstrapped Fixed Effects", "p-values" )
colnames(sim_results) <- c("Intercept", "Sw(Tm)", "On(Sat)", "Interaction")

print(sim_results)