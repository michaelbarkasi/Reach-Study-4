# rTable outlier check: Remove all reaches/trials which are outliers
cat("\nFinding outlier reaches in the rTable and removing them.\n")

# Helper functions
cat("\nLoading helper functions.\n")
find_outlier_rows <- function(pr_array) {
  
  O_P <- pr_array[,1]
  O_R <- pr_array[,2]
  
  outlier__rows <- c()
  for ( n in 1:length(O_P) ) {
    outlier__rows <- c(outlier__rows, 
                       which(rTable$Participants==O_P[n] & rTable$ReachNum==O_R[n])
    )
  }
  return(outlier__rows)
  
}
find_outliers__ <- function(dep_v, cf) {
  
  participant_c <- c()
  reach_c <- c()
  
  for( p in pTable$number ) {
    
    out_rows <- which(is_outlier(dep_v[which(rTable$Participants==p)], coef = cf))
    r_cut <- rTable$ReachNum[which(rTable$Participants==p)]
    
    participant_c <- c( participant_c, rep(p, length(out_rows)))
    reach_c <- c( reach_c, r_cut[out_rows] )
    
  }
  
  # Santiy Check
  if ( length(participant_c) != length(reach_c) ) {
    cat("\nWarning! Something went wrong finding outliers; unequal length to participant and reach columns.\n")
  }
  
  return(cbind(participant_c,reach_c)) # First column will be participant, second will be reach number
  
}

# Find starting system difference in jerk, target error, and path error
target_diff <- function() {
  running_diff <- mean(abs(rTable$RotationTargetErrorNorm - rTable$SpatialTargetErrorNorm),na.rm=TRUE)*100
  return(running_diff) 
}
path_diff <- function() {
  running_diff <- mean(abs(rTable$RotationPathErrorNorm - rTable$SpatialPathErrorNorm),na.rm=TRUE)*100
  return(running_diff) 
}
original_jerk_diff <- jerk_diff()
original_target_diff <- target_diff()
original_path_diff <- path_diff()
cat("\nMean normalized jerk measurement difference (percentage points) between systems before outlier removal: ", original_jerk_diff, "\n")
cat("Mean normalized target error measurement difference (percentage points) between systems before outlier removal: ", original_target_diff, "\n")
cat("Mean normalized path error measurement difference (percentage points) between systems before outlier removal: ", original_path_diff, "\n")

# Remove outliers
remove_outliers <- function(cf_) {
  
  # Find outliers
  cat("\nFirst identify outlier reaches for each participant.\n")
  outlier_SUj <- find_outliers__(rTable$SUjerkNorm,cf_)
  outlier_OTj <- find_outliers__(rTable$OTjerkNorm,cf_)
  outlier_SPE <- find_outliers__(rTable$SpatialPathErrorNorm,cf_)
  outlier_STE <- find_outliers__(rTable$SpatialTargetErrorNorm,cf_)
  outlier_RPE <- find_outliers__(rTable$RotationPathErrorNorm,cf_)
  outlier_RTE <- find_outliers__(rTable$RotationTargetErrorNorm,cf_)
  
  # Find rTable row numbers
  cat("\nNext identify the rows of those outliers in the rTable.\n")
  outlier_rows_SUj <- find_outlier_rows(outlier_SUj)
  outlier_rows_OTj <- find_outlier_rows(outlier_OTj)
  outlier_rows_SPE <- find_outlier_rows(outlier_SPE)
  outlier_rows_STE <- find_outlier_rows(outlier_STE)
  outlier_rows_RPE <- find_outlier_rows(outlier_RPE)
  outlier_rows_RTE <- find_outlier_rows(outlier_RTE)
  
  # Make row vector
  outlier_rows_total <- c(outlier_rows_SUj,
                          outlier_rows_OTj,
                          outlier_rows_SPE,
                          outlier_rows_STE,
                          outlier_rows_RPE,
                          outlier_rows_RTE)
  
  outlier_rows_total <- unique(outlier_rows_total)
  
  num_of_outliers_found <<- length(outlier_rows_total)
  num_of_reaches_total <<- length(rTable$Participants)
  outliers_as_percent_of_total <<- (length(outlier_rows_total)/length(rTable$Participants))*100
  num_of_reaches_after_shorts_and_floats_removed <<- length(rTable$Participants[which(!is.na(rTable$SpatialPathError))])
  outliers_as_percent_of_total_after_shorts_and_floats_removed <<- 
    (length(outlier_rows_total)/length(rTable$Participants[which(!is.na(rTable$SpatialPathError))]))*100
  
  cat("\nNumber of outliers found:", length(outlier_rows_total), "\n")
  cat("Total number of reaches:", length(rTable$Participants), "\n")
  cat("Outliers as percent of total:", (length(outlier_rows_total)/length(rTable$Participants))*100, "\n")
  cat("Total number of reaches after shorts and floats removed:", length(rTable$Participants[which(!is.na(rTable$SpatialPathError))]), "\n")
  cat("Outliers as percent of total after shorts and floats removed:", 
      (length(outlier_rows_total)/length(rTable$Participants[which(!is.na(rTable$SpatialPathError))]))*100, "\n")
  
  # Save rTable before removing rows
  cat("\nSaving original rTable before removing outlier rows (rTable_saved).\n")
  rTable_saved <<- rTable
  
  # Fill outlier rows with NA in the relevant columns
  cat("\nFilling outlier rows with NA in the relevant columns.\n")
  rTable$SUjerkNorm[outlier_rows_total] <<- NA
  rTable$OTjerkNorm[outlier_rows_total] <<- NA
  rTable$SpatialPathErrorNorm[outlier_rows_total] <<- NA
  rTable$SpatialTargetErrorNorm[outlier_rows_total] <<- NA
  rTable$RotationPathErrorNorm[outlier_rows_total] <<- NA
  rTable$RotationTargetErrorNorm[outlier_rows_total] <<- NA
  
}
remove_outliers(cf_ = outlierIQR_cutoff)

# Count remaining reaches per participant
cat("\nCounting remaining reaches per participant.\n")
count_remaining_reaches_per_participant <- function() {
  
  remaining_reaches <- c()
  for ( p in pTable$number ) {
    remaining_reaches <- c(remaining_reaches, sum(!is.na(rTable$SUjerkNorm[which(rTable$Participants==p)])))
  }
  
  mean_remaining_reaches_per_participant <<- mean(remaining_reaches)
  cat("Mean remaining reaches per participant after outlier removal:", mean(remaining_reaches),"\n")
  
  return(remaining_reaches)
  
}
re_reaches <- count_remaining_reaches_per_participant()

# Plot remaining reaches per participants
cat("\nScatter plot and histograph of remaining reaches per participant.\n")
fig_num <- fig_num + 1
plot_title1 <- "Reaches after Outlier Removal"
plot(re_reaches,
     main=paste("Fig ", sec_num, ".", fig_num, ": ", plot_title1, sep = ""),
     ylab="Remaining Reaches",
     xlab="Participant Number")
fig_num <- fig_num + 1
plot_title2 <- "Distribution, Remaining Reaches"
hist(re_reaches,
     main=paste("Fig ", sec_num, ".", fig_num, ": ", plot_title2, sep = ""),
     ylab="Number of Participants",
     xlab="Remaining Reaches")

remaining_reaches_pp_plot <- ggplot(data.frame(x = 1:length(re_reaches), y = re_reaches)) +
  geom_point(aes(x = x, y = y), color = "blue", size = 3) +
  labs(title = plot_title1,
       y = "Remaining Reaches",
       x = "Participant Number") +
  theme_minimal() +
  theme(aspect.ratio = 1)
  
remaining_reaches_pp_hist <- ggplot(data.frame(x = re_reaches)) +
  geom_histogram(aes(x = x), fill = "green", color = "black", bins = 15) +
  labs(title = plot_title2,
       y = "Number of Participants",
       x = "Remaining Reaches") +
  theme_minimal() +
  theme(aspect.ratio = 1)

afterremoval_jerk_diff <- jerk_diff()
afterremoval_target_diff <- target_diff()
afterremoval_path_diff <- path_diff()
cat("\nMean normalized jerk measurement difference (percentage points) between systems after outlier removal: ", afterremoval_jerk_diff, "\n")
cat("Mean normalized target error measurement difference (percentage points) between systems after outlier removal: ", afterremoval_target_diff, "\n")
cat("Mean normalized path error measurement difference (percentage points) between systems after outlier removal: ", afterremoval_path_diff, "\n")

# rTable outlier check: Remove all reaches/trials which are outliers
cat("\nFinding reaches in the rTable with high system disagreement and removing them.\n")

remove_disagreement <- function() {
  
  dis_rows_J <- which(abs(rTable$SUjerkNorm-rTable$OTjerkNorm) > dis_limit)
  dis_rows_T <- which(abs(rTable$RotationTargetErrorNorm - rTable$SpatialTargetErrorNorm) > dis_limit)
  dis_rows_P <- which(abs(rTable$RotationPathErrorNorm - rTable$SpatialPathErrorNorm) > dis_limit)
  dis_rows <- c( dis_rows_J, dis_rows_T, dis_rows_P)
  dis_rows <- unique(dis_rows)
  
  num_of_ExD_found <<- length(dis_rows)
  num_of_reaches_total <<- length(rTable$Participants)
  ExD_as_percent_of_total <<- (length(dis_rows)/length(rTable$Participants))*100
  num_of_reaches_after_OSF_removed <<- length(rTable$Participants[which(!is.na(rTable$SpatialPathErrorNorm))])
  ExD_as_percent_of_total_after_OSF_removed <<- 
    (length(dis_rows)/length(rTable$Participants[which(!is.na(rTable$SpatialPathErrorNorm))]))*100
  
  cat("\nNumber of reaches found with high disagreement:", length(dis_rows),"\n")
  cat("Total number of reaches:", length(rTable$Participants), "\n")
  cat("High disagreement as percent of total:", (length(dis_rows)/length(rTable$Participants))*100, "\n")
  cat("Total number of reaches after outliers, shorts and floats removed:", length(rTable$Participants[which(!is.na(rTable$SpatialPathErrorNorm))]), "\n")
  cat("Outliers as percent of total after shorts and floats removed:", 
      (length(dis_rows)/length(rTable$Participants[which(!is.na(rTable$SpatialPathErrorNorm))]))*100, "\n")
  
  # Save rTable before removing rows
  cat("\nSaving original rTable before removing extreme-disagreeing rows (rTable_saved2).\n")
  rTable_saved2 <<- rTable
  
  # Fill outlier rows with NA in the relevant columns
  cat("\nFilling outlier rows with NA in the relevant columns.\n")
  rTable$SUjerkNorm[dis_rows] <<- NA
  rTable$OTjerkNorm[dis_rows] <<- NA
  rTable$SpatialPathErrorNorm[dis_rows] <<- NA
  rTable$SpatialTargetErrorNorm[dis_rows] <<- NA
  rTable$RotationPathErrorNorm[dis_rows] <<- NA
  rTable$RotationTargetErrorNorm[dis_rows] <<- NA
  
}
remove_disagreement()

fig_num <- fig_num + 1
plot_title3 <- "Norm. Jerk"
plot_title4 <- "Norm. Target Error"
plot_title5 <- "Norm. Path Error"

hist(abs(rTable$SUjerkNorm-rTable$OTjerkNorm),
     main=paste("Fig ", sec_num, ".", fig_num, ": ", plot_title3, sep = ""),
     ylab="Number of Reaches",
     xlab="System Difference")

system_disagreement_histJ <- ggplot(rTable, aes(x = abs(SUjerkNorm - OTjerkNorm))) +
  geom_histogram(fill = "orange", color = "black", bins = 15) +
  labs(title = plot_title3,
       y = "Number of Reaches",
       x = "System Difference") +
  theme_minimal()

system_disagreement_histT <- ggplot(rTable, aes(x = abs(RotationTargetErrorNorm - SpatialTargetErrorNorm))) +
  geom_histogram(fill = "orange3", color = "black", bins = 15) +
  labs(title = plot_title4,
       y = "Number of Reaches",
       x = "System Difference") +
  theme_minimal()

system_disagreement_histP <- ggplot(rTable, aes(x = abs(RotationPathErrorNorm - SpatialPathErrorNorm))) +
  geom_histogram(fill = "orange4", color = "black", bins = 15) +
  labs(title = plot_title5,
       y = "Number of Reaches",
       x = "System Difference") +
  theme_minimal()

final_jerk_diff <- jerk_diff()
final_target_diff <- target_diff()
final_path_diff <- path_diff()
cat("\nMean normalized jerk measurement difference (percentage points) between systems after extreme disagreement removal: ", final_jerk_diff, "\n")
cat("Mean normalized target error measurement difference (percentage points) between systems after extreme disagreement removal: ", final_target_diff, "\n")
cat("Mean normalized path error measurement difference (percentage points) between systems after extreme disagreement removal: ", final_path_diff, "\n")

count_remaining_reaches <- function() {
  
  SPEn <- length(which(!is.na(rTable$SpatialPathErrorNorm)))
  STEn <- length(which(!is.na(rTable$SpatialTargetErrorNorm)))
  RPEn <- length(which(!is.na(rTable$RotationPathErrorNorm)))
  RTEn <- length(which(!is.na(rTable$RotationTargetErrorNorm)))
  SJn <- length(which(!is.na(rTable$SUjerkNorm)))
  RJn <- length(which(!is.na(rTable$OTjerkNorm)))
  
  if (1!=length(unique(c(SPEn,STEn,RPEn,RTEn,SJn,RJn)))) {
    cat("\nWarning! unequal number of reaches remaining.\n")
    cat("SPEn:", SPEn, "\n")
    cat("STEn:", STEn, "\n")
    cat("RPEn:", RPEn, "\n")
    cat("RTEn:", RTEn, "\n")
    cat("SJn:", SJn, "\n")
    cat("RJn:", RJn, "\n")
  }
  
  final_remaining_reaches <<- SPEn
  final_remaining_reaches_as_percent_of_total <<- (SPEn/length(rTable$Participants))*100
  cat("\nNumber of reaches remaining:", SPEn, "\n")
  cat("Remaining as percent of total:", (SPEn/length(rTable$Participants))*100, "\n")
  
}
count_remaining_reaches()


