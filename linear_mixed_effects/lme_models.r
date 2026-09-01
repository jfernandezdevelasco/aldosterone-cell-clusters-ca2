library(lme4)
library(lmerTest)
library(nlme)
library(broom.mixed)
library(ggplot2)
library(effects)
library(reticulate)

# Load clustering results based on distance threshold

load_cluster_data <- function(filter_method, distance_threshold, data_dir = ".") {
  
  pickle <- import("pickle")
  builtins <- import_builtins()
  pkl_filename <- file.path(
    data_dir, 
    sprintf("clustering_filter_%s_distance_%s.pkl", filter_method, as.character(distance_threshold))
  )
  py_file <- builtins$open(pkl_filename, "rb")
  results <- pickle$load(py_file)
  py_file$close()
  dat1_list <- list()
  dat2_list <- list()
  
  # Iterate through nested dictionary: results[patient_id][sample_name]
  valid_ids <- setdiff(names(results), "config")
  for (patient_id in valid_ids) {
    patient_dict <- results[[patient_id]]
    
    for (sample_name in names(patient_dict)) {
      sample_dict <- patient_dict[[sample_name]]
      
      if ("summary" %in% names(sample_dict)) {
        summary_df <- as.data.frame(sample_dict[["summary"]])
        summary_df$patient_id <- patient_id
        summary_df$sample_name <- sample_name
        dat1_list[[length(dat1_list) + 1]] <- summary_df
      }
      
      # Extract 'number_of_cells_in_predicted_clusters' for dat2
      if ("number_of_cells_in_predicted_clusters" %in% names(sample_dict)) {
        cells_data <- sample_dict[["number_of_cells_in_predicted_clusters"]]
        
        # Convert list/vector to data frame if needed
        cells_df <- as.data.frame(cells_data)
        cells_df$patient_id <- patient_id
        cells_df$sample_name <- sample_name
        dat2_list[[length(dat2_list) + 1]] <- cells_df
      }
    }
  }
  dat1 <- do.call(rbind, dat1_list)
  dat2 <- do.call(rbind, dat2_list)
  
  return(list(matrix1 = dat1, matrix2 = dat2))
}

# Run diagnostic plots for lme model assumptions check

check_model_assumptions <- function(model, data, response_var, label = "Response") {
  
  par(mfrow = c(2, 3))  
  # Random effects distribution check
  re_patient <- ranef(model)$patient_id[[1]]
  hist(
    re_patient, 
    main = paste("Random Effects:", label), 
    xlab = "Random Intercepts"
  )
  qqnorm(re_patient, main = "Q-Q Plot: Random Effects")
  qqline(re_patient, col = "red")
  
  # Residuals distribution check
  resids <- residuals(model)
  hist(
    resids, 
    main = paste("Residuals:", label), 
    xlab = "Residuals"
  )
  qqnorm(resids, main = "Q-Q Plot: Residuals")
  qqline(resids, col = "red")
  
  # Residual autocorrelation & fitted values
  plot(
    resids, 
    type = "p", 
    main = "Residual Plot", 
    xlab = "Index", 
    ylab = "Residuals", 
    pch = 19, 
    col = "blue"
  )
  abline(h = 0, col = "red", lty = 2)
  
  plot(
    response_var, 
    resids, 
    main = "Residuals vs Response", 
    xlab = label, 
    ylab = "Residuals", 
    pch = 20, 
    col = "darkgreen"
  )
  abline(h = 0, col = "red", lty = 2)
  
  par(mfrow = c(1, 1))
}

######################### Data ingestion 
# distance <- 0.05,0.10,0.15,0.2,0.25,0.3
# filter_method <- savgol_filter,mean_filter
# data_dir <- # directory 
data_list <- load_cluster_data(filter_method,distance,data_dir)
dat1 <- data_list$matrix1
dat2 <- data_list$matrix2

dat1$sn_apa <- as.factor(dat1$sn_apa)
dat2$sn_apa <- as.factor(dat2$sn_apa)

######################### Model fitting and summaries

# Model 1: Predict subclusters in functional region
mix1 <- lmer(
  number_pred_sub_clust_in_fun ~ sn_apa + num_cells_in_fun + (1 | patient_id), 
  data = dat1
)
summary(mix1)
fe_table <- tidy(mix1, effects = "fixed")
re_table <- glance(mix1, effects = "ran_coefs")
anova(mix1)

# Model 2: Number of cells in functional region ~ sn_apa
mix2 <- lmer(
  num_cells_in_fun ~ sn_apa + (1 | patient_id), 
  data = dat1
)
summary(mix2)
anova(mix2)

# Model 3: Number of cells in predicted subclusters ~ size + sn_apa
mix3 <- lmer(
  num_cells_pred_sub_in_fun ~ num_cells_in_fun + sn_apa + (1 | patient_id), 
  data = dat2
)
summary(mix3)
anova(mix3)

####################################### Visualization

ggplot(dat1, aes(x = num_cells_in_fun, y = number_pred_sub_clust_in_fun, color = sn_apa)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(
    x = "Number of Cells in Functional Cluster",
    y = "Predicted Subclusters",
    color = "SN/APA"
  )

re <- ranef(mix3)$patient_id
fe <- fixef(mix1)
clr <- rainbow(nrow(re))

par(mfrow = c(1, 2))

plot(
  number_pred_sub_clust_in_fun ~ sn_apa + num_cells_in_fun, 
  data = dat1, 
  col = clr[as.numeric(patient_id)], 
  main = "Pred w/ points"
)
lapply(seq_len(nrow(re)), function(x) {
  abline(fe[1] + re[x, 1], fe[2] + re[x, 2], col = clr[x])
})

plot(
  Reaction ~ Days, 
  data = sleepstudy, 
  col = clr[as.numeric(Subject)], 
  main = "Pred w/o points", 
  type = "n"
)
lapply(seq_len(nrow(re)), function(x) {
  abline(fe[1] + re[x, 1], fe[2] + re[x, 2], col = clr[x])
})

par(mfrow = c(1, 1))

# Effects plots
model_effects <- allEffects(mix1)
plot(model_effects)

######################### Model diagnostics / assumption verification

# Diagnostics for Model 1
check_model_assumptions(
  model = mix1, 
  data = dat1, 
  response_var = dat1$number_pred_sub_clust_in_fun, 
  label = "Model 1: Subclusters"
)

# Diagnostics for Model 2
check_model_assumptions(
  model = mix2, 
  data = dat1, 
  response_var = dat1$num_cells_in_fun, 
  label = "Model 2: Cells in Fun"
)

# Diagnostics for Model 3
check_model_assumptions(
  model = mix3, 
  data = dat2, 
  response_var = dat2$num_cells_pred_sub_in_fun, 
  label = "Model 3: Pred Sub Cells"
)