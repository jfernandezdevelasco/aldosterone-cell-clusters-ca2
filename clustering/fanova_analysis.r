library(fdANOVA)
library(fda)

##################### Helper functions

prepare_data <- function(labels_df, data_df, rows_to_keep, columns_to_delete) {
  columns_to_remove <- c("X", "X.1", "X.2")
  labels_df <- labels_df[, !(names(labels_df) %in% columns_to_remove)]
  
  data_df <- data_df[1:rows_to_keep, -1]
  data_df <- data_df[, !(names(data_df) %in% columns_to_delete)]
  
  data_df <- apply(data_df, c(1, 2), function(x) as.numeric(sub(",", ".", x)))
  data_df <- t(t(data_df) / data_df[1, ])
  
  label_col <- labels_df[, 2]
  duplicate_mask <- duplicated(label_col) | duplicated(label_col, fromLast = TRUE)
  filtered_labels <- label_col[duplicate_mask]
  
  return(list(data = data_df, labels = filtered_labels))
}

################### Configuration & Data Loading

base_dir <- "."
labels_dir <- paste0(base_dir, ".")
data_dir <- paste0(base_dir, ".")

load_patient_data <- function(file_name, rows_to_keep, cols_to_delete) {
  labels_df <- read.csv(paste0(labels_dir, file_name), sep = ';')
  data_df <- read.csv(paste0(data_dir, file_name), sep = ';')
  prepare_data(labels_df, data_df, rows_to_keep, cols_to_delete)
}

# Load all datasets into a named list
datasets <- list(
  SN  = load_patient_data(".", 1, c('X')),
  APA = load_patient_data(".", 1, c('X'))
)

# Extract matrices and labels
x <- lapply(datasets, function(ds) as.matrix(ds$data))
label_real <- lapply(datasets, `[[`, "labels")
label_rand <- lapply(label_real, sample)

###################### Plotting

plot_fanova_variants <- function(data_matrix, group_labels) {
  labels_char <- as.character(group_labels)
  plotFANOVA(x = data_matrix)
  plotFANOVA(x = data_matrix, group.label = labels_char)
  plotFANOVA(x = data_matrix, group.label = labels_char, separately = TRUE)
  plotFANOVA(x = data_matrix, group.label = labels_char, means = TRUE)
}

# Plot variants for SN (Real vs Random)
plot_fanova_variants(x$SN, label_real$SN)
plot_fanova_variants(x$SN, label_rand$SN)

##################### Statistical Testing

set.seed(123)
library(devtools)

source_url("https://raw.githubusercontent.com/jfernandezdevelasco/Fb-Heteroscedastic-Test/main/F_b_hetero.R")

# Wrapper to run the bootstrap tests across the list
run_tests <- function(data_list, label_list) {
  lapply(seq_along(data_list), function(i) {
    hetero.F.type.boot(x = data_list[[i]], group.label = label_list[[i]], nrFb = 1000)
  })
}

# Execute tests for real and random labels
results_real <- run_tests(x, label_real)
results_rand <- run_tests(x, label_rand)

#####################.  Export Results

# Aggregate F statistics and P values
stat_F_b <- c(sapply(results_real, `[[`, "statF"), sapply(results_rand, `[[`, "statF"))
p_val_F_b <- c(sapply(results_real, `[[`, "pvalueFb"), sapply(results_rand, `[[`, "pvalueFb"))

sample_names <- c('SN', 'APA','SN_r_l', 'APA_r_l')

results_table <- data.frame(
  row.names = sample_names,
  Test_Statistic = stat_F_b,
  P_Value = p_val_F_b
)

output_path <- "."
write.csv(results_table, output_path, row.names = TRUE)

print(results_table)