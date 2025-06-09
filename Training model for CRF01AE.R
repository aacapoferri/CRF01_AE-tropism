### Training tropism for CRF01AE###

# Load libraries
library(tidyverse)
library(tidymodels)
library(xgboost)
library(textrecipes)
library(themis)
library(Biostrings)

# Step 1: Feature extraction function
extract_features <- function(sequences, reference_seq) {
  ref_chars <- unlist(strsplit(reference_seq, ""))
  
  features <- map_dfr(sequences, function(seq) {
    query_chars <- unlist(strsplit(seq, ""))
    aligned_ref <- ref_chars
    aligned_query <- query_chars
    len <- min(length(aligned_ref), length(aligned_query))
    
    # Positional mapping (non-gap positions)
    ref_to_query_pos <- which(aligned_ref != "-" & aligned_query != "-")
    v3_query_pos <- aligned_query[ref_to_query_pos]
    
    # Glycosylation motif detection (NNT pattern)
    glyco_status <- "Absent"
    ref_glyco_pos <- which(sapply(1:(length(ref_chars) - 2), function(i) {
      paste0(ref_chars[i:(i + 2)], collapse = "") == "NNT"
    }))
    if (length(ref_glyco_pos) > 0 && (ref_glyco_pos + 2) <= length(query_chars)) {
      glyco_motif <- paste0(query_chars[ref_glyco_pos:(ref_glyco_pos + 2)], collapse = "")
      if (substr(glyco_motif, 1, 2) == "NN") {
        third <- substr(glyco_motif, 3, 3)
        glyco_status <- ifelse(third %in% c("I", "M"), "Lost", "Present")
      }
    }
    
    # Net charge
    pos_charge <- sum(v3_query_pos %in% c("K", "R", "H"))
    neg_charge <- sum(v3_query_pos %in% c("D", "E"))
    net_charge <- pos_charge - neg_charge
    
    # X4-like features
    features_x4 <- FALSE
    check_pos <- function(pos, accepted) {
      aa <- query_chars[ref_to_query_pos[pos]]
      is.na(aa) || aa %in% accepted
    }
    if (
      check_pos(8, c("V", "I", "-")) ||
      check_pos(10, c("I")) ||
      check_pos(11, c("R", "K")) ||
      check_pos(12, c("F")) ||
      check_pos(24, c("E")) ||
      any(sapply(32:35, function(p) is.na(query_chars[ref_to_query_pos[p]]) || query_chars[ref_to_query_pos[p]] == "-"))
    ) {
      features_x4 <- TRUE
    }
    
    # Position 25 logic
    pos25_aa <- query_chars[ref_to_query_pos[25]]
    
    tibble(
      glyco_status = glyco_status,
      net_charge = net_charge,
      features_x4 = features_x4,
      pos25 = pos25_aa
    )
  })
  
  return(features)
}

# Step 2: Load labeled training data
training_data <- read_tsv("training_new.tsv", show_col_types = FALSE)

# Step 3: Extract features from sequences using reference
reference_seq <- "CTRPSNNTRTSITIGPGQVFYRTGDIIGDIRKAYC"
feature_data <- extract_features(training_data$sequence, reference_seq)
training_set <- bind_cols(training_data %>% select(seq_name, tropism), feature_data)

# Step 4: Define preprocessing recipe
rec <- recipe(tropism ~ ., data = training_set) %>%
  update_role(seq_name, new_role = "id") %>%
  step_mutate(features_x4 = as.factor(features_x4)) %>%
  step_impute_mode(all_nominal_predictors()) %>%
  step_impute_mean(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_smote(tropism)

# Step 5: Define XGBoost model
xgb_model <- boost_tree(trees = 500, mtry = 3, learn_rate = 0.1) %>%
  set_mode("classification") %>%
  set_engine("xgboost")

# Step 6: Create workflow and fit
wf <- workflow() %>%
  add_model(xgb_model) %>%
  add_recipe(rec)

set.seed(123)
xgb_fit <- fit(wf, data = training_set)

#Prepare final object names expected by prediction script
final_fit <- xgb_fit
rec_prep <- workflows::extract_recipe(final_fit)  # safely extract prepared recipe

# Save trained model and recipe
# Extract prepared recipe for saving
rec <- workflows::extract_recipe(xgb_fit)

# Save model and recipe (must be named xgb_fit and rec for predict_from_fasta)
save(xgb_fit, rec, file = "CRF01tropism_xgb_fit.RData")

# Step 7: Predict class labels and probabilities
# Step 7: Predict tropism and assign confidence levels using the trained model

library(tidyverse)

# Predict class labels
pred_class <- predict(xgb_fit, new_data = training_set, type = "class")

# Rename the class column to "predicted"
colnames(pred_class) <- "predicted"

# Predict class probabilities
pred_probs <- predict(xgb_fit, new_data = training_set, type = "prob")

# DEBUG: Check structure of predictions (optional)
# print("Class predictions:")
# print(glimpse(pred_class))
# print("Probability predictions:")
# print(glimpse(pred_probs))

# Combine predictions with training data and assign confidence levels
results <- bind_cols(training_set, pred_class, pred_probs) %>%
  rowwise() %>%
  mutate(
    max_prob = max(c_across(starts_with(".pred_")), na.rm = TRUE),
    confidence = case_when(
      max_prob >= 0.85 ~ "High",
      max_prob >= 0.65 ~ "Moderate",
      TRUE ~ "Low"
    )
  ) %>%
  ungroup()

# View final predictions and confidence levels
results %>%
  select(seq_name, tropism, predicted, starts_with(".pred_"), confidence) %>%
  print(n = 20)

# Export results to CSV
write_csv(
  results %>%
    select(seq_name, tropism, predicted, starts_with(".pred_"), confidence),
  "CRF01tropism_xgb_fit.csv"
)
