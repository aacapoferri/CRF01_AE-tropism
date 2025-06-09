predict_from_fasta <- function(fasta_path, model_path, output_csv) {
  library(Biostrings)
  library(tidyverse)
  library(tidymodels)
  
  # Load trained model and recipe
  load(model_path)  # loads: final_fit, rec_prep
  
  # Feature extraction function
  extract_features <- function(sequences, reference_seq) {
    ref_chars <- unlist(strsplit(reference_seq, ""))
    features <- map2_dfr(names(sequences), as.character(sequences), function(id, seq) {
      query_chars <- unlist(strsplit(seq, ""))
      aligned_ref <- ref_chars
      aligned_query <- query_chars
      len <- min(length(aligned_ref), length(aligned_query))
      ref_to_query_pos <- which(aligned_ref != "-" & aligned_query != "-")
      v3_query_pos <- aligned_query[ref_to_query_pos]
      
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
      
      pos_charge <- sum(v3_query_pos %in% c("K", "R", "H"))
      neg_charge <- sum(v3_query_pos %in% c("D", "E"))
      net_charge <- pos_charge - neg_charge
      
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
      
      pos25_aa <- query_chars[ref_to_query_pos[25]]
      
      tibble(
        seq_name = id,
        glyco_status = glyco_status,
        net_charge = net_charge,
        features_x4 = features_x4,
        pos25 = pos25_aa
      )
    })
    return(features)
  }
  
  # Step 1: Load sequences
  seqs <- readAAStringSet(fasta_path)
  
  # Step 2: Extract features
  reference_seq <- "CTRPSNNTRTSITIGPGQVFYRTGDIIGDIRKAYC"
  feature_data <- extract_features(seqs, reference_seq)
  
  # ✅ FIX: Use var_info to get required pre-baking variables
  required_vars <- rec_prep$var_info %>%
    filter(role %in% c("predictor", "id")) %>%
    pull(variable)
  
  # Subset features
  feature_data <- feature_data %>% select(all_of(required_vars))
  
  # Step 3: Predict
  pred_class <- predict(final_fit, new_data = feature_data, type = "class")
  pred_probs <- predict(final_fit, new_data = feature_data, type = "prob")
  
  # Step 4: Combine and assign confidence
  results <- bind_cols(feature_data, pred_class, pred_probs) %>%
    rowwise() %>%
    mutate(
      max_prob = max(c_across(where(is.numeric) & starts_with(".pred_")), na.rm = TRUE),
      confidence = case_when(
        max_prob >= 0.85 ~ "High",
        max_prob >= 0.65 ~ "Moderate",
        TRUE ~ "Low"
      )
    ) %>%
    ungroup()
  
  # Step 5: Export
  write_csv(results, output_csv)
  message("✅ Prediction complete: ", output_csv)
}
