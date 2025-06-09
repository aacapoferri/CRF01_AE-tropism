library(tidyverse)
library(caret)
library(pROC)

# Load predictions
df <- read_csv("FPR2.csv", show_col_types = FALSE)

# Filter to X4/R5 only and ensure factors are set
df <- df %>%
  filter(true_label %in% c("R5", "X4"), model_call %in% c("R5", "X4")) %>%
  mutate(
    true_label = factor(true_label, levels = c("R5", "X4")),
    model_call = factor(model_call, levels = c("R5", "X4"))
  )

# ---- Function to calculate metrics per class ----
get_metrics <- function(pred_class, truth, positive_class) {
  cm <- confusionMatrix(pred_class, truth, positive = positive_class)
  tibble(
    Class = positive_class,
    Accuracy = cm$overall["Accuracy"],
    Sensitivity_TPR = cm$byClass["Sensitivity"],
    Specificity_TNR = cm$byClass["Specificity"],
    Precision_PPV = cm$byClass["Precision"],
    FPR = 1 - cm$byClass["Specificity"],
    FNR = 1 - cm$byClass["Sensitivity"],
    F1_Score = cm$byClass["F1"]
  )
}

# ---- Metrics for X4 ----
metrics_x4 <- get_metrics(df$model_call, df$true_label, "X4")

# ---- Metrics for R5 (invert labels) ----
df_flipped <- df %>%
  mutate(
    model_call = fct_rev(model_call),
    true_label = fct_rev(true_label)
  )
metrics_r5 <- get_metrics(df_flipped$model_call, df_flipped$true_label, "R5")

# ---- Combine and export metrics ----
all_metrics <- bind_rows(metrics_x4, metrics_r5)
write_csv(all_metrics, "tropism_classification_metrics.csv")

# ---- ROC curves (base R plotting) ----
roc_x4 <- roc(response = df$true_label, predictor = df$.pred_X4, levels = c("R5", "X4"))
roc_r5 <- roc(response = df$true_label, predictor = df$.pred_R5, levels = c("X4", "R5"))  # flipped

plot(roc_x4, col = "red", lwd = 2, main = "ROC Curves for Tropism Prediction", legacy.axes = TRUE)
lines(roc_r5, col = "blue", lwd = 2)
abline(a = 0, b = 1, lty = 3, col = "gray")
legend("bottomright", legend = c(
  paste0("X4 AUC = ", round(auc(roc_x4), 3)),
  paste0("R5 AUC = ", round(auc(roc_r5), 3))
), col = c("red", "blue"), lwd = 2, lty = c(1, 2))
