# ============================================
# Imputation & Evaluation 
# Author: Paulina Skurzak
# ============================================

#=========================
# 1. Load libraries and set working directory 
#=========================
library(dplyr)
library(mice)
library(VIM)
library(Metrics)
library(ggplot2)
library(readr)

setwd("~/Documents/GRI Research Fellow")

#=========================
# 2. Load cleaned dataset
#=========================
analysis_df <- read_csv("outputs/analysis_clean.csv")
names(analysis_df) <- make.unique(names(analysis_df))

#=========================
# 3. Keep only numeric columns of interest
#=========================
impute_df <- analysis_df %>%
  select(where(is.numeric)) %>%
  select(any_of(c("health_exp_pc", "gdp_pc", "unemployment",
                  "school_enroll", "female_labor", "oda_pc",
                  "agriculture_gdp")))

#=========================
# 4. Define evaluation function (real missingness)
#=========================
evaluate_imputation <- function(imputed_df, original_df, method_name) {
  common_cols <- intersect(names(imputed_df), names(original_df))
  rmse_vals <- c()
  mae_vals  <- c()
  
  for (col in common_cols) {
    obs <- original_df[[col]][!is.na(original_df[[col]])]
    imp <- imputed_df[[col]][!is.na(imputed_df[[col]])]
    
    if (length(obs) > 5 && length(imp) > 5) {
      rmse_vals <- c(rmse_vals, abs(mean(imp, na.rm = TRUE) - mean(obs, na.rm = TRUE)))
      mae_vals  <- c(mae_vals, abs(sd(imp, na.rm = TRUE) - sd(obs, na.rm = TRUE)))
    }
  }
  
  data.frame(
    Method = method_name,
    RMSE = mean(rmse_vals, na.rm = TRUE),
    MAE  = mean(mae_vals,  na.rm = TRUE)
  )
}

#=========================
# 5. Apply Imputation Methods (on scaled data)
#=========================

scaled_df <- as.data.frame(scale(impute_df))
colnames(scaled_df) <- names(impute_df)

# Mean imputation
mean_imp <- scaled_df
for (col in names(mean_imp)) {
  mean_imp[[col]][is.na(mean_imp[[col]])] <- mean(mean_imp[[col]], na.rm = TRUE)
}
result_mean <- evaluate_imputation(mean_imp, scaled_df, "Mean")

# KNN (on scaled data)
set.seed(42)
knn_imp <- VIM::kNN(scaled_df, k = 5, imp_var = FALSE)[, names(scaled_df)]
result_knn <- evaluate_imputation(knn_imp, scaled_df, "KNN")

# Regression imputation
set.seed(42)
reg_imp <- mice(scaled_df, method = "norm.predict", m = 1, maxit = 5, printFlag = FALSE)
reg_complete <- complete(reg_imp)
result_reg <- evaluate_imputation(reg_complete, scaled_df, "Regression")

# MICE (Predictive mean matching)
set.seed(42)
mice_imp <- mice(scaled_df, method = "pmm", m = 5, maxit = 5, printFlag = FALSE)
mice_complete <- complete(mice_imp)
result_mice <- evaluate_imputation(mice_complete, scaled_df, "MICE")

# Random Forest
set.seed(42)
rf_imp <- mice(scaled_df, method = "rf", m = 1, maxit = 5, printFlag = FALSE)
rf_complete <- complete(rf_imp)
result_rf <- evaluate_imputation(rf_complete, scaled_df, "Random Forest")

#=========================
# 6. Combine and display all results
#=========================
results_summary <- bind_rows(result_mean, result_knn, result_reg, result_mice, result_rf)
print(results_summary)

#=========================
# 7. Baseline and Post-Imputation Models (on scaled data)
#=========================
complete_cases <- na.omit(scaled_df)

baseline_model <- lm(health_exp_pc ~ gdp_pc + unemployment + school_enroll +
                       female_labor + oda_pc + agriculture_gdp,
                     data = complete_cases)

model_mean <- lm(health_exp_pc ~ gdp_pc + unemployment + school_enroll +
                   female_labor + oda_pc + agriculture_gdp, data = mean_imp)
model_knn <- lm(health_exp_pc ~ gdp_pc + unemployment + school_enroll +
                  female_labor + oda_pc + agriculture_gdp, data = knn_imp)
model_reg <- lm(health_exp_pc ~ gdp_pc + unemployment + school_enroll +
                  female_labor + oda_pc + agriculture_gdp, data = reg_complete)
model_mice <- lm(health_exp_pc ~ gdp_pc + unemployment + school_enroll +
                   female_labor + oda_pc + agriculture_gdp, data = mice_complete)
model_rf <- lm(health_exp_pc ~ gdp_pc + unemployment + school_enroll +
                 female_labor + oda_pc + agriculture_gdp, data = rf_complete)

model_comparison <- data.frame(
  Method = c("Baseline (Complete Cases)", "Mean", "KNN", "Regression", "MICE", "Random Forest"),
  R2 = c(summary(baseline_model)$r.squared,
         summary(model_mean)$r.squared,
         summary(model_knn)$r.squared,
         summary(model_reg)$r.squared,
         summary(model_mice)$r.squared,
         summary(model_rf)$r.squared)
)
print(model_comparison)
write_csv(model_comparison, "outputs/model_comparison_r2.csv")

#=========================
# 8. Compare regression coefficients
#=========================
extract_coefs <- function(model, method_name) {
  coefs <- summary(model)$coefficients
  data.frame(Method = method_name,
             Variable = rownames(coefs),
             Estimate = coefs[, "Estimate"],
             StdError = coefs[, "Std. Error"],
             PValue = coefs[, "Pr(>|t|)"])
}

coef_table <- bind_rows(
  extract_coefs(baseline_model, "Baseline (Complete Cases)"),
  extract_coefs(model_mean, "Mean"),
  extract_coefs(model_knn, "KNN"),
  extract_coefs(model_reg, "Regression"),
  extract_coefs(model_mice, "MICE"),
  extract_coefs(model_rf, "Random Forest")
)

write_csv(coef_table, "outputs/model_coefficients_comparison.csv")

#=========================
# 9. Plot RMSE comparison 
#=========================
ggplot(results_summary, aes(x = reorder(Method, RMSE), y = RMSE, fill = Method)) +
  geom_col() +
  geom_text(aes(label = round(RMSE, 3)), vjust = -0.5) +
  labs(title = "Imputation Method Comparison",
       x = "Method", y = "RMSE") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("outputs/imputation_results_real_rmse.png", width = 6, height = 4)

#=========================
# 10. Save baseline model summaries
#=========================

baseline_summary <- as.data.frame(summary(baseline_model)$coefficients)
baseline_summary$Variable <- rownames(baseline_summary)
rownames(baseline_summary) <- NULL

baseline_metrics <- data.frame(
  Metric = c("R2", "Adj_R2", "Residual_SE"),
  Value = c(summary(baseline_model)$r.squared,
            summary(baseline_model)$adj.r.squared,
            summary(baseline_model)$sigma)
)

write_csv(baseline_summary, "outputs/baseline_model_coefficients.csv")
write_csv(baseline_metrics, "outputs/baseline_model_metrics.csv")

print(baseline_summary)
print(baseline_metrics)

cat("\n Imputation Comparison (Real Missingness) Completed Successfully.\n")
