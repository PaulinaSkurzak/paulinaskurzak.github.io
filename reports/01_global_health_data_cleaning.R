# ============================================
# Data Cleaning & Preparation 
# Paulina Skurzak
# ============================================

#=========================
# 1.Load libraries and set working directory 
#=========================
library(wbstats)
library(dplyr)
library(naniar)
library(ggplot2)
library(readr)
library(tidyr)
library(MissMech)

setwd("~/Documents/GRI Research Fellow")
#=========================
# 2. Load merged dataset
#=========================
analysis_df <- read_csv("outputs/merged_analysis.csv")

#=========================
# 3. Add GDP data
#=========================
gdp_data <- wb_data(
  indicator   = "NY.GDP.PCAP.CD",
  start_date  = 2000,
  end_date    = 2025,
  return_wide = TRUE
) %>%
  select(iso3c, date, NY.GDP.PCAP.CD)

analysis_df <- left_join(analysis_df, gdp_data, by = c("iso3c", "date"))

#=========================
# 4. Clean column names
#=========================
names(analysis_df) <- make.unique(names(analysis_df))

#=========================
# 5. Match and rename columns
#=========================
get_first_match <- function(pattern) grep(pattern, names(analysis_df), value = TRUE)[1]

col_map <- c(
  health_exp_pc   = get_first_match("SH\\.XPD\\.CHEX"),
  gdp_pc          = get_first_match("NY\\.GDP\\.PCAP\\.CD"),
  unemployment    = get_first_match("SL\\.UEM\\.TOTL"),
  agriculture_gdp = get_first_match("NV\\.AGR"),
  school_enroll   = get_first_match("SE\\.SEC\\.ENRR"),
  female_labor    = get_first_match("SL\\.TLF\\.PART\\.FE\\.ZS"),
  oda_pc          = get_first_match("DT\\.ODA\\.ODAT")
)

cat("\n===== Column Mapping Used =====\n")
print(col_map)

#=========================
# 6. Create clean dataset
#=========================
analysis_clean <- analysis_df %>%
  select(iso3c, country, date, all_of(unlist(col_map))) %>%
  rename_with(~ names(col_map), -c(iso3c, country, date))

#=========================
# 7. Summarize missingness
#=========================
cat("\n===== Missing Data Summary =====\n")
miss_summary <- miss_var_summary(analysis_clean)
print(miss_summary)

#=========================
# 8. Test for MCAR
#=========================
cat("\n===== Testing MCAR =====\n")

numeric_data <- analysis_clean %>%
  select(where(is.numeric))

tryCatch({
  mcar_test <- MissMech::TestMCARNormality(numeric_data)
  print(mcar_test)
}, error = function(e) {
  cat("MissMech failed due to covariance singularity, using correlation fallback.\n")
  miss_corr <- cor(sapply(numeric_data, is.na), use = "pairwise.complete.obs")
  print(round(miss_corr, 2))
})

#=========================
# 9. Save clean dataset
#=========================
write_csv(analysis_clean, "outputs/analysis_clean.csv")
cat("\n Data Preparation Completed \n")
