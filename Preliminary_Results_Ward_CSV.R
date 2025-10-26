# ============================================================
# ESB — Minimal Loader (District-level-data.csv)
# Goal: clean, tiny, and ready for analysis
# ============================================================

# ---- Packages (install once if needed) ----
# install.packages(c("tidyverse","janitor","stringr","readr"))
suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(stringr)
  library(readr)
})

# ---- 0) Config ----
DIST_CSV <- "ESB_adoption_dataset_v9_update_june_2025_district_level.csv"          # <- your extracted CSV sheet
OUT_RDS  <- "esb_df_clean.rds"                 # local cache for reuse
OUT_CSV  <- "esb_df_clean.csv"                 # optional small, tidy export

stopifnot(file.exists(DIST_CSV))

# ---- 1) Helpers ----
to_num <- function(x) {
  if (is.numeric(x)) x else suppressWarnings(readr::parse_number(as.character(x)))
}
yesno_to_binary <- function(x) {
  xx <- as.character(x)
  case_when(
    xx %in% c("Yes","YES","yes","Y","1","true","TRUE") ~ 1L,
    xx %in% c("No","NO","no","N","0","false","FALSE")  ~ 0L,
    TRUE ~ NA_integer_
  )
}

# ---- 2) Read CSV & clean column names ----
dist_raw <- read_csv(DIST_CSV, guess_max = 100000, show_col_types = FALSE) %>%
  clean_names()

# ---- 3) Build tidy, typed dataset with outcome and predictors ----
esb_df <- dist_raw %>%
  mutate(
    # Target (robust): prefer explicit yes/no, fallback to count > 0
    has_comm_text = str_squish(str_to_lower(as.character(x0a_has_committed_es_bs))),
    committed_from_text = case_when(
      has_comm_text %in% c("yes","y","true","1","committed","has committed","has_committed") ~ 1L,
      has_comm_text %in% c("no","n","false","0","not committed","uncommitted") ~ 0L,
      has_comm_text %in% c("",".","na","n/a","null","none","missing") ~ NA_integer_,
      TRUE ~ NA_integer_
    ),
    esb_committed_cnt = to_num(x3a_number_of_es_bs_committed),
    committed_from_count = if_else(!is.na(esb_committed_cnt), as.integer(esb_committed_cnt > 0), NA_integer_),
    committed = coalesce(committed_from_text, committed_from_count)
  ) %>%
  transmute(
    # Identifiers
    lea_id        = x1c_lea_id,
    state         = x1a_state,
    district_name = x1b_local_education_agency_lea_or_entity_name,
    
    # Outcome
    esb_committed = esb_committed_cnt,
    committed     = as.integer(committed),
    
    # Geography / Operations
    total_buses = to_num(x2a_total_number_of_buses),
    contractor  = yesno_to_binary(x2b_contractor_used_for_some_or_all_of_buses),
    locale      = factor(x1p_locale_broad_type_name),
    region      = factor(x1q_census_region),
    
    # Demographics / Economics
    pct_frl          = to_num(x4e_percentage_of_students_in_district_eligible_for_free_or_reduced_price_lunch),
    pct_poverty      = to_num(x4g_percent_of_population_below_the_poverty_level),
    income_med       = to_num(x4f_median_household_income),
    pct_nonwhite_hisp= to_num(x5b_percent_non_white_and_or_hispanic),
    
    # Health / EJ
    pm25      = to_num(x5f_pm2_5_concentration),
    ozone     = to_num(x5h_ozone_concentration),
    asthma18p = to_num(x5l_average_rate_of_asthma_among_adults_aged_18_and_older),
    
    # Funding / Policy
    pod                = yesno_to_binary(x5q_wri_priority_outreach_district_pod),
    applied_not_awarded = yesno_to_binary(x6e_applied_for_esb_funding_but_not_awarded),
    arp_qualified       = yesno_to_binary(x5n_qualified_for_american_rescue_plan_funding),
    epa2022_prioritized = yesno_to_binary(x5o_epa_2022_clean_school_bus_rebate_program_prioritized_school_district),
    epa2023_prioritized = yesno_to_binary(x5p_epa_2023_clean_school_bus_grant_rebate_programs_prioritized_school_district),
    
    # Adoption intensity
    pct_esb_fleet = to_num(x3i_percent_of_fleet_that_is_electric)
  ) %>%
  filter(!is.na(lea_id)) %>%
  distinct(lea_id, .keep_all = TRUE)

# ---- 4) Save tidy local copies for quick reuse ----
saveRDS(esb_df, OUT_RDS)
readr::write_csv(esb_df, OUT_CSV)

# ---- 5) Quick sanity check ----
cat("\nSaved:\n -", OUT_RDS, "\n -", OUT_CSV, "\n")
cat("Rows:", nrow(esb_df), " | Cols:", ncol(esb_df), "\n\n")

esb_df %>%
  count(committed, drop = FALSE) %>%
  mutate(percent = round(100 * n / sum(n), 2)) %>%
  print()
