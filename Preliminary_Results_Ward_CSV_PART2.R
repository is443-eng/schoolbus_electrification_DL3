# ============================================================
# ESB — Family Summaries + Ultra-Robust Logistic Models
# Starts from: esb_df_clean.rds (created by your CSV loader)
# ============================================================

# ---- Packages (install once if needed) ----
# install.packages(c("tidyverse","broom"))
suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
})

# ---- 0) Load clean dataset ----
esb_df <- readRDS("esb_df_clean.rds")
stopifnot(is.data.frame(esb_df), "committed" %in% names(esb_df))

cat("Rows:", nrow(esb_df), " | Cols:", ncol(esb_df), "\n")
esb_df %>% count(committed, drop = FALSE) %>%
  mutate(pct = round(100*n/sum(n), 2)) %>% print()

# ---- 1) Measure families (from Charter/Research Design) ----
funding_vars <- c("applied_not_awarded","pod","arp_qualified","epa2022_prioritized","epa2023_prioritized")
demo_vars    <- c("pct_frl","pct_poverty","income_med","pct_nonwhite_hisp")
health_vars  <- c("pm25","ozone","asthma18p")
geo_vars     <- c("total_buses","contractor","locale","region","pct_esb_fleet")

# keep only vars that actually exist
funding_vars <- intersect(funding_vars, names(esb_df))
demo_vars    <- intersect(demo_vars, names(esb_df))
health_vars  <- intersect(health_vars, names(esb_df))
geo_vars     <- intersect(geo_vars, names(esb_df))

# ---- 2) Descriptive summaries by family (fast + readable) ----
summarize_family <- function(df, vars, title = "Measures") {
  present <- intersect(vars, names(df))
  if (!length(present)) { cat("\n", title, ": no variables present.\n"); return(invisible(NULL)) }
  
  num_vars <- present[vapply(df[present], is.numeric, logical(1))]
  fac_vars <- setdiff(present, num_vars)
  bin_vars <- num_vars[vapply(df[num_vars], function(x) all(na.omit(unique(x)) %in% c(0,1)), logical(1))]
  num_vars <- setdiff(num_vars, bin_vars)
  
  if (length(num_vars)) {
    cat("\n---", title, ": Numeric ---\n")
    df %>%
      summarize(across(all_of(num_vars),
                       list(n = ~sum(!is.na(.x)),
                            mean = ~mean(.x, na.rm = TRUE),
                            sd = ~sd(.x, na.rm = TRUE),
                            median = ~median(.x, na.rm = TRUE),
                            p25 = ~quantile(.x, .25, na.rm = TRUE),
                            p75 = ~quantile(.x, .75, na.rm = TRUE)
                       ), .names = "{.col}_{.fn}"
      )) %>%
      pivot_longer(everything(), names_to = c("variable",".value"), names_sep = "_(?=[^_]+$)") %>%
      arrange(variable) %>% print(n = Inf)
  }
  
  if (length(bin_vars)) {
    cat("\n---", title, ": Binary ---\n")
    df %>%
      summarize(across(all_of(bin_vars),
                       ~mean(.x == 1, na.rm = TRUE), .names = "pct_yes_{.col}"
      )) %>%
      pivot_longer(everything(), names_to = "variable", values_to = "pct_yes") %>%
      mutate(pct_yes = scales::percent(pct_yes, accuracy = 0.1)) %>%
      arrange(variable) %>% print(n = Inf)
  }
  
  if (length(fac_vars)) {
    cat("\n---", title, ": Factors (commitment rate by level) ---\n")
    invisible(lapply(fac_vars, function(v){
      cat("\n", v, "\n", sep = "")
      df %>%
        filter(!is.na(.data[[v]])) %>%
        group_by(.data[[v]]) %>%
        summarize(n = n(),
                  commit_rate = mean(committed == 1, na.rm = TRUE),
                  .groups = "drop") %>%
        arrange(desc(commit_rate), desc(n)) %>%
        mutate(commit_rate = scales::percent(commit_rate, accuracy = 0.1)) %>%
        print(n = Inf)
    }))
  }
}

summarize_family(esb_df, funding_vars, "Funding & Policy")
summarize_family(esb_df, demo_vars,    "Demographics & Economics")
summarize_family(esb_df, health_vars,  "Health & Environment")
summarize_family(esb_df, geo_vars,     "Geography & Operations")

# ============================================================
# 3) Ultra-robust logistic modeling (crash-proof helpers)
# ============================================================

# Identify binary (0/1) numerics
is_binary01 <- function(x) {
  if (!is.numeric(x)) return(FALSE)
  ux <- sort(unique(na.omit(x)))
  length(ux) <= 2 && all(ux %in% c(0,1))
}

# Build clean complete-case dataset and prune weak predictors
prep_model_data <- function(df, vars, family_name = "Model",
                            min_rows = 200,    # tweak for your data size
                            min_level_n = 20)  # min support per class/level
{
  vars <- intersect(vars, names(df))
  if (!length(vars)) { message("[", family_name, "] No variables present."); return(NULL) }
  
  dfm <- df |> dplyr::select(committed, dplyr::all_of(vars)) |> tidyr::drop_na()
  n0  <- nrow(dfm)
  message("[", family_name, "] complete cases: ", n0)
  if (n0 == 0) { message("  -> 0 rows after drop_na on: ", paste(vars, collapse=", ")); return(NULL) }
  
  y <- dfm$committed
  if (!is.numeric(y)) y <- as.numeric(y)
  if (!all(y %in% c(0,1))) { message("  -> Outcome not coded 0/1. Unique: ", paste(unique(y), collapse=", ")); return(NULL) }
  if (length(unique(y)) < 2) { message("  -> Outcome has one class after filtering."); return(NULL) }
  
  keep <- vapply(dfm[names(dfm) != "committed"], function(x){
    if (is.factor(x)) length(unique(x)) > 1 else stats::sd(as.numeric(x), na.rm=TRUE) > 0
  }, logical(1))
  dropped_zv <- names(dfm)[-1][!keep]
  if (length(dropped_zv)) message("  -> Dropped zero-variance: ", paste(dropped_zv, collapse=", "))
  dfm <- dplyr::select(dfm, committed, dplyr::all_of(names(keep)[keep]))
  
  bin_vars <- names(dfm)[-1][vapply(dfm[-1], is_binary01, logical(1))]
  bad_bin  <- c()
  for (v in bin_vars) {
    tab <- table(dfm[[v]])
    if (length(tab) < 2 || any(tab < min_level_n)) bad_bin <- c(bad_bin, v)
  }
  if (length(bad_bin)) {
    message("  -> Dropped sparse/all-one-level binaries: ", paste(bad_bin, collapse=", "))
    dfm <- dplyr::select(dfm, -dplyr::all_of(bad_bin))
  }
  
  fac_vars <- names(dfm)[-1][vapply(dfm[-1], is.factor, logical(1))]
  bad_fac  <- c()
  for (v in fac_vars) {
    tab <- table(dfm[[v]])
    if (length(tab) < 2 || any(tab < min_level_n)) bad_fac <- c(bad_fac, v)
  }
  if (length(bad_fac)) {
    message("  -> Dropped sparse factor vars: ", paste(bad_fac, collapse=", "))
    dfm <- dplyr::select(dfm, -dplyr::all_of(bad_fac))
  }
  
  if (nrow(dfm) < min_rows) message("  -> Note: only ", nrow(dfm), " rows meet criteria (min_rows=", min_rows, ").")
  if (ncol(dfm) <= 1) { message("  -> No predictors remain after pruning."); return(NULL) }
  
  droplevels(dfm)
}

fit_logit_guarded <- function(df, vars, family_name = "Model") {
  df_cc <- prep_model_data(df, vars, family_name)
  if (is.null(df_cc)) return(NULL)
  
  kept_vars <- setdiff(names(df_cc), "committed")
  fml <- as.formula(paste("committed ~", paste(kept_vars, collapse = " + ")))
  
  df_cc$committed <- as.numeric(df_cc$committed)
  
  mf <- tryCatch(stats::model.frame(fml, data = df_cc, na.action = stats::na.fail),
                 error = function(e) { message("  -> model.frame error: ", e$message); NULL })
  if (is.null(mf) || nrow(mf) == 0 || ncol(mf) <= 1) {
    message("  -> Model frame empty after encoding."); return(NULL)
  }
  
  suppressWarnings(
    tryCatch(
      stats::glm(fml, data = df_cc, family = stats::binomial(), model = FALSE),
      error = function(e) { message("  -> GLM error: ", e$message); NULL }
    )
  )
}

tidy_or <- function(model) {
  if (is.null(model)) return(tibble(term=character(), OR=numeric(), CI_low=numeric(), CI_high=numeric(), p.value=numeric()))
  broom::tidy(model, conf.int = TRUE, exponentiate = TRUE) |>
    dplyr::mutate(dplyr::across(c(estimate, conf.low, conf.high, p.value), ~round(.x, 4))) |>
    dplyr::rename(OR = estimate, CI_low = conf.low, CI_high = conf.high) |>
    dplyr::arrange(dplyr::desc(abs(log(OR + 1e-9))))
}

# ---- 4) Run models per family ----
cat("\n=== Logistic: Funding & Policy ===\n")
m_funding  <- fit_logit_guarded(esb_df, funding_vars, "Funding & Policy")
print(tidy_or(m_funding), n = Inf)

cat("\n=== Logistic: Demographics & Economics ===\n")
m_demo     <- fit_logit_guarded(esb_df, demo_vars,    "Demographics & Economics")
print(tidy_or(m_demo), n = Inf)

cat("\n=== Logistic: Health & Environment ===\n")
m_health   <- fit_logit_guarded(esb_df, health_vars,  "Health & Environment")
print(tidy_or(m_health), n = Inf)

cat("\n=== Logistic: Geography & Operations ===\n")
m_geo      <- fit_logit_guarded(esb_df, geo_vars,     "Geography & Operations")
print(tidy_or(m_geo), n = Inf)

# ---- 5) Compact combined model ----
combined_vars <- intersect(c(
  "total_buses","contractor","pct_frl","pct_poverty","income_med","pct_nonwhite_hisp",
  "pm25","asthma18p","pod","applied_not_awarded","arp_qualified",
  "epa2022_prioritized","epa2023_prioritized","region","locale"
), names(esb_df))

cat("\n=== Logistic: Compact Combined ===\n")
m_combined <- fit_logit_guarded(esb_df, combined_vars, "Combined")
print(tidy_or(m_combined), n = Inf)

if (!is.null(m_combined)) {
  pred <- predict(m_combined, type = "response")
  acc  <- mean((pred >= 0.5) == (esb_df$committed == 1), na.rm = TRUE)
  cat("Combined model in-sample accuracy @0.5:", round(acc, 3), "\n")
}
