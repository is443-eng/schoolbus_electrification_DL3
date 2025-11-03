# ===============================================================
# Electric School Bus — Modeling + Policy Simulation w/ 95% CIs
# + Scenario composition + Savings over 5 years (incl. Current)
# ===============================================================

library(tidyverse)
library(ggeffects)
library(wesanderson)
library(ggrepel)

set.seed(123)

# ---------- Helpers ----------
parse_percent <- function(x) {
  vals <- readr::parse_number(as.character(x))
  if (is.finite(max(vals, na.rm = TRUE)) && max(vals, na.rm = TRUE) > 1.5) vals <- vals / 100
  vals
}
safe_weighted_mean <- function(x, w) {
  idx <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w) & (w > 0)
  if (!any(idx)) return(NA_real_)
  sum(x[idx] * w[idx]) / sum(w[idx])
}

# ---------- Data load & clean ----------
df <- read_csv("esb_df_clean.csv")

df_clean <- df %>%
  mutate(
    epa2022_prioritized = ifelse(epa2022_prioritized %in% c("Yes","Y",1,TRUE), 1L, 0L),
    epa2023_prioritized = ifelse(epa2023_prioritized %in% c("Yes","Y",1,TRUE), 1L, 0L),
    prioritized_any      = ifelse(epa2022_prioritized == 1 | epa2023_prioritized == 1, 1L, 0L),
    pct_esb_fleet_prop   = parse_percent(pct_esb_fleet),
    pct_poverty          = as.numeric(pct_poverty),
    asthma18p            = as.numeric(asthma18p),
    total_buses          = as.numeric(total_buses),
    locale               = as.character(locale)
  ) %>%
  # Drop Unknown / NA locales
  filter(!is.na(locale), locale != "Unknown") %>%
  mutate(locale = droplevels(as.factor(locale))) %>%
  # keep valid outcome
  filter(!is.na(pct_esb_fleet_prop), pct_esb_fleet_prop >= 0, pct_esb_fleet_prop <= 1)

# ---------- Current (status-quo) prioritization by locale ----------
current_dist <- df_clean %>%
  group_by(locale) %>%
  summarise(n_districts = n(),
            pct_prioritized = mean(prioritized_any, na.rm = TRUE) * 100,
            .groups = "drop") %>%
  arrange(desc(pct_prioritized))
current_dist

# ---------- GLM with interaction ----------
glm_any_locale <- glm(
  pct_esb_fleet_prop ~ prioritized_any * locale + pct_poverty + asthma18p,
  data = df_clean,
  family = quasibinomial(link = "logit")
)
summary(glm_any_locale)

# ---------- Predicted electrification by locale × prioritization ----------
pred_locale_df <- as.data.frame(
  ggeffect(glm_any_locale, terms = c("locale", "prioritized_any"))
) %>% rename(locale = x, prioritized_any = group)

effect_by_locale <- pred_locale_df %>%
  select(locale, prioritized_any, predicted) %>%
  mutate(prioritized_any = as.numeric(as.character(prioritized_any))) %>%
  pivot_wider(names_from = prioritized_any, values_from = predicted, names_prefix = "prioritized_") %>%
  mutate(
    difference_pp      = (prioritized_1 - prioritized_0) * 100,
    pct_nonprioritized = prioritized_0 * 100,
    pct_prioritized    = prioritized_1 * 100
  ) %>% arrange(desc(difference_pp))
effect_by_locale

# ---------- Fleet-weighted national averages (+ bootstrap CI for diff) ----------
national_avg <- df_clean %>%
  filter(!is.na(total_buses), total_buses > 0) %>%
  group_by(prioritized_any) %>%
  summarise(weighted_mean = safe_weighted_mean(pct_esb_fleet_prop, total_buses),
            .groups = "drop") %>%
  mutate(label = ifelse(prioritized_any == 1, "Prioritized Avg", "Not Prioritized Avg"))
national_avg

boot_diff <- function(df, n_boot = 5000) {
  boot_est <- replicate(n_boot, {
    samp <- df[sample(nrow(df), replace = TRUE), ]
    wmeans <- samp %>%
      group_by(prioritized_any) %>%
      summarise(wm = safe_weighted_mean(pct_esb_fleet_prop, total_buses), .groups = "drop") %>%
      arrange(prioritized_any)
    wmeans$wm[2] - wmeans$wm[1]
  })
  tibble(
    diff_mean = mean(boot_est),
    diff_low  = quantile(boot_est, 0.025),
    diff_high = quantile(boot_est, 0.975),
    diff_pct_points = diff_mean * 100,
    diff_low_pp     = diff_low  * 100,
    diff_high_pp    = diff_high * 100
  )
}
national_test <- boot_diff(df_clean %>% filter(!is.na(total_buses), total_buses > 0))
national_test

# ---------- Visualization (Darjeeling palette; labels far left) ----------
palette_darjeeling <- wes_palette("Darjeeling1", 2, type = "discrete")
leftmost_label_x <- -0.5

ggplot(pred_locale_df,
       aes(x = locale, y = predicted * 100,
           group = prioritized_any, color = as.factor(prioritized_any))) +
  geom_point(aes(shape = as.factor(prioritized_any)), size = 3) +
  geom_line(aes(linetype = as.factor(prioritized_any)), linewidth = 1) +
  geom_errorbar(aes(ymin = conf.low * 100, ymax = conf.high * 100),
                width = 0.2, alpha = 0.6) +
  geom_hline(
    data = national_avg,
    aes(yintercept = weighted_mean * 100, color = as.factor(prioritized_any)),
    linetype = "dashed", linewidth = 0.9, show.legend = FALSE
  ) +
  geom_text(
    data = national_avg,
    aes(x = leftmost_label_x,
        y = weighted_mean * 100,
        label = paste0(label, ": ", round(weighted_mean * 100, 2), "%"),
        color = as.factor(prioritized_any)),
    hjust = 0, vjust = -0.4, size = 4, fontface = "bold", show.legend = FALSE
  ) +
  scale_color_manual(values = palette_darjeeling, labels = c("Not Prioritized","Prioritized")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Predicted Electrification by Locale and EPA Prioritization",
    subtitle = paste0(
      "Controls: poverty rate & asthma prevalence\n",
      "Dashed lines = fleet-weighted national averages; ",
      "bootstrap 95% CI: ",
      round(national_test$diff_pct_points, 2), " pp [",
      round(national_test$diff_low_pp, 2), ", ",
      round(national_test$diff_high_pp, 2), "]"
    ),
    x = "Locale", y = "Predicted % of Fleet Electrified",
    color = "EPA Status", shape = "EPA Status", linetype = "EPA Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title  = element_text(face = "bold", hjust = 0.5),
    plot.margin = margin(10, 60, 10, 110)
  ) +
  coord_cartesian(clip = "off")

# ===============================================================
# Policy Shift Simulation: add 95% CIs for scenario outcomes
# ===============================================================

set.seed(42)

# overall rate in data
overall_prior_rate <- mean(df_clean$prioritized_any, na.rm = TRUE)
locale_sizes <- df_clean %>% count(locale, name = "n")

rescale_policy_to_overall <- function(policy_tbl, locale_sizes, target_overall) {
  tmp <- policy_tbl %>% left_join(locale_sizes, by = "locale")
  current_overall <- sum(tmp$target_rate * tmp$n) / sum(tmp$n)
  k0 <- if (is.finite(current_overall) && current_overall > 0) target_overall / current_overall else 1
  policy_rescaled <- tmp %>%
    mutate(target_rate_adj = pmin(pmax(target_rate * k0, 0), 1)) %>%
    select(locale, target_rate_adj)
  final_overall <- sum(policy_rescaled$target_rate_adj * tmp$n) / sum(tmp$n)
  policy_rescaled %>%
    mutate(target_rate_adj = pmin(pmax(target_rate_adj * (target_overall / final_overall), 0), 1))
}

# scenarios (pre-rescale)
policy_rural <- tibble(
  locale = levels(df_clean$locale),
  target_rate = case_when(
    locale == "Rural"    ~ 0.60,
    locale == "Town"     ~ 0.40,
    locale == "Suburban" ~ 0.30,
    locale == "Urban"    ~ 0.20,
    TRUE ~ 0.40
  )
)
policy_balanced <- tibble(locale = levels(df_clean$locale), target_rate = 0.40)
policy_urban <- tibble(
  locale = levels(df_clean$locale),
  target_rate = case_when(
    locale == "Rural"    ~ 0.20,
    locale == "Town"     ~ 0.30,
    locale == "Suburban" ~ 0.40,
    locale == "Urban"    ~ 0.60,
    TRUE ~ 0.40
  )
)

# rescale to keep overall rate constant
policy_rural_adj    <- rescale_policy_to_overall(policy_rural,    locale_sizes, overall_prior_rate)
policy_balanced_adj <- rescale_policy_to_overall(policy_balanced, locale_sizes, overall_prior_rate)
policy_urban_adj    <- rescale_policy_to_overall(policy_urban,    locale_sizes, overall_prior_rate)

# -------- Bootstrap functions for scenario CIs --------
boot_statusquo_pred <- function(model, df, n_boot = 2000) {
  est <- replicate(n_boot, {
    samp <- df[sample(nrow(df), replace = TRUE), ]
    preds <- predict(model, newdata = samp, type = "response")
    safe_weighted_mean(preds, samp$total_buses)
  })
  tibble(mean = mean(est), low = quantile(est, .025), high = quantile(est, .975))
}
boot_current_observed <- function(df, n_boot = 2000) {
  est <- replicate(n_boot, {
    samp <- df[sample(nrow(df), replace = TRUE), ]
    safe_weighted_mean(samp$pct_esb_fleet_prop, samp$total_buses)
  })
  tibble(mean = mean(est), low = quantile(est, .025), high = quantile(est, .975))
}
boot_policy_pred <- function(model, df, policy_adj, n_boot = 2000) {
  est <- replicate(n_boot, {
    samp <- df[sample(nrow(df), replace = TRUE), ] %>%
      left_join(policy_adj, by = "locale") %>%
      mutate(
        prioritized_any = rbinom(n(), 1, p = target_rate_adj),
        locale = factor(locale, levels = levels(df_clean$locale))
      )
    preds <- predict(model, newdata = samp, type = "response")
    safe_weighted_mean(preds, samp$total_buses)
  })
  tibble(mean = mean(est), low = quantile(est, .025), high = quantile(est, .975))
}

# -------- Compute scenario estimates + CIs --------
ci_current  <- boot_current_observed(df_clean)
ci_statusqo <- boot_statusquo_pred(glm_any_locale, df_clean)
ci_rural    <- boot_policy_pred(glm_any_locale, df_clean, policy_rural_adj)
ci_bal      <- boot_policy_pred(glm_any_locale, df_clean, policy_balanced_adj)
ci_urban    <- boot_policy_pred(glm_any_locale, df_clean, policy_urban_adj)

scenario_results <- tribble(
  ~scenario,                   ~mean,            ~low,             ~high,
  "Current Funding Pattern",    ci_current$mean,  ci_current$low,  ci_current$high,
  "Status Quo (Predicted)",     ci_statusqo$mean,ci_statusqo$low, ci_statusqo$high,
  "Rural-focused",              ci_rural$mean,   ci_rural$low,    ci_rural$high,
  "Balanced",                   ci_bal$mean,     ci_bal$low,      ci_bal$high,
  "Urban-focused",              ci_urban$mean,   ci_urban$low,    ci_urban$high
) %>%
  mutate(across(c(mean, low, high), ~ . * 100)) %>%
  rename(weighted_electrification_pct = mean,
         ci_low_pct = low,
         ci_high_pct = high)

scenario_results

# -------- Plot scenario comparison with 95% CIs --------
pal5 <- wes_palette("Darjeeling1", 5, type = "discrete")

ggplot(scenario_results,
       aes(x = scenario, y = weighted_electrification_pct, fill = scenario)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = ci_low_pct, ymax = ci_high_pct), width = 0.15, linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f%%", weighted_electrification_pct),
                y = ci_low_pct - 0.08 * (max(weighted_electrification_pct) - min(ci_low_pct))),
            vjust = 1, size = 5, fontface = "bold") +
  scale_fill_manual(values = pal5) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand = expansion(mult = c(0.02, 0.12))) +
  labs(
    title = "National Fleet Electrification: Current vs. Alternative Prioritization Policies",
    subtitle = paste0(
      "Bars show fleet-weighted % electrified; error bars = bootstrap 95% CI (n=2000)\n",
      "Simulated policies are rescaled to keep overall prioritization rate at ",
      round(100 * overall_prior_rate, 1), "%"
    ),
    x = NULL, y = "Fleet-weighted % of Buses Electrified"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# ===============================================================
# Scenario composition — 100% stacked bars (clear alternative to pies)
# ===============================================================
scenario_alloc <- function(policy_adj, scen_name) {
  locale_sizes %>%
    left_join(policy_adj, by = "locale") %>%
    mutate(expected_prioritized = n * target_rate_adj) %>%
    summarise(locale, share = expected_prioritized / sum(expected_prioritized),
              .groups = "drop") %>%
    mutate(scenario = scen_name)
}
current_alloc <- df_clean %>%
  group_by(locale) %>%
  summarise(prioritized = sum(prioritized_any, na.rm = TRUE), .groups = "drop") %>%
  mutate(share = prioritized / sum(prioritized),
         scenario = "Current Funding Pattern") %>%
  select(locale, share, scenario)

shares_df <- bind_rows(
  current_alloc,
  scenario_alloc(policy_rural_adj,    "Rural-focused"),
  scenario_alloc(policy_balanced_adj, "Balanced"),
  scenario_alloc(policy_urban_adj,    "Urban-focused")
) %>%
  mutate(
    locale   = factor(locale, levels = levels(df_clean$locale)),
    scenario = factor(
      scenario,
      levels = c("Current Funding Pattern", "Rural-focused", "Balanced", "Urban-focused")
    ),
    pct_lab = scales::percent(share, accuracy = 0.1)
  )

palette_locales <- c("Rural"="#e41a1c", "Suburban"="#1b9e77",
                     "Town"="#ffb000", "Urban"="#d95f02")

labels_inside <- shares_df %>% filter(share >= 0.06)

ggplot(shares_df, aes(x = scenario, y = share, fill = locale)) +
  geom_col(width = 0.7, color = "white") +
  geom_text(
    data = labels_inside,
    aes(label = paste0(locale, ": ", pct_lab)),
    position = position_stack(vjust = 0.5),
    size = 3.8, fontface = "bold", color = "black"
  ) +
  scale_fill_manual(values = palette_locales) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "How Each Scenario Allocates Prioritization Across Locales",
    subtitle = "Each bar sums to 100%: share of prioritized districts by locale",
    x = NULL, y = "Share of prioritized districts",
    fill = "Locale"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# ===============================================================
# LIFETIME SAVINGS (Parametric bootstrap; aligned complete cases)
# Deterministic policies; assumes $100,000 lifetime savings / ESB
# ===============================================================

library(MASS)     # mvrnorm
library(scales)
library(ggplot2)
library(wesanderson)
library(dplyr)

lifetime_savings_per_bus <- 100000
B <- 5000  # increase for smoother CIs

# ---- Total buses held fixed (from raw df) ----
total_buses_total <- df %>%
  mutate(total_buses = as.numeric(total_buses)) %>%
  summarise(total_buses_total = sum(total_buses, na.rm = TRUE)) %>%
  pull(total_buses_total)

cat("Total school buses in service:", format(total_buses_total, big.mark=","), "\n")

# ---- Use a single complete-case mask consistent with the GLM ----
vars_needed <- c("prioritized_any","locale","pct_poverty","asthma18p","total_buses")
keep <- complete.cases(df_clean[, vars_needed])

dat_cc <- df_clean[keep, , drop = FALSE]  # filtered data used everywhere below
dat_cc$locale <- droplevels(dat_cc$locale)

# ---- Policy target rates matched to filtered rows ----
r_rural    <- dat_cc %>% left_join(policy_rural_adj,    by = "locale") %>% pull(target_rate_adj)
r_balanced <- dat_cc %>% left_join(policy_balanced_adj, by = "locale") %>% pull(target_rate_adj)
r_urban    <- dat_cc %>% left_join(policy_urban_adj,    by = "locale") %>% pull(target_rate_adj)

# ---- Model matrices built on the SAME filtered data ----
form <- formula(glm_any_locale)
X_curr <- model.matrix(form, dat_cc)

dat_p1 <- dat_cc; dat_p1$prioritized_any <- 1L
dat_p0 <- dat_cc; dat_p0$prioritized_any <- 0L
X_p1   <- model.matrix(form, dat_p1)
X_p0   <- model.matrix(form, dat_p0)

# ---- Fleet weights aligned to filtered rows ----
w <- as.numeric(dat_cc$total_buses)

# ---- Safe weighted mean ----
safe_wmean <- function(p, w) {
  ok <- is.finite(p) & is.finite(w) & !is.na(p) & !is.na(w) & w > 0
  sum(p[ok] * w[ok]) / sum(w[ok])
}

# ---- Parametric bootstrap draws of coefficients ----
beta_hat <- coef(glm_any_locale)
V_hat    <- vcov(glm_any_locale)
betas    <- MASS::mvrnorm(n = B, mu = beta_hat, Sigma = V_hat)

# ---- Storage
sav_current   <- numeric(B)
sav_statusquo <- numeric(B)
sav_rural     <- numeric(B)
sav_balanced  <- numeric(B)
sav_urban     <- numeric(B)

factor_dollars <- total_buses_total * lifetime_savings_per_bus

for (b in seq_len(B)) {
  bvec <- betas[b, ]
  
  # Probabilities (as numeric vectors; drop matrix dims)
  p_curr <- plogis(as.numeric(X_curr %*% bvec))
  p1     <- plogis(as.numeric(X_p1   %*% bvec))
  p0     <- plogis(as.numeric(X_p0   %*% bvec))
  
  # Deterministic policy expectations (same length as weights)
  p_rural    <- r_rural    * p1 + (1 - r_rural)    * p0
  p_bal      <- r_balanced * p1 + (1 - r_balanced) * p0
  p_urban    <- r_urban    * p1 + (1 - r_urban)    * p0
  
  # Fleet-weighted national proportions
  pr_curr <- safe_wmean(p_curr, w)
  pr_sq   <- pr_curr
  pr_rur  <- safe_wmean(p_rural, w)
  pr_bal  <- safe_wmean(p_bal,   w)
  pr_urb  <- safe_wmean(p_urban, w)
  
  # Convert to lifetime savings in dollars
  sav_current[b]   <- pr_curr * factor_dollars
  sav_statusquo[b] <- pr_sq   * factor_dollars
  sav_rural[b]     <- pr_rur  * factor_dollars
  sav_balanced[b]  <- pr_bal  * factor_dollars
  sav_urban[b]     <- pr_urb  * factor_dollars
}

# ---- Summaries
sumr <- function(x) tibble(mean = mean(x), low = quantile(x, .025), high = quantile(x, .975))

life_current   <- sumr(sav_current)
life_statusquo <- sumr(sav_statusquo)
life_rural     <- sumr(sav_rural)
life_balanced  <- sumr(sav_balanced)
life_urban     <- sumr(sav_urban)

lifetime_results <- tribble(
  ~scenario,                   ~mean,              ~low,               ~high,
  "Current Funding Pattern (Model-based)", life_current$mean, life_current$low, life_current$high,
  "Status Quo (Predicted)",               life_statusquo$mean, life_statusquo$low, life_statusquo$high,
  "Rural-focused",                        life_rural$mean,     life_rural$low,     life_rural$high,
  "Balanced",                             life_balanced$mean,  life_balanced$low,  life_balanced$high,
  "Urban-focused",                        life_urban$mean,     life_urban$low,     life_urban$high
) %>%
  mutate(scenario = factor(scenario,
                           levels = c("Current Funding Pattern (Model-based)","Status Quo (Predicted)",
                                      "Rural-focused","Balanced","Urban-focused")))

# ---- Plot A: Lifetime savings with parametric 95% CIs ----
pal5 <- wes_palette("Darjeeling1", 5, type = "discrete")

ggplot(lifetime_results, aes(x = scenario, y = mean, fill = scenario)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.15, linewidth = 0.8) +
  geom_text(aes(label = dollar(mean)), vjust = -0.5, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = pal5) +
  scale_y_continuous(labels = dollar_format(),
                     expand = expansion(mult = c(0.02, 0.15))) +
  labs(
    title = "Lifetime Savings by Scenario (Parametric 95% CIs)",
    subtitle = "GLM-coefficient uncertainty only; deterministic policy expectations",
    x = NULL, y = "Total Lifetime Savings (USD)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5))

# ---- Plot B: Δ vs Current (paired parametric bootstrap) ----
d_statusquo <- sav_statusquo - sav_current
d_rural     <- sav_rural     - sav_current
d_balanced  <- sav_balanced  - sav_current
d_urban     <- sav_urban     - sav_current

delta_results <- tribble(
  ~scenario,                ~mean,                    ~low,                         ~high,
  "Status Quo (Predicted)", mean(d_statusquo), quantile(d_statusquo, .025), quantile(d_statusquo, .975),
  "Rural-focused",          mean(d_rural),     quantile(d_rural,    .025), quantile(d_rural,    .975),
  "Balanced",               mean(d_balanced),  quantile(d_balanced, .025), quantile(d_balanced, .975),
  "Urban-focused",          mean(d_urban),     quantile(d_urban,    .025), quantile(d_urban,    .975)
) %>%
  mutate(scenario = factor(scenario,
                           levels = c("Status Quo (Predicted)","Rural-focused","Balanced","Urban-focused")))

ggplot(delta_results, aes(x = scenario, y = mean,
                          fill = ifelse(mean > 0, "Increase", "Decrease"))) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.15, linewidth = 0.8) +
  geom_text(aes(label = ifelse(mean > 0, paste0("+", dollar(mean)), dollar(mean))),
            vjust = ifelse(delta_results$mean >= 0, -0.5, 1.2),
            size = 4.2, fontface = "bold") +
  scale_fill_manual(values = c("Increase" = "#1b9e77", "Decrease" = "#d95f02")) +
  scale_y_continuous(labels = dollar_format(), expand = expansion(mult = c(0.05, 0.18))) +
  labs(
    title = "Δ Lifetime Savings vs Current (Parametric 95% CIs)",
    subtitle = "Paired draws from GLM coefficient uncertainty only",
    x = NULL, y = "Δ Lifetime Savings (USD)"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

# ===============================================================
# Sensitivity table: shift focus from Rural → Urban while
# keeping the overall prioritization rate constant
# ===============================================================

library(dplyr)
library(tidyr)
library(scales)
library(ggplot2)
library(purrr)

lifetime_savings_per_bus <- 100000

# Total buses (held constant)
total_buses_total <- df %>%
  mutate(total_buses = as.numeric(total_buses)) %>%
  summarise(total_buses_total = sum(total_buses, na.rm = TRUE)) %>%
  pull(total_buses_total)

# consistent complete cases with the fitted GLM
vars_needed <- c("prioritized_any","locale","pct_poverty","asthma18p","total_buses")
keep <- complete.cases(df_clean[, vars_needed])

dat_cc <- df_clean[keep, , drop = FALSE] %>%
  mutate(locale = droplevels(locale))

# design matrices once
form <- formula(glm_any_locale)
X_p1 <- model.matrix(form, transform(dat_cc, prioritized_any = 1L))
X_p0 <- model.matrix(form, transform(dat_cc, prioritized_any = 0L))
beta_hat <- coef(glm_any_locale)

# model-based predictions
p1 <- plogis(as.numeric(X_p1 %*% beta_hat))
p0 <- plogis(as.numeric(X_p0 %*% beta_hat))
w  <- as.numeric(dat_cc$total_buses)

safe_wmean <- function(p, w) {
  ok <- is.finite(p) & is.finite(w) & !is.na(p) & !is.na(w) & w > 0
  sum(p[ok] * w[ok]) / sum(w[ok])
}

# helper: rescale to keep national rate the same
rescale_policy_to_overall <- function(policy_tbl, locale_sizes, target_overall) {
  tmp <- policy_tbl %>% left_join(locale_sizes, by = "locale")
  current_overall <- sum(tmp$target_rate * tmp$n) / sum(tmp$n)
  k0 <- if (is.finite(current_overall) && current_overall > 0) target_overall / current_overall else 1
  policy_rescaled <- tmp %>%
    mutate(target_rate_adj = pmin(pmax(target_rate * k0, 0), 1))
  final_overall <- sum(policy_rescaled$target_rate_adj * tmp$n) / sum(tmp$n)
  policy_rescaled %>%
    mutate(target_rate_adj = pmin(pmax(target_rate_adj * (target_overall / final_overall), 0), 1)) %>%
    dplyr::select(locale, target_rate_adj)
}

# base endpoints
base_rural <- tibble(
  locale = levels(dat_cc$locale),
  target_rate = case_when(
    locale == "Rural"    ~ 0.60,
    locale == "Town"     ~ 0.40,
    locale == "Suburban" ~ 0.30,
    locale == "Urban"    ~ 0.20,
    TRUE ~ 0.40
  )
)
base_urban <- tibble(
  locale = levels(dat_cc$locale),
  target_rate = case_when(
    locale == "Rural"    ~ 0.20,
    locale == "Town"     ~ 0.30,
    locale == "Suburban" ~ 0.40,
    locale == "Urban"    ~ 0.60,
    TRUE ~ 0.40
  )
)

# linear "focus" blend: 0 = Rural-focused, 1 = Urban-focused
focus_grid <- tibble(focus = seq(0, 1, by = 0.05))

blend_policy <- function(focus) {
  base_rural %>%
    rename(tr_rural = target_rate) %>%
    left_join(rename(base_urban, tr_urban = target_rate), by = "locale") %>%
    mutate(target_rate = (1 - focus) * tr_rural + focus * tr_urban) %>%
    dplyr::select(locale, target_rate)
}

# compute grid of results (using map_dfr instead of do)
sens_rows <- map_dfr(focus_grid$focus, function(f) {
  pol_raw <- blend_policy(f)
  pol_adj <- rescale_policy_to_overall(pol_raw, locale_sizes, overall_prior_rate)
  r_vec <- dat_cc %>% left_join(pol_adj, by = "locale") %>% pull(target_rate_adj)
  p_exp <- r_vec * p1 + (1 - r_vec) * p0
  prop  <- safe_wmean(p_exp, w)
  tibble(
    focus = f,
    predicted_electrification_pct = 100 * prop,
    expected_electric_buses = round(prop * total_buses_total),
    lifetime_savings = prop * total_buses_total * lifetime_savings_per_bus
  )
})

sens_rows <- sens_rows %>%
  mutate(
    lifetime_savings_fmt = dollar(lifetime_savings),
    scenario_label = case_when(
      focus == 0   ~ "Rural-focused",
      focus == 0.5 ~ "Balanced",
      focus == 1   ~ "Urban-focused",
      TRUE ~ NA_character_
    )
  )

# include observed Current Funding Pattern for reference
obs_prop <- safe_wmean(dat_cc$pct_esb_fleet_prop, w)
current_row <- tibble(
  focus = NA_real_,
  predicted_electrification_pct = 100 * obs_prop,
  expected_electric_buses = round(obs_prop * total_buses_total),
  lifetime_savings = obs_prop * total_buses_total * lifetime_savings_per_bus,
  lifetime_savings_fmt = dollar(obs_prop * total_buses_total * lifetime_savings_per_bus),
  scenario_label = "Current Funding Pattern"
)

sensitivity_table <- bind_rows(current_row, sens_rows) %>% arrange(focus)

# ---- View first few rows ----
print(sensitivity_table %>%
        mutate(across(where(is.numeric), ~ round(., 3))) %>%
        head(15))

# ---- Plot sensitivity sweep ----
ggplot(sens_rows, aes(x = focus, y = predicted_electrification_pct)) +
  geom_line(linewidth = 1.1, color = "#1b9e77") +
  geom_point(data = subset(sens_rows, focus %in% c(0, 0.5, 1)), size = 3, color = "#d95f02") +
  geom_hline(yintercept = 100 * obs_prop, linetype = "dashed") +
  annotate("text", x = 0.02, y = 100 * obs_prop,
           label = paste0("Current Funding Pattern: ",
                          round(100 * obs_prop, 2), "%"),
           hjust = 0, vjust = -0.6) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("Rural-focused", "Balanced", "Urban-focused")
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Sensitivity of National Electrification to Policy Focus Mix",
    subtitle = "Focus = 0 (Rural) → 1 (Urban). Overall prioritization rate held constant.",
    x = "Policy Focus Mix",
    y = "Predicted % of Fleet Electrified (fleet-weighted)"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))



# ===============================================================
# Confidence bands for the sensitivity curve (parametric bootstrap)
#  - Uncertainty from GLM coefficients only (fast & tight CIs)
#  - Uses the same objects defined in the sensitivity section:
#    dat_cc, X_p1, X_p0, w, beta_hat, safe_wmean,
#    blend_policy(), rescale_policy_to_overall(), focus_grid,
#    locale_sizes, overall_prior_rate
# ===============================================================

library(MASS)   # mvrnorm
library(purrr)

# 1) Precompute adjusted target-rate vectors for each focus (aligned to dat_cc rows)
focuses <- focus_grid$focus
rvec_list <- setNames(
  lapply(focuses, function(f) {
    pol_raw <- blend_policy(f)
    pol_adj <- rescale_policy_to_overall(pol_raw, locale_sizes, overall_prior_rate)
    dat_cc %>% left_join(pol_adj, by = "locale") %>% pull(target_rate_adj)
  }),
  nm = focuses
)

# 2) Parametric bootstrap draws of GLM coefficients
B <- 2000  # increase to 5000+ for even smoother ribbons
V_hat <- vcov(glm_any_locale)
betas <- MASS::mvrnorm(n = B, mu = beta_hat, Sigma = V_hat)

# 3) For each draw, compute p1, p0 once; then expected prop for each focus
props_mat <- matrix(NA_real_, nrow = B, ncol = length(focuses))

for (b in seq_len(B)) {
  bvec <- betas[b, ]
  
  p1_b <- plogis(as.numeric(X_p1 %*% bvec))
  p0_b <- plogis(as.numeric(X_p0 %*% bvec))
  
  # sweep across focuses
  for (j in seq_along(focuses)) {
    rj <- rvec_list[[j]]
    p_exp <- rj * p1_b + (1 - rj) * p0_b
    props_mat[b, j] <- safe_wmean(p_exp, w)
  }
}

# 4) Summarize to mean and 95% CI for each focus
sens_ci <- tibble(
  focus = focuses,
  mean  = apply(props_mat, 2, mean),
  low   = apply(props_mat, 2, quantile, 0.025),
  high  = apply(props_mat, 2, quantile, 0.975)
) %>%
  mutate(across(c(mean, low, high), ~ . * 100))   # convert to %

# 5) Plot: ribbon + your original point estimate line/anchors
ggplot() +
  geom_ribbon(data = sens_ci,
              aes(x = focus, ymin = low, ymax = high),
              fill = "#1b9e77", alpha = 0.18) +
  geom_line(data = sens_rows,
            aes(x = focus, y = predicted_electrification_pct),
            linewidth = 1.2, color = "#1b9e77") +
  geom_point(data = subset(sens_rows, focus %in% c(0, 0.5, 1)),
             aes(x = focus, y = predicted_electrification_pct),
             size = 3, color = "#d95f02") +
  # Current Funding Pattern line/label
  { 
    obs_prop <- safe_wmean(dat_cc$pct_esb_fleet_prop, w)
    list(
      geom_hline(yintercept = 100 * obs_prop, linetype = "dashed"),
      annotate("text", x = 0.02, y = 100 * obs_prop,
               label = paste0("Current Funding Pattern: ",
                              round(100 * obs_prop, 2), "%"),
               hjust = 0, vjust = -0.6, size = 4)
    )
  } +
  scale_x_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("Rural-focused", "Balanced", "Urban-focused")
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Sensitivity of National Electrification to Policy Focus Mix",
    subtitle = "Line = point estimate; shaded band = 95% CI from GLM coefficient uncertainty",
    x = "Policy Focus Mix (0 = Rural, 1 = Urban)",
    y = "Predicted % of Fleet Electrified (fleet-weighted)"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

