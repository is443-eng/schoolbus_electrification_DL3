# ===============================================================
# Electric School Bus — Modeling + Policy Simulation w/ 95% CIs
# ===============================================================

library(tidyverse)
library(ggeffects)
library(wesanderson)

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
# Status Quo (Predicted under current prioritization)
boot_statusquo_pred <- function(model, df, n_boot = 2000) {
  est <- replicate(n_boot, {
    samp <- df[sample(nrow(df), replace = TRUE), ]
    preds <- predict(model, newdata = samp, type = "response")
    safe_weighted_mean(preds, samp$total_buses)
  })
  tibble(mean = mean(est), low = quantile(est, .025), high = quantile(est, .975))
}

# Empirical Current Funding Pattern (observed, no model)
boot_current_observed <- function(df, n_boot = 2000) {
  est <- replicate(n_boot, {
    samp <- df[sample(nrow(df), replace = TRUE), ]
    safe_weighted_mean(samp$pct_esb_fleet_prop, samp$total_buses)
  })
  tibble(mean = mean(est), low = quantile(est, .025), high = quantile(est, .975))
}

# Policy scenario predicted (random assignment per target_rate_adj within locale)
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
  geom_text(aes(label = sprintf("%.2f%%", weighted_electrification_pct)),
            vjust = -0.6, fontface = "bold") +
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



