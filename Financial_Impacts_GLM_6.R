# ===============================================================
# Electric School Bus — Final Modeling + Policy Simulation
# Publication-ready version with clean, readable charts
# ===============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(wesanderson)
  library(MASS)        # mvrnorm
  library(ggeffects)
  library(ggrepel)
  if (requireNamespace("conflicted", quietly = TRUE)) {
    conflicted::conflict_prefer("select", "dplyr", quiet = TRUE)
    conflicted::conflict_prefer("filter", "dplyr", quiet = TRUE)
    conflicted::conflict_prefer("mutate", "dplyr", quiet = TRUE)
    conflicted::conflict_prefer("rename", "dplyr", quiet = TRUE)
    conflicted::conflict_prefer("summarise", "dplyr", quiet = TRUE)
  }
})

# ------------------------ HELPERS ------------------------------
parse_percent <- function(x) {
  vals <- readr::parse_number(as.character(x))
  if (is.finite(max(vals, na.rm = TRUE)) && max(vals, na.rm = TRUE) > 1.5)
    vals <- vals / 100
  vals
}
safe_weighted_mean <- function(x, w) {
  idx <- is.finite(x) & is.finite(w) & !is.na(x) & !is.na(w) & w > 0
  if (!any(idx)) return(NA_real_)
  sum(x[idx] * w[idx]) / sum(w[idx])
}

# --------------------- LOAD & CLEAN ----------------------------
df <- read_csv("esb_df_clean.csv")

df_clean <- df %>%
  mutate(
    epa2022_prioritized = ifelse(epa2022_prioritized %in% c("Yes","Y",1,TRUE), 1L, 0L),
    epa2023_prioritized = ifelse(epa2023_prioritized %in% c("Yes","Y",1,TRUE), 1L, 0L),
    prioritized_any = ifelse(epa2022_prioritized == 1 | epa2023_prioritized == 1, 1L, 0L),
    pct_esb_fleet_prop = parse_percent(pct_esb_fleet),
    pct_poverty = as.numeric(pct_poverty),
    asthma18p = as.numeric(asthma18p),
    total_buses = as.numeric(total_buses),
    locale = as.character(locale)
  ) %>%
  filter(!is.na(locale), locale != "Unknown") %>%
  mutate(locale = droplevels(as.factor(locale))) %>%
  filter(!is.na(pct_esb_fleet_prop),
         pct_esb_fleet_prop >= 0, pct_esb_fleet_prop <= 1)

# Complete cases only
vars_needed <- c("prioritized_any","locale","pct_poverty","asthma18p","total_buses")
keep <- complete.cases(df_clean[, vars_needed])
df_cc <- df_clean[keep, , drop = FALSE] %>%
  mutate(locale = droplevels(locale))

message("Rows kept for modeling: ", nrow(df_cc))

# ------------------------ GLM MODEL ----------------------------
glm_any_locale <- glm(
  pct_esb_fleet_prop ~ prioritized_any * locale + pct_poverty + asthma18p,
  data = df_cc,
  family = quasibinomial(link = "logit")
)
summary(glm_any_locale)

# ===============================================================
# PLOT 1 — Predicted Electrification by Locale × Prioritization
# ===============================================================
pred_locale_df <- as.data.frame(
  ggeffect(glm_any_locale, terms = c("locale", "prioritized_any"))
) %>% rename(locale = x, prioritized_any = group)

national_avg <- df_cc %>%
  group_by(prioritized_any) %>%
  summarise(weighted_mean = safe_weighted_mean(pct_esb_fleet_prop, total_buses),
            .groups = "drop") %>%
  mutate(label = ifelse(prioritized_any == 1, "Prioritized Avg", "Not Prioritized Avg"))

palette_darjeeling <- wes_palette("Darjeeling1", 2, type = "discrete")

ggplot(pred_locale_df,
       aes(x = locale, y = predicted * 100,
           group = as.factor(prioritized_any),
           color = as.factor(prioritized_any))) +
  geom_point(size = 3) +
  geom_line(linewidth = 1.1,
            aes(linetype = as.factor(prioritized_any))) +
  geom_errorbar(aes(ymin = conf.low * 100, ymax = conf.high * 100),
                width = 0.2, alpha = 0.6) +
  geom_hline(
    data = national_avg,
    aes(yintercept = weighted_mean * 100, color = as.factor(prioritized_any)),
    linetype = "dashed", linewidth = 0.9, show.legend = FALSE
  ) +
  geom_text(
    data = subset(national_avg, prioritized_any == 1),
    aes(label = paste0(label, ": ", round(weighted_mean * 100, 2), "%"),
        y = weighted_mean * 100, color = as.factor(prioritized_any)),
    x = 0.6, hjust = 0, vjust = -0.6, size = 4.2, fontface = "bold"
  ) +
  geom_text(
    data = subset(national_avg, prioritized_any == 0),
    aes(label = paste0(label, ": ", round(weighted_mean * 100, 2), "%"),
        y = weighted_mean * 100, color = as.factor(prioritized_any)),
    x = 4.3, hjust = 0, vjust = -0.6, size = 4.2, fontface = "bold"
  ) +
  scale_color_manual(values = palette_darjeeling,
                     labels = c("Not Prioritized", "Prioritized")) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Predicted Electrification by Locale and EPA Prioritization",
    subtitle = "Controls: poverty rate & asthma prevalence\nDashed lines = fleet-weighted national averages",
    x = "Locale", y = "Predicted % of Fleet Electrified",
    color = "EPA Status", linetype = "EPA Status"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.key.width = unit(1.2, "cm"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

# ===============================================================
# PARAMETRIC SCENARIO SIMULATION
# ===============================================================
B_param <- 400
set.seed(123)

beta_hat <- coef(glm_any_locale)
V_hat <- vcov(glm_any_locale)
betas <- MASS::mvrnorm(n = B_param, mu = beta_hat, Sigma = V_hat)

form <- formula(glm_any_locale)
X_p1 <- model.matrix(form, transform(df_cc, prioritized_any = 1L))
X_p0 <- model.matrix(form, transform(df_cc, prioritized_any = 0L))
w <- as.numeric(df_cc$total_buses)

overall_prior_rate <- mean(df_cc$prioritized_any, na.rm = TRUE)
locale_sizes <- df_cc %>% count(locale, name = "n")

rescale_policy_to_overall <- function(policy_tbl, locale_sizes, target_overall) {
  tmp <- policy_tbl %>% left_join(locale_sizes, by = "locale")
  k0 <- target_overall / (sum(tmp$target_rate * tmp$n) / sum(tmp$n))
  policy_rescaled <- tmp %>% mutate(target_rate_adj = pmin(pmax(target_rate * k0, 0), 1))
  policy_rescaled %>% select(locale, target_rate_adj)
}

policy_rural <- tibble(
  locale = levels(df_cc$locale),
  target_rate = case_when(
    locale == "Rural" ~ 0.60,
    locale == "Town" ~ 0.40,
    locale == "Suburban" ~ 0.30,
    locale == "Urban" ~ 0.20,
    TRUE ~ 0.40
  )
)
policy_balanced <- tibble(locale = levels(df_cc$locale), target_rate = 0.40)
policy_urban <- tibble(
  locale = levels(df_cc$locale),
  target_rate = case_when(
    locale == "Rural" ~ 0.20,
    locale == "Town" ~ 0.30,
    locale == "Suburban" ~ 0.40,
    locale == "Urban" ~ 0.60,
    TRUE ~ 0.40
  )
)

policy_rural_adj    <- rescale_policy_to_overall(policy_rural,    locale_sizes, overall_prior_rate)
policy_balanced_adj <- rescale_policy_to_overall(policy_balanced, locale_sizes, overall_prior_rate)
policy_urban_adj    <- rescale_policy_to_overall(policy_urban,    locale_sizes, overall_prior_rate)

r_current  <- df_cc$prioritized_any
r_rural    <- df_cc %>% left_join(policy_rural_adj,    by = "locale") %>% pull(target_rate_adj)
r_balanced <- df_cc %>% left_join(policy_balanced_adj, by = "locale") %>% pull(target_rate_adj)
r_urban    <- df_cc %>% left_join(policy_urban_adj,    by = "locale") %>% pull(target_rate_adj)

mw_curr <- mw_rur <- mw_bal <- mw_urb <- numeric(B_param)

for (b in seq_len(B_param)) {
  bvec <- betas[b, ]
  p1 <- plogis(as.numeric(X_p1 %*% bvec))
  p0 <- plogis(as.numeric(X_p0 %*% bvec))
  p_curr <- r_current * p1 + (1 - r_current) * p0
  p_rur  <- r_rural * p1 + (1 - r_rural) * p0
  p_bal  <- r_balanced * p1 + (1 - r_balanced) * p0
  p_urb  <- r_urban * p1 + (1 - r_urban) * p0
  mw_curr[b] <- safe_weighted_mean(p_curr, w)
  mw_rur[b]  <- safe_weighted_mean(p_rur,  w)
  mw_bal[b]  <- safe_weighted_mean(p_bal,  w)
  mw_urb[b]  <- safe_weighted_mean(p_urb,  w)
}

# ===============================================================
# PLOT 2 — Scenario Composition (100% stacked)
# ===============================================================
locale_sizes <- df_cc %>% count(locale, name = "n")

scenario_alloc <- function(policy_adj, scen_name) {
  locale_sizes %>%
    left_join(policy_adj, by = "locale") %>%
    mutate(expected_prioritized = n * target_rate_adj) %>%
    summarise(locale, share = expected_prioritized / sum(expected_prioritized),
              .groups = "drop") %>%
    mutate(scenario = scen_name)
}
current_alloc <- df_cc %>%
  group_by(locale) %>%
  summarise(prioritized = sum(prioritized_any, na.rm = TRUE), .groups = "drop") %>%
  mutate(share = prioritized / sum(prioritized),
         scenario = "Current Funding Pattern") %>%
  select(locale, share, scenario)

shares_df <- bind_rows(
  current_alloc,
  scenario_alloc(policy_rural_adj, "Rural-focused"),
  scenario_alloc(policy_balanced_adj, "Balanced"),
  scenario_alloc(policy_urban_adj, "Urban-focused")
) %>%
  mutate(
    locale = factor(locale, levels = levels(df_cc$locale)),
    scenario = factor(scenario,
                      levels = c("Current Funding Pattern", "Rural-focused", "Balanced", "Urban-focused")),
    pct_lab = scales::percent(share, accuracy = 0.1)
  )

palette_locales <- c("Rural"="#e41a1c", "Suburban"="#1b9e77", "Town"="#ffb000", "Urban"="#d95f02")
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
# LIFETIME SAVINGS (PARAMETRIC) + Δ VS CURRENT
# ===============================================================
lifetime_savings_per_bus <- 100000
total_buses_total <- sum(df_cc$total_buses, na.rm = TRUE)
factor_dollars <- total_buses_total * lifetime_savings_per_bus

sav_curr <- mw_curr * factor_dollars
sav_rur  <- mw_rur  * factor_dollars
sav_bal  <- mw_bal  * factor_dollars
sav_urb  <- mw_urb  * factor_dollars

life_tbl <- function(x) tibble(mean = mean(x), low = quantile(x, .025), high = quantile(x, .975))
life_current <- life_tbl(sav_curr)
life_rural   <- life_tbl(sav_rur)
life_bal     <- life_tbl(sav_bal)
life_urban   <- life_tbl(sav_urb)

lifetime_results <- tribble(
  ~scenario, ~mean, ~low, ~high,
  "Current Funding Pattern", life_current$mean, life_current$low, life_current$high,
  "Rural-focused", life_rural$mean, life_rural$low, life_rural$high,
  "Balanced", life_bal$mean, life_bal$low, life_bal$high,
  "Urban-focused", life_urban$mean, life_urban$low, life_urban$high
)

# --- FIX: Proper factor order so "Current Funding Pattern" shows on the left
lifetime_results <- lifetime_results %>%
  mutate(
    scenario = forcats::fct_relevel(scenario,
                                    "Current Funding Pattern",
                                    "Rural-focused",
                                    "Balanced",
                                    "Urban-focused")
  )

# Helper for billions formatting
dollar_billions <- function(x) paste0("$", formatC(x / 1e9, format = "f", digits = 3), "B")

# Consistent fill colors
fill_colors <- c("Current Funding Pattern" = "#0072B2",
                 "Rural-focused" = "#e41a1c",
                 "Balanced" = "#1b9e77",
                 "Urban-focused" = "#ffb000")

# ===============================================================
# PLOT 1 — Total Lifetime Savings by Policy Scenario
# ===============================================================
ggplot(lifetime_results, aes(x = scenario, y = mean, fill = scenario)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = low, ymax = high),
                width = 0.15, linewidth = 0.8, color = "black") +
  geom_text(aes(label = dollar_billions(mean)),
            position = position_nudge(x = 0.25),
            vjust = 0.5, size = 4.5, fontface = "bold", color = "black") +
  scale_fill_manual(values = c(
    "Current Funding Pattern" = "#0072B2",  # Blue for current funding
    "Rural-focused" = "#e41a1c",
    "Balanced" = "#1b9e77",
    "Urban-focused" = "#ffb000"
  )) +
  scale_y_continuous(labels = dollar_billions,
                     expand = expansion(mult = c(0.05, 0.12))) +
  labs(
    title = "Total Lifetime Savings by Policy Scenario",
    subtitle = "Assumes $100,000 lifetime savings per ESB; Parametric 95% CIs",
    x = NULL, y = "Total Lifetime Savings (USD, billions)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)  # Center subtitle
  )


# ===============================================================
# PLOT 2 — Δ Lifetime Savings vs Current Policy
# ===============================================================

# Reorder to place Balanced in the middle
delta_results <- delta_results %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c("Rural-focused", "Balanced", "Urban-focused")
    )
  )

ggplot(delta_results, aes(x = scenario, y = mean, fill = scenario)) +
  geom_hline(yintercept = 0, color = "gray60") +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_errorbar(aes(ymin = low, ymax = high),
                width = 0.15, linewidth = 0.8, color = "black") +
  geom_text(aes(label = ifelse(mean > 0,
                               paste0("+", dollar_billions(mean)),
                               dollar_billions(mean))),
            position = position_nudge(x = 0.25),
            vjust = 0.5, size = 4.3, fontface = "bold", color = "black") +
  scale_fill_manual(values = c(
    "Rural-focused" = "#e41a1c",
    "Balanced" = "#1b9e77",
    "Urban-focused" = "#ffb000"
  )) +
  scale_y_continuous(labels = dollar_billions,
                     expand = expansion(mult = c(0.06, 0.18))) +
  labs(
    title = "Δ Lifetime Savings vs Current Policy",
    subtitle = "Positive values indicate additional lifetime savings relative to the current funding pattern",
    x = NULL, y = "Change in Total Lifetime Savings (USD, billions)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)  # Center subtitle
  )

