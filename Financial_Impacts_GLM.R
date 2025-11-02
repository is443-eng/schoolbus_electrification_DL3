library(tidyverse)

df <- read_csv("esb_df_clean.csv")

#xxxxxxxxxxxxxxxxxxxx Making the model xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
# --- robust parse for the percentage column ---
parse_percent <- function(x) {
  # strip non-numeric characters (like % or commas) then numeric
  vals <- readr::parse_number(as.character(x))
  # if it looks like 0–100 scale, convert to 0–1
  if (is.finite(max(vals, na.rm = TRUE)) && max(vals, na.rm = TRUE) > 1.5) {
    vals <- vals / 100
  }
  vals
}

df_clean <- df %>%
  mutate(
    # EPA flags to 0/1
    epa2022_prioritized = ifelse(epa2022_prioritized %in% c("Yes","Y",1,TRUE), 1L, 0L),
    epa2023_prioritized = ifelse(epa2023_prioritized %in% c("Yes","Y",1,TRUE), 1L, 0L),
    prioritized_any = ifelse(epa2022_prioritized == 1 | epa2023_prioritized == 1, 1, 0),
    # parse and scale percentage to proportion
    pct_esb_fleet_prop = parse_percent(pct_esb_fleet)
  ) %>%
  # keep only valid proportions for binomial GLM
  filter(!is.na(pct_esb_fleet_prop),
         pct_esb_fleet_prop >= 0,
         pct_esb_fleet_prop <= 1)

# quick sanity checks
summary(df_clean$pct_esb_fleet_prop)
table(df_clean$epa2022_prioritized, useNA = "ifany")
table(df_clean$epa2023_prioritized, useNA = "ifany")

# --- GLM: binomial with logit link ---
glm_model <- glm(
  pct_esb_fleet_prop ~ epa2023_prioritized + epa2022_prioritized,
  data = df_clean,
  family = binomial(link = "logit")
)

summary(glm_model)

# Odds ratios for easier reading
exp(coef(glm_model))

#xxxxxxxxxxxxxxxxxxxxxxxx TRYING TO GET QUANTITATIVE xxxxxxxxxxxxxxxxxxxxxxxxx
#how much increase is expected in the fleet electrification % with EPA prior.
newdata <- expand.grid(
  epa2022_prioritized = c(0, 1),
  epa2023_prioritized = c(0, 1)
)

newdata$predicted_prob <- predict(glm_model, newdata, type = "response")
newdata


#xxxxxxxxxxxxxxxxxxxxxxxxx PULL THE LEVER xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#xxxxxxxxxxxxxxxxxxxxxxxxx percentage of districts prioritized xxxxxxxxxxxxxx
# Count and percentage
prioritization_summary <- df_clean %>%
  summarize(
    total_districts = n(),
    pct_2022_prioritized = mean(epa2022_prioritized, na.rm = TRUE) * 100,
    pct_2023_prioritized = mean(epa2023_prioritized, na.rm = TRUE) * 100,
    pct_any_prioritized = mean(prioritized_any, na.rm = TRUE) * 100
  )

prioritization_summary


# ===============================================================
# Policy Shift Simulation — add after your current visualization
# Reuses: df_clean, glm_any_locale (quasibinomial), total_buses
# ===============================================================

library(tidyverse)
library(wesanderson)
set.seed(42)

# ---------- Helpers ----------
# 1) Status-quo overall prioritization rate
overall_prior_rate <- mean(df_clean$prioritized_any, na.rm = TRUE)

# 2) Count by locale, for rescaling
locale_sizes <- df_clean %>% count(locale, name = "n")

# 3) Rescale target rates so the *overall* prioritization rate matches status quo
rescale_policy_to_overall <- function(policy_tbl, locale_sizes, target_overall) {
  # Merge sizes to compute weighted average
  tmp <- policy_tbl %>%
    left_join(locale_sizes, by = "locale")
  
  # function to compute overall given multiplier k
  overall_given_k <- function(k) {
    rates <- pmin(pmax(tmp$target_rate * k, 0), 1)  # clip to [0,1]
    sum(rates * tmp$n) / sum(tmp$n)
  }
  
  # quick solve for k by ratio, then refine and clip
  k0 <- target_overall / (sum(tmp$target_rate * tmp$n) / sum(tmp$n))
  # small safeguard if denominator is ~0
  if (!is.finite(k0) || k0 <= 0) k0 <- 1
  
  # apply and clip
  policy_rescaled <- tmp %>%
    mutate(target_rate_adj = pmin(pmax(target_rate * k0, 0), 1)) %>%
    select(locale, target_rate_adj)
  
  # final small correction if clipping changed the overall rate notably
  final_overall <- sum(policy_rescaled$target_rate_adj * tmp$n) / sum(tmp$n)
  policy_rescaled$target_rate_adj <- policy_rescaled$target_rate_adj * (target_overall / final_overall)
  
  policy_rescaled %>%
    mutate(target_rate_adj = pmin(pmax(target_rate_adj, 0), 1))
}

# 4) Simulate a scenario and return fleet-weighted national electrification
simulate_policy <- function(df, policy_adj, model) {
  df_sim <- df %>%
    left_join(policy_adj, by = "locale") %>%
    mutate(
      # random assignment within locale to match target_rate_adj
      prioritized_any = rbinom(n(), 1, p = target_rate_adj),
      pred = predict(model, newdata = cur_data_all(), type = "response")
    )
  
  tibble(
    weighted_electrification = weighted.mean(df_sim$pred, df_sim$total_buses, na.rm = TRUE)
  )
}

# ---------- Define scenarios (pre-rescale) ----------
policy_rural <- tibble(
  locale = levels(df_clean$locale),
  target_rate = case_when(
    locale == "Rural"    ~ 0.60,
    locale == "Town"     ~ 0.40,
    locale == "Suburban" ~ 0.30,
    locale == "Urban"    ~ 0.20,
    TRUE ~ 0.30
  )
)

policy_balanced <- tibble(
  locale = levels(df_clean$locale),
  target_rate = 0.40  # same target in every locale
)

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

# ---------- Rescale to keep overall rate equal to status quo ----------
policy_rural_adj   <- rescale_policy_to_overall(policy_rural,   locale_sizes, overall_prior_rate)
policy_balanced_adj<- rescale_policy_to_overall(policy_balanced,locale_sizes, overall_prior_rate)
policy_urban_adj   <- rescale_policy_to_overall(policy_urban,   locale_sizes, overall_prior_rate)

# For transparency, show the final locale targets used
target_tables <- bind_rows(
  policy_rural_adj  %>% mutate(scenario = "Rural-focused"),
  policy_balanced_adj %>% mutate(scenario = "Balanced"),
  policy_urban_adj  %>% mutate(scenario = "Urban-focused")
) %>%
  relocate(scenario) %>%
  arrange(scenario, locale) %>%
  mutate(target_rate_pct = round(100 * target_rate_adj, 1))

target_tables  # <- inspect these to see the actual per-locale target %s

# ---------- Compute status-quo (using current prioritized_any) ----------
df_statusquo <- df_clean %>%
  mutate(pred = predict(glm_any_locale, newdata = cur_data_all(), type = "response"))

statusquo_weighted <- weighted.mean(df_statusquo$pred, df_statusquo$total_buses, na.rm = TRUE)

# ---------- Simulate scenarios ----------
sim_rural   <- simulate_policy(df_clean, policy_rural_adj,    glm_any_locale)
sim_bal     <- simulate_policy(df_clean, policy_balanced_adj, glm_any_locale)
sim_urban   <- simulate_policy(df_clean, policy_urban_adj,    glm_any_locale)

scenario_results <- tibble(
  scenario = c("Status Quo", "Rural-focused", "Balanced", "Urban-focused"),
  weighted_electrification = c(statusquo_weighted, sim_rural$weighted_electrification,
                               sim_bal$weighted_electrification, sim_urban$weighted_electrification)
) %>%
  mutate(weighted_electrification_pct = 100 * weighted_electrification)

scenario_results

# ---------- Plot comparison ----------
pal <- wes_palette("Darjeeling1", 4, type = "discrete")

ggplot(scenario_results, aes(x = scenario, y = weighted_electrification_pct, fill = scenario)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(round(weighted_electrification_pct, 2), "%")),
            vjust = -0.6, fontface = "bold") +
  scale_fill_manual(values = pal) +
  scale_y_continuous(labels = scales::percent_format(scale = 1), expand = expansion(mult = c(0.02, 0.10))) +
  labs(
    title = "National Fleet Electrification Under Alternative Prioritization Focus",
    subtitle = paste0("All scenarios rescaled to maintain overall prioritization rate of ",
                      round(100 * overall_prior_rate, 1), "%"),
    x = NULL, y = "Predicted % of Fleet Electrified (Fleet-weighted)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

