# ===============================================================
# EPA Prioritization × Locale Analysis
# ===============================================================

library(tidyverse)
library(ggeffects)

# --- Ensure combined prioritization variable exists ---
df_clean <- df_clean %>%
  mutate(
    prioritized_any = ifelse(epa2022_prioritized == 1 | epa2023_prioritized == 1, 1L, 0L),
    locale = droplevels(as.factor(locale))
  )

# ===============================================================
# --- Fit GLM with interaction between prioritization and locale ---
# ===============================================================

glm_any_locale <- glm(
  pct_esb_fleet_prop ~ prioritized_any * locale + pct_poverty + asthma18p,
  data = df_clean,
  family = quasibinomial(link = "logit")   # robust for proportion outcomes
)

summary(glm_any_locale)

# ===============================================================
# --- Predicted electrification by locale and prioritization ---
# ===============================================================

# Generate predicted probabilities for each locale × prioritization combo
pred_locale_df <- as.data.frame(
  ggeffect(glm_any_locale, terms = c("locale", "prioritized_any"))
)

# Rename columns for clarity
pred_locale_df <- pred_locale_df %>%
  rename(
    locale = x,
    prioritized_any = group
  )

# ===============================================================
# --- Summarize differences (Prioritized vs. Not) by locale ---
# ===============================================================

effect_by_locale <- pred_locale_df %>%
  select(locale, prioritized_any, predicted) %>%
  mutate(prioritized_any = as.numeric(as.character(prioritized_any))) %>%
  pivot_wider(
    names_from = prioritized_any,
    values_from = predicted,
    names_prefix = "prioritized_"
  ) %>%
  mutate(
    difference = prioritized_1 - prioritized_0,
    pct_nonprioritized = prioritized_0 * 100,
    pct_prioritized = prioritized_1 * 100,
    diff_pct_points = difference * 100
  ) %>%
  arrange(desc(diff_pct_points))

# View clean results table
effect_by_locale

# ===============================================================
# --- Compute weighted national averages for prioritized vs. not ---
# ===============================================================

national_avg <- df_clean %>%
  group_by(prioritized_any) %>%
  summarise(weighted_mean = mean(pct_esb_fleet_prop, na.rm = TRUE)) %>%
  mutate(label = ifelse(prioritized_any == 1, "Prioritized Avg", "Not Prioritized Avg"))

# ===============================================================
# --- Visualization ---
# ===============================================================

ggplot(pred_locale_df,
       aes(x = locale, y = predicted * 100,
           group = prioritized_any, color = as.factor(prioritized_any))) +
  geom_point(aes(shape = as.factor(prioritized_any)), size = 3) +
  geom_line(aes(linetype = as.factor(prioritized_any)), linewidth = 1) +
  geom_errorbar(aes(ymin = conf.low * 100, ymax = conf.high * 100),
                width = 0.2, alpha = 0.6) +
  # Add national average line
  geom_hline(data = national_avg,
             aes(yintercept = weighted_mean * 100, color = as.factor(prioritized_any)),
             linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
  geom_text(
    data = national_avg,
    aes(x = Inf, y = weighted_mean * 100,
        label = paste0(label, ": ", round(weighted_mean * 100, 2), "%"),
        color = as.factor(prioritized_any)),
    hjust = 1.1, vjust = -0.5, size = 4, fontface = "bold", show.legend = FALSE
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Predicted Electrification by Locale and EPA Prioritization",
    subtitle = "Model controls for poverty rate and asthma prevalence\nDashed lines show national averages by prioritization status",
    x = "Locale",
    y = "Predicted % of Fleet Electrified",
    color = "EPA Prioritized",
    shape = "EPA Prioritized",
    linetype = "EPA Prioritized"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

