# ===============================================================
# EPA Prioritization × Locale Analysis
# Fleet-weighted averages + bootstrap test + custom palette
# ===============================================================

library(tidyverse)
library(ggeffects)
library(wesanderson)

set.seed(123)  # for reproducibility

# ===============================================================
# --- Prepare Data ---
# ===============================================================
df_clean <- df_clean %>%
  mutate(
    prioritized_any = ifelse(epa2022_prioritized == 1 | epa2023_prioritized == 1, 1L, 0L),
    locale = droplevels(as.factor(locale)),
    total_buses = as.numeric(total_buses)
  )

# ===============================================================
# --- GLM: interaction of prioritization × locale ---
# ===============================================================
glm_any_locale <- glm(
  pct_esb_fleet_prop ~ prioritized_any * locale + pct_poverty + asthma18p,
  data = df_clean,
  family = quasibinomial(link = "logit")
)

# ===============================================================
# --- Predicted electrification by locale and prioritization ---
# ===============================================================
pred_locale_df <- as.data.frame(
  ggeffect(glm_any_locale, terms = c("locale", "prioritized_any"))
) %>%
  rename(locale = x, prioritized_any = group)

# ===============================================================
# --- Fleet-weighted national averages ---
# ===============================================================
national_avg <- df_clean %>%
  filter(!is.na(total_buses), total_buses > 0) %>%
  group_by(prioritized_any) %>%
  summarise(
    weighted_mean = weighted.mean(pct_esb_fleet_prop, total_buses, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(label = ifelse(prioritized_any == 1, "Prioritized Avg", "Not Prioritized Avg"))

# ===============================================================
# --- Weighted bootstrap test for national difference ---
# ===============================================================
boot_diff <- function(df, n_boot = 5000) {
  boot_est <- replicate(n_boot, {
    samp <- df[sample(nrow(df), replace = TRUE), ]
    wmeans <- samp %>%
      group_by(prioritized_any) %>%
      summarise(wm = weighted.mean(pct_esb_fleet_prop, total_buses, na.rm = TRUE),
                .groups = "drop")
    diff <- diff(wmeans$wm)  # prioritized(1) - not(0)
    return(diff)
  })
  tibble(
    diff_mean = mean(boot_est),
    diff_low = quantile(boot_est, 0.025),
    diff_high = quantile(boot_est, 0.975)
  )
}

national_test <- boot_diff(df_clean %>% filter(!is.na(total_buses), total_buses > 0)) %>%
  mutate(
    diff_pct_points = diff_mean * 100,
    diff_low_pp = diff_low * 100,
    diff_high_pp = diff_high * 100
  )

# ===============================================================
# --- Visualization: far-left labels + Wes Anderson palette ---
# ===============================================================

palette_darjeeling <- wes_palette("Darjeeling1", 2, type = "discrete")

# Determine an x-position just outside the first locale
leftmost_label_x <- -0.5

ggplot(pred_locale_df,
       aes(x = locale, y = predicted * 100,
           group = prioritized_any, color = as.factor(prioritized_any))) +
  geom_point(aes(shape = as.factor(prioritized_any)), size = 3) +
  geom_line(aes(linetype = as.factor(prioritized_any)), linewidth = 1) +
  geom_errorbar(aes(ymin = conf.low * 100, ymax = conf.high * 100),
                width = 0.2, alpha = 0.6) +
  # Fleet-weighted national average lines
  geom_hline(
    data = national_avg,
    aes(yintercept = weighted_mean * 100, color = as.factor(prioritized_any)),
    linetype = "dashed", linewidth = 0.9, show.legend = FALSE
  ) +
  # Labels moved FAR LEFT (outside plotting area)
  geom_text(
    data = national_avg,
    aes(
      x = leftmost_label_x,
      y = weighted_mean * 100,
      label = paste0(label, ": ", round(weighted_mean * 100, 2), "%"),
      color = as.factor(prioritized_any)
    ),
    hjust = 0, vjust = -0.4, size = 4, fontface = "bold", show.legend = FALSE
  ) +
  scale_color_manual(values = palette_darjeeling) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(
    title = "Predicted Electrification by Locale and EPA Prioritization",
    subtitle = paste0(
      "Model controls for poverty rate and asthma prevalence\n",
      "Dashed lines = fleet-weighted national averages; ",
      "bootstrap 95% CI for national difference: ",
      round(national_test$diff_pct_points, 2), " pp [",
      round(national_test$diff_low_pp, 2), ", ",
      round(national_test$diff_high_pp, 2), "]"
    ),
    x = "Locale",
    y = "Predicted % of Fleet Electrified",
    color = "EPA Prioritized",
    shape = "EPA Prioritized",
    linetype = "EPA Prioritized"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.margin = margin(10, 60, 10, 100)  # more room on left for labels
  ) +
  coord_cartesian(clip = "off")  # allow labels outside plot

