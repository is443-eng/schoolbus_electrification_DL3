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




