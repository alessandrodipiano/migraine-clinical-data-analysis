folder <- "C:/Users/aless/Desktop/medical applications/data/imputed/dropout_not_imputed_longitudinal"

full_list <- readRDS(
  file.path(folder, "full_longitudinal_all_imputations.rds")
)

length(full_list)

anyNA(full_list[[1]])








library(dplyr)
library(mice)
library(lme4)
library(broom.mixed)




calc_diffs <- function(df) {
  df %>%
    group_by(SUBJECT_ID, CYCLE) %>%
    summarise(
      mmd_m1   = MMDs[MONTH == 1][1],
      mmd_m12  = MMDs[MONTH == 12][1],
      diff_mmd = mmd_m1 - mmd_m12,
      .groups = "drop"
    ) %>%
    filter(!is.na(diff_mmd))
}

# ---- 3. The "List-to-Pool" Workflow
# Since with(imp_long, ...) fails on a list, we use lapply + as.mira

# A. Model: Overall Mean Reduction
fit_overall_list <- lapply(full_list, function(df) {
  d <- calc_diffs(df)
  lm(diff_mmd ~ 1, data = d)
})
pooled_overall <- pool(as.mira(fit_overall_list))

# B. Model: Per-Cycle Contrast (Testing if Cycle X differs from Cycle 1)
fit_cycle_test_list <- lapply(full_list, function(df) {
  d <- calc_diffs(df)
  lm(diff_mmd ~ factor(CYCLE), data = d)
})
pooled_cycle_test <- pool(as.mira(fit_cycle_test_list))

# C. Model: Baseline Trend (Mixed Effects)
fit_lme_list <- lapply(full_list, function(df) {
  df_base <- df %>% filter(MONTH == 1)
  # Using lmer for longitudinal subjects
  lmer(MMDs ~ as.numeric(CYCLE) + (1 | SUBJECT_ID), data = df_base)
})
pooled_lme <- pool(as.mira(fit_lme_list))

# ---- 4. Results Extraction
summary(pooled_overall, conf.int = TRUE)
summary(pooled_cycle_test, conf.int = TRUE)
summary(pooled_lme, conf.int = TRUE)



library(dplyr)

cycle_response <- function(df) {
  df %>%
    filter(MONTH %in% c(1, 12)) %>%
    group_by(SUBJECT_ID, CYCLE) %>%
    summarise(
      mmd_m1  = MMDs[MONTH == 1][1],
      mmd_m12 = MMDs[MONTH == 12][1],
      .groups = "drop"
    ) %>%
    filter(!is.na(mmd_m1) & !is.na(mmd_m12)) %>%
    mutate(
      perc_reduction = (mmd_m1 - mmd_m12) / mmd_m1,
      resp_30 = perc_reduction >= 0.30,
      resp_50 = perc_reduction >= 0.50
    )
}


# 1) Compute responder summaries within each imputation
resp_by_imp <- lapply(full_list, function(dfj) {
  d <- cycle_response(dfj)
  
  d %>%
    group_by(CYCLE) %>%
    summarise(
      prop_30 = mean(resp_30),
      prop_50 = mean(resp_50),
      n = n(),
      .groups = "drop"
    )
})

# 2) Pool cycle-specific proportions across imputations (robust pooling)
resp_all <- bind_rows(resp_by_imp, .id = "imp")   # keep imp as character; no coercion needed

resp_pooled_fix <- resp_all %>%
  group_by(CYCLE) %>%
  summarise(
    prop_30 = mean(prop_30, na.rm = TRUE),
    prop_50 = mean(prop_50, na.rm = TRUE),
    n = round(mean(n, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(CYCLE)

resp_pooled_fix

# 3) Full treatment period: responder if achieved threshold in ANY cycle
full_period_resp <- lapply(full_list, function(dfj) {
  d <- cycle_response(dfj)
  
  d %>%
    group_by(SUBJECT_ID) %>%
    summarise(
      resp_30_any = any(resp_30),
      resp_50_any = any(resp_50),
      .groups = "drop"
    ) %>%
    summarise(
      prop_30_any = mean(resp_30_any),
      prop_50_any = mean(resp_50_any)
    )
})

full_pooled <- bind_rows(full_period_resp) %>%
  summarise(
    prop_30_any = mean(prop_30_any, na.rm = TRUE),
    prop_50_any = mean(prop_50_any, na.rm = TRUE)
  )

full_pooled

# 4) Stability check: number of evaluable (SUBJECT_ID, CYCLE) pairs per imputation
sapply(full_list, function(dfj) nrow(cycle_response(dfj))) |> summary()




resp_table <- resp_pooled_fix %>%
  mutate(
    resp30_pct = round(100 * prop_30, 1),
    resp50_pct = round(100 * prop_50, 1)
  ) %>%
  select(CYCLE, n, resp30_pct, resp50_pct)

resp_table






