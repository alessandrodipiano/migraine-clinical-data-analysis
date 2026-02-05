library(dplyr)

# ---- Load baseline (adjust path)
df_base <- read.csv(
  "C:/Users/aless/Desktop/medical applications/data/cleaned/df_baseline_clean.csv",
  header = TRUE
)

View(df_base)


# ---- Coerce key columns
df_base$SUBJECT_ID <- as.integer(df_base$SUBJECT_ID)

# IMPORTANT: adjust these names if your baseline uses different spellings
# Expecting: Suspension (0=completed, 1=discontinued), MONTHS_OF_TREAT
df_base$Suspension      <- as.integer(df_base$Suspension)
df_base$MONTHS_OF_TREAT <- as.integer(df_base$MONTHS_OF_TREAT)

df_long <- read.csv("C:/Users/aless/Desktop/medical applications/data/cleaned/df_longitudinal_clean.csv", header=TRUE)
View(df_long)







# ---- Join baseline onto longitudinal + define absolute treatment time
df_long <- df_long %>%
  left_join(df_base %>% select(SUBJECT_ID, Suspension, MONTHS_OF_TREAT),
            by = "SUBJECT_ID") %>%
  mutate(
    t_abs = 12 * (CYCLE - 1) + MONTH
  ) %>%
  mutate(
    # RULE:
    # - completers (Suspension==0): all scheduled times are eligible (in_window=TRUE)
    # - non-completers (Suspension==1): eligible only up to MONTHS_OF_TREAT
    in_window = case_when(
      Suspension == 0 ~ TRUE,
      Suspension == 1 ~ (!is.na(MONTHS_OF_TREAT) & t_abs <= MONTHS_OF_TREAT),
      TRUE            ~ TRUE  # if Suspension missing/unknown, default to not truncating
    )
  )

# Optional sanity checks
table(df_long$Suspension, useNA = "ifany")
summary(df_long$MONTHS_OF_TREAT)



df_long_work <- df_long %>%
  filter(in_window) %>%
  select(-Suspension, -MONTHS_OF_TREAT, -t_abs, -in_window)  # keep clean for imputation


library(lme4)
# Model 1: Fixed Trend (Your current setup)
fit1 <- lmer(MMDs ~ CYCLE + (1 | SUBJECT_ID), data = df_long_work)
# Model 2: Random Trend (Proposed setup)
fit2 <- lmer(MMDs ~ CYCLE + (CYCLE | SUBJECT_ID), data = df_long_work)

anova(fit1, fit2) # Lower AIC/BIC indicates the better method
df_long$SUBJECT_ID <- as.integer(df_long$SUBJECT_ID)
df_long$CYCLE <- as.integer(df_long$CYCLE)
df_long$MONTH <- as.integer(df_long$MONTH)



library(mice)
library(miceadds)

# Drop derived variable (correct)
if ("RED_MMD_VST01" %in% names(df_long_work)) {
  df_long_work <- subset(df_long_work, select = -RED_MMD_VST01)
}

cluster_var <- "SUBJECT_ID"
time_vars   <- c("CYCLE","MONTH")

cont_vars <- c("MMDs","DOSE","GGFAR","HADSA","HADSD","HIT6","INT","MIDAS")
cont_vars <- intersect(cont_vars, names(df_long_work))

# Ensure numeric
for (v in cont_vars) df_long_work[[v]] <- as.numeric(df_long_work[[v]])

# Init
ini  <- mice(df_long_work, maxit=0, printFlag=FALSE)
meth <- ini$method
pred <- ini$predictorMatrix

# Methods
meth[cont_vars] <- "2l.pmm"
meth[c(cluster_var, time_vars)] <- ""  # do not impute ID/time

# Predictor matrix
pred[,] <- 0
pred[, cluster_var] <- -2              # random intercept per subject
pred[cont_vars, time_vars] <- 1       # time predicts outcomes
pred[cont_vars, cont_vars] <- 1        # outcomes predict each other
diag(pred) <- 0
pred[time_vars, ] <- 0                 # time vars not targets

# Sanity check
to_impute <- names(meth)[meth != ""]
bad <- to_impute[rowSums(pred[to_impute, , drop=FALSE]) == 0]
print(bad)
stopifnot(length(bad) == 0)

# Run
imp_long <- mice(
  df_long_work,
  method = meth,
  predictorMatrix = pred,
  m = 20, maxit = 20, seed = 123,
  printFlag = TRUE
)


library(tidyr)
library(zoo)

df_locf <- df_long_work %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID, CYCLE) %>%
  fill(all_of(cont_vars), .direction = "down") %>%
  ungroup()

df_nocb <- df_long_work %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID, CYCLE) %>%
  fill(all_of(cont_vars), .direction = "up") %>%
  ungroup()

interp_vars <- intersect(c("MMDs","HIT6","MIDAS"), names(df_long_work))

df_interp <- df_long_work %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID, CYCLE) %>%
  mutate(across(
    all_of(interp_vars),
    ~ zoo::na.approx(.x, x = MONTH, na.rm = FALSE)
  )) %>%
  ungroup()
# After building df_long_work, no subject should have any rows with t_abs > MONTHS_OF_TREAT if non-completer.
check_post_dropout <- df_long %>%
  filter(Suspension == 1, !is.na(MONTHS_OF_TREAT)) %>%
  summarise(
    any_post = any(t_abs > MONTHS_OF_TREAT & in_window == TRUE),
    .groups = "drop"
  )

print(check_post_dropout)

df_long_imp1 <- complete(imp_long, 1)
## -------------------------
## 3) Quick distribution comparison (one variable)
## -------------------------
compare_var <- function(var) {
  rbind(
    original = summary(df_long_work[[var]]),   # CHANGED (was df_long)
    locf     = summary(df_locf[[var]]),
    nocb     = summary(df_nocb[[var]]),
    interp   = summary(df_interp[[var]]),
    mice     = summary(df_long_imp1[[var]])
  )
}

compare_var("MMDs")





## -------------------------
## 4) Masking-based RMSE benchmark (properly aligned)
## -------------------------
calc_rmse <- function(true, imputed) {
  sqrt(mean((true - imputed)^2, na.rm = TRUE))
}

# Build truth set with a stable row_id (from dropout-safe data)
truth_set <- df_long_work %>%                    # CHANGED
  filter(!is.na(MMDs)) %>%
  mutate(row_id = row_number())

set.seed(42)
n_mask   <- floor(0.2 * nrow(truth_set))
mask_ids <- sample(truth_set$row_id, n_mask)

test_data <- truth_set %>%
  mutate(MMDs_masked = ifelse(row_id %in% mask_ids, NA, MMDs))

true_vals <- truth_set %>%
  filter(row_id %in% mask_ids) %>%
  pull(MMDs)

## ---- LOCF (within SUBJECT_ID, CYCLE)
locf_filled <- test_data %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID, CYCLE) %>%
  mutate(MMDs_f = zoo::na.locf(MMDs_masked, na.rm = FALSE)) %>%
  ungroup() %>%
  filter(row_id %in% mask_ids) %>%
  pull(MMDs_f)

## ---- NOCB (within SUBJECT_ID, CYCLE)
nocb_filled <- test_data %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID, CYCLE) %>%
  mutate(MMDs_f = zoo::na.locf(MMDs_masked, fromLast = TRUE, na.rm = FALSE)) %>%
  ungroup() %>%
  filter(row_id %in% mask_ids) %>%
  pull(MMDs_f)

## ---- Interpolation (within SUBJECT_ID, CYCLE)
interp_filled <- test_data %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID, CYCLE) %>%
  mutate(MMDs_f = zoo::na.approx(MMDs_masked, x = MONTH, na.rm = FALSE)) %>%
  ungroup() %>%
  filter(row_id %in% mask_ids) %>%
  pull(MMDs_f)


m_imp <- 20

mice_test_df <- test_data %>%
  transmute(
    row_id,
    SUBJECT_ID,
    CYCLE,
    MONTH,
    MMDs = MMDs_masked
  )

ini2  <- mice(mice_test_df, maxit = 0, printFlag = FALSE)
meth2 <- ini2$method
pred2 <- ini2$predictorMatrix

meth2["MMDs"] <- "2l.pmm"
meth2[c("row_id","SUBJECT_ID","CYCLE","MONTH")] <- ""

pred2[,] <- 0
pred2["MMDs", "SUBJECT_ID"] <- -2
pred2["MMDs", c("CYCLE","MONTH")] <- 1
pred2["MMDs", "row_id"] <- 0

set.seed(123)
imp_test <- mice(
  mice_test_df,
  method = meth2,
  predictorMatrix = pred2,
  m = m_imp,
  maxit = 10,
  printFlag = FALSE
)

calc_rmse <- function(true, imputed) sqrt(mean((true - imputed)^2, na.rm = TRUE))

rmse_mice_each <- sapply(seq_len(m_imp), function(j) {
  comp_j <- complete(imp_test, j)
  
  filled_j <- comp_j %>%
    filter(row_id %in% mask_ids) %>%
    pull(MMDs)
  
  calc_rmse(true_vals, filled_j)
})

rmse_mice_mean <- mean(rmse_mice_each)
rmse_mice_sd   <- sd(rmse_mice_each)

filled_matrix <- sapply(seq_len(m_imp), function(j) {
  complete(imp_test, j) %>%
    filter(row_id %in% mask_ids) %>%
    pull(MMDs)
})

mice_pooled_mean <- rowMeans(filled_matrix, na.rm = TRUE)
rmse_mice_pooled <- calc_rmse(true_vals, mice_pooled_mean)

cat("MICE RMSE across imputations:\n")
cat("  mean =", rmse_mice_mean, "\n")
cat("  sd   =", rmse_mice_sd, "\n\n")
cat("MICE pooled-mean RMSE =", rmse_mice_pooled, "\n")



df_long_work %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID, CYCLE) %>%
  summarise(
    n = n(),
    na_total = sum(is.na(MMDs)),
    na_first = is.na(MMDs[which.min(MONTH)]),
    na_last  = is.na(MMDs[which.max(MONTH)]),
    .groups = "drop"
  ) %>%
  summarise(
    cycles_with_any_na = sum(na_total > 0),
    cycles_na_at_first = sum(na_first, na.rm=TRUE),
    cycles_na_at_last  = sum(na_last,  na.rm=TRUE)
  )





library(dplyr)

first_minus_last_sched <- function(df) {
  df %>%
    group_by(SUBJECT_ID, CYCLE) %>%
    summarise(
      mmd_m1  = MMDs[MONTH == 1][1],
      mmd_m12 = MMDs[MONTH == 12][1],
      diff_mmd = mmd_m1 - mmd_m12,
      .groups = "drop"
    ) %>%
    filter(!is.na(mmd_m1) & !is.na(mmd_m12))
}

  
fit_overall <- with(imp_long, {
  dfk <- data.frame(SUBJECT_ID, CYCLE, MONTH, MMDs)
  d   <- first_minus_last_sched(dfk)
  lm(diff_mmd ~ 1, data = d)
})

summary(pool(fit_overall), conf.int = TRUE)

  
fit_cycle <- with(imp_long, {
  dfk <- data.frame(SUBJECT_ID, CYCLE, MONTH, MMDs)
  d   <- first_minus_last_sched(dfk)
  lm(diff_mmd ~ factor(CYCLE) - 1, data = d)
})

summary(pool(fit_cycle), conf.int = TRUE)

library(lme4)

fit_baseline_trend <- with(imp_long, {
  dfk <- data.frame(SUBJECT_ID, CYCLE, MONTH, MMDs)
  dfb <- dfk %>% filter(MONTH == 1)
  lmer(MMDs ~ as.numeric(CYCLE) + (1 | SUBJECT_ID), data = dfb)
})

summary(pool(fit_baseline_trend), conf.int = TRUE)

fit_cycle_test <- with(imp_long, {
  dfk <- data.frame(SUBJECT_ID, CYCLE, MONTH, MMDs)
  d   <- first_minus_last_sched(dfk)
  lm(diff_mmd ~ factor(CYCLE), data = d)
})

summary(pool(fit_cycle_test), conf.int = TRUE)
  
  
  
  
  
  
  
  
  
  
  
  check_range <- function(x, lo = -Inf, hi = Inf) {
    !is.na(x) & (x < lo | x > hi)
  }
  
  df_long_imp1 <- complete(imp_long, 1)
  
  bad_hit6_imp  <- check_range(df_long_imp1$HIT6,  36, 78)
  bad_hadsa_imp <- check_range(df_long_imp1$HADSA,  0, 21)
  bad_hadsd_imp <- check_range(df_long_imp1$HADSD,  0, 21)
  
  sum(bad_hit6_imp)
  sum(bad_hadsa_imp)
  sum(bad_hadsd_imp)
  
  df_long_imp1[bad_hadsa_imp, c("SUBJECT_ID","CYCLE","MONTH","HADSA")]
  







count_missing_not_lfu <- function(df) {
  df %>%
    filter(in_window) %>%                 # exclude post-LFU rows
    summarise(
      total_missing = sum(is.na(.)),
      .groups = "drop"
    )
}

count_missing_not_lfu(df_long)

missing_by_variable <- df_long %>%
  filter(in_window) %>%                   # only valid follow-up
  summarise(across(
    everything(),
    ~ sum(is.na(.)),
    .names = "n_missing_{.col}"
  )) %>%
  pivot_longer(
    everything(),
    names_to = "variable",
    values_to = "n_missing"
  ) %>%
  arrange(desc(n_missing))

missing_by_variable
















#-------------------------------------------

# -------- 1) Choose output folder
folder <- "C:/Users/aless/Desktop/medical applications/data/imputed/dropout_not_imputed_longitudinal"
dir.create(folder, showWarnings = FALSE, recursive = TRUE)

# -------- 2) Define what the "pure longitudinal" columns are
# Use the columns present in a completed imputed dataset as the template.
long_cols <- names(complete(imp_long, 1))

# -------- 3) Extract post-LFU rows (NOT imputed), keeping ONLY longitudinal columns
df_post_lfu <- df_long %>%
  filter(!in_window) %>%              # rows excluded from imputation
  select(all_of(long_cols))           # keep ONLY longitudinal columns

# -------- 4) Build and save FULL datasets for each imputation
for (j in 1:imp_long$m) {
  df_in_window_imp <- complete(imp_long, j) %>%
    select(all_of(long_cols))
  
  df_full_j <- bind_rows(df_in_window_imp, df_post_lfu) %>%
    arrange(SUBJECT_ID, CYCLE, MONTH)
  
  out_name <- sprintf("imp_dropout_long_%02d.csv", j)
  write.csv(df_full_j, file = file.path(folder, out_name), row.names = FALSE)
}

# -------- 5) Also save one RDS list (optional but convenient)
full_list <- lapply(1:imp_long$m, function(j) {
  bind_rows(
    complete(imp_long, j) %>% select(all_of(long_cols)),
    df_post_lfu
  ) %>% arrange(SUBJECT_ID, CYCLE, MONTH)
})

names(full_list) <- sprintf("imp_%02d", 1:imp_long$m)

saveRDS(
  full_list,
  file = file.path(folder, "full_longitudinal_all_imputations.rds")
)

cat("Saved full longitudinal datasets to:\n", folder, "\n")




























#----------------------------------------------------
# Check baseline sizes
n_full      <- nrow(df_long)
n_in_window <- nrow(df_long_work)

# Check one completed in-window imputation
n_imp_in_window <- nrow(complete(imp_long, 1))

# Check one full longitudinal CSV
test_full <- read.csv(file.path(folder, "imp_dropout_long_01.csv"))
n_full_csv <- nrow(test_full)

c(
  full_longitudinal_augmented = n_full,
  in_window_only              = n_in_window,
  imputed_in_window           = n_imp_in_window,
  full_longitudinal_saved     = n_full_csv
)








