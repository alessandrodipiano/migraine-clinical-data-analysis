library(mice)
library(dplyr)
library(zoo)
library(tidyr)
library(miceadds)   
exists("mice.impute.2l.pmm")

df_long <- read.csv("C:/Users/aless/Desktop/medical applications/data/cleaned/df_longitudinal_clean.csv", header=TRUE)
View(df_long)


df_long$SUBJECT_ID <- as.integer(df_long$SUBJECT_ID)
df_long$CYCLE <- as.integer(df_long$CYCLE)
df_long$MONTH <- as.integer(df_long$MONTH)

# Drop derived variable (correct)
if ("RED_MMD_VST01" %in% names(df_long)) df_long <- subset(df_long, select = -RED_MMD_VST01)

cluster_var <- "SUBJECT_ID"
time_vars <- c("CYCLE","MONTH")

cont_vars <- c("MMDs","DOSE","GGFAR","HADSA","HADSD","HIT6","INT","MIDAS")
cont_vars <- intersect(cont_vars, names(df_long))

# Ensure numeric
for (v in cont_vars) df_long[[v]] <- as.numeric(df_long[[v]])

# Init
ini  <- mice(df_long, maxit=0, printFlag=FALSE)
meth <- ini$method
pred <- ini$predictorMatrix

# Methods
meth[cont_vars] <- "2l.pmm"
meth[c(cluster_var, time_vars)] <- ""   # do not impute ID/time

# Predictor matrix
pred[,] <- 0
pred[, cluster_var] <- -2               # random intercept per subject

# Time predicts outcomes
pred[cont_vars, time_vars] <- 1

# Outcomes predict each other (same timepoint)
pred[cont_vars, cont_vars] <- 1
diag(pred) <- 0

# time vars are not targets
pred[time_vars, ] <- 0

# Sanity check
to_impute <- names(meth)[meth != ""]
bad <- to_impute[rowSums(pred[to_impute, , drop=FALSE]) == 0]
print(bad)
stopifnot(length(bad) == 0)

# Run
imp_long <- mice(
  df_long,
  method = meth,
  predictorMatrix = pred,
  m = 20, maxit = 20, seed = 123,
  printFlag = TRUE
)


plot(imp_long)
saveRDS(imp_long, file = "imp_long.rds")



# One completed dataset for quick descriptive comparison
df_long_imp1 <- complete(imp_long, 1)


## -------------------------
df_locf <- df_long %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID, CYCLE) %>%
  fill(all_of(cont_vars), .direction = "down") %>%
  ungroup()

df_nocb <- df_long %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID, CYCLE) %>%
  fill(all_of(cont_vars), .direction = "up") %>%
  ungroup()

# Interpolate only continuous-like variables where it makes sense
interp_vars <- intersect(c("MMDs","HIT6","MIDAS"), names(df_long))

df_interp <- df_long %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID, CYCLE) %>%
  mutate(across(
    all_of(interp_vars),
    ~ zoo::na.approx(.x, x = MONTH, na.rm = FALSE)
  )) %>%
  ungroup()

## -------------------------
## 3) Quick distribution comparison (one variable)
## -------------------------
compare_var <- function(var) {
  rbind(
    original = summary(df_long[[var]]),
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

# Build truth set with a stable row_id
truth_set <- df_long %>%
  filter(!is.na(MMDs)) %>%
  mutate(row_id = row_number())

set.seed(42)
n_mask  <- floor(0.2 * nrow(truth_set))
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

## ---- MICE benchmark (multilevel, consistent with longitudinal structure)
mice_test_df <- test_data %>%
  transmute(SUBJECT_ID, CYCLE, MONTH, MMDs = MMDs_masked)

ini2  <- mice(mice_test_df, maxit = 0, printFlag = FALSE)
meth2 <- ini2$method
pred2 <- ini2$predictorMatrix

meth2["MMDs"] <- "2l.pmm"
meth2[c("SUBJECT_ID","CYCLE","MONTH")] <- ""

pred2[,] <- 0
pred2["MMDs", "SUBJECT_ID"] <- -2
pred2["MMDs", c("CYCLE","MONTH")] <- 1

set.seed(123)
imp_test <- mice(
  mice_test_df,
  method = meth2,
  predictorMatrix = pred2,
  m = 1, maxit = 10,
  printFlag = FALSE
)

mice_completed <- complete(imp_test, 1)

# Reattach row_id in original row order (mice keeps row order)
mice_filled <- mice_completed %>%
  mutate(row_id = test_data$row_id) %>%
  filter(row_id %in% mask_ids) %>%
  pull(MMDs)

## ---- RMSE summary table
results <- data.frame(
  Method = c("MICE_2l_pmm", "LOCF", "NOCB", "Interpolation"),
  RMSE = c(
    calc_rmse(true_vals, mice_filled),
    calc_rmse(true_vals, locf_filled),
    calc_rmse(true_vals, nocb_filled),
    calc_rmse(true_vals, interp_filled)
  )
)

print(results)

library(mice)
library(dplyr)

m_imp <- 20

# Build MICE test data WITH row_id included
mice_test_df <- test_data %>%
  transmute(
    row_id,                  # keep it!
    SUBJECT_ID,
    CYCLE,
    MONTH,
    MMDs = MMDs_masked
  )

ini2  <- mice(mice_test_df, maxit = 0, printFlag = FALSE)
meth2 <- ini2$method
pred2 <- ini2$predictorMatrix

# Impute only MMDs
meth2["MMDs"] <- "2l.pmm"
meth2[c("row_id","SUBJECT_ID","CYCLE","MONTH")] <- ""

# Predictor matrix
pred2[,] <- 0
pred2["MMDs", "SUBJECT_ID"] <- -2
pred2["MMDs", c("CYCLE","MONTH")] <- 1

# Ensure row_id is not used as a predictor (optional but clean)
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

# RMSE per imputation
rmse_mice_each <- sapply(seq_len(m_imp), function(j) {
  comp_j <- complete(imp_test, j)
  
  filled_j <- comp_j %>%
    filter(row_id %in% mask_ids) %>%
    pull(MMDs)
  
  calc_rmse(true_vals, filled_j)
})

rmse_mice_mean <- mean(rmse_mice_each)
rmse_mice_sd   <- sd(rmse_mice_each)

# Pooled (posterior-mean) RMSE
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








first_minus_last <- function(df) {
  df %>%
    group_by(SUBJECT_ID, CYCLE) %>%
    summarise(
      diff_mmd = MMDs[which.min(MONTH)] - MMDs[which.max(MONTH)],
      .groups = "drop"
    )
}

fit <- with(imp_long, {
  dfk <- data.frame(SUBJECT_ID, CYCLE, MONTH, MMDs)  
  d   <- first_minus_last(dfk)
  lm(diff_mmd ~ 1, data = d)                       
})

pooled <- pool(fit)
summary(pooled, conf.int = TRUE)

fit_cycle <- with(imp_long, {
  dfk <- data.frame(SUBJECT_ID, CYCLE, MONTH, MMDs)
  d   <- first_minus_last(dfk)
  lm(diff_mmd ~ factor(CYCLE) - 1, data = d)
})

summary(pool(fit_cycle), conf.int = TRUE)



library(mice)
library(lme4)

fit_baseline_trend <- with(imp_long, {
  dfk <- data.frame(SUBJECT_ID, CYCLE, MONTH, MMDs)
  dfb <- dfk[dfk$MONTH == 1, ]  # baseline within each cycle
  
  # Treat cycle as ordered numeric (tests linear trend C1 -> C2 -> C3)
  lmer(MMDs ~ as.numeric(CYCLE) + (1 | SUBJECT_ID), data = dfb)
})

pooled_baseline_trend <- pool(fit_baseline_trend)
summary(pooled_baseline_trend, conf.int = TRUE)





check_range <- function(x, lo = -Inf, hi = Inf) {
  !is.na(x) & (x < lo | x > hi)
}

bad_hit6 <- check_range(df_long$HIT6, 36, 78)
bad_hadsa <- check_range(df_long$HADSA, 0, 21)
bad_hadsd <- check_range(df_long$HADSD, 0, 21)

sum(bad_hit6)
sum(bad_hadsa)
sum(bad_hadsd)

df_long[bad_hadsa, c("SUBJECT_ID","CYCLE","MONTH","HADSA")]





folder <- "C:/Users/aless/Desktop/medical applications/data/imputed/imputed_longitudinal"

files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)
files <- sort(files)  # ensure imp_01, imp_02, ... order

completed_list <- lapply(files, read.csv, stringsAsFactors = FALSE)

imp_long_reloaded <- miceadds::datlist2mids(completed_list)


fit <- with(imp_long_reloaded, lm(MMDs ~ CYCLE + MONTH))
pool(fit)

