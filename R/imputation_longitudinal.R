library(mice)

library(miceadds)   
exists("mice.impute.2l.pmm")




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






library(dplyr)
library(tidyr)

df_locf <- df_long %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID) %>%
  fill(MMDs, DOSE, GGFAR, HADSA, HADSD, HIT6, INT, MIDAS, .direction = "down") %>%
  ungroup()



df_nocb <- df_long %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID) %>%
  fill(MMDs, DOSE, GGFAR, HADSA, HADSD, HIT6, INT, MIDAS, .direction = "up") %>%
  ungroup()

install.packages('zoo')

library(zoo)

df_interp <- df_long %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID) %>%
  mutate(
    MMDs = na.approx(MMDs, x = MONTH, na.rm = FALSE),
    HIT6 = na.approx(HIT6, x = MONTH, na.rm = FALSE),
    MIDAS = na.approx(MIDAS, x = MONTH, na.rm = FALSE)
  ) %>%
  ungroup()

df_long_imp1 <- complete(imp_long, 1)
compare <- function(var) {
  rbind(
    original = summary(df_long[[var]]),
    locf     = summary(df_locf[[var]]),
    nocb     = summary(df_nocb[[var]]),
    interp   = summary(df_interp[[var]]),
    mice     = summary(df_long_imp1[[var]])
  )
}
compare("MMDs")


# 1. Select only rows where MMDs is OBSERVED (the "Ground Truth" set)
truth_set <- df_long[!is.na(df_long$MMDs), ]

# 2. Artificially mask (delete) 20% of these known values
set.seed(42)
n_mask <- floor(0.2 * nrow(truth_set))
mask_indices <- sample(1:nrow(truth_set), n_mask)

test_data <- truth_set
test_data$MMDs[mask_indices] <- NA  # Create artificial NAs

# --- METHOD A: MICE (PMM) ---
# We use a simplified run for speed in this test
imp_test <- mice(test_data, method="pmm", m=1, maxit=5, printFlag=FALSE)
mice_filled <- complete(imp_test)$MMDs[mask_indices]

# --- METHOD B: LOCF ---
# Note: Requires sorting by Subject/Time
locf_filled <- test_data %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID) %>%
  fill(MMDs, .direction = "down") %>%
  ungroup() %>%
  slice(mask_indices) %>%
  pull(MMDs)

# --- METHOD C: Linear Interpolation ---
interp_filled <- test_data %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID) %>%
  mutate(MMDs = na.approx(MMDs, x=MONTH, na.rm=FALSE)) %>%
  ungroup() %>%
  slice(mask_indices) %>%
  pull(MMDs)

# --- EVALUATION (RMSE) ---
calc_rmse <- function(true, imputed) {
  sqrt(mean((true - imputed)^2, na.rm = TRUE))
}


# --- METHOD D: NOCB (Next Observation Carried Backward) ---
# We use .direction = "up" to pull future values backward
nocb_filled <- test_data %>%
  arrange(SUBJECT_ID, CYCLE, MONTH) %>%
  group_by(SUBJECT_ID) %>%
  fill(MMDs, .direction = "up") %>% # The crucial difference
  ungroup() %>%
  slice(mask_indices) %>%
  pull(MMDs)

# Calculate RMSE for NOCB
rmse_nocb <- calc_rmse(truth_set$MMDs[mask_indices], nocb_filled)



results <- data.frame(
  Method = c("MICE", "LOCF", "Interpolation"),
  RMSE = c(
    calc_rmse(truth_set$MMDs[mask_indices], mice_filled),
    calc_rmse(truth_set$MMDs[mask_indices], locf_filled),
    calc_rmse(truth_set$MMDs[mask_indices], interp_filled)
  )
)
# Update the results table
results <- rbind(results, data.frame(Method = "NOCB", RMSE = rmse_nocb))


print(results)





