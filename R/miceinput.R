
library(moments)
library(mice)

df_baseline_clean <- read.csv("C:/Users/aless/Desktop/medical applications/data/cleaned/df_baseline_clean.csv", header=FALSE)
View(df_baseline_clean)

results <- data.frame(
  variable = character(),
  n_obs = integer(),
  shapiro_p = numeric(),
  skewness = numeric(),
  kurtosis = numeric(),
  stringsAsFactors = FALSE
)

for (col in names(df_baseline_clean)) {
  
  # Check if column has missing values
  if (any(is.na(df_baseline_clean[[col]]))) {
    
    x <- df_baseline_clean[[col]]
    
    # Keep only numeric columns
    if (is.numeric(x)) {
      
      x_obs <- x[!is.na(x)]
      
      # Shapiro-Wilk only valid for n <= 5000
      shapiro_p <- if (length(x_obs) >= 3 && length(x_obs) <= 5000) {
        shapiro.test(x_obs)$p.value
      } else {
        NA
      }
      
      results <- rbind(results, data.frame(
        variable = col,
        n_obs = length(x_obs),
        shapiro_p = shapiro_p,
        skewness = skewness(x_obs),
        kurtosis = kurtosis(x_obs)
      ))
    }
  }
}

print(results)



bin_vars <- c(
  "SEX", "Suspension", "FAMILIARITY", "PULSATING", "PAIN_MOVMENT", "Aura",
  "Bbloc", "Caant", "Tricyclic", "Antiepil", "SSRISNRI", "Antiipnt",
  "Botulin", "DETOXPRE", "Psycopathological", "Hypertension",
  "Sleep_Disorders", "PREV_T0"
)


multi_vars <- c("DIAGNOSIS", "SIDE", "T0_SYMPT_TREATMENT", "ANTIBODY", "TREATMENT_DISC")


df_baseline_clean[bin_vars]   <- lapply(df_baseline_clean[bin_vars], factor)
df_baseline_clean[multi_vars] <- lapply(df_baseline_clean[multi_vars], factor)


ini  <- mice(df_baseline_clean, maxit = 0)
meth <- ini$method
pred <- ini$predictorMatrix




if ("SUBJECT_ID" %in% names(meth)) meth["SUBJECT_ID"] <- ""


cont_vars <- c(
  "AGE", "MONTHS_OF_TREAT", "AGE_OF_ONSET", "AGE_W_CHRONICMIGRAINE",
  "WEIGTH", "HEIGTH",
  "GGCEF_T0",   # monthly headache days at baseline (count-like, PMM robust)
  "INT_T0"      # severity 0-10 (bounded, PMM robust)
)
cont_vars <- intersect(cont_vars, names(meth))
meth[cont_vars] <- "pmm"

# (C) Binary -> logistic regression
meth[bin_vars] <- "logreg"

# (D) Nominal multicategory -> multinomial regression
meth[multi_vars] <- "polyreg"



visit_order <- c(
  # Layer 1: anchors
  "AGE", "SEX", "DIAGNOSIS", "AGE_OF_ONSET", "AGE_W_CHRONICMIGRAINE", "FAMILIARITY",
  
  # Layer 2: anthropometrics
  "WEIGTH", "HEIGTH",
  
  # Layer 3: special clinical history
  "Sleep_Disorders", "PREV_T0", "INT_T0",
  
  # Layer 4: symptoms
  "SIDE", "PULSATING", "PAIN_MOVMENT", "Aura", "T0_SYMPT_TREATMENT",
  
  # Layer 5: treatments
  "Bbloc", "Caant", "Tricyclic", "Antiepil", "SSRISNRI", "Antiipnt",
  "Pizotifene", "Botulin", "DETOXPRE"
)


# Keep only columns that actually have missing values
if (length(visit_order) > 0) {
  visit_order <- visit_order[colMeans(is.na(df_baseline_clean[visit_order])) > 0]
}



#PRED MATRIX 


# Never use SUBJECT_ID as predictor; never impute it
if ("SUBJECT_ID" %in% colnames(pred)) pred[, "SUBJECT_ID"] <- 0
if ("SUBJECT_ID" %in% rownames(pred)) pred["SUBJECT_ID", ] <- 0

anchor_vars <- c("AGE","SEX","DIAGNOSIS","FAMILIARITY",
                 "AGE_OF_ONSET","AGE_W_CHRONICMIGRAINE")

clinical_vars <- c(
  "WEIGTH","HEIGTH","BMI",
  "SIDE","PULSATING","PAIN_MOVMENT","Aura",
  "T0_SYMPT_TREATMENT",
  "Psycopathological","Hypertension"
)

treat_vars <- c(
  "Bbloc","Caant","Tricyclic","Antiepil","SSRISNRI","Antiipnt",
  "Pizotifene","Botulin","DETOXPRE"
)

other_vars <- c("Sleep_Disorders","PREV_T0","INT_T0")

anchor_vars   <- anchor_vars[anchor_vars %in% colnames(pred)]
clinical_vars <- clinical_vars[clinical_vars %in% colnames(pred)]
treat_vars    <- treat_vars[treat_vars %in% colnames(pred)]
other_vars    <- other_vars[other_vars %in% colnames(pred)]

pred[,] <- 0
pred[clinical_vars, anchor_vars] <- 1
pred[treat_vars, anchor_vars] <- 1
pred[other_vars, anchor_vars] <- 1
pred[anchor_vars, ] <- 0

diag(pred) <- 0


pred[clinical_vars, clinical_vars] <- 1
diag(pred) <- 0

pred[treat_vars, clinical_vars] <- 1

pred[clinical_vars, other_vars] <- 1
pred[treat_vars, other_vars] <- 1
pred[treat_vars, treat_vars] <- 0

pred[anchor_vars, treat_vars] <- 0
pred[clinical_vars, treat_vars] <- 0




# Never impute BMI and never use it as predictor (deterministic)
if ("BMI" %in% names(meth)) meth["BMI"] <- ""

if ("BMI" %in% colnames(pred)) pred[, "BMI"] <- 0
if ("BMI" %in% rownames(pred)) pred["BMI", ] <- 0



# ---------------------------
# 6) Run MICE
# ---------------------------
imp <- mice(
  df_baseline_clean,
  method = meth,
  predictorMatrix = pred,
  visitSequence = visit_order,   # your custom order
  m = 10,
  maxit = 30,
  seed = 123
)

# ---------------------------
# 7) Diagnostics
# ---------------------------
plot(imp)                         # convergence (all variables)
# densityplot(imp, ~ INT_T0)      # example density check
# densityplot(imp, ~ BMI)

# ---------------------------
# 8) Extract a completed dataset (one of the m imputations)
# ---------------------------
df_baseline_imp1 <- complete(imp, 1)

# ---------------------------
# 9) (Optional) enforce measurement constraints after imputation
# ---------------------------

# INT_T0: clamp to [0, 10] and make integer
if ("INT_T0" %in% names(df_baseline_imp1)) {
  df_baseline_imp1$INT_T0 <- pmin(pmax(round(df_baseline_imp1$INT_T0), 0), 10)
}

# GGCEF_T0: non-negative integer count
if ("GGCEF_T0" %in% names(df_baseline_imp1)) {
  df_baseline_imp1$GGCEF_T0 <- pmax(round(df_baseline_imp1$GGCEF_T0), 0)

colSums(is.na(df_baseline_imp1))


vars_imputed <- names(which(colMeans(is.na(df_baseline_clean)) > 0))
vars_imputed

convergence_stats <- data.frame(variable = character(),
                                mean_sd = numeric(),
                                sd_sd = numeric(),
                                stringsAsFactors = FALSE)

for (v in vars_imputed) {
  imp_values <- sapply(1:imp$m, function(i) complete(imp, i)[[v]])
  
  chain_means <- apply(imp_values, 2, mean, na.rm = TRUE)
  chain_sds   <- apply(imp_values, 2, sd, na.rm = TRUE)
  
  convergence_stats <- rbind(convergence_stats,
                             data.frame(variable = v,
                                        mean_sd = sd(chain_means),
                                        sd_sd   = sd(chain_sds)))
}

convergence_stats

colnames(pred)[colSums(pred) > 0]
