
library(moments)
library(mice)

df_baseline_clean <- read.csv("C:/Users/aless/Desktop/medical applications/data/cleaned/df_baseline_clean.csv", header=TRUE)
View(df_baseline_clean)


df_baseline_clean[] <- lapply(df_baseline_clean, function(x) {
  # trim whitespace
  x <- trimws(as.character(x))
  # convert "" to NA
  x[x == ""] <- NA
  
  # if it looks numeric, convert
  suppressWarnings({
    xn <- as.numeric(x)
  })
  if (all(is.na(xn) == is.na(x))) xn else x
})


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


#----------------------------------------------
cont_vars <- c("AGE","MONTHS_OF_TREAT","AGE_OF_ONSET","AGE_W_CHRONICMIGRAINE",
               "WEIGTH","HEIGTH","GGCEF_T0","INT_T0")

bin_vars <- c("SEX","Suspension","FAMILIARITY","PULSATING","PAIN_MOVMENT","Aura",
              "Bbloc","Caant","Tricyclic","Antiepil","SSRISNRI","Antiipnt",
              "Botulin","DETOXPRE","Psycopathological","Hypertension",
              "Sleep_Disorders","PREV_T0")

multi_vars <- c("DIAGNOSIS","SIDE","T0_SYMPT_TREATMENT","ANTIBODY","TREATMENT_DISC")


df_baseline_clean[bin_vars]   <- lapply(df_baseline_clean[bin_vars], factor)
df_baseline_clean[multi_vars] <- lapply(df_baseline_clean[multi_vars], factor)
#--------------------------------------------------------------------------------

anchor_vars <- c("AGE","SEX","DIAGNOSIS","FAMILIARITY","AGE_OF_ONSET","AGE_W_CHRONICMIGRAINE")
anthro_vars <- c("WEIGTH","HEIGTH","BMI")
baseline_vars <- c("Sleep_Disorders","Hypertension","Psycopathological",
                   "SIDE","PULSATING","PAIN_MOVMENT","Aura",
                   "GGCEF_T0","INT_T0","PREV_T0","T0_SYMPT_TREATMENT")
treat_vars <- c("ANTIBODY","Bbloc","Caant","Tricyclic","Antiepil","SSRISNRI","Antiipnt",
                "Pizotifene","Botulin","DETOXPRE","MONTHS_OF_TREAT","Suspension","TREATMENT_DISC")


all_cols <- colnames(df_baseline_clean)
anchor_vars   <- intersect(anchor_vars,   all_cols)
anthro_vars   <- intersect(anthro_vars,   all_cols)
baseline_vars <- intersect(baseline_vars, all_cols)
treat_vars    <- intersect(treat_vars,    all_cols)

cont_vars  <- intersect(cont_vars,  all_cols)
bin_vars   <- intersect(bin_vars,   all_cols)
multi_vars <- intersect(multi_vars, all_cols)

#-------------------------------------------------------------
ini  <- mice(df_baseline_clean, maxit = 0, printFlag = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix


meth[cont_vars]  <- "pmm"
meth[bin_vars]   <- "logreg"
meth[multi_vars] <- "polyreg"

if ("SUBJECT_ID" %in% names(meth)) meth["SUBJECT_ID"] <- ""
if ("SUBJECT_ID" %in% colnames(pred)) pred[, "SUBJECT_ID"] <- 0
if ("SUBJECT_ID" %in% rownames(pred)) pred["SUBJECT_ID", ] <- 0

if ("BMI" %in% names(meth)) meth["BMI"] <- ""
if ("BMI" %in% colnames(pred)) pred[, "BMI"] <- 0
if ("BMI" %in% rownames(pred)) pred["BMI", ] <- 0

# Anchors not imputed
meth[anchor_vars] <- ""

# --- 6) Predictor matrix causal direction
pred[,] <- 0
pred[anthro_vars,   anchor_vars] <- 1
pred[baseline_vars, c(anchor_vars, anthro_vars)] <- 1
pred[treat_vars,    c(anchor_vars, anthro_vars, baseline_vars)] <- 1

# No downstream -> upstream
pred[baseline_vars, treat_vars] <- 0
pred[anthro_vars,   treat_vars] <- 0
pred[anchor_vars,   treat_vars] <- 0

# baseline helps baseline
pred[baseline_vars, baseline_vars] <- 1
diag(pred) <- 0

# treatments don't predict treatments
pred[treat_vars, treat_vars] <- 0

# --- 7) Sanity check
to_impute <- names(meth)[meth != ""]
bad <- to_impute[rowSums(pred[to_impute, , drop = FALSE]) == 0]
print(bad)
stopifnot(length(bad) == 0)

# variables that will actually be imputed
to_impute <- names(meth)[meth != ""]

# define visit order = causal order, restricted to imputed vars with missingness
visit_order <- c(
  # Anchors are excluded automatically (meth == "")
  
  # Anthropometrics
  anthro_vars,
  
  # Baseline clinical state
  baseline_vars,
  
  # Treatments (downstream)
  treat_vars
)

# keep only variables that are imputed and have missing values
visit_order <- visit_order[visit_order %in% to_impute]
visit_order <- visit_order[colMeans(is.na(df_baseline_clean[visit_order, drop = FALSE])) > 0]

# run mice
imp <- mice(
  df_baseline_clean,
  method           = meth,
  predictorMatrix  = pred,
  visitSequence    = visit_order,
  m                = 20,      # recommended for your block-missingness
  maxit            = 10,
  seed             = 123,
  printFlag        = TRUE
)

plot(imp)







imp_default_order <- mice(
  df_baseline_clean,
  method = meth,
  predictorMatrix = pred,
  m = 20, maxit = 10, seed = 123
)

imp_reverse_order <- mice(
  df_baseline_clean,
  method = meth,
  predictorMatrix = pred,
  visitSequence = rev(visit_order),
  m = 20, maxit = 10, seed = 123
)


vars_num <- names(meth)[meth == "pmm"]

mean_by_imp <- function(imp_obj, v) {
  sapply(1:imp_obj$m, function(k) mean(complete(imp_obj, k)[[v]], na.rm = TRUE))
}

cmp <- lapply(vars_num, function(v) {
  c(
    v = v,
    mean_current = mean(mean_by_imp(imp, v)),
    mean_default = mean(mean_by_imp(imp_default_order, v)),
    mean_reverse = mean(mean_by_imp(imp_reverse_order, v))
  )
})

do.call(rbind, cmp)



# categorical variables that were imputed
cat_vars <- names(meth)[meth %in% c("logreg", "polyreg")]

cat_vars <- intersect(cat_vars, names(imp$imp))

prop_stability <- lapply(cat_vars, function(v) {
  props <- sapply(1:imp$m, function(k) {
    tab <- table(complete(imp, k)[[v]])
    tab / sum(tab)
  })
  data.frame(
    variable = v,
    category = rownames(props),
    mean_prop = rowMeans(props),
    sd_prop   = apply(props, 1, sd)
  )
})

prop_stability <- do.call(rbind, prop_stability)
prop_stability


compare_props <- function(imp_obj, var) {
  sapply(1:imp_obj$m, function(k) {
    prop.table(table(complete(imp_obj, k)[[var]]))
  })
}

compare_props(imp, "T0_SYMPT_TREATMENT")
compare_props(imp_reverse_order, "T0_SYMPT_TREATMENT")


# Create a folder
dir.create("imputed_data", showWarnings = FALSE)

# Export each completed dataset
for (k in 1:imp$m) {
  df_k <- complete(imp, k)
  write.csv(
    df_k,
    file = sprintf("imputed_data/imputed_%02d.csv", k),
    row.names = FALSE
  )
}

