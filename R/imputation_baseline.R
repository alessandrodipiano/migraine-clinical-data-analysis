
library(moments)
library(mice)

df_baseline_clean <- read.csv("C:/Users/aless/Desktop/medical applications/data/cleaned/df_baseline_clean.csv", header=TRUE)
View(df_baseline_clean)
names(df_baseline_clean)[names(df_baseline_clean) == "WEIGTH"] <- "WEIGHT"
names(df_baseline_clean)[names(df_baseline_clean) == "HEIGTH"] <- "HEIGHT"
df_baseline_clean$BMI <- NULL


#GROUPS 
#----------------------------------------------
cont_vars <- c("AGE","MONTHS_OF_TREAT","AGE_OF_ONSET","AGE_W_CHRONICMIGRAINE",
               "WEIGHT","HEIGHT","GGCEF_T0","INT_T0")

bin_vars <- c("SEX","Suspension","FAMILIARITY","PULSATING","PAIN_MOVMENT","Aura",
              "Bbloc","Caant","Tricyclic","Antiepil","SSRISNRI","Antiipnt",
              "Botulin","DETOXPRE","Hypertension",
              "Sleep_Disorders","PREV_T0", "Pizotifene")

multi_vars <- c("DIAGNOSIS","SIDE","T0_SYMPT_TREATMENT","ANTIBODY","TREATMENT_DISC","Psycopathological")


df_baseline_clean[bin_vars]   <- lapply(df_baseline_clean[bin_vars], factor)
df_baseline_clean[multi_vars] <- lapply(df_baseline_clean[multi_vars], factor)

true_anchor <- c("AGE","SEX","DIAGNOSIS")
upstream_vars <- c("AGE_OF_ONSET","AGE_W_CHRONICMIGRAINE","FAMILIARITY")
anthro_vars <- c("WEIGHT","HEIGHT")
baseline_vars <- c("Sleep_Disorders","Hypertension","Psycopathological",
                   "SIDE","PULSATING","PAIN_MOVMENT","Aura",
                   "GGCEF_T0","INT_T0","PREV_T0","T0_SYMPT_TREATMENT")
treat_vars <- c("ANTIBODY","Bbloc","Caant","Tricyclic","Antiepil","SSRISNRI","Antiipnt",
                "Pizotifene","Botulin","DETOXPRE","MONTHS_OF_TREAT","Suspension","TREATMENT_DISC")


all_cols <- colnames(df_baseline_clean)

# Intersect groups with existing columns
true_anchor   <- intersect(true_anchor,   all_cols)
upstream_vars <- intersect(upstream_vars, all_cols)
anthro_vars   <- intersect(anthro_vars,   all_cols)
baseline_vars <- intersect(baseline_vars, all_cols)
treat_vars    <- intersect(treat_vars,    all_cols)

cont_vars  <- intersect(cont_vars,  all_cols)
bin_vars   <- intersect(bin_vars,   all_cols)
multi_vars <- intersect(multi_vars, all_cols)


#-----------------------------------------


#initialization
#-----------------------
ini  <- mice(df_baseline_clean, maxit = 0, printFlag = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix


meth[cont_vars]  <- "pmm"
meth[bin_vars]   <- "logreg"
meth[multi_vars] <- "polyreg"


if ("SUBJECT_ID" %in% names(meth)) meth["SUBJECT_ID"] <- ""
if ("SUBJECT_ID" %in% colnames(pred)) pred[, "SUBJECT_ID"] <- 0
if ("SUBJECT_ID" %in% rownames(pred)) pred["SUBJECT_ID", ] <- 0



if ("AGE_OF_ONSET" %in% upstream_vars)          meth["AGE_OF_ONSET"] <- "pmm"
if ("AGE_W_CHRONICMIGRAINE" %in% upstream_vars) meth["AGE_W_CHRONICMIGRAINE"] <- "pmm"
if ("FAMILIARITY" %in% upstream_vars)           meth["FAMILIARITY"] <- "logreg" 


meth[true_anchor] <- ""

#PRED
#-------------------------------------------

pred[,] <- 0

# Upstream predicted by true anchors
pred[upstream_vars, true_anchor] <- 1

# Anthropometrics predicted by true anchors + upstream
pred[anthro_vars, c(true_anchor, upstream_vars)] <- 1

# Baseline predicted by true anchors + upstream + anthropometrics
pred[baseline_vars, c(true_anchor, upstream_vars, anthro_vars)] <- 1

# Treatments predicted by everything before them
pred[treat_vars, c(true_anchor, upstream_vars, anthro_vars, baseline_vars)] <- 1

# Block downstream -> upstream
pred[upstream_vars, c(anthro_vars, baseline_vars, treat_vars)] <- 0
pred[anthro_vars,   treat_vars] <- 0
pred[baseline_vars, treat_vars] <- 0

# prevent baseline↔baseline circularity 
pred[baseline_vars, baseline_vars] <- 0
diag(pred) <- 0

# Safety: every imputed var must have >=1 predictor
to_impute <- names(meth)[meth != ""]
bad <- to_impute[rowSums(pred[to_impute, , drop = FALSE]) == 0]
print(bad)
stopifnot(length(bad) == 0)

# -----------------------------
# 5) TWO-STAGE IMPUTATION
# -----------------------------
# Stage 1: initialize predictors so categorical models don't see NA covariates
stage1_vars <- unique(intersect(
  c(upstream_vars, anthro_vars, "GGCEF_T0", "INT_T0", "MONTHS_OF_TREAT"),
  names(meth)
))

meth1 <- meth
meth1[setdiff(names(meth1), stage1_vars)] <- ""  # only impute stage1 vars
meth1[true_anchor] <- ""                         # still never impute anchors

pred1 <- pred
pred1[,] <- 0
pred1[stage1_vars, true_anchor] <- 1

# allow stage1 vars to use each other ONLY if you want (usually not needed)
# pred1[stage1_vars, stage1_vars] <- 0
diag(pred1) <- 0

where1 <- is.na(df_baseline_clean)

imp_stage1 <- mice(
  df_baseline_clean,
  method          = meth1,
  predictorMatrix = pred1,
  where           = where1,
  m               = 1,
  maxit           = 10,
  seed            = 123,
  printFlag       = FALSE
)

df_stage1 <- complete(imp_stage1, 1)

# Stage 2: full imputation on top of initialized predictors
where2 <- is.na(df_stage1)

visit2 <- unique(c(upstream_vars, anthro_vars, baseline_vars, treat_vars, to_impute))
visit2 <- visit2[visit2 %in% to_impute]

imp <- mice(
  df_stage1,
  method          = meth,
  predictorMatrix = pred,
  where           = where2,
  visitSequence   = visit2,
  m               = 10,
  maxit           = 20,
  seed            = 123,
  printFlag       = TRUE
)

# -----------------------------
# 6) Validate: imputed vars should have no NA
# -----------------------------
df1 <- complete(imp, 1)
should_be_complete <- names(imp$method)[imp$method != ""]
na_imputed <- colMeans(is.na(df1[, should_be_complete, drop = FALSE]))
print(na_imputed[na_imputed > 0])   # should be numeric(0)

# Optional diagnostics
plot(imp)





#TEST numerical

#--------------------------------------------------------
to_impute <- names(meth)[meth != ""]

visit_all <- to_impute



imp_default_order <- mice(
  df_baseline_clean,
  method = meth,
  predictorMatrix = pred,
  where = is.na(df_baseline_clean),
  m = 10, maxit = 20, seed = 123,
  printFlag = FALSE
)

imp_reverse_order <- mice(
  df_baseline_clean,
  method = meth,
  predictorMatrix = pred,
  where = is.na(df_baseline_clean),
  visitSequence = rev(visit_all),
  m = 10, maxit = 20, seed = 123,
  printFlag = FALSE
)

vars_num <- names(meth)[meth == "pmm"]

mean_by_imp <- function(imp_obj, v) {
  sapply(1:imp_obj$m, function(k) mean(complete(imp_obj, k)[[v]], na.rm = TRUE))
}

sd_by_imp <- function(imp_obj, v) {
  sapply(1:imp_obj$m, function(k) sd(complete(imp_obj, k)[[v]], na.rm = TRUE))
}

cmp_num <- do.call(rbind, lapply(vars_num, function(v) {
  mc <- mean_by_imp(imp, v)
  md <- mean_by_imp(imp_default_order, v)
  mr <- mean_by_imp(imp_reverse_order, v)
  
  sc <- sd_by_imp(imp, v)
  sd_ <- sd_by_imp(imp_default_order, v)
  sr <- sd_by_imp(imp_reverse_order, v)
  
  data.frame(
    variable = v,
    
    mean_current = mean(mc), mean_default = mean(md), mean_reverse = mean(mr),
    absdiff_cur_def = abs(mean(mc) - mean(md)),
    absdiff_cur_rev = abs(mean(mc) - mean(mr)),
    
    sd_current = mean(sc), sd_default = mean(sd_), sd_reverse = mean(sr),
    absdiff_sd_cur_def = abs(mean(sc) - mean(sd_)),
    absdiff_sd_cur_rev = abs(mean(sc) - mean(sr)),
    
    stringsAsFactors = FALSE
  )
}))

cmp_num




















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

fit1 <- with(imp, glm(Suspension ~ AGE + SEX + ANTIBODY + ..., family=binomial))
p1 <- pool(fit1)

fit2 <- with(imp_reverse_order, glm(Suspension ~ AGE + SEX + ANTIBODY + ..., family=binomial))
p2 <- pool(fit2)

summary(p1)
summary(p2)

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


should_be_complete <- names(imp$method)[imp$method != ""]
na_imputed <- colMeans(is.na(df1[, should_be_complete, drop = FALSE]))
na_imputed[na_imputed > 0]










