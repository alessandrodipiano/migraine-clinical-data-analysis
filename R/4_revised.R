# ==============================================================================
# TASK 4: SURVIVAL ANALYSIS (MULTIPLE IMPUTATION) - SCRIPT R STANDARD
# ==============================================================================

# --- 1. INSTALLAZIONE E CARICAMENTO LIBRERIE ---
required_packages <- c("tidyverse", "survival", "survminer", "mice", "miceadds", "gtsummary", "broom")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(tidyverse)
library(survival)
library(survminer)
library(mice)
library(miceadds)
library(broom)
library(gtsummary)

# Imposta un tema grafico pulito
theme_set(theme_minimal() + theme(text = element_text(size = 12)))

print("--- LIBRERIE CARICATE ---")

# ==============================================================================
# --- 2. CARICAMENTO DATI ---
# ==============================================================================

# !!! ATTENZIONE: MODIFICA QUESTI PERCORSI CON I TUOI !!!
path_longitudinal <- "/Users/alessandro.distefano/Desktop/AI/3rd year/Medical/migraine-clinical-data-analysis/data/imputed/imputed_longitudinal"
path_baseline     <- "/Users/alessandro.distefano/Desktop/AI/3rd year/Medical/migraine-clinical-data-analysis/data/imputed/imputed_data_baseline"

# Lettura file
files_long <- sort(list.files(path_longitudinal, pattern="\\.csv$", full.names=TRUE))
files_base <- sort(list.files(path_baseline, pattern="\\.csv$", full.names=TRUE))

# Controllo sicurezza
if(length(files_long) == 0) stop("Errore: Nessun file trovato. Controlla i percorsi!")
if(length(files_long) != length(files_base)) stop("Errore: Numero diverso di file baseline e longitudinali.")

baseline_list     <- map(files_base, ~ read_csv(.x, show_col_types = FALSE))
longitudinal_list <- map(files_long, ~ read_csv(.x, show_col_types = FALSE))

cat(paste("Caricati correttamente", length(baseline_list), "dataset imputati.\n"))

# ==============================================================================
# --- 3. PREPARAZIONE DATI (TIME-TO-EVENT) ---
# ==============================================================================

make_surv_df <- function(baseline_i, long_i) {
  
  # Join e calcolo dei tempi
  long0 <- long_i %>%
    inner_join(baseline_i %>% select(SUBJECT_ID, AGE, SEX, DIAGNOSIS, ANTIBODY, GGCEF_T0), 
               by = "SUBJECT_ID") %>%
    mutate(
      # Calcolo tempo in mesi (Verifica se la logica (CYCLE-1)*15 è corretta per il tuo studio)
      time_months = (CYCLE - 1) * 15 + MONTH,
      
      # Definizione evento: MMDs <= 50% del baseline
      responder = if_else(!is.na(MMDs) & !is.na(GGCEF_T0) & (MMDs <= 0.5 * GGCEF_T0), 1L, 0L)
    ) %>%
    arrange(SUBJECT_ID, time_months)
  
  # Collasso a livello di soggetto (Survival format)
  surv_df <- long0 %>%
    group_by(SUBJECT_ID) %>%
    summarise(
      time = {
        t_event <- time_months[responder == 1L]
        if (length(t_event) > 0) min(t_event) # Tempo al primo evento
        else max(time_months, na.rm = TRUE)   # Censura all'ultimo follow-up
      },
      event = as.integer(any(responder == 1L, na.rm = TRUE)),
      # Covariate
      AGE = first(AGE),
      SEX = first(SEX),
      DIAGNOSIS = first(DIAGNOSIS),
      ANTIBODY = first(ANTIBODY),
      GGCEF_T0 = first(GGCEF_T0),
      .groups = "drop"
    ) %>%
    filter(!is.na(time)) %>%
    mutate(
      SEX = factor(SEX, levels = c(1,2), labels = c("Female","Male")),
      DIAGNOSIS = factor(DIAGNOSIS, levels = c(1,2,3), labels = c("CM","MOH","HFEM")),
      ANTIBODY = factor(ANTIBODY, levels = c(1,2,3), labels = c("Erenumab","Galcanezumab","Fremanezumab"))
    )
  
  return(surv_df)
}

# Creazione lista dataset survival
surv_list <- map2(baseline_list, longitudinal_list, make_surv_df)

print("--- DATI SURVIVAL CREATI ---")
cat("Event rate nel primo dataset:", mean(surv_list[[1]]$event), "\n")

# ==============================================================================
# --- 4. TABLE 1: CARATTERISTICHE BASALI ---
# ==============================================================================
print("--- GENERAZIONE TABLE 1 (Vedi pannello Viewer) ---")

# Usiamo il primo dataset imputato come rappresentativo
t1 <- surv_list[[1]] %>%
  select(AGE, SEX, DIAGNOSIS, ANTIBODY, GGCEF_T0) %>%
  tbl_summary(
    by = ANTIBODY,
    statistic = list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"),
    label = list(AGE ~ "Età", SEX ~ "Sesso", GGCEF_T0 ~ "MMD Basali")
  ) %>%
  add_p() %>%
  add_overall() %>%
  bold_labels()

print(t1) # Questo apparirà nel Viewer o verrà stampato come testo

# ==============================================================================
# --- 5. CURVE DI KAPLAN-MEIER ---
# ==============================================================================
print("--- GENERAZIONE KM PLOT (Vedi pannello Plots) ---")

df1 <- surv_list[[1]] %>%
  dplyr::filter(!is.na(time), !is.na(event), !is.na(ANTIBODY)) %>%
  dplyr::mutate(
    ANTIBODY = droplevels(ANTIBODY)
  )

km_fit <- survival::survfit(survival::Surv(time, event) ~ ANTIBODY, data = df1)

p_km <- survminer::ggsurvplot(
  fit = km_fit,
  data = df1,
  
  # curve
  size = 1.1,
  linetype = "strata",
  conf.int = TRUE,
  conf.int.alpha = 0.18,
  
  # assi e tick
  xlab = "Tempo (mesi)",
  ylab = "Probabilità di NON aver ancora raggiunto risposta ≥50%",
  break.time.by = 6,              # tick ogni 6 mesi (cambia se vuoi)
  xlim = c(0, max(df1$time, na.rm = TRUE)),
  
  # p-value
  pval = TRUE,
  pval.method = TRUE,
  pval.coord = c(0, 0.08),        # sposta il p in basso a sinistra (modifica se serve)
  
  # risk table (pulita)
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.y.text = FALSE,      # evita label lunghissime nella tabella
  risk.table.fontsize = 3.2,
  
  # legenda
  legend.title = "Anticorpo",
  legend.labs = levels(df1$ANTIBODY),
  legend = c(0.8, 0.85),          # posizione dentro il grafico (x,y)
  
  # look & feel
  ggtheme = ggplot2::theme_classic(base_size = 13),
  tables.theme = ggplot2::theme_classic(base_size = 11),
  title = "Time to First 50% Response (Kaplan–Meier)"
)

print(p_km)

# ==============================================================================
# --- 6. MODELLO DI COX (POOLED) ---
# ==============================================================================
print("--- ESECUZIONE POOLED COX MODEL ---")

# 0) (Consigliato) ripulisci eventuali livelli vuoti e assicurati che il riferimento sia quello giusto
surv_list2 <- purrr::map(surv_list, ~ .x %>%
                           dplyr::filter(!is.na(time), !is.na(event), !is.na(ANTIBODY)) %>%
                           dplyr::mutate(
                             ANTIBODY  = droplevels(ANTIBODY),
                             SEX       = droplevels(SEX),
                             DIAGNOSIS = droplevels(DIAGNOSIS),
                             # scegli il riferimento (modifica se vuoi un altro baseline)
                             ANTIBODY  = stats::relevel(ANTIBODY, ref = levels(ANTIBODY)[1])
                           )
)

# 1) Converti lista in oggetto 'mids'
mids_surv <- miceadds::datlist2mids(surv_list2)

# 2) Fitta il modello su ogni dataset
cox_mira <- with(
  mids_surv,
  survival::coxph(survival::Surv(time, event) ~ ANTIBODY + AGE + SEX + DIAGNOSIS + GGCEF_T0)
)

# 3) Pooling dei risultati (Rubin)
cox_pool <- mice::pool(cox_mira)

# 4) Tabella risultati (HR + CI + p-value), ordinata e leggibile
res <- summary(cox_pool, conf.int = TRUE, exponentiate = TRUE)

out <- res %>%
  dplyr::transmute(
    Variabile = term,
    HR = estimate,
    CI95 = sprintf("%.3f – %.3f", conf.low, conf.high),
    p = p.value
  ) %>%
  dplyr::arrange(p)

print(out, row.names = FALSE)

# (Opzionale) salva anche su CSV per report
# readr::write_csv(out, "cox_pooled_results.csv")


# ==============================================================================
# --- 7. FOREST PLOT (MIGLIORATO) ---
# ==============================================================================
print("--- GENERAZIONE FOREST PLOT (MIGLIORATO) ---")

library(ggplot2)
library(dplyr)
library(stringr)

# 1) prendo la tabella pooled (HR già exponentiated)
res_plot <- res_summary %>%
  as_tibble() %>%
  filter(term != "(Intercept)") %>%
  mutate(
    # 2) etichette più umane (adatta se i tuoi nomi sono diversi)
    term_nice = case_when(
      term == "AGE" ~ "Età (per +1 anno)",
      term == "GGCEF_T0" ~ "MMD baseline (per +1)",
      term == "SEXmale" | term == "SEXMale" ~ "Sesso: Maschio vs Femmina",
      term == "DIAGNOSISMOH" ~ "Diagnosi: MOH vs Emicrania cronica",
      term %in% c("DIAGNOSIShigh_freq_ep_migraine", "DIAGNOSISHFEM") ~ "Diagnosi: HFEM vs Emicrania cronica",
      term %in% c("ANTIBODYgalcanezumab", "ANTIBODYGalcanezumab") ~ "Anticorpo: Galcanezumab vs Erenumab",
      term %in% c("ANTIBODYfremanezumab", "ANTIBODYFremanezumab") ~ "Anticorpo: Fremanezumab vs Erenumab",
      TRUE ~ term
    ),
    signif = ifelse(p.value < 0.05, "p < 0.05", "n.s.")
  ) %>%
  # 3) ordine (qui per HR; se vuoi per p: arrange(p.value))
  arrange(estimate) %>%
  mutate(term_nice = factor(term_nice, levels = term_nice))

# 4) forest plot
p_forest <- ggplot(res_plot, aes(x = estimate, y = term_nice)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.18) +
  geom_point(aes(shape = signif), size = 3) +
  scale_x_log10(
    breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4),
    labels = c("0.25","0.5","0.75","1","1.5","2","3","4")
  ) +
  labs(
    title = "Forest plot (Cox pooled su imputazioni multiple)",
    subtitle = "HR > 1 = risposta più rapida; HR < 1 = risposta più lenta",
    x = "Hazard Ratio (scala log)",
    y = NULL,
    shape = ""
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 11)
  )

print(p_forest)

# ==============================================================================
# --- 8. DIAGNOSTICA (SCHOENFELD / PH) ---
# ==============================================================================
print("--- CONTROLLO ASSUNZIONI (Proportional Hazards) ---")

# 1) Fit Cox su tutte le imputazioni (coerente con il tuo modello pooled)
cox_fits_list <- purrr::map(
  surv_list,
  ~ survival::coxph(survival::Surv(time, event) ~ ANTIBODY + AGE + SEX + DIAGNOSIS + GGCEF_T0, data = .x)
)

# 2) Schoenfeld test su tutte le imputazioni
zph_list <- purrr::map(cox_fits_list, survival::cox.zph)

# 3) P-value globale su tutte le imputazioni
global_p <- purrr::map_dbl(zph_list, ~ .x$table["GLOBAL", "p"])
print("P-value GLOBAL (PH) su 20 imputazioni:")
print(summary(global_p))

# 4) Quali covariate violano più spesso? (riassunto p per covariata)
# prendo i p-value della tabella (escludo GLOBAL)
get_pvals <- function(z) {
  tab <- as.data.frame(z$table)
  tab$term <- rownames(tab)
  tab %>%
    dplyr::filter(term != "GLOBAL") %>%
    dplyr::select(term, p = `p`)
}

p_by_imp <- purrr::map_dfr(seq_along(zph_list), function(i) {
  get_pvals(zph_list[[i]]) %>% dplyr::mutate(imp = i)
})

ph_summary <- p_by_imp %>%
  dplyr::group_by(term) %>%
  dplyr::summarise(
    median_p = median(p, na.rm = TRUE),
    min_p = min(p, na.rm = TRUE),
    prop_p_lt_005 = mean(p < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(median_p)

print("Riassunto PH per covariata (su 20 imputazioni):")
print(ph_summary)

# 5) Plot diagnostico SOLO per imputazione 1 (per report/visivo)
print("--- Plot Schoenfeld (imputazione 1) ---")
plot(zph_list[[1]])   # base R: chiaro e leggero

# Se preferisci ggcoxzph (più 'grafico'):
# survminer::ggcoxzph(zph_list[[1]])

print("--- ANALISI COMPLETATA ---")