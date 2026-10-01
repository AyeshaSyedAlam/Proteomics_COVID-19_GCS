################################################################################
## Immune-metabolic plasma proteomic signatures distinguish asymptomatic from
## symptomatic SARS-CoV-2 infection and predict post-COVID sequelae
##
## Analysis code
## Gutenberg COVID-19 Study (GCS) / Gutenberg Post-COVID Study (GPC)
##
## Contents
##   0  Setup: paths, packages, helper functions
##   1  Protein signature: ridge coefficients (Figure 2A)
##   2  Protein score by infection severity (Figure 2, Supplementary Figure 3)
##   3  Pre- and post-infection protein dynamics (Supplementary Figures 4, 5)
##   4  Protein score and acute symptoms (Figure 3A, Figure 4)
##   5  Symptom prevalence threshold: sensitivity analysis (Suppl. Figure 6)
##   6  Protein score and post-COVID sequelae (Figure 3B)
##   7  Protein score and long-term symptoms (Supplementary Figure 9)
##   8  Protein score and cardiovascular risk factors
##   9  Per-protein statistics across severity groups (Suppl. Figure 10)
##  10  Revision analyses
##        10a  Serologically confirmed infections
##        10b  Baseline comparison: with vs without proteomic measurements
##        10c  Export of analysis samples
################################################################################


################################################################################
## 0  SETUP
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(smplot2)
library(ggbeeswarm)
library(forplo)
library(sandwich)
library(openxlsx)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)

## ---- paths -------------------------------------------------------------
files <- "<path to data files>"
out   <- "<path to output folder>"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

f_prot   <- file.path(files, "data_ensemeble_feature_selection_proteins_and_peptides2026-01-22.csv")
f_scores <- file.path(files, "ensemble_feature_selection_scores_2025-08-06.csv")
f_covid  <- file.path(files, "Covid_Variables.csv")
f_sig    <- file.path(files, "relative_protein_signature_new_20260123.csv")
f_seq    <- file.path(files, "6_month_seq.csv")

stopifnot(all(file.exists(c(f_prot, f_scores, f_covid, f_sig, f_seq))))

## Cohort data frames expected in the workspace:
##   d        baseline examination of the GCS
##   d_cati   computer-assisted interviews (severity, symptoms, covariates)

## ---- symptom labels ----------------------------------------------------
symp_vars <- paste0("symp_ynu", 1:61)

symp_names <- c('Fatigue','Malaise','Fever','Chills','Cough','Sore throat',
                'Hoarseness','Runny nose','Headache','Shortness of breath',
                'Pain when breathing','Wheezing','Chest pain','Tachycardia',
                'Olfactory disturbances','Taste disturbances','Visual disturbances',
                'Conjunctivitis','Hearing problems','Tinnitus','Memory problems',
                'Confusion','Concentration problems','Sleep disturbances',
                'Less sleep','More sleep','Dizziness','Fainting spells',
                'Joint pain or swelling','Swollen ankles','Muscle pain',
                'Muscle stiffness','Body aches','Twiching of the limbs',
                'Weakness of the limbs','Tremor','Numbness or tingling',
                'Mood swings','Depression','Loss of interest or pleasure',
                'Anxiety','Abdominal pain','Diarrhea','Constipation',
                'Nausea or vomiting','Loss of appetite','Weight loss',
                'Skin rash','Hair loss','Seizures','Problems urinating',
                'Problems swallowing','Problems with balance',
                'Trouble walking or falling','Slow movements',
                'Painful menstrual bleeding','Erectile dysfunction',
                'Behavioural changes','Hallucinations','COVID toes',
                'Partial face or body paralysis')
symp_names <- gsub(" ", "_", symp_names)

severity_labels <- c("Asymptomatic", "Mild", "Moderate", "Severe")

## ---- helper functions --------------------------------------------------

## post-infection protein scores with age, sex and severity
get_scores <- function(timepoint = "post", vars = c("age", "sex", "cov_severity_a")) {
  s <- read.csv2(f_scores)
  s <- subset(s, pre_post_infection == timepoint)
  s$cov_prot_score_ensemble <- as.numeric(s$cov_prot_score_ensemble)
  merge(s, d_cati[c("sid", vars)], by = "sid", all.x = TRUE)
}

## signature protein levels (post-infection), columns renamed to gene symbols
get_proteins <- function(timepoint = "post") {
  p <- read.csv2(f_prot)[, 1:53]
  if (!is.null(timepoint)) p <- subset(p, pre_post_infection == timepoint)
  colnames(p)[5:53] <- read.csv2(f_sig)$gene_symbol
  p
}

## recode symptom variables: 99 and any non-binary value to NA
recode_symptoms <- function(df, vars) {
  for (v in vars) {
    df[[v]][df[[v]] == 99] <- NA
    df[[v]][!df[[v]] %in% c(0, 1) & !is.na(df[[v]])] <- NA
    df[[v]] <- as.numeric(df[[v]])
  }
  df
}

## linear models of the protein score on a set of binary predictors
score_on <- function(df, predictors, covariates = c("age", "sex")) {
  res <- lapply(predictors, function(v) {
    m <- lm(reformulate(c(v, covariates), response = "cov_prot_score_ensemble"),
            data = df)
    co <- summary(m)$coefficients[2, ]
    data.frame(predictor = v, beta = co[1], se = co[2], p = co[4],
               n = nobs(m), row.names = NULL)
  })
  bind_rows(res) %>%
    mutate(LCI = beta - 1.96 * se,
           UCI = beta + 1.96 * se)
}

## prevalence ratio from Poisson regression with robust (HC0) standard errors
prevalence_ratio <- function(outcome, data, covariates = c("age", "sex"), label = outcome) {
  form <- reformulate(c("cov_prot_score_ensemble", covariates), response = outcome)
  m    <- glm(form, data = data, family = poisson(link = "log"))
  se   <- sqrt(diag(vcovHC(m, type = "HC0")))
  est  <- coef(m)[2]; s <- se[2]
  data.frame(outcome = label,
             PR  = exp(est),
             LCI = exp(est - 1.96 * s),
             UCI = exp(est + 1.96 * s),
             p   = 2 * (1 - pnorm(abs(est / s))),
             n   = nobs(m), row.names = NULL)
}


################################################################################
## 1  PROTEIN SIGNATURE: RIDGE COEFFICIENTS (Figure 2A)
################################################################################

sig <- read.csv2(f_sig)
sig$Log_Odds_Ratio <- log(sig$OR.per.SD)
sig <- sig[order(sig$Log_Odds_Ratio, decreasing = FALSE), ]

signature_panel <- ggplot(sig, aes(x = Log_Odds_Ratio,
                                   y = factor(gene_name, levels = sig$gene_name),
                                   fill = Direction)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("+" = "#c98fa4", "-" = "#8eb56e")) +
  xlab("Log Odds Ratio") + ylab("Gene name") +
  theme_classic() +
  theme(axis.text  = element_text(size = 14),
        axis.title = element_text(size = 18),
        legend.position = "none",
        panel.grid = element_blank())

ggsave(file.path(out, "signature_log_odds_ratios.png"), signature_panel,
       width = 13, height = 11)


################################################################################
## 2  PROTEIN SCORE BY INFECTION SEVERITY
##    post-infection: Supplementary Figure 3A; pre-infection: Figure 3B
################################################################################

plot_score_by_severity <- function(timepoint, adjusted, title, fills, filename) {

  dat <- na.omit(get_scores(timepoint))
  dat$sex <- as.factor(dat$sex)

  if (adjusted) {
    dat$y    <- residuals(lm(cov_prot_score_ensemble ~ age + sex, data = dat))
    y_label  <- "Protein Score (adjusted for age and sex)"
  } else {
    dat$y    <- dat$cov_prot_score_ensemble
    y_label  <- "Protein Score for Symptomatic Infection"
  }

  dat$cov_severity_a <- factor(dat$cov_severity_a, levels = 1:4, labels = severity_labels)

  p <- ggplot(dat, aes(x = cov_severity_a, y = y,
                       color = cov_severity_a, fill = cov_severity_a)) +
    sm_raincloud() +
    xlab("Course of SARS-CoV-2 Infection") + ylab(y_label) +
    scale_fill_manual(values = fills) +
    scale_color_manual(values = fills) +
    stat_summary(fun.data = function(x) data.frame(y = Inf, label = paste0("n = ", length(x))),
                 geom = "text", hjust = 0.5, vjust = 1.5, col = "#262222", size = 5) +
    ylim(-3.5, 5.5) +
    theme_bw() +
    theme(axis.text  = element_text(size = 18),
          axis.title = element_text(size = 18),
          plot.title = element_text(size = 20),
          panel.grid.minor = element_blank()) +
    stat_compare_means(comparisons = list(c("Asymptomatic", "Mild"),
                                          c("Asymptomatic", "Moderate"),
                                          c("Asymptomatic", "Severe")),
                       method = "t.test", position = "identity", vjust = -0.2, size = 5) +
    ggtitle(title) +
    guides(fill = "none", color = "none")

  ggsave(file.path(out, filename), p, width = 12, height = 8)
  p
}

post_cols <- c("#cf5389", "#6379c2", "#a581c9", "#94ad5c")
pre_cols  <- c("#9c3e86", "#1e8fa6", "#714694", "#499c3e")

plot_score_by_severity("post", FALSE, "Post-infection samples", post_cols,
                       "post_infection_severity.png")
plot_score_by_severity("post", TRUE,  "Post-infection samples", post_cols,
                       "post_infection_severity_adjusted.png")
plot_score_by_severity("pre",  FALSE, "Pre-infection samples",  pre_cols,
                       "pre_infection_severity.png")
plot_score_by_severity("pre",  TRUE,  "Pre-infection samples",  pre_cols,
                       "pre_infection_severity_adjusted.png")


################################################################################
## 3  PRE- AND POST-INFECTION PROTEIN DYNAMICS
##    paired t-tests in the 189 individuals sampled at both timepoints
##    Supplementary Figure 4
################################################################################

prot_all <- get_proteins(timepoint = NULL)
gene_names <- read.csv2(f_sig)$gene_symbol

paired_sids <- intersect(prot_all$sid[prot_all$pre_post_infection == "pre"],
                         prot_all$sid[prot_all$pre_post_infection == "post"])
prot_paired <- prot_all %>% filter(sid %in% paired_sids)

pre_df  <- prot_paired %>% filter(pre_post_infection == "pre")  %>% arrange(sid)
post_df <- prot_paired %>% filter(pre_post_infection == "post") %>% arrange(sid)
stopifnot(all(pre_df$sid == post_df$sid))

paired_res <- data.frame(
  gene_name = gene_names,
  p_value   = sapply(gene_names, function(g) t.test(pre_df[[g]], post_df[[g]], paired = TRUE)$p.value),
  row.names = NULL
)
paired_res$fdr <- p.adjust(paired_res$p_value, method = "fdr")
paired_res <- paired_res[order(paired_res$p_value), ]

sig_proteins <- paired_res %>% filter(fdr < 0.05)
cat("Proteins tested:", nrow(paired_res),
    "| FDR < 0.05:", nrow(sig_proteins), "\n")

plot_df <- prot_paired %>%
  select(sid, pre_post_infection, all_of(sig_proteins$gene_name)) %>%
  pivot_longer(cols = all_of(sig_proteins$gene_name),
               names_to = "gene_name", values_to = "expression") %>%
  mutate(pre_post_infection = factor(pre_post_infection, levels = c("pre", "post")))

fdr_lab <- sig_proteins %>% mutate(lab = paste0("FDR p = ", signif(fdr, 2)))

pre_post_proteins <- ggplot(plot_df, aes(x = pre_post_infection, y = expression,
                                         fill = pre_post_infection, color = pre_post_infection)) +
  geom_violin(trim = FALSE, alpha = 1, size = 1, show.legend = FALSE) +
  geom_boxplot(width = 0.25, outlier.shape = NA, size = 1) +
  geom_beeswarm(dodge.width = 1, alpha = 0.5, cex = 1.5, priority = "ascending") +
  scale_fill_manual(values = c("pre" = "#f5dbb3", "post" = "#f5cbdf")) +
  scale_color_manual(values = c("pre" = "#db8704", "post" = "#d11d73")) +
  scale_x_discrete(labels = c("pre" = "Pre-infection", "post" = "Post-infection")) +
  xlab("") + ylab("Protein Expression") +
  facet_wrap(~ gene_name, scales = "free_y") +
  geom_text(data = fdr_lab, aes(x = 1.5, y = Inf, label = lab),
            inherit.aes = FALSE, vjust = 1.8, size = 5, colour = "black") +
  theme_bw() +
  theme(axis.text = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 15),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

ggsave(file.path(out, "pre_post_proteins.png"), pre_post_proteins,
       width = 18, height = 7, dpi = 150)


################################################################################
## 4  PROTEIN SCORE AND ACUTE SYMPTOMS
##    4a  protein score ~ symptom (Figure 3A)
##    4b  individual proteins ~ symptom (Figure 4, heatmap)
################################################################################

## ---- 4a  protein score ~ symptom ---------------------------------------

symp <- merge(get_scores("post", c("age", "sex")),
              d_cati[c("sid", symp_vars)], by = "sid", all.x = TRUE)
names(symp)[match(symp_vars, names(symp))] <- symp_names
symp <- recode_symptoms(symp, symp_names)

## symptoms reported by at least 5% of the post-infection sample
counts <- sapply(symp[symp_names], function(x) sum(x == 1, na.rm = TRUE))
reduced_symps <- symp_names[counts >= 24]

res_symptoms <- score_on(symp, reduced_symps)

forplo(res_symptoms[, c("beta", "LCI", "UCI")],
       row.labels = gsub("_", " ", res_symptoms$predictor),
       em = "\u03b2", linreg = TRUE,
       col = "black", ci.edge = FALSE, left.bar = FALSE,
       left.align = TRUE, sort = TRUE, size = 1.5,
       margin.left = 10, margin.right = 25, column.spacing = 4,
       save = TRUE, save.path = out,
       save.name = "forplo_symptoms_score", save.type = "png",
       save.width = 10, save.height = 10)

## ---- 4b  individual proteins ~ symptom ---------------------------------

prot_post <- get_proteins("post")
df_prot_symp <- merge(prot_post[, c("sid", gene_names)],
                      symp[, c("sid", "age", "sex", reduced_symps)],
                      by = "sid")

estmat  <- matrix(NA, nrow = length(gene_names), ncol = length(reduced_symps),
                  dimnames = list(gene_names, reduced_symps))
pvalmat <- estmat

for (i in seq_along(gene_names)) {
  for (j in seq_along(reduced_symps)) {
    m  <- glm(reformulate(c(paste0("scale(`", gene_names[i], "`)"), "age", "sex"),
                          response = reduced_symps[j]),
              family = binomial, data = df_prot_symp)
    co <- summary(m)$coefficients[2, ]
    estmat[i, j]  <- co[1]
    pvalmat[i, j] <- co[4]
  }
}

## signed -log10(p), capped at 10 and set to zero below p = 0.05
signed_sig <- pmin(-log10(pvalmat), 10)
signed_sig[signed_sig < 1.3] <- 0
signed_sig <- signed_sig * sign(estmat)
colnames(signed_sig) <- gsub("_", " ", colnames(signed_sig))

heatmap_matrix <- t(signed_sig)

## label proteins and symptoms with many associations
sig_threshold <- 1.3
proteins_per_symptom <- colSums(abs(signed_sig) > sig_threshold, na.rm = TRUE)
symptoms_per_protein <- rowSums(abs(signed_sig) > sig_threshold, na.rm = TRUE)
cut_sym <- 10; cut_prot <- 10

row_hl <- proteins_per_symptom[rownames(heatmap_matrix)] >= cut_sym
col_hl <- symptoms_per_protein[colnames(heatmap_matrix)] >= cut_prot

row_labels <- ifelse(row_hl,
                     paste0(rownames(heatmap_matrix), " (",
                            proteins_per_symptom[rownames(heatmap_matrix)], ")"),
                     rownames(heatmap_matrix))
col_labels <- ifelse(col_hl,
                     paste0(colnames(heatmap_matrix), " (",
                            symptoms_per_protein[colnames(heatmap_matrix)], ")"),
                     colnames(heatmap_matrix))

ht <- Heatmap(
  heatmap_matrix,
  name = "-log10(p)",
  col  = colorRamp2(c(-3, 0, 3), c("#0439b5", "white", "#027302")),
  rect_gp = gpar(col = "grey", lwd = 0.5),
  row_labels = row_labels, column_labels = col_labels,
  row_names_gp = gpar(fontsize = 8, fontface = ifelse(row_hl, "bold", "plain"),
                      col = ifelse(row_hl, "#2c8f40", "black")),
  column_names_gp = gpar(fontsize = 8, fontface = ifelse(col_hl, "bold", "plain"),
                         col = ifelse(col_hl, "#2c8f40", "black"), rot = 60),
  row_gap = unit(2, "mm"), column_gap = unit(2, "mm")
)

png(file.path(out, "heatmap_symptoms_proteins.png"),
    width = 2400, height = 1800, res = 300)
draw(ht, heatmap_legend_side = "right")
dev.off()


################################################################################
## 5  SYMPTOM PREVALENCE THRESHOLD: SENSITIVITY ANALYSIS
##    5% (>= 24 individuals) vs 10% (>= 50 individuals)
##    Supplementary Figure 6
################################################################################

run_threshold <- function(min_n) {
  score_on(symp, symp_names[counts >= min_n]) %>%
    mutate(n_cases = counts[predictor])
}

thr5  <- 24
thr10 <- ceiling(0.10 * nrow(symp))

res5  <- run_threshold(thr5)
res10 <- run_threshold(thr10)

cat("Symptoms at 5% threshold: ", nrow(res5),
    "| at 10% threshold:", nrow(res10), "\n")
cat("p < 0.05 at 5%: ", sum(res5$p  < 0.05), "of", nrow(res5),  "\n")
cat("p < 0.05 at 10%:", sum(res10$p < 0.05), "of", nrow(res10), "\n")

dropped <- setdiff(res5$predictor, res10$predictor)
cat("Symptoms dropped at the 10% threshold:", length(dropped),
    "| of these with p < 0.05:", sum(res5$p[res5$predictor %in% dropped] < 0.05), "\n")

plot_dat <- res5 %>%
  mutate(retained = ifelse(n_cases >= thr10, 1, 2),
         label = gsub("_", " ", predictor)) %>%
  arrange(desc(beta))

forplo(plot_dat[, c("beta", "LCI", "UCI")],
       row.labels = plot_dat$label,
       em = "\u03b2", linreg = TRUE,
       ci.edge = FALSE, left.bar = FALSE, right.bar = FALSE,
       left.align = TRUE, sort = FALSE, size = 2,
       margin.left = 12, margin.right = 22, column.spacing = 4,
       fill.by = plot_dat$retained,
       fill.colors = c("black", "#c1121f"),
       fill.labs = c("Retained at 10% threshold",
                     "Excluded at 10% threshold, included at 5%"),
       legend = TRUE, legend.spacing = 2, legend.vadj = 1, legend.hadj = -2.5,
       save = TRUE, save.path = out,
       save.name = "forplo_symptom_threshold_sensitivity",
       save.type = "png", save.width = 10, save.height = 10)


################################################################################
## 6  PROTEIN SCORE AND POST-COVID SEQUELAE (Figure 3B)
##    prevalence ratios from Poisson regression with robust standard errors
################################################################################

covid_data <- read.csv(f_covid)

pcs <- merge(get_scores("post", c("age", "sex")),
             covid_data[c("ID", "Sequelae.nach.3.Monaten", "Sequelae.nach.6.Monaten")],
             by.x = "sid", by.y = "ID", all.x = TRUE)

pcs$pcs3 <- ifelse(pcs$Sequelae.nach.3.Monaten == "yes", 1,
                   ifelse(pcs$Sequelae.nach.3.Monaten == "no", 0, NA))
pcs$pcs6 <- ifelse(pcs$Sequelae.nach.6.Monaten == "yes", 1,
                   ifelse(pcs$Sequelae.nach.6.Monaten == "no", 0, NA))

table(pcs$pcs3, useNA = "ifany")
table(pcs$pcs6, useNA = "ifany")

res_pcs <- rbind(
  prevalence_ratio("pcs3", pcs, label = "Post-COVID sequelae > 3 months"),
  prevalence_ratio("pcs6", pcs, label = "Post-COVID sequelae > 6 months")
)
print(res_pcs)

forplo(res_pcs[, c("PR", "LCI", "UCI")],
       row.labels = res_pcs$outcome,
       em = "PR", linreg = FALSE,
       col = "#80644f", ci.edge = FALSE,
       left.bar = FALSE, right.bar = FALSE, left.align = TRUE,
       xlim = c(0.9, 1.5), margin.left = 12, margin.right = 10,
       legend = FALSE,
       save = TRUE, save.path = out,
       save.name = "forplo_post_covid_sequelae", save.type = "png",
       save.width = 7, save.height = 4)


################################################################################
## 7  PROTEIN SCORE AND LONG-TERM SYMPTOMS (> 6 months)
##    Supplementary Figure 9
################################################################################

seq_data <- read.csv2(f_seq)
long_symp <- merge(get_scores("post", c("age", "sex")), seq_data,
                   by.x = "sid", by.y = "X", all.x = TRUE)
names(long_symp)[9:69] <- symp_names
long_symp <- recode_symptoms(long_symp, symp_names)

long_counts <- sapply(long_symp[symp_names], function(x) sum(x == 1, na.rm = TRUE))
long_reduced <- symp_names[long_counts >= 24]

res_long <- lapply(long_reduced, function(s) {
  m  <- glm(reformulate(c("cov_prot_score_ensemble", "age", "sex"), response = s),
            data = long_symp, family = binomial)
  co <- summary(m)$coefficients[2, ]
  data.frame(symptom = s, OR = exp(co[1]),
             LCI = exp(co[1] - 1.96 * co[2]),
             UCI = exp(co[1] + 1.96 * co[2]),
             p = co[4], row.names = NULL)
}) %>% bind_rows()

forplo(res_long[, c("OR", "LCI", "UCI")],
       row.labels = gsub("_", " ", res_long$symptom),
       col = "black", ci.edge = FALSE, left.bar = FALSE,
       left.align = TRUE, sort = TRUE, size = 1.5,
       margin.left = 10, margin.right = 25, column.spacing = 4,
       save = TRUE, save.path = out,
       save.name = "forplo_long_term_symptoms", save.type = "png",
       save.width = 10, save.height = 10)


################################################################################
## 8  PROTEIN SCORE AND CARDIOVASCULAR RISK FACTORS
##    (the analysis of lipid markers measured before infection, Figure 5,
##     was performed separately and is not part of this script)
################################################################################

cvrf_vars <- c("hyper", "diab", "dyslip", "cancer", "mi", "adipos", "nic",
               "stroke", "afib", "cad", "chf", "ckd", "cld", "copd", "cvd",
               "pad", "pulmd", "vte")

cvrf <- merge(get_scores("post", c("age", "sex")),
              d[c("sid", cvrf_vars)], by = "sid", all.x = TRUE)
cvrf[cvrf_vars] <- lapply(cvrf[cvrf_vars], function(x) {
  x[x == 99] <- NA
  ifelse(x == "yes", 1, ifelse(is.na(x), NA, 0))
})

res_cvrf <- score_on(cvrf, cvrf_vars)

cvrf_labels <- c("Arterial Hypertension", "Diabetes Mellitus Type 2", "Dyslipidemia",
                 "History of Cancer", "History of Myocardial Infarction", "Obesity",
                 "Smoking", "History of Stroke", "Atrial Fibrillation",
                 "Coronary Artery Disease", "Chronic Heart Failure",
                 "Chronic Kidney Disease", "Chronic Liver Disease",
                 "Chronic Obstructive Pulmonary Disease", "Cardiovascular Disease",
                 "Peripheral Artery Disease", "Pulmonary Disease",
                 "History of Venous Thromboembolism")

forplo(res_cvrf[, c("beta", "LCI", "UCI")],
       row.labels = cvrf_labels,
       em = "\u03b2", linreg = TRUE,
       col = "black", ci.edge = FALSE, left.bar = FALSE,
       left.align = TRUE, sort = TRUE,
       margin.left = 10, margin.right = 25, column.spacing = 4,
       save = TRUE, save.path = out,
       save.name = "forplo_score_cvrf", save.type = "png",
       save.width = 12, save.height = 6)

################################################################################
## 9  PER-PROTEIN STATISTICS ACROSS SEVERITY GROUPS
##    group means relative to asymptomatic (Supplementary Figure 5)
##    per-protein table (Supplementary Figure 10)
################################################################################

ens <- merge(get_proteins("post"), d_cati[c("sid", "cov_severity_a")],
             by = "sid", all.x = TRUE)

dat <- ens %>%
  filter(!is.na(cov_severity_a)) %>%
  mutate(severity = factor(cov_severity_a, levels = 1:4, labels = severity_labels))

X  <- apply(as.matrix(dat[, gene_names]), 2, as.numeric)
Xz <- scale(X)

cat("N with severity category:", nrow(dat), "\n"); print(table(dat$severity))

## ---- group means, expressed as the difference from asymptomatic --------

grp_mean <- as.data.frame(Xz) %>%
  mutate(severity = dat$severity) %>%
  group_by(severity) %>%
  summarise(across(all_of(gene_names), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

asymp_row <- grp_mean %>% filter(severity == "Asymptomatic") %>%
  select(all_of(gene_names)) %>% as.numeric()

mat_diff <- t(sweep(as.matrix(grp_mean[, gene_names]), 2, asymp_row, "-"))
colnames(mat_diff) <- as.character(grp_mean$severity)
mat_diff <- mat_diff[, colnames(mat_diff) != "Asymptomatic", drop = FALSE]

ht_diff <- Heatmap(
  mat_diff, name = "Z-score",
  col = colorRamp2(c(-0.8, 0, 0.8), c("#1d3557", "white", "#e63946")),
  cluster_columns = FALSE, cluster_rows = TRUE,
  column_names_rot = 45, column_names_gp = gpar(fontsize = 9),
  column_title = "Course of SARS-CoV-2 infection", column_title_side = "bottom",
  row_title = "Signature proteins", row_names_gp = gpar(fontsize = 8),
  width = unit(3.6, "cm"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", mat_diff[i, j]), x, y, gp = gpar(fontsize = 6.5, col = "grey20"))
  }
)

png(file.path(out, "heatmap_group_means_vs_asymptomatic.png"),
    width = 1400, height = 2600, res = 300)
draw(ht_diff, heatmap_legend_side = "right")
dev.off()

## ---- per-protein comparisons with BH adjustment -------------------------

long_df <- dat %>%
  mutate(severity_group = factor(case_when(
    cov_severity_a == 1 ~ "Asymptomatic",
    cov_severity_a == 2 ~ "Mild",
    cov_severity_a %in% c(3, 4) ~ "Moderate/Severe"),
    levels = c("Asymptomatic", "Mild", "Moderate/Severe"))) %>%
  pivot_longer(cols = all_of(gene_names), names_to = "gene", values_to = "expression") %>%
  group_by(gene) %>%
  mutate(expression_z = as.numeric(scale(expression))) %>%
  ungroup()

stats_tab <- long_df %>%
  group_by(gene) %>%
  summarise(
    n_asymp  = sum(severity_group == "Asymptomatic"    & !is.na(expression)),
    n_mild   = sum(severity_group == "Mild"            & !is.na(expression)),
    n_modsev = sum(severity_group == "Moderate/Severe" & !is.na(expression)),
    smd_mild   = mean(expression_z[severity_group == "Mild"], na.rm = TRUE) -
                 mean(expression_z[severity_group == "Asymptomatic"], na.rm = TRUE),
    smd_modsev = mean(expression_z[severity_group == "Moderate/Severe"], na.rm = TRUE) -
                 mean(expression_z[severity_group == "Asymptomatic"], na.rm = TRUE),
    p_mild   = t.test(expression[severity_group == "Asymptomatic"],
                      expression[severity_group == "Mild"])$p.value,
    p_modsev = t.test(expression[severity_group == "Asymptomatic"],
                      expression[severity_group == "Moderate/Severe"])$p.value,
    .groups = "drop") %>%
  mutate(padj_mild   = p.adjust(p_mild,   method = "BH"),
         padj_modsev = p.adjust(p_modsev, method = "BH"),
         min_padj    = pmin(padj_mild, padj_modsev, na.rm = TRUE))

cat("Proteins with BH-adjusted p < 0.05 in at least one comparison:",
    sum(stats_tab$min_padj < 0.05, na.rm = TRUE), "\n")

write.csv(stats_tab, file.path(out, "per_protein_statistics.csv"), row.names = FALSE)


################################################################################
## 10  REVISION ANALYSES
################################################################################

## -----------------------------------------------------------------------
## 10a  Serologically confirmed infections
##      nucleocapsid antibodies above the assay cut-off
##      (Abbott IgG >= 1.4 RLU or Roche >= 0.8 U/ml, baseline or follow-up)
## -----------------------------------------------------------------------

sero <- merge(get_scores("post"), d[c("sid", "igg", "roche", "igg_fu", "roche_fu")],
              by = "sid", all.x = TRUE)

sero$any_measured <- !is.na(sero$igg) | !is.na(sero$roche) |
                     !is.na(sero$igg_fu) | !is.na(sero$roche_fu)

sero$sero_pos <- (sero$igg      >= 1.4) %in% TRUE |
                 (sero$roche    >= 0.8) %in% TRUE |
                 (sero$igg_fu   >= 1.4) %in% TRUE |
                 (sero$roche_fu >= 0.8) %in% TRUE

table(sero$any_measured)
table(sero$sero_pos, useNA = "ifany")
table(sero$cov_severity_a[sero$sero_pos], useNA = "ifany")

sero_sub <- subset(sero, sero_pos & !is.na(cov_severity_a))
sero_sub$symptomatic <- ifelse(sero_sub$cov_severity_a == 1, 0, 1)

cat("N in serologically confirmed subgroup:", nrow(sero_sub), "\n")
print(table(sero_sub$cov_severity_a))

summary(lm(cov_prot_score_ensemble ~ symptomatic + age + sex, data = sero_sub))
summary(lm(cov_prot_score_ensemble ~ factor(cov_severity_a) + age + sex, data = sero_sub))


## -----------------------------------------------------------------------
## 10b  Baseline characteristics: individuals with vs without proteomics
## -----------------------------------------------------------------------

prot_ids <- unique(read.csv2(f_scores)$sid)

base_vars <- c("age", "sex", "bmi", "nic", "hyper", "diab", "dyslip",
               "adipos", "cancer", "mi", "stroke", "cvd")

base <- d[!duplicated(d$sid), c("sid", base_vars)]
base$group <- ifelse(base$sid %in% prot_ids, "proteomics", "no proteomics")
table(base$group)

bin_vars <- c("nic", "hyper", "diab", "dyslip", "adipos", "cancer", "mi", "stroke", "cvd")
base[bin_vars] <- lapply(base[bin_vars], function(x) {
  x[x %in% c(99, "unknown", "")] <- NA
  ifelse(x %in% c("yes", 1, "1"), 1, ifelse(is.na(x), NA, 0))
})
base$female <- ifelse(base$sex == "Women", 1, 0)

row_cont <- function(v, label) {
  a <- base[[v]][base$group == "proteomics"]
  b <- base[[v]][base$group == "no proteomics"]
  data.frame(Characteristic = label,
             Proteomics    = sprintf("%.1f \u00b1 %.1f", mean(a, na.rm = TRUE), sd(a, na.rm = TRUE)),
             No_proteomics = sprintf("%.1f \u00b1 %.1f", mean(b, na.rm = TRUE), sd(b, na.rm = TRUE)),
             p = signif(t.test(a, b)$p.value, 3))
}

row_bin <- function(v, label) {
  pct <- function(g) sprintf("%.1f (%d)",
                             100 * mean(base[[v]][base$group == g], na.rm = TRUE),
                             sum(base[[v]][base$group == g] == 1, na.rm = TRUE))
  data.frame(Characteristic = label,
             Proteomics = pct("proteomics"),
             No_proteomics = pct("no proteomics"),
             p = signif(chisq.test(table(base$group, base[[v]]))$p.value, 3))
}

baseline_table <- rbind(
  row_cont("age", "Age [years]"),
  row_bin("female", "Female sex"),
  row_cont("bmi", "BMI [kg/m2]"),
  row_bin("nic",    "Smoking"),
  row_bin("hyper",  "Arterial hypertension"),
  row_bin("dyslip", "Dyslipidemia"),
  row_bin("adipos", "Obesity"),
  row_bin("diab",   "Diabetes mellitus type 2"),
  row_bin("cvd",    "Cardiovascular disease"),
  row_bin("cancer", "History of cancer"),
  row_bin("mi",     "History of myocardial infarction"),
  row_bin("stroke", "History of stroke")
)

print(baseline_table)
write.xlsx(baseline_table, file.path(out, "baseline_comparison.xlsx"), rowNames = FALSE)


## -----------------------------------------------------------------------
## 10c  Export of the analysis samples
##      (shared with the statisticians for the additional revision analyses)
## -----------------------------------------------------------------------

post_export <- merge(get_scores("post"),
                     covid_data[c("ID", "Sequelae.nach.3.Monaten", "Sequelae.nach.6.Monaten")],
                     by.x = "sid", by.y = "ID", all.x = TRUE)
post_export$symptomatic <- ifelse(post_export$cov_severity_a == 1, 0,
                                  ifelse(post_export$cov_severity_a %in% 2:4, 1, NA))
post_export$pcs_3m <- ifelse(post_export$Sequelae.nach.3.Monaten == "yes", 1, 0)
post_export$pcs_6m <- ifelse(post_export$Sequelae.nach.6.Monaten == "yes", 1, 0)
post_export <- post_export[, c("sid", "cov_prot_score_ensemble", "age", "sex",
                               "cov_severity_a", "symptomatic", "pcs_3m", "pcs_6m")]

pre_export <- get_scores("pre")
pre_export$symptomatic <- ifelse(pre_export$cov_severity_a == 1, 0,
                                 ifelse(pre_export$cov_severity_a %in% 2:4, 1, NA))
pre_export <- pre_export[, c("sid", "cov_prot_score_ensemble", "age", "sex",
                             "cov_severity_a", "symptomatic")]

write.xlsx(list(post_infection = post_export, pre_infection = pre_export),
           file.path(out, "analysis_samples.xlsx"), rowNames = FALSE)

################################################################################
## End of script
################################################################################
