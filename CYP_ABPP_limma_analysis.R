## ===============================
## Setup
## ===============================

setwd("data/CYP/males")

library(tidyverse)
library(limma)
library(ggrepel)
library(gridExtra)
library(writexl)
library(imputeLCMD)
library(impute)

analysis_name <- "8weeks_males_CYP"

plot_dir <- "QC_and_Results_Plots"
dir.create(plot_dir, showWarnings = FALSE)

## ===============================
## Load data
## ===============================

prot <- read_tsv("protein_crosstab.tsv", show_col_types = FALSE)

meta <- read_tsv("metadata.tsv", show_col_types = FALSE) %>%
  transmute(
    sample_id = sample,
    treatment_group = .data[["group"]],
    probe = treatment
  ) %>%
  distinct(sample_id, .keep_all = TRUE) %>%
  mutate(
    treatment_group = factor(treatment_group, levels = c("control", "dosed")),
    probe = factor(probe, levels = c("negative", "positive")),
    group_probe = factor(
      paste(treatment_group, probe, sep = "_"),
      levels = c(
        "control_negative",
        "control_positive",
        "dosed_negative",
        "dosed_positive"
      )
    )
  )

meta_ordered <- meta %>%
  arrange(group_probe, sample_id)

annot <- prot %>%
  select(Protein, Function, Gene, Peptides, `Peptide Number`, `Unique Peptide(s)`)

expr <- prot %>%
  select(Protein, all_of(meta$sample_id)) %>%
  column_to_rownames("Protein") %>%
  as.data.frame()

expr[] <- lapply(expr, function(x) as.numeric(as.character(x)))
expr <- as.matrix(expr)

expr[expr == 0] <- NA

expr_log <- log2(expr)
expr_log[!is.finite(expr_log)] <- NA

meta <- meta %>%
  arrange(match(sample_id, colnames(expr_log)))

stopifnot(all(meta$sample_id == colnames(expr_log)))

sample_order <- meta_ordered$sample_id

## ===============================
## Raw QC plots
## ===============================

expr_long_raw <- expr_log %>%
  as.data.frame() %>%
  rownames_to_column("Protein") %>%
  pivot_longer(-Protein, names_to = "Sample", values_to = "Intensity") %>%
  left_join(meta, by = c("Sample" = "sample_id")) %>%
  mutate(Sample = factor(Sample, levels = sample_order))

p_raw_density <- ggplot(expr_long_raw, aes(Intensity, color = Sample)) +
  geom_density(na.rm = TRUE) +
  theme_minimal() +
  labs(title = "Density — raw log2")

p_raw_box <- ggplot(expr_long_raw, aes(Sample, Intensity, fill = group_probe)) +
  geom_boxplot(outlier.size = 0.3, na.rm = TRUE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot — raw log2")

## ===============================
## Median normalization within group_probe
## ===============================

expr_norm <- expr_log

for (gp in levels(meta$group_probe)) {
  
  gp_samples <- meta$sample_id[meta$group_probe == gp]
  
  gp_medians <- apply(
    expr_norm[, gp_samples, drop = FALSE],
    2,
    median,
    na.rm = TRUE
  )
  
  gp_grand_median <- median(gp_medians, na.rm = TRUE)
  
  expr_norm[, gp_samples] <- sweep(
    expr_norm[, gp_samples, drop = FALSE],
    2,
    gp_medians - gp_grand_median,
    FUN = "-"
  )
}

## ===============================
## Normalized-only long data
## ===============================

expr_long_norm_only <- expr_norm %>%
  as.data.frame() %>%
  rownames_to_column("Protein") %>%
  pivot_longer(-Protein, names_to = "Sample", values_to = "Intensity") %>%
  left_join(meta, by = c("Sample" = "sample_id")) %>%
  mutate(Sample = factor(Sample, levels = sample_order))

## ===============================
## Hybrid imputation
## ===============================

group_sample_list <- split(meta$sample_id, meta$group_probe)

detected_prop_by_group <- sapply(group_sample_list, function(samples) {
  rowMeans(!is.na(expr_norm[, samples, drop = FALSE]))
}) %>%
  as.data.frame()

shared_proteins <- apply(detected_prop_by_group >= 0.50, 1, all)
unique_proteins <- !shared_proteins

cat("\n==============================\n")
cat("Hybrid imputation summary\n")
cat("==============================\n")
cat("KNN-safe shared proteins:        ", sum(shared_proteins), "\n")
cat("Sparse/unique proteins MinDet:   ", sum(unique_proteins), "\n")
cat("==============================\n\n")

expr_shared <- expr_norm[shared_proteins, , drop = FALSE]

if (nrow(expr_shared) > 0) {
  for (gp in names(group_sample_list)) {
    
    gp_samples <- group_sample_list[[gp]]
    
    cat("KNN imputation for shared proteins |", gp, "\n")
    
    expr_shared[, gp_samples] <-
      impute.knn(
        as.matrix(expr_shared[, gp_samples, drop = FALSE]),
        k = 2
      )$data
  }
}

expr_unique <- expr_norm[unique_proteins, , drop = FALSE]

if (nrow(expr_unique) > 0) {
  
  expr_unique <- as.matrix(expr_unique)
  
  for (gp in names(group_sample_list)) {
    
    gp_samples <- group_sample_list[[gp]]
    
    cat("MinDet imputation for sparse/unique proteins |", gp, "\n")
    
    expr_unique[, gp_samples] <-
      impute.MinDet(
        as.matrix(expr_unique[, gp_samples, drop = FALSE]),
        q = 0.05
      )
  }
}

expr_imp <- rbind(expr_shared, expr_unique)
expr_imp <- expr_imp[rownames(expr_norm), ]
mode(expr_imp) <- "numeric"

cat("\nRemaining missing values after imputation:", sum(is.na(expr_imp)), "\n")

if (sum(is.na(expr_imp)) > 0) {
  stop("There are still missing values after imputation.")
}

saveRDS(
  expr_imp,
  file = file.path(plot_dir, paste0("expr_imp_", analysis_name, ".rds"))
)

write_xlsx(
  expr_imp %>%
    as.data.frame() %>%
    rownames_to_column("Protein"),
  file.path(plot_dir, paste0("expr_imp_", analysis_name, ".xlsx"))
)

## ===============================
## Limma design for enrichment-gate testing
## ===============================

design <- model.matrix(~0 + group_probe, data = meta)
colnames(design) <- levels(meta$group_probe)

## First fit all imputed proteins to test probe enrichment
fit_gate <- lmFit(expr_imp, design)

contrast_gate <- makeContrasts(
  control_probe_enrichment = control_positive - control_negative,
  dosed_probe_enrichment   = dosed_positive - dosed_negative,
  levels = design
)

fit_gate2 <- contrasts.fit(fit_gate, contrast_gate)
fit_gate2 <- eBayes(fit_gate2, trend = TRUE, robust = TRUE)

control_gate <- topTable(
  fit_gate2,
  coef = "control_probe_enrichment",
  number = Inf,
  adjust.method = "BH"
) %>%
  rownames_to_column("Protein") %>%
  transmute(
    Protein,
    control_probe_logFC = logFC,
    control_probe_p_value = P.Value,
    control_probe_adj_P = adj.P.Val,
    control_probe_t = t
  )

dosed_gate <- topTable(
  fit_gate2,
  coef = "dosed_probe_enrichment",
  number = Inf,
  adjust.method = "BH"
) %>%
  rownames_to_column("Protein") %>%
  transmute(
    Protein,
    dosed_probe_logFC = logFC,
    dosed_probe_p_value = P.Value,
    dosed_probe_adj_P = adj.P.Val,
    dosed_probe_t = t
  )

probe_enrichment_df <- control_gate %>%
  full_join(dosed_gate, by = "Protein") %>%
  mutate(
    keep_probe_enriched =
      control_probe_logFC > 0 &
      control_probe_p_value < 0.05 &
      dosed_probe_logFC > 0 &
      dosed_probe_p_value < 0.05,
    
    probe_enrichment_class = case_when(
      control_probe_logFC > 0 & control_probe_p_value < 0.05 &
        dosed_probe_logFC > 0 & dosed_probe_p_value < 0.05 ~
        "Significantly enriched in control and dosed",
      
      control_probe_logFC > 0 & control_probe_p_value < 0.05 &
        !(dosed_probe_logFC > 0 & dosed_probe_p_value < 0.05) ~
        "Significantly enriched only in control",
      
      !(control_probe_logFC > 0 & control_probe_p_value < 0.05) &
        dosed_probe_logFC > 0 & dosed_probe_p_value < 0.05 ~
        "Significantly enriched only in dosed",
      
      control_probe_logFC > 0 & dosed_probe_logFC > 0 ~
        "Positive enrichment in both but not significant in both",
      
      TRUE ~ "Not consistently probe-enriched"
    )
  ) %>%
  left_join(annot, by = "Protein")

write_xlsx(
  probe_enrichment_df,
  file.path(plot_dir, paste0("probe_enrichment_", analysis_name, ".xlsx"))
)

proteins_keep <- probe_enrichment_df %>%
  filter(keep_probe_enriched) %>%
  pull(Protein)

expr_imp_filtered <- expr_imp[proteins_keep, , drop = FALSE]

cat("\n==============================\n")
cat("Probe enrichment raw p-value filtering\n")
cat("==============================\n")
cat("Proteins before enrichment filter: ", nrow(expr_imp), "\n")
cat("Proteins kept for limma:           ", nrow(expr_imp_filtered), "\n")
cat("Proteins removed:                  ", nrow(expr_imp) - nrow(expr_imp_filtered), "\n")
cat("==============================\n\n")

print(table(probe_enrichment_df$probe_enrichment_class))

if (nrow(expr_imp_filtered) == 0) {
  stop("No proteins passed the raw p-value probe-enrichment gate.")
}

saveRDS(
  expr_imp_filtered,
  file = file.path(plot_dir, paste0("expr_imp_probe_enriched_", analysis_name, ".rds"))
)

write_xlsx(
  expr_imp_filtered %>%
    as.data.frame() %>%
    rownames_to_column("Protein"),
  file.path(plot_dir, paste0("expr_imp_probe_enriched_", analysis_name, ".xlsx"))
)

## ===============================
## Normalized QC plots
## ===============================

expr_long_norm <- expr_imp_filtered %>%
  as.data.frame() %>%
  rownames_to_column("Protein") %>%
  pivot_longer(-Protein, names_to = "Sample", values_to = "Intensity") %>%
  left_join(meta, by = c("Sample" = "sample_id")) %>%
  mutate(Sample = factor(Sample, levels = sample_order))

p_norm_density <- ggplot(expr_long_norm_only, aes(Intensity, color = Sample)) +
  geom_density(na.rm = TRUE) +
  theme_minimal() +
  labs(title = "Density — normalized log2")

p_norm_box <- ggplot(expr_long_norm_only, aes(Sample, Intensity, fill = group_probe)) +
  geom_boxplot(outlier.size = 0.3, na.rm = TRUE) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot — normalized log2")

## ===============================
## Final limma interaction test
## Only significantly probe-enriched proteins enter this model
## ===============================

contrast_final <- makeContrasts(
  probe_adjusted_dosed_vs_control =
    (dosed_positive - dosed_negative) -
    (control_positive - control_negative),
  levels = design
)

fit_final <- lmFit(expr_imp_filtered, design)
fit_final2 <- contrasts.fit(fit_final, contrast_final)
fit_final2 <- eBayes(fit_final2, trend = TRUE, robust = TRUE)

png(
  file.path(plot_dir, paste0("mean_variance_limma_", analysis_name, ".png")),
  width = 800,
  height = 600,
  res = 120
)
plotSA(fit_final2, main = "Mean–variance trend")
dev.off()

## ===============================
## PCA
## ===============================

pca <- prcomp(t(expr_imp_filtered))

pca_df <- data.frame(
  Sample = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
) %>%
  left_join(meta, by = c("Sample" = "sample_id"))

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = treatment_group, shape = probe)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA — significantly probe-enriched proteins only")

ggsave(
  file.path(plot_dir, paste0("PCA_probe_enriched_rawP_only_", analysis_name, ".png")),
  p_pca,
  width = 7,
  height = 6,
  dpi = 300
)

## ===============================
## Final results
## ===============================

res <- topTable(
  fit_final2,
  coef = "probe_adjusted_dosed_vs_control",
  number = Inf,
  adjust.method = "BH"
) %>%
  rownames_to_column("Protein") %>%
  left_join(annot, by = "Protein") %>%
  left_join(
    probe_enrichment_df %>%
      select(
        Protein,
        control_probe_logFC,
        control_probe_p_value,
        control_probe_adj_P,
        dosed_probe_logFC,
        dosed_probe_p_value,
        dosed_probe_adj_P,
        keep_probe_enriched,
        probe_enrichment_class
      ),
    by = "Protein"
  ) %>%
  mutate(
    Gene = as.character(Gene),
    Gene_clean = if_else(is.na(Gene) | Gene == "", Protein, Gene),
    is_cyp_gene = str_detect(Gene_clean, regex("cyp", ignore_case = TRUE)),
    Significance = case_when(
      adj.P.Val < 0.05 & logFC > 2  ~ "Up",
      adj.P.Val < 0.05 & logFC < -2 ~ "Down",
      TRUE ~ "NS"
    ),
    neg_log10_adj_P = -log10(adj.P.Val)
  )

volcano <- res

## ===============================
## Volcano plot
## ===============================

p_volcano <- ggplot(volcano, aes(logFC, neg_log10_adj_P, color = Significance)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c(Up = "red", Down = "blue", NS = "grey")) +
  theme_minimal() +
  geom_text_repel(
    data = volcano %>%
      filter(
        is_cyp_gene,
        adj.P.Val < 0.05,
        abs(logFC) > 2
      ),
    aes(label = Gene_clean),
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  labs(
    title = "Volcano: Significantly probe-enriched proteins only",
    subtitle = "Gate: positive enrichment with raw p < 0.05 in both control and dosed",
    x = "Probe-adjusted interaction log2FC",
    y = "-log10(adj P-value)"
  )

ggsave(
  file.path(plot_dir, paste0("volcano_probe_enriched_", analysis_name, ".png")),
  p_volcano,
  width = 8,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(plot_dir, paste0("volcano_probe_enriched_", analysis_name, ".pdf")),
  p_volcano,
  width = 8,
  height = 7
)

## ===============================
## QC comparison plots
## ===============================

pdf(
  file.path(plot_dir, paste0("density_comparison_", analysis_name, ".pdf")),
  width = 12,
  height = 6
)

grid.arrange(
  p_raw_density + labs(title = "Density — raw log2"),
  p_norm_density + labs(title = "Density — normalized"),
  ncol = 2
)

dev.off()

pdf(
  file.path(plot_dir, paste0("boxplot_comparison_", analysis_name, ".pdf")),
  width = 12,
  height = 6
)

grid.arrange(
  p_raw_box + labs(title = "Boxplot — raw log2"),
  p_norm_box + labs(title = "Boxplot — normalized"),
  ncol = 2
)

dev.off()

## ===============================
## Save all-protein, CYP-only, and non-CYP Excel files
## ===============================

all_proteins_export <- res %>%
  transmute(
    Protein,
    Gene = Gene_clean,
    Peptides,
    `Peptide Number`,
    `Unique Peptide(s)`,
    Function,
    
    control_probe_logFC,
    control_probe_p_value,
    control_probe_adj_P,
    dosed_probe_logFC,
    dosed_probe_p_value,
    dosed_probe_adj_P,
    keep_probe_enriched,
    probe_enrichment_class,
    
    interaction_log2FC = logFC,
    AveExpr,
    t_statistic = t,
    p_value = P.Value,
    adj_P_value = adj.P.Val,
    B,
    neg_log10_adj_P,
    is_cyp_gene,
    Significance
  ) %>%
  arrange(adj_P_value)

cyp_only_export <- all_proteins_export %>%
  filter(is_cyp_gene) %>%
  arrange(adj_P_value)

non_cyp_export <- all_proteins_export %>%
  filter(!is_cyp_gene) %>%
  arrange(adj_P_value)

write_xlsx(
  all_proteins_export,
  path = file.path(plot_dir, paste0("ALL_proteins_", analysis_name, ".xlsx"))
)

write_xlsx(
  cyp_only_export,
  path = file.path(plot_dir, paste0("CYP_only_", analysis_name, ".xlsx"))
)

write_xlsx(
  non_cyp_export,
  path = file.path(plot_dir, paste0("Others_only_", analysis_name, ".xlsx"))
)

## ===============================
## CYP-only volcano
## ===============================

p_cyp_volcano <- ggplot(
  volcano %>% filter(is_cyp_gene),
  aes(logFC, neg_log10_adj_P, color = Significance)
) +
  geom_point(alpha = 0.9, size = 2.8) +
  scale_color_manual(values = c(Up = "red", Down = "blue", NS = "grey")) +
  theme_minimal() +
  geom_text_repel(
    aes(label = Gene_clean),
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3
  ) +
  labs(
    title = "CYP-only volcano: Significantly probe-enriched proteins only",
    subtitle = "Gate: positive enrichment with raw p < 0.05 in both control and dosed",
    x = "Probe-adjusted interaction log2FC",
    y = "-log10(adj P-value)"
  )

ggsave(
  file.path(plot_dir, paste0("CYP_only_", analysis_name, ".png")),
  p_cyp_volcano,
  width = 8,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(plot_dir, paste0("CYP_only_", analysis_name, ".pdf")),
  p_cyp_volcano,
  width = 8,
  height = 7
)

##Code for making paneled volcano plot across time

## ===============================
## CYP volcano plots
## Grouped panel + individual plots
## Same axis limits across all timepoints
## 5W | 8W | 12W females
## ===============================

setwd("##Insert_File_Location##")

library(tidyverse)
library(readxl)
library(ggrepel)

## ===============================
## Input files
## ===============================

files <- tibble(
  Timepoint = c("5W", "8W", "12W"),
  File = c(
    "CYP_only_5wks_females_CYP.xlsx",
    "CYP_only_8wks_females_CYP.xlsx",
    "CYP_only_12wks_females_CYP.xlsx"
  )
)

## ===============================
## Read and combine CYP files
## ===============================

cyp_all <- files %>%
  mutate(data = map(File, read_excel)) %>%
  unnest(data) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("5W", "8W", "12W")),
    Gene = as.character(Gene),
    interaction_log2FC = as.numeric(interaction_log2FC),
    adj_P_value = as.numeric(adj_P_value),
    neg_log10_adj_P = -log10(adj_P_value),
    
    Volcano_Status = case_when(
      adj_P_value < 0.05 & interaction_log2FC > 2 ~ "Up",
      adj_P_value < 0.05 & interaction_log2FC < -2 ~ "Down",
      TRUE ~ "NS"
    )
  )

## ===============================
## Safety check
## ===============================

cat("\nRows per timepoint:\n")
print(table(cyp_all$Timepoint))

cat("\nVolcano status counts:\n")
print(table(cyp_all$Timepoint, cyp_all$Volcano_Status))

## ===============================
## Output directory
## ===============================

out_dir <- file.path(getwd(), "CYP_volcano_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ===============================
## Global axis limits
## Same x/y scale for grouped and individual plots
## ===============================

x_limits <- range(cyp_all$interaction_log2FC, na.rm = TRUE)
y_limits <- range(cyp_all$neg_log10_adj_P, na.rm = TRUE)

x_limits <- c(
  floor(x_limits[1]) - 0.5,
  ceiling(x_limits[2]) + 0.5
)

y_limits <- c(
  0,
  ceiling(y_limits[2]) + 0.5
)

cat("\nGlobal x-axis limits:\n")
print(x_limits)

cat("\nGlobal y-axis limits:\n")
print(y_limits)

## ===============================
## Shared volcano plotting function
## ===============================

make_cyp_volcano <- function(df, plot_title, plot_subtitle, use_facet = FALSE) {
  
  p <- ggplot(
    df,
    aes(
      x = interaction_log2FC,
      y = neg_log10_adj_P
    )
  ) +
    geom_point(
      aes(color = Volcano_Status),
      size = 2.6,
      alpha = 0.85
    ) +
    geom_text_repel(
      data = df %>% filter(Volcano_Status != "NS"),
      aes(label = Gene),
      size = 3.0,
      fontface = "bold",
      max.overlaps = Inf,
      box.padding = 1.2,
      point.padding = 0.55,
      min.segment.length = 0,
      segment.color = "grey50",
      segment.linewidth = 0.25,
      force = 15,
      force_pull = 0.05,
      max.time = 5,
      max.iter = 50000,
      direction = "both",
      xlim = x_limits,
      ylim = y_limits
    ) +
    geom_vline(
      xintercept = c(-0.5, 0.5),
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.4
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dotted",
      color = "grey45",
      linewidth = 0.4
    ) +
    scale_color_manual(
      values = c(
        "Down" = "steelblue",
        "NS" = "grey75",
        "Up" = "firebrick"
      )
    ) +
    coord_cartesian(
      xlim = x_limits,
      ylim = y_limits
    ) +
    labs(
      x = "log2FC (Dosed vs Control)",
      y = "-log10 adjusted p-value"
    ) +
    theme_bw(base_size = 13) +
    theme(
      legend.position = "none",
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "grey85", color = "grey30"),
      panel.grid.minor = element_blank()
    )
  
  if (use_facet) {
    p <- p + facet_wrap(~ Timepoint, nrow = 1)
  }
  
  return(p)
}

## ===============================
## 1. Grouped paneled plot
## ===============================

cyp_paneled_volcano <- make_cyp_volcano(
  df = cyp_all,
  plot_title = "CYP protein differential abundance | females",
  plot_subtitle = "PAH-dosed vs Control | 5W, 8W, and 12W",
  use_facet = TRUE
)

print(cyp_paneled_volcano)

ggsave(
  filename = file.path(out_dir, "CYP_paneled_volcano_females_5W_8W_12W.png"),
  plot = cyp_paneled_volcano,
  width = 14,
  height = 5.5,
  dpi = 300
)

ggsave(
  filename = file.path(out_dir, "CYP_paneled_volcano_females_5W_8W_12W.tiff"),
  plot = cyp_paneled_volcano,
  width = 14,
  height = 5.5,
  dpi = 300,
  compression = "lzw"
)

## ===============================
## 2. Individual volcano plots
## Same grey strip style as grouped plot
## Same axis limits across all timepoints
## ===============================

for (tp in levels(cyp_all$Timepoint)) {
  
  df_tp <- cyp_all %>%
    filter(Timepoint == tp) %>%
    mutate(Timepoint = factor(Timepoint, levels = c("5W", "8W", "12W")))
  
  p_tp <- make_cyp_volcano(
    df = df_tp,
    plot_title = "CYP protein differential abundance | females",
    plot_subtitle = "PAH-dosed vs Control",
    use_facet = TRUE
  )
  
  print(p_tp)
  
  ggsave(
    filename = file.path(out_dir, paste0("CYP_volcano_females_", tp, ".png")),
    plot = p_tp,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(out_dir, paste0("CYP_volcano_females_", tp, ".tiff")),
    plot = p_tp,
    width = 7,
    height = 6,
    dpi = 300,
    compression = "lzw"
  )
}
