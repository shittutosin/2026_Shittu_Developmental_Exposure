## ===============================
## Setup
## ===============================

setwd("data/proteomics/limma")

library(tidyverse)
library(limma)
library(impute)
library(ggrepel)
library(gridExtra)
library(writexl)
library(imputeLCMD)   # for MinDet


plot_dir <- "QC_and_Results_Plots"
dir.create(plot_dir, showWarnings = FALSE)

## ===============================
## Load data
## ===============================

prot <- read_tsv("protein_crosstab.tsv")
meta <- read_tsv("metadata.tsv")

meta_ordered <- meta %>% arrange(group, sample_id)

annot <- prot %>%
  select(Protein, Function, Gene, Peptides, `Peptide Number`, `Unique Peptide(s)`)

expr <- prot %>%
  select(Protein, all_of(meta$sample_id)) %>%
  column_to_rownames("Protein") %>%
  as.data.frame()

# Force numeric (VERY IMPORTANT)
expr[] <- lapply(expr, function(x) as.numeric(as.character(x)))

# Convert to matrix
expr <- as.matrix(expr)

# Replace zeros, then log2
expr[expr == 0] <- NA
expr_log <- log2(expr)

# Clean non-finite values
expr_log[!is.finite(expr_log)] <- NA

meta <- meta %>%
  mutate(group = factor(group, levels = c("control", "dosed"))) %>%
  arrange(match(sample_id, colnames(expr_log)))

stopifnot(all(meta$sample_id == colnames(expr_log)))

# Define sample order ONCE
sample_order <- c(
  meta$sample_id[meta$group == "control"],
  meta$sample_id[meta$group == "dosed"]
)

# Raw long data
expr_long_raw <- expr_log %>%
  as.data.frame() %>%
  pivot_longer(everything(),
               names_to = "Sample",
               values_to = "Intensity") %>%
  left_join(meta, by = c("Sample" = "sample_id")) %>%
  mutate(Sample = factor(Sample, levels = sample_order))


## ===============================
## Raw QC plots
## ===============================

p_raw_density <- ggplot(expr_long_raw, aes(Intensity, color = Sample)) +
  geom_density(na.rm = TRUE) +
  theme_minimal() +
  labs(title = "Density — raw log2")

p_raw_box <- ggplot(expr_long_raw, aes(Sample, Intensity, fill = group)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot — raw log2")


## ===============================
## Median normalization (per sample, within group)
## ===============================

ctrl_samples  <- meta$sample_id[meta$group == "control"]
dosed_samples <- meta$sample_id[meta$group == "dosed"]

expr_norm <- expr_log

## ---- Control group ----
ctrl_medians <- apply(expr_norm[, ctrl_samples, drop = FALSE],
                      2, median, na.rm = TRUE)
ctrl_grand_median <- median(ctrl_medians, na.rm = TRUE)

expr_norm[, ctrl_samples] <-
  sweep(expr_norm[, ctrl_samples, drop = FALSE],
        2,
        ctrl_medians - ctrl_grand_median,
        FUN = "-")

## ---- Dosed group ----
dosed_medians <- apply(expr_norm[, dosed_samples, drop = FALSE],
                       2, median, na.rm = TRUE)
dosed_grand_median <- median(dosed_medians, na.rm = TRUE)

expr_norm[, dosed_samples] <-
  sweep(expr_norm[, dosed_samples, drop = FALSE],
        2,
        dosed_medians - dosed_grand_median,
        FUN = "-")


## ===============================
## Identify shared proteins
## ===============================

detected_ctrl  <- rowSums(!is.na(expr_norm[, ctrl_samples])) > 0
detected_dosed <- rowSums(!is.na(expr_norm[, dosed_samples])) > 0

shared_proteins <- detected_ctrl & detected_dosed
unique_proteins <- !shared_proteins

expr_long_norm_only <- expr_norm %>%
  as.data.frame() %>%
  pivot_longer(everything(),
               names_to = "Sample",
               values_to = "Intensity") %>%
  left_join(meta_ordered, by = c("Sample" = "sample_id")) %>%
  mutate(Sample = factor(Sample, levels = meta_ordered$sample_id))


## ===============================
## KNN imputation on NORMALIZED data
## ===============================

expr_shared <- expr_norm[shared_proteins, , drop = FALSE]

expr_shared[, ctrl_samples] <-
  impute.knn(as.matrix(expr_shared[, ctrl_samples]), k = 2)$data

expr_shared[, dosed_samples] <-
  impute.knn(as.matrix(expr_shared[, dosed_samples]), k = 2)$data


## ===============================
## MinDet imputation for UNIQUE proteins
## ===============================

expr_unique <- expr_norm[unique_proteins, , drop = FALSE]
expr_unique <- as.matrix(expr_unique)

## ---- Control group ----
expr_unique_ctrl <- expr_unique[, ctrl_samples, drop = FALSE]

expr_unique_ctrl_imp <- impute.MinDet(
  expr_unique_ctrl,
  q = 0.05
)

## ---- Dosed group ----
expr_unique_dosed <- expr_unique[, dosed_samples, drop = FALSE]

expr_unique_dosed_imp <- impute.MinDet(
  expr_unique_dosed,
  q = 0.05
)

## ---- Recombine ----
expr_unique_imp <- expr_unique

expr_unique_imp[, ctrl_samples]  <- expr_unique_ctrl_imp
expr_unique_imp[, dosed_samples] <- expr_unique_dosed_imp


## ===============================
## Combine imputed data
## ===============================

expr_imp <- expr_norm
expr_imp[shared_proteins, ] <- expr_shared
expr_imp[unique_proteins, ] <- expr_unique_imp

mode(expr_imp) <- "numeric"





## ===============================
## Save normalized + imputed matrices
## ===============================

saveRDS(expr_imp,
        file = file.path(plot_dir, "expr_imp_12wks_F.rds"))

write_xlsx(
  expr_imp %>% 
    as.data.frame() %>% 
    rownames_to_column("Protein"),
  file.path(plot_dir, "expr_imp_12wks_F.xlsx")
)



## ===============================
## Long data after imputation
## ===============================

expr_long_raw <- expr_log %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "Intensity") %>%
  left_join(meta_ordered, by = c("Sample" = "sample_id")) %>%
  mutate(Sample = factor(Sample, levels = meta_ordered$sample_id))

expr_long_norm <- expr_imp %>%
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "Sample", values_to = "Intensity") %>%
  left_join(meta_ordered, by = c("Sample" = "sample_id")) %>%
  mutate(Sample = factor(Sample, levels = meta_ordered$sample_id))


## ===============================
## Normalized QC plots
## ===============================

p_norm_density <- ggplot(expr_long_norm_only,
                         aes(Intensity, color = Sample)) +
  geom_density(na.rm = TRUE) +
  theme_minimal() +
  labs(title = "Density — raw vs normalized")

p_norm_box <- ggplot(expr_long_norm_only,
                     aes(Sample, Intensity, fill = group)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Boxplot — raw vs normalized")


## ===============================
## Limma
## ===============================

design <- model.matrix(~0 + meta$group)
colnames(design) <- levels(meta$group)

contrast <- makeContrasts(
  dosed_vs_control = dosed - control,
  levels = design
)

fit <- lmFit(expr_imp, design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)

png(file.path(plot_dir, "mean_variance_limma_12wks_F.png"),
    width = 800, height = 600, res = 120)
plotSA(fit2, main = "Mean–variance trend (limma)")
dev.off()


## ===============================
## PCA
## ===============================

pca <- prcomp(t(expr_imp))
pca_df <- data.frame(Sample = rownames(pca$x),
                     PC1 = pca$x[,1],
                     PC2 = pca$x[,2]) %>%
  left_join(meta, by = c("Sample" = "sample_id"))

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(title = "PCA — normalized & imputed")

ggsave(file.path(plot_dir, "PCA_normalized_imputed_12wks_F.png"),
       p_pca, width = 7, height = 6, dpi = 300)


## ===============================
## Results + volcano
## ===============================

res <- topTable(fit2, coef = "dosed_vs_control",
                number = Inf, adjust.method = "BH") %>%
  rownames_to_column("Protein") %>%
  left_join(annot, by = "Protein")

volcano <- res %>%
  mutate(signif = case_when(
    adj.P.Val < 0.05 & logFC > 2  ~ "Up",
    adj.P.Val < 0.05 & logFC < -2 ~ "Down",
    TRUE                         ~ "NS"
  ))

p_volcano <- ggplot(volcano, aes(logFC, -log10(adj.P.Val), color = signif)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c(Up = "red", Down = "blue", NS = "grey")) +
  theme_minimal() +
  geom_text_repel(
    data = subset(volcano, adj.P.Val < 0.05 & abs(logFC) > 2),
    aes(label = Protein),
    max.overlaps = 10
  ) +
  labs(title = "Volcano: Dosed vs Control",
       x = "log2 Fold Change",
       y = "-log10(adj P-value)")

ggsave(file.path(plot_dir, "volcano_dosed_vs_control_12wks_F.png"),
       p_volcano, width = 8, height = 7, dpi = 300)


## ===============================
## QC comparison PDFs
## ===============================

pdf(file.path(plot_dir, "density_comparison_12wks_F.pdf"), width = 12, height = 6)

grid.arrange(
  p_raw_density + labs(title = "Density — raw log2"),
  p_norm_density + labs(title = "Density — normalized"),
  ncol = 2
)

dev.off()

pdf(file.path(plot_dir, "boxplot_comparison_12wks_F.pdf"), width = 12, height = 6)

grid.arrange(
  p_raw_box + labs(title = "Boxplot — raw log2"),
  p_norm_box + labs(title = "Boxplot — normalized"),
  ncol = 2
)

dev.off()


## ===============================
## Save volcano-ready data to Excel
## ===============================

volcano_export <- res %>%
  transmute(
    Protein,
    Peptides,
    Gene,
    Function,
    log2FC = logFC,
    t_statistic = t,
    p_value = P.Value,
    adj_P_value = adj.P.Val,
    neg_log10_adj_P = -log10(adj.P.Val),
    Significance = case_when(
      adj.P.Val < 0.05 & logFC > 2  ~ "Up",
      adj.P.Val < 0.05 & logFC < -2 ~ "Down",
      TRUE                           ~ "NS"
    )
  ) %>%
  arrange(adj_P_value)

write_xlsx(
  volcano_export,
  path = file.path(plot_dir, "volcano_ready_data_dosed_vs_control_12wks_F.xlsx")
)


## ===============================
## Summary statistics
## ===============================

n_total_proteins <- nrow(expr)
n_tested <- nrow(res)

n_up <- sum(volcano$signif == "Up", na.rm = TRUE)
n_down <- sum(volcano$signif == "Down", na.rm = TRUE)
n_ns <- sum(volcano$signif == "NS", na.rm = TRUE)

cat("\n==============================\n")
cat("Proteomics differential expression summary\n")
cat("==============================\n")
cat("Total proteins quantified:      ", n_total_proteins, "\n")
cat("Proteins tested by limma:       ", n_tested, "\n")
cat("Significantly increased (Up):   ", n_up, "\n")
cat("Significantly decreased (Down): ", n_down, "\n")
cat("Not significant:               ", n_ns, "\n")
cat("==============================\n\n")

table(volcano$signif, unique_proteins[match(volcano$Protein, rownames(expr_norm))])


## ===============================
## Supplementary information workbook
## ===============================

supp_file <- file.path(plot_dir, "Supplementary_Proteomics_Processing_12wks_F.xlsx")

raw_qc_summary <- expr_log %>%
  as.data.frame() %>%
  summarise(across(
    everything(),
    list(
      median = ~ median(.x, na.rm = TRUE),
      mean = ~ mean(.x, na.rm = TRUE),
      sd = ~ sd(.x, na.rm = TRUE),
      missing_values = ~ sum(is.na(.x)),
      detected_proteins = ~ sum(!is.na(.x))
    )
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("Sample", ".value"),
    names_pattern = "(.+)_(median|mean|sd|missing_values|detected_proteins)"
  ) %>%
  left_join(meta, by = c("Sample" = "sample_id"))

norm_qc_summary <- expr_norm %>%
  as.data.frame() %>%
  summarise(across(
    everything(),
    list(
      median = ~ median(.x, na.rm = TRUE),
      mean = ~ mean(.x, na.rm = TRUE),
      sd = ~ sd(.x, na.rm = TRUE),
      missing_values = ~ sum(is.na(.x)),
      detected_proteins = ~ sum(!is.na(.x))
    )
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("Sample", ".value"),
    names_pattern = "(.+)_(median|mean|sd|missing_values|detected_proteins)"
  ) %>%
  left_join(meta, by = c("Sample" = "sample_id"))

normalization_factors <- tibble(
  Sample = c(names(ctrl_medians), names(dosed_medians)),
  group = c(rep("control", length(ctrl_medians)),
            rep("dosed", length(dosed_medians))),
  sample_median = c(ctrl_medians, dosed_medians),
  group_grand_median = c(
    rep(ctrl_grand_median, length(ctrl_medians)),
    rep(dosed_grand_median, length(dosed_medians))
  ),
  median_shift_applied = sample_median - group_grand_median
)

protein_detection_status <- tibble(
  Protein = rownames(expr_norm),
  detected_control = detected_ctrl,
  detected_dosed = detected_dosed,
  detection_status = case_when(
    detected_control & detected_dosed ~ "Shared",
    detected_control & !detected_dosed ~ "Control_only",
    !detected_control & detected_dosed ~ "Dosed_only",
    TRUE ~ "Not_detected"
  )
)

imputation_summary <- tibble(
  Sample = colnames(expr_norm),
  missing_before_imputation = colSums(is.na(expr_norm)),
  missing_after_imputation = colSums(is.na(expr_imp))
) %>%
  left_join(meta, by = c("Sample" = "sample_id"))

supplementary_summary <- tibble(
  Metric = c(
    "Total proteins quantified",
    "Shared proteins",
    "Unique proteins",
    "Control-only proteins",
    "Dosed-only proteins",
    "Proteins tested by limma",
    "Significantly increased",
    "Significantly decreased",
    "Not significant"
  ),
  Value = c(
    nrow(expr),
    sum(shared_proteins),
    sum(unique_proteins),
    sum(protein_detection_status$detection_status == "Control_only"),
    sum(protein_detection_status$detection_status == "Dosed_only"),
    nrow(res),
    n_up,
    n_down,
    n_ns
  )
)

write_xlsx(
  list(
    Metadata = meta,
    Raw_Log2_Matrix = expr_log %>%
      as.data.frame() %>%
      rownames_to_column("Protein"),
    Raw_QC_Summary = raw_qc_summary,
    Normalization_Factors = normalization_factors,
    Normalized_Matrix = expr_norm %>%
      as.data.frame() %>%
      rownames_to_column("Protein"),
    Normalized_QC_Summary = norm_qc_summary,
    Protein_Detection_Status = protein_detection_status,
    Shared_Proteins_KNN_Imputed = expr_shared %>%
      as.data.frame() %>%
      rownames_to_column("Protein"),
    Unique_Proteins_MinDet_Imputed = expr_unique_imp %>%
      as.data.frame() %>%
      rownames_to_column("Protein"),
    Final_Imputed_Matrix = expr_imp %>%
      as.data.frame() %>%
      rownames_to_column("Protein"),
    Limma_Results_All = res,
    Volcano_Ready_Data = volcano_export,
    Differential_Summary = supplementary_summary,
    Imputation_Summary = imputation_summary
  ),
  path = supp_file
)