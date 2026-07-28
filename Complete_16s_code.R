###############################################################################
# 16S ANALYSIS PIPELINE
# Beta diversity: relative abundance
# Alpha diversity: rarefied counts
# Organism-level DA: relative abundance
# Taxa filter: within-treatment abundance filter
# Targeted single taxon report: Limosilactobacillus rodentium only
# Presence/absence export: all taxa + taxa similar to L. rodentium
###############################################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(tidyverse)
  library(RColorBrewer)
  library(writexl)
})

###############################################################################
# USER SETTINGS
###############################################################################

base_dir <- "data/16s_analysis/rel-abundance_files"

timepoint_label <- "8weeks_16s_relative_abundance_files"
sex_keep <- "M"
tissue_label <- "Colon"

metadata_file <- "Metadata.txt"
read_summary_file <- "read_summary.txt"

out_dir <- file.path(base_dir, "16S_outputs_files")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

setwd(base_dir)

###############################################################################
# LOAD INPUT FILES
###############################################################################

files <- list.files(base_dir, pattern = "_rel-abundance.tsv$", full.names = TRUE)
samples <- gsub("_rel-abundance.tsv", "", basename(files))

meta <- read.table(metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
read_counts <- read.table(read_summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

cat("\nFiles loaded:\n")
print(basename(files))

cat("\nMetadata samples:\n")
print(meta$SampleID)

###############################################################################
# BUILD LONG TABLE
###############################################################################

long_df <- map2_dfr(files, samples, function(file, samp) {
  read_tsv(file, col_types = cols()) %>%
    dplyr::select(tax_id, lineage, abundance, `estimated counts`) %>%
    dplyr::rename(est_counts = `estimated counts`) %>%
    dplyr::mutate(sample = samp)
})

long_df <- long_df %>%
  filter(!tax_id %in% c("mapped_unclassified", "unmapped")) %>%
  mutate(
    abundance = replace_na(abundance, 0),
    est_counts = replace_na(est_counts, 0)
  ) %>%
  group_by(tax_id, sample) %>%
  summarise(
    abundance = abundance[1],
    est_counts = est_counts[1],
    lineage = lineage[1],
    .groups = "drop"
  )
clean_collapsed_reads <- long_df %>%
  left_join(
    meta %>% select(SampleID, Sex, Treatment),
    by = c("sample" = "SampleID")
  ) %>%
  mutate(
    Timepoint = timepoint_label,
    Tissue = tissue_label
  ) %>%
  select(
    sample,
    Sex,
    Treatment,
    tax_id,
    lineage,
    abundance,
    est_counts,
    Timepoint,
    Tissue
  )
###############################################################################
# COUNT MATRIX
###############################################################################

otu_counts <- long_df %>%
  dplyr::select(tax_id, sample, est_counts) %>%
  pivot_wider(
    names_from = sample,
    values_from = est_counts,
    values_fill = 0
  ) %>%
  arrange(tax_id)

count_mat <- as.matrix(otu_counts[, -1])
rownames(count_mat) <- otu_counts$tax_id

###############################################################################
# TAXONOMY TABLE
###############################################################################

ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

tax_mat <- long_df %>%
  dplyr::select(tax_id, lineage) %>%
  distinct() %>%
  mutate(lineage_split = strsplit(lineage, ";")) %>%
  rowwise() %>%
  mutate(
    lineage_split = list({
      pieces <- lineage_split[lineage_split != ""]
      pieces <- trimws(pieces)
      vals <- tail(pieces, 7)
      if (length(vals) < 7) vals <- c(rep(NA, 7 - length(vals)), vals)
      rev(vals)
    })
  ) %>%
  ungroup() %>%
  unnest_wider(lineage_split, names_sep = "_") %>%
  rename_with(~ ranks, starts_with("lineage_split_")) %>%
  mutate(across(all_of(ranks), trimws))

tax_mat <- tax_mat %>%
  mutate(
    Genus = na_if(Genus, ""),
    Species = na_if(Species, ""),
    Species = ifelse(
      is.na(Species) |
        Species == "" |
        Species == Genus |
        Species == word(Genus, 1) |
        Species %in% c("sp", "sp.", "bacterium", "uncultured bacterium"),
      paste(Genus, "sp."),
      Species
    )
  )

tax_mat <- tax_mat %>%
  group_by(Species) %>%
  mutate(
    Species = ifelse(
      n() == 1,
      Species,
      paste0(Species, "_tax", row_number())
    )
  ) %>%
  ungroup()

tax_mat <- tax_mat %>%
  column_to_rownames("tax_id")

tax_mat <- tax_mat[, ranks, drop = FALSE]

###############################################################################
# SUBSET METADATA
###############################################################################

keep_samples <- meta %>%
  filter(
    Sex == sex_keep,
    Treatment %in% c("Control", "Dosed"),
    SampleID %in% colnames(count_mat)
  ) %>%
  pull(SampleID)

count_mat_sub <- count_mat[, keep_samples, drop = FALSE]

meta_ps <- meta %>%
  filter(SampleID %in% keep_samples) %>%
  arrange(match(SampleID, keep_samples)) %>%
  column_to_rownames("SampleID")

cat("\nSamples used:\n")
print(keep_samples)

cat("\nTreatment counts:\n")
print(table(meta_ps$Treatment))

###############################################################################
# BUILD PHYLOSEQ OBJECT
###############################################################################

OTU <- otu_table(count_mat_sub, taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax_mat[rownames(count_mat_sub), , drop = FALSE]))
META <- sample_data(meta_ps)

ps <- phyloseq(OTU, TAX, META)

###############################################################################
# RELATIVE ABUNDANCE OBJECT
###############################################################################

ps_rel <- transform_sample_counts(ps, function(x) {
  if (sum(x) == 0) x else x / sum(x)
})

###############################################################################
# RELATIVE ABUNDANCE LONG TABLE
###############################################################################

rel_long_taxid <- as.data.frame(otu_table(ps_rel))

if (!taxa_are_rows(ps_rel)) {
  rel_long_taxid <- t(rel_long_taxid)
}

rel_long_taxid <- rel_long_taxid %>%
  as.data.frame() %>%
  rownames_to_column("tax_id") %>%
  pivot_longer(
    cols = -tax_id,
    names_to = "SampleID",
    values_to = "RelAbund"
  ) %>%
  left_join(
    meta_ps %>% rownames_to_column("SampleID"),
    by = "SampleID"
  ) %>%
  left_join(
    tax_mat %>% rownames_to_column("tax_id"),
    by = "tax_id"
  )

###############################################################################
# WITHIN-TREATMENT TAXA FILTER FOR DA
###############################################################################

###############################################################################
# WITHIN-TREATMENT TAXA FILTER FOR DA
# Keep taxa if:
# Mean relative abundance is >= 0.1% in Control OR Dosed
###############################################################################

taxa_filter_by_treatment <- rel_long_taxid %>%
  group_by(tax_id, Species, Treatment) %>%
  summarise(
    mean_rel_abund = mean(RelAbund, na.rm = TRUE),
    mean_rel_abund_percent = mean_rel_abund * 100,
    n_samples_in_treatment = n(),
    n_present_in_treatment = sum(RelAbund > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    treatment_passes_0.1_percent_filter = mean_rel_abund >= 0.001
  )

taxa_filter_summary <- taxa_filter_by_treatment %>%
  group_by(tax_id, Species) %>%
  summarise(
    control_mean_rel_abund = mean_rel_abund[Treatment == "Control"],
    dosed_mean_rel_abund = mean_rel_abund[Treatment == "Dosed"],
    control_mean_percent = mean_rel_abund_percent[Treatment == "Control"],
    dosed_mean_percent = mean_rel_abund_percent[Treatment == "Dosed"],
    control_passes_filter = treatment_passes_0.1_percent_filter[Treatment == "Control"],
    dosed_passes_filter = treatment_passes_0.1_percent_filter[Treatment == "Dosed"],
    passed_filter = any(treatment_passes_0.1_percent_filter, na.rm = TRUE),
    .groups = "drop"
  )

taxa_keep <- taxa_filter_summary %>%
  filter(passed_filter) %>%
  pull(tax_id)

cat("\nNumber of taxa retained for DA using within-treatment abundance filter only:\n")
print(length(taxa_keep))

cat("\nL. rodentium filter status:\n")
print(
  taxa_filter_summary %>%
    filter(str_detect(Species, regex("rodentium", ignore_case = TRUE)))
)
###############################################################################
# BETA DIVERSITY — BRAY-CURTIS ON RELATIVE ABUNDANCE
###############################################################################

bray_dist <- phyloseq::distance(ps_rel, method = "bray")
ord_bray <- ordinate(ps_rel, method = "PCoA", distance = bray_dist)

bray_df <- plot_ordination(ps_rel, ord_bray, type = "samples", justDF = TRUE) %>%
  rownames_to_column("SampleID") %>%
  left_join(
    meta_ps %>% rownames_to_column("SampleID"),
    by = "SampleID"
  ) %>%
  mutate(
    Treatment = Treatment.x,
    Sex = Sex.x
  ) %>%
  dplyr::select(SampleID, Axis.1, Axis.2, Sex, Treatment)

eig <- ord_bray$values$Eigenvalues
var_explained <- eig / sum(eig)

axis1_pct <- round(var_explained[1] * 100, 1)
axis2_pct <- round(var_explained[2] * 100, 1)

set.seed(123)

permanova <- adonis2(
  bray_dist ~ Treatment,
  data = meta_ps,
  permutations = 9999
)

r2 <- permanova$R2[1]
pval <- permanova$`Pr(>F)`[1]
label_text <- paste0(
  "p = ",
  signif(pval, 3)
)
p_beta <- ggplot(bray_df, aes(Axis.1, Axis.2, color = Treatment)) +
  stat_ellipse(type = "t", linewidth = 1, level = 0.85) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c(
    "Control" = "steelblue",
    "Dosed" = "firebrick"
  )) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = label_text,
    hjust = 1.1,
    vjust = 1.5,
    size = 5,
    fontface = "bold"
  ) +
  labs(
    x = paste0("Axis 1 (", axis1_pct, "%)"),
    y = paste0("Axis 2 (", axis2_pct, "%)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    
    legend.title = element_text(size = 15, face = "bold"),
    legend.text  = element_text(size = 14, face = "bold"),
    
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 13, face = "bold"),
    
    plot.title = element_blank()
  )

ggsave(
  filename = file.path(
    out_dir,
    paste0(
      "beta_diversity_pcoa_",
      timepoint_label,
      "_",
      sex_keep,
      ".tiff"
    )
  ),
  plot = p_beta,
  device = "tiff",
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

###############################################################################
# ALPHA DIVERSITY — RAREFIED COUNTS
###############################################################################

count_mat_sub_round <- round(count_mat_sub)

ps_round <- phyloseq(
  otu_table(count_mat_sub_round, taxa_are_rows = TRUE),
  TAX,
  META
)

set.seed(123)

ps_rare <- rarefy_even_depth(
  ps_round,
  sample.size = 24000,
  rngseed = 123,
  replace = FALSE,
  trimOTUs = FALSE,
  verbose = FALSE
)

n_before <- nsamples(ps_round)
n_after <- nsamples(ps_rare)
dropped <- setdiff(sample_names(ps_round), sample_names(ps_rare))

cat("\nSamples before rarefaction:", n_before, "\n")
cat("Samples after rarefaction:", n_after, "\n")
cat("Dropped samples:", ifelse(length(dropped) == 0, "None", paste(dropped, collapse = ", ")), "\n")

alpha_df <- estimate_richness(
  ps_rare,
  measures = c("Observed", "Shannon", "Simpson")
) %>%
  rownames_to_column("SampleID") %>%
  left_join(
    meta_ps %>% rownames_to_column("SampleID"),
    by = "SampleID"
  )

pvals_alpha <- alpha_df %>%
  summarise(
    Observed_p = wilcox.test(Observed ~ Treatment)$p.value,
    Shannon_p = wilcox.test(Shannon ~ Treatment)$p.value,
    Simpson_p = wilcox.test(Simpson ~ Treatment)$p.value
  ) %>%
  pivot_longer(
    everything(),
    names_to = "Metric",
    values_to = "p_value"
  ) %>%
  mutate(Metric = gsub("_p", "", Metric))

alpha_long <- alpha_df %>%
  pivot_longer(
    cols = c("Observed", "Shannon", "Simpson"),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  left_join(pvals_alpha, by = "Metric")

p_label_df <- alpha_long %>%
  group_by(Metric) %>%
  summarise(p_value = first(p_value), .groups = "drop")

p_alpha <- ggplot(alpha_long, aes(x = Treatment, y = Value, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2.5) +
  scale_fill_manual(values = c(
    "Control" = "steelblue",
    "Dosed" = "firebrick"
  )) +
  facet_wrap(~ Metric, scales = "free_y") +
  geom_text(
    data = p_label_df,
    aes(x = 1.5, y = Inf, label = paste0("p = ", signif(p_value, 3))),
    vjust = 1.5,
    size = 4,
    inherit.aes = FALSE
  ) +
  theme_bw(base_size = 12) +
  theme(
    
    legend.title = element_text(size = 15, face = "bold"),
    legend.text  = element_text(size = 14, face = "bold"),
    
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 13, face = "bold"),
    
    strip.text = element_text(size = 15, face = "bold"),
    
    plot.title = element_blank()
    
  ) +
  labs(
    x = NULL,
    y = "Diversity Value"
  )

ggsave(
  filename = file.path(
    out_dir,
    paste0(
      "alpha_diversity_",
      timepoint_label,
      "_",
      sex_keep,
      ".tiff"
    )
  ),
  plot = p_alpha,
  device = "tiff",
  width = 10,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

alpha_stats_summary <- alpha_long %>%
  group_by(Metric, Treatment) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    sd_value = sd(Value, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Treatment,
    values_from = c(mean_value, sd_value, n_samples)
  ) %>%
  left_join(
    pvals_alpha,
    by = "Metric"
  ) %>%
  mutate(
    Sex = sex_keep,
    Timepoint = timepoint_label,
    Tissue = tissue_label
  )
###############################################################################
# DIFFERENTIAL ORGANISM ABUNDANCE — RELATIVE ABUNDANCE
###############################################################################

compute_diff_rel <- function(df, sex, group1 = "Dosed", group2 = "Control") {
  
  df_sub <- df %>%
    filter(
      Sex == sex,
      Treatment %in% c(group1, group2)
    )
  
  res <- df_sub %>%
    group_by(tax_id, Species) %>%
    summarise(
      mean_group1 = mean(RelAbund[Treatment == group1], na.rm = TRUE),
      mean_group2 = mean(RelAbund[Treatment == group2], na.rm = TRUE),
      mean_group1_percent = mean_group1 * 100,
      mean_group2_percent = mean_group2 * 100,
      prevalence_group1 = sum(RelAbund[Treatment == group1] > 0, na.rm = TRUE),
      prevalence_group2 = sum(RelAbund[Treatment == group2] > 0, na.rm = TRUE),
      n_group1 = sum(Treatment == group1),
      n_group2 = sum(Treatment == group2),
      pval = tryCatch(
        wilcox.test(
          RelAbund[Treatment == group1],
          RelAbund[Treatment == group2],
          exact = FALSE
        )$p.value,
        error = function(e) NA_real_
      ),
      log2FC = log2((mean_group1 + 1e-6) / (mean_group2 + 1e-6)),
      .groups = "drop"
    ) %>%
    filter(!is.na(pval)) %>%
    mutate(
      padj = p.adjust(pval, method = "BH"),
      Sex = sex,
      Comparison = paste(group1, "vs", group2),
      direction = case_when(
        log2FC > 0 ~ group1,
        log2FC < 0 ~ group2,
        TRUE ~ "No change"
      ),
      significant = padj < 0.05
    ) %>%
    arrange(padj, pval)
  
  return(res)
}

df_da <- rel_long_taxid %>%
  filter(tax_id %in% taxa_keep)

res_da <- compute_diff_rel(
  df_da,
  sex = sex_keep,
  group1 = "Dosed",
  group2 = "Control"
)

sig_da <- res_da %>%
  filter(padj < 0.05)

cat("\nMinimum raw Wilcoxon p-value:\n")
print(min(res_da$pval, na.rm = TRUE))

cat("\nMinimum Wilcoxon adjusted p-value:\n")
print(min(res_da$padj, na.rm = TRUE))

cat("\nNumber significant after Wilcoxon BH FDR:\n")
print(nrow(sig_da))

###############################################################################
# PRESENCE / ABSENCE ANALYSIS FOR ALL TESTED TAXA
###############################################################################

presence_absence_df <- df_da %>%
  mutate(Present = RelAbund > 0)

presence_absence_stats <- presence_absence_df %>%
  group_by(tax_id, Species) %>%
  summarise(
    present_dosed = sum(Present[Treatment == "Dosed"], na.rm = TRUE),
    absent_dosed = sum(!Present[Treatment == "Dosed"], na.rm = TRUE),
    present_control = sum(Present[Treatment == "Control"], na.rm = TRUE),
    absent_control = sum(!Present[Treatment == "Control"], na.rm = TRUE),
    fisher_p = tryCatch(
      fisher.test(
        matrix(
          c(
            present_dosed, absent_dosed,
            present_control, absent_control
          ),
          nrow = 2,
          byrow = TRUE
        )
      )$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  filter(!is.na(fisher_p)) %>%
  mutate(
    fisher_p_adj = p.adjust(fisher_p, method = "BH"),
    prevalence_dosed_percent =
      present_dosed / (present_dosed + absent_dosed) * 100,
    prevalence_control_percent =
      present_control / (present_control + absent_control) * 100,
    prevalence_difference_percent =
      prevalence_dosed_percent - prevalence_control_percent,
    presence_direction = case_when(
      prevalence_difference_percent > 0 ~ "Dosed",
      prevalence_difference_percent < 0 ~ "Control",
      TRUE ~ "No difference"
    ),
    fisher_significant = fisher_p_adj < 0.05
  ) %>%
  arrange(fisher_p, fisher_p_adj)

presence_absence_ranked <- presence_absence_stats %>%
  arrange(
    desc(abs(prevalence_difference_percent)),
    fisher_p
  )

###############################################################################
# DA PLOT — ONLY IF SIGNIFICANT TAXA EXIST
###############################################################################

if (nrow(sig_da) > 0) {
  
  p_da <- ggplot(
    sig_da,
    aes(
      x = log2FC,
      y = reorder(Species, log2FC),
      fill = factor(direction, levels = c("Control", "Dosed"))
    )
  ) +
    geom_col() +
    scale_fill_manual(
      values = c("Control" = "steelblue", "Dosed" = "firebrick"),
      labels = c("Control" = "Depleted", "Dosed" = "Enriched")
    ) +
    theme_bw(base_size = 14) +
    labs(
      title = paste0("Significant Differential Abundance: ", timepoint_label, " ", sex_keep),
      x = "Log2 Fold Change: Dosed vs Control",
      y = "Species",
      fill = "Direction"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(
        "DA_significant_species_",
        timepoint_label,
        "_",
        sex_keep,
        ".tiff"
      )
    ),
    plot = p_da,
    device = "tiff",
    width = 10,
    height = 6,
    units = "in",
    dpi = 300,
    compression = "lzw"
  )
  
} else {
  cat("\nNo significant taxa after Wilcoxon FDR correction, so no DA plot was saved.\n")
}


###############################################################################
# BARPLOTS
###############################################################################

ps_family <- tax_glom(ps_rel, taxrank = "Family", NArm = TRUE)
ps_genus <- tax_glom(ps_rel, taxrank = "Genus", NArm = TRUE)
ps_species <- tax_glom(ps_rel, taxrank = "Species", NArm = TRUE)

plot_bar_rank <- function(ps_obj, rank, topN = 20) {
  
  df <- psmelt(ps_obj)
  
  top_taxa <- df %>%
    group_by(.data[[rank]]) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    arrange(desc(Abundance)) %>%
    slice_head(n = topN) %>%
    pull(.data[[rank]])
  
  df <- df %>%
    mutate(Taxon = ifelse(.data[[rank]] %in% top_taxa, .data[[rank]], "Other"))
  
  taxa_levels <- unique(df$Taxon)
  
  base_cols <- c(
    brewer.pal(12, "Paired"),
    brewer.pal(8, "Dark2"),
    brewer.pal(8, "Set3")
  )
  
  base_cols <- base_cols[seq_along(taxa_levels)]
  names(base_cols) <- taxa_levels
  base_cols["Other"] <- "grey70"
  
  ggplot(df, aes(x = Sample, y = Abundance, fill = Taxon)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = base_cols) +
    theme_bw(base_size = 12) +
    theme(
      
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 12,
        face = "bold"
      ),
      
      axis.text.y = element_text(
        size = 12,
        face = "bold"
      ),
      
      axis.title = element_text(
        size = 16,
        face = "bold"
      ),
      
      legend.title = element_text(
        size = 15,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 13,
        face = "bold"
      ),
      
      plot.title = element_blank(),
      
      legend.position = "right"
      
    ) +
    labs(
      x = "Sample",
      y = "Relative Abundance"
    )
}

p_family <- plot_bar_rank(ps_family, "Family", topN = 20)
p_genus <- plot_bar_rank(ps_genus, "Genus", topN = 20)
p_species <- plot_bar_rank(ps_species, "Species", topN = 20)

ggsave(
  filename = file.path(
    out_dir,
    paste0(
      "family_barplot_",
      timepoint_label,
      "_",
      sex_keep,
      ".tiff"
    )
  ),
  plot = p_family,
  device = "tiff",
  width = 10,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    out_dir,
    paste0(
      "genus_barplot_",
      timepoint_label,
      "_",
      sex_keep,
      ".tiff"
    )
  ),
  plot = p_genus,
  device = "tiff",
  width = 10,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    out_dir,
    paste0(
      "species_barplot_",
      timepoint_label,
      "_",
      sex_keep,
      ".tiff"
    )
  ),
  plot = p_species,
  device = "tiff",
  width = 10,
  height = 6,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

###############################################################################
# EXPORT TABLES
###############################################################################

permanova_export <- as.data.frame(permanova) %>%
  rownames_to_column("Term") %>%
  mutate(
    Sex = sex_keep,
    Timepoint = timepoint_label,
    Tissue = tissue_label
  )

beta_coordinates <- bray_df %>%
  mutate(
    Timepoint = timepoint_label,
    Tissue = tissue_label
  )

beta_axis_summary <- tibble(
  Sex = sex_keep,
  Timepoint = timepoint_label,
  Tissue = tissue_label,
  Axis = c("Axis.1", "Axis.2"),
  Variance_Explained_Percent = c(axis1_pct, axis2_pct)
)

alpha_export <- alpha_long %>%
  mutate(
    Sex = sex_keep,
    Timepoint = timepoint_label,
    Tissue = tissue_label
  ) %>%
  select(
    SampleID,
    Sex,
    Treatment,
    Metric,
    Value,
    p_value,
    Timepoint,
    Tissue
  )

rarefaction_summary <- tibble(
  Sex = sex_keep,
  Timepoint = timepoint_label,
  Tissue = tissue_label,
  samples_before_rarefaction = n_before,
  samples_after_rarefaction = n_after,
  dropped_samples = ifelse(length(dropped) == 0, "None", paste(dropped, collapse = "; "))
)

###############################################################################
# CLEAN COLLAPSED READS — FILTER TO CURRENT SEX FOR THIS WORKBOOK
###############################################################################

clean_collapsed_reads_export <- clean_collapsed_reads %>%
  filter(Sex == sex_keep)

###############################################################################
# SPECIES-LEVEL RELATIVE ABUNDANCE
###############################################################################

relative_abundance_species <- rel_long_taxid %>%
  mutate(
    Relative_Abundance_Percent = RelAbund * 100,
    Timepoint = timepoint_label,
    Tissue = tissue_label
  ) %>%
  select(
    tax_id,
    SampleID,
    Sex,
    Treatment,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    Relative_Abundance = RelAbund,
    Relative_Abundance_Percent,
    Timepoint,
    Tissue
  )

###############################################################################
# PHYLUM / FAMILY / GENUS RELATIVE ABUNDANCE
###############################################################################

make_relative_abundance_rank <- function(relative_abundance_species, rank_col) {
  
  relative_abundance_species %>%
    filter(!is.na(.data[[rank_col]]), .data[[rank_col]] != "") %>%
    group_by(
      SampleID,
      Sex,
      Treatment,
      Timepoint,
      Tissue,
      .data[[rank_col]]
    ) %>%
    summarise(
      Relative_Abundance = sum(Relative_Abundance, na.rm = TRUE),
      Relative_Abundance_Percent = sum(Relative_Abundance_Percent, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(Taxon = .data[[rank_col]]) %>%
    mutate(Taxonomic_Level = rank_col) %>%
    select(
      SampleID,
      Sex,
      Treatment,
      Timepoint,
      Tissue,
      Taxonomic_Level,
      Taxon,
      Relative_Abundance,
      Relative_Abundance_Percent
    )
}

relative_abundance_phylum <- make_relative_abundance_rank(relative_abundance_species, "Phylum")
relative_abundance_family <- make_relative_abundance_rank(relative_abundance_species, "Family")
relative_abundance_genus  <- make_relative_abundance_rank(relative_abundance_species, "Genus")

###############################################################################
# DA EXPORT
###############################################################################

da_all_export <- res_da %>%
  mutate(
    Timepoint = timepoint_label,
    Tissue = tissue_label
  ) %>%
  rename(
    mean_dosed = mean_group1,
    mean_control = mean_group2,
    mean_dosed_percent = mean_group1_percent,
    mean_control_percent = mean_group2_percent,
    prevalence_dosed = prevalence_group1,
    prevalence_control = prevalence_group2,
    n_dosed = n_group1,
    n_control = n_group2,
    wilcoxon_p = pval,
    wilcoxon_p_adj = padj
  )

###############################################################################
# GROUP MEAN ABUNDANCE — SPECIES LEVEL
###############################################################################

group_mean_abundance_species <- relative_abundance_species %>%
  group_by(
    tax_id,
    Timepoint,
    Tissue,
    Sex,
    Treatment,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species
  ) %>%
  summarise(
    mean_relative_abundance = mean(Relative_Abundance, na.rm = TRUE),
    sd_relative_abundance = sd(Relative_Abundance, na.rm = TRUE),
    mean_relative_abundance_percent = mean(Relative_Abundance_Percent, na.rm = TRUE),
    sd_relative_abundance_percent = sd(Relative_Abundance_Percent, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  left_join(
    da_all_export %>%
      select(
        tax_id,
        Species,
        DA_log2FC = log2FC,
        DA_p_value = wilcoxon_p,
        DA_p_adj = wilcoxon_p_adj,
        DA_direction = direction,
        DA_significant = significant
      ),
    by = c("tax_id", "Species")
  )

###############################################################################
# L. RODENTIUM ONLY REPORTS
###############################################################################

target_species_pattern <- "rodentium"
target_species_label <- "Limosilactobacillus rodentium"

l_rodentium_taxa <- rel_long_taxid %>%
  filter(str_detect(Species, regex(target_species_pattern, ignore_case = TRUE))) %>%
  distinct(tax_id, Species)

l_rodentium_sample_report <- rel_long_taxid %>%
  semi_join(l_rodentium_taxa, by = c("tax_id", "Species")) %>%
  mutate(
    Relative_Abundance_Percent = RelAbund * 100,
    Present = RelAbund > 0,
    Timepoint = timepoint_label,
    Tissue = tissue_label
  ) %>%
  select(
    tax_id,
    Species,
    SampleID,
    Sex,
    Treatment,
    Relative_Abundance = RelAbund,
    Relative_Abundance_Percent,
    Present,
    Timepoint,
    Tissue
  )

l_rodentium_wilcoxon_report <- res_da %>%
  filter(str_detect(Species, regex(target_species_pattern, ignore_case = TRUE))) %>%
  mutate(
    Timepoint = timepoint_label,
    Tissue = tissue_label
  )

l_rodentium_presence_report <- presence_absence_stats %>%
  filter(str_detect(Species, regex(target_species_pattern, ignore_case = TRUE))) %>%
  mutate(
    Sex = sex_keep,
    Timepoint = timepoint_label,
    Tissue = tissue_label
  )

l_rodentium_combined_report <- l_rodentium_sample_report %>%
  distinct(tax_id, Species, Sex, Timepoint, Tissue) %>%
  left_join(
    l_rodentium_wilcoxon_report %>%
      select(
        tax_id,
        Species,
        Comparison,
        mean_dosed = mean_group1,
        mean_control = mean_group2,
        mean_dosed_percent = mean_group1_percent,
        mean_control_percent = mean_group2_percent,
        prevalence_dosed_abundance = prevalence_group1,
        prevalence_control_abundance = prevalence_group2,
        log2FC,
        wilcoxon_p = pval,
        wilcoxon_p_adj = padj,
        wilcoxon_significant = significant
      ),
    by = c("tax_id", "Species")
  ) %>%
  left_join(
    l_rodentium_presence_report %>%
      select(
        tax_id,
        Species,
        present_dosed,
        absent_dosed,
        present_control,
        absent_control,
        prevalence_dosed_percent,
        prevalence_control_percent,
        prevalence_difference_percent,
        fisher_p,
        fisher_p_adj,
        fisher_significant
      ),
    by = c("tax_id", "Species")
  )

###############################################################################
# L. REUTERI STATS REPORT
###############################################################################

l_reuteri_taxa <- rel_long_taxid %>%
  filter(
    str_detect(
      Species,
      regex(
        "reuteri",
        ignore_case = TRUE
      )
    )
  ) %>%
  distinct(
    tax_id,
    Species
  )

l_reuteri_sample_report <- rel_long_taxid %>%
  semi_join(
    l_reuteri_taxa,
    by = c(
      "tax_id",
      "Species"
    )
  ) %>%
  mutate(
    Relative_Abundance_Percent = RelAbund * 100,
    Present = RelAbund > 0,
    Timepoint = timepoint_label,
    Tissue = tissue_label
  )

l_reuteri_wilcoxon_report <- res_da %>%
  filter(
    str_detect(
      Species,
      regex(
        "reuteri",
        ignore_case = TRUE
      )
    )
  ) %>%
  mutate(
    Timepoint = timepoint_label,
    Tissue = tissue_label
  )

l_reuteri_presence_report <- presence_absence_stats %>%
  filter(
    str_detect(
      Species,
      regex(
        "reuteri",
        ignore_case = TRUE
      )
    )
  ) %>%
  mutate(
    Sex = sex_keep,
    Timepoint = timepoint_label,
    Tissue = tissue_label
  )

l_reuteri_combined_report <- l_reuteri_sample_report %>%
  distinct(
    tax_id,
    Species,
    Sex,
    Timepoint,
    Tissue
  ) %>%
  left_join(
    l_reuteri_wilcoxon_report %>%
      select(
        tax_id,
        Species,
        Comparison,
        mean_dosed = mean_group1,
        mean_control = mean_group2,
        mean_dosed_percent = mean_group1_percent,
        mean_control_percent = mean_group2_percent,
        prevalence_dosed_abundance = prevalence_group1,
        prevalence_control_abundance = prevalence_group2,
        log2FC,
        wilcoxon_p = pval,
        wilcoxon_p_adj = padj,
        wilcoxon_significant = significant
      ),
    by = c(
      "tax_id",
      "Species"
    )
  ) %>%
  left_join(
    l_reuteri_presence_report %>%
      select(
        tax_id,
        Species,
        present_dosed,
        absent_dosed,
        present_control,
        absent_control,
        prevalence_dosed_percent,
        prevalence_control_percent,
        prevalence_difference_percent,
        fisher_p,
        fisher_p_adj,
        fisher_significant
      ),
    by = c(
      "tax_id",
      "Species"
    )
  )
###############################################################################
# L. RODENTIUM PLOT EXTRACTION TABLES
###############################################################################

l_rodentium_plot_sample_values <- l_rodentium_sample_report %>%
  arrange(Timepoint, Tissue, Sex, Treatment, SampleID)

l_rodentium_plot_summary <- l_rodentium_sample_report %>%
  group_by(
    tax_id,
    Species,
    Timepoint,
    Tissue,
    Sex,
    Treatment
  ) %>%
  summarise(
    n_samples = n(),
    n_present = sum(Present, na.rm = TRUE),
    prevalence_percent = n_present / n_samples * 100,
    mean_relative_abundance = mean(Relative_Abundance, na.rm = TRUE),
    sd_relative_abundance = sd(Relative_Abundance, na.rm = TRUE),
    mean_relative_abundance_percent = mean(Relative_Abundance_Percent, na.rm = TRUE),
    sd_relative_abundance_percent = sd(Relative_Abundance_Percent, na.rm = TRUE),
    .groups = "drop"
  )

l_rodentium_plot_long <- l_rodentium_plot_summary %>%
  select(
    Species,
    Timepoint,
    Tissue,
    Sex,
    Treatment,
    mean_relative_abundance_percent,
    sd_relative_abundance_percent,
    prevalence_percent
  ) %>%
  pivot_longer(
    cols = c(mean_relative_abundance_percent, prevalence_percent),
    names_to = "Plot_Metric",
    values_to = "Plot_Value"
  )
###############################################################################
# FINAL SUPPLEMENTARY WORKBOOK
###############################################################################

supplementary_list <- list(
  "read_metrics" = read_counts,
  "clean_collapsed_reads" = clean_collapsed_reads_export,
  "relative_abundance_species" = relative_abundance_species,
  "relative_abundance_genus" = relative_abundance_genus,
  "relative_abundance_phylum" = relative_abundance_phylum,
  "relative_abundance_family" = relative_abundance_family,
  "group_mean_abundance_species" = group_mean_abundance_species,
  "alpha_sample_values" = alpha_export,
  "alpha_stats_summary" = alpha_stats_summary,
  "rarefaction_summary" = rarefaction_summary,
  "PERMANOVA" = permanova_export,
  "beta_coordinates" = beta_coordinates,
  "beta_axis_summary" = beta_axis_summary,
  "DA_all_relative_abundance" = da_all_export,
  "presence_absence_stats_all_taxa" = presence_absence_stats,
  "L_reuteri_stats_report" = l_reuteri_combined_report,
  "L_rodentium_sample_report" = l_rodentium_sample_report,
  "L_rodentium_stats_report" = l_rodentium_combined_report,
  "L_rodentium_plot_sample_values" = l_rodentium_plot_sample_values,
  "L_rodentium_plot_summary" = l_rodentium_plot_summary,
  "L_rodentium_plot_long" = l_rodentium_plot_long
)

write_xlsx(
  supplementary_list,
  path = file.path(
    out_dir,
    paste0("Supplementary_16S_outputs_", timepoint_label, "_", sex_keep, ".xlsx")
  )
)

###############################################################################
# SAVE EXTRA CSV FILES
###############################################################################

write.csv(
  da_all_export,
  file.path(out_dir, paste0("DA_all_relative_abundance_", timepoint_label, "_", sex_keep, ".csv")),
  row.names = FALSE
)

write.csv(
  presence_absence_stats,
  file.path(out_dir, paste0("presence_absence_stats_all_taxa_", timepoint_label, "_", sex_keep, ".csv")),
  row.names = FALSE
)

write.csv(
  l_rodentium_combined_report,
  file.path(out_dir, paste0("L_rodentium_stats_report_", timepoint_label, "_", sex_keep, ".csv")),
  row.names = FALSE
)

write.csv(
  l_reuteri_combined_report,
  file.path(
    out_dir,
    paste0(
      "L_reuteri_stats_report_",
      timepoint_label,
      "_",
      sex_keep,
      ".csv"
    )
  ),
  row.names = FALSE
)
write_xlsx(
  list(
    "sample_values" = l_rodentium_plot_sample_values,
    "plot_summary" = l_rodentium_plot_summary,
    "plot_long" = l_rodentium_plot_long,
    "stats_report" = l_rodentium_combined_report
  ),
  path = file.path(
    out_dir,
    paste0("L_rodentium_plot_extraction_", timepoint_label, "_", sex_keep, ".xlsx")
  )
)
###############################################################################
# DONE
###############################################################################

cat("\nDONE. Outputs saved to:\n")
cat(out_dir, "\n")

cat("\nFinal workbook tabs exported:\n")
print(names(supplementary_list))