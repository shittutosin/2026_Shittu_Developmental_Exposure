###############################################################################
# 16S BETA DIVERSITY: CONTROL AND DOSED MALE VS FEMALE
#
#
#
# Each output is a three-panel PCoA figure: 5W | 8W | 12W
# Figures are saved as 600-dpi TIFF files with LZW compression.
# PERMANOVA uses Bray-Curtis distances and 9,999 reproducible permutations.
###############################################################################

suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(tidyverse)
  library(RColorBrewer)
  library(writexl)
})

# =============================================================================
# 1. SET PATHS
# =============================================================================

base_dir <- "data/16s_analysis/rel-abundance_files"
setwd(base_dir)

out_dir <- file.path(base_dir, "beta_diversity_outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 2. FIND INPUT FILES AND READ METADATA
# =============================================================================

files <- list.files(
  base_dir,
  pattern = "_rel-abundance_(5W|8W|12W)\\.tsv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No matching relative-abundance TSV files were found in base_dir.")
}

samples <- gsub("\\.tsv$", "", basename(files))

meta <- read.table(
  "metadata.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_meta_cols <- c("SampleID", "Sex", "Treatment", "Timepoint")
missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))

if (length(missing_meta_cols) > 0) {
  stop(
    "Metadata is missing required columns: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

meta <- meta %>%
  mutate(
    SampleID  = trimws(as.character(SampleID)),
    Sex       = factor(trimws(as.character(Sex)), levels = c("M", "F")),
    Treatment = factor(
      trimws(as.character(Treatment)),
      levels = c("Control", "Dosed")
    ),
    Timepoint = factor(
      trimws(as.character(Timepoint)),
      levels = c("5W", "8W", "12W")
    )
  )

# =============================================================================
# 3. READ ABUNDANCE FILES
# =============================================================================

long_df <- map2_dfr(files, samples, function(file, samp) {
  
  dat <- read_tsv(
    file,
    col_types = cols(),
    show_col_types = FALSE
  )
  
  required_cols <- c(
    "tax_id",
    "lineage",
    "abundance",
    "estimated counts"
  )
  
  missing_cols <- setdiff(required_cols, colnames(dat))
  
  if (length(missing_cols) > 0) {
    stop(
      "File is missing required columns: ",
      basename(file),
      "\nMissing: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  dat %>%
    dplyr::select(
      tax_id,
      lineage,
      abundance,
      `estimated counts`
    ) %>%
    dplyr::rename(
      est_counts = `estimated counts`
    ) %>%
    mutate(sample = samp)
})

long_df <- long_df %>%
  filter(!tax_id %in% c("mapped_unclassified", "unmapped")) %>%
  mutate(
    abundance = replace_na(as.numeric(abundance), 0),
    est_counts = replace_na(as.numeric(est_counts), 0)
  ) %>%
  group_by(tax_id, sample) %>%
  summarise(
    abundance  = abundance[1],
    est_counts = est_counts[1],
    lineage    = lineage[1],
    .groups = "drop"
  )

# =============================================================================
# 4. CREATE COUNT MATRIX
# =============================================================================

otu_counts <- long_df %>%
  dplyr::select(tax_id, sample, est_counts) %>%
  pivot_wider(
    names_from = sample,
    values_from = est_counts,
    values_fill = 0
  ) %>%
  arrange(tax_id)

count_mat <- as.matrix(otu_counts[, -1, drop = FALSE])
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- otu_counts$tax_id

cat("\nSamples in files but not metadata:\n")
print(setdiff(colnames(count_mat), meta$SampleID))

cat("\nSamples in metadata but not files:\n")
print(setdiff(meta$SampleID, colnames(count_mat)))

keep_samples <- intersect(meta$SampleID, colnames(count_mat))

if (length(keep_samples) == 0) {
  stop("No sample names match between the metadata and abundance files.")
}

count_mat_sub <- count_mat[, keep_samples, drop = FALSE]

meta_ps <- meta %>%
  filter(SampleID %in% keep_samples) %>%
  arrange(match(SampleID, keep_samples)) %>%
  column_to_rownames("SampleID")

# =============================================================================
# 5. CREATE TAXONOMY TABLE
# =============================================================================

ranks <- c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

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
      
      if (length(vals) < 7) {
        vals <- c(rep(NA_character_, 7 - length(vals)), vals)
      }
      
      rev(vals)
    })
  ) %>%
  ungroup() %>%
  unnest_wider(lineage_split, names_sep = "_") %>%
  rename_with(~ ranks, starts_with("lineage_split_")) %>%
  mutate(across(all_of(ranks), trimws)) %>%
  mutate(
    Species = trimws(Species),
    Species = na_if(Species, ""),
    Species = ifelse(
      is.na(Species) |
        Species == "" |
        Species == Genus |
        Species == word(Genus, 1) |
        Species %in% c(
          "sp",
          "sp.",
          "bacterium",
          "uncultured bacterium"
        ),
      paste(Genus, "sp."),
      Species
    )
  ) %>%
  group_by(Species) %>%
  mutate(
    Species = ifelse(
      n() == 1,
      Species,
      paste0(Species, "_", row_number())
    )
  ) %>%
  ungroup() %>%
  column_to_rownames("tax_id")

tax_mat <- tax_mat[, ranks, drop = FALSE]

# =============================================================================
# 6. CREATE PHYLOSEQ OBJECT
# =============================================================================

OTU <- otu_table(count_mat_sub, taxa_are_rows = TRUE)

TAX <- tax_table(
  as.matrix(
    tax_mat[
      rownames(count_mat_sub),
      ,
      drop = FALSE
    ]
  )
)

META <- sample_data(meta_ps)

ps <- phyloseq(OTU, TAX, META)

ps_rel <- transform_sample_counts(
  ps,
  function(x) {
    if (sum(x) == 0) {
      x
    } else {
      x / sum(x)
    }
  }
)

ps_rel <- prune_samples(sample_sums(ps_rel) > 0, ps_rel)
ps_rel <- prune_taxa(taxa_sums(ps_rel) > 0, ps_rel)

# =============================================================================
# 7. SIGNIFICANCE LABELS AND COLORS
# =============================================================================

sig_label <- function(p) {
  case_when(
    is.na(p)     ~ "",
    p <= 0.0001  ~ "****",
    p <= 0.001   ~ "***",
    p <= 0.01    ~ "**",
    p <= 0.05    ~ "*",
    TRUE         ~ "ns"
  )
}

timepoint_cols <- brewer.pal(3, "Dark2")
names(timepoint_cols) <- c("5W", "8W", "12W")


# =============================================================================
# 8. REPRODUCIBLE PERMANOVA: MALE VS FEMALE WITHIN EACH TIMEPOINT
# =============================================================================


sex_within_timepoint_permanova <- function(ps_obj) {
  
  meta_df <- data.frame(sample_data(ps_obj))
  bray_dist <- phyloseq::distance(ps_obj, method = "bray")
  bray_mat <- as.matrix(bray_dist)
  
  results <- list()
  summary_results <- list()
  
  timepoints <- c("5W", "8W", "12W")
  
  for (i in seq_along(timepoints)) {
    
    tp <- timepoints[i]
    
    sub_meta <- meta_df %>%
      filter(
        Timepoint == tp,
        Sex %in% c("M", "F")
      )
    
    sub_meta$Sex <- droplevels(
      factor(
        sub_meta$Sex,
        levels = c("M", "F")
      )
    )
    
    if (nlevels(sub_meta$Sex) < 2) {
      warning("Both sexes are not present at ", tp, ".")
      next
    }
    
    sub_dist <- as.dist(
      bray_mat[
        rownames(sub_meta),
        rownames(sub_meta),
        drop = FALSE
      ]
    )
    
    # A fixed, timepoint-specific seed makes each p-value reproducible.
    set.seed(3000 + i)
    
    perm <- adonis2(
      sub_dist ~ Sex,
      data = sub_meta,
      permutations = 9999
    )
    
    comparison_name <- paste0(tp, ": M vs F")
    
    summary_results[[tp]] <- data.frame(
      Comparison = comparison_name,
      n_group1 = sum(sub_meta$Sex == "M"),
      n_group2 = sum(sub_meta$Sex == "F"),
      F_statistic = perm$F[1],
      R2 = perm$R2[1],
      p_value = perm$`Pr(>F)`[1],
      permutations = 9999,
      Timepoint = tp
    )
    
    perm_table <- as.data.frame(perm) %>%
      rownames_to_column("Term") %>%
      mutate(
        Term = ifelse(Term == "Sex", "Model", Term),
        Sex = "M vs F",
        Timepoint = tp,
        Tissue = "Colon",
        Comparison = comparison_name
      ) %>%
      dplyr::select(
        Term,
        Df,
        SumOfSqs,
        R2,
        F,
        `Pr(>F)`,
        Sex,
        Timepoint,
        Tissue,
        Comparison
      )
    
    results[[tp]] <- perm_table
  }
  
  list(
    summary = bind_rows(summary_results) %>%
      mutate(
        Timepoint = factor(
          Timepoint,
          levels = c("5W", "8W", "12W")
        ),
        significance = sig_label(p_value),
        label = paste0(
          "p = ",
          signif(p_value, 3),
          " (",
          significance,
          ")"
        )
      ),
    full = bind_rows(results)
  )
}

# =============================================================================
# 9. THREE-PANEL PLOTTING FUNCTION
# =============================================================================

make_faceted_sex_beta_plot <- function(
    ps_obj,
    title_text,
    file_name
) {
  
  meta_df <- data.frame(sample_data(ps_obj)) %>%
    rownames_to_column("SampleID")
  
  bray_dist <- phyloseq::distance(
    ps_obj,
    method = "bray"
  )
  
  ord_bray <- ordinate(
    ps_obj,
    method = "PCoA",
    distance = bray_dist
  )
  
  bray_df <- plot_ordination(
    ps_obj,
    ord_bray,
    type = "samples",
    justDF = TRUE
  ) %>%
    rownames_to_column("SampleID") %>%
    dplyr::select(
      SampleID,
      Axis.1,
      Axis.2
    ) %>%
    left_join(
      meta_df,
      by = "SampleID"
    ) %>%
    mutate(
      Timepoint = factor(
        Timepoint,
        levels = c("5W", "8W", "12W")
      ),
      Sex = factor(
        Sex,
        levels = c("M", "F")
      ),
      Timepoint_Sex = interaction(
        Timepoint,
        Sex,
        drop = TRUE,
        sep = "_"
      )
    )
  
  eig <- ord_bray$values$Eigenvalues
  var_explained <- eig / sum(eig)
  
  axis1_pct <- round(var_explained[1] * 100, 1)
  axis2_pct <- round(var_explained[2] * 100, 1)
  
  perm_results <- sex_within_timepoint_permanova(ps_obj)
  
  # One p-value row per facet. x = Inf and y = Inf place the label
  # in the upper-right corner instead of over the sample cloud.
  pval_df <- perm_results$summary %>%
    dplyr::select(
      Timepoint,
      label
    )
  
  p <- ggplot(
    bray_df,
    aes(
      x = Axis.1,
      y = Axis.2,
      color = Timepoint
    )
  ) +
    stat_ellipse(
      aes(
        group = Timepoint_Sex,
        linetype = Sex
      ),
      type = "t",
      linewidth = 1,
      level = 0.85,
      na.rm = TRUE
    ) +
    geom_point(
      aes(shape = Sex),
      size = 4,
      alpha = 0.9
    ) +
    facet_wrap(
      ~ Timepoint,
      nrow = 1,
      scales = "fixed"
    ) +
    geom_label(
      data = pval_df,
      aes(
        x = Inf,
        y = Inf,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = 1.05,
      vjust = 1.15,
      size = 4.5,
      fontface = "bold",
      color = "black",
      fill = "white",
      label.size = 0.25,
      label.padding = unit(0.22, "lines")
    ) +
    scale_color_manual(
      values = timepoint_cols,
      drop = FALSE
    ) +
    scale_shape_manual(
      values = c(
        "M" = 16,
        "F" = 17
      ),
      labels = c(
        "M" = "Male",
        "F" = "Female"
      ),
      drop = FALSE
    ) +
    scale_linetype_manual(
      values = c(
        "M" = "solid",
        "F" = "dashed"
      ),
      labels = c(
        "M" = "Male",
        "F" = "Female"
      ),
      drop = FALSE
    ) +
    labs(
      title = NULL,
      x = paste0("Axis 1 (", axis1_pct, "%)"),
      y = paste0("Axis 2 (", axis2_pct, "%)"),
      color = "Timepoint",
      shape = "Sex",
      linetype = "Sex"
    ) +
    theme_bw(base_size = 14) +
    theme(
      strip.background = element_rect(
        fill = "white",
        color = "black",
        linewidth = 0.7
      ),
      
      strip.text = element_text(
        face = "bold",
        size = 17
      ),
      
      axis.title.x = element_text(
        face = "bold",
        size = 17
      ),
      
      axis.title.y = element_text(
        face = "bold",
        size = 17
      ),
      
      axis.text.x = element_text(
        face = "bold",
        size = 12
      ),
      
      axis.text.y = element_text(
        face = "bold",
        size = 12
      ),
      
      legend.position = "bottom",
      
      legend.title = element_text(
        face = "bold",
        size = 14
      ),
      
      legend.text = element_text(
        face = "bold",
        size = 12
      ),
      
      panel.spacing = unit(1, "lines"),
      
      plot.margin = margin(10, 15, 10, 10)
    )
  
  ggsave(
    filename = file.path(out_dir, file_name),
    plot = p,
    device = "tiff",
    width = 16,
    height = 7,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
  )
  
  list(
    plot = p,
    permanova_summary = perm_results$summary,
    permanova_full = perm_results$full
  )
}


# =============================================================================
# 11. CREATE CONTROL AND DOSED SUBSETS
# =============================================================================

ps_controls_all <- subset_samples(
  ps_rel,
  Treatment == "Control" & Sex %in% c("M", "F")
)

ps_dosed_all <- subset_samples(
  ps_rel,
  Treatment == "Dosed" & Sex %in% c("M", "F")
)

# Remove taxa that are empty after subsetting.
ps_controls_all <- prune_taxa(
  taxa_sums(ps_controls_all) > 0,
  ps_controls_all
)

ps_dosed_all <- prune_taxa(
  taxa_sums(ps_dosed_all) > 0,
  ps_dosed_all
)



# =============================================================================
# 12. CREATE THE TWO MANUSCRIPT FIGURES
# =============================================================================

# Control: male versus female within 5W, 8W, and 12W
res_control <- make_faceted_sex_beta_plot(
  ps_obj = ps_controls_all,
  title_text = "Colon Microbiome: Controls | Male vs Female",
  file_name = "beta_controls_male_vs_female_three_panels.tiff"
)

# Dosed: male versus female within 5W, 8W, and 12W
res_dosed <- make_faceted_sex_beta_plot(
  ps_obj = ps_dosed_all,
  title_text = "Colon Microbiome: Dosed | Male vs Female",
  file_name = "beta_dosed_male_vs_female_three_panels.tiff"
)

# Display both figures in RStudio.
print(res_control$plot)
print(res_dosed$plot)



# =============================================================================
# 13. EXPORT PERMANOVA RESULTS
# =============================================================================

PERMANOVA_summary <- bind_rows(
  res_control$permanova_summary %>% mutate(Treatment = "Control"),
  res_dosed$permanova_summary %>% mutate(Treatment = "Dosed")
) %>%
  dplyr::select(
    Treatment,
    Timepoint,
    Comparison,
    n_group1,
    n_group2,
    F_statistic,
    R2,
    p_value,
    significance,
    permutations
  )

PERMANOVA_full <- bind_rows(
  res_control$permanova_full %>% mutate(Treatment = "Control"),
  res_dosed$permanova_full %>% mutate(Treatment = "Dosed")
) %>%
  dplyr::select(
    Treatment,
    Term,
    Df,
    SumOfSqs,
    R2,
    F,
    `Pr(>F)`,
    Sex,
    Timepoint,
    Tissue,
    Comparison
  )

write_xlsx(
  list(
    "PERMANOVA_summary" = PERMANOVA_summary,
    "PERMANOVA_full_adonis2" = PERMANOVA_full
  ),
  path = file.path(
    out_dir,
    "Supplementary_16S_beta_diversity_control_dosed.xlsx"
  )
)

cat("\nDONE.\n")
cat("Two three-panel TIFF figures and the PERMANOVA workbook were saved to:\n")
cat(out_dir, "\n")
