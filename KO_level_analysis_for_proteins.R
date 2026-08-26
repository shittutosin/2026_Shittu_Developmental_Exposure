# ============================================================
# KO-LEVEL STACKED BAR PLOTS FROM expr_imp FILES
# STACK BY TAXA
# SAMPLE-LEVEL NORMALIZATION TO TOTAL MEASURED PROTEOME
#
# BARS:
#   For each biological sample, summed linear KO abundance is divided by
#   the summed linear abundance of ALL measured proteins in that same sample.
#   Control/Dosed bars are then the mean sample-level relative KO abundance.
#   Taxa segments retain their proportional contribution to the KO bar.
#
# STATS:
#   PATHWAY-STYLE KO-LEVEL WILCOXON
#   Uses protein-level log2_intensity values
#   Control vs Dosed within each Timepoint and Sex
#   BH-adjusted within Pathway x Sex
#
# INPUTS USED
#   - expr_imp_* files (sample-level log2 intensities)
#   - KO_map_file.xlsx
#   - Pathway_intermediate_tables/Protein_KO_membership_all_sources.csv
#   - six limma summary .xlsx files (for organism extraction via Function)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(readxl)
  library(writexl)
  library(ggplot2)
  library(forcats)
  library(purrr)
  library(scales)
})

# ============================================================
# 0. SETTINGS
# ============================================================

setwd("/data/proteomics/KO_level_analysis")

out_dir <- "KO_selected_abundance_stacked_by_taxa_relative_proteome"
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)

# ---- Relative abundance stacked taxa plot folders ----
dir.create(
  file.path(out_dir, "plots", "relative_proteome_stacked_taxa"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(out_dir, "plots", "relative_proteome_stacked_taxa", "tiff"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(out_dir, "plots", "relative_proteome_stacked_taxa", "png"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(out_dir, "plots", "relative_proteome_stacked_taxa", "pdf"),
  recursive = TRUE,
  showWarnings = FALSE
)
plot_width  <- 9
plot_height <- 6
plot_dpi    <- 600

top_n_taxa <- 12

# ============================================================
# 1. FILES
# ============================================================

expr_files <- c(
  "expr_imp_5wks_M_8_2_3_0.05.xlsx",
  "expr_imp_8wks_M_8_2_3_0.05.xlsx",
  "expr_imp_12wks_M_8_2_3_0.05.xlsx",
  "expr_imp_5wks_F_8_2_3_0.05.xlsx",
  "expr_imp_8wks_F_8_2_3_0.05.xlsx",
  "expr_imp_12wks_F_8_2_3_0.05.xlsx"
)

limma_files <- c(
  "5wks_M_8_2_3_0.05.xlsx",
  "8wks_M_8_2_3_0.05.xlsx",
  "12wks_M_8_2_3_0.05.xlsx",
  "5wks_F_8_2_3_0.05.xlsx",
  "8wks_F_8_2_3_0.05.xlsx",
  "12wks_F_8_2_3_0.05.xlsx"
)

ko_map_file <- "KO_map_file.xlsx"
protein_ko_map_file <- file.path("Pathway_intermediate_tables", "Protein_KO_membership_all_sources.csv")

# ============================================================
# 2. HELPERS
# ============================================================

clean_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

extract_timepoint_from_expr <- function(x) {
  case_when(
    str_detect(x, "5wks")  ~ "5W",
    str_detect(x, "8wks")  ~ "8W",
    str_detect(x, "12wks") ~ "12W",
    TRUE ~ NA_character_
  )
}

extract_sex_from_expr <- function(x) {
  case_when(
    str_detect(x, "_M_") ~ "Males",
    str_detect(x, "_F_") ~ "Females",
    TRUE ~ NA_character_
  )
}

extract_timepoint_from_limma <- function(x) {
  case_when(
    str_detect(x, "^5wks_")  ~ "5W",
    str_detect(x, "^8wks_")  ~ "8W",
    str_detect(x, "^12wks_") ~ "12W",
    TRUE ~ NA_character_
  )
}

extract_sex_from_limma <- function(x) {
  case_when(
    str_detect(x, "_M_") ~ "Males",
    str_detect(x, "_F_") ~ "Females",
    TRUE ~ NA_character_
  )
}

extract_organism <- function(x) {
  out <- rep(NA_character_, length(x))
  
  idx_bracket <- str_detect(x, "\\[[^\\]]+\\]")
  out[idx_bracket] <- str_match(x[idx_bracket], "\\[([^\\]]+)\\]")[, 2]
  
  idx_tax <- is.na(out) & str_detect(x, "Tax=")
  out[idx_tax] <- str_match(
    x[idx_tax],
    "Tax=([^\\[]+?)(?:\\s+TaxID=|\\s+RepID=|$)"
  )[, 2]
  
  out <- str_squish(out)
  out[out == ""] <- NA_character_
  out
}

p_to_stars <- function(p) {
  case_when(
    is.na(p)   ~ "ns",
    p < 1e-4   ~ "****",
    p < 1e-3   ~ "***",
    p < 1e-2   ~ "**",
    p < 5e-2   ~ "*",
    TRUE       ~ "ns"
  )
}

make_taxa_palette <- function(taxa_levels) {
  base_cols <- c(
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
    "#66a61e", "#e6ab02", "#a6761d", "#6a3d9a",
    "#1f78b4", "#b2df8a", "#fb9a99", "#fdbf6f",
    "#cab2d6", "#ffff99"
  )
  
  special_taxa <- c("Other", "Unknown organism")
  regular_taxa <- setdiff(taxa_levels, special_taxa)
  
  n_regular <- length(regular_taxa)
  
  if (n_regular <= length(base_cols)) {
    cols_regular <- base_cols[seq_len(n_regular)]
  } else {
    cols_regular <- colorRampPalette(base_cols)(n_regular)
  }
  
  pal <- setNames(cols_regular, regular_taxa)
  
  if ("Other" %in% taxa_levels) {
    pal["Other"] <- "grey70"
  }
  if ("Unknown organism" %in% taxa_levels) {
    pal["Unknown organism"] <- "grey40"
  }
  
  pal[taxa_levels]
}

compute_ko_stats <- function(df_long) {
  stats_tbl <- df_long %>%
    group_by(KO, name, Pathway, Sex, Timepoint) %>%
    group_modify(~{
      dat <- .x %>%
        filter(!is.na(log2_intensity), !is.na(Group)) %>%
        filter(Group %in% c("Control", "Dosed"))
      
      n_control <- sum(dat$Group == "Control")
      n_dosed   <- sum(dat$Group == "Dosed")
      
      if (n_control < 1 || n_dosed < 1) {
        return(tibble(
          n_control = n_control,
          n_dosed = n_dosed,
          p_value = NA_real_
        ))
      }
      
      wt <- tryCatch(
        wilcox.test(log2_intensity ~ Group, data = dat, exact = FALSE),
        error = function(e) NULL
      )
      
      tibble(
        n_control = n_control,
        n_dosed = n_dosed,
        p_value = if (is.null(wt)) NA_real_ else wt$p.value
      )
    }) %>%
    ungroup()
  
  stats_tbl %>%
    group_by(Pathway, Sex) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      sig_label = p_to_stars(p_adj)
    ) %>%
    ungroup()
}

# ============================================================
# 3. READ KO MAP FILE
# ============================================================

ko_map <- read_xlsx(ko_map_file, sheet = "Sheet1") %>%
  rename_with(~ str_trim(.x)) %>%
  mutate(across(where(is.character), str_squish))

required_ko_cols <- c("KO", "name", "Pathway")
missing_ko_cols <- setdiff(required_ko_cols, colnames(ko_map))
if (length(missing_ko_cols) > 0) {
  stop("KO_map_file.xlsx is missing columns: ", paste(missing_ko_cols, collapse = ", "))
}

ko_map <- ko_map %>%
  filter(!is.na(KO), KO != "") %>%
  distinct(KO, .keep_all = TRUE)

target_kos <- unique(ko_map$KO)

# ============================================================
# 4. READ PROTEIN -> KO MAPPING
# ============================================================

protein_ko_map <- read_csv(protein_ko_map_file, show_col_types = FALSE) %>%
  mutate(
    Protein = as.character(Protein),
    KO = as.character(KO)
  ) %>%
  filter(!is.na(Protein), Protein != "", !is.na(KO), KO != "") %>%
  filter(KO %in% target_kos) %>%
  distinct(Protein, KO)

if (nrow(protein_ko_map) == 0) {
  stop("No proteins from Protein_KO_membership_all_sources.csv matched the KOs in KO_map_file.xlsx")
}

# ============================================================
# 5. READ LIMMA FILES TO EXTRACT ORGANISM NAMES
# ============================================================

read_one_limma_file <- function(f) {
  read_xlsx(f) %>%
    mutate(
      source_file = basename(f),
      Timepoint = extract_timepoint_from_limma(basename(f)),
      Sex = extract_sex_from_limma(basename(f))
    ) %>%
    select(any_of(c("Protein", "Function", "Timepoint", "Sex"))) %>%
    mutate(
      Protein = as.character(Protein),
      Function = as.character(Function),
      Organism = extract_organism(Function),
      Organism = str_replace(Organism, "\\s+DSM\\s+\\d+$", ""),
      Organism = ifelse(is.na(Organism) | Organism == "", "Unknown organism", Organism)
    ) %>%
    distinct(Protein, Timepoint, Sex, Organism, .keep_all = TRUE)
}

organism_lookup <- map_dfr(limma_files, read_one_limma_file) %>%
  distinct(Protein, Timepoint, Sex, Organism)

organism_lookup_global <- organism_lookup %>%
  distinct(Protein, Organism) %>%
  group_by(Protein) %>%
  summarise(
    Organism = first(Organism[!is.na(Organism)]),
    .groups = "drop"
  ) %>%
  mutate(
    Organism = ifelse(is.na(Organism) | Organism == "", "Unknown organism", Organism)
  )

# ============================================================
# 6. READ expr_imp FILES
# ============================================================

read_one_expr_file <- function(f) {
  df <- read_xlsx(f)
  sample_cols <- grep("^(Control|Dosed)_", colnames(df), value = TRUE)
  
  if (!("Protein" %in% colnames(df))) {
    stop("File is missing Protein column: ", f)
  }
  if (length(sample_cols) == 0) {
    stop("No sample columns like Control_* or Dosed_* found in: ", f)
  }
  
  df %>%
    mutate(
      Protein = as.character(Protein),
      Timepoint = extract_timepoint_from_expr(basename(f)),
      Sex = extract_sex_from_expr(basename(f)),
      source_file = basename(f)
    ) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "Sample",
      values_to = "log2_intensity"
    ) %>%
    mutate(
      log2_intensity = as.numeric(log2_intensity),
      Group = case_when(
        str_detect(Sample, "^Control_") ~ "Control",
        str_detect(Sample, "^Dosed_")   ~ "Dosed",
        TRUE ~ NA_character_
      ),
      Sample_ID = paste(Timepoint, Sex, Sample, sep = "_")
    )
}

expr_long <- map_dfr(expr_files, read_one_expr_file)

# ============================================================
# 7. CALCULATE TOTAL MEASURED PROTEOME PER SAMPLE
#    BEFORE RESTRICTING TO SELECTED KOs
# ============================================================

# IMPORTANT:
# The denominator is calculated from expr_long, which contains ALL proteins
# in each expr_imp file. Therefore, total_proteome_linear is NOT restricted
# to proteins present in KO_map_file.xlsx or protein_ko_map_file.
#
# For each biological sample:
#   total_proteome_linear = sum(2^log2_intensity) across ALL measured proteins

sample_total_proteome <- expr_long %>%
  mutate(
    linear_abundance = 2^log2_intensity
  ) %>%
  group_by(
    Sex, Timepoint, Sample, Sample_ID, Group
  ) %>%
  summarise(
    total_proteome_linear = sum(linear_abundance, na.rm = TRUE),
    n_proteins_total = n_distinct(Protein[!is.na(log2_intensity)]),
    .groups = "drop"
  )

# Safety check: each sample should have a finite, positive proteome denominator.
bad_denominators <- sample_total_proteome %>%
  filter(!is.finite(total_proteome_linear) | total_proteome_linear <= 0)

if (nrow(bad_denominators) > 0) {
  warning(
    nrow(bad_denominators),
    " sample(s) have non-finite or non-positive total proteome abundance. " ,
    "Their relative KO abundance will be set to NA."
  )
}


# ============================================================
# 8. JOIN EXPRESSION TO KO MAP + ORGANISM
# ============================================================

# Only now do we restrict the numerator to proteins assigned to selected KOs.
# The total-proteome denominator was already calculated using ALL proteins.

expr_ko <- expr_long %>%
  inner_join(protein_ko_map, by = "Protein") %>%
  left_join(ko_map, by = "KO") %>%
  left_join(organism_lookup, by = c("Protein", "Timepoint", "Sex")) %>%
  left_join(
    organism_lookup_global %>% rename(Organism_global = Organism),
    by = "Protein"
  ) %>%
  mutate(
    Organism = coalesce(Organism, Organism_global, "Unknown organism"),
    Organism = ifelse(
      is.na(Organism) | Organism == "",
      "Unknown organism",
      Organism
    ),
    linear_abundance = 2^log2_intensity
  ) %>%
  left_join(
    sample_total_proteome,
    by = c("Sex", "Timepoint", "Sample", "Sample_ID", "Group")
  ) %>%
  select(
    Protein, KO, name, Pathway, Timepoint, Sex, Sample, Sample_ID, Group,
    Organism, log2_intensity, linear_abundance,
    total_proteome_linear, n_proteins_total
  )

if (nrow(expr_ko) == 0) {
  stop("After joining expression data to protein-KO mapping, no rows remained.")
}


# ============================================================
# 9. SUM WITHIN EACH SAMPLE: KO x SAMPLE x ORGANISM
# ============================================================

ko_org_sample <- expr_ko %>%
  group_by(
    KO, name, Pathway, Sex, Timepoint, Sample, Sample_ID, Group,
    Organism, total_proteome_linear
  ) %>%
  summarise(
    ko_org_linear = sum(linear_abundance, na.rm = TRUE),
    n_proteins = n_distinct(Protein),
    .groups = "drop"
  ) %>%
  mutate(
    relative_ko_abundance_taxon = ifelse(
      is.finite(total_proteome_linear) & total_proteome_linear > 0,
      ko_org_linear / total_proteome_linear,
      NA_real_
    )
  )

# Total abundance of each KO within each biological sample.
ko_sample_total <- ko_org_sample %>%
  group_by(
    KO, name, Pathway, Sex, Timepoint, Sample, Sample_ID, Group,
    total_proteome_linear
  ) %>%
  summarise(
    ko_linear_total = sum(ko_org_linear, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    relative_ko_abundance = ifelse(
      is.finite(total_proteome_linear) & total_proteome_linear > 0,
      ko_linear_total / total_proteome_linear,
      NA_real_
    )
  )


# ============================================================
# 10. TAXA COLLAPSE + GROUP SUMMARY FOR STACKED BARS
# ============================================================

# Keep the top taxa per KO x Sex based on summed linear KO abundance.
top_taxa_by_ko_sex <- ko_org_sample %>%
  group_by(KO, Sex, Organism) %>%
  summarise(
    total_abundance = sum(ko_org_linear, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(KO, Sex) %>%
  arrange(desc(total_abundance), .by_group = TRUE) %>%
  mutate(
    rank_taxa = row_number(),
    Organism_collapsed = ifelse(
      rank_taxa <= top_n_taxa,
      Organism,
      "Other"
    )
  ) %>%
  ungroup() %>%
  select(KO, Sex, Organism, Organism_collapsed)

ko_org_sample_collapsed <- ko_org_sample %>%
  left_join(
    top_taxa_by_ko_sex,
    by = c("KO", "Sex", "Organism")
  ) %>%
  mutate(
    Organism_collapsed = coalesce(Organism_collapsed, "Other")
  ) %>%
  group_by(
    KO, name, Pathway, Sex, Timepoint, Sample, Sample_ID, Group,
    total_proteome_linear,
    Organism = Organism_collapsed
  ) %>%
  summarise(
    ko_org_linear = sum(ko_org_linear, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    relative_ko_abundance_taxon = ifelse(
      is.finite(total_proteome_linear) & total_proteome_linear > 0,
      ko_org_linear / total_proteome_linear,
      NA_real_
    )
  )

# ------------------------------------------------------------
# Group-level taxonomic contributions
# ------------------------------------------------------------
# Each taxon segment is first normalized within each biological sample,
# then averaged across biological samples in that Timepoint x Group.

ko_org_group_summary <- ko_org_sample_collapsed %>%
  group_by(KO, name, Pathway, Sex, Timepoint, Group, Organism) %>%
  summarise(
    mean_linear = mean(ko_org_linear, na.rm = TRUE),
    median_linear = median(ko_org_linear, na.rm = TRUE),
    mean_relative_taxon = mean(relative_ko_abundance_taxon, na.rm = TRUE),
    median_relative_taxon = median(relative_ko_abundance_taxon, na.rm = TRUE),
    n_samples = n_distinct(Sample_ID),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Group-level total KO relative abundance
# ------------------------------------------------------------
# This is explicitly the mean of the SAMPLE-LEVEL normalized KO values.

ko_group_total <- ko_sample_total %>%
  group_by(KO, name, Pathway, Sex, Timepoint, Group) %>%
  summarise(
    relative_ko_abundance = mean(relative_ko_abundance, na.rm = TRUE),
    median_relative_ko_abundance = median(relative_ko_abundance, na.rm = TRUE),
    mean_ko_linear = mean(ko_linear_total, na.rm = TRUE),
    n_samples_total = n_distinct(Sample_ID),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Join bar height back to taxonomic contributions
# ------------------------------------------------------------
# Taxa fractions are calculated from the mean sample-level relative
# contributions, and plot_height is scaled so all taxon segments sum
# exactly to the mean sample-level relative KO abundance for that bar.

ko_org_group_summary <- ko_org_group_summary %>%
  left_join(
    ko_group_total,
    by = c("KO", "name", "Pathway", "Sex", "Timepoint", "Group")
  ) %>%
  group_by(KO, name, Pathway, Sex, Timepoint, Group) %>%
  mutate(
    total_mean_relative_taxon = sum(mean_relative_taxon, na.rm = TRUE),
    taxa_fraction = ifelse(
      total_mean_relative_taxon > 0,
      mean_relative_taxon / total_mean_relative_taxon,
      0
    ),
    plot_height = taxa_fraction * relative_ko_abundance
  ) %>%
  ungroup()

# ------------------------------------------------------------
# INTERNAL SCALING CHECK
# ------------------------------------------------------------
# For every KO x Sex x Timepoint x Group:
#   sum(plot_height across taxa)
# must equal
#   mean sample-level relative KO abundance
# apart from floating-point rounding.

ko_relative_scaling_check <- ko_org_group_summary %>%
  group_by(KO, name, Pathway, Sex, Timepoint, Group) %>%
  summarise(
    stacked_bar_height = sum(plot_height, na.rm = TRUE),
    expected_relative_ko_abundance = first(relative_ko_abundance),
    difference = stacked_bar_height - expected_relative_ko_abundance,
    .groups = "drop"
  )

# ============================================================
# 11. PANEL-LEVEL STATS
# ============================================================

panel_stats <- compute_ko_stats(expr_ko)

# ============================================================
# ============================================================
# 12. SUPPLEMENTARY INFORMATION TABLES
# ============================================================

ko_protein_table <- expr_ko %>%
  distinct(KO, name, Pathway, Protein, Organism) %>%
  arrange(Pathway, KO, Organism, Protein)

ko_taxa_contributor_rank <- ko_org_group_summary %>%
  group_by(KO, name, Pathway, Sex, Timepoint, Group) %>%
  arrange(desc(taxa_fraction), .by_group = TRUE) %>%
  mutate(
    contributor_rank = row_number(),
    percent_contribution = taxa_fraction * 100
  ) %>%
  ungroup() %>%
  select(
    KO, name, Pathway, Sex, Timepoint, Group,
    contributor_rank,
    Organism,
    mean_linear,
    median_linear,
    mean_relative_taxon,
    median_relative_taxon,
    relative_ko_abundance,
    median_relative_ko_abundance,
    mean_ko_linear,
    taxa_fraction,
    percent_contribution,
    plot_height,
    n_samples
  )

ko_coverage <- expr_ko %>%
  group_by(KO, name, Pathway, Sex, Timepoint) %>%
  summarise(
    n_proteins = n_distinct(Protein),
    n_organisms = n_distinct(Organism),
    n_samples = n_distinct(Sample_ID),
    .groups = "drop"
  ) %>%
  arrange(Pathway, KO, Sex, Timepoint)

analysis_parameters <- tibble(
  Parameter = c(
    "Input expression files",
    "Input KO map",
    "Protein-KO mapping file",
    "Abundance scale",
    "Linear abundance calculation",
    "Proteome denominator",
    "Stacked bar y-axis",
    "Normalization scope",
    "Cross-sex interpretation",
    "Taxa collapse rule",
    "Statistical test",
    "Comparison",
    "Multiple testing correction",
    "FDR grouping",
    "Significance labels"
  ),
  Value = c(
    paste(expr_files, collapse = "; "),
    ko_map_file,
    protein_ko_map_file,
    "Protein-level log2 intensity from expr_imp files",
    "linear_abundance = 2^log2_intensity",
    "All measured proteins in each expr_imp sample, calculated before KO filtering",
    "Relative KO abundance",
    "Each biological sample is normalized independently: summed linear abundance of proteins assigned to the KO divided by summed linear abundance of ALL measured proteins in that same sample; group means are calculated after sample-level normalization",
    "Relative KO abundance is directly comparable across experimental conditions because every sample is normalized to its own total measured proteome; interpretation across sexes still depends on comparable proteomic acquisition and preprocessing",
    paste0("Top ", top_n_taxa, " taxa per KO x Sex retained; remaining taxa collapsed as Other"),
    "Wilcoxon rank-sum test",
    "Control vs Dosed within each KO x Sex x Timepoint",
    "Benjamini-Hochberg",
    "Adjusted within Pathway x Sex",
    "ns, *, **, ***, **** based on BH-adjusted p-value"
  )
)

write_csv(
  ko_protein_table,
  file.path(out_dir, "tables", "selected_KO_contributing_proteins.csv")
)

write_csv(
  sample_total_proteome,
  file.path(out_dir, "tables", "sample_total_measured_proteome.csv")
)

write_csv(
  ko_org_sample_collapsed,
  file.path(out_dir, "tables", "selected_KO_taxa_linear_abundance_per_sample.csv")
)

write_csv(
  ko_org_group_summary,
  file.path(out_dir, "tables", "selected_KO_taxa_group_mean_relative_abundance.csv")
)

write_csv(
  ko_relative_scaling_check,
  file.path(out_dir, "tables", "selected_KO_relative_scaling_check.csv")
)

write_csv(
  ko_sample_total,
  file.path(out_dir, "tables", "selected_KO_total_linear_abundance_per_sample.csv")
)

write_csv(
  panel_stats,
  file.path(out_dir, "tables", "selected_KO_panel_stats.csv")
)

write_csv(
  ko_taxa_contributor_rank,
  file.path(out_dir, "tables", "selected_KO_taxa_contributor_rank_percent.csv")
)

write_csv(
  ko_coverage,
  file.path(out_dir, "tables", "selected_KO_coverage_summary.csv")
)

write_csv(
  analysis_parameters,
  file.path(out_dir, "tables", "selected_KO_analysis_parameters.csv")
)

write_xlsx(
  list(
    README_analysis_parameters = analysis_parameters,
    KO_map = ko_map,
    Protein_KO_mapping_used = protein_ko_map,
    KO_contributing_proteins = ko_protein_table,
    Sample_total_proteome = sample_total_proteome,
    KO_total_linear_per_sample = ko_sample_total,
    KO_taxa_linear_per_sample = ko_org_sample_collapsed,
    KO_taxa_group_mean_relative = ko_org_group_summary,
    KO_relative_scaling_check = ko_relative_scaling_check,
    KO_taxa_contributor_rank = ko_taxa_contributor_rank,
    KO_panel_stats = panel_stats,
    KO_coverage_summary = ko_coverage
  ),
  file.path(out_dir, "tables", "selected_KO_Supplementary_Information.xlsx")
)
# ============================================================
# 13. GLOBAL TAXA PALETTE
# ============================================================

all_taxa_levels <- ko_org_group_summary %>%
  group_by(Organism) %>%
  summarise(total_abundance = sum(mean_relative_taxon, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  pull(Organism)

taxa_palette <- make_taxa_palette(all_taxa_levels)

write_csv(
  tibble(Organism = names(taxa_palette), Color = unname(taxa_palette)),
  file.path(out_dir, "tables", "taxa_color_key.csv")
)

# ============================================================
# ============================================================
# ============================================================
# 14. PLOT FUNCTION
# Generates a plot even when the KO is completely absent
# from the selected sex
# ============================================================

plot_one_ko_sex <- function(target_ko, target_sex,
                            df_plot,
                            df_stats,
                            out_dir_tiff,
                            out_dir_png,
                            out_dir_pdf) {
  
  timepoint_levels <- c("5W", "8W", "12W")
  group_levels <- c("Control", "Dosed")
  
  # ----------------------------------------------------------
  # Obtain KO metadata from the original KO map.
  # This allows the title to be generated even when the KO
  # is completely absent from the selected sex.
  # ----------------------------------------------------------
  
  ko_metadata <- ko_map %>%
    filter(KO == target_ko) %>%
    slice(1)
  
  if (nrow(ko_metadata) == 0) {
    warning("No KO metadata found for: ", target_ko)
    return(NULL)
  }
  
  ko_name <- ko_metadata$name[1]
  ko_pathway <- ko_metadata$Pathway[1]
  
  # ----------------------------------------------------------
  # Select abundance data for this KO and sex
  # ----------------------------------------------------------
  
  one <- df_plot %>%
    filter(KO == target_ko, Sex == target_sex) %>%
    mutate(
      Timepoint = factor(Timepoint, levels = timepoint_levels),
      Group = factor(Group, levels = group_levels)
    )
  
  # ----------------------------------------------------------
  # Select statistical results
  # ----------------------------------------------------------
  
  one_stats <- df_stats %>%
    filter(KO == target_ko, Sex == target_sex) %>%
    mutate(
      Timepoint = factor(Timepoint, levels = timepoint_levels)
    )
  
  plot_title <- paste0(ko_name, " (", target_ko, ")")
  
  microbiome_label <- case_when(
    target_sex == "Males"   ~ "Male Colon Microbiome",
    target_sex == "Females" ~ "Female Colon Microbiome",
    TRUE                    ~ target_sex
  )
  
  plot_subtitle <- microbiome_label
  
  # ----------------------------------------------------------
  # Determine which timepoints contain detectable abundance
  # ----------------------------------------------------------
  
  detected_timepoints <- one %>%
    filter(
      !is.na(relative_ko_abundance),
      relative_ko_abundance > 0
    ) %>%
    distinct(Timepoint) %>%
    pull(Timepoint) %>%
    as.character()
  
  missing_timepoints <- setdiff(
    timepoint_levels,
    detected_timepoints
  )
  
  # Labels for missing panels
  missing_labels <- tibble(
    Timepoint = factor(
      missing_timepoints,
      levels = timepoint_levels
    ),
    Group = factor(
      rep("Control", length(missing_timepoints)),
      levels = group_levels
    ),
    x = 1.5,
    y = 0,
    label = "Not detected"
  )
  
  # ----------------------------------------------------------
  # Determine maximum relative abundance in each panel
  # ----------------------------------------------------------
  
  if (nrow(one) > 0) {
    
    panel_max <- one %>%
      distinct(Timepoint, Group, relative_ko_abundance) %>%
      group_by(Timepoint) %>%
      summarise(
        ymax = ifelse(
          all(is.na(relative_ko_abundance)),
          0,
          max(relative_ko_abundance, na.rm = TRUE)
        ),
        .groups = "drop"
      ) %>%
      complete(
        Timepoint = factor(
          timepoint_levels,
          levels = timepoint_levels
        ),
        fill = list(ymax = 0)
      )
    
  } else {
    
    # Completely absent KO: create three empty panel records
    panel_max <- tibble(
      Timepoint = factor(
        timepoint_levels,
        levels = timepoint_levels
      ),
      ymax = 0
    )
  }
  
  # ----------------------------------------------------------
  # Add statistical annotation only to detected panels
  # ----------------------------------------------------------
  
  one_stats <- one_stats %>%
    left_join(panel_max, by = "Timepoint") %>%
    filter(
      !is.na(ymax),
      ymax > 0,
      !is.na(p_value)
    ) %>%
    mutate(
      bracket_y = ymax * 1.10,
      bracket_y_low = ymax * 1.05,
      label_y = ymax * 1.23,
      x_start = 1,
      x_end = 2,
      label_x = 1.5,
      stat_text = sig_label
    )
  
  # ----------------------------------------------------------

  # Set a suitable KO-specific y-axis
  # ----------------------------------------------------------
  
  abundance_values <- one$relative_ko_abundance
  abundance_values <- abundance_values[
    is.finite(abundance_values)
  ]
  
  stat_values <- one_stats$label_y
  stat_values <- stat_values[
    is.finite(stat_values)
  ]
  
  y_top <- max(
    c(abundance_values, stat_values),
    na.rm = TRUE
  )
  
  # Add ~10% headroom above the highest bar/annotation
  if (!is.finite(y_top) || y_top <= 0) {
    y_top <- 0.10
  } else {
    y_top <- y_top * 1.10
  }

  # Place "Not detected" labels inside the actual y-axis range.
  missing_labels <- missing_labels %>%
    mutate(y = y_top * 0.45)
  
  # ----------------------------------------------------------
  # Construct plot
  # ----------------------------------------------------------
  
  p <- ggplot(
    one,
    aes(
      x = Group,
      y = plot_height,
      fill = Organism
    )
  ) +
    
    geom_col(
      width = 0.78,
      color = "black",
      linewidth = 0.2
    ) +
    
    geom_segment(
      data = one_stats,
      aes(
        x = x_start,
        xend = x_end,
        y = bracket_y,
        yend = bracket_y
      ),
      inherit.aes = FALSE,
      linewidth = 0.5
    ) +
    
    geom_segment(
      data = one_stats,
      aes(
        x = x_start,
        xend = x_start,
        y = bracket_y_low,
        yend = bracket_y
      ),
      inherit.aes = FALSE,
      linewidth = 0.5
    ) +
    
    geom_segment(
      data = one_stats,
      aes(
        x = x_end,
        xend = x_end,
        y = bracket_y_low,
        yend = bracket_y
      ),
      inherit.aes = FALSE,
      linewidth = 0.5
    ) +
    
    geom_text(
      data = one_stats,
      aes(
        x = label_x,
        y = label_y,
        label = stat_text
      ),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    ) +
    
    geom_text(
      data = missing_labels,
      aes(
        x = x,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      size = 4,
      fontface = "bold.italic",
      color = "grey35"
    ) +
    
    facet_wrap(
      ~Timepoint,
      nrow = 1,
      scales = "fixed",
      drop = FALSE
    ) +
    
    scale_x_discrete(
      limits = group_levels,
      drop = FALSE
    ) +
    
    scale_fill_manual(
      values = taxa_palette,
      drop = FALSE
    ) +
    
    scale_y_continuous(
      limits = c(0, y_top),
      breaks = scales::pretty_breaks(n = 5),
      labels = scales::label_percent(
        accuracy = 0.001
      ),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = NULL,
      y = "Relative KO abundance (%)",
      fill = "Taxa"
    ) +
    
    theme_bw(base_size = 12) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 17
      ),
      plot.subtitle = element_text(
        face = "bold",
        size = 12
      ),
      strip.text = element_text(
        face = "bold",
        size = 12
      ),
      strip.background = element_rect(
        fill = "grey85",
        color = "grey50",
        linewidth = 0.6
      ), 
      axis.text.x = element_text(
        face = "bold",
        size = 10
      ),
      axis.text.y = element_text(
        size = 9
      ),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(
        color = "grey88",
        linewidth = 0.35
      ),
      legend.position = if (nrow(one) > 0) "right" else "none",
      legend.title = element_text(
        face = "bold",
        size = 14
      ),
      legend.text = element_text(
        face = "bold",
        size = 13
      ),
      plot.margin = margin(12, 12, 12, 12)
    )
  
  # ----------------------------------------------------------
  # Save plot
  # ----------------------------------------------------------
  
  file_stub <- paste(
    clean_filename(ko_name),
    clean_filename(target_ko),
    clean_filename(target_sex),
    "stacked_taxa_relative_proteome",
    sep = "_"
  )
  
  ggsave(
    file.path(
      out_dir_tiff,
      paste0(file_stub, ".tiff")
    ),
    p,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi,
    compression = "lzw"
  )
  
  ggsave(
    file.path(
      out_dir_png,
      paste0(file_stub, ".png")
    ),
    p,
    width = plot_width,
    height = plot_height,
    dpi = 300
  )
  
  ggsave(
    file.path(
      out_dir_pdf,
      paste0(file_stub, ".pdf")
    ),
    p,
    width = plot_width,
    height = plot_height
  )
  
  return(p)
}
# 15. GENERATE ALL PLOTS
# ============================================================

plot_keys <- expand.grid(
  KO = unique(ko_map$KO),
  Sex = c("Males", "Females"),
  stringsAsFactors = FALSE
) %>%
  as_tibble() %>%
  left_join(ko_map %>% select(KO, name, Pathway), by = "KO") %>%
  distinct(KO, name, Pathway, Sex)

plot_list <- vector("list", nrow(plot_keys))
names(plot_list) <- paste(plot_keys$KO, plot_keys$Sex, sep = " | ")

for (i in seq_len(nrow(plot_keys))) {
  this_ko <- plot_keys$KO[i]
  this_sex <- plot_keys$Sex[i]
  
  plot_list[[i]] <- plot_one_ko_sex(
    target_ko = this_ko,
    target_sex = this_sex,
    df_plot = ko_org_group_summary,
    df_stats = panel_stats,
    out_dir_tiff = file.path(out_dir, "plots", "relative_proteome_stacked_taxa", "tiff"),
    out_dir_png  = file.path(out_dir, "plots", "relative_proteome_stacked_taxa", "png"),
    out_dir_pdf  = file.path(out_dir, "plots", "relative_proteome_stacked_taxa", "pdf")
  )
}

write_csv(
  expr_ko,
  file.path(
    out_dir,
    "tables",
    "selected_KO_protein_level_values_used_for_statistics.csv"
  )
)
message("Done.")

ko_org_group_summary %>%
  filter(
    KO == "K01186",
    Sex == "Females"
  ) %>%
  distinct(
    Timepoint,
    Group,
    relative_ko_abundance,
    median_relative_ko_abundance,
    mean_ko_linear,
    n_samples_total
  ) %>%
  arrange(Timepoint, Group)


