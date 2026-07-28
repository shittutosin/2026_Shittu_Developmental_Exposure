# ============================================================
# KO-LEVEL STACKED BAR PLOTS FROM expr_imp FILES
# STACK BY TAXA
#
# BARS:
#   Log10 mean summed KO protein abundance by group
#   (stacked by taxa contribution)
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

setwd("data/proteomics/KO_analysis_for_proteins")

out_dir <- "KO_level_analysis"
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)

# ---- Relative abundance stacked taxa plot folders ----
dir.create(
  file.path(out_dir, "plots", "relative_abundance_stacked_taxa"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(out_dir, "plots", "relative_abundance_stacked_taxa", "tiff"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(out_dir, "plots", "relative_abundance_stacked_taxa", "png"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(out_dir, "plots", "relative_abundance_stacked_taxa", "pdf"),
  recursive = TRUE,
  showWarnings = FALSE
)
plot_width  <- 13
plot_height <- 7.5
plot_dpi    <- 600

top_n_taxa <- 20

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
      
      if (n_control < 2 || n_dosed < 2) {
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
# 7. JOIN EXPRESSION TO KO MAP + ORGANISM
# ============================================================

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
    Organism = ifelse(is.na(Organism) | Organism == "", "Unknown organism", Organism),
    linear_abundance = 2^log2_intensity
  ) %>%
  select(
    Protein, KO, name, Pathway, Timepoint, Sex, Sample, Sample_ID, Group,
    Organism, log2_intensity, linear_abundance
  )

if (nrow(expr_ko) == 0) {
  stop("After joining expression data to protein-KO mapping, no rows remained.")
}

# ============================================================
# 8. SUM WITHIN SAMPLE: KO x SAMPLE x ORGANISM
# ============================================================

ko_org_sample <- expr_ko %>%
  group_by(KO, name, Pathway, Sex, Timepoint, Sample, Sample_ID, Group, Organism) %>%
  summarise(
    ko_org_linear = sum(linear_abundance, na.rm = TRUE),
    n_proteins = n_distinct(Protein),
    .groups = "drop"
  )

ko_sample_total <- ko_org_sample %>%
  group_by(KO, name, Pathway, Sex, Timepoint, Sample, Sample_ID, Group) %>%
  summarise(
    ko_linear_total = sum(ko_org_linear, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# 9. TAXA COLLAPSE: KEEP TOP 20 PER KO x SEX OVERALL
# ============================================================

top_taxa_by_ko_sex <- ko_org_sample %>%
  group_by(KO, Sex, Organism) %>%
  summarise(
    total_abundance = sum(ko_org_linear, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(KO, Sex) %>%
  arrange(desc(total_abundance), .by_group = TRUE) %>%
  mutate(rank_taxa = row_number()) %>%
  mutate(Organism_collapsed = ifelse(rank_taxa <= top_n_taxa, Organism, "Other")) %>%
  ungroup() %>%
  select(KO, Sex, Organism, Organism_collapsed)

ko_org_sample_collapsed <- ko_org_sample %>%
  left_join(top_taxa_by_ko_sex, by = c("KO", "Sex", "Organism")) %>%
  mutate(
    Organism_collapsed = coalesce(Organism_collapsed, "Other")
  ) %>%
  group_by(
    KO, name, Pathway, Sex, Timepoint, Sample, Sample_ID, Group,
    Organism = Organism_collapsed
  ) %>%
  summarise(
    ko_org_linear = sum(ko_org_linear, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# ============================================================
# 10. SUMMARY FOR STACKED BARS
# Bar total height = log10(total KO abundance)
# Taxa colors show proportional contribution
# ============================================================

ko_org_group_summary <- ko_org_sample_collapsed %>%
  group_by(KO, name, Pathway, Sex, Timepoint, Group, Organism) %>%
  summarise(
    mean_linear = mean(ko_org_linear, na.rm = TRUE),
    median_linear = median(ko_org_linear, na.rm = TRUE),
    n_samples = n_distinct(Sample_ID),
    .groups = "drop"
  ) %>%
  group_by(KO, name, Pathway, Sex, Timepoint, Group) %>%
  mutate(
    total_mean_linear = sum(mean_linear, na.rm = TRUE),
    taxa_fraction = ifelse(total_mean_linear > 0, mean_linear / total_mean_linear, 0),
    total_log10_abundance = log10(total_mean_linear + 1),
    plot_height = taxa_fraction * total_log10_abundance
  ) %>%
  ungroup()
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
    total_mean_linear,
    taxa_fraction,
    percent_contribution,
    total_log10_abundance,
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
    "Stacked bar y-axis",
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
    "log10(total mean KO linear abundance + 1)",
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
  ko_org_sample_collapsed,
  file.path(out_dir, "tables", "selected_KO_taxa_linear_abundance_per_sample.csv")
)

write_csv(
  ko_org_group_summary,
  file.path(out_dir, "tables", "selected_KO_taxa_group_mean_log10_abundance.csv")
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
    KO_total_linear_per_sample = ko_sample_total,
    KO_taxa_linear_per_sample = ko_org_sample_collapsed,
    KO_taxa_group_mean_log10 = ko_org_group_summary,
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
  summarise(total_abundance = sum(mean_linear, na.rm = TRUE), .groups = "drop") %>%
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
      !is.na(total_log10_abundance),
      total_log10_abundance > 0
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
    y = 0.5,
    label = "Not detected"
  )
  
  # ----------------------------------------------------------
  # Determine maximum abundance in each panel
  # ----------------------------------------------------------
  
  if (nrow(one) > 0) {
    
    panel_max <- one %>%
      distinct(Timepoint, Group, total_log10_abundance) %>%
      group_by(Timepoint) %>%
      summarise(
        ymax = ifelse(
          all(is.na(total_log10_abundance)),
          0,
          max(total_log10_abundance, na.rm = TRUE)
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
      bracket_y = ymax + 0.35,
      bracket_y_low = ymax + 0.20,
      label_y = ymax + 0.48,
      x_start = 1,
      x_end = 2,
      label_x = 1.5,
      stat_text = sig_label
    )
  
  # ----------------------------------------------------------
  # Set a suitable common y-axis
  # ----------------------------------------------------------
  
  abundance_values <- one$total_log10_abundance
  abundance_values <- abundance_values[
    is.finite(abundance_values)
  ]
  
  stat_values <- one_stats$label_y
  stat_values <- stat_values[
    is.finite(stat_values)
  ]
  
  y_top <- max(
    c(abundance_values, stat_values, 1),
    na.rm = TRUE
  ) + 0.35
  
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
      size = 5.1,
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
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = NULL,
      y = "Log10 mean KO protein abundance",
      fill = "Taxa"
    ) +
    
    theme_bw(base_size = 12) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 15
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
        size = 11
      ),
      axis.text.y = element_text(
        size = 10
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
        size = 9
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
    "stacked_taxa_log10_abundance",
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
    out_dir_tiff = file.path(out_dir, "plots", "relative_abundance_stacked_taxa", "tiff"),
    out_dir_png  = file.path(out_dir, "plots", "relative_abundance_stacked_taxa", "png"),
    out_dir_pdf  = file.path(out_dir, "plots", "relative_abundance_stacked_taxa", "pdf")
  )
}

message("Done.")

