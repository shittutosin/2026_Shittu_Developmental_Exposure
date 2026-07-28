# ============================================================
# HADEG FUNCTIONAL-NODE PROTEIN-LEVEL BOXPLOTS
# STATS MATCH KEGG VERSION:
#   - Uses protein-level log2_intensity
#   - Wilcoxon Control vs Dosed per node x sex x timepoint
#   - BH correction by HADEG_node x Sex
#   - Y-axis = log2 protein intensity
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readxl)
  library(readr)
  library(writexl)
  library(ggplot2)
  library(purrr)
  library(forcats)
})

# ============================================================
# 0. SETTINGS
# ============================================================

setwd("data/proteomics/functional_summary_HADEG_proteins")

out_dir <- "HADEG_functional_node_boxplots"

dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots", "tiff"), recursive = TRUE, showWarnings = FALSE)

timepoint_levels <- c("5W", "8W", "12W")
group_levels <- c("Control", "Dosed")
sex_levels <- c("Males", "Females")

# Journal-ready figure dimensions (inches)
plot_width  <- 12
plot_height <- 8
plot_dpi    <- 600

hadeg_file <- "02_HADEG_best_hit_per_protein_for_manual_review_STRICT.csv"

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

# ============================================================
# 1. TARGET HADEG FUNCTIONS
# ============================================================

target_hadeg_functions <- tribble(
  ~HADEG_gene, ~HADEG_function, ~HADEG_class,
  "pcaF", "Beta-ketoadipyl-CoA thiolase", "Aromatic carbon funneling into central metabolism",
  "pcaI", "3-oxoadipate CoA-transferase subunit A", "Aromatic carbon funneling into central metabolism",
  "ahpC", "Alkyl hydroperoxide reductase C", "Oxidative stress response",
  "est", "PBAT ALS54749", "Esterase"
) %>%
  mutate(
    HADEG_gene_clean = str_to_lower(str_squish(HADEG_gene)),
    HADEG_function_clean = str_to_lower(str_squish(HADEG_function)),
    HADEG_node = paste0(HADEG_gene, " | ", HADEG_function)
  )

# ============================================================
# 2. HELPERS
# ============================================================

clean_filename <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

safe_plot_text <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "α", "alpha")
  x <- str_replace_all(x, "β", "beta")
  x <- str_replace_all(x, "γ", "gamma")
  x <- str_replace_all(x, "δ", "delta")
  x <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT")
  x
}

p_to_label <- function(p) {
  case_when(
    is.na(p)       ~ "NA",
    p <= 0.0001    ~ "****",
    p <= 0.001     ~ "***",
    p <= 0.01      ~ "**",
    p <= 0.05      ~ "*",
    TRUE           ~ "ns"
  )
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

pick_col <- function(df, candidates, required = TRUE, label = "column") {
  found <- candidates[candidates %in% colnames(df)][1]
  if (is.na(found) && required) {
    stop(
      "Could not find ", label, ". Tried: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns are:\n",
      paste(colnames(df), collapse = ", ")
    )
  }
  found
}

# ============================================================
# 3. READ HADEG BEST-HIT FILE
# ============================================================

hadeg_raw <- read_csv(hadeg_file, show_col_types = FALSE) %>%
  rename_with(~ str_squish(.x)) %>%
  mutate(across(where(is.character), str_squish))

protein_col <- pick_col(
  hadeg_raw,
  c("Protein", "protein", "Query", "query", "query_id", "Query_ID", "Protein_ID", "protein_id"),
  label = "protein ID column"
)

gene_col <- "HADEG_gene"
function_col <- "HADEG_description"
review_col <- "Include_in_heatmap"

hadeg_map <- hadeg_raw %>%
  mutate(
    Protein = as.character(.data[[protein_col]]),
    HADEG_gene_raw = as.character(.data[[gene_col]]),
    HADEG_function_raw = as.character(.data[[function_col]]),
    Review_status = as.character(.data[[review_col]]),
    
    HADEG_gene_clean = str_to_lower(
      str_squish(HADEG_gene_raw)
    ),
    
    HADEG_function_clean = str_to_lower(
      str_squish(HADEG_function_raw)
    )
  ) %>%
  filter(
    !is.na(Protein),
    Protein != ""
  ) %>%
  filter(
    # Retain proteins already approved for inclusion
    str_to_lower(
      coalesce(Review_status, "")
    ) %in% c(
      "yes",
      "y",
      "include",
      "included",
      "keep",
      "true"
    ) |
      
      # Also retain the exact est | PBAT ALS54749 class
      (
        HADEG_gene_clean == "est" &
          HADEG_function_clean == "pbat als54749"
      )
  )

target_hadeg_functions_join <- target_hadeg_functions %>%
  select(
    target_HADEG_gene = HADEG_gene,
    target_HADEG_function = HADEG_function,
    target_HADEG_class = HADEG_class,
    HADEG_gene_clean,
    HADEG_function_clean,
    target_HADEG_node = HADEG_node
  )

hadeg_map_target_gene <- hadeg_map %>%
  left_join(
    target_hadeg_functions_join,
    by = "HADEG_gene_clean"
  ) %>%
  filter(!is.na(target_HADEG_node)) %>%
  mutate(match_by = "gene")

hadeg_map_target_function <- hadeg_map %>%
  anti_join(
    hadeg_map_target_gene %>% distinct(Protein),
    by = "Protein"
  ) %>%
  left_join(
    target_hadeg_functions_join,
    by = "HADEG_function_clean"
  ) %>%
  filter(!is.na(target_HADEG_node)) %>%
  mutate(match_by = "function")

hadeg_map_target <- bind_rows(
  hadeg_map_target_gene,
  hadeg_map_target_function
) %>%
  transmute(
    Protein,
    HADEG_gene = target_HADEG_gene,
    HADEG_function = target_HADEG_function,
    HADEG_class = target_HADEG_class,
    HADEG_node = target_HADEG_node,
    HADEG_gene_raw,
    HADEG_function_raw,
    Review_status,
    match_by
  ) %>%
  distinct(
    Protein,
    HADEG_gene,
    HADEG_function,
    HADEG_class,
    HADEG_node,
    .keep_all = TRUE
  )

if (nrow(hadeg_map_target) == 0) {
  stop("No proteins matched your selected HADEG functions. Check gene/function column names and labels.")
}

write_csv(
  hadeg_map_target,
  file.path(out_dir, "tables", "HADEG_selected_function_protein_mapping_STRICT.csv")
)

# ============================================================
# 4. READ LIMMA FILES FOR ORGANISM LOOKUP
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
# 5. READ expr_imp FILES
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
# 6. JOIN EXPRESSION TO HADEG FUNCTIONS
# ============================================================

expr_hadeg <- expr_long %>%
  inner_join(
    hadeg_map_target %>%
      select(Protein, HADEG_gene, HADEG_function, HADEG_class, HADEG_node),
    by = "Protein"
  ) %>%
  left_join(organism_lookup, by = c("Protein", "Timepoint", "Sex")) %>%
  left_join(
    organism_lookup_global %>% rename(Organism_global = Organism),
    by = "Protein"
  ) %>%
  mutate(
    Organism = coalesce(Organism, Organism_global, "Unknown organism"),
    Organism = ifelse(is.na(Organism) | Organism == "", "Unknown organism", Organism),
    Timepoint = factor(Timepoint, levels = timepoint_levels),
    Sex = factor(Sex, levels = sex_levels),
    Group = factor(Group, levels = group_levels)
  ) %>%
  filter(
    !is.na(log2_intensity),
    !is.na(Group),
    !is.na(Timepoint),
    !is.na(Sex)
  )

if (nrow(expr_hadeg) == 0) {
  stop("After joining expression data to HADEG mapping, no rows remained.")
}

write_csv(
  expr_hadeg,
  file.path(out_dir, "tables", "HADEG_protein_level_expression_long_with_metadata_STRICT.csv")
)

# ============================================================
# 7. SCAFFOLD TO FORCE ALL NODE x SEX x TIMEPOINT PANELS
# ============================================================

hadeg_scaffold <- hadeg_map_target %>%
  distinct(HADEG_gene, HADEG_function, HADEG_class, HADEG_node) %>%
  tidyr::crossing(
    Sex = factor(sex_levels, levels = sex_levels),
    Timepoint = factor(timepoint_levels, levels = timepoint_levels)
  )

# ============================================================
# 8. STATS: KEGG-STYLE WILCOXON ON log2_intensity
# ============================================================

hadeg_panel_stats <- hadeg_scaffold %>%
  left_join(
    expr_hadeg,
    by = c(
      "HADEG_gene",
      "HADEG_function",
      "HADEG_class",
      "HADEG_node",
      "Sex",
      "Timepoint"
    )
  ) %>%
  group_by(HADEG_gene, HADEG_function, HADEG_class, HADEG_node, Sex, Timepoint) %>%
  group_modify(~{
    dat <- .x %>%
      filter(!is.na(log2_intensity), Group %in% c("Control", "Dosed"))
    
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
  ungroup() %>%
  group_by(HADEG_node, Sex) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    signif_label = p_to_label(p_adj)
  ) %>%
  ungroup()

# ---- KEGG-style y-position per node x sex, shared across timepoints ----
y_positions <- expr_hadeg %>%
  group_by(HADEG_node, Sex) %>%
  summarise(
    y_max = max(log2_intensity, na.rm = TRUE),
    y_min = min(log2_intensity, na.rm = TRUE),
    y_range = y_max - y_min,
    y_pos = y_max + ifelse(y_range == 0, 0.3, y_range * 0.12),
    missing_y = y_min + ifelse(y_range == 0, 0.3, y_range * 0.50),
    .groups = "drop"
  )

hadeg_panel_stats <- hadeg_panel_stats %>%
  left_join(
    y_positions,
    by = c("HADEG_node", "Sex")
  ) %>%
  mutate(
    y_pos = ifelse(is.na(y_pos), 1, y_pos),
    missing_y = ifelse(is.na(missing_y), 1, missing_y)
  )

write_csv(
  hadeg_panel_stats,
  file.path(out_dir, "tables", "HADEG_functional_node_protein_level_wilcoxon_BH_stats_STRICT.csv")
)

# ============================================================
# 9. SUMMARY TABLES
# ============================================================

hadeg_plot_summary <- expr_hadeg %>%
  group_by(HADEG_gene, HADEG_function, HADEG_class, HADEG_node, Sex, Timepoint, Group) %>%
  summarise(
    n_points = n(),
    n_unique_proteins = n_distinct(Protein),
    n_unique_organisms = n_distinct(Organism),
    median_intensity = median(log2_intensity, na.rm = TRUE),
    mean_intensity = mean(log2_intensity, na.rm = TRUE),
    proteins = paste(sort(unique(Protein)), collapse = "; "),
    organisms = paste(sort(unique(Organism)), collapse = "; "),
    .groups = "drop"
  )

write_csv(
  hadeg_plot_summary,
  file.path(out_dir, "tables", "HADEG_functional_node_protein_level_plot_summary_STRICT.csv")
)

hadeg_coverage <- expr_hadeg %>%
  distinct(
    HADEG_gene, HADEG_function, HADEG_class, HADEG_node,
    Sex, Timepoint, Protein, Organism
  ) %>%
  group_by(HADEG_gene, HADEG_function, HADEG_class, HADEG_node, Sex, Timepoint) %>%
  summarise(
    n_proteins = n_distinct(Protein),
    n_organisms = n_distinct(Organism),
    organisms = paste(sort(unique(Organism)), collapse = "; "),
    proteins = paste(sort(unique(Protein)), collapse = "; "),
    .groups = "drop"
  )

write_csv(
  hadeg_coverage,
  file.path(out_dir, "tables", "HADEG_functional_node_coverage_STRICT.csv")
)


# ============================================================
# 9B. SUPPLEMENTARY INFORMATION FOR ALL PLOTTED HADEG FUNCTIONS
# ============================================================

# This section exports supplementary information for every HADEG
# function included in target_hadeg_functions:
# pcaF, pcaI, ahpC, and est.

all_hadeg_analysis_parameters <- tibble(
  Parameter = c(
    "Figures",
    "Included HADEG genes",
    "Included HADEG functions",
    "Analysis level",
    "Input expression files",
    "Expression value",
    "Comparison",
    "Statistical test",
    "Multiple testing correction",
    "Correction grouping",
    "Significance labels",
    "Y-axis"
  ),
  Value = c(
    "HADEG protein-level boxplots for all selected functions",
    paste(target_hadeg_functions$HADEG_gene, collapse = "; "),
    paste(
      paste0(
        target_hadeg_functions$HADEG_gene,
        " | ",
        target_hadeg_functions$HADEG_function
      ),
      collapse = "; "
    ),
    "Protein-level metaproteomics",
    paste(expr_files, collapse = "; "),
    "log2_intensity from expr_imp files",
    "Control vs Dosed within each Sex x Timepoint",
    "Wilcoxon rank-sum test",
    "Benjamini-Hochberg",
    "Adjusted within HADEG_node x Sex across timepoints",
    "ns, *, **, ***, **** based on BH-adjusted p-value",
    "log2 protein intensity"
  )
)

# Protein-to-HADEG mapping for every plotted function
all_hadeg_mapping_si <- hadeg_map_target %>%
  select(
    Protein,
    HADEG_gene,
    HADEG_function,
    HADEG_class,
    HADEG_node,
    HADEG_gene_raw,
    HADEG_function_raw,
    Review_status,
    match_by
  ) %>%
  arrange(HADEG_gene, Protein)

# Complete source data used to draw every plotted boxplot
all_hadeg_boxplot_source_data_si <- expr_hadeg %>%
  select(
    Protein,
    Organism,
    HADEG_gene,
    HADEG_function,
    HADEG_class,
    HADEG_node,
    Sex,
    Timepoint,
    Group,
    Sample,
    Sample_ID,
    log2_intensity,
    source_file
  ) %>%
  arrange(
    HADEG_gene,
    Sex,
    Timepoint,
    Group,
    Protein,
    Sample_ID
  )

# Wilcoxon results and BH-adjusted p-values for every plotted function
all_hadeg_stats_si <- hadeg_panel_stats %>%
  select(
    HADEG_gene,
    HADEG_function,
    HADEG_class,
    HADEG_node,
    Sex,
    Timepoint,
    n_control,
    n_dosed,
    p_value,
    p_adj,
    signif_label
  ) %>%
  arrange(HADEG_gene, Sex, Timepoint)

# Group-level values summarized from the plotted protein intensities
all_hadeg_plot_summary_si <- hadeg_plot_summary %>%
  arrange(HADEG_gene, Sex, Timepoint, Group)

# Number and identity of proteins and organisms represented in each panel
all_hadeg_coverage_si <- hadeg_coverage %>%
  arrange(HADEG_gene, Sex, Timepoint)

# Write separate CSV files
write_csv(
  all_hadeg_mapping_si,
  file.path(
    out_dir,
    "tables",
    "All_HADEG_boxplot_protein_mapping_STRICT_SI.csv"
  )
)

write_csv(
  all_hadeg_boxplot_source_data_si,
  file.path(
    out_dir,
    "tables",
    "All_HADEG_boxplot_source_data_STRICT_SI.csv"
  )
)

write_csv(
  all_hadeg_stats_si,
  file.path(
    out_dir,
    "tables",
    "All_HADEG_boxplot_wilcoxon_BH_stats_STRICT_SI.csv"
  )
)

write_csv(
  all_hadeg_plot_summary_si,
  file.path(
    out_dir,
    "tables",
    "All_HADEG_boxplot_plot_summary_STRICT_SI.csv"
  )
)

write_csv(
  all_hadeg_coverage_si,
  file.path(
    out_dir,
    "tables",
    "All_HADEG_boxplot_coverage_STRICT_SI.csv"
  )
)

# Write one combined supplementary-information workbook
write_xlsx(
  list(
    README = all_hadeg_analysis_parameters,
    HADEG_mapping = all_hadeg_mapping_si,
    Boxplot_source_data = all_hadeg_boxplot_source_data_si,
    Wilcoxon_BH_stats = all_hadeg_stats_si,
    Plot_summary = all_hadeg_plot_summary_si,
    Coverage = all_hadeg_coverage_si
  ),
  file.path(
    out_dir,
    "tables",
    "All_HADEG_Boxplots_Supplementary_Information_STRICT.xlsx"
  )
)

# ============================================================
# 10. PLOT FUNCTION: KEGG-STYLE
# ============================================================

plot_one_hadeg_node <- function(target_node, target_sex) {
  
  df <- expr_hadeg %>%
    filter(
      HADEG_node == target_node,
      Sex == target_sex
    ) %>%
    mutate(
      Timepoint = factor(Timepoint, levels = timepoint_levels),
      Group = factor(Group, levels = group_levels)
    )
  
  scaffold_one <- hadeg_scaffold %>%
    filter(
      HADEG_node == target_node,
      Sex == target_sex
    )
  
  stat_df <- hadeg_panel_stats %>%
    filter(
      HADEG_node == target_node,
      Sex == target_sex
    ) %>%
    mutate(Timepoint = factor(Timepoint, levels = timepoint_levels))
  
  if (nrow(scaffold_one) == 0) {
    message("No scaffold found for: ", target_node, " | ", target_sex)
    return(NULL)
  }
  
  node_gene <- unique(scaffold_one$HADEG_gene)[1]
  node_function <- unique(scaffold_one$HADEG_function)[1]
  node_class <- unique(scaffold_one$HADEG_class)[1]
  
  # Add the short HADEG gene name to every plotting layer so it can
  # be displayed in a grey facet-style strip at the top of the plot.
  df <- df %>%
    mutate(Plot_label = node_gene)
  
  stat_df <- stat_df %>%
    mutate(Plot_label = node_gene)
  
  
  missing_df <- stat_df %>%
    filter(n_control == 0 & n_dosed == 0) %>%
    mutate(label = "Not detected")
  
  p <- ggplot() +
    geom_boxplot(
      data = df,
      aes(x = Timepoint, y = log2_intensity, fill = Group),
      position = position_dodge(width = 0.78),
      width = 0.64,
      alpha = 0.80,
      outlier.shape = NA,
      linewidth = 0.8
    ) +
    geom_jitter(
      data = df,
      aes(x = Timepoint, y = log2_intensity, color = Group),
      position = position_jitterdodge(
        jitter.width = 0.10,
        dodge.width = 0.78
      ),
      size = 3.0,
      alpha = 0.80
    ) +
    stat_summary(
      data = df,
      aes(x = Timepoint, y = log2_intensity, group = Group),
      fun = median,
      geom = "point",
      position = position_dodge(width = 0.78),
      shape = 23,
      size = 4.5,
      stroke = 0.9,
      fill = "white",
      color = "black"
    ) +
    geom_text(
      data = stat_df %>% filter(!(n_control == 0 & n_dosed == 0)),
      aes(x = Timepoint, y = y_pos, label = signif_label),
      inherit.aes = FALSE,
      size = 6,
      fontface = "bold"
    ) +
    geom_text(
      data = missing_df,
      aes(x = Timepoint, y = missing_y, label = label),
      inherit.aes = FALSE,
      size = 5.2,
      fontface = "bold.italic",
      color = "grey30"
    ) +
    scale_fill_manual(
      values = c("Control" = "steelblue", "Dosed" = "firebrick")
    ) +
    scale_color_manual(
      values = c("Control" = "steelblue", "Dosed" = "firebrick")
    ) +
    scale_x_discrete(drop = FALSE) +
    facet_wrap(
      ~ Plot_label,
      ncol = 1,
      labeller = label_value
    ) +
    labs(
      x = "Timepoint",
      y = expression(log[2]~protein~intensity),
      fill = "Treatment",
      color = "Treatment"
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      
      # Grey facet-style strip containing pcaF, pcaI, ahpC, or est
      strip.background = element_rect(
        fill = "grey85",
        color = "black",
        linewidth = 0.8
      ),
      strip.text = element_text(
        size = 22,
        face = "bold",
        color = "black",
        margin = margin(t = 8, r = 8, b = 8, l = 8)
      ),
      
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      legend.key.width = grid::unit(1.2, "cm"),
      legend.key.height = grid::unit(0.7, "cm"),
      
      axis.title.x = element_text(
        size = 19,
        face = "bold",
        color = "black",
        margin = margin(t = 12)
      ),
      axis.title.y = element_text(
        size = 19,
        face = "bold",
        color = "black",
        margin = margin(r = 12)
      ),
      axis.text.x = element_text(
        size = 17,
        face = "bold",
        color = "black"
      ),
      axis.text.y = element_text(
        size = 15,
        face = "bold",
        color = "black"
      ),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 1
      ),
      axis.ticks = element_line(
        linewidth = 0.8,
        color = "black"
      ),
      axis.ticks.length = grid::unit(0.22, "cm"),
      # Major horizontal grid lines
      panel.grid.major.y = element_line(
        colour = "grey90",
        linewidth = 0.5
      ),
      
      # Major vertical grid lines
      panel.grid.major.x = element_line(
        colour = "grey90",
        linewidth = 0.5
      ),
      
      # No minor grid lines
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(
        t = 14,
        r = 18,
        b = 14,
        l = 18
      )
    )
  
  file_stub <- paste(
    clean_filename(node_gene),
    clean_filename(target_sex),
    "HADEG_protein_level_boxplot_journal_ready_STRICT",
    sep = "_"
  )
  
  ggsave(
    filename = file.path(
      out_dir,
      "plots",
      "tiff",
      paste0(file_stub, ".tiff")
    ),
    plot = p,
    device = "tiff",
    width = plot_width,
    height = plot_height,
    units = "in",
    dpi = plot_dpi,
    compression = "lzw",
    bg = "white"
  )
  
  return(p)
}

# ============================================================
# 11. GENERATE ALL HADEG NODE x SEX PLOTS
# ============================================================

plot_keys <- hadeg_scaffold %>%
  distinct(HADEG_node, Sex) %>%
  arrange(HADEG_node, Sex)

hadeg_plots <- vector("list", nrow(plot_keys))
names(hadeg_plots) <- paste(plot_keys$HADEG_node, plot_keys$Sex, sep = " | ")

for (i in seq_len(nrow(plot_keys))) {
  hadeg_plots[[i]] <- plot_one_hadeg_node(
    target_node = plot_keys$HADEG_node[i],
    target_sex = plot_keys$Sex[i]
  )
}

# ============================================================
# 12. SAVE PLOT MANIFEST
# ============================================================

plot_manifest <- plot_keys %>%
  left_join(
    hadeg_scaffold %>%
      distinct(HADEG_node, HADEG_gene, HADEG_function, HADEG_class),
    by = "HADEG_node"
  ) %>%
  mutate(
    file_stub = paste(
      clean_filename(HADEG_gene),
      clean_filename(Sex),
      "HADEG_protein_level_boxplot_journal_ready_STRICT",
      sep = "_"
    ),
    tiff_file = file.path(
      "plots",
      "tiff",
      paste0(file_stub, ".tiff")
    )
  )

write_csv(
  plot_manifest,
  file.path(out_dir, "tables", "Generated_HADEG_protein_level_plot_manifest_STRICT.csv")
)

message("Done making journal-ready HADEG protein-level TIFF boxplots.")



hadeg_raw %>%
  filter(str_detect(
    str_to_lower(coalesce(as.character(.data[[function_col]]), "")),
    "pbat|als54749|esterase"
  )) %>%
  select(
    Protein = all_of(protein_col),
    HADEG_gene = all_of(gene_col),
    HADEG_description = all_of(function_col),
    Include_in_heatmap = all_of(review_col)
  )