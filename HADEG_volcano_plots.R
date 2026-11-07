

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(writexl)
  library(grid)
})

# ============================================================
# 0. USER SETTINGS
# ============================================================

root <- "/root_dir"
hadeg_dir <- file.path(root, "HADEG")
limma_dir <- file.path(root, "Litter_analysis", "Final_Litter_Blocked_Limma")

# Use the SAME HADEG annotation file you identified.
hadeg_review_file <- file.path(
  hadeg_dir,
  "02_HADEG_best_hit_per_protein_for_manual_review_STRICT.csv"
)

expected <- c("5W_M", "8W_M", "12W_M", "5W_F", "8W_F", "12W_F")

out_dir <- file.path(
  root,
  "Litter_analysis",
  "HADEG_protein_volcano_plots_Litter_Blocked_Limma"
)
plot_dir  <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")

dir.create(plot_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

plot_dpi <- 600
adj_p_cutoff <- 0.05
abs_logfc_cutoff <- 1.5

# ============================================================
# 1. HELPER FUNCTIONS
# ============================================================

clean_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

parse_dataset <- function(id) {
  tibble(
    Dataset = id,
    Timepoint = sub("_.*$", "", id),
    Sex_code = sub("^.*_", "", id),
    Sex_label = ifelse(sub("^.*_", "", id) == "M", "Males", "Females")
  )
}

clean_text <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\u00A0", " ")
  x <- str_replace_all(x, "_", " ")
  x <- str_replace_all(x, "\\s+", " ")
  str_squish(x)
}

clean_hadeg_description <- function(x) {
  x <- clean_text(x)
  x <- str_replace(x, "\\s+Bacillus subtilis\\s+RB14$", "")
  x <- str_replace(x, "\\s+Bacillus subtilis$", "")
  x <- str_replace(x, "^[A-Z0-9]+_[A-Z0-9]+\\s+", "")
  str_squish(x)
}

make_distinct_colors <- function(n) {
  if (n == 0) return(character(0))
  base_colors <- c(
    "#0072B2", "#D55E00", "#009E73", "#CC79A7",
    "#E69F00", "#56B4E9", "#000000", "#882255",
    "#44AA99", "#AA4499", "#117733", "#332288",
    "#DDCC77", "#88CCEE", "#999933", "#661100",
    "#6699CC", "#AA4466", "#4477AA", "#228833",
    "#EE6677", "#CCBB44", "#66CCEE", "#AA3377"
  )
  if (n <= length(base_colors)) return(base_colors[seq_len(n)])
  n_extra <- n - length(base_colors)
  extra_colors <- grDevices::hcl(
    h = seq(0, 360, length.out = n_extra + 1)[seq_len(n_extra)],
    c = 75,
    l = 55
  )
  c(base_colors, extra_colors)
}

safe_select <- function(df, requested_columns) {
  available_columns <- intersect(requested_columns, names(df))
  df %>% select(all_of(available_columns))
}

read_one_litter_limma <- function(id) {
  fp <- file.path(limma_dir, paste0(id, "_Litter_Blocked_Protein_Results.xlsx"))
  if (!file.exists(fp)) {
    stop("Missing litter-blocked limma workbook: ", fp)
  }
  sheets <- readxl::excel_sheets(fp)
  if (!("Protein_Statistics" %in% sheets)) {
    stop(id, ": workbook is missing Protein_Statistics sheet")
  }

  meta <- parse_dataset(id)

  df <- read_xlsx(fp, sheet = "Protein_Statistics")
  required <- c("Protein", "log2FC", "p_value", "BH_FDR")
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop(id, ": Protein_Statistics missing columns: ", paste(missing_cols, collapse = ", "))
  }

  df %>%
    mutate(
      Protein = str_squish(as.character(Protein)),
      log2FC = suppressWarnings(as.numeric(log2FC)),
      p_value = suppressWarnings(as.numeric(p_value)),
      adj_P_value = suppressWarnings(as.numeric(BH_FDR)),
      source_file = fp,
      Dataset = id,
      Timepoint = meta$Timepoint[[1]],
      Sex = meta$Sex_code[[1]],
      Sex_label = meta$Sex_label[[1]]
    ) %>%
    mutate(
      Significance = case_when(
        is.finite(adj_P_value) & adj_P_value < adj_p_cutoff &
          is.finite(log2FC) & log2FC >= abs_logfc_cutoff ~ "Up",
        is.finite(adj_P_value) & adj_P_value < adj_p_cutoff &
          is.finite(log2FC) & log2FC <= -abs_logfc_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )
}

# ============================================================
# 2. CHECK INPUT FILES
# ============================================================

missing_workbooks <- expected[!file.exists(file.path(limma_dir, paste0(expected, "_Litter_Blocked_Protein_Results.xlsx")))]
if (length(missing_workbooks) > 0) {
  stop(
    "The following litter-blocked workbooks were not found:\n",
    paste(missing_workbooks, collapse = "\n")
  )
}

if (!file.exists(hadeg_review_file)) {
  stop("The HADEG annotation file was not found:\n", hadeg_review_file)
}

# ============================================================
# 3. READ SAME HADEG ANNOTATIONS
# ============================================================

hadeg_review_raw <- read_csv(hadeg_review_file, show_col_types = FALSE)

required_hadeg_columns <- c(
  "Protein", "Manual_Broad_Class", "Manual_Subclass",
  "HADEG_gene", "HADEG_description"
)
missing_hadeg_columns <- setdiff(required_hadeg_columns, names(hadeg_review_raw))
if (length(missing_hadeg_columns) > 0) {
  stop(
    "The HADEG annotation file is missing required columns: ",
    paste(missing_hadeg_columns, collapse = ", ")
  )
}

hadeg_review <- hadeg_review_raw %>%
  mutate(
    Protein = str_squish(as.character(Protein)),
    HADEG_gene = case_when(
      is.na(HADEG_gene) | str_squish(as.character(HADEG_gene)) == "" ~ "Unknown gene",
      TRUE ~ str_squish(as.character(HADEG_gene))
    ),
    HADEG_description = case_when(
      is.na(HADEG_description) | str_squish(as.character(HADEG_description)) == "" ~ "Unknown description",
      TRUE ~ clean_hadeg_description(HADEG_description)
    ),
    Manual_Broad_Class = case_when(
      is.na(Manual_Broad_Class) | str_squish(as.character(Manual_Broad_Class)) == "" ~ "Needs manual review",
      TRUE ~ clean_text(Manual_Broad_Class)
    ),
    Manual_Subclass = case_when(
      is.na(Manual_Subclass) | str_squish(as.character(Manual_Subclass)) == "" ~ "Needs manual review",
      TRUE ~ clean_text(Manual_Subclass)
    ),
    HADEG_functional_class = paste0(HADEG_gene, " | ", HADEG_description),
    HADEG_functional_class = str_squish(str_replace_all(HADEG_functional_class, "\\s+", " ")),
    HADEG_functional_label = HADEG_functional_class,
    Annotation_review_status = case_when(
      str_to_lower(Manual_Broad_Class) == "needs manual review" |
        str_to_lower(Manual_Subclass) == "needs manual review" ~ "Needs manual review",
      TRUE ~ "Reviewed"
    )
  ) %>%
  filter(!is.na(Protein), Protein != "")

hadeg_duplicate_check <- hadeg_review %>%
  count(Protein, name = "n_HADEG_rows") %>%
  filter(n_HADEG_rows > 1)

if (nrow(hadeg_duplicate_check) > 0) {
  write_csv(hadeg_duplicate_check, file.path(table_dir, "HADEG_duplicate_annotation_check.csv"))
  message(
    nrow(hadeg_duplicate_check),
    " proteins had more than one HADEG annotation row. The first best-hit row was retained."
  )
}

hadeg_review <- hadeg_review %>% distinct(Protein, .keep_all = TRUE)
if (nrow(hadeg_review) == 0) stop("No HADEG proteins remained after reading the annotation file.")

write_csv(hadeg_review, file.path(table_dir, "HADEG_membership_ALL_proteins.csv"))

# ============================================================
# 4. READ NEW LITTER-BLOCKED LIMMA RESULTS
# ============================================================

limma_all <- map_dfr(expected, read_one_litter_limma)

# ============================================================
# 5. JOIN LIMMA STATISTICS TO HADEG ANNOTATIONS
# ============================================================

hadeg_volcano_df <- limma_all %>%
  inner_join(hadeg_review, by = "Protein", suffix = c("", "_HADEG")) %>%
  mutate(
    is_significant = Significance %in% c("Up", "Down"),
    adj_P_value_plot = case_when(
      is.na(adj_P_value) ~ NA_real_,
      adj_P_value <= 0 ~ .Machine$double.xmin,
      TRUE ~ adj_P_value
    ),
    neg_log10_adj_p = -log10(adj_P_value_plot),
    Timepoint = factor(Timepoint, levels = c("5W", "8W", "12W")),
    Sex = factor(Sex, levels = c("M", "F")),
    Sex_label = factor(Sex_label, levels = c("Males", "Females")),
    Volcano_Color_Class = if_else(is_significant, HADEG_functional_class, "Not significant")
  ) %>%
  filter(!is.na(log2FC), !is.na(neg_log10_adj_p))

if (nrow(hadeg_volcano_df) == 0) {
  stop("No proteins matched between the new litter-blocked limma results and the HADEG annotation file.")
}

write_csv(hadeg_volcano_df, file.path(table_dir, "HADEG_volcano_ALL_mapped_proteins.csv"))
write_xlsx(hadeg_volcano_df, file.path(table_dir, "HADEG_volcano_ALL_mapped_proteins.xlsx"))

# QC to show what would happen if significance were based on BH_FDR alone.
significance_check <- hadeg_volcano_df %>%
  mutate(significant_from_adj_p_only = adj_P_value < adj_p_cutoff,
         differs_from_plot_rule = is_significant != significant_from_adj_p_only) %>%
  filter(differs_from_plot_rule)
write_csv(significance_check, file.path(table_dir, "CHECK_plot_rule_vs_adjP_only.csv"))

# ============================================================
# 6. STABLE COLORS FOR SIGNIFICANT FUNCTIONAL CLASSES
# ============================================================

significant_classes <- hadeg_volcano_df %>%
  filter(is_significant) %>%
  distinct(HADEG_functional_class) %>%
  arrange(HADEG_functional_class) %>%
  pull(HADEG_functional_class)

class_colors <- make_distinct_colors(length(significant_classes))
names(class_colors) <- significant_classes
color_values <- c(class_colors, "Not significant" = "grey75")

write_csv(
  tibble(HADEG_functional_class = names(class_colors), Hex_Color = unname(class_colors)),
  file.path(table_dir, "HADEG_functional_class_color_key.csv")
)

# ============================================================
# 7. COMMON AXIS LIMITS
# ============================================================

global_x_limit <- ceiling(max(abs(hadeg_volcano_df$log2FC), na.rm = TRUE) * 1.10)
global_y_limit <- ceiling(max(hadeg_volcano_df$neg_log10_adj_p, na.rm = TRUE) * 1.10)
if (!is.finite(global_x_limit) || global_x_limit <= 0) global_x_limit <- 1
if (!is.finite(global_y_limit) || global_y_limit <= 0) global_y_limit <- 1

# ============================================================
# 8. VOLCANO PLOT FUNCTION
# ============================================================

plot_one_sex_volcano <- function(target_sex, df) {
  plot_df <- df %>% filter(Sex == target_sex)
  if (nrow(plot_df) == 0) {
    message("No data were found for sex: ", target_sex)
    return(NULL)
  }

  sex_label <- unique(as.character(plot_df$Sex_label))[1]

  classes_this_plot <- plot_df %>%
    filter(is_significant) %>%
    distinct(HADEG_functional_class) %>%
    arrange(HADEG_functional_class) %>%
    pull(HADEG_functional_class)

  legend_breaks <- c(classes_this_plot, "Not significant")

  p <- ggplot(plot_df, aes(x = log2FC, y = neg_log10_adj_p)) +
    geom_point(
      data = plot_df %>% filter(!is_significant),
      aes(color = Volcano_Color_Class),
      size = 2.7, alpha = 0.45, shape = 16
    ) +
    geom_point(
      data = plot_df %>% filter(is_significant),
      aes(color = Volcano_Color_Class),
      size = 5.0, alpha = 1, shape = 16
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey35", linewidth = 0.55) +
    geom_hline(yintercept = -log10(adj_p_cutoff), linetype = "dotted", color = "grey35", linewidth = 0.55) +
    facet_wrap(~ Timepoint, nrow = 1, scales = "fixed", drop = FALSE) +
    scale_x_continuous(
      limits = c(-global_x_limit, global_x_limit),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, global_y_limit),
      expand = expansion(mult = c(0, 0.03))
    ) +
    scale_color_manual(
      values = color_values,
      breaks = legend_breaks,
      drop = FALSE,
      name = "HADEG functional class"
    ) +
    guides(color = guide_legend(ncol = 1, byrow = TRUE,
                                override.aes = list(size = 4.5, alpha = 1))) +
    labs(
      x = expression(log[2] * " fold change (Dosed vs Control)"),
      y = expression(-log[10] * " adjusted p-value")
    ) +
    theme_bw(base_size = 16) +
    theme(
      plot.title = element_blank(),
      strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.55),
      strip.text = element_text(face = "bold", size = 16, color = "black"),
      panel.spacing = unit(0.9, "lines"),
      panel.border = element_rect(color = "black", linewidth = 0.7),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.30),
      axis.title.x = element_text(face = "bold", size = 20, color = "black", margin = margin(t = 10)),
      axis.title.y = element_text(face = "bold", size = 20, color = "black", margin = margin(r = 10)),
      axis.text = element_text(face = "bold", size = 12, color = "black"),
      axis.ticks = element_line(color = "black", linewidth = 0.60),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 14, color = "black"),
      legend.text = element_text(face = "bold", size = 12, color = "black", lineheight = 1),
      legend.key.height = unit(0.50, "cm"),
      legend.key.width = unit(0.50, "cm"),
      legend.spacing.y = unit(0.04, "cm"),
      legend.box.spacing = unit(0.25, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      plot.margin = margin(t = 8, r = 8, b = 8, l = 8)
    )

  file_stub <- paste0(clean_filename(sex_label), "_HADEG_volcano_ALL_proteins_litter_blocked_limma")

  ggsave(file.path(plot_dir, paste0(file_stub, ".pdf")), p, width = 13, height = 7, units = "in")
  ggsave(file.path(plot_dir, paste0(file_stub, ".tiff")), p, width = 13, height = 7, units = "in", dpi = plot_dpi, compression = "lzw")
  ggsave(file.path(plot_dir, paste0(file_stub, ".png")), p, width = 22, height = 8.2, units = "in", dpi = 300)

  p
}

# ============================================================
# 9. MAKE MALE AND FEMALE VOLCANO PLOTS
# ============================================================

volcano_plot_males <- plot_one_sex_volcano("M", hadeg_volcano_df)
volcano_plot_females <- plot_one_sex_volcano("F", hadeg_volcano_df)

# ============================================================
# 10. DESCRIPTIVE SUPPLEMENTARY TABLES
# ============================================================

analysis_parameters <- tibble(
  Parameter = c(
    "Analysis type",
    "Database",
    "HADEG annotation file",
    "Input differential-expression files",
    "Included HADEG proteins",
    "Protein-level statistics source",
    "Adjusted p-value source",
    "Significance rule used for plotting",
    "Adjusted p-value reference line",
    "Comparison",
    "Analysis level",
    "Volcano x-axis",
    "Volcano y-axis",
    "Color coding",
    "Panels"
  ),
  Value = c(
    "HADEG protein-level volcano visualization",
    "HADEG",
    hadeg_review_file,
    paste(file.path(limma_dir, paste0(expected, "_Litter_Blocked_Protein_Results.xlsx")), collapse = "; "),
    "All HADEG-matched proteins from the strict HADEG annotation file",
    "Existing litter-blocked protein-level limma workbooks",
    "BH_FDR from Protein_Statistics sheet",
    paste0("BH_FDR < ", adj_p_cutoff, " and |log2FC| >= ", abs_logfc_cutoff),
    paste0("Adjusted p-value = ", adj_p_cutoff),
    "Dosed vs Control within each Sex x Timepoint",
    "Protein level",
    "limma log2FC",
    "-log10 limma adjusted p-value",
    "Significant proteins colored by HADEG gene and description; non-significant proteins grey",
    "5W, 8W, and 12W"
  )
)

membership_columns <- c(
  "Protein", "HADEG_hit", "HADEG_accession", "HADEG_gene", "HADEG_description",
  "HADEG_functional_label", "HADEG_functional_class", "Manual_Broad_Class",
  "Manual_Subclass", "Annotation_review_status", "Include_in_heatmap",
  "pident", "evalue", "bitscore", "qcovhsp", "scovhsp"
)

hadeg_membership_si <- safe_select(hadeg_review, membership_columns) %>%
  arrange(HADEG_functional_class, HADEG_gene, Protein)

source_data_columns <- c(
  "Protein", "Gene", "Function", "HADEG_gene", "HADEG_description",
  "HADEG_functional_label", "HADEG_functional_class", "Manual_Broad_Class",
  "Manual_Subclass", "Annotation_review_status", "Sex_label", "Timepoint",
  "log2FC", "p_value", "adj_P_value", "neg_log10_adj_p", "Significance",
  "is_significant", "Volcano_Color_Class", "source_file"
)

hadeg_volcano_source_data_si <- safe_select(hadeg_volcano_df, source_data_columns) %>%
  arrange(Sex_label, Timepoint, desc(is_significant), HADEG_functional_class, desc(abs(log2FC)))

significant_hadeg_proteins_si <- hadeg_volcano_source_data_si %>%
  filter(is_significant) %>%
  arrange(Sex_label, Timepoint, Significance, HADEG_functional_class, desc(abs(log2FC)))

panel_summary <- hadeg_volcano_df %>%
  group_by(Sex_label, Timepoint) %>%
  summarise(
    n_total_HADEG_mapped_proteins = n_distinct(Protein),
    n_significant_HADEG_proteins = n_distinct(Protein[is_significant]),
    n_upregulated_HADEG_proteins = n_distinct(Protein[Significance == "Up"]),
    n_downregulated_HADEG_proteins = n_distinct(Protein[Significance == "Down"]),
    n_total_hadeg_functional_classes = n_distinct(HADEG_functional_class),
    n_hadeg_functional_classes_with_significant_proteins = n_distinct(HADEG_functional_class[is_significant]),
    .groups = "drop"
  ) %>%
  arrange(Sex_label, Timepoint)

hadeg_functional_class_descriptive_summary <- hadeg_volcano_df %>%
  group_by(Sex_label, Timepoint, HADEG_functional_class) %>%
  summarise(
    n_total_proteins = n_distinct(Protein),
    n_significant_proteins = n_distinct(Protein[is_significant]),
    n_upregulated_proteins = n_distinct(Protein[Significance == "Up"]),
    n_downregulated_proteins = n_distinct(Protein[Significance == "Down"]),
    mean_log2FC = mean(log2FC, na.rm = TRUE),
    median_log2FC = median(log2FC, na.rm = TRUE),
    minimum_adj_P_value = min(adj_P_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Sex_label, Timepoint, desc(n_significant_proteins), HADEG_functional_class)

write_csv(analysis_parameters, file.path(table_dir, "HADEG_volcano_analysis_parameters_ALL_proteins.csv"))
write_csv(hadeg_membership_si, file.path(table_dir, "HADEG_membership_ALL_proteins_SI.csv"))
write_csv(hadeg_volcano_source_data_si, file.path(table_dir, "HADEG_volcano_plot_source_data_ALL_proteins_SI.csv"))
write_csv(significant_hadeg_proteins_si, file.path(table_dir, "HADEG_volcano_significant_proteins_ALL_proteins_SI.csv"))
write_csv(panel_summary, file.path(table_dir, "HADEG_volcano_panel_summary_ALL_proteins_SI.csv"))
write_csv(hadeg_functional_class_descriptive_summary, file.path(table_dir, "HADEG_functional_class_DESCRIPTIVE_summary_SI.csv"))

write_xlsx(
  list(
    README_analysis_parameters = analysis_parameters,
    HADEG_membership = hadeg_membership_si,
    HADEG_volcano_source_data = hadeg_volcano_source_data_si,
    Significant_HADEG_proteins = significant_hadeg_proteins_si,
    Panel_summary = panel_summary,
    HADEG_functional_class_descriptive_summary = hadeg_functional_class_descriptive_summary,
    Color_key = tibble(HADEG_functional_class = names(class_colors), Hex_Color = unname(class_colors))
  ),
  file.path(table_dir, "HADEG_Volcano_Litter_Blocked_Supplementary_Information.xlsx")
)

message("Done making HADEG volcano plots from litter-blocked limma statistics.")
message("Results are in: ", out_dir)
