# ============================================================
# HADEG PROTEIN-LEVEL VOLCANO PLOTS - STRICT VERSION
# Colors = cleaned HADEG functional classes
# Panels = 5W | 8W | 12W
# One plot for males, one plot for females
# Uses Significance column: Up, Down, NS
# ============================================================

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

setwd("/data/proteomics/HADEG/volcano_plots")

#Load files obtained from the Global_Proteomics_ABPP_limma_analysis_code

files <- list(
  "5W_M"  = "5wks_M.xlsx",
  "8W_M"  = "8wks_M.xlsx",
  "12W_M" = "12wks_M.xlsx",
  "5W_F"  = "5wks_F.xlsx",
  "8W_F"  = "8wks_F.xlsx",
  "12W_F" = "12wks_F.xlsx"
)

hadeg_review_file <- file.path(
  "results",
  "HADEG_inspection_files",
  "02_HADEG_best_hit_per_protein_for_manual_review.csv"
)

out_dir <- "HADEG_protein_volcano_plots"
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)

plot_dpi <- 600
adj_p_cutoff <- 0.05

# ============================================================
# 1. HELPER FUNCTIONS
# ============================================================

clean_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

parse_file_label <- function(file_label) {
  parts <- strsplit(file_label, "_")[[1]]
  tibble(
    Timepoint = parts[1],
    Sex = parts[2]
  )
}

clean_hadeg_description <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_trim(x)
  
  x <- str_replace(x, "\\s+Bacillus subtilis\\s+RB14$", "")
  x <- str_replace(x, "\\s+Bacillus subtilis$", "")
  x <- str_replace(x, "^[A-Z0-9]+_[A-Z0-9]+\\s+", "")
  
  str_squish(x)
}

clean_functional_class <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "\u00A0", " ")
  x <- str_replace_all(x, "_", " ")
  x <- str_replace_all(x, "\\s+", " ")
  x <- str_squish(x)
  
  x <- case_when(
    is.na(x) | x == "" ~ "Other HADEG function",
    str_to_lower(x) %in% c("needs manual review", "review", "unknown", "na") ~ "Other HADEG function",
    TRUE ~ x
  )
  
  x
}

make_distinct_colors <- function(n) {
  
  okabe_ito <- c(
    "#E69F00", "#56B4E9", "#009E73", "#0072B2",
    "#D55E00", "#CC79A7", "#000000"
  )
  
  dark2 <- c(
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D"
  )
  
  set1 <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#A65628", "#F781BF"
  )
  
  paired_strong <- c(
    "#1F78B4", "#33A02C", "#E31A1C", "#FF7F00",
    "#6A3D9A", "#B15928", "#FB9A99", "#FDBF6F", "#CAB2D6"
  )
  
  all_cols <- unique(c(okabe_ito, dark2, set1, paired_strong))
  
  if (n > length(all_cols)) {
    extra_cols <- grDevices::hcl(
      h = seq(15, 375, length.out = n - length(all_cols) + 1)[1:(n - length(all_cols))],
      c = 100,
      l = 45
    )
    all_cols <- c(all_cols, extra_cols)
  }
  
  all_cols[1:n]
}

standardize_limma_cols <- function(df) {
  
  if ("adj.P.Val" %in% names(df) && !"adj_P_value" %in% names(df)) {
    df <- df %>% rename(adj_P_value = adj.P.Val)
  }
  
  if ("P.Value" %in% names(df) && !"p_value" %in% names(df)) {
    df <- df %>% rename(p_value = P.Value)
  }
  
  if ("logFC" %in% names(df) && !"log2FC" %in% names(df)) {
    df <- df %>% rename(log2FC = logFC)
  }
  
  df
}

# ============================================================
# 2. READ STRICT HADEG REVIEW FILE
# ============================================================

hadeg_review <- read_csv(
  hadeg_review_file,
  show_col_types = FALSE
) %>%
  mutate(
    Protein = as.character(Protein),
    Manual_Broad_Class = as.character(Manual_Broad_Class),
    Manual_Subclass = as.character(Manual_Subclass),
    Include_in_heatmap = str_to_upper(str_trim(as.character(Include_in_heatmap))),
    HADEG_gene = as.character(HADEG_gene),
    HADEG_description = clean_hadeg_description(HADEG_description),
    
    HADEG_functional_class = paste0(
      ifelse(
        is.na(HADEG_gene) | HADEG_gene == "",
        "unknown_gene",
        str_squish(HADEG_gene)
      ),
      " | ",
      ifelse(
        is.na(HADEG_description) | HADEG_description == "",
        "unknown_description",
        clean_hadeg_description(HADEG_description)
      )
    ),
    
    HADEG_functional_class = str_replace_all(
      HADEG_functional_class,
      "\\s+",
      " "
    ),
    
    HADEG_functional_class = str_squish(HADEG_functional_class),
    
    HADEG_label = paste0(
      ifelse(is.na(HADEG_gene) | HADEG_gene == "", "unknown_gene", HADEG_gene),
      " | ",
      ifelse(is.na(HADEG_description) | HADEG_description == "", "unknown_description", HADEG_description)
    )
  ) %>%
  filter(
    !is.na(Protein),
    Protein != "",
    Include_in_heatmap == "YES",
    !is.na(Manual_Broad_Class),
    Manual_Broad_Class != "Needs manual review",
    !is.na(Manual_Subclass),
    Manual_Subclass != "Needs manual review"
  ) %>%
  distinct(
    Protein,
    HADEG_functional_class,
    Manual_Broad_Class,
    Manual_Subclass,
    HADEG_gene,
    HADEG_description,
    HADEG_label,
    .keep_all = TRUE
  )

if (nrow(hadeg_review) == 0) {
  stop("No HADEG proteins remain after filtering the review file.")
}

write_csv(
  hadeg_review,
  file.path(out_dir, "tables", "HADEG_membership_used_for_volcano.csv")
)

# ============================================================
# 3. READ ALL LIMMA/VOLCANO FILES
# ============================================================

limma_all <- map2_dfr(names(files), files, function(file_label, file_path) {
  
  message("Reading: ", file_label, " | ", file_path)
  
  parsed <- parse_file_label(file_label)
  
  df <- read_xlsx(file_path) %>%
    standardize_limma_cols() %>%
    mutate(
      source_file = file_path,
      File = file_label,
      Timepoint = parsed$Timepoint,
      Sex = parsed$Sex,
      Protein = as.character(Protein)
    )
  
  required_cols <- c("Protein", "log2FC", "adj_P_value", "Significance")
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop(
      "File ", file_path, " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df
})


# ============================================================
# 4. JOIN LIMMA RESULTS WITH STRICT HADEG CLASSES
# ============================================================

hadeg_volcano_df <- limma_all %>%
  inner_join(hadeg_review, by = "Protein") %>%
  mutate(
    log2FC = suppressWarnings(as.numeric(log2FC)),
    adj_P_value = suppressWarnings(as.numeric(adj_P_value)),
    
    Significance = case_when(
      str_to_lower(Significance) == "up" ~ "Up",
      str_to_lower(Significance) == "down" ~ "Down",
      TRUE ~ "NS"
    ),
    
    is_significant = Significance %in% c("Up", "Down"),
    
    adj_P_value_plot = case_when(
      is.na(adj_P_value) ~ NA_real_,
      adj_P_value <= 0 ~ .Machine$double.xmin,
      TRUE ~ adj_P_value
    ),
    
    neg_log10_adj_p = -log10(adj_P_value_plot),
    
    Timepoint = factor(Timepoint, levels = c("5W", "8W", "12W")),
    Sex = factor(Sex, levels = c("M", "F")),
    
    Sex_label = recode(
      as.character(Sex),
      "M" = "Males",
      "F" = "Females"
    ),
    
    volcano_color_class = ifelse(
      is_significant,
      HADEG_functional_class,
      "Not significant"
    )
  ) %>%
  filter(
    !is.na(log2FC),
    !is.na(neg_log10_adj_p)
  )

if (nrow(hadeg_volcano_df) == 0) {
  stop("No HADEG proteins matched between the limma files and the strict review file.")
}

write_csv(
  hadeg_volcano_df,
  file.path(out_dir, "tables", "HADEG_volcano_all_mapped_proteins.csv")
)

write_xlsx(
  hadeg_volcano_df,
  file.path(out_dir, "tables", "HADEG_volcano_all_mapped_proteins.xlsx")
)

# ============================================================
# 5. GLOBAL COLORS
# ============================================================

sig_classes <- hadeg_volcano_df %>%
  filter(is_significant) %>%
  distinct(HADEG_functional_class) %>%
  arrange(HADEG_functional_class) %>%
  pull(HADEG_functional_class)

class_colors <- make_distinct_colors(length(sig_classes))
names(class_colors) <- sig_classes

color_values <- c(
  class_colors,
  "Not significant" = "grey75"
)

# ============================================================
# 6. PLOT FUNCTION
# ============================================================

plot_one_sex_volcano <- function(target_sex, df) {
  
  plot_df <- df %>%
    filter(Sex == target_sex)
  
  if (nrow(plot_df) == 0) {
    message("No data for sex: ", target_sex)
    return(NULL)
  }
  
  sex_label <- unique(plot_df$Sex_label)[1]
  
  sig_classes_this_plot <- plot_df %>%
    filter(is_significant) %>%
    distinct(HADEG_functional_class) %>%
    arrange(HADEG_functional_class) %>%
    pull(HADEG_functional_class)
  
  max_abs_x <- max(abs(plot_df$log2FC), na.rm = TRUE)
  max_y <- max(plot_df$neg_log10_adj_p, na.rm = TRUE)
  
  x_limit <- ceiling(max_abs_x * 1.10)
  y_limit <- ceiling(max_y * 1.10)
  
  if (is.na(x_limit) || x_limit == 0) x_limit <- 1
  if (is.na(y_limit) || y_limit == 0) y_limit <- 1
  
  p <- ggplot(
    plot_df,
    aes(
      x = log2FC,
      y = neg_log10_adj_p
    )
  ) +
    geom_point(
      data = plot_df %>% filter(!is_significant),
      aes(color = volcano_color_class),
      size = 2.0,
      alpha = 0.35
    ) +
    geom_point(
      data = plot_df %>% filter(is_significant),
      aes(color = volcano_color_class),
      size = 4.2,
      alpha = 1
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.45
    ) +
    geom_hline(
      yintercept = -log10(adj_p_cutoff),
      linetype = "dotted",
      color = "grey40",
      linewidth = 0.45
    ) +
    facet_wrap(
      ~ Timepoint,
      nrow = 1,
      scales = "fixed"
    ) +
    scale_x_continuous(
      limits = c(-x_limit, x_limit)
    ) +
    scale_y_continuous(
      limits = c(0, y_limit)
    ) +
    scale_color_manual(
      values = color_values,
      breaks = c(sig_classes_this_plot, "Not significant"),
      name = "HADEG functional class"
    ) +
    guides(
      color = guide_legend(
        ncol = 1,
        override.aes = list(size = 2.6, alpha = 1)
      )
    ) +
    labs(
      x = "log2FC (Dosed vs Control)",
      y = "-log10 adjusted p-value"
    ) +
    theme_bw(base_size = 15) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_text(size = 9.5),
      strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = 14),
      panel.spacing = unit(0.8, "lines"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(face = "bold", size = 12, lineheight = 1),
      legend.key.size = unit(0.28, "cm"),
      legend.spacing.y = unit(0.01, "cm"),
      legend.box.spacing = unit(0.1, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      plot.margin = margin(t = 6, r = 4, b = 6, l = 6)
    )
  
  file_stub <- paste0(
    clean_filename(sex_label),
    "_HADEG_protein_level_volcano_functional_classes"
  )
  
  ggsave(
    file.path(out_dir, "plots", paste0(file_stub, ".pdf")),
    p,
    width = 21,
    height = 7.5,
    units = "in"
  )
  
  ggsave(
    file.path(out_dir, "plots", paste0(file_stub, ".png")),
    p,
    width = 17,
    height = 5.8,
    units = "in",
    dpi = plot_dpi
  )
  
  ggsave(
    file.path(out_dir, "plots", paste0(file_stub, ".tiff")),
    p,
    width = 21,
    height = 7.5,
    units = "in",
    dpi = plot_dpi,
    compression = "lzw"
  )
  
  return(p)
}

# ============================================================
# 7. MAKE MALE AND FEMALE VOLCANO PLOTS
# ============================================================

volcano_plot_males <- plot_one_sex_volcano("M", hadeg_volcano_df)
volcano_plot_females <- plot_one_sex_volcano("F", hadeg_volcano_df)

# ============================================================
# ============================================================
# 8. SUPPLEMENTARY INFORMATION TABLES
# ============================================================

hadeg_volcano_analysis_parameters <- tibble(
  Parameter = c(
    "Analysis type",
    "Database",
    "HADEG review file",
    "Input limma files",
    "Strict e-value cutoff",
    "Strict bitscore cutoff",
    "Strict percent identity cutoff",
    "Strict query coverage cutoff",
    "Strict subject coverage cutoff",
    "Significance source",
    "Significance cutoff",
    "Comparison",
    "Analysis level",
    "Volcano x-axis",
    "Volcano y-axis",
    "Color coding",
    "Panels"
  ),
  Value = c(
    "HADEG protein-level volcano analysis",
    "HADEG",
    hadeg_review_file,
    paste(files, collapse = "; "),
    "evalue <= 1e-20",
    "bitscore >= 100",
    "pident >= 40",
    "qcovhsp >= 70",
    "scovhsp >= 50",
    "Significance column from limma input files",
    paste0("Adjusted p-value < ", adj_p_cutoff),
    "Dosed vs Control within each Sex x Timepoint",
    "Protein level",
    "log2FC",
    "-log10 adjusted p-value",
    "Significant HADEG-matched proteins colored by HADEG gene/description; non-significant proteins grey",
    "5W, 8W, and 12W"
  )
)

hadeg_membership_si <- hadeg_review %>%
  select(
    Protein,
    HADEG_hit,
    HADEG_accession,
    HADEG_gene,
    HADEG_description,
    HADEG_functional_class,
    Manual_Broad_Class,
    Manual_Subclass,
    Include_in_heatmap,
    pident,
    evalue,
    bitscore,
    qcovhsp,
    scovhsp
  ) %>%
  arrange(Manual_Broad_Class, Manual_Subclass, HADEG_gene, Protein)

hadeg_volcano_source_data_si <- hadeg_volcano_df %>%
  select(
    Protein,
    Peptides,
    Gene,
    Function,
    HADEG_gene,
    HADEG_description,
    HADEG_functional_class,
    Manual_Broad_Class,
    Manual_Subclass,
    Sex_label,
    Timepoint,
    log2FC,
    t_statistic,
    p_value,
    adj_P_value,
    neg_log10_adj_p,
    Significance,
    is_significant,
    volcano_color_class,
    source_file
  ) %>%
  arrange(Sex_label, Timepoint, desc(is_significant), HADEG_functional_class, desc(abs(log2FC)))

significant_hadeg_proteins_si <- hadeg_volcano_source_data_si %>%
  filter(is_significant) %>%
  arrange(Sex_label, Timepoint, Significance, HADEG_functional_class, desc(abs(log2FC)))

# Optional but useful: compact panel-level summary
hadeg_volcano_panel_summary_si <- hadeg_volcano_df %>%
  group_by(
    Sex_label,
    Timepoint
  ) %>%
  summarise(
    n_total_HADEG_mapped_proteins = n_distinct(Protein),
    n_significant_HADEG_proteins = n_distinct(Protein[is_significant]),
    n_upregulated_HADEG_proteins = n_distinct(Protein[Significance == "Up"]),
    n_downregulated_HADEG_proteins = n_distinct(Protein[Significance == "Down"]),
    n_HADEG_functional_classes = n_distinct(HADEG_functional_class),
    n_significant_HADEG_functional_classes = n_distinct(HADEG_functional_class[is_significant]),
    .groups = "drop"
  ) %>%
  arrange(Sex_label, Timepoint)

# Save individual CSVs
write_csv(
  hadeg_volcano_analysis_parameters,
  file.path(out_dir, "tables", "HADEG_volcano_analysis_parameters.csv")
)

write_csv(
  hadeg_membership_si,
  file.path(out_dir, "tables", "HADEG_membership_used_for_volcano.csv")
)

write_csv(
  hadeg_volcano_source_data_si,
  file.path(out_dir, "tables", "HADEG_volcano_plot_source_data.csv")
)

write_csv(
  significant_hadeg_proteins_si,
  file.path(out_dir, "tables", "HADEG_volcano_significant_proteins.csv")
)

write_csv(
  hadeg_volcano_panel_summary_si,
  file.path(out_dir, "tables", "HADEG_volcano_panel_summary.csv")
)

# Main supplementary workbook
write_xlsx(
  list(
    README_analysis_parameters = hadeg_volcano_analysis_parameters,
    HADEG_membership_used = hadeg_membership_si,
    HADEG_volcano_source_data = hadeg_volcano_source_data_si,
    Significant_HADEG_proteins = significant_hadeg_proteins_si,
    Panel_summary = hadeg_volcano_panel_summary_si
  ),
  file.path(out_dir, "tables", "HADEG_Volcano_Supplementary_Information.xlsx")
)

message("Done making HADEG volcano plots and supplementary information.")