# ============================================================
# CARBOHYDRATE METABOLISM VOLCANO PLOTS
# KOfam-first KO function labels
# Fixed x/y axes across timepoint panels within each plot
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(readxl)
  library(ggplot2)
  library(purrr)
  library(grid)
  library(writexl)
})

# ============================================================
# 0. SETTINGS
# ============================================================

setwd("data/proteomics/carbohydrate_metabolism_pathways")

out_dir <- "Carbohydrate_metabolism_volcano_plots"
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)

plot_dpi <- 600
adj_p_cutoff <- 0.05
logfc_cutoff <- 0

# ============================================================
# 1. LIMMA FILES
# ============================================================

limma_files <- c(
  "5wks_M_8_2_3_0.05.xlsx",
  "8wks_M_8_2_3_0.05.xlsx",
  "12wks_M_8_2_3_0.05.xlsx",
  "5wks_F_8_2_3_0.05.xlsx",
  "8wks_F_8_2_3_0.05.xlsx",
  "12wks_F_8_2_3_0.05.xlsx"
)

# ============================================================
# 2. CARBOHYDRATE PATHWAYS
# ============================================================

carbohydrate_pathways <- tibble(
  pathway_group = "Carbohydrate metabolism",
  pathway_id = c("map00051", "map00052", "map00500"),
  pathway_short = c(
    "Fructose and mannose metabolism",
    "Galactose metabolism",
    "Starch and sucrose metabolism"
  )
)

# ============================================================
# 3. HELPER FUNCTIONS
# ============================================================

clean_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

extract_timepoint <- function(x) {
  case_when(
    str_detect(x, "^5wks_")  ~ "5W",
    str_detect(x, "^8wks_")  ~ "8W",
    str_detect(x, "^12wks_") ~ "12W",
    TRUE ~ NA_character_
  )
}

extract_sex <- function(x) {
  case_when(
    str_detect(x, "_M_") ~ "Males",
    str_detect(x, "_F_") ~ "Females",
    TRUE ~ NA_character_
  )
}

clean_kofam_eggnog_function <- function(x) {
  x %>%
    as.character() %>%
    str_replace_all("\u00A0", " ") %>%
    str_replace("^\\s*[0-9.]+e[-+]?[0-9]+\\s+", "") %>%
    str_replace("^\\s*[0-9]+\\s+", "") %>%
    str_replace("\\s*\\[EC:[^\\]]+\\]", "") %>%
    str_replace("\\s*;.*$", "") %>%
    str_replace("^\\s*Belongs to the\\s+", "") %>%
    str_replace("\\s+family\\.?$", " family") %>%
    str_replace_all("\\s+", " ") %>%
    str_squish()
}

make_distinct_colors <- function(n) {
  
  okabe_ito <- c(
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#000000"
  )
  
  dark2 <- c(
    "#1B9E77",
    "#D95F02",
    "#7570B3",
    "#E7298A",
    "#66A61E",
    "#E6AB02",
    "#A6761D"
  )
  
  set1 <- c(
    "#E41A1C",
    "#377EB8",
    "#4DAF4A",
    "#984EA3",
    "#FF7F00",
    "#A65628",
    "#F781BF"
  )
  
  paired_strong <- c(
    "#1F78B4",
    "#33A02C",
    "#E31A1C",
    "#FF7F00",
    "#6A3D9A",
    "#B15928",
    "#FB9A99",
    "#FDBF6F",
    "#CAB2D6"
  )
  
  all_colors <- unique(c(okabe_ito, dark2, set1, paired_strong))
  
  if (n > length(all_colors)) {
    extra_cols <- grDevices::hcl(
      h = seq(15, 375, length.out = n - length(all_colors) + 1)[1:(n - length(all_colors))],
      c = 100,
      l = 45
    )
    all_colors <- c(all_colors, extra_cols)
  }
  
  all_colors[1:n]
}

# ============================================================
# 4. LOAD PATHWAY MEMBERSHIP
# ============================================================

expr_with_pathway <- read_csv(
  "Pathway_intermediate_tables/Expression_with_pathway_membership.csv",
  show_col_types = FALSE
)

pathway_info <- read_csv(
  "Pathway_intermediate_tables/Observed_pathway_metadata.csv",
  show_col_types = FALSE
)

if (!("pathway_name" %in% colnames(expr_with_pathway))) {
  expr_with_pathway <- expr_with_pathway %>%
    left_join(
      pathway_info %>%
        distinct(pathway_id, pathway_name, broad_category, subcategory),
      by = "pathway_id"
    )
}

pathway_membership_raw <- expr_with_pathway %>%
  filter(!is.na(pathway_id), !is.na(pathway_name)) %>%
  mutate(
    KO = str_squish(as.character(KO)),
    Source = str_squish(as.character(Source)),
    Function_from_annotation = str_squish(as.character(Function_from_annotation)),
    KO_function_clean = clean_kofam_eggnog_function(Function_from_annotation)
  ) %>%
  filter(!is.na(KO), KO != "")

ko_preferred_function <- pathway_membership_raw %>%
  filter(
    !is.na(KO_function_clean),
    KO_function_clean != "",
    !str_detect(str_to_lower(KO_function_clean), "^psort location"),
    !str_detect(str_to_lower(KO_function_clean), "^belongs to")
  ) %>%
  mutate(
    source_priority = case_when(
      str_detect(str_to_lower(Source), "kofam") ~ 1,
      str_detect(str_to_lower(Source), "eggnog") ~ 2,
      TRUE ~ 3
    )
  ) %>%
  arrange(KO, source_priority, KO_function_clean) %>%
  group_by(KO) %>%
  summarise(
    preferred_function = first(KO_function_clean),
    preferred_source = first(Source),
    all_functions_seen = paste(sort(unique(KO_function_clean)), collapse = " /// "),
    .groups = "drop"
  )

pathway_membership <- pathway_membership_raw %>%
  left_join(ko_preferred_function, by = "KO") %>%
  mutate(
    preferred_function = case_when(
      !is.na(preferred_function) & preferred_function != "" ~ preferred_function,
      !is.na(KO_function_clean) & KO_function_clean != "" ~ KO_function_clean,
      TRUE ~ "Unknown function"
    )
  ) %>%
  distinct(
    Protein,
    KO,
    Function_from_annotation,
    KO_function_clean,
    preferred_function,
    preferred_source,
    all_functions_seen,
    Source,
    pathway_id,
    pathway_name,
    broad_category,
    subcategory
  )

write_csv(
  ko_preferred_function,
  file.path(out_dir, "tables", "KO_preferred_function_labels_KOfam_first.csv")
)

write_xlsx(
  ko_preferred_function,
  file.path(out_dir, "tables", "KO_preferred_function_labels_KOfam_first.xlsx")
)

# ============================================================
# 5. READ ALL LIMMA RESULTS
# ============================================================

read_one_limma_file <- function(f) {
  read_xlsx(f) %>%
    mutate(
      source_file = f,
      Timepoint = extract_timepoint(basename(f)),
      Sex = extract_sex(basename(f))
    )
}

limma_all <- map_dfr(limma_files, read_one_limma_file)

# ============================================================
# 6. JOIN LIMMA RESULTS TO CARBOHYDRATE PATHWAYS
# ============================================================

carb_volcano_df <- limma_all %>%
  mutate(
    Significance = as.character(Significance),
    adj_P_value = suppressWarnings(as.numeric(adj_P_value)),
    log2FC = suppressWarnings(as.numeric(log2FC)),
    
    is_significant = case_when(
      !is.na(Significance) & Significance %in% c("Up", "Down") ~ TRUE,
      !is.na(adj_P_value) & adj_P_value < adj_p_cutoff & abs(log2FC) >= logfc_cutoff ~ TRUE,
      TRUE ~ FALSE
    ),
    
    direction = case_when(
      is_significant & log2FC > 0 ~ "Up",
      is_significant & log2FC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  left_join(pathway_membership, by = "Protein") %>%
  filter(!is.na(pathway_id)) %>%
  inner_join(carbohydrate_pathways, by = "pathway_id") %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("5W", "8W", "12W")),
    Sex = factor(Sex, levels = c("Males", "Females")),
    
    adj_P_value_plot = ifelse(
      is.na(adj_P_value) | adj_P_value <= 0,
      NA_real_,
      adj_P_value
    ),
    
    neg_log10_adj_p = -log10(adj_P_value_plot),
    
    preferred_function = case_when(
      !is.na(preferred_function) & preferred_function != "" ~ preferred_function,
      TRUE ~ "Unknown function"
    ),
    
    KO_plot = case_when(
      is_significant & !is.na(KO) & KO != "" ~
        paste0(str_squish(KO), " | ", preferred_function),
      is_significant ~ "Significant | Unknown KO",
      TRUE ~ "Not significant"
    ),
    
    pathway_label = paste0(pathway_name, " (", pathway_id, ")")
  ) %>%
  filter(
    !is.na(log2FC),
    !is.na(neg_log10_adj_p)
  ) %>%
  distinct(
    Protein, KO, pathway_id, pathway_name, Sex, Timepoint,
    log2FC, adj_P_value, neg_log10_adj_p,
    is_significant, direction, KO_plot,
    Function_from_annotation, KO_function_clean,
    preferred_function, preferred_source, all_functions_seen,
    .keep_all = TRUE
  )

write_csv(
  carb_volcano_df,
  file.path(out_dir, "tables", "carbohydrate_volcano_all_mapped_proteins.csv")
)

write_xlsx(
  carb_volcano_df,
  file.path(out_dir, "tables", "carbohydrate_volcano_all_mapped_proteins.xlsx")
)

# ============================================================
# 7. GLOBAL KO COLORS
# ============================================================

all_sig_kos <- carb_volcano_df %>%
  filter(is_significant) %>%
  distinct(KO_plot) %>%
  arrange(KO_plot) %>%
  pull(KO_plot)

ko_colors <- make_distinct_colors(length(all_sig_kos))
names(ko_colors) <- all_sig_kos

color_values <- c(
  ko_colors,
  "Not significant" = "grey75"
)

# ============================================================
# 8. PLOT FUNCTION
# Fixed axes across 5W, 8W, 12W within each pathway/sex plot
# ============================================================

plot_one_pathway_sex <- function(target_pathway_id, target_sex, df) {
  
  plot_df <- df %>%
    filter(
      pathway_id == target_pathway_id,
      Sex == target_sex
    )
  
  if (nrow(plot_df) == 0) {
    message("No proteins for ", target_pathway_id, " | ", target_sex)
    return(NULL)
  }
  
  pathway_nm <- unique(plot_df$pathway_name)[1]
  
  sig_kos_this_plot <- plot_df %>%
    filter(is_significant) %>%
    distinct(KO_plot) %>%
    arrange(KO_plot) %>%
    pull(KO_plot)
  
  # ----------------------------------------------------------
  # Fixed axis limits for this plot only
  # The biggest change across 5W/8W/12W sets the axis
  # ----------------------------------------------------------
  
  max_abs_x <- max(abs(plot_df$log2FC), na.rm = TRUE)
  max_y <- max(plot_df$neg_log10_adj_p, na.rm = TRUE)
  
  x_limit <- ceiling(max_abs_x * 1.10)
  y_limit <- ceiling(max_y * 1.10)
  
  if (is.na(x_limit) | x_limit == 0) x_limit <- 1
  if (is.na(y_limit) | y_limit == 0) y_limit <- 1
  
  p <- ggplot(
    plot_df,
    aes(
      x = log2FC,
      y = neg_log10_adj_p
    )
  ) +
    geom_point(
      data = plot_df %>% filter(!is_significant),
      aes(color = KO_plot),
      size = 2.1,
      alpha = 0.25
    ) +
    geom_point(
      data = plot_df %>% filter(is_significant),
      aes(color = KO_plot),
      size = 2.8,
      alpha = 0.95
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
      breaks = c(sig_kos_this_plot, "Not significant"),
      name = "KO | function"
    ) +
    guides(
      color = guide_legend(
        ncol = 1,
        override.aes = list(size = 2.4, alpha = 1)
      )
    ) +
    labs(
      title = paste0(pathway_nm, " (", target_pathway_id, ") | ", target_sex),
      x = "log2FC (Dosed vs Control)",
      y = "-log10 adjusted p-value"
    ) +
    theme_bw(base_size = 10.5) +
    theme(
      plot.title = element_text(face = "bold", size = 10.5),
      strip.background = element_rect(fill = "grey92", color = "black", linewidth = 0.35),
      strip.text = element_text(face = "bold", size = 8.5),
      panel.spacing = unit(0.7, "lines"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 7.2),
      legend.text = element_text(size = 4.7, lineheight = 0.78),
      legend.key.size = unit(0.24, "cm"),
      legend.spacing.y = unit(0.005, "cm"),
      legend.box.spacing = unit(0.08, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.title = element_text(face = "bold", size = 9.5),
      axis.text = element_text(size = 7.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      plot.margin = margin(t = 5, r = 3, b = 5, l = 5)
    )
  
  file_stub <- paste0(
    clean_filename(pathway_nm), "_",
    clean_filename(target_pathway_id), "_",
    clean_filename(target_sex),
    "_volcano_KOfam_first_fixed_axes"
  )
  
  ggsave(file.path(out_dir, "plots", paste0(file_stub, ".pdf")), p,
         width = 17, height = 5.5, units = "in")
  
  ggsave(file.path(out_dir, "plots", paste0(file_stub, ".png")), p,
         width = 17, height = 5.5, units = "in", dpi = plot_dpi)
  
  ggsave(file.path(out_dir, "plots", paste0(file_stub, ".tiff")), p,
         width = 17, height = 5.5, units = "in", dpi = plot_dpi,
         compression = "lzw")
  
  return(p)
}

# ============================================================
# 9. GENERATE PLOTS
# ============================================================

plot_keys <- carb_volcano_df %>%
  distinct(pathway_id, Sex) %>%
  arrange(pathway_id, Sex)

volcano_plot_list <- list()

for (i in seq_len(nrow(plot_keys))) {
  pid <- plot_keys$pathway_id[i]
  sx  <- plot_keys$Sex[i]
  
  key_name <- paste(pid, sx, sep = " | ")
  
  volcano_plot_list[[key_name]] <- plot_one_pathway_sex(
    target_pathway_id = pid,
    target_sex = sx,
    df = carb_volcano_df
  )
}

# ============================================================
# ============================================================
# 10. SUPPLEMENTARY INFORMATION TABLES
# ============================================================

volcano_analysis_parameters <- tibble(
  Parameter = c(
    "Analysis type",
    "Input limma files",
    "Pathways included",
    "Pathway IDs",
    "Significance cutoff",
    "log2FC cutoff",
    "Comparison",
    "Adjusted p-value source",
    "Function-label priority",
    "Volcano x-axis",
    "Volcano y-axis",
    "Fixed axes"
  ),
  Value = c(
    "Protein-level carbohydrate metabolism volcano analysis",
    paste(limma_files, collapse = "; "),
    paste(unique(carbohydrate_pathways$pathway_short), collapse = "; "),
    paste(unique(carbohydrate_pathways$pathway_id), collapse = "; "),
    paste0("adj_P_value < ", adj_p_cutoff),
    paste0("absolute log2FC >= ", logfc_cutoff),
    "Dosed vs Control within each Sex x Timepoint",
    "limma adjusted p-value from input files",
    "KOfam first, then eggNOG/other annotation when needed",
    "log2FC",
    "-log10 adjusted p-value",
    "Axes fixed across 5W, 8W, and 12W within each pathway x sex plot"
  )
)

volcano_plot_source_data <- carb_volcano_df %>%
  select(
    Protein,
    KO,
    pathway_id,
    pathway_name,
    pathway_group,
    pathway_short,
    Sex,
    Timepoint,
    log2FC,
    adj_P_value,
    neg_log10_adj_p,
    is_significant,
    direction,
    KO_plot,
    Function_from_annotation,
    KO_function_clean,
    preferred_function,
    preferred_source,
    all_functions_seen,
    source_file
  ) %>%
  arrange(pathway_name, Sex, Timepoint, desc(is_significant), KO, Protein)

significant_KO_summary <- carb_volcano_df %>%
  filter(is_significant) %>%
  filter(!is.na(KO), KO != "") %>%
  group_by(
    KO,
    KO_plot,
    preferred_function,
    preferred_source
  ) %>%
  summarise(
    n_significant_proteins = n_distinct(Protein),
    significant_proteins = paste(sort(unique(Protein)), collapse = "; "),
    pathways_found = paste(sort(unique(paste0(pathway_name, " (", pathway_id, ")"))), collapse = " /// "),
    sexes_found = paste(sort(unique(as.character(Sex))), collapse = " /// "),
    timepoints_found = paste(sort(unique(as.character(Timepoint))), collapse = " /// "),
    directions_found = paste(sort(unique(direction)), collapse = " /// "),
    min_adj_P_value = min(adj_P_value, na.rm = TRUE),
    max_abs_log2FC = max(abs(log2FC), na.rm = TRUE),
    all_functions_seen = paste(sort(unique(all_functions_seen)), collapse = " /// "),
    curated_function_name = "",
    .groups = "drop"
  ) %>%
  arrange(KO)

volcano_summary_by_panel <- carb_volcano_df %>%
  group_by(
    pathway_id,
    pathway_name,
    pathway_short,
    Sex,
    Timepoint
  ) %>%
  summarise(
    n_total_mapped_proteins = n_distinct(Protein),
    n_significant_proteins = n_distinct(Protein[is_significant]),
    n_upregulated_proteins = n_distinct(Protein[direction == "Up"]),
    n_downregulated_proteins = n_distinct(Protein[direction == "Down"]),
    n_significant_KOs = n_distinct(KO[is_significant & !is.na(KO) & KO != ""]),
    .groups = "drop"
  ) %>%
  arrange(pathway_name, Sex, Timepoint)


# Keep individual files too
write_csv(
  ko_preferred_function,
  file.path(out_dir, "tables", "KO_preferred_function_labels_KOfam_first.csv")
)

write_csv(
  volcano_plot_source_data,
  file.path(out_dir, "tables", "carbohydrate_volcano_plot_source_data.csv")
)

write_csv(
  significant_KO_summary,
  file.path(out_dir, "tables", "carbohydrate_volcano_significant_KO_summary.csv")
)

write_csv(
  volcano_summary_by_panel,
  file.path(out_dir, "tables", "carbohydrate_volcano_summary_by_panel.csv")
)

write_csv(
  volcano_analysis_parameters,
  file.path(out_dir, "tables", "carbohydrate_volcano_analysis_parameters.csv")
)

# Main supplementary workbook
write_xlsx(
  list(
    README_analysis_parameters = volcano_analysis_parameters,
    KO_preferred_function_labels = ko_preferred_function,
    volcano_plot_source_data = volcano_plot_source_data,
    significant_KO_summary = significant_KO_summary,
    volcano_summary_by_panel = volcano_summary_by_panel
  ),
  file.path(out_dir, "tables", "Carbohydrate_Volcano_Supplementary_Information.xlsx")
)

message("Done making volcano plots and supplementary information.")