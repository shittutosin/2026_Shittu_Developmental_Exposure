# ============================================================
# SELECTED KEGG PATHWAY (MAP-ID BASED) PROTEIN-LEVEL PLOTS
# WITH STATS
#
# CHANGE ADDED:
#   - If a selected pathway is not detected at a timepoint,
#     the pathway/timepoint is still shown and labeled:
#       "Not detected"
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(ggplot2)
})

# ============================================================
# 0. SETTINGS
# ============================================================

setwd("data/proteomics/selected_pathway_boxplots")

#Load pathway summary files obtained from the Pathway_summary_proteomics_heatmap code

dir.create("Selected_pathway_protein_level", showWarnings = FALSE)
dir.create(file.path("Selected_pathway_protein_level", "plots"), showWarnings = FALSE)
dir.create(file.path("Selected_pathway_protein_level", "plots", "grouped"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("Selected_pathway_protein_level", "plots", "individual"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("Selected_pathway_protein_level", "tables"), showWarnings = FALSE)

timepoint_levels <- c("5W", "8W", "12W")
group_levels <- c("Control", "Dosed")
sex_levels <- c("Males", "Females")


plot_width_grouped <- 16
plot_height_grouped <- 9

plot_width_individual <- 20
plot_height_individual <- 18

plot_dpi <- 600
# ============================================================
# 1. SELECTED PATHWAY GROUPS (MAP IDs)
# ============================================================

selected_pathways <- list(
  "SCFA metabolism" = c(
    "map00650",
    "map00640"
  ),
  
  "Tryptophan metabolism" = c(
    "map00380"
  ),
  
  "Secondary bile acids" = c(
    "map00120",
    "map04976",
    "map00121"
  ),
  
  "Xenobiotic metabolism / PAH-related" = c(
    "map00624",
    "map00626",
    "map00362",
    "map00980",
    "map00982",
    "map00983"
  ),
  
  "Carbohydrate metabolism" = c(
    "map00051",
    "map00052",
    "map00500"
  ),
  
  "Host glycan degradation (mucin-related)" = c(
    "map00511",
    "map00531",
    "map00520"
  )
)

# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

clean_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
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

# ============================================================
# 3. LOAD DATA
# ============================================================

expr_with_pathway <- read_csv(
  "Pathway_intermediate_tables/Expression_with_pathway_membership.csv",
  show_col_types = FALSE
)

pathway_info <- read_csv(
  "Pathway_intermediate_tables/Observed_pathway_metadata.csv",
  show_col_types = FALSE
)

# ============================================================
# 4. PREP DATA
# ============================================================

if (!("pathway_name" %in% colnames(expr_with_pathway))) {
  expr_with_pathway <- expr_with_pathway %>%
    left_join(
      pathway_info %>%
        distinct(pathway_id, pathway_name, broad_category, subcategory),
      by = "pathway_id"
    )
}

expr_with_pathway <- expr_with_pathway %>%
  mutate(
    Timepoint = factor(Timepoint, levels = timepoint_levels),
    Group = factor(Group, levels = group_levels),
    Sex = recode(as.character(Sex), "M" = "Males", "F" = "Females")
  ) %>%
  filter(
    !is.na(pathway_id),
    !is.na(pathway_name),
    !is.na(log2_intensity),
    Sex %in% sex_levels
  )

# ============================================================
# 5. BUILD SELECTED DATASET
# ============================================================

selected_pathway_table <- imap_dfr(selected_pathways, function(path_ids, group_name) {
  tibble(
    selected_group = group_name,
    pathway_id = path_ids
  )
})

selected_pathway_metadata <- selected_pathway_table %>%
  left_join(
    pathway_info %>%
      distinct(pathway_id, pathway_name, broad_category, subcategory),
    by = "pathway_id"
  ) %>%
  mutate(
    pathway_name = ifelse(is.na(pathway_name), pathway_id, pathway_name),
    pathway_label = pathway_name,
    selected_group = factor(selected_group, levels = names(selected_pathways))
  )

selected_expr <- expr_with_pathway %>%
  inner_join(selected_pathway_metadata, by = "pathway_id") %>%
  mutate(
    pathway_label = pathway_name.y,
    pathway_name = pathway_name.y,
    broad_category = broad_category.y,
    subcategory = subcategory.y,
    selected_group = factor(selected_group, levels = names(selected_pathways))
  ) %>%
  select(
    -pathway_name.x, -pathway_name.y,
    -broad_category.x, -broad_category.y,
    -subcategory.x, -subcategory.y
  )

# This scaffold forces every selected pathway x sex x timepoint to exist.
selected_scaffold <- selected_pathway_metadata %>%
  tidyr::crossing(
    Sex = factor(sex_levels, levels = sex_levels),
    Timepoint = factor(timepoint_levels, levels = timepoint_levels)
  )

write_csv(
  selected_expr,
  file.path("Selected_pathway_protein_level", "tables", "Selected_pathway_protein_level_long.csv")
)

selected_pathway_membership <- selected_pathway_metadata %>%
  distinct(selected_group, pathway_id, pathway_name, broad_category, subcategory) %>%
  arrange(selected_group, pathway_id)

write_csv(
  selected_pathway_membership,
  file.path("Selected_pathway_protein_level", "tables", "Selected_pathway_membership.csv")
)

# ============================================================
# 6. STATS FUNCTION
# ============================================================

compute_pathway_stats <- function(df, scaffold_df) {
  
  stats_tbl <- scaffold_df %>%
    left_join(
      df,
      by = c(
        "selected_group",
        "pathway_id",
        "pathway_name",
        "broad_category",
        "subcategory",
        "pathway_label",
        "Sex",
        "Timepoint"
      )
    ) %>%
    group_by(selected_group, Sex, pathway_id, pathway_name, pathway_label, Timepoint) %>%
    group_modify(~{
      dat <- .x %>%
        filter(!is.na(log2_intensity), !is.na(Group))
      
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
  
  stats_tbl <- stats_tbl %>%
    group_by(selected_group, Sex) %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      signif_label = p_to_label(p_adj)
    ) %>%
    ungroup()
  
  y_positions <- df %>%
    group_by(selected_group, Sex, pathway_id, pathway_name, pathway_label) %>%
    summarise(
      y_max = max(log2_intensity, na.rm = TRUE),
      y_min = min(log2_intensity, na.rm = TRUE),
      y_range = y_max - y_min,
      y_pos = y_max + ifelse(y_range == 0, 0.3, y_range * 0.12),
      missing_y = y_min + ifelse(y_range == 0, 0.3, y_range * 0.50),
      .groups = "drop"
    )
  
  stats_tbl %>%
    left_join(
      y_positions,
      by = c("selected_group", "Sex", "pathway_id", "pathway_name", "pathway_label")
    ) %>%
    mutate(
      y_pos = ifelse(is.na(y_pos), 1, y_pos),
      missing_y = ifelse(is.na(missing_y), 1, missing_y)
    )
}

all_stats <- compute_pathway_stats(selected_expr, selected_scaffold)

write_csv(
  all_stats,
  file.path("Selected_pathway_protein_level", "tables", "Selected_pathway_stats_wilcoxon.csv")
)

# ============================================================
# 7. OPTIONAL PLOT SUMMARY TABLE
# ============================================================

plot_summary <- selected_expr %>%
  group_by(selected_group, Sex, pathway_id, pathway_name, Timepoint, Group) %>%
  summarise(
    n_points = n(),
    n_unique_proteins = n_distinct(Protein),
    median_intensity = median(log2_intensity, na.rm = TRUE),
    mean_intensity = mean(log2_intensity, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  plot_summary,
  file.path("Selected_pathway_protein_level", "tables", "Selected_pathway_plot_summary.csv")
)

# ============================================================
# 8. GROUPED PLOT FUNCTION
# ============================================================

plot_group_sex <- function(group_name, target_sex, expr_df, stats_df, scaffold_df) {
  
  df <- expr_df %>%
    filter(
      selected_group == group_name,
      Sex == target_sex
    ) %>%
    mutate(
      Timepoint = factor(Timepoint, levels = c("5W", "8W", "12W")),
      Group = factor(Group, levels = c("Control", "Dosed"))
    )
  
  scaffold_one <- scaffold_df %>%
    filter(
      selected_group == group_name,
      Sex == target_sex
    )
  
  stat_df <- stats_df %>%
    filter(
      selected_group == group_name,
      Sex == target_sex
    ) %>%
    mutate(Timepoint = factor(Timepoint, levels = timepoint_levels))
  
  missing_df <- stat_df %>%
    filter(n_control == 0 & n_dosed == 0) %>%
    mutate(label = "Not detected")
  
  p <- ggplot() +
    geom_boxplot(
      data = df,
      aes(x = Timepoint, y = log2_intensity, fill = Group),
      position = position_dodge(width = 0.75),
      width = 0.60,
      alpha = 0.75,
      outlier.shape = NA
    ) +
    geom_jitter(
      data = df,
      aes(x = Timepoint, y = log2_intensity, color = Group),
      position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
      size = 1.8,
      alpha = 0.65
    ) +
    stat_summary(
      data = df,
      aes(x = Timepoint, y = log2_intensity, group = Group),
      fun = median,
      geom = "point",
      position = position_dodge(width = 0.75),
      shape = 23,
      size = 3,
      fill = "white",
      color = "black"
    ) +
    geom_text(
      data = stat_df %>% filter(!(n_control == 0 & n_dosed == 0)),
      aes(x = Timepoint, y = y_pos, label = signif_label),
      inherit.aes = FALSE,
      size = 8,
      fontface = "bold"
    ) +
    geom_text(
      data = missing_df,
      aes(x = Timepoint, y = missing_y, label = label),
      inherit.aes = FALSE,
      size = 6,
      fontface = "bold",
      color = "grey30"
    ) +
    facet_wrap(
      ~ pathway_label,
      scales = "free_y",
      drop = FALSE
    ) +
    
    scale_fill_manual(
      values = c(
        "Control" = "steelblue",
        "Dosed" = "firebrick"
      )
    ) +
    
    scale_color_manual(
      values = c(
        "Control" = "steelblue",
        "Dosed" = "firebrick"
      )
    ) +
    
    scale_x_discrete(drop = FALSE) +
    
    labs(
      title = NULL,
      x = "Timepoint",
      y = expression(Log[2]~protein~intensity),
      fill = NULL,
      color = NULL
    ) +
    
    theme_bw(base_size = 16) +
    
    theme(
      # Remove the main title such as:
      # "SCFA metabolism | Males"
      plot.title = element_blank(),
      
      # Individual pathway titles, such as:
      # "Naphthalene degradation"
      strip.background = element_rect(
        fill = "grey92",
        color = "black",
        linewidth = 0.8
      ),
      
      strip.text = element_text(
        face = "bold",
        size = 22,
        margin = margin(
          t = 7,
          r = 7,
          b = 7,
          l = 7
        )
      ),
      
      # Legend
      legend.position = "top",
      legend.title = element_blank(),
      
      legend.text = element_text(
        face = "bold",
        size = 21
      ),
      
      legend.key.size = unit(0.9, "cm"),
      
      # Axis titles
      axis.title.x = element_text(
        face = "bold",
        size = 16,
        margin = margin(t = 12)
      ),
      
      axis.title.y = element_text(
        face = "bold",
        size = 16,
        margin = margin(r = 12)
      ),
      
      # Axis tick labels
      axis.text.x = element_text(
        face = "bold",
        size = 14,
        color = "black"
      ),
      
      axis.text.y = element_text(
        face = "bold",
        size = 14,
        color = "black"
      ),
      
      # Stronger panel and axis definition
      axis.line = element_line(
        color = "black",
        linewidth = 0.7
      ),
      
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.7
      ),
      
      axis.ticks.length = unit(0.20, "cm"),
      
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      
      panel.grid.major = element_line(
        color = "grey88",
        linewidth = 0.35
      ),
      
      panel.grid.minor = element_blank(),
      
      panel.spacing = unit(1.1, "lines"),
      
      plot.margin = margin(
        t = 12,
        r = 18,
        b = 15,
        l = 18
      )
    )
  file_stub <- paste0(clean_filename(group_name), "_", clean_filename(target_sex))
  
  ggsave(
    file.path("Selected_pathway_protein_level", "plots", "grouped", paste0(file_stub, ".pdf")),
    p, width = plot_width_grouped, height = plot_height_grouped, units = "in"
  )
  
  ggsave(
    file.path("Selected_pathway_protein_level", "plots", "grouped", paste0(file_stub, ".png")),
    p, width = plot_width_grouped, height = plot_height_grouped, units = "in", dpi = plot_dpi
  )
  
  ggsave(
    file.path("Selected_pathway_protein_level", "plots", "grouped", paste0(file_stub, ".tiff")),
    p, width = plot_width_grouped, height = plot_height_grouped, units = "in", dpi = plot_dpi, compression = "lzw"
  )
  
  return(p)
}

# ============================================================
# 9. INDIVIDUAL PATHWAY PLOT FUNCTION
# ============================================================

plot_single_pathway_sex <- function(target_group, target_pathway_id, target_sex, expr_df, stats_df, scaffold_df) {
  
  df <- expr_df %>%
    filter(
      selected_group == target_group,
      pathway_id == target_pathway_id,
      Sex == target_sex
    ) %>%
    mutate(
      Timepoint = factor(Timepoint, levels = c("5W", "8W", "12W")),
      Group = factor(Group, levels = c("Control", "Dosed"))
    )
  
  scaffold_one <- scaffold_df %>%
    filter(
      selected_group == target_group,
      pathway_id == target_pathway_id,
      Sex == target_sex
    )
  
  stat_df <- stats_df %>%
    filter(
      selected_group == target_group,
      pathway_id == target_pathway_id,
      Sex == target_sex
    ) %>%
    mutate(Timepoint = factor(Timepoint, levels = timepoint_levels))
  
  pathway_nm <- unique(scaffold_one$pathway_name)[1]
  
  missing_df <- stat_df %>%
    filter(n_control == 0 & n_dosed == 0) %>%
    mutate(label = "Not detected")
  
  p <- ggplot() +
    geom_boxplot(
      data = df,
      aes(x = Timepoint, y = log2_intensity, fill = Group),
      position = position_dodge(width = 0.75),
      width = 0.60,
      alpha = 0.75,
      outlier.shape = NA
    ) +
    geom_jitter(
      data = df,
      aes(x = Timepoint, y = log2_intensity, color = Group),
      position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
      size = 2.3,
      alpha = 0.65
    ) +
    stat_summary(
      data = df,
      aes(x = Timepoint, y = log2_intensity, group = Group),
      fun = median,
      geom = "point",
      position = position_dodge(width = 0.75),
      shape = 23,
      size = 3.8,
      fill = "white",
      color = "black"
    ) +
    geom_text(
      data = stat_df %>% filter(!(n_control == 0 & n_dosed == 0)),
      aes(x = Timepoint, y = y_pos, label = signif_label),
      inherit.aes = FALSE,
      size = 8,
      fontface = "bold"
    ) +
    geom_text(
      data = missing_df,
      aes(x = Timepoint, y = missing_y, label = label),
      inherit.aes = FALSE,
      size = 6,
      fontface = "bold",
      color = "grey30"
    ) +
    facet_wrap(
      ~ pathway_label,
      scales = "free_y",
      drop = FALSE
    ) +
    
    scale_fill_manual(
      values = c(
        "Control" = "steelblue",
        "Dosed" = "firebrick"
      )
    ) +
    
    scale_color_manual(
      values = c(
        "Control" = "steelblue",
        "Dosed" = "firebrick"
      )
    ) +
    
    scale_x_discrete(drop = FALSE) +
    
    labs(
      title = NULL,
      x = "Timepoint",
      y = expression(Log[2]~protein~intensity),
      fill = NULL,
      color = NULL
    ) +
    
    theme_bw(base_size = 16) +
    
    theme(
      # Remove the broad pathway-group title
      plot.title = element_blank(),
      
      # Keep and enlarge the actual pathway title
      strip.background = element_rect(
        fill = "grey92",
        color = "black",
        linewidth = 0.8
      ),
      
      strip.text = element_text(
        face = "bold",
        size = 25,
        margin = margin(
          t = 8,
          r = 8,
          b = 8,
          l = 8
        )
      ),
      
      # Legend
      legend.position = "top",
      legend.title = element_blank(),
      
      legend.text = element_text(
        face = "bold",
        size = 21
      ),
      
      legend.key.size = unit(0.9, "cm"),
      
      # Axis titles
      axis.title.x = element_text(
        face = "bold",
        size = 22,
        margin = margin(t = 12)
      ),
      
      axis.title.y = element_text(
        face = "bold",
        size = 22,
        margin = margin(r = 15)
      ),
      
      # Axis tick labels
      axis.text.x = element_text(
        face = "bold",
        size = 20,
        color = "black"
      ),
      
      axis.text.y = element_text(
        face = "bold",
        size = 20,
        color = "black"
      ),
      
      axis.line = element_line(
        color = "black",
        linewidth = 0.7
      ),
      
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.7
      ),
      
      axis.ticks.length = unit(0.20, "cm"),
      
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      
      panel.grid.major = element_line(
        color = "grey88",
        linewidth = 0.35
      ),
      
      panel.grid.minor = element_blank(),
      
      plot.margin = margin(
        t = 15,
        r = 22,
        b = 18,
        l = 22
      )
    )
  file_stub <- paste0(
    clean_filename(target_group), "_",
    clean_filename(target_pathway_id), "_",
    clean_filename(pathway_nm), "_",
    clean_filename(target_sex)
  )
  
  ggsave(
    file.path("Selected_pathway_protein_level", "plots", "individual", paste0(file_stub, ".pdf")),
    p, width = plot_width_individual, height = plot_height_individual, units = "in"
  )
  
  ggsave(
    file.path("Selected_pathway_protein_level", "plots", "individual", paste0(file_stub, ".png")),
    p, width = plot_width_individual, height = plot_height_individual, units = "in", dpi = plot_dpi
  )
  
  ggsave(
    file.path("Selected_pathway_protein_level", "plots", "individual", paste0(file_stub, ".tiff")),
    p, width = plot_width_individual, height = plot_height_individual, units = "in", dpi = plot_dpi, compression = "lzw"
  )
  
  return(p)
}

# ============================================================
# 10. GENERATE ALL GROUPED PLOTS
# ============================================================

grouped_plot_list <- list()

for (grp in names(selected_pathways)) {
  for (sx in sex_levels) {
    grouped_plot_list[[paste(grp, sx, sep = " | ")]] <- plot_group_sex(
      group_name = grp,
      target_sex = sx,
      expr_df = selected_expr,
      stats_df = all_stats,
      scaffold_df = selected_scaffold
    )
  }
}

# ============================================================
# 11. GENERATE ALL INDIVIDUAL PATHWAY PLOTS
# ============================================================

individual_plot_list <- list()

for (grp in names(selected_pathways)) {
  for (pid in selected_pathways[[grp]]) {
    for (sx in sex_levels) {
      individual_plot_list[[paste(grp, pid, sx, sep = " | ")]] <- plot_single_pathway_sex(
        target_group = grp,
        target_pathway_id = pid,
        target_sex = sx,
        expr_df = selected_expr,
        stats_df = all_stats,
        scaffold_df = selected_scaffold
      )
    }
  }
}

# ============================================================
# 12. SAVE MANIFEST OF GENERATED GROUPED PLOTS
# ============================================================

generated_grouped_plot_table <- expand.grid(
  selected_group = names(selected_pathways),
  Sex = sex_levels,
  stringsAsFactors = FALSE
) %>%
  mutate(
    file_stub = paste0(clean_filename(selected_group), "_", clean_filename(Sex)),
    pdf_file = paste0(file_stub, ".pdf"),
    png_file = paste0(file_stub, ".png"),
    tiff_file = paste0(file_stub, ".tiff")
  )

write_csv(
  generated_grouped_plot_table,
  file.path("Selected_pathway_protein_level", "tables", "Generated_grouped_plot_file_manifest.csv")
)

# ============================================================
# 13. SAVE MANIFEST OF GENERATED INDIVIDUAL PLOTS
# ============================================================

individual_plot_manifest <- selected_pathway_metadata %>%
  tidyr::crossing(Sex = sex_levels) %>%
  mutate(
    file_stub = paste0(
      clean_filename(selected_group), "_",
      clean_filename(pathway_id), "_",
      clean_filename(pathway_name), "_",
      clean_filename(Sex)
    ),
    pdf_file = paste0(file_stub, ".pdf"),
    png_file = paste0(file_stub, ".png"),
    tiff_file = paste0(file_stub, ".tiff")
  ) %>%
  arrange(selected_group, pathway_id, Sex)

write_csv(
  individual_plot_manifest,
  file.path("Selected_pathway_protein_level", "tables", "Generated_individual_plot_file_manifest.csv")
)

message("Done.")

library(openxlsx)

# ============================================================
# 14. SAVE SUPPLEMENTARY INFORMATION WORKBOOK
# ============================================================

dir.create(file.path("Selected_pathway_protein_level", "supplementary"), showWarnings = FALSE)

supp_file <- file.path(
  "Selected_pathway_protein_level",
  "supplementary",
  "Supplementary_Selected_Pathway_Protein_Level_Tables.xlsx"
)

# -----------------------------
# Expanded stats table
# -----------------------------
expanded_stats <- selected_expr %>%
  group_by(selected_group, pathway_id, pathway_name, Sex, Timepoint, Group) %>%
  summarise(
    n_points = n(),
    n_unique_proteins = n_distinct(Protein),
    median_intensity = median(log2_intensity, na.rm = TRUE),
    mean_intensity = mean(log2_intensity, na.rm = TRUE),
    sd_intensity = sd(log2_intensity, na.rm = TRUE),
    min_intensity = min(log2_intensity, na.rm = TRUE),
    max_intensity = max(log2_intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(
      n_points,
      n_unique_proteins,
      median_intensity,
      mean_intensity,
      sd_intensity,
      min_intensity,
      max_intensity
    )
  ) %>%
  left_join(
    all_stats %>%
      select(
        selected_group, pathway_id, pathway_name, Sex, Timepoint,
        n_control, n_dosed, p_value, p_adj, signif_label
      ),
    by = c("selected_group", "pathway_id", "pathway_name", "Sex", "Timepoint")
  ) %>%
  mutate(
    difference_median_dosed_minus_control =
      median_intensity_Dosed - median_intensity_Control,
    difference_mean_dosed_minus_control =
      mean_intensity_Dosed - mean_intensity_Control,
    test_used = "Wilcoxon rank-sum test",
    p_adjust_method = "Benjamini-Hochberg FDR"
  ) %>%
  select(
    selected_group, pathway_id, pathway_name, Sex, Timepoint,
    n_control, n_dosed,
    n_unique_proteins_Control, n_unique_proteins_Dosed,
    median_intensity_Control, median_intensity_Dosed,
    mean_intensity_Control, mean_intensity_Dosed,
    sd_intensity_Control, sd_intensity_Dosed,
    min_intensity_Control, min_intensity_Dosed,
    max_intensity_Control, max_intensity_Dosed,
    difference_median_dosed_minus_control,
    difference_mean_dosed_minus_control,
    p_value, p_adj, signif_label,
    test_used, p_adjust_method
  )

# -----------------------------
# Detection / missingness table
# -----------------------------
detection_summary <- selected_scaffold %>%
  left_join(
    selected_expr %>%
      group_by(selected_group, pathway_id, pathway_name, Sex, Timepoint, Group) %>%
      summarise(
        n_detected = n(),
        n_unique_proteins = n_distinct(Protein),
        .groups = "drop"
      ),
    by = c("selected_group", "pathway_id", "pathway_name", "Sex", "Timepoint")
  ) %>%
  mutate(
    Group = as.character(Group),
    n_detected = replace_na(n_detected, 0),
    n_unique_proteins = replace_na(n_unique_proteins, 0)
  ) %>%
  select(
    selected_group, pathway_id, pathway_name, Sex, Timepoint,
    Group, n_detected, n_unique_proteins
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(n_detected, n_unique_proteins),
    values_fill = 0
  ) %>%
  mutate(
    detected_control = n_detected_Control > 0,
    detected_dosed = n_detected_Dosed > 0,
    status = case_when(
      detected_control & detected_dosed ~ "Detected in both",
      detected_control & !detected_dosed ~ "Control only",
      !detected_control & detected_dosed ~ "Dosed only",
      TRUE ~ "Not detected"
    )
  ) %>%
  select(
    selected_group, pathway_id, pathway_name, Sex, Timepoint,
    detected_control, detected_dosed,
    n_detected_Control, n_detected_Dosed,
    n_unique_proteins_Control, n_unique_proteins_Dosed,
    status
  )

# -----------------------------
# Protein contribution table
# -----------------------------
protein_contribution_table <- selected_expr %>%
  group_by(
    selected_group, pathway_id, pathway_name,
    Sex, Timepoint, Group,
    Protein, KO, Function_from_annotation
  ) %>%
  summarise(
    median_log2_intensity = median(log2_intensity, na.rm = TRUE),
    mean_log2_intensity = mean(log2_intensity, na.rm = TRUE),
    sd_log2_intensity = sd(log2_intensity, na.rm = TRUE),
    n_samples_detected = n_distinct(Sample),
    .groups = "drop"
  ) %>%
  arrange(selected_group, pathway_id, Sex, Timepoint, Group, desc(median_log2_intensity))

# -----------------------------
# Analysis settings table
# -----------------------------
analysis_settings <- tibble(
  Setting = c(
    "Input file",
    "Selected pathway source",
    "Expression value",
    "Comparison",
    "Timepoints",
    "Sexes",
    "Statistical test",
    "P-value adjustment",
    "Significance labels",
    "Missing pathway label"
  ),
  Value = c(
    "Pathway_intermediate_tables/Expression_with_pathway_membership.csv",
    "User-selected KEGG MAP IDs",
    "Protein-level log2 intensity",
    "Control vs Dosed within each pathway, sex, and timepoint",
    paste(timepoint_levels, collapse = ", "),
    paste(sex_levels, collapse = ", "),
    "Wilcoxon rank-sum test",
    "Benjamini-Hochberg FDR within selected_group and Sex",
    "**** <=0.0001; *** <=0.001; ** <=0.01; * <=0.05; ns >0.05",
    "Not detected"
  )
)

# -----------------------------
# Selected pathway group table
# -----------------------------
selected_pathway_groups <- imap_dfr(selected_pathways, function(path_ids, group_name) {
  tibble(
    selected_group = group_name,
    pathway_id = path_ids
  )
}) %>%
  left_join(
    pathway_info %>%
      distinct(pathway_id, pathway_name, broad_category, subcategory),
    by = "pathway_id"
  )

# -----------------------------
# Save Excel workbook
# -----------------------------
wb <- createWorkbook()

si_tables <- list(
  "SI_Settings" = analysis_settings,
  "Selected_pathway_groups" = selected_pathway_groups,
  "Selected_pathway_membership" = selected_pathway_membership,
  "Selected_expr_long" = selected_expr,
  "Expanded_stats" = expanded_stats,
  "Detection_summary" = detection_summary,
  "Protein_contributions" = protein_contribution_table,
  "Plot_summary" = plot_summary,
  "Grouped_plot_manifest" = generated_grouped_plot_table,
  "Individual_plot_manifest" = individual_plot_manifest
)

for (sheet_name in names(si_tables)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, si_tables[[sheet_name]])
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = 1:ncol(si_tables[[sheet_name]]), widths = "auto")
}

saveWorkbook(wb, supp_file, overwrite = TRUE)


message("Supplementary workbook saved: ", supp_file)


