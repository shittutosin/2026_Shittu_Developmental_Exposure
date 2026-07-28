# ============================================================
# LOAD PACKAGES
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(writexl)
library(tibble)
library(ggtext)
library(patchwork)
library(cowplot)
library(grid)

# ============================================================
# SET PATH
# ============================================================

path <- "data/proteomics/protein_extraction_from_L_rodentium"

# ============================================================
# FILE PAIRS, file obtained from limma analysis for all time points (colon and small intestine microbiome)
# ============================================================

file_map <- tribble(
  ~de, ~expr,
  "5wks_F.xlsx",      "expr_imp_5wks_F.xlsx",
  "5wks_M.xlsx",      "expr_imp_5wks_M.xlsx",
  "5wks_small_F.xlsx",     "expr_imp_5wks_small_F.xlsx",
  "5wks_small_M.xlsx",     "expr_imp_5wks_small_M.xlsx",
  "8wks_F.xlsx",      "expr_imp_8wks_F.xlsx",
  "8wks_M.xlsx",      "expr_imp_8wks_M.xlsx",
  "8wks_small_F.xlsx",     "expr_imp_8wks_small_F.xlsx",
  "8wks_small_M.xlsx",     "expr_imp_8wks_small_M.xlsx",
  "12wks_F.xlsx",     "expr_imp_12wks_F.xlsx",
  "12wks_M.xlsx",     "expr_imp_12wks_M.xlsx",
  "12wks_small_F.xlsx",    "expr_imp_12wks_small_F.xlsx",
  "12wks_small_M.xlsx",    "expr_imp_12wks_small_M.xlsx"
) %>%
  mutate(
    de_path = file.path(path, de),
    expr_path = file.path(path, expr)
  )

missing_files <- file_map %>%
  filter(!file.exists(de_path) | !file.exists(expr_path))

if (nrow(missing_files) > 0) {
  stop(
    "These file pairs are missing:\n",
    paste0(missing_files$de, " <-> ", missing_files$expr, collapse = "\n")
  )
}

# ============================================================
# HELPER FUNCTIONS
# ============================================================

extract_organism <- function(x) {
  x <- as.character(x)
  
  out <- case_when(
    str_detect(x, "\\[.*?\\]") ~ str_extract(x, "\\[.*?\\]"),
    str_detect(x, "Tax=.*?TaxID=") ~ str_extract(x, "Tax=.*?TaxID="),
    TRUE ~ NA_character_
  )
  
  out <- str_remove_all(out, "^\\[|\\]$")
  out <- str_remove(out, "^Tax=")
  out <- str_remove(out, "TaxID=$")
  out <- str_remove(out, " TaxID=.*$")
  out <- str_squish(out)
  out
}

clean_function <- function(x) {
  x <- as.character(x)
  x <- str_remove(x, "\\s*\\[.*?\\]\\s*$")
  x <- str_remove(x, "\\s*n=1\\s*Tax=.*$")
  x <- str_remove(x, "\\s*Tax=.*$")
  x <- str_remove(x, "^MAG:\\s*")
  x <- str_squish(x)
  x
}

get_meta <- function(filename) {
  f <- basename(filename)
  raw_tp <- str_extract(f, "^[0-9]+wks")
  
  tibble(
    Timepoint = recode(raw_tp, "5wks" = "5W", "8wks" = "8W", "12wks" = "12W"),
    Sex = case_when(
      str_detect(f, "_M_") ~ "Male",
      str_detect(f, "_F_") ~ "Female",
      TRUE ~ NA_character_
    ),
    Tissue = case_when(
      str_detect(f, "_small_") ~ "Small intestine Microbiome",
      TRUE ~ "Colon Microbiome"
    )
  )
}

make_named_palette <- function(labels) {
  labels <- unique(labels)
  n <- length(labels)
  
  cols <- colorRampPalette(c(
    "#4E6E7A", "#D39200", "#CC79A7", "#009E73",
    "#F0E442", "#0072B2", "#D55E00", "#999999",
    "#882255", "#44AA99", "#AA4499", "#DDCC77",
    "#117733", "#88CCEE", "#CC6677", "#332288"
  ))(n)
  
  names(cols) <- labels
  cols
}

# ============================================================
# LOAD ALL DATA
# ============================================================

all_data <- vector("list", nrow(file_map))

for (i in seq_len(nrow(file_map))) {
  
  de_file <- file_map$de_path[i]
  expr_file <- file_map$expr_path[i]
  
  message("Processing: ", basename(de_file))
  
  de_tbl <- read_xlsx(de_file) %>%
    mutate(Protein = as.character(Protein))
  
  expr_tbl <- read_xlsx(expr_file) %>%
    mutate(Protein = as.character(Protein))
  
  meta <- get_meta(de_file)
  
  expr_sample_cols <- c(
    "Control_1", "Control_2", "Control_3", "Control_4",
    "Dosed_1", "Dosed_2", "Dosed_3", "Dosed_4"
  )
  
  expr_sample_cols <- expr_sample_cols[expr_sample_cols %in% colnames(expr_tbl)]
  
  if (length(expr_sample_cols) == 0) {
    warning("No expected sample columns found in: ", basename(expr_file))
    next
  }
  
  merged <- de_tbl %>%
    mutate(
      Organism = extract_organism(Function),
      Protein_label = clean_function(Function)
    ) %>%
    left_join(expr_tbl, by = "Protein")
  
  long_tbl <- merged %>%
    pivot_longer(
      cols = all_of(expr_sample_cols),
      names_to = "Sample",
      values_to = "log2_intensity"
    ) %>%
    mutate(
      Group = ifelse(str_detect(Sample, "^Control"), "Control", "Dosed"),
      Timepoint = meta$Timepoint,
      Sex = meta$Sex,
      Tissue = meta$Tissue,
      Dataset = basename(de_file)
    )
  
  all_data[[i]] <- long_tbl
}

df_all <- bind_rows(all_data)

# ============================================================
# TARGET ORGANISM
# ============================================================

target <- "Limosilactobacillus reuteri subsp. rodentium"

df_org <- df_all %>%
  filter(!is.na(Organism)) %>%
  filter(str_to_lower(str_squish(Organism)) == str_to_lower(str_squish(target)))

if (nrow(df_org) == 0) {
  stop("No rows found for target organism: ", target)
}

# ============================================================
# OUTPUT FOLDER
# ============================================================

out_dir <- file.path(
  path,
  paste0(
    str_replace_all(target, "[^A-Za-z0-9]+", "_"),
    "_fixed_panel_linear_abundance_plots"
  )
)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ============================================================
# L. RODENTIUM PROTEIN SUMMARY TABLE
# Supplementary Information
# Values are derived from imputed log2 protein intensities
# ============================================================

lmp_summary <- df_org %>%
  mutate(
    linear_intensity = 2^log2_intensity
  ) %>%
  group_by(
    Timepoint,
    Sex,
    Tissue,
    Group,
    Protein_label
  ) %>%
  summarise(
    n_values = sum(!is.na(log2_intensity)),
    
    mean_log2_intensity = mean(
      log2_intensity,
      na.rm = TRUE
    ),
    
    sd_log2_intensity = sd(
      log2_intensity,
      na.rm = TRUE
    ),
    
    median_log2_intensity = median(
      log2_intensity,
      na.rm = TRUE
    ),
    
    mean_linearized_intensity = mean(
      linear_intensity,
      na.rm = TRUE
    ),
    
    sd_linearized_intensity = sd(
      linear_intensity,
      na.rm = TRUE
    ),
    
    median_linearized_intensity = median(
      linear_intensity,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    across(
      c(
        mean_log2_intensity,
        sd_log2_intensity,
        median_log2_intensity,
        mean_linearized_intensity,
        sd_linearized_intensity,
        median_linearized_intensity
      ),
      ~ ifelse(is.nan(.x), NA_real_, .x)
    )
  ) %>%
  arrange(
    Tissue,
    Sex,
    factor(Timepoint, levels = c("5W", "8W", "12W")),
    Protein_label,
    factor(Group, levels = c("Control", "Dosed"))
  )

write_xlsx(
  lmp_summary,
  file.path(
    out_dir,
    "Supplementary_Table_L_rodentium_protein_summary.xlsx"
  )
)
# ============================================================
# SUMMARISE ABUNDANCE
# ============================================================

df_abs <- df_org %>%
  mutate(linear_intensity = 2^log2_intensity) %>%
  group_by(Timepoint, Sex, Tissue, Group, Protein_label) %>%
  summarise(
    abundance_linear = mean(linear_intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(abundance_linear), !is.na(abundance_linear), abundance_linear > 0)

# ============================================================
# SETTINGS
# ============================================================

timepoint_levels <- c("5W", "8W", "12W")
group_levels <- c("Control", "Dosed")

df_abs <- df_abs %>%
  mutate(
    Timepoint = factor(Timepoint, levels = timepoint_levels),
    Group = factor(Group, levels = group_levels),
    Sex = factor(Sex, levels = c("Male", "Female"))
  )

# ============================================================
# OUTPUT FOLDER
# ============================================================

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

write_xlsx(
  df_abs %>% arrange(Sex, Tissue, Timepoint, Group, desc(abundance_linear)),
  file.path(out_dir, "protein_contribution_summary_mean_linearized_abundance.xlsx")
)

# ============================================================
# PANEL FUNCTIONS
# ============================================================

make_detected_panel <- function(df_one, tp, combo_levels, combo_palette) {
  
  df_one <- df_one %>%
    mutate(
      Protein_label = factor(Protein_label, levels = combo_levels),
      Group = factor(Group, levels = group_levels)
    )
  
  ggplot(df_one, aes(x = Group, y = abundance_linear, fill = Protein_label)) +
    geom_col(width = 0.72, color = "black", linewidth = 0.15) +
    scale_fill_manual(values = combo_palette, limits = combo_levels, breaks = combo_levels, drop = TRUE) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels = scientific_format(digits = 2), expand = expansion(mult = c(0, 0.08))) +
    labs(title = tp, x = NULL, y = "Mean linearized protein abundance", fill = "Protein") +
    theme_bw(base_size = 12) +
    theme(
      # 5W, 8W, 12W
      plot.title = element_text(
        size = 22,
        face = "bold",
        hjust = 0.5
      ),
      
      # Control and Dosed
      axis.text.x = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      
      # Y-axis title
      axis.title.y = element_text(
        size = 15,
        face = "bold",
        color = "black"
      ),
      
      # Y-axis scale numbers
      axis.text.y = element_text(
        size = 15,
        face = "bold",
        color = "black"
      ),
      
      panel.grid.major.x = element_blank(),
      legend.position = "none"
    )
}

make_not_detected_panel <- function(tp) {
  ggplot() +
    annotate("text", x = 1, y = 1, label = "Not detected", fontface = "bold.italic", size = 5) +
    labs(title = tp, x = NULL, y = NULL) +
    xlim(0, 2) +
    ylim(0, 2) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none"
    )
}

make_legend_plot <- function(
    combo_levels,
    combo_palette
) {
  
  legend_df <- tibble(
    Protein_label = factor(
      combo_levels,
      levels = combo_levels
    ),
    x = 1,
    y = seq_along(combo_levels)
  )
  
  ggplot(
    legend_df,
    aes(
      x = x,
      y = y,
      fill = Protein_label
    )
  ) +
    
    geom_tile() +
    
    scale_fill_manual(
      values = combo_palette,
      limits = combo_levels,
      breaks = combo_levels,
      drop = TRUE
    ) +
    
    labs(
      fill = "Protein"
    ) +
    
    guides(
      fill = guide_legend(
        ncol = 1,
        byrow = TRUE,
        override.aes = list(
          color = "black",
          linewidth = 0.3
        )
      )
    ) +
    
    theme_void() +
    
    theme(
      legend.position = "right",
      
      legend.title = element_text(
        size = 15,
        face = "bold",
        color = "black"
      ),
      
      # Smaller protein names
      legend.text = element_text(
        size = 11.5,
        face = "bold",
        color = "black",
        lineheight = 0.95
      ),
      
      legend.key.height = unit(0.45, "cm"),
      legend.key.width = unit(0.45, "cm"),
      
      legend.spacing.y = unit(0.08, "cm"),
      
      legend.margin = margin(
        t = 2,
        r = 2,
        b = 2,
        l = 2
      )
    )
}

# ============================================================
# MAKE PLOTS
# ============================================================

plot_index <- df_abs %>%
  distinct(Sex, Tissue) %>%
  filter(!is.na(Sex), !is.na(Tissue)) %>%
  arrange(Sex, Tissue)

for (i in seq_len(nrow(plot_index))) {
  
  sx <- plot_index$Sex[i]
  ts <- plot_index$Tissue[i]
  
  df_combo <- df_abs %>%
    filter(Sex == sx, Tissue == ts)
  
  if (nrow(df_combo) == 0) next
  
  combo_levels <- df_combo %>%
    group_by(Protein_label) %>%
    summarise(total = sum(abundance_linear, na.rm = TRUE), .groups = "drop") %>%
    filter(total > 0) %>%
    arrange(desc(total)) %>%
    pull(Protein_label)
  
  combo_palette <- make_named_palette(combo_levels)
  
  panel_list <- list()
  
  for (tp in timepoint_levels) {
    df_one <- df_combo %>% filter(as.character(Timepoint) == tp)
    
    if (nrow(df_one) > 0) {
      panel_list[[tp]] <- make_detected_panel(df_one, tp, combo_levels, combo_palette)
    } else {
      panel_list[[tp]] <- make_not_detected_panel(tp)
    }
  }
  
  panel_list <- panel_list[timepoint_levels]
  
  panel_plot <- wrap_plots(
    panel_list,
    nrow = 1,
    ncol = 3,
    widths = c(1, 1, 1)
  )
  
  legend_plot <- make_legend_plot(combo_levels, combo_palette)
  legend_grob <- cowplot::get_legend(legend_plot)
  
  plot_title <- NULL
  
  final_plot <- wrap_plots(
    panel_plot,
    patchwork::wrap_elements(legend_grob),
    nrow = 1,
    ncol = 2,
    widths = c(3, 1.2)
  ) +
    plot_annotation(
      title = plot_title,
      theme = theme(
        plot.title = ggtext::element_markdown(size = 14, hjust = 0.5)
      )
    )
  
  print(final_plot)
  
  base_name <- paste0(
    str_replace_all(target, "[^A-Za-z0-9]+", "_"), "_",
    str_replace_all(sx, "[^A-Za-z0-9]+", "_"), "_",
    str_replace_all(ts, "[^A-Za-z0-9]+", "_"),
    "_5W_8W_12W_panel_mean_linearized_abundance"
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(base_name, ".png")),
    plot = final_plot,
    width = 18,
    height = 8,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(base_name, ".tiff")),
    plot = final_plot,
    width = 18,
    height = 8,
    dpi = 300,
    compression = "lzw"
  )
}


# ============================================================
# SELECTED SINGLE PROTEIN PLOTS OVER TIME
# Uses original log2 intensities
# Saves PNG and TIFF
# ============================================================

single_protein_dir <- file.path(out_dir, "selected_single_protein_over_time_plots")
if (!dir.exists(single_protein_dir)) dir.create(single_protein_dir, recursive = TRUE)

plot_single_protein_across_time <- function(protein_name_clean,
                                            data_long = df_org,
                                            output_path = single_protein_dir,
                                            save_plot = TRUE) {
  
  df_sub <- data_long %>%
    filter(Protein_label == protein_name_clean) %>%
    mutate(
      Timepoint = factor(Timepoint, levels = timepoint_levels),
      Group = factor(Group, levels = group_levels),
      Sex = factor(Sex, levels = c("Male", "Female"))
    )
  
  if (nrow(df_sub) == 0) {
    stop("Protein label not found: ", protein_name_clean)
  }
  
  combos <- df_sub %>%
    distinct(Sex, Tissue) %>%
    filter(!is.na(Sex), !is.na(Tissue)) %>%
    arrange(Sex, Tissue)
  
  for (j in seq_len(nrow(combos))) {
    
    sx <- combos$Sex[j]
    ts <- combos$Tissue[j]
    
    one_raw <- df_sub %>%
      filter(Sex == sx, Tissue == ts)
    
    if (nrow(one_raw) == 0) next
    
    panel_df <- tidyr::expand_grid(
      Timepoint = factor(timepoint_levels, levels = timepoint_levels),
      Group = factor(group_levels, levels = group_levels)
    )
    
    one <- one_raw %>%
      right_join(panel_df, by = c("Timepoint", "Group")) %>%
      mutate(
        Sex = sx,
        Tissue = ts,
        Protein_label = protein_name_clean
      )
    
    detected_timepoints <- one_raw %>%
      filter(!is.na(log2_intensity)) %>%
      distinct(Timepoint) %>%
      pull(Timepoint) %>%
      as.character()
    
    missing_timepoints <- setdiff(timepoint_levels, detected_timepoints)
    
    p_col <- intersect(c("adj_P_value", "adj.P.Val", "adj.P.Value"), colnames(one_raw))
    
    if (length(p_col) > 0) {
      p_df <- one_raw %>%
        group_by(Timepoint) %>%
        summarise(
          p_adj = {
            vals <- unique(.data[[p_col[1]]])
            vals <- vals[!is.na(vals)]
            if (length(vals) == 0) NA_real_ else vals[1]
          },
          .groups = "drop"
        ) %>%
        mutate(
          signif_label = case_when(
            is.na(p_adj)    ~ "NA",
            p_adj <= 0.0001 ~ "****",
            p_adj <= 0.001  ~ "***",
            p_adj <= 0.01   ~ "**",
            p_adj <= 0.05   ~ "*",
            TRUE            ~ "ns"
          )
        )
    } else {
      p_df <- tibble(
        Timepoint = factor(detected_timepoints, levels = timepoint_levels),
        p_adj = NA_real_,
        signif_label = "NA"
      )
    }
    
    y_max <- max(one_raw$log2_intensity, na.rm = TRUE)
    y_min <- min(one_raw$log2_intensity, na.rm = TRUE)
    y_range <- y_max - y_min
    
    if (!is.finite(y_max) || !is.finite(y_min)) {
      y_min <- 0
      y_max <- 1
      y_range <- 1
    }
    
    y_sig <- y_max + ifelse(y_range == 0, 0.3, y_range * 0.12)
    y_top <- y_sig + ifelse(y_range == 0, 0.3, y_range * 0.10)
    
    sig_labels <- p_df %>%
      mutate(
        Timepoint = factor(Timepoint, levels = timepoint_levels),
        x = 1.5,
        y = y_sig
      )
    
    missing_labels <- tibble(
      Timepoint = factor(missing_timepoints, levels = timepoint_levels),
      x = 1.5,
      y = y_min + ((y_top - y_min) * 0.5),
      label = "Not detected"
    )
    
    p <- ggplot(one, aes(x = Group, y = log2_intensity, fill = Group)) +
      geom_boxplot(alpha = 0.75, width = 0.65, outlier.shape = NA, na.rm = TRUE) +
      geom_jitter(width = 0.12, size = 2.3, alpha = 0.9, na.rm = TRUE) +
      geom_text(
        data = sig_labels,
        aes(x = x, y = y, label = signif_label),
        inherit.aes = FALSE,
        size = 6,
        fontface = "bold"
      ) +
      geom_text(
        data = missing_labels,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        size = 5,
        fontface = "bold.italic"
      ) +
      scale_fill_manual(values = c(Control = "steelblue", Dosed = "firebrick")) +
      scale_x_discrete(drop = FALSE) +
      facet_wrap(~ Timepoint, nrow = 1, drop = FALSE) +
      labs(
        title = NULL,
        x = NULL,
        y = expression(log[2]~"protein intensity")
      ) +
      coord_cartesian(ylim = c(y_min, y_top)) +
      theme_bw(base_size = 16) +
      theme(
        
        # Remove the protein/sex/tissue title
        plot.title = element_blank(),
        
        # 5W, 8W and 12W facet labels
        strip.background = element_rect(
          fill = "grey92",
          color = "black",
          linewidth = 0.8
        ),
        
        strip.text = element_text(
          size = 20,
          face = "bold",
          color = "black",
          margin = margin(
            t = 8,
            r = 8,
            b = 8,
            l = 8
          )
        ),
        
        # Control and Dosed labels
        axis.text.x = element_text(
          size = 17,
          face = "bold",
          color = "black"
        ),
        
        # Y-axis title
        axis.title.y = element_text(
          size = 20,
          face = "bold",
          color = "black",
          margin = margin(r = 14)
        ),
        
        # Y-axis numbers
        axis.text.y = element_text(
          size = 16,
          face = "bold",
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
        
        axis.ticks.length = unit(0.18, "cm"),
        
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
        
        panel.spacing = unit(1, "lines"),
        
        legend.position = "none",
        
        plot.margin = margin(
          t = 12,
          r = 18,
          b = 15,
          l = 18
        )
      )
    print(p)
    
    if (save_plot) {
      base_name <- paste0(
        str_replace_all(target, "[^A-Za-z0-9]+", "_"), "_",
        str_replace_all(protein_name_clean, "[^A-Za-z0-9]+", "_"), "_",
        str_replace_all(sx, "[^A-Za-z0-9]+", "_"), "_",
        str_replace_all(ts, "[^A-Za-z0-9]+", "_"),
        "_single_protein_across_time"
      )
      
      ggsave(
        filename = file.path(output_path, paste0(base_name, ".png")),
        plot = p,
        width = 14,
        height = 8,
        dpi = 300
      )
      
      ggsave(
        filename = file.path(output_path, paste0(base_name, ".tiff")),
        plot = p,
        width = 14,
        height = 8,
        dpi = 300,
        compression = "lzw"
      )
    }
  }
}

# ============================================================
# EXAMPLE CALLS
# ============================================================

plot_single_protein_across_time("LPXTG-motif cell wall anchor domain protein")