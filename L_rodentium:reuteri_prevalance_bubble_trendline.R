
# ============================================================
# Limosilactobacillus presence + abundance trendline
# Line graph = mean relative abundance (%)
# Bubble strip = prevalence out of 4 mice
# Produces plots for Sex x Tissue
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(readr)
})

# ============================================================
# 1. USER SETTINGS
# ============================================================

input_file <- "/data/L_rodentium_stats_summary.xlsx"

out_dir <- file.path(
  dirname(input_file),
  "L_rodentium_prevalence_bubble_trendline"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "png"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tiff"), recursive = TRUE, showWarnings = FALSE)

target_species <- "Limosilactobacillus rodentium"

timepoint_levels <- c("5W", "8W", "12W")

# ============================================================
# 2. READ DATA
# ============================================================

available_sheets <- excel_sheets(input_file)

message("Available sheets:")
print(available_sheets)

# Uses the first worksheet by default.
df <- read_xlsx(
  input_file,
  sheet = available_sheets[1]
)

required_cols <- c(
  "tax_id",
  "Species",
  "Sex",
  "Timepoint",
  "Tissue",
  "mean_dosed",
  "mean_control",
  "mean_dosed_percent",
  "mean_control_percent",
  "present_dosed",
  "absent_dosed",
  "present_control",
  "absent_control",
  "prevalence_dosed_percent",
  "prevalence_control_percent"
)

missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

df <- df %>%
  mutate(
    tax_id = as.character(tax_id),
    Species = as.character(Species),
    
    Sex = case_when(
      Sex %in% c("M", "Male", "male") ~ "M",
      Sex %in% c("F", "Female", "female") ~ "F",
      TRUE ~ as.character(Sex)
    ),
    
    Timepoint = case_when(
      Timepoint %in% c("5weeks", "5wks", "5wk", "5") ~ "5W",
      Timepoint %in% c("8weeks", "8wks", "8wk", "8") ~ "8W",
      Timepoint %in% c("12weeks", "12wks", "12wk", "12") ~ "12W",
      TRUE ~ as.character(Timepoint)
    ),
    
    Timepoint = factor(
      Timepoint,
      levels = timepoint_levels
    ),
    
    Tissue = as.character(Tissue),
    
    mean_dosed_percent = as.numeric(mean_dosed_percent),
    mean_control_percent = as.numeric(mean_control_percent),
    
    present_dosed = as.numeric(present_dosed),
    absent_dosed = as.numeric(absent_dosed),
    present_control = as.numeric(present_control),
    absent_control = as.numeric(absent_control),
    
    prevalence_dosed_percent = as.numeric(prevalence_dosed_percent),
    prevalence_control_percent = as.numeric(prevalence_control_percent)
  )

# ============================================================
# 3. HELPER FUNCTIONS
# ============================================================

clean_filename <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

make_presence_bubbles <- function(present, total = 4) {
  
  present <- ifelse(is.na(present), 0, present)
  present <- max(0, min(total, round(present)))
  
  paste0(
    paste(rep("●", present), collapse = ""),
    paste(rep("○", total - present), collapse = "")
  )
}

# ============================================================
# 4. RESHAPE DATA FOR PLOTTING
# ============================================================

plot_df <- df %>%
  select(
    tax_id,
    Species,
    Sex,
    Timepoint,
    Tissue,
    mean_dosed_percent,
    mean_control_percent,
    present_dosed,
    present_control,
    prevalence_dosed_percent,
    prevalence_control_percent
  ) %>%
  
  pivot_longer(
    cols = c(
      mean_dosed_percent,
      mean_control_percent
    ),
    names_to = "mean_type",
    values_to = "mean_abundance_percent"
  ) %>%
  
  mutate(
    Treatment = case_when(
      mean_type == "mean_control_percent" ~ "Control",
      mean_type == "mean_dosed_percent" ~ "Dosed"
    ),
    
    present_n = case_when(
      Treatment == "Control" ~ present_control,
      Treatment == "Dosed" ~ present_dosed
    ),
    
    prevalence_percent = case_when(
      Treatment == "Control" ~ prevalence_control_percent,
      Treatment == "Dosed" ~ prevalence_dosed_percent
    ),
    
    Treatment = factor(
      Treatment,
      levels = c("Control", "Dosed")
    )
  ) %>%
  
  complete(
    Sex,
    Tissue,
    
    Timepoint = factor(
      timepoint_levels,
      levels = timepoint_levels
    ),
    
    Treatment = factor(
      c("Control", "Dosed"),
      levels = c("Control", "Dosed")
    ),
    
    fill = list(
      mean_abundance_percent = 0,
      present_n = 0,
      prevalence_percent = 0
    )
  ) %>%
  
  mutate(
    bubble_label = mapply(
      make_presence_bubbles,
      present_n
    ),
    
    bubble_text = paste0(
      bubble_label,
      "  ",
      round(prevalence_percent),
      "%"
    )
  )

write_csv(
  plot_df,
  file.path(
    out_dir,
    "L_rodentium_plotting_table_with_prevalence_bubbles.csv"
  )
)

# ============================================================
# 5. PLOT FUNCTION
# ============================================================

plot_one_panel <- function(target_sex, target_tissue) {
  
  sex_label <- ifelse(
    target_sex == "M",
    "Male",
    "Female"
  )
  
  one <- plot_df %>%
    filter(
      Sex == target_sex,
      Tissue == target_tissue
    ) %>%
    
    mutate(
      Timepoint = factor(
        Timepoint,
        levels = timepoint_levels
      ),
      
      Treatment = factor(
        Treatment,
        levels = c("Control", "Dosed")
      )
    )
  
  if (nrow(one) == 0) {
    warning(
      "No data for ",
      target_sex,
      " / ",
      target_tissue
    )
    
    return(NULL)
  }
  
  y_max_data <- max(
    one$mean_abundance_percent,
    na.rm = TRUE
  )
  
  if (!is.finite(y_max_data) || y_max_data == 0) {
    y_max_data <- 0.1
  }
  
  # Position the dosed prevalence strip above the control strip.
  bubble_y_dosed <- y_max_data * 1.60
  bubble_y_control <- y_max_data * 1.38
  y_top <- y_max_data * 1.90
  
  bubble_df <- one %>%
    mutate(
      bubble_y = case_when(
        Treatment == "Dosed" ~ bubble_y_dosed,
        Treatment == "Control" ~ bubble_y_control
      )
    )
  
  p <- ggplot(
    one,
    aes(
      x = Timepoint,
      y = mean_abundance_percent,
      color = Treatment,
      group = Treatment
    )
  ) +
    
    # Abundance trendlines
    geom_line(
      linewidth = 1.4
    ) +
    
    geom_point(
      size = 4.2,
      stroke = 1
    ) +
    
    # Prevalence bubble strips
    geom_text(
      data = bubble_df,
      aes(
        x = Timepoint,
        y = bubble_y,
        label = bubble_text,
        color = Treatment
      ),
      inherit.aes = FALSE,
      size = 5,
      fontface = "bold"
    ) +
    
    # Treatment labels beside bubble strips
    annotate(
      "text",
      x = 0.73,
      y = bubble_y_dosed,
      label = "Dosed",
      hjust = 1,
      size = 4.3,
      fontface = "bold",
      color = "firebrick"
    ) +
    
    annotate(
      "text",
      x = 0.73,
      y = bubble_y_control,
      label = "Control",
      hjust = 1,
      size = 4.3,
      fontface = "bold",
      color = "steelblue"
    ) +
    
    scale_color_manual(
      values = c(
        "Control" = "steelblue",
        "Dosed" = "firebrick"
      )
    ) +
    
    scale_x_discrete(
      drop = FALSE
    ) +
    
    scale_y_continuous(
      limits = c(0, y_top),
      expand = expansion(
        mult = c(0.03, 0.05)
      )
    ) +
    
    coord_cartesian(
      clip = "off"
    ) +
    
    # Titles and subtitles have been removed.
    labs(
      x = "Timepoint",
      y = "Mean relative abundance (%)",
      color = "Treatment",
      caption = "● = present in one mouse; ○ = absent in one mouse"
    ) +
    
    theme_bw(
      base_size = 15
    ) +
    
    theme(
      # Bold axis titles
      axis.title.x = element_text(
        size = 17,
        face = "bold",
        margin = margin(t = 12)
      ),
      
      axis.title.y = element_text(
        size = 17,
        face = "bold",
        margin = margin(r = 12)
      ),
      
      # Larger, bold axis values
      axis.text.x = element_text(
        size = 15,
        face = "bold",
        color = "black"
      ),
      
      axis.text.y = element_text(
        size = 14,
        face = "bold",
        color = "black"
      ),
      
      # Thicker axis lines and tick marks
      axis.line = element_line(
        linewidth = 0.9,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.9,
        color = "black"
      ),
      
      axis.ticks.length = unit(
        0.22,
        "cm"
      ),
      
      # Larger, bold legend
      legend.position = "right",
      
      legend.title = element_text(
        size = 15,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 14,
        face = "bold"
      ),
      
      legend.key.height = unit(
        0.8,
        "cm"
      ),
      
      # Larger, bold prevalence explanation
      plot.caption = element_text(
        size = 13,
        face = "bold",
        hjust = 0,
        margin = margin(t = 16)
      ),
      
      # Cleaner journal-style panel
      panel.grid.major = element_line(
        linewidth = 0.35,
        color = "grey85"
      ),
      
      panel.grid.minor = element_blank(),
      
      panel.border = element_rect(
        linewidth = 0.9,
        color = "black",
        fill = NA
      ),
      
      # More room around the plot
      plot.margin = margin(
        t = 20,
        r = 35,
        b = 20,
        l = 35
      )
    )
  
  file_stub <- paste0(
    clean_filename(target_species),
    "_",
    clean_filename(sex_label),
    "_",
    clean_filename(target_tissue),
    "_abundance_prevalence_bubble_trendline"
  )
  
  # Larger PNG
  ggsave(
    filename = file.path(
      out_dir,
      "png",
      paste0(file_stub, ".png")
    ),
    plot = p,
    width = 11,
    height = 7,
    units = "in",
    dpi = 600,
    bg = "white"
  )
  
  # Vector PDF
  ggsave(
    filename = file.path(
      out_dir,
      "pdf",
      paste0(file_stub, ".pdf")
    ),
    plot = p,
    width = 11,
    height = 7,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
  
  # Publication-quality TIFF
  ggsave(
    filename = file.path(
      out_dir,
      "tiff",
      paste0(file_stub, ".tiff")
    ),
    plot = p,
    width = 11,
    height = 7,
    units = "in",
    dpi = 600,
    compression = "lzw",
    bg = "white"
  )
  
  return(p)
}

# ============================================================
# 6. MAKE ALL PLOTS
# ============================================================

all_panels <- plot_df %>%
  distinct(
    Sex,
    Tissue
  ) %>%
  arrange(
    Tissue,
    Sex
  )

plots <- list()

for (i in seq_len(nrow(all_panels))) {
  
  sx <- all_panels$Sex[i]
  ti <- all_panels$Tissue[i]
  
  message(
    "Plotting: ",
    sx,
    " / ",
    ti
  )
  
  plots[[paste(sx, ti, sep = "_")]] <-
    plot_one_panel(
      target_sex = sx,
      target_tissue = ti
    )
}

# Display plots in the RStudio viewer.
plots

message(
  "Done. Files saved to: ",
  out_dir
)
