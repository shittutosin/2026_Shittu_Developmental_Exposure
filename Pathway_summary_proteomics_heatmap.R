# ============================================================
# PATHWAY-SUMMARY HEATMAPS
# expr_imp files + KOfamScan + eggNOG + KEGG REST
#
# FINAL VERSION:
#   - rows = KEGG pathway summaries
#   - pathways grouped by broad KEGG category + subcategory
#   - row selection/order = within-panel Control vs Dosed difference
#   - actual heatmap values remain sample-level, not group-averaged
#   - z-score = row z-score within each dataset
#   - sample names hidden
#   - saves PDF + TIFF + PNG
# ============================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(tibble)
  library(KEGGREST)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ============================================================
# 0. USER SETTINGS
# ============================================================
setwd("data/proteomics/pathway_analysis")

#Load expr_imp files obtained from the Global_Proteomics_ABPP_limma_analysis_code
files <- list(
  "5W_M"  = "expr_imp_5wks_M.xlsx",
  "8W_M"  = "expr_imp_8wks_M.xlsx",
  "12W_M" = "expr_imp_12wks_M.xlsx",
  "5W_F"  = "expr_imp_5wks_F.xlsx",
  "8W_F"  = "expr_imp_8wks_F.xlsx",
  "12W_F" = "expr_imp_12wks_F.xlsx"
)

# -----------------------------
# Z-score choice
# Recommended for your figure:
# FALSE = row z-score within dataset
# TRUE  = global row z-score across all datasets
# -----------------------------
use_global_zscore <- FALSE

# -----------------------------
# KOfam relaxed criteria
# -----------------------------
relaxed_evalue_cutoff <- 1e-5
relaxed_score_fraction <- 0.5

# -----------------------------
# Filtering / clutter control
# -----------------------------
min_proteins_per_pathway <- 2
min_datasets_per_pathway <- 2

# Keep only top pathways per heatmap based on
# abs(mean(Dosed) - mean(Control)) within that heatmap
# Set to NULL to keep all
max_pathways_per_heatmap <- 80

# -----------------------------
# KEGG cache
# -----------------------------
use_kegg_cache <- TRUE
kegg_cache_dir <- "KEGG_cache"
dir.create(kegg_cache_dir, showWarnings = FALSE)

ko_pathway_cache_file   <- file.path(kegg_cache_dir, "KO_to_pathway_mapping.csv")
pathway_meta_cache_file <- file.path(kegg_cache_dir, "Pathway_metadata.csv")

# -----------------------------
# Output folders
# -----------------------------
dir.create("Pathway_heatmaps", showWarnings = FALSE)
dir.create("Pathway_intermediate_tables", showWarnings = FALSE)

# -----------------------------
# Plot export settings
# -----------------------------
pdf_width  <- 14
pdf_height <- 18

png_width_px  <- 4200
png_height_px <- 5400
png_res <- 300

tiff_width_in  <- 14
tiff_height_in <- 18
tiff_res <- 300
# ============================================================
# 1. HELPER FUNCTIONS
# ============================================================

parse_file_label <- function(file_label) {
  parts <- strsplit(file_label, "_")[[1]]
  list(
    Timepoint = parts[1],
    Sex = parts[2]
  )
}

get_sample_columns <- function(df) {
  grep("^(Control|Dosed)_", colnames(df), value = TRUE)
}

get_group_from_sample <- function(sample_names) {
  case_when(
    str_detect(sample_names, "^Control_") ~ "Control",
    str_detect(sample_names, "^Dosed_")   ~ "Dosed",
    TRUE ~ NA_character_
  )
}

zscore_row <- function(x) {
  mu <- mean(x, na.rm = TRUE)
  sdv <- stats::sd(x, na.rm = TRUE)
  if (is.na(sdv) || sdv == 0) {
    return(rep(0, length(x)))
  }
  (x - mu) / sdv
}

row_zscore_matrix <- function(mat) {
  t(apply(mat, 1, zscore_row))
}

parse_kegg_class <- function(class_vec) {
  if (length(class_vec) == 0 || all(is.na(class_vec))) {
    return(list(
      broad_category = "Unclassified",
      subcategory = "Unclassified"
    ))
  }
  
  class_line <- class_vec[!is.na(class_vec) & class_vec != ""][1]
  parts <- str_split(class_line, ";")[[1]]
  parts <- str_trim(parts)
  
  broad_category <- ifelse(length(parts) >= 1, parts[1], "Unclassified")
  subcategory    <- ifelse(length(parts) >= 2, parts[2], "Unclassified")
  
  if (is.na(broad_category) || broad_category == "") broad_category <- "Unclassified"
  if (is.na(subcategory) || subcategory == "") subcategory <- "Unclassified"
  
  list(
    broad_category = broad_category,
    subcategory = subcategory
  )
}

chunk_vector <- function(x, chunk_size = 10) {
  split(x, ceiling(seq_along(x) / chunk_size))
}

safe_keggGet <- function(ids, retries = 3, sleep_sec = 1) {
  for (i in seq_len(retries)) {
    out <- tryCatch(
      KEGGREST::keggGet(ids),
      error = function(e) NULL
    )
    if (!is.null(out)) return(out)
    Sys.sleep(sleep_sec * i)
  }
  return(NULL)
}

select_top_separating_pathways <- function(mat_raw, group_vector, n_keep = NULL) {
  if (is.null(n_keep) || nrow(mat_raw) <= n_keep) {
    return(rownames(mat_raw))
  }
  
  control_cols <- which(group_vector == "Control")
  dosed_cols   <- which(group_vector == "Dosed")
  
  if (length(control_cols) == 0 || length(dosed_cols) == 0) {
    return(rownames(mat_raw))
  }
  
  effect_size <- apply(mat_raw, 1, function(x) {
    control_mean <- mean(x[control_cols], na.rm = TRUE)
    dosed_mean   <- mean(x[dosed_cols], na.rm = TRUE)
    abs(dosed_mean - control_mean)
  })
  
  effect_size[is.na(effect_size)] <- -Inf
  
  names(sort(effect_size, decreasing = TRUE))[seq_len(min(n_keep, length(effect_size)))]
}

# ============================================================
# 2. LOAD eggNOG -> ALL PROTEIN-KO MEMBERSHIPS
# ============================================================

raw_lines <- readLines("eggnog.emapper.annotations")

header_line <- sub("^#", "", raw_lines[grepl("^#query", raw_lines)])
header <- strsplit(header_line, "\t")[[1]]

eggnog <- read.delim(
  "eggnog.emapper.annotations",
  sep = "\t",
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = ""
)

colnames(eggnog) <- header

eggnog_evalue_col <- intersect(
  c("evalue", "e_value", "Evalue", "seed_ortholog_evalue"),
  colnames(eggnog)
)

if (length(eggnog_evalue_col) == 0) {
  stop("Could not find an eggNOG e-value column in merged.emapper.annotations")
}
eggnog_evalue_col <- eggnog_evalue_col[1]

if (!"Description" %in% colnames(eggnog)) {
  eggnog$Description <- NA_character_
}

eggnog_all_kos <- eggnog %>%
  mutate(
    query = as.character(query),
    KEGG_ko = as.character(KEGG_ko),
    eggnog_evalue = suppressWarnings(as.numeric(.data[[eggnog_evalue_col]])),
    Description = as.character(Description)
  ) %>%
  filter(!is.na(KEGG_ko), KEGG_ko != "") %>%
  mutate(KEGG_ko = str_replace_all(KEGG_ko, "ko:", "")) %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  mutate(KEGG_ko = str_trim(KEGG_ko)) %>%
  filter(KEGG_ko != "") %>%
  transmute(
    Protein = query,
    KO = KEGG_ko,
    Function_from_annotation = Description,
    Source = "eggNOG",
    eggnog_evalue = eggnog_evalue,
    kofamscan_evalue = NA_real_,
    kofamscan_score = NA_real_
  ) %>%
  distinct()

# ============================================================
# 3. LOAD KOfam -> ALL PROTEIN-KO MEMBERSHIPS
# ============================================================

kofam_raw <- read.table(
  "kofam_annotation_output.txt",
  header = FALSE,
  comment.char = "#",
  stringsAsFactors = FALSE,
  sep = "",
  fill = TRUE,
  quote = ""
)

kofam_fixed <- kofam_raw
star_rows <- kofam_raw$V1 == "*"

kofam_fixed[star_rows, ] <- cbind(
  gene = kofam_raw$V2[star_rows],
  KO = kofam_raw$V3[star_rows],
  threshold_flag = "*",
  score = kofam_raw$V4[star_rows],
  evalue = kofam_raw$V5[star_rows],
  kofam_raw[star_rows, 6:ncol(kofam_raw), drop = FALSE]
)

kofam_fixed[!star_rows, 1:5] <- kofam_raw[!star_rows, 1:5]

kofam <- kofam_fixed %>%
  mutate(
    definition = apply(
      kofam_fixed[, 6:ncol(kofam_fixed), drop = FALSE],
      1,
      function(x) paste(na.omit(x[x != ""]), collapse = " ")
    )
  ) %>%
  select(1:5, definition)

colnames(kofam) <- c("gene", "KO", "threshold_flag", "score", "evalue", "definition")

kofam <- kofam %>%
  mutate(
    gene = as.character(gene),
    KO = as.character(KO),
    threshold_flag = as.character(threshold_flag),
    score = suppressWarnings(as.numeric(score)),
    evalue = suppressWarnings(as.numeric(evalue)),
    definition = as.character(definition)
  )

kofam_strict_all <- kofam %>%
  filter(
    !is.na(gene), gene != "",
    !is.na(KO), KO != "",
    threshold_flag == "*"
  ) %>%
  transmute(
    Protein = gene,
    KO = KO,
    Function_from_annotation = definition,
    Source = "KOfam",
    eggnog_evalue = NA_real_,
    kofamscan_evalue = evalue,
    kofamscan_score = score
  ) %>%
  distinct()

kofam_relaxed_input <- kofam_fixed %>%
  transmute(
    gene    = as.character(V2),
    KO      = as.character(V3),
    thrshld = suppressWarnings(as.numeric(V4)),
    score   = suppressWarnings(as.numeric(V5)),
    evalue  = suppressWarnings(as.numeric(V6))
  ) %>%
  filter(!is.na(gene), gene != "", !is.na(KO), KO != "")

kofam_definition_lookup <- kofam %>%
  transmute(
    gene = gene,
    KO = KO,
    definition = definition
  ) %>%
  distinct()

kofam_relaxed_all <- kofam_relaxed_input %>%
  left_join(kofam_definition_lookup, by = c("gene", "KO")) %>%
  filter(
    !is.na(evalue),
    evalue <= relaxed_evalue_cutoff,
    !is.na(score),
    !is.na(thrshld),
    score >= relaxed_score_fraction * thrshld
  ) %>%
  transmute(
    Protein = gene,
    KO = KO,
    Function_from_annotation = definition,
    Source = "KOfam_relaxed",
    eggnog_evalue = NA_real_,
    kofamscan_evalue = evalue,
    kofamscan_score = score
  ) %>%
  distinct()

protein_ko_membership <- bind_rows(
  kofam_strict_all,
  eggnog_all_kos,
  kofam_relaxed_all
) %>%
  mutate(
    Function_from_annotation = ifelse(
      is.na(Function_from_annotation) | Function_from_annotation == "",
      KO,
      Function_from_annotation
    )
  ) %>%
  distinct(Protein, KO, Source, .keep_all = TRUE)

write_csv(
  protein_ko_membership,
  file.path("Pathway_intermediate_tables", "Protein_KO_membership_all_sources.csv")
)

protein_ko_unique <- protein_ko_membership %>%
  distinct(Protein, KO, .keep_all = TRUE)

write_csv(
  protein_ko_unique,
  file.path("Pathway_intermediate_tables", "Protein_KO_unique_for_pathway_mapping.csv")
)

# ============================================================
# 4. KEGG REST: KO -> PATHWAY
# ============================================================

if (use_kegg_cache && file.exists(ko_pathway_cache_file)) {
  message("Loading cached KO -> pathway mapping...")
  ko_to_pathway <- read_csv(ko_pathway_cache_file, show_col_types = FALSE)
} else {
  message("Downloading KO -> pathway mapping from KEGG...")
  raw_link <- KEGGREST::keggLink("pathway", "ko")
  
  ko_to_pathway <- tibble(
    KO_raw = names(raw_link),
    pathway_raw = unname(raw_link)
  ) %>%
    transmute(
      KO = str_replace(KO_raw, "^ko:", ""),
      pathway_id = str_replace(pathway_raw, "^path:", "")
    ) %>%
    filter(str_detect(pathway_id, "^map")) %>%
    distinct()
  
  if (use_kegg_cache) {
    write_csv(ko_to_pathway, ko_pathway_cache_file)
  }
}

observed_kos <- unique(protein_ko_unique$KO)

ko_to_pathway_obs <- ko_to_pathway %>%
  filter(KO %in% observed_kos) %>%
  distinct()

if (nrow(ko_to_pathway_obs) == 0) {
  stop("No observed KOs could be linked to KEGG pathways.")
}

write_csv(
  ko_to_pathway_obs,
  file.path("Pathway_intermediate_tables", "Observed_KO_to_pathway_mapping.csv")
)

# ============================================================
# 5. KEGG REST: PATHWAY METADATA
# ============================================================

observed_pathways <- sort(unique(ko_to_pathway_obs$pathway_id))

if (use_kegg_cache && file.exists(pathway_meta_cache_file)) {
  message("Loading cached pathway metadata...")
  pathway_meta <- read_csv(pathway_meta_cache_file, show_col_types = FALSE)
  
  missing_paths <- setdiff(observed_pathways, pathway_meta$pathway_id)
  
  if (length(missing_paths) > 0) {
    message("Fetching metadata for newly needed pathways...")
    path_chunks <- chunk_vector(paste0("path:", missing_paths), chunk_size = 10)
    
    newly_fetched <- map_dfr(path_chunks, function(ids) {
      res <- safe_keggGet(ids, retries = 3, sleep_sec = 1)
      if (is.null(res)) return(tibble())
      
      map_dfr(res, function(entry) {
        nm <- if ("NAME" %in% names(entry)) entry$NAME[1] else NA_character_
        class_info <- parse_kegg_class(entry$CLASS)
        
        tibble(
          pathway_id = str_replace(entry$ENTRY[1], "\\s+.*$", ""),
          pathway_name = str_replace(nm, " - .*", ""),
          broad_category = class_info$broad_category,
          subcategory = class_info$subcategory
        )
      })
    })
    
    pathway_meta <- bind_rows(pathway_meta, newly_fetched) %>%
      distinct(pathway_id, .keep_all = TRUE)
    
    if (use_kegg_cache) {
      write_csv(pathway_meta, pathway_meta_cache_file)
    }
  }
} else {
  message("Downloading pathway metadata from KEGG...")
  path_chunks <- chunk_vector(paste0("path:", observed_pathways), chunk_size = 10)
  
  pathway_meta <- map_dfr(path_chunks, function(ids) {
    res <- safe_keggGet(ids, retries = 3, sleep_sec = 1)
    if (is.null(res)) return(tibble())
    
    map_dfr(res, function(entry) {
      nm <- if ("NAME" %in% names(entry)) entry$NAME[1] else NA_character_
      class_info <- parse_kegg_class(entry$CLASS)
      
      tibble(
        pathway_id = str_replace(entry$ENTRY[1], "\\s+.*$", ""),
        pathway_name = str_replace(nm, " - .*", ""),
        broad_category = class_info$broad_category,
        subcategory = class_info$subcategory
      )
    })
  }) %>%
    distinct(pathway_id, .keep_all = TRUE)
  
  if (use_kegg_cache) {
    write_csv(pathway_meta, pathway_meta_cache_file)
  }
}

pathway_meta_obs <- pathway_meta %>%
  filter(pathway_id %in% observed_pathways) %>%
  mutate(
    pathway_name = ifelse(is.na(pathway_name) | pathway_name == "", pathway_id, pathway_name),
    broad_category = ifelse(is.na(broad_category) | broad_category == "", "Unclassified", broad_category),
    subcategory = ifelse(is.na(subcategory) | subcategory == "", "Unclassified", subcategory)
  ) %>%
  distinct()

write_csv(
  pathway_meta_obs,
  file.path("Pathway_intermediate_tables", "Observed_pathway_metadata.csv")
)

# ============================================================
# 6. PROTEIN -> KO -> PATHWAY MEMBERSHIP
# ============================================================

protein_ko_pathway <- protein_ko_unique %>%
  inner_join(ko_to_pathway_obs, by = "KO") %>%
  left_join(pathway_meta_obs, by = "pathway_id") %>%
  distinct(Protein, KO, pathway_id, pathway_name, broad_category, subcategory, .keep_all = TRUE)

if (nrow(protein_ko_pathway) == 0) {
  stop("Protein -> KO -> pathway membership table is empty.")
}

write_csv(
  protein_ko_pathway,
  file.path("Pathway_intermediate_tables", "Protein_KO_pathway_membership.csv")
)

pathway_ko_summary <- protein_ko_pathway %>%
  distinct(pathway_id, pathway_name, broad_category, subcategory, KO) %>%
  group_by(pathway_id, pathway_name, broad_category, subcategory) %>%
  summarise(
    KOs_in_pathway = paste(sort(unique(KO)), collapse = ", "),
    n_KOs_in_pathway = n_distinct(KO),
    .groups = "drop"
  )

write_csv(
  pathway_ko_summary,
  file.path("Pathway_intermediate_tables", "Pathway_KO_summary.csv")
)
# ============================================================
# 6B. PROTEIN -> KO -> PATHWAY MEMBERSHIP
#     ONLY FOR PROTEINS PRESENT IN expr_imp / CROSSTAB FILES
# ============================================================

proteins_in_crosstab <- map_dfr(files, function(file_path) {
  df <- read_xlsx(file_path)
  
  if (!"Protein" %in% colnames(df)) {
    stop("File ", file_path, " does not contain a 'Protein' column.")
  }
  
  df %>%
    transmute(Protein = as.character(Protein)) %>%
    distinct()
}) %>%
  distinct()

protein_ko_pathway_crosstab_only <- protein_ko_pathway %>%
  semi_join(proteins_in_crosstab, by = "Protein") %>%
  distinct()

write_csv(
  protein_ko_pathway_crosstab_only,
  file.path(
    "Pathway_intermediate_tables",
    "Protein_KO_pathway_membership_crosstab_only.csv"
  )
)

message(
  "Protein_KO_pathway_membership_crosstab_only.csv written with ",
  nrow(protein_ko_pathway_crosstab_only),
  " protein-KO-pathway rows from proteins present in expr_imp/crosstab files."
)
# ============================================================
# 7. READ expr_imp FILES -> LONG SAMPLE TABLE
# ============================================================

all_expr_long <- map2_dfr(names(files), files, function(file_label, file_path) {
  message("Reading expression file: ", file_label)
  
  parsed <- parse_file_label(file_label)
  timepoint <- parsed$Timepoint
  sex <- parsed$Sex
  
  df <- read_xlsx(file_path)
  
  if (!"Protein" %in% colnames(df)) {
    stop("File ", file_path, " does not contain a 'Protein' column.")
  }
  
  df <- df %>%
    mutate(Protein = as.character(Protein))
  
  sample_cols <- get_sample_columns(df)
  
  if (length(sample_cols) < 2) {
    stop("File ", file_path, " does not contain expected sample columns like Control_1, Dosed_1, etc.")
  }
  
  df %>%
    select(Protein, all_of(sample_cols)) %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "Sample",
      values_to = "log2_intensity"
    ) %>%
    mutate(
      Group = get_group_from_sample(Sample),
      File = file_label,
      Sex = sex,
      Timepoint = timepoint,
      Sample_full = paste(file_label, Sample, sep = "__")
    ) %>%
    filter(!is.na(Group))
})

write_csv(
  all_expr_long,
  file.path("Pathway_intermediate_tables", "All_expression_long.csv")
)

# ============================================================
# 8. JOIN EXPRESSION WITH PATHWAY MEMBERSHIP
#    deduplicate Protein x Pathway x Sample
# ============================================================

expr_with_pathway <- all_expr_long %>%
  inner_join(protein_ko_pathway, by = "Protein") %>%
  distinct(
    Protein, pathway_id, pathway_name, broad_category, subcategory,
    File, Sex, Timepoint, Sample, Sample_full, Group,
    log2_intensity,
    .keep_all = TRUE
  )

if (nrow(expr_with_pathway) == 0) {
  stop("No proteins from the expression files could be mapped to KEGG pathways.")
}

write_csv(
  expr_with_pathway,
  file.path("Pathway_intermediate_tables", "Expression_with_pathway_membership.csv")
)

# ============================================================
# ============================================================
# 9. PATHWAY x SAMPLE SUMMARY
#    actual values used in heatmap come from here
# ============================================================

pathway_sample_long <- expr_with_pathway %>%
  group_by(
    pathway_id, pathway_name, broad_category, subcategory,
    File, Sex, Timepoint, Sample, Sample_full, Group
  ) %>%
  summarise(
    pathway_log2_intensity = median(log2_intensity, na.rm = TRUE),
    n_proteins_in_pathway = n_distinct(Protein),
    .groups = "drop"
  ) %>%
  filter(n_proteins_in_pathway >= min_proteins_per_pathway) %>%
  # EXCLUDE Human Diseases pathways
  filter(broad_category != "Human Diseases")

pathway_dataset_counts <- pathway_sample_long %>%
  distinct(pathway_id, File) %>%
  count(pathway_id, name = "n_datasets_present")

keep_pathways <- pathway_dataset_counts %>%
  filter(n_datasets_present >= min_datasets_per_pathway) %>%
  pull(pathway_id)

pathway_sample_long <- pathway_sample_long %>%
  filter(pathway_id %in% keep_pathways)

write_csv(
  pathway_sample_long,
  file.path("Pathway_intermediate_tables", "Pathway_sample_long_raw.csv")
)

pathway_info <- pathway_sample_long %>%
  distinct(pathway_id, pathway_name, broad_category, subcategory)

write_csv(
  pathway_info,
  file.path("Pathway_intermediate_tables", "Pathway_info_used_in_heatmaps.csv")
)

# ============================================================
# 10. WIDE MATRIX OF RAW PATHWAY VALUES
# ============================================================

pathway_sample_wide <- pathway_sample_long %>%
  select(pathway_id, pathway_name, broad_category, subcategory, Sample_full, pathway_log2_intensity) %>%
  pivot_wider(
    names_from = Sample_full,
    values_from = pathway_log2_intensity
  )

write_csv(
  pathway_sample_wide,
  file.path("Pathway_intermediate_tables", "Pathway_sample_matrix_raw.csv")
)

matrix_cols <- setdiff(
  colnames(pathway_sample_wide),
  c("pathway_id", "pathway_name", "broad_category", "subcategory")
)

mat_raw <- as.matrix(pathway_sample_wide[, matrix_cols, drop = FALSE])
rownames(mat_raw) <- pathway_sample_wide$pathway_id
mode(mat_raw) <- "numeric"

keep_rows <- rowSums(!is.na(mat_raw)) > 0
mat_raw <- mat_raw[keep_rows, , drop = FALSE]

pathway_info_mat <- pathway_sample_wide %>%
  filter(pathway_id %in% rownames(mat_raw)) %>%
  distinct(pathway_id, pathway_name, broad_category, subcategory) %>%
  slice(match(rownames(mat_raw), pathway_id))

if (use_global_zscore) {
  mat_global_z <- row_zscore_matrix(mat_raw)
  
  write_csv(
    as_tibble(mat_global_z, rownames = "pathway_id") %>%
      left_join(pathway_info_mat, by = "pathway_id") %>%
      relocate(pathway_id, pathway_name, broad_category, subcategory),
    file.path("Pathway_intermediate_tables", "Pathway_sample_matrix_global_zscore.csv")
  )
}

# ============================================================
# 11. SAMPLE ANNOTATION
# ============================================================

sample_annotation <- tibble(
  Sample_full = colnames(mat_raw)
) %>%
  separate(Sample_full, into = c("File", "Sample"), sep = "__", remove = FALSE) %>%
  mutate(
    Group = get_group_from_sample(Sample),
    Timepoint = str_extract(File, "^(5W|8W|12W)"),
    Sex = str_extract(File, "(M|F)$"),
    Sex = recode(Sex, "M" = "Male", "F" = "Female")
  )

write_csv(
  sample_annotation,
  file.path("Pathway_intermediate_tables", "Sample_annotation_all_columns.csv")
)

# ============================================================
# 12. FUNCTION TO PLOT ONE HEATMAP
# ============================================================

plot_one_heatmap <- function(target_file_label,
                             mat_raw,
                             pathway_info_mat,
                             sample_annotation,
                             use_global_zscore = FALSE,
                             max_pathways_per_heatmap = NULL) {
  
  cols_this <- sample_annotation %>%
    filter(File == target_file_label) %>%
    pull(Sample_full)
  
  if (length(cols_this) == 0) {
    message("No columns found for ", target_file_label)
    return(NULL)
  }
  
  mat_sub_raw <- mat_raw[, cols_this, drop = FALSE]
  
  keep_rows <- rowSums(!is.na(mat_sub_raw)) > 0
  mat_sub_raw <- mat_sub_raw[keep_rows, , drop = FALSE]
  
  info_sub <- pathway_info_mat %>%
    filter(pathway_id %in% rownames(mat_sub_raw)) %>%
    distinct() %>%
    slice(match(rownames(mat_sub_raw), pathway_id))
  
  col_anno_df <- sample_annotation %>%
    filter(Sample_full %in% cols_this) %>%
    slice(match(cols_this, Sample_full)) %>%
    select(Group) %>%
    as.data.frame()
  
  rownames(col_anno_df) <- cols_this
  group_vector <- col_anno_df$Group
  
  # -----------------------------
  # Row selection is based on Control vs Dosed difference
  # This does NOT change heatmap values.
  # -----------------------------
  keep_ids <- select_top_separating_pathways(
    mat_raw = mat_sub_raw,
    group_vector = group_vector,
    n_keep = max_pathways_per_heatmap
  )
  
  mat_sub_raw <- mat_sub_raw[keep_ids, , drop = FALSE]
  info_sub <- info_sub %>%
    filter(pathway_id %in% rownames(mat_sub_raw)) %>%
    slice(match(rownames(mat_sub_raw), pathway_id))
  
  # z-score for display values
  if (use_global_zscore) {
    mat_sub <- mat_global_z[rownames(mat_sub_raw), cols_this, drop = FALSE]
  } else {
    mat_sub <- row_zscore_matrix(mat_sub_raw)
  }
  
  # order rows inside category blocks by strongest group separation
  control_cols <- which(group_vector == "Control")
  dosed_cols   <- which(group_vector == "Dosed")
  
  effect_size <- apply(mat_sub_raw, 1, function(x) {
    control_mean <- mean(x[control_cols], na.rm = TRUE)
    dosed_mean   <- mean(x[dosed_cols], na.rm = TRUE)
    abs(dosed_mean - control_mean)
  })
  
  info_sub$effect_size <- effect_size[info_sub$pathway_id]
  
  ord <- order(
    info_sub$broad_category,
    -info_sub$effect_size,
    info_sub$pathway_name
  )
  
  mat_sub <- mat_sub[ord, , drop = FALSE]
  info_sub <- info_sub[ord, , drop = FALSE]
  
  row_labels <- info_sub$pathway_name
  
  # IMPORTANT:
  # Use a 2-level data frame so broad categories appear properly
  row_split_df <- data.frame(
    Broad_Category = info_sub$broad_category,
    stringsAsFactors = FALSE
  )
  
  legend_title <- "Z-score"
  
  top_ha <- HeatmapAnnotation(
    Group = col_anno_df$Group,
    
    col = list(
      Group = c(
        "Control" = "steelblue",
        "Dosed" = "firebrick"
      )
    ),
    
    annotation_name_gp = gpar(
      fontsize = 13,
      fontface = "bold"
    ),
    
    annotation_name_side = "right",
    
    annotation_legend_param = list(
      Group = list(
        title = "Group",
        
        title_gp = gpar(
          fontsize = 13,
          fontface = "bold"
        ),
        
        labels_gp = gpar(
          fontsize = 12,
          fontface = "bold"
        )
      )
    ),
    
    simple_anno_size = unit(0.42, "in"),
    gap = unit(3, "mm")
  )
  
  pretty_title <- case_when(
    target_file_label == "5W_M"  ~ "5W_Males Pathway Summary",
    target_file_label == "8W_M"  ~ "8W_Males Pathway Summary",
    target_file_label == "12W_M" ~ "12W_Males Pathway Summary",
    target_file_label == "5W_F"  ~ "5W_Females Pathway Summary",
    target_file_label == "8W_F"  ~ "8W_Females Pathway Summary",
    target_file_label == "12W_F" ~ "12W_Females Pathway Summary",
    TRUE ~ paste0(target_file_label, " Pathway Summary")
  )
  
  ht <- Heatmap(
    mat_sub,
    
    # Heatmap legend
    name = legend_title,
    
    # Heatmap colors
    col = colorRamp2(
      c(-2, 0, 2),
      c("#4575B4", "#FFFFBF", "#D73027")
    ),
    
    na_col = "grey90",
    
    # Sample-group annotation
    top_annotation = top_ha,
    
    # Do not cluster rows or columns
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    # ----------------------------------------------------------
    # ACTUAL KEGG PATHWAY NAMES
    # These are now larger and bold
    # ----------------------------------------------------------
    show_row_names = TRUE,
    row_names_side = "left",
    row_labels = row_labels,
    row_names_gp = gpar(
      fontsize = 12,
      fontface = "bold"
    ),
    row_names_max_width = unit(9, "in"),
    
    # Sample names remain hidden
    show_column_names = FALSE,
    
    # ----------------------------------------------------------
    # BROAD KEGG CATEGORY LABELS
    # These remain larger and bold
    # ----------------------------------------------------------
    row_split = row_split_df,
    row_title_side = "left",
    row_title_rot = 0,
    row_title_gp = gpar(
      fontsize = 12,
      fontface = "bold"
    ),
    
    # More separation between pathway-category blocks
    row_gap = unit(5, "mm"),
    
    # ----------------------------------------------------------
    # REMOVE FIGURE TITLE
    # ----------------------------------------------------------
    column_title = NULL,
    
    # ----------------------------------------------------------
    # BOLDER HEATMAP LEGEND / AXIS
    # ----------------------------------------------------------
    heatmap_legend_param = list(
      title = "Z-score",
      
      title_gp = gpar(
        fontsize = 13,
        fontface = "bold"
      ),
      
      labels_gp = gpar(
        fontsize = 12,
        fontface = "bold"
      ),
      
      at = c(-2, -1, 0, 1, 2),
      
      grid_width = unit(7, "mm"),
      grid_height = unit(7, "mm"),
      
      legend_height = unit(50, "mm")
    ),
    
    # Thin border around the heatmap improves definition
    border = TRUE,
    
    # Rasterization helps preserve quality in large heatmaps
    use_raster = FALSE
  )
  
  base_name <- if (use_global_zscore) {
    paste0(target_file_label, "_pathway_heatmap_top_control_vs_dosed_global_zscore")
  } else {
    paste0(target_file_label, "_pathway_heatmap_top_control_vs_dosed_within_dataset_zscore")
  }
  
  out_pdf  <- file.path("Pathway_heatmaps", paste0(base_name, ".pdf"))
  out_png  <- file.path("Pathway_heatmaps", paste0(base_name, ".png"))
  out_tiff <- file.path("Pathway_heatmaps", paste0(base_name, ".tiff"))
  
  # ------------------------------------------------------------
  # SAVE PDF
  # ------------------------------------------------------------
  pdf(
    out_pdf,
    width = pdf_width,
    height = pdf_height,
    onefile = FALSE,
    useDingbats = FALSE
  )
  
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(10, 25, 10, 10), "mm")
  )
  
  dev.off()
  
  
  # ------------------------------------------------------------
  # SAVE PNG
  # ------------------------------------------------------------
  png(
    out_png,
    width = png_width_px,
    height = png_height_px,
    res = png_res
  )
  
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(10, 25, 10, 10), "mm")
  )
  
  dev.off()
  
  
  # ------------------------------------------------------------
  # SAVE TIFF
  # ------------------------------------------------------------
  tiff(
    out_tiff,
    width = tiff_width_in,
    height = tiff_height_in,
    units = "in",
    res = tiff_res,
    compression = "lzw"
  )
  
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    padding = unit(c(10, 25, 10, 10), "mm")
  )
  
  dev.off()
  
  invisible(ht)
}

# ============================================================
# 13. PLOT ALL 6 HEATMAPS
# ============================================================

ordered_files <- c("5W_M", "8W_M", "12W_M", "5W_F", "8W_F", "12W_F")

heatmap_list <- map(
  ordered_files,
  ~plot_one_heatmap(
    target_file_label = .x,
    mat_raw = mat_raw,
    pathway_info_mat = pathway_info_mat,
    sample_annotation = sample_annotation,
    use_global_zscore = use_global_zscore,
    max_pathways_per_heatmap = max_pathways_per_heatmap
  )
)

names(heatmap_list) <- ordered_files

# ============================================================
# 14. EXTRA SUMMARY TABLES
# ============================================================

pathway_presence_summary <- pathway_sample_long %>%
  distinct(pathway_id, pathway_name, broad_category, subcategory, File, Sample_full, n_proteins_in_pathway) %>%
  group_by(pathway_id, pathway_name, broad_category, subcategory, File) %>%
  summarise(
    mean_n_proteins_across_samples = mean(n_proteins_in_pathway, na.rm = TRUE),
    max_n_proteins_across_samples = max(n_proteins_in_pathway, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  pathway_presence_summary,
  file.path("Pathway_intermediate_tables", "Pathway_presence_summary_by_file.csv")
)

pathway_sample_ko_summary <- expr_with_pathway %>%
  distinct(pathway_id, pathway_name, broad_category, subcategory, Protein, KO) %>%
  group_by(pathway_id, pathway_name, broad_category, subcategory) %>%
  summarise(
    n_proteins = n_distinct(Protein),
    n_KOs = n_distinct(KO),
    KO_list = paste(sort(unique(KO)), collapse = ", "),
    .groups = "drop"
  )

write_csv(
  pathway_sample_ko_summary,
  file.path("Pathway_intermediate_tables", "Pathway_summary_with_KOs_and_protein_counts.csv")
)

message("Done.")

library(openxlsx)

# ============================================================
# 15. SAVE SUPPLEMENTARY INFORMATION WORKBOOK
# ============================================================

dir.create("Supplementary_tables", showWarnings = FALSE)

si_workbook_file <- file.path(
  "Supplementary_tables",
  "Supplementary_Pathway_Heatmap_Tables.xlsx"
)

wb <- createWorkbook()

# -----------------------------
# Sheet 1: Analysis settings
# -----------------------------
analysis_settings <- tibble(
  Setting = c(
    "Z-score method",
    "Relaxed KOfam e-value cutoff",
    "Relaxed KOfam score fraction",
    "Minimum proteins per pathway",
    "Minimum datasets per pathway",
    "Maximum pathways per heatmap",
    "Human Diseases pathways excluded",
    "KEGG cache used"
  ),
  Value = c(
    ifelse(use_global_zscore, "Global row z-score across all datasets", "Row z-score within each dataset"),
    relaxed_evalue_cutoff,
    relaxed_score_fraction,
    min_proteins_per_pathway,
    min_datasets_per_pathway,
    ifelse(is.null(max_pathways_per_heatmap), "All pathways kept", max_pathways_per_heatmap),
    "Yes",
    use_kegg_cache
  )
)

addWorksheet(wb, "SI_Settings")
writeData(wb, "SI_Settings", analysis_settings)

# -----------------------------
# Main supplementary tables
# -----------------------------
si_tables <- list(
  "Protein_KO_all_sources" = protein_ko_membership,
  "Protein_KO_unique" = protein_ko_unique,
  "Observed_KO_to_pathway" = ko_to_pathway_obs,
  "Observed_pathway_metadata" = pathway_meta_obs,
  "Protein_KO_pathway" = protein_ko_pathway,
  "Protein_KO_pathway_expr_only" = protein_ko_pathway_crosstab_only,
  "Pathway_KO_summary" = pathway_ko_summary,
  "All_expression_long" = all_expr_long,
  "Expression_with_pathway" = expr_with_pathway,
  "Pathway_sample_long_raw" = pathway_sample_long,
  "Pathway_sample_matrix_raw" = pathway_sample_wide,
  "Sample_annotation" = sample_annotation,
  "Pathway_presence_by_file" = pathway_presence_summary,
  "Pathway_summary_KO_counts" = pathway_sample_ko_summary
)

for (sheet_name in names(si_tables)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, si_tables[[sheet_name]])
}

# -----------------------------
# Heatmap row selection summary
# -----------------------------
heatmap_row_selection_summary <- map_dfr(ordered_files, function(target_file_label) {
  
  cols_this <- sample_annotation %>%
    filter(File == target_file_label) %>%
    pull(Sample_full)
  
  mat_sub_raw <- mat_raw[, cols_this, drop = FALSE]
  keep_rows <- rowSums(!is.na(mat_sub_raw)) > 0
  mat_sub_raw <- mat_sub_raw[keep_rows, , drop = FALSE]
  
  col_anno_df <- sample_annotation %>%
    filter(Sample_full %in% cols_this) %>%
    slice(match(cols_this, Sample_full))
  
  group_vector <- col_anno_df$Group
  
  control_cols <- which(group_vector == "Control")
  dosed_cols   <- which(group_vector == "Dosed")
  
  row_stats <- tibble(
    pathway_id = rownames(mat_sub_raw),
    control_mean = apply(mat_sub_raw, 1, function(x) mean(x[control_cols], na.rm = TRUE)),
    dosed_mean = apply(mat_sub_raw, 1, function(x) mean(x[dosed_cols], na.rm = TRUE))
  ) %>%
    mutate(
      difference_dosed_minus_control = dosed_mean - control_mean,
      abs_difference = abs(difference_dosed_minus_control)
    ) %>%
    arrange(desc(abs_difference)) %>%
    mutate(
      Heatmap = target_file_label,
      Rank = row_number(),
      Included_in_heatmap = ifelse(
        is.null(max_pathways_per_heatmap),
        TRUE,
        Rank <= max_pathways_per_heatmap
      )
    ) %>%
    left_join(pathway_info_mat, by = "pathway_id") %>%
    select(
      Heatmap, Rank, Included_in_heatmap,
      pathway_id, pathway_name, broad_category, subcategory,
      control_mean, dosed_mean,
      difference_dosed_minus_control, abs_difference
    )
  
  row_stats
})

addWorksheet(wb, "Heatmap_row_selection")
writeData(wb, "Heatmap_row_selection", heatmap_row_selection_summary)

# -----------------------------
# Optional z-score matrix sheet
# -----------------------------
if (use_global_zscore) {
  pathway_sample_matrix_zscore <- as_tibble(mat_global_z, rownames = "pathway_id") %>%
    left_join(pathway_info_mat, by = "pathway_id") %>%
    relocate(pathway_id, pathway_name, broad_category, subcategory)
} else {
  pathway_sample_matrix_zscore <- as_tibble(row_zscore_matrix(mat_raw), rownames = "pathway_id") %>%
    left_join(pathway_info_mat, by = "pathway_id") %>%
    relocate(pathway_id, pathway_name, broad_category, subcategory)
}

addWorksheet(wb, "Pathway_matrix_zscore")
writeData(wb, "Pathway_matrix_zscore", pathway_sample_matrix_zscore)

# -----------------------------
# Save workbook
# -----------------------------
saveWorkbook(wb, si_workbook_file, overwrite = TRUE)

message("Supplementary workbook saved: ", si_workbook_file)