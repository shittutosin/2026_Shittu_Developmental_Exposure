


suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(tidyverse)
  library(RColorBrewer)
  library(writexl)
  library(nlme)
})

###############################################################################
# USER SETTINGS
###############################################################################

base_dir <- "base?dir"

timepoint_label <- "5weeks"
sex_keep <- "M"
tissue_label <- "Colon"

metadata_file <- "5weeks_metadata.txt"
read_summary_file <- "5weeks_read_summary.txt"


litter_analysis_root <- "/litter_ID"
out_dir <- file.path(litter_analysis_root, "16S_Analysis",
                     paste0(timepoint_label, "_", sex_keep, "_", tissue_label))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(out_dir)) stop("Cannot create output folder: ", out_dir)


litter_lookup_file <- file.path(litter_analysis_root,
                                "Supplemental Data 1-Litter ID.xlsx")
litter_lookup_sheet <- "Table S1"

permutation_target <- 9999L
exact_max_assignments <- 20000L
set.seed(123)

setwd(base_dir)

###############################################################################
# LOAD INPUT FILES
###############################################################################

files <- list.files(base_dir, pattern = "_rel-abundance.tsv$", full.names = TRUE)
samples <- gsub("_rel-abundance.tsv", "", basename(files))

meta <- read.table(metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
read_counts <- read.table(read_summary_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)


###############################################################################
# CREATE ANIMAL IDs TO MATCH ORIGINAL LITTER WORKBOOK
###############################################################################

week_code <- switch(
  tolower(timepoint_label),
  "5weeks"  = "5W",
  "8weeks"  = "8W",
  "12weeks" = "12W",
  "5w"      = "5W",
  "8w"      = "8W",
  "12w"     = "12W",
  stop("Unrecognized timepoint_label")
)

meta$AnimalID <- paste0(
  week_code, "_", toupper(trimws(meta$SampleID))
)

cat("\n16S Sample-to-Animal Matching:\n")
print(meta[, c("SampleID", "AnimalID", "Sex", "Treatment")])



if (!all(c("SampleID", "Sex", "Treatment", "AnimalID") %in% names(meta))) {
  stop("Metadata needs SampleID, Sex, Treatment, and the constructed AnimalID.")
}
meta$SampleID <- trimws(as.character(meta$SampleID))
meta$AnimalID <- trimws(as.character(meta$AnimalID))
if (anyNA(meta$SampleID) || any(meta$SampleID == "") ||
    anyDuplicated(meta$SampleID)) {
  stop("Missing or duplicated 16S SampleID values in metadata.")
}
if (!file.exists(litter_lookup_file)) {
  stop("Original litter workbook not found: ", litter_lookup_file)
}
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Install readxl: install.packages('readxl')")
}

# Identical ID normalization to the earlier developmental proteomics pipeline:
# remove underscores, hyphens, whitespace, and other punctuation, then uppercase.
normalize_animal_id <- function(x) {
  toupper(gsub("[^[:alnum:]]", "", trimws(as.character(x))))
}
normalize_sex <- function(x) {
  x <- tolower(trimws(as.character(x)))
  dplyr::case_when(x %in% c("m", "male", "males") ~ "M",
                   x %in% c("f", "female", "females") ~ "F",
                   TRUE ~ NA_character_)
}
normalize_treatment <- function(x) {
  x <- tolower(trimws(as.character(x)))
  dplyr::case_when(x %in% c("control", "c") ~ "Control",
                   x %in% c("dosed", "dose", "d") ~ "Dosed",
                   TRUE ~ NA_character_)
}
normalize_timepoint <- function(x) {
  x <- tolower(gsub("[^[:alnum:]]", "", as.character(x)))
  x <- sub("weeks$", "w", x)
  x <- sub("week$", "w", x)
  toupper(x)
}

raw_litter_lookup <- readxl::read_excel(
  litter_lookup_file,
  sheet = litter_lookup_sheet,
  .name_repair = "minimal"
)
required_lookup <- c("Animal ID", "Litter ID", "Sex", "Treatment", "Timepoint")
if (!all(required_lookup %in% names(raw_litter_lookup))) {
  stop("Litter workbook Table S1 is missing: ",
       paste(setdiff(required_lookup, names(raw_litter_lookup)), collapse = ", "))
}

litter_lookup <- raw_litter_lookup %>%
  transmute(
    SourceAnimalID = as.character(.data[["Animal ID"]]),
    JoinAnimalID = normalize_animal_id(.data[["Animal ID"]]),
    LitterID = trimws(as.character(.data[["Litter ID"]])),
    SourceSex = as.character(Sex),
    SourceTreatment = as.character(Treatment),
    SourceTimepoint = as.character(Timepoint)
  ) %>%
  filter(!is.na(JoinAnimalID), JoinAnimalID != "", JoinAnimalID != "NA")

if (anyDuplicated(litter_lookup$JoinAnimalID)) {
  duplicates <- unique(litter_lookup$JoinAnimalID[duplicated(litter_lookup$JoinAnimalID)])
  stop("Duplicated animal IDs in Table S1 after ID normalization: ",
       paste(duplicates, collapse = ", "))
}


old_litter_col <- names(meta)[
  tolower(gsub("[^a-z0-9]", "", tolower(names(meta)))) == "litterid"
]
if (length(old_litter_col) > 1L) {
  stop("More than one metadata LitterID-like column; review input metadata.")
}
if (length(old_litter_col) == 1L) {
  meta$MetadataLitterID <- trimws(as.character(meta[[old_litter_col]]))
  meta[[old_litter_col]] <- NULL
} else {
  meta$MetadataLitterID <- NA_character_
}

meta$JoinAnimalID <- normalize_animal_id(meta$AnimalID)
if (anyDuplicated(meta$JoinAnimalID)) {
  stop("Multiple 16S samples would match the same animal ID. Check metadata.")
}
meta <- meta %>%
  left_join(litter_lookup, by = "JoinAnimalID", relationship = "many-to-one")

# Audit the current run only, as metadata may include both sexes.
audit_current <- meta %>%
  filter(normalize_sex(Sex) == toupper(sex_keep),
         normalize_treatment(Treatment) %in% c("Control", "Dosed")) %>%
  mutate(
    MatchStatus = dplyr::case_when(
      is.na(SourceAnimalID) ~ "NO MATCH: animal not present in Table S1",
      is.na(LitterID) | LitterID == "" | tolower(LitterID) == "na" ~
        "MATCHED ANIMAL, BUT MISSING LITTER ID",
      is.na(normalize_sex(SourceSex)) ~ "MISSING/INVALID SOURCE SEX",
      is.na(normalize_treatment(SourceTreatment)) ~ "MISSING/INVALID SOURCE TREATMENT",
      is.na(SourceTimepoint) | trimws(SourceTimepoint) == "" ~
        "MISSING SOURCE TIMEPOINT",
      normalize_sex(Sex) != normalize_sex(SourceSex) ~ "SEX MISMATCH",
      normalize_treatment(Treatment) != normalize_treatment(SourceTreatment) ~
        "TREATMENT MISMATCH",
      normalize_timepoint(timepoint_label) != normalize_timepoint(SourceTimepoint) ~
        "TIMEPOINT MISMATCH",
      !is.na(MetadataLitterID) & MetadataLitterID != "" &
        tolower(MetadataLitterID) != "na" & MetadataLitterID != LitterID ~
        "METADATA LITTER CONFLICT",
      TRUE ~ "PASS"
    )
  ) %>%
  select(SampleID, AnimalID, JoinAnimalID, SourceAnimalID, LitterID,
         MetadataLitterID, Sex, Treatment, SourceSex,
         SourceTreatment, SourceTimepoint, MatchStatus)

# Export the audit BEFORE stopping, so ID problems are visible in Excel.
audit_file <- file.path(
  out_dir,
  paste0("16S_Litter_Match_Audit_", timepoint_label, "_", sex_keep, ".csv")
)
readr::write_csv(audit_current, audit_file, na = "")
cat("\n16S -> original Table S1 litter assignments:\n")
print(
  tibble::as_tibble(audit_current),
  n = Inf,
  width = Inf
)
cat("\nLitter-match audit written to: ", audit_file, "\n", sep = "")

if (!nrow(audit_current) || anyNA(audit_current$MatchStatus) ||
    any(audit_current$MatchStatus != "PASS")) {
  cat("\nOriginal Table S1 animal IDs for the current age/sex:\n")
  print(litter_lookup %>%
          filter(normalize_sex(SourceSex) == toupper(sex_keep),
                 normalize_timepoint(SourceTimepoint) ==
                   normalize_timepoint(timepoint_label)) %>%
          select(SourceAnimalID, LitterID, SourceTreatment), n = Inf)
  stop("Litter matching failed. Open the audit CSV and check MatchStatus. ",
       "No litter IDs were invented or substituted.")
}
cat("\nPASS: All selected mouse identities and Table S1 litter assignments agree.\n")



cat("\nFiles loaded:\n")
print(basename(files))

cat("\nMetadata samples:\n")
print(meta$SampleID)



long_df <- map2_dfr(files, samples, function(file, samp) {
  read_tsv(file, col_types = cols()) %>%
    dplyr::select(tax_id, lineage, abundance, `estimated counts`) %>%
    dplyr::rename(est_counts = `estimated counts`) %>%
    dplyr::mutate(sample = samp)
})

long_df <- long_df %>%
  filter(!tax_id %in% c("mapped_unclassified", "unmapped")) %>%
  mutate(
    abundance = replace_na(abundance, 0),
    est_counts = replace_na(est_counts, 0)
  ) %>%
  group_by(tax_id, sample) %>%
  summarise(
    abundance = abundance[1],
    est_counts = est_counts[1],
    lineage = lineage[1],
    .groups = "drop"
  )
clean_collapsed_reads <- long_df %>%
  left_join(
    meta %>% select(SampleID, Sex, Treatment, LitterID),
    by = c("sample" = "SampleID")
  ) %>%
  mutate(
    Timepoint = timepoint_label,
    Tissue = tissue_label
  ) %>%
  select(
    sample,
    Sex,
    Treatment,
    LitterID,
    tax_id,
    lineage,
    abundance,
    est_counts,
    Timepoint,
    Tissue
  )
###############################################################################
# COUNT MATRIX
###############################################################################

otu_counts <- long_df %>%
  dplyr::select(tax_id, sample, est_counts) %>%
  pivot_wider(
    names_from = sample,
    values_from = est_counts,
    values_fill = 0
  ) %>%
  arrange(tax_id)

count_mat <- as.matrix(otu_counts[, -1])
rownames(count_mat) <- otu_counts$tax_id

###############################################################################
# TAXONOMY TABLE
###############################################################################

ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

tax_mat <- long_df %>%
  dplyr::select(tax_id, lineage) %>%
  distinct() %>%
  mutate(lineage_split = strsplit(lineage, ";")) %>%
  rowwise() %>%
  mutate(
    lineage_split = list({
      pieces <- lineage_split[lineage_split != ""]
      pieces <- trimws(pieces)
      vals <- tail(pieces, 7)
      if (length(vals) < 7) vals <- c(rep(NA, 7 - length(vals)), vals)
      rev(vals)
    })
  ) %>%
  ungroup() %>%
  unnest_wider(lineage_split, names_sep = "_") %>%
  rename_with(~ ranks, starts_with("lineage_split_")) %>%
  mutate(across(all_of(ranks), trimws))

tax_mat <- tax_mat %>%
  mutate(
    Genus = na_if(Genus, ""),
    Species = na_if(Species, ""),
    Species = ifelse(
      is.na(Species) |
        Species == "" |
        Species == Genus |
        Species == word(Genus, 1) |
        Species %in% c("sp", "sp.", "bacterium", "uncultured bacterium"),
      paste(Genus, "sp."),
      Species
    )
  )

tax_mat <- tax_mat %>%
  group_by(Species) %>%
  mutate(
    Species = ifelse(
      n() == 1,
      Species,
      paste0(Species, "_tax", row_number())
    )
  ) %>%
  ungroup()

tax_mat <- tax_mat %>%
  column_to_rownames("tax_id")

tax_mat <- tax_mat[, ranks, drop = FALSE]

###############################################################################
# SUBSET METADATA
###############################################################################

keep_samples <- meta %>%
  filter(
    Sex == sex_keep,
    Treatment %in% c("Control", "Dosed"),
    SampleID %in% colnames(count_mat)
  ) %>% pull(SampleID)

if (!length(keep_samples)) {
  stop("No samples matched the selected sex and timepoint input files.")
}
count_mat_sub <- count_mat[, keep_samples, drop = FALSE]
meta_ps <- meta %>%
  filter(SampleID %in% keep_samples) %>%
  arrange(match(SampleID, keep_samples)) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Control", "Dosed")),
    LitterID = trimws(as.character(LitterID))
  )

if (any(is.na(meta_ps$LitterID) | meta_ps$LitterID == "" | 
        tolower(meta_ps$LitterID) == "na")) {
  missing_ids <- meta_ps$SampleID[
    is.na(meta_ps$LitterID) | meta_ps$LitterID == "" |
      tolower(meta_ps$LitterID) == "na"
  ]
  stop("Missing litter IDs. Correct the IDENTIFIER JOIN; do not infer litters: ",
       paste(missing_ids, collapse = ", "))
}
if (anyDuplicated(meta_ps$SampleID)) stop("Duplicate retained SampleID values.")
if (!all(colnames(count_mat_sub) == meta_ps$SampleID)) {
  stop("Count matrix / metadata order mismatch.")
}

litter_treatment_check <- meta_ps %>%
  distinct(LitterID, Treatment) %>% count(LitterID, name = "n_treatments")
if (any(litter_treatment_check$n_treatments != 1L)) {
  stop("A LitterID appears in both Control and Dosed. Check sample mapping.")
}

litter_design <- meta_ps %>%
  count(Treatment, LitterID, name = "n_mice") %>%
  arrange(Treatment, LitterID)
litter_group_counts <- litter_design %>%
  count(Treatment, name = "n_independent_litters")

cat("\nMouse and independent-litter distribution:\n")
print(litter_design)
print(litter_group_counts)
if (!all(c("Control", "Dosed") %in% as.character(litter_group_counts$Treatment))) {
  stop("Need BOTH Control and Dosed after filtering.")
}
if (any(litter_group_counts$n_independent_litters < 2L)) {
  warning("Fewer than TWO independent litters in a treatment group: ",
          "inferential litter-adjusted p-values will be NA.")
}

meta_ps <- meta_ps %>% column_to_rownames("SampleID")

cat("\nSamples used:\n")
print(keep_samples)
cat("\nMouse counts by treatment:\n")
print(table(meta_ps$Treatment))

###############################################################################
# BUILD PHYLOSEQ OBJECT
###############################################################################

OTU <- otu_table(count_mat_sub, taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax_mat[rownames(count_mat_sub), , drop = FALSE]))
META <- sample_data(meta_ps)

ps <- phyloseq(OTU, TAX, META)

###############################################################################
# RELATIVE ABUNDANCE OBJECT
###############################################################################

ps_rel <- transform_sample_counts(ps, function(x) {
  if (sum(x) == 0) x else x / sum(x)
})

###############################################################################
# RELATIVE ABUNDANCE LONG TABLE
###############################################################################

rel_long_taxid <- as.data.frame(otu_table(ps_rel))

if (!taxa_are_rows(ps_rel)) {
  rel_long_taxid <- t(rel_long_taxid)
}

rel_long_taxid <- rel_long_taxid %>%
  as.data.frame() %>%
  rownames_to_column("tax_id") %>%
  pivot_longer(
    cols = -tax_id,
    names_to = "SampleID",
    values_to = "RelAbund"
  ) %>%
  left_join(
    meta_ps %>% rownames_to_column("SampleID"),
    by = "SampleID"
  ) %>%
  left_join(
    tax_mat %>% rownames_to_column("tax_id"),
    by = "tax_id"
  )

###############################################################################
# WITHIN-TREATMENT TAXA FILTER FOR DA
###############################################################################

###############################################################################
# WITHIN-TREATMENT TAXA FILTER FOR DA
# Keep taxa if:
# Mean relative abundance is >= 0.1% in Control OR Dosed
###############################################################################

taxa_filter_by_treatment <- rel_long_taxid %>%
  group_by(tax_id, Species, Treatment) %>%
  summarise(
    mean_rel_abund = mean(RelAbund, na.rm = TRUE),
    mean_rel_abund_percent = mean_rel_abund * 100,
    n_samples_in_treatment = n(),
    n_present_in_treatment = sum(RelAbund > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    treatment_passes_0.1_percent_filter = mean_rel_abund >= 0.001
  )

taxa_filter_summary <- taxa_filter_by_treatment %>%
  group_by(tax_id, Species) %>%
  summarise(
    control_mean_rel_abund = mean_rel_abund[Treatment == "Control"],
    dosed_mean_rel_abund = mean_rel_abund[Treatment == "Dosed"],
    control_mean_percent = mean_rel_abund_percent[Treatment == "Control"],
    dosed_mean_percent = mean_rel_abund_percent[Treatment == "Dosed"],
    control_passes_filter = treatment_passes_0.1_percent_filter[Treatment == "Control"],
    dosed_passes_filter = treatment_passes_0.1_percent_filter[Treatment == "Dosed"],
    passed_filter = any(treatment_passes_0.1_percent_filter, na.rm = TRUE),
    .groups = "drop"
  )

taxa_keep <- taxa_filter_summary %>%
  filter(passed_filter) %>%
  pull(tax_id)

cat("\nNumber of taxa retained for DA using within-treatment abundance filter only:\n")
print(length(taxa_keep))

cat("\nL. rodentium filter status:\n")
print(
  taxa_filter_summary %>%
    filter(str_detect(Species, regex("rodentium", ignore_case = TRUE)))
)
###############################################################################
# BETA DIVERSITY -- ORIGINAL INDIVIDUAL MICE, WHOLE-LITTER RANDOMIZATION
###############################################################################


bray_dist <- phyloseq::distance(ps_rel, method = "bray")
ord_bray <- ordinate(ps_rel, method = "PCoA", distance = bray_dist)

bray_df <- as.data.frame(ord_bray$vectors[, 1:2, drop = FALSE]) %>%
  setNames(c("Axis.1", "Axis.2")) %>%
  rownames_to_column("SampleID") %>%
  left_join(meta_ps %>% rownames_to_column("SampleID"), by = "SampleID") %>%
  select(SampleID, Axis.1, Axis.2, Sex, Treatment, LitterID)

eig <- as.numeric(ord_bray$values$Eigenvalues)
pos_eig <- eig[is.finite(eig) & eig > 0]
if (length(pos_eig) < 2L) stop("PCoA has fewer than two positive axes.")
axis1_pct <- round(100 * eig[1] / sum(pos_eig), 1)
axis2_pct <- round(100 * eig[2] / sum(pos_eig), 1)

# PERMANOVA pseudo-F from Bray-Curtis squared distances, equivalent to
# vegan::adonis2(dist ~ Treatment) for this one-factor design.
D2 <- as.matrix(bray_dist)^2
stopifnot(identical(rownames(D2), rownames(meta_ps)))
permanova_F_R2 <- function(labels, D2) {
  labels <- as.character(labels)
  n <- length(labels)
  lev <- unique(labels)
  ng <- length(lev)
  if (ng < 2L || any(table(labels) < 1L)) return(c(F = NA_real_, R2 = NA_real_))
  ss_total <- sum(D2[upper.tri(D2)]) / n
  ss_within <- 0
  for (g in lev) {
    ii <- which(labels == g)
    if (length(ii) > 1L) {
      dg <- D2[ii, ii, drop = FALSE]
      ss_within <- ss_within + sum(dg[upper.tri(dg)]) / length(ii)
    }
  }
  ss_between <- ss_total - ss_within
  if (!is.finite(ss_within) || ss_within <= 0 || ss_total <= 0 || n <= ng) {
    return(c(F = NA_real_, R2 = NA_real_))
  }
  c(F = (ss_between / (ng - 1)) / (ss_within / (n - ng)),
    R2 = ss_between / ss_total)
}

# Generate treatment labels for ALL MICE, retaining each mouse's identity.
# Under the null, treatment labels are exchangeable between entire litters.
lit_ids <- unique(as.character(meta_ps$LitterID))
lit_treatment <- vapply(lit_ids, function(id) {
  as.character(unique(meta_ps$Treatment[meta_ps$LitterID == id]))
}, character(1))
n_control_litters <- sum(lit_treatment == "Control")
n_dosed_litters <- sum(lit_treatment == "Dosed")
n_total_litters <- length(lit_ids)

observed_stats <- permanova_F_R2(meta_ps$Treatment, D2)
permanova_method <- "Unestimated: insufficient independent litters"
perm_dosed <- matrix(FALSE, nrow = 0L, ncol = nrow(meta_ps))

if (min(n_control_litters, n_dosed_litters) >= 2L &&
    is.finite(observed_stats[["F"]])) {
  n_assignments <- choose(n_total_litters, n_control_litters)
  if (is.finite(n_assignments) && n_assignments <= exact_max_assignments) {
    allocation <- combn(seq_len(n_total_litters), n_control_litters)
    perm_dosed <- matrix(FALSE, nrow = ncol(allocation), ncol = nrow(meta_ps))
    for (i in seq_len(ncol(allocation))) {
      assigned_control <- lit_ids[allocation[, i]]
      perm_dosed[i, ] <- !(meta_ps$LitterID %in% assigned_control)
    }
    permanova_method <- "Exact whole-litter treatment-label enumeration"
  } else {
    perm_dosed <- matrix(FALSE, nrow = permutation_target, ncol = nrow(meta_ps))
    for (i in seq_len(permutation_target)) {
      assigned_control <- sample(lit_ids, n_control_litters, replace = FALSE)
      perm_dosed[i, ] <- !(meta_ps$LitterID %in% assigned_control)
    }
    permanova_method <- "Monte Carlo whole-litter treatment-label permutation"
  }
  if (nrow(perm_dosed) < 20L) {
    warning("Only ", nrow(perm_dosed),
            " distinct whole-litter treatment assignments are possible. ",
            "The exact permutation p-value has coarse resolution; p < 0.05 ",
            "may be impossible for this dataset.")
  }
  perm_F <- apply(perm_dosed, 1, function(is_dosed) {
    labels <- ifelse(is_dosed, "Dosed", "Control")
    permanova_F_R2(labels, D2)[["F"]]
  })
  if (grepl("Exact", permanova_method, fixed = TRUE)) {
    permanova_p <- mean(perm_F >= observed_stats[["F"]] - 1e-10, na.rm = TRUE)
  } else {
    permanova_p <- (1 + sum(perm_F >= observed_stats[["F"]] - 1e-10,
                            na.rm = TRUE)) / (length(perm_F) + 1)
  }
} else {
  permanova_p <- NA_real_
}

# Sanity check the observed F and R2 against vegan for this design.
vegan_check <- vegan::adonis2(bray_dist ~ Treatment,
                              data = as.data.frame(meta_ps), permutations = 0)
if (is.finite(observed_stats[["F"]]) &&
    (abs(vegan_check$F[1] - observed_stats[["F"]]) > 1e-6 ||
     abs(vegan_check$R2[1] - observed_stats[["R2"]]) > 1e-6)) {
  stop("Observed PERMANOVA F / R2 disagrees with vegan; inspect calculation.")
}

# Unrestricted PERMANOVA is included ONLY as a historical comparison.
set.seed(123)
original_permanova <- vegan::adonis2(bray_dist ~ Treatment,
                                     data = as.data.frame(meta_ps),
                                     permutations = 9999)

permanova_export <- tibble(
  Term = "Treatment",
  Sex = sex_keep,
  Timepoint = timepoint_label,
  Tissue = tissue_label,
  n_mice = nrow(meta_ps),
  n_control_litters = n_control_litters,
  n_dosed_litters = n_dosed_litters,
  F = as.numeric(observed_stats[["F"]]),
  R2 = as.numeric(observed_stats[["R2"]]),
  litter_adjusted_p = permanova_p,
  n_permutations_or_assignments = nrow(perm_dosed),
  permutation_method = permanova_method,
  original_unrestricted_p = as.numeric(original_permanova$`Pr(>F)`[1]),
  note = "Individual mouse Bray distances; labels permuted at whole-litter level"
)

label_text <- if (is.finite(permanova_p)) {
  paste0("p = ", signif(permanova_p, 3))
} else {
  "p = NA"
}

p_beta <- ggplot(bray_df, aes(Axis.1, Axis.2, color = Treatment)) +
  {if (all(table(bray_df$Treatment) >= 3))
    stat_ellipse(type = "t", linewidth = 1, level = 0.85)} +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(values = c("Control" = "steelblue", "Dosed" = "firebrick")) +
  annotate("text", x = Inf, y = Inf, label = label_text,
           hjust = 1.1, vjust = 1.5, size = 5, fontface = "bold") +
  labs(x = paste0("Axis 1 (", axis1_pct, "%)"),
       y = paste0("Axis 2 (", axis2_pct, "%)")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"),
        plot.title = element_blank())

ggsave(file.path(out_dir,
                 paste0("beta_diversity_pcoa_", timepoint_label, "_", sex_keep, ".tiff")),
       p_beta, device = "tiff", width = 5, height = 4, units = "in",
       dpi = 300, compression = "lzw")


count_mat_sub_round <- round(count_mat_sub)

ps_round <- phyloseq(
  otu_table(count_mat_sub_round, taxa_are_rows = TRUE),
  TAX,
  META
)

set.seed(123)

ps_rare <- rarefy_even_depth(
  ps_round,
  sample.size = 24000,
  rngseed = 123,
  replace = FALSE,
  trimOTUs = FALSE,
  verbose = FALSE
)

n_before <- nsamples(ps_round)
n_after <- nsamples(ps_rare)
dropped <- setdiff(sample_names(ps_round), sample_names(ps_rare))

cat("\nSamples before rarefaction:", n_before, "\n")
cat("Samples after rarefaction:", n_after, "\n")
cat("Dropped samples:", ifelse(length(dropped) == 0, "None", paste(dropped, collapse = ", ")), "\n")

alpha_df <- estimate_richness(
  ps_rare,
  measures = c("Observed", "Shannon", "Simpson")
) %>%
  rownames_to_column("SampleID") %>%
  left_join(
    meta_ps %>% rownames_to_column("SampleID"),
    by = "SampleID"
  )


fit_litter_continuous <- function(dat, response, transformation = "identity") {
  empty <- tibble(p_value = NA_real_, effect = NA_real_,
                  effect_se = NA_real_, model_method = NA_character_,
                  model_status = NA_character_)
  z <- dat %>%
    select(Treatment, LitterID, all_of(response)) %>%
    rename(OriginalValue = all_of(response)) %>%
    filter(!is.na(OriginalValue), is.finite(OriginalValue))
  z$Treatment <- factor(z$Treatment, levels = c("Control", "Dosed"))
  n_lit <- tapply(z$LitterID, z$Treatment, function(x) length(unique(x)))
  if (length(n_lit) < 2L || any(is.na(n_lit)) || any(n_lit < 2L)) {
    empty$model_status <- "Fewer than 2 litters per treatment"
    return(empty)
  }
  if (transformation == "asin_sqrt") {
    z$Response <- asin(sqrt(pmin(1, pmax(0, z$OriginalValue))))
  } else if (transformation == "log1p") {
    z$Response <- log1p(z$OriginalValue)
  } else {
    z$Response <- z$OriginalValue
  }
  if (length(unique(z$Response)) < 2L) {
    empty$model_status <- "Constant response"
    return(empty)
  }
  repeat_within_litter <- any(table(z$LitterID) > 1L)
  ans <- tryCatch({
    if (repeat_within_litter) {
      fit <- nlme::lme(Response ~ Treatment, random = ~1 | LitterID,
                       data = z, method = "REML",
                       control = nlme::lmeControl(msMaxIter = 200,
                                                  msMaxEval = 300,
                                                  returnObject = TRUE))
      tt <- summary(fit)$tTable
      est <- unname(tt["TreatmentDosed", "Value"])
      se <- unname(tt["TreatmentDosed", "Std.Error"])
      p <- unname(tt["TreatmentDosed", "p-value"])
      method <- "nlme::lme; random litter intercept"
    } else {
      fit <- lm(Response ~ Treatment, data = z)
      tt <- summary(fit)$coefficients
      est <- unname(tt["TreatmentDosed", "Estimate"])
      se <- unname(tt["TreatmentDosed", "Std. Error"])
      p <- unname(tt["TreatmentDosed", "Pr(>|t|)"])
      method <- "lm; one mouse per independent litter"
    }
    tibble(p_value = p, effect = est, effect_se = se,
           model_method = method, model_status = "OK")
  }, error = function(e) {
    empty$model_status <- paste0("Model failed: ", conditionMessage(e))
    empty
  })
  ans
}

pvals_alpha <- bind_rows(lapply(c("Observed", "Shannon", "Simpson"), function(m) {
  transform_m <- if (m == "Observed") "log1p" else "identity"
  fit_litter_continuous(alpha_df, m, transform_m) %>%
    mutate(Metric = m, model_scale = transform_m)
})) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"))

alpha_original <- tibble(Metric = c("Observed", "Shannon", "Simpson")) %>%
  mutate(original_wilcoxon_p = vapply(Metric, function(m) {
    tryCatch(wilcox.test(alpha_df[[m]] ~ alpha_df$Treatment,
                         exact = FALSE)$p.value,
             error = function(e) NA_real_)
  }, numeric(1)))
pvals_alpha <- left_join(pvals_alpha, alpha_original, by = "Metric")

alpha_long <- alpha_df %>%
  pivot_longer(cols = c("Observed", "Shannon", "Simpson"),
               names_to = "Metric", values_to = "Value") %>%
  left_join(pvals_alpha, by = "Metric")

p_label_df <- alpha_long %>%
  group_by(Metric) %>%
  summarise(p_value = first(p_value), .groups = "drop") %>%
  mutate(label = ifelse(is.finite(p_value),
                        paste0("p = ", signif(p_value, 3)),
                        "p = NA"))

p_alpha <- ggplot(alpha_long, aes(x = Treatment, y = Value, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2.5) +
  scale_fill_manual(values = c("Control" = "steelblue", "Dosed" = "firebrick")) +
  facet_wrap(~ Metric, scales = "free_y") +
  geom_text(data = p_label_df,
            aes(x = 1.5, y = Inf, label = label),
            vjust = 1.5, size = 4, fontface = "bold",
            inherit.aes = FALSE) +
  theme_bw(base_size = 12) +
  theme(legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 13, face = "bold"),
        strip.text = element_text(size = 15, face = "bold"),
        plot.title = element_blank()) +
  labs(x = NULL, y = "Diversity Value")

ggsave(file.path(out_dir,
                 paste0("alpha_diversity_", timepoint_label, "_", sex_keep, ".tiff")),
       p_alpha, device = "tiff", width = 8, height = 5,
       units = "in", dpi = 300, compression = "lzw")

alpha_stats_summary <- alpha_long %>%
  group_by(Metric, Treatment) %>%
  summarise(
    mean_value = mean(Value, na.rm = TRUE),
    sd_value = sd(Value, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Treatment,
    values_from = c(mean_value, sd_value, n_samples)
  ) %>%
  left_join(
    pvals_alpha,
    by = "Metric"
  ) %>%
  mutate(
    Sex = sex_keep,
    Timepoint = timepoint_label,
    Tissue = tissue_label
  )


compute_diff_rel <- function(df, sex, group1 = "Dosed", group2 = "Control") {
  df_sub <- df %>% filter(Sex == sex, Treatment %in% c(group1, group2))
  taxa <- df_sub %>% distinct(tax_id, Species)
  if (nrow(taxa) == 0L) {
    stop("No taxa passed the 0.1% within-treatment abundance filter.")
  }
  res <- purrr::map2_dfr(taxa$tax_id, taxa$Species, function(tx, sp) {
    dd <- df_sub %>% filter(tax_id == tx)
    a <- dd$RelAbund[dd$Treatment == group1]
    b <- dd$RelAbund[dd$Treatment == group2]
    original_p <- tryCatch(wilcox.test(a, b, exact = FALSE)$p.value,
                           error = function(e) NA_real_)
    fit <- fit_litter_continuous(dd, "RelAbund", "asin_sqrt")
    tibble(
      tax_id = tx, Species = sp,
      mean_group1 = mean(a, na.rm = TRUE),
      mean_group2 = mean(b, na.rm = TRUE),
      mean_group1_percent = 100 * mean_group1,
      mean_group2_percent = 100 * mean_group2,
      prevalence_group1 = sum(a > 0, na.rm = TRUE),
      prevalence_group2 = sum(b > 0, na.rm = TRUE),
      n_group1 = length(a), n_group2 = length(b),
      n_litters_group1 = length(unique(dd$LitterID[dd$Treatment == group1])),
      n_litters_group2 = length(unique(dd$LitterID[dd$Treatment == group2])),
      original_wilcoxon_p = original_p,
      pval = fit$p_value,
      model_effect_asin_sqrt = fit$effect,
      model_effect_se = fit$effect_se,
      model_method = fit$model_method,
      model_status = fit$model_status,
      log2FC = log2((mean_group1 + 1e-6) / (mean_group2 + 1e-6))
    )
  }) %>%
    mutate(
      original_wilcoxon_BH = p.adjust(original_wilcoxon_p, method = "BH"),
      padj = p.adjust(pval, method = "BH"),
      Sex = sex,
      Comparison = paste(group1, "vs", group2),
      direction = case_when(log2FC > 0 ~ group1,
                            log2FC < 0 ~ group2,
                            TRUE ~ "No change"),
      significant = is.finite(padj) & padj < 0.05
    ) %>% arrange(padj, pval)
  res
}

df_da <- rel_long_taxid %>% filter(tax_id %in% taxa_keep)
cat("\nLitter-aware differential-abundance tests for ",
    length(unique(df_da$tax_id)), " taxa...\n", sep = "")
res_da <- compute_diff_rel(df_da, sex = sex_keep,
                           group1 = "Dosed", group2 = "Control")
sig_da <- res_da %>% filter(significant)
cat("\nNumber significant after litter-aware BH FDR:", nrow(sig_da), "\n")
if (nrow(res_da)) {
  cat("Minimum litter-adjusted p:", suppressWarnings(min(res_da$pval, na.rm = TRUE)), "\n")
  cat("Minimum litter-adjusted q:", suppressWarnings(min(res_da$padj, na.rm = TRUE)), "\n")
}

###############################################################################
# PRESENCE / ABSENCE ANALYSIS FOR ALL TESTED TAXA
###############################################################################

presence_absence_df <- df_da %>%
  mutate(Present = RelAbund > 0)

presence_absence_stats <- presence_absence_df %>%
  group_by(tax_id, Species) %>%
  summarise(
    present_dosed = sum(Present[Treatment == "Dosed"], na.rm = TRUE),
    absent_dosed = sum(!Present[Treatment == "Dosed"], na.rm = TRUE),
    present_control = sum(Present[Treatment == "Control"], na.rm = TRUE),
    absent_control = sum(!Present[Treatment == "Control"], na.rm = TRUE),
    fisher_p = tryCatch(
      fisher.test(
        matrix(
          c(
            present_dosed, absent_dosed,
            present_control, absent_control
          ),
          nrow = 2,
          byrow = TRUE
        )
      )$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  filter(!is.na(fisher_p)) %>%
  mutate(
    fisher_p_adj = p.adjust(fisher_p, method = "BH"),
    prevalence_dosed_percent =
      present_dosed / (present_dosed + absent_dosed) * 100,
    prevalence_control_percent =
      present_control / (present_control + absent_control) * 100,
    prevalence_difference_percent =
      prevalence_dosed_percent - prevalence_control_percent,
    presence_direction = case_when(
      prevalence_difference_percent > 0 ~ "Dosed",
      prevalence_difference_percent < 0 ~ "Control",
      TRUE ~ "No difference"
    ),
    fisher_significant = fisher_p_adj < 0.05
  ) %>%
  arrange(fisher_p, fisher_p_adj)

presence_absence_ranked <- presence_absence_stats %>%
  arrange(
    desc(abs(prevalence_difference_percent)),
    fisher_p
  )

# Secondary inferential presence/absence test: permute WHOLE-LITTER
# treatment labels. Original 0/1 status is retained for every mouse.
# Original Fisher p is kept as an unadjusted-for-litter reference.
presence_litter_test <- function(dx) {
  dx <- dx[match(rownames(meta_ps), dx$SampleID), , drop = FALSE]
  if (nrow(perm_dosed) == 0L || anyNA(dx$Present)) return(NA_real_)
  y <- as.numeric(dx$Present)
  group <- as.character(meta_ps$Treatment)
  obs <- mean(y[group == "Dosed"]) - mean(y[group == "Control"])
  n_dose <- rowSums(perm_dosed)
  n_ctl <- length(y) - n_dose
  sums_d <- as.numeric(perm_dosed %*% y)
  diffs <- sums_d / n_dose - (sum(y) - sums_d) / n_ctl
  if (grepl("Exact", permanova_method, fixed = TRUE)) {
    mean(abs(diffs) >= abs(obs) - 1e-10)
  } else {
    (1 + sum(abs(diffs) >= abs(obs) - 1e-10)) / (length(diffs) + 1)
  }
}

presence_litter_p <- presence_absence_df %>%
  group_by(tax_id, Species) %>%
  group_modify(~ tibble(litter_permutation_p = presence_litter_test(.x))) %>%
  ungroup() %>%
  mutate(litter_permutation_q = p.adjust(litter_permutation_p, method = "BH"))

presence_absence_stats <- presence_absence_stats %>%
  left_join(presence_litter_p, by = c("tax_id", "Species")) %>%
  mutate(n_control_litters = n_control_litters,
         n_dosed_litters = n_dosed_litters,
         litter_aware_significant = is.finite(litter_permutation_q) &
           litter_permutation_q < 0.05)

presence_absence_ranked <- presence_absence_stats %>%
  arrange(litter_permutation_q, desc(abs(prevalence_difference_percent)))

###############################################################################
# DA PLOT — ONLY IF SIGNIFICANT TAXA EXIST
###############################################################################

if (nrow(sig_da) > 0) {
  
  p_da <- ggplot(
    sig_da,
    aes(
      x = log2FC,
      y = reorder(Species, log2FC),
      fill = factor(direction, levels = c("Control", "Dosed"))
    )
  ) +
    geom_col() +
    scale_fill_manual(
      values = c("Control" = "steelblue", "Dosed" = "firebrick"),
      labels = c("Control" = "Depleted", "Dosed" = "Enriched")
    ) +
    theme_bw(base_size = 14) +
    labs(
      title = paste0("Significant Differential Abundance: ", timepoint_label, " ", sex_keep),
      x = "Log2 Fold Change: Dosed vs Control",
      y = "Species",
      fill = "Direction"
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    )
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(
        "DA_significant_species_",
        timepoint_label,
        "_",
        sex_keep,
        ".tiff"
      )
    ),
    plot = p_da,
    device = "tiff",
    width = 10,
    height = 6,
    units = "in",
    dpi = 300,
    compression = "lzw"
  )
  
} else {
  cat("\nNo significant taxa after BH correction of the mouse-level model; no DA plot was saved.\n")
}


###############################################################################
# BARPLOTS
###############################################################################

ps_family <- tax_glom(ps_rel, taxrank = "Family", NArm = TRUE)
ps_genus <- tax_glom(ps_rel, taxrank = "Genus", NArm = TRUE)
ps_species <- tax_glom(ps_rel, taxrank = "Species", NArm = TRUE)

plot_bar_rank <- function(ps_obj, rank, topN = 20) {
  
  df <- psmelt(ps_obj)
  
  top_taxa <- df %>%
    group_by(.data[[rank]]) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    arrange(desc(Abundance)) %>%
    slice_head(n = topN) %>%
    pull(.data[[rank]])
  
  df <- df %>%
    mutate(Taxon = ifelse(.data[[rank]] %in% top_taxa, .data[[rank]], "Other"))
  
  taxa_levels <- unique(df$Taxon)
  
  base_cols <- c(
    brewer.pal(12, "Paired"),
    brewer.pal(8, "Dark2"),
    brewer.pal(8, "Set3")
  )
  
  base_cols <- base_cols[seq_along(taxa_levels)]
  names(base_cols) <- taxa_levels
  base_cols["Other"] <- "grey70"
  
  ggplot(df, aes(x = Sample, y = Abundance, fill = Taxon)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = base_cols) +
    theme_bw(base_size = 12) +
    theme(
      
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 12,
        face = "bold"
      ),
      
      axis.text.y = element_text(
        size = 12,
        face = "bold"
      ),
      
      axis.title = element_text(
        size = 16,
        face = "bold"
      ),
      
      legend.title = element_text(
        size = 15,
        face = "bold"
      ),
      
      legend.text = element_text(
        size = 13,
        face = "bold"
      ),
      
      plot.title = element_blank(),
      
      legend.position = "right"
      
    ) +
    labs(
      x = "Sample",
      y = "Relative Abundance"
    )
}

p_family <- plot_bar_rank(ps_family, "Family", topN = 20)
p_genus <- plot_bar_rank(ps_genus, "Genus", topN = 20)
p_species <- plot_bar_rank(ps_species, "Species", topN = 20)

ggsave(
  filename = file.path(
    out_dir,
    paste0(
      "family_barplot_",
      timepoint_label,
      "_",
      sex_keep,
      ".tiff"
    )
  ),
  plot = p_family,
  device = "tiff",
  width = 8,
  height = 5,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    out_dir,
    paste0(
      "genus_barplot_",
      timepoint_label,
      "_",
      sex_keep,
      ".tiff"
    )
  ),
  plot = p_genus,
  device = "tiff",
  width = 8,
  height = 5,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  filename = file.path(
    out_dir,
    paste0(
      "species_barplot_",
      timepoint_label,
      "_",
      sex_keep,
      ".tiff"
    )
  ),
  plot = p_species,
  device = "tiff",
  width = 8,
  height = 5,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

###############################################################################
# EXPORT TABLES
###############################################################################

# PERMANOVA export was computed above using whole-litter permutations.

beta_coordinates <- bray_df %>%
  mutate(
    Timepoint = timepoint_label,
    Tissue = tissue_label
  )

beta_axis_summary <- tibble(
  Sex = sex_keep,
  Timepoint = timepoint_label,
  Tissue = tissue_label,
  Axis = c("Axis.1", "Axis.2"),
  Variance_Explained_Percent = c(axis1_pct, axis2_pct)
)

alpha_export <- alpha_long %>%
  mutate(
    Sex = sex_keep,
    Timepoint = timepoint_label,
    Tissue = tissue_label
  ) %>%
  select(
    SampleID,
    Sex,
    Treatment,
    LitterID,
    Metric,
    Value,
    p_value,
    p_adj_BH,
    model_method,
    model_status,
    Timepoint,
    Tissue
  )

rarefaction_summary <- tibble(
  Sex = sex_keep,
  Timepoint = timepoint_label,
  Tissue = tissue_label,
  samples_before_rarefaction = n_before,
  samples_after_rarefaction = n_after,
  dropped_samples = ifelse(length(dropped) == 0, "None", paste(dropped, collapse = "; "))
)

###############################################################################
# CLEAN COLLAPSED READS — FILTER TO CURRENT SEX FOR THIS WORKBOOK
###############################################################################

clean_collapsed_reads_export <- clean_collapsed_reads %>%
  filter(Sex == sex_keep)

###############################################################################
# SPECIES-LEVEL RELATIVE ABUNDANCE
###############################################################################

relative_abundance_species <- rel_long_taxid %>%
  mutate(
    Relative_Abundance_Percent = RelAbund * 100,
    Timepoint = timepoint_label,
    Tissue = tissue_label
  ) %>%
  select(
    tax_id,
    SampleID,
    Sex,
    Treatment,
    LitterID,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species,
    Relative_Abundance = RelAbund,
    Relative_Abundance_Percent,
    Timepoint,
    Tissue
  )

###############################################################################
# PHYLUM / FAMILY / GENUS RELATIVE ABUNDANCE
###############################################################################

make_relative_abundance_rank <- function(relative_abundance_species, rank_col) {
  
  relative_abundance_species %>%
    filter(!is.na(.data[[rank_col]]), .data[[rank_col]] != "") %>%
    group_by(
      SampleID,
      Sex,
      Treatment,
      LitterID,
      Timepoint,
      Tissue,
      .data[[rank_col]]
    ) %>%
    summarise(
      Relative_Abundance = sum(Relative_Abundance, na.rm = TRUE),
      Relative_Abundance_Percent = sum(Relative_Abundance_Percent, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(Taxon = .data[[rank_col]]) %>%
    mutate(Taxonomic_Level = rank_col) %>%
    select(
      SampleID,
      Sex,
      Treatment,
      LitterID,
      Timepoint,
      Tissue,
      Taxonomic_Level,
      Taxon,
      Relative_Abundance,
      Relative_Abundance_Percent
    )
}

relative_abundance_phylum <- make_relative_abundance_rank(relative_abundance_species, "Phylum")
relative_abundance_family <- make_relative_abundance_rank(relative_abundance_species, "Family")
relative_abundance_genus  <- make_relative_abundance_rank(relative_abundance_species, "Genus")

###############################################################################
# DA EXPORT
###############################################################################

da_all_export <- res_da %>%
  mutate(
    Timepoint = timepoint_label,
    Tissue = tissue_label
  ) %>%
  rename(
    mean_dosed = mean_group1,
    mean_control = mean_group2,
    mean_dosed_percent = mean_group1_percent,
    mean_control_percent = mean_group2_percent,
    prevalence_dosed = prevalence_group1,
    prevalence_control = prevalence_group2,
    n_dosed = n_group1,
    n_control = n_group2,
    litter_adjusted_p = pval,
    litter_adjusted_q = padj
  )

###############################################################################
# GROUP MEAN ABUNDANCE — SPECIES LEVEL
###############################################################################

group_mean_abundance_species <- relative_abundance_species %>%
  group_by(
    tax_id,
    Timepoint,
    Tissue,
    Sex,
    Treatment,
    Kingdom,
    Phylum,
    Class,
    Order,
    Family,
    Genus,
    Species
  ) %>%
  summarise(
    mean_relative_abundance = mean(Relative_Abundance, na.rm = TRUE),
    sd_relative_abundance = sd(Relative_Abundance, na.rm = TRUE),
    mean_relative_abundance_percent = mean(Relative_Abundance_Percent, na.rm = TRUE),
    sd_relative_abundance_percent = sd(Relative_Abundance_Percent, na.rm = TRUE),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  left_join(
    da_all_export %>%
      select(
        tax_id,
        Species,
        DA_log2FC = log2FC,
        DA_litter_p = litter_adjusted_p,
        DA_litter_q = litter_adjusted_q,
        DA_direction = direction,
        DA_significant = significant
      ),
    by = c("tax_id", "Species")
  )

###############################################################################
# L. RODENTIUM ONLY REPORTS
###############################################################################

target_species_pattern <- "rodentium"
target_species_label <- "Limosilactobacillus rodentium"

l_rodentium_taxa <- rel_long_taxid %>%
  filter(str_detect(Species, regex(target_species_pattern, ignore_case = TRUE))) %>%
  distinct(tax_id, Species)

l_rodentium_sample_report <- rel_long_taxid %>%
  semi_join(l_rodentium_taxa, by = c("tax_id", "Species")) %>%
  mutate(
    Relative_Abundance_Percent = RelAbund * 100,
    Present = RelAbund > 0,
    Timepoint = timepoint_label,
    Tissue = tissue_label
  ) %>%
  select(
    tax_id,
    Species,
    SampleID,
    Sex,
    Treatment,
    LitterID,
    Relative_Abundance = RelAbund,
    Relative_Abundance_Percent,
    Present,
    Timepoint,
    Tissue
  )

l_rodentium_wilcoxon_report <- res_da %>%
  filter(str_detect(Species, regex(target_species_pattern, ignore_case = TRUE))) %>%
  mutate(
    Timepoint = timepoint_label,
    Tissue = tissue_label
  )

l_rodentium_presence_report <- presence_absence_stats %>%
  filter(str_detect(Species, regex(target_species_pattern, ignore_case = TRUE))) %>%
  mutate(
    Sex = sex_keep,
    Timepoint = timepoint_label,
    Tissue = tissue_label
  )

l_rodentium_combined_report <- l_rodentium_sample_report %>%
  distinct(tax_id, Species, Sex, Timepoint, Tissue) %>%
  left_join(
    l_rodentium_wilcoxon_report %>%
      select(
        tax_id,
        Species,
        Comparison,
        mean_dosed = mean_group1,
        mean_control = mean_group2,
        mean_dosed_percent = mean_group1_percent,
        mean_control_percent = mean_group2_percent,
        prevalence_dosed_abundance = prevalence_group1,
        prevalence_control_abundance = prevalence_group2,
        log2FC,
        litter_model_p = pval,
        litter_model_q = padj,
        original_wilcoxon_p,
        litter_model_significant = significant
      ),
    by = c("tax_id", "Species")
  ) %>%
  left_join(
    l_rodentium_presence_report %>%
      select(
        tax_id,
        Species,
        present_dosed,
        absent_dosed,
        present_control,
        absent_control,
        prevalence_dosed_percent,
        prevalence_control_percent,
        prevalence_difference_percent,
        fisher_p,
        fisher_p_adj,
        fisher_significant,
        litter_permutation_p,
        litter_permutation_q,
        litter_aware_significant
      ),
    by = c("tax_id", "Species")
  )

###############################################################################
# L. REUTERI STATS REPORT
###############################################################################

l_reuteri_taxa <- rel_long_taxid %>%
  filter(
    str_detect(
      Species,
      regex(
        "reuteri",
        ignore_case = TRUE
      )
    )
  ) %>%
  distinct(
    tax_id,
    Species
  )

l_reuteri_sample_report <- rel_long_taxid %>%
  semi_join(
    l_reuteri_taxa,
    by = c(
      "tax_id",
      "Species"
    )
  ) %>%
  mutate(
    Relative_Abundance_Percent = RelAbund * 100,
    Present = RelAbund > 0,
    Timepoint = timepoint_label,
    Tissue = tissue_label
  )

l_reuteri_wilcoxon_report <- res_da %>%
  filter(
    str_detect(
      Species,
      regex(
        "reuteri",
        ignore_case = TRUE
      )
    )
  ) %>%
  mutate(
    Timepoint = timepoint_label,
    Tissue = tissue_label
  )

l_reuteri_presence_report <- presence_absence_stats %>%
  filter(
    str_detect(
      Species,
      regex(
        "reuteri",
        ignore_case = TRUE
      )
    )
  ) %>%
  mutate(
    Sex = sex_keep,
    Timepoint = timepoint_label,
    Tissue = tissue_label
  )

l_reuteri_combined_report <- l_reuteri_sample_report %>%
  distinct(
    tax_id,
    Species,
    Sex,
    Timepoint,
    Tissue
  ) %>%
  left_join(
    l_reuteri_wilcoxon_report %>%
      select(
        tax_id,
        Species,
        Comparison,
        mean_dosed = mean_group1,
        mean_control = mean_group2,
        mean_dosed_percent = mean_group1_percent,
        mean_control_percent = mean_group2_percent,
        prevalence_dosed_abundance = prevalence_group1,
        prevalence_control_abundance = prevalence_group2,
        log2FC,
        litter_model_p = pval,
        litter_model_q = padj,
        original_wilcoxon_p,
        litter_model_significant = significant
      ),
    by = c(
      "tax_id",
      "Species"
    )
  ) %>%
  left_join(
    l_reuteri_presence_report %>%
      select(
        tax_id,
        Species,
        present_dosed,
        absent_dosed,
        present_control,
        absent_control,
        prevalence_dosed_percent,
        prevalence_control_percent,
        prevalence_difference_percent,
        fisher_p,
        fisher_p_adj,
        fisher_significant,
        litter_permutation_p,
        litter_permutation_q,
        litter_aware_significant
      ),
    by = c(
      "tax_id",
      "Species"
    )
  )
###############################################################################
# L. RODENTIUM PLOT EXTRACTION TABLES
###############################################################################

l_rodentium_plot_sample_values <- l_rodentium_sample_report %>%
  arrange(Timepoint, Tissue, Sex, Treatment, SampleID)

l_rodentium_plot_summary <- l_rodentium_sample_report %>%
  group_by(
    tax_id,
    Species,
    Timepoint,
    Tissue,
    Sex,
    Treatment
  ) %>%
  summarise(
    n_samples = n(),
    n_present = sum(Present, na.rm = TRUE),
    prevalence_percent = n_present / n_samples * 100,
    mean_relative_abundance = mean(Relative_Abundance, na.rm = TRUE),
    sd_relative_abundance = sd(Relative_Abundance, na.rm = TRUE),
    mean_relative_abundance_percent = mean(Relative_Abundance_Percent, na.rm = TRUE),
    sd_relative_abundance_percent = sd(Relative_Abundance_Percent, na.rm = TRUE),
    .groups = "drop"
  )

l_rodentium_plot_long <- l_rodentium_plot_summary %>%
  select(
    Species,
    Timepoint,
    Tissue,
    Sex,
    Treatment,
    mean_relative_abundance_percent,
    sd_relative_abundance_percent,
    prevalence_percent
  ) %>%
  pivot_longer(
    cols = c(mean_relative_abundance_percent, prevalence_percent),
    names_to = "Plot_Metric",
    values_to = "Plot_Value"
  )
###############################################################################
# FINAL SUPPLEMENTARY WORKBOOK
###############################################################################

supplementary_list <- list(

  "litter_design" = litter_design,
  "litter_group_counts" = litter_group_counts,
  "sample_litter_assignments" = meta_ps %>%
    rownames_to_column("SampleID") %>%
    select(SampleID, Sex, Treatment, LitterID),
  "alpha_litter_model_details" = pvals_alpha,
  "analysis_settings" = tibble(
    Setting = c("Timepoint", "Sex", "Tissue", "DA relative abundance transform",
                "Random effect", "PERMANOVA permutations",
                "Rarefaction depth", "DA abundance filter"),
    Value = c(timepoint_label, sex_keep, tissue_label, "asin(sqrt(p))",
              "LitterID", permanova_method,
              "24000", "Mean >=0.1% in either treatment")),
  "read_metrics" = read_counts,
  "clean_collapsed_reads" = clean_collapsed_reads_export,
  "relative_abundance_species" = relative_abundance_species,
  "relative_abundance_genus" = relative_abundance_genus,
  "relative_abundance_phylum" = relative_abundance_phylum,
  "relative_abundance_family" = relative_abundance_family,
  "group_mean_abundance_species" = group_mean_abundance_species,
  "alpha_sample_values" = alpha_export,
  "alpha_stats_summary" = alpha_stats_summary,
  "rarefaction_summary" = rarefaction_summary,
  "PERMANOVA" = permanova_export,
  "beta_coordinates" = beta_coordinates,
  "beta_axis_summary" = beta_axis_summary,
  "DA_all_relative_abundance" = da_all_export,
  "presence_absence_stats_all_taxa" = presence_absence_stats,
  "L_reuteri_stats_report" = l_reuteri_combined_report,
  "L_rodentium_sample_report" = l_rodentium_sample_report,
  "L_rodentium_stats_report" = l_rodentium_combined_report,
  "L_rodentium_plot_sample_values" = l_rodentium_plot_sample_values,
  "L_rodentium_plot_summary" = l_rodentium_plot_summary,
  "L_rodentium_plot_long" = l_rodentium_plot_long
)

write_xlsx(
  supplementary_list,
  path = file.path(
    out_dir,
    paste0("Supplementary_16S_outputs_", timepoint_label, "_", sex_keep, ".xlsx")
  )
)

###############################################################################
# SAVE EXTRA CSV FILES
###############################################################################

write.csv(
  da_all_export,
  file.path(out_dir, paste0("DA_all_relative_abundance_", timepoint_label, "_", sex_keep, ".csv")),
  row.names = FALSE
)

write.csv(
  presence_absence_stats,
  file.path(out_dir, paste0("presence_absence_stats_all_taxa_", timepoint_label, "_", sex_keep, ".csv")),
  row.names = FALSE
)

write.csv(
  l_rodentium_combined_report,
  file.path(out_dir, paste0("L_rodentium_stats_report_", timepoint_label, "_", sex_keep, ".csv")),
  row.names = FALSE
)

write.csv(
  l_reuteri_combined_report,
  file.path(
    out_dir,
    paste0(
      "L_reuteri_stats_report_",
      timepoint_label,
      "_",
      sex_keep,
      ".csv"
    )
  ),
  row.names = FALSE
)
write_xlsx(
  list(
    "sample_values" = l_rodentium_plot_sample_values,
    "plot_summary" = l_rodentium_plot_summary,
    "plot_long" = l_rodentium_plot_long,
    "stats_report" = l_rodentium_combined_report
  ),
  path = file.path(
    out_dir,
    paste0("L_rodentium_plot_extraction_", timepoint_label, "_", sex_keep, ".xlsx")
  )
)
###############################################################################
# DONE
###############################################################################

cat("\nDONE. Outputs saved to:\n")
cat(out_dir, "\n")

cat("\nFinal workbook tabs exported:\n")
print(names(supplementary_list))
cat("\nLitter analysis complete. Each mouse was kept separate.\n")
cat("Litter-adjusted PERMANOVA p:", permanova_p, "\n")
cat("Permutations/assignments:", nrow(perm_dosed), "\n")
cat("DA taxa fitted:", sum(is.finite(res_da$pval)), "/", nrow(res_da), "\n")
cat("Outputs:", normalizePath(out_dir, mustWork = FALSE), "\n")
