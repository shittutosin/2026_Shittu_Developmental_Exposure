
if (!requireNamespace("statmod", quietly = TRUE)) {
  stop("Package 'statmod' is required by limma::duplicateCorrelation().")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(ggplot2)
  library(limma)
  library(openxlsx)
})


# ============================================================
# 1. SETTINGS
# ============================================================

# UPDATED input root: the output from your successfully completed
# developmental pathway-heatmap + KO-mapping script.
base_dir <- "/base_dir"
source_dir <- file.path(base_dir, "Litter_Pathways_Original_Annotation")
input_dir <- file.path(source_dir, "Pathway_intermediate_tables")
if (!dir.exists(input_dir)) stop("Cannot find input directory: ", input_dir)

# Original filtered protein tables are used ONLY for peptide grouping.
# All abundance values remain from the existing normalized/imputed pathway input.
filtered_dir <- file.path(base_dir, "filtered_stuff")
if (!dir.exists(filtered_dir)) stop("Filtered TSV folder not found: ", filtered_dir)

timepoint_levels <- c("5W", "8W", "12W")
group_levels     <- c("Control", "Dosed")
sex_levels       <- c("Males", "Females")

# Require at least 3 biological samples in BOTH groups
# for a pathway to be tested.
min_bio_reps <- 3
min_litters_per_group <- 2

plot_width_grouped    <- 16
plot_height_grouped   <- 9

plot_width_individual  <- 10
plot_height_individual <- 8

plot_dpi <- 600


# ============================================================
# 2. OUTPUT DIRECTORIES
# ============================================================

out_dir <- file.path(base_dir, "Selected_pathway_ExactPeptides_ExactRawValues_Litter_Blocked_LIMMA")

dir.create(
  out_dir,
  showWarnings = FALSE
)

dir.create(
  file.path(out_dir, "tables"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(out_dir, "plots", "grouped"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(out_dir, "plots", "individual"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(out_dir, "supplementary"),
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. SELECTED PATHWAY GROUPS
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
  
  "Glutathione metabolism / redox" = c(
    "map00480"
  ),
  
  "LPS biosynthesis / endotoxin pathway" = c(
    "map00540"
  ),
  
  "Polyamine biosynthesis" = c(
    "map00330",
    "map00270",
    "map00220"
  ),
  
  "Glucuronate metabolism" = c(
    "map00040",
    "map00053"
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
  
  "Inflammatory signaling" = c(
    "map04010",
    "map04657",
    "map04668",
    "map04064",
    "map04620",
    "map04621",
    "map04060",
    "map04062",
    "map04630"
  ),
  
  "Oxidative stress / ROS" = c(
    "map00190",
    "map04146"
  ),
  
  "Lipid metabolism" = c(
    "map00071",
    "map00564",
    "map00061",
    "map00561"
  ),
  
  "Host glycan degradation (mucin-related)" = c(
    "map00511",
    "map00531",
    "map00520"
  )
)


# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================

clean_filename <- function(x) {
  
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}


p_to_label <- function(p) {
  
  case_when(
    is.na(p)    ~ NA_character_,
    p <= 0.0001 ~ "****",
    p <= 0.001  ~ "***",
    p <= 0.01   ~ "**",
    p <= 0.05   ~ "*",
    TRUE        ~ "ns"
  )
}


first_nonmissing <- function(x) {
  
  x <- as.character(x)
  
  x <- x[
    !is.na(x) &
      x != ""
  ]
  
  if (length(x) == 0) {
    return(NA_character_)
  }
  
  x[1]
}



expr_with_pathway <- read_csv(
  file.path(input_dir, "Expression_with_pathway_membership.csv"),
  show_col_types = FALSE
)
pathway_info <- read_csv(
  file.path(input_dir, "Observed_pathway_metadata.csv"),
  show_col_types = FALSE
)

# Current output columns: MouseID, Sample_full, LitterID, Dataset, Sex.
# Original script used the name Sample; preserve it for downstream tables.
expr_with_pathway <- expr_with_pathway %>%
  mutate(Sample = as.character(MouseID))


# ============================================================
# 6. REQUIRED COLUMN CHECK
# ============================================================

required_cols <- c(
  "Sample",
  "Sample_full",
  "MouseID",
  "LitterID",
  "Dataset",
  "Sex",
  "Group",
  "Timepoint",
  "Protein",
  "pathway_id",
  "log2_intensity"
)

missing_cols <- setdiff(
  required_cols,
  colnames(expr_with_pathway)
)

if (length(missing_cols) > 0) {
  
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}


# ============================================================
# 7. CANONICAL PATHWAY METADATA
#
# Important:
# exactly ONE metadata row per pathway_id.
# ============================================================

pathway_info_one <- pathway_info %>%
  group_by(pathway_id) %>%
  summarise(
    
    pathway_name =
      first_nonmissing(pathway_name),
    
    broad_category =
      first_nonmissing(broad_category),
    
    subcategory =
      first_nonmissing(subcategory),
    
    .groups = "drop"
  )


# ============================================================
# 8. MAKE SURE METADATA COLUMNS EXIST IN EXPRESSION TABLE
# ============================================================

if (!("pathway_name" %in% names(expr_with_pathway))) {
  expr_with_pathway$pathway_name <- NA_character_
}

if (!("broad_category" %in% names(expr_with_pathway))) {
  expr_with_pathway$broad_category <- NA_character_
}

if (!("subcategory" %in% names(expr_with_pathway))) {
  expr_with_pathway$subcategory <- NA_character_
}


# ============================================================
# 9. ADD CANONICAL PATHWAY METADATA
# ============================================================

expr_with_pathway <- expr_with_pathway %>%
  
  rename(
    original_pathway_name = pathway_name,
    original_broad_category = broad_category,
    original_subcategory = subcategory
  ) %>%
  
  left_join(
    pathway_info_one %>%
      rename(
        metadata_pathway_name = pathway_name,
        metadata_broad_category = broad_category,
        metadata_subcategory = subcategory
      ),
    by = "pathway_id"
  ) %>%
  
  mutate(
    
    pathway_name = coalesce(
      metadata_pathway_name,
      original_pathway_name,
      pathway_id
    ),
    
    broad_category = coalesce(
      metadata_broad_category,
      original_broad_category
    ),
    
    subcategory = coalesce(
      metadata_subcategory,
      original_subcategory
    )
  ) %>%
  
  select(
    -original_pathway_name,
    -original_broad_category,
    -original_subcategory,
    -metadata_pathway_name,
    -metadata_broad_category,
    -metadata_subcategory
  )


# ============================================================
# 10. CLEAN EXPERIMENTAL VARIABLES
# ============================================================

expr_with_pathway <- expr_with_pathway %>%
  mutate(
    
    Timepoint = factor(
      as.character(Timepoint),
      levels = timepoint_levels
    ),
    
    Group = factor(
      as.character(Group),
      levels = group_levels
    ),
    
    Sex = dplyr::recode(
      as.character(Sex),
      "M" = "Males",
      "F" = "Females",
      "Male" = "Males",
      "Female" = "Females"
    ),
    
    Sex = factor(
      Sex,
      levels = sex_levels
    )
  ) %>%
  filter(
    !is.na(Sample),
    !is.na(Sample_full),
    !is.na(LitterID),
    !is.na(Protein),
    !is.na(pathway_id),
    !is.na(log2_intensity),
    !is.na(Sex),
    !is.na(Group),
    !is.na(Timepoint)
  )


# ============================================================
# 11. USE THE NEW GLOBALLY UNIQUE MOUSE IDENTIFIER

expr_with_pathway <- expr_with_pathway %>%
  mutate(
    BioSampleID = as.character(Sample_full),
    LitterID = as.character(LitterID)
  )

# ============================================================
# 12. SAMPLE METADATA AND DAM TREATMENT CHECK
# ============================================================
sample_metadata_check <- expr_with_pathway %>%
  distinct(BioSampleID, Sample, Dataset, Sex, Timepoint, Group, LitterID)

if (anyDuplicated(sample_metadata_check$BioSampleID)) {
  stop("Each BioSampleID must have exactly one metadata assignment.")
}
if (anyNA(sample_metadata_check$LitterID) ||
    any(sample_metadata_check$LitterID == "")) {
  stop("Missing litter IDs; please check the sample mapping.")
}

# A dam must never appear in both treatment arms.
mixed_treatment_litters <- sample_metadata_check %>%
  group_by(LitterID) %>%
  summarise(n_treatment = n_distinct(Group), .groups = "drop") %>%
  filter(n_treatment != 1L)
if (nrow(mixed_treatment_litters)) {
  stop("Some litter IDs occur in both Control and Dosed.")
}

# Each sex x timepoint dataset should contain 4 mice/group, 2 litters/group.
litter_balance <- sample_metadata_check %>%
  group_by(Sex, Timepoint, Group) %>%
  summarise(n_mice = n_distinct(BioSampleID),
            n_litters = n_distinct(LitterID), .groups = "drop")
if (any(litter_balance$n_mice != 4L) ||
    any(litter_balance$n_litters != 2L)) {
  print(litter_balance)
  stop("Unexpected mice or litter count (expected 4 mice/2 litters in each arm).")
}

write_csv(sample_metadata_check,
          file.path(out_dir, "tables", "Biological_sample_metadata.csv"))
write_csv(litter_balance,
          file.path(out_dir, "tables", "Litter_balance_by_dataset.csv"))



canonical_peptides <- function(x) {
  if (length(x) != 1L || is.na(x) ||
      !nzchar(trimws(as.character(x))) ||
      toupper(trimws(as.character(x))) %in% c("-", "NA", "N/A", "0")) {
    return(NA_character_)
  }
  parts <- unlist(strsplit(trimws(as.character(x)), "[,;|[:space:]]+"),
                  use.names = FALSE)
  parts <- parts[nzchar(parts)]
  # Never guess how to parse modified, flanked, or count-only fields.
  if (!length(parts) || any(!grepl("^[A-Z]+$", parts))) {
    return(NA_character_)
  }
  paste(sort(unique(parts)), collapse = ";")
}

extract_taxon_filtered <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  vapply(x, function(s) {
    match_os <- stringr::str_match(
      s, stringr::regex("\\b(?:Tax|OS)\\s*=\\s*(.*?)(?=\\s+(?:TaxID|OX|RepID|GN|PE|SV)\\s*=|$)",
                        ignore_case = TRUE))
    if (!is.na(match_os[1, 2])) return(trimws(match_os[1, 2]))
    br <- stringr::str_extract_all(s, "\\[[^\\]]+\\]")[[1]]
    if (length(br)) return(gsub("^\\[|\\]$", "", tail(br, 1)))
    NA_character_
  }, character(1))
}

read_original_filtered <- function(dataset) {
  parts <- strsplit(dataset, "_", fixed = TRUE)[[1]]
  week <- sub("W$", "", parts[1])
  filename <- paste0(week, "wks_", parts[2], "_8_2_filtered_3.tsv")
  f <- file.path(filtered_dir, filename)
  if (!file.exists(f)) stop(dataset, ": missing original filtered TSV: ", f)
  message("Original peptide source: ", dataset, " <- ", filename)
  raw <- readr::read_tsv(
    f, col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE, progress = FALSE
  )
  required <- c("Protein", "Peptides")
  absent <- setdiff(required, names(raw))
  if (length(absent)) stop(filename, ": missing columns: ", paste(absent, collapse = ", "))
  raw$Protein <- trimws(as.character(raw$Protein))
  if (anyNA(raw$Protein) || any(!nzchar(raw$Protein)) || anyDuplicated(raw$Protein)) {
    stop(filename, ": duplicate, empty, or missing Protein IDs")
  }
  sample_cols <- grep("^[0-9]+weeks_[0-9]+_LG$", names(raw),
                      value = TRUE, ignore.case = TRUE)
  if (length(sample_cols) != 8L) {
    stop(filename, ": expected 8 original mouse intensity columns; found ",
         paste(sample_cols, collapse = ", "))
  }
  
  sample_cols <- sort(sample_cols)
  observed_matrix <- as.matrix(raw[, sample_cols, drop = FALSE])
  raw_profiles <- apply(observed_matrix, 1, function(z) {
    zz <- trimws(as.character(z))
    missing_token <- is.na(z) | zz %in% c("", "NA", "N/A")
    v <- suppressWarnings(as.numeric(z))
    if (any(!missing_token & !is.finite(v))) {
      stop(filename, ": invalid/non-finite observed intensity in original TSV")
    }
    if (any(missing_token & !is.na(v))) {
      stop(filename, ": inconsistent missing-value parsing in original TSV")
    }
    
    vals <- ifelse(missing_token, "<MISSING>",
                   ifelse(v == 0, "0", sprintf("%a", v)))
    paste0(sample_cols, "=", vals, collapse = "|")
  })
  # Preserve source annotation fields as evidence, NOT grouping criteria.
  uniq <- if ("Unique Peptide(s)" %in% names(raw))
    as.character(raw[["Unique Peptide(s)"]]) else rep(NA_character_, nrow(raw))
  cluster <- if ("Cluster" %in% names(raw))
    as.character(raw[["Cluster"]]) else rep(NA_character_, nrow(raw))
  data.frame(
    Dataset = dataset,
    Protein = raw$Protein,
    Peptides_Original = raw$Peptides,
    PeptideSet = vapply(raw$Peptides, canonical_peptides, character(1)),
    Unique_Peptide_Field = uniq,
    Source_Cluster = cluster,
    OrganismName = if ("Function" %in% names(raw))
      extract_taxon_filtered(raw$Function) else NA_character_,
    OriginalObservedProfile = raw_profiles,
    SourceFile = filename,
    stringsAsFactors = FALSE
  ) %>% as_tibble()
}

all_datasets <- sort(unique(as.character(expr_with_pathway$Dataset)))
filtered_protein_sources <- purrr::map_dfr(all_datasets, read_original_filtered)
readr::write_csv(filtered_protein_sources,
                 file.path(out_dir, "tables", "Original_filtered_peptide_source_audit.csv"))

# Read the ORIGINAL KO reference used by the established pathway workflow.
old_ko_reference <- file.path(input_dir,
                              "Original_KO_pathway_annotations_reference_ONLY.csv")
original_ko_reference <- file.path(
  base_dir, "..", "Pathway_Enrichment_specific", "high_confidence_hits",
  "KO_summary", "Pathway_intermediate_tables",
  "Protein_KO_pathway_membership.csv"
)
ko_candidates <- c(old_ko_reference, original_ko_reference)
ko_candidates <- ko_candidates[file.exists(ko_candidates)]
if (!length(ko_candidates)) {
  stop("Original KO/pathway annotation reference not found. Needed for safe KO-specific grouping.")
}
message("Original KO source: ", ko_candidates[[1L]])
ko_source <- readr::read_csv(ko_candidates[[1L]], show_col_types = FALSE)
if (!all(c("Protein", "pathway_id", "KO") %in% names(ko_source))) {
  stop("Original KO reference requires Protein, pathway_id and KO columns")
}
if (!"Function_from_annotation" %in% names(ko_source)) {
  ko_source$Function_from_annotation <- NA_character_
}
ko_source <- ko_source %>%
  mutate(
    Protein = trimws(as.character(Protein)),
    pathway_id = trimws(as.character(pathway_id)),
    KO = trimws(as.character(KO)),
    KO = if_else(grepl("^K[0-9]{5}$", KO), KO, NA_character_),
    Function_from_annotation = as.character(Function_from_annotation)
  ) %>%
  filter(!is.na(Protein), nzchar(Protein),
         !is.na(pathway_id), nzchar(pathway_id))

# If KO mapping has no Dataset column, expand by protein IDs present in each
# of the six original filtered files; never infer organism identity.
if ("Dataset" %in% names(ko_source)) {
  ko_ref <- ko_source %>%
    mutate(Dataset = trimws(as.character(Dataset))) %>%
    filter(Dataset %in% all_datasets) %>%
    semi_join(filtered_protein_sources %>% select(Dataset, Protein),
              by = c("Dataset", "Protein"))
} else {
  ko_ref <- filtered_protein_sources %>%
    select(Dataset, Protein) %>%
    inner_join(ko_source, by = "Protein", relationship = "many-to-many")
}
ko_ref <- ko_ref %>%
  distinct(Dataset, Protein, pathway_id, KO, Function_from_annotation)
if (!nrow(ko_ref)) stop("No KO annotations could be matched to filtered proteins")

# One KO per accession-pathway membership, or NA if missing/ambiguous.
member_key <- expr_with_pathway %>% distinct(Dataset, pathway_id, Protein)
missing_original_accessions <- member_key %>%
  anti_join(filtered_protein_sources %>% select(Dataset, Protein),
            by = c("Dataset", "Protein"))
if (nrow(missing_original_accessions) > 0L) {
  readr::write_csv(missing_original_accessions,
                   file.path(out_dir, "tables", "ERROR_missing_original_filtered_accessions.csv"))
  stop("Original filtered TSV files do not cover all pathway accession IDs (",
       nrow(missing_original_accessions), " missing). See ERROR file; refusing to group incomplete evidence.")
}
member_ko <- member_key %>%
  left_join(ko_ref %>% select(Dataset, pathway_id, Protein, KO),
            by = c("Dataset", "pathway_id", "Protein"),
            relationship = "many-to-many") %>%
  group_by(Dataset, pathway_id, Protein) %>%
  summarise(
    GroupKO = {
      labels <- unique(KO[!is.na(KO)])
      if (length(labels) == 1L) labels else NA_character_
    },
    KO_Assignment_Status = {
      k <- n_distinct(KO[!is.na(KO)])
      if (k == 0L) "No unambiguous KO"
      else if (k == 1L) "Single KO" else "Multiple KOs; no collapse"
    },
    .groups = "drop"
  )

# Build ONE common accession->KO->peptide-group lookup. This can also be
# reused by your separate protein-level KO plotting pipeline.
ko_group_lookup <- ko_ref %>%
  filter(!is.na(KO)) %>%
  distinct(Dataset, Protein, KO) %>%
  inner_join(filtered_protein_sources,
             by = c("Dataset", "Protein"), relationship = "many-to-one") %>%
  mutate(
    EvidenceStatus = if_else(is.na(PeptideSet),
                             "Unresolved peptide list: retained separately",
                             "Original peptides + original raw intensity profile parsed"),
    GroupingSignature = if_else(
      is.na(PeptideSet), paste0("UNRESOLVED:", Protein),
      paste0("EXACT_PEPTIDES:", PeptideSet,
             "||EXACT_ORIGINAL_8_MOUSE_VALUES:", OriginalObservedProfile)
    )
  ) %>%
  group_by(Dataset, KO, GroupingSignature) %>%
  mutate(
    Representative = sort(unique(Protein))[1],
    GroupSize = n_distinct(Protein),
    DistinctObservedProfiles = n_distinct(OriginalObservedProfile),
    GroupOrganism = if (all(!is.na(OrganismName)) &&
                        n_distinct(OrganismName) == 1L)
      first(OrganismName) else "Taxon unresolved"
  ) %>%
  ungroup() %>%
  group_by(Dataset) %>%
  mutate(ProteinGroupID = paste0(Dataset, "_PG_",
                                 sprintf("%06d", dense_rank(paste(KO, GroupingSignature))))) %>%
  ungroup()

if (anyDuplicated(ko_group_lookup[c("Dataset", "Protein", "KO")])) {
  stop("KO group mapping is not one-to-one per Dataset x Protein x KO")
}

# This lookup is safe to join to KO-level data by Dataset, Protein, KO.
readr::write_csv(
  ko_group_lookup %>% select(Dataset, Protein, KO, ProteinGroupID,
                             Representative, GroupSize, EvidenceStatus,
                             PeptideSet, Peptides_Original, Unique_Peptide_Field,
                             OriginalObservedProfile, GroupingSignature,
                             DistinctObservedProfiles,
                             OrganismName, GroupOrganism, SourceFile),
  file.path(out_dir, "tables", "KO_protein_group_lookup_FROM_FILTERED.csv")
)

peptide_matches_different_raw_values <- ko_group_lookup %>%
  filter(!is.na(PeptideSet)) %>%
  group_by(Dataset, KO, PeptideSet) %>%
  summarise(
    N_Accessions = n_distinct(Protein),
    N_Distinct_Original_Profiles = n_distinct(OriginalObservedProfile),
    N_Retained_Groups = n_distinct(ProteinGroupID),
    Accessions = paste(sort(unique(Protein)), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(N_Accessions > 1L, N_Distinct_Original_Profiles > 1L)
readr::write_csv(
  peptide_matches_different_raw_values,
  file.path(out_dir, "tables", "Same_peptides_DIFFERENT_original_values_NOT_collapsed.csv")
)


identical_evidence_different_KOs <- ko_group_lookup %>%
  filter(!is.na(PeptideSet)) %>%
  group_by(Dataset, PeptideSet, OriginalObservedProfile) %>%
  summarise(
    N_KOs = n_distinct(KO),
    KOs = paste(sort(unique(KO)), collapse = "; "),
    Accessions = paste(sort(unique(Protein)), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(N_KOs > 1L)
readr::write_csv(
  identical_evidence_different_KOs,
  file.path(out_dir, "tables", "Same_peptides_same_original_values_DIFFERENT_KO_review.csv")
)

# Safety check: EVERY multi-accession group must be exact on BOTH criteria.
invalid_groups <- ko_group_lookup %>%
  group_by(Dataset, KO, ProteinGroupID) %>%
  summarise(
    n_accessions = n_distinct(Protein),
    n_peptide_sets = n_distinct(PeptideSet),
    n_original_profiles = n_distinct(OriginalObservedProfile),
    .groups = "drop"
  ) %>%
  filter(n_accessions > 1L &
           (n_peptide_sets != 1L | n_original_profiles != 1L))
if (nrow(invalid_groups)) {
  readr::write_csv(invalid_groups,
                   file.path(out_dir, "tables", "ERROR_invalid_merged_groups.csv"))
  stop("Invalid collapse: accession group disagrees on original peptides or raw intensities")
}
message("STRICT ORIGINAL-EVIDENCE GROUP CHECK PASSED. ",
        "Same-peptide sets held separate due to different raw values: ",
        nrow(peptide_matches_different_raw_values),
        "; identical-evidence cross-KO conflicts for review: ",
        nrow(identical_evidence_different_KOs))

# Restrict grouping to the ORIGINAL pathway membership universe. No new
# pathway memberships are inferred or imported.
member_coverage <- expr_with_pathway %>%
  distinct(Dataset, pathway_id, Protein, BioSampleID) %>%
  group_by(Dataset, pathway_id, Protein) %>%
  summarise(PathwayMice = paste(sort(BioSampleID), collapse = "|"),
            .groups = "drop")

members <- member_ko %>%
  left_join(filtered_protein_sources,
            by = c("Dataset", "Protein"), relationship = "many-to-one") %>%
  left_join(
    ko_group_lookup %>% select(Dataset, Protein, KO,
                               ProteinGroupID, Representative, GroupSize,
                               DistinctObservedProfiles, GroupOrganism),
    by = c("Dataset", "Protein", "GroupKO" = "KO"),
    relationship = "many-to-one"
  ) %>%
  left_join(member_coverage,
            by = c("Dataset", "pathway_id", "Protein"),
            relationship = "one-to-one") %>%
  mutate(
    # The globally defined KO representative may not belong to EVERY
    # pathway. Preserve it for traceability but choose an accession that
    # is actually present in this specific pathway as the plot representative.
    GlobalRepresentative = coalesce(Representative, Protein),
    ProteinGroupID = coalesce(ProteinGroupID,
                              paste0(Dataset, "_UNRESOLVED_", Protein)),
    GroupSize = coalesce(GroupSize, 1L),
    GroupOrganism = coalesce(GroupOrganism, "Taxon unresolved")
  ) %>%
  group_by(Dataset, pathway_id, ProteinGroupID) %>%
  mutate(
    Representative = sort(unique(Protein))[1],
    PathwayGroupSize = n_distinct(Protein)
  ) %>%
  ungroup()

stopifnot(nrow(members) == nrow(member_key))
if (anyDuplicated(members[c("Dataset", "pathway_id", "Protein")])) {
  stop("Peptide group mapping unexpectedly duplicates pathway accession keys")
}
readr::write_csv(members,
                 file.path(out_dir, "tables", "Peptide_group_members_FROM_FILTERED.csv"))


raw_group_audit <- members %>%
  group_by(Dataset, pathway_id, ProteinGroupID, Representative, GroupKO) %>%
  summarise(
    N_Accessions_in_Pathway = n_distinct(Protein),
    N_Different_Observed_Profiles = n_distinct(OriginalObservedProfile),
    Original_Profiles_Differ = N_Different_Observed_Profiles > 1L,
    .groups = "drop"
  ) %>%
  filter(N_Accessions_in_Pathway > 1L)
readr::write_csv(raw_group_audit,
                 file.path(out_dir, "tables", "Multimember_group_raw_intensity_audit.csv"))
if (any(raw_group_audit$Original_Profiles_Differ)) {
  stop("Failed strict raw-profile invariant. See Multimember_group_raw_intensity_audit.csv")
}


n_before <- nrow(expr_with_pathway)
expr_with_pathway <- expr_with_pathway %>%
  left_join(
    members %>% select(Dataset, pathway_id, Protein,
                       Representative, ProteinGroupID, GroupKO),
    by = c("Dataset", "pathway_id", "Protein"),
    relationship = "many-to-one"
  ) %>%
  mutate(OriginalProtein = Protein, Protein = Representative) %>%
  select(-Representative)
if (nrow(expr_with_pathway) != n_before ||
    anyNA(expr_with_pathway$ProteinGroupID)) {
  stop("Protein group join failed: changed row count or missing group ID")
}

# Quantitative disagreement is an AUDIT, not a grouping criterion.
imputed_group_audit <- expr_with_pathway %>%
  group_by(Dataset, pathway_id, ProteinGroupID, BioSampleID) %>%
  summarise(
    N_Original_Accessions = n_distinct(OriginalProtein),
    Imputed_Range = max(log2_intensity) - min(log2_intensity),
    .groups = "drop"
  ) %>%
  filter(N_Original_Accessions > 1L, Imputed_Range > 1e-8)
readr::write_csv(imputed_group_audit,
                 file.path(out_dir, "tables", "Multimember_group_imputed_intensity_disagreement.csv"))

# Group-level median of EXISTING imputed accession values; NO new imputation.
protein_sample_all <- expr_with_pathway %>%
  group_by(BioSampleID, Sample, LitterID, Sex, Timepoint,
           Group, pathway_id, ProteinGroupID, Protein) %>%
  summarise(
    protein_log2_intensity = median(log2_intensity, na.rm = TRUE),
    original_row_count = n(),
    n_accessions_contributing = n_distinct(OriginalProtein),
    .groups = "drop"
  )

if (anyDuplicated(protein_sample_all[c("BioSampleID", "pathway_id", "ProteinGroupID")])) {
  stop("More than one value per biological sample x pathway x protein group")
}



pathway_sample_all <- protein_sample_all %>%
  
  group_by(
    BioSampleID,
    Sample,
    LitterID,
    Sex,
    Timepoint,
    Group,
    pathway_id
  ) %>%
  
  summarise(
    
    pathway_sample_median_log2 =
      median(
        protein_log2_intensity,
        na.rm = TRUE
      ),
    
    pathway_sample_mean_log2 =
      mean(
        protein_log2_intensity,
        na.rm = TRUE
      ),
    
    n_unique_proteins =
      n_distinct(ProteinGroupID),
    
    .groups = "drop"
  ) %>%
  
  left_join(
    pathway_info_one,
    by = "pathway_id"
  ) %>%
  
  mutate(
    pathway_name =
      coalesce(
        pathway_name,
        pathway_id
      )
  )




duplicate_check <- pathway_sample_all %>%
  
  count(
    BioSampleID,
    pathway_id
  ) %>%
  
  filter(n > 1)


if (nrow(duplicate_check) > 0) {
  
  print(
    duplicate_check,
    n = Inf
  )
  
  stop(
    "TRUE duplicate pathway values remain."
  )
  
} else {
  
  message(
    "PASS: one pathway value per biological sample x pathway."
  )
}



write_csv(
  pathway_sample_all,
  file.path(
    out_dir,
    "tables",
    "ALL_pathways_sample_level_scores.csv"
  )
)



run_pathway_limma <- function(pathway_df, target_sex, target_timepoint,
                              min_n = 3, min_litters = 2) {
  dat <- pathway_df %>%
    filter(Sex == target_sex, Timepoint == target_timepoint,
           Group %in% c("Control", "Dosed"))
  if (!nrow(dat)) return(tibble())
  
  # First obtain the mouse and independently-assigned litter counts.
  counts_mouse <- dat %>%
    group_by(pathway_id, Group) %>%
    summarise(n_bio = n_distinct(BioSampleID), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = n_bio, values_fill = 0,
                names_prefix = "n_")
  counts_litter <- dat %>%
    group_by(pathway_id, Group) %>%
    summarise(n_litters = n_distinct(LitterID), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = n_litters, values_fill = 0,
                names_prefix = "n_litters_")
  counts <- left_join(counts_mouse, counts_litter, by = "pathway_id")
  for (nm in c("n_Control", "n_Dosed",
               "n_litters_Control", "n_litters_Dosed")) {
    if (!nm %in% names(counts)) counts[[nm]] <- 0L
  }
  for (nm in c("n_Control", "n_Dosed",
               "n_litters_Control", "n_litters_Dosed")) {
    counts[[nm]][is.na(counts[[nm]])] <- 0L
  }
  
  group_summary <- dat %>%
    group_by(pathway_id, Group) %>%
    summarise(group_mean = mean(pathway_sample_median_log2, na.rm = TRUE),
              group_median = median(pathway_sample_median_log2, na.rm = TRUE),
              group_sd = sd(pathway_sample_median_log2, na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = Group,
                values_from = c(group_mean, group_median, group_sd))
  counts <- counts %>% left_join(group_summary, by = "pathway_id")
  
  names_tbl <- dat %>%
    distinct(pathway_id, pathway_name) %>%
    group_by(pathway_id) %>%
    summarise(pathway_name = first_nonmissing(pathway_name), .groups = "drop")
  
  testable_ids <- counts %>%
    filter(n_Control >= min_n, n_Dosed >= min_n,
           n_litters_Control >= min_litters,
           n_litters_Dosed >= min_litters) %>%
    pull(pathway_id)
  
  # Output all pathways even when insufficient biological replication.
  result <- tibble(
    pathway_id = character(), logFC = numeric(), AveExpr = numeric(),
    t = numeric(), p_value = numeric(), p_adj = numeric(), B = numeric(),
    Consensus_Litter_Correlation = numeric()
  )
  
  if (length(testable_ids)) {
    # Exactly ONE pathway score per mouse; never use individual proteins as
    # repeated observations in this model.
    pathway_matrix_df <- dat %>%
      filter(pathway_id %in% testable_ids) %>%
      select(pathway_id, BioSampleID, pathway_sample_median_log2) %>%
      distinct() %>%
      pivot_wider(names_from = BioSampleID,
                  values_from = pathway_sample_median_log2)
    pathway_matrix <- as.matrix(pathway_matrix_df[, -1, drop = FALSE])
    storage.mode(pathway_matrix) <- "double"
    rownames(pathway_matrix) <- pathway_matrix_df$pathway_id
    
    sample_metadata <- dat %>%
      distinct(BioSampleID, Sample, Group, LitterID)
    if (anyDuplicated(sample_metadata$BioSampleID)) {
      stop("Conflicting metadata for an individual mouse.")
    }
    sample_metadata <- sample_metadata[
      match(colnames(pathway_matrix), sample_metadata$BioSampleID), , drop = FALSE
    ]
    if (anyNA(sample_metadata$BioSampleID) || anyNA(sample_metadata$LitterID)) {
      stop("Cannot align expression columns with litter metadata.")
    }
    sample_metadata$Group <- factor(sample_metadata$Group,
                                    levels = c("Control", "Dosed"))
    if (anyNA(sample_metadata$Group)) stop("Invalid treatment group.")
    
    design <- model.matrix(~ 0 + Group, data = sample_metadata)
    colnames(design) <- sub("^Group", "", colnames(design))
    contrast_matrix <- limma::makeContrasts(
      Dosed_vs_Control = Dosed - Control, levels = design
    )
    block <- as.character(sample_metadata$LitterID)
    if (any(tapply(as.character(sample_metadata$Group), block,
                   function(x) length(unique(x)) != 1L))) {
      stop("One litter has mice in both treatment arms.")
    }
    if (!any(duplicated(block))) {
      stop("No repeated observations within litter. Cannot estimate correlation.")
    }
    corfit <- limma::duplicateCorrelation(pathway_matrix, design,
                                          block = block)
    rho <- as.numeric(corfit$consensus.correlation)
    if (length(rho) != 1L || !is.finite(rho)) {
      stop("No finite consensus litter correlation for ", target_sex,
           " at ", target_timepoint)
    }
    fit <- limma::lmFit(pathway_matrix, design, block = block,
                        correlation = rho)
    fit2 <- limma::contrasts.fit(fit, contrast_matrix)
    fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
    result <- limma::topTable(
      fit2, coef = "Dosed_vs_Control", number = Inf,
      adjust.method = "BH", sort.by = "none"
    ) %>%
      tibble::rownames_to_column("pathway_id") %>%
      transmute(pathway_id, logFC, AveExpr, t,
                p_value = P.Value, p_adj = adj.P.Val, B,
                Consensus_Litter_Correlation = rho)
  }
  
  counts %>%
    left_join(names_tbl, by = "pathway_id") %>%
    left_join(result, by = "pathway_id") %>%
    mutate(Sex = target_sex, Timepoint = target_timepoint,
           test_status = case_when(
             n_Control < min_n | n_Dosed < min_n ~
               "Insufficient biological replicates",
             n_litters_Control < min_litters |
               n_litters_Dosed < min_litters ~
               "Insufficient independent litters",
             is.na(p_value) ~ "Model not estimable",
             TRUE ~ "Tested"
           ))
}

# ============================================================
# 18. RUN LIMMA FOR ALL SIX SEX x TIMEPOINT COMPARISONS
# ============================================================

all_pathway_stats <- map_dfr(
  
  sex_levels,
  
  function(sx) {
    
    map_dfr(
      
      timepoint_levels,
      
      function(tp) {
        
        run_pathway_limma(
          
          pathway_df =
            pathway_sample_all,
          
          target_sex =
            sx,
          
          target_timepoint =
            tp,
          
          min_n =
            min_bio_reps,
          min_litters = min_litters_per_group
        )
      }
    )
  }
)


all_pathway_stats <- all_pathway_stats %>%
  
  mutate(
    
    Sex = factor(
      Sex,
      levels =
        sex_levels
    ),
    
    Timepoint = factor(
      Timepoint,
      levels =
        timepoint_levels
    ),
    
    signif_label = case_when(
      
      test_status !=
        "Tested" ~
        
        NA_character_,
      
      TRUE ~
        
        p_to_label(
          p_adj
        )
    ),
    
    direction = case_when(
      
      test_status !=
        "Tested" ~
        
        NA_character_,
      
      p_adj < 0.05 &
        logFC > 0 ~
        
        "Higher in Dosed",
      
      p_adj < 0.05 &
        logFC < 0 ~
        
        "Higher in Control",
      
      TRUE ~
        
        "Not significant"
    )
  )


# ============================================================
# 20. SAVE ALL-PATHWAY LIMMA RESULTS
# ============================================================

write_csv(
  all_pathway_stats,
  file.path(
    out_dir,
    "tables",
    "ALL_pathways_LIMMA_results.csv"
  )
)


# ============================================================
# 21. BUILD SELECTED-PATHWAY METADATA
# ============================================================

selected_pathway_table <- imap_dfr(
  
  selected_pathways,
  
  function(path_ids, group_name) {
    
    tibble(
      
      selected_group =
        group_name,
      
      pathway_id =
        path_ids
    )
  }
)


selected_pathway_metadata <- selected_pathway_table %>%
  
  left_join(
    pathway_info_one,
    by =
      "pathway_id"
  ) %>%
  
  mutate(
    
    pathway_name =
      coalesce(
        pathway_name,
        pathway_id
      ),
    
    pathway_label =
      pathway_name,
    
    selected_group = factor(
      selected_group,
      levels =
        names(
          selected_pathways
        )
    )
  )


# ============================================================
# 22. SELECTED-PATHWAY SCAFFOLD
#
# Ensures missing pathways/timepoints still appear.
# ============================================================

selected_scaffold <- selected_pathway_metadata %>%
  
  crossing(
    
    Sex = factor(
      sex_levels,
      levels =
        sex_levels
    ),
    
    Timepoint = factor(
      timepoint_levels,
      levels =
        timepoint_levels
    )
  )


# ============================================================
# 23. EXTRACT SAMPLE-LEVEL VALUES FOR SELECTED PATHWAYS
# ============================================================

selected_sample_values <- pathway_sample_all %>%
  
  inner_join(
    
    selected_pathway_metadata %>%
      
      select(
        selected_group,
        pathway_id,
        pathway_label
      ),
    
    by =
      "pathway_id"
  ) %>%
  
  mutate(
    
    selected_group = factor(
      selected_group,
      levels =
        names(
          selected_pathways
        )
    )
  )



selected_stats <- selected_scaffold %>%
  
  left_join(
    
    all_pathway_stats %>%
      
      select(
        -pathway_name
      ),
    
    by = c(
      "pathway_id",
      "Sex",
      "Timepoint"
    )
  ) %>%
  
  mutate(
    
    n_Control =
      replace_na(
        n_Control,
        0L
      ),
    
    n_Dosed =
      replace_na(
        n_Dosed,
        0L
      ),
    
    test_status = case_when(
      
      n_Control == 0 &
        n_Dosed == 0 ~
        
        "Not detected",
      
      is.na(
        test_status
      ) ~
        
        "Insufficient biological replicates",
      
      TRUE ~
        
        test_status
    ),
    
    signif_label = case_when(
      
      test_status ==
        "Tested" ~
        
        p_to_label(
          p_adj
        ),
      
      TRUE ~
        
        NA_character_
    )
  )


# ============================================================
# 25. Y POSITIONS FOR PLOT LABELS
#
# Based on SAMPLE-LEVEL pathway values.
# ============================================================

y_positions <- selected_sample_values %>%
  
  group_by(
    selected_group,
    Sex,
    pathway_id,
    pathway_label
  ) %>%
  
  summarise(
    
    y_max =
      max(
        pathway_sample_median_log2,
        na.rm = TRUE
      ),
    
    y_min =
      min(
        pathway_sample_median_log2,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  ) %>%
  
  mutate(
    
    y_range =
      y_max - y_min,
    
    y_pos =
      y_max +
      ifelse(
        y_range == 0,
        0.3,
        y_range * 0.15
      ),
    
    status_y =
      y_min +
      ifelse(
        y_range == 0,
        0.15,
        y_range * 0.50
      )
  )


selected_stats <- selected_stats %>%
  
  left_join(
    
    y_positions,
    
    by = c(
      "selected_group",
      "Sex",
      "pathway_id",
      "pathway_label"
    )
  ) %>%
  
  mutate(
    
    y_pos =
      ifelse(
        is.na(y_pos),
        1,
        y_pos
      ),
    
    status_y =
      ifelse(
        is.na(status_y),
        1,
        status_y
      )
  )


# ============================================================
# 26. SAVE SELECTED PATHWAY TABLES
# ============================================================

write_csv(
  selected_sample_values,
  file.path(
    out_dir,
    "tables",
    "Selected_pathways_sample_level_values.csv"
  )
)


write_csv(
  selected_stats,
  file.path(
    out_dir,
    "tables",
    "Selected_pathways_LIMMA_results.csv"
  )
)


write_csv(
  selected_pathway_metadata,
  file.path(
    out_dir,
    "tables",
    "Selected_pathway_membership.csv"
  )
)


# ============================================================
# 27. SAMPLE-LEVEL DESCRIPTIVE SUMMARY
# ============================================================

selected_pathway_summary <- selected_sample_values %>%
  
  group_by(
    selected_group,
    pathway_id,
    pathway_name,
    pathway_label,
    Sex,
    Timepoint,
    Group
  ) %>%
  
  summarise(
    
    n_biological_samples =
      n_distinct(
        BioSampleID
      ),
    
    mean_pathway_score =
      mean(
        pathway_sample_median_log2,
        na.rm = TRUE
      ),
    
    median_pathway_score =
      median(
        pathway_sample_median_log2,
        na.rm = TRUE
      ),
    
    sd_pathway_score =
      sd(
        pathway_sample_median_log2,
        na.rm = TRUE
      ),
    
    mean_unique_proteins_per_sample =
      mean(
        n_unique_proteins,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )


write_csv(
  selected_pathway_summary,
  file.path(
    out_dir,
    "tables",
    "Selected_pathways_descriptive_summary.csv"
  )
)



plot_group_sex <- function(
    group_name,
    target_sex,
    sample_df,
    stats_df) {
  
  
  df <- sample_df %>%
    
    filter(
      selected_group ==
        group_name,
      
      Sex ==
        target_sex
    ) %>%
    
    mutate(
      
      Timepoint = factor(
        Timepoint,
        levels =
          timepoint_levels
      ),
      
      Group = factor(
        Group,
        levels =
          group_levels
      )
    )
  
  
  stat_df <- stats_df %>%
    
    filter(
      selected_group ==
        group_name,
      
      Sex ==
        target_sex
    ) %>%
    
    mutate(
      
      Timepoint = factor(
        Timepoint,
        levels =
          timepoint_levels
      ),
      
      status_label = case_when(
        
        test_status ==
          "Not detected" ~
          
          "Not detected",
        
        test_status ==
          "Insufficient biological replicates" ~
          
          "Insufficient n",
        
        test_status ==
          "Insufficient independent litters" ~
          "Insufficient litters",
        
        test_status ==
          "Model not estimable" ~
          
          "Not estimable",
        
        TRUE ~
          
          NA_character_
      )
    )
  
  
  p <- ggplot() +
    
    geom_boxplot(
      
      data =
        df,
      
      aes(
        x =
          Timepoint,
        
        y =
          pathway_sample_median_log2,
        
        fill =
          Group
      ),
      
      position =
        position_dodge(
          width =
            0.75
        ),
      
      width =
        0.60,
      
      alpha =
        0.75,
      
      outlier.shape =
        NA
    ) +
    
    geom_jitter(
      
      data =
        df,
      
      aes(
        x =
          Timepoint,
        
        y =
          pathway_sample_median_log2,
        
        color =
          Group
      ),
      
      position =
        position_jitterdodge(
          jitter.width =
            0.08,
          
          dodge.width =
            0.75
        ),
      
      size =
        2.8,
      
      alpha =
        0.85
    ) +
    
    
    
    stat_summary(
      
      data =
        df,
      
      aes(
        x =
          Timepoint,
        
        y =
          pathway_sample_median_log2,
        
        group =
          Group
      ),
      
      fun =
        mean,
      
      geom =
        "point",
      
      position =
        position_dodge(
          width =
            0.75
        ),
      
      shape =
        23,
      
      size =
        3.6,
      
      fill =
        "white",
      
      color =
        "black"
    ) +
    
    geom_text(
      
      data =
        stat_df %>%
        
        filter(
          test_status ==
            "Tested"
        ),
      
      aes(
        x =
          Timepoint,
        
        y =
          y_pos,
        
        label =
          signif_label
      ),
      
      inherit.aes =
        FALSE,
      
      size =
        8,
      
      fontface =
        "bold"
    ) +
    
    geom_text(
      
      data =
        stat_df %>%
        
        filter(
          test_status !=
            "Tested"
        ),
      
      aes(
        x =
          Timepoint,
        
        y =
          status_y,
        
        label =
          status_label
      ),
      
      inherit.aes =
        FALSE,
      
      size =
        4.5,
      
      fontface =
        "bold",
      
      color =
        "grey30"
    ) +
    
    facet_wrap(
      
      ~ pathway_label,
      
      scales =
        "free_y",
      
      drop =
        FALSE
    ) +
    
    scale_fill_manual(
      
      values = c(
        "Control" =
          "steelblue",
        
        "Dosed" =
          "firebrick"
      )
    ) +
    
    scale_color_manual(
      
      values = c(
        "Control" =
          "steelblue",
        
        "Dosed" =
          "firebrick"
      )
    ) +
    
    scale_x_discrete(
      drop =
        FALSE
    ) +
    
    labs(
      
      title =
        NULL,
      
      x =
        "Timepoint",
      
      y =
        expression(
          Pathway~score~
            "(" *
            median~log[2]~
            protein~intensity *
            ")"
        ),
      
      fill =
        NULL,
      
      color =
        NULL
    ) +
    
    theme_bw(
      base_size =
        16
    ) +
    
    theme(
      
      plot.title =
        element_blank(),
      
      strip.background =
        element_rect(
          fill =
            "grey92",
          
          color =
            "black",
          
          linewidth =
            0.8
        ),
      
      strip.text =
        element_text(
          face =
            "bold",
          
          size =
            22,
          
          margin =
            margin(
              7,
              7,
              7,
              7
            )
        ),
      
      legend.position =
        "top",
      
      legend.title =
        element_blank(),
      
      legend.text =
        element_text(
          face =
            "bold",
          
          size =
            20
        ),
      
      axis.title.x =
        element_text(
          face =
            "bold",
          
          size =
            16
        ),
      
      axis.title.y =
        element_text(
          face =
            "bold",
          
          size =
            16
        ),
      
      axis.text.x =
        element_text(
          face =
            "bold",
          
          size =
            14,
          
          color =
            "black"
        ),
      
      axis.text.y =
        element_text(
          face =
            "bold",
          
          size =
            14,
          
          color =
            "black"
        ),
      
      panel.border =
        element_rect(
          color =
            "black",
          
          fill =
            NA,
          
          linewidth =
            0.8
        ),
      
      panel.grid.minor =
        element_blank(),
      
      panel.spacing =
        unit(
          1.1,
          "lines"
        )
    )
  
  
  file_stub <- paste0(
    
    clean_filename(
      group_name
    ),
    
    "_",
    
    clean_filename(
      target_sex
    )
  )
  
  
  ggsave(
    
    file.path(
      out_dir,
      "plots",
      "grouped",
      paste0(
        file_stub,
        ".pdf"
      )
    ),
    
    p,
    
    width =
      plot_width_grouped,
    
    height =
      plot_height_grouped,
    
    units =
      "in"
  )
  
  
  ggsave(
    
    file.path(
      out_dir,
      "plots",
      "grouped",
      paste0(
        file_stub,
        ".png"
      )
    ),
    
    p,
    
    width =
      plot_width_grouped,
    
    height =
      plot_height_grouped,
    
    units =
      "in",
    
    dpi =
      plot_dpi
  )
  
  
  ggsave(
    
    file.path(
      out_dir,
      "plots",
      "grouped",
      paste0(
        file_stub,
        ".tiff"
      )
    ),
    
    p,
    
    width =
      plot_width_grouped,
    
    height =
      plot_height_grouped,
    
    units =
      "in",
    
    dpi =
      plot_dpi,
    
    compression =
      "lzw"
  )
  
  
  p
}


# ============================================================
# 29. INDIVIDUAL PATHWAY PLOT
# ============================================================

plot_single_pathway_sex <- function(
    target_group,
    target_pathway_id,
    target_sex,
    sample_df,
    stats_df) {
  
  
  df <- sample_df %>%
    
    filter(
      selected_group ==
        target_group,
      
      pathway_id ==
        target_pathway_id,
      
      Sex ==
        target_sex
    ) %>%
    
    mutate(
      
      Timepoint = factor(
        Timepoint,
        levels =
          timepoint_levels
      ),
      
      Group = factor(
        Group,
        levels =
          group_levels
      )
    )
  
  
  stat_df <- stats_df %>%
    
    filter(
      selected_group ==
        target_group,
      
      pathway_id ==
        target_pathway_id,
      
      Sex ==
        target_sex
    ) %>%
    
    mutate(
      
      Timepoint = factor(
        Timepoint,
        levels =
          timepoint_levels
      ),
      
      status_label = case_when(
        
        test_status ==
          "Not detected" ~
          
          "Not detected",
        
        test_status ==
          "Insufficient biological replicates" ~
          
          "Insufficient n",
        
        test_status ==
          "Insufficient independent litters" ~
          "Insufficient litters",
        
        test_status ==
          "Model not estimable" ~
          
          "Not estimable",
        
        TRUE ~
          
          NA_character_
      )
    )
  
  
  pathway_nm <- selected_pathway_metadata %>%
    
    filter(
      selected_group ==
        target_group,
      
      pathway_id ==
        target_pathway_id
    ) %>%
    
    pull(
      pathway_name
    ) %>%
    
    first()
  
  
  p <- ggplot() +
    
    geom_boxplot(
      
      data =
        df,
      
      aes(
        x =
          Timepoint,
        
        y =
          pathway_sample_median_log2,
        
        fill =
          Group
      ),
      
      position =
        position_dodge(
          width =
            0.75
        ),
      
      width =
        0.60,
      
      alpha =
        0.75,
      
      outlier.shape =
        NA
    ) +
    
    geom_jitter(
      
      data =
        df,
      
      aes(
        x =
          Timepoint,
        
        y =
          pathway_sample_median_log2,
        
        color =
          Group
      ),
      
      position =
        position_jitterdodge(
          jitter.width =
            0.08,
          
          dodge.width =
            0.75
        ),
      
      size =
        3,
      
      alpha =
        0.85
    ) +
    
    stat_summary(
      
      data =
        df,
      
      aes(
        x =
          Timepoint,
        
        y =
          pathway_sample_median_log2,
        
        group =
          Group
      ),
      
      fun =
        mean,
      
      geom =
        "point",
      
      position =
        position_dodge(
          width =
            0.75
        ),
      
      shape =
        23,
      
      size =
        4,
      
      fill =
        "white",
      
      color =
        "black"
    ) +
    
    geom_text(
      
      data =
        stat_df %>%
        
        filter(
          test_status ==
            "Tested"
        ),
      
      aes(
        x =
          Timepoint,
        
        y =
          y_pos,
        
        label =
          signif_label
      ),
      
      inherit.aes =
        FALSE,
      
      size =
        9,
      
      fontface =
        "bold"
    ) +
    
    geom_text(
      
      data =
        stat_df %>%
        
        filter(
          test_status !=
            "Tested"
        ),
      
      aes(
        x =
          Timepoint,
        
        y =
          status_y,
        
        label =
          status_label
      ),
      
      inherit.aes =
        FALSE,
      
      size =
        5,
      
      fontface =
        "bold",
      
      color =
        "grey30"
    ) +
    
    facet_wrap(
      
      ~ pathway_label,
      
      scales =
        "free_y",
      
      drop =
        FALSE
    ) +
    
    scale_fill_manual(
      
      values = c(
        "Control" =
          "steelblue",
        
        "Dosed" =
          "firebrick"
      )
    ) +
    
    scale_color_manual(
      
      values = c(
        "Control" =
          "steelblue",
        
        "Dosed" =
          "firebrick"
      )
    ) +
    
    scale_x_discrete(
      drop =
        FALSE
    ) +
    
    labs(
      
      title =
        NULL,
      
      x =
        "Timepoint",
      
      y =
        expression(
          Pathway~score~
            "(" *
            median~log[2]~
            protein~intensity *
            ")"
        ),
      
      fill =
        NULL,
      
      color =
        NULL
    ) +
    
    theme_bw(
      base_size =
        16
    ) +
    
    theme(
      
      plot.title =
        element_blank(),
      
      strip.background =
        element_rect(
          fill =
            "grey92",
          
          color =
            "black",
          
          linewidth =
            0.8
        ),
      
      strip.text =
        element_text(
          face =
            "bold",
          
          size =
            24
        ),
      
      legend.position =
        "top",
      
      legend.text =
        element_text(
          face =
            "bold",
          
          size =
            20
        ),
      
      axis.title.x =
        element_text(
          face =
            "bold",
          
          size =
            20
        ),
      
      axis.title.y =
        element_text(
          face =
            "bold",
          
          size =
            20
        ),
      
      axis.text.x =
        element_text(
          face =
            "bold",
          
          size =
            18,
          
          color =
            "black"
        ),
      
      axis.text.y =
        element_text(
          face =
            "bold",
          
          size =
            18,
          
          color =
            "black"
        ),
      
      panel.border =
        element_rect(
          color =
            "black",
          
          fill =
            NA,
          
          linewidth =
            0.8
        ),
      
      panel.grid.minor =
        element_blank()
    )
  
  
  file_stub <- paste0(
    
    clean_filename(
      target_group
    ),
    
    "_",
    
    clean_filename(
      target_pathway_id
    ),
    
    "_",
    
    clean_filename(
      pathway_nm
    ),
    
    "_",
    
    clean_filename(
      target_sex
    )
  )
  
  
  ggsave(
    
    file.path(
      out_dir,
      "plots",
      "individual",
      paste0(
        file_stub,
        ".pdf"
      )
    ),
    
    p,
    
    width =
      plot_width_individual,
    
    height =
      plot_height_individual,
    
    units =
      "in"
  )
  
  
  ggsave(
    
    file.path(
      out_dir,
      "plots",
      "individual",
      paste0(
        file_stub,
        ".png"
      )
    ),
    
    p,
    
    width =
      plot_width_individual,
    
    height =
      plot_height_individual,
    
    units =
      "in",
    
    dpi =
      plot_dpi
  )
  
  
  ggsave(
    
    file.path(
      out_dir,
      "plots",
      "individual",
      paste0(
        file_stub,
        ".tiff"
      )
    ),
    
    p,
    
    width =
      plot_width_individual,
    
    height =
      plot_height_individual,
    
    units =
      "in",
    
    dpi =
      plot_dpi,
    
    compression =
      "lzw"
  )
  
  
  p
}


# ============================================================
# 30. GENERATE GROUPED PLOTS
# ============================================================

grouped_plot_list <- list()

for (grp in names(selected_pathways)) {
  
  for (sx in sex_levels) {
    
    plot_name <- paste(
      grp,
      sx,
      sep = " | "
    )
    
    grouped_plot_list[[plot_name]] <- plot_group_sex(
      
      group_name = grp,
      
      target_sex = sx,
      
      sample_df = selected_sample_values,
      
      stats_df = selected_stats
    )
  }
}

message("Grouped pathway plots generated successfully.")

# ============================================================
# 31. GENERATE INDIVIDUAL PATHWAY PLOTS
# ============================================================

individual_plot_list <- list()

for (grp in names(selected_pathways)) {
  
  for (pid in selected_pathways[[grp]]) {
    
    for (sx in sex_levels) {
      
      plot_name <- paste(
        grp,
        pid,
        sx,
        sep = " | "
      )
      
      individual_plot_list[[plot_name]] <- plot_single_pathway_sex(
        
        target_group = grp,
        
        target_pathway_id = pid,
        
        target_sex = sx,
        
        sample_df = selected_sample_values,
        
        stats_df = selected_stats
      )
    }
  }
}

message("Individual pathway plots generated successfully.")


ko_lookup <- members %>%
  semi_join(selected_pathway_metadata %>% distinct(pathway_id),
            by = "pathway_id") %>%
  left_join(
    ko_ref %>% select(Dataset, Protein, pathway_id, KO,
                      Function_from_annotation),
    by = c("Dataset", "Protein", "pathway_id"),
    relationship = "many-to-many"
  ) %>%
  group_by(Dataset, Protein = Representative, pathway_id) %>%
  summarise(
    KO = paste(sort(unique(KO[!is.na(KO) & nzchar(KO)])),
               collapse = "; "),
    Function_from_annotation = paste(sort(unique(
      Function_from_annotation[!is.na(Function_from_annotation) &
                                 nzchar(Function_from_annotation)])),
      collapse = "; "),
    .groups = "drop"
  )

expr_with_pathway <- expr_with_pathway %>%
  left_join(ko_lookup, by = c("Dataset", "Protein", "pathway_id"),
            relationship = "many-to-one")
readr::write_csv(ko_lookup,
                 file.path(out_dir, "tables", "Selected_pathway_protein_KO_annotation.csv"))

protein_annotation_cols <- intersect(
  
  c(
    "KO",
    "Function_from_annotation"
  ),
  
  colnames(
    expr_with_pathway
  )
)


protein_group_cols <- c(
  
  "BioSampleID",
  "Sample",
  "Sex",
  "Timepoint",
  "Group",
  "pathway_id",
  "ProteinGroupID",
  "Protein",
  
  protein_annotation_cols
)


protein_contribution_sample <- expr_with_pathway %>%
  
  semi_join(
    selected_pathway_metadata,
    by =
      "pathway_id"
  ) %>%
  
  group_by(
    across(
      all_of(
        protein_group_cols
      )
    )
  ) %>%
  
  summarise(
    
    protein_log2_intensity =
      median(
        log2_intensity,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  ) %>%
  
  left_join(
    
    selected_pathway_metadata %>%
      
      select(
        selected_group,
        pathway_id,
        pathway_name
      ),
    
    by =
      "pathway_id"
  )


protein_contribution_summary <- protein_contribution_sample %>%
  
  group_by(
    
    selected_group,
    pathway_id,
    pathway_name,
    Sex,
    Timepoint,
    Group,
    ProteinGroupID,
    Protein,
    
    across(
      all_of(
        protein_annotation_cols
      )
    )
  ) %>%
  
  summarise(
    
    n_biological_samples =
      n_distinct(
        BioSampleID
      ),
    
    mean_protein_log2 =
      mean(
        protein_log2_intensity,
        na.rm = TRUE
      ),
    
    median_protein_log2 =
      median(
        protein_log2_intensity,
        na.rm = TRUE
      ),
    
    sd_protein_log2 =
      sd(
        protein_log2_intensity,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  )


write_csv(
  protein_contribution_summary,
  file.path(
    out_dir,
    "tables",
    "Selected_pathway_protein_contributions.csv"
  )
)


# ============================================================
# 33. DETECTION SUMMARY
# ============================================================

detection_summary <- selected_scaffold %>%
  
  crossing(
    
    Group = factor(
      group_levels,
      levels =
        group_levels
    )
  ) %>%
  
  left_join(
    
    selected_sample_values %>%
      
      group_by(
        selected_group,
        pathway_id,
        Sex,
        Timepoint,
        Group
      ) %>%
      
      summarise(
        
        n_samples_detected =
          n_distinct(
            BioSampleID
          ),
        
        mean_unique_proteins =
          mean(
            n_unique_proteins,
            na.rm = TRUE
          ),
        
        .groups =
          "drop"
      ),
    
    by = c(
      "selected_group",
      "pathway_id",
      "Sex",
      "Timepoint",
      "Group"
    )
  ) %>%
  
  mutate(
    
    n_samples_detected =
      replace_na(
        n_samples_detected,
        0L
      )
  )


write_csv(
  detection_summary,
  file.path(
    out_dir,
    "tables",
    "Selected_pathways_detection_summary.csv"
  )
)


# ============================================================
# 34. ANALYSIS SETTINGS TABLE
# ============================================================

analysis_settings <- tibble(
  
  Setting = c(
    
    "Input file",
    
    "Experimental unit",
    
    "Protein deduplication",
    
    "Pathway value per biological sample",
    
    "Statistical method",
    
    "Contrast",
    
    "Minimum biological replicates",
    "Minimum independent litters per treatment",
    "Litter blocking",
    "Limited independent treatment replication",
    "Multiple-testing family",
    
    "P-value adjustment",
    
    "Positive logFC",
    
    "Negative logFC",
    
    "Plot dot meaning"
  ),
  
  Value = c(
    
    file.path(input_dir, "Expression_with_pathway_membership.csv"),
    
    "One score/mouse/pathway; dams define treatment-level clusters",
    
    "Exact original peptide-set grouping within cohort and unambiguous KO (independent of imputation); median of existing imputed accession values per mouse/group; ambiguities retained separately",
    
    "Median log2 intensity across peptide-defined protein groups within each biological sample",
    
    "limma duplicateCorrelation + blocked lmFit + robust empirical Bayes",
    
    "Dosed - Control",
    
    paste0(min_bio_reps, " mice per group"),
    paste0(min_litters_per_group, " litters per group"),
    "LitterID (within-litter correlation estimated by duplicateCorrelation)",
    "2 independent dams per treatment within each sex x timepoint: interpret FDR cautiously",
    
    "All testable KEGG pathways within each Sex x Timepoint comparison",
    
    "Benjamini-Hochberg FDR",
    
    "Higher pathway score in Dosed",
    
    "Higher pathway score in Control",
    
    "One dot = one biological sample"
  )
)


# ============================================================
# 35. SUPPLEMENTARY WORKBOOK
# ============================================================

supp_file <- file.path(
  out_dir,
  "supplementary",
  "Supplementary_Selected_Pathway_Sample_Level_LIMMA.xlsx"
)

wb <- openxlsx::createWorkbook()

supp_tables <- list(
  "Analysis_settings" = analysis_settings,
  "Sample_metadata" = sample_metadata_check,
  "Litter_balance" = litter_balance,
  "Selected_pathway_membership" = selected_pathway_metadata,
  "Selected_sample_values" = selected_sample_values,
  "Selected_LIMMA_results" = selected_stats,
  "ALL_pathway_LIMMA_results" = all_pathway_stats,
  "Selected_summary" = selected_pathway_summary,
  "Detection_summary" = detection_summary,
  "Protein_contributions" = protein_contribution_summary,
  "Group_raw_profile_QC" = raw_group_audit,
  "Group_imputed_discord" = imputed_group_audit
)

for (sheet_name in names(supp_tables)) {
  
  openxlsx::addWorksheet(
    wb,
    sheet_name
  )
  
  openxlsx::writeData(
    wb,
    sheet = sheet_name,
    x = supp_tables[[sheet_name]]
  )
  
}

openxlsx::saveWorkbook(
  wb,
  supp_file,
  overwrite = TRUE
)

message("Supplementary workbook saved: ", supp_file)
