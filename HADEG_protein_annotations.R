# ============================================================
# HADEG INSPECTION FILES
# Build Protein -> HADEG -> broad functional class table
# BEFORE making heatmaps
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(stringr)
})

# ============================================================
# 0. USER SETTINGS
# ============================================================

hadeg_dir <- "data/proteomics/HADEG_annotation_mapping"

hadeg_file <- file.path(
  hadeg_dir,
  "large_HADEG_matches_filtered_annotated.tsv"
)

out_dir <- file.path(hadeg_dir, "HADEG_inspection_files")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Optional filters
min_pident <- 30
max_evalue <- 1e-5

# ============================================================
# 1. READ HADEG FILE
# ============================================================

hadeg_raw <- read_tsv(
  hadeg_file,
  col_names = TRUE,
  show_col_types = FALSE,
  progress = FALSE
)

# Remove completely empty columns if present
hadeg_raw <- hadeg_raw %>%
  select(where(~ !all(is.na(.x))))

message("Columns found:")
print(colnames(hadeg_raw))

# ============================================================
# 2. STANDARDIZE COLUMN NAMES
# ============================================================

# This assumes your first columns are:
# Protein, HADEG_hit, pident, length, mismatch, gapopen,
# qstart, qend, sstart, send, evalue, bitscore,
# HADEG_accession, HADEG_gene, HADEG_description

if (!"Protein" %in% colnames(hadeg_raw)) {
  stop("No column named Protein was found.")
}

# If your last annotation columns came in unnamed/weird, rename by position
if (ncol(hadeg_raw) >= 15) {
  colnames(hadeg_raw)[1:15] <- c(
    "Protein",
    "HADEG_hit",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "HADEG_accession",
    "HADEG_gene",
    "HADEG_description"
  )
}

hadeg <- hadeg_raw %>%
  mutate(
    Protein = as.character(Protein),
    HADEG_hit = as.character(HADEG_hit),
    HADEG_accession = as.character(HADEG_accession),
    HADEG_gene = as.character(HADEG_gene),
    HADEG_description = as.character(HADEG_description),
    pident = suppressWarnings(as.numeric(pident)),
    evalue = suppressWarnings(as.numeric(evalue)),
    bitscore = suppressWarnings(as.numeric(bitscore)),
    search_text = str_to_lower(
      paste(
        HADEG_hit,
        HADEG_gene,
        HADEG_description,
        sep = " "
      )
    )
  )

# ============================================================
# 3. FILTER TO REASONABLE MATCHES
# ============================================================

hadeg_filtered <- hadeg %>%
  filter(
    !is.na(Protein),
    Protein != "",
    !is.na(evalue),
    evalue <= max_evalue,
    !is.na(pident),
    pident >= min_pident
  )

# ============================================================
# 4. FUNCTION TO ASSIGN BROAD HADEG CLASSES
# ============================================================

assign_hadeg_class <- function(x) {
  case_when(
    str_detect(x, "ahpc|ahpf|alkyl hydroperoxide|hydroperoxide|peroxiredoxin|catalase|thioredoxin|rubrerythrin|superoxide|oxidative stress") ~
      "Oxidative stress / ROS defense",
    
    str_detect(x, "monooxygenase|dioxygenase|hydroxylase|oxygenase|laccase|oxidase") ~
      "Hydrocarbon activation / oxygenation",
    
    str_detect(x, "catechol|protocatechuate|muconate|benzoate|benzaldehyde|benzoyl|phthalate|terephthalate|gentisate|homogentisate|salicylate|aromatic|cyclohexadiene|dihydroxy.*diene") ~
      "Aromatic compound degradation",
    
    str_detect(x, "naphthalene|phenanthrene|anthracene|fluorene|pyrene|biphenyl|pah|polycyclic") ~
      "PAH-related degradation",
    
    str_detect(x, "alkane|alkan|cyclohexanol|cyclohexanone|cyclohex|chn") ~
      "Alkane / cycloalkane metabolism",
    
    str_detect(x, "aldehyde dehydrogenase|alcohol dehydrogenase|acyl-coa|acetyl-coa|malonyl-coa|thiolase|transferase|transacylase|dehydrogenase|acetoacetate|beta-oxidation|fatty acid") ~
      "Central carbon rerouting",
    
    str_detect(x, "transporter|permease|efflux|abc|outer membrane|porin|uptake") ~
      "Transport / uptake / efflux",
    
    str_detect(x, "biosurfactant|surfactin|rhamnolipid|iturin|emulsan") ~
      "Biosurfactant / hydrocarbon access",
    
    TRUE ~ "Needs manual review"
  )
}

assign_hadeg_subclass <- function(x) {
  case_when(
    str_detect(x, "ahpc|ahpf|alkyl hydroperoxide|hydroperoxide") ~
      "Alkyl hydroperoxide detoxification",
    
    str_detect(x, "peroxiredoxin") ~
      "Peroxiredoxin system",
    
    str_detect(x, "catalase") ~
      "Catalase peroxide detoxification",
    
    str_detect(x, "thioredoxin") ~
      "Thioredoxin redox system",
    
    str_detect(x, "monooxygenase") ~
      "Monooxygenase-mediated oxidation",
    
    str_detect(x, "dioxygenase") ~
      "Dioxygenase-mediated oxidation",
    
    str_detect(x, "catechol") ~
      "Catechol ring-cleavage route",
    
    str_detect(x, "protocatechuate") ~
      "Protocatechuate route",
    
    str_detect(x, "benzaldehyde|benzoate|benzoyl") ~
      "Benzoate / benzaldehyde route",
    
    str_detect(x, "naphthalene") ~
      "Naphthalene-related route",
    
    str_detect(x, "phenanthrene") ~
      "Phenanthrene-related route",
    
    str_detect(x, "cyclohexanol|cyclohexanone|cyclohex|chn") ~
      "Cyclohexane / cyclohexanol route",
    
    str_detect(x, "aldehyde dehydrogenase") ~
      "Aldehyde oxidation",
    
    str_detect(x, "alcohol dehydrogenase") ~
      "Alcohol oxidation",
    
    str_detect(x, "acyl-coa|acetyl-coa|malonyl-coa|thiolase|transacylase|transferase") ~
      "CoA / acyl-transfer metabolism",
    
    TRUE ~ "Needs manual review"
  )
}

# ============================================================
# 5. CREATE HADEG MEMBERSHIP TABLE
# ============================================================

hadeg_membership <- hadeg_filtered %>%
  mutate(
    Broad_Class_auto = assign_hadeg_class(search_text),
    Subclass_auto = assign_hadeg_subclass(search_text),
    Manual_Broad_Class = Broad_Class_auto,
    Manual_Subclass = Subclass_auto,
    Include_in_heatmap = ifelse(Broad_Class_auto == "Needs manual review", "REVIEW", "YES"),
    Reviewer_notes = ""
  ) %>%
  select(
    Protein,
    HADEG_hit,
    HADEG_accession,
    HADEG_gene,
    HADEG_description,
    Broad_Class_auto,
    Subclass_auto,
    Manual_Broad_Class,
    Manual_Subclass,
    Include_in_heatmap,
    Reviewer_notes,
    pident,
    length,
    evalue,
    bitscore,
    everything(),
    -search_text
  )

# ============================================================
# 6. KEEP BEST HADEG HIT PER PROTEIN
#    Highest bitscore, then lowest evalue
# ============================================================

hadeg_best_per_protein <- hadeg_membership %>%
  arrange(Protein, desc(bitscore), evalue) %>%
  group_by(Protein) %>%
  slice(1) %>%
  ungroup()

# ============================================================
# 7. SUMMARY TABLES FOR INSPECTION
# ============================================================

hadeg_class_summary <- hadeg_best_per_protein %>%
  count(
    Broad_Class_auto,
    Subclass_auto,
    name = "n_proteins"
  ) %>%
  arrange(Broad_Class_auto, desc(n_proteins))

hadeg_gene_summary <- hadeg_best_per_protein %>%
  count(
    Broad_Class_auto,
    HADEG_gene,
    HADEG_description,
    name = "n_proteins"
  ) %>%
  arrange(Broad_Class_auto, desc(n_proteins))

hadeg_manual_review_only <- hadeg_best_per_protein %>%
  filter(
    Broad_Class_auto == "Needs manual review" |
      Subclass_auto == "Needs manual review"
  )

# ============================================================
# 8. WRITE OUTPUT FILES
# ============================================================

write_csv(
  hadeg_membership,
  file.path(out_dir, "01_HADEG_all_hits_auto_classified_for_review.csv")
)

write_csv(
  hadeg_best_per_protein,
  file.path(out_dir, "02_HADEG_best_hit_per_protein_for_manual_review.csv")
)

write_csv(
  hadeg_class_summary,
  file.path(out_dir, "03_HADEG_class_summary_counts.csv")
)

write_csv(
  hadeg_gene_summary,
  file.path(out_dir, "04_HADEG_gene_description_summary_counts.csv")
)

write_csv(
  hadeg_manual_review_only,
  file.path(out_dir, "05_HADEG_needs_manual_review_only.csv")
)

message("Done.")
message("Inspection files saved to: ", out_dir)