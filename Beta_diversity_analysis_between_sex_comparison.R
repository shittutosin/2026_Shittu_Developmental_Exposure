
suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(tidyverse)
  library(RColorBrewer)
  library(writexl)
  library(readxl)
  library(corpcor)
})


# =============================================================================
# 1. PATHS AND SETTINGS
# =============================================================================
base_dir <- "/base_dira"
metadata_file <- file.path(base_dir, "metadata.txt")
litter_analysis_root <- "litter_Id_folder"
litter_file <- file.path(litter_analysis_root, "Litter ID.xlsx")
litter_sheet <- "Table S1"
out_dir <- file.path(litter_analysis_root, "16S_Analysis", "Combined_Beta_Diversity",
                     "Schaefer_Strimmer_Contrast_Specific_Calibrated")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop("Cannot create output folder: ", out_dir)

set.seed(123)
timepoint_levels <- c("5W", "8W", "12W")

# =============================================================================
# 2. LOAD AND CHECK METADATA; EXPLICIT MOUSE-TO-LITTER MATCHING
# =============================================================================
if (!file.exists(metadata_file)) stop("Combined 16S metadata missing: ", metadata_file)
if (!file.exists(litter_file)) stop("Original litter workbook missing: ", litter_file)

meta <- read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
required_meta <- c("SampleID", "Sex", "Treatment", "Timepoint")
if (!all(required_meta %in% names(meta))) {
  stop("16S metadata missing columns: ", paste(setdiff(required_meta, names(meta)), collapse = ", "))
}

normalize_animal_id <- function(x) toupper(gsub("[^[:alnum:]]", "", trimws(as.character(x))))
normalize_sex <- function(x) {
  y <- tolower(trimws(as.character(x)))
  dplyr::case_when(y %in% c("m", "male", "males") ~ "M",
                   y %in% c("f", "female", "females") ~ "F",
                   TRUE ~ NA_character_)
}
normalize_treatment <- function(x) {
  y <- tolower(trimws(as.character(x)))
  dplyr::case_when(y %in% c("control", "c") ~ "Control",
                   y %in% c("dosed", "dose", "d") ~ "Dosed",
                   TRUE ~ NA_character_)
}
normalize_week <- function(x) {
  y <- tolower(gsub("[^[:alnum:]]", "", as.character(x)))
  y <- sub("weeks$", "w", y)
  y <- sub("week$", "w", y)
  toupper(y)
}

meta <- meta %>%
  mutate(
    SampleID = trimws(as.character(SampleID)),
    Sex = normalize_sex(Sex),
    Treatment = normalize_treatment(Treatment),
    Timepoint = normalize_week(Timepoint)
  )
if (anyNA(meta$SampleID) || any(meta$SampleID == "") || anyDuplicated(meta$SampleID)) {
  stop("Missing/duplicate 16S SampleID values in combined metadata.")
}
if (anyNA(meta$Sex) || anyNA(meta$Treatment) ||
    anyNA(meta$Timepoint) || any(!meta$Timepoint %in% timepoint_levels)) {
  stop("Combined metadata contains invalid Sex, Treatment, or Timepoint values.")
}

# The combined 16S files are named, for example:
#   FC1_rel-abundance_5W.tsv  ->  SampleID FC1_rel-abundance_5W
#                              ->  workbook animal ID 5W_FC1
sample_pattern <- "^(.+)_rel-abundance_(5W|8W|12W)$"
if (any(!grepl(sample_pattern, meta$SampleID, ignore.case = TRUE))) {
  stop("Unexpected combined SampleID format. Expected FC1_rel-abundance_5W etc. Bad: ",
       paste(meta$SampleID[!grepl(sample_pattern, meta$SampleID, ignore.case = TRUE)], collapse = ", "))
}
meta$SampleShortID <- sub(sample_pattern, "\\1", meta$SampleID, ignore.case = TRUE)
meta$TimepointFromID <- toupper(sub(sample_pattern, "\\2", meta$SampleID, ignore.case = TRUE))
if (any(meta$Timepoint != meta$TimepointFromID)) {
  stop("Timepoint in a 16S SampleID conflicts with its metadata Timepoint.")
}
meta$AnimalID <- paste0(meta$Timepoint, "_", toupper(meta$SampleShortID))
meta$JoinAnimalID <- normalize_animal_id(meta$AnimalID)
if (anyDuplicated(meta$JoinAnimalID)) {
  stop("More than one combined 16S SampleID maps to the same mouse.")
}

lit_raw <- readxl::read_excel(litter_file, sheet = litter_sheet, .name_repair = "minimal")
required_lit <- c("Animal ID", "Litter ID", "Sex", "Treatment", "Timepoint")
if (!all(required_lit %in% names(lit_raw))) {
  stop("Original Table S1 missing columns: ", paste(setdiff(required_lit, names(lit_raw)), collapse = ", "))
}

lit_lookup <- lit_raw %>%
  transmute(
    SourceAnimalID = as.character(.data[["Animal ID"]]),
    JoinAnimalID = normalize_animal_id(.data[["Animal ID"]]),
    LitterID = trimws(as.character(.data[["Litter ID"]])),
    SourceSex = normalize_sex(.data[["Sex"]]),
    SourceTreatment = normalize_treatment(.data[["Treatment"]]),
    SourceTimepoint = normalize_week(.data[["Timepoint"]])
  ) %>%
  filter(!is.na(JoinAnimalID), JoinAnimalID != "", JoinAnimalID != "NA")
if (anyDuplicated(lit_lookup$JoinAnimalID)) {
  stop("Duplicate animal IDs in original Table S1 after normalization.")
}


existing_lit <- names(meta)[tolower(gsub("[^a-z0-9]", "", tolower(names(meta)))) == "litterid"]
if (length(existing_lit) > 1L) stop("Multiple metadata litter-ID columns detected.")
if (length(existing_lit) == 1L) {
  meta$MetadataLitterID <- as.character(meta[[existing_lit]])
  meta[[existing_lit]] <- NULL
} else {
  meta$MetadataLitterID <- NA_character_
}

meta <- meta %>% left_join(lit_lookup, by = "JoinAnimalID", relationship = "many-to-one") %>%
  mutate(
    MatchStatus = case_when(
      is.na(SourceAnimalID) ~ "NO WORKBOOK MOUSE MATCH",
      is.na(LitterID) | LitterID == "" | tolower(LitterID) == "na" ~ "MISSING LITTER ID",
      is.na(SourceSex) ~ "MISSING SOURCE SEX",
      is.na(SourceTreatment) ~ "MISSING SOURCE TREATMENT",
      is.na(SourceTimepoint) | SourceTimepoint == "" ~ "MISSING SOURCE TIMEPOINT",
      Sex != SourceSex ~ "SEX MISMATCH",
      Treatment != SourceTreatment ~ "TREATMENT MISMATCH",
      Timepoint != SourceTimepoint ~ "TIMEPOINT MISMATCH",
      !is.na(MetadataLitterID) & trimws(MetadataLitterID) != "" &
        tolower(trimws(MetadataLitterID)) != "na" &
        trimws(MetadataLitterID) != LitterID ~ "METADATA LITTER CONFLICT",
      TRUE ~ "PASS"
    )
  )

readr::write_csv(
  meta %>% select(SampleID, SampleShortID, AnimalID, LitterID, Sex, Treatment,
                  Timepoint, SourceAnimalID, SourceSex, SourceTreatment,
                  SourceTimepoint, MetadataLitterID, MatchStatus),
  file.path(out_dir, "combined_16S_mouse_litter_matching_audit.csv"), na = ""
)
cat("\nMouse-to-litter verification status:\n")
print(table(meta$MatchStatus, useNA = "ifany"))
if (anyNA(meta$MatchStatus) || any(meta$MatchStatus != "PASS")) {
  stop("Combined 16S litter matching failed. Review matching audit CSV in: ", out_dir)
}
if (any(meta %>% distinct(LitterID, Treatment) %>% count(LitterID) %>% pull(n) > 1L)) {
  stop("A litter appears in both Control and Dosed. Verify Table S1.")
}

# =============================================================================
# 3. LOAD THE ORIGINAL COMBINED 16S ABUNDANCE TABLES
# =============================================================================
files <- list.files(base_dir, pattern = "_rel-abundance_(5W|8W|12W)\\.tsv$",
                    full.names = TRUE)
if (!length(files)) stop("No combined 16S relative-abundance TSV files in: ", base_dir)
samples <- sub("\\.tsv$", "", basename(files))
if (anyDuplicated(samples)) stop("Duplicate abundance filenames after removing .tsv")
if (!setequal(samples, meta$SampleID)) {
  cat("\nFILES NOT IN METADATA:\n")
  print(setdiff(samples, meta$SampleID))
  cat("\nMETADATA NOT IN FILES:\n")
  print(setdiff(meta$SampleID, samples))
  stop("Combined 16S metadata and abundance file IDs do not match. No samples silently dropped.")
}

long_df <- purrr::map2_dfr(files, samples, function(f, id) {
  d <- readr::read_tsv(f, col_types = readr::cols(), show_col_types = FALSE)
  if (!all(c("tax_id", "estimated counts") %in% names(d))) {
    stop("Missing tax_id or estimated counts in: ", basename(f))
  }
  d %>%
    transmute(tax_id = as.character(tax_id),
              est_counts = as.numeric(.data[["estimated counts"]]), sample = id)
}) %>%
  filter(!is.na(tax_id), !tax_id %in% c("mapped_unclassified", "unmapped")) %>%
  mutate(est_counts = replace_na(est_counts, 0)) %>%
  group_by(tax_id, sample) %>%
  summarise(est_counts = est_counts[1], .groups = "drop")

otu_counts <- long_df %>%
  tidyr::pivot_wider(names_from = sample, values_from = est_counts, values_fill = 0) %>%
  arrange(tax_id)
count_mat <- as.matrix(otu_counts[, -1, drop = FALSE])
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- otu_counts$tax_id
count_mat <- count_mat[, meta$SampleID, drop = FALSE]
if (!identical(colnames(count_mat), meta$SampleID)) stop("Sample alignment error.")
if (any(!is.finite(count_mat)) || any(count_mat < 0)) stop("Bad abundance counts.")
if (any(colSums(count_mat) <= 0)) stop("At least one sample has zero total counts.")

meta_ps <- meta %>%
  select(SampleID, Sex, Treatment, Timepoint, LitterID) %>%
  mutate(
    Sex = factor(Sex, levels = c("M", "F")),
    Treatment = factor(Treatment, levels = c("Control", "Dosed")),
    Timepoint = factor(Timepoint, levels = timepoint_levels)
  ) %>%
  tibble::column_to_rownames("SampleID")

ps <- phyloseq::phyloseq(
  phyloseq::otu_table(count_mat, taxa_are_rows = TRUE),
  phyloseq::sample_data(meta_ps)
)
ps_rel <- phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
ps_rel <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps_rel) > 0, ps_rel)

cat("\nSequenced-mouse distribution by timepoint, group, sex, and litter:\n")
litter_design <- meta %>% count(Timepoint, Treatment, Sex, LitterID, name = "n_mice") %>%
  arrange(Timepoint, Treatment, LitterID, Sex)
print(tibble::as_tibble(litter_design), n = Inf)


select_phyloseq_metadata <- function(ps_obj, column, value) {
  md <- as.data.frame(phyloseq::sample_data(ps_obj))
  if (!column %in% colnames(md)) {
    stop("Sample metadata column not found: ", column)
  }
  wanted <- !is.na(md[[column]]) & as.character(md[[column]]) == as.character(value)
  sample_ids <- rownames(md)[wanted]
  if (!length(sample_ids)) {
    stop("No samples found for ", column, " = ", value)
  }
  phyloseq::prune_samples(sample_ids, ps_obj)
}



ps_names <- phyloseq::sample_names(ps_rel)


if (!all(ps_names %in% rownames(meta_ps))) {
  stop("One or more phyloseq sample IDs are absent from verified metadata.")
}
md_source <- meta_ps[ps_names, c("Sex", "Treatment", "Timepoint", "LitterID"),
                     drop = FALSE]
md <- base::data.frame(
  Sex = as.character(md_source[["Sex"]]),
  Treatment = as.character(md_source[["Treatment"]]),
  Timepoint = as.character(md_source[["Timepoint"]]),
  LitterID = as.character(md_source[["LitterID"]]),
  row.names = ps_names,
  stringsAsFactors = FALSE
)
if (!identical(rownames(md), ps_names) || !identical(class(md), "data.frame")) {
  stop("Sample metadata must be an ordinary data.frame aligned to sample IDs.")
}

md$Sex <- factor(as.character(md$Sex), levels = c("M", "F"))
md$Treatment <- factor(as.character(md$Treatment), levels = c("Control", "Dosed"))
md$Timepoint <- factor(as.character(md$Timepoint), levels = timepoint_levels)
md$LitterID <- factor(as.character(md$LitterID))
if (anyNA(md[, c("Sex", "Treatment", "Timepoint", "LitterID")])) {
  stop("A required analysis variable is missing after litter matching.")
}

n <- nrow(md)
if (n < 15L) stop("Too few sequenced mice for the 12-cell multivariate model.")
D <- as.matrix(phyloseq::distance(ps_rel, method = "bray"))
D <- D[ps_names, ps_names, drop = FALSE]
if (any(!is.finite(D))) stop("Nonfinite Bray-Curtis distances.")
if (!isTRUE(all.equal(D, t(D), tolerance = 1e-9))) {
  stop("Bray-Curtis distance matrix is not symmetric.")
}


pc <- stats::cmdscale(stats::as.dist(D), eig = TRUE, add = TRUE, k = n - 1L)
if (is.null(pc$points)) stop("Full-dimensional PCoA could not be calculated.")
Y <- as.matrix(pc$points)
if (is.null(dim(Y))) Y <- matrix(Y, ncol = 1L)
if (is.null(rownames(Y))) rownames(Y) <- ps_names
Y <- Y[ps_names, , drop = FALSE]
if (any(!is.finite(Y)) || ncol(Y) < 2L) {
  stop("Nonfinite or insufficient full PCoA coordinates.")
}
K <- ncol(Y)
cat("\nFull Bray-Curtis coordinate dimensions used for the joint sex tests: ", K, "\n")
cat("Cailliez additive correction constant: ", pc$ac, "\n")

axis_eigen <- as.numeric(pc$eig)
axis_eigen <- axis_eigen[is.finite(axis_eigen) & axis_eigen > 1e-10]
if (length(axis_eigen) < 2L) stop("Too few positive Bray-Curtis PCoA axes.")
axis1_pct <- round(100 * axis_eigen[1] / sum(axis_eigen), 1)
axis2_pct <- round(100 * axis_eigen[2] / sum(axis_eigen), 1)

coords <- tibble::tibble(
  SampleID = ps_names,
  Axis.1 = Y[, 1], Axis.2 = Y[, 2],
  Sex = md$Sex, Treatment = md$Treatment,
  Timepoint = md$Timepoint, LitterID = md$LitterID
)


md$Cell <- factor(paste(md$Treatment, md$Timepoint, md$Sex, sep = "_"))
cell_levels <- levels(md$Cell)
X <- vapply(cell_levels, function(z) as.numeric(as.character(md$Cell) == z),
            numeric(nrow(md)))
colnames(X) <- paste0("Cell", cell_levels)
rownames(X) <- rownames(md)
if (any(rowSums(X) != 1L) || any(colSums(X) == 0L) || qr(X)$rank != ncol(X))
  stop("Bad 12-cell mean-model design")


valid_axis <- which(is.finite(pc$eig[seq_len(ncol(Y))]) &
   pc$eig[seq_len(ncol(Y))] > 1e-10 * max(pc$eig, na.rm = TRUE))
Y <- Y[, valid_axis, drop = FALSE]
K <- ncol(Y)
cat("\nNumber of strictly positive corrected Bray axes:", K, "\n")

# Recreate display axes using the kept coordinates (same first two axes).
coords$Axis.1 <- Y[, 1L]
coords$Axis.2 <- Y[, 2L]

treatment_levels <- c("Control", "Dosed")
idx_by_tr <- setNames(lapply(treatment_levels, function(z)
  which(as.character(md$Treatment) == z)), treatment_levels)
x_by_tr <- setNames(lapply(treatment_levels, function(z) {
  xx <- X[idx_by_tr[[z]], , drop = FALSE]
  xx[, colSums(xx) > 0L, drop = FALSE]
}), treatment_levels)
litter_by_tr <- setNames(lapply(treatment_levels, function(z)
  as.character(md$LitterID[idx_by_tr[[z]]])), treatment_levels)

fit_one_treatment <- function(Yg, Xg, lit) {
  ng <- nrow(Yg); pg <- ncol(Xg); kg <- ncol(Yg)
  if (ng <= pg || any(colSums(Xg) < 1) || qr(Xg)$rank != pg)
    stop("Treatment-specific design not estimable")
  ZZt <- outer(lit, lit, "==") * 1.0
  df <- ng - pg
  objective <- function(lambda, details = FALSE) {
    V0 <- diag(ng) + lambda * ZZt
    ch <- tryCatch(chol(V0), error = function(e) NULL)
    if (is.null(ch)) return(if (details) NULL else Inf)
    low <- t(ch)
    Xw <- forwardsolve(low, Xg)
    Yw <- forwardsolve(low, Yg)
    qq <- qr(Xw)
    if (qq$rank != pg) return(if (details) NULL else Inf)
    E <- qr.resid(qq, Yw)
    rss <- sum(E * E)
    if (!is.finite(rss) || rss <= 0) return(if (details) NULL else Inf)
    ldV <- 2 * sum(log(diag(ch)))
    ldXtVX <- 2 * sum(log(abs(diag(qr.R(qq)))))
    nll <- kg * (ldV + ldXtVX) + kg * df * log(rss / (kg * df))
    if (!details) return(nll)
    se2 <- rss / (kg * df)
    list(lambda = lambda, residual_var = se2,
         litter_var = lambda * se2, ICC = lambda/(1 + lambda),
         Xw = Xw, Yw = Yw, E = E, Q = qr.Q(qq),
         beta = qr.coef(qq, Yw), RSS = rss,
         n = ng, p = pg, nll = nll)
  }
  opt <- optimize(function(t) objective(exp(t)), c(-13, 10), tol = 1e-5)
  cand <- c(0, exp(opt$minimum))
  vals <- vapply(cand, objective, numeric(1))
  ans <- objective(cand[which.min(vals)], TRUE)
  if (is.null(ans)) stop("Litter REML fit failed")
  ans
}

fit_all <- function(Y_input) {
  setNames(lapply(treatment_levels, function(z) {
    ii <- idx_by_tr[[z]]
    fit_one_treatment(Y_input[ii, , drop = FALSE],
                      x_by_tr[[z]], litter_by_tr[[z]])
  }), treatment_levels)
}

fits <- fit_all(Y)
model_var <- dplyr::bind_rows(lapply(treatment_levels, function(z) {
  f <- fits[[z]]
  tibble::tibble(Treatment = z, n_mice = length(idx_by_tr[[z]]),
    n_independent_litters = length(unique(litter_by_tr[[z]])),
    n_axes = K, LitterVariance = f$litter_var,
    ResidualVariance = f$residual_var, EstimatedLitterICC = f$ICC,
    LitterToResidualRatio = f$lambda,
    Assumption = "Treatment-specific scalar litter variance in corrected Bray coordinate space")
}))
cat("\nTREATMENT-SPECIFIC WORKING VARIANCE FITS:\n")
print(model_var, n = Inf)

contrast_components <- function(fit, Xg, tr, tp) {
  wanted <- paste0("Cell", tr, "_", tp)
  mi <- match(paste0(wanted, "_M"), colnames(Xg))
  fi <- match(paste0(wanted, "_F"), colnames(Xg))
  if (is.na(mi) || is.na(fi)) stop("Cannot identify sex contrast: ", wanted)
  X0 <- fit$Xw
  X0[, mi] <- X0[, mi] + X0[, fi]
  X0 <- X0[, -fi, drop = FALSE]
  q0 <- qr(X0)
  if (q0$rank != ncol(X0)) stop("Reduced mean model lost rank: ", wanted)
  Q0 <- qr.Q(q0)
  u <- fit$Xw[, fi, drop = FALSE] -
       Q0 %*% crossprod(Q0, fit$Xw[, fi, drop = FALSE])
  norm_u <- sqrt(sum(u * u))
  if (norm_u < 1e-12) stop("Sex contrast not estimable")
  u <- as.numeric(u / norm_u)
  E0 <- fit$Yw - Q0 %*% crossprod(Q0, fit$Yw)
  ss <- sum(crossprod(u, E0)^2)
  full_RSS <- sum((fit$Yw - fit$Q %*% crossprod(fit$Q, fit$Yw))^2)
  list(F = ss / (full_RSS / (fit$n - fit$p)),
       Partial_R2 = ss / sum(E0 * E0))
}

panel_grid <- expand.grid(Timepoint = timepoint_levels,
                          Treatment = treatment_levels,
                          stringsAsFactors = FALSE)
panel_grid <- panel_grid[order(match(panel_grid$Treatment, treatment_levels),
                               match(panel_grid$Timepoint, timepoint_levels)), ]
rownames(panel_grid) <- NULL
panel_grid$Key <- paste(panel_grid$Treatment, panel_grid$Timepoint, sep = "_")

scores_from_fit <- function(ff) {
  vapply(seq_len(nrow(panel_grid)), function(j) {
    tr <- panel_grid$Treatment[j]; tp <- panel_grid$Timepoint[j]
    contrast_components(ff[[tr]], x_by_tr[[tr]], tr, tp)$F
  }, numeric(1))
}
observed_F <- scores_from_fit(fits)
observed_R2 <- vapply(seq_len(nrow(panel_grid)), function(j) {
  tr <- panel_grid$Treatment[j]; tp <- panel_grid$Timepoint[j]
  contrast_components(fits[[tr]], x_by_tr[[tr]], tr, tp)$Partial_R2
}, numeric(1))


analysis_profile <- "pilot"  # "pilot" or "full"
if (!analysis_profile %in% c("pilot", "full")) stop("Invalid analysis_profile")
if (analysis_profile == "pilot") {
  B_reference <- 999L
  B_calibration_outer <- 30L
  B_calibration_inner <- 99L
  stress_scenarios <- c("ss_correlated", "cell_heterogeneous")
} else {
  B_reference <- 4999L
  B_calibration_outer <- 200L
  B_calibration_inner <- 499L
  stress_scenarios <- c("ss_correlated", "cell_heterogeneous",
                        "litter_axis_mismatch", "heavy_tailed")
}
cat("\nBootstrap profile:", analysis_profile,
    "; main replicates:", B_reference,
    "; calibration outer:", B_calibration_outer,
    "; calibration inner:", B_calibration_inner, "\n")

panel_grid <- expand.grid(Timepoint = timepoint_levels,
                          Treatment = treatment_levels,
                          stringsAsFactors = FALSE)
panel_grid <- panel_grid[order(match(panel_grid$Treatment, treatment_levels),
                               match(panel_grid$Timepoint, timepoint_levels)), ]
rownames(panel_grid) <- NULL
panel_grid$Key <- paste(panel_grid$Treatment, panel_grid$Timepoint, sep = "_")

observed_F <- vapply(seq_len(nrow(panel_grid)), function(j) {
  z <- panel_grid[j, ]
  contrast_components(fits[[z$Treatment]], x_by_tr[[z$Treatment]],
                      z$Treatment, z$Timepoint)$F
}, numeric(1))
observed_R2 <- vapply(seq_len(nrow(panel_grid)), function(j) {
  z <- panel_grid[j, ]
  contrast_components(fits[[z$Treatment]], x_by_tr[[z$Treatment]],
                      z$Treatment, z$Timepoint)$Partial_R2
}, numeric(1))


make_reduced_design <- function(tr, tp) {
  XX <- x_by_tr[[tr]]
  mi <- match(paste0("Cell", tr, "_", tp, "_M"), colnames(XX))
  fi <- match(paste0("Cell", tr, "_", tp, "_F"), colnames(XX))
  if (is.na(mi) || is.na(fi)) stop("Missing target sex cells: ", tr, " ", tp)
  XX[, mi] <- XX[, mi] + XX[, fi]
  XX <- XX[, -fi, drop = FALSE]
  if (qr(XX)$rank != ncol(XX)) stop("Reduced design rank failure")
  XX
}


estimate_ss_covariance <- function(f, tr, tp, quiet = TRUE) {
  q_all <- qr.Q(qr(f$Xw), complete = TRUE)
  q_res <- q_all[, seq.int(f$p + 1L, f$n), drop = FALSE]
  r <- crossprod(q_res, f$Yw)
  if (nrow(r) != f$n-f$p || ncol(r) != K)
    stop("Incorrect residual contrast dimensions for ", tr, " ", tp)
  S <- stats::cov(r)
  S <- (S + t(S))/2
  variances <- diag(S)
  if (any(!is.finite(variances)) || min(variances) <= 1e-14*max(variances))
    stop("Zero/near-zero residual-axis variance for ", tr, " ", tp,
         "; inspect the effective PCoA rank before analysis")
  cs <- corpcor::cor.shrink(r, verbose = !quiet)
  lam <- as.numeric(attr(cs, "lambda"))
  if (length(lam) != 1L || !is.finite(lam) || lam < 0 || lam > 1)
    stop("Invalid Schafer-Strimmer shrinkage estimate")
  C_ss <- as.matrix(cs)*outer(sqrt(variances),sqrt(variances))
  C_diag <- diag(variances, nrow = K)
  C_legacy <- .65*C_diag + .35*S
  normal <- function(C) {
    C <- (C+t(C))/2
    C <- C/mean(diag(C))
    if (any(!is.finite(C)) || inherits(try(chol(C), silent = TRUE), "try-error"))
      stop("Covariance not positive definite for ", tr, " ", tp)
    C
  }
  list(
    shape = list(ss_correlated = normal(C_ss),
                 diagonal = normal(C_diag),
                 legacy35 = normal(C_legacy)),
    lambda = lam, residual_contrasts = nrow(r),
    axis_variances = variances)
}


fit_contrast_null <- function(y_group, tr, tp) {
  X0 <- make_reduced_design(tr, tp)
  f0 <- fit_one_treatment(y_group, X0, litter_by_tr[[tr]])
  sh <- estimate_ss_covariance(f0, tr, tp)
  null_mean <- X0 %*% f0$beta
  if (!identical(dim(null_mean), dim(y_group)))
    stop("Reduced-model null mean dimension mismatch")
  list(fit = f0, mean = null_mean, cov = sh, X0 = X0)
}


cell_scale <- list()
for (tr in treatment_levels) {
  f <- fits[[tr]]
  XX <- x_by_tr[[tr]]
  yy <- Y[idx_by_tr[[tr]], , drop = FALSE]
  rr <- rowSums((yy - XX %*% f$beta)^2)
  cell_means <- vapply(seq_len(ncol(XX)), function(k)
    mean(rr[XX[, k] == 1]), numeric(1))
  scalar <- pmin(1.7, pmax(.5, .5 + .5*cell_means/mean(cell_means)))
  scalar <- scalar/weighted.mean(scalar, colSums(XX))
  cell_scale[[tr]] <- as.vector(XX%*%scalar)
}


simulate_under_reduced <- function(null_obj, tr, scenario = "ss_correlated") {
  f0 <- null_obj$fit
  lit <- litter_by_tr[[tr]]
  uniq <- unique(lit)
  Cmouse <- null_obj$cov$shape$ss_correlated
  Clitter <- Cmouse
  if (scenario == "diagonal") {
    Cmouse <- null_obj$cov$shape$diagonal
    Clitter <- Cmouse
  } else if (scenario == "litter_axis_mismatch") {
    Clitter <- null_obj$cov$shape$diagonal
  } else if (scenario == "legacy35") {
    Cmouse <- null_obj$cov$shape$legacy35
    Clitter <- Cmouse
  } else if (!scenario %in% c("ss_correlated", "cell_heterogeneous", "heavy_tailed")) {
    stop("Unknown simulation scenario: ", scenario)
  }
  sample_innovations <- function(nr, mode, covar) {
    if (mode == "heavy_tailed") {
      # Standardized t(4): variance one, but more extreme observations.
      z <- matrix(stats::rt(nr*K, df = 4)/sqrt(2), nr, K)
    } else {
      z <- matrix(stats::rnorm(nr*K), nr, K)
    }
    z %*% chol(covar)
  }
  u <- sample_innovations(length(uniq), scenario, Clitter)
  e <- sample_innovations(length(lit), scenario, Cmouse)
  if (scenario == "cell_heterogeneous")
    e <- sweep(e, 1L, sqrt(cell_scale[[tr]]), "*")
  ystar <- null_obj$mean +
    sqrt(f0$litter_var)*u[match(lit, uniq), , drop = FALSE] +
    sqrt(f0$residual_var)*e
  if (any(!is.finite(ystar))) stop("Nonfinite synthetic PCoA scores")
  ystar
}

# Fitted full-model sex contrast for one simulated dataset. Both litter and
# mouse residual variance are re-estimated in each replicate.
refit_full_F <- function(ystar, tr, tp) {
  ff <- fit_one_treatment(ystar, x_by_tr[[tr]], litter_by_tr[[tr]])
  contrast_components(ff, x_by_tr[[tr]], tr, tp)$F
}

bootstrap_reference_F <- function(null_obj, tr, tp, B, seed,
                                  scenario = "ss_correlated", progress = FALSE) {
  set.seed(seed)
  ans <- numeric(B)
  for (b in seq_len(B)) {
    if (progress && b %% 250L == 0L)
      message(tr, " ", tp, "/", scenario, ": ", b, "/", B)
    ans[b] <- tryCatch(
      refit_full_F(simulate_under_reduced(null_obj, tr, scenario), tr, tp),
      error = function(e) stop("Bootstrap failure at ", tr, " ", tp,
                               " iteration ", b, ": ", conditionMessage(e)))
  }
  if (any(!is.finite(ans))) stop("Nonfinite simulated full-model F")
  ans
}
bootstrap_tail_p <- function(observed, reference) {
  (1 + sum(reference >= observed - 1e-10))/(length(reference)+1)
}

# ============================ Primary six tests ===============================
null_objects <- vector("list", nrow(panel_grid))
reference_primary <- matrix(NA_real_, nrow = B_reference,
                            ncol = nrow(panel_grid),
                            dimnames = list(NULL, panel_grid$Key))
result_rows <- list()
shrinkage_rows <- list()
axis_variance_rows <- list()
null_mean_checks <- list()
for (j in seq_len(nrow(panel_grid))) {
  tr <- panel_grid$Treatment[j]
  tp <- panel_grid$Timepoint[j]
  ii <- idx_by_tr[[tr]]
  null_obj <- fit_contrast_null(Y[ii, , drop = FALSE], tr, tp)
  null_objects[[j]] <- null_obj
  # Independently verify the target means are constrained, while OTHER
  # male/female timepoint means remain freely represented in the null design.
  for (t2 in timepoint_levels) {
    sub <- md[idx_by_tr[[tr]], , drop=FALSE]
    sub_m <- as.character(sub$Timepoint)==t2 & as.character(sub$Sex)=="M"
    sub_f <- as.character(sub$Timepoint)==t2 & as.character(sub$Sex)=="F"
    offset <- max(abs(colMeans(null_obj$mean[sub_m,,drop=FALSE]) -
                      colMeans(null_obj$mean[sub_f,,drop=FALSE])))
    if (t2==tp && offset>1e-9)
      stop("Null model failed to constrain sex means in ", tr, " ", tp)
    null_mean_checks[[length(null_mean_checks)+1L]] <- tibble::tibble(
      TestedTreatment=tr, TestedTimepoint=tp, CheckedTimepoint=t2,
      IsTestedContrast=(tp==t2),
      MaxAbsoluteMaleFemaleNullMeanDifference=offset)
  }
  message("CONTRAST-SPECIFIC null bootstrap: ", tr, " ", tp)
  nullF <- bootstrap_reference_F(null_obj, tr, tp, B_reference,
                                seed = 72000L + 100L*j, progress = TRUE)
  # Fixed, NONPRIMARY sensitivity comparisons; never select the method by p.
  nullF_diagonal <- bootstrap_reference_F(null_obj, tr, tp, B_reference,
                                seed = 82000L + 100L*j, scenario = "diagonal")
  nullF_legacy <- bootstrap_reference_F(null_obj, tr, tp, B_reference,
                                seed = 92000L + 100L*j, scenario = "legacy35")
  reference_primary[,j] <- nullF
  p_boot <- bootstrap_tail_p(observed_F[j], nullF)
  p_diag <- bootstrap_tail_p(observed_F[j], nullF_diagonal)
  p_legacy <- bootstrap_tail_p(observed_F[j], nullF_legacy)
  n_exc <- sum(nullF >= observed_F[j]-1e-10)
  mc <- stats::binom.test(n_exc, length(nullF))$conf.int
  one_md <- md[as.character(md$Treatment) == tr &
               as.character(md$Timepoint) == tp, , drop=FALSE]
  tab <- table(one_md$LitterID, one_md$Sex)
  n_mixed <- sum(rowSums(tab > 0) == 2L)
  result_rows[[j]] <- tibble::tibble(
    Treatment = tr, Timepoint = tp,
    n_males = sum(one_md$Sex == "M"), n_females = sum(one_md$Sex == "F"),
    n_litters_this_panel = length(unique(one_md$LitterID)),
    n_mixed_sex_litters = n_mixed,
    n_independent_litters_treatment = length(unique(litter_by_tr[[tr]])),
    n_axes = K, pseudo_F_GLS = observed_F[j], Partial_R2_GLS = observed_R2[j],
    p_value = p_boot,
    p_diagonal_sensitivity = p_diag,
    p_previous35_sensitivity = p_legacy,
    n_null_bootstraps = B_reference,
    p_bootstrap_MC_interval_lower = mc[1],
    p_bootstrap_MC_interval_upper = mc[2],
    NullLitterVariance = null_obj$fit$litter_var,
    NullMouseResidualVariance = null_obj$fit$residual_var,
    NullLitterICC = null_obj$fit$ICC,
    NullSchaeferStrimmerLambda = null_obj$cov$lambda,
    NullCorrelationsRetained = 1-null_obj$cov$lambda,
    NullResidualContrastDF = null_obj$cov$residual_contrasts,
    InferenceStatus = if (n_mixed == 0L)
       "EXPLORATORY model-based: sex/litter confounding; no paired-sex litters"
     else "EXPLORATORY model-based: few independent mixed-sex litters")
  shrinkage_rows[[j]] <- tibble::tibble(
    Treatment = tr, Timepoint = tp, Model = "CONTRAST_SPECIFIC_NULL",
    EstimatedLambda_CorrelationShrunkToZero = null_obj$cov$lambda,
    RetainedEmpiricalCorrelationFraction = 1-null_obj$cov$lambda,
    n_orthogonal_residual_contrasts = null_obj$cov$residual_contrasts,
    n_litters = length(unique(litter_by_tr[[tr]])), n_axes = K)
  axis_variance_rows[[j]] <- tibble::tibble(
    Treatment = tr, Timepoint = tp, Axis = seq_len(K),
    ResidualContrastAxisVariance = null_obj$cov$axis_variances)
}
joint_results <- dplyr::bind_rows(result_rows)
joint_results$p_BH_all_six <- stats::p.adjust(joint_results$p_value, "BH")
joint_results <- joint_results %>%
  dplyr::group_by(Treatment) %>%
  dplyr::mutate(p_BH_within_treatment = stats::p.adjust(p_value, "BH")) %>%
  dplyr::ungroup()


bh_significance <- function(p) {
  if (!is.finite(p)) return(NA_character_)
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}
format_bh_p <- function(p) {
  if (!is.finite(p)) return(NA_character_)
  if (p < 0.0001) return("<0.0001")
  trimws(format.pval(p, digits = 3, eps = 0.0001))
}
joint_results$Significance_BH_all_six <- vapply(
  joint_results$p_BH_all_six, bh_significance, character(1))
joint_results$PlotLabel_BH_all_six <- vapply(seq_len(nrow(joint_results)), function(i) {
  p <- joint_results$p_BH_all_six[i]
  if (!is.finite(p)) return(NA_character_)
  paste0("adj. p = ", format_bh_p(p), " (",
         joint_results$Significance_BH_all_six[i], ")")
}, character(1))

shrinkage_results <- dplyr::bind_rows(shrinkage_rows)
axis_diagnostic <- dplyr::bind_rows(axis_variance_rows)
reference_ss <- reference_primary
null_mean_checks <- dplyr::bind_rows(null_mean_checks)
cat("\n===== CONTRAST-SPECIFIC SCHAEFER-STRIMMER P VALUES =====\n")
print(joint_results, n=Inf, width=Inf)


readr::write_csv(joint_results,
  file.path(out_dir, "SS_contrast_specific_sex_pvalues.csv"), na="")
readr::write_csv(joint_results %>%
  dplyr::select(Treatment, Timepoint, p_value, p_BH_all_six,
                Significance_BH_all_six, PlotLabel_BH_all_six),
  file.path(out_dir, "SS_plot_annotations_BH_adjusted.csv"), na="")
readr::write_csv(shrinkage_results,
  file.path(out_dir, "SS_contrast_specific_shrinkage.csv"), na="")


calibration_rows <- list()
calibration_run_details <- list()
kk <- 0L
for (mode_i in seq_along(stress_scenarios)) {
  mode <- stress_scenarios[mode_i]
  for (j in seq_len(nrow(panel_grid))) {
    tr <- panel_grid$Treatment[j]
    tp <- panel_grid$Timepoint[j]
    base_null <- null_objects[[j]]
    set.seed(100000L + mode_i*10000L + j*100L)
    rejections <- logical(B_calibration_outer)
    outer_p <- numeric(B_calibration_outer)
    outer_lambdas <- numeric(B_calibration_outer)
    for (b in seq_len(B_calibration_outer)) {
      if (b %% 10L == 0L)
        message("NESTED calibration ", mode, " ", tr, " ", tp,
                ": ", b, "/", B_calibration_outer)
      # We use a local seed for every outer/inner pair for reproducibility.
      set.seed(200000L + mode_i*100000L + j*10000L + b)
      y_outer <- simulate_under_reduced(base_null, tr, mode)
      f_outer <- fit_one_treatment(y_outer, x_by_tr[[tr]], litter_by_tr[[tr]])
      F_outer <- contrast_components(f_outer, x_by_tr[[tr]], tr, tp)$F
      h0_outer <- fit_contrast_null(y_outer, tr, tp)
      outer_lambdas[b] <- h0_outer$cov$lambda
      inner_F <- bootstrap_reference_F(h0_outer, tr, tp, B_calibration_inner,
         seed = 1000000L + mode_i*100000L + j*10000L + b)
      outer_p[b] <- bootstrap_tail_p(F_outer, inner_F)
      rejections[b] <- outer_p[b] <= 0.05
    }
    n_bad <- sum(rejections)
    ci <- stats::binom.test(n_bad, B_calibration_outer)$conf.int
    kk <- kk + 1L
    calibration_rows[[kk]] <- tibble::tibble(
      TrueNullScenario = mode, Treatment = tr, Timepoint = tp,
      n_outer_null_datasets = B_calibration_outer,
      n_inner_bootstraps_per_outer_dataset = B_calibration_inner,
      NumberOfFalsePositives = n_bad,
      EmpiricalFalsePositiveRate = n_bad/B_calibration_outer,
      Binomial95CI_Lower = ci[1], Binomial95CI_Upper = ci[2],
      MedianRefittedShrinkageLambda = stats::median(outer_lambdas),
      ReferenceMethod = "Panel-specific reduced-model SS parametric bootstrap",
      CalibrationStatus = if (ci[1] > 0.05)
        "EVIDENCE OF TYPE-I INFLATION (95% CI excludes 0.05)"
      else if (ci[2] < 0.05)
        "POTENTIALLY CONSERVATIVE (95% CI excludes 0.05)"
      else "Not distinguishable from 0.05 at current Monte Carlo precision")
    calibration_run_details[[kk]] <- tibble::tibble(
      TrueNullScenario = mode, Treatment = tr, Timepoint = tp,
      OuterSimulation = seq_len(B_calibration_outer),
      p_value = outer_p, RejectedAt005 = rejections,
      RefittedLambda = outer_lambdas)
  }
  # checkpoint after each scenario
  readr::write_csv(dplyr::bind_rows(calibration_rows),
                   file.path(out_dir, "SS_nested_null_calibration_CHECKPOINT.csv"))
}
calibration_results <- dplyr::bind_rows(calibration_rows)
nested_calibration <- calibration_results
calibration_replicates <- dplyr::bind_rows(calibration_run_details)
readr::write_csv(calibration_results,
                 file.path(out_dir, "SS_nested_null_calibration.csv"))
readr::write_csv(calibration_replicates,
                 file.path(out_dir, "SS_nested_calibration_individual_pvalues.csv"))
if (any(calibration_results$Binomial95CI_Lower > 0.05))
 
cat("\n===== COMPLETED NESTED TYPE-I ERROR CALIBRATION =====\n")
print(calibration_results, n=Inf, width=Inf)

influence_rows <- list(); jj <- 0L
for (tr in treatment_levels) {
  full_i <- idx_by_tr[[tr]]
  for (drop_lit in unique(litter_by_tr[[tr]])) {
    keep <- litter_by_tr[[tr]] != drop_lit
    partial <- tryCatch(fit_one_treatment(Y[full_i[keep],,drop=FALSE],
            x_by_tr[[tr]][keep,,drop=FALSE], litter_by_tr[[tr]][keep]),
            error=function(e) NULL)
    for (tp in timepoint_levels) {
      jj <- jj+1L
      val <- if (is.null(partial)) NA_real_ else
        tryCatch(contrast_components(partial, x_by_tr[[tr]][keep,,drop=FALSE],
                                     tr, tp)$F, error=function(e) NA_real_)
      obs <- observed_F[panel_grid$Treatment == tr & panel_grid$Timepoint == tp]
      influence_rows[[jj]] <- tibble::tibble(
        Treatment=tr, Timepoint=tp, ExcludedLitter=drop_lit,
        FullModelF=obs, AfterExcludingLitterF=val,
        F_RatioToFull=val/obs,
        Note="Descriptive leave-one-litter-out F; no p-value from full-design reference")
    }
  }
}
influence_results <- dplyr::bind_rows(influence_rows)


exact_single_sex_litter_check <- function(tr, tp) {
  ii <- which(as.character(md$Treatment)==tr & as.character(md$Timepoint)==tp)
  lit <- as.character(md$LitterID[ii]); sx <- as.character(md$Sex[ii])
  groups <- lapply(unique(lit), function(x) {
    rows <- which(lit==x)
    list(sex=unique(sx[rows]), centroid=colMeans(Y[ii[rows],,drop=FALSE]))
  })
  if (any(vapply(groups,function(x) length(x$sex)!=1L, logical(1))))
    return(tibble::tibble(Treatment=tr, Timepoint=tp,
              exact_litter_centroid_p=NA_real_, n_assignments=NA_integer_,
              Status="Not applicable: mixed-sex litter present"))
  ss <- vapply(groups, function(x)x$sex, character(1))
  zz <- do.call(rbind, lapply(groups, function(x)x$centroid))
  nm <- sum(ss=="M")
  if (nm<1 || nm>=nrow(zz))
    stop("One sex absent from litter-centroid check")
  distance_of_means <- function(idx_m) {
    a <- colMeans(zz[idx_m,,drop=FALSE]); b <- colMeans(zz[-idx_m,,drop=FALSE])
    sum((a-b)^2)
  }
  assignments <- utils::combn(seq_len(nrow(zz)),nm,simplify=FALSE)
  reference <- vapply(assignments, distance_of_means, numeric(1))
  actual <- distance_of_means(which(ss=="M"))
  tibble::tibble(Treatment=tr,Timepoint=tp,
     exact_litter_centroid_p=mean(reference >= actual-1e-12),
     n_assignments=length(assignments),
     Status="Sensitivity only: independent litters assumed exchangeable; coarse p")
}
centroid_check <- dplyr::bind_rows(lapply(treatment_levels, function(tr)
  dplyr::bind_rows(lapply(timepoint_levels, function(tp)
    exact_single_sex_litter_check(tr,tp)))))

# =============================================================================
# 7. PLOTS 5 & 6: BH-ADJUSTED p (six tests) + SIGNIFICANCE SYMBOL
# =============================================================================
timepoint_cols <- RColorBrewer::brewer.pal(3, "Dark2")
names(timepoint_cols) <- timepoint_levels

make_joint_plot <- function(treatment_name, file_name) {
  this_coords <- coords %>% filter(Treatment == treatment_name) %>%
    mutate(Timepoint = factor(Timepoint, levels = timepoint_levels),
           Timepoint_Sex = interaction(Timepoint, Sex, sep = "_", drop = TRUE))
  labels <- joint_results %>% filter(Treatment == treatment_name) %>%
    mutate(Timepoint = factor(Timepoint, levels = timepoint_levels),
           label = PlotLabel_BH_all_six) %>%
    filter(!is.na(label), label != "") %>% select(Timepoint, label)

  g <- ggplot2::ggplot(this_coords,
                       ggplot2::aes(x = Axis.1, y = Axis.2, color = Timepoint)) +
    ggplot2::stat_ellipse(
      ggplot2::aes(group = Timepoint_Sex, linetype = Sex),
      type = "t", linewidth = 1, level = 0.85, na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(shape = Sex), size = 4, alpha = 0.9) +
    ggplot2::facet_wrap(~Timepoint, nrow = 1, scales = "fixed") +
    ggplot2::geom_label(data = labels,
                        ggplot2::aes(x = Inf, y = Inf, label = label),
                        inherit.aes = FALSE, hjust = 1.05, vjust = 1.15,
                        size = 3.6, fontface = "bold", color = "black",
                        fill = "white", label.size = 0.25,
                        label.padding = grid::unit(0.22, "lines")) +
    ggplot2::scale_color_manual(values = timepoint_cols, drop = FALSE) +
    ggplot2::scale_shape_manual(values = c(M = 16, F = 17),
                                labels = c(M = "Male", F = "Female"), drop = FALSE) +
    ggplot2::scale_linetype_manual(values = c(M = "solid", F = "dashed"),
                                   labels = c(M = "Male", F = "Female"), drop = FALSE) +
    ggplot2::labs(title = NULL,
                  x = paste0("Axis 1 (", axis1_pct, "%)"),
                  y = paste0("Axis 2 (", axis2_pct, "%)"),
                  color = "Timepoint", shape = "Sex", linetype = "Sex") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.7),
      strip.text = ggplot2::element_text(face = "bold", size = 17),
      axis.title.x = ggplot2::element_text(face = "bold", size = 17),
      axis.title.y = ggplot2::element_text(face = "bold", size = 17),
      axis.text.x = ggplot2::element_text(face = "bold", size = 12),
      axis.text.y = ggplot2::element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold", size = 14),
      legend.text = ggplot2::element_text(face = "bold", size = 12),
      panel.spacing = grid::unit(1, "lines"),
      plot.margin = ggplot2::margin(10, 15, 10, 10)
    )
  ggplot2::ggsave(file.path(out_dir, file_name), plot = g,
                  device = "tiff", width = 9, height = 6.8, units = "in",
                  dpi = 600, compression = "lzw", bg = "white", limitsize = FALSE)
  g
}

p5 <- make_joint_plot("Control", "beta_5_controls_SS_contrast_specific_BH_adjusted.tiff")
p6 <- make_joint_plot("Dosed", "beta_6_dosed_SS_contrast_specific_BH_adjusted.tiff")
print(p5)
print(p6)

# =============================================================================
# 8. COMPLETE SCIENTIFIC AUDIT AND OUTPUTS
# =============================================================================
litter_design <- coords %>%
  count(Treatment, Timepoint, LitterID, Sex, name = "n_mice") %>%
  arrange(Treatment, Timepoint, LitterID, Sex)
cell_summary <- coords %>%
  group_by(Treatment, Timepoint, LitterID) %>%
  summarise(n_males = sum(Sex == "M"),
            n_females = sum(Sex == "F"),
            n_mice = dplyr::n(),
            .groups = "drop") %>%
  arrange(Treatment, Timepoint, LitterID)

writexl::write_xlsx(
  list(
    "SS_contrast_specific_sex_tests" = joint_results,
    "SS_null_shrinkage" = shrinkage_results,
    "Null_model_constraint_checks" = null_mean_checks,
    "Nested_calibration" = calibration_results,
    "Per_outer_calibration_p" = calibration_replicates,
    "Axis_variance_diagnostics" = axis_diagnostic,
    "Leave_one_litter_out_F" = influence_results,
    "Treatment_specific_full_ICC" = model_var,
    "Exact_litter_centroid_checks" = centroid_check,
    "Sequenced_litter_design" = litter_design,
    "Each_litter_sex_counts" = cell_summary,
    "Mouse_PCoA_coordinates" = coords,
    "Mouse_litter_mapping" = meta %>%
       select(SampleID, AnimalID, LitterID, Treatment,
              Timepoint, Sex, MatchStatus)
  ),
  path = file.path(out_dir, "Supplementary_16S_SS_Contrast_Specific.xlsx")
)



# EXPORT INPUTS FOR VALIDATING JOINT MULTIVARIATE MODEL

validation_dir <- file.path(
  out_dir, "Model_validation_inputs"
)

dir.create(
  validation_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 1. Export individual mouse metadata
metadata_export <- data.frame(
  SampleID = rownames(md),
  Sex = as.character(md$Sex),
  Treatment = as.character(md$Treatment),
  Timepoint = as.character(md$Timepoint),
  LitterID = as.character(md$LitterID)
)

write.csv(
  metadata_export,
  file.path(validation_dir, "16S_mouse_metadata.csv"),
  row.names = FALSE
)

# 2. Export the Bray-Curtis distance matrix
distance_export <- data.frame(
  SampleID = rownames(D),
  D,
  check.names = FALSE
)

write.csv(
  distance_export,
  file.path(validation_dir, "16S_Bray_Curtis_matrix.csv"),
  row.names = FALSE
)

cat("\nValidation files exported to:\n", validation_dir)

