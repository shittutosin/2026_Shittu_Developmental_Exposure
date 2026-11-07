
needed <- c('dplyr','tidyr','stringr','readr','readxl','ggplot2',
            'purrr','tibble','writexl')
missing_pkgs <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_pkgs))
  stop('Install these packages first: ',paste(missing_pkgs,collapse=', '))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(readr)
  library(readxl); library(ggplot2); library(purrr); library(tibble)
})

# 1. SETTINGS -----------------------------------------------------------------
root <- paste0('root_folder/',
               'root_folder')

# Same locations as your NEW individual-protein litter-blocked KO script.
limma_dir <- file.path(root,'Litter_analysis','Final_Litter_Blocked_Limma')
original_dir <- file.path(root,'Pathway_Enrichment_specific',
                          'high_confidence_hits','KO_summary',
                          'Pathway_intermediate_tables')
original_ko_file <- file.path(original_dir,'Protein_KO_pathway_membership.csv')
original_expr_file <- file.path(original_dir,'Expression_with_pathway_membership.csv')
# Restored output from your NEW litter analysis (new mouse intensities, original annotations).
litter_annotation_dir <- file.path(root,'Litter_analysis',
                                   'Litter_Pathways_Original_Annotation',
                                   'Pathway_intermediate_tables')
litter_expr_file <- file.path(litter_annotation_dir,
                              'Expression_with_pathway_membership.csv')


out_dir <- file.path(root,'Litter_analysis',
                     'Carbohydrate_Volcano_ExactRaw_Peptide_Groups_Litter_Blocked')
plots_dir <- file.path(out_dir,'plots')
tables_dir <- file.path(out_dir,'tables')
dir.create(plots_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(tables_dir,recursive=TRUE,showWarnings=FALSE)

fdr_cutoff <- 0.05
abs_logfc_cutoff <- 1.5
plot_dpi <- 600

expected <- c('5W_M','8W_M','12W_M','5W_F','8W_F','12W_F')
time_levels <- c('5W','8W','12W')
sex_levels <- c('Males','Females')

carbohydrate_pathways <- tibble(
  pathway_group='Carbohydrate metabolism',
  pathway_id=c('map00051','map00052','map00500'),
  pathway_short=c('Fructose and mannose metabolism',
                  'Galactose metabolism',
                  'Starch and sucrose metabolism')
)

# 2. HELPERS ------------------------------------------------------------------
check_columns <- function(x,required,source) {
  absent <- setdiff(required,names(x))
  if(length(absent)) stop(source,': missing columns: ',paste(absent,collapse=', '))
}

clean_filename <- function(x) {
  x <- str_replace_all(as.character(x),'[^A-Za-z0-9]+','_')
  str_replace_all(x,'^_+|_+$','')
}

normalize_sex <- function(z) {
  z <- toupper(trimws(as.character(z)))
  case_when(z %in% c('M','MALE','MALES') ~ 'M',
            z %in% c('F','FEMALE','FEMALES') ~ 'F',
            TRUE ~ NA_character_)
}

normalize_timepoint <- function(z) {
  z <- toupper(gsub('[^0-9A-Z]','',trimws(as.character(z))))
  case_when(z %in% c('5W','5WK','5WKS','5WEEK','5WEEKS') ~ '5W',
            z %in% c('8W','8WK','8WKS','8WEEK','8WEEKS') ~ '8W',
            z %in% c('12W','12WK','12WKS','12WEEK','12WEEKS') ~ '12W',
            TRUE ~ NA_character_)
}

# Same cleaning logic as your ORIGINAL KOfam-first volcano script.
clean_kofam_eggnog_function <- function(x) {
  x %>% as.character() %>%
    str_replace_all('\u00A0',' ') %>%
    str_replace('^\\s*[0-9.]+e[-+]?[0-9]+\\s+','') %>%
    str_replace('^\\s*[0-9]+\\s+','') %>%
    str_replace('\\s*\\[EC:[^\\]]+\\]','') %>%
    str_replace('\\s*;.*$','') %>%
    str_replace('^\\s*Belongs to the\\s+','') %>%
    str_replace('\\s+family\\.?$',' family') %>%
    str_replace_all('\\s+',' ') %>%
    str_squish()
}

# Same palette logic as the original plot, with a safe n=0 case.
make_distinct_colors <- function(n) {
  if(n < 1L) return(character())
  okabe_ito <- c('#E69F00','#56B4E9','#009E73','#0072B2',
                  '#D55E00','#CC79A7','#000000')
  dark2 <- c('#1B9E77','#D95F02','#7570B3','#E7298A',
             '#66A61E','#E6AB02','#A6761D')
  set1 <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3',
            '#FF7F00','#A65628','#F781BF')
  paired_strong <- c('#1F78B4','#33A02C','#E31A1C','#FF7F00',
                     '#6A3D9A','#B15928','#FB9A99','#FDBF6F','#CAB2D6')
  all_cols <- unique(c(okabe_ito,dark2,set1,paired_strong))
  if(n > length(all_cols)) {
    extra_n <- n-length(all_cols)
    extra <- grDevices::hcl(
      h=seq(15,375,length.out=extra_n+1L)[seq_len(extra_n)],
      c=100,l=45)
    all_cols <- c(all_cols,extra)
  }
  all_cols[seq_len(n)]
}

# 3. READ ORIGINAL KO/LABEL FILES AND NEW RESTORED MEMBERSHIP ---------------
for(fp in c(original_ko_file,original_expr_file,litter_expr_file))
  if(!file.exists(fp)) stop('Missing REQUIRED annotation file: ',fp)

map_raw <- readr::read_csv(original_ko_file,show_col_types=FALSE)
check_columns(map_raw,c('Protein','KO','pathway_id'),basename(original_ko_file))
if(!'Function_from_annotation' %in% names(map_raw))
  map_raw$Function_from_annotation <- NA_character_

# Original, confirmed PROTEIN -> KO -> PATHWAY triplets (no remapping).
original_ko_map <- map_raw %>%
  transmute(Protein=trimws(as.character(Protein)),
            KO=trimws(as.character(KO)),
            pathway_id=trimws(as.character(pathway_id)),
            Function_from_annotation=as.character(Function_from_annotation)) %>%
  filter(!is.na(Protein),nzchar(Protein),
         !is.na(KO),grepl('^K[0-9]{5}$',KO),
         !is.na(pathway_id),nzchar(pathway_id)) %>%
  distinct(Protein,KO,pathway_id,.keep_all=TRUE)

original_carb_map <- original_ko_map %>%
  semi_join(carbohydrate_pathways,by='pathway_id')
if(!nrow(original_carb_map))
  stop('No carbohydrate Protein-KO-pathway triplets in original mapping.')

expr_raw <- readr::read_csv(original_expr_file,show_col_types=FALSE)
check_columns(expr_raw,c('Protein','Sex','Timepoint','pathway_id'),
              basename(original_expr_file))
if(!'KO' %in% names(expr_raw)) expr_raw$KO <- NA_character_
if(!'Source' %in% names(expr_raw)) expr_raw$Source <- NA_character_
if(!'Function_from_annotation' %in% names(expr_raw))
  expr_raw$Function_from_annotation <- NA_character_

# The original dataset membership is retained to audit the NEW restored mapping.
# The actual plot uses litter_present (below), not this original presence table.
original_present <- expr_raw %>%
  transmute(Protein=trimws(as.character(Protein)),
            pathway_id=trimws(as.character(pathway_id)),
            SexCode=normalize_sex(Sex),
            Timepoint=normalize_timepoint(Timepoint)) %>%
  mutate(Dataset=paste(Timepoint,SexCode,sep='_')) %>%
  filter(Dataset %in% expected,!is.na(Protein),nzchar(Protein),
         !is.na(pathway_id),nzchar(pathway_id)) %>%
  distinct(Dataset,Protein,pathway_id)
if(!nrow(original_present))
  stop('No original dataset-specific pathway assignments were found.')


litter_expr_raw <- readr::read_csv(litter_expr_file,show_col_types=FALSE)
check_columns(litter_expr_raw,c('Dataset','Protein','pathway_id'),
              basename(litter_expr_file))
litter_present <- litter_expr_raw %>%
  transmute(Dataset=trimws(as.character(Dataset)),
            Protein=trimws(as.character(Protein)),
            pathway_id=trimws(as.character(pathway_id))) %>%
  filter(Dataset %in% expected,!is.na(Protein),nzchar(Protein),
         !is.na(pathway_id),nzchar(pathway_id)) %>%
  distinct(Dataset,Protein,pathway_id)
if(!nrow(litter_present))
  stop('No sex/timepoint-specific pathway membership in NEW litter-restored file.')

litter_only <- dplyr::anti_join(litter_present,original_present,
                               by=c('Dataset','Protein','pathway_id'))
readr::write_csv(litter_only,
  file.path(tables_dir,'QC_restored_memberships_missing_from_original.csv'))
if(nrow(litter_only))
  stop('NEW restored membership contains ',nrow(litter_only),
       ' unexpected dataset/protein/pathway assignments. See QC CSV; stopped.')


original_only <- dplyr::anti_join(original_present,litter_present,
                                 by=c('Dataset','Protein','pathway_id'))
readr::write_csv(original_only,
  file.path(tables_dir,'QC_original_memberships_absent_from_restored.csv'))
if(nrow(original_only))
  message('Note: ',nrow(original_only),
          ' original dataset/protein/pathway memberships were absent in ',
          'the restored litter-expression table (see QC CSV).')


expr_label_candidates <- expr_raw %>%
  transmute(Protein=trimws(as.character(Protein)),
            KO=trimws(as.character(KO)),
            pathway_id=trimws(as.character(pathway_id)),
            Source=as.character(Source),
            Function_from_annotation=as.character(Function_from_annotation)) %>%
  semi_join(original_carb_map,by=c('Protein','KO','pathway_id'))


map_label_candidates <- original_carb_map %>%
  transmute(Protein,KO,pathway_id,
            Source='Original protein-KO mapping (source unspecified)',
            Function_from_annotation=Function_from_annotation)

ko_preferred_function <- bind_rows(expr_label_candidates,map_label_candidates) %>%
  mutate(KO_function_clean=clean_kofam_eggnog_function(Function_from_annotation),
         Source=coalesce(Source,'Unknown original annotation source')) %>%
  filter(!is.na(KO_function_clean),nzchar(KO_function_clean),
         !str_detect(str_to_lower(KO_function_clean),'^psort location'),
         !str_detect(str_to_lower(KO_function_clean),'^belongs to')) %>%
  mutate(source_priority=case_when(
    str_detect(str_to_lower(Source),'kofam') ~ 1L,
    str_detect(str_to_lower(Source),'eggnog') ~ 2L,
    TRUE ~ 3L)) %>%
  arrange(KO,source_priority,KO_function_clean) %>%
  group_by(KO) %>%
  summarise(preferred_function=first(KO_function_clean),
            preferred_source=first(Source),
            all_functions_seen=paste(sort(unique(KO_function_clean)),collapse=' /// '),
            .groups='drop')

# One row per KO, including those with no usable descriptive annotation.
ko_preferred_function <- original_carb_map %>% distinct(KO) %>%
  left_join(ko_preferred_function,by='KO') %>%
  mutate(preferred_function=coalesce(na_if(preferred_function,''),'Unknown function'),
         preferred_source=coalesce(preferred_source,'No original annotation available')) %>%
  arrange(KO)

readr::write_csv(ko_preferred_function,
                 file.path(tables_dir,'KO_preferred_function_labels_KOfam_first.csv'))

# 5. LOAD NEW LITTER-BLOCKED PROTEIN-LEVEL LIMMA STATS ------------------------
read_litter_blocked <- function(dataset) {
  fp <- file.path(limma_dir,
                  paste0(dataset,'_Litter_Blocked_Protein_Results.xlsx'))
  if(!file.exists(fp)) stop('Missing litter-blocked limma file: ',fp)
  tabs <- readxl::excel_sheets(fp)
  if(!'Protein_Statistics' %in% tabs)
    stop(basename(fp),': missing Protein_Statistics sheet')
  tab <- readxl::read_excel(fp,sheet='Protein_Statistics')
  check_columns(tab,c('Protein','log2FC','p_value','BH_FDR'),basename(fp))
  tab$Protein <- trimws(as.character(tab$Protein))
  if(anyNA(tab$Protein) || any(!nzchar(tab$Protein)))
    stop(dataset,': missing protein identifiers in Protein_Statistics')
  if(anyDuplicated(tab$Protein))
    stop(dataset,': duplicated proteins in Protein_Statistics')
  parts <- strsplit(dataset,'_',fixed=TRUE)[[1]]
  tab %>% transmute(
    Dataset=dataset, Protein,
    Timepoint=parts[1],
    Sex=if(parts[2]=='M') 'Males' else 'Females',
    log2FC=suppressWarnings(as.numeric(log2FC)),
    p_value=suppressWarnings(as.numeric(p_value)),
    BH_FDR=suppressWarnings(as.numeric(BH_FDR)),
    source_file=basename(fp))
}

limma_all <- purrr::map_dfr(expected,read_litter_blocked)
if(!nrow(limma_all)) stop('No litter-blocked protein statistics loaded.')
if(anyDuplicated(limma_all[c('Dataset','Protein')]))
  stop('Duplicate dataset/protein combinations in limma input.')


joined_before_restriction <- limma_all %>%
  inner_join(original_carb_map,by='Protein',relationship='many-to-many')

matched_all <- joined_before_restriction %>%
  semi_join(litter_present,by=c('Dataset','Protein','pathway_id')) %>%
  inner_join(carbohydrate_pathways,by='pathway_id') %>%
  left_join(ko_preferred_function,by='KO') %>%
  distinct(Dataset,Protein,KO,pathway_id,.keep_all=TRUE) %>%
  mutate(
    # EXACT decision rule in your NEW litter-blocked KO matrix script.
    is_significant=is.finite(BH_FDR) & BH_FDR>=0 & BH_FDR<fdr_cutoff &
      is.finite(log2FC) & abs(log2FC)>=abs_logfc_cutoff,
    direction=case_when(
      is_significant & log2FC>0 ~ 'Up',
      is_significant & log2FC<0 ~ 'Down',
      TRUE ~ 'NS'),
   
    adj_P_value_plot=if_else(is.finite(BH_FDR) & BH_FDR>=0 & BH_FDR<=1,
                             pmax(BH_FDR,.Machine$double.xmin),NA_real_),
    neg_log10_adj_p=-log10(adj_P_value_plot),
    is_plot_valid=is.finite(log2FC) & is.finite(neg_log10_adj_p),
    KO_plot=if_else(is_significant,
                    paste0(KO,' | ',preferred_function),
                    'Not significant'),
    pathway_name=pathway_short,
    pathway_label=paste0(pathway_name,' (',pathway_id,')'),
    Timepoint=factor(Timepoint,levels=time_levels),
    Sex=factor(Sex,levels=sex_levels))

if(!nrow(matched_all))
  stop('No matched carbohydrate results. Check original dataset-specific pathway mapping.')
if(anyDuplicated(matched_all[c('Dataset','Protein','KO','pathway_id')]))
  stop('Duplicate dataset/protein/KO/pathway after matching!')


strict_lookup_file <- file.path(
  root, 'Litter_analysis',
  'Selected_pathway_ExactPeptides_ExactRawValues_Litter_Blocked_LIMMA',
  'tables', 'KO_protein_group_lookup_FROM_FILTERED.csv'
)
if (!file.exists(strict_lookup_file)) {
  stop('Missing STRICT pre-imputation KO grouping lookup: ', strict_lookup_file,
       '\nRun Selected_Pathways_Exact_Peptides_Exact_Filtered_Intensities_Litter_Blocked.R first.')
}

lookup <- readr::read_csv(
  strict_lookup_file,
  col_types=readr::cols(.default=readr::col_character()),
  show_col_types=FALSE
)
lookup_required <- c('Dataset','KO','Protein','ProteinGroupID','Representative',
                     'GroupSize','PeptideSet','OriginalObservedProfile')
check_columns(lookup, lookup_required, strict_lookup_file)
if (!'GroupOrganism' %in% names(lookup)) lookup$GroupOrganism <- NA_character_
lookup <- lookup %>%
  mutate(across(all_of(c('Dataset','KO','Protein','ProteinGroupID')),
                ~trimws(as.character(.x)))) %>%
  filter(Dataset %in% expected, !is.na(KO), nzchar(KO))
lookup$GroupSize <- suppressWarnings(as.integer(lookup$GroupSize))
if (anyNA(lookup$GroupSize) || any(lookup$GroupSize < 1L))
  stop('Invalid GroupSize values in strict lookup.')
profile_sample_counts <- vapply(strsplit(lookup$OriginalObservedProfile,
  '|',fixed=TRUE),length,integer(1))
if (anyNA(lookup$OriginalObservedProfile) ||
    any(profile_sample_counts!=8L))
  stop('Strict lookup contains an invalid original eight-sample profile.')
if (anyDuplicated(lookup[c('Dataset','KO','Protein')]))
  stop('Strict lookup has duplicate Dataset + KO + Protein keys.')
if (anyNA(lookup$ProteinGroupID) || any(!nzchar(lookup$ProteinGroupID)))
  stop('Strict lookup has missing protein-group IDs.')

# Revalidate each global group. Missing peptides must never collapse.
group_evidence_qc <- lookup %>%
  group_by(Dataset,KO,ProteinGroupID) %>%
  summarise(
    n_accessions=n_distinct(Protein),
    n_peptide_sets=n_distinct(PeptideSet),
    n_original_profiles=n_distinct(OriginalObservedProfile),
    missing_peptide=any(is.na(PeptideSet) | !nzchar(PeptideSet)),
    missing_original_profile=any(is.na(OriginalObservedProfile) |
                                 !nzchar(OriginalObservedProfile)),
    n_reported_group_sizes=n_distinct(GroupSize),
    reported_group_size=first(GroupSize),
    .groups='drop'
  )
readr::write_csv(group_evidence_qc,
  file.path(tables_dir,'QC_strict_original_peptide_intensity_groups.csv'))
bad_groups <- group_evidence_qc %>%
  filter(n_reported_group_sizes != 1L | reported_group_size != n_accessions |
         (n_accessions>1L &
          (n_peptide_sets!=1L | n_original_profiles!=1L |
           missing_peptide | missing_original_profile)))
if (nrow(bad_groups)) {
  readr::write_csv(bad_groups,
    file.path(tables_dir,'ERROR_conflicting_original_evidence_groups.csv'))
  stop('Strict grouping lookup contains inconsistent group sizes or original evidence.')
}

missing_group_lookup <- matched_all %>%
  anti_join(lookup, by=c('Dataset','KO','Protein'))
readr::write_csv(missing_group_lookup %>%
  select(Dataset,KO,Protein,pathway_id),
  file.path(tables_dir,'QC_carbohydrate_accessions_missing_group_lookup.csv'))
if (nrow(missing_group_lookup)) {
  stop(nrow(missing_group_lookup),
       ' carbohydrate accession/KO rows lack strict group evidence; ',
       'see QC_carbohydrate_accessions_missing_group_lookup.csv')
}


carb_members <- matched_all %>%
  left_join(lookup %>%
              select(Dataset,KO,Protein,ProteinGroupID,Representative,
                     GroupSize,PeptideSet,OriginalObservedProfile,GroupOrganism),
            by=c('Dataset','KO','Protein'), relationship='many-to-one')
if (nrow(carb_members)!=nrow(matched_all) ||
    anyNA(carb_members$ProteinGroupID)) {
  stop('Strict grouping join duplicated or lost carbohydrate accession records.')
}
if (anyDuplicated(carb_members[c('Dataset','KO','pathway_id','Protein')]))
  stop('Duplicate carbohydrate pathway accession after strict lookup join.')


plotted_group_keys <- carb_members %>%
  distinct(Dataset,KO,ProteinGroupID)
all_plotted_members <- lookup %>%
  semi_join(plotted_group_keys,by=c('Dataset','KO','ProteinGroupID')) %>%
  left_join(limma_all %>% select(Dataset,Protein,log2FC,p_value,BH_FDR),
            by=c('Dataset','Protein'),relationship='many-to-one')

finite_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) median(x) else NA_real_
}
max_range_or_inf <- function(x) {
  if (!length(x) || any(!is.finite(x))) return(Inf)
  max(x)-min(x)
}

# The original KO-level strict-group plots use 1e-8 for numerical agreement.
# Also require the identical significance classification among group members.
all_plotted_members <- all_plotted_members %>%
  mutate(Member_significant=is.finite(BH_FDR) & BH_FDR<fdr_cutoff &
           is.finite(log2FC) & abs(log2FC)>=abs_logfc_cutoff)
strict_group_stats <- all_plotted_members %>%
  group_by(Dataset,KO,ProteinGroupID) %>%
  summarise(
    N_Members_Full_Group=n_distinct(Protein),
    N_Missing_or_Invalid_Stats=sum(!is.finite(log2FC) |
                                    !is.finite(p_value) | !is.finite(BH_FDR)),
    MemberAccessions=paste(sort(unique(Protein)),collapse='; '),
    Max_FC_Difference=max_range_or_inf(log2FC),
    Max_P_Difference=max_range_or_inf(p_value),
    Max_FDR_Difference=max_range_or_inf(BH_FDR),
    N_Significance_Classes=n_distinct(Member_significant),
    Group_log2FC_descriptive=finite_median(log2FC),
    .groups='drop'
  ) %>%
  mutate(
    StatsAgree=N_Missing_or_Invalid_Stats==0L &
      Max_FC_Difference<=1e-8 & Max_P_Difference<=1e-8 &
      Max_FDR_Difference<=1e-8 & N_Significance_Classes==1L,
    NeedsReview=!StatsAgree,
    StatisticalStatus=case_when(
      N_Missing_or_Invalid_Stats>0L ~
        'Some group accessions lack valid existing limma results',
      StatsAgree ~ 'Existing accession-level limma results agree',
      TRUE ~ 'Group accession limma results differ; significance withheld'
    )
  )


# Save one row per original accession in the pathway for transparent auditing.
carb_group_keys <- c('Dataset','pathway_id','KO','ProteinGroupID')
carb_members <- carb_members %>%
  group_by(across(all_of(carb_group_keys))) %>%
  mutate(
    N_Members_This_Pathway=n_distinct(Protein),
    PathwayMemberAccessions=paste(sort(unique(Protein)),collapse='; ')
  ) %>% ungroup()
readr::write_csv(carb_members,
  file.path(tables_dir,'Carbohydrate_strict_group_accession_audit.csv'))


matched_all <- carb_members %>%
  arrange(Protein) %>%
  group_by(across(all_of(carb_group_keys))) %>%
  slice(1L) %>% ungroup() %>%
  rename(PlotMemberAccession=Protein) %>%
  left_join(strict_group_stats %>%
    select(Dataset,KO,ProteinGroupID,N_Members_Full_Group,MemberAccessions,
           StatsAgree,NeedsReview,StatisticalStatus,Group_log2FC_descriptive),
    by=c('Dataset','KO','ProteinGroupID'),relationship='many-to-one') %>%
  mutate(
    # No new statistical test: use existing p/FDR ONLY if ALL members agree.
    log2FC=if_else(StatsAgree,log2FC,Group_log2FC_descriptive),
    p_value=if_else(StatsAgree,p_value,NA_real_),
    BH_FDR=if_else(StatsAgree,BH_FDR,NA_real_),
    is_significant=is.finite(BH_FDR) & BH_FDR>=0 & BH_FDR<fdr_cutoff &
      is.finite(log2FC) & abs(log2FC)>=abs_logfc_cutoff,
    direction=case_when(
      NeedsReview ~ 'Needs review - no group-level inference',
      is_significant & log2FC>0 ~ 'Up',
      is_significant & log2FC<0 ~ 'Down',
      TRUE ~ 'NS'),
    adj_P_value_plot=if_else(is.finite(BH_FDR) & BH_FDR>=0 & BH_FDR<=1,
                             pmax(BH_FDR,.Machine$double.xmin),NA_real_),
    neg_log10_adj_p=-log10(adj_P_value_plot),
    is_plot_valid=is.finite(log2FC) & is.finite(neg_log10_adj_p),
    KO_plot=if_else(is_significant,
                    paste0(KO,' | ',preferred_function),'Not significant')
  )
if (anyNA(matched_all$StatsAgree)) stop('Some plotted groups lack stats-audit rows.')
if (anyDuplicated(matched_all[carb_group_keys]))
  stop('More than one plotted row remains per Dataset + pathway + KO + exact group.')
if (nrow(matched_all)!=nrow(distinct(carb_members,across(all_of(carb_group_keys)))))
  stop('Unexpected retained strict-group count.')


carb_volcano_df <- matched_all %>% filter(is_plot_valid)
message('Strict original TSV grouping: ',nrow(carb_members),
        ' carbohydrate accession/KO/pathway rows -> ',nrow(matched_all),
        ' group/KO/pathway rows; ',sum(strict_group_stats$NeedsReview),
        ' unique KO groups with unresolved statistics.')
if(!nrow(carb_volcano_df))
  stop('Matched carbohydrate data were found but no finite plot coordinates exist.')

# 7. MAPPING & INPUT QC ------------------------------------------------------
# Compare original vs restored pathway membership (including absent panels).
qc_grid <- tidyr::expand_grid(Dataset=expected,
                              pathway_id=carbohydrate_pathways$pathway_id)

qc_original <- original_carb_map %>%
  inner_join(original_present,by=c('Protein','pathway_id'),
             relationship='many-to-many') %>%
  group_by(Dataset,pathway_id) %>%
  summarise(n_original_protein_KO_memberships=n(),
            n_original_distinct_proteins=n_distinct(Protein),.groups='drop')

qc_restored <- original_carb_map %>%
  inner_join(litter_present,by=c('Protein','pathway_id'),
             relationship='many-to-many') %>%
  group_by(Dataset,pathway_id) %>%
  summarise(n_restored_protein_KO_memberships=n(),
            n_restored_distinct_proteins=n_distinct(Protein),.groups='drop')

qc_global <- joined_before_restriction %>%
  group_by(Dataset,pathway_id) %>%
  summarise(n_limma_results_with_global_KO_mapping=n_distinct(Protein),
            .groups='drop')

qc_matched <- matched_all %>%
  group_by(Dataset,pathway_id) %>%
  summarise(n_matched_distinct_groups=n_distinct(ProteinGroupID),
            n_matched_KO_group_rows=n(),
            n_plottable_groups=n_distinct(ProteinGroupID[is_plot_valid]),
            n_significant_groups=n_distinct(ProteinGroupID[is_significant & is_plot_valid]),
            .groups='drop')

mapping_qc <- qc_grid %>%
  left_join(qc_original,by=c('Dataset','pathway_id')) %>%
  left_join(qc_restored,by=c('Dataset','pathway_id')) %>%
  left_join(qc_global,by=c('Dataset','pathway_id')) %>%
  left_join(qc_matched,by=c('Dataset','pathway_id')) %>%
  mutate(across(starts_with('n_'),~tidyr::replace_na(.x,0L)),
         originally_absent=(n_original_protein_KO_memberships==0L)) %>%
  arrange(match(Dataset,expected),pathway_id)

if(any(mapping_qc$originally_absent & mapping_qc$n_matched_KO_group_rows>0L))
  stop('An originally ABSENT dataset/pathway was reintroduced! Check matching.')


message('Original dataset/pathway membership check (18 panels):')
print(mapping_qc %>%
        select(Dataset,pathway_id,n_original_distinct_proteins,
               n_restored_distinct_proteins,n_matched_distinct_groups,n_significant_groups,
               originally_absent),n=Inf)

# 8. GLOBAL KO COLORS (SAME COLOR ACROSS ALL FIGURES) -------------------------
sig_ko_labels <- carb_volcano_df %>%
  filter(is_significant) %>%
  distinct(KO_plot) %>% arrange(KO_plot) %>% pull(KO_plot)
ko_colors <- make_distinct_colors(length(sig_ko_labels))
names(ko_colors) <- sig_ko_labels
color_values <- c(ko_colors,'Not significant'='grey75')

# 9. PLOT: FIXED X/Y AXES ACROSS 5W, 8W, 12W WITHIN PATHWAY x SEX ----------
plot_one_pathway_sex <- function(target_pathway_id,target_sex,df) {
  plot_df <- df %>%
    filter(pathway_id==target_pathway_id,as.character(Sex)==target_sex)
  if(!nrow(plot_df)) {
    message('No plottable proteins for ',target_pathway_id,' / ',target_sex,
            '; skipping figure (see QC table).')
    return(NULL)
  }

  pathway_nm <- carbohydrate_pathways$pathway_short[
    match(target_pathway_id,carbohydrate_pathways$pathway_id)]
  sig_this_plot <- plot_df %>% filter(is_significant) %>%
    distinct(KO_plot) %>% arrange(KO_plot) %>% pull(KO_plot)

  max_abs_x <- max(abs(plot_df$log2FC),na.rm=TRUE)
  max_y <- max(plot_df$neg_log10_adj_p,na.rm=TRUE)
  x_limit <- max(ceiling(max_abs_x*1.10),ceiling(abs_logfc_cutoff*1.2),1)
  y_limit <- max(ceiling(max_y*1.10),ceiling(-log10(fdr_cutoff)*1.2),1)

  p <- ggplot(plot_df,aes(x=log2FC,y=neg_log10_adj_p)) +
    geom_point(data=plot_df %>% filter(!is_significant),
               aes(color=KO_plot),size=2.1,alpha=0.25) +
    geom_point(data=plot_df %>% filter(is_significant),
               aes(color=KO_plot),size=2.8,alpha=0.95) +
    geom_vline(xintercept=0,linetype='dashed',
               color='grey40',linewidth=0.45) +
    # These light reference lines make the NEW fold-change criterion visible.
    geom_vline(xintercept=c(-abs_logfc_cutoff,abs_logfc_cutoff),
               linetype='dotted',color='grey65',linewidth=0.35) +
    geom_hline(yintercept=-log10(fdr_cutoff),linetype='dotted',
               color='grey40',linewidth=0.45) +
    facet_wrap(~Timepoint,nrow=1,scales='fixed',drop=FALSE) +
    scale_x_continuous(limits=c(-x_limit,x_limit)) +
    scale_y_continuous(limits=c(0,y_limit)) +
    scale_color_manual(values=color_values,
                       breaks=c(sig_this_plot,'Not significant'),
                       name='KO | function') +
    guides(color=guide_legend(ncol=1,
             override.aes=list(size=2.4,alpha=1))) +
    labs(title=paste0(pathway_nm,' (',target_pathway_id,') | ',target_sex),
         x='log2FC (Dosed vs Control)',
         y='-log10 adjusted p-value (BH-FDR)') +
    theme_bw(base_size=10.5) +
    theme(
      plot.title=element_text(face='bold',size=10.5),
      strip.background=element_rect(fill='grey92',color='black',linewidth=0.35),
      strip.text=element_text(face='bold',size=8.5),
      panel.spacing=grid::unit(0.7,'lines'),
      legend.position='right',
      legend.title=element_text(face='bold',size=7.2),
      legend.text=element_text(size=6.2,face='bold',lineheight=0.78),
      legend.key.size=grid::unit(0.24,'cm'),
      legend.spacing.y=grid::unit(0.005,'cm'),
      legend.box.spacing=grid::unit(0.08,'cm'),
      legend.margin=margin(0,0,0,0),
      axis.title=element_text(face='bold',size=9.5),
      axis.text=element_text(size=7.5),
      panel.grid.minor=element_blank(),
      panel.grid.major=element_line(color='grey90',linewidth=0.25),
      plot.margin=margin(t=5,r=3,b=5,l=5)
    )

  stub <- paste0(clean_filename(pathway_nm),'_',
                 clean_filename(target_pathway_id),'_',
                 clean_filename(target_sex),
                 '_volcano_Litter_Blocked_KOfam_first_fixed_axes')

  ggsave(file.path(plots_dir,paste0(stub,'.pdf')),p,
         width=17,height=5.5,units='in',bg='white')
  ggsave(file.path(plots_dir,paste0(stub,'.png')),p,
         width=17,height=5.5,units='in',dpi=plot_dpi,bg='white')
  ggsave(file.path(plots_dir,paste0(stub,'.tiff')),p,
         width=10,height=5,units='in',dpi=plot_dpi,
         device='tiff',compression='lzw',bg='white')

  tibble(pathway_id=target_pathway_id,pathway_name=pathway_nm,
         Sex=target_sex,n_plotted_rows=nrow(plot_df),
         n_significant_points=sum(plot_df$is_significant),
         x_limit=x_limit,y_limit=y_limit,
         pdf_file=file.path(plots_dir,paste0(stub,'.pdf')))
}

plot_keys <- tidyr::expand_grid(pathway_id=carbohydrate_pathways$pathway_id,
                                Sex=sex_levels)
plot_manifest <- purrr::map2_dfr(plot_keys$pathway_id,plot_keys$Sex,
                                 ~plot_one_pathway_sex(.x,.y,carb_volcano_df))
if(!nrow(plot_manifest)) stop('No carbohydrate volcano plots generated.')
readr::write_csv(plot_manifest,file.path(tables_dir,'Plot_manifest.csv'))

# 10. SUPPLEMENTARY TABLES ---------------------------------------------------
parameters <- tibble(
  Parameter=c('Analysis','Statistics input folder','Input files','Input worksheet',
              'Protein fold change','Adjusted p-value','Significance rule',
              'Comparison','Original Protein-KO-pathway membership',
              'NEW restored dataset-specific pathway membership',
              'Original expression source for KO function labels',
              'Function annotation','Pathways','KO plot colors',
               'Axes','Volcano unit','Statistical interpretation',
               'Protein grouping lookup','Grouping criteria','Unresolved groups'),
  Value=c(
    'Protein-level carbohydrate volcano figures from precomputed litter-blocked limma',
    limma_dir,
    paste(paste0(expected,'_Litter_Blocked_Protein_Results.xlsx'),collapse='; '),
    'Protein_Statistics','log2FC','BH_FDR',
    'BH_FDR < 0.05 AND abs(log2FC) >= 1.5',
    'Dosed vs Control within each sex x timepoint',
    original_ko_file,litter_expr_file,original_expr_file,
    'Existing annotation text only, KOfam first; eggNOG second; original map fallback',
    paste(carbohydrate_pathways$pathway_id,collapse='; '),
    'Consistent globally for significant KO | function; grey for non-significant',
    'Fixed across 5W, 8W, 12W within each pathway x sex figure',
    'One point per same KO + exact peptide set + exact original eight-mouse intensity profile group within dataset/pathway',
    'Reads existing litter-blocked protein-level results only when every group member has agreeing statistics; no refitting',
    strict_lookup_file,
    'Same dataset + same KO + exact full peptide set + exact original filtered-TSV intensity and NA pattern across all 8 mice',
    'Excluded from volcano coordinates if group members lack agreeing existing limma results; exported for audit'
  )
)

volcano_plot_source_data <- carb_volcano_df %>%
  select(Dataset,ProteinGroupID,Representative,PlotMemberAccession,GroupSize,
         N_Members_This_Pathway,N_Members_Full_Group,MemberAccessions,
         PathwayMemberAccessions,StatsAgree,NeedsReview,StatisticalStatus,
         KO,pathway_group,pathway_id,pathway_name,
         pathway_short,Sex,Timepoint,log2FC,p_value,BH_FDR,
         adj_P_value_plot,neg_log10_adj_p,is_significant,direction,
         KO_plot,preferred_function,preferred_source,all_functions_seen,
         source_file) %>%
  arrange(pathway_name,Sex,Timepoint,desc(is_significant),KO,ProteinGroupID)

significant_KO_summary <- carb_volcano_df %>%
  filter(is_significant) %>%
  group_by(KO,KO_plot,preferred_function,preferred_source) %>%
  summarise(
    n_significant_group_occurrences=n(),
    member_sets=paste(sort(unique(MemberAccessions)),collapse=" /// "),
    pathways_found=paste(sort(unique(paste0(pathway_name,' (',pathway_id,')'))),
                         collapse=' /// '),
    sexes_found=paste(sort(unique(as.character(Sex))),collapse=' /// '),
    timepoints_found=paste(sort(unique(as.character(Timepoint))),collapse=' /// '),
    directions_found=paste(sort(unique(direction)),collapse=' /// '),
    min_BH_FDR=min(BH_FDR,na.rm=TRUE),
    max_abs_log2FC=max(abs(log2FC),na.rm=TRUE),
    all_functions_seen=paste(sort(unique(na.omit(all_functions_seen))),
                             collapse=' /// '),
    curated_function_name='',.groups='drop') %>%
  arrange(KO)

volcano_summary_by_panel <- matched_all %>%
  group_by(Dataset,pathway_id,pathway_name,pathway_short,Sex,Timepoint) %>%
  summarise(
    n_total_mapped_groups=n_distinct(ProteinGroupID),
    n_plottable_groups=n_distinct(ProteinGroupID[is_plot_valid]),
    n_significant_groups=n_distinct(ProteinGroupID[is_significant & is_plot_valid]),
    n_upregulated_groups=n_distinct(ProteinGroupID[direction=='Up' & is_plot_valid]),
    n_downregulated_groups=n_distinct(ProteinGroupID[direction=='Down' & is_plot_valid]),
    n_significant_KOs=n_distinct(KO[is_significant & is_plot_valid]),
    n_groups_needing_review=n_distinct(ProteinGroupID[NeedsReview]),
    .groups='drop') %>%
  arrange(pathway_name,Sex,Timepoint)

readr::write_csv(parameters,file.path(tables_dir,'carbohydrate_volcano_analysis_parameters.csv'))
readr::write_csv(volcano_plot_source_data,
                 file.path(tables_dir,'carbohydrate_volcano_plot_source_data.csv'))
readr::write_csv(significant_KO_summary,
                 file.path(tables_dir,'carbohydrate_volcano_significant_KO_summary.csv'))
readr::write_csv(volcano_summary_by_panel,
                 file.path(tables_dir,'carbohydrate_volcano_summary_by_panel.csv'))


writexl::write_xlsx(list(
  README_analysis_parameters=parameters,
  KO_preferred_function_labels=ko_preferred_function,
  volcano_plot_source_data=volcano_plot_source_data,
  significant_KO_summary=significant_KO_summary,
  volcano_summary_by_panel=volcano_summary_by_panel,
  QC_dataset_pathway_mapping=mapping_qc,
  Plot_manifest=plot_manifest
),file.path(tables_dir,'Carbohydrate_Volcano_Litter_Blocked_Supplementary.xlsx'))


