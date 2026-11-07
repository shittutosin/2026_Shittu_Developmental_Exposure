

required_packages <- c('readxl','readr','dplyr','tidyr','tibble','stringr',
                       'purrr','limma','statmod','ggplot2','writexl')
missing_packages <- required_packages[!vapply(
  required_packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing_packages)) stop(
  'Install missing R packages before running:\n',
  paste(missing_packages,collapse=', '))
suppressPackageStartupMessages({
  library(readxl); library(readr); library(dplyr); library(tidyr)
  library(tibble); library(stringr); library(purrr); library(ggplot2)
})

# 0. SETTINGS ----------------------------------------------------------------
root <- '/file_root'
limma_dir <- file.path(root,'Litter_analysis','Litter_Blocked_Limma_files')
hadeg_file <- file.path(root,'HADEG',
                       '02_HADEG_best_hit_per_protein_for_manual_review_STRICT.csv')
out_dir <- file.path(root,'Litter_analysis',
                     'PCAF_tRENDLINE_Litter_Blocked_Limma')

include_unreviewed <- TRUE

expected_datasets <- c('5W_M','8W_M','12W_M','5W_F','8W_F','12W_F')
time_levels <- c('5W','8W','12W')
sex_levels <- c('Males','Females')
fdr_cutoff <- 0.05
abs_logfc_cutoff <- 1.5          
min_litters_per_treatment <- 2L  
min_features_for_ebayes <- 3L    
plot_dpi <- 450

for(d in c(out_dir,file.path(out_dir,'tables'),file.path(out_dir,'plots')))
  dir.create(d,recursive=TRUE,showWarnings=FALSE)
table_dir <- file.path(out_dir,'tables')
plot_dir <- file.path(out_dir,'plots')
write_tab <- function(x,name) readr::write_csv(
  x,file.path(table_dir,paste0(name,'.csv')),na='')

clean_chr <- function(z) stringr::str_squish(as.character(z))
valid_chr <- function(z) !is.na(z) & nzchar(clean_chr(z)) &
  !(tolower(clean_chr(z)) %in% c('na','none','-','unknown','unknown gene'))
normalize_treatment <- function(z) {
  z <- tolower(clean_chr(z))
  dplyr::case_when(z=='control'~'Control',z=='dosed'~'Dosed',TRUE~NA_character_)
}
clean_description <- function(z) {
  z <- clean_chr(z)
  z <- stringr::str_replace_all(z,'_', ' ')
  z <- stringr::str_replace(z,'\\s+Bacillus subtilis\\s+RB14$','')
  z <- stringr::str_replace(z,'\\s+Bacillus subtilis$','')
  z <- stringr::str_replace(z,'^[A-Z0-9]+_[A-Z0-9]+\\s+','')
  stringr::str_squish(z)
}


log2_sum_linear <- function(z) {
  if(!length(z) || any(!is.finite(z))) return(NA_real_)
  m <- max(z)
  m + log2(sum(2^(z-m)))
}


if(!file.exists(hadeg_file)) {
  alternative <- file.path(root,'HADEG','results_strict',
    'HADEG_inspection_files_large_strict',
    '02_HADEG_best_hit_per_protein_for_manual_review_STRICT.csv')
  if(file.exists(alternative)) {
    hadeg_file <- alternative
    message('Using original strict HADEG annotation at: ',hadeg_file)
  } else stop('Missing strict HADEG annotation. Checked:\n',hadeg_file,
              '\nand\n',alternative)
}

anno_raw <- readr::read_csv(hadeg_file,show_col_types=FALSE)
names(anno_raw) <- clean_chr(names(anno_raw))
required_anno <- c('Protein','HADEG_gene','HADEG_description')
if(!all(required_anno %in% names(anno_raw)))
  stop('HADEG annotation missing: ',
       paste(setdiff(required_anno,names(anno_raw)),collapse=', '))
if(!'Include_in_heatmap' %in% names(anno_raw)) anno_raw$Include_in_heatmap <- NA_character_
if(!'Manual_Broad_Class' %in% names(anno_raw)) anno_raw$Manual_Broad_Class <- NA_character_
if(!'Manual_Subclass' %in% names(anno_raw)) anno_raw$Manual_Subclass <- NA_character_

approved_terms <- c('yes','y','include','included','keep','true','1')
anno <- anno_raw %>%
  mutate(Protein=clean_chr(Protein),
         HADEG_gene=clean_chr(HADEG_gene),
         HADEG_description=clean_description(HADEG_description),
         HADEG_gene_key=tolower(HADEG_gene),
         Review_raw=clean_chr(Include_in_heatmap),
         Reviewed=tolower(Review_raw) %in% approved_terms,
         HADEG_functional_class=paste0(HADEG_gene,' | ',HADEG_description)) %>%
  filter(valid_chr(Protein))


annotation_duplicates <- anno %>%
  group_by(Protein) %>%
  summarise(n_annotation_rows=n(),
            n_distinct_gene=n_distinct(HADEG_gene_key),
            n_distinct_description=n_distinct(HADEG_description),
            gene_values=paste(unique(HADEG_gene),collapse=' | '),
            .groups='drop') %>%
  filter(n_annotation_rows>1)
write_tab(annotation_duplicates,'QC_duplicate_HADEG_annotations')
if(nrow(annotation_duplicates)) warning(
  nrow(annotation_duplicates),' duplicate HADEG protein records detected; ',
  'first best hit retained as in original volcano. Inspect QC CSV.')

anno <- anno %>% distinct(Protein,.keep_all=TRUE)
write_tab(anno,'HADEG_original_best_hit_annotation')
if(!include_unreviewed) anno <- anno %>% filter(Reviewed)
anno <- anno %>% filter(valid_chr(HADEG_gene_key))
if(!nrow(anno)) stop('No HADEG proteins have usable gene labels after annotation filtering.')


gene_description_qc <- anno %>%
  group_by(HADEG_gene_key) %>%
  summarise(n_proteins=n_distinct(Protein),
            n_descriptions=n_distinct(HADEG_description),
            descriptions=paste(sort(unique(HADEG_description)),collapse=' | '),
            n_unreviewed_proteins=n_distinct(Protein[!Reviewed]),
            .groups='drop')
write_tab(gene_description_qc,'QC_HADEG_gene_multiple_descriptions')
if(any(gene_description_qc$n_descriptions>1)) warning(
  'Some HADEG genes have multiple descriptions. Check QC_HADEG_gene_multiple_descriptions.csv; ',
  'the separate gene+description sensitivity analysis avoids merging them.')


read_one_dataset <- function(id) {
  fp <- file.path(limma_dir,paste0(id,'_Litter_Blocked_Protein_Results.xlsx'))
  if(!file.exists(fp)) stop('Missing NEW litter-blocked workbook: ',fp)
  sheets <- readxl::excel_sheets(fp)
  required_sheets <- c('Individual_Mouse_Log2','Mouse_Litter_Mapping')
  if(!all(required_sheets %in% sheets)) stop(
    id,': missing sheets: ',paste(setdiff(required_sheets,sheets),collapse=', '))

  expr <- readxl::read_excel(fp,sheet='Individual_Mouse_Log2')
  meta <- readxl::read_excel(fp,sheet='Mouse_Litter_Mapping')
  needed_meta <- c('MouseID','OriginalColumn','LitterID','Treatment')
  if(!'Protein' %in% names(expr)) stop(id,': expression missing Protein column')
  if(!all(needed_meta %in% names(meta))) stop(
    id,': metadata missing: ',paste(setdiff(needed_meta,names(meta)),collapse=', '))
  if(nrow(meta)!=8L) stop(id,': expected 8 individual mice, found ',nrow(meta))

  meta <- meta %>% mutate(
    Dataset=id,
    MouseID=clean_chr(MouseID),
    OriginalColumn=clean_chr(OriginalColumn),
    LitterID=clean_chr(LitterID),
    Treatment=normalize_treatment(Treatment),
    Sex=if_else(grepl('_M$',id),'Males','Females'),
    Timepoint=sub('_.*$','',id))

  if(anyDuplicated(meta$MouseID) || anyDuplicated(meta$OriginalColumn))
    stop(id,': duplicate mouse IDs or original sample columns')
  if(any(!valid_chr(meta$MouseID)) || any(!valid_chr(meta$LitterID)) ||
     any(is.na(meta$Treatment)))
    stop(id,': invalid mouse ID, treatment, or litter ID')
  if(!setequal(meta$OriginalColumn,
               c(paste0('Control_',1:4),paste0('Dosed_',1:4))))
    stop(id,': original columns do not correspond to Control_1..4 and Dosed_1..4')
  if(any(!startsWith(meta$OriginalColumn[meta$Treatment=='Control'],'Control_')) ||
     any(!startsWith(meta$OriginalColumn[meta$Treatment=='Dosed'],'Dosed_')))
    stop(id,': treatment assignments do not match OriginalColumn')
  if(!all(meta$MouseID %in% names(expr)))
    stop(id,': expression workbook missing mice: ',
         paste(setdiff(meta$MouseID,names(expr)),collapse=', '))

  # In a dam-dosed experiment, each litter belongs to one treatment.
  if(any(meta %>% count(LitterID,Treatment) %>% count(LitterID) %>%
         pull(n) > 1L)) stop(
    id,': at least one litter ID occurs in both treatments. ',
    'Check litter IDs before fitting models.')

  expr <- expr %>% mutate(Protein=clean_chr(Protein))
  if(anyDuplicated(expr$Protein)) stop(id,': duplicated protein accession in expression')
  expr_used <- expr %>% inner_join(
    anno %>% select(Protein,HADEG_gene,HADEG_gene_key,
                    HADEG_description,HADEG_functional_class,Reviewed),
    by='Protein',relationship='many-to-one')

  message('Loaded ',id,': ',nrow(expr_used),
          ' HADEG-matched proteins; ',n_distinct(meta$LitterID),' litters')
  if(!nrow(expr_used)) return(list(
    data=tibble(),meta=meta,
    qc=tibble(Dataset=id,n_mice=nrow(meta),n_litters=n_distinct(meta$LitterID),
              n_HADEG_proteins=0L,n_HADEG_genes=0L,Workbook=fp)))

  selected <- expr_used %>%
    select(Protein,HADEG_gene,HADEG_gene_key,HADEG_description,
           HADEG_functional_class,Reviewed,all_of(meta$MouseID)) %>%
    pivot_longer(cols=all_of(meta$MouseID),names_to='MouseID',
                 values_to='log2_intensity') %>%
    mutate(log2_intensity=suppressWarnings(as.numeric(log2_intensity))) %>%
    left_join(meta %>% select(Dataset,MouseID,OriginalColumn,LitterID,
                              Treatment,Sex,Timepoint),
              by='MouseID',relationship='many-to-one')

  list(
    data=selected,
    meta=meta,
    qc=tibble(Dataset=id,n_mice=nrow(meta),n_litters=n_distinct(meta$LitterID),
              n_HADEG_proteins=n_distinct(expr_used$Protein),
              n_HADEG_genes=n_distinct(expr_used$HADEG_gene_key),
              Workbook=fp))
}

loaded <- purrr::map(expected_datasets,read_one_dataset)
all_source <- purrr::map_dfr(loaded,'data')
all_meta <- purrr::map_dfr(loaded,'meta')
source_qc <- purrr::map_dfr(loaded,'qc')
write_tab(all_meta,'All_Mouse_Litter_Metadata')
write_tab(source_qc,'QC_HADEG_matching_by_dataset')
if(!nrow(all_source)) stop('No annotated HADEG proteins matched any mouse-expression workbook.')
if(anyDuplicated(all_source[c('Dataset','MouseID','Protein')]))
  stop('Duplicate protein-by-mouse observations after the annotation join')
write_tab(all_source,'HADEG_Matched_Protein_Mouse_Intensities')

make_mouse_abundance <- function(mode) {
  if(mode=='gene') {
    dat <- all_source %>% mutate(
      GroupID=HADEG_gene_key,
      GroupLabel=HADEG_gene)
  } else if(mode=='gene_description') {
    dat <- all_source %>% mutate(
      GroupID=paste(HADEG_gene_key,HADEG_description,sep=' || '),
      GroupLabel=paste(HADEG_gene,HADEG_description,sep=' | '))
  } else stop('Unknown aggregation mode: ',mode)

  members <- dat %>%
    distinct(Dataset,Sex,Timepoint,GroupID,GroupLabel,
             Protein,HADEG_gene,HADEG_description,Reviewed) %>%
    arrange(Dataset,GroupLabel,Protein)

  mouse <- dat %>%
    group_by(Dataset,Sex,Timepoint,GroupID,GroupLabel,
             MouseID,OriginalColumn,Treatment,LitterID) %>%
    summarise(
      n_mapped_proteins=n_distinct(Protein),
      n_finite_proteins=sum(is.finite(log2_intensity)),
      any_nonfinite=any(!is.finite(log2_intensity)),
      log2_sum_abundance=log2_sum_linear(log2_intensity),
      .groups='drop')

  completeness <- mouse %>%
    group_by(Dataset,Sex,Timepoint,GroupID,GroupLabel) %>%
    summarise(n_mouse_rows=n_distinct(MouseID),
              n_finite_mouse_scores=sum(is.finite(log2_sum_abundance)),
              min_proteins_in_mouse=min(n_mapped_proteins),
              max_proteins_in_mouse=max(n_mapped_proteins),
              .groups='drop') %>%
    left_join(all_meta %>% count(Dataset,name='expected_n_mice'),by='Dataset') %>%
    mutate(Complete=n_mouse_rows==expected_n_mice &
             n_finite_mouse_scores==expected_n_mice)

  list(mouse=mouse,members=members,completeness=completeness)
}

# 4. RUN A LITTER-BLOCKED LIMMA MODEL PER SEX x TIMEPOINT --------------------
fit_dataset <- function(id,mouse,completeness,mode) {
  meta <- all_meta %>% filter(Dataset==id) %>%
    arrange(match(OriginalColumn,
                  c(paste0('Control_',1:4),paste0('Dosed_',1:4)))) %>%
    mutate(Treatment=factor(Treatment,levels=c('Control','Dosed')))
  if(!nrow(meta)) stop('Missing metadata for ',id)
  available <- completeness %>% filter(Dataset==id)
  if(!nrow(available)) return(tibble())

  # Retain descriptions/IDs of all features, including missing scores.
  status <- available %>%
    select(Dataset,Sex,Timepoint,GroupID,GroupLabel,Complete,
           n_mouse_rows,n_finite_mouse_scores,
           min_proteins_in_mouse,max_proteins_in_mouse)
  valid <- status %>% filter(Complete)
  if(!nrow(valid)) {
    return(status %>% mutate(Model_status='INCOMPLETE_MOUSE_SCORES',
                             log2FC=NA_real_,AveExpr=NA_real_,t=NA_real_,
                             p_value=NA_real_,BH_FDR=NA_real_,B=NA_real_,
                             Litter_correlation=NA_real_,
                             Correlation_method='not_tested'))
  }

  dat <- mouse %>% filter(Dataset==id,GroupID %in% valid$GroupID) %>%
    select(GroupID,MouseID,log2_sum_abundance)
  y_wide <- dat %>% pivot_wider(names_from=MouseID,
                                values_from=log2_sum_abundance)
  y_wide <- y_wide[match(valid$GroupID,y_wide$GroupID),,drop=FALSE]
  if(anyNA(y_wide$GroupID)) stop(id,': unmatched HADEG group ID in expression matrix')
  mat <- as.matrix(y_wide[,meta$MouseID,drop=FALSE]);
  storage.mode(mat) <- 'double'
  rownames(mat) <- y_wide$GroupID
  if(any(!is.finite(mat))) stop(id,': non-finite HADEG matrix despite completeness QC')

  litters_per_treatment <- meta %>%
    group_by(Treatment) %>%
    summarise(n_litters=n_distinct(LitterID),n_mice=n(),.groups='drop')
  enough_litters <- all(litters_per_treatment$n_litters>=min_litters_per_treatment)
  status$Independent_Litters_Control <- litters_per_treatment$n_litters[
    as.character(litters_per_treatment$Treatment)=='Control']
  status$Independent_Litters_Dosed <- litters_per_treatment$n_litters[
    as.character(litters_per_treatment$Treatment)=='Dosed']

  reason <- NULL
  if(!enough_litters) reason <- 'INSUFFICIENT_INDEPENDENT_LITTERS'
  if(nrow(mat)<min_features_for_ebayes) reason <- 'TOO_FEW_HADEG_FEATURES_FOR_EBAYES'
  if(!is.null(reason)) {
    warning(id,' / ',mode,': ',reason,'; no significance was inferred.')
    return(status %>% mutate(
      Model_status=if_else(Complete,reason,'INCOMPLETE_MOUSE_SCORES'),
      log2FC=NA_real_,AveExpr=NA_real_,t=NA_real_,p_value=NA_real_,
      BH_FDR=NA_real_,B=NA_real_,Litter_correlation=NA_real_,
      Correlation_method='not_tested'))
  }

  design <- model.matrix(~Treatment,data=meta)
  if(qr(design)$rank<ncol(design)) stop(id,': non-estimable treatment design')
  has_repeated_litters <- anyDuplicated(meta$LitterID)>0
  correlation <- NA_real_
  method <- 'unique_litters_no_within_litter_repeats'
  if(has_repeated_litters) {
    dc <- tryCatch(
      limma::duplicateCorrelation(mat,design=design,block=meta$LitterID),
      error=function(e) e)
    if(inherits(dc,'error') ||
       !is.finite(if(inherits(dc,'error')) NA_real_ else dc$consensus.correlation)) {
      warning(id,' / ',mode,': duplicateCorrelation failed. No unblocked fallback.')
      return(status %>% mutate(
        Model_status=if_else(Complete,'LITTER_CORRELATION_FAILED',
                             'INCOMPLETE_MOUSE_SCORES'),
        log2FC=NA_real_,AveExpr=NA_real_,t=NA_real_,p_value=NA_real_,
        BH_FDR=NA_real_,B=NA_real_,Litter_correlation=NA_real_,
        Correlation_method='failed'))
    }
    correlation <- as.numeric(dc$consensus.correlation)
    method <- 'limma_duplicateCorrelation_litter_block'
    fit0 <- tryCatch(limma::lmFit(mat,design=design,
                                 block=meta$LitterID,correlation=correlation),
                     error=function(e) e)
  } else {
    fit0 <- tryCatch(limma::lmFit(mat,design=design),error=function(e)e)
  }

  if(inherits(fit0,'error')) {
    warning(id,' / ',mode,': limma lmFit failed: ',conditionMessage(fit0))
    return(status %>% mutate(
      Model_status=if_else(Complete,'LIMMA_FIT_FAILED','INCOMPLETE_MOUSE_SCORES'),
      log2FC=NA_real_,AveExpr=NA_real_,t=NA_real_,p_value=NA_real_,
      BH_FDR=NA_real_,B=NA_real_,Litter_correlation=correlation,
      Correlation_method=method))
  }

  fit <- tryCatch(limma::eBayes(fit0),error=function(e)e)
  if(inherits(fit,'error')) {
    warning(id,' / ',mode,': limma eBayes failed: ',conditionMessage(fit))
    return(status %>% mutate(
      Model_status=if_else(Complete,'EBAYES_FAILED','INCOMPLETE_MOUSE_SCORES'),
      log2FC=NA_real_,AveExpr=NA_real_,t=NA_real_,p_value=NA_real_,
      BH_FDR=NA_real_,B=NA_real_,Litter_correlation=correlation,
      Correlation_method=method))
  }

  tt <- limma::topTable(fit,coef='TreatmentDosed',number=Inf,
                        adjust.method='BH',sort.by='none')
  tt$GroupID <- rownames(tt)
  infer <- tibble(GroupID=as.character(tt$GroupID),
                  log2FC=as.numeric(tt$logFC),
                  AveExpr=as.numeric(tt$AveExpr),
                  t=as.numeric(tt$t),
                  p_value=as.numeric(tt$P.Value),
                  BH_FDR=as.numeric(tt$adj.P.Val),
                  B=if('B' %in% names(tt)) as.numeric(tt$B)
                    else rep(NA_real_,nrow(tt)),
                  Model_status='TESTED',
                  Litter_correlation=correlation,
                  Correlation_method=method)
  status %>% left_join(infer,by='GroupID') %>%
    mutate(Model_status=if_else(!Complete,'INCOMPLETE_MOUSE_SCORES',Model_status),
           Correlation_method=if_else(!Complete,'not_tested',Correlation_method))
}

# 5. ANALYSIS FOR EACH GROUPING MODE -----------------------------------------
run_mode <- function(mode) {
  message('\n============== Analysis mode: ',mode,' ==============')
  s <- make_mouse_abundance(mode)
  prefix <- if(mode=='gene') 'HADEG_GENE' else 'HADEG_GENE_DESCRIPTION'
  write_tab(s$members,paste0(prefix,'_Original_Member_Proteins'))
  write_tab(s$mouse,paste0(prefix,'_Mouse_Pooled_Abundance'))
  write_tab(s$completeness,paste0(prefix,'_Completeness_QC'))

  results <- purrr::map_dfr(
    expected_datasets,~fit_dataset(.x,s$mouse,s$completeness,mode))
  if(!nrow(results)) {
    warning('No HADEG features were available for ',mode)
    return(list(mode=mode,results=results,mouse=s$mouse,members=s$members))
  }

  group_summary <- s$mouse %>%
    group_by(Dataset,GroupID,Treatment) %>%
    summarise(n_mice=n_distinct(MouseID),n_litters=n_distinct(LitterID),
              mean_log2_abundance=mean(log2_sum_abundance,na.rm=TRUE),
              median_log2_abundance=median(log2_sum_abundance,na.rm=TRUE),
              .groups='drop') %>%
    pivot_wider(names_from=Treatment,
                values_from=c(n_mice,n_litters,
                              mean_log2_abundance,median_log2_abundance)) %>%
    mutate(
      Descriptive_mean_Dosed_minus_Control=
        mean_log2_abundance_Dosed-mean_log2_abundance_Control,
      Descriptive_median_Dosed_minus_Control=
        median_log2_abundance_Dosed-median_log2_abundance_Control)

  results <- results %>%
    left_join(group_summary,by=c('Dataset','GroupID')) %>%
    mutate(
      Significant_FDR=is.finite(BH_FDR) & BH_FDR<fdr_cutoff,
      Significant_FDR_and_FC=Significant_FDR &
        is.finite(log2FC) & abs(log2FC)>=abs_logfc_cutoff,
      Direction=case_when(
        Significant_FDR & log2FC>0~'Higher in Dosed',
        Significant_FDR & log2FC<0~'Lower in Dosed',
        !is.na(BH_FDR)~'Not significant',
        TRUE~'Not tested'),
      Plot_status=case_when(
        Significant_FDR_and_FC & log2FC>0~'Higher in Dosed',
        Significant_FDR_and_FC & log2FC<0~'Lower in Dosed',
        !is.na(BH_FDR)~'Not significant / below FC threshold',
        TRUE~'Not tested')) %>%
    arrange(Dataset,BH_FDR,GroupID)

  write_tab(results,paste0(prefix,'_Litter_Blocked_Limma_Results'))
  write_tab(filter(results,Significant_FDR),
            paste0(prefix,'_FDR_Significant'))
  write_tab(filter(results,Significant_FDR_and_FC),
            paste0(prefix,'_FDR_and_FC_Significant'))

  
  wr <- list(
    Limma_results=results,
    Member_proteins=s$members,
    Per_mouse_abundance=s$mouse,
    Completeness_QC=s$completeness,
    Source_dataset_QC=source_qc,
    Gene_description_QC=gene_description_qc)
  writexl::write_xlsx(wr,file.path(table_dir,paste0(prefix,'_Supplementary.xlsx')))

  list(mode=mode,results=results,mouse=s$mouse,members=s$members)
}

gene <- run_mode('gene')
gene_description <- run_mode('gene_description')


member_counts <- gene$members %>%
  group_by(Dataset,Sex,Timepoint,GroupID) %>%
  summarise(n_member_proteins=n_distinct(Protein),
            n_reviewed_proteins=n_distinct(Protein[Reviewed]),
            n_unreviewed_proteins=n_distinct(Protein[!Reviewed]),
            n_distinct_HADEG_descriptions=n_distinct(HADEG_description),
            descriptions=paste(sort(unique(HADEG_description)),collapse=' | '),
            .groups='drop')

mouse_score_counts <- gene$mouse %>%
  group_by(Dataset,GroupID,Treatment) %>%
  summarise(n_mouse_scores=sum(is.finite(log2_sum_abundance)),
            n_mice_present=n_distinct(MouseID),
            .groups='drop') %>%
  tidyr::pivot_wider(names_from=Treatment,
                     values_from=c(n_mouse_scores,n_mice_present))


litter_counts <- all_meta %>%
  group_by(Dataset,Treatment) %>%
  summarise(Metadata_n_litters=n_distinct(LitterID),
            Metadata_n_mice=n_distinct(MouseID),.groups='drop') %>%
  tidyr::pivot_wider(names_from=Treatment,
                     values_from=c(Metadata_n_litters,Metadata_n_mice))


expected_litter_fields <- c('Dataset','Metadata_n_litters_Control',
                             'Metadata_n_litters_Dosed')
if(!all(expected_litter_fields %in% names(litter_counts)))
  stop('Litter count audit is missing columns: ',
       paste(setdiff(expected_litter_fields,names(litter_counts)),collapse=', '),
       '. Inspect Treatment values in all_meta.')

eligibility_master <- gene$results %>%
  left_join(member_counts,by=c('Dataset','Sex','Timepoint','GroupID')) %>%
  left_join(mouse_score_counts,by=c('Dataset','GroupID')) %>%
  left_join(litter_counts,by='Dataset') %>%
  mutate(
    Sufficient_litters_each_group=is.finite(Metadata_n_litters_Control) &
      is.finite(Metadata_n_litters_Dosed) &
      Metadata_n_litters_Control>=min_litters_per_treatment &
      Metadata_n_litters_Dosed>=min_litters_per_treatment,
    Eligibility_status=case_when(
      !Complete ~ 'INCOMPLETE_MOUSE_ABUNDANCE',
      !Sufficient_litters_each_group ~ 'INSUFFICIENT_LITTER_REPLICATION',
      Model_status=='TESTED' ~ 'TESTED',
      TRUE ~ coalesce(Model_status,'NOT_TESTED')
    ),
    Check_multiple_descriptions=n_distinct_HADEG_descriptions>1,
    Check_unreviewed_proteins=n_unreviewed_proteins>0
  ) %>%
  arrange(Dataset,Eligibility_status,GroupID)

write_tab(eligibility_master,'HADEG_ALL_Genes_Eligibility_By_Sex_Timepoint')
write_tab(eligibility_master %>% filter(Eligibility_status=='TESTED'),
          'HADEG_Eligible_Tested_Gene_List')

eligibility_summary <- eligibility_master %>%
  group_by(Dataset,Sex,Timepoint) %>%
  summarise(n_HADEG_genes_mapped=n_distinct(GroupID),
            n_HADEG_genes_tested=sum(Eligibility_status=='TESTED'),
            n_not_tested=sum(Eligibility_status!='TESTED'),
            n_with_multiple_descriptions=sum(Check_multiple_descriptions,na.rm=TRUE),
            n_with_unreviewed_proteins=sum(Check_unreviewed_proteins,na.rm=TRUE),
            .groups='drop') %>%
  arrange(Sex,Timepoint)
write_tab(eligibility_summary,'QC_HADEG_Function_Eligibility_Counts')
print(eligibility_summary,n=Inf)


old_style_target_genes <- c('pcaf','pcai')
old_style_width <- 8.0
old_style_height <- 4.8
old_style_dpi <- 600

old_style_gene_result <- gene$results
if (!all(c('GroupID','Sex','Timepoint','Model_status','log2FC','BH_FDR') %in%
         names(old_style_gene_result))) {
  stop('Cannot draw old-style plots: required fields missing from gene$results')
}

old_style_values <- old_style_gene_result %>%
  filter(tolower(as.character(GroupID)) %in% old_style_target_genes,
         Model_status == 'TESTED',is.finite(log2FC)) %>%
  pull(log2FC)
old_style_max <- if(length(old_style_values)) max(abs(old_style_values)) else 1
if(!is.finite(old_style_max) || old_style_max <= 0) old_style_max <- 1
old_style_y_limit <- max(1,ceiling(old_style_max * 1.35))
old_style_offset <- old_style_y_limit * 0.065

old_style_stars <- function(q) {
  dplyr::case_when(
    !is.finite(q) ~ '',
    q < 0.0001 ~ '****',
    q < 0.001  ~ '***',
    q < 0.01   ~ '**',
    q < 0.05   ~ '*',
    TRUE ~ ''
  )
}

old_style_colors <- c('Higher in Control'='steelblue',
                      'Higher in Dosed'='firebrick',
                      'Not significant'='grey65')

make_old_style_hadeg_plot <- function(target_gene, result_df) {
  current <- result_df %>%
    filter(tolower(as.character(GroupID)) == tolower(target_gene)) %>%
    mutate(Sex=as.character(Sex),Timepoint=as.character(Timepoint))
  if(anyDuplicated(current[c('Sex','Timepoint')]))
    stop('Duplicate sex/timepoint rows for ',target_gene)

  # Scaffold ensures all six original positions are always represented.
  panel <- tidyr::expand_grid(Sex=sex_levels,Timepoint=time_levels) %>%
    left_join(current,by=c('Sex','Timepoint')) %>%
    mutate(
      Timepoint=factor(Timepoint,levels=time_levels),
      Sex=factor(Sex,levels=sex_levels),
      Tested=(!is.na(Model_status) & Model_status=='TESTED' &
              is.finite(log2FC) & is.finite(BH_FDR)),
      Effect=if_else(Tested,as.numeric(log2FC),NA_real_),
      Plot_group=case_when(
        Tested & BH_FDR < fdr_cutoff & Effect > 0 ~ 'Higher in Dosed',
        Tested & BH_FDR < fdr_cutoff & Effect < 0 ~ 'Higher in Control',
        Tested ~ 'Not significant',
        TRUE ~ NA_character_
      ),
      Plot_group=factor(Plot_group,levels=names(old_style_colors)),
      Missing_label=case_when(
        is.na(GroupID) ~ 'Not detected',
        !Tested ~ 'Not tested',
        TRUE ~ NA_character_
      ),
      Star=old_style_stars(BH_FDR),
      Star=if_else(Tested & BH_FDR < fdr_cutoff,Star,''),
      Star_y=case_when(
        is.finite(Effect) & Effect>=0 ~ Effect+old_style_offset,
        is.finite(Effect) & Effect<0 ~ Effect-old_style_offset,
        TRUE ~ NA_real_
      )
    )

  
  line_segments <- panel %>%
    mutate(Time_index=match(as.character(Timepoint),time_levels)) %>%
    arrange(Sex,Time_index) %>%
    group_by(Sex) %>%
    mutate(Next_index=lead(Time_index),Next_effect=lead(Effect),
           Next_timepoint=factor(time_levels[Next_index],levels=time_levels)) %>%
    ungroup() %>%
    filter(is.finite(Effect),is.finite(Next_effect),
           !is.na(Next_index),Next_index==Time_index+1L)

  p <- ggplot(panel,aes(x=Timepoint,y=Effect)) +
    geom_hline(yintercept=0,linetype='dashed',linewidth=0.8,color='grey45') +
    geom_segment(data=line_segments,
                 aes(x=Timepoint,xend=Next_timepoint,y=Effect,yend=Next_effect),
                 inherit.aes=FALSE,color='grey35',linewidth=1.4) +
    geom_point(data=panel %>% filter(Tested),
               aes(color=Plot_group),size=5,show.legend=TRUE) +
    geom_text(data=panel %>% filter(Tested,Star!=''),
              aes(y=Star_y,label=Star,color=Plot_group),
              size=5,fontface='italic',show.legend=FALSE) +
    geom_text(data=panel %>% filter(!Tested),
              aes(x=Timepoint,y=0,label=Missing_label),inherit.aes=FALSE,
              size=3.1,fontface='bold.italic',color='grey50') +
    facet_wrap(~Sex,nrow=1,drop=FALSE) +
    scale_x_discrete(limits=time_levels,drop=FALSE) +
    scale_y_continuous(limits=c(-old_style_y_limit,old_style_y_limit),
                       breaks=scales::pretty_breaks(n=5),
                       expand=expansion(mult=c(0.04,0.04))) +
    scale_color_manual(name=NULL,values=old_style_colors,
                       limits=names(old_style_colors),drop=FALSE) +
    guides(color=guide_legend(nrow=1,byrow=TRUE,
                              override.aes=list(size=5,alpha=1))) +
    labs(title=NULL,x='Timepoint',
         y='Log2FC of summed abundance\n(Dosed - Control)') +
    coord_cartesian(clip='off') +
    theme_classic(base_size=17) +
    theme(
      plot.title=element_blank(),
      strip.background=element_rect(fill='grey90',color='black',linewidth=0.8),
      strip.text=element_text(size=17,face='bold',color='black',
                              margin=margin(7,7,7,7)),
      axis.title.x=element_text(size=17,face='bold',color='black',margin=margin(t=10)),
      axis.title.y=element_text(size=12,face='bold',color='black',margin=margin(r=10)),
      axis.text.x=element_text(size=15,face='bold',color='black'),
      axis.text.y=element_text(size=14,face='bold',color='black'),
      axis.ticks=element_line(linewidth=0.8,color='black'),
      panel.border=element_rect(color='black',fill=NA,linewidth=0.8),
      panel.grid.major.y=element_line(color='grey92',linewidth=0.5),
      panel.grid.major.x=element_blank(),
      panel.grid.minor=element_blank(),
      panel.spacing=grid::unit(1.2,'lines'),
      legend.position='top',legend.direction='horizontal',
      legend.text=element_text(face='bold',size=11,color='black'),
      legend.key.width=grid::unit(0.7,'cm'),
      legend.key.height=grid::unit(0.7,'cm'),
      legend.spacing.x=grid::unit(0.25,'cm'),
      plot.margin=margin(t=15,r=18,b=15,l=18)
    )

  display_gene <- if(tolower(target_gene)=='pcaf') 'pcaF' else if(tolower(target_gene)=='pcai') 'pcaI' else target_gene
  file_stub <- paste0(display_gene,'_HADEG_gene_level_litter_limma_old_style')
  ggsave(file.path(plot_dir,paste0(file_stub,'.pdf')),p,
         width=old_style_width,height=old_style_height,bg='white')
  ggsave(file.path(plot_dir,paste0(file_stub,'.tiff')),p,
         width=old_style_width,height=old_style_height,
         dpi=old_style_dpi,compression='lzw',bg='white')
  ggsave(file.path(plot_dir,paste0(file_stub,'.png')),p,
         width=old_style_width,height=old_style_height,
         dpi=300,bg='white')
  write_tab(panel %>% mutate(across(where(is.factor),as.character)),
            paste0(display_gene,'_Old_Style_Plot_Source_Data'))
  message('Saved old-style HADEG plot: ',target_gene)
  p
}

old_style_plots <- setNames(lapply(old_style_target_genes,
  make_old_style_hadeg_plot,result_df=old_style_gene_result),
  old_style_target_genes)

=
parameters <- tibble(
  Parameter=c('Analysis name','Annotation file','Include unreviewed hits',
              'Expression input','Metadata input','Unit tested',
              'Aggregation','Comparison','Sex','Timepoint',
              'Statistical model','Blocking method','Independent litter rule',
              'BH correction family','Dual-threshold plot rule',
              'Gene+description analysis','Unmatched proteins',
              'Old-style plot significance criterion','Old-style plot y-axis'),
  Value=c(
    'HADEG gene-level / HADEG gene+description functional-level (NOT KEGG KO)',
    hadeg_file,
    as.character(include_unreviewed),
    'Individual_Mouse_Log2 sheet of the six new litter-blocked workbooks',
    'Mouse_Litter_Mapping sheet of the same workbooks',
    'One summed HADEG abundance per independent mouse',
    'log2(sum(2^individual-protein log2 intensity)) within each mouse',
    'Dosed minus Control',
    'Separate male and female analyses',
    'Separate 5W, 8W and 12W analyses',
    'limma::lmFit + limma::eBayes',
    'limma::duplicateCorrelation using LitterID when litters repeat',
    paste0('At least ',min_litters_per_treatment,
           ' distinct litters in each treatment group'),
    'BH adjustment across all testable HADEG features within each sex/timepoint',
    paste0('BH_FDR < ',fdr_cutoff,' and |log2FC| >= ',abs_logfc_cutoff,
           '; FDR-only significance is also reported'),
    'Run as a separate sensitivity analysis; separate FDR family',
    'Not imputed as zero; features lacking all mouse scores are not tested',
    'BH_FDR < 0.05; the optional log2FC cutoff is NOT used for plot stars/colors',
    'Litter-aware limma log2FC of per-mouse summed HADEG protein abundance'))
write_tab(parameters,'README_Analysis_Parameters')

