

needed <- c('readxl', 'dplyr', 'tidyr', 'tibble', 'stringr', 'readr',
            'ggplot2', 'scales', 'writexl', 'patchwork', 'cowplot')
missing_pkgs <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_pkgs)) {
  stop('Please install missing packages: ', paste(missing_pkgs, collapse=', '))
}
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(ggplot2); library(patchwork)
})

# =========================== 1. USER SETTINGS ================================
root <- 'root_folder'
litter_root <- file.path(root, 'Litter_analysis')
workbook_dirs <- c(
  'Colon Microbiome' = file.path(litter_root, 'Final_Litter_Blocked_Limma'),
  'Small intestine Microbiome' = file.path(litter_root, 'small', 'Final_Litter_Blocked_Limma')
)

target <- 'Limosilactobacillus reuteri subsp. rodentium'
timepoint_levels <- c('5W', '8W', '12W')
group_levels <- c('Control', 'Dosed')
sex_levels <- c('Male', 'Female')
expected <- c('5W_M', '8W_M', '12W_M', '5W_F', '8W_F', '12W_F')

# Function names carried over from your previous organism-specific plotting code.
selected_protein_labels <- c(
  'LPXTG-motif cell wall anchor domain protein',
  'multi-sensor signal transduction histidine kinase'
)

# Create NEW output folder. Neither original workbooks nor older plots change.
out_dir <- file.path(
  litter_root,
  'Litter_Blocked_Organism_Plots',
  'Limosilactobacillus_reuteri_subsp_rodentium'
)
bar_dir <- file.path(out_dir, 'stacked_protein_contribution_plots')
single_dir <- file.path(out_dir, 'selected_single_protein_over_time_plots')
organism_dir <- file.path(out_dir, 'organism_mouse_total_plots')
for (p in c(out_dir, bar_dir, single_dir, organism_dir)) {
  dir.create(p, recursive=TRUE, showWarnings=FALSE)
}

# ========================== 2. HELPERS ======================================
extract_organism <- function(x) {
  # Same convention as original plotting code: [organism] or Tax=...TaxID=
  x <- as.character(x)
  bracket <- stringr::str_extract(x, '\\[[^\\]]+\\]')
  tax <- stringr::str_match(x, 'Tax=(.*?)\\s*TaxID=')[, 2]
  out <- dplyr::coalesce(bracket, tax)
  out <- stringr::str_remove_all(out, '^\\[|\\]$')
  stringr::str_squish(out)
}

clean_function <- function(x) {
  x <- as.character(x)
  x <- stringr::str_remove(x, '\\s*\\[.*?\\]\\s*$')
  x <- stringr::str_remove(x, '\\s*n=1\\s*Tax=.*$')
  x <- stringr::str_remove(x, '\\s*Tax=.*$')
  x <- stringr::str_remove(x, '^MAG:\\s*')
  stringr::str_squish(x)
}

safe_file <- function(x) {
  stringr::str_replace_all(x, '[^A-Za-z0-9]+', '_')
}

make_named_palette <- function(labels) {
  labels <- unique(as.character(labels))
  cols <- grDevices::colorRampPalette(c(
    '#4E6E7A', '#D39200', '#CC79A7', '#009E73',
    '#F0E442', '#0072B2', '#D55E00', '#999999',
    '#882255', '#44AA99', '#AA4499', '#DDCC77',
    '#117733', '#88CCEE', '#CC6677', '#332288'
  ))(length(labels))
  stats::setNames(cols, labels)
}

save_plot <- function(p, folder, stem, width=15, height=8) {
  ggplot2::ggsave(file.path(folder, paste0(stem,'.png')),
                  plot=p, width=width, height=height, dpi=300)
  ggplot2::ggsave(file.path(folder, paste0(stem,'.tiff')),
                  plot=p, width=width, height=height, dpi=300, compression='lzw')
  ggplot2::ggsave(file.path(folder, paste0(stem,'.pdf')),
                  plot=p, width=width, height=height, device='pdf')
}

# ========================= 3. PREFLIGHT =====================================
source_files <- tidyr::expand_grid(
  Tissue=names(workbook_dirs), Dataset=expected
) %>%
  mutate(
    File=file.path(unname(workbook_dirs[Tissue]),
                   paste0(Dataset, '_Litter_Blocked_Protein_Results.xlsx')),
    Exists=file.exists(File)
  )

if (any(!source_files$Exists)) {
  print(source_files %>% filter(!Exists), n=Inf, width=Inf)
  stop('One or more litter-blocked workbooks are missing. Verify colon and small intestine output folders.')
}

# ======================== 4. READ 12 WORKBOOKS ===============================
all_data <- list()
all_mouse_metadata <- list()
read_qc <- list()

for (i in seq_len(nrow(source_files))) {
  file <- source_files$File[[i]]
  ds <- source_files$Dataset[[i]]
  tissue <- source_files$Tissue[[i]]
  key <- paste(tissue, ds, sep='__')
  message('Reading: ', tissue, ' / ', ds)

  sh <- readxl::excel_sheets(file)
  required_sheets <- c('Protein_Annotated_Results', 'Individual_Mouse_Log2',
                       'Mouse_Litter_Mapping')
  if (!all(required_sheets %in% sh)) {
    stop('Required sheet(s) missing in ', file, ': ',
         paste(setdiff(required_sheets, sh), collapse=', '))
  }

  ann <- readxl::read_excel(file, sheet='Protein_Annotated_Results')
  # Restore the litter-blocked log2FC column name
  if ("log2FC.x" %in% names(ann)) {
    ann <- ann %>%
      rename(log2FC = log2FC.x)
  }
  expr <- readxl::read_excel(file, sheet='Individual_Mouse_Log2')
  meta <- readxl::read_excel(file, sheet='Mouse_Litter_Mapping')

  required_ann <- c('Protein','Function','log2FC','BH_FDR')
  required_meta <- c('MouseID','OriginalColumn','Treatment','LitterID')
  if (!all(required_ann %in% names(ann))) {
    stop(ds, ' / ', tissue, ': missing annotated results columns: ',
         paste(setdiff(required_ann, names(ann)), collapse=', '),
         '\nAvailable columns: ', paste(names(ann),collapse=', '))
  }
  if (!'Protein' %in% names(expr)) stop(ds, ': no Protein column in Individual_Mouse_Log2')
  if (!all(required_meta %in% names(meta))) {
    stop(ds, ': missing mouse mapping columns: ',
         paste(setdiff(required_meta,names(meta)),collapse=', '))
  }

  meta <- meta %>%
    transmute(
      Dataset=ds, Tissue=tissue,
      MouseID=as.character(MouseID),
      OriginalColumn=as.character(OriginalColumn),
      Treatment=tolower(trimws(as.character(Treatment))),
      LitterID=as.character(LitterID),
      Group=case_when(Treatment=='control'~'Control',
                      Treatment=='dosed'~'Dosed', TRUE~NA_character_),
      Timepoint=strsplit(ds, '_', fixed=TRUE)[[1]][1],
      Sex=ifelse(endsWith(ds,'_M'), 'Male','Female')
    )

  expect_cols <- c(paste0('Control_',1:4),paste0('Dosed_',1:4))
  if (nrow(meta)!=8L || anyDuplicated(meta$MouseID) ||
      anyDuplicated(meta$OriginalColumn) ||
      !setequal(meta$OriginalColumn, expect_cols) ||
      anyNA(meta$LitterID) || any(meta$LitterID=='') ||
      anyNA(meta$Group) || !all(meta$MouseID %in% names(expr)) ||
      any(meta$Group != ifelse(startsWith(meta$OriginalColumn,'Control_'), 'Control','Dosed'))) {
    print(meta)
    stop(ds, ' / ', tissue, ': mismatched individual mouse columns or litter metadata.')
  }
  if ('MatchStatus' %in% names(readxl::read_excel(file, sheet='Mouse_Litter_Mapping'))) {
    original_map <- readxl::read_excel(file, sheet='Mouse_Litter_Mapping')
    if (any(as.character(original_map$MatchStatus) != 'MATCHED')) {
      stop(ds, ' / ', tissue, ': mouse-to-litter MatchStatus contains unmatched records.')
    }
  }

  ann <- ann %>%
    mutate(Protein=trimws(as.character(Protein)),
           Function=as.character(Function),
           log2FC=as.numeric(log2FC),
           BH_FDR=as.numeric(BH_FDR),
           Organism=extract_organism(Function),
           Protein_label=clean_function(Function)) %>%
    select(Protein,Function,Organism,Protein_label,log2FC,BH_FDR)
  expr <- expr %>% mutate(Protein=trimws(as.character(Protein)))
  if (anyNA(ann$Protein) || anyDuplicated(ann$Protein) ||
      anyNA(expr$Protein) || anyDuplicated(expr$Protein)) {
    stop(ds, ' / ', tissue, ': missing/duplicated protein accessions.')
  }
  if (!all(vapply(expr[meta$MouseID],is.numeric,logical(1)))) {
    stop(ds, ' / ', tissue, ': non-numeric mouse intensities.')
  }
  if (any(!is.finite(as.matrix(expr[meta$MouseID])))) {
    stop(ds, ' / ', tissue, ': missing or nonfinite mouse intensities.')
  }

  # Match organism names exactly as in your previous plotting workflow.
  ann_target <- ann %>%
    filter(!is.na(Organism),
           tolower(str_squish(Organism))==tolower(str_squish(target)))
  read_qc[[key]] <- tibble(
    Tissue=tissue,Dataset=ds,
    N_annotated_proteins=nrow(ann),
    N_target_organism_proteins=nrow(ann_target),
    N_litters_control=n_distinct(meta$LitterID[meta$Group=='Control']),
    N_litters_dosed=n_distinct(meta$LitterID[meta$Group=='Dosed'])
  )
  all_mouse_metadata[[key]] <- meta
  if (!nrow(ann_target)) {
    message('  Target organism not detected in ', ds, ' / ', tissue)
    next
  }

  matching_expr <- expr %>%
    filter(Protein %in% ann_target$Protein) %>%
    select(Protein, all_of(meta$MouseID))
  if (nrow(matching_expr)!=nrow(ann_target)) {
    stop(ds,' / ',tissue, ': target annotations and expression protein IDs do not match.')
  }

  all_data[[key]] <- matching_expr %>%
    pivot_longer(cols=all_of(meta$MouseID),
                 names_to='MouseID',values_to='log2_intensity') %>%
    left_join(ann_target, by='Protein', relationship='many-to-one') %>%
    left_join(meta %>% select(MouseID,OriginalColumn,LitterID,Group,
                              Dataset,Tissue,Sex,Timepoint),
              by='MouseID', relationship='many-to-one')
}

source_qc <- bind_rows(read_qc)
readr::write_csv(source_qc, file.path(out_dir,'Source_coverage_and_litter_QC.csv'))
all_meta <- bind_rows(all_mouse_metadata)
readr::write_csv(all_meta, file.path(out_dir,'Mouse_litter_mapping_both_tissues.csv'))

# The same mouse ID across the two tissues must have the same assignment.
assignment_qc <- all_meta %>%
  group_by(Dataset,MouseID) %>%
  summarise(N_litters=n_distinct(LitterID), N_groups=n_distinct(Group),
            N_tissues=n_distinct(Tissue),.groups='drop')
if (any(assignment_qc$N_litters!=1L | assignment_qc$N_groups!=1L)) {
  print(assignment_qc %>% filter(N_litters!=1L | N_groups!=1L), n=Inf)
  stop('The same mouse has conflicting tissue-specific litter/group metadata.')
}

if (!length(all_data)) {
  stop('Target organism was not found in any workbook. Check the Function organism labels and target spelling.')
}
df_org <- bind_rows(all_data) %>%
  mutate(Group=factor(Group,levels=group_levels),
         Sex=factor(Sex,levels=sex_levels),
         Timepoint=factor(Timepoint,levels=timepoint_levels),
         linear_intensity=2^log2_intensity)
if (any(!is.finite(df_org$linear_intensity))) {
  stop('Nonfinite linear intensities after 2^log2 conversion.')
}
readr::write_csv(df_org, file.path(out_dir,'Target_organism_mouse_protein_long.csv'))

# ======================= 5. SUPPLEMENTARY TABLES =============================
# Mirrors the original protein-label summary; n_mice is counted separately.
protein_summary <- df_org %>%
  group_by(Tissue,Sex,Timepoint,Group,Protein_label) %>%
  summarise(
    n_protein_accessions=n_distinct(Protein),
    n_mice=n_distinct(MouseID),
    n_values=n(),
    mean_log2_intensity=mean(log2_intensity),
    sd_log2_intensity=sd(log2_intensity),
    median_log2_intensity=median(log2_intensity),
    mean_linearized_intensity=mean(linear_intensity),
    sd_linearized_intensity=sd(linear_intensity),
    median_linearized_intensity=median(linear_intensity),
    .groups='drop'
  ) %>%
  arrange(Tissue,Sex,Timepoint,Protein_label,Group)

# Each stacked bar element is mean(2^log2 intensity) for that protein label
# across protein accessions and individual mice. Same as the old plotting code.
bar_data <- protein_summary %>%
  transmute(Tissue,Sex,Timepoint,Group,Protein_label,
            abundance_linear=mean_linearized_intensity,
            n_protein_accessions,n_mice)

# ONE organism-total abundance value per individual mouse (descriptive).
organism_mouse <- df_org %>%
  group_by(Tissue,Sex,Timepoint,Dataset,Group,MouseID,LitterID) %>%
  summarise(
    n_proteins=n_distinct(Protein),
    total_linear_abundance=sum(linear_intensity),
    log2_total_abundance=log2(total_linear_abundance),
    .groups='drop'
  )

# ONE mouse measurement per functional protein label, even if multiple
# proteins carry the same function text. This avoids protein pseudoreplication.
label_mouse <- df_org %>%
  group_by(Tissue,Sex,Timepoint,Dataset,Group,MouseID,LitterID,Protein_label) %>%
  summarise(
    log2_intensity=median(log2_intensity),
    n_protein_accessions=n_distinct(Protein),
    .groups='drop'
  )


label_stats <- df_org %>%
  distinct(Tissue,Sex,Timepoint,Dataset,Protein_label,Protein,log2FC,BH_FDR) %>%
  group_by(Tissue,Sex,Timepoint,Dataset,Protein_label) %>%
  summarise(
    n_accessions=n_distinct(Protein),
    accessions=paste(sort(unique(Protein)),collapse='; '),
    blocked_protein_log2FC=if (n_distinct(Protein)==1L) first(log2FC) else NA_real_,
    blocked_protein_BH_FDR=if (n_distinct(Protein)==1L) first(BH_FDR) else NA_real_,
    statistic_applicability=ifelse(n_distinct(Protein)==1L,
      'Single accession: protein-level FDR may be displayed (exploratory)',
      'Multiple accessions: no protein-level FDR for aggregated label'),
    .groups='drop'
  )

writexl::write_xlsx(protein_summary,
  file.path(out_dir,'Supplementary_Table_L_rodentium_protein_summary_LITTER.xlsx'))
writexl::write_xlsx(bar_data,
  file.path(out_dir,'protein_contribution_summary_mean_linearized_abundance_LITTER.xlsx'))
writexl::write_xlsx(list(
  Mouse_Litter_Mapping=all_meta,
  Source_Coverage_QC=source_qc,
  Organism_Mouse_Abundance=organism_mouse,
  Label_Mouse_Median=label_mouse,
  Label_Statistics=label_stats
),file.path(out_dir,'L_rodentium_mouse_level_and_statistical_reference.xlsx'))

# ====================== 6. STACKED BAR GRAPHS ================================
for (tissue in names(workbook_dirs)) {
  for (sx in sex_levels) {
    sub <- bar_data %>% filter(Tissue==tissue,as.character(Sex)==sx)
    if (!nrow(sub)) next
    levels_protein <- sub %>%
      group_by(Protein_label) %>%
      summarise(total=sum(abundance_linear),.groups='drop') %>%
      arrange(desc(total)) %>% pull(Protein_label)
    palette <- make_named_palette(levels_protein)
    sub <- sub %>%
      mutate(Protein_label=factor(Protein_label,levels=levels_protein),
             Timepoint=factor(Timepoint,levels=timepoint_levels))

    empty_timepoints <- setdiff(timepoint_levels, as.character(unique(sub$Timepoint)))
    nd <- tibble(Timepoint=factor(empty_timepoints,levels=timepoint_levels),
                 x=rep(1.5,length(empty_timepoints)),
                 y=rep(max(sub$abundance_linear)*0.5,length(empty_timepoints)),
                 label=rep('Not detected',length(empty_timepoints)))

    p <- ggplot(sub,aes(x=Group,y=abundance_linear,fill=Protein_label)) +
      geom_col(width=0.72,color='black',linewidth=0.15) +
      facet_wrap(~Timepoint,nrow=1,drop=FALSE) +
      geom_text(data=nd,aes(x=x,y=y,label=label),inherit.aes=FALSE,
                fontface='bold.italic',size=5) +
      scale_fill_manual(values=palette,limits=levels_protein,
                        breaks=levels_protein,drop=TRUE) +
      scale_x_discrete(drop=FALSE) +
      scale_y_continuous(labels=scales::scientific_format(digits=2),
                         expand=expansion(mult=c(0,0.08))) +
      labs(x=NULL,y='Mean linearized protein abundance',fill='Protein') +
      theme_bw(base_size=12) +
      theme(
        strip.text=element_text(size=22,face='bold'),
        axis.text.x=element_text(size=16,face='bold',color='black'),
        axis.title.y=element_text(size=15,face='bold'),
        axis.text.y=element_text(size=15,face='bold',color='black'),
        panel.grid.major.x=element_blank(),
        legend.position='right',
        legend.title=element_text(size=15,face='bold'),
        legend.text=element_text(size=11,face='bold'),
        legend.key.height=grid::unit(0.45,'cm'),
        legend.key.width=grid::unit(0.45,'cm')
      ) +
      guides(fill=guide_legend(ncol=1,byrow=TRUE))

    legend <- cowplot::get_legend(p)
    if (!is.null(legend)) {
      final_plot <- patchwork::wrap_plots(
        p+theme(legend.position='none'),
        patchwork::wrap_elements(legend),
        ncol=2,widths=c(3,1.2))
    } else {
      final_plot <- p
    }
    stem <- paste0(safe_file(target),'_',safe_file(sx),'_',safe_file(tissue),
                   '_5W_8W_12W_mean_linearized_abundance_LITTER')
    save_plot(final_plot,bar_dir,stem,width=15,height=8)
    message('Saved stacked bars: ',sx,' / ',tissue)
  }
}

# ======================= 7. SINGLE FUNCTION PLOTS ============================
# Boxplots and points are ONE ROW PER MOUSE. Stars use new blocked protein BH
# FDR only when the exact functional label represents one protein accession.
plot_single_protein_across_time <- function(protein_name_clean,
                                            save_plots=TRUE) {
  if (!protein_name_clean %in% label_mouse$Protein_label) {
    warning('Protein label not detected: ',protein_name_clean)
    return(invisible(NULL))
  }
  for (tissue in names(workbook_dirs)) {
    for (sx in sex_levels) {
      one <- label_mouse %>%
        filter(Tissue==tissue,as.character(Sex)==sx,
               Protein_label==protein_name_clean) %>%
        mutate(Timepoint=factor(Timepoint,levels=timepoint_levels))
      if (!nrow(one)) next
      one_stats <- label_stats %>%
        filter(Tissue==tissue,as.character(Sex)==sx,
               Protein_label==protein_name_clean) %>%
        mutate(Timepoint=factor(Timepoint,levels=timepoint_levels),
               signif_label=case_when(
                 n_accessions!=1L | is.na(blocked_protein_BH_FDR) ~ '',
                 blocked_protein_BH_FDR<=0.0001 ~ '****',
                 blocked_protein_BH_FDR<=0.001 ~ '***',
                 blocked_protein_BH_FDR<=0.01 ~ '**',
                 blocked_protein_BH_FDR<=0.05 ~ '*',
                 TRUE ~ 'ns'
               ))

      ymin <- min(one$log2_intensity)
      ymax <- max(one$log2_intensity)
      yrange <- ymax-ymin
      if (!is.finite(yrange) || yrange==0) yrange <- 1
      y_star <- ymax+yrange*0.12
      y_top <- ymax+yrange*0.28
      one_stats <- one_stats %>% mutate(y=y_star, x=1.5)
      empty_tp <- setdiff(timepoint_levels,as.character(unique(one$Timepoint)))
      missing_df <- tibble(Timepoint=factor(empty_tp,levels=timepoint_levels),
                           x=rep(1.5,length(empty_tp)),
                           y=rep(ymin+(y_top-ymin)*0.5,length(empty_tp)),
                           label=rep('Not detected',length(empty_tp)))

      p <- ggplot(one,aes(x=Group,y=log2_intensity,fill=Group)) +
        geom_boxplot(alpha=0.75,width=0.65,outlier.shape=NA,na.rm=TRUE) +
        geom_point(position=position_jitter(width=0.12,height=0,seed=123),
                   size=2.3,alpha=0.9,na.rm=TRUE) +
        geom_text(data=one_stats %>% filter(signif_label!=''),
                  aes(x=x,y=y,label=signif_label),inherit.aes=FALSE,
                  size=6,fontface='bold') +
        geom_text(data=missing_df,aes(x=x,y=y,label=label),
                  inherit.aes=FALSE,size=5,fontface='bold.italic') +
        scale_fill_manual(values=c(Control='steelblue',Dosed='firebrick')) +
        scale_x_discrete(drop=FALSE) +
        facet_wrap(~Timepoint,nrow=1,drop=FALSE) +
        labs(x=NULL,y='Mouse-level median log2 protein intensity',
          caption='Stars: litter-blocked protein BH FDR only for single-accession labels; exploratory (two litters per treatment).') +
        coord_cartesian(ylim=c(ymin,y_top)) +
        theme_bw(base_size=16) +
        theme(
          strip.background=element_rect(fill='grey92',color='black',linewidth=0.8),
          strip.text=element_text(size=20,face='bold',color='black'),
          axis.text.x=element_text(size=17,face='bold',color='black'),
          axis.title.y=element_text(size=17,face='bold',margin=margin(r=14)),
          axis.text.y=element_text(size=16,face='bold',color='black'),
          axis.line=element_line(color='black',linewidth=0.7),
          panel.border=element_rect(color='black',fill=NA,linewidth=0.8),
          panel.grid.minor=element_blank(),
          panel.spacing=grid::unit(1,'lines'),legend.position='none',
          plot.caption=element_text(size=10,hjust=0)
        )

      if (save_plots) {
        stem <- paste0(safe_file(target),'_',safe_file(protein_name_clean),'_',
                       safe_file(sx),'_',safe_file(tissue),'_single_protein_LITTER')
        save_plot(p,single_dir,stem,width=14,height=8)
      }
      print(p)
    }
  }
  invisible(TRUE)
}

for (protein_label in selected_protein_labels) {
  plot_single_protein_across_time(protein_label)
}


for (tissue in names(workbook_dirs)) {
  for (sx in sex_levels) {
    one <- organism_mouse %>%
      filter(Tissue==tissue,as.character(Sex)==sx) %>%
      mutate(Timepoint=factor(Timepoint,levels=timepoint_levels))
    if (!nrow(one)) next
    p <- ggplot(one,aes(Group,log2_total_abundance,fill=Group)) +
      geom_boxplot(alpha=0.75,width=0.65,outlier.shape=NA) +
      geom_point(position=position_jitter(width=0.12,height=0,seed=321),size=2.5) +
      facet_wrap(~Timepoint,nrow=1,drop=FALSE) +
      scale_fill_manual(values=c(Control='steelblue',Dosed='firebrick')) +
      scale_x_discrete(drop=FALSE) +
      labs(x=NULL,y='Log2 total organism protein abundance per mouse',
           caption='Descriptive: sum of linearized intensities across unique target-organism proteins per mouse; no aggregate-level statistical test.') +
      theme_bw(base_size=16) +
      theme(strip.text=element_text(size=18,face='bold'),
            axis.text.x=element_text(size=16,face='bold'),
            legend.position='none',
            plot.caption=element_text(size=10,hjust=0))
    stem <- paste0(safe_file(target),'_',safe_file(sx),'_',safe_file(tissue),
                   '_organism_total_by_mouse_LITTER')
    save_plot(p,organism_dir,stem,width=14,height=8)
  }
}


