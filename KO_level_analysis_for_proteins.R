
needed <- c('readxl','readr','dplyr','tidyr','tibble','stringr',
            'purrr','ggplot2')
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop('Install these packages first: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({
  library(readxl); library(readr); library(dplyr); library(tidyr)
  library(tibble); library(stringr); library(purrr)
  library(ggplot2)
})

# 1. User settings ----------------------------------------------------------
root <- paste0('root/',
               'file_path')
limma_dir <- file.path(root,'Litter_analysis','Final_Litter_Blocked_Limma')
old_dir <- file.path(root,'Pathway_Enrichment_specific',
                     'high_confidence_hits','KO_summary','Pathway_intermediate_tables')
old_ko_file <- file.path(old_dir,'Protein_KO_pathway_membership.csv')
old_membership_file <- file.path(old_dir,'Expression_with_pathway_membership.csv')
old_names_file <- file.path(old_dir,'Selected_KO_name_lookup.csv')
# Additional KO roster: may contain KOs without any significant protein.
ko_reference_file <- file.path(root,'Litter_analysis','KO_map_file.xlsx')
if(!file.exists(ko_reference_file))
  stop('Missing requested KO roster workbook: ',ko_reference_file)

# Keep previous KO plots unchanged; save the new figures separately.
out_dir <- file.path(root, 'Litter_analysis',
                     'KO_Map_OkabeIto_ContrastOptimized_TIFFs')
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
audit_dir <- file.path(out_dir, 'Audit_Tables')
dir.create(audit_dir, recursive=TRUE, showWarnings=FALSE)
# Reuse the precise group definition exported by the UPDATED pathway script.
# The source filtered TSVs were inspected in that script, not here.
lookup_file <- file.path(root, 'Litter_analysis',
                         'Selected_pathway_ExactPeptides_ExactRawValues_Litter_Blocked_LIMMA',
                         'tables', 'KO_protein_group_lookup_FROM_FILTERED.csv')
if (!file.exists(lookup_file)) stop(
  'Missing exact-original-evidence lookup. Run the REVISED selected-pathway script first:\n',
  lookup_file)


# SIGNIFICANCE: inherited only for groups with agreeing existing litter-aware stats.
fdr_cutoff <- 0.05
abs_logfc_cutoff <- 1.5
# All KOs in KO_map_file.xlsx will be plotted; no significance-based KO exclusions.

# Export settings
plot_dpi <- 600
protein_max_chars <- 23
organism_max_chars <- 34


expected <- c('5W_M','8W_M','12W_M','5W_F','8W_F','12W_F')
time_levels <- c('5W','8W','12W')
sex_levels <- c('Males','Females')

# EXACT original pathway list (from previous individual-protein KO code).
selected_pathways <- tibble::tribble(
  ~analysis_theme, ~pathway_id, ~preferred_pathway_name,
  'SCFA metabolism','map00650','Butanoate metabolism',
  'SCFA metabolism','map00640','Propanoate metabolism',
  'Aromatic amino acids','map00380','Tryptophan metabolism',
  'Aromatic amino acids','map00360','Phenylalanine metabolism',
  'Aromatic amino acids','map00350','Tyrosine metabolism',
  'Aromatic amino acids','map00400','Phenylalanine, tyrosine and tryptophan biosynthesis',
  'Bile acids / host interaction','map00120','Primary bile acid biosynthesis',
  'Bile acids / host interaction','map00121','Secondary bile acid biosynthesis',
  'Bile acids / host interaction','map04976','Bile secretion',
  'Chemical stress / redox','map00480','Glutathione metabolism',
  'Chemical stress / redox','map02010','ABC transporters',
  'Chemical stress / redox','map00920','Sulfur metabolism',
  'Host-microbe interface','map00540','Lipopolysaccharide biosynthesis',
  'Host-microbe interface','map00511','Other glycan degradation',
  'Host-microbe interface','map00531','Glycosaminoglycan degradation',
  'Host-microbe interface','map00520','Amino sugar and nucleotide sugar metabolism',
  'Polyamine-related source pathways','map00330','Arginine and proline metabolism',
  'Polyamine-related source pathways','map00270','Cysteine and methionine metabolism',
  'Polyamine-related source pathways','map00220','Arginine biosynthesis',
  'Glucuronate-related metabolism','map00040','Pentose and glucuronate interconversions',
  'Glucuronate-related metabolism','map00053','Ascorbate and aldarate metabolism',
  'PAH / xenobiotic-related','map00624','Polycyclic aromatic hydrocarbon degradation',
  'PAH / xenobiotic-related','map00626','Naphthalene degradation',
  'PAH / xenobiotic-related','map00362','Benzoate degradation',
  'PAH / xenobiotic-related','map00627','Aminobenzoate degradation',
  'PAH / xenobiotic-related','map00980','Metabolism of xenobiotics by cytochrome P450',
  'PAH / xenobiotic-related','map00982','Drug metabolism - cytochrome P450',
  'PAH / xenobiotic-related','map00983','Drug metabolism - other enzymes',
  'Oxidative stress / ROS-related','map00190','Oxidative phosphorylation',
  'Oxidative stress / ROS-related','map04146','Peroxisome',
  'Lipid metabolism','map00071','Fatty acid degradation',
  'Lipid metabolism','map00564','Glycerophospholipid metabolism'
) %>% distinct(pathway_id,.keep_all=TRUE)

# 2. Helpers ----------------------------------------------------------------
clean_file <- function(s) {
  x <- stringr::str_replace_all(as.character(s),'[^A-Za-z0-9]+','_')
  x <- stringr::str_replace_all(x,'^_+|_+$','')
  substr(ifelse(is.na(x)|!nzchar(x),'Unknown',x),1,85)
}
first_valid <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x <- x[!is.na(x)&nzchar(x)&!x %in% c('NA','None','-')]
  if(length(x)) x[[1]] else NA_character_
}
most_common <- function(x) {
  x <- stringr::str_squish(as.character(x))
  x <- x[!is.na(x)&nzchar(x)&!x %in% c('NA','None','Unknown taxon','Unknown taxa','Taxon unresolved','Unknown organism','-')]
  if(!length(x)) return(NA_character_)
  tab <- sort(table(x),decreasing=TRUE)
  names(tab)[[1]]
}
pick_col <- function(nms,candidates) {
  i <- match(tolower(candidates),tolower(nms))
  i <- i[!is.na(i)]
  if(length(i)) nms[[i[[1]]]] else NA_character_
}
normalize_sex <- function(z) {
  z <- toupper(trimws(as.character(z)))
  dplyr::case_when(z %in% c('M','MALE','MALES')~'M',
                   z %in% c('F','FEMALE','FEMALES')~'F',TRUE~NA_character_)
}

# Use descriptive metadata only; never edit identifiers used for data joins.
clean_protein_label <- function(x) {
  x <- trimws(as.character(x))
  x <- sub('^>+','',x)
  x <- sub('^(sp|tr|ref|gb|emb|dbj)\\|([^|]+)\\|.*$','\\2',x,perl=TRUE)
  x <- sub('^UniRef(50|90|100)_','',x,ignore.case=TRUE)
  x <- sub('\\s+.*$','',x)
  x[is.na(x)|!nzchar(x)] <- 'Unknown protein'
  x
}

# If a text column contains Tax= / OS= / [organism], extract a taxon conservatively.
extract_taxon <- function(z) {
  vapply(as.character(z),function(s) {
    if(is.na(s)||!nzchar(trimws(s))) return(NA_character_)
    s <- stringr::str_squish(s)
    m <- stringr::str_match(s,stringr::regex('\\bTax\\s*=\\s*(.*?)(?=\\s+(?:TaxID|RepID|OX|GN|PE|SV)\\s*=|$)',ignore_case=TRUE))
    if(!is.na(m[1,2])) return(m[1,2])
    m <- stringr::str_match(s,stringr::regex('\\bOS\\s*=\\s*(.*?)(?=\\s+(?:TaxID|OX|GN|PE|SV)\\s*=|$)',ignore_case=TRUE))
    if(!is.na(m[1,2])) return(m[1,2])
    b <- stringr::str_extract_all(s,'\\[[^\\]]+\\]')[[1]]
    if(length(b)) {
      b <- gsub('^\\[|\\]$','',b)
      b <- b[!grepl('^(EC\\s*:|KO\\s*:|K[0-9]{5}$)',b,ignore.case=TRUE)]
      if(length(b)) return(tail(b,1))
    }
    NA_character_
  },character(1))
}
clean_taxon <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x,'^\\[|\\]$','')
  x <- stringr::str_replace(x,stringr::regex('^.*?\\b(?:Tax|OS)\\s*=\\s*',ignore_case=TRUE),'')
  x <- stringr::str_replace(x,stringr::regex('\\s+(?:TaxID|RepID|OX|GN|PE|SV)\\s*=.*$',ignore_case=TRUE),'')
  # Suppress strain/isolate codes for display only; retain full taxon in output.
  x <- stringr::str_replace(x,stringr::regex('\\s+(?:subsp\\.|strain|str\\.|isolate|ATCC|DSM|JCM|NCTC|NBRC)\\b.*$',ignore_case=TRUE),'')
  x <- stringr::str_replace_all(x,'[\\[\\]]','')
  x <- stringr::str_squish(x)
  x[is.na(x)|!nzchar(x)|
      grepl('^(EC\\s*:|K[0-9]{5}$|NA$)', x, ignore.case=TRUE) |
      tolower(x) %in% c('unknown taxon','unknown taxa','taxon unresolved',
                        'unknown organism','unclassified','-')] <- 'Unknown taxa'
  x
}

clean_ko_name <- function(x) {
  x <- as.character(x)
  x <- stringr::str_remove_all(x,'\\s*\\[EC[^\\]]*\\]')
  stringr::str_squish(x)
}


read_ko_reference <- function(fp) {
  tabs <- readxl::excel_sheets(fp)
  gathered <- list()
  n <- 0L
  ko_pattern <- '(?<![A-Za-z0-9])K[0-9]{5}(?![0-9])'
  
  for (sh in tabs) {
    tab <- readxl::read_excel(fp, sheet=sh, col_types='text',
                              .name_repair='unique')
    if (!nrow(tab)) next
    headers <- names(tab)
    canon <- tolower(gsub('[^A-Za-z0-9]+','',headers))
    
    # Locate the KO-ID column by both its header and actual K##### contents.
    has_ko <- vapply(tab, function(v) {
      any(stringr::str_detect(toupper(as.character(v)),ko_pattern),na.rm=TRUE)
    }, logical(1))
    id_priority <- c('ko','koid','konumber','keggid','keggko',
                     'keggortholog','keggorthologid','orthologid',
                     'ortholog','kofamid','eggnogko')
    id_named <- which(canon %in% id_priority)
    id_named <- id_named[order(match(canon[id_named],id_priority))]
    id_other <- which(grepl('ko|kegg|ortholog|kofam|eggnog',canon) &
                        !grepl('name|desc|definition|function|annotation|label',canon))
    id_candidates <- unique(c(id_named,id_other,which(has_ko)))
    id_candidates <- id_candidates[has_ko[id_candidates]]
    if (!length(id_candidates)) {
      warning('No KO IDs found in sheet: ',sh)
      next
    }
    id_col <- id_candidates[[1]]
    
    # A named KO description wins over a generic protein function/annotation.
    name_priority <- c('koname','kopreferredname','kodisplayname',
                       'kodescription','kodefinition','kofunction',
                       'orthologname','orthologdescription',
                       'keggname','keggdescription','preferredname',
                       'displayname','enzymename','name','description',
                       'definition','function','enzyme','product',
                       'genename','gene','annotation')
    exact <- which(canon %in% name_priority & seq_along(canon)!=id_col)
    exact <- exact[order(match(canon[exact],name_priority))]
    other <- which(grepl('name|description|definition|function|enzyme|product|label|gene',canon) &
                     seq_along(canon)!=id_col)
    name_candidates <- unique(c(exact,other))
    
    ko_txt <- toupper(as.character(tab[[id_col]]))
    ko_names_in_sheet <- rep(NA_character_,nrow(tab))
    for (ci in name_candidates) {
      v <- stringr::str_squish(as.character(tab[[ci]]))
      v[tolower(v) %in% c('','na','n/a','none','-','unknown')] <- NA_character_
      fill <- is.na(ko_names_in_sheet) & !is.na(v)
      ko_names_in_sheet[fill] <- v[fill]
    }
    
    # Optional fallback: an ID cell such as "K01034 | butyryl-CoA transferase".
    # A bare ID or a cell holding multiple IDs is NOT treated as a KO name.
    for (r in seq_len(nrow(tab))) {
      if (is.na(ko_txt[r])) next
      ids <- unique(unlist(stringr::str_extract_all(ko_txt[r],ko_pattern)))
      ids <- ids[!is.na(ids)]
      if (!length(ids)) next
      ko_name <- ko_names_in_sheet[r]
      if (is.na(ko_name) && length(ids)==1L) {
        raw <- stringr::str_squish(as.character(tab[[id_col]][r]))
        # Only accept a description AFTER a leading KO ID; "KO: K01034"
        # is merely an ID cell, not a name.
        suffix <- stringr::str_match(raw,
                                     stringr::regex('^K[0-9]{5}\\s*(?:[|:;,_=-]\\s*|\\s+)(.+)$',
                                                    ignore_case=TRUE))[,2]
        if (!is.na(suffix) && nzchar(stringr::str_squish(suffix)))
          ko_name <- stringr::str_squish(suffix)
      }
      n <- n+1L
      gathered[[n]] <- tibble::tibble(
        KO=ids,KO_map_name=rep(ko_name,length(ids)),
        Source_sheet=sh,Source_row=r+1L,
        KO_ID_column=headers[id_col],
        KO_name_column=if(length(name_candidates))
          paste(headers[name_candidates],collapse=' | ') else NA_character_)
    }
    message('KO map sheet "',sh,'": ID column = ',headers[id_col],
            '; preferred name column = ',
            if(length(name_candidates)) headers[name_candidates[1]] else '[none detected]')
  }
  
  if (!length(gathered)) stop('No K##### identifiers found in: ',fp)
  all_rows <- dplyr::bind_rows(gathered) %>%
    dplyr::mutate(KO_map_name=dplyr::na_if(
      stringr::str_squish(KO_map_name),''))
  
  # Flag contradictory names rather than silently hiding them.
  conflicts <- all_rows %>%
    dplyr::filter(!is.na(KO_map_name)) %>%
    dplyr::group_by(KO) %>%
    dplyr::summarise(n_distinct_names=dplyr::n_distinct(KO_map_name),
                     names_in_workbook=paste(unique(KO_map_name),collapse=' | '),
                     .groups='drop') %>%
    dplyr::filter(n_distinct_names>1L)
  
  if (nrow(conflicts)) warning(nrow(conflicts),
                               ' KO IDs have different names in the workbook; see QC_KO_map_conflicting_names.csv')
  
  # One KO = one preferred label. Prefer the first NONEMPTY workbook name.
  all_rows %>%
    dplyr::arrange(KO,dplyr::desc(!is.na(KO_map_name))) %>%
    dplyr::distinct(KO,.keep_all=TRUE)
}
ko_reference <- read_ko_reference(ko_reference_file)
message('Read ',nrow(ko_reference),' KO identifiers from KO_map_file.xlsx; ',
        sum(!is.na(ko_reference$KO_map_name)),' with workbook names')

if (all(is.na(ko_reference$KO_map_name)))
  stop('Found KO IDs in KO_map_file.xlsx but NO KO names. ',
       'Check that its name/description column was recognized. ',
       'Review the workbook headers and the KO map QC CSV.')

# 3. Read ONLY original KO and per-dataset pathway membership ---------------
for(f in c(old_ko_file,old_membership_file))
  if(!file.exists(f)) stop('Required ORIGINAL annotation file not found: ',f)

map0 <- readr::read_csv(old_ko_file,show_col_types=FALSE)
req <- c('Protein','KO','pathway_id')
if(!all(req %in% names(map0))) stop('Original KO mapping missing: ',paste(setdiff(req,names(map0)),collapse=', '))
for(col in c('pathway_name','broad_category','subcategory','Function_from_annotation'))
  if(!col %in% names(map0)) map0[[col]] <- NA_character_
map0 <- map0 %>%
  transmute(Protein=trimws(as.character(Protein)),KO=trimws(as.character(KO)),
            pathway_id=trimws(as.character(pathway_id)),
            pathway_name=as.character(pathway_name),
            broad_category=as.character(broad_category),
            subcategory=as.character(subcategory),
            Function_from_annotation=as.character(Function_from_annotation)) %>%
  filter(!is.na(Protein),nzchar(Protein),grepl('^K[0-9]{5}$',KO),
         !is.na(pathway_id),nzchar(pathway_id)) %>%
  distinct(Protein,KO,pathway_id,.keep_all=TRUE)
old_expr <- readr::read_csv(old_membership_file,show_col_types=FALSE)
req <- c('Protein','Sex','Timepoint','pathway_id')
if(!all(req %in% names(old_expr))) stop('Original expression/pathway mapping missing: ',paste(setdiff(req,names(old_expr)),collapse=', '))
old_present <- old_expr %>%
  transmute(Protein=trimws(as.character(Protein)),
            pathway_id=trimws(as.character(pathway_id)),
            SexCode=normalize_sex(Sex),
            Timepoint=toupper(trimws(as.character(Timepoint)))) %>%
  mutate(Dataset=paste(Timepoint,SexCode,sep='_')) %>%
  filter(Dataset %in% expected,!is.na(Protein),nzchar(Protein),
         !is.na(pathway_id),nzchar(pathway_id)) %>%
  distinct(Dataset,Protein,pathway_id)
if(!nrow(old_present)) stop('No original sex/timepoint-specific pathway memberships found.')


original_mapping <- map0 %>%
  left_join(selected_pathways,by='pathway_id') %>%
  filter(!is.na(analysis_theme) | KO %in% ko_reference$KO) %>%
  mutate(analysis_theme=coalesce(analysis_theme,'Additional KOs from KO_map_file'),
         pathway_name=coalesce(na_if(pathway_name,''),preferred_pathway_name,pathway_id),
         broad_category=coalesce(broad_category,'Unclassified'),
         subcategory=coalesce(subcategory,'Unclassified')) %>%
  select(-preferred_pathway_name) %>%
  distinct(Protein,KO,pathway_id,.keep_all=TRUE)
if(!nrow(original_mapping)) stop('No KO assignments matched original mapping.')


ko_reference_qc <- ko_reference %>%
  mutate(In_original_protein_KO_map=KO %in% map0$KO,
         In_original_dataset_pathway=KO %in%
           (map0 %>% semi_join(old_present,by=c('Protein','pathway_id')) %>% pull(KO)))

message('KO roster with original dataset-specific mapping: ',
        sum(ko_reference_qc$In_original_dataset_pathway),' / ',nrow(ko_reference_qc))

# 4. Load new litter-BLOCKED protein-level limma results ---------------------
read_one <- function(id) {
  fp <- file.path(limma_dir,paste0(id,'_Litter_Blocked_Protein_Results.xlsx'))
  if(!file.exists(fp)) stop('Missing NEW litter-blocked workbook: ',fp)
  tabs <- readxl::excel_sheets(fp)
  if(!'Protein_Statistics' %in% tabs) stop(id,': missing Protein_Statistics sheet')
  dat <- readxl::read_excel(fp,sheet='Protein_Statistics')
  needed_stats <- c('Protein','log2FC','p_value','BH_FDR')
  if(!all(needed_stats %in% names(dat)))
    stop(id,': Protein_Statistics missing: ',paste(setdiff(needed_stats,names(dat)),collapse=', '))
  dat$Protein <- trimws(as.character(dat$Protein))
  if(anyDuplicated(dat$Protein)) stop(id,': duplicate IDs in Protein_Statistics')
  
  # Prefer explicit original function/organism fields; never infer taxon from KO.
  ann_tab <- if('Function_Annotation' %in% tabs) 'Function_Annotation' else
    if('Protein_Annotated_Results' %in% tabs) 'Protein_Annotated_Results' else NA_character_
  if(is.na(ann_tab)) {
    warning(id,': no annotation sheet; taxonomy will be Unknown taxa')
    anno <- tibble(Protein=character(),Function=character(),Taxon_raw=character())
  } else {
    aa <- readxl::read_excel(fp,sheet=ann_tab)
    if(!'Protein' %in% names(aa)) stop(id,': annotation sheet lacks Protein')
    aa$Protein <- trimws(as.character(aa$Protein))
    if(anyDuplicated(aa$Protein)) stop(id,': duplicate IDs in annotation sheet')
    fcol <- pick_col(names(aa),c('Function','Description','Protein.Description',
                                 'Protein description','Protein_Function','Annotation'))
    tcol <- pick_col(names(aa),c('Organism','organism_name','Organism Name',
                                 'Taxon','Taxon name','Species','Taxonomy',
                                 'Organism/Taxon','Taxonomic assignment'))
    fun <- if(!is.na(fcol)) as.character(aa[[fcol]]) else rep(NA_character_,nrow(aa))
    tx <- if(!is.na(tcol)) as.character(aa[[tcol]]) else rep(NA_character_,nrow(aa))
    parsed <- extract_taxon(fun)
    tx[is.na(tx)|!nzchar(trimws(tx))] <- parsed[is.na(tx)|!nzchar(trimws(tx))]
    anno <- tibble(Protein=aa$Protein,Function=fun,Taxon_raw=tx)
  }
  parts <- strsplit(id,'_',fixed=TRUE)[[1]]
  dat %>% transmute(
    Dataset=id,Protein,
    Sex=if(parts[2]=='M') 'Males' else 'Females',
    Timepoint=parts[1],
    logFC=as.numeric(log2FC),
    P.Value=as.numeric(p_value),
    adj.P.Val=as.numeric(BH_FDR),
    limma_source_file=basename(fp)
  ) %>%
    left_join(anno,by='Protein',relationship='many-to-one') %>%
    mutate(Taxon_raw=as.character(Taxon_raw),
           Organism=clean_taxon(Taxon_raw),
           Function=as.character(Function))
}
limma_all <- purrr::map_dfr(expected,read_one)
if(!nrow(limma_all)) stop('No protein-level limma rows loaded')
if(anyDuplicated(limma_all[c('Dataset','Protein')])) stop('Duplicate dataset/protein limma statistics')


ko_protein <- limma_all %>%
  inner_join(original_mapping,by='Protein',relationship='many-to-many') %>%
  semi_join(old_present,by=c('Dataset','Protein','pathway_id')) %>%
  distinct(Dataset,Protein,KO,pathway_id,.keep_all=TRUE) %>%
  mutate(Significant=is.finite(adj.P.Val)&adj.P.Val<fdr_cutoff&
           is.finite(logFC)&abs(logFC)>=abs_logfc_cutoff,
         Direction=case_when(Significant & logFC>0 ~ 'Enriched in dosed',
                             Significant & logFC<0 ~ 'Depleted in dosed',
                             TRUE ~ 'N/S'),
         neg_log10_adj_p=ifelse(!is.na(adj.P.Val)&adj.P.Val>0,
                                -log10(adj.P.Val),NA_real_))
if(!nrow(ko_protein)) stop('No proteins matched ORIGINAL KO and dataset/pathway membership.')
if(anyDuplicated(ko_protein[c('Dataset','Protein','KO','pathway_id')]))
  stop('Duplicate dataset-protein-KO-pathway results detected')


mapping_qc <- limma_all %>%
  inner_join(original_mapping,by='Protein',relationship='many-to-many') %>%
  mutate(Present_original_dataset_pathway=
           paste(Dataset,Protein,pathway_id,sep='||') %in%
           paste(old_present$Dataset,old_present$Protein,old_present$pathway_id,sep='||')) %>%
  group_by(Dataset,pathway_id) %>%
  summarise(Global_KO_proteins=n_distinct(Protein),
            Originally_present_proteins=n_distinct(Protein[Present_original_dataset_pathway]),
            Excluded_global_only=n_distinct(Protein[!Present_original_dataset_pathway]),
            .groups='drop')

naph_qc <- tibble(Dataset=expected) %>%
  left_join(old_present %>% filter(pathway_id=='map00626') %>%
              count(Dataset,name='Original_protein_memberships'),by='Dataset') %>%
  left_join(ko_protein %>% filter(pathway_id=='map00626') %>%
              group_by(Dataset) %>%
              summarise(Matched_proteins=n_distinct(Protein),.groups='drop'),by='Dataset') %>%
  mutate(across(c(Original_protein_memberships,Matched_proteins),~replace_na(.x,0L)))

print(naph_qc,n=Inf)
if(any(naph_qc$Original_protein_memberships==0 & naph_qc$Matched_proteins>0))
  stop('Naphthalene membership mismatch: an originally absent pathway appeared.')

ko_cache <- if (file.exists(old_names_file))
  readr::read_csv(old_names_file,show_col_types=FALSE) else
    tibble(KO=character(),KO_name=character())
if (!'KO' %in% names(ko_cache))
  ko_cache <- tibble(KO=character(),KO_name=character())
if (!'KO_name' %in% names(ko_cache)) ko_cache$KO_name <- NA_character_
ko_cache <- ko_cache %>%
  transmute(KO=as.character(KO),KO_name=as.character(KO_name)) %>%
  distinct(KO,.keep_all=TRUE)

local_name <- original_mapping %>%
  group_by(KO) %>%
  summarise(Local_name=most_common(Function_from_annotation),.groups='drop')

ko_names <- ko_protein %>% distinct(KO) %>%
  left_join(ko_reference %>% select(KO,KO_map_name),by='KO') %>%
  left_join(ko_cache,by='KO') %>%
  left_join(local_name,by='KO') %>%
  mutate(
    KO_map_name=na_if(stringr::str_squish(KO_map_name),''),
    KO_name=na_if(stringr::str_squish(KO_name),''),
    Local_name=na_if(stringr::str_squish(Local_name),''),
    KO_name_source=case_when(
      !is.na(KO_map_name) ~ 'KO_map_file.xlsx',
      !is.na(KO_name) ~ 'Old cached name',
      !is.na(Local_name) ~ 'Original protein annotation',
      TRUE ~ 'KO identifier'),
    KO_display_name=coalesce(KO_map_name,KO_name,Local_name,KO),
    # Workbook labels stay intact (including any [EC ...] in the Excel name).
    KO_display_name=if_else(KO_name_source=='KO_map_file.xlsx',
                            KO_display_name,clean_ko_name(KO_display_name))
  )
ko_protein <- ko_protein %>%
  left_join(ko_names %>% select(KO,KO_display_name,KO_name_source),by='KO')

message('KO figure label sources:')
print(ko_names %>% count(KO_name_source,name='Number_of_KOs'),n=Inf)


needed_lookup <- c('Dataset','Protein','KO','ProteinGroupID','Representative',
                   'GroupSize','PeptideSet','OriginalObservedProfile','GroupOrganism')
lookup <- readr::read_csv(lookup_file,
                          col_types=readr::cols(.default=readr::col_character()),
                          show_col_types=FALSE)
missing_lookup <- setdiff(needed_lookup,names(lookup))
if (length(missing_lookup)) stop('KO grouping lookup missing columns: ',
                                 paste(missing_lookup,collapse=', '))
lookup <- lookup %>%
  mutate(across(all_of(c('Dataset','Protein','KO','ProteinGroupID')),
                ~trimws(as.character(.x)))) %>%
  filter(Dataset %in% expected, !is.na(KO),nzchar(KO)) %>%
  distinct(Dataset,KO,Protein,.keep_all=TRUE)
if(anyDuplicated(lookup[c('Dataset','KO','Protein')]))
  stop('Duplicate dataset-KO-protein keys in group lookup.')
if(anyNA(lookup$ProteinGroupID) || any(!nzchar(lookup$ProteinGroupID)))
  stop('Blank group IDs in strict lookup.')

# For every collapsed group, all three original-evidence conditions must hold.
# Missing peptide lists should NEVER be part of a multiaccession group.
lookup_qc <- lookup %>% group_by(Dataset,KO,ProteinGroupID) %>%
  summarise(N_Accession=n_distinct(Protein),
            N_Peptide_Sets=n_distinct(PeptideSet),
            N_Original_8_Mouse_Profiles=n_distinct(OriginalObservedProfile),
            Missing_Peptide=any(is.na(PeptideSet) | PeptideSet==''),
            Missing_Raw_Profile=any(is.na(OriginalObservedProfile) |
                                      OriginalObservedProfile==''),
            .groups='drop')
bad_lookup <- lookup_qc %>% filter(N_Accession>1L &
                                     (N_Peptide_Sets!=1L | N_Original_8_Mouse_Profiles!=1L |
                                        Missing_Peptide | Missing_Raw_Profile))
if(nrow(bad_lookup)) {
  readr::write_csv(bad_lookup,file.path(audit_dir,'ERROR_GroupEvidenceConflicts.csv'))
  stop('Invalid multiaccession group: peptide or original TSV profiles differ. See audit.')
}

# Accession-level stats from the ORIGINAL litter-aware run are retained.
# No protein/KO/pathway mappings are added or changed.
eligible <- ko_protein %>%
  filter(KO %in% ko_reference$KO) %>%
  distinct(Dataset,KO,Protein,.keep_all=TRUE)
missing_evidence <- eligible %>%
  anti_join(lookup,by=c('Dataset','KO','Protein'))
if(nrow(missing_evidence)) {
  readr::write_csv(missing_evidence %>%
                     select(Dataset,KO,Protein),file.path(audit_dir,'ERROR_Missing_GroupEvidence.csv'))
  stop(nrow(missing_evidence),' plotted KO/accessions lack the new strict group lookup. ',
       'Check the selected-pathway script and its original annotation sources.')
}

ko_group_members <- eligible %>%
  left_join(lookup %>%
              select(Dataset,KO,Protein,ProteinGroupID,Representative,GroupSize,
                     PeptideSet,OriginalObservedProfile,GroupOrganism),
            by=c('Dataset','KO','Protein'),relationship='many-to-one')
if(anyNA(ko_group_members$ProteinGroupID)) stop('Group join failed.')
if(anyDuplicated(ko_group_members[c('Dataset','KO','Protein')]))
  stop('Group join unexpectedly duplicated accession statistics.')

max_diff <- function(z) {
  if (!length(z) || any(!is.finite(z))) return(Inf)
  max(z)-min(z)
}
finite_median <- function(z) {
  z <- z[is.finite(z)]
  if(length(z)) median(z) else NA_real_
}


plotted_group_ids <- ko_group_members %>%
  distinct(Dataset,KO,ProteinGroupID)
missing_limma_members <- lookup %>%
  anti_join(limma_all %>% select(Dataset,Protein),by=c('Dataset','Protein')) %>%
  mutate(Affects_Plotted_Group = paste(Dataset,KO,ProteinGroupID,sep='||') %in%
           paste(plotted_group_ids$Dataset,plotted_group_ids$KO,
                 plotted_group_ids$ProteinGroupID,sep='||'))
readr::write_csv(missing_limma_members %>%
                   select(Dataset,KO,ProteinGroupID,Protein,Affects_Plotted_Group),
                 file.path(audit_dir,'QC_Lookup_Accessions_Missing_From_Limma.csv'))
if(nrow(missing_limma_members)) {
  message('Lookup accessions missing from original limma: ',
          nrow(missing_limma_members),
          '; affecting groups actually plotted: ',
          sum(missing_limma_members$Affects_Plotted_Group),
          '. See QC_Lookup_Accessions_Missing_From_Limma.csv')
}

all_group_members <- lookup %>%
  left_join(limma_all %>% select(Dataset,Protein,logFC,P.Value,adj.P.Val),
            by=c('Dataset','Protein'), relationship='many-to-one')

# A grouped p-value is displayed ONLY if ALL the original accession-level
# model results agree to numerical precision. We do NOT aggregate p-values.
ko_group_stats <- all_group_members %>%
  semi_join(plotted_group_ids,by=c('Dataset','KO','ProteinGroupID')) %>%
  group_by(Dataset,KO,ProteinGroupID) %>%
  summarise(
    N_Members_In_Group=n_distinct(Protein),
    Missing_or_Invalid_Limma_Statistics=sum(
      !is.finite(logFC) | !is.finite(P.Value) | !is.finite(adj.P.Val)),
    MemberAccessions=paste(sort(unique(Protein)),collapse='; '),
    Representative=sort(unique(Protein))[1],
    Max_logFC_difference=max_diff(logFC),
    Max_P_difference=max_diff(P.Value),
    Max_FDR_difference=max_diff(adj.P.Val),
    Group_logFC_descriptive=finite_median(logFC),
    GroupTaxon= {
      t <- unique(GroupOrganism[!is.na(GroupOrganism) & nzchar(GroupOrganism)])
      if(length(t)==1L && !t %in% c('Taxon unresolved','Unknown taxa','Unknown taxon'))
        clean_taxon(t) else 'Unknown taxa'
    },
    .groups='drop'
  ) %>%
  mutate(
    MemberStatisticsAgree=(Missing_or_Invalid_Limma_Statistics==0L &
                             Max_logFC_difference<=1e-8 &
                             Max_P_difference<=1e-8 &
                             Max_FDR_difference<=1e-8),
    NeedsReview=!MemberStatisticsAgree,
    Statistical_Status=case_when(
      Missing_or_Invalid_Limma_Statistics>0 ~
        'Group member missing or invalid limma statistics; significance unresolved',
      MemberStatisticsAgree ~ 'Existing accession-level limma results agree',
      TRUE ~ 'Member limma statistics differ; group significance unresolved')
  )


ko_plot_source <- ko_group_members %>%
  left_join(ko_group_stats %>%
              select(Dataset,KO,ProteinGroupID,N_Members_In_Group,
                     Group_logFC_descriptive,GroupTaxon,MemberStatisticsAgree,
                     NeedsReview,Statistical_Status),
            by=c('Dataset','KO','ProteinGroupID'),relationship='many-to-one') %>%
  group_by(Dataset,KO,ProteinGroupID) %>%
  arrange(Protein,.by_group=TRUE) %>%
  slice(1L) %>% ungroup() %>%
  mutate(
    logFC=if_else(MemberStatisticsAgree,logFC,Group_logFC_descriptive),
    P.Value=if_else(MemberStatisticsAgree,P.Value,NA_real_),
    adj.P.Val=if_else(MemberStatisticsAgree,adj.P.Val,NA_real_),
    Significant=is.finite(adj.P.Val) & adj.P.Val<fdr_cutoff &
      is.finite(logFC) & abs(logFC)>=abs_logfc_cutoff,
    Direction=case_when(
      NeedsReview ~ 'Group statistics need review',
      Significant & logFC>0 ~ 'Enriched in dosed',
      Significant & logFC<0 ~ 'Depleted in dosed',
      TRUE ~ 'N/S'),
    Taxon_cleaned=GroupTaxon
  )
if(anyDuplicated(ko_plot_source[c('Dataset','KO','ProteinGroupID')]))
  stop('Multiple plot rows remain per KO protein group.')
stopifnot(nrow(ko_plot_source)==nrow(ko_group_stats),
          nrow(ko_plot_source)<=nrow(eligible))

# Exports preserve every accession and every unresolved case for auditing.
readr::write_csv(lookup_qc,file.path(audit_dir,'Original_Evidence_Group_Check.csv'))
readr::write_csv(ko_group_members %>%
                   select(Dataset,KO,ProteinGroupID,Protein,Representative,GroupSize,
                          PeptideSet,OriginalObservedProfile,GroupOrganism,
                          logFC,P.Value,adj.P.Val,Significant,Sex,Timepoint),
                 file.path(audit_dir,'KO_Plot_Eligible_Accessions.csv'))
readr::write_csv(all_group_members %>%
                   semi_join(plotted_group_ids,by=c('Dataset','KO','ProteinGroupID')),
                 file.path(audit_dir,'All_Strict_Group_Members_With_Limma_Statistics.csv'))
readr::write_csv(ko_group_stats,
                 file.path(audit_dir,'All_KO_Group_Stats_Agreement.csv'))
readr::write_csv(ko_group_stats %>% filter(NeedsReview),
                 file.path(audit_dir,'REVIEW_Group_Missing_or_Differing_Limma_Stats.csv'))
readr::write_csv(ko_plot_source,
                 file.path(audit_dir,'One_Dot_Per_Protein_Group_Plot_Data.csv'))
message('Strict original-evidence grouping: ',nrow(eligible),
        ' accession/KO entries -> ',nrow(ko_plot_source),' group/KO dots; ',
        sum(ko_group_stats$NeedsReview),' groups requiring statistical review.')

ko_summary <- ko_plot_source %>%
  group_by(KO,KO_display_name,Sex,Timepoint) %>%
  summarise(n_groups=n_distinct(ProteinGroupID),
            n_significant=sum(Significant,na.rm=TRUE),
            n_enriched=sum(Significant & logFC>0,na.rm=TRUE),
            n_depleted=sum(Significant & logFC<0,na.rm=TRUE),
            n_review=sum(NeedsReview,na.rm=TRUE),
            .groups='drop')
readr::write_csv(ko_summary,file.path(audit_dir,'KO_Sex_Timepoint_Summary.csv'))



compact_plot_width <- 4.0
compact_plot_height <- 2.0
min_plot_panel_width <- 3.0
legend_title_size <- 6
legend_text_size <- 5
point_size_sig <- 1.8
point_size_ns <- 1.4
jitter_width <- 0.25   # current preferred spacing; adjust if desired
# Nonsignificant points remain light gray; significant points show taxon color.
show_timepoint_vlines <- TRUE
vline_color <- 'grey35'
vline_size <- 0.28

# Create ONE stable palette from every displayed taxon in the whole analysis.
# Taxa sharing a KO/sex plot are assigned perceptually separated colors.

# ------------------------------------------------------------
# TAXON COLORS: OKABE-ITO + DARK2 + SET1 + PAIRED
# One reproducible color per taxon, optimized for taxa that
# appear together as significant points in the same KO/sex plot.
# ------------------------------------------------------------
make_taxa_palette <- function(taxa_levels, plot_data) {
  
  taxa_levels <- sort(unique(as.character(taxa_levels)))
  taxa_levels <- taxa_levels[!is.na(taxa_levels) & nzchar(taxa_levels)]
  named_taxa <- setdiff(taxa_levels, 'Unknown taxa')
  
  # Strong, established qualitative colors; avoid the pale yellow.
  okabe_ito <- c('#E69F00', '#56B4E9', '#009E73', '#0072B2',
                 '#D55E00', '#CC79A7', '#000000')
  dark2 <- c('#1B9E77', '#D95F02', '#7570B3', '#E7298A',
             '#66A61E', '#E6AB02', '#A6761D')
  set1 <- c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3',
            '#FF7F00', '#A65628', '#F781BF')
  paired <- c('#1F78B4', '#33A02C', '#E31A1C',
              '#FF7F00', '#6A3D9A', '#B15928')
  original_palette <- unique(c(okabe_ito, dark2, set1, paired))
  
  to_lab <- function(hex) {
    grDevices::convertColor(
      t(grDevices::col2rgb(hex)) / 255,
      from = 'sRGB', to = 'Lab', scale.in = 1
    )
  }
  
  # Supply extra colors only when the original palette is not large enough.
  candidates <- original_palette
  if (length(named_taxa) > length(original_palette)) {
    extras <- unique(grDevices::hcl.colors(
      max(240L, length(named_taxa) * 5L), palette = 'Dynamic'
    ))
    extra_lab <- to_lab(extras)
    extra_chroma <- sqrt(extra_lab[, 2]^2 + extra_lab[, 3]^2)
    vivid <- extras[extra_lab[, 1] >= 35 & extra_lab[, 1] <= 75 &
                      extra_chroma >= 38]
    candidates <- unique(c(original_palette, vivid, extras))
  }
  candidates <- setdiff(candidates, c('#707070', '#BFC3C7'))
  if (length(named_taxa) > length(candidates)) {
    stop('Not enough unique colors for the number of taxa.')
  }
  
  result <- stats::setNames(character(length(named_taxa)), named_taxa)
  if (length(named_taxa)) {
    # Only significant taxon points use taxon color. N/S stays gray.
    active <- plot_data %>%
      dplyr::filter(
        is.finite(logFC),
        dplyr::coalesce(Significant, FALSE) |
          dplyr::coalesce(NeedsReview, FALSE),
        Taxon_cleaned %in% named_taxa
      ) %>%
      dplyr::distinct(KO, Sex, Taxon_cleaned)
    
    # Build a list of taxa displayed together in each KO/sex plot.
    groups <- split(active$Taxon_cleaned,
                    paste(active$KO, active$Sex, sep = '||'))
    neighbors <- stats::setNames(
      rep(list(character()), length(named_taxa)), named_taxa
    )
    for (g in groups) {
      g <- unique(g)
      for (taxon in g) {
        neighbors[[taxon]] <- union(neighbors[[taxon]],
                                    setdiff(g, taxon))
      }
    }
    
    # Prioritize taxa that share figures with the most other taxa.
    degree <- lengths(neighbors[named_taxa])
    frequency <- vapply(named_taxa, function(taxon) {
      sum(active$Taxon_cleaned == taxon)
    }, integer(1))
    assign_order <- named_taxa[order(-degree, -frequency, named_taxa)]
    
    # Calculate perceptual color distances once (CIELAB / Delta E 1976).
    lab <- to_lab(candidates)
    distances <- as.matrix(stats::dist(lab))
    assigned <- stats::setNames(rep(NA_integer_, length(named_taxa)),
                                named_taxa)
    
    for (taxon in assign_order) {
      used <- unname(assigned[!is.na(assigned)])
      available <- setdiff(seq_along(candidates), used)
      
      if (length(used)) {
        distance_all <- apply(
          distances[available, used, drop = FALSE], 1, min
        )
      } else {
        distance_all <- rep(100, length(available))
      }
      
      assigned_peers <- intersect(
        neighbors[[taxon]], names(assigned)[!is.na(assigned)]
      )
      if (length(assigned_peers)) {
        peer_idx <- unname(assigned[assigned_peers])
        distance_peers <- apply(
          distances[available, peer_idx, drop = FALSE], 1, min
        )
      } else {
        distance_peers <- distance_all
      }
      
      # Maximize separation among taxa occurring in the SAME figure;
      # also encourage separation overall and favor the original palettes.
      original_bonus <- ifelse(available <= length(original_palette), 3, 0)
      score <- distance_peers + 0.15 * distance_all + original_bonus
      assigned[[taxon]] <- available[which.max(score)]
    }
    result <- stats::setNames(candidates[unname(assigned[named_taxa])],
                              named_taxa)
  }
  
  if ('Unknown taxa' %in% taxa_levels) {
    result <- c(result, 'Unknown taxa' = '#707070')
  }
  if (anyDuplicated(unname(result))) {
    stop('Taxon colors must be globally unique.')
  }
  result
}

# Standardize taxon names, including Unknown taxa, BEFORE creating the palette.
ko_plot_source <- ko_plot_source %>%
  mutate(Taxon_cleaned = clean_taxon(Taxon_cleaned))

global_taxa <- ko_plot_source %>%
  filter(is.finite(logFC), is.finite(adj.P.Val) | NeedsReview) %>%
  distinct(Taxon_cleaned) %>%
  pull(Taxon_cleaned)
global_taxa_colors <- make_taxa_palette(global_taxa, ko_plot_source)

# Export explicit color lookup to document taxon identity across all figures.
taxon_color_manifest <- ko_plot_source %>%
  filter(Taxon_cleaned %in% names(global_taxa_colors)) %>%
  group_by(Taxon_cleaned) %>%
  summarise(N_KOs=n_distinct(KO), N_Protein_Groups=n_distinct(paste(Dataset, KO, ProteinGroupID)),
            .groups='drop') %>%
  mutate(Color_Hex=unname(global_taxa_colors[Taxon_cleaned])) %>%
  arrange(Taxon_cleaned)
readr::write_csv(taxon_color_manifest, file.path(audit_dir, 'GLOBAL_Taxon_Color_Lookup.csv'))
message('Global taxon palette: ',length(global_taxa_colors),
        ' unique colors, optimized for taxa appearing together in each KO/sex plot.')

plot_one_ko <- function(ko) {
  name_value <- ko_reference$KO_map_name[match(ko,ko_reference$KO)]
  if (!length(name_value) || is.na(name_value) || !nzchar(name_value)) {
    hit <- ko_names$KO_display_name[match(ko,ko_names$KO)]
    name_value <- if (length(hit) && !is.na(hit) && nzchar(hit)) hit else ko
  }
  
  dat <- ko_plot_source %>% filter(KO == ko)
  plot_dat <- dat %>%
    filter(is.finite(logFC), is.finite(adj.P.Val) | NeedsReview) %>%
    mutate(
      Sex=as.character(Sex),
      Timepoint=factor(as.character(Timepoint), levels=time_levels),
      Taxon_cleaned=clean_taxon(Taxon_cleaned)
    )
  
  # Significant taxa have globally consistent colors; N/S points are light gray.
  # All reviewed groups retain their taxon color and use a cross marker.
  taxa_this_ko <- sort(unique(plot_dat$Taxon_cleaned))
  taxa_this_ko <- c(setdiff(taxa_this_ko, 'Unknown taxa'),
                    intersect(taxa_this_ko, 'Unknown taxa'))
  organism_colors <- global_taxa_colors[taxa_this_ko]
  if (anyNA(organism_colors)) stop('Missing global taxon color for KO ', ko)
  
  # The two sexes share symmetric limits for the KO; safe for empty KO data.
  y_max <- if(nrow(plot_dat)) {
    max(1, ceiling(max(abs(plot_dat$logFC), na.rm=TRUE)))
  } else 1
  
  # Use all original mapped GROUPS to distinguish missing data from N/S.
  coverage <- tidyr::expand_grid(Sex=sex_levels,Timepoint=time_levels) %>%
    left_join(dat %>%
                mutate(Timepoint=as.character(Timepoint)) %>%
                count(Sex,Timepoint,name='n_matched_groups'),
              by=c('Sex','Timepoint')) %>%
    left_join(plot_dat %>%
                mutate(Timepoint=as.character(Timepoint)) %>%
                count(Sex,Timepoint,name='n_plotted_groups'),
              by=c('Sex','Timepoint')) %>%
    mutate(
      n_matched_groups=replace_na(n_matched_groups,0L),
      n_plotted_groups=replace_na(n_plotted_groups,0L),
      Status=case_when(
        n_plotted_groups>0 ~ 'Plotted',
        n_matched_groups>0 ~ 'No valid statistics',
        TRUE ~ 'Not detected'
      )
    )
  
  make_sex <- function(sex_name) {
    sex_dat <- plot_dat %>% filter(Sex == sex_name)
    missing_labels <- coverage %>%
      filter(Sex == sex_name, Status != 'Plotted') %>%
      mutate(Timepoint=factor(Timepoint, levels=time_levels),
             y=0.55*y_max, label=Status)
    sex_legend_keys <- taxa_this_ko[
      taxa_this_ko %in% unique(
        sex_dat$Taxon_cleaned[sex_dat$Significant | sex_dat$NeedsReview]
      )
    ]
    
    # Include N/S whenever nonsignificant points are plotted
    if (any(!sex_dat$Significant & !sex_dat$NeedsReview)) {
      sex_legend_keys <- c(sex_legend_keys, "N/S")
    }
    
    # Dynamic page size for long real taxa labels; NO empty legend categories.
    max_chars <- if(length(sex_legend_keys)) {
      max(nchar(sex_legend_keys, type='width'))
    } else 0
    legend_inches <- max_chars * legend_text_size * 0.57 / 72
    fig_width <- if (length(sex_legend_keys)) {
      max(compact_plot_width,min_plot_panel_width+0.38+legend_inches)
    } else compact_plot_width
    fig_height <- if (length(sex_legend_keys)) {
      max(compact_plot_height,0.86+length(sex_legend_keys)*0.17)
    } else compact_plot_height
    fig_width <- ceiling(fig_width*20)/20
    fig_height <- ceiling(fig_height*20)/20
    
    p <- ggplot(sex_dat, aes(x=Timepoint,y=logFC)) +
      {if (show_timepoint_vlines)
        geom_vline(xintercept=c(1.5,2.5),color=vline_color,
                   linewidth=vline_size)
        else NULL} +
      geom_hline(yintercept=0,color='#292929',linewidth=0.55) +
      geom_point(data=filter(sex_dat,!Significant & !NeedsReview),
                 aes(color="N/S"),
                 size=point_size_ns,shape=16,alpha=0.75,
                 position=position_jitter(width=jitter_width,height=0,seed=42),
                 show.legend=TRUE) +
      geom_point(data=filter(sex_dat,NeedsReview),
                 aes(color=Taxon_cleaned),
                 size=point_size_sig,shape=4,stroke=0.75,alpha=0.96,
                 position=position_jitter(width=jitter_width,height=0,seed=42),
                 show.legend=TRUE) +
      geom_point(data=filter(sex_dat,Significant),
                 aes(color=Taxon_cleaned),
                 size=point_size_sig,shape=16,alpha=0.98,
                 position=position_jitter(width=jitter_width,height=0,seed=42),
                 show.legend=TRUE) +
      geom_text(data=missing_labels,
                aes(x=Timepoint,y=y,label=label),
                inherit.aes=FALSE,
                color='black',
                fontface='bold',
                size=2.3) +
      scale_x_discrete(limits=time_levels,drop=FALSE,
                       expand=expansion(add=0.32)) +
      scale_y_continuous(breaks=pretty(c(-y_max,y_max),n=3),
                         expand=expansion(mult=c(0.04,0.04))) +
      coord_cartesian(ylim=c(-y_max,y_max)) +
      scale_color_manual(
        name = 'Taxa / status',
        values = c(global_taxa_colors, 'N/S' = '#BFC3C7'),
        limits = sex_legend_keys,
        breaks = sex_legend_keys,
        drop = FALSE,
        na.translate = FALSE
      ) +
      labs(title=paste0(ko,' | ',name_value,' | ',sex_name),
           x=NULL,y='log2FC (Dosed - Control)',
           caption=NULL) +
      theme_bw(base_size=9) +
      theme(
        plot.title=element_text(face='bold',size=7.2,margin=margin(b=2)),
        axis.text.x=element_text(size=8.2,color='black',face='bold'),
        axis.text.y=element_text(size=6.5,color='black'),
        axis.title.y=element_text(size=5.8,face='bold',margin=margin(r=2)),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_line(color='grey93',linewidth=0.25),
        panel.grid.minor=element_blank(),
        panel.border=element_rect(color='grey58',fill=NA,linewidth=0.35),
        legend.position=if(length(sex_legend_keys)) 'right' else 'none',
        legend.justification='center',
        legend.title=element_text(size=legend_title_size,face='bold'),
        legend.text=element_text(size=legend_text_size,face='bold'),
        plot.caption=element_text(size=5.5,color='grey35',hjust=0),
        legend.key.height=grid::unit(0.18,'cm'),
        legend.key.width=grid::unit(0.20,'cm'),
        legend.spacing.y=grid::unit(0,'cm'),
        legend.box.spacing=grid::unit(0.03,'cm'),
        legend.margin=margin(0,0,0,0),
        plot.margin=margin(4,4,3,4)
      ) +
      guides(color=guide_legend(ncol=1L,byrow=TRUE,
                                override.aes=list(shape=16,size=2.3,alpha=1,stroke=0)))
    
    stem <- paste0(clean_file(ko),'_',clean_file(name_value),'_',
                   clean_file(sex_name),'_litter_blocked.tiff')
    fp <- file.path(out_dir,stem)
    ggplot2::ggsave(fp,plot=p,width=fig_width,height=fig_height,
                    units='in',dpi=plot_dpi,device='tiff',
                    compression='lzw',bg='white',limitsize=FALSE)
    message('Saved: ',basename(fp),' | ',nrow(sex_dat),' dots; ',
            sum(missing_labels$Status=='Not detected'),' ND timepoints; ',
            sum(missing_labels$Status=='No valid statistics'),' untestable')
    invisible(fp)
  }
  invisible(lapply(sex_levels,make_sex))
}

# 9. SAVE DIRECTLY INTO ONE FOLDER: 2 TIFF FILES PER EXCEL KO ------------
requested_kos <- unique(as.character(ko_reference$KO))
stopifnot(length(requested_kos)>0L)
message('KO IDs requested in KO_map_file.xlsx: ',length(requested_kos))
message('KO IDs with matched original dataset-specific protein results: ',
        length(intersect(requested_kos,ko_plot_source$KO)))
message('Writing TIFF plots to: ',out_dir)
for (ko_now in requested_kos) {
  plot_one_ko(ko_now)
}
message('COMPLETE: ',length(requested_kos)*2L,
        ' TIFF figures (one per requested KO x sex).')
message('Audit CSVs were saved in the Audit_Tables subfolder.')
