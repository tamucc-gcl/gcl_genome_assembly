# Main report presentation. All inputs are explicit current-run artifacts.
suppressPackageStartupMessages({library(tidyverse);library(jsonlite);library(patchwork)})
argv <- commandArgs(trailingOnly=TRUE)
if(length(argv)!=1) stop('Supply report_spec.json')
spec <- fromJSON(argv[[1]],simplifyVector=FALSE)
out <- 'reports/main'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
chars <- function(x) unlist(x,use.names=FALSE)
optional_tsv <- function(path) {
    if(is.null(path) || !file.exists(path) || grepl('^NO_',basename(path))) return(tibble())
    read_delim(path,delim=if(grepl('\\.csv$',path)) ',' else '\t',
               col_types=cols(.default=col_character()),show_col_types=FALSE)
}
cell <- function(x) {
    x <- ifelse(is.na(x),'Not measured',as.character(x))
    str_replace_all(str_replace_all(x,'[\r\n]+',' '),fixed('|'),'\\|')
}
md_table <- function(x) {
    if(nrow(x)==0) return(character())
    c(paste0('| ',paste(cell(names(x)),collapse=' | '),' |'),
      paste0('| ',paste(rep('---',ncol(x)),collapse=' | '),' |'),
      apply(x,1,function(r) paste0('| ',paste(cell(r),collapse=' | '),' |')))
}
link <- function(label,path) paste0('[',label,'](',str_replace_all(path,' ','%20'),')')
report_label <- function(paths) {
    labels <- c(clipping_report.md='Graph sequence accounting',source_report.md='Assembly representation and review flags',
        sequence_graph_audit.md='Source-coordinate validation',comparison_report.md='Haplotype similarity and comparisons',
        regional_report.md='Chromosome-level sharing',regional_figures.md='Chromosome sharing figures',
        sharing_report.md='Core and accessory sequence sharing',sharing_figures.md='Sharing and growth figures',
        variant_report.md='Variant export audit',pangenome_input_report.md='Graph inputs and eligibility',
        pangenome_build_report.md='Graph construction',pangenome_report.md='Graph summary')
    result <- unname(labels[basename(paths)])
    ifelse(is.na(result),basename(paths),result)
}
artifact_rows <- if(length(spec$artifacts)) bind_rows(lapply(spec$artifacts,as_tibble)) else tibble(type=character(),id=character(),path=character(),note=character())
if(!'note' %in% names(artifact_rows)) artifact_rows$note <- ''
artifacts <- artifact_rows |> distinct(type,id,path,.keep_all=TRUE) |> arrange(type,id,path)
write_tsv(artifacts |> select(type,id,path,note),file.path(out,'artifact_registry.tsv'))
elig <- fromJSON(spec$eligibility,simplifyVector=FALSE)
assemblies <- map_dfr(elig$assemblies,function(x) tibble(id=x$id,sample=x$sample,taxid=as.character(x$taxid),
    species=if(is.null(x$species)) '' else x$species,
    haplotype=x$haplotype,completion=x$completion,graph_member=isTRUE(x$graph_member),fasta=x$fasta)) |>
    filter(completion=='finalized') |> arrange(taxid,sample,haplotype)
if(anyDuplicated(assemblies$id)) stop('Duplicate finalized identity in report')
species_label <- function(taxid) {
    labels <- unique(assemblies$species[assemblies$taxid==taxid & nzchar(assemblies$species)])
    paste0(if(length(labels)>0) paste0(paste(labels,collapse=' / '),' · ') else '', 'taxid ',taxid)
}
names_data <- map_dfr(chars(spec$name_maps),function(path) optional_tsv(path) |>
    mutate(assembly=sub('\\.name_map\\.tsv$','',basename(path)),length=as.numeric(length)))
if(anyDuplicated(names_data[c('assembly','new_name')])) stop('Duplicate final scaffold identity')
if(!all(names_data$assembly %in% assemblies$id) || any(!is.finite(names_data$length) | names_data$length<=0)) stop('Invalid final name-map inventory')
write_tsv(names_data,file.path(out,'scaffold_inventory.tsv'))
assembly_lengths <- names_data |> group_by(assembly) |>
    summarise(assembly_bp=sum(length),chromosome_assigned_bp=sum(length[str_detect(new_name,'^chr')]),.groups='drop') |>
    left_join(assemblies |> select(assembly=id,sample,haplotype),by='assembly') |>
    mutate(chromosome_assigned_percent=100*chromosome_assigned_bp/assembly_bp)
write_tsv(assembly_lengths,file.path(out,'assembly_lengths.tsv'))

md <- c('# Genome assembly results','',paste0(nrow(assemblies),' finalized assemblies from ',n_distinct(assemblies$sample),' individuals.'),'')
attention <- tibble(Area=character(),Item=character(),Issue=character(),Details=character())
sample_capabilities <- if(length(spec$run$samples)) bind_rows(lapply(spec$run$samples,as_tibble)) else tibble(sample=character(),hifi=logical(),hic=logical(),shortread=logical())
qc <- optional_tsv(spec$qc)
if(nrow(qc)>0) {
    hap_columns <- grep('^hap[0-9]+$',names(qc),value=TRUE)
    roles <- assemblies |> transmute(sample_id=sample,
        qc_haplotype=if_else(haplotype=='primary','hap1',haplotype),haplotype)
    kmer_title <- 'K-mer completeness (%)'
    q <- qc |> filter(stage=='final') |>
        pivot_longer(all_of(hap_columns),names_to='qc_haplotype',values_to='value') |>
        mutate(value=as.numeric(value)) |>
        inner_join(roles,by=c('sample_id','qc_haplotype'))
    if(anyDuplicated(q[c('sample_id','haplotype','analysis','metric')])) stop('Duplicate final QC measurement')
    totals <- q |> filter(analysis=='busco',metric=='total_busco') |> select(sample_id,haplotype,total=value)
    q <- q |> filter((analysis=='busco' & metric=='complete') |
        (analysis=='merqury' & metric %in% c('qv','kmer_completeness')) |
        (analysis=='quast' & metric=='N50')) |>
        left_join(totals,by=c('sample_id','haplotype')) |>
        mutate(measurement=case_when(metric=='complete' & !is.na(total) & total>0 ~ 100*value/total,
                                    metric=='complete' ~ NA_real_,metric=='N50' ~ value/1e6,TRUE ~ value),
               title=case_when(metric=='complete' ~ 'BUSCO complete (%)',metric=='qv' ~ 'QV',
                              metric=='kmer_completeness' ~ kmer_title,TRUE ~ 'Scaffold N50 (Mb)'),
               display=if_else(is.na(measurement),'Not measured',sprintf('%.1f',measurement)))
    write_tsv(q,file.path(out,'quality_values.tsv'))
    grid <- crossing(assemblies |> select(sample_id=sample,haplotype),
                     title=c('BUSCO complete (%)','QV',kmer_title,'Scaffold N50 (Mb)')) |>
        left_join(q |> select(sample_id,haplotype,title,display),by=c('sample_id','haplotype','title')) |>
        mutate(display=replace_na(display,'Not measured')) |> arrange(sample_id,haplotype) |>
        group_by(sample_id,title) |> summarise(value=paste(display,collapse=' / '),.groups='drop') |>
        pivot_wider(names_from=title,values_from=value) |>
        select(Individual=sample_id,all_of(c('BUSCO complete (%)','QV',kmer_title,'Scaffold N50 (Mb)')))
    length_display <- assembly_lengths |> arrange(sample,haplotype) |> group_by(sample) |>
        summarise(`Assembly length (Mb)`=paste(sprintf('%.1f',assembly_bp/1e6),collapse=' / '),
                  `Chromosome-assigned (%)`=paste(sprintf('%.1f',chromosome_assigned_percent),collapse=' / '),.groups='drop')
    grid <- grid |> left_join(length_display,by=c('Individual'='sample'))
    combined <- if('both' %in% names(qc)) qc |> filter(stage=='final',analysis=='merqury',metric=='kmer_completeness') |>
        transmute(sample_id,combined=as.numeric(both)) else tibble(sample_id=character(),combined=numeric())
    if(anyDuplicated(combined$sample_id)) stop('Duplicate combined k-mer completeness measurement')
    combined <- assemblies |> group_by(sample) |> summarise(copies=n(),primary=all(haplotype=='primary'),.groups='drop') |>
        left_join(combined,by=c('sample'='sample_id'))
    write_tsv(combined,file.path(out,'combined_kmer_values.tsv'))
    grid <- grid |> left_join(combined,by=c('Individual'='sample')) |>
        mutate(!!kmer_title:=if_else(copies>1,
            paste0(if_else(is.na(combined),'Not measured',sprintf('%.1f',combined)),' (',.data[[kmer_title]],')'),.data[[kmer_title]])) |>
        select(-copies,-primary,-combined)
    md <- c(md,'## Assembly quality','',
        'Paired values are **hap1 / hap2** (additional copies follow in haplotype order). K-mer completeness shows the **combined diploid value (hap1 / hap2)**; for more than two copies it shows their combined value, and a primary/collapsed assembly has one value. Combined completeness is measured directly, not calculated by adding haplotype percentages. Chromosome-assigned percentage is the fraction of total scaffold length with a chromosome name in the final harmonized name map, including composite names; it is not the fraction painted by alignment.','',md_table(grid),'',
        link('All QC measurements','qc/assembly/assembly_qc_metrics.csv'))
    interactive <- artifacts |> filter(type=='interactive_qc')
    if(nrow(interactive)>0) md <- c(md,paste0(' · ',link('Interactive QC rankings',interactive$path[[1]])))
    md <- c(md,'')
}

read_values <- tibble()
if(nrow(qc)>0) {
    read_values <- qc |> filter(stage=='final',
        (analysis=='mapped_reads' & metric=='coverage') |
        (analysis=='hic_contact' & metric %in% c('retention_pct','trans_to_cis'))) |>
        pivot_longer(all_of(hap_columns),names_to='qc_haplotype',values_to='value') |>
        inner_join(roles,by=c('sample_id','qc_haplotype')) |> mutate(value=as.numeric(value))
    if(anyDuplicated(read_values[c('sample_id','haplotype','analysis','metric')])) stop('Duplicate final read-support measurement')
    write_tsv(read_values,file.path(out,'final_read_support.tsv'))
    if(any(is.finite(read_values$value))) {
        metrics <- read_values |> filter(is.finite(value)) |> distinct(metric) |> pull(metric)
        labels <- c(coverage='Mean mapped-read depth (×)',retention_pct='Hi-C retention (%)',trans_to_cis='Hi-C trans:cis')
        read_display <- crossing(assemblies |> select(sample_id=sample,haplotype),metric=metrics) |>
            left_join(read_values |> select(sample_id,haplotype,metric,value),by=c('sample_id','haplotype','metric')) |>
            mutate(display=if_else(is.na(value),'Not measured',if_else(metric=='trans_to_cis',sprintf('%.3f',value),sprintf('%.1f',value)))) |>
            arrange(sample_id,haplotype) |> group_by(sample_id,metric) |>
            summarise(value=paste(display,collapse=' / '),.groups='drop') |>
            mutate(metric=unname(labels[metric])) |> pivot_wider(names_from=metric,values_from=value) |> rename(Individual=sample_id)
        if(nrow(sample_capabilities)>0) {
            read_display <- read_display |> left_join(sample_capabilities |> transmute(Individual=sample,
                `QC reads`=case_when(hifi ~ 'HiFi',shortread ~ 'Shotgun',TRUE ~ 'Not specified')),by='Individual') |> relocate(`QC reads`,.after=Individual)
        }
        if(any(artifacts$type=='contact_map')) read_display <- read_display |> mutate(`Contact maps`=map_chr(Individual,function(s) {
            ids <- assemblies$id[assemblies$sample==s]
            maps <- artifacts |> filter(type=='contact_map',id %in% ids) |> mutate(resolution=as.numeric(str_match(path,'_([0-9]+)bp_')[,2])) |> arrange(id,desc(resolution),path) |> group_by(id) |> slice_head(n=1) |> ungroup()
            if(nrow(maps)==0) return('Not available')
            paste(link(assemblies$haplotype[match(maps$id,assemblies$id)],maps$path),collapse=' · ')
        }))
        md <- c(md,'## Final read support','',md_table(read_display),'',
            'Measurements map QC reads to the finalized assemblies; paired values are hap1 / hap2. Mean depth is length-weighted samtools coverage using HiFi or shotgun reads as indicated, and is not a unique-molecule count. Hi-C retention is the reported retained-pair percentage of mapped pairs; trans:cis is the retained trans-pair count divided by cis pairs. Cis/trans use scaffold identity, so fragmentation itself changes the ratio. Neither depth nor a genome-wide contact ratio establishes support for an individual join. Not measured means unavailable, rather than zero.','',
            link('Source measurements','reports/main/final_read_support.tsv'),'')
    }
}

# Consecutive available measurements, with explicit stages and units, not causal claims.
stage_changes_present <- nrow(qc)>0 && n_distinct(qc$stage)>1
if(stage_changes_present) {
    stage_order <- c('ctg.base','ctg.org','ctg.purged','ctg.polished','ctg.cor','ctg.deco',
        'scaf.base','scaf.cor','scaf2','chim.brk','gap_fill','teloclip','final')
    measurements <- qc |> pivot_longer(all_of(hap_columns),names_to='qc_haplotype',values_to='raw') |>
        inner_join(roles,by=c('sample_id','qc_haplotype')) |> mutate(raw=as.numeric(raw)) |>
        group_by(sample_id,haplotype,stage,analysis) |>
        filter(!(analysis=='quast' & metric=='Total length' & any(metric=='Total length (>= 0 bp)'))) |> ungroup()
    busco_totals <- measurements |> filter(analysis=='busco',metric=='total_busco') |>
        select(sample_id,haplotype,stage,busco_total=raw)
    measurements <- measurements |> left_join(busco_totals,by=c('sample_id','haplotype','stage')) |>
        mutate(measure=case_when(analysis=='quast' & metric %in% c('Total length','Total length (>= 0 bp)') ~ 'Length (Mb)',
            analysis=='quast' & metric=='N50' ~ 'N50 (Mb)', analysis=='busco' & metric=='complete' ~ 'BUSCO (pp)',
            analysis=='merqury' & metric=='qv' ~ 'QV',analysis=='merqury' & metric=='kmer_completeness' ~ 'K-mer completeness (pp)',
            analysis=='hic_contact' & metric=='retention_pct' ~ 'Hi-C retention (pp)',
            analysis=='hic_contact' & metric=='trans_to_cis' ~ 'Hi-C trans:cis',TRUE ~ NA_character_),
            value=case_when(measure %in% c('Length (Mb)','N50 (Mb)') ~ raw/1e6,
                measure=='BUSCO (pp)' & busco_total>0 ~ 100*raw/busco_total,
                measure=='BUSCO (pp)' ~ NA_real_,TRUE ~ raw)) |>
        filter(!is.na(measure),is.finite(value))
    if(anyDuplicated(measurements[c('sample_id','haplotype','stage','measure')])) stop('Duplicate stage-change measurement')
    changes <- measurements |> mutate(stage_rank=match(stage,c(stage_order,sort(setdiff(unique(stage),stage_order))))) |>
        arrange(sample_id,haplotype,measure,stage_rank) |> group_by(sample_id,haplotype,measure) |>
        mutate(previous_stage=lag(stage),previous_metric=lag(metric),previous_value=lag(value),change=value-previous_value) |>
        ungroup() |> filter(!is.na(previous_stage))
    write_tsv(changes |> select(sample_id,haplotype,previous_stage,stage,analysis,previous_metric,metric,measure,previous_value,value,change),file.path(out,'stage_changes.tsv'))
    stage_md <- c('# Changes across assembly stages','',
        '[Main report](../../assembly_report.md) · [All QC measurements](../../qc/assembly/assembly_qc_metrics.csv)','',
        'Values below are signed differences between consecutive available measurements for each metric. Missing stages are not imputed. Changes describe the measured outputs and do not establish their biological cause. pp means percentage points. N50 measures scaffold contiguity and may rise after a poor join; interpret it with completeness and read support.','',
        'Length changes use the recorded QUAST total-length metric and its sequence filters. Exact all-scaffold final lengths are provided in the main report inventory. Hi-C retention and trans:cis changes may also reflect different mapping/filtering contexts between stages; these are not pure measures of assembly improvement.','',
        'Stage keys: ctg.base = initial contigs; ctg.org = organelle-filtered; ctg.purged = duplicate-purged; ctg.polished = polished; ctg.cor = corrected contigs; ctg.deco = decontaminated contigs; scaf.base = first scaffolding; scaf.cor = corrected scaffolds; scaf2 = second scaffolding; chim.brk = post-cut; gap_fill = gap-filled; teloclip = telomere extension; final = finalized. Unknown stage keys are retained verbatim.','')
    for(sample in unique(changes$sample_id)) {
        display <- changes |> filter(sample_id==.env$sample) |> arrange(haplotype,stage_rank,measure) |>
            transmute(Haplotype=haplotype,Comparison=paste(previous_stage,'→',stage),measure,
                value=sprintf('%+.2f',change)) |> pivot_wider(names_from=measure,values_from=value,values_fill='Not measured')
        stage_md <- c(stage_md,paste0('## ',sample),'',md_table(display),'')
    }
    plots <- artifacts |> filter(type=='stage_qc')
    if(nrow(plots)>0) stage_md <- c(stage_md,'## Trend plots','',map_chr(plots$path,function(p) paste0('![',basename(p),'](../../',p,')')),'')
    writeLines(stage_md,file.path(out,'assembly_stage_changes.md'))
}

# Combined scaffold-length and chromosome-assignment view, per species.
paint_files <- chars(spec$paintings)
painting <- map_dfr(paint_files[grepl('\\.tsv$',paint_files)],function(path) optional_tsv(path) |>
    mutate(assembly=sub('\\.painting\\.tsv$','',basename(path)),
           length=as.numeric(length),start_0based=as.numeric(start_0based),end_exclusive=as.numeric(end_exclusive)))
expected_painting <- assemblies$id[assemblies$taxid %in% map_chr(elig$linear_references,function(x) as.character(x$taxid))]
if(isTRUE(spec$run$painting_enabled) &&
   !setequal(expected_painting,if(nrow(painting)>0) unique(painting$assembly) else character()))
    stop('Missing or unexpected final assembly painting outputs')
if(nrow(painting)>0) {
    checked <- painting |> distinct(assembly,sequence,length) |>
        full_join(names_data |> select(assembly,sequence=new_name,map_length=length),by=c('assembly','sequence'))
    painted_ids <- unique(painting$assembly)
    if(any(is.na(checked$length[checked$assembly %in% painted_ids])) ||
       any(is.na(checked$map_length[checked$assembly %in% painted_ids])) ||
       any(checked$length != checked$map_length,na.rm=TRUE)) stop('Painting and final name-map dictionaries differ')
    if(any(!is.finite(painting$start_0based) | painting$start_0based<0 |
           painting$end_exclusive>painting$length | painting$end_exclusive<=painting$start_0based)) stop('Invalid painting coordinates')
    write_tsv(painting,file.path(out,'painting_blocks.tsv'))
    chromosome_labels <- unique(painting$chromosome)
    chromosome_labels <- chromosome_labels[!chromosome_labels %in% c('Unassigned','Ambiguous')]
    chromosome_labels <- chromosome_labels[order(suppressWarnings(as.numeric(str_remove(chromosome_labels,'^chr'))),chromosome_labels,na.last=TRUE)]
    palette <- c(setNames(grDevices::hcl.colors(max(1,length(chromosome_labels)),'Dark 3')[seq_along(chromosome_labels)],chromosome_labels),
                 Unassigned='#D9D9D9',Ambiguous='#666666')
    scope <- tibble()
    md <- c(md,'## Chromosome structure','',
        'Bar length shows scaffold size. Scaffold labels come from the final harmonized name maps; colours come from alignment to one selected final assembly per species, rather than a combined cross-assembly consensus. Grey is unassigned and dark grey is conflicting assignment. Colours summarize qualifying mapping-chain spans, not base-level coverage or proof that a join is supported.','')
    for(taxid in unique(assemblies$taxid)) {
        members <- assemblies |> filter(.data$taxid==.env$taxid,id %in% painted_ids)
        if(nrow(members)==0) next
        inventory <- names_data |> filter(assembly %in% members$id)
        # Named chromosomes first; unnamed-only assemblies retain an individual-scaffold view.
        top <- inventory |> group_by(assembly) |>
            filter(if(any(str_detect(new_name,'^chr'))) str_detect(new_name,'^chr') else TRUE) |>
            slice_max(length,n=as.integer(spec$run$top_scaffolds),with_ties=FALSE) |> ungroup() |>
            mutate(first_chrom=str_remove(str_split_i(new_name,fixed('+'),1),'_[0-9]+$'),
                   chrom_number=suppressWarnings(as.numeric(str_extract(first_chrom,'[0-9]+$'))),
                   piece=suppressWarnings(as.numeric(str_extract(str_split_i(new_name,fixed('+'),1),'(?<=_)[0-9]+$')))) |>
            arrange(assembly,desc(chrom_number),desc(first_chrom),desc(piece),desc(new_name)) |>
            mutate(row=paste(assembly,new_name,sep='|'),row=factor(row,levels=unique(row)))
        unplaced <- inventory |> group_by(assembly) |>
            summarise(unplaced_pieces=sum(!str_detect(new_name,'^chr')),
                      unplaced_bp=sum(length[!str_detect(new_name,'^chr')]),.groups='drop') |>
            mutate(assembly=factor(assembly,levels=rev(members$id)),
                   text=paste0(sprintf('%.1f Mb',unplaced_bp/1e6),' · ',format(unplaced_pieces,big.mark=','),' pieces'))
        blocks <- painting |> inner_join(top |> select(assembly,new_name,row),by=c('assembly','sequence'='new_name'))
        refs <- keep(elig$linear_references,function(x) as.character(x$taxid)==taxid)
        reference <- if(length(refs)>0) refs[[1]]$id else 'Unavailable'
        panel_labels <- top |> group_by(assembly) |> summarise(shown=sum(length),.groups='drop') |>
            left_join(assembly_lengths |> select(assembly,assembly_bp),by='assembly') |>
            mutate(label=paste0(assembly,'\nIndividual bars: ',sprintf('%.1f%%',100*shown/assembly_bp),' of assembly'))
        panel_labels <- setNames(panel_labels$label,panel_labels$assembly)
        blocks <- blocks |> mutate(panel=factor(assembly,levels=members$id))
        ncols <- max(1,min(4,n_distinct(members$haplotype)))
        p <- ggplot(blocks,aes(x=(start_0based+end_exclusive)/2e6,y=row,
                              width=(end_exclusive-start_0based)/1e6,fill=chromosome)) +
            geom_tile(height=.7) + facet_wrap(~panel,ncol=ncols,scales='free_y',labeller=labeller(panel=as_labeller(panel_labels))) +
            scale_y_discrete(labels=function(x) sub('^[^|]+\\|','',x)) +
            scale_x_continuous(limits=c(0,max(top$length)/1e6),expand=expansion(mult=c(0,.02))) +
            scale_fill_manual(values=palette,breaks=c(chromosome_labels,'Unassigned','Ambiguous'),name='Reference chromosome',drop=FALSE) +
            labs(title=paste0('Scaffold lengths and chromosome assignments · ',species_label(taxid)),
                 subtitle=paste0('Up to ',spec$run$top_scaffolds,' longest chromosome-named scaffolds per assembly, in chromosome order'),
                 x='Position along scaffold (Mb)',y=NULL,
                 caption=paste0('Painting source: ',reference,' (single final assembly). Scaffold labels: final harmonized name maps.\n',
                     'Reference panel uses identity labels, not independent alignment evidence. Unnamed-only assemblies show their longest pieces.')) +
            theme_minimal(base_size=10) + theme(panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
                 plot.background=element_rect(fill='white',colour=NA),legend.position='bottom',axis.text.y=element_text(size=7),
                 strip.text=element_text(size=9)) + guides(fill=guide_legend(nrow=2,byrow=TRUE))
        summary_plot <- ggplot(unplaced,aes(unplaced_bp/1e6,assembly)) +
            geom_col(fill='#44546A',width=.55) + geom_text(aes(label=text),hjust=-.05,size=2.7) +
            scale_x_continuous(expand=expansion(mult=c(0,.4))) +
            labs(title='All unplaced sequence · separate length scale',
                 subtitle='Each bar sums independent pieces; it is not a joined scaffold',
                 x='Summed unplaced length (Mb)',y=NULL) +
            theme_minimal(base_size=10) + theme(panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
                plot.background=element_rect(fill='white',colour=NA),axis.text.y=element_text(size=8))
        panel_height <- max(4,ceiling(nrow(members)/ncols)*3)
        summary_height <- max(2,1.2+nrow(members)*.26)
        combined <- p / summary_plot + plot_layout(heights=c(panel_height,summary_height))
        filename <- paste0('chromosome_structure_',taxid)
        ggsave(file.path(out,paste0(filename,'.png')),combined,width=14,height=panel_height+summary_height,dpi=160,bg='white')
        ggsave(file.path(out,paste0(filename,'.pdf')),combined,width=14,height=panel_height+summary_height,bg='white',device=cairo_pdf)
        md <- c(md,paste0('**',species_label(taxid),' · Painting reference: ',reference,'.** Other assemblies are aligned to this final assembly. ',
                    'Its own panel uses identity labels, not independent alignment evidence. This painting does not summarize the multi-assembly evidence used to assess chimeras.'),'',
                paste0('![Scaffold lengths, chromosome assignments and unplaced totals](reports/main/',filename,'.png)'),'',
                link('PDF',paste0('reports/main/',filename,'.pdf')),'',
                'The lower bars sum **all unplaced pieces** on a separate scale. Unplaced is a final naming category; unassigned (grey above) means no painting label passed the filters. The views are not additive: in unnamed-only assemblies the lower total includes the individual pieces shown above.','')
        total <- inventory |> group_by(assembly) |>
            summarise(total_scaffolds=n(),total_bp=sum(length),.groups='drop')
        scope <- bind_rows(scope,top |> group_by(assembly) |> summarise(shown_scaffolds=n(),shown_bp=sum(length),.groups='drop') |>
            left_join(total,by='assembly') |> left_join(unplaced |> mutate(assembly=as.character(assembly)) |>
                select(assembly,unplaced_pieces,unplaced_bp),by='assembly') |> mutate(taxid=.env$taxid))
        write_tsv(top |> select(assembly,new_name,length,row),file.path(out,paste0('display_rows_',taxid,'.tsv')))
    }
    if(anyDuplicated(scope$assembly)) stop('Assembly displayed in multiple species panels')
    write_tsv(scope,file.path(out,'display_scope.tsv'))
    scope_display <- scope |> transmute(Assembly=assembly,`Scaffolds shown / total`=paste(shown_scaffolds,total_scaffolds,sep=' / '),
                                      `Sequence shown (%)`=sprintf('%.1f',100*shown_bp/total_bp))
    md <- c(md,'<details><summary>Figure scope and underlying evidence</summary>','',md_table(scope_display),'',
        link('Complete scaffold inventory','reports/main/scaffold_inventory.tsv'),' · ',
        link('Painting blocks','reports/main/painting_blocks.tsv'),'',
        'Raw mapping chains, filters, coordinate metadata and software versions are retained under reports/chromosome_painting/. These display assignments do not alter assembly names or chimera decisions.','','</details>','')
}

# Reconcile a detected event with recorded choices and verified applied actions.
review <- optional_tsv(spec$decisions)
events <- tibble()
if(nrow(review)>0) {
    review <- review |> mutate(event_id=if_else(is.na(event_id) | event_id=='',id,event_id))
    verified <- map_dfr(chars(spec$verifications),function(path) {
        if(grepl('^NO_',basename(path))) return(tibble())
        v <- fromJSON(path)
        if(!isTRUE(v$exact_parent_reconstruction) || length(v$selected_actions)==0) return(tibble())
        as_tibble(v$selected_actions) |> select(assembly,id,cut_bp,assessment_sha256)
    })
    if(nrow(verified)>0) review <- review |> left_join(verified |> distinct() |> mutate(verified=TRUE),by=c('assembly','id','cut_bp','assessment_sha256'))
    else review$verified <- FALSE
    events <- review |> group_by(assembly,event_id) |> group_modify(function(rows,key) {
        cuts <- rows |> filter(selected=='YES',review_disposition=='CUT')
        retained <- rows |> filter(review_disposition=='RETAIN')
        deferred <- rows |> filter(review_disposition=='DEFER')
        selected <- if(nrow(cuts)>0) cuts[1,] else if(nrow(retained)>0) retained[1,] else if(nrow(deferred)>0) deferred[1,] else rows[1,]
        state <- if(nrow(cuts)>0) {
            if(all(cuts$verified %in% TRUE)) 'Cut applied and verified' else 'Cut selected; verification needed'
        } else if(nrow(retained)>0) 'Retained by review' else if(nrow(deferred)>0) 'Review deferred' else 'Review pending'
        tibble(chromosomes=selected$chromosomes,scaffold=selected$scaffold,status=state,
               cut_bp=if(nrow(cuts)>0) paste(cuts$cut_bp,collapse=', ') else '',reason=selected$reason)
    }) |> ungroup()
    write_tsv(events,file.path(out,'assembly_decisions.tsv'))
    pending <- events |> filter(!status %in% c('Cut applied and verified','Retained by review'))
    if(nrow(pending)>0) attention <- bind_rows(attention,pending |> transmute(Area='Assembly join decisions',Item=paste(assembly,event_id),Issue=status,
        Details=link(chromosomes,'reports/main/decisions.md')))
    md <- c(md,'## Assembly decisions and changes','',
        'These decisions describe joins in the original pre-finishing assemblies. A retained join is a recorded review choice, not proof of a biological rearrangement.','',
        md_table(events |> transmute(Assembly=assembly,Event=event_id,Transition=chromosomes,Outcome=status,
                                     `Original cut location`=if_else(cut_bp=='','—',paste0(scaffold,':',
                                         map_chr(str_split(cut_bp,', '),function(x) paste(format(as.numeric(x),big.mark=',',scientific=FALSE,trim=TRUE),collapse='; ')))))),'',
        link('Decision reasons and evidence','reports/main/decisions.md'),' · ',
        link('Chimera assessment across all assemblies','assembly/chimeras/README.md'),'')
    details <- c('# Assembly decisions','',md_table(events |> select(Assembly=assembly,Event=event_id,Status=status)),'')
    for(i in seq_len(nrow(events))) {
        e <- events[i,]
        details <- c(details,paste0('## ',e$assembly,' · ',e$event_id),'',paste0('**',e$status,'** · ',e$chromosomes),'',
            e$reason,'',link('Original evidence',paste0('../../assembly/chimeras/assemblies/',e$assembly,'/report.md#',tolower(e$event_id))),'')
    }
    writeLines(details,file.path(out,'decisions.md'))
} else if(isTRUE(spec$run$chimera_enabled)) {
    md <- c(md,'## Chimera evidence','',link('Assessment and manual review','assembly/chimeras/README.md'),'')
}

expected <- bind_rows(lapply(spec$expected,as_tibble))
if(nrow(expected)>0) {
    observed <- assemblies |> count(sample,name='observed')
    outcomes <- expected |> left_join(observed,by='sample') |> mutate(observed=replace_na(observed,0L),
        outcome=case_when(status!='accepted' ~ status,observed==as.integer(expected) ~ 'Finalized',TRUE ~ 'Incomplete'))
    outcomes <- outcomes |> arrange(sample)
    exceptions <- outcomes |> filter(outcome!='Finalized')
    incomplete <- exceptions |> filter(status=='accepted')
    if(nrow(incomplete)>0) attention <- bind_rows(attention,incomplete |> transmute(Area='Assembly completion',Item=sample,Issue='Assembly incomplete',Details=reason))
    md <- c(md,'## Assembly outcomes and files','',
        paste0(nrow(outcomes),' supplied individuals: ',sum(outcomes$outcome=='Finalized'),' finalized',
               if(nrow(exceptions)>0) paste0('; ',nrow(exceptions),' skipped or incomplete') else '', '.'),'')
    if(nrow(exceptions)>0) md <- c(md,md_table(exceptions |> transmute(Individual=sample,Outcome=outcome,
        Reason=if_else(str_detect(reason,'Hi-C-only'),'Hi-C reads only; no assembly reads supplied.',reason))),'')
    md <- c(md,'<details><summary>All sample outcomes</summary>','',md_table(outcomes |> transmute(Individual=sample,Outcome=outcome,
        `Finalized / expected`=if_else(as.integer(expected)>0,paste(observed,expected,sep=' / '),'Not applicable'))),'','</details>','')
}
downloads <- assemblies |> mutate(files=paste0(link('FASTA',paste0('assembly/final/',basename(fasta))), ' · ',
    link('Name map',paste0('assembly/final/',id,'.name_map.tsv')))) |>
    select(Individual=sample,haplotype,files) |> pivot_wider(names_from=haplotype,values_from=files,values_fill='—')
mitogenomes <- artifacts |> filter(type=='organelle_annotation',str_detect(path,'\\.(gb|gbk|genbank)$')) |>
    group_by(id) |> summarise(Mitogenome=paste(link('GenBank',path),collapse=' · '),.groups='drop')
if(nrow(mitogenomes)>0) downloads <- downloads |> left_join(mitogenomes,by=c('Individual'='id')) |>
    mutate(Mitogenome=replace_na(Mitogenome,'Not available'))
downloads_md <- c('## Final assemblies, chromosome name maps and mitogenomes','',md_table(downloads),'')

# Embed existing final figures directly; their evidence and renderers are unchanged.
snails <- artifacts |> filter(type=='snail') |> inner_join(assemblies |> select(id,sample,haplotype),by='id')
mito_maps <- artifacts |> filter(type=='organelle',str_detect(path,'_mito_circular\\.png$'))
thumbnail <- function(path,label) paste0('[<img src="',str_replace_all(path,' ','%20'),'" width="220" alt="',label,'" />](',str_replace_all(path,' ','%20'),')')
if(nrow(snails)>0 || nrow(mito_maps)>0) {
    figures <- tibble(Individual=sort(unique(c(snails$sample,mito_maps$id))))
    if(nrow(snails)>0) {
        if(anyDuplicated(snails$id)) stop('Multiple final snail plots for an assembly')
        figure_snails <- snails |> transmute(Individual=sample,haplotype,
            Figure=thumbnail(path,paste0(id,' final snail plot'))) |>
            pivot_wider(names_from=haplotype,values_from=Figure,names_glue='{haplotype} snail',values_fill='Not available')
        figures <- figures |> left_join(figure_snails,by='Individual')
    }
    if(nrow(mito_maps)>0) {
        figure_mito <- mito_maps |> group_by(id) |> summarise(`Mitogenome map`=paste(thumbnail(path,paste0(id,' mitogenome map')),collapse=' '),.groups='drop')
        figures <- figures |> left_join(figure_mito,by=c('Individual'='id'))
    }
    figures <- figures |> mutate(across(everything(),function(x) replace_na(x,'Not available')))
    md <- c(md,'## Assembly and mitogenome plots','','Final assembly snail plots summarize scaffold lengths, base composition and BUSCO results. Mitogenome maps show annotated feature positions on the organelle sequence; a circular drawing alone does not establish circularity. Colours follow each figure legend. Select an image to open the full figure.','',md_table(figures),'')
}
graph <- artifacts |> filter(type=='pangenome')
if(nrow(graph)>0) {
    md <- c(md,'## Pangenome','')
    for(taxid in unique(graph$id)) {
        rows <- graph |> filter(id==.env$taxid)
        review_notes <- rows |> filter(!is.na(note),str_detect(note,'REVIEW|FAIL|ERROR'))
        review_rows <- artifacts |> filter(type=='graph_review_rows',id==.env$taxid)
        review_audit <- artifacts |> filter(type=='graph_review_audit',id==.env$taxid)
        flagged <- if(nrow(review_rows)>0) read_tsv(I(review_rows$content[[1]]),show_col_types=FALSE) else tibble()
        if(nrow(flagged)>0) {
            flagged <- flagged |> filter(retention_review=='yes')
            cutoff <- if(nrow(review_audit)>0) fromJSON(review_audit$content[[1]])$min_retained_fraction_screen else NULL
            cutoff_text <- if(length(cutoff)==1) paste0(' (review screen: <',sprintf('%.0f%%',100*cutoff),')') else ' (flagged by source audit)'
            attention <- bind_rows(attention,flagged |> transmute(Area='Pangenome graph',Item=paste(assembly,sequence,sep=' / '),
                Issue=paste0('Graph path-length screen: ',if_else(is.na(retained_path_fraction),'length unavailable',sprintf('%.1f%% of input sequence length',100*retained_path_fraction)),cutoff_text),
                Details=link('Review graph representation and mapping evidence',paste0('pangenome/',taxid,'/representation/source_report.md'))))
        } else if(nrow(review_notes)>0) attention <- bind_rows(attention,tibble(Area='Pangenome graph',Item=species_label(taxid),
            Issue='Graph source audit reports a finding requiring review',Details=paste0('Inspect the listed source/graph representation discrepancies. ',
                link('Source-audit findings',review_notes$path[[1]]))))
        ref <- keep(elig$linear_references,function(x) as.character(x$taxid)==taxid)
        ref_name <- if(length(ref)) ref[[1]]$id else 'Not selected'
        md <- c(md,paste0('**',species_label(taxid),'** · Reference: ',ref_name,' · ',sum(assemblies$taxid==taxid & assemblies$graph_member),' contributing assemblies.'),'')
        summaries <- artifacts |> filter(type=='pangenome_summary',id==.env$taxid)
        sharing <- tibble()
        if(nrow(summaries)>0) {
            sharing <- read_tsv(I(summaries$content[[1]]),show_col_types=FALSE)
            sharing <- sharing |> filter(unit=='individual',class %in% c('core','softcore','shell','cloud')) |>
                mutate(class=factor(class,levels=c('core','softcore','shell','cloud')),
                    lower=case_when(class=='core' ~ core_min,class=='softcore' ~ softcore_min,class=='shell' ~ shell_min,TRUE ~ 1),
                    upper=case_when(class=='core' ~ denominator,class=='softcore' ~ core_min-1,
                        class=='shell' ~ pmin(softcore_min,core_min)-1,TRUE ~ pmin(shell_min,softcore_min,core_min)-1),
                    frequency=pmap_chr(list(lower,upper,denominator),function(lo,hi,n) {
                        if(lo>hi) return('No separate frequency tier')
                        if(lo==n && hi==n) return(paste0('All ',n,' individuals'))
                        paste0(if(lo==hi) lo else paste0(lo,'–',hi),' of ',n,' individuals')
                    })) |> arrange(class)
            if(nrow(sharing)>0) {
                write_tsv(sharing,file.path(out,paste0('pangenome_sharing_',taxid,'.tsv')))
                md <- c(md,md_table(sharing |> transmute(`Sharing across individuals`=frequency,
                    `Distinct graph sequence (Mb)`=sprintf('%.1f',bp/1e6),
                    `Graph sequence (%)`=sprintf('%.1f',percent_represented_graph_bp))),'',
                    'Presence means found in at least one included haplotype per individual. These are distinct represented graph sequences, not assembly lengths or proven biological novelty; assembly quality and graph representation affect apparent absence.','')
            }
        }
        # A single entry point retains every enabled biological and technical report.
        hub <- c(paste0('# Pangenome report · ',species_label(taxid)),'',
            paste0('Reference: **',ref_name,'**. See the linked analyses for their measurement definitions and limitations.'),'')
        if(nrow(sharing)>0) hub <- c(hub,'## Sharing across individuals','',
            md_table(sharing |> transmute(`Category`=as.character(class),`Found in`=frequency,
                `Distinct graph sequence (Mb)`=sprintf('%.1f',bp/1e6),
                `Graph sequence (%)`=sprintf('%.1f',percent_represented_graph_bp))),'',
            'An individual carries sequence when it is represented in either included haplotype. These categories describe graph sequence, not proven biological novelty. Apparent absence can reflect assembly quality or graph representation. Chromosome-level unions may overlap and must not be summed.','')
        biological <- basename(rows$path) %in% c('sharing_report.md','sharing_figures.md','regional_report.md','regional_figures.md','comparison_report.md')
        for(section in c('Biological results','Technical validation and provenance')) {
            selected <- rows[if(section=='Biological results') biological else !biological,]
            if(nrow(selected)>0) hub <- c(hub,paste0('## ',section),'',
                if(section!='Biological results') rows$note[!is.na(rows$note) & rows$note!=''] else character(),'',
                paste0('- ',map_chr(selected$path,function(p) link(report_label(p),paste0('../../',p)))),'')
        }
        hub_name <- paste0('pangenome_',taxid,'.md')
        writeLines(hub,file.path(out,hub_name))
        md <- c(md,link('Complete pangenome report: biological results and technical validation',paste0('reports/main/',hub_name)),'')
    }
}
psmc <- artifacts |> filter(type=='psmc')
if(nrow(psmc)>0) {
    md <- c(md,'## Demographic history (PSMC)','',link('PSMC histories, input scope and callable-region evidence',psmc$path[[1]]),'')
    records <- artifacts |> filter(type=='psmc_summary')
    if(nrow(records)!=1) stop('PSMC report requires exactly one current-run summary')
    tab <- read_tsv(I(records$content[[1]]),col_types=cols(.default=col_character()),show_col_types=FALSE)
    md <- c(md,md_table(tab |> transmute(Individual=sample,Mode=mode,Result=status,
        `Callable bp retained for fitting`=retained_callable_bp,`Retained scaffolds`=retained_scaffolds,
        Scaling=scaling,Details=if_else(!is.na(report) & nzchar(report),link('Report',report),reason))),
        '', 'Scaling values are explicit inputs. Unscaled histories do not assume a mutation rate or generation time. Long-scaffold scope may include sex-linked sequence.','')
}
md <- c(md,'## Detailed evidence','',link('All available figures and reports','reports/main/evidence_index.md'),'')
if(any(artifacts$type=='chimera_registry')) md <- c(md,link('Chimera screening scope and coverage','reports/main/assessment_coverage.md'),'')
if(stage_changes_present) md <- c(md,link('What changed across assembly stages','reports/main/assembly_stage_changes.md'),'')
for(type in c('stage_qc','contact_map','telomere_extension')) {
    rows <- artifacts |> filter(.data$type==.env$type)
    if(nrow(rows)>0) {
        title <- c(stage_qc='QC across assembly stages',contact_map='Final Hi-C contact maps',telomere_extension='Telomere extension results')[[type]]
        md <- c(md,paste0('- ',link(title,paste0('reports/main/evidence_index.md#',str_replace_all(type,'_','-')))))
    }
}
organelle <- artifacts |> filter(type=='organelle_stats')
recovery <- artifacts |> filter(type=='organelle_recovery')
if(nrow(organelle)>0 || nrow(recovery)>0) md <- c(md,'','## Organelle results','')
if(nrow(organelle)>0) {
    org_table <- map_dfr(organelle$content,function(text) read_tsv(I(text),col_types=cols(.default=col_character()),show_col_types=FALSE))
    for(column in c('sample_id','mitogenome_length','circular')) if(!column %in% names(org_table)) org_table[[column]] <- NA_character_
    org_display <- org_table |> transmute(Individual=sample_id,`Length (bp)`=mitogenome_length,
        Recovery=case_when(tolower(circular) %in% c('yes','true') ~ 'Circular',
                           tolower(circular) %in% c('no','false') ~ 'Not circular',TRUE ~ 'Not measured'))
    md <- c(md,md_table(org_display),'',
        link('Organelle measurements','reports/main/evidence_index.md#organelle-stats'),'')
    for(kind in c('organelle','organelle_annotation')) {
        if(any(artifacts$type==kind)) md <- c(md,link(if(kind=='organelle') 'Organelle maps' else 'Organelle annotations',
            paste0('reports/main/evidence_index.md#',str_replace_all(kind,'_','-'))),'')
    }
}
annotation_counts <- artifacts |> filter(type=='annotation_counts')
if(nrow(annotation_counts)>0) {
    counts <- map_dfr(annotation_counts$content,function(text) read_tsv(I(text),show_col_types=FALSE))
    write_tsv(counts,file.path(out,'annotation_counts.tsv'))
    audit_md <- c('# Mitogenome annotation feature counts','',
        '[Main report](../../assembly_report.md)','',
        'Counts use unique annotated locations within each GenBank record. A gene feature and a CDS/tRNA/rRNA feature at the same location count as one locus. Protein-coding, tRNA and rRNA columns count their own feature types. These are annotation counts, not a completeness assessment; unusual or missing features remain visible rather than being forced to an expected total.','',
        md_table(counts |> transmute(Individual=sample_id,`Annotation length (bp)`=annotation_length_bp,
            `Annotated loci`=annotated_loci,`Protein-coding loci`=protein_coding_loci,`tRNA loci`=trna_loci,`rRNA loci`=rrna_loci,`Other gene loci`=other_gene_loci)),'',
        'Older source summary fields named gene_count/trna_count/rrna_count were text-occurrence tallies. They are not used for the biological feature counts above. The GenBank annotations themselves are unchanged.','',
        '## Source annotations and locus audit','')
    audit_rows <- artifacts |> filter(type %in% c('organelle_annotation','annotation_loci')) |>
        transmute(Individual=id,File=link(if_else(type=='annotation_loci','Locus audit','GenBank'),paste0('../../',path)))
    audit_md <- c(audit_md,md_table(audit_rows),'')
    writeLines(audit_md,file.path(out,'organelle_annotations.md'))
    md <- c(md,link('Audited mitochondrial feature counts','reports/main/organelle_annotations.md'),'')
}
if(nrow(recovery)>0) {
    recovery_table <- map_dfr(seq_len(nrow(recovery)),function(i) {
        record <- recovery[i,]
        values <- read_tsv(I(record$content),col_types=cols(.default=col_character()),show_col_types=FALSE)
        if(!'status' %in% names(values)) stop('Organelle recovery record has no status')
        tibble(Individual=record$id,Organelle=record$organelle,Status=str_replace_all(values$status,'_',' '))
    })
    if(nrow(organelle)==0 || any(recovery_table$Status=='no recovery')) md <- c(md,md_table(recovery_table),'')
    md <- c(md,link('All recovery outcomes','reports/main/evidence_index.md#organelle-recovery'),'')
}
legacy <- artifacts |> filter(type=='assembly_details')
if(nrow(legacy)>0) md <- c(md,link('Detailed assembly, organelle, telomere and stage QC report',legacy$path[[1]]),'')
index <- c('# Evidence and detailed reports','',
    'Only artifacts supplied by this invocation are listed. Coordinate stages and measurement definitions remain with their evidence reports.','')
themes <- list('Assembly quality'=c('interactive_qc','snail','stage_qc'),
    'Chromosome structure and comparisons'=c('dotplot','ribbon','contact_map','telomere_extension','chimera_registry'),
    'Organelle assemblies'=c('organelle','organelle_annotation','annotation_counts','annotation_loci','organelle_stats','organelle_recovery'),
    'Pangenome'=c('pangenome','pangenome_summary','graph_review_rows','graph_review_audit'),
    'Demographic history'=c('psmc','psmc_summary'),
    'Methods and provenance'=c('provenance'),'Additional reports'=c('assembly_details'))
themes[['Additional reports']] <- c(themes[['Additional reports']],setdiff(unique(artifacts$type),unlist(themes)))
for(theme in names(themes)) {
    types <- themes[[theme]][themes[[theme]] %in% artifacts$type]
    if(!length(types)) next
    index <- c(index,paste0('## ',theme),'')
    for(type in types) {
    rows <- artifacts |> filter(.data$type==.env$type)
    heading <- c(stage_qc='QC across assembly stages',contact_map='Hi-C contact maps',organelle='Organelle maps',
        organelle_annotation='Organelle annotations',organelle_stats='Organelle measurements',organelle_recovery='Organelle recovery',
        telomere_extension='Telomere extension',pangenome='Pangenome analyses',pangenome_summary='Pangenome source measurements',
        annotation_counts='Annotation feature counts',annotation_loci='Annotation locus audit',chimera_registry='Chimera assessment coverage',provenance='Software and run provenance',
        graph_review_rows='Source sequence representation measurements',graph_review_audit='Graph source-audit settings and outcome',
        interactive_qc='Interactive QC rankings',snail='Final assembly snail plots',dotplot='Dotplots',ribbon='Ribbon plots',assembly_details='Detailed assembly report')
    title <- if(type %in% names(heading)) heading[[type]] else str_to_sentence(str_replace_all(type,'_',' '))
    index <- c(index,paste0('<a id="',str_replace_all(type,'_','-'),'"></a>'),'',paste0('### ',title),'')
    if(type=='contact_map') {
        grouped <- rows |> mutate(resolution=suppressWarnings(as.numeric(str_match(basename(path),'_([0-9]+)bp_')[,2]))) |>
            arrange(id,resolution,path) |> group_by(id) |> summarise(Reports=paste(link(
                if_else(is.na(resolution),basename(path),paste0(format(resolution/1000,trim=TRUE),' kb')),paste0('../../',path)),collapse=' · '),.groups='drop')
        index <- c(index,md_table(grouped |> rename(Assembly=id)),'')
    } else index <- c(index,md_table(rows |> transmute(Item=id,Report=link(report_label(path),paste0('../../',path)))),'')
    }
}
writeLines(index,file.path(out,'evidence_index.md'))
md <- c(md,'## Methods and reproducibility','',
    paste0('Run: **',spec$run$run_name,'** · Nextflow ',spec$run$nextflow,' · Profile ',spec$run$profile,'.'),'',
    if(is.character(spec$run$commit) && length(spec$run$commit)==1 && !is.na(spec$run$commit) && nzchar(spec$run$commit)) paste0('Pipeline commit: `',spec$run$commit,'`.') else
        'Pipeline executed from a working tree; source fingerprints are recorded in run provenance.','',
    link('Software versions','provenance/software_versions.tsv'),' · ',link('Version metadata issues','provenance/software_version_issues.tsv'),' · ',
    link('Run parameters and source fingerprints','provenance/run_provenance.json'),' · ',link('Environment snapshots','provenance/runtime_environments.json'),' · ',
    link('Report renderer versions','reports/main/reporting_versions.tsv'),'')
# Report expectations use the enabled controls AND the sample's actual inputs.
checks <- tibble(Feature=character(),Entity=character(),Status=character(),Reason=character())
check_output <- function(feature,entity,enabled,present,applicable=TRUE,reason='') {
    checks <<- bind_rows(checks,tibble(Feature=feature,Entity=entity,
        Status=if(!enabled) 'DISABLED' else if(!applicable) 'NOT_APPLICABLE' else if(present) 'PRESENT' else 'MISSING',Reason=reason))
}
has_artifact <- function(type,id) any(artifacts$type==type & artifacts$id==id)
registry <- artifacts |> filter(type=='chimera_registry')
if(nrow(registry)>0) {
    coverage <- read_tsv(I(registry$content[[1]]),col_types=cols(.default=col_character()),show_col_types=FALSE)
    coverage_display <- coverage |> transmute(Assembly=assembly,
        Screen=case_when(assessment_status=='assessed' ~ 'Assessed',assessment_status=='no_chromosome_scope' ~ 'Not applied: no confirmed chromosome-scale scope',TRUE ~ 'Assessment unavailable'),
        Explanation=if_else(assessment_status=='no_chromosome_scope',
            'No confirmed chromosome-scale scaffold set was available. This screen limitation does not identify a misjoin or require a cut.',result))
    writeLines(c('# Chimera screening scope and coverage','',
        'Chromosome names can identify sequence associations in fragments without establishing chromosome-scale scaffolds. Assemblies without confirmed chromosome-scale scope remain here as a screen limitation, rather than an actionable warning. Original detection summaries below are distinct from completed review decisions.','',
        md_table(coverage_display),'',link('Recorded decisions','../../assembly_report.md#assembly-decisions-and-changes')),
        file.path(out,'assessment_coverage.md'))
    # Absence of chromosome-scale scope is a screen limitation, not a repair request.
    unresolved <- coverage |> filter(!assessment_status %in% c('assessed','no_chromosome_scope'))
    if(nrow(unresolved)>0) attention <- bind_rows(attention,unresolved |> transmute(Area='Chimera assessment',Item=assembly,
        Issue=paste0('Chimera assessment incomplete: ',result),Details=link('Assessment scope','assembly/chimeras/README.md')))
}
if(isTRUE(spec$run$contract_enabled)) {
    qc_enabled <- !is.null(spec$run$qc_mode) && spec$run$qc_mode!='none'
    for(i in seq_len(nrow(assemblies))) {
        a <- assemblies[i,];cap <- sample_capabilities |> filter(sample==a$sample)
        has_hifi <- nrow(cap)>0 && isTRUE(cap$hifi[[1]])
        has_hic <- nrow(cap)>0 && isTRUE(cap$hic[[1]])
        has_reads <- has_hifi || (nrow(cap)>0 && isTRUE(cap$shortread[[1]]))
        final_qc <- if(nrow(qc)>0) qc |> filter(sample_id==a$sample,stage=='final') else tibble()
        slot <- if(a$haplotype=='primary') 'hap1' else a$haplotype
        measured <- function(analysis_name,metric_name,allow_na=FALSE) {
            if(nrow(final_qc)==0 || !slot %in% names(final_qc)) return(FALSE)
            values <- final_qc[[slot]][final_qc$analysis==analysis_name & final_qc$metric==metric_name]
            length(values)>0 && (allow_na || any(!is.na(suppressWarnings(as.numeric(values)))))
        }
        check_output('Final QC',a$id,qc_enabled,measured('quast','N50'))
        check_output('Final BUSCO',a$id,qc_enabled,measured('busco','complete') && measured('busco','total_busco'))
        check_output('Final k-mer QC',a$id,qc_enabled,measured('merqury','qv') && measured('merqury','kmer_completeness'),has_reads,'Requires QC reads and their k-mer database')
        check_output('Final snail plot',a$id,qc_enabled,has_artifact('snail',a$id))
        check_output('Final mapped-read depth',a$id,qc_enabled,measured('mapped_reads','coverage'),has_reads,'Requires QC reads')
        check_output('Final Hi-C metrics',a$id,qc_enabled,all(map_lgl(c('retention_pct','trans_to_cis','cis_pairs','trans_pairs'),function(m) measured('hic_contact',m,allow_na=m %in% c('retention_pct','trans_to_cis')))),has_hic,'Requires Hi-C reads; undefined ratios are not treated as missing output files')
        for(resolution in chars(str_split(spec$run$contact_resolutions %||% '1000000,500000,250000,100000',','))) {
            found <- any(artifacts$type=='contact_map' & artifacts$id==a$id & str_detect(artifacts$path,fixed(paste0('_',resolution,'bp_'))))
            check_output(paste0('Contact map ',resolution,' bp'),a$id,qc_enabled && isTRUE(spec$run$contact_maps_enabled),found,has_hic,'Requires final QC and Hi-C reads')
        }
        check_output('Chimera assessment',a$id,isTRUE(spec$run$chimera_enabled),nrow(registry)>0 && a$id %in% coverage$assembly)
    }
    pairs <- crossing(left=assemblies$id,right=assemblies$id) |> filter(left<right) |>
        filter(assemblies$taxid[match(left,assemblies$id)]==assemblies$taxid[match(right,assemblies$id)])
    if((spec$run$comparison_mode %||% 'within_sample')=='within_sample') pairs <- pairs |>
        filter(assemblies$sample[match(left,assemblies$id)]==assemblies$sample[match(right,assemblies$id)])
    if(nrow(pairs)==0) check_output('Pairwise figures','cohort',isTRUE(spec$run$comparisons_enabled),FALSE,FALSE,'No eligible pairs')
    for(i in seq_len(nrow(pairs))) for(type in c('dotplot','ribbon')) {
        id <- paste(pairs$left[i],pairs$right[i],sep=' / ')
        check_output(type,id,isTRUE(spec$run$comparisons_enabled),has_artifact(type,id))
    }
    for(id in unique(artifacts$id[artifacts$type=='organelle_annotation'])) {
        check_output('Annotation feature audit',id,TRUE,has_artifact('annotation_counts',id) && has_artifact('annotation_loci',id))
    }
    for(i in which(artifacts$type=='organelle_recovery')) {
        recovery_status <- read_tsv(I(artifacts$content[[i]]),show_col_types=FALSE)
        recovered <- any(recovery_status$status %in% c('recovered','circular','linear'))
        check_output('Recovered mitogenome GenBank',artifacts$id[[i]],TRUE,has_artifact('organelle_annotation',artifacts$id[[i]]),
            recovered,'Requires a recovered mitogenome; no recovery remains an explicit outcome')
    }
    for(cohort in elig$cohorts) {
        id <- as.character(cohort$taxid)
        for(name in c('pangenome_build_report.md','pangenome_input_report.md',if(!isTRUE(spec$run$graph_validate_only)) 'source_report.md'))
            check_output(name,id,isTRUE(spec$run$graph_enabled),any(artifacts$type=='pangenome' & artifacts$id==id & basename(artifacts$path)==name),isTRUE(cohort$ready),'Requires an admitted graph cohort')
    }
} else check_output('Run-control contract','run',TRUE,FALSE,FALSE,'Original run controls unavailable; output expectations not inferred')
check_output('PSMC report','run',isTRUE(spec$run$psmc_enabled),nrow(psmc)>0)
write_tsv(checks,file.path(out,'report_completeness.tsv'))
missing <- checks |> filter(Status=='MISSING')
if(nrow(missing)>0) attention <- bind_rows(attention,missing |> transmute(Area='Report completeness',Item=Entity,Issue=paste('Expected output missing:',Feature),Details=link('Completeness checks','reports/main/report_completeness.md')))
writeLines(c('# Report completeness','',
    'Expected outputs are evaluated using run controls, eligible assemblies and available read types. Disabled and inapplicable analyses are distinguished from missing expected outputs. This checks report inputs; published files and local links are checked after Nextflow finishes publishing.','',md_table(checks)),file.path(out,'report_completeness.md'))
md <- c(md,link('Report completeness checks','reports/main/report_completeness.md'),'')
if(nrow(attention)>0) {
    attention <- attention |> distinct()
    md <- append(md,c('## Needs attention','',
        if(any(attention$Area=='Pangenome graph')) 'Pangenome graph entries concern how the graph represents the input assemblies. Review whether shorter paths reflect graph representation or sequence omitted from the graph. These checks do not identify new assembly misjoins or recommend assembly edits. Path-length ratios are not measurements of missing DNA.' else character(),'',md_table(attention),''),after=4)
}
write_tsv(attention,file.path(out,'needs_attention.tsv'))
first_section <- which(str_detect(md,'^## '))[1]
md <- append(md,downloads_md,after=if(is.na(first_section)) length(md) else first_section-1)
sections <- md[str_detect(md,'^## ')]
anchor <- function(h) h |> str_remove('^## ') |> str_to_lower() |> str_remove_all('[^a-z0-9 -]') |> str_replace_all(' +','-')
toc <- c('## Contents','',map_chr(sections,function(h) paste0('- ',link(str_remove(h,'^## '),paste0('#',anchor(h))))),'')
md <- append(md,toc,after=4)
writeLines(md,'assembly_report.md')
write_json(list(format_version=1,run=spec$run,assemblies=assemblies,
    qc_present=nrow(qc)>0,painting_required=isTRUE(spec$run$painting_enabled),
    expected_painting=as.list(expected_painting),painted_assemblies=as.list(if(nrow(painting)>0) unique(painting$assembly) else character()),
    events=events),file.path(out,'report_manifest.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
write_tsv(tibble(tool=c('R','ggplot2','tidyverse','jsonlite','patchwork'),version=c(as.character(getRversion()),
    map_chr(c('ggplot2','tidyverse','jsonlite','patchwork'),function(x) as.character(packageVersion(x))))),file.path(out,'reporting_versions.tsv'))
