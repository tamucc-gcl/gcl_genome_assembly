# Main report presentation. All inputs are explicit current-run artifacts.
suppressPackageStartupMessages({library(tidyverse);library(jsonlite)})
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
    haplotype=x$haplotype,completion=x$completion,graph_member=isTRUE(x$graph_member),fasta=x$fasta)) |>
    filter(completion=='finalized') |> arrange(taxid,sample,haplotype)
if(anyDuplicated(assemblies$id)) stop('Duplicate finalized identity in report')
names_data <- map_dfr(chars(spec$name_maps),function(path) optional_tsv(path) |>
    mutate(assembly=sub('\\.name_map\\.tsv$','',basename(path)),length=as.numeric(length)))
if(anyDuplicated(names_data[c('assembly','new_name')])) stop('Duplicate final scaffold identity')
if(!all(names_data$assembly %in% assemblies$id) || any(!is.finite(names_data$length) | names_data$length<=0)) stop('Invalid final name-map inventory')
write_tsv(names_data,file.path(out,'scaffold_inventory.tsv'))

md <- c('# Genome assembly results','',paste0(nrow(assemblies),' finalized assemblies from ',n_distinct(assemblies$sample),' individuals.'),'')
qc <- optional_tsv(spec$qc)
if(nrow(qc)>0) {
    hap_columns <- grep('^hap[0-9]+$',names(qc),value=TRUE)
    roles <- assemblies |> transmute(sample_id=sample,
        qc_haplotype=if_else(haplotype=='primary','hap1',haplotype),haplotype)
    kmer_title <- if(any(assemblies$haplotype=='primary')) 'Assembly k-mer completeness (%)' else 'Haplotype k-mer completeness (%)'
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
    md <- c(md,'## Assembly quality','',
        'Each cell lists the assembled copies in order; a primary/collapsed assembly has one value. K-mer completeness is per assembly, not combined diploid completeness.','',md_table(grid),'',
        link('All QC measurements','qc/assembly/assembly_qc_metrics.csv'))
    interactive <- artifacts |> filter(type=='interactive_qc')
    if(nrow(interactive)>0) md <- c(md,paste0(' · ',link('Interactive QC rankings',interactive$path[[1]])))
    md <- c(md,'')
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
        'Bar length shows scaffold size; internal colours show chromosome assignments from the shared final reference. Grey is unassigned and dark grey is conflicting assignment. The reference panel uses its own chromosome labels.','')
    for(taxid in unique(assemblies$taxid)) {
        members <- assemblies |> filter(.data$taxid==.env$taxid,id %in% painted_ids)
        if(nrow(members)==0) next
        top <- names_data |> filter(assembly %in% members$id) |> group_by(assembly) |>
            slice_max(length,n=as.integer(spec$run$top_scaffolds),with_ties=FALSE) |>
            arrange(length,.by_group=TRUE) |> ungroup() |>
            mutate(row=paste(assembly,new_name,sep='|'),row=factor(row,levels=unique(row)))
        blocks <- painting |> inner_join(top |> select(assembly,new_name,row),by=c('assembly','sequence'='new_name'))
        p <- ggplot(blocks,aes(xmin=start_0based/1e6,xmax=end_exclusive/1e6,
                              ymin=as.numeric(row)-.35,ymax=as.numeric(row)+.35,fill=chromosome)) +
            geom_rect() + facet_wrap(~factor(assembly,levels=members$id),ncol=max(1,min(4,n_distinct(members$haplotype))),scales='free_y') +
            scale_y_continuous(breaks=seq_along(levels(top$row)),labels=sub('^[^|]+\\|','',levels(top$row))) +
            scale_x_continuous(limits=c(0,max(top$length)/1e6),expand=expansion(mult=c(0,.02))) +
            scale_fill_manual(values=palette,name='Reference chromosome',drop=FALSE) +
            labs(title=paste0('Scaffold lengths and chromosome assignments · taxid ',taxid),
                 subtitle=paste0('Longest ',spec$run$top_scaffolds,' scaffolds per assembly; final coordinates'),
                 x='Position along scaffold (Mb)',y=NULL,
                 caption='Colours summarize primary mapping-chain spans, not base-level coverage or proof of a supported join. All panels use a shared length scale.') +
            theme_minimal(base_size=10) + theme(panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
                 plot.background=element_rect(fill='white',colour=NA),legend.position='bottom',axis.text.y=element_text(size=7),
                 strip.text=element_text(size=9)) + guides(fill=guide_legend(nrow=2))
        ncols <- max(1,min(4,n_distinct(members$haplotype)))
        filename <- paste0('chromosome_structure_',taxid)
        height <- max(4,ceiling(nrow(members)/ncols)*3)
        ggsave(file.path(out,paste0(filename,'.png')),p,width=14,height=height,dpi=160,bg='white')
        ggsave(file.path(out,paste0(filename,'.pdf')),p,width=14,height=height,bg='white',device=cairo_pdf)
        md <- c(md,paste0('![Scaffold lengths and chromosome assignments](reports/main/',filename,'.png)'),'',
                link('PDF',paste0('reports/main/',filename,'.pdf')),'')
        total <- names_data |> filter(assembly %in% members$id) |> group_by(assembly) |>
            summarise(total_scaffolds=n(),total_bp=sum(length),.groups='drop')
        scope <- bind_rows(scope,top |> group_by(assembly) |> summarise(shown_scaffolds=n(),shown_bp=sum(length),.groups='drop') |>
            left_join(total,by='assembly') |> mutate(taxid=.env$taxid))
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
        selected <- if(nrow(cuts)>0) cuts[1,] else if(nrow(retained)>0) retained[1,] else rows[1,]
        state <- if(nrow(cuts)>0) {
            if(all(cuts$verified %in% TRUE)) 'Cut applied and verified' else 'Cut selected; verification needed'
        } else if(nrow(retained)>0) 'Retained by review' else 'Review pending'
        tibble(chromosomes=selected$chromosomes,scaffold=selected$scaffold,status=state,
               cut_bp=if(nrow(cuts)>0) paste(cuts$cut_bp,collapse=', ') else '',reason=selected$reason)
    }) |> ungroup()
    write_tsv(events,file.path(out,'assembly_decisions.tsv'))
    pending <- events |> filter(!status %in% c('Cut applied and verified','Retained by review'))
    if(nrow(pending)>0) md <- c(md,'## Decisions needing attention','',md_table(pending |> select(Assembly=assembly,Event=event_id,Transition=chromosomes,Status=status)),'')
    md <- c(md,'## Assembly decisions and changes','',
        'These decisions describe joins in the original pre-finishing assemblies. A retained join is a recorded review choice, not proof of a biological rearrangement.','',
        md_table(events |> transmute(Assembly=assembly,Event=event_id,Transition=chromosomes,Outcome=status,
                                     `Original cut (bp)`=if_else(cut_bp=='','—',cut_bp))),'',
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
    md <- c(md,'## Assembly outcomes and files','',md_table(outcomes |> transmute(Individual=sample,Outcome=outcome,
        `Finalized / expected`=paste(observed,expected,sep=' / '),Note=reason)),'')
}
downloads <- assemblies |> mutate(files=paste0(link('FASTA',paste0('assembly/final/',basename(fasta))), ' · ',
    link('Name map',paste0('assembly/final/',id,'.name_map.tsv')))) |> select(Assembly=id,Files=files)
md <- c(md,'<details><summary>Final assemblies and chromosome name maps</summary>','',md_table(downloads),'','</details>','')
graph <- artifacts |> filter(type=='pangenome')
if(nrow(graph)>0) {
    md <- c(md,'## Pangenome','')
    for(taxid in unique(graph$id)) {
        rows <- graph |> filter(id==.env$taxid)
        ref <- keep(elig$linear_references,function(x) as.character(x$taxid)==taxid)
        ref_name <- if(length(ref)) ref[[1]]$id else 'Not selected'
        md <- c(md,paste0('**Taxid ',taxid,'** · Reference: ',ref_name,' · ',sum(assemblies$taxid==taxid & assemblies$graph_member),' contributing assemblies.'),'')
        notes <- rows$note[!is.na(rows$note) & rows$note!='']
        if(length(notes)>0) md <- c(md,notes,'')
        md <- c(md,paste0('- ',vapply(seq_len(nrow(rows)),function(i) link(report_label(rows$path[[i]]),rows$path[[i]]),character(1))),'')
    }
}
md <- c(md,'## Detailed evidence','',link('All available figures and reports','reports/main/evidence_index.md'),'')
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
        link('Maps, annotations and detailed evidence','reports/main/evidence_index.md#organelle-stats'),'')
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
for(type in unique(artifacts$type)) {
    rows <- artifacts |> filter(.data$type==.env$type)
    index <- c(index,paste0('## ',str_replace_all(type,'_',' ')),'',
        md_table(rows |> transmute(Item=id,Report=link(report_label(path),paste0('../../',path)))),'')
}
writeLines(index,file.path(out,'evidence_index.md'))
md <- c(md,'## Methods and reproducibility','',
    paste0('Run: **',spec$run$run_name,'** · Nextflow ',spec$run$nextflow,' · Profile ',spec$run$profile,'.'),'',
    if(!is.null(spec$run$commit) && nzchar(spec$run$commit)) paste0('Pipeline commit: `',spec$run$commit,'`.') else
        'Pipeline executed from a working tree; source fingerprints are recorded in run provenance.','',
    link('Software versions','provenance/software_versions.tsv'),' · ',link('Version metadata issues','provenance/software_version_issues.tsv'),' · ',
    link('Run parameters and source fingerprints','provenance/run_provenance.json'),' · ',link('Environment snapshots','provenance/runtime_environments.json'),' · ',
    link('Report renderer versions','reports/main/reporting_versions.tsv'),'')
writeLines(md,'assembly_report.md')
write_json(list(format_version=1,run=spec$run,assemblies=assemblies,
    qc_present=nrow(qc)>0,painting_required=isTRUE(spec$run$painting_enabled),
    expected_painting=as.list(expected_painting),painted_assemblies=as.list(if(nrow(painting)>0) unique(painting$assembly) else character()),
    events=events),file.path(out,'report_manifest.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
write_tsv(tibble(tool=c('R','ggplot2','tidyverse','jsonlite'),version=c(as.character(getRversion()),
    map_chr(c('ggplot2','tidyverse','jsonlite'),function(x) as.character(packageVersion(x))))),file.path(out,'reporting_versions.tsv'))
