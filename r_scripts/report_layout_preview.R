# Presentation prototypes from retained outputs; does not alter pipeline evidence.
# Run with R 4.3, tidyverse 2.0.0 (the tested summary-report environment).
suppressPackageStartupMessages(library(tidyverse))
argv <- commandArgs(trailingOnly=TRUE)
if (length(argv)!=1) stop('Supply the accepted results directory')
root <- normalizePath(argv[[1]], mustWork=TRUE)
out <- file.path(root,'reports','layout-preview')
dir.create(out, recursive=TRUE, showWarnings=FALSE)
theme_set(theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank(),
    plot.background=element_rect(fill='white',colour=NA),legend.position='bottom'))
save_plot <- function(name,p,width,height) {
    ggsave(file.path(out,paste0(name,'.png')),p,width=width,height=height,dpi=160,bg='white')
    ggsave(file.path(out,paste0(name,'.pdf')),p,width=width,height=height,bg='white',device=cairo_pdf)
}
qc <- read_csv(file.path(root,'qc/assembly/assembly_qc_metrics.csv'),show_col_types=FALSE) |>
    filter(stage=='final') |>
    pivot_longer(c(hap1,hap2),names_to='haplotype',values_to='value') |>
    mutate(value=as.numeric(value),assembly=paste(sample_id,haplotype,sep='_'))
totals <- qc |> filter(analysis=='busco',metric=='total_busco') |>
    select(assembly,total=value)
if(anyDuplicated(totals$assembly)) stop('Ambiguous BUSCO denominator')
quality <- qc |>
    filter((analysis=='busco' & metric=='complete') |
           (analysis=='merqury' & metric %in% c('qv','kmer_completeness')) |
           (analysis=='quast' & metric=='N50')) |>
    left_join(totals,by='assembly') |>
    mutate(display_value=case_when(metric=='complete' & total>0 ~ 100*value/total,
                                  metric=='N50' ~ value/1e6,TRUE ~ value),
        panel=case_when(metric=='complete' ~ 'BUSCO complete (%)',metric=='qv' ~ 'QV',
                       metric=='kmer_completeness' ~ 'Haplotype k-mer completeness (%)',
                       metric=='N50' ~ 'Scaffold N50 (Mb)'),
        panel=factor(panel,levels=c('BUSCO complete (%)','QV','Haplotype k-mer completeness (%)','Scaffold N50 (Mb)')),
        assembly=factor(assembly,levels=rev(sort(unique(assembly)))),
        label=if_else(is.na(display_value),'Not measured',sprintf('%.1f',display_value)))
if(anyDuplicated(quality[c('assembly','panel')])) stop('Duplicate final QC measurement')
write_tsv(quality,file.path(out,'quality_values.tsv'))

maps <- list.files(file.path(root,'assembly/final'),pattern='\\.name_map\\.tsv$',full.names=TRUE)
inventory <- map_dfr(maps,function(path) read_tsv(path,show_col_types=FALSE) |>
    mutate(assembly=str_remove(basename(path),'\\.name_map\\.tsv$'))) |>
    mutate(label_type=case_when(str_detect(new_name,fixed('+')) ~ 'Composite chromosome label',
                              str_detect(new_name,'^chr') ~ 'Single chromosome label',
                              TRUE ~ 'Other / unplaced'))
write_tsv(inventory,file.path(out,'scaffold_inventory.tsv'))
top <- inventory |> group_by(assembly) |> slice_max(length,n=20,with_ties=FALSE) |>
    arrange(assembly,length) |> ungroup() |>
    mutate(row=paste(assembly,new_name,sep='|'),row=factor(row,levels=unique(row)))
p <- ggplot(top,aes(length/1e6,row,fill=label_type)) + geom_col(width=.75) +
    facet_wrap(~assembly,ncol=2,scales='free_y') +
    scale_y_discrete(labels=function(x) sub('^[^|]+\\|','',x)) +
    scale_fill_manual(values=c('Single chromosome label'='#0072B2','Composite chromosome label'='#D55E00','Other / unplaced'='#969696'),name=NULL) +
    labs(title='A · Scaffold lengths and chromosome names',
         subtitle='Real final assemblies · longest 20 scaffolds per assembly · shared length scale',
         x='Scaffold length (Mb)',y=NULL,
         caption='Whole-scaffold names describe assigned groups. A composite label does not localize its internal chromosome boundary. Smaller scaffolds remain in the linked inventory.') +
    theme(axis.text.y=element_text(size=7),strip.text=element_text(size=10))
save_plot('structure_A_inventory',p,13,17)

# Deliberately fictional coordinates: compare a proposed layout, never infer biology.
fiction <- crossing(specimen=LETTERS[1:5],haplotype=1:2,scaffold=1:4) |>
    mutate(assembly=paste('Example',specimen,paste0('hap',haplotype)),
           end=c(80,65,55,42)[scaffold],start=0,chromosome=paste0('chr',scaffold))
special <- fiction |> filter(specimen=='E',scaffold==2) |>
    mutate(start=32,end=65,chromosome='chr4')
fiction <- fiction |> mutate(end=if_else(specimen=='E' & scaffold==2,32,end)) |>
    bind_rows(special)
p <- ggplot(fiction,aes(xmin=start,xmax=end,ymin=scaffold-.24,ymax=scaffold+.24,fill=chromosome)) +
    geom_rect(colour='white',linewidth=.2) + facet_wrap(~assembly,ncol=2) +
    scale_fill_manual(values=c(chr1='#0072B2',chr2='#D55E00',chr3='#009E73',chr4='#CC79A7'),name='Aligned chromosome') +
    scale_y_reverse(breaks=1:4,labels=paste('Scaffold',1:4)) +
    labs(title='B · Chromosome painting — LAYOUT EXAMPLE ONLY',
         subtitle='Fictional specimens and blocks · these are not measurements from this run',
         x='Position along each scaffold (illustrative)',y=NULL,
         caption='A production painting would use retained, coordinate-matched alignment blocks; unassigned sequence would remain grey. Final name maps alone cannot supply these blocks.') +
    theme(axis.text.x=element_blank(),panel.grid=element_blank(),strip.text=element_text(size=10))
save_plot('structure_B_painting_concept',p,13,10)

# Combined A+B design: scaffold length/name and internal painting on the same bar.
# Entire figure is fictional so illustrative block boundaries cannot be mistaken
# for measurements on the accepted assemblies.
combined_scaffolds <- crossing(specimen=LETTERS[1:5],haplotype=1:2,rank=1:6) |>
    mutate(assembly=paste('Example',specimen,paste0('hap',haplotype)),
           scaffold=c('chr1_1','chr2_1','chr3_1','chr4_1','chr4_2','unplaced_1')[rank],
           chromosome=c('chr1','chr2','chr3','chr4','chr4','Unassigned')[rank],
           size=c(90,70,60,35,12,8)[rank] * (1+(haplotype-1)*.03),
           size=if_else(specimen=='E' & rank==2,size+35,size),
           scaffold=if_else(specimen=='E' & rank==2,'chr2_1+chr4_1',scaffold),
           row=paste(assembly,scaffold,sep='|')) |>
    arrange(assembly,desc(rank)) |>
    mutate(row=factor(row,levels=unique(row)))
combined_blocks <- pmap_dfr(combined_scaffolds,function(specimen,haplotype,rank,assembly,scaffold,chromosome,size,row) {
    if(specimen=='E' && rank==2) {
        tibble(assembly=assembly,row=row,start=c(0,size-35),end=c(size-35,size),chromosome=c('chr2','chr4'))
    } else if(specimen=='D' && rank==3) {
        tibble(assembly=assembly,row=row,start=c(0,20,25),end=c(20,25,size),chromosome=c('chr3','Unassigned','chr3'))
    } else tibble(assembly=assembly,row=row,start=0,end=size,chromosome=chromosome)
}) |> mutate(row=factor(row,levels=levels(combined_scaffolds$row)),y=as.numeric(row))
combined_scaffolds <- combined_scaffolds |> mutate(y=as.numeric(row))
p <- ggplot() +
    geom_rect(data=combined_blocks,aes(xmin=start,xmax=end,ymin=y-.31,ymax=y+.31,fill=chromosome)) +
    geom_rect(data=combined_scaffolds,aes(xmin=0,xmax=size,ymin=y-.31,ymax=y+.31),fill=NA,colour='#555555',linewidth=.2) +
    geom_text(data=combined_scaffolds,aes(size+2,y,label=sprintf('%.0f',size)),hjust=0,size=2.8) +
    facet_wrap(~assembly,ncol=2,scales='free_y') +
    scale_y_continuous(breaks=combined_scaffolds$y,labels=combined_scaffolds$scaffold) +
    scale_x_continuous(limits=c(0,120),breaks=c(0,30,60,90,120),expand=expansion(mult=c(0,.04))) +
    scale_fill_manual(values=c(chr1='#0072B2',chr2='#D55E00',chr3='#009E73',chr4='#CC79A7',Unassigned='#D9D9D9'),name='Chromosome assignment') +
    labs(title='Scaffold lengths and chromosome painting - combined layout',
         subtitle='FICTIONAL EXAMPLE: names, lengths and internal blocks illustrate the design',
         x='Position along scaffold / scaffold length (Mb, illustrative)',y=NULL,
         caption='One bar per scaffold; colour changes show chromosome transitions. Grey marks unassigned intervals. Haplotype 1 and 2 are side by side; all panels use the same length scale.') +
    theme(panel.grid.major.y=element_blank(),axis.text.y=element_text(size=8),strip.text=element_text(size=10),plot.caption=element_text(size=9))
save_plot('structure_AB_combined_concept',p,14,12)

md_table <- function(x) {
    c(paste0('| ',paste(names(x),collapse=' | '),' |'),
      paste0('| ',paste(rep('---',ncol(x)),collapse=' | '),' |'),
      apply(x,1,function(row) paste0('| ',paste(replace(row,is.na(row),'Not measured'),collapse=' | '),' |')))
}
review_path <- file.path(root,'assembly/chimeras/review.tsv')
events <- tibble()
if(file.exists(review_path)) {
    review <- read_tsv(review_path,col_types=cols(.default=col_character()))
    verification <- list.files(file.path(root,'assembly/chimeras'),pattern='verification\\.json$',full.names=TRUE)
    verified <- map_dfr(verification,function(path) {
        v <- jsonlite::fromJSON(path)
        if(!isTRUE(v$exact_parent_reconstruction) || is.null(v$selected_actions)) return(tibble())
        as_tibble(v$selected_actions) |> select(assembly,id,cut_bp,assessment_sha256)
    }) |> distinct() |> mutate(verified=TRUE)
    review <- review |> left_join(verified,by=c('assembly','id','cut_bp','assessment_sha256'))
    events <- review |> group_by(assembly,event_id) |> group_modify(function(rows,key) {
        cuts <- rows |> filter(selected=='YES',review_disposition=='CUT')
        keep <- rows |> filter(review_disposition=='RETAIN')
        chosen <- if(nrow(cuts)>0) cuts[1,] else if(nrow(keep)>0) keep[1,] else rows[1,]
        state <- if(nrow(cuts)>0) {
            if(all(cuts$verified %in% TRUE)) 'Cut applied and verified' else 'Cut selected; verification needed'
        } else if(nrow(keep)>0) 'Retained by review' else 'Review pending'
        tibble(chromosomes=chosen$chromosomes,scaffold=chosen$scaffold,status=state,
               cut_bp=if(nrow(cuts)>0) paste(cuts$cut_bp,collapse=', ') else '',reason=chosen$reason)
    }) |> ungroup()
    write_tsv(events,file.path(out,'decision_values.tsv'))
}
pair_plan <- read_tsv(file.path(root,'pairwise_alignments/comparison_plan.tsv'),show_col_types=FALSE)
candidate_pair <- if(nrow(events)>0) pair_plan |> filter(ref_id %in% events$assembly | query_id %in% events$assembly) else pair_plan
if(nrow(candidate_pair)==0) candidate_pair <- pair_plan
pair <- candidate_pair[1,]
ribbon <- paste0('pairwise_alignments/',pair$ref_id,'_vs_',pair$query_id,'_riparian.png')
invisible(file.copy(file.path(root,ribbon),file.path(out,'structure_C_ribbon.png'),overwrite=TRUE))
numbers <- quality |> select(sample_id,panel,haplotype,display_value) |>
    mutate(display_value=if_else(is.na(display_value),'Not measured',sprintf('%.1f',display_value))) |>
    pivot_wider(names_from=haplotype,values_from=display_value) |>
    mutate(value=paste(hap1,hap2,sep=' / ')) |>
    select(sample_id,panel,value) |>
    pivot_wider(names_from=panel,values_from=value) |> arrange(sample_id)
names(numbers)[1] <- 'Individual'
quality_note <- 'Each cell shows **haplotype 1 / haplotype 2**. K-mer completeness is per haplotype; combined diploid completeness is a different measurement.'
writeLines(c('# Final assembly quality','',quality_note,'',md_table(numbers)),
           file.path(root,'quality-overview.md'))
writeLines(c('# Quality overview','',
    'The quality figures have been dropped. [Open the native Markdown table](quality-overview.md).'),
    file.path(root,'quality-layout-comparison.md'))

decisions <- events |> transmute(Assembly=assembly,Transition=chromosomes,
    Outcome=status,`Original scaffold`=scaffold,`Cut (bp)`=if_else(cut_bp=='','—',cut_bp))
decision_detail <- c('# Chimera decisions — layout preview','',
    'Recorded decisions and verified actions from this run. Coordinates refer to the original pre-finishing scaffolds.','',
    md_table(decisions),'',
    'Retained by review records a reviewer choice; it does not establish that a rearrangement is biologically proven.','')
for(i in seq_len(nrow(events))) {
    e <- events[i,]
    decision_detail <- c(decision_detail,paste0('## ',e$assembly,' · ',e$event_id),'',
        paste0('**',e$status,'** · ',e$chromosomes),'',e$reason,'',
        paste0('[Original evidence](../../assembly/chimeras/assemblies/',e$assembly,'/report.md#',tolower(e$event_id),')'),'')
}
writeLines(decision_detail,file.path(out,'chimera-decisions.md'))
graph_audits <- list.files(file.path(root,'pangenome'),pattern='source_audit.json$',recursive=TRUE,full.names=TRUE)
graph_lines <- map(graph_audits,function(path) {
    a <- jsonlite::fromJSON(path)
    taxid <- basename(dirname(dirname(path)))
    base <- paste0('pangenome/',taxid)
    c(paste0('**Taxid ',taxid,'**: ',a$input_sequences,' source sequences; ',
        a$chromosome_scale_paths_requiring_review,' chromosome-scale paths flagged for representation review.'),'',
      paste0('[Graph construction and inputs](',base,'/pangenome_build_report.md) · [Representation review](',base,
             '/representation/source_report.md) · [Sharing results](',base,'/sharing/sharing_report.md)'), '')
}) |> unlist()
preview <- c('# Genome assembly results','',
    '> Layout preview using the accepted integration run. Quality and recorded decisions are real. The combined chromosome-structure figure is fictional and illustrates the proposed design.','',
    '## Assembly quality','',
    paste0(n_distinct(quality$assembly),' finalized assemblies from ',n_distinct(qc$sample_id),
        ' individuals. Compare the measured metrics directly; workflow completion is separate from biological quality.'),'',
    quality_note,'',md_table(numbers),'',
    '[Interactive QC ranking](reports/assembly_qc_report.html) · [All QC measurements](qc/assembly/assembly_qc_metrics.csv)','',
    '## Chromosome structure','',
    '![Scaffold lengths with chromosome painting](reports/layout-preview/structure_AB_combined_concept.png)','',
    '**Illustrative data only.** Each bar combines the scaffold name and length with chromosome assignments along it. Haplotype 1 and 2 appear side by side. The final figure will use measured, coordinate-matched blocks; no boundaries have been inferred from scaffold names.','',
    '[Real scaffold-length inventory from this run](reports/layout-preview/structure_A_inventory.png)','',
    '<details><summary>Haplotype correspondence: ribbon example</summary>','',
    '![Haplotype ribbon plot](reports/layout-preview/structure_C_ribbon.png)','',
    paste0('Real retained comparison: ',pair$ref_id,' and ',pair$query_id,'. All pairs remain accessible below.'),'','</details>','',
    '## Recorded assembly decisions','',md_table(decisions),'',
    paste0(sum(events$status=='Cut applied and verified'),' cuts were applied and sequence-verified. ',
        sum(events$status=='Retained by review'),' transitions were retained by review. These rows describe the original assessed joins and what happened to them.'),'',
    '[Decision reasons and original evidence](reports/layout-preview/chimera-decisions.md)','',
    '## Pangenome','',
    graph_lines,
    '## Detailed results','',
    '[Assembly and input outcomes](assembly_run_summary.md) · [QC ranking](reports/assembly_qc_report.html) · [Software and environment provenance](provenance/run_provenance.json)','',
    '<details><summary>All haplotype comparisons</summary>','')
for(i in seq_len(nrow(pair_plan))) {
    pair <- pair_plan[i,];prefix <- paste0('pairwise_alignments/',pair$ref_id,'_vs_',pair$query_id)
    preview <- c(preview,paste0('- **',pair$ref_sample,'**: [dotplot](',prefix,'_dotplot.png) · [ribbon](',prefix,'_riparian.png)'))
}
preview <- c(preview,'','</details>','',
    '<details><summary>Organelle, telomere and contact-map evidence</summary>','',
    'These existing figure sets will be linked per sample in the detailed report; this preview focuses on the opening layout.','','</details>')
writeLines(preview,file.path(root,'report-layout-preview.md'))
writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
message('Preview figures written to ',out)
