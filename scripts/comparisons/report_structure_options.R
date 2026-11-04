# Read-only layout comparison using retained reporting tables; no new mapping.
suppressPackageStartupMessages({library(tidyverse);library(patchwork)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1) stop('Supply the results directory to assess')
root <- normalizePath(args[[1]],mustWork=TRUE)
out <- file.path(root,'reports','layout-options')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
inventory <- read_tsv(file.path(root,'reports/main/scaffold_inventory.tsv'),show_col_types=FALSE)
blocks <- read_tsv(file.path(root,'reports/main/painting_blocks.tsv'),show_col_types=FALSE)
manifest <- jsonlite::fromJSON(file.path(root,'reports/main/report_manifest.json'))
members <- as_tibble(manifest$assemblies) |> arrange(taxid,sample,haplotype)
top <- inventory |> filter(str_detect(new_name,'^chr')) |> group_by(assembly) |>
    slice_max(length,n=20,with_ties=FALSE) |> ungroup() |>
    mutate(first_chrom=str_remove(str_split_i(new_name,fixed('+'),1),'_[0-9]+$'),
           chrom_number=as.numeric(str_extract(first_chrom,'[0-9]+$')),
           piece=as.numeric(str_extract(str_split_i(new_name,fixed('+'),1),'(?<=_)[0-9]+$')),
           aggregate=FALSE)
unplaced <- inventory |> filter(!str_detect(new_name,'^chr')) |> group_by(assembly) |>
    summarise(length=sum(length),pieces=n(),.groups='drop') |>
    mutate(new_name=paste0('Unplaced total (',format(pieces,big.mark=','),' pieces)'),
           first_chrom='Unplaced',chrom_number=Inf,piece=Inf,aggregate=TRUE)
write_tsv(unplaced,file.path(out,'unplaced_totals.tsv'))
chromosomes <- unique(blocks$chromosome)
chromosomes <- chromosomes[!chromosomes %in% c('Unassigned','Ambiguous')]
chromosomes <- chromosomes[order(as.numeric(str_extract(chromosomes,'[0-9]+$')),chromosomes,na.last=TRUE)]
palette <- c(setNames(hcl.colors(length(chromosomes),'Dark 3'),chromosomes),
             Ambiguous='#666666',Unassigned='#D9D9D9',`Unplaced sum`='#44546A')
theme_set(theme_minimal(base_size=10)+theme(panel.grid.minor=element_blank(),
    panel.grid.major.y=element_blank(),plot.background=element_rect(fill='white',colour=NA),
    legend.position='bottom',axis.text.y=element_text(size=7)))

structure_plot <- function(include_unplaced) {
    rows <- if(include_unplaced) bind_rows(top,unplaced) else top
    rows <- rows |> arrange(assembly,desc(aggregate),desc(chrom_number),desc(first_chrom),desc(piece),desc(new_name)) |>
        mutate(row=paste(assembly,new_name,sep='|'),row=factor(row,levels=unique(row)))
    write_tsv(rows,file.path(out,if(include_unplaced) 'same_scale_rows.tsv' else 'separate_scale_rows.tsv'))
    draw <- blocks |> inner_join(rows |> filter(!aggregate) |> select(assembly,new_name,row),by=c('assembly','sequence'='new_name'))
    if(include_unplaced) draw <- bind_rows(draw,rows |> filter(aggregate) |>
        transmute(assembly,sequence=new_name,length,start_0based=0,end_exclusive=length,chromosome='Unplaced sum',row))
    maximum <- max(rows$length)/1e6
    ggplot(draw,aes(x=(start_0based+end_exclusive)/2e6,y=row,width=(end_exclusive-start_0based)/1e6,fill=chromosome)) +
        geom_tile(height=.7) +
        facet_wrap(~factor(assembly,levels=members$id),ncol=2,scales='free_y') +
        scale_y_discrete(labels=function(x) sub('^[^|]+\\|','',x)) +
        scale_x_continuous(limits=c(0,maximum),expand=expansion(mult=c(0,.02))) +
        scale_fill_manual(values=palette,breaks=c(chromosomes,'Unassigned','Ambiguous',if(include_unplaced) 'Unplaced sum'),drop=FALSE,name='Reference chromosome') +
        guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
        labs(title=if(include_unplaced) 'Option 1: chromosome order and unplaced sum on one scale' else 'Option 2: chromosome order; unplaced sum shown separately',
             subtitle='Real final assemblies; individual scaffolds sorted by chromosome, then piece',
             x='Scaffold length / position (Mb)',y=NULL,
             caption=if(include_unplaced) 'The unplaced bar is the sum of independent pieces, not a joined scaffold. All panels share one scale.' else
                 'Longest 20 named scaffolds per assembly. Smaller named pieces remain in the linked inventory. Reference identity panel: Sde-CMat_203_hap1.')
}
save <- function(name,p,height) {
    ggsave(file.path(out,paste0(name,'.png')),p,width=14,height=height,dpi=160,bg='white')
    ggsave(file.path(out,paste0(name,'.pdf')),p,width=14,height=height,device=cairo_pdf,bg='white')
}
save('sorted_same_scale',structure_plot(TRUE),17)
save('sorted_table_option',structure_plot(FALSE),17)
totals <- unplaced |> mutate(assembly=factor(assembly,levels=rev(members$id)),
    text=paste0(sprintf('%.1f Mb',length/1e6),' · ',format(pieces,big.mark=','),' pieces'))
sum_plot <- ggplot(totals,aes(length/1e6,assembly)) + geom_col(fill='#44546A',width=.55) +
    geom_text(aes(label=text),hjust=-.05,size=2.7) +
    scale_x_continuous(expand=expansion(mult=c(0,.4))) +
    labs(title='All unplaced sequence: summed length, with a separate scale',
         subtitle='Each bar pools independent pieces; colours do not represent chromosome assignment here',
         x='Summed unplaced length (Mb)',y=NULL) + theme(legend.position='none')
save('unplaced_summary',sum_plot,4.3)
save('sorted_separate_scale',structure_plot(FALSE)/sum_plot + plot_layout(heights=c(4,1.3)),20)
writeLines(c('# Structure plot options','',
    'Both options use the same retained final sequence lengths and painting blocks. Sorting and presentation change; evidence does not.','',
    '## Option 1: one shared length scale','',
    '![Sorted with unplaced sum](sorted_same_scale.png)','',
    '## Option 2: separate summary scale','',
    '![Sorted with separate unplaced summary](sorted_separate_scale.png)','',
    'Unplaced is a final naming category. Unassigned is an interval with no chromosome label passing the painting filters. These are different: an unplaced piece can still align to a reference chromosome.','',
    '[Unplaced lengths and piece counts](unplaced_totals.tsv)','',
    '## Option 3: chromosome plot and native Markdown table','',
    '![Sorted chromosome plot](sorted_table_option.png)','',
    '| Assembly | Unplaced pieces | Total unplaced (Mb) |','| --- | ---: | ---: |',
    paste0('| ',unplaced$assembly,' | ',format(unplaced$pieces,big.mark=','),' | ',sprintf('%.1f',unplaced$length/1e6),' |')),
    file.path(out,'README.md'))
for(option in 1:3) {
    body <- c(paste0('# Option ',option),'',
        c('One shared scale: individual scaffolds and the unplaced sum.',
          'Chromosome panels retain their scale; the unplaced summary has its own scale.',
          'Chromosome panels plus a native Markdown table of unplaced totals.')[[option]],'',
        paste0('![Structure option ',option,'](',c('sorted_same_scale.png','sorted_separate_scale.png','sorted_table_option.png')[[option]],')'),'')
    if(option==3) body <- c(body,'| Assembly | Unplaced pieces | Total unplaced (Mb) |','| --- | ---: | ---: |',
        paste0('| ',unplaced$assembly,' | ',format(unplaced$pieces,big.mark=','),' | ',sprintf('%.1f',unplaced$length/1e6),' |'))
    writeLines(body,file.path(out,paste0('option-',option,'.md')))
}
