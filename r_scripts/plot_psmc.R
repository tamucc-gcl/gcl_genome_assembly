# Presentation consumes validated, scaled tables; it never calls variants or refits PSMC.
suppressPackageStartupMessages({library(ggplot2); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=4L) stop('Expected context.json input_directory curves_directory output_directory')
context <- fromJSON(args[[1]])
qc <- fromJSON(file.path(args[[2]],'input_qc.json'))
out <- args[[4]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
ready <- identical(qc$status,'ready')
scaled <- !is.null(context$mutation_rate)
summary <- list(sample=context$sample,taxid=context$taxid,mode=context$mode,
    reference=context$reference_id,status=if(ready) 'completed' else 'not_fitted',
    reason=qc$reason,scope_bp=qc$scope_bp,callable_bp=qc$callable_bp,heterozygous_snps=qc$heterozygous_snps,
    retained_callable_bp=qc$psmc_retained_callable_bp,retained_scaffolds=qc$segments,
    retained_contigs=paste(qc$retained_contigs,collapse=', '),missing_bin_fraction=qc$missing_bin_fraction,
    scaling=if(scaled) 'years / Ne' else 'unscaled',
    report=paste0('post_assembly/psmc/',context$taxid,'/',context$sample,'/report/report.md'))
md <- c(paste0('# PSMC: ',context$sample),'',
    paste0('Species: ',context$species,' (taxid ',context$taxid,'). Input mode: **',context$mode,'**.'),'',
    paste0('Reference: **',context$reference_id,'**. SHA256: `',context$reference_sha256,'`.'),'',
    paste0('Reference selection: `',context$reference_selection,'`. [Own reference provenance](../reference/reference.json).'),'',
    '| Measurement | Value |','| --- | --- |',
    paste0('| Outcome | ',summary$status,' |'),
    paste0('| Assembly representation | ',context$representation,' |'),
    paste0('| Callable sequence before conversion | ',format(qc$callable_bp,scientific=FALSE,trim=TRUE),' bp |'),
    paste0('| Callable fraction of selected non-gap scope | ',sprintf('%.2f%%',100*qc$callable_fraction),' |'),
    paste0('| Heterozygous SNPs | ',format(qc$heterozygous_snps,scientific=FALSE,trim=TRUE),' |'),
    paste0('| Heterozygosity per callable bp | ',if(is.null(qc$heterozygosity)) 'Not measurable' else signif(qc$heterozygosity,5),' |'),
    paste0('| PSMC segments | ',qc$segments,' |'),
    paste0('| Callable bp on retained PSMC sequences | ',format(qc$psmc_retained_callable_bp,scientific=FALSE,trim=TRUE),' |'),
    paste0('| Unknown PSMC bins | ',sprintf('%.2f%%',100*qc$missing_bin_fraction),' |'),
    paste0('| Bin size | ',qc$bin_size,' bp |'),
    paste0('| Scaling | ',summary$scaling,' |'),'')
if(context$mode=='assemblies') md <- c(md,
    '| Assembly-calling setting | Value |','| --- | --- |',
    paste0('| Minimap2 preset | ',context$assembly_preset,' |'),
    paste0('| Minimum alignment for masks and variant filters | ',context$assembly_min_alignment_bp,' bp |'),'',
    '[Exact caller settings](../genotypes/dipcall_settings.json) · [Generated-command changes](../genotypes/dipcall.makefile.diff)','')
if(qc$psmc_retained_callable_bp < qc$callable_bp) md <- c(md,
    paste0('**Input scope warning:** conversion retained ',qc$segments,' scaffolds and ',
      format(qc$psmc_retained_callable_bp,scientific=FALSE,trim=TRUE),' callable bp of ',
      format(qc$callable_bp,scientific=FALSE,trim=TRUE),' available before conversion. Completion does not establish representative genomic scope.'),'')
if(length(qc$retained_contigs)) md <- c(md,paste0('Retained scaffolds: ',paste(qc$retained_contigs,collapse=', ')), '')
if(!is.null(qc$mask_audit_totals)) {
    a <- qc$mask_audit_totals
    md <- c(md,'## Where sequence was lost','',
       '| Stage | Callable bp remaining |','| --- | ---: |',
       paste0('| Selected non-gap scope | ',qc$scope_bp,' |'),
       if(!is.null(a$confident_scope_bp)) paste0('| Inside Dipcall confident regions | ',a$confident_scope_bp,' |'),
       paste0('| After simple genotype checks | ',a$after_simple_genotypes_bp,' |'),
       paste0('| After complex-variant masking | ',a$after_simple_genotypes_bp-a$complex_mask_loss_bp,' |'),
       paste0('| After duplicate-position masking | ',qc$callable_bp,' |'),
       paste0('| On sequences retained by fq2psmcfa | ',qc$psmc_retained_callable_bp,' |'),'',
       '[Per-scaffold mask losses](../input/mask_audit.tsv) · [Conversion retention and rejection reasons](../input/conversion_audit.tsv)','',
       'fq2psmcfa requires at least 20% callable bases across each whole input scaffold, in addition to its minimum good-base count. Unknown intervals retain their original lengths.','')
}
diagnostics <- read.delim(file.path(out,'psmc_qc.tsv'),colClasses='character',na.strings='')
diagnostic_value <- function(value, unit, status) {
    if(is.na(value)) return(if(status=='open-ended') 'Open-ended' else 'Unavailable')
    number <- suppressWarnings(as.numeric(value))
    if(identical(unit,'fraction')) return(sprintf('%.3g%%',100*number))
    if(!is.na(number)) value <- format(signif(number,6),trim=TRUE)
    if(is.na(unit)) value else paste(value,unit)
}
values <- mapply(diagnostic_value,diagnostics$value,diagnostics$unit,diagnostics$status)
md <- c(md,'## PSMC diagnostics','',
    '| Category | Metric | Value | Screening status |','| --- | --- | ---: | --- |',
    paste0('| ',diagnostics$category,' | ',diagnostics$metric,' | ',values,' | ',diagnostics$status,' |'),'',
    '[Diagnostic metrics (TSV)](psmc_qc.tsv)','',
    'Input flags use 50% selected-scope callability, 50% masked retained bins and a heterozygosity review range of 0.01–2%. These are screening thresholds. Usable bins can contain masked bases.','')
if(ready) {
    intervals <- read.delim(file.path(out,'psmc_intervals.tsv'))
    iterations <- read.delim(file.path(out,'psmc_iterations.tsv'))
    support_colors <- c(insufficient='#B73D31',caution='#A46200',reasonable='#0072B2')
    support_scale <- scale_color_manual(values=support_colors,limits=names(support_colors),
        labels=c('Below 10','10 to below 20','20 or more'),name='Expected recombinations')
    positive_times <- intervals$time_start[intervals$time_start>0]
    time_scale <- if(length(positive_times)) min(positive_times) else 1
    intervals$plot_start <- intervals$time_start/time_scale
    intervals$plot_end <- intervals$time_end/time_scale
    time_ticks <- expm1(pretty(log1p(range(intervals$plot_start)),n=5))
    time_axis <- scale_x_continuous(trans='log1p',breaks=time_ticks,minor_breaks=NULL,
        labels=function(x) format(x*time_scale,digits=2,trim=TRUE))
    time_label <- if(scaled) 'Time before present (years)' else 'Scaled divergence'
    peaks <- intervals[intervals$recent_peak=='yes',]
    subtitle <- paste(context$sample,'· observed fit, iteration',unique(intervals$iteration))
    history_qc <- ggplot(intervals,aes(x=plot_start,y=population,color=support)) +
        geom_segment(data=intervals[!is.na(intervals$plot_end),],aes(xend=plot_end,yend=population),linewidth=1.2) +
        geom_point(size=1.7) + geom_point(data=peaks,shape=21,fill='white',color='black',size=3,stroke=1) +
        support_scale + time_axis + expand_limits(y=0) + theme_bw(base_size=11) + theme(legend.position='bottom') +
        labs(title='Observed history by interval support',subtitle=subtitle,x=time_label,
             y=if(scaled) 'Effective population size (Ne)' else 'Scaled size',
             caption='Open circles flag recent peaks for review. The terminal interval is shown at its start.')
    recombination_qc <- ggplot(intervals,aes(plot_start,expected_recombinations,color=support)) +
        geom_hline(yintercept=c(10,20),linetype='dashed',color='grey45') +
        geom_point(size=2.3) + support_scale + time_axis +
        scale_y_continuous(trans='log1p',breaks=function(limits)
            c(0,1,10,20,10^seq.int(2,max(2,floor(log10(max(limits,20,na.rm=TRUE))))))) +
        theme_bw(base_size=11) + theme(legend.position='bottom') +
        labs(title='Recombination support across time',subtitle=subtitle,x=paste(time_label,'at interval start'),
             y='Expected recombinations',caption='Each point represents one free interval. Count axis uses log(1 + count).')
    observed_likelihood <- iterations[is.finite(iterations$log_likelihood),]
    likelihood_qc <- ggplot(observed_likelihood,aes(iteration,log_likelihood)) +
        geom_point(color='#0072B2',size=2) +
        theme_bw(base_size=11) + labs(title='Observed log-likelihood by iteration',subtitle=context$sample,
            x='Iteration',y='Log-likelihood',caption=if(nrow(observed_likelihood))
                'PSMC reports likelihood before the parameter update recorded in each iteration.' else
                'No fitted-iteration likelihood values were available.')
    if(nrow(observed_likelihood)>1L) likelihood_qc <- likelihood_qc + geom_line(color='#0072B2')
    figures <- list(psmc_history_support=history_qc,psmc_recombination_support=recombination_qc,
                    psmc_likelihood=likelihood_qc)
    for(name in names(figures))
        ggsave(file.path(out,paste0(name,'.png')),figures[[name]],width=8,height=4.8,dpi=160,bg='white')
    md <- c(md,'[Free intervals](psmc_intervals.tsv) · [Iteration history](psmc_iterations.tsv)','',
        '![Observed history by interval support](psmc_history_support.png)','',
        '![Recombination support across time](psmc_recombination_support.png)','',
        '![Observed log-likelihood by iteration](psmc_likelihood.png)','',
        'Temporal endpoints describe the first and last groups with at least 20 expected recombinations. Groups between them can have lower support. Size changes compare corresponding groups across iterations. Review convergence from the changes and likelihood trajectory.','',
        'The peak screen checks for a group in the first five that is at least twice as high as each immediate neighbor. A flagged peak warrants parameter sensitivity checks. Population structure can also produce peak-and-decline patterns.','')
    curves <- read.delim(file.path(args[[3]],'curves.tsv'),check.names=FALSE)
    if(!all(c('fit','kind','time','population','scaling') %in% names(curves)) ||
       any(!is.finite(curves$time)) || any(!is.finite(curves$population))) stop('Invalid exported curve')
    # Time zero is retained in the download, and omitted only from the log-axis display.
    display <- curves[curves$time>0 & curves$population>0,]
    if(!nrow(display) || !any(display$kind=='observed')) stop('No positive observed curve to plot')
    figure <- ggplot() +
        geom_step(data=display[display$kind=='bootstrap',],aes(time,population,group=fit),
                  color='grey65',alpha=.45,linewidth=.4) +
        geom_step(data=display[display$kind=='observed',],aes(time,population),
                  color='#0072B2',linewidth=.9) +
        scale_x_log10() + scale_y_continuous(limits=c(0,NA),expand=expansion(mult=c(0,.06))) +
        theme_bw(base_size=12) +
        labs(title=context$sample,subtitle=paste0(context$mode,'; ',context$bootstrap_count,
             ' block bootstrap replicates; iteration ',context$iteration),
             x=if(scaled) 'Time before present (years)' else 'Scaled divergence (t × theta per base)',
             y=if(scaled) 'Effective population size (Ne)' else 'Scaled size (lambda × theta per base)',
             caption='Blue: observed fit. Grey: block-resampled fits. Bootstrap spread is not parameter or model certainty.')
    ggsave(file.path(out,'history.png'),figure,width=8,height=5,dpi=180,bg='white')
    ggsave(file.path(out,'history.pdf'),figure,width=8,height=5,bg='white')
    for(name in c('curves.tsv','fit_manifest.json','interval_qc.tsv'))
        if(!file.copy(file.path(args[[3]],name),file.path(out,name),overwrite=TRUE)) stop('Failed to retain curve evidence')
    md <- c(md,'![PSMC history](history.png)','',
        '[PDF](history.pdf) · [Curve table](curves.tsv) · [Fit/iteration/scaling manifest](fit_manifest.json)','',
        '[Observed fit](../fits/observed.psmc) · [Bootstrap fits and resampling records](../fits/)','')
    interval_qc <- read.delim(file.path(out,'interval_qc.tsv'))
    low <- sum(interval_qc$fit=='observed' & interval_qc$below_upstream_advisory_ten=='yes')
    md <- c(md,paste0('Observed fit: **',low,' free-parameter groups below the upstream advisory of approximately ten expected recombinations**. This is a model-information diagnostic, not an automatic fit rejection.'),
        '', '[Per-group interval information](interval_qc.tsv)','')
    if(scaled) md <- c(md,paste0('Mutation rate: ',context$mutation_rate,
        ' per base per generation; generation time: ',context$generation_time,' years.'),'')
    else md <- c(md,'Mutation rate and generation time are unset. These axes are explicitly unscaled; no species-specific values were assumed.','')
} else md <- c(md,paste0('**No fit produced:** ',qc$reason),
    'The masked consensus and QC remain available; this outcome is not a demographic estimate.','')
md <- c(md,'## Scope and limitations','',qc$scope_note,'',
    'Only long reference scaffolds are included. This does not establish autosomal identity. Gaps, excluded intervals, complex-variant flanks and uncertain calls remain unknown. No unknown position is filled as homozygous reference.','',
    if(context$mode=='assemblies')
      'Dipcall confident regions define the assembly-based callability mask. Assembly phasing and collapsed sequence can affect heterozygosity.' else
      'Complete same-individual HiFi reads were mapped and called with the general bcftools diploid model. This route needs biological sensitivity validation; mapping/calling thresholds do not guarantee unbiased heterozygote detection.',
    '', 'The default pattern estimates the four youngest atomic intervals separately. Review recombination support when changing the pattern. The main history plot omits the terminal free-parameter block by default.',
    '', '[Per-contig QC](../input/contig_qc.tsv) · [Input QC and filters](../input/input_qc.json)',
    if(qc$callable_bp>0) '[Callable BED](../input/callable.bed)' else 'Callable BED is empty: no callable positions.',
    if(file.info(file.path(args[[2]],'input.psmcfa'))$size>0) '[PSMC input](../input/input.psmcfa)' else 'PSMC input is empty: no sequence passed the input filters.',
    '', '[Masked diploid consensus](../input/diploid.fa.gz) · [Genotype/reference provenance](../input/provenance.json)')
writeLines(md,file.path(out,'report.md'))
write_json(summary,file.path(out,'summary.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
