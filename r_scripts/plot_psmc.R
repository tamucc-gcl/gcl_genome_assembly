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
if(ready) {
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
    '', 'Time-interval settings are configurable and require review of inferred recombination information. The default pattern is an upstream example, not an organism-specific optimization. The terminal free-parameter block is omitted from plots by default.',
    '', '[Per-contig QC](../input/contig_qc.tsv) · [Input QC and filters](../input/input_qc.json)',
    if(qc$callable_bp>0) '[Callable BED](../input/callable.bed)' else 'Callable BED is empty: no callable positions.',
    if(file.info(file.path(args[[2]],'input.psmcfa'))$size>0) '[PSMC input](../input/input.psmcfa)' else 'PSMC input is empty: no sequence passed the input filters.',
    '', '[Masked diploid consensus](../input/diploid.fa.gz) · [Genotype/reference provenance](../input/provenance.json)')
writeLines(md,file.path(out,'report.md'))
write_json(summary,file.path(out,'summary.json'),auto_unbox=TRUE,pretty=TRUE,null='null')
