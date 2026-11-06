# Presentation only: combine already exported curves without changing evidence or fits.
suppressPackageStartupMessages({library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop('Expected PSMC index output directory')
out <- args[[1]]
path <- file.path(out,'combined_curves.tsv')
if(!file.exists(path)) quit(status=0) # All individuals skipped or not fitted.
d <- read.delim(path,check.names=FALSE)
s <- read.delim(file.path(out,'summary.tsv'),check.names=FALSE)
required <- c('sample','taxid','reference','kind','fit','time','population','scaling')
if(!all(required %in% names(d)) || any(!is.finite(d$time)) || any(!is.finite(d$population)))
    stop('Invalid combined curve evidence')
# PSMC histories are individual summaries, not reference-coordinate plots.
# Different own-reference IDs must not split comparable individuals into panels.
d$comparison <- paste0('Taxid ',d$taxid,' | ',d$scaling)
display <- d[d$time>0 & d$population>0,]
observed <- display[display$kind=='observed',]
if(!nrow(observed)) stop('No positive observed curves in completed results')
not_fitted <- s$sample[s$status!='completed']
omitted <- if(length(not_fitted)) paste0('No inferred curve: ',paste(not_fitted,collapse=', '),'. ') else ''
unscaled <- unique(observed$scaling)
xlab <- if(identical(unscaled,'unscaled')) 'Scaled divergence (t × theta per base)' else
    if(identical(unscaled,'years_Ne')) 'Time before present (years)' else 'Time (units shown in panel strip)'
ylab <- if(identical(unscaled,'unscaled')) 'Scaled size (lambda × theta per base)' else
    if(identical(unscaled,'years_Ne')) 'Effective population size (Ne)' else 'Population size (units shown in panel strip)'
common <- theme_bw(base_size=11) + theme(legend.position='bottom',
    plot.caption=element_text(hjust=0,size=9),strip.text=element_text(size=9))
caption <- paste(strwrap(paste0(omitted,
    'Exploratory histories: check callable scope and model information. Zero-time rows remain in the evidence table.'),width=115),collapse='\n')
overlay <- ggplot(observed,aes(time,population,color=sample,group=sample)) +
    geom_step(linewidth=.85) + scale_color_viridis_d(end=.85,name='Individual') +
    scale_x_log10() + scale_y_log10() + facet_wrap(~comparison,scales='free') + common +
    labs(title='PSMC histories across individuals',subtitle='Observed fits; both axes logarithmic',
         x=xlab,y=ylab,caption=caption)
ggsave(file.path(out,'combined_history.png'),overlay,width=10,height=6,dpi=180,bg='white')
ggsave(file.path(out,'combined_history.pdf'),overlay,width=10,height=6,bg='white')
display$individual <- paste(display$sample,display$comparison,sep='\n')
panels <- ggplot(display,aes(time,population,group=interaction(sample,fit))) +
    geom_step(data=display[display$kind=='bootstrap',],color='grey60',alpha=.35,linewidth=.35) +
    geom_step(data=display[display$kind=='observed',],aes(color=sample),linewidth=.85) +
    scale_color_viridis_d(end=.85,guide='none') + scale_x_log10() +
    scale_y_log10() +
    facet_wrap(~individual,ncol=2,scales='free') + common +
    labs(title='PSMC histories and block bootstraps',
         subtitle='Coloured: observed; grey: bootstraps. Both axes logarithmic; panel scales differ.',
         x=xlab,y=ylab,caption=caption)
height <- 2.8*ceiling(length(unique(display$individual))/2)+1.6
ggsave(file.path(out,'individual_histories.png'),panels,width=11,height=height,dpi=180,bg='white')
ggsave(file.path(out,'individual_histories.pdf'),panels,width=11,height=height,bg='white')
