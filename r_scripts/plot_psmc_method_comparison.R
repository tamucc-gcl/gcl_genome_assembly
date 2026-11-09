#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(jsonlite)})
root <- commandArgs(trailingOnly=TRUE)[[1]]
context <- fromJSON(file.path(root,'comparison.json'))
d <- read.delim(file.path(root,'method_curves.tsv'))
d <- d[d$kind=='observed' & d$time>0,]
if(length(unique(d$scaling))!=1L) stop('Cannot overlay methods with different scaling conventions')
scaled <- !identical(unique(d$scaling),'unscaled')
p <- ggplot(d,aes(time,population,colour=run)) + geom_step(linewidth=.8) + scale_x_log10() +
 scale_y_continuous(limits=c(0,max(.0001,max(d$population)*1.1)),expand=c(0,0)) +
 scale_colour_manual(values=c(cohort='#009E73',control_20kb='#D55E00',hifi='#0072B2')) +
 labs(title=paste(context$sample,'· PSMC method sensitivity'),
 subtitle='Observed fits on the callable sequence retained by each method',
 x=if(scaled) 'Time before present (years)' else 'Scaled divergence (t × theta per base)',
 y=if(scaled) 'Effective population size (Ne)' else 'Scaled size (lambda × theta per base)',colour=NULL,
 caption='Different masks and genotype ascertainment can change histories. This is not a matched-scope fit comparison.') +
 theme_bw(base_size=12) + theme(legend.position='top',plot.margin=margin(12,16,12,24))
ggsave(file.path(root,'method_histories.png'),p,width=10,height=6,dpi=160,bg='white')
ggsave(file.path(root,'method_histories.pdf'),p,width=10,height=6,bg='white')
write.table(data.frame(tool=c('R','ggplot2','jsonlite'),version=c(as.character(getRversion()),as.character(packageVersion('ggplot2')),as.character(packageVersion('jsonlite')))),
 file.path(root,'plot_versions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
