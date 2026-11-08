#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(jsonlite))
root <- commandArgs(trailingOnly=TRUE)[[1]]
context <- fromJSON(file.path(root,'comparison.json'))
d <- read.delim(file.path(root,'curves.tsv'))
d <- d[d$kind=='observed' & d$time>0,]
d$run <- factor(d$run,levels=c('asm20_50kb','asm20_20kb','asm20_10kb'))
p <- ggplot(d,aes(time,population,colour=run)) + geom_step(linewidth=.8) + scale_x_log10() +
 scale_y_continuous(limits=c(0,max(.0001,max(d$population)*1.1)),expand=c(0,0)) +
 scale_colour_manual(values=c(asm20_50kb='#0072B2',asm20_20kb='#D55E00',asm20_10kb='#009E73'),drop=FALSE) +
 theme_bw(base_size=12) + theme(legend.position='top',plot.margin=margin(12,16,12,24)) +
 labs(title=paste(context$sample,'· PSMC sensitivity to minimum alignment length'),
 subtitle='asm20 throughout; observed fits; confidence masks and variant filters changed together',
 x='Scaled divergence (t × theta per base)',y='Scaled size (lambda × theta per base)',colour=NULL,
 caption='Callable scope differs between runs. Read-support and genotype-agreement audits are linked separately.')
ggsave(file.path(root,'length_histories.png'),p,width=10,height=6,dpi=160,bg='white')
ggsave(file.path(root,'length_histories.pdf'),p,width=10,height=6,bg='white')
write.table(data.frame(tool=c('R','ggplot2','jsonlite'),version=c(as.character(getRversion()),as.character(packageVersion('ggplot2')),as.character(packageVersion('jsonlite')))),
 file.path(root,'plot_versions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
