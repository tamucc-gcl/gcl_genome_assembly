#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(jsonlite))
root <- commandArgs(trailingOnly=TRUE)[[1]]
context <- fromJSON(file.path(root,'summary.json'))
d <- read.delim(file.path(root,'windows.tsv'))
shared <- d[d$region=='shared' & d$bp>0,]
shared$discordance <- 100*shared$discordant/shared$bp
shared$position <- (shared$start+shared$end)/2/1e6
shared$scaffold <- factor(shared$scaffold,levels=unique(shared$scaffold))
shared$window_count <- ave(shared$start,shared$scaffold,FUN=length)
p <- ggplot(shared,aes(position,discordance)) +
 geom_line(data=shared[shared$window_count>1,],colour='#0072B2',linewidth=.5) +
 geom_point(colour='#0072B2',size=.6) +
 scale_y_continuous(limits=c(0,max(.01,max(shared$discordance)*1.05)),expand=c(0,0)) +
 facet_wrap(~scaffold,scales='free_x',ncol=3) + theme_bw(base_size=11) +
 labs(title='Genotype disagreement within shared callable sequence',
      subtitle=paste(context$labels$left,'versus',context$labels$right,'· masked PSMC input · 1 Mb windows'),
      x='Scaffold position (Mb)',y='Discordant calls / shared callable bases (%)')
ggsave(file.path(root,'reliability.png'),p,width=10,height=8,dpi=160,bg='white')
ggsave(file.path(root,'reliability.pdf'),p,width=10,height=8,bg='white')
write.table(data.frame(tool=c('R','ggplot2','jsonlite'),version=c(as.character(getRversion()),as.character(packageVersion('ggplot2')),as.character(packageVersion('jsonlite')))),
 file.path(root,'plot_versions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
