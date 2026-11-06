# Diagnostic presentation only; all coverage/filter decisions come from retained tables.
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop('Expected reports directory')
out <- args[[1]]
d <- read.delim(file.path(out,'classification_summary.tsv'),check.names=FALSE)
if(!nrow(d)) quit(status=0)
d$label <- factor(d$label,levels=unique(d$label))
d$haplotype <- factor(d$haplotype,levels=c('hap2','hap1'))
palette <- c('#BDBDBD','#969696','#756BB1','#E6AB02','#A6761D',
             '#E7298A','#D95F02','#66A61E','#1B9E77')
p <- ggplot(d,aes(haplotype,percent,fill=label)) +
    geom_col(width=.65) + coord_flip(ylim=c(0,100)) + facet_wrap(~sample,ncol=1) +
    scale_fill_manual(values=palette,name=NULL,drop=FALSE,
        labels=function(v) vapply(v,function(x) paste(strwrap(x,width=30),collapse='\n'),character(1))) +
    scale_y_continuous(breaks=seq(0,100,20),expand=expansion(mult=c(0,0))) +
    theme_bw(base_size=11) + theme(legend.position='bottom',
        legend.text=element_text(size=9),plot.caption=element_text(hjust=0,size=9),
        plot.margin=margin(12,22,12,12),
        panel.grid.major.y=element_blank(),strip.text=element_text(hjust=0)) +
    guides(fill=guide_legend(ncol=2,byrow=TRUE)) +
    labs(title='Why Dipcall alignment spans do or do not support a confident mask',
         subtitle='Hap1 is the reference-to-itself control; categories partition selected non-gap scope',
         x=NULL,y='Selected scope (%)',
         caption=paste(strwrap('Overlap and low-quality categories describe observed alignments; they do not establish biological absence or repeat identity.',width=110),collapse='\n'))
height <- 1.3*length(unique(d$sample))+4
ggsave(file.path(out,'mask_reasons.png'),p,width=11,height=height,dpi=180,bg='white')
ggsave(file.path(out,'mask_reasons.pdf'),p,width=11,height=height,bg='white')
write.table(data.frame(tool=c('R','ggplot2'),version=c(R.version.string,as.character(packageVersion('ggplot2')))),
    file.path(out,'plot_versions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
