#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2);library(pafr);library(gggenomes);library(argparse)})
script <- sub("^--file=","",grep("^--file=",commandArgs(FALSE),value=TRUE)[1])
source(file.path(dirname(normalizePath(script)),"pairwise_plot_helpers.R"))
p <- ArgumentParser(description="Final-assembly pafr dotplot and gggenomes riparian plot")
for (name in c("paf","ref-fai","query-fai","ref","query","prefix")) p$add_argument(paste0("--",name),required=TRUE)
for (name in c("dot-min-seq","ribbon-min-seq","ribbon-min-aln","upstream-min-aln"))
  p$add_argument(paste0("--",name),type="double",default=if(grepl("seq",name))1000000 else 150000)
p$add_argument("--upstream-min-mapq",type="integer",default=5)
for (name in c("dot-width","dot-height","ribbon-width","ribbon-height"))
  p$add_argument(paste0("--",name),type="double",default=10)
p$add_argument("--ribbon-alpha",type="double",default=.45)
a <- p$parse_args()
tryCatch({
  if (min(a$dot_min_seq,a$ribbon_min_seq,a$ribbon_min_aln,a$upstream_min_aln)<0 ||
      min(a$dot_width,a$dot_height,a$ribbon_width,a$ribbon_height)<=0 ||
      a$ribbon_alpha<0 || a$ribbon_alpha>1 || a$upstream_min_mapq<0 || a$upstream_min_mapq>254)
    stop("Invalid plotting thresholds/dimensions")
  ref <- read_final_fai(a$ref_fai);query <- read_final_fai(a$query_fai)
  d <- read_final_paf(a$paf,ref,query)
  rows <- list()
  for (kind in c("dotplot","ribbon")) {
    minimum <- if(kind=="dotplot")a$dot_min_seq else a$ribbon_min_seq
    ml <- if(kind=="dotplot")a$upstream_min_aln else max(a$upstream_min_aln,a$ribbon_min_aln)
    r <- sequence_layout(ref,minimum);q <- sequence_layout(query,minimum)
    keep_mq <- if(a$upstream_min_mapq==0) rep(TRUE,nrow(d)) else d$mapq>=a$upstream_min_mapq & d$mapq!=255
    chosen <- d[d$tname %in% r$seq_id & d$qname %in% q$seq_id & d$alen>=ml & keep_mq,,drop=FALSE]
    subtitle <- paste0(format(nrow(chosen),big.mark=",")," mapping chains; sequences >= ",format(minimum,big.mark=","),
                       " bp; effective chain length >= ",format(ml,big.mark=",")," bp; MAPQ >= ",a$upstream_min_mapq)
    title <- paste(a$query,"vs",a$ref)
    if(kind=="dotplot") {
      plot <- dotplot_final(chosen,r,q,title,subtitle)+labs(x=paste0(a$query," (final sequences; numeric coordinates in Mb)"),
                                                          y=paste0(a$ref," (final sequences; numeric coordinates in Mb)"))
      output <- paste0(a$prefix,"_dotplot.png");width <- a$dot_width;height <- a$dot_height
    } else {
      plot <- riparian_final(chosen,r,q,a$ref,a$query,title,subtitle,a$ribbon_alpha)
      output <- paste0(a$prefix,"_riparian.png");width <- a$ribbon_width;height <- a$ribbon_height
    }
    ggsave(output,plot,width=width,height=height,dpi=300,bg="white")
    rows[[kind]] <- plot_status(kind,chosen,r,q)
  }
  write.table(do.call(rbind,rows),paste0(a$prefix,".plot_status.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
},error=function(e){message(conditionMessage(e));quit(status=65)})
