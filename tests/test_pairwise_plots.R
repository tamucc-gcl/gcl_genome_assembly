#!/usr/bin/env Rscript
# Real pafr/ggplot2/gggenomes regression checks; no plotting executable doubles.
suppressPackageStartupMessages({library(ggplot2);library(pafr);library(gggenomes)})
args <- commandArgs(TRUE)
source(file.path(args[1],"r_scripts","pairwise_plot_helpers.R"))
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
ref <- data.frame(seq_id=c("chr1","chr2","unaligned"),length=c(1e7,1e7,1e7))
query <- data.frame(seq_id=c("chr1","chr2"),length=c(1e7,1e7))
path <- file.path(args[2],"fixture.paf.gz")
con <- gzfile(path,"wt")
writeLines(c("chr1\t10000000\t0\t5000000\t+\tchr1\t10000000\t0\t5000000\t4900000\t5000000\t60",
             "chr2\t10000000\t1000000\t9000000\t-\tchr2\t10000000\t1000000\t9000000\t7800000\t8000000\t60\tdv:f:0.02"),con)
close(con)
d <- read_final_paf(path,ref,query);stopifnot(nrow(d)==2)
r <- sequence_layout(ref,0);q <- sequence_layout(query,0)
p <- dotplot_final(d,r,q,"Forward and inverted mappings","Synthetic layout validation")
b <- ggplot_build(p)
stopifnot(max(b$layout$panel_params[[1]]$y.range)>=3e7,
          max(b$layout$panel_params[[1]]$x.range)>=2e7,
          b$data[[2]]$x[1]>b$data[[2]]$xend[1],
          b$data[[2]]$y[1]<b$data[[2]]$yend[1])
ggsave(file.path(args[2],"dotplot.png"),p,width=9,height=8,dpi=180,bg="white")
p <- riparian_final(d,r,q,"reference","query","Forward and inverted mappings","Synthetic layout validation",.45)
b <- ggplot_build(p)
links <- Filter(function(x) all(c("xmin","xmax","fill") %in% names(x)),b$data)[[1]]
stopifnot(nrow(links)==2,all(is.finite(links$x)),any(links$xmin>links$xmax))
ggsave(file.path(args[2],"riparian.png"),p,width=12,height=6,dpi=180,bg="white")
stopifnot(plot_status("dotplot",d[FALSE,],r,q)$status=="no_alignments_passing_plot_filters")
ggsave(file.path(args[2],"empty.png"),empty_comparison_plot("no_alignments_passing_plot_filters","Empty","Measured empty"),width=5,height=4)
bad <- ref;bad$length[1] <- 999
stopifnot(inherits(try(read_final_paf(path,bad,query),silent=TRUE),"try-error"))
con <- gzfile(path,"wt");writeLines("bad PAF",con);close(con)
stopifnot(inherits(try(read_final_paf(path,ref,query),silent=TRUE),"try-error"))
writeBin(as.raw(c(31,139,8,0,0,0)),path)
stopifnot(inherits(try(read_final_paf(path,ref,query),silent=TRUE),"try-error"))
# Many scaffolds must render without hundreds of text labels or grid boundaries.
many <- data.frame(seq_id=paste0("scaffold_",seq_len(60)),length=rep(1e6,60))
line <- "scaffold_1\t1000000\t0\t500000\t+\tscaffold_1\t1000000\t0\t500000\t490000\t500000\t60"
con <- gzfile(path,"wt");writeLines(line,con);close(con)
d <- read_final_paf(path,many,many);m <- sequence_layout(many,0)
ggsave(file.path(args[2],"fragmented-dotplot.png"),dotplot_final(d,m,m,"Fragmented","Synthetic layout validation"),width=8,height=8)
ggsave(file.path(args[2],"fragmented-riparian.png"),riparian_final(d,m,m,"reference","query","Fragmented","Synthetic layout validation",.45),width=12,height=6)
cat("PASS: real R plots, full FAI axes, inversion, empty/malformed/corrupt input, fragmented layouts\n")
