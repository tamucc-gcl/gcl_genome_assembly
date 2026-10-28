# Final-assembly dictionaries and strict PAF parsing shared by the R plots.
read_final_fai <- function(path) {
  x <- read.delim(path, header=FALSE, sep="\t", quote="", comment.char="",
                  colClasses="character", check.names=FALSE)
  if (ncol(x)<2 || !nrow(x)) stop("Empty or invalid final assembly index")
  d <- data.frame(seq_id=x[[1]], length=suppressWarnings(as.numeric(x[[2]])))
  if (any(!nzchar(d$seq_id)) || anyDuplicated(d$seq_id) ||
      any(!is.finite(d$length) | d$length<=0 | d$length!=floor(d$length)))
    stop("Invalid or duplicate final assembly index entry")
  d
}

read_final_paf <- function(path, ref, query) {
  con <- if (grepl("\\.gz$",path)) gzfile(path,"rt") else file(path,"rt")
  on.exit(close(con))
  # Corrupt gzip warnings must fail rather than becoming an empty measurement.
  lines <- withCallingHandlers(readLines(con,warn=FALSE), warning=function(w) stop(conditionMessage(w)))
  lines <- lines[nzchar(trimws(lines))]
  fields <- strsplit(lines,"\t",fixed=TRUE)
  if (any(lengths(fields)<12)) stop("PAF record has fewer than 12 columns")
  nm <- c("qname","qlen","qstart","qend","strand","tname","tlen","tstart","tend","nmatch","alen","mapq")
  if (!length(fields)) {
    d <- as.data.frame(setNames(replicate(12,character(),simplify=FALSE),nm))
  } else {
    d <- as.data.frame(do.call(rbind,lapply(fields,head,12)),stringsAsFactors=FALSE)
    names(d) <- nm
  }
  for (key in setdiff(nm,c("qname","tname","strand"))) {
    d[[key]] <- suppressWarnings(as.numeric(d[[key]]))
    if (any(!is.finite(d[[key]]) | d[[key]]!=floor(d[[key]]))) stop("Invalid PAF integer")
  }
  if (any(!d$strand %in% c("+","-")) || any(!nzchar(d$qname) | !nzchar(d$tname)) ||
      any(d$qstart<0 | d$qend<=d$qstart | d$qend>d$qlen |
          d$tstart<0 | d$tend<=d$tstart | d$tend>d$tlen |
          d$alen<=0 | d$nmatch<0 | d$nmatch>d$alen | d$mapq<0 | d$mapq>255)) stop("Invalid PAF coordinates/counts")
  ri <- match(d$tname,ref$seq_id);qi <- match(d$qname,query$seq_id)
  if (anyNA(ri) || anyNA(qi) || any(ref$length[ri]!=d$tlen | query$length[qi]!=d$qlen))
    stop("PAF sequence names/lengths disagree with final assembly indexes")
  d
}

sequence_layout <- function(d,minimum) {
  d <- d[d$length>=minimum,,drop=FALSE]
  chr <- suppressWarnings(as.numeric(sub("^chr([0-9]+).*$","\\1",d$seq_id)))
  chr[!grepl("^chr[0-9]+",d$seq_id)] <- Inf
  d <- d[order(chr,-d$length,d$seq_id),,drop=FALSE]
  d$offset <- cumsum(d$length)-d$length
  d$midpoint <- d$offset+d$length/2
  d
}

plot_status <- function(kind,chosen,ref,query) {
  state <- if (!nrow(ref) || !nrow(query)) "no_sequences_above_threshold" else
    if (!nrow(chosen)) "no_alignments_passing_plot_filters" else "plotted"
  data.frame(plot=kind,status=state,alignments=nrow(chosen),ref_sequences=nrow(ref),
             query_sequences=nrow(query),ref_axis_bp=sum(ref$length),query_axis_bp=sum(query$length))
}

empty_comparison_plot <- function(state,title,subtitle) {
  ggplot2::ggplot()+ggplot2::annotate("text",x=.5,y=.5,label=gsub("_"," ",state))+
    ggplot2::theme_void()+ggplot2::labs(title=title,subtitle=subtitle)
}

dotplot_final <- function(d,ref,query,title,subtitle) {
  if (!nrow(d)) return(empty_comparison_plot(plot_status("dotplot",d,ref,query)$status,title,subtitle))
  # pafr draws the real mapping records. Replace its aligned-only concatenation
  # with offsets from the full final FAI, without fabricating alignment records.
  qord <- query$seq_id[query$seq_id %in% d$qname]
  rord <- ref$seq_id[ref$seq_id %in% d$tname]
  p <- pafr::dotplot(d,order_by="provided",ordering=list(qord,rord),
                    label_seqs=FALSE,dashes=FALSE,line_size=.4)
  for (i in seq_along(p$layers)) {
    x <- p$layers[[i]]$data
    if (!is.data.frame(x) || !all(c("qname","tname","concat_qstart") %in% names(x))) next
    x$concat_qstart <- x$qstart+query$offset[match(x$qname,query$seq_id)]
    x$concat_qend <- x$qend+query$offset[match(x$qname,query$seq_id)]
    x$concat_tstart <- x$tstart+ref$offset[match(x$tname,ref$seq_id)]
    x$concat_tend <- x$tend+ref$offset[match(x$tname,ref$seq_id)]
    p$layers[[i]]$data <- x
  }
  # Dense fragmented inputs get numeric axes, not hundreds of grid lines/labels.
  if (nrow(query)<=40) {
    p <- p+ggplot2::geom_vline(xintercept=query$offset,colour="grey85",linewidth=.2)+
      ggplot2::scale_x_continuous(breaks=query$midpoint,labels=query$seq_id,limits=c(0,sum(query$length)),expand=c(0,0))
  } else p <- p+ggplot2::scale_x_continuous(labels=function(x)x/1e6,limits=c(0,sum(query$length)),expand=c(0,0))
  if (nrow(ref)<=40) {
    p <- p+ggplot2::geom_hline(yintercept=ref$offset,colour="grey85",linewidth=.2)+
      ggplot2::scale_y_continuous(breaks=ref$midpoint,labels=ref$seq_id,limits=c(0,sum(ref$length)),expand=c(0,0))
  } else p <- p+ggplot2::scale_y_continuous(labels=function(x)x/1e6,limits=c(0,sum(ref$length)),expand=c(0,0))
  p+ggplot2::theme_classic(base_size=11)+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=60,hjust=1,size=7),
      axis.text.y=ggplot2::element_text(size=7),panel.border=ggplot2::element_rect(fill=NA,colour="black"))+
    ggplot2::labs(title=title,subtitle=subtitle,x="Query final sequences",y="Reference final sequences",
      caption="Mapping-chain synteny overview; unaligned final sequences remain on the axes.")
}

riparian_final <- function(d,ref,query,ref_id,query_id,title,subtitle,alpha) {
  if (!nrow(d)) return(empty_comparison_plot(plot_status("ribbon",d,ref,query)$status,title,subtitle))
  # Prefix sequence identifiers internally: scaffold names can repeat between bins.
  # Display labels retain the original final assembly names.
  seqs <- rbind(data.frame(bin_id=ref_id,seq_id=paste0("ref:",ref$seq_id),length=ref$length,label=ref$seq_id),
                data.frame(bin_id=query_id,seq_id=paste0("qry:",query$seq_id),length=query$length,label=query$seq_id))
  links <- data.frame(bin_id=ref_id,seq_id=paste0("ref:",d$tname),start=d$tstart+1,end=d$tend,
                      bin_id2=query_id,seq_id2=paste0("qry:",d$qname),
                      start2=ifelse(d$strand=="+",d$qstart+1,d$qend),
                      end2=ifelse(d$strand=="+",d$qend,d$qstart+1),ref_chr=d$tname)
  palette <- setNames(grDevices::hcl.colors(nrow(ref),"Dark 3"),ref$seq_id)
  p <- gggenomes::gggenomes(seqs=tibble::as_tibble(seqs),links=tibble::as_tibble(links),spacing=0)+
    gggenomes::geom_link(ggplot2::aes(fill=ref_chr),alpha=alpha,colour="grey40",linewidth=.05)+
    gggenomes::geom_seq(linewidth=2,colour="grey55")+
    gggenomes::geom_bin_label(size=3)+
    ggplot2::scale_fill_manual(values=palette,name="Reference sequence",breaks=ref$seq_id,
                              guide=if(nrow(ref)<=30)ggplot2::guide_legend(ncol=2) else "none")
  if (max(nrow(ref),nrow(query))<=40)
    p <- p+gggenomes::geom_seq_label(ggplot2::aes(label=label),size=2.2,vjust=-.8,check_overlap=TRUE)
  p+ggplot2::labs(title=title,subtitle=subtitle,
      caption=paste0("minimap2 mapping chains; gggenomes ribbons. ",nrow(ref)," reference / ",nrow(query),
                     " query sequences. Neutral bars do not imply chromosome identity."))+
    ggplot2::theme(plot.background=ggplot2::element_rect(fill="white",colour=NA),
                    panel.background=ggplot2::element_rect(fill="white",colour=NA),legend.text=ggplot2::element_text(size=7))
}
