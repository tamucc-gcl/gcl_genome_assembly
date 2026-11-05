#!/usr/bin/env Rscript
# Presentation only: all tracks and profiles are retained assay measurements.
suppressPackageStartupMessages({library(ggplot2);library(jsonlite);library(patchwork)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3) stop('Expected packet, event ID and selected-cut JSON')
packet <- args[1];id <- args[2]
data <- fromJSON(file.path(packet,'evidence.json'),simplifyVector=FALSE)
e <- Filter(function(x) x$id==id,data$events)[[1]]
cuts <- fromJSON(args[3],simplifyVector=FALSE)
# Explicit converters keep strings and null counts separate.
numeric_rows <- function(values,fields) {
  as.data.frame(setNames(lapply(fields,function(key) vapply(values,function(x)
    if(is.null(x[[key]])) NA_real_ else as.numeric(x[[key]]),numeric(1))),fields))
}
tracks <- do.call(rbind,lapply(seq_along(e$peer_tracks),function(i) {
  blocks <- Filter(function(b) !is.null(b$chrom) && nzchar(b$chrom),e$peer_tracks[[i]]$blocks)
  if(!length(blocks)) return(NULL)
  data.frame(row=i,chrom=vapply(blocks,`[[`,character(1),'chrom'),
    lo=vapply(blocks,`[[`,numeric(1),'lo'),hi=vapply(blocks,`[[`,numeric(1),'hi'))
}))
if(is.null(tracks)) tracks <- data.frame(row=numeric(),chrom=character(),lo=numeric(),hi=numeric())
regions <- do.call(rbind,lapply(e$chromosome_regions,function(r) data.frame(chrom=r$chrom,lo=r$lo,hi=r$hi)))
if(is.null(regions)) regions <- data.frame(chrom=character(),lo=numeric(),hi=numeric())
chrom_order <- function(values) {
  number <- suppressWarnings(as.numeric(sub('^chr','',sub('_.*$','',values))))
  values[order(is.na(number),number,values)]
}
focal <- chrom_order(unique(c(e$left,e$right)))
chroms <- c(focal,setdiff(chrom_order(unique(c(tracks$chrom,regions$chrom))),focal))
seeds <- suppressWarnings(as.numeric(sub('^chr','',sub('_.*$','',chroms))))
for(i in which(is.na(seeds))) seeds[i] <- sum(utf8ToInt(chroms[i]))
colours <- setNames(hsv((seeds*.61803398875)%%1,.25,.65),chroms)
colours[focal] <- c('#0072B2','#D55E00')[seq_along(focal)]
fill_scale <- function() scale_fill_manual(values=colours,limits=chroms,drop=FALSE,
  name='Chromosome assignments (focal pair first)',guide=guide_legend(nrow=2,byrow=TRUE))
profile <- function(name) {
  value <- e[[paste0(name,'_profile')]]
  if(is.null(value)) value <- e[[name]]
  numeric_rows(value$probes,c('position','molecules'))
}
normal <- profile('hifi');strict <- profile('hifi_strict');allprimary <- profile('hifi_all_primary')
depth <- numeric_rows(e$depth,c('start','end','depth'))
depth_all <- numeric_rows(e$depth_all_primary,c('start','end','depth'))
contacts <- do.call(rbind,lapply(names(e$contact_profiles),function(lib) {
  x <- numeric_rows(e$contact_profiles[[lib]],c('position','ratio'))
  if(nrow(x)) x$library <- lib
  x
}))
if(is.null(contacts)) contacts <- data.frame(position=numeric(),ratio=numeric(),library=character())
libraries <- sort(unique(contacts$library))
library_colours <- c('#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9','#000000','#882255')
if(length(libraries)>length(library_colours)) library_colours <- hcl.colors(length(libraries),'Dark 3')
library_colours <- setNames(library_colours[seq_along(libraries)],libraries)
telomeres <- numeric_rows(e$telomeres,c('start','count'))
alternatives <- Filter(function(x) is.null(x$site_verdict)||x$site_verdict!='EXCLUDE',e$alternatives)
alternatives <- head(alternatives,6)
gaps <- Filter(function(x) is.null(x$site_verdict)||x$site_verdict!='EXCLUDE',e$cut_options)
markers <- rbind(
  if(length(alternatives)) data.frame(point=vapply(alternatives,`[[`,numeric(1),'probe_bp'),
    label=vapply(alternatives,`[[`,character(1),'id'),kind='R assay (not a cut)'),
  if(length(gaps)) data.frame(point=vapply(gaps,`[[`,numeric(1),'cut_bp'),
    label=vapply(gaps,`[[`,character(1),'id'),kind='Exact gap option'),
  if(length(cuts)) data.frame(point=vapply(cuts,function(x) as.numeric(x$cut_bp),numeric(1)),
    label='Selected cut',kind='Submitted cut'))
if(is.null(markers)) markers <- data.frame(point=numeric(),label=character(),kind=character())
marker_colours <- c('R assay (not a cut)'='#762A83','Exact gap option'='#B2182B','Submitted cut'='#007A3D')
base <- theme_bw(base_size=10)+theme(panel.grid.minor=element_blank(),
  legend.position='top',plot.background=element_rect(fill='white',colour=NA),
  plot.title=element_text(size=10),legend.text=element_text(size=8),legend.title=element_text(size=8))
save_figure <- function(figure,stem,width,height) {
  ggsave(paste0(stem,'.png'),figure,width=width,height=height,dpi=150,bg='white',limitsize=FALSE)
  ggsave(paste0(stem,'.pdf'),figure,width=width,height=height,bg='white',device=cairo_pdf,limitsize=FALSE)
}
detail_window <- function() {
  candidates <- alternatives
  focus <- e$recommended_cut
  if(is.null(focus)) focus <- e$preferred_region
  if(is.null(focus)) focus <- if(length(candidates)) candidates[[1]] else e
  selected <- vapply(cuts,function(x) as.numeric(x$cut_bp),numeric(1))
  point <- if(length(selected)) selected[1] else if(!is.null(focus$cut_bp)&&focus$cut_bp!=0) focus$cut_bp else
    if(!is.null(focus$probe_bp)&&focus$probe_bp!=0) focus$probe_bp else (focus$lo+focus$hi)%/%2
  lo <- hi <- point
  if(!length(selected)&&is.null(focus$cut_bp)&&focus$hi-focus$lo<=200000) {lo<-focus$lo;hi<-focus$hi}
  c(max(0,lo-100000),hi+100000)
}
axis_pair <- function(primary,secondary) {
  left <- max(1,primary,na.rm=TRUE);right <- max(1,secondary,na.rm=TRUE)
  list(factor=left/right,limit=left*1.08)
}
span_scale <- axis_pair(c(normal$molecules,strict$molecules),allprimary$molecules)
depth_scale <- axis_pair(depth$depth,depth_all$depth)
for(detail in c(FALSE,TRUE)) {
  window <- if(detail) detail_window() else c(max(0,e$lo-500000),e$hi+500000)
  decorate <- function(plot,shade=TRUE) {
    if(shade&&nrow(regions)) plot <- plot+geom_rect(data=regions,
      aes(xmin=lo/1e6,xmax=hi/1e6,ymin=-Inf,ymax=Inf,fill=chrom),inherit.aes=FALSE,alpha=.1,show.legend=FALSE)
    plot <- plot+annotate('rect',xmin=e$lo/1e6,xmax=e$hi/1e6,ymin=-Inf,ymax=Inf,
      fill='#FFA500',alpha=if(e$hi-e$lo>200000) .025 else .12)
    for(location in alternatives) plot <- plot+annotate('rect',xmin=location$lo/1e6,xmax=location$hi/1e6,
      ymin=-Inf,ymax=Inf,fill='#762A83',alpha=.08)
    plot+fill_scale()+coord_cartesian(xlim=window/1e6,clip='on')+base
  }
  mark <- function(plot) {
    # Fixed colours do not replace each panel's measurement legend.
    for(kind in names(marker_colours)) {
      x <- markers[markers$kind==kind,,drop=FALSE]
      if(nrow(x)) plot <- plot+geom_vline(data=x,aes(xintercept=point/1e6),
        colour=marker_colours[kind],linetype=if(kind=='Submitted cut') 'solid' else 'dashed',linewidth=if(kind=='Submitted cut') .8 else .5)
    }
    plot
  }
  peer <- decorate(ggplot(),FALSE)+geom_rect(data=tracks,
    aes(xmin=lo/1e6,xmax=hi/1e6,ymin=row-.35,ymax=row+.35,fill=chrom))+
    scale_y_reverse(breaks=seq_along(e$peer_tracks),labels=vapply(e$peer_tracks,`[[`,character(1),'peer'),
      limits=c(max(1,length(e$peer_tracks))+.7,-.8))+labs(x=NULL,y=NULL,
      title='Chromosome assignments; blank regions have no usable assignment')
  # Group nearly coincident labels as before; individual sites remain in the tables.
  label_rows <- markers[markers$kind=='R assay (not a cut)'&markers$point>=window[1]&markers$point<=window[2],,drop=FALSE]
  if(nrow(label_rows)) {
    label_rows <- label_rows[order(label_rows$point),]
    groups <- cumsum(c(TRUE,diff(label_rows$point)>=.03*(e$hi-e$lo+1000000)))
    label_rows <- do.call(rbind,lapply(split(label_rows,groups),function(g)
      data.frame(point=mean(g$point),label=paste(sub(paste0('^',id),'',g$label),collapse='/'))))
    peer <- peer+geom_text(data=label_rows,aes(x=point/1e6,y=-.45,label=label),colour='#762A83',size=2.6)
  }
  span <- decorate(ggplot())+geom_hline(yintercept=2,linetype='dotted',colour='grey50')+
    geom_line(data=normal,aes(position/1e6,molecules,colour='MAPQ ≥20'),linewidth=.45)+
    geom_line(data=strict,aes(position/1e6,molecules,colour='MAPQ ≥30'),linewidth=.45)+
    geom_line(data=allprimary,aes(position/1e6,molecules*span_scale$factor,colour='All primary · right axis'),linetype='dotted',linewidth=.5)+
    scale_colour_manual(values=c('MAPQ ≥20'='#2166AC','MAPQ ≥30'='#B35806','All primary · right axis'='#888888'),
      limits=c('MAPQ ≥20','MAPQ ≥30','All primary · right axis'),name=NULL)+
    scale_y_continuous(limits=c(0,span_scale$limit),sec.axis=sec_axis(~./span_scale$factor,name='All primary molecules'))+
    labs(x=NULL,y='HiFi spanning molecules')
  dep <- decorate(ggplot())+geom_line(data=depth,aes((start+end)/2e6,depth,colour='MAPQ ≥20'),linewidth=.45)+
    geom_line(data=depth_all,aes((start+end)/2e6,depth*depth_scale$factor,colour='All primary · right axis'),linetype='dotted',linewidth=.5)+
    scale_colour_manual(values=c('MAPQ ≥20'='#238B45','All primary · right axis'='#888888'),
      limits=c('MAPQ ≥20','All primary · right axis'),name=NULL)+
    scale_y_continuous(limits=c(0,depth_scale$limit),sec.axis=sec_axis(~./depth_scale$factor,name='All primary depth'))+
    labs(x=NULL,y='Mean HiFi depth')
  contact <- decorate(ggplot())+geom_line(data=contacts,aes(position/1e6,ratio,colour=library),linewidth=.5,na.rm=TRUE)+
    scale_colour_manual(values=library_colours,limits=libraries,name='Hi-C library')+scale_y_continuous(limits=c(0,NA))+
    labs(x=NULL,y='Hi-C contact support\nCross / sqrt(within L × R)',
      title='100 kb bins, 300 kb flanks · libraries shown separately; gaps mean unavailable')
  telo <- decorate(ggplot())+geom_line(data=telomeres,aes(start/1e6,count),colour='#762A83',linewidth=.35)+
    scale_y_continuous(limits=c(0,NA))+labs(x='Original pre-finishing scaffold coordinate (Mb)',y='Telomere motif count\nper 1 kb')
  missing_panel <- function(plot,unavailable) {
    if(unavailable) plot+annotate('text',x=mean(window)/1e6,y=.5,label='Evidence unavailable',size=3) else plot
  }
  span <- missing_panel(span,!nrow(normal)&&!nrow(strict)&&!nrow(allprimary))
  dep <- missing_panel(dep,!nrow(depth)&&!nrow(depth_all))
  contact <- missing_panel(contact,!nrow(contacts));telo <- missing_panel(telo,!nrow(telomeres))
  panels <- lapply(list(peer,span,dep,contact,telo),mark)
  for(i in 1:4) panels[[i]] <- panels[[i]]+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
  figure <- wrap_plots(panels,ncol=1,heights=c(max(2,length(e$peer_tracks)*.25),1.5,1.1,1.4,1))+
    plot_annotation(title=paste(data$assembly,id,e$scaffold,paste(e$left,'→',e$right),sep=' · '),
      caption=paste('Pale backgrounds: retained independent-individual consensus; ambiguous regions stay neutral. Yellow: detected transition bracket.',
        'Purple ranges/dashes: R assays, not cut recommendations. Red dashes: gap options. Solid green: submitted cuts.',
        'Spanning molecules are HiFi reads. Right axes show all primary placements on separate scales. Hi-C curves are descriptive, not calibrated cut tests.',sep='\n'),
      theme=theme(plot.caption=element_text(hjust=0,size=8.5),plot.title=element_text(size=13)))
  save_figure(figure,file.path(packet,'plots',paste0(id,if(detail) '.detail' else '')),13,max(12,9+.32*length(e$peer_tracks)))
}

# Reconstruct retained count matrices exactly from their sparse presentation tables.
spec_path <- file.path(packet,'junction',paste0(id,'.contact_plot_spec.json'))
if(file.exists(spec_path)) {
  specs <- fromJSON(spec_path,simplifyVector=FALSE)
  sites <- e$junction_evidence$sites
  for(s in specs) {
    path <- file.path(packet,'junction',s$file)
    cells <- read.delim(gzfile(path))
    if(nrow(cells)!=s$nonzero_cells||sum(cells$count)!=s$total||
       s$zero_cells+s$nonzero_cells!=s$bins^2) stop('Contact export changed retained counts: ',path)
    cells$x <- (s$lo+(cells$column+.5)*s$resolution)/1e6
    cells$y <- (s$lo+(cells$row+.5)*s$resolution)/1e6
    lo <- s$lo/1e6;hi <- (s$lo+s$bins*s$resolution)/1e6
    maximum <- max(1,vapply(Filter(function(x) x$resolution==s$resolution,specs),`[[`,numeric(1),'maximum'))
    colour_breaks <- unique(c(0,1,10^seq_len(max(1,floor(log10(maximum))))))
    colour_breaks <- colour_breaks[colour_breaks<=maximum]
    contact_map <- ggplot(cells,aes(x,y,fill=count))+
      annotate('rect',xmin=lo,xmax=hi,ymin=lo,ymax=hi,fill='#000004')+
      geom_tile(width=s$resolution/1e6,height=s$resolution/1e6)+
      scale_fill_viridis_c(option='magma',trans='log1p',limits=c(0,maximum),
        breaks=colour_breaks,
        name='Retained pairs\nlog(1 + count)',guide=guide_colourbar(barheight=grid::unit(45,'mm')))+
      coord_fixed(xlim=c(lo,hi),ylim=c(lo,hi),expand=FALSE)+theme_bw(base_size=11)+
      labs(title=paste(data$assembly,id,s$library,paste0(s$resolution/1000,' kb'),sep=' · '),
        subtitle='Raw counts; shared colour scale across libraries at this resolution',
        x='Original scaffold coordinate (Mb)',y='Original scaffold coordinate (Mb)',
        caption='Dark cells are measured zeros, not missing maps. Raw library depths differ; a zero alone does not establish a misjoin.')+
      theme(legend.position='right',plot.caption=element_text(size=8,hjust=0))
    visible <- Filter(function(x) x$point>=s$lo&&x$point<=s$lo+s$bins*s$resolution,sites)
    if(length(visible)) {
      probes <- data.frame(point=vapply(visible,`[[`,numeric(1),'point')/1e6,
        label=vapply(visible,`[[`,character(1),'id'))
      # Axis tick labels provide an exact ID/coordinate key without overlapping map labels.
      key <- paste(paste0(probes$label,' = ',format(probes$point,digits=7),' Mb'),collapse='; ')
      key <- paste(strwrap(key,width=105),collapse='\n')
      contact_map <- contact_map+geom_vline(data=probes,aes(xintercept=point),colour='cyan',alpha=.45,linewidth=.3)+
        geom_hline(data=probes,aes(yintercept=point),colour='cyan',alpha=.45,linewidth=.3)+
        labs(caption=paste(key,'Dark cells = measured zeros. Cyan = assay probes; green = submitted cuts. Raw counts do not normalize library depth.',sep='\n'))
    }
    if(length(cuts)) for(cut in cuts) contact_map <- contact_map+
      geom_vline(xintercept=as.numeric(cut$cut_bp)/1e6,colour='#00CC66',linewidth=.6)+
      geom_hline(yintercept=as.numeric(cut$cut_bp)/1e6,colour='#00CC66',linewidth=.6)
    save_figure(contact_map,file.path(packet,'junction',paste0(id,'.',s$library,'.',s$resolution)),9,8.5)
  }
}
write.table(data.frame(track=c('assigned_peer_intervals','consensus_intervals','hifi20','hifi30','hifi_all_primary','depth20','depth_all_primary','contact_profile','telomeres'),
  rows=c(nrow(tracks),nrow(regions),nrow(normal),nrow(strict),nrow(allprimary),nrow(depth),nrow(depth_all),nrow(contacts),nrow(telomeres))),
  file.path(packet,'plots',paste0(id,'.render_audit.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(packet,'plots','evidence_plot_R_session.txt'))
write_json(list(event=id,overview_window_bp=c(max(0,e$lo-500000),e$hi+500000),
  detail_window_bp=detail_window(),focal_colours=as.list(colours[focal]),
  spanning_scale=span_scale,depth_scale=depth_scale,
  markers=markers),file.path(packet,'plots',paste0(id,'.render_spec.json')),auto_unbox=TRUE,pretty=TRUE)
write.table(data.frame(process='MISASSEMBLY_REVIEW',tool=c('R','ggplot2','jsonlite','patchwork'),
  version=c(as.character(getRversion()),vapply(c('ggplot2','jsonlite','patchwork'),function(p) as.character(packageVersion(p)),character(1)))),
  file.path(packet,'plots','evidence_plot_versions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
