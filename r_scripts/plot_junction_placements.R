#!/usr/bin/env Rscript
# Render retained evidence only. Classification and mapping remain in Python assays.
suppressPackageStartupMessages({library(ggplot2); library(jsonlite); library(patchwork)})
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3) stop('Usage: plot_junction_placements.R packet event_id selected_cuts.json')
packet <- args[1]; id <- args[2]
evidence <- fromJSON(file.path(packet,'evidence.json'),simplifyVector=FALSE)
event <- Filter(function(e) e$id==id,evidence$events)[[1]]
root <- file.path(packet,'junction')
counts <- read.delim(file.path(root,paste0(id,'.placement_counts.tsv')),check.names=FALSE)
profile <- read.delim(file.path(root,'hifi/spanning_profile.tsv'))
profile <- profile[profile$event==id,]
cuts <- fromJSON(args[3],simplifyVector=FALSE)
palette <- c('Primary, any MAPQ'='#777777','Primary, MAPQ ≥30'='#0072B2',
             'Credible after competition'='#009E73','Reaches both unique arms'='#CC79A7',
             'Same-haplotype alternative'='#D55E00','Other-haplotype match (context)'='#6A51A3')
long <- do.call(rbind,lapply(seq_len(6),function(i) data.frame(
  site=counts$site,value=counts[[c('primary_all','primary_high','credible','unique_arms','competing','other_haplotype')[i]]],
  measurement=names(palette)[i])))
long$measurement <- factor(long$measurement,levels=names(palette))
long$site <- factor(long$site,levels=counts$site[order(counts$point)])
labels <- setNames(paste0(counts$site,'\n',format(counts$point/1e6,digits=6),' Mb\n',
  ifelse(is.na(counts$local_reads_long_enough),'Reach unknown',paste0('Long enough: ',counts$local_reads_long_enough))),counts$site)
base <- theme_bw(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position='bottom',
  plot.background=element_rect(fill='white',colour=NA))
# Keep exact marker coordinates, but spread their labels in coordinate order.
# Leader lines make even neighbouring gap/probe positions individually traceable.
markers <- counts[order(counts$point,counts$site),c('site','point')]
xrange <- range(profile$position/1e6)
xspan <- max(diff(xrange),.001)
markers$label_x <- markers$point/1e6
spacing <- min(.09*xspan,.92*xspan/max(1,nrow(markers)-1))
if(nrow(markers)>1) for(i in 2:nrow(markers))
  markers$label_x[i] <- max(markers$label_x[i],markers$label_x[i-1]+spacing)
markers$label_x <- markers$label_x-max(0,max(markers$label_x)-(xrange[2]-.04*xspan))
if(min(markers$label_x)<xrange[1]+.04*xspan)
  markers$label_x <- markers$label_x+(xrange[1]+.04*xspan-min(markers$label_x))
ymax <- max(1,profile$mapq0,profile$mapq30,na.rm=TRUE)
p <- ggplot(profile,aes(position/1e6))+geom_line(aes(y=mapq0,colour='Primary, any MAPQ'),linewidth=.5)+
  geom_line(aes(y=mapq30,colour='Primary, MAPQ ≥30'),linewidth=.6)+
  geom_segment(data=markers,aes(x=point/1e6,xend=point/1e6,y=0,yend=1.03*ymax),
    inherit.aes=FALSE,linetype='dashed',colour='#8E44AD',linewidth=.4)+
  geom_segment(data=markers,aes(x=point/1e6,xend=label_x,y=1.03*ymax,yend=1.17*ymax),
    inherit.aes=FALSE,colour='#8E44AD',linewidth=.4)+
  geom_label(data=markers,aes(x=label_x,y=1.21*ymax,label=site),inherit.aes=FALSE,
    colour='#71368A',fill='white',label.size=.25,size=3,label.padding=grid::unit(.15,'lines'))+
  coord_cartesian(ylim=c(0,1.29*ymax))+
  scale_colour_manual(values=palette,name=NULL)+labs(x='Original pre-finishing scaffold coordinate (Mb)',
  y='Primary HiFi spans',title='Regional context: primary placements before competitive checks')+base
if(length(cuts)) for(cut in cuts) p <- p+geom_vline(xintercept=as.numeric(cut$cut_bp)/1e6,colour='#006D2C',linewidth=.8)
q <- ggplot(long,aes(site,value,fill=measurement))+geom_col(position=position_dodge(width=.85),width=.8,na.rm=TRUE)+
  geom_text(aes(y=ifelse(is.na(value),0,value),label=ifelse(is.na(value),'NA',value)),
    position=position_dodge(width=.85),vjust=-.4,size=2.7)+
  scale_fill_manual(values=palette,name=NULL,drop=FALSE,guide=guide_legend(nrow=2,byrow=TRUE))+scale_x_discrete(labels=labels)+
  scale_y_continuous(expand=expansion(mult=c(0,.18)))+labs(x='Assay probe (not a proposed cut)',y='Molecules at probe',
  title='Placement checks at measured probes: compare counts, do not add them')+base
figure <- (p/q)+plot_layout(heights=c(1,1.2))+plot_annotation(
  title=paste(evidence$assembly,id,event$scaffold,paste(event$left,'→',event$right),sep=' · '),
  caption=paste('Spans require 1 kb on each side. Competitive checks are measured at probes only, not along the whole curve.',
    'Counts overlap. Same-haplotype alternatives may explain apparent spans; other-haplotype similarity alone is not a veto.',
    'NA = not measured (including uninterrupted assays at N gaps). Zero is an observed count, not a cut decision.',
    'Matching site IDs connect top markers to bottom bar groups; leaders preserve exact positions for neighbouring probes.',
    'Dashed purple: probes. Solid green: submitted cuts. Long enough: local reads whose length can reach both unique anchors; not proof they span.',sep='\n'),
  theme=theme(plot.caption=element_text(hjust=0,size=9),plot.title=element_text(size=13)))
figure_width <- max(13,nrow(markers)*1.25)
ggsave(file.path(root,paste0(id,'.placements.png')),figure,width=figure_width,height=9.5,dpi=160,bg='white')
ggsave(file.path(root,paste0(id,'.placements.pdf')),figure,width=figure_width,height=9.5,bg='white',device=cairo_pdf)
writeLines(capture.output(sessionInfo()),file.path(root,'placement_plot_R_session.txt'))
write.table(data.frame(process='MISASSEMBLY_REVIEW',tool=c('R','ggplot2','jsonlite','patchwork'),
  version=c(as.character(getRversion()),vapply(c('ggplot2','jsonlite','patchwork'),function(p) as.character(packageVersion(p)),character(1)))),
  file.path(root,'placement_plot_versions.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
