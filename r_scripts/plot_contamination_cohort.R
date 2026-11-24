# Presentation only: retain the accepted sequence evidence and taxonomic assignments.
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly=TRUE)
if(!length(args) %in% c(2,3)) stop('Supply sequence evidence, output directory and optional fcs mode')
fcs <- length(args)==3 && args[3]=='fcs'
d <- read.delim(args[1], check.names=FALSE, stringsAsFactors=FALSE)
required <- c('cohort_assembly','sequence','length_bp','gc_fraction','mean_depth','plot_taxon')
if(!all(required %in% names(d))) stop('Required GC/coverage evidence columns missing')
dir.create(args[2],recursive=TRUE,showWarnings=FALSE)
if(fcs && all(paste0('cat_',c('primary','repeat_expected','low_expected','repeat_other','inconclusive','taxonomy_flag','unreported'),'_bp') %in% names(d))) {
    source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'fcs_plot_common.R'))
    fcs_pies(d,file.path(args[2],'cohort_gc_coverage'))
    quit(save='no')
}
for(k in c('length_bp','gc_fraction','mean_depth')) d[[k]] <- as.numeric(d[[k]])
if(any(d$length_bp<0,na.rm=TRUE) || any(d$mean_depth<0,na.rm=TRUE)) stop('Negative sequence length or depth')
d$gc_pct <- 100*d$gc_fraction
d$plottable <- is.finite(d$gc_pct) & d$gc_pct>=0 & d$gc_pct<=100 & is.finite(d$mean_depth)
assemblies <- sort(unique(d$cohort_assembly))
scope <- do.call(rbind,lapply(assemblies,function(a) {
    x <- d[d$cohort_assembly==a,]
    data.frame(assembly=a,contigs=nrow(x),plotted=sum(x$plottable),
        missing_gc=sum(!is.finite(x$gc_pct)),missing_depth=sum(!is.finite(x$mean_depth)),
        zero_depth=sum(x$mean_depth==0,na.rm=TRUE),omitted=sum(!x$plottable))
}))
write.table(scope,file.path(args[2],'cohort_gc_coverage_scope.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
special <- c('no-hit'='#B1B7BD','unresolved'='#666E76','conflicting'='#A03B83',
             'not-searched'='#D8DCE0','unsearched'='#D8DCE0')
if(fcs) special <- c('unassigned'='#B1B7BD','unresolved: repeat'='#92999F',
                    'unresolved: low-coverage'='#BFC5C9','unresolved: inconclusive'='#666E76',
                    'mixed FCS classifications'='#A03B83')
span <- sort(tapply(d$length_bp,d$plot_taxon,sum,na.rm=TRUE),decreasing=TRUE)
taxa <- head(setdiff(names(span),names(special)),10)
palette <- c(setNames(grDevices::hcl.colors(max(1,length(taxa)),'Dark 3')[seq_along(taxa)],taxa),
             special,'other taxa'='#C98942')
d$display_taxon <- ifelse(d$plot_taxon %in% names(palette),d$plot_taxon,'other taxa')
labels <- setNames(sprintf('%s\n%s / %s contigs plotted; %s omitted',scope$assembly,
                          scope$plotted,scope$contigs,scope$omitted),scope$assembly)
points <- d[d$plottable,]
visible <- intersect(names(palette),unique(points$display_taxon))
gc_limits <- if(nrow(points)) c(max(0,5*floor((min(points$gc_pct)-2)/5)),
                               min(100,5*ceiling((max(points$gc_pct)+2)/5))) else c(0,100)
positive <- points$mean_depth[points$mean_depth>0]
breaks <- if(length(positive)) unique(c(0,10^seq(max(-3,floor(log10(min(positive)))),
                                                      max(-3,ceiling(log10(max(positive))))))) else c(0,1,10)
p <- ggplot(points,aes(gc_pct,mean_depth,color=display_taxon)) +
    geom_blank(data=data.frame(cohort_assembly=assemblies,gc_pct=50,mean_depth=1),
               aes(gc_pct,mean_depth),inherit.aes=FALSE) +
    geom_point(size=.5,alpha=.6) + geom_point(aes(size=length_bp),alpha=.55) +
    facet_wrap(~cohort_assembly,ncol=2,labeller=as_labeller(labels),drop=FALSE) +
    scale_x_continuous(limits=gc_limits,breaks=scales::breaks_pretty(n=6)) +
    scale_y_continuous(trans=scales::pseudo_log_trans(base=10,sigma=.1),breaks=breaks,
                       labels=scales::label_number()) +
    scale_color_manual(values=palette,limits=visible,drop=FALSE) +
    scale_size_area(max_size=7,labels=scales::label_number(scale_cut=scales::cut_short_scale())) +
    guides(color=guide_legend(ncol=if(fcs) 2 else 4,title.position='top',override.aes=list(size=3,alpha=1)),
           size=guide_legend(nrow=1,title.position='top')) +
    labs(title='GC and read coverage across assemblies',x='GC (%)',y='Mean read depth',
         color=if(fcs) 'FCS-GX division / unresolved status' else 'Taxonomic group / search status',size='Contig length (bp)',
         caption=paste('One dot per contig; shared scales and colours. Assemblies are shown separately; homologous contigs are not deduplicated.',
             if(fcs) 'Colours use FCS-GX classifications on the exact screened FASTA. Inconclusive and mixed assignments remain explicit.' else
             'Taxonomy uses retained search evidence. Unsearched sequence is not a searched no-hit.',sep='\n')) +
    theme_bw(base_size=11) + theme(legend.position='bottom',legend.box='vertical',strip.text=element_text(size=10),
                                 plot.caption=element_text(hjust=0),panel.grid.minor=element_blank())
height <- max(5,ceiling(length(assemblies)/2)*3+2)
ggsave(file.path(args[2],'cohort_gc_coverage.png'),p,width=12,height=height,dpi=160,bg='white',limitsize=FALSE)
ggsave(file.path(args[2],'cohort_gc_coverage.pdf'),p,width=12,height=height,bg='white',limitsize=FALSE)
