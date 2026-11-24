# Shared presentation categories. No inference or alignment changes are made here.
fcs_keys <- c('primary','repeat_expected','low_expected','repeat_other','inconclusive','taxonomy_flag','unreported')
fcs_labels <- c('Primary division','Repeat; expected division reported','Low GX coverage; expected division reported','Repeat / low GX coverage; other division reported','Other inconclusive results','Taxonomy contaminant / chimeric','Not reported')
fcs_colors <- setNames(c('#38734B','#237D94','#A0CDDA','#E5A644','#A8A9AB','#C54C54','#EEEEEE'),fcs_keys)
fcs_theme <- theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),legend.position='bottom',legend.title=element_blank(),plot.caption=element_text(hjust=0,size=9),strip.text=element_text(face='bold'),panel.spacing.x=grid::unit(1,'cm'),plot.margin=margin(12,18,12,12))
fcs_fill <- function() scale_fill_manual(values=fcs_colors,breaks=fcs_keys,labels=fcs_labels,drop=FALSE)
fcs_legend <- function() guides(fill=guide_legend(ncol=2,byrow=TRUE))

fcs_pies <- function(d,prefix) {
  for(k in c('length_bp','gc_fraction','mean_depth',paste0('cat_',fcs_keys,'_bp'))) d[[k]]<-as.numeric(d[[k]])
  if(any(abs(rowSums(d[paste0('cat_',fcs_keys,'_bp')])-d$length_bp)>.5)) stop('Category fractions do not cover each record')
  if(any(d$length_bp<=0)|any(d$mean_depth<0,na.rm=TRUE)) stop('Invalid record metrics')
  d$plottable<-is.finite(d$gc_fraction)&is.finite(d$mean_depth)&d$gc_fraction>=0&d$gc_fraction<=1
  assemblies<-sort(unique(d$cohort_assembly))
  scope<-do.call(rbind,lapply(assemblies,function(a){x<-d[d$cohort_assembly==a,];data.frame(assembly=a,contigs=nrow(x),plotted=sum(x$plottable),missing_gc=sum(!is.finite(x$gc_fraction)),missing_depth=sum(!is.finite(x$mean_depth)),zero_depth=sum(x$mean_depth==0,na.rm=TRUE),omitted=sum(!x$plottable))}))
  write.table(scope,paste0(prefix,'_scope.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  x<-d[d$plottable,];x<-x[order(-x$length_bp),]
  trans<-scales::pseudo_log_trans(base=10,sigma=.1)
  gx<-if(nrow(x))range(100*x$gc_fraction) else c(0,100)
  gy<-if(nrow(x))range(trans$transform(x$mean_depth)) else c(0,1)
  dx<-max(5,diff(gx));dy<-max(.3,diff(gy))
  xrange<-gx+c(-.07,.07)*dx;yrange<-gy+c(-.07,.07)*dy
  verts<-list();v<-0L
  if(nrow(x))for(i in seq_len(nrow(x))) {
    radius<-.025*sqrt(x$length_bp[i]/max(x$length_bp));angle<-0
    for(key in fcs_keys) {
      fraction<-x[[paste0('cat_',key,'_bp')]][i]/x$length_bp[i]
      if(fraction<=0)next
      angles<-seq(angle,angle+2*pi*fraction,length.out=max(3,ceiling(48*fraction)))
      v<-v+1L
      verts[[v]]<-data.frame(cohort_assembly=x$cohort_assembly[i],group=v,category=key,
        x=100*x$gc_fraction[i]+c(0,cos(angles),0)*radius*dx,
        y=trans$transform(x$mean_depth[i])+c(0,sin(angles),0)*radius*dy)
      angle<-angle+2*pi*fraction
    }
  }
  polygons<-if(length(verts))do.call(rbind,verts) else data.frame(cohort_assembly=character(),group=integer(),category=character(),x=numeric(),y=numeric())
  labels<-setNames(sprintf('%s\n%d / %d records plotted; %d omitted',scope$assembly,scope$plotted,scope$contigs,scope$omitted),scope$assembly)
  positive<-x$mean_depth[x$mean_depth>0]
  ticks<-if(length(positive))unique(c(0,10^seq(floor(log10(min(positive))),ceiling(log10(max(positive)))))) else c(0,1,10)
  ticks<-ticks[trans$transform(ticks)>=yrange[1]&trans$transform(ticks)<=yrange[2]]
  p<-ggplot(polygons,aes(x,y,group=group,fill=category))+
    geom_point(data=x,aes(100*gc_fraction,trans$transform(mean_depth)),inherit.aes=FALSE,size=.45,color='#5B666E')+geom_polygon(alpha=.9)+
    geom_blank(data=data.frame(cohort_assembly=assemblies,x=mean(xrange),y=mean(yrange)),aes(x,y),inherit.aes=FALSE)+
    facet_wrap(~cohort_assembly,ncol=2,labeller=as_labeller(labels),drop=FALSE)+fcs_fill()+fcs_legend()+
    scale_x_continuous(limits=xrange,breaks=scales::breaks_pretty(n=6))+
    scale_y_continuous(limits=yrange,breaks=trans$transform(ticks),labels=scales::label_number()(ticks))+
    coord_fixed(ratio=dx/dy)+labs(title='GC and read coverage: FCS evidence composition',x='GC (%)',y='Mean read depth',
    caption='One pie per complete record; pie area scales with record length. Wedges show base-pair proportions in the shared FCS evidence categories.\nSmall neutral centre marks keep tiny records visible. Expected division is the GX asserted division.\nRepeat / low GX coverage remain inconclusive; GX coverage is separate from read depth.')+fcs_theme
  height<-max(7,ceiling(length(assemblies)/2)*3.3+3)
  ggsave(paste0(prefix,'.png'),p,width=13,height=height,dpi=160,bg='white',limitsize=FALSE)
  ggsave(paste0(prefix,'.pdf'),p,width=13,height=height,bg='white',limitsize=FALSE)
}

fcs_composition <- function(d,file,detail=FALSE) {
  categories<-paste0('cat_',fcs_keys,'_bp')
  if(detail) {
    alln<-sum(d$is_mixed==1);d<-d[d$is_mixed==1,];d<-head(d[order(-d$length_bp,d$sequence),],30)
    d$label<-sprintf('%s | %.2f Mb',d$sequence,d$length_bp/1e6)
    totals<-data.frame(label=d$label,denominator=d$length_bp)
    if(!nrow(d)){p<-ggplot()+annotate('text',x=0,y=0,label='No mixed records in this assembly')+theme_void();ggsave(paste0(file,'.png'),p,width=10,height=3,bg='white');return(invisible(NULL))}
    data<-do.call(rbind,lapply(seq_len(nrow(d)),function(i)data.frame(label=d$label[i],stage=d$coordinate_stage[i],category=fcs_keys,pct=100*as.numeric(d[i,categories])/d$length_bp[i])))
    data$label<-factor(data$label,levels=rev(d$label))
    subtitle<-sprintf('Showing %d of %d mixed records; limit 30, longest first',nrow(d),alln)
    height<-max(5,nrow(d)*.25+3)
  } else {
    mixed<-d[d$is_mixed==1,]
    groups<-unique(d[c('cohort_assembly','coordinate_stage')])
    data<-do.call(rbind,lapply(seq_len(nrow(groups)),function(i){g<-groups[i,];z<-mixed[mixed$cohort_assembly==g$cohort_assembly&mixed$coordinate_stage==g$coordinate_stage,];denom<-sum(z$length_bp);data.frame(label=g$cohort_assembly,stage=g$coordinate_stage,category=fcs_keys,pct=if(denom>0)100*colSums(z[categories])/denom else rep(0,length(fcs_keys)))}))
    data$label<-factor(data$label,levels=rev(sort(unique(data$label))))
    subtitle<-'One row per assembly; bars partition mixed-record base pairs. Empty bars indicate no mixed records.'
    height<-max(6,length(unique(data$label))*.35+4)
  }
  data$category<-factor(data$category,levels=fcs_keys)
  data$stage<-factor(data$stage,levels=c('contig','scaffold','final'),labels=c('Contigs','Scaffolds','Final assemblies'))
  p<-ggplot(data,aes(label,pct,fill=category))+geom_col(width=.7)+coord_flip()+fcs_fill()+fcs_legend()+
    facet_grid(.~stage)+scale_y_continuous(limits=c(0,100),breaks=seq(0,100,25),labels=function(x)paste0(x,'%'),expand=expansion(mult=c(0,.02)))+
    labs(title=if(detail)'FCS evidence inside mixed records' else 'FCS evidence composition across assemblies',subtitle=subtitle,x=NULL,y=if(detail)'Base pairs (% of each complete record)' else 'Base pairs (% of mixed-record span)',caption='Expected-division repeat / low GX coverage remain inconclusive. Taxonomy flags are not cleaning recommendations.\nThe complete evidence table retains all records; display limits do not filter biological evidence.')+fcs_theme
  if(!detail) p<-p+geom_text(aes(group=interaction(label,category),label=ifelse(pct>=3,sprintf('%.1f%%',pct),''),color=category=='repeat_expected'),position=position_stack(vjust=.5),size=3)+scale_color_manual(values=c('FALSE'='#172A38','TRUE'='white'),guide='none')
  ggsave(paste0(file,'.png'),p,width=14,height=height,dpi=150,bg='white',limitsize=FALSE)
  ggsave(paste0(file,'.pdf'),p,width=14,height=height,bg='white',limitsize=FALSE)
}
