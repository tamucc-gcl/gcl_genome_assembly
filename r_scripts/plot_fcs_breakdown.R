suppressPackageStartupMessages(library(ggplot2))
source(file.path(dirname(sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])),'fcs_plot_common.R'))
args<-commandArgs(trailingOnly=TRUE)
d<-read.delim(args[1],check.names=FALSE)
s<-read.delim(args[2],check.names=FALSE)
out<-args[3]
fcs_composition(d,file.path(out,'cohort_fcs_composition'))
assemblies<-sort(unique(d$cohort_assembly))
for(i in seq_along(assemblies)) {
  a<-assemblies[i];x<-d[d$cohort_assembly==a,];prefix<-file.path(out,sprintf('assembly_%04d',i))
  fcs_pies(x,paste0(prefix,'_gc_coverage'))
  fcs_composition(x,paste0(prefix,'_composition'),TRUE)
  z<-s[s$cohort_assembly==a&s$is_mixed==1&!is.na(s$alignment_bp),]
  fields<-c('alignment_bp','transposon_mask_bp','low_complexity_mask_bp','conserved_mask_bp')
  labels<-c('GX database alignment coverage','GX transposon masking','GX low-complexity masking','GX conserved-sequence masking')
  if(nrow(z)) {
    den<-sum(z$span_bp);value<-vapply(fields,function(k)if(any(is.na(z[[k]])))NA_real_ else 100*sum(z[[k]])/den,numeric(1))
    q<-data.frame(metric=factor(labels,levels=rev(labels)),pct=value)
    p<-ggplot(q,aes(metric,pct))+geom_col(fill='#237D94',width=.6,na.rm=TRUE)+geom_text(aes(label=sprintf('%.1f%%',pct)),hjust=-.1,na.rm=TRUE)+coord_flip()+scale_y_continuous(limits=c(0,100))+labs(title=paste(a,'— GX alignment and masking'),x=NULL,y='Percent of reported GX segment span in mixed records',caption='Independent measurements: masking tracks may overlap and must not be summed. Missing measurements remain unavailable.\nGX coverage is separate from read depth; GX masking is not a separate repeat annotation.')+fcs_theme
  } else p<-ggplot()+annotate('text',x=0,y=0,label='No mixed-record alignment / masking measurements available')+theme_void()
  ggsave(paste0(prefix,'_diagnostics.png'),p,width=12,height=5,dpi=150,bg='white')
  ggsave(paste0(prefix,'_diagnostics.pdf'),p,width=12,height=5,bg='white')
}
