rm(list=ls())
dat <- read.table('./data/TRF2/Table-2A-figure.txt',sep='\t',stringsAsFactors=FALSE,header=TRUE)
p_dat <- read.table('./data/TRF2/Table-2A-p-values.txt',sep='\t',stringsAsFactors=FALSE,header=TRUE)
sig <- rep('NS',nrow(p_dat))
sig[p_dat$p.value < 0.05] <- '*'
sig[p_dat$p.value < 0.01] <- '**'
sig[p_dat$p.value < 0.001] <- '***'
sig <- paste(as.character(p_dat$chr),sig,sep=' | ')
names(sig) <- p_dat$chr

dat_df <- plyr::dlply(dat,'chr')
df_proper <- function(x){
	des = c(x[,2],x[,3])
	var = rep(c('TRF2\noverlapping','Extra\nTRF2'),each=nrow(x))
	chr = rep(x[1,1],length(var))
	df_o = data.frame(chr=chr,des=des,var=var,stringsAsFactors=FALSE)
	return(df_o)
}
dat_df_prop <- lapply(dat_df,df_proper)
for(i in 1:length(dat_df_prop)){
	dat_df_prop[[i]]$chr.sig <- sig[names(dat_df_prop)[i]]
}

dat_p <- plyr::ldply(dat_df_prop)[,-1]
chr_ord <- paste('chr',c(1:22,'X','Y'),sep='')
dat_p$chr.sig <- factor(dat_p$chr.sig, levels=sig[chr_ord])
dat_p$des.log2 <- log2(dat_p$des+1)
dat_p$var <- factor(dat_p$var, levels=c('TRF2\noverlapping','Extra\nTRF2'))
library(ggplot2)
p <- ggplot(dat_p,aes(y=des.log2,x=var)) +
	geom_boxplot(aes(fill=var)) +
	facet_wrap('chr.sig', scale='free_y') +
	ylab('log2(DNase enrichment)') +
	scale_fill_manual(values=c('#80808094','#80808094'),name='') +
	xlab('') +
	theme(axis.text.x=element_text(angle=45,hjust=1),legend.position='none')
png('./figures/Supp-Figure-2-boxplot.png',width=1000,height=1400,res=200)
print(p)
dev.off()

pdf('./figures/Supp-Figure-2-boxplot.pdf')
print(p)
dev.off()


