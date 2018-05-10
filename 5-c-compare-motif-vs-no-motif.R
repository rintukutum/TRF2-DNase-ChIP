#############
rm(list=ls())
source('./util-funcs.R')
flank_reads <- read_bed('./data/TRF2/flank_50bp_TRF2_reads.bed')
getBox <- function(
	idx,
	flank_reads
){
	reads.df <- data.frame(
		matrix(
			data = NA,
			nrow = nrow(flank_reads),
			ncol = 2
		),
		stringsAsFactors=FALSE
	)
	colnames(reads.df) <- c('cond','tag.reads')
	reads.df$tag.reads <- flank_reads$V4
	reads.df$cond <- 'no-motif'
	reads.df$cond[idx] <- 'motif'
	library(ggplot2)
	p <- ggplot(reads.df,aes(x=as.factor(cond),y=log2(tag.reads+1)))+
		geom_boxplot()
	print(p)	
}
##############
###### 1 mismatch
motif_mismatch <- read.csv(
	'./data/TRF2/motif-TTAGGG/motif-TTAGGG-1-mismatch-TRF2-summit-peak.csv',
	header = TRUE,
	stringsAsFactors = FALSE
)
motif_no_mismatch <- read.csv(
	'./data/TRF2/motif-TTAGGG/motif-TTAGGG-NO-mismatch-TRF2-summit-peak.csv',
	header = TRUE,
	stringsAsFactors = FALSE
)
getBox(unique(motif_mismatch$idx),flank_reads)
getBox(unique(motif_no_mismatch$idx),flank_reads)

#########
motif_mapped <- read.csv(
'./data/TRF2/motif-TTAGGG/motif-TTAGGG-1-mismatch-TRF2-summit-peak.csv',
stringsAsFactors=FALSE
)
coord_hit <- motif_no_mismatch
idx <- coord_hit$hit >= 4
coord_hit_idx <- plyr::dlply(coord_hit[idx,],'idx')
source('./util-funcs.R')
hits.status <- lapply(coord_hit_idx,getConsHits)
hits.status.df <- na.omit(plyr::ldply(hits.status))[,-1]
idx_cons <- unique(hits.status.df$idx)
getBox(idx_cons,flank_reads)
