#----------------------
rm(list=ls())
source('./util-funcs.R')
bed.file <- './data/TRF2/all_common_peaks_pankaj_sorted.bed'
rep1.bam <- '../DNase-ENCODE-common-peaks/sorted_rep1_dedup.bam'
rep2.bam <- '../DNase-ENCODE-common-peaks/sorted_rep2_dedup.bam'
####
# rep1 coverage
getCoverageSamtools(
	bed.file = bed.file,
	bam.file = rep1.bam,
	output.name = './data/TRF2/TRF2_seq_rep1.bed'
)
####
# rep2 coverage
getCoverageSamtools(
	bed.file = bed.file,
	bam.file = rep2.bam,
	output.name = './data/TRF2/TRF2_seq_rep2.bed'
)
#############
rm(list=ls())
source('./util-funcs.R')
rep1_reads <- read_bed('./data/TRF2/TRF2_seq_rep1.bed',)
rep2_reads <- read_bed('./data/TRF2/TRF2_seq_rep2.bed')
reads_mean <- ((rep1_reads[,'V4'])/23 + (rep2_reads[,'V4'])/24)/2
#############
# TRF2 reads
trf2_reads <- data.frame(rep1_reads[,1:3],mean=reads_mean)
write.table(
	x = trf2_reads,
	file = './data/TRF2/TRF2_seq_reads.bed',
	col.names = FALSE,
	row.names = FALSE,
	sep = '\t',
	quote = FALSE
)
rm(list=ls())
source('./util-funcs.R')
trf2_reads <- read_bed('./data/TRF2/TRF2_seq_reads.bed')
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
trf2_seq <- read_bed('./data/TRF2/TRF2_peaks_coord_seq.bed')
motif_mismatch <- getMotifMaps(
	motif = 'TTAGGG',
	bed = trf2_seq,
	seq = trf2_seq$V4,
	max.mismatch = 1,
	flank = FALSE
)

getBox(unique(motif_mismatch$idx),trf2_reads)

coord_hit <- motif_mismatch
idx <- coord_hit$hit >= 4
coord_hit_idx <- plyr::dlply(coord_hit[idx,],'idx')
hits.status <- lapply(coord_hit_idx,getConsHits)
hits.status.df <- na.omit(plyr::ldply(hits.status))[,-1]
idx_cons <- unique(hits.status.df$idx)
getBox(idx_cons,flank_reads)

reads.df <- data.frame(
	matrix(
		data = NA,
		nrow = nrow(trf2_reads),
		ncol = 2
	),
	stringsAsFactors=FALSE
)
colnames(reads.df) <- c('cond','tag.reads')
reads.df$tag.reads <- trf2_reads$V4
reads.df$cond <- 'no-motif'
reads.df$cond[idx_cons] <- 'motif'
library(ggplot2)
p <- ggplot(reads.df,aes(x=as.factor(cond),y=log2(tag.reads+1)))+
	geom_boxplot() +
	xlab('')
pdf('./figures/motif-TRF2-boxplot.pdf',width=2.7,height=4)
p
dev.off()

