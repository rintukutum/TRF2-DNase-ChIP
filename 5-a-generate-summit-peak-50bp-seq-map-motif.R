#############
# author : Rintu Kutum
#############
rm(list=ls())
trf2_bed <- read.table(
	'./data/TRF2/all_common_peaks_pankaj_sorted.bed',
	header=FALSE,
	sep='\t',
	stringsAsFactors=FALSE
)
### create the flanking regions from the mid peaks
getFlankingRegions <- function(
	bed,
	bp_span
){
	library('BSgenome.Hsapiens.UCSC.hg19')
	flank_bed <- bed
	seq <- c()
	for(i in 1:nrow(bed)){
		x <- bed[i,]
		print(i)
		chr <- as.character(x[1])
		start <- as.numeric(x[2])
		end <- as.numeric(x[3])
		mp <- round((end + start)/2)
		s <- mp - bp_span
		flank_bed[i,2] <- s
		e <- mp + bp_span
		flank_bed[i,3] <- e
		seq[i] <- as.character(
			getSeq(
				Hsapiens,
				names=chr,
				start=s,
				end=e
			)
		)
	}
	flank_bed$seq <- seq
	return(flank_bed)
}
######
flank_bed <- getFlankingRegions(
	trf2_bed,
	25
)
write.table(
	flank_bed,
	'./data/TRF2/flank_50bp_TRF2_summit.bed',
	col.names=FALSE,
	sep='\t',
	row.names=FALSE,
	quote=FALSE
)
rm(list=ls())
library(Biostrings)
trf2_flank <- read.table(
	'./data/TRF2/flank_50bp_TRF2_summit.bed',
	header=FALSE,
	sep= '\t',
	stringsAsFactors=FALSE
)
query_motif <- 'TTAGGG'
flank_seq <- trf2_flank$V4
mapped <- c()
for(i in 1:length(flank_seq)){
	print(i)
	mapped[i] <- countPattern(
		query_motif,
		flank_seq[i],
		max.mismatch=1
	)
}
idx <- which(mapped != 0)
test <- matchPattern(
	query_motif,
		flank_seq[i],
		max.mismatch=1
)
output <- list()
for(i in 1:length(idx)){
	hit <- matchPattern(
		query_motif,
		flank_seq[idx[i]],
		max.mismatch=1
	)
	print(i)
	output[[i]] <- cbind(
		data.frame(hit@ranges),
		query = query_motif,
		flank.seq = as.character(hit@subject),
		chr = trf2_flank[idx[i],1],
		flank.start = trf2_flank[idx[i],2],
		flank.end = trf2_flank[idx[i],3],
		hit = length(hit@ranges)
	)
}

output.df <- plyr::ldply(output)
dir.create('./data/TRF2/motif-TTAGGG/',showWarnings=FALSE)
write.csv(
output.df,
'./data/TRF2/motif-TTAGGG/motif-TTAGGG-1-mismatch-TRF2-summit-peak.csv',
row.names=FALSE
)

