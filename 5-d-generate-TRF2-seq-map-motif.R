rm(list=ls())
trf2_bed <- read.table(
	'./data/TRF2/all_common_peaks_pankaj_sorted.bed',
	header=FALSE,
	sep='\t',
	stringsAsFactors=FALSE
)

### create the flanking regions from the mid peaks
getSeqRegions <- function(
	bed
){
	library('BSgenome.Hsapiens.UCSC.hg19')
	seq_bed <- bed
	seq <- c()
	for(i in 1:nrow(bed)){
		x <- bed[i,]
		print(i)
		chr <- as.character(x[1])
		start <- as.numeric(x[2])
		end <- as.numeric(x[3])
		seq[i] <- as.character(
			getSeq(
				Hsapiens,
				names=chr,
				start=start,
				end=end
			)
		)
	}
	seq_bed$seq <- seq
	return(seq_bed)
}
######
trf2_seq_bed <- getSeqRegions(trf2_bed)
write.table(
	trf2_seq_bed,
	'./data/TRF2/TRF2_seq.bed',
	col.names=FALSE,
	sep='\t',
	row.names=FALSE,
	quote=FALSE
)
rm(list=ls())
trf2_flank <- read.table(
	'./data/TRF2/TRF2_seq.bed',
	header=FALSE,
	sep= '\t',
	stringsAsFactors=FALSE
)
query_motif <- paste(rep('TTAGGG',4),collapse='')
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
		hit = length(hit@ranges),
		idx = i
	)
}

output.df <- plyr::ldply(output)
dir.create('./data/TRF2/motif-TTAGGG/',showWarnings=FALSE)
output.df <- output.df[,colnames(output.df)[c(6,1:4,7:10,5)],]
write.csv(
output.df,
'./data/TRF2/motif-TTAGGG/motif-TTAGGG-1-mismatch-TRF2-seq.csv',
row.names=FALSE
)
