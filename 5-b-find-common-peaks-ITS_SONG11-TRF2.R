#######################
# author : Rintu Kutum
#######################
rm(list=ls())
its_bed <- read.table(
	'./data/TRF2/ITS.bed',
	header=FALSE,
	sep='\t',
	stringsAsFactors=FALSE)
song11_bed <- read.table(
	'./data/TRF2/SONG11.bed',
	header=FALSE,
	sep='\t',
	stringsAsFactors=FALSE)
trf2_bed <- read.table(
	'./data/TRF2/all_common_peaks_pankaj_sorted.bed',
	header=FALSE,
	sep='\t',
	stringsAsFactors=FALSE
)
checkSEswap <- function(s,e){
	if(s > e){
		sd <- e
		ed <- s
		# swap
		s <- sd
		e <- ed
	}
	return(c(s=s,e=e))
}
findCommonPeaks.TRF2 <- function(
	query,
	search
){
	# format
	colnames(query) <- colnames(search) <- c('chr','start','end')
	# sorting
	for(i in 1:nrow(query)){
		s <- query[i,2]
		e <- query[i,3]
		x <- checkSEswap(s,e)
		query[i,2:3] <- x
	}
	query.s <- plyr::ddply(query, 'chr', function(x){x[order(x$start),]})
	search.s <- plyr::ddply(search, 'chr', function(x){x[order(x$start),]})
	# genomic ranges class
	library('GenomicRanges')
	q.gr <- makeGRangesFromDataFrame(query.s)
	s.gr <- makeGRangesFromDataFrame(search.s)
	sq_overlap <- findOverlaps(s.gr,q.gr)
	sq_count <- countOverlaps(s.gr,q.gr)
	#
	q.coord <- q.gr[sq_overlap@to,]
	q.ov <- data.frame(q.coord)[,1:4]
	colnames(q.ov)[1] <- 'chr'
	colnames(q.ov)[-1] <- paste(colnames(q.ov)[-1],'-q',sep='')
	s.coord <- s.gr[sq_overlap@from,]
	s.ov <- data.frame(s.coord)[,2:4]
	colnames(s.ov) <- paste(colnames(s.ov),'-s',sep='')
	if(nrow(s.ov) == 0){
		s.ov[1,] <- NA
		q.ov[1,] <- NA
	}
	output <- data.frame(
		q.ov,
		s.ov
	)
	return(output)
}
dir.create('./data/intersect-ITS-SONG11',showWarnings=FALSE)
its_common <- findCommonPeaks.TRF2(its_bed,trf2_bed)
song11_common <- findCommonPeaks.TRF2(song11_bed,trf2_bed)

write.csv(
	its_common,
	'./data/intersect-ITS-SONG11/ITS-TRF2-intersect.csv',
	row.names=FALSE
	)

write.csv(
	song11_common,
	'./data/intersect-ITS-SONG11/SONG11-TRF2-intersect.csv',
	row.names=FALSE
	)
#-----------

